//
//  S2Loop.swift
//  S2Geometry
//
//  Created by Alex Studnicka on 7/1/16.
//  Copyright © 2016 Alex Studnicka. MIT License.
//

#if os(Linux)
	import Glibc
#else
	import Darwin.C
#endif

class S2LoopEdgeIndex: S2EdgeIndex {
  
  let vertices: [S2Point]
  /// Maps each S2Point to its order in the loop, from 1 to numVertices.
  let vertexToIndex: [S2Point: Int]
  /// The index (into "vertices") of the vertex that comes first in the total ordering of all vertices in this loop.
  let firstLogicalVertex: Int

  init?(vertices: [S2Point]) {
    guard vertices.count > 2 else { return nil }
    self.vertices = vertices
    guard let min = vertices.min() else { return nil }
    guard let first = vertices.index(of: min) else { return nil }
    firstLogicalVertex = first
    // make sure we don't have duplicates, which would be an invalid loop
    guard vertices.count == Set(vertices).count else { return nil }
    let s = zip(vertices, 1 ... vertices.count)
    let p = Dictionary(uniqueKeysWithValues: s)
    vertexToIndex = p
    super.init()
  }

  /// Return the index of a vertex at point "p", or -1 if not found. The return
  /// value is in the range 1..num_vertices_ if found.
  func find(vertex p: S2Point) -> Int? {
    return vertexToIndex[p]
  }
  
  /// Return the index of a vertex at point "p", or -1 if not found. The return
  /// value is in the range 1..num_vertices_ if found.
  func contains(vertex p: S2Point) -> Bool {
    return vertexToIndex[p] != nil
  }
  
  public override var numEdges: Int {
    return vertices.count
  }
  
  public override func edgeFrom(index: Int) -> S2Point {
    return vertices[index % vertices.count]
  }
  
  public override func edgeTo(index: Int) -> S2Point {
    return vertices[(index + 1) % vertices.count]
  }
  
}

/// An S2Loop represents a simple spherical polygon. It consists of a single
/// chain of vertices where the first vertex is implicitly connected to the last.
/// All loops are defined to have a CCW orientation, i.e. the interior of the
/// polygon is on the left side of the edges. This implies that a clockwise loop
/// enclosing a small area is interpreted to be a CCW loop enclosing a very large area.
///
/// Loops are not allowed to have any duplicate vertices (whether adjacent or
/// not), and non-adjacent edges are not allowed to intersect. Loops must have at
/// least 3 vertices. Although these restrictions are not enforced in optimized
/// code, you may get unexpected results if they are violated.
///
/// Point containment is defined such that if the sphere is subdivided into
/// faces (loops), every point is contained by exactly one face. This implies
/// that loops do not necessarily contain all (or any) of their vertices An
/// S2LatLngRect represents a latitude-longitude rectangle. It is capable of
/// representing the empty and full rectangles as well as single points.
public struct S2Loop {
	
	/// Max angle that intersections can be off by and yet still be considered colinear.
	public static let maxIntersectionError = 1e-15
	
	public let vertices: [S2Point]
  let originInside: Bool
  let bound: S2Rect

	/// Edge index used for performance-critical operations. For example,
  /// contains() can determine whether a point is inside a loop in nearly
  /// constant time, whereas without an edge index it is forced to compare the
  /// query point against every edge in the loop.
	let index: S2LoopEdgeIndex
  
	/// The depth of a loop is defined as its nesting level within its containing
  /// polygon. "Outer shell" loops have depth 0, holes within those loops have
  /// depth 1, shells within those holes have depth 2, etc. This field is only
  /// used by the S2Polygon implementation.
	public var depth: Int = 0
	
  // MARK: inits / factory
  
  init?(vertices: [S2Point], originInside: Bool, bound: S2Rect) {
    self.vertices = vertices
    self.originInside = originInside
    self.bound = bound
    guard let index = S2LoopEdgeIndex(vertices: vertices) else { return nil }
    self.index = index
  }
  
  static func contains(point p: S2Point, vertices: [S2Point], originInside: Bool) -> Bool {
    // Empty and full loops don't need a special case, but invalid loops with
    // zero vertices do, so we might as well handle them all at once.
    if vertices.count < 3 {
      return originInside
    }
    var inside = originInside
    guard let last = vertices.last else { fatalError() }
    var crosser = EdgeCrosser(a: S2Point.origin, b: p, c: last)
    // The s2edgeindex library is not optimized yet for long edges,
    // so the tradeoff to using it comes with larger loops.
    for v in vertices {
      inside = inside != crosser.isEdgeOrVertexCrossing(point: v)
    }
    return inside
  }
  
  static func computeOrigin(vertices: [S2Point]) -> Bool {
    switch vertices.count {
    case 1:
      // this is the special empty or full loop
      // the origin depends on if the vertex is in the southern hemisphere or not.
      return vertices[0].z < 0
    case 0, 2:
      // these are incomplete loops
      return false
    default:
      // Point containment testing is done by counting edge crossings starting
      // at a fixed point on the sphere (OriginPoint). We need to know whether
      // the reference point (OriginPoint) is inside or outside the loop before
      // we can construct the S2ShapeIndex. We do this by first guessing that
      // it is outside, and then seeing whether we get the correct containment
      // result for vertex 1. If the result is incorrect, the origin must be
      // inside the loop.
      // A loop with consecutive vertices A,B,C contains vertex B if and only if
      // the fixed vector R = B.Ortho is contained by the wedge ABC. The
      // wedge is closed at A and open at C, i.e. the point B is inside the loop
      // if A = R but not if C = R. This convention is required for compatibility
      // with VertexCrossing. (Note that we can't use OriginPoint
      // as the fixed vector because of the possibility that B == OriginPoint.)
      let v1Inside = S2Point.orderedCCW(a: vertices[1].ortho, b: vertices[0], c: vertices[2], o: vertices[1])
      return v1Inside != S2Loop.contains(point: vertices[1], vertices: vertices, originInside: false)
    }
  }

  static func computeBound(vertices: [S2Point], originInside: Bool) -> S2Rect {
    //
    switch vertices.count {
    case 1:
      // this is the special empty or full loop
      return originInside ? .full : .empty
    case 0, 2:
      // these are incomplete loops
      return .empty
    default: break
    }
    // We *must* call initBound before initIndex, because initBound calls
    // ContainsPoint(s2.Point), and ContainsPoint(s2.Point) does a bounds check whenever the
    // index is not fresh (i.e., the loop has been added to the index but the
    // index has not been updated yet).
    // The bounding rectangle of a loop is not necessarily the same as the
    // bounding rectangle of its vertices. First, the maximal latitude may be
    // attained along the interior of an edge. Second, the loop may wrap
    // entirely around the sphere (e.g. a loop that defines two revolutions of a
    // candy-cane stripe). Third, the loop may include one or both poles.
    // Note that a small clockwise loop near the equator contains both poles.
    var bounder = RectBounder()
    for p in vertices {
      bounder.add(point: p)
    }
    bounder.add(point: vertices[0])
    var b = bounder.bound
    if S2Loop.contains(point: S2Point(x: 0, y: 0, z: 1), vertices: vertices, originInside: originInside) {
      b = S2Rect(lat: R1Interval(lo: b.lat.lo, hi: .pi / 2.0), lng: S1Interval.full)
    }
    // If a loop contains the south pole, then either it wraps entirely
    // around the sphere (full longitude range), or it also contains the
    // north pole in which case b.Lng.isFull due to the test above.
    // Either way, we only need to do the south pole containment test if
    // b.Lng.isFull.
    if b.lng.isFull && S2Loop.contains(point: S2Point(x: 0, y: 0, z: -1), vertices: vertices, originInside: originInside) {
      b = S2Rect(lat: R1Interval(lo: -.pi / 2.0, hi: b.lat.hi), lng: S1Interval.full)
    }
    return b
  }
  
  /// Constructs a loop from the given points
  /// 1 point is a special empty or full, depending on it being in the northern hemisphere or not
  /// 0 or 2 points are incomplete loops
  public init?(points: [S2Point]) {
    let originInside = S2Loop.computeOrigin(vertices: points)
    let bound = S2Loop.computeBound(vertices: points, originInside: originInside)
    self.init(vertices: points, originInside: originInside, bound: bound)
  }
  
  /// Constructs a loop corresponding to the given cell.
  /// Note that the loop and cell *do not* contain exactly the same set of
  /// points, because Loop and Cell have slightly different definitions of
  /// point containment. For example, a Cell vertex is contained by all
  /// four neighboring Cells, but it is contained by exactly one of four
  /// Loops constructed from those cells. As another example, the cell
  /// coverings of cell and LoopFromCell(cell) will be different, because the
  /// loop contains points on its boundary that actually belong to other cells
  /// (i.e., the covering will include a layer of neighboring cells).
  public init(cell: S2Cell, bound: S2Rect? = nil) {
    let bound = bound ?? cell.rectBound
    let vertices = (0..<4).map { cell.vertex($0) }
    let originInside = S2Loop.computeOrigin(vertices: vertices)
    self.init(vertices: vertices, originInside: originInside, bound: bound)!
  }

  // MARK: tests
  
	/// Return true if this loop represents a hole in its containing polygon.
	public var isHole: Bool {
		return (depth & 1) != 0
	}
	
  /// Return true if the loop area is at most 2*Pi.
  public var isNormalized: Bool {
    // We allow a bit of error so that exact hemispheres are considered normalized.
    return area <= 2 * .pi + 1e-14
  }
  
  // MARK: computed members
  
	/// The sign of a loop is -1 if the loop represents a hole in its containing polygon, and +1 otherwise.
	public var sign: Int {
		return isHole ? -1 : 1
	}
	
	public var numVertices: Int {
		return vertices.count
	}
	
	/// For convenience, we make two entire copies of the vertex list available:
  /// vertex(n..2*n-1) is mapped to vertex(0..n-1), where n == numVertices().
	public func vertex(_ i: Int) -> S2Point {
		return vertices[i >= vertices.count ? i - vertices.count : i]
	}
	
	/// Invert the loop if necessary so that the area enclosed by the loop is at most 2*Pi.
	public func normalized() -> S2Loop {
		if isNormalized {
      return self
    }
    return inverted()
	}
	
	/// Reverse the order of the loop vertices, effectively complementing the region represented by the loop.
	public func inverted() -> S2Loop {
    let newVertices = Array(vertices.reversed())
		let newOriginInside = !originInside
    let newBound: S2Rect
		if bound.lat.lo > -0.5 * .pi && bound.lat.hi < 0.5 * .pi {
			// The complement of this loop contains both poles.
			newBound = .full
		} else {
      newBound = S2Loop.computeBound(vertices: newVertices, originInside: newOriginInside)
		}
    return S2Loop(vertices: newVertices, originInside: newOriginInside, bound: newBound)! // , newIndex, newFirstLogicalVertex, newVertexToIndex)
	}

	/// Helper method to get area and optionally centroid.
  func getAreaCentroid(doCentroid: Bool = true) -> (area: Double, centroid: S2Point) {
		var centroid = S2Point()
		// Don't crash even if loop is not well-defined.
		if numVertices < 3 {
			return (area: 0, centroid: centroid)
		}
		// The triangle area calculation becomes numerically unstable as the length
		// of any edge approaches 180 degrees. However, a loop may contain vertices
		// that are 180 degrees apart and still be valid, e.g. a loop that defines
		// the northern hemisphere using four points. We handle this case by using
		// triangles centered around an origin that is slightly displaced from the
		// first vertex. The amount of displacement is enough to get plenty of
		// accuracy for antipodal points, but small enough so that we still get
		// accurate areas for very tiny triangles.
		//
		// Of course, if the loop contains a point that is exactly antipodal from
		// our slightly displaced vertex, the area will still be unstable, but we
		// expect this case to be very unlikely (i.e. a polygon with two vertices on
		// opposite sides of the Earth with one of them displaced by about 2mm in
		// exactly the right direction). Note that the approximate point resolution
		// using the E7 or S2CellId representation is only about 1cm.
		var origin = vertex(0)
		let axis = (origin.largestAbsComponent + 1) % 3
		let slightlyDisplaced = origin.get(axis: axis) + M_E * 1e-10
		origin = S2Point(x: (axis == 0) ? slightlyDisplaced : origin.x, y: (axis == 1) ? slightlyDisplaced : origin.y, z: (axis == 2) ? slightlyDisplaced : origin.z)
		origin = S2Point.normalize(point: origin)
		var areaSum: Double = 0
		var centroidSum = S2Point()
		for i in 1 ... numVertices {
			areaSum += S2Point.signedArea(a: origin, b: vertex(i - 1), c: vertex(i));
			if doCentroid {
				// The true centroid is already premultiplied by the triangle area.
				let trueCentroid = S2Point.trueCentroid(a: origin, b: vertex(i - 1), c: vertex(i))
				centroidSum = centroidSum + trueCentroid
			}
		}
		// The calculated area at this point should be between -4*Pi and 4*Pi,
		// although it may be slightly larger or smaller than this due to
		// numerical errors.
		// assert (Math.abs(areaSum) <= 4 * S2..pi + 1e-12);
		if (areaSum < 0) {
			// If the area is negative, we have computed the area to the right of the
			// loop. The area to the left is 4*Pi - (-area). Amazingly, the centroid
			// does not need to be changed, since it is the negative of the integral
			// of position over the region to the right of the loop. This is the same
			// as the integral of position over the region to the left of the loop,
			// since the integral of position over the entire sphere is (0, 0, 0).
			areaSum += 4 * .pi
		}
		// The loop's sign() does not affect the return result and should be taken
		// into account by the caller.
		if doCentroid {
			centroid = centroidSum
		}
		return (area: areaSum, centroid: centroid)
	}

	/// Return the area of the loop interior, i.e. the region on the left side of
	/// the loop. The return value is between 0 and 4*Pi and the true centroid of
  /// the loop multiplied by the area of the loop (see S2.java for details on
  /// centroids). Note that the centroid may not be contained by the loop.
	public var areaAndCentroid: (area: Double, centroid: S2Point) {
		return getAreaCentroid()
	}
	
	/// Return the area of the polygon interior, i.e. the region on the left side
  /// of an odd number of loops. The return value is between 0 and 4*Pi.
	public var area: Double {
		return getAreaCentroid(doCentroid: false).area
	}
	
	/// Return the true centroid of the polygon multiplied by the area of the
  /// polygon (see {@link S2} for details on centroids). Note that the centroid
  /// may not be contained by the polygon.
	public var centroid: S2Point {
		return getAreaCentroid().centroid
	}
	
  // MARK: contains/ intersects
  
	// The following are the possible relationships between two loops A and B:
	//
	// (1) A and B do not intersect.
	// (2) A contains B.
	// (3) B contains A.
	// (4) The boundaries of A and B cross (i.e. the boundary of A
	// intersects the interior and exterior of B and vice versa).
	// (5) (A union B) is the entire sphere (i.e. A contains the
	// complement of B and vice versa).
	//
	// More than one of these may be true at the same time, for example if
	// A == B or A == Complement(B).
	
  /// Return true if the region contained by this loop is a superset of the region contained by the given other loop.
  public func contains(loop b: S2Loop) -> Bool {
    // For this loop A to contains the given loop B, all of the following must
    // be true:
    // (1) There are no edge crossings between A and B except at vertices.
    // (2) At every vertex that is shared between A and B, the local edge
    //     ordering implies that A contains B.
    // (3) If there are no shared vertices, then A must contain a vertex of B
    //     and B must not contain a vertex of A. (An arbitrary vertex may be
    //     chosen in each case.)
    // The second part of (3) is necessary to detect the case of two loops whose
    // union is the entire sphere, i.e. two loops that contains each other's
    // boundaries but not each other's interiors.
    if (!bound.contains(rect: b.rectBound)) { return false }
    // Unless there are shared vertices, we need to check whether A contains a
    // vertex of B. Since shared vertices are rare, it is more efficient to do
    // this test up front as a quick rejection test.
    if !contains(point: b.vertex(0)) && !index.contains(vertex: b.vertex(0)) { return false }
    // Now check whether there are any edge crossings, and also check the loop
    // relationship at any shared vertices.
    if (checkEdgeCrossings(b: b, relation: WedgeContains.self) <= 0) { return false }
    // At this point we know that the boundaries of A and B do not intersect,
    // and that A contains a vertex of B. However we still need to check for
    // the case mentioned above, where (A union B) is the entire sphere.
    // Normally this check is very cheap due to the bounding box precondition.
    if bound.union(rect: b.rectBound).isFull {
      if b.contains(point: vertex(0)) && !b.index.contains(vertex: vertex(0)) { return false }
    }
    return true
  }
  
	/// Return true if the region contained by this loop intersects the region contained by the given other loop.
	public func intersects(loop b: S2Loop) -> Bool {
		// a->Intersects(b) if and only if !a->Complement()->Contains(b).
		// This code is similar to Contains(), but is optimized for the case
		// where both loops enclose less than half of the sphere.
		if !bound.intersects(rect: b.rectBound) { return false }
		// Normalize the arguments so that B has a smaller longitude span than A.
		// This makes intersection tests much more efficient in the case where
		// longitude pruning is used (see CheckEdgeCrossings).
		if b.rectBound.lng.length > bound.lng.length { return b.intersects(loop: self) }
		// Unless there are shared vertices, we need to check whether A contains a
		// vertex of B. Since shared vertices are rare, it is more efficient to do
		// this test up front as a quick acceptance test.
		if contains(point: b.vertex(0)) && !index.contains(vertex: b.vertex(0)) { return true }
		// Now check whether there are any edge crossings, and also check the loop
		// relationship at any shared vertices.
		if (checkEdgeCrossings(b: b, relation: WedgeIntersects.self) < 0) { return true }
		// We know that A does not contain a vertex of B, and that there are no edge
		// crossings. Therefore the only way that A can intersect B is if B
		// entirely contains A. We can check this by testing whether B contains an
		// arbitrary non-shared vertex of A. Note that this check is cheap because
		// of the bounding box precondition and the fact that we normalized the
		// arguments so that A's longitude span is at least as long as B's.
		if b.rectBound.contains(rect: bound) {
      if b.contains(point: vertex(0)) && !b.index.contains(vertex: vertex(0)) { return true }
		}
		return false
	}
	
	/// Given two loops of a polygon, return true if A contains B. This version of
	/// contains() is much cheaper since it does not need to check whether the boundaries of the two loops cross.
	public func containsNested(loop b: S2Loop) -> Bool {
		if !bound.contains(rect: b.rectBound) { return false }
		// We are given that A and B do not share any edges, and that either one
		// loop contains the other or they do not intersect.
		guard let m = index.find(vertex: b.vertex(1)) else {
			// Since b->vertex(1) is not shared, we can check whether A contains it.
			return contains(point: b.vertex(1))
		}
		// Check whether the edge order around b->vertex(1) is compatible with
		// A containin B.
		return WedgeContains.test(a0: vertex(m - 1), ab1: vertex(m), a2: vertex(m + 1), b0: b.vertex(0), b2: b.vertex(2)) > 0
	}
	
	/// Return +1 if A contains B (i.e. the interior of B is a subset of the
  /// interior of A), -1 if the boundaries of A and B cross, and 0 otherwise.
  /// Requires that A does not properly contain the complement of B, i.e. A and B
  /// do not contain each other's boundaries. This method is used for testing
  /// whether multi-loop polygons contain each other.
	public func containsOrCrosses(loop b: S2Loop) -> Int {
		// There can be containment or crossing only if the bounds intersect.
		if !bound.intersects(rect: b.rectBound) { return 0 }
		// Now check whether there are any edge crossings, and also check the loop
		// relationship at any shared vertices. Note that unlike Contains() or
		// Intersects(), we can't do a point containment test as a shortcut because
		// we need to detect whether there are any edge crossings.
		let result = checkEdgeCrossings(b: b, relation: WedgeContainsOrCrosses.self)
		// If there was an edge crossing or a shared vertex, we know the result
		// already. (This is true even if the result is 1, but since we don't
		// bother keeping track of whether a shared vertex was seen, we handle this
		// case below.)
		if result <= 0 { return result }
		// At this point we know that the boundaries do not intersect, and we are
		// given that (A union B) is a proper subset of the sphere. Furthermore
		// either A contains B, or there are no shared vertices (due to the check
		// above). So now we just need to distinguish the case where A contains B
		// from the case where B contains A or the two loops are disjoint.
		if !bound.contains(rect: b.rectBound) { return 0 }
		if !contains(point: b.vertex(0)) && !index.contains(vertex: b.vertex(0)) { return 0 }
		return 1
	}
	
	/// Returns true if two loops have the same boundary except for vertex
  /// perturbations. More precisely, the vertices in the two loops must be in the
  /// same cyclic order, and corresponding vertex pairs must be separated by no
  /// more than maxError. Note: This method mostly useful only for testing
  /// purposes.
	internal func boundaryApproxEquals(boundary b: S2Loop, maxError: Double) -> Bool {
		if numVertices != b.numVertices { return false }
		let maxVertices = numVertices
		var iThis = index.firstLogicalVertex
		var iOther = b.index.firstLogicalVertex
		for _ in 0 ..< maxVertices {
			if !S2Point.approxEquals(vertex(iThis), b.vertex(iOther), maxError: maxError) { return false }
			iThis += 1
			iOther += 1
		}
		return true
	}

	/// The point 'p' does not need to be normalized.
	public func contains(point p: S2Point) -> Bool {
		if !bound.contains(point: p) {
			return false
		}
    var inside = originInside
		let origin = S2Point.origin
		var crosser = EdgeCrosser(a: origin, b: p, c: vertices[numVertices - 1])
		// The s2edgeindex library is not optimized yet for long edges,
		// so the tradeoff to using it comes with larger loops.
		if numVertices < 2000 {
			for i in 0 ..< numVertices {
				inside = inside != crosser.isEdgeOrVertexCrossing(point: vertices[i])
			}
		} else {
			var previousIndex = -2
			var it = getEdgeIterator(expectedQueries: numVertices)
			it.getCandidates(a: origin, b: p)
			for ai in it {
				if previousIndex != ai - 1 {
					crosser.restartAt(point: vertices[ai])
				}
				previousIndex = ai
				inside = inside != crosser.isEdgeOrVertexCrossing(point: vertex(ai + 1))
			}
		}
		return inside
	}

	/// Returns the shortest distance from a point P to this loop, given as the
	/// angle formed between P, the origin and the nearest point on the loop to P.
  /// This angle in radians is equivalent to the arclength along the unit sphere.
	public func getDistance(to p: S2Point) -> S1Angle {
		let normalized = S2Point.normalize(point: p)
		// The furthest point from p on the sphere is its antipode, which is an
		// angle of PI radians. This is an upper bound on the angle.
		var minDistance = S1Angle(radians: .pi)
		for i in 0 ..< numVertices {
			minDistance = min(minDistance, EdgeUtil.getDistance(x: normalized, a: vertex(i), b: vertex(i + 1)))
		}
		return minDistance
	}
	
	/// Creates an edge index over the vertices, which by itself takes no time.
  /// Then the expected number of queries is used to determine whether brute
  /// force lookups are likely to be slower than really creating an index, and if
  /// so, we do so. Finally an iterator is returned that can be used to perform
  /// edge lookups.
	func getEdgeIterator(expectedQueries: Int) -> S2EdgeIndex.DataEdgeIterator {
		index.predictAdditionalCalls(n: expectedQueries)
		return S2EdgeIndex.DataEdgeIterator(edgeIndex: index)
	}
	
	/// Return true if this loop is valid.
	public var isValid: Bool {
		if numVertices < 3 {
			return false
		}
		// All vertices must be unit length.
		for i in 0 ..< numVertices {
			if !vertex(i).isUnitLength {
				return false
			}
		}
		// Loops are not allowed to have any duplicate vertices.
		var vmap: Set<S2Point> = []
		for i in 0 ..< numVertices {
			if !vmap.insert(vertex(i)).inserted {
				return false
			}
		}
		// Non-adjacent edges are not allowed to intersect.
		var crosses = false
		var it = getEdgeIterator(expectedQueries: numVertices)
		for a1 in 0 ..< numVertices {
			let a2 = (a1 + 1) % numVertices
			var crosser = EdgeCrosser(a: vertex(a1), b: vertex(a2), c: vertex(0))
			var previousIndex = -2
			it.getCandidates(a: vertex(a1), b: vertex(a2))
			for b1 in it {
				let b2 = (b1 + 1) % numVertices
				// If either 'a' index equals either 'b' index, then these two edges
				// share a vertex. If a1==b1 then it must be the case that a2==b2, e.g.
				// the two edges are the same. In that case, we skip the test, since we
				// don't want to test an edge against itself. If a1==b2 or b1==a2 then
				// we have one edge ending at the start of the other, or in other words,
				// the edges share a vertex -- and in S2 space, where edges are always
				// great circle segments on a sphere, edges can only intersect at most
				// once, so we don't need to do further checks in that case either.
				if (a1 != b2 && a2 != b1 && a1 != b1) {
					// WORKAROUND(shakusa, ericv): S2.robustCCW() currently
					// requires arbitrary-precision arithmetic to be truly robust. That
					// means it can give the wrong answers in cases where we are trying
					// to determine edge intersections. The workaround is to ignore
					// intersections between edge pairs where all four points are
					// nearly colinear.
					let abc = S2Point.angle(a: vertex(a1), b: vertex(a2), c: vertex(b1))
					let abcNearlyLinear = S2Point.approxEquals(abc, 0, maxError: S2Loop.maxIntersectionError) || S2Point.approxEquals(abc, .pi, maxError: S2Loop.maxIntersectionError)
					let abd = S2Point.angle(a: vertex(a1), b: vertex(a2), c: vertex(b2));
					let abdNearlyLinear = S2Point.approxEquals(abd, 0, maxError: S2Loop.maxIntersectionError) || S2Point.approxEquals(abd, .pi, maxError: S2Loop.maxIntersectionError)
					if abcNearlyLinear && abdNearlyLinear { continue }
					if previousIndex != b1 {
						crosser.restartAt(point: vertex(b1))
					}
					// Beware, this may return the loop is valid if there is a
					// "vertex crossing".
					// TODO(user): Fix that.
					crosses = crosser.robustCrossing(point: vertex(b2)) > 0
					previousIndex = b2;
					if crosses {
//						log.info("Edges " + a1 + " and " + b1 + " cross");
//						log.info(String.format("Edge locations in degrees: " + "%s-%s and %s-%s",
//						new S2LatLng(vertex(a1)).toStringDegrees(),
//						new S2LatLng(vertex(a2)).toStringDegrees(),
//						new S2LatLng(vertex(b1)).toStringDegrees(),
//						new S2LatLng(vertex(b2)).toStringDegrees()));
						return false
					}
				}
			}
		}
		return true
	}

	/// This method encapsulates the common code for loop containment and
  /// intersection tests. It is used in three slightly different variations to
  /// implement contains(), intersects(), and containsOrCrosses().
  /// In a nutshell, this method checks all the edges of this loop (A) for
  /// intersection with all the edges of B. It returns -1 immediately if any edge
  /// intersections are found. Otherwise, if there are any shared vertices, it
  /// returns the minimum value of the given WedgeRelation for all such vertices
  /// (returning immediately if any wedge returns -1). Returns +1 if there are no
  /// intersections and no shared vertices.
	func checkEdgeCrossings(b: S2Loop, relation: WedgeRelation.Type) -> Int {
		var it = getEdgeIterator(expectedQueries: b.numVertices)
		var result = 1
		// since 'this' usually has many more vertices than 'b', use the index on
		// 'this' and loop over 'b'
		for j in 0 ..< b.numVertices {
			var crosser = EdgeCrosser(a: b.vertex(j), b: b.vertex(j + 1), c: vertex(0))
			var previousIndex = -2
			it.getCandidates(a: b.vertex(j), b: b.vertex(j + 1))
			for i in it {
				if previousIndex != i - 1 {
					crosser.restartAt(point: vertex(i))
				}
				previousIndex = i
				let crossing = crosser.robustCrossing(point: vertex(i + 1))
				if crossing < 0 { continue }
				if crossing > 0 { return -1 } // There is a proper edge crossing.
				if vertex(i + 1) == (b.vertex(j + 1)) {
					result = min(result, relation.test(a0: vertex(i), ab1: vertex(i + 1), a2: vertex(i + 2), b0: b.vertex(j), b2: b.vertex(j + 2)))
					if result < 0 { return result }
				}
			}
		}
		return result
	}
	
}

extension S2Loop: S2Region {
  
	// MARK: S2Region
	
	/// Return a bounding spherical cap.
	public var capBound: S2Cap {
		return bound.capBound
	}
	
	/// Return a bounding latitude-longitude rectangle.
	public var rectBound: S2Rect {
		return bound
	}
	
	/// If this method returns true, the region completely contains the given cell.
  /// Otherwise, either the region does not contain the cell or the containment
  /// relationship could not be determined.
	public func contains(cell: S2Cell) -> Bool {
		// It is faster to construct a bounding rectangle for an S2Cell than for
		// a general polygon. A future optimization could also take advantage of
		// the fact than an S2Cell is convex.
		let cellBound = cell.rectBound
		if !bound.contains(rect: cellBound) {
			return false
		}
		let cellLoop = S2Loop(cell: cell, bound: cellBound)
		return contains(loop: cellLoop)
	}
	
	/// If this method returns false, the region does not intersect the given cell.
	/// Otherwise, either region intersects the cell, or the intersection
  /// relationship could not be determined.
	public func mayIntersect(cell: S2Cell) -> Bool {
		// It is faster to construct a bounding rectangle for an S2Cell than for
		// a general polygon. A future optimization could also take advantage of
		// the fact than an S2Cell is convex.
		let cellBound = cell.rectBound
		if !bound.intersects(rect: cellBound) {
			return false
		}
		let cellLoop = S2Loop(cell: cell, bound: cellBound)
		return cellLoop.intersects(loop: self)
	}
	
}

extension S2Loop: Equatable, Comparable {

  public static func ==(lhs: S2Loop, rhs: S2Loop) -> Bool {
    guard lhs.numVertices == rhs.numVertices else { return false }
    // Compare the two loops' vertices, starting with each loop's
    // firstLogicalVertex. This allows us to always catch cases where logically
    // identical loops have different vertex orderings (e.g. ABCD and BCDA).
    let maxVertices = lhs.numVertices
    var iThis = lhs.index.firstLogicalVertex
    var iOther = rhs.index.firstLogicalVertex
    for _ in 0 ..< maxVertices {
      if lhs.vertex(iThis) != rhs.vertex(iOther) { return false }
      iThis += 1
      iOther += 1
    }
    return true
  }
  
  public static func <(lhs: S2Loop, rhs: S2Loop) -> Bool {
    if lhs.numVertices != rhs.numVertices {
      return lhs.numVertices < rhs.numVertices
    }
    // Compare the two loops' vertices, starting with each loop's
    // firstLogicalVertex. This allows us to always catch cases where logically
    // identical loops have different vertex orderings (e.g. ABCD and BCDA).
    let maxVertices = lhs.numVertices
    var iThis = lhs.index.firstLogicalVertex
    var iOther = rhs.index.firstLogicalVertex
    for _ in 0 ..< maxVertices {
      if lhs.vertex(iThis) != rhs.vertex(iOther) {
        return lhs.vertex(iThis) < rhs.vertex(iOther)
      }
      iThis += 1
      iOther += 1
    }
    return false
  }

}
