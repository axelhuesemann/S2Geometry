//
//  S2Point.swift
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

/// An S2Point represents a point on the unit sphere as a 3D vector. Usually
/// points are normalized to be unit length, but some methods do not require this.
public struct S2Point {

	public let x: Double
	public let y: Double
	public let z: Double
	
  // MARK: initialize
  
	public init(x: Double = 0, y: Double = 0, z: Double = 0) {
		self.x = x
		self.y = y
		self.z = z
	}
	
  /// Return a unique "origin" on the sphere for operations that need a fixed
  /// reference point. It should *not* be a point that is commonly used in edge
  /// tests in order to avoid triggering code to handle degenerate cases. (This
  /// rules out the north and south poles.)
  public static let origin = S2Point(x: 0, y: 1, z: 0)

  // MARK: computed members
  
	public var norm2: Double {
		return x * x + y * y + z * z
	}
	
	public var norm: Double {
		return sqrt(norm2)
	}
	
	/// Return a vector orthogonal to this one
	public var ortho: S2Point {
		let k = largestAbsComponent
		var temp: S2Point
		if k == 1 {
			temp = S2Point(x: 1, y: 0, z: 0)
		} else if k == 2 {
			temp = S2Point(x: 0, y: 1, z: 0)
		} else {
			temp = S2Point(x: 0, y: 0, z: 1)
		}
		return S2Point.normalize(point: crossProd(temp))
	}
	
	/// Return the index of the largest component fabs
	public var largestAbsComponent: Int {
		let temp = S2Point(x: abs(x), y: abs(y), z: abs(z))
		if temp.x > temp.y {
			if temp.x > temp.z {
				return 0
			} else {
				return 2
			}
		} else {
			if temp.y > temp.z {
				return 1
			} else {
				return 2
			}
		}
	}
	
	public static func normalize(point p: S2Point) -> S2Point {
		var norm = p.norm
		if norm != 0 {
			norm = 1.0 / norm
		}
		return p * norm
	}
	
	public func get(axis: Int) -> Double {
		return (axis == 0) ? x : (axis == 1) ? y : z
	}
	
	/// Return the angle between two vectors in radians
	public func angle(to x: S2Point) -> Double {
		return atan2((self × x).norm, self ⋅ x)
	}
	
	/// Compare two vectors, return true if all their components are within a difference of margin.
  public func approxEquals(_ point: S2Point, accuracy: Double = 1e-15) -> Bool {
    return (Swift.abs(x - point.x) < accuracy) && (Swift.abs(y - point.y) < accuracy) && (Swift.abs(z - point.z) < accuracy)
	}
	
  // MARK: arithmetic
	
	public func dotProd(_ b: S2Point) -> Double {
		return x * b.x + y * b.y + z * b.z
	}
	
	public func crossProd(_ b: S2Point) -> S2Point {
		return S2Point(x: y * b.z - z * b.y, y: z * b.x - x * b.z, z: x * b.y - y * b.x)
	}
	
}

extension S2Point: Equatable, Comparable, Hashable {
  
  // MARK: Equatable, Comparable, Hashable
  
  public static func ==(lhs: S2Point, rhs: S2Point) -> Bool {
    return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z
  }
  
  public static func <(lhs: S2Point, rhs: S2Point) -> Bool {
    if lhs.x < rhs.x {
      return true
    }
    if rhs.x < lhs.x {
      return false
    }
    if lhs.y < rhs.y {
      return true
    }
    if rhs.y < lhs.y {
      return false
    }
    if lhs.z < rhs.z {
      return true
    }
    return false
  }
  
  public var hashValue: Int {
    var value: UInt64 = 17
    let vx = value.addingReportingOverflow(Swift.abs(x).bitPattern).0
    let vy = value.addingReportingOverflow(Swift.abs(y).bitPattern).0
    let vz = value.addingReportingOverflow(Swift.abs(z).bitPattern).0
    value = value.addingReportingOverflow(vx.multipliedReportingOverflow(by: 37).0).0
    value = value.addingReportingOverflow(vy.multipliedReportingOverflow(by: 37).0).0
    value = value.addingReportingOverflow(vz.multipliedReportingOverflow(by: 37).0).0
    value ^= (value >> 32)
    return Int(Int64(bitPattern: value))
  }
  
  // MARK: arithmetic
  
  public prefix static func -(x: S2Point) -> S2Point {
    return S2Point(x: -x.x, y: -x.y, z: -x.z)
  }
  
  public static func +(lhs: S2Point, rhs: S2Point) -> S2Point {
    return S2Point(x: lhs.x + rhs.x, y: lhs.y + rhs.y, z: lhs.z + rhs.z)
  }
  
  public static func -(lhs: S2Point, rhs: S2Point) -> S2Point {
    return S2Point(x: lhs.x - rhs.x, y: lhs.y - rhs.y, z: lhs.z - rhs.z)
  }
  
  public static func *(lhs: S2Point, rhs: Double) -> S2Point {
    return S2Point(x: lhs.x * rhs, y: lhs.y * rhs, z: lhs.z * rhs)
  }
  
  public static func /(lhs: S2Point, rhs: Double) -> S2Point {
    return S2Point(x: lhs.x / rhs, y: lhs.y / rhs, z: lhs.z / rhs)
  }
  
  public static func ⋅(lhs: S2Point, rhs: S2Point) -> Double {
    return lhs.dotProd(rhs)
  }
  
  public static func ×(lhs: S2Point, rhs: S2Point) -> S2Point {
    return lhs.crossProd(rhs)
  }

}

extension S2Point {
  
  // MARK: geometry methods
  /// Return true if the given point is approximately unit length (this is mainly useful for assertions).
  public var isUnitLength: Bool {
    return Swift.abs(norm2 - 1) <= 1e-15
  }
  
  /**
   Return true if edge AB crosses CD at a point that is interior to both
   edges. Properties:
   
   (1) SimpleCrossing(b,a,c,d) == SimpleCrossing(a,b,c,d)
   (2) SimpleCrossing(c,d,a,b) == SimpleCrossing(a,b,c,d)
   */
  public static func simpleCrossing(a: S2Point, b: S2Point, c: S2Point, d: S2Point) -> Bool {
    // We compute SimpleCCW() for triangles ACB, CBD, BDA, and DAC. All
    // of these triangles need to have the same orientation (CW or CCW)
    // for an intersection to exist. Note that this is slightly more
    // restrictive than the corresponding definition for planar edges,
    // since we need to exclude pairs of line segments that would
    // otherwise "intersect" by crossing two antipodal points.
    let ab = a.crossProd(b)
    let cd = c.crossProd(d)
    let acb = -ab.dotProd(c)
    let cbd = -cd.dotProd(b)
    let bda = ab.dotProd(d)
    let dac = cd.dotProd(a)
    return (acb * cbd > 0) && (cbd * bda > 0) && (bda * dac > 0)
  }
  
  /**
   Return a vector "c" that is orthogonal to the given unit-length vectors "a"
   and "b". This function is similar to a.CrossProd(b) except that it does a
   better job of ensuring orthogonality when "a" is nearly parallel to "b",
   and it returns a non-zero result even when a == b or a == -b.
   
   It satisfies the following properties (RCP == RobustCrossProd):
   
   1. RCP(a,b) != 0 for all a, b
   2. RCP(b,a) == -RCP(a,b) unless a == b or a == -b
   3. RCP(-a,b) == -RCP(a,b) unless a == b or a == -b
   4. RCP(a,-b) == -RCP(a,b) unless a == b or a == -b
   */
  public static func robustCrossProd(a: S2Point, b: S2Point) -> S2Point {
    // The direction of a.CrossProd(b) becomes unstable as (a + b) or (a - b)
    // approaches zero. This leads to situations where a.CrossProd(b) is not
    // very orthogonal to "a" and/or "b". We could fix this using Gram-Schmidt,
    // but we also want b.RobustCrossProd(a) == -b.RobustCrossProd(a).
    //
    // The easiest fix is to just compute the cross product of (b+a) and (b-a).
    // Given that "a" and "b" are unit-length, this has good orthogonality to
    // "a" and "b" even if they differ only in the lowest bit of one component.
    // assert (isUnitLength(a) && isUnitLength(b));
    let x = (b + a).crossProd(b - a)
    if x != S2Point() {
      return x
    }
    // The only result that makes sense mathematically is to return zero, but
    // we find it more convenient to return an arbitrary orthogonal vector.
    return a.ortho
  }
  
  /**
   Return the area of triangle ABC. The method used is about twice as
   expensive as Girard's formula, but it is numerically stable for both large
   and very small triangles. The points do not need to be normalized. The area
   is always positive.
   
   The triangle area is undefined if it contains two antipodal points, and
   becomes numerically unstable as the length of any edge approaches 180 degrees.
   */
  static func area(a: S2Point, b: S2Point, c: S2Point) -> Double {
    // This method is based on l'Huilier's theorem,
    //
    // tan(E/4) = sqrt(tan(s/2) tan((s-a)/2) tan((s-b)/2) tan((s-c)/2))
    //
    // where E is the spherical excess of the triangle (i.e. its area),
    // a, b, c, are the side lengths, and
    // s is the semiperimeter (a + b + c) / 2 .
    //
    // The only significant source of error using l'Huilier's method is the
    // cancellation error of the terms (s-a), (s-b), (s-c). This leads to a
    // *relative* error of about 1e-16 * s / min(s-a, s-b, s-c). This compares
    // to a relative error of about 1e-15 / E using Girard's formula, where E is
    // the true area of the triangle. Girard's formula can be even worse than
    // this for very small triangles, e.g. a triangle with a true area of 1e-30
    // might evaluate to 1e-5.
    //
    // So, we prefer l'Huilier's formula unless dmin < s * (0.1 * E), where
    // dmin = min(s-a, s-b, s-c). This basically includes all triangles
    // except for extremely long and skinny ones.
    //
    // Since we don't know E, we would like a conservative upper bound on
    // the triangle area in terms of s and dmin. It's possible to show that
    // E <= k1 * s * sqrt(s * dmin), where k1 = 2*sqrt(3)/Pi (about 1).
    // Using this, it's easy to show that we should always use l'Huilier's
    // method if dmin >= k2 * s^5, where k2 is about 1e-2. Furthermore,
    // if dmin < k2 * s^5, the triangle area is at most k3 * s^4, where
    // k3 is about 0.1. Since the best case error using Girard's formula
    // is about 1e-15, this means that we shouldn't even consider it unless
    // s >= 3e-4 or so.
    // We use volatile doubles to force the compiler to truncate all of these
    // quantities to 64 bits. Otherwise it may compute a value of dmin > 0
    // simply because it chose to spill one of the intermediate values to
    // memory but not one of the others.
    let sa = b.angle(to: c)
    let sb = c.angle(to: a)
    let sc = a.angle(to: b)
    let s = 0.5 * (sa + sb + sc)
    if s >= 3e-4 {
      // Consider whether Girard's formula might be more accurate.
      let s2 = s * s
      let dmin = s - max(sa, max(sb, sc))
      if dmin < 1e-2 * s * s2 * s2 {
        // This triangle is skinny enough to consider Girard's formula.
        let area = girardArea(a: a, b: b, c: c)
        if dmin < s * (0.1 * area) {
          return area
        }
      }
    }
    // Use l'Huilier's formula.
    return 4 * atan(sqrt(max(0.0, tan(0.5 * s) * tan(0.5 * (s - sa)) * tan(0.5 * (s - sb)) * tan(0.5 * (s - sc)))))
  }
  
  /**
   * Return the area of the triangle computed using Girard's formula. This is
   * slightly faster than the Area() method above is not accurate for very small
   * triangles.
   */
  public static func girardArea(a: S2Point, b: S2Point, c: S2Point) -> Double {
    // This is equivalent to the usual Girard's formula but is slightly
    // more accurate, faster to compute, and handles a == b == c without
    // a special case.
    let ab = a.crossProd(b)
    let bc = b.crossProd(c)
    let ac = a.crossProd(c)
    return max(0.0, ab.angle(to: ac) - ab.angle(to: bc) + bc.angle(to: ac))
  }
  
  /// Like Area(), but returns a positive value for counterclockwise triangles and a negative value otherwise.
  public static func signedArea(a: S2Point, b: S2Point, c: S2Point) -> Double {
    return area(a: a, b: b, c: c) * Double(robustCCW(a: a, b: b, c: c))
  }
  
  // About centroids:
  // ----------------
  //
  // There are several notions of the "centroid" of a triangle. First, there
  // // is the planar centroid, which is simply the centroid of the ordinary
  // (non-spherical) triangle defined by the three vertices. Second, there is
  // the surface centroid, which is defined as the intersection of the three
  // medians of the spherical triangle. It is possible to show that this
  // point is simply the planar centroid projected to the surface of the
  // sphere. Finally, there is the true centroid (mass centroid), which is
  // defined as the area integral over the spherical triangle of (x,y,z)
  // divided by the triangle area. This is the point that the triangle would
  // rotate around if it was spinning in empty space.
  //
  // The best centroid for most purposes is the true centroid. Unlike the
  // planar and surface centroids, the true centroid behaves linearly as
  // regions are added or subtracted. That is, if you split a triangle into
  // pieces and compute the average of their centroids (weighted by triangle
  // area), the result equals the centroid of the original triangle. This is
  // not true of the other centroids.
  //
  // Also note that the surface centroid may be nowhere near the intuitive
  // "center" of a spherical triangle. For example, consider the triangle
  // with vertices A=(1,eps,0), B=(0,0,1), C=(-1,eps,0) (a quarter-sphere).
  // The surface centroid of this triangle is at S=(0, 2*eps, 1), which is
  // within a distance of 2*eps of the vertex B. Note that the median from A
  // (the segment connecting A to the midpoint of BC) passes through S, since
  // this is the shortest path connecting the two endpoints. On the other
  // hand, the true centroid is at M=(0, 0.5, 0.5), which when projected onto
  // the surface is a much more reasonable interpretation of the "center" of
  // this triangle.
  
  /**
   Return the centroid of the planar triangle ABC. This can be normalized to
   unit length to obtain the "surface centroid" of the corresponding spherical
   triangle, i.e. the intersection of the three medians. However, note that
   for large spherical triangles the surface centroid may be nowhere near the
   intuitive "center" (see example above).
   */
  public static func planarCentroid(a: S2Point, b: S2Point, c: S2Point) -> S2Point {
    return S2Point(x: (a.x + b.x + c.x) / 3.0, y: (a.y + b.y + c.y) / 3.0, z: (a.z + b.z + c.z) / 3.0)
  }
  
  /**
   Returns the true centroid of the spherical triangle ABC multiplied by the
   signed area of spherical triangle ABC. The reasons for multiplying by the
   signed area are (1) this is the quantity that needs to be summed to compute
   the centroid of a union or difference of triangles, and (2) it's actually
   easier to calculate this way.
   */
  public static func trueCentroid(a: S2Point, b: S2Point, c: S2Point) -> S2Point {
    // I couldn't find any references for computing the true centroid of a
    // spherical triangle... I have a truly marvellous demonstration of this
    // formula which this margin is too narrow to contain :)
    // assert (isUnitLength(a) && isUnitLength(b) && isUnitLength(c));
    let sina = b.crossProd(c).norm
    let sinb = c.crossProd(a).norm
    let sinc = a.crossProd(b).norm
    let ra = (sina == 0) ? 1 : (asin(sina) / sina)
    let rb = (sinb == 0) ? 1 : (asin(sinb) / sinb)
    let rc = (sinc == 0) ? 1 : (asin(sinc) / sinc)
    // Now compute a point M such that M.X = rX * det(ABC) / 2 for X in A,B,C.
    let x = S2Point(x: a.x, y: b.x, z: c.x)
    let y = S2Point(x: a.y, y: b.y, z: c.y)
    let z = S2Point(x: a.z, y: b.z, z: c.z)
    let r = S2Point(x: ra, y: rb, z: rc)
    return S2Point(x: 0.5 * y.crossProd(z).dotProd(r), y: 0.5 * z.crossProd(x).dotProd(r), z: 0.5 * x.crossProd(y).dotProd(r))
  }
  
  /**
   Return true if the points A, B, C are strictly counterclockwise. Return
   false if the points are clockwise or colinear (i.e. if they are all
   contained on some great circle).
   
   Due to numerical errors, situations may arise that are mathematically
   impossible, e.g. ABC may be considered strictly CCW while BCA is not.
   However, the implementation guarantees the following:
   
   If SimpleCCW(a,b,c), then !SimpleCCW(c,b,a) for all a,b,c.
   
   In other words, ABC and CBA are guaranteed not to be both CCW
   */
  public static func simpleCCW(a: S2Point, b: S2Point, c: S2Point) -> Bool {
    // We compute the signed volume of the parallelepiped ABC. The usual
    // formula for this is (AxB).C, but we compute it here using (CxA).B
    // in order to ensure that ABC and CBA are not both CCW. This follows
    // from the following identities (which are true numerically, not just
    // mathematically):
    //
    // (1) x.CrossProd(y) == -(y.CrossProd(x))
    // (2) (-x).DotProd(y) == -(x.DotProd(y))
    return c.crossProd(a).dotProd(b) > 0
  }
  
  /**
   WARNING! This requires arbitrary precision arithmetic to be truly robust.
   This means that for nearly colinear AB and AC, this function may return the
   wrong answer.
   
   Like SimpleCCW(), but returns +1 if the points are counterclockwise and -1
   if the points are clockwise. It satisfies the following conditions:
   
   (1) RobustCCW(a,b,c) == 0 if and only if a == b, b == c, or c == a (2)
   RobustCCW(b,c,a) == RobustCCW(a,b,c) for all a,b,c (3) RobustCCW(c,b,a)
   ==-RobustCCW(a,b,c) for all a,b,c
   
   In other words:
   
   (1) The result is zero if and only if two points are the same. (2)
   Rotating the order of the arguments does not affect the result. (3)
   Exchanging any two arguments inverts the result.
   
   This function is essentially like taking the sign of the determinant of
   a,b,c, except that it has additional logic to make sure that the above
   properties hold even when the three points are coplanar, and to deal with
   the limitations of floating-point arithmetic.
   
   Note: a, b and c are expected to be of unit length. Otherwise, the results
   are undefined.
   */
  public static func robustCCW(a: S2Point, b: S2Point, c: S2Point) -> Int {
    return robustCCW(a: a, b: b, c: c, aCrossB: a.crossProd(b))
  }
  
  /// A more efficient version of RobustCCW that allows the precomputed
  /// cross-product of A and B to be specified.
  /// Note: a, b and c are expected to be of unit length. Otherwise, the results
  /// are undefined
  public static func robustCCW(a: S2Point, b: S2Point, c: S2Point, aCrossB: S2Point) -> Int {
    // assert (isUnitLength(a) && isUnitLength(b) && isUnitLength(c));
    // There are 14 multiplications and additions to compute the determinant
    // below. Since all three points are normalized, it is possible to show
    // that the average rounding error per operation does not exceed 2**-54,
    // the maximum rounding error for an operation whose result magnitude is in
    // the range [0.5,1). Therefore, if the absolute value of the determinant
    // is greater than 2*14*(2**-54), the determinant will have the same sign
    // even if the arguments are rotated (which produces a mathematically
    // equivalent result but with potentially different rounding errors).
    let kMinAbsValue = 1.6e-15 // 2 * 14 * 2**-54
    let det = aCrossB.dotProd(c)
    // Double-check borderline cases in debug mode.
    // assert ((Math.abs(det) < kMinAbsValue) || (Math.abs(det) > 1000 * kMinAbsValue)
    //    || (det * expensiveCCW(a, b, c) > 0));
    if (det > kMinAbsValue) {
      return 1
    }
    if (det < -kMinAbsValue) {
      return -1
    }
    return expensiveCCW(a: a, b: b, c: c)
  }
  
  /// A relatively expensive calculation invoked by RobustCCW() if the sign of the determinant is uncertain.
  private static func expensiveCCW(a: S2Point, b: S2Point, c: S2Point) -> Int {
    // Return zero if and only if two points are the same. This ensures (1).
    if a == b || b == c || c == a {
      return 0
    }
    // Now compute the determinant in a stable way. Since all three points are
    // unit length and we know that the determinant is very close to zero, this
    // means that points are very nearly colinear. Furthermore, the most common
    // situation is where two points are nearly identical or nearly antipodal.
    // To get the best accuracy in this situation, it is important to
    // immediately reduce the magnitude of the arguments by computing either
    // A+B or A-B for each pair of points. Note that even if A and B differ
    // only in their low bits, A-B can be computed very accurately. On the
    // other hand we can't accurately represent an arbitrary linear combination
    // of two vectors as would be required for Gaussian elimination. The code
    // below chooses the vertex opposite the longest edge as the "origin" for
    // the calculation, and computes the different vectors to the other two
    // vertices. This minimizes the sum of the lengths of these vectors.
    //
    // This implementation is very stable numerically, but it still does not
    // return consistent results in all cases. For example, if three points are
    // spaced far apart from each other along a great circle, the sign of the
    // result will basically be random (although it will still satisfy the
    // conditions documented in the header file). The only way to return
    // consistent results in all cases is to compute the result using
    // arbitrary-precision arithmetic. I considered using the Gnu MP library,
    // but this would be very expensive (up to 2000 bits of precision may be
    // needed to store the intermediate results) and seems like overkill for
    // this problem. The MP library is apparently also quite particular about
    // compilers and compilation options and would be a pain to maintain.
    // We want to handle the case of nearby points and nearly antipodal points
    // accurately, so determine whether A+B or A-B is smaller in each case.
    let sab: Double = (a.dotProd(b) > 0) ? -1 : 1
    let sbc: Double = (b.dotProd(c) > 0) ? -1 : 1
    let sca: Double = (c.dotProd(a) > 0) ? -1 : 1
    let vab = a + (b * sab)
    let vbc = b + (c * sbc)
    let vca = c + (a * sca)
    let dab = vab.norm2
    let dbc = vbc.norm2
    let dca = vca.norm2
    // Sort the difference vectors to find the longest edge, and use the
    // opposite vertex as the origin. If two difference vectors are the same
    // length, we break ties deterministically to ensure that the symmetry
    // properties guaranteed in the header file will be true.
    var sign: Double
    if dca < dbc || (dca == dbc && a < b) {
      if dab < dbc || (dab == dbc && a < c) {
        // The "sab" factor converts A +/- B into B +/- A.
        sign = vab.crossProd(vca).dotProd(a) * sab // BC is longest edge
      } else {
        sign = vca.crossProd(vbc).dotProd(c) * sca // AB is longest edge
      }
    } else {
      if dab < dca || (dab == dca && b < c) {
        sign = vbc.crossProd(vab).dotProd(b) * sbc // CA is longest edge
      } else {
        sign = vca.crossProd(vbc).dotProd(c) * sca // AB is longest edge
      }
    }
    if sign > 0 {
      return 1
    }
    if sign < 0 {
      return -1
    }
    // The points A, B, and C are numerically indistinguishable from coplanar.
    // This may be due to roundoff error, or the points may in fact be exactly
    // coplanar. We handle this situation by perturbing all of the points by a
    // vector (eps, eps**2, eps**3) where "eps" is an infinitesmally small
    // positive number (e.g. 1 divided by a googolplex). The perturbation is
    // done symbolically, i.e. we compute what would happen if the points were
    // perturbed by this amount. It turns out that this is equivalent to
    // checking whether the points are ordered CCW around the origin first in
    // the Y-Z plane, then in the Z-X plane, and then in the X-Y plane.
    var ccw = R2Vector.planarOrderedCCW(a: R2Vector(x: a.y, y: a.z), b: R2Vector(x: b.y, y: b.z), c: R2Vector(x: c.y, y: c.z))
    if (ccw == 0) {
      ccw = R2Vector.planarOrderedCCW(a: R2Vector(x: a.z, y: a.x), b: R2Vector(x: b.z, y: b.x), c: R2Vector(x: c.z, y: c.x))
      if (ccw == 0) {
        ccw = R2Vector.planarOrderedCCW(a: R2Vector(x: a.x, y: a.y), b: R2Vector(x: b.x, y: b.y), c: R2Vector(x: c.x, y: c.y))
        // assert (ccw != 0);
      }
    }
    return ccw
  }
  
  /**
   Return true if the edges OA, OB, and OC are encountered in that order while
   sweeping CCW around the point O. You can think of this as testing whether
   A <= B <= C with respect to a continuous CCW ordering around O.
   
   Properties:
   - If orderedCCW(a,b,c,o) && orderedCCW(b,a,c,o), then a == b
   - If orderedCCW(a,b,c,o) && orderedCCW(a,c,b,o), then b == c
   - If orderedCCW(a,b,c,o) && orderedCCW(c,b,a,o), then a == b == c
   - If a == b or b == c, then orderedCCW(a,b,c,o) is true
   - Otherwise if a == c, then orderedCCW(a,b,c,o) is false
   */
  public static func orderedCCW(a: S2Point, b: S2Point, c: S2Point, o: S2Point) -> Bool {
    // The last inequality below is ">" rather than ">=" so that we return true
    // if A == B or B == C, and otherwise false if A == C. Recall that
    // RobustCCW(x,y,z) == -RobustCCW(z,y,x) for all x,y,z.
    var sum = 0
    if robustCCW(a: b, b: o, c: a) >= 0 {
      sum += 1
    }
    if robustCCW(a: c, b: o, c: b) >= 0 {
      sum += 1
    }
    if robustCCW(a: a, b: o, c: c) > 0 {
      sum += 1
    }
    return sum >= 2
  }
  
  /**
   * Return the angle at the vertex B in the triangle ABC. The return value is
   * always in the range [0, Pi]. The points do not need to be normalized.
   * Ensures that Angle(a,b,c) == Angle(c,b,a) for all a,b,c.
   *
   *  The angle is undefined if A or C is diametrically opposite from B, and
   * becomes numerically unstable as the length of edge AB or BC approaches 180
   * degrees.
   */
  public static func angle(a: S2Point, b: S2Point, c: S2Point) -> Double {
    return a.crossProd(b).angle(to: c.crossProd(b))
  }
  
  /**
   Return the exterior angle at the vertex B in the triangle ABC. The return
   value is positive if ABC is counterclockwise and negative otherwise. If you
   imagine an ant walking from A to B to C, this is the angle that the ant
   turns at vertex B (positive = left, negative = right). Ensures that
   TurnAngle(a,b,c) == -TurnAngle(c,b,a) for all a,b,c.
   
   - Returns: the exterior angle at the vertex B in the triangle ABC
   */
  public static func turnAngle(a: S2Point, b: S2Point, c: S2Point) -> Double {
    // This is a bit less efficient because we compute all 3 cross products, but
    // it ensures that turnAngle(a,b,c) == -turnAngle(c,b,a) for all a,b,c.
    let outAngle = b.crossProd(a).angle(to: c.crossProd(b))
    return (robustCCW(a: a, b: b, c: c) > 0) ? outAngle : -outAngle
  }
  
  /// Return true if two points are within the given distance of each other (mainly useful for testing).
  public static func approxEquals(_ a: S2Point, _ b: S2Point, maxError: Double = 1e-15) -> Bool {
    return a.angle(to: b) <= maxError
  }
  
  /// Return true if two points are within the given distance of each other (mainly useful for testing).
  public static func approxEquals(_ a: Double, _ b: Double, maxError: Double = 1e-15) -> Bool {
    return Swift.abs(a - b) <= maxError
  }
  
}
