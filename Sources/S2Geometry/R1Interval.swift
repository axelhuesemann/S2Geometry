//
//  R1Interval.swift
//  S2Geometry
//
//  Created by Alex Studnicka on 6/30/16.
//  Copyright © 2016 Alex Studnicka. MIT License.
//

/// An R1Interval represents a closed, bounded interval on the real line. It is
/// capable of representing the empty interval (containing no points) and
/// zero-length intervals (containing a single point).
public struct R1Interval {
	
	public let lo: Double
	public let hi: Double
	
  // MARK: init / factory
  
	/// Interval constructor. If lo > hi, the interval is empty.
	public init(lo: Double, hi: Double) {
		self.lo = lo
		self.hi = hi
	}
	
	/// Convenience method to construct an interval containing a single point.
	public init(point p: Double) {
		self.init(lo: p, hi: p)
	}
	
	/// Convenience method to construct the minimal interval containing the two
  /// given points. This is equivalent to starting with an empty interval and
  /// calling AddPoint() twice, but it is more efficient.
	public init(p1: Double, p2: Double) {
		if p1 <= p2 {
			self.init(lo: p1, hi: p2)
		} else {
			self.init(lo: p2, hi: p1)
		}
	}
	
  /// Returns an empty interval. (Any interval where lo > hi is considered empty.)
  public static let empty = R1Interval(lo: 1, hi: 0)
  
  // MARK: tests
  
	/// Return true if the interval is empty, i.e. it contains no points.
	public var isEmpty: Bool {
		return lo > hi
	}
	
  // MARK: computed members
  
	/// Return the center of the interval. For empty intervals, the result is arbitrary.
	public var center: Double {
		return 0.5 * (lo + hi)
	}
	
	/// Return the length of the interval. The length of an empty interval is negative.
	public var length: Double {
		return hi - lo
	}

  // MARK: contains / intersects
  
	public func contains(point p: Double) -> Bool {
		return p >= lo && p <= hi
	}
	
	public func interiorContains(point p: Double) -> Bool {
		return p > lo && p < hi
	}
	
	/// Return true if this interval contains the interval 'y'.
	public func contains(interval y: R1Interval) -> Bool {
		guard !y.isEmpty else { return true }
		return y.lo >= lo && y.hi <= hi
	}
	
	/// Return true if the interior of this interval contains the entire interval 'y' (including its boundary).
	public func interiorContains(interval y: R1Interval) -> Bool {
		guard !y.isEmpty else { return true }
		return y.lo > lo && y.hi < hi
	}
	
	/// Return true if this interval intersects the given interval, i.e. if they have any points in common.
	public func intersects(interval y: R1Interval) -> Bool {
		if lo <= y.lo {
			return y.lo <= hi && y.lo <= y.hi
		} else {
			return lo <= y.hi && lo <= hi
		}
	}
	
	/// Return true if the interior of this interval intersects any point of the given interval (including its boundary).
	public func interiorIntersects(interval y: R1Interval) -> Bool {
		return y.lo < hi && lo < y.hi && lo < hi && y.lo <= y.hi
	}
	
  // MARK: arithmetic
  
	/// Expand the interval so that it contains the given point "p".
	public func add(point p: Double) -> R1Interval {
		if isEmpty {
			return R1Interval(point: p)
		} else if p < lo {
			return R1Interval(lo: p, hi: hi)
		} else if p > hi {
			return R1Interval(lo: lo, hi: p)
		} else {
			return R1Interval(lo: lo, hi: hi)
		}
	}
	
	/// Return an interval that contains all points with a distance "radius" of a
  /// point in this interval. Note that the expansion of an empty interval is always empty.
	public func expanded(radius: Double) -> R1Interval {
		guard !isEmpty else { return self }
		return R1Interval(lo: lo - radius, hi: hi + radius)
	}
	
	/// Return the smallest interval that contains this interval and the given interval "y".
	public func union(interval y: R1Interval) -> R1Interval {
		guard !isEmpty else { return y }
		guard !y.isEmpty else { return self }
		return R1Interval(lo: min(lo, y.lo), hi: max(hi, y.hi))
	}
	
	/// Return the intersection of this interval with the given interval. Empty intervals do not need to be special-cased.
	public func intersection(interval y: R1Interval) -> R1Interval {
		return R1Interval(lo: max(lo, y.lo), hi: min(hi, y.hi))
	}
	
  // MARK: compare
  
	/// Return true if length of the symmetric difference between the two intervals is at most the given tolerance.
	public func approxEquals(_ other: R1Interval, maxError: Double = 1e-15) -> Bool {
		guard !isEmpty else { return other.length <= maxError }
		guard !other.isEmpty else { return length <= maxError }
		return abs(other.lo - lo) + abs(other.hi - hi) <= maxError
	}
	
}

extension R1Interval: Equatable, Hashable {

  public static func ==(lhs: R1Interval, rhs: R1Interval) -> Bool {
    return (lhs.lo == rhs.lo && lhs.hi == rhs.hi) || (lhs.isEmpty && rhs.isEmpty)
  }

  public var hashValue: Int {
    guard !isEmpty else { return 17 }
    var value: Int64 = 17
    value = 37 * value + Int64(bitPattern: lo.bitPattern)
    value = 37 * value + Int64(bitPattern: hi.bitPattern)
    return Int(value ^ (value >> 32))
  }
  
}
