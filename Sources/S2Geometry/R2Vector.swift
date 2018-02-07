//
//  R2Vector.swift
//  S2Geometry
//
//  Created by Alex Studnicka on 6/30/16.
//  Copyright © 2016 Alex Studnicka. MIT License.
//

/// Represents a vector in the two-dimensional space. It defines the
/// basic geometrical operations for 2D vectors, e.g. cross product, addition, norm, comparison etc.
public struct R2Vector {
	
	public let x: Double
	public let y: Double
	
  // MARK: init
  
	public init(x: Double = 0, y: Double = 0) {
		self.x = x
		self.y = y
	}
	
  // MARK: computed members
  
	public func get(index: Int) -> Double {
		return index == 0 ? x : y
	}
	
	public var norm2: Double {
		return (x * x) + (y * y)
	}
	
  // MARK: arithmetic
	
	public func dotProd(_ b: R2Vector) -> Double {
		return x * b.x + y * b.y
	}
	
	public func crossProd(_ b: R2Vector) -> Double {
		return x * b.y - y * b.x
	}
	
  // MARK: ccw
  
  public static func planarCCW(a: R2Vector, b: R2Vector) -> Int {
    // Return +1 if the edge AB is CCW around the origin, etc.
    let sab: Double = (a.dotProd(b) > 0) ? -1 : 1
    let vab = a + (b * sab)
    let da = a.norm2
    let db = b.norm2
    let sign: Double
    if da < db || (da == db && a < b) {
      sign = a.crossProd(vab) * sab
    } else {
      sign = vab.crossProd(b)
    }
    if sign > 0 {
      return 1
    }
    if sign < 0 {
      return -1
    }
    return 0
  }
  
  public static func planarOrderedCCW(a: R2Vector, b: R2Vector, c: R2Vector) -> Int {
    let sum = planarCCW(a: a, b: b) + planarCCW(a: b, b: c) + planarCCW(a: c, b: a)
    if sum > 0 {
      return 1
    }
    if sum < 0 {
      return -1
    }
    return 0
  }

}

extension R2Vector: Equatable, Comparable {
  
  public static func ==(lhs: R2Vector, rhs: R2Vector) -> Bool {
    return lhs.x == rhs.x && lhs.y == rhs.y
  }
  public static func <(lhs: R2Vector, rhs: R2Vector) -> Bool {
    if lhs.x < rhs.x {
      return true
    }
    if rhs.x < lhs.x {
      return false
    }
    if lhs.y < rhs.y {
      return true
    }
    return false
  }
  
  public static func +(lhs: R2Vector, rhs: R2Vector) -> R2Vector {
    return R2Vector(x: lhs.x + rhs.x, y: lhs.y + rhs.y)
  }
  
  public static func -(lhs: R2Vector, rhs: R2Vector) -> R2Vector {
    return R2Vector(x: lhs.x - rhs.x, y: lhs.y - rhs.y)
  }
  
  public static func *(lhs: R2Vector, m: Double) -> R2Vector {
    return R2Vector(x: lhs.x * m, y: lhs.y * m)
  }
  
  public static func ⋅(lhs: R2Vector, rhs: R2Vector) -> Double {
    return lhs.dotProd(rhs)
  }

  public static func ×(lhs: R2Vector, rhs: R2Vector) -> Double {
    return lhs.crossProd(rhs)
  }
  
}

precedencegroup DotProductPrecedence {
	associativity: left
	higherThan: MultiplicationPrecedence
}
infix operator ⋅: DotProductPrecedence

precedencegroup CrossProductPrecedence {
	associativity: left
	higherThan: DotProductPrecedence
}
infix operator ×: CrossProductPrecedence
