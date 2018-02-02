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
	
	public init(x: Double = 0, y: Double = 0) {
		self.x = x
		self.y = y
	}
	
	public func get(index: Int) -> Double {
		return index == 0 ? x : y
	}
	
	public var norm2: Double {
		return (x * x) + (y * y)
	}
	
  // MARK:
	
	public func dotProd(_ b: R2Vector) -> Double {
		return x * b.x + y * b.y
	}
	
	public func crossProd(_ b: R2Vector) -> Double {
		return x * b.y - y * b.x
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
