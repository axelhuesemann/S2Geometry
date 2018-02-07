//
//  S2Metric.swift
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

/// Defines an area or a length cell metric.
public struct S2Metric {
  
  /// The deriv value of a metric is a derivative, and must be multiplied by a length or area in (s,t)-space to get a useful value.
  public let deriv: Double
  public let dim: Int
  
  /// Defines a cell metric of the given dimension (1 == length, 2 == area).
  public init(dim: Int, deriv: Double) {
    self.deriv = deriv
    self.dim = dim
  }
  
  /// Return the value of a metric for cells at the given level.
  public func getValue(level: Int) -> Double {
    return scalb(deriv, Double(dim) * (1 - Double(level)))
  }
  
  /// Return the level at which the metric has approximately the given value.
  /// For example, S2::kAvgEdge.GetClosestLevel(0.1) returns the level at which
  /// the average cell edge length is approximately 0.1. The return value is
  /// always a valid level.
  public func getClosestLevel(value: Double) -> Int {
    return getMinLevel(value: 2.squareRoot() * value)
  }
  
  /// Return the minimum level such that the metric is at most the given value,
  /// or S2CellId::kMaxLevel if there is no such level. For example,
  /// S2::kMaxDiag.GetMinLevel(0.1) returns the minimum level such that all
  /// cell diagonal lengths are 0.1 or smaller. The return value is always a
  /// valid level.
  public func getMinLevel(value: Double) -> Int {
    guard value > 0 else { return S2CellId.maxLevel }
    // This code is equivalent to computing a floating-point "level" value and rounding up.
    #if os(Linux)
      let exponent = Glibc.exp(value / (Double(1 << dim) * deriv))
    #else
      let exponent = Darwin.exp(value / (Double(1 << dim) * deriv))
    #endif
    let level = max(0, min(S2CellId.maxLevel, -((Int(exponent.bitPattern) - 1) >> (dim - 1))))
    // assert (level == S2CellId.MAX_LEVEL || getValue(level) <= value);
    // assert (level == 0 || getValue(level - 1) > value);
    return level
  }
  
  /// Return the maximum level such that the metric is at least the given
  /// value, or zero if there is no such level. For example,
  /// S2.kMinWidth.GetMaxLevel(0.1) returns the maximum level such that all
  /// cells have a minimum width of 0.1 or larger. The return value is always a
  /// valid level.
  public func getMaxLevel(value: Double) -> Int {
    guard value > 0 else { return S2CellId.maxLevel }
    // This code is equivalent to computing a floating-point "level" value and rounding down.
    #if os(Linux)
      let exponent = Glibc.exp(Double(1 << dim) * deriv / value)
    #else
      let exponent = Darwin.exp(Double(1 << dim) * deriv / value)
    #endif
    let level = max(0, min(S2CellId.maxLevel, (Int(exponent.bitPattern - 1) >> (dim - 1))))
    // assert (level == 0 || getValue(level) >= value);
    // assert (level == S2CellId.MAX_LEVEL || getValue(level + 1) < value);
    return level
  }
  
}
