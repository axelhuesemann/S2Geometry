// swift-tools-version:4.0

import PackageDescription

let package = Package(
  name: "S2Geometry",
  dependencies: [],
  targets: [
    .target(
      name: "S2Geometry",
      dependencies: []
    ),
    .testTarget(name: "S2GeometryTests", dependencies: ["S2Geometry"]),
  ]
)
