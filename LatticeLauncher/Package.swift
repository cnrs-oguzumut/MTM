// swift-tools-version: 5.9
import PackageDescription

let package = Package(
    name: "LatticeLauncher",
    platforms: [
        .macOS(.v14)
    ],
    targets: [
        .executableTarget(
            name: "LatticeLauncher",
            path: "Sources"
        )
    ]
)
