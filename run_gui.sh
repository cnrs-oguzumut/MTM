#!/bin/bash
# Launches the Swift GUI for Lattice Triangulation on Mac
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR/LatticeLauncher"

# Build and run with Swift
swift run LatticeLauncher
