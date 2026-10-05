#!/bin/bash
set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
APP_DIR="${DIR}/Lattice Launcher.app"
CONTENTS="${APP_DIR}/Contents"
MACOS="${CONTENTS}/MacOS"
RESOURCES="${CONTENTS}/Resources"

echo "Building LatticeLauncher (Release)..."
cd "$DIR"
swift build -c release

BIN_PATH="$(swift build -c release --show-bin-path)/LatticeLauncher"

mkdir -p "$MACOS"
mkdir -p "$RESOURCES"
cp "$BIN_PATH" "$MACOS/LatticeLauncher"

# Create Info.plist
cat <<EOF > "${CONTENTS}/Info.plist"
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
    <key>CFBundleExecutable</key>
    <string>LatticeLauncher</string>
    <key>CFBundleIdentifier</key>
    <string>com.usalman.LatticeLauncher</string>
    <key>CFBundleName</key>
    <string>Lattice Launcher</string>
    <key>CFBundlePackageType</key>
    <string>APPL</string>
    <key>CFBundleShortVersionString</key>
    <string>1.0</string>
    <key>CFBundleVersion</key>
    <string>1</string>
    <key>NSHighResolutionCapable</key>
    <true/>
</dict>
</plist>
EOF

echo "App bundle created: ${APP_DIR}"
