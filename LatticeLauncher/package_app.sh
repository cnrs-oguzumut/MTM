#!/bin/bash
set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
APP_DIR="${DIR}/Lattice Launcher.app"
CONTENTS="${APP_DIR}/Contents"
MACOS="${CONTENTS}/MacOS"
RESOURCES="${CONTENTS}/Resources"

echo "Building LatticeLauncher (Release)..."
cd "$DIR"

mkdir -p "$MACOS"
mkdir -p "$RESOURCES"

XCODE_SWIFTC="/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/swiftc"
XCODE_SDK="/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk"
XCODE_PLUGINS="/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/usr/lib/swift/host/plugins"

if [ -x "$XCODE_SWIFTC" ] && [ -d "$XCODE_PLUGINS" ]; then
    "$XCODE_SWIFTC" -O -sdk "$XCODE_SDK" -plugin-path "$XCODE_PLUGINS" -module-cache-path /tmp/swift-module-cache Sources/*.swift -o "$MACOS/LatticeLauncher"
else
    swift build -c release --cache-path /tmp/swift-cache
    BIN_PATH="$(swift build -c release --show-bin-path)/LatticeLauncher"
    cp "$BIN_PATH" "$MACOS/LatticeLauncher"
fi

# Copy AppIcon if present
if [ -f "${DIR}/AppIcon.icns" ]; then
    cp "${DIR}/AppIcon.icns" "${RESOURCES}/AppIcon.icns"
fi

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
    <key>CFBundleIconFile</key>
    <string>AppIcon</string>
    <key>NSHighResolutionCapable</key>
    <true/>
</dict>
</plist>
EOF

echo "App bundle created: ${APP_DIR}"
