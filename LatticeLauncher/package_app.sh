#!/bin/bash
set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${DIR}/.." && pwd)"
APP_DIR="${DIR}/Lattice Launcher.app"
CONTENTS="${APP_DIR}/Contents"
MACOS="${CONTENTS}/MacOS"
RESOURCES="${CONTENTS}/Resources"
FRAMEWORKS="${CONTENTS}/Frameworks"

echo "Building LatticeLauncher (Release)..."
cd "$DIR"

mkdir -p "$MACOS"
mkdir -p "$RESOURCES"
mkdir -p "$FRAMEWORKS"

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

# Bundle C++ solver executable and all non-system dylibs for standalone distribution
CPP_BIN="${ROOT_DIR}/build/lattice_triangulation"
if [ -f "$CPP_BIN" ]; then
    echo "Bundling C++ solver executable and dynamic libraries..."
    cp "$CPP_BIN" "$MACOS/lattice_triangulation"
    chmod +x "$MACOS/lattice_triangulation"

    python3 - "$MACOS/lattice_triangulation" "$FRAMEWORKS" << 'PYEOF'
import os, shutil, subprocess, sys

exe_path = sys.argv[1]
fw_dir = sys.argv[2]
os.makedirs(fw_dir, exist_ok=True)

search_dirs = [
    "/opt/homebrew/lib",
    "/opt/homebrew/opt/gcc/lib/gcc/current",
    "/opt/homebrew/opt/libomp/lib",
    "/opt/homebrew/opt/suite-sparse/lib",
    "/opt/homebrew/opt/armadillo/lib",
    "/opt/homebrew/opt/arpack/lib",
    "/opt/homebrew/opt/openblas/lib",
    "/opt/homebrew/opt/metis/lib",
    "/opt/homebrew/opt/gmp/lib",
    "/opt/homebrew/opt/mpfr/lib",
]

def get_deps(binary):
    out = subprocess.check_output(["otool", "-L", binary], text=True)
    lines = out.strip().splitlines()[1:]
    deps = []
    for line in lines:
        dep = line.strip().split(" (compatibility")[0].strip()
        deps.append(dep)
    return deps

def resolve_dep(dep, loader_dir):
    if dep.startswith("/System/") or dep.startswith("/usr/lib/"):
        return None
    if dep.startswith("/"):
        if os.path.exists(dep):
            return os.path.realpath(dep), os.path.basename(dep)
    if dep.startswith("@rpath/") or dep.startswith("@loader_path/"):
        name = dep.split("/", 1)[1]
        cand = os.path.join(loader_dir, name)
        if os.path.exists(cand):
            return os.path.realpath(cand), os.path.basename(name)
        for sdir in search_dirs:
            cand2 = os.path.join(sdir, os.path.basename(name))
            if os.path.exists(cand2):
                return os.path.realpath(cand2), os.path.basename(name)
    return None

queue = [(exe_path, os.path.dirname(exe_path), True)]
visited = set()

while queue:
    target, loader_dir, is_exe = queue.pop(0)
    deps = get_deps(target)
    for dep in deps:
        res = resolve_dep(dep, loader_dir)
        if not res:
            continue
        real_src, base_name = res
        dst = os.path.join(fw_dir, base_name)
        if not is_exe and os.path.basename(target) == base_name:
            subprocess.run(["install_name_tool", "-id", f"@loader_path/{base_name}", target], check=False)
            continue
        if base_name not in visited:
            visited.add(base_name)
            shutil.copy2(real_src, dst)
            os.chmod(dst, 0o755)
            subprocess.run(["install_name_tool", "-id", f"@loader_path/{base_name}", dst], check=False)
            queue.append((dst, os.path.dirname(real_src), False))
        new_ref = f"@executable_path/../Frameworks/{base_name}" if is_exe else f"@loader_path/{base_name}"
        subprocess.run(["install_name_tool", "-change", dep, new_ref, target], check=False)

print(f"Bundled {len(visited)} dynamic libraries into Contents/Frameworks.")
PYEOF
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

# Ad-hoc codesign Frameworks and App Bundle so macOS Gatekeeper / dyld accepts rewritten load commands
for dylib in "${FRAMEWORKS}"/*.dylib; do
    if [ -f "$dylib" ]; then
        codesign --force --sign - "$dylib"
    fi
done
if [ -f "${MACOS}/lattice_triangulation" ]; then
    codesign --force --sign - "${MACOS}/lattice_triangulation"
fi
codesign --force --deep --sign - "${APP_DIR}"

echo "App bundle created and signed: ${APP_DIR}"

# Build distributable DMG
DMG_STAGING="/tmp/LatticeLauncher_dmg_staging"
DMG_PATH="${ROOT_DIR}/LatticeLauncher.dmg"
rm -rf "$DMG_STAGING" "$DMG_PATH"
mkdir -p "$DMG_STAGING"
cp -R "${APP_DIR}" "${DMG_STAGING}/Lattice Launcher.app"
ln -s /Applications "${DMG_STAGING}/Applications"
if [ -f "${ROOT_DIR}/docs/discrete_burgers_circuit_method.pdf" ]; then
    cp "${ROOT_DIR}/docs/discrete_burgers_circuit_method.pdf" "${DMG_STAGING}/Discrete_Burgers_Circuit_Guide.pdf"
fi

if command -v hdiutil >/dev/null 2>&1; then
    hdiutil create -volname "Lattice Launcher" -srcfolder "$DMG_STAGING" -ov -format UDZO "$DMG_PATH" || true
    if [ -f "$DMG_PATH" ]; then
        echo "Created standalone DMG installer: ${DMG_PATH}"
    fi
fi
rm -rf "$DMG_STAGING"
