#!/bin/bash
set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${DIR}/.." && pwd)"
APP_DIR="${DIR}/Lattice Launcher.app"
ROOT_APP_DIR="${ROOT_DIR}/Lattice Launcher.app"
CONTENTS="${APP_DIR}/Contents"
MACOS="${CONTENTS}/MacOS"
RESOURCES="${CONTENTS}/Resources"
DOCS_OUT="${RESOURCES}/Docs"
FRAMEWORKS="${CONTENTS}/Frameworks"

# Universal Apple Silicon baseline: M1 / M2 / M3 / M4 / M5+ on macOS 14.0 (Sonoma) or newer
MIN_MACOS="14.0"
SWIFT_TARGET="arm64-apple-macos${MIN_MACOS}"

echo "Building LatticeLauncher GUI (${SWIFT_TARGET}, Release)..."
cd "$DIR"

mkdir -p "$MACOS" "$RESOURCES" "$DOCS_OUT"
rm -rf "$FRAMEWORKS"
mkdir -p "$FRAMEWORKS"

XCODE_SWIFTC="/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/swiftc"
XCODE_SDK="/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk"
XCODE_PLUGINS="/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/usr/lib/swift/host/plugins"

if [ -x "$XCODE_SWIFTC" ] && [ -d "$XCODE_PLUGINS" ]; then
    "$XCODE_SWIFTC" -O -target "$SWIFT_TARGET" -sdk "$XCODE_SDK" -plugin-path "$XCODE_PLUGINS" -module-cache-path /tmp/swift-module-cache Sources/*.swift -o "$MACOS/LatticeLauncher"
else
    swift build -c release --arch arm64 -Xswiftc -target -Xswiftc "$SWIFT_TARGET" --cache-path /tmp/swift-cache
    BIN_PATH="$(swift build -c release --arch arm64 --show-bin-path)/LatticeLauncher"
    cp "$BIN_PATH" "$MACOS/LatticeLauncher"
fi
chmod +x "$MACOS/LatticeLauncher"

# Bundle C++ solver executable and all transitive non-system dylibs for standalone M-series distribution
CPP_BIN="${ROOT_DIR}/build/lattice_triangulation"
if [ ! -f "$CPP_BIN" ]; then
    echo "ERROR: C++ solver binary not found at ${CPP_BIN}. Build CMake target lattice_triangulation first." >&2
    exit 1
fi

echo "Bundling C++ solver executable and dynamic libraries..."
cp "$CPP_BIN" "$MACOS/lattice_triangulation"
chmod +x "$MACOS/lattice_triangulation"

SDK_VER="$(xcrun --show-sdk-version 2>/dev/null || echo "14.0")"

python3 - "$MACOS/lattice_triangulation" "$FRAMEWORKS" "$MIN_MACOS" "$SDK_VER" << 'PYEOF'
import os, shutil, subprocess, sys

exe_path = sys.argv[1]
fw_dir = sys.argv[2]
min_macos = sys.argv[3]
sdk_ver = sys.argv[4]
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
    if dep.startswith("@rpath/") or dep.startswith("@loader_path/") or dep.startswith("@executable_path/"):
        name = os.path.basename(dep)
        cand = os.path.join(loader_dir, name)
        if os.path.exists(cand):
            return os.path.realpath(cand), name
        for sdir in search_dirs:
            cand2 = os.path.join(sdir, name)
            if os.path.exists(cand2):
                return os.path.realpath(cand2), name
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
            subprocess.run(["install_name_tool", "-id", f"@loader_path/{base_name}", target], check=True)
            continue
        if base_name not in visited:
            visited.add(base_name)
            shutil.copy2(real_src, dst)
            os.chmod(dst, 0o755)
            subprocess.run(["install_name_tool", "-id", f"@loader_path/{base_name}", dst], check=True)
            queue.append((dst, os.path.dirname(real_src), False))
        new_ref = f"@executable_path/../Frameworks/{base_name}" if is_exe else f"@loader_path/{base_name}"
        subprocess.run(["install_name_tool", "-change", dep, new_ref, target], check=True)

# Normalize LC_BUILD_VERSION minos on bundled dylibs so macOS 14+ dyld accepts all libraries cleanly
vtool_bin = shutil.which("vtool")
if vtool_bin:
    for base_name in sorted(visited):
        dst = os.path.join(fw_dir, base_name)
        subprocess.run(
            [vtool_bin, "-set-build-version", "macos", min_macos, sdk_ver, "-replace", "-output", dst, dst],
            check=True
        )

# Verify zero external /opt/homebrew, /usr/local, or /Users load paths remain
all_bins = [exe_path] + [os.path.join(fw_dir, b) for b in sorted(visited)]
leaks = []
for b in all_bins:
    for dep in get_deps(b):
        if any(dep.startswith(p) for p in ("/opt/homebrew", "/usr/local", "/Users/")):
            leaks.append((os.path.basename(b), dep))

if leaks:
    print("ERROR: Non-portable dynamic library references remain:", file=sys.stderr)
    for bname, dep in leaks:
        print(f"  {bname} -> {dep}", file=sys.stderr)
    sys.exit(1)

print(f"Bundled and verified {len(visited)} portable arm64 dynamic libraries into Contents/Frameworks (minos {min_macos}).")
PYEOF

# Copy AppIcon if present
if [ -f "${DIR}/AppIcon.icns" ]; then
    cp "${DIR}/AppIcon.icns" "${RESOURCES}/AppIcon.icns"
fi

# Bundle documentation guides (.md and .pdf) into Contents/Resources/Docs
rm -rf "${DOCS_OUT}"
mkdir -p "${DOCS_OUT}"
if [ -d "${ROOT_DIR}/docs" ]; then
    for docfile in "${ROOT_DIR}/docs"/*.md "${ROOT_DIR}/docs"/*.pdf; do
        if [ -f "$docfile" ]; then
            cp "$docfile" "${DOCS_OUT}/"
        fi
    done
fi

# Create Info.plist with explicit macOS 14.0 minimum system version & arm64 priority
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
    <key>CFBundleDisplayName</key>
    <string>Lattice Launcher</string>
    <key>CFBundlePackageType</key>
    <string>APPL</string>
    <key>CFBundleShortVersionString</key>
    <string>1.2</string>
    <key>CFBundleVersion</key>
    <string>3</string>
    <key>CFBundleIconFile</key>
    <string>AppIcon</string>
    <key>LSMinimumSystemVersion</key>
    <string>${MIN_MACOS}</string>
    <key>LSArchitecturePriority</key>
    <array>
        <string>arm64</string>
    </array>
    <key>NSHighResolutionCapable</key>
    <true/>
</dict>
</plist>
EOF

# Strip extended attributes/quarantine and ad-hoc codesign all bundled Mach-O files and app bundle
xattr -cr "${APP_DIR}" || true
for dylib in "${FRAMEWORKS}"/*.dylib; do
    if [ -f "$dylib" ]; then
        codesign --force --sign - "$dylib"
    fi
done
codesign --force --sign - "${MACOS}/lattice_triangulation"
codesign --force --sign - "${MACOS}/LatticeLauncher"
codesign --force --deep --sign - "${APP_DIR}"
codesign --verify --deep --strict "${APP_DIR}"

# Sync root-level Lattice Launcher.app copy if kept in repository root
rm -rf "${ROOT_APP_DIR}"
cp -R "${APP_DIR}" "${ROOT_APP_DIR}"

echo "App bundle created, verified, and signed: ${APP_DIR}"

# Build distributable DMG for Apple Silicon (M1 / M2 / M3 / M4+) Macs
DMG_STAGING="$(mktemp -d /tmp/LatticeLauncher_dmg_staging.XXXXXX)"
DMG_PATH="${ROOT_DIR}/LatticeLauncher.dmg"
rm -f "$DMG_PATH"

cp -R "${APP_DIR}" "${DMG_STAGING}/Lattice Launcher.app"
ln -s /Applications "${DMG_STAGING}/Applications"

if [ -f "${DIR}/AppIcon.icns" ]; then
    cp "${DIR}/AppIcon.icns" "${DMG_STAGING}/.VolumeIcon.icns"
fi

mkdir -p "${DMG_STAGING}/Documentation"
for pdf in \
    "discrete_burgers_circuit_method.pdf" \
    "plasticity_memory_remeshing.pdf" \
    "aligned_deformation_gradient.pdf" \
    "special_grain_boundaries_and_bcs_math.pdf" \
    "stress_and_traction_control.pdf" \
    "preconditioned_lbfgs.pdf"; do
    if [ -f "${ROOT_DIR}/docs/${pdf}" ]; then
        cp "${ROOT_DIR}/docs/${pdf}" "${DMG_STAGING}/Documentation/${pdf}"
    fi
done

cat << 'READMEEOF' > "${DMG_STAGING}/READ_ME_FIRST.txt"
========================================================================
  Lattice Launcher — Mesoscopic Tensorial Model (MTM) Crystal Studio
  Standalone Apple Silicon Edition (M1 / M2 / M3 / M4+ Macs)
========================================================================

SYSTEM REQUIREMENTS
  • Any Mac with Apple Silicon (M1, M2, M3, M4, or newer)
  • macOS 14.0 (Sonoma), macOS 15 (Sequoia), macOS 26+ (Tahoe), or newer
  • Zero external dependencies required — the C++ solver and all math/BLAS
    libraries (OpenMP, SuiteSparse/CHOLMOD, Armadillo, OpenBLAS, ARPACK,
    GMP/MPFR) are bundled directly inside Lattice Launcher.app.

INSTALLATION
  1. Drag "Lattice Launcher.app" into the "Applications" folder shortcut.
  2. Open "Lattice Launcher" from /Applications.
  3. By default, standalone simulation runs and VTK/CSV outputs are saved
     in: ~/Documents/LatticeRuns

FIRST-TIME OPEN NOTE (macOS Gatekeeper)
  Because this research distribution is ad-hoc signed rather than Apple
  Developer ID notarized, macOS may show a prompt on the very first launch
  when downloaded from a browser or email ("Apple could not verify...").

  To allow opening on first launch, use EITHER option once:

  Option A — System Settings (No Terminal needed):
    1. Try to open Lattice Launcher once, then click "Done".
    2. Open System Settings -> Privacy & Security.
    3. Scroll down to the Security section and click "Open Anyway"
       next to "Lattice Launcher".

  Option B — One-Line Terminal Command:
    xattr -cr "/Applications/Lattice Launcher.app"

========================================================================
READMEEOF

xattr -cr "${DMG_STAGING}" || true

SETFILE_BIN="$(xcrun -f SetFile 2>/dev/null || true)"
if [ -n "$SETFILE_BIN" ] && [ -f "${DMG_STAGING}/.VolumeIcon.icns" ]; then
    "$SETFILE_BIN" -c icnC "${DMG_STAGING}/.VolumeIcon.icns" 2>/dev/null || true
    "$SETFILE_BIN" -a C "${DMG_STAGING}" 2>/dev/null || true
fi

if command -v hdiutil >/dev/null 2>&1; then
    hdiutil create -volname "Lattice Launcher" -srcfolder "$DMG_STAGING" -ov -format UDZO -imagekey zlib-level=9 "$DMG_PATH"
    xattr -cr "$DMG_PATH" || true
    echo "Created standalone Apple Silicon DMG installer: ${DMG_PATH}"
fi
rm -rf "$DMG_STAGING"
