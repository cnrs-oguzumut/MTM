#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

echo "Checking for disallowed direct lagrange::reduce / reduceAndCount calls..."

DISALLOWED_REDUCE=$(grep -rn "lagrange::reduce" src/ include/ \
    | grep -v "^include/reductions/LagrangeReduction\.h:" \
    | grep -v "^include/lattice_energy/ReferenceFrame\.h:" \
    | grep -v "^src/lagrange_reduction\.cpp:" \
    | grep -v "^src/FEMHessianAssembler\.cpp:" \
    | grep -v "^src/energy_spline_interpolator\.cpp:" \
    | grep -v "^src/experiments/acoustic_studies\.cpp:" \
    | grep -v ":[0-9]*:[[:space:]]*//" \
    | grep -v ":[0-9]*:[[:space:]]*\*" \
    || true)

if [[ -n "$DISALLOWED_REDUCE" ]]; then
    echo "ERROR: Found direct lagrange::reduce calls outside allow-listed files:"
    echo "$DISALLOWED_REDUCE"
    exit 1
fi

echo "Checking for disallowed local p2(h / calculateShapeDerivatives calls in active drivers..."

DISALLOWED_P2H=$(grep -rn "p2(h" src/ include/ \
    | grep -v "^src/ReferenceFrame\.cpp:" \
    | grep -v "^src/experiments/dislocation_study\.cpp:40:" \
    | grep -v "^src/experiments/zanzotto_examples\.cpp:31:" \
    | grep -v "^src/experiments/data_analysis\.cpp:" \
    | grep -v "^src/experiments/shifted_crystal_study\.cpp:" \
    | grep -v "^src/experiments/stress_controlled_examples\.cpp:" \
    | grep -v "^src/experiments/shifting_examples\.cpp:" \
    | grep -v "^src/experiments/shifted_dislocation_study\.cpp:" \
    | grep -v "^src/experiments/dislocation_indentation\.cpp:" \
    | grep -v ":[0-9]*:[[:space:]]*//" \
    || true)

if [[ -n "$DISALLOWED_P2H" ]]; then
    echo "ERROR: Found local p2(h...) shape derivative setup outside ReferenceFrame.cpp:"
    echo "$DISALLOWED_P2H"
    exit 1
fi

DISALLOWED_DNDX=$(grep -rn "calculateShapeDerivatives(" src/ include/ \
    | grep -v "^include/experiments/common_simulation_helpers\.h:" \
    | grep -v "^src/experiments/common_simulation_helpers\.cpp:" \
    | grep -v "^src/ReferenceFrame\.cpp:" \
    | grep -v "^src/experiments/dislocation_study\.cpp:43:" \
    | grep -v "^src/experiments/zanzotto_examples\.cpp:34:" \
    | grep -v "^src/experiments/data_analysis\.cpp:" \
    | grep -v "^src/experiments/shifted_crystal_study\.cpp:" \
    | grep -v "^src/experiments/stress_controlled_examples\.cpp:" \
    | grep -v "^src/experiments/shifting_examples\.cpp:" \
    | grep -v "^src/experiments/shifted_dislocation_study\.cpp:" \
    | grep -v "^src/experiments/dislocation_indentation\.cpp:" \
    | grep -v ":[0-9]*:[[:space:]]*//" \
    || true)

if [[ -n "$DISALLOWED_DNDX" ]]; then
    echo "ERROR: Found local calculateShapeDerivatives calls outside make_reference_setup:"
    echo "$DISALLOWED_DNDX"
    exit 1
fi

echo "OK: All reduction and dndx setup sites are properly routed through ReferenceFrame."
