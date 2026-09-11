#!/bin/bash
# ==================================================================
# Module loader script for Tchenla cluster (Paris 13 University)
# Usage:
#   source load_modules_tchenla.sh
# ==================================================================

echo "Purging existing modules..."
module purge

echo "Loading Tchenla OpenHPC modules..."
# 1. Compiler + MPI
module load gnu14/14.2.0
module load openmpi5/5.0.7

# 2. Linear Algebra & Math stack
module load openblas/0.3.29
module load boost/1.88.0
# CGAL (Header-only. Prefer local ~/required_libraries/cgal if downloaded)
if [ -d "$HOME/required_libraries/cgal/include" ]; then
    export CGAL_DIR="$HOME/required_libraries/cgal"
elif compgen -G "$HOME/required_libraries/CGAL-*/include" > /dev/null; then
    export CGAL_DIR=$(ls -d $HOME/required_libraries/CGAL-* | head -n 1)
else
    export CGAL_DIR=/opt/software/libs/cgal/gnu14/5.6.2
fi
export CPATH=$CGAL_DIR/include:${CPATH}
module load armadillo/14.4.0
module load eigen/3.4.0

# 3. Python stack for analysis & VTK scripts (optional)
module load py3-numpy/1.26.4 2>/dev/null || true
module load py3-scipy/1.5.4 2>/dev/null || true

echo ""
echo "=========================================================="
echo "Loaded modules on Tchenla:"
module list
echo "=========================================================="
echo "To build, run:"
echo "  make -f Makefile.tchenla -j 8"
echo "=========================================================="
