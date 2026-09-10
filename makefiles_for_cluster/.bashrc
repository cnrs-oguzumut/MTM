# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
    . /etc/bashrc
fi

if [ -f /etc/profile.d/modules.sh ]; then
    export MODULES_AUTO_HANDLING=1
    source /etc/profile.d/modules.sh
    # Load SLURM (new Magi environment Debian Trixie / fallback to legacy)
    if [ -d /softs/trixie/modules/x86-64/slurm ] || [ -f /softs/trixie/modules/x86-64/slurm/default ]; then
        module load /softs/trixie/modules/x86-64/slurm/default
    else
        module load /softs/modules/slurm/default 2>/dev/null || module load slurm 2>/dev/null
    fi
    module use /softs/modules/oneapi/2024.1.0.596
    module load compiler/2024.1.0
    module load mkl/2024.1
    module load mpi/2021.12
fi

# --- Library paths (append, don't overwrite!) ---
export LD_LIBRARY_PATH=/softs/armadillo-12.8.2/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/softs/boost_1_85_0/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/softs/superlu-5.2.2/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/dist/umut.salman/codes/fadeRelease_v2.03/lib_ubuntu16.04_x86_64:$LD_LIBRARY_PATH

# Intel runtime (for libiomp5.so and MKL)
export LD_LIBRARY_PATH=/softs/oneapi/2024.1.0.596/compiler/2024.1/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/softs/oneapi/2024.1.0.596/mkl/2024.1/lib/intel64:$LD_LIBRARY_PATH

# --- Other paths ---
export LIBRARY_PATH=/home/dist/umut.salman/codes/armadillo-11.1.1/lib:/softs/superlu-5.2.2/lib:$LIBRARY_PATH
export C_INCLUDE_PATH=/home/dist/umut.salman/codes/armadillo-11.1.1/include:/softs/superlu-5.2.2/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=/home/dist/umut.salman/codes/armadillo-11.1.1/include:/softs/superlu-5.2.2/include:$CPLUS_INCLUDE_PATH
export CPATH=/home/dist/umut.salman/codes/armadillo-11.1.1/include:/softs/superlu-5.2.2/include:$CPATH
export OBJC_INCLUDE_PATH=/home/dist/umut.salman/codes/armadillo-11.1.1/include:/softs/superlu-5.2.2/include:$OBJC_INCLUDE_PATH

# Python virtual environment
source ~/venv/bin/activate
export PATH=$HOME/.local/bin:$PATH
export PATH=/home/dist/umut.salman/.local/bin:$PATH
export LD_LIBRARY_PATH=/home/dist/umut.salman/.local/lib:$LD_LIBRARY_PATH
