#!/bin/bash
#SBATCH --job-name=setupbuii
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --mem=30G
#SBATCH --time=00:30:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh1126
#SBATCH --output=./setupbuii_out.%j.out
#SBATCH --error=./setupbuii_err.%j.out

### ----- You need to edit these lines to set your ----- ###
### ----- default compiler and python environment   ---- ###
### ----  and paths for CLEO and build directories  ---- ###
module load gcc/11.2.0-gcc-11.2.0
module load python3/2022.01-gcc-11.2.0
module load nvhpc/23.7-gcc-11.2.0
spack load cmake@3.23.1%gcc
source activate /work/mh1126/m300950/condaenvs/superdropsenv

path2CLEO=${HOME}/CLEO/
path2src=${HOME}/breakup_partii/src/
path2build=/work/mh1126/m300950/droplet_breakup_partii/build/
orig_configfile=${HOME}/breakup_partii/src/src/buii_config.txt 
tmp_configfile=${path2build}/tmp/buii_config.txt 

python=/work/mh1126/m300950/condaenvs/superdropsenv/bin/python
gxx="g++"
gcc="gcc"
cuda="nvc++"
### ---------------------------------------------------- ###

### ------------ choose Kokkos configuration ----------- ###
kokkosflags="-DKokkos_ARCH_NATIVE=ON -DKokkos_ARCH_AMPERE80=ON -DKokkos_ENABLE_SERIAL=ON" # serial kokkos
kokkoshost="-DKokkos_ENABLE_OPENMP=ON"                                          # flags for host parallelism (e.g. using OpenMP)
kokkosdevice=""           # flags for device parallelism (e.g. using CUDA)
### ---------------------------------------------------- ###

### ------------------------ build --------------------- ###
### build CLEO using cmake (with openMP thread parallelism through Kokkos)
buildcmd="CXX=${gxx} CC=${gcc} CUDA=${cuda} cmake -S ${path2src} -B ${path2build} ${kokkosflags} ${kokkoshost} ${kokkosdevice}"
echo ${buildcmd}
CXX=${gxx} CC=${gcc} CUDA=${cuda} cmake -S ${path2src} -B ${path2build} ${kokkosflags} ${kokkoshost} ${kokkosdevice} -DCLEOLIBS_SOURCE_DIR:STRING=${path2CLEO}libs

export OMP_PROC_BIND=spread
export OMP_PLACES=threads
### ---------------------------------------------------- ###

### --------------------- compile ---------------------- ###
### generate input files and compile 1-D rainshaft for list of executables
${python} build_buii.py ${path2CLEO} ${path2build} ${orig_configfile} ${tmp_configfile} 
### ---------------------------------------------------- ###