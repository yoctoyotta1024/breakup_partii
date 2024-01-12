#!/bin/bash
#SBATCH --job-name=runbuii
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --mem=30G
#SBATCH --time=00:30:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh1126
#SBATCH --output=./runbuii_out.%j.out
#SBATCH --error=./runbuii_err.%j.out

### ----- You need to edit these lines to set your ----- ###
### ----- default compiler and python environment   ---- ###
### ----  and paths for CLEO and build directories  ---- ###
module load python3/2022.01-gcc-11.2.0
source activate /work/mh1126/m300950/condaenvs/superdropsenv

executable=buii_${1}

# path2build=/work/mh1126/m300950/droplet_breakup_partii/build/
path2build=${HOME}/breakup_partii/build/
configfile=${HOME}/breakup_partii/src/src/buii_config.txt 

python=/work/mh1126/m300950/condaenvs/superdropsenv/bin/python
### ---------------------------------------------------- ###

### ------------------- prepare for run ---------------- ###
export OMP_PROC_BIND=spread
export OMP_PLACES=threads

### ensure these directories exist (it's a good idea for later use)
mkdir ${path2build}bin
mkdir ${path2build}share
### ---------------------------------------------------- ###

### ----------------------- run ------------------------ ###
### run 1-D rainshaft 
${python} run_buii.py ${path2build} ${configfile} ${executable}
### ---------------------------------------------------- ###