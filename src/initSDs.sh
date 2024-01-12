#!/bin/bash
#SBATCH --job-name=initSDs
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --mem=30G
#SBATCH --time=00:30:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh1126
#SBATCH --output=./initSDs_out.%j.out
#SBATCH --error=./initSDs_err.%j.out

### ----- You need to edit these lines to set your ----- ###
### ----- default compiler and python environment   ---- ###
### ----  and paths for CLEO and build directories  ---- ###
module load python3/2022.01-gcc-11.2.0
source activate /work/mh1126/m300950/condaenvs/superdropsenv

path2CLEO=${HOME}/CLEO/
# path2build=/work/mh1126/m300950/droplet_breakup_partii/build/
path2build=${HOME}/breakup_partii/build/
configfile=${HOME}/breakup_partii/src/src/buii_config.txt 

python=/work/mh1126/m300950/condaenvs/superdropsenv/bin/python
### ---------------------------------------------------- ###

### --------------------- initSDs ---------------------- ###
### make initial SD conditions for 1-D rainshaft
${python} initSDs.py ${path2CLEO} ${path2build} ${configfile}
### ---------------------------------------------------- ###