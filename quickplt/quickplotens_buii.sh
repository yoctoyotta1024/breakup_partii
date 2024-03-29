#!/bin/bash
#SBATCH --job-name=plotens
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --mem=30G
#SBATCH --time=00:30:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh1126
#SBATCH --output=./plotens_out.%j.out
#SBATCH --error=./plotens_err.%j.out

### ----- You need to edit these lines to set your ----- ###
### ----- default compiler and python environment   ---- ###
### ----  and paths for CLEO and build directories  ---- ###
module load python3/2022.01-gcc-11.2.0
source activate /work/mh1126/m300950/condaenvs/superdropsenv

path2CLEO=${HOME}/CLEO/
path2build=/work/mh1126/m300950/droplet_breakup_partii/build/

datapath=${path2build}bin/${1}/ensemb/

python=/work/mh1126/m300950/condaenvs/superdropsenv/bin/python
### ---------------------------------------------------- ###

### ----------------------- plots ------------------------ ###
### run 1-D rainshaft quick plotting script
${python} quickplotens_buii.py ${path2CLEO} ${path2build} ${datapath}
### ---------------------------------------------------- ###