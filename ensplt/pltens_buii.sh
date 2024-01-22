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
savefigpath="/home/m/m300950/breakup_partii/plots/"

python=/work/mh1126/m300950/condaenvs/superdropsenv/bin/python

what2plot=${1}
### ---------------------------------------------------- ###

### ----------------------- plots ------------------------ ###
if [ "${what2plot}" != "massmoms" ] && [ "${what2plot}" != "dists" ] && [ "${what2plot}" != "probs" ];
then
  echo "please specify 'massmoms', 'dists' or 'probs' for plotting"
else
  ${python} pltens_buii.py ${path2CLEO} ${path2build} ${savefigpath} ${what2plot}
fi
### ---------------------------------------------------- ###