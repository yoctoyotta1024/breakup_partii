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

wetdry=${1}

path2CLEO=${HOME}/CLEO/
path2build=/work/mh1126/m300950/droplet_breakup_partii/build/
tmp_configfile=${path2build}/tmp/buii_config.txt 

python=/work/mh1126/m300950/condaenvs/superdropsenv/bin/python
### ---------------------------------------------------- ###

### --------------------- initSDs ---------------------- ###
### make initial SD conditions for 1-D rainshaft
if [ "${wetdry}" != "dry" ] && [ "${wetdry}" != "wet" ];
then
  echo "please specify 'wet' or 'dry' initial superdroplets"

elif [ "${wetdry}" == "dry" ];
then
  ${python} initSDs_dry.py ${path2CLEO} ${path2build} ${tmp_configfile}

elif [ "${wetdry}" == "wet" ];
then
  ${python} initSDs_wet.py ${path2CLEO} ${path2build} ${tmp_configfile}
fi
### ---------------------------------------------------- ###