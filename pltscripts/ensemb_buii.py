'''
----- CLEO -----
File: ensemb_buii.py
Project: breakup_partii
Created Date: Friday 17th November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Thursday 11th January 2024
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
Script writes mean and stddev of 
ensemble of datasets from
buii 1-D rainshaft output
'''

import sys
from pathlib import Path

path2CLEO = sys.argv[1]
path2build = sys.argv[2]

sys.path.append(path2CLEO)  # for imports from pySD package
from pySD.sdmout_src.ensembzarr import write_ensemble

### ---------------------------------------------------------------- ###
### ----------------------- INPUT PARAMETERS ----------------------- ###
### ---------------------------------------------------------------- ###
# label and path for each ensembles of datasets 
labels = ["coalbure", "coalonly", "coalbreakup"]
binpath = path2build+"/bin/"      # path before directory called "label" containing zarr datasets
ensembzarr = "sol_ensemb.zarr"      # name of ensemble dataset
ensembsetuptxt = "setup_ensemb.txt" # name of ensemble dataset

# runs in each ensemble
runs = [0, 1, 2]
runnums = {
  "coalbure" : runs,
  "coalonly" : runs,
  "coalbreakup" : runs,
}

# variables in datasets to create ensemble dataset for
vars4ensemb = ["nsupers", "massmom0", "massmom1", "massmom2"]

# path and filenames for plotting functions
constsfile    = path2CLEO+"/libs/cleoconstants.hpp"
gridfile      = path2build+"/share/buii_dimlessGBxboundaries.dat"

### ------------------------------------------------------------ ###
### ------------ CREATE ENSEMBLE OF DATASETS RESULTS ----------- ###
### ------------------------------------------------------------ ###
for lab in labels:

  # directories for making ensemble dataset
  datapath = binpath+"/"+lab+"/runs/"   # directory of datasets
  ensembdatapath = binpath+"/"+lab+"/ensemb/" # directory to write ensem dataset to

  if path2CLEO == ensembdatapath:
    raise ValueError("ensemble mean directory cannot be CLEO")
  else:
    Path(ensembdatapath).mkdir(exist_ok=True) 

  datasets = [] 
  setupfile = datapath+"/setup_run"+str(runnums[lab][0])+".txt"
  for n in runnums[lab]:
    # setup and zarr for run[n] of ensemble
    runstr = "run"+str(n)
    dataset = datapath+"/sol_"+runstr+".zarr"
    datasets.append(dataset)
  
  ensembdataset = ensembdatapath+ensembzarr
  ensembsetupfile = ensembdatapath+ensembsetuptxt
  write_ensemble(ensembdataset, ensembsetupfile,
                vars4ensemb, setupfile, datasets)