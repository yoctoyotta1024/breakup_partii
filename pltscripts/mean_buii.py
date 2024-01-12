'''
----- CLEO -----
File: quickplot_buii.py
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
Script plots mean of ensemble of data
from buii 1-D rainshaft output
'''

import os
import sys
import numpy as np
import random 
from pathlib import Path
from matplotlib.colors import LogNorm, Normalize

path2CLEO = sys.argv[1]
path2build = sys.argv[2]

sys.path.append(path2CLEO)  # for imports from pySD package
sys.path.append(path2CLEO+"/examples/exampleplotting/") # for imports from example plotting package

from plotssrc import pltsds, pltmoms, animations
from pySD.sdmout_src import *

### ---------------------------------------------------------------- ###
### ----------------------- INPUT PARAMETERS ----------------------- ###
### ---------------------------------------------------------------- ###
# label and path for each ensembles of datasets 
labels = ["coalbure", "coalonly"]
binpath = path2build+"/bin/" # path before directory called "label" containing zarr datasets
meanzarr = "sol_ensemb.zarr" # name of ensemble dataset

# runs in each ensemble
runs = [0, 1, 2, 3]
runnums = {
  "coalbure" : runs,
  "coalonly" : runs
}

# path and filenames for plotting functions
constsfile    = path2CLEO+"/libs/cleoconstants.hpp"
gridfile      = path2build+"/share/buii_dimlessGBxboundaries.dat"

### ------------------------------------------------------------ ###
### ------------ CREATE ENSEMBLE OF DATASETS RESULTS ----------- ###
### ------------------------------------------------------------ ###

### ------------------------------------------------------------ ###
### ------------------------------------------------------------ ###                                
def write_ensemble_dataset(meandataset, lab, datapath, runnums):
  
  for n in runnums[lab]:
    # setup and zarr for run[n] of ensemble
    runstr = "run"+str(n)
    setupfile     = datapath+"/setup_"+runstr+".txt"
    dataset       = datapath+"/sol_"+runstr+".zarr"

    print(dataset, setupfile)
  print(meandataset)
  
### ------------------------------------------------------------ ###
### ------------ CREATE ENSEMBLE OF DATASETS RESULTS ----------- ###
### ------------------------------------------------------------ ###
for lab in labels:

  # directories for making ensemble dataset
  datapath = binpath+"/"+lab+"/runs/"   # directory of datasets
  meandatapath = binpath+"/"+lab+"/ensemb/" # directory to write ensem dataset to
  meandataset = meandatapath+meanzarr

  if path2CLEO == meandatapath:
    raise ValueError("ensemble mean directory cannot be CLEO")
  else:
    Path(meandatapath).mkdir(exist_ok=True) 

  write_ensemble_dataset(meandataset, lab, datapath, runnums) 