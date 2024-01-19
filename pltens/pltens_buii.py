'''
----- CLEO -----
File: pltens_buii.py
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
Script plots some data from ensemble
dataset of buii 1-D rainshaft output
'''

import os
import sys
import numpy as np
import random 
from pathlib import Path
from matplotlib.colors import LogNorm, Normalize

path2CLEO = sys.argv[1]
path2build = sys.argv[2]
savefigpath = sys.argv[3]

sys.path.append(path2CLEO)  # for imports from pySD package
sys.path.append(path2CLEO+"/examples/exampleplotting/") # for imports from example plotting package

from src.pltens_src import *
from pySD.sdmout_src import *

### ---------------------------------------------------------------- ###
### ----------------------- INPUT PARAMETERS ----------------------- ###
### ---------------------------------------------------------------- ###
### --- essential paths and filenames --- ###
labels = ["coalbure", "coalonly", "coalbreakup", "coalnobure"]

constsfile    = path2CLEO+"/libs/cleoconstants.hpp"
gridfile      = path2build+"/share/buii_dimlessGBxboundaries.dat"

### ------------------------------------------------------------ ###
### ----------------------- PLOT RESULTS ----------------------- ###
### ------------------------------------------------------------ ###
if path2CLEO == savefigpath:
  raise ValueError("plots directory cannot be CLEO")
else:
  Path(savefigpath).mkdir(parents=True, exist_ok=True) 
  

for label in labels:
  
  ### ----- load data to plot ----- ###
  # path and file names for plotting results
  datapath = path2build+"/bin/"+label+"/ensemb/"
  setupfile = datapath+"/setup_ensemb.txt"
  dataset = datapath+"/sol_ensemb.zarr"

  # read in constants and data
  config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=True)
  consts = pysetuptxt.get_consts(setupfile, isprint=True)
  gbxs = pygbxsdat.get_gridboxes(gridfile, consts["COORD0"], isprint=True)

  time = pyzarr.get_time(dataset)
  massmoms = pyzarr.get_massmoms(dataset, config["ntime"], gbxs["ndims"])

  ### ----- plot figures ----- ###
  savename = savefigpath + "domainmassmoms_"+label+".png"
  pltmoms.plot_domainmassmoments(time, massmoms, savename=savename)