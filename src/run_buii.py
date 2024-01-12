'''
----- CLEO -----
File: run_buii.py 
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
Script runs CLEO buii_[exectuable] for 1-D rainshaft
'''

import os
import sys
import numpy as np
import random 
from pathlib import Path
from matplotlib.colors import LogNorm, Normalize

path2CLEO = sys.argv[1]
path2build = sys.argv[2]
tmp_configfile = sys.argv[3]

sys.path.append(path2CLEO)  # for imports from pySD package
from pySD import editconfigfile

### ---------------------------------------------------------------- ###
### ----------------------- INPUT PARAMETERS ----------------------- ###
### ---------------------------------------------------------------- ###
# labels for model compilations
labels = ["coalbure", "coalonly"]

# name of executables and data output
executables = {
  "coalbure" : "buii_coalbure",
  "coalonly" : "buii_coalonly"
}

# initSDs run numbers to use
runnums = [0, 1, 2, 3]
executables = {
  "coalbure" : runnums,
  "coalonly" : runnums
}

binpath       = path2build+"/bin/"
sharepath     = path2build+"/share/"
tmppath       = path2build+"/tmp/"
initSDspath   = sharepath+"/buii_dimlessSDsinits/"
### ---------------------------------------------------------------- ###
### ---------------------------- RUN CLEO -------------------------- ###
### ---------------------------------------------------------------- ###
if path2CLEO == path2build:
  raise ValueError("build directory cannot be CLEO")

for n in runnums:
  runstr = "run"+str(n)

  for lab in labels:
  
    ### ----- copy tmp_config to config then edit ----- ###
    configfile = tmppath+"/config_"+lab+"_"+runstr+".txt"
    os.system('cp '+tmp_configfile+" "+configfile)
    params = {
      "initsupers_filename" : initSDspath+"/"+runstr+".dat",
      "setuptxt" : binpath+"/"+lab+"/setup_"+runstr+".txt",
      "zarrbasedir" : binpath+"/"+lab+"/sol_"+runstr+".zarr",
      "stats_filename" : binpath+"/"+lab+"/stats_"+runstr+".txt",
      }
    editconfigfile.edit_config_params(configfile, params)

    ### ----- run executabel with new configfile ----- ### 
    os.chdir(path2build)
    os.system('pwd')
    executable = path2build+'/src/'+executables[lab]
    os.system(executable + ' ' + configfile)
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###