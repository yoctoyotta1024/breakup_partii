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
configfile_draft = sys.argv[3]

sys.path.append(path2CLEO)  # for imports from pySD package
from pySD.initsuperdropsbinary_src import *
from pySD.initsuperdropsbinary_src import create_initsuperdrops as csupers 
from pySD.initsuperdropsbinary_src import read_initsuperdrops as rsupers 

### ---------------------------------------------------------------- ###
### ----------------------- INPUT PARAMETERS ----------------------- ###
### ---------------------------------------------------------------- ###
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

### ---------------------------------------------------------------- ###
### ---------------------------- RUN CLEO -------------------------- ###
### ---------------------------------------------------------------- ###
os.chdir(path2build)
os.system('pwd')
for exec in executables:
  executable = path2build+'/src/'+exec

  for n in runnums:
    runstr = "run"+str(n)
    
    configfile = 
    os.system(executable + ' ' + configfile)
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###