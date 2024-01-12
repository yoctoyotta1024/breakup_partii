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
configfile = sys.argv[3]
executable = sys.argv[4]

sys.path.append(path2CLEO)  # for imports from pySD package
from pySD.initsuperdropsbinary_src import *
from pySD.initsuperdropsbinary_src import create_initsuperdrops as csupers 
from pySD.initsuperdropsbinary_src import read_initsuperdrops as rsupers 

### ---------------------------------------------------------------- ###
### ----------------------- INPUT PARAMETERS ----------------------- ###
### ---------------------------------------------------------------- ###
### --- essential paths and filenames --- ###
# path and filenames for creating initial SD conditions
constsfile    = path2CLEO+"/libs/cleoconstants.hpp"
sharepath     = path2build+"/share/"
gridfile      = sharepath+"buii_dimlessGBxboundaries.dat"
initSDsfile   = sharepath+"buii_dimlessSDsinit.dat"

### --- plotting initialisation figures --- ###
isfigures   = [True, True] # booleans for [making, saving] initialisation figures
savefigpath = path2build+"/bin/" # directory for saving init figures
SDgbxs2plt  = list(range(39, 99))
SDgbxs2plt  = [random.choice(SDgbxs2plt)] # choose random gbx from list to plot 

### --- settings for initial superdroplets --- ###
# initial superdroplet coordinates
zlim        = 800       # min z coord of superdroplets [m]
npergbx     = 256       # number of superdroplets per gridbox 

# initial superdroplet radii (and implicitly solute masses)
rspan       = [3e-9, 5e-5]                      # min and max range of radii to sample [m]
dryr_sf     = 1.0                               # dryradii are 1/sf of radii [m]

# settings for initial superdroplet multiplicies
geomeans             = [0.02e-6, 0.2e-6, 3.5e-6]               
geosigs              = [1.55, 2.3, 2]                    
scalefacs            = [1e6, 0.3e6, 0.025e6]   
numconc = np.sum(scalefacs) * 100
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###

### ---------------------------------------------------------------- ###
### ---------------------------- RUN CLEO -------------------------- ###
### ---------------------------------------------------------------- ###
os.chdir(path2build)
os.system('pwd')
executable = path2build+'/src/'+executable
os.system(executable + ' ' + configfile)
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###