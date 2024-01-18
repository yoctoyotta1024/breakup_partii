'''
----- CLEO -----
File: initSDs.py 
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
Script generated initial superdrop
conditions for 1-D rainshaft
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
initSDspath   = sharepath+"buii_dimlessSDsinits/"
runs          = [0, 1, 2, 3]

### --- plotting initialisation figures --- ###
isfigures   = [True, True] # booleans for [making, saving] initialisation figures
savefigpath = path2build+"/bin/initSDs/" # directory for saving init figures
SDgbxs2plt  = 0

### --- settings for initial superdroplets --- ###
# initial superdroplet coordinates
nsupers = 16384                                 # number of superdroplets per gridbox 

# initial superdroplet radii (and implicitly solute masses)
rspan   = [3e-7, 2e-3]                          # max and min range of radii to sample [m]
dryr_sf = 1e3                                   # dryradii are 1/sf of radii [m]

# settings for initial superdroplet multiplicies
numconc               = 1e9                     # droplet number concentration [m^3]      
scalefacs             = [10000, 1]              # Ratio [Cloud : Raindrop] disribution
 
reff                 = 7e-6                     # effective radius [m] (Hansen Cloud droplet Gamma Distribution )
nueff                = 0.08                     # effective variance (Hansen Cloud droplet Gamma Distribution )

nrain                = 3000                     # raindrop concentration [m^-3] (Geoffroy Raindrop Gamma Distribution)
qrain                = 0.9                      # rainwater content [g/m^3] (Geoffroy Raindrop Gamma Distribution)
dvol                 = 8e-4                     # mean volume diameter [m] (Geoffroy Raindrop Gamma Distribution)


### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###

### ---------------------------------------------------------------- ###
### ------------------- BINARY FILES GENERATION--------------------- ###
### ---------------------------------------------------------------- ###
### --- ensure build, share and bin directories exist --- ###
if path2CLEO == path2build:
  raise ValueError("build directory cannot be CLEO")
else:
  Path(path2build).mkdir(exist_ok=True) 
  Path(sharepath).mkdir(exist_ok=True) 
  Path(initSDspath).mkdir(exist_ok=True) 

### ----- create initial superdroplets generator ----- ###
coord3gen = None                        # do not generate superdroplet coord3s
coord1gen = None                        # do not generate superdroplet coord1s
coord2gen = None                        # do not generate superdroplet coord2s
rdist1 = probdists.ClouddropsHansenGamma(reff, nueff)
rdist2 = probdists.RaindropsGeoffroyGamma(nrain, qrain, dvol)
distribs = [rdist1, rdist2]
xiprobdist = probdists.CombinedRadiiProbDistribs(distribs, scalefacs)
radiigen = rgens.SampleLog10RadiiGen(rspan)
dryradiigen =  dryrgens.ScaledRadiiGen(dryr_sf)

initattrsgen = attrsgen.AttrsGenerator(radiigen, dryradiigen, xiprobdist,
                                        coord3gen, coord1gen, coord2gen)

for n in runs:
  runstr = "run"+str(n)
  initSDsfile = initSDspath+"/"+runstr+".dat"
  savefigstem = savefigpath+runstr+"_"

  ### ----- write initial superdroplets binary ----- ###
  csupers.write_initsuperdrops_binary(initSDsfile, initattrsgen, 
                                      tmp_configfile, constsfile,
                                      gridfile, nsupers, numconc)

  ### ----- show (and save) plots of binary file data ----- ###
  if isfigures[0]:
    if isfigures[1]:
      Path(savefigpath).mkdir(exist_ok=True) 
    rsupers.plot_initGBxs_distribs(tmp_configfile, constsfile, initSDsfile,
                                   gridfile, savefigstem, isfigures[1],
                                   SDgbxs2plt) 
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###