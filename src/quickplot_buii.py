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
Script plots some raw data from 
buii 1-D rainshaft output
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
### --- essential paths and filenames --- ###
# path and filenames for creating initial SD conditions
constsfile    = path2CLEO+"/libs/cleoconstants.hpp"
gridfile      = path2build+"/share/buii_dimlessGBxboundaries.dat"

# path and file names for plotting results
datapath       = path2build+"/bin/"
setupfile     = datapath+"buiisetup.txt"
dataset       = datapath+"buiisol.zarr"

# directory for saving figures and animations
savefigpath = path2build+"/bin/"

### ------------------------------------------------------------ ###
### ----------------------- PLOT RESULTS ----------------------- ###
### ------------------------------------------------------------ ###
# read in constants and intial setup from setup .txt file
config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=True)
consts = pysetuptxt.get_consts(setupfile, isprint=True)
gbxs = pygbxsdat.get_gridboxes(gridfile, consts["COORD0"], isprint=True)

### ----- load data to plot ----- ###
time = pyzarr.get_time(dataset)
sddata = pyzarr.get_supers(dataset, consts)
totnsupers = pyzarr.get_totnsupers(dataset)
massmoms = pyzarr.get_massmoms(dataset, config["ntime"], gbxs["ndims"])

### ----- plot figures ----- ###
savename = savefigpath + "buii_totnsupers.png"
pltmoms.plot_totnsupers(time, totnsupers, savename=savename)

savename = savefigpath + "buii_domainmassmoms.png"
pltmoms.plot_domainmassmoments(time, massmoms, savename=savename)

nsample = 500
savename = savefigpath + "buii_randomsample.png"
pltsds.plot_randomsample_superdrops(time, sddata,
                                        config["totnsupers"],
                                        nsample,
                                        savename=savename)

### ----- plot 1-D .gif animations ----- ###
nframes = len(time.mins)
mom2ani = np.sum(massmoms.nsupers, axis=(1,2))
xlims = [0, np.amax(mom2ani)]
xlabel = "number of super-droplets"
savename=savefigpath+"buii_nsupers1d"
animations.animate1dprofile(gbxs, mom2ani, time.mins, nframes,
                            xlabel=xlabel, xlims=xlims,
                            color="green", saveani=True,
                            savename=savename, fps=5)   

nframes = len(time.mins)
norm = gbxs["gbxvols"] * 1e6 # volume [cm^3]
mom2ani = np.sum(massmoms.mom0 / norm[None,:], axis=(1,2))
xlims = [0, np.amax(mom2ani)]
xlabel = "number concentration /cm$^{-3}$"
savename=savefigpath+"buii_numconc1d"
animations.animate1dprofile(gbxs, mom2ani, time.mins, nframes,
                            xlabel=xlabel, xlims=xlims,
                            color="green", saveani=True,
                            savename=savename, fps=5)

nframes = len(time.mins)
norm = gbxs["gbxvols"] # volume [m^3]
mom2ani = np.sum(massmoms.mom1/ norm[None,:], axis=(1,2))
xlims = [0, np.amax(mom2ani)]
xlabel = "mass concentration /g m$^{-3}$"
savename=savefigpath+"buii_massconc1d"
animations.animate1dprofile(gbxs, mom2ani, time.mins, nframes,
                            xlabel=xlabel, xlims=xlims,
                            color="green", saveani=True,
                            savename=savename, fps=5)                        
### ------------------------------------------------------------ ###
### ------------------------------------------------------------ ###                                