'''
----- CLEO -----
File: quickplotens_buii.py
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

import sys
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

path2CLEO = sys.argv[1]
path2build = sys.argv[2]
datapath = sys.argv[3]

sys.path.append(path2CLEO)  # for imports from pySD package
sys.path.append(path2CLEO+"/examples/exampleplotting/") # for imports from example plotting package

from plotssrc import pltmoms, animations
from pySD.sdmout_src import *             # pyzarr, pysetuptxt & pygbxsdat

### ---------------------------------------------------------------- ###
### ----------------------- INPUT PARAMETERS ----------------------- ###
### ---------------------------------------------------------------- ###
### --- essential paths and filenames --- ###
# path and filenames for plotting functions
constsfile    = path2CLEO+"/libs/cleoconstants.hpp"
gridfile      = path2build+"/share/buii_dimlessGBxboundaries.dat"

# path and file names for plotting results
setupfile     = datapath+"/setup_ensemb.txt"
dataset       = datapath+"/sol_ensemb.zarr"

# directory for saving figures and animations
pltgifs = False # plot gifs or not
savefigpath = datapath+"/plots/"

### ------------------------------------------------------------ ###
### --------------------- EXTRA PLOT FUNCS --------------------- ###
### ------------------------------------------------------------ ###
def savefig(fig, savename):
  fig.savefig(savename, dpi=400, bbox_inches="tight",
              facecolor='w', format="png")
  print("Figure .png saved as: "+savename)

def plot_distrib(ax, redges, mean, std, logy=False):

  ntime = mean.shape[0] # number of timesthat distribution has
  n2plt = 10 # number of distirbutiosn to plot
  ts2plt = list(range(0, ntime, ntime//n2plt)) # index of times to plot
  colors = plt.cm.plasma(np.linspace(0.2, 0.8, n2plt))
  
  for i, c in zip(ts2plt, colors):
    ax.step(redges[:-1], mean[i, :], where='pre', color=c)

    ax.step(redges[:-1], (mean-std)[i, :],
            where='pre', color=c, linestyle="--")
    ax.step(redges[:-1], (mean+std)[i, :],
            where='pre', color=c, linestyle="--") 
  
  if logy:
    ax.set_yscale("log")
  ax.set_xscale("log")

def plot_domaindistribs(dataset, savename=""):
  
  fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(8,12),
                          sharex=True)

  ds = pyzarr.get_rawdataset(dataset)
  redges = ds["h_redges"]

  mean, std = ds["h_numconc"], ds["h_numconcstd"]
  plot_distrib(axs[0], redges, mean, std, logy=True)
  axs[0].set_ylabel("droplet number concentration /cm$^{-3}$")

  mean, std = ds["h_watermass"], ds["h_watermassstd"]
  plot_distrib(axs[1], redges, mean, std, logy=True)
  axs[1].set_ylabel("droplet mass concentration /g m$^{-3}$")

  mean, std = ds["h_refproxy"], ds["h_refproxystd"]
  plot_distrib(axs[2], redges, mean, std, logy=True)
  axs[2].set_ylabel("reflectivity proxy /m$^{6}$")

  axs[-1].set_xlabel("radius /\u03BCm")

  fig.tight_layout()
  if savename != "":
    savefig(fig, savename)

### ------------------------------------------------------------ ###
### ------------------------------------------------------------ ###                                

### ------------------------------------------------------------ ###
### ----------------------- PLOT RESULTS ----------------------- ###
### ------------------------------------------------------------ ###
if path2CLEO == savefigpath:
  raise ValueError("plots directory cannot be CLEO")
else:
  Path(savefigpath).mkdir(parents=True, exist_ok=True) 
  
# read in constants and intial setup from setup .txt file
config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=True)
consts = pysetuptxt.get_consts(setupfile, isprint=True)
gbxs = pygbxsdat.get_gridboxes(gridfile, consts["COORD0"], isprint=True)

### ----- load data to plot ----- ###
time = pyzarr.get_time(dataset)
massmoms = pyzarr.get_massmoms(dataset, config["ntime"], gbxs["ndims"])

### ----- plot figures ----- ###
savename = savefigpath + "domainmassmoms.png"
pltmoms.plot_domainmassmoments(time, massmoms, savename=savename)

savename = savefigpath + "domaindistribs.png"
plot_domaindistribs(dataset, savename=savename)

### ----- plot 1-D .gif animations ----- ###
if pltgifs:
  nframes = len(time.mins)
  mom2ani = np.sum(massmoms.nsupers, axis=(1,2))
  xlims = [0, np.amax(mom2ani)]
  xlabel = "number of super-droplets"
  savename=savefigpath+"nsupers1d"
  animations.animate1dprofile(gbxs, mom2ani, time.mins, nframes,
                              xlabel=xlabel, xlims=xlims,
                              color="green", saveani=True,
                              savename=savename, fps=5)   

  nframes = len(time.mins)
  norm = gbxs["gbxvols"] * 1e6 # volume [cm^3]
  mom2ani = np.sum(massmoms.mom0 / norm[None,:], axis=(1,2))
  xlims = [0, np.amax(mom2ani)]
  xlabel = "number concentration /cm$^{-3}$"
  savename=savefigpath+"numconc1d"
  animations.animate1dprofile(gbxs, mom2ani, time.mins, nframes,
                              xlabel=xlabel, xlims=xlims,
                              color="green", saveani=True,
                              savename=savename, fps=5)

  nframes = len(time.mins)
  norm = gbxs["gbxvols"] # volume [m^3]
  mom2ani = np.sum(massmoms.mom1/ norm[None,:], axis=(1,2))
  xlims = [0, np.amax(mom2ani)]
  xlabel = "mass concentration /g m$^{-3}$"
  savename=savefigpath+"massconc1d"
  animations.animate1dprofile(gbxs, mom2ani, time.mins, nframes,
                              xlabel=xlabel, xlims=xlims,
                              color="green", saveani=True,
                              savename=savename, fps=5)                        
### ------------------------------------------------------------ ###
### ------------------------------------------------------------ ###                                