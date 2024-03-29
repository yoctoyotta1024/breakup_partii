'''
----- CLEO -----
File: quickplotrun_buii.py
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
import matplotlib.pyplot as plt
from pathlib import Path

path2CLEO = sys.argv[1]
path2build = sys.argv[2]
datapath = sys.argv[3]
runstr = "run"+sys.argv[4]

sys.path.append(path2CLEO)  # for imports from pySD package
sys.path.append(path2CLEO+"/examples/exampleplotting/") # for imports from example plotting package

from plotssrc import pltdist, pltmoms, animations
from pySD.sdmout_src import *             # pyzarr, pysetuptxt & pygbxsdat
from pySD.sdmout_src import sdtracing

### ---------------------------------------------------------------- ###
### ----------------------- INPUT PARAMETERS ----------------------- ###
### ---------------------------------------------------------------- ###
### --- essential paths and filenames --- ###
# path and filenames for plotting functions
constsfile    = path2CLEO+"/libs/cleoconstants.hpp"
gridfile      = path2build+"/share/buii_dimlessGBxboundaries.dat"

# path and file names for plotting results
setupfile     = datapath+"/setup_"+runstr+".txt"
dataset       = datapath+"/sol_"+runstr+".zarr"

# directory for saving figures and animations
pltgifs = False # plot gifs or not
savefigpath = datapath+"/plots/"+runstr+"/"

### ------------------------------------------------------------ ###
### ------------------- EXTRA PLOT FUNCTIONS ------------------- ###
### ------------------------------------------------------------ ###

def savefig(fig, savename):
  fig.savefig(savename, dpi=400, bbox_inches="tight",
              facecolor='w', format="png")
  print("Figure .png saved as: "+savename)

def sample_data(npopln, nsample, sddata,
                attrs=["xi", "radius", "coord3"]):

  data = sdtracing.attrs_for_superdroplets_sample(sddata,
                                                  attrs,
                                                  ndrops2sample=nsample,
                                                  minid=0,
                                                  maxid=int(npopln))
  return data

def plot_sample(npopln, nsample, time, sddata, savename):
  ''' takes random sample of 'nsample' superdroplets from total 
  'npopln' population and plots their radius and coord3 evolution'''

  attrs = ["xi", "radius", "msol", "coord3"]
  sample = sample_data(npopln, nsample, sddata, attrs=attrs)

  fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12,8))    
  axs = axs.flatten()

  fig.suptitle("Random Sample of Superdroplets")
  
  diam = sample["radius"] * 2 / 1e4 #[cm]
  axs[0].plot(time.mins, diam, linewidth=0.8)
  axs[0].set_ylabel("diameter /cm")
  axs[0].set_yscale("log")

  crd3 = sample["coord3"] / 1000 # [km]
  axs[1].plot(time.mins, crd3, linewidth=0.8)
  axs[1].set_ylabel("z /km")

  axs[2].plot(time.mins, sample["msol"], linewidth=0.8)
  axs[2].set_ylabel("solute mass / g")
  axs[2].set_yscale("log")
  axs[2].set_xlabel("time /mins")

  axs[3].plot(time.mins, sample["xi"], linewidth=0.8)
  axs[3].set_ylabel("multiplicity, xi")
  axs[3].set_yscale("log")
  axs[3].set_xlabel("time /mins")

  fig.tight_layout()

  savefig(fig, savename) 
  plt.show()

  return fig, axs

def plot_radius_coord3(time, sample, savename):
  ''' takes random sample of 'nsample' superdroplets from total 
  'npopln' population and plots their radius and coord3 evolution'''

  fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(12,6))    

  fig.suptitle("Random Sample of Superdroplets")
  
  diam = sample["radius"] * 2 / 1e4 #[cm]
  axs[0].plot(time.mins, diam, linewidth=0.8)
  axs[0].set_xlabel("time /mins")
  axs[0].set_ylabel("diameter /cm")
  axs[0].set_yscale("log")

  crd3 = sample["coord3"] / 1000 # [km]
  axs[1].plot(time.mins, crd3, linewidth=0.8)
  axs[1].set_xlabel("time /mins")
  axs[1].set_ylabel("z /km")

  axs[2].plot(diam, crd3, linewidth=0.8)
  axs[2].set_xlabel("diameter /cm")
  axs[2].set_xscale("log")
  axs[2].set_ylabel("z /km")

  fig.tight_layout()

  savefig(fig, savename)
  plt.show()

  return fig, axs


def plot_xi_radius_coord3(time, sample, savename):
  ''' takes random sample of 'nsample' superdroplets from total 
  'npopln' population and plots their radius and coord3 evolution'''

  fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(12,6))    

  fig.suptitle("Random Sample of Superdroplets")
  
  xi = sample["xi"]
  axs[0].plot(time.mins, xi, linewidth=0.8)
  axs[0].set_xlabel("time /mins")
  axs[0].set_ylabel("multiplicity")
  axs[0].set_yscale("log")

  crd3 = sample["coord3"] / 1000 # [km]
  axs[1].plot(xi, crd3, linewidth=0.8)
  axs[1].set_ylabel("z /km")
  axs[1].set_xlabel("multiplicity")
  axs[1].set_xscale("log")

  diam = sample["radius"] * 2 / 1e4 #[cm]
  axs[2].plot(diam, xi, linewidth=0.8)
  axs[2].set_xlabel("diameter /cm")
  axs[2].set_xscale("log")
  axs[2].set_yscale("log")

  fig.tight_layout()

  savefig(fig, savename)
  plt.show()

  return fig, axs

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
sddata = pyzarr.get_supers(dataset, consts)
totnsupers = pyzarr.get_totnsupers(dataset)
massmoms = pyzarr.get_massmoms(dataset, config["ntime"], gbxs["ndims"])

### ----- plot domain figures ----- ###
savename = savefigpath + "domainmassmoms.png"
pltmoms.plot_domainmassmoments(time, massmoms, savename=savename)

t2plts = np.linspace(0, time.secs[-1], 15)
rspan = [np.nanmin(sddata["radius"]), np.nanmax(sddata["radius"])]
nbins = 200
smoothsig = False

savename = savefigpath + "domainnconcdist.png"
fig, ax = pltdist.plot_domainnumconc_distribs(time.secs, sddata, t2plts, 
                                     gbxs["domainvol"], rspan, nbins,
                                     smoothsig=smoothsig,
                                     perlogR=False,
                                     ylog=True) 
savefig(fig, savename)

savename = savefigpath + "domainmassdist.png"
fig, ax = pltdist.plot_domainmass_distribs(time.secs, sddata, t2plts, 
                                     gbxs["domainvol"], rspan, nbins,
                                     smoothsig=smoothsig,
                                     perlogR=False,
                                     ylog=True)
savefig(fig, savename)

savename = savefigpath + "domainnsupersdist.png"
fig, ax = pltdist.plot_domainnsupers_distribs(time.secs, sddata, t2plts, 
                                     gbxs["domainvol"], rspan, nbins,
                                     smoothsig=smoothsig,
                                     perlogR=True,
                                     ylog=False)
savefig(fig, savename)

### ----- plot random sample figures ----- ###
nsample = 250
savename = savefigpath + "randomsample.png"
plot_sample(totnsupers[0], nsample, time, sddata, savename)

nsample = 250
sample = sample_data(totnsupers[0], nsample, sddata)
savename = savefigpath + "randomsample_radiuscoord3.png"
plot_radius_coord3(time, sample, savename)

savename = savefigpath + "randomsample_radiusxi.png"
plot_xi_radius_coord3(time, sample, savename)

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