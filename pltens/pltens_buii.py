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

import sys
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

path2CLEO = sys.argv[1]
path2build = sys.argv[2]
savefigpath = sys.argv[3]

import src.pltens_src as src

### ---------------------------------------------------------------- ###
### ----------------------- INPUT PARAMETERS ----------------------- ###
### ---------------------------------------------------------------- ###
### --- essential paths and filenames --- ###
datalabs = ["coalbure", "coalonly", "coalbreakup", "coalnobure"]
labels = {
  "coalbure": "CoalBuRe", 
  "coalonly": "Orig", 
  "coalbreakup": "CoalBu",
  "coalnobure": "Coal",
}
colors = {
  "coalbure": "C0", 
  "coalonly": "C1", 
  "coalbreakup": "C2",
  "coalnobure": "C3",
}

zgbx          = 0 # z index of gridbox to plot for mass moments
constsfile    = path2CLEO+"/libs/cleoconstants.hpp"
gridfile      = path2build+"/share/buii_dimlessGBxboundaries.dat"
### ------------------------------------------------------------ ###
### ----------------------- PLOT RESULTS ----------------------- ###
### ------------------------------------------------------------ ###
if path2CLEO == savefigpath:
  raise ValueError("plots directory cannot be CLEO")
else:
  Path(savefigpath).mkdir(parents=True, exist_ok=True) 
  
def plot_all_on_axs(fig, axs, plotfunc, datalabs, savename=""):
  
  handles, handlelabs = [], []
  for datalab in datalabs:

    ### ----- load data to plot ----- ###
    # path and file names for plotting results
    datapath = path2build+"/bin/"+datalab+"/ensemb/"
    label = labels[datalab]
    color = colors[datalab]

    ### ----- plot data ----- ###
    line0 = plotfunc(axs, zgbx, datapath, gridfile, color=color)
    handles.append(line0)
    handlelabs.append(label)
  
  try:
    axs[0].legend(handles, handlelabs)
  except:
    axs.legend(handles, handlelabs)

  if savename != "":
    src.savefig(fig, savename, show=False)
  
### ----- plot domain mass moments ----- ###
def plot_all_massmoments(datalabs, savename=""):
  fig, axs = plt.subplots(nrows=5, ncols=1, figsize=(6,8), sharex=True)
  fig.suptitle("Total Mass Moments Over Domain")
  plotfunc = src.plot_gbxmassmoments
  plot_all_on_axs(fig, axs, plotfunc, datalabs, savename=savename) 
  
savename = savefigpath + "massmoments.png"
plot_all_massmoments(datalabs, savename=savename)

### ----- plot domain number concentration ----- ###
def plot_all_numconc(datalabs, savename=""):
  fig, axs = plt.subplots(figsize=(6,8))
  plotfunc = src.plot_gbxnumconc
  plot_all_on_axs(fig, axs, plotfunc, datalabs, savename=savename) 

savename = savefigpath + "numconc.png"
plot_all_numconc(datalabs, savename=savename)

### ----- plot domain reflectivity ----- ###
def plot_all_reflectivity(datalabs, savename=""):
  fig, axs = plt.subplots(figsize=(6,8))
  plotfunc = src.plot_gbxreflectivity
  plot_all_on_axs(fig, axs, plotfunc, datalabs, savename=savename) 
  
savename = savefigpath + "reflectivity.png"
plot_all_reflectivity(datalabs, savename=savename)