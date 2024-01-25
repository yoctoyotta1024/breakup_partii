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
what2plot = sys.argv[4]

import src.pltens_src as src

### ---------------------------------------------------------------- ###
### ----------------------- INPUT PARAMETERS ----------------------- ###
### ---------------------------------------------------------------- ###
### --- essential paths and filenames --- ###
nfrags = ["2p6", "4", "8", "16", "32", "128", "256", "512", "2048"]
datalabs_1 = [lab+"/coalbure" for lab in nfrags]
datalabs = ["2p6/coalre"] + datalabs_1
datalabs = ["bin_nfrags"+lab for lab in datalabs]

print("datalabs: ", datalabs)

labels = {
  "bin_nfrags2p6/coalre": "no breakup",
  "bin_nfrags2p6/coalbure": "\u03A9 = 2.6", 
  "bin_nfrags4/coalbure": "\u03A9 = 4", 
  "bin_nfrags8/coalbure": "\u03A9 = 8", 
  "bin_nfrags16/coalbure": "\u03A9 = 16", 
  "bin_nfrags32/coalbure": "\u03A9 = 32", 
  "bin_nfrags128/coalbure": "\u03A9 = 128",
  "bin_nfrags256/coalbure": "\u03A9 = 256", 
  "bin_nfrags512/coalbure": "\u03A9 = 512", 
  "bin_nfrags2048/coalbure": "\u03A9 = 2048", 
}

colors = {
  "bin_nfrags2p6/coalre": "grey",
  "bin_nfrags2p6/coalbure": "fuchsia",
  "bin_nfrags4/coalbure": "purple",
  "bin_nfrags8/coalbure": "steelblue", 
  "bin_nfrags16/coalbure": "deepskyblue", 
  "bin_nfrags32/coalbure": "green",
  "bin_nfrags128/coalbure": "olive", 
  "bin_nfrags256/coalbure": "orange", 
  "bin_nfrags512/coalbure": "red", 
  "bin_nfrags2048/coalbure": "brown", 
}

linestyles = {}
for datalab in datalabs:
  if "coalbure" in datalab:
    linestyles[datalab] = "solid"
  else:
    linestyles[datalab] = "solid"

constsfile    = path2CLEO+"/libs/cleoconstants.hpp"
gridfile      = path2build+"/share/buii_dimlessGBxboundaries.dat"
### ------------------------------------------------------------ ###
### ----------------------- PLOT RESULTS ----------------------- ###
### ------------------------------------------------------------ ###
if path2CLEO == savefigpath:
  raise ValueError("plots directory cannot be CLEO")
else:
  Path(savefigpath).mkdir(parents=True, exist_ok=True) 
  
if what2plot == "massmoms":

  ### ----- plot domain mass moments ----- ###
  def plot_all_massmoments(datalabs, savename=""):
    fig, axs = plt.subplots(nrows=5, ncols=1, figsize=(8,12), sharex=True)
    fig.suptitle("Total Mass Moments Over Domain")
    plotfunc = src.plot_gbxmassmoments
    args = [gridfile]
    src.plot_all_on_axs(path2build, plotfunc, args, fig, axs,
                        datalabs, labels, colors, linestyles,
                        savename=savename) 
    
  savename = savefigpath + "massmoments.png"
  plot_all_massmoments(datalabs, savename=savename)

  ### ----- plot domain number concentration ----- ###
  def plot_all_numconc(datalabs, savename=""):
    fig, axs = plt.subplots(figsize=(6,8))
    plotfunc = src.plot_gbxnumconc
    args = [gridfile]
    src.plot_all_on_axs(path2build, plotfunc, args, fig, axs,
                        datalabs, labels, colors, linestyles,
                        savename=savename) 

  savename = savefigpath + "numconc.png"
  plot_all_numconc(datalabs, savename=savename)

  ### ----- plot domain reflectivity ----- ###
  def plot_all_reflectivity(datalabs, savename=""):
    fig, axs = plt.subplots(figsize=(6,8))
    plotfunc = src.plot_gbxreflectivity
    args = [gridfile]
    src.plot_all_on_axs(path2build, plotfunc, args, fig, axs,
                        datalabs, labels, colors, linestyles,
                        savename=savename) 
    
  savename = savefigpath + "reflectivity.png"
  plot_all_reflectivity(datalabs, savename=savename)


if what2plot == "dists":
  ### --- plot domain droplet distibutions --- ###
  def plot_all_distrib(datalabs, plotfunc, trange, savename=""):

    fig, axs = plt.subplots(nrows=4, ncols=2, figsize=(16,12))
    axs = axs.flatten()  
    t2plts = np.arange(trange[0], trange[1]+trange[2], trange[2]) # [s]
    print("t2plts: ", t2plts)
    args = [t2plts]
    src.plot_all_on_axs(path2build, plotfunc, args, fig, axs,
                        datalabs, labels, colors, linestyles,
                        savename=savename) 

  trange = [0, 1400, 200] #[s]

  plotfunc = src.plot_domainnumconc_dist
  savename = savefigpath + "dist_numconc.png"
  plot_all_distrib(datalabs, plotfunc, trange, savename=savename)

  plotfunc = src.plot_domainwatermass_dist
  savename = savefigpath + "dist_watermass.png"
  plot_all_distrib(datalabs, plotfunc, trange, savename=savename)

  plotfunc = src.plot_domainreflectivity_dist
  savename = savefigpath + "dist_reflectivity.png"
  plot_all_distrib(datalabs, plotfunc, trange, savename=savename)

if what2plot == "probs":
  ### --- plot domain collision probabilities --- ###
  def plot_all_prob(datalabs, probcalc, trange, levels, savename=""):

    fig, axs = plt.subplots(nrows=7, ncols=3, figsize=(12,18))
    plotfunc = src.plot_collisions_overtime
    t2plts = np.arange(trange[0], trange[1]+trange[2], trange[2]) # [s]
    args = [t2plts, probcalc, levels]
    src.plot_all_on_fig(path2build, plotfunc, args, fig, axs,
                        datalabs, labels, colors, 
                        savename=savename) 
  
  trange = [0, 1200, 200] #[s]
  levels = np.linspace(-16, 5, 50)

  savename = savefigpath + "prob_colls.png"
  probcalc = src.log10_collprob
  plot_all_prob(datalabs, probcalc, trange, levels, savename=savename)

  savename = savefigpath + "prob_collcoal.png"
  probcalc = src.log10_collcoal_prob
  plot_all_prob(datalabs, probcalc, trange, levels, savename=savename)

  savename = savefigpath + "prob_collbreakup.png"
  probcalc = src.log10_collbreakup_prob
  plot_all_prob(datalabs, probcalc, trange, levels, savename=savename)

  savename = savefigpath + "prob_collrebound.png"
  probcalc = src.log10_collrebound_prob
  plot_all_prob(datalabs, probcalc, trange, levels, savename=savename)
