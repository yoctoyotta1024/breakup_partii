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
datalabs = ["coalre", "coalbu", "coalbure"]
labels = {
  "coalbure": "Coal+Bu+Re", 
  "coalbu": "Coal+Bu",
  "coalre": "Coal+Re",
}
colors = {
  "coalbure": "C0", 
  "coalbu": "C2",
  "coalre": "C3",
}

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
    fig, axs = plt.subplots(nrows=5, ncols=1, figsize=(6,8), sharex=True)
    fig.suptitle("Total Mass Moments Over Domain")
    plotfunc = src.plot_gbxmassmoments
    args = [gridfile]
    src.plot_all_on_axs(path2build, plotfunc, args, fig, axs,
                        datalabs, labels, colors,
                        savename=savename) 
    
  savename = savefigpath + "massmoments.png"
  plot_all_massmoments(datalabs, savename=savename)

  ### ----- plot domain number concentration ----- ###
  def plot_all_numconc(datalabs, savename=""):
    fig, axs = plt.subplots(figsize=(6,8))
    plotfunc = src.plot_gbxnumconc
    args = [gridfile]
    src.plot_all_on_axs(path2build, plotfunc, args, fig, axs,
                        datalabs, labels, colors,
                        savename=savename) 

  savename = savefigpath + "numconc.png"
  plot_all_numconc(datalabs, savename=savename)

  ### ----- plot domain reflectivity ----- ###
  def plot_all_reflectivity(datalabs, savename=""):
    fig, axs = plt.subplots(figsize=(6,8))
    plotfunc = src.plot_gbxreflectivity
    args = [gridfile]
    src.plot_all_on_axs(path2build, plotfunc, args, fig, axs,
                        datalabs, labels, colors,
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
                        datalabs, labels, colors,
                        savename=savename) 

  trange = [0, 70, 10] #[s]

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
  
  trange = [0, 60, 10] #[s]

  savename = savefigpath + "prob_colls.png"
  probcalc = src.log10_collprob
  levels = np.linspace(-10, 5, 50)
  plot_all_prob(datalabs, probcalc, trange, levels, savename=savename)

  savename = savefigpath + "prob_collcoal.png"
  probcalc = src.relative_collcoal_probability
  levels = np.linspace(-20, -5, 50)
  plot_all_prob(datalabs, probcalc, trange, levels, savename=savename)

  savename = savefigpath + "prob_collbreakup.png"
  probcalc = src.relative_collbreakup_probability
  levels = np.linspace(-16, -5, 50)
  plot_all_prob(datalabs, probcalc, trange, levels, savename=savename)

  savename = savefigpath + "prob_collrebound.png"
  probcalc = src.relative_collrebound_probability
  levels = np.linspace(-16, -5, 50)
  plot_all_prob(datalabs, probcalc, trange, levels, savename=savename)
