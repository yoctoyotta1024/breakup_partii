'''
----- CLEO -----
File: pltens_src.py
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
import matplotlib.pyplot as plt
from pathlib import Path

sys.path.append(str(Path.home())+"/CLEO/")  # path2CLEO for imports from pySD package
from pySD.sdmout_src import *

def savefig(fig, savename, show=True):
  
  fig.tight_layout()
  
  fig.savefig(savename, dpi=400, bbox_inches="tight",
              facecolor='w', format="png")
  print("Figure .png saved as: "+savename)

  if show:
    plt.show()

def get_gbxs(datapath, gridfile):

  setupfile = datapath+"/setup_ensemb.txt"
  consts = pysetuptxt.get_consts(setupfile, isprint=True) 
  gbxs = pygbxsdat.get_gridboxes(gridfile, consts["COORD0"], isprint=True)
  
  return gbxs

def get_massmoms(datapath, gridfile):
  
  setupfile = datapath+"/setup_ensemb.txt"
  dataset = datapath+"/sol_ensemb.zarr"

  # read in constants and data
  config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=True)
  consts = pysetuptxt.get_consts(setupfile, isprint=True)
  gbxs = pygbxsdat.get_gridboxes(gridfile, consts["COORD0"], isprint=True)

  time = pyzarr.get_time(dataset)
  massmoms = pyzarr.get_massmoms(dataset, config["ntime"], gbxs["ndims"])

  return time, massmoms

def plot_all_on_axs(path2build, gridfile,
                    fig, axs, plotfunc, datalabs,
                    labels, colors, savename=""):
  
  handles, handlelabs = [], []
  for datalab in datalabs:

    ### ----- load data to plot ----- ###
    # path and file names for plotting results
    datapath = path2build+"/bin/"+datalab+"/ensemb/"
    label = labels[datalab]
    color = colors[datalab]

    ### ----- plot data ----- ###
    line0 = plotfunc(axs, datapath, gridfile, color=color)
    handles.append(line0)
    handlelabs.append(label)
  
  try:
    axs[0].legend(handles, handlelabs)
  except:
    axs.legend(handles, handlelabs)

  if savename != "":
    savefig(fig, savename, show=False)
    
def plot_gbxmassmoments(axs, datapath, gridfile, color="k"):
  ''' plot mass moments 0th gridbox in domain '''

  zgbx=0
  time, massmoms = get_massmoms(datapath, gridfile)

  line0 = axs[0].plot(time.mins, massmoms.nsupers[:,0,0,zgbx], color=color)
  axs[1].plot(time.mins, massmoms.mom0[:,0,0,zgbx], color=color)
  axs[2].plot(time.mins, massmoms.mom1[:,0,0,zgbx], color=color)
  axs[3].plot(time.mins, massmoms.mom2[:,0,0,zgbx], color=color)
  axs[4].plot(time.mins, massmoms.effmass[:,0,0,zgbx], color=color)

  axs[0].set_ylabel("number of\nsuperdroplets")
  axs[1].set_ylabel("$\u03BB^{m}_{0}$, number\nof droplets")
  axs[2].set_ylabel("$\u03BB^{m}_{1}$, droplet\nmass /g")
  axs[3].set_ylabel("$\u03BB^{m}_{2}$\n~reflectivity /g$^2$")
  ylab4 = "mean effective\ndroplet mass,\n<$\u03BB^{m}_{2}$/$\u03BB^{m}_{1}>$ /g"
  axs[4].set_ylabel(ylab4)

  for ax in [axs[1], axs[2]]:
    ax.set_yscale("log")
  axs[-1].set_xlabel("time /min")

  return line0[0]

def plot_gbxnumconc(ax, datapath, gridfile, color="k"):
  ''' plot number concentration of 0th gridbox in domain '''

  zgbx=0
  time, massmoms = get_massmoms(datapath, gridfile)
  gbxs = get_gbxs(datapath, gridfile)
  volcm3 = gbxs["gbxvols"][0,0,zgbx] * 1e6 # [cm^3]
  numconc = massmoms.mom0[:,0,0,zgbx] / volcm3 # [cm^-3]
  
  line = ax.plot(time.mins, numconc, color=color)
  ax.set_ylabel("number concentration /cm$^{-3}$")
  ax.set_yscale("log")
  ax.set_xlabel("time /min")

  return line[0]

def plot_gbxreflectivity(ax, datapath, gridfile, color="k"):
  ''' plot reflectivity proxy of 0th gridbox in domain '''

  zgbx=0
  time, massmoms = get_massmoms(datapath, gridfile)
  refproxy = massmoms.mom2[:,0,0,zgbx] * 1e25 # 10^25 g^2

  line = ax.plot(time.mins, refproxy, color=color)
  ax.set_ylabel("reflectivity proxy /10$^{25}$ g$^2$")
  ax.set_xlabel("time /min")

  return line[0]