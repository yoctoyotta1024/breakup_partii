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
from pySD.sdmout_src import *             # pyzarr, pysetuptxt & pygbxsdat
from pySD.sdmout_src import massmoms
from .probcalcs import *

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

def get_ntime_ndims(setupfile, gridfile):
  
  config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=True)
  consts = pysetuptxt.get_consts(setupfile, isprint=True)
  gbxs = pygbxsdat.get_gridboxes(gridfile, consts["COORD0"], isprint=True)

  return config["ntime"], gbxs["ndims"]

def get_time_massmoms(datapath, gridfile):
  
  setupfile = datapath+"/setup_ensemb.txt"
  dataset = datapath+"/sol_ensemb.zarr"

  # read in data
  ntime, ndims = get_ntime_ndims(setupfile, gridfile)
  time = pyzarr.get_time(dataset)
  massmoms = pyzarr.get_massmoms(dataset, ntime, ndims)

  return time, massmoms

def get_massmomstds(datapath, gridfile):
  
  setupfile = datapath+"/setup_ensemb.txt"
  dataset = datapath+"/sol_ensemb.zarr"
  ntime, ndims = get_ntime_ndims(setupfile, gridfile)

  return massmoms.MassMoms(dataset, ntime, ndims, lab="_std")

def get_time_dist(datapath, distname, std=False, rcens=False):

  dataset = datapath+"/sol_ensemb.zarr"

  ds = pyzarr.get_rawdataset(dataset)
  
  time = pyzarr.get_time(dataset)
  mean = ds["h_"+distname]

  if rcens:
    rbins = ds["h_rcens"]
  else:
    rbins = ds["h_redges"]
  
  if std:
    std = ds["h_"+distname+"std"]
    return time, rbins, mean, std
  
  else:
    return time, rbins, mean


def plot_all_on_axs(path2build, plotfunc, args,
                    fig, axs, datalabs, labels,
                    colors, savename=""):
  ''' load dataset, label and color for each datalab
  in datalabs then use plotfunc to plot each dataset
  on figure axis/axes '''

  handles, handlelabs = [], []
  for datalab in datalabs:

    ### ----- load data to plot ----- ###
    # path and file names for plotting results
    datapath = path2build+"/bin/"+datalab+"/ensemb/"
    label = labels[datalab]
    color = colors[datalab]

    ### ----- plot data ----- ###
    line0 = plotfunc(axs, datapath, color, *args)
    handles.append(line0)
    handlelabs.append(label)
  
  try:
    axs[0].legend(handles, handlelabs)
  except:
    axs.legend(handles, handlelabs)

  if savename != "":
    savefig(fig, savename, show=False)

def plot_all_on_fig(path2build, plotfunc, args,
                    fig, axs, datalabs, labels,
                    colors, savename=""):
  ''' load dataset, label and color for each datalab
  in datalabs then use plotfunc to plot each dataset
  on axis/axes in seperate columns of figure '''
    
  if axs.shape[1] != len(datalabs):
    raise ValueError("number of columns to figure"+\
                     " must equal number of datasets")

  handles, handlelabs = [], []
  for d, datalab in enumerate(datalabs):
    axs_d = axs[:,d]

    ### ----- load data to plot ----- ###
    # path and file names for plotting results
    datapath = path2build+"/bin/"+datalab+"/ensemb/"

    ### ----- plot data ----- ###
    plotfunc(axs_d, datapath, datalab, *args)
    title = labels[datalab]+"\n"+axs_d[0].get_title()
    axs_d[0].set_title(title, color=colors[datalab])

  if savename != "":
    savefig(fig, savename, show=False)

def plot_gbxmassmoments(axs, datapath, color, gridfile):
  ''' plot mass moments 0th gridbox in domain '''

  zgbx=0
  time, massmoms = get_time_massmoms(datapath, gridfile)

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

def plot_gbxnumconc(ax, datapath, color, gridfile):
  ''' plot number concentration of 0th gridbox in domain '''

  zgbx=0
  gbxs = get_gbxs(datapath, gridfile)
  volcm3 = gbxs["gbxvols"][0,0,zgbx] * 1e6 # [cm^3]
  
  time, massmoms = get_time_massmoms(datapath, gridfile)
  numconc = massmoms.mom0[:,0,0,zgbx] / volcm3 # [cm^-3]
  
  massmomstds = get_massmomstds(datapath, gridfile)
  std = massmomstds.mom0[:,0,0,zgbx] / volcm3 # [cm^-3] 
  
  line = ax.plot(time.mins, numconc, color=color)
  ax.fill_between(time.mins, numconc-std, numconc+std,
                  color=color, alpha=0.2)

  ax.set_ylabel("number concentration /cm$^{-3}$")
  ax.set_yscale("log")
  ax.set_xlabel("time /min")

  return line[0]

def plot_gbxreflectivity(ax, datapath, color, gridfile):
  ''' plot reflectivity proxy of 0th gridbox in domain '''

  zgbx=0
  time, massmoms = get_time_massmoms(datapath, gridfile)
  refproxy = massmoms.mom2[:,0,0,zgbx] * 1e25 # 10^25 g^2

  massmomstds = get_massmomstds(datapath, gridfile)
  std = massmomstds.mom2[:,0,0,zgbx] * 1e25 # 10^25 g^2

  line = ax.plot(time.mins, refproxy, color=color)
  ax.fill_between(time.mins, refproxy-std, refproxy+std,
                  color=color, alpha=0.2)
  ax.set_ylabel("reflectivity proxy /10$^{25}$ g$^2$")
  ax.set_yscale("log")
  ax.set_xlabel("time /min")

  return line[0]

def plot_distribs_overtime(axs, t2plts, time, redges, mean, std,
                           color="k", ylab=None, logy=False):

  if len(axs) != len(t2plts):
    raise ValueError("number of times to plot != number of axes")
  
  for n in range(len(t2plts)):
    ax = axs[n]
    idx = np.argmin(abs(time.secs-t2plts[n])) # index of time to plot
    t2plt = time.secs[idx] # [min]
    tlab = "t = {:.1f}secs".format(t2plt)

    line = ax.step(redges[:-1], mean[idx, :], where='pre',
                   color=color, linewidth=0.8)

    ax.step(redges[:-1], (mean-std)[idx, :], where='pre',
            color=color, linewidth=0.8, linestyle="--")
    ax.step(redges[:-1], (mean+std)[idx, :], where='pre',
            color=color, linewidth=0.8, linestyle="--") 

    ax.set_title(tlab)
    ax.set_ylabel(ylab)
    if logy:
      ax.set_yscale("log")

    ax.set_xlabel("radius /\u03BCm")
    ax.set_xscale("log")

  return line[0]

def plot_domainnumconc_dist(axs, datapath, color, t2plts):
  ''' plots seperate distribution for each time in
  t2plts [s] on each axis in axs '''

  time, redges, mean, std = get_time_dist(datapath, "numconc",
                                          std=True, rcens=False)

  ylab = "number concentration /cm$^{-3}$"
  line = plot_distribs_overtime(axs, t2plts, time, redges, mean, std,
                                color=color, ylab=ylab, logy=True)
  
  for ax in axs:
    ax.set_ylim([1e-10, 175])
  
  return line

def plot_domainwatermass_dist(axs, datapath, color, t2plts):
  ''' plots seperate distribution for each time in
  t2plts [s] on each axis in axs '''

  time, redges, mean, std = get_time_dist(datapath, "watermass",
                                          std=True, rcens=False)

  ylab = "mass concentration /g m$^{-3}$"
  line = plot_distribs_overtime(axs, t2plts, time, redges, mean, std,
                                color=color, ylab=ylab, logy=True)
  
  for ax in axs:
    ax.set_ylim([1e-13, 2.5])
  
  return line

def plot_domainreflectivity_dist(axs, datapath, color, t2plts):
  ''' plots seperate distribution for each time in
  t2plts [s] on each axis in axs '''

  time, redges, mean, std = get_time_dist(datapath, "refproxy",
                                          std=True, rcens=False)

  ylab = "reflectivity proxy /m$^{6}$"
  line = plot_distribs_overtime(axs, t2plts, time, redges, mean, std,
                                color=color, ylab=ylab, logy=True)
  
  for ax in axs:
    ax.set_ylim([1e-14, 2e-8])
  
  return line

def plot_collisions_overtime(axs, datapath, datalab, t2plts, probcalc, levels):
  ''' plot probability for each time in t2plts given probcalc function'''

  time, rcens, numconc = get_time_dist(datapath, "numconc",
                                        std=False, rcens=True)

  if len(axs) != len(t2plts):
    raise ValueError("number of times to plot != number of axes")
  
  for n in range(len(t2plts)):
    ax = axs[n]
    idx = np.argmin(abs(time.secs-t2plts[n])) # index of time to plot
    t2plt = time.secs[idx] # [min]
    tlab = "t = {:.1f}secs".format(t2plt)

    rr1, rr2, prob = probcalc(datalab, rcens, numconc[idx, :])
    ax.contourf(rr1, rr2, prob, levels=levels, extend="both")

    ax.set_xlim([np.amin(rcens), np.amax(rcens)])
    ax.set_ylim([np.amin(rcens), 1e3])

    fill_r2_greaterthan_r1(ax)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_aspect("equal")               
    ax.set_title(tlab)

def fill_r2_greaterthan_r1(ax):
  ''' shade over in grey area of plot where
  yaxis radii >= xaxis radii'''

  xlims = ax.get_xlim()
  ylims = ax.get_ylim()

  r1 = np.logspace(np.log10(xlims[0]), np.log10(xlims[1]), 500)

  ax.fill_between(r1, np.full(r1.shape, ylims[1]), r1,
                  step="post", color="grey")