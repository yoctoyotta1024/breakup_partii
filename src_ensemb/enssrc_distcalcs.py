'''
----- CLEO -----
File: ensemb_buii.py
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
functions for writing out
ensemble data in zarr format
'''

import numpy as np

import sys
from pathlib import Path
sys.path.append("/home/m/m300950/CLEO/")  # path2CLEO for imports from pySD package
from pySD.sdmout_src import *             # pyzarr, pysetuptxt & pygbxsdat

def get_log10redgs(rspan, nbins):
  ''' returns edges of log10(r) bins for 'nbins'
  evenly spaced from rspan[0] to rspan[1]'''

  log10redgs = np.linspace(np.log10(rspan[0]),
                             np.log10(rspan[1]),
                             nbins+1)  
  
  return log10redgs # log10(r /microns)

def log10r_histogram(log10redgs, log10radius, wghts):
  ''' returns (weighted) frequency in each log(r) bin
  with edges given by log10r_hedges'''

   # get (weighted) number frequency in each bin
  return np.histogram(log10radius, bins=log10redgs,
                      weights=wghts, density=None)[0]

def get_redges_rcens(log10redgs):

  redges = 10**log10redgs                                             # radius edges of bins
  rcens = (10**(log10redgs[1:]) + 10**(log10redgs[:-1])) / 2        # radius centres of bins

  return redges, rcens # [microns]

def get_domainvol(setupfile, gridfile):

  consts = pysetuptxt.get_consts(setupfile, isprint=True) 
  gbxs = pygbxsdat.get_gridboxes(gridfile, consts["COORD0"], isprint=True)
    
  return gbxs["domainvol"] #[m^3]

def watermass(radius):
  ''' mass [g] as if entire droplet is
  water with density RHO_L
  (neglecting volume of solute) 
  given radius in microns '''
  
  RHO_L = 998.203 #[Kg/m^3]
  rcubed = radius * radius * radius / 1e18 # [m]

  return 4.0 / 3.0 * np.pi * RHO_L * rcubed * 1000 #[g]

def log10r_distrib(rspan, nbins, radius, wghts, perlog10r=False):
  ''' get distribution of data with weights 'wghts' against
  log10(r). Uses np.histogram to get frequency of a particular
  value of data that falls into bins evenly spaced in log10(r) '''

  # create edges of log10(r) histogram bins (evenly spaced in log10(r))
  log10redgs = get_log10redgs(rspan, nbins) 
  log10r = np.log10(radius)
  
  hist = log10r_histogram(log10redgs, log10r, wghts)
  if perlog10r == True: # histogram frequency / delta_log10(r)
    log10rwdths = log10redgs[1:]- log10redgs[:-1]                 # ln10(r) histogram bin widths
    hist = hist/log10rwdths 
 
  redges, rcens = get_redges_rcens(log10redgs)

  return hist, redges, rcens # units of hedgs and hcens = units of rspan (usually [microns])

def numconc_distrib(dataset, log10redgs, gbxidx, vol):
  '''calculate the real droplet number concentration
  for a gridbox with volume 'vol' and index 'gbxidx'.
  If gbxidx=="domain", all superdroplets in dataset
  are used, so 'vol' should be domain volume) '''

  numconc = [] # array dims [time, nbins]
  if gbxidx == "domain":
    
    radius = pyzarr.get_rawdata4raggedkey(dataset, "radius")
    xi = pyzarr.get_rawdata4raggedkey(dataset, "xi")

    log10r = np.log10(radius)
    wghts = xi / vol / 1e6          # real droplets [/cm^3]
    for t in range(len(radius)): # for each timestep
      hist = log10r_histogram(log10redgs, log10r[t], wghts[t])
      numconc.append(hist)
  
  numconc = np.asarray(numconc) # array dims [time, nbins]
  
  return numconc # units: [cm^-3]

def watermass_distrib(dataset, log10redgs, gbxidx, vol):
  '''calculate the real droplet mass concentration
  for a gridbox with volume 'vol' and index 'gbxidx'.
  If gbxidx=="domain", all superdroplets in dataset
  are used, so 'vol' should be domain volume) '''

  massconc = [] # array dims [time, nbins]
  if gbxidx == "domain":
    
    radius = pyzarr.get_rawdata4raggedkey(dataset, "radius")
    xi = pyzarr.get_rawdata4raggedkey(dataset, "xi")
    wghts = xi * watermass(radius) / vol # real droplets mass [g/m^3]

    log10r = np.log10(radius)
    for t in range(len(radius)): # for each timestep
      hist = log10r_histogram(log10redgs, log10r[t], wghts[t])
      massconc.append(hist)
  
  massconc = np.asarray(massconc) # array dims [time, nbins]
  
  return massconc # units: [g/m^3]