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
import awkward as ak

import sys
from pathlib import Path
sys.path.append("/home/m/m300950/CLEO/")  # path2CLEO for imports from pySD package
from pySD.sdmout_src import *             # pyzarr, pysetuptxt & pygbxsdat
import pySD.sdmout_src.ensembzarr as enszarr

def write_ensemble_domaindists(ensembdataset, ensembsetupfile,
                               setupfile, gridfile, datasets,
                               distparams):
  ''' write number, mass and mass^2 droplet 
  distributions for entire domain to ensemble
  zarr. parametrs for distirubtions given by
  distparams={nbins, rspan} dictionary'''
  refset = datasets[0] # reference dataset
  for dataset in datasets:
    enszarr.check_dataset_for_ensemb(dataset, refset)
  
  ensemble_domainnumconc_distrib(ensembdataset, ensembsetupfile,
                                 setupfile, gridfile, datasets,
                                 distparams)
  
  write_domaindistrib_to_zarr()

def get_domainvol(setupfile, gridfile):

  consts = pysetuptxt.get_consts(setupfile, isprint=True) 
  gbxs = pygbxsdat.get_gridboxes(gridfile, consts["COORD0"], isprint=True)
    
  return gbxs["domainvol"]

def ensemble_domainnumconc_distrib(ensembdataset,
                                   ensembsetupfile,
                                   setupfile,
                                   gridfile,
                                   datasets,
                                   distparams):
  ''' take mean of real droplet number
  concentration distributions.
  parametrs for distirubtions given by
  distparams={nbins, rspan} dictionary'''

  for dataset in datasets:
    domainvol = get_domainvol(setupfile, gridfile) 
    numconc_distrib(dataset, domainvol, "domain",
                    distparams["nbins"], distparams["rspan"])
  
  ensemble_distrib()

def log10r_distribution(rspan, nbins, radius, wghts, perlog10r=False):
  ''' get distribution of data with weights 'wghts' against
  log10(r). Uses np.histogram to get frequency of a particular
  value of data that falls into bins evenly spaced in log10(r) '''

  # create edges of log10(r) histogram bins (evenly spaced in log10(r))
  log10r_hedgs = np.linspace(np.log10(rspan[0]), np.log10(rspan[1]), nbins+1)  
  
  # get (weighted) number frequency in each bin
  hist = np.histogram(np.log10(radius), bins=log10r_hedgs,
                      weights=wghts, density=None)[0]

  if perlog10r == True: # histogram frequency / delta_log10(r)
    log10r_hwdths = log10r_hedgs[1:]- log10r_hedgs[:-1]                 # ln10(r) histogram bin widths
    hist = hist/log10r_hwdths

  redges = 10**log10r_hedgs                                             # radius edges of bins
  rcens = (10**(log10r_hedgs[1:]) + 10**(log10r_hedgs[:-1])) / 2        # radius centres of bins

  return hist, redges, rcens # units of hedgs and hcens = units of rspan (usually [microns])

def numconc_distrib(dataset, vol, gbxidx, nbins, rspan):
  '''calculate the real droplet number concentration
  for a gridbox with volume 'vol' and index 'gbxidx'.
  If gbxidx=="domain", all superdroplets in dataset
  are used, so 'vol' should be domain volume) '''

  

  print("---- WIP -----")
  numconc = [] # array dims [time, nbins]
  if gbxidx == "domain":
    xi = pyzarr.get_rawdata4raggedkey(dataset, "xi")

    for 
  print(xi, len(xi), len(xi[0]))
  # print(ak.shape(xi))
  
  
  print("---- WIP -----")

def ensemble_distrib():

  print("now take mean and std distrib for ensemb")

def write_domaindistrib_to_zarr():

  print("now write to zarr")