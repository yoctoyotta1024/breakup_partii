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

def get_log10redgs(rspan, nbins):
  ''' returns edges of log10(r) bins for 'nbins'
  evenly spaced from rspan[0] to rspan[1]'''

  log10redgs = np.linspace(np.log10(rspan[0]),
                             np.log10(rspan[1]),
                             nbins+1)  
  
  return log10redgs

def log10r_histogram(log10redgs, log10radius, wghts):
  ''' returns (weighted) frequency in each log(r) bin
  with edges given by log10r_hedges'''

   # get (weighted) number frequency in each bin
  return np.histogram(log10radius, bins=log10redgs,
                      weights=wghts, density=None)[0]

def get_redges_rcens(log10redgs):

  redges = 10**log10redgs                                             # radius edges of bins
  rcens = (10**(log10redgs[1:]) + 10**(log10redgs[:-1])) / 2        # radius centres of bins

  return redges, rcens

def get_domainvol(setupfile, gridfile):

  consts = pysetuptxt.get_consts(setupfile, isprint=True) 
  gbxs = pygbxsdat.get_gridboxes(gridfile, consts["COORD0"], isprint=True)
    
  return gbxs["domainvol"]

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
  
  log10redgs = get_log10redgs(distparams["rspan"],
                              distparams["nbins"])
  redges, rcens = get_redges_rcens(log10redgs)

  write_ensemble_domainnumconc_distrib(ensembdataset,
                                             ensembsetupfile,
                                             setupfile, gridfile,
                                             datasets, log10redgs) 
  
  # ensemble_domainmass_distribs(ensembdataset, ensembsetupfile,
  #                             setupfile, gridfile, datasets,
  #                             distparams) 
  # write_domaindistrib_to_zarr("massdistribs")

  write_redges_rcens()

def write_redges_rcens():

  print("now write redges rcens")

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

def calc_dists_for_ensemb(calc_distrib_func, datasets,
                          log10redgs, args):
  ''' take mean of real droplet number
  concentration distributions.
  parametrs for distirubtions given by
  distparams={nbins, rspan} dictionary'''
  
  dists = []
  for dataset in datasets:
    dist = calc_distrib_func(dataset, log10redgs, *args)
    dists.append(dist)
  dists = np.asarray(dists)

  return dists

def write_ensemble_domainnumconc_distrib(ensembdataset,
                                   ensembsetupfile,
                                   setupfile,
                                   gridfile,
                                   datasets,
                                   log10redgs):
  ''' take mean of real droplet number
  concentration distributions.
  parametrs for distirubtions given by
  distparams={nbins, rspan} dictionary'''
  
  domainvol = get_domainvol(setupfile, gridfile) 
  numconc_dists = calc_dists_for_ensemb(numconc_distrib,
                                        datasets, log10redgs,
                                        ["domain", domainvol])
  
  meandist, stddist = ensemble_distrib(numconc_dists)

  import matplotlib.pyplot as plt
  redges, rcens = get_redges_rcens(log10redgs)
  plt.step(redges[:-1], meandist.T[:,::30], where='pre')
  plt.step(redges[:-1], (meandist+stddist).T[:,::30],
           where='pre', linestyle="--")
  plt.step(redges[:-1], (meandist-stddist).T[:,::30],
           where='pre', linestyle="--")
  plt.yscale("log")
  plt.xscale("log")
  plt.savefig("histt_test.png")

  write_domaindistrib_to_zarr("numconc", meandist, stddist)

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
  
  return numconc

def ensemble_distrib(ensemb_dists):
  ''' returns mean of ensemble, and its standard 
   deviation for 'nruns' of distributions
  with dims [nruns, time, nbins] '''

  nruns = ensemb_dists.shape[0]
  meandist = np.mean(ensemb_dists, axis=0)
  stddist = np.std(ensemb_dists, axis=0) / np.sqrt(nruns) # sigma/sqrt(N) 

  return meandist, stddist

def write_domaindistrib_to_zarr(name, meandist, stddist):

  print("now write to zarr: "+name)