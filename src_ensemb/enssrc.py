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
import zarr

import enssrc_distcalcs as distcalcs

import sys
from pathlib import Path
sys.path.append("/home/m/m300950/CLEO/")  # path2CLEO for imports from pySD package
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
  
  log10redgs = distcalcs.get_log10redgs(distparams["rspan"],
                              distparams["nbins"])
  
  redges, rcens = distcalcs.get_redges_rcens(log10redgs)
  write_redges_rcens(ensembdataset, redges, rcens)

  write_ensemble_domainnumconc_distrib(ensembdataset,
                                       ensembsetupfile,
                                       gridfile,
                                       datasets,
                                       log10redgs) 

  write_ensemble_domainwatermass_distrib(ensembdataset,
                                       ensembsetupfile,
                                       gridfile,
                                       datasets,
                                       log10redgs)  
  
  write_ensemble_domainreflectivity_distrib(ensembdataset,
                                            datasets,
                                            log10redgs) 

def write_domaindistrib_to_zarr(ensembdataset, name, meandist, stddist):

  print("TODO: write mean and std of "+name+" to zarr: "+ensembdataset)
  
def write_zarrarray(array, zarrarrayname, shape, 
               chunks=(1250000), dtype='<f8'):
  
  z1 = zarr.open(zarrarrayname, mode='w', shape=shape,
                 chunks=chunks, dtype=dtype, order='C',
                 compressor=None, fill_value=None, 
                 filters=None, zarr_format=2)

  z1[:] = array

def write_zattrs_metadata(arraydims, units=" ", scale_factor=1.0):

  metadata = '{\n"_ARRAY_DIMENSIONS": '+arraydims+',\n'+\
    '"units": '+units+',\n'+\
      '"scale_factor": '+str(scale_factor)+\
        '\n}'
  print(metadata)

def write_redges_rcens(ensembdataset, redges, rcens):
  
  zarrdataset = "/home/m/m300950/breakup_partii/plots/testzarr.zarr/" # TODO 

  sf = 1e6
  arrayname = zarrdataset+"/redges"
  write_zarrarray(redges/sf, arrayname, redges.shape)

  arrayname = zarrdataset+"/rcens"
  write_zarrarray(rcens/sf, arrayname, rcens.shape)

  write_zattrs_metadata('["rcens"]', units="micro m", scale_factor=sf)
  write_zattrs_metadata('["redges"]', units="micro m", scale_factor=sf)

  print("TODO: write redges rcens")

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

def ensemble_distrib_mean_std(ensemb_dists):
  ''' returns mean of ensemble, and its standard 
   deviation for 'nruns' of distributions
  with dims [nruns, time, nbins] '''

  nruns = ensemb_dists.shape[0]
  meandist = np.mean(ensemb_dists, axis=0)
  stddist = np.std(ensemb_dists, axis=0) / np.sqrt(nruns) # sigma/sqrt(N) 

  return meandist, stddist

def write_ensemble_domainnumconc_distrib(ensembdataset,
                                         ensembsetupfile,
                                         gridfile, datasets,
                                         log10redgs):
  ''' write mean and standard deviation real droplet number
  concentration distributions for domain over datasets
  of ensemble to zarr 'ensembdataset' where bins of
  distribution are defined by log10redges [microns] '''
  
  domainvol = distcalcs.get_domainvol(ensembsetupfile, gridfile) # [m^3]
  numconc_dists = calc_dists_for_ensemb(distcalcs.numconc_distrib,
                                        datasets, log10redgs,
                                        ["domain", domainvol])
  
  meandist, stddist = ensemble_distrib_mean_std(numconc_dists)
  
  write_domaindistrib_to_zarr(ensembdataset, "dist_num", meandist, stddist)

def write_ensemble_domainwatermass_distrib(ensembdataset,
                                           ensembsetupfile,
                                           gridfile,
                                           datasets,
                                           log10redgs): 
  ''' write mean and standard deviation of real droplet
  mass (as if water) distribution for domain over
  datasets of ensemble to zarr 'ensembdataset' where bins
  of distribution are defined by log10redges [microns]'''
    
  domainvol = distcalcs.get_domainvol(ensembsetupfile, gridfile) # [m^3]
  watermass_dists = calc_dists_for_ensemb(distcalcs.watermass_distrib,
                                        datasets, log10redgs,
                                        ["domain", domainvol])
  
  meandist, stddist = ensemble_distrib_mean_std(watermass_dists)
  
  write_domaindistrib_to_zarr(ensembdataset, "dist_watermass",
                              meandist, stddist)
  
def write_ensemble_domainreflectivity_distrib(ensembdataset,
                                              datasets,
                                              log10redgs):
  ''' write mean and standard deviation of real droplet
  reflectivity proxy (6th moment of radius) distribution
  for domain over datasets of ensemble to zarr 'ensembdataset'
  where bins of distribution are defined by log10redges [microns] '''

  refproxy_dists = calc_dists_for_ensemb(distcalcs.reflectproxy_distrib,
                                          datasets, log10redgs,
                                          ["domain"])
  
  meandist, stddist = ensemble_distrib_mean_std(refproxy_dists)
  
  write_domaindistrib_to_zarr(ensembdataset, "dist_refproxy",
                              meandist, stddist)