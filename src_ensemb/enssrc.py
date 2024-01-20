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
import pySD.sdmout_src.ensembzarr as enszarr

def write_ensemble_distributions(ensembdataset,
                                 ensembsetupfile,
                                 setupfile,
                                 datasets):
  
  refset = datasets[0] # reference dataset
  for dataset in datasets:
    enszarr.check_dataset_for_ensemb(dataset, refset)
  
  ensemble_numconc_distrib(ensembdataset, ensembsetupfile,
                           setupfile, datasets)
  
  write_distrib_to_zarr()

def ensemble_numconc_distrib(ensembdataset, ensembsetupfile,
                             setupfile, datasets):
  
  for dataset in datasets:
    numconc_distrib()
  
  ensemble_distrib()

def numconc_distrib():
  print("now calc numconc distrib for run")

def ensemble_distrib():

  print("now take mean and std distrib for ensemb")

def write_distrib_to_zarr():

  print("now write to zarr")