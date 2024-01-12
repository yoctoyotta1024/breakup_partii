'''
----- CLEO -----
File: write_ensemble_dataset.py
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
functions to write a new zarr dataset from 
mean of an ensemble of datasets
'''

import os
import sys
import numpy as np

def import_pyzarr(path2CLEO):
  
  sys.path.append(path2CLEO)  # for imports from pySD package
  from pySD import pyzarr, editconfigfile

def write_ensemb_setupfile(ensembsetuptxt, setupfile):
  

  print(ensembsetuptxt, setupfile)
  # os.system('cp '+setupfile+" "+ensembsetuptxt)
  # params = {
  #   "initsupers_filename" : "ensemble (see below)",
  #   "setuptxt" : "ensemble (see below)",
  #   "zarrbasedir" : "ensemble (see below)",
  #   "stats_filename" "ensemble (see below)",
  #   }
  # editconfigfile.edit_config_params(ensembsetuptxt, params)

def write_ensemb_zarr(ensembdataset, vars4ensemb, datasets):

  for dataset in datasets:
    ds = pyzarr.get_rawdataset(dataset)
    print(dataset)

def write_ensemble_dataset(path2CLEO, ensembdataset, ensembsetuptxt,
                           vars4ensemb, setupfile, datasets):

  import_pyzarr(path2CLEO) # TODO: move import away

  write_ensemb_setupfile(ensembsetuptxt, setupfile)

  write_ensemb_zarr(ensembdataset, vars4ensemb, datasets)