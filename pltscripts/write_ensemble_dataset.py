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
import zarr

path2pySD = "/home/m/m300950/CLEO/" # TODO: remove
sys.path.append(path2pySD)  # for imports from pySD package
from pySD.sdmout_src import pyzarr
from pySD import editconfigfile

def write_ensemble_info(ensembsetupfile, setupfile, datasets):
  header = "// ----------------------------- //" +\
            "\n// --------- ensembsetupfile --------- //" +\
            "\n// ----------------------------- //\n"
  footer = "// ----------------------------- //"

  datasets_str = "\ndatasets in ensemble: \n     "+\
                "\n     ".join(datasets) + "\n"
        
  setup_str = "\nsetup copied from: \n     "+setupfile+"\n"

  with open(ensembsetupfile, "a") as file:
    file.write(header)
    file.write(setup_str)
    file.write(datasets_str)
    file.write(footer)

def write_ensemb_setupfile(ensembsetupfile, setupfile, datasets):
  
  os.system('cp '+setupfile+" "+ensembsetupfile)
  params = {
    "initsupers_filename" : "[ensemble, see below]",
    "setuptxt" : "[ensemble, see below]", 
    "zarrbasedir" : "[ensemble, see below]",
    "stats_filename" : "[ensemble, see below]"
    }
  editconfigfile.edit_config_params(ensembsetupfile, params)

  write_ensemble_info(ensembsetupfile, setupfile, datasets)

def get_datasets_for_ensemb(datasets, refset):
  '''returns opened datasets after checking that time
  in each is consistent with the refset'''
  
  dss = []
  time = pyzarr.get_rawdataset(refset)["time"].values
  for dataset in datasets:
    d = pyzarr.get_rawdataset(dataset)
  
    if np.any(d["time"].values != time):
      raise ValueError("data for time in datasets must be the same")

    dss.append(d)

  return dss

def write_time_to_ensembzarr(ensembdataset, dataset):
  ''' create or replace time group in ensembdataset
  with time group copied from dataset '''

  src = zarr.open(dataset)["time"]
  dest = zarr.open(ensembdataset)
  zarr.copy(src, dest, log=sys.stdout, if_exists='replace')

def write_ensemb_zarr(ensembdataset, vars4ensemb, datasets):

  refset = datasets[0] # reference dataset
  dss = get_datasets_for_ensemb(datasets, refset)
  
  write_time_to_ensembzarr(ensembdataset, refset)
  
def write_ensemble_dataset(ensembdataset, ensembsetupfile,
                           vars4ensemb, setupfile, datasets):

  write_ensemb_setupfile(ensembsetupfile, setupfile, datasets)

  write_ensemb_zarr(ensembdataset, vars4ensemb, datasets)
  time = pyzarr.get_rawdataset(ensembdataset)["time"].values
  print("time: ", time)