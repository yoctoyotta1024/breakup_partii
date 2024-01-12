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

def check_dataset_for_ensemb(dataset, refset):
  '''returns opened dataset after checking that
  time is consistent with the refset'''
  
  time = pyzarr.get_rawdataset(refset)["time"].values
  ds = pyzarr.get_rawdataset(dataset)

  if np.any(ds["time"].values != time):
    raise ValueError("data for time in datasets must be the same")

def write_time_to_ensembzarr(ensembdataset, dataset):
  ''' create or replace time group in ensembdataset
  with time group copied from dataset '''

  src = zarr.open(dataset)["time"]
  dest = zarr.open(ensembdataset)
  zarr.copy(src, dest, log=sys.stdout, if_exists='replace')

def ensemble_data(func, datasets, var):
  ''' call func on ensemble of var from datasets'''
  
  data4ensemb = []
  for dataset in datasets:
    v = pyzarr.get_rawdataset(dataset)[var].values
    data4ensemb.append(v)
  data4ensemb = np.asarray(data4ensemb)
  
  return func(data4ensemb)

def  write_meanvar_to_array(arrayname, array, refvar):
  ''' write array to zarr store under 'arrayname'
  using same metadata as refvar '''
  
  z = zarr.open(arrayname, mode='w',
                  shape=refvar.shape,
                  chunks=refvar.chunks,
                  compressor=refvar.compressor,
                  dtype=refvar.dtype, 
                  fill_value=refvar.fill_value,
                  filters=refvar.filters,
                  order=refvar.order)
  z[:] = array

def write_meanvars_to_ensembzarr(ensembdataset, vars4ensemb,
                                 datasets, refset):
  ''' write mean over ensemble of datasets for
    each var in vars4ensemb into ensembdataset zarr'''
  
  for var in vars4ensemb:
    meanname = ensembdataset+"/"+var
    meanvar = ensemble_data(lambda x : np.mean(x, axis=0),
                            datasets, var)
    write_meanvar_to_array(meanname, meanvar, zarr.open(refset)[var])

def write_ensemb_zarr(ensembdataset, vars4ensemb, datasets):

  refset = datasets[0] # reference dataset
  for dataset in datasets:
    check_dataset_for_ensemb(dataset, refset)

  write_time_to_ensembzarr(ensembdataset, refset)
  write_meanvars_to_ensembzarr(ensembdataset, vars4ensemb,
                               datasets, refset)
      
def write_ensemble_dataset(ensembdataset, ensembsetupfile,
                           vars4ensemb, setupfile, datasets):

  write_ensemb_setupfile(ensembsetupfile, setupfile, datasets)
  write_ensemb_zarr(ensembdataset, vars4ensemb, datasets)