'''
----- CLEO -----
File: build_buii.py
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
Script compiles CLEO and generates thermodynamics 
and gridfiles for a constant 1-D rainshaft
'''

import os
import sys
import numpy as np
import random 
from pathlib import Path

path2CLEO = sys.argv[1]
path2build = sys.argv[2]
orig_configfile = sys.argv[3]
tmp_configfile = sys.argv[4]

sys.path.append(path2CLEO)  # for imports from pySD package
from pySD.gbxboundariesbinary_src import read_gbxboundaries as rgrid
from pySD.gbxboundariesbinary_src import create_gbxboundaries as cgrid
from pySD.thermobinary_src import thermogen
from pySD.thermobinary_src import create_thermodynamics as cthermo
from pySD.thermobinary_src import read_thermodynamics as rthermo

### ---------------------------------------------------------------- ###
### ----------------------- INPUT PARAMETERS ----------------------- ###
### ---------------------------------------------------------------- ###
# list of names of executables to compile
executables = ["buii_coalbure"]

### --- essential paths and filenames --- ###
# path and filenames for creating initial conditions
constsfile    = path2CLEO+"/libs/cleoconstants.hpp"
binpath       = path2build+"/bin/"
sharepath     = path2build+"/share/"
tmppath       = path2build+"/tmp/"
gridfile      = sharepath+"buii_dimlessGBxboundaries.dat"
thermofile    = sharepath+"buii_dimlessthermo.dat"
### --- plotting initialisation figures --- ###
isfigures   = [True, True] # booleans for [making, saving] initialisation figures
savefigpath = path2build+"/bin/" # directory for saving init figures

### --- settings for 1-D gridbox boundaries --- ###
zgrid       = [0, 20, 20]        # evenly spaced zhalf coords [zmin, zmax, zdelta] [m]
xgrid       = np.array([0, 20])  # array of xhalf coords [m]
ygrid       = np.array([0, 20])  # array of yhalf coords [m]

### --- settings for 1-D Thermodynamics --- ###
PRESS0      = 101315                # [Pa]
TEMP0       = 297.9                 # [K]
qvap0       = 0.016                 # [Kg/Kg]
Zbase       = 800                   # [m]
TEMPlapses  = [9.8, 6.5]            # -dT/dz [K/km]
qvaplapses  = [2.97, "saturated"]   # -dvap/dz [g/Kg km^-1]
qcond       = 0.0                   # [Kg/Kg]
WMAX        = 0.0                   # [m/s]
Wlength     = 0.0                   # [m] use constant W (Wlength=0.0), or sinusoidal 1-D profile below cloud base

### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###


### ---------------------------------------------------------------- ###
### ------------------- BINARY FILES GENERATION--------------------- ###
### ---------------------------------------------------------------- ###
### --- ensure build, share and bin directories exist --- ###
if path2CLEO == path2build:
  raise ValueError("build directory cannot be CLEO")
else:
  Path(path2build).mkdir(exist_ok=True) 
  Path(sharepath).mkdir(exist_ok=True) 
  Path(binpath).mkdir(exist_ok=True)
  Path(tmppath).mkdir(exist_ok=True)

### ----- copy root config file to tmp_config file ----- ###
os.system('cp '+orig_configfile+" "+tmp_configfile)

### ----- write gridbox boundaries binary ----- ###
cgrid.write_gridboxboundaries_binary(gridfile, zgrid, xgrid, ygrid, constsfile)
rgrid.print_domain_info(constsfile, gridfile)

### ----- write thermodynamics binaries ----- ###
thermodyngen = thermogen.ConstHydrostaticLapseRates(tmp_configfile, constsfile,
                                                    PRESS0, TEMP0, qvap0,
                                                    Zbase, TEMPlapses,
                                                    qvaplapses, qcond,
                                                    WMAX, None, None,
                                                    Wlength)
cthermo.write_thermodynamics_binary(thermofile, thermodyngen, 
                                    tmp_configfile, constsfile,
                                    gridfile)

### ----- show (and save) plots of binary file data ----- ###
if isfigures[0]:
  if isfigures[1]:
    Path(savefigpath).mkdir(exist_ok=True) 
  rgrid.plot_gridboxboundaries(constsfile, gridfile,
                               savefigpath, isfigures[1])
  rthermo.plot_thermodynamics(constsfile, tmp_configfile, gridfile,
                              thermofile, savefigpath, isfigures[1])
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###

### ---------------------------------------------------------------- ###
### ------------------------- COMPILE CLEO ------------------------- ###
### ---------------------------------------------------------------- ###
# 2. compile and the run model
os.chdir(path2build)
os.system('pwd')
for exec in executables:
  os.system('make -j 64 '+exec)
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###