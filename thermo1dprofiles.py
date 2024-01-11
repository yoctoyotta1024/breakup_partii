'''
----- CLEO -----
File: temp.py
Project: breakup_partii
Created Date: Tuesday 9th January 2024
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Wednesday 10th January 2024
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2024 MPI-M, Clara Bayley
-----
'''

import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate

sys.path.append("../CLEO") # path to pySD (same as to CLEO)
from pySD.thermobinary_src import thermogen

GRAVG = 9.8
RGAS_DRY = 287.04 
CP_DRY = 1004.64
Mr_ratio = 0.01801528 / 0.028966216
configfile = "./src/src/buii_config.txt"
constsfile = "../CLEO/libs/cleoconstants.hpp"

def supersaturation(press, temp, qvap, Mr_ratio):
  
  pv = qvap*press/(Mr_ratio + qvap) # vapour pressure
  psat = thermogen.saturation_press(temp)
  
  qsat = Mr_ratio * psat/(press-pv) 
  supersat = qvap/qsat - 1
  
  return supersat

### --- settings for 1-D Thermodynamics --- ###
PRESS0 = 101315         # [Pa]
TEMP0 = 297.9           # [K]
qvap0 = 0.016           # [Kg/Kg]
Zbase = 800             # [m]
TEMPlapses = [9.8, 6.5]  # -dT/dz [K/km]
qvaplapses = [2.97, "saturated"] # -dvap/dz [g/Kg km^-1]
qcond = 0.0             # [Kg/Kg]
WVEL = 0.0              # [m/s]

thermo = thermogen.ConstHydrostaticLapseRates(configfile, constsfile,
                                              PRESS0, TEMP0, qvap0,
                                              Zbase, TEMPlapses,
                                              qvaplapses, qcond,
                                              WVEL, None, None)

zfulls = np.arange(0,2500,100)
temp, press, qvap = thermo.hydrostatic_lapserates_thermo(zfulls)
theta = temp * (PRESS0 / press) ** (RGAS_DRY/CP_DRY)
supersat = supersaturation(press, temp, qvap, Mr_ratio)

fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(8, 8))
axs = axs.flatten()

axs[0].plot(temp, zfulls/1000)
axs[0].set_xlabel("temp / K")

axs[1].plot(press, zfulls/1000)
axs[1].set_xlabel("press / Pa")

axs[2].plot(theta, zfulls/1000)
axs[2].set_xlabel("theta / K")

axs[3].plot(qvap*1000, zfulls/1000)
axs[3].set_xlabel("qvap / g/Kg")

axs[4].plot(supersat, zfulls/1000)
axs[4].set_xlabel("supersaturation")

qcond = np.full(zfulls.shape, qcond)
axs[5].plot(qcond, zfulls/1000)
axs[5].set_xlabel("qcond / g/Kg")

for ax in axs:
  ax.set_ylabel("z /km")

fig.tight_layout()
plt.savefig("./initprofs.png")