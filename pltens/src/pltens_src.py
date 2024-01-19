'''
----- CLEO -----
File: pltens_src.py
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
Script plots some data from ensemble
dataset of buii 1-D rainshaft output
'''

import numpy as np
import matplotlib.pyplot as plt

def savefig(fig, savename, show=True):
  
  fig.tight_layout()
  
  fig.savefig(savename, dpi=400, bbox_inches="tight",
              facecolor='w', format="png")
  print("Figure .png saved as: "+savename)

  if show:
    plt.show()

def plot_gbxmassmoments(axs, zgbx, time, massmoms, color="k"):

  line0 = axs[0].plot(time.mins, massmoms.nsupers[:,0,0,zgbx], color=color)
  axs[1].plot(time.mins, massmoms.mom0[:,0,0,zgbx], color=color)
  axs[2].plot(time.mins, massmoms.mom1[:,0,0,zgbx], color=color)
  axs[3].plot(time.mins, massmoms.mom2[:,0,0,zgbx], color=color)
  axs[4].plot(time.mins, massmoms.effmass[:,0,0,zgbx], color=color)

  axs[0].set_ylabel("number of\nsuperdroplets")
  axs[1].set_ylabel("$\u03BB^{m}_{0}$, number\nof droplets")
  axs[2].set_ylabel("$\u03BB^{m}_{1}$, droplet\nmass /g")
  axs[3].set_ylabel("$\u03BB^{m}_{2}$\n~reflectivity /g$^2$")
  ylab4 = "mean effective\ndroplet mass,\n<$\u03BB^{m}_{2}$/$\u03BB^{m}_{1}>$ /g"
  axs[4].set_ylabel(ylab4)

  for ax in [axs[1], axs[2]]:
    ax.set_yscale("log")
  axs[-1].set_xlabel("time /min")

  return line0[0]
