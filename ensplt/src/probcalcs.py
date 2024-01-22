'''
----- CLEO -----
File: probcalcs.py
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
calculations for probabilities of collisions
given radii and droplet distributions
'''

import numpy as np

def collision_probability(rcens, numconc):
  ''' calculate probability of collision using
  Long's hydrodynamic kernel according
  to Simmel et al. 2002'''

  rr1, rr2 = np.meshgrid(rcens, rcens)
  
  prob = np.outer(numconc, numconc)

  prob = np.where(rr1 <= rr2, np.nan, prob) / np.nanmax(prob) #normalise and remove data where rr2 > rr1

  return rr1, rr2, prob    