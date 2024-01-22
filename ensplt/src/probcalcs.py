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

def watermass(diam):
  ''' mass as if water droplet given diameter [microns]'''

  radius = diam / 2e6 # convert diameter [microns] to radius [m]
  rho = 998.203 # density [kg/m^3]

  return 4.0/3 * np.pi * rho * radius**3 #[kg]

def simmel_terminalv(radius):
  ''' returns terminal velocity [m/s] of droplets
  given array of their radii [microns] according to
  simmel et al. 2002 '''

  diam = 2.0*radius # [microns]

  alpha = np.full(diam.shape, 9.17)
  beta = np.full(diam.shape, 0.0)

  alpha = np.where(diam < 3477.84, 17.32, alpha)
  beta = np.where(diam < 3477.84, 1.0/6, beta)

  alpha = np.where(diam < 1511.64, 49.62, alpha)
  beta = np.where(diam < 1511.64, 1.0/3, beta)

  alpha = np.where(diam < 134.43, 4579.5, alpha)
  beta = np.where(diam < 134.43, 2.0/3, beta)

  mass = watermass(diam) * 1000 # [grams]
  terminalv = alpha * mass**beta

  return terminalv # [m/s]

def hydrodyanmic_kernel(rr1, rr2, terminalv, eff=1.0):
  ''' returns kernel K(drop1, drop2) (proportional to probability)
  for a pair of droplets colliding according to the hydrodynamic,
  i.e. gravitational, collision kernel.
  Probability is given by prob_jk = K(drop1, drop2) * delta_t/delta_vol,
  (see Shima 2009 eqn 3) where the kernel,
  K(drop1, drop2) := eff * pi * (r1 + r2)^2 * |v1−v2|,
  given the efficiency factor eff = eff(drop1, drop2) '''
  
  sumrsqrd = (rr1 + rr2)*(rr1 + rr2)
  vdiff = terminalv(rr1) - terminalv(rr2)
  hydro_kernel = np.pi * eff * sumrsqrd * vdiff

  return hydro_kernel

def relative_collision_probability(rcens, numconc):
  ''' calculate probability of collision using
  Long's hydrodynamic kernel according
  to Simmel et al. 2002'''

  rr1, rr2 = np.meshgrid(rcens, rcens)
  kernel = hydrodyanmic_kernel(rr1, rr2, simmel_terminalv, eff=1.0)
  numdens = np.outer(numconc, numconc).T
  
  relprob = kernel * numdens # proportional to probability

  # remove data where rr2 > rr1
  relprob = np.where(rr1 <= rr2, np.nan, relprob)

  return rr1, rr2, relprob

def log10_collprob(datalab, rcens, numconc):
  ''' calculate probability of collision using
  Long's hydrodynamic kernel according
  to Simmel et al. 2002'''

  rr1, rr2, relprob = relative_collision_probability(rcens, numconc)
  
  # log10(relprob)
  relprob = np.where(relprob == 0.0, np.nan, relprob)
  relprob = np.log10(relprob)

  return rr1, rr2, relprob

def relative_outcome_probability(datalab, rcens, numconc):
  ''' calculate probability of collision-outcome
  using Long's hydrodynamic kernel according
  to Simmel et al. 2002'''

  rr1, rr2, relprob = relative_collision_probability(rcens, numconc)
  
  if datalab == "coalbure":
    outcome = relative_outcome_probability_coalbure(relprob)
  elif datalab == "coalbu":
    outcome = relative_outcome_probability_coalbu(relprob)
  elif datalab == "coalre":
    outcome = relative_outcome_probability_coalre(relprob)
  
  return rr1, rr2, outcome

def relative_outcome_probability_coalbure(relprob):
  
  print("TODO: coalbure")
 
  outcome = {
    "coal" : np.full(relprob.shape, 0.0), # TODO 
    "bu" : np.full(relprob.shape, 0.0), # TODO 
    "re" : np.full(relprob.shape, 0.0), # TODO 
  }

  return outcome

def relative_outcome_probability_coalbu(relprob):
  
  print("TODO: coalbu")
  
  outcome = {
    "coal" : np.full(relprob.shape, 1.0) * relprob, # TODO 
    "bu" : np.full(relprob.shape, 1.0) * relprob, # TODO 
    "re" : np.full(relprob.shape, 0.0),
  }

  return outcome

def relative_outcome_probability_coalre(relprob):

  print("TODO: coalre")

  outcome = {
    "coal" : np.full(relprob.shape, 1.0) * relprob, # TODO 
    "bu" : np.full(relprob.shape, 0.0),
    "re" : np.full(relprob.shape, 1.0) * relprob, # TODO 
  }

  return outcome

def relative_collcoal_probability(datalab, rcens, numconc):
  ''' calculate probability of collision-coalescence
  using Long's hydrodynamic kernel according
  to Simmel et al. 2002'''

  rr1, rr2, outcome = relative_outcome_probability(datalab, rcens, numconc)

  return rr1, rr2, outcome["coal"]

def relative_collbreakup_probability(datalab, rcens, numconc):
  ''' calculate probability of collision-coalescence
  using Long's hydrodynamic kernel according
  to Simmel et al. 2002'''

  rr1, rr2, outcome = relative_outcome_probability(datalab, rcens, numconc)

  return rr1, rr2, outcome["bu"]

def relative_collrebound_probability(datalab, rcens, numconc):
  ''' calculate probability of collision-coalescence
  using Long's hydrodynamic kernel according
  to Simmel et al. 2002'''

  rr1, rr2, outcome = relative_outcome_probability(datalab, rcens, numconc)

  return rr1, rr2, outcome["re"]