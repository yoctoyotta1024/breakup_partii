/*
 * ----- CLEO -----
 * File: coalbreakup.hpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Monday 15th January 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * functionality to enact collision-
 * coalescence or breakup events
 * in SDM analagous to to Shima et al. 2009.
 * CoalBreakup struct satisfies PairEnactX
 * concept used in Collisions struct
 */

#ifndef COALBREAKUP_HPP
#define COALBREAKUP_HPP

#include <functional>
#include <concepts>

#include <Kokkos_Core.hpp>

#include "superdrops/breakup.hpp"
#include "superdrops/breakup_nfrags.hpp"
#include "superdrops/coalbure_flag.hpp"
#include "superdrops/coalescence.hpp"
#include "superdrops/collisions.hpp"
#include "superdrops/collisionkinetics.hpp"
#include "superdrops/microphysicalprocess.hpp"
#include "superdrops/superdrop.hpp"

template <NFragments NFrags, CoalBuReFlag Flag>
struct DoCoalBreakup
/* ie. DoCoalescenceBreakup */
{
private:
  DoCoalescence coal;
  DoBreakup<NFrags> bu;
  Flag coalbure_flag;

  KOKKOS_FUNCTION
  unsigned long long collision_gamma(const unsigned long long xi1,
                                     const unsigned long long xi2,
                                     const double prob,
                                     const double phi) const
  /* calculates value of gamma factor in Monte Carlo
  collision as in Shima et al. 2009 given probability of
  collision. Note: probability is probability of
  collision *NOT* collision-coalescence! */
  {
    return coal.coalescence_gamma(xi1, xi2, prob, phi);
  }

  KOKKOS_FUNCTION
  bool coalesce_or_breakup(const unsigned long long gamma,
                           const double phi,
                           Superdrop &drop1,
                           Superdrop &drop2) const;
  /* function enacts rebound or coalescence or breakup
  depending on value of flag. If flag = 2 -> breakup.
  Otherwise -> coalescence. */

public:
  DoCoalBreakup(const NFrags nfrags, const Flag flag)
      : bu(nfrags), coalbure_flag(flag) {}

  KOKKOS_INLINE_FUNCTION
  bool operator()(Superdrop &drop1, Superdrop &drop2,
                  const double prob, const double phi) const;
  /* this operator is used as an "adaptor" for
  using DoCoalBreakup for collision - coalescence,
  or breakup as a function in DoCollisions
  that satistfies the PairEnactX concept */
};

template <PairProbability Probability,
          NFragments NFrags,
          CoalBuReFlag Flag>
inline MicrophysicalProcess auto
CoalBreakup(const unsigned int interval,
            const std::function<double(unsigned int)> int2realtime,
            const Probability collprob,
            const NFrags nfrags,
            const Flag coalbure_flag)
/* constructs Microphysical Process for collision-
coalscence, breakup or rebound of superdroplets with
a constant timestep 'interval' and probability
of collision determined by 'collprob' */
{
  const auto DELT = double{int2realtime(interval)};

  const DoCoalBreakup<NFrags, Flag>
      coalbreakup(nfrags, coalbure_flag);
  const DoCollisions<Probability, DoCoalBreakup<NFrags, Flag>>
      colls(DELT, collprob, coalbreakup);

  return ConstTstepMicrophysics(interval, colls);
}

template <NFragments NFrags, CoalBuReFlag Flag>
KOKKOS_FUNCTION bool
DoCoalBreakup<NFrags, Flag>::operator()(Superdrop &drop1,
                                        Superdrop &drop2,
                                        const double prob,
                                        const double phi) const
/* this operator is used as an "adaptor" for
using DoCoalBreakup for collision - coalescence,
breakup or rebound as a function in DoCollisions
that satistfies the PairEnactX concept */
{
  /* 1. calculate gamma factor for collision  */
  const auto xi1 = drop1.get_xi();
  const auto xi2 = drop2.get_xi();
  const auto gamma = collision_gamma(xi1, xi2, prob, phi);

  /* 2. enact collision between pair
  of superdroplets if gamma is not zero */
  if (gamma != 0)
  {
    return coalesce_or_breakup(gamma, phi, drop1, drop2);
  }

  return 0;
}

template <NFragments NFrags, CoalBuReFlag Flag>
KOKKOS_FUNCTION bool
DoCoalBreakup<NFrags, Flag>::
    coalesce_or_breakup(const unsigned long long gamma,
                        const double phi,
                        Superdrop &drop1,
                        Superdrop &drop2) const
/* function enacts rebound or coalescence or breakup
depending on value of flag. If flag = 2 -> breakup.
Otherwise -> coalescence.
So if flag=0 or 1 -> coalescence (flag for rebound
(flag=0) enacts coalescence instead) */
{
  const auto flag = coalbure_flag(phi, drop1, drop2);

  if (flag == 2) // breakup
  {
    bu.breakup_superdroplet_pair(drop1, drop2);
    return 0;
  }
  else // coalescence
  {
    auto is_null = coal.coalesce_superdroplet_pair(gamma, drop1, drop2);
    return is_null;
  }
}

#endif // COALBREAKUP_HPP
