/*
 * ----- CLEO -----
 * File: coalrebound.hpp
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
 * coalescence/rebound events
 * in SDM analagous to to Shima et al. 2009.
 * CoalRe struct satisfies PairEnactX
 * concept used in Collisions struct
 */

#ifndef COALREBOUND_HPP
#define COALREBOUND_HPP

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
struct DoCoalRebound
/* ie. DoCoalescenceBreakupRebound without Breakup */
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
  bool coalesce_or_rebound(const unsigned long long gamma,
                           const double phi,
                           Superdrop &drop1,
                           Superdrop &drop2) const;
  /* function enacts coalescence if flag = 1 -> coalescence.
  Otherwise -> rebound. */

public:
  DoCoalRebound(const NFrags nfrags, const Flag flag)
      : bu(nfrags), coalbure_flag(flag) {}

  KOKKOS_INLINE_FUNCTION
  bool operator()(Superdrop &drop1, Superdrop &drop2,
                  const double prob, const double phi) const;
  /* this operator is used as an "adaptor" for
  using DoCoaNoBuRe for collision - coalescence, without
  breakup or rebound as a function in DoCollisions
  that satistfies the PairEnactX concept */
};

template <PairProbability Probability,
          NFragments NFrags,
          CoalBuReFlag Flag>
inline MicrophysicalProcess auto
CoalRebound(const unsigned int interval,
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

  const DoCoalRebound<NFrags, Flag>
      coalrebound(nfrags, coalbure_flag);
  const DoCollisions<Probability, DoCoalRebound<NFrags, Flag>>
      colls(DELT, collprob, coalrebound);

  return ConstTstepMicrophysics(interval, colls);
}

template <NFragments NFrags, CoalBuReFlag Flag>
KOKKOS_FUNCTION bool
DoCoalRebound<NFrags, Flag>::operator()(Superdrop &drop1, Superdrop &drop2,
                                     const double prob, const double phi) const
/* this operator is used as an "adaptor" for
using DoCoalRebound for collision - coalescence/rebound,
without breakup as a function in DoCollisions
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
    return coalesce_or_rebound(gamma, phi, drop1, drop2);
  }

  return 0;
}

template <NFragments NFrags, CoalBuReFlag Flag>
KOKKOS_FUNCTION bool
DoCoalRebound<NFrags, Flag>::
    coalesce_or_rebound(const unsigned long long gamma,
                                const double phi,
                                Superdrop &drop1,
                                Superdrop &drop2) const
/*  function enacts rebound if flag = 0 , otherwise coalescence,
so flag = 1 or 2 -> coalescence (ie. breakup flag (flag=2) 
does coalescence instead of breakup) */
{
  const auto flag = coalbure_flag(phi, drop1, drop2);

  if (flag != 0) // coalescence 
  {
    const auto is_null = coal.coalesce_superdroplet_pair(gamma, drop1, drop2);
    return is_null;
  }
  else
  {
    return 0;
  }
}

#endif // COALREBOUND_HPP
