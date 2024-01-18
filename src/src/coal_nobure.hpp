/*
 * ----- CLEO -----
 * File: coal_nobure.hpp
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
 * coalescence, breakup or rebound events
 * in SDM analagous to to Shima et al. 2009.
 * CoalNoBuRe struct satisfies PairEnactX
 * concept used in Collisions struct
 */

#ifndef COALNOBURE_HPP
#define COALNOBURE_HPP

#include <functional>
#include <concepts>

#include <Kokkos_Core.hpp>

#include "./breakup.hpp"
#include "./breakup_nfrags.hpp"
#include "./coalbure_flag.hpp"
#include "./coalescence.hpp"
#include "./collisions.hpp"
#include "./collisionkinetics.hpp"
#include "./microphysicalprocess.hpp"
#include "./superdrop.hpp"

template <NFragments NFrags, CoalBuReFlag Flag>
struct DoCoalNoBuRe
/* ie. DoCoalescenceBreakupRebound without BreakupRebound */
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
  bool if_coalesce(const unsigned long long gamma,
                   const double phi,
                   Superdrop &drop1,
                   Superdrop &drop2) const;
  /* function enacts coalescence if flag = 1 -> coalescence.
  Otherwise -> rebound. */

public:
  DoCoalNoBuRe(const NFrags nfrags, const Flag flag)
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
CoalNoBuRe(const unsigned int interval,
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

  const DoCoalNoBuRe<NFrags, Flag> coalnobure(nfrags, coalbure_flag);
  const DoCollisions<Probability, DoCoalNoBuRe<NFrags, Flag>> colls(DELT,
                                                                  collprob,
                                                                  coalnobure);

  return ConstTstepMicrophysics(interval, colls);
}

template <NFragments NFrags, CoalBuReFlag Flag>
KOKKOS_FUNCTION bool
DoCoalNoBuRe<NFrags, Flag>::operator()(Superdrop &drop1, Superdrop &drop2,
                                     const double prob, const double phi) const
/* this operator is used as an "adaptor" for
using DoCoalNoBuRe for collision - coalescence,
without breakup or rebound as a function in DoCollisions
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
    return coalesce_breakup_or_rebound(gamma, phi, drop1, drop2);
  }

  return 0;
}

template <NFragments NFrags, CoalBuReFlag Flag>
KOKKOS_FUNCTION bool
DoCoalNoBuRe<NFrags, Flag>::
    if_coalesce(const unsigned long long gamma,
                                const double phi,
                                Superdrop &drop1,
                                Superdrop &drop2) const
/*  function enacts coalescence if flag = 1 -> coalescence.
Otherwise -> rebound. */
{
  const auto flag = coalbure_flag(phi, drop1, drop2);

  if (flag == 1) // coalescence
  {
    const auto is_null = coal.coalesce_superdroplet_pair(gamma, drop1, drop2);
    return is_null;
  }
  else
  {
    return 0;
  }
}

#endif // COALNOBURE_HPP
