/*
 * ----- CLEO -----
 * File: config_collisions.hpp
 * Project: src
 * Created Date: Thursday 11th January 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Thursday 11th January 2024
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * CLEO super-droplet model (SDM) for 1-D
 * rainshaft investigation of droplet breakup
 */

#ifndef CONFIG_COLLISIONS_HPP 
#define CONFIG_COLLISIONS_HPP 

#include <concepts>

#include "./coalbreakup.hpp"
#include "./coalrebound.hpp"

#include "initialise/config.hpp"
#include "initialise/timesteps.hpp"

#include "superdrops/breakup_nfrags.hpp"
#include "superdrops/breakup.hpp"
#include "superdrops/coalbure_flag.hpp"
#include "superdrops/coalbure.hpp"
#include "superdrops/coalescence.hpp"
#include "superdrops/collisionprobs/longhydroprob.hpp"
#include "superdrops/terminalvelocity.hpp"

struct ConfigCollisions_CoalBuRe
{
  inline MicrophysicalProcess auto
  operator()(const Config &config, const Timesteps &tsteps) const
  {
    const double bu_nfrags = 8;

    const PairProbability auto collprob = LongHydroProb(1.0);
    const NFragments auto nfrags = ConstNFrags(bu_nfrags);
    const CoalBuReFlag auto coalbure_flag = TSCoalBuReFlag{};
    const MicrophysicalProcess auto colls = CoalBuRe(tsteps.get_collstep(),
                                                     &step2realtime,
                                                     collprob,
                                                     nfrags,
                                                     coalbure_flag);

    return colls;
  }
};

struct ConfigCollisions_CoalBreakup
{
  inline MicrophysicalProcess auto
  operator()(const Config &config, const Timesteps &tsteps) const
  {
    const double bu_nfrags = 8;

    const PairProbability auto collprob = LongHydroProb(1.0);
    const NFragments auto nfrags = ConstNFrags(bu_nfrags);
    const CoalBuReFlag auto coalbure_flag = TSCoalBuReFlag{};
    const MicrophysicalProcess auto colls = CoalBreakup(tsteps.get_collstep(),
                                                        &step2realtime,
                                                        collprob,
                                                        nfrags,
                                                        coalbure_flag);

    return colls;
  }
};

struct ConfigCollisions_CoalRebound
{
  inline MicrophysicalProcess auto
  operator()(const Config &config, const Timesteps &tsteps) const
  {
    const PairProbability auto collprob = LongHydroProb(1.0);
    const CoalBuReFlag auto coalbure_flag = TSCoalBuReFlag{};
    const MicrophysicalProcess auto colls = CoalRebound(tsteps.get_collstep(),
                                                        &step2realtime,
                                                        collprob,
                                                        coalbure_flag);

    return colls;
  }
};

struct ConfigCollisions_CoalOnly
{
  inline MicrophysicalProcess auto
  operator()(const Config &config, const Timesteps &tsteps) const
  {
    const PairProbability auto coalprob = LongHydroProb(1.0);
    const MicrophysicalProcess auto colls = CollCoal(tsteps.get_collstep(),
                                                    &step2realtime,
                                                      coalprob);

    return colls;
  }
};

#endif // CONFIG_COLLISIONS_HPP 