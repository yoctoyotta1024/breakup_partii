/*
 * ----- CLEO -----
 * File: main_supplement.hpp
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

#ifndef MAIN_SUPPLEMENT_HPP
#define MAIN_SUPPLEMENT_HPP

#include <iostream>
#include <stdexcept>
#include <string_view>
#include <concepts>
#include <cmath>
#include <array>

#include <Kokkos_Core.hpp>

#include "cartesiandomain/cartesianmaps.hpp"
#include "cartesiandomain/cartesianmotion.hpp"
#include "cartesiandomain/createcartesianmaps.hpp"

#include "coupldyn_fromfile/fromfilecomms.hpp"
#include "coupldyn_fromfile/fromfile_cartesian_dynamics.hpp"

#include "gridboxes/gridboxmaps.hpp"

#include "initialise/config.hpp"
#include "initialise/timesteps.hpp"
#include "initialise/initgbxs_null.hpp"
#include "initialise/initsupers_frombinary.hpp"

#include "observers/massmomentsobs.hpp"
#include "observers/nsupersobs.hpp"
#include "observers/observers.hpp"
#include "observers/printobs.hpp"
#include "observers/timeobs.hpp"
#include "observers/supersattrsobs.hpp"

#include "runcleo/coupleddynamics.hpp"
#include "runcleo/couplingcomms.hpp"
#include "runcleo/initialconditions.hpp"
#include "runcleo/runcleo.hpp"
#include "runcleo/sdmmethods.hpp"

#include "superdrops/motion.hpp"
#include "superdrops/terminalvelocity.hpp"

#include "zarr/fsstore.hpp"
#include "zarr/superdropattrsbuffers.hpp"
#include "zarr/superdropsbuffers.hpp"

inline CoupledDynamics auto
create_coupldyn(const Config &config,
                const CartesianMaps &gbxmaps,
                const unsigned int couplstep,
                const unsigned int t_end)

{
  const auto h_ndims(gbxmaps.ndims_hostcopy());
  const std::array<size_t, 3>
      ndims({h_ndims(0), h_ndims(1), h_ndims(2)});

  const auto nsteps = (unsigned int)(std::ceil(t_end / couplstep) + 1);

  return FromFileDynamics(config, couplstep, ndims, nsteps);
}

inline InitialConditions auto
create_initconds(const Config &config)
{
  const InitSupersFromBinary initsupers(config);
  const InitGbxsNull initgbxs(config);

  return InitConds(initsupers, initgbxs);
}

inline GridboxMaps auto
create_gbxmaps(const Config &config)
{
  const auto gbxmaps = create_cartesian_maps(config.ngbxs,
                                             config.nspacedims,
                                             config.grid_filename);
  return gbxmaps;
}

inline Motion<CartesianMaps> auto
create_motion(const unsigned int motionstep)
{
  const auto terminalv = SimmelTerminalVelocity{};
  
  return CartesianMotion(motionstep,
                         &step2dimlesstime,
                         terminalv);                                                                            
}

inline Observer auto
create_supersattrs_observer(const unsigned int interval,
                            FSStore &store,
                            const int maxchunk)
{
  SuperdropsBuffers auto buffers = SdIdBuffer() >>
                                   XiBuffer() >>
                                   MsolBuffer() >>
                                   RadiusBuffer() >>
                                   Coord3Buffer() >>
                                   SdgbxindexBuffer();
  return SupersAttrsObserver(interval, store, maxchunk, buffers);
}

inline Observer auto
create_observer(const Config &config,
                const Timesteps &tsteps,
                FSStore &store)
{
  const auto obsstep = (unsigned int)tsteps.get_obsstep();
  const auto maxchunk = int{config.maxchunk};

  const Observer auto obs1 = PrintObserver(obsstep * 10, &step2realtime);

  const Observer auto obs2 = TimeObserver(obsstep, store, maxchunk,
                                          &step2dimlesstime);

  const Observer auto obs3 = NsupersObserver(obsstep, store, maxchunk,
                                             config.ngbxs);

  const Observer auto obs4 = MassMomentsObserver(obsstep, store, maxchunk,
                                                 config.ngbxs);

  const Observer auto obs5 = create_supersattrs_observer(obsstep, store,
                                                          maxchunk);

  return obs1 >> obs2 >> obs3 >> obs4 >> obs5;
}

#endif // MAIN_SUPPLEMENT_HPP