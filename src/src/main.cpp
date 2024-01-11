/*
 * ----- CLEO -----
 * File: main.cpp
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

#include "main_supplement.hpp"
#include "coalbure_supplement.hpp"

inline MicrophysicalProcess auto
create_microphysics(const Config &config, const Timesteps &tsteps)
{
  const MicrophysicalProcess auto colls = config_collisions();

  const MicrophysicalProcess auto cond = Condensation(tsteps.get_condstep(),
                                                      config.doAlterThermo,
                                                      config.cond_iters,
                                                      &step2dimlesstime,
                                                      config.cond_rtol,
                                                      config.cond_atol,
                                                      config.cond_SUBTSTEP,
                                                      &realtime2dimless);

  return cond >> colls;
}

inline auto create_sdm(const Config &config,
                       const Timesteps &tsteps,
                       FSStore &store)
{
  const auto couplstep = (unsigned int)tsteps.get_couplstep();
  const GridboxMaps auto gbxmaps(create_gbxmaps(config));
  const MicrophysicalProcess auto microphys(create_microphysics(config, tsteps));
  const Motion<CartesianMaps> auto movesupers(create_motion(tsteps.get_motionstep()));
  const Observer auto obs(create_observer(config, tsteps, store));

  return SDMMethods(couplstep, gbxmaps,
                    microphys, movesupers, obs);
}

int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    throw std::invalid_argument("configuration file(s) not specified");
  }

  Kokkos::Timer kokkostimer;

  /* Read input parameters from configuration file(s) */
  const std::string_view config_filename(argv[1]); // path to configuration file
  const Config config(config_filename);
  const Timesteps tsteps(config); // timesteps for model (e.g. coupling and end time)

  /* Create zarr store for writing output to storage */
  FSStore fsstore(config.zarrbasedir);

  /* Initial conditions for CLEO run */
  const InitialConditions auto initconds = create_initconds(config);

  /* Initialise Kokkos parallel environment */
  Kokkos::initialize(argc, argv);
  {
    /* CLEO Super-Droplet Model (excluding coupled dynamics solver) */
    const SDMMethods sdm(create_sdm(config, tsteps, fsstore));

    /* Solver of dynamics coupled to CLEO SDM */
    CoupledDynamics auto coupldyn(
        create_coupldyn(config, sdm.gbxmaps,
                        tsteps.get_couplstep(),
                        tsteps.get_t_end()));

    /* coupling between coupldyn and SDM */
    const CouplingComms<FromFileDynamics> auto comms = FromFileComms{};
    
    /* Run CLEO (SDM coupled to dynamics solver) */
    const RunCLEO runcleo(sdm, coupldyn, comms);
    runcleo(initconds, tsteps.get_t_end());
  }
  Kokkos::finalize();

  const auto ttot = double{kokkostimer.seconds()};
  std::cout << "-----\n Total Program Duration: "
            << ttot << "s \n-----\n";

  return 0;
}