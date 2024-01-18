/*
 * ----- CLEO -----
 * File: main_coalbreakup.cpp
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
#include "config_collisions.hpp"

int main(int argc, char *argv[])
{
  return main_supplement(argc, argv, ConfigCollisions_CoalBreakup{});
}