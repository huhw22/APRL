#ifndef TRAJECTORY_RADIATION_CONFIGURATION_H
#define TRAJECTORY_RADIATION_CONFIGURATION_H

#include <string>

#include "radiation_types.h"

namespace radiation
{
  RadiationConfig loadConfiguration(const std::string& filename);
}

#endif
