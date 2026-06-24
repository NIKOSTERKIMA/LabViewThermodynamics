#pragma once

#include <impl/lvt_internal.h>

#include <optional>
#include <string>

std::optional<std::string> lvt_compute_temperatures_naive_1d(SpatialGrid1D &grid, double time_step,
                                                             double simulation_time);
std::optional<std::string> lvt_compute_temperatures_naive_2d(SpatialGrid2D &grid, double time_step,
                                                             double simulation_time);
