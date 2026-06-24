#include "lvt_naive.h"
#include "lvt_internal.h"

#include <cmath>
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

static constexpr double STEFAN_BOLTZMANN = 5.670374419e-8;

static double radiation_face_delta(double temperature, const RadiationParams &radiation, double spatial_step,
                                   double c_vol, double time_step) {
  const double dx = spatial_step;
  const double volume = dx * dx * radiation.plate_thickness;
  const double face_area = dx * radiation.plate_thickness;
  const double q_rad =
      radiation.emissivity * STEFAN_BOLTZMANN * (std::pow(temperature, 4) - std::pow(radiation.ambient_temperature, 4));
  return -q_rad * face_area * time_step / (c_vol * volume);
}

std::optional<std::string> lvt_compute_temperatures_naive_1d(SpatialGrid1D &grid, double time_step,
                                                             double simulation_time) {
  const size_t cell_count = grid.cells.size();
  std::vector<double> current(cell_count);
  std::vector<double> next(cell_count);

  for (size_t i = 0; i < cell_count; i++) {
    current[i] = grid.get_cell(i).temperature;
  }

  const double dx2 = grid.spatial_step * grid.spatial_step;

  for (double t = 0.0; t < simulation_time; t += time_step) {
    for (size_t i = 0; i < cell_count; i++) {
      const Cell &center = grid.get_cell(i);
      if (center.heat_source) {
        next[i] = current[i];
        continue;
      }

      double flux = 0.0;
      if (i > 0) {
        const Cell &left = grid.get_cell(i - 1);
        flux += 2.0 * left.lambda * center.lambda / (left.lambda + center.lambda) * (current[i - 1] - current[i]);
      }
      if (i + 1 < cell_count) {
        const Cell &right = grid.get_cell(i + 1);
        flux += 2.0 * center.lambda * right.lambda / (center.lambda + right.lambda) * (current[i + 1] - current[i]);
      }

      next[i] = current[i] + flux * time_step / (center.c_vol * dx2);
    }

    current.swap(next);
  }

  for (size_t i = 0; i < cell_count; i++) {
    grid.get_cell(i).temperature = current[i];
  }

  return std::nullopt;
}

std::optional<std::string> lvt_compute_temperatures_naive_2d(SpatialGrid2D &grid, double time_step,
                                                             double simulation_time) {
  const size_t cell_count = grid.cells.size();
  std::vector<double> current(cell_count);
  std::vector<double> next(cell_count);

  for (size_t row = 0; row < grid.row_count; row++) {
    for (size_t col = 0; col < grid.col_count; col++) {
      current[row * grid.col_count + col] = grid.get_cell(row, col).temperature;
    }
  }

  const double dx2 = grid.spatial_step * grid.spatial_step;

  for (double t = 0.0; t < simulation_time; t += time_step) {
    for (size_t row = 0; row < grid.row_count; row++) {
      for (size_t col = 0; col < grid.col_count; col++) {
        const size_t index = row * grid.col_count + col;
        const Cell &center = grid.get_cell(row, col);
        if (center.heat_source) {
          next[index] = current[index];
          continue;
        }

        double flux_x = 0.0;
        if (col > 0) {
          const Cell &left = grid.get_cell(row, col - 1);
          flux_x +=
              2.0 * left.lambda * center.lambda / (left.lambda + center.lambda) * (current[index - 1] - current[index]);
        }
        if (col + 1 < grid.col_count) {
          const Cell &right = grid.get_cell(row, col + 1);
          flux_x += 2.0 * center.lambda * right.lambda / (center.lambda + right.lambda) *
                    (current[index + 1] - current[index]);
        }

        double flux_y = 0.0;
        if (row > 0) {
          const Cell &bottom = grid.get_cell(row - 1, col);
          flux_y += 2.0 * bottom.lambda * center.lambda / (bottom.lambda + center.lambda) *
                    (current[index - grid.col_count] - current[index]);
        }
        if (row + 1 < grid.row_count) {
          const Cell &top = grid.get_cell(row + 1, col);
          flux_y += 2.0 * center.lambda * top.lambda / (center.lambda + top.lambda) *
                    (current[index + grid.col_count] - current[index]);
        }

        double radiation_delta = 0.0;
        if (grid.radiation.has_value()) {
          const RadiationParams &radiation = *grid.radiation;
          if (col == 0) {
            radiation_delta +=
                radiation_face_delta(current[index], radiation, grid.spatial_step, center.c_vol, time_step);
          }
          if (col + 1 == grid.col_count) {
            radiation_delta +=
                radiation_face_delta(current[index], radiation, grid.spatial_step, center.c_vol, time_step);
          }
          if (row == 0) {
            radiation_delta +=
                radiation_face_delta(current[index], radiation, grid.spatial_step, center.c_vol, time_step);
          }
          if (row + 1 == grid.row_count) {
            radiation_delta +=
                radiation_face_delta(current[index], radiation, grid.spatial_step, center.c_vol, time_step);
          }
        }

        next[index] = current[index] + (flux_x + flux_y) * time_step / (center.c_vol * dx2) + radiation_delta;
      }
    }

    current.swap(next);
  }

  for (size_t row = 0; row < grid.row_count; row++) {
    for (size_t col = 0; col < grid.col_count; col++) {
      grid.get_cell(row, col).temperature = current[row * grid.col_count + col];
    }
  }

  return std::nullopt;
}
