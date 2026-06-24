#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

struct Region1DInputs {
  const double *x;
  const double *length;
  const double *lambda;
  const double *c_vol;
  const double *initial_temperature;
  const uint8_t *heat_source;
  uint32_t count;
};

struct Region1D {
  double x;
  double length;
  double lambda;
  double c_vol;
  double initial_temperature;
  bool heat_source;
};

Region1D assemble_region_1d(const Region1DInputs &inputs, uint32_t index);

struct Region2DInputs {
  const double *x;
  const double *y;
  const double *width;
  const double *height;
  const double *lambda;
  const double *c_vol;
  const double *initial_temperature;
  const uint8_t *heat_source;
  uint32_t count;
};

struct Region2D {
  double x;
  double y;
  double width;
  double height;
  double lambda;
  double c_vol;
  double initial_temperature;
  bool heat_source;
};

Region2D assemble_region_2d(const Region2DInputs &inputs, uint32_t index);

struct Cell {
  double temperature = 0.0;
  double lambda = 0.0;
  double c_vol = 1.0;      // avoid division by zero by default
  bool heat_source = true; // avoid redundant updates by default
};

struct SpatialGrid1D {
  std::vector<Cell> cells;
  double spatial_step = 0.0;
  double origin_x = 0.0;

  SpatialGrid1D(double spatial_step, size_t cell_count, double origin_x);

  Cell &get_cell(size_t index);
  const Cell &get_cell(size_t index) const;
};

struct RadiationParams {
  double ambient_temperature;
  double emissivity;
  double plate_thickness;
};

struct SpatialGrid2D {
  std::vector<Cell> cells;
  size_t col_count = 0;
  size_t row_count = 0;
  double spatial_step = 0.0;
  double origin_x = 0.0;
  double origin_y = 0.0;
  std::optional<RadiationParams> radiation;

  SpatialGrid2D(double spatial_step, size_t col_count, size_t row_count, double origin_x, double origin_y);

  Cell &get_cell(size_t row, size_t col);
  const Cell &get_cell(size_t row, size_t col) const;
};
