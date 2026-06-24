#pragma once

#include <lvt.h>

#include <cstddef>
#include <vector>

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

struct SpatialGrid2D {
  std::vector<Cell> cells;
  size_t col_count = 0;
  size_t row_count = 0;
  double spatial_step = 0.0;
  double origin_x = 0.0;
  double origin_y = 0.0;

  SpatialGrid2D(double spatial_step, size_t col_count, size_t row_count, double origin_x, double origin_y);

  Cell &get_cell(size_t row, size_t col);
  const Cell &get_cell(size_t row, size_t col) const;
};
