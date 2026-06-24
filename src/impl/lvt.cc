#include "lvt.h"
#include "backend/lvt_naive.h"
#include "lvt_internal.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <expected>
#include <optional>
#include <string>

using Grid1DResult = std::expected<SpatialGrid1D, std::string>;
using Grid2DResult = std::expected<SpatialGrid2D, std::string>;

SpatialGrid1D::SpatialGrid1D(double spatial_step, size_t cell_count, double origin_x)
    : cells(cell_count), spatial_step(spatial_step), origin_x(origin_x) {}

Cell &SpatialGrid1D::get_cell(size_t index) { return cells[index]; }

const Cell &SpatialGrid1D::get_cell(size_t index) const { return cells[index]; }

SpatialGrid2D::SpatialGrid2D(double spatial_step, size_t col_count, size_t row_count, double origin_x, double origin_y)
    : cells(col_count * row_count), col_count(col_count), row_count(row_count), spatial_step(spatial_step),
      origin_x(origin_x), origin_y(origin_y) {}

Cell &SpatialGrid2D::get_cell(size_t row, size_t col) { return cells[row * col_count + col]; }

const Cell &SpatialGrid2D::get_cell(size_t row, size_t col) const { return cells[row * col_count + col]; }

static void write_error(char *buf, uint32_t len, const char *msg) {
  if (buf != nullptr && len > 0)
    std::snprintf(buf, len, "%s", msg);
}

static void write_success_timing(char *buf, uint32_t len, double elapsed_ms) {
  if (buf != nullptr && len > 0)
    std::snprintf(buf, len, "ran successful in %.3f ms", elapsed_ms);
}

static std::optional<std::string> validate_compute_args(const void *regions, uint32_t region_count, double spatial_step,
                                                        double time_step, double simulation_time) {

  if (region_count == 0) {
    return std::string("region_count must be positive");
  }

  if (regions == nullptr) {
    return std::string("regions must be non-null");
  }

  if (spatial_step <= 0.0) {
    return std::string("spatial_step must be positive");
  }

  if (time_step <= 0.0) {
    return std::string("time_step must be positive");
  }

  if (simulation_time < 0.0) {
    return std::string("simulation_time must be non-negative");
  }

  return std::nullopt;
}

static std::optional<std::string> dispatch_compute_1d(SpatialGrid1D &grid, double time_step, double simulation_time,
                                                      int8_t backend_dispatch) {
  switch (backend_dispatch) {
  case 0:
    return lvt_compute_temperatures_naive_1d(grid, time_step, simulation_time);
  case 1:
    return std::string("SIMD backend is not implemented");
  case 2:
    return std::string("multithreaded SIMD backend is not implemented");
  case 3:
    return std::string("naive CUDA backend is not implemented");
  case 4:
    return std::string("optimized CUDA backend is not implemented");
  default:
    return std::string("unsupported backend");
  }
}

static std::optional<std::string> dispatch_compute_2d(SpatialGrid2D &grid, double time_step, double simulation_time,
                                                      int8_t backend_dispatch) {
  switch (backend_dispatch) {
  case 0:
    return lvt_compute_temperatures_naive_2d(grid, time_step, simulation_time);
  case 1:
    return std::string("SIMD backend is not implemented");
  case 2:
    return std::string("multithreaded SIMD backend is not implemented");
  case 3:
    return std::string("naive CUDA backend is not implemented");
  case 4:
    return std::string("optimized CUDA backend is not implemented");
  default:
    return std::string("unsupported backend");
  }
}

static std::optional<std::string> validate_output_buffers_1d(double *x_out, double *temperatures_out,
                                                             uint32_t point_count) {
  if (point_count == 0) {
    return std::string("point_count must be positive");
  }

  if (x_out == nullptr || temperatures_out == nullptr) {
    return std::string("x_out and temperatures_out must be non-null");
  }

  return std::nullopt;
}

static std::optional<std::string> validate_output_buffers_2d(double *x_out, double *y_out, double *temperatures_out,
                                                             uint32_t sample_col_count, uint32_t sample_row_count) {
  if (sample_col_count == 0 || sample_row_count == 0) {
    return std::string("sample_col_count and sample_row_count must be positive");
  }

  if (x_out == nullptr || y_out == nullptr || temperatures_out == nullptr) {
    return std::string("x_out, y_out and temperatures_out must be non-null");
  }

  return std::nullopt;
}

static size_t cell_index(double coord, double origin, double step) {
  return static_cast<size_t>((coord - origin) / step);
}

static size_t grid_extent(double span, double step) { return static_cast<size_t>(std::ceil(span / step)); }

static bool overlap_1d(const Region1D &a, const Region1D &b, double step) {
  const double overlap = std::fmax(0.0, std::fmin(a.x + a.length, b.x + b.length) - std::fmax(a.x, b.x));
  return overlap >= step;
}

static bool overlap_2d(const Region2D &a, const Region2D &b, double step) {
  const double overlap_x = std::fmax(0.0, std::fmin(a.x + a.width, b.x + b.width) - std::fmax(a.x, b.x));
  const double overlap_y = std::fmax(0.0, std::fmin(a.y + a.height, b.y + b.height) - std::fmax(a.y, b.y));
  return overlap_x >= step && overlap_y >= step;
}

static Grid1DResult prepare_grid_1d(const Region1D *regions, uint32_t region_count, double spatial_step) {
  double x_min = regions[0].x;
  double x_max = regions[0].x + regions[0].length;
  for (uint32_t i = 1; i < region_count; i++) {
    if (regions[i].x < x_min) {
      x_min = regions[i].x;
    }
    if (regions[i].x + regions[i].length > x_max) {
      x_max = regions[i].x + regions[i].length;
    }
  }

  for (uint32_t i = 0; i < region_count; i++) {
    if (regions[i].length < spatial_step) {
      char msg[128];
      std::snprintf(msg, sizeof(msg), "region %u length %.6f is smaller than spatial_step %.6f", i, regions[i].length,
                    spatial_step);
      return std::unexpected(std::string(msg));
    }
  }

  for (uint32_t i = 0; i < region_count; i++) {
    for (uint32_t j = i + 1; j < region_count; j++) {
      if (overlap_1d(regions[i], regions[j], spatial_step)) {
        char msg[128];
        std::snprintf(msg, sizeof(msg), "regions %u and %u overlap by at least spatial_step %.6f", i, j, spatial_step);
        return std::unexpected(std::string(msg));
      }
    }
  }

  const size_t cell_count = grid_extent(x_max - x_min, spatial_step);
  SpatialGrid1D grid(spatial_step, cell_count, x_min);

  for (uint32_t i = 0; i < region_count; i++) {
    for (size_t j = cell_index(regions[i].x, x_min, spatial_step);
         j < cell_index(regions[i].x + regions[i].length, x_min, spatial_step); j++) {
      Cell &cell = grid.get_cell(j);
      cell.temperature = regions[i].initial_temperature;
      cell.lambda = regions[i].lambda;
      cell.c_vol = regions[i].c_vol;
      cell.heat_source = regions[i].heat_source != 0;
    }
  }

  return grid;
}

static Grid2DResult prepare_grid_2d(const Region2D *regions, uint32_t region_count, double spatial_step) {
  double x_min = regions[0].x;
  double x_max = regions[0].x + regions[0].width;
  double y_min = regions[0].y;
  double y_max = regions[0].y + regions[0].height;
  for (uint32_t i = 1; i < region_count; i++) {
    if (regions[i].x < x_min) {
      x_min = regions[i].x;
    }
    if (regions[i].x + regions[i].width > x_max) {
      x_max = regions[i].x + regions[i].width;
    }
    if (regions[i].y < y_min) {
      y_min = regions[i].y;
    }
    if (regions[i].y + regions[i].height > y_max) {
      y_max = regions[i].y + regions[i].height;
    }
  }

  for (uint32_t i = 0; i < region_count; i++) {
    if (regions[i].width < spatial_step || regions[i].height < spatial_step) {
      char msg[128];
      std::snprintf(msg, sizeof(msg), "region %u width %.6f or height %.6f is smaller than spatial_step %.6f", i,
                    regions[i].width, regions[i].height, spatial_step);
      return std::unexpected(std::string(msg));
    }
  }

  for (uint32_t i = 0; i < region_count; i++) {
    for (uint32_t j = i + 1; j < region_count; j++) {
      if (overlap_2d(regions[i], regions[j], spatial_step)) {
        char msg[128];
        std::snprintf(msg, sizeof(msg), "regions %u and %u overlap by at least spatial_step %.6f", i, j, spatial_step);
        return std::unexpected(std::string(msg));
      }
    }
  }

  const size_t col_count = grid_extent(x_max - x_min, spatial_step);
  const size_t row_count = grid_extent(y_max - y_min, spatial_step);
  SpatialGrid2D grid(spatial_step, col_count, row_count, x_min, y_min);

  for (uint32_t i = 0; i < region_count; i++) {
    for (size_t row = cell_index(regions[i].y, y_min, spatial_step);
         row < cell_index(regions[i].y + regions[i].height, y_min, spatial_step); row++) {
      for (size_t col = cell_index(regions[i].x, x_min, spatial_step);
           col < cell_index(regions[i].x + regions[i].width, x_min, spatial_step); col++) {
        Cell &cell = grid.get_cell(row, col);
        cell.temperature = regions[i].initial_temperature;
        cell.lambda = regions[i].lambda;
        cell.c_vol = regions[i].c_vol;
        cell.heat_source = regions[i].heat_source != 0;
      }
    }
  }

  return grid;
}

static void sample_grid_1d(const SpatialGrid1D &grid, double *x_out, double *temperatures_out, uint32_t point_count) {
  const double span = static_cast<double>(grid.cells.size()) * grid.spatial_step;
  for (uint32_t i = 0; i < point_count; i++) {
    const double x = grid.origin_x + (static_cast<double>(i) + 0.5) * span / static_cast<double>(point_count);
    x_out[i] = x;

    size_t cell = cell_index(x, grid.origin_x, grid.spatial_step);
    if (cell >= grid.cells.size()) {
      cell = grid.cells.size() - 1;
    }
    temperatures_out[i] = grid.get_cell(cell).temperature;
  }
}

static void sample_grid_2d(const SpatialGrid2D &grid, double *x_out, double *y_out, double *temperatures_out,
                           uint32_t sample_col_count, uint32_t sample_row_count) {
  const double x_span = static_cast<double>(grid.col_count) * grid.spatial_step;
  const double y_span = static_cast<double>(grid.row_count) * grid.spatial_step;

  for (uint32_t row = 0; row < sample_row_count; row++) {
    for (uint32_t col = 0; col < sample_col_count; col++) {
      const uint32_t i = row * sample_col_count + col;
      const double x =
          grid.origin_x + (static_cast<double>(col) + 0.5) * x_span / static_cast<double>(sample_col_count);
      const double y =
          grid.origin_y + (static_cast<double>(row) + 0.5) * y_span / static_cast<double>(sample_row_count);
      x_out[i] = x;
      y_out[i] = y;

      size_t cell_row = cell_index(y, grid.origin_y, grid.spatial_step);
      size_t cell_col = cell_index(x, grid.origin_x, grid.spatial_step);
      if (cell_row >= grid.row_count) {
        cell_row = grid.row_count - 1;
      }
      if (cell_col >= grid.col_count) {
        cell_col = grid.col_count - 1;
      }
      temperatures_out[i] = grid.get_cell(cell_row, cell_col).temperature;
    }
  }
}

LVT_API uint8_t LVT_CALL lvt_compute_temperatures_1d(Region1D *regions, uint32_t region_count, double spatial_step,
                                                     double time_step, double simulation_time, double *x_out,
                                                     double *temperatures_out, uint32_t point_count,
                                                     char *error_message, uint32_t error_message_length,
                                                     int8_t backend_dispatch) {
  if (const std::optional<std::string> validation_error =
          validate_compute_args(regions, region_count, spatial_step, time_step, simulation_time)) {
    write_error(error_message, error_message_length, validation_error->c_str());
    return 1;
  }

  if (const std::optional<std::string> output_error =
          validate_output_buffers_1d(x_out, temperatures_out, point_count)) {
    write_error(error_message, error_message_length, output_error->c_str());
    return 1;
  }

  const Grid1DResult grid_result = prepare_grid_1d(regions, region_count, spatial_step);
  if (!grid_result) {
    write_error(error_message, error_message_length, grid_result.error().c_str());
    return 1;
  }

  SpatialGrid1D grid = std::move(*grid_result);
  const auto backend_start = std::chrono::steady_clock::now();
  if (const std::optional<std::string> backend_error =
          dispatch_compute_1d(grid, time_step, simulation_time, backend_dispatch)) {
    write_error(error_message, error_message_length, backend_error->c_str());
    return 1;
  }
  const double backend_ms =
      std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - backend_start).count();
  sample_grid_1d(grid, x_out, temperatures_out, point_count);
  write_success_timing(error_message, error_message_length, backend_ms);
  return 0;
}

LVT_API uint8_t LVT_CALL lvt_compute_temperatures_2d(Region2D *regions, uint32_t region_count, double spatial_step,
                                                     double time_step, double simulation_time, double *x_out,
                                                     double *y_out, double *temperatures_out, uint32_t sample_col_count,
                                                     uint32_t sample_row_count, char *error_message,
                                                     uint32_t error_message_length, int8_t backend_dispatch) {
  if (const std::optional<std::string> validation_error =
          validate_compute_args(regions, region_count, spatial_step, time_step, simulation_time)) {
    write_error(error_message, error_message_length, validation_error->c_str());
    return 1;
  }

  if (const std::optional<std::string> output_error =
          validate_output_buffers_2d(x_out, y_out, temperatures_out, sample_col_count, sample_row_count)) {
    write_error(error_message, error_message_length, output_error->c_str());
    return 1;
  }

  const Grid2DResult grid_result = prepare_grid_2d(regions, region_count, spatial_step);
  if (!grid_result) {
    write_error(error_message, error_message_length, grid_result.error().c_str());
    return 1;
  }

  SpatialGrid2D grid = std::move(*grid_result);
  const auto backend_start = std::chrono::steady_clock::now();
  if (const std::optional<std::string> backend_error =
          dispatch_compute_2d(grid, time_step, simulation_time, backend_dispatch)) {
    write_error(error_message, error_message_length, backend_error->c_str());
    return 1;
  }
  const double backend_ms =
      std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - backend_start).count();
  sample_grid_2d(grid, x_out, y_out, temperatures_out, sample_col_count, sample_row_count);
  write_success_timing(error_message, error_message_length, backend_ms);
  return 0;
}
