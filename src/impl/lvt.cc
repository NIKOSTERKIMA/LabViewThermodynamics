#include "lvt.h"
#include "backend/lvt_naive.h"
#include "lvt_internal.h"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <expected>
#include <optional>
#include <ratio>
#include <string>
#include <utility>

using Grid1DResult = std::expected<SpatialGrid1D, std::string>;
using Grid2DResult = std::expected<SpatialGrid2D, std::string>;

Region1D assemble_region_1d(const Region1DInputs &inputs, uint32_t index) {
  return Region1D{.x = inputs.x[index],
                  .length = inputs.length[index],
                  .lambda = inputs.lambda[index],
                  .c_vol = inputs.c_vol[index],
                  .initial_temperature = inputs.initial_temperature[index],
                  .heat_source = inputs.heat_source[index] != 0};
}

Region2D assemble_region_2d(const Region2DInputs &inputs, uint32_t index) {
  return Region2D{.x = inputs.x[index],
                  .y = inputs.y[index],
                  .width = inputs.width[index],
                  .height = inputs.height[index],
                  .lambda = inputs.lambda[index],
                  .c_vol = inputs.c_vol[index],
                  .initial_temperature = inputs.initial_temperature[index],
                  .heat_source = inputs.heat_source[index] != 0};
}

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

static std::optional<std::string> validate_region_inputs_1d(const Region1DInputs &inputs) {
  if (inputs.count == 0) {
    return std::string("region_count must be positive");
  }

  if (inputs.x == nullptr || inputs.length == nullptr || inputs.lambda == nullptr || inputs.c_vol == nullptr ||
      inputs.initial_temperature == nullptr || inputs.heat_source == nullptr) {
    return std::string("region input arrays must be non-null");
  }

  return std::nullopt;
}

static std::optional<std::string> validate_region_inputs_2d(const Region2DInputs &inputs) {
  if (inputs.count == 0) {
    return std::string("region_count must be positive");
  }

  if (inputs.x == nullptr || inputs.y == nullptr || inputs.width == nullptr || inputs.height == nullptr ||
      inputs.lambda == nullptr || inputs.c_vol == nullptr || inputs.initial_temperature == nullptr ||
      inputs.heat_source == nullptr) {
    return std::string("region input arrays must be non-null");
  }

  return std::nullopt;
}

static std::optional<std::string> validate_compute_args(double spatial_step, double time_step, double simulation_time) {
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

static std::optional<std::string> validate_radiation_args(double emissivity, double plate_thickness) {
  if (emissivity < 0.0 || emissivity > 1.0) {
    return std::string("emissivity must be between 0.0 and 1.0");
  }

  if (plate_thickness <= 0.0) {
    return std::string("plate_thickness must be positive");
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

static Grid1DResult prepare_grid_1d(const Region1DInputs &inputs, double spatial_step) {
  const Region1D first = assemble_region_1d(inputs, 0);
  double x_min = first.x;
  double x_max = first.x + first.length;
  for (uint32_t i = 1; i < inputs.count; i++) {
    const Region1D region = assemble_region_1d(inputs, i);
    if (region.x < x_min) {
      x_min = region.x;
    }
    if (region.x + region.length > x_max) {
      x_max = region.x + region.length;
    }
  }

  for (uint32_t i = 0; i < inputs.count; i++) {
    const Region1D region = assemble_region_1d(inputs, i);
    if (region.length < spatial_step) {
      char msg[128];
      std::snprintf(msg, sizeof(msg), "region %u length %.6f is smaller than spatial_step %.6f", i, region.length,
                    spatial_step);
      return std::unexpected(std::string(msg));
    }
  }

  for (uint32_t i = 0; i < inputs.count; i++) {
    for (uint32_t j = i + 1; j < inputs.count; j++) {
      if (overlap_1d(assemble_region_1d(inputs, i), assemble_region_1d(inputs, j), spatial_step)) {
        char msg[128];
        std::snprintf(msg, sizeof(msg), "regions %u and %u overlap by at least spatial_step %.6f", i, j, spatial_step);
        return std::unexpected(std::string(msg));
      }
    }
  }

  const size_t cell_count = grid_extent(x_max - x_min, spatial_step);
  SpatialGrid1D grid(spatial_step, cell_count, x_min);

  for (uint32_t i = 0; i < inputs.count; i++) {
    const Region1D region = assemble_region_1d(inputs, i);
    for (size_t j = cell_index(region.x, x_min, spatial_step);
         j < cell_index(region.x + region.length, x_min, spatial_step); j++) {
      Cell &cell = grid.get_cell(j);
      cell.temperature = region.initial_temperature;
      cell.lambda = region.lambda;
      cell.c_vol = region.c_vol;
      cell.heat_source = region.heat_source;
    }
  }

  return grid;
}

static Grid2DResult prepare_grid_2d(const Region2DInputs &inputs, double spatial_step) {
  const Region2D first = assemble_region_2d(inputs, 0);
  double x_min = first.x;
  double x_max = first.x + first.width;
  double y_min = first.y;
  double y_max = first.y + first.height;
  for (uint32_t i = 1; i < inputs.count; i++) {
    const Region2D region = assemble_region_2d(inputs, i);
    if (region.x < x_min) {
      x_min = region.x;
    }
    if (region.x + region.width > x_max) {
      x_max = region.x + region.width;
    }
    if (region.y < y_min) {
      y_min = region.y;
    }
    if (region.y + region.height > y_max) {
      y_max = region.y + region.height;
    }
  }

  for (uint32_t i = 0; i < inputs.count; i++) {
    const Region2D region = assemble_region_2d(inputs, i);
    if (region.width < spatial_step || region.height < spatial_step) {
      char msg[128];
      std::snprintf(msg, sizeof(msg), "region %u width %.6f or height %.6f is smaller than spatial_step %.6f", i,
                    region.width, region.height, spatial_step);
      return std::unexpected(std::string(msg));
    }
  }

  for (uint32_t i = 0; i < inputs.count; i++) {
    for (uint32_t j = i + 1; j < inputs.count; j++) {
      if (overlap_2d(assemble_region_2d(inputs, i), assemble_region_2d(inputs, j), spatial_step)) {
        char msg[128];
        std::snprintf(msg, sizeof(msg), "regions %u and %u overlap by at least spatial_step %.6f", i, j, spatial_step);
        return std::unexpected(std::string(msg));
      }
    }
  }

  const size_t col_count = grid_extent(x_max - x_min, spatial_step);
  const size_t row_count = grid_extent(y_max - y_min, spatial_step);
  SpatialGrid2D grid(spatial_step, col_count, row_count, x_min, y_min);

  for (uint32_t i = 0; i < inputs.count; i++) {
    const Region2D region = assemble_region_2d(inputs, i);
    for (size_t row = cell_index(region.y, y_min, spatial_step);
         row < cell_index(region.y + region.height, y_min, spatial_step); row++) {
      for (size_t col = cell_index(region.x, x_min, spatial_step);
           col < cell_index(region.x + region.width, x_min, spatial_step); col++) {
        Cell &cell = grid.get_cell(row, col);
        cell.temperature = region.initial_temperature;
        cell.lambda = region.lambda;
        cell.c_vol = region.c_vol;
        cell.heat_source = region.heat_source;
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

LVT_API uint8_t LVT_CALL lvt_compute_temperatures_1d(const double *region_x, const double *region_length,
                                                     const double *region_lambda, const double *region_c_vol,
                                                     const double *region_initial_temperature,
                                                     const uint8_t *region_heat_source, uint32_t region_count,
                                                     double spatial_step, double time_step, double simulation_time,
                                                     double *x_out, double *temperatures_out, uint32_t point_count,
                                                     char *error_message, uint32_t error_message_length) {
  const Region1DInputs inputs{.x = region_x,
                              .length = region_length,
                              .lambda = region_lambda,
                              .c_vol = region_c_vol,
                              .initial_temperature = region_initial_temperature,
                              .heat_source = region_heat_source,
                              .count = region_count};

  if (const std::optional<std::string> validation_error = validate_region_inputs_1d(inputs)) {
    write_error(error_message, error_message_length, validation_error->c_str());
    return 1;
  }

  if (const std::optional<std::string> validation_error =
          validate_compute_args(spatial_step, time_step, simulation_time)) {
    write_error(error_message, error_message_length, validation_error->c_str());
    return 1;
  }

  if (const std::optional<std::string> output_error =
          validate_output_buffers_1d(x_out, temperatures_out, point_count)) {
    write_error(error_message, error_message_length, output_error->c_str());
    return 1;
  }

  const Grid1DResult grid_result = prepare_grid_1d(inputs, spatial_step);
  if (!grid_result) {
    write_error(error_message, error_message_length, grid_result.error().c_str());
    return 1;
  }

  SpatialGrid1D grid = std::move(*grid_result);
  const auto backend_start = std::chrono::steady_clock::now();
  if (const std::optional<std::string> backend_error =
          lvt_compute_temperatures_naive_1d(grid, time_step, simulation_time)) {
    write_error(error_message, error_message_length, backend_error->c_str());
    return 1;
  }
  const double backend_ms =
      std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - backend_start).count();
  sample_grid_1d(grid, x_out, temperatures_out, point_count);
  write_success_timing(error_message, error_message_length, backend_ms);
  return 0;
}

LVT_API uint8_t LVT_CALL lvt_compute_temperatures_2d(
    const double *region_x, const double *region_y, const double *region_width, const double *region_height,
    const double *region_lambda, const double *region_c_vol, const double *region_initial_temperature,
    const uint8_t *region_heat_source, uint32_t region_count, double spatial_step, double time_step,
    double simulation_time, double *x_out, double *y_out, double *temperatures_out, uint32_t sample_col_count,
    uint32_t sample_row_count, char *error_message, uint32_t error_message_length) {
  const Region2DInputs inputs{.x = region_x,
                              .y = region_y,
                              .width = region_width,
                              .height = region_height,
                              .lambda = region_lambda,
                              .c_vol = region_c_vol,
                              .initial_temperature = region_initial_temperature,
                              .heat_source = region_heat_source,
                              .count = region_count};

  if (const std::optional<std::string> validation_error = validate_region_inputs_2d(inputs)) {
    write_error(error_message, error_message_length, validation_error->c_str());
    return 1;
  }

  if (const std::optional<std::string> validation_error =
          validate_compute_args(spatial_step, time_step, simulation_time)) {
    write_error(error_message, error_message_length, validation_error->c_str());
    return 1;
  }

  if (const std::optional<std::string> output_error =
          validate_output_buffers_2d(x_out, y_out, temperatures_out, sample_col_count, sample_row_count)) {
    write_error(error_message, error_message_length, output_error->c_str());
    return 1;
  }

  const Grid2DResult grid_result = prepare_grid_2d(inputs, spatial_step);
  if (!grid_result) {
    write_error(error_message, error_message_length, grid_result.error().c_str());
    return 1;
  }

  SpatialGrid2D grid = std::move(*grid_result);
  const auto backend_start = std::chrono::steady_clock::now();
  if (const std::optional<std::string> backend_error =
          lvt_compute_temperatures_naive_2d(grid, time_step, simulation_time)) {
    write_error(error_message, error_message_length, backend_error->c_str());
    return 1;
  }
  const double backend_ms =
      std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - backend_start).count();
  sample_grid_2d(grid, x_out, y_out, temperatures_out, sample_col_count, sample_row_count);
  write_success_timing(error_message, error_message_length, backend_ms);
  return 0;
}

LVT_API uint8_t LVT_CALL lvt_compute_temperatures_2d_with_radiation(
    const double *region_x, const double *region_y, const double *region_width, const double *region_height,
    const double *region_lambda, const double *region_c_vol, const double *region_initial_temperature,
    const uint8_t *region_heat_source, uint32_t region_count, double spatial_step, double time_step,
    double simulation_time, double ambient_temperature, double emissivity, double plate_thickness, double *x_out,
    double *y_out, double *temperatures_out, uint32_t sample_col_count, uint32_t sample_row_count, char *error_message,
    uint32_t error_message_length) {
  const Region2DInputs inputs{.x = region_x,
                              .y = region_y,
                              .width = region_width,
                              .height = region_height,
                              .lambda = region_lambda,
                              .c_vol = region_c_vol,
                              .initial_temperature = region_initial_temperature,
                              .heat_source = region_heat_source,
                              .count = region_count};

  if (const std::optional<std::string> validation_error = validate_region_inputs_2d(inputs)) {
    write_error(error_message, error_message_length, validation_error->c_str());
    return 1;
  }

  if (const std::optional<std::string> validation_error =
          validate_compute_args(spatial_step, time_step, simulation_time)) {
    write_error(error_message, error_message_length, validation_error->c_str());
    return 1;
  }

  if (const std::optional<std::string> radiation_error = validate_radiation_args(emissivity, plate_thickness)) {
    write_error(error_message, error_message_length, radiation_error->c_str());
    return 1;
  }

  if (const std::optional<std::string> output_error =
          validate_output_buffers_2d(x_out, y_out, temperatures_out, sample_col_count, sample_row_count)) {
    write_error(error_message, error_message_length, output_error->c_str());
    return 1;
  }

  const Grid2DResult grid_result = prepare_grid_2d(inputs, spatial_step);
  if (!grid_result) {
    write_error(error_message, error_message_length, grid_result.error().c_str());
    return 1;
  }

  SpatialGrid2D grid = std::move(*grid_result);
  grid.radiation = RadiationParams{
      .ambient_temperature = ambient_temperature, .emissivity = emissivity, .plate_thickness = plate_thickness};

  const auto backend_start = std::chrono::steady_clock::now();
  if (const std::optional<std::string> backend_error =
          lvt_compute_temperatures_naive_2d(grid, time_step, simulation_time)) {
    write_error(error_message, error_message_length, backend_error->c_str());
    return 1;
  }
  const double backend_ms =
      std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - backend_start).count();
  sample_grid_2d(grid, x_out, y_out, temperatures_out, sample_col_count, sample_row_count);
  write_success_timing(error_message, error_message_length, backend_ms);
  return 0;
}
