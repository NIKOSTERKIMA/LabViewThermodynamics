#include "lvt.h"

#include <cstdint>
#include <cstdio>
#include <vector>

int main() {
  const double region_x[] = {0.0, 0.50, 0.0, 0.50};
  const double region_y[] = {0.0, 0.0, 0.50, 0.50};
  const double region_width[] = {0.50, 0.50, 0.50, 0.50};
  const double region_height[] = {0.50, 0.50, 0.50, 0.50};
  const double region_lambda[] = {400.0, 50.0, 0.04, 50.0};
  const double region_c_vol[] = {3.45e6, 3.5e6, 1.0e6, 3.5e6};
  const double region_initial_temperature[] = {20.0, 200.0, 20.0, 300.0};
  const uint8_t region_heat_source[] = {0, 1, 0, 1};

  char error[256] = {};
  constexpr uint32_t region_count = static_cast<uint32_t>(sizeof(region_x) / sizeof(region_x[0]));
  constexpr double spatial_step = 0.02;
  constexpr double time_step = 0.05;
  constexpr double simulation_time = 60.0;
  constexpr uint32_t sample_col_count = 20;
  constexpr uint32_t sample_row_count = 20;
  constexpr uint32_t sample_count = sample_col_count * sample_row_count;

  std::vector<double> x(sample_count);
  std::vector<double> y(sample_count);
  std::vector<double> temperatures(sample_count);

  if (lvt_compute_temperatures_2d(region_x, region_y, region_width, region_height, region_lambda, region_c_vol,
                                  region_initial_temperature, region_heat_source, region_count, spatial_step, time_step,
                                  simulation_time, x.data(), y.data(), temperatures.data(), sample_col_count,
                                  sample_row_count, error, sizeof(error)) != 0) {
    std::fprintf(stderr, "compute failed: %s\n", error);
    return 1;
  }

  std::fprintf(stderr, "%s\n", error);

  std::printf("x_m,y_m,temperature\n");
  for (uint32_t i = 0; i < sample_count; i++) {
    std::printf("%.6f,%.6f,%.6f\n", x[i], y[i], temperatures[i]);
  }

  return 0;
}
