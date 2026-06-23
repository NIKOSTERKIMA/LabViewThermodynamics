#include "lvt.h"

#include <cstdio>
#include <vector>

int main() {
  // 1 m x 1 m plate in four quadrants: copper | heater / insulation | heater
  Region2D regions[] = {
      {.x = 0.0,
       .y = 0.0,
       .width = 0.50,
       .height = 0.50,
       .lambda = 400.0,
       .c_vol = 3.45e6,
       .initial_temperature = 20.0,
       .heat_source = 0},
      {.x = 0.50,
       .y = 0.0,
       .width = 0.50,
       .height = 0.50,
       .lambda = 50.0,
       .c_vol = 3.5e6,
       .initial_temperature = 200.0,
       .heat_source = 1},
      {.x = 0.0,
       .y = 0.50,
       .width = 0.50,
       .height = 0.50,
       .lambda = 0.04,
       .c_vol = 1.0e6,
       .initial_temperature = 20.0,
       .heat_source = 0},
      {.x = 0.50,
       .y = 0.50,
       .width = 0.50,
       .height = 0.50,
       .lambda = 50.0,
       .c_vol = 3.5e6,
       .initial_temperature = 300.0,
       .heat_source = 1},
  };

  char error[256] = {};
  constexpr uint32_t region_count = static_cast<uint32_t>(sizeof(regions) / sizeof(regions[0]));
  constexpr double spatial_step = 0.02;
  constexpr double time_step = 0.05;
  constexpr double simulation_time = 60.0;
  constexpr uint32_t sample_col_count = 20;
  constexpr uint32_t sample_row_count = 20;
  constexpr uint32_t sample_count = sample_col_count * sample_row_count;

  std::vector<double> x(sample_count);
  std::vector<double> y(sample_count);
  std::vector<double> temperatures(sample_count);

  if (lvt_compute_temperatures_2d(regions, region_count, spatial_step, time_step, simulation_time, x.data(), y.data(),
                                  temperatures.data(), sample_col_count, sample_row_count, error, sizeof(error),
                                  0) != 0) {
    std::fprintf(stderr, "compute failed: %s\n", error);
    return 1;
  }

  std::fprintf(stderr, "%s\n", error);

  std::printf("x_m,y_m,temperature_C\n");
  for (uint32_t i = 0; i < sample_count; i++) {
    std::printf("%.6f,%.6f,%.6f\n", x[i], y[i], temperatures[i]);
  }

  return 0;
}
