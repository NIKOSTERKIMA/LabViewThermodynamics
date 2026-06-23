#include "lvt.h"

#include <cstdio>
#include <vector>

int main() {
  Region1D regions[] = {
      {.x = 0.0,
       .length = 0.40,
       .lambda = 400.0,
       .c_vol = 3.45e6,
       .initial_temperature = 20.0,
       .heat_source = 0},
      {.x = 0.40,
       .length = 0.20,
       .lambda = 50.0,
       .c_vol = 3.5e6,
       .initial_temperature = 200.0,
       .heat_source = 1},
      {.x = 0.60,
       .length = 0.25,
       .lambda = 0.04,
       .c_vol = 1.0e6,
       .initial_temperature = 20.0,
       .heat_source = 0},
      {.x = 0.85,
       .length = 0.15,
       .lambda = 50.0,
       .c_vol = 3.5e6,
       .initial_temperature = 300.0,
       .heat_source = 1},
  };

  char error[256] = {};
  constexpr uint32_t region_count = static_cast<uint32_t>(sizeof(regions) / sizeof(regions[0]));
  constexpr double spatial_step = 0.01;
  constexpr double time_step = 0.0005;
  constexpr double simulation_time = 3000.0;
  constexpr uint32_t sample_count = 100;

  std::vector<double> x(sample_count);
  std::vector<double> temperatures(sample_count);

  if (lvt_compute_temperatures_1d(regions, region_count, spatial_step, time_step, simulation_time, x.data(),
                                  temperatures.data(), sample_count, error, sizeof(error), 0) != 0) {
    std::fprintf(stderr, "compute failed: %s\n", error);
    return 1;
  }

  std::fprintf(stderr, "%s\n", error);

  std::printf("x_m,temperature_C\n");
  for (uint32_t i = 0; i < sample_count; i++) {
    std::printf("%.6f,%.6f\n", x[i], temperatures[i]);
  }

  return 0;
}
