#include "lvt.h"

#include <cstdint>
#include <cstdio>
#include <vector>

int main() {
  const double region_x[] = {0.0, 0.40, 0.60, 0.85};
  const double region_length[] = {0.40, 0.20, 0.25, 0.15};
  const double region_lambda[] = {400.0, 50.0, 0.04, 50.0};
  const double region_c_vol[] = {3.45e6, 3.5e6, 1.0e6, 3.5e6};
  const double region_initial_temperature[] = {20.0, 200.0, 20.0, 300.0};
  const uint8_t region_heat_source[] = {0, 1, 0, 1};

  char error[256] = {};
  constexpr uint32_t region_count = static_cast<uint32_t>(sizeof(region_x) / sizeof(region_x[0]));
  constexpr double spatial_step = 0.01;
  constexpr double time_step = 0.0005;
  constexpr double simulation_time = 3000.0;
  constexpr uint32_t sample_count = 100;

  std::vector<double> x(sample_count);
  std::vector<double> temperatures(sample_count);

  if (lvt_compute_temperatures_1d(region_x, region_length, region_lambda, region_c_vol, region_initial_temperature,
                                  region_heat_source, region_count, spatial_step, time_step, simulation_time, x.data(),
                                  temperatures.data(), sample_count, error, sizeof(error)) != 0) {
    std::fprintf(stderr, "compute failed: %s\n", error);
    return 1;
  }

  std::fprintf(stderr, "%s\n", error);

  std::printf("x_m,temperature\n");
  for (uint32_t i = 0; i < sample_count; i++) {
    std::printf("%.6f,%.6f\n", x[i], temperatures[i]);
  }

  return 0;
}
