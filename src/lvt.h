#ifndef LVT_H
#define LVT_H

#include <cstdint>

#ifdef _WIN32
#ifdef LVT_BUILD
#define LVT_API __declspec(dllexport)
#else
#define LVT_API __declspec(dllimport)
#endif
#define LVT_CALL __cdecl
#else
#ifdef LVT_BUILD
#define LVT_API __attribute__((visibility("default")))
#else
#define LVT_API
#endif
#define LVT_CALL
#endif

#ifdef __cplusplus
extern "C" {
#endif

#pragma pack(push, 1)

typedef struct {
  double x;
  double length;
  double lambda;
  double c_vol;
  double initial_temperature;
  uint8_t heat_source;
} Region1D;

typedef struct {
  double x;
  double y;
  double width;
  double height;
  double lambda;
  double c_vol;
  double initial_temperature;
  uint8_t heat_source;
} Region2D;

#pragma pack(pop)

// Shared LabView facing dll interface. Backend dispatch argument chooses concrete computational engine/implementation.
// 0 - naive, 1 - SIMD, 2 - naive CUDA, 3 - optimized CUDA
// Returns 0 if successful, 1 if error.
LVT_API uint8_t LVT_CALL lvt_compute_temperatures_1d(Region1D *regions, uint32_t region_count, double spatial_step,
                                                     double time_step, double simulation_time, double *x_out,
                                                     double *temperatures_out, uint32_t point_count,
                                                     char *error_message, uint32_t error_message_length,
                                                     int8_t backend_dispatch);
LVT_API uint8_t LVT_CALL lvt_compute_temperatures_2d(Region2D *regions, uint32_t region_count, double spatial_step,
                                                     double time_step, double simulation_time, double *x_out,
                                                     double *y_out, double *temperatures_out, uint32_t sample_col_count,
                                                     uint32_t sample_row_count, char *error_message,
                                                     uint32_t error_message_length, int8_t backend_dispatch);

#ifdef __cplusplus
}
#endif

#endif
