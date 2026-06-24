#ifndef LVT_H
#define LVT_H

#ifdef __cplusplus
#include <cstdint>
#else
typedef signed char int8_t;
typedef unsigned char uint8_t;
typedef signed short int16_t;
typedef unsigned short uint16_t;
typedef signed int int32_t;
typedef unsigned int uint32_t;
typedef signed long long int64_t;
typedef unsigned long long uint64_t;
#endif

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

// Returns 0 if successful, 1 if error.
LVT_API uint8_t LVT_CALL lvt_compute_temperatures_1d(const double *region_x, const double *region_length,
                                                     const double *region_lambda, const double *region_c_vol,
                                                     const double *region_initial_temperature,
                                                     const uint8_t *region_heat_source, uint32_t region_count,
                                                     double spatial_step, double time_step, double simulation_time,
                                                     double *x_out, double *temperatures_out, uint32_t point_count,
                                                     char *error_message, uint32_t error_message_length);

LVT_API uint8_t LVT_CALL lvt_compute_temperatures_2d(
    const double *region_x, const double *region_y, const double *region_width, const double *region_height,
    const double *region_lambda, const double *region_c_vol, const double *region_initial_temperature,
    const uint8_t *region_heat_source, uint32_t region_count, double spatial_step, double time_step,
    double simulation_time, double *x_out, double *y_out, double *temperatures_out, uint32_t sample_col_count,
    uint32_t sample_row_count, char *error_message, uint32_t error_message_length);

LVT_API uint8_t LVT_CALL lvt_compute_temperatures_2d_with_radiation(
    const double *region_x, const double *region_y, const double *region_width, const double *region_height,
    const double *region_lambda, const double *region_c_vol, const double *region_initial_temperature,
    const uint8_t *region_heat_source, uint32_t region_count, double spatial_step, double time_step,
    double simulation_time, double ambient_temperature, double emissivity, double plate_thickness, double *x_out,
    double *y_out, double *temperatures_out, uint32_t sample_col_count, uint32_t sample_row_count, char *error_message,
    uint32_t error_message_length);

#ifdef __cplusplus
}
#endif

#endif
