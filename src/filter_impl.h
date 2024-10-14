#pragma once

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

void filter_impl(uint8_t* buffer, int width, int height, int plane_stride, int pixel_stride,
                    const char *bg, int opening_size, int th_low, int th_high,
                    int bg_sampling_rate, int bg_number_frame, struct timespec time);

#ifdef __cplusplus
}
#endif
