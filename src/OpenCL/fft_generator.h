#ifndef FFT_GENERATOR_H
#define FFT_GENERATOR_H
#include <CL/cl.h>
#include "OpenCL_wrappers.h"
typedef struct {
  char *code;
  double *cossin;
} fft_code;
#ifdef __cplusplus 
  extern "C" fft_code * generate_fft_program(cl_uint fft_size, struct bigdft_device_infos * infos);
  extern "C" fft_code * generate_fft_program_no_shared(cl_uint fft_size, struct bigdft_device_infos * infos);
#else
  fft_code * generate_fft_program_no_shared(cl_uint fft_size, struct bigdft_device_infos * infos);
  fft_code * generate_fft_program(cl_uint fft_size, struct bigdft_device_infos * infos);
#endif
#endif
