#ifndef FFT_GENERATOR_H
#define FFT_GENERATOR_H
#include <CL/cl.h>
typedef struct {
  char *code;
  double *cossin;
} fft_code;
#ifdef __cplusplus 
  extern "C" fft_code * generate_fft_program(cl_uint fft_size);
  extern "C" fft_code * generate_fft_program_no_shared(cl_uint fft_size);
#else
  fft_code * generate_fft_program_no_shared(cl_uint fft_size);
  fft_code * generate_fft_program(cl_uint fft_size);
#endif
#endif
