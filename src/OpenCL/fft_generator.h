#ifndef FFT_GENERATOR_H
#define FFT_GENERATOR_H
#include <CL/cl.h>
#ifdef __cplusplus 
  extern "C" char * generate_fft_program(cl_uint fft_size);
#else
  char * generate_fft_program(cl_uint fft_size);
#endif
#endif
