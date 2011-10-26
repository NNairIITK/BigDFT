#ifndef FFT_GENERATOR_H
#define FFT_GENERATOR_H
#include <CL/cl.h>
#include "OpenCL_wrappers.h"
#define TWOPI 3.14159265358979323846264338327950288419716939937510*2
typedef struct {
  char *code;
  double *cossin;
} fft_code;
extern cl_uint use_constant_memory;
#ifdef __cplusplus 
  void sincos_access(std::stringstream &program, const char* index);
  void generate_header(std::stringstream &program);
  void generate_cosin_tables(std::stringstream &program, cl_uint fft_size, fft_code* output);
  void generate_radixes(cl_uint fft_size, std::list<unsigned int> &available_radixes, std::list<unsigned int> &radixes, std::list <unsigned int> &uniq_radixes);
  extern "C" fft_code * generate_fft_program(cl_uint fft_size, struct bigdft_device_infos * infos);
  extern "C" fft_code * generate_fft_program_no_shared(cl_uint fft_size, struct bigdft_device_infos * infos);
#else
  fft_code * generate_fft_program_no_shared(cl_uint fft_size, struct bigdft_device_infos * infos);
  fft_code * generate_fft_program(cl_uint fft_size, struct bigdft_device_infos * infos);
#endif
#endif
