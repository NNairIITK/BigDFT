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

inline static void fft_compute_sizes(struct bigdft_device_infos * infos, cl_uint fft_size, size_t *block_size_i, size_t *block_size_j, size_t *elem_per_thread){
  assert(fft_size <= infos->LOCAL_MEM_SIZE/16);
  cl_uint shared_size_used;
  *elem_per_thread = 1;
  if(fft_size == infos->MAX_WORK_GROUP_SIZE){
    *block_size_i = fft_size;
    *block_size_j = 1;
  } else if (fft_size > infos->MAX_WORK_GROUP_SIZE){
    unsigned int ratio =(unsigned int)ceil((double)fft_size/(double)infos->MAX_WORK_GROUP_SIZE);
    size_t reduced_size = (size_t)ceil((double)fft_size/ratio);
    *block_size_i = (reduced_size / 16) * 16 + (reduced_size % 16 ? 16 : 0);
    *block_size_j = 1;
    *elem_per_thread = ratio;
  } else {
    if(fft_size <= 64)
       shared_size_used=512;
    else
       shared_size_used=1024;
    int FFT_LENGTH = (fft_size / 16) * 16 + (fft_size % 16 ? 16 : 0);
    int LINE_LIMIT = ( (shared_size_used/FFT_LENGTH) & (~3) ) <= infos->MAX_WORK_GROUP_SIZE/FFT_LENGTH ? (shared_size_used/FFT_LENGTH) & (~3) : infos->MAX_WORK_GROUP_SIZE/FFT_LENGTH ;
    int LINE_NUMBER = 1;
    while( LINE_LIMIT /= 2 )
    LINE_NUMBER *= 2;
    if( LINE_NUMBER > FFT_LENGTH )
       LINE_NUMBER = FFT_LENGTH;
    if( LINE_NUMBER < 1 )
       LINE_NUMBER=1;
    *block_size_i=FFT_LENGTH, *block_size_j=LINE_NUMBER;
  }
}

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
