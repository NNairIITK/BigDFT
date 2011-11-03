#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>
#include <list>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include "fft_generator.h"


static void generate_radix_macro(std::stringstream &program, unsigned int radix_size){
  program<<"#define radix"<<radix_size<<"m(il,jl,N,A,B,in,out,sign,div)\
  { \
  double2 tmp,val,w;\
  int a,b,p,r,id,idd;\
  b = jl / ("<<radix_size<<"*A);\
  r = jl % ("<<radix_size<<"*A);\
  a = r % A;\
  tmp = in[A*b+a][il];\
  val = in[(N/"<<radix_size<<")+A*b+a][il];\
  idd = id = r*B;";
  if(use_constant_memory)
    program<<"  w.x = cosar[idd];\
  w.y = sinar[idd];";
  else
    program<<"  w = as_double2(read_imageui(cosat,smplr,(int2)(r*(N/("<<radix_size<<"*A)),0)));";
  
  for(unsigned int i=2; i<radix_size;i++){
    program<<"  tmp.x += val.x * w.x;\
  tmp.x += sign val.y * w.y;\
  tmp.y += sign - val.x * w.y;\
  tmp.y += val.y * w.x;\
  val = in[(N/"<<radix_size<<")*"<<i<<"+A*b+a][il];\
  idd += id;\
  idd>=N ? idd-=N : idd;";
    if(use_constant_memory)
      program<<"  w.x = cosar[idd];\
  w.y = sinar[idd];";
    else
      program<<"  w = as_double2(read_imageui(cosat,smplr,(int2)(r*("<<i<<"*N/("<<radix_size<<"*A))%N,0)));";
  }

  program<<"  tmp.x += val.x * w.x;\
  tmp.x += sign val.y * w.y;\
  tmp.x div;\
  tmp.y += sign - val.x * w.y;\
  tmp.y += val.y * w.x;\
  tmp.y div;\
  out[jl][il]=tmp;\
}\n";
}

static void generate_buffer_size(std::stringstream &program, cl_uint fft_size, struct bigdft_device_infos * infos){
  unsigned int buffer_length;
  unsigned int buffer_width;
  unsigned int buffer_limit;
  cl_uint shared_size_used=0;

  if(fft_size <= 64)
    shared_size_used=512;
  else if(fft_size <= 256)
    shared_size_used=1024;
  buffer_length = (fft_size / 16) * 16 + (fft_size % 16 ? 16 : 0);
  buffer_limit = (shared_size_used/buffer_length);
  buffer_limit = ( (shared_size_used/buffer_length) & (~3) ) <= infos->MAX_WORK_GROUP_SIZE/buffer_length ? (shared_size_used/buffer_length) & (~3) : infos->MAX_WORK_GROUP_SIZE/buffer_length ;
  buffer_width = 1;
  while( buffer_limit /= 2 )
    buffer_width *= 2;
  if( buffer_width > buffer_length)
    buffer_width = buffer_length;
  if( buffer_width < 1 )
    buffer_width = 1;

  program<<"#undef FFT_LENGTH\n\
#define FFT_LENGTH "<<buffer_length<<"\n\
#undef LINE_NUMBER\n\
#define LINE_NUMBER "<<buffer_width<<"\n\
#undef BUFFER_DEPTH\n\
#define BUFFER_DEPTH LINE_NUMBER+"<<(buffer_width>1?1:0)<<"\n";
}

static void generate_kernel(std::stringstream &program, cl_uint fft_size, std::list<unsigned int> &radixes, bool reverse){
  if( reverse )
    program<<"__kernel void fftKernel_"<<fft_size<<"_r_d(uint n, uint ndat, __global const double2 *psi, __global double2 *out";
  else
    program<<"__kernel void fftKernel_"<<fft_size<<"_d(uint n, uint ndat, __global const double2 *psi, __global double2 *out";
  if(!use_constant_memory)
    program<<", __read_only image2d_t cosat";
  program<<"){\n\
__local double2 tmp1[FFT_LENGTH][BUFFER_DEPTH];\n\
__local double2 tmp2[FFT_LENGTH][BUFFER_DEPTH];\n\
  size_t il = get_local_id(0);\n\
  size_t jl = get_local_id(1);\n\
  size_t jg = get_global_id(1);\n\
  size_t ilt = jl+(il/LINE_NUMBER)*LINE_NUMBER;\n\
  size_t jlt = il%LINE_NUMBER;\n\
  ptrdiff_t jgt = get_group_id(1);\n\
  jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1) - ndat ) : jg;\n\
  jgt = jg - jl + jlt;\n\
  tmp1[ilt][jlt] = ilt < "<<fft_size<<" ? psi[jgt + ( ilt ) * ndat] : 0.0;\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n";

  unsigned int A=1,B=fft_size;
  std::string in="tmp1", out="tmp2", tmp;
  std::list<unsigned int>::iterator it;

  for( it = radixes.begin(); it != radixes.end(); it++){
     B/=*it;
     program<<"  radix"<<*it<<"m(jlt, ilt, "<<fft_size<<", "<<A<<", "<<B<<", "<<in<<", "<<out;
     if( reverse ){
       if(it == --(radixes.end()))
         program<<",-,*="<<(double)1/(double)fft_size<<");\n";
       else
         program<<",-,);\n";
     } else
       program<<",+,);\n";
     program<<"  barrier(CLK_LOCAL_MEM_FENCE);\n";
     A*=*it;
     tmp=out; out=in; in=tmp;
  }
  program<<"  if(il<"<<fft_size<<")\n\
    out[jg*"<<fft_size<<"+il] = "<<in<<"[il][jl];\n\
}\n";
}

static void generate_kernel_k(std::stringstream &program, cl_uint fft_size, std::list<unsigned int> &radixes, bool reverse){
  if( reverse )
    program<<"__kernel void fftKernel_k_"<<fft_size<<"_r_d(uint n, uint ndat, __global const double2 *psi, __global double2 *out, __global const double *k";
  else
    program<<"__kernel void fftKernel_k_"<<fft_size<<"_d(uint n, uint ndat, __global const double2 *psi, __global double2 *out, __global const double *k";
  if(!use_constant_memory)
    program<<", __read_only image2d_t cosat";
  program<<"){\n\
__local double2 tmp1[FFT_LENGTH][BUFFER_DEPTH];\n\
__local double2 tmp2[FFT_LENGTH][BUFFER_DEPTH];\n\
  size_t il = get_local_id(0);\n\
  size_t jl = get_local_id(1);\n\
  size_t jg = get_global_id(1);\n\
  size_t ilt = jl+(il/LINE_NUMBER)*LINE_NUMBER;\n\
  size_t jlt = il%LINE_NUMBER;\n\
  ptrdiff_t jgt = get_group_id(1);\n\
  jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1) - ndat ) : jg;\n\
  jgt = jg - jl + jlt;\n\
  tmp1[ilt][jlt] = ilt < "<<fft_size<<" ? psi[jgt + ( ilt ) * ndat] : 0.0;\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n";

  unsigned int A=1,B=fft_size;
  std::string in="tmp1", out="tmp2", tmp;
  std::list<unsigned int>::iterator it;

  for( it = radixes.begin(); it != radixes.end(); it++){
     B/=*it;
     program<<"  radix"<<*it<<"m(jlt, ilt, "<<fft_size<<", "<<A<<", "<<B<<", "<<in<<", "<<out;
     if( reverse ){
       if(it == --(radixes.end()))
         program<<",-,*="<<(double)1/(double)fft_size<<");\n";
       else
         program<<",-,);\n";
     } else
       program<<",+,);\n";
     program<<"  barrier(CLK_LOCAL_MEM_FENCE);\n";
     A*=*it;
     tmp=out; out=in; in=tmp;
  }
  program<<"  if(il<"<<fft_size<<")\n\
    out[jg*"<<fft_size<<"+il] = "<<in<<"[il][jl]*k[jg*"<<fft_size<<"+il];\n\
}\n";
}


extern "C" fft_code * generate_fft_program(cl_uint fft_size, struct bigdft_device_infos * infos){
  unsigned int available_rad[] = {2,3,5,7,11,13,17,19,23,29,31};
  std::list<unsigned int> available_radixes (available_rad, available_rad + sizeof(available_rad) / sizeof(unsigned int) );
  std::list<unsigned int> radixes;
  std::list<unsigned int> uniq_radixes;
  std::list<unsigned int>::iterator it;
  std::stringstream program;
  fft_code* output = (fft_code*)malloc(sizeof(fft_code));

  generate_radixes(fft_size, available_radixes, radixes, uniq_radixes);

  generate_header(program);
  generate_cosin_tables(program, fft_size, output);

  for( it = uniq_radixes.begin(); it != uniq_radixes.end(); it++ )
    generate_radix_macro(program,*it);

  generate_buffer_size(program,fft_size,infos);

  generate_kernel(program,fft_size,radixes,false);
  generate_kernel(program,fft_size,radixes,true);
  generate_kernel_k(program,fft_size,radixes,false);
 
  output->code = (char *)malloc((program.str().size()+1)*sizeof(char));
  strcpy(output->code, program.str().c_str());
  return output;
}

