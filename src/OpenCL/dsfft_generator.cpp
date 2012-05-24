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

/*
static void print_index(unsigned int jl, unsigned int N, unsigned int radix_size, unsigned int A) {
  unsigned int b,r,a;
  b = jl / ( radix_size * A );
  r = jl % ( radix_size * A );
  a = r  %  A;
  std::cout<<"{"<< A*b+a;
  std::cout<<","<< N/radix_size + A*b+a;
  for(unsigned int i=2; i<radix_size;i++){
     std::cout<<","<< N*i/radix_size + A*b+a;
  }
  std::cout<<"} ";
}

static void print_indexes( cl_uint fft_size, std::list<unsigned int> &radixes) {
  unsigned int A=1,B=fft_size;
  std::list<unsigned int>::iterator it;
  for( it = radixes.begin(); it != radixes.end(); it++){
     B/=*it;
     for(unsigned int i=0; i<fft_size; i++){
       print_index(i,fft_size,*it,A);
     }
     std::cout<<std::endl;
     A*=*it;
  }
  std::cout<<std::endl;
}*/

static void generate_buffer_size(std::stringstream &program, cl_uint fft_size, struct bigdft_device_infos * infos){
  size_t block_size_i, block_size_j, elem_per_thread;
  fft_compute_sizes(infos, fft_size, &block_size_i, &block_size_j, &elem_per_thread);
  unsigned int buffer_length = block_size_i * elem_per_thread;
  unsigned int buffer_width = block_size_j;

  program<<"#undef FFT_LENGTH\n\
#define FFT_LENGTH "<<buffer_length<<"\n\
#undef LINE_NUMBER\n\
#define LINE_NUMBER "<<buffer_width<<"\n\
#undef BUFFER_DEPTH\n\
#define BUFFER_DEPTH LINE_NUMBER+"<<(buffer_width>1?1:0)<<"\n";
}

static void generate_kernel(std::stringstream &program, cl_uint fft_size, std::list<unsigned int> &radixes, bool reverse, bool k, bool r2c, bool c2r, bool free, struct bigdft_device_infos * infos){
  size_t block_size_i, block_size_j, elem_per_thread;
  fft_compute_sizes(infos, fft_size, &block_size_i, &block_size_j, &elem_per_thread);
  program<<"__kernel void fftKernel_";
  if( k ) program<<"k_";
  program<<fft_size;
  if( reverse ) program<<"_r";
  if( r2c ) program<<"_r2c";
  if( c2r ) program<<"_c2r";
  if( free ) program<<"_f";
  program<<"_d(uint n, uint ndat,";
  if( r2c ) program<<" __global const double *psi,";
  else program<<" __global const double2 *psi,";
  if( c2r ) program<<" __global double *out";
  else program<<" __global double2 *out";
  if( k ) program<<", __global const double *k";
  if(!use_constant_memory) program<<", __read_only image2d_t cosat";
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
  jgt = jg - jl + jlt;\n";
  for(unsigned int i=0; i<elem_per_thread; i++){
    program<<"  tmp1[ilt+"<<i<<"*get_local_size(0)][jlt] = ilt + "<<i<<"*get_local_size(0) < "<<((free && !reverse) ? fft_size/2 : fft_size)<<" ? ";
    if( r2c ) program<<"(double2)( ";
    program<<"psi[jgt + ( ilt + "<<i<<"*get_local_size(0) ) * ndat]";
    if( r2c ) program<<", 0.0)";
    program<<" : (double2)(0.0, 0.0);\n";
  }
  program<<"  barrier(CLK_LOCAL_MEM_FENCE);\n";

  unsigned int A=1,B=fft_size;
  std::string in="tmp1", out="tmp2", tmp;
  std::list<unsigned int>::iterator it;
  
  for( it = radixes.begin(); it != radixes.end(); it++){
     B/=*it;
     for(unsigned int i=0; i<elem_per_thread; i++){
       program<<"  radix"<<(*it)<<"m(jlt, (ilt+"<<i<<"*get_local_size(0)), "<<fft_size<<", "<<A<<", "<<B<<", "<<in<<", "<<out;
       if( reverse ){
         if(it == --(radixes.end()))
           program<<",-,*="<<(double)1/(double)fft_size<<");\n";
         else
           program<<",-,);\n";
       } else
         program<<",+,);\n";
     }
     program<<"  barrier(CLK_LOCAL_MEM_FENCE);\n";
     A*=*it;
     tmp=out; out=in; in=tmp;
  }
  for(unsigned int i=0; i<elem_per_thread; i++){
    program<<"  if(get_local_size(0)*"<<i<<"+il < "<<((free && reverse) ? fft_size/2 : fft_size)<<")\n";
    program<<"    out[jg*"<<fft_size<<"+get_local_size(0)*"<<i<<"+il] = "<<in<<"[get_local_size(0)*"<<i<<"+il][jl]";
    if( c2r ) program<<".x";
    if( k ) program<<"*k[jg*"<<fft_size<<"+get_local_size(0)*"<<i<<"+il]";
    program<<";\n";
  }
  program<<"}\n";
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
//  print_indexes(fft_size,radixes);
  generate_kernel(program,fft_size,radixes,false,false,false,false,false,infos);
  generate_kernel(program,fft_size,radixes,false,false,true,false,false,infos);
  generate_kernel(program,fft_size,radixes,true,false,false,false,false,infos);
  generate_kernel(program,fft_size,radixes,true,false,false,true,false,infos);
  generate_kernel(program,fft_size,radixes,false,true,false,false,false,infos);
  generate_kernel(program,fft_size,radixes,false,false,false,false,true,infos);
  generate_kernel(program,fft_size,radixes,false,false,true,false,true,infos);
  generate_kernel(program,fft_size,radixes,true,false,false,false,true,infos);
  generate_kernel(program,fft_size,radixes,true,false,false,true,true,infos);
  generate_kernel(program,fft_size,radixes,false,true,false,false,true,infos);
 
  output->code = (char *)malloc((program.str().size()+1)*sizeof(char));
  strcpy(output->code, program.str().c_str());
  return output;
}

