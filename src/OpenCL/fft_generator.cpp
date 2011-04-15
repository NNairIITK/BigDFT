#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>
#include <list>
#include <malloc.h>
#include <math.h>
#include "fft_generator.h"
#define TWOPI 3.14159265358979323846264338327950288419716939937510*2

extern cl_uint use_constant_memory;

void sincos_access(std::stringstream &program, const char* index){
  if(use_constant_memory)
    program<<"  w.x = cosar["<<index<<"];\
  w.y = sinar["<<index<<"];";
  else
    program<<"  w = as_double2(read_imageui(cosat,smplr,(int2)("<<index<<",0)));";
}

void generate_radix(std::stringstream &program, unsigned int radix_size){
  program<<"#define radix"<<radix_size<<"m(il,jl,N,A,B,in,out,sign,div)\
  { \
  double2 tmp,val,w;\
  int a,b,p,r;\
  b = jl / ("<<radix_size<<"*A);\
  r = jl % ("<<radix_size<<"*A);\
  a = r % A;\
  val = in[A*b+a][il];\
  tmp.x = val.x;\
  tmp.y = val.y;\
  val = in[(N/"<<radix_size<<")+A*b+a][il];";
  if(use_constant_memory)
    program<<"  w.x = cosar[r*(N/("<<radix_size<<"*A))];\
  w.y = sinar[r*(N/("<<radix_size<<"*A))];";
  else
    program<<"  w = as_double2(read_imageui(cosat,smplr,(int2)(r*(N/("<<radix_size<<"*A)),0)));";
  
  for(int i=2; i<radix_size;i++){
    program<<"  tmp.x += val.x * w.x;\
  tmp.x += sign val.y * w.y;\
  tmp.y += sign - val.x * w.y;\
  tmp.y += val.y * w.x;\
  val = in[(N/"<<radix_size<<")*"<<i<<"+A*b+a][il];";
    if(use_constant_memory)
      program<<"  w.x = cosar[r*("<<i<<"*N/("<<radix_size<<"*A))%N];\
  w.y = sinar[r*("<<i<<"*N/("<<radix_size<<"*A))%N];";
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

extern "C" fft_code * generate_fft_program(cl_uint fft_size){
  unsigned int available_rad[] = {2,3,5,7,11,13,17,19,23,29,31};
  std::list<unsigned int> available_radixes (available_rad, available_rad + sizeof(available_rad) / sizeof(unsigned int) );
  std::list<unsigned int> radixes;
  std::vector<double> sines, cosines;
  std::stringstream program;
  fft_code* output = (fft_code*)malloc(sizeof(fft_code));
  unsigned int buffer_length;
  unsigned int buffer_width;
  unsigned int buffer_limit;
  cl_uint i;
  cl_uint fft_size_o=fft_size;
  for(i=0; i<fft_size; i++){
    sines.push_back(sin(TWOPI*i/(double)fft_size));
    cosines.push_back(cos(TWOPI*i/(double)fft_size));
  }
  std::list<unsigned int>::iterator it;
  for( it = available_radixes.begin(); it != available_radixes.end(); it++ ){
    while(fft_size_o % *it == 0){
      fft_size_o /= *it;
      radixes.push_back(*it);
    }
  }
  if(fft_size_o != 1){
    std::cerr<<"Invalid FFT size : "<<fft_size_o<<" is irreductible!"<<std::endl;
    return NULL;
  }
  std::list <unsigned int> uniq_radixes = radixes;
  uniq_radixes.unique();

  program<<"\
#pragma OPENCL EXTENSION cl_khr_fp64: enable \n\
#define TWOPI 3.14159265358979323846264338327950288419716939937510*2\n\
";
  program<< std::showpoint<<std::scientific;
  program.precision(40);
  if(use_constant_memory){
    output->cossin=NULL;
    program<<"__constant double sinar["<<fft_size<<"] = { "<<sines[0];
    for(i=1; i<fft_size; i++){
      program<<" ,\n"<<sines[i];
    }
    program<<"\n};\n";

    program<<"__constant double cosar["<<fft_size<<"] = { "<<cosines[0];
    for(i=1; i<fft_size; i++){
      program<<" ,\n"<<cosines[i];
    }
    program<<"\n};\n";
  } else {
    output->cossin = (double *)malloc(fft_size*sizeof(double)*2);
    for(i=0; i<fft_size; i++){
      output->cossin[2*i]=cosines[i];
      output->cossin[2*i+1]=sines[i];
    }
    program<<"const sampler_t smplr = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST;\n";
  }

  for( it = uniq_radixes.begin(); it != uniq_radixes.end(); it++ )
    generate_radix(program,*it);

  cl_uint shared_size_used=0;
  if(fft_size <= 64)
    shared_size_used=512;
  else if(fft_size <= 256)
    shared_size_used=1024;
  buffer_length = (fft_size / 16) * 16 + (fft_size % 16 ? 16 : 0);
  buffer_limit = (shared_size_used/buffer_length);
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
#define BUFFER_DEPTH LINE_NUMBER+"<<(buffer_width>1?1:0)<<"\n\
__kernel void fftKernel_"<<fft_size<<"_d(uint n, uint ndat, __global const double2 *psi, __global double2 *out";
  if(!use_constant_memory)
    program<<", __read_only image2d_t cosat";
  program<<"){\n\
\n\
__local double2 tmp1[FFT_LENGTH][BUFFER_DEPTH];\n\
__local double2 tmp2[FFT_LENGTH][BUFFER_DEPTH];\n\
\n";
  program<<"  size_t il = get_local_id(0);\n\
  size_t jl = get_local_id(1);\n\
  size_t jg = get_global_id(1);\n\
  size_t ilt = jl+(il/LINE_NUMBER)*LINE_NUMBER;\n\
  size_t jlt = il%LINE_NUMBER;\n\
  ptrdiff_t jgt = get_group_id(1);\n\
  jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1) - ndat ) : jg;\n\
  jgt = jg - jl + jlt;\n\
  tmp1[ilt][jlt] = ilt < "<<fft_size<<" ? psi[jgt + ( ilt ) * ndat] : 0.0;\n\
  \n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
";

  unsigned int A,B;

  A=1;
  B=fft_size;


  
  std::string in, out, tmp;
  in = "tmp1";
  out = "tmp2";

  for( it = radixes.begin(); it != radixes.end(); it++){
     B/=*it;
     program<<"  radix"<<*it<<"m(jlt, ilt, "<<fft_size<<", "<<A<<", "<<B<<", "<<in<<", "<<out<<",+,);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
";
     A*=*it;
     tmp=out; out=in; in=tmp;
  }

  program<<"\n\
  if(il<"<<fft_size<<")\n\
    out[jg*"<<fft_size<<"+il] = "<<in<<"[il][jl];\n\
}\n\
"; 
  program<<"__kernel void fftKernel_"<<fft_size<<"_r_d(uint n, uint ndat, __global const double2 *psi, __global double2 *out";
  if(!use_constant_memory)
    program<<", __read_only image2d_t cosat";
  program<<"){\n\
\n\
__local double2 tmp1[FFT_LENGTH][BUFFER_DEPTH];\n\
__local double2 tmp2[FFT_LENGTH][BUFFER_DEPTH];\n\
\n";
  program<<"  size_t il = get_local_id(0);\n\
  size_t jl = get_local_id(1);\n\
  size_t jg = get_global_id(1);\n\
  size_t ilt = jl+(il/LINE_NUMBER)*LINE_NUMBER;\n\
  size_t jlt = il%LINE_NUMBER;\n\
  ptrdiff_t jgt = get_group_id(1);\n\
  jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1) - ndat ) : jg;\n\
  jgt = jg - jl + jlt;\n\
  tmp1[ilt][jlt] = ilt < "<<fft_size<<" ? psi[jgt + ( ilt ) * ndat] : 0.0;\n\
  \n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
";

  A=1;
  B=fft_size;

  in = "tmp1";
  out = "tmp2";

  for( it = radixes.begin(); it != radixes.end(); it++){
     B/=*it;
     if(it == --(radixes.end()))
       program<<"  radix"<<*it<<"m(jlt, ilt, "<<fft_size<<", "<<A<<", "<<B<<", "<<in<<", "<<out<<",-,*="<<(double)1/(double)fft_size<<");\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
";
     else
       program<<"  radix"<<*it<<"m(jlt, ilt, "<<fft_size<<", "<<A<<", "<<B<<", "<<in<<", "<<out<<",-,);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
";
     A*=*it;
     tmp=out; out=in; in=tmp;
  }

  program<<"\n\
  if(il<"<<fft_size<<")\n\
    out[jg*"<<fft_size<<"+il] = "<<in<<"[il][jl];\n\
}\n\
";

  output->code = (char *)malloc((program.str().size()+1)*sizeof(char));
  strcpy(output->code, program.str().c_str());
  return output;
}
