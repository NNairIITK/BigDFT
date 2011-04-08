#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>
#include <malloc.h>
#include "fft_generator.h"

extern "C" char * generate_fft_program(cl_uint fft_size){
  unsigned int available_rad[] = {2,3,5};
  std::vector<unsigned int> available_radixes (available_rad, available_rad + sizeof(available_rad) / sizeof(unsigned int) );
  std::vector<unsigned int> radixes;
  std::stringstream program;
  char* output;
  unsigned int buffer_length;
  unsigned int buffer_width;
  program<<"\
#pragma OPENCL EXTENSION cl_khr_fp64: enable \n\
#define TWOPI 3.14159265358979323846264338327950288419716939937510*2\n\
#define radix2m(il,jl,N,A,B,in,out) \
{ \
  double2 tmp,val,w;\
  int a,b,p,r;\
  b = jl / (2*A);\
  r = jl % (2*A);\
/*  p = r / A;*/\
  a = r % A;\
  val = in[A*b+a][il];\
  tmp.x = val.x;\
  tmp.y = val.y;\
  val = in[(N/2)+A*b+a][il];\
  w.x = cos((TWOPI/(A*2))*r);\
  w.y = sin((TWOPI/(A*2))*r);\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  out[jl][il]=tmp;\
}\n\
#define radix3m(il,jl,N,A,B,in,out) \
{ \
  double2 tmp,val,w;\
  int a,b,p,r;\
  b = jl / (3*A);\
  r = jl % (3*A);\
  a = r % A;\
  val = in[A*b+a][il];\
  tmp.x = val.x;\
  tmp.y = val.y;\
  val = in[(N/3)+A*b+a][il];\
  w.x = cos((TWOPI/(A*3))*r);\
  w.y = sin((TWOPI/(A*3))*r);\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  val = in[(N/3)*2+A*b+a][il];\
  w.x = cos((TWOPI/(A*3))*r*2);\
  w.y = sin((TWOPI/(A*3))*r*2);\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  out[jl][il]=tmp;\
}\n\
#define radix5m(il,jl,N,A,B,in,out) \
{ \
  double2 tmp,val,w;\
  int a,b,p,r;\
  b = jl / (5*A);\
  r = jl % (5*A);\
  a = r % A;\
  val = in[A*b+a][il];\
  tmp.x = val.x;\
  tmp.y = val.y;\
  val = in[(N/5)+A*b+a][il];\
  w.x = cos((TWOPI/(A*5))*r);\
  w.y = sin((TWOPI/(A*5))*r);\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  val = in[(N/5)*2+A*b+a][il];\
  w.x = cos((TWOPI/(A*5))*r*2);\
  w.y = sin((TWOPI/(A*5))*r*2);\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  val = in[(N/5)*3+A*b+a][il];\
  w.x = cos((TWOPI/(A*5))*r*3);\
  w.y = sin((TWOPI/(A*5))*r*3);\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  val = in[(N/5)*4+A*b+a][il];\
  w.x = cos((TWOPI/(A*5))*r*4);\
  w.y = sin((TWOPI/(A*5))*r*4);\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  out[jl][il]=tmp;\
}\n\
";

  buffer_length = (fft_size / 16) * 16 + (fft_size % 16 ? 16 : 0);
//  std::cout<<"buffer length : "<<buffer_length<<std::endl;
  buffer_width = (128/buffer_length) & (~1);
//  std::cout<<"buffer width : "<<buffer_width<<std::endl;

  program<<"#undef FFT_LENGTH\n\
#define FFT_LENGTH "<<buffer_length<<"\n\
#undef LINE_NUMBER\n\
#define LINE_NUMBER "<<buffer_width<<"\n\
#undef BUFFER_DEPTH\n\
#define BUFFER_DEPTH LINE_NUMBER+"<<(buffer_width>1?1:0)<<"\n\
__kernel void fftKernel_"<<fft_size<<"_d(uint n, uint ndat, __global const double2 *psi, __global double2 *out){\n\
\n\
__local double2 tmp1[FFT_LENGTH][BUFFER_DEPTH];\n\
__local double2 tmp2[FFT_LENGTH][BUFFER_DEPTH];\n\
\n\
  const uint fftl="<<fft_size<<";\n\
  size_t il = get_local_id(0);\n\
  size_t jl = get_local_id(1);\n\
  size_t jg = get_global_id(1);\n\
  size_t ilt = jl+(il/LINE_NUMBER)*LINE_NUMBER;\n\
  size_t jlt = il%LINE_NUMBER;\n\
  ptrdiff_t jgt = get_group_id(1);\n\
  jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1) - ndat ) : jg;\n\
  jgt = jg - jl + jlt;\n\
  tmp1[ilt][jlt] = ilt < fftl ? psi[jgt + ( ilt ) * ndat] : 0.0;\n\
  \n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
";

  unsigned int A,B;

  A=1;
  B=fft_size;

  std::vector<unsigned int>::iterator it;
//  std::cout<<"radixes : "<<std::endl;
  for( it = available_radixes.begin(); it < available_radixes.end(); it++ ){
    while(fft_size % *it == 0){
      fft_size /= *it;
      radixes.push_back(*it);
//      std::cout<<*it<<std::endl;
    }
  }
  if(fft_size != 1){
    std::cerr<<"Invalid FFT size : "<<fft_size<<" is irreducible!"<<std::endl;
    return NULL;
  }
  
  std::string in, out, tmp;
  in = "tmp1";
  out = "tmp2";

  for( it = radixes.begin(); it < radixes.end(); it++){
     B/=*it;
     program<<"  radix"<<*it<<"m(jlt, ilt, fftl, "<<A<<", "<<B<<", "<<in<<", "<<out<<");\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
";
     A*=*it;
     tmp=out; out=in; in=tmp;
  }

  program<<"\n\
  if(il<fftl)\n\
    out[jg*n+il] = "<<in<<"[il][jl];\n\
}\n\
"; 
 
  output = (char *)malloc((program.str().size()+1)*sizeof(char));
  strcpy(output, program.str().c_str());
  return output;
}
