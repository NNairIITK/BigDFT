#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>
#include <malloc.h>
#include <math.h>
#include "fft_generator.h"
#define TWOPI 3.14159265358979323846264338327950288419716939937510*2

extern "C" char * generate_fft_program(cl_uint fft_size){
  unsigned int available_rad[] = {2,3,5,7,11};
  std::vector<unsigned int> available_radixes (available_rad, available_rad + sizeof(available_rad) / sizeof(unsigned int) );
  std::vector<unsigned int> radixes;
  std::vector<double> sines, cosines;
  std::stringstream program;
  char* output;
  unsigned int buffer_length;
  unsigned int buffer_width;
  unsigned int buffer_limit;
  cl_uint i;
  cl_uint fft_size_o=fft_size;
  for(i=0; i<fft_size; i++){
    sines.push_back(sin(TWOPI*i/(double)fft_size));
    cosines.push_back(cos(TWOPI*i/(double)fft_size));
  }

  program<<"\
#pragma OPENCL EXTENSION cl_khr_fp64: enable \n\
#define TWOPI 3.14159265358979323846264338327950288419716939937510*2\n\
";

  program<< std::showpoint<<std::scientific;
  program.precision(40);
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

  program<<"\
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
  w.x = cosar[r*(N/(2*A))];\
  w.y = sinar[r*(N/(2*A))];\
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
  w.x = cosar[r*(N/(3*A))];\
  w.y = sinar[r*(N/(3*A))];\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  val = in[(N/3)*2+A*b+a][il];\
  w.x = cosar[(r*(2*N/(3*A)))%N];\
  w.y = sinar[(r*(2*N/(3*A)))%N];\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  out[jl][il]=tmp;\
}\n\
#define radix4m(il,jl,N,A,B,in,out) \
{ \
  double2 tmp,val,w;\
  int a,b,p,r;\
  b = jl / (4*A);\
  r = jl % (4*A);\
  a = r % A;\
  val = in[A*b+a][il];\
  tmp.x = val.x;\
  tmp.y = val.y;\
  val = in[(N/4)+A*b+a][il];\
  w.x = cosar[r*(N/(4*A))];\
  w.y = sinar[r*(N/(4*A))];\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  val = in[(N/4)*2+A*b+a][il];\
  w.x = cosar[r*(2*N/(4*A))%N];\
  w.y = sinar[r*(2*N/(4*A))%N];\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  val = in[(N/4)*3+A*b+a][il];\
  w.x = cosar[r*(3*N/(4*A))%N];\
  w.y = sinar[r*(3*N/(4*A))%N];\
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
  w.x = cosar[r*(N/(5*A))];\
  w.y = sinar[r*(N/(5*A))];\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  val = in[(N/5)*2+A*b+a][il];\
  w.x = cosar[r*(2*N/(5*A))%N];\
  w.y = sinar[r*(2*N/(5*A))%N];\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  val = in[(N/5)*3+A*b+a][il];\
  w.x = cosar[r*(3*N/(5*A))%N];\
  w.y = sinar[r*(3*N/(5*A))%N];\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  val = in[(N/5)*4+A*b+a][il];\
  w.x = cosar[r*(4*N/(5*A))%N];\
  w.y = sinar[r*(4*N/(5*A))%N];\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  out[jl][il]=tmp;\
}\n\
#define radix7m(il,jl,N,A,B,in,out) \
{ \
  double2 tmp,val,w;\
  int a,b,p,r;\
  b = jl / (7*A);\
  r = jl % (7*A);\
  a = r % A;\
  val = in[A*b+a][il];\
  tmp.x = val.x;\
  tmp.y = val.y;\
  val = in[(N/7)+A*b+a][il];\
  w.x = cosar[r*(N/(7*A))];\
  w.y = sinar[r*(N/(7*A))];\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  val = in[(N/7)*2+A*b+a][il];\
  w.x = cosar[r*(2*N/(7*A))%N];\
  w.y = sinar[r*(2*N/(7*A))%N];\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  val = in[(N/7)*3+A*b+a][il];\
  w.x = cosar[r*(3*N/(7*A))%N];\
  w.y = sinar[r*(3*N/(7*A))%N];\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  val = in[(N/7)*4+A*b+a][il];\
  w.x = cosar[r*(4*N/(7*A))%N];\
  w.y = sinar[r*(4*N/(7*A))%N];\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  val = in[(N/7)*5+A*b+a][il];\
  w.x = cosar[r*(5*N/(7*A))%N];\
  w.y = sinar[r*(5*N/(7*A))%N];\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  val = in[(N/7)*6+A*b+a][il];\
  w.x = cosar[r*(6*N/(7*A))%N];\
  w.y = sinar[r*(6*N/(7*A))%N];\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  out[jl][il]=tmp;\
}\n\
#define radix11m(il,jl,N,A,B,in,out) \
{ \
  double2 tmp,val,w;\
  int a,b,p,r;\
  b = jl / (11*A);\
  r = jl % (11*A);\
  a = r % A;\
  val = in[A*b+a][il];\
  tmp.x = val.x;\
  tmp.y = val.y;\
  val = in[(N/11)+A*b+a][il];\
  w.x = cosar[r*(N/(11*A))];\
  w.y = sinar[r*(N/(11*A))];\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  val = in[(N/11)*2+A*b+a][il];\
  w.x = cosar[r*(2*N/(11*A))%N];\
  w.y = sinar[r*(2*N/(11*A))%N];\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  val = in[(N/11)*3+A*b+a][il];\
  w.x = cosar[r*(3*N/(11*A))%N];\
  w.y = sinar[r*(3*N/(11*A))%N];\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  val = in[(N/11)*4+A*b+a][il];\
  w.x = cosar[r*(4*N/(11*A))%N];\
  w.y = sinar[r*(4*N/(11*A))%N];\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  val = in[(N/11)*5+A*b+a][il];\
  w.x = cosar[r*(5*N/(11*A))%N];\
  w.y = sinar[r*(5*N/(11*A))%N];\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  val = in[(N/11)*6+A*b+a][il];\
  w.x = cosar[r*(6*N/(11*A))%N];\
  w.y = sinar[r*(6*N/(11*A))%N];\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  val = in[(N/11)*7+A*b+a][il];\
  w.x = cosar[r*(7*N/(11*A))%N];\
  w.y = sinar[r*(7*N/(11*A))%N];\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  val = in[(N/11)*8+A*b+a][il];\
  w.x = cosar[r*(8*N/(11*A))%N];\
  w.y = sinar[r*(8*N/(11*A))%N];\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  val = in[(N/11)*9+A*b+a][il];\
  w.x = cosar[r*(9*N/(11*A))%N];\
  w.y = sinar[r*(9*N/(11*A))%N];\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  val = in[(N/11)*10+A*b+a][il];\
  w.x = cosar[r*(10*N/(11*A))%N];\
  w.y = sinar[r*(10*N/(11*A))%N];\
  tmp.x += val.x * w.x;\
  tmp.x += val.y * w.y;\
  tmp.y -= val.x * w.y;\
  tmp.y += val.y * w.x;\
  out[jl][il]=tmp;\
}\n\
";
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
__kernel void fftKernel_"<<fft_size<<"_d(uint n, uint ndat, __global const double2 *psi, __global double2 *out){\n\
\n\
__local double2 tmp1[FFT_LENGTH][BUFFER_DEPTH];\n\
__local double2 tmp2[FFT_LENGTH][BUFFER_DEPTH];\n\
\n\
  size_t il = get_local_id(0);\n\
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

  std::vector<unsigned int>::iterator it;
//  std::cout<<"radixes : "<<std::endl;
  for( it = available_radixes.begin(); it < available_radixes.end(); it++ ){
    while(fft_size_o % *it == 0){
      fft_size_o /= *it;
      radixes.push_back(*it);
//      std::cout<<*it<<std::endl;
    }
  }
  if(fft_size_o != 1){
    std::cerr<<"Invalid FFT size : "<<fft_size_o<<" is irreducible!"<<std::endl;
    return NULL;
  }
  
  std::string in, out, tmp;
  in = "tmp1";
  out = "tmp2";

  for( it = radixes.begin(); it < radixes.end(); it++){
     B/=*it;
     program<<"  radix"<<*it<<"m(jlt, ilt, "<<fft_size<<", "<<A<<", "<<B<<", "<<in<<", "<<out<<");\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
";
     A*=*it;
     tmp=out; out=in; in=tmp;
  }

  program<<"\n\
  if(il<"<<fft_size<<")\n\
    out[jg*n+il] = "<<in<<"[il][jl];\n\
}\n\
"; 
 
  output = (char *)malloc((program.str().size()+1)*sizeof(char));
  strcpy(output, program.str().c_str());
  return output;
}
