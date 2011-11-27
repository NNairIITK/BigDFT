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

void sincos_access(std::stringstream &program, const char* index){
  if(use_constant_memory)
    program<<"  w.x = cosar["<<index<<"];\
  w.y = sinar["<<index<<"];";
  else
    program<<"  w = as_double2(read_imageui(cosat,smplr,(int2)("<<index<<",0)));";
}

void generate_header(std::stringstream &program) {
  program<<"#ifdef cl_khr_fp64\n\
#pragma OPENCL EXTENSION cl_khr_fp64: enable \n\
#elif defined (cl_amd_fp64)\n\
#pragma OPENCL EXTENSION cl_amd_fp64: enable \n\
#endif\n";
}

void generate_cosin_tables(std::stringstream &program, cl_uint fft_size, fft_code* output){
  std::vector<double> sines, cosines;
  for(unsigned int i=0; i<fft_size; i++){
    sines.push_back(sin(TWOPI*i/(double)fft_size));
    cosines.push_back(cos(TWOPI*i/(double)fft_size));
  }
  program<< std::showpoint<<std::scientific;
  program.precision(20);
  if(use_constant_memory){
    output->cossin=NULL;
    program<<"__constant double sinar["<<fft_size<<"] = { "<<sines[0];
    for(unsigned int i=1; i<fft_size; i++){
      program<<" ,\n"<<sines[i];
    }
    program<<"\n};\n";
    program<<"__constant double cosar["<<fft_size<<"] = { "<<cosines[0];
    for(unsigned int i=1; i<fft_size; i++){
      program<<" ,\n"<<cosines[i];
    }
    program<<"\n};\n";
  } else {
    output->cossin = (double *)malloc(fft_size*sizeof(double)*2);
    for(unsigned int i=0; i<fft_size; i++){
      output->cossin[2*i]=cosines[i];
      output->cossin[2*i+1]=sines[i];
    }
    program<<"const sampler_t smplr = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST;\n";
  }
}

void generate_radixes(cl_uint fft_size, std::list<unsigned int> &available_radixes, std::list<unsigned int> &radixes, std::list <unsigned int> &uniq_radixes){
  cl_uint fft_size_o=fft_size;
  std::list<unsigned int>::iterator it;

  for( it = available_radixes.begin(); it != available_radixes.end(); it++ ){
    while(fft_size_o % *it == 0){
      fft_size_o /= *it;
      radixes.push_back(*it);
    }
  }
  if(fft_size_o != 1){
    std::cerr<<"Invalid FFT size : "<<fft_size_o<<" is irreductible!"<<std::endl;
    exit(1);
  }
  uniq_radixes = radixes;
  uniq_radixes.unique();
}


