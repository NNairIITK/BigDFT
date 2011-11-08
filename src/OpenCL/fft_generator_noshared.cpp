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

static void decompose_radix(unsigned int radix_size, std::list<unsigned int> &sub_radixes){
  std::list<unsigned int>::iterator it;
  unsigned int available_rad[] = {2,3,5};
  std::list<unsigned int> available_radixes (available_rad, available_rad + sizeof(available_rad) / sizeof(unsigned int) );

  for( it = available_radixes.begin(); it != available_radixes.end(); it++ ){
    while(radix_size % *it == 0){
      radix_size /= *it;
      sub_radixes.push_back(*it);
    }
  }
  if(radix_size != 1){
    std::cerr<<"Invalid FFT size : "<<radix_size<<" is irreductible!"<<std::endl;
    exit(1);
  }

}

/*
static void generate_radix_no_shared(std::stringstream &program, unsigned int radix_size, unsigned int fft_size, std::string &sign, std::stringstream &div){
  unsigned int A=1,B;
  unsigned int order1[16]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
  unsigned int order2[16];
  unsigned int *order,*order_in,*order_out;
  unsigned int index;
  order_out = order2;
  order_in = order1;
  order = order1;
  program<<"  double2 tmp1,tmp2,val,w;\n\
  double2 tmp["<<radix_size<<"];\n";
  for(A=1,B=16/2;A<radix_size;A*=2,B/=2){
    for(int j=0;j<2;j++){
      for(int i=0; i<8;i++){
       order_out[j+i*2] = order_in[j*8+i];
      }
    }
    order = order_out; order_out = order_in; order_in = order;
    for(unsigned int i=0; i<radix_size;i+=2){
      if(A==1)
        program<<"  tmp1 = tmp2 = psi["<<order[i]<<"*ndat];\n\
  val = psi["<<order[i+1]<<"*ndat];\n";
      else
        program<<"  tmp1 = tmp2 = tmp["<<order[i]<<"];\n\
  val = tmp["<<order[i+1]<<"];\n";
      program<<"  tmp1 += val;\n";
      if(A*2==radix_size) program<<"  tmp1.x "<<div.str()<<";\n";
      if(A*2==radix_size) program<<"  tmp1.y "<<div.str()<<";\n";
      index = order[i+1]%(B)*(fft_size/(B*2));
      if(index == 0) {
        program<<"  tmp2 -= val;\n";
      } else if(index == fft_size/4) {
        program<<"  val = tmp2 - val;\n";
        program<<"  tmp2.x = "<<sign<<" val.y;\n";
        program<<"  tmp2.y = "<<sign<<" - val.x;\n";
      } else{
        program<<"  val = tmp2 - val;\n";
        program<<"  w.x = cosar["<<index<<"];\n\
  w.y = sinar["<<index<<"];\n\
  tmp2.x = val.x * w.x;\n\
  tmp2.x += "<<sign<<" val.y * w.y;\n\
  tmp2.y = "<<sign<<" - val.x * w.y;\n\
  tmp2.y += val.y * w.x;\n";
      }
      if(A*2==radix_size) program<<"  tmp2.x "<<div.str()<<";\n";
      if(A*2==radix_size) program<<"  tmp2.y "<<div.str()<<";\n";
      program<<"  tmp["<<order[i]<<"]=tmp1;\n\
  tmp["<<order[i+1]<<"]=tmp2;\n";
    }
  }

  for(int i=0; i<2; i++){
    for(int j=0; j<2; j++){
      for(int k=0; k<2; k++){
        for(int l=0; l<2; l++){
          program<<"  out[jg*"<<fft_size<<"+"<<((l*2+k)*2+j)*2+i<<"]=tmp["<<((i*2+j)*2+k)*2+l<<"];\n";
          //program<<"  out["<<((l*2+k)*2+j)*2+i<<"*ndat]=tmp["<<((i*2+j)*2+k)*2+l<<"];\n";
        }
      }
    }
  }
}*/

static void generate_radix_no_shared_generic(std::stringstream &program, unsigned int radix_size, unsigned int fft_size, unsigned int &A, unsigned int &B, std::string &in, std::string &out, bool transpose, unsigned int stride, std::string &sign, std::stringstream &div){
  std::list<unsigned int> sub_radixes;
  unsigned int *order1 = new unsigned int[radix_size];
  unsigned int *order2 = new unsigned int[radix_size];
  unsigned int *order,*order_in,*order_out;
  unsigned int index;
  std::list<unsigned int>::iterator it;

  order_out = order2; order_in = order1; order = order1;

  decompose_radix(radix_size, sub_radixes);

  for(unsigned int i=0; i<radix_size; i++)
    order1[i] = i;
  program<<"{  double2 cossin,val,t;\n\
  double2 tmp["<<radix_size<<"];\n\
  double2 tmp_val["<<sub_radixes.back()<<"];\n";
  for(unsigned int i=0; i < radix_size; i++)
    program<<"  tmp["<<i<<"] = "<<in<<"[("<<i<<"+ig*"<<radix_size<<")*ndat];\n";
  for(it = sub_radixes.begin(); it != sub_radixes.end(); it++, A *= *it ){
    B /= *it;
    for(unsigned int j=0; j<*it; j++){
      for(unsigned int i=0; i<radix_size/(*it); i++){
        order_out[j + i * *it] = order_in[j*radix_size/(*it) + i];
      }
    }
    order = order_out; order_out = order_in; order_in = order;
    for(unsigned int i=0; i < radix_size; i += *it){
      for(unsigned int j=0; j < *it; j++){
        program<<"  t = tmp["<<order[i]<<"];\n";
        for(unsigned int k=1; k < *it; k++){
          index = (k*j*fft_size/(*it))%fft_size;
          program<<"  cossin.x = cosar["<<index<<"];\n";
          program<<"  cossin.y = sinar["<<index<<"];\n";
          program<<"  val = tmp["<<order[i+k]<<"];\n";
          program<<"  t.x += val.x * cossin.x;\n";
          program<<"  t.x += "<<sign<<" val.y * cossin.y;\n";
          program<<"  t.y += "<<sign<<" - val.x * cossin.y;\n";
          program<<"  t.y += val.y * cossin.x;\n";
        }
        index = (order[i+j]%(B))*j*(fft_size/(B*(*it)));
        program<<"  cossin.x = cosar["<<index<<"];\n";
        program<<"  cossin.y = sinar["<<index<<"];\n";
        program<<"  tmp_val["<<j<<"].x = t.x * cossin.x;\n";
        program<<"  tmp_val["<<j<<"].x += "<<sign<<" t.y * cossin.y;\n";
        program<<"  tmp_val["<<j<<"].y = "<<sign<<" - t.x * cossin.y;\n";
        program<<"  tmp_val["<<j<<"].y += t.y * cossin.x;\n";
      }
      for(unsigned int j=0; j < *it; j++){
        if(B==1) program<<"  tmp_val["<<j<<"].x "<<div.str()<<";\n";
        if(B==1) program<<"  tmp_val["<<j<<"].y "<<div.str()<<";\n";
        program<<"  tmp["<<order[i+j]<<"] = tmp_val["<<j<<"];\n";
      }
    }
  }
  unsigned int *digits = new unsigned int[sub_radixes.size()];
  std::list<unsigned int>::reverse_iterator rit;
  for(unsigned int i=0; i < sub_radixes.size(); i++)
    digits[i] = 0;
  for(unsigned int i=0; i < radix_size; i++){
    unsigned int j;
    order[i] = 0;
    unsigned int radix_remainder = radix_size;
    for(rit = sub_radixes.rbegin(), j=0; rit != sub_radixes.rend(); rit++, j++){
      radix_remainder /= *rit;
      order[i] += digits[j]*radix_remainder;
    }
    for(rit = sub_radixes.rbegin(), j=0; rit != sub_radixes.rend(); rit++, j++){
      digits[j]++;
      if(digits[j] != *rit){
        break;
      }
      else
        digits[j]=0;
    }
  }
  for(unsigned int i=0;i <radix_size; i++){
    if( transpose )
      program<<"  "<<out<<"[jg*"<<fft_size<<"+ig+"<<order[i]*B<<"] = tmp["<<i<<"];\n";
    else
      program<<"  "<<out<<"[jg+(ig+"<<order[i]*B<<")*ndat] = tmp["<<i<<"];\n";
  }
  program<<"  }\n";
  delete[] order1;
  delete[] order2;
  delete[] digits;
}

static void generate_kernel_no_shared(std::stringstream &program, cl_uint fft_size, std::list<unsigned int> &radixes, bool reverse){
  if( reverse )
    program<<"__kernel void fftKernel_"<<fft_size<<"_r_d(uint n, uint ndat, __global const double2 *psi, __global double2 *out";
  else
    program<<"__kernel void fftKernel_"<<fft_size<<"_d(uint n, uint ndat, __global const double2 *psi, __global double2 *out";
  if(!use_constant_memory)
    program<<", __read_only image2d_t cosat";
  program<<"){\n\
  size_t jg = get_global_id(1);\n\
  size_t ig = get_global_id(0);\n\
  jg  = get_group_id(1) == get_num_groups(1) - 1 ? jg - ( get_global_size(1) - ndat ) : jg;\n\
  psi = &psi[jg];\n\
  //out = &out[jg];\n";
  std::string sign;
  std::stringstream div;
  std::string in="psi";
  std::string out="out";
  std::list<unsigned int>::iterator it;
  div<<std::showpoint<<std::scientific;
  div.precision(20);
  if( reverse ){
    sign = "-";
  } else {
    sign = "+";
  }
  div<<"";
  unsigned int A=1, B=fft_size;
  bool transpose = false;
  for( it=radixes.begin(); it!=radixes.end(); it++){
    if (it == --(radixes.end())){
      if( reverse)
        div<<"*="<<(double)1/(double)fft_size;
      transpose = true;
    }
    generate_radix_no_shared_generic(program, *it, fft_size, A, B, in, out, transpose, 0, sign, div);
  }
//  generate_radix_no_shared(program, 16, fft_size, sign, div);
  program<<"}\n";
}



extern "C" fft_code * generate_fft_program_no_shared(cl_uint fft_size, struct bigdft_device_infos * infos){
  unsigned int available_rad[] = {15,16};
  std::list<unsigned int> available_radixes (available_rad, available_rad + sizeof(available_rad) / sizeof(unsigned int) );
  std::list<unsigned int> radixes;
  std::list<unsigned int> uniq_radixes;
  std::list<unsigned int>::iterator it;
  std::stringstream program;
  fft_code* output = (fft_code*)malloc(sizeof(fft_code));

  generate_radixes(fft_size, available_radixes, radixes, uniq_radixes);

  generate_header(program);
  generate_cosin_tables(program, fft_size, output);

  generate_kernel_no_shared(program,fft_size,radixes,false);
  generate_kernel_no_shared(program,fft_size,radixes,true);
 
  output->code = (char *)malloc((program.str().size()+1)*sizeof(char));
  strcpy(output->code, program.str().c_str());
  return output;
}
