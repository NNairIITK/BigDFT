//! @file
//!  Magic filter (??)
//!
//! @author
//!    Copyright (C) 2009-2011 BigDFT group 
//!    This file is distributed under the terms of the
//!    GNU General Public License, see ~/COPYING file
//!    or http://www.gnu.org/copyleft/gpl.txt .
//!    For the list of contributors, see ~/AUTHORS 


#include <stdlib.h>
#include <stdio.h>
#include <thread_engine.h>
#include <time.h>
#define CACHE_SIZE 6291456
#define NB_SYB 2
#define NB_CORE 8

///void rdtsc_(long long unsigned int * t) {
//  rdtscll(*t);
//}

void nanosec_(long long unsigned int * t){
  struct timespec time;
  clock_gettime(CLOCK_REALTIME, &time);
  *t = time.tv_sec;
  *t *= 1000000000;
  *t += time.tv_nsec;
}

const double filt[] = { 8.4334247333529341094733325815816e-7,
                       -0.1290557201342060969516786758559028e-4,
                        0.8762984476210559564689161894116397e-4,
                       -0.30158038132690463167163703826169879e-3,
                        0.174723713672993903449447812749852942e-2,
                       -0.942047030201080385922711540948195075e-2,
                        0.2373821463724942397566389712597274535e-1,
                        0.612625895831207982195380597e-1,
                        0.9940415697834003993178616713,
                       -0.604895289196983516002834636e-1,
                       -0.2103025160930381434955489412839065067e-1,
                        0.1337263414854794752733423467013220997e-1,
                       -0.344128144493493857280881509686821861e-2,
                        0.49443227688689919192282259476750972e-3,
                       -0.5185986881173432922848639136911487e-4,
                        2.72734492911979659657715313017228e-6};

void magicfilter_partial(unsigned int n1, unsigned int start, unsigned int ndat, unsigned int ld, const double *x, double *y) {
  int i,j;
  double result00,result10,result20,result30;//,result4,result5,result6,result7,result8,result9,result10,result11,result12,result13,result14,result15;
  double result01,result11,result21,result31;//,result4,result5,result6,result7,result8,result9,result10,result11,result12,result13,result14,result15;
  double result0;
  double result1;
  const double *x_t0, *x_t1;
  for(j=start; j<start+ndat; j+=2 ) {
    x_t0 = x+j*n1;
    x_t1 = x+(j+1)*n1;
    for(i=0; i<n1; i+=4) {
      result00 = 0;
      result10 = 0;
      result20 = 0;
      result30 = 0;

      result01 = 0;
      result11 = 0;
      result21 = 0;
      result31 = 0;

      result00 += x_t0[(i+n1-8)%n1]*filt[0];
      result10 += x_t0[(i+n1-7)%n1]*filt[0];
      result20 += x_t0[(i+n1-6)%n1]*filt[0];
      result30 += x_t0[(i+n1-5)%n1]*filt[0];

      result01 += x_t1[(i+n1-8)%n1]*filt[0];
      result11 += x_t1[(i+n1-7)%n1]*filt[0];
      result21 += x_t1[(i+n1-6)%n1]*filt[0];
      result31 += x_t1[(i+n1-5)%n1]*filt[0];

      result00 += x_t0[(i+n1-7)%n1]*filt[1];
      result10 += x_t0[(i+n1-6)%n1]*filt[1];
      result20 += x_t0[(i+n1-5)%n1]*filt[1];
      result30 += x_t0[(i+n1-4)%n1]*filt[1];

      result01 += x_t1[(i+n1-7)%n1]*filt[1];
      result11 += x_t1[(i+n1-6)%n1]*filt[1];
      result21 += x_t1[(i+n1-5)%n1]*filt[1];
      result31 += x_t1[(i+n1-4)%n1]*filt[1];

      result00 += x_t0[(i+n1-6)%n1]*filt[2];
      result10 += x_t0[(i+n1-5)%n1]*filt[2];
      result20 += x_t0[(i+n1-4)%n1]*filt[2];
      result30 += x_t0[(i+n1-3)%n1]*filt[2];

      result01 += x_t1[(i+n1-6)%n1]*filt[2];
      result11 += x_t1[(i+n1-5)%n1]*filt[2];
      result21 += x_t1[(i+n1-4)%n1]*filt[2];
      result31 += x_t1[(i+n1-3)%n1]*filt[2];

      result00 += x_t0[(i+n1-5)%n1]*filt[3];
      result10 += x_t0[(i+n1-4)%n1]*filt[3];
      result20 += x_t0[(i+n1-3)%n1]*filt[3];
      result30 += x_t0[(i+n1-2)%n1]*filt[3];

      result01 += x_t1[(i+n1-5)%n1]*filt[3];
      result11 += x_t1[(i+n1-4)%n1]*filt[3];
      result21 += x_t1[(i+n1-3)%n1]*filt[3];
      result31 += x_t1[(i+n1-2)%n1]*filt[3];

      result00 += x_t0[(i+n1-4)%n1]*filt[4];
      result10 += x_t0[(i+n1-3)%n1]*filt[4];
      result20 += x_t0[(i+n1-2)%n1]*filt[4];
      result30 += x_t0[(i+n1-1)%n1]*filt[4];

      result01 += x_t1[(i+n1-4)%n1]*filt[4];
      result11 += x_t1[(i+n1-3)%n1]*filt[4];
      result21 += x_t1[(i+n1-2)%n1]*filt[4];
      result31 += x_t1[(i+n1-1)%n1]*filt[4];

      result00 += x_t0[(i+n1-3)%n1]*filt[5];
      result10 += x_t0[(i+n1-2)%n1]*filt[5];
      result20 += x_t0[(i+n1-1)%n1]*filt[5];
      result30 += x_t0[ i         ]*filt[5];

      result01 += x_t1[(i+n1-3)%n1]*filt[5];
      result11 += x_t1[(i+n1-2)%n1]*filt[5];
      result21 += x_t1[(i+n1-1)%n1]*filt[5];
      result31 += x_t1[ i         ]*filt[5];

      result00 += x_t0[(i+n1-2)%n1]*filt[6];
      result10 += x_t0[(i+n1-1)%n1]*filt[6];
      result20 += x_t0[ i         ]*filt[6];
      result30 += x_t0[(i   +1)%n1]*filt[6];

      result01 += x_t1[(i+n1-2)%n1]*filt[6];
      result11 += x_t1[(i+n1-1)%n1]*filt[6];
      result21 += x_t1[ i         ]*filt[6];
      result31 += x_t1[(i   +1)%n1]*filt[6];

      result00 += x_t0[(i+n1-1)%n1]*filt[7];
      result10 += x_t0[ i         ]*filt[7];
      result20 += x_t0[(i   +1)%n1]*filt[7];
      result30 += x_t0[(i   +2)%n1]*filt[7];

      result01 += x_t1[(i+n1-1)%n1]*filt[7];
      result11 += x_t1[ i         ]*filt[7];
      result21 += x_t1[(i   +1)%n1]*filt[7];
      result31 += x_t1[(i   +2)%n1]*filt[7];

      result00 += x_t0[ i         ]*filt[8];
      result10 += x_t0[(i   +1)%n1]*filt[8];
      result20 += x_t0[(i   +2)%n1]*filt[8];
      result30 += x_t0[(i   +3)%n1]*filt[8];

      result01 += x_t1[ i         ]*filt[8];
      result11 += x_t1[(i   +1)%n1]*filt[8];
      result21 += x_t1[(i   +2)%n1]*filt[8];
      result31 += x_t1[(i   +3)%n1]*filt[8];

      result00 += x_t0[(i   +1)%n1]*filt[9];
      result10 += x_t0[(i   +2)%n1]*filt[9];
      result20 += x_t0[(i   +3)%n1]*filt[9];
      result30 += x_t0[(i   +4)%n1]*filt[9];

      result01 += x_t1[(i   +1)%n1]*filt[9];
      result11 += x_t1[(i   +2)%n1]*filt[9];
      result21 += x_t1[(i   +3)%n1]*filt[9];
      result31 += x_t1[(i   +4)%n1]*filt[9];

      result00 += x_t0[(i   +2)%n1]*filt[10];
      result10 += x_t0[(i   +3)%n1]*filt[10];
      result20 += x_t0[(i   +4)%n1]*filt[10];
      result30 += x_t0[(i   +5)%n1]*filt[10];

      result01 += x_t1[(i   +2)%n1]*filt[10];
      result11 += x_t1[(i   +3)%n1]*filt[10];
      result21 += x_t1[(i   +4)%n1]*filt[10];
      result31 += x_t1[(i   +5)%n1]*filt[10];

      result00 += x_t0[(i   +3)%n1]*filt[11];
      result10 += x_t0[(i   +4)%n1]*filt[11];
      result20 += x_t0[(i   +5)%n1]*filt[11];
      result30 += x_t0[(i   +6)%n1]*filt[11];

      result01 += x_t1[(i   +3)%n1]*filt[11];
      result11 += x_t1[(i   +4)%n1]*filt[11];
      result21 += x_t1[(i   +5)%n1]*filt[11];
      result31 += x_t1[(i   +6)%n1]*filt[11];

      result00 += x_t0[(i   +4)%n1]*filt[12];
      result10 += x_t0[(i   +5)%n1]*filt[12];
      result20 += x_t0[(i   +6)%n1]*filt[12];
      result30 += x_t0[(i   +7)%n1]*filt[12];

      result01 += x_t1[(i   +4)%n1]*filt[12];
      result11 += x_t1[(i   +5)%n1]*filt[12];
      result21 += x_t1[(i   +6)%n1]*filt[12];
      result31 += x_t1[(i   +7)%n1]*filt[12];

      result00 += x_t0[(i   +5)%n1]*filt[13];
      result10 += x_t0[(i   +6)%n1]*filt[13];
      result20 += x_t0[(i   +7)%n1]*filt[13];
      result30 += x_t0[(i   +8)%n1]*filt[13];

      result01 += x_t1[(i   +5)%n1]*filt[13];
      result11 += x_t1[(i   +6)%n1]*filt[13];
      result21 += x_t1[(i   +7)%n1]*filt[13];
      result31 += x_t1[(i   +8)%n1]*filt[13];

      result00 += x_t0[(i   +6)%n1]*filt[14];
      result10 += x_t0[(i   +7)%n1]*filt[14];
      result20 += x_t0[(i   +8)%n1]*filt[14];
      result30 += x_t0[(i   +9)%n1]*filt[14];

      result01 += x_t1[(i   +6)%n1]*filt[14];
      result11 += x_t1[(i   +7)%n1]*filt[14];
      result21 += x_t1[(i   +8)%n1]*filt[14];
      result31 += x_t1[(i   +9)%n1]*filt[14];

      result00 += x_t0[(i   +7)%n1]*filt[15];
      result10 += x_t0[(i   +8)%n1]*filt[15];
      result20 += x_t0[(i   +9)%n1]*filt[15];
      result30 += x_t0[(i   +10)%n1]*filt[15];

      result01 += x_t1[(i   +7)%n1]*filt[15];
      result11 += x_t1[(i   +8)%n1]*filt[15];
      result21 += x_t1[(i   +9)%n1]*filt[15];
      result31 += x_t1[(i   +10)%n1]*filt[15];

      y[i*ld + j] = result00;
      y[(i+1)*ld + j] = result10;
      y[(i+2)*ld + j] = result20;
      y[(i+3)*ld + j] = result30;
      y[i*ld + j+1] = result01;
      y[(i+1)*ld + j+1] = result11;
      y[(i+2)*ld + j+1] = result21;
      y[(i+3)*ld + j+1] = result31;
    }
    for(i=(n1/4)*4; i<n1; i++) {
      result0 = 0;
      result1 = 0;
      result0 += x_t0[(i+n1-8)%n1]*filt[0];
      result1 += x_t1[(i+n1-8)%n1]*filt[0];
      result0 += x_t0[(i+n1-7)%n1]*filt[1];
      result1 += x_t1[(i+n1-7)%n1]*filt[1];
      result0 += x_t0[(i+n1-6)%n1]*filt[2];
      result1 += x_t1[(i+n1-6)%n1]*filt[2];
      result0 += x_t0[(i+n1-5)%n1]*filt[3];
      result1 += x_t1[(i+n1-5)%n1]*filt[3];
      result0 += x_t0[(i+n1-4)%n1]*filt[4];
      result1 += x_t1[(i+n1-4)%n1]*filt[4];
      result0 += x_t0[(i+n1-3)%n1]*filt[5];
      result1 += x_t1[(i+n1-3)%n1]*filt[5];
      result0 += x_t0[(i+n1-2)%n1]*filt[6];
      result1 += x_t1[(i+n1-2)%n1]*filt[6];
      result0 += x_t0[(i+n1-1)%n1]*filt[7];
      result1 += x_t1[(i+n1-1)%n1]*filt[7];
      result0 += x_t0[ i         ]*filt[8];
      result1 += x_t1[ i         ]*filt[8];
      result0 += x_t0[(i   +1)%n1]*filt[9];
      result1 += x_t1[(i   +1)%n1]*filt[9];
      result0 += x_t0[(i   +2)%n1]*filt[10];
      result1 += x_t1[(i   +2)%n1]*filt[10];
      result0 += x_t0[(i   +3)%n1]*filt[11];
      result1 += x_t1[(i   +3)%n1]*filt[11];
      result0 += x_t0[(i   +4)%n1]*filt[12];
      result1 += x_t1[(i   +5)%n1]*filt[13];
      result0 += x_t0[(i   +6)%n1]*filt[14];
      result1 += x_t1[(i   +6)%n1]*filt[14];
      result0 += x_t0[(i   +7)%n1]*filt[15];
      result1 += x_t1[(i   +7)%n1]*filt[15];
      y[i*ld + j] = result0;
      y[i*ld + j+1] = result1;
    }
  }
  for(j=start + (ndat/2)*4; j<start+ndat; j++) {
    x_t0 = x+j*n1;
    for(i=0; i<n1; i+=4) {
      result00 = 0;
      result10 = 0;
      result20 = 0;
      result30 = 0;

      result00 += x_t0[(i+n1-8)%n1]*filt[0];
      result10 += x_t0[(i+n1-7)%n1]*filt[0];
      result20 += x_t0[(i+n1-6)%n1]*filt[0];
      result30 += x_t0[(i+n1-5)%n1]*filt[0];

      result00 += x_t0[(i+n1-7)%n1]*filt[1];
      result10 += x_t0[(i+n1-6)%n1]*filt[1];
      result20 += x_t0[(i+n1-5)%n1]*filt[1];
      result30 += x_t0[(i+n1-4)%n1]*filt[1];

      result00 += x_t0[(i+n1-6)%n1]*filt[2];
      result10 += x_t0[(i+n1-5)%n1]*filt[2];
      result20 += x_t0[(i+n1-4)%n1]*filt[2];
      result30 += x_t0[(i+n1-3)%n1]*filt[2];

      result00 += x_t0[(i+n1-5)%n1]*filt[3];
      result10 += x_t0[(i+n1-4)%n1]*filt[3];
      result20 += x_t0[(i+n1-3)%n1]*filt[3];
      result30 += x_t0[(i+n1-2)%n1]*filt[3];

      result00 += x_t0[(i+n1-4)%n1]*filt[4];
      result10 += x_t0[(i+n1-3)%n1]*filt[4];
      result20 += x_t0[(i+n1-2)%n1]*filt[4];
      result30 += x_t0[(i+n1-1)%n1]*filt[4];

      result00 += x_t0[(i+n1-3)%n1]*filt[5];
      result10 += x_t0[(i+n1-2)%n1]*filt[5];
      result20 += x_t0[(i+n1-1)%n1]*filt[5];
      result30 += x_t0[ i         ]*filt[5];

      result00 += x_t0[(i+n1-2)%n1]*filt[6];
      result10 += x_t0[(i+n1-1)%n1]*filt[6];
      result20 += x_t0[ i         ]*filt[6];
      result30 += x_t0[(i   +1)%n1]*filt[6];

      result00 += x_t0[(i+n1-1)%n1]*filt[7];
      result10 += x_t0[ i         ]*filt[7];
      result20 += x_t0[(i   +1)%n1]*filt[7];
      result30 += x_t0[(i   +2)%n1]*filt[7];

      result00 += x_t0[ i         ]*filt[8];
      result10 += x_t0[(i   +1)%n1]*filt[8];
      result20 += x_t0[(i   +2)%n1]*filt[8];
      result30 += x_t0[(i   +3)%n1]*filt[8];

      result00 += x_t0[(i   +1)%n1]*filt[9];
      result10 += x_t0[(i   +2)%n1]*filt[9];
      result20 += x_t0[(i   +3)%n1]*filt[9];
      result30 += x_t0[(i   +4)%n1]*filt[9];

      result00 += x_t0[(i   +2)%n1]*filt[10];
      result10 += x_t0[(i   +3)%n1]*filt[10];
      result20 += x_t0[(i   +4)%n1]*filt[10];
      result30 += x_t0[(i   +5)%n1]*filt[10];

      result00 += x_t0[(i   +3)%n1]*filt[11];
      result10 += x_t0[(i   +4)%n1]*filt[11];
      result20 += x_t0[(i   +5)%n1]*filt[11];
      result30 += x_t0[(i   +6)%n1]*filt[11];

      result00 += x_t0[(i   +4)%n1]*filt[12];
      result10 += x_t0[(i   +5)%n1]*filt[12];
      result20 += x_t0[(i   +6)%n1]*filt[12];
      result30 += x_t0[(i   +7)%n1]*filt[12];

      result00 += x_t0[(i   +5)%n1]*filt[13];
      result10 += x_t0[(i   +6)%n1]*filt[13];
      result20 += x_t0[(i   +7)%n1]*filt[13];
      result30 += x_t0[(i   +8)%n1]*filt[13];

      result00 += x_t0[(i   +6)%n1]*filt[14];
      result10 += x_t0[(i   +7)%n1]*filt[14];
      result20 += x_t0[(i   +8)%n1]*filt[14];
      result30 += x_t0[(i   +9)%n1]*filt[14];

      result00 += x_t0[(i   +7)%n1]*filt[15];
      result10 += x_t0[(i   +8)%n1]*filt[15];
      result20 += x_t0[(i   +9)%n1]*filt[15];
      result30 += x_t0[(i   +10)%n1]*filt[15];

      y[i*ld + j] = result00;
      y[(i+1)*ld + j] = result10;
      y[(i+2)*ld + j] = result20;
      y[(i+3)*ld + j] = result30;
    }
    for(i=(n1/4)*4; i<n1; i++) {
      result0 = 0;
      result0 += x_t0[(i+n1-8)%n1]*filt[0];
      result0 += x_t0[(i+n1-7)%n1]*filt[1];
      result0 += x_t0[(i+n1-6)%n1]*filt[2];
      result0 += x_t0[(i+n1-5)%n1]*filt[3];
      result0 += x_t0[(i+n1-4)%n1]*filt[4];
      result0 += x_t0[(i+n1-3)%n1]*filt[5];
      result0 += x_t0[(i+n1-2)%n1]*filt[6];
      result0 += x_t0[(i+n1-1)%n1]*filt[7];
      result0 += x_t0[ i         ]*filt[8];
      result0 += x_t0[(i   +1)%n1]*filt[9];
      result0 += x_t0[(i   +2)%n1]*filt[10];
      result0 += x_t0[(i   +3)%n1]*filt[11];
      result0 += x_t0[(i   +4)%n1]*filt[12];
      result0 += x_t0[(i   +6)%n1]*filt[14];
      result0 += x_t0[(i   +7)%n1]*filt[15];
      y[i*ld + j] = result0;
    }
  }
}



void magicfilter_partial2(unsigned int n1, unsigned int start, unsigned int ndat, unsigned int ld, const double *x, double *y) {
  int i,j,k;
  double result0,result1,result2,result3,result4,result5,result6,result7;//,result8;//,result9,result10,result11;//,result12,result13,result14,result15;
//  double result01,result11,result21,result31;//,result4,result5,result6,result7,result8,result9,result10,result11,result12,result13,result14,result15;
  const double *x_t0, *x_t1, *x_t2, *x_t3,*x_t4, *x_t5, *x_t6, *x_t7;//,*x_t8;//, *x_t9, *x_t10, *x_t11;
  double *y_t;
  x_t7 = x-n1;
  const double *fil_t = filt+8;
  for(j=start; j+7<start+ndat; j+=8 ) {
    x_t0=x_t7+n1;
    x_t1=x_t0+n1;
    x_t2=x_t1+n1;
    x_t3=x_t2+n1;
    x_t4=x_t3+n1;
    x_t5=x_t4+n1;
    x_t6=x_t5+n1;
    x_t7=x_t6+n1;
    y_t = y + j - ld + 8;

    for(i=0;i<8;i++){
      unsigned int i_t = i+n1;
      result0=0;
      result1=0;
      result2=0;
      result3=0;
      result4=0;
      result5=0;
      result6=0;
      result7=0;
      y_t += ld-8;
      for(k=-8; k<8; k++) {
        size_t o = (i_t+k)%n1;
        result0+=x_t0[o]*fil_t[k];
        result1+=x_t1[o]*fil_t[k];
        result2+=x_t2[o]*fil_t[k];
        result3+=x_t3[o]*fil_t[k];
        result4+=x_t4[o]*fil_t[k];
        result5+=x_t5[o]*fil_t[k];
        result6+=x_t6[o]*fil_t[k];
        result7+=x_t7[o]*fil_t[k];
      }
      *y_t++ = result0;
      *y_t++ = result1;
      *y_t++ = result2;
      *y_t++ = result3;
      *y_t++ = result4;
      *y_t++ = result5;
      *y_t++ = result6;
      *y_t++ = result7;
    }

    for(i=8; i<n1-8; i++) {
      result0=0;
      result1=0;
      result2=0;
      result3=0;
      result4=0;
      result5=0;
      result6=0;
      result7=0;
      y_t += ld-8;
      for(k=-8; k<8; k++) {
        size_t o = i+k;
        result0+=x_t0[o]*fil_t[k];
        result1+=x_t1[o]*fil_t[k];
        result2+=x_t2[o]*fil_t[k];
        result3+=x_t3[o]*fil_t[k];
        result4+=x_t4[o]*fil_t[k];
        result5+=x_t5[o]*fil_t[k];
        result6+=x_t6[o]*fil_t[k];
        result7+=x_t7[o]*fil_t[k];
      }
      *y_t++ = result0;
      *y_t++ = result1;
      *y_t++ = result2;
      *y_t++ = result3;
      *y_t++ = result4;
      *y_t++ = result5;
      *y_t++ = result6;
      *y_t++ = result7;
    }

    for(i=n1-8;i<n1;i++){
      unsigned int i_t = i+n1;
      result0=0;
      result1=0;
      result2=0;
      result3=0;
      result4=0;
      result5=0;
      result6=0;
      result7=0;
      y_t += ld-8;
      for(k=-8; k<8; k++) {
        size_t o = (i_t+k)%n1;
        result0+=x_t0[o]*fil_t[k];
        result1+=x_t1[o]*fil_t[k];
        result2+=x_t2[o]*fil_t[k];
        result3+=x_t3[o]*fil_t[k];
        result4+=x_t4[o]*fil_t[k];
        result5+=x_t5[o]*fil_t[k];
        result6+=x_t6[o]*fil_t[k];
        result7+=x_t7[o]*fil_t[k];
      }
      *y_t++ = result0;
      *y_t++ = result1;
      *y_t++ = result2;
      *y_t++ = result3;
      *y_t++ = result4;
      *y_t++ = result5;
      *y_t++ = result6;
      *y_t++ = result7;
    }
  }
  for(j=start + (ndat/8)*8; j<start+ndat; j++) {
//    printf("plop\n");
    x_t0 = x+j*n1;
    for(i=0; i<n1; i++) {
      result0 = 0;
      unsigned int i_t = i+n1;

      for(k=-8; k<8; k++) {
        result0 += x_t0[(i_t+k)%n1]*fil_t[k];
      }

      y[i*ld + j] = result0;
    }
  }
}

struct magicfilter_params {
  unsigned int n1;
  unsigned int start;
  unsigned int ndat;
  unsigned int ld;
  unsigned long long int start_date;
  unsigned long long int stop_date;
  const double *x;
  double *y;
};


struct magicfilter_params ** params;

void init_thread_engine_(){
  char *argv[]={"conv_check", "-t", "topology_muscade.cfg", "-p", "placement_muscade.cfg", "-v", "1"};
  int argc = 7;
  int i;

  get_engine_opt( argc, argv);

  init_thread_engine( );

  params = (struct magicfilter_params **)malloc( (1 + engine_params->thread_number) * sizeof(struct magicfilter_params *));
  for(i=0; i < 1 + engine_params->thread_number; i++) {
    params[i] = (struct magicfilter_params *)malloc( sizeof(struct magicfilter_params) );
  }
}

void * magicfilter_part_wrapper(void *p) {
  struct magicfilter_params * params;
  params = ( struct magicfilter_params *) p;
//  rdtscll(params->start_date);
  magicfilter_partial2(params->n1, params->start, params->ndat, params->ld, params->x, params->y);
//  rdtscll(params->stop_date);
  return NULL;
}

void magicfilter1d_d_par_(unsigned int *n1, unsigned int *ndat, double *x, double *y) {
  int i,j;
  size_t slice;
  unsigned int steps = 1;
  slice = *ndat/(1 + engine_params->thread_number);
  while( (slice* *n1 *sizeof(double)*NB_SYB) > CACHE_SIZE ) {
    slice /= 2;
    steps++;
  }
  slice = *ndat/(steps*(1 + engine_params->thread_number));
//  printf("steps : %d, slice : %lu\n",steps,slice);
  for(i=0; i<steps-1; i++) {
    for(j=0; j<(1 + engine_params->thread_number); j++) {
      params[j]->n1 = *n1;
      params[j]->start = (i*(1 + engine_params->thread_number)+j)*slice;
      params[j]->ndat = slice;
      params[j]->ld = *ndat;
      params[j]->x = x;
      params[j]->y = y;
    }
    run_bench(magicfilter_part_wrapper,magicfilter_part_wrapper, (void **)params );
  }
  for(j=0; j<engine_params->thread_number; j++) {
    params[j]->n1 = *n1;
    params[j]->start = ((steps-1)*(1 + engine_params->thread_number)+j)*slice;
    params[j]->ndat = slice;
    params[j]->ld = *ndat;
    params[j]->x = x;
    params[j]->y = y;
  }
  params[engine_params->thread_number]->n1 = *n1;
  params[engine_params->thread_number]->start = ((steps-1)*(1 + engine_params->thread_number)+engine_params->thread_number)*slice;
  params[engine_params->thread_number]->ndat = slice + (*ndat % slice);
  params[engine_params->thread_number]->ld = *ndat;
  params[engine_params->thread_number]->x = x;
  params[engine_params->thread_number]->y = y;
//  unsigned long long int start, stop;
//  rdtscll(start);
  run_bench(magicfilter_part_wrapper,magicfilter_part_wrapper, (void **)params );
//  rdtscll(stop);
//  printf("start : %llu, stop : %llu\n",start,stop-start);
//  for(j=0; j<(1 + engine_params->thread_number); j++)
//    printf("  start : %llu, stop %llu\n", params[j]->start_date-start, params[j]->stop_date-start);
}

void magicfilter1d_d_seq_(unsigned int *n1, unsigned int *ndat, double *x, double *y) {
   int i,j;
   size_t slice;
   unsigned int steps = 1;
   slice = *ndat/NB_CORE;
   while( (slice* *n1 *sizeof(double)*NB_SYB) > CACHE_SIZE ) {
     slice /= 2;
     steps++;
   }
   slice = *ndat/(steps*NB_CORE);
   printf("steps : %d, slice : %lu\n",steps,slice);
   for(i=0; i<steps-1; i++) {
     for(j=0; j<NB_CORE; j++) {
       magicfilter_partial2(*n1, (i*NB_CORE+j)*slice, slice, *ndat, x, y);
     }
   }
   for(j=0; j<NB_CORE-1; j++) {
     magicfilter_partial2(*n1, ((steps-1)*NB_CORE+j)*slice, slice, *ndat, x, y);
   }
   magicfilter_partial2(*n1, ((steps-1)*NB_CORE+NB_CORE-1)*slice, slice+ (*ndat % slice), *ndat, x, y);
}
