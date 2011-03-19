//! @file
//! ???
//!
//! @author
//!    Copyright (C) 2009-2011 BigDFT group 
//!    This file is distributed under the terms of the
//!    GNU General Public License, see ~/COPYING file
//!    or http://www.gnu.org/copyleft/gpl.txt .
//!    For the list of contributors, see ~/AUTHORS 


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <emmintrin.h>
#include <pmmintrin.h>
#include <assert.h>
#include <math.h>


const double filt[] __attribute__ ((aligned (16))) = { 8.4334247333529341094733325815816e-7,
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

void nanosec(unsigned long long int * t){
  struct timespec time;
  clock_gettime(CLOCK_REALTIME, &time);
  *t = time.tv_sec;
  *t *= 1000000000;
  *t += time.tv_nsec;
}
#define FILTER_SIZE 16
#define BUFFER_WIDTH 1024
#define BUFFER_DEPTH 64


void conv_interleaved_ref(unsigned int n, double const * source, double * dest){
  unsigned int i,j;
  double tmp0,tmp1;
  for(j=0;j<n;j++,source+=2) {
    tmp0 = 0.0;
    tmp1 = 0.0;
    for(i=0;i<FILTER_SIZE;i++) {
      tmp0 += source[2*i]*filt[i];
      tmp1 += source[2*i+1]*filt[i];
    }
    dest[2*j] = tmp0;
    dest[2*j+1] = tmp1;
  }
}

void interleave(unsigned int n, double const * source, double * dest){
  int i,j;
  int ndat=2;
  for(j=0;j<ndat;j++) {
    for(i=0;i<n;i++){
      dest[i*ndat] = *source++;
    }
    dest++;
  }
}
#define microloop_nostore(source,D,F,B0,B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,B12,B13,B14,B15) \
  D = _mm_load_pd(source);\
  F = _mm_load1_pd(filt); \
  B0 = _mm_mul_pd(F,D);\
  F = _mm_load1_pd(filt+1);\
  B15 = _mm_add_pd(B15,_mm_mul_pd(F,D));\
  F = _mm_load1_pd(filt+2);\
  B14 = _mm_add_pd(B14,_mm_mul_pd(F,D));\
  F = _mm_load1_pd(filt+3);\
  B13 = _mm_add_pd(B13,_mm_mul_pd(F,D));\
  F = _mm_load1_pd(filt+4);\
  B12 = _mm_add_pd(B12,_mm_mul_pd(F,D));\
  F = _mm_load1_pd(filt+5);\
  B11 = _mm_add_pd(B11,_mm_mul_pd(F,D));\
  F = _mm_load1_pd(filt+6);\
  B10 = _mm_add_pd(B10,_mm_mul_pd(F,D));\
  F = _mm_load1_pd(filt+7);\
  B9 = _mm_add_pd(B9,_mm_mul_pd(F,D));\
  F = _mm_load1_pd(filt+8);\
  B8 = _mm_add_pd(B8,_mm_mul_pd(F,D));\
  F = _mm_load1_pd(filt+9);\
  B7 = _mm_add_pd(B7,_mm_mul_pd(F,D));\
  F = _mm_load1_pd(filt+10);\
  B6 = _mm_add_pd(B6,_mm_mul_pd(F,D));\
  F = _mm_load1_pd(filt+11);\
  B5 = _mm_add_pd(B5,_mm_mul_pd(F,D));\
  F = _mm_load1_pd(filt+12);\
  B4 = _mm_add_pd(B4,_mm_mul_pd(F,D));\
  F = _mm_load1_pd(filt+13);\
  B3 = _mm_add_pd(B3,_mm_mul_pd(F,D));\
  F = _mm_load1_pd(filt+14);\
  B2 = _mm_add_pd(B2,_mm_mul_pd(F,D));\
  F = _mm_load1_pd(filt+15);\
  B1 = _mm_add_pd(B1,_mm_mul_pd(F,D));\

inline void conv_interleaved(unsigned int n, double const * source, double * dest){
  int i,j;
  __v2df S0,S1,S2,S3,S4,S5,S6,S7;
  __v2df S8,S9,S10,S11,S12,S13,S14,S15;
  __v2df F;
  __v2df D;

  microloop_nostore(source,D,F,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15);
  microloop_nostore(source+2,D,F,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0);
  microloop_nostore(source+4,D,F,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1);
  microloop_nostore(source+6,D,F,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2);
  microloop_nostore(source+8,D,F,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3);
  microloop_nostore(source+10,D,F,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4);
  microloop_nostore(source+12,D,F,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5);
  microloop_nostore(source+14,D,F,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6);
  microloop_nostore(source+16,D,F,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7);
  microloop_nostore(source+18,D,F,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8);
  microloop_nostore(source+20,D,F,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9);
  microloop_nostore(source+22,D,F,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10);
  microloop_nostore(source+24,D,F,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11);
  microloop_nostore(source+26,D,F,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12);
  microloop_nostore(source+28,D,F,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13);
  microloop_nostore(source+30,D,F,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14);
  source+=32;
  i+=16;
  while(i<n) {
    _mm_store_pd(dest,S0);
    microloop_nostore(source,D,F,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15);
    _mm_store_pd(dest+2,S1);
    microloop_nostore(source+2,D,F,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0);
    _mm_store_pd(dest+4,S2);
    microloop_nostore(source+4,D,F,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1);
    _mm_store_pd(dest+6,S3);
    microloop_nostore(source+6,D,F,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2);
    _mm_store_pd(dest+8,S4);
    microloop_nostore(source+8,D,F,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3);
    _mm_store_pd(dest+10,S5);
    microloop_nostore(source+10,D,F,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4);
    _mm_store_pd(dest+12,S6);
    microloop_nostore(source+12,D,F,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5);
    _mm_store_pd(dest+14,S7);
    microloop_nostore(source+14,D,F,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6);
    _mm_store_pd(dest+16,S8);
    microloop_nostore(source+16,D,F,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7);
    _mm_store_pd(dest+18,S9);
    microloop_nostore(source+18,D,F,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8);
    _mm_store_pd(dest+20,S10);
    microloop_nostore(source+20,D,F,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9);
    _mm_store_pd(dest+22,S11);
    microloop_nostore(source+22,D,F,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10);
    _mm_store_pd(dest+24,S12);
    microloop_nostore(source+24,D,F,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11);
    _mm_store_pd(dest+26,S13);
    microloop_nostore(source+26,D,F,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12);
    _mm_store_pd(dest+28,S14);
    microloop_nostore(source+28,D,F,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13);
    _mm_store_pd(dest+30,S15);
    microloop_nostore(source+30,D,F,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14);
    source+=32;
    dest+=32;
    i += 16;
  }
  _mm_store_pd(dest,S0);
  microloop_nostore(source,D,F,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15);
  _mm_store_pd(dest+2,S1);
  microloop_nostore(source+2,D,F,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0);
  _mm_store_pd(dest+4,S2);
  microloop_nostore(source+4,D,F,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1);
  _mm_store_pd(dest+6,S3);
  microloop_nostore(source+6,D,F,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2);
  _mm_store_pd(dest+8,S4);
  microloop_nostore(source+8,D,F,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3);
  _mm_store_pd(dest+10,S5);
  microloop_nostore(source+10,D,F,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4);
  _mm_store_pd(dest+12,S6);
  microloop_nostore(source+12,D,F,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5);
  _mm_store_pd(dest+14,S7);
  microloop_nostore(source+14,D,F,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6);
  _mm_store_pd(dest+16,S8);
  microloop_nostore(source+16,D,F,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7);
  _mm_store_pd(dest+18,S9);
  microloop_nostore(source+18,D,F,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8);
  _mm_store_pd(dest+20,S10);
  microloop_nostore(source+20,D,F,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9);
  _mm_store_pd(dest+22,S11);
  microloop_nostore(source+22,D,F,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10);
  _mm_store_pd(dest+24,S12);
  microloop_nostore(source+24,D,F,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11);
  _mm_store_pd(dest+26,S13);
  microloop_nostore(source+26,D,F,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12);
  _mm_store_pd(dest+28,S14);
  microloop_nostore(source+28,D,F,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13);
  _mm_store_pd(dest+30,S15);
}

#define FLOP (BUFFER_WIDTH*BUFFER_DEPTH)*4*FILTER_SIZE
int main(void) {
  double a[BUFFER_WIDTH][2*(BUFFER_DEPTH+FILTER_SIZE)] __attribute__ ((aligned (16)));
  double b[BUFFER_WIDTH][2*BUFFER_DEPTH] __attribute__ ((aligned (16)));
  double c[BUFFER_WIDTH][2*BUFFER_DEPTH] __attribute__ ((aligned (16)));
  unsigned int i,j,k;
  for(i=0; i<BUFFER_WIDTH;i++) {
    for(j=0;j<2*(BUFFER_DEPTH+FILTER_SIZE);j++){
      a[i][j] = rand()/(double)RAND_MAX;
    }
  }
  for(i=0; i<BUFFER_WIDTH;i++) {
    for(j=0;j<2*BUFFER_DEPTH;j++){
      b[i][j] =0.0;
    }
  }
  unsigned long long int t1,t2;
  nanosec(&t1);
  for(i=0; i<BUFFER_WIDTH;i++) {
    conv_interleaved(BUFFER_DEPTH,&(a[i][0]),&(b[i][0]));
  }
  nanosec(&t2);
  for(i=0; i<BUFFER_WIDTH;i++) {
    conv_interleaved_ref(BUFFER_DEPTH,&(a[i][0]),&(c[i][0]));
  }
  double br=0;
  for(i=0; i<BUFFER_WIDTH;i++) {
    for(j=0;j<2*BUFFER_DEPTH;j++) {
      br+=b[i][j];
      if(fabs(b[i][j]-c[i][j])>1e-10){
        printf("error %u %u: %lf != %lf!\n",i,j , b[i][j], c[i][j]);
      }
    }
  }
  printf("result %lf, duration %llu ns, FLOP %d, GFLOPS %lf\n", br, t2-t1, FLOP, (double)FLOP/(float)(t2-t1));
}
