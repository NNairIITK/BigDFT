//! @file
//!  Convolutions SSE (vectorisation)
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
#include <string.h>
#include <time.h>
#include <papi.h>

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
const double filt_u[] __attribute__ ((aligned (16))) = { 2.72734492911979659657715313017228e-6,
                        8.4334247333529341094733325815816e-7,
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
#define BUFFER_WIDTH (128*132)
#define BUFFER_DEPTH 128

#define store2_reverse(dest,offset,S) \
    _mm_storel_pd(dest+offset,S);\
    _mm_storeh_pd(dest+offset+n,S);

#define store2(dest,offset,S) \
    _mm_store_pd(dest+ndat*offset,S);

#define microloop_nostore(source0,source1,offset,D,F,B0,B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,B12,B13,B14,B15) \
  D = _mm_set_pd(*(source1+offset),*(source0+offset));\
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

#define microloop_nostore_reverse(source,offset,D,F,B0,B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,B12,B13,B14,B15) \
  D = _mm_load_pd(source+offset*ndat);\
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


inline void conv_gliding(unsigned int n,unsigned int ndat, double const * source0, double const * source1, double * dest){
  int i;
  __m128d S0,S1,S2,S3,S4,S5,S6,S7;
  __m128d S8,S9,S10,S11,S12,S13,S14,S15;
  __m128d F;
  __m128d D;
  S0=S1=S2=S3=S4=S5=S6=S7=S8=S9=S10=S11=S12=S13=S14=S15=_mm_setzero_pd();

  microloop_nostore(source0,source1,0,D,F,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15);
  microloop_nostore(source0,source1,1,D,F,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0);
  microloop_nostore(source0,source1,2,D,F,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1);
  microloop_nostore(source0,source1,3,D,F,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2);
  microloop_nostore(source0,source1,4,D,F,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3);
  microloop_nostore(source0,source1,5,D,F,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4);
  microloop_nostore(source0,source1,6,D,F,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5);
  microloop_nostore(source0,source1,7,D,F,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6);
  microloop_nostore(source0,source1,8,D,F,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7);
  microloop_nostore(source0,source1,9,D,F,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8);
  microloop_nostore(source0,source1,10,D,F,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9);
  microloop_nostore(source0,source1,11,D,F,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10);
  microloop_nostore(source0,source1,12,D,F,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11);
  microloop_nostore(source0,source1,13,D,F,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12);
  microloop_nostore(source0,source1,14,D,F,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13);
  microloop_nostore(source0,source1,15,D,F,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14);
  source0 += 16;
  source1 += 16;
  i=16;
  while(i<n) {
    store2(dest,0,S0);
    microloop_nostore(source0,source1,0,D,F,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15);
    store2(dest,1,S1);
    microloop_nostore(source0,source1,1,D,F,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0);
    store2(dest,2,S2);
    microloop_nostore(source0,source1,2,D,F,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1);
    store2(dest,3,S3);
    microloop_nostore(source0,source1,3,D,F,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2);
    store2(dest,4,S4);
    microloop_nostore(source0,source1,4,D,F,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3);
    store2(dest,5,S5);
    microloop_nostore(source0,source1,5,D,F,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4);
    store2(dest,6,S6);
    microloop_nostore(source0,source1,6,D,F,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5);
    store2(dest,7,S7);
    microloop_nostore(source0,source1,7,D,F,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6);
    store2(dest,8,S8);
    microloop_nostore(source0,source1,8,D,F,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7);
    store2(dest,9,S9);
    microloop_nostore(source0,source1,9,D,F,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8);
    store2(dest,10,S10);
    microloop_nostore(source0,source1,10,D,F,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9);
    store2(dest,11,S11);
    microloop_nostore(source0,source1,11,D,F,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10);
    store2(dest,12,S12);
    microloop_nostore(source0,source1,12,D,F,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11);
    store2(dest,13,S13);
    microloop_nostore(source0,source1,13,D,F,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12);
    store2(dest,14,S14);
    microloop_nostore(source0,source1,14,D,F,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13);
    store2(dest,15,S15);
    microloop_nostore(source0,source1,15,D,F,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14);
    source0+=16;
    dest+=16*ndat;
    source1+=16;
    i += 16;
  }
  store2(dest,0,S0);
  microloop_nostore(source0,source1,0,D,F,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15);
  store2(dest,1,S1);
  microloop_nostore(source0,source1,1,D,F,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0);
  store2(dest,2,S2);
  microloop_nostore(source0,source1,2,D,F,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1);
  store2(dest,3,S3);
  microloop_nostore(source0,source1,3,D,F,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2);
  store2(dest,4,S4);
  microloop_nostore(source0,source1,4,D,F,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3);
  store2(dest,5,S5);
  microloop_nostore(source0,source1,5,D,F,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4);
  store2(dest,6,S6);
  microloop_nostore(source0,source1,6,D,F,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5);
  store2(dest,7,S7);
  microloop_nostore(source0,source1,7,D,F,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6);
  store2(dest,8,S8);
  microloop_nostore(source0,source1,8,D,F,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7);
  store2(dest,9,S9);
  microloop_nostore(source0,source1,9,D,F,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8);
  store2(dest,10,S10);
  microloop_nostore(source0,source1,10,D,F,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9);
  store2(dest,11,S11);
  microloop_nostore(source0,source1,11,D,F,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10);
  store2(dest,12,S12);
  microloop_nostore(source0,source1,12,D,F,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11);
  store2(dest,13,S13);
  microloop_nostore(source0,source1,13,D,F,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12);
  store2(dest,14,S14);
  microloop_nostore(source0,source1,14,D,F,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13);
  store2(dest,15,S15);
}

inline void conv_gliding_reverse(unsigned int n,unsigned int ndat, double const * source, double * dest){
  int i;
  __m128d S0,S1,S2,S3,S4,S5,S6,S7;
  __m128d S8,S9,S10,S11,S12,S13,S14,S15;
  __m128d F;
  __m128d D;
  S0=S1=S2=S3=S4=S5=S6=S7=S8=S9=S10=S11=S12=S13=S14=S15=_mm_setzero_pd();

  microloop_nostore_reverse(source,0,D,F,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15);
  microloop_nostore_reverse(source,1,D,F,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0);
  microloop_nostore_reverse(source,2,D,F,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1);
  microloop_nostore_reverse(source,3,D,F,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2);
  microloop_nostore_reverse(source,4,D,F,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3);
  microloop_nostore_reverse(source,5,D,F,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4);
  microloop_nostore_reverse(source,6,D,F,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5);
  microloop_nostore_reverse(source,7,D,F,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6);
  microloop_nostore_reverse(source,8,D,F,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7);
  microloop_nostore_reverse(source,9,D,F,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8);
  microloop_nostore_reverse(source,10,D,F,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9);
  microloop_nostore_reverse(source,11,D,F,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10);
  microloop_nostore_reverse(source,12,D,F,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11);
  microloop_nostore_reverse(source,13,D,F,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12);
  microloop_nostore_reverse(source,14,D,F,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13);
  microloop_nostore_reverse(source,15,D,F,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14);
  source += 16*ndat;
  i=16;
  while(i<n) {
    store2_reverse(dest,0,S0);
    microloop_nostore_reverse(source,0,D,F,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15);
    store2_reverse(dest,1,S1);
    microloop_nostore_reverse(source,1,D,F,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0);
    store2_reverse(dest,2,S2);
    microloop_nostore_reverse(source,2,D,F,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1);
    store2_reverse(dest,3,S3);
    microloop_nostore_reverse(source,3,D,F,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2);
    store2_reverse(dest,4,S4);
    microloop_nostore_reverse(source,4,D,F,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3);
    store2_reverse(dest,5,S5);
    microloop_nostore_reverse(source,5,D,F,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4);
    store2_reverse(dest,6,S6);
    microloop_nostore_reverse(source,6,D,F,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5);
    store2_reverse(dest,7,S7);
    microloop_nostore_reverse(source,7,D,F,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6);
    store2_reverse(dest,8,S8);
    microloop_nostore_reverse(source,8,D,F,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7);
    store2_reverse(dest,9,S9);
    microloop_nostore_reverse(source,9,D,F,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8);
    store2_reverse(dest,10,S10);
    microloop_nostore_reverse(source,10,D,F,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9);
    store2_reverse(dest,11,S11);
    microloop_nostore_reverse(source,11,D,F,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10);
    store2_reverse(dest,12,S12);
    microloop_nostore_reverse(source,12,D,F,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11);
    store2_reverse(dest,13,S13);
    microloop_nostore_reverse(source,13,D,F,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12);
    store2_reverse(dest,14,S14);
    microloop_nostore_reverse(source,14,D,F,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13);
    store2_reverse(dest,15,S15);
    microloop_nostore_reverse(source,15,D,F,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14);
    source+=16*ndat;
    dest+=16;
    i += 16;
  }
  store2_reverse(dest,0,S0);
  microloop_nostore_reverse(source,0,D,F,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15);
  store2_reverse(dest,1,S1);
  microloop_nostore_reverse(source,1,D,F,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0);
  store2_reverse(dest,2,S2);
  microloop_nostore_reverse(source,2,D,F,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1);
  store2_reverse(dest,3,S3);
  microloop_nostore_reverse(source,3,D,F,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2);
  store2_reverse(dest,4,S4);
  microloop_nostore_reverse(source,4,D,F,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3);
  store2_reverse(dest,5,S5);
  microloop_nostore_reverse(source,5,D,F,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4);
  store2_reverse(dest,6,S6);
  microloop_nostore_reverse(source,6,D,F,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5);
  store2_reverse(dest,7,S7);
  microloop_nostore_reverse(source,7,D,F,S7,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6);
  store2_reverse(dest,8,S8);
  microloop_nostore_reverse(source,8,D,F,S8,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7);
  store2_reverse(dest,9,S9);
  microloop_nostore_reverse(source,9,D,F,S9,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8);
  store2_reverse(dest,10,S10);
  microloop_nostore_reverse(source,10,D,F,S10,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9);
  store2_reverse(dest,11,S11);
  microloop_nostore_reverse(source,11,D,F,S11,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10);
  store2_reverse(dest,12,S12);
  microloop_nostore_reverse(source,12,D,F,S12,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11);
  store2_reverse(dest,13,S13);
  microloop_nostore_reverse(source,13,D,F,S13,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12);
  store2_reverse(dest,14,S14);
  microloop_nostore_reverse(source,14,D,F,S14,S15,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13);
  store2_reverse(dest,15,S15);
}


inline void conv_gliding_t(size_t n,size_t ndat, double const * source, double * dest){
  unsigned int i=0;
  do{
    conv_gliding(n,ndat,source,source+n+FILTER_SIZE,dest);
    source+=FILTER_SIZE*2+2*n;
    i+=2;
    dest+=2;
  } while (i<ndat);
}

inline void conv_gliding_t_reverse(size_t n,size_t ndat, double const * source, double * dest){
  unsigned int i=0;
  do{
    conv_gliding_reverse(n,ndat,source,dest);
    source+=2;
    i+=2;
    dest+=2*n;
  } while (i<ndat);
}
#define conv_2x2_block_fused(offset_filter,offset_source,d00,d10) \
FA = _mm_load_pd(filt+offset_filter);\
d00 = _mm_load_pd(source0+offset_source);\
S00 = _mm_add_pd(S00,_mm_mul_pd(d00,FA));\
FU = _mm_load_pd(filt_u+offset_filter);\
d10 = _mm_load_pd(source1+offset_source);\
S01 = _mm_add_pd(S01,_mm_mul_pd(d00,FU));\
S11 = _mm_add_pd(S11,_mm_mul_pd(d10,FU));\
S10 = _mm_add_pd(S10,_mm_mul_pd(d10,FA));

#define conv_2_block_fused(offset_filter,offset_source,d0) \
FA = _mm_load_pd(filt+offset_filter);\
d0 = _mm_load_pd(source+offset_source);\
S0 = _mm_add_pd(S0,_mm_mul_pd(d0,FA));\
FU = _mm_load_pd(filt_u+offset_filter);\
S1 = _mm_add_pd(S1,_mm_mul_pd(d0,FU));

#define conv_2_block(filter,offset_source,d0) \
F = _mm_load_pd(filter);\
d0 = _mm_load_pd(source+offset_source);\
S0 = _mm_add_pd(S0,_mm_mul_pd(d0,F));

#define conv_4_block_fused(offset_filter,offset_source,d0,d1) \
FA = _mm_load_pd(filt+offset_filter);\
S0 = _mm_add_pd(S0,_mm_mul_pd(d0,FA));\
FU = _mm_load_pd(filt_u+offset_filter);\
S1 = _mm_add_pd(S1,_mm_mul_pd(d0,FU));\
d1 = _mm_load_pd(source+offset_source);\
S2 = _mm_add_pd(S2,_mm_mul_pd(d1,FA));\
S3 = _mm_add_pd(S3,_mm_mul_pd(d1,FU));

#define conv_4_block(filter,offset_source,d0,d1) \
F = _mm_load_pd(filter);\
S0 = _mm_add_pd(S0,_mm_mul_pd(d0,F));\
d1 = _mm_load_pd(source+offset_source);\
S1 = _mm_add_pd(S1,_mm_mul_pd(d1,F));

#define conv_6_block_fused(offset_filter,offset_source,d0,d1,d2) \
FA = _mm_load_pd(filt+offset_filter);\
S0 = _mm_add_pd(S0,_mm_mul_pd(d0,FA));\
FU = _mm_load_pd(filt_u+offset_filter);\
S1 = _mm_add_pd(S1,_mm_mul_pd(d0,FU));\
S2 = _mm_add_pd(S2,_mm_mul_pd(d1,FA));\
S3 = _mm_add_pd(S3,_mm_mul_pd(d1,FU));\
d2 = _mm_load_pd(source+offset_source);\
S4 = _mm_add_pd(S4,_mm_mul_pd(d2,FA));\
S5 = _mm_add_pd(S5,_mm_mul_pd(d2,FU));

#define conv_6_block(filter,offset_source,d0,d1,d2) \
F = _mm_load_pd(filter);\
S0 = _mm_add_pd(S0,_mm_mul_pd(d0,F));\
S1 = _mm_add_pd(S1,_mm_mul_pd(d1,F));\
d2 = _mm_load_pd(source+offset_source);\
S2 = _mm_add_pd(S2,_mm_mul_pd(d2,F));

#define conv_4x2_block_fused(offset_filter,offset_source,d00,d10,d20,d30) \
FA = _mm_load_pd(filt+offset_filter);\
d00 = _mm_load_pd(source0+offset_source);\
S00 = _mm_add_pd(S00,_mm_mul_pd(d00,FA));\
FU = _mm_load_pd(filt_u+offset_filter);\
d10 = _mm_load_pd(source1+offset_source);\
S01 = _mm_add_pd(S01,_mm_mul_pd(d00,FU));\
S11 = _mm_add_pd(S11,_mm_mul_pd(d10,FU));\
d20 = _mm_load_pd(source2+offset_source);\
S10 = _mm_add_pd(S10,_mm_mul_pd(d10,FA));\
d30 = _mm_load_pd(source3+offset_source);\
S20 = _mm_add_pd(S20,_mm_mul_pd(d20,FA));\
S30 = _mm_add_pd(S30,_mm_mul_pd(d30,FA));\
S31 = _mm_add_pd(S31,_mm_mul_pd(d30,FU));\
S21 = _mm_add_pd(S21,_mm_mul_pd(d20,FU));

#define conv_2x4_block_fused(offset_filter,offset_source,d00,d01,d10,d11) \
FA = _mm_load_pd(filt+offset_filter);\
S00 = _mm_add_pd(S00,_mm_mul_pd(d00,FA));\
FU = _mm_load_pd(filt_u+offset_filter);\
S10 = _mm_add_pd(S10,_mm_mul_pd(d10,FA));\
S11 = _mm_add_pd(S11,_mm_mul_pd(d10,FU));\
d01 = _mm_load_pd(source0+offset_source);\
S01 = _mm_add_pd(S01,_mm_mul_pd(d00,FU));\
d11 = _mm_load_pd(source1+offset_source);\
S03 = _mm_add_pd(S03,_mm_mul_pd(d01,FU));\
S13 = _mm_add_pd(S13,_mm_mul_pd(d11,FU));\
S02 = _mm_add_pd(S02,_mm_mul_pd(d01,FA));\
S12 = _mm_add_pd(S12,_mm_mul_pd(d11,FA));

#define conv_8_block_fused(offset_filter,offset_source,d0,d1,d2,d3) \
FA = _mm_load_pd(filt+offset_filter);\
S0 = _mm_add_pd(S0,_mm_mul_pd(d0,FA));\
S2 = _mm_add_pd(S2,_mm_mul_pd(d1,FA));\
FU = _mm_load_pd(filt_u+offset_filter);\
S4 = _mm_add_pd(S4,_mm_mul_pd(d2,FA));\
S1 = _mm_add_pd(S1,_mm_mul_pd(d0,FU));\
d3 = _mm_load_pd(source+offset_source);\
S3 = _mm_add_pd(S3,_mm_mul_pd(d1,FU));\
S5 = _mm_add_pd(S5,_mm_mul_pd(d2,FU));\
S7 = _mm_add_pd(S7,_mm_mul_pd(d3,FU));\
S6 = _mm_add_pd(S6,_mm_mul_pd(d3,FA));

#define conv_8_block_fused_t(offset_filter,offset_source,d0,d1,d2,d3) \
FA = _mm_load_pd(filt+offset_filter);\
S0 = _mm_add_pd(S0,_mm_mul_pd(d0,FA));\
S2 = _mm_add_pd(S2,_mm_mul_pd(d1,FA));\
FU = _mm_load_pd(filt_u+offset_filter);\
S4 = _mm_add_pd(S4,_mm_mul_pd(d2,FA));\
S1 = _mm_add_pd(S1,_mm_mul_pd(d0,FU));\
d3 = _mm_set_pd(*(source+(offset_source+1)*ndat),*(source+offset_source*ndat));\
S3 = _mm_add_pd(S3,_mm_mul_pd(d1,FU));\
S5 = _mm_add_pd(S5,_mm_mul_pd(d2,FU));\
S7 = _mm_add_pd(S7,_mm_mul_pd(d3,FU));\
S6 = _mm_add_pd(S6,_mm_mul_pd(d3,FA));

#define conv_8_block(filter,offset_source,d0,d1,d2,d3) \
F = _mm_load_pd(filter);\
S0 = _mm_add_pd(S0,_mm_mul_pd(d0,F));\
S1 = _mm_add_pd(S1,_mm_mul_pd(d1,F));\
S2 = _mm_add_pd(S2,_mm_mul_pd(d2,F));\
d3 = _mm_load_pd(source+offset_source);\
S3 = _mm_add_pd(S3,_mm_mul_pd(d3,F));

#define conv_10_block_fused(offset_filter,offset_source,d0,d1,d2,d3,d4) \
FA = _mm_load_pd(filt+offset_filter);\
S0 = _mm_add_pd(S0,_mm_mul_pd(d0,FA));\
S2 = _mm_add_pd(S2,_mm_mul_pd(d1,FA));\
FU = _mm_load_pd(filt_u+offset_filter);\
S4 = _mm_add_pd(S4,_mm_mul_pd(d2,FA));\
S6 = _mm_add_pd(S6,_mm_mul_pd(d3,FA));\
S1 = _mm_add_pd(S1,_mm_mul_pd(d0,FU));\
d4 = _mm_load_pd(source+offset_source);\
S3 = _mm_add_pd(S3,_mm_mul_pd(d1,FU));\
S5 = _mm_add_pd(S5,_mm_mul_pd(d2,FU));\
S7 = _mm_add_pd(S7,_mm_mul_pd(d3,FU));\
S9 = _mm_add_pd(S9,_mm_mul_pd(d4,FU));\
S8 = _mm_add_pd(S8,_mm_mul_pd(d4,FA));

#define conv_10_block(filter,offset_source,d0,d1,d2,d3,d4) \
F = _mm_load_pd(filter);\
S0 = _mm_add_pd(S0,_mm_mul_pd(d0,F));\
S1 = _mm_add_pd(S1,_mm_mul_pd(d1,F));\
S2 = _mm_add_pd(S2,_mm_mul_pd(d2,F));\
S3 = _mm_add_pd(S3,_mm_mul_pd(d3,F));\
d4 = _mm_load_pd(source+offset_source);\
S4 = _mm_add_pd(S4,_mm_mul_pd(d4,F));

#define conv_6x2_block_fused(offset_filter,offset_source,d00,d10,d20,d30,d40,d50) \
FA = _mm_load_pd(filt+offset_filter);\
d00 = _mm_load_pd(source0+offset_source);\
S00 = _mm_add_pd(S00,_mm_mul_pd(d00,FA));\
FU = _mm_load_pd(filt_u+offset_filter);\
d10 = _mm_load_pd(source1+offset_source);\
S01 = _mm_add_pd(S01,_mm_mul_pd(d00,FU));\
S11 = _mm_add_pd(S11,_mm_mul_pd(d10,FU));\
d20 = _mm_load_pd(source2+offset_source);\
S10 = _mm_add_pd(S10,_mm_mul_pd(d10,FA));\
d30 = _mm_load_pd(source3+offset_source);\
S20 = _mm_add_pd(S20,_mm_mul_pd(d20,FA));\
S30 = _mm_add_pd(S30,_mm_mul_pd(d30,FA));\
S31 = _mm_add_pd(S31,_mm_mul_pd(d30,FU));\
d40 = _mm_load_pd(source4+offset_source);\
S21 = _mm_add_pd(S21,_mm_mul_pd(d20,FU));\
d50 = _mm_load_pd(source5+offset_source);\
S41 = _mm_add_pd(S41,_mm_mul_pd(d40,FU));\
S51 = _mm_add_pd(S51,_mm_mul_pd(d50,FU));\
S50 = _mm_add_pd(S50,_mm_mul_pd(d50,FA));\
S40 = _mm_add_pd(S40,_mm_mul_pd(d40,FA));

#define conv_12_block_fused(offset_filter,offset_source,d0,d1,d2,d3,d4,d5) \
FA = _mm_load_pd(filt+offset_filter);\
FU = _mm_load_pd(filt_u+offset_filter);\
S0 = _mm_add_pd(S0,_mm_mul_pd(d0,FA));\
S1 = _mm_add_pd(S1,_mm_mul_pd(d0,FU));\
S2 = _mm_add_pd(S2,_mm_mul_pd(d1,FA));\
S3 = _mm_add_pd(S3,_mm_mul_pd(d1,FU));\
d5 = _mm_load_pd(source+offset_source);\
S4 = _mm_add_pd(S4,_mm_mul_pd(d2,FA));\
S5 = _mm_add_pd(S5,_mm_mul_pd(d2,FU));\
S6 = _mm_add_pd(S6,_mm_mul_pd(d3,FA));\
S7 = _mm_add_pd(S7,_mm_mul_pd(d3,FU));\
S8 = _mm_add_pd(S8,_mm_mul_pd(d4,FA));\
S9 = _mm_add_pd(S9,_mm_mul_pd(d4,FU));\
S10 = _mm_add_pd(S10,_mm_mul_pd(d5,FA));\
S11 = _mm_add_pd(S11,_mm_mul_pd(d5,FU));

#define conv_12_block(filter,offset_source,d0,d1,d2,d3,d4,d5) \
F = _mm_load_pd(filter);\
S0 = _mm_add_pd(S0,_mm_mul_pd(d0,F));\
S1 = _mm_add_pd(S1,_mm_mul_pd(d1,F));\
S2 = _mm_add_pd(S2,_mm_mul_pd(d2,F));\
S3 = _mm_add_pd(S3,_mm_mul_pd(d3,F));\
S4 = _mm_add_pd(S4,_mm_mul_pd(d4,F));\
d5 = _mm_load_pd(source+offset_source);\
S5 = _mm_add_pd(S5,_mm_mul_pd(d5,F));

#define conv_14_block(filter,offset_source,d0,d1,d2,d3,d4,d5,d6) \
F = _mm_load_pd(filter);\
S0 = _mm_add_pd(S0,_mm_mul_pd(d0,F));\
S1 = _mm_add_pd(S1,_mm_mul_pd(d1,F));\
S2 = _mm_add_pd(S2,_mm_mul_pd(d2,F));\
S3 = _mm_add_pd(S3,_mm_mul_pd(d3,F));\
S4 = _mm_add_pd(S4,_mm_mul_pd(d4,F));\
S5 = _mm_add_pd(S5,_mm_mul_pd(d5,F));\
d6 = _mm_load_pd(source+offset_source);\
S6 = _mm_add_pd(S6,_mm_mul_pd(d6,F));

#define conv_16_block(filter,offset_source,d0,d1,d2,d3,d4,d5,d6,d7) \
F = _mm_load_pd(filter);\
S0 = _mm_add_pd(S0,_mm_mul_pd(d0,F));\
S1 = _mm_add_pd(S1,_mm_mul_pd(d1,F));\
S2 = _mm_add_pd(S2,_mm_mul_pd(d2,F));\
S3 = _mm_add_pd(S3,_mm_mul_pd(d3,F));\
S4 = _mm_add_pd(S4,_mm_mul_pd(d4,F));\
S5 = _mm_add_pd(S5,_mm_mul_pd(d5,F));\
S6 = _mm_add_pd(S6,_mm_mul_pd(d6,F));\
d7 = _mm_load_pd(source+offset_source);\
S7 = _mm_add_pd(S7,_mm_mul_pd(d7,F));

#define conv_18_block(filter,offset_source,d0,d1,d2,d3,d4,d5,d6,d7,d8) \
F = _mm_load_pd(filter);\
S0 = _mm_add_pd(S0,_mm_mul_pd(d0,F));\
S1 = _mm_add_pd(S1,_mm_mul_pd(d1,F));\
S2 = _mm_add_pd(S2,_mm_mul_pd(d2,F));\
S3 = _mm_add_pd(S3,_mm_mul_pd(d3,F));\
S4 = _mm_add_pd(S4,_mm_mul_pd(d4,F));\
S5 = _mm_add_pd(S5,_mm_mul_pd(d5,F));\
S6 = _mm_add_pd(S6,_mm_mul_pd(d6,F));\
S7 = _mm_add_pd(S7,_mm_mul_pd(d7,F));\
d8 = _mm_load_pd(source+offset_source);\
S8 = _mm_add_pd(S8,_mm_mul_pd(d8,F));


inline void conv_2x2_fused(double const * source0,double const * source1, double * dest0, double * dest1){
  __m128d S00,S01,S10,S11;
  __m128d FA,FU;
  __m128d D00,D10;
  FA = _mm_load_pd(filt);
  D00 = _mm_load_pd(source0);
  S00 = _mm_mul_pd(D00,FA);
  D10 = _mm_load_pd(source1);
  S10 = _mm_mul_pd(D10,FA);

  FU = _mm_set_pd(filt[0],filt[15]);
  S01 = _mm_set_pd(*(source0+1),*(source0+16));
  S01 = _mm_mul_pd(S01,FU);
  S11 = _mm_set_pd(*(source1+1),*(source1+16));
  S11 = _mm_mul_pd(S11,FU);

  conv_2x2_block_fused(2,2,D00,D10);  
  conv_2x2_block_fused(4,4,D00,D10);
  conv_2x2_block_fused(6,6,D00,D10);
  conv_2x2_block_fused(8,8,D00,D10);
  conv_2x2_block_fused(10,10,D00,D10);
  conv_2x2_block_fused(12,12,D00,D10);
  conv_2x2_block_fused(14,14,D00,D10);

  _mm_store_pd(dest0,_mm_hadd_pd(S00,S01));
  _mm_store_pd(dest1,_mm_hadd_pd(S10,S11));
}

inline void conv_2_fused(double const * source, double * dest){
  __m128d S0,S1;
  __m128d FA,FU;
  __m128d D0;
  FA = _mm_load_pd(filt);
  D0 = _mm_load_pd(source);
  S0 = _mm_mul_pd(D0,FA);

  FU = _mm_set_pd(filt[0],filt[15]);
  S1 = _mm_set_pd(*(source+1),*(source+16));
  S1 = _mm_mul_pd(S1,FU);

  conv_2_block_fused(2,2,D0);  
  conv_2_block_fused(4,4,D0);
  conv_2_block_fused(6,6,D0);
  conv_2_block_fused(8,8,D0);
  conv_2_block_fused(10,10,D0);
  conv_2_block_fused(12,12,D0);
  conv_2_block_fused(14,14,D0);

  _mm_store_pd(dest,_mm_hadd_pd(S0,S1));
}

inline void conv_2_even(double const * source, double * dest){
  __m128d S0;
  __m128d F;
  __m128d D0;
  F = _mm_load_pd(filt);
  D0 = _mm_load_pd(source);
  S0 = _mm_mul_pd(D0,F);

  conv_2_block(filt+2,2,D0);
  conv_2_block(filt+4,4,D0);
  conv_2_block(filt+6,6,D0);
  conv_2_block(filt+8,8,D0);
  conv_2_block(filt+10,10,D0);
  conv_2_block(filt+12,12,D0);
  conv_2_block(filt+14,14,D0);

  _mm_storel_pd(dest,_mm_hadd_pd(S0,S0));

} 

inline void conv_2_odd(double const * source, double * dest){
  __m128d S0;
  __m128d F;
  __m128d D0;
  F = _mm_set_pd(filt[0],filt[15]);
  S0 = _mm_set_pd(*(source+1),*(source+16));
  S0 = _mm_mul_pd(S0,F);

  F = _mm_load_pd(filt_u+2);
  D0 = _mm_load_pd(source+2);
  S0 = _mm_add_pd(S0,_mm_mul_pd(D0,F));

  conv_2_block(filt_u+4,4,D0);
  conv_2_block(filt_u+6,6,D0);
  conv_2_block(filt_u+8,8,D0);
  conv_2_block(filt_u+10,10,D0);
  conv_2_block(filt_u+12,12,D0);
  conv_2_block(filt_u+14,14,D0);

  _mm_storel_pd(dest+1,_mm_hadd_pd(S0,S0));
}

inline void conv_2x2_line_fused(size_t n,  double const * source0, double const * source1, double * dest0, double * dest1){
  unsigned int j=0;
  do {
    conv_2x2_fused(&source0[j],&source1[j],&dest0[j],&dest1[j]);
    j+=2;
  } while(j<n);
}

inline void conv_2_line_fused(size_t n,  double const * source, double * dest){
  unsigned int j=0;
  do {
    conv_2_fused(&source[j],&dest[j]);
    j+=2;
  } while(j<n);
}

inline void conv_2_line(size_t n,  double const * source, double * dest){
  unsigned int j=0;
  do {
    conv_2_even(&source[j],&dest[j]);
    conv_2_odd(&source[j],&dest[j]);
    j+=2;
  } while(j<n);
}

inline void conv_4x2_fused(size_t ndat, double const * source0, double const * source1, double const * source2, double const * source3,
                           double       * dest){
  __m128d S00,S01,S10,S11,S20,S21,S30,S31;
  __m128d FA,FU;
  __m128d D00,D10,D20,D30;

  FA = _mm_load_pd(filt);
  D00 = _mm_load_pd(source0);
  S00 = _mm_mul_pd(D00,FA);
  D10 = _mm_load_pd(source1);
  S10 = _mm_mul_pd(D10,FA);
  D20 = _mm_load_pd(source2);
  S20 = _mm_mul_pd(D20,FA);
  D30 = _mm_load_pd(source3);
  S30 = _mm_mul_pd(D30,FA);
  
  FU = _mm_load_pd(filt_u);
  S01 = _mm_loadl_pd(D00,source0+16);
  S01 = _mm_mul_pd(S01,FU);
  S11 = _mm_loadl_pd(D10,source1+16);
  S11 = _mm_mul_pd(S11,FU);
  S21 = _mm_loadl_pd(D20,source2+16);
  S21 = _mm_mul_pd(S21,FU);
  S31 = _mm_loadl_pd(D30,source3+16);
  S31 = _mm_mul_pd(S31,FU);

  conv_4x2_block_fused(2,2,D00,D10,D20,D30);
  conv_4x2_block_fused(4,4,D00,D10,D20,D30);
  conv_4x2_block_fused(6,6,D00,D10,D20,D30);
  conv_4x2_block_fused(8,8,D00,D10,D20,D30);
  conv_4x2_block_fused(10,10,D00,D10,D20,D30);
  conv_4x2_block_fused(12,12,D00,D10,D20,D30);
  conv_4x2_block_fused(14,14,D00,D10,D20,D30);

  _mm_store_pd(dest,_mm_hadd_pd(S00,S10));
  _mm_store_pd(dest+2,_mm_hadd_pd(S20,S30));
  _mm_store_pd(dest+ndat,_mm_hadd_pd(S01,S11));
  _mm_store_pd(dest+2+ndat,_mm_hadd_pd(S21,S31));
}

inline void conv_2x4_fused(size_t ndat, double const * source0, double const * source1, double * dest){
  __m128d S00,S01,S02,S03,S10,S11,S12,S13;
  __m128d FA,FU;
  __m128d D00,D01,D10,D11;

  FA = _mm_load_pd(filt);
  D00 = _mm_load_pd(source0);
  S00 = _mm_mul_pd(D00,FA);
  D01 = _mm_load_pd(source0+2);
  S02 = _mm_mul_pd(D01,FA);
  D10 = _mm_load_pd(source1);
  S10 = _mm_mul_pd(D10,FA);
  D11 = _mm_load_pd(source1+2);
  S12 = _mm_mul_pd(D11,FA);

  FU = _mm_load_pd(filt_u);
  S01 = _mm_loadl_pd(D00,source0+16);
  S01 = _mm_mul_pd(S01,FU);
  S03 = _mm_loadl_pd(D01,source0+18);
  S03 = _mm_mul_pd(S03,FU);
  S11 = _mm_loadl_pd(D10,source1+16);
  S11 = _mm_mul_pd(S11,FU);
  S13 = _mm_loadl_pd(D11,source1+18);
  S13 = _mm_mul_pd(S13,FU);

  conv_2x4_block_fused(2,4,D01,D00,D11,D10);
  conv_2x4_block_fused(4,6,D00,D01,D10,D11);
  conv_2x4_block_fused(6,8,D01,D00,D11,D10);
  conv_2x4_block_fused(8,10,D00,D01,D10,D11);
  conv_2x4_block_fused(10,12,D01,D00,D11,D10);
  conv_2x4_block_fused(12,14,D00,D01,D10,D11);
  conv_2x4_block_fused(14,16,D01,D00,D11,D10);

  _mm_store_pd(dest,_mm_hadd_pd(S00,S10));
  _mm_store_pd(dest+ndat,_mm_hadd_pd(S01,S11));
  _mm_store_pd(dest+2*ndat,_mm_hadd_pd(S02,S12));
  _mm_store_pd(dest+3*ndat,_mm_hadd_pd(S03,S13));
}

inline void conv_4_fused(double const * source, double * dest){
  __m128d S0,S1,S2,S3;
  __m128d FA,FU;
  __m128d D0,D1;

  FA = _mm_load_pd(filt);
  D0 = _mm_load_pd(source);
  S0 = _mm_mul_pd(D0,FA);
  D1 = _mm_load_pd(source+2);
  S2 = _mm_mul_pd(D1,FA);

  FU = _mm_set_pd(filt[0],filt[15]);
  S1 = _mm_set_pd(*(source+1),*(source+16));
  S1 = _mm_mul_pd(S1,FU);
  S3 = _mm_set_pd(*(source+3),*(source+18));
  S3 = _mm_mul_pd(S3,FU);

  conv_4_block_fused(2,4,D1,D0);
  conv_4_block_fused(4,6,D0,D1);
  conv_4_block_fused(6,8,D1,D0);
  conv_4_block_fused(8,10,D0,D1);
  conv_4_block_fused(10,12,D1,D0);
  conv_4_block_fused(12,14,D0,D1);
  conv_4_block_fused(14,16,D1,D0);

  _mm_store_pd(dest,_mm_hadd_pd(S0,S1));
  _mm_store_pd(dest+2,_mm_hadd_pd(S2,S3));
}

inline void conv_4_even(double const * source, double * dest){
  __m128d S0,S1;
  __m128d F;
  __m128d D0,D1;
  F = _mm_load_pd(filt);
  D0 = _mm_load_pd(source);
  S0 = _mm_mul_pd(D0,F);
  D1 = _mm_load_pd(source+2);
  S1 = _mm_mul_pd(D1,F);

  conv_4_block(filt+2,4,D1,D0);
  conv_4_block(filt+4,6,D0,D1);
  conv_4_block(filt+6,8,D1,D0);
  conv_4_block(filt+8,10,D0,D1);
  conv_4_block(filt+10,12,D1,D0);
  conv_4_block(filt+12,14,D0,D1);
  conv_4_block(filt+14,16,D1,D0);

  _mm_storel_pd(dest,_mm_hadd_pd(S0,S0));
  _mm_storel_pd(dest+2,_mm_hadd_pd(S1,S1));

} 

inline void conv_4_odd(double const * source, double * dest){
  __m128d S0,S1;
  __m128d F;
  __m128d D0,D1;
  F = _mm_set_pd(filt[0],filt[15]);
  S0 = _mm_set_pd(*(source+1),*(source+16));
  S0 = _mm_mul_pd(S0,F);
  S1 = _mm_set_pd(*(source+3),*(source+18));
  S1 = _mm_mul_pd(S1,F);

  F = _mm_load_pd(filt_u+2);
  D0 = _mm_load_pd(source+2);
  S0 = _mm_add_pd(S0,_mm_mul_pd(D0,F));
  D1 = _mm_load_pd(source+4);
  S1 = _mm_add_pd(S1,_mm_mul_pd(D1,F));

  conv_4_block(filt_u+4,6,D1,D0);
  conv_4_block(filt_u+6,8,D0,D1);
  conv_4_block(filt_u+8,10,D1,D0);
  conv_4_block(filt_u+10,12,D0,D1);
  conv_4_block(filt_u+12,14,D1,D0);
  conv_4_block(filt_u+14,16,D0,D1);

  _mm_storel_pd(dest+1,_mm_hadd_pd(S0,S0));
  _mm_storel_pd(dest+3,_mm_hadd_pd(S1,S1));
}

inline void conv_4_line(size_t n, double const * source, double * dest){
  unsigned int j=0;
  do {
    conv_4_even(&source[j],&dest[j]);
    conv_4_odd(&source[j],&dest[j]);
    j+=4;
  } while(j<n);
}

inline void conv_4_line_fused(size_t n, double const * source, double * dest){
  unsigned int j=0;
  do {
    conv_4_fused(&source[j],&dest[j]);
    j+=4;
  } while(j<n);
}

inline void conv_2x4_line_fused_t(size_t n,size_t ndat, double const * source, double * dest){
  double * dest_t;
  unsigned int i=0;
  do{
    dest_t=dest;
    unsigned int j=0;
    do {
      conv_2x4_fused(ndat,source,source+n+FILTER_SIZE,dest_t);
      source+=4;
      j+=4;
      dest_t+=4*ndat;
    } while(j<n);
    source+=FILTER_SIZE*2+n;
    i+=2;
    dest+=2;
  } while (i<ndat);
}

inline void conv_4x2_line_fused_t(size_t n,size_t ndat, double const * source, double * dest){
  double * dest_t;
  unsigned int i=0;
  do{
    dest_t=dest;
    unsigned int j=0;
    do {
      conv_4x2_fused(ndat,source,source+n+FILTER_SIZE,source+2*n+2*FILTER_SIZE,source+3*n+3*FILTER_SIZE,dest_t);
      source+=2;
      j+=2;
      dest_t+=2*ndat;
    } while(j<n);
    source+=FILTER_SIZE*4+3*n;
    i+=4;
    dest+=4;
  } while (i<ndat);
}


inline void copy_buffer_16x16(size_t ndat, double *dest, double *buff) {
  unsigned int i;
  for(i=0;i<16;i++,dest+=ndat,buff+=16){
    _mm_store_pd(dest,_mm_load_pd(buff));
    _mm_store_pd(dest+2,_mm_load_pd(buff+2));
    _mm_store_pd(dest+4,_mm_load_pd(buff+4));
    _mm_store_pd(dest+6,_mm_load_pd(buff+6));
    _mm_store_pd(dest+8,_mm_load_pd(buff+8));
    _mm_store_pd(dest+10,_mm_load_pd(buff+10));
    _mm_store_pd(dest+12,_mm_load_pd(buff+12));
    _mm_store_pd(dest+14,_mm_load_pd(buff+14));
  }
}

inline void copy_buffer_8x16(size_t ndat, double *dest, double *buff) {
  unsigned int i;
  for(i=0;i<8;i++,dest+=ndat,buff+=16){
    _mm_store_pd(dest,_mm_load_pd(buff));
    _mm_store_pd(dest+2,_mm_load_pd(buff+2));
    _mm_store_pd(dest+4,_mm_load_pd(buff+4));
    _mm_store_pd(dest+6,_mm_load_pd(buff+6));
    _mm_store_pd(dest+8,_mm_load_pd(buff+8));
    _mm_store_pd(dest+10,_mm_load_pd(buff+10));
    _mm_store_pd(dest+12,_mm_load_pd(buff+12));
    _mm_store_pd(dest+14,_mm_load_pd(buff+14));
  }
}

inline void copy_buffer_8x32(size_t ndat, double *dest, double *buff) {
  unsigned int i;
  for(i=0;i<8;i++,dest+=ndat,buff+=32){
    memcpy(dest,buff,32*sizeof(double));
  }
}

inline void copy_buffer_4x32(size_t ndat, double *dest, double *buff) {
  unsigned int i;
  for(i=0;i<4;i++,dest+=ndat,buff+=32){
    memcpy(dest,buff,32*sizeof(double));
  }
}

inline void copy_buffer_4x64(size_t ndat, double *dest, double *buff) {
  unsigned int i;
  for(i=0;i<4;i++,dest+=ndat,buff+=64){
    memcpy(dest,buff,64*sizeof(double));
  }
}

inline void copy_buffer_2x128(size_t ndat, double *dest, double *buff) {
  unsigned int i;
  for(i=0;i<2;i++,dest+=ndat,buff+=128){
    memcpy(dest,buff,128*sizeof(double));
  }
}

inline void copy_buffer_2x256(size_t ndat, double *dest, double *buff) {
  unsigned int i;
  for(i=0;i<2;i++,dest+=ndat,buff+=256){
    memcpy(dest,buff,256*sizeof(double));
  }
}

inline void fill_buffer_16x16_t(size_t n, double const *src, double *buff) {
  unsigned int i;
  __m128d D0,D1;
  for(i=0;i<16;i+=2,src+=2*n){
    D0 = _mm_load_pd(src);
    D1 = _mm_load_pd(src+n);
    _mm_store_pd(buff+i, _mm_unpacklo_pd( D0, D1));
    _mm_store_pd(buff+i+16*1, _mm_unpackhi_pd( D0, D1));
    D0 = _mm_load_pd(src+2);
    D1 = _mm_load_pd(src+2+n);
    _mm_store_pd(buff+i+16*2, _mm_unpacklo_pd( D0, D1));
    _mm_store_pd(buff+i+16*3, _mm_unpackhi_pd( D0, D1));
    D0 = _mm_load_pd(src+4);
    D1 = _mm_load_pd(src+4+n);
    _mm_store_pd(buff+i+16*4, _mm_unpacklo_pd( D0, D1));
    _mm_store_pd(buff+i+16*5, _mm_unpackhi_pd( D0, D1));
    D0 = _mm_load_pd(src+6);
    D1 = _mm_load_pd(src+6+n);
    _mm_store_pd(buff+i+16*6, _mm_unpacklo_pd( D0, D1));
    _mm_store_pd(buff+i+16*7, _mm_unpackhi_pd( D0, D1));
    D0 = _mm_load_pd(src+8);
    D1 = _mm_load_pd(src+8+n);
    _mm_store_pd(buff+i+16*8, _mm_unpacklo_pd( D0, D1));
    _mm_store_pd(buff+i+16*9, _mm_unpackhi_pd( D0, D1));
    D0 = _mm_load_pd(src+10);
    D1 = _mm_load_pd(src+10+n);
    _mm_store_pd(buff+i+16*10, _mm_unpacklo_pd( D0, D1));
    _mm_store_pd(buff+i+16*11, _mm_unpackhi_pd( D0, D1));
    D0 = _mm_load_pd(src+12);
    D1 = _mm_load_pd(src+12+n);
    _mm_store_pd(buff+i+16*12, _mm_unpacklo_pd( D0, D1));
    _mm_store_pd(buff+i+16*13, _mm_unpackhi_pd( D0, D1));
    D0 = _mm_load_pd(src+14);
    D1 = _mm_load_pd(src+14+n);
    _mm_store_pd(buff+i+16*14, _mm_unpacklo_pd( D0, D1));
    _mm_store_pd(buff+i+16*15, _mm_unpackhi_pd( D0, D1));
  }
}

inline void fill_buffer_8x16_t(size_t n, double const *src, double *buff) {
  unsigned int i;
  __m128d D0,D1;
  for(i=0;i<16;i+=2,src+=2*n){
    D0 = _mm_load_pd(src);
    D1 = _mm_load_pd(src+n);
    _mm_store_pd(buff+i, _mm_unpacklo_pd( D0, D1));
    _mm_store_pd(buff+i+16*1, _mm_unpackhi_pd( D0, D1));
    D0 = _mm_load_pd(src+2);
    D1 = _mm_load_pd(src+2+n);
    _mm_store_pd(buff+i+16*2, _mm_unpacklo_pd( D0, D1));
    _mm_store_pd(buff+i+16*3, _mm_unpackhi_pd( D0, D1));
    D0 = _mm_load_pd(src+4);
    D1 = _mm_load_pd(src+4+n);
    _mm_store_pd(buff+i+16*4, _mm_unpacklo_pd( D0, D1));
    _mm_store_pd(buff+i+16*5, _mm_unpackhi_pd( D0, D1));
    D0 = _mm_load_pd(src+6);
    D1 = _mm_load_pd(src+6+n);
    _mm_store_pd(buff+i+16*6, _mm_unpacklo_pd( D0, D1));
    _mm_store_pd(buff+i+16*7, _mm_unpackhi_pd( D0, D1));
  }
}

inline void fill_buffer_8x32_t(size_t n, double const *src, double *buff) {
  unsigned int i;
  __m128d D0,D1;
  for(i=0;i<32;i+=2,src+=2*n){
    D0 = _mm_load_pd(src);
    D1 = _mm_load_pd(src+n);
    _mm_store_pd(buff+i, _mm_unpacklo_pd( D0, D1));
    _mm_store_pd(buff+i+32*1, _mm_unpackhi_pd( D0, D1));
    D0 = _mm_load_pd(src+2);
    D1 = _mm_load_pd(src+2+n);
    _mm_store_pd(buff+i+32*2, _mm_unpacklo_pd( D0, D1));
    _mm_store_pd(buff+i+32*3, _mm_unpackhi_pd( D0, D1));
    D0 = _mm_load_pd(src+4);
    D1 = _mm_load_pd(src+4+n);
    _mm_store_pd(buff+i+32*4, _mm_unpacklo_pd( D0, D1));
    _mm_store_pd(buff+i+32*5, _mm_unpackhi_pd( D0, D1));
    D0 = _mm_load_pd(src+6);
    D1 = _mm_load_pd(src+6+n);
    _mm_store_pd(buff+i+32*6, _mm_unpacklo_pd( D0, D1));
    _mm_store_pd(buff+i+32*7, _mm_unpackhi_pd( D0, D1));
  }
}


inline void fill_buffer_4x64_t(size_t n, double const *src, double *buff) {
  unsigned int i;
  __m128d D0,D1;
  for(i=0;i<64;i+=2,src+=2*n){
    D0 = _mm_load_pd(src);
    D1 = _mm_load_pd(src+n);
    _mm_store_pd(buff+i, _mm_unpacklo_pd( D0, D1));
    _mm_store_pd(buff+i+64*1, _mm_unpackhi_pd( D0, D1));
    D0 = _mm_load_pd(src+2);
    D1 = _mm_load_pd(src+2+n);
    _mm_store_pd(buff+i+64*2, _mm_unpacklo_pd( D0, D1));
    _mm_store_pd(buff+i+64*3, _mm_unpackhi_pd( D0, D1));
  }
}

inline void fill_buffer_2x128_t(size_t n, double const *src, double *buff) {
  unsigned int i;
  __m128d D0,D1;
  for(i=0;i<128;i+=2,src+=2*n){
    D0 = _mm_load_pd(src);
    D1 = _mm_load_pd(src+n);
    _mm_store_pd(buff+i, _mm_unpacklo_pd( D0, D1));
    _mm_store_pd(buff+i+128*1, _mm_unpackhi_pd( D0, D1));
  }
}

inline void fill_buffer_2x256_t(size_t n, double const *src, double *buff) {
  unsigned int i;
  __m128d D0,D1;
  for(i=0;i<256;i+=2,src+=2*n){
    D0 = _mm_load_pd(src);
    D1 = _mm_load_pd(src+n);
    _mm_store_pd(buff+i, _mm_unpacklo_pd( D0, D1));
    _mm_store_pd(buff+i+256, _mm_unpackhi_pd( D0, D1));
  }
}

void print_first_block_16x16(size_t ndat, double const * source){
  int i,j;
  for(i=0;i<16;i++){
    for(j=0;j<16;j++){
      printf("%lf ",source[i*ndat+j]);
    }
    printf("\n");
  }
}

void zero_buffer_16x16(double *buff){
  int i,j;
  for(i=0;i<16;i++){
    for(j=0;j<16;j++){
      buff[i*16+j]=0.0;
    }
  }

}

void print_buffer_16x16(double * buff){
  int i,j;
  for(i=0;i<16;i++){
    for(j=0;j<16;j++){
      printf("%lf ",buff[i*16+j]);
    }
    printf("\n");
  }
}

inline void transpose_8x32(size_t n, size_t ndat, double const * source, double * dest){
  double buff[8*32] __attribute__ ((aligned (16)));
  double * dest_t;
  unsigned int i=0;
  do{
    dest_t=dest;
    unsigned int j=0;
    do{
      fill_buffer_8x32_t(n+FILTER_SIZE, source,buff);
      copy_buffer_8x32(ndat,dest_t,buff);
      source+=8;
      j+=8;
      dest_t+=8*ndat;
    } while(j<n);
    source+=FILTER_SIZE*32+31*n;
    i+=32;
    dest+=32;
  } while (i<ndat);
}

inline void transpose_4x64(size_t n, size_t ndat, double const * source, double * dest){
  double buff[4*64] __attribute__ ((aligned (16)));
  double * dest_t;
  unsigned int i=0;
  do{
    dest_t=dest;
    unsigned int j=0;
    do{
      fill_buffer_4x64_t(n+FILTER_SIZE, source,buff);
      copy_buffer_4x64(ndat,dest_t,buff);
      source+=4;
      j+=4;
      dest_t+=4*ndat;
    } while(j<n);
    source+=FILTER_SIZE*64+63*n;
    i+=64;
    dest+=64;
  } while (i<ndat);
}

inline void transpose_2x128(size_t n, size_t ndat, double const * source, double * dest){
  double buff[2*128] __attribute__ ((aligned (16)));
  double * dest_t;
  unsigned int i=0;
  do{
    dest_t=dest;
    unsigned int j=0;
    do{
      fill_buffer_2x128_t(n+FILTER_SIZE, source,buff);
      copy_buffer_2x128(ndat,dest_t,buff);
      source+=2;
      j+=2;
      dest_t+=2*ndat;
    } while(j<n);
    source+=FILTER_SIZE*128+127*n;
    i+=128;
    dest+=128;
  } while (i<ndat);
}

inline void transpose_2x256(size_t n, size_t ndat, double const * source, double * dest){
  double buff[2*256] __attribute__ ((aligned (16)));
  double * dest_t;
  unsigned int i=0;
  do{
    dest_t=dest;
    unsigned int j=0;
    do{
      fill_buffer_2x256_t(n+FILTER_SIZE, source,buff);
      copy_buffer_2x256(ndat,dest_t,buff);
      source+=2;
      j+=2;
      dest_t+=2*ndat;
    } while(j<n);
    source+=FILTER_SIZE*256+255*n;
    i+=256;
    dest+=256;
  } while (i<ndat);
}

inline void transpose_8x16(size_t n, size_t ndat, double const * source, double * dest){
  double buff[8*16] __attribute__ ((aligned (16)));
  double * dest_t;
  unsigned int i=0;
  do{
    dest_t=dest;
    unsigned int j=0;
    do{
      fill_buffer_8x16_t(n+FILTER_SIZE, source,buff);
      copy_buffer_8x16(ndat,dest_t,buff);
      source+=8;
      j+=8;
      dest_t+=8*ndat;
    } while(j<n);
    source+=FILTER_SIZE*16+15*n;
    i+=16;
    dest+=16;
  } while (i<ndat);
}

inline void copy(size_t n, size_t ndat, double const * source, double * dest){
  unsigned int i=0;
  do{
    memcpy(dest,source,n*sizeof(double));
    source+=n+FILTER_SIZE;
    dest+=n;
    i+=1;
  } while (i<ndat);
}

inline void copy_v2(size_t n, size_t ndat, double const * source, double * dest){
  unsigned int i=0;
  do{
    unsigned int j=0;
    do{
      _mm_store_pd(dest+j,_mm_load_pd(source+j));
      j+=2;
    } while(j<n);
    source+=n+FILTER_SIZE;
    dest+=n;
    i+=1;
  } while (i<ndat);
}

inline void copy_v16(size_t n, size_t ndat, double const * source, double * dest){
  unsigned int i=0;
  do{
    unsigned int j=0;
    do{
      _mm_store_pd(dest,_mm_load_pd(source));
      _mm_store_pd(dest+2,_mm_load_pd(source+2));
      _mm_store_pd(dest+4,_mm_load_pd(source+4));
      _mm_store_pd(dest+6,_mm_load_pd(source+6));
      _mm_store_pd(dest+8,_mm_load_pd(source+8));
      _mm_store_pd(dest+10,_mm_load_pd(source+10));
      _mm_store_pd(dest+12,_mm_load_pd(source+12));
      _mm_store_pd(dest+14,_mm_load_pd(source+14));
      dest+=16;
      source+=16;
      j+=16;
    } while(j<n);
    source+=FILTER_SIZE;
    i+=1;
  } while (i<ndat);
}


inline void transpose_16x16(size_t n, size_t ndat, double const * source, double * dest){
  double buff[16*16] __attribute__ ((aligned (16)));
  double * dest_t;
  unsigned int i=0;
  do{
    dest_t=dest;
    unsigned int j=0;
    do{
      fill_buffer_16x16_t(n+FILTER_SIZE, source,buff);
      copy_buffer_16x16(ndat,dest_t,buff);
      source+=16;
      j+=16;
      dest_t+=16*ndat;
    } while(j<n);
    source+=FILTER_SIZE*16+15*n;
    i+=16;
    dest+=16;
  } while (i<ndat);
}

inline void transpose_ref(size_t n, size_t ndat, double const * source, double * dest){
  unsigned int i=0;
  do{
    unsigned int j=0;
    do{
      dest[j*ndat+i]=source[i*(n+FILTER_SIZE)+j];
      j+=1;
    } while(j<n);
    i+=1;
  } while (i<ndat);
}

inline void transpose_ref2(size_t n, size_t ndat, double const * source, double * dest){
  unsigned int j=0;
  do{
    unsigned int i=0;
    do{
      dest[j*ndat+i]=source[i*(n+FILTER_SIZE)+j];
      i+=1;
    } while (i<ndat);
    j+=1;
  } while(j<n);
}

inline void conv_4x2_line_fused_2x128_tb(size_t n,size_t ndat, double const * source, double * dest){
  double buff[2*128] __attribute__ ((aligned (16)));
  double * dest_t;
  double const * source_t0, *source_t1, *source_t2, *source_t3;
  unsigned int i=0;
  do{
    dest_t=dest;
    unsigned int j=0;
    do {
      unsigned int k=0;
      source_t0 = source;
      source_t1 = source_t0 + n+FILTER_SIZE;
      source_t2 = source_t1 + n+FILTER_SIZE;
      source_t3 = source_t2 + n+FILTER_SIZE;
      do {
        conv_4x2_fused(128,source_t0,source_t1,source_t2,source_t3,buff+k);
        source_t0 = source_t3 + n+FILTER_SIZE;
        source_t1 = source_t0 + n+FILTER_SIZE;
        source_t2 = source_t1 + n+FILTER_SIZE;
        source_t3 = source_t2 + n+FILTER_SIZE;
        k+=4;
      } while(k<128);


      copy_buffer_2x128(ndat,dest_t,buff);
      source+=2;
      j+=2;
      dest_t+=2*ndat;
    } while(j<n);
    source+=FILTER_SIZE*128+127*n;
    i+=128;
    dest+=128;
  } while (i<ndat);
}

inline void conv_4x2_line_fused_4x64_tb(size_t n,size_t ndat, double const * source, double * dest){
  double buff[4*64] __attribute__ ((aligned (16)));
  double * dest_t;
  double const * source_t0, *source_t1, *source_t2, *source_t3;
  unsigned int i=0;
  do{
    dest_t=dest;
    unsigned int j=0;
    do {
      unsigned int k=0;
      source_t0 = source;
      source_t1 = source_t0 + n+FILTER_SIZE;
      source_t2 = source_t1 + n+FILTER_SIZE;
      source_t3 = source_t2 + n+FILTER_SIZE;
      do {
        conv_4x2_fused(64,source_t0,source_t1,source_t2,source_t3,buff+k);
        conv_4x2_fused(64,source_t0+2,source_t1+2,source_t2+2,source_t3+2,buff+k+2*64);
        source_t0 = source_t3 + n+FILTER_SIZE;
        source_t1 = source_t0 + n+FILTER_SIZE;
        source_t2 = source_t1 + n+FILTER_SIZE;
        source_t3 = source_t2 + n+FILTER_SIZE;
        k+=4;
      } while(k<64);
      copy_buffer_4x64(ndat,dest_t,buff);
      source+=4;
      j+=4;
      dest_t+=4*ndat;
    } while(j<n);
    source+=FILTER_SIZE*64+63*n;
    i+=64;
    dest+=64;
  } while (i<ndat);
}

inline void conv_4x2_line_fused_8x32_tb(size_t n,size_t ndat, double const * source, double * dest){
  double buff[8*32] __attribute__ ((aligned (16)));
  double * dest_t;
  double const * source_t0, *source_t1, *source_t2, *source_t3;
  unsigned int i=0;
  do{
    dest_t=dest;
    unsigned int j=0;
    do {
      unsigned int k=0;
      source_t0 = source;
      source_t1 = source_t0 + n+FILTER_SIZE;
      source_t2 = source_t1 + n+FILTER_SIZE;
      source_t3 = source_t2 + n+FILTER_SIZE;
      do {
        conv_4x2_fused(32,source_t0,source_t1,source_t2,source_t3,buff+k);
        conv_4x2_fused(32,source_t0+2,source_t1+2,source_t2+2,source_t3+2,buff+k+2*32);
        conv_4x2_fused(32,source_t0+4,source_t1+4,source_t2+4,source_t3+4,buff+k+4*32);
        conv_4x2_fused(32,source_t0+6,source_t1+6,source_t2+6,source_t3+6,buff+k+6*32);
        source_t0 = source_t3 + n+FILTER_SIZE;
        source_t1 = source_t0 + n+FILTER_SIZE;
        source_t2 = source_t1 + n+FILTER_SIZE;
        source_t3 = source_t2 + n+FILTER_SIZE;
        k+=4;
      } while(k<32);
      copy_buffer_8x32(ndat,dest_t,buff);
      source+=8;
      j+=8;
      dest_t+=8*ndat;
    } while(j<n);
    source+=FILTER_SIZE*32+31*n;
    i+=32;
    dest+=32;
  } while (i<ndat);
}

inline void conv_4x2_line_fused_4x32_tb(size_t n,size_t ndat, double const * source, double * dest){
  double buff[4*32] __attribute__ ((aligned (16)));
  double * dest_t;
  double const * source_t0, *source_t1, *source_t2, *source_t3;
  unsigned int i=0;
  do{
    dest_t=dest;
    unsigned int j=0;
    do {
      unsigned int k=0;
      source_t0 = source;
      source_t1 = source_t0 + n+FILTER_SIZE;
      source_t2 = source_t1 + n+FILTER_SIZE;
      source_t3 = source_t2 + n+FILTER_SIZE;
      do {
        conv_4x2_fused(32,source_t0,source_t1,source_t2,source_t3,buff+k);
        conv_4x2_fused(32,source_t0+2,source_t1+2,source_t2+2,source_t3+2,buff+k+2*32);
        source_t0 = source_t3 + n+FILTER_SIZE;
        source_t1 = source_t0 + n+FILTER_SIZE;
        source_t2 = source_t1 + n+FILTER_SIZE;
        source_t3 = source_t2 + n+FILTER_SIZE;
        k+=4;
      } while(k<32);
      copy_buffer_4x32(ndat,dest_t,buff);
      source+=4;
      j+=4;
      dest_t+=4*ndat;
    } while(j<n);
    source+=FILTER_SIZE*32+31*n;
    i+=32;
    dest+=32;
  } while (i<ndat);
}


inline void conv_4x2_line_fused_16x16_tb(size_t n,size_t ndat, double const * source, double * dest){
  double buff[16*16] __attribute__ ((aligned (16)));
  double * dest_t;
  double const * source_t0, *source_t1, *source_t2, *source_t3;
  unsigned int i=0;
  do{
    dest_t=dest;
    unsigned int j=0;
    do {
      source_t0 = source;
      source_t1 = source_t0+n+FILTER_SIZE;
      source_t2 = source_t1+n+FILTER_SIZE;
      source_t3 = source_t2+n+FILTER_SIZE;
      conv_4x2_fused(16,source_t0,source_t1,source_t2,source_t3,buff);
      conv_4x2_fused(16,source_t0+2,source_t1+2,source_t2+2,source_t3+2,buff+2*16);
      conv_4x2_fused(16,source_t0+4,source_t1+4,source_t2+4,source_t3+4,buff+4*16);
      conv_4x2_fused(16,source_t0+6,source_t1+6,source_t2+6,source_t3+6,buff+6*16);
      conv_4x2_fused(16,source_t0+8,source_t1+8,source_t2+8,source_t3+8,buff+8*16);
      conv_4x2_fused(16,source_t0+10,source_t1+10,source_t2+10,source_t3+10,buff+10*16);
      conv_4x2_fused(16,source_t0+12,source_t1+12,source_t2+12,source_t3+12,buff+12*16);
      conv_4x2_fused(16,source_t0+14,source_t1+14,source_t2+14,source_t3+14,buff+14*16);
      source_t0 = source_t3 + n+FILTER_SIZE;
      source_t1 = source_t0 + n+FILTER_SIZE;
      source_t2 = source_t1 + n+FILTER_SIZE;
      source_t3 = source_t2 + n+FILTER_SIZE;
      conv_4x2_fused(16,source_t0,source_t1,source_t2,source_t3,buff+4);
      conv_4x2_fused(16,source_t0+2,source_t1+2,source_t2+2,source_t3+2,buff+4+2*16);
      conv_4x2_fused(16,source_t0+4,source_t1+4,source_t2+4,source_t3+4,buff+4+4*16);
      conv_4x2_fused(16,source_t0+6,source_t1+6,source_t2+6,source_t3+6,buff+4+6*16);
      conv_4x2_fused(16,source_t0+8,source_t1+8,source_t2+8,source_t3+8,buff+4+8*16);
      conv_4x2_fused(16,source_t0+10,source_t1+10,source_t2+10,source_t3+10,buff+4+10*16);
      conv_4x2_fused(16,source_t0+12,source_t1+12,source_t2+12,source_t3+12,buff+4+12*16);
      conv_4x2_fused(16,source_t0+14,source_t1+14,source_t2+14,source_t3+14,buff+4+14*16);
      source_t0 = source_t3 + n+FILTER_SIZE;
      source_t1 = source_t0 + n+FILTER_SIZE;
      source_t2 = source_t1 + n+FILTER_SIZE;
      source_t3 = source_t2 + n+FILTER_SIZE;
      conv_4x2_fused(16,source_t0,source_t1,source_t2,source_t3,buff+8);
      conv_4x2_fused(16,source_t0+2,source_t1+2,source_t2+2,source_t3+2,buff+8+2*16);
      conv_4x2_fused(16,source_t0+4,source_t1+4,source_t2+4,source_t3+4,buff+8+4*16);
      conv_4x2_fused(16,source_t0+6,source_t1+6,source_t2+6,source_t3+6,buff+8+6*16);
      conv_4x2_fused(16,source_t0+8,source_t1+8,source_t2+8,source_t3+8,buff+8+8*16);
      conv_4x2_fused(16,source_t0+10,source_t1+10,source_t2+10,source_t3+10,buff+8+10*16);
      conv_4x2_fused(16,source_t0+12,source_t1+12,source_t2+12,source_t3+12,buff+8+12*16);
      conv_4x2_fused(16,source_t0+14,source_t1+14,source_t2+14,source_t3+14,buff+8+14*16);
      source_t0 = source_t3 + n+FILTER_SIZE;
      source_t1 = source_t0 + n+FILTER_SIZE;
      source_t2 = source_t1 + n+FILTER_SIZE;
      source_t3 = source_t2 + n+FILTER_SIZE;
      conv_4x2_fused(16,source_t0,source_t1,source_t2,source_t3,buff+12);
      conv_4x2_fused(16,source_t0+2,source_t1+2,source_t2+2,source_t3+2,buff+12+2*16);
      conv_4x2_fused(16,source_t0+4,source_t1+4,source_t2+4,source_t3+4,buff+12+4*16);
      conv_4x2_fused(16,source_t0+6,source_t1+6,source_t2+6,source_t3+6,buff+12+6*16);
      conv_4x2_fused(16,source_t0+8,source_t1+8,source_t2+8,source_t3+8,buff+12+8*16);
      conv_4x2_fused(16,source_t0+10,source_t1+10,source_t2+10,source_t3+10,buff+12+10*16);
      conv_4x2_fused(16,source_t0+12,source_t1+12,source_t2+12,source_t3+12,buff+12+12*16);
      conv_4x2_fused(16,source_t0+14,source_t1+14,source_t2+14,source_t3+14,buff+12+14*16);
      copy_buffer_16x16(ndat,dest_t,buff);
      source+=16;
      j+=16;
      dest_t+=16*ndat;
    } while(j<n);
    source+=FILTER_SIZE*16+15*n;
    i+=16;
    dest+=16;
  } while (i<ndat);
}

inline void conv_6_fused(double const * source, double * dest){
  __m128d S0,S1,S2,S3,S4,S5;
  __m128d FA,FU;
  __m128d D0,D1,D2;

  FA = _mm_load_pd(filt);
  D0 = _mm_load_pd(source);
  S0 = _mm_mul_pd(D0,FA);
  D1 = _mm_load_pd(source+2);
  S2 = _mm_mul_pd(D1,FA);
  D2 = _mm_load_pd(source+4);
  S4 = _mm_mul_pd(D2,FA);

  FU = _mm_set_pd(filt[0],filt[15]);
  S1 = _mm_set_pd(*(source+1),*(source+16));
  S1 = _mm_mul_pd(S1,FU);
  S3 = _mm_set_pd(*(source+3),*(source+18));
  S3 = _mm_mul_pd(S3,FU);
  S5 = _mm_set_pd(*(source+5),*(source+20));
  S5 = _mm_mul_pd(S5,FU);

  conv_6_block_fused(2,6,D1,D2,D0);
  conv_6_block_fused(4,8,D2,D0,D1);
  conv_6_block_fused(6,10,D0,D1,D2);
  conv_6_block_fused(8,12,D1,D2,D0);
  conv_6_block_fused(10,14,D2,D0,D1);
  conv_6_block_fused(12,16,D0,D1,D2);
  conv_6_block_fused(14,18,D1,D2,D0);

  _mm_store_pd(dest,_mm_hadd_pd(S0,S1));
  _mm_store_pd(dest+2,_mm_hadd_pd(S2,S3));
  _mm_store_pd(dest+4,_mm_hadd_pd(S4,S5));
}

inline void conv_6_even(double const * source, double * dest){
  __m128d S0,S1,S2;
  __m128d F;
  __m128d D0,D1,D2;
  F = _mm_load_pd(filt);
  D0 = _mm_load_pd(source);
  S0 = _mm_mul_pd(D0,F);
  D1 = _mm_load_pd(source+2);
  S1 = _mm_mul_pd(D1,F);
  D2 = _mm_load_pd(source+4);
  S2 = _mm_mul_pd(D2,F);

  conv_6_block(filt+2,6,D1,D2,D0);
  conv_6_block(filt+4,8,D2,D0,D1);
  conv_6_block(filt+6,10,D0,D1,D2);
  conv_6_block(filt+8,12,D1,D2,D0);
  conv_6_block(filt+10,14,D2,D0,D1);
  conv_6_block(filt+12,16,D0,D1,D2);
  conv_6_block(filt+14,18,D1,D2,D0);

  _mm_storel_pd(dest,_mm_hadd_pd(S0,S0));
  _mm_storel_pd(dest+2,_mm_hadd_pd(S1,S1));
  _mm_storel_pd(dest+4,_mm_hadd_pd(S2,S2));

} 

inline void conv_6_odd(double const * source, double * dest){
  __m128d S0,S1,S2;
  __m128d F;
  __m128d D0,D1,D2;
  F = _mm_set_pd(filt[0],filt[15]);
  S0 = _mm_set_pd(*(source+1),*(source+16));
  S0 = _mm_mul_pd(S0,F);
  S1 = _mm_set_pd(*(source+3),*(source+18));
  S1 = _mm_mul_pd(S1,F);
  S2 = _mm_set_pd(*(source+5),*(source+20));
  S2 = _mm_mul_pd(S2,F);

  F = _mm_load_pd(filt_u+2);
  D0 = _mm_load_pd(source+2);
  S0 = _mm_add_pd(S0,_mm_mul_pd(D0,F));
  D1 = _mm_load_pd(source+4);
  S1 = _mm_add_pd(S1,_mm_mul_pd(D1,F));
  D2 = _mm_load_pd(source+6);
  S2 = _mm_add_pd(S2,_mm_mul_pd(D2,F));

  conv_6_block(filt_u+4,8,D1,D2,D0);
  conv_6_block(filt_u+6,10,D2,D0,D1);
  conv_6_block(filt_u+8,12,D0,D1,D2);
  conv_6_block(filt_u+10,14,D1,D2,D0);
  conv_6_block(filt_u+12,16,D2,D0,D1);
  conv_6_block(filt_u+14,18,D0,D1,D2);

  _mm_storel_pd(dest+1,_mm_hadd_pd(S0,S0));
  _mm_storel_pd(dest+3,_mm_hadd_pd(S1,S1));
  _mm_storel_pd(dest+5,_mm_hadd_pd(S2,S2));
}

inline void conv_6_line(size_t n,  double const * source, double * dest){
  unsigned int j=0;
  do {
    conv_6_even(&source[j],&dest[j]);
    conv_6_odd(&source[j],&dest[j]);
    j+=6;
  } while(j<n);
}

inline void conv_6_line_fused(size_t n,  double const * source, double * dest){
  unsigned int j=0;
  do {
    conv_6_fused(&source[j],&dest[j]);
    j+=6;
  } while(j<n);
}

inline void conv_8_even_init(size_t n, double const * source, double * dest){
  __m128d S0,S1,S2,S3;
  __m128d F;
  __m128d D0,D1,D2,D3;
  F = _mm_load_pd(filt);
  D0 = _mm_load_pd(source+n-8);
  S0 = _mm_mul_pd(D0,F);
  D1 = _mm_load_pd(source+n-6);
  S1 = _mm_mul_pd(D1,F);
  D2 = _mm_load_pd(source+n-4);
  S2 = _mm_mul_pd(D2,F);
  D3 = _mm_load_pd(source+n-2);
  S3 = _mm_mul_pd(D3,F);

  conv_8_block(filt+2,0,D1,D2,D3,D0);
  conv_8_block(filt+4,2,D2,D3,D0,D1);
  conv_8_block(filt+6,4,D3,D0,D1,D2);
  conv_8_block(filt+8,6,D0,D1,D2,D3);
  conv_8_block(filt+10,8,D1,D2,D3,D0);
  conv_8_block(filt+12,10,D2,D3,D0,D1);
  conv_8_block(filt+14,12,D3,D0,D1,D2);

  _mm_storel_pd(dest,_mm_hadd_pd(S0,S0));
  _mm_storel_pd(dest+2,_mm_hadd_pd(S1,S1));
  _mm_storel_pd(dest+4,_mm_hadd_pd(S2,S2));
  _mm_storel_pd(dest+6,_mm_hadd_pd(S3,S3));

}

inline void conv_8_fused_t(unsigned int ndat, double const * source, double * dest){
  __m128d S0,S1,S2,S3,S4,S5,S6,S7;
  __m128d FA,FU;
  __m128d D0,D1,D2,D3;

  FA = _mm_load_pd(filt);
  D0 = _mm_set_pd(*(source+ndat),*(source));
  S0 = _mm_mul_pd(D0,FA);
  D1 = _mm_set_pd(*(source+3*ndat),*(source+2*ndat));
  S2 = _mm_mul_pd(D1,FA);
  D2 = _mm_set_pd(*(source+5*ndat),*(source+4*ndat));
  S4 = _mm_mul_pd(D2,FA);
  D3 = _mm_set_pd(*(source+7*ndat),*(source+6*ndat));
  S6 = _mm_mul_pd(D3,FA);

  FU = _mm_load_pd(filt_u);
  S1 = _mm_loadl_pd(D0,source+16*ndat);
  S1 = _mm_mul_pd(S1,FU);
  S3 = _mm_loadl_pd(D1,source+18*ndat);
  S3 = _mm_mul_pd(S3,FU);
  S5 = _mm_loadl_pd(D2,source+20*ndat);
  S5 = _mm_mul_pd(S5,FU);
  S7 = _mm_loadl_pd(D3,source+22*ndat);
  S7 = _mm_mul_pd(S7,FU);

  conv_8_block_fused_t(2,8,D1,D2,D3,D0);
  conv_8_block_fused_t(4,10,D2,D3,D0,D1);
  conv_8_block_fused_t(6,12,D3,D0,D1,D2);
  conv_8_block_fused_t(8,14,D0,D1,D2,D3);
  conv_8_block_fused_t(10,16,D1,D2,D3,D0);
  conv_8_block_fused_t(12,18,D2,D3,D0,D1);
  conv_8_block_fused_t(14,20,D3,D0,D1,D2);

  _mm_store_pd(dest,_mm_hadd_pd(S0,S1));
  _mm_store_pd(dest+2,_mm_hadd_pd(S2,S3));
  _mm_store_pd(dest+4,_mm_hadd_pd(S4,S5));
  _mm_store_pd(dest+6,_mm_hadd_pd(S6,S7));
}

inline void conv_8_fused(double const * source, double * dest){
  __m128d S0,S1,S2,S3,S4,S5,S6,S7;
  __m128d FA,FU;
  __m128d D0,D1,D2,D3;

  FA = _mm_load_pd(filt);
  D0 = _mm_load_pd(source);
  S0 = _mm_mul_pd(D0,FA);
  D1 = _mm_load_pd(source+2);
  S2 = _mm_mul_pd(D1,FA);
  D2 = _mm_load_pd(source+4);
  S4 = _mm_mul_pd(D2,FA);
  D3 = _mm_load_pd(source+6);
  S6 = _mm_mul_pd(D3,FA);

  FU = _mm_load_pd(filt_u);
  S1 = _mm_loadl_pd(D0,source+16);
  S1 = _mm_mul_pd(S1,FU);
  S3 = _mm_loadl_pd(D1,source+18);
  S3 = _mm_mul_pd(S3,FU);
  S5 = _mm_loadl_pd(D2,source+20);
  S5 = _mm_mul_pd(S5,FU);
  S7 = _mm_loadl_pd(D3,source+22);
  S7 = _mm_mul_pd(S7,FU);

  conv_8_block_fused(2,8,D1,D2,D3,D0);
  conv_8_block_fused(4,10,D2,D3,D0,D1);
  conv_8_block_fused(6,12,D3,D0,D1,D2);
  conv_8_block_fused(8,14,D0,D1,D2,D3);
  conv_8_block_fused(10,16,D1,D2,D3,D0);
  conv_8_block_fused(12,18,D2,D3,D0,D1);
  conv_8_block_fused(14,20,D3,D0,D1,D2);

  _mm_store_pd(dest,_mm_hadd_pd(S0,S1));
  _mm_store_pd(dest+2,_mm_hadd_pd(S2,S3));
  _mm_store_pd(dest+4,_mm_hadd_pd(S4,S5));
  _mm_store_pd(dest+6,_mm_hadd_pd(S6,S7));
}

inline void conv_8_even(double const * source, double * dest){
  __m128d S0,S1,S2,S3;
  __m128d F;
  __m128d D0,D1,D2,D3;
  F = _mm_load_pd(filt);
  D0 = _mm_load_pd(source);
  S0 = _mm_mul_pd(D0,F);
  D1 = _mm_load_pd(source+2);
  S1 = _mm_mul_pd(D1,F);
  D2 = _mm_load_pd(source+4);
  S2 = _mm_mul_pd(D2,F);
  D3 = _mm_load_pd(source+6);
  S3 = _mm_mul_pd(D3,F);

  conv_8_block(filt+2,8,D1,D2,D3,D0);
  conv_8_block(filt+4,10,D2,D3,D0,D1);
  conv_8_block(filt+6,12,D3,D0,D1,D2);
  conv_8_block(filt+8,14,D0,D1,D2,D3);
  conv_8_block(filt+10,16,D1,D2,D3,D0);
  conv_8_block(filt+12,18,D2,D3,D0,D1);
  conv_8_block(filt+14,20,D3,D0,D1,D2);

  _mm_storel_pd(dest,_mm_hadd_pd(S0,S0));
  _mm_storel_pd(dest+2,_mm_hadd_pd(S1,S1));
  _mm_storel_pd(dest+4,_mm_hadd_pd(S2,S2));
  _mm_storel_pd(dest+6,_mm_hadd_pd(S3,S3));

} 

inline void conv_8_odd(double const * source, double * dest){
  __m128d S0,S1,S2,S3;
  __m128d F;
  __m128d D0,D1,D2,D3;
  F = _mm_set_pd(filt[0],filt[15]);
  S0 = _mm_set_pd(*(source+1),*(source+16));
  S0 = _mm_mul_pd(S0,F);
  S1 = _mm_set_pd(*(source+3),*(source+18));
  S1 = _mm_mul_pd(S1,F);
  S2 = _mm_set_pd(*(source+5),*(source+20));
  S2 = _mm_mul_pd(S2,F);
  S3 = _mm_set_pd(*(source+7),*(source+22));
  S3 = _mm_mul_pd(S3,F);

  F = _mm_load_pd(filt_u+2);
  D0 = _mm_load_pd(source+2);
  S0 = _mm_add_pd(S0,_mm_mul_pd(D0,F));
  D1 = _mm_load_pd(source+4);
  S1 = _mm_add_pd(S1,_mm_mul_pd(D1,F));
  D2 = _mm_load_pd(source+6);
  S2 = _mm_add_pd(S2,_mm_mul_pd(D2,F));
  D3 = _mm_load_pd(source+8);
  S3 = _mm_add_pd(S3,_mm_mul_pd(D3,F));

  conv_8_block(filt_u+4,10,D1,D2,D3,D0);
  conv_8_block(filt_u+6,12,D2,D3,D0,D1);
  conv_8_block(filt_u+8,14,D3,D0,D1,D2);
  conv_8_block(filt_u+10,16,D0,D1,D2,D3);
  conv_8_block(filt_u+12,18,D1,D2,D3,D0);
  conv_8_block(filt_u+14,20,D2,D3,D0,D1);

  _mm_storel_pd(dest+1,_mm_hadd_pd(S0,S0));
  _mm_storel_pd(dest+3,_mm_hadd_pd(S1,S1));
  _mm_storel_pd(dest+5,_mm_hadd_pd(S2,S2));
  _mm_storel_pd(dest+7,_mm_hadd_pd(S3,S3));
}

inline void conv_8_line_fused(size_t n,  double const * source, double * dest){
  unsigned int j=0;
  do {
    conv_8_fused(&source[j],&dest[j]);
    j+=8;
  } while(j<n);
}

inline void conv_8_fused_tr(size_t n, size_t ndat,  double const * source, double * dest){
  unsigned int i=0;
  do {
    double const * source_t=source+i;
    unsigned int j=0;
    do {
      conv_8_fused_t(ndat,source_t+j*ndat,&dest[j]);
      j+=8;
    } while(j<n);
    i += 1;
    dest += n; 
  } while(i<ndat);
}


inline void conv_8_line(size_t n,  double const * source, double * dest){
  unsigned int j=0;
  do {
    conv_8_even(&source[j],&dest[j]);
    conv_8_odd(&source[j],&dest[j]);
    j+=8;
  } while(j<n);
}

inline void conv_10_fused(double const * source, double * dest){
  __m128d S0,S1,S2,S3,S4,S5,S6,S7,S8,S9;
  __m128d FA,FU;
  __m128d D0,D1,D2,D3,D4;

  FA = _mm_load_pd(filt);
  D0 = _mm_load_pd(source);
  S0 = _mm_mul_pd(D0,FA);
  D1 = _mm_load_pd(source+2);
  S2 = _mm_mul_pd(D1,FA);
  D2 = _mm_load_pd(source+4);
  S4 = _mm_mul_pd(D2,FA);
  D3 = _mm_load_pd(source+6);
  S6 = _mm_mul_pd(D3,FA);
  D4 = _mm_load_pd(source+8);
  S8 = _mm_mul_pd(D4,FA);

  FU = _mm_load_pd(filt_u);
  S1 = _mm_loadl_pd(D0,source+16);
  S1 = _mm_mul_pd(S1,FU);
  S3 = _mm_loadl_pd(D1,source+18);
  S3 = _mm_mul_pd(S3,FU);
  S5 = _mm_loadl_pd(D2,source+20);
  S5 = _mm_mul_pd(S5,FU);
  S7 = _mm_loadl_pd(D3,source+22);
  S7 = _mm_mul_pd(S7,FU);
  S9 = _mm_loadl_pd(D4,source+24);
  S9 = _mm_mul_pd(S9,FU);

  conv_10_block_fused(2,10,D1,D2,D3,D4,D0);
  conv_10_block_fused(4,12,D2,D3,D4,D0,D1);
  conv_10_block_fused(6,14,D3,D4,D0,D1,D2);
  conv_10_block_fused(8,16,D4,D0,D1,D2,D3);
  conv_10_block_fused(10,18,D0,D1,D2,D3,D4);
  conv_10_block_fused(12,20,D1,D2,D3,D4,D0);
  conv_10_block_fused(14,22,D2,D3,D4,D0,D1);

  _mm_store_pd(dest,_mm_hadd_pd(S0,S1));
  _mm_store_pd(dest+2,_mm_hadd_pd(S2,S3));
  _mm_store_pd(dest+4,_mm_hadd_pd(S4,S5));
  _mm_store_pd(dest+6,_mm_hadd_pd(S6,S7));
  _mm_store_pd(dest+8,_mm_hadd_pd(S8,S9));
}

inline void conv_10_even(double const * source, double * dest){
  __m128d S0,S1,S2,S3,S4;
  __m128d F;
  __m128d D0,D1,D2,D3,D4;
  F = _mm_load_pd(filt);
  D0 = _mm_load_pd(source);
  S0 = _mm_mul_pd(D0,F);
  D1 = _mm_load_pd(source+2);
  S1 = _mm_mul_pd(D1,F);
  D2 = _mm_load_pd(source+4);
  S2 = _mm_mul_pd(D2,F);
  D3 = _mm_load_pd(source+6);
  S3 = _mm_mul_pd(D3,F);
  D4 = _mm_load_pd(source+8);
  S4 = _mm_mul_pd(D4,F);

  conv_10_block(filt+2,10,D1,D2,D3,D4,D0);
  conv_10_block(filt+4,12,D2,D3,D4,D0,D1);
  conv_10_block(filt+6,14,D3,D4,D0,D1,D2);
  conv_10_block(filt+8,16,D4,D0,D1,D2,D3);
  conv_10_block(filt+10,18,D0,D1,D2,D3,D4);
  conv_10_block(filt+12,20,D1,D2,D3,D4,D0);
  conv_10_block(filt+14,22,D2,D3,D4,D0,D1);

  _mm_storel_pd(dest,_mm_hadd_pd(S0,S0));
  _mm_storel_pd(dest+2,_mm_hadd_pd(S1,S1));
  _mm_storel_pd(dest+4,_mm_hadd_pd(S2,S2));
  _mm_storel_pd(dest+6,_mm_hadd_pd(S3,S3));
  _mm_storel_pd(dest+8,_mm_hadd_pd(S4,S4));

} 

inline void conv_10_odd(double const * source, double * dest){
  __m128d S0,S1,S2,S3,S4;
  __m128d F;
  __m128d D0,D1,D2,D3,D4;
  F = _mm_set_pd(filt[0],filt[15]);
  S0 = _mm_set_pd(*(source+1),*(source+16));
  S0 = _mm_mul_pd(S0,F);
  S1 = _mm_set_pd(*(source+3),*(source+18));
  S1 = _mm_mul_pd(S1,F);
  S2 = _mm_set_pd(*(source+5),*(source+20));
  S2 = _mm_mul_pd(S2,F);
  S3 = _mm_set_pd(*(source+7),*(source+22));
  S3 = _mm_mul_pd(S3,F);
  S4 = _mm_set_pd(*(source+9),*(source+24));
  S4 = _mm_mul_pd(S4,F);

  F = _mm_load_pd(filt_u+2);
  D0 = _mm_load_pd(source+2);
  S0 = _mm_add_pd(S0,_mm_mul_pd(D0,F));
  D1 = _mm_load_pd(source+4);
  S1 = _mm_add_pd(S1,_mm_mul_pd(D1,F));
  D2 = _mm_load_pd(source+6);
  S2 = _mm_add_pd(S2,_mm_mul_pd(D2,F));
  D3 = _mm_load_pd(source+8);
  S3 = _mm_add_pd(S3,_mm_mul_pd(D3,F));
  D4 = _mm_load_pd(source+10);
  S4 = _mm_add_pd(S4,_mm_mul_pd(D4,F));

  conv_10_block(filt_u+4,12,D1,D2,D3,D4,D0);
  conv_10_block(filt_u+6,14,D2,D3,D4,D0,D1);
  conv_10_block(filt_u+8,16,D3,D4,D0,D1,D2);
  conv_10_block(filt_u+10,18,D4,D0,D1,D2,D3);
  conv_10_block(filt_u+12,20,D0,D1,D2,D3,D4);
  conv_10_block(filt_u+14,22,D1,D2,D3,D4,D0);

  _mm_storel_pd(dest+1,_mm_hadd_pd(S0,S0));
  _mm_storel_pd(dest+3,_mm_hadd_pd(S1,S1));
  _mm_storel_pd(dest+5,_mm_hadd_pd(S2,S2));
  _mm_storel_pd(dest+7,_mm_hadd_pd(S3,S3));
  _mm_storel_pd(dest+9,_mm_hadd_pd(S4,S4));
}
 
inline void conv_10_line(size_t n,  double const * source, double * dest){
  unsigned int j=0;
  do {
    conv_10_even(&source[j],&dest[j]);
    conv_10_odd(&source[j],&dest[j]);
    j+=10;
  } while(j<n);
}

inline void conv_10_line_fused(size_t n,  double const * source, double * dest){
  unsigned int j=0;
  do {
    conv_10_fused(&source[j],&dest[j]);
    j+=10;
  } while(j<n);
}

inline void conv_12_fused(double const * source, double * dest){
  __m128d S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11;
  __m128d FA,FU;
  __m128d D0,D1,D2,D3,D4,D5;

  FA = _mm_load_pd(filt);
  D0 = _mm_load_pd(source);
  S0 = _mm_mul_pd(D0,FA);
  D1 = _mm_load_pd(source+2);
  S2 = _mm_mul_pd(D1,FA);
  D2 = _mm_load_pd(source+4);
  S4 = _mm_mul_pd(D2,FA);
  D3 = _mm_load_pd(source+6);
  S6 = _mm_mul_pd(D3,FA);
  D4 = _mm_load_pd(source+8);
  S8 = _mm_mul_pd(D4,FA);
  D5 = _mm_load_pd(source+10);
  S10 = _mm_mul_pd(D5,FA);

  FU = _mm_set_pd(filt[0],filt[15]);
  S1 = _mm_set_pd(*(source+1),*(source+16));
  S1 = _mm_mul_pd(S1,FU);
  S3 = _mm_set_pd(*(source+3),*(source+18));
  S3 = _mm_mul_pd(S3,FU);
  S5 = _mm_set_pd(*(source+5),*(source+20));
  S5 = _mm_mul_pd(S5,FU);
  S7 = _mm_set_pd(*(source+7),*(source+22));
  S7 = _mm_mul_pd(S7,FU);
  S9 = _mm_set_pd(*(source+9),*(source+24));
  S9 = _mm_mul_pd(S9,FU);
  S11 = _mm_set_pd(*(source+11),*(source+26));
  S11 = _mm_mul_pd(S11,FU);

  conv_12_block_fused(2,12,D1,D2,D3,D4,D5,D0);
  conv_12_block_fused(4,14,D2,D3,D4,D5,D0,D1);
  conv_12_block_fused(6,16,D3,D4,D5,D0,D1,D2);
  conv_12_block_fused(8,18,D4,D5,D0,D1,D2,D3);
  conv_12_block_fused(10,20,D5,D0,D1,D2,D3,D4);
  conv_12_block_fused(12,22,D0,D1,D2,D3,D4,D5);
  conv_12_block_fused(14,24,D1,D2,D3,D4,D5,D0);

  _mm_store_pd(dest,_mm_hadd_pd(S0,S1));
  _mm_store_pd(dest+2,_mm_hadd_pd(S2,S3));
  _mm_store_pd(dest+4,_mm_hadd_pd(S4,S5));
  _mm_store_pd(dest+6,_mm_hadd_pd(S6,S7));
  _mm_store_pd(dest+8,_mm_hadd_pd(S8,S9));
  _mm_store_pd(dest+10,_mm_hadd_pd(S10,S11));
}

inline void conv_6x2_fused(size_t ndat, double const * source0, double const * source1,
                                        double const * source2, double const * source3,
                                        double const * source4, double const * source5,
                                        double       * dest){
  __m128d S00,S01,S10,S11,S20,S21,S30,S31,S40,S41,S50,S51;
  __m128d FA,FU;
  __m128d D00,D10,D20,D30,D40,D50;

  FA = _mm_load_pd(filt);
  D00 = _mm_load_pd(source0);
  S00 = _mm_mul_pd(D00,FA);
  D10 = _mm_load_pd(source1);
  S10 = _mm_mul_pd(D10,FA);
  D20 = _mm_load_pd(source2);
  S20 = _mm_mul_pd(D20,FA);
  D30 = _mm_load_pd(source3);
  S30 = _mm_mul_pd(D30,FA);
  D40 = _mm_load_pd(source4);
  S40 = _mm_mul_pd(D40,FA);
  D50 = _mm_load_pd(source5);
  S50 = _mm_mul_pd(D50,FA);
  
  FU = _mm_load_pd(filt_u);
  S01 = _mm_loadl_pd(D00,source0+16);
  S01 = _mm_mul_pd(S01,FU);
  S11 = _mm_loadl_pd(D10,source1+16);
  S11 = _mm_mul_pd(S11,FU);
  S21 = _mm_loadl_pd(D20,source2+16);
  S21 = _mm_mul_pd(S21,FU);
  S31 = _mm_loadl_pd(D30,source3+16);
  S31 = _mm_mul_pd(S31,FU);
  S41 = _mm_loadl_pd(D40,source4+16);
  S41 = _mm_mul_pd(S41,FU);
  S51 = _mm_loadl_pd(D50,source5+16);
  S51 = _mm_mul_pd(S51,FU);

  conv_6x2_block_fused(2,2,D00,D10,D20,D30,D40,D50);
  conv_6x2_block_fused(4,4,D00,D10,D20,D30,D40,D50);
  conv_6x2_block_fused(6,6,D00,D10,D20,D30,D40,D50);
  conv_6x2_block_fused(8,8,D00,D10,D20,D30,D40,D50);
  conv_6x2_block_fused(10,10,D00,D10,D20,D30,D40,D50);
  conv_6x2_block_fused(12,12,D00,D10,D20,D30,D40,D50);
  conv_6x2_block_fused(14,14,D00,D10,D20,D30,D40,D50);

  _mm_store_pd(dest,_mm_hadd_pd(S00,S10));
  _mm_store_pd(dest+2,_mm_hadd_pd(S20,S30));
  _mm_store_pd(dest+4,_mm_hadd_pd(S40,S50));
  _mm_store_pd(dest+ndat,_mm_hadd_pd(S01,S11));
  _mm_store_pd(dest+2+ndat,_mm_hadd_pd(S21,S31));
  _mm_store_pd(dest+4+ndat,_mm_hadd_pd(S41,S51));
}

inline void conv_12_even(double const * source, double * dest){
  __m128d S0,S1,S2,S3,S4,S5;
  __m128d F;
  __m128d D0,D1,D2,D3,D4,D5;
  F = _mm_load_pd(filt);
  D0 = _mm_load_pd(source);
  S0 = _mm_mul_pd(D0,F);
  D1 = _mm_load_pd(source+2);
  S1 = _mm_mul_pd(D1,F);
  D2 = _mm_load_pd(source+4);
  S2 = _mm_mul_pd(D2,F);
  D3 = _mm_load_pd(source+6);
  S3 = _mm_mul_pd(D3,F);
  D4 = _mm_load_pd(source+8);
  S4 = _mm_mul_pd(D4,F);
  D5 = _mm_load_pd(source+10);
  S5 = _mm_mul_pd(D5,F);

  conv_12_block(filt+2,12,D1,D2,D3,D4,D5,D0);
  conv_12_block(filt+4,14,D2,D3,D4,D5,D0,D1);
  conv_12_block(filt+6,16,D3,D4,D5,D0,D1,D2);
  conv_12_block(filt+8,18,D4,D5,D0,D1,D2,D3);
  conv_12_block(filt+10,20,D5,D0,D1,D2,D3,D4);
  conv_12_block(filt+12,22,D0,D1,D2,D3,D4,D5);
  conv_12_block(filt+14,24,D1,D2,D3,D4,D5,D0);

  _mm_storel_pd(dest,_mm_hadd_pd(S0,S0));
  _mm_storel_pd(dest+2,_mm_hadd_pd(S1,S1));
  _mm_storel_pd(dest+4,_mm_hadd_pd(S2,S2));
  _mm_storel_pd(dest+6,_mm_hadd_pd(S3,S3));
  _mm_storel_pd(dest+8,_mm_hadd_pd(S4,S4));
  _mm_storel_pd(dest+10,_mm_hadd_pd(S5,S5));

} 

inline void conv_12_odd(double const * source, double * dest){
  __m128d S0,S1,S2,S3,S4,S5;
  __m128d F;
  __m128d D0,D1,D2,D3,D4,D5;
  F = _mm_set_pd(filt[0],filt[15]);
  S0 = _mm_set_pd(*(source+1),*(source+16));
  S0 = _mm_mul_pd(S0,F);
  S1 = _mm_set_pd(*(source+3),*(source+18));
  S1 = _mm_mul_pd(S1,F);
  S2 = _mm_set_pd(*(source+5),*(source+20));
  S2 = _mm_mul_pd(S2,F);
  S3 = _mm_set_pd(*(source+7),*(source+22));
  S3 = _mm_mul_pd(S3,F);
  S4 = _mm_set_pd(*(source+9),*(source+24));
  S4 = _mm_mul_pd(S4,F);
  S5 = _mm_set_pd(*(source+11),*(source+26));
  S5 = _mm_mul_pd(S5,F);

  F = _mm_load_pd(filt_u+2);
  D0 = _mm_load_pd(source+2);
  S0 = _mm_add_pd(S0,_mm_mul_pd(D0,F));
  D1 = _mm_load_pd(source+4);
  S1 = _mm_add_pd(S1,_mm_mul_pd(D1,F));
  D2 = _mm_load_pd(source+6);
  S2 = _mm_add_pd(S2,_mm_mul_pd(D2,F));
  D3 = _mm_load_pd(source+8);
  S3 = _mm_add_pd(S3,_mm_mul_pd(D3,F));
  D4 = _mm_load_pd(source+10);
  S4 = _mm_add_pd(S4,_mm_mul_pd(D4,F));
  D5 = _mm_load_pd(source+12);
  S5 = _mm_add_pd(S5,_mm_mul_pd(D5,F));

  conv_12_block(filt_u+4,14,D1,D2,D3,D4,D5,D0);
  conv_12_block(filt_u+6,16,D2,D3,D4,D5,D0,D1);
  conv_12_block(filt_u+8,18,D3,D4,D5,D0,D1,D2);
  conv_12_block(filt_u+10,20,D4,D5,D0,D1,D2,D3);
  conv_12_block(filt_u+12,22,D5,D0,D1,D2,D3,D4);
  conv_12_block(filt_u+14,24,D0,D1,D2,D3,D4,D5);

  _mm_storel_pd(dest+1,_mm_hadd_pd(S0,S0));
  _mm_storel_pd(dest+3,_mm_hadd_pd(S1,S1));
  _mm_storel_pd(dest+5,_mm_hadd_pd(S2,S2));
  _mm_storel_pd(dest+7,_mm_hadd_pd(S3,S3));
  _mm_storel_pd(dest+9,_mm_hadd_pd(S4,S4));
  _mm_storel_pd(dest+11,_mm_hadd_pd(S5,S5));
}
 
inline void conv_12_line(size_t n,  double const * source, double * dest){
  unsigned int j=0;
  do {
    conv_12_even(&source[j],&dest[j]);
    conv_12_odd(&source[j],&dest[j]);
    j+=12;
  } while(j<n);
}

inline void conv_12_line_fused(size_t n,  double const * source, double * dest){
  unsigned int j=0;
  do {
    conv_12_fused(&source[j],&dest[j]);
    j+=12;
  } while(j<n);
}

inline void conv_6x2_line_fused_t(size_t n,size_t ndat, double const * source, double * dest){
  double * dest_t;
  unsigned int i=0;
  do{
    dest_t=dest;
    unsigned int j=0;
    do {
      conv_6x2_fused(ndat,source,                  source+n+FILTER_SIZE,
                          source+2*n+2*FILTER_SIZE, source+3*n+3*FILTER_SIZE,
                          source+4*n+4*FILTER_SIZE, source+5*n+5*FILTER_SIZE,
                          dest_t);
      source+=2;
      j+=2;
      dest_t+=2*ndat;
    } while(j<n);
    source+=FILTER_SIZE*6+5*n;
    i+=6;
    dest+=6;
  } while (i<ndat);
}
inline void conv_14_even(double const * source, double * dest){
  __m128d S0,S1,S2,S3,S4,S5,S6;
  __m128d F;
  __m128d D0,D1,D2,D3,D4,D5,D6;
  F = _mm_load_pd(filt);
  D0 = _mm_load_pd(source);
  S0 = _mm_mul_pd(D0,F);
  D1 = _mm_load_pd(source+2);
  S1 = _mm_mul_pd(D1,F);
  D2 = _mm_load_pd(source+4);
  S2 = _mm_mul_pd(D2,F);
  D3 = _mm_load_pd(source+6);
  S3 = _mm_mul_pd(D3,F);
  D4 = _mm_load_pd(source+8);
  S4 = _mm_mul_pd(D4,F);
  D5 = _mm_load_pd(source+10);
  S5 = _mm_mul_pd(D5,F);
  D6 = _mm_load_pd(source+12);
  S6 = _mm_mul_pd(D6,F);

  conv_14_block(filt+2,14,D1,D2,D3,D4,D5,D6,D0);
  conv_14_block(filt+4,16,D2,D3,D4,D5,D6,D0,D1);
  conv_14_block(filt+6,18,D3,D4,D5,D6,D0,D1,D2);
  conv_14_block(filt+8,20,D4,D5,D6,D0,D1,D2,D3);
  conv_14_block(filt+10,22,D5,D6,D0,D1,D2,D3,D4);
  conv_14_block(filt+12,24,D6,D0,D1,D2,D3,D4,D5);
  conv_14_block(filt+14,26,D0,D1,D2,D3,D4,D5,D6);

  _mm_storel_pd(dest,_mm_hadd_pd(S0,S0));
  _mm_storel_pd(dest+2,_mm_hadd_pd(S1,S1));
  _mm_storel_pd(dest+4,_mm_hadd_pd(S2,S2));
  _mm_storel_pd(dest+6,_mm_hadd_pd(S3,S3));
  _mm_storel_pd(dest+8,_mm_hadd_pd(S4,S4));
  _mm_storel_pd(dest+10,_mm_hadd_pd(S5,S5));
  _mm_storel_pd(dest+12,_mm_hadd_pd(S6,S6));

} 

inline void conv_14_odd(double const * source, double * dest){
  __m128d S0,S1,S2,S3,S4,S5,S6;
  __m128d F;
  __m128d D0,D1,D2,D3,D4,D5,D6;
  F = _mm_set_pd(filt[0],filt[15]);
  S0 = _mm_set_pd(*(source+1),*(source+16));
  S0 = _mm_mul_pd(S0,F);
  S1 = _mm_set_pd(*(source+3),*(source+18));
  S1 = _mm_mul_pd(S1,F);
  S2 = _mm_set_pd(*(source+5),*(source+20));
  S2 = _mm_mul_pd(S2,F);
  S3 = _mm_set_pd(*(source+7),*(source+22));
  S3 = _mm_mul_pd(S3,F);
  S4 = _mm_set_pd(*(source+9),*(source+24));
  S4 = _mm_mul_pd(S4,F);
  S5 = _mm_set_pd(*(source+11),*(source+26));
  S5 = _mm_mul_pd(S5,F);
  S6 = _mm_set_pd(*(source+13),*(source+28));
  S6 = _mm_mul_pd(S6,F);

  F = _mm_load_pd(filt_u+2);
  D0 = _mm_load_pd(source+2);
  S0 = _mm_add_pd(S0,_mm_mul_pd(D0,F));
  D1 = _mm_load_pd(source+4);
  S1 = _mm_add_pd(S1,_mm_mul_pd(D1,F));
  D2 = _mm_load_pd(source+6);
  S2 = _mm_add_pd(S2,_mm_mul_pd(D2,F));
  D3 = _mm_load_pd(source+8);
  S3 = _mm_add_pd(S3,_mm_mul_pd(D3,F));
  D4 = _mm_load_pd(source+10);
  S4 = _mm_add_pd(S4,_mm_mul_pd(D4,F));
  D5 = _mm_load_pd(source+12);
  S5 = _mm_add_pd(S5,_mm_mul_pd(D5,F));
  D6 = _mm_load_pd(source+14);
  S6 = _mm_add_pd(S6,_mm_mul_pd(D6,F));

  conv_14_block(filt_u+4,16,D1,D2,D3,D4,D5,D6,D0);
  conv_14_block(filt_u+6,18,D2,D3,D4,D5,D6,D0,D1);
  conv_14_block(filt_u+8,20,D3,D4,D5,D6,D0,D1,D2);
  conv_14_block(filt_u+10,22,D4,D5,D6,D0,D1,D2,D3);
  conv_14_block(filt_u+12,24,D5,D6,D0,D1,D2,D3,D4);
  conv_14_block(filt_u+14,26,D6,D0,D1,D2,D3,D4,D5);

  _mm_storel_pd(dest+1,_mm_hadd_pd(S0,S0));
  _mm_storel_pd(dest+3,_mm_hadd_pd(S1,S1));
  _mm_storel_pd(dest+5,_mm_hadd_pd(S2,S2));
  _mm_storel_pd(dest+7,_mm_hadd_pd(S3,S3));
  _mm_storel_pd(dest+9,_mm_hadd_pd(S4,S4));
  _mm_storel_pd(dest+11,_mm_hadd_pd(S5,S5));
  _mm_storel_pd(dest+13,_mm_hadd_pd(S6,S6));
}
 
inline void conv_14_line(size_t n,  double const * source, double * dest){
  unsigned int j=0;
  do {
    conv_14_even(&source[j],&dest[j]);
    conv_14_odd(&source[j],&dest[j]);
    j+=14;
  } while(j<n);
}

inline void conv_16_even(double const * source, double * dest){
  __m128d S0,S1,S2,S3,S4,S5,S6,S7;
  __m128d F;
  __m128d D0,D1,D2,D3,D4,D5,D6,D7;
  F = _mm_load_pd(filt);
  D0 = _mm_load_pd(source);
  S0 = _mm_mul_pd(D0,F);
  D1 = _mm_load_pd(source+2);
  S1 = _mm_mul_pd(D1,F);
  D2 = _mm_load_pd(source+4);
  S2 = _mm_mul_pd(D2,F);
  D3 = _mm_load_pd(source+6);
  S3 = _mm_mul_pd(D3,F);
  D4 = _mm_load_pd(source+8);
  S4 = _mm_mul_pd(D4,F);
  D5 = _mm_load_pd(source+10);
  S5 = _mm_mul_pd(D5,F);
  D6 = _mm_load_pd(source+12);
  S6 = _mm_mul_pd(D6,F);
  D7 = _mm_load_pd(source+14);
  S7 = _mm_mul_pd(D7,F);

  conv_16_block(filt+2,16,D1,D2,D3,D4,D5,D6,D7,D0);
  conv_16_block(filt+4,18,D2,D3,D4,D5,D6,D7,D0,D1);
  conv_16_block(filt+6,20,D3,D4,D5,D6,D7,D0,D1,D2);
  conv_16_block(filt+8,22,D4,D5,D6,D7,D0,D1,D2,D3);
  conv_16_block(filt+10,24,D5,D6,D7,D0,D1,D2,D3,D4);
  conv_16_block(filt+12,26,D6,D7,D0,D1,D2,D3,D4,D5);
  conv_16_block(filt+14,28,D7,D0,D1,D2,D3,D4,D5,D6);

  _mm_storel_pd(dest,_mm_hadd_pd(S0,S0));
  _mm_storel_pd(dest+2,_mm_hadd_pd(S1,S1));
  _mm_storel_pd(dest+4,_mm_hadd_pd(S2,S2));
  _mm_storel_pd(dest+6,_mm_hadd_pd(S3,S3));
  _mm_storel_pd(dest+8,_mm_hadd_pd(S4,S4));
  _mm_storel_pd(dest+10,_mm_hadd_pd(S5,S5));
  _mm_storel_pd(dest+12,_mm_hadd_pd(S6,S6));
  _mm_storel_pd(dest+14,_mm_hadd_pd(S7,S7));

} 

inline void conv_16_odd(double const * source, double * dest){
  __m128d S0,S1,S2,S3,S4,S5,S6,S7;
  __m128d F;
  __m128d D0,D1,D2,D3,D4,D5,D6,D7;
  F = _mm_set_pd(filt[0],filt[15]);
  S0 = _mm_set_pd(*(source+1),*(source+16));
  S0 = _mm_mul_pd(S0,F);
  S1 = _mm_set_pd(*(source+3),*(source+18));
  S1 = _mm_mul_pd(S1,F);
  S2 = _mm_set_pd(*(source+5),*(source+20));
  S2 = _mm_mul_pd(S2,F);
  S3 = _mm_set_pd(*(source+7),*(source+22));
  S3 = _mm_mul_pd(S3,F);
  S4 = _mm_set_pd(*(source+9),*(source+24));
  S4 = _mm_mul_pd(S4,F);
  S5 = _mm_set_pd(*(source+11),*(source+26));
  S5 = _mm_mul_pd(S5,F);
  S6 = _mm_set_pd(*(source+13),*(source+28));
  S6 = _mm_mul_pd(S6,F);
  S7 = _mm_set_pd(*(source+15),*(source+30));
  S7 = _mm_mul_pd(S7,F);

  F = _mm_load_pd(filt_u+2);
  D0 = _mm_load_pd(source+2);
  S0 = _mm_add_pd(S0,_mm_mul_pd(D0,F));
  D1 = _mm_load_pd(source+4);
  S1 = _mm_add_pd(S1,_mm_mul_pd(D1,F));
  D2 = _mm_load_pd(source+6);
  S2 = _mm_add_pd(S2,_mm_mul_pd(D2,F));
  D3 = _mm_load_pd(source+8);
  S3 = _mm_add_pd(S3,_mm_mul_pd(D3,F));
  D4 = _mm_load_pd(source+10);
  S4 = _mm_add_pd(S4,_mm_mul_pd(D4,F));
  D5 = _mm_load_pd(source+12);
  S5 = _mm_add_pd(S5,_mm_mul_pd(D5,F));
  D6 = _mm_load_pd(source+14);
  S6 = _mm_add_pd(S6,_mm_mul_pd(D6,F));
  D7 = _mm_load_pd(source+16);
  S7 = _mm_add_pd(S7,_mm_mul_pd(D7,F));

  conv_16_block(filt_u+4,18,D1,D2,D3,D4,D5,D6,D7,D0);
  conv_16_block(filt_u+6,20,D2,D3,D4,D5,D6,D7,D0,D1);
  conv_16_block(filt_u+8,22,D3,D4,D5,D6,D7,D0,D1,D2);
  conv_16_block(filt_u+10,24,D4,D5,D6,D7,D0,D1,D2,D3);
  conv_16_block(filt_u+12,26,D5,D6,D7,D0,D1,D2,D3,D4);
  conv_16_block(filt_u+14,28,D6,D7,D0,D1,D2,D3,D4,D5);

  _mm_storel_pd(dest+1,_mm_hadd_pd(S0,S0));
  _mm_storel_pd(dest+3,_mm_hadd_pd(S1,S1));
  _mm_storel_pd(dest+5,_mm_hadd_pd(S2,S2));
  _mm_storel_pd(dest+7,_mm_hadd_pd(S3,S3));
  _mm_storel_pd(dest+9,_mm_hadd_pd(S4,S4));
  _mm_storel_pd(dest+11,_mm_hadd_pd(S5,S5));
  _mm_storel_pd(dest+13,_mm_hadd_pd(S6,S6));
  _mm_storel_pd(dest+15,_mm_hadd_pd(S7,S7));
}
 
inline void conv_16_line(size_t n,  double const * source, double * dest){
  unsigned int j=0;
  do {
    conv_16_even(&source[j],&dest[j]);
    conv_16_odd(&source[j],&dest[j]);
    j+=16;
  } while(j<n);
}

inline void conv_18_even(double const * source, double * dest){
  __m128d S0,S1,S2,S3,S4,S5,S6,S7,S8;
  __m128d F;
  __m128d D0,D1,D2,D3,D4,D5,D6,D7,D8;
  F = _mm_load_pd(filt);
  D0 = _mm_load_pd(source);
  S0 = _mm_mul_pd(D0,F);
  D1 = _mm_load_pd(source+2);
  S1 = _mm_mul_pd(D1,F);
  D2 = _mm_load_pd(source+4);
  S2 = _mm_mul_pd(D2,F);
  D3 = _mm_load_pd(source+6);
  S3 = _mm_mul_pd(D3,F);
  D4 = _mm_load_pd(source+8);
  S4 = _mm_mul_pd(D4,F);
  D5 = _mm_load_pd(source+10);
  S5 = _mm_mul_pd(D5,F);
  D6 = _mm_load_pd(source+12);
  S6 = _mm_mul_pd(D6,F);
  D7 = _mm_load_pd(source+14);
  S7 = _mm_mul_pd(D7,F);
  D8 = _mm_load_pd(source+16);
  S8 = _mm_mul_pd(D8,F);

  conv_18_block(filt+2,18,D1,D2,D3,D4,D5,D6,D7,D8,D0);
  conv_18_block(filt+4,20,D2,D3,D4,D5,D6,D7,D8,D0,D1);
  conv_18_block(filt+6,22,D3,D4,D5,D6,D7,D8,D0,D1,D2);
  conv_18_block(filt+8,24,D4,D5,D6,D7,D8,D0,D1,D2,D3);
  conv_18_block(filt+10,26,D5,D6,D7,D8,D0,D1,D2,D3,D4);
  conv_18_block(filt+12,28,D6,D7,D8,D0,D1,D2,D3,D4,D5);
  conv_18_block(filt+14,30,D7,D8,D0,D1,D2,D3,D4,D5,D6);

  _mm_storel_pd(dest,_mm_hadd_pd(S0,S0));
  _mm_storel_pd(dest+2,_mm_hadd_pd(S1,S1));
  _mm_storel_pd(dest+4,_mm_hadd_pd(S2,S2));
  _mm_storel_pd(dest+6,_mm_hadd_pd(S3,S3));
  _mm_storel_pd(dest+8,_mm_hadd_pd(S4,S4));
  _mm_storel_pd(dest+10,_mm_hadd_pd(S5,S5));
  _mm_storel_pd(dest+12,_mm_hadd_pd(S6,S6));
  _mm_storel_pd(dest+14,_mm_hadd_pd(S7,S7));
  _mm_storel_pd(dest+16,_mm_hadd_pd(S8,S8));

} 

inline void conv_18_odd(double const * source, double * dest){
  __m128d S0,S1,S2,S3,S4,S5,S6,S7,S8;
  __m128d F;
  __m128d D0,D1,D2,D3,D4,D5,D6,D7,D8;
  F = _mm_set_pd(filt[0],filt[15]);
  S0 = _mm_set_pd(*(source+1),*(source+16));
  S0 = _mm_mul_pd(S0,F);
  S1 = _mm_set_pd(*(source+3),*(source+18));
  S1 = _mm_mul_pd(S1,F);
  S2 = _mm_set_pd(*(source+5),*(source+20));
  S2 = _mm_mul_pd(S2,F);
  S3 = _mm_set_pd(*(source+7),*(source+22));
  S3 = _mm_mul_pd(S3,F);
  S4 = _mm_set_pd(*(source+9),*(source+24));
  S4 = _mm_mul_pd(S4,F);
  S5 = _mm_set_pd(*(source+11),*(source+26));
  S5 = _mm_mul_pd(S5,F);
  S6 = _mm_set_pd(*(source+13),*(source+28));
  S6 = _mm_mul_pd(S6,F);
  S7 = _mm_set_pd(*(source+15),*(source+30));
  S7 = _mm_mul_pd(S7,F);
  S8 = _mm_set_pd(*(source+17),*(source+32));
  S8 = _mm_mul_pd(S8,F);

  F = _mm_load_pd(filt_u+2);
  D0 = _mm_load_pd(source+2);
  S0 = _mm_add_pd(S0,_mm_mul_pd(D0,F));
  D1 = _mm_load_pd(source+4);
  S1 = _mm_add_pd(S1,_mm_mul_pd(D1,F));
  D2 = _mm_load_pd(source+6);
  S2 = _mm_add_pd(S2,_mm_mul_pd(D2,F));
  D3 = _mm_load_pd(source+8);
  S3 = _mm_add_pd(S3,_mm_mul_pd(D3,F));
  D4 = _mm_load_pd(source+10);
  S4 = _mm_add_pd(S4,_mm_mul_pd(D4,F));
  D5 = _mm_load_pd(source+12);
  S5 = _mm_add_pd(S5,_mm_mul_pd(D5,F));
  D6 = _mm_load_pd(source+14);
  S6 = _mm_add_pd(S6,_mm_mul_pd(D6,F));
  D7 = _mm_load_pd(source+16);
  S7 = _mm_add_pd(S7,_mm_mul_pd(D7,F));
  D8 = _mm_load_pd(source+18);
  S8 = _mm_add_pd(S8,_mm_mul_pd(D8,F));

  conv_18_block(filt_u+4,20,D1,D2,D3,D4,D5,D6,D7,D8,D0);
  conv_18_block(filt_u+6,22,D2,D3,D4,D5,D6,D7,D8,D0,D1);
  conv_18_block(filt_u+8,24,D3,D4,D5,D6,D7,D8,D0,D1,D2);
  conv_18_block(filt_u+10,26,D4,D5,D6,D7,D8,D0,D1,D2,D3);
  conv_18_block(filt_u+12,28,D5,D6,D7,D8,D0,D1,D2,D3,D4);
  conv_18_block(filt_u+14,30,D6,D7,D8,D0,D1,D2,D3,D4,D5);

  _mm_storel_pd(dest+1,_mm_hadd_pd(S0,S0));
  _mm_storel_pd(dest+3,_mm_hadd_pd(S1,S1));
  _mm_storel_pd(dest+5,_mm_hadd_pd(S2,S2));
  _mm_storel_pd(dest+7,_mm_hadd_pd(S3,S3));
  _mm_storel_pd(dest+9,_mm_hadd_pd(S4,S4));
  _mm_storel_pd(dest+11,_mm_hadd_pd(S5,S5));
  _mm_storel_pd(dest+13,_mm_hadd_pd(S6,S6));
  _mm_storel_pd(dest+15,_mm_hadd_pd(S7,S7));
  _mm_storel_pd(dest+17,_mm_hadd_pd(S8,S8));
}
 
inline void conv_18_line(size_t n,  double const * source, double * dest){
  unsigned int j=0;
  do {
    conv_18_even(&source[j],&dest[j]);
    conv_18_odd(&source[j],&dest[j]);
    j+=18;
  } while(j<n);
}


inline void conv_16_even_bis(double const * source, double * dest){
  __m128d S0,S1,S2,S3,S4,S5,S6,S7;
  __m128d F;
  __m128d D;
  S0 = _mm_setzero_pd();
  S1 = _mm_setzero_pd();
  S2 = _mm_setzero_pd();
  S3 = _mm_setzero_pd();
  S4 = _mm_setzero_pd();
  S5 = _mm_setzero_pd();
  S6 = _mm_setzero_pd();
  S7 = _mm_setzero_pd();
  unsigned int i = 0;
  do {
    F = _mm_load_pd(filt+i);
    D = _mm_load_pd(source+i);
    S0 = _mm_add_pd(S0,_mm_mul_pd(D,F));
    D = _mm_load_pd(source+2+i);
    S1 = _mm_add_pd(S1,_mm_mul_pd(D,F));
    D = _mm_load_pd(source+4+i);
    S2 = _mm_add_pd(S2,_mm_mul_pd(D,F));
    D = _mm_load_pd(source+6+i);
    S3 = _mm_add_pd(S3,_mm_mul_pd(D,F));
    D = _mm_load_pd(source+8+i);
    S4 = _mm_add_pd(S4,_mm_mul_pd(D,F));
    D = _mm_load_pd(source+10+i);
    S5 = _mm_add_pd(S5,_mm_mul_pd(D,F));
    D = _mm_load_pd(source+12+i);
    S6 = _mm_add_pd(S6,_mm_mul_pd(D,F));
    D = _mm_load_pd(source+14+i);
    S7 = _mm_add_pd(S7,_mm_mul_pd(D,F));
//    source+=2;
    i+=2;
  } while (i<16);
  _mm_storel_pd(dest,_mm_hadd_pd(S0,S0));
  _mm_storel_pd(dest+2,_mm_hadd_pd(S1,S1));
  _mm_storel_pd(dest+4,_mm_hadd_pd(S2,S2));
  _mm_storel_pd(dest+6,_mm_hadd_pd(S3,S3));
  _mm_storel_pd(dest+8,_mm_hadd_pd(S4,S4));
  _mm_storel_pd(dest+10,_mm_hadd_pd(S5,S5));
  _mm_storel_pd(dest+12,_mm_hadd_pd(S6,S6));
  _mm_storel_pd(dest+14,_mm_hadd_pd(S7,S7));
  
}

inline void conv_16_odd_bis(double const * source, double * dest){
  __m128d S0,S1,S2,S3,S4,S5,S6,S7;
  __m128d F;
  __m128d D;
  F = _mm_set_pd(filt[0],filt[15]);
  S0 = _mm_set_pd(*(source+1),*(source+16));
  S0 = _mm_mul_pd(S0,F);
  S1 = _mm_set_pd(*(source+3),*(source+18));
  S1 = _mm_mul_pd(S1,F);
  S2 = _mm_set_pd(*(source+5),*(source+20));
  S2 = _mm_mul_pd(S2,F);
  S3 = _mm_set_pd(*(source+7),*(source+22));
  S3 = _mm_mul_pd(S3,F);
  S4 = _mm_set_pd(*(source+9),*(source+24));
  S4 = _mm_mul_pd(S4,F);
  S5 = _mm_set_pd(*(source+11),*(source+26));
  S5 = _mm_mul_pd(S5,F);
  S6 = _mm_set_pd(*(source+13),*(source+28));
  S6 = _mm_mul_pd(S6,F);
  S7 = _mm_set_pd(*(source+15),*(source+30));
  S7 = _mm_mul_pd(S7,F);
  
  unsigned int i = 2;
  do {
    F = _mm_load_pd(filt_u+i);
    D = _mm_load_pd(source+i);
    S0 = _mm_add_pd(S0,_mm_mul_pd(D,F));
    D = _mm_load_pd(source+2+i);
    S1 = _mm_add_pd(S1,_mm_mul_pd(D,F));
    D = _mm_load_pd(source+4+i);
    S2 = _mm_add_pd(S2,_mm_mul_pd(D,F));
    D = _mm_load_pd(source+6+i);
    S3 = _mm_add_pd(S3,_mm_mul_pd(D,F));
    D = _mm_load_pd(source+8+i);
    S4 = _mm_add_pd(S4,_mm_mul_pd(D,F));
    D = _mm_load_pd(source+10+i);
    S5 = _mm_add_pd(S5,_mm_mul_pd(D,F));
    D = _mm_load_pd(source+12+i);
    S6 = _mm_add_pd(S6,_mm_mul_pd(D,F));
    D = _mm_load_pd(source+14+i);
    S7 = _mm_add_pd(S7,_mm_mul_pd(D,F));
//    source += 2;
    i += 2;
  } while (i<16);
  _mm_storel_pd(dest+1,_mm_hadd_pd(S0,S0));
  _mm_storel_pd(dest+3,_mm_hadd_pd(S1,S1));
  _mm_storel_pd(dest+5,_mm_hadd_pd(S2,S2));
  _mm_storel_pd(dest+7,_mm_hadd_pd(S3,S3));
  _mm_storel_pd(dest+9,_mm_hadd_pd(S4,S4));
  _mm_storel_pd(dest+11,_mm_hadd_pd(S5,S5));
  _mm_storel_pd(dest+13,_mm_hadd_pd(S6,S6));
  _mm_storel_pd(dest+15,_mm_hadd_pd(S7,S7));
  
}

inline void conv_16_line_bis(size_t n,  double const * source, double * dest){
  unsigned int j=0;
  do {
    conv_16_even_bis(&source[j],&dest[j]);
    conv_16_odd_bis(&source[j],&dest[j]);
    j+=16;
  } while(j<n);
}

void conv_ref(size_t n, size_t ndat, double const * source, double * dest){
  unsigned int i,j,k;
  double tmp;
  double const * source_t;
  for(k=0;k<ndat;k++) {
    source_t = source + k*(n+FILTER_SIZE);
    for(j=0;j<n;j++,source_t++) {
      tmp = 0.0;
      for(i=0;i<FILTER_SIZE;i++) {
        tmp += source_t[i]*filt[i];
      }
      dest[j*ndat+k] = tmp;
    }
  }
}

void conv_ref_reverse(size_t n, size_t ndat, double const * source, double * dest){
  unsigned int i,j,k;
  double tmp;
  double const * source_t;
  for(k=0;k<ndat;k++) {
    source_t = source + k;
    for(j=0;j<n;j++,source_t+=ndat) {
      tmp = 0.0;
      for(i=0;i<FILTER_SIZE;i++) {
        tmp += source_t[i*ndat]*filt[i];
      }
      dest[j+k*n] = tmp;
      
    }
  }
}
#define FLOP (BUFFER_WIDTH*BUFFER_DEPTH)*2*FILTER_SIZE
#define MOP (BUFFER_WIDTH*BUFFER_DEPTH)*2*8

#define EVENT_NUMBER 7
int read_print_stop_counters( long long int *counters, char **event_name, int event_number){
    int i;
    int error;

    error = PAPI_stop_counters(counters,event_number);
    for(i=0; i<event_number; i++)
      printf("%s : %lld\n",event_name[i],counters[i]);
    return error;
}


inline void check_conv(void (*func) (size_t, size_t,  double const *, double *), double const * source, double * dest, double * check, char const * description) {
    int event_number = EVENT_NUMBER;
    long long int counters[EVENT_NUMBER];
    int events[EVENT_NUMBER] = {PAPI_L1_TCM,PAPI_L2_TCM,PAPI_L3_TCM,
                                PAPI_L2_TCA,PAPI_L3_TCA,
                                PAPI_TOT_CYC,PAPI_TOT_INS};
    char * event_name[] = { "PAPI_L1_TCM","PAPI_L2_TCM","PAPI_L3_TCM",
                            "PAPI_L2_TCA","PAPI_L3_TCA",
                            "PAPI_TOT_CYC","PAPI_TOT_INS"};
/*    int events[EVENT_NUMBER] = {PAPI_RES_STL, PAPI_BR_MSP, PAPI_BR_PRC, PAPI_BR_CN,
                                  PAPI_TOT_INS, PAPI_TOT_CYC};
      char * event_name[] = { "PAPI_RES_STL", "PAPI_BR_MSP", "PAPI_BR_PRC", "PAPI_BR_CN", 
                              "PAPI_TOT_INS", "PAPI_TOT_CYC"};*/

    unsigned long long int t1,t2;
//    memset(dest, 0, BUFFER_DEPTH*BUFFER_WIDTH*sizeof(double));

    PAPI_start_counters(events,event_number);
    nanosec(&t1);
    (*func)(BUFFER_DEPTH,BUFFER_WIDTH,source,dest);
    nanosec(&t2);
    read_print_stop_counters(counters, event_name, event_number);


    double br;
    br=0;
    unsigned int i,j;
    for(i=0; i<BUFFER_DEPTH;i++) {
      for(j=0;j<BUFFER_WIDTH;j++) {
        br+=dest[i*BUFFER_WIDTH+j];
        if(fabs(dest[i*BUFFER_WIDTH+j]-check[i*BUFFER_WIDTH+j])>1e-10){
          printf("error %u %u: %1.15lf != %1.15lf (error %1.15lf)!\n",i,j , dest[i*BUFFER_WIDTH+j], check[i*BUFFER_WIDTH+j], fabs(dest[i*BUFFER_WIDTH+j]-check[i*BUFFER_WIDTH+j]));
        }
      }
    }
    printf("result %10s %le, duration %10llu ns, FLOP %d, GFLOPS %lf\n", description, br, t2-t1, FLOP, (double)FLOP/(float)(t2-t1));
}

inline void check_trans(void (*func) (size_t, size_t,  double const *, double *), double const * source, double * dest, double * check, char const * description) {
    unsigned long long int t1,t2;
    int event_number = EVENT_NUMBER;
    long long int counters[EVENT_NUMBER];
    int events[EVENT_NUMBER] = {PAPI_L1_TCM,PAPI_L2_TCM,PAPI_L3_TCM,
                                PAPI_L2_TCA,PAPI_L3_TCA,
                                PAPI_TOT_CYC,PAPI_TOT_INS};
    char * event_name[] = { "PAPI_L1_TCM","PAPI_L2_TCM","PAPI_L3_TCM",
                            "PAPI_L2_TCA","PAPI_L3_TCA",
                            "PAPI_TOT_CYC","PAPI_TOT_INS"};
/*    int events[EVENT_NUMBER] = {PAPI_RES_STL, PAPI_BR_MSP, PAPI_BR_PRC, PAPI_BR_CN,
                                PAPI_TOT_INS, PAPI_TOT_CYC};
    char * event_name[] = { "PAPI_RES_STL", "PAPI_BR_MSP", "PAPI_BR_PRC", "PAPI_BR_CN", 
                            "PAPI_TOT_INS", "PAPI_TOT_CYC"};*/

    PAPI_start_counters(events,event_number);
    nanosec(&t1);
    (*func)(BUFFER_DEPTH,BUFFER_WIDTH,source,dest);
    nanosec(&t2);
    read_print_stop_counters(counters, event_name, event_number);


    double br;
    br=0;
    unsigned int i,j;
    for(i=0; i<BUFFER_DEPTH;i++) {
      for(j=0;j<BUFFER_WIDTH;j++) {
        br+=dest[i*BUFFER_WIDTH+j];
        if(fabs(dest[i*BUFFER_WIDTH+j]-check[i*BUFFER_WIDTH+j])>1e-10){
          printf("error %u %u: %1.15lf != %1.15lf (error %1.15lf)!\n",i,j , dest[i*BUFFER_WIDTH+j], check[i*BUFFER_WIDTH+j], fabs(dest[i*BUFFER_WIDTH+j]-check[i*BUFFER_WIDTH+j]));
        }
      }
    }
    printf("result %10s %le, duration %10llu ns, B %d, GB/s %lf (conv limit : %lf)\n",description, br, t2-t1, MOP, (double)MOP/(float)(t2-t1),(double)FLOP/(float)(t2-t1));
}


int main(void) {

  void * t;
  double * a;
  int err = posix_memalign(&t,16,BUFFER_WIDTH*(BUFFER_DEPTH+FILTER_SIZE)*sizeof(double));
  if(err){
    printf("memalign error : %i\n",err);
    exit(0);
  }
  a=t;
  double * b;
  err = posix_memalign(&t,16,BUFFER_DEPTH*BUFFER_WIDTH*sizeof(double));
  if(err){
    printf("memalign error : %i\n",err);
    exit(0);
  }
  b=t;
  double * c;
  err = posix_memalign(&t,16,BUFFER_DEPTH*BUFFER_WIDTH*sizeof(double));
  if(err){
    printf("memalign error : %i\n",err);
    exit(0);
  }
  c=t;
//  printf("a : %p, b : %p, c : %p\n",a,b,c);
  unsigned int i,j;
  for(i=0; i<BUFFER_WIDTH;i++) {
    for(j=0;j<BUFFER_DEPTH+FILTER_SIZE;j++){
      a[i*(BUFFER_DEPTH+FILTER_SIZE)+j] = rand()/(double)RAND_MAX;
    }
  }
  for(i=0; i<BUFFER_DEPTH;i++) {
    for(j=0;j<BUFFER_WIDTH;j++){
      b[i*BUFFER_WIDTH+j] =0.0;
    }
  }
  conv_ref(BUFFER_DEPTH,BUFFER_WIDTH,a,c);
  unsigned long long int t1,t2;
  double br;

  check_conv(conv_2x4_line_fused_t,a,b,c,"2x4f");
  check_conv(conv_4x2_line_fused_t,a,b,c,"4x2f");
  check_conv(conv_4x2_line_fused_16x16_tb,a,b,c,"4x2fb16x16");
  check_conv(conv_4x2_line_fused_8x32_tb,a,b,c,"4x2fb8x32");
  check_conv(conv_4x2_line_fused_4x32_tb,a,b,c,"4x2fb4x32");
  check_conv(conv_4x2_line_fused_4x64_tb,a,b,c,"4x2fb4x64");
  check_conv(conv_4x2_line_fused_2x128_tb,a,b,c,"4x2fb2x128");
  check_conv(conv_gliding_t,a,b,c,"gliding_t");
  conv_ref_reverse(BUFFER_DEPTH,BUFFER_WIDTH,a,c);
  check_conv(conv_gliding_t_reverse,a,b,c,"gliding_r");
  check_conv(conv_8_fused_tr,a,b,c,"8fr");

  transpose_ref(BUFFER_DEPTH,BUFFER_WIDTH,a,c);

  check_trans(transpose_ref,a,b,c,"t");
  check_trans(transpose_ref2,a,b,c,"t(2)");
  check_trans(transpose_16x16,a,b,c,"t 16x16");
  check_trans(transpose_8x16,a,b,c,"t 8x16");
  check_trans(transpose_8x32,a,b,c,"t 8x32");
  check_trans(transpose_4x64,a,b,c,"t 4x64");
  check_trans(transpose_2x128,a,b,c,"t 2x128");
  check_trans(transpose_2x256,a,b,c,"t 2x256");


  copy(BUFFER_DEPTH,BUFFER_WIDTH,a,c);
  check_trans(copy,a,b,c,"c");
  check_trans(copy_v2,a,b,c,"cv2");
  check_trans(copy_v16,a,b,c,"cv16");

  nanosec(&t1);
  memcpy(b,a,BUFFER_DEPTH*BUFFER_WIDTH*sizeof(double));
  nanosec(&t2);
  br=0;
  for(i=0; i<BUFFER_DEPTH;i++) {
    for(j=0;j<BUFFER_WIDTH;j++) {
      br+=b[i*BUFFER_WIDTH+j];
    }
  }
  printf("result %10s %le, duration %10llu ns, B %d, GB/s %lf (conv limit : %lf)\n","memcpy", br, t2-t1, MOP, (double)MOP/(float)(t2-t1),(double)FLOP/(float)(t2-t1));

  return 0;
}
