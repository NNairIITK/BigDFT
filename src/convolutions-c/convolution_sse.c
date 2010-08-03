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
#define BUFFER_WIDTH 32
#define BUFFER_DEPTH 8*3*2*5*7*3

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

#define conv_8_block_fused(offset_filter,offset_source,d0,d1,d2,d3) \
FA = _mm_load_pd(filt+offset_filter);\
FU = _mm_load_pd(filt_u+offset_filter);\
S0 = _mm_add_pd(S0,_mm_mul_pd(d0,FA));\
S1 = _mm_add_pd(S1,_mm_mul_pd(d0,FU));\
S2 = _mm_add_pd(S2,_mm_mul_pd(d1,FA));\
S3 = _mm_add_pd(S3,_mm_mul_pd(d1,FU));\
d3 = _mm_load_pd(source+offset_source);\
S4 = _mm_add_pd(S4,_mm_mul_pd(d2,FA));\
S5 = _mm_add_pd(S5,_mm_mul_pd(d2,FU));\
S6 = _mm_add_pd(S6,_mm_mul_pd(d3,FA));\
S7 = _mm_add_pd(S7,_mm_mul_pd(d3,FU));

#define conv_8_block(filter,offset_source,d0,d1,d2,d3) \
F = _mm_load_pd(filter);\
S0 = _mm_add_pd(S0,_mm_mul_pd(d0,F));\
S1 = _mm_add_pd(S1,_mm_mul_pd(d1,F));\
S2 = _mm_add_pd(S2,_mm_mul_pd(d2,F));\
d3 = _mm_load_pd(source+offset_source);\
S3 = _mm_add_pd(S3,_mm_mul_pd(d3,F));

#define conv_10_block_fused(offset_filter,offset_source,d0,d1,d2,d3,d4) \
FA = _mm_load_pd(filt+offset_filter);\
FU = _mm_load_pd(filt_u+offset_filter);\
S0 = _mm_add_pd(S0,_mm_mul_pd(d0,FA));\
S1 = _mm_add_pd(S1,_mm_mul_pd(d0,FU));\
S2 = _mm_add_pd(S2,_mm_mul_pd(d1,FA));\
S3 = _mm_add_pd(S3,_mm_mul_pd(d1,FU));\
d4 = _mm_load_pd(source+offset_source);\
S4 = _mm_add_pd(S4,_mm_mul_pd(d2,FA));\
S5 = _mm_add_pd(S5,_mm_mul_pd(d2,FU));\
S6 = _mm_add_pd(S6,_mm_mul_pd(d3,FA));\
S7 = _mm_add_pd(S7,_mm_mul_pd(d3,FU));\
S8 = _mm_add_pd(S8,_mm_mul_pd(d4,FA));\
S9 = _mm_add_pd(S9,_mm_mul_pd(d4,FU));

#define conv_10_block(filter,offset_source,d0,d1,d2,d3,d4) \
F = _mm_load_pd(filter);\
S0 = _mm_add_pd(S0,_mm_mul_pd(d0,F));\
S1 = _mm_add_pd(S1,_mm_mul_pd(d1,F));\
S2 = _mm_add_pd(S2,_mm_mul_pd(d2,F));\
S3 = _mm_add_pd(S3,_mm_mul_pd(d3,F));\
d4 = _mm_load_pd(source+offset_source);\
S4 = _mm_add_pd(S4,_mm_mul_pd(d4,F));

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

inline void conv_4_line(size_t n,  double const * source, double * dest){
  unsigned int j=0;
  do {
    conv_4_even(&source[j],&dest[j]);
    conv_4_odd(&source[j],&dest[j]);
    j+=4;
  } while(j<n);
}

inline void conv_4_line_fused(size_t n,  double const * source, double * dest){
  unsigned int j=0;
  do {
    conv_4_fused(&source[j],&dest[j]);
    j+=4;
  } while(j<n);
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

void conv_ref(double const * source, double * dest){
  unsigned int i,j;
  double tmp;
  for(j=0;j<BUFFER_DEPTH;j++,source++) {
    tmp = 0.0;
    for(i=0;i<FILTER_SIZE;i++) {
     tmp += source[i]*filt[i];
    }
    dest[j] = tmp;
  }
}

#define FLOP (BUFFER_WIDTH*BUFFER_DEPTH)*2*FILTER_SIZE
int main(void) {
  double a[BUFFER_WIDTH][BUFFER_DEPTH+FILTER_SIZE] __attribute__ ((aligned (16)));
  double b[BUFFER_WIDTH][BUFFER_DEPTH] __attribute__ ((aligned (16)));
  double c[BUFFER_WIDTH][BUFFER_DEPTH] __attribute__ ((aligned (16)));
  unsigned int i,j,k;
  for(i=1; i<BUFFER_WIDTH;i++) {
    for(j=0;j<BUFFER_DEPTH+FILTER_SIZE;j++){
      a[i][j] = rand()/(double)RAND_MAX;
    }
  }
  for(i=0; i<BUFFER_WIDTH;i++) {
    for(j=0;j<BUFFER_DEPTH;j++){
      b[i][j] =0.0;
    }
  }
  for(i=0; i<BUFFER_WIDTH;i++) {
      conv_ref(a[i],c[i]);
  }
  unsigned long long int t1,t2;
  nanosec(&t1);
  for(i=0; i<BUFFER_WIDTH;i++) {
      conv_16_line_bis(BUFFER_DEPTH,a[i],b[i]);
  }
  nanosec(&t2);
  double br=0;
  for(i=0; i<BUFFER_WIDTH;i++) {
    for(j=0;j<BUFFER_DEPTH;j++) {
      br+=b[i][j];
      if(fabs(b[i][j]-c[i][j])>1e-10){
        printf("error %u %u: %1.15lf != %1.15lf (error %1.15lf)!\n",i,j , b[i][j], c[i][j], fabs(b[i][j]-c[i][j]));
     }
    }
  }
  printf("result  16b %lf, duration %llu ns, FLOP %d, GFLOPS %lf\n", br, t2-t1, FLOP, (double)FLOP/(float)(t2-t1));

  nanosec(&t1);
  for(i=0; i<BUFFER_WIDTH;i++) {
      conv_2_line(BUFFER_DEPTH,a[i],b[i]);
  }
  nanosec(&t2);
  br=0;
  for(i=0; i<BUFFER_WIDTH;i++) {
    for(j=0;j<BUFFER_DEPTH;j++) {
      br+=b[i][j];
      if(fabs(b[i][j]-c[i][j])>1e-10){
        printf("error %u %u: %1.15lf != %1.15lf (error %1.15lf)!\n",i,j , b[i][j], c[i][j], fabs(b[i][j]-c[i][j]));
     }
    }
  }
  printf("result  2eo %lf, duration %llu ns, FLOP %d, GFLOPS %lf\n", br, t2-t1, FLOP, (double)FLOP/(float)(t2-t1));

  nanosec(&t1);
  for(i=0; i<BUFFER_WIDTH;i++) {
      conv_2_line_fused(BUFFER_DEPTH,a[i],b[i]);
  }
  nanosec(&t2);
  br=0;
  for(i=0; i<BUFFER_WIDTH;i++) {
    for(j=0;j<BUFFER_DEPTH;j++) {
      br+=b[i][j];
      if(fabs(b[i][j]-c[i][j])>1e-10){
        printf("error %u %u: %1.15lf != %1.15lf (error %1.15lf)!\n",i,j , b[i][j], c[i][j], fabs(b[i][j]-c[i][j]));
     }
    }
  }
  printf("result   2f %lf, duration %llu ns, FLOP %d, GFLOPS %lf\n", br, t2-t1, FLOP, (double)FLOP/(float)(t2-t1));

  nanosec(&t1);
  for(i=0; i<BUFFER_WIDTH;i++) {
      conv_4_line(BUFFER_DEPTH,a[i],b[i]);
  }
  nanosec(&t2);
  br=0;
  for(i=0; i<BUFFER_WIDTH;i++) {
    for(j=0;j<BUFFER_DEPTH;j++) {
      br+=b[i][j];
      if(fabs(b[i][j]-c[i][j])>1e-10){
        printf("error %u %u: %1.15lf != %1.15lf (error %1.15lf)!\n",i,j , b[i][j], c[i][j], fabs(b[i][j]-c[i][j]));
     }
    }
  }
  printf("result  4eo %lf, duration %llu ns, FLOP %d, GFLOPS %lf\n", br, t2-t1, FLOP, (double)FLOP/(float)(t2-t1));

  nanosec(&t1);
  for(i=0; i<BUFFER_WIDTH;i++) {
      conv_4_line_fused(BUFFER_DEPTH,a[i],b[i]);
  }
  nanosec(&t2);
  br=0;
  for(i=0; i<BUFFER_WIDTH;i++) {
    for(j=0;j<BUFFER_DEPTH;j++) {
      br+=b[i][j];
      if(fabs(b[i][j]-c[i][j])>1e-10){
        printf("error %u %u: %1.15lf != %1.15lf (error %1.15lf)!\n",i,j , b[i][j], c[i][j], fabs(b[i][j]-c[i][j]));
     }
    }
  }
  printf("result   4f %lf, duration %llu ns, FLOP %d, GFLOPS %lf\n", br, t2-t1, FLOP, (double)FLOP/(float)(t2-t1));

  nanosec(&t1);
  for(i=0; i<BUFFER_WIDTH;i++) {
      conv_6_line(BUFFER_DEPTH,a[i],b[i]);
  }
  nanosec(&t2);
  br=0;
  for(i=0; i<BUFFER_WIDTH;i++) {
    for(j=0;j<BUFFER_DEPTH;j++) {
      br+=b[i][j];
      if(fabs(b[i][j]-c[i][j])>1e-10){
        printf("error %u %u: %1.15lf != %1.15lf (error %1.15lf)!\n",i,j , b[i][j], c[i][j], fabs(b[i][j]-c[i][j]));
     }
    }
  }
  printf("result  6eo %lf, duration %llu ns, FLOP %d, GFLOPS %lf\n", br, t2-t1, FLOP, (double)FLOP/(float)(t2-t1));

  nanosec(&t1);
  for(i=0; i<BUFFER_WIDTH;i++) {
      conv_6_line_fused(BUFFER_DEPTH,a[i],b[i]);
  }
  nanosec(&t2);
  br=0;
  for(i=0; i<BUFFER_WIDTH;i++) {
    for(j=0;j<BUFFER_DEPTH;j++) {
      br+=b[i][j];
      if(fabs(b[i][j]-c[i][j])>1e-10){
        printf("error %u %u: %1.15lf != %1.15lf (error %1.15lf)!\n",i,j , b[i][j], c[i][j], fabs(b[i][j]-c[i][j]));
     }
    }
  }
  printf("result   6f %lf, duration %llu ns, FLOP %d, GFLOPS %lf\n", br, t2-t1, FLOP, (double)FLOP/(float)(t2-t1));

  nanosec(&t1);
  for(i=0; i<BUFFER_WIDTH;i++) {
      conv_8_line(BUFFER_DEPTH,a[i],b[i]);
  }
  nanosec(&t2);
  br=0;
  for(i=0; i<BUFFER_WIDTH;i++) {
    for(j=0;j<BUFFER_DEPTH;j++) {
      br+=b[i][j];
      if(fabs(b[i][j]-c[i][j])>1e-10){
        printf("error %u %u: %1.15lf != %1.15lf (error %1.15lf)!\n",i,j , b[i][j], c[i][j], fabs(b[i][j]-c[i][j]));
     }
    }
  }
  printf("result  8eo %lf, duration %llu ns, FLOP %d, GFLOPS %lf\n", br, t2-t1, FLOP, (double)FLOP/(float)(t2-t1));

  nanosec(&t1);
  for(i=0; i<BUFFER_WIDTH;i++) {
      conv_8_line_fused(BUFFER_DEPTH,a[i],b[i]);
  }
  nanosec(&t2);
  br=0;
  for(i=0; i<BUFFER_WIDTH;i++) {
    for(j=0;j<BUFFER_DEPTH;j++) {
      br+=b[i][j];
      if(fabs(b[i][j]-c[i][j])>1e-10){
        printf("error %u %u: %1.15lf != %1.15lf (error %1.15lf)!\n",i,j , b[i][j], c[i][j], fabs(b[i][j]-c[i][j]));
     }
    }
  }
  printf("result   8f %lf, duration %llu ns, FLOP %d, GFLOPS %lf\n", br, t2-t1, FLOP, (double)FLOP/(float)(t2-t1));

  nanosec(&t1);
  for(i=0; i<BUFFER_WIDTH;i++) {
      conv_10_line(BUFFER_DEPTH,a[i],b[i]);
  }
  nanosec(&t2);
  br=0;
  for(i=0; i<BUFFER_WIDTH;i++) {
    for(j=0;j<BUFFER_DEPTH;j++) {
      br+=b[i][j];
      if(fabs(b[i][j]-c[i][j])>1e-10){
        printf("error %u %u: %1.15lf != %1.15lf (error %1.15lf)!\n",i,j , b[i][j], c[i][j], fabs(b[i][j]-c[i][j]));
     }
    }
  }
  printf("result 10eo %lf, duration %llu ns, FLOP %d, GFLOPS %lf\n", br, t2-t1, FLOP, (double)FLOP/(float)(t2-t1));

  nanosec(&t1);
  for(i=0; i<BUFFER_WIDTH;i++) {
      conv_10_line_fused(BUFFER_DEPTH,a[i],b[i]);
  }
  nanosec(&t2);
  br=0;
  for(i=0; i<BUFFER_WIDTH;i++) {
    for(j=0;j<BUFFER_DEPTH;j++) {
      br+=b[i][j];
      if(fabs(b[i][j]-c[i][j])>1e-10){
        printf("error %u %u: %1.15lf != %1.15lf (error %1.15lf)!\n",i,j , b[i][j], c[i][j], fabs(b[i][j]-c[i][j]));
     }
    }
  }
  printf("result  10f %lf, duration %llu ns, FLOP %d, GFLOPS %lf\n", br, t2-t1, FLOP, (double)FLOP/(float)(t2-t1));

  nanosec(&t1);
  for(i=0; i<BUFFER_WIDTH;i++) {
      conv_12_line(BUFFER_DEPTH,a[i],b[i]);
  }
  nanosec(&t2);
  br=0;
  for(i=0; i<BUFFER_WIDTH;i++) {
    for(j=0;j<BUFFER_DEPTH;j++) {
      br+=b[i][j];
      if(fabs(b[i][j]-c[i][j])>1e-10){
        printf("error %u %u: %1.15lf != %1.15lf (error %1.15lf)!\n",i,j , b[i][j], c[i][j], fabs(b[i][j]-c[i][j]));
     }
    }
  }
  printf("result 12eo %lf, duration %llu ns, FLOP %d, GFLOPS %lf\n", br, t2-t1, FLOP, (double)FLOP/(float)(t2-t1));

  nanosec(&t1);
  for(i=0; i<BUFFER_WIDTH;i++) {
      conv_12_line_fused(BUFFER_DEPTH,a[i],b[i]);
  }
  nanosec(&t2);
  br=0;
  for(i=0; i<BUFFER_WIDTH;i++) {
    for(j=0;j<BUFFER_DEPTH;j++) {
      br+=b[i][j];
      if(fabs(b[i][j]-c[i][j])>1e-10){
        printf("error %u %u: %1.15lf != %1.15lf (error %1.15lf)!\n",i,j , b[i][j], c[i][j], fabs(b[i][j]-c[i][j]));
     }
    }
  }
  printf("result  12f %lf, duration %llu ns, FLOP %d, GFLOPS %lf\n", br, t2-t1, FLOP, (double)FLOP/(float)(t2-t1));

  nanosec(&t1);
  for(i=0; i<BUFFER_WIDTH;i++) {
      conv_14_line(BUFFER_DEPTH,a[i],b[i]);
  }
  nanosec(&t2);
  br=0;
  for(i=0; i<BUFFER_WIDTH;i++) {
    for(j=0;j<BUFFER_DEPTH;j++) {
      br+=b[i][j];
      if(fabs(b[i][j]-c[i][j])>1e-10){
        printf("error %u %u: %1.15lf != %1.15lf (error %1.15lf)!\n",i,j , b[i][j], c[i][j], fabs(b[i][j]-c[i][j]));
     }
    }
  }
  printf("result 14eo %lf, duration %llu ns, FLOP %d, GFLOPS %lf\n", br, t2-t1, FLOP, (double)FLOP/(float)(t2-t1));

  nanosec(&t1);
  for(i=0; i<BUFFER_WIDTH;i++) {
      conv_16_line(BUFFER_DEPTH,a[i],b[i]);
  }
  nanosec(&t2);
  br=0;
  for(i=0; i<BUFFER_WIDTH;i++) {
    for(j=0;j<BUFFER_DEPTH;j++) {
      br+=b[i][j];
      if(fabs(b[i][j]-c[i][j])>1e-10){
        printf("error %u %u: %1.15lf != %1.15lf (error %1.15lf)!\n",i,j , b[i][j], c[i][j], fabs(b[i][j]-c[i][j]));
     }
    }
  }
  printf("result 16eo %lf, duration %llu ns, FLOP %d, GFLOPS %lf\n", br, t2-t1, FLOP, (double)FLOP/(float)(t2-t1));

  nanosec(&t1);
  for(i=0; i<BUFFER_WIDTH;i++) {
      conv_18_line(BUFFER_DEPTH,a[i],b[i]);
  }
  nanosec(&t2);
  br=0;
  for(i=0; i<BUFFER_WIDTH;i++) {
    for(j=0;j<BUFFER_DEPTH;j++) {
      br+=b[i][j];
      if(fabs(b[i][j]-c[i][j])>1e-10){
        printf("error %u %u: %1.15lf != %1.15lf (error %1.15lf)!\n",i,j , b[i][j], c[i][j], fabs(b[i][j]-c[i][j]));
     }
    }
  }
  printf("result 18eo %lf, duration %llu ns, FLOP %d, GFLOPS %lf\n", br, t2-t1, FLOP, (double)FLOP/(float)(t2-t1));


}
