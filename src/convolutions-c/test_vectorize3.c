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
const double filt_u[] __attribute__ ((aligned (16))) = { 0.0,
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
#define BUFFER_WIDTH 512
#define BUFFER_DEPTH 256

inline void conv_16_even(double const * source, double * dest){
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

inline void conv_16_odd(double const * source, double * dest){
  __m128d S0,S1,S2,S3,S4,S5,S6,S7;
  __m128d F;
  __m128d D;
  F = _mm_set_pd(filt[0],filt[15]);
  S0 = _mm_set_pd(*(source+1),*(source+16));
  S0 = _mm_mul_pd(S0,F);
  S1 = _mm_set_pd(*(source+2),*(source+18));
  S1 = _mm_mul_pd(S1,F);
  S2 = _mm_set_pd(*(source+3),*(source+20));
  S2 = _mm_mul_pd(S2,F);
  S3 = _mm_set_pd(*(source+4),*(source+22));
  S3 = _mm_mul_pd(S3,F);
  S4 = _mm_set_pd(*(source+5),*(source+24));
  S4 = _mm_mul_pd(S4,F);
  S5 = _mm_set_pd(*(source+6),*(source+26));
  S5 = _mm_mul_pd(S5,F);
  S6 = _mm_set_pd(*(source+7),*(source+28));
  S6 = _mm_mul_pd(S6,F);
  S7 = _mm_set_pd(*(source+8),*(source+30));
  S7 = _mm_mul_pd(S7,F);
  
/*  S0 = _mm_set_sd(*(source+1)*filt[0]+*(source+16)*filt[15]);
  S1 = _mm_set_sd(*(source+3)*filt[0]+*(source+18)*filt[15]);
  S2 = _mm_set_sd(*(source+5)*filt[0]+*(source+20)*filt[15]);
  S3 = _mm_set_sd(*(source+7)*filt[0]+*(source+22)*filt[15]);
  S4 = _mm_set_sd(*(source+9)*filt[0]+*(source+24)*filt[15]);
  S5 = _mm_set_sd(*(source+11)*filt[0]+*(source+26)*filt[15]);
  S6 = _mm_set_sd(*(source+13)*filt[0]+*(source+28)*filt[15]);
  S7 = _mm_set_sd(*(source+15)*filt[0]+*(source+30)*filt[15]);*/
//  source+=1;
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
    for(j=0; j<BUFFER_DEPTH; j+=16){
      conv_ref(&(a[i][j]),&(c[i][j]));
    }
  }
  unsigned long long int t1,t2;
  nanosec(&t1);
  for(i=0; i<BUFFER_WIDTH;i++) {
    for(j=0; j<BUFFER_DEPTH; j+=16){
      conv_16_even(&(a[i][j]),&(b[i][j]));
      conv_16_odd(&(a[i][j]),&(b[i][j]));
    }
  }
  nanosec(&t2);
  double br=0;
  for(i=0; i<BUFFER_WIDTH;i++) {
    for(j=0;j<BUFFER_DEPTH;j++) {
      br+=b[i][j];
      if(fabs(b[i][j]-c[i][j])>1e-10){
//        printf("error %u %u: %lf != %lf!\n",i,j , b[i][j], c[i][j]);
     }
    }
  }
  printf("result %lf, duration %llu ns, FLOP %d, GFLOPS %lf\n", br, t2-t1, FLOP, (double)FLOP/(float)(t2-t1));
}
