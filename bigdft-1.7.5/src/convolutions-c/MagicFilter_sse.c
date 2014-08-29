//! @file
//!  Magic filter optimised for sse (x86 vectorisation)
//!
//! @author
//!    Copyright (C) 2009-2011 BigDFT group 
//!    This file is distributed under the terms of the
//!    GNU General Public License, see ~/COPYING file
//!    or http://www.gnu.org/copyleft/gpl.txt .
//!    For the list of contributors, see ~/AUTHORS 


#include <emmintrin.h>
#include <pmmintrin.h>
#include <time.h>


const double filter[] __attribute__ ((aligned (16))) = { 8.4334247333529341094733325815816e-7,
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
const double filter_u[] __attribute__ ((aligned (16))) = { 2.72734492911979659657715313017228e-6,
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

const double filter_reverse[] __attribute__ ((aligned (16))) = {
                        2.72734492911979659657715313017228e-6,
                       -0.5185986881173432922848639136911487e-4,
                        0.49443227688689919192282259476750972e-3,
                       -0.344128144493493857280881509686821861e-2,
                        0.1337263414854794752733423467013220997e-1,
                       -0.2103025160930381434955489412839065067e-1,
                       -0.604895289196983516002834636e-1,
                        0.9940415697834003993178616713,
                        0.612625895831207982195380597e-1,
                        0.2373821463724942397566389712597274535e-1,
                       -0.942047030201080385922711540948195075e-2,
                        0.174723713672993903449447812749852942e-2,
                       -0.30158038132690463167163703826169879e-3,
                        0.8762984476210559564689161894116397e-4,
                       -0.1290557201342060969516786758559028e-4,
                        8.4334247333529341094733325815816e-7
};

const double filter_reverse_u[] __attribute__ ((aligned (16))) = {
                        8.4334247333529341094733325815816e-7,
                        2.72734492911979659657715313017228e-6,
                       -0.5185986881173432922848639136911487e-4,
                        0.49443227688689919192282259476750972e-3,
                       -0.344128144493493857280881509686821861e-2,
                        0.1337263414854794752733423467013220997e-1,
                       -0.2103025160930381434955489412839065067e-1,
                       -0.604895289196983516002834636e-1,
                        0.9940415697834003993178616713,
                        0.612625895831207982195380597e-1,
                        0.2373821463724942397566389712597274535e-1,
                       -0.942047030201080385922711540948195075e-2,
                        0.174723713672993903449447812749852942e-2,
                       -0.30158038132690463167163703826169879e-3,
                        0.8762984476210559564689161894116397e-4,
                       -0.1290557201342060969516786758559028e-4,
                        8.4334247333529341094733325815816e-7
};


#define conv_4x2_block_fused(offset_filter,offset_source,d00,d10,d20,d30) \
FA = _mm_load_pd(filter+offset_filter);\
d00 = _mm_load_pd(source0+offset_source);\
S00 = _mm_add_pd(S00,_mm_mul_pd(d00,FA));\
FU = _mm_load_pd(filter_u+offset_filter);\
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

inline void conv_4x2_fused_init(size_t n, size_t ndat, double const * source0, double const * source1, double const * source2, double const * source3,
                           double       * dest){
  __m128d S00,S01,S10,S11,S20,S21,S30,S31;
  __m128d FA,FU;
  __m128d D00,D10,D20,D30;

  FA = _mm_load_pd(filter);
  D00 = _mm_load_pd(source0+n-8);
  S00 = _mm_mul_pd(D00,FA);
  D10 = _mm_load_pd(source1+n-8);
  S10 = _mm_mul_pd(D10,FA);
  D20 = _mm_load_pd(source2+n-8);
  S20 = _mm_mul_pd(D20,FA);
  D30 = _mm_load_pd(source3+n-8);
  S30 = _mm_mul_pd(D30,FA);
  
  FU = _mm_load_pd(filter_u);
  S01 = _mm_loadl_pd(D00,source0+8);
  S01 = _mm_mul_pd(S01,FU);
  S11 = _mm_loadl_pd(D10,source1+8);
  S11 = _mm_mul_pd(S11,FU);
  S21 = _mm_loadl_pd(D20,source2+8);
  S21 = _mm_mul_pd(S21,FU);
  S31 = _mm_loadl_pd(D30,source3+8);
  S31 = _mm_mul_pd(S31,FU);

  conv_4x2_block_fused(2,n-6,D00,D10,D20,D30);
  conv_4x2_block_fused(4,n-4,D00,D10,D20,D30);
  conv_4x2_block_fused(6,n-2,D00,D10,D20,D30);
  conv_4x2_block_fused(8,0,D00,D10,D20,D30);
  conv_4x2_block_fused(10,2,D00,D10,D20,D30);
  conv_4x2_block_fused(12,4,D00,D10,D20,D30);
  conv_4x2_block_fused(14,6,D00,D10,D20,D30);

  _mm_store_pd(dest,_mm_hadd_pd(S00,S10));
  _mm_store_pd(dest+2,_mm_hadd_pd(S20,S30));
  _mm_store_pd(dest+ndat,_mm_hadd_pd(S01,S11));
  _mm_store_pd(dest+2+ndat,_mm_hadd_pd(S21,S31));
 
  dest+=2*ndat;

  FA = _mm_load_pd(filter);
  D00 = _mm_load_pd(source0+n-6);
  S00 = _mm_mul_pd(D00,FA);
  D10 = _mm_load_pd(source1+n-6);
  S10 = _mm_mul_pd(D10,FA);
  D20 = _mm_load_pd(source2+n-6);
  S20 = _mm_mul_pd(D20,FA);
  D30 = _mm_load_pd(source3+n-6);
  S30 = _mm_mul_pd(D30,FA);
  
  FU = _mm_load_pd(filter_u);
  S01 = _mm_loadl_pd(D00,source0+10);
  S01 = _mm_mul_pd(S01,FU);
  S11 = _mm_loadl_pd(D10,source1+10);
  S11 = _mm_mul_pd(S11,FU);
  S21 = _mm_loadl_pd(D20,source2+10);
  S21 = _mm_mul_pd(S21,FU);
  S31 = _mm_loadl_pd(D30,source3+10);
  S31 = _mm_mul_pd(S31,FU);

  conv_4x2_block_fused(2,n-4,D00,D10,D20,D30);
  conv_4x2_block_fused(4,n-2,D00,D10,D20,D30);
  conv_4x2_block_fused(6,0,D00,D10,D20,D30);
  conv_4x2_block_fused(8,2,D00,D10,D20,D30);
  conv_4x2_block_fused(10,4,D00,D10,D20,D30);
  conv_4x2_block_fused(12,6,D00,D10,D20,D30);
  conv_4x2_block_fused(14,8,D00,D10,D20,D30);

  _mm_store_pd(dest,_mm_hadd_pd(S00,S10));
  _mm_store_pd(dest+2,_mm_hadd_pd(S20,S30));
  _mm_store_pd(dest+ndat,_mm_hadd_pd(S01,S11));
  _mm_store_pd(dest+2+ndat,_mm_hadd_pd(S21,S31));

  dest+=2*ndat;

  FA = _mm_load_pd(filter);
  D00 = _mm_load_pd(source0+n-4);
  S00 = _mm_mul_pd(D00,FA);
  D10 = _mm_load_pd(source1+n-4);
  S10 = _mm_mul_pd(D10,FA);
  D20 = _mm_load_pd(source2+n-4);
  S20 = _mm_mul_pd(D20,FA);
  D30 = _mm_load_pd(source3+n-4);
  S30 = _mm_mul_pd(D30,FA);
  
  FU = _mm_load_pd(filter_u);
  S01 = _mm_loadl_pd(D00,source0+12);
  S01 = _mm_mul_pd(S01,FU);
  S11 = _mm_loadl_pd(D10,source1+12);
  S11 = _mm_mul_pd(S11,FU);
  S21 = _mm_loadl_pd(D20,source2+12);
  S21 = _mm_mul_pd(S21,FU);
  S31 = _mm_loadl_pd(D30,source3+12);
  S31 = _mm_mul_pd(S31,FU);

  conv_4x2_block_fused(2,n-2,D00,D10,D20,D30);
  conv_4x2_block_fused(4,0,D00,D10,D20,D30);
  conv_4x2_block_fused(6,2,D00,D10,D20,D30);
  conv_4x2_block_fused(8,4,D00,D10,D20,D30);
  conv_4x2_block_fused(10,6,D00,D10,D20,D30);
  conv_4x2_block_fused(12,8,D00,D10,D20,D30);
  conv_4x2_block_fused(14,10,D00,D10,D20,D30);

  _mm_store_pd(dest,_mm_hadd_pd(S00,S10));
  _mm_store_pd(dest+2,_mm_hadd_pd(S20,S30));
  _mm_store_pd(dest+ndat,_mm_hadd_pd(S01,S11));
  _mm_store_pd(dest+2+ndat,_mm_hadd_pd(S21,S31));

  dest+=2*ndat;

  FA = _mm_load_pd(filter);
  D00 = _mm_load_pd(source0+n-2);
  S00 = _mm_mul_pd(D00,FA);
  D10 = _mm_load_pd(source1+n-2);
  S10 = _mm_mul_pd(D10,FA);
  D20 = _mm_load_pd(source2+n-2);
  S20 = _mm_mul_pd(D20,FA);
  D30 = _mm_load_pd(source3+n-2);
  S30 = _mm_mul_pd(D30,FA);
  
  FU = _mm_load_pd(filter_u);
  S01 = _mm_loadl_pd(D00,source0+14);
  S01 = _mm_mul_pd(S01,FU);
  S11 = _mm_loadl_pd(D10,source1+14);
  S11 = _mm_mul_pd(S11,FU);
  S21 = _mm_loadl_pd(D20,source2+14);
  S21 = _mm_mul_pd(S21,FU);
  S31 = _mm_loadl_pd(D30,source3+14);
  S31 = _mm_mul_pd(S31,FU);

  conv_4x2_block_fused(2,0,D00,D10,D20,D30);
  conv_4x2_block_fused(4,2,D00,D10,D20,D30);
  conv_4x2_block_fused(6,4,D00,D10,D20,D30);
  conv_4x2_block_fused(8,6,D00,D10,D20,D30);
  conv_4x2_block_fused(10,8,D00,D10,D20,D30);
  conv_4x2_block_fused(12,10,D00,D10,D20,D30);
  conv_4x2_block_fused(14,12,D00,D10,D20,D30);

  _mm_store_pd(dest,_mm_hadd_pd(S00,S10));
  _mm_store_pd(dest+2,_mm_hadd_pd(S20,S30));
  _mm_store_pd(dest+ndat,_mm_hadd_pd(S01,S11));
  _mm_store_pd(dest+2+ndat,_mm_hadd_pd(S21,S31));


}

inline void conv_4x2_fused_finish(size_t n, size_t ndat, double const * source0, double const * source1, double const * source2, double const * source3,
                           double       * dest){
  __m128d S00,S01,S10,S11,S20,S21,S30,S31;
  __m128d FA,FU;
  __m128d D00,D10,D20,D30;

  FA = _mm_load_pd(filter);
  D00 = _mm_load_pd(source0);
  S00 = _mm_mul_pd(D00,FA);
  D10 = _mm_load_pd(source1);
  S10 = _mm_mul_pd(D10,FA);
  D20 = _mm_load_pd(source2);
  S20 = _mm_mul_pd(D20,FA);
  D30 = _mm_load_pd(source3);
  S30 = _mm_mul_pd(D30,FA);
  
  FU = _mm_load_pd(filter_u);
  S01 = _mm_loadl_pd(D00,source0-n+16);
  S01 = _mm_mul_pd(S01,FU);
  S11 = _mm_loadl_pd(D10,source1-n+16);
  S11 = _mm_mul_pd(S11,FU);
  S21 = _mm_loadl_pd(D20,source2-n+16);
  S21 = _mm_mul_pd(S21,FU);
  S31 = _mm_loadl_pd(D30,source3-n+16);
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

  dest +=2*ndat;

  FA = _mm_load_pd(filter);
  D00 = _mm_load_pd(source0+2);
  S00 = _mm_mul_pd(D00,FA);
  D10 = _mm_load_pd(source1+2);
  S10 = _mm_mul_pd(D10,FA);
  D20 = _mm_load_pd(source2+2);
  S20 = _mm_mul_pd(D20,FA);
  D30 = _mm_load_pd(source3+2);
  S30 = _mm_mul_pd(D30,FA);
  
  FU = _mm_load_pd(filter_u);
  S01 = _mm_loadl_pd(D00,source0-n+18);
  S01 = _mm_mul_pd(S01,FU);
  S11 = _mm_loadl_pd(D10,source1-n+18);
  S11 = _mm_mul_pd(S11,FU);
  S21 = _mm_loadl_pd(D20,source2-n+18);
  S21 = _mm_mul_pd(S21,FU);
  S31 = _mm_loadl_pd(D30,source3-n+18);
  S31 = _mm_mul_pd(S31,FU);

  conv_4x2_block_fused(2,4,D00,D10,D20,D30);
  conv_4x2_block_fused(4,6,D00,D10,D20,D30);
  conv_4x2_block_fused(6,8,D00,D10,D20,D30);
  conv_4x2_block_fused(8,10,D00,D10,D20,D30);
  conv_4x2_block_fused(10,12,D00,D10,D20,D30);
  conv_4x2_block_fused(12,14,D00,D10,D20,D30);
  conv_4x2_block_fused(14,-n+16,D00,D10,D20,D30);

  _mm_store_pd(dest,_mm_hadd_pd(S00,S10));
  _mm_store_pd(dest+2,_mm_hadd_pd(S20,S30));
  _mm_store_pd(dest+ndat,_mm_hadd_pd(S01,S11));
  _mm_store_pd(dest+2+ndat,_mm_hadd_pd(S21,S31));

  dest +=2*ndat;

  FA = _mm_load_pd(filter);
  D00 = _mm_load_pd(source0+4);
  S00 = _mm_mul_pd(D00,FA);
  D10 = _mm_load_pd(source1+4);
  S10 = _mm_mul_pd(D10,FA);
  D20 = _mm_load_pd(source2+4);
  S20 = _mm_mul_pd(D20,FA);
  D30 = _mm_load_pd(source3+4);
  S30 = _mm_mul_pd(D30,FA);
  
  FU = _mm_load_pd(filter_u);
  S01 = _mm_loadl_pd(D00,source0-n+20);
  S01 = _mm_mul_pd(S01,FU);
  S11 = _mm_loadl_pd(D10,source1-n+20);
  S11 = _mm_mul_pd(S11,FU);
  S21 = _mm_loadl_pd(D20,source2-n+20);
  S21 = _mm_mul_pd(S21,FU);
  S31 = _mm_loadl_pd(D30,source3-n+20);
  S31 = _mm_mul_pd(S31,FU);

  conv_4x2_block_fused(2,6,D00,D10,D20,D30);
  conv_4x2_block_fused(4,8,D00,D10,D20,D30);
  conv_4x2_block_fused(6,10,D00,D10,D20,D30);
  conv_4x2_block_fused(8,12,D00,D10,D20,D30);
  conv_4x2_block_fused(10,14,D00,D10,D20,D30);
  conv_4x2_block_fused(12,-n+16,D00,D10,D20,D30);
  conv_4x2_block_fused(14,-n+18,D00,D10,D20,D30);

  _mm_store_pd(dest,_mm_hadd_pd(S00,S10));
  _mm_store_pd(dest+2,_mm_hadd_pd(S20,S30));
  _mm_store_pd(dest+ndat,_mm_hadd_pd(S01,S11));
  _mm_store_pd(dest+2+ndat,_mm_hadd_pd(S21,S31));

  dest +=2*ndat;

  FA = _mm_load_pd(filter);
  D00 = _mm_load_pd(source0+6);
  S00 = _mm_mul_pd(D00,FA);
  D10 = _mm_load_pd(source1+6);
  S10 = _mm_mul_pd(D10,FA);
  D20 = _mm_load_pd(source2+6);
  S20 = _mm_mul_pd(D20,FA);
  D30 = _mm_load_pd(source3+6);
  S30 = _mm_mul_pd(D30,FA);
  
  FU = _mm_load_pd(filter_u);
  S01 = _mm_loadl_pd(D00,source0-n+22);
  S01 = _mm_mul_pd(S01,FU);
  S11 = _mm_loadl_pd(D10,source1-n+22);
  S11 = _mm_mul_pd(S11,FU);
  S21 = _mm_loadl_pd(D20,source2-n+22);
  S21 = _mm_mul_pd(S21,FU);
  S31 = _mm_loadl_pd(D30,source3-n+22);
  S31 = _mm_mul_pd(S31,FU);

  conv_4x2_block_fused(2,8,D00,D10,D20,D30);
  conv_4x2_block_fused(4,10,D00,D10,D20,D30);
  conv_4x2_block_fused(6,12,D00,D10,D20,D30);
  conv_4x2_block_fused(8,14,D00,D10,D20,D30);
  conv_4x2_block_fused(10,-n+16,D00,D10,D20,D30);
  conv_4x2_block_fused(12,-n+18,D00,D10,D20,D30);
  conv_4x2_block_fused(14,-n+20,D00,D10,D20,D30);

  _mm_store_pd(dest,_mm_hadd_pd(S00,S10));
  _mm_store_pd(dest+2,_mm_hadd_pd(S20,S30));
  _mm_store_pd(dest+ndat,_mm_hadd_pd(S01,S11));
  _mm_store_pd(dest+2+ndat,_mm_hadd_pd(S21,S31));

}


inline void conv_4x2_fused(size_t ndat, double const * source0, double const * source1, double const * source2, double const * source3,
                           double       * dest){
  __m128d S00,S01,S10,S11,S20,S21,S30,S31;
  __m128d FA,FU;
  __m128d D00,D10,D20,D30;

  FA = _mm_load_pd(filter);
  D00 = _mm_load_pd(source0);
  S00 = _mm_mul_pd(D00,FA);
  D10 = _mm_load_pd(source1);
  S10 = _mm_mul_pd(D10,FA);
  D20 = _mm_load_pd(source2);
  S20 = _mm_mul_pd(D20,FA);
  D30 = _mm_load_pd(source3);
  S30 = _mm_mul_pd(D30,FA);
  
  FU = _mm_load_pd(filter_u);
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

#define conv_4x2_block_fused_t(offset_filter,offset_source,d00,d10,d20,d30) \
FA = _mm_load_pd(filter_reverse+offset_filter);\
d00 = _mm_load_pd(source0+offset_source);\
S00 = _mm_add_pd(S00,_mm_mul_pd(d00,FA));\
FU = _mm_load_pd(filter_reverse_u+offset_filter);\
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

inline void conv_4x2_fused_init_t(size_t n, size_t ndat, double const * source0, double const * source1, double const * source2, double const * source3,
                           double       * dest){
  __m128d S00,S01,S10,S11,S20,S21,S30,S31;
  __m128d FA,FU;
  __m128d D00,D10,D20,D30;

  FA = _mm_load_pd(filter_reverse);
  D00 = _mm_load_pd(source0+n-8);
  S00 = _mm_mul_pd(D00,FA);
  D10 = _mm_load_pd(source1+n-8);
  S10 = _mm_mul_pd(D10,FA);
  D20 = _mm_load_pd(source2+n-8);
  S20 = _mm_mul_pd(D20,FA);
  D30 = _mm_load_pd(source3+n-8);
  S30 = _mm_mul_pd(D30,FA);
  
  FU = _mm_load_pd(filter_reverse_u);
  S01 = _mm_loadl_pd(D00,source0+8);
  S01 = _mm_mul_pd(S01,FU);
  S11 = _mm_loadl_pd(D10,source1+8);
  S11 = _mm_mul_pd(S11,FU);
  S21 = _mm_loadl_pd(D20,source2+8);
  S21 = _mm_mul_pd(S21,FU);
  S31 = _mm_loadl_pd(D30,source3+8);
  S31 = _mm_mul_pd(S31,FU);

  conv_4x2_block_fused_t(2,n-6,D00,D10,D20,D30);
  conv_4x2_block_fused_t(4,n-4,D00,D10,D20,D30);
  conv_4x2_block_fused_t(6,n-2,D00,D10,D20,D30);
  conv_4x2_block_fused_t(8,0,D00,D10,D20,D30);
  conv_4x2_block_fused_t(10,2,D00,D10,D20,D30);
  conv_4x2_block_fused_t(12,4,D00,D10,D20,D30);
  conv_4x2_block_fused_t(14,6,D00,D10,D20,D30);

  _mm_store_pd(dest+(n-1)*ndat,_mm_hadd_pd(S00,S10));
  _mm_store_pd(dest+(n-1)*ndat+2,_mm_hadd_pd(S20,S30));
  _mm_store_pd(dest,_mm_hadd_pd(S01,S11));
  _mm_store_pd(dest+2,_mm_hadd_pd(S21,S31));
 
  dest+=ndat;

  FA = _mm_load_pd(filter_reverse);
  D00 = _mm_load_pd(source0+n-6);
  S00 = _mm_mul_pd(D00,FA);
  D10 = _mm_load_pd(source1+n-6);
  S10 = _mm_mul_pd(D10,FA);
  D20 = _mm_load_pd(source2+n-6);
  S20 = _mm_mul_pd(D20,FA);
  D30 = _mm_load_pd(source3+n-6);
  S30 = _mm_mul_pd(D30,FA);
  
  FU = _mm_load_pd(filter_reverse_u);
  S01 = _mm_loadl_pd(D00,source0+10);
  S01 = _mm_mul_pd(S01,FU);
  S11 = _mm_loadl_pd(D10,source1+10);
  S11 = _mm_mul_pd(S11,FU);
  S21 = _mm_loadl_pd(D20,source2+10);
  S21 = _mm_mul_pd(S21,FU);
  S31 = _mm_loadl_pd(D30,source3+10);
  S31 = _mm_mul_pd(S31,FU);

  conv_4x2_block_fused_t(2,n-4,D00,D10,D20,D30);
  conv_4x2_block_fused_t(4,n-2,D00,D10,D20,D30);
  conv_4x2_block_fused_t(6,0,D00,D10,D20,D30);
  conv_4x2_block_fused_t(8,2,D00,D10,D20,D30);
  conv_4x2_block_fused_t(10,4,D00,D10,D20,D30);
  conv_4x2_block_fused_t(12,6,D00,D10,D20,D30);
  conv_4x2_block_fused_t(14,8,D00,D10,D20,D30);

  _mm_store_pd(dest,_mm_hadd_pd(S00,S10));
  _mm_store_pd(dest+2,_mm_hadd_pd(S20,S30));
  _mm_store_pd(dest+ndat,_mm_hadd_pd(S01,S11));
  _mm_store_pd(dest+2+ndat,_mm_hadd_pd(S21,S31));

  dest+=2*ndat;

  FA = _mm_load_pd(filter_reverse);
  D00 = _mm_load_pd(source0+n-4);
  S00 = _mm_mul_pd(D00,FA);
  D10 = _mm_load_pd(source1+n-4);
  S10 = _mm_mul_pd(D10,FA);
  D20 = _mm_load_pd(source2+n-4);
  S20 = _mm_mul_pd(D20,FA);
  D30 = _mm_load_pd(source3+n-4);
  S30 = _mm_mul_pd(D30,FA);
  
  FU = _mm_load_pd(filter_reverse_u);
  S01 = _mm_loadl_pd(D00,source0+12);
  S01 = _mm_mul_pd(S01,FU);
  S11 = _mm_loadl_pd(D10,source1+12);
  S11 = _mm_mul_pd(S11,FU);
  S21 = _mm_loadl_pd(D20,source2+12);
  S21 = _mm_mul_pd(S21,FU);
  S31 = _mm_loadl_pd(D30,source3+12);
  S31 = _mm_mul_pd(S31,FU);

  conv_4x2_block_fused_t(2,n-2,D00,D10,D20,D30);
  conv_4x2_block_fused_t(4,0,D00,D10,D20,D30);
  conv_4x2_block_fused_t(6,2,D00,D10,D20,D30);
  conv_4x2_block_fused_t(8,4,D00,D10,D20,D30);
  conv_4x2_block_fused_t(10,6,D00,D10,D20,D30);
  conv_4x2_block_fused_t(12,8,D00,D10,D20,D30);
  conv_4x2_block_fused_t(14,10,D00,D10,D20,D30);

  _mm_store_pd(dest,_mm_hadd_pd(S00,S10));
  _mm_store_pd(dest+2,_mm_hadd_pd(S20,S30));
  _mm_store_pd(dest+ndat,_mm_hadd_pd(S01,S11));
  _mm_store_pd(dest+2+ndat,_mm_hadd_pd(S21,S31));

  dest+=2*ndat;

  FA = _mm_load_pd(filter_reverse);
  D00 = _mm_load_pd(source0+n-2);
  S00 = _mm_mul_pd(D00,FA);
  D10 = _mm_load_pd(source1+n-2);
  S10 = _mm_mul_pd(D10,FA);
  D20 = _mm_load_pd(source2+n-2);
  S20 = _mm_mul_pd(D20,FA);
  D30 = _mm_load_pd(source3+n-2);
  S30 = _mm_mul_pd(D30,FA);
  
  FU = _mm_load_pd(filter_reverse_u);
  S01 = _mm_loadl_pd(D00,source0+14);
  S01 = _mm_mul_pd(S01,FU);
  S11 = _mm_loadl_pd(D10,source1+14);
  S11 = _mm_mul_pd(S11,FU);
  S21 = _mm_loadl_pd(D20,source2+14);
  S21 = _mm_mul_pd(S21,FU);
  S31 = _mm_loadl_pd(D30,source3+14);
  S31 = _mm_mul_pd(S31,FU);

  conv_4x2_block_fused_t(2,0,D00,D10,D20,D30);
  conv_4x2_block_fused_t(4,2,D00,D10,D20,D30);
  conv_4x2_block_fused_t(6,4,D00,D10,D20,D30);
  conv_4x2_block_fused_t(8,6,D00,D10,D20,D30);
  conv_4x2_block_fused_t(10,8,D00,D10,D20,D30);
  conv_4x2_block_fused_t(12,10,D00,D10,D20,D30);
  conv_4x2_block_fused_t(14,12,D00,D10,D20,D30);

  _mm_store_pd(dest,_mm_hadd_pd(S00,S10));
  _mm_store_pd(dest+2,_mm_hadd_pd(S20,S30));
  _mm_store_pd(dest+ndat,_mm_hadd_pd(S01,S11));
  _mm_store_pd(dest+2+ndat,_mm_hadd_pd(S21,S31));


}

inline void conv_4x2_fused_finish_t(size_t n, size_t ndat, double const * source0, double const * source1, double const * source2, double const * source3,
                           double       * dest){
  __m128d S00,S01,S10,S11,S20,S21,S30,S31;
  __m128d FA,FU;
  __m128d D00,D10,D20,D30;

  FA = _mm_load_pd(filter_reverse);
  D00 = _mm_load_pd(source0);
  S00 = _mm_mul_pd(D00,FA);
  D10 = _mm_load_pd(source1);
  S10 = _mm_mul_pd(D10,FA);
  D20 = _mm_load_pd(source2);
  S20 = _mm_mul_pd(D20,FA);
  D30 = _mm_load_pd(source3);
  S30 = _mm_mul_pd(D30,FA);
  
  FU = _mm_load_pd(filter_reverse_u);
  S01 = _mm_loadl_pd(D00,source0-n+16);
  S01 = _mm_mul_pd(S01,FU);
  S11 = _mm_loadl_pd(D10,source1-n+16);
  S11 = _mm_mul_pd(S11,FU);
  S21 = _mm_loadl_pd(D20,source2-n+16);
  S21 = _mm_mul_pd(S21,FU);
  S31 = _mm_loadl_pd(D30,source3-n+16);
  S31 = _mm_mul_pd(S31,FU);

  conv_4x2_block_fused_t(2,2,D00,D10,D20,D30);
  conv_4x2_block_fused_t(4,4,D00,D10,D20,D30);
  conv_4x2_block_fused_t(6,6,D00,D10,D20,D30);
  conv_4x2_block_fused_t(8,8,D00,D10,D20,D30);
  conv_4x2_block_fused_t(10,10,D00,D10,D20,D30);
  conv_4x2_block_fused_t(12,12,D00,D10,D20,D30);
  conv_4x2_block_fused_t(14,14,D00,D10,D20,D30);

  _mm_store_pd(dest,_mm_hadd_pd(S00,S10));
  _mm_store_pd(dest+2,_mm_hadd_pd(S20,S30));
  _mm_store_pd(dest+ndat,_mm_hadd_pd(S01,S11));
  _mm_store_pd(dest+2+ndat,_mm_hadd_pd(S21,S31));

  dest +=2*ndat;

  FA = _mm_load_pd(filter_reverse);
  D00 = _mm_load_pd(source0+2);
  S00 = _mm_mul_pd(D00,FA);
  D10 = _mm_load_pd(source1+2);
  S10 = _mm_mul_pd(D10,FA);
  D20 = _mm_load_pd(source2+2);
  S20 = _mm_mul_pd(D20,FA);
  D30 = _mm_load_pd(source3+2);
  S30 = _mm_mul_pd(D30,FA);
  
  FU = _mm_load_pd(filter_reverse_u);
  S01 = _mm_loadl_pd(D00,source0-n+18);
  S01 = _mm_mul_pd(S01,FU);
  S11 = _mm_loadl_pd(D10,source1-n+18);
  S11 = _mm_mul_pd(S11,FU);
  S21 = _mm_loadl_pd(D20,source2-n+18);
  S21 = _mm_mul_pd(S21,FU);
  S31 = _mm_loadl_pd(D30,source3-n+18);
  S31 = _mm_mul_pd(S31,FU);

  conv_4x2_block_fused_t(2,4,D00,D10,D20,D30);
  conv_4x2_block_fused_t(4,6,D00,D10,D20,D30);
  conv_4x2_block_fused_t(6,8,D00,D10,D20,D30);
  conv_4x2_block_fused_t(8,10,D00,D10,D20,D30);
  conv_4x2_block_fused_t(10,12,D00,D10,D20,D30);
  conv_4x2_block_fused_t(12,14,D00,D10,D20,D30);
  conv_4x2_block_fused_t(14,-n+16,D00,D10,D20,D30);

  _mm_store_pd(dest,_mm_hadd_pd(S00,S10));
  _mm_store_pd(dest+2,_mm_hadd_pd(S20,S30));
  _mm_store_pd(dest+ndat,_mm_hadd_pd(S01,S11));
  _mm_store_pd(dest+2+ndat,_mm_hadd_pd(S21,S31));

  dest +=2*ndat;

  FA = _mm_load_pd(filter_reverse);
  D00 = _mm_load_pd(source0+4);
  S00 = _mm_mul_pd(D00,FA);
  D10 = _mm_load_pd(source1+4);
  S10 = _mm_mul_pd(D10,FA);
  D20 = _mm_load_pd(source2+4);
  S20 = _mm_mul_pd(D20,FA);
  D30 = _mm_load_pd(source3+4);
  S30 = _mm_mul_pd(D30,FA);
  
  FU = _mm_load_pd(filter_reverse_u);
  S01 = _mm_loadl_pd(D00,source0-n+20);
  S01 = _mm_mul_pd(S01,FU);
  S11 = _mm_loadl_pd(D10,source1-n+20);
  S11 = _mm_mul_pd(S11,FU);
  S21 = _mm_loadl_pd(D20,source2-n+20);
  S21 = _mm_mul_pd(S21,FU);
  S31 = _mm_loadl_pd(D30,source3-n+20);
  S31 = _mm_mul_pd(S31,FU);

  conv_4x2_block_fused_t(2,6,D00,D10,D20,D30);
  conv_4x2_block_fused_t(4,8,D00,D10,D20,D30);
  conv_4x2_block_fused_t(6,10,D00,D10,D20,D30);
  conv_4x2_block_fused_t(8,12,D00,D10,D20,D30);
  conv_4x2_block_fused_t(10,14,D00,D10,D20,D30);
  conv_4x2_block_fused_t(12,-n+16,D00,D10,D20,D30);
  conv_4x2_block_fused_t(14,-n+18,D00,D10,D20,D30);

  _mm_store_pd(dest,_mm_hadd_pd(S00,S10));
  _mm_store_pd(dest+2,_mm_hadd_pd(S20,S30));
  _mm_store_pd(dest+ndat,_mm_hadd_pd(S01,S11));
  _mm_store_pd(dest+2+ndat,_mm_hadd_pd(S21,S31));

  dest +=2*ndat;

  FA = _mm_load_pd(filter_reverse);
  D00 = _mm_load_pd(source0+6);
  S00 = _mm_mul_pd(D00,FA);
  D10 = _mm_load_pd(source1+6);
  S10 = _mm_mul_pd(D10,FA);
  D20 = _mm_load_pd(source2+6);
  S20 = _mm_mul_pd(D20,FA);
  D30 = _mm_load_pd(source3+6);
  S30 = _mm_mul_pd(D30,FA);
  
  FU = _mm_load_pd(filter_reverse_u);
  S01 = _mm_loadl_pd(D00,source0-n+22);
  S01 = _mm_mul_pd(S01,FU);
  S11 = _mm_loadl_pd(D10,source1-n+22);
  S11 = _mm_mul_pd(S11,FU);
  S21 = _mm_loadl_pd(D20,source2-n+22);
  S21 = _mm_mul_pd(S21,FU);
  S31 = _mm_loadl_pd(D30,source3-n+22);
  S31 = _mm_mul_pd(S31,FU);

  conv_4x2_block_fused_t(2,8,D00,D10,D20,D30);
  conv_4x2_block_fused_t(4,10,D00,D10,D20,D30);
  conv_4x2_block_fused_t(6,12,D00,D10,D20,D30);
  conv_4x2_block_fused_t(8,14,D00,D10,D20,D30);
  conv_4x2_block_fused_t(10,-n+16,D00,D10,D20,D30);
  conv_4x2_block_fused_t(12,-n+18,D00,D10,D20,D30);
  conv_4x2_block_fused_t(14,-n+20,D00,D10,D20,D30);

  _mm_store_pd(dest,_mm_hadd_pd(S00,S10));
  _mm_store_pd(dest+2,_mm_hadd_pd(S20,S30));
  _mm_store_pd(dest+ndat,_mm_hadd_pd(S01,S11));
  _mm_store_pd(dest+2+ndat,_mm_hadd_pd(S21,S31));

}


inline void conv_4x2_fused_t(size_t ndat, double const * source0, double const * source1, double const * source2, double const * source3,
                           double       * dest){
  __m128d S00,S01,S10,S11,S20,S21,S30,S31;
  __m128d FA,FU;
  __m128d D00,D10,D20,D30;

  FA = _mm_load_pd(filter_reverse);
  D00 = _mm_load_pd(source0);
  S00 = _mm_mul_pd(D00,FA);
  D10 = _mm_load_pd(source1);
  S10 = _mm_mul_pd(D10,FA);
  D20 = _mm_load_pd(source2);
  S20 = _mm_mul_pd(D20,FA);
  D30 = _mm_load_pd(source3);
  S30 = _mm_mul_pd(D30,FA);
  
  FU = _mm_load_pd(filter_reverse_u);
  S01 = _mm_loadl_pd(D00,source0+16);
  S01 = _mm_mul_pd(S01,FU);
  S11 = _mm_loadl_pd(D10,source1+16);
  S11 = _mm_mul_pd(S11,FU);
  S21 = _mm_loadl_pd(D20,source2+16);
  S21 = _mm_mul_pd(S21,FU);
  S31 = _mm_loadl_pd(D30,source3+16);
  S31 = _mm_mul_pd(S31,FU);

  conv_4x2_block_fused_t(2,2,D00,D10,D20,D30);
  conv_4x2_block_fused_t(4,4,D00,D10,D20,D30);
  conv_4x2_block_fused_t(6,6,D00,D10,D20,D30);
  conv_4x2_block_fused_t(8,8,D00,D10,D20,D30);
  conv_4x2_block_fused_t(10,10,D00,D10,D20,D30);
  conv_4x2_block_fused_t(12,12,D00,D10,D20,D30);
  conv_4x2_block_fused_t(14,14,D00,D10,D20,D30);

  _mm_store_pd(dest,_mm_hadd_pd(S00,S10));
  _mm_store_pd(dest+2,_mm_hadd_pd(S20,S30));
  _mm_store_pd(dest+ndat,_mm_hadd_pd(S01,S11));
  _mm_store_pd(dest+2+ndat,_mm_hadd_pd(S21,S31));
}

void magicfilter1d_naive_(unsigned int *n, unsigned int *ndat, double const *source, double *dest) {
  double tmp;
  unsigned int i,j,k;
  for(i=0;i<(*ndat);i++){
    for(j=0;j<(*n);j++) {
      tmp=0;
      for(k=0;k<16;k++){
        tmp+=source[(j-8+k+(*n))%(*n)]*filter[k];
      }
      dest[j*(*ndat)]=tmp;
    }
    dest += 1;
    source += (*n);
  } 
}

void magicfilter1d_t_naive_(unsigned int *n, unsigned int *ndat, double const *source, double *dest) {
  double tmp;
  unsigned int i,j,k;
  for(i=0;i<(*ndat);i++){
    for(j=0;j<(*n);j++) {
      tmp=0;
      for(k=0;k<16;k++){
        tmp+=source[(j-7+k+(*n))%(*n)]*filter_reverse[k];
      }
      dest[j*(*ndat)]=tmp;
    }
    dest += 1;
    source += (*n);
  } 
}

void magicfilter1d_sse_(unsigned int *n, unsigned int *ndat, double const *source, double *dest) {
  double * dest_t;
  unsigned int i=0;
  do{
    dest_t=dest;
    conv_4x2_fused_init(*n,*ndat,source,source+*n,source+2*(*n),source+3*(*n),dest_t);
    dest_t += 8*(*ndat);
    unsigned int j=8;
    do {
      conv_4x2_fused(*ndat,source,source+*n,source+2*(*n),source+3*(*n),dest_t);
      source+=2;
      j+=2;
      dest_t+=2*(*ndat);
    } while(j<(*n-8));
    conv_4x2_fused_finish(*n,*ndat,source,source+*n,source+2*(*n),source+3*(*n),dest_t);
    source+=16+3*(*n);
    i+=4;
    dest+=4;
  } while (i<(*ndat));
}

void magicfilter1d_t_sse_(unsigned int *n, unsigned int *ndat, double const *source, double *dest) {
  double * dest_t;
  unsigned int i=0;
  do{
    dest_t=dest;
    conv_4x2_fused_init_t(*n,*ndat,source,source+*n,source+2*(*n),source+3*(*n),dest_t);
    dest_t += 7*(*ndat);
    unsigned int j=8;
    do {
      conv_4x2_fused_t(*ndat,source,source+*n,source+2*(*n),source+3*(*n),dest_t);
      source+=2;
      j+=2;
      dest_t+=2*(*ndat);
    } while(j<(*n-8));
    conv_4x2_fused_finish_t(*n,*ndat,source,source+*n,source+2*(*n),source+3*(*n),dest_t);
    source+=16+3*(*n);
    i+=4;
    dest+=4;
  } while (i<(*ndat));
}
