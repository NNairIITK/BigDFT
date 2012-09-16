#ifndef BENCH_LIB_H
#define BENCH_LIB_H
#include "OpenCL_wrappers.h"

extern bigdft_context context;
extern bigdft_command_queue queue;

void init_key(cl_int *keyg, cl_int *keyv, cl_int nseg, cl_int *nvctr_cf);
void init_random(double * data, size_t size);
void bench_magicfilter1d(cl_uint n1, cl_uint n2, cl_uint n3, double * in, double * out);
void bench_magicfilter1d_straight(cl_uint n1, cl_uint n2, cl_uint n3, double * in, double * out);
void bench_magicfilter1d_block(cl_uint n1, cl_uint n2, cl_uint n3, double * in, double * out);
void bench_magicfiltershrink1d(cl_uint n1, cl_uint n2, cl_uint n3, double * in, double * out);
void bench_magicfiltergrow1d(cl_uint n1, cl_uint n2, cl_uint n3, double * in, double * out);
void bench_kinetic1d(cl_uint n1, cl_uint n2, cl_uint n3, double * in, double * out);
void bench_ana1d(cl_uint n1, cl_uint n2, cl_uint n3, double * in, double * out);
void bench_ana1d_block(cl_uint n1, cl_uint n2, cl_uint n3, double * in, double * out);
void bench_anashrink1d(cl_uint n1, cl_uint n2, cl_uint n3, double * in, double * out);
void bench_syn1d(cl_uint n1, cl_uint n2, cl_uint n3, double * in, double * out);
void bench_syngrow1d(cl_uint n1, cl_uint n2, cl_uint n3, double * in, double * out);
void bench_uncompress(cl_uint n1, cl_uint n2, cl_uint n3, cl_uint nseg, cl_uint nvctr_cf, cl_uint * keyg, cl_uint * keyv, double * psi_in, double * psi_out);
void bench_gemm(cl_uint n1, cl_uint n2, cl_uint n3, double * in1, double * in2, double * out);
void bench_zgemm(cl_uint n1, cl_uint n2, cl_uint n3, double * in1, double * in2, double * out);
void bench_zgemmd(cl_uint n1, cl_uint n2, cl_uint n3, double * in1, double * in2, double * out);

#endif
