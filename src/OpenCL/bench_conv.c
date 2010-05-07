#include "bench_lib.h"


#define MAX_N1 256
#define MAX_N2 256
#define MAX_N3 256
#define N1_STEP 64
#define N2_STEP 64
#define N3_STEP 64
#define N1_URANGE 16 
#define N2_URANGE 16 
#define N3_URANGE 0
#define N1_USTEP 8
#define N2_USTEP 1
#define N3_USTEP 1

int main(){
  cl_uint n1,n2,n3;
  cl_uint un1,un2,un3;
  double *in, *out;
  cl_uint size = (MAX_N1+N1_URANGE)*(MAX_N2+N2_URANGE)*(MAX_N3+N3_URANGE);


  in = (double*) malloc(size*sizeof(double));
  out = (double*) malloc(size*sizeof(double));

  init_random(in, size);

  ocl_create_gpu_context_(&context);
  ocl_create_command_queue_(&queue,&context);
  ocl_build_kernels_(&context);
  init_event_list_();

  for( n1 = N1_STEP; n1 <= MAX_N1; n1 += N1_STEP ){
    for( un1 = n1 - N1_URANGE; un1 <= n1 + N1_URANGE; un1 += N1_USTEP){
      printf("%u\n",un1);
      for( n2 = N2_STEP; n2 <= MAX_N2; n2 += N2_STEP ){
        for( un2 = n2 - N2_URANGE; un2 <= n2 + N2_URANGE; un2 += N2_USTEP){
          for( n3 = n2; n3 <= MAX_N3; n3 += N3_STEP ){
            for( un3 = n3 - N3_URANGE; un3 <= n3 + N3_URANGE; un3 += N3_USTEP){
              bench_magicfilter1d(un1,un2,un3,in,out);
              bench_magicfilter1d_straight(un1,un2,un3,in,out);
              bench_magicfiltergrow1d(un1,un2,un3,in,out);
              bench_magicfiltershrink1d(un1,un2,un3,in,out);
              bench_kinetic1d(un1,un2,un3,in,out);
              bench_ana1d(un1,un2,un3,in,out);
              bench_anashrink1d(un1,un2,un3,in,out);
              bench_syn1d(un1,un2,un3,in,out);
              bench_syngrow1d(un1,un2,un3,in,out);
            }
          }
        }
      }
    }
  }

  print_event_list_();
  ocl_clean_(&queue, &context);
  return 0;
}
