#include <time.h>
#include "config.h"
void FC_FUNC_(nanosec_cuda,NANOSEC_CUDA)(unsigned long long int * t){
  struct timespec time;
  clock_gettime(CLOCK_REALTIME, &time);
  *t = time.tv_sec;
  *t *= 1000000000;
  *t += time.tv_nsec;
}
