#ifndef _lochamh_
#define _lochamh_

#define max(a,b) (a > b ? a : b)
#define min(a,b) (a < b ? a : b)

//maximum size of the shared memory array
//conceived for maximize occupancy on a hardware of compute
//capability 1.2 and higher (1024 threads at same time on a given multiprocessor)
#define MAX_SHARED_SIZE 3072 //16*256 4 kB (should be =~ 3.9 kB, try also 3072)
#define HALF_WARP_SIZE 16 // for all architectures
#define NUM_LINES 8 
#define HW_ELEM 2 //this is HALF_WARP_SIZE/NUM_LINES


void correctSequence(int thds,int elem,int * tab);

typedef struct  _parGPU
{
  unsigned int ElementsPerBlock;

  int thline[HALF_WARP_SIZE]; //line considered by a thread within the half-warp
  int thelem[HALF_WARP_SIZE]; //elements considered by a thread within the half-warp
  int hwelem_calc[16]; //maximum number of half warps
  int hwelem_copy[16]; //maximum number of half-warps
  int hwoffset_calc[16]; //maximum number of half warps
  int hwoffset_copy[16]; //maximum number of half-warps

} parGPU_t;

extern "C" 
void gpulocham_(int *n1,int *n2, int *n3,
		double *h1,double *h2,double *h3,
		double **psi,double **out,double **pot,			     
		double **work,
		double **work2,
		double *epot,double *ekinpot);

#endif
