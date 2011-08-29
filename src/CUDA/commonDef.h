#ifndef _commonDefh_
#define _commonDefh_

#define max(a,b) (a > b ? a : b)
#define min(a,b) (a < b ? a : b)

//parameters for convolutions
//maximum size of the shared memory array
//conceived for maximize occupancy on a hardware of compute
//capability 1.2 and higher (1024 threads at same time on a given multiprocessor)
#define MAX_SHARED_SIZE 3072 //16*256 4 kB (should be =~ 3.9 kB, try also 3072)
#define HALF_WARP_SIZE 16 // for all architectures
#define NUM_LINES 8 
#define HW_ELEM 2 //this is HALF_WARP_SIZE/NUM_LINES


//parameters for compression-decompression
#define MAX_CONSTANT_SIZE 32768 //The overall available memory is 64KB
#define ELEMS_BLOCK 256  // maximum size of a segment for de-compression
#define NBLOCKS_MAX 65536 //for all architectures


void correctSequence(int thds,int elem,int * tab);

//structure for convolutions
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

//structure for de-compression
typedef struct  _keysGPU
{
  unsigned int nsegs;

  int keys[MAX_CONSTANT_SIZE/sizeof(int)]; //array containing the grid information 

} keysGPU_t;


#endif
