/****u* CUDA/newcase.cu
**
** 
** AUTHOR
**  Matthieu Ospici
**
** SOURCE
*/

//#include <stdio.h>
//#include <cutil.h>
//#include <multithreading.h>
//#include <pthread.h>
//#include <semaphore.h>
//#include <sched.h>

//#include <time.h>
//#include <sys/time.h>

//#include "convolution.h"

//#include "conv_pot_shared_kernel.cu"
//#include "deffct.h" //convseria

//#define NUM_MULTIPROCESSORS 30 // valid for GTX280
//#define BLOCK_SIZE 256 //must be at least 64

#define max(a,b) (a > b ? a : b)
#define min(a,b) (a < b ? a : b)

__global__ void conv_lg(unsigned int n,unsigned int ndat,float *psi_in,float *psi_out,
			int lowfil,int lupfil)
{

  const unsigned int tid = threadIdx.x; //thread id
  const unsigned int bid = blockIdx.x;  //block id (can be 2-dim)

  const unsigned int numthreads = blockDim.x;

  //number of parts in which the line is separated
  int linecuts = ((n-1)/numthreads + 1);

  //line cut in which the block is placed
  int icut = bid % linecuts; // to be controlled if modulo operator is
			     // not too slow

  //start point in n axis of the line treated by the block
  int istart = icut * (n/linecuts);

  //end point
  int iend = min((icut+1)*(n/linecuts),n);

  //number of elements
  int num_elem = iend-istart+1;

  //actual column treated by the block, in ndat axis
  int idat = bid / linecuts;

  //starting element in the global memory for the OUTPUT array
  int startelem = n*idat + istart;

  //distance from the halfwarp base element to guarantee coalesced
  //access in global memory
  int startd = startelem & 15;  // bitwise modulo 16
  
  //the first element which should be written in the output array
  //should be processed by thread startd
  //dimension of the array in shared memory (to be half-warp aligned)
  const unsigned int dim_shared = max(numthreads,num_elem+startd+lowfil+lupfil);

  __shared__ float psi_sh[dim_shared];  
  
  //in case the line cut is not the first, correct the starting point
  //such that the left elements are copied also
  //NOTE: it is assumed that for non-first segments the starting
  //points is far enough for the filter to be contained
  int startdcorr=(icut > 1 ? lowfil : 0); 

  //copy of elements in shared memory
  //if (tid >= startd - startdcorr && tid <= startd + num_elem)
  //  {
      //bank conflicts for psi_sh: for each half-warp there is 
      //a linear addressing with stride 1, so ok.

      //coalesced accesses for psi_in: the access are completely
      //uncoalesced, must pass through texture fetches

  psi_sh[tid+lowfil]=psi_in[ndat*(istart+tid-startd)+idat];
      //  }
  
  /* this part cane be neglected for the moment
  //copy the other values in shared memory
  //first start with the values on left side, if not already been copied
  if (icut == 0 && tid)
    {
      
    }

  //next the values of the right side
  if (

  */

  //end shared memory copy
  __syncthreads();


  //perform convolution in shared memory and write results in the
  //lowfil-scaled address

  psi_sh[tid]=psi_sh[tid+lowfil];

  /*
  //write in global memory by taking care of the coalescing
  if (tid >= startd && tid <= startd + num_elem)
    {
      //psi_out[startelem-startd+tid]=psi_sh[tid];
    }
  */

}
/****/
