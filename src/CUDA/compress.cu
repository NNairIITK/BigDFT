/*!
 @author
    Copyright (C) 2010-2011 BigDFT group (LG)
    This file is distributed under the terms of the
    GNU General Public License, see ~/COPYING file
    or http://www.gnu.org/copyleft/gpl.txt .
    For the list of contributors, see ~/AUTHORS 

*/

#include <stdio.h>
#include "commonDef.h"


//parameters for the kernel (global variables)

#include "kernels_compress.hcu"


template<typename T>
int uncompressgpu(int n1, int n2, int n3,
		  T *psicf,T *psig, int *keys)
{


  //decide also the number of threads and block of the grid
  dim3  gridC(nblocksC, 1, 1);  
  dim3  threadsC(ELEMS_BLOCK, nseg_blockC , 1);

  //set the value of the psig array to zero
  cudaMemset((void*) psig,0,8*n1*n2*n3*sizeof(T));

  uncompresscoarsefine<T> <<< gridC, threadsC >>>(n1,n2,n3,psicf,psig,keys);
  cudaThreadSynchronize();

  return 0;

}

template<typename T>
int compressgpu(int n1, int n2, int n3, 
		T *psig,T *psicf, int *keys)
{


  //decide also the number of threads and block of the grid
  dim3  gridC(nblocksC, 1, 1);  
  dim3  threadsC(ELEMS_BLOCK, nseg_blockC , 1);

  compresscoarsefine<T> <<< gridC, threadsC >>>(n1,n2,n3,psig,psicf,keys);
  cudaThreadSynchronize();

  return 0;

}
/****/
