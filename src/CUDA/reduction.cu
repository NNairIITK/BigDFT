/****c* CUDA/reduction.h
** 
** AUTHOR
**  Luigi Genovese
** 
** CHANGELOG
**
** SOURCE
*/

#include <stdio.h>

#include "reduction.hcu"

//values from which the reduction is done on CPU
#define  COPYVALUES 4096




template<typename T>
float reducearrays(int n,
		   int ndat,
		   T *psi,
		   T *vpsi,
		   T *epot)
{

  float epots[COPYVALUES];

  //first reduce the two arrays in the input-output form
  int threads=256; //first value, fixed
  int blocks=(n*ndat-1)/256 + 1; //ntot-1/256 +1

  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  //allocate on the GPU the memory to host the reduction arrays
  T* wrkred[2];
  size_t size0=blocks*sizeof(T);
  size_t size1=(blocks/2+1)*sizeof(T);

  if(cudaMalloc( (void**) &(wrkred[0]), size0) != 0)
    {
      printf("reducearrays:GPU allocation error 1 \n");
      return 0.f;
    }
  if(cudaMalloc( (void**) &(wrkred[1]), size1) != 0)
    {
      printf("reducearrays:GPU allocation error 2\n");
      return 0.f;
    }

  switch (threads)
    {
    case 512:
      reducefirst<T, 512><<< dimGrid, dimBlock >>>(n,ndat,psi,vpsi,wrkred[0]); break;
    case 256:	 
      reducefirst<T, 256><<< dimGrid, dimBlock >>>(n,ndat,psi,vpsi,wrkred[0]); break;
    case 128:	 
      reducefirst<T, 128><<< dimGrid, dimBlock >>>(n,ndat,psi,vpsi,wrkred[0]); break;
    case 64:	 
      reducefirst<T, 64><<< dimGrid, dimBlock >>>(n,ndat,psi,vpsi,wrkred[0]); break;
    case 32:	 
      reducefirst<T, 32><<< dimGrid, dimBlock >>>(n,ndat,psi,vpsi,wrkred[0]); break;
    case 16:	 
      reducefirst<T, 16><<< dimGrid, dimBlock >>>(n,ndat,psi,vpsi,wrkred[0]); break;
    case 8:	 
      reducefirst<T,  8><<< dimGrid, dimBlock >>>(n,ndat,psi,vpsi,wrkred[0]); break;
    case 4:	 
      reducefirst<T,  4><<< dimGrid, dimBlock >>>(n,ndat,psi,vpsi,wrkred[0]); break;
    case 2:	 
      reducefirst<T,  2><<< dimGrid, dimBlock >>>(n,ndat,psi,vpsi,wrkred[0]); break;
    case 1:	 
      reducefirst<T,  1><<< dimGrid, dimBlock >>>(n,ndat,psi,vpsi,wrkred[0]); break;
    default:
      exit(1);
    }

  int ntot=blocks;
  threads=512;
  int iin=0;
  int iout=1;

  //then pass to the reduction case until only one element is left

  while (blocks > COPYVALUES)
    {
      while(ntot < 2*threads){threads/=2;}

      blocks=(ntot-1)/(2*threads) +1;

      dim3 dimBlock(threads, 1, 1);
      dim3 dimGrid(blocks, 1, 1);

      switch (threads)
	{
	case 512:
	  reducethen<T, 512><<< dimGrid, dimBlock >>>(ntot,wrkred[iin],wrkred[iout]); break;
	case 256:
	  reducethen<T, 256><<< dimGrid, dimBlock >>>(ntot,wrkred[iin],wrkred[iout]); break;
	case 128:
	  reducethen<T, 128><<< dimGrid, dimBlock >>>(ntot,wrkred[iin],wrkred[iout]); break;
	case 64:
	  reducethen<T, 64><<< dimGrid, dimBlock >>>(ntot,wrkred[iin],wrkred[iout]); break;
	case 32:
	  reducethen<T, 32><<< dimGrid, dimBlock >>>(ntot,wrkred[iin],wrkred[iout]); break;
	case 16:
	  reducethen<T, 16><<< dimGrid, dimBlock >>>(ntot,wrkred[iin],wrkred[iout]); break;
	case 8:
	  reducethen<T,  8><<< dimGrid, dimBlock >>>(ntot,wrkred[iin],wrkred[iout]); break;
	case 4:
	  reducethen<T,  4><<< dimGrid, dimBlock >>>(ntot,wrkred[iin],wrkred[iout]); break;
	case 2:
	  reducethen<T,  2><<< dimGrid, dimBlock >>>(ntot,wrkred[iin],wrkred[iout]); break;
	case 1:
	  reducethen<T,  1><<< dimGrid, dimBlock >>>(ntot,wrkred[iin],wrkred[iout]); break;
	default:
	  exit(1);
	}

      ntot=blocks;
      iin=1-iin;
      iout=1-iout;
      //printf("ntot,blocks,iin,iout %i %i %i %i\n",ntot,blocks,iin,iout); 

    }

  //printf("ntot,blocks,iin,iout %i %i %i %i\n",ntot,blocks,iin,iout); 
  cudaFree(wrkred[iout]);


  if(cudaMemcpy(epots,wrkred[iin], blocks*sizeof(T), cudaMemcpyDeviceToHost)  != 0)
    {
      printf("reducearrays: DeviceToHost Memcpy error \n");
      return 0.f;
    }


  //cudaFree(wrkred);
  cudaFree(wrkred[iin]);

  
  register T  result=epots[0];

    for( int i = 1 ; i < blocks; i++ )  
      {
      result += epots[i];
      //printf ("%f",epots);
    }
  
    *epot=result;

  return result;
  
}


/****/

//little hack : specification of double and float reducearray in order to compil only once the template
//Of course, this fct is never called
void fctHackTemplate()
{
  float yo = reducearrays<float>(NULL,NULL,NULL,NULL,NULL);
  float yo2 = reducearrays<double>(NULL,NULL,NULL,NULL,NULL);
}
