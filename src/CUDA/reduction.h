/****c* CUDA/reduction.h
** 
** AUTHOR
**  Luigi Genovese
** 
** CHANGELOG
**
** SOURCE
*/

#ifndef   	REDUCTION_H_
# define   	REDUCTION_H_

//values from which the reduction is done on CPU
#define  COPYVALUES 1024



float reducearrays(int n,
		   int ndat,
		   float *psi,
		   float *vpsi);


template <unsigned int blockSize>
__global__ void reducefirst(int n,int ndat, float *psi, float *vpsi, float *psivpsi)
{

  const unsigned int tid = threadIdx.x; //from 1 to 256
  const unsigned int ivp = blockIdx.x*blockSize + tid; //index of vpsi array
						 //(coalesced)

  __shared__ float pvp_sh[blockSize];
  
  //index of psi array (uncoalesced), may benefit from textures
  const unsigned int ix=ivp/n;
  const unsigned int iy=ivp - n*ix;
  const unsigned int ip=ix + ndat*iy;

  //copy the two arrays in shared memory
  pvp_sh[tid]=0.f;
  
  if(ivp < n*ndat)
    {
      pvp_sh[tid]=psi[ip]*vpsi[ivp];
    }

  //end shared memory copy
  __syncthreads();

  //reduction of the array in one element
  if (blockSize >= 512) { if (tid < 256) { pvp_sh[tid] += pvp_sh[tid + 256]; } __syncthreads(); }
  if (blockSize >= 256) { if (tid < 128) { pvp_sh[tid] += pvp_sh[tid + 128]; } __syncthreads(); }
  if (blockSize >= 128) { if (tid < 64) { pvp_sh[tid] += pvp_sh[tid + 64]; } __syncthreads(); }

  if (tid < 32) {
    if (blockSize >= 64) pvp_sh[tid] += pvp_sh[tid + 32];
    if (blockSize >= 32) pvp_sh[tid] += pvp_sh[tid + 16];
    if (blockSize >= 16) pvp_sh[tid] += pvp_sh[tid + 8];
    if (blockSize >= 8) pvp_sh[tid] += pvp_sh[tid + 4];
    if (blockSize >= 4) pvp_sh[tid] += pvp_sh[tid + 2];
    if (blockSize >= 2) pvp_sh[tid] += pvp_sh[tid + 1];
  }
  if (tid == 0) psivpsi[blockIdx.x] = pvp_sh[0];
}


template <unsigned int blockSize>
__global__ void reducethen(int n, float *in,float *out)
{

  const unsigned int tid = threadIdx.x; //from 1 to 256
  unsigned int i = blockIdx.x*blockSize*2 + tid; //(coalesced)
  const unsigned int gridSize = blockSize*2*gridDim.x;

  __shared__ float pvp_sh[blockSize];
  
  pvp_sh[tid]=0.f;
  
  while(i < n)
    {
      pvp_sh[tid]+=in[i]+in[i+blockSize];
      i+=gridSize;
    }

  //end shared memory copy
  __syncthreads();

  //reduction of the array in one element
  if (blockSize >= 512) { if (tid < 256) { pvp_sh[tid] += pvp_sh[tid + 256]; } __syncthreads(); }
  if (blockSize >= 256) { if (tid < 128) { pvp_sh[tid] += pvp_sh[tid + 128]; } __syncthreads(); }
  if (blockSize >= 128) { if (tid < 64) { pvp_sh[tid] += pvp_sh[tid + 64]; } __syncthreads(); }

  if (tid < 32) {
    if (blockSize >= 64) pvp_sh[tid] += pvp_sh[tid + 32];
    if (blockSize >= 32) pvp_sh[tid] += pvp_sh[tid + 16];
    if (blockSize >= 16) pvp_sh[tid] += pvp_sh[tid + 8];
    if (blockSize >= 8) pvp_sh[tid] += pvp_sh[tid + 4];
    if (blockSize >= 4) pvp_sh[tid] += pvp_sh[tid + 2];
    if (blockSize >= 2) pvp_sh[tid] += pvp_sh[tid + 1];
  }
  if (tid == 0) out[blockIdx.x] = pvp_sh[0];
}

float reducearrays(int n,
		   int ndat,
		   float *psi,
		   float *vpsi)
{

  float epot[COPYVALUES];

  //first reduce the two arrays in the input-output form
  int threads=256; //first value, fixed
  int blocks=(n*ndat-1) >> 8 + 1; //ntot-1/256 +1

  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  //allocate on the GPU the memory to host the reduction arrays
  float* wrkred[2];
  size_t size0=blocks*sizeof(float);
  size_t size1=(blocks/2+1)*sizeof(float);

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
      reducefirst<512><<< dimGrid, dimBlock >>>(n,ndat,psi,vpsi,wrkred[0]); break;
    case 256:
      reducefirst<256><<< dimGrid, dimBlock >>>(n,ndat,psi,vpsi,wrkred[0]); break;
    case 128:
      reducefirst<128><<< dimGrid, dimBlock >>>(n,ndat,psi,vpsi,wrkred[0]); break;
    case 64:
      reducefirst< 64><<< dimGrid, dimBlock >>>(n,ndat,psi,vpsi,wrkred[0]); break;
    case 32:
      reducefirst< 32><<< dimGrid, dimBlock >>>(n,ndat,psi,vpsi,wrkred[0]); break;
    case 16:
      reducefirst< 16><<< dimGrid, dimBlock >>>(n,ndat,psi,vpsi,wrkred[0]); break;
    case 8:
      reducefirst<  8><<< dimGrid, dimBlock >>>(n,ndat,psi,vpsi,wrkred[0]); break;
    case 4:
      reducefirst<  4><<< dimGrid, dimBlock >>>(n,ndat,psi,vpsi,wrkred[0]); break;
    case 2:
      reducefirst<  2><<< dimGrid, dimBlock >>>(n,ndat,psi,vpsi,wrkred[0]); break;
    case 1:
      reducefirst<  1><<< dimGrid, dimBlock >>>(n,ndat,psi,vpsi,wrkred[0]); break;
    default:
      exit(1);
    }

  int ntot=n*ndat;
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
	  reducethen<512><<< dimGrid, dimBlock >>>(ntot,wrkred[iin],wrkred[iout]); break;
	case 256:
	  reducethen<256><<< dimGrid, dimBlock >>>(ntot,wrkred[iin],wrkred[iout]); break;
	case 128:
	  reducethen<128><<< dimGrid, dimBlock >>>(ntot,wrkred[iin],wrkred[iout]); break;
	case 64:
	  reducethen< 64><<< dimGrid, dimBlock >>>(ntot,wrkred[iin],wrkred[iout]); break;
	case 32:
	  reducethen< 32><<< dimGrid, dimBlock >>>(ntot,wrkred[iin],wrkred[iout]); break;
	case 16:
	  reducethen< 16><<< dimGrid, dimBlock >>>(ntot,wrkred[iin],wrkred[iout]); break;
	case 8:
	  reducethen<  8><<< dimGrid, dimBlock >>>(ntot,wrkred[iin],wrkred[iout]); break;
	case 4:
	  reducethen<  4><<< dimGrid, dimBlock >>>(ntot,wrkred[iin],wrkred[iout]); break;
	case 2:
	  reducethen<  2><<< dimGrid, dimBlock >>>(ntot,wrkred[iin],wrkred[iout]); break;
	case 1:
	  reducethen<  1><<< dimGrid, dimBlock >>>(ntot,wrkred[iin],wrkred[iout]); break;
	default:
	  exit(1);
	}

      ntot=blocks;
      iin=1-iin;
      iout=1-iout;

      printf("ntot,blocks,iin,iout %i %i %i %i\n",ntot,blocks,iin,iout); 

    }

  if(cudaMemcpy(&epot,&(wrkred[iin]), blocks*sizeof(float), cudaMemcpyDeviceToHost)  != 0)
    {
      printf("reducearrays: DeviceToHost Memcpy error \n");
      return 0.f;
    }

  cudaFree(wrkred[1]);
  cudaFree(wrkred[0]);

  register float result=epot[0];

    for( int i = 1 ; i < blocks; i++ ) 
      {
      result += epot[i];
    }
  
  return result;

}
#endif 	    /* !REDUCTION_H_ */

/****/
