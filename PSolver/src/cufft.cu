#include <iostream>
#include <stdio.h>
#include "cufft.h"
#include "cuda.h"
#include "cuda_runtime_api.h"
#include "config.h"
 
#define DOUBLE

#ifdef DOUBLE
#define Complex  cufftDoubleComplex
#define Real double
#define Transform CUFFT_Z2Z
#define TransformExec cufftExecZ2Z
#else
#define Complex  cufftComplex
#define Real float
#define Transform CUFFT_C2C
#define TransformExec cufftExecC2C
#endif

#define TILE_DIM  8


#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

// synchronize blocks
extern "C" void synchronize_() {
 
  cudaThreadSynchronize();
}

// allocate device memory
extern "C" void FC_FUNC(cudamalloc, CUDAMALLOC) (int *size, Real **d_data,int *ierr) {

  *ierr = cudaMalloc((void**)d_data, sizeof(Real)*(*size));
  if( cudaGetLastError() != cudaSuccess)
      printf("allocate error\n");
}

// allocate device memory
extern "C" void FC_FUNC(cudamemset, CUDAMEMSET) (Real **d_data, int* value, int* size,int *ierr) {

  *ierr = cudaMemset((void*)*d_data, *value, sizeof(Real)*(*size));
  if( cudaGetLastError() != cudaSuccess)
      printf("cudamemset error\n");
}

extern "C" void FC_FUNC(cudafree, CUDAFREE) (Real **d_data) {

  cudaFree(*d_data);
}

extern "C" void FC_FUNC(cufftdestroy, CUFFTDESTROY) (cufftHandle *plan) {

  cufftDestroy(*plan);
}

// set device memory
extern "C" void FC_FUNC_(reset_gpu_data, RESET_GPU_DATA)(int *size, Real* h_data, Real **d_data){

 cudaMemcpy(*d_data, h_data, sizeof(Real)*(*size),
         cudaMemcpyHostToDevice);
 if( cudaGetLastError() != cudaSuccess)
      printf("transfer error\n");

}

// read device memory
extern "C" void FC_FUNC_(get_gpu_data, GET_GPU_DATA)(int *size, Real *h_data, Real **d_data) {

 cudaMemcpy(h_data, *d_data, sizeof(Real)*(*size),
         cudaMemcpyDeviceToHost);
 if (cudaGetLastError() != cudaSuccess)
        printf("transfer back error\n");
}


// transpose
__global__ void transpose(Complex *idata, Complex *odata,
        int width, int height)
{
  __shared__ Complex tile[TILE_DIM][TILE_DIM+1];

  int xIndex = blockIdx.x * TILE_DIM + threadIdx.x;
  int yIndex = blockIdx.y * TILE_DIM + threadIdx.y;
  int index_in = xIndex + (yIndex)*(width);
  int xIndex1 = blockIdx.y * TILE_DIM + threadIdx.x;
  int yIndex1 = blockIdx.x * TILE_DIM + threadIdx.y;
  int index_out = xIndex1 + (yIndex1)*height;

  if (xIndex < width && yIndex < height)
      tile[threadIdx.y][threadIdx.x] = idata[index_in];
    __syncthreads();

  if (xIndex1 < height && yIndex1 < width) {
      odata[index_out] = tile[threadIdx.x][threadIdx.y];
  }
}

// transpose together with spread operation
__global__ void transpose_spread(Complex *idata, Complex *odata, 
	int width, int height, int bign_h)
{
  __shared__ Complex tile[TILE_DIM][TILE_DIM+1];

  int xIndex = blockIdx.x * TILE_DIM + threadIdx.x;
  int yIndex = blockIdx.y * TILE_DIM + threadIdx.y;
  int index_in = xIndex + (yIndex)*(width);
  int xIndex1 = blockIdx.y * TILE_DIM + threadIdx.x;
  int yIndex1 = blockIdx.x * TILE_DIM + threadIdx.y;
  int index_out = xIndex1 + (yIndex1)*height;
  int div = index_out / bign_h;
  int mod = index_out % bign_h;
  index_out = div * (bign_h << 1) + mod+bign_h;
  int plus = -bign_h;

  if (xIndex < width && yIndex < height)
      tile[threadIdx.y][threadIdx.x] = idata[index_in];
    __syncthreads();

  if (xIndex1 < height && yIndex1 < width) {
      odata[index_out] = tile[threadIdx.x][threadIdx.y];
    #ifdef DOUBLE
      odata[index_out + plus] = make_double2(0., 0.);
    #else
      odata[index_out + plus] = make_float2(0.f, 0.f);
    #endif
  }
}

// transpose together with inverse spread operation
__global__ void transpose_spread_i(Complex *idata, Complex *odata,
        int width, int height, int bign_h)
{
  __shared__ Complex tile[TILE_DIM][TILE_DIM+1];

  int xIndex = blockIdx.x * TILE_DIM + threadIdx.x;
  int yIndex = blockIdx.y * TILE_DIM + threadIdx.y;
  int index_in = xIndex + (yIndex)*(width);
  int xIndex1 = blockIdx.y * TILE_DIM + threadIdx.x;
  int yIndex1 = blockIdx.x * TILE_DIM + threadIdx.y;
  int index_out = xIndex1 + (yIndex1)*height;
  int div = index_in / bign_h;
  int mod = index_in % bign_h;
  index_in = div * (bign_h << 1) + mod;

  if (xIndex < width && yIndex < height)
      tile[threadIdx.y][threadIdx.x] = idata[index_in];
    __syncthreads();

  if (xIndex1 < height && yIndex1 < width)
      odata[index_out] = tile[threadIdx.x][threadIdx.y];
}

// spread operation
__global__ void spread(Real* src, unsigned int spitch, Real* dst, unsigned int dpitch)
{
   unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x;
   unsigned int tid = threadIdx.x;
 
   Real res = (tid >= spitch) ? src[bid * spitch + tid-spitch] : 0.0;
   if( tid < dpitch) {
	dst[bid * dpitch + tid] = res;
   }
}

// inverse spread operation
__global__ void spread_i(Real* src, unsigned int spitch, Real* dst, unsigned int dpitch)
{
   unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x;
   unsigned int tid = threadIdx.x;

   Real res = src[bid * dpitch + tid];
   if( tid < dpitch) dst[bid * spitch + tid] = res;
}

// spread operation for 2nd dim
__global__ void spread_y(Complex* src, Complex* dst)
{
   unsigned int tid = (blockIdx.y * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;
   unsigned int tid1 = (blockIdx.y * gridDim.x * 2 + blockIdx.x) * blockDim.x + threadIdx.x;

   Complex res =  src[tid];
   dst[tid1 + blockDim.x*gridDim.x] = res;
#ifdef DOUBLE
   dst[tid1] = make_double2(0., 0.);
#else
   dst[tid1] = make_float2(0.f, 0.f);
#endif
}

__global__ void spread_y_r(Real* src, Real* dst)
{
   unsigned int tid = (blockIdx.y * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;
   unsigned int tid1 = (blockIdx.y * gridDim.x * 2 + blockIdx.x) * blockDim.x + threadIdx.x;

   Real res =  src[tid];
   dst[tid1 + blockDim.x*gridDim.x] = res;
#ifdef DOUBLE
   dst[tid1] = 0.;
#else
   dst[tid1] = 0.f;
#endif
}

__global__ void spread_z(Real* src, Real* dst)
{
   unsigned int tid = (blockIdx.y * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;
   //unsigned int tid1 = (blockIdx.y * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;

   Real res =  src[tid];
   src[tid + (gridDim.y * gridDim.x) * blockDim.x] = res;
#ifdef DOUBLE
   src[tid] = 0.0;
#else
   src[tid] = 0.f;
#endif
}


// inverse spread operation for 2nd dim
__global__ void spread_y_i(Complex* src, Complex* dst)
{
   unsigned int tid = (blockIdx.y * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;
   unsigned int tid1 = (blockIdx.y * gridDim.x * 2 + blockIdx.x) * blockDim.x + threadIdx.x;

   Complex res =  src[tid1];
   dst[tid] = res;
}

__global__ void spread_y_i_r(Real* src, Real* dst)
{
   unsigned int tid = (blockIdx.y * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;
   unsigned int tid1 = (blockIdx.y * gridDim.x * 2 + blockIdx.x) * blockDim.x + threadIdx.x;

   Real res =  src[tid1];
   dst[tid] = res;
}


// multiply with potential
__global__ void multiply_kernel(int nx, int ny, int nz, Complex *d_data, Real *d_kernel, Real scal) {

 int tj = threadIdx.x;
 int td = blockDim.x;

 int blockData = (nx*ny*nz)/(gridDim.x*gridDim.y);

 int jj = (blockIdx.y*gridDim.x + blockIdx.x)*blockData;

 for (int k=0; k<blockData/td; k++) {
     d_data[jj + tj+ k*td].x *= d_kernel[jj + tj+ k*td]*scal;
     d_data[jj + tj+ k*td].y *= d_kernel[jj + tj+ k*td]*scal;
 }

}

// scale
__global__ void scale_kernel(int nx, int ny, int nz, Complex *d_data, Real mult) {

 int tj = threadIdx.x;
 int td = blockDim.x;

 int blockData = (nx*ny*nz)/(gridDim.x*gridDim.y);

 int jj = (blockIdx.y*gridDim.x + blockIdx.x)*blockData;

 for (int k=0; k<blockData/td; k++) {
     d_data[jj + tj+ k*td].x *= mult;
     d_data[jj + tj+ k*td].y *= mult;
 }

}

__global__ void zero(int nx, int ny, int nz, Real *z) {

        int tj = threadIdx.x;
        int td = blockDim.x;

	int blockData = (nx*ny*nz)/(gridDim.x*gridDim.y);

        int jj = ((blockIdx.y)*gridDim.x + (blockIdx.x))*blockData;

        for (int k=0; k<blockData/td; k++) {
        	z[jj + tj+ k*td] = 0.0;
        }
}

__global__ void copy_0(int nx, int ny, int nz, Real *in, Real *out) {

	int tj = threadIdx.x;
        int td = blockDim.x;

        int jj =  (blockIdx.y*nx*ny/4 + blockIdx.x*nx/2);
        int jj1 =  ((blockIdx.y+nz/2)*nx*ny + (blockIdx.x+ny/2)*nx);


        out[jj1+tj+td] = in[jj+tj];

}

__global__ void copy(int nx,int ny,int nz, Real *in, Real *out) {

        int tj = threadIdx.x;
        //int td = blockDim.x;

        int jj =  (blockIdx.y*nx*ny/4 + blockIdx.x*nx/2);
        int jj1 =  ((blockIdx.y)*nx*ny + (blockIdx.x)*nx);

        out[jj+tj] = in[jj1+tj];
}

/************ 1D transform *************/

extern "C" void cuda_1d_plan_(int *NX_p, int *Nbatch_p,
                 cufftHandle *plan) {

 int NX = *NX_p;
 int Nbatch = *Nbatch_p;

 int n1d[3]= {NX, 1, 1};

 if(cufftPlanMany(plan,  1, n1d,
              NULL, 1, NX,
              NULL, 1, NX, Transform, Nbatch) != CUFFT_SUCCESS)
      printf("Error creating plan\n");

 //cufftPlan1d(plan, NX, Transform, Nbatch );

}

extern "C" void cuda_1d_forward_(cufftHandle *plan,
                Complex **d_data, Complex **d_data2) {

   if( TransformExec(*plan, *d_data, *d_data2, CUFFT_FORWARD)!= CUFFT_SUCCESS){
      printf("error in 1D forward transform\n");
   }

}

extern "C" void cuda_1d_inverse_(cufftHandle *plan,
                Complex **d_data, Complex **d_data2) {

   if( TransformExec(*plan, *d_data, *d_data2, CUFFT_INVERSE)!= CUFFT_SUCCESS){
      printf("error in 1D inverse transform\n");
   }

}

/************ 2D transform *************/

extern "C" void cuda_2d_plan_(int *NX_p, int *NY_p, int *Nbatch_p,
                 cufftHandle *plan) {

 int NX = *NX_p;
 int NY = *NY_p;
 int Nbatch = *Nbatch_p;

 int n1d[3]= {NX, NY, 1};

 if(cufftPlanMany(plan,  1, n1d,
              NULL, 1, NX*NY,
              NULL, 1, NX*NY, Transform, Nbatch) != CUFFT_SUCCESS)
      printf("Error creating plan\n");

}

extern "C" void cuda_2d_forward_(cufftHandle *plan,
                Complex **d_data, Complex **d_data2) {

   if( TransformExec(*plan, *d_data, *d_data2, CUFFT_FORWARD)!= CUFFT_SUCCESS){
      printf("error in 2D forward transform\n");
   }

}

extern "C" void cuda_2d_inverse_(cufftHandle *plan,
                Complex **d_data, Complex **d_data2) {

   if( TransformExec(*plan, *d_data, *d_data2, CUFFT_INVERSE)!= CUFFT_SUCCESS){
      printf("error in 2D inverse transform\n");
   }

}

/************ 3D transform *************/

extern "C" void cuda_3d_plan_(int *NX_p, int *NY_p, int *NZ_p,
                 cufftHandle *plan) {

 int NX = *NX_p;
 int NY = *NY_p;
 int NZ = *NZ_p;

 int n[3] = { NZ, NY, NX };
 if(cufftPlanMany(plan, 3, n,
              NULL, 1, NX*NY*NZ,
              NULL, 1, NX*NY*NZ, Transform, 1) != CUFFT_SUCCESS)
      printf("Error creating plan\n");
}

extern "C" void cuda_3d_forward_(cufftHandle *plan,
                Complex **d_data, Complex **d_data2) {

   if( TransformExec(*plan, *d_data, *d_data2, CUFFT_FORWARD)!= CUFFT_SUCCESS){
      printf("error in 3D forward transform\n");
   }

}

extern "C" void cuda_3d_inverse_(int *NX_p, int *NY_p, int *NZ_p ,cufftHandle *plan,
                Complex **d_data, Complex **d_data2) {

   int NX = *NX_p;
   int NY = *NY_p;
   int NZ = *NZ_p;

   if( TransformExec(*plan, *d_data, *d_data2, CUFFT_INVERSE)!= CUFFT_SUCCESS){
      printf("error in 3D inverse transform\n");
   }

   // scale kernel paramters
   int nThreads = NX;
   dim3 nBlocks(NY,NZ,1);

   scale_kernel <<< nBlocks, nThreads >>> (NX,NY,NZ,*d_data2, 1.0/double(NX*NY*NZ));
}

/************ 3D Poisson Solver for periodic boundary *************/

extern "C" void cuda_3d_psolver_cufft3d_plan_(int *NX_p, int *NY_p, int *NZ_p,
                 cufftHandle *plan, cufftHandle *plan1) {

 int NX = *NX_p;
 int NY = *NY_p;
 int NZ = *NZ_p;

 int n[3] = { NZ, NY, NX };
 if(cufftPlanMany(plan, 3, n,
              NULL, 1, NX*NY*NZ,
              NULL, 1, NX*NY*NZ, CUFFT_D2Z, 1) != CUFFT_SUCCESS)
      printf("Error creating plan\n");

 if(cufftPlanMany(plan1, 3, n,
              NULL, 1, NX*NY*NZ,
              NULL, 1, NX*NY*NZ, CUFFT_Z2D, 1) != CUFFT_SUCCESS)
      printf("Error creating plan\n");

}


extern "C" void cuda_3d_psolver_cufft3d_(int *NX_p, int *NY_p, int *NZ_p,cufftHandle *plan,
             cufftHandle *plan1, Complex **d_data, Complex **d_data2, Real **d_kernel, Real *scal_p,
	     int *geo1_p, int *geo2_p, int *geo3_p) {

 int NX = *NX_p;
 int NY = *NY_p;
 int NZ = *NZ_p;

 int geo1 = *geo1_p;
 int geo2 = *geo2_p;
 int geo3 = *geo3_p;

 int ysize = NY/2 + geo2 * NY/2;
 int zsize = NZ/2 + geo3 * NZ/2;

 Real scal = *scal_p;

 // multiply kernel paramters
 int nThreads = NX/2+1;
 dim3 nBlocks(NY,NZ,1);

 // copy kernel paramters
 int nthreads = NX/2;
 dim3 nblocks(NY/2,NZ/2,1);

 // spread kernel parameters
 dim3 nblocks_s(zsize,ysize,1);

 Complex* dst = *d_data;
 Complex* src = *d_data2;

   if (geo1==0 && geo2==0 && geo3==0) {
    src = *d_data;
    dst = *d_data2;
    zero <<< nblocks, nthreads >>> (NX,NY,NZ, (Real*)dst);
    copy_0 <<< nblocks, nthreads  >>> (NX,NY,NZ, (Real*)src, (Real*)dst);
   }
   else {
    if (geo1==0) {
     src = *d_data;
     dst = *d_data2;
     spread<<<nblocks_s, NX>>>((Real*)src, NX/2, (Real*)dst, NX);
    }
    if (geo2==0) {
      if (geo1==0) {
        src = *d_data2;
        dst = *d_data;
      } else {
        src = *d_data;
        dst = *d_data2;
      }
      nblocks_s.x=ysize;
      nblocks_s.y=zsize;
      spread_y_r<<<nblocks_s, NX>>>((Real*)src, (Real*)dst);
    }
    if (geo3==0) {
      nblocks_s.x=NY;
      nblocks_s.y=zsize;
      spread_z<<<nblocks_s, NX>>>((Real*)dst, (Real*)src);
    }
   }

   // Forward FFT

   if( cufftExecD2Z(*plan, (Real*)dst, src)!= CUFFT_SUCCESS){
      printf("error in PSper forward transform\n");
   }

   // multiply with kernel

   multiply_kernel <<< nBlocks, nThreads >>> (NX/2+1,NY,NZ,src,*d_kernel,scal);

   // Inverse FFT

   if( cufftExecZ2D(*plan1, src, (Real*)dst)!= CUFFT_SUCCESS){
      printf("error in PSper inverse transform\n");
   }

   if (geo1==0 && geo2==0 && geo3==0)
     copy <<< nblocks, nthreads >>> (NX,NY,NZ, (Real*)dst, (Real*)src);
   else { 
    if (geo2==0) {
       nblocks_s.x=ysize;
       nblocks_s.y=zsize;
       spread_y_i_r<<<nblocks_s, NX>>>((Real*)dst, (Real*)src);
    }
    if (geo1==0) {
       if (geo2==0) {
        Complex* tmp = src;
        src = dst;
        dst = tmp;
       }
      nblocks_s.x=zsize;
      nblocks_s.y=ysize; 
      spread_i<<<nblocks_s, NX/2>>>((Real*)dst, NX/2, (Real*)src, NX);
    }
   }
}

/************ 3D Poisson Solver for general boundary *************/

extern "C" void FC_FUNC_(cuda_3d_psolver_general_plan, CUDA_3D_PSOLVER_GENERAL_PLAN)(int *N,
                 cufftHandle *plan, int *switch_alg,
		 int *geo) {

 int NX = N[0];
 int NY = N[1];
 int NZ = N[2];

 //int geo1 = geo[0];
 int geo2 = geo[1];
 int geo3 = geo[2];

 int n1d[3]= {1, 1, 1};

 int ysize = NY/2 + geo2 * NY/2;
 int zsize = NZ/2 + geo3 * NZ/2;

 n1d[0] = NX;
 if(cufftPlanMany(plan,  1, n1d,
              NULL, 1, NX,
              NULL, 1, NX, CUFFT_D2Z, ysize*zsize) != CUFFT_SUCCESS)
      printf("Error creating plan\n");

 if(cufftPlanMany(plan+1,  1, n1d,
              NULL, 1, NX,
              NULL, 1, NX, CUFFT_Z2D, ysize*zsize) != CUFFT_SUCCESS)
      printf("Error creating plan\n");

 n1d[0] = NY;
 if(cufftPlanMany(plan+2,  1, n1d,
              NULL, 1, NY,
              NULL, 1, NY, Transform, (NX/2+1)*zsize) != CUFFT_SUCCESS)
      printf("Error creating plan\n");

 n1d[0] = NZ;
 if(cufftPlanMany(plan+3,  1, n1d,
              NULL, 1, NZ,
              NULL, 1, NZ, Transform, (NX/2+1)*NY) != CUFFT_SUCCESS)
      printf("Error creating plan\n");

 *switch_alg = 0;

 /*int nPrimeSize = 17;
 int primeSize[] = {92,104,116,124,136,148,152,164,172,184,188,204,208,220,228,232,248};

 for (int p=0; p<nPrimeSize; p++)
   if (NZ == primeSize[p]) {
     *switch_alg = 1;
     break;
   }

 n1d[0] = NZ;

 int inembed[1];
 int onembed[1];
 inembed[0] = 1;
 onembed[0] = 1;
 if(cufftPlanMany(plan+4,  1, n1d,
              inembed, NY, 1,
              onembed, NY, 1, Transform, NY) != CUFFT_SUCCESS)
      printf("Error creating plan\n");*/

}

extern "C" void FC_FUNC_(cuda_3d_psolver_general, CUDA_3D_PSOLVER_GENERAL)(int *N,
          cufftHandle *plan,
          Complex **d_data, Complex **d_data2, Real **d_kernel, int *switch_alg,
          int *geo, Real *scal_p) {

 int NX = N[0];
 int NY = N[1];
 int NZ = N[2];

 Real scal = *scal_p;

 int geo1 = geo[0];
 int geo2 = geo[1];
 int geo3 = geo[2];

 int ysize=NY/2+geo2*NY/2;
 int zsize=NZ/2+geo3*NZ/2;

 // transpose kernel parameters
 dim3 grid((NX/2+1+TILE_DIM-1)/TILE_DIM,(ysize*zsize+TILE_DIM-1)/TILE_DIM,1);
 dim3 threads(TILE_DIM,TILE_DIM,1);

 // spread kernel parameters
 dim3 nblocks(zsize,ysize,1);

 // multiply kernel paramters
 int nThreads = NX/2+1;
 dim3 nBlocks(NZ,NY,1);

 Complex* dst = *d_data;
 Complex* src = *d_data2;

 // X transform 

   if (geo1==0) {
     src = *d_data;
     dst = *d_data2;
     spread<<<nblocks, NX>>>((Real*)src, NX/2, (Real*)dst, NX);
   }

   if( cufftExecD2Z(plan[0], (Real*)dst, src)!= CUFFT_SUCCESS){
      printf("error in PSolver forward transform 1\n");
   }

   if (geo2==0) {
     transpose_spread <<< grid, threads >>>(src, dst,NX/2+1,ysize*zsize,NY/2);
   } else {
     transpose <<< grid, threads >>>(src, dst,NX/2+1,ysize*zsize);
   }

   // Y transform
   if( TransformExec(plan[2], dst, src, CUFFT_FORWARD)!= CUFFT_SUCCESS){
      printf("error in PSolver forward transform 2\n");
   }

  // Z transform, on entire cube
  if (!(*switch_alg)) {
   grid.x = (NY+TILE_DIM-1)/TILE_DIM;
   grid.y = ((NX/2+1)*zsize+TILE_DIM-1)/TILE_DIM;

   if (geo3==0) {
     transpose_spread <<< grid, threads >>>(src,dst,NY,(NX/2+1)*NZ/2,NZ/2);
   } else {
     transpose <<< grid, threads >>>(src, dst,NY,(NX/2+1)*NZ);
   }

   if( TransformExec(plan[3], dst, src, CUFFT_FORWARD)!= CUFFT_SUCCESS){
      printf("error in PSolver forward transform 3\n");
   }
  }
  else {
   if (geo3==0) {
      nblocks.x=zsize;
      nblocks.y=NX;
      spread_y<<<nblocks, NY>>>(src, dst);
   }

   for(int k=0; k<NX; ++k){
     if( TransformExec(plan[4], dst, src, CUFFT_FORWARD)!= CUFFT_SUCCESS){
      printf("error in PSolver forward transform 3\n");
     }
     src += NY*NZ;
     dst += NY*NZ;
   }

   src -= NX*NY*NZ;
   dst -= NX*NY*NZ;
  }

  // multiply with kernel

  multiply_kernel <<< nBlocks, nThreads >>> (NX/2+1,NY,NZ,src,*d_kernel,scal);

  // inverse transform

  // Z transform, on entire cube 
  if (!(*switch_alg)) {
   if( TransformExec(plan[3], src, dst, CUFFT_INVERSE)!= CUFFT_SUCCESS){
      printf("error in PSolver inverse transform 1\n");
   }

   grid.x = (zsize*(NX/2+1)+TILE_DIM-1)/TILE_DIM;
   grid.y = (NY+TILE_DIM-1)/TILE_DIM;

   if (geo3==0) {
     transpose_spread_i <<< grid, threads >>>(dst, src,NZ/2*(NX/2+1),NY,NZ/2);
   } else {
     transpose <<< grid, threads >>>(dst, src,NZ*(NX/2+1),NY);
   }

  }
  else {

   for(int k=0; k<NX; ++k){
     if( TransformExec(plan[4], src, dst, CUFFT_INVERSE)!= CUFFT_SUCCESS){
      printf("error in PSolver inverse transform 3\n");
     }
     src += NY*NZ;
     dst += NY*NZ;
   }

   src -= NX*NY*NZ;
   dst -= NX*NY*NZ;

   if (geo3==0)
      spread_y_i<<<nblocks, NY>>>(dst, src);
  }

  // Y transform

   if( TransformExec(plan[2], src, dst, CUFFT_INVERSE)!= CUFFT_SUCCESS){
      printf("error in PSolver inverse transform 2\n");
   }

   grid.x = (ysize*zsize+TILE_DIM-1)/TILE_DIM;
   grid.y = (NX/2+1+TILE_DIM-1)/TILE_DIM;

   if (geo2==0) {
      transpose_spread_i <<< grid, threads >>>(dst, src,ysize*zsize,NX/2+1, NY/2);
   } else
      transpose <<< grid, threads >>>(dst, src,ysize*zsize,NX/2+1);

   // X transform

   if( cufftExecZ2D(plan[1], src, (Real*)dst)!= CUFFT_SUCCESS){
      printf("error in PSolver inverse transform 3\n");
   }

   nblocks.x=zsize;
   nblocks.y=ysize;
   if (geo1==0) {
      spread_i<<<nblocks, NX/2>>>((Real*)dst,NX/2, (Real*)src, NX);
   }
}


extern "C" void FC_FUNC_(cuda_3d_psolver_plangeneral, CUDA_3D_PSOLVER_PLANGENERAL)(int *N,
          Complex **d_data, Complex **d_data2, Real **d_kernel,
          int *geo, Real *scal_p) {

 cufftHandle plan;

 int NX = N[0];
 int NY = N[1];
 int NZ = N[2];

 Real scal = *scal_p;

 int geo1 = geo[0];
 int geo2 = geo[1];
 int geo3 = geo[2];

 int ysize=NY/2+geo2*NY/2;
 int zsize=NZ/2+geo3*NZ/2;

 // transpose kernel parameters
 dim3 grid((NX/2+1+TILE_DIM-1)/TILE_DIM,(ysize*zsize+TILE_DIM-1)/TILE_DIM,1);
 dim3 threads(TILE_DIM,TILE_DIM,1);

 // spread kernel parameters
 dim3 nblocks(zsize,ysize,1);

 // multiply kernel paramters
 int nThreads = NX/2+1;
 dim3 nBlocks(NZ,NY,1);

 Complex* dst = *d_data;
 Complex* src = *d_data2;

 int n1d[3]= {1, 1, 1};

 n1d[0] = NX;
 if(cufftPlanMany(&plan,  1, n1d,
              NULL, 1, NX,
              NULL, 1, NX, CUFFT_D2Z, ysize*zsize) != CUFFT_SUCCESS)
      printf("Error creating plan 1\n");

 // X transform 

   if (geo1==0) {
     src = *d_data;
     dst = *d_data2;
     spread<<<nblocks, NX>>>((Real*)src, NX/2, (Real*)dst, NX);
   }

   if( cufftExecD2Z(plan, (Real*)dst, src)!= CUFFT_SUCCESS){
      printf("error in PSolver forward transform 1\n");
   }

   if (geo2==0) {
     transpose_spread <<< grid, threads >>>(src, dst,NX/2+1,ysize*zsize,NY/2);
   } else {
     transpose <<< grid, threads >>>(src, dst,NX/2+1,ysize*zsize);
   }

   cufftDestroy(plan);

   n1d[0] = NY;
   if(cufftPlanMany(&plan,  1, n1d,
              NULL, 1, NY,
              NULL, 1, NY, Transform, (NX/2+1)*zsize) != CUFFT_SUCCESS)
      printf("Error creating plan 2\n");

   // Y transform
   if( TransformExec(plan, dst, src, CUFFT_FORWARD)!= CUFFT_SUCCESS){
      printf("error in PSolver forward transform 2\n");
   }

  // Z transform, on entire cube
   grid.x = (NY+TILE_DIM-1)/TILE_DIM;
   grid.y = ((NX/2+1)*zsize+TILE_DIM-1)/TILE_DIM;

   if (geo3==0) {
     transpose_spread <<< grid, threads >>>(src, dst,NY,(NX/2+1)*NZ/2,NZ/2);
   } else {
     transpose <<< grid, threads >>>(src, dst,NY,(NX/2+1)*NZ);
   }

   cufftDestroy(plan);
   n1d[0] = NZ;
   if(cufftPlanMany(&plan,  1, n1d,
              NULL, 1, NZ,
              NULL, 1, NZ, Transform, (NX/2+1)*NY) != CUFFT_SUCCESS)
      printf("Error creating plan 3\n");

   if( TransformExec(plan, dst, src, CUFFT_FORWARD)!= CUFFT_SUCCESS){
      printf("error in PSolver forward transform 3\n");
   }

  // multiply with kernel

  multiply_kernel <<< nBlocks, nThreads >>> (NX/2+1,NY,NZ,src,*d_kernel,scal);

  // inverse transform

  // Z transform, on entire cube 
   if( TransformExec(plan, src, dst, CUFFT_INVERSE)!= CUFFT_SUCCESS){
      printf("error in PSolver inverse transform 1\n");
   }

   grid.x = (zsize*(NX/2+1)+TILE_DIM-1)/TILE_DIM;
   grid.y = (NY+TILE_DIM-1)/TILE_DIM;

   if (geo3==0) {
     transpose_spread_i <<< grid, threads >>>(dst, src,NZ/2*(NX/2+1),NY,NZ/2);
   } else {
     transpose <<< grid, threads >>>(dst, src,NZ*(NX/2+1),NY);
   }

  // Y transform

   cufftDestroy(plan);
   n1d[0] = NY;
   if(cufftPlanMany(&plan,  1, n1d,
              NULL, 1, NY,
              NULL, 1, NY, Transform, (NX/2+1)*zsize) != CUFFT_SUCCESS)
      printf("Error creating plan 4\n");

   if( TransformExec(plan, src, dst,CUFFT_INVERSE)!= CUFFT_SUCCESS){
      printf("error in PSolver inverse transform 2\n");
   }

   grid.x = (ysize*zsize+TILE_DIM-1)/TILE_DIM;
   grid.y = (NX/2+1+TILE_DIM-1)/TILE_DIM;

   if (geo2==0) {
      transpose_spread_i <<< grid, threads >>>(dst,src,ysize*zsize,NX/2+1, NY/2);
   } else
      transpose <<< grid, threads >>>(dst, src,ysize*zsize,NX/2+1);

   // X transform

   cufftDestroy(plan);
   n1d[0] = NX;
   if(cufftPlanMany(&plan,  1, n1d,
              NULL, 1, NX,
              NULL, 1, NX, CUFFT_Z2D, ysize*zsize) != CUFFT_SUCCESS)
      printf("Error creating plan 5\n");

   if( cufftExecZ2D(plan, src, (Real*)dst)!= CUFFT_SUCCESS){
      printf("error in PSolver inverse transform 3\n");
   }

   nblocks.x=zsize;
   nblocks.y=ysize;
   if (geo1==0) {
      spread_i<<<nblocks, NX/2>>>((Real*)dst,NX/2, (Real*)src, NX);
   }

   cufftDestroy(plan);
}


//Specialization of the computation part for each reduction kernel.
//the kern1_red itself is useless as it is the same for all 3 reductions
//keeping it, as we may want to use another someday

typedef void(*comp_and_red_op)(int, Real*, Real*, Real*, Real*, Real*, Real*, Real*, Real*, Real*, Real*, Real*, Real*, Real*);
typedef void(*red_op)(int, Real*, Real*, Real*, Real*, Real*, Real*, Real*, Real*, Real*, Real*, Real*, Real*, Real*);

__device__
void kern1_comp_and_red (int i , Real* p_GPU, Real* q_GPU, Real* r_GPU, Real* x_GPU, Real* z_GPU, Real* corr_GPU, Real* oneoeps_GPU, Real* alpha_GPU, Real* beta_GPU, Real* beta0_GPU, Real* kappa_GPU, Real* g_odata, Real* sum){
  Real zeta=z_GPU[i]*oneoeps_GPU[i];
  z_GPU[i]=zeta;
  *sum+= (r_GPU[i]*zeta);
}

__device__
void kern1_red (int i , Real* p_GPU, Real* q_GPU, Real* r_GPU, Real* x_GPU, Real* z_GPU, Real* corr_GPU, Real* oneoeps_GPU, Real* alpha_GPU, Real* beta_GPU, Real* beta0_GPU, Real* kappa_GPU, Real* g_odata, Real* sum){
  *sum+= (g_odata[i]);
}


__device__
void kern2_comp_and_red (int i , Real* p_GPU, Real* q_GPU, Real* r_GPU, Real* x_GPU, Real* z_GPU, Real* corr_GPU, Real* oneoeps_GPU, Real* alpha_GPU, Real* beta_GPU, Real* beta0_GPU, Real* kappa_GPU, Real* g_odata, Real* sum){
  Real zeta=z_GPU[i];
  Real pval = zeta+(*beta_GPU / *beta0_GPU)*p_GPU[i];
  Real qval = zeta*corr_GPU[i]+r_GPU[i]+(*beta_GPU / *beta0_GPU)*q_GPU[i];
  p_GPU[i] = pval;
  q_GPU[i] = qval;
  *sum+= (pval*qval);
}

__device__
void kern3_comp_and_red (int i , Real* p_GPU, Real* q_GPU, Real* r_GPU, Real* x_GPU, Real* z_GPU, Real* corr_GPU, Real* oneoeps_GPU, Real* alpha_GPU, Real* beta_GPU, Real* beta0_GPU, Real* kappa_GPU, Real* g_odata, Real* sum){
  x_GPU[i] = x_GPU[i] + *alpha_GPU*p_GPU[i];
  r_GPU[i] = r_GPU[i] - *alpha_GPU*q_GPU[i];
  z_GPU[i] = r_GPU[i] * oneoeps_GPU[i];
  *sum+=(r_GPU[i]*r_GPU[i]);
}




//helper functions for the reduction (reduction taken from NVIDIA cuda samples)
template<class T>
struct SharedMemory
{
    __device__ inline operator       T *()
    {
        extern __shared__ int __smem[];
        return (T *)__smem;
    }

    __device__ inline operator const T *() const
    {
        extern __shared__ int __smem[];
        return (T *)__smem;
    }
};
// specialize for double to avoid unaligned memory
// access compile errors
template<>
struct SharedMemory<double>
{
    __device__ inline operator       double *()
    {
        extern __shared__ double __smem_d[];
        return (double *)__smem_d;
    }

    __device__ inline operator const double *() const
    {
        extern __shared__ double __smem_d[];
        return (double *)__smem_d;
    }
};



/*actual kernel call for the reduction, that is specialized with 2 template 
subkernels, one for computation, the other for specific reduction part.
Result is written in g_odata array in GPU memory. So this must be called several 
times to actually reduce to a single element.
*/
template <unsigned int blockSize, bool nIsPow2, comp_and_red_op op1, red_op op2>
__global__ void
reduce_kernel(int n, int reduceOnly, Real* p_GPU, Real* q_GPU, Real* r_GPU, Real* x_GPU, Real* z_GPU, Real* corr_GPU, Real* oneoeps_GPU, Real* alpha_GPU, Real* beta_GPU, Real* beta0_GPU, Real* kappa_GPU, Real* g_odata)
{
    Real *sdata = SharedMemory<Real>();

    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*blockSize*2 + threadIdx.x;
    unsigned int gridSize = blockSize*2*gridDim.x;

    Real mySum = 0;
    // we reduce multiple elements per thread.  The number is determined by the
    // number of active thread blocks (via gridDim).  More blocks will result
    // in a larger gridSize and therefore fewer elements per thread
    while (i < n)
    {
    if(!reduceOnly){
        op1 (i , p_GPU, q_GPU, r_GPU, x_GPU, z_GPU, corr_GPU, oneoeps_GPU, alpha_GPU, beta_GPU, beta0_GPU, kappa_GPU, g_odata, &mySum);
        // ensure we don't read out of bounds -- this is optimized away for powerOf2 sized arrays
        if (nIsPow2 || i + blockSize < n){
            op1 (i + blockSize, p_GPU, q_GPU, r_GPU, x_GPU, z_GPU, corr_GPU, oneoeps_GPU, alpha_GPU, beta_GPU, beta0_GPU, kappa_GPU, g_odata, &mySum);
        }
    }else{

        //subsequent calls to the kernel after the first one don't have to perform 
        // the computations
        op2 (i, p_GPU, q_GPU, r_GPU, x_GPU, z_GPU, corr_GPU, oneoeps_GPU, alpha_GPU, beta_GPU, beta0_GPU, kappa_GPU, g_odata, &mySum);

        // ensure we don't read out of bounds -- this is optimized away for powerOf2 sized arrays
        if (nIsPow2 || i + blockSize < n)
            op2 (i + blockSize, p_GPU, q_GPU, r_GPU, x_GPU, z_GPU, corr_GPU, oneoeps_GPU, alpha_GPU, beta_GPU, beta0_GPU, kappa_GPU, g_odata, &mySum);
    }

        i += gridSize;
    }

    // each thread puts its local sum into shared memory
    sdata[tid] = mySum;
    __syncthreads();


    // do reduction in shared mem
    if ((blockSize >= 512) && (tid < 256))
    {
        sdata[tid] = mySum = mySum + sdata[tid + 256];
    }

    __syncthreads();

    if ((blockSize >= 256) &&(tid < 128))
    {
            sdata[tid] = mySum = mySum + sdata[tid + 128];
    }

     __syncthreads();

    if ((blockSize >= 128) && (tid <  64))
    {
       sdata[tid] = mySum = mySum + sdata[tid +  64];
    }

    __syncthreads();

#if (__CUDA_ARCH__ >= 300 )
    if ( tid < 32 )
    {
        // Fetch final intermediate sum from 2nd warp
        if (blockSize >=  64) mySum += sdata[tid + 32];
        // Reduce final warp using shuffle
        for (int offset = warpSize/2; offset > 0; offset /= 2) 
        {
            mySum += __shfl_down(mySum, offset);
        }
    }
#else
    // fully unroll reduction within a single warp
    if ((blockSize >=  64) && (tid < 32))
    {
        sdata[tid] = mySum = mySum + sdata[tid + 32];
    }

    __syncthreads();

    if ((blockSize >=  32) && (tid < 16))
    {
        sdata[tid] = mySum = mySum + sdata[tid + 16];
    }

    __syncthreads();

    if ((blockSize >=  16) && (tid <  8))
    {
        sdata[tid] = mySum = mySum + sdata[tid +  8];
    }

    __syncthreads();

    if ((blockSize >=   8) && (tid <  4))
    {
        sdata[tid] = mySum = mySum + sdata[tid +  4];
    }

    __syncthreads();

    if ((blockSize >=   4) && (tid <  2))
    {
        sdata[tid] = mySum = mySum + sdata[tid +  2];
    }

    __syncthreads();

    if ((blockSize >=   2) && ( tid <  1))
    {
        sdata[tid] = mySum = mySum + sdata[tid +  1];
    }

    __syncthreads();
#endif

    // write result for this block to global mem
    if (tid == 0) g_odata[blockIdx.x] = mySum;
}


//wrapper for templated kernel
template <comp_and_red_op op1, red_op op2>
void reduce_step(int s, int threads, int blocks, int reduceOnly,  Real* p_GPU, Real* q_GPU, Real* r_GPU, Real* x_GPU, Real* z_GPU, Real* corr_GPU, Real* oneoeps_GPU, Real* alpha_GPU, Real* beta_GPU, Real* beta0_GPU, Real* kappa_GPU, Real* d_odata){
    //TODO : 2D
    dim3 dimBlock(threads, 1, 1);
    dim3 dimGrid(blocks, 1, 1);
    int smemSize = (threads <= 32) ? 2 * threads * sizeof(Real) : threads * sizeof(Real);

    if (((s&(s-1))==0))//pow2
    {
        switch (threads)
        {
            case 512:
                reduce_kernel<512, true, op1, op2><<< dimGrid, dimBlock, smemSize >>>(s, reduceOnly,  p_GPU, q_GPU, r_GPU, x_GPU, z_GPU, corr_GPU, oneoeps_GPU, alpha_GPU, beta_GPU, beta0_GPU, kappa_GPU, d_odata);
                break;
            case 256:
                reduce_kernel<256, true, op1, op2><<< dimGrid, dimBlock, smemSize >>>(s, reduceOnly,  p_GPU, q_GPU, r_GPU, x_GPU, z_GPU, corr_GPU, oneoeps_GPU, alpha_GPU, beta_GPU, beta0_GPU, kappa_GPU, d_odata);
                break;

            case 128:
                reduce_kernel<128, true, op1, op2><<< dimGrid, dimBlock, smemSize >>>(s, reduceOnly,  p_GPU, q_GPU, r_GPU, x_GPU, z_GPU, corr_GPU, oneoeps_GPU, alpha_GPU, beta_GPU, beta0_GPU, kappa_GPU, d_odata);
                break;
            case 64:
                reduce_kernel<64, true, op1, op2><<< dimGrid, dimBlock, smemSize >>>(s, reduceOnly,  p_GPU, q_GPU, r_GPU, x_GPU, z_GPU, corr_GPU, oneoeps_GPU, alpha_GPU, beta_GPU, beta0_GPU, kappa_GPU, d_odata);
                break;
            case 32:
                reduce_kernel<32, true, op1, op2><<< dimGrid, dimBlock, smemSize >>>(s, reduceOnly,  p_GPU, q_GPU, r_GPU, x_GPU, z_GPU, corr_GPU, oneoeps_GPU, alpha_GPU, beta_GPU, beta0_GPU, kappa_GPU, d_odata);
                break;
            case 16:
                reduce_kernel<16, true, op1, op2><<< dimGrid, dimBlock, smemSize >>>(s, reduceOnly,  p_GPU, q_GPU, r_GPU, x_GPU, z_GPU, corr_GPU, oneoeps_GPU, alpha_GPU, beta_GPU, beta0_GPU, kappa_GPU, d_odata);
                break;
            case  8:
                reduce_kernel<8, true, op1, op2><<< dimGrid, dimBlock, smemSize >>>(s, reduceOnly,  p_GPU, q_GPU, r_GPU, x_GPU, z_GPU, corr_GPU, oneoeps_GPU, alpha_GPU, beta_GPU, beta0_GPU, kappa_GPU, d_odata);
                break;
            case  4:
                reduce_kernel<4, true, op1, op2><<< dimGrid, dimBlock, smemSize >>>(s, reduceOnly,  p_GPU, q_GPU, r_GPU, x_GPU, z_GPU, corr_GPU, oneoeps_GPU, alpha_GPU, beta_GPU, beta0_GPU, kappa_GPU, d_odata);
                break;
            case  2:
                reduce_kernel<2, true, op1, op2><<< dimGrid, dimBlock, smemSize >>>(s, reduceOnly,  p_GPU, q_GPU, r_GPU, x_GPU, z_GPU, corr_GPU, oneoeps_GPU, alpha_GPU, beta_GPU, beta0_GPU, kappa_GPU, d_odata);
                break;
            case  1:
                reduce_kernel<1, true, op1, op2><<< dimGrid, dimBlock, smemSize >>>(s, reduceOnly,  p_GPU, q_GPU, r_GPU, x_GPU, z_GPU, corr_GPU, oneoeps_GPU, alpha_GPU, beta_GPU, beta0_GPU, kappa_GPU, d_odata);
                break;
        }
    }
    else
    {
        switch (threads)
        {
            case 512:
                reduce_kernel<512, false, op1, op2><<< dimGrid, dimBlock, smemSize >>>(s, reduceOnly,  p_GPU, q_GPU, r_GPU, x_GPU, z_GPU, corr_GPU, oneoeps_GPU, alpha_GPU, beta_GPU, beta0_GPU, kappa_GPU, d_odata);
                break;
            case 256:
                reduce_kernel<256, false, op1, op2><<< dimGrid, dimBlock, smemSize >>>(s, reduceOnly,  p_GPU, q_GPU, r_GPU, x_GPU, z_GPU, corr_GPU, oneoeps_GPU, alpha_GPU, beta_GPU, beta0_GPU, kappa_GPU, d_odata);
                break;
            case 128:
                reduce_kernel<128, false, op1, op2><<< dimGrid, dimBlock, smemSize >>>(s, reduceOnly,  p_GPU, q_GPU, r_GPU, x_GPU, z_GPU, corr_GPU, oneoeps_GPU, alpha_GPU, beta_GPU, beta0_GPU, kappa_GPU, d_odata);
                break;
            case 64:
                reduce_kernel<64, false, op1, op2><<< dimGrid, dimBlock, smemSize >>>(s, reduceOnly,  p_GPU, q_GPU, r_GPU, x_GPU, z_GPU, corr_GPU, oneoeps_GPU, alpha_GPU, beta_GPU, beta0_GPU, kappa_GPU, d_odata);
                break;
            case 32:
                reduce_kernel<32, false, op1, op2><<< dimGrid, dimBlock, smemSize >>>(s, reduceOnly,  p_GPU, q_GPU, r_GPU, x_GPU, z_GPU, corr_GPU, oneoeps_GPU, alpha_GPU, beta_GPU, beta0_GPU, kappa_GPU, d_odata);
                break;
            case 16:
                reduce_kernel<16, false, op1, op2><<< dimGrid, dimBlock, smemSize >>>(s, reduceOnly,  p_GPU, q_GPU, r_GPU, x_GPU, z_GPU, corr_GPU, oneoeps_GPU, alpha_GPU, beta_GPU, beta0_GPU, kappa_GPU, d_odata);
                break;
            case  8:
                reduce_kernel<8, false, op1, op2><<< dimGrid, dimBlock, smemSize >>>(s, reduceOnly,  p_GPU, q_GPU, r_GPU, x_GPU, z_GPU, corr_GPU, oneoeps_GPU, alpha_GPU, beta_GPU, beta0_GPU, kappa_GPU, d_odata);
                break;
            case  4:
                reduce_kernel<4, false, op1, op2><<< dimGrid, dimBlock, smemSize >>>(s, reduceOnly,  p_GPU, q_GPU, r_GPU, x_GPU, z_GPU, corr_GPU, oneoeps_GPU, alpha_GPU, beta_GPU, beta0_GPU, kappa_GPU, d_odata);
                break;
            case  2:
                reduce_kernel<2, false, op1, op2><<< dimGrid, dimBlock, smemSize >>>(s, reduceOnly,  p_GPU, q_GPU, r_GPU, x_GPU, z_GPU, corr_GPU, oneoeps_GPU, alpha_GPU, beta_GPU, beta0_GPU, kappa_GPU, d_odata);
                break;
            case  1:
                reduce_kernel<1, false, op1, op2><<< dimGrid, dimBlock, smemSize >>>(s, reduceOnly,  p_GPU, q_GPU, r_GPU, x_GPU, z_GPU, corr_GPU, oneoeps_GPU, alpha_GPU, beta_GPU, beta0_GPU, kappa_GPU, d_odata);
                break;
        }
    }

gpuErrchk( cudaPeekAtLastError() );
gpuErrchk( cudaDeviceSynchronize() );

//  cudaError_t error = cudaGetLastError();
//  if(error != cudaSuccess)
//  {
//    // print the CUDA error message and exit
//    printf("CUDA error in kernel : %s", cudaGetErrorString(error));
//exit(0);
//  }
}

unsigned int nextPow2(unsigned int x)
{
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return ++x;
}

/*this performs some calculations to chose the size of the blocks we want to use 
for reduction, while limiting their number for efficiency purposes, as each kernel
 will handle several elements in this version (see last version of reduction in
reduction sample code from Nvidia)
*/
template <comp_and_red_op op1, red_op op2>
void apply_reduction(int n,
          Real* p_GPU, Real* q_GPU, Real* r_GPU, Real* x_GPU, Real* z_GPU, Real* corr_GPU, Real* oneoeps_GPU, Real* alpha_GPU, Real* beta_GPU, Real* beta0_GPU, Real* kappa_GPU, Real* result) {
    int maxThreads=256;
    int maxBlocks=64;
    int blocks=0;
    int threads=0;

    //get device capability, to avoid block/grid size excceed the upbound
    cudaDeviceProp prop;
    int device;
    cudaGetDevice(&device);
    cudaGetDeviceProperties(&prop, device);

    threads = (n < maxThreads*2) ? nextPow2((n + 1)/ 2) : maxThreads;
    blocks = (n + (threads * 2 - 1)) / (threads * 2);

    if ((Real)threads*blocks > (Real)prop.maxGridSize[0] * prop.maxThreadsPerBlock)
    {
        printf("n is too large, please choose a smaller number!\n");
    }

    if (blocks > prop.maxGridSize[0])
    {
        printf("Grid size <%d> excceeds the device capability <%d>, set block size as %d (original %d)\n",
               blocks, prop.maxGridSize[0], threads*2, threads);

        blocks /= 2;
        threads *= 2;
    }

    //we will only use maxblocks blocks, and make each thread work on more data
    blocks = min(maxBlocks, blocks);


    Real *d_odata = NULL;
    cudaMalloc((void **) &d_odata, blocks*sizeof(Real));
  if( cudaGetLastError() != cudaSuccess)
      printf("allocate error\n");
    //first reduction
    cudaDeviceSynchronize();
    reduce_step<op1, op2>(n, threads, blocks, 0,  p_GPU, q_GPU, r_GPU, x_GPU, z_GPU, corr_GPU, oneoeps_GPU, alpha_GPU, beta_GPU, beta0_GPU, kappa_GPU, d_odata);
  
  if( cudaGetLastError() != cudaSuccess)
    printf("result during first reduction kernel\n");

    int s=blocks;
    //loop and perform as many reductions steps as necessary
    while (s > 1)
    {
        threads = (s < maxThreads*2) ? nextPow2((s + 1)/ 2) : maxThreads;
        blocks = (s + (threads * 2 - 1)) / (threads * 2);
        if (blocks > prop.maxGridSize[0])
        {
            printf("Grid size <%d> excceeds the device capability <%d>, set block size as %d (original %d)\n",
            blocks, prop.maxGridSize[0], threads*2, threads);

            blocks /= 2;
            threads *= 2;
        }
        blocks = min(maxBlocks, blocks);

        reduce_step<op1, op2>(s, threads, blocks, 1, p_GPU, q_GPU, r_GPU, x_GPU, z_GPU, corr_GPU, oneoeps_GPU, alpha_GPU, beta_GPU, beta0_GPU, kappa_GPU, d_odata);
        if( cudaGetLastError() != cudaSuccess)
            printf("result during reduction kernel\n");
        s = (s + (threads*2-1)) / (threads*2);
    }

  //TODO: move result copy to user code ?
  cudaMemcpy(result, d_odata, sizeof(Real), cudaMemcpyDeviceToHost);
  if( cudaGetLastError() != cudaSuccess)
    printf("result retrieve error\n");
  cudaFree(d_odata);
  cudaDeviceSynchronize();
}

//these will be called from fortran, and apply the reduction with the right subkernels

extern "C" void FC_FUNC_(first_reduction_kernel, FIRST_REDUCTION_KERNEL)(int* n1, int* n23,
          Real** p_GPU, Real** q_GPU, Real** r_GPU, Real** x_GPU, Real** z_GPU, Real** corr_GPU, Real** oneoeps_GPU, Real** alpha_GPU, Real** beta_GPU, Real** beta0_GPU, Real** kappa_GPU, Real* result) {

    int n=(*n1) * (*n23);
    apply_reduction<kern1_comp_and_red, kern1_red>(n, *p_GPU, *q_GPU, *r_GPU, *x_GPU, *z_GPU, *corr_GPU, *oneoeps_GPU, *alpha_GPU, *beta_GPU, *beta0_GPU, *kappa_GPU, result);

}

extern "C" void FC_FUNC_(second_reduction_kernel, SECOND_REDUCTION_KERNEL)(int* n1, int* n23,
          Real** p_GPU, Real** q_GPU, Real** r_GPU, Real** x_GPU, Real** z_GPU, Real** corr_GPU, Real** oneoeps_GPU, Real** alpha_GPU, Real** beta_GPU, Real** beta0_GPU, Real** kappa_GPU, Real* result) {

    int n=(*n1) * (*n23);
    apply_reduction<kern2_comp_and_red, kern1_red>(n, *p_GPU, *q_GPU, *r_GPU, *x_GPU, *z_GPU, *corr_GPU, *oneoeps_GPU, *alpha_GPU, *beta_GPU, *beta0_GPU, *kappa_GPU, result);

}

extern "C" void FC_FUNC_(third_reduction_kernel, THIRD_REDUCTION_KERNEL)(int* n1, int* n23,
          Real** p_GPU, Real** q_GPU, Real** r_GPU, Real** x_GPU, Real** z_GPU, Real** corr_GPU, Real** oneoeps_GPU, Real** alpha_GPU, Real** beta_GPU, Real** beta0_GPU, Real** kappa_GPU, Real* result) {

    int n=(*n1) * (*n23);
    apply_reduction<kern3_comp_and_red, kern1_red>(n, *p_GPU, *q_GPU, *r_GPU, *x_GPU, *z_GPU, *corr_GPU, *oneoeps_GPU, *alpha_GPU, *beta_GPU, *beta0_GPU, *kappa_GPU, result);

}
