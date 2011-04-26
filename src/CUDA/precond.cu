//! @file
//! Cuda
//!
//! @author
//!    Copyright (C) 2009-2011 BigDFT group 
//!    This file is distributed under the terms of the
//!    GNU General Public License, see ~/COPYING file
//!    or http://www.gnu.org/copyleft/gpl.txt .
//!    For the list of contributors, see ~/AUTHORS 


#define CUERR { cudaError_t err; \
 if ((err = cudaGetLastError()) != cudaSuccess) { \
 printf("CUDA error: %s, line %d\n", cudaGetErrorString(err), __LINE__); }}



template<typename T>
int gpuapply_hp(int n1,int n2, int n3,
		T h1,T h2,T h3,T c,
		T *in,T *out,int *keys,
		T *work1,T *work2,T *cpsifscf)
{

  
  dim3  gridWT3(linecutsWTz,  numBlocksWTz, 1);  
  dim3  threadsWT3(HALF_WARP_SIZE, num_halfwarpsWTz , 1);
  dim3  gridWT2(linecutsWTy,  numBlocksWTy, 1);  
  dim3  threadsWT2(HALF_WARP_SIZE, num_halfwarpsWTy , 1);
  dim3  gridWT1(linecutsWTx,  numBlocksWTx, 1);  
  dim3  threadsWT1(HALF_WARP_SIZE, num_halfwarpsWTx , 1);

  dim3  gridK3(linecutsKz,  numBlocksKz, 1);  
  dim3  threadsK3(HALF_WARP_SIZE, num_halfwarpsKz , 1);
  dim3  gridK2(linecutsKy,  numBlocksKy, 1);  
  dim3  threadsK2(HALF_WARP_SIZE, num_halfwarpsKy , 1);
  dim3  gridK1(linecutsKx,  numBlocksKx, 1);  
  dim3  threadsK1(HALF_WARP_SIZE, num_halfwarpsKx , 1);

  dim3  gridC(nblocksC, 1, 1);  
  dim3  threadsC(ELEMS_BLOCK, nseg_blockC , 1);

  //uncompression
  //set the value of the psig array to zero
  //do_blasScal(8*n1*n2*n3,(T) 0.,psi); //alternative way
  cudaMemset((void*) work1,0,8*n1*n2*n3*sizeof(T));
  cudaThreadSynchronize();

  //decompress wavefunctions with scaling factor
  uncompresscoarsefinescal<T> <<< gridC, threadsC >>>(n1,n2,n3,h1,h2,h3,c,in,work1,keys);
  cudaThreadSynchronize();
  

  //wavelet transformation
  waveletstosf<T>(n1,n2,n3,gridWT1,gridWT2,gridWT3,
		  threadsWT1,threadsWT2,threadsWT3,
		  work1,work2);

  //worky array should be initialised to c*x
  c_initialize<T> <<< gridK3, threadsK3 >>>(2*n3,4*n1*n2,work2,cpsifscf,c,3);
  cudaThreadSynchronize();

  kineticapplication<T>(n1,n2,n3,h1,h2,h3,
			gridK1,gridK2,gridK3,
			threadsK1,threadsK2,threadsK3,
			work2,cpsifscf,work1,
			out); //work array

  sftowavelets<T>(n1,n2,n3,gridWT1,gridWT2,gridWT3,
		  threadsWT1,threadsWT2,threadsWT3,
		  work1,work2);


  //recompress
  compresscoarsefinescal<T> <<< gridC, threadsC >>>(n1,n2,n3,h1,h2,h3,c,work2,out,keys);
  cudaThreadSynchronize();

  return 0;

}

//preconditioning loop
//reproduces the precong_per routine without prec_fft_c
template<typename T>
int gpucg_precong(int n1,int n2, int n3,int npsi,int ncong,
		  T h1,T h2,T h3,T c,T *x,int *keys,T *r,T *b,T *d,
		  T *work1,T *work2,T *work3, T *gnrm)
{

  dim3  gridC(nblocksC, 1, 1);  
  dim3  threadsC(ELEMS_BLOCK, nseg_blockC , 1);


  //calculate the residue for the orbital
  *gnrm=do_blasDot(npsi, x, x);
  cudaThreadSynchronize();

  wscalgpu<T> <<< gridC, threadsC >>>(x,h1,h2,h3,c,keys);
  cudaThreadSynchronize();

  gpuapply_hp<T>(n1,n2,n3,h1,h2,h3,c,
		 x,d,keys,
		 work1,work2,work3);

  //change the sign of the 0-th step such as to use axpy calls
  //d=d-x
  do_blasAxpy(npsi,(T)(-1.),x,d);
  cudaThreadSynchronize();

  //r=d
  cudaMemcpy(r,d,npsi*sizeof(T), cudaMemcpyDeviceToDevice);
  cudaThreadSynchronize();

  T rmr_new=do_blasDot(npsi, r, r);
  cudaThreadSynchronize();

  T alpha,beta,rmr_old;

  for(int i=0;i < ncong;++i)
    {
      gpuapply_hp<T>(n1,n2,n3,h1,h2,h3,c,
		     d,b,keys,
		     work1,work2,work3);

      alpha=rmr_new/do_blasDot(npsi, d, b);
      cudaThreadSynchronize();

      //here the sign is inverted because of the
      //mapping d -> -d => b -> -b
      
      //x=x-alpha*d
      do_blasAxpy(npsi,-alpha,d,x);
      cudaThreadSynchronize();

      if (i != ncong-1)
	{
	  //r=r-alpha*b (r does not change sign since also b is opposite)
	  do_blasAxpy(npsi,-alpha,b,r);
	  cudaThreadSynchronize();

	  rmr_old=rmr_new;
	  rmr_new=do_blasDot(npsi, r, r);
	  cudaThreadSynchronize();

	  beta=rmr_new/rmr_old;

	  //d=d+1/beta*r
	  do_blasAxpy(npsi,(T)(1)/beta,r,d);
	  cudaThreadSynchronize();

	  //d=beta*d
	  do_blasScal(npsi,beta,d);
	  cudaThreadSynchronize();
	}
  
	
    }


  wscalgpu<T> <<< gridC, threadsC >>>(x,h1,h2,h3,c,keys);
  cudaThreadSynchronize();


  return 0;

}

//preconditioning loop
//reproduces the precong_per after the precondition_preconditioner
template<typename T>
int gpucg_intprecong(int n1,int n2, int n3,int npsi,int ncong,
		     T h1,T h2,T h3,T c,T *x,int *keys,T *r,T *b,T *d,
		     T *work1,T *work2,T *work3)
{

  dim3  gridC(nblocksC, 1, 1);  
  dim3  threadsC(ELEMS_BLOCK, nseg_blockC , 1);

  gpuapply_hp<T>(n1,n2,n3,h1,h2,h3,c,
		 x,d,keys,
		 work1,work2,work3);

  //change the sign of the 0-th step such as to use axpy calls
  //d=d-b
  do_blasAxpy(npsi,(T)(-1.),b,d);
  cudaThreadSynchronize();

  //r=d
  cudaMemcpy(r,d,npsi*sizeof(T), cudaMemcpyDeviceToDevice);
  cudaThreadSynchronize();

  T rmr_new=do_blasDot(npsi, r, r);
  cudaThreadSynchronize();

  T alpha,beta,rmr_old;

  for(int i=0;i < ncong;++i)
    {
      gpuapply_hp<T>(n1,n2,n3,h1,h2,h3,c,
		     d,b,keys,
		     work1,work2,work3);

      alpha=rmr_new/do_blasDot(npsi, d, b);
      cudaThreadSynchronize();

      //here the sign is inverted because of the
      //mapping d -> -d => b -> -b
      
      //x=x-alpha*d
      do_blasAxpy(npsi,-alpha,d,x);
      cudaThreadSynchronize();

      if (i != ncong-1)
	{
	  //r=r-alpha*b (r does not change sign since also b is opposite)
	  do_blasAxpy(npsi,-alpha,b,r);
	  cudaThreadSynchronize();

	  rmr_old=rmr_new;
	  rmr_new=do_blasDot(npsi, r, r);
	  cudaThreadSynchronize();

	  beta=rmr_new/rmr_old;

	  //d=d+1/beta*r
	  do_blasAxpy(npsi,(T)(1)/beta,r,d);
	  cudaThreadSynchronize();

	  //d=beta*d
	  do_blasScal(npsi,beta,d);
	  cudaThreadSynchronize();
	}
  
	
    }


  wscalgpu<T> <<< gridC, threadsC >>>(x,h1,h2,h3,c,keys);
  cudaThreadSynchronize();


  return 0;

}


extern "C" 
void gpuprecond_(int *n1,int *n2, int *n3,int *npsi,
		 double *h1,double *h2,double *h3,
		 double **x,int **keys, 
		 double **r,double **b,double **d,
		 double **work1,double **work2,double **work3,
		 double *c,int *ncong, double *gnrm)
{
  
  /* printf("%i, %i, %i, %i,  %lf, %lf, %lf, %p, %p, %p, %p, %p %p, %p, %p, %lf, %i, %lf\n",*n1+1,*n2+1,*n3+1,*npsi,*ncong,
	 *h1,*h2,*h3,*c,
	 *x,*keys,*r,*b,*d,*work1,*work2,*work3,*gnrm);*/

  if(gpucg_precong<double>(*n1+1,*n2+1,*n3+1,*npsi,*ncong,
			   *h1,*h2,*h3,*c,
			   *x,*keys,*r,*b,*d,*work1,*work2,*work3,gnrm)!= 0)
    {
      printf("ERROR: GPU fulllocalhamiltonian\n ");
      return;
    } 

  CUERR;
  return; 
}

extern "C" 
void gpuintprecond_(int *n1,int *n2, int *n3,int *npsi,
		    double *h1,double *h2,double *h3,
		    double **x,int **keys, 
		    double **r,double **b,double **d,
		    double **work1,double **work2,double **work3,
		    double *c,int *ncong)
{
  
  /* printf("%i, %i, %i, %i,  %lf, %lf, %lf, %p, %p, %p, %p, %p %p, %p, %p, %lf, %i, %lf\n",*n1+1,*n2+1,*n3+1,*npsi,*ncong,
	 *h1,*h2,*h3,*c,
	 *x,*keys,*r,*b,*d,*work1,*work2,*work3,*gnrm);*/

  if(gpucg_intprecong<double>(*n1+1,*n2+1,*n3+1,*npsi,*ncong,
			      *h1,*h2,*h3,*c,
			      *x,*keys,*r,*b,*d,*work1,*work2,*work3)!= 0)
    {
      printf("ERROR: GPU fulllocalhamiltonian\n ");
      return;
    } 

  CUERR;
  return; 
}
