/*!
  accumulate the partial density of a decompressed, fscf psi
  add a work array to do not overwrite psifscf, since it is needed for
  the hamiltonian
 @author
    Copyright (C) 2010-2011 BigDFT group
    This file is distributed under the terms of the
    GNU General Public License, see ~/COPYING file
    or http://www.gnu.org/copyleft/gpl.txt .
    For the list of contributors, see ~/AUTHORS 

*/
template<typename T>
int densityaccumulation(int n1,int n2, int n3,
			dim3 gridMF1, dim3 gridMF2, dim3 gridMF3,
			dim3 threadsMF1, dim3 threadsMF2, dim3 threadsMF3,
			T *psifscf,T *psisq,T *work,T hfac,T *rhop)
{

  //calculate the MF transformation
  magicfilter1d<T> <<< gridMF3, threadsMF3 >>>(2*n3,4*n1*n2,psifscf,psisq,6);
  cudaThreadSynchronize();

  magicfilter1d<T> <<< gridMF2, threadsMF2 >>>(2*n2,4*n1*n3,psisq,work,7);
  cudaThreadSynchronize();

  magicfilter1d_den<T> <<< gridMF1, threadsMF1 >>>(2*n1,4*n2*n3,work,psisq,8);
  cudaThreadSynchronize();
 
  //accumulate the density with cuBlas
  //rhop=rhop+hfac*psisq
  do_blasAxpy(8*n1*n2*n3,hfac,psisq,rhop);

  return 0;

}

template<typename T>
int localpartialdensity(int n1,int n2, int n3, int nspin, T norbp,
			T h1,T h2,T h3,
			T *occup, T *spinsgn,
			T **psi,T *rhop,int *keys,
			T * work1,T *work2)
{

  
  dim3  gridWT3(linecutsWTz,  numBlocksWTz, 1);  
  dim3  threadsWT3(HALF_WARP_SIZE, num_halfwarpsWTz , 1);
  dim3  gridWT2(linecutsWTy,  numBlocksWTy, 1);  
  dim3  threadsWT2(HALF_WARP_SIZE, num_halfwarpsWTy , 1);
  dim3  gridWT1(linecutsWTx,  numBlocksWTx, 1);  
  dim3  threadsWT1(HALF_WARP_SIZE, num_halfwarpsWTx , 1);

  dim3  gridMF3(linecutsMFz,  numBlocksMFz, 1);  
  dim3  threadsMF3(HALF_WARP_SIZE, num_halfwarpsMFz , 1);
  dim3  gridMF2(linecutsMFy,  numBlocksMFy, 1);  
  dim3  threadsMF2(HALF_WARP_SIZE, num_halfwarpsMFy , 1);
  dim3  gridMF1(linecutsMFx,  numBlocksMFx, 1);  
  dim3  threadsMF1(HALF_WARP_SIZE, num_halfwarpsMFx , 1);

  dim3  gridC(nblocksC, 1, 1);  
  dim3  threadsC(ELEMS_BLOCK, nseg_blockC , 1);

  //set the density value to zero (eventually it should be 10^-20)
  cudaMemset((void*) rhop,0,8*n1*n2*n3*nspin*sizeof(T));
  cudaThreadSynchronize();


  T hfac;

  int iaddjmp;

  for(int iorb=0;iorb< norbp;++iorb)
    {

      hfac=occup[iorb]/(h1*h2*h3);

      iaddjmp=(spinsgn[iorb] > 0. ? 0 : 8*n1*n2*n3);

      //uncompression
      //set the value of the psig array to zero
      cudaMemset((void*) work1,0,8*n1*n2*n3*sizeof(T));
      cudaThreadSynchronize();

      uncompresscoarsefine<T> <<< gridC, threadsC >>>(n1,n2,n3,psi[iorb],work1,keys);
      cudaThreadSynchronize();
  
      //wavelet transformation.
      //overwrite the psi array with psifscf
      waveletstosf<T>(n1,n2,n3,gridWT1,gridWT2,gridWT3,
		      threadsWT1,threadsWT2,threadsWT3,
		      work1,psi[iorb]);


      densityaccumulation<T>(n1,n2,n3,gridMF1,gridMF2,gridMF3,
			     threadsMF1,threadsMF2,threadsMF3,
			     psi[iorb],work1,work2,hfac,&rhop[iaddjmp]);

    }

  return 0;

}

extern "C" 
void gpulocden_(int *n1,int *n2, int *n3,int *norbp,int *nspin,
		double *h1,double *h2,double *h3,
		double *occup,double *spinsgn,
		double **psi,int **keys, 
		double **work1,double **work2,
		double **rho)
{
  
  //  printf("%i %i %i %i %i %lf %lf %lf %lf %lf %p %p %p %p %p",
  if(localpartialdensity<double>(*n1+1,*n2+1,*n3+1,*nspin,*norbp,
				 *h1,*h2,*h3,occup,spinsgn,
				 psi,*rho,*keys,
				 *work1,*work2)!= 0)
    {
      printf("ERROR: GPU localpartialdensity\n ");
      return;
    } 


  return; 
}
