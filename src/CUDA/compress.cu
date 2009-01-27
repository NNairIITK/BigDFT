 /****u* CUDA/1Dconv_new.cu
**
** 
** AUTHOR
**  Luigi Genovese
**
** SOURCE
*/
  
#include <stdio.h>
#include "commonDef.h"


//parameters for the kernel (global variables)
//__constant__ int keys[MAX_CONSTANT_SIZE/sizeof(int)];

#include "kernels_compress.hcu"

/*





// copy the keys parameter into the GPU memory
void GPUgridkeys(int* keys,
		 int nseg_c,int nseg_f,
		 int* keyv_c,int* keyg_c,int* keyv_f,int* keyg_f)
{

  //printf("Anseg_block %i,nblocks %i,nseg_c %i \n",
  //nseg_block,nblocks,nseg_c);


  int nseg_tot=-1;
  int segment_elements;
  //assign the keys values, coarse part
  for(int i=0; i<nseg_c; ++i)
    {
      //printf("i %i,nblocks %i,nseg_c %i,segelems %i,elems_block %i \n",
      //i,nblocks,nseg_c,keyg_c[2*i+1]-keyg_c[2*i],ELEMS_BLOCK);
  
      //control if the number of elements of a segment is less than
      //the number of elements per block
      segment_elements=keyg_c[2*i+1]-keyg_c[2*i];

      if (segment_elements <= ELEMS_BLOCK)
	{

	  printf("i %i,j %i,nseg_tot %i \n",i,i,nseg_tot);
	  nseg_tot+=1;
	  //put the values of the keyg and keyv arrays  
	  keys[4*nseg_tot]=1; //which means coarse segment
	  keys[4*nseg_tot+1]=segment_elements;
	  keys[4*nseg_tot+2]=keyv_c[i];
	  keys[4*nseg_tot+3]=keyg_c[2*i];
	  
	}
      else
	{
	  //printf("i %i,segelems/elemsblock %i,nseg_tot %i \n",
	  //i,segment_elements/ELEMS_BLOCK,nseg_tot);
	  nseg_tot+=1;
	  //put the values of the keyg and keyv arrays  
	  keys[4*nseg_tot]=1; //which means coarse segment
	  keys[4*nseg_tot+1]=ELEMS_BLOCK;
	  keys[4*nseg_tot+2]=keyv_c[i];
	  keys[4*nseg_tot+3]=keyg_c[2*i];	  

	  for(int j=1;j < segment_elements/ELEMS_BLOCK;++j)
	    {
	      //printf("i %i,j %i,nseg_tot %i \n",i,j,nseg_tot);

	      nseg_tot+=1;
	      //put the values of the keyg and keyv arrays  
	      keys[4*nseg_tot]=1; //which means coarse segment
	      keys[4*nseg_tot+1]=ELEMS_BLOCK;
	      keys[4*nseg_tot+2]=keyv_c[i]+j*ELEMS_BLOCK+1;
	      keys[4*nseg_tot+3]=keyg_c[2*i]+j*ELEMS_BLOCK+1;
	    }
	  int j=segment_elements/ELEMS_BLOCK;
	  nseg_tot+=1;
	  //put the values of the keyg and keyv arrays  
	  printf("i %i,jAA %i,nseg_tot %i,maxsize %i \n",
		 i,j,nseg_tot,MAX_CONSTANT_SIZE/sizeof(int));
	  keys[4*nseg_tot]=1; //which means coarse segment
	  keys[4*nseg_tot+1]=
	    segment_elements-j*ELEMS_BLOCK;
	  keys[4*nseg_tot+2]=keyv_c[i]+j*ELEMS_BLOCK+1;
	  keys[4*nseg_tot+3]=keyg_c[2*i]+j*ELEMS_BLOCK+1;
	}


    }

  printf("nseg_block %i,nblocks %i,nseg_tot %i \n",
	 nseg_block,nblocks,nseg_tot);


  //assign the keys values, fine part
  for(int i=0; i<nseg_f; ++i)
    {
      //control if the number of elements of a segment is less than
      //the number of elements per block
      segment_elements=keyg_f[2*i+1]-keyg_f[2*i];
      printf("i %i,j %i,nseg_tot %i \n",i,i,nseg_tot);

      if (segment_elements <= ELEMS_BLOCK)
	{
	  nseg_tot+=1;
	  //put the values of the keyg and keyv arrays  
	  keys[4*nseg_tot]=7; //which means fine segment
	  keys[4*nseg_tot+1]=segment_elements;
	  keys[4*nseg_tot+2]=keyv_f[i];
	  keys[4*nseg_tot+3]=keyg_f[2*i];
	  
	}
      else
	{
	  nseg_tot+=1;
	  //put the values of the keyg and keyv arrays  
	  keys[4*nseg_tot]=7; //which means fine segment
	  keys[4*nseg_tot+1]=ELEMS_BLOCK;
	  keys[4*nseg_tot+2]=keyv_f[i];
	  keys[4*nseg_tot+3]=keyg_f[2*i];	  

	  for(int j=1;j < segment_elements/ELEMS_BLOCK;++j)
	    {
	      printf("i %i,j %i,nseg_tot %i \n",i,j,nseg_tot);
	      nseg_tot+=1;
	      //put the values of the keyg and keyv arrays  
	      keys[4*nseg_tot]=7; //which means fine segment
	      keys[4*nseg_tot+1]=ELEMS_BLOCK;
	      keys[4*nseg_tot+2]=keyv_f[i]+j*ELEMS_BLOCK+1;
	      keys[4*nseg_tot+3]=keyg_f[2*i]+j*ELEMS_BLOCK+1;
	    }
	  int j=segment_elements/ELEMS_BLOCK;
	  nseg_tot+=1;
	  //put the values of the keyg and keyv arrays  
	  keys[4*nseg_tot]=7; //which means fine segment
	  keys[4*nseg_tot+1]=
	    segment_elements-j*ELEMS_BLOCK;
	  keys[4*nseg_tot+2]=keyv_f[i]+j*ELEMS_BLOCK+1;
	  keys[4*nseg_tot+3]=keyg_f[2*i]+j*ELEMS_BLOCK+1;
	}
    }

  //number of segments treated by each block;
  nseg_block =(nseg_tot+1)/NBLOCKS_MAX +1;

  nblocks=min(NBLOCKS_MAX,nseg_tot+1);

  if(4*nblocks*nseg_block > MAX_CONSTANT_SIZE/sizeof(int))
    {
      printf("ERROR: Too many segments for Constant Memory!\n");

      return;
    } 

  //refill the last part of the segments
  for(int i=nseg_tot+1; i< nblocks*nseg_block; ++i)
    {
      keys[4*i]=0; //which means empty segment
      keys[4*i+1]=0;
      keys[4*i+2]=0;
      keys[4*i+3]=0;
    }

    

}

extern "C"
void copygpukeys_(int *nseg_c,int *nseg_f,
		int *keyv_c,
		int *keyg_c,
		int *keyv_f,
		int *keyg_f)
{

  //create the keys
  int keysCPU[MAX_CONSTANT_SIZE/sizeof(int)];

  GPUgridkeys(keysCPU,*nseg_c,*nseg_f,keyv_c,keyg_c,keyv_f,keyg_f);
  //send them to constant memory, once and for all
  //here the size can be reduced according to the actual number of
  //segments

  //if(cudaMemcpyToSymbol(*keys,&keysCPU,4*nseg_block*nblocks*sizeof(int)) != 0)
  //  {
  //    printf("MemcpyToSymbol error\n");
  //
  //return;
  //}

}

*/

extern "C"
void uncompressgpu_(int* n1, int* n2,int *n3 ,
		    double **psi, double **psig, int **keys)
{
  if(uncompressgpu<double>(*n1+1,*n2+1,*n3+1,*psi,*psig,*keys) != 0)
    {
      printf("MemcpyToSymbol error\n");

      return;
    }

}

extern "C"
void compressgpu_(int* n1, int* n2,int *n3 ,
		    double **psig, double **psi, int **keys)
{
  if(compressgpu<double>(*n1+1,*n2+1,*n3+1,*psig,*psi,*keys) != 0)
    {
      printf("MemcpyToSymbol error\n");

      return;
    }

}

extern "C"
void copynseginfo_(int *nblocks, int *nseg_block)
{
  nseg_blockC = *nseg_block;  
  nblocksC = *nblocks;

  //printf("nsegC %i, nblC %i,\n",nseg_blockC,nblocksC);
}


extern "C"
void readkeysinfo_(int *elems_block, int *nblocks_max)
{
  *elems_block=ELEMS_BLOCK;
  *nblocks_max=NBLOCKS_MAX;
}

