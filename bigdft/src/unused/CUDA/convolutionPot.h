/****c* CUDA/convolutionPot.h
** 
** AUTHOR
**  Matthieu Ospici
** 
** CHANGELOG
**  Started on  Tue Mar 11 14:55:07 2008 Matthieu Ospici
**
**  Last update Tue Mar 11 14:55:07 2008 Matthieu Ospici
**
** SOURCE
*/

#ifndef   	CONVOLUTIONPOT_H_
# define   	CONVOLUTIONPOT_H_

#include "structUtil.h"


int conv3dGPU_pot(int dn1,
		  int dn2,
		  int dn3,

		  float *t_psi_GPU,      		  
		  float *t_V_GPU,
		  float *t_Work_GPU,
		  
		  
		  const float *f_data1,
		  const float *f_data2,

		  int fsize1,
		  int lowfil1,
		  int lupfil1,

		  int fsize2,
		  int lowfil2,
		  int lupfil2,

		  unsigned int num_threads,
		  evalPerfGPU_t *evPerf);



void doMalloc(multiTab_t *t_m_dataIn,
	      multiTab_t *t_m_dataOut,
	      float **t_f_data,
	      int fsize,
	      unsigned int num_threads,
	      unsigned int num_GPU,
	      evalPerfGPU_t *t_evPerf);



//******** PRIVATE MEMBER, do not try to call *******
void createParam(param_t* param,
		 unsigned int num_threads,
		 unsigned int num_elements_dim1,
		 unsigned int num_elements_dim2,
		 unsigned int *numBlock,
		 unsigned int *numBlockDim2,
		 const float *f_data,
		 int fsize,
		 int lowfil, 
		 int lupfil);



void uniformiseTab( int* tab,unsigned int *offsetTab,int sizeTab,int lineSize);


void createTabs(int tabSeq[2][40],
		unsigned int currPosTab[2][40],
		unsigned int numThread,
		unsigned int sizeLine,
		short step);

#endif 	    /* !CONVOLUTION_H_ */

/****/
