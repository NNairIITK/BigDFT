/****c* CUDA/convolution.h
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

#ifndef   	CONVOLUTION_H_
# define   	CONVOLUTION_H_

#include "structUtil.h"


void conv3dCPU(multiTab_t* m_dataIn,
	       multiTab_t* m_dataOut,
	       //  float *dataOut,
	       const double *f_data,
	       int fsize);


void conv3dCPU_s(multiTab_t* m_dataIn,
		 multiTab_t* m_dataOut,
		 //  float *dataOut,
		 const float *f_data,
		 int fsize);


void conv3dCPU_sd(multiTab_t* m_dataIn,
		  multiTab_t* m_dataOut,
		  const double *f_data,
		  int fsize);

/*int conv3dGPU(multiTab_t* m_dataIn,
	      multiTab_t* m_dataOut,
	      //     float *dataOut,
	      const float *f_data,
	      int fsize,
	      unsigned int num_threads,
	      evalPerfGPU_t *evPerf);*/

int conv3dGPU(multiTab_t* m_dataIn,
	      multiTab_t* m_dataOut,
	      //    float *dataOut,
	      const float *f_data,
	      int fsize,
	      int lowfil,
	      int lupfil,
	      unsigned int num_threads,
	      evalPerfGPU_t *evPerf);

int conv1dGPU(multiTab_t* m_dataIn,
	      multiTab_t* m_dataOut,
	      //    float *dataOut,
	      const float *f_data,
	      int fsize,
	      int lowfil,
	      int lupfil,
	      unsigned int num_threads,
	      evalPerfGPU_t *evPerf);

int conv3d_multi_GPU(multiTab_t *t_m_dataIn,
		     multiTab_t *t_m_dataOut,
		     float **t_f_data,
		     int fsize,
		     unsigned int num_threads,
		     unsigned int num_GPU,
		     evalPerfGPU_t *t_evPerf);

void doMalloc(multiTab_t *t_m_dataIn,
	      multiTab_t *t_m_dataOut,
	      float **t_f_data,
	      int fsize,
	      unsigned int num_threads,
	      unsigned int num_GPU,
	      evalPerfGPU_t *t_evPerf);




void freeMultiTab(multiTab_t* m);


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
