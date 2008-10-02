/****c* CUDA/structUtil.h
**
** AUTHOR
**  Matthieu Ospici
** 
** CHANGELOG
**  Started on  Mon Feb 11 10:28:17 2008 Matthieu Ospici
**
**  Last update Mon Feb 11 10:28:17 2008 Matthieu Ospici
**
** SOURCE
*/

#ifndef   	STRUCTUTIL_H_
# define   	STRUCTUTIL_H_

//#include <stdbool.h>
//#define SIZE_SHARED_TOTAL  4080
#define MAX_LINE_SIZE 100
#define SIZE_SHARED_TOTAL  16*MAX_LINE_SIZE






typedef struct _evalPerfGPU
{
  float GPU_calc;
  float GPU_trsf_CPUGPU;
  float GPU_trsf_GPUCPU;
  float cuda_malloc;
  float host_malloc;
  float GPU_total_1;
} evalPerfGPU_t;





//#define fSIZE 40
typedef struct  multiTab
{
  float *tab;

  double *tab_d;
  unsigned int n1; //first dim
  unsigned int n2; //second dimension
  unsigned int n3;

  //  unsigned int n1_real;
  // unsigned int n2_real;
  

  float *GPU_data;

 
} multiTab_t;


typedef struct _threadParam
{
  multiTab_t *m_dataIn;
  multiTab_t *m_dataOut;
  const float *f_data;
  int fsize;
  unsigned int num_threads;
  evalPerfGPU_t *evPerf;
  unsigned int gpuID;
  unsigned int num_GPU;
} threadParam;



/*typedef struct _threadParamMalloc
{
  multiTab_t *m_dataIn;
  multiTab_t *m_dataOut;
  unsigned int gpuID;

} threadParamMalloc;*/

typedef struct _fetchTabs
{
  int tabSeq[2][40]; //for exmaple : 16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,15,
  unsigned  int currPosTab[2][40]; //for exmaple 0,16,32,48,64,80,96,112,128,144,160,176,192,208,224,240

} fetchTabs_t;

typedef struct _calcTabs
{
 int tabSeq[2][40]; 
  unsigned  int currPosTab[2][40];

} calcTabs_t;

typedef struct   _param
{
  unsigned int SIZE_SHARED_1 ;
  int SIZE_SHARED_2;
 

  fetchTabs_t fetchTabs;
  calcTabs_t calcTabs;

  
  unsigned  int sizeLineFirstBlocks;
  unsigned  int sizeLineLastBlock;


  int lowfil, lupfil; //structure of f_fct
  float f_fct[100];

  unsigned int lineLastBlock;


 

} param_t;



void printTab(const multiTab_t *toPrint,int ib =0,int jb=0);

bool compareTab(const multiTab_t *a, const multiTab_t *b,int *i,int *j,int ib =0,int jb=0); //true = equal

/*
inline float getTab(const multiTab_t* m,unsigned int i,unsigned int j)
{
  return m->tab[i*m->n2_real + j];
}

inline void setTab(multiTab_t* m,float val,unsigned int i,unsigned int j)
{
  m->tab[i*m->n2_real + j] = val;
}*/

#endif 	    /* !STRUCTUTIL_H_ */

/****/
