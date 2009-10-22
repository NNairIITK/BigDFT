#include <stdio.h>

#include <S_GPU/include/sg_common_def.h>



//============= PRECOND BINDING ====================
typedef struct sg_param_precond
{
  int n1, n2,  n3, npsi;
  double h1, h2, h3;
  double **x;
  int **keys;
  double **r, **b, **d;
  double **work1, **work2, **work3;
  double c;
  int ncong;
  double *gnrm;

  
} sg_param_precond_t;




void sg_callback_precond(void *param)
{
  sg_param_precond_t *locpar = ((sg_param_precond_t*)(param));
 
  double gnrmToAdd;
  
  gpuprecond_(&(locpar->n1),&(locpar->n2), &(locpar->n3),&(locpar->npsi),
	      &(locpar->h1),&(locpar->h2),&(locpar->h3),
	      locpar->x,locpar->keys, 
	      locpar->r,locpar->b,locpar->d,
	      locpar->work1,locpar->work2,locpar->work3,
	      &(locpar->c),&(locpar->ncong), &gnrmToAdd);

  *(locpar->gnrm) += gnrmToAdd;
}

void sg_precond_adapter__(int *n1,int *n2, int *n3,int *npsi,
			  double *h1,double *h2,double *h3,
			  double **x,int **keys, 
			  double **r,double **b,double **d,
			  double **work1,double **work2,double **work3,
			  double *c,int *ncong, double *gnrm,
			  sg_stream_ptr *stream)
{
  sg_param_precond_t param;

  param.n1 = *n1;
  param.n2 = *n2;
  param.n3 = *n3;
  param.npsi = *npsi;
  param.h1 = *h1;
  param.h2 = *h2;
  param.h3 = *h3;
  param.x = x;
  param.keys = keys;
  param.r = r;
  param.b = b;
  param.d = d;
  param.work1 = work1;
  param.work2 = work2;
  param.work3 = work3;
  param.c = *c;
  param.ncong = *ncong;
  param.gnrm = gnrm;
  

  sg_calc(&sg_callback_precond,&param,sizeof(sg_param_precond_t),*stream);
}

//============= END PRECOND BINDING ====================

//============= LOCHAM BINDING ====================

typedef struct sg_param_locham
{

  int n1, n2, n3;
  double h1, h2, h3;
  double **psi, **pot;
  int **keys; 
  double **work1, **work2, **work3;
  double *epot_sum, *ekin_sum;
  double occup_gpu;
} sg_param_locham_t;



void sg_callback_locham(void *param)
{
  sg_param_locham_t *locpar = ((sg_param_locham_t*)(param));
 


  double epotToAdd,ekinToAdd;
  gpulocham_(&locpar->n1,&locpar->n2, &locpar->n3,
	     &locpar->h1,&locpar->h2,&locpar->h3,
	     locpar->psi,locpar->pot,locpar->keys,			     
	     locpar->work1, locpar->work2, locpar->work3,
	     &epotToAdd,&ekinToAdd);


  *locpar->ekin_sum += locpar->occup_gpu*ekinToAdd;
 *locpar->epot_sum += locpar->occup_gpu*epotToAdd;
}

void sg_callback_fulllocham(void *param)
{
  sg_param_locham_t *locpar = ((sg_param_locham_t*)(param));
 


  double epotToAdd,ekinToAdd;
  gpufulllocham_(&locpar->n1,&locpar->n2, &locpar->n3,
	     &locpar->h1,&locpar->h2,&locpar->h3,
	     locpar->psi,locpar->pot,locpar->keys,			     
	     locpar->work1, locpar->work2, locpar->work3,
	     &epotToAdd,&ekinToAdd);


  *locpar->ekin_sum += locpar->occup_gpu*ekinToAdd;
 *locpar->epot_sum += locpar->occup_gpu*epotToAdd;
}


void sg_locham_adapter__(int *n1,int *n2, int *n3,
			double *h1,double *h2,double *h3,
			double **psi,double **pot,int **keys, 
			double **work1,double **work2,double **work3,
			double *epot_sum,double *ekin_sum,
			double *occup_gpu,
			sg_stream_ptr *stream)
{
  sg_param_locham_t param;

  
  param.n1 = *n1;
  param.n2 = *n2;
  param.n3 = *n3;

  param.h1 = *h1;
  param.h2 = *h2;
  param.h3 = *h3;
  param.psi = psi;
  param.pot = pot;
  param.keys = keys;
  param.work1 = work1;
  param.work2 = work2;
  param.work3 = work3;
  param.epot_sum = epot_sum;
  param.ekin_sum = ekin_sum;
  param.occup_gpu = *occup_gpu;

  
  sg_calc(&sg_callback_locham,&param,sizeof(sg_param_locham_t),*stream);
}


void sg_fulllocham_adapter__(int *n1,int *n2, int *n3,
			double *h1,double *h2,double *h3,
			double **psi,double **pot,int **keys, 
			double **work1,double **work2,double **work3,
			double *epot_sum,double *ekin_sum,
			double *occup_gpu,
			sg_stream_ptr *stream)
{
  sg_param_locham_t param;

  
  param.n1 = *n1;
  param.n2 = *n2;
  param.n3 = *n3;

  param.h1 = *h1;
  param.h2 = *h2;
  param.h3 = *h3;
  param.psi = psi;
  param.pot = pot;
  param.keys = keys;
  param.work1 = work1;
  param.work2 = work2;
  param.work3 = work3;
  param.epot_sum = epot_sum;
  param.ekin_sum = ekin_sum;
  param.occup_gpu = *occup_gpu;

  
  sg_calc(&sg_callback_fulllocham,&param,sizeof(sg_param_locham_t),*stream);
}


//============= END LOCHAM BINDING ====================

//=============  LOCDEN BINDING ====================

typedef struct sg_param_locden
{
  int n1,n2, n3,norbp,nspin;
  double h1, h2, h3;
  double *occup, *spinsgn;
  double **psi;
  int **keys; 
  double **work1,**work2;
  double **rho;
} sg_param_locden_t;



void sg_callback_locden(void *param)
{
  sg_param_locden_t *locpar = ((sg_param_locden_t*)(param));
 

  gpulocden_(&locpar->n1,&locpar->n2,&locpar->n3,&locpar->norbp,&locpar->nspin,
	     &locpar->h1,&locpar->h2,&locpar->h3,
	     locpar->occup,locpar->spinsgn,
	     locpar->psi,locpar->keys, 
	     locpar->work1,locpar->work2,
	     locpar->rho);
}


void  sg_locden_adapter__(int *n1,int *n2, int *n3,int *norbp,int *nspin,
			  double *h1,double *h2,double *h3,
			  double *occup,double *spinsgn,
			  double **psi,int **keys, 
			  double **work1,double **work2,
			  double **rho,
			  sg_stream_ptr *stream)
{
  
  sg_param_locden_t param;

  param.n1 = *n1;
  param.n2 = *n2;
  param.n3 = *n3;
  param.norbp = *norbp;
  param.nspin = *nspin;
  param.h1 = *h1;
  param.h2 = *h2;
  param.h3 = *h3;
  param.occup = occup;
  param.spinsgn = spinsgn;
  param.psi = psi;
  param.keys = keys;
  param.work1 = work1;
  param.work2 = work2;
  param.rho = rho;

sg_calc(&sg_callback_locden,&param,sizeof(sg_param_locden_t),*stream);
}
