//! @file
//! binding s_gpu
//!
//! @author
//!    Copyright (C) 2009-2011 BigDFT group 
//!    This file is distributed under the terms of the
//!    GNU General Public License, see ~/COPYING file
//!    or http://www.gnu.org/copyleft/gpl.txt .
//!    For the list of contributors, see ~/AUTHORS 


#include <stdio.h>
#include <config.h>

#include <s_gpu.h>
#include <sg_common_def.h>



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

typedef struct sg_param_precondprecond
{
  int n1,n2,n3,nvctr_c,nvctr_f,nseg_c,nseg_f,ncplx,hybrid_on;
  double h1, h2, h3, cprecr;
  int **modul1,**modul2,**modul3,**keyg,**keyv;
  double **x, **scal, **b, **af, **bf, **cf, **ef, **kern_k1, **kern_k2, **kern_k3, **z1, **z3, **x_c, **psifscf,** ww;
  
} sg_param_precondprecond_t;


typedef struct sg_param_intprecond
{
  int n1, n2,  n3, npsi;
  double h1, h2, h3;
  double **x;
  int **keys;
  double **r, **b, **d;
  double **work1, **work2, **work3;
  double c;
  int ncong;
  
} sg_param_intprecond_t;

void gpuprecond_(int *n1,int *n2, int *n3,int *npsi,
		 double *h1,double *h2,double *h3,
		 double **x,int **keys, 
		 double **r,double **b,double **d,
		 double **work1,double **work2,double **work3,
		 double *c,int *ncong, double *gnrm);

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

void FC_FUNC_(precond_preconditioner_wrapper,PRECOND_PRECONDITIONER_WRAPPER)
     (int *hybrid_on, int *n1, int *n2, int *n3, int *nvctr_c, int *nvctr_f, int *nseg_c, int *nseg_f, int *ncplx,
     double *cprecr, double *hx, double *hy, double *hz, double **scal, int **keyg, int **keyv, int **modul1, int **modul2, int **modul3,
      double **af, double **bf, double **cf, double **ef, double **kern_k1, double **kern_k2, double **kern_k3, double **z1, double **z3, double **x_c, double **psifscf, double **ww, double **x, double **b);

void sg_callback_precondprecond(void *param)
{
  sg_param_precondprecond_t *locpar = ((sg_param_precondprecond_t*)(param));
 
  FC_FUNC_(precond_preconditioner_wrapper,PRECOND_PRECONDITIONER_WRAPPER)(&(locpar->hybrid_on),
				 &(locpar->n1),&(locpar->n2),&(locpar->n3),
				 &(locpar->nvctr_c),&(locpar->nvctr_f),
				 &(locpar->nseg_c),&(locpar->nseg_f),
				 &(locpar->ncplx),&(locpar->cprecr),
				 &(locpar->h1),&(locpar->h2),&(locpar->h3),
				 locpar->scal,locpar->keyg,locpar->keyv,
				 locpar->modul1,locpar->modul2,locpar->modul3,
				 locpar->af,locpar->bf,locpar->cf,locpar->ef,
				 locpar->kern_k1,locpar->kern_k2,locpar->kern_k3,
				 locpar->z1,locpar->z3,locpar->x_c,locpar->psifscf,
				 locpar->ww,
				 locpar->x,locpar->b);

}

void gpuintprecond_(int *n1,int *n2, int *n3,int *npsi,
		    double *h1,double *h2,double *h3,
		    double **x,int **keys, 
		    double **r,double **b,double **d,
		    double **work1,double **work2,double **work3,
		    double *c,int *ncong);

void sg_callback_intprecond(void *param)
{
  sg_param_intprecond_t *locpar = ((sg_param_intprecond_t*)(param));
 
  gpuintprecond_(&(locpar->n1),&(locpar->n2), &(locpar->n3),&(locpar->npsi),
		 &(locpar->h1),&(locpar->h2),&(locpar->h3),
		 locpar->x,locpar->keys, 
		 locpar->r,locpar->b,locpar->d,
		 locpar->work1,locpar->work2,locpar->work3,
		 &(locpar->c),&(locpar->ncong));
}


void FC_FUNC_(sg_precond_adapter,SG_PRECOND_ADAPTER)(int *n1,int *n2, int *n3,int *npsi,
			  double *h1,double *h2,double *h3,
			  double **x,int **keys, 
			  double **r,double **b,double **d,
			  double **work1,double **work2,double **work3,
			  double *c,int *ncong, double *gnrm,
			  sg_stream_ptr_t *stream)
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
  

  sg_gpu_send_calc(&sg_callback_precond,&param,sizeof(sg_param_precond_t),*stream);
}

void FC_FUNC_(sg_intprecond_adapter,SG_INTPRECOND_ADAPTER)(int *n1,int *n2, int *n3,int *npsi,
			     double *h1,double *h2,double *h3,
			     double **x,int **keys, 
			     double **r,double **b,double **d,
			     double **work1,double **work2,double **work3,
			     double *c,int *ncong,
			     sg_stream_ptr_t *stream)
{
  sg_param_intprecond_t param;

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

  sg_gpu_send_calc(&sg_callback_intprecond,&param,sizeof(sg_param_intprecond_t),*stream);
}

void FC_FUNC_(sg_precond_preconditioner_adapter,SG_PRECOND_PRECONDITIONER_ADAPTER)(int *hybrid_on,int *n1, int *n2,int *n3,
					 int *nvctr_c,int *nvctr_f,
					 int *nseg_c,int *nseg_f,int *ncplx,
					 double *cprecr,
					 double *h1,double *h2, double *h3,
					 double **scal,
					 int **keyg,int **keyv,
					 int **modul1, int **modul2,int **modul3,
					 double **af,double **bf,double **cf,
					 double **ef,double **kern_k1,double **kern_k2,
					 double **kern_k3,double **z1,double **z3,
					 double **x_c,double **psifscf,double **ww,
					 double **x,double **b,
					 sg_stream_ptr_t *stream)
{
  sg_param_precondprecond_t param;

  param.hybrid_on = *hybrid_on;
  param.n1 = *n1;
  param.n2 = *n2;
  param.n3 = *n3;
  param.nvctr_c = *nvctr_c;
  param.nvctr_f = *nvctr_f;
  param.nseg_c = *nseg_c;
  param.nseg_f = *nseg_f;
  param.ncplx = *ncplx;
  param.h1 = *h1;
  param.h2 = *h2;
  param.h3 = *h3;
  param.scal = scal;
  param.cprecr = *cprecr;
  param.keyg = keyg;
  param.keyv = keyv;
  param.modul1 = modul1;
  param.modul2 = modul2;
  param.modul3 = modul3;
  param.af = af;
  param.bf = bf;
  param.cf = cf;
  param.ef = ef;
  param.kern_k1 = kern_k1;
  param.kern_k2 = kern_k2;
  param.kern_k3 = kern_k3;
  param.z1 = z1;
  param.z3 = z3;
  param.psifscf = psifscf;
  param.ww = ww;
  param.x = x;
  param.b = b;
  sg_gpu_send_calc(&sg_callback_precondprecond,&param,sizeof(sg_param_precondprecond_t),*stream);
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


void gpulocham_(int *n1,int *n2, int *n3,
		double *h1,double *h2,double *h3,
		double **psi,double **pot,int **keys, 
		double **work1,double **work2,double **work3,
		double *epot,double *ekin);

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

void gpufulllocham_(int *n1,int *n2, int *n3,
		    double *h1,double *h2,double *h3,
		    double **psiw,double **pot,int **keys, 
		    double **psi,double **out,
		    double **work,
		    double *epot,double *ekinpot);

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


void FC_FUNC_(sg_locham_adapter,SG_LOCHAM_ADAPTER)(int *n1,int *n2, int *n3,
			double *h1,double *h2,double *h3,
			double **psi,double **pot,int **keys, 
			double **work1,double **work2,double **work3,
			double *epot_sum,double *ekin_sum,
			double *occup_gpu,
			sg_stream_ptr_t *stream)
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

  
  sg_gpu_send_calc(&sg_callback_locham,&param,sizeof(sg_param_locham_t),*stream);
}


void FC_FUNC_(sg_fulllocham_adapter,SG_FULLLOCHAM_ADAPTER)(int *n1,int *n2, int *n3,
			double *h1,double *h2,double *h3,
			double **psi,double **pot,int **keys, 
			double **work1,double **work2,double **work3,
			double *epot_sum,double *ekin_sum,
			double *occup_gpu,
			sg_stream_ptr_t *stream)
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

  
  sg_gpu_send_calc(&sg_callback_fulllocham,&param,sizeof(sg_param_locham_t),*stream);
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


void gpulocden_(int *n1,int *n2, int *n3,int *norbp,int *nspin,
		double *h1,double *h2,double *h3,
		double *occup,double *spinsgn,
		double **psi,int **keys, 
		double **work1,double **work2,
		double **rho);

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


void  FC_FUNC_(sg_locden_adapter,SG_LOCDEN_ADAPTER)(int *n1,int *n2, int *n3,int *norbp,int *nspin,
			  double *h1,double *h2,double *h3,
			  double *occup,double *spinsgn,
			  double **psi,int **keys, 
			  double **work1,double **work2,
			  double **rho,
			  sg_stream_ptr_t *stream)
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

  sg_gpu_send_calc(&sg_callback_locden,&param,sizeof(sg_param_locden_t),*stream);
}
