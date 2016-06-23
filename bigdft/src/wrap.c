#include <config.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
//#include "/home/ghasemi/Programs/amber11/AmberTools/src/nab/nabcode.h"
#include "nabcode.h"
//#include "/home/ghasemi/Programs/nab-5.1.2/include/nabcode.h"
//#include "/home/physik/roysh/MHOP/NAB/nab-5.1.2/src/nabcode.h"
//#include "static_in"
extern char NAB_rsbuf[];


static MOLECULE_T *m;
static PARMSTRUCT_T *prm1; //=m->m_prm;

static REAL_T m_xyz[2000];

static INT_T ier  ; //, mytaskid , numtasks;

//static REAL_T *m_xyz,  *f_xyz,  *v;

static REAL_T dgrad, fret, dummy[2];
//static REAL_T dfr = 0.0001 ;
//static INT_T maxit = 1000 ;  
//static REAL_T dgrad = 0.001 ;
//static INT_T __gdab0001__;
//static INT_T __gdab0002__;
//static INT_T __gdab0003__;
//static INT_T __it0001__;
//static INT_T __it0002__;
//static INT_T __it0003__;
//static REAL_T __ft0001__;
//static STRING_T *__st0001__ = NULL;
//static STRING_T *__st0002__ = NULL;
//static STRING_T *__st0003__ = NULL;
#define NUL '\0'
//int nab_init_(INT_T *natoms, REAL_T *rxyz, REAL_T  *fxyz ,  char *filename_t , int *ilenn,int *l_sat,char *sat)
void FC_FUNC_(nab_init,NAB_INIT)(char *args)
// , INT_T *iproc , INT_T *nproc , char *argv[])// , char *filename)
{
//void  nab_pdbread_(REAL_T * , char * , int * );
nabout = stdout; /*default*/

m = getpdb("init.pdb", NULL );
readparm( m, "struct.top");

setxyz_from_mol(  &m, NULL, m_xyz );

//mm_options( "cut=25.0, ntpr=10, nsnb=999, gamma_ln=5.0" );
//mm_options( "ntpr=10, gb=1, surften=0.00500,  kappa=0.10395, rgbmax=99., cut=99.0, diel=C ");
//mm_options( "ntpr=1000, gb=1, gbsa=1 , rgbmax=99.,cut=99.0, diel=C ");
mm_options( args);
//mm_options( "ntpr=1000, gb=0, gbsa=0 , rgbmax=99.,cut=99.0, diel=C ");
//mm_options( "cut=999., ntpr=10, nsnb=100, diel=R" );
mme_init( m, NULL, "::Z", m_xyz, NULL );
prm1=m->m_prm;

return ;
} //end of void  nab_init_
//******************************************************************************
void  FC_FUNC_(nab_finalize,NAB_FINALIZE)()
{
ier = 	mpifinalize(  ) ; 
}


void  FC_FUNC_(call_nab_gradient,CALL_NAB_GRADIENT)(REAL_T *rxyz , REAL_T *fxyz, REAL_T  *e_pos, INT_T *icount){

int i , j ; 
//e_pos = mme( rxyz, fxyz, ITEMP( __it0001__, 1));

*e_pos = mme( rxyz, fxyz, icount );
//printf( "...energy is %8.3f\n", *e_pos );
//j =3* *( NAB_mri( m, "natoms" ) );
//for (i=0;i<j ; i++){
//printf( "%16.10f\n", fxyz[i]);
//}
}

/*
void  call_nab_hessian_(REAL_T *rxyz,REAL_T *gxyz,REAL_T *hess,REAL_T  *amass) {
	REAL_T *grad=NULL;
	INT_T *descF_PxQ=NULL,*descF_1x1=NULL,*descG_PxQ=NULL,context_PxQ;
	INT_T context_1x1,context_Nx1,gridim,natom,iter;
	iter=0;
	//mme2(rxyz,gxyz,hess,amass,grad,descF_PxQ,descF_1x1,descG_PxQ,context_PxQ, \ 
	//	context_1x1,context_Nx1,gridim,natom,iter);
	//mme2(rxyz,gxyz,hess,amass,grad,descG_PxQ,descF_1x1,descG_PxQ,&context_PxQ,&context_1x1,&context_Nx1,&gridim,&natom,&iter);
//#ifdef SCALAPACK
//	printf("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n");
//#endif
}
*/
//**************************************************************************************************
void FC_FUNC_(call_nab_hessian,CALL_NAB_HESSIAN)(REAL_T *x,int *n,REAL_T *h,REAL_T *amass,REAL_T *g) {
	nmode_hess(x,n,mme2,h,amass,g);
	return;
}
//**************************************************************************************************
#define DTYPE_ (0)
#define CTXT_ (1)
#define M_ (2)
#define N_ (3)
#define MB_ (4)
#define NB_ (5)
#define RSRC_ (6)
#define CSRC_ (7)
#define LLD_ (8)
#define DLEN_ (9)

int nmode_hess( REAL_T *x, int *n, 
	   REAL_T ( *func)( REAL_T*, REAL_T*, REAL_T*, REAL_T*,
			    REAL_T*, int*, int*, int*,
			    int*, int*, int*, int*,
			    int*, int* ),REAL_T *h,REAL_T *m,REAL_T *g)

{
  REAL_T *work=NULL, *grad=NULL;
  REAL_T sumg, energy, t1, t2;
  int niter, i, j, k, mcopy, natom, ret_val, lwork, info;
  int mytaskid, numtasks, gridim;
  int context_PxQ=-1, context_1x1=-1, context_Nx1=-1;
  int descH_PxQ[DLEN_], descG_PxQ[DLEN_], descG_1x1[DLEN_];
  FILE *vfp;
  REAL_T pressure, temperature;
  char uplo, jobz;
  size_t ncopy;

#ifdef SCALAPACK
  int zero=0, one=1;
  int myrow, mycol, nprow, npcol;
  int myrowC, mycolC, nprowC, npcolC;
  int bs3, ierror, lld;
  int np, nq, nb, sizesytrd, nprocs, contextC, nrc, ldc, qrmem;
  int lldH_PxQ, lldG_PxQ, lldZ_PxQ, lldG_1x1;
  size_t locpH_PxQ, locpG_PxQ, locpG_1x1;
  size_t locqH_PxQ, locqG_PxQ, locqG_1x1;
  size_t sizeH_PxQ, sizeG_PxQ, sizeG_1x1;
  size_t adr, sizemqrleft;
  REAL_T *ptr, *reductarr=NULL, *eigvecrow=NULL;
#endif

  /* If PRINT_NM_TIMES is defined print some calculation times. */

#undef PRINT_NM_TIMES

  /* Get mytaskid and, if SCALAPACK is defined, numtasks. */

  mytaskid = get_mytaskid();

#ifdef SCALAPACK
  numtasks = get_numtasks();
#endif

  /* Allocate some dynamic vectors and matrices. */

  mcopy = *n;
  ncopy = *n;
  //m = vector( 1, ncopy );
  //v = vector( 1, ncopy );

#ifndef SCALAPACK

  /* If SCALAPACK is not defined, allocate full copies of g and h. */

  //g = vector( 1, ncopy );
  //h = vector( 0, ncopy*ncopy );

#else

  /*
   * If SCALAPACK is defined, allocate distributed copies of g and h,
   * as well as a single copy of grad.
   *
   * Create a general context.  Although context_PxQ does comprise all
   * of the processes, it appears that the general context must be
   * distinct from the context(s) of the matrices that participate
   * in the redistribution via pdgemr2d.
   *
   * The topologic layout of the general context does not appear
   * to be important, so this general context is row cyclic.
   *
   * It appears to be important that the most general context be created
   * first, followed by successively less general contexts.  For this
   * code the correct order of creation is Nx1 then PxQ then 1x1.  Each
   * context is subsumed by the earlier created contexts.  I don't know
   * what to do about overlapping contexts where one does not sumsume
   * the other.  Failure to the most general context first leads to
   * synchronization errors.
   */

  sl_init_(&context_Nx1, &numtasks, &one);

  /* Calculate the dimensions of the largest possible square process grid. */

  gridim = (int)(sqrt((REAL_T)numtasks)+0.5);
  if (gridim*gridim > numtasks) {
    gridim -= 1;
  }
  if (gridim == 0) {
    gridim = 1;
  }

  /*
   * Initialize the process grid for block cyclic distribution of matrices.
   * Note that a returned context of -1 indicates that the task is not
   * active on the process grid.
   */

  sl_init_(&context_PxQ, &gridim, &gridim);

  /* Initialize the process grid for a single (1x1) process. */

  sl_init_(&context_1x1, &one, &one);

  /*
   * Get the number of rows and columns on the block cyclic (PxQ) process grid,
   * as well as this task's row and column on the grid.
   */

  blacs_gridinfo_(&context_PxQ, &nprow, &npcol, &myrow, &mycol);

  /*
   * Get the blocksize for a square block.  Because in egb2 the
   * loop index i selects three rows of the Hessian, multiply
   * the blocksize by three.
   */

  bs3 = 3*get_blocksize();

  /*
   * If this task is on the process grid, set up the array descriptors.
   * If this task isn't on the process grid, set descZ_PxQ[CTXT_],
   * descG_PxQ[CTXT_] and descH_PxQ[CTXT_] to -1.  These values will
   * be used by pdgemr2d to determine activity on the grid.
   */

  if (context_PxQ >= 0) {

    /*
     * Allocate then distribute the Hessian matrix h, the eigenvector
     * matrix z and the gradient vector g on the block cyclic process grid.
     * The gradient vector g does not need to be allocated for columns other
     * than column 0, but the descinit_ function must be called for all columns.
     *
     * The numroc_ function is used to calculate the number of matrix
     * elements that are distributed across a PxQ processor grid.
     */

    locpG_PxQ = numroc_(&mcopy, &bs3, &myrow, &zero, &nprow);
    locqG_PxQ = numroc_(&one, &bs3, &mycol, &zero, &npcol);
    sizeG_PxQ = locpG_PxQ * locqG_PxQ;
    lldG_PxQ = locpG_PxQ;
    descinit_(descG_PxQ, &mcopy, &one, &bs3, &bs3,
	      &zero, &zero, &context_PxQ, &lldG_PxQ, &info);
    g = vector( 1, sizeG_PxQ );

    locpH_PxQ = numroc_(&mcopy, &bs3, &myrow, &zero, &nprow);
    locqH_PxQ = numroc_(&mcopy, &bs3, &mycol, &zero, &npcol);
    sizeH_PxQ = locpH_PxQ * locqH_PxQ;
    lldH_PxQ = locpH_PxQ;
    descinit_(descH_PxQ, &mcopy, &mcopy, &bs3, &bs3,
	      &zero, &zero, &context_PxQ, &lldH_PxQ, &info);
    h = vector( 0, sizeH_PxQ );
  
  } else {
    descG_PxQ[CTXT_] = -1;
    descH_PxQ[CTXT_] = -1;
  }

  /*
   * Get the number of rows and columns on the single process grid,
   * as well as this task's row and column on the grid.
   */

  blacs_gridinfo_(&context_1x1, &nprow, &npcol, &myrow, &mycol);

  /*
   * If this task is on the process grid, set up the array descriptors.
   * If this task isn't on the process grid, set descG_1x1[CTXT_] ] to -1.
   * This value will be used by pdgemr2d to determine activity on the grid.
   */

  if (context_1x1 >= 0) {

    /*
     * Allocate then distribute the gradient vector grad on the single
     * process grid.  The descinit_ function is called for this array
     * for only the task that is active on the grid.
     *
     * Also, for the other tasks that are not active, the vector grad
     * is allocated because the copy of grad from the single process
     * grid is broadcast to the other copies.  The size of this vector
     * for these other tasks is determined by the ncopy variable.
     *
     * The numroc_ function is used to calculate the number of matrix
     * elements that are distributed across a 1x1 processor grid.
     */

    locpG_1x1 = numroc_(&mcopy, &bs3, &myrow, &zero, &nprow);
    locqG_1x1 = numroc_(&one, &bs3, &mycol, &zero, &npcol);
    sizeG_1x1 = locpG_1x1 * locqG_1x1;
    lldG_1x1 = locpG_1x1;
    descinit_(descG_1x1, &mcopy, &one, &bs3, &bs3,
	      &zero, &zero, &context_1x1, &lldG_1x1, &info);
    grad = vector( 1, sizeG_1x1 );

  } else {
    descG_1x1[CTXT_] = -1;
    grad = vector( 1, ncopy );
  }

  /* Allocate the eigenvector row and reduction arrays. */

  eigvecrow = vector( 0, ncopy );
  reductarr = vector( 0, ncopy );
  
#endif /* ifndef SCALAPACK */

  niter = 1;

  t1 = seconds();

  /*
   * For non-ScaLAPACK execution, set some variables to values that
   * will select the proper sections of code below.
   */

#ifndef SCALAPACK

  gridim = 1;
  context_PxQ = 0;

#endif

  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

  /*
   * Compute the function value in "f",
   * its gradient (with respect to x) in "g", and
   * its hessian (with respect to x) in "h".
   * The atomic masses are returned in "m".
   * The number of atoms is returned in "natom".
   *
   * The grad, descG_PxQ, descG_1x1, descH_PxQ,
   * gridim, context_PxQ, context_1x1 and context_Nx1
   * calling parameters supply ScaLAPACK information,
   * or are dummy arguments for non-ScaLAPACK execution.
   *
   * The &g[1] and &m[1] calling parameters
   * map from 1..n indexing in this newton function
   * to 0..n-1 indexing in *func (which defaults to mme2).
   * This technique is not used for h because it must be
   * indexed from 1 to n.
   */

  energy=( *func )( x, &g[ 0 ], h, &m[ 0 ],
		      &grad[ 1 ], descG_PxQ, descG_1x1, descH_PxQ,
		      &context_PxQ, &context_1x1, &context_Nx1,
		      &gridim, &natom, &niter );
  return 1;
}//end nmode
//**************************************************************************************************
//void  call_nab_conj_(REAL_T *rxyz , REAL_T  *e_pos , REAL_T *dgradd , REAL_T *dfrr, INT_T *maxitt )
//{
//ier = conjgrad( rxyz, ITEMP( __it0001__, 3 *  *( NAB_mri( m, "natoms" ) ) ),e_pos, mme,  dgradd, dfrr , maxitt) ;
////ier = conjgrad( rxyz,&__gdab0001__,e_pos, mme,  &dgrad, &dfr , &maxit) ;
//}
//**************************************************************************************************
void  FC_FUNC_(nab_pdbwrite,NAB_PDBWRITE)(REAL_T *rxyz , char *filename ,int *ii, char *remline)
{
char *trim(char *) ; 
char *file1 ; 
int i, iii ; 
iii = *ii ; 
filename[iii]='\0' ;
remline[100]='\0' ;
setmol_from_xyz(  &m, NULL, rxyz );
 putpdb( filename, m, remline );
}


 void  FC_FUNC_(nab_pdbread,NAB_PDBREAD)(REAL_T *rxyz , char *filename , int *ii )
{
char *trim(char *) ; 
char *file1 ; 
int i, iii ;
iii = *ii ;
filename[iii]='\0' ;
file1 = trim(filename) ;
m = getpdb( filename, NULL );
setxyz_from_mol(  &m, NULL, rxyz );
}


 void  FC_FUNC_(nab_bonds_angles,NAB_BONDS_ANGLES)(INT_T *nbonds,INT_T *nangles ,INT_T *tbonds1, INT_T *tbonds2,INT_T *tangles1 ,  INT_T *tangles2, INT_T *tangles3) {
int i , ii ,iii ;
int *ik ;

	
*nbonds = m->m_prm->Nbonh ; 
*nbonds += m->m_prm->Nbona ; 
 *nangles = m->m_prm->Ntheth ;
 *nangles += m->m_prm->Ntheta ;
	ii = m->m_prm->Nbonh ;
	iii = m->m_prm->Nbona ;
for (i=0 ; i < ii  ; i ++){
tbonds1[i] = m->m_prm->BondHAt1[i] ;
tbonds2[i] = m->m_prm->BondHAt2[i] ;
}
for (i= 0 ; i < iii  ; i ++){
tbonds1[ii+i] = m->m_prm->BondAt1[i] ;
tbonds2[ii+i] = m->m_prm->BondAt2[i] ;
}
	ii =  m->m_prm->Ntheth ;
	iii = m->m_prm->Ntheta ;
for (i=0 ; i < ii  ; i ++){
tangles1[i] = m->m_prm->AngleHAt1[i] ;
tangles2[i] = m->m_prm->AngleHAt2[i] ;
tangles3[i] = m->m_prm->AngleHAt3[i] ;
}
for (i= 0 ; i < iii  ; i ++){
tangles1[ii+i] = m->m_prm->AngleAt1[i] ;
tangles2[ii+i] = m->m_prm->AngleAt2[i] ;
tangles3[ii+i] = m->m_prm->AngleAt3[i] ;
}
//
//	!ndHAt1 = (int *) get(sizeof(int)* prm->Nbonh);
//	!  prm->BondHAt2 = (int *) get(sizeof(int)* prm->Nbonh);
//	!  prm->BondHNum = (int *) get(sizeof(int)* prm->Nbonh);
//	!  prm->BondAt1 = (int *) get(sizeof(int)* prm->Nbona);
//	!  prm->BondAt2 = (int *) get(sizeof(int)* prm->Nbona);
//	!  prm->BondNum = (int *) get(sizeof(int)* prm->Nbona);
//	!  prm->AngleHAt1 = (int *) get(sizeof(int)* prm->Ntheth);
//	!  prm->AngleHAt2 = (int *) get(sizeof(int)* prm->Ntheth);
//	!  prm->AngleHAt3 = (int *) get(sizeof(int)* prm->Ntheth);
//	!  prm->AngleHNum = (int *) get(sizeof(int)* prm->Ntheth);
//	!  prm->AngleAt1 = (int *) get(sizeof(int)* prm->Ntheta);
//	!  prm->AngleAt2 = (int *) get(sizeof(int)*prm->Ntheta);
//	!  prm->AngleAt3 = (int *) get(sizeof(int)*prm->Ntheta);
//	!  prm->AngleNum = (int *) get(sizeof(int)*prm->Ntheta);
//	!
//
//
}

char *trim(char *str)
{
  	char *ibuf = str, *obuf = str;
      int i = 0, cnt = 0;

      /*
 *       **  Trap NULL
 *             */

      if (str)
      {
            /*
 *             **  Remove leading spaces (from RMLEAD.C)
 *                         */

            for (ibuf = str; *ibuf && isspace(*ibuf); ++ibuf)
                  ;
            if (str != ibuf)
                  memmove(str, ibuf, ibuf - str);

            /*
 *             **  Collapse embedded spaces (from LV1WS.C)
 *                         */

            while (*ibuf)
            {
                  if (isspace(*ibuf) && cnt)
                        ibuf++;
                  else
                  {
                        if (!isspace(*ibuf))
                              cnt = 0;
                        else
                        {
                              *ibuf = ' ';
                              cnt = 1;
                        }
                        obuf[i++] = *ibuf++;
                  }
            }
            obuf[i] = NUL;

            /*
 *             **  Remove trailing spaces (from RMTRAIL.C)
 *                         */

            while (--i >= 0)
            {
                  if (!isspace(obuf[i]))
                        break;
            }
            obuf[++i] = NUL;
      }
      return str;
}

