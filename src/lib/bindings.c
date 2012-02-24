#include <config.h>

#ifdef HAVE_GLIB
#include <glib-object.h>
#endif

#include "bigdft.h"
#include "bindings.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>


 /* Duplicate functions in C due to multiple interface definition in Fortran */
void FC_FUNC_(inquire_pointer1, INQUIRE_POINTER1)(void *pt, void *add, int *size)
{
  memcpy(pt, add, sizeof(void*) * *size);
}
void FC_FUNC_(inquire_pointer2, INQUIRE_POINTER2)(void *pt, void *add, int *size)
{
  memcpy(pt, add, sizeof(void*) * *size);
}
void FC_FUNC_(inquire_pointer3, INQUIRE_POINTER3)(void *pt, void *add, int *size)
{
  memcpy(pt, add, sizeof(void*) * *size);
}
void FC_FUNC_(inquire_pointer4, INQUIRE_POINTER4)(void *pt, void *add, int *size)
{
  memcpy(pt, add, sizeof(void*) * *size);
}
void FC_FUNC_(inquire_pointer5, INQUIRE_POINTER5)(void *pt, void *add, int *size)
{
  memcpy(pt, add, sizeof(void*) * *size);
}


void FC_FUNC_(inputs_new, INPUTS_NEW)(void *in);
void FC_FUNC_(inputs_free, INPUTS_FREE)(void *in);
void FC_FUNC_(inputs_set_radical, INPUTS_SET_RADICAL)(void *in,
                                                      const gchar *rad, int *len);
void FC_FUNC_(inputs_get_dft, INPUTS_GET_DFT)(const void *in, double *hx, double *hy, double *hz, double *crmult, double *frmult, int *ixc, int *ncharge, double *elecfield, int *nspin, int *mpol, double *gnrm_cv, int *itermax, int *nrepmax, int *ncong, int *idsx, int *dispersion, int *inputPsiId, int *output_wf_format, int *output_grid, double *rbuf, int *ncongt, int *norbv, int *nvirt, int *nplot, int *disableSym);
void FC_FUNC_(inputs_get_mix, INPUTS_GET_MIX)(void *in, int *iscf, int *itrpmax,
                                              int *norbsempty, int *occopt, double *alphamix,
                                              double *rpnrm_cv, double *gnrm_startmix,
                                              double *Tel, double *alphadiis);
void FC_FUNC_(inputs_get_geopt, INPUTS_GET_GEOPT)(void *in, char *geopt_approach,
                                                  int *ncount_cluster_x, double *frac_fluct,
                                                  double *forcemax, double *randdis,
                                                  double *betax, int *history, int *ionmov,
                                                  double *dtion, double *strtarget,
                                                  f90_pointer_double *qmass);
void FC_FUNC_(inputs_parse_params, INPUTS_PARSE_PARAMS)(void *in,
                                                        int *iproc, int *dump);
void FC_FUNC_(inputs_get_files, INPUTS_GET_FILES)(const void *in, int *files);
void FC_FUNC_(inputs_parse_add, INPUTS_PARSE_ADD)(void *in, const void *sym,
                                                  const gchar *geocode, const double *alat,
                                                  int *iproc, int *dump);

static BigDFT_Inputs* bigdft_inputs_init()
{
  BigDFT_Inputs *in;

  in = g_malloc(sizeof(BigDFT_Inputs));
  memset(in, 0, sizeof(BigDFT_Inputs));
  in->data = (void*)0;
  F90_1D_POINTER_INIT(&in->qmass);

  return in;
}
static void bigdft_inputs_dispose(BigDFT_Inputs *in)
{
  g_free(in->data);
  g_free(in);
}

BigDFT_Inputs* bigdft_inputs_new(const gchar *naming)
{
  BigDFT_Inputs *in;
  int iproc = 0, len, dump = 0;

  in = bigdft_inputs_init();
  FC_FUNC_(inputs_new, INPUTS_NEW)(&in->data);
  if (naming && naming[0])
    {
      len = strlen(naming);
      FC_FUNC_(inputs_set_radical, INPUTS_SET_RADICAL)(in->data, naming, &len);
    }
  else
    {
      len = 0;
      FC_FUNC_(inputs_set_radical, INPUTS_SET_RADICAL)(in->data, " ", &len);
    }
  FC_FUNC_(inputs_parse_params, INPUTS_PARSE_PARAMS)(in->data, &iproc, &dump);

  FC_FUNC_(inputs_get_files, INPUTS_GET_FILES)(in->data, &in->files);
  FC_FUNC_(inputs_get_dft, INPUTS_GET_DFT)(in->data, in->h, in->h + 1, in->h + 2,
                                           &in->crmult, &in->frmult, &in->ixc,
                                           &in->ncharge, in->elecfield, &in->nspin,
                                           &in->mpol, &in->gnrm_cv, &in->itermax,
                                           &in->nrepmax, &in->ncong, &in->idsx,
                                           &in->dispersion, &in->inputPsiId,
                                           &in->output_wf_format, &in->output_grid,
                                           &in->rbuf, &in->ncongt, &in->norbv, &in->nvirt,
                                           &in->nplot, &in->disableSym);
  FC_FUNC_(inputs_get_mix, INPUTS_GET_MIX)(in->data, &in->iscf, &in->itrpmax,
                                           &in->norbsempty, (int*)(&in->occopt), &in->alphamix,
                                           &in->rpnrm_cv, &in->gnrm_startmix, &in->Tel,
                                           &in->alphadiis);
  FC_FUNC_(inputs_get_geopt, INPUTS_GET_GEOPT)(in->data, in->geopt_approach,
                                               &in->ncount_cluster_x, &in->frac_fluct,
                                               &in->forcemax, &in->randdis, &in->betax,
                                               &in->history, &in->ionmov, &in->dtion,
                                               in->strtarget, &in->qmass);
  /* FC_FUNC_(inputs_get_sic, INPUTS_GET_SIC)(); */
  /* FC_FUNC_(inputs_get_tddft, INPUTS_GET_TDDFT)(); */
  
  return in;
}
void bigdft_inputs_free(BigDFT_Inputs *in)
{
  FC_FUNC_(inputs_free, INPUTS_FREE)(&in->data);
  bigdft_inputs_dispose(in);
}
void bigdft_inputs_parse_additional(BigDFT_Inputs *in, BigDFT_Atoms *atoms)
{
  int iproc = 0, dump = 0;

  FC_FUNC_(inputs_parse_add, INPUTS_PARSE_ADD)(in->data, atoms->sym, &atoms->geocode,
                                               atoms->alat, &iproc, &dump);
  FC_FUNC_(inputs_get_files, INPUTS_GET_FILES)(in->data, &in->files);

  /* FC_FUNC_(inputs_get_kpt, INPUTS_GET_KPT)(); */
}

/*******************************/
/* BigDFT_Proj data structure. */
/*******************************/
void FC_FUNC_(proj_new, PROJ_NEW)(void *nlpspd);
void FC_FUNC_(proj_free, PROJ_FREE)(void *nlpspd, f90_pointer_double *proj);
void FC_FUNC(createprojectorsarrays, CREATEPROJECTORSARRAYS)
    (int *iproc, const void *lr, double *rxyz, void *atoms,
     void *orbs, double *radii, double *cpmult, double *fpmult, const double *h1,
     const double *h2, const double *h3, void *nlpspd, f90_pointer_double *proj);
void FC_FUNC_(proj_get_dimensions, PROJ_GET_DIMENSIONS)(void *nlpspd,
                                                        guint *nproj, guint *nprojel);

static BigDFT_Proj* bigdft_proj_init()
{
  BigDFT_Proj *proj;

  proj = g_malloc(sizeof(BigDFT_Proj));
  memset(proj, 0, sizeof(BigDFT_Proj));
  F90_1D_POINTER_INIT(&proj->proj);

  return proj;
}
static void bigdft_proj_dispose(BigDFT_Proj *proj)
{
  g_free(proj);
}

BigDFT_Proj* bigdft_proj_new(const BigDFT_LocReg *glr, const BigDFT_Orbs *orbs, double frmult)
{
  BigDFT_Proj *proj;
  int iproc = 1;

  proj = bigdft_proj_init();
  FC_FUNC_(proj_new, PROJ_NEW)(&proj->nlpspd);
  FC_FUNC(createprojectorsarrays, CREATEPROJECTORSARRAYS)
    (&iproc, glr->data, glr->atoms->rxyz.data,
     glr->atoms->data, orbs->data, glr->radii, &frmult, &frmult,
     glr->h, glr->h + 1, glr->h + 2, proj->nlpspd, &proj->proj);
  FC_FUNC_(proj_get_dimensions, PROJ_GET_DIMENSIONS)(proj->nlpspd, &proj->nproj,
                                                     &proj->nprojel);

  return proj;
}
void bigdft_proj_free(BigDFT_Proj *proj)
{
  FC_FUNC_(proj_free, PROJ_FREE)(&proj->nlpspd, &proj->proj);
  bigdft_proj_dispose(proj);
}

/******************/
/* Miscellaneous. */
/******************/
void FC_FUNC(memoryestimator, MEMORYESTIMATOR)(const int *nproc, const int *idsx,
                                               const void *lr,
                                               const int *nat, const int *norb,
                                               const int *nspinor, const int *nkpt,
                                               const guint *nprojel, const int *nspin,
                                               const int *itrpmax, const int *iscf,
                                               double *peak);
double bigdft_memory_get_peak(int nproc, const BigDFT_LocReg *lr, const BigDFT_Inputs *in,
                              const BigDFT_Orbs *orbs, const BigDFT_Proj *proj)
{
  double peak;
  int nat = -1;

  FC_FUNC(memoryestimator, MEMORYESTIMATOR)(&nproc, &in->idsx, lr->data, &nat,
                                            &orbs->norb, &orbs->nspinor, &orbs->nkpts,
                                            &proj->nprojel, &orbs->nspin, &in->itrpmax,
                                            &in->iscf, &peak);
  return peak;
}
