#include <config.h>

#ifdef HAVE_GLIB
#include <glib-object.h>
#endif

#include "bigdft.h"
#include "bindings.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

void FC_FUNC_(read_wave_to_isf, READ_WAVE_TO_ISF)
     (int *lstat, const char* filename, int *ln, int *iorbp,
      double *hx, double *hy, double *hz,
      int *n1, int *n2, int *n3, int *nspinor, f90_pointer_double_4D *psiscf);
void FC_FUNC_(free_wave_to_isf, FREE_WAVE_TO_ISF)(f90_pointer_double_4D *psiscf);

void FC_FUNC_(read_wave_descr, READ_WAVE_DESCR)
     (int *lstat, const char* filename, int *ln, int *norbu,
      int *norbd, int *iorb, int *ispin, int *nkpt, int *ikpt, int *nspinor, int *ispinor);

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



f90_pointer_double_4D* bigdft_read_wave_to_isf(const char *filename, int iorbp,
                                               double h[3], int n[3], int *nspinor)
{
  int ln, lstat;
  f90_pointer_double_4D *psiscf;

  psiscf = g_malloc(sizeof(f90_pointer_double_4D));
  F90_4D_POINTER_INIT(psiscf);
  
  ln = strlen(filename);
  FC_FUNC_(read_wave_to_isf, READ_WAVE_TO_ISF)
    (&lstat, filename, &ln, &iorbp, h, h + 1, h + 2, n, n + 1, n + 2, nspinor, psiscf);
  if (!lstat)
    {
      g_free(psiscf);
      psiscf = (f90_pointer_double_4D*)0;
    }

  DBG_MEM(psiscf, f90_pointer_double_4D);

  return psiscf;
}
void bigdft_free_wave_to_isf(f90_pointer_double_4D *psiscf)
{
  FC_FUNC_(free_wave_to_isf, FREE_WAVE_TO_ISF)(psiscf);
  g_free(psiscf);
}

gboolean bigdft_read_wave_descr(const char *filename, int *norbu,
                                int *norbd, int *nkpt, int *nspinor,
                                int *iorb, int *ispin, int *ikpt, int *ispinor)
{
  int ln, lstat, norbu_, norbd_, nkpt_, nspinor_;
  int iorb_, ispin_, ikpt_, ispinor_;
  
  ln = strlen(filename);
  FC_FUNC_(read_wave_descr, READ_WAVE_DESCR)
    (&lstat, filename, &ln, &norbu_, &norbd_, &iorb_, &ispin_,
     &nkpt_, &ikpt_, &nspinor_, &ispinor_);
  if (!lstat)
    return FALSE;

  if (norbu)   *norbu   = norbu_;
  if (norbd)   *norbd   = norbd_;
  if (nkpt)    *nkpt    = nkpt_;
  if (nspinor) *nspinor = nspinor_;

  if (iorb)    *iorb    = iorb_;
  if (ispin)   *ispin   = ispin_;
  if (ikpt)    *ikpt    = ikpt_;
  if (ispinor) *ispinor = ispinor_;

  return TRUE;
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

/* Orbital descriptor part. */
struct f90_pointer_orbs_
{
  void *orbs;
  /* void *info[F90_POINTER_SIZE]; */
};

void FC_FUNC_(orbs_new, ORBS_NEW)(f90_pointer_orbs *orbs);
void FC_FUNC_(orbs_free, ORBS_FREE)(f90_pointer_orbs *orbs);
void FC_FUNC_(read_orbital_variables, READ_ORBITAL_VARIABLES)(int *iproc, int *nproc,
                                                              int *verb, void *in, void *atoms,
                                                              void *orbs, int *nelec);
void FC_FUNC_(orbs_comm, ORBS_COMM)(void *orbs, void *glr, const int *iproc, const int *nproc);
void FC_FUNC_(orbs_get_dimensions, ORBS_GET_DIMENSIONS)(void *orbs, int *norb,
                                                        int *norbp, int *norbu,
                                                        int *norbd, int *nspin,
                                                        int *nspinor, int *npsidim,
                                                        int *nkpts, int *nkptsp,
                                                        int *isorb, int *iskpts);

static BigDFT_Orbs* bigdft_orbs_init()
{
  BigDFT_Orbs *orbs;

  orbs = g_malloc(sizeof(BigDFT_Orbs));
  memset(orbs, 0, sizeof(BigDFT_Orbs));
  orbs->data = g_malloc(sizeof(f90_pointer_orbs));
  memset(orbs->data, 0, sizeof(f90_pointer_orbs));

  return orbs;
}
static void bigdft_orbs_dispose(BigDFT_Orbs *orbs)
{
  g_free(orbs->data);
  g_free(orbs);
}

BigDFT_Orbs* bigdft_orbs_new(const BigDFT_Atoms *atoms, const BigDFT_Inputs *in,
                             const BigDFT_LocReg *glr, int iproc, int nproc, guint *nelec)
{
  BigDFT_Orbs *orbs;
  int nelec_, verb = 0;

  orbs = bigdft_orbs_init();
  FC_FUNC_(orbs_new, ORBS_NEW)(orbs->data);
  FC_FUNC_(read_orbital_variables, READ_ORBITAL_VARIABLES)(&iproc, &nproc, &verb, in->data,
                                                           atoms->data,
                                                           orbs->data->orbs, &nelec_);
  if (nelec)
    *nelec = nelec_;
  FC_FUNC_(orbs_comm, ORBS_COMM)(orbs->data->orbs, glr->data, &iproc, &nproc);
  
  FC_FUNC_(orbs_get_dimensions, ORBS_GET_DIMENSIONS)(orbs->data->orbs, &orbs->norb,
                                                     &orbs->norbp, &orbs->norbu,
                                                     &orbs->norbd, &orbs->nspin,
                                                     &orbs->nspinor, &orbs->npsidim,
                                                     &orbs->nkpts, &orbs->nkptsp,
                                                     &orbs->isorb, &orbs->iskpts);
  
  return orbs;
}
void bigdft_orbs_free(BigDFT_Orbs *orbs)
{
  FC_FUNC_(orbs_free, ORBS_FREE)(orbs->data);
  bigdft_orbs_dispose(orbs);
}

/*******************************/
/* BigDFT_Proj data structure. */
/*******************************/
struct f90_pointer_nlpspd_
{
  void *proj;
};

void FC_FUNC_(proj_new, PROJ_NEW)(f90_pointer_nlpspd *nlpspd);
void FC_FUNC_(proj_free, PROJ_FREE)(f90_pointer_nlpspd *nlpspd, f90_pointer_double *proj);
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
  proj->nlpspd = g_malloc(sizeof(f90_pointer_nlpspd));
  memset(proj->nlpspd, 0, sizeof(f90_pointer_nlpspd));
  F90_1D_POINTER_INIT(&proj->proj);

  return proj;
}
static void bigdft_proj_dispose(BigDFT_Proj *proj)
{
  g_free(proj->nlpspd);
  g_free(proj);
}

BigDFT_Proj* bigdft_proj_new(const BigDFT_Atoms *atoms, const BigDFT_LocReg *glr,
                             const BigDFT_Orbs *orbs, double *radii, double frmult)
{
  BigDFT_Proj *proj;
  int iproc = 1;

  proj = bigdft_proj_init();
  FC_FUNC_(proj_new, PROJ_NEW)(proj->nlpspd);
  FC_FUNC(createprojectorsarrays, CREATEPROJECTORSARRAYS)
    (&iproc, glr->data, atoms->rxyz.data,
     atoms->data, orbs->data->orbs, radii, &frmult, &frmult,
     glr->h, glr->h + 1, glr->h + 2, proj->nlpspd->proj, &proj->proj);
  FC_FUNC_(proj_get_dimensions, PROJ_GET_DIMENSIONS)(proj->nlpspd->proj, &proj->nproj,
                                                     &proj->nprojel);

  return proj;
}
void bigdft_proj_free(BigDFT_Proj *proj)
{
  FC_FUNC_(proj_free, PROJ_FREE)(proj->nlpspd, &proj->proj);
  bigdft_proj_dispose(proj);
}

/******************/
/* Miscellaneous. */
/******************/
void FC_FUNC(memoryestimator, MEMORYESTIMATOR)(int *nproc, int *idsx, void *lr, int *nat,
                                               int *norb, int *nspinor, int *nkpt,
                                               guint *nprojel, int *nspin, int *itrpmax,
                                               int *iscf, double *peak);
double bigdft_memory_peak(int nproc, BigDFT_LocReg *lr, BigDFT_Inputs *in,
                          BigDFT_Orbs *orbs, BigDFT_Proj *proj)
{
  double peak;
  int nat = -1;

  FC_FUNC(memoryestimator, MEMORYESTIMATOR)(&nproc, &in->idsx, lr->data, &nat,
                                            &orbs->norb, &orbs->nspinor, &orbs->nkpts,
                                            &proj->nprojel, &orbs->nspin, &in->itrpmax,
                                            &in->iscf, &peak);
  return peak;
}
