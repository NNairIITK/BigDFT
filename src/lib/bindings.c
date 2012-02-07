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

/* Wavefunction descriptor part. */
struct f90_pointer_glr_
{
  void *glr;
  /* void *info[F90_POINTER_SIZE]; */
};

void FC_FUNC_(glr_new, GLR_NEW)(f90_pointer_glr *glr);
void FC_FUNC_(system_size, SYSTEM_SIZE)(int *iproc, void *atoms, double *rxyz,
                                        double *radii_cf, double *crmult, double *frmult,
                                        double *hx, double *hy, double *hz,
                                        f90_pointer_glr *glr, double *shift);
void FC_FUNC_(glr_get_dimensions, GLR_GET_DIMENSIONS)(void *glr, char *geocode,
                                                      int *n, int *ni);
void FC_FUNC_(glr_free, GLR_FREE)(f90_pointer_glr *glr);
void FC_FUNC_(glr_set_wave_descriptors,
             GLR_SET_WAVE_DESCRIPTORS)(int *iproc, double *hx, double *hy,
                                       double *hz, void *atoms, double *rxyz, double *radii,
                                       double *crmult, double *frmult, void *glr);

static BigDFT_Glr* bigdft_glr_init(BigDFT_Atoms *atoms, double *radii, double h[3],
                            double crmult, double frmult)
{
  BigDFT_Glr *glr;
  int iproc = 1;

  glr = g_malloc(sizeof(BigDFT_Glr));
  memset(glr, 0, sizeof(BigDFT_Glr));
  glr->data = g_malloc(sizeof(f90_pointer_glr));
  memset(glr->data, 0, sizeof(f90_pointer_glr));

  FC_FUNC_(glr_new, GLR_NEW)(glr->data);
  FC_FUNC_(system_size, SYSTEM_SIZE)(&iproc, atoms->data, atoms->rxyz.data, radii,
                                     &crmult, &frmult, h, h + 1, h + 2, glr->data->glr,
                                     atoms->shift);
  glr->h[0] = h[0];
  glr->h[1] = h[1];
  glr->h[2] = h[2];
  FC_FUNC_(glr_get_dimensions, GLR_GET_DIMENSIONS)(glr->data->glr, &glr->geocode,
                                                   (int*)glr->n, (int*)glr->ni);
  FC_FUNC_(atoms_copy_alat, ATOMS_COPY_ALAT)(atoms->data, atoms->alat,
                                             atoms->alat + 1, atoms->alat + 2);

  return glr;
}
static void bigdft_glr_dispose(BigDFT_Glr *glr)
{
  g_free(glr->data);
  g_free(glr);
}

BigDFT_Glr* bigdft_glr_new(BigDFT_Atoms *atoms, double *radii, double h[3],
                           double crmult, double frmult)
{
  BigDFT_Glr *glr;

  glr = bigdft_glr_init(atoms, radii, h, crmult, frmult);
  
  return glr;
}
BigDFT_Glr* bigdft_glr_new_with_wave_descriptors(BigDFT_Atoms *atoms, double *radii,
                                                 double h[3], double crmult, double frmult)
{
  BigDFT_Glr *glr;

  glr = bigdft_glr_init(atoms, radii, h, crmult, frmult);
  bigdft_glr_set_wave_descriptors(glr, atoms, radii, crmult, frmult);
  
  return glr;
}
void bigdft_glr_free(BigDFT_Glr *glr)
{
  FC_FUNC_(glr_free, GLR_FREE)(glr->data);
  bigdft_glr_dispose(glr);
}
void bigdft_glr_set_wave_descriptors(BigDFT_Glr *glr, BigDFT_Atoms *atoms, double *radii,
                                     double crmult, double frmult)
{
  int iproc = 1;

  FC_FUNC_(glr_set_wave_descriptors,
          GLR_SET_WAVE_DESCRIPTORS)(&iproc, glr->h, glr->h + 1, glr->h + 2,
                                    atoms->data, atoms->rxyz.data, radii,
                                    &crmult, &frmult, glr->data->glr);
}

void FC_FUNC_(fill_logrid, FILL_LOGRID)(char *geocode, int *n1, int *n2, int *n3,
                                        int *nl1, int *nu1, int *nl2, int *nu2, int *nl3, int *nu3,
                                        int *orig, int *nat, int *ntypes, int *iatype,
                                        double *rxyz, double *radii, double *mult,
                                        double *hx, double *hy, double *hz, int *grid);
guint* bigdft_fill_logrid(BigDFT_Atoms *atoms, guint n[3], double *radii,
                          double mult, double h[3])
{
  guint *grid;
  int orig = 0;

  grid = g_malloc(sizeof(guint) * (n[0] + 1) * (n[1] + 1) * (n[2] + 1));

  FC_FUNC_(fill_logrid, FILL_LOGRID)(&atoms->geocode, (int*)n, (int*)(n + 1), (int*)(n + 2),
                                     &orig, (int*)n, &orig, (int*)(n + 1), &orig, (int*)(n + 2),
                                     &orig, (int*)(&atoms->nat), (int*)(&atoms->ntypes), (int*)atoms->iatype,
                                     atoms->rxyz.data, radii, &mult, h, h + 1, h + 2, (int*)grid);

  return grid;
}


/* Wavefunction descriptor part. */
struct f90_pointer_inputs_
{
  void *in;
  /* void *info[F90_POINTER_SIZE]; */
};

void FC_FUNC_(inputs_new, INPUTS_NEW)(f90_pointer_inputs *in);
void FC_FUNC_(inputs_free, INPUTS_FREE)(f90_pointer_inputs *in);
void FC_FUNC_(inputs_set_radical, INPUTS_SET_RADICAL)(f90_pointer_inputs *in,
                                                      const gchar *rad, int *len);
void FC_FUNC_(inputs_get_dft, INPUTS_GET_DFT)(const f90_pointer_inputs *in, double *hx, double *hy, double *hz, double *crmult, double *frmult, int *ixc, int *ncharge, double *elecfield, int *nspin, int *mpol, double *gnrm_cv, int *itermax, int *nrepmax, int *ncong, int *idsx, int *dispersion, int *inputPsiId, int *output_wf_format, int *output_grid, double *rbuf, int *ncongt, int *norbv, int *nvirt, int *nplot, int *disableSym);
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
void FC_FUNC_(inputs_parse_params, INPUTS_PARSE_PARAMS)(f90_pointer_inputs *in,
                                                        int *iproc, int *dump);
void FC_FUNC_(inputs_get_files, INPUTS_GET_FILES)(const f90_pointer_inputs *in, int *files);
void FC_FUNC_(inputs_parse_add, INPUTS_PARSE_ADD)(void *in, void *atoms,
                                                  int *iproc, int *dump);

static BigDFT_Inputs* bigdft_inputs_init()
{
  BigDFT_Inputs *in;

  in = g_malloc(sizeof(BigDFT_Inputs));
  memset(in, 0, sizeof(BigDFT_Inputs));
  in->data = g_malloc(sizeof(f90_pointer_inputs));
  memset(in->data, 0, sizeof(f90_pointer_inputs));
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
  FC_FUNC_(inputs_new, INPUTS_NEW)(in->data);
  if (naming && naming[0])
    {
      len = strlen(naming);
      FC_FUNC_(inputs_set_radical, INPUTS_SET_RADICAL)(in->data->in, naming, &len);
    }
  else
    {
      len = 0;
      FC_FUNC_(inputs_set_radical, INPUTS_SET_RADICAL)(in->data->in, " ", &len);
    }
  FC_FUNC_(inputs_parse_params, INPUTS_PARSE_PARAMS)(in->data->in, &iproc, &dump);

  FC_FUNC_(inputs_get_files, INPUTS_GET_FILES)(in->data->in, &in->files);
  FC_FUNC_(inputs_get_dft, INPUTS_GET_DFT)(in->data->in, in->h, in->h + 1, in->h + 2,
                                           &in->crmult, &in->frmult, &in->ixc,
                                           &in->ncharge, in->elecfield, &in->nspin,
                                           &in->mpol, &in->gnrm_cv, &in->itermax,
                                           &in->nrepmax, &in->ncong, &in->idsx,
                                           &in->dispersion, &in->inputPsiId,
                                           &in->output_wf_format, &in->output_grid,
                                           &in->rbuf, &in->ncongt, &in->norbv, &in->nvirt,
                                           &in->nplot, &in->disableSym);
  FC_FUNC_(inputs_get_mix, INPUTS_GET_MIX)(in->data->in, &in->iscf, &in->itrpmax,
                                           &in->norbsempty, (int*)(&in->occopt), &in->alphamix,
                                           &in->rpnrm_cv, &in->gnrm_startmix, &in->Tel,
                                           &in->alphadiis);
  FC_FUNC_(inputs_get_geopt, INPUTS_GET_GEOPT)(in->data->in, in->geopt_approach,
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
  FC_FUNC_(inputs_free, INPUTS_FREE)(in->data);
  bigdft_inputs_dispose(in);
}
void bigdft_inputs_parse_additional(BigDFT_Inputs *in, BigDFT_Atoms *atoms)
{
  int iproc = 0, dump = 0;

  FC_FUNC_(inputs_parse_add, INPUTS_PARSE_ADD)(in->data->in, atoms->data, &iproc, &dump);
  FC_FUNC_(inputs_get_files, INPUTS_GET_FILES)(in->data->in, &in->files);

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
                             const BigDFT_Glr *glr, int iproc, int nproc, guint *nelec)
{
  BigDFT_Orbs *orbs;
  int nelec_, verb = 0;

  orbs = bigdft_orbs_init();
  FC_FUNC_(orbs_new, ORBS_NEW)(orbs->data);
  FC_FUNC_(read_orbital_variables, READ_ORBITAL_VARIABLES)(&iproc, &nproc, &verb, in->data->in,
                                                           atoms->data,
                                                           orbs->data->orbs, &nelec_);
  if (nelec)
    *nelec = nelec_;
  FC_FUNC_(orbs_comm, ORBS_COMM)(orbs->data->orbs, glr->data->glr, &iproc, &nproc);
  
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

BigDFT_Proj* bigdft_proj_new(const BigDFT_Atoms *atoms, const BigDFT_Glr *glr,
                             const BigDFT_Orbs *orbs, double *radii, double frmult)
{
  BigDFT_Proj *proj;
  int iproc = 1;

  proj = bigdft_proj_init();
  FC_FUNC_(proj_new, PROJ_NEW)(proj->nlpspd);
  FC_FUNC(createprojectorsarrays, CREATEPROJECTORSARRAYS)
    (&iproc, glr->data->glr, atoms->rxyz.data,
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

/**********************************/
/* BigDFT_DensPot data structure. */
/**********************************/
struct f90_pointer_rhodsc_
{
  void *rhodsc;
};
struct f90_pointer_denspotd_
{
  void *denspotd;
};

void FC_FUNC_(denspot_new, DENSPOT_NEW)(f90_pointer_denspotd *denspotd,
                                        f90_pointer_rhodsc *rhodsc);
void FC_FUNC_(denspot_free, DENSPOT_FREE)(f90_pointer_denspotd *denspotd,
                                          f90_pointer_rhodsc *rhodsc,
                                          f90_pointer_double *pot_ion,
                                          f90_pointer_double *rhopot,
                                          f90_pointer_double *rhocore,
                                          f90_pointer_double_4D *potxc);
void FC_FUNC(allocaterhopot, ALLOCATERHOPOT)(const guint *iproc, const guint *nproc,
                                             const void *glr, const double *hxh,
                                             const double *hyh, const double *hzh,
                                             const void *in, const void *atoms,
                                             const double *rxyz,
                                             const double *radii,
                                             void *denspotd, void *rhodsc,
                                             f90_pointer_double *rhopot,
                                             f90_pointer_double *pot_ion,
                                             f90_pointer_double_4D *potxc,
                                             f90_pointer_double *rhocore);

static BigDFT_DensPot* bigdft_denspot_init()
{
  BigDFT_DensPot *denspot;

  denspot = g_malloc(sizeof(BigDFT_DensPot));
  memset(denspot, 0, sizeof(BigDFT_DensPot));

  denspot->rhodsc = g_malloc(sizeof(f90_pointer_rhodsc));
  memset(denspot->rhodsc, 0, sizeof(f90_pointer_rhodsc));
  denspot->denspotd = g_malloc(sizeof(f90_pointer_denspotd));
  memset(denspot->denspotd, 0, sizeof(f90_pointer_denspotd));

  return denspot;
}
static void bigdft_denspot_dispose(BigDFT_DensPot *denspot)
{
  g_free(denspot->rhodsc);
  g_free(denspot->denspotd);
  g_free(denspot);
}

BigDFT_DensPot* bigdft_denspot_new (const BigDFT_Atoms *atoms, const BigDFT_Glr *glr,
                                    const BigDFT_Inputs *in, const double *radii,
                                    guint iproc, guint nproc)
{
  BigDFT_DensPot *denspot;
  double hh[3];

  denspot = bigdft_denspot_init();
  FC_FUNC_(denspot_new, DENSPOT_NEW)(denspot->denspotd, denspot->rhodsc);
  hh[0] = glr->h[0] * 0.5;
  hh[1] = glr->h[1] * 0.5;
  hh[2] = glr->h[2] * 0.5;
  FC_FUNC(allocaterhopot, ALLOCATERHOPOT)(&iproc, &nproc, glr->data->glr,
                                          hh, hh + 1, hh + 2, in->data->in,
                                          atoms->data, atoms->rxyz.data,
                                          radii, denspot->denspotd->denspotd,
                                          denspot->rhodsc->rhodsc, &denspot->rhopot,
                                          &denspot->pot_ion, &denspot->potxc,
                                          &denspot->rhocore);

  return denspot;
}
void bigdft_denspot_free(BigDFT_DensPot *denspotd)
{
  FC_FUNC_(denspot_free, DENSPOT_FREE)(denspotd->denspotd, denspotd->rhodsc,
                                       &denspotd->pot_ion, &denspotd->rhopot,
                                       &denspotd->rhocore, &denspotd->potxc);
  bigdft_denspot_dispose(denspotd);
}

/******************/
/* Miscellaneous. */
/******************/
void FC_FUNC(memoryestimator, MEMORYESTIMATOR)(int *nproc, int *idsx, void *lr, int *nat,
                                               int *norb, int *nspinor, int *nkpt,
                                               guint *nprojel, int *nspin, int *itrpmax,
                                               int *iscf, double *peak);
double bigdft_memory_peak(int nproc, BigDFT_Glr *lr, BigDFT_Inputs *in,
                          BigDFT_Orbs *orbs, BigDFT_Proj *proj)
{
  double peak;
  int nat = -1;

  FC_FUNC(memoryestimator, MEMORYESTIMATOR)(&nproc, &in->idsx, lr->data->glr, &nat,
                                            &orbs->norb, &orbs->nspinor, &orbs->nkpts,
                                            &proj->nprojel, &orbs->nspin, &in->itrpmax,
                                            &in->iscf, &peak);
  return peak;
}

void FC_FUNC(createkernel, CREATEKERNEL)(const guint *iproc, const guint *nproc,
                                         const gchar *geocode,
                                         const guint *n1i, const guint *n2i, const guint *n3i,
                                         const double *hxh, const double *hyh,
                                         const double *hzh, const guint *ndegree_ip,
                                         f90_pointer_double *pkernel, const guint *verb);
f90_pointer_double* bigdft_psolver_create_kernel(const BigDFT_Glr *glr, guint iproc,
                                                 guint nproc)
{
  f90_pointer_double *pkernel;
  guint ndegree_ip = 16, verb = 0;
  double hh[3];

  pkernel = g_malloc(sizeof(f90_pointer_double));
  F90_1D_POINTER_INIT(pkernel);

  hh[0] = glr->h[0] * 0.5;
  hh[1] = glr->h[1] * 0.5;
  hh[2] = glr->h[2] * 0.5;
  FC_FUNC(createkernel, CREATEKERNEL)(&iproc, &nproc, &glr->geocode, glr->ni,
                                      glr->ni + 1, glr->ni + 2, hh, hh + 1, hh + 2,
                                      &ndegree_ip, pkernel, &verb);
  
  return pkernel;
}
void bigdft_psolver_free_kernel(f90_pointer_double *pkernel)
{
  FC_FUNC_(deallocate_double_1d, DEALLOCATE_DOUBLE_1D)(pkernel);
  g_free(pkernel);
}
