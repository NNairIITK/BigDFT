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

void FC_FUNC_(inquire_address1, INQUIRE_ADDRESS1)(double *add, void *pt)
{
  double *val = (double*)(&pt);
  *add = val[0];
}
void FC_FUNC_(inquire_address2, INQUIRE_ADDRESS2)(double *add, void *pt)
{
  double *val = (double*)(&pt);
  *add = val[0];
}


void FC_FUNC_(inputs_new, INPUTS_NEW)(void *in);
void FC_FUNC_(inputs_free, INPUTS_FREE)(void *in);
void FC_FUNC_(inputs_set_radical, INPUTS_SET_RADICAL)(void *in, int *nproc,
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
  int iproc = 0, len, dump = 0, nproc = 1;

  in = bigdft_inputs_init();
  FC_FUNC_(inputs_new, INPUTS_NEW)(&in->data);
  if (naming && naming[0])
    {
      len = strlen(naming);
      FC_FUNC_(inputs_set_radical, INPUTS_SET_RADICAL)(in->data, &nproc, naming, &len);
    }
  else
    {
      len = 0;
      FC_FUNC_(inputs_set_radical, INPUTS_SET_RADICAL)(in->data, &nproc, " ", &len);
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

#ifdef HAVE_GLIB
G_DEFINE_TYPE(BigDFT_Proj, bigdft_proj, G_TYPE_OBJECT)

static void bigdft_proj_dispose(GObject *proj);
static void bigdft_proj_finalize(GObject *proj);

static void bigdft_proj_class_init(BigDFT_ProjClass *klass)
{
  /* Connect the overloading methods. */
  G_OBJECT_CLASS(klass)->dispose      = bigdft_proj_dispose;
  G_OBJECT_CLASS(klass)->finalize     = bigdft_proj_finalize;
  /* G_OBJECT_CLASS(klass)->set_property = visu_data_set_property; */
  /* G_OBJECT_CLASS(klass)->get_property = visu_data_get_property; */
}
#endif

static void bigdft_proj_init(BigDFT_Proj *obj)
{
#ifdef HAVE_GLIB
  memset((void*)((char*)obj + sizeof(GObject)), 0, sizeof(BigDFT_Proj) - sizeof(GObject));
#else
  memset(obj, 0, sizeof(BigDFT_Proj));
#endif
  FC_FUNC_(proj_new, PROJ_NEW)(&obj->nlpspd);
  F90_1D_POINTER_INIT(&obj->proj);
}
static void bigdft_proj_dispose(GObject *obj)
{
#ifdef HAVE_GLIB
  BigDFT_Proj *proj = BIGDFT_PROJ(obj);

  if (proj->dispose_has_run)
    return;
  proj->dispose_has_run = TRUE;

  /* Chain up to the parent class */
  G_OBJECT_CLASS(bigdft_proj_parent_class)->dispose(obj);
#endif
}
static void bigdft_proj_finalize(GObject *obj)
{
  BigDFT_Proj *proj = BIGDFT_PROJ(obj);

  FC_FUNC_(proj_free, PROJ_FREE)(&proj->nlpspd, &proj->proj);

#ifdef HAVE_GLIB
  G_OBJECT_CLASS(bigdft_proj_parent_class)->finalize(obj);
  /* g_debug("Freeing proj object %p done.\n", obj); */
#endif
}
BigDFT_Proj* bigdft_proj_new(const BigDFT_LocReg *glr, const BigDFT_Orbs *orbs, double frmult)
{
  BigDFT_Proj *proj;
  int iproc = 1;

#ifdef HAVE_GLIB
  proj = BIGDFT_PROJ(g_object_new(BIGDFT_PROJ_TYPE, NULL));
#else
  proj = g_malloc(sizeof(BigDFT_Proj));
  bigdft_proj_init(proj);
#endif

  FC_FUNC(createprojectorsarrays, CREATEPROJECTORSARRAYS)
    (&iproc, glr->data, glr->parent.rxyz.data,
     glr->parent.data, orbs->data, glr->radii, &frmult, &frmult,
     glr->h, glr->h + 1, glr->h + 2, proj->nlpspd, &proj->proj);
  FC_FUNC_(proj_get_dimensions, PROJ_GET_DIMENSIONS)(proj->nlpspd, &proj->nproj,
                                                     &proj->nprojel);

  return proj;
}
void bigdft_proj_free(BigDFT_Proj *proj)
{
#ifdef HAVE_GLIB
  g_object_unref(G_OBJECT(proj));
#else
  bigdft_proj_finalize(proj);
  g_free(proj);
#endif
}

/********************************/
/* BigDFT_Energs data structure */
/********************************/
#ifdef HAVE_GLIB
enum {
  EKS_READY_SIGNAL,
  LAST_SIGNAL
};

G_DEFINE_TYPE(BigDFT_Energs, bigdft_energs, G_TYPE_OBJECT)

static guint bigdft_energs_signals[LAST_SIGNAL] = { 0 };

static void bigdft_energs_dispose(GObject *energs);
static void bigdft_energs_finalize(GObject *energs);

static void bigdft_energs_class_init(BigDFT_EnergsClass *klass)
{
  /* Connect the overloading methods. */
  G_OBJECT_CLASS(klass)->dispose      = bigdft_energs_dispose;
  G_OBJECT_CLASS(klass)->finalize     = bigdft_energs_finalize;
  /* G_OBJECT_CLASS(klass)->set_property = visu_data_set_property; */
  /* G_OBJECT_CLASS(klass)->get_property = visu_data_get_property; */

  bigdft_energs_signals[EKS_READY_SIGNAL] =
    g_signal_new("eks-ready", G_TYPE_FROM_CLASS(klass),
                 G_SIGNAL_RUN_LAST | G_SIGNAL_NO_RECURSE | G_SIGNAL_NO_HOOKS,
		 0, NULL, NULL, g_cclosure_marshal_VOID__UINT,
                 G_TYPE_NONE, 1, G_TYPE_UINT, NULL);
}
#endif

static void bigdft_energs_init(BigDFT_Energs *obj)
{
  double self;

#ifdef HAVE_GLIB
  memset((void*)((char*)obj + sizeof(GObject)), 0, sizeof(BigDFT_Energs) - sizeof(GObject));
#else
  memset(obj, 0, sizeof(BigDFT_Energs));
#endif
  self = *((double*)&obj);
  FC_FUNC_(energs_new, ENERGS_NEW)(&self, &obj->data);
}
static void bigdft_energs_dispose(GObject *obj)
{
#ifdef HAVE_GLIB
  BigDFT_Energs *energs = BIGDFT_ENERGS(obj);

  if (energs->dispose_has_run)
    return;
  energs->dispose_has_run = TRUE;

  /* Chain up to the parent class */
  G_OBJECT_CLASS(bigdft_energs_parent_class)->dispose(obj);
#endif
}
static void bigdft_energs_finalize(GObject *obj)
{
  BigDFT_Energs *energs = BIGDFT_ENERGS(obj);

  FC_FUNC_(energs_free, ENERGS_FREE)(&energs->data);

#ifdef HAVE_GLIB
  G_OBJECT_CLASS(bigdft_energs_parent_class)->finalize(obj);
#endif
}
BigDFT_Energs* bigdft_energs_new()
{
  BigDFT_Energs *energs;

#ifdef HAVE_GLIB
  energs = BIGDFT_ENERGS(g_object_new(BIGDFT_ENERGS_TYPE, NULL));
#else
  energs = g_malloc(sizeof(BigDFT_Energs));
  bigdft_energs_init(energs);
#endif

  return energs;
}
void bigdft_energs_free(BigDFT_Energs *energs)
{
#ifdef HAVE_GLIB
  g_object_unref(G_OBJECT(energs));
#else
  bigdft_energs_finalize(energs);
  g_free(energs);
#endif
}
void FC_FUNC_(energs_emit, ENERGS_EMIT)(BigDFT_Energs **obj, guint *istep,
                                        BigDFT_EnergsKind *kind)
{
  BigDFT_Energs *energs = *obj;

#ifdef HAVE_GLIB
  FC_FUNC_(energs_copy_data, ENERGS_COPY_DATA)
    (energs->data, &energs->eh, &energs->exc,
     &energs->evxc, &energs->eion, &energs->edisp,
     &energs->ekin, &energs->epot, &energs->eproj,
     &energs->eexctX, &energs->ebs, &energs->eKS,
     &energs->trH, &energs->evsum, &energs->evsic);
  switch (*kind)
    {
    case BIGDFT_E_KS:
      g_signal_emit(G_OBJECT(energs), bigdft_energs_signals[EKS_READY_SIGNAL],
                    0 /* details */, *istep, NULL);
      break;
    default:
      break;
    }
#endif  
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

