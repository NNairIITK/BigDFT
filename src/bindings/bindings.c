#include <config.h>

#ifdef HAVE_GLIB
#include <glib-object.h>
#include <gio/gio.h>
#endif

#include "bigdft.h"
#include "bindings.h"
#include "bindings_api.h"

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

static BigDFT_Inputs* bigdft_inputs_init()
{
  BigDFT_Inputs *in;

  in = g_malloc(sizeof(BigDFT_Inputs));
  memset(in, 0, sizeof(BigDFT_Inputs));
  in->data = (void*)0;
  in->refCount = 1;
  F90_1D_POINTER_INIT(&in->qmass);

  return in;
}
static void bigdft_inputs_dispose(BigDFT_Inputs *in)
{
  g_free(in->data);
  g_free(in);
}
/**
 * bigdft_inputs_new:
 * @naming: (allow-none): a naming scheme, or none.
 *
 * Create a new #BigDFT_Inputs structure.
 * 
 * Returns: (transfer full): a new structure.
 */
BigDFT_Inputs* bigdft_inputs_new(const gchar *naming)
{
  BigDFT_Inputs *in;
  int iproc = 0, len, nproc = 1;
  int dump = 0;

  in = bigdft_inputs_init();
  FC_FUNC_(inputs_new, INPUTS_NEW)(&in->data);
  if (naming && naming[0])
    {
      len = strlen(naming);
      FC_FUNC_(inputs_set_radical, INPUTS_SET_RADICAL)(in->data, &nproc, naming, &len, len);
    }
  else
    {
      len = 0;
      FC_FUNC_(inputs_set_radical, INPUTS_SET_RADICAL)(in->data, &nproc, " ", &len, 1);
    }
  FC_FUNC_(inputs_parse_params, INPUTS_PARSE_PARAMS)(in->data, &iproc, &dump);

  FC_FUNC_(inputs_get_files, INPUTS_GET_FILES)(in->data, &in->files);
  FC_FUNC_(inputs_get_dft, INPUTS_GET_DFT)(in->data, in->h, in->h + 1, in->h + 2,
                                           &in->crmult, &in->frmult, &in->ixc,
                                           &in->ncharge, in->elecfield, &in->nspin,
                                           &in->mpol, &in->gnrm_cv, (int*)&in->itermax,
                                           (int*)&in->nrepmax, &in->ncong, (int*)&in->idsx,
                                           &in->dispersion, &in->inputPsiId,
                                           &in->output_wf_format, &in->output_grid,
                                           &in->rbuf, &in->ncongt, &in->norbv, &in->nvirt,
                                           &in->nplot, &in->disableSym);
  FC_FUNC_(inputs_get_mix, INPUTS_GET_MIX)(in->data, (int*)&in->iscf, (int*)&in->itrpmax,
                                           (int*)&in->norbsempty, (int*)(&in->occopt),
                                           &in->alphamix,
                                           &in->rpnrm_cv, &in->gnrm_startmix, &in->Tel,
                                           &in->alphadiis);
  FC_FUNC_(inputs_get_geopt, INPUTS_GET_GEOPT)(in->data, in->geopt_approach,
                                               &in->ncount_cluster_x, &in->frac_fluct,
                                               &in->forcemax, &in->randdis, &in->betax,
                                               &in->history, &in->ionmov, &in->dtion,
                                               in->strtarget, &in->qmass, 10);
  /* FC_FUNC_(inputs_get_sic, INPUTS_GET_SIC)(); */
  /* FC_FUNC_(inputs_get_tddft, INPUTS_GET_TDDFT)(); */
  FC_FUNC_(inputs_get_perf, INPUTS_GET_PERF)(in->data, (int*)&in->linear);
  
  return in;
}
void bigdft_inputs_free(BigDFT_Inputs *in)
{
  if (in->data)
    FC_FUNC_(inputs_free, INPUTS_FREE)(&in->data);
  bigdft_inputs_dispose(in);
}
BigDFT_Inputs* bigdft_inputs_ref(BigDFT_Inputs *in)
{
  in->refCount += 1;
  return in;
}
void bigdft_inputs_unref(BigDFT_Inputs *in)
{
  in->refCount -= 1;
  if (!in->refCount)
    bigdft_inputs_free(in);
}
#ifdef GLIB_MAJOR_VERSION
GType bigdft_inputs_get_type(void)
{
  static GType g_define_type_id = 0;

  if (g_define_type_id == 0)
    g_define_type_id =
      g_boxed_type_register_static("BigDFT_Inputs", 
                                   (GBoxedCopyFunc)bigdft_inputs_ref,
                                   (GBoxedFreeFunc)bigdft_inputs_unref);
  return g_define_type_id;
}
#endif
void bigdft_inputs_parse_additional(BigDFT_Inputs *in, BigDFT_Atoms *atoms)
{
  int iproc = 0, dump = 0;

  FC_FUNC_(inputs_parse_add, INPUTS_PARSE_ADD)(in->data, atoms->data, &iproc, &dump);
  FC_FUNC_(inputs_get_files, INPUTS_GET_FILES)(in->data, &in->files);

  /* FC_FUNC_(inputs_get_kpt, INPUTS_GET_KPT)(); */
}

/*******************************/
/* BigDFT_Proj data structure. */
/*******************************/
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
BigDFT_Proj* bigdft_proj_new(const BigDFT_Locreg *glr, const BigDFT_Orbs *orbs, double frmult)
{
  BigDFT_Proj *proj;
  int iproc = 1;
  _gaussian_basis *gauss;

#ifdef HAVE_GLIB
  proj = BIGDFT_PROJ(g_object_new(BIGDFT_PROJ_TYPE, NULL));
#else
  proj = g_malloc(sizeof(BigDFT_Proj));
  bigdft_proj_init(proj);
#endif

  FC_FUNC(createprojectorsarrays, CREATEPROJECTORSARRAYS)
    (&iproc, glr->data, glr->parent.rxyz.data,
     glr->parent.data, orbs->data, &g_array_index(glr->radii, double, 0), &frmult, &frmult,
     glr->h, glr->h + 1, glr->h + 2, proj->nlpspd, gauss, &proj->proj);
  FC_FUNC_(proj_get_dimensions, PROJ_GET_DIMENSIONS)(proj->nlpspd, (int*)&proj->nproj,
                                                     (int*)&proj->nprojel);

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
#ifdef HAVE_GLIB
  memset((void*)((char*)obj + sizeof(GObject)), 0, sizeof(BigDFT_Energs) - sizeof(GObject));
#else
  memset(obj, 0, sizeof(BigDFT_Energs));
#endif
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

  if (energs->data)
    FC_FUNC_(energs_free, ENERGS_FREE)(&energs->data);

#ifdef HAVE_GLIB
  G_OBJECT_CLASS(bigdft_energs_parent_class)->finalize(obj);
#endif
}
BigDFT_Energs* bigdft_energs_new()
{
  BigDFT_Energs *energs;
  long self;

#ifdef HAVE_GLIB
  energs = BIGDFT_ENERGS(g_object_new(BIGDFT_ENERGS_TYPE, NULL));
#else
  energs = g_malloc(sizeof(BigDFT_Energs));
  bigdft_energs_init(energs);
#endif
  self = *((long*)&energs);
  FC_FUNC_(energs_new, ENERGS_NEW)(&self, &energs->data);

  return energs;
}
void FC_FUNC_(energs_new_wrapper, ENERGS_NEW_WRAPPER)(double *self, void *obj)
{
  BigDFT_Energs *energs;

  energs = bigdft_energs_new_from_fortran(obj);
  *self = *((double*)&energs);
}
BigDFT_Energs* bigdft_energs_new_from_fortran(void *obj)
{
  BigDFT_Energs *energs;

#ifdef HAVE_GLIB
  energs = BIGDFT_ENERGS(g_object_new(BIGDFT_ENERGS_TYPE, NULL));
#else
  energs = g_malloc(sizeof(BigDFT_Energs));
  bigdft_energs_init(energs);
#endif
  energs->data = obj;

  return energs;
}
void FC_FUNC_(energs_free_wrapper, ENERGS_FREE_WRAPPER)(gpointer *obj)
{
  BigDFT_Energs *energs = BIGDFT_ENERGS(*obj);

  energs->data = (gpointer)0;
  bigdft_energs_free(energs);
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
                                        BigDFT_EnergsIds *kind)
{
  BigDFT_Energs *energs = BIGDFT_ENERGS(*obj);

  FC_FUNC_(energs_copy_data, ENERGS_COPY_DATA)
    (energs->data, &energs->eh, &energs->exc,
     &energs->evxc, &energs->eion, &energs->edisp,
     &energs->ekin, &energs->epot, &energs->eproj,
     &energs->eexctX, &energs->ebs, &energs->eKS,
     &energs->trH, &energs->evsum, &energs->evsic);
  bigdft_energs_emit(*obj, *istep, *kind);
}
void bigdft_energs_emit(BigDFT_Energs *energs, guint istep, BigDFT_EnergsIds kind)
{
#ifdef HAVE_GLIB
  switch (kind)
    {
    case BIGDFT_ENERGS_EKS:
      g_signal_emit(G_OBJECT(energs), bigdft_energs_signals[EKS_READY_SIGNAL],
                    0 /* details */, istep, NULL);
      break;
    default:
      break;
    }
#endif  
}

/*********************************/
/* BigDFT_Restart data structure */
/*********************************/
#ifdef HAVE_GLIB
G_DEFINE_TYPE(BigDFT_Restart, bigdft_restart, G_TYPE_OBJECT)

static void bigdft_restart_dispose(GObject *restart);
static void bigdft_restart_finalize(GObject *restart);

static void bigdft_restart_class_init(BigDFT_RestartClass *klass)
{
  /* Connect the overloading methods. */
  G_OBJECT_CLASS(klass)->dispose      = bigdft_restart_dispose;
  G_OBJECT_CLASS(klass)->finalize     = bigdft_restart_finalize;
  /* G_OBJECT_CLASS(klass)->set_property = visu_data_set_property; */
  /* G_OBJECT_CLASS(klass)->get_property = visu_data_get_property; */
}
#endif

static void bigdft_restart_init(BigDFT_Restart *obj)
{
#ifdef HAVE_GLIB
  memset((void*)((char*)obj + sizeof(GObject)), 0, sizeof(BigDFT_Restart) - sizeof(GObject));
#else
  memset(obj, 0, sizeof(BigDFT_Restart));
#endif
}
static void bigdft_restart_dispose(GObject *obj)
{
#ifdef HAVE_GLIB
  BigDFT_Restart *restart = BIGDFT_RESTART(obj);

  if (restart->dispose_has_run)
    return;
  restart->dispose_has_run = TRUE;

  /* Chain up to the parent class */
  G_OBJECT_CLASS(bigdft_restart_parent_class)->dispose(obj);
#endif
}
static void bigdft_restart_finalize(GObject *obj)
{
  BigDFT_Restart *restart = BIGDFT_RESTART(obj);

  if (restart->data)
    FC_FUNC_(rst_free, RST_FREE)(&restart->data);
  if (restart->in)
    bigdft_inputs_unref(restart->in);

#ifdef HAVE_GLIB
  G_OBJECT_CLASS(bigdft_restart_parent_class)->finalize(obj);
#endif
}
BigDFT_Restart* bigdft_restart_new(BigDFT_Atoms *atoms, BigDFT_Inputs *in, guint iproc)
{
  BigDFT_Restart *restart;
  long self;

#ifdef HAVE_GLIB
  restart = BIGDFT_RESTART(g_object_new(BIGDFT_RESTART_TYPE, NULL));
#else
  restart = g_malloc(sizeof(BigDFT_Restart));
  bigdft_restart_init(restart);
#endif
  self = *((long*)&restart);
  FC_FUNC_(rst_new, RST_NEW)(&self, &restart->data);
  FC_FUNC_(rst_init, RST_INIT)(restart->data, &iproc, atoms->data, in->data);
  restart->inputPsiId = BIGDFT_RESTART_LCAO;
  restart->in = in;
  bigdft_inputs_ref(in);

  return restart;
}
void FC_FUNC_(restart_new_wrapper, RESTART_NEW_WRAPPER)(double *self, void *obj)
{
  BigDFT_Restart *restart;

  restart = bigdft_restart_new_from_fortran(obj);
  *self = *((double*)&restart);
}
BigDFT_Restart* bigdft_restart_new_from_fortran(void *obj)
{
  BigDFT_Restart *restart;

#ifdef HAVE_GLIB
  restart = BIGDFT_RESTART(g_object_new(BIGDFT_RESTART_TYPE, NULL));
#else
  restart = g_malloc(sizeof(BigDFT_Restart));
  bigdft_restart_init(restart);
#endif
  restart->data = obj;

  return restart;
}
void FC_FUNC_(restart_free_wrapper, RESTART_FREE_WRAPPER)(gpointer *obj)
{
  BigDFT_Restart *restart = BIGDFT_RESTART(*obj);

  restart->data = (gpointer)0;
  bigdft_restart_free(restart);
}
void bigdft_restart_free(BigDFT_Restart *restart)
{
#ifdef HAVE_GLIB
  g_object_unref(G_OBJECT(restart));
#else
  bigdft_restart_finalize(restart);
  g_free(restart);
#endif
}
void bigdft_restart_set(BigDFT_Restart *restart, BigDFT_RestartIds id)
{
  int inputPsiId[] = {0, 1};

  restart->inputPsiId = id;
  FC_FUNC_(inputs_set_restart, INPUTS_SET_RESTART)(restart->in->data, inputPsiId + id);
}

/*********************************/
/* BigDFT_OptLoop data structure */
/*********************************/
#ifdef HAVE_GLIB
enum {
  SYNC_FORTRAN,
  ITER_HAMILTONIAN,
  ITER_SUBSPACE,
  ITER_WAVEFUNCTIONS,
  DONE_HAMILTONIAN,
  DONE_SUBSPACE,
  DONE_WAVEFUNCTIONS,
  LAST_SIGNAL_OPTLOOP
};

G_DEFINE_TYPE(BigDFT_OptLoop, bigdft_optloop, G_TYPE_OBJECT)

static guint bigdft_optloop_signals[LAST_SIGNAL_OPTLOOP] = { 0 };

static void bigdft_optloop_dispose(GObject *optloop);
static void bigdft_optloop_finalize(GObject *optloop);

static void bigdft_optloop_class_init(BigDFT_OptLoopClass *klass)
{
  /* Connect the overloading methods. */
  G_OBJECT_CLASS(klass)->dispose      = bigdft_optloop_dispose;
  G_OBJECT_CLASS(klass)->finalize     = bigdft_optloop_finalize;
  /* G_OBJECT_CLASS(klass)->set_property = visu_data_set_property; */
  /* G_OBJECT_CLASS(klass)->get_property = visu_data_get_property; */

  bigdft_optloop_signals[SYNC_FORTRAN] =
    g_signal_new("sync-fortran", G_TYPE_FROM_CLASS(klass),
                 G_SIGNAL_RUN_LAST | G_SIGNAL_NO_RECURSE | G_SIGNAL_NO_HOOKS,
        	 0, NULL, NULL, g_cclosure_marshal_VOID__VOID,
                 G_TYPE_NONE, 0, NULL);

  bigdft_optloop_signals[ITER_HAMILTONIAN] =
    g_signal_new("iter-hamiltonian", G_TYPE_FROM_CLASS(klass),
                 G_SIGNAL_RUN_LAST | G_SIGNAL_NO_RECURSE | G_SIGNAL_NO_HOOKS,
        	 0, NULL, NULL, g_cclosure_marshal_VOID__OBJECT,
                 G_TYPE_NONE, 1, G_TYPE_OBJECT, NULL);

  bigdft_optloop_signals[ITER_SUBSPACE] =
    g_signal_new("iter-subspace", G_TYPE_FROM_CLASS(klass),
                 G_SIGNAL_RUN_LAST | G_SIGNAL_NO_RECURSE | G_SIGNAL_NO_HOOKS,
        	 0, NULL, NULL, g_cclosure_marshal_VOID__OBJECT,
                 G_TYPE_NONE, 1, G_TYPE_OBJECT, NULL);

  bigdft_optloop_signals[ITER_WAVEFUNCTIONS] =
    g_signal_new("iter-wavefunctions", G_TYPE_FROM_CLASS(klass),
                 G_SIGNAL_RUN_LAST | G_SIGNAL_NO_RECURSE | G_SIGNAL_NO_HOOKS,
        	 0, NULL, NULL, g_cclosure_marshal_VOID__OBJECT,
                 G_TYPE_NONE, 1, G_TYPE_OBJECT, NULL);

  bigdft_optloop_signals[DONE_HAMILTONIAN] =
    g_signal_new("done-hamiltonian", G_TYPE_FROM_CLASS(klass),
                 G_SIGNAL_RUN_LAST | G_SIGNAL_NO_RECURSE | G_SIGNAL_NO_HOOKS,
        	 0, NULL, NULL, g_cclosure_marshal_VOID__OBJECT,
                 G_TYPE_NONE, 1, G_TYPE_OBJECT, NULL);

  bigdft_optloop_signals[DONE_SUBSPACE] =
    g_signal_new("done-subspace", G_TYPE_FROM_CLASS(klass),
                 G_SIGNAL_RUN_LAST | G_SIGNAL_NO_RECURSE | G_SIGNAL_NO_HOOKS,
        	 0, NULL, NULL, g_cclosure_marshal_VOID__OBJECT,
                 G_TYPE_NONE, 1, G_TYPE_OBJECT, NULL);

  bigdft_optloop_signals[DONE_WAVEFUNCTIONS] =
    g_signal_new("done-wavefunctions", G_TYPE_FROM_CLASS(klass),
                 G_SIGNAL_RUN_LAST | G_SIGNAL_NO_RECURSE | G_SIGNAL_NO_HOOKS,
        	 0, NULL, NULL, g_cclosure_marshal_VOID__OBJECT,
                 G_TYPE_NONE, 1, G_TYPE_OBJECT, NULL);
}
#endif

static void bigdft_optloop_init(BigDFT_OptLoop *obj)
{
#ifdef HAVE_GLIB
  memset((void*)((char*)obj + sizeof(GObject)), 0, sizeof(BigDFT_OptLoop) - sizeof(GObject));
#else
  memset(obj, 0, sizeof(BigDFT_OptLoop));
#endif
  obj->iscf = 0;
  obj->itrpmax = 1;
  obj->nrepmax = 1;
  obj->itermax = 50;
  obj->gnrm_cv = 1e-4;
  obj->rpnrm_cv = 1e-4;
  obj->gnrm_startmix = 0.;
}
static void bigdft_optloop_dispose(GObject *obj)
{
#ifdef HAVE_GLIB
  BigDFT_OptLoop *optloop = BIGDFT_OPTLOOP(obj);

  if (optloop->dispose_has_run)
    return;
  optloop->dispose_has_run = TRUE;

  /* Chain up to the parent class */
  G_OBJECT_CLASS(bigdft_optloop_parent_class)->dispose(obj);
#endif
}
static void bigdft_optloop_finalize(GObject *obj)
{
  BigDFT_OptLoop *optloop = BIGDFT_OPTLOOP(obj);

  if (optloop->data)
    FC_FUNC_(optloop_free, OPTLOOP_FREE)(&optloop->data);

#ifdef HAVE_GLIB
  G_OBJECT_CLASS(bigdft_optloop_parent_class)->finalize(obj);
#endif
}
BigDFT_OptLoop* bigdft_optloop_new()
{
  BigDFT_OptLoop *optloop;
  long self;

#ifdef HAVE_GLIB
  optloop = BIGDFT_OPTLOOP(g_object_new(BIGDFT_OPTLOOP_TYPE, NULL));
#else
  optloop = g_malloc(sizeof(BigDFT_OptLoop));
  bigdft_optloop_init(optloop);
#endif
  self = *((long*)&optloop);
  FC_FUNC_(optloop_new, OPTLOOP_NEW)(&self, &optloop->data);
  FC_FUNC_(optloop_sync_data, OPTLOOP_SYNC_DATA)
    (optloop->data, &optloop->gnrm_cv, &optloop->rpnrm_cv, &optloop->gnrm_startmix,
     &optloop->gnrm, &optloop->rpnrm,
     (int*)&optloop->itrpmax, (int*)&optloop->nrepmax, (int*)&optloop->itermax,
     (int*)&optloop->itrp, (int*)&optloop->itrep, (int*)&optloop->iter,
     &optloop->iscf, &optloop->infocode);

  return optloop;
}
void FC_FUNC_(optloop_new_wrapper, OPTLOOP_NEW_WRAPPER)(double *self, void *obj)
{
  BigDFT_OptLoop *optloop;

  optloop = bigdft_optloop_new_from_fortran(obj);
  *self = *((double*)&optloop);
}
BigDFT_OptLoop* bigdft_optloop_new_from_fortran(void *obj)
{
  BigDFT_OptLoop *optloop;

#ifdef HAVE_GLIB
  optloop = BIGDFT_OPTLOOP(g_object_new(BIGDFT_OPTLOOP_TYPE, NULL));
#else
  optloop = g_malloc(sizeof(BigDFT_OptLoop));
  bigdft_optloop_init(optloop);
#endif
  optloop->data = obj;

  return optloop;
}
void FC_FUNC_(optloop_free_wrapper, OPTLOOP_FREE_WRAPPER)(gpointer *obj)
{
  BigDFT_OptLoop *optloop = BIGDFT_OPTLOOP(*obj);

  optloop->data = (gpointer)0;
  bigdft_optloop_free(optloop);
}
void bigdft_optloop_free(BigDFT_OptLoop *optloop)
{
#ifdef HAVE_GLIB
  g_object_unref(G_OBJECT(optloop));
#else
  bigdft_optloop_finalize(optloop);
  g_free(optloop);
#endif
}
void FC_FUNC_(optloop_emit, OPTLOOP_EMIT)(BigDFT_OptLoop **obj, BigDFT_OptLoopIds *kind,
                                          BigDFT_Energs **energs)
{
  BigDFT_OptLoop *optloop = BIGDFT_OPTLOOP(*obj);

  bigdft_optloop_copy_from_fortran(optloop);
  bigdft_optloop_emit(optloop, *kind, *energs);
}
void bigdft_optloop_emit(BigDFT_OptLoop *optloop, BigDFT_OptLoopIds kind, BigDFT_Energs *energs)
{
#ifdef HAVE_GLIB
  switch (kind)
    {
    case BIGDFT_OPTLOOP_ITER_HAMILTONIAN:
      g_signal_emit(G_OBJECT(optloop), bigdft_optloop_signals[ITER_HAMILTONIAN],
                    0 /* details */, energs, NULL);
      break;
    case BIGDFT_OPTLOOP_ITER_SUBSPACE:
      g_signal_emit(G_OBJECT(optloop), bigdft_optloop_signals[ITER_SUBSPACE],
                    0 /* details */, energs, NULL);
      break;
    case BIGDFT_OPTLOOP_ITER_WAVEFUNCTIONS:
      g_signal_emit(G_OBJECT(optloop), bigdft_optloop_signals[ITER_WAVEFUNCTIONS],
                    0 /* details */, energs, NULL);
      break;
    case BIGDFT_OPTLOOP_DONE_HAMILTONIAN:
      g_signal_emit(G_OBJECT(optloop), bigdft_optloop_signals[DONE_HAMILTONIAN],
                    0 /* details */, energs, NULL);
      break;
    case BIGDFT_OPTLOOP_DONE_SUBSPACE:
      g_signal_emit(G_OBJECT(optloop), bigdft_optloop_signals[DONE_SUBSPACE],
                    0 /* details */, energs, NULL);
      break;
    case BIGDFT_OPTLOOP_DONE_WAVEFUNCTIONS:
      g_signal_emit(G_OBJECT(optloop), bigdft_optloop_signals[DONE_WAVEFUNCTIONS],
                    0 /* details */, energs, NULL);
      break;
    default:
      break;
    }
#endif
}
void bigdft_optloop_copy_from_fortran(BigDFT_OptLoop *optloop)
{
  FC_FUNC_(optloop_copy_data, OPTLOOP_COPY_DATA)
    (optloop->data, &optloop->gnrm_cv, &optloop->rpnrm_cv, &optloop->gnrm_startmix,
     &optloop->gnrm, &optloop->rpnrm,
     (int*)&optloop->itrpmax, (int*)&optloop->nrepmax, (int*)&optloop->itermax,
     (int*)&optloop->itrp, (int*)&optloop->itrep, (int*)&optloop->iter,
     &optloop->iscf, &optloop->infocode);
}
void bigdft_optloop_sync_to_fortran(BigDFT_OptLoop *optloop)
{
  FC_FUNC_(optloop_sync_data, OPTLOOP_SYNC_DATA)
    (optloop->data, &optloop->gnrm_cv, &optloop->rpnrm_cv, &optloop->gnrm_startmix,
     &optloop->gnrm, &optloop->rpnrm,
     (int*)&optloop->itrpmax, (int*)&optloop->nrepmax, (int*)&optloop->itermax,
     (int*)&optloop->itrp, (int*)&optloop->itrep, (int*)&optloop->iter,
     &optloop->iscf, &optloop->infocode);

#ifdef HAVE_GLIB
  g_signal_emit(G_OBJECT(optloop), bigdft_optloop_signals[SYNC_FORTRAN],
                0 /* details */, NULL);
#endif
}


/******************/
/* Miscellaneous. */
/******************/
double bigdft_memory_get_peak(guint nproc, const BigDFT_Locreg *lr, const BigDFT_Inputs *in,
                              const BigDFT_Orbs *orbs, const BigDFT_Proj *proj)
{
  double peak;

  FC_FUNC(memoryestimator, MEMORYESTIMATOR)((int*)&nproc, (int*)&in->idsx, lr->data,
                                            (int*)&lr->parent.nat, (int*)&orbs->norb,
                                            (int*)&orbs->nspinor, (int*)&orbs->nkpts,
                                            (int*)&proj->nprojel, (int*)&orbs->nspin,
                                            (int*)&in->itrpmax, (int*)&in->iscf, &peak);
  return peak;
}

/**
 * bigdft_init:
 * @mpi_iproc: (out):
 * @mpi_nproc: (out):
 * @mpi_igroup: (out):
 * @mpi_ngroup: (out):
 * @mpi_groupsize: 0 to use all MPI resources.
 *
 * Setup MPI and other variables.
 *
 * Returns: 
 **/
int bigdft_init(guint *mpi_iproc, guint *mpi_nproc, guint *mpi_igroup, guint *mpi_ngroup,
                guint mpi_groupsize)
{
  int ierr;
  int info[4];

  FC_FUNC_(bigdft_mpi_init, BIGDFT_MPI_INIT)(&ierr);
  FC_FUNC_(bigdft_init_mpi_env, BIGDFT_INIT_MPI_ENV)(info, &mpi_groupsize, &ierr);
  if (mpi_iproc)
    *mpi_iproc = (guint)info[0];
  if (mpi_nproc)
    *mpi_nproc = (guint)info[1];
  if (mpi_igroup)
    *mpi_igroup = (guint)info[2];
  if (mpi_ngroup)
    *mpi_ngroup = (guint)info[3];

  return ierr;
}
int bigdft_finalize()
{
  int ierr;

  FC_FUNC_(bigdft_finalize, BIGDFT_FINALIZE)(&ierr);
  return ierr;
}
/**
 * bigdft_set_input:
 * @radical: 
 * @posinp: 
 * @atoms: (out) (transfer full):
 *
 * Pouet.
 *
 * Returns: (transfer full):
 **/
BigDFT_Inputs* bigdft_set_input(const gchar *radical, const gchar *posinp, BigDFT_Atoms **atoms)
{
  BigDFT_Atoms *at;
  BigDFT_Inputs *in;

  at = bigdft_atoms_new();
  in = bigdft_inputs_init();
  FC_FUNC_(inputs_new, INPUTS_NEW)(&in->data);
  FC_FUNC_(bigdft_set_input, BIGDFT_SET_INPUT)(radical, posinp, &at->rxyz,
                                               in->data, at->data, strlen(radical), strlen(posinp));

  bigdft_atoms_copy_from_fortran(at);
  *atoms = at;
  return in;
}

/**
 * bigdft_eval_forces:
 * @atoms: 
 * @in: 
 * @rst:
 * @iproc: 
 * @nproc: 
 *
 * Pouet again.
 *
 * Returns: (transfer full):
 **/
BigDFT_Energs* bigdft_eval_forces(BigDFT_Atoms *atoms, BigDFT_Inputs *in, BigDFT_Restart *rst,
                                  guint iproc, guint nproc)
{
  int infocode;
  BigDFT_Energs *en;
  

  en = bigdft_energs_new();
  en->nat = atoms->nat;
  en->fxyz = g_malloc(sizeof(double) * atoms->nat * 3); 
  FC_FUNC_(call_bigdft, CALL_BIGDFT)(&nproc, &iproc, atoms->data, atoms->rxyz.data, in->data,
                                     &en->etot, en->fxyz, en->strten, &en->fnoise, rst->data,
                                     &infocode);
  FC_FUNC_(energs_copy_data, ENERGS_COPY_DATA)
    (en->data, &en->eh, &en->exc, &en->evxc, &en->eion, &en->edisp,
     &en->ekin, &en->epot, &en->eproj, &en->eexctX, &en->ebs, &en->eKS,
     &en->trH, &en->evsum, &en->evsic);
  
  return en;
}
