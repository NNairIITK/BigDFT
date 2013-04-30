#include <config.h>

#include "bigdft.h"
#include "bindings.h"
#include "bindings_api.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

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
