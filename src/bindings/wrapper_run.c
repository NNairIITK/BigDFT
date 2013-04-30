#include <config.h>

#include "bigdft.h"
#include "bindings.h"
#include "bindings_api.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

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
#ifdef HAVE_GLIB
static void bigdft_energs_dispose(GObject *obj)
{
  BigDFT_Energs *energs = BIGDFT_ENERGS(obj);

  if (energs->dispose_has_run)
    return;
  energs->dispose_has_run = TRUE;

  /* Chain up to the parent class */
  G_OBJECT_CLASS(bigdft_energs_parent_class)->dispose(obj);
}
#endif
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
void bigdft_restart_set_mode(BigDFT_Restart *restart, BigDFT_RestartModes id)
{
  int inputPsiId[] = {0, 1};

  restart->inputPsiId = id;
  FC_FUNC_(inputs_set_restart, INPUTS_SET_RESTART)(restart->in->data, inputPsiId + id);
}

/**
 * bigdft_run_run:
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
BigDFT_Energs* bigdft_run_run(BigDFT_Atoms *atoms, BigDFT_Inputs *in, BigDFT_Restart *rst,
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
