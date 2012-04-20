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


static void bigdft_orbs_dispose(GObject *orbs);
static void bigdft_orbs_finalize(GObject *orbs);

#ifdef HAVE_GLIB
G_DEFINE_TYPE(BigDFT_Orbs, bigdft_orbs, G_TYPE_OBJECT)

static void bigdft_orbs_class_init(BigDFT_OrbsClass *klass)
{
  /* Connect the overloading methods. */
  G_OBJECT_CLASS(klass)->dispose      = bigdft_orbs_dispose;
  G_OBJECT_CLASS(klass)->finalize     = bigdft_orbs_finalize;
  /* G_OBJECT_CLASS(klass)->set_property = visu_data_set_property; */
  /* G_OBJECT_CLASS(klass)->get_property = visu_data_get_property; */
}
#endif

static void bigdft_orbs_init(BigDFT_Orbs *obj)
{
#ifdef HAVE_GLIB
  memset((void*)((char*)obj + sizeof(GObject)), 0, sizeof(BigDFT_Orbs) - sizeof(GObject));
#else
  memset(obj, 0, sizeof(BigDFT_Orbs));
#endif
}
static void bigdft_orbs_dispose(GObject *obj)
{
#ifdef HAVE_GLIB
  BigDFT_Orbs *orbs = BIGDFT_ORBS(obj);

  if (orbs->dispose_has_run)
    return;
  orbs->dispose_has_run = TRUE;

  if (orbs->data)
    FC_FUNC_(orbs_empty, ORBS_EMPTY)(orbs->data);
  if (orbs->comm)
    FC_FUNC_(orbs_comm_empty, ORBS_COMM_EMPTY)(orbs->comm);
  
  /* Chain up to the parent class */
  G_OBJECT_CLASS(bigdft_orbs_parent_class)->dispose(obj);
#endif
}
static void bigdft_orbs_finalize(GObject *obj)
{
  BigDFT_Orbs *orbs = BIGDFT_ORBS(obj);

  if (orbs->data)
    FC_FUNC_(orbs_free, ORBS_FREE)(&orbs->data);
  if (orbs->comm)
    FC_FUNC_(orbs_comm_free, ORBS_COMM_FREE)(&orbs->comm);

#ifdef HAVE_GLIB
  G_OBJECT_CLASS(bigdft_orbs_parent_class)->finalize(obj);
#endif
}

BigDFT_Orbs* bigdft_orbs_new()
{
  BigDFT_Orbs *orbs;

#ifdef HAVE_GLIB
  orbs = BIGDFT_ORBS(g_object_new(BIGDFT_ORBS_TYPE, NULL));
#else
  orbs = g_malloc(sizeof(BigDFT_Orbs));
  bigdft_orbs_init(orbs);
#endif

  FC_FUNC_(orbs_new, ORBS_NEW)(&orbs->data);
  FC_FUNC_(orbs_init, ORBS_INIT)(orbs->data);
  
  return orbs;
}
void bigdft_orbs_free(BigDFT_Orbs *orbs)
{
#ifdef HAVE_GLIB
  g_object_unref(G_OBJECT(orbs));
#else
  bigdft_orbs_finalize(orbs);
  g_free(orbs);
#endif
}
guint bigdft_orbs_define(BigDFT_Orbs *orbs, const BigDFT_Lzd *lzd, const BigDFT_Inputs *in,
                         guint iproc, guint nproc)
{
  int nelec_, verb = 0;

  orbs->in = in;
  FC_FUNC_(orbs_empty, ORBS_EMPTY)(orbs->data);
  FC_FUNC_(read_orbital_variables, READ_ORBITAL_VARIABLES)(&iproc, &nproc, &verb, in->data,
                                                           lzd->parent.parent.data,
                                                           orbs->data, &nelec_);
  if (!orbs->comm)
    FC_FUNC_(orbs_comm_new, ORBS_COMM_NEW)(&orbs->comm);
  else
    FC_FUNC_(orbs_comm_empty, ORBS_COMM_EMPTY)(orbs->comm);
  FC_FUNC_(orbs_comm_init, ORBS_COMM_INIT)(orbs->comm, orbs->data,
                                           BIGDFT_LOCREG(lzd)->data, &iproc, &nproc);
  
  FC_FUNC_(orbs_get_dimensions, ORBS_GET_DIMENSIONS)(orbs->data, &orbs->norb,
                                                     &orbs->norbp, &orbs->norbu,
                                                     &orbs->norbd, &orbs->nspin,
                                                     &orbs->nspinor, &orbs->npsidim,
                                                     &orbs->nkpts, &orbs->nkptsp,
                                                     &orbs->isorb, &orbs->iskpts);
  GET_ATTR_DBL   (orbs, ORBS, occup, OCCUP);
  GET_ATTR_DBL   (orbs, ORBS, kwgts, KWGTS);
  GET_ATTR_DBL_2D(orbs, ORBS, kpts,  KPTS);

  return nelec_;
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

/******************************************/
/* Devel. version of a wavefunction type. */
/******************************************/
static void bigdft_wf_dispose(GObject *atoms);
static void bigdft_wf_finalize(GObject *atoms);

#ifdef HAVE_GLIB
enum {
  PSI_READY_SIGNAL,
  ONE_WAVE_READY_SIGNAL,
  LAST_SIGNAL
};

G_DEFINE_TYPE(BigDFT_Wf, bigdft_wf, BIGDFT_ORBS_TYPE)

static guint bigdft_wf_signals[LAST_SIGNAL] = { 0 };

static void g_cclosure_marshal_ONE_WAVE(GClosure *closure,
                                        GValue *return_value,
                                        guint n_param_values,
                                        const GValue *param_values,
                                        gpointer invocation_hint,
                                        gpointer marshal_data)
{
  typedef void (*callbackFunc)(gpointer data1, guint iter, GArray *arg_psi, guint arg_kpt,
                               guint arg_orb, guint arg_spin, gpointer data2);
  register callbackFunc callback;
  register GCClosure *cc = (GCClosure*)closure;
  register gpointer data1, data2;

  g_return_if_fail(n_param_values == 6);

  if (G_CCLOSURE_SWAP_DATA(closure))
    {
      data1 = closure->data;
      data2 = g_value_peek_pointer(param_values + 0);
    }
  else
    {
      data1 = g_value_peek_pointer(param_values + 0);
      data2 = closure->data;
    }
  callback = (callbackFunc)(size_t)(marshal_data ? marshal_data : cc->callback);

  callback(data1, g_value_get_uint(param_values + 1),
           (GArray*)g_value_get_boxed(param_values + 2),
           g_value_get_uint(param_values + 3), 
           g_value_get_uint(param_values + 4), 
           g_value_get_uint(param_values + 5), data2);
}

static void bigdft_wf_class_init(BigDFT_WfClass *klass)
{
  /* Connect the overloading methods. */
  G_OBJECT_CLASS(klass)->dispose      = bigdft_wf_dispose;
  G_OBJECT_CLASS(klass)->finalize     = bigdft_wf_finalize;
  /* G_OBJECT_CLASS(klass)->set_property = visu_data_set_property; */
  /* G_OBJECT_CLASS(klass)->get_property = visu_data_get_property; */

  bigdft_wf_signals[PSI_READY_SIGNAL] =
    g_signal_new("psi-ready", G_TYPE_FROM_CLASS(klass),
                 G_SIGNAL_RUN_LAST | G_SIGNAL_NO_RECURSE | G_SIGNAL_NO_HOOKS,
		 0, NULL, NULL, g_cclosure_marshal_VOID__UINT,
                 G_TYPE_NONE, 1, G_TYPE_UINT, NULL);

  bigdft_wf_signals[ONE_WAVE_READY_SIGNAL] =
    g_signal_new("one-wave-ready", G_TYPE_FROM_CLASS(klass),
                 G_SIGNAL_RUN_LAST | G_SIGNAL_NO_RECURSE | G_SIGNAL_NO_HOOKS | G_SIGNAL_DETAILED,
		 0, NULL, NULL, g_cclosure_marshal_ONE_WAVE,
                 G_TYPE_NONE, 5, G_TYPE_UINT, G_TYPE_ARRAY, G_TYPE_UINT,
                 G_TYPE_UINT, G_TYPE_UINT, NULL);
}
#endif

static void bigdft_wf_init(BigDFT_Wf *obj)
{
#ifdef HAVE_GLIB
  memset((void*)((char*)obj + sizeof(BigDFT_Orbs)), 0, sizeof(BigDFT_Wf) - sizeof(BigDFT_Orbs));
#else
  memset(obj, 0, sizeof(BigDFT_Wf));
#endif
}
static void bigdft_wf_dispose(GObject *obj)
{
#ifdef HAVE_GLIB
  BigDFT_Wf *wf = BIGDFT_WF(obj);

  if (wf->dispose_has_run)
    return;
  wf->dispose_has_run = TRUE;

  wf->parent.data = (void*)0;
  wf->parent.comm = (void*)0;

  /* Destroy only the C wrappers. */
  wf->lzd->data = (void*)0;
  g_object_unref(G_OBJECT(wf->lzd));

  if (wf->data)
    FC_FUNC_(wf_empty, WF_EMPTY)(wf->data);

  /* Chain up to the parent class */
  G_OBJECT_CLASS(bigdft_wf_parent_class)->dispose(obj);
#endif
}
static void bigdft_wf_finalize(GObject *obj)
{
  BigDFT_Wf *wf = BIGDFT_WF(obj);

  if (wf->data)
    FC_FUNC_(wf_free, WF_FREE)(&wf->data);

#ifdef HAVE_GLIB
  G_OBJECT_CLASS(bigdft_wf_parent_class)->finalize(obj);
  /* g_debug("Freeing wf object %p done.\n", obj); */
#endif
}
void FC_FUNC_(wf_emit_psi, WF_EMIT_PSI)(BigDFT_Wf **wf, guint *istep)
{
#ifdef HAVE_GLIB
  g_signal_emit(G_OBJECT(*wf), bigdft_wf_signals[PSI_READY_SIGNAL],
                0 /* details */, *istep, NULL);
#endif  
}
void bigdft_wf_emit_one_wave(BigDFT_Wf *wf, guint iter, GArray *psic,
                             GQuark quark, guint ikpt, guint iorb, guint ispin)
{
#ifdef HAVE_GLIB
  g_signal_emit(G_OBJECT(wf), bigdft_wf_signals[ONE_WAVE_READY_SIGNAL],
                quark, iter, psic, ikpt, iorb, ispin, NULL);
#endif  
}

BigDFT_Wf* bigdft_wf_new()
{
  double self;
  BigDFT_Wf *wf;

#ifdef HAVE_GLIB
  wf = BIGDFT_WF(g_object_new(BIGDFT_WF_TYPE, NULL));
#else
  wf = g_malloc(sizeof(BigDFT_Wf));
  bigdft_wf_init(wf);
#endif
  self = *((double*)&wf);
  FC_FUNC_(wf_new, WF_NEW)(&self, &wf->data, &wf->parent.data, &wf->parent.comm,
                           &wf->data_lzd);
  FC_FUNC_(orbs_init, ORBS_INIT)(wf->parent.data);
  FC_FUNC_(wf_get_psi, WF_GET_PSI)(wf->data, &wf->psi);

  wf->lzd = bigdft_lzd_new_with_fortran(wf->data_lzd);

  return wf;
}
void FC_FUNC_(wf_new_wrapper, WF_NEW_WRAPPER)(double *self, void *obj)
{
  BigDFT_Wf *wf;

  wf = bigdft_wf_new_from_fortran(obj);
  *self = *((double*)&wf);
}
BigDFT_Wf* bigdft_wf_new_from_fortran(void *obj)
{
  BigDFT_Wf *wf;

#ifdef HAVE_GLIB
  wf = BIGDFT_WF(g_object_new(BIGDFT_WF_TYPE, NULL));
#else
  wf = g_malloc(sizeof(BigDFT_Wf));
  bigdft_wf_init(wf);
#endif
  wf->data = obj;
  FC_FUNC_(wf_get_data, WF_GET_DATA)(wf->data, &wf->parent.data, &wf->parent.comm,
                                     &wf->data_lzd);
  FC_FUNC_(wf_get_psi, WF_GET_PSI)(wf->data, &wf->psi);

  wf->lzd = bigdft_lzd_new_from_fortran(wf->data_lzd);

  return wf;
}
void FC_FUNC_(wf_free_wrapper, WF_FREE_WRAPPER)(gpointer *obj)
{
  BigDFT_Wf *wf = BIGDFT_WF(*obj);

  wf->data = (gpointer)0;
  bigdft_wf_free(wf);
}
void bigdft_wf_free(BigDFT_Wf *wf)
{
#ifdef HAVE_GLIB
  g_object_unref(G_OBJECT(wf));
#else
  bigdft_wf_finalize(wf);
  g_free(wf);
#endif
}
void FC_FUNC_(wf_copy_from_fortran, WF_COPY_FROM_FORTRAN)
     (gpointer *self, const double *radii, const double *crmult, const double *frmult)
{
  BigDFT_Wf *wf = BIGDFT_WF(*self);

  bigdft_lzd_copy_from_fortran(wf->lzd, radii, *crmult, *frmult);
}
guint bigdft_wf_define(BigDFT_Wf *wf, const BigDFT_Inputs *in, guint iproc, guint nproc)
{
  int nelec;

  nelec = bigdft_orbs_define(&wf->parent, wf->lzd, in, iproc, nproc);
  
  FC_FUNC_(lzd_empty, LZD_EMPTY)(wf->lzd->data);
  FC_FUNC_(check_linear_and_create_lzd, CHECK_LINEAR_AND_CREATE_LZD)
    (&iproc, &nproc, in->data, wf->lzd->data, wf->lzd->parent.parent.data,
     wf->parent.data, wf->lzd->parent.parent.rxyz.data);

  FC_FUNC_(wf_empty, WF_EMPTY)(wf->data);

  return nelec;
}
void bigdft_wf_calculate_psi0(BigDFT_Wf *wf, BigDFT_LocalFields *denspot,
                              BigDFT_Proj *proj, BigDFT_Energs *energs,
                              guint iproc, guint nproc)
{
  int inputpsi;
  guint norbv;
  void *GPU;
  BigDFT_Orbs *orbs;

  FC_FUNC_(gpu_new, GPU_NEW)(&GPU);
  FC_FUNC_(input_wf, INPUT_WF)(&iproc, &nproc, wf->parent.in->data, GPU,
                               BIGDFT_ATOMS(wf->lzd)->data,
                               BIGDFT_ATOMS(wf->lzd)->rxyz.data,
                               denspot->data, proj->nlpspd, &proj->proj,
                               wf->data, energs->data, &inputpsi, &norbv,
                               (void*)0, (void*)0, (void*)0, (void*)0, (void*)0,
                               (void*)0, (void*)0);
  FC_FUNC_(gpu_free, GPU_FREE)(&GPU);
  orbs = &wf->parent;
  GET_ATTR_DBL(orbs, ORBS, eval,  EVAL);
}
guint bigdft_wf_optimization_loop(BigDFT_Wf *wf, BigDFT_LocalFields *denspot,
                                  BigDFT_Proj *proj, BigDFT_Energs *energs,
                                  guint iproc, guint nproc, BigDFT_optLoopParams *params)
{
  guint infocode, itrp, icycle, iter;
  guint inputpsi = 0;
  double xcstr[6];
  void *GPU;
  BigDFT_optLoopParams p;

  if (params)
    p = *params;
  else
    bigdft_optloopparams_init(&p);
  FC_FUNC_(gpu_new, GPU_NEW)(&GPU);
  FC_FUNC_(kswfn_optimization_loop, KSWFN_OPTIMIZATION_LOOP)
    (&infocode, &itrp, &icycle, &iter, &iproc, &nproc,
     &p.iscf, &p.itrpmax, &p.nrepmax, &p.itermax, &p.gnrm_cv, &p.rpnrm_cv,
     &p.gnrm_startmix, &p.alphamix, &p.idsx, &inputpsi,
     wf->data, denspot->data, proj->nlpspd, &proj->proj,
     energs->data, BIGDFT_ATOMS(wf->lzd)->data, BIGDFT_ATOMS(wf->lzd)->rxyz.data,
     GPU, xcstr, wf->parent.in->data);
  FC_FUNC_(energs_copy_data, ENERGS_COPY_DATA)
    (energs->data, &energs->eh, &energs->exc,
     &energs->evxc, &energs->eion, &energs->edisp,
     &energs->ekin, &energs->epot, &energs->eproj,
     &energs->eexctX, &energs->ebs, &energs->eKS,
     &energs->trH, &energs->evsum, &energs->evsic);
  FC_FUNC_(gpu_free, GPU_FREE)(&GPU);

  return infocode;
}
const double* bigdft_wf_get_psi_compress(const BigDFT_Wf *wf, guint ikpt, guint iorb,
                                         BigDFT_Spin ispin, BigDFT_Spinor ispinor,
                                         guint *psiSize, guint iproc)
{
  guint ispinor_, orbSize;
  int iorbp, jproc;

  *psiSize = 0;
  if (ispin == BIGDFT_SPIN_DOWN && wf->parent.norbd == 0)
    return (double*)0;
  if (ispinor == BIGDFT_IMAG && wf->parent.nspinor == 1)
    return (double*)0;

  /* Get the shift to apply on wf->psi to get the right orbital. */
  ispinor_  = (ispinor != BIGDFT_PARTIAL_DENSITY)?ispinor:BIGDFT_REAL;
  ispinor_ += 1;
  ispin    += 1;
  FC_FUNC_(orbs_get_iorbp, ORBS_GET_IORBP)(wf->parent.data, &iorbp, &jproc,
                                           &ikpt, &iorb, &ispin, &ispinor_);
  if (iorbp < 0)
    return (double*)0;
  FC_FUNC_(glr_get_psi_size, GLR_GET_PSI_SIZE)(wf->lzd->parent.data, &orbSize);
  *psiSize = orbSize;
  if (ispinor == BIGDFT_PARTIAL_DENSITY && wf->parent.nspinor == 2)
    *psiSize *= 2;
  
  return (iproc == jproc)?wf->psi->data + (iorbp * orbSize):(double*)0;
}
gboolean bigdft_wf_copy_psi_compress(const BigDFT_Wf *wf, guint ikpt, guint iorb,
                                     BigDFT_Spin ispin, BigDFT_Spinor ispinor,
                                     guint iproc, double *psic, guint psiSize)
{
  guint ispinor_, orbSize;
  int iorbp, jproc;

  if (ispin == BIGDFT_SPIN_DOWN && wf->parent.norbd == 0)
    return FALSE;
  if (ispinor == BIGDFT_IMAG && wf->parent.nspinor == 1)
    return FALSE;

  /* Get the shift to apply on wf->psi to get the right orbital. */
  ispinor_  = (ispinor != BIGDFT_PARTIAL_DENSITY)?ispinor:BIGDFT_REAL;
  ispinor_ += 1;
  ispin    += 1;
  FC_FUNC_(orbs_get_iorbp, ORBS_GET_IORBP)(wf->parent.data, &iorbp, &jproc,
                                           &ikpt, &iorb, &ispin, &ispinor_);
  if (iorbp < 0)
    return FALSE;
  FC_FUNC_(glr_get_psi_size, GLR_GET_PSI_SIZE)(wf->lzd->parent.data, &orbSize);
  if (ispinor == BIGDFT_PARTIAL_DENSITY && wf->parent.nspinor == 2)
    {
      if (psiSize != 2 * orbSize)
        return FALSE;
    }
  else
    {
      if (psiSize != orbSize)
        return FALSE;
    }
  if (iproc == jproc)
    memcpy(psic, wf->psi->data + (iorbp * orbSize), sizeof(double) * psiSize);
  else
    FC_FUNC_(kswfn_mpi_copy, KSWFN_MPI_COPY)(psic, &jproc, &iorbp, &psiSize);
  
  return TRUE;
}
double* bigdft_wf_convert_to_isf(const BigDFT_Wf *wf, guint ikpt, guint iorb,
                                 BigDFT_Spin ispin, BigDFT_Spinor ispinor, guint iproc)
{
  guint psiSize, i, n;
  const double *psic;
  double *psir, *psii;

  psic = bigdft_wf_get_psi_compress(wf, ikpt, iorb, ispin, ispinor, &psiSize, iproc);
  if (!psic)
    return (double *)0;
  
  psir = bigdft_locreg_convert_to_isf(&wf->lzd->parent, psic);
  if (ispinor == BIGDFT_PARTIAL_DENSITY)
    {
      n = wf->lzd->parent.ni[0] * wf->lzd->parent.ni[1] * wf->lzd->parent.ni[2];
      if (wf->parent.nspinor == 2)
        {
          psii = bigdft_locreg_convert_to_isf(&wf->lzd->parent, psic + psiSize / 2);
          for(i = 0; i < n; i++)
            psir[i] = psir[i] * psir[i] + psii[i] * psii[i];
          g_free(psii);
        }
      else
        for(i = 0; i < n; i++)
          psir[i] *= psir[i];;
    }

  return psir;
}
typedef struct bigdft_data
{
  guint                iproc, nproc;
  const BigDFT_Inputs *in;
  BigDFT_Proj         *proj;
  BigDFT_LocalFields  *denspot;
  BigDFT_Wf           *wf;
  BigDFT_Energs       *energs;
} BigDFT_Data;
static gpointer wf_optimization_thread(gpointer data)
{
  BigDFT_Data *ct = (BigDFT_Data*)data;
  
  bigdft_localfields_create_poisson_kernels(ct->denspot, ct->wf->lzd,
                                            ct->in, ct->iproc, ct->nproc);
  bigdft_localfields_create_effective_ionic_pot(ct->denspot, ct->wf->lzd,
                                                ct->in, ct->iproc, ct->nproc);
  bigdft_wf_calculate_psi0(ct->wf, ct->denspot, ct->proj, ct->energs, ct->iproc, ct->nproc);
  bigdft_wf_optimization_loop(ct->wf, ct->denspot, ct->proj, ct->energs,
                              ct->iproc, ct->nproc, (BigDFT_optLoopParams*)0);
#ifdef HAVE_GLIB
  g_object_unref(G_OBJECT(ct->wf));
  g_object_unref(G_OBJECT(ct->denspot));
  g_object_unref(G_OBJECT(ct->energs));
  g_object_unref(G_OBJECT(ct->proj));
#endif
  g_free(ct);

  return (gpointer)0;
}
void bigdft_wf_optimization(BigDFT_Wf *wf, BigDFT_Proj *proj, BigDFT_LocalFields *denspot,
                            BigDFT_Energs *energs, const BigDFT_Inputs *in,
                            gboolean threaded, guint iproc, guint nproc)
{
  BigDFT_Data *ct;
#ifdef HAVE_GLIB
  GThread *ld_thread;
  GError *error = (GError*)0;
#endif

  ct = g_malloc(sizeof(BigDFT_Data));
  ct->iproc   = iproc;
  ct->nproc   = nproc;
  ct->denspot = denspot;
  ct->in      = in;
  ct->proj    = proj;
  ct->wf      = wf;
  ct->energs  = energs;
#ifdef HAVE_GLIB
  g_object_ref(G_OBJECT(wf));
  g_object_ref(G_OBJECT(denspot));
  g_object_ref(G_OBJECT(proj));
  g_object_ref(G_OBJECT(energs));
#endif
#ifdef G_THREADS_ENABLED
  if (threaded)
    ld_thread = g_thread_create(wf_optimization_thread, ct, FALSE, &error);
  else
    wf_optimization_thread(ct);
#else
  wf_optimization_thread(ct);
#endif
}

void bigdft_optloopparams_init(BigDFT_optLoopParams *params)
{
  params->iscf = 0;
  params->itrpmax = 1;
  params->nrepmax = 1;
  params->itermax = 50;
  params->gnrm_cv = 1e-4;
  params->rpnrm_cv = 1e-4;
  params->gnrm_startmix = 0.;
  params->alphamix = 0.;
  params->idsx = 6;
}
