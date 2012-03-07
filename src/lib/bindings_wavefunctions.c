#include <config.h>

#ifdef HAVE_GLIB
#include <glib-object.h>
#endif

#include "bigdft.h"
#include "bindings.h"
#include "bindings_api.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>


static void bigdft_orbs_dispose(GObject *orbs);
static void bigdft_orbs_finalize(GObject *orbs);
static guint bigdft_orbs_define(BigDFT_Orbs *orbs,
                                const BigDFT_LocReg *glr, const BigDFT_Inputs *in,
                                guint iproc, guint nproc);

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
  FC_FUNC_(orbs_new, ORBS_NEW)(&obj->data);
}
static void bigdft_orbs_dispose(GObject *obj)
{
#ifdef HAVE_GLIB
  BigDFT_Orbs *orbs = BIGDFT_ORBS(obj);

  if (orbs->dispose_has_run)
    return;
  orbs->dispose_has_run = TRUE;

  g_object_unref(G_OBJECT(orbs->glr));

  /* Chain up to the parent class */
  G_OBJECT_CLASS(bigdft_orbs_parent_class)->dispose(obj);
#endif
}
static void bigdft_orbs_finalize(GObject *obj)
{
  BigDFT_Orbs *orbs = BIGDFT_ORBS(obj);

  FC_FUNC_(orbs_free, ORBS_FREE)(&orbs->data);
  FC_FUNC_(orbs_comm_free, ORBS_COMM_FREE)(&orbs->comm);

#ifdef HAVE_GLIB
  G_OBJECT_CLASS(bigdft_orbs_parent_class)->finalize(obj);
#endif
}

BigDFT_Orbs* bigdft_orbs_new(const BigDFT_LocReg *glr, const BigDFT_Inputs *in,
                             guint iproc, guint nproc, guint *nelec)
{
  BigDFT_Orbs *orbs;
  int nelec_, verb = 0;

#ifdef HAVE_GLIB
  orbs = BIGDFT_ORBS(g_object_new(BIGDFT_ORBS_TYPE, NULL));
#else
  orbs = g_malloc(sizeof(BigDFT_Orbs));
  bigdft_orbs_init(orbs);
#endif

  nelec_ = bigdft_orbs_define(orbs, glr, in, iproc, nproc);
  if (nelec)
    *nelec = nelec_;
  
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
static guint bigdft_orbs_define(BigDFT_Orbs *orbs,
                                const BigDFT_LocReg *glr, const BigDFT_Inputs *in,
                                guint iproc, guint nproc)
{
  int nelec_, verb = 0;

  orbs->in = in;
#ifdef HAVE_GLIB
  if (orbs->glr)
    g_object_unref(G_OBJECT(orbs->glr));
  g_object_ref(G_OBJECT(glr));
#endif
  orbs->glr = glr;
  FC_FUNC_(read_orbital_variables, READ_ORBITAL_VARIABLES)(&iproc, &nproc, &verb, in->data,
                                                           glr->parent.data,
                                                           orbs->data, &nelec_);
  FC_FUNC_(orbs_comm, ORBS_COMM)(&orbs->comm, orbs->data, glr->data, &iproc, &nproc);
  
  FC_FUNC_(orbs_get_dimensions, ORBS_GET_DIMENSIONS)(orbs->data, &orbs->norb,
                                                     &orbs->norbp, &orbs->norbu,
                                                     &orbs->norbd, &orbs->nspin,
                                                     &orbs->nspinor, &orbs->npsidim,
                                                     &orbs->nkpts, &orbs->nkptsp,
                                                     &orbs->isorb, &orbs->iskpts);
  
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
G_DEFINE_TYPE(BigDFT_Wf, bigdft_wf, BIGDFT_ORBS_TYPE)

static void bigdft_wf_class_init(BigDFT_WfClass *klass)
{
  /* Connect the overloading methods. */
  G_OBJECT_CLASS(klass)->dispose      = bigdft_wf_dispose;
  G_OBJECT_CLASS(klass)->finalize     = bigdft_wf_finalize;
  /* G_OBJECT_CLASS(klass)->set_property = visu_data_set_property; */
  /* G_OBJECT_CLASS(klass)->get_property = visu_data_get_property; */
}
#endif

static void bigdft_wf_init(BigDFT_Wf *obj)
{
#ifdef HAVE_GLIB
  memset((void*)((char*)obj + sizeof(BigDFT_Orbs)), 0, sizeof(BigDFT_Wf) - sizeof(BigDFT_Orbs));
#else
  memset(obj, 0, sizeof(BigDFT_Wf));
#endif
  
  F90_1D_POINTER_INIT(&obj->psi);
  F90_1D_POINTER_INIT(&obj->hpsi);
  F90_1D_POINTER_INIT(&obj->psit);
}
static void bigdft_wf_dispose(GObject *obj)
{
#ifdef HAVE_GLIB
  BigDFT_Wf *wf = BIGDFT_WF(obj);

  if (wf->dispose_has_run)
    return;
  wf->dispose_has_run = TRUE;

  g_object_unref(G_OBJECT(wf->lzd));

  /* Chain up to the parent class */
  G_OBJECT_CLASS(bigdft_wf_parent_class)->dispose(obj);
#endif
}
static void bigdft_wf_finalize(GObject *obj)
{
  BigDFT_Wf *wf = BIGDFT_WF(obj);

  FC_FUNC_(deallocate_double_1d, DEALLOCATE_DOUBLE_1D)(&wf->psi);
  FC_FUNC_(deallocate_double_1d, DEALLOCATE_DOUBLE_1D)(&wf->hpsi);
  FC_FUNC_(deallocate_double_1d, DEALLOCATE_DOUBLE_1D)(&wf->psit);

#ifdef HAVE_GLIB
  G_OBJECT_CLASS(bigdft_wf_parent_class)->finalize(obj);
#endif
}

BigDFT_Wf* bigdft_wf_new(BigDFT_Lzd *lzd, BigDFT_Inputs *in,
                         guint iproc, guint nproc, guint *nelec)
{
  BigDFT_Wf *wf;
  guint nelec_;

#ifdef HAVE_GLIB
  wf = BIGDFT_WF(g_object_new(BIGDFT_WF_TYPE, NULL));
#else
  wf = g_malloc(sizeof(BigDFT_Wf));
  bigdft_wf_init(wf);
#endif

  nelec_ = bigdft_orbs_define(&wf->parent, BIGDFT_LOCREG(lzd), in, iproc, nproc);
  if (nelec)
    *nelec = nelec_;
  
#ifdef HAVE_GLIB
  g_object_ref(G_OBJECT(lzd));
#endif
  wf->lzd = lzd;

  return wf;
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
void bigdft_wf_calculate_psi0(BigDFT_Wf *wf, BigDFT_LocalFields *denspot, BigDFT_Proj *proj,
                              guint iproc, guint nproc)
{
  int inputpsi, norbv;
  void *GPU;

  FC_FUNC_(gpu_new, GPU_NEW)(&GPU);
  FC_FUNC_(input_wf, INPUT_WF)(&iproc, &nproc, wf->parent.in->data, GPU,
                               BIGDFT_ATOMS(wf->lzd)->data,
                               BIGDFT_ATOMS(wf->lzd)->rxyz.data,
                               wf->lzd->data, BIGDFT_LOCREG(wf->lzd)->h,
                               BIGDFT_LOCREG(wf->lzd)->h + 1, BIGDFT_LOCREG(wf->lzd)->h + 2,
                               denspot->data, proj->nlpspd, &proj->proj,
                               wf->parent.data, wf->parent.comm,
                               &wf->psi, &wf->hpsi, &wf->psit,
                               &inputpsi, &norbv,
                               (void*)0, (void*)0, (void*)0, (void*)0, (void*)0,
                               (void*)0, (void*)0, (void*)0, (void*)0
                               );
  FC_FUNC_(gpu_free, GPU_FREE)(&GPU);
}
