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

static void bigdft_locreg_dispose(GObject *atoms);
static void bigdft_locreg_finalize(GObject *atoms);

#ifdef HAVE_GLIB
G_DEFINE_TYPE(BigDFT_LocReg, bigdft_locreg, G_TYPE_OBJECT)

static void bigdft_locreg_class_init(BigDFT_LocRegClass *klass)
{
  /* Connect the overloading methods. */
  G_OBJECT_CLASS(klass)->dispose      = bigdft_locreg_dispose;
  G_OBJECT_CLASS(klass)->finalize     = bigdft_locreg_finalize;
  /* G_OBJECT_CLASS(klass)->set_property = visu_data_set_property; */
  /* G_OBJECT_CLASS(klass)->get_property = visu_data_get_property; */
}
#endif

static void bigdft_locreg_init(BigDFT_LocReg *glr)
{
#ifdef HAVE_GLIB
  memset(glr + sizeof(GObject), 0, sizeof(BigDFT_LocReg) - sizeof(GObject));
#else
  memset(glr, 0, sizeof(BigDFT_LocReg));
#endif
}
static void bigdft_locreg_dispose(GObject *obj)
{
#ifdef HAVE_GLIB
  BigDFT_LocReg *glr = BIGDFT_LOCREG(obj);

  if (glr->dispose_has_run)
    return;
  glr->dispose_has_run = TRUE;
  
  g_object_unref(G_OBJECT(glr->atoms));
  if (glr->data)
    FC_FUNC_(glr_empty, GLR_EMPTY)(glr->data);

  /* Chain up to the parent class */
  G_OBJECT_CLASS(bigdft_locreg_parent_class)->dispose(obj);
#endif
}
static void bigdft_locreg_finalize(GObject *obj)
{
  BigDFT_LocReg *glr = BIGDFT_LOCREG(obj);
  guint i;

  if (glr->data)
    FC_FUNC_(glr_free, GLR_FREE)(&glr->data);

#ifdef HAVE_GLIB
  G_OBJECT_CLASS(bigdft_locreg_parent_class)->finalize(obj);
#endif
}

BigDFT_LocReg* bigdft_locreg_new(BigDFT_Atoms *atoms, double *radii, double h[3],
                                 double crmult, double frmult)
{
  BigDFT_LocReg *glr;

#ifdef HAVE_GLIB
  glr = BIGDFT_LOCREG(g_object_new(BIGDFT_LOCREG_TYPE, NULL));
#else
  glr = g_malloc(sizeof(BigDFT_LocReg));
  bigdft_locreg_init(glr);
#endif

  FC_FUNC_(glr_new, GLR_NEW)(&glr->data);
  FC_FUNC_(glr_init, GLR_INIT)(glr->data, &glr->d);
  bigdft_locreg_define(glr, atoms, radii, h, crmult, frmult);

  return glr;
}
BigDFT_LocReg* bigdft_locreg_new_with_wave_descriptors(BigDFT_Atoms *atoms, double *radii,
                                                       double h[3], double crmult, double frmult)
{
  BigDFT_LocReg *glr;

  glr = bigdft_locreg_new(atoms, radii, h, crmult, frmult);
  bigdft_locreg_set_wave_descriptors(glr);
  
  return glr;
}
void bigdft_locreg_define(BigDFT_LocReg *glr, BigDFT_Atoms *atoms, double *radii,
                          double h[3], double crmult, double frmult)
{
  int iproc = 1;

  /* Deallocate all previously allocated data. */
  FC_FUNC_(glr_empty, GLR_EMPTY)(glr->data);

  /* Set the new size. */
  FC_FUNC_(system_size, SYSTEM_SIZE)(&iproc, atoms->data, atoms->rxyz.data, radii,
                                     &crmult, &frmult, h, h + 1, h + 2, glr->data,
                                     atoms->shift);
  /* Assign values. */
  glr->h[0] = h[0];
  glr->h[1] = h[1];
  glr->h[2] = h[2];
  FC_FUNC_(glr_get_dimensions, GLR_GET_DIMENSIONS)(glr->data, &glr->geocode,
                                                   (int*)glr->n, (int*)glr->ni);
#ifdef HAVE_GLIB
  if (glr->atoms)
    g_object_unref(G_OBJECT(glr->atoms));
  g_object_ref(G_OBJECT(atoms));
#endif
  glr->atoms  = atoms;
  glr->crmult = crmult;
  glr->frmult = frmult;
  if (glr->radii)
    g_free(glr->radii);
  glr->radii = g_malloc(sizeof(double) * atoms->ntypes * 3);
  memcpy(glr->radii, radii, sizeof(double) * atoms->ntypes * 3);
  
  /* Update Atoms accordingly. */
  FC_FUNC_(atoms_copy_alat, ATOMS_COPY_ALAT)(atoms->data, atoms->alat,
                                             atoms->alat + 1, atoms->alat + 2);
}
void bigdft_locreg_free(BigDFT_LocReg *glr)
{
#ifdef HAVE_GLIB
  g_object_unref(G_OBJECT(glr));
#else
  bigdft_locreg_finalize(glr);
  g_free(glr);
#endif
}
void bigdft_locreg_set_wave_descriptors(BigDFT_LocReg *glr)
{
  int iproc = 1;

  FC_FUNC_(glr_set_wave_descriptors,
           GLR_SET_WAVE_DESCRIPTORS)(&iproc, glr->h, glr->h + 1, glr->h + 2,
                                     glr->atoms->data, glr->atoms->rxyz.data, glr->radii,
                                     &glr->crmult, &glr->frmult, glr->data);
}


/********************************/
/* Bind the Lzd data structure. */
/********************************/
static void bigdft_lzd_dispose(GObject *atoms);
static void bigdft_lzd_finalize(GObject *atoms);

#ifdef HAVE_GLIB
G_DEFINE_TYPE(BigDFT_Lzd, bigdft_lzd, BIGDFT_LOCREG_TYPE)

static void bigdft_lzd_class_init(BigDFT_LzdClass *klass)
{
  /* Connect the overloading methods. */
  G_OBJECT_CLASS(klass)->dispose      = bigdft_lzd_dispose;
  G_OBJECT_CLASS(klass)->finalize     = bigdft_lzd_finalize;
  /* G_OBJECT_CLASS(klass)->set_property = visu_data_set_property; */
  /* G_OBJECT_CLASS(klass)->get_property = visu_data_get_property; */
}
#endif

static void bigdft_lzd_init(BigDFT_Lzd *lzd)
{
#ifdef HAVE_GLIB
  memset(lzd + sizeof(GObject), 0, sizeof(BigDFT_Lzd) - sizeof(GObject));
#else
  memset(lzd, 0, sizeof(BigDFT_Lzd));
#endif
  
  FC_FUNC_(lzd_new, LZD_NEW)(&lzd->data, &lzd->parent.data);
  FC_FUNC_(glr_init, GLR_INIT)(lzd->parent.data, &lzd->parent.d);
}
static void bigdft_lzd_dispose(GObject *obj)
{
#ifdef HAVE_GLIB
  BigDFT_Lzd *lzd = BIGDFT_LZD(obj);

  if (lzd->dispose_has_run)
    return;
  lzd->dispose_has_run = TRUE;
  
  lzd->parent.data = (void*)0;

  /* Chain up to the parent class */
  G_OBJECT_CLASS(bigdft_lzd_parent_class)->dispose(obj);
#endif
}
static void bigdft_lzd_finalize(GObject *obj)
{
  BigDFT_Lzd *lzd = BIGDFT_LZD(obj);

  FC_FUNC_(lzd_free, LZD_FREE)(&lzd->data);

#ifdef HAVE_GLIB
  G_OBJECT_CLASS(bigdft_lzd_parent_class)->finalize(obj);
#endif
}

BigDFT_Lzd* bigdft_lzd_new(BigDFT_Atoms *atoms, double *radii, double h[3],
                           double crmult, double frmult)
{
  BigDFT_Lzd *lzd;
  int iproc = 1;

#ifdef HAVE_GLIB
  lzd = BIGDFT_LZD(g_object_new(BIGDFT_LZD_TYPE, NULL));
#else
  lzd = g_malloc(sizeof(BigDFT_Lzd));
  bigdft_lzd_init(lzd);
#endif

  bigdft_locreg_define(BIGDFT_LOCREG(lzd), atoms, radii, h, crmult, frmult);

  return lzd;
}
void bigdft_lzd_free(BigDFT_Lzd *lzd)
{
#ifdef HAVE_GLIB
  g_object_unref(G_OBJECT(lzd));
#else
  bigdft_lzd_finalize(lzd);
  g_free(lzd);
#endif
}
