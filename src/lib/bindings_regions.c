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


static void bigdft_locreg_dispose(GObject *atoms);
static void bigdft_locreg_finalize(GObject *atoms);

#ifdef HAVE_GLIB
G_DEFINE_TYPE(BigDFT_LocReg, bigdft_locreg, BIGDFT_ATOMS_TYPE)

static void bigdft_locreg_class_init(BigDFT_LocRegClass *klass)
{
  /* Connect the overloading methods. */
  G_OBJECT_CLASS(klass)->dispose      = bigdft_locreg_dispose;
  G_OBJECT_CLASS(klass)->finalize     = bigdft_locreg_finalize;
  /* G_OBJECT_CLASS(klass)->set_property = visu_data_set_property; */
  /* G_OBJECT_CLASS(klass)->get_property = visu_data_get_property; */
}
#endif

static void bigdft_locreg_init(BigDFT_LocReg *obj)
{
#ifdef HAVE_GLIB
  memset((void*)((char*)obj + sizeof(BigDFT_Atoms)), 0, sizeof(BigDFT_LocReg) - sizeof(BigDFT_Atoms));
#else
  memset(obj, 0, sizeof(BigDFT_LocReg));
#endif
}
static void bigdft_locreg_dispose(GObject *obj)
{
#ifdef HAVE_GLIB
  BigDFT_LocReg *glr = BIGDFT_LOCREG(obj);

  if (glr->dispose_has_run)
    return;
  glr->dispose_has_run = TRUE;

  if (glr->data)
    FC_FUNC_(glr_empty, GLR_EMPTY)(glr->data);

  /* Chain up to the parent class */
  G_OBJECT_CLASS(bigdft_locreg_parent_class)->dispose(obj);
#endif
}
static void bigdft_locreg_finalize(GObject *obj)
{
  BigDFT_LocReg *glr = BIGDFT_LOCREG(obj);

  if (glr->data)
    FC_FUNC_(glr_free, GLR_FREE)(&glr->data);

#ifdef HAVE_GLIB
  G_OBJECT_CLASS(bigdft_locreg_parent_class)->finalize(obj);
#endif
}

BigDFT_LocReg* bigdft_locreg_new()
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

  return glr;
}
void bigdft_locreg_set_radii(BigDFT_LocReg *glr, const double *radii)
{
  if (glr->radii)
    g_free(glr->radii);
  glr->radii = g_malloc(sizeof(double) * glr->parent.ntypes * 3);
  memcpy(glr->radii, radii, sizeof(double) * glr->parent.ntypes * 3);
}
void bigdft_locreg_set_size(BigDFT_LocReg *glr, const double h[3],
                            double crmult, double frmult)
{
  int iproc = 1;

  /* Deallocate all previously allocated data. */
  FC_FUNC_(glr_empty, GLR_EMPTY)(glr->data);

  /* Assign values. */
  glr->h[0] = h[0];
  glr->h[1] = h[1];
  glr->h[2] = h[2];
  /* Set the new size. */
  FC_FUNC_(system_size, SYSTEM_SIZE)(&iproc, glr->parent.data, glr->parent.rxyz.data,
                                     glr->radii, &crmult, &frmult, glr->h, glr->h + 1,
                                     glr->h + 2, glr->data, glr->parent.shift);
  FC_FUNC_(glr_get_dimensions, GLR_GET_DIMENSIONS)(glr->data, (int*)glr->n, (int*)glr->ni);
  glr->crmult = crmult;
  glr->frmult = frmult;

  /* Update Atoms accordingly. */
  FC_FUNC_(atoms_copy_alat, ATOMS_COPY_ALAT)(glr->parent.data, glr->parent.alat,
                                             glr->parent.alat + 1, glr->parent.alat + 2);
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
static void bigdft_locreg_copy_data(BigDFT_LocReg *glr,
                                    const double *radii, const double h[3],
                                    double crmult, double frmult)
{
  FC_FUNC_(glr_get_dimensions, GLR_GET_DIMENSIONS)(glr->data, (int*)glr->n, (int*)glr->ni);
  bigdft_locreg_set_radii(glr, radii);
  glr->h[0] = h[0];
  glr->h[1] = h[1];
  glr->h[2] = h[2];
  glr->crmult = crmult;
  glr->frmult = frmult;
}
void bigdft_locreg_set_wave_descriptors(BigDFT_LocReg *glr)
{
  int iproc = 1;

  FC_FUNC_(glr_set_wave_descriptors,
           GLR_SET_WAVE_DESCRIPTORS)(&iproc, glr->h, glr->h + 1, glr->h + 2,
                                     glr->parent.data, glr->parent.rxyz.data, glr->radii,
                                     &glr->crmult, &glr->frmult, glr->data);
}
gboolean* bigdft_locreg_get_grid(const BigDFT_LocReg *glr, BigDFT_Grid gridType)
{
  gboolean *grid;
  guint orig = 0;
  double h_[3];
  double mult;
  BigDFT_Atoms *atoms = BIGDFT_ATOMS(glr);
  double *radii;

  grid = g_malloc(sizeof(gboolean) * (glr->n[0] + 1) * (glr->n[1] + 1) * (glr->n[2] + 1));
  if (atoms->geocode == 'F')
    h_[0] = atoms->alat[0] / glr->n[0];
  else
    h_[0] = atoms->alat[0] / (glr->n[0] + 1);
  if (atoms->geocode == 'F' || atoms->geocode == 'S')
    h_[1] = atoms->alat[1] / glr->n[1];
  else
    h_[1] = atoms->alat[1] / (glr->n[1] + 1);
  if (atoms->geocode == 'F')
    h_[2] = atoms->alat[2] / glr->n[2];
  else
    h_[2] = atoms->alat[2] / (glr->n[2] + 1);
  mult = (gridType == GRID_COARSE)?glr->crmult:glr->frmult;
  radii = (gridType == GRID_COARSE)?glr->radii:glr->radii + atoms->ntypes;

  FC_FUNC_(fill_logrid, FILL_LOGRID)(&atoms->geocode, glr->n, glr->n + 1, glr->n + 2,
                                     &orig, glr->n, &orig, glr->n + 1, &orig, glr->n + 2,
                                     &orig, &atoms->nat, &atoms->ntypes, atoms->iatype,
                                     atoms->rxyz.data, radii, &mult, h_, h_ + 1, h_ + 2,
                                     (int*)grid);

  return grid;
}
double* bigdft_locreg_convert_to_isf(const BigDFT_LocReg *glr, const double *psic)
{
  guint n;
  double *psir;

  n = glr->ni[0] * glr->ni[1] * glr->ni[2];
  psir = g_malloc(sizeof(double) * n);
  
  FC_FUNC_(wf_iorbp_to_psi, WF_IORBP_TO_PSI)(psir, psic, glr->data);
  
  return psir;
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

static void bigdft_lzd_init(BigDFT_Lzd *obj)
{
#ifdef HAVE_GLIB
  memset((void*)((char*)obj + sizeof(BigDFT_LocReg)), 0, sizeof(BigDFT_Lzd) - sizeof(BigDFT_LocReg));
#else
  memset(obj, 0, sizeof(BigDFT_Lzd));
#endif
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

  if (lzd->data)
    FC_FUNC_(lzd_free, LZD_FREE)(&lzd->data);

#ifdef HAVE_GLIB
  G_OBJECT_CLASS(bigdft_lzd_parent_class)->finalize(obj);
#endif
}

BigDFT_Lzd* bigdft_lzd_new()
{
  BigDFT_Lzd *lzd;

#ifdef HAVE_GLIB
  lzd = BIGDFT_LZD(g_object_new(BIGDFT_LZD_TYPE, NULL));
#else
  lzd = g_malloc(sizeof(BigDFT_Lzd));
  bigdft_lzd_init(lzd);
#endif
   
  FC_FUNC_(lzd_new, LZD_NEW)(&lzd->data);
  FC_FUNC_(lzd_init, LZD_INIT)(lzd->data, &lzd->parent.data);
  FC_FUNC_(glr_init, GLR_INIT)(lzd->parent.data, &lzd->parent.d);

  return lzd;
}
BigDFT_Lzd* bigdft_lzd_new_with_fortran(void *fortran_lzd)
{
  BigDFT_Lzd *lzd;

#ifdef HAVE_GLIB
  lzd = BIGDFT_LZD(g_object_new(BIGDFT_LZD_TYPE, NULL));
#else
  lzd = g_malloc(sizeof(BigDFT_Lzd));
  bigdft_lzd_init(lzd);
#endif

  lzd->data = fortran_lzd;
  FC_FUNC_(lzd_init, LZD_INIT)(lzd->data, &lzd->parent.data);
  FC_FUNC_(glr_init, GLR_INIT)(lzd->parent.data, &lzd->parent.d);

  return lzd;
}
BigDFT_Lzd* bigdft_lzd_new_from_fortran(void *fortran_lzd)
{
  BigDFT_Lzd *lzd;

#ifdef HAVE_GLIB
  lzd = BIGDFT_LZD(g_object_new(BIGDFT_LZD_TYPE, NULL));
#else
  lzd = g_malloc(sizeof(BigDFT_Lzd));
  bigdft_lzd_init(lzd);
#endif

  lzd->data = fortran_lzd;
  FC_FUNC_(lzd_get_data, LZD_GET_DATA)(lzd->data, &lzd->parent.data);
  FC_FUNC_(glr_get_data, GLR_GET_DATA)(lzd->parent.data, &lzd->parent.d);

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
void bigdft_lzd_set_size(BigDFT_Lzd *lzd, const double h[3],
                         double crmult, double frmult)
{
  bigdft_locreg_set_size(&lzd->parent, h, crmult, frmult);
  FC_FUNC_(lzd_set_hgrids, LZD_SET_HGRIDS)(lzd->data, lzd->parent.h);
  lzd->h[0] = lzd->parent.h[0];
  lzd->h[1] = lzd->parent.h[1];
  lzd->h[2] = lzd->parent.h[2];
}
void bigdft_lzd_copy_from_fortran(BigDFT_Lzd *lzd, const double *radii,
                                  double crmult, double frmult)
{
  FC_FUNC_(lzd_get_hgrids, LZD_GET_HGRIDS)(lzd->data, lzd->h);
  bigdft_locreg_copy_data(&lzd->parent, radii, lzd->h, crmult, frmult);
}
