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

  FC_FUNC_(atoms_new, ATOMS_NEW)(&glr->parent.data, &glr->parent.sym);
  FC_FUNC_(glr_new, GLR_NEW)(&glr->data);
  FC_FUNC_(glr_init, GLR_INIT)(glr->data, &glr->d, &glr->wfd);

  return glr;
}
BigDFT_LocReg* bigdft_locreg_new_from_fortran(void *glr_fortran)
{
  BigDFT_LocReg *glr;

#ifdef HAVE_GLIB
  glr = BIGDFT_LOCREG(g_object_new(BIGDFT_LOCREG_TYPE, NULL));
#else
  glr = g_malloc(sizeof(BigDFT_LocReg));
  bigdft_locreg_init(glr);
#endif
  glr->data = glr_fortran;

  FC_FUNC_(glr_get_data, GLR_GET_DATA)(glr->data, &glr->d, &glr->wfd);

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
  FC_FUNC_(glr_get_dimensions, GLR_GET_DIMENSIONS)(glr->data, glr->n, glr->ni,
                                                   glr->ns, glr->nsi, &glr->norb);
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
  FC_FUNC_(glr_get_dimensions, GLR_GET_DIMENSIONS)(glr->data, glr->n, glr->ni,
                                                   glr->ns, glr->nsi, &glr->norb);
  bigdft_locreg_set_radii(glr, radii);
  glr->h[0] = h[0];
  glr->h[1] = h[1];
  glr->h[2] = h[2];
  glr->crmult = crmult;
  glr->frmult = frmult;
}
static void bigdft_locreg_copy_wfd(BigDFT_LocReg *glr)
{
  f90_pointer_int    tmp, tmp2;
  f90_pointer_int_2D tmp2D, tmp2D2;

  memset(&tmp2D, 0, sizeof(f90_pointer_int_2D));
  memset(&tmp2D2, 0, sizeof(f90_pointer_int_2D));
  memset(&tmp, 0, sizeof(f90_pointer_int));
  memset(&tmp2, 0, sizeof(f90_pointer_int));
  FC_FUNC_(glr_wfd_get_data, GLR_WFD_GET_DATA)(glr->wfd, &glr->nvctr_c, &glr->nvctr_f,
                                               &glr->nseg_c, &glr->nseg_f, &tmp2D, &tmp2D2,
                                               &tmp, &tmp2);
  glr->keyglob  = (guint*)tmp2D.data;
  glr->keygloc  = (guint*)tmp2D2.data;
  glr->keyvglob = (guint*)tmp.data;
  glr->keyvloc  = (guint*)tmp2.data;
}
void bigdft_locreg_set_wave_descriptors(BigDFT_LocReg *glr)
{
  int iproc = 1;

  FC_FUNC_(glr_set_wave_descriptors,
           GLR_SET_WAVE_DESCRIPTORS)(&iproc, glr->h, glr->h + 1, glr->h + 2,
                                     glr->parent.data, glr->parent.rxyz.data, glr->radii,
                                     &glr->crmult, &glr->frmult, glr->data);
  bigdft_locreg_copy_wfd(glr);
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
gboolean bigdft_locreg_iter_new(const BigDFT_LocReg *glr, BigDFT_LocRegIter *iter, BigDFT_Grid gridType)
{
  memset(iter, 0, sizeof(BigDFT_LocRegIter));
  iter->nseg = (gridType == GRID_COARSE)?glr->nseg_c:glr->nseg_c + glr->nseg_f;
  iter->iseg = (gridType == GRID_COARSE)?(guint)-1:glr->nseg_c - 1;
  iter->glr  = glr;

  switch (glr->parent.geocode)
    {
    case 'S':
      iter->grid[0] = glr->n[0] + 1;
      iter->grid[1] = glr->n[1];
      iter->grid[2] = glr->n[2] + 1;
      break;
    case 'F':
      iter->grid[0] = glr->n[0];
      iter->grid[1] = glr->n[1];
      iter->grid[2] = glr->n[2];
      break;
    default:
      fprintf(stderr, "WARNING: unknown geocode.\n");
    case 'P':
      iter->grid[0] = glr->n[0] + 1;
      iter->grid[1] = glr->n[1] + 1;
      iter->grid[2] = glr->n[2] + 1;
      break;
    }

  return bigdft_locreg_iter_next(iter);
}
gboolean bigdft_locreg_iter_next(BigDFT_LocRegIter *iter)
{
  guint jj, j0, j1, ii;

  iter->iseg += 1;
  if (iter->iseg >= iter->nseg)
    return FALSE;

  jj = iter->glr->keyvloc[iter->iseg];
  j0 = iter->glr->keygloc[iter->iseg * 2];
  j1 = iter->glr->keygloc[iter->iseg * 2 + 1];
      
  ii = j0 - 1;
  iter->i3 = ii / ((iter->glr->n[0] + 1) * (iter->glr->n[1] + 1));
  ii = ii - iter->i3 * (iter->glr->n[0] + 1) * (iter->glr->n[1] + 1);
  iter->i2 = ii / (iter->glr->n[0] + 1);
  iter->i0 = ii - iter->i2 * (iter->glr->n[0] + 1);
  iter->i1 = iter->i0 + j1 - j0;

  iter->z  = (double)iter->i3 / (double)iter->grid[2];
  iter->y  = (double)iter->i2 / (double)iter->grid[1];
  iter->x0 = (double)iter->i0 / (double)iter->grid[0];
  iter->x1 = (double)iter->i1 / (double)iter->grid[0];

  return TRUE;
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
  guint i;

  if (lzd->data)
    FC_FUNC_(lzd_free, LZD_FREE)(&lzd->data);

  if (lzd->Llr)
    {
      for (i = 0; i < lzd->nlr; i++)
        {
          /* We free only the C wrapper. */
          lzd->Llr[i]->data = (gpointer)0;
          lzd->Llr[i]->parent.data = (gpointer)0;
          bigdft_locreg_free(lzd->Llr[i]);
        }
      g_free(lzd->Llr);
    }

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
  FC_FUNC_(glr_init, GLR_INIT)(lzd->parent.data, &lzd->parent.d, &lzd->parent.wfd);
  FC_FUNC_(atoms_new, ATOMS_NEW)(&lzd->parent.parent.data, &lzd->parent.parent.sym);

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
  FC_FUNC_(glr_init, GLR_INIT)(lzd->parent.data, &lzd->parent.d, &lzd->parent.wfd);
  FC_FUNC_(atoms_new, ATOMS_NEW)(&lzd->parent.parent.data, &lzd->parent.parent.sym);

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
  FC_FUNC_(glr_get_data, GLR_GET_DATA)(lzd->parent.data, &lzd->parent.d, &lzd->parent.wfd);
  FC_FUNC_(lzd_copy_data, LZD_COPY_DATA)(lzd->data, &lzd->nlr);

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
void bigdft_lzd_set_irreductible_zone(BigDFT_Lzd *lzd, guint nspin)
{
  FC_FUNC_(symmetry_set_irreductible_zone, SYMMETRY_SET_IRREDUCTIBLE_ZONE)
    (lzd->parent.parent.sym, &lzd->parent.parent.geocode,
     lzd->parent.ni, lzd->parent.ni + 1, lzd->parent.ni + 2, &nspin);
}
void bigdft_lzd_copy_from_fortran(BigDFT_Lzd *lzd, const double *radii,
                                  double crmult, double frmult)
{
  FC_FUNC_(lzd_get_hgrids, LZD_GET_HGRIDS)(lzd->data, lzd->h);
  bigdft_locreg_copy_data(&lzd->parent, radii, lzd->h, crmult, frmult);
}
void bigdft_lzd_define(BigDFT_Lzd *lzd, guint type,
                       BigDFT_Orbs *orbs, guint iproc, guint nproc)
{
  guint i, j;
  gpointer llr;
  guint withderorbs = 0;
  void *dorbs;

  FC_FUNC_(lzd_empty, LZD_EMPTY)(lzd->data);

  FC_FUNC_(check_linear_and_create_lzd, CHECK_LINEAR_AND_CREATE_LZD)
    (&iproc, &nproc, &type, lzd->data, lzd->parent.parent.data,
     orbs->data, &orbs->nspin, lzd->parent.parent.rxyz.data);
  if (bigdft_orbs_get_linear(orbs))
    {
      FC_FUNC_(lzd_empty, LZD_EMPTY)(lzd->data);
      FC_FUNC_(orbs_new, ORBS_NEW)(&dorbs);
      FC_FUNC_(lzd_init_llr, LZD_INIT_LLR)
        (&iproc, &nproc, orbs->in->data, lzd->parent.parent.data, lzd->parent.parent.rxyz.data,
         orbs->data, dorbs, &withderorbs, lzd->data);
      GET_ATTR_UINT(orbs, ORBS, inwhichlocreg, INWHICHLOCREG);
      GET_ATTR_UINT(orbs, ORBS, onwhichmpi,    ONWHICHMPI);
      GET_ATTR_UINT(orbs, ORBS, onwhichatom,   ONWHICHATOM);
      FC_FUNC_(orbs_free, ORBS_FREE)(&dorbs);
    }

  /* Get the llr array. */
  FC_FUNC_(lzd_copy_data, LZD_COPY_DATA)(lzd->data, &lzd->nlr);
  if (lzd->nlr > 0)
    {
      lzd->Llr = g_malloc(sizeof(BigDFT_LocReg*) * lzd->nlr);
      for (i = 0; i < lzd->nlr; i++)
        {
          j = i + 1;
          FC_FUNC_(lzd_get_llr, LZD_GET_LLR)(lzd->data, &j, &llr);
          lzd->Llr[i] = bigdft_locreg_new_from_fortran(llr);
          lzd->Llr[i]->parent.data = lzd->parent.parent.data;
          bigdft_atoms_copy_from_fortran(&lzd->Llr[i]->parent);
          bigdft_locreg_copy_data(lzd->Llr[i], lzd->parent.radii, lzd->h,
                                  lzd->parent.crmult, lzd->parent.frmult);
          bigdft_locreg_copy_wfd(lzd->Llr[i]);
        }
    }
}
gboolean bigdft_lzd_iter_new(const BigDFT_Lzd *lzd, BigDFT_LocRegIter *iter,
                             BigDFT_Grid gridType, guint ilr)
{
  memset(iter, 0, sizeof(BigDFT_LocRegIter));
  iter->glr  = (ilr)?lzd->Llr[ilr - 1]:&lzd->parent;
  iter->nseg = (gridType == GRID_COARSE)?iter->glr->nseg_c:iter->glr->nseg_c + iter->glr->nseg_f;
  iter->iseg = (gridType == GRID_COARSE)?(guint)-1:iter->glr->nseg_c - 1;

  switch (lzd->parent.parent.geocode)
    {
    case 'S':
      iter->grid[0] = lzd->parent.n[0] + 1;
      iter->grid[1] = lzd->parent.n[1];
      iter->grid[2] = lzd->parent.n[2] + 1;
      break;
    case 'F':
      iter->grid[0] = lzd->parent.n[0];
      iter->grid[1] = lzd->parent.n[1];
      iter->grid[2] = lzd->parent.n[2];
      break;
    default:
      fprintf(stderr, "WARNING: unknown geocode.\n");
    case 'P':
      iter->grid[0] = lzd->parent.n[0] + 1;
      iter->grid[1] = lzd->parent.n[1] + 1;
      iter->grid[2] = lzd->parent.n[2] + 1;
      break;
    }

  return bigdft_lzd_iter_next(iter);
}
gboolean bigdft_lzd_iter_next(BigDFT_LocRegIter *iter)
{
  if (!bigdft_locreg_iter_next(iter))
    return FALSE;
  
  iter->z  += (double)iter->glr->ns[2] / (double)iter->grid[2];
  iter->y  += (double)iter->glr->ns[1] / (double)iter->grid[1];
  iter->x0 += (double)iter->glr->ns[0] / (double)iter->grid[0];
  iter->x1 += (double)iter->glr->ns[0] / (double)iter->grid[0];

  return TRUE;
}
