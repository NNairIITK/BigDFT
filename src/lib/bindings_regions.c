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

static gboolean _print_error(const gchar *message)
{
  fprintf(stderr, "Error on %s.\n", message);
  return FALSE;
}

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
  /* Assign values. */
  glr->h[0] = h[0];
  glr->h[1] = h[1];
  glr->h[2] = h[2];

  glr->crmult = crmult;
  glr->frmult = frmult;
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
static void _locreg_copy_d(BigDFT_LocReg *glr)
{
  FC_FUNC_(glr_get_dimensions, GLR_GET_DIMENSIONS)(glr->data, glr->n, glr->ni, glr->ns, glr->nsi,
                                                   glr->nfl, glr->nfu, &glr->norb);
  FC_FUNC_(glr_get_locreg_data, GLR_GET_LOCREG_DATA)(glr->data, &glr->locrad, glr->locregCenter);
}
static void _locreg_copy_wfd(BigDFT_LocReg *glr)
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
void bigdft_locreg_set_d_dims(BigDFT_LocReg *lr, guint n[3], guint ni[3], guint ns[3],
                              guint nsi[3], guint nfl[3], guint nfu[3])
{
  FC_FUNC_(glr_set_dimensions, GLR_SET_DIMENSIONS)(lr->data, n, ni, ns, nsi, nfl, nfu);
  _locreg_copy_d(lr);
}
void bigdft_locreg_set_wfd_dims(BigDFT_LocReg *lr, guint nseg_c, guint nseg_f,
                                guint nvctr_c, guint nvctr_f)
{
  FC_FUNC_(glr_set_wfd_dims, GLR_SET_WFD_DIMS)(lr->data, &nseg_c, &nseg_f, &nvctr_c, &nvctr_f);
  _locreg_copy_wfd(lr);
}
void bigdft_locreg_init_d(BigDFT_LocReg *glr)
{
  int iproc = 1;

  /* Deallocate all previously allocated data. */
  FC_FUNC_(glr_empty, GLR_EMPTY)(glr->data);

  /* Set the new size. */
  FC_FUNC_(system_size, SYSTEM_SIZE)(&iproc, glr->parent.data, glr->parent.rxyz.data,
                                     glr->radii, &glr->crmult, &glr->frmult, glr->h, glr->h + 1,
                                     glr->h + 2, glr->data, glr->parent.shift);
  _locreg_copy_d(glr);
  /* Update Atoms accordingly. */
  FC_FUNC_(atoms_copy_alat, ATOMS_COPY_ALAT)(glr->parent.data, glr->parent.alat,
                                             glr->parent.alat + 1, glr->parent.alat + 2);
}
void bigdft_locreg_init_wfd(BigDFT_LocReg *glr)
{
  int iproc = 1;

  FC_FUNC_(glr_set_wave_descriptors,
           GLR_SET_WAVE_DESCRIPTORS)(&iproc, glr->h, glr->h + 1, glr->h + 2,
                                     glr->parent.data, glr->parent.rxyz.data, glr->radii,
                                     &glr->crmult, &glr->frmult, glr->data);
  _locreg_copy_wfd(glr);
}
void bigdft_locreg_init_bounds(BigDFT_LocReg *lr)
{
  FC_FUNC_(glr_set_bounds, GLR_SET_BOUNDS)(lr->data);
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
void bigdft_locreg_write_psi_compress(const BigDFT_LocReg *lr,
                                      guint unitwf, BigDFT_WfFileFormats format,
                                      gboolean linear, guint iorb, const guint n[3], const double *psic)
{
  guint useFormattedOutput = (format == BIGDFT_WF_FORMAT_PLAIN), confPotOrder = 4;
  double eval = 0., confPotprefac = 0.;

  if (linear)
    FC_FUNC_(writeonewave_linear, WRITEONEWAVE_LINEAR)
      (&unitwf, &useFormattedOutput, &iorb, n, n + 1, n + 2, lr->h, lr->h + 1, lr->h + 2,
       &lr->locregCenter, &lr->locrad, &confPotOrder, &confPotprefac,
       &lr->parent.nat, lr->parent.rxyz.data,
       &lr->nseg_c, &lr->nvctr_c, lr->keyglob, lr->keyvglob,
       &lr->nseg_f, &lr->nvctr_f, lr->keyglob + 2 * lr->nseg_c, lr->keyvglob + lr->nseg_c,
       psic, psic + lr->nvctr_c, &eval);
  else
    FC_FUNC(writeonewave, WRITEONEWAVE)
      (&unitwf, &useFormattedOutput, &iorb, n, n + 1, n + 2,
       lr->h, lr->h + 1, lr->h + 2, &lr->parent.nat, lr->parent.rxyz.data,
       &lr->nseg_c, &lr->nvctr_c, lr->keyglob, lr->keyvglob,
       &lr->nseg_f, &lr->nvctr_f, lr->keyglob + 2 * lr->nseg_c, lr->keyvglob + lr->nseg_c,
       psic, psic + lr->nvctr_c, &eval);
}

gboolean bigdft_locreg_check(const BigDFT_LocReg *glr)
{
  guint n[3], ni[3], ns[3], nsi[3], nfl[3], nfu[3], norb;
  BigDFT_LocReg ref;

  /* Check the consistency of everything with Fortran values. */
  FC_FUNC_(glr_get_dimensions, GLR_GET_DIMENSIONS)(glr->data, n, ni, ns, nsi, nfl, nfu, &norb);
  if (n[0] != glr->n[0] || n[1] != glr->n[1] || n[2] != glr->n[2])
    return _print_error("n");
  if (ni[0] != glr->ni[0] || ni[1] != glr->ni[1] || ni[2] != glr->ni[2])
    return _print_error("ni");
  if (ns[0] != glr->ns[0] || ns[1] != glr->ns[1] || ns[2] != glr->ns[2])
    return _print_error("ns");
  if (nsi[0] != glr->nsi[0] || nsi[1] != glr->nsi[1] || nsi[2] != glr->nsi[2])
    return _print_error("nsi");
  if (nfl[0] != glr->nfl[0] || nfl[1] != glr->nfl[1] || nfl[2] != glr->nfl[2])
    return _print_error("nfl");
  if (nfu[0] != glr->nfu[0] || nfu[1] != glr->nfu[1] || nfu[2] != glr->nfu[2])
    return _print_error("nfu");

  ref.wfd = glr->wfd;
  _locreg_copy_wfd(&ref);
  if (ref.nvctr_c != glr->nvctr_c || ref.nvctr_f != glr->nvctr_f)
    return _print_error("nvctr");
  if (ref.nseg_c != glr->nseg_c || ref.nseg_f != glr->nseg_f)
    return _print_error("nseg");
  if (ref.keygloc != glr->keygloc || ref.keyvloc != glr->keyvloc)
    return _print_error("keyloc");
  if (ref.keyglob != glr->keyglob || ref.keyvglob != glr->keyvglob)
    return _print_error("keyglob");

  return TRUE;
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
enum {
  LZD_DEFINED_SIGNAL,
  LAST_LZD_SIGNAL
};

G_DEFINE_TYPE(BigDFT_Lzd, bigdft_lzd, BIGDFT_LOCREG_TYPE)

static guint bigdft_lzd_signals[LAST_LZD_SIGNAL] = { 0 };

static void bigdft_lzd_class_init(BigDFT_LzdClass *klass)
{
  /* Connect the overloading methods. */
  G_OBJECT_CLASS(klass)->dispose      = bigdft_lzd_dispose;
  G_OBJECT_CLASS(klass)->finalize     = bigdft_lzd_finalize;
  /* G_OBJECT_CLASS(klass)->set_property = visu_data_set_property; */
  /* G_OBJECT_CLASS(klass)->get_property = visu_data_get_property; */

  bigdft_lzd_signals[LZD_DEFINED_SIGNAL] =
    g_signal_new("defined", G_TYPE_FROM_CLASS(klass),
                 G_SIGNAL_RUN_LAST | G_SIGNAL_NO_RECURSE | G_SIGNAL_NO_HOOKS,
		 0, NULL, NULL, g_cclosure_marshal_VOID__VOID,
                 G_TYPE_NONE, 0, NULL);
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
static void _free_llr(BigDFT_Lzd *lzd, gboolean separate)
{
  guint i;

  if (lzd->Llr)
    {
      /* We create a Fortran copy of the nderlying structure for
         each locreg and dec ref the C wrapper. Caller is responsible
         to destroy original Fortran underlying structure. */
      for (i = 0; i < lzd->nlr; i++)
        {
          if (separate)
            {
              FC_FUNC_(glr_copy, GLR_COPY)(&lzd->Llr[i]->data, &lzd->Llr[i]->d,
                                           &lzd->Llr[i]->wfd, lzd->Llr[i]->data);
              _locreg_copy_wfd(lzd->Llr[i]);
            }
          else
            {
              lzd->Llr[i]->data = (gpointer)0;
            }
          /* Currently the underlying Atom structure is not copied. */
          F90_2D_POINTER_INIT(&lzd->Llr[i]->parent.rxyz);
          lzd->Llr[i]->parent.data = (gpointer)0;
          /* Finally, we unref the locreg. */
          bigdft_locreg_free(lzd->Llr[i]);
        }
      g_free(lzd->Llr);
    }
  lzd->nlr = 0;
}
static void bigdft_lzd_finalize(GObject *obj)
{
  BigDFT_Lzd *lzd = BIGDFT_LZD(obj);

  if (lzd->data)
    FC_FUNC_(lzd_free, LZD_FREE)(&lzd->data);

  _free_llr(lzd, FALSE);

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
static void _allocate_llr(BigDFT_Lzd *lzd)
{
  guint i, j;
  gpointer llr;

  if (lzd->nlr > 0)
    {
      lzd->Llr = g_malloc(sizeof(BigDFT_LocReg*) * lzd->nlr);
      for (i = 0; i < lzd->nlr; i++)
        {
          j = i + 1;
          FC_FUNC_(lzd_get_llr, LZD_GET_LLR)(lzd->data, &j, &llr);
          lzd->Llr[i] = bigdft_locreg_new_from_fortran(llr);
          /* Copy the atoms data. */
          lzd->Llr[i]->parent.data = lzd->parent.parent.data;
          if (lzd->Llr[i]->parent.data)
            bigdft_atoms_copy_from_fortran(&lzd->Llr[i]->parent);
          memcpy(&lzd->Llr[i]->parent.rxyz, &lzd->parent.parent.rxyz,
                 sizeof(f90_pointer_double_2D) + 3 * sizeof(double));
          /* Copy some building values. */
          bigdft_locreg_set_radii(lzd->Llr[i], lzd->parent.radii);
          bigdft_locreg_set_size(lzd->Llr[i], lzd->parent.h, lzd->parent.crmult, lzd->parent.frmult);
        }
    }
}
static void _lzd_wrap_llr(BigDFT_Lzd *lzd)
{
  guint i;

  /* Get the llr array. */
  FC_FUNC_(lzd_copy_data, LZD_COPY_DATA)(lzd->data, &lzd->nlr);
  _allocate_llr(lzd);
  for (i = 0; i < lzd->nlr; i++)
    {
      /* Copy the locreg data. */
      _locreg_copy_d(lzd->Llr[i]);
      _locreg_copy_wfd(lzd->Llr[i]);
    }
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
  _lzd_wrap_llr(lzd);

#ifdef HAVE_GLIB
  g_signal_emit(G_OBJECT(lzd), bigdft_lzd_signals[LZD_DEFINED_SIGNAL],
                0 /* details */, NULL);
#endif

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
void bigdft_lzd_init_d(BigDFT_Lzd *lzd)
{
  bigdft_locreg_init_d(&lzd->parent);
  FC_FUNC_(lzd_set_hgrids, LZD_SET_HGRIDS)(lzd->data, lzd->parent.h);
}
void bigdft_lzd_set_irreductible_zone(BigDFT_Lzd *lzd, guint nspin)
{
  FC_FUNC_(symmetry_set_irreductible_zone, SYMMETRY_SET_IRREDUCTIBLE_ZONE)
    (lzd->parent.parent.sym, &lzd->parent.parent.geocode,
     lzd->parent.ni, lzd->parent.ni + 1, lzd->parent.ni + 2, &nspin, 1);
}
void bigdft_lzd_copy_from_fortran(BigDFT_Lzd *lzd, const double *radii,
                                  double crmult, double frmult)
{
  FC_FUNC_(lzd_get_hgrids, LZD_GET_HGRIDS)(lzd->data, lzd->parent.h);
  bigdft_locreg_set_radii(&lzd->parent, radii);
  bigdft_locreg_set_size(&lzd->parent, lzd->parent.h, crmult, frmult);
  _locreg_copy_d(&lzd->parent);
}
void bigdft_lzd_define(BigDFT_Lzd *lzd, guint type,
                       BigDFT_Orbs *orbs, guint iproc, guint nproc)
{
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
  _lzd_wrap_llr(lzd);

#ifdef HAVE_GLIB
  g_signal_emit(G_OBJECT(lzd), bigdft_lzd_signals[LZD_DEFINED_SIGNAL],
                0 /* details */, NULL);
#endif
}
void bigdft_lzd_emit_defined(BigDFT_Lzd *lzd)
{
  guint ilr, j;

  /* First, sync. */
  FC_FUNC_(lzd_copy_data, LZD_COPY_DATA)(lzd->data, &lzd->nlr);
  _locreg_copy_d(&lzd->parent);
  _locreg_copy_wfd(&lzd->parent);
  for (ilr = 0; ilr < lzd->nlr; ilr++)
    {
      j = ilr + 1;
      FC_FUNC_(lzd_get_llr, LZD_GET_LLR)(lzd->data, &j, &lzd->Llr[ilr]->data);
      FC_FUNC_(glr_get_data, GLR_GET_DATA)(lzd->Llr[ilr]->data, &lzd->Llr[ilr]->d, &lzd->Llr[ilr]->wfd);
      _locreg_copy_d(lzd->Llr[ilr]);
      _locreg_copy_wfd(lzd->Llr[ilr]);
    }
  
#ifdef HAVE_GLIB
  g_signal_emit(G_OBJECT(lzd), bigdft_lzd_signals[LZD_DEFINED_SIGNAL],
                0 /* details */, NULL);
#endif
}
void bigdft_lzd_set_n_locreg(BigDFT_Lzd *lzd, guint nlr)
{
  _free_llr(lzd, TRUE);

  lzd->nlr = nlr;
  FC_FUNC_(lzd_set_nlr, LZD_SET_NLR)(lzd->data, &nlr, &lzd->parent.parent.geocode, 1);
  _allocate_llr(lzd);
}

gboolean bigdft_lzd_check(const BigDFT_Lzd *lzd)
{
  guint nlr, ilr;

  /* Check the consistency of everything with Fortran values. */

  /* nlr. */
  FC_FUNC_(lzd_copy_data, LZD_COPY_DATA)(lzd->data, &nlr);
  if (nlr != lzd->nlr)
    return _print_error("nlr");

  for (ilr = 0; ilr < lzd->nlr; ilr++)
    if (!bigdft_locreg_check(lzd->Llr[ilr]))
      return FALSE;
  
  return TRUE;
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
