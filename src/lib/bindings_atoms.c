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

static void bigdft_atoms_dispose(GObject *atoms);
static void bigdft_atoms_finalize(GObject *atoms);

#ifdef HAVE_GLIB
G_DEFINE_TYPE(BigDFT_Atoms, bigdft_atoms, G_TYPE_OBJECT)

static void bigdft_atoms_class_init(BigDFT_AtomsClass *klass)
{
  /* Connect the overloading methods. */
  G_OBJECT_CLASS(klass)->dispose      = bigdft_atoms_dispose;
  G_OBJECT_CLASS(klass)->finalize     = bigdft_atoms_finalize;
  /* G_OBJECT_CLASS(klass)->set_property = visu_data_set_property; */
  /* G_OBJECT_CLASS(klass)->get_property = visu_data_get_property; */
}
#endif

static void bigdft_atoms_free_additional(BigDFT_Atoms *atoms)
{
  guint i;

  FC_FUNC_(deallocate_double_2d, DEALLOCATE_DOUBLE_2D)(&atoms->rxyz);
  if (atoms->atomnames)
    {
      for (i = 0; i < atoms->ntypes; i++)
        if (atoms->atomnames[i])
          g_free(atoms->atomnames[i]);
      g_free(atoms->atomnames);
    }
}
static void bigdft_atoms_init(BigDFT_Atoms *obj)
{
#ifdef HAVE_GLIB
  memset((void*)((char*)obj + sizeof(GObject)), 0, sizeof(BigDFT_Atoms) - sizeof(GObject));
#else
  memset(obj, 0, sizeof(BigDFT_Atoms));
#endif
  F90_2D_POINTER_INIT(&obj->rxyz);
  FC_FUNC_(atoms_new, ATOMS_NEW)(&obj->data, &obj->sym);
}
static void bigdft_atoms_dispose(GObject *obj)
{
#ifdef HAVE_GLIB
  BigDFT_Atoms *atoms = BIGDFT_ATOMS(obj);

  if (atoms->dispose_has_run)
    return;
  atoms->dispose_has_run = TRUE;

  /* Chain up to the parent class */
  G_OBJECT_CLASS(bigdft_atoms_parent_class)->dispose(obj);
#endif
}
static void bigdft_atoms_finalize(GObject *obj)
{
  BigDFT_Atoms *atoms = BIGDFT_ATOMS(obj);

  FC_FUNC_(atoms_free, ATOMS_FREE)(&atoms->data);
  bigdft_atoms_free_additional(atoms);

#ifdef HAVE_GLIB
  G_OBJECT_CLASS(bigdft_atoms_parent_class)->finalize(obj);
#endif
}
static void bigdft_atoms_get_nat_arrays(BigDFT_Atoms *atoms)
{
  GET_ATTR_UINT  (atoms, ATOMS, iatype,   IATYPE);
  GET_ATTR_UINT  (atoms, ATOMS, iasctype, IASCTYPE);
  GET_ATTR_UINT  (atoms, ATOMS, natpol,   NATPOL);
  GET_ATTR_INT   (atoms, ATOMS, ifrztyp,  IFRZTYP);
  GET_ATTR_DBL   (atoms, ATOMS, amu,      AMU);
  GET_ATTR_DBL_2D(atoms, ATOMS, aocc,     AOCC);
}
static void bigdft_atoms_get_ntypes_arrays(BigDFT_Atoms *atoms)
{
  GET_ATTR_UINT  (atoms, ATOMS, nelpsp,     NELPSP);
  GET_ATTR_UINT  (atoms, ATOMS, npspcode,   NPSPCODE);
  GET_ATTR_UINT  (atoms, ATOMS, nzatom,     NZATOM);
  GET_ATTR_INT   (atoms, ATOMS, nlcc_ngv,   NLCC_NGV);
  GET_ATTR_INT   (atoms, ATOMS, nlcc_ngc,   NLCC_NGC);
  GET_ATTR_INT   (atoms, ATOMS, ixcpsp,     IXCPSP);
  GET_ATTR_DBL_2D(atoms, ATOMS, radii_cf,   RADII_CF);
  GET_ATTR_DBL_3D(atoms, ATOMS, psppar,     PSPPAR);
  GET_ATTR_DBL_2D(atoms, ATOMS, nlccpar,    NLCCPAR);
  GET_ATTR_DBL_2D(atoms, ATOMS, ig_nlccpar, IG_NLCCPAR);
}

BigDFT_Atoms* bigdft_atoms_new()
{
  BigDFT_Atoms *atoms;

#ifdef HAVE_GLIB
  atoms = BIGDFT_ATOMS(g_object_new(BIGDFT_ATOMS_TYPE, NULL));
#else
  atoms = g_malloc(sizeof(BigDFT_Atoms));
  bigdft_atoms_init(atoms);
#endif

  return atoms;
}
void bigdft_atoms_free(BigDFT_Atoms *atoms)
{
#ifdef HAVE_GLIB
  g_object_unref(G_OBJECT(atoms));
#else
  bigdft_atoms_finalize(atoms);
  g_free(atoms);
#endif
}
BigDFT_Atoms* bigdft_atoms_new_from_file(const gchar *filename)
{
  BigDFT_Atoms *atoms;
  guint lstat, ln, i, j;
  gchar str[20];

#ifdef HAVE_GLIB
  atoms = BIGDFT_ATOMS(g_object_new(BIGDFT_ATOMS_TYPE, NULL));
#else
  atoms = g_malloc(sizeof(BigDFT_Atoms));
  bigdft_atoms_init(atoms);
#endif
  
  if (!bigdft_atoms_set_structure_from_file(atoms, filename))
    {
      bigdft_atoms_free(atoms);
      atoms = (BigDFT_Atoms*)0;
    }

  return atoms;
}
void bigdft_atoms_set_n_atoms(BigDFT_Atoms *atoms, guint nat)
{
  FC_FUNC_(atoms_set_n_atoms, ATOMS_SET_N_ATOMS)(atoms->data, &atoms->rxyz,
                                                 (int*)(&nat));
  atoms->nat = nat;
  bigdft_atoms_get_nat_arrays(atoms);
}
void bigdft_atoms_set_n_types(BigDFT_Atoms *atoms, guint ntypes)
{
  FC_FUNC_(atoms_set_n_types, ATOMS_SET_N_TYPES)(atoms->data, (int*)(&ntypes));
  atoms->ntypes = ntypes;
  bigdft_atoms_get_ntypes_arrays(atoms);
  atoms->atomnames = g_malloc(sizeof(gchar*) * ntypes);
  memset(atoms->atomnames, 0, sizeof(gchar*) * ntypes);
}
void bigdft_atoms_sync(BigDFT_Atoms *atoms)
{
  gchar format[5], units[20];
  guint i, j;

  memset(format, ' ', 5);
  if (atoms->format)
    memcpy(format, atoms->format, strlen(atoms->format));
  memset(units, ' ', 20);
  if (atoms->units)
    memcpy(units, atoms->units, strlen(atoms->units));
  FC_FUNC_(atoms_sync, ATOMS_SYNC)(atoms->data, atoms->alat, atoms->alat + 1, atoms->alat + 2,
                                   &atoms->geocode, format, units);
  for (i = 0; i < atoms->ntypes; i++)
    {
      j = i + 1;
      memset(units, ' ', 20);
      if (atoms->atomnames[i])
        memcpy(units, atoms->atomnames[i], strlen(atoms->atomnames[i]));
      FC_FUNC_(atoms_set_name, ATOMS_SET_NAME)(atoms->data, (int*)(&j), units);
    }
}
gboolean bigdft_atoms_set_structure_from_file(BigDFT_Atoms *atoms, const gchar *filename)
{
  guint lstat, ln, i, j;
  gchar str[20];

  FC_FUNC_(atoms_empty, ATOMS_EMPTY)(atoms->data);
  bigdft_atoms_free_additional(atoms);

  ln = strlen(filename);
  FC_FUNC_(atoms_set_from_file, ATOMS_SET_FROM_FILE)((int*)(&lstat), atoms->data,
                                                     &atoms->rxyz, filename, (int*)(&ln));

  if (!lstat)
    return FALSE;
  else
    {
      FC_FUNC_(atoms_copy_geometry_data, ATOMS_COPY_GEOMETRY_DATA)
        (atoms->data, &atoms->geocode, atoms->format, atoms->units);
      FC_FUNC_(atoms_copy_alat, ATOMS_COPY_ALAT)(atoms->data, atoms->alat,
                                                 atoms->alat + 1, atoms->alat + 2);
      FC_FUNC_(atoms_copy_nat, ATOMS_COPY_NAT)(atoms->data, (int*)(&atoms->nat));
      bigdft_atoms_get_nat_arrays(atoms);
      FC_FUNC_(atoms_copy_ntypes, ATOMS_COPY_NTYPES)(atoms->data,
                                                     (int*)(&atoms->ntypes));
      bigdft_atoms_get_ntypes_arrays(atoms);
      atoms->atomnames = g_malloc(sizeof(gchar*) * atoms->ntypes);
      for (i = 0; i < atoms->ntypes; i++)
        {
          j = i + 1;
          FC_FUNC_(atoms_copy_name, ATOMS_COPY_NAME)(atoms->data,
                                                     (int*)(&j), str, (int*)(&ln));
          atoms->atomnames[i] = g_malloc(sizeof(gchar) * (ln + 1));
          memcpy(atoms->atomnames[i], str, sizeof(gchar) * ln);
          atoms->atomnames[i][ln] = '\0';
        }
    }
  return TRUE;
}

void bigdft_atoms_set_psp(BigDFT_Atoms *atoms, int ixc, guint nspin, const gchar *occup)
{
  int verb = 0;
  guint ln;

  FC_FUNC_(init_atomic_values, INIT_ATOMIC_VALUES)(&verb, atoms->data, &ixc);
  FC_FUNC_(atoms_copy_psp_data, ATOMS_COPY_PSP_DATA)
    (atoms->data, (int*)(&atoms->natsc), (int*)(&atoms->donlcc));
  ln = (occup)?strlen(occup):0;
  FC_FUNC_(atoms_read_variables, ATOMS_READ_VARIABLES)(atoms->data, &nspin, occup, &ln);
}

void bigdft_atoms_set_symmetries(BigDFT_Atoms *atoms, gboolean active,
                                 double tol, double elecfield[3])
{
  int disable;

  disable = (!active);
  FC_FUNC_(atoms_set_symmetries, ATOMS_SET_SYMMETRIES)(atoms->data, atoms->rxyz.data,
                                                       &disable, &tol, elecfield);
}

void bigdft_atoms_set_displacement(BigDFT_Atoms *atoms, double randdis)
{
  FC_FUNC_(atoms_set_displacement, ATOMS_SET_DISPLACEMENT)(atoms->data, atoms->rxyz.data, &randdis);
}

double* bigdft_atoms_get_radii(const BigDFT_Atoms *atoms, double crmult,
                               double frmult, double projrad)
{
  double *radii_cf;
  double crmult_, frmult_, projrad_;

  radii_cf = g_malloc(sizeof(double) * 3 * atoms->ntypes);
  crmult_  = (crmult <= 0.)?5.:crmult;
  frmult_  = (frmult <= 0.)?8.:frmult;
  projrad_ = (projrad <= 0.)?15.:projrad;
  FC_FUNC_(read_radii_variables, READ_RADII_VARIABLES)(atoms->data, radii_cf,
                                                       &crmult_, &frmult_, &projrad_);
  return radii_cf;
}

void bigdft_atoms_write(const BigDFT_Atoms *atoms, const gchar *filename)
{
  gchar *comment;
  int ln, ln2;
  f90_pointer_double forces;

  if (atoms->comment)
    {
      comment = atoms->comment;
      ln = strlen(atoms->comment);
    }
  else
    {
      comment = " ";
      ln = 1;
    }
  F90_1D_POINTER_INIT(&forces);
  ln2 = strlen(filename);
  FC_FUNC_(atoms_write, ATOMS_WRITE)(atoms->data, filename, &ln2, atoms->rxyz.data,
                                     &forces, &atoms->energy, comment, &ln);
}
