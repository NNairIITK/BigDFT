#include <config.h>

#include "bigdft.h"
#include "bindings.h"
#include "bindings_api.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

static void bigdft_atoms_dispose(GObject *atoms);
static void bigdft_atoms_finalize(GObject *atoms);

/**
 * BigDFT_Atoms:
 * @parent: its #GObject parent.
 * @dispose_has_run: a private internal flag.
 * @geocode: a flag for the boundary conditions.
 * @format: the format to represent these atoms on disk.
 * @units: the unit of the stored data on disk.
 * @ntypes:
 * @nat:
 * @natsc:
 * @atomnames: (array zero-terminated=1):
 * @donlcc:
 * @iatype:
 */

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
  G_OBJECT(obj)->ref_count = 1;
#endif
}
static void bigdft_atoms_dispose(GObject *obj)
{
  BigDFT_Atoms *atoms = BIGDFT_ATOMS(obj);

  if (atoms->dispose_has_run)
    return;
  atoms->dispose_has_run = TRUE;

#ifdef HAVE_GLIB
  /* Chain up to the parent class */
  G_OBJECT_CLASS(bigdft_atoms_parent_class)->dispose(obj);
#endif
}
static void bigdft_atoms_finalize(GObject *obj)
{
  BigDFT_Atoms *atoms = BIGDFT_ATOMS(obj);

  if (F_TYPE(atoms->data))
    FC_FUNC_(atoms_free, ATOMS_FREE)(&atoms->data);
  bigdft_atoms_free_additional(atoms);

#ifdef HAVE_GLIB
  G_OBJECT_CLASS(bigdft_atoms_parent_class)->finalize(obj);
#endif
}
void bigdft_atoms_get_nat_arrays(BigDFT_Atoms *atoms)
{
  GET_ATTR_UINT  (atoms, ATOMS, iatype,   IATYPE);
  GET_ATTR_UINT  (atoms, ATOMS, iasctype, IASCTYPE);
  GET_ATTR_UINT  (atoms, ATOMS, natpol,   NATPOL);
  GET_ATTR_INT   (atoms, ATOMS, ifrztyp,  IFRZTYP);
  GET_ATTR_DBL   (atoms, ATOMS, amu,      AMU);
  GET_ATTR_DBL_2D(atoms, ATOMS, aocc,     AOCC);
  GET_ATTR_DBL_2D(atoms, ATOMS, rxyz,     RXYZ);
}
void bigdft_atoms_get_ntypes_arrays(BigDFT_Atoms *atoms)
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
  // GET_ATTR_DBL_2D(atoms, ATOMS, ig_nlccpar, IG_NLCCPAR);
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
  FC_FUNC_(atoms_new, ATOMS_NEW)(&atoms->data);
  FC_FUNC_(atoms_get, ATOMS_GET)(F_TYPE(atoms->data), &atoms->astruct, &atoms->sym);
  bigdft_atoms_copy_from_fortran(atoms);

  return atoms;
}
BigDFT_Atoms* bigdft_atoms_new_from_fortran(_atoms_data_pointer at)
{
  BigDFT_Atoms *atoms;

#ifdef HAVE_GLIB
  atoms = BIGDFT_ATOMS(g_object_new(BIGDFT_ATOMS_TYPE, NULL));
#else
  atoms = g_malloc(sizeof(BigDFT_Atoms));
  bigdft_atoms_init(atoms);
#endif
  atoms->data = at;
  FC_FUNC_(atoms_get, ATOMS_get)(F_TYPE(atoms->data), &atoms->astruct, &atoms->sym);
  bigdft_atoms_copy_from_fortran(atoms);

  return atoms;
}
void bigdft_atoms_unref(BigDFT_Atoms *atoms)
{
  g_object_unref(G_OBJECT(atoms));
#ifdef HAVE_GLIB
#else
  if (G_OBJECT(atoms)->ref_count <= 0)
    {
      bigdft_atoms_dispose(G_OBJECT(atoms));
      bigdft_atoms_finalize(G_OBJECT(atoms));
      g_free(atoms);
    }
#endif
}
BigDFT_Atoms* bigdft_atoms_new_from_file(const gchar *filename)
{
  BigDFT_Atoms *atoms;

#ifdef HAVE_GLIB
  atoms = BIGDFT_ATOMS(g_object_new(BIGDFT_ATOMS_TYPE, NULL));
#else
  atoms = g_malloc(sizeof(BigDFT_Atoms));
  bigdft_atoms_init(atoms);
#endif
  FC_FUNC_(atoms_new, ATOMS_NEW)(&atoms->data);
  FC_FUNC_(atoms_get, ATOMS_GET)(F_TYPE(atoms->data), &atoms->astruct, &atoms->sym);
  
  if (!bigdft_atoms_set_structure_from_file(atoms, filename))
    {
      bigdft_atoms_unref(atoms);
      atoms = (BigDFT_Atoms*)0;
    }

  return atoms;
}

void bigdft_atoms_set_n_atoms(BigDFT_Atoms *atoms, guint nat)
{
  const gchar subname[] = "bigdft_atoms_set_n_atoms";

  FC_FUNC_(astruct_set_n_atoms, ASTRUCT_SET_N_ATOMS)(F_TYPE(atoms->astruct), (int*)(&nat),
                                                     subname, strlen(subname));
  atoms->nat = nat;
  bigdft_atoms_get_nat_arrays(atoms);
}
static void _sync_atomnames(BigDFT_Atoms *atoms)
{
  gchar name[20];
  guint i, j;

  for (i = 0; i < atoms->ntypes; i++)
    {
      j = i + 1;
      memset(name, ' ', 20);
      if (atoms->atomnames[i])
        memcpy(name, atoms->atomnames[i], strlen(atoms->atomnames[i]));
      FC_FUNC_(atoms_set_name, ATOMS_SET_NAME)(F_TYPE(atoms->data), (int*)(&j), name, 20);
    }
}
/**
 * bigdft_atoms_set_types:
 * @atoms: 
 * @names: (array zero-terminated=1):
 *
 * Pouet.
 **/
void bigdft_atoms_set_types(BigDFT_Atoms *atoms, const gchar **names)
{
  guint i, ntypes;
  const gchar subname[] = "bigdft_atoms_set_types";

  for (ntypes = 0; names[ntypes]; ntypes++);
  FC_FUNC_(astruct_set_n_types, ASTRUCT_SET_N_TYPES)(F_TYPE(atoms->astruct),
                                                     (int*)(&ntypes),
                                                     subname, strlen(subname));
  atoms->ntypes = ntypes;
  bigdft_atoms_get_ntypes_arrays(atoms);
  atoms->atomnames = g_malloc(sizeof(gchar*) * (ntypes + 1));
  for (i = 0; i < ntypes; i++)
    atoms->atomnames[i] = g_strdup(names[i]);
  atoms->atomnames[ntypes] = (gchar*)0;
  _sync_atomnames(atoms);
}
static void _sync_geometry(BigDFT_Atoms *atoms)
{
  gchar format[5], units[20];

  memset(format, ' ', 5);
  if (atoms->format)
    memcpy(format, atoms->format, strlen(atoms->format));
  memset(units, ' ', 20);
  if (atoms->units)
    memcpy(units, atoms->units, strlen(atoms->units));
  FC_FUNC_(astruct_set_geometry, ASTRUCT_SET_GEOMETRY)(F_TYPE(atoms->astruct), atoms->alat,
                                                       &atoms->geocode, format, units, 1, 5, 20);
}
/**
 * bigdft_atoms_set_geometry:
 * @atoms: 
 * @geocode: 
 * @alat: (array fixed-size=3):
 * @units: 
 *
 * 
 **/
void bigdft_atoms_set_geometry(BigDFT_Atoms *atoms, gchar geocode, double alat[3], const gchar *units)
{
  atoms->geocode = geocode;
  atoms->alat[0] = alat[0];
  atoms->alat[1] = alat[1];
  atoms->alat[2] = alat[2];
  strncpy(atoms->units, units, 20);
  _sync_geometry(atoms);
}
void bigdft_atoms_set_default_file_format(BigDFT_Atoms *atoms, const gchar *format)
{
  strncpy(atoms->format, format, 6);
  _sync_geometry(atoms);
}
void bigdft_atoms_copy_from_fortran(BigDFT_Atoms *atoms)
{
  guint ln, i, j;
  gchar str[20];
  int nat, ntypes;

  FC_FUNC_(astruct_copy_geometry_data, ASTRUCT_COPY_GEOMETRY_DATA)
    (F_TYPE(atoms->astruct), &atoms->geocode, atoms->format, atoms->units, 1, 5, 20);
  FC_FUNC_(astruct_copy_alat, ASTRUCT_COPY_ALAT)(F_TYPE(atoms->astruct), atoms->alat);
  FC_FUNC_(astruct_copy_nat, ASTRUCT_COPY_NAT)(F_TYPE(atoms->astruct), &nat);
  if (nat > 0)
    {
      atoms->nat = (guint)nat;
      bigdft_atoms_get_nat_arrays(atoms);
    }
  else
    atoms->nat = 0;
  FC_FUNC_(astruct_copy_ntypes, ASTRUCT_COPY_NTYPES)(F_TYPE(atoms->astruct), &ntypes);
  if (atoms->ntypes > 0)
    {
      atoms->ntypes = (guint)ntypes;
      bigdft_atoms_get_ntypes_arrays(atoms);
    }
  else
    atoms->ntypes = 0;
  atoms->atomnames = g_malloc(sizeof(gchar*) * (atoms->ntypes + 1));
  for (i = 0; i < atoms->ntypes; i++)
    {
      j = i + 1;
      FC_FUNC_(astruct_copy_name, ASTRUCT_COPY_NAME)(F_TYPE(atoms->astruct),
                                                     (int*)(&j), str, (int*)(&ln), 20);
      fprintf(stderr, "%d\n", ln);
      atoms->atomnames[i] = g_malloc(sizeof(gchar) * (ln + 1));
      memcpy(atoms->atomnames[i], str, sizeof(gchar) * ln);
      atoms->atomnames[i][ln] = '\0';
    }
  atoms->atomnames[atoms->ntypes] = (gchar*)0;
}
gboolean bigdft_atoms_set_structure_from_file(BigDFT_Atoms *atoms, const gchar *filename)
{
  guint lstat;

  FC_FUNC_(atoms_empty, ATOMS_EMPTY)(F_TYPE(atoms->data));
  bigdft_atoms_free_additional(atoms);

  FC_FUNC_(astruct_set_from_file, ASTRUCT_SET_FROM_FILE)((int*)(&lstat),
                                                         F_TYPE(atoms->astruct),
                                                         filename, strlen(filename));

  if (!lstat)
    return FALSE;
  else
    bigdft_atoms_copy_from_fortran(atoms);

  return TRUE;
}
/**
 * bigdft_atoms_set_symmetries:
 * @atoms:
 * @active:
 * @tol:
 * @elecfield: (array fixed-size=3):
 *
 */
void bigdft_atoms_set_symmetries(BigDFT_Atoms *atoms, gboolean active,
                                 double tol, double elecfield[3], guint nspin)
{
  int disable;

  disable = (!active);
  FC_FUNC_(astruct_set_symmetries, ASTRUCT_SET_SYMMETRIES)(F_TYPE(atoms->astruct), &disable, &tol, elecfield, (int*)&nspin);
}

void bigdft_atoms_set_displacement(BigDFT_Atoms *atoms, double randdis)
{
  FC_FUNC_(astruct_set_displacement, ASTRUCT_SET_DISPLACEMENT)(F_TYPE(atoms->astruct), &randdis);
}
/**
 * bigdft_atoms_set_psp:
 * @atoms:
 * @ixc:
 * @nspin:
 * @occup: (allow-none):
 *
 */
void bigdft_atoms_set_psp(BigDFT_Atoms *atoms, int ixc, guint nspin, const gchar *occup)
{
  /* int verb = 0; */
  /* int ln; */
  const gchar subname[] = "bigdft_atoms_set_psp";

  /* Allocate the atomic arrays if not already. */
  if (!atoms->psppar)
    {
      FC_FUNC_(allocate_atoms_nat, ALLOCATE_ATOMS_NAT)(F_TYPE(atoms->data), subname, strlen(subname));
      bigdft_atoms_get_nat_arrays(atoms);
      FC_FUNC_(allocate_atoms_ntypes, ALLOCATE_ATOMS_NTYPES)(F_TYPE(atoms->data), subname, strlen(subname));
      bigdft_atoms_get_ntypes_arrays(atoms);
    }

  /* FC_FUNC_(init_atomic_values, INIT_ATOMIC_VALUES)(&verb, atoms->data, &ixc); */
  /* FC_FUNC_(atoms_copy_psp_data, ATOMS_COPY_PSP_DATA) */
  /*   (F_TYPE(atoms->data), (int*)(&atoms->natsc), (int*)(&atoms->donlcc)); */
  /* ln = (occup)?strlen(occup):0; */
  /* FC_FUNC_(atoms_read_variables, ATOMS_READ_VARIABLES)(F_TYPE(atoms->data), (int*)&nspin, */
  /*                                                      occup, &ln, ln); */
}
/**
 * bigdft_atoms_get_radii:
 * @atoms:
 * @crmult:
 * @frmult:
 * @projrad:
 *
 * Returns: (transfer full) (element-type double):
 */
GArray* bigdft_atoms_get_radii(const BigDFT_Atoms *atoms, double crmult,
                               double frmult, double projrad)
{
  GArray *arr;
  double crmult_, frmult_, projrad_;

  arr = g_array_sized_new(FALSE, FALSE, sizeof(double), 3 * atoms->ntypes);
  g_array_set_size(arr, 3 * atoms->ntypes);
  crmult_  = (crmult <= 0.)?5.:crmult;
  frmult_  = (frmult <= 0.)?8.:frmult;
  projrad_ = (projrad <= 0.)?15.:projrad;
  FC_FUNC_(read_radii_variables, READ_RADII_VARIABLES)(F_TYPE(atoms->data), (double*)arr->data,
                                                       &crmult_, &frmult_, &projrad_);
  return arr;
}

/**
 * bigdft_atoms_write:
 * @atoms: 
 * @filename: 
 * @format: (allow-none):
 *
 * 
 **/
void bigdft_atoms_write(BigDFT_Atoms *atoms, const gchar *filename, const gchar *format)
{
  gchar *comment;
  f90_pointer_double_2D forces;

  /* Update format if necessary. */
  if (format)
    {
      strncpy(atoms->format, format, 6);
      _sync_geometry(atoms);
    }
  if (atoms->comment)
    comment = atoms->comment;
  else
    comment = " ";

  F90_2D_POINTER_INIT(&forces);
  FC_FUNC_(atoms_write, ATOMS_WRITE)(F_TYPE(atoms->data), filename, &forces, &atoms->energy, comment,
                                     strlen(filename), strlen(comment));
}

gchar* bigdft_atoms_get_extra_as_label(const BigDFT_Atoms *atoms, guint iat)
{
  gchar buf[50], *ret;
  guint i, j;

  FC_FUNC_(write_extra_info, WRITE_EXTRA_INFO)(buf, (int*)atoms->natpol + iat,
                                               atoms->ifrztyp + iat, 50);
  for (i = 50; i > 0 && buf[i - 1] == ' '; i--);
  if (i == 0)
    return (gchar*)0;

  for (j = 0; buf[j] == ' '; j++);
  ret = g_malloc(sizeof(gchar) * (i - j + 1));
  memcpy(ret, buf + j, sizeof(gchar) * (i - j));
  ret[i - j] = '\0';

  return ret;
}
