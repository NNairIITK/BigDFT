#include "bigdft.h"

#include <config.h>

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

void FC_FUNC_(deallocate_double, DEALLOCATE_DOUBLE)(f90_pointer_double *array);

void FC_FUNC_(read_wave_to_isf, READ_WAVE_TO_ISF)
     (int *lstat, const char* filename, int *ln, int *iorbp,
      double *hx, double *hy, double *hz,
      int *n1, int *n2, int *n3, int *nspinor, f90_pointer_double *psiscf);
void FC_FUNC_(free_wave_to_isf, FREE_WAVE_TO_ISF)(f90_pointer_double *psiscf);

void FC_FUNC_(read_wave_descr, READ_WAVE_DESCR)
     (int *lstat, const char* filename, int *ln, int *norbu,
      int *norbd, int *iorb, int *ispin, int *nkpt, int *ikpt, int *nspinor, int *ispinor);

f90_pointer_double* bigdft_read_wave_to_isf(const char *filename, int iorbp,
                                            double h[3], int n[3], int *nspinor)
{
  int ln, lstat;
  f90_pointer_double *psiscf;

  psiscf = g_malloc(sizeof(f90_pointer_double));
  memset(psiscf, 0, sizeof(f90_pointer_double));

  ln = strlen(filename);
  FC_FUNC_(read_wave_to_isf, READ_WAVE_TO_ISF)
    (&lstat, filename, &ln, &iorbp, h, h + 1, h + 2, n, n + 1, n + 2, nspinor, psiscf);
  if (!lstat)
    {
      g_free(psiscf);
      psiscf = (f90_pointer_double*)0;
    }

  return psiscf;
}
void bigdft_free_wave_to_isf(f90_pointer_double *psiscf)
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

/* The atoms_data structure binding. */
struct f90_pointer_atoms_
{
  void *atoms;
  /* void *info[F90_POINTER_SIZE]; */
};

void FC_FUNC_(atoms_new, ATOMS_NEW)(f90_pointer_atoms *atoms);
void FC_FUNC_(atoms_free, ATOMS_FREE)(f90_pointer_atoms *atoms);
void FC_FUNC_(atoms_set_n_atoms, ATOMS_SET_N_ATOMS)(f90_pointer_atoms *atoms, int *nat);
void FC_FUNC_(atoms_set_n_types, ATOMS_SET_N_TYPES)(f90_pointer_atoms *atoms, int *ntypes);
void FC_FUNC_(atoms_get_geocode, ATOMS_GET_GEOCODE)(f90_pointer_atoms *atoms, char *geocode);
void FC_FUNC_(atoms_get_nat, ATOMS_GET_NAT)(f90_pointer_atoms *atoms, int *nat);
void FC_FUNC_(atoms_get_ntypes, ATOMS_GET_NTYPES)(f90_pointer_atoms *atoms, int *ntypes);
void FC_FUNC_(atoms_get_iatype, ATOMS_GET_IATYPE)(f90_pointer_atoms *atoms,
                                                  f90_pointer_int *iatype);
void FC_FUNC_(atoms_new_from_file, ATOMS_NEW_FROM_FILE)(int *lstat, f90_pointer_atoms *atoms,
                                                        f90_pointer_double *rxyz,
                                                        const gchar *filename, int *ln);
void FC_FUNC_(init_atomic_values, INIT_ATOMIC_VALUES)(int *iproc, f90_pointer_atoms *atoms,
                                                      int *ixc);
void FC_FUNC_(read_radii_variables, READ_RADII_VARIABLES)(f90_pointer_atoms *atoms,
                                                          double *radii_cf);

#define GET_ATTR_INT(name,NAME) { \
  f90_pointer_int tmp; \
  memset(&tmp, 0, sizeof(f90_pointer_int)); \
  FC_FUNC_(atoms_get_ ## name, ATOMS_GET_ ## NAME)(atoms->data->atoms, &tmp); \
  atoms->name = (guint*)tmp.data; \
  }

BigDFT_Atoms* bigdft_atoms_init()
{
  BigDFT_Atoms *atoms;

  atoms = g_malloc(sizeof(BigDFT_Atoms));
  memset(atoms, 0, sizeof(BigDFT_Atoms));
  atoms->data = g_malloc(sizeof(f90_pointer_atoms));
  memset(atoms->data, 0, sizeof(f90_pointer_atoms));
  memset(&atoms->rxyz, 0, sizeof(f90_pointer_double));

  return atoms;
}
void bigdft_atoms_dispose(BigDFT_Atoms *atoms)
{
  g_free(atoms->data);
  FC_FUNC_(deallocate_double, DEALLOCATE_DOUBLE)(&atoms->rxyz);
  g_free(atoms);
}

BigDFT_Atoms* bigdft_atoms_new()
{
  BigDFT_Atoms *atoms;

  atoms = bigdft_atoms_init();
  FC_FUNC_(atoms_new, ATOMS_NEW)(atoms->data);

  return atoms;
}
void bigdft_atoms_free(BigDFT_Atoms *atoms)
{
  FC_FUNC_(atoms_free, ATOMS_FREE)(atoms->data);
  bigdft_atoms_dispose(atoms);
}
BigDFT_Atoms* bigdft_atoms_new_from_file(const gchar *filename)
{
  BigDFT_Atoms *atoms;
  guint lstat, ln;

  atoms = bigdft_atoms_init();
  ln = strlen(filename);
  FC_FUNC_(atoms_new_from_file, ATOMS_NEW_FROM_FILE)((int*)(&lstat), atoms->data,
                                                     &atoms->rxyz, filename, (int*)(&ln));

  if (!lstat)
    {
      bigdft_atoms_dispose(atoms);
      atoms = (BigDFT_Atoms*)0;
    }
  else
    {
      FC_FUNC_(atoms_get_geocode, ATOMS_GET_GEOCODE)(atoms->data->atoms, (gchar*)(&atoms->geocode));
      FC_FUNC_(atoms_get_nat, ATOMS_GET_NAT)(atoms->data->atoms, (int*)(&atoms->nat));
      FC_FUNC_(atoms_get_ntypes, ATOMS_GET_NTYPES)(atoms->data->atoms, (int*)(&atoms->ntypes));
      GET_ATTR_INT(iatype,IATYPE);
    }

  return atoms;
}
void bigdft_atoms_set_n_atoms(BigDFT_Atoms *atoms, guint nat)
{
  FC_FUNC_(atoms_set_n_atoms, ATOMS_SET_N_ATOMS)(atoms->data->atoms, (int*)(&nat));
  atoms->nat = nat;

  GET_ATTR_INT(iatype,IATYPE);
}
void bigdft_atoms_set_n_types(BigDFT_Atoms *atoms, guint ntypes)
{
  FC_FUNC_(atoms_set_n_types, ATOMS_SET_N_TYPES)(atoms->data, (int*)(&ntypes));
  atoms->ntypes = ntypes;
}

void bigdft_atoms_set_psp(BigDFT_Atoms *atoms, int ixc)
{
  int iproc = 0;
  FC_FUNC_(init_atomic_values, INIT_ATOMIC_VALUES)(&iproc, atoms->data->atoms, &ixc);
}

double* bigdft_atoms_get_radii(BigDFT_Atoms *atoms)
{
  double *radii_cf;

  radii_cf = g_malloc(sizeof(double) * 3 * atoms->ntypes);
  FC_FUNC_(read_radii_variables, READ_RADII_VARIABLES)(atoms->data->atoms, radii_cf);
  return radii_cf;
}


/* Wavefunction descriptor part. */
struct f90_pointer_glr_
{
  void *glr;
  /* void *info[F90_POINTER_SIZE]; */
};

void FC_FUNC_(glr_new, GLR_NEW)(f90_pointer_glr *glr);
void FC_FUNC_(system_size, SYSTEM_SIZE)(int *iproc, f90_pointer_atoms *atoms, double *rxyz,
                                        double *radii_cf, double *crmult, double *frmult,
                                        double *hx, double *hy, double *hz,
                                        f90_pointer_glr *glr, double *shift);
void FC_FUNC_(glr_get_n, GLR_GET_N)(f90_pointer_glr *glr, int *n);
void FC_FUNC_(glr_free, GLR_FREE)(f90_pointer_glr *glr);

BigDFT_Glr* bigdft_glr_init()
{
  BigDFT_Glr *glr;

  glr = g_malloc(sizeof(BigDFT_Glr));
  memset(glr, 0, sizeof(BigDFT_Glr));
  glr->data = g_malloc(sizeof(f90_pointer_glr));
  memset(glr->data, 0, sizeof(f90_pointer_glr));

  return glr;
}
void bigdft_glr_dispose(BigDFT_Glr *glr)
{
  g_free(glr->data);
  g_free(glr);
}

BigDFT_Glr* bigdft_glr_new(BigDFT_Atoms *atoms, double *radii, double h[3],
                           double crmult, double frmult)
{
  BigDFT_Glr *glr;
  int iproc = 0;

  glr = bigdft_glr_init();
  FC_FUNC_(glr_new, GLR_NEW)(glr->data);
  FC_FUNC_(system_size, SYSTEM_SIZE)(&iproc, atoms->data->atoms, atoms->rxyz.data, radii,
                                     &crmult, &frmult, h, h + 1, h + 2, glr->data->glr,
                                     atoms->shift);
  glr->h[0] = h[0];
  glr->h[1] = h[1];
  glr->h[2] = h[2];
  FC_FUNC_(glr_get_n, GLR_GET_N)(glr->data->glr, glr->n);
  
  return glr;
}
void bigdft_glr_free(BigDFT_Glr *glr)
{
  FC_FUNC_(glr_free, GLR_FREE)(glr->data);
  bigdft_glr_dispose(glr);
}

void FC_FUNC_(fill_logrid, FILL_LOGRID)(char *geocode, int *n1, int *n2, int *n3,
                                        int *nl1, int *nu1, int *nl2, int *nu2, int *nl3, int *nu3,
                                        int *orig, int *nat, int *ntypes, int *iatype,
                                        double *rxyz, double *radii, double *mult,
                                        double *hx, double *hy, double *hz, int *grid);
guint* bigdft_fill_logrid(BigDFT_Atoms *atoms, guint n[3], double *radii,
                          double mult, double h[3])
{
  guint *grid;
  int orig = 0;

  grid = g_malloc(sizeof(guint) * (n[0] + 1) * (n[1] + 1) * (n[2] + 1));

  FC_FUNC_(fill_logrid, FILL_LOGRID)(&atoms->geocode, (int*)n, (int*)(n + 1), (int*)(n + 2),
                                     &orig, (int*)n, &orig, (int*)(n + 1), &orig, (int*)(n + 2),
                                     &orig, (int*)(&atoms->nat), (int*)(&atoms->ntypes), (int*)atoms->iatype,
                                     atoms->rxyz.data, radii, &mult, h, h + 1, h + 2, (int*)grid);

  return grid;
}
