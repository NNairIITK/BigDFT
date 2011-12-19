#ifndef BIGDFT_H
#define BIGDFT_H

#ifndef GLIB_MAJOR_VERSION
#include <stdlib.h>
#define TRUE 1
#define FALSE 0
#define gboolean int
#define g_malloc(A) malloc(A)
#define g_free(A)   free(A)
#define guint unsigned int
#define gchar char
#endif

#define F90_POINTER_SIZE 16

typedef struct f90_pointer_double_
{
  double *data;
  void *info[F90_POINTER_SIZE];
} f90_pointer_double;
typedef struct f90_pointer_int_
{
  int *data;
  void *info[F90_POINTER_SIZE];
} f90_pointer_int;

f90_pointer_double* bigdft_read_wave_to_isf(const gchar *filename, int iorbp,
                                            double h[3], int n[3], int *nspinor);
void bigdft_free_wave_to_isf(f90_pointer_double *psiscf);
gboolean bigdft_read_wave_descr(const gchar *filename, int *norbu,
                                int *norbd, int *nkpt, int *nspinor,
                                int *iorb, int *ispin, int *ikpt, int *ispinor);

typedef struct f90_pointer_atoms_ f90_pointer_atoms;
typedef struct BigDFT_Atoms_
{
  /* TODO: bindings to values... */
  gchar geocode;
  guint  ntypes, nat;
  guint *iatype;
  

  /* Coordinates. */
  f90_pointer_double rxyz;
  double shift[3];

  /* Private. */
  f90_pointer_atoms *data;
} BigDFT_Atoms;

BigDFT_Atoms* bigdft_atoms_new();
void bigdft_atoms_free(BigDFT_Atoms *atoms);
BigDFT_Atoms* bigdft_atoms_new_from_file(const gchar *filename);
void bigdft_atoms_set_n_atoms(BigDFT_Atoms *atoms, guint nat);
void bigdft_atoms_set_n_types(BigDFT_Atoms *atoms, guint ntypes);
void bigdft_atoms_set_psp(BigDFT_Atoms *atoms, int ixc);
double* bigdft_atoms_get_radii(BigDFT_Atoms *atoms);

typedef struct f90_pointer_glr_ f90_pointer_glr;
typedef struct BigDFT_glr_
{
  double h[3];
  int n[3];
  
  /* TODO: bindings to values... */

  /* Private. */
  f90_pointer_glr *data;
} BigDFT_Glr;

BigDFT_Glr* bigdft_glr_new(BigDFT_Atoms *atoms, double *radii, double h[3],
                           double crmult, double frmult);
void bigdft_glr_free(BigDFT_Glr *glr);

guint* bigdft_fill_logrid(BigDFT_Atoms *atoms, guint n[3], double *radii,
                          double mult, double h[3]);

#endif
