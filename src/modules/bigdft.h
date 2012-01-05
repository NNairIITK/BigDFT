/*> @file 
     Header for the public BigDFT API.
    @author
     Copyright (C) 2011-2012 BigDFT group 
     This file is distributed under the terms of the
     GNU General Public License, see ~/COPYING file
     or http://www.gnu.org/copyleft/gpl.txt .
     For the list of contributors, see ~/AUTHORS 
*/

#ifndef BIGDFT_H
#define BIGDFT_H

#include <bigdft_cst.h>

/***************************/
/* Generic pointer arrays. */
/***************************/
typedef struct f90_pointer_double_
{
  void *shift[POINTER_SHIFT_SIZE];
  double *data;
  void *info[F90_POINTER_SIZE];
} f90_pointer_double;
typedef struct f90_pointer_int_
{
  void *shift[POINTER_SHIFT_SIZE];
  int *data;
  void *info[F90_POINTER_SIZE];
} f90_pointer_int;

/********************************/
/* BigDFT_Atoms data structure. */
/********************************/
typedef struct f90_pointer_atoms_ f90_pointer_atoms;
typedef struct BigDFT_Atoms_
{
  /* TODO: bindings to values... */
  gchar geocode;
  guint  ntypes, nat;
  guint *iatype;
  gchar **atomnames;
  double alat[3];

  /* Coordinates. */
  f90_pointer_double rxyz;
  double shift[3];

  /* Private. */
  f90_pointer_atoms *data;
} BigDFT_Atoms;

BigDFT_Atoms* bigdft_atoms_new();
BigDFT_Atoms* bigdft_atoms_new_from_file(const gchar *filename);
void          bigdft_atoms_free       (BigDFT_Atoms *atoms);
void          bigdft_atoms_set_n_atoms(BigDFT_Atoms *atoms, guint nat);
void          bigdft_atoms_set_n_types(BigDFT_Atoms *atoms, guint ntypes);
void          bigdft_atoms_set_psp    (BigDFT_Atoms *atoms, int ixc);
double*       bigdft_atoms_get_radii  (BigDFT_Atoms *atoms);

/******************************/
/* BigDFT_Glr data structure. */
/******************************/
typedef struct f90_pointer_glr_ f90_pointer_glr;
typedef struct BigDFT_glr_
{
  double h[3];
  guint n[3];
  
  /* TODO: bindings to values... */

  /* Private. */
  f90_pointer_glr *data;
} BigDFT_Glr;

BigDFT_Glr* bigdft_glr_new (BigDFT_Atoms *atoms, double *radii, double h[3],
                            double crmult, double frmult);
void        bigdft_glr_free(BigDFT_Glr *glr);

/*********************************/
/* BigDFT_Inputs data structure. */
/*********************************/
typedef struct f90_pointer_inputs_ f90_pointer_inputs;
typedef struct BigDFT_Inputs_
{
  /* TODO: bindings to values... */
  int files;
  int ixc, ncharge, nspin, mpol, itermax, nrepmax, ncong, idsx,
    dispersion, inputPsiId, output_wf_format, output_grid, ncongt, norbv, nvirt,
    nplot, disableSym;
  double crmult, frmult, gnrm_cv, rbuf;
  double h[3], elecfield[3];

  /* Private. */
  f90_pointer_inputs *data;
} BigDFT_Inputs;

#define BIGDFT_INPUTS_NONE    0
#define BIGDFT_INPUTS_DFT     1
#define BIGDFT_INPUTS_GEOPT   2
#define BIGDFT_INPUTS_PERF    4
#define BIGDFT_INPUTS_KPT     8
#define BIGDFT_INPUTS_MIX    16
#define BIGDFT_INPUTS_TDDFT  32
#define BIGDFT_INPUTS_SIC    64
#define BIGDFT_INPUTS_FREQ  128

BigDFT_Inputs* bigdft_inputs_new (const gchar *radical);
void           bigdft_inputs_free(BigDFT_Inputs *in);

/******************/
/* Miscellaneous. */
/******************/
guint* bigdft_fill_logrid(BigDFT_Atoms *atoms, guint n[3], double *radii,
                          double mult, double h[3]);
f90_pointer_double* bigdft_read_wave_to_isf(const gchar *filename, int iorbp,
                                            double h[3], int n[3], int *nspinor);
void bigdft_free_wave_to_isf(f90_pointer_double *psiscf);
gboolean bigdft_read_wave_descr(const gchar *filename, int *norbu,
                                int *norbd, int *nkpt, int *nspinor,
                                int *iorb, int *ispin, int *ikpt, int *ispinor);


#endif
