#include "bigdft.h"

#include <config.h>

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

void FC_FUNC_(read_wave_to_isf_etsf, READ_WAVE_TO_ISF_ETSF)
     (int *lstat, const char* filename, int *ln, int *iorbp,
      double *hx, double *hy, double *hz,
      int *n1, int *n2, int *n3, int *nspinor, f90_pointer_double *psiscf);
void FC_FUNC_(free_wave_to_isf_etsf, FREE_WAVE_TO_ISF_ETSF)(f90_pointer_double *psiscf);

void FC_FUNC_(read_wave_descr_etsf, READ_WAVE_DESCR_ETSF)
     (int *lstat, const char* filename, int *ln, int *norbu,
      int *norbd, int *nkpt, int *nspinor);

f90_pointer_double* bigdft_read_wave_to_isf_etsf(const char *filename, int iorbp,
                                                 double h[3], int n[3], int *nspinor)
{
  int ln, lstat;
  char *str;
  f90_pointer_double *psiscf;

  psiscf = g_malloc(sizeof(f90_pointer_double));
  memset(psiscf, 0, sizeof(f90_pointer_double));

  ln = strlen(filename);
  str = g_malloc(sizeof(char) * (ln + 1));
  memcpy(str, filename, sizeof(char) * ln);
  str[ln] = ' ';
  FC_FUNC_(read_wave_to_isf_etsf, READ_WAVE_TO_ISF_ETSF)
    (&lstat, str, &ln, &iorbp, h, h + 1, h + 2, n, n + 1, n + 2, nspinor, psiscf);
  g_free(str);
  if (!lstat)
    {
      g_free(psiscf);
      psiscf = (f90_pointer_double*)0;
    }

  return psiscf;
}
void bigdft_free_wave_to_isf_etsf(f90_pointer_double *psiscf)
{
  FC_FUNC_(free_wave_to_isf_etsf, FREE_WAVE_TO_ISF_ETSF)(psiscf);
  g_free(psiscf);
}

gboolean bigdft_read_wave_descr_etsf(const char *filename, int *norbu,
                                     int *norbd, int *nkpt, int *nspinor)
{
  int ln, lstat, norbu_, norbd_, nkpt_, nspinor_;
  char *str;
  
  ln = strlen(filename);
  str = g_malloc(sizeof(char) * (ln + 1));
  memcpy(str, filename, sizeof(char) * ln);
  str[ln] = ' ';
  FC_FUNC_(read_wave_descr_etsf, READ_WAVE_DESCR_ETSF)
    (&lstat, str, &ln, &norbu_, &norbd_, &nkpt_, &nspinor_);
  g_free(str);
  if (!lstat)
    return FALSE;

  if (norbu)   *norbu   = norbu_;
  if (norbd)   *norbd   = norbd_;
  if (nkpt)    *nkpt    = nkpt_;
  if (nspinor) *nspinor = nspinor_;

  return TRUE;
}
