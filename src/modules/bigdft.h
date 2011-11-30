#ifndef BIGDFT_H
#define BIGDFT_H

#define F90_POINTER_SIZE 16

typedef struct f90_pointer_double_
{
  double *data;
  void *info[F90_POINTER_SIZE];
} f90_pointer_double;

f90_pointer_double* bigdft_read_wave_to_isf_etsf(const char *filename, int iorbp,
                                                 double h[3], int n[3], int *nspinor);
void bigdft_free_wave_to_isf_etsf(f90_pointer_double *psiscf);

#endif
