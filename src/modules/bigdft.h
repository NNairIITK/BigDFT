#ifndef BIGDFT_H
#define BIGDFT_H

#ifndef GLIB_MAJOR_VERSION
#define TRUE 1
#define FALSE 0
#define gboolean int
#define g_malloc(A) malloc(A)
#define g_free(A)   free(A)
#endif

#define F90_POINTER_SIZE 16

typedef struct f90_pointer_double_
{
  double *data;
  void *info[F90_POINTER_SIZE];
} f90_pointer_double;

f90_pointer_double* bigdft_read_wave_to_isf(const char *filename, int iorbp,
                                            double h[3], int n[3], int *nspinor);
void bigdft_free_wave_to_isf(f90_pointer_double *psiscf);
gboolean bigdft_read_wave_descr(const char *filename, int *norbu,
                                int *norbd, int *nkpt, int *nspinor,
                                int *iorb, int *ispin, int *ikpt, int *ispinor);

#endif
