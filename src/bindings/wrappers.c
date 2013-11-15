#include <config.h>

#include "bigdft.h"
#include "bindings.h"
#include "bindings_api.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

 /* Duplicate functions in C due to multiple interface definition in Fortran */
void FC_FUNC_(inquire_pointer1, INQUIRE_POINTER1)(void *pt, void *add, int *size)
{
  memcpy(pt, add, sizeof(void*) * *size);
}
void FC_FUNC_(inquire_pointer2, INQUIRE_POINTER2)(void *pt, void *add, int *size)
{
  memcpy(pt, add, sizeof(void*) * *size);
}
void FC_FUNC_(inquire_pointer3, INQUIRE_POINTER3)(void *pt, void *add, int *size)
{
  memcpy(pt, add, sizeof(void*) * *size);
}
void FC_FUNC_(inquire_pointer4, INQUIRE_POINTER4)(void *pt, void *add, int *size)
{
  memcpy(pt, add, sizeof(void*) * *size);
}
void FC_FUNC_(inquire_pointer5, INQUIRE_POINTER5)(void *pt, void *add, int *size)
{
  memcpy(pt, add, sizeof(void*) * *size);
}

void FC_FUNC_(inquire_address1, INQUIRE_ADDRESS1)(double *add, void *pt)
{
  double *val = (double*)(&pt);
  *add = val[0];
}
void FC_FUNC_(inquire_address2, INQUIRE_ADDRESS2)(double *add, void *pt)
{
  double *val = (double*)(&pt);
  *add = val[0];
}

/**
 * bigdft_init_f_lib:
 *
 */
void bigdft_init_f_lib()
{
  FC_FUNC_(f_lib_initialize, F_LIB_INITIALIZE)();
}

/**
 * bigdft_init:
 * @mpi_iproc: (out):
 * @mpi_nproc: (out):
 * @mpi_igroup: (out):
 * @mpi_ngroup: (out):
 * @mpi_groupsize: 0 to use all MPI resources.
 *
 * Setup MPI and other variables. Implicitly call bigdft_mpi_set_distribution().
 *
 * Returns: 
 **/
int bigdft_init(guint *mpi_iproc, guint *mpi_nproc, guint *mpi_igroup, guint *mpi_ngroup,
                guint mpi_groupsize)
{
  int ierr;

  FC_FUNC_(f_lib_initialize, F_LIB_INITIALIZE)();
  FC_FUNC_(bigdft_mpi_init, BIGDFT_MPI_INIT)(&ierr);
  ierr = bigdft_mpi_set_distribution(mpi_iproc, mpi_nproc, mpi_igroup, mpi_ngroup, mpi_groupsize);

#ifdef HAVE_GLIB
  g_type_init();
#endif

  return ierr;
}
/**
 * bigdft_mpi_set_distribution:
 * @mpi_iproc: (out):
 * @mpi_nproc: (out):
 * @mpi_igroup: (out):
 * @mpi_ngroup: (out):
 * @mpi_groupsize: 0 to use all MPI resources.
 *
 * 
 *
 * Returns: 
 **/
int bigdft_mpi_set_distribution(guint *mpi_iproc, guint *mpi_nproc,
                                guint *mpi_igroup, guint *mpi_ngroup,
                                guint mpi_groupsize)
{
  int ierr;
  int info[4];

  FC_FUNC_(bigdft_init_mpi_env, BIGDFT_INIT_MPI_ENV)(info, (int*)&mpi_groupsize, &ierr);
  if (mpi_iproc)
    *mpi_iproc = (guint)info[0];
  if (mpi_nproc)
    *mpi_nproc = (guint)info[1];
  if (mpi_igroup)
    *mpi_igroup = (guint)info[2];
  if (mpi_ngroup)
    *mpi_ngroup = (guint)info[3];

  return ierr;
}
void bigdft_mpi_force_group(guint igroup, guint ngroup)
{
  FC_FUNC_(bigdft_init_mpi_force, BIGDFT_INIT_MPI_FORCE)((int*)&igroup, (int*)&ngroup);
}
int bigdft_finalize()
{
  int ierr;

  FC_FUNC_(bigdft_finalize, BIGDFT_FINALIZE)(&ierr);
  FC_FUNC_(f_lib_finalize, F_LIB_FINALIZE)();
  return ierr;
}

guint bigdft_get_count(GObject *obj)
{
  return G_OBJECT(obj)->ref_count;
}

gchar* _get_c_string(const gchar *fstr, guint len)
{
  guint i;
  
  for (i = len; i > 0 && fstr[i - 1] == ' '; i--);
  if (i == 0)
    return (gchar*)0;

  return g_strndup(fstr, i);
}

#ifndef GLIB_MAJOR_VERSION
GArray* g_array_sized_new(gboolean zero, gboolean nullify, guint ele_size, guint n_ele)
{
  GArray *arr;

  arr = g_malloc(sizeof(GArray));
  arr->data = g_malloc(ele_size * n_ele);
  arr->len = 0;
  return arr;
}
#endif
