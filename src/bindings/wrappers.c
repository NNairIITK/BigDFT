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
 * bigdft_init:
 * @mpi_iproc: (out):
 * @mpi_nproc: (out):
 * @mpi_igroup: (out):
 * @mpi_ngroup: (out):
 * @mpi_groupsize: 0 to use all MPI resources.
 *
 * Setup MPI and other variables.
 *
 * Returns: 
 **/
int bigdft_init(guint *mpi_iproc, guint *mpi_nproc, guint *mpi_igroup, guint *mpi_ngroup,
                guint mpi_groupsize)
{
  int ierr;
  int info[4];

  FC_FUNC_(bigdft_mpi_init, BIGDFT_MPI_INIT)(&ierr);
  FC_FUNC_(bigdft_init_mpi_env, BIGDFT_INIT_MPI_ENV)(info, (int*)&mpi_groupsize, &ierr);
  if (mpi_iproc)
    *mpi_iproc = (guint)info[0];
  if (mpi_nproc)
    *mpi_nproc = (guint)info[1];
  if (mpi_igroup)
    *mpi_igroup = (guint)info[2];
  if (mpi_ngroup)
    *mpi_ngroup = (guint)info[3];

#ifdef HAVE_GLIB
  g_type_init();
#endif

  return ierr;
}
int bigdft_finalize()
{
  int ierr;

  FC_FUNC_(bigdft_finalize, BIGDFT_FINALIZE)(&ierr);
  return ierr;
}
