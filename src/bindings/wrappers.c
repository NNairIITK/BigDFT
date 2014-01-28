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

/******************************/
/* BigDFT_Dict data structure */
/******************************/
#ifdef GLIB_MAJOR_VERSION
G_DEFINE_TYPE(BigDFT_Dict, bigdft_dict, G_TYPE_OBJECT)

static void bigdft_dict_dispose(GObject *dict);
static void bigdft_dict_finalize(GObject *dict);

static void bigdft_dict_class_init(BigDFT_DictClass *klass)
{
  /* Connect the overloading methods. */
  G_OBJECT_CLASS(klass)->dispose      = bigdft_dict_dispose;
  G_OBJECT_CLASS(klass)->finalize     = bigdft_dict_finalize;
  /* G_OBJECT_CLASS(klass)->set_property = visu_data_set_property; */
  /* G_OBJECT_CLASS(klass)->get_property = visu_data_get_property; */
}
#endif

static void bigdft_dict_init(BigDFT_Dict *obj)
{
#ifdef HAVE_GLIB
  memset((void*)((char*)obj + sizeof(GObject)), 0, sizeof(BigDFT_Dict) - sizeof(GObject));
#else
  memset(obj, 0, sizeof(BigDFT_Dict));
  G_OBJECT(obj)->ref_count = 1;
#endif
}
#ifdef HAVE_GLIB
static void bigdft_dict_dispose(GObject *obj)
{
  BigDFT_Dict *dict = BIGDFT_DICT(obj);

  if (dict->dispose_has_run)
    return;
  dict->dispose_has_run = TRUE;

  /* Chain up to the parent class */
  G_OBJECT_CLASS(bigdft_dict_parent_class)->dispose(obj);
}
#endif
static void bigdft_dict_finalize(GObject *obj)
{
  BigDFT_Dict *dict = BIGDFT_DICT(obj);

  if (dict->root)
    FC_FUNC_(dict_free, DICT_FREE)(&dict->root);

#ifdef HAVE_GLIB
  G_OBJECT_CLASS(bigdft_dict_parent_class)->finalize(obj);
#endif
}
void bigdft_dict_unref(BigDFT_Dict *dict)
{
  g_object_unref(G_OBJECT(dict));
#ifdef HAVE_GLIB
#else
  if (G_OBJECT(dict)->ref_count <= 0)
    {
      bigdft_dict_finalize(G_OBJECT(dict));
      g_free(dict);
    }
#endif
}

/**
 * bigdft_dict_new:
 * @root: (allow-none) (out) (caller-allocates):
 *
 * Pouet.
 *
 * Returns: (transfer full):
 **/
BigDFT_Dict *bigdft_dict_new(BigDFT_DictIter *root)
{
  BigDFT_Dict *dict;

#ifdef HAVE_GLIB
  dict = BIGDFT_DICT(g_object_new(BIGDFT_DICT_TYPE, NULL));
#else
  dict = g_malloc(sizeof(BigDFT_Dict));
  bigdft_dict_init(dict);
#endif
  FC_FUNC_(dict_new, DICT_NEW)(&dict->root);
  dict->current = dict->root;

  if (root)
    {
      (*root).dict = dict;
      (*root).pointer = dict->root;
    } 

  return dict;
}
/**
 * bigdft_dict_new_from_yaml:
 * @buf: 
 *
 * Pouet.
 *
 * Returns: (transfer full):
 **/
BigDFT_Dict *bigdft_dict_new_from_yaml(const gchar *buf)
{
  BigDFT_Dict *dict;

  dict = bigdft_dict_new(NULL);
  FC_FUNC_(dict_parse, DICT_PARSE)(&dict->root, buf, strlen(buf));
  dict->current = dict->root;
  return dict;
}
gboolean bigdft_dict_move_to(BigDFT_Dict *dict, BigDFT_DictIter *iter)
{
  if (iter->dict != dict)
    return FALSE;
  dict->current = iter->pointer;
  return TRUE;
}
/**
 * bigdft_dict_insert:
 * @dict: 
 * @id: 
 * @iter: (out) (caller-allocates) (allow-none):
 *
 * Pouet.
 **/
void bigdft_dict_insert(BigDFT_Dict *dict, const gchar *id, BigDFT_DictIter *iter)
{
  FC_FUNC_(dict_insert, DICT_INSERT)(&dict->current, id, strlen(id));
  if (iter)
    {
      (*iter).dict = dict;
      (*iter).pointer = dict->current;
    }
}
/**
 * bigdft_dict_append:
 * @dict: 
 * @iter: (out) (caller-allocates) (allow-none):
 *
 * Pouet.
 **/
void bigdft_dict_append(BigDFT_Dict *dict, BigDFT_DictIter *iter)
{
  FC_FUNC_(dict_append, DICT_APPEND)(&dict->current);
  if (iter)
    {
      (*iter).dict = dict;
      (*iter).pointer = dict->current;
    }
}
/**
 * bigdft_dict_set:
 * @dict: 
 * @id: (allow-none): 
 * @value: 
 *
 * Pouet.
 **/
void  bigdft_dict_set(BigDFT_Dict *dict, const gchar *id, const gchar *value)
{
  _dictionary *root;
  
  root = dict->current;
  if (id)
    FC_FUNC_(dict_insert, DICT_INSERT)(&dict->current, id, strlen(id));
  FC_FUNC_(dict_put, DICT_PUT)(&dict->current, value, strlen(value));
  dict->current = root;
}
/**
 * bigdft_dict_set_array:
 * @in: 
 * @id: (allow-none):
 * @value: (array zero-terminated=1):
 *
 * 
 **/
void  bigdft_dict_set_array(BigDFT_Dict *dict, const gchar *id, const gchar **value)
{
  guint i;
  _dictionary *root, *key;

  root = dict->current;
  if (id)
    FC_FUNC_(dict_insert, DICT_INSERT)(&dict->current, id, strlen(id));
  key = dict->current;
  for (i = 0; value[i]; i++)
    {
      dict->current = key;
      FC_FUNC_(dict_append, DICT_APPEND)(&dict->current);
      FC_FUNC_(dict_put, DICT_PUT)(&dict->current, value[i], strlen(value[i]));
    }
  dict->current = root;
}
void bigdft_dict_dump(BigDFT_Dict *dict)
{
  FC_FUNC_(dict_dump, DICT_DUMP)(&dict->root);
}
/*********************************/
