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
 * bigdft_lib_init:
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
int bigdft_lib_init(guint *mpi_iproc, guint *mpi_nproc, guint *mpi_igroup, guint *mpi_ngroup,
                guint mpi_groupsize)
{
  int ierr;

  FC_FUNC_(f_lib_initialize, F_LIB_INITIALIZE)();
  FC_FUNC_(bigdft_mpi_init, BIGDFT_MPI_INIT)(&ierr);
  ierr = bigdft_mpi_set_distribution(mpi_iproc, mpi_nproc, mpi_igroup, mpi_ngroup, mpi_groupsize);

#ifdef HAVE_GLIB
#if GLIB_MINOR_VERSION < 36
  g_type_init();
#endif
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
int bigdft_lib_finalize()
{
  int ierr;

  FC_FUNC_(bigdft_finalize, BIGDFT_FINALIZE)(&ierr);
  FC_FUNC_(f_lib_finalize, F_LIB_FINALIZE)();
  return ierr;
}
/**
 * bigdft_lib_err_severe_override:
 * @func: (allow-none) (scope call): a routine.
 *
 * Change the default callback for severe errors.
 **/
void bigdft_lib_err_severe_override(BigdftErrorCallback func)
{
  FC_FUNC_(call_external_c_fromadd, CALL_EXTERNAL_C_FROMADD)(&func);
  /*FC_FUNC_(err_severe_override, ERR_SEVERE_OVERRIDE)(func);*/
  fprintf(stderr, "%ld\n", (long long int)func);
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

  /* g_message("New dict %p.", (gpointer)obj); */
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

  /* g_message("Killing %p (%p).", (gpointer)obj, dict->root); */
  if (F_TYPE(dict->root))
    FC_FUNC_(dict_free_binding, DICT_FREE_BINDING)(&dict->root);

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
 * @root: (allow-none) (out caller-allocates):
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
  FC_FUNC_(dict_init_binding, DICT_INIT_BINDING)(&dict->root);
  dict->current = dict->root;

  if (root)
    {
      root->dict = dict;
      root->pointer = dict->root;
    }

  return dict;
}
/**
 * bigdft_dict_new_from_yaml:
 * @buf: 
 * @root: (allow-none) (out caller-allocates):
 *
 * Pouet.
 *
 * Returns: (transfer full):
 **/
BigDFT_Dict *bigdft_dict_new_from_yaml(const gchar *buf, BigDFT_DictIter *root)
{
  BigDFT_Dict *dict;

#ifdef HAVE_GLIB
  dict = BIGDFT_DICT(g_object_new(BIGDFT_DICT_TYPE, NULL));
#else
  dict = g_malloc(sizeof(BigDFT_Dict));
  bigdft_dict_init(dict);
#endif
  FC_FUNC_(dict_parse, DICT_PARSE)(&dict->root, buf, strlen(buf));
  dict->current = dict->root;

  if (root)
    {
      root->dict = dict;
      root->pointer = dict->root;
    }

  return dict;
}
BigDFT_Dict *bigdft_dict_new_from_fortran(f90_dictionary_pointer dictf)
{
  BigDFT_Dict *dict;

#ifdef HAVE_GLIB
  dict = BIGDFT_DICT(g_object_new(BIGDFT_DICT_TYPE, NULL));
#else
  dict = g_malloc(sizeof(BigDFT_Dict));
  bigdft_dict_init(dict);
#endif
  dict->root = dictf;
  dict->current = dict->root;

  return dict;
}
/**
 * bigdft_dict_get_current:
 * @dict: 
 * @iter: (out caller-allocates)
 *
 * 
 **/
void bigdft_dict_get_current(BigDFT_Dict *dict, BigDFT_DictIter *iter)
{
  iter->dict = dict;
  iter->pointer = dict->current;
}
/**
 * bigdft_dict_move_to:
 * @dict: 
 * @iter: (allow-none):
 *
 * 
 *
 * Returns: 
 **/
gboolean bigdft_dict_move_to(BigDFT_Dict *dict, BigDFT_DictIter *iter)
{
  if (iter && iter->dict != dict)
    return FALSE;

  dict->current = (iter)?iter->pointer:dict->root;
  return TRUE;
}
/**
 * bigdft_dict_move_to_key:
 * @dict: 
 * @iter: (out caller-allocates) (allow-none):
 * @key:
 *
 * 
 *
 * Returns: TRUE, if @key exists.
 **/
gboolean bigdft_dict_move_to_key(BigDFT_Dict *dict, const gchar *key, BigDFT_DictIter *iter)
{
  int exists;

  FC_FUNC_(dict_move_to_key, DICT_MOVE_TO_KEY)(&dict->current, &exists, key, strlen(key));
  if (iter)
    {
      iter->dict = dict;
      iter->pointer = dict->current;
    }
  return (gboolean)exists;
}
/**
 * bigdft_dict_move_to_item:
 * @dict: 
 * @iter: (out caller-allocates) (allow-none):
 * @id:
 *
 * 
 *
 * Returns: TRUE, if @key exists.
 **/
gboolean bigdft_dict_move_to_item(BigDFT_Dict *dict, guint id, BigDFT_DictIter *iter)
{
  int exists;

  FC_FUNC_(dict_move_to_item, DICT_MOVE_TO_ITEM)(&dict->current, &exists, (int*)&id);
  if (iter)
    {
      iter->dict = dict;
      iter->pointer = dict->current;
    }
  return (gboolean)exists;
}
/**
 * bigdft_dict_iter:
 * @dict: 
 * @iter: (out caller-allocates) (allow-none):
 *
 * 
 *
 * Returns: 
 **/
gboolean bigdft_dict_iter(BigDFT_Dict *dict, BigDFT_DictIter *iter)
{
  int exists;

  FC_FUNC_(dict_iter, DICT_ITER)(&dict->current, &exists);
  if (iter && exists)
    {
      iter->dict = dict;
      iter->pointer = dict->current;
    }  
  return (gboolean)exists;
}
/**
 * bigdft_dict_next:
 * @dict: 
 * @iter: (out caller-allocates) (allow-none):
 *
 * 
 *
 * Returns: 
 **/
gboolean bigdft_dict_next(BigDFT_Dict *dict, BigDFT_DictIter *iter)
{
  int exists;

  FC_FUNC_(dict_next, DICT_NEXT)(&dict->current, &exists);
  if (iter && exists)
    {
      iter->dict = dict;
      iter->pointer = dict->current;
    }  
  return (gboolean)exists;
}
/**
 * bigdft_dict_insert:
 * @dict: 
 * @key: 
 * @iter: (out caller-allocates) (allow-none):
 *
 * Pouet.
 **/
void bigdft_dict_insert(BigDFT_Dict *dict, const gchar *key, BigDFT_DictIter *iter)
{
  FC_FUNC_(dict_insert, DICT_INSERT)(&dict->current, key, strlen(key));
  if (iter)
    {
      iter->dict = dict;
      iter->pointer = dict->current;
    }
}
/**
 * bigdft_dict_append:
 * @dict: 
 * @iter: (out caller-allocates) (allow-none):
 *
 * Pouet.
 **/
void bigdft_dict_append(BigDFT_Dict *dict, BigDFT_DictIter *iter)
{
  FC_FUNC_(dict_append, DICT_APPEND)(&dict->current);
  if (iter)
    {
      iter->dict = dict;
      iter->pointer = dict->current;
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
  f90_dictionary_pointer root;
  
  root = dict->current;
  if (id)
    FC_FUNC_(dict_insert, DICT_INSERT)(&dict->current, id, strlen(id));
  FC_FUNC_(dict_put, DICT_PUT)(&dict->current, value, strlen(value));
  dict->current = root;
}
/**
 * bigdft_dict_set_array:
 * @dict: 
 * @id: (allow-none):
 * @value: (array zero-terminated=1):
 *
 * 
 **/
void  bigdft_dict_set_array(BigDFT_Dict *dict, const gchar *id, const gchar **value)
{
  guint i;
  f90_dictionary_pointer root, key;

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
/**
 * bigdft_dict_set_dict:
 * @dict: 
 * @id: (allow-none):
 * @value:
 *
 * 
 **/
void  bigdft_dict_set_dict(BigDFT_Dict *dict, const gchar *id, const BigDFT_Dict *value)
{
  f90_dictionary_pointer root;
  
  root = dict->current;
  if (id)
    FC_FUNC_(dict_insert, DICT_INSERT)(&dict->current, id, strlen(id));
  FC_FUNC_(dict_update_binding, DICT_UPDATE_BINDING)(&dict->current, &value->current);
  dict->current = root;
}
gboolean bigdft_dict_pop(BigDFT_Dict *dict, const gchar *key)
{
  int exists;

  FC_FUNC_(dict_pop, DICT_POP)(&dict->current, &exists, key, strlen(key));
  return (gboolean)exists;
}
/**
 * bigdft_dict_value:
 * @dict: 
 *
 * 
 *
 * Returns: (transfer full):
 **/
gchar* bigdft_dict_value(BigDFT_Dict *dict)
{
#define max_field_length 256
  char buf[max_field_length + 1];
  guint i, ln;
  gchar *out;
  
  buf[max_field_length] = ' ';
  FC_FUNC_(dict_value_binding, DICT_VALUE_BINDING)(&dict->current, buf, max_field_length);
  for (i = max_field_length; i > 0 && buf[i] == ' '; i--)
    buf[i] = '\0';
  ln = max_field_length - i;
  out = g_malloc(sizeof(gchar) * (ln + 1));
  memcpy(out, buf, sizeof(gchar) * (ln + 1));
  return out;
}
/**
 * bigdft_dict_key:
 * @dict: 
 *
 * 
 *
 * Returns: (transfer full):
 **/
gchar* bigdft_dict_key(BigDFT_Dict *dict)
{
#define max_field_length 256
  char buf[max_field_length + 1];
  guint i, ln;
  gchar *out;
  
  buf[max_field_length] = ' ';
  FC_FUNC_(dict_key, DICT_KEY)(&dict->current, buf, max_field_length);
  for (i = max_field_length; i >= 0 && buf[i] == ' '; i--)
    buf[i] = '\0';
  ln = max_field_length - i;
  out = g_malloc(sizeof(gchar) * (ln + 1));
  memcpy(out, buf, sizeof(gchar) * (ln + 1));
  return out;
}
guint bigdft_dict_len(BigDFT_Dict *dict)
{
  int ln;

  FC_FUNC_(dict_len, DICT_LEN)(&dict->current, &ln);
  if (ln >= 0)
    return (guint)ln;
  FC_FUNC_(dict_size, DICT_SIZE)(&dict->current, &ln);
  if (ln >= 0)
    return (guint)ln;
  return 1;
}
void bigdft_dict_dump(BigDFT_Dict *dict, gint unit)
{
  FC_FUNC_(dict_dump, DICT_DUMP)(&dict->root, &unit);
}
void bigdft_dict_dump_to_file(BigDFT_Dict *dict, const gchar *filename)
{
  FC_FUNC_(dict_dump_to_file, DICT_DUMP_TO_FILE)(&dict->root, filename, strlen(filename));
}
/*********************************/
