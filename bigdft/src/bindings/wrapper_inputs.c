/* @file
 * Bindings for the BigDFT package
 * @author
 * Copyright (C) 2013-2013 BigDFT group
 * This file is distributed under the terms of the
 * GNU General Public License, see ~/COPYING file
 * or http://www.gnu.org/copyleft/gpl.txt .
 * For the list of contributors, see ~/AUTHORS
**/


#include <config.h>

#include "bigdft.h"
#include "bindings.h"
#include "bindings_api.h"

#include <string.h>
#include <stdio.h>




static void _free_output(BigDFT_Inputs *in)
{
  g_free(in->dir_output);
}


void _inputs_sync(BigDFT_Inputs *in)

{
  gchar dir_output[100];

  FC_FUNC_(inputs_get_dft, INPUTS_GET_DFT)(F_TYPE(in->data), in->h, in->h + 1, in->h + 2,
                                           &in->crmult, &in->frmult, &in->ixc,
                                           &in->ncharge, in->elecfield, &in->nspin,
                                           &in->mpol, &in->gnrm_cv, (int*)&in->itermax,
                                           (int*)&in->nrepmax, &in->ncong, (int*)&in->idsx,
                                           &in->dispersion, &in->inputPsiId,
                                           &in->output_wf_format, &in->output_grid,
                                           &in->rbuf, &in->ncongt, &in->norbv, &in->nvirt,
                                           &in->nplot, &in->disableSym, &in->last_run);
  FC_FUNC_(inputs_get_mix, INPUTS_GET_MIX)(F_TYPE(in->data), (int*)&in->iscf, (int*)&in->itrpmax,
                                           (int*)&in->norbsempty, (int*)(&in->occopt),
                                           &in->alphamix,
                                           &in->rpnrm_cv, &in->gnrm_startmix, &in->Tel,
                                           &in->alphadiis);
  FC_FUNC_(inputs_get_geopt, INPUTS_GET_GEOPT)(F_TYPE(in->data), in->geopt_approach,
                                               &in->ncount_cluster_x, &in->frac_fluct,
                                               &in->forcemax, &in->randdis, &in->betax,
                                               &in->history, &in->ionmov, &in->dtion,
                                               in->strtarget, &in->qmass, 10);
  /* FC_FUNC_(inputs_get_sic, INPUTS_GET_SIC)(); */
  /* FC_FUNC_(inputs_get_tddft, INPUTS_GET_TDDFT)(); */
  FC_FUNC_(inputs_get_perf, INPUTS_GET_PERF)(F_TYPE(in->data), (int*)&in->linear);

  FC_FUNC_(inputs_get_output, INPUTS_GET_OUTPUT)(F_TYPE(in->data), dir_output, 100);
  in->dir_output = _get_c_string(dir_output, 100);
}


static BigDFT_Inputs* bigdft_inputs_init()
{
  BigDFT_Inputs *in;

  in = g_malloc(sizeof(BigDFT_Inputs));
  memset(in, 0, sizeof(BigDFT_Inputs));
  in->refCount = 1;
  F90_1D_POINTER_INIT(&in->qmass);

  return in;
}


static void bigdft_inputs_dispose(BigDFT_Inputs *in)
{
  if (F_TYPE(in->data))
    FC_FUNC_(inputs_free, INPUTS_FREE)(&in->data);

  _free_output(in);

  g_free(in);
}
BigDFT_Inputs* bigdft_inputs_new_from_fortran(f90_input_variables_pointer inputs)






{
  BigDFT_Inputs *in;

  in = bigdft_inputs_init();
  in->data = inputs;

  _inputs_sync(in);

  return in;
}
void bigdft_inputs_free(BigDFT_Inputs *in)
{
  bigdft_inputs_dispose(in);
}
BigDFT_Inputs* bigdft_inputs_ref(BigDFT_Inputs *in)
{
  in->refCount += 1;
  return in;
}
void bigdft_inputs_unref(BigDFT_Inputs *in)
{
  in->refCount -= 1;
  if (!in->refCount)
    bigdft_inputs_free(in);
}
#ifdef GLIB_MAJOR_VERSION
GType bigdft_inputs_get_type(void)
{
  static GType g_define_type_id = 0;

  if (g_define_type_id == 0)
    g_define_type_id =
      g_boxed_type_register_static("BigDFT_Inputs", 
                                   (GBoxedCopyFunc)bigdft_inputs_ref,
                                   (GBoxedFreeFunc)bigdft_inputs_unref);
  return g_define_type_id;
}
#endif


/* Wrappers on dictionaries, for the input variables. */
void bigdft_inputs_set(BigDFT_Inputs *in, const gchar *level,
                       const gchar *id, const gchar *value)
{
  BigDFT_Dict *dict;

  dict = bigdft_dict_new(NULL);
  bigdft_dict_set(dict, id, value);  
  FC_FUNC_(inputs_set_dict, INPUTS_SET_DICT)(F_TYPE(in->data), level, &dict->root, strlen(level));
  g_object_unref(G_OBJECT(dict));

  _inputs_sync(in);
}
/**
 * bigdft_inputs_set_array:
 * @in: 
 * @id: 
 * @value: (array zero-terminated=1):
 *
 * 
 **/
void bigdft_inputs_set_array(BigDFT_Inputs *in, const gchar *level,
                             const gchar *id, const gchar **value)
{
  BigDFT_Dict *dict;

  dict = bigdft_dict_new(NULL);
  bigdft_dict_set_array(dict, id, value);  
  FC_FUNC_(inputs_set_dict, INPUTS_SET_DICT)(F_TYPE(in->data), level, &dict->root, strlen(level));
  g_object_unref(G_OBJECT(dict));

  _inputs_sync(in);
}
/**
 * bigdft_inputs_set_array_at:
 * @in: 
 * @n_row: 
 * @n_cols:
 * @value: (array zero-terminated=1):
 *
 * 
 **/
void bigdft_inputs_set_matrix(BigDFT_Inputs *in, const gchar *id,
                              guint n_row, guint n_cols, const gchar **value)
{
  /* _dictionary *dict; */
  /* guint i; */

  /* FC_FUNC_(dict_new, DICT_NEW)(&dict); */
  
  /* for (i = 0; value[i]; i++) */
  /*   FC_FUNC_(dict_set_at, DICT_SET_AT)(&dict, id, (int*)&i, value[i], */
  /*                                      strlen(id), strlen(value[i])); */
  /* FC_FUNC_(inputs_set_dict, INPUTS_SET_DICT)(F_TYPE(in->data), &dict); */

  /* FC_FUNC_(dict_free, DICT_FREE)(&dict); */

  /* _inputs_sync(in); */
}
