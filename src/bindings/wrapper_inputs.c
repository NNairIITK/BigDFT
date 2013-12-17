#include <config.h>

#include "bigdft.h"
#include "bindings.h"
#include "bindings_api.h"

#include <string.h>
#include <stdio.h>

static void _free_names(BigDFT_Inputs *in)
{
  g_free(in->run_name);

  g_free(in->file_occnum);
  g_free(in->file_igpop);
  g_free(in->file_lin);
}
static void _free_output(BigDFT_Inputs *in)
{
  g_free(in->dir_output);
  g_free(in->writing_directory);
}
static void _sync_output(BigDFT_Inputs *in)
{
  gchar dir_output[100], writing_directory[500];

  FC_FUNC_(inputs_get_output, INPUTS_GET_OUTPUT)(in->data, dir_output, writing_directory, 100, 500);
  in->dir_output = _get_c_string(dir_output, 100);
  in->writing_directory = _get_c_string(writing_directory, 500);
}
static void _sync(BigDFT_Inputs *in)
{
  FC_FUNC_(inputs_get_files, INPUTS_GET_FILES)(in->data, &in->files);
  FC_FUNC_(inputs_get_dft, INPUTS_GET_DFT)(in->data, in->h, in->h + 1, in->h + 2,
                                           &in->crmult, &in->frmult, &in->ixc,
                                           &in->ncharge, in->elecfield, &in->nspin,
                                           &in->mpol, &in->gnrm_cv, (int*)&in->itermax,
                                           (int*)&in->nrepmax, &in->ncong, (int*)&in->idsx,
                                           &in->dispersion, &in->inputPsiId,
                                           &in->output_wf_format, &in->output_grid,
                                           &in->rbuf, &in->ncongt, &in->norbv, &in->nvirt,
                                           &in->nplot, &in->disableSym);
  FC_FUNC_(inputs_get_mix, INPUTS_GET_MIX)(in->data, (int*)&in->iscf, (int*)&in->itrpmax,
                                           (int*)&in->norbsempty, (int*)(&in->occopt),
                                           &in->alphamix,
                                           &in->rpnrm_cv, &in->gnrm_startmix, &in->Tel,
                                           &in->alphadiis);
  FC_FUNC_(inputs_get_geopt, INPUTS_GET_GEOPT)(in->data, in->geopt_approach,
                                               &in->ncount_cluster_x, &in->frac_fluct,
                                               &in->forcemax, &in->randdis, &in->betax,
                                               &in->history, &in->ionmov, &in->dtion,
                                               in->strtarget, &in->qmass, 10);
  /* FC_FUNC_(inputs_get_sic, INPUTS_GET_SIC)(); */
  /* FC_FUNC_(inputs_get_tddft, INPUTS_GET_TDDFT)(); */
  FC_FUNC_(inputs_get_perf, INPUTS_GET_PERF)(in->data, (int*)&in->linear);
}
static void _sync_add(BigDFT_Inputs *in)
{
  FC_FUNC_(inputs_get_files, INPUTS_GET_FILES)(in->data, &in->files);
  /* FC_FUNC_(inputs_get_kpt, INPUTS_GET_KPT)(); */
  _sync_output(in);
}

static BigDFT_Inputs* bigdft_inputs_init()
{
  BigDFT_Inputs *in;

  in = g_malloc(sizeof(BigDFT_Inputs));
  memset(in, 0, sizeof(BigDFT_Inputs));
  in->refCount = 1;
  F90_1D_POINTER_INIT(&in->qmass);
  in->files = BIGDFT_INPUTS_UNPARSED;

  return in;
}
static void bigdft_inputs_dispose(BigDFT_Inputs *in)
{
  if (in->data)
    FC_FUNC_(inputs_free, INPUTS_FREE)(&in->data);
  if (in->input_values)
    FC_FUNC_(dict_free, DICT_FREE)(&in->input_values);

  _free_names(in);
  _free_output(in);

  g_free(in);
}
/**
 * bigdft_inputs_new:
 * @naming: (allow-none): a naming scheme, or none.
 *
 * Create a new #BigDFT_Inputs structure, empty.
 * 
 * Returns: (transfer full): a new structure.
 */
BigDFT_Inputs* bigdft_inputs_new(const gchar *naming)
{
  BigDFT_Inputs *in;
  gchar file_occnum[100], file_igpop[100], file_lin[100], run_name[100];

  in = bigdft_inputs_init();
  FC_FUNC_(inputs_new, INPUTS_NEW)(&in->data);
  FC_FUNC_(dict_new, DICT_NEW)(&in->input_values);
  
  if (naming && naming[0])
    FC_FUNC_(standard_inputfile_names, STANDARD_INPUTFILE_NAMES)(in->data, naming, strlen(naming));
  else
    FC_FUNC_(standard_inputfile_names, STANDARD_INPUTFILE_NAMES)(in->data, " ", 1);
  /* Get naming schemes. */
  _free_names(in);
  FC_FUNC_(inputs_get_naming, INPUTS_GET_NAMING)(in->data, run_name, file_occnum, file_igpop,
                                                 file_lin, 100, 100, 100, 100);
  in->run_name = _get_c_string(run_name, 100);
  in->file_occnum = _get_c_string(file_occnum, 100);
  in->file_igpop = _get_c_string(file_igpop, 100);
  in->file_lin = _get_c_string(file_lin, 100);

  return in;
}
/**
 * bigdft_inputs_new_from_files:
 * @naming: (allow-none): a naming scheme, or none.
 * @iproc:
 *
 * Create a new #BigDFT_Inputs structure, parsing files with the given
 * naming scheme. Use bigdft_inputs_analyse() to actually analyse the
 * values and transfer them into variables.
 * 
 * Returns: (transfer full): a new structure.
 */
BigDFT_Inputs* bigdft_inputs_new_from_files(const gchar *naming, guint iproc)
{
  BigDFT_Inputs *in;

  in = bigdft_inputs_new(naming);

  if (naming && naming[0])
    FC_FUNC_(inputs_set_from_file, INPUTS_SET_FROM_FILE)
      (&in->input_values, naming, strlen(naming));
  else
    FC_FUNC_(inputs_set_from_file, INPUTS_SET_FROM_FILE)
      (&in->input_values, " ", 1);
  
  return in;
}
BigDFT_Inputs* bigdft_inputs_new_from_fortran(_input_variables *inputs)
{
  BigDFT_Inputs *in;

  in = bigdft_inputs_init();
  in->data = inputs;

  _sync(in);
  _sync_add(in);

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
void bigdft_inputs_analyse(BigDFT_Inputs *in, BigDFT_Atoms *atoms, gboolean dump)
{
  FC_FUNC_(inputs_from_dict, INPUTS_FROM_DICT)(in->data, atoms->data, &in->input_values, (gint*)&dump);
  _sync(in);
  _sync_add(in);
  /* To be removed later, currently, this allocates atoms also. */
  bigdft_atoms_get_nat_arrays(atoms);
  bigdft_atoms_get_ntypes_arrays(atoms);
}
void bigdft_inputs_create_dir_output(BigDFT_Inputs *in, guint iproc)
{
  FC_FUNC_(create_dir_output, CREATE_DIR_OUTPUT)((int*)&iproc, in->data);
  _sync_output(in);
}
gboolean bigdft_inputs_dump(BigDFT_Inputs *in, const gchar *filename, gboolean useronly)
{
  int iostat;

  FC_FUNC_(inputs_dump_to_file, INPUTS_DUMP_TO_FILE)(&iostat,
                                                     &in->input_values,
                                                     filename, &useronly,
                                                     strlen(filename));
  return (iostat != 0);
}

/**
 * bigdft_set_input:
 * @radical: 
 * @posinp: 
 * @atoms: (out) (transfer full):
 *
 * Pouet.
 *
 * Returns: (transfer full):
 **/
BigDFT_Inputs* bigdft_set_input(const gchar *radical, const gchar *posinp, BigDFT_Atoms **atoms)
{
  BigDFT_Atoms *at;
  BigDFT_Inputs *in;

  at = bigdft_atoms_new();
  in = bigdft_inputs_init();
  FC_FUNC_(inputs_new, INPUTS_NEW)(&in->data);
  FC_FUNC_(bigdft_set_input, BIGDFT_SET_INPUT)(radical, posinp, in->data, at->data,
                                               strlen(radical), strlen(posinp));
  _sync(in);
  _sync_add(in);
  bigdft_atoms_copy_from_fortran(at);
  *atoms = at;
  return in;
}

/* Wrappers on dictionaries, for the input variables. */
#include "input_keys.h"
void bigdft_inputs_set(BigDFT_Inputs *in, BigDFT_InputsKeyIds id, const gchar *value)
{
  const gchar *name, *file;

  name = _input_keys[id];
  file = _input_keys[_input_files[id]];
  FC_FUNC_(inputs_set, INPUTS_SET)(&in->input_values, file, name, value,
                                   strlen(file), strlen(name), strlen(value));
}
/**
 * bigdft_inputs_set_array:
 * @in: 
 * @id: 
 * @value: (array zero-terminated=1):
 *
 * 
 **/
void bigdft_inputs_set_array(BigDFT_Inputs *in, BigDFT_InputsKeyIds id,
                             const gchar **value)
{
  const gchar *name, *file;
  guint i;

  name = _input_keys[id];
  file = _input_keys[_input_files[id]];
  for (i = 0; value[i]; i++)
    FC_FUNC_(inputs_set_at, INPUTS_SET_AT)(&in->input_values, file, name, (gint*)&i, value[i],
                                           strlen(file), strlen(name), strlen(value[i]));
}
/**
 * bigdft_inputs_set_array_at:
 * @in: 
 * @id: 
 * @at:
 * @value: (array zero-terminated=1):
 *
 * 
 **/
void bigdft_inputs_set_array_at(BigDFT_Inputs *in, BigDFT_InputsKeyIds id,
                                guint at, const gchar **value)
{
  const gchar *name, *file;
  guint i;

  name = _input_keys[id];
  file = _input_keys[_input_files[id]];
  for (i = 0; value[i]; i++)
    FC_FUNC_(inputs_set_at2, INPUTS_SET_AT2)(&in->input_values, file, name,
                                             (gint*)&at, (gint*)&i, value[i],
                                             strlen(file), strlen(name), strlen(value[i]));
}
