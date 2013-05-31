#include <config.h>

#include "bigdft.h"
#include "bindings.h"
#include "bindings_api.h"

#include <string.h>

static void _free_names(BigDFT_Inputs *in)
{
  g_free(in->run_name);

  g_free(in->file_dft);
  g_free(in->file_geopt);
  g_free(in->file_kpt);
  g_free(in->file_perf);
  g_free(in->file_tddft);
  g_free(in->file_mix);
  g_free(in->file_sic);
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

  _free_names(in);
  _free_output(in);

  g_free(in);
}
/**
 * bigdft_inputs_new:
 * @naming: (allow-none): a naming scheme, or none.
 * @iproc:
 *
 * Create a new #BigDFT_Inputs structure, issuing into default values.
 * 
 * Returns: (transfer full): a new structure.
 */
BigDFT_Inputs* bigdft_inputs_new(const gchar *naming, guint iproc)
{
  BigDFT_Inputs *in;
  int len;
  gchar file_dft[100], file_geopt[100], file_kpt[100], file_perf[100], file_tddft[100], file_mix[100],
    file_sic[100], file_occnum[100], file_igpop[100], file_lin[100], run_name[100];

  in = bigdft_inputs_init();
  FC_FUNC_(inputs_new, INPUTS_NEW)(&in->data);
  if (naming && naming[0])
    {
      len = strlen(naming);
      FC_FUNC_(standard_inputfile_names, STANDARD_INPUTFILE_NAMES)(in->data, naming, (int*)&iproc, len);
    }
  else
    {
      len = 1;
      FC_FUNC_(standard_inputfile_names, STANDARD_INPUTFILE_NAMES)(in->data, " ", (int*)&iproc, 1);
    }
  /* Get naming schemes. */
  _free_names(in);
  FC_FUNC_(inputs_get_naming, INPUTS_GET_NAMING)(in->data, run_name, file_dft, file_geopt, file_kpt,
                                                 file_perf, file_tddft, file_mix, file_sic,
                                                 file_occnum, file_igpop, file_lin, 100, 100, 100,
                                                 100, 100, 100, 100, 100, 100, 100, 100);
  in->run_name = _get_c_string(run_name, 100);
  in->file_dft = _get_c_string(file_dft, 100);
  in->file_geopt = _get_c_string(file_geopt, 100);
  in->file_kpt = _get_c_string(file_kpt, 100);
  in->file_perf = _get_c_string(file_perf, 100);
  in->file_tddft = _get_c_string(file_tddft, 100);
  in->file_mix = _get_c_string(file_mix, 100);
  in->file_sic = _get_c_string(file_sic, 100);
  in->file_occnum = _get_c_string(file_occnum, 100);
  in->file_igpop = _get_c_string(file_igpop, 100);
  in->file_lin = _get_c_string(file_lin, 100);
  
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
void bigdft_inputs_parse(BigDFT_Inputs *in, guint iproc, gboolean dump)
{
  FC_FUNC_(inputs_parse_params, INPUTS_PARSE_PARAMS)(in->data, (int*)&iproc, &dump);
  _sync(in);
}
void bigdft_inputs_parse_additional(BigDFT_Inputs *in, BigDFT_Atoms *atoms, guint iproc, gboolean dump)
{
  FC_FUNC_(inputs_parse_add, INPUTS_PARSE_ADD)(in->data, atoms->data, (int*)&iproc, &dump);
  _sync_add(in);
}
void bigdft_inputs_create_dir_output(BigDFT_Inputs *in, guint iproc)
{
  FC_FUNC_(create_dir_output, CREATE_DIR_OUTPUT)((int*)&iproc, in->data);
  _sync_output(in);
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

  bigdft_atoms_copy_from_fortran(at);
  *atoms = at;
  return in;
}

