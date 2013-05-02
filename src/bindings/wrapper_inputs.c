#include <config.h>

#include "bigdft.h"
#include "bindings.h"
#include "bindings_api.h"

#include <string.h>

static void _inputs_sync(BigDFT_Inputs *in)
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
static void _inputs_sync_add(BigDFT_Inputs *in)
{
  FC_FUNC_(inputs_get_files, INPUTS_GET_FILES)(in->data, &in->files);
  /* FC_FUNC_(inputs_get_kpt, INPUTS_GET_KPT)(); */
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
  if (in->data)
    FC_FUNC_(inputs_free, INPUTS_FREE)(&in->data);
  g_free(in);
}
/**
 * bigdft_inputs_new:
 * @naming: (allow-none): a naming scheme, or none.
 *
 * Create a new #BigDFT_Inputs structure, issuing into default values.
 * 
 * Returns: (transfer full): a new structure.
 */
BigDFT_Inputs* bigdft_inputs_new(const gchar *naming)
{
  BigDFT_Inputs *in;
  int len;

  in = bigdft_inputs_init();
  FC_FUNC_(inputs_new, INPUTS_NEW)(&in->data);
  if (naming && naming[0])
    {
      len = strlen(naming);
      FC_FUNC_(inputs_set_radical, INPUTS_SET_RADICAL)(in->data, naming, &len, len);
    }
  else
    {
      len = 1;
      FC_FUNC_(inputs_set_radical, INPUTS_SET_RADICAL)(in->data, " ", &len, 1);
    }
  
  return in;
}
BigDFT_Inputs* bigdft_inputs_new_from_fortran(_input_variables *inputs)
{
  BigDFT_Inputs *in;

  in = bigdft_inputs_init();
  in->data = inputs;

  _inputs_sync(in);
  _inputs_sync_add(in);

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
void bigdft_inputs_parse(BigDFT_Inputs *in, guint iproc)
{
  int dump = 0;

  FC_FUNC_(inputs_parse_params, INPUTS_PARSE_PARAMS)(in->data, (int*)&iproc, &dump);
  _inputs_sync(in);
}
void bigdft_inputs_parse_additional(BigDFT_Inputs *in, BigDFT_Atoms *atoms, guint iproc)
{
  int dump = 0;

  FC_FUNC_(inputs_parse_add, INPUTS_PARSE_ADD)(in->data, atoms->data, (int*)&iproc, &dump);
  _inputs_sync_add(in);
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
  FC_FUNC_(bigdft_set_input, BIGDFT_SET_INPUT)(radical, posinp, &at->rxyz,
                                               in->data, at->data, strlen(radical), strlen(posinp));

  bigdft_atoms_copy_from_fortran(at);
  *atoms = at;
  return in;
}

