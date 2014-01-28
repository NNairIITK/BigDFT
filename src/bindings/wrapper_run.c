#include <config.h>

#include "bigdft.h"
#include "bindings.h"
#include "bindings_api.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

/*********************************/
/* BigDFT_Goutput data structure */
/*********************************/
static void bigdft_goutput_dispose(GObject *energs);
static void bigdft_goutput_finalize(GObject *energs);
#ifdef HAVE_GLIB
enum {
  EKS_READY_SIGNAL,
  LAST_SIGNAL
};

G_DEFINE_TYPE(BigDFT_Goutput, bigdft_goutput, G_TYPE_OBJECT)

static guint bigdft_goutput_signals[LAST_SIGNAL] = { 0 };

static void bigdft_goutput_class_init(BigDFT_GoutputClass *klass)
{
  /* Connect the overloading methods. */
  G_OBJECT_CLASS(klass)->dispose      = bigdft_goutput_dispose;
  G_OBJECT_CLASS(klass)->finalize     = bigdft_goutput_finalize;
  /* G_OBJECT_CLASS(klass)->set_property = visu_data_set_property; */
  /* G_OBJECT_CLASS(klass)->get_property = visu_data_get_property; */

  bigdft_goutput_signals[EKS_READY_SIGNAL] =
    g_signal_new("eks-ready", G_TYPE_FROM_CLASS(klass),
                 G_SIGNAL_RUN_LAST | G_SIGNAL_NO_RECURSE | G_SIGNAL_NO_HOOKS,
		 0, NULL, NULL, g_cclosure_marshal_VOID__UINT,
                 G_TYPE_NONE, 1, G_TYPE_UINT, NULL);
}
#endif

static void bigdft_goutput_init(BigDFT_Goutput *obj)
{
#ifdef HAVE_GLIB
  memset((void*)((char*)obj + sizeof(GObject)), 0, sizeof(BigDFT_Goutput) - sizeof(GObject));
#else
  memset(obj, 0, sizeof(BigDFT_Goutput));
  G_OBJECT(obj)->ref_count = 1;
#endif
}
static void bigdft_goutput_dispose(GObject *obj)
{
  BigDFT_Goutput *energs = BIGDFT_GOUTPUT(obj);

  if (energs->dispose_has_run)
    return;
  energs->dispose_has_run = TRUE;

#ifdef HAVE_GLIB
  /* Chain up to the parent class */
  G_OBJECT_CLASS(bigdft_goutput_parent_class)->dispose(obj);
#endif
}
static void bigdft_goutput_finalize(GObject *obj)
{
  BigDFT_Goutput *outs = BIGDFT_GOUTPUT(obj);

  if (outs->data)
    FC_FUNC_(global_output_free, GLOBAL_OUTPUT_FREE)(&outs->data);

#ifdef HAVE_GLIB
  G_OBJECT_CLASS(bigdft_goutput_parent_class)->finalize(obj);
#endif
}
BigDFT_Goutput* bigdft_goutput_new(guint nat)
{
  BigDFT_Goutput *outs;
  long self;
  f90_pointer_double_2D fxyz;

#ifdef HAVE_GLIB
  outs = BIGDFT_GOUTPUT(g_object_new(BIGDFT_GOUTPUT_TYPE, NULL));
#else
  outs = g_malloc(sizeof(BigDFT_Goutput));
  bigdft_goutput_init(outs);
#endif
  self = *((long*)&outs);
  F90_2D_POINTER_INIT(&fxyz);
  FC_FUNC_(global_output_new, GLOBAL_OUTPUT_NEW)(&self, &outs->data, &outs->energs,
                                                 &fxyz, (int*)&nat);
  outs->fdim = nat;
  outs->fxyz = fxyz.data;

  return outs;
}
void FC_FUNC_(energs_new_wrapper, ENERGS_NEW_WRAPPER)(double *self, _DFT_global_output *obj)
{
  BigDFT_Goutput *outs;

  outs = bigdft_goutput_new_from_fortran(obj);
  *self = *((double*)&outs);
}
static void _sync_energs(BigDFT_Goutput *outs)
{
  FC_FUNC_(energs_copy_data, ENERGS_COPY_DATA)
    (outs->energs, &outs->eh, &outs->exc,
     &outs->evxc, &outs->eion, &outs->edisp,
     &outs->ekin, &outs->epot, &outs->eproj,
     &outs->eexctX, &outs->ebs, &outs->eKS,
     &outs->trH, &outs->evsum, &outs->evsic);
}
static void _sync_outs(BigDFT_Goutput *outs)
{
  f90_pointer_double_2D fxyz;

  F90_2D_POINTER_INIT(&fxyz);
  FC_FUNC_(global_output_get, GLOBAL_OUTPUT_GET)(outs->data, &outs->energs, &fxyz,
                                                 (int*)&outs->fdim,
                                                 &outs->fnoise, &outs->pressure,
                                                 outs->strten, &outs->etot);
  outs->fxyz = fxyz.data;
  _sync_energs(outs);
}
BigDFT_Goutput* bigdft_goutput_new_from_fortran(_DFT_global_output *obj)
{
  BigDFT_Goutput *outs;

#ifdef HAVE_GLIB
  outs = BIGDFT_GOUTPUT(g_object_new(BIGDFT_GOUTPUT_TYPE, NULL));
#else
  outs = g_malloc(sizeof(BigDFT_Goutput));
  bigdft_goutput_init(outs);
#endif
  outs->data = obj;
  _sync_outs(outs);

  return outs;
}
void FC_FUNC_(energs_free_wrapper, ENERGS_FREE_WRAPPER)(gpointer *obj)
{
  BigDFT_Goutput *outs = BIGDFT_GOUTPUT(*obj);

  outs->data = (_DFT_global_output*)0;
  outs->energs = (_energy_terms*)0;
  bigdft_goutput_unref(outs);
}
void bigdft_goutput_unref(BigDFT_Goutput *outs)
{
  g_object_unref(G_OBJECT(outs));
#ifndef HAVE_GLIB
  if (G_OBJECT(outs)->ref_count <= 0)
    {
      bigdft_goutput_dispose(G_OBJECT(outs));
      bigdft_goutput_finalize(G_OBJECT(outs));
      g_free(outs);
    }
#endif
}
void FC_FUNC_(energs_emit, ENERGS_EMIT)(BigDFT_Goutput **obj, guint *istep,
                                        BigDFT_EnergsIds *kind)
{
  BigDFT_Goutput *outs = BIGDFT_GOUTPUT(*obj);

  _sync_energs(outs);
  bigdft_goutput_emit_energs(*obj, *istep, *kind);
}
void bigdft_goutput_emit_energs(BigDFT_Goutput *outs, guint istep, BigDFT_EnergsIds kind)
{
#ifdef HAVE_GLIB
  switch (kind)
    {
    case BIGDFT_ENERGS_EKS:
      g_signal_emit(G_OBJECT(outs), bigdft_goutput_signals[EKS_READY_SIGNAL],
                    0 /* details */, istep, NULL);
      break;
    default:
      break;
    }
#endif  
}

/*********************************/
/* BigDFT_Restart data structure */
/*********************************/
static void bigdft_restart_dispose(GObject *restart);
static void bigdft_restart_finalize(GObject *restart);
#ifdef HAVE_GLIB
G_DEFINE_TYPE(BigDFT_Restart, bigdft_restart, G_TYPE_OBJECT)

static void bigdft_restart_class_init(BigDFT_RestartClass *klass)
{
  /* Connect the overloading methods. */
  G_OBJECT_CLASS(klass)->dispose      = bigdft_restart_dispose;
  G_OBJECT_CLASS(klass)->finalize     = bigdft_restart_finalize;
  /* G_OBJECT_CLASS(klass)->set_property = visu_data_set_property; */
  /* G_OBJECT_CLASS(klass)->get_property = visu_data_get_property; */
}
#endif

static void bigdft_restart_init(BigDFT_Restart *obj)
{
#ifdef HAVE_GLIB
  memset((void*)((char*)obj + sizeof(GObject)), 0, sizeof(BigDFT_Restart) - sizeof(GObject));
#else
  memset(obj, 0, sizeof(BigDFT_Restart));
  G_OBJECT(obj)->ref_count = 1;
#endif
}
static void bigdft_restart_dispose(GObject *obj)
{
  BigDFT_Restart *restart = BIGDFT_RESTART(obj);

  if (restart->dispose_has_run)
    return;
  restart->dispose_has_run = TRUE;

#ifdef HAVE_GLIB
  /* Chain up to the parent class */
  G_OBJECT_CLASS(bigdft_restart_parent_class)->dispose(obj);
#endif
}
static void bigdft_restart_finalize(GObject *obj)
{
  BigDFT_Restart *restart = BIGDFT_RESTART(obj);

  if (restart->data)
    FC_FUNC_(rst_free, RST_FREE)(&restart->data);
  if (restart->in)
    bigdft_inputs_unref(restart->in);

#ifdef HAVE_GLIB
  G_OBJECT_CLASS(bigdft_restart_parent_class)->finalize(obj);
#endif
}
BigDFT_Restart* bigdft_restart_new(BigDFT_Atoms *atoms, BigDFT_Inputs *in, guint iproc)
{
  BigDFT_Restart *restart;
  long self;

#ifdef HAVE_GLIB
  restart = BIGDFT_RESTART(g_object_new(BIGDFT_RESTART_TYPE, NULL));
#else
  restart = g_malloc(sizeof(BigDFT_Restart));
  bigdft_restart_init(restart);
#endif
  self = *((long*)&restart);
  FC_FUNC_(rst_new, RST_NEW)(&self, &restart->data);
  FC_FUNC_(rst_init, RST_INIT)(restart->data, (int*)&iproc, atoms->data, in->data);
  restart->inputPsiId = BIGDFT_RESTART_LCAO;
  restart->in = in;
  bigdft_inputs_ref(in);

  return restart;
}
void FC_FUNC_(restart_new_wrapper, RESTART_NEW_WRAPPER)(double *self, void *obj)
{
  BigDFT_Restart *restart;

  restart = bigdft_restart_new_from_fortran(obj);
  *self = *((double*)&restart);
}
BigDFT_Restart* bigdft_restart_new_from_fortran(_restart_objects *obj)
{
  BigDFT_Restart *restart;

#ifdef HAVE_GLIB
  restart = BIGDFT_RESTART(g_object_new(BIGDFT_RESTART_TYPE, NULL));
#else
  restart = g_malloc(sizeof(BigDFT_Restart));
  bigdft_restart_init(restart);
#endif
  restart->data = obj;

  return restart;
}
void FC_FUNC_(restart_free_wrapper, RESTART_FREE_WRAPPER)(gpointer *obj)
{
  BigDFT_Restart *restart = BIGDFT_RESTART(*obj);

  restart->data = (gpointer)0;
  bigdft_restart_unref(restart);
}
void bigdft_restart_unref(BigDFT_Restart *restart)
{
  g_object_unref(G_OBJECT(restart));
#ifndef HAVE_GLIB
  if (G_OBJECT(restart)->ref_count <= 0)
    {
      bigdft_restart_dispose(G_OBJECT(restart));
      bigdft_restart_finalize(G_OBJECT(restart));
      g_free(restart);
    }
#endif
}
void bigdft_restart_set_mode(BigDFT_Restart *restart, BigDFT_RestartModes id)
{
  int inputPsiId[] = {0, 1};

  restart->inputPsiId = id;
  FC_FUNC_(inputs_set_restart, INPUTS_SET_RESTART)(restart->in->data, inputPsiId + id);
}

static BigDFT_Memory* memory_ref(BigDFT_Memory *boxed)
{
  boxed->ref += 1;
  return boxed;
}
static void memory_unref(BigDFT_Memory *boxed)
{
  boxed->ref -= 1;
  if (!boxed->ref)
    {
      FC_FUNC_(mem_destroy, MEM_DESTROY)(&boxed->data);
      g_free(boxed);
    }
}
#ifdef HAVE_GLIB
/**
 * bigdft_memory_get_type:
 *
 * Plop.
 *
 * Returns: a new type for #BigDFT_Memory structures.
 */
GType bigdft_memory_get_type(void)
{
  static GType g_define_type_id = 0;

  if (g_define_type_id == 0)
    g_define_type_id = g_boxed_type_register_static("BigDFT_Memory", 
                                                    (GBoxedCopyFunc)memory_ref,
                                                    (GBoxedFreeFunc)memory_unref);
  return g_define_type_id;
}
#endif

/*****************************/
/* BigDFT_Run data structure */
/*****************************/
static void bigdft_run_dispose(GObject *run);
static void bigdft_run_finalize(GObject *run);
#ifdef HAVE_GLIB
G_DEFINE_TYPE(BigDFT_Run, bigdft_run, G_TYPE_OBJECT)

static void bigdft_run_class_init(BigDFT_RunClass *klass)
{
  /* Connect the overloading methods. */
  G_OBJECT_CLASS(klass)->dispose      = bigdft_run_dispose;
  G_OBJECT_CLASS(klass)->finalize     = bigdft_run_finalize;
  /* G_OBJECT_CLASS(klass)->set_property = visu_data_set_property; */
  /* G_OBJECT_CLASS(klass)->get_property = visu_data_get_property; */
}
#endif

static void bigdft_run_init(BigDFT_Run *obj)
{
#ifdef HAVE_GLIB
  memset((void*)((char*)obj + sizeof(GObject)), 0, sizeof(BigDFT_Run) - sizeof(GObject));
#else
  memset(obj, 0, sizeof(BigDFT_Run));
  G_OBJECT(obj)->ref_count = 1;
#endif
}
static void bigdft_run_dispose(GObject *obj)
{
  BigDFT_Run *run = BIGDFT_RUN(obj);

  if (run->dispose_has_run)
    return;
  run->dispose_has_run = TRUE;

  bigdft_inputs_unref(run->inputs);
  bigdft_atoms_unref(run->atoms);
  bigdft_restart_unref(run->restart);

  if (run->dict)
    {
      /* Release ownership of dict. */
      FC_FUNC_(run_objects_nullify_dict, RUN_OBJECTS_NULLIFY_DICT)(run->data);
      bigdft_dict_unref(run->dict);
    }

#ifdef HAVE_GLIB
  /* Chain up to the parent class */
  G_OBJECT_CLASS(bigdft_run_parent_class)->dispose(obj);
#endif
}
static void bigdft_run_finalize(GObject *obj)
{
  BigDFT_Run *run = BIGDFT_RUN(obj);

  if (run->data)
    FC_FUNC_(run_objects_destroy, RUN_OBJECTS_DESTROY)(&run->data);

#ifdef HAVE_GLIB
  G_OBJECT_CLASS(bigdft_run_parent_class)->finalize(obj);
#endif
}
BigDFT_Run* bigdft_run_new()
{
  BigDFT_Run *run;
  /* long self; */

#ifdef HAVE_GLIB
  run = BIGDFT_RUN(g_object_new(BIGDFT_RUN_TYPE, NULL));
#else
  run = g_malloc(sizeof(BigDFT_Run));
  bigdft_run_init(run);
#endif
  /* self = *((long*)&run); */
  FC_FUNC_(run_objects_new, RUN_OBJECTS_NEW)(&run->data);

  return run;
}
static void _attributes_from_fortran(BigDFT_Run *run)
{
  _atoms_data *atoms;
  _input_variables *inputs;
  _restart_objects *rst;

  /* Create C wrappers for Fortran objects. */
  FC_FUNC_(run_objects_get, RUN_OBJECTS_GET)(run->data, &inputs, &atoms, &rst);
  run->inputs = bigdft_inputs_new_from_fortran(inputs);
  run->atoms  = bigdft_atoms_new_from_fortran(atoms);
  run->restart = bigdft_restart_new_from_fortran(rst);
  run->restart->in = run->inputs;
  bigdft_inputs_ref(run->inputs);
}
BigDFT_Run* bigdft_run_new_from_files(const gchar *radical, const gchar *posinp)
{
  BigDFT_Run *run;

  run = bigdft_run_new();

  /* Call the creation routine of Fortran. */
  FC_FUNC_(run_objects_init_from_files, RUN_OBJECTS_SET_FROM_FILES)
    (run->data, radical, posinp, strlen(radical), strlen(posinp));

  _attributes_from_fortran(run);

  return run;
}
/**
 * bigdft_run_new_from_dict:
 * @inputs: 
 * @atoms: 
 * @rst: (allow-none):
 * @iproc: 
 *
 * Pouet.
 *
 * Returns: (transfer full):
 **/
BigDFT_Run* bigdft_run_new_from_dict(BigDFT_Dict *dict, gboolean dump)
{
  BigDFT_Run *run;

  run = bigdft_run_new();

  /* Taking a pointer on the dictionary. */
  g_object_ref(G_OBJECT(dict));
  run->dict = dict;

  /* Associate the dictionary and parse it. */
  FC_FUNC_(run_objects_set_dict, RUN_OBJECTS_SET_DICT)(run->data, &dict->root);
  FC_FUNC_(run_objects_parse, RUN_OBJECTS_PARSE)(run->data, (gint*)&dump);

  _attributes_from_fortran(run);

  return run;
}

/* void FC_FUNC_(run_new_wrapper, RUN_NEW_WRAPPER)(double *self, void *obj) */
/* { */
/*   BigDFT_Run *run; */

/*   run = bigdft_run_new_from_fortran(obj); */
/*   *self = *((double*)&run); */
/* } */
BigDFT_Run* bigdft_run_new_from_fortran(_run_objects *obj, gboolean create_wrappers)
{
  BigDFT_Run *run;

#ifdef HAVE_GLIB
  run = BIGDFT_RUN(g_object_new(BIGDFT_RUN_TYPE, NULL));
#else
  run = g_malloc(sizeof(BigDFT_Run));
  bigdft_run_init(run);
#endif
  run->data = obj;

  if (create_wrappers)
    _attributes_from_fortran(run);

  return run;
}
/* void FC_FUNC_(run_free_wrapper, RUN_FREE_WRAPPER)(gpointer *obj) */
/* { */
/*   BigDFT_Run *run = BIGDFT_RUN(*obj); */

/*   run->data = (gpointer)0; */
/*   bigdft_run_free(run); */
/* } */
void bigdft_run_unref(BigDFT_Run *run)
{
  g_object_unref(G_OBJECT(run));
#ifndef HAVE_GLIB
  if (G_OBJECT(run)->ref_count <= 0)
    {
      bigdft_run_dispose(G_OBJECT(run));
      bigdft_run_finalize(G_OBJECT(run));
      g_free(run);
    }
#endif
}
/**
 * bigdft_run_dump:
 * @run: 
 * @filename: (type filename):
 * @full: 
 *
 * Plop.
 *
 * Returns: 
 **/
gboolean bigdft_run_dump(BigDFT_Run *run, const gchar *filename, gboolean full)
{
  int iostat;
  int userOnly = !full;

  FC_FUNC_(run_objects_dump_to_file, RUN_OBJECTS_DUMP_TO_FILE)
    (&iostat, &run->dict->root, filename, &userOnly, strlen(filename));

  return (iostat == 0);
}
/**
 * bigdft_run_memoryEstimation:
 * @run: 
 * @iproc: 
 * @nproc: 
 *
 * Pouet here also.
 *
 * Returns: (transfer full) (type BigDFT_Memory*):
 **/
BigDFT_Memory* bigdft_run_memoryEstimation(BigDFT_Run *run, guint iproc, guint nproc)
{
  double shift[3];
  double *rxyz;
  BigDFT_Memory *mem;

  mem = g_malloc(sizeof(BigDFT_Memory));
  mem->ref = 1;
  FC_FUNC_(mem_new, MEM_NEW)(&mem->data);

  rxyz = g_malloc(sizeof(double) * 3 * run->atoms->nat);
  FC_FUNC_(run_objects_system_setup, RUN_OBJECTS_SYSTEM_SETUP)
    (run->data, (int*)&iproc, (int*)&nproc, rxyz, shift, mem->data);
  g_free(rxyz);
  FC_FUNC_(mem_to_c, MEM_TO_C)(mem->data, &mem->submat, &mem->ncomponents,
                               &mem->norb, &mem->norbp, &mem->oneorb,
                               &mem->allpsi_mpi, &mem->psistorage,
                               &mem->projarr, &mem->grid, &mem->workarr,
                               &mem->kernel, &mem->density, &mem->psolver,
                               &mem->ham, &mem->peak);
  return mem;
}
/**
 * bigdft_run_calculate:
 * @run: 
 * @iproc: 
 * @nproc: 
 *
 * Pouet again.
 *
 * Returns: (transfer full):
 **/
BigDFT_Goutput* bigdft_run_calculate(BigDFT_Run *run, guint iproc, guint nproc)
{
  int infocode;
  BigDFT_Goutput *outs;

  outs = bigdft_goutput_new(run->atoms->nat);
  FC_FUNC_(call_bigdft, CALL_BIGDFT)(run->data, outs->data, (int*)&nproc, (int*)&iproc,
                                     &infocode);
  _sync_outs(outs);
  
  return outs;
}

/**
 * bigdft_run_get_atoms:
 * @run: 
 *
 * Pouet.
 *
 * Returns: (transfer full):
 **/
BigDFT_Atoms* bigdft_run_get_atoms(BigDFT_Run *run)
{
  g_object_ref(G_OBJECT(run->atoms));
  return run->atoms;
}
/**
 * bigdft_run_get_inputs:
 * @run: 
 *
 * Pouet.
 *
 * Returns: (transfer full):
 **/
BigDFT_Inputs* bigdft_run_get_inputs(BigDFT_Run *run)
{
  bigdft_inputs_ref(run->inputs);
  return run->inputs;
}
/**
 * bigdft_run_get_restart:
 * @run: 
 *
 * Pouet.
 *
 * Returns: (transfer full):
 **/
BigDFT_Restart* bigdft_run_get_restart(BigDFT_Run *run)
{
  g_object_ref(G_OBJECT(run->restart));
  return run->restart;
}
