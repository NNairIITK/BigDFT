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
  BigDFT_Goutput *outs = BIGDFT_GOUTPUT(obj);
  long self = 0;

  if (outs->dispose_has_run)
    return;
  outs->dispose_has_run = TRUE;
  /* Detach the wrapper from the Fortran object. */
  if (F_TYPE(outs->data))
    FC_FUNC_(state_properties_set_c_obj, STATE_PROPERTIES_SET_C_OBJ)
      (F_TYPE(outs->data), &self);

#ifdef HAVE_GLIB
  /* Chain up to the parent class */
  G_OBJECT_CLASS(bigdft_goutput_parent_class)->dispose(obj);
#endif
}
static void bigdft_goutput_finalize(GObject *obj)
{
  BigDFT_Goutput *outs = BIGDFT_GOUTPUT(obj);

  if (F_TYPE(outs->data))
    FC_FUNC_(state_properties_delete, STATE_PROPERTIES_DELETE)(&outs->data);

#ifdef HAVE_GLIB
  G_OBJECT_CLASS(bigdft_goutput_parent_class)->finalize(obj);
#endif
}
static void _sync_energs(BigDFT_Goutput *outs)
{
  FC_FUNC_(energs_copy_data, ENERGS_COPY_DATA)
    (F_TYPE(outs->energs), &outs->eh, &outs->exc,
     &outs->evxc, &outs->eion, &outs->edisp,
     &outs->ekin, &outs->epot, &outs->eproj,
     &outs->eexctX, &outs->ebs, &outs->eKS,
     &outs->trH, &outs->evsum, &outs->evsic);
}
static void _sync_outs(BigDFT_Goutput *outs)
{
  f90_pointer_double_2D fxyz;
  long self;

  F90_2D_POINTER_INIT(&fxyz);
  FC_FUNC_(state_properties_get, STATE_PROPERTIES_GET)(F_TYPE(outs->data),
                                                       &outs->energs, &fxyz,
                                                       (int*)&outs->fdim,
                                                       &outs->fnoise,
                                                       &outs->pressure,
                                                       outs->strten,
                                                       &outs->etot);
  outs->fxyz = fxyz.data;
  _sync_energs(outs);

  self = *((long*)&outs);
  FC_FUNC_(state_properties_set_c_obj, STATE_PROPERTIES_SET_C_OBJ)
    (F_TYPE(outs->data), &self);
}
static BigDFT_Goutput* bigdft_goutput_from_fortran(f90_state_properties *fobj)
{
  BigDFT_Goutput *outs;
  gpointer cobj;

  FC_FUNC_(state_properties_c_obj, STATE_PROPERTIES_C_OBJ)(fobj, (long*)&cobj);
  if (cobj)
    {
      /* Read existing wrapper. */
      outs = BIGDFT_GOUTPUT(cobj);
      _sync_outs(outs);
      g_object_ref(G_OBJECT(outs));
    }
  else
    {
      /* Create a new wrapper. */
#ifdef HAVE_GLIB
      outs = BIGDFT_GOUTPUT(g_object_new(BIGDFT_GOUTPUT_TYPE, NULL));
#else
      outs = g_malloc(sizeof(BigDFT_Goutput));
      bigdft_goutput_init(outs);
#endif
      F_TYPE(outs->data) = fobj;
      _sync_outs(outs);
    }

  return outs;
}
BigDFT_Goutput* bigdft_goutput_new(guint nat)
{
  f90_state_properties_pointer fobj;
  BigDFT_Goutput *outs;
  
  FC_FUNC_(state_properties_alloc, STATE_PROPERTIES_ALLOC)(&fobj, (int*)&nat);
  outs = bigdft_goutput_from_fortran(F_TYPE(fobj));
  outs->data = fobj; /* Need to copy the full structure here for later
                        deallocate. */
  return outs;
}
/**
 * bigdft_goutput_new_from_fortran:
 * @fadd: 
 *
 * Create a C wrapper around a Fortran object.
 *
 * Returns: (transfer none):
 **/
BigDFT_Goutput* bigdft_goutput_new_from_fortran(long fadd)
{
  gpointer fobj;

  fobj = fadd;
  return bigdft_goutput_from_fortran((f90_state_properties*)fobj);
}
void FC_FUNC_(state_properties_to_python, STATE_PROPERTIES_TO_PYTHON)(long outs, const char *lbl, int *len)
{
#ifdef HAVE_PYTHON
  char buf[64];

  memcpy(buf, lbl, sizeof(char) * *len);
  sprintf(buf + *len, " = BigDFT.Goutput.new_from_fortran(%ld)\n", outs);
  PyRun_SimpleString(buf);
#endif
}
void FC_FUNC_(state_properties_wrapper_detach, STATE_PROPERTIES_WRAPPER_DETACH)(gpointer *obj, int *claim)
{
  BigDFT_Goutput *outs = BIGDFT_GOUTPUT(*obj);
  long self = 0;

  claim = bigdft_get_count(G_OBJECT(outs)) - 1;
  /* Detach the wrapper from the Fortran object. */
  if (!claim)
    {
      FC_FUNC_(state_properties_set_c_obj, STATE_PROPERTIES_SET_C_OBJ)
        (F_TYPE(outs->data), &self);
      F_TYPE(outs->data) = (f90_state_properties*)0;
    }
  bigdft_goutput_unref(outs);
}
void FC_FUNC_(state_properties_wrapper_attach, STATE_PROPERTIES_WRAPPER_ATTACH)(gpointer *obj, f90_state_properties *fobj)
{
  BigDFT_Goutput *outs = BIGDFT_GOUTPUT(*obj);

  FC_FUNC_(state_properties_copy, STATE_PROPERTIES_COPY)(&outs->data, fobj);
  _sync_outs(outs);
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

/**
 * bigdft_memory_ref:
 * @boxed: (transfer full):
 *
 * returns: (transfer full):
 */
BigDFT_Memory* bigdft_memory_ref(BigDFT_Memory *boxed)
{
  boxed->ref += 1;
  return boxed;
}
/**
 * bigdft_memory_unref:
 * @boxed: (transfer full):
 *
 */
void bigdft_memory_unref(BigDFT_Memory *boxed)
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
                                                    (GBoxedCopyFunc)bigdft_memory_ref,
                                                    (GBoxedFreeFunc)bigdft_memory_unref);
  return g_define_type_id;
}
#endif
void bigdft_memory_dump(BigDFT_Memory *mem)
{
  FC_FUNC_(print_memory_estimation, PRINT_MEMORY_ESTIMATION)(F_TYPE(mem->data));

}

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
  long self = 0;

  if (run->dispose_has_run)
    return;
  run->dispose_has_run = TRUE;

  /* Release ownership of atoms and inputs. */
  FC_FUNC_(run_objects_nullify_volatile, RUN_OBJECTS_NULLIFY_VOLATILE)
    (F_TYPE(run->data));
  bigdft_inputs_unref(run->inputs);
  bigdft_atoms_unref(run->atoms);

  /* Release ownership of dict. */
  FC_FUNC_(run_objects_nullify_dict, RUN_OBJECTS_NULLIFY_DICT)(F_TYPE(run->data));
  bigdft_dict_unref(run->dict);

  /* Detach the wrapper from the Fortran object. */
  if (F_TYPE(run->data))
    FC_FUNC_(run_objects_set_c_obj, RUN_OBJECTS_SET_C_OBJ)
      (F_TYPE(run->data), &self);

#ifdef HAVE_GLIB
  /* Chain up to the parent class */
  G_OBJECT_CLASS(bigdft_run_parent_class)->dispose(obj);
#endif
}
static void bigdft_run_finalize(GObject *obj)
{
  BigDFT_Run *run = BIGDFT_RUN(obj);

  if (F_TYPE(run->data))
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
  f90_atoms_data_pointer atoms;
  f90_input_variables_pointer inputs;
  f90_dictionary_pointer dict;

  /* Create C wrappers for Fortran objects. */
  FC_FUNC_(run_objects_get, RUN_OBJECTS_GET)(F_TYPE(run->data),
                                             &dict, &inputs, &atoms);
  run->dict   = bigdft_dict_new_from_fortran(dict);
  run->inputs = bigdft_inputs_new_from_fortran(inputs);
  run->atoms  = bigdft_atoms_new_from_fortran(atoms);
}
BigDFT_Run* bigdft_run_new_from_files(const gchar *radical, const gchar *posinp)
{
  BigDFT_Run *run;

  run = bigdft_run_new();

  /* Call the creation routine of Fortran. */
  FC_FUNC_(run_objects_init_from_run_name, RUN_OBJECTS_INIT_FROM_RUN_NAME)
    (F_TYPE(run->data), radical, posinp, strlen(radical), strlen(posinp));

  _attributes_from_fortran(run);

  return run;
}
/**
 * bigdft_run_new_from_dict:
 * @dict: 
 *
 * Pouet.
 *
 * Returns: (transfer full):
 **/
BigDFT_Run* bigdft_run_new_from_dict(BigDFT_Dict *dict)
{
  BigDFT_Run *run;

  run = bigdft_run_new();

  /* Associate the dictionary and parse it. */
  FC_FUNC_(run_objects_update, RUN_OBJECTS_UPDATE)(F_TYPE(run->data), &dict->root);
  _attributes_from_fortran(run);
  return run;
}

/**
 * bigdft_run_new_from_fortran:
 * @fadd: 
 *
 * Create a C wrapper around a Fortran object.
 *
 * Returns: (transfer none):
 **/
BigDFT_Run* bigdft_run_new_from_fortran(long fadd)
{
  gpointer fobj, cobj;
  BigDFT_Run *run;

  fobj = fadd;
  FC_FUNC_(run_objects_c_obj, RUN_OBJECTS_C_OBJ)(fobj, (long*)&cobj);
  if (cobj)
    {
      /* Read existing wrapper. */
      run = BIGDFT_RUN(cobj);
      g_object_ref(G_OBJECT(run));
    }
  else
    {
      /* Create a new wrapper. */
#ifdef HAVE_GLIB
      run = BIGDFT_RUN(g_object_new(BIGDFT_RUN_TYPE, NULL));
#else
      run = g_malloc(sizeof(BigDFT_Run));
      bigdft_run_init(run);
#endif
      F_TYPE(run->data) = fobj;
      _attributes_from_fortran(run);
    }

  return run;
}
void FC_FUNC_(run_objects_to_python, RUN_OBJECTS_TO_PYTHON)(long run, const char *lbl, int *len)
{
#ifdef HAVE_PYTHON
  char buf[64];

  memcpy(buf, lbl, sizeof(char) * *len);
  sprintf(buf + *len, " = BigDFT.Run.new_from_fortran(%ld)\n", run);
  PyRun_SimpleString(buf);
#endif
}
void FC_FUNC_(run_objects_wrapper_detach, RUN_OBJECTS_WRAPPER_DETACH)(gpointer *obj, int *claim)
{
  BigDFT_Run *outs = BIGDFT_RUN(*obj);
  long self = 0;

  claim = bigdft_get_count(G_OBJECT(outs)) - 1;
  /* Detach the wrapper from the Fortran object. */
  if (!claim)
    {
      FC_FUNC_(run_objects_set_c_obj, RUN_OBJECTS_SET_C_OBJ)
        (F_TYPE(outs->data), &self);
      F_TYPE(outs->data) = (f90_run_objects*)0;
    }
  bigdft_run_unref(outs);
}
void FC_FUNC_(run_objects_wrapper_attach, RUN_OBJECTS_WRAPPER_ATTACH)(gpointer *obj, f90_run_objects *fobj)
{
  BigDFT_Run *run = BIGDFT_RUN(*obj);

  FC_FUNC_(run_objects_copy, RUN_OBJECTS_COPY)(&run->data, fobj);
  _attributes_from_fortran(run);
}
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
void bigdft_run_update(BigDFT_Run *run, BigDFT_Dict *dict)
{
  f90_atoms_data_pointer atoms;
  f90_input_variables_pointer inputs;
  f90_dictionary_pointer dictf;

  if (run->inputs && run->atoms)
    {
      /* Internal atoms, inputs structures will change, need to update
         their containers. */
      FC_FUNC_(run_objects_nullify_volatile, RUN_OBJECTS_NULLIFY_VOLATILE)
        (F_TYPE(run->data));
      bigdft_inputs_unref(run->inputs);
      bigdft_atoms_unref(run->atoms);
    }
  FC_FUNC_(run_objects_update, RUN_OBJECTS_UPDATE)(F_TYPE(run->data), &dict->root);
  /* Reassociate atoms, inputs structures. */
  FC_FUNC_(run_objects_get, RUN_OBJECTS_GET)(F_TYPE(run->data),
                                             &dictf, &inputs, &atoms);
  run->inputs = bigdft_inputs_new_from_fortran(inputs);
  run->atoms  = bigdft_atoms_new_from_fortran(atoms);
  if (!run->dict)
    run->dict = bigdft_dict_new_from_fortran(dictf);
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
  int ln;
  ln=strlen(filename);
  /*fprintf(stdout,"length %i\n",ln);*/
  FC_FUNC_(run_objects_dump_to_file, RUN_OBJECTS_DUMP_TO_FILE)
    (&iostat, &run->dict->root, filename, &userOnly, &ln, ln);
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
    (F_TYPE(run->data), (int*)&iproc, (int*)&nproc, rxyz, shift, F_TYPE(mem->data));
  g_free(rxyz);
  FC_FUNC_(mem_to_c, MEM_TO_C)(F_TYPE(mem->data), &mem->submat, &mem->ncomponents,
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
  int infocode, ncount_bigdft;
  BigDFT_Goutput *outs;

  outs = bigdft_goutput_new(run->atoms->nat);
  FC_FUNC_(bigdft_exec, BIGDFT_EXEC)(F_TYPE(run->data), F_TYPE(outs->data),&infocode);
  _inputs_sync(run->inputs);

  if (run->inputs->ncount_cluster_x > 1)
    {
      FC_FUNC(geopt, GEOPT)(F_TYPE(run->data), F_TYPE(outs->data),
                            (int*)&nproc, (int*)&iproc, &ncount_bigdft);
      _inputs_sync(run->inputs);
    }

  /* if there is a last run to be performed do it now before stopping */
  if (run->inputs->last_run == -1)
    {
      FC_FUNC_(bigdft_exec, BIGDFT_EXEC)(F_TYPE(run->data), F_TYPE(outs->data),
                                         &infocode);
      _inputs_sync(run->inputs);
    }

  _sync_outs(outs);
  
  return outs;
}

/**
 * bigdft_run_get_dict:
 * @run: 
 *
 * Pouet.
 *
 * Returns: (transfer full):
 **/
BigDFT_Dict* bigdft_run_get_dict(BigDFT_Run *run)
{
  g_object_ref(G_OBJECT(run->dict));
  return run->dict;
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
