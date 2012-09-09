/**************************************/
/* BigDFT_LocalFields data structure. */
/**************************************/

#include <config.h>

#ifdef HAVE_GLIB
#include <glib-object.h>
#include <gio/gio.h>
#endif

#include "bigdft.h"
#include "bindings.h"
#include "bindings_api.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

static void bigdft_localfields_dispose(GObject *localfields);
static void bigdft_localfields_finalize(GObject *localfields);

#ifdef HAVE_GLIB
enum {
  V_EXT_READY_SIGNAL,
  DENSITY_READY_SIGNAL,
  LAST_SIGNAL
};

G_DEFINE_TYPE(BigDFT_LocalFields, bigdft_localfields, G_TYPE_OBJECT)

static guint bigdft_localfields_signals[LAST_SIGNAL] = { 0 };

static void bigdft_localfields_class_init(BigDFT_LocalFieldsClass *klass)
{
  /* Connect the overloading methods. */
  G_OBJECT_CLASS(klass)->dispose      = bigdft_localfields_dispose;
  G_OBJECT_CLASS(klass)->finalize     = bigdft_localfields_finalize;
  /* G_OBJECT_CLASS(klass)->set_property = visu_data_set_property; */
  /* G_OBJECT_CLASS(klass)->get_property = visu_data_get_property; */

  bigdft_localfields_signals[V_EXT_READY_SIGNAL] =
    g_signal_new("v-ext-ready", G_TYPE_FROM_CLASS(klass),
                 G_SIGNAL_RUN_LAST | G_SIGNAL_NO_RECURSE | G_SIGNAL_NO_HOOKS,
		 0, NULL, NULL, g_cclosure_marshal_VOID__VOID,
                 G_TYPE_NONE, 0, NULL);
  bigdft_localfields_signals[DENSITY_READY_SIGNAL] =
    g_signal_new("density-ready", G_TYPE_FROM_CLASS(klass),
                 G_SIGNAL_RUN_LAST | G_SIGNAL_NO_RECURSE | G_SIGNAL_NO_HOOKS,
		 0, NULL, NULL, g_cclosure_marshal_VOID__UINT,
                 G_TYPE_NONE, 1, G_TYPE_UINT, NULL);
}
#endif

static void bigdft_localfields_init(BigDFT_LocalFields *obj)
{
#ifdef HAVE_GLIB
  memset((void*)((char*)obj + sizeof(GObject)), 0, sizeof(BigDFT_LocalFields) - sizeof(GObject));
#else
  memset(obj, 0, sizeof(BigDFT_LocalFields));
#endif
}
static void bigdft_localfields_dispose(GObject *obj)
{
#ifdef HAVE_GLIB
  BigDFT_LocalFields *denspot = BIGDFT_LOCALFIELDS(obj);

  if (denspot->dispose_has_run)
    return;
  denspot->dispose_has_run = TRUE;

  /* Chain up to the parent class */
  G_OBJECT_CLASS(bigdft_localfields_parent_class)->dispose(obj);
#endif
}
static void bigdft_localfields_finalize(GObject *obj)
{
  BigDFT_LocalFields *denspotd = BIGDFT_LOCALFIELDS(obj);

  if (denspotd->data)
    FC_FUNC_(localfields_free, LOCALFIELDS_FREE)(&denspotd->data);

#ifdef HAVE_GLIB
  G_OBJECT_CLASS(bigdft_localfields_parent_class)->finalize(obj);
  /* g_debug("Freeing localfields object %p done.\n", obj); */
#endif
}
void FC_FUNC_(localfields_emit_rhov, LOCALFIELDS_EMIT_RHOV)(BigDFT_LocalFields **denspot,
                                                            guint *istep)
{
  FC_FUNC_(localfields_copy_metadata, LOCALFIELDS_COPY_METADATA)
    ((*denspot)->data, &(*denspot)->rhov_is, (*denspot)->h,
     (*denspot)->ni, &(*denspot)->psoffset);
  bigdft_localfields_emit_rhov(*denspot, *istep);
}
void bigdft_localfields_emit_rhov(BigDFT_LocalFields *denspot, guint istep)
{
#ifdef HAVE_GLIB
  switch (denspot->rhov_is)
    {
    case BIGDFT_RHO_IS_ELECTRONIC_DENSITY:
      g_signal_emit(G_OBJECT(denspot), bigdft_localfields_signals[DENSITY_READY_SIGNAL],
                    0 /* details */, istep, NULL);
      break;
    default:
      break;
    }
#endif  
}
void FC_FUNC_(localfields_emit_v_ext, LOCALFIELDS_EMIT_V_EXT)(BigDFT_LocalFields **denspot)
{
  bigdft_localfields_emit_v_ext(*denspot);
}
void bigdft_localfields_emit_v_ext(BigDFT_LocalFields *denspot)
{
#ifdef HAVE_GLIB
  g_signal_emit(G_OBJECT(denspot), bigdft_localfields_signals[V_EXT_READY_SIGNAL],
		0 /* details */, NULL);
#endif
}

BigDFT_LocalFields* bigdft_localfields_new(const BigDFT_Lzd *lzd,
                                           const BigDFT_Inputs *in,
                                           guint iproc, guint nproc)
{
  BigDFT_LocalFields *localfields;
  double self;
  gchar SICapproach[4] = "NO  ", rho_commun[3] = "DBL";
  guint verb = 0;
  guint iproc_world, nproc_world, iproc_grp, nproc_grp, mpi_comm;

#ifdef HAVE_GLIB
  localfields = BIGDFT_LOCALFIELDS(g_object_new(BIGDFT_LOCALFIELDS_TYPE, NULL));
#else
  localfields = g_malloc(sizeof(BigDFT_LocalFields));
  memset(localfields, 0, sizeof(BigDFT_LocalFields));
  bigdft_localfields_init(localfields);
#endif
  self = *((double*)&localfields);
  FC_FUNC_(localfields_new, LOCALFIELDS_NEW)(&self, &localfields->data,
                                             &localfields->rhod, &localfields->dpbox);
  FC_FUNC_(initialize_dft_local_fields, INITIALIZE_DFT_LOCAL_FIELDS)(localfields->data);

  FC_FUNC_(dpbox_set_box, DPBOX_SET_BOX)(localfields->dpbox, lzd->data);
  /* We get hgrids and ni. */
  FC_FUNC_(localfields_copy_metadata, LOCALFIELDS_COPY_METADATA)
    (localfields->data, &localfields->rhov_is, localfields->h, localfields->ni,
     &localfields->psoffset);

  FC_FUNC_(system_initkernels, SYSTEM_INITKERNELS)(&verb, &iproc, &nproc,
                                                   &lzd->parent.parent.geocode,
                                                   in->data, localfields->data, 1);
  GET_ATTR_DBL(localfields, LOCALFIELDS, pkernel,    PKERNEL);
  GET_ATTR_DBL(localfields, LOCALFIELDS, pkernelseq, PKERNELSEQ);
  FC_FUNC_(kernel_get_comm, KERNEL_GET_COMM)(localfields->pkernel, &iproc_world,
                                             &nproc_world, &iproc_grp, &nproc_grp,
                                             &mpi_comm);
  FC_FUNC_(denspot_communications, DENSPOT_COMMUNICATIONS)
    (&iproc_world, &nproc_world, &iproc_grp, &nproc_grp, &mpi_comm,
     &in->ixc, &in->nspin, &lzd->parent.parent.geocode,
     SICapproach, localfields->dpbox, 1, 4);
  FC_FUNC_(density_descriptors, DENSITY_DESCRIPTORS)(&iproc, &nproc, &in->nspin,
                                                     &in->crmult, &in->frmult,
                                                     lzd->parent.parent.data,
                                                     localfields->dpbox, rho_commun,
                                                     lzd->parent.parent.rxyz.data,
                                                     lzd->parent.radii, localfields->rhod, 3);
  FC_FUNC(allocaterhopot, ALLOCATERHOPOT)(&iproc, lzd->parent.data,
                                          &in->nspin, lzd->parent.parent.data,
					  lzd->parent.parent.rxyz.data,
                                          localfields->data);
  FC_FUNC_(localfields_copy_metadata, LOCALFIELDS_COPY_METADATA)
    (localfields->data, &localfields->rhov_is, localfields->h, localfields->ni,
     &localfields->psoffset);
  GET_ATTR_DBL   (localfields, LOCALFIELDS, rhov,  RHOV);
  GET_ATTR_DBL_4D(localfields, LOCALFIELDS, v_ext, V_EXT);
  GET_ATTR_DBL_4D(localfields, LOCALFIELDS, v_xc,  V_XC);
  
  return localfields;
}
void FC_FUNC_(localfields_new_wrapper, LOCALFIELDS_NEW_WRAPPER)(double *self, void *obj)
{
  BigDFT_LocalFields *localfields;

  localfields = bigdft_localfields_new_from_fortran(obj);
  *self = *((double*)&localfields);
}
BigDFT_LocalFields* bigdft_localfields_new_from_fortran(void *obj)
{
  BigDFT_LocalFields *localfields;

#ifdef HAVE_GLIB
  localfields = BIGDFT_LOCALFIELDS(g_object_new(BIGDFT_LOCALFIELDS_TYPE, NULL));
#else
  localfields = g_malloc(sizeof(BigDFT_LocalFields));
  bigdft_localfields_init(localfields);
#endif
  localfields->data = obj;
  FC_FUNC_(localfields_get_data, LOCALFIELDS_GET_DATA)
    (localfields->data, &localfields->rhod, &localfields->dpbox);

  FC_FUNC_(localfields_copy_metadata, LOCALFIELDS_COPY_METADATA)
    (localfields->data, &localfields->rhov_is, localfields->h, localfields->ni,
     &localfields->psoffset);
  GET_ATTR_DBL   (localfields, LOCALFIELDS, rhov,  RHOV);
  GET_ATTR_DBL_4D(localfields, LOCALFIELDS, v_ext, V_EXT);
  GET_ATTR_DBL_4D(localfields, LOCALFIELDS, v_xc,  V_XC);

  GET_ATTR_DBL   (localfields, LOCALFIELDS, pkernel,    PKERNEL);
  GET_ATTR_DBL   (localfields, LOCALFIELDS, pkernelseq, PKERNELSEQ);

  return localfields;
}
void FC_FUNC_(localfields_free_wrapper, LOCALFIELDS_FREE_WRAPPER)(gpointer *obj)
{
  BigDFT_LocalFields *localfields = BIGDFT_LOCALFIELDS(*obj);

  localfields->data = (gpointer)0;
  bigdft_localfields_free(localfields);
}
void bigdft_localfields_free(BigDFT_LocalFields *denspot)
{
#ifdef HAVE_GLIB
  g_object_unref(G_OBJECT(denspot));
#else
  bigdft_localfields_finalize(denspot);
  g_free(denspot);
#endif
}
void bigdft_localfields_create_poisson_kernels(BigDFT_LocalFields *localfields)
{
  guint verb = 0;
  
  FC_FUNC_(system_createkernels, SYSTEM_CREATEKERNELS)(localfields->data, &verb);
}
void bigdft_localfields_create_effective_ionic_pot(BigDFT_LocalFields *denspot,
                                                   const BigDFT_Lzd *lzd,
                                                   const BigDFT_Inputs *in,
                                                   guint iproc, guint nproc)
{
  int verb = 0;
  
  FC_FUNC(createeffectiveionicpotential, CREATEEFFECTIVEIONICPOTENTIAL)
    (&iproc, &nproc, &verb, in->data, lzd->parent.parent.data, lzd->parent.parent.rxyz.data,
     lzd->parent.parent.shift, lzd->parent.data, denspot->h, denspot->h + 1, denspot->h + 2,
     denspot->dpbox, denspot->pkernel, denspot->v_ext, in->elecfield,
     &denspot->psoffset);
  
  bigdft_localfields_emit_v_ext(denspot);
}
