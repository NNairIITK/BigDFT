/**************************************/
/* BigDFT_LocalFields data structure. */
/**************************************/

#include <config.h>

#ifdef HAVE_GLIB
#include <glib-object.h>
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
G_DEFINE_TYPE(BigDFT_LocalFields, bigdft_localfields, G_TYPE_OBJECT)

static void bigdft_localfields_class_init(BigDFT_LocalFieldsClass *klass)
{
  /* Connect the overloading methods. */
  G_OBJECT_CLASS(klass)->dispose      = bigdft_localfields_dispose;
  G_OBJECT_CLASS(klass)->finalize     = bigdft_localfields_finalize;
  /* G_OBJECT_CLASS(klass)->set_property = visu_data_set_property; */
  /* G_OBJECT_CLASS(klass)->get_property = visu_data_get_property; */
}
#endif

static void bigdft_localfields_init(BigDFT_LocalFields *obj)
{
#ifdef HAVE_GLIB
  memset(obj + sizeof(GObject), 0, sizeof(BigDFT_LocalFields) - sizeof(GObject));
#else
  memset(obj, 0, sizeof(BigDFT_LocalFields));
#endif
  FC_FUNC_(localfields_new, LOCALFIELDS_NEW)(&obj->data, &obj->rhod, &obj->dpcom);
  FC_FUNC_(initialize_dft_local_fields, INITIALIZE_DFT_LOCAL_FIELDS)(obj->data);
}
static void bigdft_localfields_dispose(GObject *obj)
{
  BigDFT_LocalFields *denspot = BIGDFT_LOCALFIELDS(obj);

  if (denspot->dispose_has_run)
    return;
  denspot->dispose_has_run = TRUE;

#ifdef HAVE_GLIB
  g_object_unref(G_OBJECT(denspot->atoms));
  /* Chain up to the parent class */
  G_OBJECT_CLASS(bigdft_localfields_parent_class)->dispose(obj);
#endif
}
static void bigdft_localfields_finalize(GObject *obj)
{
  BigDFT_LocalFields *denspotd = BIGDFT_LOCALFIELDS(obj);

  FC_FUNC_(localfields_free, LOCALFIELDS_FREE)(&denspotd->data);

#ifdef HAVE_GLIB
  G_OBJECT_CLASS(bigdft_localfields_parent_class)->finalize(obj);
#endif
}

BigDFT_LocalFields* bigdft_localfields_new (const BigDFT_Atoms *atoms,
                                            const BigDFT_LocReg *glr,
                                            const BigDFT_Inputs *in,
                                            const double *radii,
                                            guint iproc, guint nproc)
{
  BigDFT_LocalFields *localfields;
  double hh[3];
  guint ndegree_ip = 16, verb = 0;

#ifdef HAVE_GLIB
  localfields = BIGDFT_LOCALFIELDS(g_object_new(BIGDFT_LOCALFIELDS_TYPE, NULL));
  g_object_ref(G_OBJECT(atoms));
#else
  localfields = g_malloc(sizeof(BigDFT_LocalFields));
  memset(localfields, 0, sizeof(BigDFT_LocalFields));
  bigdft_localfields_init(localfields);
#endif
  localfields->atoms = atoms;
  localfields->glr   = glr;

  hh[0] = glr->h[0] * 0.5;
  hh[1] = glr->h[1] * 0.5;
  hh[2] = glr->h[2] * 0.5;
  FC_FUNC_(denspot_communications, DENSPOT_COMMUNICATIONS)
    (&iproc, &nproc, glr->d, hh, hh + 1, hh + 2, in->data,
     atoms->data, atoms->rxyz.data, radii, localfields->dpcom, localfields->rhod);
  FC_FUNC(allocaterhopot, ALLOCATERHOPOT)(&iproc, glr->data,
                                          hh, hh + 1, hh + 2, in->data,
                                          atoms->data, atoms->rxyz.data,
                                          localfields->data);
  FC_FUNC_(localfields_copy_metadata, LOCALFIELDS_COPY_METADATA)
    (localfields->data, &localfields->rhov_is, localfields->h,
     &localfields->psoffset);
  GET_ATTR_DBL   (localfields, LOCALFIELDS, rhov,  RHOV);
  GET_ATTR_DBL_4D(localfields, LOCALFIELDS, v_ext, V_EXT);
  GET_ATTR_DBL_4D(localfields, LOCALFIELDS, v_xc,  V_XC);
  
  FC_FUNC_(system_createkernels, SYSTEM_CREATEKERNELS)
    (&iproc, &nproc, &verb, &glr->geocode, glr->d, hh,
     in->data, localfields->data);
  GET_ATTR_DBL   (localfields, LOCALFIELDS, pkernel,    PKERNEL);
  GET_ATTR_DBL   (localfields, LOCALFIELDS, pkernelseq, PKERNELSEQ);

  return localfields;
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
void bigdft_localfields_create_effective_ionic_pot(BigDFT_LocalFields *denspot,
                                                   const BigDFT_Inputs *in,
                                                   guint iproc, guint nproc)
{
  double hh[3];
  
  hh[0] = denspot->glr->h[0] * 0.5;
  hh[1] = denspot->glr->h[1] * 0.5;
  hh[2] = denspot->glr->h[2] * 0.5;

  FC_FUNC(createeffectiveionicpotential, CREATEEFFECTIVEIONICPOTENTIAL)
    (&iproc, &nproc, in->data, denspot->atoms->data, denspot->atoms->rxyz.data,
     denspot->atoms->shift, denspot->glr->data, hh, hh + 1, hh + 2,
     denspot->dpcom, denspot->pkernel, denspot->v_ext, in->elecfield,
     &denspot->psoffset);
}
