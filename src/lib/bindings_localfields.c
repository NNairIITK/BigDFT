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
  obj->rhod  = (void*)0;
  obj->dpcom = (void*)0;
  obj->data  = (void*)0;
  F90_1D_POINTER_INIT(&obj->rhov);
  F90_1D_POINTER_INIT(&obj->rho_full);
  F90_1D_POINTER_INIT(&obj->pot_full);
  F90_1D_POINTER_INIT(&obj->pkernel);
  F90_1D_POINTER_INIT(&obj->pkernelseq);
  F90_2D_POINTER_INIT(&obj->rho_psi);
  F90_4D_POINTER_INIT(&obj->rho_c);
  F90_4D_POINTER_INIT(&obj->v_ext);
  F90_4D_POINTER_INIT(&obj->v_xc);
  F90_4D_POINTER_INIT(&obj->vloc_ks);
  F90_4D_POINTER_INIT(&obj->f_xc);
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
                                            const BigDFT_Glr *glr,
                                            const BigDFT_Inputs *in,
                                            const double *radii,
                                            guint iproc, guint nproc)
{
  BigDFT_LocalFields *denspot;
  double hh[3];
  guint ndegree_ip = 16, verb = 0;

#ifdef HAVE_GLIB
  denspot = BIGDFT_LOCALFIELDS(g_object_new(BIGDFT_LOCALFIELDS_TYPE, NULL));
#else
  denspot = g_malloc(sizeof(BigDFT_LocalFields));
  memset(denspot, 0, sizeof(BigDFT_LocalFields));
  bigdft_localfields_init(denspot);
#endif
  hh[0] = glr->h[0] * 0.5;
  hh[1] = glr->h[1] * 0.5;
  hh[2] = glr->h[2] * 0.5;
  FC_FUNC_(denspot_communications, DENSPOT_COMMUNICATIONS)
    (&iproc, &nproc, glr->d, hh, hh + 1, hh + 2, in->data,
     atoms->data, atoms->rxyz.data, radii, denspot->dpcom, denspot->rhod);
  FC_FUNC(allocaterhopot, ALLOCATERHOPOT)(&iproc, glr->data,
                                          hh, hh + 1, hh + 2, in->data,
                                          atoms->data, atoms->rxyz.data,
                                          denspot->data);
  FC_FUNC_(localfields_copy_metadata, LOCALFIELDS_COPY_METADATA)
    (denspot->data, &denspot->rhov_is, denspot->h, &denspot->psoffset);
  
  FC_FUNC_(system_createkernels, SYSTEM_CREATEKERNELS)
    (&iproc, &nproc, &verb, &glr->geocode, glr->d, hh,
     in->data, denspot->data);

  return denspot;
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
