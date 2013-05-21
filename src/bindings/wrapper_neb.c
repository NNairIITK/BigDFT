#include <config.h>

#include "bigdft.h"
#include "bindings.h"
#include "bindings_api.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

/*****************************/
/* BigDFT_Neb data structure */
/*****************************/
#ifdef HAVE_GLIB
G_DEFINE_TYPE(BigDFT_Neb, bigdft_neb, G_TYPE_OBJECT)

static void bigdft_neb_dispose(GObject *neb);
static void bigdft_neb_finalize(GObject *neb);

static void bigdft_neb_class_init(BigDFT_NebClass *klass)
{
  /* Connect the overloading methods. */
  G_OBJECT_CLASS(klass)->dispose      = bigdft_neb_dispose;
  G_OBJECT_CLASS(klass)->finalize     = bigdft_neb_finalize;
  /* G_OBJECT_CLASS(klass)->set_property = visu_data_set_property; */
  /* G_OBJECT_CLASS(klass)->get_property = visu_data_get_property; */
}
#endif

static void bigdft_neb_init(BigDFT_Neb *obj)
{
#ifdef HAVE_GLIB
  memset((void*)((char*)obj + sizeof(GObject)), 0, sizeof(BigDFT_Neb) - sizeof(GObject));
#else
  memset(obj, 0, sizeof(BigDFT_Neb));
  G_OBJECT(obj)->ref_count = 1;
#endif
}
#ifdef HAVE_GLIB
static void bigdft_neb_dispose(GObject *obj)
{
  BigDFT_Neb *neb = BIGDFT_NEB(obj);

  if (neb->dispose_has_neb)
    return;
  neb->dispose_has_neb = TRUE;

  /* Chain up to the parent class */
  G_OBJECT_CLASS(bigdft_neb_parent_class)->dispose(obj);
}
#endif
static void bigdft_neb_finalize(GObject *obj)
{
  BigDFT_Neb *neb = BIGDFT_NEB(obj);

  if (neb->data)
    FC_FUNC_(neb_free, NEB_FREE)(&neb->data);

#ifdef HAVE_GLIB
  G_OBJECT_CLASS(bigdft_neb_parent_class)->finalize(obj);
#endif
}
BigDFT_Neb* bigdft_neb_new(BigDFT_Atoms *atoms, guint nimages, BigDFT_NebAlgo algo)
{
  BigDFT_Neb *neb;
  /* long self; */

#ifdef HAVE_GLIB
  neb = BIGDFT_NEB(g_object_new(BIGDFT_NEB_TYPE, NULL));
#else
  neb = g_malloc(sizeof(BigDFT_Neb));
  bigdft_neb_init(neb);
#endif
  /* self = *((long*)&neb); */
  FC_FUNC_(neb_new, NEB_NEW)(&neb->data, atoms->nat, (int*)&nimages, (int*)&algo);

  return neb;
}
/* void FC_FUNC_(neb_new_wrapper, NEB_NEW_WRAPPER)(double *self, void *obj) */
/* { */
/*   BigDFT_Neb *neb; */

/*   neb = bigdft_neb_new_from_fortran(obj); */
/*   *self = *((double*)&neb); */
/* } */
/* BigDFT_Neb* bigdft_neb_new_from_fortran(void *obj) */
/* { */
/*   BigDFT_Neb *neb; */

/* #ifdef HAVE_GLIB */
/*   neb = BIGDFT_NEB(g_object_new(BIGDFT_NEB_TYPE, NULL)); */
/* #else */
/*   neb = g_malloc(sizeof(BigDFT_Neb)); */
/*   bigdft_neb_init(neb); */
/* #endif */
/*   neb->data = obj; */

/*   return neb; */
/* } */
/* void FC_FUNC_(neb_free_wrapper, NEB_FREE_WRAPPER)(gpointer *obj) */
/* { */
/*   BigDFT_Neb *neb = BIGDFT_NEB(*obj); */

/*   neb->data = (gpointer)0; */
/*   bigdft_neb_free(neb); */
/* } */
void bigdft_neb_unref(BigDFT_Neb *neb)
{
  g_object_unref(G_OBJECT(neb));
#ifdef HAVE_GLIB
#else
  if (G_OBJECT(neb)->ref_count <= 0)
    {
      bigdft_neb_finalize(G_OBJECT(neb));
      g_free(neb);
    }
#endif
}
