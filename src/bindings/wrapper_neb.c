#include <config.h>

#include "bigdft.h"
#include "bindings.h"
#include "bindings_api.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

/*******************************/
/* BigDFT_Image data structure */
/*******************************/
#ifdef HAVE_GLIB
G_DEFINE_TYPE(BigDFT_Image, bigdft_image, G_TYPE_OBJECT)

static void bigdft_image_dispose(GObject *image);
static void bigdft_image_finalize(GObject *image);

static void bigdft_image_class_init(BigDFT_ImageClass *klass)
{
  /* Connect the overloading methods. */
  G_OBJECT_CLASS(klass)->dispose      = bigdft_image_dispose;
  G_OBJECT_CLASS(klass)->finalize     = bigdft_image_finalize;
  /* G_OBJECT_CLASS(klass)->set_property = visu_data_set_property; */
  /* G_OBJECT_CLASS(klass)->get_property = visu_data_get_property; */
}
#endif

static void bigdft_image_init(BigDFT_Image *obj)
{
#ifdef HAVE_GLIB
  memset((void*)((char*)obj + sizeof(GObject)), 0, sizeof(BigDFT_Image) - sizeof(GObject));
#else
  memset(obj, 0, sizeof(BigDFT_Image));
  G_OBJECT(obj)->ref_count = 1;
#endif
}
#ifdef HAVE_GLIB
static void bigdft_image_dispose(GObject *obj)
{
  BigDFT_Image *image = BIGDFT_IMAGE(obj);

  if (image->dispose_has_image)
    return;
  image->dispose_has_image = TRUE;

  /* Chain up to the parent class */
  G_OBJECT_CLASS(bigdft_image_parent_class)->dispose(obj);
}
#endif
static void bigdft_image_finalize(GObject *obj)
{
  BigDFT_Image *image = BIGDFT_IMAGE(obj);
  f90_run_objects_pointer run;
  f90_DFT_global_output_pointer outs;

  if (F_TYPE(image->data))
    {
      FC_FUNC_(image_free, IMAGE_FREE)(&image->data, &run, &outs);
      image->run->data  = run;
      image->outs->data = outs;
    }
  bigdft_run_unref(image->run);
  bigdft_goutput_unref(image->outs);

#ifdef HAVE_GLIB
  G_OBJECT_CLASS(bigdft_image_parent_class)->finalize(obj);
#endif
}
BigDFT_Image* bigdft_image_new(BigDFT_Atoms *atoms, BigDFT_Inputs *ins, BigDFT_Restart *rst, BigDFT_ImageAlgo algo)
{
  BigDFT_Image *image;
  f90_run_objects_pointer run;
  f90_DFT_global_output_pointer outs;
  int algo_;
  /* long self; */

#ifdef HAVE_GLIB
  image = BIGDFT_IMAGE(g_object_new(BIGDFT_IMAGE_TYPE, NULL));
#else
  image = g_malloc(sizeof(BigDFT_Image));
  bigdft_image_init(image);
#endif
  /* self = *((long*)&image); */
  algo_ = algo + 1;
  FC_FUNC_(image_new, IMAGE_NEW)(&image->data, &run, &outs,
                                 F_TYPE(atoms->data), F_TYPE(ins->data),
                                 F_TYPE(rst->data), (int*)&algo_);
  image->run = bigdft_run_new_from_fortran(run, FALSE);
  image->run->atoms = atoms;
  g_object_ref(G_OBJECT(atoms));
  image->run->inputs = ins;
  bigdft_inputs_ref(ins);
  image->run->restart = rst;
  g_object_ref(G_OBJECT(rst));
  image->outs = bigdft_goutput_new_from_fortran(outs);
  FC_FUNC_(image_get_attributes, IMAGE_GET_ATTRIBUTES)(F_TYPE(image->data), &image->error, &image->F,
                                                       (int*)&image->id);

  return image;
}
/* void FC_FUNC_(image_new_wrapper, IMAGE_NEW_WRAPPER)(double *self, void *obj) */
/* { */
/*   BigDFT_Image *image; */

/*   image = bigdft_image_new_from_fortran(obj); */
/*   *self = *((double*)&image); */
/* } */
/* BigDFT_Image* bigdft_image_new_from_fortran(void *obj) */
/* { */
/*   BigDFT_Image *image; */

/* #ifdef HAVE_GLIB */
/*   image = BIGDFT_IMAGE(g_object_new(BIGDFT_IMAGE_TYPE, NULL)); */
/* #else */
/*   image = g_malloc(sizeof(BigDFT_Image)); */
/*   bigdft_image_init(image); */
/* #endif */
/*   image->data = obj; */

/*   return image; */
/* } */
/* void FC_FUNC_(image_free_wrapper, IMAGE_FREE_WRAPPER)(gpointer *obj) */
/* { */
/*   BigDFT_Image *image = BIGDFT_IMAGE(*obj); */

/*   image->data = (gpointer)0; */
/*   bigdft_image_free(image); */
/* } */
void bigdft_image_unref(BigDFT_Image *image)
{
  g_object_unref(G_OBJECT(image));
#ifdef HAVE_GLIB
#else
  if (G_OBJECT(image)->ref_count <= 0)
    {
      bigdft_image_finalize(G_OBJECT(image));
      g_free(image);
    }
#endif
}
/**
 * bigdft_image_update_pos:
 * @image: 
 * @iteration: 
 * @imgm1: (allow-none):
 * @imgp1: (allow-none):
 * @k_before: 
 * @k_after: 
 * @climbing: 
 *
 * 
 **/
void bigdft_image_update_pos(BigDFT_Image *image, guint iteration,
                             const BigDFT_Image *imgm1, const BigDFT_Image *imgp1,
                             double k_before, double k_after, gboolean climbing)
{
  const BigDFT_Image *m1, *p1;
  gboolean optimization;
  f90_NEB_data_pointer neb;

  m1 = (imgm1)?imgm1:image;
  p1 = (imgp1)?imgp1:image;
  optimization = (m1 == image) || (p1 == image);
  FC_FUNC_(neb_new, NEB_NEW)(&neb);
  FC_FUNC_(image_update_pos, IMAGE_UPDATE_POS)(F_TYPE(image->data), (int*)&iteration,
                                               m1->run->atoms->rxyz, p1->run->atoms->rxyz,
                                               &m1->outs->etot, &p1->outs->etot,
                                               &k_before, &k_after,
                                               &optimization, &climbing, F_TYPE(neb));
  FC_FUNC_(neb_free, NEB_FREE)(&neb);
  FC_FUNC_(image_get_attributes, IMAGE_GET_ATTRIBUTES)(F_TYPE(image->data), &image->error, &image->F,
                                                       (int*)&image->id);
}
/**
 * bigdft_image_update_pos_from_file:
 * @image: 
 * @iteration: 
 * @filem1: (allow-none):
 * @filep1: (allow-none):
 * @k_before: 
 * @k_after: 
 * @climbing: 
 *
 * 
 *
 * Returns: 
 **/
gboolean bigdft_image_update_pos_from_file(BigDFT_Image *image, guint iteration,
                                           const gchar *filem1, const gchar *filep1,
                                           double k_before, double k_after, gboolean climbing)
{
  f90_NEB_data_pointer neb;
  gchar *filem1_, *filep1_;

  filem1_ = g_strdup((filem1)?filem1:" ");
  filep1_ = g_strdup((filep1)?filep1:" ");
  FC_FUNC_(neb_new, NEB_NEW)(&neb);
  FC_FUNC_(image_update_pos_from_file, IMAGE_UPDATE_POS_FROM_FILE)
    (F_TYPE(image->data), (int*)&iteration, filem1_, filep1_, &k_before, &k_after, &climbing, F_TYPE(neb),
     strlen(filem1_), strlen(filep1_));
  FC_FUNC_(neb_free, NEB_FREE)(&neb);
  FC_FUNC_(image_get_attributes, IMAGE_GET_ATTRIBUTES)(F_TYPE(image->data), &image->error, &image->F,
                                                       (int*)&image->id);

  return (image->error >= 0.);
}
void bigdft_image_calculate(BigDFT_Image *image, guint iteration, guint id)
{
  FC_FUNC_(image_calculate, IMAGE_CALCULATE)(F_TYPE(image->data), (int*)&iteration, (int*)&id);
  FC_FUNC_(image_get_attributes, IMAGE_GET_ATTRIBUTES)(F_TYPE(image->data), &image->error, &image->F,
                                                       (int*)&image->id);
}

/**
 * bigdft_image_set_distribute:
 * @update: (array length=nimages):
 * @nimages: 
 * @ngroup: 
 *
 *
 * Returns: (transfer full) (element-type int):
 **/
GArray* bigdft_image_set_distribute(gboolean *update, guint nimages, guint ngroup)
{
  GArray *igroup;

  igroup = g_array_sized_new(FALSE, FALSE, sizeof(int), nimages);
  g_array_set_size(igroup, nimages);
  FC_FUNC_(images_distribute_tasks, IMAGES_DISTRIBUTE_TASKS)((int*)igroup->data, (int*)update,
                                                             (int*)&nimages, (int*)&ngroup);
  return igroup;
}

/**
 * bigdft_image_get_run:
 * @image: 
 *
 * 
 *
 * Returns: (transfer full):
 **/
BigDFT_Run* bigdft_image_get_run(BigDFT_Image *image)
{
  g_object_ref(G_OBJECT(image->run));
  return image->run;
}
/**
 * bigdft_image_get_outs:
 * @image: 
 *
 * 
 *
 * Returns: (transfer full):
 **/
BigDFT_Goutput* bigdft_image_get_outs(BigDFT_Image *image)
{
  g_object_ref(G_OBJECT(image->outs));
  return image->outs;
}
