#ifndef BIGDFT_SPRED_H
#define BIGDFT_SPRED_H

#include <bigdft.h>

G_BEGIN_DECLS

/*******************************/
/* BigDFT_Image data structure */
/*******************************/
typedef enum
  {
    BIGDFT_IMAGE_STEEPEST_DESCENT,
    BIGDFT_IMAGE_FLETCHER_REEVES,
    BIGDFT_IMAGE_POLAK_RIBIERE,
    BIGDFT_IMAGE_QUICK_MIN,
    BIGDFT_IMAGE_DAMPED_VERLET,
    BIGDFT_IMAGE_SIM_ANNEALING
  } BigDFT_ImageAlgo;

#ifdef GLIB_MAJOR_VERSION
#define BIGDFT_IMAGE_TYPE    (bigdft_image_get_type())
#define BIGDFT_IMAGE(obj)                                               \
  (G_TYPE_CHECK_INSTANCE_CAST(obj, BIGDFT_IMAGE_TYPE, BigDFT_Image))
typedef struct _BigDFT_ImageClass BigDFT_ImageClass;
struct _BigDFT_ImageClass
{
  GObjectClass parent;
};
GType bigdft_image_get_type(void);
#else
#define BIGDFT_IMAGE_TYPE    (999)
#define BIGDFT_IMAGE(obj)    ((BigDFT_Image*)obj)
#endif
typedef struct _BigDFT_Image BigDFT_Image;
struct _BigDFT_Image
{
  GObject parent;
  gboolean dispose_has_image;

  /* Bind attributes. */
  BigDFT_Run *run;
  BigDFT_Goutput *outs;
  double error, F;
  guint id;

  /* Private. */
  f90_run_image_pointer data;
};
/*********************************/
/* BigDFT_Image* bigdft_image_new       (BigDFT_Atoms *atoms, BigDFT_Inputs *ins, */
/*                                       BigDFT_Restart *rst, BigDFT_ImageAlgo algo); */
void          bigdft_image_unref     (BigDFT_Image *image);
void          bigdft_image_update_pos(BigDFT_Image *image, guint iteration,
                                      const BigDFT_Image *imgm1, const BigDFT_Image *imgp1,
                                      double k_before, double k_after,
                                      gboolean climbing);
gboolean      bigdft_image_update_pos_from_file(BigDFT_Image *image, guint iteration,
                                                const gchar *filem1, const gchar *filep1,
                                                double k_before, double k_after, gboolean climbing);
void          bigdft_image_calculate (BigDFT_Image *image, guint iteration, guint id);
BigDFT_Run*   bigdft_image_get_run(BigDFT_Image *image);
BigDFT_Goutput* bigdft_image_get_outs(BigDFT_Image *image);

GArray* bigdft_image_set_distribute(gboolean *update, guint nimages, guint ngroup);
/*********************************/

G_END_DECLS

#endif
