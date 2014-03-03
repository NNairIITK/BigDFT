#ifndef BINDINGS_H
#define BINDINGS_H

#ifdef HAVE_DEBUG
#define DBG_MEM(A, T) {int i__; for (i__ = 0; i__ < sizeof(T) / sizeof(void*); i__++) \
                                  fprintf(stderr, "DBG (%2d) -> %p\n", i__, ((void**)A)[i__]); }
#else
#define DBG_MEM(A, T)
#endif

#define GET_ATTR_UINT(obj,OBJ,name,NAME) {       \
  f90_pointer_int tmp; \
  memset(&tmp, 0, sizeof(f90_pointer_int)); \
  FC_FUNC_(obj ## _get_ ## name, OBJ ## _GET_ ## NAME)(F_TYPE(obj->data), &tmp); \
  obj->name = (guint*)tmp.data; \
  }
#define GET_ATTR_UINT_2D(obj,OBJ,name,NAME) {       \
  f90_pointer_int tmp; \
  memset(&tmp, 0, sizeof(f90_pointer_int_2D)); \
  FC_FUNC_(obj ## _get_ ## name, OBJ ## _GET_ ## NAME)(F_TYPE(obj->data), &tmp); \
  obj->name = (guint*)tmp.data; \
  }
#define GET_ATTR_INT(obj,OBJ,name,NAME) {       \
  f90_pointer_int tmp; \
  memset(&tmp, 0, sizeof(f90_pointer_int)); \
  FC_FUNC_(obj ## _get_ ## name, OBJ ## _GET_ ## NAME)(F_TYPE(obj->data), &tmp); \
  obj->name = (int*)tmp.data; \
  }
#define GET_ATTR_DBL(obj,OBJ,name,NAME) {       \
  f90_pointer_double tmp; \
  memset(&tmp, 0, sizeof(f90_pointer_double)); \
  FC_FUNC_(obj ## _get_ ## name, OBJ ## _GET_ ## NAME)(F_TYPE(obj->data), &tmp); \
  obj->name = (double*)tmp.data; \
  }
#define GET_ATTR_DBL_2D(obj,OBJ,name,NAME) {       \
  f90_pointer_double_2D tmp; \
  memset(&tmp, 0, sizeof(f90_pointer_double_2D)); \
  FC_FUNC_(obj ## _get_ ## name, OBJ ## _GET_ ## NAME)(F_TYPE(obj->data), &tmp); \
  obj->name = (double*)tmp.data; \
  }
#define GET_ATTR_DBL_3D(obj,OBJ,name,NAME) {       \
  f90_pointer_double_3D tmp; \
  memset(&tmp, 0, sizeof(f90_pointer_double_3D)); \
  FC_FUNC_(obj ## _get_ ## name, OBJ ## _GET_ ## NAME)(F_TYPE(obj->data), &tmp); \
  obj->name = (double*)tmp.data; \
  }
#define GET_ATTR_DBL_4D(obj,OBJ,name,NAME) {       \
  f90_pointer_double_4D tmp; \
  memset(&tmp, 0, sizeof(f90_pointer_double_4D)); \
  FC_FUNC_(obj ## _get_ ## name, OBJ ## _GET_ ## NAME)(F_TYPE(obj->data), &tmp); \
  obj->name = (double*)tmp.data; \
  }

void FC_FUNC_(f90_pointer_1d_init, F90_POINTER_1D_INIT)(f90_pointer_double *pt, guint *s);
void FC_FUNC_(f90_pointer_2d_init, F90_POINTER_2D_INIT)(f90_pointer_double_2D *pt, guint *s);
void FC_FUNC_(f90_pointer_3d_init, F90_POINTER_3D_INIT)(f90_pointer_double_3D *pt, guint *s);
void FC_FUNC_(f90_pointer_4d_init, F90_POINTER_4D_INIT)(f90_pointer_double_4D *pt, guint *s);
void FC_FUNC_(f90_pointer_5d_init, F90_POINTER_5D_INIT)(f90_pointer_double_5D *pt, guint *s);

#define F90_1D_POINTER_INIT(pt) {                                       \
  guint size_ = F90_1D_POINTER_SIZE;                                    \
  FC_FUNC_(f90_pointer_1d_init, F90_POINTER_1D_INIT)(pt, &size_);       \
  }
#define F90_2D_POINTER_INIT(pt) {                                       \
  guint size_ = F90_2D_POINTER_SIZE;                                    \
  FC_FUNC_(f90_pointer_2d_init, F90_POINTER_2D_INIT)(pt, &size_);       \
  }
#define F90_3D_POINTER_INIT(pt) {                                       \
  guint size_ = F90_3D_POINTER_SIZE;                                    \
  FC_FUNC_(f90_pointer_3d_init, F90_POINTER_3D_INIT)(pt, &size_);       \
  }
#define F90_4D_POINTER_INIT(pt) {                                       \
  guint size_ = F90_4D_POINTER_SIZE;                                    \
  FC_FUNC_(f90_pointer_4d_init, F90_POINTER_4D_INIT)(pt, &size_);       \
  }
#define F90_5D_POINTER_INIT(pt) {                                       \
  guint size_ = F90_5D_POINTER_SIZE;                                    \
  FC_FUNC_(f90_pointer_5d_init, F90_POINTER_5D_INIT)(pt, &size_);       \
  }

/* Constructors of C wrappers around already built Fortran objects. */
BigDFT_Atoms*   bigdft_atoms_new_from_fortran  (f90_atoms_data_pointer at);
BigDFT_Inputs*  bigdft_inputs_new_from_fortran (f90_input_variables_pointer inputs);
BigDFT_Restart* bigdft_restart_new_from_fortran(f90_restart_objects_pointer obj);
BigDFT_Run*     bigdft_run_new_from_fortran    (f90_run_objects_pointer obj,
                                                gboolean create_wrappers);
BigDFT_Goutput* bigdft_goutput_new_from_fortran(f90_DFT_global_output_pointer obj);

/* Additional private methods. */
void _inputs_sync(BigDFT_Inputs *in);
void _inputs_sync_add(BigDFT_Inputs *in);

/*  Generic tools. */
gchar* _get_c_string(const gchar *fstr, guint len);

#endif
