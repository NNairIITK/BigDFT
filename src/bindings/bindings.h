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
  FC_FUNC_(obj ## _get_ ## name, OBJ ## _GET_ ## NAME)(obj->data, &tmp); \
  obj->name = (guint*)tmp.data; \
  }
#define GET_ATTR_UINT_2D(obj,OBJ,name,NAME) {       \
  f90_pointer_int tmp; \
  memset(&tmp, 0, sizeof(f90_pointer_int_2D)); \
  FC_FUNC_(obj ## _get_ ## name, OBJ ## _GET_ ## NAME)(obj->data, &tmp); \
  obj->name = (guint*)tmp.data; \
  }
#define GET_ATTR_INT(obj,OBJ,name,NAME) {       \
  f90_pointer_int tmp; \
  memset(&tmp, 0, sizeof(f90_pointer_int)); \
  FC_FUNC_(obj ## _get_ ## name, OBJ ## _GET_ ## NAME)(obj->data, &tmp); \
  obj->name = (int*)tmp.data; \
  }
#define GET_ATTR_DBL(obj,OBJ,name,NAME) {       \
  f90_pointer_double tmp; \
  memset(&tmp, 0, sizeof(f90_pointer_double)); \
  FC_FUNC_(obj ## _get_ ## name, OBJ ## _GET_ ## NAME)(obj->data, &tmp); \
  obj->name = (double*)tmp.data; \
  }
#define GET_ATTR_DBL_2D(obj,OBJ,name,NAME) {       \
  f90_pointer_double_2D tmp; \
  memset(&tmp, 0, sizeof(f90_pointer_double_2D)); \
  FC_FUNC_(obj ## _get_ ## name, OBJ ## _GET_ ## NAME)(obj->data, &tmp); \
  obj->name = (double*)tmp.data; \
  }
#define GET_ATTR_DBL_3D(obj,OBJ,name,NAME) {       \
  f90_pointer_double_3D tmp; \
  memset(&tmp, 0, sizeof(f90_pointer_double_3D)); \
  FC_FUNC_(obj ## _get_ ## name, OBJ ## _GET_ ## NAME)(obj->data, &tmp); \
  obj->name = (double*)tmp.data; \
  }
#define GET_ATTR_DBL_4D(obj,OBJ,name,NAME) {       \
  f90_pointer_double_4D tmp; \
  memset(&tmp, 0, sizeof(f90_pointer_double_4D)); \
  FC_FUNC_(obj ## _get_ ## name, OBJ ## _GET_ ## NAME)(obj->data, &tmp); \
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

void FC_FUNC_(deallocate_double_1d, DEALLOCATE_DOUBLE_1D)(f90_pointer_double *array);
void FC_FUNC_(deallocate_double_2d, DEALLOCATE_DOUBLE_2D)(f90_pointer_double_2D *array);

/* Constructors of C wrappers around already built Fortran objects. */
BigDFT_Atoms*   bigdft_atoms_new_from_fortran  (_atoms_data *at, f90_pointer_double_2D *rxyz);
BigDFT_Inputs*  bigdft_inputs_new_from_fortran (_input_variables *inputs);
BigDFT_Energs*  bigdft_energs_new_from_fortran (_energy_terms *obj);
BigDFT_Restart* bigdft_restart_new_from_fortran(_restart_objects *obj);

#endif
