/*
 * Copyright 1993-2007 NVIDIA Corporation.  All rights reserved.
 *
 * NOTICE TO USER:   
 *
 * This source code is subject to NVIDIA ownership rights under U.S. and 
 * international Copyright laws.  Users and possessors of this source code 
 * are hereby granted a nonexclusive, royalty-free license to use this code 
 * in individual and commercial software.
 *
 * NVIDIA MAKES NO REPRESENTATION ABOUT THE SUITABILITY OF THIS SOURCE 
 * CODE FOR ANY PURPOSE.  IT IS PROVIDED "AS IS" WITHOUT EXPRESS OR 
 * IMPLIED WARRANTY OF ANY KIND.  NVIDIA DISCLAIMS ALL WARRANTIES WITH 
 * REGARD TO THIS SOURCE CODE, INCLUDING ALL IMPLIED WARRANTIES OF 
 * MERCHANTABILITY, NONINFRINGEMENT, AND FITNESS FOR A PARTICULAR PURPOSE.
 * IN NO EVENT SHALL NVIDIA BE LIABLE FOR ANY SPECIAL, INDIRECT, INCIDENTAL, 
 * OR CONSEQUENTIAL DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS 
 * OF USE, DATA OR PROFITS,  WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE 
 * OR OTHER TORTIOUS ACTION,  ARISING OUT OF OR IN CONNECTION WITH THE USE 
 * OR PERFORMANCE OF THIS SOURCE CODE.  
 *
 * U.S. Government End Users.   This source code is a "commercial item" as 
 * that term is defined at  48 C.F.R. 2.101 (OCT 1995), consisting  of 
 * "commercial computer  software"  and "commercial computer software 
 * documentation" as such terms are  used in 48 C.F.R. 12.212 (SEPT 1995) 
 * and is provided to the U.S. Government only as a commercial end item.  
 * Consistent with 48 C.F.R.12.212 and 48 C.F.R. 227.7202-1 through 
 * 227.7202-4 (JUNE 1995), all U.S. Government End Users acquire the 
 * source code with only those rights set forth herein. 
 *
 * Any use of this source code in individual and commercial software must 
 * include, in the user documentation and internal comments to the code,
 * the above Disclaimer and U.S. Government End Users Notice.
 */

/*
 * This file contains example Fortran bindings for the CUBLAS library, These
 * bindings have been tested with Intel Fortran 9.0 on 32-bit and 64-bit 
 * Windows, and with g77 3.4.5 on 32-bit and 64-bit Linux. They will likely
 * have to be adjusted for other Fortran compilers and platforms.
 */

#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#if defined(__GNUC__)
#include <stdint.h>
#endif /* __GNUC__ */
#include "cublas.h"   /* CUBLAS public header file  */

#define imin(a,b) (((a)<(b))?(a):(b))
#define imax(a,b) (((a)<(b))?(b):(a))

#define CUBLAS_G77              1
#define CUBLAS_INTEL_FORTRAN    2

/* Default to g77 on Linux, and Intel Fortran on Win32 */
#if defined(_WIN32)
#define CUBLAS_FORTRAN_COMPILER CUBLAS_INTEL_FORTRAN
#elif defined(__linux)
#define CUBLAS_FORTRAN_COMPILER CUBLAS_G77
#else
#error unsupported platform
#endif

#if CUBLAS_FORTRAN_COMPILER==CUBLAS_G77
/* NOTE: Must use -fno-second-underscore when building Fortran source with g77
 *       g77 invocation may not use -fno-f2c, which forces different return 
 *       type conventions than the one used below
 */
#define CUBLAS_INIT             cublas_init_
#define CUBLAS_SHUTDOWN         cublas_shutdown_
#define CUBLAS_ALLOC            cublas_alloc_
#define CUBLAS_FREE             cublas_free_
#define CUBLAS_SET_VECTOR       cublas_set_vector_
#define CUBLAS_GET_VECTOR       cublas_get_vector_
#define CUBLAS_SET_MATRIX       cublas_set_matrix_
#define CUBLAS_GET_MATRIX       cublas_get_matrix_
#define CUBLAS_GET_ERROR        cublas_get_error_
#define CUBLAS_XERBLA           cublas_xerbla_
#define CUBLAS_ISAMAX           cublas_isamax_
#define CUBLAS_ISAMIN           cublas_isamin_
#define CUBLAS_SASUM            cublas_sasum_
#define CUBLAS_SAXPY            cublas_saxpy_
#define CUBLAS_SCOPY            cublas_scopy_
#define CUBLAS_SDOT             cublas_sdot_
#define CUBLAS_SNRM2            cublas_snrm2_
#define CUBLAS_SROT             cublas_srot_
#define CUBLAS_SROTG            cublas_srotg_
#define CUBLAS_SROTM            cublas_srotm_
#define CUBLAS_SROTMG           cublas_srotmg_
#define CUBLAS_SSCAL            cublas_sscal_
#define CUBLAS_SSWAP            cublas_sswap_
#define CUBLAS_CAXPY            cublas_caxpy_
#define CUBLAS_CCOPY            cublas_ccopy_
#define CUBLAS_CROT             cublas_crot_
#define CUBLAS_CROTG            cublas_crotg_
#define CUBLAS_CSCAL            cublas_cscal_
#define CUBLAS_CSROT            cublas_csrot_
#define CUBLAS_CSSCAL           cublas_csscal_
#define CUBLAS_CSWAP            cublas_cswap_
#define CUBLAS_CDOTU            cublas_cdotu_
#define CUBLAS_CDOTC            cublas_cdotc_
#define CUBLAS_ICAMAX           cublas_icamax_
#define CUBLAS_SCASUM           cublas_scasum_
#define CUBLAS_SCNRM2           cublas_scnrm2_
#define CUBLAS_SGBMV            cublas_sgbmv_
#define CUBLAS_SGEMV            cublas_sgemv_
#define CUBLAS_SGER             cublas_sger_
#define CUBLAS_SSBMV            cublas_ssbmv_
#define CUBLAS_SSPMV            cublas_sspmv_
#define CUBLAS_SSPR             cublas_sspr_
#define CUBLAS_SSPR2            cublas_sspr2_
#define CUBLAS_SSYMV            cublas_ssymv_
#define CUBLAS_SSYR             cublas_ssyr_
#define CUBLAS_SSYR2            cublas_ssyr2_
#define CUBLAS_STBMV            cublas_stbmv_
#define CUBLAS_STBSV            cublas_stbsv_
#define CUBLAS_STPMV            cublas_stpmv_
#define CUBLAS_STPSV            cublas_stpsv_
#define CUBLAS_STRMV            cublas_strmv_
#define CUBLAS_STRSV            cublas_strsv_
#define CUBLAS_SGEMM            cublas_sgemm_
#define CUBLAS_SSYMM            cublas_ssymm_
#define CUBLAS_SSYR2K           cublas_ssyr2k_
#define CUBLAS_SSYRK            cublas_ssyrk_
#define CUBLAS_STRMM            cublas_strmm_
#define CUBLAS_STRSM            cublas_strsm_
#define CUBLAS_CGEMM            cublas_cgemm_
#define CUBLAS_CHEMM            cublas_chemm_
#define CUBLAS_CSYMM            cublas_csymm_
#define CUBLAS_CTRMM            cublas_ctrmm_
#define CUBLAS_CTRSM            cublas_ctrsm_
#define CUBLAS_CHERK            cublas_cherk_
#define CUBLAS_CSYRK            cublas_csyrk_
#define CUBLAS_CHER2K           cublas_cher2k_
#define CUBLAS_CSYR2K           cublas_csyr2k_
#elif CUBLAS_FORTRAN_COMPILER==CUBLAS_INTEL_FORTRAN
#define CUBLAS_INIT             CUBLAS_INIT 
#define CUBLAS_SHUTDOWN         CUBLAS_SHUTDOWN
#define CUBLAS_ALLOC            CUBLAS_ALLOC
#define CUBLAS_FREE             CUBLAS_FREE
#define CUBLAS_SET_VECTOR       CUBLAS_SET_VECTOR
#define CUBLAS_GET_VECTOR       CUBLAS_GET_VECTOR
#define CUBLAS_SET_MATRIX       CUBLAS_SET_MATRIX
#define CUBLAS_GET_MATRIX       CUBLAS_GET_MATRIX
#define CUBLAS_GET_ERROR        CUBLAS_GET_ERROR
#define CUBLAS_XERBLA           CUBLAS_XERBLA
#define CUBLAS_ISAMAX           CUBLAS_ISAMAX
#define CUBLAS_ISAMIN           CUBLAS_ISAMIN
#define CUBLAS_SASUM            CUBLAS_SASUM
#define CUBLAS_SAXPY            CUBLAS_SAXPY
#define CUBLAS_SCOPY            CUBLAS_SCOPY
#define CUBLAS_SDOT             CUBLAS_SDOT
#define CUBLAS_SNRM2            CUBLAS_SNRM2
#define CUBLAS_SROT             CUBLAS_SROT
#define CUBLAS_SROTG            CUBLAS_SROTG
#define CUBLAS_SROTM            CUBLAS_SROTM
#define CUBLAS_SROTMG           CUBLAS_SROTMG
#define CUBLAS_SSCAL            CUBLAS_SSCAL
#define CUBLAS_SSWAP            CUBLAS_SSWAP
#define CUBLAS_CAXPY            CUBLAS_CAXPY
#define CUBLAS_CCOPY            CUBLAS_CCOPY
#define CUBLAS_CROT             CUBLAS_CROT
#define CUBLAS_CROTG            CUBLAS_CROTG
#define CUBLAS_CSCAL            CUBLAS_CSCAL
#define CUBLAS_CSROT            CUBLAS_CSROT
#define CUBLAS_CSSCAL           CUBLAS_CSSCAL
#define CUBLAS_CSWAP            CUBLAS_CSWAP 
#define CUBLAS_CDOTU            CUBLAS_CDOTU
#define CUBLAS_CDOTC            CUBLAS_CDOTC
#define CUBLAS_ICAMAX           CUBLAS_ICAMAX
#define CUBLAS_SCASUM           CUBLAS_SCASUM
#define CUBLAS_SCNRM2           CUBLAS_SCNRM2
#define CUBLAS_SGBMV            CUBLAS_SGBMV
#define CUBLAS_SGEMV            CUBLAS_SGEMV
#define CUBLAS_SGER             CUBLAS_SGER
#define CUBLAS_SSBMV            CUBLAS_SSBMV
#define CUBLAS_SSPMV            CUBLAS_SSPMV
#define CUBLAS_SSPR             CUBLAS_SSPR
#define CUBLAS_SSPR2            CUBLAS_SSPR2
#define CUBLAS_SSYMV            CUBLAS_SSYMV
#define CUBLAS_SSYR             CUBLAS_SSYR
#define CUBLAS_SSYR2            CUBLAS_SSYR2
#define CUBLAS_STBMV            CUBLAS_STBMV
#define CUBLAS_STBSV            CUBLAS_STBSV
#define CUBLAS_STPMV            CUBLAS_STPMV
#define CUBLAS_STPSV            CUBLAS_STPSV
#define CUBLAS_STRMV            CUBLAS_STRMV
#define CUBLAS_STRSV            CUBLAS_STRSV
#define CUBLAS_SGEMM            CUBLAS_SGEMM
#define CUBLAS_SSYMM            CUBLAS_SSYMM
#define CUBLAS_SSYR2K           CUBLAS_SSYR2K
#define CUBLAS_SSYRK            CUBLAS_SSYRK
#define CUBLAS_STRMM            CUBLAS_STRMM
#define CUBLAS_STRSM            CUBLAS_STRSM
#define CUBLAS_CGEMM            CUBLAS_CGEMM
#define CUBLAS_CHEMM            CUBLAS_CHEMM
#define CUBLAS_CSYMM            CUBLAS_CSYMM
#define CUBLAS_CTRMM            CUBLAS_CTRMM
#define CUBLAS_CTRSM            CUBLAS_CTRSM
#define CUBLAS_CHERK            CUBLAS_CHERK
#define CUBLAS_CSYRK            CUBLAS_CSYRK
#define CUBLAS_CHER2K           CUBLAS_CHER2K
#define CUBLAS_CSYR2K           CUBLAS_CSYR2K
#else
#error unsupported Fortran compiler
#endif

/* For now, the GPU only supports a 32-bit address space, so device pointers
   can be represented as INTEGER*4 in Fortran. In the future, device pointers
   may become 64-bit pointers, and will have to be represented as INTEGER*8 in
   Fortran, at which point devptr_t needs to be typedef'ed as long long.
*/
typedef int devptr_t;

#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus */
int CUBLAS_INIT (void);
int CUBLAS_SHUTDOWN (void);
int CUBLAS_ALLOC (const int *n, const int *elemSize, devptr_t *devicePtr);
int CUBLAS_FREE (const devptr_t *devicePtr);
int CUBLAS_SET_VECTOR (const int *n, const int *elemSize, const void *x,
                       const int *incx, const devptr_t *y, const int *incy);
int CUBLAS_GET_VECTOR (const int *n, const int *elemSize, const devptr_t *x,
                       const int *incx, void *y, const int *incy);
int CUBLAS_SET_MATRIX (const int *rows, const int *cols, const int *elemSize,
                       const void *A, const int *lda, const int *B, 
                       const int *ldb);
int CUBLAS_GET_MATRIX (const int *rows, const int *cols, const int *elemSize,
                       const int *A, const int *lda, void *B, const int *ldb);

/* BLAS util */
void CUBLAS_XERBLA (const char *srName, int *info);
#if defined(__cplusplus)
}
#endif /* __cplusplus */

int CUBLAS_INIT (void) 
{
    return (int)cublasInit ();
}

int CUBLAS_SHUTDOWN (void) 
{
    return (int)cublasShutdown ();
}

int CUBLAS_ALLOC (const int *n, const int *elemSize, devptr_t *devicePtr)
{    
    void *tPtr;
    int retVal;
    retVal = (int)cublasAlloc (*n, *elemSize, &tPtr);
    *devicePtr = (devptr_t)(uintptr_t)tPtr;
    return retVal;
}

int CUBLAS_FREE (const devptr_t *devicePtr)
{
    void *tPtr;
    tPtr = (void *)(uintptr_t)(*devicePtr);
    return (int)cublasFree (tPtr);
}

int CUBLAS_SET_VECTOR (const int *n, const int *elemSize, const void *x,
                       const int *incx, const devptr_t *y, const int *incy)
{
    void *tPtr = (void *)(uintptr_t)(*y);
    return (int)cublasSetVector (*n, *elemSize, x, *incx, tPtr, *incy);
}

int CUBLAS_GET_VECTOR (const int *n, const int *elemSize, const devptr_t *x,
                       const int *incx, void *y, const int *incy)
{
    const void *tPtr = (const void *)(uintptr_t)(*x);
    return (int)cublasGetVector (*n, *elemSize, tPtr, *incx, y, *incy);
}

int CUBLAS_SET_MATRIX (const int *rows, const int *cols, const int *elemSize,
                       const void *A, const int *lda, const devptr_t *B, 
                       const int *ldb)
{
    void *tPtr = (void *)(uintptr_t)(*B);
    return (int)cublasSetMatrix (*rows, *cols, *elemSize, A, *lda, tPtr,*ldb);
}

int CUBLAS_GET_MATRIX (const int *rows, const int *cols, const int *elemSize,
                       const devptr_t *A, const int *lda, void *B, 
                       const int *ldb)
{
    const void *tPtr = (const void *)(uintptr_t)(*A);
    return (int)cublasGetMatrix (*rows, *cols, *elemSize, tPtr, *lda, B, *ldb);
}

int CUBLAS_GET_ERROR (void)
{
    return (int)cublasGetError();
}

void CUBLAS_XERBLA (const char *srName, int *info)
{
    cublasXerbla (srName, *info);
}

/*
 *  Fortran callable BLAS functions that include GPU memory allocation and
 *  copy-up and copy-down code. These can be called from unmodified Fortran 
 *  code, but they are inefficient due to the data constantly bouncing back 
 *  and forth between CPU and GPU.
 */
#if defined(CUBLAS_USE_THUNKING)

#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus */
/* BLAS1 */
#if CUBLAS_FORTRAN_COMPILER==CUBLAS_G77
double CUBLAS_SDOT (const int *n, const float *x, const int *incx, float *y, 
                    const int *incy);
double CUBLAS_SASUM (const int *n, const float *x, const int *incx);
double CUBLAS_SNRM2 (const int *n, const float *x, const int *incx);
double CUBLAS_SCASUM (const int *n, const cuComplex *x, const int *incx);
double CUBLAS_SCNRM2 (const int *n, const cuComplex *x, const int *incx);
#else
float CUBLAS_SDOT (const int *n, const float *x, const int *incx, float *y, 
                   const int *incy);
float CUBLAS_SASUM (const int *n, const float *x, const int *incx);
float CUBLAS_SNRM2 (const int *n, const float *x, const int *incx);
float CUBLAS_SCASUM (const int *n, const cuComplex *x, const int *incx);
float CUBLAS_SCNRM2 (const int *n, const cuComplex *x, const int *incx);
#endif
int CUBLAS_ISAMAX (const int *n, const float *x, const int *incx);
int CUBLAS_ISAMIN (const int *n, const float *x, const int *incx);
void CUBLAS_SAXPY (const int *n, const float *alpha, const float *x, 
                   const int *incx, float *y, const int *incy);
void CUBLAS_SCOPY (const int *n, const float *x, const int *incx, float *y, 
                   const int *incy);
void CUBLAS_SROT (const int *n, float *x, const int *incx, float *y, 
                  const int *incy, const float *sc, const float *ss);
void CUBLAS_SROTG (float *sa, float *sb, float *sc, float *ss);
void CUBLAS_SROTM (const int *n, float *x, const int *incx, float *y,
                   const int *incy, const float* sparam);
void CUBLAS_SROTMG (float *sd1, float *sd2, float *sx1, const float *sy1, 
                    float* sparam);
void CUBLAS_SSCAL (const int *n, const float *alpha, float *x,
                   const int *incx);
void CUBLAS_SSWAP(const int *n, float *x, const int *incx, float *y,
                  const int *incy);

void CUBLAS_CAXPY (const int *n, const cuComplex *alpha, const cuComplex *x, 
                   const int *incx, cuComplex *y, const int *incy);
void CUBLAS_CCOPY (const int *n, const cuComplex *x, const int *incx, 
                   cuComplex *y, const int *incy);
void CUBLAS_CROT (const int *n, cuComplex *x, const int *incx, cuComplex *y, 
                  const int *incy, const float *sc, const cuComplex *cs);
void CUBLAS_CROTG (cuComplex *ca, const cuComplex *cb, float *sc, 
                   cuComplex *cs);
void CUBLAS_CSCAL (const int *n, const cuComplex *alpha, cuComplex *x, 
                   const int *incx);
void CUBLAS_CSROT (const int *n, cuComplex *x, const int *incx, cuComplex *y, 
                   const int *incy, const float *sc, const float *ss);
void CUBLAS_CSSCAL (const int *n, const float *alpha, cuComplex *x,
                    const int *incx);
void CUBLAS_CSWAP (const int *n, cuComplex *x, const int *incx, cuComplex *y,
                   const int *incy);
void CUBLAS_CDOTU (cuComplex *retVal, const int *n, const cuComplex *x, 
                   const int *incx, const cuComplex *y, const int *incy);
void CUBLAS_CDOTC (cuComplex *retVal,const int *n, const cuComplex *x, 
                   const int *incx, const cuComplex *y, const int *incy);
int CUBLAS_ICAMAX (const int *n, const cuComplex *x, const int *incx);
int CUBLAS_ICAMIN (const int *n, const cuComplex *x, const int *incx);

/* BLAS2 */
void CUBLAS_SGBMV (const char *trans, const int *m, const int *n, 
                   const int *kl, const int *ku, const float *alpha, 
                   const float *A, const int *lda, const float *x,
                   const int *incx, const float *beta, float *y, 
                   const int *incy);
void CUBLAS_SGEMV (const char *trans, const int *m, const int *n,
                   const float *alpha, const float *A, const int *lda,
                   const float *x, const int *incx, const float *beta, 
                   float *y, const int *incy);
void CUBLAS_SGER (const int *m, const int *n, const float *alpha,
                  const float *x, const int *incx, const float *y,
                  const int *incy, float *A, const int *lda);
void CUBLAS_SSBMV (const char *uplo, const int *n, const int *k, 
                   const float *alpha, const float *A, const int *lda,
                   const float *x, const int *incx, const float *beta,
                   float *y, const int *incy);
void CUBLAS_SSPMV (const char *uplo, const int *n, const float *alpha, 
                   const float *AP, const float *x, const int *incx, 
                   const float *beta, float *y, const int *incy);
void CUBLAS_SSPR (const char *uplo, const int *n, const float *alpha,
                  const float *x, const int *incx, float *AP);
void CUBLAS_SSPR2 (const char *uplo, const int *n, const float *alpha,
                   const float *x, const int *incx, const float *y,
                   const int *incy, float *AP);
void CUBLAS_SSYMV (const char *uplo, const int *n, const float *alpha,
                   const float *A, const int *lda, const float *x,
                   const int *incx, const float *beta, float *y, 
                   const int *incy);
void CUBLAS_SSYR (const char *uplo, const int *n, const float *alpha,
                  const float *x, const int *incx, float *A, const int *lda);
void CUBLAS_SSYR2 (const char *uplo, const int *n, const float *alpha,
                   const float *x, const int *incx, const float *y,
                   const int *incy, float *A, const int *lda);
void CUBLAS_STBMV (const char *uplo, const char *trans, const char *diag, 
                   const int *n, const int *k, const float *A, const int *lda, 
                   float *x, const int *incx);
void CUBLAS_STBSV (const char *uplo, const char *trans, const char *diag, 
                   const int *n, const int *k, const float *A, const int *lda, 
                   float *x, const int *incx);
void CUBLAS_STPMV (const char *uplo, const char *trans, const char *diag, 
                   const int *n, const float *AP, float *x, const int *incx);
void CUBLAS_STPSV (const char *uplo, const char *trans, const char *diag, 
                   const int *n, const float *AP, float *x, const int *incx);
void CUBLAS_STRMV (const char *uplo, const char *trans, const char *diag, 
                   const int *n, const float *A, const int *lda, float *x, 
                   const int *incx);
void CUBLAS_STRSV (const char *uplo, const char *trans, const char *diag, 
                   const int *n, const float *A, const int *lda, float *x, 
                   const int *incx);

/* BLAS3 */
void CUBLAS_SGEMM (const char *transa, const char *transb, const int *m,
                   const int *n, const int *k, const float *alpha, 
                   const float *A, const int *lda, const float *B, 
                   const int *ldb, const float *beta, float *C, 
                   const int *ldc);
void CUBLAS_SSYMM (const char *side, const char *uplo, const int *m,
                   const int *n, const float *alpha, const float *A,
                   const int *lda, const float *B, const int *ldb, 
                   const float *beta, float *C, const int *ldc);
void CUBLAS_SSYR2K (const char *uplo, const char *trans, const int *n, 
                    const int *k, const float *alpha, const float *A,
                    const int *lda,const float *B, const int *ldb, 
                    const float *beta, float *C, const int *ldc);
void CUBLAS_SSYRK (const char *uplo, const char *trans, const int *n,
                   const int *k, const float *alpha, const float *A,
                   const int *lda, const float *beta, float *C,
                   const int *ldc);
void CUBLAS_STRMM (const char *side, const char *uplo, const char *transa, 
                   const char *diag, const int *m, const int *n,
                   const float *alpha, const float *A, const int *lda,
                   float *B, const int *ldb);
void CUBLAS_STRSM (const char *side, const char *uplo, const char *transa, 
                   const char *diag, const int *m, const int *n,
                   const float *alpha, const float *A, const int *lda,
                   float *B, const int *ldb);

void CUBLAS_CGEMM (const char *transa, const char *transb, const int *m,
                   const int *n, const int *k, const cuComplex *alpha,
                   const cuComplex *A, const int *lda, const cuComplex *B,
                   const int *ldb, const cuComplex *beta, cuComplex *C, 
                   const int *ldc);
#if defined(__cplusplus)
}
#endif /* __cplusplus */

/*---------------------------------------------------------------------------*/
/*---------------------------------- BLAS1 ----------------------------------*/
/*---------------------------------------------------------------------------*/

int CUBLAS_ISAMAX (const int *n, const float *x, const int *incx)
{
    float *devPtrx;
    int retVal;

    cublasAlloc (*n * abs(*incx), sizeof(x[0]), (void**)&devPtrx);
    cublasSetVector (*n, sizeof(x[0]), x, *incx, devPtrx, *incx);
    retVal = cublasIsamax (*n, devPtrx, *incx);
    cublasFree (devPtrx);
    return retVal;
}

int CUBLAS_ISAMIN (const int *n, const float *x, const int *incx)
{
    float *devPtrx;
    int retVal;

    cublasAlloc (*n * abs(*incx), sizeof(x[0]), (void**)&devPtrx);
    cublasSetVector (*n, sizeof(x[0]), x, *incx, devPtrx, *incx);
    retVal = cublasIsamin (*n, devPtrx, *incx);
    cublasFree (devPtrx);
    return retVal;
}

#if CUBLAS_FORTRAN_COMPILER==CUBLAS_G77
double CUBLAS_SASUM (const int *n, const float *x, const int *incx)
#else
float CUBLAS_SASUM (const int *n, const float *x, const int *incx)
#endif
{
    float *devPtrx;
    float retVal;
    
    cublasAlloc (*n * abs(*incx), sizeof(x[0]), (void**)&devPtrx);
    cublasSetVector (*n, sizeof(x[0]), x, *incx, devPtrx, *incx);
    retVal = cublasSasum (*n, devPtrx, *incx);
    cublasFree (devPtrx);
    return retVal;
}

void CUBLAS_SAXPY (const int *n, const float *alpha, const float *x, 
                   const int *incx, float *y, const int *incy)
{
    float *devPtrx, *devPtry;

    cublasAlloc (*n * abs(*incx), sizeof(x[0]), (void**)&devPtrx);
    cublasAlloc (*n * abs(*incy), sizeof(y[0]), (void**)&devPtry);
    cublasSetVector (*n, sizeof(x[0]), x, abs(*incx), devPtrx, abs(*incx));
    cublasSetVector (*n, sizeof(y[0]), y, abs(*incy), devPtry, abs(*incy));
    cublasSaxpy (*n, *alpha, devPtrx, *incx, devPtry, *incy);
    cublasGetVector (*n, sizeof(y[0]), devPtry, abs(*incy), y, abs(*incy));
    cublasFree (devPtrx);
    cublasFree (devPtry);
}

void CUBLAS_SCOPY (const int *n, const float *x, const int *incx, float *y,
                   const int *incy)
{
    float *devPtrx, *devPtry;

    cublasAlloc (*n * abs(*incx), sizeof(x[0]), (void**)&devPtrx);
    cublasAlloc (*n * abs(*incy), sizeof(y[0]), (void**)&devPtry);
    cublasSetVector (*n, sizeof(x[0]), x, abs(*incx), devPtrx, abs(*incx));
    cublasSetVector (*n, sizeof(y[0]), y, abs(*incy), devPtry, abs(*incy));
    cublasScopy (*n, devPtrx, *incx, devPtry, *incy);
    cublasGetVector (*n, sizeof(y[0]), devPtry, abs(*incy), y, abs(*incy));
    cublasFree (devPtrx);
    cublasFree (devPtry);
}

#if CUBLAS_FORTRAN_COMPILER==CUBLAS_G77
double CUBLAS_SDOT (const int *n, const float *x, const int *incx, float *y,
                    const int *incy)
#else
float CUBLAS_SDOT (const int *n, const float *x, const int *incx, float *y,
                   const int *incy)
#endif
{
    float *devPtrx, *devPtry, retVal;

    cublasAlloc (*n * abs(*incx), sizeof(x[0]), (void**)&devPtrx);
    cublasAlloc (*n * abs(*incy), sizeof(y[0]), (void**)&devPtry);
    cublasSetVector (*n, sizeof(x[0]), x, abs(*incx), devPtrx, abs(*incx));
    cublasSetVector (*n, sizeof(y[0]), y, abs(*incy), devPtry, abs(*incy));
    retVal = cublasSdot (*n, devPtrx, *incx, devPtry, *incy);
    cublasFree (devPtrx);
    cublasFree (devPtry);
    return retVal;
}

#if CUBLAS_FORTRAN_COMPILER==CUBLAS_G77
double CUBLAS_SNRM2 (const int *n, const float *x, const int *incx)
#else
float CUBLAS_SNRM2 (const int *n, const float *x, const int *incx)
#endif
{
    float *devPtrx;
    float retVal;
    
    cublasAlloc (*n * abs(*incx), sizeof(x[0]), (void**)&devPtrx);
    cublasSetVector (*n, sizeof(x[0]), x, *incx, devPtrx, *incx);
    retVal = cublasSnrm2 (*n, devPtrx, *incx);
    cublasFree (devPtrx);
    return retVal;
}

void CUBLAS_SROT (const int *n, float *x, const int *incx, float *y, 
                  const int *incy, const float *sc, const float *ss)
{
    float *devPtrx, *devPtry;

    cublasAlloc (*n * abs(*incx), sizeof(x[0]), (void**)&devPtrx);
    cublasAlloc (*n * abs(*incy), sizeof(y[0]), (void**)&devPtry);
    cublasSetVector (*n, sizeof(x[0]), x, abs(*incx), devPtrx, abs(*incx));
    cublasSetVector (*n, sizeof(y[0]), y, abs(*incy), devPtry, abs(*incy));
    cublasSrot (*n, devPtrx, *incx, devPtry, *incy, *sc, *ss);
    cublasGetVector (*n, sizeof(x[0]), devPtrx, abs(*incx), x, abs(*incx));
    cublasGetVector (*n, sizeof(y[0]), devPtry, abs(*incy), y, abs(*incy));
    cublasFree (devPtrx);
    cublasFree (devPtry);
}

void CUBLAS_SROTG (float *sa, float *sb, float *sc, float *ss)
{
    cublasSrotg (sa, sb, sc, ss);
}

void CUBLAS_SROTM (const int *n, float *x, const int *incx, float *y, 
                   const int *incy, const float* sparam)
{
    float *devPtrx, *devPtry;

    cublasAlloc (*n * abs(*incx), sizeof(x[0]), (void**)&devPtrx);
    cublasAlloc (*n * abs(*incy), sizeof(y[0]), (void**)&devPtry);
    cublasSetVector (*n, sizeof(x[0]), x, abs(*incx), devPtrx, abs(*incx));
    cublasSetVector (*n, sizeof(y[0]), y, abs(*incy), devPtry, abs(*incy));
    cublasSrotm (*n, devPtrx, *incx, devPtry, *incy, sparam);
    cublasGetVector (*n, sizeof(x[0]), devPtrx, abs(*incx), x, abs(*incx));
    cublasGetVector (*n, sizeof(y[0]), devPtry, abs(*incy), y, abs(*incy));
    cublasFree (devPtrx);
    cublasFree (devPtry);
}

void CUBLAS_SROTMG (float *sd1, float *sd2, float *sx1, const float *sy1,
                    float* sparam)
{
    cublasSrotmg (sd1, sd2, sx1, sy1, sparam);
}

void CUBLAS_SSCAL (const int *n, const float *alpha, float *x, const int *incx)
{
    float *devPtrx;
    
    cublasAlloc (*n * abs(*incx), sizeof(x[0]), (void**)&devPtrx);
    cublasSetVector (*n, sizeof(x[0]), x, *incx, devPtrx, *incx);
    cublasSscal (*n, *alpha, devPtrx, *incx);
    cublasGetVector (*n, sizeof(x[0]), devPtrx, *incx, x, *incx);
    cublasFree (devPtrx); 
}

void CUBLAS_SSWAP (const int *n, float *x, const int *incx, float *y, 
                   const int *incy)
{
    float *devPtrx, *devPtry;

    cublasAlloc (*n * abs(*incx), sizeof(x[0]), (void**)&devPtrx);
    cublasAlloc (*n * abs(*incy), sizeof(y[0]), (void**)&devPtry);
    cublasSetVector (*n, sizeof(x[0]), x, abs(*incx), devPtrx, abs(*incx));
    cublasSetVector (*n, sizeof(y[0]), y, abs(*incy), devPtry, abs(*incy));
    cublasSswap (*n, devPtrx, *incx, devPtry, *incy);
    cublasGetVector (*n, sizeof(x[0]), devPtrx, abs(*incx), x, abs(*incx));
    cublasGetVector (*n, sizeof(y[0]), devPtry, abs(*incy), y, abs(*incy));
    cublasFree (devPtrx);
    cublasFree (devPtry);
}

void CUBLAS_CAXPY (const int *n, const cuComplex *alpha, const cuComplex *x, 
                   const int *incx, cuComplex *y, const int *incy)
{
    cuComplex *devPtrx, *devPtry;

    cublasAlloc (*n * abs(*incx), sizeof(x[0]), (void**)&devPtrx);
    cublasAlloc (*n * abs(*incy), sizeof(y[0]), (void**)&devPtry);
    cublasSetVector (*n, sizeof(x[0]), x, abs(*incx), devPtrx, abs(*incx));
    cublasSetVector (*n, sizeof(y[0]), y, abs(*incy), devPtry, abs(*incy));
    cublasCaxpy (*n, *alpha, devPtrx, *incx, devPtry, *incy);
    cublasGetVector (*n, sizeof(y[0]), devPtry, abs(*incy), y, abs(*incy));
    cublasFree (devPtrx);
    cublasFree (devPtry);
}

void CUBLAS_CCOPY (const int *n, const cuComplex *x, const int *incx, 
                   cuComplex *y, const int *incy)
{
    cuComplex *devPtrx, *devPtry;

    cublasAlloc (*n * abs(*incx), sizeof(x[0]), (void**)&devPtrx);
    cublasAlloc (*n * abs(*incy), sizeof(y[0]), (void**)&devPtry);
    cublasSetVector (*n, sizeof(x[0]), x, abs(*incx), devPtrx, abs(*incx));
    cublasSetVector (*n, sizeof(y[0]), y, abs(*incy), devPtry, abs(*incy));
    cublasCcopy (*n, devPtrx, *incx, devPtry, *incy);
    cublasGetVector (*n, sizeof(y[0]), devPtry, abs(*incy), y, abs(*incy));
    cublasFree (devPtrx);
    cublasFree (devPtry);
}

void CUBLAS_CROT (const int *n, cuComplex *x, const int *incx, cuComplex *y, 
                   const int *incy, const float *sc, const cuComplex *cs)
{
    cuComplex *devPtrx, *devPtry;

    cublasAlloc (*n * abs(*incx), sizeof(x[0]), (void**)&devPtrx);
    cublasAlloc (*n * abs(*incy), sizeof(y[0]), (void**)&devPtry);
    cublasSetVector (*n, sizeof(x[0]), x, abs(*incx), devPtrx, abs(*incx));
    cublasSetVector (*n, sizeof(y[0]), y, abs(*incy), devPtry, abs(*incy));
    cublasCrot (*n, devPtrx, *incx, devPtry, *incy, *sc, *cs);
    cublasGetVector (*n, sizeof(x[0]), devPtrx, abs(*incx), x, abs(*incx));
    cublasGetVector (*n, sizeof(y[0]), devPtry, abs(*incy), y, abs(*incy));
    cublasFree (devPtrx);
    cublasFree (devPtry);
}

void CUBLAS_CROTG (cuComplex *ca, const cuComplex *cb, float *sc,
                   cuComplex *cs)
{
    cublasCrotg (ca, *cb, sc, cs);
}

void CUBLAS_CSCAL (const int *n, const cuComplex *alpha, cuComplex *x, 
                   const int *incx)
{
    cuComplex *devPtrx;
    
    cublasAlloc (*n * abs(*incx), sizeof(x[0]), (void**)&devPtrx);
    cublasSetVector (*n, sizeof(x[0]), x, *incx, devPtrx, *incx);
    cublasCscal (*n, *alpha, devPtrx, *incx);
    cublasGetVector (*n, sizeof(x[0]), devPtrx, *incx, x, *incx);
    cublasFree (devPtrx); 
}

void CUBLAS_CSROT (const int *n, cuComplex *x, const int *incx, cuComplex *y, 
                   const int *incy, const float *sc, const float *ss)
{
    cuComplex *devPtrx, *devPtry;

    cublasAlloc (*n * abs(*incx), sizeof(x[0]), (void**)&devPtrx);
    cublasAlloc (*n * abs(*incy), sizeof(y[0]), (void**)&devPtry);
    cublasSetVector (*n, sizeof(x[0]), x, abs(*incx), devPtrx, abs(*incx));
    cublasSetVector (*n, sizeof(y[0]), y, abs(*incy), devPtry, abs(*incy));
    cublasCsrot (*n, devPtrx, *incx, devPtry, *incy, *sc, *ss);
    cublasGetVector (*n, sizeof(x[0]), devPtrx, abs(*incx), x, abs(*incx));
    cublasGetVector (*n, sizeof(y[0]), devPtry, abs(*incy), y, abs(*incy));
    cublasFree (devPtrx);
    cublasFree (devPtry);
}

void CUBLAS_CSSCAL (const int *n, const float *alpha, cuComplex *x, 
                    const int *incx)
{
    cuComplex *devPtrx;
    
    cublasAlloc (*n * abs(*incx), sizeof(x[0]), (void**)&devPtrx);
    cublasSetVector (*n, sizeof(x[0]), x, *incx, devPtrx, *incx);
    cublasCsscal (*n, *alpha, devPtrx, *incx);
    cublasGetVector (*n, sizeof(x[0]), devPtrx, *incx, x, *incx);
    cublasFree (devPtrx); 
}

void CUBLAS_CSWAP (const int *n, cuComplex *x, const int *incx, cuComplex *y,
                   const int *incy)
{
    cuComplex *devPtrx, *devPtry;

    cublasAlloc (*n * abs(*incx), sizeof(x[0]), (void**)&devPtrx);
    cublasAlloc (*n * abs(*incy), sizeof(y[0]), (void**)&devPtry);
    cublasSetVector (*n, sizeof(x[0]), x, abs(*incx), devPtrx, abs(*incx));
    cublasSetVector (*n, sizeof(y[0]), y, abs(*incy), devPtry, abs(*incy));
    cublasCswap (*n, devPtrx, *incx, devPtry, *incy);
    cublasGetVector (*n, sizeof(x[0]), devPtrx, abs(*incx), x, abs(*incx));
    cublasGetVector (*n, sizeof(y[0]), devPtry, abs(*incy), y, abs(*incy));
    cublasFree (devPtrx);
    cublasFree (devPtry);
}

void CUBLAS_CDOTU (cuComplex *retVal, const int *n, const cuComplex *x, 
                   const int *incx, const cuComplex *y, const int *incy)
{
    cuComplex *devPtrx, *devPtry;

    cublasAlloc (*n * abs(*incx), sizeof(x[0]), (void**)&devPtrx);
    cublasAlloc (*n * abs(*incy), sizeof(y[0]), (void**)&devPtry);
    cublasSetVector (*n, sizeof(x[0]), x, abs(*incx), devPtrx, abs(*incx));
    cublasSetVector (*n, sizeof(y[0]), y, abs(*incy), devPtry, abs(*incy));
    *retVal = cublasCdotu (*n, devPtrx, *incx, devPtry, *incy);
    cublasFree (devPtrx);
    cublasFree (devPtry);
}

void CUBLAS_CDOTC (cuComplex *retVal, const int *n, const cuComplex *x, 
                   const int *incx, const cuComplex *y, const int *incy)
{
    cuComplex *devPtrx, *devPtry;

    cublasAlloc (*n * abs(*incx), sizeof(x[0]), (void**)&devPtrx);
    cublasAlloc (*n * abs(*incy), sizeof(y[0]), (void**)&devPtry);
    cublasSetVector (*n, sizeof(x[0]), x, abs(*incx), devPtrx, abs(*incx));
    cublasSetVector (*n, sizeof(y[0]), y, abs(*incy), devPtry, abs(*incy));
    *retVal = cublasCdotc (*n, devPtrx, *incx, devPtry, *incy);
    cublasFree (devPtrx);
    cublasFree (devPtry);
}

int CUBLAS_ICAMAX (const int *n, const cuComplex *x, const int *incx)
{
    cuComplex *devPtrx;
    int retVal;

    cublasAlloc (*n * abs(*incx), sizeof(x[0]), (void**)&devPtrx);
    cublasSetVector (*n, sizeof(x[0]), x, *incx, devPtrx, *incx);
    retVal = cublasIcamax (*n, devPtrx, *incx);
    cublasFree (devPtrx);
    return retVal;
}

int CUBLAS_ICAMIN (const int *n, const cuComplex *x, const int *incx)
{
    cuComplex *devPtrx;
    int retVal;

    cublasAlloc (*n * abs(*incx), sizeof(x[0]), (void**)&devPtrx);
    cublasSetVector (*n, sizeof(x[0]), x, *incx, devPtrx, *incx);
    retVal = cublasIcamin (*n, devPtrx, *incx);
    cublasFree (devPtrx);
    return retVal;
}

#if CUBLAS_FORTRAN_COMPILER==CUBLAS_G77
double CUBLAS_SCASUM (const int *n, const cuComplex *x, const int *incx)
#else
float CUBLAS_SCASUM (const int *n, const cuComplex *x, const int *incx)
#endif
{
    cuComplex *devPtrx;
    float retVal;
    
    cublasAlloc (*n * abs(*incx), sizeof(x[0]), (void**)&devPtrx);
    cublasSetVector (*n, sizeof(x[0]), x, *incx, devPtrx, *incx);
    retVal = cublasScasum (*n, devPtrx, *incx);
    cublasFree (devPtrx);
    return retVal;
}

#if CUBLAS_FORTRAN_COMPILER==CUBLAS_G77
double CUBLAS_SCNRM2 (const int *n, const cuComplex *x, const int *incx)
#else
float CUBLAS_SCNRM2 (const int *n, const cuComplex *x, const int *incx)
#endif
{
    cuComplex *devPtrx;
    float retVal;
    
    cublasAlloc (*n * abs(*incx), sizeof(x[0]), (void**)&devPtrx);
    cublasSetVector (*n, sizeof(x[0]), x, *incx, devPtrx, *incx);
    retVal = cublasScnrm2 (*n, devPtrx, *incx);
    cublasFree (devPtrx);
    return retVal;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------- BLAS2 ----------------------------------*/
/*---------------------------------------------------------------------------*/

void CUBLAS_SGBMV (const char *trans, const int *m, const int *n, 
                   const int *kl, const int *ku, const float *alpha, 
                   const float *A, const int *lda, const float *x, 
                   const int *incx, const float *beta, float *y, 
                   const int *incy)
{
    float *devPtrx, *devPtry, *devPtrA;

    /*  X      - REAL             array of DIMENSION at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
     *           and at least
     *           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
     */
    if (toupper(trans[0]) == 'N') {
        cublasAlloc (*n * abs(*incx), sizeof(x[0]), (void**)&devPtrx);
        cublasSetVector (*n, sizeof(x[0]), x, abs(*incx), devPtrx, abs(*incx));
    } else {
        cublasAlloc (*m * abs(*incx), sizeof(x[0]), (void**)&devPtrx);
        cublasSetVector (*m, sizeof(x[0]), x, abs(*incx), devPtrx, abs(*incx));
    }
    /*  Y      - REAL             array of DIMENSION at least
     *           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
     *           and at least
     *           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
     */
    if (toupper(trans[0]) == 'N') {
        cublasAlloc (*m * abs(*incy), sizeof(y[0]), (void**)&devPtry);
        cublasSetVector (*m, sizeof(y[0]), y, abs(*incy), devPtry, abs(*incy));
    } else {
        cublasAlloc (*n * abs(*incy), sizeof(y[0]), (void**)&devPtry);
        cublasSetVector (*n, sizeof(y[0]), y, abs(*incy), devPtry, abs(*incy));
    }       
    /* A      - REAL             array of DIMENSION ( LDA, n ). 
     * Before entry, the leading ( kl + ku + 1 ) by n part of the
     * array A must contain the matrix of coefficients, supplied
     */
    cublasAlloc ((*lda) * (*n), sizeof(devPtrA[0]), (void**)&devPtrA);
    cublasSetMatrix (imin(*kl+*ku+1,*lda), *n, sizeof(A[0]), A, *lda, devPtrA, 
                     *lda);
    cublasSgbmv (trans[0], *m, *n, *kl, *ku, *alpha, devPtrA, *lda, devPtrx, 
                 *incx, *beta, devPtry, *incy);
    if (toupper(trans[0]) == 'N') {
        cublasGetVector (*m, sizeof(y[0]), devPtry, abs(*incy), y, abs(*incy));
    } else {
        cublasGetVector (*n, sizeof(y[0]), devPtry, abs(*incy), y, abs(*incy));
    }       
    cublasFree (devPtrA);
    cublasFree (devPtry);
    cublasFree (devPtrx);
}

void CUBLAS_SGEMV (const char *trans, const int *m, const int *n,
                   const float *alpha, const float *A, const int *lda,
                   const float *x, const int *incx, const float *beta,
                   float *y, const int *incy)
{
    float *devPtrA, *devPtrx, *devPtry;
    
    /*  X      - REAL             array of DIMENSION at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
     *           and at least
     *           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
     */
    if (toupper(trans[0]) == 'N') {
        cublasAlloc (*n * abs(*incx), sizeof(x[0]), (void**)&devPtrx);
        cublasSetVector (*n, sizeof(x[0]), x, abs(*incx), devPtrx, abs(*incx));
    } else {
        cublasAlloc (*m * abs(*incx), sizeof(x[0]), (void**)&devPtrx);
        cublasSetVector (*m, sizeof(x[0]), x, abs(*incx), devPtrx, abs(*incx));
    }
    /*  Y      - REAL             array of DIMENSION at least
     *           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
     *           and at least
     *           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
     */
    if (toupper(trans[0]) == 'N') {
        cublasAlloc (*m * abs(*incy), sizeof(y[0]), (void**)&devPtry);
        cublasSetVector (*m, sizeof(y[0]), y, abs(*incy), devPtry, abs(*incy));
    } else {
        cublasAlloc (*n * abs(*incy), sizeof(y[0]), (void**)&devPtry);
        cublasSetVector (*n, sizeof(y[0]), y, abs(*incy), devPtry, abs(*incy));
    }       
    /*  A      - REAL             array of DIMENSION ( LDA, n ).
     *           Before entry, the leading m by n part of the array A must
     *           contain the matrix of coefficients.
     */
    cublasAlloc ((*lda) * (*n), sizeof(devPtrA[0]), (void**)&devPtrA);
    cublasSetMatrix (imin(*m,*lda), *n, sizeof(A[0]), A, *lda, devPtrA, *lda);
    cublasSgemv (trans[0], *m, *n, *alpha, devPtrA, *lda, devPtrx, *incx,
                 *beta, devPtry, *incy);
    if (toupper(trans[0]) == 'N') {
        cublasGetVector (*m, sizeof(y[0]), devPtry, abs(*incy), y, abs(*incy));
    } else {
        cublasGetVector (*n, sizeof(y[0]), devPtry, abs(*incy), y, abs(*incy));
    }       
    cublasFree (devPtrA);
    cublasFree (devPtry);
    cublasFree (devPtrx);
}

void CUBLAS_SGER (const int *m, const int *n, const float *alpha, 
                  const float *x, const int *incx, const float *y,
                  const int *incy, float *A, const int *lda)
{
    float *devPtrA, *devPtrx, *devPtry;

    cublasAlloc (*m * abs(*incx), sizeof(devPtrx[0]), (void**)&devPtrx);
    cublasAlloc (*n * abs(*incy), sizeof(devPtry[0]), (void**)&devPtry);
    cublasSetVector (*m* abs(*incx), sizeof(x[0]), x, 1, devPtrx, 1);
    cublasSetVector (*n* abs(*incy), sizeof(y[0]), y, 1, devPtry, 1);

    // REAL array of DIMENSION ( LDA, n ).
    //      Before entry, the leading m by n part of the array A must
    //      contain the matrix of coefficients. On exit, A is
    cublasAlloc ((*lda) * (*n), sizeof(devPtrA[0]), (void**)&devPtrA);
    cublasSetMatrix (imin(*m,*lda), *n, sizeof(A[0]), A, *lda, devPtrA, *lda);
    cublasSger (*m, *n, *alpha, devPtrx, *incx, devPtry, *incy, devPtrA, *lda);
    cublasGetMatrix (imin(*m,*lda), *n, sizeof(A[0]), devPtrA, *lda, A, *lda);
    cublasFree (devPtrA);
    cublasFree (devPtrx);
    cublasFree (devPtry);
}

void CUBLAS_SSBMV (const char *uplo, const int *n, const int *k, 
                   const float *alpha, const float *A, const int *lda,
                   const float *x, const int *incx, const float *beta, 
                   float *y, const int *incy)
{
    float *devPtrA, *devPtrx, *devPtry;
    
    /*  X      - REAL             array of DIMENSION at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     */
    cublasAlloc (*n * abs(*incx), sizeof(devPtrx[0]), (void**)&devPtrx);
    cublasSetVector (*n* abs(*incx), sizeof(x[0]), x, 1, devPtrx, 1);
    /*  Y      - REAL             array of DIMENSION at least
     *           ( 1 + ( n - 1 )*abs( INCY ) ).
     */
    cublasAlloc (*n * abs(*incy), sizeof(devPtry[0]), (void**)&devPtry);
    cublasSetVector (*n* abs(*incy), sizeof(y[0]), y, 1, devPtry, 1);    
    /*  A      - REAL             array of DIMENSION ( LDA, n ).
     *           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
     *           by n part of the array A must contain the upper triangular
     */
    cublasAlloc ((*lda) * (*n), sizeof(devPtrA[0]), (void**)&devPtrA);
    cublasSetMatrix (imin(*k+1,*lda), *n, sizeof(A[0]), A, *lda, devPtrA, *lda);
    cublasSsbmv (uplo[0], *n, *k, *alpha, devPtrA, *lda, devPtrx, *incx, *beta,
                 devPtry, *incy);
    cublasGetVector (*n* abs(*incy), sizeof(y[0]), devPtry, 1, y, 1);
    cublasFree (devPtrA);
    cublasFree (devPtrx);
    cublasFree (devPtry);
}

void CUBLAS_SSPMV (const char *uplo, const int *n, const float *alpha,
                   const float *AP, const float *x, const int *incx, 
                   const float *beta, float *y, const int *incy)
{
    float *devPtrAP, *devPtrx, *devPtry;

    /*  X      - REAL             array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     */
    cublasAlloc (*n * abs(*incx), sizeof(devPtrx[0]), (void**)&devPtrx);
    cublasSetVector (*n* abs(*incx), sizeof(x[0]), x, 1, devPtrx, 1);
    /*  Y      - REAL             array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCY ) ).
     */
    cublasAlloc (*n * abs(*incy), sizeof(devPtry[0]), (void**)&devPtry);
    cublasSetVector (*n* abs(*incy), sizeof(y[0]), y, 1, devPtry, 1);
    /*  AP     - REAL             array of DIMENSION at least
     *           ( ( n*( n + 1 ) )/2 ).
     *           Before entry with UPLO = 'U' or 'u', the array AP must
     *           contain the upper triangular part of the symmetric matrix
     */
    cublasAlloc (((*n)*(*n+1))/2, sizeof(devPtrAP[0]), (void**)&devPtrAP);
    cublasSetVector (((*n)*(*n+1))/2, sizeof(AP[0]), AP, 1, devPtrAP, 1);
    cublasSspmv (*uplo, *n, *alpha, devPtrAP, devPtrx, *incx, *beta, devPtry,
                 *incy);
    cublasGetVector ((*n) * abs(*incy), sizeof(y[0]), devPtry, 1, y, 1);
    cublasFree (devPtrAP);
    cublasFree (devPtrx);
    cublasFree (devPtry);
}

void CUBLAS_SSPR (const char *uplo, const int *n, const float *alpha, 
                    const float *x, const int *incx, float *AP)
{
    float *devPtrAP, *devPtrx;

    /*  X      - REAL             array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     */
    cublasAlloc (*n * abs(*incx), sizeof(devPtrx[0]), (void**)&devPtrx);
    cublasSetVector (*n* abs(*incx), sizeof(x[0]), x, 1, devPtrx, 1);
    /*  AP     - REAL             array of DIMENSION at least
     *           ( ( n*( n + 1 ) )/2 ).
     */
    cublasAlloc (((*n) * (*n+1))/2, sizeof(devPtrAP[0]), (void**)&devPtrAP);
    cublasSetVector (((*n) * (*n+1))/2, sizeof(AP[0]), AP, 1, devPtrAP, 1);
    cublasSspr (uplo[0], *n, *alpha, devPtrx, *incx, devPtrAP);
    cublasGetVector (((*n) * (*n+1))/2, sizeof(AP[0]), devPtrAP, 1, AP, 1);
    cublasFree (devPtrAP);
    cublasFree (devPtrx);
}

void CUBLAS_SSPR2 (const char *uplo, const int *n, const float *alpha,
                   const float *x, const int *incx, const float *y, 
                   const int *incy, float *AP)
{
    float *devPtrAP, *devPtrx, *devPtry;

    /*  X      - REAL             array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     */
    cublasAlloc (*n * abs(*incx), sizeof(devPtrx[0]), (void**)&devPtrx);
    cublasSetVector (*n* abs(*incx), sizeof(x[0]), x, 1, devPtrx, 1);
    /*  Y      - REAL             array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCY ) ).
     */
    cublasAlloc (*n * abs(*incy), sizeof(devPtry[0]), (void**)&devPtry);
    cublasSetVector (*n* abs(*incy), sizeof(y[0]), y, 1, devPtry, 1);
   /*  AP     - REAL             array of DIMENSION at least
    *           ( ( n*( n + 1 ) )/2 ).
    *           Before entry with  UPLO = 'U' or 'u', the array AP must
    *           contain the upper triangular part of the symmetric matrix
    */
    cublasAlloc (((*n) * (*n+1))/2, sizeof(devPtrAP[0]), (void**)&devPtrAP);
    cublasSetVector (((*n) * (*n+1))/2, sizeof(AP[0]), AP, 1, devPtrAP, 1);
    cublasSspr2 (uplo[0], *n, *alpha, devPtrx, *incx, devPtry, *incy,devPtrAP);
    cublasGetVector (((*n) * (*n+1))/2, sizeof(AP[0]), devPtrAP, 1, AP, 1);
    cublasFree (devPtrAP);
    cublasFree (devPtrx);
    cublasFree (devPtry);
}

void CUBLAS_SSYMV (const char *uplo, const int *n, const float *alpha,
                   const float *A, const int *lda, const float *x, 
                   const int *incx, const float *beta, float *y, 
                   const int *incy)
{
    float *devPtrA, *devPtrx, *devPtry;
    
    /*  X      - REAL             array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     */
    cublasAlloc (*n * abs(*incx), sizeof(devPtrx[0]), (void**)&devPtrx);
    cublasSetVector (*n* abs(*incx), sizeof(x[0]), x, 1, devPtrx, 1);
    /*  Y      - REAL             array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCY ) ).
     */
    cublasAlloc (*n * abs(*incy), sizeof(devPtry[0]), (void**)&devPtry);
    cublasSetVector (*n* abs(*incy), sizeof(y[0]), y, 1, devPtry, 1);
    /*  A      - REAL             array of DIMENSION ( LDA, n ).
     *           Before entry with  UPLO = 'U' or 'u', the leading n by n
     *           upper triangular part of the array A must contain the upper
     */
    cublasAlloc ((*lda) * (*n), sizeof(devPtrA[0]), (void**)&devPtrA);
    cublasSetMatrix (imin(*n,*lda), *n, sizeof(A[0]), A, *lda, devPtrA, *lda);
    cublasSsymv (uplo[0], *n, *alpha, devPtrA, *lda, devPtrx, *incx, *beta,
                 devPtry, *incy);
    cublasGetVector (*n * abs(*incy), sizeof(y[0]), devPtry, 1, y, 1);
    cublasFree (devPtrA);
    cublasFree (devPtrx);
    cublasFree (devPtry);
}

void CUBLAS_SSYR (const char *uplo, const int *n, const float *alpha, 
                  const float *x, const int *incx, float *A, const int *lda)
{
    float *devPtrA, *devPtrx;

    /*  X      - REAL             array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     */
    cublasAlloc (*n * abs(*incx), sizeof(devPtrx[0]), (void**)&devPtrx);
    cublasSetVector (*n* abs(*incx), sizeof(x[0]), x, 1, devPtrx, 1);
    /*  A      - REAL             array of DIMENSION ( LDA, n ).
     *           Before entry with  UPLO = 'U' or 'u', the leading n by n
     *           upper triangular part of the array A must contain the upper
     */
    cublasAlloc ((*lda) * (*n), sizeof(devPtrA[0]), (void**)&devPtrA);
    cublasSetMatrix (imin(*n,*lda), *n, sizeof(A[0]), A, *lda, devPtrA, *lda);
    cublasSsyr (uplo[0], *n, *alpha, devPtrx, *incx, devPtrA, *lda);
    cublasGetMatrix (imin(*n,*lda), *n, sizeof(A[0]), devPtrA, *lda, A, *lda);
    cublasFree (devPtrA);
    cublasFree (devPtrx);
}

void CUBLAS_SSYR2 (const char *uplo, const int *n, const float *alpha,
                   const float *x, const int *incx, const float *y,
                   const int *incy, float *A, const int *lda)
{
    float *devPtrA, *devPtrx, *devPtry;

    /*  X      - REAL             array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     */ 
    cublasAlloc (*n * abs(*incx), sizeof(devPtrx[0]), (void**)&devPtrx);
    cublasSetVector (*n * abs(*incx), sizeof(x[0]), x, 1, devPtrx, 1);
    /*
     *  Y      - REAL             array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCY ) ).
     */
    cublasAlloc (*n * abs(*incy), sizeof(devPtry[0]), (void**)&devPtry);
    cublasSetVector (*n* abs(*incy), sizeof(y[0]), y, 1, devPtry, 1);
    /*  A      - REAL             array of DIMENSION ( LDA, n ).
     *           Before entry with  UPLO = 'U' or 'u', the leading n by n
     *           upper triangular part of the array A must contain the upper
     */
    cublasAlloc ((*lda) * (*n), sizeof(devPtrA[0]), (void**)&devPtrA);
    cublasSetMatrix (imin(*n,*lda), *n, sizeof(A[0]), A, *lda, devPtrA, *lda);
    cublasSsyr2 (uplo[0], *n, *alpha, devPtrx, *incx, devPtry, *incy, devPtrA,
                 *lda);
    cublasGetMatrix (imin(*n,*lda), *n, sizeof(A[0]), devPtrA, *lda, A, *lda);
    cublasFree (devPtrA);
    cublasFree (devPtrx);
    cublasFree (devPtry);
}

void CUBLAS_STBMV (const char *uplo, const char *trans, const char *diag,
                   const int *n, const int *k, const float *A, const int *lda,
                   float *x, const int *incx)
{
    float *devPtrA, *devPtrx;
    
    /*  X      - REAL             array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     */
    cublasAlloc (*n * abs(*incx), sizeof(devPtrx[0]), (void**)&devPtrx);
    cublasSetVector (*n* abs(*incx), sizeof(x[0]), x, 1, devPtrx, 1);
    /*  A      - REAL             array of DIMENSION ( LDA, n ).
     *           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
     *           by n part of the array A must contain the upper triangular
     */
    cublasAlloc ((*lda) * (*n), sizeof(devPtrA[0]), (void**)&devPtrA);
    cublasSetMatrix (imin(*k+1,*lda), *n, sizeof(A[0]), A, *lda, devPtrA, *lda);
    cublasStbmv (uplo[0], trans[0], diag[0], *n, *k, devPtrA, *lda, devPtrx, 
                 *incx);
    cublasGetVector (*n* abs(*incx), sizeof(x[0]), devPtrx, 1, x, 1);
    cublasFree (devPtrA);
    cublasFree (devPtrx);
}

void CUBLAS_STBSV (const char *uplo, const char *trans, const char *diag,
                   const int *n, const int *k, const float *A, const int *lda,
                   float *x, const int *incx)
{
    float *devPtrA, *devPtrx;

    /*  X      - REAL             array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     */
    cublasAlloc (*n * abs(*incx), sizeof(devPtrx[0]), (void**)&devPtrx);
    cublasSetVector (*n * abs(*incx), sizeof(x[0]), x, 1, devPtrx, 1);
    /*  A      - REAL             array of DIMENSION ( LDA, n ).
     *           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
     *           by n part of the array A must contain the upper triangular
     */
    cublasAlloc ((*lda) * (*n), sizeof(devPtrA[0]), (void**)&devPtrA);
    cublasSetMatrix (imin(*k+1,*lda), *n, sizeof(A[0]), A, *lda, devPtrA, *lda);
    cublasStbsv (uplo[0], trans[0], diag[0], *n, *k, devPtrA, *lda, devPtrx,
                 *incx);
    cublasGetVector (*n * abs(*incx), sizeof(x[0]), devPtrx, 1, x, 1);
    cublasFree (devPtrA);
    cublasFree (devPtrx);
}

void CUBLAS_STPMV (const char *uplo, const char *trans, const char *diag,
                   const int *n, const float *AP, float *x, const int *incx)
{
    float *devPtrAP, *devPtrx;

    /*  X      - REAL             array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     */
    cublasAlloc (*n * abs(*incx), sizeof(devPtrx[0]), (void**)&devPtrx);
    cublasSetVector (*n* abs(*incx), sizeof(x[0]), x, 1, devPtrx, 1);
    /*  AP     - REAL             array of DIMENSION at least
     *           ( ( n*( n + 1 ) )/2 ).
     *           Before entry with  UPLO = 'U' or 'u', the array AP must
     */
    cublasAlloc (((*n) * (*n+1))/2, sizeof(devPtrAP[0]), (void**)&devPtrAP);
    cublasSetVector (((*n) * (*n+1))/2, sizeof(AP[0]), AP, 1, devPtrAP, 1);
    cublasStpmv (uplo[0], trans[0], diag[0], *n, devPtrAP, devPtrx, *incx);
    cublasGetVector (*n* abs(*incx), sizeof(x[0]), devPtrx, 1, x, 1);
    cublasFree (devPtrAP);
    cublasFree (devPtrx);
}

void CUBLAS_STPSV (const char *uplo, const char *trans, const char *diag,
                   const int *n, const float *AP, float *x, const int *incx)
{
    float *devPtrAP, *devPtrx;

     /*  X      - REAL             array of dimension at least
      *           ( 1 + ( n - 1 )*abs( INCX ) ).
      */
    cublasAlloc (*n * abs(*incx), sizeof(devPtrx[0]), (void**)&devPtrx);
    cublasSetVector (*n* abs(*incx), sizeof(x[0]), x, 1, devPtrx, 1);
    /*  AP     - REAL             array of DIMENSION at least
     *           ( ( n*( n + 1 ) )/2 ).
     *           Before entry with  UPLO = 'U' or 'u', the array AP must
     */
    cublasAlloc (((*n) * (*n+1))/2, sizeof(devPtrAP[0]), (void**)&devPtrAP);
    cublasSetVector (((*n) * (*n+1))/2, sizeof(AP[0]), AP, 1, devPtrAP, 1);
    cublasStpsv (uplo[0], trans[0], diag[0], *n, devPtrAP, devPtrx, *incx);
    cublasGetVector (*n* abs(*incx), sizeof(x[0]), devPtrx, 1, x, 1);
    cublasFree (devPtrAP);
    cublasFree (devPtrx);
}

void CUBLAS_STRMV (const char *uplo, const char *trans,
                            const char *diag, const int *n, const float *A,
                            const int *lda, float *x, const int *incx)
{
    float *devPtrA, *devPtrx;
    
    /*  X      - REAL             array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     */
    cublasAlloc (*n * abs(*incx), sizeof(devPtrx[0]), (void**)&devPtrx);
    cublasSetVector (*n* abs(*incx), sizeof(x[0]), x, 1, devPtrx, 1);
    /*  A      - REAL             array of DIMENSION ( LDA, n ).
     *           Before entry with  UPLO = 'U' or 'u', the leading n by n
     *           upper triangular part of the array A must contain the upper
     */
    cublasAlloc ((*lda) * (*n), sizeof(devPtrA[0]), (void**)&devPtrA);
    cublasSetMatrix (imin(*n,*lda), *n, sizeof(A[0]), A, *lda, devPtrA, *lda);
    cublasStrmv (uplo[0], trans[0], diag[0], *n, devPtrA, *lda, devPtrx,*incx);
    cublasGetVector (*n* abs(*incx), sizeof(x[0]), devPtrx, 1, x, 1);
    cublasFree (devPtrA);
    cublasFree (devPtrx);
}

void CUBLAS_STRSV (const char *uplo, const char *trans, const char *diag, 
                   const int *n, const float *A, const int *lda, float *x, 
                   const int *incx)
{
    float *devPtrA, *devPtrx;

    /*  X      - REAL             array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     */
    cublasAlloc (*n * abs(*incx), sizeof(devPtrx[0]), (void**)&devPtrx);
    cublasSetVector (*n* abs(*incx), sizeof(x[0]), x, 1, devPtrx, 1);
    /*  A      - REAL             array of DIMENSION ( LDA, n ).
     *           Before entry with  UPLO = 'U' or 'u', the leading n by n
     *           upper triangular part of the array A must contain the upper
     */
    cublasAlloc ((*lda) * (*n), sizeof(devPtrA[0]), (void**)&devPtrA);
    cublasSetMatrix (imin(*n,*lda), *n, sizeof(A[0]), A, *lda, devPtrA, *lda);
    cublasStrsv (uplo[0], trans[0], diag[0], *n, devPtrA, *lda, devPtrx,*incx);
    cublasGetVector (*n* abs(*incx), sizeof(x[0]), devPtrx, 1, x, 1);
    cublasFree (devPtrA);
    cublasFree (devPtrx);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------- BLAS3 ----------------------------------*/
/*---------------------------------------------------------------------------*/

void CUBLAS_SGEMM (const char *transa, const char *transb, const int *m,
                   const int *n, const int *k, const float *alpha,
                   const float *A, const int *lda, const float *B,
                   const int *ldb, const float *beta, float *C, const int *ldc)
{
    int ka, kb;
    float *devPtrA, *devPtrB, *devPtrC;

    /*  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is
     *           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
     *           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
     *           part of the array  A  must contain the matrix  A,  otherwise
     *           the leading  k by m  part of the array  A  must contain  the
     *           matrix A.
     */
    ka = (toupper(transa[0]) == 'N') ? *k : *m;
    cublasAlloc (*lda * ka, sizeof(devPtrA[0]), (void**)&devPtrA);
    if (toupper(transa[0]) == 'N') {
        cublasSetMatrix (imin(*m,*lda), *k, sizeof(A[0]), A, *lda, devPtrA, 
                         *lda);
    } else {
        cublasSetMatrix (imin(*k,*lda), *m, sizeof(A[0]), A, *lda, devPtrA, 
                         *lda);
    }

    /*  B      - REAL             array of DIMENSION ( LDB, kb ), where kb is
     *           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
     *           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
     *           part of the array  B  must contain the matrix  B,  otherwise
     *           the leading  n by k  part of the array  B  must contain  the
     *           matrix B.
     */
    kb = (toupper(transb[0]) == 'N') ? *n : *k;
    cublasAlloc (*ldb * kb, sizeof(devPtrB[0]), (void**)&devPtrB);
    if (toupper(transb[0]) == 'N') {
        cublasSetMatrix (imin(*k,*ldb), *n, sizeof(B[0]), B, *ldb, devPtrB, 
                         *ldb);
    } else {
        cublasSetMatrix (imin(*n,*ldb), *k, sizeof(B[0]), B, *ldb, devPtrB,
                         *ldb);
    }
    
    /*  C      - REAL             array of DIMENSION ( LDC, n ).
     *           Before entry, the leading  m by n  part of the array  C must
     *           contain the matrix  C,  except when  beta  is zero, in which
     *           case C need not be set on entry.
     *           On exit, the array  C  is overwritten by the  m by n  matrix
     */
    cublasAlloc ((*ldc) * (*n), sizeof(devPtrC[0]), (void**)&devPtrC);
    cublasSetMatrix (imin(*m,*ldc), *n, sizeof(C[0]), C, *ldc, devPtrC, *ldc);

    cublasSgemm (transa[0], transb[0], *m, *n, *k, *alpha, devPtrA, *lda, 
                 devPtrB, *ldb, *beta, devPtrC, *ldc);

    cublasGetMatrix (imin(*m,*ldc), *n, sizeof(C[0]), devPtrC, *ldc, C, *ldc);
    cublasFree (devPtrA);
    cublasFree (devPtrB);
    cublasFree (devPtrC);
}

void CUBLAS_SSYMM (const char *side, const char *uplo, const int *m, 
                   const int *n, const float *alpha, const float *A, 
                   const int *lda, const float *B, const int *ldb, 
                   const float *beta, float *C, const int *ldc)
{
    int ka;
    float *devPtrA, *devPtrB, *devPtrC;
    
    /*  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is
     *           m  when  SIDE = 'L' or 'l'  and is  n otherwise.
     *           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of
     *           the array  A  must contain the  symmetric matrix, [..]
     *           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of
     *           the array  A  must contain the  symmetric matrix, [..]
     */
    ka = (toupper(side[0]) == 'L') ? *m : *n;
    cublasAlloc ((*lda) * ka, sizeof(devPtrA[0]), (void**)&devPtrA);
    if (toupper(side[0]) == 'L') {
        cublasSetMatrix (imin(*m,*lda), *m, sizeof(A[0]), A, *lda, devPtrA, 
                         *lda);
    } else {
        cublasSetMatrix (imin(*n,*lda), *n, sizeof(A[0]), A, *lda, devPtrA,
                         *lda);
    }

    /*  B      - REAL             array of DIMENSION ( LDB, n ).
     *           Before entry, the leading  m by n part of the array  B  must
     *           contain the matrix B.
     */
    cublasAlloc ((*ldb) * (*n), sizeof(devPtrB[0]), (void**)&devPtrB);
    cublasSetMatrix (imin(*m,*ldb), *n, sizeof(B[0]), B, *ldb, devPtrB, *ldb);

    /*  C      - REAL             array of DIMENSION ( LDC, n ).
     *           Before entry, the leading  m by n  part of the array  C must
     *           contain the matrix  C,  except when  beta  is zero, in which
     *           case C need not be set on entry.
     */
    cublasAlloc ((*ldc) * (*n), sizeof(devPtrC[0]), (void**)&devPtrC);
    cublasSetMatrix (imin(*m,*ldc), *n, sizeof(C[0]), C, *ldc, devPtrC, *ldc);
    
    cublasSsymm (side[0], uplo[0], *m, *n, *alpha, devPtrA, *lda, devPtrB,
                 *ldb, *beta, devPtrC, *ldc);

    cublasGetMatrix (imin(*m,*ldc), *n, sizeof(C[0]), devPtrC, *ldc, C, *ldc);
    cublasFree (devPtrA);
    cublasFree (devPtrB);
    cublasFree (devPtrC);
}

void CUBLAS_SSYR2K (const char *uplo, const char *trans, const int *n,
                    const int *k, const float *alpha, const float *A, 
                    const int *lda, const float *B, const int *ldb, 
                    const float *beta, float *C, const int *ldc)
{
    int ka, kb;
    float *devPtrA, *devPtrB, *devPtrC;

    /*  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is
     *           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
     *           Before entry with  TRANS = 'N' or 'n',  the leading  n by k
     *           part of the array  A  must contain the matrix  A,  otherwise
     *           the leading  k by n  part of the array  A  must contain  the
     *           matrix A.
     */
    ka = (toupper(trans[0]) == 'N') ? *k : *n;
    cublasAlloc (*lda * ka, sizeof(devPtrA[0]), (void**)&devPtrA);
    if (toupper(trans[0]) == 'N') {
        cublasSetMatrix (imin(*n,*lda), *k, sizeof(A[0]), A, *lda, devPtrA,
                         *lda);
    } else {
        cublasSetMatrix (imin(*k,*lda), *n, sizeof(A[0]), A, *lda, devPtrA,
                         *lda);
    }

    /*  B      - REAL             array of DIMENSION ( LDB, kb ), where kb is
     *           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
     *           Before entry with  TRANS = 'N' or 'n',  the leading  n by k
     *           part of the array  B  must contain the matrix  B,  otherwise
     *           the leading  k by n  part of the array  B  must contain  the
     *           matrix B.
     */
    kb = (toupper(trans[0]) == 'N') ? *k : *n;
    cublasAlloc ((*ldb) * kb, sizeof(devPtrB[0]), (void**)&devPtrB);
    if (toupper(trans[0]) == 'N') {
        cublasSetMatrix (imin(*n,*ldb), *k, sizeof(B[0]), B, *ldb, devPtrB,
                         *ldb);
    } else {
        cublasSetMatrix (imin(*k,*ldb), *n, sizeof(B[0]), B, *ldb, devPtrB,
                         *ldb);
    }

    /* C      single precision array of dimensions (ldc, n). If uplo == 'U' or
     *        'u', the leading n x n triangular part of the array C must 
     *        contain the upper triangular part of the symmetric matrix C and 
     *        the strictly lower triangular part of C is not referenced. On 
     *        exit, the upper 
     */
    cublasAlloc ((*ldc) * (*n), sizeof(devPtrC[0]), (void**)&devPtrC);
    cublasSetMatrix (imin(*n,*ldc), *n, sizeof(C[0]), C, *ldc, devPtrC, *ldc);

    cublasSsyr2k (uplo[0], trans[0], *n, *k, *alpha, devPtrA, *lda, devPtrB, 
           *ldb, *beta, devPtrC, *ldc);

    cublasGetMatrix (imin(*n,*ldc), *n, sizeof(C[0]), devPtrC, *ldc, C, *ldc);
    cublasFree (devPtrA);
    cublasFree (devPtrB);
    cublasFree (devPtrC);
}

void CUBLAS_SSYRK (const char *uplo, const char *trans, const int *n, 
                   const int *k, const float *alpha, const float *A, 
                   const int *lda, const float *beta, float *C, const int *ldc)
{
    int ka;
    float *devPtrA, *devPtrC;

    /* A      single precision array of dimensions (lda, ka), where ka is k 
     *        when trans == 'N' or 'n', and is n otherwise. When trans == 'N' 
     *        or 'n', the leading n x k part of array A must contain the matrix
     *        A, otherwise the leading k x n part of the array must contain the
     *        matrix A.
     */
    ka = (toupper(trans[0]) == 'N') ? *k : *n;
    cublasAlloc (*lda * ka, sizeof(devPtrA[0]), (void**)&devPtrA);
    if (toupper(trans[0]) == 'N') {
        cublasSetMatrix (imin(*n,*lda), *k, sizeof(A[0]), A, *lda, devPtrA,
                         *lda);
    } else {
        cublasSetMatrix (imin(*k,*lda), *n, sizeof(A[0]), A, *lda, devPtrA,
                         *lda);
    }
    
    /* C      single precision array of dimensions (ldc, n). If uplo='U'or'u',
     *        the leading n x n triangular part of the array C must contain the
     *        upper triangular part of the symmetric matrix C and the strictly 
     */
    cublasAlloc ((*ldc) * (*n), sizeof(devPtrC[0]), (void**)&devPtrC);
    cublasSetMatrix (imin(*n,*ldc), *n, sizeof(C[0]), C, *ldc, devPtrC, *ldc);
    
    cublasSsyrk (uplo[0], trans[0], *n, *k, *alpha, devPtrA, *lda, *beta,
                 devPtrC, *ldc);

    cublasGetMatrix (imin(*n,*ldc), *n, sizeof(C[0]), devPtrC, *ldc, C, *ldc);
    cublasFree (devPtrA);
    cublasFree (devPtrC);
}

void CUBLAS_STRMM (const char *side, const char *uplo, const char *transa,
                   const char *diag, const int *m, const int *n,
                   const float *alpha, const float *A, const int *lda,
                   float *B, const int *ldb)
{
    int k;
    float *devPtrA, *devPtrB;

    /* A      single precision array of dimensions (lda, k). k = m if side =
     *        'L' or 'l', k = n if side = 'R' or 'r'. If uplo = 'U' or 'u'
     *        the leading k x k upper triangular part of the array A must
     *        contain the upper triangular matrix, and the strictly lower
     *        triangular part of A is not referenced. If uplo = 'L' or 'l'
     *        the leading k x k lower triangular part of the array A must
     *        contain the lower triangular matrix, and the strictly upper
     */
    k = (toupper(side[0]) == 'L') ? *m : *n;
    cublasAlloc (*lda * k, sizeof(devPtrA[0]), (void**)&devPtrA);
    if (toupper(side[0]) == 'L') {
        cublasSetMatrix (imin(k,*lda), k, sizeof(A[0]), A, *lda, devPtrA, *lda);
    } else {
        cublasSetMatrix (imin(k,*lda), k, sizeof(A[0]), A, *lda, devPtrA, *lda);
    }

    /* B      single precision array of dimensions (ldb, n). On entry, the 
     *        leading m x n part of the array contains the matrix B. It is
     *        overwritten with the transformed matrix on exit.
     */
    cublasAlloc ((*ldb) * (*n), sizeof(devPtrB[0]), (void**)&devPtrB);
    cublasSetMatrix (imin(*m,*ldb), *n, sizeof(B[0]), B, *ldb, devPtrB, *ldb);

    cublasStrmm (side[0], uplo[0], transa[0], diag[0], *m, *n, *alpha, devPtrA,
           *lda, devPtrB, *ldb);

    cublasGetMatrix (imin(*m,*ldb), *n, sizeof(B[0]), devPtrB, *ldb, B, *ldb);
    
    cublasFree (devPtrA);
    cublasFree (devPtrB);
}

void CUBLAS_STRSM (const char *side, const char *uplo, const char *transa, 
                   const char *diag, const int *m, const int *n, 
                   const float *alpha, const float *A, const int *lda,
                   float *B, const int *ldb)
{
    float *devPtrA, *devPtrB;
    int k;

    //  A      - REAL             array of DIMENSION ( LDA, k ), where k is m
    //           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
    //           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
    //           upper triangular part of the array  A must contain the upper
    k = (toupper(side[0]) == 'L') ? *m : *n;
    cublasAlloc (*lda * k, sizeof(devPtrA[0]), (void**)&devPtrA);
    cublasSetMatrix (imin(k,*lda), k, sizeof(A[0]), A, *lda, devPtrA, *lda);

    //  B      - REAL             array of DIMENSION ( LDB, n ).
    //           Before entry,  the leading  m by n part of the array  B must
    //           contain  the  right-hand  side  matrix  B,  and  on exit  is
    cublasAlloc ((*ldb) * (*n), sizeof(devPtrB[0]), (void**)&devPtrB);
    cublasSetMatrix (imin(*m,*ldb), *n, sizeof(B[0]), B, *ldb, devPtrB, *ldb);
    cublasStrsm (side[0], uplo[0], transa[0], diag[0], *m, *n, *alpha, devPtrA,
                 *lda, devPtrB, *ldb);
    cublasGetMatrix (imin(*m,*ldb), *n, sizeof(B[0]), devPtrB, *ldb, B, *ldb);
    cublasFree (devPtrA);
    cublasFree (devPtrB);
}

void CUBLAS_CGEMM (const char *transa, const char *transb, const int *m,
                   const int *n, const int *k, const cuComplex *alpha,
                   const cuComplex *A, const int *lda, const cuComplex *B,
                   const int *ldb, const cuComplex *beta, cuComplex *C, 
                   const int *ldc)
{
    int ka, kb;
    cuComplex *devPtrA, *devPtrB, *devPtrC;

    /*  A      - COMPLEX          array of DIMENSION ( LDA, ka ), where ka is
     *           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
     *           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
     *           part of the array  A  must contain the matrix  A,  otherwise
     *           the leading  k by m  part of the array  A  must contain  the
     *           matrix A.
     */

    ka = (toupper(transa[0]) == 'N') ? *k : *m;
    cublasAlloc (*lda * ka, sizeof(devPtrA[0]), (void**)&devPtrA);
    if (toupper(transa[0]) == 'N') {
        cublasSetMatrix (imin(*m,*lda), *k, sizeof(A[0]), A, *lda, devPtrA, 
                         *lda);
    } else {
        cublasSetMatrix (imin(*k,*lda), *m, sizeof(A[0]), A, *lda, devPtrA, 
                         *lda);
    }

    /*  B      - COMPLEX          array of DIMENSION ( LDB, kb ), where kb is
     *           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
     *           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
     *           part of the array  B  must contain the matrix  B,  otherwise
     *           the leading  n by k  part of the array  B  must contain  the
     *           matrix B.
     */

    kb = (toupper(transb[0]) == 'N') ? *n : *k;
    cublasAlloc (*ldb * kb, sizeof(devPtrB[0]), (void**)&devPtrB);
    if (toupper(transb[0]) == 'N') {
        cublasSetMatrix (imin(*k,*ldb), *n, sizeof(B[0]), B, *ldb, devPtrB, 
                         *ldb);
    } else {
        cublasSetMatrix (imin(*n,*ldb), *k, sizeof(B[0]), B, *ldb, devPtrB,
                         *ldb);
    }

    /*  C      - COMPLEX          array of DIMENSION ( LDC, n ).
     *           Before entry, the leading  m by n  part of the array  C must
     *           contain the matrix  C,  except when  beta  is zero, in which
     *           case C need not be set on entry.
     *           On exit, the array  C  is overwritten by the  m by n  matrix
     */

    cublasAlloc ((*ldc) * (*n), sizeof(devPtrC[0]), (void**)&devPtrC);
    cublasSetMatrix (imin(*m,*ldc), *n, sizeof(C[0]), C, *ldc, devPtrC, *ldc);

    cublasCgemm (transa[0], transb[0], *m, *n, *k, *alpha, devPtrA, *lda, 
                 devPtrB, *ldb, *beta, devPtrC, *ldc);

    cublasGetMatrix (imin(*m,*ldc), *n, sizeof(C[0]), devPtrC, *ldc, C, *ldc);
    cublasFree (devPtrA);
    cublasFree (devPtrB);
    cublasFree (devPtrC);
}

void CUBLAS_CHEMM (void)
{
    printf ("CUBLAS_CHEMM stub\n");
}
void CUBLAS_CSYMM (void)
{
    printf ("CUBLAS_CSYMM stub\n");
}
void CUBLAS_CTRMM (void)
{
    printf ("CUBLAS_CTRMM stub\n");
}
void CUBLAS_CTRSM (void)
{
    printf ("CUBLAS_CTRSM stub\n");
}
void CUBLAS_CHERK (void)
{
    printf ("CUBLAS_CHERK stub\n");
}
void CUBLAS_CSYRK (void)
{
    printf ("CUBLAS_CSYRK stub\n");
}
void CUBLAS_CHER2K (void)
{
    printf ("CUBLAS_CHER2K stub\n");
}
void CUBLAS_CSYR2K (void)
{
    printf ("CUBLAS_CSYR2K stub\n");
}

#else /* defined(CUBLAS_USE_THUNKING) */

/*
 * Fortran callable thin wrappers. Fortran application must allocate and
 * deallocate GPU memory, and copy data up and down.
 */
#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus */
#if CUBLAS_FORTRAN_COMPILER==CUBLAS_G77
double CUBLAS_SDOT (const int *n, const devptr_t *devPtrx, const int *incx, 
                    const devptr_t *devPtry, const int *incy);
double CUBLAS_SASUM (const int *n, const devptr_t *devPtrx, const int *incx);
double CUBLAS_SNRM2 (const int *n, const devptr_t *devPtrx, const int *incx);
double CUBLAS_SCASUM (const int *n, const devptr_t *devPtrx, const int *incx);
double CUBLAS_SCNRM2 (const int *n, const devptr_t *devPtrx, const int *incx);
#else
float CUBLAS_SDOT (const int *n, const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy);
float CUBLAS_SASUM (const int *n, const devptr_t *devPtrx, const int *incx);
float CUBLAS_SNRM2 (const int *n, const devptr_t *devPtrx, const int *incx);
float CUBLAS_SCASUM (const int *n, const devptr_t *devPtrx, const int *incx);
float CUBLAS_SCNRM2 (const int *n, const devptr_t *devPtrx, const int *incx);
#endif

int CUBLAS_ISAMAX (const int *n, const devptr_t *devPtrx, const int *incx);
int CUBLAS_ISAMIN (const int *n, const devptr_t *devPtrx, const int *incx);
void CUBLAS_SAXPY (const int *n, const float *alpha, const devptr_t *devPtrx, 
                   const int *incx, const devptr_t *devPtry, const int *incy);
void CUBLAS_SCOPY (const int *n, const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy);
void CUBLAS_SROT (const int *n, const devptr_t *devPtrX, const int *incx, 
                  const devptr_t *devPtrY, const int *incy, const float *sc, 
                  const float *ss);
void CUBLAS_SROTG (float *sa, float *sb, float *sc, float *ss);
void CUBLAS_SROTM (const int *n, const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy, 
                   const float* sparam);
void CUBLAS_SROTMG (float *sd1, float *sd2, float *sx1, const float *sy1, 
                    float* sparam);
void CUBLAS_SSCAL (const int *n, const float *alpha, const devptr_t *devPtrx,
                   const int *incx);
void CUBLAS_SSWAP (const int *n, const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy);

void CUBLAS_CAXPY (const int *n, const cuComplex *alpha,
                   const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy);
void CUBLAS_CCOPY (const int *n, const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy);
void CUBLAS_CROT (const int *n, const devptr_t *devPtrX, const int *incx, 
                  const devptr_t *devPtrY, const int *incy, const float *sc, 
                  const cuComplex *cs);
void CUBLAS_CROTG (cuComplex *ca, const cuComplex *cb, float *sc,
                   cuComplex *cs);
void CUBLAS_CSCAL (const int *n, const cuComplex *alpha, 
                   const devptr_t *devPtrx, const int *incx);
void CUBLAS_CSROT (const int *n, const devptr_t *devPtrX, const int *incx, 
                   const devptr_t *devPtrY, const int *incy, const float *sc, 
                   const float *ss);
void CUBLAS_CSSCAL (const int *n, const float *alpha, const devptr_t *devPtrx, 
                    const int *incx);
void CUBLAS_CSWAP (const int *n, const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy);
void CUBLAS_CDOTU (cuComplex *retVal, const int *n, const devptr_t *devPtrx, 
                   const int *incx, const devptr_t *devPtry, const int *incy);
void CUBLAS_CDOTC (cuComplex *retVal, const int *n, const devptr_t *devPtrx, 
                   const int *incx, const devptr_t *devPtry, const int *incy);
int CUBLAS_ICAMAX (const int *n, const devptr_t *devPtrx, const int *incx);
int CUBLAS_ICAMIN (const int *n, const devptr_t *devPtrx, const int *incx);

/* BLAS2 */
void CUBLAS_SGBMV (const char *trans, const int *m, const int *n,
                   const int *kl, const int *ku, const float *alpha, 
                   const devptr_t *devPtrA, const int *lda, 
                   const devptr_t *devPtrx, const int *incx, const float *beta,
                   const devptr_t *devPtry, const int *incy);
void CUBLAS_SGEMV (const char *trans, const int *m, const int *n,
                   const float *alpha, const devptr_t *devPtrA, const int *lda,
                   const devptr_t *devPtrx, const int *incx, const float *beta,
                   const devptr_t *devPtry, const int *incy);
void CUBLAS_SGER (const int *m, const int *n, const float *alpha, 
                  const devptr_t *devPtrx, const int *incx, 
                  const devptr_t *devPtry, const int *incy, 
                  const devptr_t *devPtrA, const int *lda);
void CUBLAS_SSBMV (const char *uplo, const int *n, const int *k, 
                   const float *alpha, const devptr_t *devPtrA, const int *lda, 
                   const devptr_t *devPtrx, const int *incx, const float *beta,
                   const devptr_t *devPtry, const int *incy);
void CUBLAS_SSPMV (const char *uplo, const int *n, const float *alpha, 
                   const devptr_t *devPtrAP, const devptr_t *devPtrx, 
                   const int *incx, const float *beta, const devptr_t *devPtry,
                   const int *incy);
void CUBLAS_SSPR (const char *uplo, const int *n, const float *alpha,
                  const devptr_t *devPtrx, const int *incx,
                  const devptr_t *devPtrAP);
void CUBLAS_SSPR2 (const char *uplo, const int *n, const float *alpha, 
                   const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy, 
                   const devptr_t *devPtrAP);
void CUBLAS_SSYMV (const char *uplo, const int *n, const float *alpha, 
                   const devptr_t *devPtrA, const int *lda, 
                   const devptr_t *devPtrx, const int *incx, const float *beta,
                   const devptr_t *devPtry, const int *incy);
void CUBLAS_SSYR (const char *uplo, const int *n, const float *alpha,
                  const devptr_t *devPtrx, const int *incx, 
                  const devptr_t *devPtrA, const int *lda);
void CUBLAS_SSYR2 (const char *uplo, const int *n, const float *alpha, 
                   const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy, 
                   const devptr_t *devPtrA, const int *lda);
void CUBLAS_STBMV (const char *uplo, const char *trans, const char *diag, 
                   const int *n, const int *k, const devptr_t *devPtrA,
                   const int *lda, const devptr_t *devPtrx, const int *incx);
void CUBLAS_STBSV (const char *uplo, const char *trans, const char *diag, 
                   const int *n, const int *k, const devptr_t *devPtrA, 
                   const int *lda, const devptr_t *devPtrx, const int *incx);
void CUBLAS_STPMV (const char *uplo, const char *trans, const char *diag, 
                   const int *n, const devptr_t *devPtrAP, 
                   const devptr_t *devPtrx, const int *incx);
void CUBLAS_STPSV (const char *uplo, const char *trans, const char *diag, 
                   const int *n, const devptr_t *devPtrAP, 
                   const devptr_t *devPtrx, const int *incx);
void CUBLAS_STRMV (const char *uplo, const char *trans, const char *diag, 
                   const int *n, const devptr_t *devPtrA, const int *lda,
                   const devptr_t *devPtrx, const int *incx);
void CUBLAS_STRSV (const char *uplo, const char *trans, const char *diag, 
                   const int *n, const devptr_t *devPtrA, const int *lda,
                   const devptr_t *devPtrx, const int *incx);

/* BLAS 3 */
void CUBLAS_SGEMM (const char *transa, const char *transb, const int *m,
                   const int *n, const int *k, const float *alpha, 
                   const devptr_t *A, const int *lda, const devptr_t *B, 
                   const int *ldb, const float *beta, const devptr_t *C, 
                   const int *ldc);
void CUBLAS_SSYMM (const char *side, const char *uplo, const int *m,
                   const int *n, const float *alpha, const devptr_t *devPtrA,
                   const int *lda, const devptr_t *devPtrB, const int *ldb,
                   const float *beta, const devptr_t *devPtrC, const int *ldc);
void CUBLAS_SSYR2K (const char *uplo, const char *trans, const int *n,
                    const int *k, const float *alpha, const devptr_t *devPtrA,
                    const int *lda, const devptr_t *devPtrB, const int *ldb,
                    const float *beta, const devptr_t *devPtrC, const int *ldc);
void CUBLAS_SSYRK (const char *uplo, const char *trans, const int *n,
                   const int *k, const float *alpha, const devptr_t *devPtrA,
                   const int *lda, const float *beta, const devptr_t *devPtrC,
                   const int *ldc);
void CUBLAS_STRMM (const char *side, const char *uplo, const char *transa, 
                   const char *diag, const int *m, const int *n,
                   const float *alpha, const devptr_t *devPtrA, const int *lda,
                   const devptr_t *devPtrB, const int *ldb);
void CUBLAS_STRSM (const char *side, const char *uplo, const char *transa, 
                   const char *diag, const int *m, const int *n, 
                   const float *alpha, const devptr_t *devPtrA, const int *lda, 
                   const devptr_t *devPtrB, const int *ldb);

void CUBLAS_CGEMM (const char *transa, const char *transb, const int *m,
                   const int *n, const int *k, const cuComplex *alpha,
                   const devptr_t *devPtrA, const int *lda,
                   const devptr_t *devPtrB, const int *ldb, 
                   const cuComplex *beta, const devptr_t *devPtrC,
                   const int *ldc);
#if defined(__cplusplus)
}
#endif /* __cplusplus */

/*---------------------------------------------------------------------------*/
/*---------------------------------- BLAS1 ----------------------------------*/
/*---------------------------------------------------------------------------*/

int CUBLAS_ISAMAX (const int *n, const devptr_t *devPtrx, const int *incx)
{
    float *x = (float *)(uintptr_t)(*devPtrx);
    int retVal;
    retVal = cublasIsamax (*n, x, *incx);
    return retVal;
}

int CUBLAS_ISAMIN (const int *n, const devptr_t *devPtrx, const int *incx)
{
    float *x = (float *)(uintptr_t)(*devPtrx);
    int retVal;
    retVal = cublasIsamin (*n, x, *incx);
    return retVal;
}

#if CUBLAS_FORTRAN_COMPILER==CUBLAS_G77
double CUBLAS_SASUM (const int *n, const devptr_t *devPtrx, const int *incx)
#else
float CUBLAS_SASUM (const int *n, const devptr_t *devPtrx, const int *incx)
#endif
{
    float *x = (float *)(uintptr_t)(*devPtrx);
    float retVal;
    retVal = cublasSasum (*n, x, *incx);
    return retVal;
}

void CUBLAS_SAXPY (const int *n, const float *alpha, const devptr_t *devPtrx, 
                   const int *incx, const devptr_t *devPtry, const int *incy)
{
    float *x = (float *)(uintptr_t)(*devPtrx);
    float *y = (float *)(uintptr_t)(*devPtry);
    cublasSaxpy (*n, *alpha, x, *incx, y, *incy);
}

void CUBLAS_SCOPY (const int *n, const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy)
{
    float *x = (float *)(uintptr_t)(*devPtrx);
    float *y = (float *)(uintptr_t)(*devPtry);
    cublasScopy (*n, x, *incx, y, *incy);
}

#if CUBLAS_FORTRAN_COMPILER==CUBLAS_G77
double CUBLAS_SDOT (const int *n, const devptr_t *devPtrx, const int *incx, 
                    const devptr_t *devPtry, const int *incy)
#else
float CUBLAS_SDOT (const int *n, const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy)
#endif
{
    float *x = (float *)(uintptr_t)(*devPtrx);
    float *y = (float *)(uintptr_t)(*devPtry);
    return cublasSdot (*n, x, *incx, y, *incy);
}

#if CUBLAS_FORTRAN_COMPILER==CUBLAS_G77
double CUBLAS_SNRM2 (const int *n, const devptr_t *devPtrx, const int *incx)
#else
float CUBLAS_SNRM2 (const int *n, const devptr_t *devPtrx, const int *incx)
#endif
{
    float *x = (float *)(uintptr_t)(*devPtrx);
    return cublasSnrm2 (*n, x, *incx);
}

void CUBLAS_SROT (const int *n, const devptr_t *devPtrx, const int *incx, 
                  const devptr_t *devPtry, const int *incy, const float *sc, 
                  const float *ss)
{
    float *x = (float *)(uintptr_t)(*devPtrx);
    float *y = (float *)(uintptr_t)(*devPtry);
    cublasSrot (*n, x, *incx, y, *incy, *sc, *ss);
}

void CUBLAS_SROTG (float *sa, float *sb, float *sc, float *ss)
{
    cublasSrotg (sa, sb, sc, ss);
}

void CUBLAS_SROTM (const int *n, const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy, 
                   const float* sparam) 
{
    float *x = (float *)(uintptr_t)(*devPtrx);
    float *y = (float *)(uintptr_t)(*devPtry);
    cublasSrotm (*n, x, *incx, y, *incy, sparam);
}

void CUBLAS_SROTMG (float *sd1, float *sd2, float *sx1, const float *sy1,
                    float* sparam)
{
    cublasSrotmg (sd1, sd2, sx1, sy1, sparam);
}

void CUBLAS_SSCAL (const int *n, const float *alpha, const devptr_t *devPtrx,
                   const int *incx)
{
    float *x = (float *)(uintptr_t)(*devPtrx);
    cublasSscal (*n, *alpha, x, *incx);
}

void CUBLAS_SSWAP (const int *n, const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy)
{
    float *x = (float *)(uintptr_t)(*devPtrx);
    float *y = (float *)(uintptr_t)(*devPtry);
    cublasSswap (*n, x, *incx, y, *incy);
}

void CUBLAS_CAXPY (const int *n, const cuComplex *alpha, 
                   const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy)
{
    cuComplex *x = (cuComplex *)(uintptr_t)(*devPtrx);
    cuComplex *y = (cuComplex *)(uintptr_t)(*devPtry);
    cublasCaxpy (*n, *alpha, x, *incx, y, *incy);
}

void CUBLAS_CCOPY (const int *n, const devptr_t *devPtrx, const int *incx,
                   const devptr_t *devPtry, const int *incy)
{
    cuComplex *x = (cuComplex *)(uintptr_t)(*devPtrx);
    cuComplex *y = (cuComplex *)(uintptr_t)(*devPtry);
    cublasCcopy (*n, x, *incx, y, *incy);
}

void CUBLAS_CROT (const int *n, const devptr_t *devPtrx, const int *incx, 
                  const devptr_t *devPtry, const int *incy, const float *sc, 
                  const cuComplex *cs)
{
    cuComplex *x = (cuComplex *)(uintptr_t)(*devPtrx);
    cuComplex *y = (cuComplex *)(uintptr_t)(*devPtry);
    cublasCrot (*n, x, *incx, y, *incy, *sc, *cs);
}

void CUBLAS_CROTG (cuComplex *ca, const cuComplex *cb, float *sc,
                   cuComplex *cs)
{
    cublasCrotg (ca, *cb, sc, cs);
}

void CUBLAS_CSCAL (const int *n, const cuComplex *alpha, 
                   const devptr_t *devPtrx, const int *incx)
{
    cuComplex *x = (cuComplex *)(uintptr_t)(*devPtrx);
    cublasCscal (*n, *alpha, x, *incx);
}

void CUBLAS_CSROT (const int *n, const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy, const float *sc, 
                   const float *ss)
{
    cuComplex *x = (cuComplex *)(uintptr_t)(*devPtrx);
    cuComplex *y = (cuComplex *)(uintptr_t)(*devPtry);
    cublasCsrot (*n, x, *incx, y, *incy, *sc, *ss);
}

void CUBLAS_CSSCAL (const int *n, const float *alpha, const devptr_t *devPtrx,
                    const int *incx)
{
    cuComplex *x = (cuComplex *)(uintptr_t)(*devPtrx);
    cublasCsscal (*n, *alpha, x, *incx);
}

void CUBLAS_CSWAP (const int *n, const devptr_t *devPtrx, const int *incx,
                   const devptr_t *devPtry, const int *incy)
{
    cuComplex *x = (cuComplex *)(uintptr_t)(*devPtrx);
    cuComplex *y = (cuComplex *)(uintptr_t)(*devPtry);
    cublasCswap (*n, x, *incx, y, *incy);
}

void CUBLAS_CDOTU (cuComplex *retVal, const int *n, const devptr_t *devPtrx,
                   const int *incx, const devptr_t *devPtry,const int *incy)
{
    cuComplex *x = (cuComplex *)(uintptr_t)(*devPtrx);
    cuComplex *y = (cuComplex *)(uintptr_t)(*devPtry);
    *retVal = cublasCdotu (*n, x, *incx, y, *incy);
}

void CUBLAS_CDOTC (cuComplex *retVal, const int *n, const devptr_t *devPtrx,
                   const int *incx, const devptr_t *devPtry, const int *incy)
{
    cuComplex *x = (cuComplex *)(uintptr_t)(*devPtrx);
    cuComplex *y = (cuComplex *)(uintptr_t)(*devPtry);
    *retVal = cublasCdotc (*n, x, *incx, y, *incy);
}

int CUBLAS_ICAMAX (const int *n, const devptr_t *devPtrx, const int *incx)
{
    cuComplex *x = (cuComplex *)(uintptr_t)(*devPtrx);
    return cublasIcamax (*n, x, *incx);
}

int CUBLAS_ICAMIN (const int *n, const devptr_t *devPtrx, const int *incx)
{
    cuComplex *x = (cuComplex *)(uintptr_t)(*devPtrx);
    return cublasIcamin (*n, x, *incx);
}

#if CUBLAS_FORTRAN_COMPILER==CUBLAS_G77
double CUBLAS_SCASUM (const int *n, const devptr_t *devPtrx, const int *incx)
#else
float CUBLAS_SCASUM (const int *n, const devptr_t *devPtrx, const int *incx)
#endif
{
    cuComplex *x = (cuComplex *)(uintptr_t)(*devPtrx);
    return cublasScasum (*n, x, *incx);
}

#if CUBLAS_FORTRAN_COMPILER==CUBLAS_G77
double CUBLAS_SCNRM2 (const int *n, const devptr_t *devPtrx, const int *incx)
#else
float CUBLAS_SCNRM2 (const int *n, const devptr_t *devPtrx, const int *incx)
#endif
{
    cuComplex *x = (cuComplex *)(uintptr_t)(*devPtrx);
    return cublasScnrm2 (*n, x, *incx);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------- BLAS2 ----------------------------------*/
/*---------------------------------------------------------------------------*/

void CUBLAS_SGBMV (const char *trans, const int *m, const int *n, 
                   const int *kl, const int *ku, const float *alpha,
                   const devptr_t *devPtrA, const int *lda, 
                   const devptr_t *devPtrx, const int *incx, const float *beta,
                   const devptr_t *devPtry, const int *incy)
{
    float *A = (float *)(uintptr_t)(*devPtrA);
    float *x = (float *)(uintptr_t)(*devPtrx);
    float *y = (float *)(uintptr_t)(*devPtry);
    cublasSgbmv (trans[0], *m, *n, *kl, *ku, *alpha, A, *lda, x, *incx, *beta,
                 y, *incy);
}

void CUBLAS_SGEMV (const char *trans, const int *m, const int *n, 
                   const float *alpha, const devptr_t *devPtrA, const int *lda,
                   const devptr_t *devPtrx, const int *incx, const float *beta,
                   const devptr_t *devPtry, const int *incy)
{
    float *A = (float *)(uintptr_t)(*devPtrA);
    float *x = (float *)(uintptr_t)(*devPtrx);
    float *y = (float *)(uintptr_t)(*devPtry);
    cublasSgemv (trans[0], *m, *n, *alpha, A, *lda, x, *incx, *beta, y, *incy);
}

void CUBLAS_SGER (const int *m, const int *n, const float *alpha, 
                  const devptr_t *devPtrx, const int *incx,
                  const devptr_t *devPtry, const int *incy,
                  const devptr_t *devPtrA, const int *lda)
{
    float *A = (float *)(uintptr_t)(*devPtrA);
    float *x = (float *)(uintptr_t)(*devPtrx);
    float *y = (float *)(uintptr_t)(*devPtry);    
    cublasSger (*m, *n, *alpha, x, *incx, y, *incy, A, *lda);
}

void CUBLAS_SSBMV (const char *uplo, const int *n, const int *k,
                   const float *alpha, const devptr_t *devPtrA, const int *lda,
                   const devptr_t *devPtrx, const int *incx, const float *beta,
                   const devptr_t *devPtry, const int *incy)
{
    float *A = (float *)(uintptr_t)(*devPtrA);
    float *x = (float *)(uintptr_t)(*devPtrx);
    float *y = (float *)(uintptr_t)(*devPtry);    
    cublasSsbmv (uplo[0], *n, *k, *alpha, A, *lda, x, *incx, *beta, y, *incy);
}

void CUBLAS_SSPMV (const char *uplo, const int *n, const float *alpha,
                   const devptr_t *devPtrAP, const devptr_t *devPtrx,
                   const int *incx, const float *beta, const devptr_t *devPtry,
                   const int *incy)
{
    float *AP = (float *)(uintptr_t)(*devPtrAP);
    float *x = (float *)(uintptr_t)(*devPtrx);
    float *y = (float *)(uintptr_t)(*devPtry);    
    cublasSspmv (uplo[0], *n, *alpha, AP, x, *incx, *beta, y, *incy);
}

void CUBLAS_SSPR (const char *uplo, const int *n, const float *alpha, 
                  const devptr_t *devPtrx, const int *incx,
                  const devptr_t *devPtrAP)
{
    float *AP = (float *)(uintptr_t)(*devPtrAP);
    float *x = (float *)(uintptr_t)(*devPtrx);
    cublasSspr (uplo[0], *n, *alpha, x, *incx, AP);
}

void CUBLAS_SSPR2 (const char *uplo, const int *n, const float *alpha,
                   const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy,
                   const devptr_t *devPtrAP)
{
    float *AP = (float *)(uintptr_t)(*devPtrAP);
    float *x = (float *)(uintptr_t)(*devPtrx);
    float *y = (float *)(uintptr_t)(*devPtry);    
    cublasSspr2 (uplo[0], *n, *alpha, x, *incx, y, *incy, AP);
}

void CUBLAS_SSYMV (const char *uplo, const int *n, const float *alpha,
                   const devptr_t *devPtrA, const int *lda, 
                   const devptr_t *devPtrx, const int *incx, const float *beta,
                   const devptr_t *devPtry,
                   const int *incy)
{
    float *A = (float *)(uintptr_t)(*devPtrA);
    float *x = (float *)(uintptr_t)(*devPtrx);
    float *y = (float *)(uintptr_t)(*devPtry);    
    cublasSsymv (uplo[0], *n, *alpha, A, *lda, x, *incx, *beta, y, *incy);
}

void CUBLAS_SSYR (const char *uplo, const int *n, const float *alpha,
                  const devptr_t *devPtrx, const int *incx,
                  const devptr_t *devPtrA, const int *lda)
{
    float *A = (float *)(uintptr_t)(*devPtrA);
    float *x = (float *)(uintptr_t)(*devPtrx);    
    cublasSsyr (uplo[0], *n, *alpha, x, *incx, A, *lda);
}

void CUBLAS_SSYR2 (const char *uplo, const int *n, const float *alpha,
                   const devptr_t *devPtrx, const int *incx,
                   const devptr_t *devPtry, const int *incy, 
                   const devptr_t *devPtrA, const int *lda)
{
    float *A = (float *)(uintptr_t)(*devPtrA);
    float *x = (float *)(uintptr_t)(*devPtrx);
    float *y = (float *)(uintptr_t)(*devPtry);    
    cublasSsyr2 (uplo[0], *n, *alpha, x, *incx, y, *incy, A, *lda);
}

void CUBLAS_STBMV (const char *uplo, const char *trans, const char *diag,
                   const int *n, const int *k, const devptr_t *devPtrA, 
                   const int *lda, const devptr_t *devPtrx, const int *incx)
{
    float *A = (float *)(uintptr_t)(*devPtrA);
    float *x = (float *)(uintptr_t)(*devPtrx);    
    cublasStbmv (uplo[0], trans[0], diag[0], *n, *k, A, *lda, x, *incx);
}

void CUBLAS_STBSV (const char *uplo, const char *trans, const char *diag,
                   const int *n, const int *k, const devptr_t *devPtrA, 
                   const int *lda, const devptr_t *devPtrx, const int *incx)
{
    float *A = (float *)(uintptr_t)(*devPtrA);
    float *x = (float *)(uintptr_t)(*devPtrx);       
    cublasStbsv (uplo[0], trans[0], diag[0], *n, *k, A, *lda, x, *incx);
}

void CUBLAS_STPMV (const char *uplo, const char *trans, const char *diag,
                   const int *n,  const devptr_t *devPtrAP, 
                   const devptr_t *devPtrx, const int *incx)
{
    float *AP = (float *)(uintptr_t)(*devPtrAP);
    float *x = (float *)(uintptr_t)(*devPtrx);       
    cublasStpmv (uplo[0], trans[0], diag[0], *n, AP, x, *incx);
}

void CUBLAS_STPSV (const char *uplo, const char *trans, const char *diag,
                   const int *n, const devptr_t *devPtrAP, 
                   const devptr_t *devPtrx, const int *incx)
{
    float *AP = (float *)(uintptr_t)(*devPtrAP);
    float *x = (float *)(uintptr_t)(*devPtrx);       
    cublasStpsv (uplo[0], trans[0], diag[0], *n, AP, x, *incx);
}

void CUBLAS_STRMV (const char *uplo, const char *trans, const char *diag,
                   const int *n, const devptr_t *devPtrA, const int *lda,
                   const devptr_t *devPtrx, const int *incx)
{
    float *A = (float *)(uintptr_t)(*devPtrA);
    float *x = (float *)(uintptr_t)(*devPtrx);       
    cublasStrmv (uplo[0], trans[0], diag[0], *n, A, *lda, x, *incx);
}

void CUBLAS_STRSV (const char *uplo, const char *trans, const char *diag,
                   const int *n, const devptr_t *devPtrA, const int *lda,
                   const devptr_t *devPtrx, const int *incx)
{
    float *A = (float *)(uintptr_t)(*devPtrA);
    float *x = (float *)(uintptr_t)(*devPtrx);       
    cublasStrsv (uplo[0], trans[0], diag[0], *n, A, *lda, x, *incx);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------- BLAS3 ----------------------------------*/
/*---------------------------------------------------------------------------*/

void CUBLAS_SGEMM (const char *transa, const char *transb, const int *m,
                   const int *n, const int *k, const float *alpha,
                   const devptr_t *devPtrA, const int *lda, 
                   const devptr_t *devPtrB, const int *ldb, const float *beta,
                   const devptr_t *devPtrC, const int *ldc)
{
    float *A = (float *)(uintptr_t)(*devPtrA);
    float *B = (float *)(uintptr_t)(*devPtrB);
    float *C = (float *)(uintptr_t)(*devPtrC);
    cublasSgemm (transa[0], transb[0], *m, *n, *k, *alpha, A, *lda, 
                 B, *ldb, *beta, C, *ldc);
}

void CUBLAS_SSYMM (const char *side, const char *uplo, const int *m, 
                   const int *n, const float *alpha, const devptr_t *devPtrA,
                   const int *lda, const devptr_t *devPtrB, const int *ldb, 
                   const float *beta, const devptr_t *devPtrC, const int *ldc)
{
    float *A = (float *)(uintptr_t)(*devPtrA);
    float *B = (float *)(uintptr_t)(*devPtrB);
    float *C = (float *)(uintptr_t)(*devPtrC);
    cublasSsymm (*side, *uplo, *m, *m, *alpha, A, *lda, B, *ldb, *beta, C,
                 *ldc);
}

void CUBLAS_SSYR2K (const char *uplo, const char *trans, const int *n,
                    const int *k, const float *alpha, const devptr_t *devPtrA,
                    const int *lda, const devptr_t *devPtrB, const int *ldb, 
                    const float *beta, const devptr_t *devPtrC, const int *ldc)
{
    float *A = (float *)(uintptr_t)(*devPtrA);
    float *B = (float *)(uintptr_t)(*devPtrB);
    float *C = (float *)(uintptr_t)(*devPtrC);
    cublasSsyr2k (*uplo, *trans, *n, *k, *alpha, A, *lda, B, *ldb, *beta, 
                  C, *ldc);
}

void CUBLAS_SSYRK (const char *uplo, const char *trans, const int *n, 
                   const int *k, const float *alpha, const devptr_t *devPtrA, 
                   const int *lda, const float *beta, const devptr_t *devPtrC,
                   const int *ldc)
{
    float *A = (float *)(uintptr_t)(*devPtrA);
    float *C = (float *)(uintptr_t)(*devPtrC);
    cublasSsyrk (*uplo, *trans, *n, *k, *alpha, A, *lda, *beta, C, *ldc);
}

void CUBLAS_STRMM (const char *side, const char *uplo, const char *transa,
                   const char *diag, const int *m, const int *n,
                   const float *alpha, const devptr_t *devPtrA, const int *lda,
                   const devptr_t *devPtrB, const int *ldb)
{
    float *A = (float *)(uintptr_t)(*devPtrA);
    float *B = (float *)(uintptr_t)(*devPtrB);
    cublasStrmm (*side, *uplo, *transa, *diag, *m, *n, *alpha, A, *lda, B,
                 *ldb);
}

void CUBLAS_STRSM (const char *side, const char *uplo, const char *transa,
                   const char *diag, const int *m, const int *n, 
                   const float *alpha, const devptr_t *devPtrA, const int *lda,
                   const devptr_t *devPtrB, const int *ldb)
{
    float *A = (float *)(uintptr_t)*devPtrA;
    float *B = (float *)(uintptr_t)*devPtrB;
    cublasStrsm (side[0], uplo[0], transa[0], diag[0], *m, *n, *alpha,
                 A, *lda, B, *ldb);
}

void CUBLAS_CGEMM (const char *transa, const char *transb, const int *m,
                   const int *n, const int *k, const cuComplex *alpha,
                   const devptr_t *devPtrA, const int *lda,
                   const devptr_t *devPtrB, const int *ldb, 
                   const cuComplex *beta, const devptr_t *devPtrC,
                   const int *ldc)
{
    cuComplex *A = (cuComplex *)(uintptr_t)*devPtrA;
    cuComplex *B = (cuComplex *)(uintptr_t)*devPtrB;
    cuComplex *C = (cuComplex *)(uintptr_t)*devPtrC;    
    cublasCgemm (transa[0], transb[0], *m, *n, *k, *alpha, A, *lda, B, *ldb, 
                 *beta, C, *ldc);
}

void CUBLAS_CHEMM (void)
{
    printf ("CUBLAS_CHEMM stub\n");
}
void CUBLAS_CSYMM (void)
{
    printf ("CUBLAS_CSYMM stub\n");
}
void CUBLAS_CTRMM (void)
{
    printf ("CUBLAS_CTRMM stub\n");
}
void CUBLAS_CTRSM (void)
{
    printf ("CUBLAS_CTRSM stub\n");
}
void CUBLAS_CHERK (void)
{
    printf ("CUBLAS_CHERK stub\n");
}
void CUBLAS_CSYRK (void)
{
    printf ("CUBLAS_CSYRK stub\n");
}
void CUBLAS_CHER2K (void)
{
    printf ("CUBLAS_CHER2K stub\n");
}
void CUBLAS_CSYR2K (void)
{
    printf ("CUBLAS_CSYR2K stub\n");
}

#endif  /* defined(CUBLAS_USE_THUNKING) */

