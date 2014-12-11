/* abi_common.h */

/*
 * Copyright (C) 2008-2014 ABINIT Group (MG)
 *
 * This file is part of the ABINIT software package. For license information,
 * please see the COPYING file in the top-level directory of the ABINIT source
 * distribution.
 *
 */

#ifndef _ABINIT_COMMON_H
#define _ABINIT_COMMON_H

/*
 * Language standards requires the existance of pre-defined macros
 * Microsoft Visual C++ does not define __STDC__,
 * Sun Workshop 4.2 supports C94 without setting __STDC_VERSION__ to the proper value
 */

#if defined (__STDC__)
# define PREDEF_STANDARD_C_1989    /** ANSI X3.159-1989 **/
# if defined (__STDC_VERSION__)
#  define PREDEF_STANDARD_C_1990   /** ISO/IEC 9899:1990 **/
#  if (__STDC_VERSION__ >= 199409L)
#   define PREDEF_STANDARD_C_1994  /** ISO/IEC 9899-1:1994 **/
#  endif
#  if (__STDC_VERSION__ >= 199901L)
#   define PREDEF_STANDARD_C_1999  /** ISO/IEC 9899:1999 **/
#  endif
# endif
#endif

/** #define DEBUG_MODE   **/

#if defined(HAVE_FC_LONG_LINES) || defined(__INTEL_COMPILER) || defined(FC_NAG) || !defined(HAVE_FC_MACRO_NEWLINE)
# define NEWLINE ;
#else
# define NEWLINE \newline
#endif

/**
#ifdef FC_NAG
#  undef HAVE_FC_LONG_LINES
#  define HAVE_FC_LONG_LINES
#endif
**/


/** this does not work with gfort, pgi, **/
#if defined (FC_GNU) || defined(FC_G95) || defined (FC_PGI)
#define QUOTE(x)     'x'
#else
#define QUOTE(x)     #x
#endif


/** Token concatenation
 1) ANSI CPP
 2) Traditional CPP
**/

#if defined (FC_INTEL)
#define CONCAT(x,y) x ## y
#else
#define CONCAT(x,y) x/**/y
#endif

#define BYTE_SIZE(array)  PRODUCT(SHAPE(array)) * DBLE(KIND(array))

/*
 * ABI_  abinit macros.
 * DBG_  macros for debugging. Defined only if abinit is compiled in DEBUG_MODE.
 * MSG_  macros for logging.
*/

#ifdef HAVE_FC_LONG_LINES
#  define ABI_CHECK(expr,str) if (.not.(expr)) call assert(.FALSE.,str,__FILE__,__LINE__)
#else
#  define ABI_CHECK(expr,str) if (.not.(expr)) call assert(.FALSE.,str)
#endif

#if defined HAVE_FC_LONG_LINES
#  define ABI_CHECK_MPI(ierr,msg)      call check_mpi_ierr(ierr,msg,"PERS",__FILE__,__LINE__)
#  define ABI_CHECK_MPI_PERS(ierr,msg) call check_mpi_ierr(ierr,msg,"PERS",__FILE__,__LINE__)
#else
#  define ABI_CHECK_MPI(ierr,msg) call check_mpi_ierr(ierr,msg,"PERS")
#  define ABI_CHECK_MPI_PERS(ierr,msg) call check_mpi_ierr(ierr,msg,"PERS")
#endif

/* This macro is used in low-level MPI wrappers to return mpierr to the caller
 * so that one can locate the problematic call. Use it wisely since it may cause memory leaks.
 * #define ABI_HANDLE_MPIERR(mpierr) ABI_CHECK_MPI(mpierr,"ABI_HANDLE: Fatal error")
 */
#define ABI_HANDLE_MPIERR(mpierr) if (mpierr/=MPI_SUCCESS) RETURN

/* Macros for memory checking and profiling
 *
 *   ABI_ALLOCATE   --> Allocate memory for intrinsic datatypes (real, integer, complex).
 *   ABI_CALLOC     --> Clear alloc: same as ABI_ALLOCATE but initializes memory with zeros
 *   ABI_DEALLOCATE --> Free memory allocated by ABI_ALLOCATE
 *
 * To allocate/deallocate arrays of user-defined type use:
 *
 *   ABI_DT_ALLOCATE
 *   ABI_DT_DEALLOCATE
 *
 * Note for programmers using OpenMP:
 *   Do not use stat=ABI_ALLOC_STAT_ABI inside an OpenMP parallel region
 *   ABI_ALLOC_STAT_ABI is a global variable defined in m_profiling_abi.F90
 *   and this can have a detrimental effect if the allocation is done inside an OpenMP parallel region.
 *
 #define HAVE_MEM_PROFILING
*/

#ifdef HAVE_MEM_PROFILING

#  define ABI_ALLOCATE(ARR,SIZE) \
   allocate(ARR SIZE,stat=ABI_ALLOC_STAT_ABI) NEWLINE \
   call memocc_abi(ABI_ALLOC_STAT_ABI,product(int(shape(ARR),kind=8))*kind(ARR),QUOTE(ARR),ABI_FUNC)

#  define ABI_CALLOC(ARR,SIZE) \
   allocate(ARR SIZE,stat=ABI_ALLOC_STAT_ABI) NEWLINE \
   call memocc_abi(ABI_ALLOC_STAT_ABI,product(int(shape(ARR),kind=8))*kind(ARR),QUOTE(ARR),ABI_FUNC) NEWLINE ARR = zero

#  define ABI_DEALLOCATE(ARR) \
   ABI_ALLOC_SIZE_ABI=-product(int(shape(ARR),kind=8))*kind(ARR) NEWLINE \
   deallocate(ARR,stat=ABI_ALLOC_STAT_ABI) NEWLINE \
   call memocc_abi(ABI_ALLOC_STAT_ABI,ABI_ALLOC_SIZE_ABI,QUOTE(ARR),ABI_FUNC)

#else
#  ifdef HAVE_MY_PROFILE
#    define ABI_ALLOCATE(ARR,SIZE)  allocate(ARR SIZE, stat=ABI_ALLOC_STAT_ABI) NEWLINE \
        write(HAVE_MY_PROFILE,*) "+1  ",__FILE__,__LINE__, product(int(shape(ARR),kind=8))*kind(ARR), "ARR"
#    define ABI_DEALLOCATE(ARR)   write(HAVE_MY_PROFILE,*) "-1  ",__FILE__,__LINE__, product(int(shape(ARR),kind=8))*kind(ARR),"ARR" NEWLINE \
  deallocate(ARR, stat=ABI_ALLOC_STAT_ABI)
#    define ABI_CALLOC(ARR,SIZE) ABI_ALLOCATE(ARR,SIZE) NEWLINE ARR = zero
#  else
#    define ABI_ALLOCATE(ARR,SIZE) allocate(ARR SIZE, stat=ABI_ALLOC_STAT_ABI)
#    define ABI_CALLOC(ARR,SIZE) allocate(ARR SIZE, stat=ABI_ALLOC_STAT_ABI) NEWLINE ARR = zero
#    define ABI_DEALLOCATE(ARR)  deallocate(ARR, stat=ABI_ALLOC_STAT_ABI)
#  endif
#endif

#define ABI_DATATYPE_ALLOCATE(ARR,SIZE)  allocate(ARR SIZE, stat=ABI_ALLOC_STAT_ABI)
#define ABI_DATATYPE_DEALLOCATE(ARR)   deallocate(ARR, stat=ABI_ALLOC_STAT_ABI)

/* Shorthand versions for lazy programmers that are familiar with C */
#define ABI_MALLOC(ARR,SIZE) ABI_ALLOCATE(ARR,SIZE)
#define ABI_FREE(ARR) ABI_DEALLOCATE(ARR)

#define ABI_DT_MALLOC(ARR,SIZE)  ABI_DATATYPE_ALLOCATE(ARR,SIZE)
#define ABI_DT_FREE(ARR)         ABI_DATATYPE_DEALLOCATE(ARR)

/*
 * Macros to allocate/deallocate depending on the status of the entity
 * Caveat: pointers must use ABI_PTR_FREE_IF
 *
#define ABI_MALLOC_IFNOT(ARR) if (.not.allocated(ARR)) ABI_MALLOC(ARR)
#define ABI_FREE_IF(ARR) if (allocated(ARR)) ABI_FREE(ARR)
#define ABI_PTR_FREE_IF(PTR) if (associated(PTR)) ABI_FREE(PTR)
*/

#ifdef DEBUG_MODE
/* Macros used in debug mode */
#  ifdef HAVE_FC_LONG_LINES
#    define ASSERT(expr) if (.not.expr) call assert((expr), "Assertion failed", __FILE__,__LINE__)
#    define ASSERT_IF(condition, expr) if (condition) call assert((expr), "Assertion failed", __FILE__,__LINE__)
#    define DBG_CHECK(expr,str) if (.not.expr) call assert((expr), str,__FILE__,__LINE__)
#    define DBG_CHKPT(value) write(std_out,*)__FILE__,":",__LINE__,":",value
#    define DBG_ENTER(mode) call sentinel(1,mode,__FILE__,ABI_FUNC,__LINE__)
#    define DBG_EXIT(mode)  call sentinel(2,mode,__FILE__,ABI_FUNC,__LINE__)
#  else
#    define ASSERT(expr) if (.not.expr) call assert((expr), "Assertion failed")
#    define ASSERT_IF(condition, expr) if (condition) call assert((expr), "Assertion failed")
#    define DBG_CHECK(expr,str) if (.not.expr) call assert((expr),str)
#    define DBG_CHKPT(value) write(std_out,*)value
#    define DBG_ENTER(mode) call sentinel(1,mode)
#    define DBG_EXIT(mode)  call sentinel(2,mode)
#  endif

/* Stop if two arrays have different shape */
#  define DBG_EQSHAPE(arr1, arr2) if (any(shape(arr1)/=shape(arr2))) MSG_ERROR("Different shape")

#else
/* nops */
#  define ASSERT(expr)
#  define ASSERT_IF(condition, expr)
#  define DBG_CHECK(expr,str)
#  define DBG_CHKPT(value)
#  define DBG_ENTER(mode)
#  define DBG_EXIT(mode)
#  define DBG_EQSHAPE(arr1, arr2)
#endif

/* Macro for basic messages */
#ifdef HAVE_FC_LONG_LINES

#  define MSG_COMMENT(msg) call msg_hndl(msg,"COMMENT","PERS",__FILE__,__LINE__)
#  define MSG_WARNING(msg) call msg_hndl(msg,"WARNING","PERS",__FILE__,__LINE__)
#  define MSG_ERROR(msg)   call msg_hndl(msg,"ERROR"  ,"PERS",__FILE__,__LINE__)
#  define MSG_BUG(msg)     call msg_hndl(msg,"BUG"    ,"PERS",__FILE__,__LINE__)

#  define MSG_ERROR_NODUMP(msg) call msg_hndl(msg,"ERROR","PERS",__FILE__,__LINE__,NODUMP=.TRUE.)
#  define MSG_ERROR_NOSTOP(msg,ierr) \
   ierr=ierr+1;call msg_hndl(msg,"ERROR","PERS",__FILE__,__LINE__,NOSTOP=.TRUE.)

#  define ETSF_CHECK_ERROR(lstat,Error_data)   if (.not. lstat) call abietsf_msg_hndl(lstat,Error_data,"PERS",__FILE__,__LINE__)
#  define ETSF_WARN(lstat,Error_data) call abietsf_warn(lstat,Error_data,"PERS",__FILE__,__LINE__)

#  define NCF_CHECK(ncerr,msg) if (ncerr/=nf90_noerr) call netcdf_check(ncerr,msg,__FILE__,__LINE__)
#else
/*
 * Safe macros for emergency cases!
 * Useful if __FILE__ expands to the full path name exceeding
 * the max number of Fortran columns. ISO doesn't define any standard!
 */
#  define MSG_COMMENT(msg) call msg_hndl(msg,"COMMENT","PERS")
#  define MSG_WARNING(msg) call msg_hndl(msg,"WARNING","PERS")
#  define MSG_ERROR(msg)   call msg_hndl(msg,"ERROR"  ,"PERS")
#  define MSG_BUG(msg)     call msg_hndl(msg,"BUG"    ,"PERS")

#  define MSG_ERROR_NODUMP(msg) call msg_hndl(msg,"ERROR","PERS",NODUMP=.TRUE.)
#  define MSG_ERROR_NOSTOP(msg,ierr) \
   ierr=ierr+1;call msg_hndl(msg,"ERROR","PERS",NOSTOP=.TRUE.)

#  define ETSF_CHECK_ERROR(lstat,Error_data)   call abietsf_msg_hndl(lstat,Error_data,"PERS")
#  define ETSF_WARN(lstat,Error_data) call abietsf_warn(lstat,Error_data,"PERS")

#  define NCF_CHECK(ncerr,msg) call netcdf_check(ncerr,msg)
#endif

/* Macro for clean exit */
/*
#define ABI_EXIT(exit_status) call leave_new("COLL",exit_status=exit_status,print_config=.False.)
*/

/* Macro for checking whether allocation was successful */
#define ABI_CHECK_ALLOC(msg) if (ABI_ALLOC_STAT_ABI/=0) MSG_ERROR(msg)

/* Macros used for stopping the code if external libraries have not been enabled */
#ifdef HAVE_FC_LONG_LINES
#define BIGDFT_NOTENABLED_ERROR() call bigdft_lib_error(__FILE__, __LINE__)
#else
#define BIGDFT_NOTENABLED_ERROR() call bigdft_lib_error()
#endif

/* Write a warning if msg is not empty */
#define MSG_WARNING_IF(msg) if (len_trim(msg)/=0) MSG_WARNING(msg)

/* Dummy use of unused arguments to silence compiler warnings */
#define ABI_UNUSED(var) if (.FALSE.) call unused_var(var)

#ifdef HAVE_TIMER_PAPI
#ifdef HAVE_FC_LONG_LINES
#  define XPAPI_CHECK(check,msg) if (check/=PAPI_OK) call xpapi_handle_error(check,msg,__FILE__,__LINE__)
#else
#  define XPAPI_CHECK(check,msg) if (check/=PAPI_OK) call xpapi_handle_error(check,msg)
#endif
#else
#  define XPAPI_CHECK(check,msg)
#endif

/* Portable support for /dev/null */
#if defined HAVE_OS_WINDOWS
#define NULL_FILE "NUL"
#else
#define NULL_FILE "/dev/null"
#endif

/* Macros for Fortran IO (requires m_io_tools module) */
#ifdef DEBUG_MODE
/* #  define ABI_FCLOSE(fort_unit, msg) if (close_unit(fort_unit, msg)/=0) MSG_ERROR(msg) */
#  define ABI_FCLOSE(fort_unit, msg) close(fort_unit)
#else
/* #  define ABI_FCLOSE(fort_unit, msg) if (close_unit(fort_unit, msg)/=0) MSG_WARNING(msg) */
#  define ABI_FCLOSE(fort_unit, msg) close(fort_unit)
#endif

/* F2003 support  */
#define ABI_CHECK_CNULL(cptr,msg) if (.not.C_ASSOCIATED(cptr)) MSG_ERROR(msg)

#ifdef HAVE_FC_ASYNC
#define ABI_ASYNC asynchronous,
#else
#define ABI_ASYNC 
#endif

#ifdef HAVE_FC_PRIVATE
#define ABI_PRIVATE ,private
#else
#define ABI_PRIVATE
#endif

#ifdef HAVE_FC_PROTECTED
#define ABI_PROTECTED ,protected
#else
#define ABI_PROTECTED
#endif


/* F2008 support  */
#ifdef HAVE_FC_CONTIGUOUS
#define ABI_CONTIGUOUS contiguous,
#else
#define ABI_CONTIGUOUS
#endif

/* 
 * Temporary trick used to pass contiguous array descriptors to F90 routines 
 * Mainly used to avoid copy-in copy-out in m_abi_linalg
#define DEV_CONTARRD contiguous,
*/
#define DEV_CONTARRD

/* OpenMP support */
#ifndef HAVE_OMP_COLLAPSE
#define COLLAPSE(x)
#endif

/* DFTI macros (should be declared in m_dfti but build-sys tests complain */
#ifdef HAVE_FC_LONG_LINES
#  define DFTI_CHECK(status) if (status/=0) call dfti_check_status(status,__FILE__,__LINE__)
#else
#  define DFTI_CHECK(status) if (status/=0) call dfti_check_status(status)
#endif


/* Macros used in the GW code */
#ifdef HAVE_GW_DPC
#  define GWPC_CONJG(cvar)  DCONJG(cvar)
#  define GWPC_CMPLX(re,im) DCMPLX(re,im)
#else
#  define GWPC_CONJG(cvar)  CONJG(cvar)
#  define GWPC_CMPLX(re,im) CMPLX(re,im)
#endif

/*
#define yaml_out ab_out
*/
#define yaml_out std_out


/* Macros used for particular slaves of the Abinit tests farm */

/* 
 * IBM6 has a very problematic IO. Flushing the IO buffer helps avoid segmentation fault 
 * To disable these tricks, comment the three lines below.
 * */
#ifdef FC_IBM
#define HAVE_IBM6
#endif
#undef HAVE_IBM6

#ifdef HAVE_IBM6
#define _IBM6(message) call wrtout(std_out,message,"COLL",do_flush=.True.)
#else
#define _IBM6(message)
#endif

#endif
/* _ABINIT_COMMON_H */
