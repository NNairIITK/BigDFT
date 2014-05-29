#ifndef _ABINIT_COMMON_H
#define _ABINIT_COMMON_H

#if defined (__STDC__)
# define PREDEF_STANDARD_C_1989
# if defined (__STDC_VERSION__)
#  define PREDEF_STANDARD_C_1990
#  if (__STDC_VERSION__ >= 199409L)
#   define PREDEF_STANDARD_C_1994
#  endif
#  if (__STDC_VERSION__ >= 199901L)
#   define PREDEF_STANDARD_C_1999
#  endif
# endif
#endif

#if defined(HAVE_FC_LONG_LINES) || defined(__INTEL_COMPILER) || defined(FC_NAG) || !defined(HAVE_FC_MACRO_NEWLINE)
# define NEWLINE ;
#else
# define NEWLINE \newline
#endif

#if defined (FC_GNU) || defined(FC_G95) || defined (FC_PGI)
#define QUOTE(x)     'x'
#else
#define QUOTE(x)     #x
#endif


#if defined (FC_INTEL)
#define CONCAT(x,y) x ## y
#else
#define CONCAT(x,y) x/**/y
#endif

#define BYTE_SIZE(array)  PRODUCT(SHAPE(array)) * DBLE(KIND(array))

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

#define ABI_HANDLE_MPIERR(mpierr) if (mpierr/=MPI_SUCCESS) RETURN

#ifdef HAVE_MEM_PROFILING

#  define ABI_ALLOCATE(ARR,SIZE) \
   allocate(ARR SIZE,stat=ABI_ALLOC_STAT) NEWLINE \
   call memocc(ABI_ALLOC_STAT,product(shape(ARR))*kind(ARR),QUOTE(ARR),ABI_FUNC)

#  define ABI_CALLOC(ARR,SIZE) \
   allocate(ARR SIZE,stat=ABI_ALLOC_STAT) NEWLINE \
   call memocc(ABI_ALLOC_STAT,product(shape(ARR))*kind(ARR),QUOTE(ARR),ABI_FUNC) NEWLINE ARR = zero

#  define ABI_DEALLOCATE(ARR) \
   ABI_ALLOC_SIZE=-product(shape(ARR))*kind(ARR) NEWLINE \
   deallocate(ARR,stat=ABI_ALLOC_STAT) NEWLINE \
   call memocc(ABI_ALLOC_STAT,ABI_ALLOC_SIZE,QUOTE(ARR),ABI_FUNC)

#else
#  ifdef HAVE_MY_PROFILE
#    define ABI_ALLOCATE(ARR,SIZE)  allocate(ARR SIZE, stat=ABI_ALLOC_STAT) NEWLINE \
        write(HAVE_MY_PROFILE,*) "+1  ",__FILE__,__LINE__, product(shape(ARR))*kind(ARR), "ARR"
#    define ABI_DEALLOCATE(ARR)   write(HAVE_MY_PROFILE,*) "-1  ",__FILE__,__LINE__, product(shape(ARR))*kind(ARR),"ARR" NEWLINE \
  deallocate(ARR, stat=ABI_ALLOC_STAT)
#    define ABI_CALLOC(ARR,SIZE) ABI_ALLOCATE(ARR,SIZE) NEWLINE ARR = zero
#  else
#    define ABI_ALLOCATE(ARR,SIZE) allocate(ARR SIZE, stat=ABI_ALLOC_STAT)
#    define ABI_CALLOC(ARR,SIZE) allocate(ARR SIZE, stat=ABI_ALLOC_STAT) NEWLINE ARR = zero
#    define ABI_DEALLOCATE(ARR)  deallocate(ARR, stat=ABI_ALLOC_STAT)
#  endif
#endif

#define ABI_DATATYPE_ALLOCATE(ARR,SIZE)  allocate(ARR SIZE, stat=ABI_ALLOC_STAT)
#define ABI_DATATYPE_DEALLOCATE(ARR)   deallocate(ARR, stat=ABI_ALLOC_STAT)

#define ABI_MALLOC(ARR,SIZE) ABI_ALLOCATE(ARR,SIZE)
#define ABI_FREE(ARR) ABI_DEALLOCATE(ARR)

#define ABI_DT_MALLOC(ARR,SIZE)  ABI_DATATYPE_ALLOCATE(ARR,SIZE)
#define ABI_DT_FREE(ARR)         ABI_DATATYPE_DEALLOCATE(ARR)

#ifdef DEBUG_MODE
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

#  define DBG_EQSHAPE(arr1, arr2) if (any(shape(arr1)/=shape(arr2))) MSG_ERROR("Different shape")

#else
#  define ASSERT(expr)
#  define ASSERT_IF(condition, expr)
#  define DBG_CHECK(expr,str)
#  define DBG_CHKPT(value)
#  define DBG_ENTER(mode)
#  define DBG_EXIT(mode)
#  define DBG_EQSHAPE(arr1, arr2)
#endif

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

#define ABI_CHECK_ALLOC(msg) if (ABI_ALLOC_STAT/=0) MSG_ERROR(msg)

#ifdef HAVE_FC_LONG_LINES
#define BIGDFT_NOTENABLED_ERROR() call bigdft_lib_error(__FILE__, __LINE__)
#else
#define BIGDFT_NOTENABLED_ERROR() call bigdft_lib_error()
#endif

#define MSG_WARNING_IF(msg) if (len_trim(msg)/=0) MSG_WARNING(msg)

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

#if defined HAVE_OS_WINDOWS
#define NULL_FILE "NUL"
#else
#define NULL_FILE "/dev/null"
#endif

#ifdef DEBUG_MODE
#  define ABI_FCLOSE(fort_unit, msg) close(fort_unit)
#else
#  define ABI_FCLOSE(fort_unit, msg) close(fort_unit)
#endif

#define ABI_CHECK_CNULL(cptr,msg) if (.not.C_ASSOCIATED(cptr)) MSG_ERROR(msg)

#ifdef HAVE_FC_ASYNC
#define ABI_ASYNC asynchronous,
#else
#define ABI_ASYNC 
#endif


#ifdef HAVE_FC_CONTIGUOUS
#define ABI_CONTIGUOUS contiguous,
#else
#define ABI_CONTIGUOUS
#endif


#ifndef HAVE_OMP_COLLAPSE
#define COLLAPSE(x)
#endif

#ifdef HAVE_FC_LONG_LINES
#  define DFTI_CHECK(status) if (status/=0) call dfti_check_status(status,__FILE__,__LINE__)
#else
#  define DFTI_CHECK(status) if (status/=0) call dfti_check_status(status)
#endif


#ifdef HAVE_GW_DPC
#  define GWPC_CONJG(cvar)  DCONJG(cvar)
#  define GWPC_CMPLX(re,im) CMPLX(re,im)
#else
#  define GWPC_CONJG(cvar)  CONJG(cvar)
#  define GWPC_CMPLX(re,im) DCMPLX(re,im)
#endif

#define yaml_out std_out


#ifdef HAVE_IBM
#define _IBM6(message) call wrtout(std_out,message,"COLL",do_flush=.True.)
#else
#define _IBM6(message)
#endif


#endif
