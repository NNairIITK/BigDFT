

/*
 * This file is part of the libPAW library.
 * It has to be customized according to the host code.
 * For the time being there are 2 known host codes:
 * ABINIT (www.abinit.org) and BigDFT (bigdft.org).
 */

/*
 * Copyright (C) 2014-2014 ABINIT Group (MT)
 * This file is part of the ABINIT software package. For license information,
 * please see the COPYING file in the top-level directory of the ABINIT source
 * distribution.
 */


#if defined HAVE_CONFIG_H
#include "config.h"
#endif

/* =============================
 * ========= ABINIT ============
 * ============================= */
#if defined HAVE_LIBPAW_ABINIT


#  include "abi_common.h"


#  define USE_DEFS use defs_basis


#  define USE_MPI_WRAPPERS use m_xmpi



#  define USE_MSG_HANDLING use m_errors, only : msg_hndl
#  undef  HAVE_YAML


#  define USE_MEMORY_PROFILING use m_profiling_abi

#  define LIBPAW_ALLOCATE(ARR,SIZE) ABI_ALLOCATE(ARR,SIZE)
#  define LIBPAW_DEALLOCATE(ARR) ABI_DEALLOCATE(ARR)

#  define LIBPAW_POINTER_ALLOCATE(ARR,SIZE) ABI_ALLOCATE(ARR,SIZE)
#  define LIBPAW_POINTER_DEALLOCATE(ARR) ABI_DEALLOCATE(ARR)

#  define LIBPAW_DATATYPE_ALLOCATE(ARR,SIZE) ABI_DATATYPE_ALLOCATE(ARR,SIZE)
#  define LIBPAW_DATATYPE_DEALLOCATE(ARR) ABI_DATATYPE_DEALLOCATE(ARR)

#  define LIBPAW_BOUND1_ALLOCATE(ARR,BND1) ABI_ALLOCATE(ARR,(BND1))
#  define LIBPAW_BOUND2_ALLOCATE(ARR,BND1,BND2) ABI_ALLOCATE(ARR,(BND1,BND2))
#  define BOUNDS(LBND,UBND) LBND : UBND 


#  if defined HAVE_DFT_LIBXC
#    define LIBPAW_HAVE_LIBXC HAVE_DFT_LIBXC
#  else
#    undef LIBPAW_HAVE_LIBXC
#  endif


#  define LIBPAW_CONTIGUOUS ABI_CONTIGUOUS

/* =============================
 * ========= BIGDFT ============
 * ============================= */
#elif HAVE_LIBPAW_BIGDFT


#  define USE_DEFS use m_libpaw_defs


#  define USE_MPI_WRAPPERS use m_libpaw_mpi


#  define USE_MSG_HANDLING use m_libpaw_tools, only : wrtout => libpaw_wrtout, libpaw_msg_hndl
#  define MSG_COMMENT(msg) call libpaw_msg_hndl(msg,"COMMENT","PERS")
#  define MSG_WARNING(msg) call libpaw_msg_hndl(msg,"WARNING","PERS")
#  define MSG_ERROR(msg)   call libpaw_msg_hndl(msg,"ERROR"  ,"PERS")
#  define MSG_BUG(msg)     call libpaw_msg_hndl(msg,"BUG"    ,"PERS")


#  define HAVE_YAML


#  define USE_MEMORY_PROFILING use dynamic_memory

#  define LIBPAW_ALLOCATE(ARR,SIZE) ARR=f_malloc(to_array SIZE ) 
#  define LIBPAW_DEALLOCATE(ARR) call f_free(ARR)

#  define LIBPAW_POINTER_ALLOCATE(ARR,SIZE) ARR=f_malloc_ptr(to_array SIZE ) 
#  define LIBPAW_POINTER_DEALLOCATE(ARR) call f_free_ptr(ARR)

#  define LIBPAW_DATATYPE_ALLOCATE(ARR,SIZE) allocate(ARR SIZE)
#  define LIBPAW_DATATYPE_DEALLOCATE(ARR) deallocate(ARR)

#  define LIBPAW_BOUND1_ALLOCATE(ARR,BND1) ARR=f_malloc((/ BND1 /))
#  define LIBPAW_BOUND2_ALLOCATE(ARR,BND1,BND2) ARR=f_malloc((/ BND1 , BND2 /))
#  define BOUNDS(LBND,UBND) LBND .to. UBND 


#  define LIBPAW_HAVE_LIBXC


#  define LIBPAW_CONTIGUOUS 

/* =============================
 * ========= DEFAULT ===========
 * ============================= */
#else


#  define USE_DEFS use m_libpaw_defs


#  define USE_MPI_WRAPPERS use m_libpaw_mpi


#  define USE_MSG_HANDLING use m_libpaw_tools, only : wrtout => libpaw_wrtout, libpaw_msg_hndl
#  define MSG_COMMENT(msg) call libpaw_msg_hndl(msg,"COMMENT","PERS")
#  define MSG_WARNING(msg) call libpaw_msg_hndl(msg,"WARNING","PERS")
#  define MSG_ERROR(msg)   call libpaw_msg_hndl(msg,"ERROR"  ,"PERS")
#  define MSG_BUG(msg)     call libpaw_msg_hndl(msg,"BUG"    ,"PERS")
#  undef  HAVE_YAML


#  define USE_MEMORY_PROFILING

#  define LIBPAW_ALLOCATE(ARR,SIZE) allocate(ARR SIZE)
#  define LIBPAW_DEALLOCATE(ARR) deallocate(ARR)

#  define LIBPAW_POINTER_ALLOCATE(ARR,SIZE) allocate(ARR SIZE)
#  define LIBPAW_POINTER_DEALLOCATE(ARR) deallocate(ARR)

#  define LIBPAW_DATATYPE_ALLOCATE(ARR,SIZE) allocate(ARR SIZE)
#  define LIBPAW_DATATYPE_DEALLOCATE(ARR) deallocate(ARR)

#  define LIBPAW_BOUND1_ALLOCATE(ARR,BND1) allocate(ARR(BND1))
#  define LIBPAW_BOUND2_ALLOCATE(ARR,BND1,BND2) allocate(ARR(BND1,BND2))
#  define BOUNDS(LBND,UBND) LBND : UBND 


#  undef LIBPAW_HAVE_LIBXC


#  define LIBPAW_CONTIGUOUS 

/* =============================
 * =========== END =============
 * ============================= */
#endif
