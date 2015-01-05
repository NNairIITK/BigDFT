/* libpaw.h */

/*
 * This part is part of the libPAW library.
 * It has to be customized according to the host code.
 * For the time being there are 2 known host codes:
 * ABINIT (www.abinit.org) and BigDFT (bigdft.org).
 */

/*
 * Copyright (C) 2014-2014 ABINIT Group (MT)
 *
 * This file is part of the ABINIT software package. For license information,
 * please see the COPYING file in the top-level directory of the ABINIT source
 * distribution.
 *
 */

/* config.h should contain all preprocessing directives */
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

/* =============================
 * ========= ABINIT ============
 * ============================= */
#if defined HAVE_LIBPAW_ABINIT

/* ABINIT specific macros */
#include "abi_common.h"

/* Constants and defs */
#define USE_DEFS use defs_basis

/* MPI wrappers */
#define USE_MPI_WRAPPERS use m_xmpi

/* Messages, errors */
#define USE_MSG_HANDLING use m_errors, only : msg_hndl
/* Other macros already defined in abi_common.h */

/* Allocation/deallocation with memory profiling */
#define USE_MEMORY_PROFILING use m_profiling_abi
#define LIBPAW_ALLOCATE(ARR,SIZE) ABI_ALLOCATE(ARR,SIZE)
#define LIBPAW_DEALLOCATE(ARR)    ABI_DEALLOCATE(ARR)
#define LIBPAW_POINTER_ALLOCATE(ARR,SIZE) ABI_ALLOCATE(ARR,SIZE)
#define LIBPAW_POINTER_DEALLOCATE(ARR)    ABI_DEALLOCATE(ARR)
#define LIBPAW_DATATYPE_ALLOCATE(ARR,SIZE) ABI_DATATYPE_ALLOCATE(ARR,SIZE)
#define LIBPAW_DATATYPE_DEALLOCATE(ARR)    ABI_DATATYPE_DEALLOCATE(ARR)


/* =============================
 * ========= BIGDFT ============
 * ============================= */
#elif HAVE_LIBPAW_BIGDFT

/* Constants and defs */
#define USE_DEFS use m_libpaw_defs

/* MPI wrappers */
#define USE_MPI_WRAPPERS use m_libpaw_mpi

/* Messages, errors */
#define USE_MSG_HANDLING use m_libpaw_tools, only : wrtout => libpaw_wrtout, libpaw_msg_hndl
#define MSG_COMMENT(msg) call libpaw_msg_hndl(msg,"COMMENT","PERS",__FILE__,__LINE__)
#define MSG_WARNING(msg) call libpaw_msg_hndl(msg,"WARNING","PERS",__FILE__,__LINE__)
#define MSG_ERROR(msg)   call libpaw_msg_hndl(msg,"ERROR"  ,"PERS",__FILE__,__LINE__)
#define MSG_BUG(msg)     call libpaw_msg_hndl(msg,"BUG"    ,"PERS",__FILE__,__LINE__)

/* Allocation/deallocation with memory profiling */
#define USE_MEMORY_PROFILING use dynamic_memory
#define LIBPAW_ALLOCATE(ARR,SIZE) ARR=f_malloc(to_array SIZE ) 
#define LIBPAW_DEALLOCATE(ARR)    call f_free(ARR)
#define LIBPAW_POINTER_ALLOCATE(ARR,SIZE) ARR=f_malloc_ptr(to_array SIZE ) 
#define LIBPAW_POINTER_DEALLOCATE(ARR)    call f_free_ptr(ARR)
#define LIBPAW_DATATYPE_ALLOCATE(ARR,SIZE) allocate(ARR SIZE)
#define LIBPAW_DATATYPE_DEALLOCATE(ARR)    deallocate(ARR)


/* =============================
 * ========= DEFAULT ===========
 * ============================= */
#else

/* Constants and defs */
#define USE_DEFS use m_libpaw_defs

/* MPI wrappers */
#define USE_MPI_WRAPPERS use m_libpaw_mpi

/* Messages, errors */
#define USE_MSG_HANDLING use m_libpaw_tools, only : wrtout => libpaw_wrtout, libpaw_msg_hndl
#define MSG_COMMENT(msg) call libpaw_msg_hndl(msg,"COMMENT","PERS")
#define MSG_WARNING(msg) call libpaw_msg_hndl(msg,"WARNING","PERS")
#define MSG_ERROR(msg)   call libpaw_msg_hndl(msg,"ERROR"  ,"PERS")
#define MSG_BUG(msg)     call libpaw_msg_hndl(msg,"BUG"    ,"PERS")

/* Allocation/deallocation */
#define USE_MEMORY_PROFILING
#define LIBPAW_ALLOCATE(ARR,SIZE) allocate(ARR SIZE)
#define LIBPAW_DEALLOCATE(ARR)    deallocate(ARR)
#define LIBPAW_POINTER_ALLOCATE(ARR,SIZE) allocate(ARR SIZE)
#define LIBPAW_POINTER_DEALLOCATE(ARR)    deallocate(ARR)
#define LIBPAW_DATATYPE_ALLOCATE(ARR,SIZE) allocate(ARR SIZE)
#define LIBPAW_DATATYPE_DEALLOCATE(ARR)    deallocate(ARR)


/* =============================
 * =========== END =============
 * ============================= */
#endif
