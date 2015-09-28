!> @file
!!    Routines to bind fake argumentf for cuda solver
!! @author
!!    Copyright (C) 2002-2015 BigDFT group  (LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
 subroutine cuda_estimate_memory_needs_cu(iproc,n,&
      geo,plansSize,maxPlanSize,freeGPUSize, totalGPUSize )
   use dictionaries
   call f_err_throw('We should not enter into the cuda_estimation_routine')
 end subroutine cuda_estimate_memory_needs_cu

 subroutine pad_data()
   use dictionaries
   call f_err_throw('We should not enter into the pad_data routine')
 end subroutine pad_data

 subroutine unpad_data()
   use dictionaries
   call f_err_throw('We should not enter into the unpad_data routine')
 end subroutine unpad_data

 subroutine finalize_reduction_kernel()
   use dictionaries
   call f_err_throw('We should not enter into the finalize_reduction_kernel routine')
 end subroutine finalize_reduction_kernel
