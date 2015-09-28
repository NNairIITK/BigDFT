!> @file
!! BigDFT package performing ab initio calculation based on wavelets
!! @author
!!    Copyright (C) 2007-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Experimental program to exchange forces after configuration for images
program NEB_images

   use BigDFT_API

   implicit none     !< As a general policy, we will have "implicit none" by assuming the same
   integer :: iproc,nproc,ierr,infocode
   integer :: ncount_bigdft
   !input variables
   type(run_objects) :: runObj
   !output variables
   type(DFT_global_output) :: outs
   character(len=60), dimension(:), allocatable :: arr_posinp,arr_radical
   character(len=60) :: run_id
   !information for mpi_initalization
   integer, dimension(4) :: mpi_info
   integer :: iconfig,nconfig,istat,group_size,ngroups,igroup

   !-finds the number of taskgroup size
   !-initializes the mpi_environment for each group
   !-decides the radical name for each run
   call bigdft_init(mpi_info,nconfig,run_id,ierr)

   !just for backward compatibility
   iproc=mpi_info(1)
   nproc=mpi_info(2)
   igroup=mpi_info(3)
   !number of groups simultaneously running
   ngroups=mpi_info(4)
  
   !allocate arrays of run ids
   allocate(arr_radical(abs(nconfig)))
   allocate(arr_posinp(abs(nconfig)))

   !Here we call  a routine which reads a possible radical format argument.
   call bigdft_get_run_ids(nconfig,trim(run_id),arr_radical,arr_posinp,ierr)

   do iconfig=1,abs(nconfig)
      if (modulo(iconfig-1,ngroups)==igroup) then
         !print *,'iconfig,arr_radical(iconfig),arr_posinp(iconfig)',arr_radical(iconfig),arr_posinp(iconfig),iconfig,igroup
         ! Read all input files. This should be the sole routine which is called to initialize the run.
         call run_objects_init_from_files(runObj, arr_radical(iconfig),arr_posinp(iconfig))
         call init_global_output(outs, runObj%atoms%nat)
         call call_bigdft(runObj, outs, bigdft_mpi%nproc,bigdft_mpi%iproc,infocode)

         call call_bigdft(runObj, outs, bigdft_mpi%nproc,bigdft_mpi%iproc,infocode)

         if (runObj%inputs%ncount_cluster_x > 1) then
            filename=trim('final_'//trim(arr_posinp(iconfig)))
            if (bigdft_mpi%iproc == 0) call write_atomic_file(filename,outs%energy,&
                 & runObj%rxyz,runObj%ixyz_int,runObj%atoms,'FINAL CONFIGURATION',forces=outs%fxyz)
         else
            filename=trim('forces_'//trim(arr_posinp(iconfig)))
            if (bigdft_mpi%iproc == 0) call write_atomic_file(filename,outs%energy,&
                 & runObj%rxyz,runObj%ixyz_int,runObj%atoms,'Geometry + metaData forces',forces=outs%fxyz)
         end if

         ! Deallocations.
         call deallocate_global_output(outs)
         call run_objects_free(runObj)

      end if
   enddo !loop over iconfig

   deallocate(arr_posinp,arr_radical)

   call bigdft_finalize(ierr)

 END PROGRAM NEB_images
