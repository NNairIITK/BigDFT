!> @file
!! BigDFT package performing ab initio calculation based on wavelets
!! @author
!!    Copyright (C) 2007-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Experimental program to exchange forces after configuration for images
program NEB_images

   use BigDFT_API

   implicit none     !< As a general policy, we will have "implicit none" by assuming the same
   integer :: iproc,nproc,i_stat,i_all,ierr,infocode
   integer :: ncount_bigdft
   real(gp) :: etot,fnoise
   !input variables
   type(run_objects) :: runObj
   character(len=60), dimension(:), allocatable :: arr_posinp,arr_radical
   character(len=60) :: run_id
   ! atomic coordinates, forces, strten
   !information for mpi_initalization
   integer, dimension(4) :: mpi_info
   real(gp), dimension(6) :: strten
   real(gp), dimension(:,:), allocatable :: fxyz
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

   !here we call  a routine which
   ! Read a possible radical format argument.
   call bigdft_get_run_ids(nconfig,trim(run_id),arr_radical,arr_posinp,ierr)

   do iconfig=1,abs(nconfig)
      if (modulo(iconfig-1,ngroups)==igroup) then
         !print *,'iconfig,arr_radical(iconfig),arr_posinp(iconfig)',arr_radical(iconfig),arr_posinp(iconfig),iconfig,igroup
         ! Read all input files. This should be the sole routine which is called to initialize the run.
         call run_objects_set_from_files(runObj, arr_radical(iconfig),arr_posinp(iconfig))

         !here we should define a routine to extract the number of atoms and the positions, and allocate forces array

         allocate(fxyz(3,runObj%atoms%nat+ndebug),stat=i_stat)
         call memocc(i_stat,fxyz,'fxyz',subname)

         call call_bigdft(runObj,bigdft_mpi%nproc,bigdft_mpi%iproc,etot,fxyz,strten,fnoise,infocode)


         call call_bigdft(runObj,bigdft_mpi%nproc,bigdft_mpi%iproc,etot,fxyz,strten,fnoise,infocode)

         if (runObj%inputs%ncount_cluster_x > 1) then
            filename=trim('final_'//trim(arr_posinp(iconfig)))
            if (bigdft_mpi%iproc == 0) call write_atomic_file(filename,etot,&
                 & runObj%rxyz,runObj%atoms,'FINAL CONFIGURATION',forces=fxyz)
         else
            filename=trim('forces_'//trim(arr_posinp(iconfig)))
            if (bigdft_mpi%iproc == 0) call write_atomic_file(filename,etot,&
                 & runObj%rxyz,runObj%atoms,'Geometry + metaData forces',forces=fxyz)
         end if

         i_all=-product(shape(fxyz))*kind(fxyz)
         deallocate(fxyz,stat=i_stat)
         call memocc(i_stat,i_all,'fxyz',subname)

         call run_objects_free(runObj)

      end if
   enddo !loop over iconfig

   deallocate(arr_posinp,arr_radical)

   call bigdft_finalize(ierr)

 END PROGRAM NEB_images
