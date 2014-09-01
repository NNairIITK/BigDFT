!> @file
!! BigDFT package performing ab initio calculation based on wavelets
!! @author
!!    Copyright (C) 2007-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Main program to calculate electronic structures
program BigDFT
   use module_base
   use bigdft_run!module_types
   use module_interfaces, only: write_atomic_file
   use yaml_strings, only: f_strcpy
   use yaml_output, only: yaml_map
   use yaml_parse
   !use internal_coordinates, only : get_neighbors

   implicit none     !< As a general policy, we will have "implicit none" by assuming the same

   character(len=*), parameter :: subname='BigDFT' !< Used by memocc routine (timing)
   integer :: ierr,infocode!,iproc,nproc
   integer :: ncount_bigdft
   !input variables
   type(run_objects) :: runObj
   !output variables
   type(DFT_global_output) :: outs
   character(len=60), dimension(:), allocatable :: arr_posinp,arr_radical
   character(len=60) :: filename, run_id, posinp_id
   !information for mpi_initalization
   integer, dimension(4) :: mpi_info
   integer :: iconfig,nconfig!,ngroups,igroup
   real(kind=8),dimension(:,:),allocatable :: fxyz
   integer :: iat
   logical :: file_exists
   integer,dimension(:),allocatable :: atoms_ref
   type(yaml_cl_parse) :: parser !< command line parser
   type(dictionary), pointer :: run

   call f_lib_initialize()


   !define command-line options
   parser=yaml_cl_parse_null()
   !between these lines, for another executable using BigDFT as a blackbox,
   !other command line options can be specified
   !then the bigdft options can be specified
   call bigdft_command_line_options(parser)
   !parse command line
   call yaml_cl_parse_cmd_line(parser)

   !-finds the number of taskgroup size
   !-initializes the mpi_environment for each group
   !-decides the radical name for each run
   call bigdft_init_new(parser,mpi_info)

!!$   call bigdft_init(mpi_info,nconfig,run_id,ierr)

!!$   !just for backward compatibility
!!$   iproc=mpi_info(1)
!!$   nproc=mpi_info(2)
!!$
!!$   igroup=mpi_info(3)
!!$   !number of groups
!!$   ngroups=mpi_info(4)
   
!!$   !allocate arrays of run ids
!!$   allocate(arr_radical(abs(nconfig)))
!!$   allocate(arr_posinp(abs(nconfig)))
!!$
!!$   !here we call a routine which reads a possible radical format argument.
!!$   call bigdft_get_run_ids(nconfig,trim(run_id),arr_radical,arr_posinp,ierr)
!!$
!!$   do iconfig=1,abs(nconfig)
!!$      if (modulo(iconfig-1,bigdft_mpi%ngroup)==bigdft_mpi%igroup) then
!!$         !print *,'iconfig,arr_radical(iconfig),arr_posinp(iconfig)',arr_radical(iconfig),arr_posinp(iconfig),iconfig,igroup
!!$         ! Read all input files.
!!$         call f_strcpy(src=arr_radical(iconfig),dest=run_id)
!!$         call f_strcpy(src=arr_posinp(iconfig) ,dest=posinp_id)

   !case with parser information
   !this key will contain the runs which are associated to the current BigDFT instance
   run => dict_iter(parser%args .get. 'BigDFT')
   do while(associated(run))
      run_id = run // 'name'
      posinp_id = run // 'posinp' ! the position id has to be reconsidered

         call run_objects_init_from_files(runObj, run_id, posinp_id)
         call init_global_output(outs, runObj%atoms%astruct%nat)

         call call_bigdft(runObj,outs,bigdft_mpi%nproc,bigdft_mpi%iproc,infocode)

         if (runObj%inputs%ncount_cluster_x > 1) then
            if (bigdft_mpi%iproc ==0 ) call yaml_map('Wavefunction Optimization Finished, exit signal',infocode)
            ! geometry optimization
            call geopt(runObj,outs,bigdft_mpi%nproc,bigdft_mpi%iproc,ncount_bigdft)
         end if

         !if there is a last run to be performed do it now before stopping
         if (runObj%inputs%last_run == -1) then
            runObj%inputs%last_run = 1
            call call_bigdft(runObj, outs, bigdft_mpi%nproc,bigdft_mpi%iproc,infocode)
         end if

         if (runObj%inputs%ncount_cluster_x > 1) then
            !filename='final_'//trim(posinp_id)
            call f_strcpy(src='final_'//trim(posinp_id),dest=filename)
            if (bigdft_mpi%iproc == 0) call write_atomic_file(filename,outs%energy,runObj%atoms%astruct%rxyz, &
                 & runObj%atoms%astruct%ixyz_int, runObj%atoms,'FINAL CONFIGURATION',forces=outs%fxyz)
         else
            !filename='forces_'//trim(arr_posinp(iconfig))
            call f_strcpy(src='forces_'//trim(posinp_id),dest=filename)
            if (bigdft_mpi%iproc == 0) call write_atomic_file(filename,outs%energy,runObj%atoms%astruct%rxyz, &
                 & runObj%atoms%astruct%ixyz_int, runObj%atoms,'Geometry + metaData forces',forces=outs%fxyz)
         end if

         ! Deallocations.
         call deallocate_global_output(outs)
         call run_objects_free(runObj)
         !temporary
         !call f_malloc_dump_status()
!!$      end if
      run => dict_next(run)
   end do !loop over iconfig

!!$   deallocate(arr_posinp,arr_radical)

   call bigdft_finalize(ierr)

   !free command line parser information
   call yaml_cl_parse_free(parser)

   call f_lib_finalize()

END PROGRAM BigDFT
