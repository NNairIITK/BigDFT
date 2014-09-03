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
   character(len=60) :: filename,posinp_id!, run_id
   integer :: iconfig,nconfig!,ngroups,igroup
   real(kind=8),dimension(:,:),allocatable :: fxyz
   integer :: iat
   logical :: file_exists
   integer,dimension(:),allocatable :: atoms_ref
   type(dictionary), pointer :: run,options

   call f_lib_initialize()

   !call test_dictionaries0()
   !call test_error_handling()
   !call test_timing()
   call bigdft_command_line_options(options)

   !-finds the number of taskgroup size
   !-initializes the mpi_environment for each group
   !-decides the radical name for each run
   call bigdft_init(options)

   !case with parser information
   !this key will contain the runs which are associated to the current BigDFT instance
   run => dict_iter(options .get. 'BigDFT')
   do while(associated(run))
      call run_objects_init(runObj,run)
      call init_global_output(outs,runObj%atoms%astruct%nat)

      posinp_id = run // 'posinp' 
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
         run => dict_next(run)
   end do !loop over iconfig

   call dict_free(options)
   call bigdft_finalize(ierr)


   call f_lib_finalize()

END PROGRAM BigDFT



subroutine test_dictionaries0()
  use yaml_output                                                     !contains the routines for the yaml output
  use dictionaries                                                    !contains the dictionary routines
  implicit none
  type(dictionary),pointer :: d1, d2, d3
 
  call dict_init(d1)                                                  !initialize the dictionary ''d1''
  call set(d1//'toto',1)                                              !add the key ''toto'' to it and assign it the value 1
  call set(d1//'titi',1.d0)                                           !add the key ''titi'' to it and assign it the value 1.d0
  call set(d1//'tutu',(/ '1', '2' /))                                 !add the key ''tutu'' to it and assign it the array [1,2]

  call dict_init(d2)                                                  !initialize the array dictionary ''d2''
  call set(d2//'a',0)                                                 !add the key ''a'' and assign it the value 0
  call set(d1//'List',list_new((/.item.d2, .item.'4', .item.'1.0'/))) !create a list from ''d2'' and the values 4 and 1.0

  call yaml_dict_dump(d1)                                             !output the content of ''d1'' in the yaml format

  d3 => d1//'New key'                                                 !points to ''d1["New key"]'' (and creates it at the same time)
  call set(d3//'Example',4)                                           !add the key ''Example'' to ''d3'' and assign it the value 4
  call yaml_dict_dump(d3)                                             !output the content of ''d2'' in the yaml format

  call yaml_map('List length',dict_len(d1//'List'))                   !print the length of the key ''List'' in dictionary ''d1''
  call yaml_map('Dictionary size',dict_size(d1))                      !print the size of the dictionary ''d1''
  call dict_free(d1)                                                  !destroy the dictionary ''d1''
  
end subroutine test_dictionaries0


subroutine test_error_handling()
  use yaml_output
  use dictionaries
  implicit none
  integer :: ERR_TEST
  external :: abort_test

  call yaml_comment('Error Handling Module Test',hfill='~')                !just a comment

  call f_err_define(err_name='ERR_TEST',&                                  !define the error
       err_msg='This is the error message for the error "ERR_TEST" and'//&
       ' it is written extensively on purpose to see whether the yaml'//&
       ' module can still handle it',&
       err_action='For this error, contact the routine developer',&
       err_id=ERR_TEST,callback=abort_test)

  call yaml_map("Raising the TEST error, errcode",ERR_TEST)                !print that the error will now be triggered
  if (f_err_raise(.true.,'Extra message added',err_id=ERR_TEST)) return    !rasie the error and return

end subroutine test_error_handling

subroutine abort_test()
  use yaml_output
  implicit none
  call yaml_comment('printing error informations',hfill='_')               !just a comment indicating that the error informations will now be written
  call f_dump_last_error()                                                 !print the error information
end subroutine abort_test



subroutine test_timing()
  use dynamic_memory
  use time_profiling
  use dictionaries
  use yaml_output
  type(dictionary), pointer :: dict_tmp

  call sub1()
  call f_malloc_dump_status(dict_summary=dict_tmp)
  call yaml_dict_dump(dict_tmp)

end subroutine test_timing

subroutine sub1()
  use dynamic_memory
  call f_routine('sub1')
  call sub2()
  call sub2()
  call sub3()
  call f_release_routine()
end subroutine sub1

subroutine sub2()
  use dynamic_memory
  call f_routine('sub2')
  call waste_time()
  call f_release_routine()
end subroutine sub2

subroutine sub3()
  use dynamic_memory
  call f_routine('sub3')
  call waste_time()
  call f_release_routine()
end subroutine sub3

subroutine waste_time()
  implicit none
  integer :: i
  real(kind=8) :: tt
  tt=0.d0
  do i=1,1000000
      tt = tt + sin(real(i,kind=8))
  end do
end subroutine waste_time

