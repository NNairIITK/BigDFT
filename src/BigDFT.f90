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
   use module_types
   use module_interfaces
   use yaml_output
   use internal_coordinates, only : get_neighbors

   implicit none     !< As a general policy, we will have "implicit none" by assuming the same

   character(len=*), parameter :: subname='BigDFT' !< Used by memocc routine (timing)
   integer :: iproc,nproc,ierr,infocode
   integer :: ncount_bigdft
   !input variables
   type(run_objects) :: runObj
   !output variables
   type(DFT_global_output) :: outs
   character(len=60), dimension(:), allocatable :: arr_posinp,arr_radical
   character(len=60) :: filename, run_id
   !information for mpi_initalization
   integer, dimension(4) :: mpi_info
   integer :: iconfig,nconfig,ngroups,igroup
   real(kind=8),dimension(:,:),allocatable :: fxyz
   integer :: iat
   logical :: file_exists
   integer,dimension(:),allocatable :: atoms_ref

   call f_lib_initialize()

!call test_yaml()
!call test_dictionaries0()
!call test_error_handling()
!call test_timing()

   !-finds the number of taskgroup size
   !-initializes the mpi_environment for each group
   !-decides the radical name for each run
   call bigdft_init(mpi_info,nconfig,run_id,ierr)

   !just for backward compatibility
   iproc=mpi_info(1)
   nproc=mpi_info(2)

   igroup=mpi_info(3)
   !number of groups
   ngroups=mpi_info(4)
   
   !allocate arrays of run ids
   allocate(arr_radical(abs(nconfig)))
   allocate(arr_posinp(abs(nconfig)))

   !here we call a routine which reads a possible radical format argument.
   call bigdft_get_run_ids(nconfig,trim(run_id),arr_radical,arr_posinp,ierr)

   do iconfig=1,abs(nconfig)
      if (modulo(iconfig-1,ngroups)==igroup) then
         !print *,'iconfig,arr_radical(iconfig),arr_posinp(iconfig)',arr_radical(iconfig),arr_posinp(iconfig),iconfig,igroup
         ! Read all input files.
         call run_objects_init_from_files(runObj, arr_radical(iconfig), arr_posinp(iconfig))
         call init_global_output(outs, runObj%atoms%astruct%nat)


         !!! THIS IS TEMPORARY, SHOULD BE DONE IN A BETTER WAY ####################
         !!inquire(file='posinp.fix',exist=file_exists)
         !!if (file_exists) then
         !!    atoms_ref = f_malloc(runObj%atoms%astruct%nat,id='atoms_ref')
         !!    open(unit=123,file='posinp.fix')
         !!    do iat=1,runObj%atoms%astruct%nat
         !!        read(123,*) atoms_ref(iat), runObj%atoms%astruct%fix_int(1:3,iat)
         !!    end do
         !!    close(unit=123)
         !!    if (iproc==0) call yaml_map('before: runObj%atoms%astruct%ixyz_int',runObj%atoms%astruct%ixyz_int)
         !!    !!call get_neighbors(runObj%atoms%astruct%rxyz, runObj%atoms%astruct%nat, runObj%atoms%astruct%ixyz_int(1,:), &
         !!    !!     runObj%atoms%astruct%ixyz_int(2,:), runObj%atoms%astruct%ixyz_int(3,:),atoms_ref)
         !!    call f_free(atoms_ref)
         !!    if (iproc==0) call yaml_map('after: runObj%atoms%astruct%ixyz_int',runObj%atoms%astruct%ixyz_int)
         !!else
         !!    call get_neighbors(runObj%atoms%astruct%rxyz, runObj%atoms%astruct%nat, runObj%atoms%astruct%ixyz_int(1,:), &
         !!         runObj%atoms%astruct%ixyz_int(2,:), runObj%atoms%astruct%ixyz_int(3,:))
         !!end if
         !!! ######################################################################

         call call_bigdft(runObj,outs,bigdft_mpi%nproc,bigdft_mpi%iproc,infocode)

         if (runObj%inputs%ncount_cluster_x > 1) then
            if (iproc ==0 ) call yaml_map('Wavefunction Optimization Finished, exit signal',infocode)
            ! geometry optimization
            call geopt(runObj,outs,bigdft_mpi%nproc,bigdft_mpi%iproc,ncount_bigdft)
         end if

         !if there is a last run to be performed do it now before stopping
         if (runObj%inputs%last_run == -1) then
            runObj%inputs%last_run = 1
            call call_bigdft(runObj, outs, bigdft_mpi%nproc,bigdft_mpi%iproc,infocode)
         end if

         if (runObj%inputs%ncount_cluster_x > 1) then
            filename=trim('final_'//trim(arr_posinp(iconfig)))
            if (bigdft_mpi%iproc == 0) call write_atomic_file(filename,outs%energy,runObj%atoms%astruct%rxyz, &
                 & runObj%atoms%astruct%ixyz_int, runObj%atoms,'FINAL CONFIGURATION',forces=outs%fxyz)
         else
            filename=trim('forces_'//trim(arr_posinp(iconfig)))
            if (bigdft_mpi%iproc == 0) call write_atomic_file(filename,outs%energy,runObj%atoms%astruct%rxyz, &
                 & runObj%atoms%astruct%ixyz_int, runObj%atoms,'Geometry + metaData forces',forces=outs%fxyz)
         end if

         ! Deallocations.
         call deallocate_global_output(outs)
         call run_objects_free(runObj, subname)
         !temporary
         !call f_malloc_dump_status()
      end if
   enddo !loop over iconfig

   deallocate(arr_posinp,arr_radical)

   call bigdft_finalize(ierr)

   call f_lib_finalize()

END PROGRAM BigDFT

   subroutine test_yaml()
     use yaml_output

     call yaml_comment('Yaml Invoice Example',hfill='-')
     call yaml_map('invoice',34843)
     call yaml_map('date',trim(yaml_date_toa()))
     call yaml_mapping_open('bill-to',label='id001')
       call yaml_map('given','Chris')
       call yaml_mapping_open('address')
         call yaml_mapping_open('lines')
           call yaml_scalar('458 Walkman Dr.')
           call yaml_scalar('Suite #292')
         call yaml_mapping_close()
       call yaml_mapping_close()
     call yaml_mapping_close()
     call yaml_map('ship_to','*id001')
     call yaml_open_sequence('product')
       call yaml_sequence(advance='no')
       call yaml_map('sku','BL394D')
       call yaml_map('quantity',4)
       call yaml_map('description','Basketball')
       call yaml_map('price',450.,fmt='(f6.2)')
       call yaml_map('parcel dimensions',(/30,32,35/))
       call yaml_sequence(advance='no')
       call yaml_mapping_open(flow=.true.)
         call yaml_map('sku','BL4438H')
         call yaml_map('quantity',1)
         call yaml_newline()
         call yaml_map('description','Super Hoop')
         call yaml_map('price',2392.,fmt='(f8.2)')
         call yaml_map('parcel dimensions',(/120,20,15/))
       call yaml_mapping_close()
     call yaml_sequence_close()
     call yaml_map('tax',251.42,fmt='(f6.2)')
     call yaml_map('total',4443.52d0,fmt='(f6.2)') !wrong format on purpose
     call yaml_map('comments','Late afternoon is best. Backup contact is Nancy Billsmer @ 338-4338.')

   end subroutine test_yaml


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

  d3 => d1//'New key'                                                 !point to ''d1''?
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

