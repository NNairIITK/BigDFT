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

   implicit none     !< As a general policy, we will have "implicit none" by assuming the same

   character(len=*), parameter :: subname='BigDFT' !< Used by memocc routine (timing)
   integer :: iproc,nproc,i_stat,i_all,ierr,infocode
   integer :: ncount_bigdft
   real(gp) :: etot,fnoise
   !input variables
   type(atoms_data) :: atoms
   type(input_variables) :: inputs
   type(restart_objects) :: rst
   character(len=60), dimension(:), allocatable :: arr_posinp,arr_radical
   character(len=60) :: filename, run_id
   ! atomic coordinates, forces, strten
   !information for mpi_initalization
   integer, dimension(4) :: mpi_info
   real(gp), dimension(6) :: strten
   real(gp), dimension(:,:), allocatable :: fxyz
   real(gp), dimension(:,:), pointer :: rxyz
   integer :: iconfig,nconfig,ngroups,igroup

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
   
   
   !Nullify LZD for cubic version (new input guess)
   call nullify_local_zone_descriptors(rst%tmb%lzd)


  
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

         call bigdft_set_input(arr_radical(iconfig),arr_posinp(iconfig),rxyz,inputs,atoms)

         !here we should define a routine to extract the number of atoms and the positions, and allocate forces array

         allocate(fxyz(3,atoms%astruct%nat+ndebug),stat=i_stat)
         call memocc(i_stat,fxyz,'fxyz',subname)
         call init_restart_objects(bigdft_mpi%iproc,inputs,atoms,rst,subname)

         call call_bigdft(bigdft_mpi%nproc,bigdft_mpi%iproc,atoms,rxyz,inputs,etot,fxyz,strten,fnoise,rst,infocode)

         if (inputs%ncount_cluster_x > 1) then
            filename=trim(inputs%dir_output)//'geopt.mon'
            open(unit=16,file=filename,status='unknown',position='append')
            if (iproc ==0 ) write(16,*) '----------------------------------------------------------------------------'
            if (iproc ==0 ) call yaml_map('Wavefunction Optimization Finished, exit signal',infocode)
            !if (iproc ==0 ) write(*,"(1x,a,2i5)") 'Wavefunction Optimization Finished, exit signal=',infocode
            ! geometry optimization
            call geopt(bigdft_mpi%nproc,bigdft_mpi%iproc,rxyz,atoms,fxyz,strten,etot,rst,inputs,ncount_bigdft)
            close(16)
         end if

         !if there is a last run to be performed do it now before stopping
         if (inputs%last_run == -1) then
            inputs%last_run = 1
            call call_bigdft(bigdft_mpi%nproc,bigdft_mpi%iproc,atoms,rxyz,inputs,etot,fxyz,strten,fnoise,rst,infocode)
         end if

         if (inputs%ncount_cluster_x > 1) then
            filename=trim('final_'//trim(arr_posinp(iconfig)))
            if (bigdft_mpi%iproc == 0) call write_atomic_file(filename,etot,rxyz,atoms,'FINAL CONFIGURATION',forces=fxyz)
         else
            filename=trim('forces_'//trim(arr_posinp(iconfig)))
            if (bigdft_mpi%iproc == 0) call write_atomic_file(filename,etot,rxyz,atoms,'Geometry + metaData forces',forces=fxyz)
         end if


     !always deallocate lzd for new input guess
     !call deallocate_lzd(rst%tmb%lzd, subname)
     ! Modified by SM
     call deallocate_local_zone_descriptors(rst%tmb%lzd, subname)

         i_all=-product(shape(rxyz))*kind(rxyz)
         deallocate(rxyz,stat=i_stat)
         call memocc(i_stat,i_all,'rxyz',subname)
         i_all=-product(shape(fxyz))*kind(fxyz)
         deallocate(fxyz,stat=i_stat)
         call memocc(i_stat,i_all,'fxyz',subname)

         call free_restart_objects(rst,subname)

         call deallocate_atoms(atoms,subname) 
         
         !temporary
         !call f_malloc_dump_status()
         call bigdft_free_input(inputs)

      end if
   enddo !loop over iconfig

   deallocate(arr_posinp,arr_radical)

   call bigdft_finalize(ierr)

END PROGRAM BigDFT

