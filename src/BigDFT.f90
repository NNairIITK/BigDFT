!> @file
!! BigDFT package performing ab initio calculation based on wavelets
!! @author
!!    Copyright (C) 2007-2011 BigDFT group
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
   logical :: exist_list
   !input variables
   type(atoms_data) :: atoms
   type(input_variables) :: inputs
   type(restart_objects) :: rst
   character(len=50), dimension(:), allocatable :: arr_posinp
   character(len=60) :: filename, radical
   ! atomic coordinates, forces, strten
   real(gp), dimension(6) :: strten
   real(gp), dimension(:,:), allocatable :: fxyz
   real(gp), dimension(:,:), pointer :: rxyz
   integer :: iconfig,nconfig,istat,group_size

   ! Start MPI in parallel version
   !in the case of MPIfake libraries the number of processors is automatically adjusted
   call bigdft_mpi_init(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

   ! GLOBAL TASK GROUPS: modify here to make multiple runs
   group_size=nproc
   call mpi_environment_set(bigdft_mpi,iproc,nproc,MPI_COMM_WORLD,group_size)
   
   call memocc_set_memory_limit(memorylimit)


   ! Read a possible radical format argument.
   call get_command_argument(1, value = radical, status = istat)
   if (istat > 0) then
      write(radical, "(A)") "input"   
   end if

   ! find out which input files will be used
   inquire(file="list_posinp",exist=exist_list)
   if (exist_list) then
      open(54,file="list_posinp")
      read(54,*) nconfig
      if (nconfig > 0) then 
         !allocation not referenced since memocc count not initialised
         allocate(arr_posinp(1:nconfig))

         do iconfig=1,nconfig
            read(54,*) arr_posinp(iconfig)
         enddo
      else
         !normal case
         nconfig=1
         allocate(arr_posinp(1:1))
         if (istat > 0) then
            arr_posinp(1)='posinp'// trim(bigdft_run_id_toa())
         else
            arr_posinp(1)=trim(radical)
         end if
      endif
      close(54)
   else
      nconfig=1
      allocate(arr_posinp(1:1))
         if (istat > 0) then
            arr_posinp='posinp' // trim(bigdft_run_id_toa())
         else
            arr_posinp(1)=trim(radical)
         end if
   end if

   do iconfig=1,nconfig

      ! Read all input files.
      !standard names
      call standard_inputfile_names(inputs, radical, bigdft_mpi%nproc)
      call read_input_variables(bigdft_mpi%iproc,trim(arr_posinp(iconfig)),inputs, atoms, rxyz)
      
      ! Decide whether we use the cubic or the linear version
      select case (inputs%inputpsiid)
      case (INPUT_PSI_EMPTY, INPUT_PSI_RANDOM, INPUT_PSI_CP2K, INPUT_PSI_LCAO, INPUT_PSI_MEMORY_WVL, &
            INPUT_PSI_DISK_WVL, INPUT_PSI_LCAO_GAUSS, INPUT_PSI_MEMORY_GAUSS, INPUT_PSI_DISK_GAUSS)
          rst%version = CUBIC_VERSION
      case (INPUT_PSI_LINEAR_AO, INPUT_PSI_MEMORY_LINEAR, INPUT_PSI_DISK_LINEAR)
          rst%version = LINEAR_VERSION
      end select

      !initialize memory counting
      !call memocc(0,iproc,'count','start')

      allocate(fxyz(3,atoms%nat+ndebug),stat=i_stat)
      call memocc(i_stat,fxyz,'fxyz',subname)

      call init_restart_objects(bigdft_mpi%iproc,inputs%matacc,atoms,rst,subname)

      if (bigdft_mpi%iproc == 0) then
         call print_general_parameters(bigdft_mpi%nproc,inputs,atoms)
         !call write_input_parameters(inputs,atoms)
      end if

      !if other steps are supposed to be done leave the last_run to minus one
      !otherwise put it to one
      if (inputs%last_run == -1 .and. inputs%ncount_cluster_x <=1 .or. inputs%ncount_cluster_x <= 1) then
         inputs%last_run = 1
      end if

      call call_bigdft(bigdft_mpi%nproc,bigdft_mpi%iproc,atoms,rxyz,inputs,etot,fxyz,strten,fnoise,rst,infocode)

      if (inputs%ncount_cluster_x > 1) then
         filename=trim(inputs%dir_output)//'geopt.mon'
         open(unit=16,file=filename,status='unknown',position='append')
         if (iproc ==0 ) write(16,*) '----------------------------------------------------------------------------'
         if (iproc ==0 ) write(*,"(1x,a,2i5)") 'Wavefunction Optimization Finished, exit signal=',infocode
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

      call deallocate_atoms(atoms,subname) 

      if (inputs%inputPsiId==INPUT_PSI_LINEAR_AO .or. inputs%inputPsiId==INPUT_PSI_MEMORY_LINEAR &
          .or. inputs%inputPsiId==INPUT_PSI_DISK_LINEAR) then
          call destroy_DFT_wavefunction(rst%tmb)
          call deallocate_local_zone_descriptors(rst%tmb%lzd, subname)
      end if


      if(inputs%linear /= INPUT_IG_OFF .and. inputs%linear /= INPUT_IG_LIG) &
           & call deallocateBasicArraysInput(inputs%lin)

      call free_restart_objects(rst,subname)

      i_all=-product(shape(rxyz))*kind(rxyz)
      deallocate(rxyz,stat=i_stat)
      call memocc(i_stat,i_all,'rxyz',subname)
      i_all=-product(shape(fxyz))*kind(fxyz)
      deallocate(fxyz,stat=i_stat)
      call memocc(i_stat,i_all,'fxyz',subname)

      call free_input_variables(inputs)

      !finalize memory counting
      call memocc(0,0,'count','stop')

   enddo !loop over iconfig

   deallocate(arr_posinp)

   call mpi_environment_free(bigdft_mpi)
   !wait all processes before finalisation
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)

   call MPI_FINALIZE(ierr)

END PROGRAM BigDFT

