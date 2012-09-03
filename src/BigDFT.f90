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
   integer :: iproc,nproc,iat,j,i_stat,i_all,ierr,infocode
   integer :: ncount_bigdft
   real(gp) :: etot,sumx,sumy,sumz,fnoise
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
   integer :: iconfig,nconfig,istat

   ! Start MPI in parallel version
   !in the case of MPIfake libraries the number of processors is automatically adjusted
   call bigdft_mpi_init(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

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
            arr_posinp(1)='posinp'
         else
            arr_posinp(1)=trim(radical)
         end if
      endif
      close(54)
   else
      nconfig=1
      allocate(arr_posinp(1:1))
         if (istat > 0) then
            arr_posinp(1)='posinp'
         else
            arr_posinp(1)=trim(radical)
         end if
   end if

   do iconfig=1,nconfig

      ! Read all input files.
      !standard names
      call standard_inputfile_names(inputs, radical, nproc)
      call read_input_variables(iproc,trim(arr_posinp(iconfig)),inputs, atoms, rxyz)
      if (iproc == 0) then
         call print_general_parameters(nproc,inputs,atoms)
         !call write_input_parameters(inputs,atoms)
      end if

      !initialize memory counting
      !call memocc(0,iproc,'count','start')

      allocate(fxyz(3,atoms%nat+ndebug),stat=i_stat)
      call memocc(i_stat,fxyz,'fxyz',subname)

      call init_restart_objects(iproc,inputs%iacceleration,atoms,rst,subname)

      !if other steps are supposed to be done leave the last_run to minus one
      !otherwise put it to one
      if (inputs%last_run == -1 .and. inputs%ncount_cluster_x <=1 .or. inputs%ncount_cluster_x <= 1) then
         inputs%last_run = 1
      end if

      call call_bigdft(nproc,iproc,atoms,rxyz,inputs,etot,fxyz,strten,fnoise,rst,infocode)

      if (inputs%ncount_cluster_x > 1) then
         open(unit=16,file='geopt.mon',status='unknown',position='append')
         if (iproc ==0 ) write(16,*) '----------------------------------------------------------------------------'
         if (iproc ==0 ) write(*,"(1x,a,2i5)") 'Wavefunction Optimization Finished, exit signal=',infocode
         ! geometry optimization
         call geopt(nproc,iproc,rxyz,atoms,fxyz,strten,etot,rst,inputs,ncount_bigdft)
         close(16)
      end if

      !if there is a last run to be performed do it now before stopping
      if (inputs%last_run == -1) then
         inputs%last_run = 1
         call call_bigdft(nproc,iproc,atoms,rxyz,inputs,etot,fxyz,strten,fnoise,rst,infocode)
      end if

      if (inputs%ncount_cluster_x > 1) then
         filename=trim('final_'//trim(arr_posinp(iconfig)))
         if (iproc == 0) call write_atomic_file(filename,etot,rxyz,atoms,'FINAL CONFIGURATION',forces=fxyz)
      else
         filename=trim('forces_'//trim(arr_posinp(iconfig)))
         if (iproc == 0) call write_atomic_file(filename,etot,rxyz,atoms,'Geometry + metaData forces',forces=fxyz)
      end if

      if (iproc == 0) then
         sumx=0.d0
         sumy=0.d0
         sumz=0.d0
         write(*,'(1x,a,19x,a)') 'Final values of the Forces for each atom'
         do iat=1,atoms%nat
            write(*,'(1x,i5,1x,a6,3(1x,1pe12.5))') &
            iat,trim(atoms%atomnames(atoms%iatype(iat))),(fxyz(j,iat),j=1,3)
            sumx=sumx+fxyz(1,iat)
            sumy=sumy+fxyz(2,iat)
            sumz=sumz+fxyz(3,iat)
         enddo
         !$$        if (.not. inputs%gaussian_help .or. .true.) then !zero of the forces calculated
         !$$           write(*,'(1x,a)')'the sum of the forces is'
         !$$           write(*,'(1x,a16,3x,1pe16.8)')'x direction',sumx
         !$$           write(*,'(1x,a16,3x,1pe16.8)')'y direction',sumy
         !$$           write(*,'(1x,a16,3x,1pe16.8)')'z direction',sumz
         !$$        end if
      endif

      call deallocate_atoms(atoms,subname) 

      if (inputs%inputPsiId==INPUT_PSI_LINEAR_AO .or. inputs%inputPsiId==INPUT_PSI_MEMORY_LINEAR &
          .or. inputs%inputPsiId==INPUT_PSI_DISK_LINEAR) then
          call destroy_DFT_wavefunction(rst%tmb)
          call deallocate_local_zone_descriptors(rst%tmb%lzd, subname)
      end if


!      call deallocate_lr(rst%Lzd%Glr,subname)    
!      call deallocate_local_zone_descriptors(rst%Lzd, subname)
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

   !wait all processes before finalisation
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)

   call MPI_FINALIZE(ierr)

END PROGRAM BigDFT

