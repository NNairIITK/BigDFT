!> @file 
!!   Routines to test atomic forces
!! @author
!!   Copyright (C) 2005-2011 BigDFT group 
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 


!> Runs BigDFT and test whether the forces are the 
!! derivative of the energy.
!! Performs the integration of the calculated forces over
!! some random displacement and compare the result with the 
!! difference of the energy between the final and the initial 
!! position
!! @warning
!!    Date: 10/07; THIS PROGRAM MUST BE COMPLETELY CHANGED
!!    Date: Feb 2011:  This program was modified and updated by Ali Sadeghi
program test_forces

   use module_base
   use module_types
   use module_interfaces
   use m_ab6_symmetry
   use yaml_output

   implicit none
   character(len=*), parameter :: subname='test_forces'
   integer :: iproc,nproc,iat,i_stat,i_all,ierr,infocode,istat
   real(gp) :: etot,fnoise
   logical :: exist_list
   !input variables
   type(atoms_data) :: atoms
   type(input_variables) :: inputs
   type(restart_objects) :: rst
   character(len=50), dimension(:), allocatable :: arr_posinp
   character(len=60), parameter :: filename="list_posinp"
   ! atomic coordinates, forces
   real(gp), dimension(:,:), allocatable :: fxyz
   real(gp), dimension(:,:), pointer :: rxyz,drxyz
   integer :: iconfig,nconfig
   integer :: ipath,npath
   real(gp):: dx,etot0,path,fdr
   !parameter (dx=1.d-2 , npath=2*16+1)  ! npath = 2*n+1 where n=2,4,6,8,...
   parameter (dx=1.d-2 , npath=5)
   real(gp) :: simpson(1:npath)
   character(len=60) :: radical
   real(gp), dimension(6) :: strten

   ! Start MPI in parallel version
   !in the case of MPIfake libraries the number of processors is automatically adjusted
   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
   call mpi_environment_set(bigdft_mpi,iproc,nproc,MPI_COMM_WORLD,0)

   call memocc_set_memory_limit(memorylimit)

   ! Read a possible radical format argument.
   call get_command_argument(1, value = radical, status = istat)
   if (istat > 0) then
      write(radical, "(A)") "input"
   end if


   ! find out which input files will be used
   inquire(file=filename,exist=exist_list)
   if (exist_list) then
      open(54,file=filename)
      read(54,*) nconfig
      if (nconfig > 0) then 
         !allocation not referenced since memocc count not initialised
         allocate(arr_posinp(1:nconfig))

         do iconfig=1,nconfig
            read(54,*) arr_posinp(iconfig)
         enddo
      else
         nconfig=1
         allocate(arr_posinp(1:1))
         arr_posinp(1)='posinp'
      endif
      close(54)
   else
      nconfig=1
      allocate(arr_posinp(1:1))
      arr_posinp(1)='posinp'
   end if

   !prepare the array of the correct Simpson's rule weigths for the integration
   if (mod(npath,2).ne.1) stop 'the number of iteration steps has to be odd'
   simpson(1)=1.d0/3.d0
   simpson(2)=4.d0/3.d0
   do ipath=3,npath-2,2
      simpson(ipath)=2.d0/3.d0
      simpson(ipath+1)=4.d0/3.d0
   enddo
   simpson(npath)=1.d0/3.d0



   do iconfig=1,nconfig

      !welcome screen
      if (iproc==0) call print_logo()
      if (iproc==0) then 
         print*
         print*
         print*,'*******************************************************************************************************'
         print*,"This is a test program to verify  whether the force are the derivative of the energy: F=-dE/dR"
         print '(a,i3,a,f7.4,a,f6.4,a)', "It performs the integration of the calculated forces over " , npath,  &
            &   " random displacement (in the range [", -dx, ",",dx,"] a.u.)"
         print*, "and compares the result with the difference of the energy  between the final and the initial position:"// &
            &   " E2-E1 = -Integral F.dR" 
         print*," The advantage is two fold: 1) avoiding cancellation error in finite difference derivative," // & 
         " 2) considernig the forces over all atoms "
         print*,'*********************************************************************************************************'
         print*
         print*
      endif

      ! Read all input files.
      !standard names
      call standard_inputfile_names(inputs,radical,nproc)
      call read_input_variables(iproc,trim(arr_posinp(iconfig)),inputs, atoms, rxyz)
      !     if (iproc == 0) then
      !       call print_general_parameters(nproc,inputs,atoms)
      !    end if


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

      call init_restart_objects(iproc,inputs%matacc,atoms,rst,subname)

      !if other steps are supposed to be done leave the last_run to minus one
      !otherwise put it to one
      if (inputs%last_run == -1 .and. inputs%ncount_cluster_x <=1 .or. inputs%ncount_cluster_x <= 1) then
         inputs%last_run = 1
      end if

      ! path integral   
      path=0.d0
      !calculate the displacement at each integration step
      !(use sin instead of random numbers)
      allocate(drxyz(1:3,1:atoms%nat))
      do iat=1,atoms%nat
         drxyz(1,iat)=dx*sin(iat+.2d0)   
         drxyz(2,iat)=dx*sin(iat+.4d0)  
         drxyz(3,iat)=dx*sin(iat+.7d0)  
      end do

      ! loop for ipath 
      do ipath=1,npath

         !update atomic positions alog the path
         if(ipath>1) then
            rxyz(:,:)=rxyz(:,:)+drxyz(:,:)
            inputs%inputPsiId=1
         end if

         if (iproc == 0) then
            call print_general_parameters(nproc,inputs,atoms) ! to know the new positions
         end if

         call call_bigdft(nproc,iproc,atoms,rxyz,inputs,etot,fxyz,strten,fnoise,rst,infocode)
         !        inputs%inputPsiId=0   ! change PsiId to 0 if you want to  generate a new Psi and not use the found one

         if (iproc == 0 ) write(*,"(1x,a,2i5)") 'Wavefunction Optimization Finished, exit signal=',infocode

         if (ipath == 1 ) etot0=etot
         !   do one step of the path integration
         if (iproc == 0) then
            !integrate forces*displacement
            !fdr=sum(fxyz(1:3,1:atoms%nat)*drxyz(1:3,1:atoms%nat))
            fdr=sum(fxyz(:,:)*drxyz(:,:))
            path=path-simpson(ipath)*fdr
            write(*,"('path iter:',i3,'   -F.dr=',e13.5,'    path integral=',e13.5 )") ipath,-fdr, path 

            !Print atomic forces
            call write_forces(atoms,fxyz)
         end if
      end do !loop over ipath


      deallocate(drxyz)

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
      if (iproc==0) then 
         write(*,*) 
         write(*,*) 'Check correctness of forces'
         write(*,*) 'Difference of total energies =',etot-etot0
         write(*,*) 'Integral force*displacement = ',path
         write(*,*) 'Difference = ',(etot-etot0)-path
         write(*,*) 
      endif

   enddo !loop over iconfig

   deallocate(arr_posinp)

   call MPI_FINALIZE(ierr)

END PROGRAM test_forces
