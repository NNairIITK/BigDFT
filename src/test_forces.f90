!> @file 
!!   Routines to test atomic forces
!! @author
!!   Copyright (C) 2005-2013 BigDFT group 
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
!!    Date: 10/07: THIS PROGRAM MUST BE COMPLETELY CHANGED
!!    Date: Feb 2011:  This program was modified and updated by Ali Sadeghi
program test_forces

   use module_base
   use module_types
   use module_interfaces
   use m_ab6_symmetry
   use yaml_output

   implicit none
   character(len=*), parameter :: subname='test_forces'
   integer :: iproc,nproc,iat,ierr,infocode!,istat
   !logical :: exist_list
   !input variables
   type(run_objects) :: runObj
   type(DFT_global_output) :: outs
   character(len=60), parameter :: filename="list_posinp"
   character(len=60), dimension(:), allocatable :: arr_posinp,arr_radical
   character(len=60) :: run_id
   ! atomic coordinates, forces
   real(gp), dimension(:,:), pointer :: drxyz
   integer :: iconfig,nconfig,igroup,ngroups
   integer :: ipath,npath
   real(gp):: dx,etot0,path,fdr
   !parameter (dx=1.d-2 , npath=2*16+1)  ! npath = 2*n+1 where n=2,4,6,8,...
   parameter (dx=1.d-2 , npath=5)
   real(gp) :: simpson(1:npath)
   !character(len=60) :: radical
   integer, dimension(4) :: mpi_info

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

   !here we call  a routine which
   ! Read a possible radical format argument.
   call bigdft_get_run_ids(nconfig,trim(run_id),arr_radical,arr_posinp,ierr)

   !prepare the array of the correct Simpson's rule weigths for the integration
   if (mod(npath,2).ne.1) stop 'the number of iteration steps has to be odd'
   simpson(1)=1.d0/3.d0
   simpson(2)=4.d0/3.d0
   do ipath=3,npath-2,2
      simpson(ipath)=2.d0/3.d0
      simpson(ipath+1)=4.d0/3.d0
   enddo
   simpson(npath)=1.d0/3.d0

   if (iproc==0) then
!!$         !start a new document in the beginning of the output, if the document is closed before
      call yaml_set_stream(record_length=95,istat=ierr)
      call yaml_new_document()
!!$         call print_logo()
      call yaml_comment('',hfill='-')
      call yaml_comment('This is a test program to verify  whether the forces are the derivative of the energy: F=-dE/dR')
      call yaml_comment('It performs the integration of the calculated forces over ' // trim(yaml_toa(npath)))
      call yaml_comment(' random displacement (in the range [", -dx, ",",dx,"] a.u.)')
      call yaml_comment(' and compares the result with the difference of the energy')
      call yaml_comment(' between the final and the initial position: E2-E1 = -Integral F.dR"')
      call yaml_comment(' The advantage is two fold:')
      call yaml_comment(' 1) avoiding cancellation error in finite difference derivative,')
      call yaml_comment(' 2) considering the forces over all atoms.')
      call yaml_comment('',hfill='-')
      !print*
      !print*
      !print*,'*******************************************************************************************************'
      !print*,"This is a test program to verify  whether the force are the derivative of the energy: F=-dE/dR"
      !print '(a,i3,a,f7.4,a,f6.4,a)', "It performs the integration of the calculated forces over " , npath,  &
      !   &   " random displacement (in the range [", -dx, ",",dx,"] a.u.)"
      !print*, "and compares the result with the difference of the energy  between the final and the initial position:"// &
      !   &   " E2-E1 = -Integral F.dR" 
      !print*," The advantage is two fold: 1) avoiding cancellation error in finite difference derivative," // & 
      !" 2) considernig the forces over all atoms "
      !print*,'*********************************************************************************************************'
      !print*
      !print*
   endif

   do iconfig=1,abs(nconfig)
      if (modulo(iconfig-1,ngroups)==igroup) then

         ! Read all input files.
         call run_objects_init_from_files(runObj, arr_radical(iconfig),arr_posinp(iconfig))

!!$      !standard names
!!$      call standard_inputfile_names(inputs,radical,nproc)
!!$      call read_input_variables(iproc,nproc,arr_posinp(iconfig),inputs, atoms, rxyz,nconfig,radical,istat)


      !initialize memory counting
      !call memocc(0,iproc,'count','start')
         call init_global_output(outs, runObj%atoms%astruct%nat)

      !     if (iproc == 0) then
      !       call print_general_parameters(nproc,inputs,runObj%atoms)
      !    end if

      !if other steps are supposed to be done leave the last_run to minus one
      !otherwise put it to one
      if (runObj%inputs%last_run == -1 .and. runObj%inputs%ncount_cluster_x <=1 .or. &
           & runObj%inputs%ncount_cluster_x <= 1) then
         runObj%inputs%last_run = 1
      end if

      ! path integral   
      path=0.d0
      !calculate the displacement at each integration step
      !(use sin instead of random numbers)
      allocate(drxyz(1:3,1:runObj%atoms%astruct%nat))
      do iat=1,runObj%atoms%astruct%nat
         drxyz(1,iat)=dx*sin(iat+.2d0)   
         drxyz(2,iat)=dx*sin(iat+.4d0)  
         drxyz(3,iat)=dx*sin(iat+.7d0)  
      end do

      ! loop for ipath 
      do ipath=1,npath

         !update atomic positions along the path
         if(ipath>1) then
            runObj%atoms%astruct%rxyz(:,:)=runObj%atoms%astruct%rxyz(:,:)+drxyz(:,:)
            runObj%inputs%inputPsiId=1
            if(runObj%rst%version == LINEAR_VERSION) runObj%inputs%inputPsiId=101
         end if

         if (iproc == 0) then
            call print_general_parameters(nproc,runObj%inputs,runObj%atoms) ! to know the new positions
         end if

         call call_bigdft(runObj, outs, nproc,iproc,infocode)
         !        inputs%inputPsiId=0   ! change PsiId to 0 if you want to  generate a new Psi and not use the found one

         if (iproc == 0 ) call yaml_map('Wavefunction Optimization Finished, exit signal',infocode)
         !if (iproc == 0 ) write(*,"(1x,a,2i5)") 'Wavefunction Optimization Finished, exit signal=',infocode

         if (ipath == 1 ) etot0=outs%energy
         !   do one step of the path integration
         if (iproc == 0) then
            !integrate forces*displacement
            !fdr=sum(fxyz(1:3,1:runObj%atoms%nat)*drxyz(1:3,1:runObj%atoms%nat))
            fdr=sum(outs%fxyz(:,:)*drxyz(:,:))
            path=path-simpson(ipath)*fdr
            call yaml_map('Path iteration',ipath)
            call yaml_map('-F.dr',-fdr,fmt='(1pe13.5)')
            call yaml_map('Path integral',path,fmt='(1pe13.5)')
            !write(*,"('path iter:',i3,'   -F.dr=',e13.5,'    path integral=',e13.5 )") ipath,-fdr, path 
            
            !Print atomic forces
            call write_forces(runObj%atoms,outs%fxyz)
         end if
      end do !loop over ipath

      deallocate(drxyz)

      if (iproc==0) then 
         write(*,*) 
         write(*,*) 'Check correctness of forces'
         write(*,*) 'Difference of total energies =',outs%energy-etot0
         write(*,*) 'Integral force*displacement = ',path
         write(*,*) 'Difference = ',(outs%energy-etot0)-path
         write(*,*) 
      endif

      call deallocate_global_output(outs)
      call run_objects_free(runObj, subname)

!!$      if (inputs%inputPsiId==INPUT_PSI_LINEAR_AO .or. inputs%inputPsiId==INPUT_PSI_MEMORY_LINEAR &
!!$          .or. inputs%inputPsiId==INPUT_PSI_DISK_LINEAR) then
!!$          call destroy_DFT_wavefunction(rst%tmb)
!!$          call deallocate_local_zone_descriptors(rst%tmb%lzd, subname)
!!$      end if
!!$
!!$      if(inputs%linear /= INPUT_IG_OFF .and. inputs%linear /= INPUT_IG_LIG) &
!!$           & call deallocateBasicArraysInput(inputs%lin)


!!$      call free_input_variables(inputs)

!!$      !finalize memory counting
!!$      call memocc(0,0,'count','stop')
   end if
enddo !loop over iconfig

   deallocate(arr_posinp,arr_radical)

   call bigdft_finalize(ierr)

!!$   call MPI_FINALIZE(ierr)

END PROGRAM test_forces
