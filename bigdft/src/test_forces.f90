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
!!    Date: Jun 2016:  This program was modified and updated by Deb Sankar De
program test_forces

   use module_base
   use bigdft_run
   !use public_enums, only: SPHERICAL!LINEAR_VERSION
   !use module_interfaces
   !use m_ab6_symmetry
   use yaml_output
   use f_random
   use module_input_keys, only: print_general_parameters
   implicit none
   character(len=*), parameter :: subname='test_forces'
   !integer :: iproc,nproc,
   integer :: iat,ierr,infocode!,istat
   !logical :: exist_list
   !input variables
   type(run_objects) :: runObj
   type(state_properties) :: outs
   !character(len=60), parameter :: filename="list_posinp"
   ! atomic coordinates, forces
   real(gp), dimension(:,:), pointer :: drxyz,dr,int0
   !integer :: iconfig,nconfig,igroup,ngroups
   integer :: ipath,npath,natx
   real(gp):: dx,etot0,path,fdr,stepsize,fact,t1,t2,t3
   !parameter (dx=1.d-2 , npath=2*16+1)  ! npath = 2*n+1 where n=2,4,6,8,...
   parameter (dx=1.d-2 , npath=500,natx=1000)
   real(gp) :: simpson(1:npath),dd(3,natx,2)
   !character(len=60) :: radical
   !integer, dimension(4) :: mpi_info
   type(dictionary), pointer :: run,options
   character(len = max_field_length) :: input_id
   real(gp) ::e_min, e_max
      

   call f_lib_initialize()

   call bigdft_command_line_options(options)
   call bigdft_init(options)

   !prepare the array of the correct Simpson's rule weigths for the integration
   if (mod(npath,2).eq.1) stop 'the number of iteration steps has to be even' !.ne.1) stop 'the number of iteration steps has to be odd'
   simpson(1)=1.d0/3.d0
   simpson(2)=4.d0/3.d0
   do ipath=3,npath-2,2
      simpson(ipath)=2.d0/3.d0
      simpson(ipath+1)=4.d0/3.d0
   enddo
   simpson(npath)=1.d0/3.d0
   stepsize=1.0d0/npath
   fact=2*pi/npath
   etot0=0.0_gp

   if (bigdft_mpi%iproc==0) then
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
   endif


   !call random_number(dd)
   call f_random_number(dd)
   dd=dd*.1d0

   run => dict_iter(options .get. 'BigDFT')
   do while(associated(run))
!!$
!!$   do iconfig=1,abs(nconfig)
!!$      if (modulo(iconfig-1,ngroups)==igroup) then

         ! Read all input files.
         call run_objects_init(runObj,run)

!!$      !standard names
!!$      call standard_inputfile_names(inputs,radical,nproc)
!!$      call read_input_variables(iproc,nproc,arr_posinp(iconfig),inputs, atoms, rxyz,nconfig,radical,istat)


      !initialize memory counting
         call init_state_properties(outs, runObj%atoms%astruct%nat)

      !     if (iproc == 0) then
      !       call print_general_parameters(inputs,atoms)
      !    end if

      !if other steps are supposed to be done leave the last_run to minus one
      !otherwise put it to one
      if (runObj%inputs%last_run == -1 .and. runObj%inputs%ncount_cluster_x <=1 .or. &
           & runObj%inputs%ncount_cluster_x <= 1) then
         runObj%inputs%last_run = 1
      end if

      ! path integral   
      path=0.d0
      e_min=1.d100
      e_max=-1.d100
      !calculate the displacement at each integration step
      !(use sin instead of random numbers)
      
      !allocate(drxyz(1:3,1:runObj%atoms%astruct%nat),dr(1:3,1:runObj%atoms%astruct%nat),int0(1:3,1:runObj%atoms%astruct%nat))
      drxyz=f_malloc_ptr([3,runObj%atoms%astruct%nat],id='drxyz')
      dr=f_malloc_ptr([3,runObj%atoms%astruct%nat],id='dr')
      int0=f_malloc_ptr([3,runObj%atoms%astruct%nat],id='int0')
!     do iat=1,runObj%atoms%astruct%nat
!        drxyz(1,iat)=  !dx*sin(iat+.2d0)   
!        drxyz(2,iat)=  !dx*sin(iat+.4d0)  
!        drxyz(3,iat)=  !dx*sin(iat+.7d0)  
!     end do



!dd(1,1,1)=9.9755959009261722E-002
!dd(2,1,1)=5.6682470761127340E-002
!dd(3,1,1)=9.6591537549612499E-002
!dd(1,1,2)=7.3754263633984520E-003
!dd(2,1,2)=5.3552292777272472E-004
!dd(3,1,2)=3.4708128851801544E-002
!dd(1,2,1)=7.4792768547143215E-002
!dd(2,2,1)=3.6739089737475576E-002
!dd(3,2,1)=4.8063689875473152E-002
!dd(1,2,2)=3.4224381607283506E-002
!dd(2,2,2)=2.1795172633847261E-002
!dd(3,2,2)=1.3316041001365931E-002



int0(:,:)=runObj%atoms%astruct%rxyz(:,:)
!DEBug write(*,*)int0 
!stop
      ! loop for ipath 
      do ipath=0,npath

         !update atomic positions along the path
!         if(ipath.ge.0) then

     do iat=1,runObj%atoms%astruct%nat
        drxyz(1,iat)=dd(1,iat,1)*cos(2*pi*stepsize*ipath)-dd(1,iat,2)*sin(2*pi*stepsize*ipath)
        drxyz(2,iat)=dd(2,iat,1)*cos(2*pi*stepsize*ipath)-dd(2,iat,2)*sin(2*pi*stepsize*ipath)   
        drxyz(3,iat)=dd(3,iat,1)*cos(2*pi*stepsize*ipath)-dd(3,iat,2)*sin(2*pi*stepsize*ipath)
           
!       dr(1,iat)=dd(1,iat,1)*cos(2*pi*stepsize*ipath)+dd(1,iat,2)*sin(2*pi*stepsize*ipath)
!       dr(2,iat)=dd(2,iat,1)*cos(2*pi*stepsize*ipath)+dd(2,iat,2)*sin(2*pi*stepsize*ipath)   
!       dr(3,iat)=dd(3,iat,1)*cos(2*pi*stepsize*ipath)+dd(3,iat,2)*sin(2*pi*stepsize*ipath)
        runObj%atoms%astruct%rxyz(1,iat)=int0(1,iat)+dd(1,iat,1)*sin(2*pi*stepsize*ipath)+dd(1,iat,2)*cos(2*pi*stepsize*ipath)
        runObj%atoms%astruct%rxyz(2,iat)=int0(2,iat)+dd(2,iat,1)*sin(2*pi*stepsize*ipath)+dd(2,iat,2)*cos(2*pi*stepsize*ipath)   
        runObj%atoms%astruct%rxyz(3,iat)=int0(3,iat)+dd(3,iat,1)*sin(2*pi*stepsize*ipath)+dd(3,iat,2)*cos(2*pi*stepsize*ipath)
     end do
!            runObj%atoms%astruct%rxyz(:,:)=int0(:,:)+dr(:,:) !runObj%atoms%astruct%rxyz(:,:)+dr(:,:)
!            write(*,*)runObj%atoms%astruct%rxyz(:,:)
!!$            runObj%inputs%inputPsiId=1
!!$            if(runObj%rst%version == LINEAR_VERSION) then
!!$               runObj%inputs%inputPsiId=101
!!$               !switch off fragment calculation after this point
!!$               !runObj%inputs%lin%fragment_calculation=.false.
!!$               !runObj%inputs%frag%nfrag=1
!!$            end if
!            runObj%atoms%astruct%rxyz(:,:)=int0(:,:)+dr(:,:)
            call bigdft_set_input_policy(INPUT_POLICY_MEMORY, runObj)
!         end if

!            runObj%atoms%astruct%rxyz(:,:)=int0(:,:)+dr(:,:)
!DEBug            write(7,*)runObj%atoms%astruct%nat
!DEBug            write(7,*)
!DEBug            do iat=1,runObj%atoms%astruct%nat
!            write(7,'(1X)')'free'
!DEBug           write(7,'(A,3(2x,E14.6))')runObj%atoms%astruct%atomnames(runObj%atoms%astruct%iatype(iat)),runObj%atoms%astruct%rxyz(1,iat),runObj%atoms%astruct%rxyz(2,iat),runObj%atoms%astruct%rxyz(3,iat)
!DEBug           enddo
            call bigdft_state(runObj, outs,infocode)



         if (bigdft_mpi%iproc == 0) then
            call bigdft_get_run_properties(run, input_id = input_id)
            call print_general_parameters(runObj%inputs,runObj%atoms,input_id) ! to know the new positions
         end if

!         call bigdft_state(runObj, outs,infocode)
         !        inputs%inputPsiId=0   ! change PsiId to 0 if you want to  generate a new Psi and not use the found one

         if (bigdft_mpi%iproc == 0 ) call yaml_map('Wavefunction Optimization Finished, exit signal',infocode)
         !if (iproc == 0 ) write(*,"(1x,a,2i5)") 'Wavefunction Optimization Finished, exit signal=',infocode

!         if (ipath == 1 ) etot0=outs%energy
         !   do one step of the path integration
!        if (bigdft_mpi%iproc == 0) then
!           !integrate forces*displacement
!           !fdr=sum(fxyz(1:3,1:runObj%atoms%nat)*drxyz(1:3,1:runObj%atoms%nat))
!           fdr=sum(outs%fxyz(:,:)*drxyz(:,:))
!           path=path-simpson(ipath)*fdr
           
!           runObj%atoms%astruct%rxyz(:,:)=int0(:,:)+dr(:,:)
!           write(*,*)runObj%atoms%astruct%rxyz(:,:)


     if (ipath.eq.0) etot0=outs%energy

      e_min=min(e_min,outs%energy)
      e_max=max(e_max,outs%energy)
!DEBug       do iat=1,runObj%atoms%astruct%nat
!DEBug!       write(8,*)outs%fxyz(1,iat),drxyz(1,iat),outs%fxyz(2,iat),drxyz(2,iat),outs%fxyz(3,iat),drxyz(3,iat)
!DEBug       write(8,*)drxyz(1,iat)*outs%fxyz(1,iat),drxyz(2,iat)*outs%fxyz(2,iat),drxyz(3,iat)*outs%fxyz(3,iat)
!DEBug       enddo
      if (ipath < npath) then
        t1=0.d0
        t2=0.d0
        t3=0.d0
        do iat=1,runObj%atoms%astruct%nat
           t1=t1+outs%fxyz(1,iat)*drxyz(1,iat)
           t2=t2+outs%fxyz(2,iat)*drxyz(2,iat)
           t3=t3+outs%fxyz(3,iat)*drxyz(3,iat)
!          write(8,*)ipath
        enddo
!DEBug        write(18,*)etot0,outs%energy
        fdr=(t1+t2+t3)*fact 
        path=path+fdr
!DEBug      write(11,'(i5,2(2x,e17.10))') ipath,outs%energy-etot0,path

            call yaml_map('(Test forces) Results',[outs%energy-etot0,path],fmt='(e17.10)')
            call yaml_map('Path iteration',ipath)
            call yaml_map('-F.dr',-fdr,fmt='(1pe13.5)')
            call yaml_map('Path integral',path,fmt='(1pe13.5)')
            !write(*,"('path iter:',i3,'   -F.dr=',e13.5,'    path integral=',e13.5 )") ipath,-fdr, path 
            
            !Print atomic forces
            call write_forces(bigdft_get_astruct_ptr(runObj),outs%fxyz)
         end if
      end do !loop over ipath

      !deallocate(drxyz,dr,int0)
      call f_free_ptr(drxyz)
      call f_free_ptr(dr)
      call f_free_ptr(int0)

      if (bigdft_mpi%iproc==0) then 
         write(*,*) 
         write(*,*) 'Check correctness of forces'
         write(*,*) 'Difference of total energies =',outs%energy-etot0
         write(*,*) 'Integral force*displacement = ',path
         write(*,*) 'Difference = ',(outs%energy-etot0)-path
         write(*,*) 
      endif

      call deallocate_state_properties(outs)
      call free_run_objects(runObj)

!!$   end if
      run => dict_next(run)
   end do !loop over iconfig

!!$   deallocate(arr_posinp,arr_radical)
   call dict_free(options)
   call bigdft_finalize(ierr)

   call f_lib_finalize()

END PROGRAM test_forces
