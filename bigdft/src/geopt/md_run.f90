!> @file
!!  Routines for the run of MD in NVE and NVT ensembles
!! @author
!!    Copyright (C) 2007-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
!FIXME: get number of SCF cycles and total CPU time for each SCF
subroutine md(run_md,outs,nproc,iproc)

   use module_base
   use bigdft_run

  IMPLICIT NONE
  TYPE(run_objects), intent(inout) :: run_md
  TYPE(state_properties), intent(inout) :: outs 
  integer, intent(in) :: nproc, iproc
  call f_routine(id='md')

!CASE BOMD (NVT,NVE)
  call bomd(run_md,outs,nproc,iproc)
!TODO other cases will be developed soon: CPMD/XL-BOMD (NVT,NVE)
!TODO NPT

  call f_release_routine()
END subroutine md

!> routine to use NH MD files inside BigDFT geometry drivers
subroutine bomd(run_md,outs,nproc,iproc)
  !use nose_hoover_chains_data
  use nose_hoover_chains
  use module_base
  use bigdft_run
  use yaml_output
!NNdbg
   use module_types
  implicit none
!
  TYPE(run_objects), intent(inout) :: run_md
  TYPE(state_properties), intent(inout) :: outs 
  integer, intent(in) :: nproc, iproc
!
  
  INTEGER ierr, iat, jatm, ii, maxsteps, istep, printfrq, ndof

  CHARACTER(LEN=50) :: comment,naming_id
  CHARACTER(LEN=60) :: run_id

  INTEGER :: natoms

  !  TYPE(nhc_type)  :: nhcdata
  !this will become an input variable
  type(NHC_data) :: NHC

  REAL(KIND=8), DIMENSION(:,:), POINTER :: rxyz, fxyz, vxyz
  REAL(KIND=8), DIMENSION(:), POINTER   :: amass
  CHARACTER(LEN=5), DIMENSION(:), POINTER :: alabel

  REAL(KIND=8) :: dist, dd(3), epe, eke, ete, T0ions, Tions, dt, &
       com(3), tcpu0, tcpu1
  REAL(KIND=8), PARAMETER :: amu_to_au=1822.888485D0

  character(len=*), parameter :: subname='bomd'
  LOGICAL :: ionode, no_translation
!NNdbg
  integer :: norbp, iwfn

  call f_routine(id='bomd')

  call cpu_time(tcpu0)
  natoms=bigdft_nat(run_md)

  !Getting a local copy of the coordinates
  rxyz => bigdft_get_rxyz_ptr(run_md)

  ionode=.false.
  IF(iproc==0)ionode=.true. 

  if(ionode) call yaml_comment('Starting BO MD',hfill='*')

  dt= run_md%inputs%dt 
  T0ions=run_md%inputs%temperature
  printfrq=run_md%inputs%md_printfrq
  no_translation=run_md%inputs%no_translation

  call initialize_NHC_data(nhc,run_md%inputs%nhc, &
                                run_md%inputs%nhnc, &
                                run_md%inputs%nsuzuki,run_md%inputs%nmultint,run_md%inputs%nosefrq) 



  if(no_translation)call get_com(natoms,com,rxyz)

  !FIXME
  if(no_translation)then
     ndof=3*natoms-3 !FIXME
  else
     ndof=3*natoms
  end if
 
  vxyz = f_malloc_ptr([3,natoms],id='vxyz')
  fxyz = f_malloc_ptr([3,natoms],id='fxyz')
  amass = f_malloc_ptr(natoms,id='amass')
  alabel = f_malloc_str_ptr(len(alabel),natoms,id='alabel')

  DO iat=1,natoms
     ii=run_md%atoms%astruct%iatype(iat)
     amass(iat)  = run_md%atoms%amu(ii)*amu_to_au
     alabel(iat) = run_md%atoms%astruct%atomnames(ii)
  END DO


  !Allocate initial velocities 

  CALL init_velocities(natoms,3,ndof,amass,T0ions,vxyz,eke)

  Tions=T0ions

  !Initialize Nose Variables
  IF(nhc%nhchain)THEN
     call nose_init(natoms,ndof,dt,T0ions,amass,vxyz,nhc)
     call nose_energy(natoms,ndof,T0ions,nhc)
  END IF

  istep=0

  !Restart 
  !if(ionode)&
  call restart_md('read',iproc,run_md%inputs%restart_pos,&
                         run_md%inputs%restart_vel, &
                         run_md%inputs%restart_nose,&
                         natoms,istep,alabel,&
                         run_md%inputs%dir_output,rxyz,vxyz,nhc)

  if(run_md%inputs%restart_vel)call temperature(natoms,3,ndof,amass,vxyz,Tions,eke)
  if(run_md%inputs%restart_nose)call nose_energy(natoms,ndof,T0ions,nhc)
  !Do SCF if coordinates are restarted
  if(run_md%inputs%restart_pos)call bigdft_state(run_md,outs,ierr)
  epe=outs%energy
  DO iat=1,natoms
     fxyz(1:3,iat)=outs%fxyz(1:3,iat)/amass(iat)
  END DO


  !Total energy
  ete=epe+eke+nhc%enose

  call cpu_time(tcpu1)
  
  !MD printout energies
  IF (ionode)THEN 
     CALL write_md_trajectory(istep,natoms,alabel,rxyz,vxyz)
     CALL write_md_energy(istep,Tions,eke,epe,ete,tcpu1-tcpu0)
     call yaml_comment('Starting MD',hfill='*')
     call yaml_map('Number of degrees of freedom',ndof)
     call yaml_flush_document()
  END IF
  
  maxsteps= run_md%inputs%mdsteps+istep 

  !setting inputpsiid=1 (after the first SCF)
  call bigdft_set_input_policy(INPUT_POLICY_MEMORY, run_md)

  !----------------------------------------------------------------------!
  MD_loop: DO !MD loop starts here

     call cpu_time(tcpu0)

     istep=istep+1
     IF(istep.gt.maxsteps)exit MD_loop

     IF(nhc%NHCHAIN)CALL NOSE_EVOLVE(natoms,ndof,T0ions,amass,vxyz,nhc)

     CALL velocity_verlet_vel(dt,natoms,rxyz,vxyz,fxyz)

     CALL velocity_verlet_pos(dt,natoms,rxyz,vxyz)
     if(no_translation)CALL shift_com(natoms,com,rxyz)



     !> SCF
     CALL bigdft_state(run_md,outs,ierr)
     !FIXME: get number of SCF cycles and total CPU time for each SCF
     epe=outs%energy
     DO iat=1,natoms
        fxyz(1:3,iat)=outs%fxyz(1:3,iat)/amass(iat)
     END DO

     CALL velocity_verlet_vel(dt,natoms,rxyz,vxyz,fxyz)

     IF(nhc%NHCHAIN)THEN
        CALL NOSE_EVOLVE(natoms,ndof,T0ions,amass,vxyz,nhc)
        CALL NOSE_ENERGY(natoms,ndof,T0ions,nhc)
     END IF

     CALL temperature(natoms,3,ndof,amass,vxyz,Tions,eke)

     ete=epe+eke+nhc%enose


     IF (ionode.and.mod(istep,printfrq)==0)& 
          CALL write_md_trajectory(istep,natoms,alabel,rxyz,vxyz)

     if(ionode)&
       call restart_md('write',iproc,run_md%inputs%restart_pos,&
                               run_md%inputs%restart_vel, &
                               run_md%inputs%restart_nose,&
                               natoms,istep,alabel,&
                               run_md%inputs%dir_output,rxyz,vxyz,nhc)
     call cpu_time(tcpu1)
     IF (ionode) & 
          CALL write_md_energy(istep,Tions,eke,epe,ete,tcpu1-tcpu0)

  END DO MD_loop !MD loop ends here
  !----------------------------------------------------------------------!

  !the deallocation of the pointers
  call finalize_NHC_data(nhc)
  call f_free_ptr(vxyz)
  call f_free_ptr(fxyz)
  call f_free_ptr(amass)
  call f_free_str_ptr(len(alabel),alabel)

  call f_release_routine()

end subroutine bomd

!>Initialize velocities for MD using Box-Muller Sampling
SUBROUTINE init_velocities(natoms,ndim,ndof,amass,T0ions,vxyz,eke)
  use numerics, only: pi,au_to_k => Ha_K
  use dynamic_memory, only: f_release_routine,f_routine
  use random, only: builtin_rand
  IMPLICIT NONE
  INTEGER, INTENT(IN):: natoms, ndof, ndim
  REAL(KIND=8), INTENT(IN) :: T0ions, amass(natoms)
  REAL(KIND=8), INTENT(OUT)::  vxyz(ndim,natoms),eke  
  !
  REAL(KIND=8) :: sigma, dum(2)
  INTEGER      :: iat, k
  INTEGER      :: idum=0
  
  call f_routine(id='init_velocities')

  !FIXME: a proper ndof has to be determined
  DO iat=1,natoms,2
     DO k=1,ndim
        dum(1)=real(builtin_rand(idum),8)
        dum(2)=real(builtin_rand(idum),8)
        sigma=SQRT(T0ions/au_to_k/amass(iat)) !Sqrt(kT/M_i) 
        vxyz(k,iat)=sigma*SQRT(-2.D0*LOG(dum(2)))*COS(2.D0*pi*dum(1))
        IF(iat+1.LE.natoms)then 
           sigma=SQRT(T0ions/au_to_k/amass(iat+1)) 
           vxyz(k,iat+1)=sigma*SQRT(-2.D0*LOG(dum(2)))*SIN(2.D0*pi*dum(1))
        END IF
     END DO
  END DO

  CALL remove_lin_momentum(natoms,vxyz)

  CALL rescale_velocities(natoms,ndim,ndof,amass,T0ions,vxyz,eke)
  !
  call f_release_routine()

END SUBROUTINE init_velocities


!> Subroutine to Rescale the velocities to temperature Ttarg
SUBROUTINE rescale_velocities(natoms,ndim,ndof,amass,Ttarg,vxyz,eke)
  IMPLICIT NONE
  !
  INTEGER :: natoms,ndof,ndim
  REAL(KIND=8)  :: Ttarg, vxyz(ndim,*), amass(*), eke
  !
  INTEGER:: iat
  REAL*8 :: Tcur,SCAL
  !
  CALL temperature(natoms,ndim,ndof,amass,vxyz,Tcur,eke)
  scal=DSQRT(Ttarg/Tcur)
  DO iat=1,natoms
     vxyz(1:ndim,iat)=vxyz(1:ndim,iat)*scal
  END DO
  eke=eke*scal*scal
END SUBROUTINE rescale_velocities


!> Subroutine to Compute Instantaneous Temperature Tinst
SUBROUTINE temperature(natoms,ndim,ndof,amass,vxyz,Tinst,eke)
  use numerics, only: au_to_k => Ha_K
  IMPLICIT NONE
  INTEGER :: natoms, ndim, ndof
  REAL(KIND=8) :: vxyz(ndim,natoms), Tinst, amass(natoms), eke 
  !
  INTEGER :: iat, k
  REAL(KIND=8) :: mv2

  !FIXME: ndof should be properly computed (constraints, bulk, gasphase etc.)
  !  ndof=3*natoms-3
  !
  mv2=0.D0
  DO iat=1,natoms
     DO k=1,ndim
        mv2=mv2+amass(iat)*vxyz(k,iat)*vxyz(k,iat)
     END DO
  END DO
  eke=0.5d0*mv2
  Tinst=mv2*au_to_k/real(ndof,kind=8)
END SUBROUTINE temperature



!>Update positions
SUBROUTINE velocity_verlet_pos(dt,natoms,rxyz,vxyz)
  IMPLICIT NONE
  !
  REAL(KIND=8) :: dt, rxyz(3,*), vxyz(3,*)
  INTEGER :: natoms
  !
  INTEGER :: iat
  !
  DO iat=1,natoms
     rxyz(1,iat)=rxyz(1,iat)+vxyz(1,iat)*dt
     rxyz(2,iat)=rxyz(2,iat)+vxyz(2,iat)*dt
     rxyz(3,iat)=rxyz(3,iat)+vxyz(3,iat)*dt
  END DO
END SUBROUTINE velocity_verlet_pos



!>Update velocities
SUBROUTINE velocity_verlet_vel(dt,natoms,rxyz,vxyz,fxyz)
  IMPLICIT NONE
  !
  REAL(KIND=8) :: dt, rxyz(3,*), vxyz(3,*), fxyz(3,*)
  INTEGER :: natoms
  !
  INTEGER :: iat
  !
  DO iat=1,natoms
     vxyz(1,iat)=vxyz(1,iat)+0.5D0*dt*fxyz(1,iat)
     vxyz(2,iat)=vxyz(2,iat)+0.5D0*dt*fxyz(2,iat)
     vxyz(3,iat)=vxyz(3,iat)+0.5D0*dt*fxyz(3,iat)
  END DO
END SUBROUTINE velocity_verlet_vel

SUBROUTINE write_md_trajectory(istep,natoms,alabel,rxyz,vxyz)
  use f_utils
  IMPLICIT NONE
  INTEGER :: istep, natoms
  REAL(KIND=8) :: rxyz(3,*), vxyz(3,*)
  CHARACTER(LEN=*) :: alabel(*)
  !
  INTEGER :: iat,unt
  unt=111
  call f_open_file(unt,FILE='Trajectory.xyz',status='UNKNOWN',position='APPEND',binary=.false.)
  WRITE(unt,*)natoms
  WRITE(unt,'(A,I16)')'Step:',istep
  DO iat=1,natoms
     WRITE(unt,'(A,6f16.6)')alabel(iat),rxyz(1:3,iat)*0.529,vxyz(1:3,iat)*0.529
  END DO
  call f_close(unt)
END SUBROUTINE write_md_trajectory

SUBROUTINE write_md_energy(istep,Tions,eke,epe,ete,tcpu)
  use yaml_output
  use yaml_strings
  use f_utils
  IMPLICIT NONE
  INTEGER :: istep
  REAL(KIND=8) :: Tions, eke, epe, ete,tcpu
  !local variables
  character(len=*), parameter :: fm='(f16.6)', fm2='(f10.2)'
  integer :: unt
  unt = 111
  !IF(istep.eq.0)PRINT "(2X,A,5A16)", "(MD)","ISTEP","TEMP.","EKE","EPE","ETE"
  !PRINT "(2X,A,I16,4F16.6)","(MD)",istep,Tions,eke,epe,ete
  call yaml_mapping_open("(MD)",flow=.true.)
  call yaml_map('istep',istep,fmt='(i6)')
  call yaml_map('T',Tions,fmt=fm2)
  call yaml_map('Eke',eke,fmt=fm)
  call yaml_map('Epe',epe,fmt=fm)
  call yaml_map('Ete',ete,fmt=fm)
  call yaml_map('tcpu',tcpu,fmt=fm2)
  call yaml_mapping_close()
  call f_open_file(unt,FILE='energy.dat',STATUS='UNKNOWN',position='APPEND',binary=.false.)
  WRITE(unt,"(I16,F10.2,3F16.6,F10.2)")istep,Tions,eke,epe,ete,tcpu
  call f_close(unt)
END SUBROUTINE write_md_energy

SUBROUTINE remove_lin_momentum(natoms,vxyz)
  IMPLICIT NONE
  INTEGER :: natoms
  REAL(KIND=8) :: vxyz(3,natoms)
  !
  REAL(KIND=8) :: rlm(natoms)
  INTEGER :: iat

  rlm(1:3)=0.d0
  DO iat=1,natoms
     rlm(1)=rlm(1) + vxyz(1,iat)
     rlm(2)=rlm(2) + vxyz(2,iat)
     rlm(3)=rlm(3) + vxyz(3,iat)
  END DO
  rlm(1:3)=rlm(1:3)/real(natoms,kind=8)
  DO iat=1,natoms
     vxyz(1,iat)=vxyz(1,iat) - rlm(1)
     vxyz(2,iat)=vxyz(2,iat) - rlm(2)
     vxyz(3,iat)=vxyz(3,iat) - rlm(3)
  END DO
END SUBROUTINE remove_lin_momentum

SUBROUTINE get_com(natoms,com,rxyz)
  IMPLICIT NONE
  INTEGER :: natoms
  REAL(KIND=8) :: rxyz(3,natoms), com(3)
  !
  INTEGER :: iat, k
  DO k=1,3
     com(k) = 0.d0
     DO iat=1,natoms
        com(k)=com(k) + rxyz(k,iat)
     END DO
     com(k)=com(k)/dble(natoms) 
  END DO
END SUBROUTINE get_com

SUBROUTINE shift_com(natoms,com,rxyz)
  IMPLICIT NONE
  INTEGER :: natoms
  REAL (KIND=8) :: com(3), rxyz(3,natoms)

  INTEGER :: iat, k
  REAL (KIND=8) :: newcom(3)

  CALL get_com(natoms,newcom,rxyz)

  DO iat=1,natoms
     DO k=1,3
        rxyz(k,iat)=rxyz(k,iat)-(newcom(k)-com(k))
     END DO
  END DO
END SUBROUTINE shift_com

SUBROUTINE restart_md(control,iproc,restart_pos,restart_vel,restart_nose,natoms,istep,&
                      alabel,dir,rxyz,vxyz,nhc)
  use f_utils
  use nose_hoover_chains
  use yaml_output
  implicit none
  integer :: istep, natoms, iproc
  real(kind=8) :: rxyz(3,*), vxyz(3,*)
  character(len=*) :: alabel(*), dir, control
  type(NHC_data), intent(inout) :: nhc
  logical :: restart_pos, restart_vel, restart_nose
  !
  integer :: iat,ierr, ios
  integer :: unt
  real (kind=8) :: dum(3)
  !real (kind=8), parameter :: bohr_to_ang=0.529d0
  character (len=20) :: dummy_char

  unt=113


  select case (trim(control))
  case ('write')
    !only done at iproc==0
    if(iproc==0)then
      call f_open_file(unt,file=trim(dir)//'md.restart',status='UNKNOWN',position='REWIND',binary=.false.)
      call yaml_map('MD Restart file opened for writing',trim(dir)//'md.restart')
      !write coordinates in xyz format for easy editing by the users
      write(unt,*)natoms
      write(unt,'(A,I16)')'Step:',istep
      do iat=1,natoms
         write(unt,'(A,6e16.8)')alabel(iat),rxyz(1:3,iat),vxyz(1:3,iat)
      end do
      call write_nhc_restart(unt,nhc)
      call f_close(unt)
    end if
!
  case ('read')
    !done by all the processors
    call f_open_file(unt,file=trim(dir)//'md.restart',status='UNKNOWN',position='REWIND',binary=.false.)
    if(iproc==0)call yaml_map('MD Restart file opened for reading',trim(dir)//'md.restart')
    if(restart_pos.and.restart_vel)then
      !read positions and velocities
      read(unt,*)iat
      if(iat/=natoms)&
      call f_err_throw('Error reading md.restart! wrong number of atoms', &
                       err_name='BIGDFT_RUNTIME_ERROR')
      read(unt,*)dummy_char,istep
      do iat=1,natoms
         read(unt,*,iostat=ios)alabel(iat),rxyz(1:3,iat),vxyz(1:3,iat)
      end do
      if(ios/=0) call f_err_throw('Error reading pos and vel from md.restart!', &
                       err_name='BIGDFT_RUNTIME_ERROR')
      if(iproc==0)then
        call yaml_map('Initial positions  restarted from step',istep)
        call yaml_map('Initial velocities restarted from step',istep)
      end if
    else if(restart_pos)then
      read(unt,*,iostat=ios)iat
      if(iat/=natoms)&
        call f_err_throw('Error reading md.restart! wrong number of atoms', &
                         err_name='BIGDFT_RUNTIME_ERROR')
      read(unt,*,iostat=ios)dummy_char,istep
      do iat=1,natoms
        read(unt,*,iostat=ios)alabel(iat),rxyz(1:3,iat)
      end do
      if(ios/=0) call f_err_throw('Error reading pos from md.restart!', &
                         err_name='BIGDFT_RUNTIME_ERROR')
      if(iproc==0)then
        call yaml_map('Initial positions restarted from step',istep)
      end if
    else if(restart_vel)then
      !read velocities 
      read(unt,*,iostat=ios)iat
      if(iat/=natoms)&
        call f_err_throw('Error reading md.restart! wrong number of atoms', &
                       err_name='BIGDFT_RUNTIME_ERROR')
      read(unt,*,iostat=ios)dummy_char,istep
      do iat=1,natoms
         read(unt,*,iostat=ios)alabel(iat),dum(1:3),vxyz(1:3,iat)
      end do
      if(ios/=0) call f_err_throw('Error reading vel from md.restart!', &
                         err_name='BIGDFT_RUNTIME_ERROR')
      if(iproc==0)&
        call yaml_map('Initial velocities restarted from step',istep)
    end if
    if(restart_nose)then
      rewind(unt)
      read(unt,*,iostat=ios)iat
      read(unt,*,iostat=ios)dummy_char,istep 
      do iat=1,natoms 
        read(unt,*,iostat=ios) ! skip atomic positions and velocities from the restart
      end do
      if(ios/=0) call f_err_throw('Error skipping data while reading nose info from md.restart!', &
                       err_name='BIGDFT_RUNTIME_ERROR')
      call read_nhc_restart(unt,nhc,ierr)
      if(ierr==0.and.iproc==0) &
        call yaml_map('Initial Nose Hoover Chains restarted from step',istep)
      if(ierr/=0.and.iproc==0)then
         call yaml_warning('Nose Hoover information is absent in the md.restart file')
         call yaml_warning('Nose Hoover information is not restarted')
      end if
    end if
    call f_close(unt)
  end select
END SUBROUTINE restart_md

