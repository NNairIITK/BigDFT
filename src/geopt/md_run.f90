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
       com(3)
  REAL(KIND=8), PARAMETER :: amu_to_au=1822.888485D0

  character(len=*), parameter :: subname='bomd'
  LOGICAL :: ionode, no_translation

  call f_routine(id='bomd')

  natoms=bigdft_nat(run_md)

  !Getting a local copy of the coordinates
  rxyz => bigdft_get_rxyz_ptr(run_md)

  ionode=.false.
  IF(iproc==0)ionode=.true. 

  maxsteps= run_md%inputs%mdsteps
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
!!not permitted directly, see http://bigdft.org/Wiki/index.php?title=Coding_Rules#Low_level_operations
!!$  ALLOCATE(vxyz(3,natoms))
!!$  ALLOCATE(fxyz(3,natoms))
!!$  ALLOCATE(amass(natoms))
!!$  ALLOCATE(alabel(natoms))

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


  epe=outs%energy
  DO iat=1,natoms
     fxyz(1:3,iat)=outs%fxyz(1:3,iat)/amass(iat)
  END DO


  !Total energy
  ete=epe+eke+nhc%enose
  !  ete=epe+eke

  istep=0


  !MD printout energies
  IF (ionode)THEN 
     CALL write_md_trajectory(istep,natoms,alabel,rxyz,vxyz)
     CALL write_md_energy(istep,Tions,eke,epe,ete)
     call yaml_comment('Starting MD',hfill='*')
     call yaml_map('Number of degrees of freedom',ndof)
!     call yaml_map('Maximum number of steps (maxsteps)',maxsteps)
!     call yaml_map('Initial Temperature (T0ions)',T0ions)
  END IF


  !----------------------------------------------------------------------!
  MD_loop: DO !MD loop starts here

     istep=istep+1

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
!     call yaml_map('enose', nhc%enose)
!     call yaml_map('istep',istep)
     !print *, "enose =", nhc%enose, "istep=",istep

     ete=epe+eke+nhc%enose
     !    ete=epe+eke

     IF (ionode) & 
          CALL write_md_energy(istep,Tions,eke,epe,ete)

     IF (ionode.and.mod(istep,printfrq)==0)& 
          CALL write_md_trajectory(istep,natoms,alabel,rxyz,vxyz)

     IF(istep+1.gt.maxsteps)exit MD_loop
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
  IMPLICIT NONE
  INTEGER, INTENT(IN):: natoms, ndof, ndim
  REAL(KIND=8), INTENT(IN) :: T0ions, amass(natoms)
  REAL(KIND=8), INTENT(OUT)::  vxyz(ndim,natoms),eke  
  !
  REAL(KIND=8) :: sigma, dum(2)
  INTEGER      :: iat, k
  REAL(KIND=4) :: builtin_rand
  INTEGER      :: idum=0

  !FIXME: a proper ndof has to be determined
  DO iat=1,natoms,2
     DO k=1,ndim
        dum(1)=real(builtin_rand(idum),8)
        dum(2)=real(builtin_rand(idum),8)
        sigma=DSQRT(T0ions/au_to_k/amass(iat)) !Sqrt(kT/M_i) 
        vxyz(k,iat)=sigma*DSQRT(-2.D0*DLOG(dum(2)))*DCOS(2.D0*pi*dum(1))
        IF(iat+1.LE.natoms)then 
           sigma=DSQRT(T0ions/au_to_k/amass(iat+1)) 
           vxyz(k,iat+1)=sigma*DSQRT(-2.D0*DLOG(dum(2)))*DSIN(2.D0*pi*dum(1))
        END IF
     END DO
  END DO

  CALL remove_lin_momentum(natoms,vxyz)

  CALL rescale_velocities(natoms,ndim,ndof,amass,T0ions,vxyz,eke)
  !
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
  Tinst=mv2*au_to_k/DFLOAT(ndof)
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

SUBROUTINE write_md_energy(istep,Tions,eke,epe,ete)
  use yaml_output
  use yaml_strings
  use f_utils
  IMPLICIT NONE
  INTEGER :: istep
  REAL(KIND=8) :: Tions, eke, epe, ete
  !local variables
  character(len=*), parameter :: fm='(f16.6)'
  integer :: unt
  unt = 111
  !IF(istep.eq.0)PRINT "(2X,A,5A16)", "(MD)","ISTEP","TEMP.","EKE","EPE","ETE"
  !PRINT "(2X,A,I16,4F16.6)","(MD)",istep,Tions,eke,epe,ete
  call yaml_mapping_open("(MD)",flow=.true.)
  call yaml_map('istep',istep,fmt='(i6)')
  call yaml_map('T',Tions,fmt=fm)
  call yaml_map('Eke',eke,fmt=fm)
  call yaml_map('Epe',epe,fmt=fm)
  call yaml_map('Ete',ete,fmt=fm)
  call yaml_mapping_close()
  call f_open_file(unt,FILE='energy.dat',STATUS='UNKNOWN',position='APPEND',binary=.false.)
  WRITE(unt,"(I16,4F16.6)")istep,Tions,eke,epe,ete
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
  rlm(1:3)=rlm(1:3)/dfloat(natoms)
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

