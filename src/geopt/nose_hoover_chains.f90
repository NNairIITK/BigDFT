!> @file
!!  Routines for MD
!! @author
!!    Nisanth Nair
!!    Copyright (C) 2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
MODULE Nose_Hoover_Chains
  use module_base
  implicit none

  type, public :: NHC_data
     LOGICAL          :: NHCHAIN, NOSPEAT 
     INTEGER          :: ETA_SIZE, NHNC, NMULTINT, NSUZUKI
     REAL (KIND=8)    :: NOSEFRQ, ENOSE
     REAL (KIND=8), DIMENSION(:), pointer :: ETA, VETA, AETA, NOSEQ, DTYS
  end type NHC_data

CONTAINS

  pure subroutine nullify_NHC_data(NHC)
    implicit none
    type(NHC_data), intent(out) :: NHC
    NHC%NHCHAIN=.false.
    NHC%NOSPEAT=.false.
    NHC%ETA_SIZE=0
    NHC%NHNC=0
    NHC%NMULTINT=0
    NHC%NSUZUKI=0
    NHC%NOSEFRQ=0.d0
    NHC%ENOSE=0.d0
    nullify(NHC%ETA)
    nullify(NHC%VETA)
    nullify(NHC%AETA)
    nullify(NHC%NOSEQ)
    nullify(NHC%DTYS)
  end subroutine nullify_NHC_data

  subroutine initialize_NHC_data(nhc,nhchain,nospeat,nhnc,nsuzuki,nmultint,nosefrq)
    implicit none
    logical, intent(in)  :: NHCHAIN !<Nose Hoover Chain thermostat 
    logical, intent(in) :: NOSPEAT !<!separate Nose hoover chain for all the atoms
    INTEGER, intent(in)  :: NHNC, NMULTINT, NSUZUKI
    REAL (KIND=8), intent(in) :: NOSEFRQ
    type(NHC_data), intent(out) :: NHC

    !perform checks
    if (all(nsuzuki /= [3,5,7])) then
       call f_err_throw('Error in the nsuzuki parameter, must be 3,5 or 7, found'//nsuzuki,&
            err_name='BIGDFT_INPUT_VARIABLES_ERROR')
    end if

    call nullify_NHC_data(NHC)
    NHC%NHCHAIN=NHCHAIN
    NHC%NOSPEAT=NOSPEAT
    NHC%NHNC=NHNC
    NHC%NSUZUKI=NSUZUKI
    NHC%NMULTINT=NMULTINT
    NHC%NOSEFRQ=NOSEFRQ

  end subroutine initialize_NHC_data

  !
  !     -----------------------------------------------------
  !>     Function : Nose-Hoover-chains are initialized 
  !!     Note     : Units used for Nose-Hoover chain integration 
  !!                is a.u. 
  !!     -----------------------------------------------------         
  !!
  SUBROUTINE NOSE_INIT(natoms,ndof,dt,T0ions,amass,vxyz,nhc)
    !use nose_hoover_chains_data
    IMPLICIT NONE 
    INTEGER :: natoms,ndof
    REAL (KIND=8)   :: dt, T0ions, amass(natoms),vxyz(3,natoms) 
    type(NHC_data), intent(inout) :: nhc
    !local variables
    INTEGER         :: I, II, INHC, IATOMS, DOF
    LOGICAL         :: TO_RESTART
    REAL (KIND=8)   :: GAUSS, DUMMY, DUM(2) 
    REAL (KIND=8)   :: KT, FKT, AKIN, sigma 
    REAL (KIND=8)   :: WYS1, WYS2, WYS3, WYS4

    !
    !
    REAL (KIND=8), PARAMETER  :: AU_TO_K   = 315774.664550534774D0, & 
         CMI_TO_AU = 7.26D-7, &
         pi        = 4.d0*ATAN(1.D0),  &
         CUBROOT   = 1.D0/3.D0

    nhc%ETA_SIZE = nhc%NHNC
    IF(nhc%NOSPEAT)nhc%ETA_SIZE = nhc%ETA_SIZE * NATOMS

    !ETA_SIZE=ETA_SIZE

    NHC%NOSEQ=f_malloc0_ptr(nhc%ETA_SIZE,id='noseq')
    NHC%ETA=f_malloc0_ptr(nhc%ETA_SIZE,id='eta')
    NHC%VETA=f_malloc0_ptr(nhc%ETA_SIZE,id='veta')
    NHC%AETA=f_malloc0_ptr(nhc%ETA_SIZE,id='aeta')
    NHC%DTYS=f_malloc0_ptr(7,id='dtys')
    !
    DOF = ndof
    IF(nhc%NOSPEAT) DOF = 3
    !
    !   ** To a.u **
    !
    KT  = T0ions / AU_TO_K
    FKT = KT * REAL(DOF,gp)
    !
    !   ** Nose freq. in a.u **
    !
    nhc%NOSEFRQ = nhc%NOSEFRQ * CMI_TO_AU
    !
    !   ** Nose masses in a.u **
    !
    print *, "done nose init state 0 "
    IF(.NOT. nhc%NOSPEAT)THEN
       nhc%NOSEQ(1) = FKT / nhc%NOSEFRQ**2    
       DO INHC = 2, nhc%NHNC
          nhc%NOSEQ(INHC) = KT / nhc%NOSEFRQ**2  
       END DO
    ELSE
       DO IATOMS = 1, NATOMS
          II = (IATOMS - 1) * nhc%NHNC + 1
          nhc%NOSEQ(II) = FKT / nhc%NOSEFRQ**2
          DO INHC = 2, nhc%NHNC
             II = (IATOMS - 1) * nhc%NHNC + INHC
             nhc%NOSEQ(II) = KT / nhc%NOSEFRQ**2
          END DO
       END DO
    END IF
    !
    !  ** Restart Nose Hoover Chain **
    !  FIXME

    !   ** Assign initial nose velocities **
    !
    IF(.NOT. nhc%NOSPEAT)THEN
       DO INHC = 1, nhc%NHNC, 2
          sigma=DSQRT(KT/nhc%NOSEQ(inhc)) !Sqrt(kT/M_i) 
          CALL RANDOM_NUMBER(dum(1))
          CALL RANDOM_NUMBER(dum(2))
          nhc%VETA(INHC) =sigma*DSQRT(-2.D0*DLOG(dum(2)))*DCOS(2.D0*pi*dum(1))
          IF(INHC+1.LE.nhc%NHNC)THEN 
             sigma=DSQRT(KT/nhc%NOSEQ(inhc+1)) !Sqrt(kT/M_i) 
             nhc%VETA(INHC+1)=sigma*DSQRT(-2.D0*DLOG(dum(2)))*DSIN(2.D0*pi*dum(1))
          END IF
       END DO
       CALL rescale_velocities(nhc%nhnc,1,nhc%nhnc,nhc%noseq,T0ions,nhc%veta,dum(1))
    ELSE
       DO IATOMS = 1, NATOMS
          DO INHC = 1, nhc%NHNC
             II = (IATOMS - 1) * nhc%NHNC + INHC
             sigma=DSQRT(KT/nhc%NOSEQ(ii)) !Sqrt(kT/M_i) 
             CALL RANDOM_NUMBER(dum(1))
             CALL RANDOM_NUMBER(dum(2))
             nhc%VETA(II) =sigma*DSQRT(-2.D0*DLOG(dum(2)))*DCOS(2.D0*pi*dum(1))
             IF(II+1.LE.NATOMS*nhc%NHNC)THEN 
                sigma=DSQRT(KT/nhc%NOSEQ(ii+1)) !Sqrt(kT/M_i) 
                nhc%VETA(II+1)=sigma*DSQRT(-2.D0*DLOG(dum(2)))*DSIN(2.D0*pi*dum(1))
             END IF
          END DO
       END DO
       CALL rescale_velocities(natoms*nhc%nhnc,1,natoms*nhc%nhnc,nhc%noseq,T0ions,nhc%veta,dum(1))
    END IF
    !
    !     ** Calculate mv2 **
    !
    AKIN = 0.d0
    DO I = 1, NATOMS
       AKIN = AKIN + amass(I)*             &
            ( vxyz(1,I) * vxyz(1,I) + &
            vxyz(2,I) * vxyz(2,I) + &
            vxyz(3,I) * vxyz(3,I) )
    END DO
    print *, "done nose init state 2 "
    !
    !   ** Calculate the nose accelerations **
    !
    IF(.NOT. nhc%NOSPEAT)THEN
       nhc%AETA(1) = ( AKIN - FKT ) / nhc%NOSEQ(1)
       DO I = 2, nhc%NHNC
          nhc%AETA(I) = ( nhc%NOSEQ(I-1) * nhc%VETA(I-1)**2 - KT) / nhc%NOSEQ(I)
       END DO
    ELSE
       DO IATOMS = 1, NATOMS
          II = (IATOMS - 1) * nhc%NHNC + 1
          AKIN = amass(IATOMS)*                       &
               ( vxyz(1,IATOMS) * vxyz(1,IATOMS) + &
               vxyz(2,IATOMS) * vxyz(2,IATOMS) + &
               vxyz(3,IATOMS) * vxyz(3,IATOMS) )
          nhc%AETA(II) = ( AKIN - FKT) / nhc%NOSEQ(II)
          DO INHC = 2, nhc%NHNC
             II = (IATOMS - 1) * nhc%NHNC + INHC
             nhc%AETA(II) = ( nhc%NOSEQ(II-1) * nhc%VETA(II-1) * nhc%VETA(II-1) &
                  - KT ) / nhc%NOSEQ(II)
          END DO
       END DO
    END IF
    print *, "done nose init state 3 "
    !
    !     ** Initializing factors for Yoshida-Suzuki integration **
    !
    IF(nhc%NSUZUKI > 7 )THEN
       WRITE(*,'(/,A)')' ** MAXIMUM YOSHIKA-SUZUKI ORDER IS 5 **'
       WRITE(*,'(A)')  ' **         NSUZUKI SET TO  7         **'
       nhc%NSUZUKI = 7
    END IF
    IF(MOD(nhc%NSUZUKI,2)==0)THEN
       WRITE(*,'(/,A)')' **    NSUZUKI  CAN ONLY BE 3,5 or 7   **'
       WRITE(*,'(A)')  ' **           NSUZUKI SET TO  7        **'
       nhc%NSUZUKI = 7
    END IF
    IF(nhc%NSUZUKI == 3)THEN
       nhc%DTYS(1) = 1.d0 / ( 2.d0 - 2.d0**CUBROOT )
       nhc%DTYS(2) = 1.d0 - 2.d0 * nhc%DTYS(1)
       nhc%DTYS(3) = nhc%DTYS(1)
       nhc%DTYS(1) = nhc%DTYS(1) * DT  / REAL(NHC%NMULTINT,gp) 
       nhc%DTYS(2) = nhc%DTYS(2) * DT  / REAL(NHC%NMULTINT,gp)
       nhc%DTYS(3) = nhc%DTYS(3) * DT  / REAL(NHC%NMULTINT,gp)
       nhc%DTYS(4) = 0.d0
       nhc%DTYS(5) = 0.d0
       nhc%DTYS(6) = 0.d0
       nhc%DTYS(7) = 0.d0
    ELSE IF(nhc%NSUZUKI == 5)THEN
       nhc%DTYS(1) = 1.d0 / ( 4.d0 - 4.d0**CUBROOT )
       nhc%DTYS(3) = 1.d0 - 4.d0 * nhc%DTYS(1)
       nhc%DTYS(1) = nhc%DTYS(1) * DT  / REAL(NHC%NMULTINT,gp)
       nhc%DTYS(3) = nhc%DTYS(3) * DT  / REAL(NHC%NMULTINT,gp)
       nhc%DTYS(2) = nhc%DTYS(1)
       nhc%DTYS(4) = nhc%DTYS(1)
       nhc%DTYS(5) = nhc%DTYS(1)
       nhc%DTYS(6) = nhc%DTYS(1)
       nhc%DTYS(7) = nhc%DTYS(1)
    ELSE IF(nhc%NSUZUKI == 7)THEN
       WYS1    =  0.784513610477560D0
       WYS2    =  0.235573213359357D0
       WYS3    = -0.117767998417887D1
       WYS4    =  1.d0 - 2.d0 * (WYS1 + WYS2 + WYS3 )
       nhc%DTYS(1) =  WYS3 * DT  / REAL(NHC%NMULTINT,gp)
       nhc%DTYS(2) =  WYS2 * DT  / REAL(NHC%NMULTINT,gp)       
       nhc%DTYS(3) =  WYS1 * DT  / REAL(NHC%NMULTINT,gp)
       nhc%DTYS(4) =  WYS4 * DT  / REAL(NHC%NMULTINT,gp)
       nhc%DTYS(5) =  WYS1 * DT  / REAL(NHC%NMULTINT,gp)
       nhc%DTYS(6) =  WYS2 * DT  / REAL(NHC%NMULTINT,gp)
       nhc%DTYS(7) =  WYS3 * DT  / REAL(NHC%NMULTINT,gp)
    END IF
    !
    print *, "done nose init state 4 "
    RETURN
  END SUBROUTINE nose_init
  !
  SUBROUTINE NOSE_EVOLVE(natoms,ndof,T0ions,amass,vxyz,nhc)
    !use nose_hoover_chains_data
    IMPLICIT NONE 
    !
    !   ------------------------------------------------------------------
    !   Function : Integration of Nose Hoover chain equations of motion
    !   ------------------------------------------------------------------
    !
    !
    INTEGER :: natoms,ndof
    REAL (KIND=8) :: T0ions, amass(*), vxyz(3,*)
    type(NHC_data), intent(inout) :: nhc
    !
    INTEGER         ::  I, J, K, M, II, INHC, IATOMS, III, DOF
    REAL (KIND=8)            ::  WFUN
    REAL (KIND=8)            ::  AKIN, EFR, DTYS2, DTYS4, DTYS8, SCAL, FKT, KT
    !
    REAL (KIND=8), PARAMETER  :: AMU_TO_AU = 1822.888485D0, &
         AU_TO_K   = 315774.664550534774D0, & 
         CMI_TO_AU = 7.26D-7

    !
    IF(.NOT. nhc%NOSPEAT)THEN
       DOF = ndof
    ELSE
       DOF = 3
    END IF
    !
    KT  = T0ions / AU_TO_K
    FKT = KT * REAL(DOF)
    !
    SCAL = 1.d0
    !
    !   ** Calculate kinetic energy **
    !
    AKIN = 0.d0
    DO I = 1, NATOMS
       AKIN = AKIN + amass(I)*        &
            ( vxyz(1,I)**2 +   &
            vxyz(2,I)**2 +   &
            vxyz(3,I)**2 )
    END DO
    !
    !
    DO K = 1, NHC%NMULTINT
       DO J = 1, nhc%NSUZUKI
          !
          !
          DTYS2 = NHC%DTYS(J) / 2.d0
          DTYS4 = NHC%DTYS(J) / 4.d0
          DTYS8 = NHC%DTYS(J) / 8.d0
          !
          IF(.NOT. NHC%NOSPEAT)THEN
             NHC%AETA(NHC%NHNC) = (NHC%NOSEQ(NHC%NHNC-1)*NHC%VETA(NHC%NHNC-1)*NHC%VETA(NHC%NHNC-1) & 
                  - KT)/ NHC%NOSEQ(NHC%NHNC)
             NHC%VETA(NHC%NHNC) = NHC%VETA(NHC%NHNC) +  DTYS4 * NHC%AETA(NHC%NHNC)
          ELSE
             DO IATOMS = 1, NATOMS
                II = IATOMS * NHC%NHNC
                NHC%AETA(II) = (nhc%NOSEQ(II-1)*NHC%VETA(II-1)*NHC%VETA(II-1) &
                     - KT) / nhc%NOSEQ(II)
                NHC%VETA(II) = NHC%VETA(II) + DTYS4 * NHC%AETA(II)
             END DO
          END IF
          !
          IF(.NOT. NHC%NOSPEAT)THEN
             DO II = 1, NHC%NHNC - 2
                INHC = NHC%NHNC - II
                EFR  = DEXP(-1.d0 * DTYS8 * NHC%VETA(INHC+1) )
                NHC%VETA(INHC) = NHC%VETA(INHC) * EFR * EFR &
                     + DTYS4 * NHC%AETA(INHC) * EFR
             END DO
          ELSE
             DO IATOMS = 1, NATOMS
                DO INHC = 1, NHC%NHNC - 2
                   II = IATOMS * NHC%NHNC - INHC
                   EFR = DEXP(-1.d0 * DTYS8 * NHC%VETA(II+1) )
                   NHC%VETA(II) = NHC%VETA(II) * EFR * EFR &
                        +  DTYS4 * NHC%AETA(II) * EFR
                END DO
             END DO
          END IF
          !
          !
          IF(.NOT. NHC%NOSPEAT)THEN
             EFR = DEXP(-1.d0 * DTYS8 * NHC%VETA(2) )
             NHC%VETA(1) = NHC%VETA(1) * EFR * EFR + DTYS4 * NHC%AETA(1) * EFR
          ELSE
             DO IATOMS = 1, NATOMS
                II = (IATOMS - 1) * NHC%NHNC + 1
                EFR = DEXP(-1.d0 * DTYS8 * NHC%VETA(II+1) )
                NHC%VETA(II) = NHC%VETA(II) * EFR * EFR + DTYS4 * NHC%AETA(II) * EFR
             END DO
          END IF
          !
          !        ** Operation on position of thermostat **
          !
          IF(.NOT. NHC%NOSPEAT)THEN
             DO INHC = 1, NHC%NHNC
                nhc%ETA(INHC) = nhc%ETA(INHC) + DTYS2 * NHC%VETA(INHC)
             END DO
          ELSE
             DO IATOMS = 1, NATOMS
                DO INHC = 1, NHC%NHNC
                   II = (IATOMS-1)*NHC%NHNC + INHC
                   nhc%ETA(II) = nhc%ETA(II) + DTYS2 * NHC%VETA(II)
                END DO
             END DO
          END IF
          !
          !        ** Operation on particle velocities **
          !
          IF(.NOT. NHC%NOSPEAT)THEN
             EFR = DEXP(-1.d0 * DTYS2 * NHC%VETA(1) )
             SCAL = SCAL * EFR
             AKIN = AKIN * EFR * EFR
          ELSE
             DO IATOMS = 1, NATOMS
                II  = (IATOMS-1) * NHC%NHNC + 1
                EFR = DEXP(-1.d0 * DTYS2 * NHC%VETA(II) )
                vxyz(1,IATOMS) = vxyz(1,IATOMS) * EFR
                vxyz(2,IATOMS) = vxyz(2,IATOMS) * EFR
                vxyz(3,IATOMS) = vxyz(3,IATOMS) * EFR 
             END DO
          END IF
          !
          !        ** Operations in reverse NHC cylcles **
          !
          !
          IF(.NOT. NHC%NOSPEAT)THEN
             EFR = DEXP(-1.0 * DTYS8 * NHC%VETA(2) )
             NHC%AETA(1) = (AKIN - FKT)/nhc%NOSEQ(1)
             NHC%VETA(1) = NHC%VETA(1) * EFR * EFR + DTYS4 * NHC%AETA(1) * EFR
          ELSE
             DO IATOMS = 1, NATOMS
                II   = (IATOMS-1)*NHC%NHNC + 1
                EFR  = DEXP(-1.d0 * DTYS8 * NHC%VETA(II+1) )
                AKIN = amass(IATOMS)*       &
                     ( vxyz(1,IATOMS)**2 + &
                     vxyz(2,IATOMS)**2 + &
                     vxyz(3,IATOMS)**2 )
                NHC%AETA(II) = (AKIN - FKT)/nhc%NOSEQ(II)
                NHC%VETA(II) = NHC%VETA(II) * EFR * EFR + DTYS4 * NHC%AETA(II) * EFR
             END DO
          END IF
          !
          !        ** Operation on V_z(j) **
          !
          IF(.NOT. NHC%NOSPEAT)THEN
             DO INHC = 2, NHC%NHNC - 1
                EFR  = DEXP(-1.d0 * DTYS8 * NHC%VETA(INHC+1) )
                NHC%AETA(INHC) = (NHC%NOSEQ(INHC-1)*NHC%VETA(INHC-1)*NHC%VETA(INHC-1) &
                     - KT)/NHC%NOSEQ(INHC)
                NHC%VETA(INHC) = NHC%VETA(INHC) * EFR * EFR & 
                     + DTYS4 * NHC%AETA(INHC) * EFR
             END DO
          ELSE 
             DO IATOMS  = 1, NATOMS
                DO INHC = 2, NHC%NHNC - 1
                   II   = (IATOMS - 1)*NHC%NHNC + INHC
                   EFR  = DEXP(-1.d0 * DTYS8 * NHC%VETA(II+1) )
                   NHC%AETA(II) = (NHC%NOSEQ(II-1)*NHC%VETA(II-1)*NHC%VETA(II-1) &
                        - KT)/NHC%NOSEQ(II)
                   NHC%VETA(II) = NHC%VETA(II) * EFR * EFR &
                        + DTYS4 * NHC%AETA(II) * EFR
                END DO
             END DO
          END IF
          !
          !          ** Operation on V_z(M)**
          !
          IF(.NOT. NHC%NOSPEAT)THEN
             NHC%AETA(NHC%NHNC) = (NHC%NOSEQ(NHC%NHNC-1)*NHC%VETA(NHC%NHNC-1)*NHC%VETA(NHC%NHNC-1) & 
                  - KT)/ NHC%NOSEQ(NHC%NHNC)
             NHC%VETA(NHC%NHNC) = NHC%VETA(NHC%NHNC) + DTYS4 * NHC%AETA(NHC%NHNC)
          ELSE
             DO IATOMS = 1, NATOMS
                II = IATOMS * NHC%NHNC
                NHC%AETA(II) = (NHC%NOSEQ(II-1)*NHC%VETA(II-1)*NHC%VETA(II-1) &
                     - KT)/ NHC%NOSEQ(II)
                NHC%VETA(II) = NHC%VETA(II) + DTYS4 * NHC%AETA(II)
             END DO
          END IF
          !
       END DO
    END DO
    !
    !     ** Scale particles with accumulated scaling factor **
    !
    IF(.NOT. NHC%NOSPEAT)THEN
       DO I = 1, NATOMS
          vxyz(1:3,I) = vxyz(1:3,I) * SCAL
       END DO
    END IF
    !
  END SUBROUTINE NOSE_EVOLVE
  !
  !     ***************************************************************
  !     ***************************************************************
  !
  SUBROUTINE NOSE_ENERGY(natoms,ndof,T0ions,nhc)
    !use nose_hoover_chains_data
    IMPLICIT NONE 
    !
    !     ---------------------------------------------------------------
    !     Function : To calculate the energy of extended system variables
    !     ---------------------------------------------------------------
    !
    INTEGER :: natoms, ndof
    REAL (KIND=8) :: T0ions
    type(NHC_data), intent(inout) :: nhc
    !
    REAL (KIND=8), PARAMETER  :: AU_TO_K   = 315774.664550534774D0

    !
    INTEGER         :: INHC, IATOMS, II, DOF
    REAL (KIND=8)            :: KT, FKT
    print *, "Entering nose energy"
    !
    IF(.NOT. NHC%NOSPEAT)THEN
       DOF = ndof
    ELSE
       DOF = 3
    END IF
    KT  = T0ions / AU_TO_K 
    FKT = KT * REAL(DOF)
    !
    nhc%ENOSE = 0.D0
    !
    !     ** Kinetic erergy + potential energy of the thermostat **
    !
    IF(.NOT. NHC%NOSPEAT)THEN
       DO INHC = 1, NHC%NHNC
          NHC%ENOSE = NHC%ENOSE + 0.5D0 * NHC%NOSEQ(INHC) * NHC%VETA(INHC) * NHC%VETA(INHC)
          IF(INHC == 1) THEN
             NHC%ENOSE = NHC%ENOSE + FKT * nhc%ETA(1)
          ELSE 
             NHC%ENOSE = NHC%ENOSE + KT  * nhc%ETA(INHC)
          END IF
       END DO
    ELSE
       DO IATOMS = 1, NATOMS
          II = (IATOMS-1)*NHC%NHNC + 1
          NHC%ENOSE = NHC%ENOSE + 0.5D0 * NHC%NOSEQ(II) * NHC%VETA(II) * NHC%VETA(II)
          NHC%ENOSE = NHC%ENOSE + FKT * nhc%ETA(II)
          DO INHC = 2, NHC%NHNC
             II = (IATOMS-1)*NHC%NHNC + INHC
             NHC%ENOSE = NHC%ENOSE + 0.5D0 * NHC%NOSEQ(II) * NHC%VETA(II) * NHC%VETA(II)
             NHC%ENOSE = NHC%ENOSE + KT  * nhc%ETA(II)
          END DO
       END DO
    END IF
    !
  END SUBROUTINE nose_energy
  !

END MODULE Nose_Hoover_Chains
