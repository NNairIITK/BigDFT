!> @file
!! Part of the pseudo program (pseudopotential generation)
!! @author
!!    Alex Willand, under the supervision of Stefan Goedecker
!!    gpu accelerated routines by Raffael Widmer
!!    parts of this program were based on the fitting program by Matthias Krack
!!    http://cvs.berlios.de/cgi-bin/viewcvs.cgi/cp2k/potentials/goedecker/pseudo/v2.2/
!!
!!    Copyright (C) 2010-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Adapted version from numerical recipies
real(kind=8) FUNCTION ZBRENT &
        (FUNC,ng,ngmx,l,lmx,xp,psi,nocc,noccmx,ispin,nsmx, &
        X1,X2,TOL)
   implicit none
   !Arguments
   integer, intent(in) :: ng,ngmx,l,lmx,nocc,noccmx,ispin,nsmx
   real(kind=8), dimension(0:ng), intent(in) :: xp
   real(kind=8), dimension(0:ngmx,noccmx,lmx,nsmx), intent(in) :: psi
   real(kind=8), intent(in) :: X1,X2,TOL
   real(kind=8), external :: FUNC
   !Local variables
   INTEGER, PARAMETER :: ITMAX=100
   REAL(KIND=8), PARAMETER :: EPS=3.D-20
   real(kind=8) :: FA,FB,FC,A,B,C,D,E,P,Q,R,S,TOL1,XM
   integer :: ITER

!   print*,'entered zbrent'
   A=X1
   B=X2
!   FA=FUNC(A)
!   FB=FUNC(B)
   FA=FUNC(ng,l,xp,psi(0,nocc,l+1,ispin),A)
   FB=FUNC(ng,l,xp,psi(0,nocc,l+1,ispin),B)
   IF(FB*FA.GT.0.d0) then
      print*,FA,FB
      STOP 'Root must be bracketed for ZBRENT.'
   endif
   FC=FB
   DO ITER=1,ITMAX
      IF(FB*FC.GT.0.d0) THEN
        C=A
        FC=FA
        D=B-A
        E=D
      ENDIF
      IF(ABS(FC).LT.ABS(FB)) THEN
        A=B
        B=C
        C=A
        FA=FB
        FB=FC
        FC=FA
      ENDIF
      TOL1=TOL
!      TOL1=2.d0*EPS*ABS(B)+0.5d0*TOL
      XM=.5d0*(C-B)
      IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0.d0)THEN
        ZBRENT=B
        RETURN
      ENDIF
      IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
        S=FB/FA
        IF(A.EQ.C) THEN
          P=2.d0*XM*S
          Q=1.d0-S
        ELSE
          Q=FA/FC
          R=FB/FC
          P=S*(2.d0*XM*Q*(Q-R)-(B-A)*(R-1.d0))
          Q=(Q-1.d0)*(R-1.d0)*(S-1.d0)
        ENDIF
        IF(P.GT.0.d0) Q=-Q
        P=ABS(P)
        IF(2.d0*P .LT. MIN(3.d0*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
          E=D
          D=P/Q
        ELSE
          D=XM
          E=D
        ENDIF
      ELSE
        D=XM
        E=D
      ENDIF
      A=B
      FA=FB
      IF(ABS(D) .GT. TOL1) THEN
        B=B+D
      ELSE
        B=B+SIGN(TOL1,XM)
      ENDIF
      FB=FUNC(ng,l,xp,psi(0,nocc,l+1,ispin),B)
   END DO
   STOP 'ZBRENT exceeding maximum iterations.'

END FUNCTION ZBRENT
