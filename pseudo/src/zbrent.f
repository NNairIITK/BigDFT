c     adapted version from numerical recepies
      real*8 FUNCTION ZBRENT
     :     (FUNC,ng,ngmx,l,lmx,xp,psi,nocc,noccmx,ispin,nsmx,
     :     X1,X2,TOL)
      implicit real*8 (a-h,o-z)
      DIMENSION xp(0:ng),psi(0:ngmx,noccmx,lmx,nsmx)
      real*8 func
      external func

      PARAMETER (ITMAX=100,EPS=3.D-20)

c      print*,'entered zbrent'
      A=X1
      B=X2
c      FA=FUNC(A)
c      FB=FUNC(B)
      FA=func(ng,l,xp,psi(0,nocc,l+1,ispin),A)
      FB=func(ng,l,xp,psi(0,nocc,l+1,ispin),B)
      IF(FB*FA.GT.0.d0) then
         print*,FA,FB
         STOP 'Root must be bracketed for ZBRENT.'
      endif
      FC=FB
      DO 11 ITER=1,ITMAX
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
c        TOL1=2.d0*EPS*ABS(B)+0.5d0*TOL
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
        FB=func(ng,l,xp,psi(0,nocc,l+1,ispin),B)
11    CONTINUE
      STOP 'ZBRENT exceeding maximum iterations.'
      ZBRENT=B
      RETURN
      END

