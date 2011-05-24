
C     ==================================================================
      SUBROUTINE XC(RHO,EX,EC,VX,VC)
C     ==--------------------------------------------------------------==
C     ==  LDA EXCHANGE AND CORRELATION FUNCTIONALS                    ==
C     ==                                                              ==
C     ==  EXCHANGE  :  SLATER alpha                                   ==
C     ==  CORRELATION : CEPERLEY & ALDER (PERDEW-ZUNGER PARAMETERS)   ==
C     ==                VOSKO, WILK & NUSSAIR                         ==
C     ==                LEE, YANG & PARR                              ==
C     ==                PERDEW & WANG                                 ==
C     ==                WIGNER                                        ==
C     ==                HEDIN & LUNDQVIST                             ==
C     ==                ORTIZ & BALLONE (PERDEW-ZUNGER FORMULA)       ==
C     ==                ORTIZ & BALLONE (PERDEW-WANG FORMULA)         ==
C     ==                HCTH/120                                      ==
C     ==--------------------------------------------------------------==
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'func.inc'
      PARAMETER (SMALL=1.D-10)
      PARAMETER (PI34= 0.75D0 / 3.141592653589793D+00
     +         ,THIRD=1.D0/3.D0)
C     ==--------------------------------------------------------------==
C..Exchange
      IF(MFXCX.EQ.1) THEN
        CALL SLATERX(RHO,EX,VX,SALPHA)
      ELSE
        EX=0.0D0
        VX=0.0D0
      ENDIF
      IF(RHO.LE.SMALL) THEN
        EC = 0.0D0
        VC = 0.0D0
        EX = 0.0D0
        VX = 0.0D0
      ELSE IF(MFXCC.EQ.1) THEN
        RS=(PI34/RHO)**THIRD
        IFLG=2
        IF(RS.LT.1.0D0) IFLG=1
        CALL PZ(RS,EC,VC,IFLG)
      ELSEIF(MFXCC.EQ.2) THEN
        RS = (PI34/RHO)**THIRD
        CALL VWN(RS,EC,VC)
      ELSEIF(MFXCC.EQ.3) THEN
        CALL LYP(RHO,EC,VC)
      ELSEIF(MFXCC.EQ.4) THEN
        RS=(PI34/RHO)**THIRD
        IFLG=2
        IF(RS.LT.0.5D0) IFLG=1
        IF(RS.GT.100.D0) IFLG=3
        CALL PW(RS,EC,VC,IFLG)
      ELSEIF(MFXCC.EQ.5) THEN
        CALL WIGNER(RHO,EC,VC)
      ELSEIF(MFXCC.EQ.6) THEN
        CALL HEDIN(RHO,EC,VC)
      ELSEIF(MFXCC.EQ.7) THEN
        RS=(PI34/RHO)**THIRD
        IFLG=2
        IF(RS.LT.1.0D0) IFLG=1
        CALL OBPZ(RS,EC,VC,IFLG)
      ELSEIF(MFXCC.EQ.8) THEN
        RS=(PI34/RHO)**THIRD
        IFLG=2
        IF(RS.LT.0.5D0) IFLG=1
        IF(RS.GT.100.D0) IFLG=3
        CALL OBPW(RS,EC,VC,IFLG)
      ELSEIF(MFXCC.EQ.9) THEN
        RS=(PI34/RHO)**THIRD
        CALL PADE(RS,EC,VC)
      ELSE
        EC=0.0D0
        VC=0.0D0
      ENDIF
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE GCXC(RHO,GRHO,SX,SC,V1X,V2X,V1C,V2C)
      use xc_b97, only: eval_b97
C     ==--------------------------------------------------------------==
C     ==  GRADIENT CORRECTIONS FOR EXCHANGE AND CORRELATION           ==
C     ==                                                              ==
C     ==  EXCHANGE  :  BECKE88                                        ==
C     ==               GGAX                                           ==
C     ==               PBEX                                           ==
C     ==               PBESX                                           ==
C     ==               revPBEX                                        ==
C     ==               HCTH/120                                       ==
C     ==               OPTX                                           ==
C     ==  CORRELATION : PERDEW86                                      ==
C     ==                LEE, YANG & PARR                              ==
C     ==                GGAC                                          ==
C     ==                PBEC                                          ==
C     ==                PBESC                                          ==
C     ==                HCTH/120                                      ==
C     ==--------------------------------------------------------------==
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (SMALL=1.D-10)
      INCLUDE 'func.inc'
C     ==--------------------------------------------------------------==
C..Exchange
      IF(RHO.LE.SMALL) THEN
        SX  = 0.0D0
        V1X = 0.0D0
        V2X = 0.0D0
      ELSEIF(MGCX.EQ.1) THEN
        CALL BECKE88(bbeta,RHO,GRHO,SX,V1X,V2X)
      ELSEIF(MGCX.EQ.2) THEN
        CALL GGAX(RHO,GRHO,SX,V1X,V2X)
      ELSEIF(MGCX.EQ.3) THEN
        CALL PBEX(RHO,GRHO,SX,V1X,V2X)
      ELSEIF(MGCX.EQ.4) THEN
        CALL revPBEX(RHO,GRHO,SX,V1X,V2X)
      ELSEIF(MGCX.EQ.5.AND.MGCC.EQ.5) THEN
        CALL HCTH120(RHO,GRHO,SX,V1X,V2X) ! x&c
        SC=0.0D0
        V1C=0.0D0 
        V2C=0.0D0
      ELSEIF(MGCX.EQ.6) THEN
        CALL OPTX(RHO,GRHO,SX,V1X,V2X)
      ELSEIF(MGCX.EQ.7) THEN
        CALL BECKE88(BBETA,RHO,GRHO,SXA,V1XA,V2XA)
        CALL GGAX(RHO,GRHO,SXB,V1XB,V2XB)
        SX=0.722D0*SXA+0.347D0*SXB
        V1X=0.722D0*V1XA+0.347D0*V1XB
        V2X=0.722D0*V2XA+0.347D0*V2XB
      ELSEIF(MGCX.EQ.8) THEN
        CALL BECKE88(BBETA,RHO,GRHO,SXA,V1XA,V2XA)
        CALL GGAX(RHO,GRHO,SXB,V1XB,V2XB)
        SX=0.542D0*SXA+0.167D0*SXB
        V1X=0.542D0*V1XA+0.167D0*V1XB
        V2X=0.542D0*V2XA+0.167D0*V2XB
      ELSEIF(MGCX.EQ.9) THEN
        CALL PBESX(RHO,GRHO,SX,V1X,V2X)
CMK Extra functionals not yet available with CPMD
      ELSEIF(MGCX.EQ.13.AND.MGCC.EQ.13) THEN
        CALL HCTH(93,RHO,SQRT(GRHO),SX,V1X,V2X)
        SC=0.0D0
        V1C=0.0D0
        V2C=0.0D0
      ELSEIF(MGCX.EQ.14.AND.MGCC.EQ.14) THEN
        CALL HCTH(120,RHO,SQRT(GRHO),SX,V1X,V2X)
        SC=0.0D0
        V1C=0.0D0
        V2C=0.0D0
      ELSEIF(MGCX.EQ.15.AND.MGCC.EQ.15) THEN
        CALL HCTH(147,RHO,SQRT(GRHO),SX,V1X,V2X)
        SC=0.0D0
        V1C=0.0D0
        V2C=0.0D0
      ELSEIF(MGCX.EQ.16.AND.MGCC.EQ.16) THEN
        CALL HCTH(407,RHO,SQRT(GRHO),SX,V1X,V2X)
        SC=0.0D0
        V1C=0.0D0
        V2C=0.0D0
      ELSEIF(MGCX.EQ.11) THEN
        CALL eval_b97(1,RHO,GRHO,SX,V1X,V2X)
        SC=0.0D0
        V1C=0.0D0
        V2C=0.0D0
      ELSEIF(MGCX.EQ.12) THEN
        CALL eval_b97(2,RHO,GRHO,SX,V1X,V2X)
        SC=0.0D0
        V1C=0.0D0
        V2C=0.0D0
      ELSE
        SX=0.0D0
        V1X=0.0D0
        V2X=0.0D0
      ENDIF
C..Correlation
      IF(RHO.LE.SMALL) THEN
        SC  = 0.0D0
        V1C = 0.0D0
        V2C = 0.0D0
      ELSEIF(MGCC.EQ.1) THEN
        CALL PERDEW86(RHO,GRHO,SC,V1C,V2C)
      ELSEIF(MGCC.EQ.2) THEN
        CALL GLYP(RHO,GRHO,SC,V1C,V2C)
      ELSEIF(MGCC.EQ.3) THEN
        CALL GGAC(RHO,GRHO,SC,V1C,V2C)
      ELSEIF(MGCC.EQ.4) THEN
        W1=1.D0
        CALL PBEC(RHO,GRHO,W1,SC,V1C,V2C)
      ELSEIF(MGCC.EQ.6) THEN
        W1=0.74D0
        CALL PBEC(RHO,GRHO,W1,SC,V1C,V2C)
      ELSEIF(MGCC.EQ.7) THEN
        CALL PBESC(RHO,GRHO,W1,SC,V1C,V2C)
      ELSE
        SC=0.0D0
        V1C=0.0D0
        V2C=0.0D0
      ENDIF
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE SLATERX(RHO,EX,VX,ALPHA)
C     ==--------------------------------------------------------------==
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (SMALL=1.D-10)
      PARAMETER (F1 = -1.10783814957303361D0)
      PARAMETER (THIRD=1.D0/3.D0,F43=4.D0/3.D0)
C     ==--------------------------------------------------------------==
      IF(RHO.LE.SMALL) THEN
        EX = 0.0D0
        VX = 0.0D0
      ELSE
        RS = RHO**THIRD
        EX = F1*ALPHA*RS
        VX = F43*F1*ALPHA*RS
      ENDIF
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PZ(RS,EPZ,VPZ,IFLG)
C     ==--------------------------------------------------------------==
C     ==  J.P. PERDEW AND ALEX ZUNGER PRB 23, 5048 (1981)             ==
C     ==--------------------------------------------------------------==
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (A=0.0311D0,B=-0.048D0,C=0.0020D0,D=-0.0116D0,
     *           GC=-0.1423D0,B1=1.0529D0,B2=0.3334D0)
C     ==--------------------------------------------------------------==
      IF(IFLG.EQ.1) THEN
C..High density formula
        XLN=LOG(RS)
        EPZ=A*XLN+B+C*RS*XLN+D*RS
        VPZ=A*XLN+(B-A/3.D0)+2.D0/3.D0*C*RS*XLN+
     *              (2.D0*D-C)/3.D0*RS
      ELSEIF(IFLG.EQ.2) THEN
C..Interpolation formula
        RS1=SQRT(RS)
        RS2=RS
        OX=1.D0+B1*RS1+B2*RS2
        DOX=1.D0+7.D0/6.D0*B1*RS1+4.D0/3.D0*B2*RS2
        EPZ=GC/OX
        VPZ=EPZ*DOX/OX
      ENDIF
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE VWN(RS,EVWN,VVWN)
C     ==--------------------------------------------------------------==
C     ==  S.H VOSKO, L.WILK, AND M. NUSAIR,                           ==
C     ==                 CAN. J. PHYS. 58 1200  (1980)                ==
C     ==--------------------------------------------------------------==
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (A=0.0310907D0,B=3.72744D0,C=12.9352D0,X0=-0.10498D0)
      PARAMETER (TWO=2.0D0)
C     ==--------------------------------------------------------------==
      Q  = SQRT(4.D0*C-B*B)
      F1 = TWO*B/Q
      F2 = B*X0/(X0*X0+B*X0+C)
      F3 = TWO*(TWO*X0+B)/Q
      X  = SQRT(RS)
      FX = X*X+B*X+C
      QX = ATAN(Q/(TWO*X+B))
      EVWN=A*(LOG(RS/FX)+F1*QX-F2*(LOG((X-X0)**2/FX)+F3*QX))
      TXPB=TWO*X+B
      TTQQ=TXPB*TXPB+Q*Q
      VVWN=EVWN - X*A/6.D0*(TWO/X-TXPB/FX-4.D0*B/TTQQ-F2*(TWO/(X-X0)
     *          -TXPB/FX-4.D0*(TWO*X0+B)/TTQQ))
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE LYP(RHO,ELYP,VLYP)
C     ==--------------------------------------------------------------==
C     ==  C. LEE, W. YANG, AND R.G. PARR, PRB 37, 785 (1988)          ==
C     ==  THIS IS ONLY THE LDA PART                                   ==
C     ==--------------------------------------------------------------==
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (A=0.04918D0,B=0.132D0,C=0.2533D0,D=0.349D0)
      PARAMETER (CF=2.87123400018819108D0)
C     ==--------------------------------------------------------------==
      RS=RHO**(-1.D0/3.D0)
      ECRS=B*CF*EXP(-C*RS)
      OX=1.D0/(1.D0+D*RS)
      ELYP=-A*OX*(1.D0+ECRS)
      VLYP=ELYP-RS/3.D0*A*OX*(D*OX+ECRS*(D*OX+C))
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PW(RS,EPWC,VPWC,IFLG)
C     ==--------------------------------------------------------------==
C     ==  J.P. PERDEW AND YUE WANG PRB 45, 13244 (1992)               ==
C     ==--------------------------------------------------------------==
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (A=0.031091D0,A1=0.21370D0,B1=7.5957D0,B2=3.5876D0,
     *           B3=1.6382D0,B4=0.49294D0,C0=A,C1=0.046644D0,
     *           C2=0.00664D0,C3=0.01043D0,D0=0.4335D0,D1=1.4408D0)
C     ==--------------------------------------------------------------==
      EPWC=0.0D0
      VPWC=0.0D0
      IF(IFLG.EQ.1) THEN
C..High density formula
        XLN=LOG(RS)
        EPWC=C0*XLN-C1+C2*RS*XLN-C3*RS
        VPWC=C0*XLN-(C1+C0/3.D0)+2.D0/3.D0*C2*RS*XLN-
     *              (2.D0*C3+C2)/3.D0*RS
      ELSEIF(IFLG.EQ.2) THEN
C..Interpolation formula
        RS1=SQRT(RS)
        RS2=RS
        RS3=RS2*RS1
        RS4=RS2*RS2
        OM=2.D0*A*(B1*RS1+B2*RS2+B3*RS3+B4*RS4)
        DOM=2.D0*A*(0.5D0*B1*RS1+B2*RS2+1.5D0*B3*RS3+2.D0*B4*RS4)
        OLOG=LOG(1.D0+1.0D0/OM)
        EPWC=-2.D0*A*(1.D0+A1*RS)*OLOG
        VPWC=-2.D0*A*(1.D0+2.D0/3.D0*A1*RS)*OLOG
     *       -2.D0/3.D0*A*(1.D0+A1*RS)*DOM/(OM*(OM+1.D0))
      ELSEIF(IFLG.EQ.3) THEN
C..Low density formula
        EPWC=-D0/RS+D1/RS**1.5D0
        VPWC=-4.D0/3.D0*D0/RS+1.5D0*D1/RS**1.5D0
      ENDIF
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE WIGNER(RHO,EXC,FXC)
      IMPLICIT REAL*8 (A-H,O-Z)
      RH=RHO
      X=RH**0.33333333333333333D0
      FXC=-X*((0.943656D0+8.8963D0*X)/(1.0D0+12.57D0*X)**2)
      EXC=-0.738D0*X*(0.959D0/(1.0D0+12.57D0*X))
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE HEDIN(RHO,ECP,FCP)
      IMPLICIT REAL*8 (A-H,O-Z)
cmb-ike   VARIABLES with SAVE attribute hinder the vectorization
cmb-ike      SAVE RH
cmb-ike      IF(RH .EQ. 0.0D0) RETURN
      RH=RHO
      RSM1=0.62035049D0*RH**(0.3333333333333333D0)
      ALN=DLOG(1.0D0 + 21.0D0*RSM1)
      X=21.0D0/RSM1
      ECP = ALN+(X**3*ALN-X*X)+X/2.0D0-1.0D0/3.0D0
      ECP = -0.0225D0*ECP
      FCP = -0.0225D0*ALN
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE OBPZ(RS,EPZ,VPZ,IFLG)
C     ==--------------------------------------------------------------==
C     ==  G.ORTIZ AND P. BALLONE PRB 50, 1391 (1994)                  ==
C     ==  PERDEW-ZUNGER FORMULA                                       ==
C     ==--------------------------------------------------------------==
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (A=0.031091D0,B=-0.046644D0,C=0.00419D0,D=-0.00983D0,
     *           GC=-0.103756D0,B1=0.56371D0,B2=0.27358D0)
C     ==--------------------------------------------------------------==
      IF(IFLG.EQ.1) THEN
C..High density formula
        XLN=LOG(RS)
        EPZ=A*XLN+B+C*RS*XLN+D*RS
        VPZ=A*XLN+(B-A/3.D0)+2.D0/3.D0*C*RS*XLN+
     *              (2.D0*D-C)/3.D0*RS
      ELSEIF(IFLG.EQ.2) THEN
C..Interpolation formula
        RS1=SQRT(RS)
        RS2=RS
        OX=1.D0+B1*RS1+B2*RS2
        DOX=1.D0+7.D0/6.D0*B1*RS1+4.D0/3.D0*B2*RS2
        EPZ=GC/OX
        VPZ=EPZ*DOX/OX
      ENDIF
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE OBPW(RS,EPWC,VPWC,IFLG)
C     ==--------------------------------------------------------------==
C     ==  G.ORTIZ AND P. BALLONE PRB 50, 1391 (1994)                  ==
C     ==  PERDEW-WANG FORMULA                                         ==
C     ==--------------------------------------------------------------==
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (A=0.031091D0,A1=0.026481D0,B1=7.5957D0,B2=3.5876D0,
     *           B3=-0.46647D0,B4=0.13354D0,C0=A,C1=0.046644D0,
     *           C2=0.00664D0,C3=0.01043D0,D0=0.4335D0,D1=1.4408D0)
C     ==--------------------------------------------------------------==
      EPWC=0.0D0
      VPWC=0.0D0
      IF(IFLG.EQ.1) THEN
C..High density formula
        XLN=LOG(RS)
        EPWC=C0*XLN-C1+C2*RS*XLN-C3*RS
        VPWC=C0*XLN-(C1+C0/3.D0)+2.D0/3.D0*C2*RS*XLN-
     *              (2.D0*C3+C2)/3.D0*RS
      ELSEIF(IFLG.EQ.2) THEN
C..Interpolation formula
        RS1=SQRT(RS)
        RS2=RS
        RS3=RS2*RS1
        RS4=RS2*RS2
        OM=2.D0*A*(B1*RS1+B2*RS2+B3*RS3+B4*RS4)
        DOM=2.D0*A*(0.5D0*B1*RS1+B2*RS2+1.5D0*B3*RS3+2.D0*B4*RS4)
        OLOG=LOG(1.D0+1.0D0/OM)
        EPWC=-2.D0*A*(1.0D0+A1*RS)*OLOG
        VPWC=-2.D0*A*(1.D0+2.D0/3.D0*A1*RS)*OLOG
     *       -2.D0/3.D0*A*(1.D0+A1*RS)*DOM/(OM*(OM+1.D0))
      ELSEIF(IFLG.EQ.3) THEN
C..Low density formula
        EPWC=-D0/RS+D1/RS**1.5D0
        VPWC=-4.D0/3.D0*D0/RS+1.5D0*D1/RS**1.5D0
      ENDIF
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PADE(RS,EC,VC)
C     ==--------------------------------------------------------------==
C     ==  PADE APPROXIMATION                                          ==
C     ==  S. GOEDECKER, M. TETER, J. HUTTER, PRB in press             ==
C     ==--------------------------------------------------------------==
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (A0=0.4581652932831429D0,A1=2.217058676663745D0,
     *           A2=0.7405551735357053D0,A3=0.01968227878617998D0)
      PARAMETER (B1=1.0000000000000000D0,B2=4.504130959426697D0,
     *           B3=1.110667363742916D0,B4=0.02359291751427506D0)
      PARAMETER (O3=1.D0/3.D0)
C     ==--------------------------------------------------------------==
      TOP=A0+RS*(A1+RS*(A2+RS*A3))
      DTOP=A1+RS*(2.D0*A2+3.D0*A3*RS)
      BOT=RS*(B1+RS*(B2+RS*(B3+RS*B4)))
      DBOT=B1+RS*(2.D0*B2+RS*(3.D0*B3+RS*4.D0*B4))
      EC=-TOP/BOT
      VC=EC+RS*O3*(DTOP/BOT-TOP*DBOT/(BOT*BOT))
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE BECKE88(B1,RHO,GRHO,SX,V1X,V2X)
C     ==--------------------------------------------------------------==
C BECKE EXCHANGE: PRA 38, 3098 (1988)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(OB3=1.D0/3.D0)
C     ==--------------------------------------------------------------==
      TWO13 = 2.0D0**(1.D0/3.D0)
      AA    = GRHO
      A     = SQRT(AA)
      BR1   = RHO**OB3
      BR2   = BR1*BR1
      BR4   = BR2*BR2
      XS    = TWO13*A/BR4
      XS2   = XS*XS
      SA2B8 = SQRT(1.0D0+XS2)
      SHM1  = LOG(XS+SA2B8)
      DD    = 1.0D0 + 6.0D0*B1*XS*SHM1
      DD2   = DD*DD
      EE    = 6.0D0*B1*XS2/SA2B8 - 1.D0
      SX    = TWO13*AA/BR4*(-B1/DD)
      V1X   = -(4.D0/3.D0)/TWO13*XS2*B1*BR1*EE/DD2
      V2X   = TWO13*B1*(EE-DD)/(BR4*DD2)
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE GGAX(RHO,GRHO,SX,V1X,V2X)
C     ==--------------------------------------------------------------==
C J.P.PERDEW ET AL. PRB 46 6671 (1992)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(F1=0.19645D0,F2=7.7956D0,F3=0.2743D0,F4=0.1508D0,
     *          F5=0.004D9,PI=3.141592653589793D0)
C     ==--------------------------------------------------------------==
      FP1   = -3.D0/(16.D0*PI)*(3.D0*PI*PI)**(-1.D0/3.D0)
      FP2   = 0.5D0*(3.D0*PI*PI)**(-1.D0/3.D0)
      AA    = GRHO
      A     = SQRT(AA)
      RR    = RHO**(-4.D0/3.D0)
      S     = FP2*A*RR
      S2    = S*S
      S3    = S2*S
      S4    = S2*S2
      EXPS  = F4*EXP(-100.D0*S2)
      AS    = F3-EXPS-F5*S2
      SA2B8 = SQRT(1.0D0+F2*F2*S2)
      SHM1  = LOG(F2*S+SA2B8)
      BS    = 1.D0+F1*S*SHM1+F5*S4
      DAS   = 200.D0*S*EXPS-2.D0*S*F5
      DBS   = F1*(SHM1+F2*S/SA2B8)+4.D0*F5*S3
      DLS   = (DAS/AS-DBS/BS)
      SX    = FP1*AA*RR*AS/BS
      V1X   = -4.D0/3.D0*SX/RHO*(1.D0+S*DLS)
      V2X   = FP1*RR*AS/BS*(2.D0+S*DLS)
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PBESX(RHO,GRHO,SX,V1X,V2X)
C     ==--------------------------------------------------------------==
C PBESol functional
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(US=0.161620459673995492D0,AX=-0.738558766382022406D0,
     *          UM=0.123456790123456789D0,UK=0.8040D0,UL=UM/UK)
C     ==--------------------------------------------------------------==
      AA    = GRHO
      RR    = RHO**(-4.D0/3.D0)
      EX    = AX/RR
      S2    = AA*RR*RR*US*US
      PO    = 1.D0/(1.D0 + UL*S2)
      FX    = UK-UK*PO
      SX    = EX*FX
      DFX   = 2.D0*UK*UL*PO*PO
      V1X   = 1.33333333333333D0*AX*RHO**0.333333333333D0*(FX-S2*DFX)
      V2X   = EX*DFX*(US*RR)**2
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PBEX(RHO,GRHO,SX,V1X,V2X)
C     ==--------------------------------------------------------------==
C J.P.PERDEW ET AL. PRL 77 3865 (1996)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(US=0.161620459673995492D0,AX=-0.738558766382022406D0,
     *          UM=0.2195149727645171D0,UK=0.8040D0,UL=UM/UK)
C     ==--------------------------------------------------------------==
      AA    = GRHO
      RR    = RHO**(-4.D0/3.D0)
      EX    = AX/RR
      S2    = AA*RR*RR*US*US
      PO    = 1.D0/(1.D0 + UL*S2)
      FX    = UK-UK*PO
      SX    = EX*FX
      DFX   = 2.D0*UK*UL*PO*PO
      V1X   = 1.33333333333333D0*AX*RHO**0.333333333333D0*(FX-S2*DFX)
      V2X   = EX*DFX*(US*RR)**2
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE revPBEX(RHO,GRHO,SX,V1X,V2X)
C     ==--------------------------------------------------------------==
C Y. ZHANG ET AL. PRL 80 890 (1998)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(US=0.161620459673995492D0,AX=-0.738558766382022406D0,
     *          UM=0.2195149727645171D0,UK=1.2450D0,UL=UM/UK)
C     ==--------------------------------------------------------------==
      AA    = GRHO
      RR    = RHO**(-4.D0/3.D0)
      EX    = AX/RR
      S2    = AA*RR*RR*US*US
      PO    = 1.D0/(1.D0 + UL*S2)
      FX    = UK-UK*PO
      SX    = EX*FX
      DFX   = 2.D0*UK*UL*PO*PO
      V1X   = 1.33333333333333D0*AX*RHO**0.333333333333D0*(FX-S2*DFX)
      V2X   = EX*DFX*(US*RR)**2
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PERDEW86(RHO,GRHO,SC,V1C,V2C)
C     ==--------------------------------------------------------------==
C PERDEW CORRELATION: PRB 33, 8822 (1986)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(P1=0.023266D0,P2=7.389D-6,P3=8.723D0,P4=0.472D0)
      PARAMETER(PC1=0.001667D0,PC2=0.002568D0,PCI=PC1+PC2)
      PARAMETER(OB3=1.D0/3.D0, FPI=4.0D0*3.141592653589793D0)
C     ==--------------------------------------------------------------==
      AA    = GRHO
      A     = SQRT(AA)
      BR1   = RHO**OB3
      BR2   = BR1*BR1
      BR4   = BR2*BR2
      RS    = (3.D0/(FPI*RHO))**OB3
      RS2   = RS*RS
      RS3   = RS*RS2
      CNA   = PC2+P1*RS+P2*RS2
      CNB   = 1.D0+P3*RS+P4*RS2+1.D4*P2*RS3
      CN    = PC1 + CNA/CNB
      DRS   = -OB3*(3.D0/FPI)**OB3 / BR4
      DCNA  = (P1+2.D0*P2*RS)*DRS
      DCNB  = (P3+2.D0*P4*RS+3.D4*P2*RS2)*DRS
      DCN   = DCNA/CNB - CNA/(CNB*CNB)*DCNB
      PHI   = 0.192D0*PCI/CN*A*RHO**(-7.D0/6.D0)
      EPHI  = EXP(-PHI)
      SC    = AA/BR4*CN*EPHI
      V1C   = SC*((1.D0+PHI)*DCN/CN -((4.D0/3.D0)-(7.D0/6.D0)*PHI)/RHO)
      V2C   = CN*EPHI/BR4*(2.D0-PHI)
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE GLYP(RHO,GRHO,SC,V1C,V2C)
C     ==--------------------------------------------------------------==
C LEE, YANG PARR: GRADIENT CORRECTION PART
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(A=0.04918D0,B=0.132D0,C=0.2533D0,D=0.349D0)
C     ==--------------------------------------------------------------==
      AA    = GRHO
      R     = RHO**(-1.d0/3.d0)
      OM    = EXP(-C*R)/(1.d0+D*R)
      R5    = R**5
      XL    = 1.d0+(7.d0/3.d0)*(C*R + D*R/(1.d0+D*R))
      FF    = A*B*AA/24.d0
      SC    = FF*R5*OM*XL
      DR5   = 5.d0*R*R*R*R
      DOM   = -OM*(C+D+C*D*R)/(1.d0+D*R)
      DXL   = (7.d0/3.d0)*(C+D+2.d0*C*D*R+C*D*D*R*R)/(1.d0+D*R)**2
      V1C   = -FF*(R*R*R*R)/3.d0*( DR5*OM*XL + R5*DOM*XL + R5*OM*DXL)
      V2C   = A*B*R5*OM*XL/12.d0
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE GGAC(RHO,GRHO,SC,V1C,V2C)
C     ==--------------------------------------------------------------==
C PERDEW & WANG GGA CORRELATION PART      
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(AL=0.09D0,PA=0.023266D0,PB=7.389D-6,PC=8.723D0,
     *          PD=0.472D0,CX=-0.001667D0,CXC0=0.002568D0,CC0=-CX+CXC0)
      PARAMETER(OB3=1.D0/3.D0, PI=3.141592653589793D0)
C     ==--------------------------------------------------------------==
      XNU   = 16.D0/PI*(3.D0*PI*PI)**OB3
      BE    = XNU*CC0
      CALL XC(RHO,EX,EC,VX,VC)
      AA    = GRHO
      A     = SQRT(AA)
      RS    = (3.D0/(4.D0*PI*RHO))**OB3
      RS2   = RS*RS
      RS3   = RS*RS2
      XKF   = (9.D0*PI/4.D0)**OB3/RS
      XKS   = SQRT(4.D0*XKF/PI)
      T     = A/(2.D0*XKS*RHO)
      EXPE  = EXP(-2.D0*AL*EC/(BE*BE))
      AF    = 2.D0*AL/BE * (1.D0/(EXPE-1.D0))
      BF    = EXPE*(VC-EC)
      Y     = AF*T*T
      XY    = (1.D0+Y)/(1.D0+Y+Y*Y)
      QY    = Y*Y*(2.D0+Y)/(1.D0+Y+Y*Y)**2
      S1    = 1.D0+2.D0*AL/BE*T*T*XY
      H0    = BE*BE/(2.D0*AL) * LOG(S1)
      DH0   = BE*T*T/S1*(-7.D0/3.D0*XY-QY*(AF*BF/BE-7.D0/3.D0))
      DDH0  = BE/(2.D0*XKS*XKS*RHO)*(XY-QY)/S1
      EE    = -100.D0*(XKS/XKF*T)**2
      CNA   = CXC0+PA*RS+PB*RS2
      DCNA  = -(PA*RS+2.D0*PB*RS2)/3.D0
      CNB   = 1.D0+PC*RS+PD*RS2+1.D4*PB*RS3
      DCNB  = -(PC*RS+2.D0*PD*RS2+3.D4*PB*RS3)/3.D0
      CN    = CNA/CNB - CX
      DCN   = DCNA/CNB - CNA*DCNB/(CNB*CNB)
      H1    = XNU*(CN-CC0-3.D0/7.D0*CX)*T*T*EXP(EE)
      DH1   = -OB3*(H1*(7.D0+8.D0*EE)+XNU*T*T*EXP(EE)*DCN)
      DDH1  = 2.D0*H1*(1.D0+EE)*RHO/AA
      SC    = RHO*(H0+H1)
      V1C   = H0+H1+DH0+DH1
      V2C   = DDH0+DDH1
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PBESC(RHO,GRHO,W1,SC,V1C,V2C)
C     ==--------------------------------------------------------------==
C PBESol Correlation functional
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(BE=0.046D0,GA=0.031090690869654895D0)
      PARAMETER(OB3=1.D0/3.D0,PI=3.141592653589793D0)
C     ==--------------------------------------------------------------==
      CALL XC(RHO,EX,EC,VX,VC)
      AA    = GRHO
      A     = SQRT(AA)
      RS    = (3.D0/(4.D0*PI*RHO))**OB3
      XKF   = (9.D0*PI/4.D0)**OB3/RS
      XKS   = SQRT(4.D0*XKF/PI)
      T     = A/(2.D0*XKS*RHO)
      EXPE  = EXP(-EC/GA)
      AF    = BE/GA * (1.D0/(EXPE-1.D0))
      Y     = AF*T*T
      XY    = (1.D0+Y)/(1.D0+Y+Y*Y)
      S1    = 1.D0+BE/GA*T*T*XY
      H0    = GA * LOG(S1)
      DTDR  = -T*7.D0/(6.D0*RHO)
      DADR  = AF*AF*EXPE/BE*(VC-EC)/RHO
      DSDA  = -BE/GA * AF * T**6 * (2.D0+Y) / (1.D0+Y+Y*Y)**2
      DSDT  = 2.D0*BE/GA * T * (1.D0+2.D0*Y) / (1.D0+Y+Y*Y)**2
      DSDR  = DSDA*DADR + DSDT*DTDR
      DHDT  = GA/S1*DSDT
      DHDR  = GA/S1*DSDR
      SC    = W1*RHO*H0
      V1C   = W1*H0+W1*DHDR*RHO
      V2C   = W1*RHO*DHDT*T/AA
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE PBEC(RHO,GRHO,W1,SC,V1C,V2C)
C     ==--------------------------------------------------------------==
C PBE Correlation functional
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(BE=0.06672455060314922D0,GA=0.031090690869654895D0)
      PARAMETER(OB3=1.D0/3.D0,PI=3.141592653589793D0)
C     ==--------------------------------------------------------------==
      CALL XC(RHO,EX,EC,VX,VC)
      AA    = GRHO
      A     = SQRT(AA)
      RS    = (3.D0/(4.D0*PI*RHO))**OB3
      XKF   = (9.D0*PI/4.D0)**OB3/RS
      XKS   = SQRT(4.D0*XKF/PI)
      T     = A/(2.D0*XKS*RHO)
      EXPE  = EXP(-EC/GA)
      AF    = BE/GA * (1.D0/(EXPE-1.D0))
      Y     = AF*T*T
      XY    = (1.D0+Y)/(1.D0+Y+Y*Y)
      S1    = 1.D0+BE/GA*T*T*XY
      H0    = GA * LOG(S1)
      DTDR  = -T*7.D0/(6.D0*RHO)
      DADR  = AF*AF*EXPE/BE*(VC-EC)/RHO
      DSDA  = -BE/GA * AF * T**6 * (2.D0+Y) / (1.D0+Y+Y*Y)**2
      DSDT  = 2.D0*BE/GA * T * (1.D0+2.D0*Y) / (1.D0+Y+Y*Y)**2
      DSDR  = DSDA*DADR + DSDT*DTDR
      DHDT  = GA/S1*DSDT
      DHDR  = GA/S1*DSDR
      SC    = W1*RHO*H0
      V1C   = W1*H0+W1*DHDR*RHO
      V2C   = W1*RHO*DHDT*T/AA
C     ==--------------------------------------------------------------==
      RETURN
      END
C     ==================================================================
      SUBROUTINE hcth120(rho,grho,sx,v1x,v2x)
C     HCTH, JCP 109, 6264 (1998)
C     Parameters set-up after N.L. Doltsisnis & M. Sprik (1999)
C     Present release: Tsukuba, 09/02/2005
c--------------------------------------------------------------------------
c     rhoa = rhob = 0.5 * rho
c     grho is the SQUARE of the gradient of rho! --> gr=sqrt(grho)
c     sx  : total exchange correlation energy at point r 
c     v1x : d(sx)/drho  (eq. dfdra = dfdrb in original)
c     v2x : 1/gr*d(sx)/d(gr) (eq. 0.5 * dfdza = 0.5 * dfdzb in original)
c--------------------------------------------------------------------------
      IMPLICIT REAL*8 (a-h,o-z)
      PARAMETER(o3=1.0d0/3.0d0,fr83=8.d0/3.d0)
      DIMENSION cg0(6),cg1(6),caa(6),cab(6),cx(6)
      r3q2=DEXP(-o3*0.69314718055994531d0)
      r3pi=DEXP(-o3*0.04611759718129048d0)
c.....coefficients for PW correlation......................................
      cg0(1)= 0.031091d0
      cg0(2)= 0.213700d0
      cg0(3)= 7.595700d0
      cg0(4)= 3.587600d0
      cg0(5)= 1.638200d0
      cg0(6)= 0.492940d0
      cg1(1)= 0.015545d0
      cg1(2)= 0.205480d0
      cg1(3)=14.118900d0
      cg1(4)= 6.197700d0
      cg1(5)= 3.366200d0
      cg1(6)= 0.625170d0
C......HCTH-19-4.....................................
      caa(1)=  0.489508D+00
      caa(2)= -0.260699D+00
      caa(3)=  0.432917D+00
      caa(4)= -0.199247D+01
      caa(5)=  0.248531D+01
      caa(6)=  0.200000D+00
      cab(1)=  0.514730D+00
      cab(2)=  0.692982D+01
      cab(3)= -0.247073D+02
      cab(4)=  0.231098D+02
      cab(5)= -0.113234D+02
      cab(6)=  0.006000D+00
      cx(1) =  0.109163D+01
      cx(2) = -0.747215D+00
      cx(3) =  0.507833D+01
      cx(4) = -0.410746D+01
      cx(5) =  0.117173D+01
      cx(6) =  0.004000D+00
c...........................................................................
      gr=DSQRT(grho)
      rho_o3=rho**(o3) 
      rho_o34=rho*rho_o3
      xa=1.25992105d0*gr/rho_o34
      xa2=xa*xa
      ra=0.781592642d0/rho_o3
      rab=r3q2*ra
      dra_drho=-0.260530881d0/rho_o34
      drab_drho=r3q2*dra_drho
      CALL pwcorr(ra,cg1,g,dg)
      era1=g
      dera1_dra=dg
      CALL pwcorr(rab,cg0,g,dg)
      erab0=g
      derab0_drab=dg
      ex=-0.75d0*r3pi*rho_o34
      dex_drho=-r3pi*rho_o3
      uaa=caa(6)*xa2
      uaa=uaa/(1.0d0+uaa)
      uab=cab(6)*xa2
      uab=uab/(1.0d0+uab)
      ux=cx(6)*xa2
      ux=ux/(1.0d0+ux)
      ffaa=rho*era1
      ffab=rho*erab0-ffaa
      dffaa_drho=era1+rho*dera1_dra*dra_drho
      dffab_drho=erab0+rho*derab0_drab*drab_drho-dffaa_drho
cmb-> i-loop removed
      denaa=1.d0/(1.0d0+caa(6)*xa2)
      denab=1.d0/(1.0d0+cab(6)*xa2)
      denx =1.d0/(1.0d0+cx(6)*xa2)
      f83rho=fr83/rho
      bygr=2.0d0/gr
      gaa=caa(1)+uaa*(caa(2)+uaa*(caa(3)+uaa*(caa(4)+uaa*caa(5))))
      gab=cab(1)+uab*(cab(2)+uab*(cab(3)+uab*(cab(4)+uab*cab(5))))
      gx=cx(1)+ux*(cx(2)+ux*(cx(3)+ux*(cx(4)+ux*cx(5))))
      taa=denaa*uaa*(caa(2)+uaa*(2.d0*caa(3)+uaa 
     &    *(3.d0*caa(4)+uaa*4.d0*caa(5))))
      tab=denab*uab*(cab(2)+uab*(2.d0*cab(3)+uab
     &    *(3.d0*cab(4)+uab*4.d0*cab(5))))
      txx=denx*ux*(cx(2)+ux*(2.d0*cx(3)+ux
     &    *(3.d0*cx(4)+ux*4.d0*cx(5))))
      dgaa_drho=-f83rho*taa
      dgab_drho=-f83rho*tab
      dgx_drho=-f83rho*txx
      dgaa_dgr=bygr*taa
      dgab_dgr=bygr*tab
      dgx_dgr=bygr*txx
cmb
      sx=ex*gx+ffaa*gaa+ffab*gab
      v1x=dex_drho*gx+ex*dgx_drho
     .   +dffaa_drho*gaa+ffaa*dgaa_drho
     .   +dffab_drho*gab+ffab*dgab_drho
      v2x=(ex*dgx_dgr+ffaa*dgaa_dgr+ffab*dgab_dgr)/gr
      RETURN
      END
C =-------------------------------------------------------------------=
      SUBROUTINE pwcorr(r,c,g,dg)
      IMPLICIT real*8 (a-h,o-z)
      DIMENSION c(6)
      r12=DSQRT(r)
      r32=r*r12
      r2=r*r
      rb=c(3)*r12+c(4)*r+c(5)*r32+c(6)*r2
      sb=1.0d0+1.0d0/(2.0d0*c(1)*rb)
      g=-2.0d0*c(1)*(1.0d0+c(2)*r)*DLOG(sb)
      drb=c(3)/(2.0d0*r12)+c(4)+1.5d0*c(5)*r12+2.0d0*c(6)*r
      dg=(1.0d0+c(2)*r)*drb/(rb*rb*sb)-2.0d0*c(1)*c(2)*DLOG(sb)
      RETURN
      END 
C     ==================================================================
      SUBROUTINE OPTX(rho,grho,sx,v1x,v2x)
C     OPTX, Handy et al. JCP 116, p. 5411 (2002) and refs. therein
C     Present release: Tsukuba, 20/6/2002
c--------------------------------------------------------------------------
c     rhoa = rhob = 0.5 * rho in LDA implementation
c     grho is the SQUARE of the gradient of rho! --> gr=sqrt(grho)
c     sx  : total exchange correlation energy at point r
c     v1x : d(sx)/drho
c     v2x : 1/gr*d(sx)/d(gr)
c--------------------------------------------------------------------------
      IMPLICIT REAL*8 (a-h,o-z)
      PARAMETER(SMALL=1.D-20,SMAL2=1.D-08)
C.......coefficients and exponents....................
      PARAMETER(o43=4.0d0/3.0d0,two13=1.259921049894873D0
     .         ,two53=3.174802103936399D0,gam=0.006D0
     .         ,a1cx=0.9784571170284421D0,a2=1.43169D0)
C.......OPTX in compact form..........................
      IF(RHO.LE.SMALL) THEN
       sx=0.0D0
       v1x=0.0D0
       v2x=0.0D0
      ELSE
       gr=DMAX1(grho,SMAL2)
       rho43=rho**o43
       xa=two13*DSQRT(gr)/rho43
       gamx2=gam*xa*xa
       uden=1.d+00/(1.d+00+gamx2)
       uu=a2*gamx2*gamx2*uden*uden
       uden=rho43*uu*uden
       sx=-rho43*(a1cx+uu)/two13
       v1x=o43*(sx+two53*uden)/rho
       v2x=-two53*uden/gr
      ENDIF
C
      RETURN
      END
C =-------------------------------------------------------------------=


