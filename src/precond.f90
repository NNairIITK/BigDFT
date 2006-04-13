

        subroutine preconditionall(iproc,nproc,norb,n1,n2,n3,nbox_c,hgrid,  & 
                   nseg,nvctr,keyg,keyv,hpsi)
! Calls the preconditioner for each orbital treated by the processor
        implicit real*8 (a-h,o-z)
        parameter(eps_mach=1.d-12)
        dimension nseg(0:2*norb),nvctr(0:2*norb)
        dimension keyg(2,nseg(2*norb)),keyv(nseg(2*norb))
        dimension hpsi(nvctr(2*norb)),nbox_c(2,3,norb)

        onem=1.d0-eps_mach
        norb_p=int(onem+dble(norb)/dble(nproc))

     do iorb=iproc*norb_p+1,min((iproc+1)*norb_p,norb)
         mseg_c=nseg(2*iorb-1)-nseg(2*iorb-2) 
         mseg_f=nseg(2*iorb  )-nseg(2*iorb-1)
         iseg_c=nseg(2*iorb-2)+1
         iseg_f=nseg(2*iorb-1)+1
         mvctr_c= nvctr(2*iorb-1)-nvctr(2*iorb-2) 
         mvctr_f=(nvctr(2*iorb  )-nvctr(2*iorb-1))/7
         ipsi_c=nvctr(2*iorb-2)+1
         ipsi_f=nvctr(2*iorb-1)+1
         nl1=nbox_c(1,1,iorb) ; nu1=nbox_c(2,1,iorb)
         nl2=nbox_c(1,2,iorb) ; nu2=nbox_c(2,2,iorb)
         nl3=nbox_c(1,3,iorb) ; nu3=nbox_c(2,3,iorb)

!       call precondition_simple(mvctr_c,mvctr_f,hpsi(ipsi_c),hpsi(ipsi_f))
!The whole box is used instead of subbox
       call precondition(n1,n2,n3,hgrid,mseg_c,mvctr_c,mvctr_f,  & 
            keyg(1,iseg_c),keyv(iseg_c),hpsi(ipsi_c),hpsi(ipsi_f))

     enddo

     return
     end


     subroutine precondition_simple(mvctr_c,mvctr_f,hpsi_c,hpsi_f)
     implicit real*8 (a-h,o-z)
     dimension hpsi_c(mvctr_c),hpsi_f(7,mvctr_f)
! Trivial preconditioner just for testing

      small=5.d-3
      do i=1,mvctr_c
        hpsi_c(i)=small*hpsi_c(i)
      enddo

      do i=1,mvctr_f
      do l=1,7
        hpsi_f(l,i)=small*hpsi_f(l,i)
      enddo
      enddo

     return
     end




        subroutine precondition(n1,n2,n3,hgrid,mseg_c,mvctr_c,mvctr_f,keyg_c,keyv_c,hpsi_c,hpsi_f)
        implicit real*8 (a-h,o-z)
        logical :: iprint = .true.
        dimension keyg_c(2,mseg_c),keyv_c(mseg_c),hpsi_c(mvctr_c),hpsi_f(7,mvctr_f)
        real*8, allocatable, dimension(:,:,:) :: hpsip
!
!       WAVELET AND SCALING FUNCTION SECOND DERIVATIVE FILTERS
        PARAMETER(B2=24.8758460293923314D0,A2=3.55369228991319019D0)

       atomic_length=2.d0
! Number of sweeps in wavelet transformation
!      THE BIGGEST SCALING FUNCTION STEP: atomic_length*FAC_LEN
!      (NOT JUST ATOMIC_LENGTH, BECAUSE SO IT IS BETTER IN PRACTICE) 
       FAC_LEN=2.D0  
       NUM_TRANS=NINT(log(atomic_length*FAC_LEN/hgrid)/log(2.d0))
! Find right leading dimensions for array


        N2_NT=2**NUM_TRANS
!       ND1M*2**NUM_TRANS=ND1>N1   
!       ND1M*N2_NT=ND1  
!       ND1M=N1/N2_NT
!       ND1=ND1M*N2_NT=MODULO(N1,N2_NT)*N2_NT
        
        ND1=CEILING( ((N1+1)*1.D0)/(N2_NT*1.D0) ) *N2_NT-1
        ND2=CEILING( ((N2+1)*1.D0)/(N2_NT*1.D0) ) *N2_NT-1
        ND3=CEILING( ((N3+1)*1.D0)/(N2_NT*1.D0) ) *N2_NT-1
 
        if (iprint) then
        write(*,*) 'NUMBER OF WAVELET TRANSFORMS (sweeps)',NUM_TRANS
        WRITE(*,*)'N1=',N1,'N2=',N2,'N3=',N3
        WRITE(*,*)'ND1=',ND1,'ND2=',ND2,'ND3=',ND3
        iprint=.false.
        endif

        allocate(hpsip(0:nd1,0:nd2,0:nd3))
! find leading dimensions that allow for a wavelet analysis
        call zero((nd1+1)*(nd2+1)*(nd3+1),hpsip)

! coarse part
	do iseg=1,mseg_c
          jj=keyv_c(iseg)
          j0=keyg_c(1,iseg)
          j1=keyg_c(2,iseg)
             ii=j0-1
             i3=ii/((n1+1)*(n2+1))
             ii=ii-i3*(n1+1)*(n2+1)
             i2=ii/(n1+1)
             i0=ii-i2*(n1+1)
             i1=i0+j1-j0
	  do i=i0,i1
            hpsip(i,i2,i3)=hpsi_c(i-i0+jj)
          enddo
         enddo

        FAC_H=1.D0/(HGRID*N2_NT)**2
        H0=    1.5D0*A2*FAC_H;    H1=(A2+B2*.5D0)*FAC_H
        H2=(A2*.5D0+B2)*FAC_H;    H3=    1.5D0*B2*FAC_H
!       H0=5.33

!        EPS=.75d0
        EPS=.5d0

 
        CALL MULTI_FORWARD(ND1,ND2,ND3,HPSIP,NUM_TRANS,NN1,NN2,NN3) 

        NNN1=NN1; NNN2=NN2; NNN3=NN3 

        CALL PRECOND_PROPER(ND1,ND2,ND3,HPSIP,NUM_TRANS,NNN1,NNN2,NNN3,H0,H1,H2,H3,EPS)

        CALL MULTI_BACKWARD(ND1,ND2,ND3,HPSIP,NUM_TRANS,NN1,NN2,NN3)

        F1=1.D0/(H1+EPS);         F2=1.D0/(H2+EPS);         F3=1.D0/(H3+EPS)

        DO I=1,MVCTR_F
          HPSI_F(1,I)=HPSI_F(1,I)*F1
          HPSI_F(2,I)=HPSI_F(2,I)*F1
          HPSI_F(4,I)=HPSI_F(4,I)*F1

          HPSI_F(3,I)=HPSI_F(3,I)*F2
          HPSI_F(5,I)=HPSI_F(5,I)*F2
          HPSI_F(6,I)=HPSI_F(6,I)*F2

          HPSI_F(7,I)=HPSI_F(7,I)*F3
        ENDDO


!	write(90,*) hpsip

! coarse part
        do iseg=1,mseg_c
          jj=keyv_c(iseg)
          j0=keyg_c(1,iseg)
          j1=keyg_c(2,iseg)
             ii=j0-1
             i3=ii/((n1+1)*(n2+1))
             ii=ii-i3*(n1+1)*(n2+1)
             i2=ii/(n1+1)
             i0=ii-i2*(n1+1)
             i1=i0+j1-j0
          do i=i0,i1
            hpsi_c(i-i0+jj)=hpsip(i,i2,i3)
          enddo
        enddo

        deallocate(hpsip)

       end         

       SUBROUTINE PRECOND_PROPER(nd1,nd2,nd3,x,NUM_TRANS,N1,N2,N3,H0,H1,H2,H3,EPS)
       implicit real*8 (a-h,o-z)
       dimension  x(0:nd1,0:nd2,0:nd3)

       DO I_TRANS=1,NUM_TRANS
         N1P=2*(N1+1)-1
         N2P=2*(N2+1)-1
         N3P=2*(N3+1)-1

         IF (N1P.GT.ND1) STOP 'N1 BEYOND BORDERS'
         IF (N2P.GT.ND2) STOP 'N2 BEYOND BORDERS'
         IF (N3P.GT.ND3) STOP 'N3 BEYOND BORDERS'

         N1PP=N1+1
         N2PP=N2+1
         N3PP=N3+1

         F1=1.D0/(H1+EPS);         F2=1.D0/(H2+EPS);         F3=1.D0/(H3+EPS)         

         IF (I_TRANS.EQ.1) THEN 

           F0=1.D0/(H0+EPS)

           DO I3=0,N3
             I3P=I3+N3PP
             DO I2=0,N2
               I2P=I2+N2PP
               DO I1=0,N1
                 I1P=I1+N1PP

                 X(I1,I2,I3)=X(I1,I2,I3)*F0

                 X(I1P,I2,I3)=X(I1P,I2,I3)*F1
                 X(I1,I2P,I3)=X(I1,I2P,I3)*F1
                 X(I1,I2,I3P)=X(I1,I2,I3P)*F1

                 X(I1P,I2P,I3)=X(I1P,I2P,I3)*F2
                 X(I1,I2P,I3P)=X(I1,I2P,I3P)*F2
                 X(I1P,I2,I3P)=X(I1P,I2,I3P)*F2

                 X(I1P,I2P,I3P)=X(I1P,I2P,I3P)*F3

               ENDDO
             ENDDO
           ENDDO

         ELSE

           DO I3=0,N3
             I3P=I3+N3PP
             DO I2=0,N2
               I2P=I2+N2PP
               DO I1=0,N1
                 I1P=I1+N1PP

                 X(I1P,I2,I3)=X(I1P,I2,I3)*F1
                 X(I1,I2P,I3)=X(I1,I2P,I3)*F1
                 X(I1,I2,I3P)=X(I1,I2,I3P)*F1

                 X(I1P,I2P,I3)=X(I1P,I2P,I3)*F2
                 X(I1,I2P,I3P)=X(I1,I2P,I3P)*F2
                 X(I1P,I2,I3P)=X(I1P,I2,I3P)*F2

                 X(I1P,I2P,I3P)=X(I1P,I2P,I3P)*F3

               ENDDO
             ENDDO
           ENDDO

         ENDIF  

         N1=N1P
         N2=N2P
         N3=N3P

         H1=H1*4.D0
         H2=H2*4.D0
         H3=H3*4.D0

       ENDDO

       END

       SUBROUTINE MULTI_BACKWARD(nd1,nd2,nd3,x,NUM_TRANS,N1,N2,N3)
       implicit real*8 (a-h,o-z)
       dimension  x(0:nd1,0:nd2,0:nd3)
       ALLOCATABLE XX(:),YY(:),WW(:)

       IF (NUM_TRANS.GE.1)  THEN

         ALLOCATE(YY((ND1+1)*(ND2+1)*(ND3+1)))
         ALLOCATE(XX((ND1+1)*(ND2+1)*(ND3+1)))

       ENDIF

       IF (NUM_TRANS.GE.2) THEN

         NN1=(ND1+1)/2-1
         NN2=(ND2+1)/2-1
         NN3=(ND3+1)/2-1

         ALLOCATE(WW((NN1+1)*(NN2+1)*(NN3+1)))

         DO I_TRANS=1,NUM_TRANS-1

           N1=2*(N1+1)-1
           N2=2*(N2+1)-1
           N3=2*(N3+1)-1

           IF (N1.GT.ND1) STOP 'N1 BEYOND BORDERS'
           IF (N2.GT.ND2) STOP 'N2 BEYOND BORDERS'
           IF (N3.GT.ND3) STOP 'N3 BEYOND BORDERS'

           I=1
           DO I3=0,N3
           DO I2=0,N2
           DO I1=0,N1
                 XX(I)=X(I1,I2,I3)
                 I=I+1
           ENDDO
           ENDDO
           ENDDO

           CALL BACKWARD_3D(N1,N2,N3,XX,YY,WW)

           I=1
           DO I3=0,N3
           DO I2=0,N2
           DO I1=0,N1
                 X(I1,I2,I3)=YY(I)
                 I=I+1
           ENDDO
           ENDDO
           ENDDO

         ENDDO

         DEALLOCATE(WW)

       ENDIF

       IF (NUM_TRANS.GE.1) THEN

         N1=2*(N1+1)-1
         N2=2*(N2+1)-1
         N3=2*(N3+1)-1

         CALL BACKWARD_3D_SELF(N1,N2,N3,X,XX,YY)
         DEALLOCATE(XX,YY)

       ENDIF

       END



       SUBROUTINE MULTI_FORWARD(nd1,nd2,nd3,x,NUM_TRANS,N1,N2,N3)
       implicit real*8 (a-h,o-z)
       dimension  x(0:nd1,0:nd2,0:nd3)
       ALLOCATABLE XX(:),YY(:),WW(:) 

       N1=ND1
       N2=ND2
       N3=ND3

       IF (NUM_TRANS.GE.1)  THEN

         ALLOCATE(YY((ND1+1)*(ND2+1)*(ND3+1))) 
         ALLOCATE(XX((ND1+1)*(ND2+1)*(ND3+1)))

         CALL FORWARD_3D_SELF(N1,N2,N3,X,YY,XX)

         N1=(N1+1)/2-1
         N2=(N2+1)/2-1
         N3=(N3+1)/2-1

       ENDIF

       IF (NUM_TRANS.GE.2) THEN

         ALLOCATE(WW((N1+1)*(N2+1)*(N3+1)))

         DO I_TRANS=2,NUM_TRANS

           I=1
           DO I3=0,N3
           DO I2=0,N2
           DO I1=0,N1
                 XX(I)=X(I1,I2,I3)
                 I=I+1
           ENDDO
           ENDDO
           ENDDO

           CALL FORWARD_3D(N1,N2,N3,XX,YY,WW)
 
           I=1
           DO I3=0,N3
           DO I2=0,N2
           DO I1=0,N1
                 X(I1,I2,I3)=YY(I)
                 I=I+1
           ENDDO
           ENDDO
           ENDDO

           N1=(N1+1)/2-1
           N2=(N2+1)/2-1
           N3=(N3+1)/2-1

         ENDDO
 
         DEALLOCATE(WW)

       ENDIF

       IF (NUM_TRANS.GE.1) THEN 
         DEALLOCATE(YY,XX)         
       ENDIF 

       END




       subroutine BACKWARD_3D(nd1,nd2,nd3,x,y,ww)
! A periodic synthesis (BACKWARD) wavelet transformation
! The input array x is not overwritten
        implicit real*8 (a-h,o-z)
        dimension  x(0:nd1,0:nd2,0:nd3)
        dimension ww(0:nd1,0:nd2,0:nd3)
        dimension  y(0:nd1,0:nd2,0:nd3)

! i1,i2,i3 -> i2,i3,I1
        nt=(nd2+1)*(nd3+1)
        call  BACKWARD_FAST(nd1,nt,x,y)
! i2,i3,I1 -> i3,I1,I2
        nt=(nd3+1)*(nd1+1)
        call  BACKWARD_FAST(nd2,nt,y,ww)
! i3,I1,I2  -> I1,I2,I3
        nt=(nd1+1)*(nd2+1)
        call  BACKWARD_FAST(nd3,nt,ww,y)

        return
        end

       subroutine BACKWARD_3D_SELF(nd1,nd2,nd3,x,y,ww)
! A periodic synthesis (BACKWARD) wavelet transformation
! The input array x is not overwritten
        implicit real*8 (a-h,o-z)
        dimension  x(0:nd1,0:nd2,0:nd3)
        dimension ww(0:nd1,0:nd2,0:nd3)
        dimension  y(0:nd1,0:nd2,0:nd3)

! i1,i2,i3 -> i2,i3,I1
        nt=(nd2+1)*(nd3+1)
        call  BACKWARD_FAST(nd1,nt,x,y)
! i2,i3,I1 -> i3,I1,I2
        nt=(nd3+1)*(nd1+1)
        call  BACKWARD_FAST(nd2,nt,y,ww)
! i3,I1,I2  -> I1,I2,I3
        nt=(nd1+1)*(nd2+1)
        call  BACKWARD_FAST(nd3,nt,ww,x)

        return
        end


        subroutine FORWARD_3D(nd1,nd2,nd3,y,x,ww)
! An analysis (FORWARD) periodic wavelet transformation
! The input array y is NOT overwritten
        implicit real*8 (a-h,o-z)
        dimension  x(0:nd1,0:nd2,0:nd3)
        dimension ww(0:nd1,0:nd2,0:nd3)
        dimension  y(0:nd1,0:nd2,0:nd3)

! I1,I2,I3 -> I2,I3,i1
        nt=(nd2+1)*(nd3+1)
        call  FORWARD_FAST(nd1,nt,y,x)
! I2,I3,i1 -> I3,i1,i2
        nt=(nd3+1)*(nd1+1)
        call  FORWARD_FAST(nd2,nt,x,ww)
! I3,i1,i2 -> i1,i2,i3
        nt=(nd1+1)*(nd2+1)
        call  FORWARD_FAST(nd3,nt,ww,x)

        return
        end

        subroutine FORWARD_3D_SELF(nd1,nd2,nd3,y,x,ww)
! An analysis (FORWARD) periodic wavelet transformation
! The input array y is NOT overwritten
        implicit real*8 (a-h,o-z)
        dimension  x(0:nd1,0:nd2,0:nd3)
        dimension ww(0:nd1,0:nd2,0:nd3)
        dimension  y(0:nd1,0:nd2,0:nd3)

! I1,I2,I3 -> I2,I3,i1
        nt=(nd2+1)*(nd3+1)
        call  FORWARD_FAST(nd1,nt,y,x)
! I2,I3,i1 -> I3,i1,i2
        nt=(nd3+1)*(nd1+1)
        call  FORWARD_FAST(nd2,nt,x,ww)
! I3,i1,i2 -> i1,i2,i3
        nt=(nd1+1)*(nd2+1)
        call  FORWARD_FAST(nd3,nt,ww,y)

        return
        end


      SUBROUTINE FORWARD_FAST(RIGHT,NT,C,CD_1)
!
!      FORWARD WAVELET TRANSFORM, ANALYSIS, PERIODIC
!
       implicit real*8 (a-h,o-z)
       INTEGER RIGHT
       DIMENSION C(0:RIGHT,NT),CD_1(NT,0:RIGHT)
       ALLOCATABLE MOD_MY(:)
       INCLUDE 'sym_16.inc'

       LENC=RIGHT+1
       LEN_2=LENC/2

       MOD_LEFT=1-M
       MOD_RIGHT=2*LEN_2-2+M

       ALLOCATE(MOD_MY(MOD_LEFT:MOD_RIGHT))

       DO I=MOD_LEFT,MOD_RIGHT
         MOD_MY(I)=MODULO(I,LENC)
       ENDDO

       DO IT=1,NT-7,8
         DO I=0,LEN_2-1
           I2=2*I

           CI_0=0.D0
           CI_1=0.D0
           CI_2=0.D0
           CI_3=0.D0

           CI_4=0.D0
           CI_5=0.D0
           CI_6=0.D0
           CI_7=0.D0

           DI_0=0.D0
           DI_1=0.D0
           DI_2=0.D0
           DI_3=0.D0

           DI_4=0.D0
           DI_5=0.D0
           DI_6=0.D0
           DI_7=0.D0

           DO J=1-M,M
             JI2=MOD_MY(J+I2)

             CHTJ=CHT(J)
             CGTJ=CGT(J) 

             CI_0=CI_0+CHTJ*C(JI2,IT+0)
             CI_1=CI_1+CHTJ*C(JI2,IT+1)
             CI_2=CI_2+CHTJ*C(JI2,IT+2)
             CI_3=CI_3+CHTJ*C(JI2,IT+3)

             CI_4=CI_4+CHTJ*C(JI2,IT+4)
             CI_5=CI_5+CHTJ*C(JI2,IT+5)
             CI_6=CI_6+CHTJ*C(JI2,IT+6)
             CI_7=CI_7+CHTJ*C(JI2,IT+7)

             DI_0=DI_0+CGTJ*C(JI2,IT+0)
             DI_1=DI_1+CGTJ*C(JI2,IT+1)
             DI_2=DI_2+CGTJ*C(JI2,IT+2)
             DI_3=DI_3+CGTJ*C(JI2,IT+3)

             DI_4=DI_4+CGTJ*C(JI2,IT+4)
             DI_5=DI_5+CGTJ*C(JI2,IT+5)
             DI_6=DI_6+CGTJ*C(JI2,IT+6)
             DI_7=DI_7+CGTJ*C(JI2,IT+7)

           ENDDO

           CD_1(IT+0,I)=CI_0
           CD_1(IT+1,I)=CI_1
           CD_1(IT+2,I)=CI_2
           CD_1(IT+3,I)=CI_3
           CD_1(IT+4,I)=CI_4
           CD_1(IT+5,I)=CI_5
           CD_1(IT+6,I)=CI_6
           CD_1(IT+7,I)=CI_7

           IL2=LEN_2+I

           CD_1(IT+0,IL2)=DI_0
           CD_1(IT+1,IL2)=DI_1
           CD_1(IT+2,IL2)=DI_2
           CD_1(IT+3,IL2)=DI_3
           CD_1(IT+4,IL2)=DI_4
           CD_1(IT+5,IL2)=DI_5
           CD_1(IT+6,IL2)=DI_6
           CD_1(IT+7,IL2)=DI_7

         ENDDO
       ENDDO

!       IT0=IT+8
       IT0=IT

       DO IT=IT0,NT
         DO I=0,LEN_2-1
           I2=2*I
           CI=0.D0
           DI=0.D0
           DO J=1-M,M
             JI2=MOD_MY(J+I2)
             CI=CI+CHT(J)*C(JI2,IT)
             DI=DI+CGT(J)*C(JI2,IT)
           ENDDO
           CD_1(IT,I)=CI
           CD_1(IT,LEN_2+I)=DI
         ENDDO
       ENDDO

       DEALLOCATE(MOD_MY) 

       END





      SUBROUTINE BACKWARD_FAST(RIGHT1,NT,CD,C1)
!
!     BACKWARD WAVELET TRANSFORM, SYNTHESIS, PERIODIC
!
      implicit real*8 (a-h,o-z)
      INTEGER RIGHT1
      DIMENSION CD(0:RIGHT1,NT),C1(NT,0:RIGHT1)
      ALLOCATABLE MOD_MY(:)
      include 'sym_16.inc'

       M_2=M/2
       LEN_2=(RIGHT1+1)/2

       MOD_LEFT=-M_2
       MOD_RIGHT=LEN_2-1+M_2

       ALLOCATE(MOD_MY(MOD_LEFT:MOD_RIGHT))

       DO I=MOD_LEFT,MOD_RIGHT
         MOD_MY(I)=MODULO(I,LEN_2)
       ENDDO

        DO IT=1,NT-7,8 
          DO I=0,LEN_2-1

            CI2_0 =0.D0
            CI2_1 =0.D0
            CI2_2 =0.D0
            CI2_3 =0.D0
            CI2_4 =0.D0
            CI2_5 =0.D0
            CI2_6 =0.D0
            CI2_7 =0.D0

            CI21_0=0.D0
            CI21_1=0.D0
            CI21_2=0.D0
            CI21_3=0.D0
            CI21_4=0.D0
            CI21_5=0.D0
            CI21_6=0.D0
            CI21_7=0.D0

            DO J=-M_2,M_2
              I_J=MOD_MY(I-J)

              I_J2=I_J+LEN_2

              J2=2*J
              J21=J2+1

              CHJ2=CH(J2)
              CGJ2=CG(J2)

              CHJ21=CH(J21)
              CGJ21=CG(J21)

              CI2_0  = CI2_0  + CHJ2*CD(I_J,IT+0) + CGJ2*CD(I_J2,IT+0)
              CI2_1  = CI2_1  + CHJ2*CD(I_J,IT+1) + CGJ2*CD(I_J2,IT+1)
              CI2_2  = CI2_2  + CHJ2*CD(I_J,IT+2) + CGJ2*CD(I_J2,IT+2)
              CI2_3  = CI2_3  + CHJ2*CD(I_J,IT+3) + CGJ2*CD(I_J2,IT+3)
              CI2_4  = CI2_4  + CHJ2*CD(I_J,IT+4) + CGJ2*CD(I_J2,IT+4)
              CI2_5  = CI2_5  + CHJ2*CD(I_J,IT+5) + CGJ2*CD(I_J2,IT+5)
              CI2_6  = CI2_6  + CHJ2*CD(I_J,IT+6) + CGJ2*CD(I_J2,IT+6)
              CI2_7  = CI2_7  + CHJ2*CD(I_J,IT+7) + CGJ2*CD(I_J2,IT+7)

              CI21_0 = CI21_0 + CHJ21*CD(I_J,IT+0) + CGJ21*CD(I_J2,IT+0)
              CI21_1 = CI21_1 + CHJ21*CD(I_J,IT+1) + CGJ21*CD(I_J2,IT+1)
              CI21_2 = CI21_2 + CHJ21*CD(I_J,IT+2) + CGJ21*CD(I_J2,IT+2)
              CI21_3 = CI21_3 + CHJ21*CD(I_J,IT+3) + CGJ21*CD(I_J2,IT+3)
              CI21_4 = CI21_4 + CHJ21*CD(I_J,IT+4) + CGJ21*CD(I_J2,IT+4)
              CI21_5 = CI21_5 + CHJ21*CD(I_J,IT+5) + CGJ21*CD(I_J2,IT+5)
              CI21_6 = CI21_6 + CHJ21*CD(I_J,IT+6) + CGJ21*CD(I_J2,IT+6)
              CI21_7 = CI21_7 + CHJ21*CD(I_J,IT+7) + CGJ21*CD(I_J2,IT+7)

            ENDDO

            I2=2*I
            I21=I2+1

            C1(IT+0,I2 ) = CI2_0
            C1(IT+1,I2 ) = CI2_1
            C1(IT+2,I2 ) = CI2_2
            C1(IT+3,I2 ) = CI2_3
            C1(IT+4,I2 ) = CI2_4
            C1(IT+5,I2 ) = CI2_5
            C1(IT+6,I2 ) = CI2_6
            C1(IT+7,I2 ) = CI2_7

            C1(IT+0,I21) = CI21_0 
            C1(IT+1,I21) = CI21_1 
            C1(IT+2,I21) = CI21_2 
            C1(IT+3,I21) = CI21_3 
            C1(IT+4,I21) = CI21_4 
            C1(IT+5,I21) = CI21_5 
            C1(IT+6,I21) = CI21_6 
            C1(IT+7,I21) = CI21_7 

          ENDDO
        ENDDO

!       IT0=IT+8
       IT0=IT

       DO IT=IT0,NT
         DO I=0,LEN_2-1
           CI2 =0.D0
           CI21=0.D0
           DO J=-M_2,M_2
             I_J=MOD_MY(I-J)
             CI2  = CI2  + CH(2*J  )*CD(I_J,IT) + CG(2*J  )*CD(I_J+LEN_2,IT)
             CI21 = CI21 + CH(2*J+1)*CD(I_J,IT) + CG(2*J+1)*CD(I_J+LEN_2,IT)
           ENDDO
           C1(IT,2*I  ) = CI2
           C1(IT,2*I+1) = CI21
         ENDDO
       ENDDO

       DEALLOCATE(MOD_MY)

      END






