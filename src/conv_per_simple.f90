       SUBROUTINE SYN_REPEATED_PER(nd1,nd2,nd3,x,NUM_TRANS,N1,N2,N3)
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

           CALL SYNTHESE_PER(N1,N2,N3,XX,YY,WW)

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

         CALL SYNTHESE_PER_SELF(N1,N2,N3,X,XX,YY)
         DEALLOCATE(XX,YY)

       ENDIF

       END



       SUBROUTINE ANA_REPEATED_PER(nd1,nd2,nd3,x,NUM_TRANS,N1,N2,N3)
       implicit real*8 (a-h,o-z)
       dimension  x(0:nd1,0:nd2,0:nd3)
       ALLOCATABLE XX(:),YY(:),WW(:) 

       N1=ND1
       N2=ND2
       N3=ND3

       IF (NUM_TRANS.GE.1)  THEN

         ALLOCATE(YY((ND1+1)*(ND2+1)*(ND3+1))) 
         ALLOCATE(XX((ND1+1)*(ND2+1)*(ND3+1)))

         CALL ANALYSE_PER_SELF(N1,N2,N3,X,YY,XX)

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

           CALL ANALYSE_PER(N1,N2,N3,XX,YY,WW)
 
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




       subroutine SYNTHESE_PER(nd1,nd2,nd3,x,y,ww)
! A periodic synthesis (BACKWARD) wavelet transformation
! The input array x is not overwritten
        implicit real*8 (a-h,o-z)
        dimension  x(0:nd1,0:nd2,0:nd3)
        dimension ww(0:nd1,0:nd2,0:nd3)
        dimension  y(0:nd1,0:nd2,0:nd3)

! i1,i2,i3 -> i2,i3,I1
        nt=(nd2+1)*(nd3+1)
        call  SYN_ROT_PER(nd1,nt,x,y)
! i2,i3,I1 -> i3,I1,I2
        nt=(nd3+1)*(nd1+1)
        call  SYN_ROT_PER(nd2,nt,y,ww)
! i3,I1,I2  -> I1,I2,I3
        nt=(nd1+1)*(nd2+1)
        call  SYN_ROT_PER(nd3,nt,ww,y)

        return
        end

       subroutine SYNTHESE_PER_SELF(nd1,nd2,nd3,x,y,ww)
! A periodic synthesis (BACKWARD) wavelet transformation
! The input array x is not overwritten
        implicit real*8 (a-h,o-z)
        dimension  x(0:nd1,0:nd2,0:nd3)
        dimension ww(0:nd1,0:nd2,0:nd3)
        dimension  y(0:nd1,0:nd2,0:nd3)

! i1,i2,i3 -> i2,i3,I1
        nt=(nd2+1)*(nd3+1)
        call  SYN_ROT_PER(nd1,nt,x,y)
! i2,i3,I1 -> i3,I1,I2
        nt=(nd3+1)*(nd1+1)
        call  SYN_ROT_PER(nd2,nt,y,ww)
! i3,I1,I2  -> I1,I2,I3
        nt=(nd1+1)*(nd2+1)
        call  SYN_ROT_PER(nd3,nt,ww,x)

        return
        end


        subroutine ANALYSE_PER(nd1,nd2,nd3,y,x,ww)
! An analysis (FORWARD) periodic wavelet transformation
! The input array y is NOT overwritten
        implicit real*8 (a-h,o-z)
        dimension  x(0:nd1,0:nd2,0:nd3)
        dimension ww(0:nd1,0:nd2,0:nd3)
        dimension  y(0:nd1,0:nd2,0:nd3)

! I1,I2,I3 -> I2,I3,i1
        nt=(nd2+1)*(nd3+1)
        call  ANA_ROT_PER(nd1,nt,y,x)
! I2,I3,i1 -> I3,i1,i2
        nt=(nd3+1)*(nd1+1)
        call  ANA_ROT_PER(nd2,nt,x,ww)
! I3,i1,i2 -> i1,i2,i3
        nt=(nd1+1)*(nd2+1)
        call  ANA_ROT_PER(nd3,nt,ww,x)

        return
        end

        subroutine ANALYSE_PER_SELF(nd1,nd2,nd3,y,x,ww)
! An analysis (FORWARD) periodic wavelet transformation
! The input array y is NOT overwritten
        implicit real*8 (a-h,o-z)
        dimension  x(0:nd1,0:nd2,0:nd3)
        dimension ww(0:nd1,0:nd2,0:nd3)
        dimension  y(0:nd1,0:nd2,0:nd3)

! I1,I2,I3 -> I2,I3,i1
        nt=(nd2+1)*(nd3+1)
        call  ANA_ROT_PER(nd1,nt,y,x)
! I2,I3,i1 -> I3,i1,i2
        nt=(nd3+1)*(nd1+1)
        call  ANA_ROT_PER(nd2,nt,x,ww)
! I3,i1,i2 -> i1,i2,i3
        nt=(nd1+1)*(nd2+1)
        call  ANA_ROT_PER(nd3,nt,ww,y)

        return
        end


      SUBROUTINE ANA_ROT_PER(RIGHT,NT,C,CD_1)
!
!      FORWARD WAVELET TRANSFORM, ANALYSIS, PERIODIC
!
       implicit real*8 (a-h,o-z)
       INTEGER RIGHT
       DIMENSION C(0:RIGHT,NT),CD_1(NT,0:RIGHT)
        parameter(m=8)
        real*8 ch(-8:9) ,cg(-8:9)
!       Daubechy S16
        data ch  /  0.d0 , -0.0033824159510050025955D0, & 
                -0.00054213233180001068935D0, 0.031695087811525991431D0, & 
                 0.0076074873249766081919D0, -0.14329423835127266284D0, & 
                -0.061273359067811077843D0, 0.48135965125905339159D0,  & 
                 0.77718575169962802862D0,0.36444189483617893676D0, &
                -0.051945838107881800736D0,-0.027219029917103486322D0, &
                 0.049137179673730286787D0,0.0038087520138944894631D0, &
                -0.014952258337062199118D0,-0.00030292051472413308126D0, &
                 0.0018899503327676891843D0 , 0.d0 /
        data cg  / 0.d0 , -0.0018899503327676891843D0, &
                -0.00030292051472413308126D0, 0.014952258337062199118D0, &
                 0.0038087520138944894631D0, -0.049137179673730286787D0, &
                -0.027219029917103486322D0, 0.051945838107881800736D0, &
                 0.36444189483617893676D0, -0.77718575169962802862D0, &
                 0.48135965125905339159D0, 0.061273359067811077843D0, &
                -0.14329423835127266284D0, -0.0076074873249766081919D0, &
                 0.031695087811525991431D0, 0.00054213233180001068935D0, &
                -0.0033824159510050025955D0 , 0.d0 /

       LENC=RIGHT+1
       LEN_2=LENC/2

       DO IT=1,NT
!      *NT       
         DO I=0,LEN_2-1
!        *LEN_2
           I2=2*I
           CI=0.D0
           DI=0.D0
           DO J=1-M,M
!          *2*M (I.E.,16)
             JI2=MODULO(J+I2,LENC)
             CI=CI+CH(J)*C(JI2,IT)
             DI=DI+CG(J)*C(JI2,IT)
!            *4: DO NOT COUNT MODULO             
           ENDDO
           CD_1(IT,I)=CI
           CD_1(IT,LEN_2+I)=DI
         ENDDO
       ENDDO
!      ANA_ROT_PER: NT*LEN_2*2*M*4 FLOPS

       END

      SUBROUTINE SYN_ROT_PER(RIGHT1,NT,CD,C1)
!
!     BACKWARD WAVELET TRANSFORM, SYNTHESIS, PERIODIC
!
      implicit real*8 (a-h,o-z)
      INTEGER RIGHT1
      DIMENSION CD(0:RIGHT1,NT),C1(NT,0:RIGHT1)
        parameter(m=8)
        real*8 ch(-8:9) ,cg(-8:9)
!       Daubechy S16
        data ch  /  0.d0 , -0.0033824159510050025955D0, & 
                -0.00054213233180001068935D0, 0.031695087811525991431D0, & 
                 0.0076074873249766081919D0, -0.14329423835127266284D0, & 
                -0.061273359067811077843D0, 0.48135965125905339159D0,  & 
                 0.77718575169962802862D0,0.36444189483617893676D0, &
                -0.051945838107881800736D0,-0.027219029917103486322D0, &
                 0.049137179673730286787D0,0.0038087520138944894631D0, &
                -0.014952258337062199118D0,-0.00030292051472413308126D0, &
                 0.0018899503327676891843D0 , 0.d0 /
        data cg  / 0.d0 , -0.0018899503327676891843D0, &
                -0.00030292051472413308126D0, 0.014952258337062199118D0, &
                 0.0038087520138944894631D0, -0.049137179673730286787D0, &
                -0.027219029917103486322D0, 0.051945838107881800736D0, &
                 0.36444189483617893676D0, -0.77718575169962802862D0, &
                 0.48135965125905339159D0, 0.061273359067811077843D0, &
                -0.14329423835127266284D0, -0.0076074873249766081919D0, &
                 0.031695087811525991431D0, 0.00054213233180001068935D0, &
                -0.0033824159510050025955D0 , 0.d0 /

        M_2=M/2
        LEN_2=(RIGHT1+1)/2

        DO IT=1,NT
!       *NT
          DO I=0,LEN_2-1
!         *LEN_2
            CI2 =0.D0
            CI21=0.D0
            DO J=-M_2,M_2
!           *(2*M_2+1)
              I_J=MODULO(I-J,LEN_2)
              CI2  = CI2  + CH(2*J  )*CD(I_J,IT) + CG(2*J  )*CD(I_J+LEN_2,IT)
              CI21 = CI21 + CH(2*J+1)*CD(I_J,IT) + CG(2*J+1)*CD(I_J+LEN_2,IT)
!             *8: DO NOT COUNT MODULO
            ENDDO
            C1(IT,2*I  ) = CI2
            C1(IT,2*I+1) = CI21
          ENDDO
        ENDDO
!       SYN_ROT_PER:  NT*LEN_2*(2*M_2+1)*8 FLOPS

      END





