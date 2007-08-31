       SUBROUTINE SYN_REPEATED_PER(nd1,nd2,nd3,x,NUM_TRANS,N1,N2,N3)
       implicit real(kind=8) (a-h,o-z)
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
       implicit real(kind=8) (a-h,o-z)
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
        implicit real(kind=8) (a-h,o-z)
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
        implicit real(kind=8) (a-h,o-z)
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
        implicit real(kind=8) (a-h,o-z)
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
        implicit real(kind=8) (a-h,o-z)
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
       implicit real(kind=8) (a-h,o-z)
       INTEGER RIGHT
       DIMENSION C(0:RIGHT,NT),CD_1(NT,0:RIGHT)
       ALLOCATABLE MOD_MY(:)
        parameter(m=8)
        real(kind=8) ch(-8:9) ,cg(-8:9)
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

       MOD_LEFT=1-M
       MOD_RIGHT=2*LEN_2-2+M

       ALLOCATE(MOD_MY(MOD_LEFT:MOD_RIGHT))

       DO I=MOD_LEFT,MOD_RIGHT
         MOD_MY(I)=MODULO(I,LENC)
       ENDDO

!       nflop=NT*LEN_2*2*M*4
!       call system_clock(ncount1,ncount_rate,ncount_max)

       DO IT=1,NT-11,12
         DO I=0,LEN_2-1
           I2=2*I

           CI_0 =0.D0
           CI_1 =0.D0
           CI_2 =0.D0
           CI_3 =0.D0
           CI_4 =0.D0
           CI_5 =0.D0
           CI_6 =0.D0
           CI_7 =0.D0
           CI_8 =0.D0
           CI_9 =0.D0
           CI_10=0.D0
           CI_11=0.D0

           DI_0 =0.D0
           DI_1 =0.D0
           DI_2 =0.D0
           DI_3 =0.D0
           DI_4 =0.D0
           DI_5 =0.D0
           DI_6 =0.D0
           DI_7 =0.D0
           DI_8 =0.D0
           DI_9 =0.D0
           DI_10=0.D0
           DI_11=0.D0

           DO J=1-M,M
             JI2=MOD_MY(J+I2)

             CHJ=CH(J)
             CGJ=CG(J) 

             CI_0 =CI_0 +CHJ*C(JI2,IT+0 )
             CI_1 =CI_1 +CHJ*C(JI2,IT+1 )
             CI_2 =CI_2 +CHJ*C(JI2,IT+2 )
             CI_3 =CI_3 +CHJ*C(JI2,IT+3 )

             CI_4 =CI_4 +CHJ*C(JI2,IT+4 )
             CI_5 =CI_5 +CHJ*C(JI2,IT+5 )
             CI_6 =CI_6 +CHJ*C(JI2,IT+6 )
             CI_7 =CI_7 +CHJ*C(JI2,IT+7 )
             CI_8 =CI_8  +CHJ*C(JI2,IT+8  )
             CI_9 =CI_9  +CHJ*C(JI2,IT+9  )
             CI_10=CI_10 +CHJ*C(JI2,IT+10 )
             CI_11=CI_11 +CHJ*C(JI2,IT+11 )

             DI_0 =DI_0 +CGJ*C(JI2,IT+0 )
             DI_1 =DI_1 +CGJ*C(JI2,IT+1 )
             DI_2 =DI_2 +CGJ*C(JI2,IT+2 )
             DI_3 =DI_3 +CGJ*C(JI2,IT+3 )

             DI_4 =DI_4 +CGJ*C(JI2,IT+4 )
             DI_5 =DI_5 +CGJ*C(JI2,IT+5 )
             DI_6 =DI_6 +CGJ*C(JI2,IT+6 )
             DI_7 =DI_7 +CGJ*C(JI2,IT+7 )
             DI_8 =DI_8 +CGJ*C(JI2,IT+8 )
             DI_9 =DI_9 +CGJ*C(JI2,IT+9 )
             DI_10=DI_10+CGJ*C(JI2,IT+10)
             DI_11=DI_11+CGJ*C(JI2,IT+11)

           ENDDO

           CD_1(IT+0,I)=CI_0
           CD_1(IT+1,I)=CI_1
           CD_1(IT+2,I)=CI_2
           CD_1(IT+3,I)=CI_3
           CD_1(IT+4,I)=CI_4
           CD_1(IT+5,I)=CI_5
           CD_1(IT+6,I)=CI_6
           CD_1(IT+7,I)=CI_7
           CD_1(IT+8 ,I)=CI_8 
           CD_1(IT+9 ,I)=CI_9 
           CD_1(IT+10,I)=CI_10
           CD_1(IT+11,I)=CI_11

           IL2=LEN_2+I

           CD_1(IT+0,IL2)=DI_0
           CD_1(IT+1,IL2)=DI_1
           CD_1(IT+2,IL2)=DI_2
           CD_1(IT+3,IL2)=DI_3
           CD_1(IT+4,IL2)=DI_4
           CD_1(IT+5,IL2)=DI_5
           CD_1(IT+6,IL2)=DI_6
           CD_1(IT+7,IL2)=DI_7
           CD_1(IT+8 ,IL2)=DI_8 
           CD_1(IT+9 ,IL2)=DI_9 
           CD_1(IT+10,IL2)=DI_10
           CD_1(IT+11,IL2)=DI_11

         ENDDO
       ENDDO

       IT0=IT

       DO IT=IT0,NT
         DO I=0,LEN_2-1
           I2=2*I
           CI=0.D0
           DI=0.D0
           DO J=1-M,M
             JI2=MOD_MY(J+I2)
             CI=CI+CH(J)*C(JI2,IT)
             DI=DI+CG(J)*C(JI2,IT)
           ENDDO
           CD_1(IT,I)=CI
           CD_1(IT,LEN_2+I)=DI
         ENDDO
       ENDDO

!        call system_clock(ncount2,ncount_rate,ncount_max)
!        tel=dble(ncount2-ncount1)/dble(ncount_rate)
!        write(95,'(a40,1x,e11.4,1x,f10.1,1x,i9)') 'ana_rot_per',tel,1.d-6*nflop/tel,nflop
       DEALLOCATE(MOD_MY) 

       END





      SUBROUTINE SYN_ROT_PER(RIGHT1,NT,CD,C1)
!
!     BACKWARD WAVELET TRANSFORM, SYNTHESIS, PERIODIC
!
      implicit real(kind=8) (a-h,o-z)
      INTEGER RIGHT1
      DIMENSION CD(0:RIGHT1,NT),C1(NT,0:RIGHT1)
      ALLOCATABLE MOD_MY(:)
        parameter(m=8)
        real(kind=8) ch(-8:9) ,cg(-8:9)
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

       MOD_LEFT=-M_2
       MOD_RIGHT=LEN_2-1+M_2

       ALLOCATE(MOD_MY(MOD_LEFT:MOD_RIGHT))

       DO I=MOD_LEFT,MOD_RIGHT
         MOD_MY(I)=MODULO(I,LEN_2)
       ENDDO

!       nflop=NT*LEN_2*(2*M_2+1)*8
!       call system_clock(ncount1,ncount_rate,ncount_max)

        DO IT=1,NT-11,12
          DO I=0,LEN_2-1

            CI2_0 =0.D0
            CI2_1 =0.D0
            CI2_2 =0.D0
            CI2_3 =0.D0
            CI2_4 =0.D0
            CI2_5 =0.D0
            CI2_6 =0.D0
            CI2_7 =0.D0
            CI2_8  =0.D0
            CI2_9  =0.D0
            CI2_10 =0.D0
            CI2_11 =0.D0

            CI21_0=0.D0
            CI21_1=0.D0
            CI21_2=0.D0
            CI21_3=0.D0
            CI21_4=0.D0
            CI21_5=0.D0
            CI21_6=0.D0
            CI21_7=0.D0
            CI21_8 =0.D0
            CI21_9 =0.D0
            CI21_10=0.D0
            CI21_11=0.D0

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
              CI2_8   = CI2_8   + CHJ2*CD(I_J,IT+8 ) + CGJ2*CD(I_J2,IT+8 )
              CI2_9   = CI2_9   + CHJ2*CD(I_J,IT+9 ) + CGJ2*CD(I_J2,IT+9 )
              CI2_10  = CI2_10  + CHJ2*CD(I_J,IT+10) + CGJ2*CD(I_J2,IT+10)
              CI2_11  = CI2_11  + CHJ2*CD(I_J,IT+11) + CGJ2*CD(I_J2,IT+11)

              CI21_0 = CI21_0 + CHJ21*CD(I_J,IT+0) + CGJ21*CD(I_J2,IT+0)
              CI21_1 = CI21_1 + CHJ21*CD(I_J,IT+1) + CGJ21*CD(I_J2,IT+1)
              CI21_2 = CI21_2 + CHJ21*CD(I_J,IT+2) + CGJ21*CD(I_J2,IT+2)
              CI21_3 = CI21_3 + CHJ21*CD(I_J,IT+3) + CGJ21*CD(I_J2,IT+3)
              CI21_4 = CI21_4 + CHJ21*CD(I_J,IT+4) + CGJ21*CD(I_J2,IT+4)
              CI21_5 = CI21_5 + CHJ21*CD(I_J,IT+5) + CGJ21*CD(I_J2,IT+5)
              CI21_6 = CI21_6 + CHJ21*CD(I_J,IT+6) + CGJ21*CD(I_J2,IT+6)
              CI21_7 = CI21_7 + CHJ21*CD(I_J,IT+7) + CGJ21*CD(I_J2,IT+7)
              CI21_8  = CI21_8  + CHJ21*CD(I_J,IT+8 ) + CGJ21*CD(I_J2,IT+8 )
              CI21_9  = CI21_9  + CHJ21*CD(I_J,IT+9 ) + CGJ21*CD(I_J2,IT+9 )
              CI21_10 = CI21_10 + CHJ21*CD(I_J,IT+10) + CGJ21*CD(I_J2,IT+10)
              CI21_11 = CI21_11 + CHJ21*CD(I_J,IT+11) + CGJ21*CD(I_J2,IT+11)

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
            C1(IT+8 ,I2 ) = CI2_8 
            C1(IT+9 ,I2 ) = CI2_9 
            C1(IT+10,I2 ) = CI2_10
            C1(IT+11,I2 ) = CI2_11

            C1(IT+0,I21) = CI21_0 
            C1(IT+1,I21) = CI21_1 
            C1(IT+2,I21) = CI21_2 
            C1(IT+3,I21) = CI21_3 
            C1(IT+4,I21) = CI21_4 
            C1(IT+5,I21) = CI21_5 
            C1(IT+6,I21) = CI21_6 
            C1(IT+7,I21) = CI21_7 
            C1(IT+8 ,I21) = CI21_8  
            C1(IT+9 ,I21) = CI21_9  
            C1(IT+10,I21) = CI21_10 
            C1(IT+11,I21) = CI21_11 

          ENDDO
        ENDDO

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

!        call system_clock(ncount2,ncount_rate,ncount_max)
!        tel=dble(ncount2-ncount1)/dble(ncount_rate)
!        write(95,'(a40,1x,e11.4,1x,f10.1,1x,i9)') 'syn_rot_per',tel,1.d-6*nflop/tel,nflop
       DEALLOCATE(MOD_MY)

      END




