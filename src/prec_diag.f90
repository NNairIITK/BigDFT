        subroutine prec_diag(n1,n2,n3,hgrid,nseg_c,nvctr_c,nvctr_f,&
         keyg_c,keyv_c,hpsi_c,hpsi_f,C,scal,A2,B2)
! 
!
        implicit real(kind=8) (a-h,o-z)
        dimension keyg_c(2,nseg_c),keyv_c(nseg_c),hpsi_c(nvctr_c),hpsi_f(7,nvctr_f)
        real(kind=8), allocatable, dimension(:,:,:) :: hpsip
        real(kind=8)::scal(0:3) 
       real(kind=8),parameter::atomic_length=2.d0,FAC_LEN=2.D0

!      Number of sweeps in wavelet transformation
!      THE BIGGEST SCALING FUNCTION STEP: atomic_length*FAC_LEN
!      (NOT JUST ATOMIC_LENGTH, BECAUSE SO IT IS BETTER IN PRACTICE) 
       NUM_TRANS=NINT(log(atomic_length*FAC_LEN/hgrid)/log(2.d0))
       N2_NT=2**NUM_TRANS
       !write(*,'(1x,a)') 'NUMBER OF WAVELET TRANSFORMS (sweeps)',NUM_TRANS

! Find right leading dimensions for array

        
!       ND1+1 IS THE MULTIPLE OF N2_N
!       WHICH IS CLOSEST TO N1+1 FROM ABOVE. 
        ND1=CEILING( ((N1+1)*1.D0)/(N2_NT*1.D0) ) *N2_NT-1
!       THE SAME FOR ND2,ND3.
        ND2=CEILING( ((N2+1)*1.D0)/(N2_NT*1.D0) ) *N2_NT-1
        ND3=CEILING( ((N3+1)*1.D0)/(N2_NT*1.D0) ) *N2_NT-1

        !write(*,'(3(1x,a,i0))')'ND1=',ND1,'ND2=',ND2,'ND3=',ND3

        allocate(hpsip(0:nd1,0:nd2,0:nd3),stat=i_stat)
        call memocc(i_stat,product(shape(hpsip))*kind(hpsip),'hpsip','prec_diag')

        HPSIP=0.D0

! coarse part
        do iseg=1,nseg_c
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

        FAC_H=1.D0/((HGRID*N2_NT)**2)

        H0=    1.5D0*A2*FAC_H;    H1=(A2+B2*.5D0)*FAC_H
        H2=(A2*.5D0+B2)*FAC_H;    H3=    1.5D0*B2*FAC_H

!       FORWARD TRANSFORM THE COARSE SCALING FUNCTIONS NUM_TRANS TIMES
        CALL ANA_REPEATED_PER(ND1,ND2,ND3,HPSIP,NUM_TRANS,NN1,NN2,NN3) 

        NNN1=NN1; NNN2=NN2; NNN3=NN3 

!       DIAGONALLY PRECONDITION THE RESULTING COARSE WAVELETS
        CALL PRECOND_PROPER(ND1,ND2,ND3,HPSIP,NUM_TRANS,NNN1,NNN2,NNN3,H0,H1,H2,H3,C)

        HPSIP=HPSIP/SCAL(0) ! apply (wscal)^(-1)

!       BACKWARD TRANSFORM THE COARSE SCALING FUNCTIONS NUM_TRANS TIMES
        CALL SYN_REPEATED_PER(ND1,ND2,ND3,HPSIP,NUM_TRANS,NN1,NN2,NN3)

!       DIAGONALLY PRECONDITION THE FINE WAVELETS
        DO I=1,NVCTR_F
          HPSI_F(1,I)=HPSI_F(1,I)*scal(1)
          HPSI_F(2,I)=HPSI_F(2,I)*scal(1)
          HPSI_F(4,I)=HPSI_F(4,I)*scal(1)

          HPSI_F(3,I)=HPSI_F(3,I)*scal(2)
          HPSI_F(5,I)=HPSI_F(5,I)*scal(2)
          HPSI_F(6,I)=HPSI_F(6,I)*scal(2)

          HPSI_F(7,I)=HPSI_F(7,I)*scal(3)
        ENDDO

! coarse part
        do iseg=1,nseg_c
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

        i_all=-product(shape(hpsip))*kind(hpsip)
        deallocate(hpsip,stat=i_stat)
        call memocc(i_stat,i_all,'hpsip','prec_diag')

       end         

       SUBROUTINE PRECOND_PROPER(nd1,nd2,nd3,x,NUM_TRANS,N1,N2,N3,H0,H1,H2,H3,EPS)
       implicit real(kind=8) (a-h,o-z)
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

           F1=1.D0/(H1+EPS); F2=1.D0/(H2+EPS);  F3=1.D0/(H3+EPS)       


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

