!> @file
!!   Old Preconditioner routines (unused) 
!! @deprecated
!! @author
!!    Copyright (C) 2010 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Calls the preconditioner for each orbital treated by the processor
subroutine preconditionall(iproc,nproc,norb,norbp,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, &
                   hgrid,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,cprec,logrid_c,logrid_f,hpsi)
  implicit real(kind=8) (a-h,o-z)
  dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
  dimension hpsi(nvctr_c+7*nvctr_f,norbp)

  call cpu_time(tr0)
  call system_clock(ncount1,ncount_rate,ncount_max)

  do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
     cr=-cprec
     call precondition(cr,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,hgrid,&
                   nseg_c,nvctr_c,keyg(1,1),keyv(1),  & 
                   nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),&
                   hpsi(1,iorb-iproc*norbp),hpsi(nvctr_c+1,iorb-iproc*norbp),IORB)

  enddo

  call cpu_time(tr1)
  call system_clock(ncount2,ncount_rate,ncount_max)
  tel=dble(ncount2-ncount1)/dble(ncount_rate)
  write(77,'(a40,i4,2(1x,e10.3))') 'PRECONDITIONING TIME',iproc,tr1-tr0,tel

END SUBROUTINE preconditionall


subroutine precondition(C,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
              hgrid,nseg_c,nvctr_c,keyg_c,&
              keyv_c,nseg_f,nvctr_f,keyg_f,&
              keyv_f,psi_c,psi_f,IORB)
  implicit real(kind=8) (a-h,o-z)
  dimension keyg_c(2,nseg_c),keyv_c(nseg_c),keyg_f(2,nseg_f),keyv_f(nseg_f)
  dimension psi_c(nvctr_c),psi_f(7,nvctr_f)
  
  real(kind=8),allocatable,dimension(:) ::   GG_C, HPSI_C, HPSI_C_OLD
  real(kind=8),allocatable,dimension(:) ::    HPSI_C_MOD, PSI_C_NEW
  real(kind=8),allocatable,dimension(:) ::    AUX,ALPHAS_PR
  real(kind=8),allocatable,dimension(:,:) :: GG_F, HPSI_F, HPSI_F_OLD
  real(kind=8),allocatable,dimension(:,:) ::  HPSI_F_MOD,   PSI_F_NEW
  
  real(kind=8), parameter :: A2=3.55369228991319019D0,ALPHAMIN=.05d0,COSMIN=.5D0

  OPEN(10,FILE='input.dat') 
    DO I=1,8
      READ(10,*)
    ENDDO
    READ(10,*) NSTEP_MIN,NSTEP_MAX
    READ(10,*) 
    READ(10,*) FAC_NORM2          
  CLOSE(10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       DETERMINE THE NUMBER OF SWEEPS ETC
!
  NUM_TRANS=NINT(  LOG( 1.5D0*A2/(HGRID**2*C) )/(2.D0*LOG(2.D0))  )
  write(*,*) 'NUMBER OF WAVELET TRANSFORMS (sweeps)',NUM_TRANS,hgrid,c
! Find right leading dimensions for array

  N2_NT=2**NUM_TRANS
  
  ND1=CEILING( ((N1+1)*1.D0)/(N2_NT*1.D0) ) *N2_NT-1
  ND2=CEILING( ((N2+1)*1.D0)/(N2_NT*1.D0) ) *N2_NT-1
  ND3=CEILING( ((N3+1)*1.D0)/(N2_NT*1.D0) ) *N2_NT-1

  WRITE(*,*)'ND1=',ND1,'ND2=',ND2,'ND3=',ND3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!       
! FIND THE SQUARE OF UNPRECONDITIONED GRADIENT 
!   

  PSI2=DOT_PRODUCT(PSI_C,PSI_C)
  DO I=1,7
    DO J=1,NVCTR_F
      PSI2=PSI2+PSI_F(I,J)**2  
    ENDDO
  ENDDO

  ALLOCATE(AUX(NSTEP_MAX),ALPHAS_PR(NSTEP_MAX)) 
  allocate(GG_c(nvctr_c),stat=i_stat)
  call memocc(i_stat,product(shape(gg_c))*kind(gg_c),'gg_c','precondition')
  allocate(GG_f(7,nvctr_f),stat=i_stat)
  call memocc(i_stat,product(shape(gg_f))*kind(gg_f),'gg_f','precondition')
  allocate(PSI_c_NEW(nvctr_c),stat=i_stat)
  call memocc(i_stat,product(shape(psi_c_new))*kind(psi_c_new),'psi_c_new','precondition')
  allocate(PSI_f_NEW(7,nvctr_f),stat=i_stat)
  call memocc(i_stat,product(shape(psi_f_new))*kind(psi_f_new),'psi_f_new','precondition')
  allocate(HPSI_c(nvctr_c),stat=i_stat)
  call memocc(i_stat,product(shape(hpsi_c))*kind(hpsi_c),'hpsi_c','precondition')
  allocate(HPSI_f(7,nvctr_f),stat=i_stat)
  call memocc(i_stat,product(shape(hpsi_f))*kind(hpsi_f),'hpsi_f','precondition')
  allocate(HPSI_c_OLD(nvctr_c),stat=i_stat)
  call memocc(i_stat,product(shape(hpsi_c_old))*kind(hpsi_c_old),'hpsi_c_old','precondition')
  allocate(HPSI_f_OLD(7,nvctr_f),stat=i_stat)
  call memocc(i_stat,product(shape(hpsi_f_old))*kind(hpsi_f_old),'hpsi_f_old','precondition')
  allocate(HPSI_c_MOD(nvctr_c),stat=i_stat)
  call memocc(i_stat,product(shape(hpsi_c_mod))*kind(hpsi_c_mod),'hpsi_c_mod','precondition')
  allocate(HPSI_f_MOD(7,nvctr_f),stat=i_stat)
  call memocc(i_stat,product(shape(hpsi_f_mod))*kind(hpsi_f_mod),'hpsi_f_mod','precondition')


  GG_C=PSI_C
  GG_F=PSI_F
  
  nl1_c=0 ; nu1_c=n1 
  nl2_c=0 ; nu2_c=n2
  nl3_c=0 ; nu3_c=n3
  nl1_f=nfl1 ; nu1_f=nfu1 
  nl2_f=nfl2 ; nu2_f=nfu2 
  nl3_f=nfl3 ; nu3_f=nfu3
  
  CALL prec_diag(n1,n2,n3,hgrid,nseg_c,nvctr_c,nvctr_f,&
            keyg_c,keyv_c,psi_c,psi_f,C,&
            NUM_TRANS,ND1,ND2,ND3,N2_NT)

  CALL CALC_GRAD_reza(n1,n2,n3,&
            nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c, & 
            nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f, &
            nseg_c,nvctr_c,keyg_c,&
            keyv_c,nseg_f,nvctr_f,&
            keyg_f,keyv_f,psi_c,psi_f, &
            hgrid,gg_c,gg_f,C, &
            hpsi_c_old,hpsi_f_old,hpsi2old)


  ALPHA=PIECELINE(C*HGRID*HGRID*.25d0)*.5d0
  WRITE(*,*)'IORB,C,ALPHA_PR:',IORB,C,ALPHA

  LOOP:DO ISTEP=1,NSTEP_MAX

!    HPSI_MOD=HPSI_OLD
    HPSI_C_MOD=HPSI_C_OLD
    HPSI_F_MOD=HPSI_F_OLD

!    CALL prec_diag_simp(n1,n2,n3,hgrid,C,hpsi_MOD)
  CALL prec_diag(n1,n2,n3,hgrid,nseg_c,nvctr_c,nvctr_f,&
            keyg_c,keyv_c,Hpsi_C_MOD,HPSI_F_MOD,C,&
            NUM_TRANS,ND1,ND2,ND3,N2_NT)

!    PSI_NEW=PSI-ALPHA*HPSI_MOD

    PSI_C_NEW=PSI_C-ALPHA*HPSI_C_MOD
    PSI_F_NEW=PSI_F-ALPHA*HPSI_F_MOD

!    CALL CALC_GRAD(n1,n2,n3,psi_new,hgrid,gg,C,hpsi,hpsi2)           

    CALL CALC_GRAD_reza(n1,n2,n3,&
                   nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c, & 
                   nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f, &
                   nseg_c,nvctr_c,keyg_c,&
                  keyv_c,nseg_f,nvctr_f,&
                  keyg_f,keyv_f,psi_c_NEW,psi_f_NEW, &
                  hgrid,gg_c,gg_f,C, &
                   hpsi_c,hpsi_f,hpsi2)

!   SCALPROD=DOT_PRODUCT(HPSI,HPSI_OLD)
    SCALPROD=DOT_PRODUCT(HPSI_C,HPSI_C_OLD)
    DO I=1,7
      DO J=1,NVCTR_F
        SCALPROD=SCALPROD+HPSI_F(I,J)*HPSI_F_OLD(I,J)
      ENDDO
    ENDDO


    COSANG=SCALPROD/SQRT(HPSI2*HPSI2OLD)
    AUX(ISTEP)=HPSI2
    ALPHAS_PR(ISTEP)=ALPHA
           
    IF (COSANG.GT.COSMIN) THEN
!      PSI=PSI_NEW
      PSI_C=PSI_C_NEW
      PSI_F=PSI_F_NEW
!
      IF ((HPSI2.LT.PSI2*FAC_NORM2).AND.(ISTEP.GT.NSTEP_MIN)) EXIT  
!      HPSI_OLD=HPSI
      HPSI_C_OLD=HPSI_C
      HPSI_F_OLD=HPSI_F
!
      HPSI2OLD=HPSI2
!
      ALPHA=ALPHA*1.1d0
    ELSE
      ALPHA=ALPHA*0.5d0
      IF (ALPHA.LT.ALPHAMIN) THEN
        WRITE(*,*)'ALPHA HIT MINIMUM'
        EXIT LOOP
     ENDIF

    ENDIF  

  ENDDO LOOP

  IF (ISTEP.EQ.NSTEP_MAX+1) ISTEP = NSTEP_MAX

  WRITE(*,*)'IORB,  # OF STEPS: ', IORB,ISTEP
!  WRITE(*,'(a,20(e9.2))')'AUX_GRAD2=',(AUX(I),I=1,ISTEP)
  WRITE(*,'(a,20(e12.5))')'AUX_GRAD2=',(AUX(I),I=1,ISTEP)
  WRITE(*,'(a,20(e9.2))')'ALPHAS_PR==',(ALPHAS_PR(I),I=1,ISTEP)

  DEALLOCATE(GG_C,GG_F,AUX,ALPHAS_PR) 
  DEALLOCATE(HPSI_C,HPSI_F) 
  DEALLOCATE(HPSI_C_OLD,HPSI_F_OLD) 
  DEALLOCATE(HPSI_C_MOD,HPSI_F_MOD) 
  DEALLOCATE(PSI_C_NEW,PSI_F_NEW) 
          
END SUBROUTINE precondition

          
SUBROUTINE CALC_GRAD_REZA(n1,n2,n3,&
       nl1,nu1,nl2,nu2,nl3,nu3, & 
       nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f, &
       nseg_c,nvctr_c,keyg_c,keyv_c,nseg_f,nvctr_f,keyg_f,keyv_f, &
       psi_c,psi_f, &
       hgrid,gg_c,gg_f,C, &
       hpsi_c,hpsi_f,hpsi2) 
        implicit real(kind=8) (a-h,o-z)        

        allocatable psig_stand1(:,:,:,:,:,:),psig_stand2(:,:,:,:,:,:)
                
        dimension keyg_c(2,nseg_c),keyv_c(nseg_c),keyg_f(2,nseg_f),keyv_f(nseg_f)
        dimension psi_c(nvctr_c),psi_f(7,nvctr_f)
        dimension GG_c(nvctr_c),GG_f(7,nvctr_f)
        dimension hpsi_c(nvctr_c),hpsi_f(7,nvctr_f)

        allocate(psig_stand1(0:n1,2,0:n2,2,0:n3,2),stat=i_stat)
        call memocc(i_stat,product(shape(psig_stand1))*kind(psig_stand1),'psig_stand1','precondition')
        allocate(psig_stand2(0:n1,2,0:n2,2,0:n3,2),stat=i_stat)
        call memocc(i_stat,product(shape(psig_stand2))*kind(psig_stand2),'psig_stand2','precondition')

        call uncompress_forstandard(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
                              nseg_c,nvctr_c,keyg_c,keyv_c,  & 
                              nseg_f,nvctr_f,keyg_f,keyv_f,  & 
                              psi_c,psi_f,psig_stand1)

          call ConvolStand(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,&
           nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f, &
           hgrid,psig_stand1,psig_stand2)
          
        call compress_forstandard(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
                            nseg_c,nvctr_c,keyg_c,keyv_c,  & 
                            nseg_f,nvctr_f,keyg_f,keyv_f,  & 
                            psig_stand2,hpsi_c,hpsi_f)

          HPSI_C=HPSI_C+C*PSI_C-GG_C   
          HPSI_F=HPSI_F+C*PSI_F-GG_F            

          HPSI2=DOT_PRODUCT(HPSI_C,HPSI_C)
          DO I=1,7
            DO J=1,NVCTR_F
              HPSI2=HPSI2+HPSI_F(I,J)**2  
            ENDDO
          ENDDO

          i_all=-product(shape(psig_stand1))*kind(psig_stand1)
          deallocate(psig_stand1,stat=i_stat)
          call memocc(i_stat,i_all,'psig_stand1','precondition')
          i_all=-product(shape(psig_stand2))*kind(psig_stand2)
          deallocate(psig_stand2,stat=i_stat)
          call memocc(i_stat,i_all,'psig_stand2','precondition')

END SUBROUTINE CALC_GRAD_REZA


subroutine prec_diag(n1,n2,n3,hgrid,nseg_c,nvctr_c,nvctr_f,&
         keyg_c,keyv_c,hpsi_c,hpsi_f,C,&
         NUM_TRANS,ND1,ND2,ND3,N2_NT)
! 
!
        implicit real(kind=8) (a-h,o-z)
        dimension keyg_c(2,nseg_c),keyv_c(nseg_c),hpsi_c(nvctr_c),hpsi_f(7,nvctr_f)
        real(kind=8), allocatable, dimension(:,:,:) :: hpsip
!
!       WAVELET AND SCALING FUNCTION SECOND DERIVATIVE FILTERS
        PARAMETER(B2=24.8758460293923314D0,A2=3.55369228991319019D0)



        allocate(hpsip(0:nd1,0:nd2,0:nd3),stat=i_stat)
        call memocc(i_stat,product(shape(hpsip))*kind(hpsip),'hpsip','prec_diag')
! find leading dimensions that allow for a wavelet analysis
        call dzero((nd1+1)*(nd2+1)*(nd3+1),hpsip)

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

          FAC_H=1.D0/(HGRID*N2_NT)**2

        H0=    1.5D0*A2*FAC_H;    H1=(A2+B2*.5D0)*FAC_H
        H2=(A2*.5D0+B2)*FAC_H;    H3=    1.5D0*B2*FAC_H
!
        CALL MULTI_FORWARD(ND1,ND2,ND3,HPSIP,NUM_TRANS,NN1,NN2,NN3) 

        NNN1=NN1; NNN2=NN2; NNN3=NN3 

        CALL PRECOND_PROPER(ND1,ND2,ND3,HPSIP,NUM_TRANS,NNN1,NNN2,NNN3,H0,H1,H2,H3,C)

        CALL MULTI_BACKWARD(ND1,ND2,ND3,HPSIP,NUM_TRANS,NN1,NN2,NN3)

        F1=1.D0/(H1+C);         F2=1.D0/(H2+C);         F3=1.D0/(H3+C)

        DO I=1,NVCTR_F
          HPSI_F(1,I)=HPSI_F(1,I)*F1
          HPSI_F(2,I)=HPSI_F(2,I)*F1
          HPSI_F(4,I)=HPSI_F(4,I)*F1

          HPSI_F(3,I)=HPSI_F(3,I)*F2
          HPSI_F(5,I)=HPSI_F(5,I)*F2
          HPSI_F(6,I)=HPSI_F(6,I)*F2

          HPSI_F(7,I)=HPSI_F(7,I)*F3
        ENDDO

!    write(90,*) hpsip

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

END SUBROUTINE


FUNCTION PIECELINE(C)
   IMPLICIT real(kind=8)(A-H,O-Z)

   PARAMETER(A=0.1086D0,B=.9895D0)
   
   WRITE(*,*)'IN PIECELINE, C=',C 
   
   IF (C.LT. .011D0) THEN 
        PIECELINE=.5D0
   ELSE 
        IF (C.LT.1.1D0) THEN 
           PIECELINE=A*LOG(C)+B
        ELSE   
           PIECELINE=1.D0   
        ENDIF        
   ENDIF
     
END 

       
!> Wavefunction in real space
!! Applies the magic filter matrix transposed ; data set shrinks
!! The input array x is overwritten
subroutine convolut_magic_t2(n1,n2,n3,ww,x,y)
        implicit real(kind=8) (a-h,o-z)
        parameter(lowfil=-8,lupfil=8)
        dimension x(-lupfil:n1-lowfil,-lupfil:n2-lowfil,-lupfil:n3-lowfil),y(0:n1,0:n2,0:n3)
        dimension ww(0:n1,-lupfil:n2-lowfil,-lupfil:n3-lowfil)
!          THE MAGIC FILTER FOR DAUBECHIES-16
           real(kind=8), parameter :: fil(lowfil:lupfil)= (/0.D0,&
         2.72734492911979659657715313017228D-6,&
       -0.5185986881173432922848639136911487D-4,&
        0.49443227688689919192282259476750972D-3,&
       -0.344128144493493857280881509686821861D-2,&
        0.1337263414854794752733423467013220997D-1,&
       -0.2103025160930381434955489412839065067D-1,&
       -0.604895289196983516002834636D-1,& 
        0.9940415697834003993178616713D0,&
        0.612625895831207982195380597D-1,&
        0.2373821463724942397566389712597274535D-1,&
       -0.942047030201080385922711540948195075D-2,&
        0.174723713672993903449447812749852942D-2,&
       -0.30158038132690463167163703826169879D-3,&
        0.8762984476210559564689161894116397D-4,&
       -0.1290557201342060969516786758559028D-4,&
        8.4334247333529341094733325815816D-7 /)

!  (I1,I2*I3) -> (I2*I3,i1)
        ndat=(n2+1+lupfil-lowfil)*(n3+1+lupfil-lowfil)
        call convrot_shrink2(lowfil,lupfil,fil,n1,ndat,x,ww)
!  (I2,I3*i1) -> (I3*i1,i2)
        ndat=(n3+1+lupfil-lowfil)*(n1+1)
        call convrot_shrink2(lowfil,lupfil,fil,n2,ndat,ww,x)
!  (I3,i1*i2) -> (i1*i2,i3)
        ndat=(n1+1)*(n2+1)
        call convrot_shrink2(lowfil,lupfil,fil,n3,ndat,x,y)

    return
    end

        subroutine convrot_shrink2(lowfil,lupfil,fil,n1,ndat,x,y)
        implicit real(kind=8) (a-h,o-z)
        dimension fil(lowfil:lupfil),x(lowfil:n1+lupfil,ndat),y(ndat,0:n1)
! the filtered output data structure has shrunk by the filter length

    do j=1,ndat
    do i=0,n1
        
        tt=0.d0
        do l=lowfil,lupfil
    tt=tt+x(i+l,j)*(fil(l))**2
        enddo
    y(j,i)=tt

        enddo
        enddo

end


!> Optimized standard version of the major convolution routines
subroutine ConvolStand(n1,n2,n3,&
               nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,&
               nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,&
               hgrid,x,y)
    !applies the kinetic energy operator onto x to get y
    implicit real(kind=8) (a-h,o-z)
!    dimension nbox_c(2,3),nbox_f(2,3)
    !dimension x(0:n1,2,0:n2,2,0:n3,2),y(0:n1,2,0:n2,2,0:n3,2)
    dimension x(nl1_c:nu1_c,2,nl2_c:nu2_c,2,nl3_c:nu3_c,2)
    dimension y(nl1_c:nu1_c,2,nl2_c:nu2_c,2,nl3_c:nu3_c,2)

    parameter(lowfil=-14,lupfil=14)
    dimension a(lowfil:lupfil),b(lowfil:lupfil),c(lowfil:lupfil),e(lowfil:lupfil)
    scale=-.5d0/hgrid**2
    !write(6,*) 'scale=',scale
    !---------------------------------------------------------------------------
    ! second derivative filters for Daubechies 16
    !<phi|D^2|phi_i>
    a(0)=   -3.5536922899131901941296809374d0*scale
    a(1)=    2.2191465938911163898794546405d0*scale
    a(2)=   -0.6156141465570069496314853949d0*scale
    a(3)=    0.2371780582153805636239247476d0*scale
    a(4)=   -0.0822663999742123340987663521d0*scale
    a(5)=    0.02207029188482255523789911295638968409d0*scale
    a(6)=   -0.409765689342633823899327051188315485d-2*scale
    a(7)=    0.45167920287502235349480037639758496d-3*scale
    a(8)=   -0.2398228524507599670405555359023135d-4*scale
    a(9)=    2.0904234952920365957922889447361d-6*scale
    a(10)=  -3.7230763047369275848791496973044d-7*scale
    a(11)=  -1.05857055496741470373494132287d-8*scale
    a(12)=  -5.813879830282540547959250667d-11*scale
    a(13)=   2.70800493626319438269856689037647576d-13*scale
    a(14)=  -6.924474940639200152025730585882d-18*scale
    do i=1,14
        a(-i)=a(i)
    enddo
    !<phi|D^2|psi_i>
    c(-14)=     -3.869102413147656535541850057188d-18*scale
    c(-13)=      1.5130616560866154733900029272077362d-13*scale
    c(-12)=     -3.2264702314010525539061647271983988409d-11*scale
    c(-11)=     -5.96264938781402337319841002642d-9*scale
    c(-10)=     -2.1656830629214041470164889350342d-7*scale
    c(-9 )=      8.7969704055286288323596890609625d-7*scale
    c(-8 )=     -0.00001133456724516819987751818232711775d0*scale
    c(-7 )=      0.00021710795484646138591610188464622454d0*scale
    c(-6 )=     -0.0021356291838797986414312219042358542d0*scale
    c(-5 )=      0.00713761218453631422925717625758502986d0*scale
    c(-4 )=     -0.0284696165863973422636410524436931061d0*scale
    c(-3 )=      0.14327329352510759457155821037742893841d0*scale
    c(-2 )=     -0.42498050943780130143385739554118569733d0*scale
    c(-1 )=      0.65703074007121357894896358254040272157d0*scale
    c( 0 )=     -0.42081655293724308770919536332797729898d0*scale
    c( 1 )=     -0.21716117505137104371463587747283267899d0*scale
    c( 2 )=      0.63457035267892488185929915286969303251d0*scale
    c( 3 )=     -0.53298223962800395684936080758073568406d0*scale
    c( 4 )=      0.23370490631751294307619384973520033236d0*scale
    c( 5 )=     -0.05657736973328755112051544344507997075d0*scale
    c( 6 )=      0.0080872029411844780634067667008050127d0*scale
    c( 7 )=     -0.00093423623304808664741804536808932984d0*scale
    c( 8 )=      0.00005075807947289728306309081261461095d0*scale
    c( 9 )=     -4.62561497463184262755416490048242d-6*scale
    c( 10)=      6.3919128513793415587294752371778d-7*scale
    c( 11)=      1.87909235155149902916133888931d-8*scale
    c( 12)=      1.04757345962781829480207861447155543883d-10*scale
    c( 13)=     -4.84665690596158959648731537084025836d-13*scale
    c( 14)=      1.2392629629188986192855777620877d-17*scale
    !<psi|D^2|phi_i>
    do i=-14,14
        b(i)=c(-i)
    enddo
    !<psi|D^2|psi_i>
    e(0)=   -24.875846029392331358907766562d0*scale
    e(1)=   -7.1440597663471719869313377994d0*scale
    e(2)=   -0.04251705323669172315864542163525830944d0*scale
    e(3)=   -0.26995931336279126953587091167128839196d0*scale
    e(4)=    0.08207454169225172612513390763444496516d0*scale
    e(5)=   -0.02207327034586634477996701627614752761d0*scale
    e(6)=    0.00409765642831595181639002667514310145d0*scale
    e(7)=   -0.00045167920287507774929432548999880117d0*scale
    e(8)=    0.00002398228524507599670405555359023135d0*scale
    e(9)=   -2.0904234952920365957922889447361d-6*scale
    e(10)=   3.7230763047369275848791496973044d-7*scale
    e(11)=   1.05857055496741470373494132287d-8*scale
    e(12)=   5.8138798302825405479592506674648873655d-11*scale
    e(13)=  -2.70800493626319438269856689037647576d-13*scale
    e(14)=   6.924474940639200152025730585882d-18*scale
    do i=1,14
        e(-i)=e(i)
    enddo
    !---------------------------------------------------------------------------


!     WRITE(*,*)'CONVOL_STAND STARTED,N1,NL1_C:',N1,NL1_C



!    nl1_c=nbox_c(1,1) ; nu1_c=nbox_c(2,1)
!    nl2_c=nbox_c(1,2) ; nu2_c=nbox_c(2,2)
!    nl3_c=nbox_c(1,3) ; nu3_c=nbox_c(2,3)
!    nl1_f=nbox_f(1,1) ; nu1_f=nbox_f(2,1)
!    nl2_f=nbox_f(1,2) ; nu2_f=nbox_f(2,2)
!    nl3_f=nbox_f(1,3) ; nu3_f=nbox_f(2,3)
    y=0.d0
    !nl1_c=0 ; nu1_c=n1
    !nl2_c=0 ; nu2_c=n2
    !nl3_c=0 ; nu3_c=n3
    !nl1_f=0 ; nu1_f=n1
    !nl2_f=0 ; nu2_f=n2
    !nl3_f=0 ; nu3_f=n3
!    write(6,*)  'n1=',n1
!    write(6,*)  'n2=',n2
!    write(6,*)  'n3=',n3
!    write(6,'(2(a,x,i5,x))')  'nl1_c=',nl1_c,'nu1_c=',nu1_c
!    write(6,'(2(a,x,i5,x))')  'nl2_c=',nl2_c,'nu2_c=',nu2_c
!    write(6,'(2(a,x,i5,x))')  'nl3_c=',nl3_c,'nu3_c=',nu3_c
!    write(6,'(2(a,x,i5,x))')  'nl1_f=',nl1_f,'nu1_f=',nu1_f
!    write(6,'(2(a,x,i5,x))')  'nl2_f=',nl2_f,'nu2_f=',nu2_f
!    write(6,'(2(a,x,i5,x))')  'nl3_f=',nl3_f,'nu3_f=',nu3_f

                                                            
    do i3=nl3_c,nu3_c
        !(1/2) d^2/dx^2
        do i2=nl2_c,nu2_c
            do i1=nl1_c,nu1_c
                t111=0.d0;t211=0.d0
                do l=max(nl1_c-i1,lowfil),min(lupfil,nu1_c-i1)
                    t111=t111 + x(i1+l,1,i2,1,i3,1)*a(l) + x(i1+l,2,i2,1,i3,1)*b(l)
                    t211=t211 + x(i1+l,1,i2,1,i3,1)*c(l) + x(i1+l,2,i2,1,i3,1)*e(l)
                enddo
                y(i1,1,i2,1,i3,1)=t111
                y(i1,2,i2,1,i3,1)=t211
            enddo
        enddo
        !+ (1/2) d^2/dy^2
        do i1=nl1_c,nu1_c
            do i2=nl2_c,nu2_c
                t111=0.d0;t121=0.d0
                do l=max(nl2_c-i2,lowfil),min(lupfil,nu2_c-i2)
                    t111=t111 + x(i1,1,i2+l,1,i3,1)*a(l) + x(i1,1,i2+l,2,i3,1)*b(l)
                    t121=t121 + x(i1,1,i2+l,1,i3,1)*c(l) + x(i1,1,i2+l,2,i3,1)*e(l)
                enddo
                y(i1,1,i2,1,i3,1)=y(i1,1,i2,1,i3,1) + t111
                y(i1,1,i2,2,i3,1)=t121
            enddo
        enddo
    enddo
    do i3=nl3_f,nu3_f
        !(1/2) d^2/dx^2
        do i2=nl2_f,nu2_f
            do i1=nl1_f,nu1_f
                t112=0.d0;t121=0.d0;t122=0.d0;t212=0.d0;t221=0.d0;t222=0.d0
                do l=max(nl1_f-i1,lowfil),min(lupfil,nu1_f-i1)
                    t112=t112 + x(i1+l,1,i2,1,i3,2)*a(l) + x(i1+l,2,i2,1,i3,2)*b(l)
                    t121=t121 + x(i1+l,1,i2,2,i3,1)*a(l) + x(i1+l,2,i2,2,i3,1)*b(l)
                    t122=t122 + x(i1+l,1,i2,2,i3,2)*a(l) + x(i1+l,2,i2,2,i3,2)*b(l)
                    t212=t212 + x(i1+l,1,i2,1,i3,2)*c(l) + x(i1+l,2,i2,1,i3,2)*e(l)
                    t221=t221 + x(i1+l,1,i2,2,i3,1)*c(l) + x(i1+l,2,i2,2,i3,1)*e(l)
                    t222=t222 + x(i1+l,1,i2,2,i3,2)*c(l) + x(i1+l,2,i2,2,i3,2)*e(l)
                enddo
                y(i1,1,i2,1,i3,2)=t112
                y(i1,1,i2,2,i3,1)=y(i1,1,i2,2,i3,1) + t121
                y(i1,1,i2,2,i3,2)=t122
                y(i1,2,i2,1,i3,2)=t212
                y(i1,2,i2,2,i3,1)=t221
                y(i1,2,i2,2,i3,2)=t222
            enddo
        enddo
        !+ (1/2) d^2/dy^2
        do i1=nl1_f,nu1_f
            do i2=nl2_f,nu2_f
                t112=0.d0;t211=0.d0;t122=0.d0;t212=0.d0;t221=0.d0;t222=0.d0
                do l=max(nl2_f-i2,lowfil),min(lupfil,nu2_f-i2)
                    t112=t112 + x(i1,1,i2+l,1,i3,2)*a(l) + x(i1,1,i2+l,2,i3,2)*b(l)
                    t211=t211 + x(i1,2,i2+l,1,i3,1)*a(l) + x(i1,2,i2+l,2,i3,1)*b(l)
                    t122=t122 + x(i1,1,i2+l,1,i3,2)*c(l) + x(i1,1,i2+l,2,i3,2)*e(l)
                    t212=t212 + x(i1,2,i2+l,1,i3,2)*a(l) + x(i1,2,i2+l,2,i3,2)*b(l)
                    t221=t221 + x(i1,2,i2+l,1,i3,1)*c(l) + x(i1,2,i2+l,2,i3,1)*e(l)
                    t222=t222 + x(i1,2,i2+l,1,i3,2)*c(l) + x(i1,2,i2+l,2,i3,2)*e(l)
                enddo
                y(i1,1,i2,1,i3,2)=y(i1,1,i2,1,i3,2) + t112
                y(i1,2,i2,1,i3,1)=y(i1,2,i2,1,i3,1) + t211
                y(i1,1,i2,2,i3,2)=y(i1,1,i2,2,i3,2) + t122
                y(i1,2,i2,1,i3,2)=y(i1,2,i2,1,i3,2) + t212
                y(i1,2,i2,2,i3,1)=y(i1,2,i2,2,i3,1) + t221
                y(i1,2,i2,2,i3,2)=y(i1,2,i2,2,i3,2) + t222
            enddo
        enddo
    enddo
    ! + (1/2) d^2/dz^2
    do i2=nl2_c,nu2_c
        do i1=nl1_c,nu1_c
            do i3=nl3_c,nu3_c
                t111=0.d0;t112=0.d0
                do l=max(nl3_c-i3,lowfil),min(lupfil,nu3_c-i3)
                    t111=t111 + x(i1,1,i2,1,i3+l,1)*a(l) + x(i1,1,i2,1,i3+l,2)*b(l)
                    t112=t112 + x(i1,1,i2,1,i3+l,1)*c(l) + x(i1,1,i2,1,i3+l,2)*e(l)
                enddo
                y(i1,1,i2,1,i3,1)=y(i1,1,i2,1,i3,1) + t111
                y(i1,1,i2,1,i3,2)=y(i1,1,i2,1,i3,2) + t112
            enddo
        enddo
    enddo
    do i2=nl2_f,nu2_f
        do i1=nl1_f,nu1_f
            do i3=nl3_f,nu3_f
                t121=0.d0;t211=0.d0;t122=0.d0;t212=0.d0;t221=0.d0;t222=0.d0
                do l=max(nl3_f-i3,lowfil),min(lupfil,nu3_f-i3)
                    t121=t121 + x(i1,1,i2,2,i3+l,1)*a(l) + x(i1,1,i2,2,i3+l,2)*b(l)
                    t211=t211 + x(i1,2,i2,1,i3+l,1)*a(l) + x(i1,2,i2,1,i3+l,2)*b(l)
                    t122=t122 + x(i1,1,i2,2,i3+l,1)*c(l) + x(i1,1,i2,2,i3+l,2)*e(l)
                    t212=t212 + x(i1,2,i2,1,i3+l,1)*c(l) + x(i1,2,i2,1,i3+l,2)*e(l)
                    t221=t221 + x(i1,2,i2,2,i3+l,1)*a(l) + x(i1,2,i2,2,i3+l,2)*b(l)
                    t222=t222 + x(i1,2,i2,2,i3+l,1)*c(l) + x(i1,2,i2,2,i3+l,2)*e(l)
                enddo
                y(i1,1,i2,2,i3,1)=y(i1,1,i2,2,i3,1) + t121
                y(i1,2,i2,1,i3,1)=y(i1,2,i2,1,i3,1) + t211
                y(i1,1,i2,2,i3,2)=y(i1,1,i2,2,i3,2) + t122
                y(i1,2,i2,1,i3,2)=y(i1,2,i2,1,i3,2) + t212
                y(i1,2,i2,2,i3,1)=y(i1,2,i2,2,i3,1) + t221
                y(i1,2,i2,2,i3,2)=y(i1,2,i2,2,i3,2) + t222
            enddo
        enddo
    enddo

!    WRITE(*,*) 'CONVOL_STAND ENDED'

    return
end
    
    
!> Expands the compressed wavefunction in vector form (psi_c,psi_f) 
!! into fine scaling functions (psifscf)
subroutine uncompress_forstandard(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
                              mseg_c,mvctr_c,keyg_c,keyv_c,  & 
                              mseg_f,mvctr_f,keyg_f,keyv_f,  & 
                              psi_c,psi_f,psig)
        implicit real(kind=8) (a-h,o-z)
        dimension keyg_c(2,mseg_c),keyv_c(mseg_c),keyg_f(2,mseg_f),keyv_f(mseg_f)
        dimension psi_c(mvctr_c),psi_f(7,mvctr_f)
        dimension psig(nl1:nu1,2,nl2:nu2,2,nl3:nu3,2)
        !real(kind=8), allocatable :: psig(:,:,:,:,:,:),ww(:)

        !allocate(psig(nl1:nu1,2,nl2:nu2,2,nl3:nu3,2),ww((2*(nu1-nl1)+16)*(2*(nu2-nl2)+16)*(2*(nu3-nl3)+16)))

        call dzero(8*(nu1-nl1+1)*(nu2-nl2+1)*(nu3-nl3+1),psig)

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
              psig(i,1,i2,1,i3,1)=psi_c(i-i0+jj)
           enddo
        enddo

! fine part
        do iseg=1,mseg_f
           jj=keyv_f(iseg)
           j0=keyg_f(1,iseg)
           j1=keyg_f(2,iseg)
           ii=j0-1
           i3=ii/((n1+1)*(n2+1))
           ii=ii-i3*(n1+1)*(n2+1)
           i2=ii/(n1+1)
           i0=ii-i2*(n1+1)
           i1=i0+j1-j0
           do i=i0,i1
              psig(i,2,i2,1,i3,1)=psi_f(1,i-i0+jj)
              psig(i,1,i2,2,i3,1)=psi_f(2,i-i0+jj)
              psig(i,2,i2,2,i3,1)=psi_f(3,i-i0+jj)
              psig(i,1,i2,1,i3,2)=psi_f(4,i-i0+jj)
              psig(i,2,i2,1,i3,2)=psi_f(5,i-i0+jj)
              psig(i,1,i2,2,i3,2)=psi_f(6,i-i0+jj)
              psig(i,2,i2,2,i3,2)=psi_f(7,i-i0+jj)
           enddo
        enddo

        ! calculate fine scaling functions.  It is not needed for standard model.
        !call synthese_grow(nu1-nl1,nu2-nl2,nu3-nl3,ww,psig,psifscf)

        !deallocate(psig,ww)

end


!> Compresses a wavefunction that is given in terms of fine scaling functions (psifscf) into 
!! the retained coarse scaling functions and wavelet coefficients (psi_c,psi_f)
subroutine compress_forstandard(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
                            mseg_c,mvctr_c,keyg_c,keyv_c,  & 
                            mseg_f,mvctr_f,keyg_f,keyv_f,  & 
                            psig,psi_c,psi_f)
        implicit real(kind=8) (a-h,o-z)
        dimension keyg_c(2,mseg_c),keyv_c(mseg_c),keyg_f(2,mseg_f),keyv_f(mseg_f)
        dimension psi_c(mvctr_c),psi_f(7,mvctr_f)
!        dimension psifscf((2*(nu1-nl1)+16)*(2*(nu2-nl2)+16)*(2*(nu3-nl3)+16))
        real(kind=8), allocatable :: ww(:)

        dimension psig(nl1:nu1,2,nl2:nu2,2,nl3:nu3,2)
        
!        allocate(ww((2*(nu1-nl1)+16)*(2*(nu2-nl2)+16)*(2*(nu3-nl3)+16)))
! decompose wavelets into coarse scaling functions and wavelets
!    call analyse_shrink(nu1-nl1,nu2-nl2,nu3-nl3,ww,psifscf,psig)

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
            psi_c(i-i0+jj)=psig(i,1,i2,1,i3,1)
          enddo
        enddo

! fine part
    do iseg=1,mseg_f
          jj=keyv_f(iseg)
          j0=keyg_f(1,iseg)
          j1=keyg_f(2,iseg)
             ii=j0-1
             i3=ii/((n1+1)*(n2+1))
             ii=ii-i3*(n1+1)*(n2+1)
             i2=ii/(n1+1)
             i0=ii-i2*(n1+1)
             i1=i0+j1-j0
      do i=i0,i1
            psi_f(1,i-i0+jj)=psig(i,2,i2,1,i3,1)
            psi_f(2,i-i0+jj)=psig(i,1,i2,2,i3,1)
            psi_f(3,i-i0+jj)=psig(i,2,i2,2,i3,1)
            psi_f(4,i-i0+jj)=psig(i,1,i2,1,i3,2)
            psi_f(5,i-i0+jj)=psig(i,2,i2,1,i3,2)
            psi_f(6,i-i0+jj)=psig(i,1,i2,2,i3,2)
            psi_f(7,i-i0+jj)=psig(i,2,i2,2,i3,2)
          enddo
        enddo

 !       deallocate(ww)

end


SUBROUTINE MULTI_BACKWARD(nd1,nd2,nd3,x,NUM_TRANS,N1,N2,N3)
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
   implicit none
   !Arguments
   integer, intent(in) :: n1,n2,n3,nd1,nd2,nd3,num_trans
   real(kind=8) :: x(0:nd1,0:nd2,0:nd3)
   !Local variables
   real(kind=8), allocatable :: XX(:),YY(:),WW(:) 
   integer :: i,i1,i2,i3,i_trans
   
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

END SUBROUTINE MULTI_FORWARD


!> A periodic synthesis (BACKWARD) wavelet transformation
!! The input array x is not overwritten
subroutine BACKWARD_3D(nd1,nd2,nd3,x,y,ww)
        implicit real(kind=8) (a-h,o-z)
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


!> A periodic synthesis (BACKWARD) wavelet transformation
!! The input array x is not overwritten
subroutine BACKWARD_3D_SELF(nd1,nd2,nd3,x,y,ww)
        implicit real(kind=8) (a-h,o-z)
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


!> An analysis (FORWARD) periodic wavelet transformation
!! The input array y is NOT overwritten
subroutine FORWARD_3D(nd1,nd2,nd3,y,x,ww)
        implicit real(kind=8) (a-h,o-z)
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


!> An analysis (FORWARD) periodic wavelet transformation
!! The input array y is NOT overwritten
subroutine FORWARD_3D_SELF(nd1,nd2,nd3,y,x,ww)
        implicit real(kind=8) (a-h,o-z)
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


!>     FORWARD WAVELET TRANSFORM, ANALYSIS, PERIODIC
SUBROUTINE FORWARD_FAST(RIGHT,NT,C,CD_1)
       implicit real(kind=8) (a-h,o-z)
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


!>  BACKWARD WAVELET TRANSFORM, SYNTHESIS, PERIODIC
SUBROUTINE BACKWARD_FAST(RIGHT1,NT,CD,C1)
      implicit real(kind=8) (a-h,o-z)
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
