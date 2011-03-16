!> @file
!!  Old gauss to daubechies (unused)
!! @deprecated
!! @author
!!    Copyright (C) 2009-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

 !       PROGRAM MAIN
 !       implicit real(kind=8) (a-h,o-z)          
 !       PARAMETER(NMAX=1000,NWORK=10000)
 !       DIMENSION C(0:NMAX,2)
!!
!!       NOW THE DIMENSION OF C IS (0:NMAX,2) INSTEAD OF (-N_INTVX:N_INTVX)
!!
 !       DIMENSION WW(0:NWORK,2)  

 !           GAU_A=1.D0
 !           GAU_CEN=20.D0
 !           N_GAU=0
 !          FACTOR=1.D0
 !
 !        HGRID=1.D0

 !       CALL GAUSS_TO_DAUB(HGRID,FACTOR,GAU_CEN,GAU_A,N_GAU,&!NO ERR, ERRSUC
 !            NMAX,N_LEFT,N_RIGHT,C,ERR_NORM,&               !NO ERR_WAV. NMAX INSTEAD OF N_INTVX
 !            WW,NWORK)                                            !ADDED WORK ARRAY WW(:,:)  


 !           WRITE(*,*)'ERROR=',ERR_NORM

 !       END


         SUBROUTINE GAUSS_TO_DAUB(HGRID,FACTOR,GAU_CEN,GAU_A,N_GAU,&!NO ERR, ERRSUC
              NMAX,N_LEFT,N_RIGHT,C,ERR_NORM,&              !NO ERR_WAV. NMAX INSTEAD OF N_INTVX
              WW,NWORK)                             !ADDED WORK ARRAYS WW WITH DIMENSION  NWORK
! Gives the expansion coefficients of exp(-(1/2)*(x/gau_a)**2)
!!  INPUT: hgrid
!          FACTOR
!!         gau_cen
!!          gau_a
!!          n_gau
!           NMAX
!! OUTPUT: N_LEFT,N_RIGHT: Intervall where the GAUSSIAN IS LARGER THAN
!than thE MACHINE PRECISION
!!         C(:,1) array of scaling function coefficients:
!!         C(:,2) array of wavelet coefficients:
!!         WW(:,1),WW(:,2): work arrays that have to be 16 times larger than C
            implicit real(kind=8) (a-h,o-z)
            INTEGER LEFTS(0:4),RIGHTS(0:4),RIGHTX,LEFTX,RIGHT_T
            DIMENSION C(0:NMAX,2)
            DIMENSION WW(0:NWORK,2)
         INCLUDE 'recs16.inc'! MAGIC FILTER  
         INCLUDE 'intots.inc'! HERE WE KEEP THE ANALYTICAL NORMS OF GAUSSIANS
         INCLUDE 'sym_16.inc'! WAVELET FILTERS

!
!            RESCALE THE PARAMETERS SO THAT HGRID GOES TO 1.D0  
!                   
             A=GAU_A/HGRID
             I0=NINT(GAU_CEN/HGRID) ! THE ARRAY IS CENTERED AT I0
             Z0=GAU_CEN/HGRID-REAL(I0,KIND=8)
         
             H=.125D0*.5d0
!
!            CALCULATE THE ARRAY SIZES;
!            AT LEVEL 0, POSITIONS SHIFTED BY I0 
!
             RIGHT_T= CEILING(15.D0*A)

!             WRITE(*,*)'RIGHT_T=',RIGHT_T,'A=',A,'HGRID=',HGRID,'NMAX=',NMAX

!             IF (2*RIGHT_T.GT.NMAX) &
             IF (2*RIGHT_T.GT.NWORK) &
!             STOP 'A GAUSSIAN IS GREATER THAN THE CELL'
             STOP 'INCREASE THE NWORK IN SUBROUTINE gauss_to_daub.f90'

             LEFTS( 0)=MAX(I0-RIGHT_T,   0)
             RIGHTS(0)=MIN(I0+RIGHT_T,NMAX)

             N_LEFT=LEFTS(0)
             N_RIGHT=RIGHTS(0)

             DO K=1,4
               RIGHTS(K)=2*RIGHTS(K-1)+M
               LEFTS( K)=2*LEFTS( K-1)-M
             ENDDO 

             LEFTX = LEFTS(4)-N
             RIGHTX=RIGHTS(4)+N  
!
!            EIGENTLICH, CALCULATE THE EXPANSION COEFFICIENTS
!            AT LEVEL 4, POSITIONS SHIFTED BY 16*I0 
!         
             !corrected for avoiding 0**0 problem
             if (n_gau == 0) then
                DO I=LEFTX,RIGHTX
                   x=REAL(I-I0*16,KIND=8)*H
                   r=x-z0
                   r2=r/a
                   r2=r2*r2
                   r2=0.5d0*r2
                   func=dexp(-r2)
                   WW(I-LEFTX,1)=func
                ENDDO
             else
                DO I=LEFTX,RIGHTX
                   x=REAL(I-I0*16,KIND=8)*H
                   r=x-z0
                   coeff=r**n_gau
                   r2=r/a
                   r2=r2*r2
                   r2=0.5d0*r2
                   func=dexp(-r2)
                   func=coeff*func
                   WW(I-LEFTX,1)=func
                   !PSI(REAL(I-I0*16,KIND=8)*H,A,Z0,N_GAU)
                ENDDO
             end if

             CALL APPLY_W(WW(:,1),WW(:,2),&
                                LEFTX   ,RIGHTX   ,LEFTS(4),RIGHTS(4),H)

             CALL FORWARD_C(WW(:,2),WW(:,1),&
                                LEFTS(4),RIGHTS(4),LEFTS(3),RIGHTS(3)) 
             CALL FORWARD_C(WW(:,1),WW(:,2),&
                                LEFTS(3),RIGHTS(3),LEFTS(2),RIGHTS(2)) 
             CALL FORWARD_C(WW(:,2),WW(:,1),&
                                LEFTS(2),RIGHTS(2),LEFTS(1),RIGHTS(1)) 

             CALL FORWARD(  WW(:,1),WW(:,2),&
                                LEFTS(1),RIGHTS(1),LEFTS(0),RIGHTS(0)) 

             C=0.D0
             LENGTH=N_RIGHT-N_LEFT+1
             DO I=0,LENGTH-1
               C(I+N_LEFT,1)=WW(I       ,2) !N_LEFT..N_RIGHT <->    0  ..  LENGTH-1
               C(I+N_LEFT,2)=WW(I+LENGTH,2) !N_LEFT..N_RIGHT <-> LENGTH..2*LENGTH-1
             ENDDO 
!
!            CALCULATE THE (RELATIVE) ERROR
!
             CN2=0.D0
             DO I=0,LENGTH*2-1
               CN2=CN2+WW(I,2)**2
             ENDDO 
             
             THEOR_NORM2=VALINTS(N_GAU)*A**(2*N_GAU+1)

             ERROR=SQRT(ABS(1.D0-CN2/THEOR_NORM2))
!
!            RESCALE BACK THE COEFFICIENTS AND THE ERROR
!
             FAC= HGRID**N_GAU*SQRT(HGRID)*FACTOR
             C=C*FAC
             ERR_NORM=ERROR*FAC
!
!            CALCULATE THE OUTPUT ARRAY DIMENSIONS
!

!             WRITE(*,*)'N_LEFT=',N_LEFT,'        N_RIGHT=',N_RIGHT

         END
         
        SUBROUTINE APPLY_W(CX,C,LEFTX,RIGHTX,LEFT,RIGHT,H)
!
!       APPLYING THE MAGIC FILTER ("SHRINK") 
!
        implicit real(kind=8) (a-h,o-z)
        INTEGER RIGHT,RIGHTX
        DIMENSION CX(LEFTX:RIGHTX),C(LEFT:RIGHT)
        INCLUDE 'recs16.inc'

        SQH=SQRT(H)
        
        DO I=LEFT,RIGHT
          CI=0.D0
          DO J=-N,N
            CI=CI+CX(I+J)*W(J)         
          ENDDO
          C(I)=CI*SQH
        ENDDO
        
        END 


      SUBROUTINE FORWARD_C(C,C_1,LEFT,RIGHT,LEFT_1,RIGHT_1)
!
!      FORWARD WAVELET TRANSFORM WITHOUT WAVELETS ("SHRINK")
!
       implicit real(kind=8) (a-h,o-z)
       INTEGER RIGHT,RIGHT_1
       DIMENSION C(LEFT:RIGHT)
       DIMENSION C_1(LEFT_1:RIGHT_1)

       INCLUDE 'sym_16.inc'

!
!      GET THE COARSE SCFUNCTIONS AND WAVELETS
!
       DO I=LEFT_1,RIGHT_1
         I2=2*I
         CI=0.D0
         DO J=-M,M
           CI=CI+CHT(J)*C(J+I2)
         ENDDO
         C_1(I)=CI
       ENDDO

       END

      SUBROUTINE FORWARD(C,CD_1,LEFT,RIGHT,LEFT_1,RIGHT_1)
!
!      CONVENTIONAL FORWARD WAVELET TRANSFORM ("SHRINK")
!
       implicit real(kind=8) (a-h,o-z)
       INTEGER RIGHT,RIGHT_1
       DIMENSION C(LEFT:RIGHT)
       DIMENSION CD_1(LEFT_1:RIGHT_1,2)

       INCLUDE 'sym_16.inc'

!
!      GET THE COARSE SCFUNCTIONS AND WAVELETS
!
       DO I=LEFT_1,RIGHT_1
         I2=2*I
         CI=0.D0
         DI=0.D0
         DO J=-M,M
           CI=CI+CHT(J)*C(J+I2)
           DI=DI+CGT(J)*C(J+I2)
         ENDDO
         CD_1(I,1)=CI
         CD_1(I,2)=DI
       ENDDO
 
       END
