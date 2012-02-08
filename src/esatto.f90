!> @file
!!    Routines to do XANES calculation
!! @author
!!    Copyright (C) 2009-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> module for XANES calculation
module esatto
  use module_base
  use module_interfaces
  implicit none

  integer ::  NLS
  private

  public :: esatto_CalcolaRiflettivita, wig, inizializza, facttable, binary_search
  
  integer  :: ngrid, lpot
  real(gp) , pointer ::  rv(:),  vv(:) , psi(:), derpsi(:)
  real(gp) :: rpot, spot, hpot,Rmts
  complex(gp), pointer :: ham2ham(:,:)
  real(gp) facttable(0:50)
contains
  
  subroutine binomial(n,m,res)
    integer n,m
    real(gp) res
    res = facttable(n)/facttable(n-m)/facttable(m)
  END SUBROUTINE binomial


  real(gp) function Wig( two_ja,  two_jb,  two_jc,&
       two_ma,  two_mb,  two_mc             )
    
    integer two_ja,  two_jb,  two_jc,two_ma,  two_mb,  two_mc
    
    integer  jca ,jcb ,jcc ,jmma ,jmmb, jmmc,jpma,jpmb,jpmc,jsum,kmin, kmax, k, sign !(c) , status
    real(gp)  sum_pos , sum_neg , norm, term
    real(gp)  bc1, bc2, bc3, bcn1, bcn2, bcd1, bcd2, bcd3, bcd4

    if(two_ja < 0 .or. two_jb < 0 .or. two_jc < 0)  then
       print *, " domain error  in Wig"
       stop
    else if(  (two_ja > two_jb+ two_jc)  .or. ((two_jb > two_jc+ two_ja)  .or. two_jc > two_ja+ two_jb)          ) then
       Wig=0.0
    else if (  two_ma+ two_mb+ two_mc /= 0 ) then
       Wig=0
    else if ( abs(two_ma) >two_ja  .or. abs(two_mb) >two_jb  .or. abs(two_mc) >two_jc ) then
       print *, " domain error  in Wig"
       stop
    else if ( mod(two_ma+two_ja,2) /=0   .or. mod(two_mb+two_jb,2) /=0   .or.  mod(two_mc+two_jc,2) /=0  ) then
       print *, " domain error  in Wig"
       stop
    else
       jca  = (-two_ja + two_jb + two_jc) / 2
       jcb  = ( two_ja - two_jb + two_jc) / 2
       jcc  = ( two_ja + two_jb - two_jc) / 2
       jmma = ( two_ja - two_ma) / 2
       jmmb = ( two_jb - two_mb) / 2
       jmmc = ( two_jc - two_mc) / 2
       jpma = ( two_ja + two_ma) / 2
       jpmb = ( two_jb + two_mb) / 2
       jpmc = ( two_jc + two_mc) / 2
       jsum = ( two_ja + two_jb + two_jc) / 2
       kmin = max(0, jpmb - jmmc, jmma - jpmc)
       kmax = min (jcc, jmma, jpmb)
       sign = (-1)**(kmin - jpma + jmmb) 
       !n(c) status = 0
       sum_pos = 0.0
       sum_neg = 0.0
       
       call binomial (two_ja, jcc , bcn1)
       call binomial(two_jb, jcc , bcn2)
       call binomial(jsum+1, jcc , bcd1)
       call binomial(two_ja, jmma, bcd2)
       call binomial(two_jb, jmmb, bcd3)
       call binomial(two_jc, jpmc, bcd4)
       
       norm = sqrt (bcn1 * bcn2) / sqrt (bcd1 * bcd2 * bcd3 * bcd4 * ( two_jc + 1.0))
       
       do k = kmin, kmax
          call binomial (jcc, k, bc1)
          call binomial (jcb, jmma - k, bc2)
          call binomial (jca, jpmb - k, bc3)
          
          
          term = bc1  * bc2  * bc3
          
          if (sign < 0)  then
             sum_neg  = sum_neg  +norm * term
          else 
             sum_pos = sum_pos +norm * term
          endif
          sign = -sign
       enddo
       Wig  = sum_pos - sum_neg
    endif
  end function Wig
  
  
  real(gp) function  ThreeYintegral(Lout,Lpot,Lin,  Mout,Mpot,Min  )
    
    integer , intent(IN) :: Lout,Lpot,Lin,  Mout,Mpot,Min 
    
    real(8) , PARAMETER :: PI=3.141592653589793D0 !n(c) ,TWOPI=2D0*PI
    
    ThreeYintegral  =  sqrt(((2*Lout+1)*(2*Lpot+1)*(2*Lin+1.0))/( 4*PI ))  *Wig(  2*Lout,2*Lpot,2*Lin,  0,0,0)&
         *Wig(  2*Lout,2*Lpot,2*Lin,  -2*Mout,2*Mpot,2*Min)
    
    if( MOD(Mout,2) == 1 ) then
       ThreeYintegral=-ThreeYintegral
    endif
    
  end function ThreeYintegral


  subroutine  inizializza(  nls_a, ngrid_A , rv_A,vv_A,   lpot_A, rpot_A, spot_A, hpot_A,    psi_A, derpsi_A, Rmts_A)
    integer  :: ngrid_A, lpot_A, nls_a
    real(gp) , target ::  rv_A(1:ngrid_A),  vv_A(1:ngrid_A) 
    real(gp) :: rpot_A, spot_A, hpot_A,Rmts_A
    real(gp), pointer :: psi_p(:),  psi_A(:),  derpsi_A(:)
    integer l,j

    NLS=nls_a
    ngrid=ngrid_A
    rv=>rv_A
    vv=>vv_A

    lpot=lpot_A
    rpot=rpot_A
    spot=spot_A
    hpot=hpot_A

    psi=>psi_A
    psi_p=>psi_A
    derpsi=>derpsi_A

    Rmts=Rmts_A

    facttable(0)=1.0
    do l=1,50
       facttable(l)=l*facttable(l-1)
    enddo

    allocate( ham2ham( 0:nls-1, 0:   nls-1    ) )
    
    do l=0, nls-1
       do j=0, nls-1

          ham2ham(l,j) = ThreeYintegral(l,lpot,j,  0,0,0   )

       enddo
    enddo
    facttable(0)=1.0
    do l=1,30
       facttable(l)=l*facttable(l-1)
    enddo

  END SUBROUTINE inizializza


  integer function binary_search(dval, dlist, len)
    !Arguments
    integer :: len
    real(8) :: dval
    real(8) :: dlist(1:len)
    !Local variables
    integer :: bottom , top , middle
    
    if (dval < dlist (1) ) then
       binary_search = 0
    else 
       bottom = 1
       top = len 
       do while (bottom < top)
          middle = (top + bottom) / 2 

          if (dlist (middle) < dval) then
             bottom = middle + 1 
          else if (dlist (middle) > dval) then
             top = middle 
          else
             binary_search=middle
             exit
          endif

       enddo
       if(top==bottom) then
          if (dlist (bottom) > dval) then
             binary_search = bottom - 1 
          else
             binary_search= bottom 
          endif
       endif
    endif
  end function binary_search


  subroutine   getKs( R,E  , Ks , quadra)
    real(gp) R, E
    complex(gp) :: Ks(0:nls-1)
    logical quadra
    integer pos
    real(gp) v0
    integer l

    if(R<=rv(1)) then
       v0= vv(1)
    else if( R>=rv(ngrid)) then
       v0= vv(ngrid)
    else
       pos = binary_search( R, rv, ngrid)
       v0 = ( vv(pos)*(rv(pos+1)-R) + vv(pos+1)*(R-rv(pos)) )/( rv(pos+1) -rv(pos)         )
    endif
    
    do l=0, nls-1
       Ks(l)= 2*E-2*v0- l*(l+1)/R/R
       !! Ks(l)=2*E
       if(.not. quadra) then
          !! print *, "Ks(l)" , Ks(l)
          Ks(l)= sqrt(Ks(l))
          !! print *, "sqrt Ks(l)" , Ks(l)
       endif
    enddo
  END SUBROUTINE getKs


  subroutine getHam2ham(R, h2h_res, fattore)
    real(gp) R, fattore
      complex(gp) :: h2h_res(0:nls-1,0:nls-1)

      real(gp) d,fatt
      integer i,j

      d=R-rpot
      d=(d*d)/(2*spot*spot)
      if(d>60) then
         fatt=0
      else
         fatt = hpot * exp(-d)
      endif
      do i=0, nls-1
         do j =0 , nls-1
            h2h_res(i,j)= fatt* ham2ham(i,j) * fattore
         enddo
      enddo
    END SUBROUTINE getHam2ham


    subroutine RiflettivitaSub( E , Ref, Trasm)
      !Arguments
      real(gp) :: E 
      complex(gp), intent(out) :: Ref(0:nls-1,0:nls-1), Trasm(0:nls-1,0:nls-1)
      !Local variables
      complex(gp), target ::  K(0:nls-1)
      integer :: i
      complex(gp) :: UIC 

      UIC=(0.0,1.0_gp)

      call getKs(RMts,E,K, .false. )

      Ref(:,:)=0.0
      
!!!      print *, " D  psi(1) " ,psi(1)

      do i=0, nls-1
!!!         print *, "per l=  ", i, " psi ", psi(i)
         Ref(i,i) = ( UIC *K(i) *psi(i) + derpsi(i) )/(  UIC *K(i) * psi(i) - derpsi(i) )
      enddo
 
      Trasm(:,:)=0.0
      do i=0, nls-1
         Trasm(i,i) = 2*UIC*K(i) /( UIC *K(i) * psi(i) - derpsi(i) )
      enddo
    END SUBROUTINE RiflettivitaSub


    subroutine Up_Down( R, Energy, E, dE_dz,  Up,  Down)
      real(gp) R, Energy
      complex(gp) ::  E(0:nls-1,0:nls-1) , dE_dz(0:nls-1,0:nls-1) ,  Up(0:nls-1,0:nls-1) , Down(0:nls-1,0:nls-1) 
      
      complex(gp) K(0:nls-1)
      integer i,j
      complex(gp) UIC 
      UIC=(0.0,1.0_gp)
      
      call getKs(R,Energy, K, .false. )
      do i=0, nls-1
         do j=0, nls-1
            Up(i,j)   = ( E(i,j) +    dE_dz(i,j)  /  ( K (i)*UIC)      )*0.5
            Down(i,j) = ( E(i,j) -    dE_dz(i,j)  /  ( K (i)*UIC)      )*0.5
         enddo
      enddo
    END SUBROUTINE Up_Down


    SUBROUTINE ZGEMM_i(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      complex(kind=8) :: ALPHA,BETA
      integer :: K,LDA,LDB,LDC,M,N
      character(len=1) :: TRANSA,TRANSB
      complex(kind=8) :: A(LDA,*),B(LDB,*),C(LDC,*)
      print *, TRANSA,TRANSB,M,N,K,ALPHA,LDA,LDB,BETA,LDC
      call ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    END SUBROUTINE ZGEMM_i


   subroutine Propaga( R1, R2,Energy, E, dE_dz)
     real(gp) :: R1,R2,Energy
     complex(gp) :: E(0:nls-1,0:nls-1), dE_dz(0:nls-1,0:nls-1)

     real(gp) hh,h6
     complex(gp) ::alpha, beta
     complex(gp) :: tf(0:nls-1,0:nls-1)
     complex(gp) :: di1(0:nls-1,0:nls-1)
     complex(gp) :: di2(0:nls-1,0:nls-1)
     complex(gp) :: dt1    (0:nls-1,0:nls-1)
     complex(gp) :: dt2 (0:nls-1,0:nls-1)
     complex(gp) :: dm1    (0:nls-1,0:nls-1)
     complex(gp) :: dm2 (0:nls-1,0:nls-1)

     complex(gp) ::  Ea(0:nls-1,0:nls-1)
     complex(gp) ::  dE_dza(0:nls-1,0:nls-1)
 
     complex(gp) :: K2(0:nls-1)
     integer :: i

     hh=(R2-R1)/2.0
     call  getHam2ham(R1 , tf , 2.0_gp ) 
     call getKs(R1 , Energy, K2 , .true. )
     do i=0, nls-1
        tf(i,i) =  tf(i,i) - K2(i) 
     enddo

     alpha=1.0_gp
     beta=0.0_gp
        
     di1 =  dE_dz 
!!     print *, " propaga 1 ", nls
     call zgemm( 'N','N',  nls,nls,nls  , alpha,    tf(0,0) , nls,  E(0,0), nls , beta,         di2(0,0), nls     ) 


     Ea     =  E    + hh*di1
     dE_dza =  dE_dz+ hh*di2

     call  getHam2ham(R1+hh  , tf , 2.0_gp ) 
     call getKs(R1+hh , Energy, K2 , .true. )
     do i=0, nls-1
        tf(i,i) =  tf(i,i) - K2(i) 
     enddo

     dt1 =  dE_dza 
!!     print *, " propaga 2 "
     call zgemm( 'N','N',  nls,nls,nls  , alpha,    tf(0,0) , nls,  Ea(0,0), nls , beta,         dt2(0,0), nls     ) 

     Ea     =  E    + hh*dt1
     dE_dza =  dE_dz+ hh*dt2

     dm1=dE_dza
     !!     print *, " propaga 3 "

     call zgemm( 'N','N',  nls,nls,nls  , alpha,    tf(0,0) , nls,  Ea(0,0), nls , beta,         dm2(0,0), nls     ) 

     Ea     =  E    + 2*hh*dm1
     dE_dza =  dE_dz+ 2*hh*dm2

     dm1 =dm1+dt1
     dm2 =dm2+dt2

     call  getHam2ham(R2  , tf , 2.0_gp ) 
     call getKs(R2 , Energy, K2 , .true. )
     do i=0, nls-1
        tf(i,i) =  tf(i,i) - K2(i) 
     enddo

     dt1 = dE_dza
     !!     print *, " propaga 4 "

     call zgemm( 'N','N',  nls,nls,nls  , alpha,    tf(0,0) , nls,  Ea(0,0), nls , beta,    dt2(0,0), nls     ) 

     h6=hh/3.0     

     E     =E+h6*((di1+dt1)+   2.0* dm1)
     dE_dz  =dE_dz+ h6*((di2+dt2)+  2.0* dm2)
   END SUBROUTINE Propaga


   function inverse(n,A)
     integer n
     complex(gp) inverse(0:n-1,0:n-1)
     complex(gp) :: A(0:n-1,0:n-1)

     integer INFO,i, IPIV(n)

     inverse=0.0
     do i=0,n-1
        inverse(i,i)=1.0
     enddo

     call ZGESV( n , n , A(0,0), n, IPIV, inverse(0,0), n, INFO )
   end function inverse


   real(gp) function esatto_CalcolaRiflettivita( ngrid_A ,rgrid, dumgrid1, nls_a, lpot_a, rpot_a,spot_a,hpot_a,y_r,d_r,&
        Rmts,    Rinf ,nsteps_coarse ,nsteps_fine, Energia , Labs)
     real(gp) rpot_a, spot_a, hpot_a,Rmts , Rinf, Energia
     integer nls_a, lpot_a, nsteps_coarse, nsteps_fine, ngrid_A
     real(gp), target :: rgrid(1:ngrid_A), dumgrid1(1:ngrid_A)
     real(gp), pointer ::  y_r(:), d_r(:)
     integer Labs

     integer iCs, iFs,i
     real(gp) R0, Rh, rf0, rfh
     complex(gp) dE_dz(  0:nls_a-1,0:nls_a-1       ), E(  0:nls_a-1,0:nls_a-1       )
     complex(gp) UU(  0:nls_a-1,0:nls_a-1       ), DU(  0:nls_a-1,0:nls_a-1       )
     complex(gp) UD(  0:nls_a-1,0:nls_a-1       ), DD(  0:nls_a-1,0:nls_a-1       )
     complex(gp) UIC 
     complex(gp) Ref(0:nls_a-1,0:nls_a-1), Transm(0:nls_a-1,0:nls_a-1), Kh(0:nls_a-1)
!!!     real(gp) , pointer :: py_r(:)

     UIC=(0.0,1.0)

!!!     print *, " B  y_r(1) " ,y_r(1)
!!!     py_r=> y_r
!!!     print *, " B  py_r(1) " ,py_r(1)

     call  inizializza(  nls_a, ngrid_A ,rgrid  ,dumgrid1 ,   lpot_A, rpot_A, spot_A, hpot_A,  y_r  ,d_r ,Rmts )

     call RiflettivitaSub(Energia, Ref, Transm)

     open(unit=22,file='runge.dat')

     do iCs=0, nsteps_coarse-1
        R0 = Rmts +(  (Rinf-Rmts )/NSteps_coarse )*iCs
        Rh = Rmts +(  (Rinf-Rmts )/NSteps_coarse )*(iCs+1)

        call getKs(Rh  ,Energia, Kh,.false. )
        dE_dz = 0.0
        E     = 0.0 
        do i=0, nls-1
           E    (i,i) = 1.0
           dE_dz(i,i) = UIC * Kh(i)  
        enddo

        do iFs=0, nsteps_fine-1
           rfh = Rh +(  (R0-Rh )/nsteps_fine  )*iFs
           rf0 = Rh +(  (R0-Rh )/nsteps_fine  )*(iFs+1)

           call Propaga( rfh, rf0 , Energia , E, dE_dz)

        enddo
        call Up_Down( R0, Energia, E, dE_dz,UU, DU )

        dE_dz = 0.0
        E     = 0.0 
        do i=0, nls-1
           E    (i,i) = 1.0
           dE_dz(i,i) = -UIC * Kh(i)  
        enddo

        do iFs=0, nsteps_fine-1
           rfh = Rh +(  (R0-Rh )/nsteps_fine  )*iFs
           rf0 = Rh +(  (R0-Rh )/nsteps_fine  )*(iFs+1)
            
           call Propaga( rfh, rf0 , Energia , E, dE_dz)

        enddo
        call Up_Down( R0, Energia, E, dE_dz,UD, DD )

        Ref = MatMul(  inverse( nls,  UU-MatMul(Ref, DU)    ), -UD+MatMul(Ref, DD) )
        
        Transm = MatMul(  Transm, DD + MatMul(DU, Ref)      )
        write(22,'(200(f20.10,1x))')  Rh, (1+Ref(1,1))/Transm(1,1), Kh(1) 
        ! print *, Rh, (1+Ref(1,1))/Transm(1,1)
     enddo

     esatto_CalcolaRiflettivita= DBLE(sum( Transm(Labs,:)*conjg(Transm(Labs,:)))     )

     close(unit=22)


   end function esatto_CalcolaRiflettivita
 end module esatto
