!> @file
!!  Old preconditioner routines
!! @deprecated
!! @author
!!    Copyright (C) 2009-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

subroutine preconditionall(iproc,nproc,norb,norbp,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,hgrid,  & 
                   ncong,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,eval,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,hpsi)
! Calls the preconditioner for each orbital treated by the processor
        implicit real(kind=8) (a-h,o-z)
        dimension ibyz_c(2,0:n2,0:n3),ibxz_c(2,0:n1,0:n3),ibxy_c(2,0:n1,0:n2)
        dimension ibyz_f(2,0:n2,0:n3),ibxz_f(2,0:n1,0:n3),ibxy_f(2,0:n1,0:n2)
        dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
        dimension hpsi(nvctr_c+7*nvctr_f,norbp),eval(norb)

       call timing(iproc,'Precondition  ','ON')

      do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)

       cprecr=-eval(iorb)
!       write(*,*) 'cprecr',iorb,cprecr
        call precong(iorb,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, &
                   nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
                   ncong,cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,hpsi(1,iorb-iproc*norbp))

      enddo

       call timing(iproc,'Precondition  ','OF')

END SUBROUTINE preconditionall


subroutine precong(iorb,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, &
                   nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
                   ncong,cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,hpsi)
! Solves (KE+cprecr*I)*xx=yy by conjugate gradient method
! hpsi is the right hand side on input and the solution on output
  implicit none
!Arguments 
  integer :: iorb,n1,n2,n3,ncong,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,nseg_c,nseg_f,nvctr_c
  integer :: nvctr_f
  real(kind=8) :: cprecr,hgrid,hpsi(nvctr_c+7*nvctr_f)
  integer :: ibxy_c(2,0:n1,0:n2),ibxy_f(2,0:n1,0:n2),ibxz_c(2,0:n1,0:n3),ibxz_f(2,0:n1,0:n3)
  integer :: ibyz_c(2,0:n2,0:n3),ibyz_f(2,0:n2,0:n3),keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
!Local variables 
  integer :: i,i_all,i_stat,icong
  real(kind=8) :: a2,alpha,alpha1,alpha2,b2,beta,beta1,beta2,fac_h,h0,h1,h2,h3,tt
  real(kind=8), allocatable :: ppsi(:)
  real(kind=8) :: residues(ncong)
  real(kind=8), allocatable :: rpsi(:)
  real(kind=8) :: scal(0:3)
  real(kind=8), allocatable :: wpsi(:)
!  real(kind=8), allocatable :: spsi(:)

        allocate(rpsi(nvctr_c+7*nvctr_f),stat=i_stat)
        call memocc(i_stat,product(shape(rpsi))*kind(rpsi),'rpsi','precong')
        allocate(ppsi(nvctr_c+7*nvctr_f),stat=i_stat)
        call memocc(i_stat,product(shape(ppsi))*kind(ppsi),'ppsi','precong')
        allocate(wpsi(nvctr_c+7*nvctr_f),stat=i_stat)
        call memocc(i_stat,product(shape(wpsi))*kind(wpsi),'wpsi','precong')

!! check initial residue of original equation
!        allocate(spsi(nvctr_c+7*nvctr_f))
!        do i=1,nvctr_c+7*nvctr_f
!        spsi(i)=hpsi(i)
!        enddo
!        do i=0,3
!        scal(i)=1.d0
!        enddo
!        call  CALC_GRAD_REZA(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, &
!              nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), &
!              scal,cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,hpsi,hpsi(nvctr_c+1),wpsi,wpsi(nvctr_c+1))
!        tt=0.d0
!    do i=1,nvctr_c+7*nvctr_f
!    tt=tt+(wpsi(i)-spsi(i))**2
!        enddo
!        write(*,*) 'initial residue',iorb,sqrt(tt)
!! checkend

!       WAVELET AND SCALING FUNCTION SECOND DERIVATIVES
        B2=24.8758460293923314D0 ; A2=3.55369228991319019D0
        FAC_H=1.d0/hgrid**2
        H0=    1.5D0*A2*FAC_H;    H1=(A2+B2*.5D0)*FAC_H
        H2=(A2*.5D0+B2)*FAC_H;    H3=    1.5D0*B2*FAC_H

        scal(0)=sqrt(1.D0/(H0+cprecr)) ;  scal(1)=sqrt(1.D0/(H1+cprecr)) 
        scal(2)=sqrt(1.D0/(H2+cprecr)) ;  scal(3)=sqrt(1.D0/(H3+cprecr))

        !print "(1x,a,4(1pe26.14))","scal  ",scal
        !print "(1x,a,1pe26.14)",   "cprecr",cprecr

!  D^{-1/2} times right hand side
        call  wscalv(nvctr_c,nvctr_f,scal,hpsi,hpsi(nvctr_c+1))

! assume as input guess x=y
        call  CALC_GRAD_REZA(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, &
              nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), &
              scal,cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,hpsi,hpsi(nvctr_c+1),wpsi,wpsi(nvctr_c+1))
        do i=1,nvctr_c+7*nvctr_f
           tt=wpsi(i)-hpsi(i)
           rpsi(i)=tt
           ppsi(i)=tt
        enddo

        do 1000 icong=2,ncong
        call  CALC_GRAD_REZA(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, &
              nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), &
              scal,cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,ppsi,ppsi(nvctr_c+1),wpsi,wpsi(nvctr_c+1))

        alpha1=0.d0 ; alpha2=0.d0
        do i=1,nvctr_c+7*nvctr_f
           alpha1=alpha1+rpsi(i)*rpsi(i)
           alpha2=alpha2+rpsi(i)*wpsi(i)
        enddo
        residues(icong)=alpha1
        alpha=alpha1/alpha2

        do i=1,nvctr_c+7*nvctr_f
           hpsi(i)=hpsi(i)-alpha*ppsi(i)
           rpsi(i)=rpsi(i)-alpha*wpsi(i)
        end do

        !print "(1x,a,i0,a,1pe26.14)","residues(",icong,")=",residues(icong)
        !print "(1x,a,3(1pe26.14))","alpha",alpha1,alpha2,alpha

!        if (alpha1.lt.tol) goto 1010
        if (icong.ge.ncong) goto 1010

        beta1=0.d0 ; beta2=0.d0
        do i=1,nvctr_c+7*nvctr_f
           beta1=beta1+rpsi(i)*wpsi(i)
           beta2=beta2+ppsi(i)*wpsi(i)
        enddo
        beta=beta1/beta2

        !print "(1x,a,3(1pe26.14))","beta ",alpha1,alpha2,alpha

        do i=1,nvctr_c+7*nvctr_f
           ppsi(i)=rpsi(i)-beta*ppsi(i)
        end do

1000 continue
1010 continue

!  D^{-1/2} times solution
        call  wscalv(nvctr_c,nvctr_f,scal,hpsi,hpsi(nvctr_c+1))

!        write(*,'(i4,(100(1x,e8.2)))') iorb,(residues(icong),icong=2,ncong)

!! check final residue of original equation
!        do i=0,3
!        scal(i)=1.d0
!        enddo
!        call  CALC_GRAD_REZA(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, &
!              nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), &
!              scal,cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,hpsi,hpsi(nvctr_c+1),wpsi,wpsi(nvctr_c+1))
!        tt=0.d0
!    do i=1,nvctr_c+7*nvctr_f
!    tt=tt+(wpsi(i)-spsi(i))**2
!        enddo
!        write(*,*) 'final residue',iorb,sqrt(tt)
!        deallocate(spsi)
!! checkend

        i_all=-product(shape(rpsi))*kind(rpsi)
        deallocate(rpsi,stat=i_stat)
        call memocc(i_stat,i_all,'rpsi','precong')
        i_all=-product(shape(ppsi))*kind(ppsi)
        deallocate(ppsi,stat=i_stat)
        call memocc(i_stat,i_all,'ppsi','precong')
        i_all=-product(shape(wpsi))*kind(wpsi)
        deallocate(wpsi,stat=i_stat)
        call memocc(i_stat,i_all,'wpsi','precong')

END SUBROUTINE precong


subroutine wscalv(nvctr_c,nvctr_f,scal,psi_c,psi_f)
! multiplies a wavefunction psi_c,psi_f (in vector form) with a scaling vector (scal)
        implicit real(kind=8) (a-h,o-z)
        dimension psi_c(nvctr_c),psi_f(7,nvctr_f),scal(0:3)

        do i=1,nvctr_c
           psi_c(i)=psi_c(i)*scal(0)           !  1 1 1
        enddo
        do i=1,nvctr_f
           psi_f(1,i)=psi_f(1,i)*scal(1)       !  2 1 1
           psi_f(2,i)=psi_f(2,i)*scal(1)       !  1 2 1
           psi_f(3,i)=psi_f(3,i)*scal(2)       !  2 2 1
           psi_f(4,i)=psi_f(4,i)*scal(1)       !  1 1 2
           psi_f(5,i)=psi_f(5,i)*scal(2)       !  2 1 2
           psi_f(6,i)=psi_f(6,i)*scal(2)       !  1 2 2
           psi_f(7,i)=psi_f(7,i)*scal(3)       !  2 2 2
        enddo

END SUBROUTINE wscalv


!> ypsi = @f$(1/2) \nabla^2 xpsi@f$
SUBROUTINE CALC_GRAD_REZA(n2,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
                  nseg_c,nvctr_c,keyg_c,keyv_c,nseg_f,nvctr_f,keyg_f,keyv_f, &
                  scal,cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,xpsi_c,xpsi_f,ypsi_c,ypsi_f)
        implicit real(kind=8) (a-h,o-z)        
        dimension keyg_c(2,nseg_c),keyv_c(nseg_c),keyg_f(2,nseg_f),keyv_f(nseg_f)
        dimension xpsi_c(nvctr_c),xpsi_f(7,nvctr_f),scal(0:3)
        dimension ypsi_c(nvctr_c),ypsi_f(7,nvctr_f)
        dimension ibyz_c(2,0:n2,0:n3),ibxz_c(2,0:n1,0:n3),ibxy_c(2,0:n1,0:n2)
        dimension ibyz_f(2,0:n2,0:n3),ibxz_f(2,0:n1,0:n3),ibxy_f(2,0:n1,0:n2)
        allocatable xpsig_c(:,:,:),ypsig_c(:,:,:)
        allocatable xpsig_fc(:,:,:,:)
        allocatable xpsig_f(:,:,:,:),ypsig_f(:,:,:,:)

        allocate(xpsig_c(0:n1,0:n2,0:n3),stat=i_stat)
        call memocc(i_stat,product(shape(xpsig_c))*kind(xpsig_c),'xpsig_c','calc_grad_reza')
        allocate(ypsig_c(0:n1,0:n2,0:n3),stat=i_stat)
        call memocc(i_stat,product(shape(ypsig_c))*kind(ypsig_c),'ypsig_c','calc_grad_reza')
        allocate(xpsig_fc(0:n1,0:n2,0:n3,3),stat=i_stat)
        call memocc(i_stat,product(shape(xpsig_fc))*kind(xpsig_fc),'xpsig_fc','calc_grad_reza')
        allocate(xpsig_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3),stat=i_stat)
        call memocc(i_stat,product(shape(xpsig_f))*kind(xpsig_f),'xpsig_f','calc_grad_reza')
        allocate(ypsig_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3),stat=i_stat)
        call memocc(i_stat,product(shape(ypsig_f))*kind(ypsig_f),'ypsig_f','calc_grad_reza')

        call uncompress_forstandard(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
                              nseg_c,nvctr_c,keyg_c,keyv_c,  & 
                              nseg_f,nvctr_f,keyg_f,keyv_f,  & 
                              scal,xpsi_c,xpsi_f,xpsig_c,xpsig_fc,xpsig_f)

        call Convolkinetic(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
               cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,xpsig_c,xpsig_fc,xpsig_f,ypsig_c,ypsig_f)
          
        call compress_forstandard(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
                            nseg_c,nvctr_c,keyg_c,keyv_c,  & 
                            nseg_f,nvctr_f,keyg_f,keyv_f,  & 
                            scal,ypsig_c,ypsig_f,ypsi_c,ypsi_f)

        i_all=-product(shape(xpsig_c))*kind(xpsig_c)
        deallocate(xpsig_c,stat=i_stat)
        call memocc(i_stat,i_all,'xpsig_c','calc_grad_reza')
        i_all=-product(shape(ypsig_c))*kind(ypsig_c)
        deallocate(ypsig_c,stat=i_stat)
        call memocc(i_stat,i_all,'ypsig_c','calc_grad_reza')
        i_all=-product(shape(xpsig_fc))*kind(xpsig_fc)
        deallocate(xpsig_fc,stat=i_stat)
        call memocc(i_stat,i_all,'xpsig_fc','calc_grad_reza')
        i_all=-product(shape(xpsig_f))*kind(xpsig_f)
        deallocate(xpsig_f,stat=i_stat)
        call memocc(i_stat,i_all,'xpsig_f','calc_grad_reza')
        i_all=-product(shape(ypsig_f))*kind(ypsig_f)
        deallocate(ypsig_f,stat=i_stat)
        call memocc(i_stat,i_all,'ypsig_f','calc_grad_reza')

        END 

    
        subroutine uncompress_forstandard(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
                              mseg_c,nvctr_c,keyg_c,keyv_c,  & 
                              mseg_f,nvctr_f,keyg_f,keyv_f,  & 
                              scal,psi_c,psi_f,psig_c,psig_fc,psig_f) 
! Expands the compressed wavefunction in vector form (psi_c,psi_f) into the psig format
        implicit real(kind=8) (a-h,o-z)
        dimension keyg_c(2,mseg_c),keyv_c(mseg_c),keyg_f(2,mseg_f),keyv_f(mseg_f)
        dimension psi_c(nvctr_c),psi_f(7,nvctr_f),scal(0:3)
        dimension psig_c(0:n1,0:n2,0:n3)
        dimension psig_fc(0:n1,0:n2,0:n3,3),psig_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3)

        call to_zero((n1+1)*(n2+1)*(n3+1),psig_c)
        call to_zero(3*(n1+1)*(n2+1)*(n3+1),psig_fc)
        call to_zero(7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1),psig_f)

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
              psig_c(i,i2,i3)=psi_c(i-i0+jj)*scal(0)
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
            psig_f(1,i,i2,i3)=psi_f(1,i-i0+jj)*scal(1)
            psig_fc(i,i2,i3,1)=psig_f(1,i,i2,i3)
            psig_f(2,i,i2,i3)=psi_f(2,i-i0+jj)*scal(1)
            psig_fc(i,i2,i3,2)=psig_f(2,i,i2,i3)
            psig_f(3,i,i2,i3)=psi_f(3,i-i0+jj)*scal(2)
            psig_f(4,i,i2,i3)=psi_f(4,i-i0+jj)*scal(1)
            psig_fc(i,i2,i3,3)=psig_f(4,i,i2,i3)
            psig_f(5,i,i2,i3)=psi_f(5,i-i0+jj)*scal(2)
            psig_f(6,i,i2,i3)=psi_f(6,i-i0+jj)*scal(2)
            psig_f(7,i,i2,i3)=psi_f(7,i-i0+jj)*scal(3)
          enddo
         enddo


        end

    
        subroutine compress_forstandard(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
                            mseg_c,nvctr_c,keyg_c,keyv_c,  & 
                            mseg_f,nvctr_f,keyg_f,keyv_f,  & 
                            scal,psig_c,psig_f,psi_c,psi_f)
! Compresses a psig wavefunction into psi_c,psi_f form
        implicit real(kind=8) (a-h,o-z)
        dimension keyg_c(2,mseg_c),keyv_c(mseg_c),keyg_f(2,mseg_f),keyv_f(mseg_f)
        dimension psi_c(nvctr_c),psi_f(7,nvctr_f),scal(0:3)
        dimension psig_c(0:n1,0:n2,0:n3)
        dimension psig_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3)
        
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
            psi_c(i-i0+jj)=psig_c(i,i2,i3)*scal(0)
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
            psi_f(1,i-i0+jj)=psig_f(1,i,i2,i3)*scal(1)
            psi_f(2,i-i0+jj)=psig_f(2,i,i2,i3)*scal(1)
            psi_f(3,i-i0+jj)=psig_f(3,i,i2,i3)*scal(2)
            psi_f(4,i-i0+jj)=psig_f(4,i,i2,i3)*scal(1)
            psi_f(5,i-i0+jj)=psig_f(5,i,i2,i3)*scal(2)
            psi_f(6,i-i0+jj)=psig_f(6,i,i2,i3)*scal(2)
            psi_f(7,i-i0+jj)=psig_f(7,i,i2,i3)*scal(3)
          enddo
        enddo

        end


