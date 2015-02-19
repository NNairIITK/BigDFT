!> @file
!!  Routines of compression and uncompression of the wavefunctions
!! @author
!!    Copyright (C) 2010-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Compresses a psig wavefunction into psi_c,psi_f form
subroutine compress_plain(n1,n2,nl1,nu1,nl2,nu2,nl3,nu3, & 
     mseg_c,mvctr_c,keyg_c,keyv_c,  & 
     mseg_f,mvctr_f,keyg_f,keyv_f,  & 
     psig,psi_c,psi_f)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,nl1,nu1,nl2,nu2,nl3,nu3
  integer, intent(in) :: mseg_c,mvctr_c,mseg_f,mvctr_f
  integer, dimension(mseg_c), intent(in) :: keyv_c
  integer, dimension(mseg_f), intent(in) :: keyv_f
  integer, dimension(2,mseg_c), intent(in) :: keyg_c
  integer, dimension(2,mseg_f), intent(in) :: keyg_f
  real(wp), dimension(nl1:nu1,2,nl2:nu2,2,nl3:nu3,2), intent(in) :: psig
  real(wp), dimension(mvctr_c), intent(out) :: psi_c
  real(wp), dimension(7,mvctr_f), intent(out) :: psi_f
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

  !$omp parallel default(private) shared(keyv_c,keyv_f,keyg_c,keyg_f,psig,psi_c,psi_f) &
  !$omp shared(n1,n2,nl1,mseg_c,mseg_f)
  ! coarse part
  !$omp do
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
  !$omp enddo
  ! fine part
  !$omp do
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
  !$omp enddo
  !$omp end parallel

END SUBROUTINE compress_plain


!> Expands the compressed wavefunction in vector form (psi_c,psi_f) 
!! into fine scaling functions (psifscf)
!! The number of operations is 
!! Read IO= nvctr_c +7*nvctr_f
!!
!! Write IO = 8*(n1+1)(n2+1)(n3+1) +  nvctr_c +7*nvctr_f
!!
!! Calculation = 0 (Only data transfer)            
subroutine uncompress(n1,n2,n3,nseg_c,nvctr_c,keyg_c,keyv_c,  & 
     nseg_f,nvctr_f,keyg_f,keyv_f,  & 
     psi_c,psi_f,psig)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(nseg_f), intent(in) :: keyv_f
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  integer, dimension(2,nseg_f), intent(in) :: keyg_f
  real(wp), dimension(nvctr_c), intent(in) :: psi_c
  real(wp), dimension(7,nvctr_f), intent(in) :: psi_f
  real(wp), dimension(0:n1,2,0:n2,2,0:n3,2), intent(out) :: psig
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

  call f_zero(psig)
  !$omp parallel default(private) &
  !$omp shared(psig,psi_c,psi_f,keyv_c,keyg_c,keyv_f,keyg_f,n1,n2,n3,nseg_c,nseg_f)
  

  ! coarse part
  !$omp do
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
        psig(i,1,i2,1,i3,1)=psi_c(i-i0+jj)
     enddo
  enddo
  !$omp enddo

  ! fine part
  !$omp do
  do iseg=1,nseg_f
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
  !$omp enddo

  !$omp end parallel
END SUBROUTINE uncompress


subroutine fill_random(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & !n(c) mvctr_c, mvctr_f (arg:11,15)
     mseg_c,keyg_c,  &
     mseg_f,keyg_f,  & 
     psig_c,psig_f)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,mseg_c,mseg_f !n(c) mvctr_c,mvctr_f 
  integer, dimension(2,mseg_c), intent(in) :: keyg_c
  integer, dimension(2,mseg_f), intent(in) :: keyg_f
  real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: psig_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(out) :: psig_f
  !local variables
  real(wp)::x
  integer :: iseg,j0,j1,ii,i1,i2,i3,i0,i,l !n(c) jj

  psig_c=0._wp
  psig_f=0._wp
  ! coarse part
  do iseg=1,mseg_c
     !n(c) jj=keyv_c(iseg)
     j0=keyg_c(1,iseg)
     j1=keyg_c(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
      call random_number(x) 
        psig_c(i,i2,i3)=x
     enddo
  enddo

  ! fine part
  do iseg=1,mseg_f
     !n(c) jj=keyv_f(iseg)
     j0=keyg_f(1,iseg)
     j1=keyg_f(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
       do l=1,7
          call random_number(x)
          psig_f(l,i,i2,i3)=x
       enddo
     enddo
  enddo

END SUBROUTINE fill_random


!> Expands the compressed wavefunction in vector form (psi_c,psi_f) into the psig format
subroutine uncompress_forstandard_short(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
     mseg_c,mvctr_c,keyg_c,keyv_c,  & 
     mseg_f,mvctr_f,keyg_f,keyv_f,  & 
     scal,psi_c,psi_f,psig_c,psig_f)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,mseg_c,mvctr_c,mseg_f,mvctr_f
  integer, dimension(mseg_c), intent(in) :: keyv_c
  integer, dimension(mseg_f), intent(in) :: keyv_f
  integer, dimension(2,mseg_c), intent(in) :: keyg_c
  integer, dimension(2,mseg_f), intent(in) :: keyg_f
  real(wp), dimension(0:3), intent(in) :: scal
  real(wp), dimension(mvctr_c), intent(in) :: psi_c
  real(wp), dimension(7,mvctr_f), intent(in) :: psi_f
  real(wp), dimension(0:n1,0:n2,0:n3), intent(inout) :: psig_c !these have been already initialised
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(inout) :: psig_f
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

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
!if(abs(psi_c(i-i0+jj)-3.4239806d-05)<1.d-12) then
!    write(950,'(a,6i8,2es16.7)') 'iseg, keyv_c(iseg), keyg_c(1,iseg), i, i2, i3, psi_c(i-i0+jj), abs(psi_c(i-i0+jj)-3.4239806d-05)', &
!        iseg, keyv_c(iseg), keyg_c(1,iseg), i, i2, i3, psi_c(i-i0+jj), abs(psi_c(i-i0+jj)-3.4239806d-05)
!end if
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
        
        psig_f(2,i,i2,i3)=psi_f(2,i-i0+jj)*scal(1)
        
        psig_f(3,i,i2,i3)=psi_f(3,i-i0+jj)*scal(2)
        psig_f(4,i,i2,i3)=psi_f(4,i-i0+jj)*scal(1)
        
        psig_f(5,i,i2,i3)=psi_f(5,i-i0+jj)*scal(2)
        psig_f(6,i,i2,i3)=psi_f(6,i-i0+jj)*scal(2)
        psig_f(7,i,i2,i3)=psi_f(7,i-i0+jj)*scal(3)
     enddo
  enddo

END SUBROUTINE uncompress_forstandard_short


!> Expands the compressed wavefunction in vector form (psi_c,psi_f) into the psig format
subroutine uncompress_forstandard(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
     mseg_c,mvctr_c,keyg_c,keyv_c,  & 
     mseg_f,mvctr_f,keyg_f,keyv_f,  & 
     scal,psi_c,psi_f,psig_c,psig_f,&
     x_f1,x_f2,x_f3)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,mseg_c,mvctr_c,mseg_f,mvctr_f
  integer, dimension(mseg_c), intent(in) :: keyv_c
  integer, dimension(mseg_f), intent(in) :: keyv_f
  integer, dimension(2,mseg_c), intent(in) :: keyg_c
  integer, dimension(2,mseg_f), intent(in) :: keyg_f
  real(wp), dimension(0:3), intent(in) :: scal
  real(wp), dimension(mvctr_c), intent(in) :: psi_c
  real(wp), dimension(7,mvctr_f), intent(in) :: psi_f
  real(wp), dimension(0:n1,0:n2,0:n3), intent(inout) :: psig_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(inout) :: psig_f
  real(wp), dimension(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(inout) :: x_f1
  real(wp), dimension(nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), intent(inout) :: x_f2
  real(wp), dimension(nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), intent(inout) :: x_f3
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

  !$omp parallel default(private) &
  !$omp shared(scal,psig_c,psig_f,x_f1,x_f2,x_f3) &
  !$omp shared(psi_c,psi_f,keyv_c,keyg_c,keyv_f,keyg_f,n1,n2,n3,mseg_c,mseg_f)
  ! coarse part
  !$omp do
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
  !$omp enddo
  ! fine part
  !$omp do
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
        x_f1(i,i2,i3)=psig_f(1,i,i2,i3)
        
        psig_f(2,i,i2,i3)=psi_f(2,i-i0+jj)*scal(1)
        x_f2(i2,i,i3)=psig_f(2,i,i2,i3)
        
        psig_f(3,i,i2,i3)=psi_f(3,i-i0+jj)*scal(2)
        psig_f(4,i,i2,i3)=psi_f(4,i-i0+jj)*scal(1)
        x_f3(i3,i,i2)=psig_f(4,i,i2,i3)
        
        psig_f(5,i,i2,i3)=psi_f(5,i-i0+jj)*scal(2)
        psig_f(6,i,i2,i3)=psi_f(6,i-i0+jj)*scal(2)
        psig_f(7,i,i2,i3)=psi_f(7,i-i0+jj)*scal(3)
     enddo
  enddo
 !$omp enddo
 !$omp end parallel

END SUBROUTINE uncompress_forstandard


subroutine uncompress_f(n1,n2,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
     mseg_f,mvctr_f,keyg_f,keyv_f,  & 
     scal,psi_f,psig_f)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,mseg_f,mvctr_f
  integer, dimension(mseg_f), intent(in) :: keyv_f
  integer, dimension(2,mseg_f), intent(in) :: keyg_f
  real(wp), dimension(0:3), intent(in) :: scal
  real(wp), dimension(7,mvctr_f), intent(in) :: psi_f
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(out) :: psig_f
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

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
        psig_f(2,i,i2,i3)=psi_f(2,i-i0+jj)*scal(1)
        psig_f(3,i,i2,i3)=psi_f(3,i-i0+jj)*scal(2)
        psig_f(4,i,i2,i3)=psi_f(4,i-i0+jj)*scal(1)
        psig_f(5,i,i2,i3)=psi_f(5,i-i0+jj)*scal(2)
        psig_f(6,i,i2,i3)=psi_f(6,i-i0+jj)*scal(2)
        psig_f(7,i,i2,i3)=psi_f(7,i-i0+jj)*scal(3)
     enddo
  enddo

END SUBROUTINE uncompress_f


subroutine compress_f(n1,n2,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
     mseg_f,mvctr_f,keyg_f,keyv_f,  & 
     scal,psig_f,psi_f)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,mseg_f,mvctr_f
  integer, dimension(mseg_f), intent(in) :: keyv_f
  integer, dimension(2,mseg_f), intent(in) :: keyg_f
  real(wp), dimension(0:3), intent(in) :: scal
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in) :: psig_f
  real(wp), dimension(7,mvctr_f), intent(out) :: psi_f
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

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

END SUBROUTINE compress_f


!> Compresses a psig wavefunction into psi_c,psi_f form
subroutine compress_forstandard(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
     mseg_c,mvctr_c,keyg_c,keyv_c,  & 
     mseg_f,mvctr_f,keyg_f,keyv_f,  & 
     scal,psig_c,psig_f,psi_c,psi_f)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,mseg_c,mvctr_c,mseg_f,mvctr_f
  integer, dimension(mseg_c), intent(in) :: keyv_c
  integer, dimension(mseg_f), intent(in) :: keyv_f
  integer, dimension(2,mseg_c), intent(in) :: keyg_c
  integer, dimension(2,mseg_f), intent(in) :: keyg_f
  real(wp), dimension(0:3), intent(in) :: scal
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: psig_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in) :: psig_f
  real(wp), dimension(mvctr_c), intent(out) :: psi_c
  real(wp), dimension(7,mvctr_f), intent(out) :: psi_f
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

  !$omp parallel default(private) &
  !$omp shared(scal,psig_c,psig_f) &
  !$omp shared(psi_c,psi_f,keyv_c,keyg_c,keyv_f,keyg_f,n1,n2,n3,mseg_c,mseg_f)
  ! coarse part
  !$omp do
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
  !$omp enddo
  ! fine part
  !$omp do
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
  !$omp enddo
  !$omp end parallel

END SUBROUTINE compress_forstandard


!> Compresses a wavefunction that is given in terms of fine scaling functions (psifscf) into 
!! the retained coarse scaling functions and wavelet coefficients (psi_c,psi_f)
subroutine compress_per(n1,n2,n3,nseg_c,nvctr_c,keyg_c,keyv_c,  & 
     nseg_f,nvctr_f,keyg_f,keyv_f,  & 
     psifscf,psi_c,psi_f,psig)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(nseg_f), intent(in) :: keyv_f
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  integer, dimension(2,nseg_f), intent(in) :: keyg_f
  real(wp), dimension((2*n1+2)*(2*n2+2)*(2*n3+2)), intent(in) :: psifscf
  real(wp), dimension(0:n1,2,0:n2,2,0:n3,2), intent(inout) :: psig
  real(wp), dimension(nvctr_c), intent(out) :: psi_c
  real(wp), dimension(7,nvctr_f), intent(out) :: psi_f
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

  ! decompose wavelets into coarse scaling functions and wavelets
  call analyse_per_self(n1,n2,n3,psifscf,psig)

  !$omp parallel default(private) &
  !$omp shared(psig,psi_c,psi_f,keyv_c,keyg_c,keyv_f,keyg_f,n1,n2,n3,nseg_c,nseg_f)
  
  ! coarse part
  !$omp do
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
        psi_c(i-i0+jj)=psig(i,1,i2,1,i3,1)
     enddo
  enddo
  !$omp enddo

  ! fine part
  !$omp do
  do iseg=1,nseg_f
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
  !$omp enddo

  !$omp end parallel
END SUBROUTINE compress_per


!> Compresses a wavefunction that is given in terms of fine scaling functions (psifscf) into 
!! the retained coarse scaling functions and wavelet coefficients (psi_c,psi_f)
subroutine compress_per_scal(n1,n2,n3,nseg_c,nvctr_c,keyg_c,keyv_c,  & 
     nseg_f,nvctr_f,keyg_f,keyv_f,  & 
     psifscf,psi_c,psi_f,psig,scal)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(nseg_f), intent(in) :: keyv_f
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  integer, dimension(2,nseg_f), intent(in) :: keyg_f
  real(wp), dimension(0:7), intent(in) :: scal
  real(wp), dimension((2*n1+2)*(2*n2+2)*(2*n3+2)), intent(in) :: psifscf
  real(wp), dimension(0:n1,2,0:n2,2,0:n3,2), intent(inout) :: psig
  real(wp), dimension(nvctr_c), intent(out) :: psi_c
  real(wp), dimension(7,nvctr_f), intent(out) :: psi_f
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

  ! decompose wavelets into coarse scaling functions and wavelets
  call analyse_per_self(n1,n2,n3,psifscf,psig)

  !$omp parallel default(private) &
  !$omp shared(psig,psi_c,psi_f,keyv_c,keyg_c,keyv_f,keyg_f,n1,n2,n3,nseg_c,nseg_f,scal)
  
  ! coarse part
  !$omp do
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
        psi_c(i-i0+jj)=psig(i,1,i2,1,i3,1)*scal(0)
     enddo
  enddo
  !$omp enddo

  ! fine part
  !$omp do
  do iseg=1,nseg_f
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
        psi_f(1,i-i0+jj)=psig(i,2,i2,1,i3,1)*scal(1)
        psi_f(2,i-i0+jj)=psig(i,1,i2,2,i3,1)*scal(2)
        psi_f(3,i-i0+jj)=psig(i,2,i2,2,i3,1)*scal(3)
        psi_f(4,i-i0+jj)=psig(i,1,i2,1,i3,2)*scal(4)
        psi_f(5,i-i0+jj)=psig(i,2,i2,1,i3,2)*scal(5)
        psi_f(6,i-i0+jj)=psig(i,1,i2,2,i3,2)*scal(6)
        psi_f(7,i-i0+jj)=psig(i,2,i2,2,i3,2)*scal(7)
     enddo
  enddo
  !$omp enddo

  !$omp end parallel
END SUBROUTINE compress_per_scal


subroutine compress_scal(n1,n2,n3,nseg_c,nvctr_c,keyg_c,keyv_c,  & 
     nseg_f,nvctr_f,keyg_f,keyv_f,  & 
     psifscf,psi_c,psi_f,scal)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(nseg_f), intent(in) :: keyv_f
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  integer, dimension(2,nseg_f), intent(in) :: keyg_f
  real(wp), dimension(0:7), intent(in) :: scal
  real(wp), dimension(0:n1,2,0:n2,2,0:n3,2), intent(in) :: psifscf
  real(wp), dimension(nvctr_c), intent(out) :: psi_c
  real(wp), dimension(7,nvctr_f), intent(out) :: psi_f
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

  !$omp parallel default(private) &
  !$omp shared(psi_c,psi_f,keyv_c,keyg_c,keyv_f,keyg_f,n1,n2,n3,nseg_c,nseg_f,scal) &
  !$omp shared(psifscf)
  
  ! coarse part
  !$omp do
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
        psi_c(i-i0+jj)=psifscf(i,1,i2,1,i3,1)*scal(0)
     enddo
  enddo
  !$omp enddo

  ! fine part
  !$omp do
  do iseg=1,nseg_f
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
        psi_f(1,i-i0+jj)=psifscf(i,2,i2,1,i3,1)*scal(1)
        psi_f(2,i-i0+jj)=psifscf(i,1,i2,2,i3,1)*scal(2)
        psi_f(3,i-i0+jj)=psifscf(i,2,i2,2,i3,1)*scal(3)
        psi_f(4,i-i0+jj)=psifscf(i,1,i2,1,i3,2)*scal(4)
        psi_f(5,i-i0+jj)=psifscf(i,2,i2,1,i3,2)*scal(5)
        psi_f(6,i-i0+jj)=psifscf(i,1,i2,2,i3,2)*scal(6)
        psi_f(7,i-i0+jj)=psifscf(i,2,i2,2,i3,2)*scal(7)
     enddo
  enddo
  !$omp enddo

  !$omp end parallel
END SUBROUTINE compress_scal


!> Expands the compressed wavefunction in vector form (psi_c,psi_f) 
!! into fine scaling functions (psifscf)
subroutine uncompress_scal(n1,n2,n3,nseg_c,nvctr_c,keyg_c,keyv_c,  & 
     nseg_f,nvctr_f,keyg_f,keyv_f,  & 
     psi_c,psi_f,psifscf,scal)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(nseg_f), intent(in) :: keyv_f
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  integer, dimension(2,nseg_f), intent(in) :: keyg_f
  real(wp), dimension(0:7), intent(in) :: scal
  real(wp), dimension(nvctr_c), intent(in) :: psi_c
  real(wp), dimension(7,nvctr_f), intent(in) :: psi_f
  real(wp), dimension(0:n1,2,0:n2,2,0:n3,2), intent(out) :: psifscf
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

  !$omp parallel default(private) &
  !$omp shared(psi_c,psi_f,keyv_c,keyg_c,keyv_f,keyg_f,n1,n2,n3,nseg_c,nseg_f,scal) &
  !$omp shared(psifscf)
  
  call f_zero(psifscf)

  ! coarse part
  !$omp do
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
        psifscf(i,1,i2,1,i3,1)=psi_c(i-i0+jj)*scal(0)
     enddo
  enddo
  !$omp enddo

  ! fine part
  !$omp do
  do iseg=1,nseg_f
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
        psifscf(i,2,i2,1,i3,1)=psi_f(1,i-i0+jj)*scal(1)
        psifscf(i,1,i2,2,i3,1)=psi_f(2,i-i0+jj)*scal(2)
        psifscf(i,2,i2,2,i3,1)=psi_f(3,i-i0+jj)*scal(3)
        psifscf(i,1,i2,1,i3,2)=psi_f(4,i-i0+jj)*scal(4)
        psifscf(i,2,i2,1,i3,2)=psi_f(5,i-i0+jj)*scal(5)
        psifscf(i,1,i2,2,i3,2)=psi_f(6,i-i0+jj)*scal(6)
        psifscf(i,2,i2,2,i3,2)=psi_f(7,i-i0+jj)*scal(7)
     enddo
  enddo
  !$omp enddo

  !$omp end parallel
  !psig=1.d0
  !psig=1.d0/sqrt(real(8*(n1+1)*(n2+1)*(n3+1),wp))


END SUBROUTINE uncompress_scal


!> Expands the compressed wavefunction in vector form (psi_c,psi_f) 
!! into fine scaling functions (psifscf)
subroutine uncompress_per_scal(n1,n2,n3,nseg_c,nvctr_c,keyg_c,keyv_c,  & 
     nseg_f,nvctr_f,keyg_f,keyv_f,  & 
     psi_c,psi_f,psifscf,psig,scal)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(nseg_f), intent(in) :: keyv_f
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  integer, dimension(2,nseg_f), intent(in) :: keyg_f
  real(wp), dimension(0:7), intent(in) :: scal
  real(wp), dimension(nvctr_c), intent(in) :: psi_c
  real(wp), dimension(7,nvctr_f), intent(in) :: psi_f
  real(wp), dimension(0:n1,2,0:n2,2,0:n3,2), intent(inout) :: psig
  real(wp), dimension((2*n1+2)*(2*n2+2)*(2*n3+2)), intent(out) :: psifscf
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

  call f_zero(psig)

  !$omp parallel default(private) &
  !$omp shared(psig,psi_c,psi_f,keyv_c,keyg_c,keyv_f,keyg_f,n1,n2,n3,nseg_c,nseg_f,scal)

  ! coarse part
  !$omp do
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
        psig(i,1,i2,1,i3,1)=psi_c(i-i0+jj)*scal(0)
     enddo
  enddo
  !$omp enddo

  ! fine part
  !$omp do
  do iseg=1,nseg_f
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
        psig(i,2,i2,1,i3,1)=psi_f(1,i-i0+jj)*scal(1)
        psig(i,1,i2,2,i3,1)=psi_f(2,i-i0+jj)*scal(2)
        psig(i,2,i2,2,i3,1)=psi_f(3,i-i0+jj)*scal(3)
        psig(i,1,i2,1,i3,2)=psi_f(4,i-i0+jj)*scal(4)
        psig(i,2,i2,1,i3,2)=psi_f(5,i-i0+jj)*scal(5)
        psig(i,1,i2,2,i3,2)=psi_f(6,i-i0+jj)*scal(6)
        psig(i,2,i2,2,i3,2)=psi_f(7,i-i0+jj)*scal(7)
     enddo
  enddo
  !$omp enddo

  !$omp end parallel
  !psig=1.d0
  !psig=1.d0/sqrt(real(8*(n1+1)*(n2+1)*(n3+1),wp))

  ! calculate fine scaling functions
  call synthese_per_self(n1,n2,n3,psig,psifscf)

END SUBROUTINE uncompress_per_scal


!> Expands the compressed wavefunction in vector form (psi_c,psi_f) 
!! into fine scaling functions (psifscf)
subroutine uncompress_per(n1,n2,n3,nseg_c,nvctr_c,keyg_c,keyv_c,  & 
     nseg_f,nvctr_f,keyg_f,keyv_f,  & 
     psi_c,psi_f,psifscf,psig)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(nseg_f), intent(in) :: keyv_f
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  integer, dimension(2,nseg_f), intent(in) :: keyg_f
  real(wp), dimension(nvctr_c), intent(in) :: psi_c
  real(wp), dimension(7,nvctr_f), intent(in) :: psi_f
  real(wp), dimension(0:n1,2,0:n2,2,0:n3,2), intent(inout) :: psig
  real(wp), dimension((2*n1+2)*(2*n2+2)*(2*n3+2)), intent(out) :: psifscf
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

  call f_zero(psig)
  !$omp parallel default(private) &
  !$omp shared(psig,psi_c,psi_f,keyv_c,keyg_c,keyv_f,keyg_f,n1,n2,n3,nseg_c,nseg_f)
  

  ! coarse part
  !$omp do
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
        psig(i,1,i2,1,i3,1)=psi_c(i-i0+jj)
     enddo
  enddo
  !$omp enddo

  ! fine part
  !$omp do
  do iseg=1,nseg_f
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
  !$omp enddo

  !$omp end parallel
  !psig=1.d0
  !psig=1.d0/sqrt(real(8*(n1+1)*(n2+1)*(n3+1),wp))

  ! calculate fine scaling functions
  call synthese_per_self(n1,n2,n3,psig,psifscf)

END SUBROUTINE uncompress_per


!> Expands the compressed wavefunction in vector form (psi_c,psi_f) 
!! into fine scaling functions (psifscf)
subroutine uncompress_sd_scal(n1,n2,n3,nseg_c,nvctr_c,keyg_c,keyv_c,  & 
     nseg_f,nvctr_f,keyg_f,keyv_f,  & 
     psi_c,psi_f,psig,scal)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(nseg_f), intent(in) :: keyv_f
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  integer, dimension(2,nseg_f), intent(in) :: keyg_f
  real(gp),intent(in)::scal(0:7)
  real(wp), dimension(nvctr_c), intent(in) :: psi_c
  real(wp), dimension(7,nvctr_f), intent(in) :: psi_f
  real(wp), dimension(0:7,0:n1,0:n2,0:n3), intent(inout) :: psig
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

  call f_zero(psig)
  !$omp parallel default(private) &
  !$omp shared(psig,psi_c,psi_f,keyv_c,keyg_c,keyv_f,keyg_f,n1,n2,n3,nseg_c,nseg_f,scal)
  
  ! coarse part

  !$omp do
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
        psig(0,i,i2,i3)=psi_c(i-i0+jj)*scal(0)
     enddo
  enddo
  !$omp enddo

  ! fine part
  !$omp do
  do iseg=1,nseg_f
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
        psig(1,i,i2,i3)=psi_f(1,i-i0+jj)*scal(1)
        psig(2,i,i2,i3)=psi_f(2,i-i0+jj)*scal(2)
        psig(3,i,i2,i3)=psi_f(3,i-i0+jj)*scal(3)
        psig(4,i,i2,i3)=psi_f(4,i-i0+jj)*scal(4)
        psig(5,i,i2,i3)=psi_f(5,i-i0+jj)*scal(5)
        psig(6,i,i2,i3)=psi_f(6,i-i0+jj)*scal(6)
        psig(7,i,i2,i3)=psi_f(7,i-i0+jj)*scal(7)
     enddo
  enddo

  !$omp enddo

  !$omp end parallel

END SUBROUTINE uncompress_sd_scal


subroutine compress_sd_scal(n1,n2,n3,nseg_c,nvctr_c,keyg_c,keyv_c,  & 
     nseg_f,nvctr_f,keyg_f,keyv_f,  & 
     psig,psi_c,psi_f,scal)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(nseg_f), intent(in) :: keyv_f
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  integer, dimension(2,nseg_f), intent(in) :: keyg_f
  real(gp),intent(in)::scal(0:7)
  real(wp), dimension(0:7,0:n1,0:n2,0:n3), intent(inout) :: psig
  real(wp), dimension(nvctr_c), intent(out) :: psi_c
  real(wp), dimension(7,nvctr_f), intent(out) :: psi_f
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

  !$omp parallel default(private) &
  !$omp shared(psig,psi_c,psi_f,keyv_c,keyg_c,keyv_f,keyg_f,n1,n2,n3,nseg_c,nseg_f,scal)
  
  ! coarse part
  !$omp do
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
        psi_c(i-i0+jj)=psig(0,i,i2,i3)*scal(0)
     enddo
  enddo
  !$omp enddo

  ! fine part
  !$omp do
  do iseg=1,nseg_f
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
        psi_f(1,i-i0+jj)=psig(1,i,i2,i3)*scal(1)
        psi_f(2,i-i0+jj)=psig(2,i,i2,i3)*scal(2)
        psi_f(3,i-i0+jj)=psig(3,i,i2,i3)*scal(3)
        psi_f(4,i-i0+jj)=psig(4,i,i2,i3)*scal(4)
        psi_f(5,i-i0+jj)=psig(5,i,i2,i3)*scal(5)
        psi_f(6,i-i0+jj)=psig(6,i,i2,i3)*scal(6)
        psi_f(7,i-i0+jj)=psig(7,i,i2,i3)*scal(7)
     enddo
  enddo
  !$omp enddo
  !$omp end parallel

END SUBROUTINE compress_sd_scal


!> Expands the compressed wavefunction in vector form (psi_c,psi_f) 
!! into fine scaling functions
subroutine uncompress_sd(n1,n2,n3,nseg_c,nvctr_c,keyg_c,keyv_c,  & 
     nseg_f,nvctr_f,keyg_f,keyv_f,  & 
     psi_c,psi_f,psig)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(nseg_f), intent(in) :: keyv_f
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  integer, dimension(2,nseg_f), intent(in) :: keyg_f
  real(wp), dimension(nvctr_c), intent(in) :: psi_c
  real(wp), dimension(7,nvctr_f), intent(in) :: psi_f
  real(wp), dimension(0:7,0:n1,0:n2,0:n3), intent(out) :: psig
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

  call f_zero(psig)

  !$omp parallel default(private) &
  !$omp shared(psig,psi_c,psi_f,keyv_c,keyg_c,keyv_f,keyg_f,n1,n2,n3,nseg_c,nseg_f)
  
  ! coarse part

  !$omp do
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
        psig(0,i,i2,i3)=psi_c(i-i0+jj)
     enddo
  enddo
  !$omp enddo

  ! fine part
  !$omp do
  do iseg=1,nseg_f
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
        psig(1,i,i2,i3)=psi_f(1,i-i0+jj)
        psig(2,i,i2,i3)=psi_f(2,i-i0+jj)
        psig(3,i,i2,i3)=psi_f(3,i-i0+jj)
        psig(4,i,i2,i3)=psi_f(4,i-i0+jj)
        psig(5,i,i2,i3)=psi_f(5,i-i0+jj)
        psig(6,i,i2,i3)=psi_f(6,i-i0+jj)
        psig(7,i,i2,i3)=psi_f(7,i-i0+jj)
     enddo
  enddo

  !$omp enddo

  !$omp end parallel

END SUBROUTINE uncompress_sd


subroutine compress_sd(n1,n2,n3,nseg_c,nvctr_c,keyg_c,keyv_c,  & 
     nseg_f,nvctr_f,keyg_f,keyv_f,  & 
     psig,psi_c,psi_f)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(nseg_f), intent(in) :: keyv_f
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  integer, dimension(2,nseg_f), intent(in) :: keyg_f
  real(wp), dimension(0:7,0:n1,0:n2,0:n3), intent(in) :: psig
  real(wp), dimension(nvctr_c), intent(out) :: psi_c
  real(wp), dimension(7,nvctr_f), intent(out) :: psi_f
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

  !$omp parallel default(private) &
  !$omp shared(psig,psi_c,psi_f,keyv_c,keyg_c,keyv_f,keyg_f,n1,n2,n3,nseg_c,nseg_f)
  
  ! coarse part
  !$omp do
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
        psi_c(i-i0+jj)=psig(0,i,i2,i3)
     enddo
  enddo
  !$omp enddo

  ! fine part
  !$omp do
  do iseg=1,nseg_f
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
        psi_f(1,i-i0+jj)=psig(1,i,i2,i3)
        psi_f(2,i-i0+jj)=psig(2,i,i2,i3)
        psi_f(3,i-i0+jj)=psig(3,i,i2,i3)
        psi_f(4,i-i0+jj)=psig(4,i,i2,i3)
        psi_f(5,i-i0+jj)=psig(5,i,i2,i3)
        psi_f(6,i-i0+jj)=psig(6,i,i2,i3)
        psi_f(7,i-i0+jj)=psig(7,i,i2,i3)
     enddo
  enddo
  !$omp enddo
  !$omp end parallel

END SUBROUTINE compress_sd


subroutine uncompress_c(hpsi,x_c,keyg_c,keyv_c,nseg_c,nvctr_c,n1,n2,n3)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  integer, intent(in) :: nseg_c,nvctr_c
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  real(wp), dimension(nvctr_c), intent(in) :: hpsi
  real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: x_c
  !local variables
  integer iseg,jj,j0,j1,ii,i3,i2,i0,i1,i

  call f_zero(x_c)
  !$omp parallel default(private) &
  !$omp shared(hpsi,x_c,keyv_c,keyg_c,n1,n2,n3,nseg_c)
  
  !$omp do
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
        x_c(i,i2,i3)=hpsi(i-i0+jj)
     enddo
  enddo
  !$omp enddo
  !$omp end parallel
  
END SUBROUTINE uncompress_c


subroutine compress_c(hpsi,y_c,keyg_c,keyv_c,nseg_c,nvctr_c,n1,n2,n3)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  integer, intent(in) :: nseg_c,nvctr_c
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: y_c
  real(wp), dimension(nvctr_c), intent(out) :: hpsi
  !local variables
  integer iseg,jj,j0,j1,ii,i3,i2,i0,i1,i

  !$omp parallel default(private) &
  !$omp shared(hpsi,y_c,keyv_c,keyg_c,n1,n2,n3,nseg_c)
  !$omp do
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
        hpsi(i-i0+jj)=y_c(i,i2,i3)
     enddo
  enddo
  !$omp enddo

  !$omp end parallel
END SUBROUTINE compress_c


!> Compresses a wavefunction that is given in terms of fine scaling functions (psifscf) into 
!! the retained coarse scaling functions and wavelet coefficients (psi_c,psi_f)
subroutine compress_slab_scal(n1,n2,n3,nseg_c,nvctr_c,keyg_c,keyv_c,  & 
     nseg_f,nvctr_f,keyg_f,keyv_f,  & 
     psifscf,psi_c,psi_f,psig,scal)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(nseg_f), intent(in) :: keyv_f
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  integer, dimension(2,nseg_f), intent(in) :: keyg_f
  real(wp), dimension(0:7), intent(in) :: scal
  real(wp), dimension((2*n1+2)*(2*n2+16)*(2*n3+2)), intent(in) :: psifscf
  real(wp), dimension(0:n1,2,0:n2,2,0:n3,2), intent(inout) :: psig
  real(wp), dimension(nvctr_c), intent(out) :: psi_c
  real(wp), dimension(7,nvctr_f), intent(out) :: psi_f
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

  ! decompose wavelets into coarse scaling functions and wavelets

  call analyse_slab_self(n1,n2,n3,psifscf,psig)

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
        psi_c(i-i0+jj)=psig(i,1,i2,1,i3,1)*scal(0)
     enddo
  enddo

  ! fine part
  do iseg=1,nseg_f
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
        psi_f(1,i-i0+jj)=psig(i,2,i2,1,i3,1)*scal(1)
        psi_f(2,i-i0+jj)=psig(i,1,i2,2,i3,1)*scal(2)
        psi_f(3,i-i0+jj)=psig(i,2,i2,2,i3,1)*scal(3)
        psi_f(4,i-i0+jj)=psig(i,1,i2,1,i3,2)*scal(4)
        psi_f(5,i-i0+jj)=psig(i,2,i2,1,i3,2)*scal(5)
        psi_f(6,i-i0+jj)=psig(i,1,i2,2,i3,2)*scal(6)
        psi_f(7,i-i0+jj)=psig(i,2,i2,2,i3,2)*scal(7)
     enddo
  enddo

END SUBROUTINE compress_slab_scal


!> Expands the compressed wavefunction in vector form (psi_c,psi_f) 
!! into fine scaling functions (psifscf)
subroutine uncompress_slab_scal(n1,n2,n3,nseg_c,nvctr_c,keyg_c,keyv_c,  & 
     nseg_f,nvctr_f,keyg_f,keyv_f,  & 
     psi_c,psi_f,psifscf,psig,scal)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(nseg_f), intent(in) :: keyv_f
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  integer, dimension(2,nseg_f), intent(in) :: keyg_f
  real(wp), dimension(0:7), intent(in) :: scal
  real(wp), dimension(nvctr_c), intent(in) :: psi_c
  real(wp), dimension(7,nvctr_f), intent(in) :: psi_f
  real(wp), dimension(0:n1,2,0:n2,2,0:n3,2), intent(inout) :: psig
  real(wp), dimension((2*n1+2)*(2*n2+16)*(2*n3+2)), intent(out) :: psifscf
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

  call f_zero(psig)

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
        psig(i,1,i2,1,i3,1)=psi_c(i-i0+jj)*scal(0)
     enddo
  enddo

  ! fine part
  do iseg=1,nseg_f
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
        psig(i,2,i2,1,i3,1)=psi_f(1,i-i0+jj)*scal(1)
        psig(i,1,i2,2,i3,1)=psi_f(2,i-i0+jj)*scal(2)
        psig(i,2,i2,2,i3,1)=psi_f(3,i-i0+jj)*scal(3)
        psig(i,1,i2,1,i3,2)=psi_f(4,i-i0+jj)*scal(4)
        psig(i,2,i2,1,i3,2)=psi_f(5,i-i0+jj)*scal(5)
        psig(i,1,i2,2,i3,2)=psi_f(6,i-i0+jj)*scal(6)
        psig(i,2,i2,2,i3,2)=psi_f(7,i-i0+jj)*scal(7)
     enddo
  enddo

  ! calculate fine scaling functions
  
  call synthese_slab_self(n1,n2,n3,psig,psifscf)

END SUBROUTINE uncompress_slab_scal


!> Compresses a wavefunction that is given in terms of fine scaling functions (psifscf) into 
!! the retained coarse scaling functions and wavelet coefficients (psi_c,psi_f)
subroutine compress_slab(n1,n2,n3,nseg_c,nvctr_c,keyg_c,keyv_c,  & 
     nseg_f,nvctr_f,keyg_f,keyv_f,  & 
     psifscf,psi_c,psi_f,psig)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(nseg_f), intent(in) :: keyv_f
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  integer, dimension(2,nseg_f), intent(in) :: keyg_f
  real(wp), dimension((2*n1+2)*(2*n2+16)*(2*n3+2)), intent(inout) :: psifscf
  real(wp), dimension(0:n1,2,0:n2,2,0:n3,2), intent(inout) :: psig
  real(wp), dimension(nvctr_c), intent(out) :: psi_c
  real(wp), dimension(7,nvctr_f), intent(out) :: psi_f
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

  ! decompose wavelets into coarse scaling functions and wavelets

  call analyse_slab_self(n1,n2,n3,psifscf,psig)

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
        psi_c(i-i0+jj)=psig(i,1,i2,1,i3,1)
     enddo
  enddo

  ! fine part
  do iseg=1,nseg_f
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

END SUBROUTINE compress_slab


!> Expands the compressed wavefunction in vector form (psi_c,psi_f) 
!! into fine scaling functions (psifscf)
subroutine uncompress_slab(n1,n2,n3,nseg_c,nvctr_c,keyg_c,keyv_c,  & 
     nseg_f,nvctr_f,keyg_f,keyv_f,  & 
     psi_c,psi_f,psifscf,psig)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(nseg_f), intent(in) :: keyv_f
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  integer, dimension(2,nseg_f), intent(in) :: keyg_f
  real(wp), dimension(nvctr_c), intent(in) :: psi_c
  real(wp), dimension(7,nvctr_f), intent(in) :: psi_f
  real(wp), dimension(0:n1,2,0:n2,2,0:n3,2), intent(inout) :: psig
  real(wp), dimension((2*n1+2)*(2*n2+16)*(2*n3+2)), intent(out) :: psifscf
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

  call f_zero(psig)

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
        psig(i,1,i2,1,i3,1)=psi_c(i-i0+jj)
     enddo
  enddo

  ! fine part
  do iseg=1,nseg_f
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

  ! calculate fine scaling functions
  
  call synthese_slab_self(n1,n2,n3,psig,psifscf)

END SUBROUTINE uncompress_slab


!> Compresses a wavefunction that is given in terms of scfunctions (y_c)
!! and wavelets (y_f) into 
!! the retained coarse scaling functions and wavelet coefficients (psi_c,psi_f)
subroutine compress_per_f(n1,n2,n3,nseg_c,nvctr_c,keyg_c,keyv_c,  & 
     nseg_f,nvctr_f,keyg_f,keyv_f,  & 
     y_c,y_f,psi_c,psi_f,min1,min2,min3,max1,max2,max3)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f
  integer, intent(in) :: min1,min2,min3,max1,max2,max3
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(nseg_f), intent(in) :: keyv_f
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  integer, dimension(2,nseg_f), intent(in) :: keyg_f
  real(wp),intent(in)::y_f(7,min1:max1,min2:max2,min3:max3)
  real(wp),intent(in)::y_c(0:n1,0:n2,0:n3)

  real(wp), dimension(nvctr_c), intent(out) :: psi_c
  real(wp), dimension(7,nvctr_f), intent(out) :: psi_f
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

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
        psi_c(i-i0+jj)=y_c(i,i2,i3)
     enddo
  enddo

  ! fine part
  do iseg=1,nseg_f
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
        psi_f(1,i-i0+jj)=y_f(1,i,i2,i3)
        psi_f(2,i-i0+jj)=y_f(2,i,i2,i3)
        psi_f(3,i-i0+jj)=y_f(3,i,i2,i3)
        psi_f(4,i-i0+jj)=y_f(4,i,i2,i3)
        psi_f(5,i-i0+jj)=y_f(5,i,i2,i3)
        psi_f(6,i-i0+jj)=y_f(6,i,i2,i3)
        psi_f(7,i-i0+jj)=y_f(7,i,i2,i3)
     enddo
  enddo

END SUBROUTINE compress_per_f


!> Expands the compressed wavefunction in vector form (psi_c,psi_f) 
!! into fine scaling functions (psifscf)
!! in addition, calculates the wavelet coefficient array x_f
subroutine uncompress_per_f(n1,n2,n3,nseg_c,nvctr_c,keyg_c,keyv_c,  & 
     nseg_f,nvctr_f,keyg_f,keyv_f,  & 
     psi_c,psi_f,x_c,x_f,x_f1,x_f2,x_f3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f
  integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(nseg_f), intent(in) :: keyv_f
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  integer, dimension(2,nseg_f), intent(in) :: keyg_f
  real(wp), dimension(nvctr_c), intent(in) :: psi_c
  real(wp), dimension(7,nvctr_f), intent(in) :: psi_f
  real(wp),intent(out)::x_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3)
  real(wp),intent(out)::x_c(0:n1,0:n2,0:n3)
  real(wp), dimension(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(out) :: x_f1
  real(wp), dimension(nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), intent(out) :: x_f2
  real(wp), dimension(nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), intent(out) :: x_f3
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

  x_c=0._wp
  x_f=0._wp
  x_f1=0._wp
  x_f2=0._wp
  x_f3=0._wp

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
        x_c(i,i2,i3)       =psi_c(i-i0+jj)
     enddo
  enddo

  ! fine part
  do iseg=1,nseg_f
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
      
      x_f(1,i,i2,i3)=psi_f(1,i-i0+jj)
      x_f(2,i,i2,i3)=psi_f(2,i-i0+jj)
      x_f(3,i,i2,i3)=psi_f(3,i-i0+jj)
      x_f(4,i,i2,i3)=psi_f(4,i-i0+jj)
      x_f(5,i,i2,i3)=psi_f(5,i-i0+jj)
      x_f(6,i,i2,i3)=psi_f(6,i-i0+jj)
      x_f(7,i,i2,i3)=psi_f(7,i-i0+jj)
      
        x_f1(i,i2,i3)=psi_f(1,i-i0+jj)
        x_f2(i2,i,i3)=psi_f(2,i-i0+jj)
        x_f3(i3,i,i2)=psi_f(4,i-i0+jj)
     enddo
  enddo

END SUBROUTINE uncompress_per_f


!> Expands the compressed wavefunction in vector form (psi_c,psi_f) 
!! into fine scaling functions (psifscf)
!! in addition, calculates the wavelet coefficient array x_f
subroutine uncompress_per_f_short(n1,n2,n3,nseg_c,nvctr_c,keyg_c,keyv_c,  & 
     nseg_f,nvctr_f,keyg_f,keyv_f,  & 
     psi_c,psi_f,x_c,x_f,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f
  integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(nseg_f), intent(in) :: keyv_f
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  integer, dimension(2,nseg_f), intent(in) :: keyg_f
  real(wp), dimension(nvctr_c), intent(in) :: psi_c
  real(wp), dimension(7,nvctr_f), intent(in) :: psi_f
  real(wp),intent(out)::x_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3)
  real(wp),intent(out)::x_c(0:n1,0:n2,0:n3)
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

  x_c=0._wp
  x_f=0._wp

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
        x_c(i,i2,i3) = psi_c(i-i0+jj)
     enddo
  enddo

  ! fine part
  do iseg=1,nseg_f
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
      
      x_f(1,i,i2,i3)=psi_f(1,i-i0+jj)
      x_f(2,i,i2,i3)=psi_f(2,i-i0+jj)
      x_f(3,i,i2,i3)=psi_f(3,i-i0+jj)
      x_f(4,i,i2,i3)=psi_f(4,i-i0+jj)
      x_f(5,i,i2,i3)=psi_f(5,i-i0+jj)
      x_f(6,i,i2,i3)=psi_f(6,i-i0+jj)
      x_f(7,i,i2,i3)=psi_f(7,i-i0+jj)
     enddo
  enddo

END SUBROUTINE uncompress_per_f_short


!> Expands the compressed wavefunction in vector form (psi_c,psi_f) into the psig format
!! note that psig should be put to zero outside segment regions
subroutine uncompress_standard_scal(grid,wfd,scal,keyv_c,keyv_f,keyg_c,keyg_f,&
     psi_c,psi_f,psig_c,psig_f)
  use module_base
  use module_types
  implicit none
  type(grid_dimensions), intent(in) :: grid
  type(wavefunctions_descriptors), intent(in) :: wfd !< to be used only for the dimensions
  real(wp), dimension(0:3), intent(in) :: scal
  integer, dimension(wfd%nseg_c), intent(in) :: keyv_c
  integer, dimension(wfd%nseg_f), intent(in) :: keyv_f
  integer, dimension(2,wfd%nseg_c), intent(in) :: keyg_c
  integer, dimension(2,wfd%nseg_f), intent(in) :: keyg_f
  real(wp), dimension(wfd%nvctr_c), intent(in) :: psi_c
  real(wp), dimension(7,wfd%nvctr_f), intent(in) :: psi_f
  real(wp), dimension(0:grid%n1,0:grid%n2,0:grid%n3), intent(inout) :: psig_c
  real(wp), dimension(7,grid%nfl1:grid%nfu1,grid%nfl2:grid%nfu2,grid%nfl3:grid%nfu3), intent(inout) :: psig_f
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

  ! coarse part
  !$omp parallel default(shared) &
  !$omp private(iseg,jj,j0,j1,ii,i1,i2,i3,i0,i)
  
  !$omp do
  do iseg=1,wfd%nseg_c
     jj=keyv_c(iseg)
     j0=keyg_c(1,iseg)
     j1=keyg_c(2,iseg)
     ii=j0-1
     i3=ii/((grid%n1+1)*(grid%n2+1))
     ii=ii-i3*(grid%n1+1)*(grid%n2+1)
     i2=ii/(grid%n1+1)
     i0=ii-i2*(grid%n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psig_c(i,i2,i3)=psi_c(i-i0+jj)*scal(0)
     enddo
  enddo
  !$omp enddo
  ! fine part, to be done only if nseg_f is nonzero
  !$omp do
  do iseg=1,wfd%nseg_f
     jj=keyv_f(iseg)
     j0=keyg_f(1,iseg)
     j1=keyg_f(2,iseg)
     ii=j0-1
     i3=ii/((grid%n1+1)*(grid%n2+1))
     ii=ii-i3*(grid%n1+1)*(grid%n2+1)
     i2=ii/(grid%n1+1)
     i0=ii-i2*(grid%n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psig_f(1,i,i2,i3)=psi_f(1,i-i0+jj)*scal(1)
        psig_f(2,i,i2,i3)=psi_f(2,i-i0+jj)*scal(1)
        psig_f(3,i,i2,i3)=psi_f(3,i-i0+jj)*scal(2)
        psig_f(4,i,i2,i3)=psi_f(4,i-i0+jj)*scal(1)
        psig_f(5,i,i2,i3)=psi_f(5,i-i0+jj)*scal(2)
        psig_f(6,i,i2,i3)=psi_f(6,i-i0+jj)*scal(2)
        psig_f(7,i,i2,i3)=psi_f(7,i-i0+jj)*scal(3)
     enddo
  enddo
 !$omp enddo
 !$omp end parallel

END SUBROUTINE uncompress_standard_scal

subroutine compress_standard_scal(grid,wfd,scal,keyv_c,keyv_f,keyg_c,keyg_f,&
     psig_c,psig_f,psi_c,psi_f)
  use module_base
  use module_types
  implicit none
  type(grid_dimensions), intent(in) :: grid
  type(wavefunctions_descriptors), intent(in) :: wfd
  real(wp), dimension(0:3), intent(in) :: scal
  integer, dimension(wfd%nseg_c), intent(in) :: keyv_c
  integer, dimension(wfd%nseg_f), intent(in) :: keyv_f
  integer, dimension(2,wfd%nseg_c), intent(in) :: keyg_c
  integer, dimension(2,wfd%nseg_f), intent(in) :: keyg_f
  real(wp), dimension(0:grid%n1,0:grid%n2,0:grid%n3), intent(in) :: psig_c
  real(wp), dimension(7,grid%nfl1:grid%nfu1,grid%nfl2:grid%nfu2,grid%nfl3:grid%nfu3), intent(in) :: psig_f
  real(wp), dimension(wfd%nvctr_c), intent(inout) :: psi_c
  real(wp), dimension(7,wfd%nvctr_f), intent(inout) :: psi_f
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

  !$omp parallel default(shared) &
  !$omp private(iseg,jj,j0,j1,ii,i1,i2,i3,i0,i)

  ! coarse part
  !$omp do
  do iseg=1,wfd%nseg_c
     jj=keyv_c(iseg)
     j0=keyg_c(1,iseg)
     j1=keyg_c(2,iseg)
     ii=j0-1
     i3=ii/((grid%n1+1)*(grid%n2+1))
     ii=ii-i3*(grid%n1+1)*(grid%n2+1)
     i2=ii/(grid%n1+1)
     i0=ii-i2*(grid%n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psi_c(i-i0+jj)=psig_c(i,i2,i3)*scal(0)
     enddo
  enddo
  !$omp enddo
  ! fine part
  !$omp do
  do iseg=1,wfd%nseg_f
     jj=keyv_f(iseg)
     j0=keyg_f(1,iseg)
     j1=keyg_f(2,iseg)
     ii=j0-1
     i3=ii/((grid%n1+1)*(grid%n2+1))
     ii=ii-i3*(grid%n1+1)*(grid%n2+1)
     i2=ii/(grid%n1+1)
     i0=ii-i2*(grid%n1+1)
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
  !$omp enddo
  !$omp end parallel

end subroutine compress_standard_scal


!> Compress the wavefunction psig and accumulate the result on the psi array
!! The wavefunction psig is distributed in the standard form (coarse and fine arrays)
subroutine compress_and_accumulate_standard(grid,wfd,&
     keyv_c,keyv_f,keyg_c,keyg_f,&
     psig_c,psig_f,psi_c,psi_f)
  use module_base
  use module_types
  implicit none
  type(grid_dimensions), intent(in) :: grid
  type(wavefunctions_descriptors), intent(in) :: wfd
  integer, dimension(wfd%nseg_c), intent(in) :: keyv_c
  integer, dimension(wfd%nseg_f), intent(in) :: keyv_f
  integer, dimension(2,wfd%nseg_c), intent(in) :: keyg_c
  integer, dimension(2,wfd%nseg_f), intent(in) :: keyg_f
  real(wp), dimension(0:grid%n1,0:grid%n2,0:grid%n3), intent(in) :: psig_c
  real(wp), dimension(7,grid%nfl1:grid%nfu1,grid%nfl2:grid%nfu2,grid%nfl3:grid%nfu3), intent(in) :: psig_f
  real(wp), dimension(wfd%nvctr_c), intent(inout) :: psi_c
  real(wp), dimension(7,wfd%nvctr_f), intent(inout) :: psi_f
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

  !$omp parallel default(shared) &
  !$omp private(iseg,jj,j0,j1,ii,i1,i2,i3,i0,i)

  ! coarse part
  !$omp do
  do iseg=1,wfd%nseg_c
     jj=keyv_c(iseg)
     j0=keyg_c(1,iseg)
     j1=keyg_c(2,iseg)
     ii=j0-1
     i3=ii/((grid%n1+1)*(grid%n2+1))
     ii=ii-i3*(grid%n1+1)*(grid%n2+1)
     i2=ii/(grid%n1+1)
     i0=ii-i2*(grid%n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psi_c(i-i0+jj)=psi_c(i-i0+jj)+psig_c(i,i2,i3)
     enddo
  enddo
  !$omp enddo
  ! fine part
  !$omp do
  do iseg=1,wfd%nseg_f
     jj=keyv_f(iseg)
     j0=keyg_f(1,iseg)
     j1=keyg_f(2,iseg)
     ii=j0-1
     i3=ii/((grid%n1+1)*(grid%n2+1))
     ii=ii-i3*(grid%n1+1)*(grid%n2+1)
     i2=ii/(grid%n1+1)
     i0=ii-i2*(grid%n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psi_f(1,i-i0+jj)=psi_f(1,i-i0+jj)+psig_f(1,i,i2,i3)
        psi_f(2,i-i0+jj)=psi_f(2,i-i0+jj)+psig_f(2,i,i2,i3)
        psi_f(3,i-i0+jj)=psi_f(3,i-i0+jj)+psig_f(3,i,i2,i3)
        psi_f(4,i-i0+jj)=psi_f(4,i-i0+jj)+psig_f(4,i,i2,i3)
        psi_f(5,i-i0+jj)=psi_f(5,i-i0+jj)+psig_f(5,i,i2,i3)
        psi_f(6,i-i0+jj)=psi_f(6,i-i0+jj)+psig_f(6,i,i2,i3)
        psi_f(7,i-i0+jj)=psi_f(7,i-i0+jj)+psig_f(7,i,i2,i3)
     enddo
  enddo
  !$omp enddo
  !$omp end parallel

end subroutine compress_and_accumulate_standard

!> Compress the wavefunction psig and accumulate the result on the psi array
!! The wavefunction psig is distributed in the mixed form (fisrt coarse then fine components in each direction)
subroutine compress_and_accumulate_mixed(grid,wfd,&
     keyv_c,keyv_f,keyg_c,keyg_f,&
     psig,psi_c,psi_f)
  use module_base
  use module_types
  implicit none
  type(grid_dimensions), intent(in) :: grid
  type(wavefunctions_descriptors), intent(in) :: wfd
  integer, dimension(wfd%nseg_c), intent(in) :: keyv_c
  integer, dimension(wfd%nseg_f), intent(in) :: keyv_f
  integer, dimension(2,wfd%nseg_c), intent(in) :: keyg_c
  integer, dimension(2,wfd%nseg_f), intent(in) :: keyg_f
  real(wp), dimension(0:grid%n1,2,0:grid%n2,2,0:grid%n3,2), intent(in) :: psig
  real(wp), dimension(wfd%nvctr_c), intent(inout) :: psi_c
  real(wp), dimension(7,wfd%nvctr_f), intent(inout) :: psi_f
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

  !$omp parallel default(shared) &
  !$omp private(iseg,jj,j0,j1,ii,i1,i2,i3,i0,i)
  
  ! coarse part
  !$omp do
  do iseg=1,wfd%nseg_c
     jj=keyv_c(iseg)
     j0=keyg_c(1,iseg)
     j1=keyg_c(2,iseg)
     ii=j0-1
     i3=ii/((grid%n1+1)*(grid%n2+1))
     ii=ii-i3*(grid%n1+1)*(grid%n2+1)
     i2=ii/(grid%n1+1)
     i0=ii-i2*(grid%n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psi_c(i-i0+jj)=psi_c(i-i0+jj)+psig(i,1,i2,1,i3,1)
     enddo
  enddo
  !$omp enddo
  ! fine part
  !$omp do
  do iseg=1,wfd%nseg_f
     jj=keyv_f(iseg)
     j0=keyg_f(1,iseg)
     j1=keyg_f(2,iseg)
     ii=j0-1
     i3=ii/((grid%n1+1)*(grid%n2+1))
     ii=ii-i3*(grid%n1+1)*(grid%n2+1)
     i2=ii/(grid%n1+1)
     i0=ii-i2*(grid%n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psi_f(1,i-i0+jj)=psi_f(1,i-i0+jj)+psig(i,2,i2,1,i3,1)
        psi_f(2,i-i0+jj)=psi_f(2,i-i0+jj)+psig(i,1,i2,2,i3,1)
        psi_f(3,i-i0+jj)=psi_f(3,i-i0+jj)+psig(i,2,i2,2,i3,1)
        psi_f(4,i-i0+jj)=psi_f(4,i-i0+jj)+psig(i,1,i2,1,i3,2)
        psi_f(5,i-i0+jj)=psi_f(5,i-i0+jj)+psig(i,2,i2,1,i3,2)
        psi_f(6,i-i0+jj)=psi_f(6,i-i0+jj)+psig(i,1,i2,2,i3,2)
        psi_f(7,i-i0+jj)=psi_f(7,i-i0+jj)+psig(i,2,i2,2,i3,2)
     enddo
  enddo
  !$omp enddo
  !$omp end parallel
  
end subroutine compress_and_accumulate_mixed
