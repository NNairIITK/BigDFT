!> @file
!!    Wrapper for simplifying the call
!! @author
!!    Copyright (C) 2010-2011 BigDFT group (LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

!> Wrapper
subroutine wnrm_wrap(ncplx,mvctr_c,mvctr_f,psi,scpr)
  use module_base
  implicit none
  integer, intent(in) :: mvctr_c,mvctr_f,ncplx
  real(wp), dimension((mvctr_c+7*mvctr_f)*ncplx), intent(in) :: psi
  real(dp), intent(out) :: scpr
  !local variables
  integer :: i_f
  real(dp) :: scalp

  i_f=min(mvctr_f,1)
 
  call wnrm(mvctr_c,mvctr_f,psi,psi(mvctr_c+i_f),scpr)

  if (ncplx ==2) then
     call wnrm(mvctr_c,mvctr_f,&
          psi(mvctr_c+7*mvctr_f+1),psi(mvctr_c+7*mvctr_f+mvctr_c+i_f),scalp)
     scpr=scpr+scalp
  end if
  
END SUBROUTINE wnrm_wrap


!> Calculates the norm SQUARED (scpr) of a wavefunction (in vector form)
subroutine wnrm(mvctr_c,mvctr_f,psi_c,psi_f,scpr)
  use module_base
  implicit none
  !Arguments
  integer, intent(in) :: mvctr_c,mvctr_f
  real(wp), dimension(mvctr_c), intent(in) :: psi_c
  real(wp), dimension(7,mvctr_f), intent(in) :: psi_f
  real(dp), intent(out) :: scpr
  !local variables
  integer :: i
  real(dp) :: scpr0,scpr1,scpr2,scpr3,scpr4,scpr5,scpr6,scpr7
!!!    integer :: ncount0,ncount2,ncount_rate,ncount_max
!!!    real(gp) :: tel

!!!  !dee
!!!    open(unit=97,file='time_wnrm',status='unknown',position='append')
!!!    call system_clock(ncount0,ncount_rate,ncount_max)

    scpr=0.0_dp
    scpr0=0.0_dp
    scpr1=0.0_dp
    scpr2=0.0_dp
    scpr3=0.0_dp
    scpr4=0.0_dp
    scpr5=0.0_dp
    scpr6=0.0_dp
    scpr7=0.0_dp

!$omp parallel default(private)&
!$omp shared(mvctr_c,mvctr_f,psi_c,psi_f,scpr) firstprivate(scpr0,scpr1,scpr2,scpr3,scpr4,scpr5,scpr6,scpr7)

!$omp do schedule(static)
 do i=1,mvctr_c
    scpr0=scpr0+real(psi_c(i),dp)**2
 enddo
!$omp end do nowait
!$omp do schedule(guided)
 do i=1,mvctr_f
    scpr1=scpr1+real(psi_f(1,i),dp)**2
    scpr2=scpr2+real(psi_f(2,i),dp)**2
    scpr3=scpr3+real(psi_f(3,i),dp)**2
    scpr4=scpr4+real(psi_f(4,i),dp)**2
    scpr5=scpr5+real(psi_f(5,i),dp)**2
    scpr6=scpr6+real(psi_f(6,i),dp)**2
    scpr7=scpr7+real(psi_f(7,i),dp)**2
enddo
!$omp end do
    scpr0=scpr0+scpr1+scpr2+scpr3+scpr4+scpr5+scpr6+scpr7

!$omp critical 
    scpr=scpr+scpr0
!$omp end critical

!$omp end parallel

!!!    call system_clock(ncount2,ncount_rate,ncount_max)
!!!    tel=dble(ncount2-ncount0)/dble(ncount_rate)
!!!    write(97,*) 'wnrm:',tel
!!!    close(97)

END SUBROUTINE wnrm


!> Wrapper for simplifying the call
subroutine wscal_wrap(mvctr_c,mvctr_f,scal,psi)
  use module_base
  implicit none
  integer, intent(in) :: mvctr_c,mvctr_f
  real(wp), intent(in) :: scal
  real(wp), dimension(mvctr_c+7*mvctr_f), intent(inout) :: psi
  !local variables
  integer :: i_f

  i_f=min(mvctr_f,1)
 
  call wscal(mvctr_c,mvctr_f,scal,psi,psi(mvctr_c+i_f))
  
END SUBROUTINE wscal_wrap


!> Multiplies a wavefunction psi_c,psi_f (in vector form) with a scalar (scal)
subroutine wscal(mvctr_c,mvctr_f,scal,psi_c,psi_f)
  use module_base
  implicit none
  integer, intent(in) :: mvctr_c,mvctr_f
  real(wp), intent(in) :: scal
  real(wp), dimension(mvctr_c), intent(inout) :: psi_c
  real(wp), dimension(7,mvctr_f), intent(inout) :: psi_f
  !local variables
  integer :: i
!!!$omp parallel default(private) shared(mvctr_c,mvctr_f,scal,psi_c,psi_f)
!!!$omp do
  do i=1,mvctr_c
     psi_c(i)=psi_c(i)*scal
  enddo
!!!$omp enddo
!!!$omp do
  do i=1,mvctr_f
     psi_f(1,i)=psi_f(1,i)*scal
     psi_f(2,i)=psi_f(2,i)*scal
     psi_f(3,i)=psi_f(3,i)*scal
     psi_f(4,i)=psi_f(4,i)*scal
     psi_f(5,i)=psi_f(5,i)*scal
     psi_f(6,i)=psi_f(6,i)*scal
     psi_f(7,i)=psi_f(7,i)*scal
  enddo
!!!$omp enddo
!!!$omp end parallel

END SUBROUTINE wscal


!> Wrapper for simplifying the call
subroutine wscalv_wrap(mvctr_c,mvctr_f,scal,psi)
  use module_base
  implicit none
  integer, intent(in) :: mvctr_c,mvctr_f
  real(wp), dimension(0:3), intent(in) :: scal
  real(wp), dimension(mvctr_c+7*mvctr_f), intent(inout) :: psi
  !local variables
  integer :: i_f

  i_f=min(mvctr_f,1)
 
  call wscalv(mvctr_c,mvctr_f,scal,psi,psi(mvctr_c+i_f))
  
END SUBROUTINE wscalv_wrap


!> Multiplies a wavefunction psi_c,psi_f (in vector form) with a scaling vector (scal)
subroutine wscalv(mvctr_c,mvctr_f,scal,psi_c,psi_f)
  use module_base
  implicit none
  integer, intent(in) :: mvctr_c,mvctr_f
  real(wp), dimension(0:3), intent(in) :: scal
  real(wp), dimension(mvctr_c), intent(inout) :: psi_c
  real(wp), dimension(7,mvctr_f), intent(inout) :: psi_f
  !local variables
  integer :: i
!!!$omp parallel default(private) shared(mvctr_c,mvctr_f,scal,psi_c,psi_f)
!!!$omp do
  do i=1,mvctr_c
     psi_c(i)=psi_c(i)*scal(0)           !  1 1 1
  enddo
!!!$omp enddo
!!!$omp do
  do i=1,mvctr_f
     psi_f(1,i)=psi_f(1,i)*scal(1)       !  2 1 1
     psi_f(2,i)=psi_f(2,i)*scal(1)       !  1 2 1
     psi_f(3,i)=psi_f(3,i)*scal(2)       !  2 2 1
     psi_f(4,i)=psi_f(4,i)*scal(1)       !  1 1 2
     psi_f(5,i)=psi_f(5,i)*scal(2)       !  2 1 2
     psi_f(6,i)=psi_f(6,i)*scal(2)       !  1 2 2
     psi_f(7,i)=psi_f(7,i)*scal(3)       !  2 2 2
  enddo
!!!$omp enddo
!!!$omp end parallel
END SUBROUTINE wscalv


!> Initializes a wavefunction to zero
subroutine wzero(mvctr_c,mvctr_f,psi_c,psi_f)
  use module_base
  implicit none
  integer, intent(in) :: mvctr_c,mvctr_f
  real(wp), dimension(mvctr_c), intent(out) :: psi_c
  real(wp), dimension(7,mvctr_f), intent(out) :: psi_f
  !local variables
  integer :: i

! seems to be never called
!write(*,*) ' i am in wzero'
!!!$omp parallel default(private) shared(mvctr_c,mvctr_f,psi_c,psi_f)
!!!$omp do
  do i=1,mvctr_c
     psi_c(i)=0.0_wp
  enddo
!!!$omp enddo
!!!$omp do
  do i=1,mvctr_f
     psi_f(1,i)=0.0_wp
     psi_f(2,i)=0.0_wp
     psi_f(3,i)=0.0_wp
     psi_f(4,i)=0.0_wp
     psi_f(5,i)=0.0_wp
     psi_f(6,i)=0.0_wp
     psi_f(7,i)=0.0_wp
  enddo
!!!$omp enddo
!!!$omp end parallel
END SUBROUTINE wzero


!> Wrapper of wpdot1 to avoid boundary problems in absence of wavelets
!! and to perform scalar product for complex wavefunctions and projectors
!! if the wavefunctions are complex, so should be also the projectors
subroutine wpdot_wrap1(ncplx,mavctr_c,mavctr_f,maseg_c,maseg_f,keyav,keyag,apsi,  &
     mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keybv,keybg,bpsi,scpr,proj_count)
  use module_base
  implicit none
  integer, intent(in) :: mavctr_c,mavctr_f,maseg_c,maseg_f,proj_count
  integer, intent(in) :: ncplx,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f
  integer, dimension(maseg_c+maseg_f), intent(in) :: keyav
  integer, dimension(mbseg_c+mbseg_f), intent(in) :: keybv
  integer, dimension(2,maseg_c+maseg_f), intent(in) :: keyag
  integer, dimension(2,mbseg_c+mbseg_f), intent(in) :: keybg
  real(wp), dimension(mavctr_c+7*mavctr_f,ncplx), intent(in) :: apsi
  real(wp), dimension(mbvctr_c+7*mbvctr_f,ncplx,proj_count), intent(in) :: bpsi
  real(dp), dimension(proj_count,ncplx), intent(out) :: scpr
  !local variables
  integer :: ia_f,ib_f,iaseg_f,ibseg_f,ia,ib
  real(dp), dimension(proj_count,ncplx,ncplx) :: scalprod 

  ia_f=min(mavctr_f,1)
  ib_f=min(mbvctr_f,1)

  iaseg_f=min(maseg_f,1)
  ibseg_f=min(mbseg_f,1)

  do ia=1,ncplx
     do ib=1,ncplx
        select case (proj_count)

        case(4)
        call wpdot_4(mavctr_c,mavctr_f,maseg_c,maseg_f,&
             keyav,keyav(maseg_c+iaseg_f),&
             keyag,keyag(1,maseg_c+iaseg_f),&
             apsi(1,ia),apsi(mavctr_c+ia_f,ia),  &
             mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
             keybv,keybv(mbseg_c+ibseg_f),&
             keybg,keybg(1,mbseg_c+ibseg_f),&
             bpsi(1,ib,1),bpsi(1,ib,2),bpsi(1,ib,3),bpsi(1,ib,4), &
             bpsi(mbvctr_c+ib_f,ib,1),bpsi(mbvctr_c+ib_f,ib,2),bpsi(mbvctr_c+ib_f,ib,3), &
             bpsi(mbvctr_c+ib_f,ib,4), &
             scalprod(1,ia,ib),&
             proj_count)

        case(5)
        call wpdot_5(mavctr_c,mavctr_f,maseg_c,maseg_f,&
             keyav,keyav(maseg_c+iaseg_f),&
             keyag,keyag(1,maseg_c+iaseg_f),&
             apsi(1,ia),apsi(mavctr_c+ia_f,ia),  &
             mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
             keybv,keybv(mbseg_c+ibseg_f),&
             keybg,keybg(1,mbseg_c+ibseg_f),&
             bpsi(1,ib,1),bpsi(1,ib,2),bpsi(1,ib,3),bpsi(1,ib,4),bpsi(1,ib,5), &
             bpsi(mbvctr_c+ib_f,ib,1),bpsi(mbvctr_c+ib_f,ib,2),bpsi(mbvctr_c+ib_f,ib,3), &
             bpsi(mbvctr_c+ib_f,ib,4),bpsi(mbvctr_c+ib_f,ib,5), &
             scalprod(1,ia,ib),&
             proj_count)

        case(8)
        call wpdot_8(mavctr_c,mavctr_f,maseg_c,maseg_f,&
             keyav,keyav(maseg_c+iaseg_f),&
             keyag,keyag(1,maseg_c+iaseg_f),&
             apsi(1,ia),apsi(mavctr_c+ia_f,ia),  &
             mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
             keybv,keybv(mbseg_c+ibseg_f),&
             keybg,keybg(1,mbseg_c+ibseg_f),&
             bpsi(1,ib,1),bpsi(1,ib,2),bpsi(1,ib,3),bpsi(1,ib,4),bpsi(1,ib,5),bpsi(1,ib,6),bpsi(1,ib,7), &
             bpsi(1,ib,8), &
             bpsi(mbvctr_c+ib_f,ib,1),bpsi(mbvctr_c+ib_f,ib,2),bpsi(mbvctr_c+ib_f,ib,3), &
             bpsi(mbvctr_c+ib_f,ib,4),bpsi(mbvctr_c+ib_f,ib,5),bpsi(mbvctr_c+ib_f,ib,6), &
             bpsi(mbvctr_c+ib_f,ib,7),bpsi(mbvctr_c+ib_f,ib,8), &
             scalprod(1,ia,ib),&
             proj_count)

        case(13)
        call wpdot_13(mavctr_c,mavctr_f,maseg_c,maseg_f,&
             keyav,keyav(maseg_c+iaseg_f),&
             keyag,keyag(1,maseg_c+iaseg_f),&
             apsi(1,ia),apsi(mavctr_c+ia_f,ia),  &
             mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
             keybv,keybv(mbseg_c+ibseg_f),&
             keybg,keybg(1,mbseg_c+ibseg_f),&
             bpsi(1,ib,1),bpsi(1,ib,2),bpsi(1,ib,3),bpsi(1,ib,4),bpsi(1,ib,5),bpsi(1,ib,6),bpsi(1,ib,7), &
             bpsi(1,ib,8),bpsi(1,ib,9),bpsi(1,ib,10),bpsi(1,ib,11),bpsi(1,ib,12),bpsi(1,ib,13), &
             bpsi(mbvctr_c+ib_f,ib,1),bpsi(mbvctr_c+ib_f,ib,2),bpsi(mbvctr_c+ib_f,ib,3), &
             bpsi(mbvctr_c+ib_f,ib,4),bpsi(mbvctr_c+ib_f,ib,5),bpsi(mbvctr_c+ib_f,ib,6), &
             bpsi(mbvctr_c+ib_f,ib,7),bpsi(mbvctr_c+ib_f,ib,8),bpsi(mbvctr_c+ib_f,ib,9), &
             bpsi(mbvctr_c+ib_f,ib,10),bpsi(mbvctr_c+ib_f,ib,11),bpsi(mbvctr_c+ib_f,ib,12), &
             bpsi(mbvctr_c+ib_f,ib,13), &
             scalprod(1,ia,ib),&
             proj_count)

        case(14)
        call wpdot_14(mavctr_c,mavctr_f,maseg_c,maseg_f,&
             keyav,keyav(maseg_c+iaseg_f),&
             keyag,keyag(1,maseg_c+iaseg_f),&
             apsi(1,ia),apsi(mavctr_c+ia_f,ia),  &
             mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
             keybv,keybv(mbseg_c+ibseg_f),&
             keybg,keybg(1,mbseg_c+ibseg_f),&
             bpsi(1,ib,1),bpsi(1,ib,2),bpsi(1,ib,3),bpsi(1,ib,4),bpsi(1,ib,5),bpsi(1,ib,6),bpsi(1,ib,7),&
             bpsi(1,ib,8),bpsi(1,ib,9),bpsi(1,ib,10),bpsi(1,ib,11),bpsi(1,ib,12),bpsi(1,ib,13),bpsi(1,ib,14),&
             bpsi(mbvctr_c+ib_f,ib,1),bpsi(mbvctr_c+ib_f,ib,2),bpsi(mbvctr_c+ib_f,ib,3), &
             bpsi(mbvctr_c+ib_f,ib,4),bpsi(mbvctr_c+ib_f,ib,5),bpsi(mbvctr_c+ib_f,ib,6), &
             bpsi(mbvctr_c+ib_f,ib,7),bpsi(mbvctr_c+ib_f,ib,8),bpsi(mbvctr_c+ib_f,ib,9), &
             bpsi(mbvctr_c+ib_f,ib,10),bpsi(mbvctr_c+ib_f,ib,11),bpsi(mbvctr_c+ib_f,ib,12), &
             bpsi(mbvctr_c+ib_f,ib,13),bpsi(mbvctr_c+ib_f,ib,14), &
             scalprod(1,ia,ib),&
             proj_count)

        case(18)
        call wpdot_18(mavctr_c,mavctr_f,maseg_c,maseg_f,&
             keyav,keyav(maseg_c+iaseg_f),&
             keyag,keyag(1,maseg_c+iaseg_f),&
             apsi(1,ia),apsi(mavctr_c+ia_f,ia),  &
             mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
             keybv,keybv(mbseg_c+ibseg_f),&
             keybg,keybg(1,mbseg_c+ibseg_f),&
             bpsi(1,ib,1),bpsi(1,ib,2),bpsi(1,ib,3),bpsi(1,ib,4),bpsi(1,ib,5),bpsi(1,ib,6),bpsi(1,ib,7),&
             bpsi(1,ib,8),bpsi(1,ib,9),bpsi(1,ib,10),bpsi(1,ib,11),bpsi(1,ib,12),bpsi(1,ib,13),bpsi(1,ib,14),&
             bpsi(1,ib,15),bpsi(1,ib,16),bpsi(1,ib,17),bpsi(1,ib,18),&
             bpsi(mbvctr_c+ib_f,ib,1),bpsi(mbvctr_c+ib_f,ib,2),bpsi(mbvctr_c+ib_f,ib,3), &
             bpsi(mbvctr_c+ib_f,ib,4),bpsi(mbvctr_c+ib_f,ib,5),bpsi(mbvctr_c+ib_f,ib,6), &
             bpsi(mbvctr_c+ib_f,ib,7),bpsi(mbvctr_c+ib_f,ib,8),bpsi(mbvctr_c+ib_f,ib,9), &
             bpsi(mbvctr_c+ib_f,ib,10),bpsi(mbvctr_c+ib_f,ib,11),bpsi(mbvctr_c+ib_f,ib,12), &
             bpsi(mbvctr_c+ib_f,ib,13),bpsi(mbvctr_c+ib_f,ib,14),bpsi(mbvctr_c+ib_f,ib,15), &
             bpsi(mbvctr_c+ib_f,ib,16),bpsi(mbvctr_c+ib_f,ib,17),bpsi(mbvctr_c+ib_f,ib,18), &
             scalprod(1,ia,ib),&
             proj_count)

        case(19)
        call wpdot_19(mavctr_c,mavctr_f,maseg_c,maseg_f,&
             keyav,keyav(maseg_c+iaseg_f),&
             keyag,keyag(1,maseg_c+iaseg_f),&
             apsi(1,ia),apsi(mavctr_c+ia_f,ia),  &
             mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
             keybv,keybv(mbseg_c+ibseg_f),&
             keybg,keybg(1,mbseg_c+ibseg_f),&
             bpsi(1,ib,1),bpsi(1,ib,2),bpsi(1,ib,3),bpsi(1,ib,4),bpsi(1,ib,5),bpsi(1,ib,6),bpsi(1,ib,7),&
             bpsi(1,ib,8),bpsi(1,ib,9),bpsi(1,ib,10),bpsi(1,ib,11),bpsi(1,ib,12),bpsi(1,ib,13),bpsi(1,ib,14),&
             bpsi(1,ib,15),bpsi(1,ib,16),bpsi(1,ib,17),bpsi(1,ib,18),bpsi(1,ib,19),&
             bpsi(mbvctr_c+ib_f,ib,1),bpsi(mbvctr_c+ib_f,ib,2),bpsi(mbvctr_c+ib_f,ib,3), &
             bpsi(mbvctr_c+ib_f,ib,4),bpsi(mbvctr_c+ib_f,ib,5),bpsi(mbvctr_c+ib_f,ib,6), &
             bpsi(mbvctr_c+ib_f,ib,7),bpsi(mbvctr_c+ib_f,ib,8),bpsi(mbvctr_c+ib_f,ib,9), &
             bpsi(mbvctr_c+ib_f,ib,10),bpsi(mbvctr_c+ib_f,ib,11),bpsi(mbvctr_c+ib_f,ib,12), &
             bpsi(mbvctr_c+ib_f,ib,13),bpsi(mbvctr_c+ib_f,ib,14),bpsi(mbvctr_c+ib_f,ib,15), &
             bpsi(mbvctr_c+ib_f,ib,16),bpsi(mbvctr_c+ib_f,ib,17),bpsi(mbvctr_c+ib_f,ib,18), &
             bpsi(mbvctr_c+ib_f,ib,19),scalprod(1,ia,ib),&
             proj_count)

        case(20)
        call wpdot_20(mavctr_c,mavctr_f,maseg_c,maseg_f,&
             keyav,keyav(maseg_c+iaseg_f),&
             keyag,keyag(1,maseg_c+iaseg_f),&
             apsi(1,ia),apsi(mavctr_c+ia_f,ia),  &
             mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
             keybv,keybv(mbseg_c+ibseg_f),&
             keybg,keybg(1,mbseg_c+ibseg_f),&
             bpsi(1,ib,1),bpsi(1,ib,2),bpsi(1,ib,3),bpsi(1,ib,4),bpsi(1,ib,5),bpsi(1,ib,6),bpsi(1,ib,7),&
             bpsi(1,ib,8),bpsi(1,ib,9),bpsi(1,ib,10),bpsi(1,ib,11),bpsi(1,ib,12),bpsi(1,ib,13),bpsi(1,ib,14),&
             bpsi(1,ib,15),bpsi(1,ib,16),bpsi(1,ib,17),bpsi(1,ib,18),bpsi(1,ib,19),bpsi(1,ib,20),&
             bpsi(mbvctr_c+ib_f,ib,1),bpsi(mbvctr_c+ib_f,ib,2),bpsi(mbvctr_c+ib_f,ib,3), &
             bpsi(mbvctr_c+ib_f,ib,4),bpsi(mbvctr_c+ib_f,ib,5),bpsi(mbvctr_c+ib_f,ib,6), &
             bpsi(mbvctr_c+ib_f,ib,7),bpsi(mbvctr_c+ib_f,ib,8),bpsi(mbvctr_c+ib_f,ib,9), &
             bpsi(mbvctr_c+ib_f,ib,10),bpsi(mbvctr_c+ib_f,ib,11),bpsi(mbvctr_c+ib_f,ib,12), &
             bpsi(mbvctr_c+ib_f,ib,13),bpsi(mbvctr_c+ib_f,ib,14),bpsi(mbvctr_c+ib_f,ib,15), &
             bpsi(mbvctr_c+ib_f,ib,16),bpsi(mbvctr_c+ib_f,ib,17),bpsi(mbvctr_c+ib_f,ib,18), &
             bpsi(mbvctr_c+ib_f,ib,19),bpsi(mbvctr_c+ib_f,ib,20),scalprod(1,ia,ib),&
             proj_count)
 
        case(22)
        call wpdot_22(mavctr_c,mavctr_f,maseg_c,maseg_f,&
             keyav,keyav(maseg_c+iaseg_f),&
             keyag,keyag(1,maseg_c+iaseg_f),&
             apsi(1,ia),apsi(mavctr_c+ia_f,ia),  &
             mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
             keybv,keybv(mbseg_c+ibseg_f),&
             keybg,keybg(1,mbseg_c+ibseg_f),&
             bpsi(1,ib,1),bpsi(1,ib,2),bpsi(1,ib,3),bpsi(1,ib,4),bpsi(1,ib,5),bpsi(1,ib,6),bpsi(1,ib,7),&
             bpsi(1,ib,8),bpsi(1,ib,9),bpsi(1,ib,10),bpsi(1,ib,11),bpsi(1,ib,12),bpsi(1,ib,13),bpsi(1,ib,14),&
             bpsi(1,ib,15),bpsi(1,ib,16),bpsi(1,ib,17),bpsi(1,ib,18),bpsi(1,ib,19),bpsi(1,ib,20),&
             bpsi(1,ib,21),bpsi(1,ib,22),&
             bpsi(mbvctr_c+ib_f,ib,1),bpsi(mbvctr_c+ib_f,ib,2),bpsi(mbvctr_c+ib_f,ib,3), &
             bpsi(mbvctr_c+ib_f,ib,4),bpsi(mbvctr_c+ib_f,ib,5),bpsi(mbvctr_c+ib_f,ib,6), &
             bpsi(mbvctr_c+ib_f,ib,7),bpsi(mbvctr_c+ib_f,ib,8),bpsi(mbvctr_c+ib_f,ib,9), &
             bpsi(mbvctr_c+ib_f,ib,10),bpsi(mbvctr_c+ib_f,ib,11),bpsi(mbvctr_c+ib_f,ib,12), &
             bpsi(mbvctr_c+ib_f,ib,13),bpsi(mbvctr_c+ib_f,ib,14),bpsi(mbvctr_c+ib_f,ib,15), &
             bpsi(mbvctr_c+ib_f,ib,16),bpsi(mbvctr_c+ib_f,ib,17),bpsi(mbvctr_c+ib_f,ib,18), &
             bpsi(mbvctr_c+ib_f,ib,19),bpsi(mbvctr_c+ib_f,ib,20), &
             bpsi(mbvctr_c+ib_f,ib,21),bpsi(mbvctr_c+ib_f,ib,22),scalprod(1,ia,ib),&
             proj_count)

         end select
     end do
  end do


  !then define the result
  if (ncplx == 1) then
     scpr(:,1)=scalprod(:,1,1)
  else if (ncplx == 2) then
     scpr(:,1)=scalprod(:,1,1)+scalprod(:,2,2)
     scpr(:,2)=scalprod(:,1,2)-scalprod(:,2,1)
  else
     write(*,*)'ERROR wpdot: ncplx not valid:',ncplx
     stop
  end if

END SUBROUTINE wpdot_wrap1

!> This function must be generalized for the linear scaling code                
!! @warning   
!! calculates the dot product between a wavefunctions apsi and a projector bpsi (both in compressed form)
!! Warning: It is assumed that the segments of bpsi are always contained inside
!! the segments of apsi.
!! This version multiplies apsi with all the projectors at once in order to increase efficiency.

subroutine wpdot_4(  &
     mavctr_c,mavctr_f,maseg_c,maseg_f,keyav_c,keyav_f,keyag_c,keyag_f,apsi_c,apsi_f,  &
     mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keybv_c,keybv_f,keybg_c,keybg_f, &
     bpsi1,bpsi2,bpsi3,bpsi4,bpsi_f1,bpsi_f2,bpsi_f3,bpsi_f4, &
     scpr,proj_count)

  use module_base
  implicit none
  integer, intent(in) :: mavctr_c,mavctr_f,maseg_c,maseg_f,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f
  integer, intent(in) :: proj_count
  integer, dimension(maseg_c), intent(in) :: keyav_c
  integer, dimension(mbseg_c), intent(in) :: keybv_c
  integer, dimension(maseg_f), intent(in) :: keyav_f
  integer, dimension(mbseg_f), intent(in) :: keybv_f
  integer, dimension(2,maseg_c), intent(in) :: keyag_c
  integer, dimension(2,mbseg_c), intent(in) :: keybg_c
  integer, dimension(2,maseg_f), intent(in) :: keyag_f
  integer, dimension(2,mbseg_f), intent(in) :: keybg_f
  real(wp), dimension(mavctr_c), intent(in) :: apsi_c
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi1,bpsi2,bpsi3,bpsi4
  real(wp), dimension(7,mavctr_f), intent(in) :: apsi_f
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f1,bpsi_f2,bpsi_f3,bpsi_f4
  real(dp), dimension(proj_count),intent(out) :: scpr
  !local variables
  integer :: ibseg,jaj,jb1,jb0,jbj,iaoff,length,i,ja0,ja1
  real(dp),dimension(proj_count) :: scpr1,scpr2,scpr3,scpr4,scpr5,scpr6,scpr7,scpr0
  integer :: iaseg0
  integer, dimension(maseg_c) :: keyag_c_lin !>linear version of second indices of keyag_c
  integer, dimension(maseg_f) :: keyag_f_lin !>linear version of second indices of keyag_f
  real(dp) :: apsi_temp
  real(dp), dimension(7) :: apsi_f_temp
!  !$ integer :: omp_get_max_threads,omp_get_thread_num,omp_get_num_threads
!!!    integer :: ncount0,ncount2,ncount_rate,ncount_max
!!!    real(gp) :: tel

!!!  !dee
!!!    open(unit=97,file='time_wpdot',status='unknown',position='append')
!!!    call system_clock(ncount0,ncount_rate,ncount_max)

  keyag_c_lin = keyag_c(1,:)!speed up access in hunt subroutine by consecutive arrangement in memory
  keyag_f_lin = keyag_f(1,:)!speed up access in hunt subroutine by consecutive arrangement in memory

  scpr=0.0_dp

!$omp parallel default(private) shared(maseg_c,keyav_c,keyag_c,keyag_c_lin) &
!$omp shared(keybg_c,mbseg_c,keybv_c,mbseg_f,maseg_f,apsi_c,bpsi1,bpsi2) &
!$omp shared(bpsi3,bpsi4,bpsi_f1,bpsi_f2,bpsi_f3,bpsi_f4) &
!$omp shared(keybv_f,keybg_f,keyag_f,keyag_f_lin,keyav_f) &
!$omp shared(apsi_f,scpr)

!  !$  print *,'AAANonLocalHamiltonian with nthread:, out to:' ,omp_get_num_threads(),omp_get_thread_num()
!!!!$omp shared (ncount0,ncount2,ncount_rate,ncount_max,tel)

  iaseg0=1 

  scpr0=0.0_dp
  scpr1=0.0_dp
  scpr2=0.0_dp
  scpr3=0.0_dp
  scpr4=0.0_dp
  scpr5=0.0_dp
  scpr6=0.0_dp
  scpr7=0.0_dp
!coarse part. Loop on the projectors segments
!$omp do schedule(static)
   do ibseg=1,mbseg_c
     jbj=keybv_c(ibseg)
     jb0=keybg_c(1,ibseg) !starting point of projector segment
     jb1=keybg_c(2,ibseg) !ending point of projector segment
!     print *,'huntenter',ibseg,jb0,jb1
 
     !find the starting point of the wavefunction segment
     !warning: hunt is assuming that the variable is always found
     !if it is not, iaseg0 is put to maseg + 1 so that the loop is disabled
     call hunt1(.true.,keyag_c_lin,maseg_c,keybg_c(1,ibseg),iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     !now pass through all the wavefunction segments until the end of the segment is 
     !still contained in projector segment
     nonconvex_loop_c: do while(iaseg0 <= maseg_c)
!     print *,'huntexit',iaseg0,maseg_c,keyag_c_lin(iaseg0),keyag_c(2,iaseg0)

        !length = jb1-jb0
        !iaoff = jb0-keyag_c_lin(iaseg0)!jb0-ja0

        ja0=keyag_c_lin(iaseg0)
        ja1=min(jb1,keyag_c(2,iaseg0)) 
        length = ja1-jb0
        iaoff = max(jb0-ja0,0) !no offset if we are already inside

        jaj=keyav_c(iaseg0)
        do i=0,length
           apsi_temp=real(apsi_c(jaj+iaoff+i),dp)

           scpr0(1)=scpr0(1)+apsi_temp *real(bpsi1(jbj+i),dp)
           scpr0(2)=scpr0(2)+apsi_temp *real(bpsi2(jbj+i),dp)
           scpr0(3)=scpr0(3)+apsi_temp *real(bpsi3(jbj+i),dp)
           scpr0(4)=scpr0(4)+apsi_temp *real(bpsi4(jbj+i),dp)
        enddo
 !       print *,'length',length,ibseg,scpr0,iaseg0,ja1,jb1

        !print *,'ibseg,mbseg_c,iaseg0,maseg_c',ibseg,mbseg_c,iaseg0,maseg_c
        !print '(a,6(i8),1pe25.17)','ja0,ja1t,ja1,jb0,jb1',&
        !     ibseg,ja0,keyag_c(2,iaseg0),ja1,jb0,jb1,scpr0
        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_c) exit nonconvex_loop_c !segment is finished
        iaseg0=iaseg0+1
        jb0=max(keybg_c(1,ibseg),keyag_c_lin(iaseg0))
        if (keyag_c_lin(iaseg0)>jb1) exit nonconvex_loop_c !segment is not covered
        jbj=jbj+max(jb0-keybg_c(1,ibseg),0)
     end do nonconvex_loop_c
     !disable loop if the end is reached
     if (iaseg0 == maseg_c .and. keybg_c(1,ibseg)> keyag_c_lin(maseg_c)) iaseg0=iaseg0+1


   enddo
!stop
!$omp end do nowait

! fine part

iaseg0=1

!$omp do schedule(static)
   do ibseg=1,mbseg_f
     jbj=keybv_f(ibseg)
     jb0=keybg_f(1,ibseg)
     jb1=keybg_f(2,ibseg)
!    print *,'huntenter',ibseg,jb0,jb1
     call hunt1(.true.,keyag_f_lin,maseg_f,keybg_f(1,ibseg),iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     nonconvex_loop_f: do while(iaseg0 <= maseg_f)
!!$     length = jb1-jb0
!!$     iaoff = jb0-keyag_f_lin(iaseg0)

        ja0=keyag_f_lin(iaseg0) !still doubts about copying in automatic array
        ja1=min(jb1,keyag_f(2,iaseg0)) 
        length = ja1-jb0
        iaoff = max(jb0-ja0,0) !no offset if we are already inside

        jaj=keyav_f(iaseg0)
        do i=0,length
           apsi_f_temp(1)=real(apsi_f(1,jaj+iaoff+i),dp)
           apsi_f_temp(2)=real(apsi_f(2,jaj+iaoff+i),dp)
           apsi_f_temp(3)=real(apsi_f(3,jaj+iaoff+i),dp)
           apsi_f_temp(4)=real(apsi_f(4,jaj+iaoff+i),dp)
           apsi_f_temp(5)=real(apsi_f(5,jaj+iaoff+i),dp)
           apsi_f_temp(6)=real(apsi_f(6,jaj+iaoff+i),dp)
           apsi_f_temp(7)=real(apsi_f(7,jaj+iaoff+i),dp)

           scpr1(1)=scpr1(1)+apsi_f_temp(1)*real(bpsi_f1(1,jbj+i),dp)
           scpr2(1)=scpr2(1)+apsi_f_temp(2)*real(bpsi_f1(2,jbj+i),dp)
           scpr3(1)=scpr3(1)+apsi_f_temp(3)*real(bpsi_f1(3,jbj+i),dp)
           scpr4(1)=scpr4(1)+apsi_f_temp(4)*real(bpsi_f1(4,jbj+i),dp)
           scpr5(1)=scpr5(1)+apsi_f_temp(5)*real(bpsi_f1(5,jbj+i),dp)
           scpr6(1)=scpr6(1)+apsi_f_temp(6)*real(bpsi_f1(6,jbj+i),dp)
           scpr7(1)=scpr7(1)+apsi_f_temp(7)*real(bpsi_f1(7,jbj+i),dp)

           scpr1(2)=scpr1(2)+apsi_f_temp(1)*real(bpsi_f2(1,jbj+i),dp)
           scpr2(2)=scpr2(2)+apsi_f_temp(2)*real(bpsi_f2(2,jbj+i),dp)
           scpr3(2)=scpr3(2)+apsi_f_temp(3)*real(bpsi_f2(3,jbj+i),dp)
           scpr4(2)=scpr4(2)+apsi_f_temp(4)*real(bpsi_f2(4,jbj+i),dp)
           scpr5(2)=scpr5(2)+apsi_f_temp(5)*real(bpsi_f2(5,jbj+i),dp)
           scpr6(2)=scpr6(2)+apsi_f_temp(6)*real(bpsi_f2(6,jbj+i),dp)
           scpr7(2)=scpr7(2)+apsi_f_temp(7)*real(bpsi_f2(7,jbj+i),dp)

           scpr1(3)=scpr1(3)+apsi_f_temp(1)*real(bpsi_f3(1,jbj+i),dp)
           scpr2(3)=scpr2(3)+apsi_f_temp(2)*real(bpsi_f3(2,jbj+i),dp)
           scpr3(3)=scpr3(3)+apsi_f_temp(3)*real(bpsi_f3(3,jbj+i),dp)
           scpr4(3)=scpr4(3)+apsi_f_temp(4)*real(bpsi_f3(4,jbj+i),dp)
           scpr5(3)=scpr5(3)+apsi_f_temp(5)*real(bpsi_f3(5,jbj+i),dp)
           scpr6(3)=scpr6(3)+apsi_f_temp(6)*real(bpsi_f3(6,jbj+i),dp)
           scpr7(3)=scpr7(3)+apsi_f_temp(7)*real(bpsi_f3(7,jbj+i),dp)

           scpr1(4)=scpr1(4)+apsi_f_temp(1)*real(bpsi_f4(1,jbj+i),dp)
           scpr2(4)=scpr2(4)+apsi_f_temp(2)*real(bpsi_f4(2,jbj+i),dp)
           scpr3(4)=scpr3(4)+apsi_f_temp(3)*real(bpsi_f4(3,jbj+i),dp)
           scpr4(4)=scpr4(4)+apsi_f_temp(4)*real(bpsi_f4(4,jbj+i),dp)
           scpr5(4)=scpr5(4)+apsi_f_temp(5)*real(bpsi_f4(5,jbj+i),dp)
           scpr6(4)=scpr6(4)+apsi_f_temp(6)*real(bpsi_f4(6,jbj+i),dp)
           scpr7(4)=scpr7(4)+apsi_f_temp(7)*real(bpsi_f4(7,jbj+i),dp)

        enddo
 !       print *,'length',length,ibseg,scpr1,iaseg0,ja1,jb1
        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_f) exit nonconvex_loop_f !segment is finished  
        iaseg0=iaseg0+1
        jb0=max(keybg_f(1,ibseg),keyag_f_lin(iaseg0))
        if (keyag_f_lin(iaseg0)>jb1) exit nonconvex_loop_f !segment is not covered 
        jbj=jbj+max(jb0-keybg_f(1,ibseg),0)
     end do nonconvex_loop_f
     !disable loop if the end is reached
     if (iaseg0 == maseg_f .and. keybg_f(1,ibseg)> keyag_f_lin(maseg_f)) iaseg0=iaseg0+1

   enddo
!$omp end do !!!implicit barrier 

   scpr0=scpr0+scpr1+scpr2+scpr3+scpr4+scpr5+scpr6+scpr7

!$omp critical 
   scpr=scpr+scpr0
!$omp end critical
  
!$omp end parallel

!!!    call system_clock(ncount2,ncount_rate,ncount_max)
!!!    tel=dble(ncount2-ncount0)/dble(ncount_rate)
!!!    write(97,*) 'wpdot:',tel
!!!    close(97)

END SUBROUTINE wpdot_4

subroutine wpdot_5(  &
     mavctr_c,mavctr_f,maseg_c,maseg_f,keyav_c,keyav_f,keyag_c,keyag_f,apsi_c,apsi_f,  &
     mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keybv_c,keybv_f,keybg_c,keybg_f, &
     bpsi1,bpsi2,bpsi3,bpsi4,bpsi5,bpsi_f1,bpsi_f2,bpsi_f3,bpsi_f4,bpsi_f5, &
     scpr,proj_count)

  use module_base
  implicit none
  integer, intent(in) :: mavctr_c,mavctr_f,maseg_c,maseg_f,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f
  integer, intent(in) :: proj_count
  integer, dimension(maseg_c), intent(in) :: keyav_c
  integer, dimension(mbseg_c), intent(in) :: keybv_c
  integer, dimension(maseg_f), intent(in) :: keyav_f
  integer, dimension(mbseg_f), intent(in) :: keybv_f
  integer, dimension(2,maseg_c), intent(in) :: keyag_c
  integer, dimension(2,mbseg_c), intent(in) :: keybg_c
  integer, dimension(2,maseg_f), intent(in) :: keyag_f
  integer, dimension(2,mbseg_f), intent(in) :: keybg_f
  real(wp), dimension(mavctr_c), intent(in) :: apsi_c
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi1,bpsi2,bpsi3,bpsi4,bpsi5
  real(wp), dimension(7,mavctr_f), intent(in) :: apsi_f
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f1,bpsi_f2,bpsi_f3,bpsi_f4,bpsi_f5
  real(dp), dimension(proj_count),intent(out) :: scpr
  !local variables
  integer :: ibseg,jaj,jb1,jb0,jbj,iaoff,length,i,ja0,ja1
  real(dp),dimension(proj_count) :: scpr1,scpr2,scpr3,scpr4,scpr5,scpr6,scpr7,scpr0
  integer :: iaseg0
  integer, dimension(maseg_c) :: keyag_c_lin !>linear version of second indices of keyag_c
  integer, dimension(maseg_f) :: keyag_f_lin !>linear version of second indices of keyag_f
  real(dp) :: apsi_temp
  real(dp), dimension(7) :: apsi_f_temp
!!!    integer :: ncount0,ncount2,ncount_rate,ncount_max
!!!    real(gp) :: tel

!!!  !dee
!!!    open(unit=97,file='time_wpdot',status='unknown',position='append')
!!!    call system_clock(ncount0,ncount_rate,ncount_max)

  keyag_c_lin = keyag_c(1,:)!speed up access in hunt subroutine by consecutive arrangement in memory
  keyag_f_lin = keyag_f(1,:)!speed up access in hunt subroutine by consecutive arrangement in memory

  scpr=0.0_dp

!$omp parallel default(private) shared(maseg_c,keyav_c,keyag_c,keyag_c_lin) &
!$omp shared(keybg_c,mbseg_c,keybv_c,mbseg_f,maseg_f,apsi_c,bpsi1,bpsi2) &
!$omp shared(bpsi3,bpsi4,bpsi5,bpsi_f1,bpsi_f2,bpsi_f3,bpsi_f4,bpsi_f5) &
!$omp shared(keybv_f,keybg_f,keyag_f,keyag_f_lin,keyav_f,apsi_f,scpr)

!!!!$omp shared (ncount0,ncount2,ncount_rate,ncount_max,tel)

  iaseg0=1 

  scpr0=0.0_dp
  scpr1=0.0_dp
  scpr2=0.0_dp
  scpr3=0.0_dp
  scpr4=0.0_dp
  scpr5=0.0_dp
  scpr6=0.0_dp
  scpr7=0.0_dp
!coarse part. Loop on the projectors segments
!$omp do schedule(static)
   do ibseg=1,mbseg_c
     jbj=keybv_c(ibseg)
     jb0=keybg_c(1,ibseg) !starting point of projector segment
     jb1=keybg_c(2,ibseg) !ending point of projector segment
!     print *,'huntenter',ibseg,jb0,jb1
 
     !find the starting point of the wavefunction segment
     !warning: hunt is assuming that the variable is always found
     !if it is not, iaseg0 is put to maseg + 1 so that the loop is disabled
     call hunt1(.true.,keyag_c_lin,maseg_c,keybg_c(1,ibseg),iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     !now pass through all the wavefunction segments until the end of the segment is 
     !still contained in projector segment
     nonconvex_loop_c: do while(iaseg0 <= maseg_c)
!     print *,'huntexit',iaseg0,maseg_c,keyag_c_lin(iaseg0),keyag_c(2,iaseg0)

        !length = jb1-jb0
        !iaoff = jb0-keyag_c_lin(iaseg0)!jb0-ja0

        ja0=keyag_c_lin(iaseg0)
        ja1=min(jb1,keyag_c(2,iaseg0)) 
        length = ja1-jb0
        iaoff = max(jb0-ja0,0) !no offset if we are already inside

        jaj=keyav_c(iaseg0)
        do i=0,length
           apsi_temp=real(apsi_c(jaj+iaoff+i),dp)

           scpr0(1)=scpr0(1)+apsi_temp *real(bpsi1(jbj+i),dp)
           scpr0(2)=scpr0(2)+apsi_temp *real(bpsi2(jbj+i),dp)
           scpr0(3)=scpr0(3)+apsi_temp *real(bpsi3(jbj+i),dp)
           scpr0(4)=scpr0(4)+apsi_temp *real(bpsi4(jbj+i),dp)
           scpr0(5)=scpr0(5)+apsi_temp *real(bpsi5(jbj+i),dp)
        enddo
 !       print *,'length',length,ibseg,scpr0,iaseg0,ja1,jb1

        !print *,'ibseg,mbseg_c,iaseg0,maseg_c',ibseg,mbseg_c,iaseg0,maseg_c
        !print '(a,6(i8),1pe25.17)','ja0,ja1t,ja1,jb0,jb1',&
        !     ibseg,ja0,keyag_c(2,iaseg0),ja1,jb0,jb1,scpr0
        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_c) exit nonconvex_loop_c !segment is finished
        iaseg0=iaseg0+1
        jb0=max(keybg_c(1,ibseg),keyag_c_lin(iaseg0))
        if (keyag_c_lin(iaseg0)>jb1) exit nonconvex_loop_c !segment is not covered
        jbj=jbj+max(jb0-keybg_c(1,ibseg),0)
     end do nonconvex_loop_c
     !disable loop if the end is reached
     if (iaseg0 == maseg_c .and. keybg_c(1,ibseg)> keyag_c_lin(maseg_c)) iaseg0=iaseg0+1


   enddo
!stop
!$omp end do nowait

! fine part

iaseg0=1

!$omp do schedule(static)
   do ibseg=1,mbseg_f
     jbj=keybv_f(ibseg)
     jb0=keybg_f(1,ibseg)
     jb1=keybg_f(2,ibseg)
!    print *,'huntenter',ibseg,jb0,jb1
     call hunt1(.true.,keyag_f_lin,maseg_f,keybg_f(1,ibseg),iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     nonconvex_loop_f: do while(iaseg0 <= maseg_f)
!!$     length = jb1-jb0
!!$     iaoff = jb0-keyag_f_lin(iaseg0)

        ja0=keyag_f_lin(iaseg0) !still doubts about copying in automatic array
        ja1=min(jb1,keyag_f(2,iaseg0)) 
        length = ja1-jb0
        iaoff = max(jb0-ja0,0) !no offset if we are already inside

        jaj=keyav_f(iaseg0)
        do i=0,length
           apsi_f_temp(1)=real(apsi_f(1,jaj+iaoff+i),dp)
           apsi_f_temp(2)=real(apsi_f(2,jaj+iaoff+i),dp)
           apsi_f_temp(3)=real(apsi_f(3,jaj+iaoff+i),dp)
           apsi_f_temp(4)=real(apsi_f(4,jaj+iaoff+i),dp)
           apsi_f_temp(5)=real(apsi_f(5,jaj+iaoff+i),dp)
           apsi_f_temp(6)=real(apsi_f(6,jaj+iaoff+i),dp)
           apsi_f_temp(7)=real(apsi_f(7,jaj+iaoff+i),dp)

           scpr1(1)=scpr1(1)+apsi_f_temp(1)*real(bpsi_f1(1,jbj+i),dp)
           scpr2(1)=scpr2(1)+apsi_f_temp(2)*real(bpsi_f1(2,jbj+i),dp)
           scpr3(1)=scpr3(1)+apsi_f_temp(3)*real(bpsi_f1(3,jbj+i),dp)
           scpr4(1)=scpr4(1)+apsi_f_temp(4)*real(bpsi_f1(4,jbj+i),dp)
           scpr5(1)=scpr5(1)+apsi_f_temp(5)*real(bpsi_f1(5,jbj+i),dp)
           scpr6(1)=scpr6(1)+apsi_f_temp(6)*real(bpsi_f1(6,jbj+i),dp)
           scpr7(1)=scpr7(1)+apsi_f_temp(7)*real(bpsi_f1(7,jbj+i),dp)

           scpr1(2)=scpr1(2)+apsi_f_temp(1)*real(bpsi_f2(1,jbj+i),dp)
           scpr2(2)=scpr2(2)+apsi_f_temp(2)*real(bpsi_f2(2,jbj+i),dp)
           scpr3(2)=scpr3(2)+apsi_f_temp(3)*real(bpsi_f2(3,jbj+i),dp)
           scpr4(2)=scpr4(2)+apsi_f_temp(4)*real(bpsi_f2(4,jbj+i),dp)
           scpr5(2)=scpr5(2)+apsi_f_temp(5)*real(bpsi_f2(5,jbj+i),dp)
           scpr6(2)=scpr6(2)+apsi_f_temp(6)*real(bpsi_f2(6,jbj+i),dp)
           scpr7(2)=scpr7(2)+apsi_f_temp(7)*real(bpsi_f2(7,jbj+i),dp)

           scpr1(3)=scpr1(3)+apsi_f_temp(1)*real(bpsi_f3(1,jbj+i),dp)
           scpr2(3)=scpr2(3)+apsi_f_temp(2)*real(bpsi_f3(2,jbj+i),dp)
           scpr3(3)=scpr3(3)+apsi_f_temp(3)*real(bpsi_f3(3,jbj+i),dp)
           scpr4(3)=scpr4(3)+apsi_f_temp(4)*real(bpsi_f3(4,jbj+i),dp)
           scpr5(3)=scpr5(3)+apsi_f_temp(5)*real(bpsi_f3(5,jbj+i),dp)
           scpr6(3)=scpr6(3)+apsi_f_temp(6)*real(bpsi_f3(6,jbj+i),dp)
           scpr7(3)=scpr7(3)+apsi_f_temp(7)*real(bpsi_f3(7,jbj+i),dp)

           scpr1(4)=scpr1(4)+apsi_f_temp(1)*real(bpsi_f4(1,jbj+i),dp)
           scpr2(4)=scpr2(4)+apsi_f_temp(2)*real(bpsi_f4(2,jbj+i),dp)
           scpr3(4)=scpr3(4)+apsi_f_temp(3)*real(bpsi_f4(3,jbj+i),dp)
           scpr4(4)=scpr4(4)+apsi_f_temp(4)*real(bpsi_f4(4,jbj+i),dp)
           scpr5(4)=scpr5(4)+apsi_f_temp(5)*real(bpsi_f4(5,jbj+i),dp)
           scpr6(4)=scpr6(4)+apsi_f_temp(6)*real(bpsi_f4(6,jbj+i),dp)
           scpr7(4)=scpr7(4)+apsi_f_temp(7)*real(bpsi_f4(7,jbj+i),dp)

           scpr1(5)=scpr1(5)+apsi_f_temp(1)*real(bpsi_f5(1,jbj+i),dp)
           scpr2(5)=scpr2(5)+apsi_f_temp(2)*real(bpsi_f5(2,jbj+i),dp)
           scpr3(5)=scpr3(5)+apsi_f_temp(3)*real(bpsi_f5(3,jbj+i),dp)
           scpr4(5)=scpr4(5)+apsi_f_temp(4)*real(bpsi_f5(4,jbj+i),dp)
           scpr5(5)=scpr5(5)+apsi_f_temp(5)*real(bpsi_f5(5,jbj+i),dp)
           scpr6(5)=scpr6(5)+apsi_f_temp(6)*real(bpsi_f5(6,jbj+i),dp)
           scpr7(5)=scpr7(5)+apsi_f_temp(7)*real(bpsi_f5(7,jbj+i),dp)

        enddo
 !       print *,'length',length,ibseg,scpr1,iaseg0,ja1,jb1
        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_f) exit nonconvex_loop_f !segment is finished  
        iaseg0=iaseg0+1
        jb0=max(keybg_f(1,ibseg),keyag_f_lin(iaseg0))
        if (keyag_f_lin(iaseg0)>jb1) exit nonconvex_loop_f !segment is not covered 
        jbj=jbj+max(jb0-keybg_f(1,ibseg),0)
     end do nonconvex_loop_f
     !disable loop if the end is reached
     if (iaseg0 == maseg_f .and. keybg_f(1,ibseg)> keyag_f_lin(maseg_f)) iaseg0=iaseg0+1

   enddo
!$omp end do !!!implicit barrier 

   scpr0=scpr0+scpr1+scpr2+scpr3+scpr4+scpr5+scpr6+scpr7

!$omp critical 
   scpr=scpr+scpr0
!$omp end critical
  
!$omp end parallel

!!!    call system_clock(ncount2,ncount_rate,ncount_max)
!!!    tel=dble(ncount2-ncount0)/dble(ncount_rate)
!!!    write(97,*) 'wpdot:',tel
!!!    close(97)

END SUBROUTINE wpdot_5

subroutine wpdot_8(  &
     mavctr_c,mavctr_f,maseg_c,maseg_f,keyav_c,keyav_f,keyag_c,keyag_f,apsi_c,apsi_f,  &
     mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keybv_c,keybv_f,keybg_c,keybg_f, &
     bpsi1,bpsi2,bpsi3,bpsi4,bpsi5,bpsi6,bpsi7,bpsi8, &
     bpsi_f1,bpsi_f2,bpsi_f3,bpsi_f4,bpsi_f5,bpsi_f6,bpsi_f7,bpsi_f8, &
     scpr,proj_count)

  use module_base
  implicit none
  integer, intent(in) :: mavctr_c,mavctr_f,maseg_c,maseg_f,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f
  integer, intent(in) :: proj_count
  integer, dimension(maseg_c), intent(in) :: keyav_c
  integer, dimension(mbseg_c), intent(in) :: keybv_c
  integer, dimension(maseg_f), intent(in) :: keyav_f
  integer, dimension(mbseg_f), intent(in) :: keybv_f
  integer, dimension(2,maseg_c), intent(in) :: keyag_c
  integer, dimension(2,mbseg_c), intent(in) :: keybg_c
  integer, dimension(2,maseg_f), intent(in) :: keyag_f
  integer, dimension(2,mbseg_f), intent(in) :: keybg_f
  real(wp), dimension(mavctr_c), intent(in) :: apsi_c
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi1,bpsi2,bpsi3,bpsi4,bpsi5,bpsi6,bpsi7,bpsi8
  real(wp), dimension(7,mavctr_f), intent(in) :: apsi_f
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f1,bpsi_f2,bpsi_f3,bpsi_f4,bpsi_f5,bpsi_f6,bpsi_f7
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f8
  real(dp), dimension(proj_count),intent(out) :: scpr
  !local variables
  integer :: ibseg,jaj,jb1,jb0,jbj,iaoff,length,i,ja0,ja1
  real(dp),dimension(proj_count) :: scpr1,scpr2,scpr3,scpr4,scpr5,scpr6,scpr7,scpr0
  integer :: iaseg0
  integer, dimension(maseg_c) :: keyag_c_lin !>linear version of second indices of keyag_c
  integer, dimension(maseg_f) :: keyag_f_lin !>linear version of second indices of keyag_f
  real(dp) :: apsi_temp
  real(dp), dimension(7) :: apsi_f_temp
!!!    integer :: ncount0,ncount2,ncount_rate,ncount_max
!!!    real(gp) :: tel

!!!  !dee
!!!    open(unit=97,file='time_wpdot',status='unknown',position='append')
!!!    call system_clock(ncount0,ncount_rate,ncount_max)

  keyag_c_lin = keyag_c(1,:)!speed up access in hunt subroutine by consecutive arrangement in memory
  keyag_f_lin = keyag_f(1,:)!speed up access in hunt subroutine by consecutive arrangement in memory

  scpr=0.0_dp

!$omp parallel default(private) shared(maseg_c,keyav_c,keyag_c,keyag_c_lin) &
!$omp shared(keybg_c,mbseg_c,keybv_c,mbseg_f,maseg_f,apsi_c,bpsi1,bpsi2) &
!$omp shared(bpsi3,bpsi4,bpsi5,bpsi6,bpsi7,bpsi8,bpsi_f1,bpsi_f2,bpsi_f3) &
!$omp shared(bpsi_f4,bpsi_f5,bpsi_f6,bpsi_f7,bpsi_f8) &
!$omp shared(keybv_f,keybg_f,keyag_f,keyag_f_lin,keyav_f) &
!$omp shared(apsi_f,scpr)

!!!!$omp shared (ncount0,ncount2,ncount_rate,ncount_max,tel)

  iaseg0=1 

  scpr0=0.0_dp
  scpr1=0.0_dp
  scpr2=0.0_dp
  scpr3=0.0_dp
  scpr4=0.0_dp
  scpr5=0.0_dp
  scpr6=0.0_dp
  scpr7=0.0_dp
!coarse part. Loop on the projectors segments
!$omp do schedule(static)
   do ibseg=1,mbseg_c
     jbj=keybv_c(ibseg)
     jb0=keybg_c(1,ibseg) !starting point of projector segment
     jb1=keybg_c(2,ibseg) !ending point of projector segment
!     print *,'huntenter',ibseg,jb0,jb1
 
     !find the starting point of the wavefunction segment
     !warning: hunt is assuming that the variable is always found
     !if it is not, iaseg0 is put to maseg + 1 so that the loop is disabled
     call hunt1(.true.,keyag_c_lin,maseg_c,keybg_c(1,ibseg),iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     !now pass through all the wavefunction segments until the end of the segment is 
     !still contained in projector segment
     nonconvex_loop_c: do while(iaseg0 <= maseg_c)
!     print *,'huntexit',iaseg0,maseg_c,keyag_c_lin(iaseg0),keyag_c(2,iaseg0)

        !length = jb1-jb0
        !iaoff = jb0-keyag_c_lin(iaseg0)!jb0-ja0

        ja0=keyag_c_lin(iaseg0)
        ja1=min(jb1,keyag_c(2,iaseg0)) 
        length = ja1-jb0
        iaoff = max(jb0-ja0,0) !no offset if we are already inside

        jaj=keyav_c(iaseg0)
        do i=0,length
           apsi_temp=real(apsi_c(jaj+iaoff+i),dp)

           scpr0(1)=scpr0(1)+apsi_temp *real(bpsi1(jbj+i),dp)
           scpr0(2)=scpr0(2)+apsi_temp *real(bpsi2(jbj+i),dp)
           scpr0(3)=scpr0(3)+apsi_temp *real(bpsi3(jbj+i),dp)
           scpr0(4)=scpr0(4)+apsi_temp *real(bpsi4(jbj+i),dp)
           scpr0(5)=scpr0(5)+apsi_temp *real(bpsi5(jbj+i),dp)
           scpr0(6)=scpr0(6)+apsi_temp *real(bpsi6(jbj+i),dp)
           scpr0(7)=scpr0(7)+apsi_temp *real(bpsi7(jbj+i),dp)
           scpr0(8)=scpr0(8)+apsi_temp *real(bpsi8(jbj+i),dp)
        enddo
 !       print *,'length',length,ibseg,scpr0,iaseg0,ja1,jb1

        !print *,'ibseg,mbseg_c,iaseg0,maseg_c',ibseg,mbseg_c,iaseg0,maseg_c
        !print '(a,6(i8),1pe25.17)','ja0,ja1t,ja1,jb0,jb1',&
        !     ibseg,ja0,keyag_c(2,iaseg0),ja1,jb0,jb1,scpr0
        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_c) exit nonconvex_loop_c !segment is finished
        iaseg0=iaseg0+1
        jb0=max(keybg_c(1,ibseg),keyag_c_lin(iaseg0))
        if (keyag_c_lin(iaseg0)>jb1) exit nonconvex_loop_c !segment is not covered
        jbj=jbj+max(jb0-keybg_c(1,ibseg),0)
     end do nonconvex_loop_c
     !disable loop if the end is reached
     if (iaseg0 == maseg_c .and. keybg_c(1,ibseg)> keyag_c_lin(maseg_c)) iaseg0=iaseg0+1


   enddo
!stop
!$omp end do nowait

! fine part

iaseg0=1

!$omp do schedule(static)
   do ibseg=1,mbseg_f
     jbj=keybv_f(ibseg)
     jb0=keybg_f(1,ibseg)
     jb1=keybg_f(2,ibseg)
!    print *,'huntenter',ibseg,jb0,jb1
     call hunt1(.true.,keyag_f_lin,maseg_f,keybg_f(1,ibseg),iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     nonconvex_loop_f: do while(iaseg0 <= maseg_f)
!!$     length = jb1-jb0
!!$     iaoff = jb0-keyag_f_lin(iaseg0)

        ja0=keyag_f_lin(iaseg0) !still doubts about copying in automatic array
        ja1=min(jb1,keyag_f(2,iaseg0)) 
        length = ja1-jb0
        iaoff = max(jb0-ja0,0) !no offset if we are already inside

        jaj=keyav_f(iaseg0)
        do i=0,length
           apsi_f_temp(1)=real(apsi_f(1,jaj+iaoff+i),dp)
           apsi_f_temp(2)=real(apsi_f(2,jaj+iaoff+i),dp)
           apsi_f_temp(3)=real(apsi_f(3,jaj+iaoff+i),dp)
           apsi_f_temp(4)=real(apsi_f(4,jaj+iaoff+i),dp)
           apsi_f_temp(5)=real(apsi_f(5,jaj+iaoff+i),dp)
           apsi_f_temp(6)=real(apsi_f(6,jaj+iaoff+i),dp)
           apsi_f_temp(7)=real(apsi_f(7,jaj+iaoff+i),dp)

           scpr1(1)=scpr1(1)+apsi_f_temp(1)*real(bpsi_f1(1,jbj+i),dp)
           scpr2(1)=scpr2(1)+apsi_f_temp(2)*real(bpsi_f1(2,jbj+i),dp)
           scpr3(1)=scpr3(1)+apsi_f_temp(3)*real(bpsi_f1(3,jbj+i),dp)
           scpr4(1)=scpr4(1)+apsi_f_temp(4)*real(bpsi_f1(4,jbj+i),dp)
           scpr5(1)=scpr5(1)+apsi_f_temp(5)*real(bpsi_f1(5,jbj+i),dp)
           scpr6(1)=scpr6(1)+apsi_f_temp(6)*real(bpsi_f1(6,jbj+i),dp)
           scpr7(1)=scpr7(1)+apsi_f_temp(7)*real(bpsi_f1(7,jbj+i),dp)

           scpr1(2)=scpr1(2)+apsi_f_temp(1)*real(bpsi_f2(1,jbj+i),dp)
           scpr2(2)=scpr2(2)+apsi_f_temp(2)*real(bpsi_f2(2,jbj+i),dp)
           scpr3(2)=scpr3(2)+apsi_f_temp(3)*real(bpsi_f2(3,jbj+i),dp)
           scpr4(2)=scpr4(2)+apsi_f_temp(4)*real(bpsi_f2(4,jbj+i),dp)
           scpr5(2)=scpr5(2)+apsi_f_temp(5)*real(bpsi_f2(5,jbj+i),dp)
           scpr6(2)=scpr6(2)+apsi_f_temp(6)*real(bpsi_f2(6,jbj+i),dp)
           scpr7(2)=scpr7(2)+apsi_f_temp(7)*real(bpsi_f2(7,jbj+i),dp)

           scpr1(3)=scpr1(3)+apsi_f_temp(1)*real(bpsi_f3(1,jbj+i),dp)
           scpr2(3)=scpr2(3)+apsi_f_temp(2)*real(bpsi_f3(2,jbj+i),dp)
           scpr3(3)=scpr3(3)+apsi_f_temp(3)*real(bpsi_f3(3,jbj+i),dp)
           scpr4(3)=scpr4(3)+apsi_f_temp(4)*real(bpsi_f3(4,jbj+i),dp)
           scpr5(3)=scpr5(3)+apsi_f_temp(5)*real(bpsi_f3(5,jbj+i),dp)
           scpr6(3)=scpr6(3)+apsi_f_temp(6)*real(bpsi_f3(6,jbj+i),dp)
           scpr7(3)=scpr7(3)+apsi_f_temp(7)*real(bpsi_f3(7,jbj+i),dp)

           scpr1(4)=scpr1(4)+apsi_f_temp(1)*real(bpsi_f4(1,jbj+i),dp)
           scpr2(4)=scpr2(4)+apsi_f_temp(2)*real(bpsi_f4(2,jbj+i),dp)
           scpr3(4)=scpr3(4)+apsi_f_temp(3)*real(bpsi_f4(3,jbj+i),dp)
           scpr4(4)=scpr4(4)+apsi_f_temp(4)*real(bpsi_f4(4,jbj+i),dp)
           scpr5(4)=scpr5(4)+apsi_f_temp(5)*real(bpsi_f4(5,jbj+i),dp)
           scpr6(4)=scpr6(4)+apsi_f_temp(6)*real(bpsi_f4(6,jbj+i),dp)
           scpr7(4)=scpr7(4)+apsi_f_temp(7)*real(bpsi_f4(7,jbj+i),dp)

           scpr1(5)=scpr1(5)+apsi_f_temp(1)*real(bpsi_f5(1,jbj+i),dp)
           scpr2(5)=scpr2(5)+apsi_f_temp(2)*real(bpsi_f5(2,jbj+i),dp)
           scpr3(5)=scpr3(5)+apsi_f_temp(3)*real(bpsi_f5(3,jbj+i),dp)
           scpr4(5)=scpr4(5)+apsi_f_temp(4)*real(bpsi_f5(4,jbj+i),dp)
           scpr5(5)=scpr5(5)+apsi_f_temp(5)*real(bpsi_f5(5,jbj+i),dp)
           scpr6(5)=scpr6(5)+apsi_f_temp(6)*real(bpsi_f5(6,jbj+i),dp)
           scpr7(5)=scpr7(5)+apsi_f_temp(7)*real(bpsi_f5(7,jbj+i),dp)

           scpr1(6)=scpr1(6)+apsi_f_temp(1)*real(bpsi_f6(1,jbj+i),dp)
           scpr2(6)=scpr2(6)+apsi_f_temp(2)*real(bpsi_f6(2,jbj+i),dp)
           scpr3(6)=scpr3(6)+apsi_f_temp(3)*real(bpsi_f6(3,jbj+i),dp)
           scpr4(6)=scpr4(6)+apsi_f_temp(4)*real(bpsi_f6(4,jbj+i),dp)
           scpr5(6)=scpr5(6)+apsi_f_temp(5)*real(bpsi_f6(5,jbj+i),dp)
           scpr6(6)=scpr6(6)+apsi_f_temp(6)*real(bpsi_f6(6,jbj+i),dp)
           scpr7(6)=scpr7(6)+apsi_f_temp(7)*real(bpsi_f6(7,jbj+i),dp)

           scpr1(7)=scpr1(7)+apsi_f_temp(1)*real(bpsi_f7(1,jbj+i),dp)
           scpr2(7)=scpr2(7)+apsi_f_temp(2)*real(bpsi_f7(2,jbj+i),dp)
           scpr3(7)=scpr3(7)+apsi_f_temp(3)*real(bpsi_f7(3,jbj+i),dp)
           scpr4(7)=scpr4(7)+apsi_f_temp(4)*real(bpsi_f7(4,jbj+i),dp)
           scpr5(7)=scpr5(7)+apsi_f_temp(5)*real(bpsi_f7(5,jbj+i),dp)
           scpr6(7)=scpr6(7)+apsi_f_temp(6)*real(bpsi_f7(6,jbj+i),dp)
           scpr7(7)=scpr7(7)+apsi_f_temp(7)*real(bpsi_f7(7,jbj+i),dp)

           scpr1(8)=scpr1(8)+apsi_f_temp(1)*real(bpsi_f8(1,jbj+i),dp)
           scpr2(8)=scpr2(8)+apsi_f_temp(2)*real(bpsi_f8(2,jbj+i),dp)
           scpr3(8)=scpr3(8)+apsi_f_temp(3)*real(bpsi_f8(3,jbj+i),dp)
           scpr4(8)=scpr4(8)+apsi_f_temp(4)*real(bpsi_f8(4,jbj+i),dp)
           scpr5(8)=scpr5(8)+apsi_f_temp(5)*real(bpsi_f8(5,jbj+i),dp)
           scpr6(8)=scpr6(8)+apsi_f_temp(6)*real(bpsi_f8(6,jbj+i),dp)
           scpr7(8)=scpr7(8)+apsi_f_temp(7)*real(bpsi_f8(7,jbj+i),dp)

        enddo
 !       print *,'length',length,ibseg,scpr1,iaseg0,ja1,jb1
        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_f) exit nonconvex_loop_f !segment is finished  
        iaseg0=iaseg0+1
        jb0=max(keybg_f(1,ibseg),keyag_f_lin(iaseg0))
        if (keyag_f_lin(iaseg0)>jb1) exit nonconvex_loop_f !segment is not covered 
        jbj=jbj+max(jb0-keybg_f(1,ibseg),0)
     end do nonconvex_loop_f
     !disable loop if the end is reached
     if (iaseg0 == maseg_f .and. keybg_f(1,ibseg)> keyag_f_lin(maseg_f)) iaseg0=iaseg0+1

   enddo
!$omp end do !!!implicit barrier 

   scpr0=scpr0+scpr1+scpr2+scpr3+scpr4+scpr5+scpr6+scpr7

!$omp critical 
   scpr=scpr+scpr0
!$omp end critical
  
!$omp end parallel

!!!    call system_clock(ncount2,ncount_rate,ncount_max)
!!!    tel=dble(ncount2-ncount0)/dble(ncount_rate)
!!!    write(97,*) 'wpdot:',tel
!!!    close(97)

END SUBROUTINE wpdot_8

subroutine wpdot_13(  &
     mavctr_c,mavctr_f,maseg_c,maseg_f,keyav_c,keyav_f,keyag_c,keyag_f,apsi_c,apsi_f,  &
     mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keybv_c,keybv_f,keybg_c,keybg_f, &
     bpsi1,bpsi2,bpsi3,bpsi4,bpsi5,bpsi6,bpsi7,bpsi8,bpsi9,bpsi10,bpsi11,bpsi12,bpsi13, &
     bpsi_f1,bpsi_f2,bpsi_f3,bpsi_f4,bpsi_f5,bpsi_f6,bpsi_f7,bpsi_f8,bpsi_f9,bpsi_f10,bpsi_f11, &
     bpsi_f12,bpsi_f13, &
     scpr,proj_count)

  use module_base
  implicit none
  integer, intent(in) :: mavctr_c,mavctr_f,maseg_c,maseg_f,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f
  integer, intent(in) :: proj_count
  integer, dimension(maseg_c), intent(in) :: keyav_c
  integer, dimension(mbseg_c), intent(in) :: keybv_c
  integer, dimension(maseg_f), intent(in) :: keyav_f
  integer, dimension(mbseg_f), intent(in) :: keybv_f
  integer, dimension(2,maseg_c), intent(in) :: keyag_c
  integer, dimension(2,mbseg_c), intent(in) :: keybg_c
  integer, dimension(2,maseg_f), intent(in) :: keyag_f
  integer, dimension(2,mbseg_f), intent(in) :: keybg_f
  real(wp), dimension(mavctr_c), intent(in) :: apsi_c
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi1,bpsi2,bpsi3,bpsi4,bpsi5,bpsi6,bpsi7
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi8,bpsi9,bpsi10,bpsi11,bpsi12,bpsi13
  real(wp), dimension(7,mavctr_f), intent(in) :: apsi_f
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f1,bpsi_f2,bpsi_f3,bpsi_f4,bpsi_f5,bpsi_f6,bpsi_f7
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f8,bpsi_f9,bpsi_f10,bpsi_f11,bpsi_f12,bpsi_f13
  real(dp), dimension(proj_count),intent(out) :: scpr
  !local variables
  integer :: ibseg,jaj,jb1,jb0,jbj,iaoff,length,i,ja0,ja1
  real(dp),dimension(proj_count) :: scpr1,scpr2,scpr3,scpr4,scpr5,scpr6,scpr7,scpr0
  integer :: iaseg0
  integer, dimension(maseg_c) :: keyag_c_lin !>linear version of second indices of keyag_c
  integer, dimension(maseg_f) :: keyag_f_lin !>linear version of second indices of keyag_f
  real(dp) :: apsi_temp
  real(dp), dimension(7) :: apsi_f_temp
!!!    integer :: ncount0,ncount2,ncount_rate,ncount_max
!!!    real(gp) :: tel

!!!  !dee
!!!    open(unit=97,file='time_wpdot',status='unknown',position='append')
!!!    call system_clock(ncount0,ncount_rate,ncount_max)

  keyag_c_lin = keyag_c(1,:)!speed up access in hunt subroutine by consecutive arrangement in memory
  keyag_f_lin = keyag_f(1,:)!speed up access in hunt subroutine by consecutive arrangement in memory

  scpr=0.0_dp

!$omp parallel default(private) shared(maseg_c,keyav_c,keyag_c,keyag_c_lin) &
!$omp shared(keybg_c,mbseg_c,keybv_c,mbseg_f,maseg_f,apsi_c,bpsi1,bpsi2) &
!$omp shared(bpsi3,bpsi4,bpsi5,bpsi6,bpsi7,bpsi8,bpsi9,bpsi10,bpsi11) &
!$omp shared(bpsi12,bpsi13,bpsi_f1,bpsi_f2,bpsi_f3,bpsi_f4,bpsi_f5) &
!$omp shared(bpsi_f6,bpsi_f7,bpsi_f8,bpsi_f9,bpsi_f10,bpsi_f11,bpsi_f12) &
!$omp shared(bpsi_f13,keybv_f,keybg_f,keyag_f,keyag_f_lin,keyav_f) &
!$omp shared(apsi_f,scpr)

!!!!$omp shared (ncount0,ncount2,ncount_rate,ncount_max,tel)

  iaseg0=1 

  scpr0=0.0_dp
  scpr1=0.0_dp
  scpr2=0.0_dp
  scpr3=0.0_dp
  scpr4=0.0_dp
  scpr5=0.0_dp
  scpr6=0.0_dp
  scpr7=0.0_dp
!coarse part. Loop on the projectors segments
!$omp do schedule(static)
   do ibseg=1,mbseg_c
     jbj=keybv_c(ibseg)
     jb0=keybg_c(1,ibseg) !starting point of projector segment
     jb1=keybg_c(2,ibseg) !ending point of projector segment
!     print *,'huntenter',ibseg,jb0,jb1
 
     !find the starting point of the wavefunction segment
     !warning: hunt is assuming that the variable is always found
     !if it is not, iaseg0 is put to maseg + 1 so that the loop is disabled
     call hunt1(.true.,keyag_c_lin,maseg_c,keybg_c(1,ibseg),iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     !now pass through all the wavefunction segments until the end of the segment is 
     !still contained in projector segment
     nonconvex_loop_c: do while(iaseg0 <= maseg_c)
!     print *,'huntexit',iaseg0,maseg_c,keyag_c_lin(iaseg0),keyag_c(2,iaseg0)

        !length = jb1-jb0
        !iaoff = jb0-keyag_c_lin(iaseg0)!jb0-ja0

        ja0=keyag_c_lin(iaseg0)
        ja1=min(jb1,keyag_c(2,iaseg0)) 
        length = ja1-jb0
        iaoff = max(jb0-ja0,0) !no offset if we are already inside

        jaj=keyav_c(iaseg0)
        do i=0,length
           apsi_temp=real(apsi_c(jaj+iaoff+i),dp)

           scpr0(1)=scpr0(1)+apsi_temp *real(bpsi1(jbj+i),dp)
           scpr0(2)=scpr0(2)+apsi_temp *real(bpsi2(jbj+i),dp)
           scpr0(3)=scpr0(3)+apsi_temp *real(bpsi3(jbj+i),dp)
           scpr0(4)=scpr0(4)+apsi_temp *real(bpsi4(jbj+i),dp)
           scpr0(5)=scpr0(5)+apsi_temp *real(bpsi5(jbj+i),dp)
           scpr0(6)=scpr0(6)+apsi_temp *real(bpsi6(jbj+i),dp)
           scpr0(7)=scpr0(7)+apsi_temp *real(bpsi7(jbj+i),dp)
           scpr0(8)=scpr0(8)+apsi_temp *real(bpsi8(jbj+i),dp)
           scpr0(9)=scpr0(9)+apsi_temp *real(bpsi9(jbj+i),dp)
           scpr0(10)=scpr0(10)+apsi_temp *real(bpsi10(jbj+i),dp)
           scpr0(11)=scpr0(11)+apsi_temp *real(bpsi11(jbj+i),dp)
           scpr0(12)=scpr0(12)+apsi_temp *real(bpsi12(jbj+i),dp)
           scpr0(13)=scpr0(13)+apsi_temp *real(bpsi13(jbj+i),dp)
        enddo
 !       print *,'length',length,ibseg,scpr0,iaseg0,ja1,jb1

        !print *,'ibseg,mbseg_c,iaseg0,maseg_c',ibseg,mbseg_c,iaseg0,maseg_c
        !print '(a,6(i8),1pe25.17)','ja0,ja1t,ja1,jb0,jb1',&
        !     ibseg,ja0,keyag_c(2,iaseg0),ja1,jb0,jb1,scpr0
        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_c) exit nonconvex_loop_c !segment is finished
        iaseg0=iaseg0+1
        jb0=max(keybg_c(1,ibseg),keyag_c_lin(iaseg0))
        if (keyag_c_lin(iaseg0)>jb1) exit nonconvex_loop_c !segment is not covered
        jbj=jbj+max(jb0-keybg_c(1,ibseg),0)
     end do nonconvex_loop_c
     !disable loop if the end is reached
     if (iaseg0 == maseg_c .and. keybg_c(1,ibseg)> keyag_c_lin(maseg_c)) iaseg0=iaseg0+1


   enddo
!stop
!$omp end do nowait

! fine part

iaseg0=1

!$omp do schedule(static)
   do ibseg=1,mbseg_f
     jbj=keybv_f(ibseg)
     jb0=keybg_f(1,ibseg)
     jb1=keybg_f(2,ibseg)
!    print *,'huntenter',ibseg,jb0,jb1
     call hunt1(.true.,keyag_f_lin,maseg_f,keybg_f(1,ibseg),iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     nonconvex_loop_f: do while(iaseg0 <= maseg_f)
!!$     length = jb1-jb0
!!$     iaoff = jb0-keyag_f_lin(iaseg0)

        ja0=keyag_f_lin(iaseg0) !still doubts about copying in automatic array
        ja1=min(jb1,keyag_f(2,iaseg0)) 
        length = ja1-jb0
        iaoff = max(jb0-ja0,0) !no offset if we are already inside

        jaj=keyav_f(iaseg0)
        do i=0,length
           apsi_f_temp(1)=real(apsi_f(1,jaj+iaoff+i),dp)
           apsi_f_temp(2)=real(apsi_f(2,jaj+iaoff+i),dp)
           apsi_f_temp(3)=real(apsi_f(3,jaj+iaoff+i),dp)
           apsi_f_temp(4)=real(apsi_f(4,jaj+iaoff+i),dp)
           apsi_f_temp(5)=real(apsi_f(5,jaj+iaoff+i),dp)
           apsi_f_temp(6)=real(apsi_f(6,jaj+iaoff+i),dp)
           apsi_f_temp(7)=real(apsi_f(7,jaj+iaoff+i),dp)

           scpr1(1)=scpr1(1)+apsi_f_temp(1)*real(bpsi_f1(1,jbj+i),dp)
           scpr2(1)=scpr2(1)+apsi_f_temp(2)*real(bpsi_f1(2,jbj+i),dp)
           scpr3(1)=scpr3(1)+apsi_f_temp(3)*real(bpsi_f1(3,jbj+i),dp)
           scpr4(1)=scpr4(1)+apsi_f_temp(4)*real(bpsi_f1(4,jbj+i),dp)
           scpr5(1)=scpr5(1)+apsi_f_temp(5)*real(bpsi_f1(5,jbj+i),dp)
           scpr6(1)=scpr6(1)+apsi_f_temp(6)*real(bpsi_f1(6,jbj+i),dp)
           scpr7(1)=scpr7(1)+apsi_f_temp(7)*real(bpsi_f1(7,jbj+i),dp)

           scpr1(2)=scpr1(2)+apsi_f_temp(1)*real(bpsi_f2(1,jbj+i),dp)
           scpr2(2)=scpr2(2)+apsi_f_temp(2)*real(bpsi_f2(2,jbj+i),dp)
           scpr3(2)=scpr3(2)+apsi_f_temp(3)*real(bpsi_f2(3,jbj+i),dp)
           scpr4(2)=scpr4(2)+apsi_f_temp(4)*real(bpsi_f2(4,jbj+i),dp)
           scpr5(2)=scpr5(2)+apsi_f_temp(5)*real(bpsi_f2(5,jbj+i),dp)
           scpr6(2)=scpr6(2)+apsi_f_temp(6)*real(bpsi_f2(6,jbj+i),dp)
           scpr7(2)=scpr7(2)+apsi_f_temp(7)*real(bpsi_f2(7,jbj+i),dp)

           scpr1(3)=scpr1(3)+apsi_f_temp(1)*real(bpsi_f3(1,jbj+i),dp)
           scpr2(3)=scpr2(3)+apsi_f_temp(2)*real(bpsi_f3(2,jbj+i),dp)
           scpr3(3)=scpr3(3)+apsi_f_temp(3)*real(bpsi_f3(3,jbj+i),dp)
           scpr4(3)=scpr4(3)+apsi_f_temp(4)*real(bpsi_f3(4,jbj+i),dp)
           scpr5(3)=scpr5(3)+apsi_f_temp(5)*real(bpsi_f3(5,jbj+i),dp)
           scpr6(3)=scpr6(3)+apsi_f_temp(6)*real(bpsi_f3(6,jbj+i),dp)
           scpr7(3)=scpr7(3)+apsi_f_temp(7)*real(bpsi_f3(7,jbj+i),dp)

           scpr1(4)=scpr1(4)+apsi_f_temp(1)*real(bpsi_f4(1,jbj+i),dp)
           scpr2(4)=scpr2(4)+apsi_f_temp(2)*real(bpsi_f4(2,jbj+i),dp)
           scpr3(4)=scpr3(4)+apsi_f_temp(3)*real(bpsi_f4(3,jbj+i),dp)
           scpr4(4)=scpr4(4)+apsi_f_temp(4)*real(bpsi_f4(4,jbj+i),dp)
           scpr5(4)=scpr5(4)+apsi_f_temp(5)*real(bpsi_f4(5,jbj+i),dp)
           scpr6(4)=scpr6(4)+apsi_f_temp(6)*real(bpsi_f4(6,jbj+i),dp)
           scpr7(4)=scpr7(4)+apsi_f_temp(7)*real(bpsi_f4(7,jbj+i),dp)

           scpr1(5)=scpr1(5)+apsi_f_temp(1)*real(bpsi_f5(1,jbj+i),dp)
           scpr2(5)=scpr2(5)+apsi_f_temp(2)*real(bpsi_f5(2,jbj+i),dp)
           scpr3(5)=scpr3(5)+apsi_f_temp(3)*real(bpsi_f5(3,jbj+i),dp)
           scpr4(5)=scpr4(5)+apsi_f_temp(4)*real(bpsi_f5(4,jbj+i),dp)
           scpr5(5)=scpr5(5)+apsi_f_temp(5)*real(bpsi_f5(5,jbj+i),dp)
           scpr6(5)=scpr6(5)+apsi_f_temp(6)*real(bpsi_f5(6,jbj+i),dp)
           scpr7(5)=scpr7(5)+apsi_f_temp(7)*real(bpsi_f5(7,jbj+i),dp)

           scpr1(6)=scpr1(6)+apsi_f_temp(1)*real(bpsi_f6(1,jbj+i),dp)
           scpr2(6)=scpr2(6)+apsi_f_temp(2)*real(bpsi_f6(2,jbj+i),dp)
           scpr3(6)=scpr3(6)+apsi_f_temp(3)*real(bpsi_f6(3,jbj+i),dp)
           scpr4(6)=scpr4(6)+apsi_f_temp(4)*real(bpsi_f6(4,jbj+i),dp)
           scpr5(6)=scpr5(6)+apsi_f_temp(5)*real(bpsi_f6(5,jbj+i),dp)
           scpr6(6)=scpr6(6)+apsi_f_temp(6)*real(bpsi_f6(6,jbj+i),dp)
           scpr7(6)=scpr7(6)+apsi_f_temp(7)*real(bpsi_f6(7,jbj+i),dp)

           scpr1(7)=scpr1(7)+apsi_f_temp(1)*real(bpsi_f7(1,jbj+i),dp)
           scpr2(7)=scpr2(7)+apsi_f_temp(2)*real(bpsi_f7(2,jbj+i),dp)
           scpr3(7)=scpr3(7)+apsi_f_temp(3)*real(bpsi_f7(3,jbj+i),dp)
           scpr4(7)=scpr4(7)+apsi_f_temp(4)*real(bpsi_f7(4,jbj+i),dp)
           scpr5(7)=scpr5(7)+apsi_f_temp(5)*real(bpsi_f7(5,jbj+i),dp)
           scpr6(7)=scpr6(7)+apsi_f_temp(6)*real(bpsi_f7(6,jbj+i),dp)
           scpr7(7)=scpr7(7)+apsi_f_temp(7)*real(bpsi_f7(7,jbj+i),dp)

           scpr1(8)=scpr1(8)+apsi_f_temp(1)*real(bpsi_f8(1,jbj+i),dp)
           scpr2(8)=scpr2(8)+apsi_f_temp(2)*real(bpsi_f8(2,jbj+i),dp)
           scpr3(8)=scpr3(8)+apsi_f_temp(3)*real(bpsi_f8(3,jbj+i),dp)
           scpr4(8)=scpr4(8)+apsi_f_temp(4)*real(bpsi_f8(4,jbj+i),dp)
           scpr5(8)=scpr5(8)+apsi_f_temp(5)*real(bpsi_f8(5,jbj+i),dp)
           scpr6(8)=scpr6(8)+apsi_f_temp(6)*real(bpsi_f8(6,jbj+i),dp)
           scpr7(8)=scpr7(8)+apsi_f_temp(7)*real(bpsi_f8(7,jbj+i),dp)

           scpr1(9)=scpr1(9)+apsi_f_temp(1)*real(bpsi_f9(1,jbj+i),dp)
           scpr2(9)=scpr2(9)+apsi_f_temp(2)*real(bpsi_f9(2,jbj+i),dp)
           scpr3(9)=scpr3(9)+apsi_f_temp(3)*real(bpsi_f9(3,jbj+i),dp)
           scpr4(9)=scpr4(9)+apsi_f_temp(4)*real(bpsi_f9(4,jbj+i),dp)
           scpr5(9)=scpr5(9)+apsi_f_temp(5)*real(bpsi_f9(5,jbj+i),dp)
           scpr6(9)=scpr6(9)+apsi_f_temp(6)*real(bpsi_f9(6,jbj+i),dp)
           scpr7(9)=scpr7(9)+apsi_f_temp(7)*real(bpsi_f9(7,jbj+i),dp)

           scpr1(10)=scpr1(10)+apsi_f_temp(1)*real(bpsi_f10(1,jbj+i),dp)
           scpr2(10)=scpr2(10)+apsi_f_temp(2)*real(bpsi_f10(2,jbj+i),dp)
           scpr3(10)=scpr3(10)+apsi_f_temp(3)*real(bpsi_f10(3,jbj+i),dp)
           scpr4(10)=scpr4(10)+apsi_f_temp(4)*real(bpsi_f10(4,jbj+i),dp)
           scpr5(10)=scpr5(10)+apsi_f_temp(5)*real(bpsi_f10(5,jbj+i),dp)
           scpr6(10)=scpr6(10)+apsi_f_temp(6)*real(bpsi_f10(6,jbj+i),dp)
           scpr7(10)=scpr7(10)+apsi_f_temp(7)*real(bpsi_f10(7,jbj+i),dp)

           scpr1(11)=scpr1(11)+apsi_f_temp(1)*real(bpsi_f11(1,jbj+i),dp)
           scpr2(11)=scpr2(11)+apsi_f_temp(2)*real(bpsi_f11(2,jbj+i),dp)
           scpr3(11)=scpr3(11)+apsi_f_temp(3)*real(bpsi_f11(3,jbj+i),dp)
           scpr4(11)=scpr4(11)+apsi_f_temp(4)*real(bpsi_f11(4,jbj+i),dp)
           scpr5(11)=scpr5(11)+apsi_f_temp(5)*real(bpsi_f11(5,jbj+i),dp)
           scpr6(11)=scpr6(11)+apsi_f_temp(6)*real(bpsi_f11(6,jbj+i),dp)
           scpr7(11)=scpr7(11)+apsi_f_temp(7)*real(bpsi_f11(7,jbj+i),dp)

           scpr1(12)=scpr1(12)+apsi_f_temp(1)*real(bpsi_f12(1,jbj+i),dp)
           scpr2(12)=scpr2(12)+apsi_f_temp(2)*real(bpsi_f12(2,jbj+i),dp)
           scpr3(12)=scpr3(12)+apsi_f_temp(3)*real(bpsi_f12(3,jbj+i),dp)
           scpr4(12)=scpr4(12)+apsi_f_temp(4)*real(bpsi_f12(4,jbj+i),dp)
           scpr5(12)=scpr5(12)+apsi_f_temp(5)*real(bpsi_f12(5,jbj+i),dp)
           scpr6(12)=scpr6(12)+apsi_f_temp(6)*real(bpsi_f12(6,jbj+i),dp)
           scpr7(12)=scpr7(12)+apsi_f_temp(7)*real(bpsi_f12(7,jbj+i),dp)

           scpr1(13)=scpr1(13)+apsi_f_temp(1)*real(bpsi_f13(1,jbj+i),dp)
           scpr2(13)=scpr2(13)+apsi_f_temp(2)*real(bpsi_f13(2,jbj+i),dp)
           scpr3(13)=scpr3(13)+apsi_f_temp(3)*real(bpsi_f13(3,jbj+i),dp)
           scpr4(13)=scpr4(13)+apsi_f_temp(4)*real(bpsi_f13(4,jbj+i),dp)
           scpr5(13)=scpr5(13)+apsi_f_temp(5)*real(bpsi_f13(5,jbj+i),dp)
           scpr6(13)=scpr6(13)+apsi_f_temp(6)*real(bpsi_f13(6,jbj+i),dp)
           scpr7(13)=scpr7(13)+apsi_f_temp(7)*real(bpsi_f13(7,jbj+i),dp)

        enddo
 !       print *,'length',length,ibseg,scpr1,iaseg0,ja1,jb1
        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_f) exit nonconvex_loop_f !segment is finished  
        iaseg0=iaseg0+1
        jb0=max(keybg_f(1,ibseg),keyag_f_lin(iaseg0))
        if (keyag_f_lin(iaseg0)>jb1) exit nonconvex_loop_f !segment is not covered 
        jbj=jbj+max(jb0-keybg_f(1,ibseg),0)
     end do nonconvex_loop_f
     !disable loop if the end is reached
     if (iaseg0 == maseg_f .and. keybg_f(1,ibseg)> keyag_f_lin(maseg_f)) iaseg0=iaseg0+1

   enddo
!$omp end do !!!implicit barrier 

   scpr0=scpr0+scpr1+scpr2+scpr3+scpr4+scpr5+scpr6+scpr7

!$omp critical 
   scpr=scpr+scpr0
!$omp end critical
  
!$omp end parallel

!!!    call system_clock(ncount2,ncount_rate,ncount_max)
!!!    tel=dble(ncount2-ncount0)/dble(ncount_rate)
!!!    write(97,*) 'wpdot:',tel
!!!    close(97)

END SUBROUTINE wpdot_13

subroutine wpdot_14(  &
     mavctr_c,mavctr_f,maseg_c,maseg_f,keyav_c,keyav_f,keyag_c,keyag_f,apsi_c,apsi_f,  &
     mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keybv_c,keybv_f,keybg_c,keybg_f, &
     bpsi1,bpsi2,bpsi3,bpsi4,bpsi5,bpsi6,bpsi7,bpsi8,bpsi9,bpsi10,bpsi11,bpsi12,bpsi13,bpsi14, &
     bpsi_f1,bpsi_f2,bpsi_f3,bpsi_f4,bpsi_f5,bpsi_f6,bpsi_f7,bpsi_f8,bpsi_f9,bpsi_f10,bpsi_f11,&
     bpsi_f12,bpsi_f13,bpsi_f14, &
     scpr,proj_count)

  use module_base
  implicit none
  integer, intent(in) :: mavctr_c,mavctr_f,maseg_c,maseg_f,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f
  integer, intent(in) :: proj_count
  integer, dimension(maseg_c), intent(in) :: keyav_c
  integer, dimension(mbseg_c), intent(in) :: keybv_c
  integer, dimension(maseg_f), intent(in) :: keyav_f
  integer, dimension(mbseg_f), intent(in) :: keybv_f
  integer, dimension(2,maseg_c), intent(in) :: keyag_c
  integer, dimension(2,mbseg_c), intent(in) :: keybg_c
  integer, dimension(2,maseg_f), intent(in) :: keyag_f
  integer, dimension(2,mbseg_f), intent(in) :: keybg_f
  real(wp), dimension(mavctr_c), intent(in) :: apsi_c
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi1,bpsi2,bpsi3,bpsi4,bpsi5,bpsi6,bpsi7
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi8,bpsi9,bpsi10,bpsi11,bpsi12,bpsi13,bpsi14
  real(wp), dimension(7,mavctr_f), intent(in) :: apsi_f
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f1,bpsi_f2,bpsi_f3,bpsi_f4,bpsi_f5,bpsi_f6,bpsi_f7
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f8,bpsi_f9,bpsi_f10,bpsi_f11,bpsi_f12,bpsi_f13,bpsi_f14
  real(dp), dimension(proj_count),intent(out) :: scpr
  !local variables
  integer :: ibseg,jaj,jb1,jb0,jbj,iaoff,length,i,ja0,ja1
  real(dp),dimension(proj_count) :: scpr1,scpr2,scpr3,scpr4,scpr5,scpr6,scpr7,scpr0
  integer :: iaseg0
  integer, dimension(maseg_c) :: keyag_c_lin !>linear version of second indices of keyag_c
  integer, dimension(maseg_f) :: keyag_f_lin !>linear version of second indices of keyag_f
  real(dp) :: apsi_temp
  real(dp), dimension(7) :: apsi_f_temp
!!!    integer :: ncount0,ncount2,ncount_rate,ncount_max
!!!    real(gp) :: tel

!!!  !dee
!!!    open(unit=97,file='time_wpdot',status='unknown',position='append')
!!!    call system_clock(ncount0,ncount_rate,ncount_max)

  keyag_c_lin = keyag_c(1,:)!speed up access in hunt subroutine by consecutive arrangement in memory
  keyag_f_lin = keyag_f(1,:)!speed up access in hunt subroutine by consecutive arrangement in memory

  scpr=0.0_dp

!$omp parallel default(private) shared(maseg_c,keyav_c,keyag_c,keyag_c_lin) &
!$omp shared(keybg_c,mbseg_c,keybv_c,mbseg_f,maseg_f,apsi_c,bpsi1,bpsi2) &
!$omp shared(bpsi3,bpsi4,bpsi5,bpsi6,bpsi7,bpsi8,bpsi9,bpsi10,bpsi11) &
!$omp shared(bpsi12,bpsi13,bpsi14,bpsi_f1,bpsi_f2,bpsi_f3,bpsi_f4,bpsi_f5) &
!$omp shared(bpsi_f6,bpsi_f7,bpsi_f8,bpsi_f9,bpsi_f10,bpsi_f11,bpsi_f12) &
!$omp shared(bpsi_f13,bpsi_f14,keybv_f,keybg_f,keyag_f,keyag_f_lin,keyav_f) &
!$omp shared(apsi_f,scpr)

!!!!$omp shared (ncount0,ncount2,ncount_rate,ncount_max,tel)

  iaseg0=1 

  scpr0=0.0_dp
  scpr1=0.0_dp
  scpr2=0.0_dp
  scpr3=0.0_dp
  scpr4=0.0_dp
  scpr5=0.0_dp
  scpr6=0.0_dp
  scpr7=0.0_dp
!coarse part. Loop on the projectors segments
!$omp do schedule(static)
   do ibseg=1,mbseg_c
     jbj=keybv_c(ibseg)
     jb0=keybg_c(1,ibseg) !starting point of projector segment
     jb1=keybg_c(2,ibseg) !ending point of projector segment
!     print *,'huntenter',ibseg,jb0,jb1
 
     !find the starting point of the wavefunction segment
     !warning: hunt is assuming that the variable is always found
     !if it is not, iaseg0 is put to maseg + 1 so that the loop is disabled
     call hunt1(.true.,keyag_c_lin,maseg_c,keybg_c(1,ibseg),iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     !now pass through all the wavefunction segments until the end of the segment is 
     !still contained in projector segment
     nonconvex_loop_c: do while(iaseg0 <= maseg_c)
!     print *,'huntexit',iaseg0,maseg_c,keyag_c_lin(iaseg0),keyag_c(2,iaseg0)

        !length = jb1-jb0
        !iaoff = jb0-keyag_c_lin(iaseg0)!jb0-ja0

        ja0=keyag_c_lin(iaseg0)
        ja1=min(jb1,keyag_c(2,iaseg0)) 
        length = ja1-jb0
        iaoff = max(jb0-ja0,0) !no offset if we are already inside

        jaj=keyav_c(iaseg0)
        do i=0,length
           apsi_temp=real(apsi_c(jaj+iaoff+i),dp)

           scpr0(1)=scpr0(1)+apsi_temp *real(bpsi1(jbj+i),dp)
           scpr0(2)=scpr0(2)+apsi_temp *real(bpsi2(jbj+i),dp)
           scpr0(3)=scpr0(3)+apsi_temp *real(bpsi3(jbj+i),dp)
           scpr0(4)=scpr0(4)+apsi_temp *real(bpsi4(jbj+i),dp)
           scpr0(5)=scpr0(5)+apsi_temp *real(bpsi5(jbj+i),dp)
           scpr0(6)=scpr0(6)+apsi_temp *real(bpsi6(jbj+i),dp)
           scpr0(7)=scpr0(7)+apsi_temp *real(bpsi7(jbj+i),dp)
           scpr0(8)=scpr0(8)+apsi_temp *real(bpsi8(jbj+i),dp)
           scpr0(9)=scpr0(9)+apsi_temp *real(bpsi9(jbj+i),dp)
           scpr0(10)=scpr0(10)+apsi_temp *real(bpsi10(jbj+i),dp)
           scpr0(11)=scpr0(11)+apsi_temp *real(bpsi11(jbj+i),dp)
           scpr0(12)=scpr0(12)+apsi_temp *real(bpsi12(jbj+i),dp)
           scpr0(13)=scpr0(13)+apsi_temp *real(bpsi13(jbj+i),dp)
           scpr0(14)=scpr0(14)+apsi_temp *real(bpsi14(jbj+i),dp)
        enddo
 !       print *,'length',length,ibseg,scpr0,iaseg0,ja1,jb1

        !print *,'ibseg,mbseg_c,iaseg0,maseg_c',ibseg,mbseg_c,iaseg0,maseg_c
        !print '(a,6(i8),1pe25.17)','ja0,ja1t,ja1,jb0,jb1',&
        !     ibseg,ja0,keyag_c(2,iaseg0),ja1,jb0,jb1,scpr0
        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_c) exit nonconvex_loop_c !segment is finished
        iaseg0=iaseg0+1
        jb0=max(keybg_c(1,ibseg),keyag_c_lin(iaseg0))
        if (keyag_c_lin(iaseg0)>jb1) exit nonconvex_loop_c !segment is not covered
        jbj=jbj+max(jb0-keybg_c(1,ibseg),0)
     end do nonconvex_loop_c
     !disable loop if the end is reached
     if (iaseg0 == maseg_c .and. keybg_c(1,ibseg)> keyag_c_lin(maseg_c)) iaseg0=iaseg0+1


   enddo
!stop
!$omp end do nowait

! fine part

iaseg0=1

!$omp do schedule(static)
   do ibseg=1,mbseg_f
     jbj=keybv_f(ibseg)
     jb0=keybg_f(1,ibseg)
     jb1=keybg_f(2,ibseg)
!    print *,'huntenter',ibseg,jb0,jb1
     call hunt1(.true.,keyag_f_lin,maseg_f,keybg_f(1,ibseg),iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     nonconvex_loop_f: do while(iaseg0 <= maseg_f)
!!$     length = jb1-jb0
!!$     iaoff = jb0-keyag_f_lin(iaseg0)

        ja0=keyag_f_lin(iaseg0) !still doubts about copying in automatic array
        ja1=min(jb1,keyag_f(2,iaseg0)) 
        length = ja1-jb0
        iaoff = max(jb0-ja0,0) !no offset if we are already inside

        jaj=keyav_f(iaseg0)
        do i=0,length
           apsi_f_temp(1)=real(apsi_f(1,jaj+iaoff+i),dp)
           apsi_f_temp(2)=real(apsi_f(2,jaj+iaoff+i),dp)
           apsi_f_temp(3)=real(apsi_f(3,jaj+iaoff+i),dp)
           apsi_f_temp(4)=real(apsi_f(4,jaj+iaoff+i),dp)
           apsi_f_temp(5)=real(apsi_f(5,jaj+iaoff+i),dp)
           apsi_f_temp(6)=real(apsi_f(6,jaj+iaoff+i),dp)
           apsi_f_temp(7)=real(apsi_f(7,jaj+iaoff+i),dp)

           scpr1(1)=scpr1(1)+apsi_f_temp(1)*real(bpsi_f1(1,jbj+i),dp)
           scpr2(1)=scpr2(1)+apsi_f_temp(2)*real(bpsi_f1(2,jbj+i),dp)
           scpr3(1)=scpr3(1)+apsi_f_temp(3)*real(bpsi_f1(3,jbj+i),dp)
           scpr4(1)=scpr4(1)+apsi_f_temp(4)*real(bpsi_f1(4,jbj+i),dp)
           scpr5(1)=scpr5(1)+apsi_f_temp(5)*real(bpsi_f1(5,jbj+i),dp)
           scpr6(1)=scpr6(1)+apsi_f_temp(6)*real(bpsi_f1(6,jbj+i),dp)
           scpr7(1)=scpr7(1)+apsi_f_temp(7)*real(bpsi_f1(7,jbj+i),dp)

           scpr1(2)=scpr1(2)+apsi_f_temp(1)*real(bpsi_f2(1,jbj+i),dp)
           scpr2(2)=scpr2(2)+apsi_f_temp(2)*real(bpsi_f2(2,jbj+i),dp)
           scpr3(2)=scpr3(2)+apsi_f_temp(3)*real(bpsi_f2(3,jbj+i),dp)
           scpr4(2)=scpr4(2)+apsi_f_temp(4)*real(bpsi_f2(4,jbj+i),dp)
           scpr5(2)=scpr5(2)+apsi_f_temp(5)*real(bpsi_f2(5,jbj+i),dp)
           scpr6(2)=scpr6(2)+apsi_f_temp(6)*real(bpsi_f2(6,jbj+i),dp)
           scpr7(2)=scpr7(2)+apsi_f_temp(7)*real(bpsi_f2(7,jbj+i),dp)

           scpr1(3)=scpr1(3)+apsi_f_temp(1)*real(bpsi_f3(1,jbj+i),dp)
           scpr2(3)=scpr2(3)+apsi_f_temp(2)*real(bpsi_f3(2,jbj+i),dp)
           scpr3(3)=scpr3(3)+apsi_f_temp(3)*real(bpsi_f3(3,jbj+i),dp)
           scpr4(3)=scpr4(3)+apsi_f_temp(4)*real(bpsi_f3(4,jbj+i),dp)
           scpr5(3)=scpr5(3)+apsi_f_temp(5)*real(bpsi_f3(5,jbj+i),dp)
           scpr6(3)=scpr6(3)+apsi_f_temp(6)*real(bpsi_f3(6,jbj+i),dp)
           scpr7(3)=scpr7(3)+apsi_f_temp(7)*real(bpsi_f3(7,jbj+i),dp)

           scpr1(4)=scpr1(4)+apsi_f_temp(1)*real(bpsi_f4(1,jbj+i),dp)
           scpr2(4)=scpr2(4)+apsi_f_temp(2)*real(bpsi_f4(2,jbj+i),dp)
           scpr3(4)=scpr3(4)+apsi_f_temp(3)*real(bpsi_f4(3,jbj+i),dp)
           scpr4(4)=scpr4(4)+apsi_f_temp(4)*real(bpsi_f4(4,jbj+i),dp)
           scpr5(4)=scpr5(4)+apsi_f_temp(5)*real(bpsi_f4(5,jbj+i),dp)
           scpr6(4)=scpr6(4)+apsi_f_temp(6)*real(bpsi_f4(6,jbj+i),dp)
           scpr7(4)=scpr7(4)+apsi_f_temp(7)*real(bpsi_f4(7,jbj+i),dp)

           scpr1(5)=scpr1(5)+apsi_f_temp(1)*real(bpsi_f5(1,jbj+i),dp)
           scpr2(5)=scpr2(5)+apsi_f_temp(2)*real(bpsi_f5(2,jbj+i),dp)
           scpr3(5)=scpr3(5)+apsi_f_temp(3)*real(bpsi_f5(3,jbj+i),dp)
           scpr4(5)=scpr4(5)+apsi_f_temp(4)*real(bpsi_f5(4,jbj+i),dp)
           scpr5(5)=scpr5(5)+apsi_f_temp(5)*real(bpsi_f5(5,jbj+i),dp)
           scpr6(5)=scpr6(5)+apsi_f_temp(6)*real(bpsi_f5(6,jbj+i),dp)
           scpr7(5)=scpr7(5)+apsi_f_temp(7)*real(bpsi_f5(7,jbj+i),dp)

           scpr1(6)=scpr1(6)+apsi_f_temp(1)*real(bpsi_f6(1,jbj+i),dp)
           scpr2(6)=scpr2(6)+apsi_f_temp(2)*real(bpsi_f6(2,jbj+i),dp)
           scpr3(6)=scpr3(6)+apsi_f_temp(3)*real(bpsi_f6(3,jbj+i),dp)
           scpr4(6)=scpr4(6)+apsi_f_temp(4)*real(bpsi_f6(4,jbj+i),dp)
           scpr5(6)=scpr5(6)+apsi_f_temp(5)*real(bpsi_f6(5,jbj+i),dp)
           scpr6(6)=scpr6(6)+apsi_f_temp(6)*real(bpsi_f6(6,jbj+i),dp)
           scpr7(6)=scpr7(6)+apsi_f_temp(7)*real(bpsi_f6(7,jbj+i),dp)

           scpr1(7)=scpr1(7)+apsi_f_temp(1)*real(bpsi_f7(1,jbj+i),dp)
           scpr2(7)=scpr2(7)+apsi_f_temp(2)*real(bpsi_f7(2,jbj+i),dp)
           scpr3(7)=scpr3(7)+apsi_f_temp(3)*real(bpsi_f7(3,jbj+i),dp)
           scpr4(7)=scpr4(7)+apsi_f_temp(4)*real(bpsi_f7(4,jbj+i),dp)
           scpr5(7)=scpr5(7)+apsi_f_temp(5)*real(bpsi_f7(5,jbj+i),dp)
           scpr6(7)=scpr6(7)+apsi_f_temp(6)*real(bpsi_f7(6,jbj+i),dp)
           scpr7(7)=scpr7(7)+apsi_f_temp(7)*real(bpsi_f7(7,jbj+i),dp)

           scpr1(8)=scpr1(8)+apsi_f_temp(1)*real(bpsi_f8(1,jbj+i),dp)
           scpr2(8)=scpr2(8)+apsi_f_temp(2)*real(bpsi_f8(2,jbj+i),dp)
           scpr3(8)=scpr3(8)+apsi_f_temp(3)*real(bpsi_f8(3,jbj+i),dp)
           scpr4(8)=scpr4(8)+apsi_f_temp(4)*real(bpsi_f8(4,jbj+i),dp)
           scpr5(8)=scpr5(8)+apsi_f_temp(5)*real(bpsi_f8(5,jbj+i),dp)
           scpr6(8)=scpr6(8)+apsi_f_temp(6)*real(bpsi_f8(6,jbj+i),dp)
           scpr7(8)=scpr7(8)+apsi_f_temp(7)*real(bpsi_f8(7,jbj+i),dp)

           scpr1(9)=scpr1(9)+apsi_f_temp(1)*real(bpsi_f9(1,jbj+i),dp)
           scpr2(9)=scpr2(9)+apsi_f_temp(2)*real(bpsi_f9(2,jbj+i),dp)
           scpr3(9)=scpr3(9)+apsi_f_temp(3)*real(bpsi_f9(3,jbj+i),dp)
           scpr4(9)=scpr4(9)+apsi_f_temp(4)*real(bpsi_f9(4,jbj+i),dp)
           scpr5(9)=scpr5(9)+apsi_f_temp(5)*real(bpsi_f9(5,jbj+i),dp)
           scpr6(9)=scpr6(9)+apsi_f_temp(6)*real(bpsi_f9(6,jbj+i),dp)
           scpr7(9)=scpr7(9)+apsi_f_temp(7)*real(bpsi_f9(7,jbj+i),dp)

           scpr1(10)=scpr1(10)+apsi_f_temp(1)*real(bpsi_f10(1,jbj+i),dp)
           scpr2(10)=scpr2(10)+apsi_f_temp(2)*real(bpsi_f10(2,jbj+i),dp)
           scpr3(10)=scpr3(10)+apsi_f_temp(3)*real(bpsi_f10(3,jbj+i),dp)
           scpr4(10)=scpr4(10)+apsi_f_temp(4)*real(bpsi_f10(4,jbj+i),dp)
           scpr5(10)=scpr5(10)+apsi_f_temp(5)*real(bpsi_f10(5,jbj+i),dp)
           scpr6(10)=scpr6(10)+apsi_f_temp(6)*real(bpsi_f10(6,jbj+i),dp)
           scpr7(10)=scpr7(10)+apsi_f_temp(7)*real(bpsi_f10(7,jbj+i),dp)

           scpr1(11)=scpr1(11)+apsi_f_temp(1)*real(bpsi_f11(1,jbj+i),dp)
           scpr2(11)=scpr2(11)+apsi_f_temp(2)*real(bpsi_f11(2,jbj+i),dp)
           scpr3(11)=scpr3(11)+apsi_f_temp(3)*real(bpsi_f11(3,jbj+i),dp)
           scpr4(11)=scpr4(11)+apsi_f_temp(4)*real(bpsi_f11(4,jbj+i),dp)
           scpr5(11)=scpr5(11)+apsi_f_temp(5)*real(bpsi_f11(5,jbj+i),dp)
           scpr6(11)=scpr6(11)+apsi_f_temp(6)*real(bpsi_f11(6,jbj+i),dp)
           scpr7(11)=scpr7(11)+apsi_f_temp(7)*real(bpsi_f11(7,jbj+i),dp)

           scpr1(12)=scpr1(12)+apsi_f_temp(1)*real(bpsi_f12(1,jbj+i),dp)
           scpr2(12)=scpr2(12)+apsi_f_temp(2)*real(bpsi_f12(2,jbj+i),dp)
           scpr3(12)=scpr3(12)+apsi_f_temp(3)*real(bpsi_f12(3,jbj+i),dp)
           scpr4(12)=scpr4(12)+apsi_f_temp(4)*real(bpsi_f12(4,jbj+i),dp)
           scpr5(12)=scpr5(12)+apsi_f_temp(5)*real(bpsi_f12(5,jbj+i),dp)
           scpr6(12)=scpr6(12)+apsi_f_temp(6)*real(bpsi_f12(6,jbj+i),dp)
           scpr7(12)=scpr7(12)+apsi_f_temp(7)*real(bpsi_f12(7,jbj+i),dp)

           scpr1(13)=scpr1(13)+apsi_f_temp(1)*real(bpsi_f13(1,jbj+i),dp)
           scpr2(13)=scpr2(13)+apsi_f_temp(2)*real(bpsi_f13(2,jbj+i),dp)
           scpr3(13)=scpr3(13)+apsi_f_temp(3)*real(bpsi_f13(3,jbj+i),dp)
           scpr4(13)=scpr4(13)+apsi_f_temp(4)*real(bpsi_f13(4,jbj+i),dp)
           scpr5(13)=scpr5(13)+apsi_f_temp(5)*real(bpsi_f13(5,jbj+i),dp)
           scpr6(13)=scpr6(13)+apsi_f_temp(6)*real(bpsi_f13(6,jbj+i),dp)
           scpr7(13)=scpr7(13)+apsi_f_temp(7)*real(bpsi_f13(7,jbj+i),dp)

           scpr1(14)=scpr1(14)+apsi_f_temp(1)*real(bpsi_f14(1,jbj+i),dp)
           scpr2(14)=scpr2(14)+apsi_f_temp(2)*real(bpsi_f14(2,jbj+i),dp)
           scpr3(14)=scpr3(14)+apsi_f_temp(3)*real(bpsi_f14(3,jbj+i),dp)
           scpr4(14)=scpr4(14)+apsi_f_temp(4)*real(bpsi_f14(4,jbj+i),dp)
           scpr5(14)=scpr5(14)+apsi_f_temp(5)*real(bpsi_f14(5,jbj+i),dp)
           scpr6(14)=scpr6(14)+apsi_f_temp(6)*real(bpsi_f14(6,jbj+i),dp)
           scpr7(14)=scpr7(14)+apsi_f_temp(7)*real(bpsi_f14(7,jbj+i),dp)
        enddo
 !       print *,'length',length,ibseg,scpr1,iaseg0,ja1,jb1
        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_f) exit nonconvex_loop_f !segment is finished  
        iaseg0=iaseg0+1
        jb0=max(keybg_f(1,ibseg),keyag_f_lin(iaseg0))
        if (keyag_f_lin(iaseg0)>jb1) exit nonconvex_loop_f !segment is not covered 
        jbj=jbj+max(jb0-keybg_f(1,ibseg),0)
     end do nonconvex_loop_f
     !disable loop if the end is reached
     if (iaseg0 == maseg_f .and. keybg_f(1,ibseg)> keyag_f_lin(maseg_f)) iaseg0=iaseg0+1

   enddo
!$omp end do !!!implicit barrier 

   scpr0=scpr0+scpr1+scpr2+scpr3+scpr4+scpr5+scpr6+scpr7

!$omp critical 
   scpr=scpr+scpr0
!$omp end critical
  
!$omp end parallel

!!!    call system_clock(ncount2,ncount_rate,ncount_max)
!!!    tel=dble(ncount2-ncount0)/dble(ncount_rate)
!!!    write(97,*) 'wpdot:',tel
!!!    close(97)

END SUBROUTINE wpdot_14

subroutine wpdot_18(  &
     mavctr_c,mavctr_f,maseg_c,maseg_f,keyav_c,keyav_f,keyag_c,keyag_f,apsi_c,apsi_f,  &
     mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keybv_c,keybv_f,keybg_c,keybg_f, &
     bpsi1,bpsi2,bpsi3,bpsi4,bpsi5,bpsi6,bpsi7,bpsi8,bpsi9,bpsi10,bpsi11,bpsi12,bpsi13,bpsi14, &
     bpsi15,bpsi16,bpsi17,bpsi18,bpsi_f1,bpsi_f2,bpsi_f3,bpsi_f4,bpsi_f5,bpsi_f6,bpsi_f7, &
     bpsi_f8,bpsi_f9,bpsi_f10,bpsi_f11,bpsi_f12,bpsi_f13,bpsi_f14,bpsi_f15,bpsi_f16,&
     bpsi_f17,bpsi_f18, &
     scpr,proj_count)

  use module_base
  implicit none
  integer, intent(in) :: mavctr_c,mavctr_f,maseg_c,maseg_f,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f
  integer, intent(in) :: proj_count
  integer, dimension(maseg_c), intent(in) :: keyav_c
  integer, dimension(mbseg_c), intent(in) :: keybv_c
  integer, dimension(maseg_f), intent(in) :: keyav_f
  integer, dimension(mbseg_f), intent(in) :: keybv_f
  integer, dimension(2,maseg_c), intent(in) :: keyag_c
  integer, dimension(2,mbseg_c), intent(in) :: keybg_c
  integer, dimension(2,maseg_f), intent(in) :: keyag_f
  integer, dimension(2,mbseg_f), intent(in) :: keybg_f
  real(wp), dimension(mavctr_c), intent(in) :: apsi_c
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi1,bpsi2,bpsi3,bpsi4,bpsi5,bpsi6,bpsi7
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi8,bpsi9,bpsi10,bpsi11,bpsi12,bpsi13,bpsi14
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi15,bpsi16,bpsi17,bpsi18
  real(wp), dimension(7,mavctr_f), intent(in) :: apsi_f
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f1,bpsi_f2,bpsi_f3,bpsi_f4,bpsi_f5,bpsi_f6,bpsi_f7
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f8,bpsi_f9,bpsi_f10,bpsi_f11,bpsi_f12,bpsi_f13,bpsi_f14
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f15,bpsi_f16,bpsi_f17,bpsi_f18
  real(dp), dimension(proj_count),intent(out) :: scpr
  !local variables
  integer :: ibseg,jaj,jb1,jb0,jbj,iaoff,length,i,ja0,ja1
  real(dp),dimension(proj_count) :: scpr1,scpr2,scpr3,scpr4,scpr5,scpr6,scpr7,scpr0
  integer :: iaseg0
  integer, dimension(maseg_c) :: keyag_c_lin !>linear version of second indices of keyag_c
  integer, dimension(maseg_f) :: keyag_f_lin !>linear version of second indices of keyag_f
  real(dp) :: apsi_temp
  real(dp), dimension(7) :: apsi_f_temp
!!!    integer :: ncount0,ncount2,ncount_rate,ncount_max
!!!    real(gp) :: tel

!!!  !dee
!!!    open(unit=97,file='time_wpdot',status='unknown',position='append')
!!!    call system_clock(ncount0,ncount_rate,ncount_max)

  keyag_c_lin = keyag_c(1,:)!speed up access in hunt subroutine by consecutive arrangement in memory
  keyag_f_lin = keyag_f(1,:)!speed up access in hunt subroutine by consecutive arrangement in memory

  scpr=0.0_dp

!$omp parallel default(private) shared(maseg_c,keyav_c,keyag_c,keyag_c_lin) &
!$omp shared(keybg_c,mbseg_c,keybv_c,mbseg_f,maseg_f,apsi_c,bpsi1,bpsi2) &
!$omp shared(bpsi3,bpsi4,bpsi5,bpsi6,bpsi7,bpsi8,bpsi9,bpsi10,bpsi11) &
!$omp shared(bpsi12,bpsi13,bpsi14,bpsi15,bpsi16,bpsi17,bpsi18) &
!$omp shared(bpsi_f1,bpsi_f2,bpsi_f3,bpsi_f4,bpsi_f5,bpsi_f6,bpsi_f7) &
!$omp shared(bpsi_f8,bpsi_f9,bpsi_f10,bpsi_f11,bpsi_f12,bpsi_f13,bpsi_f14) &
!$omp shared(bpsi_f15,bpsi_f16,bpsi_f17,bpsi_f18,keybv_f,keybg_f) &
!$omp shared(keyag_f,keyag_f_lin,keyav_f,apsi_f,scpr)

!!!!$omp shared (ncount0,ncount2,ncount_rate,ncount_max,tel)

  iaseg0=1 

  scpr0=0.0_dp
  scpr1=0.0_dp
  scpr2=0.0_dp
  scpr3=0.0_dp
  scpr4=0.0_dp
  scpr5=0.0_dp
  scpr6=0.0_dp
  scpr7=0.0_dp
!coarse part. Loop on the projectors segments
!$omp do schedule(static)
   do ibseg=1,mbseg_c
     jbj=keybv_c(ibseg)
     jb0=keybg_c(1,ibseg) !starting point of projector segment
     jb1=keybg_c(2,ibseg) !ending point of projector segment
!     print *,'huntenter',ibseg,jb0,jb1
 
     !find the starting point of the wavefunction segment
     !warning: hunt is assuming that the variable is always found
     !if it is not, iaseg0 is put to maseg + 1 so that the loop is disabled
     call hunt1(.true.,keyag_c_lin,maseg_c,keybg_c(1,ibseg),iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     !now pass through all the wavefunction segments until the end of the segment is 
     !still contained in projector segment
     nonconvex_loop_c: do while(iaseg0 <= maseg_c)
!     print *,'huntexit',iaseg0,maseg_c,keyag_c_lin(iaseg0),keyag_c(2,iaseg0)

        !length = jb1-jb0
        !iaoff = jb0-keyag_c_lin(iaseg0)!jb0-ja0

        ja0=keyag_c_lin(iaseg0)
        ja1=min(jb1,keyag_c(2,iaseg0)) 
        length = ja1-jb0
        iaoff = max(jb0-ja0,0) !no offset if we are already inside

        jaj=keyav_c(iaseg0)
        do i=0,length
           apsi_temp=real(apsi_c(jaj+iaoff+i),dp)

           scpr0(1)=scpr0(1)+apsi_temp *real(bpsi1(jbj+i),dp)
           scpr0(2)=scpr0(2)+apsi_temp *real(bpsi2(jbj+i),dp)
           scpr0(3)=scpr0(3)+apsi_temp *real(bpsi3(jbj+i),dp)
           scpr0(4)=scpr0(4)+apsi_temp *real(bpsi4(jbj+i),dp)
           scpr0(5)=scpr0(5)+apsi_temp *real(bpsi5(jbj+i),dp)
           scpr0(6)=scpr0(6)+apsi_temp *real(bpsi6(jbj+i),dp)
           scpr0(7)=scpr0(7)+apsi_temp *real(bpsi7(jbj+i),dp)
           scpr0(8)=scpr0(8)+apsi_temp *real(bpsi8(jbj+i),dp)
           scpr0(9)=scpr0(9)+apsi_temp *real(bpsi9(jbj+i),dp)
           scpr0(10)=scpr0(10)+apsi_temp *real(bpsi10(jbj+i),dp)
           scpr0(11)=scpr0(11)+apsi_temp *real(bpsi11(jbj+i),dp)
           scpr0(12)=scpr0(12)+apsi_temp *real(bpsi12(jbj+i),dp)
           scpr0(13)=scpr0(13)+apsi_temp *real(bpsi13(jbj+i),dp)
           scpr0(14)=scpr0(14)+apsi_temp *real(bpsi14(jbj+i),dp)
           scpr0(15)=scpr0(15)+apsi_temp *real(bpsi15(jbj+i),dp)
           scpr0(16)=scpr0(16)+apsi_temp *real(bpsi16(jbj+i),dp)
           scpr0(17)=scpr0(17)+apsi_temp *real(bpsi17(jbj+i),dp)
           scpr0(18)=scpr0(18)+apsi_temp *real(bpsi18(jbj+i),dp)
        enddo
 !       print *,'length',length,ibseg,scpr0,iaseg0,ja1,jb1

        !print *,'ibseg,mbseg_c,iaseg0,maseg_c',ibseg,mbseg_c,iaseg0,maseg_c
        !print '(a,6(i8),1pe25.17)','ja0,ja1t,ja1,jb0,jb1',&
        !     ibseg,ja0,keyag_c(2,iaseg0),ja1,jb0,jb1,scpr0
        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_c) exit nonconvex_loop_c !segment is finished
        iaseg0=iaseg0+1
        jb0=max(keybg_c(1,ibseg),keyag_c_lin(iaseg0))
        if (keyag_c_lin(iaseg0)>jb1) exit nonconvex_loop_c !segment is not covered
        jbj=jbj+max(jb0-keybg_c(1,ibseg),0)
     end do nonconvex_loop_c
     !disable loop if the end is reached
     if (iaseg0 == maseg_c .and. keybg_c(1,ibseg)> keyag_c_lin(maseg_c)) iaseg0=iaseg0+1


   enddo
!stop
!$omp end do nowait

! fine part

iaseg0=1

!$omp do schedule(static)
   do ibseg=1,mbseg_f
     jbj=keybv_f(ibseg)
     jb0=keybg_f(1,ibseg)
     jb1=keybg_f(2,ibseg)
!    print *,'huntenter',ibseg,jb0,jb1
     call hunt1(.true.,keyag_f_lin,maseg_f,keybg_f(1,ibseg),iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     nonconvex_loop_f: do while(iaseg0 <= maseg_f)
!!$     length = jb1-jb0
!!$     iaoff = jb0-keyag_f_lin(iaseg0)

        ja0=keyag_f_lin(iaseg0) !still doubts about copying in automatic array
        ja1=min(jb1,keyag_f(2,iaseg0)) 
        length = ja1-jb0
        iaoff = max(jb0-ja0,0) !no offset if we are already inside

        jaj=keyav_f(iaseg0)
        do i=0,length
           apsi_f_temp(1)=real(apsi_f(1,jaj+iaoff+i),dp)
           apsi_f_temp(2)=real(apsi_f(2,jaj+iaoff+i),dp)
           apsi_f_temp(3)=real(apsi_f(3,jaj+iaoff+i),dp)
           apsi_f_temp(4)=real(apsi_f(4,jaj+iaoff+i),dp)
           apsi_f_temp(5)=real(apsi_f(5,jaj+iaoff+i),dp)
           apsi_f_temp(6)=real(apsi_f(6,jaj+iaoff+i),dp)
           apsi_f_temp(7)=real(apsi_f(7,jaj+iaoff+i),dp)

           scpr1(1)=scpr1(1)+apsi_f_temp(1)*real(bpsi_f1(1,jbj+i),dp)
           scpr2(1)=scpr2(1)+apsi_f_temp(2)*real(bpsi_f1(2,jbj+i),dp)
           scpr3(1)=scpr3(1)+apsi_f_temp(3)*real(bpsi_f1(3,jbj+i),dp)
           scpr4(1)=scpr4(1)+apsi_f_temp(4)*real(bpsi_f1(4,jbj+i),dp)
           scpr5(1)=scpr5(1)+apsi_f_temp(5)*real(bpsi_f1(5,jbj+i),dp)
           scpr6(1)=scpr6(1)+apsi_f_temp(6)*real(bpsi_f1(6,jbj+i),dp)
           scpr7(1)=scpr7(1)+apsi_f_temp(7)*real(bpsi_f1(7,jbj+i),dp)

           scpr1(2)=scpr1(2)+apsi_f_temp(1)*real(bpsi_f2(1,jbj+i),dp)
           scpr2(2)=scpr2(2)+apsi_f_temp(2)*real(bpsi_f2(2,jbj+i),dp)
           scpr3(2)=scpr3(2)+apsi_f_temp(3)*real(bpsi_f2(3,jbj+i),dp)
           scpr4(2)=scpr4(2)+apsi_f_temp(4)*real(bpsi_f2(4,jbj+i),dp)
           scpr5(2)=scpr5(2)+apsi_f_temp(5)*real(bpsi_f2(5,jbj+i),dp)
           scpr6(2)=scpr6(2)+apsi_f_temp(6)*real(bpsi_f2(6,jbj+i),dp)
           scpr7(2)=scpr7(2)+apsi_f_temp(7)*real(bpsi_f2(7,jbj+i),dp)

           scpr1(3)=scpr1(3)+apsi_f_temp(1)*real(bpsi_f3(1,jbj+i),dp)
           scpr2(3)=scpr2(3)+apsi_f_temp(2)*real(bpsi_f3(2,jbj+i),dp)
           scpr3(3)=scpr3(3)+apsi_f_temp(3)*real(bpsi_f3(3,jbj+i),dp)
           scpr4(3)=scpr4(3)+apsi_f_temp(4)*real(bpsi_f3(4,jbj+i),dp)
           scpr5(3)=scpr5(3)+apsi_f_temp(5)*real(bpsi_f3(5,jbj+i),dp)
           scpr6(3)=scpr6(3)+apsi_f_temp(6)*real(bpsi_f3(6,jbj+i),dp)
           scpr7(3)=scpr7(3)+apsi_f_temp(7)*real(bpsi_f3(7,jbj+i),dp)

           scpr1(4)=scpr1(4)+apsi_f_temp(1)*real(bpsi_f4(1,jbj+i),dp)
           scpr2(4)=scpr2(4)+apsi_f_temp(2)*real(bpsi_f4(2,jbj+i),dp)
           scpr3(4)=scpr3(4)+apsi_f_temp(3)*real(bpsi_f4(3,jbj+i),dp)
           scpr4(4)=scpr4(4)+apsi_f_temp(4)*real(bpsi_f4(4,jbj+i),dp)
           scpr5(4)=scpr5(4)+apsi_f_temp(5)*real(bpsi_f4(5,jbj+i),dp)
           scpr6(4)=scpr6(4)+apsi_f_temp(6)*real(bpsi_f4(6,jbj+i),dp)
           scpr7(4)=scpr7(4)+apsi_f_temp(7)*real(bpsi_f4(7,jbj+i),dp)

           scpr1(5)=scpr1(5)+apsi_f_temp(1)*real(bpsi_f5(1,jbj+i),dp)
           scpr2(5)=scpr2(5)+apsi_f_temp(2)*real(bpsi_f5(2,jbj+i),dp)
           scpr3(5)=scpr3(5)+apsi_f_temp(3)*real(bpsi_f5(3,jbj+i),dp)
           scpr4(5)=scpr4(5)+apsi_f_temp(4)*real(bpsi_f5(4,jbj+i),dp)
           scpr5(5)=scpr5(5)+apsi_f_temp(5)*real(bpsi_f5(5,jbj+i),dp)
           scpr6(5)=scpr6(5)+apsi_f_temp(6)*real(bpsi_f5(6,jbj+i),dp)
           scpr7(5)=scpr7(5)+apsi_f_temp(7)*real(bpsi_f5(7,jbj+i),dp)

           scpr1(6)=scpr1(6)+apsi_f_temp(1)*real(bpsi_f6(1,jbj+i),dp)
           scpr2(6)=scpr2(6)+apsi_f_temp(2)*real(bpsi_f6(2,jbj+i),dp)
           scpr3(6)=scpr3(6)+apsi_f_temp(3)*real(bpsi_f6(3,jbj+i),dp)
           scpr4(6)=scpr4(6)+apsi_f_temp(4)*real(bpsi_f6(4,jbj+i),dp)
           scpr5(6)=scpr5(6)+apsi_f_temp(5)*real(bpsi_f6(5,jbj+i),dp)
           scpr6(6)=scpr6(6)+apsi_f_temp(6)*real(bpsi_f6(6,jbj+i),dp)
           scpr7(6)=scpr7(6)+apsi_f_temp(7)*real(bpsi_f6(7,jbj+i),dp)

           scpr1(7)=scpr1(7)+apsi_f_temp(1)*real(bpsi_f7(1,jbj+i),dp)
           scpr2(7)=scpr2(7)+apsi_f_temp(2)*real(bpsi_f7(2,jbj+i),dp)
           scpr3(7)=scpr3(7)+apsi_f_temp(3)*real(bpsi_f7(3,jbj+i),dp)
           scpr4(7)=scpr4(7)+apsi_f_temp(4)*real(bpsi_f7(4,jbj+i),dp)
           scpr5(7)=scpr5(7)+apsi_f_temp(5)*real(bpsi_f7(5,jbj+i),dp)
           scpr6(7)=scpr6(7)+apsi_f_temp(6)*real(bpsi_f7(6,jbj+i),dp)
           scpr7(7)=scpr7(7)+apsi_f_temp(7)*real(bpsi_f7(7,jbj+i),dp)

           scpr1(8)=scpr1(8)+apsi_f_temp(1)*real(bpsi_f8(1,jbj+i),dp)
           scpr2(8)=scpr2(8)+apsi_f_temp(2)*real(bpsi_f8(2,jbj+i),dp)
           scpr3(8)=scpr3(8)+apsi_f_temp(3)*real(bpsi_f8(3,jbj+i),dp)
           scpr4(8)=scpr4(8)+apsi_f_temp(4)*real(bpsi_f8(4,jbj+i),dp)
           scpr5(8)=scpr5(8)+apsi_f_temp(5)*real(bpsi_f8(5,jbj+i),dp)
           scpr6(8)=scpr6(8)+apsi_f_temp(6)*real(bpsi_f8(6,jbj+i),dp)
           scpr7(8)=scpr7(8)+apsi_f_temp(7)*real(bpsi_f8(7,jbj+i),dp)

           scpr1(9)=scpr1(9)+apsi_f_temp(1)*real(bpsi_f9(1,jbj+i),dp)
           scpr2(9)=scpr2(9)+apsi_f_temp(2)*real(bpsi_f9(2,jbj+i),dp)
           scpr3(9)=scpr3(9)+apsi_f_temp(3)*real(bpsi_f9(3,jbj+i),dp)
           scpr4(9)=scpr4(9)+apsi_f_temp(4)*real(bpsi_f9(4,jbj+i),dp)
           scpr5(9)=scpr5(9)+apsi_f_temp(5)*real(bpsi_f9(5,jbj+i),dp)
           scpr6(9)=scpr6(9)+apsi_f_temp(6)*real(bpsi_f9(6,jbj+i),dp)
           scpr7(9)=scpr7(9)+apsi_f_temp(7)*real(bpsi_f9(7,jbj+i),dp)

           scpr1(10)=scpr1(10)+apsi_f_temp(1)*real(bpsi_f10(1,jbj+i),dp)
           scpr2(10)=scpr2(10)+apsi_f_temp(2)*real(bpsi_f10(2,jbj+i),dp)
           scpr3(10)=scpr3(10)+apsi_f_temp(3)*real(bpsi_f10(3,jbj+i),dp)
           scpr4(10)=scpr4(10)+apsi_f_temp(4)*real(bpsi_f10(4,jbj+i),dp)
           scpr5(10)=scpr5(10)+apsi_f_temp(5)*real(bpsi_f10(5,jbj+i),dp)
           scpr6(10)=scpr6(10)+apsi_f_temp(6)*real(bpsi_f10(6,jbj+i),dp)
           scpr7(10)=scpr7(10)+apsi_f_temp(7)*real(bpsi_f10(7,jbj+i),dp)

           scpr1(11)=scpr1(11)+apsi_f_temp(1)*real(bpsi_f11(1,jbj+i),dp)
           scpr2(11)=scpr2(11)+apsi_f_temp(2)*real(bpsi_f11(2,jbj+i),dp)
           scpr3(11)=scpr3(11)+apsi_f_temp(3)*real(bpsi_f11(3,jbj+i),dp)
           scpr4(11)=scpr4(11)+apsi_f_temp(4)*real(bpsi_f11(4,jbj+i),dp)
           scpr5(11)=scpr5(11)+apsi_f_temp(5)*real(bpsi_f11(5,jbj+i),dp)
           scpr6(11)=scpr6(11)+apsi_f_temp(6)*real(bpsi_f11(6,jbj+i),dp)
           scpr7(11)=scpr7(11)+apsi_f_temp(7)*real(bpsi_f11(7,jbj+i),dp)

           scpr1(12)=scpr1(12)+apsi_f_temp(1)*real(bpsi_f12(1,jbj+i),dp)
           scpr2(12)=scpr2(12)+apsi_f_temp(2)*real(bpsi_f12(2,jbj+i),dp)
           scpr3(12)=scpr3(12)+apsi_f_temp(3)*real(bpsi_f12(3,jbj+i),dp)
           scpr4(12)=scpr4(12)+apsi_f_temp(4)*real(bpsi_f12(4,jbj+i),dp)
           scpr5(12)=scpr5(12)+apsi_f_temp(5)*real(bpsi_f12(5,jbj+i),dp)
           scpr6(12)=scpr6(12)+apsi_f_temp(6)*real(bpsi_f12(6,jbj+i),dp)
           scpr7(12)=scpr7(12)+apsi_f_temp(7)*real(bpsi_f12(7,jbj+i),dp)

           scpr1(13)=scpr1(13)+apsi_f_temp(1)*real(bpsi_f13(1,jbj+i),dp)
           scpr2(13)=scpr2(13)+apsi_f_temp(2)*real(bpsi_f13(2,jbj+i),dp)
           scpr3(13)=scpr3(13)+apsi_f_temp(3)*real(bpsi_f13(3,jbj+i),dp)
           scpr4(13)=scpr4(13)+apsi_f_temp(4)*real(bpsi_f13(4,jbj+i),dp)
           scpr5(13)=scpr5(13)+apsi_f_temp(5)*real(bpsi_f13(5,jbj+i),dp)
           scpr6(13)=scpr6(13)+apsi_f_temp(6)*real(bpsi_f13(6,jbj+i),dp)
           scpr7(13)=scpr7(13)+apsi_f_temp(7)*real(bpsi_f13(7,jbj+i),dp)

           scpr1(14)=scpr1(14)+apsi_f_temp(1)*real(bpsi_f14(1,jbj+i),dp)
           scpr2(14)=scpr2(14)+apsi_f_temp(2)*real(bpsi_f14(2,jbj+i),dp)
           scpr3(14)=scpr3(14)+apsi_f_temp(3)*real(bpsi_f14(3,jbj+i),dp)
           scpr4(14)=scpr4(14)+apsi_f_temp(4)*real(bpsi_f14(4,jbj+i),dp)
           scpr5(14)=scpr5(14)+apsi_f_temp(5)*real(bpsi_f14(5,jbj+i),dp)
           scpr6(14)=scpr6(14)+apsi_f_temp(6)*real(bpsi_f14(6,jbj+i),dp)
           scpr7(14)=scpr7(14)+apsi_f_temp(7)*real(bpsi_f14(7,jbj+i),dp)

           scpr1(15)=scpr1(15)+apsi_f_temp(1)*real(bpsi_f15(1,jbj+i),dp)
           scpr2(15)=scpr2(15)+apsi_f_temp(2)*real(bpsi_f15(2,jbj+i),dp)
           scpr3(15)=scpr3(15)+apsi_f_temp(3)*real(bpsi_f15(3,jbj+i),dp)
           scpr4(15)=scpr4(15)+apsi_f_temp(4)*real(bpsi_f15(4,jbj+i),dp)
           scpr5(15)=scpr5(15)+apsi_f_temp(5)*real(bpsi_f15(5,jbj+i),dp)
           scpr6(15)=scpr6(15)+apsi_f_temp(6)*real(bpsi_f15(6,jbj+i),dp)
           scpr7(15)=scpr7(15)+apsi_f_temp(7)*real(bpsi_f15(7,jbj+i),dp)

           scpr1(16)=scpr1(16)+apsi_f_temp(1)*real(bpsi_f16(1,jbj+i),dp)
           scpr2(16)=scpr2(16)+apsi_f_temp(2)*real(bpsi_f16(2,jbj+i),dp)
           scpr3(16)=scpr3(16)+apsi_f_temp(3)*real(bpsi_f16(3,jbj+i),dp)
           scpr4(16)=scpr4(16)+apsi_f_temp(4)*real(bpsi_f16(4,jbj+i),dp)
           scpr5(16)=scpr5(16)+apsi_f_temp(5)*real(bpsi_f16(5,jbj+i),dp)
           scpr6(16)=scpr6(16)+apsi_f_temp(6)*real(bpsi_f16(6,jbj+i),dp)
           scpr7(16)=scpr7(16)+apsi_f_temp(7)*real(bpsi_f16(7,jbj+i),dp)

           scpr1(17)=scpr1(17)+apsi_f_temp(1)*real(bpsi_f17(1,jbj+i),dp)
           scpr2(17)=scpr2(17)+apsi_f_temp(2)*real(bpsi_f17(2,jbj+i),dp)
           scpr3(17)=scpr3(17)+apsi_f_temp(3)*real(bpsi_f17(3,jbj+i),dp)
           scpr4(17)=scpr4(17)+apsi_f_temp(4)*real(bpsi_f17(4,jbj+i),dp)
           scpr5(17)=scpr5(17)+apsi_f_temp(5)*real(bpsi_f17(5,jbj+i),dp)
           scpr6(17)=scpr6(17)+apsi_f_temp(6)*real(bpsi_f17(6,jbj+i),dp)
           scpr7(17)=scpr7(17)+apsi_f_temp(7)*real(bpsi_f17(7,jbj+i),dp)

           scpr1(18)=scpr1(18)+apsi_f_temp(1)*real(bpsi_f18(1,jbj+i),dp)
           scpr2(18)=scpr2(18)+apsi_f_temp(2)*real(bpsi_f18(2,jbj+i),dp)
           scpr3(18)=scpr3(18)+apsi_f_temp(3)*real(bpsi_f18(3,jbj+i),dp)
           scpr4(18)=scpr4(18)+apsi_f_temp(4)*real(bpsi_f18(4,jbj+i),dp)
           scpr5(18)=scpr5(18)+apsi_f_temp(5)*real(bpsi_f18(5,jbj+i),dp)
           scpr6(18)=scpr6(18)+apsi_f_temp(6)*real(bpsi_f18(6,jbj+i),dp)
           scpr7(18)=scpr7(18)+apsi_f_temp(7)*real(bpsi_f18(7,jbj+i),dp)

        enddo
 !       print *,'length',length,ibseg,scpr1,iaseg0,ja1,jb1
        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_f) exit nonconvex_loop_f !segment is finished  
        iaseg0=iaseg0+1
        jb0=max(keybg_f(1,ibseg),keyag_f_lin(iaseg0))
        if (keyag_f_lin(iaseg0)>jb1) exit nonconvex_loop_f !segment is not covered 
        jbj=jbj+max(jb0-keybg_f(1,ibseg),0)
     end do nonconvex_loop_f
     !disable loop if the end is reached
     if (iaseg0 == maseg_f .and. keybg_f(1,ibseg)> keyag_f_lin(maseg_f)) iaseg0=iaseg0+1

   enddo
!$omp end do !!!implicit barrier 

   scpr0=scpr0+scpr1+scpr2+scpr3+scpr4+scpr5+scpr6+scpr7

!$omp critical 
   scpr=scpr+scpr0
!$omp end critical
  
!$omp end parallel

!!!    call system_clock(ncount2,ncount_rate,ncount_max)
!!!    tel=dble(ncount2-ncount0)/dble(ncount_rate)
!!!    write(97,*) 'wpdot:',tel
!!!    close(97)

END SUBROUTINE wpdot_18

subroutine wpdot_19(  &
     mavctr_c,mavctr_f,maseg_c,maseg_f,keyav_c,keyav_f,keyag_c,keyag_f,apsi_c,apsi_f,  &
     mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keybv_c,keybv_f,keybg_c,keybg_f, &
     bpsi1,bpsi2,bpsi3,bpsi4,bpsi5,bpsi6,bpsi7,bpsi8,bpsi9,bpsi10,bpsi11,bpsi12,bpsi13,bpsi14, &
     bpsi15,bpsi16,bpsi17,bpsi18,bpsi19,bpsi_f1,bpsi_f2,bpsi_f3,bpsi_f4,bpsi_f5,bpsi_f6,bpsi_f7, &
     bpsi_f8,bpsi_f9,bpsi_f10,bpsi_f11,bpsi_f12,bpsi_f13,bpsi_f14,bpsi_f15,bpsi_f16,&
     bpsi_f17,bpsi_f18,bpsi_f19, &
     scpr,proj_count)

  use module_base
  implicit none
  integer, intent(in) :: mavctr_c,mavctr_f,maseg_c,maseg_f,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f
  integer, intent(in) :: proj_count
  integer, dimension(maseg_c), intent(in) :: keyav_c
  integer, dimension(mbseg_c), intent(in) :: keybv_c
  integer, dimension(maseg_f), intent(in) :: keyav_f
  integer, dimension(mbseg_f), intent(in) :: keybv_f
  integer, dimension(2,maseg_c), intent(in) :: keyag_c
  integer, dimension(2,mbseg_c), intent(in) :: keybg_c
  integer, dimension(2,maseg_f), intent(in) :: keyag_f
  integer, dimension(2,mbseg_f), intent(in) :: keybg_f
  real(wp), dimension(mavctr_c), intent(in) :: apsi_c
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi1,bpsi2,bpsi3,bpsi4,bpsi5,bpsi6,bpsi7
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi8,bpsi9,bpsi10,bpsi11,bpsi12,bpsi13,bpsi14
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi15,bpsi16,bpsi17,bpsi18,bpsi19
  real(wp), dimension(7,mavctr_f), intent(in) :: apsi_f
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f1,bpsi_f2,bpsi_f3,bpsi_f4,bpsi_f5,bpsi_f6,bpsi_f7
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f8,bpsi_f9,bpsi_f10,bpsi_f11,bpsi_f12,bpsi_f13,bpsi_f14
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f15,bpsi_f16,bpsi_f17,bpsi_f18,bpsi_f19
  real(dp), dimension(proj_count),intent(out) :: scpr
  !local variables
  integer :: ibseg,jaj,jb1,jb0,jbj,iaoff,length,i,ja0,ja1
  real(dp),dimension(proj_count) :: scpr1,scpr2,scpr3,scpr4,scpr5,scpr6,scpr7,scpr0
  integer :: iaseg0
  integer, dimension(maseg_c) :: keyag_c_lin !>linear version of second indices of keyag_c
  integer, dimension(maseg_f) :: keyag_f_lin !>linear version of second indices of keyag_f
  real(dp) :: apsi_temp
  real(dp), dimension(7) :: apsi_f_temp
!!!    integer :: ncount0,ncount2,ncount_rate,ncount_max
!!!    real(gp) :: tel

!!!  !dee
!!!    open(unit=97,file='time_wpdot',status='unknown',position='append')
!!!    call system_clock(ncount0,ncount_rate,ncount_max)

  keyag_c_lin = keyag_c(1,:)!speed up access in hunt subroutine by consecutive arrangement in memory
  keyag_f_lin = keyag_f(1,:)!speed up access in hunt subroutine by consecutive arrangement in memory

  scpr=0.0_dp

!$omp parallel default(private) shared(maseg_c,keyav_c,keyag_c,keyag_c_lin) &
!$omp shared(keybg_c,mbseg_c,keybv_c,mbseg_f,maseg_f,apsi_c,bpsi1,bpsi2) &
!$omp shared(bpsi3,bpsi4,bpsi5,bpsi6,bpsi7,bpsi8,bpsi9,bpsi10,bpsi11) &
!$omp shared(bpsi12,bpsi13,bpsi14,bpsi15,bpsi16,bpsi17,bpsi18,bpsi19) &
!$omp shared(bpsi_f1,bpsi_f2,bpsi_f3,bpsi_f4,bpsi_f5,bpsi_f6,bpsi_f7) &
!$omp shared(bpsi_f8,bpsi_f9,bpsi_f10,bpsi_f11,bpsi_f12,bpsi_f13,bpsi_f14) &
!$omp shared(bpsi_f15,bpsi_f16,bpsi_f17,bpsi_f18,bpsi_f19,keybv_f,keybg_f) &
!$omp shared(keyag_f,keyag_f_lin,keyav_f,apsi_f,scpr)

!!!!$omp shared (ncount0,ncount2,ncount_rate,ncount_max,tel)

  iaseg0=1 

  scpr0=0.0_dp
  scpr1=0.0_dp
  scpr2=0.0_dp
  scpr3=0.0_dp
  scpr4=0.0_dp
  scpr5=0.0_dp
  scpr6=0.0_dp
  scpr7=0.0_dp
!coarse part. Loop on the projectors segments
!$omp do schedule(static)
   do ibseg=1,mbseg_c
     jbj=keybv_c(ibseg)
     jb0=keybg_c(1,ibseg) !starting point of projector segment
     jb1=keybg_c(2,ibseg) !ending point of projector segment
!     print *,'huntenter',ibseg,jb0,jb1
 
     !find the starting point of the wavefunction segment
     !warning: hunt is assuming that the variable is always found
     !if it is not, iaseg0 is put to maseg + 1 so that the loop is disabled
     call hunt1(.true.,keyag_c_lin,maseg_c,keybg_c(1,ibseg),iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     !now pass through all the wavefunction segments until the end of the segment is 
     !still contained in projector segment
     nonconvex_loop_c: do while(iaseg0 <= maseg_c)
!     print *,'huntexit',iaseg0,maseg_c,keyag_c_lin(iaseg0),keyag_c(2,iaseg0)

        !length = jb1-jb0
        !iaoff = jb0-keyag_c_lin(iaseg0)!jb0-ja0

        ja0=keyag_c_lin(iaseg0)
        ja1=min(jb1,keyag_c(2,iaseg0)) 
        length = ja1-jb0
        iaoff = max(jb0-ja0,0) !no offset if we are already inside

        jaj=keyav_c(iaseg0)
        do i=0,length
           apsi_temp=real(apsi_c(jaj+iaoff+i),dp)

           scpr0(1)=scpr0(1)+apsi_temp *real(bpsi1(jbj+i),dp)
           scpr0(2)=scpr0(2)+apsi_temp *real(bpsi2(jbj+i),dp)
           scpr0(3)=scpr0(3)+apsi_temp *real(bpsi3(jbj+i),dp)
           scpr0(4)=scpr0(4)+apsi_temp *real(bpsi4(jbj+i),dp)
           scpr0(5)=scpr0(5)+apsi_temp *real(bpsi5(jbj+i),dp)
           scpr0(6)=scpr0(6)+apsi_temp *real(bpsi6(jbj+i),dp)
           scpr0(7)=scpr0(7)+apsi_temp *real(bpsi7(jbj+i),dp)
           scpr0(8)=scpr0(8)+apsi_temp *real(bpsi8(jbj+i),dp)
           scpr0(9)=scpr0(9)+apsi_temp *real(bpsi9(jbj+i),dp)
           scpr0(10)=scpr0(10)+apsi_temp *real(bpsi10(jbj+i),dp)
           scpr0(11)=scpr0(11)+apsi_temp *real(bpsi11(jbj+i),dp)
           scpr0(12)=scpr0(12)+apsi_temp *real(bpsi12(jbj+i),dp)
           scpr0(13)=scpr0(13)+apsi_temp *real(bpsi13(jbj+i),dp)
           scpr0(14)=scpr0(14)+apsi_temp *real(bpsi14(jbj+i),dp)
           scpr0(15)=scpr0(15)+apsi_temp *real(bpsi15(jbj+i),dp)
           scpr0(16)=scpr0(16)+apsi_temp *real(bpsi16(jbj+i),dp)
           scpr0(17)=scpr0(17)+apsi_temp *real(bpsi17(jbj+i),dp)
           scpr0(18)=scpr0(18)+apsi_temp *real(bpsi18(jbj+i),dp)
           scpr0(19)=scpr0(19)+apsi_temp *real(bpsi19(jbj+i),dp)
        enddo
 !       print *,'length',length,ibseg,scpr0,iaseg0,ja1,jb1

        !print *,'ibseg,mbseg_c,iaseg0,maseg_c',ibseg,mbseg_c,iaseg0,maseg_c
        !print '(a,6(i8),1pe25.17)','ja0,ja1t,ja1,jb0,jb1',&
        !     ibseg,ja0,keyag_c(2,iaseg0),ja1,jb0,jb1,scpr0
        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_c) exit nonconvex_loop_c !segment is finished
        iaseg0=iaseg0+1
        jb0=max(keybg_c(1,ibseg),keyag_c_lin(iaseg0))
        if (keyag_c_lin(iaseg0)>jb1) exit nonconvex_loop_c !segment is not covered
        jbj=jbj+max(jb0-keybg_c(1,ibseg),0)
     end do nonconvex_loop_c
     !disable loop if the end is reached
     if (iaseg0 == maseg_c .and. keybg_c(1,ibseg)> keyag_c_lin(maseg_c)) iaseg0=iaseg0+1


   enddo
!stop
!$omp end do nowait

! fine part

iaseg0=1

!$omp do schedule(static)
   do ibseg=1,mbseg_f
     jbj=keybv_f(ibseg)
     jb0=keybg_f(1,ibseg)
     jb1=keybg_f(2,ibseg)
!    print *,'huntenter',ibseg,jb0,jb1
     call hunt1(.true.,keyag_f_lin,maseg_f,keybg_f(1,ibseg),iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     nonconvex_loop_f: do while(iaseg0 <= maseg_f)
!!$     length = jb1-jb0
!!$     iaoff = jb0-keyag_f_lin(iaseg0)

        ja0=keyag_f_lin(iaseg0) !still doubts about copying in automatic array
        ja1=min(jb1,keyag_f(2,iaseg0)) 
        length = ja1-jb0
        iaoff = max(jb0-ja0,0) !no offset if we are already inside

        jaj=keyav_f(iaseg0)
        do i=0,length
           apsi_f_temp(1)=real(apsi_f(1,jaj+iaoff+i),dp)
           apsi_f_temp(2)=real(apsi_f(2,jaj+iaoff+i),dp)
           apsi_f_temp(3)=real(apsi_f(3,jaj+iaoff+i),dp)
           apsi_f_temp(4)=real(apsi_f(4,jaj+iaoff+i),dp)
           apsi_f_temp(5)=real(apsi_f(5,jaj+iaoff+i),dp)
           apsi_f_temp(6)=real(apsi_f(6,jaj+iaoff+i),dp)
           apsi_f_temp(7)=real(apsi_f(7,jaj+iaoff+i),dp)

           scpr1(1)=scpr1(1)+apsi_f_temp(1)*real(bpsi_f1(1,jbj+i),dp)
           scpr2(1)=scpr2(1)+apsi_f_temp(2)*real(bpsi_f1(2,jbj+i),dp)
           scpr3(1)=scpr3(1)+apsi_f_temp(3)*real(bpsi_f1(3,jbj+i),dp)
           scpr4(1)=scpr4(1)+apsi_f_temp(4)*real(bpsi_f1(4,jbj+i),dp)
           scpr5(1)=scpr5(1)+apsi_f_temp(5)*real(bpsi_f1(5,jbj+i),dp)
           scpr6(1)=scpr6(1)+apsi_f_temp(6)*real(bpsi_f1(6,jbj+i),dp)
           scpr7(1)=scpr7(1)+apsi_f_temp(7)*real(bpsi_f1(7,jbj+i),dp)

           scpr1(2)=scpr1(2)+apsi_f_temp(1)*real(bpsi_f2(1,jbj+i),dp)
           scpr2(2)=scpr2(2)+apsi_f_temp(2)*real(bpsi_f2(2,jbj+i),dp)
           scpr3(2)=scpr3(2)+apsi_f_temp(3)*real(bpsi_f2(3,jbj+i),dp)
           scpr4(2)=scpr4(2)+apsi_f_temp(4)*real(bpsi_f2(4,jbj+i),dp)
           scpr5(2)=scpr5(2)+apsi_f_temp(5)*real(bpsi_f2(5,jbj+i),dp)
           scpr6(2)=scpr6(2)+apsi_f_temp(6)*real(bpsi_f2(6,jbj+i),dp)
           scpr7(2)=scpr7(2)+apsi_f_temp(7)*real(bpsi_f2(7,jbj+i),dp)

           scpr1(3)=scpr1(3)+apsi_f_temp(1)*real(bpsi_f3(1,jbj+i),dp)
           scpr2(3)=scpr2(3)+apsi_f_temp(2)*real(bpsi_f3(2,jbj+i),dp)
           scpr3(3)=scpr3(3)+apsi_f_temp(3)*real(bpsi_f3(3,jbj+i),dp)
           scpr4(3)=scpr4(3)+apsi_f_temp(4)*real(bpsi_f3(4,jbj+i),dp)
           scpr5(3)=scpr5(3)+apsi_f_temp(5)*real(bpsi_f3(5,jbj+i),dp)
           scpr6(3)=scpr6(3)+apsi_f_temp(6)*real(bpsi_f3(6,jbj+i),dp)
           scpr7(3)=scpr7(3)+apsi_f_temp(7)*real(bpsi_f3(7,jbj+i),dp)

           scpr1(4)=scpr1(4)+apsi_f_temp(1)*real(bpsi_f4(1,jbj+i),dp)
           scpr2(4)=scpr2(4)+apsi_f_temp(2)*real(bpsi_f4(2,jbj+i),dp)
           scpr3(4)=scpr3(4)+apsi_f_temp(3)*real(bpsi_f4(3,jbj+i),dp)
           scpr4(4)=scpr4(4)+apsi_f_temp(4)*real(bpsi_f4(4,jbj+i),dp)
           scpr5(4)=scpr5(4)+apsi_f_temp(5)*real(bpsi_f4(5,jbj+i),dp)
           scpr6(4)=scpr6(4)+apsi_f_temp(6)*real(bpsi_f4(6,jbj+i),dp)
           scpr7(4)=scpr7(4)+apsi_f_temp(7)*real(bpsi_f4(7,jbj+i),dp)

           scpr1(5)=scpr1(5)+apsi_f_temp(1)*real(bpsi_f5(1,jbj+i),dp)
           scpr2(5)=scpr2(5)+apsi_f_temp(2)*real(bpsi_f5(2,jbj+i),dp)
           scpr3(5)=scpr3(5)+apsi_f_temp(3)*real(bpsi_f5(3,jbj+i),dp)
           scpr4(5)=scpr4(5)+apsi_f_temp(4)*real(bpsi_f5(4,jbj+i),dp)
           scpr5(5)=scpr5(5)+apsi_f_temp(5)*real(bpsi_f5(5,jbj+i),dp)
           scpr6(5)=scpr6(5)+apsi_f_temp(6)*real(bpsi_f5(6,jbj+i),dp)
           scpr7(5)=scpr7(5)+apsi_f_temp(7)*real(bpsi_f5(7,jbj+i),dp)

           scpr1(6)=scpr1(6)+apsi_f_temp(1)*real(bpsi_f6(1,jbj+i),dp)
           scpr2(6)=scpr2(6)+apsi_f_temp(2)*real(bpsi_f6(2,jbj+i),dp)
           scpr3(6)=scpr3(6)+apsi_f_temp(3)*real(bpsi_f6(3,jbj+i),dp)
           scpr4(6)=scpr4(6)+apsi_f_temp(4)*real(bpsi_f6(4,jbj+i),dp)
           scpr5(6)=scpr5(6)+apsi_f_temp(5)*real(bpsi_f6(5,jbj+i),dp)
           scpr6(6)=scpr6(6)+apsi_f_temp(6)*real(bpsi_f6(6,jbj+i),dp)
           scpr7(6)=scpr7(6)+apsi_f_temp(7)*real(bpsi_f6(7,jbj+i),dp)

           scpr1(7)=scpr1(7)+apsi_f_temp(1)*real(bpsi_f7(1,jbj+i),dp)
           scpr2(7)=scpr2(7)+apsi_f_temp(2)*real(bpsi_f7(2,jbj+i),dp)
           scpr3(7)=scpr3(7)+apsi_f_temp(3)*real(bpsi_f7(3,jbj+i),dp)
           scpr4(7)=scpr4(7)+apsi_f_temp(4)*real(bpsi_f7(4,jbj+i),dp)
           scpr5(7)=scpr5(7)+apsi_f_temp(5)*real(bpsi_f7(5,jbj+i),dp)
           scpr6(7)=scpr6(7)+apsi_f_temp(6)*real(bpsi_f7(6,jbj+i),dp)
           scpr7(7)=scpr7(7)+apsi_f_temp(7)*real(bpsi_f7(7,jbj+i),dp)

           scpr1(8)=scpr1(8)+apsi_f_temp(1)*real(bpsi_f8(1,jbj+i),dp)
           scpr2(8)=scpr2(8)+apsi_f_temp(2)*real(bpsi_f8(2,jbj+i),dp)
           scpr3(8)=scpr3(8)+apsi_f_temp(3)*real(bpsi_f8(3,jbj+i),dp)
           scpr4(8)=scpr4(8)+apsi_f_temp(4)*real(bpsi_f8(4,jbj+i),dp)
           scpr5(8)=scpr5(8)+apsi_f_temp(5)*real(bpsi_f8(5,jbj+i),dp)
           scpr6(8)=scpr6(8)+apsi_f_temp(6)*real(bpsi_f8(6,jbj+i),dp)
           scpr7(8)=scpr7(8)+apsi_f_temp(7)*real(bpsi_f8(7,jbj+i),dp)

           scpr1(9)=scpr1(9)+apsi_f_temp(1)*real(bpsi_f9(1,jbj+i),dp)
           scpr2(9)=scpr2(9)+apsi_f_temp(2)*real(bpsi_f9(2,jbj+i),dp)
           scpr3(9)=scpr3(9)+apsi_f_temp(3)*real(bpsi_f9(3,jbj+i),dp)
           scpr4(9)=scpr4(9)+apsi_f_temp(4)*real(bpsi_f9(4,jbj+i),dp)
           scpr5(9)=scpr5(9)+apsi_f_temp(5)*real(bpsi_f9(5,jbj+i),dp)
           scpr6(9)=scpr6(9)+apsi_f_temp(6)*real(bpsi_f9(6,jbj+i),dp)
           scpr7(9)=scpr7(9)+apsi_f_temp(7)*real(bpsi_f9(7,jbj+i),dp)

           scpr1(10)=scpr1(10)+apsi_f_temp(1)*real(bpsi_f10(1,jbj+i),dp)
           scpr2(10)=scpr2(10)+apsi_f_temp(2)*real(bpsi_f10(2,jbj+i),dp)
           scpr3(10)=scpr3(10)+apsi_f_temp(3)*real(bpsi_f10(3,jbj+i),dp)
           scpr4(10)=scpr4(10)+apsi_f_temp(4)*real(bpsi_f10(4,jbj+i),dp)
           scpr5(10)=scpr5(10)+apsi_f_temp(5)*real(bpsi_f10(5,jbj+i),dp)
           scpr6(10)=scpr6(10)+apsi_f_temp(6)*real(bpsi_f10(6,jbj+i),dp)
           scpr7(10)=scpr7(10)+apsi_f_temp(7)*real(bpsi_f10(7,jbj+i),dp)

           scpr1(11)=scpr1(11)+apsi_f_temp(1)*real(bpsi_f11(1,jbj+i),dp)
           scpr2(11)=scpr2(11)+apsi_f_temp(2)*real(bpsi_f11(2,jbj+i),dp)
           scpr3(11)=scpr3(11)+apsi_f_temp(3)*real(bpsi_f11(3,jbj+i),dp)
           scpr4(11)=scpr4(11)+apsi_f_temp(4)*real(bpsi_f11(4,jbj+i),dp)
           scpr5(11)=scpr5(11)+apsi_f_temp(5)*real(bpsi_f11(5,jbj+i),dp)
           scpr6(11)=scpr6(11)+apsi_f_temp(6)*real(bpsi_f11(6,jbj+i),dp)
           scpr7(11)=scpr7(11)+apsi_f_temp(7)*real(bpsi_f11(7,jbj+i),dp)

           scpr1(12)=scpr1(12)+apsi_f_temp(1)*real(bpsi_f12(1,jbj+i),dp)
           scpr2(12)=scpr2(12)+apsi_f_temp(2)*real(bpsi_f12(2,jbj+i),dp)
           scpr3(12)=scpr3(12)+apsi_f_temp(3)*real(bpsi_f12(3,jbj+i),dp)
           scpr4(12)=scpr4(12)+apsi_f_temp(4)*real(bpsi_f12(4,jbj+i),dp)
           scpr5(12)=scpr5(12)+apsi_f_temp(5)*real(bpsi_f12(5,jbj+i),dp)
           scpr6(12)=scpr6(12)+apsi_f_temp(6)*real(bpsi_f12(6,jbj+i),dp)
           scpr7(12)=scpr7(12)+apsi_f_temp(7)*real(bpsi_f12(7,jbj+i),dp)

           scpr1(13)=scpr1(13)+apsi_f_temp(1)*real(bpsi_f13(1,jbj+i),dp)
           scpr2(13)=scpr2(13)+apsi_f_temp(2)*real(bpsi_f13(2,jbj+i),dp)
           scpr3(13)=scpr3(13)+apsi_f_temp(3)*real(bpsi_f13(3,jbj+i),dp)
           scpr4(13)=scpr4(13)+apsi_f_temp(4)*real(bpsi_f13(4,jbj+i),dp)
           scpr5(13)=scpr5(13)+apsi_f_temp(5)*real(bpsi_f13(5,jbj+i),dp)
           scpr6(13)=scpr6(13)+apsi_f_temp(6)*real(bpsi_f13(6,jbj+i),dp)
           scpr7(13)=scpr7(13)+apsi_f_temp(7)*real(bpsi_f13(7,jbj+i),dp)

           scpr1(14)=scpr1(14)+apsi_f_temp(1)*real(bpsi_f14(1,jbj+i),dp)
           scpr2(14)=scpr2(14)+apsi_f_temp(2)*real(bpsi_f14(2,jbj+i),dp)
           scpr3(14)=scpr3(14)+apsi_f_temp(3)*real(bpsi_f14(3,jbj+i),dp)
           scpr4(14)=scpr4(14)+apsi_f_temp(4)*real(bpsi_f14(4,jbj+i),dp)
           scpr5(14)=scpr5(14)+apsi_f_temp(5)*real(bpsi_f14(5,jbj+i),dp)
           scpr6(14)=scpr6(14)+apsi_f_temp(6)*real(bpsi_f14(6,jbj+i),dp)
           scpr7(14)=scpr7(14)+apsi_f_temp(7)*real(bpsi_f14(7,jbj+i),dp)

           scpr1(15)=scpr1(15)+apsi_f_temp(1)*real(bpsi_f15(1,jbj+i),dp)
           scpr2(15)=scpr2(15)+apsi_f_temp(2)*real(bpsi_f15(2,jbj+i),dp)
           scpr3(15)=scpr3(15)+apsi_f_temp(3)*real(bpsi_f15(3,jbj+i),dp)
           scpr4(15)=scpr4(15)+apsi_f_temp(4)*real(bpsi_f15(4,jbj+i),dp)
           scpr5(15)=scpr5(15)+apsi_f_temp(5)*real(bpsi_f15(5,jbj+i),dp)
           scpr6(15)=scpr6(15)+apsi_f_temp(6)*real(bpsi_f15(6,jbj+i),dp)
           scpr7(15)=scpr7(15)+apsi_f_temp(7)*real(bpsi_f15(7,jbj+i),dp)

           scpr1(16)=scpr1(16)+apsi_f_temp(1)*real(bpsi_f16(1,jbj+i),dp)
           scpr2(16)=scpr2(16)+apsi_f_temp(2)*real(bpsi_f16(2,jbj+i),dp)
           scpr3(16)=scpr3(16)+apsi_f_temp(3)*real(bpsi_f16(3,jbj+i),dp)
           scpr4(16)=scpr4(16)+apsi_f_temp(4)*real(bpsi_f16(4,jbj+i),dp)
           scpr5(16)=scpr5(16)+apsi_f_temp(5)*real(bpsi_f16(5,jbj+i),dp)
           scpr6(16)=scpr6(16)+apsi_f_temp(6)*real(bpsi_f16(6,jbj+i),dp)
           scpr7(16)=scpr7(16)+apsi_f_temp(7)*real(bpsi_f16(7,jbj+i),dp)

           scpr1(17)=scpr1(17)+apsi_f_temp(1)*real(bpsi_f17(1,jbj+i),dp)
           scpr2(17)=scpr2(17)+apsi_f_temp(2)*real(bpsi_f17(2,jbj+i),dp)
           scpr3(17)=scpr3(17)+apsi_f_temp(3)*real(bpsi_f17(3,jbj+i),dp)
           scpr4(17)=scpr4(17)+apsi_f_temp(4)*real(bpsi_f17(4,jbj+i),dp)
           scpr5(17)=scpr5(17)+apsi_f_temp(5)*real(bpsi_f17(5,jbj+i),dp)
           scpr6(17)=scpr6(17)+apsi_f_temp(6)*real(bpsi_f17(6,jbj+i),dp)
           scpr7(17)=scpr7(17)+apsi_f_temp(7)*real(bpsi_f17(7,jbj+i),dp)

           scpr1(18)=scpr1(18)+apsi_f_temp(1)*real(bpsi_f18(1,jbj+i),dp)
           scpr2(18)=scpr2(18)+apsi_f_temp(2)*real(bpsi_f18(2,jbj+i),dp)
           scpr3(18)=scpr3(18)+apsi_f_temp(3)*real(bpsi_f18(3,jbj+i),dp)
           scpr4(18)=scpr4(18)+apsi_f_temp(4)*real(bpsi_f18(4,jbj+i),dp)
           scpr5(18)=scpr5(18)+apsi_f_temp(5)*real(bpsi_f18(5,jbj+i),dp)
           scpr6(18)=scpr6(18)+apsi_f_temp(6)*real(bpsi_f18(6,jbj+i),dp)
           scpr7(18)=scpr7(18)+apsi_f_temp(7)*real(bpsi_f18(7,jbj+i),dp)

           scpr1(19)=scpr1(19)+apsi_f_temp(1)*real(bpsi_f19(1,jbj+i),dp)
           scpr2(19)=scpr2(19)+apsi_f_temp(2)*real(bpsi_f19(2,jbj+i),dp)
           scpr3(19)=scpr3(19)+apsi_f_temp(3)*real(bpsi_f19(3,jbj+i),dp)
           scpr4(19)=scpr4(19)+apsi_f_temp(4)*real(bpsi_f19(4,jbj+i),dp)
           scpr5(19)=scpr5(19)+apsi_f_temp(5)*real(bpsi_f19(5,jbj+i),dp)
           scpr6(19)=scpr6(19)+apsi_f_temp(6)*real(bpsi_f19(6,jbj+i),dp)
           scpr7(19)=scpr7(19)+apsi_f_temp(7)*real(bpsi_f19(7,jbj+i),dp)

        enddo
 !       print *,'length',length,ibseg,scpr1,iaseg0,ja1,jb1
        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_f) exit nonconvex_loop_f !segment is finished  
        iaseg0=iaseg0+1
        jb0=max(keybg_f(1,ibseg),keyag_f_lin(iaseg0))
        if (keyag_f_lin(iaseg0)>jb1) exit nonconvex_loop_f !segment is not covered 
        jbj=jbj+max(jb0-keybg_f(1,ibseg),0)
     end do nonconvex_loop_f
     !disable loop if the end is reached
     if (iaseg0 == maseg_f .and. keybg_f(1,ibseg)> keyag_f_lin(maseg_f)) iaseg0=iaseg0+1

   enddo
!$omp end do !!!implicit barrier 

   scpr0=scpr0+scpr1+scpr2+scpr3+scpr4+scpr5+scpr6+scpr7

!$omp critical 
   scpr=scpr+scpr0
!$omp end critical
  
!$omp end parallel

!!!    call system_clock(ncount2,ncount_rate,ncount_max)
!!!    tel=dble(ncount2-ncount0)/dble(ncount_rate)
!!!    write(97,*) 'wpdot:',tel
!!!    close(97)

END SUBROUTINE wpdot_19

subroutine wpdot_20(  &
     mavctr_c,mavctr_f,maseg_c,maseg_f,keyav_c,keyav_f,keyag_c,keyag_f,apsi_c,apsi_f,  &
     mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keybv_c,keybv_f,keybg_c,keybg_f, &
     bpsi1,bpsi2,bpsi3,bpsi4,bpsi5,bpsi6,bpsi7,bpsi8,bpsi9,bpsi10,bpsi11,bpsi12,bpsi13,bpsi14, &
     bpsi15,bpsi16,bpsi17,bpsi18,bpsi19,bpsi20,bpsi_f1,bpsi_f2,bpsi_f3,bpsi_f4,bpsi_f5,bpsi_f6,bpsi_f7, &
     bpsi_f8,bpsi_f9,bpsi_f10,bpsi_f11,bpsi_f12,bpsi_f13,bpsi_f14,bpsi_f15,bpsi_f16,&
     bpsi_f17,bpsi_f18,bpsi_f19,bpsi_f20, &
     scpr,proj_count)

  use module_base
  implicit none
  integer, intent(in) :: mavctr_c,mavctr_f,maseg_c,maseg_f,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f
  integer, intent(in) :: proj_count
  integer, dimension(maseg_c), intent(in) :: keyav_c
  integer, dimension(mbseg_c), intent(in) :: keybv_c
  integer, dimension(maseg_f), intent(in) :: keyav_f
  integer, dimension(mbseg_f), intent(in) :: keybv_f
  integer, dimension(2,maseg_c), intent(in) :: keyag_c
  integer, dimension(2,mbseg_c), intent(in) :: keybg_c
  integer, dimension(2,maseg_f), intent(in) :: keyag_f
  integer, dimension(2,mbseg_f), intent(in) :: keybg_f
  real(wp), dimension(mavctr_c), intent(in) :: apsi_c
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi1,bpsi2,bpsi3,bpsi4,bpsi5,bpsi6,bpsi7
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi8,bpsi9,bpsi10,bpsi11,bpsi12,bpsi13,bpsi14
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi15,bpsi16,bpsi17,bpsi18,bpsi19,bpsi20
  real(wp), dimension(7,mavctr_f), intent(in) :: apsi_f
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f1,bpsi_f2,bpsi_f3,bpsi_f4,bpsi_f5,bpsi_f6,bpsi_f7
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f8,bpsi_f9,bpsi_f10,bpsi_f11,bpsi_f12,bpsi_f13,bpsi_f14
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f15,bpsi_f16,bpsi_f17,bpsi_f18,bpsi_f19,bpsi_f20
  real(dp), dimension(proj_count),intent(out) :: scpr
  !local variables
  integer :: ibseg,jaj,jb1,jb0,jbj,iaoff,length,i,ja0,ja1
  real(dp),dimension(proj_count) :: scpr1,scpr2,scpr3,scpr4,scpr5,scpr6,scpr7,scpr0
  integer :: iaseg0
  integer, dimension(maseg_c) :: keyag_c_lin !>linear version of second indices of keyag_c
  integer, dimension(maseg_f) :: keyag_f_lin !>linear version of second indices of keyag_f
  real(dp) :: apsi_temp
  real(dp), dimension(7) :: apsi_f_temp
!!!    integer :: ncount0,ncount2,ncount_rate,ncount_max
!!!    real(gp) :: tel

!!!  !dee
!!!    open(unit=97,file='time_wpdot',status='unknown',position='append')
!!!    call system_clock(ncount0,ncount_rate,ncount_max)

  keyag_c_lin = keyag_c(1,:)!speed up access in hunt subroutine by consecutive arrangement in memory
  keyag_f_lin = keyag_f(1,:)!speed up access in hunt subroutine by consecutive arrangement in memory

  scpr=0.0_dp

!$omp parallel default(private) shared(maseg_c,keyav_c,keyag_c,keyag_c_lin) &
!$omp shared(keybg_c,mbseg_c,keybv_c,mbseg_f,maseg_f,apsi_c,bpsi1,bpsi2) &
!$omp shared(bpsi3,bpsi4,bpsi5,bpsi6,bpsi7,bpsi8,bpsi9,bpsi10,bpsi11) &
!$omp shared(bpsi12,bpsi13,bpsi14,bpsi15,bpsi16,bpsi17,bpsi18,bpsi19) &
!$omp shared(bpsi20,bpsi_f1,bpsi_f2,bpsi_f3,bpsi_f4,bpsi_f5,bpsi_f6,bpsi_f7) &
!$omp shared(bpsi_f8,bpsi_f9,bpsi_f10,bpsi_f11,bpsi_f12,bpsi_f13,bpsi_f14) &
!$omp shared(bpsi_f15,bpsi_f16,bpsi_f17,bpsi_f18,bpsi_f19,bpsi_f20) &
!$omp shared(keybv_f,keybg_f,keyag_f,keyag_f_lin,keyav_f,apsi_f,scpr)

!!!!$omp shared (ncount0,ncount2,ncount_rate,ncount_max,tel)

  iaseg0=1 

  scpr0=0.0_dp
  scpr1=0.0_dp
  scpr2=0.0_dp
  scpr3=0.0_dp
  scpr4=0.0_dp
  scpr5=0.0_dp
  scpr6=0.0_dp
  scpr7=0.0_dp
!coarse part. Loop on the projectors segments
!$omp do schedule(static)
   do ibseg=1,mbseg_c
     jbj=keybv_c(ibseg)
     jb0=keybg_c(1,ibseg) !starting point of projector segment
     jb1=keybg_c(2,ibseg) !ending point of projector segment
!     print *,'huntenter',ibseg,jb0,jb1
 
     !find the starting point of the wavefunction segment
     !warning: hunt is assuming that the variable is always found
     !if it is not, iaseg0 is put to maseg + 1 so that the loop is disabled
     call hunt1(.true.,keyag_c_lin,maseg_c,keybg_c(1,ibseg),iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     !now pass through all the wavefunction segments until the end of the segment is 
     !still contained in projector segment
     nonconvex_loop_c: do while(iaseg0 <= maseg_c)
!     print *,'huntexit',iaseg0,maseg_c,keyag_c_lin(iaseg0),keyag_c(2,iaseg0)

        !length = jb1-jb0
        !iaoff = jb0-keyag_c_lin(iaseg0)!jb0-ja0

        ja0=keyag_c_lin(iaseg0)
        ja1=min(jb1,keyag_c(2,iaseg0)) 
        length = ja1-jb0
        iaoff = max(jb0-ja0,0) !no offset if we are already inside

        jaj=keyav_c(iaseg0)
        do i=0,length
           apsi_temp=real(apsi_c(jaj+iaoff+i),dp)

           scpr0(1)=scpr0(1)+apsi_temp *real(bpsi1(jbj+i),dp)
           scpr0(2)=scpr0(2)+apsi_temp *real(bpsi2(jbj+i),dp)
           scpr0(3)=scpr0(3)+apsi_temp *real(bpsi3(jbj+i),dp)
           scpr0(4)=scpr0(4)+apsi_temp *real(bpsi4(jbj+i),dp)
           scpr0(5)=scpr0(5)+apsi_temp *real(bpsi5(jbj+i),dp)
           scpr0(6)=scpr0(6)+apsi_temp *real(bpsi6(jbj+i),dp)
           scpr0(7)=scpr0(7)+apsi_temp *real(bpsi7(jbj+i),dp)
           scpr0(8)=scpr0(8)+apsi_temp *real(bpsi8(jbj+i),dp)
           scpr0(9)=scpr0(9)+apsi_temp *real(bpsi9(jbj+i),dp)
           scpr0(10)=scpr0(10)+apsi_temp *real(bpsi10(jbj+i),dp)
           scpr0(11)=scpr0(11)+apsi_temp *real(bpsi11(jbj+i),dp)
           scpr0(12)=scpr0(12)+apsi_temp *real(bpsi12(jbj+i),dp)
           scpr0(13)=scpr0(13)+apsi_temp *real(bpsi13(jbj+i),dp)
           scpr0(14)=scpr0(14)+apsi_temp *real(bpsi14(jbj+i),dp)
           scpr0(15)=scpr0(15)+apsi_temp *real(bpsi15(jbj+i),dp)
           scpr0(16)=scpr0(16)+apsi_temp *real(bpsi16(jbj+i),dp)
           scpr0(17)=scpr0(17)+apsi_temp *real(bpsi17(jbj+i),dp)
           scpr0(18)=scpr0(18)+apsi_temp *real(bpsi18(jbj+i),dp)
           scpr0(19)=scpr0(19)+apsi_temp *real(bpsi19(jbj+i),dp)
           scpr0(20)=scpr0(20)+apsi_temp *real(bpsi20(jbj+i),dp)
        enddo
 !       print *,'length',length,ibseg,scpr0,iaseg0,ja1,jb1

        !print *,'ibseg,mbseg_c,iaseg0,maseg_c',ibseg,mbseg_c,iaseg0,maseg_c
        !print '(a,6(i8),1pe25.17)','ja0,ja1t,ja1,jb0,jb1',&
        !     ibseg,ja0,keyag_c(2,iaseg0),ja1,jb0,jb1,scpr0
        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_c) exit nonconvex_loop_c !segment is finished
        iaseg0=iaseg0+1
        jb0=max(keybg_c(1,ibseg),keyag_c_lin(iaseg0))
        if (keyag_c_lin(iaseg0)>jb1) exit nonconvex_loop_c !segment is not covered
        jbj=jbj+max(jb0-keybg_c(1,ibseg),0)
     end do nonconvex_loop_c
     !disable loop if the end is reached
     if (iaseg0 == maseg_c .and. keybg_c(1,ibseg)> keyag_c_lin(maseg_c)) iaseg0=iaseg0+1


   enddo
!stop
!$omp end do nowait

! fine part

iaseg0=1

!$omp do schedule(static)
   do ibseg=1,mbseg_f
     jbj=keybv_f(ibseg)
     jb0=keybg_f(1,ibseg)
     jb1=keybg_f(2,ibseg)
!    print *,'huntenter',ibseg,jb0,jb1
     call hunt1(.true.,keyag_f_lin,maseg_f,keybg_f(1,ibseg),iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     nonconvex_loop_f: do while(iaseg0 <= maseg_f)
!!$     length = jb1-jb0
!!$     iaoff = jb0-keyag_f_lin(iaseg0)

        ja0=keyag_f_lin(iaseg0) !still doubts about copying in automatic array
        ja1=min(jb1,keyag_f(2,iaseg0)) 
        length = ja1-jb0
        iaoff = max(jb0-ja0,0) !no offset if we are already inside

        jaj=keyav_f(iaseg0)
        do i=0,length
           apsi_f_temp(1)=real(apsi_f(1,jaj+iaoff+i),dp)
           apsi_f_temp(2)=real(apsi_f(2,jaj+iaoff+i),dp)
           apsi_f_temp(3)=real(apsi_f(3,jaj+iaoff+i),dp)
           apsi_f_temp(4)=real(apsi_f(4,jaj+iaoff+i),dp)
           apsi_f_temp(5)=real(apsi_f(5,jaj+iaoff+i),dp)
           apsi_f_temp(6)=real(apsi_f(6,jaj+iaoff+i),dp)
           apsi_f_temp(7)=real(apsi_f(7,jaj+iaoff+i),dp)

           scpr1(1)=scpr1(1)+apsi_f_temp(1)*real(bpsi_f1(1,jbj+i),dp)
           scpr2(1)=scpr2(1)+apsi_f_temp(2)*real(bpsi_f1(2,jbj+i),dp)
           scpr3(1)=scpr3(1)+apsi_f_temp(3)*real(bpsi_f1(3,jbj+i),dp)
           scpr4(1)=scpr4(1)+apsi_f_temp(4)*real(bpsi_f1(4,jbj+i),dp)
           scpr5(1)=scpr5(1)+apsi_f_temp(5)*real(bpsi_f1(5,jbj+i),dp)
           scpr6(1)=scpr6(1)+apsi_f_temp(6)*real(bpsi_f1(6,jbj+i),dp)
           scpr7(1)=scpr7(1)+apsi_f_temp(7)*real(bpsi_f1(7,jbj+i),dp)

           scpr1(2)=scpr1(2)+apsi_f_temp(1)*real(bpsi_f2(1,jbj+i),dp)
           scpr2(2)=scpr2(2)+apsi_f_temp(2)*real(bpsi_f2(2,jbj+i),dp)
           scpr3(2)=scpr3(2)+apsi_f_temp(3)*real(bpsi_f2(3,jbj+i),dp)
           scpr4(2)=scpr4(2)+apsi_f_temp(4)*real(bpsi_f2(4,jbj+i),dp)
           scpr5(2)=scpr5(2)+apsi_f_temp(5)*real(bpsi_f2(5,jbj+i),dp)
           scpr6(2)=scpr6(2)+apsi_f_temp(6)*real(bpsi_f2(6,jbj+i),dp)
           scpr7(2)=scpr7(2)+apsi_f_temp(7)*real(bpsi_f2(7,jbj+i),dp)

           scpr1(3)=scpr1(3)+apsi_f_temp(1)*real(bpsi_f3(1,jbj+i),dp)
           scpr2(3)=scpr2(3)+apsi_f_temp(2)*real(bpsi_f3(2,jbj+i),dp)
           scpr3(3)=scpr3(3)+apsi_f_temp(3)*real(bpsi_f3(3,jbj+i),dp)
           scpr4(3)=scpr4(3)+apsi_f_temp(4)*real(bpsi_f3(4,jbj+i),dp)
           scpr5(3)=scpr5(3)+apsi_f_temp(5)*real(bpsi_f3(5,jbj+i),dp)
           scpr6(3)=scpr6(3)+apsi_f_temp(6)*real(bpsi_f3(6,jbj+i),dp)
           scpr7(3)=scpr7(3)+apsi_f_temp(7)*real(bpsi_f3(7,jbj+i),dp)

           scpr1(4)=scpr1(4)+apsi_f_temp(1)*real(bpsi_f4(1,jbj+i),dp)
           scpr2(4)=scpr2(4)+apsi_f_temp(2)*real(bpsi_f4(2,jbj+i),dp)
           scpr3(4)=scpr3(4)+apsi_f_temp(3)*real(bpsi_f4(3,jbj+i),dp)
           scpr4(4)=scpr4(4)+apsi_f_temp(4)*real(bpsi_f4(4,jbj+i),dp)
           scpr5(4)=scpr5(4)+apsi_f_temp(5)*real(bpsi_f4(5,jbj+i),dp)
           scpr6(4)=scpr6(4)+apsi_f_temp(6)*real(bpsi_f4(6,jbj+i),dp)
           scpr7(4)=scpr7(4)+apsi_f_temp(7)*real(bpsi_f4(7,jbj+i),dp)

           scpr1(5)=scpr1(5)+apsi_f_temp(1)*real(bpsi_f5(1,jbj+i),dp)
           scpr2(5)=scpr2(5)+apsi_f_temp(2)*real(bpsi_f5(2,jbj+i),dp)
           scpr3(5)=scpr3(5)+apsi_f_temp(3)*real(bpsi_f5(3,jbj+i),dp)
           scpr4(5)=scpr4(5)+apsi_f_temp(4)*real(bpsi_f5(4,jbj+i),dp)
           scpr5(5)=scpr5(5)+apsi_f_temp(5)*real(bpsi_f5(5,jbj+i),dp)
           scpr6(5)=scpr6(5)+apsi_f_temp(6)*real(bpsi_f5(6,jbj+i),dp)
           scpr7(5)=scpr7(5)+apsi_f_temp(7)*real(bpsi_f5(7,jbj+i),dp)

           scpr1(6)=scpr1(6)+apsi_f_temp(1)*real(bpsi_f6(1,jbj+i),dp)
           scpr2(6)=scpr2(6)+apsi_f_temp(2)*real(bpsi_f6(2,jbj+i),dp)
           scpr3(6)=scpr3(6)+apsi_f_temp(3)*real(bpsi_f6(3,jbj+i),dp)
           scpr4(6)=scpr4(6)+apsi_f_temp(4)*real(bpsi_f6(4,jbj+i),dp)
           scpr5(6)=scpr5(6)+apsi_f_temp(5)*real(bpsi_f6(5,jbj+i),dp)
           scpr6(6)=scpr6(6)+apsi_f_temp(6)*real(bpsi_f6(6,jbj+i),dp)
           scpr7(6)=scpr7(6)+apsi_f_temp(7)*real(bpsi_f6(7,jbj+i),dp)

           scpr1(7)=scpr1(7)+apsi_f_temp(1)*real(bpsi_f7(1,jbj+i),dp)
           scpr2(7)=scpr2(7)+apsi_f_temp(2)*real(bpsi_f7(2,jbj+i),dp)
           scpr3(7)=scpr3(7)+apsi_f_temp(3)*real(bpsi_f7(3,jbj+i),dp)
           scpr4(7)=scpr4(7)+apsi_f_temp(4)*real(bpsi_f7(4,jbj+i),dp)
           scpr5(7)=scpr5(7)+apsi_f_temp(5)*real(bpsi_f7(5,jbj+i),dp)
           scpr6(7)=scpr6(7)+apsi_f_temp(6)*real(bpsi_f7(6,jbj+i),dp)
           scpr7(7)=scpr7(7)+apsi_f_temp(7)*real(bpsi_f7(7,jbj+i),dp)

           scpr1(8)=scpr1(8)+apsi_f_temp(1)*real(bpsi_f8(1,jbj+i),dp)
           scpr2(8)=scpr2(8)+apsi_f_temp(2)*real(bpsi_f8(2,jbj+i),dp)
           scpr3(8)=scpr3(8)+apsi_f_temp(3)*real(bpsi_f8(3,jbj+i),dp)
           scpr4(8)=scpr4(8)+apsi_f_temp(4)*real(bpsi_f8(4,jbj+i),dp)
           scpr5(8)=scpr5(8)+apsi_f_temp(5)*real(bpsi_f8(5,jbj+i),dp)
           scpr6(8)=scpr6(8)+apsi_f_temp(6)*real(bpsi_f8(6,jbj+i),dp)
           scpr7(8)=scpr7(8)+apsi_f_temp(7)*real(bpsi_f8(7,jbj+i),dp)

           scpr1(9)=scpr1(9)+apsi_f_temp(1)*real(bpsi_f9(1,jbj+i),dp)
           scpr2(9)=scpr2(9)+apsi_f_temp(2)*real(bpsi_f9(2,jbj+i),dp)
           scpr3(9)=scpr3(9)+apsi_f_temp(3)*real(bpsi_f9(3,jbj+i),dp)
           scpr4(9)=scpr4(9)+apsi_f_temp(4)*real(bpsi_f9(4,jbj+i),dp)
           scpr5(9)=scpr5(9)+apsi_f_temp(5)*real(bpsi_f9(5,jbj+i),dp)
           scpr6(9)=scpr6(9)+apsi_f_temp(6)*real(bpsi_f9(6,jbj+i),dp)
           scpr7(9)=scpr7(9)+apsi_f_temp(7)*real(bpsi_f9(7,jbj+i),dp)

           scpr1(10)=scpr1(10)+apsi_f_temp(1)*real(bpsi_f10(1,jbj+i),dp)
           scpr2(10)=scpr2(10)+apsi_f_temp(2)*real(bpsi_f10(2,jbj+i),dp)
           scpr3(10)=scpr3(10)+apsi_f_temp(3)*real(bpsi_f10(3,jbj+i),dp)
           scpr4(10)=scpr4(10)+apsi_f_temp(4)*real(bpsi_f10(4,jbj+i),dp)
           scpr5(10)=scpr5(10)+apsi_f_temp(5)*real(bpsi_f10(5,jbj+i),dp)
           scpr6(10)=scpr6(10)+apsi_f_temp(6)*real(bpsi_f10(6,jbj+i),dp)
           scpr7(10)=scpr7(10)+apsi_f_temp(7)*real(bpsi_f10(7,jbj+i),dp)

           scpr1(11)=scpr1(11)+apsi_f_temp(1)*real(bpsi_f11(1,jbj+i),dp)
           scpr2(11)=scpr2(11)+apsi_f_temp(2)*real(bpsi_f11(2,jbj+i),dp)
           scpr3(11)=scpr3(11)+apsi_f_temp(3)*real(bpsi_f11(3,jbj+i),dp)
           scpr4(11)=scpr4(11)+apsi_f_temp(4)*real(bpsi_f11(4,jbj+i),dp)
           scpr5(11)=scpr5(11)+apsi_f_temp(5)*real(bpsi_f11(5,jbj+i),dp)
           scpr6(11)=scpr6(11)+apsi_f_temp(6)*real(bpsi_f11(6,jbj+i),dp)
           scpr7(11)=scpr7(11)+apsi_f_temp(7)*real(bpsi_f11(7,jbj+i),dp)

           scpr1(12)=scpr1(12)+apsi_f_temp(1)*real(bpsi_f12(1,jbj+i),dp)
           scpr2(12)=scpr2(12)+apsi_f_temp(2)*real(bpsi_f12(2,jbj+i),dp)
           scpr3(12)=scpr3(12)+apsi_f_temp(3)*real(bpsi_f12(3,jbj+i),dp)
           scpr4(12)=scpr4(12)+apsi_f_temp(4)*real(bpsi_f12(4,jbj+i),dp)
           scpr5(12)=scpr5(12)+apsi_f_temp(5)*real(bpsi_f12(5,jbj+i),dp)
           scpr6(12)=scpr6(12)+apsi_f_temp(6)*real(bpsi_f12(6,jbj+i),dp)
           scpr7(12)=scpr7(12)+apsi_f_temp(7)*real(bpsi_f12(7,jbj+i),dp)

           scpr1(13)=scpr1(13)+apsi_f_temp(1)*real(bpsi_f13(1,jbj+i),dp)
           scpr2(13)=scpr2(13)+apsi_f_temp(2)*real(bpsi_f13(2,jbj+i),dp)
           scpr3(13)=scpr3(13)+apsi_f_temp(3)*real(bpsi_f13(3,jbj+i),dp)
           scpr4(13)=scpr4(13)+apsi_f_temp(4)*real(bpsi_f13(4,jbj+i),dp)
           scpr5(13)=scpr5(13)+apsi_f_temp(5)*real(bpsi_f13(5,jbj+i),dp)
           scpr6(13)=scpr6(13)+apsi_f_temp(6)*real(bpsi_f13(6,jbj+i),dp)
           scpr7(13)=scpr7(13)+apsi_f_temp(7)*real(bpsi_f13(7,jbj+i),dp)

           scpr1(14)=scpr1(14)+apsi_f_temp(1)*real(bpsi_f14(1,jbj+i),dp)
           scpr2(14)=scpr2(14)+apsi_f_temp(2)*real(bpsi_f14(2,jbj+i),dp)
           scpr3(14)=scpr3(14)+apsi_f_temp(3)*real(bpsi_f14(3,jbj+i),dp)
           scpr4(14)=scpr4(14)+apsi_f_temp(4)*real(bpsi_f14(4,jbj+i),dp)
           scpr5(14)=scpr5(14)+apsi_f_temp(5)*real(bpsi_f14(5,jbj+i),dp)
           scpr6(14)=scpr6(14)+apsi_f_temp(6)*real(bpsi_f14(6,jbj+i),dp)
           scpr7(14)=scpr7(14)+apsi_f_temp(7)*real(bpsi_f14(7,jbj+i),dp)

           scpr1(15)=scpr1(15)+apsi_f_temp(1)*real(bpsi_f15(1,jbj+i),dp)
           scpr2(15)=scpr2(15)+apsi_f_temp(2)*real(bpsi_f15(2,jbj+i),dp)
           scpr3(15)=scpr3(15)+apsi_f_temp(3)*real(bpsi_f15(3,jbj+i),dp)
           scpr4(15)=scpr4(15)+apsi_f_temp(4)*real(bpsi_f15(4,jbj+i),dp)
           scpr5(15)=scpr5(15)+apsi_f_temp(5)*real(bpsi_f15(5,jbj+i),dp)
           scpr6(15)=scpr6(15)+apsi_f_temp(6)*real(bpsi_f15(6,jbj+i),dp)
           scpr7(15)=scpr7(15)+apsi_f_temp(7)*real(bpsi_f15(7,jbj+i),dp)

           scpr1(16)=scpr1(16)+apsi_f_temp(1)*real(bpsi_f16(1,jbj+i),dp)
           scpr2(16)=scpr2(16)+apsi_f_temp(2)*real(bpsi_f16(2,jbj+i),dp)
           scpr3(16)=scpr3(16)+apsi_f_temp(3)*real(bpsi_f16(3,jbj+i),dp)
           scpr4(16)=scpr4(16)+apsi_f_temp(4)*real(bpsi_f16(4,jbj+i),dp)
           scpr5(16)=scpr5(16)+apsi_f_temp(5)*real(bpsi_f16(5,jbj+i),dp)
           scpr6(16)=scpr6(16)+apsi_f_temp(6)*real(bpsi_f16(6,jbj+i),dp)
           scpr7(16)=scpr7(16)+apsi_f_temp(7)*real(bpsi_f16(7,jbj+i),dp)

           scpr1(17)=scpr1(17)+apsi_f_temp(1)*real(bpsi_f17(1,jbj+i),dp)
           scpr2(17)=scpr2(17)+apsi_f_temp(2)*real(bpsi_f17(2,jbj+i),dp)
           scpr3(17)=scpr3(17)+apsi_f_temp(3)*real(bpsi_f17(3,jbj+i),dp)
           scpr4(17)=scpr4(17)+apsi_f_temp(4)*real(bpsi_f17(4,jbj+i),dp)
           scpr5(17)=scpr5(17)+apsi_f_temp(5)*real(bpsi_f17(5,jbj+i),dp)
           scpr6(17)=scpr6(17)+apsi_f_temp(6)*real(bpsi_f17(6,jbj+i),dp)
           scpr7(17)=scpr7(17)+apsi_f_temp(7)*real(bpsi_f17(7,jbj+i),dp)

           scpr1(18)=scpr1(18)+apsi_f_temp(1)*real(bpsi_f18(1,jbj+i),dp)
           scpr2(18)=scpr2(18)+apsi_f_temp(2)*real(bpsi_f18(2,jbj+i),dp)
           scpr3(18)=scpr3(18)+apsi_f_temp(3)*real(bpsi_f18(3,jbj+i),dp)
           scpr4(18)=scpr4(18)+apsi_f_temp(4)*real(bpsi_f18(4,jbj+i),dp)
           scpr5(18)=scpr5(18)+apsi_f_temp(5)*real(bpsi_f18(5,jbj+i),dp)
           scpr6(18)=scpr6(18)+apsi_f_temp(6)*real(bpsi_f18(6,jbj+i),dp)
           scpr7(18)=scpr7(18)+apsi_f_temp(7)*real(bpsi_f18(7,jbj+i),dp)

           scpr1(19)=scpr1(19)+apsi_f_temp(1)*real(bpsi_f19(1,jbj+i),dp)
           scpr2(19)=scpr2(19)+apsi_f_temp(2)*real(bpsi_f19(2,jbj+i),dp)
           scpr3(19)=scpr3(19)+apsi_f_temp(3)*real(bpsi_f19(3,jbj+i),dp)
           scpr4(19)=scpr4(19)+apsi_f_temp(4)*real(bpsi_f19(4,jbj+i),dp)
           scpr5(19)=scpr5(19)+apsi_f_temp(5)*real(bpsi_f19(5,jbj+i),dp)
           scpr6(19)=scpr6(19)+apsi_f_temp(6)*real(bpsi_f19(6,jbj+i),dp)
           scpr7(19)=scpr7(19)+apsi_f_temp(7)*real(bpsi_f19(7,jbj+i),dp)

           scpr1(20)=scpr1(20)+apsi_f_temp(1)*real(bpsi_f20(1,jbj+i),dp)
           scpr2(20)=scpr2(20)+apsi_f_temp(2)*real(bpsi_f20(2,jbj+i),dp)
           scpr3(20)=scpr3(20)+apsi_f_temp(3)*real(bpsi_f20(3,jbj+i),dp)
           scpr4(20)=scpr4(20)+apsi_f_temp(4)*real(bpsi_f20(4,jbj+i),dp)
           scpr5(20)=scpr5(20)+apsi_f_temp(5)*real(bpsi_f20(5,jbj+i),dp)
           scpr6(20)=scpr6(20)+apsi_f_temp(6)*real(bpsi_f20(6,jbj+i),dp)
           scpr7(20)=scpr7(20)+apsi_f_temp(7)*real(bpsi_f20(7,jbj+i),dp)

        enddo
 !       print *,'length',length,ibseg,scpr1,iaseg0,ja1,jb1
        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_f) exit nonconvex_loop_f !segment is finished  
        iaseg0=iaseg0+1
        jb0=max(keybg_f(1,ibseg),keyag_f_lin(iaseg0))
        if (keyag_f_lin(iaseg0)>jb1) exit nonconvex_loop_f !segment is not covered 
        jbj=jbj+max(jb0-keybg_f(1,ibseg),0)
     end do nonconvex_loop_f
     !disable loop if the end is reached
     if (iaseg0 == maseg_f .and. keybg_f(1,ibseg)> keyag_f_lin(maseg_f)) iaseg0=iaseg0+1

   enddo
!$omp end do !!!implicit barrier 

   scpr0=scpr0+scpr1+scpr2+scpr3+scpr4+scpr5+scpr6+scpr7

!$omp critical 
   scpr=scpr+scpr0
!$omp end critical
  
!$omp end parallel

!!!    call system_clock(ncount2,ncount_rate,ncount_max)
!!!    tel=dble(ncount2-ncount0)/dble(ncount_rate)
!!!    write(97,*) 'wpdot:',tel
!!!    close(97)

END SUBROUTINE wpdot_20

subroutine wpdot_22(  &
     mavctr_c,mavctr_f,maseg_c,maseg_f,keyav_c,keyav_f,keyag_c,keyag_f,apsi_c,apsi_f,  &
     mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keybv_c,keybv_f,keybg_c,keybg_f, &
     bpsi1,bpsi2,bpsi3,bpsi4,bpsi5,bpsi6,bpsi7,bpsi8,bpsi9,bpsi10,bpsi11,bpsi12,bpsi13,bpsi14, &
     bpsi15,bpsi16,bpsi17,bpsi18,bpsi19,bpsi20,bpsi21,bpsi22,bpsi_f1,bpsi_f2,bpsi_f3,bpsi_f4, &
     bpsi_f5,bpsi_f6,bpsi_f7,bpsi_f8,bpsi_f9,bpsi_f10,bpsi_f11,bpsi_f12,bpsi_f13,bpsi_f14, &
     bpsi_f15,bpsi_f16,bpsi_f17,bpsi_f18,bpsi_f19,bpsi_f20,bpsi_f21,bpsi_f22, &
     scpr,proj_count)

  use module_base
  implicit none
  integer, intent(in) :: mavctr_c,mavctr_f,maseg_c,maseg_f,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f
  integer, intent(in) :: proj_count
  integer, dimension(maseg_c), intent(in) :: keyav_c
  integer, dimension(mbseg_c), intent(in) :: keybv_c
  integer, dimension(maseg_f), intent(in) :: keyav_f
  integer, dimension(mbseg_f), intent(in) :: keybv_f
  integer, dimension(2,maseg_c), intent(in) :: keyag_c
  integer, dimension(2,mbseg_c), intent(in) :: keybg_c
  integer, dimension(2,maseg_f), intent(in) :: keyag_f
  integer, dimension(2,mbseg_f), intent(in) :: keybg_f
  real(wp), dimension(mavctr_c), intent(in) :: apsi_c
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi1,bpsi2,bpsi3,bpsi4,bpsi5,bpsi6,bpsi7
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi8,bpsi9,bpsi10,bpsi11,bpsi12,bpsi13,bpsi14
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi15,bpsi16,bpsi17,bpsi18,bpsi19,bpsi20,bpsi21,bpsi22
  real(wp), dimension(7,mavctr_f), intent(in) :: apsi_f
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f1,bpsi_f2,bpsi_f3,bpsi_f4,bpsi_f5,bpsi_f6,bpsi_f7
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f8,bpsi_f9,bpsi_f10,bpsi_f11,bpsi_f12,bpsi_f13,bpsi_f14
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f15,bpsi_f16,bpsi_f17,bpsi_f18,bpsi_f19,bpsi_f20
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f21,bpsi_f22
  real(dp), dimension(proj_count),intent(out) :: scpr
  !local variables
  integer :: ibseg,jaj,jb1,jb0,jbj,iaoff,length,i,ja0,ja1
  real(dp),dimension(proj_count) :: scpr1,scpr2,scpr3,scpr4,scpr5,scpr6,scpr7,scpr0
  integer :: iaseg0
  integer, dimension(maseg_c) :: keyag_c_lin !>linear version of second indices of keyag_c
  integer, dimension(maseg_f) :: keyag_f_lin !>linear version of second indices of keyag_f
  real(dp) :: apsi_temp
  real(dp), dimension(7) :: apsi_f_temp
!!!    integer :: ncount0,ncount2,ncount_rate,ncount_max
!!!    real(gp) :: tel

!!!  !dee
!!!    open(unit=97,file='time_wpdot',status='unknown',position='append')
!!!    call system_clock(ncount0,ncount_rate,ncount_max)

  keyag_c_lin = keyag_c(1,:)!speed up access in hunt subroutine by consecutive arrangement in memory
  keyag_f_lin = keyag_f(1,:)!speed up access in hunt subroutine by consecutive arrangement in memory

  scpr=0.0_dp

!$omp parallel default(private) shared(maseg_c,keyav_c,keyag_c,keyag_c_lin) &
!$omp shared(keybg_c,mbseg_c,keybv_c,mbseg_f,maseg_f,apsi_c,bpsi1,bpsi2) &
!$omp shared(bpsi3,bpsi4,bpsi5,bpsi6,bpsi7,bpsi8,bpsi9,bpsi10,bpsi11) &
!$omp shared(bpsi12,bpsi13,bpsi14,bpsi15,bpsi16,bpsi17,bpsi18,bpsi19) &
!$omp shared(bpsi20,bpsi21,bpsi22,bpsi_f1,bpsi_f2,bpsi_f3,bpsi_f4,bpsi_f5) &
!$omp shared(bpsi_f6,bpsi_f7,bpsi_f8,bpsi_f9,bpsi_f10,bpsi_f11,bpsi_f12) &
!$omp shared(bpsi_f13,bpsi_f14,bpsi_f15,bpsi_f16,bpsi_f17,bpsi_f18) &
!$omp shared(bpsi_f19,bpsi_f20,bpsi_f21,bpsi_f22) &
!$omp shared(keybv_f,keybg_f,keyag_f,keyag_f_lin,keyav_f,apsi_f,scpr)

!!!!$omp shared (ncount0,ncount2,ncount_rate,ncount_max,tel)

  iaseg0=1 

  scpr0=0.0_dp
  scpr1=0.0_dp
  scpr2=0.0_dp
  scpr3=0.0_dp
  scpr4=0.0_dp
  scpr5=0.0_dp
  scpr6=0.0_dp
  scpr7=0.0_dp
!coarse part. Loop on the projectors segments
!$omp do schedule(static)
   do ibseg=1,mbseg_c
     jbj=keybv_c(ibseg)
     jb0=keybg_c(1,ibseg) !starting point of projector segment
     jb1=keybg_c(2,ibseg) !ending point of projector segment
!     print *,'huntenter',ibseg,jb0,jb1
 
     !find the starting point of the wavefunction segment
     !warning: hunt is assuming that the variable is always found
     !if it is not, iaseg0 is put to maseg + 1 so that the loop is disabled
     call hunt1(.true.,keyag_c_lin,maseg_c,keybg_c(1,ibseg),iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     !now pass through all the wavefunction segments until the end of the segment is 
     !still contained in projector segment
     nonconvex_loop_c: do while(iaseg0 <= maseg_c)
!     print *,'huntexit',iaseg0,maseg_c,keyag_c_lin(iaseg0),keyag_c(2,iaseg0)

        !length = jb1-jb0
        !iaoff = jb0-keyag_c_lin(iaseg0)!jb0-ja0

        ja0=keyag_c_lin(iaseg0)
        ja1=min(jb1,keyag_c(2,iaseg0)) 
        length = ja1-jb0
        iaoff = max(jb0-ja0,0) !no offset if we are already inside

        jaj=keyav_c(iaseg0)
        do i=0,length
           apsi_temp=real(apsi_c(jaj+iaoff+i),dp)

           scpr0(1)=scpr0(1)+apsi_temp *real(bpsi1(jbj+i),dp)
           scpr0(2)=scpr0(2)+apsi_temp *real(bpsi2(jbj+i),dp)
           scpr0(3)=scpr0(3)+apsi_temp *real(bpsi3(jbj+i),dp)
           scpr0(4)=scpr0(4)+apsi_temp *real(bpsi4(jbj+i),dp)
           scpr0(5)=scpr0(5)+apsi_temp *real(bpsi5(jbj+i),dp)
           scpr0(6)=scpr0(6)+apsi_temp *real(bpsi6(jbj+i),dp)
           scpr0(7)=scpr0(7)+apsi_temp *real(bpsi7(jbj+i),dp)
           scpr0(8)=scpr0(8)+apsi_temp *real(bpsi8(jbj+i),dp)
           scpr0(9)=scpr0(9)+apsi_temp *real(bpsi9(jbj+i),dp)
           scpr0(10)=scpr0(10)+apsi_temp *real(bpsi10(jbj+i),dp)
           scpr0(11)=scpr0(11)+apsi_temp *real(bpsi11(jbj+i),dp)
           scpr0(12)=scpr0(12)+apsi_temp *real(bpsi12(jbj+i),dp)
           scpr0(13)=scpr0(13)+apsi_temp *real(bpsi13(jbj+i),dp)
           scpr0(14)=scpr0(14)+apsi_temp *real(bpsi14(jbj+i),dp)
           scpr0(15)=scpr0(15)+apsi_temp *real(bpsi15(jbj+i),dp)
           scpr0(16)=scpr0(16)+apsi_temp *real(bpsi16(jbj+i),dp)
           scpr0(17)=scpr0(17)+apsi_temp *real(bpsi17(jbj+i),dp)
           scpr0(18)=scpr0(18)+apsi_temp *real(bpsi18(jbj+i),dp)
           scpr0(19)=scpr0(19)+apsi_temp *real(bpsi19(jbj+i),dp)
           scpr0(20)=scpr0(20)+apsi_temp *real(bpsi20(jbj+i),dp)
           scpr0(21)=scpr0(21)+apsi_temp *real(bpsi21(jbj+i),dp)
           scpr0(22)=scpr0(22)+apsi_temp *real(bpsi22(jbj+i),dp)

        enddo
 !       print *,'length',length,ibseg,scpr0,iaseg0,ja1,jb1

        !print *,'ibseg,mbseg_c,iaseg0,maseg_c',ibseg,mbseg_c,iaseg0,maseg_c
        !print '(a,6(i8),1pe25.17)','ja0,ja1t,ja1,jb0,jb1',&
        !     ibseg,ja0,keyag_c(2,iaseg0),ja1,jb0,jb1,scpr0
        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_c) exit nonconvex_loop_c !segment is finished
        iaseg0=iaseg0+1
        jb0=max(keybg_c(1,ibseg),keyag_c_lin(iaseg0))
        if (keyag_c_lin(iaseg0)>jb1) exit nonconvex_loop_c !segment is not covered
        jbj=jbj+max(jb0-keybg_c(1,ibseg),0)
     end do nonconvex_loop_c
     !disable loop if the end is reached
     if (iaseg0 == maseg_c .and. keybg_c(1,ibseg)> keyag_c_lin(maseg_c)) iaseg0=iaseg0+1


   enddo
!stop
!$omp end do nowait

! fine part

iaseg0=1

!$omp do schedule(static)
   do ibseg=1,mbseg_f
     jbj=keybv_f(ibseg)
     jb0=keybg_f(1,ibseg)
     jb1=keybg_f(2,ibseg)
!    print *,'huntenter',ibseg,jb0,jb1
     call hunt1(.true.,keyag_f_lin,maseg_f,keybg_f(1,ibseg),iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     nonconvex_loop_f: do while(iaseg0 <= maseg_f)
!!$     length = jb1-jb0
!!$     iaoff = jb0-keyag_f_lin(iaseg0)

        ja0=keyag_f_lin(iaseg0) !still doubts about copying in automatic array
        ja1=min(jb1,keyag_f(2,iaseg0)) 
        length = ja1-jb0
        iaoff = max(jb0-ja0,0) !no offset if we are already inside

        jaj=keyav_f(iaseg0)
        do i=0,length
           apsi_f_temp(1)=real(apsi_f(1,jaj+iaoff+i),dp)
           apsi_f_temp(2)=real(apsi_f(2,jaj+iaoff+i),dp)
           apsi_f_temp(3)=real(apsi_f(3,jaj+iaoff+i),dp)
           apsi_f_temp(4)=real(apsi_f(4,jaj+iaoff+i),dp)
           apsi_f_temp(5)=real(apsi_f(5,jaj+iaoff+i),dp)
           apsi_f_temp(6)=real(apsi_f(6,jaj+iaoff+i),dp)
           apsi_f_temp(7)=real(apsi_f(7,jaj+iaoff+i),dp)

           scpr1(1)=scpr1(1)+apsi_f_temp(1)*real(bpsi_f1(1,jbj+i),dp)
           scpr2(1)=scpr2(1)+apsi_f_temp(2)*real(bpsi_f1(2,jbj+i),dp)
           scpr3(1)=scpr3(1)+apsi_f_temp(3)*real(bpsi_f1(3,jbj+i),dp)
           scpr4(1)=scpr4(1)+apsi_f_temp(4)*real(bpsi_f1(4,jbj+i),dp)
           scpr5(1)=scpr5(1)+apsi_f_temp(5)*real(bpsi_f1(5,jbj+i),dp)
           scpr6(1)=scpr6(1)+apsi_f_temp(6)*real(bpsi_f1(6,jbj+i),dp)
           scpr7(1)=scpr7(1)+apsi_f_temp(7)*real(bpsi_f1(7,jbj+i),dp)

           scpr1(2)=scpr1(2)+apsi_f_temp(1)*real(bpsi_f2(1,jbj+i),dp)
           scpr2(2)=scpr2(2)+apsi_f_temp(2)*real(bpsi_f2(2,jbj+i),dp)
           scpr3(2)=scpr3(2)+apsi_f_temp(3)*real(bpsi_f2(3,jbj+i),dp)
           scpr4(2)=scpr4(2)+apsi_f_temp(4)*real(bpsi_f2(4,jbj+i),dp)
           scpr5(2)=scpr5(2)+apsi_f_temp(5)*real(bpsi_f2(5,jbj+i),dp)
           scpr6(2)=scpr6(2)+apsi_f_temp(6)*real(bpsi_f2(6,jbj+i),dp)
           scpr7(2)=scpr7(2)+apsi_f_temp(7)*real(bpsi_f2(7,jbj+i),dp)

           scpr1(3)=scpr1(3)+apsi_f_temp(1)*real(bpsi_f3(1,jbj+i),dp)
           scpr2(3)=scpr2(3)+apsi_f_temp(2)*real(bpsi_f3(2,jbj+i),dp)
           scpr3(3)=scpr3(3)+apsi_f_temp(3)*real(bpsi_f3(3,jbj+i),dp)
           scpr4(3)=scpr4(3)+apsi_f_temp(4)*real(bpsi_f3(4,jbj+i),dp)
           scpr5(3)=scpr5(3)+apsi_f_temp(5)*real(bpsi_f3(5,jbj+i),dp)
           scpr6(3)=scpr6(3)+apsi_f_temp(6)*real(bpsi_f3(6,jbj+i),dp)
           scpr7(3)=scpr7(3)+apsi_f_temp(7)*real(bpsi_f3(7,jbj+i),dp)

           scpr1(4)=scpr1(4)+apsi_f_temp(1)*real(bpsi_f4(1,jbj+i),dp)
           scpr2(4)=scpr2(4)+apsi_f_temp(2)*real(bpsi_f4(2,jbj+i),dp)
           scpr3(4)=scpr3(4)+apsi_f_temp(3)*real(bpsi_f4(3,jbj+i),dp)
           scpr4(4)=scpr4(4)+apsi_f_temp(4)*real(bpsi_f4(4,jbj+i),dp)
           scpr5(4)=scpr5(4)+apsi_f_temp(5)*real(bpsi_f4(5,jbj+i),dp)
           scpr6(4)=scpr6(4)+apsi_f_temp(6)*real(bpsi_f4(6,jbj+i),dp)
           scpr7(4)=scpr7(4)+apsi_f_temp(7)*real(bpsi_f4(7,jbj+i),dp)

           scpr1(5)=scpr1(5)+apsi_f_temp(1)*real(bpsi_f5(1,jbj+i),dp)
           scpr2(5)=scpr2(5)+apsi_f_temp(2)*real(bpsi_f5(2,jbj+i),dp)
           scpr3(5)=scpr3(5)+apsi_f_temp(3)*real(bpsi_f5(3,jbj+i),dp)
           scpr4(5)=scpr4(5)+apsi_f_temp(4)*real(bpsi_f5(4,jbj+i),dp)
           scpr5(5)=scpr5(5)+apsi_f_temp(5)*real(bpsi_f5(5,jbj+i),dp)
           scpr6(5)=scpr6(5)+apsi_f_temp(6)*real(bpsi_f5(6,jbj+i),dp)
           scpr7(5)=scpr7(5)+apsi_f_temp(7)*real(bpsi_f5(7,jbj+i),dp)

           scpr1(6)=scpr1(6)+apsi_f_temp(1)*real(bpsi_f6(1,jbj+i),dp)
           scpr2(6)=scpr2(6)+apsi_f_temp(2)*real(bpsi_f6(2,jbj+i),dp)
           scpr3(6)=scpr3(6)+apsi_f_temp(3)*real(bpsi_f6(3,jbj+i),dp)
           scpr4(6)=scpr4(6)+apsi_f_temp(4)*real(bpsi_f6(4,jbj+i),dp)
           scpr5(6)=scpr5(6)+apsi_f_temp(5)*real(bpsi_f6(5,jbj+i),dp)
           scpr6(6)=scpr6(6)+apsi_f_temp(6)*real(bpsi_f6(6,jbj+i),dp)
           scpr7(6)=scpr7(6)+apsi_f_temp(7)*real(bpsi_f6(7,jbj+i),dp)

           scpr1(7)=scpr1(7)+apsi_f_temp(1)*real(bpsi_f7(1,jbj+i),dp)
           scpr2(7)=scpr2(7)+apsi_f_temp(2)*real(bpsi_f7(2,jbj+i),dp)
           scpr3(7)=scpr3(7)+apsi_f_temp(3)*real(bpsi_f7(3,jbj+i),dp)
           scpr4(7)=scpr4(7)+apsi_f_temp(4)*real(bpsi_f7(4,jbj+i),dp)
           scpr5(7)=scpr5(7)+apsi_f_temp(5)*real(bpsi_f7(5,jbj+i),dp)
           scpr6(7)=scpr6(7)+apsi_f_temp(6)*real(bpsi_f7(6,jbj+i),dp)
           scpr7(7)=scpr7(7)+apsi_f_temp(7)*real(bpsi_f7(7,jbj+i),dp)

           scpr1(8)=scpr1(8)+apsi_f_temp(1)*real(bpsi_f8(1,jbj+i),dp)
           scpr2(8)=scpr2(8)+apsi_f_temp(2)*real(bpsi_f8(2,jbj+i),dp)
           scpr3(8)=scpr3(8)+apsi_f_temp(3)*real(bpsi_f8(3,jbj+i),dp)
           scpr4(8)=scpr4(8)+apsi_f_temp(4)*real(bpsi_f8(4,jbj+i),dp)
           scpr5(8)=scpr5(8)+apsi_f_temp(5)*real(bpsi_f8(5,jbj+i),dp)
           scpr6(8)=scpr6(8)+apsi_f_temp(6)*real(bpsi_f8(6,jbj+i),dp)
           scpr7(8)=scpr7(8)+apsi_f_temp(7)*real(bpsi_f8(7,jbj+i),dp)

           scpr1(9)=scpr1(9)+apsi_f_temp(1)*real(bpsi_f9(1,jbj+i),dp)
           scpr2(9)=scpr2(9)+apsi_f_temp(2)*real(bpsi_f9(2,jbj+i),dp)
           scpr3(9)=scpr3(9)+apsi_f_temp(3)*real(bpsi_f9(3,jbj+i),dp)
           scpr4(9)=scpr4(9)+apsi_f_temp(4)*real(bpsi_f9(4,jbj+i),dp)
           scpr5(9)=scpr5(9)+apsi_f_temp(5)*real(bpsi_f9(5,jbj+i),dp)
           scpr6(9)=scpr6(9)+apsi_f_temp(6)*real(bpsi_f9(6,jbj+i),dp)
           scpr7(9)=scpr7(9)+apsi_f_temp(7)*real(bpsi_f9(7,jbj+i),dp)

           scpr1(10)=scpr1(10)+apsi_f_temp(1)*real(bpsi_f10(1,jbj+i),dp)
           scpr2(10)=scpr2(10)+apsi_f_temp(2)*real(bpsi_f10(2,jbj+i),dp)
           scpr3(10)=scpr3(10)+apsi_f_temp(3)*real(bpsi_f10(3,jbj+i),dp)
           scpr4(10)=scpr4(10)+apsi_f_temp(4)*real(bpsi_f10(4,jbj+i),dp)
           scpr5(10)=scpr5(10)+apsi_f_temp(5)*real(bpsi_f10(5,jbj+i),dp)
           scpr6(10)=scpr6(10)+apsi_f_temp(6)*real(bpsi_f10(6,jbj+i),dp)
           scpr7(10)=scpr7(10)+apsi_f_temp(7)*real(bpsi_f10(7,jbj+i),dp)

           scpr1(11)=scpr1(11)+apsi_f_temp(1)*real(bpsi_f11(1,jbj+i),dp)
           scpr2(11)=scpr2(11)+apsi_f_temp(2)*real(bpsi_f11(2,jbj+i),dp)
           scpr3(11)=scpr3(11)+apsi_f_temp(3)*real(bpsi_f11(3,jbj+i),dp)
           scpr4(11)=scpr4(11)+apsi_f_temp(4)*real(bpsi_f11(4,jbj+i),dp)
           scpr5(11)=scpr5(11)+apsi_f_temp(5)*real(bpsi_f11(5,jbj+i),dp)
           scpr6(11)=scpr6(11)+apsi_f_temp(6)*real(bpsi_f11(6,jbj+i),dp)
           scpr7(11)=scpr7(11)+apsi_f_temp(7)*real(bpsi_f11(7,jbj+i),dp)

           scpr1(12)=scpr1(12)+apsi_f_temp(1)*real(bpsi_f12(1,jbj+i),dp)
           scpr2(12)=scpr2(12)+apsi_f_temp(2)*real(bpsi_f12(2,jbj+i),dp)
           scpr3(12)=scpr3(12)+apsi_f_temp(3)*real(bpsi_f12(3,jbj+i),dp)
           scpr4(12)=scpr4(12)+apsi_f_temp(4)*real(bpsi_f12(4,jbj+i),dp)
           scpr5(12)=scpr5(12)+apsi_f_temp(5)*real(bpsi_f12(5,jbj+i),dp)
           scpr6(12)=scpr6(12)+apsi_f_temp(6)*real(bpsi_f12(6,jbj+i),dp)
           scpr7(12)=scpr7(12)+apsi_f_temp(7)*real(bpsi_f12(7,jbj+i),dp)

           scpr1(13)=scpr1(13)+apsi_f_temp(1)*real(bpsi_f13(1,jbj+i),dp)
           scpr2(13)=scpr2(13)+apsi_f_temp(2)*real(bpsi_f13(2,jbj+i),dp)
           scpr3(13)=scpr3(13)+apsi_f_temp(3)*real(bpsi_f13(3,jbj+i),dp)
           scpr4(13)=scpr4(13)+apsi_f_temp(4)*real(bpsi_f13(4,jbj+i),dp)
           scpr5(13)=scpr5(13)+apsi_f_temp(5)*real(bpsi_f13(5,jbj+i),dp)
           scpr6(13)=scpr6(13)+apsi_f_temp(6)*real(bpsi_f13(6,jbj+i),dp)
           scpr7(13)=scpr7(13)+apsi_f_temp(7)*real(bpsi_f13(7,jbj+i),dp)

           scpr1(14)=scpr1(14)+apsi_f_temp(1)*real(bpsi_f14(1,jbj+i),dp)
           scpr2(14)=scpr2(14)+apsi_f_temp(2)*real(bpsi_f14(2,jbj+i),dp)
           scpr3(14)=scpr3(14)+apsi_f_temp(3)*real(bpsi_f14(3,jbj+i),dp)
           scpr4(14)=scpr4(14)+apsi_f_temp(4)*real(bpsi_f14(4,jbj+i),dp)
           scpr5(14)=scpr5(14)+apsi_f_temp(5)*real(bpsi_f14(5,jbj+i),dp)
           scpr6(14)=scpr6(14)+apsi_f_temp(6)*real(bpsi_f14(6,jbj+i),dp)
           scpr7(14)=scpr7(14)+apsi_f_temp(7)*real(bpsi_f14(7,jbj+i),dp)

           scpr1(15)=scpr1(15)+apsi_f_temp(1)*real(bpsi_f15(1,jbj+i),dp)
           scpr2(15)=scpr2(15)+apsi_f_temp(2)*real(bpsi_f15(2,jbj+i),dp)
           scpr3(15)=scpr3(15)+apsi_f_temp(3)*real(bpsi_f15(3,jbj+i),dp)
           scpr4(15)=scpr4(15)+apsi_f_temp(4)*real(bpsi_f15(4,jbj+i),dp)
           scpr5(15)=scpr5(15)+apsi_f_temp(5)*real(bpsi_f15(5,jbj+i),dp)
           scpr6(15)=scpr6(15)+apsi_f_temp(6)*real(bpsi_f15(6,jbj+i),dp)
           scpr7(15)=scpr7(15)+apsi_f_temp(7)*real(bpsi_f15(7,jbj+i),dp)

           scpr1(16)=scpr1(16)+apsi_f_temp(1)*real(bpsi_f16(1,jbj+i),dp)
           scpr2(16)=scpr2(16)+apsi_f_temp(2)*real(bpsi_f16(2,jbj+i),dp)
           scpr3(16)=scpr3(16)+apsi_f_temp(3)*real(bpsi_f16(3,jbj+i),dp)
           scpr4(16)=scpr4(16)+apsi_f_temp(4)*real(bpsi_f16(4,jbj+i),dp)
           scpr5(16)=scpr5(16)+apsi_f_temp(5)*real(bpsi_f16(5,jbj+i),dp)
           scpr6(16)=scpr6(16)+apsi_f_temp(6)*real(bpsi_f16(6,jbj+i),dp)
           scpr7(16)=scpr7(16)+apsi_f_temp(7)*real(bpsi_f16(7,jbj+i),dp)

           scpr1(17)=scpr1(17)+apsi_f_temp(1)*real(bpsi_f17(1,jbj+i),dp)
           scpr2(17)=scpr2(17)+apsi_f_temp(2)*real(bpsi_f17(2,jbj+i),dp)
           scpr3(17)=scpr3(17)+apsi_f_temp(3)*real(bpsi_f17(3,jbj+i),dp)
           scpr4(17)=scpr4(17)+apsi_f_temp(4)*real(bpsi_f17(4,jbj+i),dp)
           scpr5(17)=scpr5(17)+apsi_f_temp(5)*real(bpsi_f17(5,jbj+i),dp)
           scpr6(17)=scpr6(17)+apsi_f_temp(6)*real(bpsi_f17(6,jbj+i),dp)
           scpr7(17)=scpr7(17)+apsi_f_temp(7)*real(bpsi_f17(7,jbj+i),dp)

           scpr1(18)=scpr1(18)+apsi_f_temp(1)*real(bpsi_f18(1,jbj+i),dp)
           scpr2(18)=scpr2(18)+apsi_f_temp(2)*real(bpsi_f18(2,jbj+i),dp)
           scpr3(18)=scpr3(18)+apsi_f_temp(3)*real(bpsi_f18(3,jbj+i),dp)
           scpr4(18)=scpr4(18)+apsi_f_temp(4)*real(bpsi_f18(4,jbj+i),dp)
           scpr5(18)=scpr5(18)+apsi_f_temp(5)*real(bpsi_f18(5,jbj+i),dp)
           scpr6(18)=scpr6(18)+apsi_f_temp(6)*real(bpsi_f18(6,jbj+i),dp)
           scpr7(18)=scpr7(18)+apsi_f_temp(7)*real(bpsi_f18(7,jbj+i),dp)

           scpr1(19)=scpr1(19)+apsi_f_temp(1)*real(bpsi_f19(1,jbj+i),dp)
           scpr2(19)=scpr2(19)+apsi_f_temp(2)*real(bpsi_f19(2,jbj+i),dp)
           scpr3(19)=scpr3(19)+apsi_f_temp(3)*real(bpsi_f19(3,jbj+i),dp)
           scpr4(19)=scpr4(19)+apsi_f_temp(4)*real(bpsi_f19(4,jbj+i),dp)
           scpr5(19)=scpr5(19)+apsi_f_temp(5)*real(bpsi_f19(5,jbj+i),dp)
           scpr6(19)=scpr6(19)+apsi_f_temp(6)*real(bpsi_f19(6,jbj+i),dp)
           scpr7(19)=scpr7(19)+apsi_f_temp(7)*real(bpsi_f19(7,jbj+i),dp)

           scpr1(20)=scpr1(20)+apsi_f_temp(1)*real(bpsi_f20(1,jbj+i),dp)
           scpr2(20)=scpr2(20)+apsi_f_temp(2)*real(bpsi_f20(2,jbj+i),dp)
           scpr3(20)=scpr3(20)+apsi_f_temp(3)*real(bpsi_f20(3,jbj+i),dp)
           scpr4(20)=scpr4(20)+apsi_f_temp(4)*real(bpsi_f20(4,jbj+i),dp)
           scpr5(20)=scpr5(20)+apsi_f_temp(5)*real(bpsi_f20(5,jbj+i),dp)
           scpr6(20)=scpr6(20)+apsi_f_temp(6)*real(bpsi_f20(6,jbj+i),dp)
           scpr7(20)=scpr7(20)+apsi_f_temp(7)*real(bpsi_f20(7,jbj+i),dp)

           scpr1(21)=scpr1(21)+apsi_f_temp(1)*real(bpsi_f21(1,jbj+i),dp)
           scpr2(21)=scpr2(21)+apsi_f_temp(2)*real(bpsi_f21(2,jbj+i),dp)
           scpr3(21)=scpr3(21)+apsi_f_temp(3)*real(bpsi_f21(3,jbj+i),dp)
           scpr4(21)=scpr4(21)+apsi_f_temp(4)*real(bpsi_f21(4,jbj+i),dp)
           scpr5(21)=scpr5(21)+apsi_f_temp(5)*real(bpsi_f21(5,jbj+i),dp)
           scpr6(21)=scpr6(21)+apsi_f_temp(6)*real(bpsi_f21(6,jbj+i),dp)
           scpr7(21)=scpr7(21)+apsi_f_temp(7)*real(bpsi_f21(7,jbj+i),dp)

           scpr1(22)=scpr1(22)+apsi_f_temp(1)*real(bpsi_f22(1,jbj+i),dp)
           scpr2(22)=scpr2(22)+apsi_f_temp(2)*real(bpsi_f22(2,jbj+i),dp)
           scpr3(22)=scpr3(22)+apsi_f_temp(3)*real(bpsi_f22(3,jbj+i),dp)
           scpr4(22)=scpr4(22)+apsi_f_temp(4)*real(bpsi_f22(4,jbj+i),dp)
           scpr5(22)=scpr5(22)+apsi_f_temp(5)*real(bpsi_f22(5,jbj+i),dp)
           scpr6(22)=scpr6(22)+apsi_f_temp(6)*real(bpsi_f22(6,jbj+i),dp)
           scpr7(22)=scpr7(22)+apsi_f_temp(7)*real(bpsi_f22(7,jbj+i),dp)

        enddo
 !       print *,'length',length,ibseg,scpr1,iaseg0,ja1,jb1
        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_f) exit nonconvex_loop_f !segment is finished  
        iaseg0=iaseg0+1
        jb0=max(keybg_f(1,ibseg),keyag_f_lin(iaseg0))
        if (keyag_f_lin(iaseg0)>jb1) exit nonconvex_loop_f !segment is not covered 
        jbj=jbj+max(jb0-keybg_f(1,ibseg),0)
     end do nonconvex_loop_f
     !disable loop if the end is reached
     if (iaseg0 == maseg_f .and. keybg_f(1,ibseg)> keyag_f_lin(maseg_f)) iaseg0=iaseg0+1

   enddo
!$omp end do !!!implicit barrier 

   scpr0=scpr0+scpr1+scpr2+scpr3+scpr4+scpr5+scpr6+scpr7

!$omp critical 
   scpr=scpr+scpr0
!$omp end critical
  
!$omp end parallel

!!!    call system_clock(ncount2,ncount_rate,ncount_max)
!!!    tel=dble(ncount2-ncount0)/dble(ncount_rate)
!!!    write(97,*) 'wpdot:',tel
!!!    close(97)

END SUBROUTINE wpdot_22

!> Wrapper of wpdot to avoid boundary problems in absence of wavelets
!! and to perform scalar product for complex wavefunctions and projectors
!! if the wavefunctions are complex, so should be also the projectors
subroutine wpdot_wrap(ncplx,mavctr_c,mavctr_f,maseg_c,maseg_f,keyav,keyag,apsi,  &
     mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keybv,keybg,bpsi,scpr)
  use module_base
  implicit none
  integer, intent(in) :: mavctr_c,mavctr_f,maseg_c,maseg_f
  integer, intent(in) :: ncplx,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f
  integer, dimension(maseg_c+maseg_f), intent(in) :: keyav
  integer, dimension(mbseg_c+mbseg_f), intent(in) :: keybv
  integer, dimension(2,maseg_c+maseg_f), intent(in) :: keyag
  integer, dimension(2,mbseg_c+mbseg_f), intent(in) :: keybg
  real(wp), dimension(mavctr_c+7*mavctr_f,ncplx), intent(in) :: apsi
  real(wp), dimension(mbvctr_c+7*mbvctr_f,ncplx), intent(in) :: bpsi
  real(dp), dimension(ncplx), intent(out) :: scpr
  !local variables
  integer :: ia_f,ib_f,iaseg_f,ibseg_f,ia,ib
  real(dp), dimension(ncplx,ncplx) :: scalprod 

  ia_f=min(mavctr_f,1)
  ib_f=min(mbvctr_f,1)

  iaseg_f=min(maseg_f,1)
  ibseg_f=min(mbseg_f,1)

  do ia=1,ncplx
     do ib=1,ncplx
        call wpdot(mavctr_c,mavctr_f,maseg_c,maseg_f,&
             keyav,keyav(maseg_c+iaseg_f),&
             keyag,keyag(1,maseg_c+iaseg_f),&
             apsi(1,ia),apsi(mavctr_c+ia_f,ia),  &
             mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
             keybv,keybv(mbseg_c+ibseg_f),&
             keybg,keybg(1,mbseg_c+ibseg_f),&
             bpsi(1,ib),bpsi(mbvctr_c+ib_f,ib),scalprod(ia,ib))
     end do
  end do

  !then define the result
  if (ncplx == 1) then
     scpr(1)=scalprod(1,1)
  else if (ncplx == 2) then
     scpr(1)=scalprod(1,1)+scalprod(2,2)
     scpr(2)=scalprod(1,2)-scalprod(2,1)
  else
     write(*,*)'ERROR wpdot: ncplx not valid:',ncplx
     stop
  end if

END SUBROUTINE wpdot_wrap

!> This function must be generalized for the linear scaling code                
!! @warning   
!! calculates the dot product between a wavefunctions apsi and a projector bpsi (both in compressed form)
!! Warning: It is assumed that the segments of bpsi are always contained inside
!! the segments of apsi.
subroutine wpdot(  &
     mavctr_c,mavctr_f,maseg_c,maseg_f,keyav_c,keyav_f,keyag_c,keyag_f,apsi_c,apsi_f,  &
     mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keybv_c,keybv_f,keybg_c,keybg_f,bpsi_c,bpsi_f,scpr)
  use module_base
  implicit none
  integer, intent(in) :: mavctr_c,mavctr_f,maseg_c,maseg_f,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f
  integer, dimension(maseg_c), intent(in) :: keyav_c
  integer, dimension(mbseg_c), intent(in) :: keybv_c
  integer, dimension(maseg_f), intent(in) :: keyav_f
  integer, dimension(mbseg_f), intent(in) :: keybv_f
  integer, dimension(2,maseg_c), intent(in) :: keyag_c
  integer, dimension(2,mbseg_c), intent(in) :: keybg_c
  integer, dimension(2,maseg_f), intent(in) :: keyag_f
  integer, dimension(2,mbseg_f), intent(in) :: keybg_f
  real(wp), dimension(mavctr_c), intent(in) :: apsi_c
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi_c
  real(wp), dimension(7,mavctr_f), intent(in) :: apsi_f
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f
  real(dp), intent(out) :: scpr
  !local variables
  integer :: ibseg,jaj,jb1,jb0,jbj,iaoff,iboff,length,i,ja0,ja1
  real(dp) :: scpr1,scpr2,scpr3,scpr4,scpr5,scpr6,scpr7,scpr0
  integer :: iaseg0
  integer, dimension(maseg_c) :: keyag_c_lin !>linear version of second indices of keyag_c
  integer, dimension(maseg_f) :: keyag_f_lin !>linear version of second indices of keyag_f
!!!    integer :: ncount0,ncount2,ncount_rate,ncount_max
!!!    real(gp) :: tel


!!!  !dee
!!!    open(unit=97,file='time_wpdot',status='unknown',position='append')
!!!    call system_clock(ncount0,ncount_rate,ncount_max)

  keyag_c_lin = keyag_c(1,:)!speed up access in hunt subroutine by consecutive arrangement in memory
  keyag_f_lin = keyag_f(1,:)!speed up access in hunt subroutine by consecutive arrangement in memory

  scpr=0.0_dp
  scpr0=0.0_dp
  scpr1=0.0_dp
  scpr2=0.0_dp
  scpr3=0.0_dp
  scpr4=0.0_dp
  scpr5=0.0_dp
  scpr6=0.0_dp
  scpr7=0.0_dp

!$omp parallel default (private) firstprivate(scpr0,scpr1,scpr2,scpr3,scpr4,scpr5,scpr6,scpr7)&
!$omp shared (maseg_c,keyav_c,keyag_c,keyag_c_lin,keybg_c,mbseg_c,keybv_c,mbseg_f,maseg_f)&
!$omp shared (apsi_c,bpsi_c,bpsi_f,keybv_f,keybg_f,keyag_f,keyag_f_lin,keyav_f)&
!$omp shared (apsi_f,scpr)
!!!!$omp shared (ncount0,ncount2,ncount_rate,ncount_max,tel)

  iaseg0=1 

!coarse part. Loop on the projectors segments
!$omp do schedule(static)
   do ibseg=1,mbseg_c
     jbj=keybv_c(ibseg)
!     jb0=keybg_c(1,ibseg) !starting point of projector segment
     jb0=max(keybg_c(1,ibseg),keyag_c_lin(1))
     jb1=keybg_c(2,ibseg) !ending point of projector segment
     iboff = max(jb0-keybg_c(1,ibseg),0)
     !print *,'huntenter',ibseg,jb0,jb1
 
     !find the starting point of the wavefunction segment
     !warning: hunt is assuming that the variable is always found
     !if it is not, iaseg0 is put to maseg + 1 so that the loop is disabled
!     call hunt1(.true.,keyag_c_lin,maseg_c,keybg_c(1,ibseg),iaseg0)
     call hunt1(.true.,keyag_c_lin,maseg_c,jb0,iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     !now pass through all the wavefunction segments until the end of the segment is 
     !still contained in projector segment
     nonconvex_loop_c: do while(iaseg0 <= maseg_c)
!     print *,'huntexit',iaseg0,maseg_c,keyag_c_lin(iaseg0),keyag_c(2,iaseg0)

        !length = jb1-jb0
        !iaoff = jb0-keyag_c_lin(iaseg0)!jb0-ja0

        ja0=keyag_c_lin(iaseg0)
        ja1=min(jb1,keyag_c(2,iaseg0)) 
        length = ja1-jb0
        iaoff = max(jb0-ja0,0) !no offset if we are already inside

        jaj=keyav_c(iaseg0)
        do i=0,length
           scpr0=scpr0+real(apsi_c(jaj+iaoff+i),dp)*real(bpsi_c(jbj+i+iboff),dp)
        enddo
       !print *,'length',length,ibseg,scpr0,iaseg0,ja1,jb1

        !print *,'ibseg,mbseg_c,iaseg0,maseg_c',ibseg,mbseg_c,iaseg0,maseg_c
        !print '(a,6(i8),1pe25.17)','ja0,ja1t,ja1,jb0,jb1',&
        !     ibseg,ja0,keyag_c(2,iaseg0),ja1,jb0,jb1,scpr0
        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_c) exit nonconvex_loop_c !segment is finished
        iaseg0=iaseg0+1
        jb0=max(keybg_c(1,ibseg),keyag_c_lin(iaseg0))
        if (keyag_c_lin(iaseg0)>jb1) exit nonconvex_loop_c !segment is not covered
        jbj=jbj+max(jb0-keybg_c(1,ibseg),0)
     end do nonconvex_loop_c
     !disable loop if the end is reached
     if (iaseg0 == maseg_c .and. keybg_c(1,ibseg)> keyag_c_lin(maseg_c)) iaseg0=iaseg0+1


   enddo
!stop
!$omp end do nowait

! fine part

iaseg0=1

!$omp do schedule(static)
   do ibseg=1,mbseg_f
     jbj=keybv_f(ibseg)
     !jb0=keybg_f(1,ibseg)
     jb0=max(keybg_f(1,ibseg),keyag_f_lin(1))
     jb1=keybg_f(2,ibseg)
     iboff = max(jb0-keybg_f(1,ibseg),0)
!    print *,'huntenter',ibseg,jb0,jb1
     !call hunt1(.true.,keyag_f_lin,maseg_f,keybg_f(1,ibseg),iaseg0)
     call hunt1(.true.,keyag_f_lin,maseg_f,jb0,iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     nonconvex_loop_f: do while(iaseg0 <= maseg_f)
!!$     length = jb1-jb0
!!$     iaoff = jb0-keyag_f_lin(iaseg0)

        ja0=keyag_f_lin(iaseg0) !still doubts about copying in automatic array
        ja1=min(jb1,keyag_f(2,iaseg0)) 
        length = ja1-jb0
        iaoff = max(jb0-ja0,0) !no offset if we are already inside

        jaj=keyav_f(iaseg0)
        do i=0,length
           scpr1=scpr1+real(apsi_f(1,jaj+iaoff+i),dp)*real(bpsi_f(1,jbj+i+iboff),dp)
           scpr2=scpr2+real(apsi_f(2,jaj+iaoff+i),dp)*real(bpsi_f(2,jbj+i+iboff),dp)
           scpr3=scpr3+real(apsi_f(3,jaj+iaoff+i),dp)*real(bpsi_f(3,jbj+i+iboff),dp)
           scpr4=scpr4+real(apsi_f(4,jaj+iaoff+i),dp)*real(bpsi_f(4,jbj+i+iboff),dp)
           scpr5=scpr5+real(apsi_f(5,jaj+iaoff+i),dp)*real(bpsi_f(5,jbj+i+iboff),dp)
           scpr6=scpr6+real(apsi_f(6,jaj+iaoff+i),dp)*real(bpsi_f(6,jbj+i+iboff),dp)
           scpr7=scpr7+real(apsi_f(7,jaj+iaoff+i),dp)*real(bpsi_f(7,jbj+i+iboff),dp)
        enddo
 !       print *,'length',length,ibseg,scpr1,iaseg0,ja1,jb1
        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_f) exit nonconvex_loop_f !segment is finished  
        iaseg0=iaseg0+1
        jb0=max(keybg_f(1,ibseg),keyag_f_lin(iaseg0))
        if (keyag_f_lin(iaseg0)>jb1) exit nonconvex_loop_f !segment is not covered 
        jbj=jbj+max(jb0-keybg_f(1,ibseg),0)
     end do nonconvex_loop_f
     !disable loop if the end is reached
     if (iaseg0 == maseg_f .and. keybg_f(1,ibseg)> keyag_f_lin(maseg_f)) iaseg0=iaseg0+1

   enddo
!$omp end do !!!implicit barrier 

   scpr0=scpr0+scpr1+scpr2+scpr3+scpr4+scpr5+scpr6+scpr7

!$omp critical 
   scpr=scpr+scpr0
!$omp end critical

!$omp end parallel

!!!    call system_clock(ncount2,ncount_rate,ncount_max)
!!!    tel=dble(ncount2-ncount0)/dble(ncount_rate)
!!!    write(97,*) 'wpdot:',tel
!!!    close(97)

END SUBROUTINE wpdot


!> Wrapper of waxpy for complex Ax+y and for no fine grid cases
!! @warning 
!!    in complex cases, it acts with y = Conj(A) *x +y, with the complex conjugate
!!    if the a function is complex, so should be the b function
subroutine waxpy_wrap(ncplx,scpr,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keybv,keybg,bpsi,&
     mavctr_c,mavctr_f,maseg_c,maseg_f,keyav,keyag,apsi)
  use module_base
  implicit none
  integer, intent(in) :: mavctr_c,mavctr_f,maseg_c,maseg_f
  integer, intent(in) :: ncplx,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f
  real(dp), dimension(ncplx), intent(in) :: scpr
  integer, dimension(maseg_c+maseg_f), intent(in) :: keyav
  integer, dimension(mbseg_c+mbseg_f), intent(in) :: keybv
  integer, dimension(2,maseg_c+maseg_f), intent(in) :: keyag
  integer, dimension(2,mbseg_c+mbseg_f), intent(in) :: keybg
  real(wp), dimension(mbvctr_c+7*mbvctr_f,ncplx), intent(in) :: bpsi
  real(wp), dimension(mavctr_c+7*mavctr_f,ncplx), intent(inout) :: apsi
  !local variables
  integer :: ia_f,ib_f,iaseg_f,ibseg_f,ia,ib,is

  ia_f=min(mavctr_f,1)
  ib_f=min(mbvctr_f,1)

  iaseg_f=min(maseg_f,1)
  ibseg_f=min(mbseg_f,1)


  ia=1
  ib=1
  is=1
  call waxpy(scpr(is),mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
       keybv,keybv(mbseg_c+ibseg_f),&
       keybg,keybg(1,mbseg_c+ibseg_f),&
       bpsi(1,ib),bpsi(mbvctr_c+ib_f,ib), &
       mavctr_c,mavctr_f,maseg_c,maseg_f,&
       keyav,keyav(maseg_c+iaseg_f),&
       keyag,keyag(1,maseg_c+iaseg_f),&
       apsi(1,ia),apsi(mavctr_c+ia_f,ia))
   if (ncplx == 2) then
      ib=2
      is=2

      call waxpy(scpr(is),mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
           keybv,keybv(mbseg_c+ibseg_f),&
           keybg,keybg(1,mbseg_c+ibseg_f),&
           bpsi(1,ib),bpsi(mbvctr_c+ib_f,ib), &
           mavctr_c,mavctr_f,maseg_c,maseg_f,&
           keyav,keyav(maseg_c+iaseg_f),&
           keyag,keyag(1,maseg_c+iaseg_f),&
           apsi(1,ia),apsi(mavctr_c+ia_f,ia))

      ia=2
      is=1
      ib=2

      call waxpy(scpr(is),mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
           keybv,keybv(mbseg_c+ibseg_f),&
           keybg,keybg(1,mbseg_c+ibseg_f),&
           bpsi(1,ib),bpsi(mbvctr_c+ib_f,ib), &
           mavctr_c,mavctr_f,maseg_c,maseg_f,&
           keyav,keyav(maseg_c+iaseg_f),&
           keyag,keyag(1,maseg_c+iaseg_f),&
           apsi(1,ia),apsi(mavctr_c+ia_f,ia))


      is=2
      ib=1
      !beware the minus sign
      call waxpy(-scpr(is),mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
           keybv,keybv(mbseg_c+ibseg_f),&
           keybg,keybg(1,mbseg_c+ibseg_f),&
           bpsi(1,ib),bpsi(mbvctr_c+ib_f,ib), &
           mavctr_c,mavctr_f,maseg_c,maseg_f,&
           keyav,keyav(maseg_c+iaseg_f),&
           keyag,keyag(1,maseg_c+iaseg_f),&
           apsi(1,ia),apsi(mavctr_c+ia_f,ia))
      
   end if



END SUBROUTINE waxpy_wrap


!> Rank 1 update of wavefunction a with wavefunction b: apsi=apsi+scpr*bpsi
!! The update is only done in the localization region of apsi
subroutine waxpy(  & 
     scpr,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keybv_c,keybv_f,keybg_c,keybg_f,bpsi_c,bpsi_f, & 
     mavctr_c,mavctr_f,maseg_c,maseg_f,keyav_c,keyav_f,keyag_c,keyag_f,apsi_c,apsi_f)
  use module_base
  implicit none
  integer, intent(in) :: mavctr_c,mavctr_f,maseg_c,maseg_f,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f
  real(dp), intent(in) :: scpr
  integer, dimension(maseg_c), intent(in) :: keyav_c
  integer, dimension(mbseg_c), intent(in) :: keybv_c
  integer, dimension(maseg_f), intent(in) :: keyav_f
  integer, dimension(mbseg_f), intent(in) :: keybv_f
  integer, dimension(2,maseg_c), intent(in) :: keyag_c
  integer, dimension(2,mbseg_c), intent(in) :: keybg_c
  integer, dimension(2,maseg_f), intent(in) :: keyag_f
  integer, dimension(2,mbseg_f), intent(in) :: keybg_f
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi_c
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f
  real(wp), dimension(mavctr_c), intent(inout) :: apsi_c
  real(wp), dimension(7,mavctr_f), intent(inout) :: apsi_f
  !local variables
  integer :: ibseg,iaseg0,jaj,jb1,jb0,jbj,iaoff,length,i,ja0,ja1
  real(wp) :: scprwp
  integer, dimension(maseg_c) :: keyag_c_lin!linear version of second inidces of keyag_c
  integer, dimension(maseg_f) :: keyag_f_lin!linear version of second inidces of keyag_f
!!!    integer :: ncount0,ncount2,ncount_rate,ncount_max
!!!    real(gp) :: tel 

!!!!dee
!!!    open(unit=97,file='time_waxpy',status='unknown',position='append')
!!!    call system_clock(ncount0,ncount_rate,ncount_max)

    keyag_c_lin = keyag_c(1,:)!speed up access in hunt subroutine by consecutive arrangement in memory
    keyag_f_lin = keyag_f(1,:)!speed up access in hunt subroutine by consecutive arrangement in memory


  scprwp=real(scpr,wp)

!$omp parallel default (private) &
!$omp shared (maseg_c,keyav_c,keyag_c,keyag_c_lin,keybg_c,mbseg_c,mbseg_f,maseg_f)&
!$omp shared (keyav_f,keyag_f,keyag_f_lin,keybg_f,keybv_f,scprwp,bpsi_c,bpsi_f)&
!$omp shared (apsi_f,apsi_c,keybv_c)
!!!!$omp shared (ncount0,ncount2,ncount_rate,ncount_max,tel)
  
iaseg0=1

! coarse part

!$omp do schedule(static) 
   do ibseg=1,mbseg_c
     jbj=keybv_c(ibseg)
     jb0=keybg_c(1,ibseg)
     jb1=keybg_c(2,ibseg)
    
     call hunt1(.true.,keyag_c_lin,maseg_c,keybg_c(1,ibseg),iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
 
     !now pass through all the wavefunction segments until the end of the segment is 
     !still contained in projector segment
     nonconvex_loop_c: do while(iaseg0 <= maseg_c)
        !length = jb1-jb0
        !iaoff = jb0-keyag_c_lin(iaseg0)!jb0-ja0

        ja0=keyag_c_lin(iaseg0) !still doubts about copying in automatic array
        ja1=min(jb1,keyag_c(2,iaseg0)) 
        length = ja1-jb0
        iaoff = max(jb0-ja0,0) !no offset if we are already inside

        jaj=keyav_c(iaseg0)

        do i=0,length
           apsi_c(jaj+iaoff+i)=apsi_c(jaj+iaoff+i)+scprwp*bpsi_c(jbj+i)
        enddo
        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_c) exit nonconvex_loop_c !segment is finished
        iaseg0=iaseg0+1
        jb0=max(keybg_c(1,ibseg),keyag_c_lin(iaseg0))
        if (keyag_c_lin(iaseg0)>jb1) exit nonconvex_loop_c
        jbj=jbj+max(jb0-keybg_c(1,ibseg),0)
     end do nonconvex_loop_c
     !disable loop if the end is reached
     if (iaseg0 == maseg_c .and. keybg_c(1,ibseg)> keyag_c_lin(maseg_c)) iaseg0=iaseg0+1

   enddo
!$omp end do nowait

! fine part

   iaseg0=1

!$omp do schedule(static)
   do ibseg=1,mbseg_f
      jbj=keybv_f(ibseg)
      jb0=keybg_f(1,ibseg)
      jb1=keybg_f(2,ibseg)

      call hunt1(.true.,keyag_f_lin,maseg_f,keybg_f(1,ibseg),iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
      nonconvex_loop_f: do while(iaseg0 <= maseg_f)
!!$     length = jb1-jb0
!!$     iaoff = jb0-keyag_f_lin(iaseg0)

         ja0=keyag_f_lin(iaseg0) !still doubts about copying in automatic array
         ja1=min(jb1,keyag_f(2,iaseg0)) 
         length = ja1-jb0
         iaoff = max(jb0-ja0,0) !no offset if we are already inside

         jaj=keyav_f(iaseg0)
         do i=0,length
            apsi_f(1,jaj+iaoff+i)=apsi_f(1,jaj+iaoff+i)+&
                 scprwp*bpsi_f(1,jbj+i)
            apsi_f(2,jaj+iaoff+i)=apsi_f(2,jaj+iaoff+i)+&
                 scprwp*bpsi_f(2,jbj+i)
            apsi_f(3,jaj+iaoff+i)=apsi_f(3,jaj+iaoff+i)+&
                 scprwp*bpsi_f(3,jbj+i)
            apsi_f(4,jaj+iaoff+i)=apsi_f(4,jaj+iaoff+i)+&
                 scprwp*bpsi_f(4,jbj+i)
            apsi_f(5,jaj+iaoff+i)=apsi_f(5,jaj+iaoff+i)+&
                 scprwp*bpsi_f(5,jbj+i)
            apsi_f(6,jaj+iaoff+i)=apsi_f(6,jaj+iaoff+i)+&
                 scprwp*bpsi_f(6,jbj+i)
            apsi_f(7,jaj+iaoff+i)=apsi_f(7,jaj+iaoff+i)+&
                 scprwp*bpsi_f(7,jbj+i)
         enddo
        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_f) exit nonconvex_loop_f !segment is finished  
        iaseg0=iaseg0+1
        jb0=max(keybg_f(1,ibseg),keyag_f_lin(iaseg0))
        if (keyag_f_lin(iaseg0)>jb1) exit nonconvex_loop_f !segment is not covered 
        jbj=jbj+max(jb0-keybg_f(1,ibseg),0)
      end do nonconvex_loop_f
      !disable loop if the end is reached
      if (iaseg0 == maseg_f .and. keybg_f(1,ibseg)> keyag_f_lin(maseg_f)) iaseg0=iaseg0+1


   enddo
!$omp end do
!$omp end parallel

!!!    call system_clock(ncount2,ncount_rate,ncount_max)
!!!    tel=dble(ncount2-ncount0)/dble(ncount_rate)
!!!    write(97,*) 'waxpy:',tel    
!!!    close(97)

END SUBROUTINE waxpy

!> Search the segments which intersect each other
!! @todo modify this routine to have also the end as result
subroutine hunt1(ascnd,xx,n,x,jlo)
  implicit none
  logical, intent(in) :: ascnd
  integer, intent(in) :: x !<starting point in grid coordinates
  integer, intent(in) :: n !<number of segments
  integer, dimension(n), intent(in) :: xx !<array of segment starting points
  integer, intent(inout) :: jlo !<input: starting segment, 
                                ! output: closest segment corresponding to x
                                ! warning: if jlo is outside range, routine is disabled
  !local variables
  integer :: inc,jhi,jm

!  print *,'jlo,n,xx(n),x',jlo,n,xx(n),x

  !check array extremes
  if (ascnd) then
     if (jlo > n) return
     !if (x > xx(n)) then
     !   print *,'that is the last'
     !   jlo=n+1
     !   return
     !end if
  else
     if (jlo < 1) return
     !if (x < xx(1)) then
     !   jlo=n+1
     !   return
     !end if
  end if

  !start searching
  if(x == xx(1))then
     jlo=1
     return
  end if
  if(x == xx(n)) then 
     jlo=n
     return
  end if

!  print *,'quickreturn'
 
  !check if the order is ascending (the sense of the ordering)
  !ascnd=xx(n) >= xx(1) now passed as an argument
  
  !the starting point is external to the array, complete search (commented out)
!!$  if(jlo <= 0 .or. jlo > n)then
!!$     jlo=0
!!$     jhi=n+1
!!$     goto 3
!!$  endif

  !increment of the segment
  inc=1
  !target is above starting point
  if ((x >= xx(jlo)) .eqv. ascnd) then
     guess_end: do
        jhi=jlo+inc
        !number of segments is over
        if(jhi > n)then
           jhi=n+1
           exit guess_end
        !increase until the target is below
        else if((x >= xx(jhi)) .eqv. ascnd)then
           jlo=jhi
           inc=inc+inc
        else
           exit guess_end
        endif
     end do guess_end
!print *,'range',jlo,inc,jhi,x,xx(jlo),xx(jhi)
  else
     !target is below, invert start and end
     jhi=jlo
     guess_start: do
        jlo=jhi-inc
        !segment are over (from below)
        if (jlo < 1) then
           jlo=0
           exit guess_start
        !decrease until the target is above
        else if((x < xx(jlo)) .eqv. ascnd)then
           jhi=jlo
           inc=inc+inc
        else
           exit guess_start
        endif
     end do guess_start
  endif

!3 continue
  binary_search: do
     !the end and the beginning are contiguous: segment number has been found
     if (jhi-jlo == 1) then
        !comment: this condition is known from the beginning, moving it
        !if(x == xx(n))jlo=n
        !if(x == xx(1))jlo=1
        exit binary_search
     endif
     !mean point (integer division, rounded towards jhi)
     jm=(jhi+jlo)/2
     !restrict search from the bottom of from the top
     if ((x >= xx(jm)) .eqv. ascnd) then
        jlo=jm
     else
        jhi=jm
     endif
  end do binary_search

  !print *,'good,jm:',jm,jlo,jhi

END SUBROUTINE hunt1








!!!!
!!!!
!!!!
!!!!
!!!!
!!!!
!!!!
!!!!
!!!!!> Wrapper of wpdot to avoid boundary problems in absence of wavelets
!!!!!! and to perform scalar product for complex wavefunctions and projectors
!!!!!! if the wavefunctions are complex, so should be also the projectors
!!!!subroutine wpdot_wrap_debug1(ncplx,mavctr_c,mavctr_f,maseg_c,maseg_f,keyav,keyag,apsi,  &
!!!!     mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keybv,keybg,bpsi,scpr)
!!!!  use module_base
!!!!  implicit none
!!!!  integer, intent(in) :: mavctr_c,mavctr_f,maseg_c,maseg_f
!!!!  integer, intent(in) :: ncplx,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f
!!!!  integer, dimension(maseg_c+maseg_f), intent(in) :: keyav
!!!!  integer, dimension(mbseg_c+mbseg_f), intent(in) :: keybv
!!!!  integer, dimension(2,maseg_c+maseg_f), intent(in) :: keyag
!!!!  integer, dimension(2,mbseg_c+mbseg_f), intent(in) :: keybg
!!!!  real(wp), dimension(mavctr_c+7*mavctr_f,ncplx), intent(in) :: apsi
!!!!  real(wp), dimension(mbvctr_c+7*mbvctr_f,ncplx), intent(in) :: bpsi
!!!!  real(dp), dimension(ncplx), intent(out) :: scpr
!!!!  !local variables
!!!!  integer :: ia_f,ib_f,iaseg_f,ibseg_f,ia,ib
!!!!  real(dp), dimension(ncplx,ncplx) :: scalprod 
!!!!
!!!!  ia_f=min(mavctr_f,1)
!!!!  ib_f=min(mbvctr_f,1)
!!!!
!!!!  iaseg_f=min(maseg_f,1)
!!!!  ibseg_f=min(mbseg_f,1)
!!!!
!!!!  do ia=1,ncplx
!!!!     do ib=1,ncplx
!!!!        call wpdot_debug1(mavctr_c,mavctr_f,maseg_c,maseg_f,&
!!!!             keyav,keyav(maseg_c+iaseg_f),&
!!!!             keyag,keyag(1,maseg_c+iaseg_f),&
!!!!             apsi(1,ia),apsi(mavctr_c+ia_f,ia),  &
!!!!             mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
!!!!             keybv,keybv(mbseg_c+ibseg_f),&
!!!!             keybg,keybg(1,mbseg_c+ibseg_f),&
!!!!             bpsi(1,ib),bpsi(mbvctr_c+ib_f,ib),scalprod(ia,ib))
!!!!     end do
!!!!  end do
!!!!
!!!!  !then define the result
!!!!  if (ncplx == 1) then
!!!!     scpr(1)=scalprod(1,1)
!!!!  else if (ncplx == 2) then
!!!!     scpr(1)=scalprod(1,1)+scalprod(2,2)
!!!!     scpr(2)=scalprod(1,2)-scalprod(2,1)
!!!!  else
!!!!     write(*,*)'ERROR wpdot: ncplx not valid:',ncplx
!!!!     stop
!!!!  end if
!!!!
!!!!END SUBROUTINE wpdot_wrap_debug1
!!!!
!!!!
!!!!!> This function must be generalized for the linear scaling code                
!!!!!! @warning   
!!!!!! calculates the dot product between a wavefunctions apsi and a projector bpsi (both in compressed form)
!!!!!! Warning: It is assumed that the segments of bpsi are always contained inside
!!!!!! the segments of apsi.
!!!!subroutine wpdot_debug1(  &
!!!!     mavctr_c,mavctr_f,maseg_c,maseg_f,keyav_c,keyav_f,keyag_c,keyag_f,apsi_c,apsi_f,  &
!!!!     mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keybv_c,keybv_f,keybg_c,keybg_f,bpsi_c,bpsi_f,scpr)
!!!!  use module_base
!!!!  implicit none
!!!!  integer, intent(in) :: mavctr_c,mavctr_f,maseg_c,maseg_f,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f
!!!!  integer, dimension(maseg_c), intent(in) :: keyav_c
!!!!  integer, dimension(mbseg_c), intent(in) :: keybv_c
!!!!  integer, dimension(maseg_f), intent(in) :: keyav_f
!!!!  integer, dimension(mbseg_f), intent(in) :: keybv_f
!!!!  integer, dimension(2,maseg_c), intent(in) :: keyag_c
!!!!  integer, dimension(2,mbseg_c), intent(in) :: keybg_c
!!!!  integer, dimension(2,maseg_f), intent(in) :: keyag_f
!!!!  integer, dimension(2,mbseg_f), intent(in) :: keybg_f
!!!!  real(wp), dimension(mavctr_c), intent(in) :: apsi_c
!!!!  real(wp), dimension(mbvctr_c), intent(in) :: bpsi_c
!!!!  real(wp), dimension(7,mavctr_f), intent(in) :: apsi_f
!!!!  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f
!!!!  real(dp), intent(out) :: scpr
!!!!  !local variables
!!!!  integer :: ibseg,jaj,jb1,jb0,jbj,iaoff,iboff,length,i,ja0,ja1
!!!!  real(dp) :: scpr1,scpr2,scpr3,scpr4,scpr5,scpr6,scpr7,scpr0
!!!!  integer :: iaseg0
!!!!  integer, dimension(maseg_c) :: keyag_c_lin !>linear version of second indices of keyag_c
!!!!  integer, dimension(maseg_f) :: keyag_f_lin !>linear version of second indices of keyag_f
!!!!!!!    integer :: ncount0,ncount2,ncount_rate,ncount_max
!!!!!!!    real(gp) :: tel
!!!!
!!!!
!!!!!!!  !dee
!!!!!!!    open(unit=97,file='time_wpdot',status='unknown',position='append')
!!!!!!!    call system_clock(ncount0,ncount_rate,ncount_max)
!!!!
!!!!  keyag_c_lin = keyag_c(1,:)!speed up access in hunt subroutine by consecutive arrangement in memory
!!!!  keyag_f_lin = keyag_f(1,:)!speed up access in hunt subroutine by consecutive arrangement in memory
!!!!
!!!!  scpr=0.0_dp
!!!!  scpr0=0.0_dp
!!!!  scpr1=0.0_dp
!!!!  scpr2=0.0_dp
!!!!  scpr3=0.0_dp
!!!!  scpr4=0.0_dp
!!!!  scpr5=0.0_dp
!!!!  scpr6=0.0_dp
!!!!  scpr7=0.0_dp
!!!!
!!!!!$omp parallel default (private) firstprivate(scpr0,scpr1,scpr2,scpr3,scpr4,scpr5,scpr6,scpr7)&
!!!!!$omp shared (maseg_c,keyav_c,keyag_c,keyag_c_lin,keybg_c,mbseg_c,keybv_c,mbseg_f,maseg_f)&
!!!!!$omp shared (apsi_c,bpsi_c,bpsi_f,keybv_f,keybg_f,keyag_f,keyag_f_lin,keyav_f)&
!!!!!$omp shared (apsi_f,scpr)
!!!!!!!!$omp shared (ncount0,ncount2,ncount_rate,ncount_max,tel)
!!!!
!!!!  iaseg0=1 
!!!!
!!!!!coarse part. Loop on the projectors segments
!!!!!$omp do schedule(static)
!!!!   do ibseg=1,mbseg_c
!!!!     jbj=keybv_c(ibseg)
!!!!!     jb0=keybg_c(1,ibseg) !starting point of projector segment
!!!!     jb0=max(keybg_c(1,ibseg),keyag_c_lin(1))
!!!!     jb1=keybg_c(2,ibseg) !ending point of projector segment
!!!!     iboff = max(jb0-keybg_c(1,ibseg),0)
!!!!     !print *,'huntenter',ibseg,jb0,jb1
!!!! 
!!!!     !find the starting point of the wavefunction segment
!!!!     !warning: hunt is assuming that the variable is always found
!!!!     !if it is not, iaseg0 is put to maseg + 1 so that the loop is disabled
!!!!!     call hunt1(.true.,keyag_c_lin,maseg_c,keybg_c(1,ibseg),iaseg0)
!!!!     call hunt1(.true.,keyag_c_lin,maseg_c,jb0,iaseg0)
!!!!     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
!!!!        iaseg0=1
!!!!        cycle     
!!!!     end if
!!!!     !now pass through all the wavefunction segments until the end of the segment is 
!!!!     !still contained in projector segment
!!!!     nonconvex_loop_c: do while(iaseg0 <= maseg_c)
!!!!!     print *,'huntexit',iaseg0,maseg_c,keyag_c_lin(iaseg0),keyag_c(2,iaseg0)
!!!!
!!!!        !length = jb1-jb0
!!!!        !iaoff = jb0-keyag_c_lin(iaseg0)!jb0-ja0
!!!!
!!!!        ja0=keyag_c_lin(iaseg0)
!!!!        ja1=min(jb1,keyag_c(2,iaseg0)) 
!!!!        length = ja1-jb0
!!!!        iaoff = max(jb0-ja0,0) !no offset if we are already inside
!!!!
!!!!        jaj=keyav_c(iaseg0)
!!!!        do i=0,length
!!!!           scpr0=scpr0+real(apsi_c(jaj+iaoff+i),dp)*real(bpsi_c(jbj+i+iboff),dp)
!!!!           write(1000,'(a,3i8,3es20.10)') 'i, jaj+iaoff+i, jbj+i+iboff, apsi_c(jaj+iaoff+i), bpsi_c(jbj+i+iboff), scpr0', i, jaj+iaoff+i, jbj+i+iboff, apsi_c(jaj+iaoff+i), bpsi_c(jbj+i+iboff), scpr0
!!!!        enddo
!!!!       !print *,'length',length,ibseg,scpr0,iaseg0,ja1,jb1
!!!!
!!!!        !print *,'ibseg,mbseg_c,iaseg0,maseg_c',ibseg,mbseg_c,iaseg0,maseg_c
!!!!        !print '(a,6(i8),1pe25.17)','ja0,ja1t,ja1,jb0,jb1',&
!!!!        !     ibseg,ja0,keyag_c(2,iaseg0),ja1,jb0,jb1,scpr0
!!!!        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_c) exit nonconvex_loop_c !segment is finished
!!!!        iaseg0=iaseg0+1
!!!!        jb0=max(keybg_c(1,ibseg),keyag_c_lin(iaseg0))
!!!!        if (keyag_c_lin(iaseg0)>jb1) exit nonconvex_loop_c !segment is not covered
!!!!        jbj=jbj+max(jb0-keybg_c(1,ibseg),0)
!!!!     end do nonconvex_loop_c
!!!!     !disable loop if the end is reached
!!!!     if (iaseg0 == maseg_c .and. keybg_c(1,ibseg)> keyag_c_lin(maseg_c)) iaseg0=iaseg0+1
!!!!
!!!!
!!!!   enddo
!!!!!stop
!!!!!$omp end do nowait
!!!!
!!!!! fine part
!!!!
!!!!iaseg0=1
!!!!
!!!!!$omp do schedule(static)
!!!!   do ibseg=1,mbseg_f
!!!!     jbj=keybv_f(ibseg)
!!!!     jb0=max(keybg_f(1,ibseg),keyag_f_lin(1))
!!!!     jb1=keybg_f(2,ibseg)
!!!!     iboff = max(jb0-keybg_f(1,ibseg),0)
!!!!!    print *,'huntenter',ibseg,jb0,jb1
!!!!     !call hunt1(.true.,keyag_f_lin,maseg_f,keybg_f(1,ibseg),iaseg0)
!!!!     call hunt1(.true.,keyag_f_lin,maseg_f,jb0,iaseg0)
!!!!     write(1011,'(a,7i8)') 'ibseg, jbj, jb0, jb1, iaseg0, keybg_f(1,ibseg), keyag_f_lin(14), ', ibseg, jbj, jb0, jb1, iaseg0, keybg_f(1,ibseg), keyag_f_lin(14)
!!!!     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
!!!!        iaseg0=1
!!!!        cycle     
!!!!     end if
!!!!     write(1010,'(a,6i8)') 'ibseg, jbj, jb0, jb1, iaseg0, maseg_f', ibseg, jbj, jb0, jb1, iaseg0, maseg_f
!!!!     nonconvex_loop_f: do while(iaseg0 <= maseg_f)
!!!!!!$     length = jb1-jb0
!!!!!!$     iaoff = jb0-keyag_f_lin(iaseg0)
!!!!
!!!!        ja0=keyag_f_lin(iaseg0) !still doubts about copying in automatic array
!!!!        ja1=min(jb1,keyag_f(2,iaseg0)) 
!!!!        length = ja1-jb0
!!!!        iaoff = max(jb0-ja0,0) !no offset if we are already inside
!!!!
!!!!        jaj=keyav_f(iaseg0)
!!!!        do i=0,length
!!!!           scpr1=scpr1+real(apsi_f(1,jaj+iaoff+i),dp)*real(bpsi_f(1,jbj+i+iboff),dp)
!!!!           scpr2=scpr2+real(apsi_f(2,jaj+iaoff+i),dp)*real(bpsi_f(2,jbj+i+iboff),dp)
!!!!           scpr3=scpr3+real(apsi_f(3,jaj+iaoff+i),dp)*real(bpsi_f(3,jbj+i+iboff),dp)
!!!!           scpr4=scpr4+real(apsi_f(4,jaj+iaoff+i),dp)*real(bpsi_f(4,jbj+i+iboff),dp)
!!!!           scpr5=scpr5+real(apsi_f(5,jaj+iaoff+i),dp)*real(bpsi_f(5,jbj+i+iboff),dp)
!!!!           scpr6=scpr6+real(apsi_f(6,jaj+iaoff+i),dp)*real(bpsi_f(6,jbj+i+iboff),dp)
!!!!           scpr7=scpr7+real(apsi_f(7,jaj+iaoff+i),dp)*real(bpsi_f(7,jbj+i+iboff),dp)
!!!!           write(1001,'(a,5i8,3es20.10)') 'ibseg, iaseg0, i, jaj+iaoff+i, jbj+i+iboff, apsi_f(1,jaj+iaoff+i), bpsi_f(1,jbj+i+iboff), scpr1', ibseg, iaseg0, i, jaj+iaoff+i, jbj+i+iboff, apsi_f(1,jaj+iaoff+i), bpsi_f(1,jbj+i+iboff), scpr1
!!!!           write(1002,'(a,5i8,3es20.10)') 'ibseg, iaseg0, i, jaj+iaoff+i, jbj+i+iboff, apsi_f(2,jaj+iaoff+i), bpsi_f(2,jbj+i+iboff), scpr2', ibseg, iaseg0, i, jaj+iaoff+i, jbj+i+iboff, apsi_f(2,jaj+iaoff+i), bpsi_f(2,jbj+i+iboff), scpr2
!!!!           write(1003,'(a,5i8,3es20.10)') 'ibseg, iaseg0, i, jaj+iaoff+i, jbj+i+iboff, apsi_f(3,jaj+iaoff+i), bpsi_f(3,jbj+i+iboff), scpr3', ibseg, iaseg0, i, jaj+iaoff+i, jbj+i+iboff, apsi_f(3,jaj+iaoff+i), bpsi_f(3,jbj+i+iboff), scpr3
!!!!           write(1004,'(a,5i8,3es20.10)') 'ibseg, iaseg0, i, jaj+iaoff+i, jbj+i+iboff, apsi_f(4,jaj+iaoff+i), bpsi_f(4,jbj+i+iboff), scpr4', ibseg, iaseg0, i, jaj+iaoff+i, jbj+i+iboff, apsi_f(4,jaj+iaoff+i), bpsi_f(4,jbj+i+iboff), scpr4
!!!!           write(1005,'(a,5i8,3es20.10)') 'ibseg, iaseg0, i, jaj+iaoff+i, jbj+i+iboff, apsi_f(5,jaj+iaoff+i), bpsi_f(5,jbj+i+iboff), scpr5', ibseg, iaseg0, i, jaj+iaoff+i, jbj+i+iboff, apsi_f(5,jaj+iaoff+i), bpsi_f(5,jbj+i+iboff), scpr5
!!!!           write(1006,'(a,5i8,3es20.10)') 'ibseg, iaseg0, i, jaj+iaoff+i, jbj+i+iboff, apsi_f(6,jaj+iaoff+i), bpsi_f(6,jbj+i+iboff), scpr6', ibseg, iaseg0, i, jaj+iaoff+i, jbj+i+iboff, apsi_f(6,jaj+iaoff+i), bpsi_f(6,jbj+i+iboff), scpr6
!!!!           write(1007,'(a,5i8,3es20.10)') 'ibseg, iaseg0, i, jaj+iaoff+i, jbj+i+iboff, apsi_f(7,jaj+iaoff+i), bpsi_f(7,jbj+i+iboff), scpr7', ibseg, iaseg0, i, jaj+iaoff+i, jbj+i+iboff, apsi_f(7,jaj+iaoff+i), bpsi_f(7,jbj+i+iboff), scpr7
!!!!        enddo
!!!! !       print *,'length',length,ibseg,scpr1,iaseg0,ja1,jb1
!!!!        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_f) exit nonconvex_loop_f !segment is finished  
!!!!        iaseg0=iaseg0+1
!!!!        jb0=max(keybg_f(1,ibseg),keyag_f_lin(iaseg0))
!!!!        if (keyag_f_lin(iaseg0)>jb1) exit nonconvex_loop_f !segment is not covered 
!!!!        jbj=jbj+max(jb0-keybg_f(1,ibseg),0)
!!!!     end do nonconvex_loop_f
!!!!     !disable loop if the end is reached
!!!!     if (iaseg0 == maseg_f .and. keybg_f(1,ibseg)> keyag_f_lin(maseg_f)) iaseg0=iaseg0+1
!!!!
!!!!   enddo
!!!!!$omp end do !!!implicit barrier 
!!!!
!!!!   scpr0=scpr0+scpr1+scpr2+scpr3+scpr4+scpr5+scpr6+scpr7
!!!!
!!!!!$omp critical 
!!!!   scpr=scpr+scpr0
!!!!!$omp end critical
!!!!
!!!!!$omp end parallel
!!!!
!!!!!!!    call system_clock(ncount2,ncount_rate,ncount_max)
!!!!!!!    tel=dble(ncount2-ncount0)/dble(ncount_rate)
!!!!!!!    write(97,*) 'wpdot:',tel
!!!!!!!    close(97)
!!!!
!!!!END SUBROUTINE wpdot_debug1
!!!!
!!!!
!!!!!> Wrapper of wpdot to avoid boundary problems in absence of wavelets
!!!!!! and to perform scalar product for complex wavefunctions and projectors
!!!!!! if the wavefunctions are complex, so should be also the projectors
!!!!subroutine wpdot_wrap_debug2(ncplx,mavctr_c,mavctr_f,maseg_c,maseg_f,keyav,keyag,apsi,  &
!!!!     mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keybv,keybg,bpsi,scpr)
!!!!  use module_base
!!!!  implicit none
!!!!  integer, intent(in) :: mavctr_c,mavctr_f,maseg_c,maseg_f
!!!!  integer, intent(in) :: ncplx,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f
!!!!  integer, dimension(maseg_c+maseg_f), intent(in) :: keyav
!!!!  integer, dimension(mbseg_c+mbseg_f), intent(in) :: keybv
!!!!  integer, dimension(2,maseg_c+maseg_f), intent(in) :: keyag
!!!!  integer, dimension(2,mbseg_c+mbseg_f), intent(in) :: keybg
!!!!  real(wp), dimension(mavctr_c+7*mavctr_f,ncplx), intent(in) :: apsi
!!!!  real(wp), dimension(mbvctr_c+7*mbvctr_f,ncplx), intent(in) :: bpsi
!!!!  real(dp), dimension(ncplx), intent(out) :: scpr
!!!!  !local variables
!!!!  integer :: ia_f,ib_f,iaseg_f,ibseg_f,ia,ib
!!!!  real(dp), dimension(ncplx,ncplx) :: scalprod 
!!!!
!!!!  ia_f=min(mavctr_f,1)
!!!!  ib_f=min(mbvctr_f,1)
!!!!
!!!!  iaseg_f=min(maseg_f,1)
!!!!  ibseg_f=min(mbseg_f,1)
!!!!
!!!!  do ia=1,ncplx
!!!!     do ib=1,ncplx
!!!!        call wpdot_debug2(mavctr_c,mavctr_f,maseg_c,maseg_f,&
!!!!             keyav,keyav(maseg_c+iaseg_f),&
!!!!             keyag,keyag(1,maseg_c+iaseg_f),&
!!!!             apsi(1,ia),apsi(mavctr_c+ia_f,ia),  &
!!!!             mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
!!!!             keybv,keybv(mbseg_c+ibseg_f),&
!!!!             keybg,keybg(1,mbseg_c+ibseg_f),&
!!!!             bpsi(1,ib),bpsi(mbvctr_c+ib_f,ib),scalprod(ia,ib))
!!!!     end do
!!!!  end do
!!!!
!!!!  !then define the result
!!!!  if (ncplx == 1) then
!!!!     scpr(1)=scalprod(1,1)
!!!!  else if (ncplx == 2) then
!!!!     scpr(1)=scalprod(1,1)+scalprod(2,2)
!!!!     scpr(2)=scalprod(1,2)-scalprod(2,1)
!!!!  else
!!!!     write(*,*)'ERROR wpdot: ncplx not valid:',ncplx
!!!!     stop
!!!!  end if
!!!!
!!!!END SUBROUTINE wpdot_wrap_debug2
!!!!
!!!!
!!!!!> This function must be generalized for the linear scaling code                
!!!!!! @warning   
!!!!!! calculates the dot product between a wavefunctions apsi and a projector bpsi (both in compressed form)
!!!!!! Warning: It is assumed that the segments of bpsi are always contained inside
!!!!!! the segments of apsi.
!!!!subroutine wpdot_debug2(  &
!!!!     mavctr_c,mavctr_f,maseg_c,maseg_f,keyav_c,keyav_f,keyag_c,keyag_f,apsi_c,apsi_f,  &
!!!!     mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keybv_c,keybv_f,keybg_c,keybg_f,bpsi_c,bpsi_f,scpr)
!!!!  use module_base
!!!!  implicit none
!!!!  integer, intent(in) :: mavctr_c,mavctr_f,maseg_c,maseg_f,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f
!!!!  integer, dimension(maseg_c), intent(in) :: keyav_c
!!!!  integer, dimension(mbseg_c), intent(in) :: keybv_c
!!!!  integer, dimension(maseg_f), intent(in) :: keyav_f
!!!!  integer, dimension(mbseg_f), intent(in) :: keybv_f
!!!!  integer, dimension(2,maseg_c), intent(in) :: keyag_c
!!!!  integer, dimension(2,mbseg_c), intent(in) :: keybg_c
!!!!  integer, dimension(2,maseg_f), intent(in) :: keyag_f
!!!!  integer, dimension(2,mbseg_f), intent(in) :: keybg_f
!!!!  real(wp), dimension(mavctr_c), intent(in) :: apsi_c
!!!!  real(wp), dimension(mbvctr_c), intent(in) :: bpsi_c
!!!!  real(wp), dimension(7,mavctr_f), intent(in) :: apsi_f
!!!!  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f
!!!!  real(dp), intent(out) :: scpr
!!!!  !local variables
!!!!  integer :: ibseg,jaj,jb1,jb0,jbj,iaoff,iboff,length,i,ja0,ja1
!!!!  real(dp) :: scpr1,scpr2,scpr3,scpr4,scpr5,scpr6,scpr7,scpr0
!!!!  integer :: iaseg0
!!!!  integer, dimension(maseg_c) :: keyag_c_lin !>linear version of second indices of keyag_c
!!!!  integer, dimension(maseg_f) :: keyag_f_lin !>linear version of second indices of keyag_f
!!!!!!!    integer :: ncount0,ncount2,ncount_rate,ncount_max
!!!!!!!    real(gp) :: tel
!!!!
!!!!
!!!!!!!  !dee
!!!!!!!    open(unit=97,file='time_wpdot',status='unknown',position='append')
!!!!!!!    call system_clock(ncount0,ncount_rate,ncount_max)
!!!!
!!!!  keyag_c_lin = keyag_c(1,:)!speed up access in hunt subroutine by consecutive arrangement in memory
!!!!  keyag_f_lin = keyag_f(1,:)!speed up access in hunt subroutine by consecutive arrangement in memory
!!!!
!!!!  scpr=0.0_dp
!!!!  scpr0=0.0_dp
!!!!  scpr1=0.0_dp
!!!!  scpr2=0.0_dp
!!!!  scpr3=0.0_dp
!!!!  scpr4=0.0_dp
!!!!  scpr5=0.0_dp
!!!!  scpr6=0.0_dp
!!!!  scpr7=0.0_dp
!!!!
!!!!!$omp parallel default (private) firstprivate(scpr0,scpr1,scpr2,scpr3,scpr4,scpr5,scpr6,scpr7)&
!!!!!$omp shared (maseg_c,keyav_c,keyag_c,keyag_c_lin,keybg_c,mbseg_c,keybv_c,mbseg_f,maseg_f)&
!!!!!$omp shared (apsi_c,bpsi_c,bpsi_f,keybv_f,keybg_f,keyag_f,keyag_f_lin,keyav_f)&
!!!!!$omp shared (apsi_f,scpr)
!!!!!!!!$omp shared (ncount0,ncount2,ncount_rate,ncount_max,tel)
!!!!
!!!!  iaseg0=1 
!!!!
!!!!!coarse part. Loop on the projectors segments
!!!!!$omp do schedule(static)
!!!!   do ibseg=1,mbseg_c
!!!!     jbj=keybv_c(ibseg)
!!!!!     jb0=keybg_c(1,ibseg) !starting point of projector segment
!!!!     jb0=max(keybg_c(1,ibseg),keyag_c_lin(1))
!!!!     jb1=keybg_c(2,ibseg) !ending point of projector segment
!!!!     iboff = max(jb0-keybg_c(1,ibseg),0)
!!!!     !print *,'huntenter',ibseg,jb0,jb1
!!!! 
!!!!     !find the starting point of the wavefunction segment
!!!!     !warning: hunt is assuming that the variable is always found
!!!!     !if it is not, iaseg0 is put to maseg + 1 so that the loop is disabled
!!!!!     call hunt1(.true.,keyag_c_lin,maseg_c,keybg_c(1,ibseg),iaseg0)
!!!!     call hunt1(.true.,keyag_c_lin,maseg_c,jb0,iaseg0)
!!!!     write(2011,'(a,7i8)') 'ibseg, jbj, jb0, jb1, iaseg0, keybg_f(1,ibseg), keyag_f_lin(14), ', ibseg, jbj, jb0, jb1, iaseg0, keybg_f(1,ibseg), keyag_f_lin(14)
!!!!     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
!!!!        iaseg0=1
!!!!        cycle     
!!!!     end if
!!!!     !now pass through all the wavefunction segments until the end of the segment is 
!!!!     !still contained in projector segment
!!!!     nonconvex_loop_c: do while(iaseg0 <= maseg_c)
!!!!!     print *,'huntexit',iaseg0,maseg_c,keyag_c_lin(iaseg0),keyag_c(2,iaseg0)
!!!!
!!!!        !length = jb1-jb0
!!!!        !iaoff = jb0-keyag_c_lin(iaseg0)!jb0-ja0
!!!!
!!!!        ja0=keyag_c_lin(iaseg0)
!!!!        ja1=min(jb1,keyag_c(2,iaseg0)) 
!!!!        length = ja1-jb0
!!!!        iaoff = max(jb0-ja0,0) !no offset if we are already inside
!!!!
!!!!        jaj=keyav_c(iaseg0)
!!!!        do i=0,length
!!!!           scpr0=scpr0+real(apsi_c(jaj+iaoff+i),dp)*real(bpsi_c(jbj+i+iboff),dp)
!!!!           write(2000,'(a,3i8,3es20.10)') 'i, jaj+iaoff+i, jbj+i+iboff, apsi_c(jaj+iaoff+i), bpsi_c(jbj+i+iboff), scpr0', i, jaj+iaoff+i, jbj+i+iboff, apsi_c(jaj+iaoff+i), bpsi_c(jbj+i+iboff), scpr0
!!!!        enddo
!!!!       !print *,'length',length,ibseg,scpr0,iaseg0,ja1,jb1
!!!!
!!!!        !print *,'ibseg,mbseg_c,iaseg0,maseg_c',ibseg,mbseg_c,iaseg0,maseg_c
!!!!        !print '(a,6(i8),1pe25.17)','ja0,ja1t,ja1,jb0,jb1',&
!!!!        !     ibseg,ja0,keyag_c(2,iaseg0),ja1,jb0,jb1,scpr0
!!!!        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_c) exit nonconvex_loop_c !segment is finished
!!!!        iaseg0=iaseg0+1
!!!!        jb0=max(keybg_c(1,ibseg),keyag_c_lin(iaseg0))
!!!!        if (keyag_c_lin(iaseg0)>jb1) exit nonconvex_loop_c !segment is not covered
!!!!        jbj=jbj+max(jb0-keybg_c(1,ibseg),0)
!!!!     end do nonconvex_loop_c
!!!!     !disable loop if the end is reached
!!!!     if (iaseg0 == maseg_c .and. keybg_c(1,ibseg)> keyag_c_lin(maseg_c)) iaseg0=iaseg0+1
!!!!
!!!!
!!!!   enddo
!!!!!stop
!!!!!$omp end do nowait
!!!!
!!!!! fine part
!!!!
!!!!iaseg0=1
!!!!
!!!!!$omp do schedule(static)
!!!!   do ibseg=1,mbseg_f
!!!!     jbj=keybv_f(ibseg)
!!!!     jb0=max(keybg_f(1,ibseg),keyag_f_lin(1))
!!!!     jb1=keybg_f(2,ibseg)
!!!!     iboff = max(jb0-keybg_f(1,ibseg),0)
!!!!!    print *,'huntenter',ibseg,jb0,jb1
!!!!     !call hunt1(.true.,keyag_f_lin,maseg_f,keybg_f(1,ibseg),iaseg0)
!!!!     call hunt1(.true.,keyag_f_lin,maseg_f,jb0,iaseg0)
!!!!     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
!!!!        iaseg0=1
!!!!        cycle     
!!!!     end if
!!!!     write(2010,'(a,6i8)') 'ibseg, jbj, jb0, jb1, iaseg0, maseg_f', ibseg, jbj, jb0, jb1, iaseg0, maseg_f
!!!!     nonconvex_loop_f: do while(iaseg0 <= maseg_f)
!!!!!!$     length = jb1-jb0
!!!!!!$     iaoff = jb0-keyag_f_lin(iaseg0)
!!!!
!!!!        ja0=keyag_f_lin(iaseg0) !still doubts about copying in automatic array
!!!!        ja1=min(jb1,keyag_f(2,iaseg0)) 
!!!!        length = ja1-jb0
!!!!        iaoff = max(jb0-ja0,0) !no offset if we are already inside
!!!!
!!!!        jaj=keyav_f(iaseg0)
!!!!        do i=0,length
!!!!           scpr1=scpr1+real(apsi_f(1,jaj+iaoff+i),dp)*real(bpsi_f(1,jbj+i+iboff),dp)
!!!!           scpr2=scpr2+real(apsi_f(2,jaj+iaoff+i),dp)*real(bpsi_f(2,jbj+i+iboff),dp)
!!!!           scpr3=scpr3+real(apsi_f(3,jaj+iaoff+i),dp)*real(bpsi_f(3,jbj+i+iboff),dp)
!!!!           scpr4=scpr4+real(apsi_f(4,jaj+iaoff+i),dp)*real(bpsi_f(4,jbj+i+iboff),dp)
!!!!           scpr5=scpr5+real(apsi_f(5,jaj+iaoff+i),dp)*real(bpsi_f(5,jbj+i+iboff),dp)
!!!!           scpr6=scpr6+real(apsi_f(6,jaj+iaoff+i),dp)*real(bpsi_f(6,jbj+i+iboff),dp)
!!!!           scpr7=scpr7+real(apsi_f(7,jaj+iaoff+i),dp)*real(bpsi_f(7,jbj+i+iboff),dp)
!!!!           write(2001,'(a,5i8,3es20.10)') 'ibseg, iaseg0, i, jaj+iaoff+i, jbj+i+iboff, apsi_f(1,jaj+iaoff+i), bpsi_f(1,jbj+i+iboff), scpr1', ibseg, iaseg0, i, jaj+iaoff+i, jbj+i+iboff, apsi_f(1,jaj+iaoff+i), bpsi_f(1,jbj+i+iboff), scpr1
!!!!           write(2002,'(a,5i8,3es20.10)') 'ibseg, iaseg0, i, jaj+iaoff+i, jbj+i+iboff, apsi_f(2,jaj+iaoff+i), bpsi_f(2,jbj+i+iboff), scpr2', ibseg, iaseg0, i, jaj+iaoff+i, jbj+i+iboff, apsi_f(2,jaj+iaoff+i), bpsi_f(2,jbj+i+iboff), scpr2
!!!!           write(2003,'(a,5i8,3es20.10)') 'ibseg, iaseg0, i, jaj+iaoff+i, jbj+i+iboff, apsi_f(3,jaj+iaoff+i), bpsi_f(3,jbj+i+iboff), scpr3', ibseg, iaseg0, i, jaj+iaoff+i, jbj+i+iboff, apsi_f(3,jaj+iaoff+i), bpsi_f(3,jbj+i+iboff), scpr3
!!!!           write(2004,'(a,5i8,3es20.10)') 'ibseg, iaseg0, i, jaj+iaoff+i, jbj+i+iboff, apsi_f(4,jaj+iaoff+i), bpsi_f(4,jbj+i+iboff), scpr4', ibseg, iaseg0, i, jaj+iaoff+i, jbj+i+iboff, apsi_f(4,jaj+iaoff+i), bpsi_f(4,jbj+i+iboff), scpr4
!!!!           write(2005,'(a,5i8,3es20.10)') 'ibseg, iaseg0, i, jaj+iaoff+i, jbj+i+iboff, apsi_f(5,jaj+iaoff+i), bpsi_f(5,jbj+i+iboff), scpr5', ibseg, iaseg0, i, jaj+iaoff+i, jbj+i+iboff, apsi_f(5,jaj+iaoff+i), bpsi_f(5,jbj+i+iboff), scpr5
!!!!           write(2006,'(a,5i8,3es20.10)') 'ibseg, iaseg0, i, jaj+iaoff+i, jbj+i+iboff, apsi_f(6,jaj+iaoff+i), bpsi_f(6,jbj+i+iboff), scpr6', ibseg, iaseg0, i, jaj+iaoff+i, jbj+i+iboff, apsi_f(6,jaj+iaoff+i), bpsi_f(6,jbj+i+iboff), scpr6
!!!!           write(2007,'(a,5i8,3es20.10)') 'ibseg, iaseg0, i, jaj+iaoff+i, jbj+i+iboff, apsi_f(7,jaj+iaoff+i), bpsi_f(7,jbj+i+iboff), scpr7', ibseg, iaseg0, i, jaj+iaoff+i, jbj+i+iboff, apsi_f(7,jaj+iaoff+i), bpsi_f(7,jbj+i+iboff), scpr7
!!!!        enddo
!!!! !       print *,'length',length,ibseg,scpr1,iaseg0,ja1,jb1
!!!!        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_f) exit nonconvex_loop_f !segment is finished  
!!!!        iaseg0=iaseg0+1
!!!!        jb0=max(keybg_f(1,ibseg),keyag_f_lin(iaseg0))
!!!!        if (keyag_f_lin(iaseg0)>jb1) exit nonconvex_loop_f !segment is not covered 
!!!!        jbj=jbj+max(jb0-keybg_f(1,ibseg),0)
!!!!     end do nonconvex_loop_f
!!!!     !disable loop if the end is reached
!!!!     if (iaseg0 == maseg_f .and. keybg_f(1,ibseg)> keyag_f_lin(maseg_f)) iaseg0=iaseg0+1
!!!!
!!!!   enddo
!!!!!$omp end do !!!implicit barrier 
!!!!
!!!!   scpr0=scpr0+scpr1+scpr2+scpr3+scpr4+scpr5+scpr6+scpr7
!!!!
!!!!!$omp critical 
!!!!   scpr=scpr+scpr0
!!!!!$omp end critical
!!!!
!!!!!$omp end parallel
!!!!
!!!!!!!    call system_clock(ncount2,ncount_rate,ncount_max)
!!!!!!!    tel=dble(ncount2-ncount0)/dble(ncount_rate)
!!!!!!!    write(97,*) 'wpdot:',tel
!!!!!!!    close(97)
!!!!
!!!!END SUBROUTINE wpdot_debug2
