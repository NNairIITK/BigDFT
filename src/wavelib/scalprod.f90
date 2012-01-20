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
  real(wp), dimension(mvctr_c+7*mvctr_f), intent(in) :: psi
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
  integer :: ibseg,jaj,jb1,jb0,jbj,iaoff,length,i,ja0,ja1
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
           scpr0=scpr0+real(apsi_c(jaj+iaoff+i),dp)*real(bpsi_c(jbj+i),dp)
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
           scpr1=scpr1+real(apsi_f(1,jaj+iaoff+i),dp)*real(bpsi_f(1,jbj+i),dp)
           scpr2=scpr2+real(apsi_f(2,jaj+iaoff+i),dp)*real(bpsi_f(2,jbj+i),dp)
           scpr3=scpr3+real(apsi_f(3,jaj+iaoff+i),dp)*real(bpsi_f(3,jbj+i),dp)
           scpr4=scpr4+real(apsi_f(4,jaj+iaoff+i),dp)*real(bpsi_f(4,jbj+i),dp)
           scpr5=scpr5+real(apsi_f(5,jaj+iaoff+i),dp)*real(bpsi_f(5,jbj+i),dp)
           scpr6=scpr6+real(apsi_f(6,jaj+iaoff+i),dp)*real(bpsi_f(6,jbj+i),dp)
           scpr7=scpr7+real(apsi_f(7,jaj+iaoff+i),dp)*real(bpsi_f(7,jbj+i),dp)
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
