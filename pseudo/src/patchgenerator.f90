!> @file
!! Paw generation in pseudo program
!! @author
!!    Alex Willand, under the supervision of Stefan Goedecker
!!    gpu accelerated routines by Raffael Widmer
!!    parts of this program were based on the fitting program by matthias krack
!!    http://cvs.berlios.de/cgi-bin/viewcvs.cgi/cp2k/potentials/goedecker/pseudo/v2.2/
!!
!!    Copyright (C) 2010-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


subroutine paw_generator(zion, lmx,  lpmx, lmax,  hsep, gpot, &
     alpz, alps, &
     ng, noccmax,noccmx, expo,&
     psi, aeval, occup, &
     Nsol, Labs, Ngrid,Ngrid_box, Egrid,  rgrid ,rw,rd,  psigrid, Npaw,&
     PAWpatch , psipsigrid, rprb, rcore, zcore , Ngrid_box_larger)

  implicit none
  !Arguments
  integer, intent(in) :: ng,noccmax,noccmx,Nsol, labs, Ngrid,  Ngrid_box, Ngrid_box_larger
  integer, intent(in)::lpmx, lmx, lmax
  integer, intent(in) :: Npaw
  real(kind=8), dimension(ng+1), intent(out) :: expo
  real(kind=8), intent(in), dimension(Ngrid) :: rgrid, rd, rw
  real(kind=8), intent(in) :: rcore, zcore
  real(kind=8), dimension(0:ng,noccmax, lmx), intent(out) :: psi
  real(kind=8), dimension(Nsol), intent(out) :: Egrid
  real(kind=8), dimension(Ngrid,Nsol) , intent(out) :: psigrid, psipsigrid
  real(kind=8), dimension(4), intent(in) :: gpot !! gpot dim e diventata 4!!!
  real(kind=8), intent(in) :: rprb, zion
  real(kind=8), dimension(noccmx,lmx), intent(in) ::  occup
  real(kind=8), dimension(noccmx,lmx), intent(out) ::  aeval
  real(kind=8), dimension(Npaw,Npaw), intent(out):: PAWpatch
  real(kind=8), dimension(6,lpmx), intent(in):: hsep
  real(kind=8), intent(in) :: alpz
  real(kind=8), dimension(*), intent(in) :: alps
 
  !Local variables
  integer, parameter :: n_int=1000
  real(kind=8), parameter :: fact=4.0d0
  real(kind=8) :: alpl
  !real(kind=8), dimension(noccmx,lmx) ::chrg,res
  real(kind=8), dimension(:), allocatable :: xp
  real(kind=8), dimension(:,:), allocatable :: vh

  real(kind=8), dimension(:,:,:,:), allocatable :: rmt
  integer :: l,i,iocc 
  real(kind=8) :: rij,a,a0,a0in,tt
  real(kind=8) :: value_at_r
  integer :: lpx

  !filename = 'psppar.'//trim(atomname)

  lpx=0

  lpx_determination: do i=1,4
     if (alps(i) == 0.0_8) then
        exit lpx_determination
     else
        lpx=i-1
     end if
  end do lpx_determination

  alpl=alpz
  
  !allocate arrays for the gatom routine
  allocate(vh(4*(ng+1)**2,4*(ng+1)**2))
  
  allocate(xp(0:ng))
  allocate(rmt(n_int,0:ng,0:ng,lmax+1))
  
  !can be switched on for debugging
  !if (iproc.eq.0) write(*,'(1x,a,a7,a9,i3,i3,a9,i3,f5.2)')&
  !     'Input Guess Generation for atom',trim(atomname),&
  !     'Z,Zion=',izatom,ielpsp,'ng,rprb=',ng+1,rprb
  
  rij=3._8
  ! exponents of gaussians
  ! print *, " ESPONENTI " 
  ! print *, " alpz " , alpz
  a0in=alpz
  a0=a0in/rij
  !       tt=sqrt(sqrt(2._8))
  tt=2._8**(.3_8)
  do i=0,ng
     a=a0*tt**i
     xp(i)=.5_8/a**2
     ! print *, " xp(", i,")", xp(i)
  end do
  
  ! initial guess
  do l=0,lmx-1
     do iocc=1,noccmax
        do i=0,ng
           psi(i,iocc,l+1)=0.0_8
        end do
     end do
  end do
  
  call crtvh_paw(ng,lmax,xp,vh,rprb,fact,n_int,rmt)
  
!!!  call gatom(rcov,rprb,lmax,lpx,noccmax,occup,&
!!!       zion,alpz,gpot,alpl,hsep,alps,vh,xp,rmt,fact,n_int,&
!!!       aeval,ng,psi,res,chrg)
  
  PAWpatch=0.0_8
  call gatom_modified(rprb,lmax,lpx,lpmx, noccmax,noccmx, occup,&
       zion,alpz,gpot,alpl,hsep,alps,vh,xp,rmt,fact,n_int,&
       aeval,ng,psi,&
       Nsol, Labs, Ngrid,Ngrid_box,Egrid,  rgrid,rw, rd,  psigrid,Npaw,  PAWpatch,&
       psipsigrid,rcore,zcore , Ngrid_box_larger   )              
 
  
  do i=1,ng+1
     expo(i)=sqrt(0.5_8/xp(i-1))
  end do
  
  do l=0,lmx-1
     do iocc=1,noccmax
        if( value_at_r(rprb, ng , expo,psi(0,iocc,l+1)).lt.0.0     ) then
           do i=0,ng
              psi(i,iocc,l+1)=-psi(i,iocc,l+1)
           enddo
        endif
     enddo
  enddo
  
END SUBROUTINE paw_generator


function value_at_r(r, ng , expo,psi     )

  implicit none
  integer, parameter :: gp=kind(1.0d0) 
 
  real(gp) , intent(in) ::  r
  integer, intent(in) :: ng
  real(gp), dimension(ng+1), intent(in) :: expo
  real(gp), dimension(0:ng), intent(in) :: psi

  ! local
  integer, parameter :: n_int=1000

  integer ig
  real(gp) sum
  real(gp) :: value_at_r

  sum=0.0
  do ig = 0,ng
     sum=sum+psi(ig)*exp( -r*r/2.0/expo(ig+1)/expo(ig+1) )
  enddo
  value_at_r=sum

end function value_at_r
 
