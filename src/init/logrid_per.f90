!> @file
!! Routines to creation localisation regions
!! @author
!!    Copyright (C) 2010-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Define the bounds of wavefunctions for periodic systems
subroutine make_bounds_per(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,bounds,wfd)
  use module_base
  use locregs
  implicit none
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(convolutions_bounds),intent(out):: bounds
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3

  logical,allocatable,dimension(:,:,:) :: logrid
  character(len=*), parameter :: subname='make_bounds_per'
  integer :: i_stat,i_all,nseg_c

  bounds%kb%ibyz_f = f_malloc_ptr((/ 1.to.2, 0.to.n2, 0.to.n3 /),id='bounds%kb%ibyz_f')
  bounds%kb%ibxz_f = f_malloc_ptr((/ 1.to.2, 0.to.n1, 0.to.n3 /),id='bounds%kb%ibxz_f')
  bounds%kb%ibxy_f = f_malloc_ptr((/ 1.to.2, 0.to.n1, 0.to.n2 /),id='bounds%kb%ibxy_f')

  bounds%gb%ibyz_ff = f_malloc_ptr((/ 1.to.2, nfl2.to.nfu2, nfl3.to.nfu3 /),id='bounds%gb%ibyz_ff')
  bounds%gb%ibzxx_f = f_malloc_ptr((/ 1.to.2, nfl3.to.nfu3, 0.to.2*n1+1 /),id='bounds%gb%ibzxx_f')
  bounds%gb%ibxxyy_f = f_malloc_ptr((/ 1.to.2, 0.to.2*n1+1, 0.to.2*n2+1 /),id='bounds%gb%ibxxyy_f')

  bounds%sb%ibxy_ff = f_malloc_ptr((/ 1.to.2, nfl1.to.nfu1, nfl2.to.nfu2 /),id='bounds%sb%ibxy_ff')
  bounds%sb%ibzzx_f = f_malloc_ptr((/ 1.to.2, 0.to.2*n3+1, nfl1.to.nfu1 /),id='bounds%sb%ibzzx_f')
  bounds%sb%ibyyzz_f = f_malloc_ptr((/ 1.to.2, 0.to.2*n2+1, 0.to.2*n3+1 /),id='bounds%sb%ibyyzz_f')

  logrid = f_malloc((/ 0.to.n1, 0.to.n2, 0.to.n3 /),id='logrid')

  nseg_c=wfd%nseg_c
  call make_logrid_f(n1,n2,n3, & 
       wfd%nseg_f,wfd%keygloc(1,nseg_c+min(1,wfd%nseg_f)),  & 
       logrid)

  call make_bounds(n1,n2,n3,logrid,bounds%kb%ibyz_f,bounds%kb%ibxz_f,bounds%kb%ibxy_f)

  call f_free(logrid)

END SUBROUTINE make_bounds_per


subroutine make_all_ib_per(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
     ibxy_f,ibxy_ff,ibzzx_f,ibyyzz_f,&
     ibyz_f,ibyz_ff,ibzxx_f,ibxxyy_f)
  !    creates complicated ib arrays    
  use module_base
  implicit none
  integer,intent(in)::n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  integer :: i1,i2,i3,i_stat,i_all !n(c) m1,m2,m3

  integer,intent(in):: ibyz_f(2,0:n2,0:n3+ndebug),ibxy_f(2,0:n1,0:n2+ndebug)

  !    for shrink:    
  integer,intent(out):: ibxy_ff(2,nfl1:nfu1,nfl2:nfu2)
  integer,intent(out):: ibzzx_f(2,0:2*n3+1,nfl1:nfu1) 
  integer,intent(out):: ibyyzz_f(2,0:2*n2+1,0:2*n3+1)

  !    for grow:    
  integer,intent(out):: ibyz_ff(2,nfl2:nfu2,nfl3:nfu3)
  integer,intent(out):: ibzxx_f(2,nfl3:nfu3,0:2*n1+1)
  integer,intent(out):: ibxxyy_f(2,0:2*n1+1,0:2*n2+1)

  character(len=*), parameter :: subname=' make_all_ib'
  logical,allocatable:: logrid_big(:)

  logrid_big = f_malloc((2*n1+2)*(2*n2+2)*(2*n3+2),id='logrid_big')

  !n(c) m1=nfu1-nfl1
  !n(c) m2=nfu2-nfl2
  !n(c) m3=nfu3-nfl3

  !   (0:n3,0:2*n1+1,0:2*n2+1) from grow
  !   (0:2*n2+1,0:2*n3+1,0:n1) from shrink

  !   for shrink:
  do i2=nfl2,nfu2
     do i1=nfl1,nfu1
        ibxy_ff(:,i1,i2)=ibxy_f(:,i1,i2)
     enddo
  enddo

  call make_ib_inv_per(logrid_big,ibxy_ff,ibzzx_f,ibyyzz_f,&
             n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3)

  !    for grow:
  do i2=nfl2,nfu2
     do i3=nfl3,nfu3
        ibyz_ff(:,i2,i3)=ibyz_f(:,i2,i3)
!      write(9,*) i2,i3,ibyz_ff(1,i2,i3),ibyz_ff(2,i2,i3)
     enddo
  enddo

  call make_ib_per(logrid_big,ibyz_ff,ibzxx_f,ibxxyy_f,n1,n2,&
       nfl2,nfu2,nfl3,nfu3)

  call f_free(logrid_big)

END SUBROUTINE make_all_ib_per


!>   This subroutine mimics the comb_grow_f one
subroutine make_ib_inv_per(logrid_big,ibxy,ibzzx,ibyyzz,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3)
  implicit none
  integer nt,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,n1,n2,n3
  integer,intent(out):: ibxy(2,nfl1:nfu1,nfl2:nfu2)
  integer,intent(out):: ibzzx(2,0:2*nfu3+1,nfl1:nfu1) 
  integer,intent(out):: ibyyzz(2,0:2*n2+1,0:2*n3+1)
  logical logrid_big(nfl3:nfu3,0:2*n1+1,0:2*n2+1)! work array

  ! I3,i1,i2 -> i1,i2,i3 
  nt=(nfu1-nfl1+1)*(nfu2-nfl2+1)
  call ib_to_logrid_inv_per(ibxy,logrid_big,n3,nt)

  ! I2,I3,i1 -> I3,i1,i2
  nt=(2*n3+2)*(nfu1-nfl1+1)
  call ib_from_logrid_inv(ibzzx,logrid_big,nfl2,nfu2,nt)
  call ib_to_logrid_inv_per(ibzzx,logrid_big,n2,nt)

  ! I1,I2,I3  -> I2,I3,i1
  nt=(2*n2+2)*(2*n3+2)
  call ib_from_logrid_inv( ibyyzz,logrid_big,nfl1,nfu1,nt)

END SUBROUTINE make_ib_inv_per


!> This one mimics the comb_rot_grow_f_loc
subroutine ib_to_logrid_inv_per(ib,logrid,n,ndat)
  implicit none
  integer ndat,l,i,n,ii
  integer ib(2,ndat)! input
  logical logrid(0:2*n+1,ndat)! output

  logrid=.false.

  do l=1,ndat
     do i = 2*ib(1,l)-14 , 2*ib(2,l)+16
      ii=modulo(i,2*n+2) 
        logrid(ii,l)=.true.
     enddo
  enddo

END SUBROUTINE ib_to_logrid_inv_per


!>    This subroutine mimics the comb_grow_f one
subroutine make_ib_per(logrid_big,ibyz,ibzxx,ibxxyy,n1,n2,nfl2,nfu2,nfl3,nfu3)
  implicit none
  integer :: nt,nfl2,nfu2,nfl3,nfu3,n1,n2
  integer :: ibyz(  2,nfl2:nfu2,nfl3:nfu3)! input
  integer :: ibzxx( 2,          nfl3:nfu3,0:2*n1+1)!output
  integer :: ibxxyy(2,                    0:2*n1+1,0:2*n2+1)!output
  logical :: logrid_big(        nfl3:nfu3,0:2*n1+1,0:2*n2+1)! work array

  ! i1,i2,i3 -> i2,i3,I1
  nt=(nfu2-nfl2+1)*(nfu3-nfl3+1)
  call ib_to_logrid_rot_per(  ibyz,logrid_big,n1,nt)

  ! i2,i3,I1 -> i3,I1,I2
  nt=(nfu3-nfl3+1)*(2*n1+2)
  call ib_from_logrid( ibzxx,logrid_big,nfl2,nfu2,nt)
  call ib_to_logrid_rot_per( ibzxx,logrid_big,n2,nt)

  ! i3,I1,I2  -> I1,I2,I3
  nt=(2*n1+2)*(2*n2+2)
  call ib_from_logrid(ibxxyy,logrid_big,nfl3,nfu3,nt)

END SUBROUTINE make_ib_per


!>   This one mimics the comb_rot_grow_f_loc
subroutine ib_to_logrid_rot_per(ib,logrid,n,ndat)
  implicit none
  integer ndat,n,l,i,ii
  integer ib(2,ndat)! input
  logical logrid(ndat,0:2*n+1)! output

  logrid=.false.

  do l=1,ndat
     do i = 2*ib(1,l)-14 , 2*ib(2,l)+16
        ii=modulo(i,2*n+2) 
        logrid(l,ii)=.true.
     enddo
  enddo

END SUBROUTINE ib_to_logrid_rot_per


subroutine make_logrid_f(n1,n2,n3, & 
     mseg_f,keyg_f,&
     logrid)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  integer, intent(in) :: mseg_f
  integer, dimension(2,mseg_f), intent(in) :: keyg_f
  logical,intent(out),dimension(0:n1,0:n2,0:n3)::logrid
  !local variables
  integer :: iseg,j0,j1,ii,i1,i2,i3,i0,i !n(c) jj

  logrid=.false.
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
        logrid(i,i2,i3)=.true. 
     enddo
  enddo

END SUBROUTINE make_logrid_f
