!! @file
!! @author Bastian Schaefer
!! @section LICENCE
!!    Copyright (C) 2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

module module_hessian
   implicit none

   private

   public :: cal_hessian_fd

contains

   !> Reza's routine for finite difference hessian
   subroutine cal_hessian_fd(iproc,nat,alat,pos,hess)
      use module_base, only: gp, f_malloc, f_free, assignment(=)
      use module_energyandforces, only: mhgpsenergyandforces
      implicit none
      integer, intent(in):: iproc, nat
      real(gp), dimension(3*nat), intent(in) :: pos
      real(gp), dimension(3), intent(in) :: alat
      real(gp), dimension(3*nat,3*nat), intent(inout) :: hess
      !local variables
      integer :: iat
      !real(gp) :: t1,t2,t3
      !real(gp), allocatable, dimension(:,:) :: hess
      real(gp), allocatable, dimension(:) :: tpos,grad,eval,workf
      real(gp) :: h,rlarge,twelfth,twothird,etot,cmx,cmy,cmz,dm,tt
      real(gp) :: s,fnoise
      integer :: i,j,k,lworkf

      !allocate(hess(3*nat,3*nat))
      tpos = f_malloc(3*nat,id='tpos')
      grad = f_malloc(3*nat,id='grad')
      eval = f_malloc(3*nat,id='eval')

      lworkf=1000*nat
      workf = f_malloc(lworkf,id='workf')

      !h=1.e-1_gp
      !h=7.5e-2_gp
      !h=5.e-2_gp
         h=1.e-3_gp
     ! h=1.e-2_gp
      !h=5.e-3_gp
      !h=2.e-2_gp
      rlarge=1._gp*1.e4_gp
      twelfth=-1._gp/(12._gp*h)
      twothird=-2._gp/(3._gp*h)
      if(iproc==0) write(*,*) '(hess) HESSIAN: h',h
      !-------------------------------------------------------
      do i=1,3*nat
         iat=(i-1)/3+1
         do k=1,3*nat
             tpos(k)=pos(k)
             grad(k)=0._gp
         enddo
         !-----------------------------------------
         tpos(i)=tpos(i)-2*h
         call mhgpsenergyandforces(nat,alat,tpos,grad,fnoise,etot)
         do j=1,3*nat
             hess(j,i)=twelfth*grad(j)
         enddo
         !if(iproc==0) write(*,*) 'ALIREZA-6',i,iat
         !-----------------------------------------
         tpos(i)=tpos(i)+h
         call mhgpsenergyandforces(nat,alat,tpos,grad,fnoise,etot)
         do j=1,3*nat
         hess(j,i)=hess(j,i)-twothird*grad(j)
         enddo
         !-----------------------------------------
         tpos(i)=tpos(i)+2*h
         call mhgpsenergyandforces(nat,alat,tpos,grad,fnoise,etot)
         do j=1,3*nat
         hess(j,i)=hess(j,i)+twothird*grad(j)
         enddo
         !-----------------------------------------
         tpos(i)=tpos(i)+h
         call mhgpsenergyandforces(nat,alat,tpos,grad,fnoise,etot)
         do j=1,3*nat
         hess(j,i)=hess(j,i)-twelfth*grad(j)
         !write(*,*) 'HESS ',j,i,hess(j,i)
         enddo
         !-----------------------------------------
      enddo
      !-------------------------------------------------------

      !check symmetry
      dm=0._gp
      do i=1,3*nat
         do j=1,i-1
            s=.5_gp*(hess(i,j)+hess(j,i))
            tt=abs(hess(i,j)-hess(j,i))/(1._gp+abs(s))
            dm=max(dm,tt)
            hess(i,j)=s
            hess(j,i)=s
         enddo
      enddo
      if (dm.gt.1.e-1_gp) write(*,*) '(hess) max dev from sym',dm

   !   do j=1,3*nat
   !   do i=1,3*nat
   !   write(*,*) '(hess) hier',nat,hess(i,j)
   !   write(499,*) hess(i,j)
   !   enddo
   !   enddo

      !-------------------------------------------------------
      !project out rotations
      cmx=0._gp ; cmy=0._gp ; cmz=0._gp
      do i=1,3*nat-2,3
         cmx=cmx+pos(i+0)
         cmy=cmy+pos(i+1)
         cmz=cmz+pos(i+2)
      enddo
      cmx=cmx/nat ; cmy=cmy/nat ; cmz=cmz/nat
     
      !x-y plane
      do i=1,3*nat-2,3
         workf(i+1)= (pos(i+0)-cmx)
         workf(i+0)=-(pos(i+1)-cmy)
      enddo
      call f_free(tpos)
      call f_free(grad)
      call f_free(eval)
      call f_free(workf)
   end subroutine cal_hessian_fd

end module module_hessian
