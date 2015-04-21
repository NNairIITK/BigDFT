!! @file
!! @author Bastian Schaefer
!! @section LICENCE
!!    Copyright (C) 2015-2015 BigDFT group
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
   subroutine cal_hessian_fd(mhgpsst,runObj,outs,pos,hess,nfree)
      use module_base, only: gp, f_malloc, f_free, assignment(=)
      use module_energyandforces, only: mhgpsenergyandforces
      use bigdft_run, only: run_objects, state_properties
      use module_mhgps_state
   use module_Atoms, only: move_this_coordinate
      implicit none
      type(mhgps_state), intent(inout) :: mhgpsst
      type(run_objects), intent(inout) :: runObj
      type(state_properties), intent(inout) :: outs
      real(gp), dimension(3*runObj%atoms%astruct%nat), intent(in) :: pos
      real(gp), dimension(3*runObj%atoms%astruct%nat,3*runObj%atoms%astruct%nat), intent(inout) :: hess
      !local variables
      integer :: iat
      !real(gp) :: t1,t2,t3
      !real(gp), allocatable, dimension(:,:) :: hess
      real(gp), allocatable, dimension(:) :: tpos,grad,eval,workf
      real(gp) :: h,twelfth,twothird,etot,cmx,cmy,cmz,dm,tt
      real(gp) :: s
      integer :: i,j,k,lworkf,infocode,idir,jat,jdir
      integer, dimension(:), allocatable :: ifrztyp0 !< To avoid to freeze the atoms for bigdft_state
      integer :: ifree, jfree,nfree

      !allocate(hess(3*runObj%atoms%astruct%nat,3*runObj%atoms%astruct%nat))
      tpos = f_malloc(3*runObj%atoms%astruct%nat,id='tpos')
      grad = f_malloc(3*runObj%atoms%astruct%nat,id='grad')
      eval = f_malloc(3*runObj%atoms%astruct%nat,id='eval')
      ifrztyp0 = f_malloc(runObj%atoms%astruct%nat, id = 'ifrztyp0')

      lworkf=1000*runObj%atoms%astruct%nat
      workf = f_malloc(lworkf,id='workf')

      ifrztyp0 = runObj%atoms%astruct%ifrztyp                             
!      runObj%atoms%astruct%ifrztyp = 0 

      !h=1.e-1_gp
      !h=7.5e-2_gp
      !h=5.e-2_gp
         h=1.e-3_gp
      !h=1.e-2_gp
      !h=5.e-3_gp
      !h=2.e-2_gp
      twelfth=-1._gp/(12._gp*h)
      twothird=-2._gp/(3._gp*h)
      hess=0.0_gp
      if(mhgpsst%iproc==0) write(*,*) '(hess) HESSIAN: h',h
      !-------------------------------------------------------
      ifree=0
      nfree=0
      do i=1,3*runObj%atoms%astruct%nat
         iat=(i-1)/3+1
         idir=i-3*(iat-1)
         if (.not.move_this_coordinate(ifrztyp0(iat),idir)) then          
            cycle                                                      
         end if 
         nfree=nfree+1
         ifree=ifree+1
         do j=1,3*runObj%atoms%astruct%nat
             tpos(j)=pos(j)
             grad(j)=0._gp
         enddo
         !-----------------------------------------
         tpos(i)=tpos(i)-2*h
         call mhgpsenergyandforces(mhgpsst,runObj,outs,tpos,grad,etot,infocode)
         jfree=0
         do j=1,3*runObj%atoms%astruct%nat
         jat=(j-1)/3+1
         jdir=j-3*(jat-1)
         if (.not.move_this_coordinate(ifrztyp0(jat),jdir)) then          
            cycle                                                      
         end if 
             jfree=jfree+1
             hess(jfree,ifree)=twelfth*grad(j)
!             hess(j,i)=twelfth*grad(j)
         enddo
         !if(mhgpsst%iproc==0) write(*,*) 'ALIREZA-6',i,iat
         !-----------------------------------------
         tpos(i)=tpos(i)+h
         call mhgpsenergyandforces(mhgpsst,runObj,outs,tpos,grad,etot,infocode)
         jfree=0
         do j=1,3*runObj%atoms%astruct%nat
         jat=(j-1)/3+1
         jdir=j-3*(jat-1)
         if (.not.move_this_coordinate(ifrztyp0(jat),jdir)) then          
            cycle                                                      
         end if 
             jfree=jfree+1
         hess(jfree,ifree)=hess(jfree,ifree)-twothird*grad(j)
!         hess(j,i)=hess(j,i)-twothird*grad(j)
         enddo
         !-----------------------------------------
         tpos(i)=tpos(i)+2*h
         call mhgpsenergyandforces(mhgpsst,runObj,outs,tpos,grad,etot,infocode)
         jfree=0
         do j=1,3*runObj%atoms%astruct%nat
         jat=(j-1)/3+1
         jdir=j-3*(jat-1)
         if (.not.move_this_coordinate(ifrztyp0(jat),jdir)) then          
            cycle                                                      
         end if 
             jfree=jfree+1
         hess(jfree,ifree)=hess(jfree,ifree)+twothird*grad(j)
!         hess(j,i)=hess(j,i)+twothird*grad(j)
         enddo
         !-----------------------------------------
         tpos(i)=tpos(i)+h
         call mhgpsenergyandforces(mhgpsst,runObj,outs,tpos,grad,etot,infocode)
         jfree=0
         do j=1,3*runObj%atoms%astruct%nat
         jat=(j-1)/3+1
         jdir=j-3*(jat-1)
         if (.not.move_this_coordinate(ifrztyp0(jat),jdir)) then          
            cycle                                                      
         end if 
             jfree=jfree+1
         hess(jfree,ifree)=hess(jfree,ifree)-twelfth*grad(j)
!         hess(j,i)=hess(j,i)-twelfth*grad(j)
         !write(*,*) 'HESS ',j,i,hess(j,i)
         enddo
         !-----------------------------------------
      enddo
      !-------------------------------------------------------

      !check symmetry
      dm=0._gp
      do i=1,nfree
!      do i=1,3*runObj%atoms%astruct%nat
         do j=1,i-1
            s=.5_gp*(hess(i,j)+hess(j,i))
            tt=abs(hess(i,j)-hess(j,i))/(1._gp+abs(s))
            dm=max(dm,tt)
            hess(i,j)=s
            hess(j,i)=s
         enddo
      enddo
      if (dm.gt.1.e-1_gp) write(*,*) '(hess) max dev from sym',dm

   !   do j=1,3*runObj%atoms%astruct%nat
   !   do i=1,3*runObj%atoms%astruct%nat
   !   write(*,*) '(hess) hier',runObj%atoms%astruct%nat,hess(i,j)
   !   write(499,*) hess(i,j)
   !   enddo
   !   enddo

      !-------------------------------------------------------
      !project out rotations
      cmx=0._gp ; cmy=0._gp ; cmz=0._gp
      do i=1,3*runObj%atoms%astruct%nat-2,3
         cmx=cmx+pos(i+0)
         cmy=cmy+pos(i+1)
         cmz=cmz+pos(i+2)
      enddo
      cmx=cmx/runObj%atoms%astruct%nat ; cmy=cmy/runObj%atoms%astruct%nat ; cmz=cmz/runObj%atoms%astruct%nat
     
      !x-y plane
      do i=1,3*runObj%atoms%astruct%nat-2,3
         workf(i+1)= (pos(i+0)-cmx)
         workf(i+0)=-(pos(i+1)-cmy)
      enddo
      runObj%atoms%astruct%ifrztyp = ifrztyp0 
      call f_free(tpos)
      call f_free(grad)
      call f_free(eval)
      call f_free(workf)
      call f_free(ifrztyp0)
   end subroutine cal_hessian_fd

end module module_hessian
