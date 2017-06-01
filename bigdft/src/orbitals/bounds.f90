module bounds
  use locregs
  implicit none

  private

  !> Public routines
  public :: locreg_bounds
  public :: wfd_to_logrids
  public :: make_bounds
  public :: make_all_ib
  public :: ext_buffers,ext_buffers_coarse
  public :: make_bounds_per
  public :: make_all_ib_per
  public :: geocode_buffers
  public :: locreg_mesh_origin,locreg_mesh_shape

  contains

    !> Calculates the bounds arrays needed for convolutions
    subroutine locreg_bounds(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,wfd,bounds)
      use module_base
      !use module_interfaces, except_this_one => locreg_bounds
      implicit none
      !Arguments
      integer, intent(in) :: n1,n2,n3
      integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
      type(wavefunctions_descriptors), intent(in) :: wfd
      type(convolutions_bounds), intent(out) :: bounds
      !Local variables
      character(len=*), parameter :: subname='locreg_bounds'
      logical, dimension(:,:,:), allocatable :: logrid_c,logrid_f
    
      call f_routine(id=subname)
    
      !define logrids
      logrid_c = f_malloc((/ 0.to.n1, 0.to.n2, 0.to.n3 /),id='logrid_c')
      logrid_f = f_malloc((/ 0.to.n1, 0.to.n2, 0.to.n3 /),id='logrid_f')
      
      call wfd_to_logrids(n1,n2,n3,wfd,logrid_c,logrid_f)
    
      !allocate and calculate kinetic bounds
      bounds%kb%ibyz_c = f_malloc_ptr((/ 1.to.2,0.to.n2,0.to.n3/),id='bounds%kb%ibyz_c')
      bounds%kb%ibxz_c = f_malloc_ptr((/ 1.to.2,0.to.n1,0.to.n3/),id='bounds%kb%ibxz_c')
      bounds%kb%ibxy_c = f_malloc_ptr((/ 1.to.2,0.to.n1,0.to.n2/),id='bounds%kb%ibxy_c')
      bounds%kb%ibyz_f = f_malloc_ptr((/ 1.to.2,0.to.n2,0.to.n3/),id='bounds%kb%ibyz_f')
      bounds%kb%ibxz_f = f_malloc_ptr((/ 1.to.2,0.to.n1,0.to.n3/),id='bounds%kb%ibxz_f')
      bounds%kb%ibxy_f = f_malloc_ptr((/ 1.to.2,0.to.n1,0.to.n2/),id='bounds%kb%ibxy_f')
    
      call make_bounds(n1,n2,n3,logrid_c,bounds%kb%ibyz_c,bounds%kb%ibxz_c,bounds%kb%ibxy_c)
      call make_bounds(n1,n2,n3,logrid_f,bounds%kb%ibyz_f,bounds%kb%ibxz_f,bounds%kb%ibxy_f)
    
      call f_free(logrid_c)
      call f_free(logrid_f)
      
      !allocate grow, shrink and real bounds
      bounds%gb%ibzxx_c = f_malloc_ptr((/ 1.to.2, 0.to.n3, -14.to.2*n1+16 /),id='bounds%gb%ibzxx_c')
      bounds%gb%ibxxyy_c = f_malloc_ptr((/ 1.to.2, -14.to.2*n1+16, -14.to.2*n2+16 /),id='bounds%gb%ibxxyy_c')
      bounds%gb%ibyz_ff = f_malloc_ptr((/ 1.to.2, nfl2.to.nfu2, nfl3.to.nfu3 /),id='bounds%gb%ibyz_ff')
      bounds%gb%ibzxx_f = f_malloc_ptr((/ 1.to.2, nfl3.to.nfu3, 2*nfl1-14.to.2*nfu1+16 /),id='bounds%gb%ibzxx_f')
      bounds%gb%ibxxyy_f = f_malloc_ptr((/ 1.to.2, 2*nfl1-14.to.2*nfu1+16, 2*nfl2-14.to.2*nfu2+16 /),id='bounds%gb%ibxxyy_f')
    
      bounds%sb%ibzzx_c = f_malloc_ptr((/ 1.to.2 , -14.to.2*n3+16 , 0.to.n1 /),id='bounds%sb%ibzzx_c')
      bounds%sb%ibyyzz_c = f_malloc_ptr((/ 1.to.2, -14.to.2*n2+16, -14.to.2*n3+16 /),id='bounds%sb%ibyyzz_c')
      bounds%sb%ibxy_ff = f_malloc_ptr((/ 1.to.2, nfl1.to.nfu1, nfl2.to.nfu2 /),id='bounds%sb%ibxy_ff')
      bounds%sb%ibzzx_f = f_malloc_ptr((/ 1.to.2, -14+2*nfl3.to.2*nfu3+16, nfl1.to.nfu1 /),id='bounds%sb%ibzzx_f')
      bounds%sb%ibyyzz_f = f_malloc_ptr((/ 1.to.2, -14+2*nfl2.to.2*nfu2+16, -14+2*nfl3.to.2*nfu3+16 /),id='bounds%sb%ibyyzz_f')
    
      bounds%ibyyzz_r = f_malloc_ptr((/ 1.to.2,-14.to.2*n2+16,-14.to.2*n3+16/),id='bounds%ibyyzz_r')
    
      call make_all_ib(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
           bounds%kb%ibxy_c,bounds%sb%ibzzx_c,bounds%sb%ibyyzz_c,&
           bounds%kb%ibxy_f,bounds%sb%ibxy_ff,bounds%sb%ibzzx_f,bounds%sb%ibyyzz_f,&
           bounds%kb%ibyz_c,bounds%gb%ibzxx_c,bounds%gb%ibxxyy_c,&
           bounds%kb%ibyz_f,bounds%gb%ibyz_ff,bounds%gb%ibzxx_f,bounds%gb%ibxxyy_f,&
           bounds%ibyyzz_r)
    
      call f_release_routine()
    
    END SUBROUTINE locreg_bounds

    subroutine wfd_to_logrids(n1,n2,n3,wfd,logrid_c,logrid_f)
      use module_base
      use locregs
      implicit none
      !Arguments
      integer, intent(in) :: n1,n2,n3
      type(wavefunctions_descriptors), intent(in) :: wfd
      logical, dimension(0:n1,0:n2,0:n3), intent(out) :: logrid_c,logrid_f
      !local variables
      integer :: iseg,j0,j1,ii,i1,i2,i3,i0,nvctr_check,i,n1p1,np
    
      n1p1=n1+1
      np=n1p1*(n2+1)
    
      !coarse part
      logrid_c(:,:,:)=.false.
      !control variable
      nvctr_check=0
      !write(*,*) 'np, n1p1', np, n1p1
      do iseg=1,wfd%nseg_c
         j0=wfd%keygloc(1,iseg)
         j1=wfd%keygloc(2,iseg)
         ii=j0-1
         i3=ii/np
         ii=ii-i3*np
         i2=ii/n1p1
         i0=ii-i2*n1p1
         i1=i0+j1-j0
         do i=i0,i1
            nvctr_check=nvctr_check+1
            logrid_c(i,i2,i3)=.true.
         end do
      end do
      !check
      if (nvctr_check /= wfd%nvctr_c) then
         write(*,'(1x,a,3(i6))')&
              'ERROR: problem in wfd_to_logrid(coarse)',nvctr_check,wfd%nvctr_c,wfd%nseg_c
         stop
      end if
      !!do i3=0,n3
      !!  do i2=0,n2
      !!    do i1=0,n1
      !!      write(700,'(a,3i9,l)') 'i1, i2, i3, logrid_c(i1,i2,i3)', i1, i2, i3, logrid_c(i1,i2,i3)
      !!    end do
      !!  end do
      !!end do
    
      !fine part
      logrid_f(:,:,:)=.false.
      !control variable
      nvctr_check=0
      do iseg=wfd%nseg_c+1,wfd%nseg_c+wfd%nseg_f
         j0=wfd%keygloc(1,iseg)
         j1=wfd%keygloc(2,iseg)
         ii=j0-1
         i3=ii/np
         ii=ii-i3*np
         i2=ii/n1p1
         i0=ii-i2*n1p1
         i1=i0+j1-j0
         do i=i0,i1
            nvctr_check=nvctr_check+1
            logrid_f(i,i2,i3)=.true.
         end do
      end do
      !check
      if (nvctr_check /= wfd%nvctr_f) then
         write(*,'(1x,a,2(i6))')&
              'ERROR: problem in wfd_to_logrid(fine)',nvctr_check,wfd%nvctr_f
         stop
      end if
    
    END SUBROUTINE wfd_to_logrids


    subroutine make_bounds(n1,n2,n3,logrid,ibyz,ibxz,ibxy)
       implicit none
       integer, intent(in) :: n1,n2,n3
       logical, dimension(0:n1,0:n2,0:n3), intent(in) :: logrid
       integer, dimension(2,0:n2,0:n3), intent(out) :: ibyz
       integer, dimension(2,0:n1,0:n3), intent(out) :: ibxz
       integer, dimension(2,0:n1,0:n2), intent(out) :: ibxy
       !local variables
       integer :: i1,i2,i3
    
       !$omp parallel default(shared) private(i3,i2,i1)
       !$omp do
       do i3=0,n3 
          do i2=0,n2 
             ibyz(1,i2,i3)= 1000
             ibyz(2,i2,i3)=-1000
    
             loop_i1s: do i1=0,n1
                if (logrid(i1,i2,i3)) then 
                   ibyz(1,i2,i3)=i1
                   exit loop_i1s
                endif
             enddo loop_i1s
    
             loop_i1e: do i1=n1,0,-1
                if (logrid(i1,i2,i3)) then 
                   ibyz(2,i2,i3)=i1
                   exit loop_i1e
                endif
             enddo loop_i1e
          end do
       end do
       !$omp end do
    
       !$omp do
       do i3=0,n3 
          do i1=0,n1
             ibxz(1,i1,i3)= 1000
             ibxz(2,i1,i3)=-1000
    
             loop_i2s: do i2=0,n2 
                if (logrid(i1,i2,i3)) then 
                   ibxz(1,i1,i3)=i2
                   exit loop_i2s
                endif
             enddo loop_i2s
    
             loop_i2e: do i2=n2,0,-1
                if (logrid(i1,i2,i3)) then 
                   ibxz(2,i1,i3)=i2
                   exit loop_i2e
                endif
             enddo loop_i2e
    
          end do
       end do
       !$omp end do
    
       !$omp do
       do i2=0,n2 
          do i1=0,n1 
             ibxy(1,i1,i2)= 1000
             ibxy(2,i1,i2)=-1000
    
             loop_i3s: do i3=0,n3
                if (logrid(i1,i2,i3)) then 
                   ibxy(1,i1,i2)=i3
                   exit loop_i3s
                endif
             enddo loop_i3s
    
             loop_i3e: do i3=n3,0,-1
                if (logrid(i1,i2,i3)) then 
                   ibxy(2,i1,i2)=i3
                   exit loop_i3e
                endif
             enddo loop_i3e
          end do
       end do
       !$omp end do
       !$omp end parallel
    
    END SUBROUTINE make_bounds


    !> Cleaned version of the logrid_old.f90 in the unused directory (with newmethod=.true.)
    !! Creates complicated ib arrays    
    subroutine make_all_ib(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
         ibxy_c,ibzzx_c,ibyyzz_c,ibxy_f,ibxy_ff,ibzzx_f,ibyyzz_f,&
         ibyz_c,ibzxx_c,ibxxyy_c,ibyz_f,ibyz_ff,ibzxx_f,ibxxyy_f,ibyyzz_r)
      use module_base
      !use module_interfaces, except_this_one => make_all_ib
      implicit none
      integer,intent(in)::n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
      integer :: i1,i2,i3 !n(c) m1,m2,m3
    
      integer,intent(in):: ibyz_c(2,0:n2,0:n3),ibxy_c(2,0:n1,0:n2)
      integer,intent(in):: ibyz_f(2,0:n2,0:n3),ibxy_f(2,0:n1,0:n2)
    
      !    for shrink:    
      integer,intent(inout):: ibzzx_c(2,-14:2*n3+16,0:n1) 
      integer,intent(out):: ibyyzz_c(2,-14:2*n2+16,-14:2*n3+16)
    
      integer,intent(out):: ibxy_ff(2,nfl1:nfu1,nfl2:nfu2)
      integer,intent(inout):: ibzzx_f(2,-14+2*nfl3:2*nfu3+16,nfl1:nfu1) 
      integer,intent(out):: ibyyzz_f(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16)
    
      !    for grow:    
      integer,intent(out):: ibzxx_c(2,0:n3,-14:2*n1+16) ! extended boundary arrays
      integer,intent(out):: ibxxyy_c(2,-14:2*n1+16,-14:2*n2+16)
    
      integer,intent(inout):: ibyz_ff(2,nfl2:nfu2,nfl3:nfu3)
      integer,intent(out):: ibzxx_f(2,nfl3:nfu3,2*nfl1-14:2*nfu1+16)
      integer,intent(out):: ibxxyy_f(2,2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16)
    
      character(len=*), parameter :: subname='make_all_ib'
      logical,allocatable:: logrid_big(:)
    
      !    for real space:
      integer,intent(out):: ibyyzz_r(2,-14:2*n2+16,-14:2*n3+16)
    
    
      logrid_big = f_malloc((2*n1+31)*(2*n2+31)*(2*n3+31),id='logrid_big')
    
      !n(c) m1=nfu1-nfl1
      !n(c) m2=nfu2-nfl2
      !n(c) m3=nfu3-nfl3
    
      !   (0:n3,-14:2*n1+16,-14:2*n2+16) from grow
      !   (-14:2*n2+16,-14:2*n3+16,0:n1) from shrink
      !   (-14:2*n1+16,-14:2*n2+16,-14:2*n3+16) from real space
    
      !	for shrink:
      do i2=nfl2,nfu2
         do i1=nfl1,nfu1
            ibxy_ff(:,i1,i2)=ibxy_f(:,i1,i2)
         enddo
      enddo
    
      call make_ib_inv(logrid_big, ibxy_c,ibzzx_c,ibyyzz_c,0,n1,0,n2,0,n3)
      call make_ib_inv(logrid_big,ibxy_ff,ibzzx_f,ibyyzz_f,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3)
    
      !for realspace:
      !-14:2*n2+16,-14:2*n3+16
      do i3=-14,2*n3+16
         do i2=-14,2*n2+16
            if (ibyyzz_c(1,i2,i3).ne.1000) then
               ibyyzz_r(1,i2,i3)=2*ibyyzz_c(1,i2,i3)
               ibyyzz_r(2,i2,i3)=2*ibyyzz_c(2,i2,i3)+30
            else
               ibyyzz_r(1,i2,i3)=1000
               ibyyzz_r(2,i2,i3)=-1000
            endif
         enddo
      enddo
      call squares(ibyyzz_r,2*n2+30,2*n3+30)
    
      !    for grow:
    
      do i2=nfl2,nfu2
         do i3=nfl3,nfu3
            ibyz_ff(:,i2,i3)=ibyz_f(:,i2,i3)
         enddo
      enddo
    
      call make_ib_c(logrid_big,ibyz_c,ibzxx_c,ibxxyy_c,n1,n2,n3)
    
      call make_ib(logrid_big,ibyz_ff,ibzxx_f,ibxxyy_f,&
           nfl1,nfu1,nfl2,nfu2,nfl3,nfu3)
    
      call squares_1d(ibxxyy_f,2*nfl1-14,2*nfu1+16,2*nfl2-14,2*nfu2+16)
    
      call f_free(logrid_big)
    
    END SUBROUTINE make_all_ib
    
    
    !> This subroutine mimics the comb_grow_f one
    subroutine make_ib_inv(logrid_big,ibxy,ibzzx,ibyyzz,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3)
      !use module_interfaces, except_this_one => make_ib_inv
      implicit none
      integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
      integer,intent(in):: ibxy(2,nfl1:nfu1,nfl2:nfu2)
      integer,intent(inout):: ibzzx(2,-14+2*nfl3:2*nfu3+16,nfl1:nfu1) 
      integer,intent(out):: ibyyzz(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16)
      logical, intent(inout) :: logrid_big(nfl3:nfu3,2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16)! work array
      integer :: nt
    
      ! I3,i1,i2 -> i1,i2,i3 
      nt=(nfu1-nfl1+1)*(nfu2-nfl2+1)
      call ib_to_logrid_inv(ibxy,logrid_big,nfl3,nfu3,nt)
    
      ! I2,I3,i1 -> I3,i1,i2
      nt=(2*(nfu3-nfl3)+31)*(nfu1-nfl1+1)
      call ib_from_logrid_inv(ibzzx,logrid_big,nfl2,nfu2,nt)
      call ib_to_logrid_inv(ibzzx,logrid_big,nfl2,nfu2,nt)
    
      ! I1,I2,I3  -> I2,I3,i1
      nt=(2*(nfu2-nfl2)+31)*(2*(nfu3-nfl3)+31)
      call ib_from_logrid_inv(ibyyzz,logrid_big,nfl1,nfu1,nt)
    
    END SUBROUTINE make_ib_inv
    
    
    !> This one mimics the comb_rot_grow_f_loc
    subroutine ib_to_logrid_inv(ib,logrid,nfl,nfu,ndat)
      implicit none
      integer, intent(in) :: ndat,nfl,nfu
      integer, intent(in) :: ib(2,ndat)! input
      logical, intent(out) :: logrid(-14+2*nfl:2*nfu+16,ndat)! output
    
      integer :: l,i
    
      logrid=.false.
    
      !$omp parallel default(shared) private(i,l)
      !$omp do
      do l=1,ndat
         do i = 2*ib(1,l)-14 , 2*ib(2,l)+16
            logrid(i,l)=.true.
         enddo
      enddo
      !$omp end do
      !$omp end parallel
    
    END SUBROUTINE ib_to_logrid_inv
    
    
    !> Mimics the bounds subroutine    
    subroutine ib_from_logrid_inv(ib,logrid,ml1,mu1,ndat)
      implicit none
      integer, intent(in) :: ml1,mu1,ndat
      integer, intent(out) :: ib(2,ndat)
      logical, intent(in) :: logrid(ndat,ml1:mu1)
    
      integer :: i,i1
    
      !$omp parallel default(shared) private(i,i1)
      !$omp do
      do i=1,ndat
         ib(1,i)= 1000
         ib(2,i)=-1000
    
         inner1:do i1=ml1,mu1
            if (logrid(i,i1)) then 
               ib(1,i)=i1
               exit inner1
            endif
         enddo inner1
    
         inner2:do i1=mu1,ml1,-1
            if (logrid(i,i1)) then 
               ib(2,i)=i1
               exit inner2
            endif
         enddo inner2
      enddo
      !$omp end do
      !$omp end parallel
    
    END SUBROUTINE ib_from_logrid_inv
    
    
    !> This subroutine mimics the comb_grow_f one
    subroutine make_ib_c(logrid_big,ibyz,ibzxx,ibxxyy,n1,n2,n3)
      !use module_interfaces, except_this_one => make_ib_c
      implicit none
      integer nt,n1,n2,n3
      integer ibyz(2,0:n2,0:n3)! input
      integer ibzxx(2,0:n3,-14:2*n1+16)!output
      integer ibxxyy(2,-14:2*n1+16,-14:2*n2+16)!output
      logical logrid_big(0:n3,-14:2*n1+16,-14:2*n2+16)! work array
    
      call squares(ibyz,n2,n3)
      ! i1,i2,i3 -> i2,i3,I1
      nt=(n2+1)*(n3+1)
      call ib_to_logrid_rot(ibyz,logrid_big,0,n1,nt)
    
      ! i2,i3,I1 -> i3,I1,I2
      nt=(n3+1)*(2*n1+31)
      call ib_from_logrid(ibzxx,logrid_big,0,n2,nt)
      call squares(ibzxx,n3,2*n1+30)
    
      call ib_to_logrid_rot(ibzxx,logrid_big,0,n2,nt)
    
      ! i3,I1,I2  -> I1,I2,I3
      nt=(2*n1+31)*(2*n2+31)
      call ib_from_logrid(ibxxyy,logrid_big,0,n3,nt)
      call squares(ibxxyy,2*n1+30,2*n2+30)
    
    END SUBROUTINE make_ib_c
    
    
    !> This subroutine mimics the comb_grow_f one
    subroutine make_ib(logrid_big,ibyz,ibzxx,ibxxyy,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3)
      !use module_interfaces, except_this_one => make_ib
      implicit none
      integer nt,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
      integer ibyz(  2,nfl2:nfu2,nfl3:nfu3)! input
      integer ibzxx( 2,          nfl3:nfu3,2*nfl1-14:2*nfu1+16)!output
      integer ibxxyy(2,                    2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16)!output
      logical logrid_big(           nfl3:nfu3,2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16)! work array
    
      ! i1,i2,i3 -> i2,i3,I1
      nt=(nfu2-nfl2+1)*(nfu3-nfl3+1)
      call ib_to_logrid_rot(  ibyz,logrid_big,nfl1,nfu1,nt)
    
      ! i2,i3,I1 -> i3,I1,I2
      nt=(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31)
      call ib_from_logrid( ibzxx,logrid_big,nfl2,nfu2,nt)
      call ib_to_logrid_rot( ibzxx,logrid_big,nfl2,nfu2,nt)
    
      ! i3,I1,I2  -> I1,I2,I3
      nt=(2*(nfu1-nfl1)+31)*(2*(nfu2-nfl2)+31)
      call ib_from_logrid(ibxxyy,logrid_big,nfl3,nfu3,nt)
    
    END SUBROUTINE make_ib
    
    
    !> This one mimics the comb_rot_grow_f_loc
    subroutine ib_to_logrid_rot(ib,logrid,nfl,nfu,ndat)
      implicit none
      integer ndat,nfl,nfu,l,i
      integer ib(2,ndat)! input
      logical logrid(ndat,-14+2*nfl:2*nfu+16)! output
    
      logrid=.false.
    
      !$omp parallel default(shared) private(l,i)
      !$omp do
      do l=1,ndat
         do i = 2*ib(1,l)-14 , 2*ib(2,l)+16
            logrid(l,i)=.true.
         enddo
      enddo
      !$omp end do
      !$omp end parallel
    
    END SUBROUTINE ib_to_logrid_rot
    
    
    !> Mimics the bounds subroutine    
    subroutine ib_from_logrid(ib,logrid,ml1,mu1,ndat)
      implicit none
      integer i,i1
      integer ml1,mu1,ndat
      integer ib(2,ndat)
      logical logrid(ml1:mu1,ndat)
    
      !$omp parallel default(shared) private(i,i1)
      !$omp do
      do i=1,ndat
         ib(1,i)= 1000
         ib(2,i)=-1000
    
         inner1:do i1=ml1,mu1
            if (logrid(i1,i)) then 
               ib(1,i)=i1
               exit inner1
            endif
         enddo inner1
    
         inner2:do i1=mu1,ml1,-1
            if (logrid(i1,i)) then 
               ib(2,i)=i1
               exit inner2
            endif
         enddo inner2
      enddo
      !$omp end do
      !$omp end parallel
    
    END SUBROUTINE ib_from_logrid
    
    
    !> Modifies the ib array
    !! so that it is made up of blocks of size 2
    !! the localization region is enlarged as a result
    !! works for even nfl2 only
    subroutine squares_1d(ib,nfl2,nfu2,nfl3,nfu3)
      implicit none
      !Arguments
      integer,intent(in) :: nfl2,nfu2,nfl3,nfu3
      integer,intent(inout) :: ib(2,nfl2:nfu2,nfl3:nfu3)
      !Local variables
      integer :: i2,i3,ii2,ibmin,ibmax
    
      do i3=nfl3,nfu3
         do i2=nfl2/2,(nfu2-1)/2
            ii2=2*i2
            ibmin=min(ib(1,ii2,i3),ib(1,ii2+1,i3))
    
            ib(1,ii2,i3)=ibmin
            ib(1,ii2+1,i3)=ibmin
    
            ibmax=max(ib(2,ii2,i3),ib(2,ii2+1,i3))
    
            ib(2,ii2,i3)=ibmax
            ib(2,ii2+1,i3)=ibmax
         enddo
      enddo
    END SUBROUTINE squares_1d
    
    
    !> Modifies the ib array 
    !! so that it is made up of squares 2x2
    !! the localization region is enlarged as a result
    subroutine squares(ib,n2,n3)
      implicit none
      !Arguments
      integer, intent(in) :: n2, n3
      integer, dimension(2,0:n2,0:n3), intent(inout) :: ib
      !Local variables
      integer :: i2,i3,ii2,ii3,ibmin,ibmax
    
      !If one dimension is zero: do nothing
      if (n2 == 0 .or. n3 == 0) return
    
      do i3=0,(n3-1)/2
         ii3=2*i3
         do i2=0,(n2-1)/2
            ii2=2*i2
            ibmin=min(ib(1,ii2,ii3),ib(1,ii2+1,ii3),&
                 ib(1,ii2,ii3+1),ib(1,ii2+1,ii3+1))
    
            ib(1,ii2,ii3)=ibmin
            ib(1,ii2+1,ii3)=ibmin
            ib(1,ii2,ii3+1)=ibmin
            ib(1,ii2+1,ii3+1)=ibmin
    
            ibmax=max(ib(2,ii2,ii3),ib(2,ii2+1,ii3),&
                 ib(2,ii2,ii3+1),ib(2,ii2+1,ii3+1))
    
            ib(2,ii2,ii3)=ibmax
            ib(2,ii2+1,ii3)=ibmax
            ib(2,ii2,ii3+1)=ibmax
            ib(2,ii2+1,ii3+1)=ibmax
         enddo
      enddo
    END SUBROUTINE squares

    pure subroutine ext_buffers(periodic,nl,nr)
      implicit none
      logical, intent(in) :: periodic
      integer, intent(out) :: nl,nr
    
      if (periodic) then
         nl=0
         nr=0
      else
         nl=14
         nr=15
      end if
    END SUBROUTINE ext_buffers

    pure subroutine ext_buffers_coarse(periodic,nb)
      implicit none
      logical, intent(in) :: periodic
      integer, intent(out) :: nb
      if (periodic) then
         nb=0
      else
         nb=7
      end if
    END SUBROUTINE ext_buffers_coarse



    !> Define the bounds of wavefunctions for periodic systems
    subroutine make_bounds_per(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,cbounds,wfd)
      use module_base
      use locregs
      implicit none
      type(wavefunctions_descriptors), intent(in) :: wfd
      type(convolutions_bounds),intent(out):: cbounds
      integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
    
      logical,allocatable,dimension(:,:,:) :: logrid
      character(len=*), parameter :: subname='make_bounds_per'
      integer :: nseg_c
    
      cbounds%kb%ibyz_f = f_malloc_ptr((/ 1.to.2, 0.to.n2, 0.to.n3 /),id='cbounds%kb%ibyz_f')
      cbounds%kb%ibxz_f = f_malloc_ptr((/ 1.to.2, 0.to.n1, 0.to.n3 /),id='cbounds%kb%ibxz_f')
      cbounds%kb%ibxy_f = f_malloc_ptr((/ 1.to.2, 0.to.n1, 0.to.n2 /),id='cbounds%kb%ibxy_f')
    
      cbounds%gb%ibyz_ff = f_malloc_ptr((/ 1.to.2, nfl2.to.nfu2, nfl3.to.nfu3 /),id='cbounds%gb%ibyz_ff')
      cbounds%gb%ibzxx_f = f_malloc_ptr((/ 1.to.2, nfl3.to.nfu3, 0.to.2*n1+1 /),id='cbounds%gb%ibzxx_f')
      cbounds%gb%ibxxyy_f = f_malloc_ptr((/ 1.to.2, 0.to.2*n1+1, 0.to.2*n2+1 /),id='cbounds%gb%ibxxyy_f')
    
      cbounds%sb%ibxy_ff = f_malloc_ptr((/ 1.to.2, nfl1.to.nfu1, nfl2.to.nfu2 /),id='cbounds%sb%ibxy_ff')
      cbounds%sb%ibzzx_f = f_malloc_ptr((/ 1.to.2, 0.to.2*n3+1, nfl1.to.nfu1 /),id='cbounds%sb%ibzzx_f')
      cbounds%sb%ibyyzz_f = f_malloc_ptr((/ 1.to.2, 0.to.2*n2+1, 0.to.2*n3+1 /),id='cbounds%sb%ibyyzz_f')
    
      logrid = f_malloc((/ 0.to.n1, 0.to.n2, 0.to.n3 /),id='logrid')
    
      nseg_c=wfd%nseg_c
      call make_logrid_f(n1,n2,n3, & 
           wfd%nseg_f,wfd%keygloc(1:,nseg_c+min(1,wfd%nseg_f):),  & 
           logrid)
    
      call make_bounds(n1,n2,n3,logrid,cbounds%kb%ibyz_f,cbounds%kb%ibxz_f,cbounds%kb%ibxy_f)
    
      call f_free(logrid)
    
    END SUBROUTINE make_bounds_per
    
    
    subroutine make_all_ib_per(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
         ibxy_f,ibxy_ff,ibzzx_f,ibyyzz_f,&
         ibyz_f,ibyz_ff,ibzxx_f,ibxxyy_f)
      !    creates complicated ib arrays    
      use module_base
      implicit none
      integer,intent(in)::n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
      integer :: i1,i2,i3 !n(c) m1,m2,m3
    
      integer,intent(in):: ibyz_f(2,0:n2,0:n3),ibxy_f(2,0:n1,0:n2)
    
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


    subroutine geocode_buffers(geocode_local, geocode_global, nl1, nl2, nl3)
      implicit none
      character(len=1), intent(in) :: geocode_local, geocode_global !< @copydoc poisson_solver::doc::geocode
      integer, intent(out) :: nl1, nl2, nl3
      !local variables
      logical :: perx_local, pery_local, perz_local
      logical :: perx_global,pery_global,perz_global
      integer :: nr1, nr2, nr3
    
      !conditions for periodicity in the three directions
      perx_local=(geocode_local /= 'F')
      pery_local=(geocode_local == 'P')
      perz_local=(geocode_local /= 'F')
      perx_global=(geocode_global /= 'F')
      pery_global=(geocode_global == 'P')
      perz_global=(geocode_global /= 'F')
    
      call ext_buffers(perx_local, nl1, nr1)
      call ext_buffers(pery_local, nl2, nr2)
      call ext_buffers(perz_local, nl3, nr3)
    
      ! If the global box has non-free boundary conditions, the shift is already
      ! contained in nsi1,nsi2,nsi3 and does not need to be subtracted.
      if (perx_global) then
          nl1 = 0
      end if
      if (pery_global) then
          nl2 = 0
      end if
      if (perz_global) then
          nl3 = 0
      end if
    
    end subroutine geocode_buffers

    pure function locreg_mesh_origin(mesh) result(or)
      use box, only: cell,cell_r,cell_periodic_dims
      use module_defs, only: gp
      implicit none
      type(cell), intent(in) :: mesh
      real(gp), dimension(3) :: or
      !local variables
      logical, dimension(3) :: peri
      integer :: nbli,nbri,i
    
      !buffers associated to the geocode
      !conditions for periodicity in the three directions
      peri=cell_periodic_dims(mesh)
      do i=1,3
         call ext_buffers(peri(i),nbli,nbri)
         or(i)=cell_r(mesh,nbli+1,dim=i)
      end do
 
    end function locreg_mesh_origin

    !> return the shapes of the localisation region
    !!useful to allocate the array psifscf
    function locreg_mesh_shape(mesh,highres) result(ndims)
      use box, only: cell,cell_periodic_dims
      implicit none
      type(cell), intent(in) :: mesh
      logical, intent(in), optional :: highres
      integer, dimension(3) :: ndims
      !local variables
      logical :: hr
      integer :: i,nb
      logical, dimension(3) :: peri
      hr=.false.
      if (present(highres)) hr=highres
      peri=cell_periodic_dims(mesh)
      do i=1,3
         call ext_buffers_coarse(peri(i),nb)
         if (hr) then
            ndims(i)=2*(mesh%ndims(i)+1+nb)
         else
            ndims(i)=mesh%ndims(i)+1
         end if
      end do
    end function locreg_mesh_shape

end module bounds
