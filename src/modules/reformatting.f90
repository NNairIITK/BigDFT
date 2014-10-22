!> @file
!!  Module to handle the rotation and translation of a scalar field
!!  according to a transformation
!!
!! @author
!!    Copyright (C) 2013-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Module handling the rotation and trasnlation of a scalar field as wavefunctions
module reformatting
  implicit none
  
  private

  public :: field_rototranslation3D,my_scaling_function4b2B,back_trans_16_reversed

contains

  !> Routine which directly applies the 3D transformation of the rototranslation
  !this routine has to be cleaned and the allocations of the work arrays have to be deplaced
  subroutine field_rototranslation3D(n_phi,nrange_phi,phi_ISF,Rmat,da,&
       centre_old,centre_new,&
       !newz,cost,sint,onemc,&
       iorder,hgrids_old,ndims_old,f_old,&
       hgrids_new,ndims_new,f_new)
    use module_base
    use yaml_output
    implicit none
    integer, intent(in) :: n_phi,nrange_phi !< number of points of ISF array and real-space range
    integer, dimension(3), intent(in) :: iorder
    real(gp), dimension(3), intent(in) :: da !<coordinates of rigid shift vector
    real(gp), dimension(3), intent(in) :: centre_old,centre_new !<centre of rotation
    real(gp), dimension(3), intent(in) :: hgrids_old,hgrids_new !<dimension of old and new box
    integer, dimension(3), intent(in) :: ndims_old,ndims_new !<dimension of old and new box
    !real(gp), dimension(3) :: newz
    !real(gp), intent(in) ::  cost,sint,onemc
    real(gp), dimension(3,3), intent(in) :: Rmat !<rotation matrix
    real(gp), dimension(n_phi,2), intent(in) :: phi_ISF
    real(gp), dimension(ndims_old(1),ndims_old(2),ndims_old(3)), intent(in) :: f_old
    real(gp), dimension(ndims_new(1),ndims_new(2),ndims_new(3)), intent(out) :: f_new
    !local variables
    integer :: m_isf,k1,i,j,k,me,ms
    real(gp) :: dt,norm!,ux,uy,uz
    !$ real(gp) :: tt
    real(gp) :: scal !<scaling factor
    !real(gp), dimension(3,3) :: Rmat !<rotation matrix
    integer, dimension(3) :: isign,irp
    real(gp), dimension(:), allocatable :: shf
    real(gp), dimension(:), allocatable :: work,work2

    !print *,'3d'
    call f_routine(id='field_rototranslation3D')
    work =f_malloc(ndims_new(1)*(maxval(ndims_old))**2,id='work')
    work2=f_malloc(ndims_new(1)*ndims_new(2)*maxval(ndims_old),id='work2')

    m_isf=nrange_phi/2
    !shf=f_malloc(-m_isf .to. m_isf,id='shf')
    !for each of the dimensions build the interpolating vector which is needed

!!$    !identify the rotation matrix elements
!!$    ux=newz(1)
!!$    uy=newz(2)
!!$    uz=newz(3)
!!$    !first row (xp)
!!$    rmat(:,1)=(/cost+onemc*ux**2   ,ux*uy*onemc-uz*sint,ux*uz*onemc+uy*sint/)
!!$    !second row (yp)
!!$    rmat(:,2)=(/ux*uy*onemc+uz*sint,cost+onemc*uy**2   ,uy*uz*onemc-ux*sint/)
!!$    !third row (zp)
!!$    rmat(:,3)=(/ux*uz*onemc-uy*sint,uy*uz*onemc+ux*sint,cost+onemc*uz**2   /)


    !first step: determine xn from a coordinate n13o=xo or yo or zo
    !f_old (nxo,nyo,nzo) -> work(n11o,n12o,nxn) !n11o and n12o are the remaining dimensions
    !second step: determine yn from n22o=n11o or n12o
    !work(n11o,n12o,nxn) -> work2(n21o,nxn,nyn)
    !third step: determine zn from n21o
    !work2(n21o,nxn,nyn) -> f_new(xn,yn,zn)

    isign=1
    if (iorder(1)<0) isign(1)=2
    if (iorder(2)<0) isign(2)=2
    if (iorder(3)<0) isign(3)=2
    irp=abs(iorder)

    norm=0.0_gp

    !first step
    select case(irp(1))
    case(1) !xn is derived from xo 
       !$omp parallel default(shared) private(scal,k,j,i,dt,k1,ms,me,shf)
       allocate(shf(-m_isf:m_isf))
       !$omp do
       do k=1,ndims_old(3)
          do j=1,ndims_old(2)
             do i=1,ndims_new(1)
                call shift_and_start(irp(1),1,2,3,i,j,k,&
                     dt,k1,ms,me,scal)

                call define_filter(dt,nrange_phi,n_phi,phi_ISF(1,isign(1)),shf)
                call redirect(m_isf,ms,me,shf,scal)
                !print *,'filter',sum(shf),dt,sum(shf**2),'vals',shf
                !work(j,k+(i-1)*ndims_old(3))
                work(j+ind(2,3,k,i))=convolve(irp(1),k1,j,k,ms,me,&
                     m_isf,shf,ndims_old(1),ndims_old(2),ndims_old(3),f_old,scal)
             end do
          end do
       end do
       !$omp end do
       deallocate(shf)
       !$omp end parallel
    case(2) !xn is derived from yo
       !$omp parallel default(shared) private(scal,k,j,i,dt,k1,ms,me,shf)
       allocate(shf(-m_isf:m_isf))
       !$omp do
       do k=1,ndims_old(3)
          do j=1,ndims_old(1)
             do i=1,ndims_new(1)
                call shift_and_start(irp(1),1,1,3,i,j,k,&
                     dt,k1,ms,me,scal)

                call define_filter(dt,nrange_phi,n_phi,phi_ISF(1,isign(1)),shf)
                call redirect(m_isf,ms,me,shf,scal)
                !work(j,k+(i-1)*ndims_old(3))
                work(j+ind(1,3,k,i))=convolve(irp(1),j,k1,k,ms,me,&
                     m_isf,shf,ndims_old(1),ndims_old(2),ndims_old(3),f_old,scal)
             end do
          end do
       end do
       !$omp end do
       deallocate(shf)
       !$omp end parallel
    case(3) !xn is derived from zo
       !$omp parallel default(shared) private(scal,k,j,i,dt,k1,ms,me,shf)
       allocate(shf(-m_isf:m_isf))
       !$omp do
       do k=1,ndims_old(2)
          do j=1,ndims_old(1)
             do i=1,ndims_new(1)
                call shift_and_start(irp(1),1,1,2,i,j,k,&
                     dt,k1,ms,me,scal)

                call define_filter(dt,nrange_phi,n_phi,phi_ISF(1,isign(1)),shf)
                call redirect(m_isf,ms,me,shf,scal)
                !work(k,j+(i-1)*ndims_old(2))
                work(j+ind(1,2,k,i))=convolve(irp(1),j,k,k1,ms,me,&
                     m_isf,shf,ndims_old(1),ndims_old(2),ndims_old(3),f_old,scal)
             end do
          end do
       end do
       !$omp end do
       deallocate(shf)
       !$omp end parallel
    end select

    !second step
    select case(irp(1)*10+irp(2))
    case(21) !yp is derived from xo (and xp has been derived from y)
       !$omp parallel default(shared) private(scal,i,k,j,dt,k1,ms,me,shf)
       allocate(shf(-m_isf:m_isf))
       !$omp do
       do i=1,ndims_new(1)
          do k=1,ndims_old(irp(3))
             do j=1,ndims_new(2)
                call shift_and_start(irp(2),2,2,irp(3),i,j,k,&
                     dt,k1,ms,me,scal)
                call define_filter(dt,nrange_phi,n_phi,phi_ISF(1,isign(2)),shf)
                call redirect(m_isf,ms,me,shf,scal)
                work2(k+ind2(irp(3),i,j))=convolve(1,k1,k,i,ms,me,m_isf,shf,&
                     ndims_old(1),ndims_old(3),ndims_new(1),work,scal)
             end do
          end do
       end do
       !$omp end do
       deallocate(shf)
       !$omp end parallel
    case(23) !yp is derived from zo (and xp has been derived from y)
       !$omp parallel default(shared) private(scal,i,k,j,dt,k1,ms,me,shf)
       allocate(shf(-m_isf:m_isf))
       !$omp do
       do i=1,ndims_new(1)
          do k=1,ndims_old(irp(3))
             do j=1,ndims_new(2)
                call shift_and_start(irp(2),2,2,irp(3),i,j,k,&
                     dt,k1,ms,me,scal)

                call define_filter(dt,nrange_phi,n_phi,phi_ISF(1,isign(2)),shf)
                call redirect(m_isf,ms,me,shf,scal)
                work2(k+ind2(irp(3),i,j))=convolve(2,k,k1,i,ms,me,m_isf,shf,&
                     ndims_old(1),ndims_old(3),ndims_new(1),work,scal)
             end do
          end do
       end do
       !$omp end do
       deallocate(shf)
       !$omp end parallel
    case(12) !yp is derived from yo (and xp has been derived from x)
       !$omp parallel default(shared) private(scal,i,k,j,dt,k1,ms,me,shf)
       allocate(shf(-m_isf:m_isf))
       !$omp do
       do i=1,ndims_new(1)
          do k=1,ndims_old(irp(3))
             do j=1,ndims_new(2)
                call shift_and_start(irp(2),2,2,irp(3),i,j,k,&
                     dt,k1,ms,me,scal)
                call define_filter(dt,nrange_phi,n_phi,phi_ISF(1,isign(2)),shf)
                call redirect(m_isf,ms,me,shf,scal)
                !work2(k,i+(j-1)*ndims_new(1))
                work2(k+ind2(irp(3),i,j))=convolve(1,k1,k,i,ms,me,m_isf,shf,&
                     ndims_old(2),ndims_old(3),ndims_new(1),work,scal)
             end do
          end do
       end do
       !$omp end do
       deallocate(shf)
       !$omp end parallel
    case(13) !yp is derived from zo (and xp has been derived from x)
       !$omp parallel default(shared) private(scal,i,k,j,dt,k1,ms,me,shf)
       allocate(shf(-m_isf:m_isf))
       !$omp do
       do i=1,ndims_new(1)
          do k=1,ndims_old(irp(3))
             do j=1,ndims_new(2)
                call shift_and_start(irp(2),2,2,irp(3),i,j,k,&
                     dt,k1,ms,me,scal)

                !              print *,'value fouund',dt,k1,j

                call define_filter(dt,nrange_phi,n_phi,phi_ISF(1,isign(2)),shf)
                call redirect(m_isf,ms,me,shf,scal)
                !work2(k,i+(j-1)*ndims_new(1))
                work2(k+ind2(irp(3),i,j))=convolve(2,k,k1,i,ms,me,m_isf,shf,&
                     ndims_old(2),ndims_old(3),ndims_new(1),work,scal)
             end do
          end do
       end do
       !$omp end do
       deallocate(shf)
       !$omp end parallel
    case(32) !yp is derived from yo (and xp has been derived from z)
       !$omp parallel default(shared) private(scal,i,k,j,dt,k1,ms,me,shf)
       allocate(shf(-m_isf:m_isf))
       !$omp do
       do i=1,ndims_new(1)
          do k=1,ndims_old(irp(3))
             do j=1,ndims_new(2)
                call shift_and_start(irp(2),2,2,irp(3),i,j,k,&
                     dt,k1,ms,me,scal)

                call define_filter(dt,nrange_phi,n_phi,phi_ISF(1,isign(2)),shf)
                call redirect(m_isf,ms,me,shf,scal)
                work2(k+ind2(irp(3),i,j))=convolve(2,k,k1,i,ms,me,m_isf,shf,&
                     ndims_old(1),ndims_old(2),ndims_new(1),work,scal)
             end do
          end do
       end do
       !$omp end do
       deallocate(shf)
       !$omp end parallel
    case(31) !yp is derived from xo (and xp has been derived from z)
       !$omp parallel default(shared) private(scal,i,k,j,dt,k1,ms,me,shf)
       allocate(shf(-m_isf:m_isf))
       !$omp do
       do i=1,ndims_new(1)
          do k=1,ndims_old(irp(3))
             do j=1,ndims_new(2)
                call shift_and_start(irp(2),2,2,irp(3),i,j,k,&
                     dt,k1,ms,me,scal)

                call define_filter(dt,nrange_phi,n_phi,phi_ISF(1,isign(2)),shf)
                call redirect(m_isf,ms,me,shf,scal)
                work2(k+ind2(irp(3),i,j))=convolve(1,k1,k,i,ms,me,m_isf,shf,&
                     ndims_old(1),ndims_old(2),ndims_new(1),work,scal)
             end do
          end do
       end do
       !$omp end do
       deallocate(shf)
       !$omp end parallel
    end select

    !third step
    !$omp parallel default(shared) private(scal,j,i,k,dt,k1,ms,me,shf,tt)
    allocate(shf(-m_isf:m_isf)) !the usage of f_malloc here can be considered
    !$omp do !!reduction(+:norm)
    do j=1,ndims_new(2)
       do i=1,ndims_new(1)
          do k=1,ndims_new(3)
             call shift_and_start(irp(3),3,2,3,i,j,k,&
                  dt,k1,ms,me,scal)

             call define_filter(dt,nrange_phi,n_phi,phi_ISF(1,isign(3)),shf)
             call redirect(m_isf,ms,me,shf,scal)
             f_new(i,j,k)=convolve(1,k1,i,j,ms,me,m_isf,shf,&
                  ndims_old(irp(3)),ndims_new(1),ndims_new(2),work2,scal)
             !tt=f_new(i,j,k)
             !tt=tt*tt
             !norm=norm+tt
          end do
       end do
    end do
    !$omp end do
    deallocate(shf)
    !$omp end parallel

    !print *,'norm of the function',sqrt(norm)

    call f_free(work)
    call f_free(work2)
    !call f_free(shf)
    call f_release_routine()

  contains

    !index of work array for step 1
    pure function ind(jc2,jc3,i2,i3)
      implicit none
      integer, intent(in) :: jc2,jc3,i2,i3
      integer :: ind

      ind=ndims_old(jc2)*(i2-1)+ndims_old(jc2)*ndims_old(jc3)*(i3-1)

    end function ind

    pure function ind2(jc3,i2,i3)
      implicit none
      integer, intent(in) :: jc3,i2,i3
      integer :: ind2

      ind2=ndims_old(jc3)*(i2-1)+ndims_old(jc3)*ndims_new(1)*(i3-1)

    end function ind2


    !> Correct the filter for the interpolation if the direction is inverted
    !! after this check put the absolute value in scal
    pure subroutine redirect(m_isf,ms,me,shf,scal)
      implicit none
      integer, intent(in) :: m_isf
      integer, intent(inout) :: ms,me
      real(gp), intent(inout) :: scal
      real(gp), dimension(-m_isf:m_isf), intent(inout) :: shf
      !local variables
      !!!! integer :: ishf
      !!!! real(gp) :: tt

      !do nothing otherwise
!!$      if (sign(1.0_gp,scal) < 0.0_gp) then
!!$         !swap array 
!!$         do ishf=1,m_isf
!!$            tt=shf(ishf)
!!$            shf(ishf)=shf(-ishf)
!!$            shf(-ishf)=tt
!!$         end do
!!$         !invert convolutions
!!$         ishf=me
!!$         me=-ms
!!$         ms=-ishf
!!$      end if
      
      scal=sqrt(abs(scal))
    end subroutine redirect


    pure subroutine shift_and_start(ntr,istep,i2,i3,j1,j2,j3,&
         dt,istart,ms,me,alpha)
      use module_base
      implicit none
      integer, intent(in) :: ntr !< id of the dimension to be transformed
      integer, intent(in) :: istep,i2,i3
      integer, intent(in) :: j1,j2,j3
      integer, intent(out) :: istart,ms,me
      real(gp), intent(out) :: dt
      real(gp), intent(out) :: alpha !< scaling to preserve the normalization
      !local variables
      integer :: ivars
      !!!! integer(kind=8) :: istart_shift
      real(gp), dimension(3) :: t
      real(gp) :: coord_old
      !!!! real(gp) :: tt

      !define the coordinates in the reference frame, which depends on the transformed variables
      t(1)=-centre_new(1)+real(j1-1,gp)*hgrids_new(1) !the first step is always the same
      if (istep >=2) then
         t(2)=-centre_new(2)+real(j2-1,gp)*hgrids_new(2)
      else
         t(2)=-centre_old(i2)+real(j2-1,gp)*hgrids_old(i2)
      end if
      if (istep ==3) then
         t(3)=-centre_new(3)+real(j3-1,gp)*hgrids_new(3)
      else
         t(3)=-centre_old(i3)+real(j3-1,gp)*hgrids_old(i3)
      end if

      !code for the coords
      ivars=1000*istep+100+10*i2+i3

      !define the value of the shift of the variable we are going to transform
      !coordinate that has to be found in the old box, including the shift
      call old_coord(ntr,ivars,rmat,t(1),t(2),t(3),coord_old,alpha)
      coord_old=coord_old-da(ntr)

      !scaling factor to preserve the normalization of the convolution
      alpha=alpha*hgrids_new(ntr)/hgrids_old(ntr)

      !now finds ibar and delta
      call ibar_and_delta(coord_old,centre_old(ntr),hgrids_old(ntr),ndims_old(ntr),&
           centre_new(ntr),hgrids_new(ntr),m_isf,istart,dt,ms,me)      

      !if (sign(1.0_gp,alpha) < 0.0_gp) dt=-dt
      !!!!seems that coord_old can be multiplied by  hold/(hnew*a), where a is the scaling between xnew and xold
      !!!!central point of the convolution rounded to the grid points
      !!!tt=(coord_old+centre_old(ntr)+hgrids_old(ntr))/hgrids_old(ntr)
      !!!
      !!!!after this points has been identified, it seems we might start finding istart
      !!!!and the corresponding delta, followed by the step sizes
      !!!if (tt > real(ndims_old(ntr),gp)) then
      !!!   istart=ndims_old(ntr)
      !!!else
      !!!   istart=min(max(1,nint(tt)),ndims_old(ntr))
      !!!end if
      !!!
      !!!!this shift brings the old point in the new reference frame and it defines the local translation.
      !!!dt=real(istart,gp)-(coord_old+centre_new(ntr)+hgrids_new(ntr))/hgrids_old(ntr)
      !!!
      !!!!test: defines the local translation without new reference:
      !!!!this shift brings the old point in the new reference frame and it defines the local translation.
      !!!!dt=real(istart,gp)-tt
      !!!
      !!!
      !!!!purify the shift to be lower than a multiple of the grid spacing
      !!!istart_shift=nint(dt,kind=8)
      !!!dt=dt-real(istart_shift,gp)
      !!!   !!$istart=istart-istart_shift
      !!!   !!$!identify extremes for the convolution
      !!!   !!$ms=-min(m_isf,istart-1)
      !!!   !!$me=min(m_isf,ndims_old(ntr)-istart)
      !!!
      !!!istart_shift=int(istart,kind=8)-istart_shift
      !!!!identify extremes for the convolution, here a scaling factor should be added
      !!!if (istart_shift - 1 > int(m_isf,kind=8)) then
      !!!   ms=-m_isf
      !!!else if (1 - istart_shift > int(m_isf,kind=8)) then !starting point is too far
      !!!   ms=m_isf+1 !loop is disabled
      !!!else
      !!!   ms= 1 - int(istart_shift)
      !!!end if
      !!!! check the other extreme
      !!!if (int(ndims_old(ntr),kind=8)-istart_shift > int(m_isf,kind=8)) then
      !!!   me=m_isf
      !!!else if (int(ndims_old(ntr),kind=8)-istart_shift < -int(m_isf,kind=8)) then
      !!!   me=-m_isf-1 !disable loop
      !!!else
      !!!   me=ndims_old(ntr)-int(istart_shift)
      !!!end if
      !!!!check ms <= me ==> 1-istart <= ndims_old - istart ==> -2*istart <= ndims_old -1 ==> 2*istart >= 1 - ndims_old
      !!!!reconvert to default integer value
      !!!istart=int(istart_shift)
    end subroutine shift_and_start

    !> identify the point from which the old data has to be interpolated
    !! and the shift with respect to this point
    !! identify also the start and end points of the convolution
    pure subroutine ibar_and_delta(xold,cold,hold,nold,cnew,hnew,m_isf,ibar,delta,ms,me)
      implicit none
      integer, intent(in) :: nold  !<number of points in of the old direction
      integer, intent(in) :: m_isf !<(half of the) size of the ISF support 
      real(gp), intent(in) :: xold !<previous coordinate in the old reference frame
      real(gp), intent(in) :: cold !<coordinate of the center of rotation in the old reference
      real(gp), intent(in) :: hold !<grid spacing of the old reference frame
      real(gp), intent(in) :: cnew !<coordinate of the center of rotation in the new reference
      real(gp), intent(in) :: hnew !<grid spacing of the new reference frame
      integer, intent(out) :: ibar !< point in the old frame center of the convolution 
      integer, intent(out) :: ms !< starting point of the convolution
      integer, intent(out) :: me !< ending point of the convolution
      real(gp), intent(out) :: delta !< shift of the convolution, defining the filter
      !local variables
      integer(kind=8) :: istart_shift
      real(gp) :: tt

      !seems that xold can be multiplied by  hold/(hnew*a), where a is the scaling between xnew and xold
      !central point of the convolution rounded to the grid points
      tt=(xold+cold+hold)/hold

      !after this point has been identified, we might start finding ibar
      !and the corresponding delta, followed by the step sizes
      if (tt > real(nold,gp)) then
         ibar=nold
      else
         ibar=min(max(1,nint(tt)),nold)
      end if

      !this shift brings the old point in the new reference frame and it defines the local translation.
      delta=real(ibar,gp)-(xold+cnew+hnew)/hold

      !purify the shift to be lower than a multiple of the grid spacing
      istart_shift=nint(delta,kind=8)
      delta=delta-real(istart_shift,gp)

      istart_shift=int(ibar,kind=8)-istart_shift
      !identify extremes for the convolution, here a scaling factor should be added
      if (istart_shift - 1 > int(m_isf,kind=8)) then
         ms=-m_isf
      else if (1 - istart_shift > int(m_isf,kind=8)) then !starting point is too far
         ms=m_isf+1 !loop is disabled
         delta=0.0_gp
      else
         ms= 1 - int(istart_shift)
      end if
      ! check the other extreme
      if (int(nold,kind=8)-istart_shift > int(m_isf,kind=8)) then
         me=m_isf
      else if (int(nold,kind=8)-istart_shift < -int(m_isf,kind=8)) then
         me=-m_isf-1 !disable loop
         delta=0.0_gp
      else
         me=nold-int(istart_shift)
      end if
      !check ms <= me ==> 1-ibar <= nold - ibar ==> -2*ibar <= nold -1 ==> 2*ibar >= 1 - nold
      !reconvert to default integer value. If the loop is disabled this conversion is useless
      ibar=int(istart_shift)
      
    end subroutine ibar_and_delta


    pure subroutine old_coord(icrd,ivars,rmat,x,y,z,coord,alpha)
      use module_defs, only: gp
      implicit none
      integer, intent(in) :: icrd !<id of the old coordinate to be retrieved
      integer, intent(in) :: ivars !< order of the variables in terms of 1000*istep+first*100+second*10+third
      real(gp), dimension(3,3), intent(in) :: rmat !< rotation matrix entries
      real(gp), intent(in) :: x,y,z !<coordinates to be used for the mapping
      real(gp), intent(out) :: coord
      real(gp), intent(out) :: alpha !<scaling between the new coordinate and the old one xold=alpha*xnew, xnew=x(istep),xold=x(icrd)
      !local variables


      !the possible cases are
      !123 (1,2,3)
      !122 (1
      !113 (2
      !121 (2
      !112 (3)

      coord=0.0_gp
      alpha=1.0_gp
      select case(icrd)
      case(1) !x coordinate
         select case(ivars)
         case(1123)!'xnyozo')
            coord=(x-rmat(2,1)*y-rmat(3,1)*z)/rmat(1,1)
            alpha=1.0_gp/rmat(1,1)
         case(2123)!'xnynzo')
            coord=(rmat(2,2)*x-rmat(2,1)*y+rmat(1,3)*z)/rmat(3,3)
            alpha=-rmat(2,1)/rmat(3,3)
         case(3123)!'xnynzn')
            coord=rmat(1,1)*x + rmat(1,2)*y  + rmat(1,3)*z
            alpha=rmat(1,3)
         case(2122)!'xnynyo')
            coord=(-rmat(3,2)*x + rmat(3,1)*y + rmat(1,3)*z)/rmat(2,3)
            alpha=rmat(3,1)/rmat(2,3)
         end select
      case(2) !y coordinate
         select case(ivars)
         case(1113)!'xnxozo')
            coord=(x-rmat(1,1)*y-rmat(3,1)*z)/rmat(2,1)
            alpha=1.0_gp/rmat(2,1)
         case(2121)!'xnynxo')
            coord=(rmat(3,2)*x-rmat(3,1)*y +rmat(2,3)*z)/rmat(1,3)
            alpha=-rmat(3,1)/rmat(1,3)
         case(2123)!'xnynzo')
            coord=(-rmat(1,2)*x +rmat(1,1)*y + rmat(2,3)*z )/rmat(3,3)
            alpha=rmat(1,1)/rmat(3,3)
         case(3123)!'xnynzn')
            coord=rmat(2,1)*x +rmat(2,2)*y +rmat(2,3)*z
            alpha=rmat(2,3)
         end select
      case(3) !z coordinate
         select case(ivars)
         case(1112)!'xnxoyo')
            coord=(x-rmat(1,1)*y - rmat(2,1)*z)/rmat(3,1)
            alpha=1.0_gp/rmat(3,1)
         case(2121)!'xnynxo')
            coord=(rmat(3,3)*z - rmat(2,2)*x + rmat(2,1)*y)/rmat(1,3)
            alpha=rmat(2,1)/rmat(1,3)
         case(2122)!'xnynyo')
            coord=(rmat(1,2)*x+rmat(3,3)*z-rmat(1,1)*y)/rmat(2,3)
            alpha=-rmat(1,1)/rmat(2,3)
         case(3123)!'xnynzn')
            coord=rmat(3,3)*z + rmat(3,2)*y + rmat(3,1)*x
            alpha=rmat(3,3)
         end select
      end select
    end subroutine old_coord

    pure function convolve(idim,i,j,k,ms,me,m_isf,shf,n1,n2,n3,f_in,alpha)
      use module_defs, only: gp
      implicit none
      integer, intent(in) :: idim !<dimension to be convolved
      integer, intent(in) :: n1,n2,n3,m_isf
      integer, intent(in) :: i,j,k !< starting point of the convolution
      integer, intent(in) :: ms,me !< extremes for the shift
      real(gp), dimension(-m_isf:m_isf), intent(in) :: shf
      real(gp), dimension(n1,n2,n3), intent(in) :: f_in
      real(gp), intent(in) :: alpha !<scaling factor to be applied after convolution
      real(gp) :: convolve
      !local variables
      integer :: l
      real(gp) :: tt

      tt=0.0_gp
      select case(idim)
      case(1)
         do l=ms,me
            tt=tt+shf(l)*f_in(i+l,j,k)
         end do
      case(2)
         do l=ms,me
            tt=tt+shf(l)*f_in(i,j+l,k)
         end do
      case(3)
         do l=ms,me
            tt=tt+shf(l)*f_in(i,j,k+l)
         end do
      end select

      !end of interpolate coefficient
      convolve=alpha*tt

    end function convolve

    pure subroutine define_filter(dt,nrange,nphi,phi,shf)
      use module_base
      implicit none
      integer, intent(in) :: nphi !< number of sampling points of the ISF function (multiple of nrange)
      integer, intent(in) :: nrange !< extension of the ISF domain in dimensionless units (even number)
      real(gp), intent(in) :: dt
      real(gp), dimension(nphi), intent(in) :: phi !< interpolating scaling function array
      real(gp), dimension(-nrange/2:nrange/2), intent(out) :: shf !< interpolating filter to be applied
      !local variables
      integer :: nunit,ish,ipos,m_isf,l,jisf

      m_isf=nrange/2
      !number of points for a unit displacement
      nunit=nphi/nrange 

      !evaluate the shift
      ish=nint(real(nunit,gp)*dt)

      !if (ish /= 0) print *,'dt',dt,ish

      !starting point in the filter definition
      ipos=ish+1
      if (ish<= 0) then
         jisf=-(abs(ish))/nunit-1
      else if (ish > 0) then
         jisf=ish/nunit+1
      else
         jisf=0
      end if
      jisf=jisf-m_isf

      !fill the filters in its nonzero coefficients
      do l=-m_isf,m_isf
         if (jisf >= -m_isf .and. jisf <= m_isf) then
            shf(l)=phi(ipos)
         else
            shf(l)=0.0_gp
         end if
         jisf=jisf+1
         ipos=ipos+nunit
      end do

      !if the sign is negative the filter can be inverted

    end subroutine define_filter


  end subroutine field_rototranslation3D


  !> Backward wavelet transform
  !! gives the anti-correlation
  subroutine back_trans_16_reversed(nd,nt,x,y)
    implicit none
    !Arguments
    integer, intent(in) :: nd !< Length of data set
    integer, intent(in) :: nt !< Length of data in data set to be transformed
    real(kind=8), intent(in) :: x(0:nd-1)  !< Input data
    real(kind=8), intent(out) :: y(0:nd-1) !< Output data
    !Local variables
    integer :: i,j,ind

    integer, parameter :: m=18
    real(kind=8), dimension(-m:m) :: ch=(/0d0, 0d0,&
         3.571912260328699082d-6, -1.1450094552100700164d-6, &
         -0.00005642629040127758254d0, 0.00002345539585568117642d0, &
         0.0004069961892884996228d0, -0.0002465534369237166607d0, &
         -0.001634776719899382798d0, 0.00259729967896342247d0, &
         0.006477427625463336123d0, -0.01262044842878062896d0, &
         -0.02535252967734825372d0, 0.02966399618206407251d0, &
         0.06485097060728547963d0, -0.0289320622117497406d0, &
         0.0185085845718848147d0, 0.5048199552943667001d0, &
         0.970046711566057329d0, 0.7212353426722887695d0, &
         0.0294258861485558961d0, -0.2797722999367705543d0, &
         -0.0990303522418633099d0, 0.07410630821538452139d0,&
         0.04680637576666147908d0, -0.011843799423550127927d0, &
         -0.0122154536585793166d0, 0.0010521128108874154748d0, &
         0.00196569149666800115d0, -0.00008582923667387588177d0, &
         -0.0002141180336992365887d0, 3.667434093271785533d-6,&
         0.000011440737665613076119d0, 0d0, 0d0, 0d0, 0d0/)
    real(kind=8), dimension(-m:m) :: cg,cht,cgt

    !******** coefficients for wavelet transform *********************
    do i=-m,m
       cht(i)=0.d0
       cg(i)=0.d0
       cgt(i)=0.d0
    enddo

    ! the normalization is chosen such that a constant function remains the same constant 
    ! on each level of the transform

    cht( 0)=1.D0

    ! g coefficients from h coefficients
    do i=-m,m-1
       cg(i+1)=cht(-i)*(-1.d0)**(i+1)
       cgt(i+1)=ch(-i)*(-1.d0)**(i+1)
    enddo


    do i=0,nt/2-1
       y(2*i+0)=0.d0
       y(2*i+1)=0.d0

       do j=-m/2,m/2-1

          ! periodically wrap index if necessary
          ind=i-j
          loop99: do
             if (ind.lt.0) then 
                ind=ind+nt/2
                cycle loop99
             end if
             if (ind.ge.nt/2) then 
                ind=ind-nt/2
                cycle loop99
             end if
             exit loop99
          end do loop99

          y(2*i+0)=y(2*i+0) + ch(2*j-0)*x(ind)+cg(2*j-0)*x(ind+nt/2)
          y(2*i+1)=y(2*i+1) + ch(2*j+1)*x(ind)+cg(2*j+1)*x(ind+nt/2)
       end do
    end do

  END SUBROUTINE back_trans_16_reversed

  subroutine my_scaling_function4b2B(itype,nd,nrange,a,x)
    use module_base
    implicit none
    !Arguments
    !Type of interpolating functions
    integer, intent(in) :: itype
    !Number of points: must be 2**nex
    integer, intent(in) :: nd
    integer, intent(out) :: nrange
    real(kind=8), dimension(0:nd), intent(out) :: a
    real(kind=8), dimension(0:nd,2), intent(out) :: x
    !Local variables
    character(len=*), parameter :: subname='scaling_function4b2B'
    real(kind=8), dimension(:), allocatable :: y
    integer :: i,nt,ni

    call f_routine(id=subname)

    !Only itype=8,14,16,20,24,30,40,50,60,100
    select case(itype)
    case(16)
       !O.K.
    case default
       print *,"Only interpolating functions 8, 14, 16, 20, 24, 30, 40, 50, 60, 100"
       stop
    end select
!!$  write(unit=*,fmt="(1x,a,i0,a)") &
!!$       "Use interpolating scaling functions of ",itype," order"

    !Give the range of the scaling function
    !from -itype to itype
    ni=2*itype
    nrange = ni

    y = f_malloc(0.to.nd,id='y')

    ! plot scaling function
    call zero(nd+1,x(0,1))
    call zero(nd+1,y)
    nt=ni
    x(nt/2,1)=1.d0
    loop1: do
       nt=2*nt
       call back_trans_16(nd,nt,x(0,1),y)
       do i=0,nt-1
          x(i,1)=y(i)
       end do
       if (nt.eq.nd) then
          exit loop1
       end if
    end do loop1

    ! plot reversed scaling function
    call zero(nd+1,x(0,2))
    call zero(nd+1,y)
    nt=ni
    x(nt/2,2)=1.d0
    loop2: do
       nt=2*nt
       call back_trans_16_reversed(nd,nt,x(0,2),y)
       do i=0,nt-1
          x(i,2)=y(i)
       end do
       if (nt.eq.nd) then
          exit loop2
       end if
    end do loop2


    !open (unit=1,file='scfunction',status='unknown')
    do i=0,nd
       a(i) = real(i*ni,kind=8)/real(nd,kind=8)-(.5d0*real(ni,kind=8)-1.d0)
       !write(1,*) a(i),x(i)
    end do
    !close(1)

    call f_free(y)
    call f_release_routine()
  END SUBROUTINE my_scaling_function4b2B

  
end module reformatting
