!> @file
!!    Modulefile for handling the read-write of a given simulation 
!!    box
!!
!! @author
!!    G. Fisicaro, L. Genovese (September 2015)
!!    Copyright (C) 2002-2015 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
module IObox
  use PSbase
  implicit none

  integer, parameter :: UNKNOWN=0
  integer, parameter :: CUBE=1
  integer, parameter :: ETSF=2
  
  private

  public :: read_field,read_field_dimensions

  contains

    pure subroutine startend_buffers(geocode,nl1,nl2,nl3,nbx,nby,nbz,nc1,nc2,nc3)
      implicit none
      character(len=1), intent(in) :: geocode
      integer, intent(out) :: nl1,nl2,nl3,nbx,nby,nbz
      integer, intent(in), optional :: nc1,nc2,nc3

      if (geocode /= 'F') then
         nl1=1
         nl3=1
         nbx = 1
         nbz = 1
      else
         nl1=15
         nl3=15
         nbx = 0
         nbz = 0
      end if
      !value of the buffer in the y direction
      if (geocode == 'P') then
         nl2=1
         nby = 1
      else
         nl2=15
         nby = 0
      end if
     
    end subroutine startend_buffers

    function get_file_format(filename,isuffix) result(fformat)
      implicit none
      character(len=*), intent(in) :: filename
      integer, intent(out) :: isuffix
      integer :: fformat


      ! Format = 1 -> cube (default)
      ! Format = 2 -> ETSF
      ! ...
      fformat = UNKNOWN
      isuffix = index(filename, ".cube", back = .true.)
      if (isuffix > 0) then
         isuffix = isuffix - 1
         fformat = CUBE
      else
         isuffix = index(filename, ".etsf", back = .true.)
         if (isuffix <= 0) isuffix = index(filename, ".etsf.nc", back = .true.)
         if (isuffix > 0) then
            isuffix = isuffix - 1
            fformat = ETSF
         else
            isuffix = len(filename)
            fformat=UNKNOWN
         end if
      end if
      
    end function get_file_format

    subroutine read_cube_header(filename,geocode,ndims,hgrids,&
         nat,rxyz, iatypes, znucl)
      use f_utils
      use dynamic_memory
      implicit none
      character(len=*), intent(in) :: filename
      character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
      integer, dimension(3), intent(out) :: ndims
      real(gp), dimension(3), intent(out) :: hgrids
      integer, intent(out) :: nat
      real(gp), dimension(:,:), pointer :: rxyz
      integer, dimension(:), pointer :: iatypes, znucl
      !local variables
      character(len=*), parameter :: subname='read_cube_header'
      integer :: n1t,n2t,n3t,n1,n2,n3,idum,iat,j,unt
      integer :: nl1,nl2,nl3,nbx,nby,nbz,n1i,n2i,n3i
      real(gp) :: dum1,dum2,dum3,hxh,hyh,hzh
      integer, dimension(:), allocatable :: znucl_

      call startend_buffers(geocode,nl1,nl2,nl3,nbx,nby,nbz)
      unt=f_get_free_unit(22)
      call f_open_file(unit=unt,file=trim(filename)//".cube",status='old')
      read(unt,*)! 'CUBE file for charge density'
      read(unt,*)! 'Case for '//trim(message)

      read(unt,'(i5,3(f12.6),a)') nat, dum1, dum2, dum3 
      read(unt,'(i5,3(f12.6))') n1t , hxh  , dum1 , dum2
      read(unt,'(i5,3(f12.6))') n2t , dum1 , hyh  , dum2
      read(unt,'(i5,3(f12.6))') n3t , dum1 , dum2 , hzh

      !grid positions
      n1=n1t/2-nbx
      n1i=2*n1+(1-nbx)+2*nl1
      n2=n2t/2-nby
      n2i=2*n2+(1-nby)+2*nl2
      n3=n3t/2-nbz
      n3i=2*n3+(1-nbz)+2*nl3

      !atomic positions
      rxyz = f_malloc_ptr((/ 3, nat /),id='rxyz')
      iatypes = f_malloc_ptr(nat,id='iatypes')
      znucl_ = f_malloc(nat,id='znucl_')
      znucl_(:) = -1

      do iat=1,nat
         read(unt,'(i5,4(f12.6))') idum , dum1 , (rxyz(j,iat),j=1,3)
         do j = 1, nat, 1
            if (znucl_(j) == idum .or. znucl_(j) == -1) then
               znucl_(j) = idum
               exit
            end if
         end do
         iatypes(iat) = j
      end do

      do j = 1, nat, 1
         if (znucl_(j) == -1) then
            exit
         end if
      end do
      znucl = f_malloc_ptr(j-1,id='znucl')
      znucl(1:j-1) = znucl_(1:j-1)

      call f_free(znucl_)

      call f_close(unt)

      ndims(1)=n1i
      ndims(2)=n2i
      ndims(3)=n3i
      hgrids(1)=hxh
      hgrids(2)=hyh
      hgrids(3)=hzh

    END SUBROUTINE read_cube_header

    !>   Read a cube field which have been plotted previously by write_cube_fields
    subroutine read_cube_field(filename,geocode,ndims,rho)
      use f_utils
      use dynamic_memory
      implicit none
      character(len=*), intent(in) :: filename
      character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
      integer, dimension(3), intent(in) :: ndims
      real(dp), dimension(ndims(1),ndims(2),ndims(3)) :: rho
      !local variables
      !n(c) character(len=*), parameter :: subname='read_cube_field'
      character(len=3) :: advancestring
      integer :: n1t,n2t,n3t,n1,n2,n3,i1,i2,i3,nat,iat,unt
      integer :: nl1,nl2,nl3,nbx,nby,nbz,icount,ind,n1i,n2i,n3i
      real(gp) :: dum1,dum2,dum3,tt

      call startend_buffers(geocode,nl1,nl2,nl3,nbx,nby,nbz)

      !aliasing
      n1i=ndims(1)
      n2i=ndims(2)
      n3i=ndims(3)

      unt=f_get_free_unit(22)
      call f_open_file(unit=unt,file=trim(filename)//'.cube',status='old')
      read(unt,*)! 'CUBE file for charge density'
      read(unt,*)! 'Case for '//trim(message)

      read(unt,'(i5,3(f12.6),a)')  nat , dum1, dum2, dum3 
      read(unt,'(i5,3(f12.6))') n1t , dum3,   dum1 ,dum2
      read(unt,'(i5,3(f12.6))') n2t ,dum1 , dum3  ,  dum2
      read(unt,'(i5,3(f12.6))') n3t ,dum1 , dum2 , dum3

      !grid positions
      n1=n1t/2-nbx
      if (n1i /= 2*n1+(1-nbx)+2*nl1) stop 'n1i not valid'
      n2=n2t/2-nby
      if (n2i /= 2*n2+(1-nby)+2*nl2) stop 'n2i not valid'
      n3=n3t/2-nbz
      if (n3i /= 2*n3+(1-nbz)+2*nl3) stop 'n3i not valid'

      !zero the buffer
      call f_zero(rho)

      do iat=1,nat
         !read(unt,'(i5,4(f12.6))')! idum , dum1 , (rxyz(j,iat),j=1,3)
         read(unt,*)! idum , dum1 , (rxyz(j,iat),j=1,3)
      end do

      !the loop is reverted for a cube file
      !charge normalised to the total charge
      do i1=0,2*(n1+nbx) - 1
         do i2=0,2*(n2+nby) - 1
            icount=0
            do i3=0,2*(n3+nbz) - 1
               icount=icount+1
               if (icount == 6 .or. i3==2*(n3+nbz) - 1) then
                  advancestring='yes'
                  icount=0
               else
                  advancestring='no'
               end if
               ind=i1+nl1+(i2+nl2-1)*n1i+(i3+nl3-1)*n1i*n2i
               read(unt,'(1x,1pe13.6)',advance=advancestring) tt !rho(ind) 
               !rho(ind)=tt
               rho(i1+nl1,i2+nl2,i3+nl3)=tt
               !           write(16,*)i1,i2,i3,ind,rho(ind)
               !read(unt,*)',advance=advancestring) rho(ind) 
            end do
         end do
      end do
      call f_close(unt)

      ! write(14,*)rho

    END SUBROUTINE read_cube_field

    !> routine to be used to estimate the dimension of the array to be allocated
    subroutine read_field_dimensions(filename,geocode,ndims,nspin)
      use dictionaries, only: f_err_throw
      use dynamic_memory
      implicit none
      character(len=*), intent(in) :: filename
      character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
      integer, intent(out) :: nspin
      integer, dimension(3), intent(out) ::  ndims
      !local variables
      integer :: nat,fformat,isuffix
      real(gp), dimension(1) :: rho_fake
      real(gp), dimension(3) :: hgrids
      real(gp), dimension(:,:), pointer   :: rxyz2
      integer, dimension(:), pointer   :: iatypes2, znucl2

      fformat=get_file_format(filename,isuffix)

      select case(fformat)
      case(CUBE)
         call read_cube(filename(1:isuffix),geocode,ndims,hgrids,nspin,1,1,rho_fake,&
              nat,rxyz2, iatypes2, znucl2,dry_run=.true.)
         call f_free_ptr(rxyz2)
         call f_free_ptr(iatypes2)
         call f_free_ptr(znucl2)
      case(ETSF)
         call f_err_throw('Size estimator not (yet) possible for ETSF format')
      end select

    end subroutine read_field_dimensions

    !> Read a density file using file format depending on the extension.
    subroutine read_field(filename,geocode,ndims,hgrids,nspin,ldrho,nrho,rho,&
         nat,rxyz,iatypes,znucl)
      use dynamic_memory
      use dictionaries, only: f_err_throw
      implicit none
      character(len=*), intent(in) :: filename
      character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
      integer, intent(in) :: ldrho,nrho
      integer, intent(out) :: nspin
      integer, dimension(3), intent(out) ::  ndims
      real(gp), dimension(3), intent(out) :: hgrids
      real(dp), dimension(ldrho,nrho), intent(inout) :: rho
      real(gp), dimension(:,:), pointer, optional :: rxyz
      integer, intent(out), optional ::  nat
      integer, dimension(:), pointer, optional :: iatypes, znucl
      !local variables
      character(len = *), parameter :: subname = "read_field"
      logical, dimension(4) :: optargs
      integer :: isuffix,fformat,nat_read
      real(gp), dimension(:,:), pointer :: rxyz_read
      integer, dimension(:), pointer :: iatypes_read, znucl_read

      optargs=[present(rxyz),present(nat),present(iatypes),present(znucl)]
      !check the arguments
      if ((.not. all(optargs)) .and. (any(optargs))) &
           call f_err_throw('Wrong usage of read_field, rxyz, znucl and iatypes should be _all_ present')

      call f_routine(id=subname)

      fformat=get_file_format(filename,isuffix)

      select case(fformat)
      case(CUBE)
         call read_cube(filename(1:isuffix),geocode,ndims,hgrids,nspin,ldrho,nrho,rho,&
              nat_read,rxyz_read, iatypes_read, znucl_read)
      case(ETSF)
         call read_etsf(filename(1:isuffix),geocode,&
              ndims(1),ndims(2),ndims(3),nspin,hgrids(1),hgrids(2),hgrids(3),rho,&
              nat_read,rxyz_read, iatypes_read, znucl_read)
         if (ldrho < product(ndims) .or. nrho < nspin) &
              call f_err_throw('Severe error, the sizes of the rho array have revealed not to be sufficient. '//&
              'The etsf reading  might have caused a boundary error')
      end select

      if (all(optargs)) then
         rxyz => rxyz_read
         iatypes => iatypes_read
         znucl => znucl_read
         nat=nat_read
      else
         call f_free_ptr(rxyz_read)
         call f_free_ptr(iatypes_read)
         call f_free_ptr(znucl_read)
      end if
      call f_release_routine()
    END SUBROUTINE read_field

    !> Read density or potential in cube format
    subroutine read_cube(filename,geocode,ndims,hgrids,nspin,ldrho,nrho,rho,&
         nat,rxyz, iatypes, znucl,dry_run)
      use PSbase
      use f_utils
      use yaml_strings, only: operator(+),f_strcpy
      use dictionaries, only: f_err_throw
      use dynamic_memory
      implicit none
      character(len=*), intent(in) :: filename
      character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
      integer, intent(in) :: ldrho,nrho !<dimensions of the rho array
      integer, intent(out) :: nspin
      integer, dimension(3), intent(out) ::  ndims
      real(gp), dimension(3), intent(out) :: hgrids
      real(dp), dimension(ldrho,nrho), intent(inout) :: rho
      real(gp), dimension(:,:), pointer   :: rxyz
      integer, intent(out)   ::  nat
      integer, dimension(:), pointer   :: iatypes, znucl
      logical, intent(in), optional :: dry_run !<only retrieve dimensions, do not fill the array
      !local variables
      !n(c) character(len=*), parameter :: subname='read_cube'
      character(len=5) :: suffix
      character(len=15) :: message
      integer, dimension(3) :: na,nb
      real(gp), dimension(3) :: ha,hb
      integer :: ia,nat2
      logical :: exists,drr
      real(gp), dimension(:,:), pointer   :: rxyz2
      integer, dimension(:), pointer   :: iatypes2, znucl2

      ! Test if we have up and down densities.
      call f_file_exists(file=trim(filename)//"-up.cube",exists=exists)
      if (exists) then
         call f_file_exists(file=trim(filename)//"-down.cube",exists=exists)
         if (.not.exists) then
            call f_err_throw("found a "+filename+"-up.cube file, but no -down.cube...")
            nspin = 1
         else
            nspin = 2
         end if
      else
         call f_file_exists(file=trim(filename)//".cube",exists=exists)
         if (.not. exists) call f_err_throw('The file '+filename+' does not exists')
         nspin = 1
      end if

      drr=.false.
      if (present(dry_run)) drr=dry_run
      
      !read the header of the files and verify it is coherent for nspin==2
      if (nspin /= 2) then
         call f_strcpy(src='',dest=suffix)
         call read_cube_header(filename//trim(suffix),geocode,ndims,hgrids,&
              nat,rxyz, iatypes, znucl)
      else
         call f_strcpy(src='-up',dest=suffix)
         call read_cube_header(filename//trim(suffix),geocode,na,ha,&
              nat,rxyz, iatypes, znucl)
         call f_strcpy(src='-down',dest=suffix)
         call read_cube_header(filename//trim(suffix),geocode,nb,hb,&
              nat2,rxyz2, iatypes2, znucl2)
         if (any(na /= nb .or. ha /= hb)) then
            call f_err_throw('Error in reading .cube file, the dimensions between spin up and down are not coherent')
         else
            ndims=na
            hgrids=ha
         end if
         call f_free_ptr(rxyz2)
         call f_free_ptr(iatypes2)
         call f_free_ptr(znucl2)
      end if

      !now the information to allocate the array has been set up
      if (drr) then
         return
      else
         !check if the given dimension for rho is compatible with the array
         if (ldrho < product(ndims) .or. nrho /= nspin) &
              call f_err_throw('The dimension of the rho array is not coherent with the file')
      end if

      !read the header of the files and verify it is coherent for nspin==2
      if (nspin /= 2) then
         call f_strcpy(src='',dest=suffix)
         call read_cube_field(filename//trim(suffix),geocode,ndims,rho)
      else
         call f_strcpy(src='-up',dest=suffix)
         call read_cube_field(filename//trim(suffix),geocode,ndims,rho(1,1))
         call f_strcpy(src='-down',dest=suffix)
         call read_cube_field(filename//trim(suffix),geocode,ndims,rho(1,2))
      end if

    END SUBROUTINE read_cube

end module IObox
