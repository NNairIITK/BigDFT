!> @file
!>  Modules which contains information about multipole structure
!! @author
!!    Copyright (C) 2013-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

module multipole_base
  use module_base
  implicit none

  private

  integer,parameter,public :: lmax=3

  type,public :: multipole
    !integer :: l
    real(dp),dimension(:),pointer :: q
  end type multipole

  type,public :: multipole_set
    real(dp),dimension(3) :: rxyz
    type(multipole),dimension(:),pointer :: qlm
  end type multipole_set

  type,public :: external_potential_descriptors
    integer :: nmpl
    type(multipole_set),dimension(:),pointer :: mpl
  end type external_potential_descriptors

  !> Public routines
  public :: multipoles_from_dict
  public :: external_potential_descriptors_null
  public :: deallocate_external_potential_descriptors

  contains


    pure function multipole_null() result(mp)
      implicit none
      type(multipole) :: mp
      call nullify_multipole(mp)
    end function multipole_null


    pure subroutine nullify_multipole(mp)
      implicit none
      type(multipole),intent(out) :: mp
      !mp%l = 0
      nullify(mp%q)
    end subroutine nullify_multipole


    pure function multipole_set_null() result(mps)
      implicit none
      type(multipole_set) :: mps
      call nullify_multipole_set(mps)
    end function multipole_set_null


    pure subroutine nullify_multipole_set(mps)
      implicit none
      type(multipole_set),intent(out) :: mps
      mps%rxyz(1:3) = 0._dp
      !mps%lmax = 0
      nullify(mps%qlm)
    end subroutine nullify_multipole_set


    pure function external_potential_descriptors_null() result(ep)
      implicit none
      type(external_potential_descriptors) :: ep
      call nullify_external_potential_descriptors(ep)
    end function external_potential_descriptors_null


    pure subroutine nullify_external_potential_descriptors(ep)
      implicit none
      type(external_potential_descriptors),intent(out) :: ep
      ep%nmpl=0
      nullify(ep%mpl)
    end subroutine nullify_external_potential_descriptors


    subroutine deallocate_multipole(mp)
      implicit none
      type(multipole),intent(inout) :: mp
      call f_free_ptr(mp%q)
    end subroutine deallocate_multipole


    subroutine deallocate_multipole_set(mps)
      implicit none
      type(multipole_set),intent(inout) :: mps
      integer :: l
      do l=0,lmax
          call deallocate_multipole(mps%qlm(l))
      end do
      deallocate(mps%qlm)
    end subroutine deallocate_multipole_set


    subroutine deallocate_external_potential_descriptors(ep)
      implicit none
      type(external_potential_descriptors),intent(inout) :: ep
      integer :: impl
      do impl=1,ep%nmpl
          call deallocate_multipole_set(ep%mpl(impl))
      end do
      deallocate(ep%mpl)
    end subroutine deallocate_external_potential_descriptors


    subroutine multipoles_from_dict(dict, ep)
      use yaml_output
      implicit none

      ! Calling arguments
      type(dictionary),pointer :: dict
      type(external_potential_descriptors),intent(out) :: ep

      ! Local variables
      integer :: l, ilen, impl
      type(dictionary),pointer :: iter
      character(len=2) :: key
      real(kind=8),dimension(3) :: rxyz
      !integer,parameter :: lmax=3
      real(kind=8),dimension(2*lmax+1) :: mp_tmp
      real(dp) :: norm, dnrm2


      ! Get the number of multipole centers
      !call yaml_map('Size of the mp list',nmpl)
      !call yaml_sequence_open('Values')
      ep = external_potential_descriptors_null()
      ep%nmpl = dict_len(dict)
      allocate(ep%mpl(ep%nmpl))
      do impl=1,ep%nmpl
          ep%mpl(impl) = multipole_set_null()
          allocate(ep%mpl(impl)%qlm(0:lmax))
          do l=0,lmax
              ep%mpl(impl)%qlm(l) = multipole_null()
          end do
      end do

      ! Get the value from the dict
      call f_zero(mp_tmp)
      do impl=0,ep%nmpl-1
          !call yaml_sequence()
          !call yaml_map('Size of element'//trim(yaml_toa(impl)),dict_size(dict//impl))
          iter => dict//impl
          ! Retrieve the position of the multipole center, which always have to be given.
          if ('r' .notin. iter) then
              call f_err_throw('No positions were provided for the multipole center no. '//trim(yaml_toa(impl+1)))
          end if
          ilen = dict_len(iter//'r')
          if (ilen/=3) then
              call f_err_throw('For the multipole center no. '//trim(yaml_toa(impl+1))//&
                   &' the number of coordinates specified are wrong ('//trim(yaml_toa(ilen))//')')
          end if
          ep%mpl(impl+1)%rxyz = iter//'r'
          !call yaml_map('rxyz',ep%mpl(impl)%rxyz)
          do l=0,lmax
              key='q'//trim(adjustl(yaml_toa(l)))
              if (key .in. iter) then
                  ilen = dict_len(iter//key)
                  if (ilen/=2*l+1) then
                      call f_err_throw('Wrong len ('//trim(yaml_toa(ilen))//') of the mutipole no.'&
                           &//trim(yaml_toa(impl+1))//' (l='//trim(yaml_toa(l))//')')
                  end if
                  mp_tmp(1:2*l+1) = iter//key
                  !call yaml_map(key,mp_tmp(1:2*l+1))
              else
                  mp_tmp(1:2*l+1) = 0._dp
              end if
              !ep%mpl(impl+1)%qlm(l)%l = l
              norm = dnrm2(2*l+1, mp_tmp(1), 1)
              if (norm>0.d0) then
                  ep%mpl(impl+1)%qlm(l)%q = f_malloc_ptr(2*l+1,id='q')
                  !call yaml_map('l',mp_tmp(1:2*l+1))
                  call vcopy(2*l+1, mp_tmp(1), 1, ep%mpl(impl+1)%qlm(l)%q(1), 1)
              end if
          end do
      end do
      !call yaml_sequence_close()

    end subroutine multipoles_from_dict

end module multipole_base
