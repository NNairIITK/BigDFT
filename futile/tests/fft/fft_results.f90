!> @file
!!   Test the results for the FFT API 
!! @copyright
!!   Copyright (C) Stefan Goedecker, CEA Grenoble, 2002, Basel University, 2009
!!   Copyright (C) 2009-2013 BigDFT group
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS
program fft_results
  use futile
  use numerics
  implicit none
  character(len=*), parameter :: inputs=&
       "- {name: ndims, shortname: n, default: 36,"//&
       "  help_string: Array of dimension 3 for the points of the simulation box,"//&
       "  help_dict: {Allowed values: list of integers}}"//f_cr//&
       "- {name: sizes, shortname: s, default: 10.0,"//&
       "  help_string: Array of dimension 3 for the size of the simulation box,"//&
       "  help_dict: {Allowed values: list of floating point numbers}}"//f_cr//&
       "- {name: pref, shortname: p, default: 2,"//&
       "  help_string: Array of dimension 3 for the reference impulses,"//&
       "  help_dict: {Allowed values: list of integers}}"//f_cr//&

       "- {name: nrep, shortname: r, default: 1,"//&
       "  help_string: Number of repetitions of the bench,"//&
       "  help_dict: {Allowed values: integer}}"


  integer :: i1,i2,i3,inzee,nrep
  real(f_double) :: f,r2
  integer, dimension(3) :: n,pref
  real(f_double), dimension(3) :: L,h,sigma,r
  type(dictionary), pointer :: options
  real(f_double), dimension(:,:,:,:,:), allocatable :: zinout

  call f_lib_initialize()

  call yaml_argparse(options,inputs)
  
  n=options//'ndims'
  L=options//'sizes'
  pref=options//'pref'
  nrep=options//'nrep'

  call dict_free(options)
  zinout=f_malloc([2,n(1),n(2),n(3),2],id='zinout')

  h=L/n
  sigma=L/10.0_f_double

  inzee=1
  !fill the inout array
  do i3=1,n(3)
     r(3)=(i3-1)!*h(3)
     do i2=1,n(2)
        r(2)=(i2-1)!*h(2)
        do i1=1,n(1)
           r(1)=(i1-1)!*h(1)
           r2=sum(r**2/sigma)
           f=safe_exp(-onehalf*r2)
           zinout(1,i1,i2,i3,inzee)=product(cos(r*twopi/n*pref))
           zinout(2,i1,i2,i3,inzee)=0.0_f_double
        end do
     end do
  end do

  !test the 3d FFT
  !out-of-place case of traditional output
  call FFT_3d(n(1),n(2),n(3),n(1),n(2),n(3),zinout,1,inzee)
  
  call yaml_map('Value at pref',&
       zinout(:,j(pref(1),n(1))+1,j(pref(2),n(2))+1,j(pref(3),n(3))+1,inzee)/product(n)*8.0_f_double)
  call yaml_map('Symmetric value',&
       zinout(:,j(-pref(1),n(1))+1,j(-pref(2),n(2))+1,j(-pref(3),n(3))+1,inzee)/product(n)*8.0_f_double)
  
  call yaml_map('Other values',sum(abs(zinout(:,:,:,:,inzee)))-product(n))
    
  !then inspect results of the modification of coordinates
  do i1=1,n(1)
     call f_assert(j(p(i1,n(1)),n(1))+1 == i1,id=trim(yaml_toa(i1,fmt='(i6)')))
     !print *,'p',i1,p(i1,n(1)),n(1)+2-i1,p(n(1)+2-i1,n(1))
     call f_assert(-p(i1,n(1)) == p(n(1)+2-i1,n(1)) .or. i1==n(1)/2+1 ,id=trim(yaml_toa(n(1)+2-i1,fmt='(i6)')))
     !print *,'j',i1,j(-p(i1,n(1)),n(1))+1,n(1)+2-i1
     call f_assert( j(-p(i1,n(1)),n(1)) == n(1)+1-i1 .or. i1==1 ,id=trim(yaml_toa(n(1)+2-i1,fmt='(i6)')))


  end do

  call take_timings(nrep,n,zinout,1,inzee)

  call f_free(zinout)
  call f_lib_finalize()

  contains

    !>impulse coordinate, from 0,...,n/2+1,-n/2+1,...-1
    pure function p(i,n)
      implicit none
      integer, intent(in) :: i,n
      integer :: p
      p=i-(i/(n/2+2))*n-1
    end function p

    !>real space coordinate, from 0,...,n-1
    pure function j(p,n)
      implicit none
      integer, intent(in) :: p,n
      integer :: j
      j=p-((p+n)/n-1)*n
    end function j


end program fft_results

subroutine take_timings(nrep,n,zinout,isign,inzee) 
  use futile
  implicit none
  integer, intent(in) :: nrep,isign
  integer, dimension(3) :: n
  integer, intent(inout) :: inzee
  real(f_double), dimension(2,n(1),n(2),n(3),2), intent(inout) :: zinout
  !local variables
  integer :: irep
  integer(f_long) :: t0,t1
  !then take timing for the number of repetitions
  t0=f_time()
  do irep=1,nrep
     !test the 3d FFT
     !out-of-place case of traditional output
     call FFT_3d(n(1),n(2),n(3),n(1),n(2),n(3),zinout,1,inzee)
  end do
  t1=f_time()
  call yaml_map('No. of repetitions',nrep)
  call yaml_map('Time per repetition per element (ns)',&
       1.e-3_f_double*(t1-t0)/nrep/product(n),fmt='(f12.7)')

end subroutine take_timings
