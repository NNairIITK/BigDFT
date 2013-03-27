!> @file
!! Include fortran file for deallocation template
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


  !local variables
  integer :: ierr
  character(len=info_length) :: address
  !local variables
  integer(kind=8) :: ilsize
  !      call timing(0,'AllocationProf','IR') 
  !profile the array allocation
  call getaddress(array,address,len(address),ierr)
  !here the size should be corrected with ndebug
  ilsize=int(product(shape(array))*kind(array),kind=8)
  !fortran deallocation
  deallocate(array,stat=ierror)
  !hopefully only address is necessary for the deallocation
  !in the include file this part can be expanded if needed
  call profile_deallocation(ierror,ilsize,address) 
  !      call timing(0,'AllocationProf','RS') 
