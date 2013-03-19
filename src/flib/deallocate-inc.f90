!include file for deallocation template

  !local variables
  integer :: istat,ierr
  character(len=info_length) :: address
  !local variables
  integer :: i_all
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
