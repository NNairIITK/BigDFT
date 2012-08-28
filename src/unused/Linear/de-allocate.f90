subroutine deallocateRecvBufferOrtho(comon, subname)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(p2pComms),intent(inout):: comon
  character(len=*),intent(in):: subname

  ! Local variables
  integer:: istat, iall

  iall = -product(shape(comon%recvBuf))*kind(comon%recvBuf)
  deallocate(comon%recvBuf, stat=istat)
  call memocc(istat, iall, 'comon%recvBuf',subname)

end subroutine deallocateRecvBufferOrtho


subroutine allocateRecvBufferOrtho(comon, subname)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(p2pComms),intent(inout):: comon
  character(len=*),intent(in):: subname

  ! Local variables
  integer:: istat

  allocate(comon%recvBuf(comon%nrecvBuf), stat=istat)
  call memocc(istat, comon%recvBuf, 'comon%recvBuf', subname)

end subroutine allocateRecvBufferOrtho


subroutine deallocateSendBufferOrtho(comon, subname)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(p2pComms),intent(inout):: comon
  character(len=*),intent(in):: subname

  ! Local variables
  integer:: istat, iall

  iall = -product(shape(comon%sendBuf))*kind(comon%sendBuf)
  deallocate(comon%sendBuf, stat=istat)
  call memocc(istat, iall, 'comon%sendBuf',subname)

end subroutine deallocateSendBufferOrtho


subroutine allocateSendBufferOrtho(comon, subname)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(p2pComms),intent(inout):: comon
  character(len=*),intent(in):: subname

  ! Local variables
  integer:: istat

  allocate(comon%sendBuf(comon%nsendBuf), stat=istat)
  call memocc(istat, comon%sendBuf, 'comon%sendBuf', subname)

end subroutine allocateSendBufferOrtho

subroutine allocateCommuncationBuffersOrtho(comon, subname)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(p2pComms),intent(inout):: comon
  character(len=*),intent(in):: subname

  ! Local variables
  integer:: istat

  allocate(comon%recvBuf(comon%nrecvBuf), stat=istat)
  call memocc(istat, comon%recvBuf, 'comon%recvBuf', subname)
  allocate(comon%sendBuf(comon%nsendBuf), stat=istat)
  call memocc(istat, comon%sendBuf, 'comon%sendBuf', subname)

end subroutine allocateCommuncationBuffersOrtho



subroutine deallocateCommuncationBuffersOrtho(comon, subname)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(p2pComms),intent(inout):: comon
  character(len=*),intent(in):: subname

  ! Local variables
  integer:: istat, iall

  iall = -product(shape(comon%recvBuf))*kind(comon%recvBuf)
  deallocate(comon%recvBuf, stat=istat)
  call memocc(istat, iall, 'comon%recvBuf',subname)
  iall = -product(shape(comon%sendBuf))*kind(comon%sendBuf)
  deallocate(comon%sendBuf, stat=istat)
  call memocc(istat, iall, 'comon%sendBuf',subname)

end subroutine deallocateCommuncationBuffersOrtho


