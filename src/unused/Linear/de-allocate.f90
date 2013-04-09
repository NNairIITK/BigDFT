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







subroutine allocateCommunicationbufferSumrho(iproc, comsr, subname)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc
type(p2pComms),intent(inout):: comsr
character(len=*),intent(in):: subname

! Local variables
integer:: istat
call timing(iproc,'allocommsumrho','ON') !lr408t
allocate(comsr%sendBuf(comsr%nsendBuf), stat=istat)
call memocc(istat, comsr%sendBuf, 'comsr%sendBuf', subname)
call razero(comsr%nSendBuf, comsr%sendBuf)

allocate(comsr%recvBuf(comsr%nrecvBuf), stat=istat)
call memocc(istat, comsr%recvBuf, 'comsr%recvBuf', subname)
call razero(comsr%nrecvBuf, comsr%recvBuf)
call timing(iproc,'allocommsumrho','OF') !lr408t
end subroutine allocateCommunicationbufferSumrho


subroutine deallocateCommunicationbufferSumrho(comsr, subname)
use module_base
use module_types
implicit none

! Calling arguments
!type(p2pCommsSumrho),intent(inout):: comsr
type(p2pComms),intent(inout):: comsr
character(len=*),intent(in):: subname

! Local variables
integer:: istat, iall

iall=-product(shape(comsr%sendBuf))*kind(comsr%sendBuf)
deallocate(comsr%sendBuf, stat=istat)
call memocc(istat, iall, 'comsr%sendBuf', subname)

iall=-product(shape(comsr%recvBuf))*kind(comsr%recvBuf)
deallocate(comsr%recvBuf, stat=istat)
call memocc(istat, iall, 'comsr%recvBuf', subname)

end subroutine deallocateCommunicationbufferSumrho
