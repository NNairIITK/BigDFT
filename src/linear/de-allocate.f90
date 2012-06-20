subroutine allocateCommunicationsBuffersPotential(comgp, subname)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  !type(p2pCommsGatherPot),intent(inout):: comgp
  type(p2pComms),intent(inout):: comgp
  character(len=*),intent(in):: subname
  
  ! Local variables
  integer:: istat
  
  allocate(comgp%recvBuf(comgp%nrecvBuf), stat=istat)
  call memocc(istat, comgp%recvBuf, 'comgp%recvBuf', subname)

end subroutine allocateCommunicationsBuffersPotential



subroutine deallocateCommunicationsBuffersPotential(comgp, subname)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  !type(p2pCommsGatherPot),intent(inout):: comgp
  type(p2pComms),intent(inout):: comgp
  character(len=*),intent(in):: subname
  
  ! Local variables
  integer:: istat, iall
  
  iall=-product(shape(comgp%recvBuf))*kind(comgp%recvBuf)
  deallocate(comgp%recvBuf, stat=istat)
  call memocc(istat, iall, 'comgp%recvBuf', subname)

end subroutine deallocateCommunicationsBuffersPotential



subroutine allocate_workarrays_quartic_convolutions(lr, subname, work)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  type(locreg_descriptors),intent(in):: lr
  character(len=*),intent(in):: subname
  type(workarrays_quartic_convolutions),intent(out):: work
  
  ! Local variables
  integer:: istat
  
  allocate(work%xx_c(0:lr%d%n1,0:lr%d%n2,0:lr%d%n3), stat=istat)
  call memocc(istat, work%xx_c, 'work%xx_c', subname)
  allocate(work%xy_c(0:lr%d%n2,0:lr%d%n1,0:lr%d%n3), stat=istat)
  call memocc(istat, work%xy_c, 'work%xy_c', subname)
  allocate(work%xz_c(0:lr%d%n3,0:lr%d%n1,0:lr%d%n2), stat=istat)
  call memocc(istat, work%xz_c, 'work%xz_c', subname)
  !!work%xx_c=0.d0
  !!work%xy_c=0.d0
  !!work%xz_c=0.d0
  call to_zero((lr%d%n1+1)*(lr%d%n2+1)*(lr%d%n3+1), work%xx_c(0,0,0))
  call to_zero((lr%d%n1+1)*(lr%d%n2+1)*(lr%d%n3+1), work%xy_c(0,0,0))
  call to_zero((lr%d%n1+1)*(lr%d%n2+1)*(lr%d%n3+1), work%xz_c(0,0,0))
  
  
  allocate(work%xx_f1(lr%d%nfl1:lr%d%nfu1,lr%d%nfl2:lr%d%nfu2,lr%d%nfl3:lr%d%nfu3), stat=istat)
  call memocc(istat, work%xx_f1, 'work%xx_f1', subname)
  !!work%xx_f1=0.d0
  call to_zero((lr%d%nfu1-lr%d%nfl1+1)*(lr%d%nfu2-lr%d%nfl2+1)*(lr%d%nfu3-lr%d%nfl3+1), &
       work%xx_f1(lr%d%nfl1,lr%d%nfl2,lr%d%nfl3))
  allocate(work%xx_f(7,lr%d%nfl1:lr%d%nfu1,lr%d%nfl2:lr%d%nfu2,lr%d%nfl3:lr%d%nfu3), stat=istat)
  call memocc(istat, work%xx_f, 'work%xx_f', subname)
  !!work%xx_f=0.d0
  call to_zero(7*(lr%d%nfu1-lr%d%nfl1+1)*(lr%d%nfu2-lr%d%nfl2+1)*(lr%d%nfu3-lr%d%nfl3+1), &
       work%xx_f(1,lr%d%nfl1,lr%d%nfl2,lr%d%nfl3))
  
  
  allocate(work%xy_f2(lr%d%nfl2:lr%d%nfu2,lr%d%nfl1:lr%d%nfu1,lr%d%nfl3:lr%d%nfu3), stat=istat)
  call memocc(istat, work%xy_f2, 'work%xy_f2', subname)
  !!work%xy_f2=0.d0
  call to_zero((lr%d%nfu1-lr%d%nfl1+1)*(lr%d%nfu2-lr%d%nfl2+1)*(lr%d%nfu3-lr%d%nfl3+1), &
       work%xy_f2(lr%d%nfl2,lr%d%nfl1,lr%d%nfl3))
  allocate(work%xy_f(7,lr%d%nfl2:lr%d%nfu2,lr%d%nfl1:lr%d%nfu1,lr%d%nfl3:lr%d%nfu3), stat=istat)
  call memocc(istat, work%xy_f, 'work%xy_f', subname)
  !!work%xy_f=0.d0
  call to_zero(7*(lr%d%nfu1-lr%d%nfl1+1)*(lr%d%nfu2-lr%d%nfl2+1)*(lr%d%nfu3-lr%d%nfl3+1), &
       work%xy_f(1,lr%d%nfl2,lr%d%nfl1,lr%d%nfl3))
  
  
  allocate(work%xz_f4(lr%d%nfl3:lr%d%nfu3,lr%d%nfl1:lr%d%nfu1,lr%d%nfl2:lr%d%nfu2), stat=istat)
  call memocc(istat, work%xz_f4, 'work%xz_f4', subname)
  !!work%xz_f4=0.d0
  call to_zero((lr%d%nfu1-lr%d%nfl1+1)*(lr%d%nfu2-lr%d%nfl2+1)*(lr%d%nfu3-lr%d%nfl3+1), &
       work%xz_f4(lr%d%nfl3,lr%d%nfl1,lr%d%nfl2))
  allocate(work%xz_f(7,lr%d%nfl3:lr%d%nfu3,lr%d%nfl1:lr%d%nfu1,lr%d%nfl2:lr%d%nfu2), stat=istat)
  call memocc(istat, work%xz_f, 'work%xz_f', subname)
  !!work%xz_f=0.d0
  call to_zero(7*(lr%d%nfu1-lr%d%nfl1+1)*(lr%d%nfu2-lr%d%nfl2+1)*(lr%d%nfu3-lr%d%nfl3+1), &
       work%xz_f(1,lr%d%nfl3,lr%d%nfl1,lr%d%nfl2))
  
  
  allocate(work%y_c(0:lr%d%n1,0:lr%d%n2,0:lr%d%n3), stat=istat)
  call memocc(istat, work%y_c, 'work%y_c', subname)
  !!work%y_c=0.d0
  call to_zero((lr%d%n1+1)*(lr%d%n2+1)*(lr%d%n3+1), work%y_c(0,0,0))
  
  allocate(work%y_f(7,lr%d%nfl1:lr%d%nfu1,lr%d%nfl2:lr%d%nfu2,lr%d%nfl3:lr%d%nfu3), stat=istat)
  call memocc(istat, work%y_f, 'work%y_f', subname)
  !!work%y_f=0.d0
  call to_zero(7*(lr%d%nfu1-lr%d%nfl1+1)*(lr%d%nfu2-lr%d%nfl2+1)*(lr%d%nfu3-lr%d%nfl3+1), &
       work%y_f(1,lr%d%nfl1,lr%d%nfl2,lr%d%nfl3))

end subroutine allocate_workarrays_quartic_convolutions



subroutine deallocate_workarrays_quartic_convolutions(lr, subname, work)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  type(locreg_descriptors),intent(in):: lr
  character(len=*),intent(in):: subname
  type(workarrays_quartic_convolutions),intent(out):: work
  
  ! Local variables
  integer:: iall, istat


  iall=-product(shape(work%xx_c))*kind(work%xx_c)
  deallocate(work%xx_c, stat=istat)
  call memocc(istat, iall, 'work%xx_c', subname)

  iall=-product(shape(work%xy_c))*kind(work%xy_c)
  deallocate(work%xy_c, stat=istat)
  call memocc(istat, iall, 'work%xy_c', subname)

  iall=-product(shape(work%xz_c))*kind(work%xz_c)
  deallocate(work%xz_c, stat=istat)
  call memocc(istat, iall, 'work%xz_c', subname)


  iall=-product(shape(work%xx_f1))*kind(work%xx_f1)
  deallocate(work%xx_f1, stat=istat)
  call memocc(istat, iall, 'work%xx_f1', subname)

  iall=-product(shape(work%xx_f))*kind(work%xx_f)
  deallocate(work%xx_f, stat=istat)
  call memocc(istat, iall, 'work%xx_f', subname)

  iall=-product(shape(work%xy_f2))*kind(work%xy_f2)
  deallocate(work%xy_f2, stat=istat)
  call memocc(istat, iall, 'work%xy_f2', subname)


  iall=-product(shape(work%xy_f))*kind(work%xy_f)
  deallocate(work%xy_f, stat=istat)
  call memocc(istat, iall, 'work%xy_f', subname)


  iall=-product(shape(work%xz_f4))*kind(work%xz_f4)
  deallocate(work%xz_f4, stat=istat)
  call memocc(istat, iall, 'work%xz_f4', subname)

  iall=-product(shape(work%xz_f))*kind(work%xz_f)
  deallocate(work%xz_f, stat=istat)
  call memocc(istat, iall, 'work%xz_f', subname)


  iall=-product(shape(work%y_c))*kind(work%y_c)
  deallocate(work%y_c, stat=istat)
  call memocc(istat, iall, 'work%y_c', subname)

  iall=-product(shape(work%y_f))*kind(work%y_f)
  deallocate(work%y_f, stat=istat)
  call memocc(istat, iall, 'work%y_f', subname)

end subroutine deallocate_workarrays_quartic_convolutions



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


