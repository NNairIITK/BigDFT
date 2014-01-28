!> @file
!! ODe-Allocation of arrays related to the linear version
!! @author
!!    Copyright (C) 2011-2012 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


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

END SUBROUTINE allocateCommunicationsBuffersPotential



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

END SUBROUTINE deallocateCommunicationsBuffersPotential



!!subroutine allocate_workarrays_quartic_convolutions(lr, subname, work)
!!  use module_base
!!  use module_types
!!  implicit none
!!  
!!  ! Calling arguments
!!  type(locreg_descriptors),intent(in):: lr
!!  character(len=*),intent(in):: subname
!!  type(workarrays_quartic_convolutions),intent(out):: work
!!  
!!  ! Local variables
!!  integer:: istat
!!  
!!  allocate(work%xx_c(0:lr%d%n1,0:lr%d%n2,0:lr%d%n3), stat=istat)
!!  call memocc(istat, work%xx_c, 'work%xx_c', subname)
!!  allocate(work%xy_c(0:lr%d%n2,0:lr%d%n1,0:lr%d%n3), stat=istat)
!!  call memocc(istat, work%xy_c, 'work%xy_c', subname)
!!  allocate(work%xz_c(0:lr%d%n3,0:lr%d%n1,0:lr%d%n2), stat=istat)
!!  call memocc(istat, work%xz_c, 'work%xz_c', subname)
!!
!!  call to_zero((lr%d%n1+1)*(lr%d%n2+1)*(lr%d%n3+1), work%xx_c(0,0,0))
!!  call to_zero((lr%d%n1+1)*(lr%d%n2+1)*(lr%d%n3+1), work%xy_c(0,0,0))
!!  call to_zero((lr%d%n1+1)*(lr%d%n2+1)*(lr%d%n3+1), work%xz_c(0,0,0))
!!  
!!  allocate(work%xx_f1(lr%d%nfl1:lr%d%nfu1,lr%d%nfl2:lr%d%nfu2,lr%d%nfl3:lr%d%nfu3), stat=istat)
!!  call memocc(istat, work%xx_f1, 'work%xx_f1', subname)
!!  call to_zero((lr%d%nfu1-lr%d%nfl1+1)*(lr%d%nfu2-lr%d%nfl2+1)*(lr%d%nfu3-lr%d%nfl3+1), &
!!       work%xx_f1(lr%d%nfl1,lr%d%nfl2,lr%d%nfl3))
!!  allocate(work%xx_f(7,lr%d%nfl1:lr%d%nfu1,lr%d%nfl2:lr%d%nfu2,lr%d%nfl3:lr%d%nfu3), stat=istat)
!!  call memocc(istat, work%xx_f, 'work%xx_f', subname)
!!  call to_zero(7*(lr%d%nfu1-lr%d%nfl1+1)*(lr%d%nfu2-lr%d%nfl2+1)*(lr%d%nfu3-lr%d%nfl3+1), &
!!       work%xx_f(1,lr%d%nfl1,lr%d%nfl2,lr%d%nfl3))
!!  
!!  
!!  allocate(work%xy_f2(lr%d%nfl2:lr%d%nfu2,lr%d%nfl1:lr%d%nfu1,lr%d%nfl3:lr%d%nfu3), stat=istat)
!!  call memocc(istat, work%xy_f2, 'work%xy_f2', subname)
!!  call to_zero((lr%d%nfu1-lr%d%nfl1+1)*(lr%d%nfu2-lr%d%nfl2+1)*(lr%d%nfu3-lr%d%nfl3+1), &
!!       work%xy_f2(lr%d%nfl2,lr%d%nfl1,lr%d%nfl3))
!!  allocate(work%xy_f(7,lr%d%nfl2:lr%d%nfu2,lr%d%nfl1:lr%d%nfu1,lr%d%nfl3:lr%d%nfu3), stat=istat)
!!  call memocc(istat, work%xy_f, 'work%xy_f', subname)
!!  call to_zero(7*(lr%d%nfu1-lr%d%nfl1+1)*(lr%d%nfu2-lr%d%nfl2+1)*(lr%d%nfu3-lr%d%nfl3+1), &
!!       work%xy_f(1,lr%d%nfl2,lr%d%nfl1,lr%d%nfl3))
!!  
!!  
!!  allocate(work%xz_f4(lr%d%nfl3:lr%d%nfu3,lr%d%nfl1:lr%d%nfu1,lr%d%nfl2:lr%d%nfu2), stat=istat)
!!  call memocc(istat, work%xz_f4, 'work%xz_f4', subname)
!!  call to_zero((lr%d%nfu1-lr%d%nfl1+1)*(lr%d%nfu2-lr%d%nfl2+1)*(lr%d%nfu3-lr%d%nfl3+1), &
!!       work%xz_f4(lr%d%nfl3,lr%d%nfl1,lr%d%nfl2))
!!  allocate(work%xz_f(7,lr%d%nfl3:lr%d%nfu3,lr%d%nfl1:lr%d%nfu1,lr%d%nfl2:lr%d%nfu2), stat=istat)
!!  call memocc(istat, work%xz_f, 'work%xz_f', subname)
!!  call to_zero(7*(lr%d%nfu1-lr%d%nfl1+1)*(lr%d%nfu2-lr%d%nfl2+1)*(lr%d%nfu3-lr%d%nfl3+1), &
!!       work%xz_f(1,lr%d%nfl3,lr%d%nfl1,lr%d%nfl2))
!!  
!!  
!!  allocate(work%y_c(0:lr%d%n1,0:lr%d%n2,0:lr%d%n3), stat=istat)
!!  call memocc(istat, work%y_c, 'work%y_c', subname)
!!  call to_zero((lr%d%n1+1)*(lr%d%n2+1)*(lr%d%n3+1), work%y_c(0,0,0))
!!  
!!  allocate(work%y_f(7,lr%d%nfl1:lr%d%nfu1,lr%d%nfl2:lr%d%nfu2,lr%d%nfl3:lr%d%nfu3), stat=istat)
!!  call memocc(istat, work%y_f, 'work%y_f', subname)
!!  call to_zero(7*(lr%d%nfu1-lr%d%nfl1+1)*(lr%d%nfu2-lr%d%nfl2+1)*(lr%d%nfu3-lr%d%nfl3+1), &
!!       work%y_f(1,lr%d%nfl1,lr%d%nfl2,lr%d%nfl3))
!!
!!end subroutine allocate_workarrays_quartic_convolutions



subroutine deallocate_workarrays_quartic_convolutions(lr, subname, work)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  type(locreg_descriptors),intent(in):: lr
  character(len=*),intent(in):: subname
  type(workarrays_quartic_convolutions),intent(inout):: work
  
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


  iall=-product(shape(work%aeff0array))*kind(work%aeff0array)
  deallocate(work%aeff0array, stat=istat)
  call memocc(istat, iall, 'work%aeff0array', subname)

  iall=-product(shape(work%beff0array))*kind(work%beff0array)
  deallocate(work%beff0array, stat=istat)
  call memocc(istat, iall, 'work%beff0array', subname)

  iall=-product(shape(work%ceff0array))*kind(work%ceff0array)
  deallocate(work%ceff0array, stat=istat)
  call memocc(istat, iall, 'work%ceff0array', subname)

  iall=-product(shape(work%eeff0array))*kind(work%eeff0array)
  deallocate(work%eeff0array, stat=istat)
  call memocc(istat, iall, 'work%eeff0array', subname)


  iall=-product(shape(work%aeff0_2array))*kind(work%aeff0_2array)
  deallocate(work%aeff0_2array, stat=istat)
  call memocc(istat, iall, 'work%aeff0_2array', subname)

  iall=-product(shape(work%beff0_2array))*kind(work%beff0_2array)
  deallocate(work%beff0_2array, stat=istat)
  call memocc(istat, iall, 'work%beff0_2array', subname)

  iall=-product(shape(work%ceff0_2array))*kind(work%ceff0_2array)
  deallocate(work%ceff0_2array, stat=istat)
  call memocc(istat, iall, 'work%ceff0_2array', subname)

  iall=-product(shape(work%eeff0_2array))*kind(work%eeff0_2array)
  deallocate(work%eeff0_2array, stat=istat)
  call memocc(istat, iall, 'work%eeff0_2array', subname)


  iall=-product(shape(work%aeff0_2auxarray))*kind(work%aeff0_2auxarray)
  deallocate(work%aeff0_2auxarray, stat=istat)
  call memocc(istat, iall, 'work%aeff0_2auxarray', subname)

  iall=-product(shape(work%beff0_2auxarray))*kind(work%beff0_2auxarray)
  deallocate(work%beff0_2auxarray, stat=istat)
  call memocc(istat, iall, 'work%beff0_2auxarray', subname)

  iall=-product(shape(work%ceff0_2auxarray))*kind(work%ceff0_2auxarray)
  deallocate(work%ceff0_2auxarray, stat=istat)
  call memocc(istat, iall, 'work%ceff0_2auxarray', subname)

  iall=-product(shape(work%eeff0_2auxarray))*kind(work%eeff0_2auxarray)
  deallocate(work%eeff0_2auxarray, stat=istat)
  call memocc(istat, iall, 'work%eeff0_2auxarray', subname)


  iall=-product(shape(work%xya_c))*kind(work%xya_c)
  deallocate(work%xya_c, stat=istat)
  call memocc(istat, iall, 'work%xya_c', subname)

  iall=-product(shape(work%xyc_c))*kind(work%xyc_c)
  deallocate(work%xyc_c, stat=istat)
  call memocc(istat, iall, 'work%xyc_c', subname)


  iall=-product(shape(work%xza_c))*kind(work%xza_c)
  deallocate(work%xza_c, stat=istat)
  call memocc(istat, iall, 'work%xza_c', subname)

  iall=-product(shape(work%xzc_c))*kind(work%xzc_c)
  deallocate(work%xzc_c, stat=istat)
  call memocc(istat, iall, 'work%xzc_c', subname)


  iall=-product(shape(work%yza_c))*kind(work%yza_c)
  deallocate(work%yza_c, stat=istat)
  call memocc(istat, iall, 'work%yza_c', subname)

  iall=-product(shape(work%yzb_c))*kind(work%yzb_c)
  deallocate(work%yzb_c, stat=istat)
  call memocc(istat, iall, 'work%yzb_c', subname)

  iall=-product(shape(work%yzc_c))*kind(work%yzc_c)
  deallocate(work%yzc_c, stat=istat)
  call memocc(istat, iall, 'work%yzc_c', subname)

  iall=-product(shape(work%yze_c))*kind(work%yze_c)
  deallocate(work%yze_c, stat=istat)
  call memocc(istat, iall, 'work%yze_c', subname)


  iall=-product(shape(work%xya_f))*kind(work%xya_f)
  deallocate(work%xya_f, stat=istat)
  call memocc(istat, iall, 'work%xya_f', subname)
  iall=-product(shape(work%xyb_f))*kind(work%xyb_f)
  deallocate(work%xyb_f, stat=istat)
  call memocc(istat, iall, 'work%xyb_f', subname)
  iall=-product(shape(work%xyc_f))*kind(work%xyc_f)
  deallocate(work%xyc_f, stat=istat)
  call memocc(istat, iall, 'work%xyc_f', subname)
  iall=-product(shape(work%xye_f))*kind(work%xye_f)
  deallocate(work%xye_f, stat=istat)
  call memocc(istat, iall, 'work%yze_f7', subname)

  iall=-product(shape(work%xza_f))*kind(work%xza_f)
  deallocate(work%xza_f, stat=istat)
  call memocc(istat, iall, 'work%xza_f', subname)
  iall=-product(shape(work%xzb_f))*kind(work%xzb_f)
  deallocate(work%xzb_f, stat=istat)
  call memocc(istat, iall, 'work%xzb_f', subname)
  iall=-product(shape(work%xzc_f))*kind(work%xzc_f)
  deallocate(work%xzc_f, stat=istat)
  call memocc(istat, iall, 'work%xzc_f', subname)
  iall=-product(shape(work%xze_f))*kind(work%xze_f)
  deallocate(work%xze_f, stat=istat)
  call memocc(istat, iall, 'zze_f7', subname)

  iall=-product(shape(work%yza_f))*kind(work%yza_f)
  deallocate(work%yza_f, stat=istat)
  call memocc(istat, iall, 'work%yza_f', subname)
  iall=-product(shape(work%yzb_f))*kind(work%yzb_f)
  deallocate(work%yzb_f, stat=istat)
  call memocc(istat, iall, 'work%yzb_f', subname)
  iall=-product(shape(work%yzc_f))*kind(work%yzc_f)
  deallocate(work%yzc_f, stat=istat)
  call memocc(istat, iall, 'work%yzc_f', subname)
  iall=-product(shape(work%yze_f))*kind(work%yze_f)
  deallocate(work%yze_f, stat=istat)
  call memocc(istat, iall, 'work%yze_f', subname)

end subroutine deallocate_workarrays_quartic_convolutions


subroutine init_local_work_arrays(n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3, with_confpot, work, subname)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in)::n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3
  logical,intent(in):: with_confpot
  type(workarrays_quartic_convolutions),intent(inout):: work
  character(len=*),intent(in):: subname

  ! Local variables
  integer:: i, istat
  integer,parameter :: lowfil=-14,lupfil=14

  allocate(work%xx_c(0:n1,0:n2,0:n3), stat=istat)
  call memocc(istat, work%xx_c, 'work%xx_c', subname)
  allocate(work%xy_c(0:n2,0:n1,0:n3), stat=istat)
  call memocc(istat, work%xy_c, 'work%xy_c', subname)
  allocate(work%xz_c(0:n3,0:n1,0:n2), stat=istat)
  call memocc(istat, work%xz_c, 'work%xz_c', subname)

  call to_zero((n1+1)*(n2+1)*(n3+1), work%xx_c(0,0,0))
  call to_zero((n1+1)*(n2+1)*(n3+1), work%xy_c(0,0,0))
  call to_zero((n1+1)*(n2+1)*(n3+1), work%xz_c(0,0,0))
  
  allocate(work%xx_f1(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), stat=istat)
  call memocc(istat, work%xx_f1, 'work%xx_f1', subname)
  call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1), &
       work%xx_f1(nfl1,nfl2,nfl3))
  allocate(work%xx_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), stat=istat)
  call memocc(istat, work%xx_f, 'work%xx_f', subname)
  call to_zero(7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1), &
       work%xx_f(1,nfl1,nfl2,nfl3))
  
  
  allocate(work%xy_f2(nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), stat=istat)
  call memocc(istat, work%xy_f2, 'work%xy_f2', subname)
  call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1), &
       work%xy_f2(nfl2,nfl1,nfl3))
  allocate(work%xy_f(7,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), stat=istat)
  call memocc(istat, work%xy_f, 'work%xy_f', subname)
  call to_zero(7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1), &
       work%xy_f(1,nfl2,nfl1,nfl3))
  
  
  allocate(work%xz_f4(nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
  call memocc(istat, work%xz_f4, 'work%xz_f4', subname)
  call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1), &
       work%xz_f4(nfl3,nfl1,nfl2))
  allocate(work%xz_f(7,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
  call memocc(istat, work%xz_f, 'work%xz_f', subname)
  call to_zero(7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1), &
       work%xz_f(1,nfl3,nfl1,nfl2))
  
  
  allocate(work%y_c(0:n1,0:n2,0:n3), stat=istat)
  call memocc(istat, work%y_c, 'work%y_c', subname)
  call to_zero((n1+1)*(n2+1)*(n3+1), work%y_c(0,0,0))
  
  allocate(work%y_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), stat=istat)
  call memocc(istat, work%y_f, 'work%y_f', subname)
  call to_zero(7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1), &
       work%y_f(1,nfl1,nfl2,nfl3))

  i=max(n1,n2,n3)
  allocate(work%aeff0array(-3+lowfil:lupfil+3,0:i), stat=istat)
  call memocc(istat, work%aeff0array, 'work%aeff0array', subname)
  allocate(work%beff0array(-3+lowfil:lupfil+3,0:i), stat=istat)
  call memocc(istat, work%beff0array, 'work%beff0array', subname)
  allocate(work%ceff0array(-3+lowfil:lupfil+3,0:i), stat=istat)
  call memocc(istat, work%ceff0array, 'work%ceff0array', subname)
  allocate(work%eeff0array(lowfil:lupfil,0:i), stat=istat)
  call memocc(istat, work%eeff0array, 'work%eeff0array', subname)
  call to_zero((i+1)*(lupfil-lowfil+7), work%aeff0array(-3+lowfil,0))
  call to_zero((i+1)*(lupfil-lowfil+7), work%beff0array(-3+lowfil,0))
  call to_zero((i+1)*(lupfil-lowfil+7), work%ceff0array(-3+lowfil,0))
  call to_zero((i+1)*(lupfil-lowfil+1), work%eeff0array(lowfil,0))
  
  allocate(work%aeff0_2array(-3+lowfil:lupfil+3,0:i), stat=istat)
  call memocc(istat, work%aeff0_2array, 'work%aeff0_2array', subname)
  allocate(work%beff0_2array(-3+lowfil:lupfil+3,0:i), stat=istat)
  call memocc(istat, work%beff0_2array, 'work%beff0_2array', subname)
  allocate(work%ceff0_2array(-3+lowfil:lupfil+3,0:i), stat=istat)
  call memocc(istat, work%ceff0_2array, 'work%ceff0_2array', subname)
  allocate(work%eeff0_2array(lowfil:lupfil,0:i), stat=istat)
  call memocc(istat, work%eeff0_2array, 'work%eeff0_2array', subname)
  call to_zero((i+1)*(lupfil-lowfil+7), work%aeff0_2array(-3+lowfil,0))
  call to_zero((i+1)*(lupfil-lowfil+7), work%beff0_2array(-3+lowfil,0))
  call to_zero((i+1)*(lupfil-lowfil+7), work%ceff0_2array(-3+lowfil,0))
  call to_zero((i+1)*(lupfil-lowfil+1), work%eeff0_2array(lowfil,0))
  
  allocate(work%aeff0_2auxarray(-3+lowfil:lupfil+3,0:i), stat=istat)
  call memocc(istat, work%aeff0_2auxarray, 'work%aeff0_2auxarray', subname)
  allocate(work%beff0_2auxarray(-3+lowfil:lupfil+3,0:i), stat=istat)
  call memocc(istat, work%beff0_2auxarray, 'work%beff0_2auxarray', subname)
  allocate(work%ceff0_2auxarray(-3+lowfil:lupfil+3,0:i), stat=istat)
  call memocc(istat, work%ceff0_2auxarray, 'work%ceff0_2auxarray', subname)
  allocate(work%eeff0_2auxarray(-3+lowfil:lupfil+3,0:i), stat=istat)
  call memocc(istat, work%eeff0_2auxarray, 'work%eeff0_2auxarray', subname)
  call to_zero((i+1)*(lupfil-lowfil+7), work%aeff0_2auxarray(-3+lowfil,0))
  call to_zero((i+1)*(lupfil-lowfil+7), work%beff0_2auxarray(-3+lowfil,0))
  call to_zero((i+1)*(lupfil-lowfil+7), work%ceff0_2auxarray(-3+lowfil,0))
  call to_zero((i+1)*(lupfil-lowfil+7), work%eeff0_2auxarray(-3+lowfil,0))
  
  allocate(work%xya_c(0:n2,0:n1,0:n3), stat=istat)
  call memocc(istat, work%xya_c, 'work%xya_c', subname)
  allocate(work%xyc_c(0:n2,0:n1,0:n3), stat=istat)
  call memocc(istat, work%xyc_c, 'work%xyc_c', subname)
  if(with_confpot) then
     call to_zero((n1+1)*(n2+1)*(n3+1), work%xya_c(0,0,0))
     call to_zero((n1+1)*(n2+1)*(n3+1), work%xyc_c(0,0,0))
  end if
  
  allocate(work%xza_c(0:n3,0:n1,0:n2), stat=istat)
  call memocc(istat, work%xza_c, 'work%xza_c', subname)
  allocate(work%xzc_c(0:n3,0:n1,0:n2), stat=istat)
  call memocc(istat, work%xzc_c, 'work%xzc_c', subname)
  if(with_confpot) then
     call to_zero((n1+1)*(n2+1)*(n3+1), work%xza_c(0,0,0))
     call to_zero((n1+1)*(n2+1)*(n3+1), work%xzc_c(0,0,0))
  end if
  
  allocate(work%yza_c(0:n3,0:n1,0:n2), stat=istat)
  call memocc(istat, work%yza_c, 'work%yza_c', subname)
  allocate(work%yzb_c(0:n3,0:n1,0:n2), stat=istat)
  call memocc(istat, work%yzb_c, 'work%yzb_c', subname)
  allocate(work%yzc_c(0:n3,0:n1,0:n2), stat=istat)
  call memocc(istat, work%yzc_c, 'work%yzc_c', subname)
  allocate(work%yze_c(0:n3,0:n1,0:n2), stat=istat)
  call memocc(istat, work%yze_c, 'work%yze_c', subname)
  if(with_confpot) then
     call to_zero((n1+1)*(n2+1)*(n3+1), work%yza_c(0,0,0))
     call to_zero((n1+1)*(n2+1)*(n3+1), work%yzb_c(0,0,0))
     call to_zero((n1+1)*(n2+1)*(n3+1), work%yzc_c(0,0,0))
     call to_zero((n1+1)*(n2+1)*(n3+1), work%yze_c(0,0,0))
  end if
  
  allocate(work%xya_f(3,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), stat=istat)
  call memocc(istat, work%xya_f, 'work%xya_f', subname)
  allocate(work%xyb_f(4,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), stat=istat)
  call memocc(istat, work%xyb_f, 'work%xyb_f', subname)
  allocate(work%xyc_f(3,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), stat=istat)
  call memocc(istat, work%xyc_f, 'work%xyc_f', subname)
  allocate(work%xye_f(4,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), stat=istat)
  call memocc(istat, work%xye_f, 'work%xye_f', subname)
  if(with_confpot) then
     call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*3, work%xya_f(1,nfl2,nfl1,nfl3))
     call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*4, work%xyb_f(1,nfl2,nfl1,nfl3))
     call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*3, work%xyc_f(1,nfl2,nfl1,nfl3))
     call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*4, work%xye_f(1,nfl2,nfl1,nfl3))
  end if
  
  allocate(work%xza_f(3,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
  call memocc(istat, work%xza_f, 'work%xza_f', subname)
  allocate(work%xzb_f(4,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
  call memocc(istat, work%xzb_f, 'work%xzb_f', subname)
  allocate(work%xzc_f(3,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
  call memocc(istat, work%xzc_f, 'work%xzc_f', subname)
  allocate(work%xze_f(4,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
  call memocc(istat, work%xze_f, 'work%xze_f', subname)
  if(with_confpot) then
     call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*3, work%xza_f(1,nfl3,nfl1,nfl2))
     call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*4, work%xzb_f(1,nfl3,nfl1,nfl2))
     call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*3, work%xzc_f(1,nfl3,nfl1,nfl2))
     call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*4, work%xze_f(1,nfl3,nfl1,nfl2))
  end if
  
  allocate(work%yza_f(3,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
  call memocc(istat, work%yza_f, 'work%yza_f', subname)
  allocate(work%yzb_f(4,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
  call memocc(istat, work%yzb_f, 'work%yzb_f', subname)
  allocate(work%yzc_f(3,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
  call memocc(istat, work%yzc_f, 'work%yzc_f', subname)
  allocate(work%yze_f(4,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
  call memocc(istat, work%yze_f, 'work%yze_f', subname)
  if(with_confpot) then
     call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*3, work%yza_f(1,nfl3,nfl1,nfl2))
     call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*4, work%yzb_f(1,nfl3,nfl1,nfl2))
     call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*3, work%yzc_f(1,nfl3,nfl1,nfl2))
     call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*4, work%yze_f(1,nfl3,nfl1,nfl2))
  end if
  
  
!!$  call to_zero(lupfil-lowfil+7, work%aeff0(-3+lowfil))
!!$  call to_zero(lupfil-lowfil+7, work%aeff1(-3+lowfil))
!!$  call to_zero(lupfil-lowfil+7, work%aeff2(-3+lowfil))
!!$  call to_zero(lupfil-lowfil+7, work%aeff3(-3+lowfil))
!!$  
!!$  call to_zero(lupfil-lowfil+7, work%beff0(-3+lowfil))
!!$  call to_zero(lupfil-lowfil+7, work%beff1(-3+lowfil))
!!$  call to_zero(lupfil-lowfil+7, work%beff2(-3+lowfil))
!!$  call to_zero(lupfil-lowfil+7, work%beff3(-3+lowfil))
!!$  
!!$  call to_zero(lupfil-lowfil+7, work%ceff0(-3+lowfil))
!!$  call to_zero(lupfil-lowfil+7, work%ceff1(-3+lowfil))
!!$  call to_zero(lupfil-lowfil+7, work%ceff2(-3+lowfil))
!!$  call to_zero(lupfil-lowfil+7, work%ceff3(-3+lowfil))
!!$  
!!$  call to_zero(lupfil-lowfil+1, work%eeff0(lowfil))
!!$  call to_zero(lupfil-lowfil+1, work%eeff1(lowfil))
!!$  call to_zero(lupfil-lowfil+1, work%eeff2(lowfil))
!!$  call to_zero(lupfil-lowfil+1, work%eeff3(lowfil))
!!$  
!!$  
!!$  call to_zero(lupfil-lowfil+7, work%aeff0_2(-3+lowfil))
!!$  call to_zero(lupfil-lowfil+7, work%aeff1_2(-3+lowfil))
!!$  call to_zero(lupfil-lowfil+7, work%aeff2_2(-3+lowfil))
!!$  call to_zero(lupfil-lowfil+7, work%aeff3_2(-3+lowfil))
!!$  
!!$  call to_zero(lupfil-lowfil+7, work%beff0_2(-3+lowfil))
!!$  call to_zero(lupfil-lowfil+7, work%beff1_2(-3+lowfil))
!!$  call to_zero(lupfil-lowfil+7, work%beff2_2(-3+lowfil))
!!$  call to_zero(lupfil-lowfil+7, work%beff3_2(-3+lowfil))
!!$  
!!$  call to_zero(lupfil-lowfil+7, work%ceff0_2(-3+lowfil))
!!$  call to_zero(lupfil-lowfil+7, work%ceff1_2(-3+lowfil))
!!$  call to_zero(lupfil-lowfil+7, work%ceff2_2(-3+lowfil))
!!$  call to_zero(lupfil-lowfil+7, work%ceff3_2(-3+lowfil))
!!$  
!!$  call to_zero(lupfil-lowfil+1, work%eeff0_2(lowfil))
!!$  call to_zero(lupfil-lowfil+1, work%eeff1_2(lowfil))
!!$  call to_zero(lupfil-lowfil+1, work%eeff2_2(lowfil))
!!$  call to_zero(lupfil-lowfil+1, work%eeff3_2(lowfil))

END SUBROUTINE init_local_work_arrays
