!> @file
!! De-Allocation of arrays related to the linear version
!! @author
!!    Copyright (C) 2011-2012 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS



subroutine deallocate_workarrays_quartic_convolutions(work)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  type(workarrays_quartic_convolutions),intent(inout):: work
  
  ! Local variables
  integer:: iall, istat


  call f_free_ptr(work%xx_c)

  call f_free_ptr(work%xy_c)

  call f_free_ptr(work%xz_c)

  call f_free_ptr(work%xx_f1)

  call f_free_ptr(work%xx_f)

  call f_free_ptr(work%xy_f2)

  call f_free_ptr(work%xy_f)

  call f_free_ptr(work%xz_f4)

  call f_free_ptr(work%xz_f)

  call f_free_ptr(work%y_c)

  call f_free_ptr(work%y_f)

  call f_free_ptr(work%aeff0array)

  call f_free_ptr(work%beff0array)

  call f_free_ptr(work%ceff0array)

  call f_free_ptr(work%eeff0array)

  call f_free_ptr(work%aeff0_2array)

  call f_free_ptr(work%beff0_2array)

  call f_free_ptr(work%ceff0_2array)

  call f_free_ptr(work%eeff0_2array)

  call f_free_ptr(work%aeff0_2auxarray)

  call f_free_ptr(work%beff0_2auxarray)

  call f_free_ptr(work%ceff0_2auxarray)

  call f_free_ptr(work%eeff0_2auxarray)

  call f_free_ptr(work%xya_c)

  call f_free_ptr(work%xyc_c)

  call f_free_ptr(work%xza_c)

  call f_free_ptr(work%xzc_c)

  call f_free_ptr(work%yza_c)

  call f_free_ptr(work%yzb_c)

  call f_free_ptr(work%yzc_c)

  call f_free_ptr(work%yze_c)

  call f_free_ptr(work%xya_f)

  call f_free_ptr(work%xyb_f)

  call f_free_ptr(work%xyc_f)

  call f_free_ptr(work%xye_f)

  call f_free_ptr(work%xza_f)

  call f_free_ptr(work%xzb_f)

  call f_free_ptr(work%xzc_f)

  call f_free_ptr(work%xze_f)

  call f_free_ptr(work%yza_f)

  call f_free_ptr(work%yzb_f)

  call f_free_ptr(work%yzc_f)

  call f_free_ptr(work%yze_f)

end subroutine deallocate_workarrays_quartic_convolutions


subroutine init_local_work_arrays(n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3, with_confpot, work)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in)::n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3
  logical,intent(in):: with_confpot
  type(workarrays_quartic_convolutions),intent(inout):: work

  ! Local variables
  integer:: i, istat
  integer,parameter :: lowfil=-14,lupfil=14

  work%xx_c = f_malloc0_ptr((/ 0.to.n1, 0.to.n2, 0.to.n3 /),id='work%xx_c')
  work%xy_c = f_malloc0_ptr((/ 0.to.n2, 0.to.n1, 0.to.n3 /),id='work%xy_c')
  work%xz_c = f_malloc0_ptr((/ 0.to.n3, 0.to.n1, 0.to.n2 /),id='work%xz_c')
  
  work%xx_f1 = f_malloc0_ptr((/ nfl1.to.nfu1, nfl2.to.nfu2, nfl3.to.nfu3 /),id='work%xx_f1')
  work%xx_f = f_malloc0_ptr((/ 1.to.7, nfl1.to.nfu1, nfl2.to.nfu2, nfl3.to.nfu3 /),id='work%xx_f')
  
  
  work%xy_f2 = f_malloc0_ptr((/ nfl2.to.nfu2, nfl1.to.nfu1, nfl3.to.nfu3 /),id='work%xy_f2')
  work%xy_f = f_malloc0_ptr((/ 1.to.7, nfl2.to.nfu2, nfl1.to.nfu1, nfl3.to.nfu3 /),id='work%xy_f')
  
  
  work%xz_f4 = f_malloc0_ptr((/ nfl3.to.nfu3, nfl1.to.nfu1, nfl2.to.nfu2 /),id='work%xz_f4')
  work%xz_f = f_malloc0_ptr((/ 1.to.7, nfl3.to.nfu3, nfl1.to.nfu1, nfl2.to.nfu2 /),id='work%xz_f')
  
  
  work%y_c = f_malloc0_ptr((/ 0.to.n1, 0.to.n2, 0.to.n3 /),id='work%y_c')
  
  work%y_f = f_malloc0_ptr((/ 1.to.7, nfl1.to.nfu1, nfl2.to.nfu2, nfl3.to.nfu3 /),id='work%y_f')

  i=max(n1,n2,n3)
  work%aeff0array = f_malloc0_ptr((/ -3+lowfil.to.lupfil+3, 0.to.i /),id='work%aeff0array')
  work%beff0array = f_malloc0_ptr((/ -3+lowfil.to.lupfil+3, 0.to.i /),id='work%beff0array')
  work%ceff0array = f_malloc0_ptr((/ -3+lowfil.to.lupfil+3, 0.to.i /),id='work%ceff0array')
  work%eeff0array = f_malloc0_ptr((/ lowfil.to.lupfil, 0.to.i /),id='work%eeff0array')
  
  work%aeff0_2array = f_malloc0_ptr((/ -3+lowfil.to.lupfil+3, 0.to.i /),id='work%aeff0_2array')
  work%beff0_2array = f_malloc0_ptr((/ -3+lowfil.to.lupfil+3, 0.to.i /),id='work%beff0_2array')
  work%ceff0_2array = f_malloc0_ptr((/ -3+lowfil.to.lupfil+3, 0.to.i /),id='work%ceff0_2array')
  work%eeff0_2array = f_malloc0_ptr((/ lowfil.to.lupfil, 0.to.i /),id='work%eeff0_2array')
  
  work%aeff0_2auxarray = f_malloc0_ptr((/ -3+lowfil.to.lupfil+3, 0.to.i /),id='work%aeff0_2auxarray')
  work%beff0_2auxarray = f_malloc0_ptr((/ -3+lowfil.to.lupfil+3, 0.to.i /),id='work%beff0_2auxarray')
  work%ceff0_2auxarray = f_malloc0_ptr((/ -3+lowfil.to.lupfil+3, 0.to.i /),id='work%ceff0_2auxarray')
  work%eeff0_2auxarray = f_malloc0_ptr((/ -3+lowfil.to.lupfil+3, 0.to.i /),id='work%eeff0_2auxarray')
  
  work%xya_c = f_malloc_ptr((/ 0.to.n2, 0.to.n1, 0.to.n3 /),id='work%xya_c')
  work%xyc_c = f_malloc_ptr((/ 0.to.n2, 0.to.n1, 0.to.n3 /),id='work%xyc_c')
  if(with_confpot) then
     call to_zero((n1+1)*(n2+1)*(n3+1), work%xya_c(0,0,0))
     call to_zero((n1+1)*(n2+1)*(n3+1), work%xyc_c(0,0,0))
  end if
  
  work%xza_c = f_malloc_ptr((/ 0.to.n3, 0.to.n1, 0.to.n2 /),id='work%xza_c')
  work%xzc_c = f_malloc_ptr((/ 0.to.n3, 0.to.n1, 0.to.n2 /),id='work%xzc_c')
  if(with_confpot) then
     call to_zero((n1+1)*(n2+1)*(n3+1), work%xza_c(0,0,0))
     call to_zero((n1+1)*(n2+1)*(n3+1), work%xzc_c(0,0,0))
  end if
  
  work%yza_c = f_malloc_ptr((/ 0.to.n3, 0.to.n1, 0.to.n2 /),id='work%yza_c')
  work%yzb_c = f_malloc_ptr((/ 0.to.n3, 0.to.n1, 0.to.n2 /),id='work%yzb_c')
  work%yzc_c = f_malloc_ptr((/ 0.to.n3, 0.to.n1, 0.to.n2 /),id='work%yzc_c')
  work%yze_c = f_malloc_ptr((/ 0.to.n3, 0.to.n1, 0.to.n2 /),id='work%yze_c')
  if(with_confpot) then
     call to_zero((n1+1)*(n2+1)*(n3+1), work%yza_c(0,0,0))
     call to_zero((n1+1)*(n2+1)*(n3+1), work%yzb_c(0,0,0))
     call to_zero((n1+1)*(n2+1)*(n3+1), work%yzc_c(0,0,0))
     call to_zero((n1+1)*(n2+1)*(n3+1), work%yze_c(0,0,0))
  end if
  
  work%xya_f = f_malloc_ptr((/ 1.to.3, nfl2.to.nfu2, nfl1.to.nfu1, nfl3.to.nfu3 /),id='work%xya_f')
  work%xyb_f = f_malloc_ptr((/ 1.to.4, nfl2.to.nfu2, nfl1.to.nfu1, nfl3.to.nfu3 /),id='work%xyb_f')
  work%xyc_f = f_malloc_ptr((/ 1.to.3, nfl2.to.nfu2, nfl1.to.nfu1, nfl3.to.nfu3 /),id='work%xyc_f')
  work%xye_f = f_malloc_ptr((/ 1.to.4, nfl2.to.nfu2, nfl1.to.nfu1, nfl3.to.nfu3 /),id='work%xye_f')
  if(with_confpot) then
     call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*3, work%xya_f(1,nfl2,nfl1,nfl3))
     call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*4, work%xyb_f(1,nfl2,nfl1,nfl3))
     call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*3, work%xyc_f(1,nfl2,nfl1,nfl3))
     call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*4, work%xye_f(1,nfl2,nfl1,nfl3))
  end if
  
  work%xza_f = f_malloc_ptr((/ 1.to.3, nfl3.to.nfu3, nfl1.to.nfu1, nfl2.to.nfu2 /),id='work%xza_f')
  work%xzb_f = f_malloc_ptr((/ 1.to.4, nfl3.to.nfu3, nfl1.to.nfu1, nfl2.to.nfu2 /),id='work%xzb_f')
  work%xzc_f = f_malloc_ptr((/ 1.to.3, nfl3.to.nfu3, nfl1.to.nfu1, nfl2.to.nfu2 /),id='work%xzc_f')
  work%xze_f = f_malloc_ptr((/ 1.to.4, nfl3.to.nfu3, nfl1.to.nfu1, nfl2.to.nfu2 /),id='work%xze_f')
  if(with_confpot) then
     call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*3, work%xza_f(1,nfl3,nfl1,nfl2))
     call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*4, work%xzb_f(1,nfl3,nfl1,nfl2))
     call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*3, work%xzc_f(1,nfl3,nfl1,nfl2))
     call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*4, work%xze_f(1,nfl3,nfl1,nfl2))
  end if
  
  work%yza_f = f_malloc_ptr((/ 1.to.3, nfl3.to.nfu3, nfl1.to.nfu1, nfl2.to.nfu2 /),id='work%yza_f')
  work%yzb_f = f_malloc_ptr((/ 1.to.4, nfl3.to.nfu3, nfl1.to.nfu1, nfl2.to.nfu2 /),id='work%yzb_f')
  work%yzc_f = f_malloc_ptr((/ 1.to.3, nfl3.to.nfu3, nfl1.to.nfu1, nfl2.to.nfu2 /),id='work%yzc_f')
  work%yze_f = f_malloc_ptr((/ 1.to.4, nfl3.to.nfu3, nfl1.to.nfu1, nfl2.to.nfu2 /),id='work%yze_f')
  if(with_confpot) then
     call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*3, work%yza_f(1,nfl3,nfl1,nfl2))
     call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*4, work%yzb_f(1,nfl3,nfl1,nfl2))
     call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*3, work%yzc_f(1,nfl3,nfl1,nfl2))
     call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*4, work%yze_f(1,nfl3,nfl1,nfl2))
  end if
  

END SUBROUTINE init_local_work_arrays


subroutine zero_local_work_arrays(n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3, with_confpot, work, subname)
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

  call f_routine(id='zero_local_work_arrays')

  call to_zero((n1+1)*(n2+1)*(n3+1),work%xx_c(0,0,0))!! = f_malloc0_ptr((/ 0.to.n1, 0.to.n2, 0.to.n3 /),id='work%xx_c')
  call to_zero((n2+1)*(n1+1)*(n3+1),work%xy_c(0,0,0))!! = f_malloc0_ptr((/ 0.to.n2, 0.to.n1, 0.to.n3 /),id='work%xy_c')
  call to_zero((n3+1)*(n1+1)*(n2+1),work%xz_c(0,0,0))!! = f_malloc0_ptr((/ 0.to.n3, 0.to.n1, 0.to.n2 /),id='work%xz_c')
  
  call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1),work%xx_f1(nfl1,nfl2,nfl3))!! = f_malloc0_ptr((/ nfl1.to.nfu1, nfl2.to.nfu2, nfl3.to.nfu3 /),id='work%xx_f1')
  call to_zero(7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1),work%xx_f(1,nfl1,nfl2,nfl3))!! = f_malloc0_ptr((/ 1.to.7, nfl1.to.nfu1, nfl2.to.nfu2, nfl3.to.nfu3 /),id='work%xx_f')
  
  
  call to_zero((nfu2-nfl2+1)*(nfu1-nfl1+1)*(nfu3-nfl3+1),work%xy_f2(nfl2,nfl1,nfl3))!! = f_malloc0_ptr((/ nfl2.to.nfu2, nfl1.to.nfu1, nfl3.to.nfu3 /),id='work%xy_f2')
  call to_zero(7*(nfu2-nfl2+1)*(nfu1-nfl1+1)*(nfu3-nfl3+1),work%xy_f(1,nfl2,nfl1,nfl3))!! = f_malloc0_ptr((/ 1.to.7, nfl2.to.nfu2, nfl1.to.nfu1, nfl3.to.nfu3 /),id='work%xy_f')
  
  
  call to_zero((nfu3-nfl3+1)*(nfu1-nfl1+1)*(nfu2-nfl2+1),work%xz_f4(nfl3,nfl1,nfl2))!! = f_malloc0_ptr((/ nfl3.to.nfu3, nfl1.to.nfu1, nfl2.to.nfu2 /),id='work%xz_f4')
  call to_zero(7*(nfu3-nfl3+1)*(nfu1-nfl1+1)*(nfu2-nfl2+1),work%xz_f(1,nfl3,nfl1,nfl2))!! = f_malloc0_ptr((/ 1.to.7, nfl3.to.nfu3, nfl1.to.nfu1, nfl2.to.nfu2 /),id='work%xz_f')
  
  
  call to_zero((n1+1)*(n2+1)*(n3+1),work%y_c(0,0,0))!! = f_malloc0_ptr((/ 0.to.n1, 0.to.n2, 0.to.n3 /),id='work%y_c')
  
  call to_zero(7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1),work%y_f(1,nfl1,nfl2,nfl3))!! = f_malloc0_ptr((/ 1.to.7, nfl1.to.nfu1, nfl2.to.nfu2, nfl3.to.nfu3 /),id='work%y_f')

  i=max(n1,n2,n3)
  call to_zero((lupfil-lowfil+7)*(i+1),work%aeff0array(-3+lowfil,0))!! = f_malloc0_ptr((/ -3+lowfil.to.lupfil+3, 0.to.i /),id='work%aeff0array')
  call to_zero((lupfil-lowfil+7)*(i+1),work%beff0array(-3+lowfil,0))!! = f_malloc0_ptr((/ -3+lowfil.to.lupfil+3, 0.to.i /),id='work%beff0array')
  call to_zero((lupfil-lowfil+7)*(i+1),work%ceff0array(-3+lowfil,0))!! = f_malloc0_ptr((/ -3+lowfil.to.lupfil+3, 0.to.i /),id='work%ceff0array')
  call to_zero((lupfil-lowfil+1)*(i+1),work%eeff0array(lowfil,0))!! = f_malloc0_ptr((/ lowfil.to.lupfil, 0.to.i /),id='work%eeff0array')
  
  call to_zero((lupfil-lowfil+7)*(i+1),work%aeff0_2array(-3+lowfil,0))!! = f_malloc0_ptr((/ -3+lowfil.to.lupfil+3, 0.to.i /),id='work%aeff0_2array')
  call to_zero((lupfil-lowfil+7)*(i+1),work%beff0_2array(-3+lowfil,0))!! = f_malloc0_ptr((/ -3+lowfil.to.lupfil+3, 0.to.i /),id='work%beff0_2array')
  call to_zero((lupfil-lowfil+7)*(i+1),work%ceff0_2array(-3+lowfil,0))!! = f_malloc0_ptr((/ -3+lowfil.to.lupfil+3, 0.to.i /),id='work%ceff0_2array')
  call to_zero((lupfil-lowfil+1)*(i+1),work%eeff0_2array(lowfil,0))!! = f_malloc0_ptr((/ lowfil.to.lupfil, 0.to.i /),id='work%eeff0_2array')
  
  call to_zero((lupfil-lowfil+7)*(i+1),work%aeff0_2auxarray(-3+lowfil,0))!! = f_malloc0_ptr((/ -3+lowfil.to.lupfil+3, 0.to.i /),id='work%aeff0_2auxarray')
  call to_zero((lupfil-lowfil+7)*(i+1),work%beff0_2auxarray(-3+lowfil,0))!! = f_malloc0_ptr((/ -3+lowfil.to.lupfil+3, 0.to.i /),id='work%beff0_2auxarray')
  call to_zero((lupfil-lowfil+7)*(i+1),work%ceff0_2auxarray(-3+lowfil,0))!! = f_malloc0_ptr((/ -3+lowfil.to.lupfil+3, 0.to.i /),id='work%ceff0_2auxarray')
  call to_zero((lupfil-lowfil+7)*(i+1),work%eeff0_2auxarray(-3+lowfil,0))!! = f_malloc0_ptr((/ -3+lowfil.to.lupfil+3, 0.to.i /),id='work%eeff0_2auxarray')
  
  !!work%xya_c = f_malloc_ptr((/ 0.to.n2, 0.to.n1, 0.to.n3 /),id='work%xya_c')
  !!work%xyc_c = f_malloc_ptr((/ 0.to.n2, 0.to.n1, 0.to.n3 /),id='work%xyc_c')
  if(with_confpot) then
     call to_zero((n1+1)*(n2+1)*(n3+1), work%xya_c(0,0,0))
     call to_zero((n1+1)*(n2+1)*(n3+1), work%xyc_c(0,0,0))
  end if
  
  !!work%xza_c = f_malloc_ptr((/ 0.to.n3, 0.to.n1, 0.to.n2 /),id='work%xza_c')
  !!work%xzc_c = f_malloc_ptr((/ 0.to.n3, 0.to.n1, 0.to.n2 /),id='work%xzc_c')
  if(with_confpot) then
     call to_zero((n1+1)*(n2+1)*(n3+1), work%xza_c(0,0,0))
     call to_zero((n1+1)*(n2+1)*(n3+1), work%xzc_c(0,0,0))
  end if
  
  !!work%yza_c = f_malloc_ptr((/ 0.to.n3, 0.to.n1, 0.to.n2 /),id='work%yza_c')
  !!work%yzb_c = f_malloc_ptr((/ 0.to.n3, 0.to.n1, 0.to.n2 /),id='work%yzb_c')
  !!work%yzc_c = f_malloc_ptr((/ 0.to.n3, 0.to.n1, 0.to.n2 /),id='work%yzc_c')
  !!work%yze_c = f_malloc_ptr((/ 0.to.n3, 0.to.n1, 0.to.n2 /),id='work%yze_c')
  if(with_confpot) then
     call to_zero((n1+1)*(n2+1)*(n3+1), work%yza_c(0,0,0))
     call to_zero((n1+1)*(n2+1)*(n3+1), work%yzb_c(0,0,0))
     call to_zero((n1+1)*(n2+1)*(n3+1), work%yzc_c(0,0,0))
     call to_zero((n1+1)*(n2+1)*(n3+1), work%yze_c(0,0,0))
  end if
  
  !!work%xya_f = f_malloc_ptr((/ 1.to.3, nfl2.to.nfu2, nfl1.to.nfu1, nfl3.to.nfu3 /),id='work%xya_f')
  !!work%xyb_f = f_malloc_ptr((/ 1.to.4, nfl2.to.nfu2, nfl1.to.nfu1, nfl3.to.nfu3 /),id='work%xyb_f')
  !!work%xyc_f = f_malloc_ptr((/ 1.to.3, nfl2.to.nfu2, nfl1.to.nfu1, nfl3.to.nfu3 /),id='work%xyc_f')
  !!work%xye_f = f_malloc_ptr((/ 1.to.4, nfl2.to.nfu2, nfl1.to.nfu1, nfl3.to.nfu3 /),id='work%xye_f')
  if(with_confpot) then
     call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*3, work%xya_f(1,nfl2,nfl1,nfl3))
     call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*4, work%xyb_f(1,nfl2,nfl1,nfl3))
     call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*3, work%xyc_f(1,nfl2,nfl1,nfl3))
     call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*4, work%xye_f(1,nfl2,nfl1,nfl3))
  end if
  
  !!work%xza_f = f_malloc_ptr((/ 1.to.3, nfl3.to.nfu3, nfl1.to.nfu1, nfl2.to.nfu2 /),id='work%xza_f')
  !!work%xzb_f = f_malloc_ptr((/ 1.to.4, nfl3.to.nfu3, nfl1.to.nfu1, nfl2.to.nfu2 /),id='work%xzb_f')
  !!work%xzc_f = f_malloc_ptr((/ 1.to.3, nfl3.to.nfu3, nfl1.to.nfu1, nfl2.to.nfu2 /),id='work%xzc_f')
  !!work%xze_f = f_malloc_ptr((/ 1.to.4, nfl3.to.nfu3, nfl1.to.nfu1, nfl2.to.nfu2 /),id='work%xze_f')
  if(with_confpot) then
     call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*3, work%xza_f(1,nfl3,nfl1,nfl2))
     call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*4, work%xzb_f(1,nfl3,nfl1,nfl2))
     call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*3, work%xzc_f(1,nfl3,nfl1,nfl2))
     call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*4, work%xze_f(1,nfl3,nfl1,nfl2))
  end if
  
  !!work%yza_f = f_malloc_ptr((/ 1.to.3, nfl3.to.nfu3, nfl1.to.nfu1, nfl2.to.nfu2 /),id='work%yza_f')
  !!work%yzb_f = f_malloc_ptr((/ 1.to.4, nfl3.to.nfu3, nfl1.to.nfu1, nfl2.to.nfu2 /),id='work%yzb_f')
  !!work%yzc_f = f_malloc_ptr((/ 1.to.3, nfl3.to.nfu3, nfl1.to.nfu1, nfl2.to.nfu2 /),id='work%yzc_f')
  !!work%yze_f = f_malloc_ptr((/ 1.to.4, nfl3.to.nfu3, nfl1.to.nfu1, nfl2.to.nfu2 /),id='work%yze_f')
  if(with_confpot) then
     call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*3, work%yza_f(1,nfl3,nfl1,nfl2))
     call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*4, work%yzb_f(1,nfl3,nfl1,nfl2))
     call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*3, work%yzc_f(1,nfl3,nfl1,nfl2))
     call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)*4, work%yze_f(1,nfl3,nfl1,nfl2))
  end if
  
  call f_release_routine()

END SUBROUTINE zero_local_work_arrays
