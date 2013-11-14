!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_build_info_fake
!! NAME
!!  m_build_info_fake
!!
!! FUNCTION
!!  Fake m_build_info module, because abilint is at present (Feb 2009)
!!  unable to treat correctly the m_build_info.F90.in file
!!  inside the m_error.F90 module, to generate the
!!  src/27_toolbox_oop/abinit.dir file.
!!
!! COPYRIGHT
!!  Copyright (C) 2009-2013 ABINIT group (XG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.inc"
#endif

!#include "abi_common.h"


module m_build_info_fake

 use defs_basis
 use m_profiling

 implicit none

contains  !===========================================================
!!***

!!****f* ABINIT/dump_config_fake
!! NAME
!!  dump_config_fake
!!
!! FUNCTION
!!  fake dump_config
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine dump_config_fake()

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dump_config_fake'
!End of the abilint section

 implicit none

!Arguments ------------------------------------

!Local variables-------------------------------
! *********************************************************************
  
 end subroutine dump_config_fake

end module m_build_info_fake
!!***
