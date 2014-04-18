!> @file
!! Include fortran file for memcpy interfaces
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  if (n1 /= n2) then
     call f_err_throw('Error in f_memcpy; the size of the source ('//trim(yaml_toa(n2))//&
          ') and of the destination buffer ('//trim(yaml_toa(n1))//&
          ') does not coincide',err_id=ERR_INVALID_COPY)
     return
  end if
  if (n1 <=0) return
  call c_memcopy(dest,f_loc(src),n1*kind(dest))
