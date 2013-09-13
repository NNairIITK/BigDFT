!> @file
!! Include fortran file for f_malloc routines
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

  if (present(id)) then
     lgt=min(len(id),namelen)
     m%array_id(1:lgt)=id(1:lgt)
  end if
  if (present(routine_id)) then
     lgt=min(len(routine_id),namelen)
     m%routine_id(1:lgt)=routine_id(1:lgt)
  else
     m%routine_id=mems(ictrl)%present_routine
  end if

  if(present(profile)) m%profile=profile
