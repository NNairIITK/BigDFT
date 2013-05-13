!> @file
!! Other fortran file for f_malloc routines
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  !guess the rank
  m%rank=0

  if (present(lbounds)) then
     m%rank=size(lbounds)
     m%lbounds(1:m%rank)=lbounds
  end if

  if (present(shape)) then
     if (m%rank == 0) then
        m%rank=size(shape)
     else if (m%rank/=size(shape)) then
        stop 'ERROR, f_malloc: shape not conformal with lbounds'
     end if
     m%shape(1:m%rank)=shape
     do i=1,m%rank
        m%ubounds(i)=m%lbounds(i)+m%shape(i)-1
     end do
     if (present(ubounds)) then
        if (m%rank/=size(ubounds)) stop &
             'ERROR, f_malloc: shape not conformal with ubounds'
        do i=1,m%rank
           if (m%ubounds(i) /=ubounds(i)) stop &
                'ERROR, f_malloc: ubounds not conformal with shape and lbounds'
        end do
     end if
  else
     if (present(ubounds)) then
        if (m%rank == 0) then
           m%rank=size(ubounds)
        else if (m%rank/=size(ubounds)) then
           stop 'ERROR, f_malloc: ubounds not conformal with lbounds'
        end if
        m%ubounds(1:m%rank)=ubounds
        do i=1,m%rank
           m%shape(i)=m%ubounds(i)-m%lbounds(i)+1
        end do
     else
        stop 'ERROR, f_malloc: at least shape or ubounds should be defined'
     end if
  end if

