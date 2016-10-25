!> @file
!!   Fake routines for BLACS/ScaLAPACK
!! @author
!!   Copyright (C) 2016 CheSS developers
!!
!!   This file is part of CheSS.
!!   
!!   CheSS is free software: you can redistribute it and/or modify
!!   it under the terms of the GNU Lesser General Public License as published by
!!   the Free Software Foundation, either version 3 of the License, or
!!   (at your option) any later version.
!!   
!!   CheSS is distributed in the hope that it will be useful,
!!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!   GNU Lesser General Public License for more details.
!!   
!!   You should have received a copy of the GNU Lesser General Public License
!!   along with CheSS.  If not, see <http://www.gnu.org/licenses/>.


subroutine blacs_get()
  implicit none
  stop 'FAKE BLACS_GET'
END SUBROUTINE blacs_get

subroutine blacs_gridinit()
  implicit none
  stop 'FAKE BLABS_GRIDINIT'
END SUBROUTINE blacs_gridinit

subroutine blacs_gridinfo()
  implicit none
  stop 'FAKE BLABS_GRIDINFO'
END SUBROUTINE blacs_gridinfo

subroutine descinit()
  implicit none
  stop 'FAKE DESCINIT'
END SUBROUTINE descinit

subroutine pdelset()
  implicit none
  stop 'FAKE PDELSET'
END SUBROUTINE pdelset

subroutine pdsygvx()
  implicit none
  stop 'FAKE PDSYGVX'
END SUBROUTINE pdsygvx

subroutine pdsyevx()
  implicit none
  stop 'FAKE PDSYEGVX'
END SUBROUTINE pdsyevx

subroutine pdelset2()
  implicit none
  stop 'FAKE PDELSET2'
END SUBROUTINE pdelset2

subroutine pdgemm()
  implicit none
  stop 'FAKE PDGEMM'
END SUBROUTINE pdgemm

subroutine pdsymm()
  implicit none
  stop 'FAKE PDSYMM'
END SUBROUTINE pdsymm

integer function numroc()
  implicit none
  numroc=1
  stop 'FAKE NUMROC'
end function numroc

subroutine pdgesv()
  implicit none
  stop 'FAKE PDGESV'
END SUBROUTINE pdgesv

subroutine pdpotrf()
  implicit none
  stop 'FAKE PDPOTRF'
END SUBROUTINE pdpotrf

subroutine pdpotri()
  implicit none
  stop 'FAKE PDPOTRI'
END SUBROUTINE pdpotri

subroutine blacs_gridexit()
  implicit none
  stop 'FAKE BLACS_GRIDEXIT'
END SUBROUTINE blacs_gridexit
