!{\src2tex{textfont=tt}}
!!****f* ABINIT/xcast_mpi
!! NAME
!!  xcast_mpi
!!
!! FUNCTION
!!  This module contains functions that calls MPI routine,
!!  if we compile the code using the MPI CPP flags.
!!  xcast_mpi is the generic function.
!!
!! COPYRIGHT
!!  Copyright (C) 2001-2011 ABINIT group (Rshaltaf,AR,XG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~ABINIT/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine xcast_mpi_intv(xval,master,spaceComm,ier)

 use defs_basis
#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 integer,intent(inout) :: xval
 integer,intent(in) :: spaceComm,master
 integer,intent(out) :: ier

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_BCAST(xval,1,MPI_INTEGER,master,spaceComm,ier)
 end if
#endif

end subroutine xcast_mpi_intv
!!***

!!****f* ABINIT/xcast_mpi_int1d
!! NAME
!!  xcast_mpi_int1d
!!
!! FUNCTION
!!  Broadcasts data from master to slaves.
!!  Target: one-dimensional integer arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
!!  master= master MPI node
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xcast_mpi_int1d(xval,master,spaceComm,ier)

 use defs_basis
#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 integer,intent(inout) :: xval(:)
 integer,intent(in) :: spaceComm,master
 integer,intent(out) :: ier
!Local variables-------------------------------
 integer :: n

! *************************************************************************

 ier=0
 n=size(xval)
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_BCAST(xval,n,MPI_INTEGER,master,spaceComm,ier)
 end if
#endif
end subroutine xcast_mpi_int1d
!!***

!!****f* ABINIT/xcast_mpi_int2d
!! NAME
!!  xcast_mpi_int2d
!!
!! FUNCTION
!!  Broadcasts data from master to slaves.
!!  Target: two-dimensional integer arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
!!  master= master MPI node
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xcast_mpi_int2d(xval,master,spaceComm,ier)

 use defs_basis
#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 integer,intent(inout) :: xval(:,:)
 integer,intent(in) :: spaceComm,master
 integer,intent(out) :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n1,n2
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   n1=size(xval,dim=1)
   n2=size(xval,dim=2)
   call MPI_BCAST(xval,n1*n2,MPI_INTEGER,master,spaceComm,ier)
 end if
#endif
end subroutine xcast_mpi_int2d
!!***

!!****f* ABINIT/xcast_mpi_int3d
!! NAME
!!  xcast_mpi_int3d
!!
!! FUNCTION
!!  Broadcasts data from master to slaves.
!!  Target: three-dimensional integer arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
!!  master= master MPI node
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xcast_mpi_int3d(xval,master,spaceComm,ier)

 use defs_basis
#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 integer,intent(inout) :: xval(:,:,:)
 integer,intent(in) :: spaceComm,master
 integer,intent(out) :: ier

!Local variables-------------------------------
#if defined HAVE_MPI
 integer :: n1,n2,n3
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   n1=size(xval,dim=1)
   n2=size(xval,dim=2)
   n3=size(xval,dim=3)
   call MPI_BCAST(xval,n1*n2*n3,MPI_INTEGER,master,spaceComm,ier)
 end if
#endif
end subroutine xcast_mpi_int3d
!!***

!!****f* ABINIT/xcast_mpi_dpv
!! NAME
!!  xcast_mpi_dpv
!!
!! FUNCTION
!!  Broadcasts data from master to slaves.
!!  Target: scalar double precisions.
!!
!! INPUTS
!!  spaceComm= MPI communicator
!!  master= master MPI node
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xcast_mpi_dpv(xval,master,spaceComm,ier)

 use defs_basis
#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 real(dp),intent(inout) :: xval
 integer ,intent(in) :: spaceComm,master
 integer ,intent(out) :: ier
!Local variables-------------------------------

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_BCAST(xval,1,MPI_DOUBLE_PRECISION,master,spaceComm,ier)
 end if
#endif

end subroutine xcast_mpi_dpv
!!***

!!****f* ABINIT/xcast_mpi_dp1d
!! NAME
!!  xcast_mpi_dp1d
!!
!! FUNCTION
!!  Broadcasts data from master to slaves.
!!  Target: double precision one-dimensional arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
!!  master= master MPI node
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xcast_mpi_dp1d(xval,master,spaceComm,ier)

 use defs_basis
#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(inout) :: xval(:)
 integer ,intent(in) :: spaceComm,master
 integer ,intent(out) :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   n=size(xval,dim=1)
   call MPI_BCAST(xval,n,MPI_DOUBLE_PRECISION,master,spaceComm,ier)
 end if
#endif

end subroutine xcast_mpi_dp1d
!!***

!!****f* ABINIT/xcast_mpi_dp2d
!! NAME
!!  xcast_mpi_dp2d
!!
!! FUNCTION
!!  Broadcasts data from master to slaves.
!!  Target: double precision two-dimensional arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
!!  master= master MPI node
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xcast_mpi_dp2d(xval,master,spaceComm,ier)

 use defs_basis
#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(inout) :: xval(:,:)
 integer ,intent(in) :: spaceComm,master
 integer ,intent(out) :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n1,n2
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   n1=size(xval,dim=1)
   n2=size(xval,dim=2)
   call MPI_BCAST(xval,n1*n2,MPI_DOUBLE_PRECISION,master,spaceComm,ier)
 end if
#endif

end subroutine xcast_mpi_dp2d
!!***

!!****f* ABINIT/xcast_mpi_dp3d
!! NAME
!!  xcast_mpi_dp3d
!!
!! FUNCTION
!!  Broadcasts data from master to slaves.
!!  Target: double precision three-dimensional arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
!!  master= master MPI node
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xcast_mpi_dp3d(xval,master,spaceComm,ier)

 use defs_basis
#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(inout) :: xval(:,:,:)
 integer ,intent(in) :: spaceComm,master
 integer ,intent(out) :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n1,n2,n3
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   n1=size(xval,dim=1)
   n2=size(xval,dim=2)
   n3=size(xval,dim=3)
   call MPI_BCAST(xval,n1*n2*n3,MPI_DOUBLE_PRECISION,master,spaceComm,ier)
 end if
#endif

end subroutine xcast_mpi_dp3d
!!***

!!****f* ABINIT/xcast_mpi_dp4d
!! NAME
!!  xcast_mpi_dp4d
!!
!! FUNCTION
!!  Broadcasts data from master to slaves.
!!  Target: double precision four-dimensional arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
!!  master= master MPI node
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xcast_mpi_dp4d(xval,master,spaceComm,ier)

 use defs_basis
#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(inout) :: xval(:,:,:,:)
 integer ,intent(in) :: spaceComm,master
 integer ,intent(out) :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n1,n2,n3,n4
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   n1=size(xval,dim=1)
   n2=size(xval,dim=2)
   n3=size(xval,dim=3)
   n4=size(xval,dim=4)
   call MPI_BCAST(xval,n1*n2*n3*n4,MPI_DOUBLE_PRECISION,master,spaceComm,ier)
 end if
#endif

end subroutine xcast_mpi_dp4d
!!***

!!****f* ABINIT/xcast_mpi_spv
!! NAME
!!  xcast_mpi_spv
!!
!! FUNCTION
!!  Broadcasts data from master to slaves.
!!  Target: scalar single precisions.
!!
!! INPUTS
!!  spaceComm= MPI communicator
!!  master= master MPI node
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xcast_mpi_spv(xval,master,spaceComm,ier)

use defs_basis
#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real,intent(inout) :: xval
 integer,intent(in) :: spaceComm,master
 integer,intent(out) :: ier

!Local variables-------------------

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_BCAST(xval,1,MPI_REAL,master,spaceComm,ier)
 end if
#endif

end subroutine xcast_mpi_spv
!!***

!!****f* ABINIT/xcast_mpi_sp1d
!! NAME
!!  xcast_mpi_sp1d
!!
!! FUNCTION
!!  Broadcasts data from master to slaves.
!!  Target: one-dimensional single precision arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
!!  master= master MPI node
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xcast_mpi_sp1d(xval,master,spaceComm,ier)

 use defs_basis
#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real,intent(inout) :: xval(:)
 integer ,intent(in) :: spaceComm,master
 integer ,intent(out) :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   n=size(xval,dim=1)
   call MPI_BCAST(xval,n,MPI_REAL,master,spaceComm,ier)
 end if
#endif

end subroutine xcast_mpi_sp1d
!!***

!!****f* ABINIT/xcast_mpi_sp2d
!! NAME
!!  xcast_mpi_sp2d
!!
!! FUNCTION
!!  Broadcasts data from master to slaves.
!!  Target: two-dimensional single precision arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
!!  master= master MPI node
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xcast_mpi_sp2d(xval,master,spaceComm,ier)

 use defs_basis
#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real,intent(inout) :: xval(:,:)
 integer ,intent(in) :: spaceComm,master
 integer ,intent(out) :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n1,n2
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   n1=size(xval,dim=1)
   n2=size(xval,dim=2)
   call MPI_BCAST(xval,n1*n2,MPI_REAL,master,spaceComm,ier)
 end if
#endif

end subroutine xcast_mpi_sp2d
!!***

!!****f* ABINIT/xcast_mpi_sp3d
!! NAME
!!  xcast_mpi_sp3d
!!
!! FUNCTION
!!  Broadcasts data from master to slaves.
!!  Target: three-dimensional single precision arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
!!  master= master MPI node
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xcast_mpi_sp3d(xval,master,spaceComm,ier)

 use defs_basis
#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real,intent(inout) :: xval(:,:,:)
 integer ,intent(in) :: spaceComm,master
 integer ,intent(out) :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n1,n2,n3
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   n1=size(xval,dim=1)
   n2=size(xval,dim=2)
   n3=size(xval,dim=3)
   call MPI_BCAST(xval,n1*n2*n3,MPI_REAL,master,spaceComm,ier)
 end if
#endif

end subroutine xcast_mpi_sp3d
!!***

!!****f* ABINIT/xcast_mpi_sp4d
!! NAME
!!  xcast_mpi_sp4d
!!
!! FUNCTION
!!  Broadcasts data from master to slaves.
!!  Target: four-dimensional single precision arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
!!  master= master MPI node
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xcast_mpi_sp4d(xval,master,spaceComm,ier)

 use defs_basis
#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real,intent(inout) :: xval(:,:,:,:)
 integer ,intent(in) :: spaceComm,master
 integer ,intent(out) :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n1,n2,n3,n4
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   n1=size(xval,dim=1)
   n2=size(xval,dim=2)
   n3=size(xval,dim=3)
   n4=size(xval,dim=4)
   call MPI_BCAST(xval,n1*n2*n3*n4,MPI_REAL,master,spaceComm,ier)
 end if
#endif

end subroutine xcast_mpi_sp4d
!!***

!!****f* ABINIT/xcast_mpi_cplxv
!! NAME
!!  xcast_mpi_cplxv
!!
!! FUNCTION
!!  Broadcasts data from master to slaves.
!!  Target: scalar complexs.
!!
!! INPUTS
!!  spaceComm= MPI communicator
!!  master= master MPI node
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xcast_mpi_cplxv(xval,master,spaceComm,ier)

 use defs_basis
#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 complex,intent(inout) :: xval
 integer ,intent(in) :: spaceComm,master
 integer ,intent(out) :: ier

!Local variables-------------------

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_BCAST(xval,1,MPI_COMPLEX,master,spaceComm,ier)
 end if
#endif

end subroutine xcast_mpi_cplxv
!!***

!!****f* ABINIT/xcast_mpi_cplx1d
!! NAME
!!  xcast_mpi_cplx1d
!!
!! FUNCTION
!!  Broadcasts data from master to slaves.
!!  Target: one-dimensional complex arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
!!  master= master MPI node
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xcast_mpi_cplx1d(xval,master,spaceComm,ier)

 use defs_basis
#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 complex,intent(inout) :: xval(:)
 integer ,intent(in) :: spaceComm,master
 integer ,intent(out) :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   n=size(xval(:))
   call MPI_BCAST(xval,n,MPI_COMPLEX,master,spaceComm,ier)
 end if
#endif

end subroutine xcast_mpi_cplx1d
!!***

!!****f* ABINIT/xcast_mpi_cplx2d
!! NAME
!!  xcast_mpi_cplx2d
!!
!! FUNCTION
!!  Broadcasts data from master to slaves.
!!  Target: two-dimensional complex arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
!!  master= master MPI node
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xcast_mpi_cplx2d(xval,master,spaceComm,ier)

 use defs_basis
#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 complex,intent(inout) :: xval(:,:)
 integer ,intent(in) :: spaceComm,master
 integer ,intent(out) :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n1,n2
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   n1=size(xval,dim=1)
   n2=size(xval,dim=2)
   call MPI_BCAST(xval,n1*n2,MPI_COMPLEX,master,spaceComm,ier)
 end if
#endif

end subroutine xcast_mpi_cplx2d
!!***

!!****f* ABINIT/xcast_mpi_cplx3d
!! NAME
!!  xcast_mpi_cplx3d
!!
!! FUNCTION
!!  Broadcasts data from master to slaves.
!!  Target: three-dimensional complex arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
!!  master= master MPI node
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xcast_mpi_cplx3d(xval,master,spaceComm,ier)

 use defs_basis
#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 complex,intent(inout) :: xval(:,:,:)
 integer ,intent(in) :: spaceComm,master
 integer ,intent(out) :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n1,n2,n3
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   n1=size(xval,dim=1)
   n2=size(xval,dim=2)
   n3=size(xval,dim=3)
   call MPI_BCAST(xval,n1*n2*n3,MPI_COMPLEX,master,spaceComm,ier)
 end if
#endif

end subroutine xcast_mpi_cplx3d
!!***

!!****f* ABINIT/xcast_mpi_cplx4d
!! NAME
!!  xcast_mpi_cplx4d
!!
!! FUNCTION
!!  Broadcasts data from master to slaves.
!!  Target: four-dimensional complex arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
!!  master= master MPI node
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xcast_mpi_cplx4d(xval,master,spaceComm,ier)

 use defs_basis
#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 complex,intent(inout) :: xval(:,:,:,:)
 integer ,intent(in) :: spaceComm,master
 integer ,intent(out) :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n1,n2,n3,n4
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   n1=size(xval,dim=1)
   n2=size(xval,dim=2)
   n3=size(xval,dim=3)
   n4=size(xval,dim=4)
   call MPI_BCAST(xval,n1*n2*n3*n4,MPI_COMPLEX,master,spaceComm,ier)
 end if
#endif

end subroutine xcast_mpi_cplx4d
!!***

!!****f* ABINIT/xcast_mpi_dcv
!! NAME
!!  xcast_mpi_dcv
!!
!! FUNCTION
!!  Broadcasts data from master to slaves.
!!  Target: scalar double complexs.
!!
!! INPUTS
!!  spaceComm= MPI communicator
!!  master= master MPI node
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xcast_mpi_dcv(xval,master,spaceComm,ier)

use defs_basis
#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 complex(dpc),intent(inout):: xval
 integer ,intent(in) :: spaceComm,master
 integer ,intent(out) :: ier

!Local variables-------------------

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_BCAST(xval,1,MPI_DOUBLE_COMPLEX,master,spaceComm,ier)
 end if
#endif

end subroutine xcast_mpi_dcv
!!***

!!****f* ABINIT/xcast_mpi_dc1d
!! NAME
!!  xcast_mpi_dc1d
!!
!! FUNCTION
!!  Broadcasts data from master to slaves.
!!  Target: one-dimensional double complex arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
!!  master= master MPI node
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xcast_mpi_dc1d(xval,master,spaceComm,ier)

use defs_basis
#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 complex(dpc),intent(inout):: xval(:)
 integer ,intent(in) :: spaceComm,master
 integer ,intent(out) :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   n=size(xval(:))
   call MPI_BCAST(xval,n,MPI_DOUBLE_COMPLEX,master,spaceComm,ier)
 end if
#endif

end subroutine xcast_mpi_dc1d
!!***

!!****f* ABINIT/xcast_mpi_dc2d
!! NAME
!!  xcast_mpi_dc2d
!!
!! FUNCTION
!!  Broadcasts data from master to slaves.
!!  Target: two-dimensional double complex arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
!!  master= master MPI node
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xcast_mpi_dc2d(xval,master,spaceComm,ier)

use defs_basis
#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 complex(dpc),intent(inout):: xval(:,:)
 integer ,intent(in) :: spaceComm,master
 integer ,intent(out) :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n1,n2
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   n1=size(xval,dim=1)
   n2=size(xval,dim=2)
   call MPI_BCAST(xval,n1*n2,MPI_DOUBLE_COMPLEX,master,spaceComm,ier)
 end if
#endif

end subroutine xcast_mpi_dc2d
!!***

!!****f* ABINIT/xcast_mpi_dc3d
!! NAME
!!  xcast_mpi_dc3d
!!
!! FUNCTION
!!  Broadcasts data from master to slaves.
!!  Target: three-dimensional double complex arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
!!  master= master MPI node
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xcast_mpi_dc3d(xval,master,spaceComm,ier)

use defs_basis
#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 complex(dpc),intent(inout):: xval(:,:,:)
 integer ,intent(in) :: spaceComm,master
 integer ,intent(out) :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n1,n2,n3
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   n1=size(xval,dim=1)
   n2=size(xval,dim=2)
   n3=size(xval,dim=3)
   call MPI_BCAST(xval,n1*n2*n3,MPI_DOUBLE_COMPLEX,master,spaceComm,ier)
 end if
#endif

end subroutine xcast_mpi_dc3d
!!***

!!****f* ABINIT/xcast_mpi_dc4d
!! NAME
!!  xcast_mpi_dc4d
!!
!! FUNCTION
!!  Broadcasts data from master to slaves.
!!  Target: four-dimensional complex arrays in double precision.
!!
!! INPUTS
!!  spaceComm= MPI communicator
!!  master= master MPI node
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xcast_mpi_dc4d(xval,master,spaceComm,ier)

 use defs_basis
#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 complex(dpc),intent(inout) :: xval(:,:,:,:)
 integer,intent(in) :: spaceComm,master
 integer,intent(out) :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n1,n2,n3,n4
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   n1=size(xval,dim=1)
   n2=size(xval,dim=2)
   n3=size(xval,dim=3)
   n4=size(xval,dim=4)
   call MPI_BCAST(xval,n1*n2*n3*n4,MPI_DOUBLE_COMPLEX,master,spaceComm,ier)
 end if
#endif

end subroutine xcast_mpi_dc4d
!!***

!!****f* ABINIT/xcast_mpi_ch0d
!! NAME
!!  xcast_mpi_ch0d
!!
!! FUNCTION
!!  Broadcasts data from master to slaves.
!!  Target: character strings.
!!
!! INPUTS
!!  spaceComm= MPI communicator
!!  master= master MPI node
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xcast_mpi_ch0d(xval,master,spaceComm,ier)

 use defs_basis

#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 character(len=*),intent(inout) :: xval
 integer,intent(in) :: spaceComm,master
 integer,intent(out) :: ier

!Local variables-------------------------------
#if defined HAVE_MPI
 integer :: nch,rank
#endif

!*************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_COMM_RANK(spaceComm,rank,ier)
   if (rank==master) nch=len_trim(xval)
   call MPI_BCAST(nch,1,MPI_INTEGER,master,spaceComm,ier)
   call MPI_BCAST(xval,nch,MPI_CHARACTER,master,spaceComm,ier)
   if (rank/=master) xval(nch+1:)=''
 end if
#endif

end subroutine xcast_mpi_ch0d
!!***

!!****f* ABINIT/xcast_mpi_ch1d
!! NAME
!!  xcast_mpi_ch1d
!!
!! FUNCTION
!!  Broadcasts data from master to slaves.
!!  Target: one-dimensional array of character stringss.
!!
!! INPUTS
!!  spaceComm= MPI communicator
!!  master= master MPI node
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xcast_mpi_ch1d(xval,master,spaceComm,ier)

 use defs_basis

#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 Character(len=*),intent(inout) :: xval(:)
 integer,intent(in) :: spaceComm,master
 integer,intent(out) :: ier

!Local variables-------------------------------
#if defined HAVE_MPI
 integer :: ii,nch
#endif

!*************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   nch=0
   do ii=1,size(xval)
     nch=nch+len(xval(ii))
   end do
   call MPI_BCAST(xval,nch,MPI_CHARACTER,master,spaceComm,ier)
 end if
#endif

end subroutine xcast_mpi_ch1d
!!***
