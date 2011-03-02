!{\src2tex{textfont=tt}}
!!****f* ABINIT/xsum_mpi
!! NAME
!!  xsum_mpi
!!
!! FUNCTION
!!  This module contains functions that calls MPI routine,
!!  if we compile the code using the MPI CPP flags.
!!  xsum_mpi is the generic function.
!!
!! COPYRIGHT
!!  Copyright (C) 2001-2011 ABINIT group (AR,XG,MB)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~ABINIT/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! TODO
!!
!! NOTES
!!  MG: The additional array xsum is not needed if MPI2 is available,
!!  MPI2 defines an option MPI_IN_PLACE to do the SUM in-place in case
!!  of intra-communicators.
!!  It would be useful to have a wrapper for the in-place option
!!  since it allows for memory saving if the array to be summed up is large.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.inc"
#endif

subroutine xsum_mpi_int(xval,spaceComm,ier)

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
 integer,intent(in) :: spaceComm
 integer,intent(out)   :: ier

!Local variables-------------------------------
#if defined HAVE_MPI
 integer :: n1, istat
 integer , allocatable :: xsum(:)
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
!  Accumulate xval on all proc. in spaceComm
   n1 = size(xval)
   allocate(xsum(n1), STAT=istat)
   if (istat /= 0) stop 'error allocating xsum in xsum_mpi_int'
   call MPI_ALLREDUCE(xval,xsum,n1,MPI_INTEGER,&
&   MPI_SUM,spaceComm,ier)
   xval (:) = xsum(:)
   deallocate(xsum)
 end if
#endif
end subroutine xsum_mpi_int
!!***

!!****f* ABINIT/xsum_mpi_intv
!! NAME
!!  xsum_mpi_intv
!!
!! FUNCTION
!!  Combines values from all processes and distribute
!!  the result back to all processes.
!!  Target: scalar integers.
!!
!! INPUTS
!!  spaceComm= MPI communicator
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
subroutine xsum_mpi_intv(xval,spaceComm,ier)

 use defs_basis

#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments----------------------
 integer,intent(inout) :: xval
 integer,intent(in) :: spaceComm
 integer,intent(out)   :: ier

!Local variables----------------
#if defined HAVE_MPI
 integer  :: xsum
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
!  Accumulate xval on all proc. in spaceComm
   call MPI_ALLREDUCE(xval,xsum,1,MPI_INTEGER,MPI_SUM,spaceComm,ier)
   xval = xsum
 end if
#endif
end subroutine xsum_mpi_intv
!!***

!!****f* ABINIT/xsum_mpi_intv2
!! NAME
!!  xsum_mpi_intv2
!!
!! FUNCTION
!!  Combines values from all processes and distribute
!!  the result back to all processes.
!!  Target: scalar integer without transfers.
!!
!! INPUTS
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!  xsum= receive buffer
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xsum_mpi_intv2(xval,xsum,spaceComm,ier)

 use defs_basis

#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments---------------------
 integer,intent(inout) :: xval,xsum
 integer,intent(in) :: spaceComm
 integer,intent(out)   :: ier

!Local variables---------------

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
!  Accumulate xval on all proc. in spaceComm
   call MPI_ALLREDUCE(xval,xsum,1,MPI_INTEGER,&
&   MPI_SUM,spaceComm,ier)
 else
#endif
   xsum=xval
#if defined HAVE_MPI
 end if
#endif
end subroutine xsum_mpi_intv2
!!***

!!****f* ABINIT/xsum_mpi_intn
!! NAME
!!  xsum_mpi_intn
!!
!! FUNCTION
!!  Combines values from all processes and distribute
!!  the result back to all processes.
!!  Target: one-dimensional integer arrays.
!!
!! INPUTS
!!  n1= first dimension of the array
!!  spaceComm= MPI communicator
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
subroutine xsum_mpi_intn(xval,n1,spaceComm,ier)

 use defs_basis

#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 integer,intent(inout) :: xval(:)
 integer,intent(in)    :: n1
 integer,intent(in) :: spaceComm
 integer,intent(out)   :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer , allocatable :: xsum(:)
 integer :: nproc_space_comm, istat
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
!Accumulate xval on all proc. in spaceComm
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
   if (nproc_space_comm /= 1) then
     allocate(xsum(n1), STAT=istat)
     if (istat /= 0) stop 'error allocating xsum in xsum_mpi_intn'
     call MPI_ALLREDUCE(xval,xsum,n1,MPI_INTEGER,&
&     MPI_SUM,spaceComm,ier)
     xval (:) = xsum(:)
     deallocate(xsum)
   end if
 end if
#endif
end subroutine xsum_mpi_intn
!!***

!!****f* ABINIT/xsum_mpi_int2t
!! NAME
!!  xsum_mpi_int2t
!!
!! FUNCTION
!!  Combines values from all processes and distribute
!!  the result back to all processes.
!!  Target: one-dimensional integer array without transfers.
!!
!! INPUTS
!!  n1= first dimension of the array
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!  xsum= receive buffer
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xsum_mpi_int2t(xval,xsum,n1,spaceComm,ier)

 use defs_basis

#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none
#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 integer ,intent(inout) :: xval(:),xsum(:)
 integer ,intent(in)    :: n1
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
!  Accumulate xval on all proc. in spaceComm
   call MPI_ALLREDUCE(xval,xsum,n1,MPI_INTEGER,&
&   MPI_SUM,spaceComm,ier)
 else
#endif
   xsum=xval
#if defined HAVE_MPI
 end if
#endif
end subroutine xsum_mpi_int2t
!!***

!!****f* ABINIT/xsum_mpi_int2d
!! NAME
!!  xsum_mpi_int2d
!!
!! FUNCTION
!!  Combines values from all processes and distribute
!!  the result back to all processes.
!!  Target: two-dimensional integer arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
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
subroutine xsum_mpi_int2d(xval,spaceComm,ier)

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
 integer,intent(in) :: spaceComm
 integer,intent(out)   :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer ::   n1,n2, istat
 integer , allocatable :: xsum(:,:)
 integer :: nproc_space_comm
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
   if (nproc_space_comm /= 1) then
!    Accumulate xval on all proc. in spaceComm
     n1 =size(xval,dim=1)
     n2 =size(xval,dim=2)
     allocate(xsum(n1,n2), STAT=istat)
     if (istat /= 0) stop 'error allocating xsum in xsum_mpi_int2d'
     call MPI_ALLREDUCE(xval,xsum,n1*n2,MPI_INTEGER,&
&     MPI_SUM,spaceComm,ier)
     xval (:,:) = xsum(:,:)
     deallocate(xsum)
   end if
 end if
#endif
end subroutine xsum_mpi_int2d
!!***

!!****f* ABINIT/xsum_mpi_int3d
!! NAME
!!  xsum_mpi_int3d
!!
!! FUNCTION
!!  Combines values from all processes and distribute
!!  the result back to all processes.
!!  Target: three-dimensional integer arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
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
subroutine xsum_mpi_int3d(xval,spaceComm,ier)

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
 integer,intent(in) :: spaceComm
 integer,intent(out)   :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer ::   n1,n2,n3, istat
 integer , allocatable :: xsum(:,:,:)
 integer :: nproc_space_comm
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
!Accumulate xval on all proc. in spaceComm
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
   if (nproc_space_comm /= 1) then
     n1 =size(xval,dim=1)
     n2 =size(xval,dim=2)
     n3 =size(xval,dim=3)
     allocate(xsum(n1,n2,n3), STAT=istat)
     if (istat /= 0) stop 'error allocating xsum in xsum_mpi_int3d'
     call MPI_ALLREDUCE(xval,xsum,n1*n2*n3,MPI_INTEGER,&
&     MPI_SUM,spaceComm,ier)
     xval (:,:,:) = xsum(:,:,:)
     deallocate(xsum)
   end if
 end if
#endif
end subroutine xsum_mpi_int3d
!!***

!!****f* ABINIT/xsum_mpi_int4d
!! NAME
!!  xsum_mpi_int4d
!!
!! FUNCTION
!!  Combines values from all processes and distribute
!!  the result back to all processes.
!!  Target: four-diemnsional integer arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
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
subroutine xsum_mpi_int4d(xval,spaceComm,ier)

 use defs_basis

#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 integer,intent(inout) :: xval(:,:,:,:)
 integer,intent(in) :: spaceComm
 integer,intent(out)   :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer ::   n1,n2,n3,n4, istat
 integer , allocatable :: xsum(:,:,:,:)
 integer :: nproc_space_comm
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
!Accumulate xval on all proc. in spaceComm
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
   if (nproc_space_comm /= 1) then
     n1 =size(xval,dim=1)
     n2 =size(xval,dim=2)
     n3 =size(xval,dim=3)
     n4 =size(xval,dim=4)
     allocate(xsum(n1,n2,n3,n4), STAT=istat)
     if (istat /= 0) stop 'error allocating xsum in xsum_mpi_int4d'
     call MPI_ALLREDUCE(xval,xsum,n1*n2*n3*n4,MPI_INTEGER,&
&     MPI_SUM,spaceComm,ier)
     xval (:,:,:,:) = xsum(:,:,:,:)
     deallocate(xsum)
   end if
 end if
#endif
end subroutine xsum_mpi_int4d
!!***

!!****f* ABINIT/xsum_mpi_dp
!! NAME
!!  xsum_mpi_dp
!!
!! FUNCTION
!!  Combines values from all processes and distribute
!!  the result back to all processes.
!!  Target: one-dimensional double precision arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
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
subroutine xsum_mpi_dp(xval,spaceComm,ier)

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
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer ::   n1, istat
 real(dp) , allocatable :: xsum(:)
 integer :: nproc_space_comm
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
!  Accumulate xval on all proc. in spaceComm
   call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
   if (nproc_space_comm /= 1) then
     n1 = size(xval)
     allocate(xsum(n1), STAT=istat)
     if (istat /= 0) stop 'error allocating xsum in xsum_mpi_dp'
     call MPI_ALLREDUCE(xval,xsum,n1,MPI_DOUBLE_PRECISION,&
&     MPI_SUM,spaceComm,ier)
     xval (:) = xsum(:)
     deallocate(xsum)
   end if
 end if
#endif

end subroutine xsum_mpi_dp
!!***

!!****f* ABINIT/xsum_mpi_dpvt
!! NAME
!!  xsum_mpi_dpvt
!!
!! FUNCTION
!!  Combines values from all processes and distribute
!!  the result back to all processes.
!!  Target: scalar double precisions.
!!
!! INPUTS
!!  xval= buffer array
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  xsum= receive buffer
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  None
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xsum_mpi_dpvt(xval,xsum,spaceComm,ier)

 use defs_basis

#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(in) :: xval
 real(dp),intent(out) :: xsum
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: nproc_space_comm
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
!Accumulate xval on all proc. in spaceComm
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
   if (nproc_space_comm /= 1) then
     call MPI_ALLREDUCE(xval,xsum,1,MPI_DOUBLE_PRECISION,&
&     MPI_SUM,spaceComm,ier)
   end if
 else
#endif
   xsum=xval
#if defined HAVE_MPI
 end if
#endif
end subroutine xsum_mpi_dpvt
!!***

!!****f* ABINIT/xsum_mpi_dpv
!! NAME
!!  xsum_mpi_dpv
!!
!! FUNCTION
!!  Combines values from all processes and distribute
!!  the result back to all processes.
!!  Target: scalar double precisions.
!!
!! INPUTS
!!  spaceComm= MPI communicator
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
subroutine xsum_mpi_dpv(xval,spaceComm,ier)

 use defs_basis

#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(inout) :: xval
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: nproc_space_comm
 real(dp)  :: xsum
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
!Accumulate xval on all proc. in spaceComm
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
   if (nproc_space_comm /= 1) then
     call MPI_ALLREDUCE(xval,xsum,1,MPI_DOUBLE_PRECISION,&
&     MPI_SUM,spaceComm,ier)
     xval  = xsum
   end if
 end if
#endif
end subroutine xsum_mpi_dpv
!!***

!!****f* ABINIT/xsum_mpi_dpn
!! NAME
!!  xsum_mpi_dpn
!!
!! FUNCTION
!!  Combines values from all processes and distribute
!!  the result back to all processes.
!!  Target: one-dimensional double precision arrays.
!!
!! INPUTS
!!  n1= first dimension of the array
!!  spaceComm= MPI communicator
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
subroutine xsum_mpi_dpn(xval,n1,spaceComm,ier)

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
 integer ,intent(in)    :: n1
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: nproc_space_comm, istat
 real(dp) , allocatable :: xsum(:)
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
!  Accumulate xval on all proc. in spaceComm
   call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
   if (nproc_space_comm /= 1) then
     allocate(xsum(n1), STAT=istat)
     if (istat /= 0) stop 'error allocating xsum in xsum_mpi_dpn'
     call MPI_ALLREDUCE(xval,xsum,n1,MPI_DOUBLE_PRECISION,&
&     MPI_SUM,spaceComm,ier)
     xval (:) = xsum(:)
     deallocate(xsum)
   end if
 end if
#endif
end subroutine xsum_mpi_dpn
!!***

!!****f* ABINIT/xsum_mpi_dp2d
!! NAME
!!  xsum_mpi_dp2d
!!
!! FUNCTION
!!  Combines values from all processes and distribute
!!  the result back to all processes.
!!  Target: double precision two-dimensional arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
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
subroutine xsum_mpi_dp2d(xval,spaceComm,ier)

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
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n1,n2, istat
 real(dp) , allocatable :: xsum(:,:)
 integer :: nproc_space_comm
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
   if (nproc_space_comm /= 1) then
     n1 = size(xval,dim=1)
     n2 = size(xval,dim=2)
!    Accumulate xval on all proc. in spaceComm
     allocate(xsum(n1,n2), STAT=istat)
     if (istat /= 0) stop 'error allocating xsum in xsum_mpi_dp2d'
     call MPI_ALLREDUCE(xval,xsum,n1*n2,MPI_DOUBLE_PRECISION,&
&     MPI_SUM,spaceComm,ier)
     xval (:,:) = xsum(:,:)
     deallocate(xsum)
   end if
 end if
#endif
end subroutine xsum_mpi_dp2d
!!***

!!****f* ABINIT/xsum_mpi_dp3d
!! NAME
!!  xsum_mpi_dp3d
!!
!! FUNCTION
!!  Combines values from all processes and distribute
!!  the result back to all processes.
!!  Target: double precision three-dimensional arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
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
subroutine xsum_mpi_dp3d(xval,spaceComm,ier)

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
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n1,n2,n3, istat
 real(dp) , allocatable :: xsum(:,:,:)
 integer :: nproc_space_comm
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
   if (nproc_space_comm /= 1) then
     n1 = size(xval,dim=1)
     n2 = size(xval,dim=2)
     n3 = size(xval,dim=3)
!    Accumulate xval on all proc. in spaceComm
     allocate(xsum(n1,n2,n3), STAT=istat)
     if (istat /= 0) stop 'error allocating xsum in xsum_mpi_dp3d'
     call MPI_ALLREDUCE(xval,xsum,n1*n2*n3,MPI_DOUBLE_PRECISION,&
&     MPI_SUM,spaceComm,ier)
     xval (:,:,:) = xsum(:,:,:)
     deallocate(xsum)
   end if
 end if
#endif
end subroutine xsum_mpi_dp3d
!!***

!!****f* ABINIT/xsum_mpi_dp4d
!! NAME
!!  xsum_mpi_dp4d
!!
!! FUNCTION
!!  Combines values from all processes and distribute
!!  the result back to all processes.
!!  Target: double precision four-dimensional arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
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
subroutine xsum_mpi_dp4d(xval,spaceComm,ier)

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
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n1,n2,n3,n4, istat
 real(dp) , allocatable :: xsum(:,:,:,:)
 integer :: nproc_space_comm
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
   if (nproc_space_comm /= 1) then
     n1 = size(xval,dim=1)
     n2 = size(xval,dim=2)
     n3 = size(xval,dim=3)
     n4 = size(xval,dim=4)
!    Accumulate xval on all proc. in spaceComm
     allocate(xsum(n1,n2,n3,n4), STAT=istat)
     if (istat /= 0) stop 'error allocating xsum in xsum_mpi_dp4d'
     call MPI_ALLREDUCE(xval,xsum,n1*n2*n3*n4,MPI_DOUBLE_PRECISION,&
&     MPI_SUM,spaceComm,ier)
     xval (:,:,:,:) = xsum(:,:,:,:)
     deallocate(xsum)
   end if
 end if
#endif
end subroutine xsum_mpi_dp4d
!!***

!!****f* ABINIT/xsum_mpi_dp5d
!! NAME
!!  xsum_mpi_dp5d
!!
!! FUNCTION
!!  Combines values from all processes and distribute
!!  the result back to all processes.
!!  Target: double precision five-dimensional arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
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
subroutine xsum_mpi_dp5d(xval,spaceComm,ier)

 use defs_basis

#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(inout) :: xval(:,:,:,:,:)
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n1,n2,n3,n4,n5, istat
 real(dp) , allocatable :: xsum(:,:,:,:,:)
 integer :: nproc_space_comm
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
   if (nproc_space_comm /= 1) then
     n1 = size(xval,dim=1)
     n2 = size(xval,dim=2)
     n3 = size(xval,dim=3)
     n4 = size(xval,dim=4)
     n5 = size(xval,dim=5)
!    Accumulate xval on all proc. in spaceComm
     allocate(xsum(n1,n2,n3,n4,n5), STAT=istat)
     if (istat /= 0) stop 'error allocating xsum in xsum_mpi_dp5d'
     call MPI_ALLREDUCE(xval,xsum,n1*n2*n3*n4*n5,MPI_DOUBLE_PRECISION,&
&     MPI_SUM,spaceComm,ier)
     xval (:,:,:,:,:) = xsum(:,:,:,:,:)
     deallocate(xsum)
   end if
 end if
#endif
end subroutine xsum_mpi_dp5d
!!***

!!****f* ABINIT/xsum_mpi_dp6d
!! NAME
!!  xsum_mpi_dp6d
!!
!! FUNCTION
!!  Combines values from all processes and distribute
!!  the result back to all processes.
!!  Target: double precision six-dimensional arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
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
subroutine xsum_mpi_dp6d(xval,spaceComm,ier)

 use defs_basis

#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(inout) :: xval(:,:,:,:,:,:)
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n1,n2,n3,n4,n5,n6, istat
 real(dp) , allocatable :: xsum(:,:,:,:,:,:)
 integer :: nproc_space_comm
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
   if (nproc_space_comm /= 1) then
     n1 = size(xval,dim=1)
     n2 = size(xval,dim=2)
     n3 = size(xval,dim=3)
     n4 = size(xval,dim=4)
     n5 = size(xval,dim=5)
     n6 = size(xval,dim=6)
!    Accumulate xval on all proc. in spaceComm
     allocate(xsum(n1,n2,n3,n4,n5,n6), STAT=istat)
     if (istat /= 0) stop 'error allocating xsum in xsum_mpi_dp6d'
     call MPI_ALLREDUCE(xval,xsum,n1*n2*n3*n4*n5*n6,MPI_DOUBLE_PRECISION,&
&     MPI_SUM,spaceComm,ier)
     xval (:,:,:,:,:,:) = xsum(:,:,:,:,:,:)
     deallocate(xsum)
   end if
 end if
#endif
end subroutine xsum_mpi_dp6d
!!***

!!****f* ABINIT/xsum_mpi_dp2t
!! NAME
!!  xsum_mpi_dp2t
!!
!! FUNCTION
!!  Combines values from all processes and distribute
!!  the result back to all processes.
!!  Target: double precision one-dimensional array without transfers.
!!
!! INPUTS
!!  n1= first dimension of the array
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!  xsum= receive buffer
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xsum_mpi_dp2t(xval,xsum,n1,spaceComm,ier)

 use defs_basis

#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(inout) :: xval(:),xsum(:)
 integer ,intent(in)    :: n1
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
!  Accumulate xval on all proc. in spaceComm
   call MPI_ALLREDUCE(xval,xsum,n1,MPI_DOUBLE_PRECISION,MPI_SUM,spaceComm,ier)
 else
#endif
   xsum=xval
#if defined HAVE_MPI
 end if
#endif
end subroutine xsum_mpi_dp2t
!!***

!!****f* ABINIT/xsum_mpi_dp3d2t
!! NAME
!!  xsum_mpi_dp3d2t
!!
!! FUNCTION
!!  Combines values from all processes and distribute
!!  the result back to all processes.
!!  Target: double precision three-dimensional array without transfers.
!!
!! INPUTS
!!  n1= first dimension of the array
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!  xsum= receive buffer
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xsum_mpi_dp3d2t(xval,xsum,n1,spaceComm,ier)

 use defs_basis

#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(inout) :: xval(:,:,:),xsum(:,:,:)
 integer ,intent(in)    :: n1
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
!  Accumulate xval on all proc. in spaceComm
   call MPI_ALLREDUCE(xval,xsum,n1,MPI_DOUBLE_PRECISION,MPI_SUM,spaceComm,ier)
 else
#endif
   xsum=xval
#if defined HAVE_MPI
 end if
#endif
end subroutine xsum_mpi_dp3d2t
!!***

!!****f* ABINIT/xsum_mpi_dp4d2t
!! NAME
!!  xsum_mpi_dp4d2t
!!
!! FUNCTION
!!  Combines values from all processes and distribute
!!  the result back to all processes.
!!  Target: double precision four-dimensional array without transfers.
!!
!! INPUTS
!!  n1= first dimension of the array
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!  xsum= receive buffer
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xsum_mpi_dp4d2t(xval,xsum,n1,spaceComm,ier)

 use defs_basis

#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(inout) :: xval(:,:,:,:),xsum(:,:,:,:)
 integer ,intent(in)    :: n1
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
!  Accumulate xval on all proc. in spaceComm
   call MPI_ALLREDUCE(xval,xsum,n1,MPI_DOUBLE_PRECISION,MPI_SUM,spaceComm,ier)
 else
#endif
   xsum=xval
#if defined HAVE_MPI
 end if
#endif
end subroutine xsum_mpi_dp4d2t
!!***

!!****f* ABINIT/xsum_mpi_c0dc
!! NAME
!!  xsum_mpi_c0dc
!!
!! FUNCTION
!!  Combines values from all processes and distribute the result back to all processes.
!!  Target: double complex scalar
!!
!! INPUTS
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= scalar to be summed.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine xsum_mpi_c0dc(xval,spaceComm,ier)

 use defs_basis

#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 complex(dpc),intent(inout) :: xval
 integer,intent(in) :: spaceComm
 integer,intent(out)   :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: nproc_space_comm
 complex(dpc) :: xsum
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
   if (nproc_space_comm /= 1) then
!    Accumulate xval on all proc. in spaceComm
     call MPI_ALLREDUCE(xval,xsum,1,MPI_DOUBLE_COMPLEX,MPI_SUM,spaceComm,ier)
     xval = xsum
   end if
 end if
#endif

end subroutine xsum_mpi_c0dc
!!***


!!****f* ABINIT/xsum_mpi_c1dc
!! NAME
!!  xsum_mpi_c1dc
!!
!! FUNCTION
!!  Combines values from all processes and distribute
!!  the result back to all processes.
!!  Target: one-dimensional double complex arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
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
subroutine xsum_mpi_c1dc(xval,spaceComm,ier)

 use defs_basis

#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 complex(dpc),intent(inout) :: xval(:)
 integer,intent(in) :: spaceComm
 integer,intent(out)   :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer ::   n1
 complex(dpc) , allocatable :: xsum(:)
 integer :: nproc_space_comm, istat
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
   if (nproc_space_comm /= 1) then
!    Accumulate xval on all proc. in spaceComm
     n1 =size(xval,dim=1)
     allocate(xsum(n1), STAT=istat)
     if (istat /= 0) stop 'error allocating xsum in xsum_mpi_c1dc'
     call MPI_ALLREDUCE(xval,xsum,n1,MPI_DOUBLE_COMPLEX,MPI_SUM,spaceComm,ier)
     xval (:) = xsum(:)
     deallocate(xsum)
   end if
 end if
#endif

end subroutine xsum_mpi_c1dc
!!***

!!****f* ABINIT/xsum_mpi_c2dc
!! NAME
!!  xsum_mpi_c2dc
!!
!! FUNCTION
!!  Combines values from all processes and distribute
!!  the result back to all processes.
!!  Target: two-dimensional double complex arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
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
subroutine xsum_mpi_c2dc(xval,spaceComm,ier)

 use defs_basis

#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 complex(dpc),intent(inout) :: xval(:,:)
 integer,intent(in) :: spaceComm
 integer,intent(out)   :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer ::   n1,n2, istat
 integer :: nproc_space_comm
 complex(dpc),allocatable :: xsum(:,:)
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
   if (nproc_space_comm /= 1) then
!    Accumulate xval on all proc. in spaceComm
     n1 =size(xval,dim=1)
     n2 =size(xval,dim=2)
     allocate(xsum(n1,n2), STAT=istat)
     if (istat /= 0) stop 'error allocating xsum in xsum_mpi_c2dc'
     call MPI_ALLREDUCE(xval,xsum,n1*n2,MPI_DOUBLE_COMPLEX,MPI_SUM,spaceComm,ier)
     xval (:,:) = xsum(:,:)
     deallocate(xsum)
   end if
 end if
#endif

end subroutine xsum_mpi_c2dc
!!***

!!****f* ABINIT/xsum_mpi_c3dc
!! NAME
!!  xsum_mpi_c3dc
!!
!! FUNCTION
!!  Combines values from all processes and distribute
!!  the result back to all processes.
!!  Target: three-dimensional double complex arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
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
subroutine xsum_mpi_c3dc(xval,spaceComm,ier)

 use defs_basis

#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 complex(dpc),intent(inout) :: xval(:,:,:)
 integer,intent(in) :: spaceComm
 integer,intent(out)   :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer ::   n1,n2,n3, istat
 complex(dpc) , allocatable :: xsum(:,:,:)
 integer :: nproc_space_comm
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
   if (nproc_space_comm /= 1) then
!    Accumulate xval on all proc. in spaceComm
     n1 =size(xval,dim=1)
     n2 =size(xval,dim=2)
     n3 =size(xval,dim=3)
     allocate(xsum(n1,n2,n3), STAT=istat)
     if (istat /= 0) stop 'error allocating xsum in xsum_mpi_c3dc'
     call MPI_ALLREDUCE(xval,xsum,n1*n2*n3,MPI_DOUBLE_COMPLEX,&
&     MPI_SUM,spaceComm,ier)
     xval (:,:,:) = xsum(:,:,:)
     deallocate(xsum)
   end if
 end if
#endif
end subroutine xsum_mpi_c3dc
!!***

!!****f* ABINIT/xsum_mpi_c4dc
!! NAME
!!  xsum_mpi_c4dc
!!
!! FUNCTION
!!  Combines values from all processes and distribute
!!  the result back to all processes.
!!  Target: four-dimensional double complex arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
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
subroutine xsum_mpi_c4dc(xval,spaceComm,ier)

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
 integer,intent(in) :: spaceComm
 integer,intent(out)   :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n1,n2,n3,n4, istat
 complex(dpc),allocatable :: xsum(:,:,:,:)
 integer :: nproc_space_comm
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
   if (nproc_space_comm /= 1) then
!    Accumulate xval on all proc. in spaceComm
     n1 =size(xval,dim=1)
     n2 =size(xval,dim=2)
     n3 =size(xval,dim=3)
     n4 =size(xval,dim=4)
     allocate(xsum(n1,n2,n3,n4), STAT=istat)
     if (istat /= 0) stop 'error allocating xsum in xsum_mpi_c4dc'
     call MPI_ALLREDUCE(xval,xsum,n1*n2*n3*n4,MPI_DOUBLE_COMPLEX,MPI_SUM,spaceComm,ier)
     xval (:,:,:,:) = xsum(:,:,:,:)
     deallocate(xsum)
   end if
 end if
#endif

end subroutine xsum_mpi_c4dc
!!***

!!****f* ABINIT/xsum_mpi_c5dc
!! NAME
!!  xsum_mpi_c5dc
!!
!! FUNCTION
!!  Combines values from all processes and distribute the result back to all processes.
!!  Target: five-dimensional double precision complex arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
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
subroutine xsum_mpi_c5dc(xval,spaceComm,ier)

 use defs_basis

#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 complex(dpc),intent(inout) :: xval(:,:,:,:,:)
 integer,intent(in) :: spaceComm
 integer,intent(out)   :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n1,n2,n3,n4,n5, istat
 complex(dpc),allocatable :: xsum(:,:,:,:,:)
 integer :: nproc_space_comm
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
   if (nproc_space_comm /= 1) then !Accumulate xval on all proc. in spaceComm
     n1 =size(xval,dim=1)
     n2 =size(xval,dim=2)
     n3 =size(xval,dim=3)
     n4 =size(xval,dim=4)
     n5 =size(xval,dim=5)
     allocate(xsum(n1,n2,n3,n4,n5), STAT=istat)
     if (istat /= 0) stop 'error allocating xsum in xsum_mpi_c5dc'
     call MPI_ALLREDUCE(xval,xsum,n1*n2*n3*n4*n5,MPI_DOUBLE_COMPLEX,MPI_SUM,spaceComm,ier)
     xval (:,:,:,:,:) = xsum(:,:,:,:,:)
     deallocate(xsum)
   end if
 end if
#endif

end subroutine xsum_mpi_c5dc
!!***

!!****f* ABINIT/xsum_mpi_c6dc
!! NAME
!!  xsum_mpi_c6dc
!!
!! FUNCTION
!!  Combines values from all processes and distribute the result back to all processes.
!!  Target: six-dimensional double precision complex arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
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
subroutine xsum_mpi_c6dc(xval,spaceComm,ier)

 use defs_basis

#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 complex(dpc),intent(inout) :: xval(:,:,:,:,:,:)
 integer,intent(in) :: spaceComm
 integer,intent(out)   :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n1,n2,n3,n4,n5,n6,istat
 complex(dpc),allocatable :: xsum(:,:,:,:,:,:)
 integer :: nproc_space_comm
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
   if (nproc_space_comm /= 1) then !Accumulate xval on all proc. in spaceComm
     n1 =size(xval,dim=1)
     n2 =size(xval,dim=2)
     n3 =size(xval,dim=3)
     n4 =size(xval,dim=4)
     n5 =size(xval,dim=5)
     n6 =size(xval,dim=6)
     allocate(xsum(n1,n2,n3,n4,n5,n6), STAT=istat)
     if (istat /= 0) stop 'error allocating xsum in xsum_mpi_c6dc'
     call MPI_ALLREDUCE(xval,xsum,n1*n2*n3*n4*n5*n6,MPI_DOUBLE_COMPLEX,MPI_SUM,spaceComm,ier)
     xval = xsum
     deallocate(xsum)
   end if
 end if
#endif

end subroutine xsum_mpi_c6dc
!!***

!!****f* ABINIT/xsum_mpi_c1cplx
!! NAME
!!  xsum_mpi_c1cplx
!!
!! FUNCTION
!!  Combines values from all processes and distribute
!!  the result back to all processes.
!!  Target: one-dimensional complex arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
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
subroutine xsum_mpi_c1cplx(xval,spaceComm,ier)

 use defs_basis

#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments----------------
 complex,intent(inout) :: xval(:)
 integer,intent(in) :: spaceComm
 integer,intent(out)   :: ier

!Local variables--------------
#if defined HAVE_MPI
 integer ::   n1
 complex, allocatable :: xsum(:)
 integer :: nproc_space_comm, istat
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
   if (nproc_space_comm /= 1) then
!    Accumulate xval on all proc. in spaceComm
     n1 =size(xval,dim=1)
     allocate(xsum(n1), STAT=istat)
     if (istat /= 0) stop 'error allocating xsum in xsum_mpi_c1cplx'
     call MPI_ALLREDUCE(xval,xsum,n1,MPI_COMPLEX,MPI_SUM,spaceComm,ier)
     xval (:) = xsum(:)
     deallocate(xsum)
   end if
 end if
#endif
end subroutine xsum_mpi_c1cplx
!!***

!!****f* ABINIT/xsum_mpi_c2cplx
!! NAME
!!  xsum_mpi_c2cplx
!!
!! FUNCTION
!!  Combines values from all processes and distribute
!!  the result back to all processes.
!!  Target: two-dimensional complex arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
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
subroutine xsum_mpi_c2cplx(xval,spaceComm,ier)

 use defs_basis

#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments----------------
 complex,intent(inout) :: xval(:,:)
 integer,intent(in) :: spaceComm
 integer,intent(out)   :: ier

!Local variables--------------
#if defined HAVE_MPI
 integer ::   n1,n2, istat
 complex, allocatable :: xsum(:,:)
 integer :: nproc_space_comm
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
   if (nproc_space_comm /= 1) then
!    Accumulate xval on all proc. in spaceComm
     n1 =size(xval,dim=1)
     n2 =size(xval,dim=2)
     allocate(xsum(n1,n2), STAT=istat)
     if (istat /= 0) stop 'error allocating xsum in xsum_mpi_c2cplx'
     call MPI_ALLREDUCE(xval,xsum,n1*n2,MPI_COMPLEX,&
&     MPI_SUM,spaceComm,ier)
     xval (:,:) = xsum(:,:)
     deallocate(xsum)
   end if
 end if
#endif
end subroutine xsum_mpi_c2cplx
!!***

!!****f* ABINIT/xsum_mpi_c3cplx
!! NAME
!!  xsum_mpi_c3cplx
!!
!! FUNCTION
!!  Combines values from all processes and distribute
!!  the result back to all processes.
!!  Target: three-dimensional complex arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
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
subroutine xsum_mpi_c3cplx(xval,spaceComm,ier)

 use defs_basis

#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments----------------
 complex,intent(inout) :: xval(:,:,:)
 integer,intent(in) :: spaceComm
 integer,intent(out)   :: ier

!Local variables--------------
#if defined HAVE_MPI
 integer ::   n1,n2,n3, istat
 complex , allocatable :: xsum(:,:,:)
 integer :: nproc_space_comm
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
   if (nproc_space_comm /= 1) then
!    Accumulate xval on all proc. in spaceComm
     n1 =size(xval,dim=1)
     n2 =size(xval,dim=2)
     n3 =size(xval,dim=3)
     allocate(xsum(n1,n2,n3), STAT=istat)
     if (istat /= 0) stop 'error allocating xsum in xsum_mpi_c3cplx'
     call MPI_ALLREDUCE(xval,xsum,n1*n2*n3,MPI_COMPLEX,&
&     MPI_SUM,spaceComm,ier)
     xval (:,:,:) = xsum(:,:,:)
     deallocate(xsum)
   end if
 end if
#endif
end subroutine xsum_mpi_c3cplx
!!***

!!****f* ABINIT/xsum_mpi_c4cplx
!! NAME
!!  xsum_mpi_c4cplx
!!
!! FUNCTION
!!  Combines values from all processes and distribute
!!  the result back to all processes.
!!  Target: four-dimensional complex arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
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
subroutine xsum_mpi_c4cplx(xval,spaceComm,ier)

 use defs_basis

#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments----------------
 complex,intent(inout) :: xval(:,:,:,:)
 integer,intent(in) :: spaceComm
 integer,intent(out)   :: ier

!Local variables--------------
#if defined HAVE_MPI
 integer ::   n1,n2,n3,n4, istat
 complex , allocatable :: xsum(:,:,:,:)
 integer :: nproc_space_comm
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
   if (nproc_space_comm /= 1) then
!    Accumulate xval on all proc. in spaceComm
     n1 =size(xval,dim=1)
     n2 =size(xval,dim=2)
     n3 =size(xval,dim=3)
     n4 =size(xval,dim=4)
     allocate(xsum(n1,n2,n3,n4), STAT=istat)
     if (istat /= 0) stop 'error allocating xsum in xsum_mpi_c4cplx'
     call MPI_ALLREDUCE(xval,xsum,n1*n2*n3*n4,MPI_COMPLEX,&
&     MPI_SUM,spaceComm,ier)
     xval (:,:,:,:) = xsum(:,:,:,:)
     deallocate(xsum)
   end if
 end if
#endif
end subroutine xsum_mpi_c4cplx
!!***

!!****f* ABINIT/xsum_mpi_c5cplx
!! NAME
!!  xsum_mpi_c5cplx
!!
!! FUNCTION
!!  Combines values from all processes and distribute
!!  the result back to all processes.
!!  Target: five-dimensional complex arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
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
subroutine xsum_mpi_c5cplx(xval,spaceComm,ier)

 use defs_basis

#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments----------------
 complex,intent(inout) :: xval(:,:,:,:,:)
 integer,intent(in) :: spaceComm
 integer,intent(out)   :: ier

!Local variables--------------
#if defined HAVE_MPI
 integer ::   n1,n2,n3,n4,n5, istat
 complex , allocatable :: xsum(:,:,:,:,:)
 integer :: nproc_space_comm
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
   if (nproc_space_comm /= 1) then
!    Accumulate xval on all proc. in spaceComm
     n1 =size(xval,dim=1)
     n2 =size(xval,dim=2)
     n3 =size(xval,dim=3)
     n4 =size(xval,dim=4)
     n5 =size(xval,dim=5)
     allocate(xsum(n1,n2,n3,n4,n5), STAT=istat)
     if (istat /= 0) stop 'error allocating xsum in xsum_mpi_c5cplx'
     call MPI_ALLREDUCE(xval,xsum,n1*n2*n3*n4*n5,MPI_COMPLEX,MPI_SUM,spaceComm,ier)
     xval (:,:,:,:,:) = xsum(:,:,:,:,:)
     deallocate(xsum)
   end if
 end if
#endif

end subroutine xsum_mpi_c5cplx
!!***

!!****f* ABINIT/xsum_mpi_c6cplx
!! NAME
!!  xsum_mpi_c6cplx
!!
!! FUNCTION
!!  Combines values from all processes and distribute the result back to all processes.
!!  Target: six-dimensional complex arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
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

subroutine xsum_mpi_c6cplx(xval,spaceComm,ier)

 use defs_basis

#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments----------------
 complex,intent(inout) :: xval(:,:,:,:,:,:)
 integer,intent(in) :: spaceComm
 integer,intent(out)   :: ier

!Local variables--------------
#if defined HAVE_MPI
 integer ::   n1,n2,n3,n4,n5,n6,istat
 complex,allocatable :: xsum(:,:,:,:,:,:)
 integer :: nproc_space_comm
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_COMM_SIZE(spaceComm,nproc_space_comm,ier)
   if (nproc_space_comm /= 1) then
!    Accumulate xval on all proc. in spaceComm
     n1 =size(xval,dim=1)
     n2 =size(xval,dim=2)
     n3 =size(xval,dim=3)
     n4 =size(xval,dim=4)
     n5 =size(xval,dim=5)
     n6 =size(xval,dim=6)
     allocate(xsum(n1,n2,n3,n4,n5,n6), STAT=istat)
     if (istat /= 0) stop 'error allocating xsum in xsum_mpi_c6cplx'
     call MPI_ALLREDUCE(xval,xsum,n1*n2*n3*n4*n5*n6,MPI_COMPLEX,MPI_SUM,spaceComm,ier)
     xval = xsum
     deallocate(xsum)
   end if
 end if
#endif

end subroutine xsum_mpi_c6cplx
!!***

!!****f* ABINIT/xsum_mpi_log1d
!! NAME
!!  xsum_mpi_log1d
!!
!! FUNCTION
!!  Combines values from all processes and distribute
!!  the result back to all processes.
!!  Target: one-dimensional logical arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
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
subroutine xsum_mpi_log1d(xval,spaceComm,ier)

 use defs_basis

#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 integer,intent(in) :: spaceComm
 integer,intent(out)   :: ier
 logical,intent(inout) :: xval(:)

!Local variables-------------------------------
#if defined HAVE_MPI
 integer :: n1, istat
 logical , allocatable :: xsum(:)
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
!  Accumulate xval on all proc. in spaceComm
   n1 = size(xval)
   allocate(xsum(n1), STAT=istat)
   if (istat /= 0) stop 'error allocating xsum in xsum_mpi_log1d'
   call MPI_ALLREDUCE(xval,xsum,n1,MPI_LOGICAL,&
&   MPI_LOR,spaceComm,ier)
   xval (:) = xsum(:)
   deallocate(xsum)
 end if
#endif
end subroutine xsum_mpi_log1d

!!***


!!****f* ABINIT/xsum_mpi_log2d
!! NAME
!!  xsum_mpi_log2d
!!
!! FUNCTION
!!  Combines values from all processes and distribute
!!  the result back to all processes.
!!  Target: two-dimensional logical arrays.
!!
!! INPUTS
!!  spaceComm= MPI communicator
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
subroutine xsum_mpi_log2d(xval,spaceComm,ier)

 use defs_basis

#if defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 integer,intent(in) :: spaceComm
 integer,intent(out)   :: ier
 logical,intent(inout) :: xval(:,:)

!Local variables-------------------------------
#if defined HAVE_MPI
 integer :: n1,n2, istat
 logical , allocatable :: xsum(:,:)
#endif

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
!  Accumulate xval on all proc. in spaceComm
   n1 = size(xval,1)
   n2 = size(xval,2)
   allocate(xsum(n1,n2), STAT=istat)
   if (istat /= 0) stop 'error allocating xsum in xsum_mpi_log2d'
   call MPI_ALLREDUCE(xval,xsum,n1*n2,MPI_LOGICAL,&
&   MPI_LOR,spaceComm,ier)
   xval (:,:) = xsum(:,:)
   deallocate(xsum)
 end if
#endif
end subroutine xsum_mpi_log2d

!!***
