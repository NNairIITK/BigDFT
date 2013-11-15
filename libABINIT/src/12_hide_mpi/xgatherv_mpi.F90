!{\src2tex{textfont=tt}}
!!****f* ABINIT/xgatherv_mpi
!! NAME
!!  xgatherv_mpi
!!
!! FUNCTION
!!  This module contains functions that calls MPI routine,
!!  if we compile the code using the MPI CPP flags.
!!  xgatherv_mpi is the generic function.
!!
!! COPYRIGHT
!!  Copyright (C) 2001-2013 ABINIT group (MT,GG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~ABINIT/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

!!***

!!****f* ABINIT/xgatherv_mpi_int
!! NAME
!!  xgatherv_mpi_int
!!
!! FUNCTION
!!  Gathers data from all tasks and delivers it to all.
!!  Target: one-dimensional integer arrays.
!!
!! INPUTS
!!  xval= buffer array
!!  recvcounts= number of received elements
!!  displs= relative offsets for incoming data
!!  nelem= number of elements
!!  root= rank of receiving process
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  recvbuf= received buffer
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_gatherv
!!
!! SOURCE
subroutine xgatherv_mpi_int(xval,nelem,recvbuf,recvcounts,displs,root,spaceComm,ier)


 use defs_basis
#if defined HAVE_MPI && defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgatherv_mpi_int'
!End of the abilint section

 implicit none

#if defined HAVE_MPI && defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 integer,intent(in) :: xval(:)
 integer,intent(inout)   :: recvbuf(:)
 integer,intent(in) :: recvcounts(:),displs(:)
 integer,intent(in) :: nelem,root,spaceComm
 integer,intent(out)   :: ier

!Local variables-------------------
 integer :: cc,dd

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_gatherV(xval,nelem,MPI_INTEGER,recvbuf,recvcounts,displs,&
&   MPI_INTEGER,root,spaceComm,ier)
 else if (spaceComm == MPI_COMM_SELF) then
#endif
   dd=0;if (size(displs)>0) dd=displs(1)
   cc=size(xval);if (size(recvcounts)>0) cc=recvcounts(1)
   recvbuf(dd+1:dd+cc)=xval(1:cc)
#if defined HAVE_MPI
 end if
#endif
end subroutine xgatherv_mpi_int
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/xgatherv_mpi_int2d
!! NAME
!!  xgatherv_mpi_int2d
!!
!! FUNCTION
!!  This module contains functions that calls MPI routine,
!!  if we compile the code using the MPI CPP flags.
!!  xgatherv_mpi is the generic function.
!!
!! INPUTS
!!  xval= buffer array
!!  recvcounts= number of received elements
!!  displs= relative offsets for incoming data
!!  nelem= number of elements
!!  root= rank of receiving process
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  recvbuf= received buffer
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_gatherv
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine xgatherv_mpi_int2d(xval,nelem,recvbuf,recvcounts,displs,root,spaceComm,ier)


 use defs_basis
#if defined HAVE_MPI && defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgatherv_mpi_int2d'
!End of the abilint section

 implicit none

#if defined HAVE_MPI && defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 integer,intent(in) :: xval(:,:)
 integer,intent(inout) :: recvbuf(:,:)
 integer,intent(in) :: recvcounts(:),displs(:)
 integer,intent(in) :: nelem,root,spaceComm
 integer,intent(out)   :: ier

!Local variables--------------
 integer :: cc,dd,sz1

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_gatherV(xval,nelem,MPI_INTEGER,recvbuf,recvcounts,displs,&
&   MPI_INTEGER,root,spaceComm,ier)
 else if (spaceComm == MPI_COMM_SELF) then
#endif
   sz1=size(xval,1)
   dd=0;if (size(displs)>0) dd=displs(1)/sz1
   cc=size(xval,2);if (size(recvcounts)>0) cc=recvcounts(1)/sz1
   recvbuf(:,dd+1:dd+cc)=xval(:,1:cc)
#if defined HAVE_MPI
 end if
#endif
end subroutine xgatherv_mpi_int2d
!!***

!!****f* ABINIT/xgatherv_mpi_dp
!! NAME
!!  xgatherv_mpi_dp
!!
!! FUNCTION
!!  Gathers data from all tasks and delivers it to all.
!!  Target: one-dimensional double precision arrays.
!!
!! INPUTS
!!  xval= buffer array
!!  recvcounts= number of received elements
!!  displs= relative offsets for incoming data
!!  nelem= number of elements
!!  root= rank of receiving process
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  recvbuf= received buffer
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_gatherv
!!
!! SOURCE
subroutine xgatherv_mpi_dp(xval,nelem,recvbuf,recvcounts,displs,root,spaceComm,ier)


 use defs_basis
#if defined HAVE_MPI && defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgatherv_mpi_dp'
!End of the abilint section

 implicit none

#if defined HAVE_MPI && defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(in) :: xval(:)
 real(dp),intent(inout)   :: recvbuf(:)
 integer,intent(in) :: recvcounts(:),displs(:)
 integer,intent(in) :: nelem,root,spaceComm
 integer,intent(out)   :: ier

!Local variables--------------
 integer :: cc,dd

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_gatherV(xval,nelem,MPI_DOUBLE_PRECISION,recvbuf,recvcounts,displs,&
&   MPI_DOUBLE_PRECISION,root,spaceComm,ier)
 else if (spaceComm == MPI_COMM_SELF) then
#endif
   dd=0;if (size(displs)>0) dd=displs(1)
   cc=size(xval);if (size(recvcounts)>0) cc=recvcounts(1)
   recvbuf(dd+1:dd+cc)=xval(1:cc)
#if defined HAVE_MPI
 end if
#endif
end subroutine xgatherv_mpi_dp
!!***

!!****f* ABINIT/xgatherv_mpi_dp2d
!! NAME
!!  xgatherv_mpi_dp2d
!!
!! FUNCTION
!!  Gathers data from all tasks and delivers it to all.
!!  Target: double precision two-dimensional arrays.
!!
!! INPUTS
!!  xval= buffer array
!!  recvcounts= number of received elements
!!  displs= relative offsets for incoming data
!!  nelem= number of elements
!!  root= rank of receiving process
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  recvbuf= received buffer
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_gatherv
!!
!! SOURCE
subroutine xgatherv_mpi_dp2d(xval,nelem,recvbuf,recvcounts,displs,root,spaceComm,ier)


 use defs_basis
#if defined HAVE_MPI && defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgatherv_mpi_dp2d'
!End of the abilint section

 implicit none

#if defined HAVE_MPI && defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(in) :: xval(:,:)
 real(dp),intent(inout) :: recvbuf(:,:)
 integer,intent(in) :: recvcounts(:),displs(:)
 integer,intent(in) :: nelem,root,spaceComm
 integer,intent(out)   :: ier

!Local variables--------------
 integer :: cc,dd,sz1

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_gatherV(xval,nelem,MPI_DOUBLE_PRECISION,recvbuf,recvcounts,displs,&
&   MPI_DOUBLE_PRECISION,root,spaceComm,ier)
 else if (spaceComm == MPI_COMM_SELF) then
#endif
   sz1=size(xval,1)
   dd=0;if (size(displs)>0) dd=displs(1)/sz1
   cc=size(xval,2);if (size(recvcounts)>0) cc=recvcounts(1)/sz1
   recvbuf(:,dd+1:dd+cc)=xval(:,1:cc)
#if defined HAVE_MPI
 end if
#endif
end subroutine xgatherv_mpi_dp2d
!!***

!!****f* ABINIT/xgatherv_mpi_dp3d
!! NAME
!!  xgatherv_mpi_dp3d
!!
!! FUNCTION
!!  Gathers data from all tasks and delivers it to all.
!!  Target: double precision three-dimensional arrays.
!!
!! INPUTS
!!  xval= buffer array
!!  recvcounts= number of received elements
!!  displs= relative offsets for incoming data
!!  nelem= number of elements
!!  root= rank of receiving process
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  recvbuf= received buffer
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_gatherv
!!
!! SOURCE
subroutine xgatherv_mpi_dp3d(xval,nelem,recvbuf,recvcounts,displs,root,spaceComm,ier)


 use defs_basis
#if defined HAVE_MPI && defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgatherv_mpi_dp3d'
!End of the abilint section

 implicit none

#if defined HAVE_MPI && defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(in) :: xval(:,:,:)
 real(dp),intent(inout)   :: recvbuf(:,:,:)
 integer,intent(in) :: recvcounts(:),displs(:)
 integer,intent(in) :: nelem,root,spaceComm
 integer,intent(out)   :: ier

!Local variables--------------
 integer :: cc,dd,sz12

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_gatherV(xval,nelem,MPI_DOUBLE_PRECISION,recvbuf,recvcounts,displs,&
&   MPI_DOUBLE_PRECISION,root,spaceComm,ier)
 else if (spaceComm == MPI_COMM_SELF) then
#endif
   sz12=size(xval,1)*size(xval,2)
   dd=0;if (size(displs)>0) dd=displs(1)/sz12
   cc=size(xval,3);if (size(recvcounts)>0) cc=recvcounts(1)/sz12
   recvbuf(:,:,dd+1:dd+cc)=xval(:,:,1:cc)
#if defined HAVE_MPI
 end if
#endif
end subroutine xgatherv_mpi_dp3d
!!***

!!****f* ABINIT/xgatherv_mpi_dp4d
!! NAME
!!  xgatherv_mpi_dp4d
!!
!! FUNCTION
!!  Gathers data from all tasks and delivers it to all.
!!  Target: double precision four-dimensional arrays.
!!
!! INPUTS
!!  xval= buffer array
!!  recvcounts= number of received elements
!!  displs= relative offsets for incoming data
!!  nelem= number of elements
!!  root= rank of receiving process
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  recvbuf= received buffer
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_gatherv
!!
!! SOURCE
subroutine xgatherv_mpi_dp4d(xval,nelem,recvbuf,recvcounts,displs,root,spaceComm,ier)


 use defs_basis
#if defined HAVE_MPI && defined HAVE_MPI2 && ! defined HAVE_MPI_INCLUDED_ONCE
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgatherv_mpi_dp4d'
!End of the abilint section

 implicit none

#if defined HAVE_MPI && defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(in) :: xval(:,:,:,:)
 real(dp),intent(inout)   :: recvbuf(:,:,:,:)
 integer,intent(in) :: recvcounts(:),displs(:)
 integer,intent(in) :: nelem,root,spaceComm
 integer,intent(out)   :: ier

!Local variables-------------------
 integer :: cc,dd,sz123

! *************************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
   call MPI_gatherV(xval,nelem,MPI_DOUBLE_PRECISION,recvbuf,recvcounts,displs,&
&   MPI_DOUBLE_PRECISION,root,spaceComm,ier)
 else if (spaceComm == MPI_COMM_SELF) then
#endif
   sz123=size(xval,1)*size(xval,2)*size(xval,3)
   dd=0;if (size(displs)>0) dd=displs(1)/sz123
   cc=size(xval,4);if (size(recvcounts)>0) cc=recvcounts(1)/sz123
   recvbuf(:,:,:,dd+1:dd+cc)=xval(:,:,:,1:cc)
#if defined HAVE_MPI
 end if
#endif
end subroutine xgatherv_mpi_dp4d
!!***
