!!****m* ABINIT/interfaces_12_hide_mpi
!! NAME
!! interfaces_12_hide_mpi
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/12_hide_mpi
!!
!! COPYRIGHT
!! Copyright (C) 2010-2011 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! THIS FILE IS GENERATED AUTOMATICALLY BY abilint.
!! To do that: config/scripts/abilint . .
!! 
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module interfaces_12_hide_mpi

 implicit none

interface
 subroutine xallgather_mpi_int(xval,recvcounts,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: spaceComm
  integer,intent(inout) :: xval
  integer,intent(inout) :: recvcounts(:)
 end subroutine xallgather_mpi_int
end interface

interface
 subroutine xallgather_mpi_char(charval,recvcounts,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: spaceComm
  character(20),intent(inout) :: charval
  character(20),intent(inout) :: recvcounts(:)
 end subroutine xallgather_mpi_char
end interface

interface
 subroutine xallgather_mpi_int1d(xval,nelem,recvcounts,spaceComm,ier)
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: nelem
  integer ,intent(in) :: spaceComm
  integer,intent(inout) :: recvcounts(:)
  integer,intent(in) :: xval(:)
 end subroutine xallgather_mpi_int1d
end interface

interface
 subroutine xallgather_mpi_dp1d(xval,nelem,recvcounts,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: nelem
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: recvcounts(:)
  real(dp),intent(in) :: xval(:)
 end subroutine xallgather_mpi_dp1d
end interface

interface
 subroutine xallgather_mpi_dp2d(xval,nelem,recvcounts,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: nelem
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: recvcounts(:,:)
  real(dp),intent(in) :: xval(:,:)
 end subroutine xallgather_mpi_dp2d
end interface

interface
 subroutine xallgather_mpi_dp3d(xval,nelem,recvcounts,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: nelem
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: recvcounts(:,:,:)
  real(dp),intent(in) :: xval(:,:,:)
 end subroutine xallgather_mpi_dp3d
end interface

interface
 subroutine xallgather_mpi_dp4d(xval,nelem,recvcounts,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: nelem
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: recvcounts(:,:,:,:)
  real(dp),intent(in) :: xval(:,:,:,:)
 end subroutine xallgather_mpi_dp4d
end interface

interface
 subroutine xallgatherv_mpi_int2d(xval,nelem,recvbuf,recvcounts,displs,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: nelem
  integer,intent(in) :: spaceComm
  integer,intent(in) :: displs(:)
  integer,intent(inout) :: recvbuf(:,:)
  integer,intent(in) :: recvcounts(:)
  integer,intent(in) :: xval(:,:)
 end subroutine xallgatherv_mpi_int2d
end interface

interface
 subroutine xallgatherv_mpi_int(xval,nelem,recvbuf,recvcounts,displs,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: nelem
  integer,intent(in) :: spaceComm
  integer,intent(in) :: displs(:)
  integer,intent(inout) :: recvbuf(:)
  integer,intent(in) :: recvcounts(:)
  integer,intent(in) :: xval(:)
 end subroutine xallgatherv_mpi_int
end interface

interface
 subroutine xallgatherv_mpi_dp(xval,nelem,recvbuf,recvcounts,displs,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: nelem
  integer,intent(in) :: spaceComm
  integer,intent(in) :: displs(:)
  integer,intent(in) :: recvcounts(:)
  real(dp),intent(inout) :: recvbuf(:)
  real(dp),intent(in) :: xval(:)
 end subroutine xallgatherv_mpi_dp
end interface

interface
 subroutine xallgatherv_mpi_dp2d(xval,nelem,recvbuf,recvcounts,displs,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: nelem
  integer,intent(in) :: spaceComm
  integer,intent(in) :: displs(:)
  integer,intent(in) :: recvcounts(:)
  real(dp),intent(inout) :: recvbuf(:,:)
  real(dp),intent(in) :: xval(:,:)
 end subroutine xallgatherv_mpi_dp2d
end interface

interface
 subroutine xallgatherv_mpi_dp3d(xval,nelem,recvbuf,recvcounts,displs,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: nelem
  integer,intent(in) :: spaceComm
  integer,intent(in) :: displs(:)
  integer,intent(in) :: recvcounts(:)
  real(dp),intent(inout) :: recvbuf(:,:,:)
  real(dp),intent(in) :: xval(:,:,:)
 end subroutine xallgatherv_mpi_dp3d
end interface

interface
 subroutine xallgatherv_mpi_dp4d(xval,nelem,recvbuf,recvcounts,displs,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: nelem
  integer,intent(in) :: spaceComm
  integer,intent(in) :: displs(:)
  integer,intent(in) :: recvcounts(:)
  real(dp),intent(inout) :: recvbuf(:,:,:,:)
  real(dp),intent(in) :: xval(:,:,:,:)
 end subroutine xallgatherv_mpi_dp4d
end interface

interface
 subroutine xalltoall_mpi_dp2d(xval, sendsize, recvbuf, recvsize, spaceComm, ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: recvsize
  integer ,intent(in) :: sendsize
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: recvbuf(:,:)
  real(dp),intent(in) :: xval(:,:)
 end subroutine xalltoall_mpi_dp2d
end interface

interface
 subroutine xalltoallv_mpi_dp2d(xval,sendcnts,sdispls,recvbuf,recvcnts,rdispls,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: spaceComm
  integer ,intent(in) :: rdispls(:)
  integer ,intent(in) :: recvcnts(:)
  integer ,intent(in) :: sdispls(:)
  integer ,intent(in) :: sendcnts(:)
  real(dp),intent(inout) :: recvbuf(:,:)
  real(dp),intent(in) :: xval(:,:)
 end subroutine xalltoallv_mpi_dp2d
end interface

interface
 subroutine xalltoallv_mpi_int2d(xval,sendcnts,sdispls,recvbuf,recvcnts,rdispls,spaceComm,ier)
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: spaceComm
  integer ,intent(in) :: rdispls(:)
  integer,intent(inout) :: recvbuf(:,:)
  integer ,intent(in) :: recvcnts(:)
  integer ,intent(in) :: sdispls(:)
  integer ,intent(in) :: sendcnts(:)
  integer,intent(in) :: xval(:,:)
 end subroutine xalltoallv_mpi_int2d
end interface

interface
 subroutine xalltoallv_mpi_dp1d(xval,sendcnts,sdispls,recvbuf,recvcnts,rdispls,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: rdispls
  integer ,intent(in) :: spaceComm
  integer ,intent(in) :: recvcnts(:)
  integer ,intent(in) :: sdispls(:)
  integer ,intent(in) :: sendcnts(:)
  real(dp),intent(inout) :: recvbuf(:)
  real(dp),intent(in) :: xval(:)
 end subroutine xalltoallv_mpi_dp1d
end interface

interface
 subroutine xalltoallv_mpi_dp1d2(xval,sendcnts,sdispls,recvbuf,recvcnts,rdispls,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: spaceComm
  integer ,intent(in) :: rdispls(:)
  integer ,intent(in) :: recvcnts(:)
  integer ,intent(in) :: sdispls(:)
  integer ,intent(in) :: sendcnts(:)
  real(dp),intent(inout) :: recvbuf(:)
  real(dp),intent(in) :: xval(:)
 end subroutine xalltoallv_mpi_dp1d2
end interface

interface
 subroutine xcast_mpi_intv(xval,master,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: master
  integer,intent(in) :: spaceComm
  integer,intent(inout) :: xval
 end subroutine xcast_mpi_intv
end interface

interface
 subroutine xcast_mpi_int1d(xval,master,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: master
  integer,intent(in) :: spaceComm
  integer,intent(inout) :: xval(:)
 end subroutine xcast_mpi_int1d
end interface

interface
 subroutine xcast_mpi_int2d(xval,master,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: master
  integer,intent(in) :: spaceComm
  integer,intent(inout) :: xval(:,:)
 end subroutine xcast_mpi_int2d
end interface

interface
 subroutine xcast_mpi_int3d(xval,master,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: master
  integer,intent(in) :: spaceComm
  integer,intent(inout) :: xval(:,:,:)
 end subroutine xcast_mpi_int3d
end interface

interface
 subroutine xcast_mpi_dpv(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval
 end subroutine xcast_mpi_dpv
end interface

interface
 subroutine xcast_mpi_dp1d(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval(:)
 end subroutine xcast_mpi_dp1d
end interface

interface
 subroutine xcast_mpi_dp2d(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval(:,:)
 end subroutine xcast_mpi_dp2d
end interface

interface
 subroutine xcast_mpi_dp3d(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval(:,:,:)
 end subroutine xcast_mpi_dp3d
end interface

interface
 subroutine xcast_mpi_dp4d(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval(:,:,:,:)
 end subroutine xcast_mpi_dp4d
end interface

interface
 subroutine xcast_mpi_spv(xval,master,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: master
  integer,intent(in) :: spaceComm
  real,intent(inout) :: xval
 end subroutine xcast_mpi_spv
end interface

interface
 subroutine xcast_mpi_sp1d(xval,master,spaceComm,ier)
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  real,intent(inout) :: xval(:)
 end subroutine xcast_mpi_sp1d
end interface

interface
 subroutine xcast_mpi_sp2d(xval,master,spaceComm,ier)
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  real,intent(inout) :: xval(:,:)
 end subroutine xcast_mpi_sp2d
end interface

interface
 subroutine xcast_mpi_sp3d(xval,master,spaceComm,ier)
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  real,intent(inout) :: xval(:,:,:)
 end subroutine xcast_mpi_sp3d
end interface

interface
 subroutine xcast_mpi_sp4d(xval,master,spaceComm,ier)
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  real,intent(inout) :: xval(:,:,:,:)
 end subroutine xcast_mpi_sp4d
end interface

interface
 subroutine xcast_mpi_cplxv(xval,master,spaceComm,ier)
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  complex,intent(inout) :: xval
 end subroutine xcast_mpi_cplxv
end interface

interface
 subroutine xcast_mpi_cplx1d(xval,master,spaceComm,ier)
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  complex,intent(inout) :: xval(:)
 end subroutine xcast_mpi_cplx1d
end interface

interface
 subroutine xcast_mpi_cplx2d(xval,master,spaceComm,ier)
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  complex,intent(inout) :: xval(:,:)
 end subroutine xcast_mpi_cplx2d
end interface

interface
 subroutine xcast_mpi_cplx3d(xval,master,spaceComm,ier)
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  complex,intent(inout) :: xval(:,:,:)
 end subroutine xcast_mpi_cplx3d
end interface

interface
 subroutine xcast_mpi_cplx4d(xval,master,spaceComm,ier)
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  complex,intent(inout) :: xval(:,:,:,:)
 end subroutine xcast_mpi_cplx4d
end interface

interface
 subroutine xcast_mpi_dcv(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  complex(dpc),intent(inout) :: xval
 end subroutine xcast_mpi_dcv
end interface

interface
 subroutine xcast_mpi_dc1d(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  complex(dpc),intent(inout) :: xval(:)
 end subroutine xcast_mpi_dc1d
end interface

interface
 subroutine xcast_mpi_dc2d(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  complex(dpc),intent(inout) :: xval(:,:)
 end subroutine xcast_mpi_dc2d
end interface

interface
 subroutine xcast_mpi_dc3d(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  complex(dpc),intent(inout) :: xval(:,:,:)
 end subroutine xcast_mpi_dc3d
end interface

interface
 subroutine xcast_mpi_dc4d(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: master
  integer,intent(in) :: spaceComm
  complex(dpc),intent(inout) :: xval(:,:,:,:)
 end subroutine xcast_mpi_dc4d
end interface

interface
 subroutine xcast_mpi_ch0d(xval,master,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: master
  integer,intent(in) :: spaceComm
  character(len=*),intent(inout) :: xval
 end subroutine xcast_mpi_ch0d
end interface

interface
 subroutine xcast_mpi_ch1d(xval,master,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: master
  integer,intent(in) :: spaceComm
  character(len=*),intent(inout) :: xval(:)
 end subroutine xcast_mpi_ch1d
end interface

interface
 subroutine xexch_mpi_intn(vsend,n1,sender,vrecv,recever,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: n1
  integer,intent(in) :: recever
  integer,intent(in) :: sender
  integer,intent(in) :: spaceComm
  integer,intent(inout) :: vrecv(:)
  integer,intent(in) :: vsend(:)
 end subroutine xexch_mpi_intn
end interface

interface
 subroutine xexch_mpi_int2d(vsend,nt,sender,vrecv,recever,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: nt
  integer,intent(in) :: recever
  integer,intent(in) :: sender
  integer,intent(in) :: spaceComm
  integer,intent(inout) :: vrecv(:,:)
  integer,intent(in) :: vsend(:,:)
 end subroutine xexch_mpi_int2d
end interface

interface
 subroutine xexch_mpi_dpn(vsend,n1,sender,vrecv,recever,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: n1
  integer,intent(in) :: recever
  integer,intent(in) :: sender
  integer,intent(in) :: spaceComm
  real(dp),intent(inout) :: vrecv(:)
  real(dp),intent(in) :: vsend(:)
 end subroutine xexch_mpi_dpn
end interface

interface
 subroutine xexch_mpi_dp2d(vsend,nt,sender,vrecv,recever,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: nt
  integer,intent(in) :: recever
  integer,intent(in) :: sender
  integer,intent(in) :: spaceComm
  real(dp),intent(inout) :: vrecv(:,:)
  real(dp),intent(in) :: vsend(:,:)
 end subroutine xexch_mpi_dp2d
end interface

interface
 subroutine xexch_mpi_dp3d(vsend,nt,sender,vrecv,recever,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: nt
  integer,intent(in) :: recever
  integer,intent(in) :: sender
  integer,intent(in) :: spaceComm
  real(dp),intent(inout) :: vrecv(:,:,:)
  real(dp),intent(in) :: vsend(:,:,:)
 end subroutine xexch_mpi_dp3d
end interface

interface
 subroutine xexch_mpi_dp4d_tag(vsend,mtag,sender,vrecv,recever,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: mtag
  integer,intent(in) :: recever
  integer,intent(in) :: sender
  integer,intent(in) :: spaceComm
  real(dp),intent(inout) :: vrecv(:,:,:,:)
  real(dp),intent(in) :: vsend(:,:,:,:)
 end subroutine xexch_mpi_dp4d_tag
end interface

interface
 subroutine xexch_mpi_dp5d_tag(vsend,mtag,sender,vrecv,recever,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: mtag
  integer,intent(in) :: recever
  integer,intent(in) :: sender
  integer,intent(in) :: spaceComm
  real(dp),intent(inout) :: vrecv(:,:,:,:,:)
  real(dp),intent(in) :: vsend(:,:,:,:,:)
 end subroutine xexch_mpi_dp5d_tag
end interface

interface
 subroutine xexch_mpi_spc_1d(vsend,n1,sender,vrecv,recever,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: n1
  integer,intent(in) :: recever
  integer,intent(in) :: sender
  integer,intent(in) :: spaceComm
  complex(spc),intent(inout) :: vrecv(:)
  complex(spc),intent(in) :: vsend(:)
 end subroutine xexch_mpi_spc_1d
end interface

interface
 subroutine xexch_mpi_dpc_1d(vsend,n1,sender,vrecv,recever,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: n1
  integer,intent(in) :: recever
  integer,intent(in) :: sender
  integer,intent(in) :: spaceComm
  complex(dpc),intent(inout) :: vrecv(:)
  complex(dpc),intent(in) :: vsend(:)
 end subroutine xexch_mpi_dpc_1d
end interface

interface
 subroutine xexch_mpi_dpc_2d(vsend,nt,sender,vrecv,recever,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: nt
  integer,intent(in) :: recever
  integer,intent(in) :: sender
  integer,intent(in) :: spaceComm
  complex(dpc),intent(inout) :: vrecv(:,:)
  complex(dpc),intent(in) :: vsend(:,:)
 end subroutine xexch_mpi_dpc_2d
end interface

interface
 subroutine xgather_mpi_int(xval,sendcount,recvbuf,recvcount,root,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer :: recvcount
  integer,intent(in) :: root
  integer :: sendcount
  integer,intent(in) :: spaceComm
  integer,intent(inout) :: recvbuf(:)
  integer,intent(in) :: xval(:)
 end subroutine xgather_mpi_int
end interface

interface
 subroutine xgather_mpi_int2d(xval,sendcount,recvbuf,recvcount,root,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer :: recvcount
  integer,intent(in) :: root
  integer :: sendcount
  integer,intent(in) :: spaceComm
  integer,intent(inout) :: recvbuf(:,:)
  integer,intent(in) :: xval(:,:)
 end subroutine xgather_mpi_int2d
end interface

interface
 subroutine xgather_mpi_dp(xval,sendcount,recvbuf,recvcount,root,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer :: recvcount
  integer,intent(in) :: root
  integer :: sendcount
  integer,intent(in) :: spaceComm
  real(dp),intent(inout) :: recvbuf(:)
  real(dp),intent(in) :: xval(:)
 end subroutine xgather_mpi_dp
end interface

interface
 subroutine xgather_mpi_dp2d(xval,sendcount,recvbuf,recvcount,root,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer :: recvcount
  integer,intent(in) :: root
  integer :: sendcount
  integer,intent(in) :: spaceComm
  real(dp),intent(inout) :: recvbuf(:,:)
  real(dp),intent(in) :: xval(:,:)
 end subroutine xgather_mpi_dp2d
end interface

interface
 subroutine xgather_mpi_dp3d(xval,sendcount,recvbuf,recvcount,root,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer :: recvcount
  integer,intent(in) :: root
  integer :: sendcount
  integer,intent(in) :: spaceComm
  real(dp),intent(inout) :: recvbuf(:,:,:)
  real(dp),intent(in) :: xval(:,:,:)
 end subroutine xgather_mpi_dp3d
end interface

interface
 subroutine xgather_mpi_dp4d(xval,sendcount,recvbuf,recvcount,root,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer :: recvcount
  integer,intent(in) :: root
  integer :: sendcount
  integer,intent(in) :: spaceComm
  real(dp),intent(inout) :: recvbuf(:,:,:,:)
  real(dp),intent(in) :: xval(:,:,:,:)
 end subroutine xgather_mpi_dp4d
end interface

interface
 subroutine xgatherv_mpi_int(xval,nelem,recvbuf,recvcounts,displs,root,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: nelem
  integer,intent(in) :: root
  integer,intent(in) :: spaceComm
  integer,intent(in) :: displs(:)
  integer,intent(inout) :: recvbuf(:)
  integer,intent(in) :: recvcounts(:)
  integer,intent(in) :: xval(:)
 end subroutine xgatherv_mpi_int
end interface

interface
 subroutine xgatherv_mpi_int2d(xval,nelem,recvbuf,recvcounts,displs,root,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: nelem
  integer,intent(in) :: root
  integer,intent(in) :: spaceComm
  integer,intent(in) :: displs(:)
  integer,intent(inout) :: recvbuf(:,:)
  integer,intent(in) :: recvcounts(:)
  integer,intent(in) :: xval(:,:)
 end subroutine xgatherv_mpi_int2d
end interface

interface
 subroutine xgatherv_mpi_dp(xval,nelem,recvbuf,recvcounts,displs,root,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: nelem
  integer,intent(in) :: root
  integer,intent(in) :: spaceComm
  integer,intent(in) :: displs(:)
  integer,intent(in) :: recvcounts(:)
  real(dp),intent(inout) :: recvbuf(:)
  real(dp),intent(in) :: xval(:)
 end subroutine xgatherv_mpi_dp
end interface

interface
 subroutine xgatherv_mpi_dp2d(xval,nelem,recvbuf,recvcounts,displs,root,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: nelem
  integer,intent(in) :: root
  integer,intent(in) :: spaceComm
  integer,intent(in) :: displs(:)
  integer,intent(in) :: recvcounts(:)
  real(dp),intent(inout) :: recvbuf(:,:)
  real(dp),intent(in) :: xval(:,:)
 end subroutine xgatherv_mpi_dp2d
end interface

interface
 subroutine xgatherv_mpi_dp3d(xval,nelem,recvbuf,recvcounts,displs,root,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: nelem
  integer,intent(in) :: root
  integer,intent(in) :: spaceComm
  integer,intent(in) :: displs(:)
  integer,intent(in) :: recvcounts(:)
  real(dp),intent(inout) :: recvbuf(:,:,:)
  real(dp),intent(in) :: xval(:,:,:)
 end subroutine xgatherv_mpi_dp3d
end interface

interface
 subroutine xgatherv_mpi_dp4d(xval,nelem,recvbuf,recvcounts,displs,root,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: nelem
  integer,intent(in) :: root
  integer,intent(in) :: spaceComm
  integer,intent(in) :: displs(:)
  integer,intent(in) :: recvcounts(:)
  real(dp),intent(inout) :: recvbuf(:,:,:,:)
  real(dp),intent(in) :: xval(:,:,:,:)
 end subroutine xgatherv_mpi_dp4d
end interface



interface
 subroutine xmax_mpi_intv(xval,xmax,spaceComm,ier)
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: spaceComm
  integer ,intent(inout) :: xmax
  integer ,intent(inout) :: xval
 end subroutine xmax_mpi_intv
end interface

interface
 subroutine xmax_mpi_dpv(xval,xmax,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xmax
  real(dp),intent(inout) :: xval
 end subroutine xmax_mpi_dpv
end interface

interface
 subroutine xmin_mpi_intv(xval,xmin,spaceComm,ier)
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: spaceComm
  integer ,intent(inout) :: xmin
  integer ,intent(inout) :: xval
 end subroutine xmin_mpi_intv
end interface

interface
 subroutine xmin_mpi_dpv(xval,xmin,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xmin
  real(dp),intent(inout) :: xval
 end subroutine xmin_mpi_dpv
end interface

interface
 subroutine xrecv_mpi_intv(xval,source,tag,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: source
  integer,intent(in) :: spaceComm
  integer,intent(in) :: tag
  integer,intent(inout) :: xval
 end subroutine xrecv_mpi_intv
end interface

interface
 subroutine xrecv_mpi_dp2d(xval,source,tag,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: source
  integer ,intent(in) :: spaceComm
  integer ,intent(in) :: tag
  real(dp),intent(inout) :: xval(:,:)
 end subroutine xrecv_mpi_dp2d
end interface

interface
 subroutine xrecv_mpi_dp3d(xval,source,tag,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: source
  integer ,intent(in) :: spaceComm
  integer ,intent(in) :: tag
  real(dp),intent(inout) :: xval(:,:,:)
 end subroutine xrecv_mpi_dp3d
end interface

interface
 subroutine xsend_mpi_intv(xval,dest,tag,spaceComm,ier)
  implicit none
  integer,intent(in) :: dest
  integer,intent(out) :: ier
  integer,intent(in) :: spaceComm
  integer,intent(in) :: tag
  integer,intent(inout) :: xval
 end subroutine xsend_mpi_intv
end interface

interface
 subroutine xsend_mpi_dp2d(xval,dest,tag,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(in) :: dest
  integer ,intent(out) :: ier
  integer ,intent(in) :: spaceComm
  integer ,intent(in) :: tag
  real(dp),intent(inout) :: xval(:,:)
 end subroutine xsend_mpi_dp2d
end interface

interface
 subroutine xsend_mpi_dp3d(xval,dest,tag,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(in) :: dest
  integer ,intent(out) :: ier
  integer ,intent(in) :: spaceComm
  integer ,intent(in) :: tag
  real(dp),intent(inout) :: xval(:,:,:)
 end subroutine xsend_mpi_dp3d
end interface

interface
 subroutine xsum_master_dp1d(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval(:)
 end subroutine xsum_master_dp1d
end interface

interface
 subroutine xsum_master_dp2d(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval(:,:)
 end subroutine xsum_master_dp2d
end interface

interface
 subroutine xsum_master_dp3d(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval(:,:,:)
 end subroutine xsum_master_dp3d
end interface

interface
 subroutine xsum_master_dp4d(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval(:,:,:,:)
 end subroutine xsum_master_dp4d
end interface

interface
 subroutine xsum_master_dp5d(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval(:,:,:,:,:)
 end subroutine xsum_master_dp5d
end interface

interface
 subroutine xsum_master_dp6d(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval(:,:,:,:,:,:)
 end subroutine xsum_master_dp6d
end interface

interface
 subroutine xsum_master_dp7d(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval(:,:,:,:,:,:,:)
 end subroutine xsum_master_dp7d
end interface

interface
 subroutine xsum_master_int4d(xval,master,spaceComm,ier)
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  integer ,intent(inout) :: xval(:,:,:,:)
 end subroutine xsum_master_int4d
end interface

interface
 subroutine xsum_master_c2cplx(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  complex(spc),intent(inout) :: xval(:,:)
 end subroutine xsum_master_c2cplx
end interface

interface
 subroutine xsum_master_c3cplx(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  complex(spc),intent(inout) :: xval(:,:,:)
 end subroutine xsum_master_c3cplx
end interface

interface
 subroutine xsum_master_c4cplx(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  complex(spc),intent(inout) :: xval(:,:,:,:)
 end subroutine xsum_master_c4cplx
end interface

interface
 subroutine xsum_master_c5cplx(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: master
  integer,intent(in) :: spaceComm
  complex(spc) ,intent(inout) :: xval(:,:,:,:,:)
 end subroutine xsum_master_c5cplx
end interface

interface
 subroutine xsum_master_c1dpc(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: master
  integer,intent(in) :: spaceComm
  complex(dpc) ,intent(inout) :: xval(:)
 end subroutine xsum_master_c1dpc
end interface

interface
 subroutine xsum_master_c2dpc(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  complex(dpc) ,intent(inout) :: xval(:,:)
 end subroutine xsum_master_c2dpc
end interface

interface
 subroutine xsum_master_c3dpc(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  complex(dpc) ,intent(inout) :: xval(:,:,:)
 end subroutine xsum_master_c3dpc
end interface

interface
 subroutine xsum_master_c4dpc(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  complex(dpc) ,intent(inout) :: xval(:,:,:,:)
 end subroutine xsum_master_c4dpc
end interface

interface
 subroutine xsum_master_c5dpc(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: master
  integer,intent(in) :: spaceComm
  complex(dpc) ,intent(inout) :: xval(:,:,:,:,:)
 end subroutine xsum_master_c5dpc
end interface

interface
 subroutine xsum_mpi_int(xval,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: spaceComm
  integer,intent(inout) :: xval(:)
 end subroutine xsum_mpi_int
end interface

interface
 subroutine xsum_mpi_intv(xval,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: spaceComm
  integer,intent(inout) :: xval
 end subroutine xsum_mpi_intv
end interface

interface
 subroutine xsum_mpi_intv2(xval,xsum,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: spaceComm
  integer,intent(inout) :: xsum
  integer,intent(inout) :: xval
 end subroutine xsum_mpi_intv2
end interface

interface
 subroutine xsum_mpi_intn(xval,n1,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: n1
  integer,intent(in) :: spaceComm
  integer,intent(inout) :: xval(:)
 end subroutine xsum_mpi_intn
end interface

interface
 subroutine xsum_mpi_int2t(xval,xsum,n1,spaceComm,ier)
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: n1
  integer ,intent(in) :: spaceComm
  integer ,intent(inout) :: xsum(:)
  integer ,intent(inout) :: xval(:)
 end subroutine xsum_mpi_int2t
end interface

interface
 subroutine xsum_mpi_int2d(xval,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: spaceComm
  integer,intent(inout) :: xval(:,:)
 end subroutine xsum_mpi_int2d
end interface

interface
 subroutine xsum_mpi_int3d(xval,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: spaceComm
  integer,intent(inout) :: xval(:,:,:)
 end subroutine xsum_mpi_int3d
end interface

interface
 subroutine xsum_mpi_int4d(xval,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: spaceComm
  integer,intent(inout) :: xval(:,:,:,:)
 end subroutine xsum_mpi_int4d
end interface

interface
 subroutine xsum_mpi_dp(xval,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval(:)
 end subroutine xsum_mpi_dp
end interface

interface
 subroutine xsum_mpi_dpvt(xval,xsum,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: spaceComm
  real(dp),intent(out) :: xsum
  real(dp),intent(in) :: xval
 end subroutine xsum_mpi_dpvt
end interface

interface
 subroutine xsum_mpi_dpv(xval,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval
 end subroutine xsum_mpi_dpv
end interface

interface
 subroutine xsum_mpi_dpn(xval,n1,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: n1
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval(:)
 end subroutine xsum_mpi_dpn
end interface

interface
 subroutine xsum_mpi_dp2d(xval,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval(:,:)
 end subroutine xsum_mpi_dp2d
end interface

interface
 subroutine xsum_mpi_dp3d(xval,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval(:,:,:)
 end subroutine xsum_mpi_dp3d
end interface

interface
 subroutine xsum_mpi_dp4d(xval,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval(:,:,:,:)
 end subroutine xsum_mpi_dp4d
end interface

interface
 subroutine xsum_mpi_dp5d(xval,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval(:,:,:,:,:)
 end subroutine xsum_mpi_dp5d
end interface

interface
 subroutine xsum_mpi_dp6d(xval,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval(:,:,:,:,:,:)
 end subroutine xsum_mpi_dp6d
end interface

interface
 subroutine xsum_mpi_dp2t(xval,xsum,n1,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: n1
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xsum(:)
  real(dp),intent(inout) :: xval(:)
 end subroutine xsum_mpi_dp2t
end interface

interface
 subroutine xsum_mpi_dp3d2t(xval,xsum,n1,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: n1
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xsum(:,:,:)
  real(dp),intent(inout) :: xval(:,:,:)
 end subroutine xsum_mpi_dp3d2t
end interface

interface
 subroutine xsum_mpi_dp4d2t(xval,xsum,n1,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: n1
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xsum(:,:,:,:)
  real(dp),intent(inout) :: xval(:,:,:,:)
 end subroutine xsum_mpi_dp4d2t
end interface

interface
 subroutine xsum_mpi_c0dc(xval,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: spaceComm
  complex(dpc),intent(inout) :: xval
 end subroutine xsum_mpi_c0dc
end interface

interface
 subroutine xsum_mpi_c1dc(xval,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: spaceComm
  complex(dpc),intent(inout) :: xval(:)
 end subroutine xsum_mpi_c1dc
end interface

interface
 subroutine xsum_mpi_c2dc(xval,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: spaceComm
  complex(dpc),intent(inout) :: xval(:,:)
 end subroutine xsum_mpi_c2dc
end interface

interface
 subroutine xsum_mpi_c3dc(xval,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: spaceComm
  complex(dpc),intent(inout) :: xval(:,:,:)
 end subroutine xsum_mpi_c3dc
end interface

interface
 subroutine xsum_mpi_c4dc(xval,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: spaceComm
  complex(dpc),intent(inout) :: xval(:,:,:,:)
 end subroutine xsum_mpi_c4dc
end interface

interface
 subroutine xsum_mpi_c5dc(xval,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: spaceComm
  complex(dpc),intent(inout) :: xval(:,:,:,:,:)
 end subroutine xsum_mpi_c5dc
end interface

interface
 subroutine xsum_mpi_c6dc(xval,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: spaceComm
  complex(dpc),intent(inout) :: xval(:,:,:,:,:,:)
 end subroutine xsum_mpi_c6dc
end interface

interface
 subroutine xsum_mpi_c1cplx(xval,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: spaceComm
  complex,intent(inout) :: xval(:)
 end subroutine xsum_mpi_c1cplx
end interface

interface
 subroutine xsum_mpi_c2cplx(xval,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: spaceComm
  complex,intent(inout) :: xval(:,:)
 end subroutine xsum_mpi_c2cplx
end interface

interface
 subroutine xsum_mpi_c3cplx(xval,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: spaceComm
  complex,intent(inout) :: xval(:,:,:)
 end subroutine xsum_mpi_c3cplx
end interface

interface
 subroutine xsum_mpi_c4cplx(xval,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: spaceComm
  complex,intent(inout) :: xval(:,:,:,:)
 end subroutine xsum_mpi_c4cplx
end interface

interface
 subroutine xsum_mpi_c5cplx(xval,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: spaceComm
  complex,intent(inout) :: xval(:,:,:,:,:)
 end subroutine xsum_mpi_c5cplx
end interface

interface
 subroutine xsum_mpi_c6cplx(xval,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: spaceComm
  complex,intent(inout) :: xval(:,:,:,:,:,:)
 end subroutine xsum_mpi_c6cplx
end interface

interface
 subroutine xsum_mpi_log1d(xval,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: spaceComm
  logical,intent(inout) :: xval(:)
 end subroutine xsum_mpi_log1d
end interface

interface
 subroutine xsum_mpi_log2d(xval,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: spaceComm
  logical,intent(inout) :: xval(:,:)
 end subroutine xsum_mpi_log2d
end interface

end module interfaces_12_hide_mpi
!!***
