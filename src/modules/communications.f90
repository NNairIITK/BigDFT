module communications
  use module_base
  use module_types
  implicit none

  private


  public :: collective_comms_null
  public :: allocate_MPI_communication_arrays

  contains

    !> Creators and destructors

    pure function collective_comms_null() result(comms)
      implicit none
      type(collective_comms) :: comms
      call nullify_collective_comms(comms)
    end function collective_comms_null

    pure subroutine nullify_collective_comms(comms)
      implicit none
      type(collective_comms),intent(inout) :: comms
      nullify(comms%nsendcounts_c)
      nullify(comms%nsenddspls_c)
      nullify(comms%nrecvcounts_c)
      nullify(comms%nrecvdspls_c)
      nullify(comms%isendbuf_c)
      nullify(comms%iextract_c)
      nullify(comms%iexpand_c)
      nullify(comms%irecvbuf_c)
      nullify(comms%norb_per_gridpoint_c)
      nullify(comms%indexrecvorbital_c)
      nullify(comms%isptsp_c)
      nullify(comms%psit_c)
      nullify(comms%nsendcounts_f)
      nullify(comms%nsenddspls_f)
      nullify(comms%nrecvcounts_f)
      nullify(comms%nrecvdspls_f)
      nullify(comms%isendbuf_f)
      nullify(comms%iextract_f)
      nullify(comms%iexpand_f)
      nullify(comms%irecvbuf_f)
      nullify(comms%norb_per_gridpoint_f)
      nullify(comms%indexrecvorbital_f)
      nullify(comms%isptsp_f)
      nullify(comms%psit_f)
      nullify(comms%nsendcounts_repartitionrho)
      nullify(comms%nrecvcounts_repartitionrho)
      nullify(comms%nsenddspls_repartitionrho)
      nullify(comms%nrecvdspls_repartitionrho)
      nullify(comms%commarr_repartitionrho)
    end subroutine nullify_collective_comms

    subroutine allocate_MPI_communication_arrays(nproc, comms)
      implicit none
      integer,intent(in) :: nproc
      type(collective_comms),intent(inout) :: comms
      integer :: istat
      character(len=*),parameter :: subname='allocate_MPI_communication_arrays'
      allocate(comms%nsendcounts_c(0:nproc-1), stat=istat)
      call memocc(istat, comms%nsendcounts_c, 'comms%nsendcounts_c', subname)
      allocate(comms%nsenddspls_c(0:nproc-1), stat=istat)
      call memocc(istat, comms%nsenddspls_c, 'comms%nsenddspls_c', subname)
      allocate(comms%nrecvcounts_c(0:nproc-1), stat=istat)
      call memocc(istat, comms%nrecvcounts_c, 'comms%nrecvcounts_c', subname)
      allocate(comms%nrecvdspls_c(0:nproc-1), stat=istat)
      call memocc(istat, comms%nrecvdspls_c, 'comms%nrecvdspls_c', subname)
      allocate(comms%nsendcounts_f(0:nproc-1), stat=istat)
      call memocc(istat, comms%nsendcounts_f, 'comms%nsendcounts_f', subname)
      allocate(comms%nsenddspls_f(0:nproc-1), stat=istat)
      call memocc(istat, comms%nsenddspls_f, 'comms%nsenddspls_f', subname)
      allocate(comms%nrecvcounts_f(0:nproc-1), stat=istat)
      call memocc(istat, comms%nrecvcounts_f, 'comms%nrecvcounts_f', subname)
      allocate(comms%nrecvdspls_f(0:nproc-1), stat=istat)
      call memocc(istat, comms%nrecvdspls_f, 'comms%nrecvdspls_f', subname)
    end subroutine allocate_MPI_communication_arrays


end module communications
