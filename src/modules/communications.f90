module communications
  use module_base
  use module_types
  implicit none

  private


  public :: collective_comms_null
  public :: allocate_MPI_communication_arrays, allocate_local_communications_arrays

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


    subroutine allocate_local_communications_arrays(comms)
      implicit none
      type(collective_comms),intent(inout) :: comms
      integer :: istat
      character(len=*),parameter :: subname='allocate_MPI_communication_arrays'

      allocate(comms%irecvbuf_c(comms%ndimpsi_c), stat=istat)
      call memocc(istat, comms%irecvbuf_c, 'comms%irecvbuf_c', subname)
      allocate(comms%indexrecvorbital_c(comms%ndimind_c), stat=istat)
      call memocc(istat, comms%indexrecvorbital_c, 'comms%indexrecvorbital_c', subname)
      allocate(comms%iextract_c(comms%ndimind_c), stat=istat)
      call memocc(istat, comms%iextract_c, 'comms%iextract_c', subname)
      allocate(comms%iexpand_c(comms%ndimind_c), stat=istat)
      call memocc(istat, comms%iexpand_c, 'comms%iexpand_c', subname)
      allocate(comms%isendbuf_c(comms%ndimpsi_c), stat=istat)
      call memocc(istat, comms%isendbuf_c, 'comms%isendbuf_c', subname)

      allocate(comms%irecvbuf_f(comms%ndimpsi_f), stat=istat)
      call memocc(istat, comms%irecvbuf_f, 'comms%irecvbuf_f', subname)
      allocate(comms%indexrecvorbital_f(comms%ndimind_f), stat=istat)
      call memocc(istat, comms%indexrecvorbital_f, 'comms%indexrecvorbital_f', subname)
      allocate(comms%iextract_f(comms%ndimind_f), stat=istat)
      call memocc(istat, comms%iextract_f, 'comms%iextract_f', subname)
      allocate(comms%iexpand_f(comms%ndimind_f), stat=istat)
      call memocc(istat, comms%iexpand_f, 'comms%iexpand_f', subname)
      allocate(comms%isendbuf_f(comms%ndimpsi_f), stat=istat)
      call memocc(istat, comms%isendbuf_f, 'comms%isendbuf_f', subname)

      allocate(comms%isptsp_c(max(comms%nptsp_c,1)), stat=istat)
      call memocc(istat, comms%isptsp_c, 'comms%isptsp_c', subname)
      allocate(comms%isptsp_f(max(comms%nptsp_f,1)), stat=istat)
      call memocc(istat, comms%isptsp_f, 'comms%isptsp_f', subname)

      allocate(comms%norb_per_gridpoint_c(comms%nptsp_c), stat=istat)
      call memocc(istat, comms%norb_per_gridpoint_c, 'comms%norb_per_gridpoint_c', subname)
      allocate(comms%norb_per_gridpoint_f(comms%nptsp_f), stat=istat)
      call memocc(istat, comms%norb_per_gridpoint_f, 'comms%norb_per_gridpoint_f', subname)

    end subroutine allocate_local_communications_arrays


end module communications
