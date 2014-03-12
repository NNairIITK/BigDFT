module communications
  use module_base
  use module_types
  implicit none

  private


  public :: collective_comms_null
  public :: allocate_MPI_communication_arrays
  public :: allocate_local_communications_arrays
  public :: allocate_MPI_communications_arrays_repartition
  public :: init_collective_comms

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
      character(len=*),parameter :: subname='allocate_local_communications_arrays'

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


    subroutine allocate_MPI_communications_arrays_repartition(nproc, comms)
      implicit none
      integer,intent(in) :: nproc
      type(collective_comms),intent(inout) :: comms
      integer :: istat
      character(len=*),parameter :: subname='allocate_MPI_communications_arrays_repartition'
      allocate(comms%nsendcounts_repartitionrho(0:nproc-1), stat=istat)
      call memocc(istat, comms%nsendcounts_repartitionrho, 'comms%nsendcounts_repartitionrho', subname)
      allocate(comms%nrecvcounts_repartitionrho(0:nproc-1), stat=istat)
      call memocc(istat, comms%nrecvcounts_repartitionrho, 'comms%nrecvcounts_repartitionrho', subname)
      allocate(comms%nsenddspls_repartitionrho(0:nproc-1), stat=istat)
      call memocc(istat, comms%nsenddspls_repartitionrho, 'comms%nsenddspls_repartitionrho', subname)
      allocate(comms%nrecvdspls_repartitionrho(0:nproc-1), stat=istat)
      call memocc(istat, comms%nrecvdspls_repartitionrho, 'comms%nrecvdspls_repartitionrho', subname)
    end subroutine allocate_MPI_communications_arrays_repartition


    subroutine init_collective_comms(iproc, nproc, npsidim_orbs, orbs, lzd, collcom)
      use module_base
      use module_types
      use module_interfaces, except_this_one => init_collective_comms
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, npsidim_orbs
      type(orbitals_data),intent(in) :: orbs
      type(local_zone_descriptors),intent(in) :: lzd
      type(collective_comms),intent(inout) :: collcom
      
      ! Local variables
      integer :: ii, iorb, iiorb, ilr, istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, ierr
      integer :: ipt, nvalp_c, nvalp_f
      real(kind=8),dimension(:,:,:),allocatable :: weight_c, weight_f
      real(kind=8) :: weight_c_tot, weight_f_tot, weightp_c, weightp_f
      integer,dimension(:,:),allocatable :: istartend_c, istartend_f
      integer,dimension(:,:,:),allocatable :: index_in_global_c, index_in_global_f
      integer,dimension(:),allocatable :: npts_par_c, npts_par_f
      
      call timing(iproc,'init_collcomm ','ON')
    
      call f_routine('init_collective_comms')
    
      weight_c=f_malloc((/0.to.lzd%glr%d%n1,0.to.lzd%glr%d%n2,0.to.lzd%glr%d%n3/))
      weight_f=f_malloc((/0.to.lzd%glr%d%n1,0.to.lzd%glr%d%n2,0.to.lzd%glr%d%n3/))
      index_in_global_c=f_malloc((/0.to.lzd%glr%d%n1,0.to.lzd%glr%d%n2,0.to.lzd%glr%d%n3/))
      index_in_global_f=f_malloc((/0.to.lzd%glr%d%n1,0.to.lzd%glr%d%n2,0.to.lzd%glr%d%n3/))
      
    
      call get_weights(iproc, nproc, orbs, lzd, weight_c, weight_f, weight_c_tot, weight_f_tot)
    
      ! Assign the grid points to the processes such that the work is equally dsitributed
      istartend_c=f_malloc((/1.to.2,0.to.nproc-1/))
      istartend_f=f_malloc((/1.to.2,0.to.nproc-1/))
      call assign_weight_to_process(iproc, nproc, lzd, weight_c, weight_f, weight_c_tot, weight_f_tot, &
           istartend_c, istartend_f, istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, &
           weightp_c, weightp_f, collcom%nptsp_c, collcom%nptsp_f, nvalp_c, nvalp_f)
    
    
      ! Determine the index of a grid point i1,i2,i3 in the compressed array
      call get_index_in_global2(lzd%glr, index_in_global_c, index_in_global_f)
    
    
      ! Determine values for mpi_alltoallv
      call allocate_MPI_communication_arrays(nproc, collcom)
      call determine_communication_arrays(iproc, nproc, npsidim_orbs, orbs, lzd, istartend_c, istartend_f, &
           index_in_global_c, index_in_global_f, nvalp_c, nvalp_f, &
           collcom%nsendcounts_c, collcom%nsenddspls_c, collcom%nrecvcounts_c, collcom%nrecvdspls_c, &
           collcom%nsendcounts_f, collcom%nsenddspls_f, collcom%nrecvcounts_f, collcom%nrecvdspls_f)
    
    
    
      !Now set some integers in the collcomm structure
      collcom%ndimind_c = sum(collcom%nrecvcounts_c)
      collcom%ndimind_f = sum(collcom%nrecvcounts_f)
    
      ! Now rearrange the data on the process to communicate them
      collcom%ndimpsi_c=0
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          ilr=orbs%inwhichlocreg(iiorb)
          collcom%ndimpsi_c=collcom%ndimpsi_c+lzd%llr(ilr)%wfd%nvctr_c
      end do
      collcom%ndimpsi_f=0
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          ilr=orbs%inwhichlocreg(iiorb)
          collcom%ndimpsi_f=collcom%ndimpsi_f+lzd%llr(ilr)%wfd%nvctr_f
      end do
    
      call allocate_local_communications_arrays(collcom)
    
      call determine_num_orbs_per_gridpoint_new(iproc, nproc, lzd, istartend_c, istartend_f, &
           istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, &
           weightp_c, weightp_f, collcom%nptsp_c, collcom%nptsp_f, weight_c, weight_f, &
           collcom%norb_per_gridpoint_c, collcom%norb_per_gridpoint_f)
    
      call f_free(weight_c)
      call f_free(weight_f)
    
    
      call get_switch_indices(iproc, nproc, orbs, lzd, collcom%ndimpsi_c, collcom%ndimpsi_f, istartend_c, istartend_f, &
           collcom%nsendcounts_c, collcom%nsenddspls_c, collcom%ndimind_c, collcom%nrecvcounts_c, collcom%nrecvdspls_c, &
           collcom%nsendcounts_f, collcom%nsenddspls_f, collcom%ndimind_f, collcom%nrecvcounts_f, collcom%nrecvdspls_f, &
           index_in_global_c, index_in_global_f, &
           weightp_c, weightp_f, collcom%isendbuf_c, collcom%irecvbuf_c, collcom%isendbuf_f, collcom%irecvbuf_f, &
           collcom%indexrecvorbital_c, collcom%iextract_c, collcom%iexpand_c, &
           collcom%indexrecvorbital_f, collcom%iextract_f, collcom%iexpand_f)
    
    
      ! These variables are used in various subroutines to speed up the code
      collcom%isptsp_c(1) = 0
      do ipt=2,collcom%nptsp_c
            collcom%isptsp_c(ipt) = collcom%isptsp_c(ipt-1) + collcom%norb_per_gridpoint_c(ipt-1)
      end do
      if (maxval(collcom%isptsp_c)>collcom%ndimind_c) stop 'maxval(collcom%isptsp_c)>collcom%ndimind_c'
    
      collcom%isptsp_f(1) = 0
      do ipt=2,collcom%nptsp_f
            collcom%isptsp_f(ipt) = collcom%isptsp_f(ipt-1) + collcom%norb_per_gridpoint_f(ipt-1)
      end do
      if (maxval(collcom%isptsp_f)>collcom%ndimind_f) stop 'maxval(collcom%isptsp_f)>collcom%ndimind_f'
    
    
      ! Not used any more, so deallocate...
      call f_free(istartend_c)
      call f_free(istartend_f)
    
      call f_free(index_in_global_c)
      call f_free(index_in_global_f)
      
    
      call f_release_routine()
      
      call timing(iproc,'init_collcomm ','OF')
      
    end subroutine init_collective_comms

end module communications
