module communications
  use module_base
  use communications_type
  implicit none

  private

  !!type,public:: collective_comms
  !!  integer :: nptsp_c, ndimpsi_c, ndimind_c, ndimind_f, nptsp_f, ndimpsi_f
  !!  integer,dimension(:),pointer :: nsendcounts_c, nsenddspls_c, nrecvcounts_c, nrecvdspls_c
  !!  integer,dimension(:),pointer :: isendbuf_c, iextract_c, iexpand_c, irecvbuf_c
  !!  integer,dimension(:),pointer :: norb_per_gridpoint_c, indexrecvorbital_c
  !!  integer,dimension(:),pointer :: nsendcounts_f, nsenddspls_f, nrecvcounts_f, nrecvdspls_f
  !!  integer,dimension(:),pointer :: isendbuf_f, iextract_f, iexpand_f, irecvbuf_f
  !!  integer,dimension(:),pointer :: norb_per_gridpoint_f, indexrecvorbital_f
  !!  integer,dimension(:),pointer :: isptsp_c, isptsp_f !<starting index of a given gridpoint (basically summation of norb_per_gridpoint_*)
  !!  real(kind=8),dimension(:),pointer :: psit_c, psit_f
  !!  integer,dimension(:),pointer :: nsendcounts_repartitionrho, nrecvcounts_repartitionrho
  !!  integer,dimension(:),pointer :: nsenddspls_repartitionrho, nrecvdspls_repartitionrho
  !!  integer :: ncomms_repartitionrho, window
  !!  integer,dimension(:,:),pointer :: commarr_repartitionrho
  !!end type collective_comms



  public :: collective_comms_null
  public :: init_collective_comms
  public :: transpose_localized
  public :: untranspose_localized
  public :: init_collective_comms_sumrho
  public :: transpose_switch_psir
  public :: transpose_communicate_psir
  public :: transpose_unswitch_psirt
  public :: deallocate_collective_comms

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

    subroutine allocate_MPI_communication_arrays(nproc, comms, only_coarse)
      implicit none
      integer,intent(in) :: nproc
      type(collective_comms),intent(inout) :: comms
      logical,intent(in),optional :: only_coarse
      logical :: allocate_fine
      if (present(only_coarse)) then
          allocate_fine=.not.only_coarse
      else
          allocate_fine=.true.
      end if
      comms%nsendcounts_c=f_malloc_ptr(0.to.nproc-1,id='comms%nsendcounts_c')
      comms%nsenddspls_c=f_malloc_ptr(0.to.nproc-1,id='comms%nsenddspls_c')
      comms%nrecvcounts_c=f_malloc_ptr(0.to.nproc-1,id='comms%nrecvcounts_c')
      comms%nrecvdspls_c=f_malloc_ptr(0.to.nproc-1,id='comms%nrecvdspls_c')
      if (allocate_fine) then
          comms%nsendcounts_f=f_malloc_ptr(0.to.nproc-1,id='comms%nsendcounts_f')
          comms%nsenddspls_f=f_malloc_ptr(0.to.nproc-1,id='comms%nsenddspls_f')
          comms%nrecvcounts_f=f_malloc_ptr(0.to.nproc-1,id='comms%nrecvcounts_f')
          comms%nrecvdspls_f=f_malloc_ptr(0.to.nproc-1,id='comms%nrecvdspls_f')
      end if
    end subroutine allocate_MPI_communication_arrays


    subroutine allocate_local_communications_arrays(comms, only_coarse)
      implicit none
      type(collective_comms),intent(inout) :: comms
      logical,intent(in),optional :: only_coarse
      logical :: allocate_fine
      if (present(only_coarse)) then
          allocate_fine=.not.only_coarse
      else
          allocate_fine=.true.
      end if
      comms%irecvbuf_c=f_malloc_ptr(comms%ndimpsi_c,id='comms%irecvbuf_c')
      comms%indexrecvorbital_c=f_malloc_ptr(comms%ndimind_c,id='comms%indexrecvorbital_c')
      comms%iextract_c=f_malloc_ptr(comms%ndimind_c,id='comms%iextract_c')
      comms%iexpand_c=f_malloc_ptr(comms%ndimind_c,id='comms%iexpand_c')
      comms%isendbuf_c=f_malloc_ptr(comms%ndimpsi_c,id='comms%isendbuf_c')
      comms%isptsp_c=f_malloc_ptr(max(comms%nptsp_c,1),id='comms%isptsp_c')
      comms%norb_per_gridpoint_c=f_malloc_ptr(comms%nptsp_c,id='comms%norb_per_gridpoint_c')
      if (allocate_fine) then
          comms%irecvbuf_f=f_malloc_ptr(comms%ndimpsi_f,id='comms%irecvbuf_f')
          comms%indexrecvorbital_f=f_malloc_ptr(comms%ndimind_f,id='comms%indexrecvorbital_f')
          comms%iextract_f=f_malloc_ptr(comms%ndimind_f,id='comms%iextract_f')
          comms%iexpand_f=f_malloc_ptr(comms%ndimind_f,id='comms%iexpand_f')
          comms%isendbuf_f=f_malloc_ptr(comms%ndimpsi_f,id='comms%isendbuf_f')
          comms%isptsp_f=f_malloc_ptr(max(comms%nptsp_f,1),id='comms%isptsp_f')
          comms%norb_per_gridpoint_f=f_malloc_ptr(comms%nptsp_f,id='comms%norb_per_gridpoint_f')
      end if
    end subroutine allocate_local_communications_arrays


    subroutine allocate_MPI_communications_arrays_repartition(nproc, comms)
      implicit none
      integer,intent(in) :: nproc
      type(collective_comms),intent(inout) :: comms
      comms%nsendcounts_repartitionrho=f_malloc_ptr(0.to.nproc-1,id='comms%nsendcounts_repartitionrho')
      comms%nrecvcounts_repartitionrho=f_malloc_ptr(0.to.nproc-1,id='comms%nrecvcounts_repartitionrho')
      comms%nsenddspls_repartitionrho=f_malloc_ptr(0.to.nproc-1,id='comms%nsenddspls_repartitionrho')
      comms%nrecvdspls_repartitionrho=f_malloc_ptr(0.to.nproc-1,id='comms%nrecvdspls_repartitionrho')
    end subroutine allocate_MPI_communications_arrays_repartition


    subroutine allocate_MPI_communications_arrays_repartitionp2p(ncommunications, commarr_repartitionrho)
      implicit none
      integer,intent(in) :: ncommunications
      integer,dimension(:,:),pointer,intent(inout) :: commarr_repartitionrho
      commarr_repartitionrho=f_malloc_ptr((/4,ncommunications/),id='commarr_repartitionrho')
    end subroutine allocate_MPI_communications_arrays_repartitionp2p


    subroutine deallocate_collective_comms(comms)
      implicit none
      type(collective_comms),intent(inout) :: comms
      call deallocate_MPI_communication_arrays(comms)
      call deallocate_local_communications_arrays(comms)
      call deallocate_MPI_communications_arrays_repartition(comms)
      if (associated(comms%psit_c)) call f_free_ptr(comms%psit_c)
      if (associated(comms%psit_f)) call f_free_ptr(comms%psit_f)
      call deallocate_MPI_communications_arrays_repartitionp2p(comms%commarr_repartitionrho)
    end subroutine deallocate_collective_comms


    subroutine deallocate_MPI_communication_arrays(comms)
      implicit none
      type(collective_comms),intent(inout) :: comms
      if (associated(comms%nsendcounts_c)) call f_free_ptr(comms%nsendcounts_c)
      if (associated(comms%nsenddspls_c)) call f_free_ptr(comms%nsenddspls_c)
      if (associated(comms%nrecvcounts_c)) call f_free_ptr(comms%nrecvcounts_c)
      if (associated(comms%nrecvdspls_c)) call f_free_ptr(comms%nrecvdspls_c)
      if (associated(comms%nsendcounts_f)) call f_free_ptr(comms%nsendcounts_f)
      if (associated(comms%nsenddspls_f)) call f_free_ptr(comms%nsenddspls_f)
      if (associated(comms%nrecvcounts_f)) call f_free_ptr(comms%nrecvcounts_f)
      if (associated(comms%nrecvdspls_f)) call f_free_ptr(comms%nrecvdspls_f)
    end subroutine deallocate_MPI_communication_arrays

    subroutine deallocate_local_communications_arrays(comms)
      implicit none
      type(collective_comms),intent(inout) :: comms
      if (associated(comms%irecvbuf_c)) call f_free_ptr(comms%irecvbuf_c)
      if (associated(comms%indexrecvorbital_c)) call f_free_ptr(comms%indexrecvorbital_c)
      if (associated(comms%iextract_c)) call f_free_ptr(comms%iextract_c)
      if (associated(comms%iexpand_c)) call f_free_ptr(comms%iexpand_c)
      if (associated(comms%isendbuf_c)) call f_free_ptr(comms%isendbuf_c)
      if (associated(comms%irecvbuf_f)) call f_free_ptr(comms%irecvbuf_f)
      if (associated(comms%indexrecvorbital_f)) call f_free_ptr(comms%indexrecvorbital_f)
      if (associated(comms%iextract_f)) call f_free_ptr(comms%iextract_f)
      if (associated(comms%iexpand_f)) call f_free_ptr(comms%iexpand_f)
      if (associated(comms%isendbuf_f)) call f_free_ptr(comms%isendbuf_f)
      if (associated(comms%isptsp_c)) call f_free_ptr(comms%isptsp_c)
      if (associated(comms%isptsp_f)) call f_free_ptr(comms%isptsp_f)
      if (associated(comms%norb_per_gridpoint_c)) call f_free_ptr(comms%norb_per_gridpoint_c)
      if (associated(comms%norb_per_gridpoint_f)) call f_free_ptr(comms%norb_per_gridpoint_f)
    end subroutine deallocate_local_communications_arrays

    subroutine deallocate_MPI_communications_arrays_repartition(comms)
      implicit none
      type(collective_comms),intent(inout) :: comms
      if (associated(comms%nsendcounts_repartitionrho)) call f_free_ptr(comms%nsendcounts_repartitionrho)
      if (associated(comms%nrecvcounts_repartitionrho)) call f_free_ptr(comms%nrecvcounts_repartitionrho)
      if (associated(comms%nsenddspls_repartitionrho)) call f_free_ptr(comms%nsenddspls_repartitionrho)
      if (associated(comms%nrecvdspls_repartitionrho)) call f_free_ptr(comms%nrecvdspls_repartitionrho)
    end subroutine deallocate_MPI_communications_arrays_repartition


    subroutine deallocate_MPI_communications_arrays_repartitionp2p(commarr_repartitionrho)
      implicit none
      integer,dimension(:,:),pointer,intent(inout) :: commarr_repartitionrho
      call f_free_ptr(commarr_repartitionrho)
    end subroutine deallocate_MPI_communications_arrays_repartitionp2p



    subroutine init_collective_comms(iproc, nproc, npsidim_orbs, orbs, lzd, collcom)
      use module_base
      use module_types
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


    subroutine get_weights(iproc, nproc, orbs, lzd, weight_c, weight_f, weight_c_tot, weight_f_tot)
      use module_base
      use module_types
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(orbitals_data),intent(in) :: orbs
      type(local_zone_descriptors),intent(in) :: lzd
      real(kind=8),dimension(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3),intent(out) :: weight_c, weight_f
      real(kind=8),intent(out) :: weight_c_tot, weight_f_tot
      
      ! Local variables
      integer :: iorb, iiorb, i0, i1, i2, i3, ii, jj, iseg, ierr, ilr, istart, iend, i, j0, j1, ii1, ii2, ii3
    
    
      ii=(lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(lzd%glr%d%n3+1)
      call to_zero(ii, weight_c(0,0,0))
      call to_zero(ii, weight_f(0,0,0))
      weight_c_tot=0.d0
      weight_f_tot=0.d0
    
      !$omp parallel default(private) &
      !$omp shared(orbs,lzd,weight_c,weight_c_tot,weight_f,weight_f_tot,ilr,iiorb)
    
    
      ! Calculate the weights for the coarse part.
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          ilr=orbs%inwhichlocreg(iiorb)
          if (lzd%llr(ilr)%wfd%nseg_c>0) then
              !$omp do reduction(+:weight_c_tot) 
              do iseg=1,lzd%llr(ilr)%wfd%nseg_c
                  jj=lzd%llr(ilr)%wfd%keyvloc(iseg)
                  j0=lzd%llr(ilr)%wfd%keygloc(1,iseg)
                  j1=lzd%llr(ilr)%wfd%keygloc(2,iseg)
                  ii=j0-1
                  i3=ii/((lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1))
                  ii=ii-i3*(lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)
                  i2=ii/(lzd%llr(ilr)%d%n1+1)
                  i0=ii-i2*(lzd%llr(ilr)%d%n1+1)
                  i1=i0+j1-j0
                  !write(*,'(a,8i8)') 'jj, ii, j0, j1, i0, i1, i2, i3',jj,ii,j0,j1,i0,i1,i2,i3
                  do i=i0,i1
                      ii1=i+lzd%llr(ilr)%ns1
                      ii2=i2+lzd%llr(ilr)%ns2
                      ii3=i3+lzd%llr(ilr)%ns3
                      weight_c(ii1,ii2,ii3)=weight_c(ii1,ii2,ii3)+1.d0
                      weight_c_tot=weight_c_tot+1.d0
                  end do
              end do
              !$omp end do
          end if
    
          ! Calculate the weights for the fine part.
          istart=lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)
          iend=istart+lzd%llr(ilr)%wfd%nseg_f-1
          if (istart<=iend) then
              !$omp do reduction(+:weight_f_tot)
              do iseg=istart,iend
                  jj=lzd%llr(ilr)%wfd%keyvloc(iseg)
                  j0=lzd%llr(ilr)%wfd%keygloc(1,iseg)
                  j1=lzd%llr(ilr)%wfd%keygloc(2,iseg)
                  ii=j0-1
                  i3=ii/((lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1))
                  ii=ii-i3*(lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)
                  i2=ii/(lzd%llr(ilr)%d%n1+1)
                  i0=ii-i2*(lzd%llr(ilr)%d%n1+1)
                  i1=i0+j1-j0
                  do i=i0,i1
                      ii1=i+lzd%llr(ilr)%ns1
                      ii2=i2+lzd%llr(ilr)%ns2
                      ii3=i3+lzd%llr(ilr)%ns3
                      weight_f(ii1,ii2,ii3)=weight_f(ii1,ii2,ii3)+1.d0
                      weight_f_tot=weight_f_tot+1.d0
                  end do
              end do
              !$omp end do
          end if
      end do
      !$omp end parallel
    
      ! Sum up among all processes.
      if(nproc>1) then
          call mpiallred(weight_c_tot, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
          call mpiallred(weight_f_tot, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
          ii=(lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(lzd%glr%d%n3+1)
          call mpiallred(weight_c(0,0,0), ii,  mpi_sum, bigdft_mpi%mpi_comm, ierr)
          call mpiallred(weight_f(0,0,0), ii,  mpi_sum, bigdft_mpi%mpi_comm, ierr)
      end if
    
      weight_c_tot=0.d0
      weight_f_tot=0.d0
      !$omp parallel default(none) &
      !$omp shared(lzd, weight_c, weight_f, weight_c_tot, weight_f_tot) private(i3, i2, i1)
      !$omp do reduction(+: weight_c_tot, weight_f_tot)
      do i3=0,lzd%glr%d%n3
          do i2=0,lzd%glr%d%n2
              do i1=0,lzd%glr%d%n1
                  weight_c(i1,i2,i3)=weight_c(i1,i2,i3)**2
                  weight_f(i1,i2,i3)=weight_f(i1,i2,i3)**2
                  weight_c_tot=weight_c_tot+weight_c(i1,i2,i3)
                  weight_f_tot=weight_f_tot+weight_f(i1,i2,i3)
              end do
          end do
      end do
      !$omp end do
      !$omp end parallel
    
    
    end subroutine get_weights



    subroutine assign_weight_to_process(iproc, nproc, lzd, weight_c, weight_f, weight_tot_c, weight_tot_f, &
               istartend_c, istartend_f, istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, &
               weightp_c, weightp_f, nptsp_c, nptsp_f, nvalp_c, nvalp_f)
      use module_base
      use module_types
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(local_zone_descriptors),intent(in) :: lzd
      real(kind=8),dimension(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3),intent(in) :: weight_c, weight_f
      real(kind=8),intent(in) :: weight_tot_c, weight_tot_f
      integer,dimension(2,0:nproc-1),intent(out) :: istartend_c, istartend_f
      integer,intent(out) :: istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f
      real(kind=8),intent(out) :: weightp_c, weightp_f
      integer,intent(out) :: nptsp_c, nptsp_f
      integer,intent(out) :: nvalp_c, nvalp_f
      
      ! Local variables
      integer :: jproc, i1, i2, i3, ii, istart, iend, jj, j0, j1, ii_c, ii_f
      !!$$integer :: ii2, iiseg, jprocdone
      integer :: i, iseg, i0, iitot, ierr, istat, iall
      real(kind=8) :: tt, tt2, weight_c_ideal, weight_f_ideal, ttt, tmp, tmp2
      real(8),dimension(:,:),allocatable :: weights_c_startend, weights_f_startend
      character(len=*),parameter :: subname='assign_weight_to_process'
    
      ! Ideal weight per process.
      weight_c_ideal=weight_tot_c/dble(nproc)
      weight_f_ideal=weight_tot_f/dble(nproc)
    
    
    
    
    allocate(weights_c_startend(2,0:nproc-1), stat=istat)
    call memocc(istat, weights_c_startend, 'weights_c_startend', subname)
    allocate(weights_f_startend(2,0:nproc-1), stat=istat)
    call memocc(istat, weights_f_startend, 'weights_f_startend', subname)
    
      tt=0.d0
      weights_c_startend(1,0)=0.d0
      do jproc=0,nproc-2
          tt=tt+weight_c_ideal
          weights_c_startend(2,jproc)=dble(floor(tt,kind=8))
          weights_c_startend(1,jproc+1)=dble(floor(tt,kind=8))+1.d0
      end do
      weights_c_startend(2,nproc-1)=weight_tot_c
    
      ! Iterate through all grid points and assign them to processes such that the
      ! load balancing is optimal.
      if (nproc==1) then
          istartend_c(1,0)=1
          istartend_c(2,0)=lzd%glr%wfd%nvctr_c
          weightp_c = weight_tot_c 
          istartp_seg_c=1
          iendp_seg_c=lzd%glr%wfd%nseg_c
          nvalp_c=0
          do i1=0,lzd%glr%d%n1
             do i2=0,lzd%glr%d%n2
                do i3=0,lzd%glr%d%n3
                   nvalp_c = nvalp_c+nint(sqrt(weight_c(i1,i2,i3)))
                end do
             end do
          end do
      else
          tt=0.d0
          tt2=0.d0
          ttt=0.d0
          iitot=0
          jproc=0
          istartend_c(1,0)=1
          if (iproc==0) then
              istartp_seg_c=1
          end if
          loop_nseg_c: do iseg=1,lzd%glr%wfd%nseg_c
              jj=lzd%glr%wfd%keyvloc(iseg)
              j0=lzd%glr%wfd%keygloc(1,iseg)
              j1=lzd%glr%wfd%keygloc(2,iseg)
              ii=j0-1
              i3=ii/((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1))
              ii=ii-i3*(lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)
              i2=ii/(lzd%glr%d%n1+1)
              i0=ii-i2*(lzd%glr%d%n1+1)
              i1=i0+j1-j0
              tmp=0.d0
              tmp2=0.d0
              do i=i0,i1
                  tt=tt+weight_c(i,i2,i3)
                  tt2=tt2+weight_c(i,i2,i3)
                  ttt=ttt+sqrt(weight_c(i,i2,i3))
                  tmp=tmp+weight_c(i,i2,i3)
                  tmp2=tmp2+sqrt(weight_c(i,i2,i3))
                  iitot=iitot+1
                  if (jproc<nproc) then
                      if (tt>weights_c_startend(1,jproc)) then
                          if (jproc>0) then
                              if (iproc==jproc) then
                                  istartp_seg_c=iseg
                              end if
                              if (iproc==jproc-1) then
                                  iendp_seg_c=iseg
                                  weightp_c=tt2
                                  nvalp_c=nint(ttt)
                              end if
                              istartend_c(1,jproc)=iitot+1
                              tt2=0.d0
                              ttt=0.d0
                          end if
                          jproc=jproc+1
                      end if
                  end if
              end do
          end do loop_nseg_c
    
          do jproc=0,nproc-2
              istartend_c(2,jproc)=istartend_c(1,jproc+1)-1
          end do
          istartend_c(2,nproc-1)=lzd%glr%wfd%nvctr_c
    
          if(iproc==nproc-1) then
              nvalp_c=nint(ttt)
              weightp_c=tt2
              iendp_seg_c=lzd%glr%wfd%nseg_c
          end if
          !ii=iendp_seg_c-istartp_seg_c+1
          !call mpiallred(ii, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
          !if (ii/=lzd%glr%wfd%nseg_c) stop 'ii/=lzd%glr%wfd%nseg_c'
      end if
    
    
      ! Same for fine region
      tt=0.d0
      weights_f_startend(1,0)=0.d0
      do jproc=0,nproc-2
          tt=tt+weight_f_ideal
          weights_f_startend(2,jproc)=dble(floor(tt,kind=8))
          weights_f_startend(1,jproc+1)=dble(floor(tt,kind=8))+1.d0
      end do
      weights_f_startend(2,nproc-1)=weight_tot_f
    
    
      istart=lzd%glr%wfd%nseg_c+min(1,lzd%glr%wfd%nseg_f)
      iend=istart+lzd%glr%wfd%nseg_f-1
      if (nproc==1) then
          istartend_f(1,0)=1
          istartend_f(2,0)=lzd%glr%wfd%nvctr_f
          weightp_f = weight_tot_f 
          istartp_seg_f=1+lzd%glr%wfd%nseg_c
          iendp_seg_f=lzd%glr%wfd%nseg_c+lzd%glr%wfd%nseg_f
          nvalp_f=0
          do i1=0,lzd%glr%d%n1
             do i2=0,lzd%glr%d%n2
                do i3=0,lzd%glr%d%n3
                   nvalp_f = nvalp_f+nint(sqrt(weight_f(i1,i2,i3)))
                end do
             end do
          end do
      else
          tt=0.d0
          tt2=0.d0
          ttt=0.d0
          iitot=0
          jproc=0
          istart=lzd%glr%wfd%nseg_c+min(1,lzd%glr%wfd%nseg_f)
          iend=istart+lzd%glr%wfd%nseg_f-1
          istartend_f(1,0)=1
          if (iproc==0) then
              istartp_seg_f=istart
          end if
          loop_nseg_f: do iseg=istart,iend
              jj=lzd%glr%wfd%keyvloc(iseg)
              j0=lzd%glr%wfd%keygloc(1,iseg)
              j1=lzd%glr%wfd%keygloc(2,iseg)
              ii=j0-1
              i3=ii/((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1))
              ii=ii-i3*(lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)
              i2=ii/(lzd%glr%d%n1+1)
              i0=ii-i2*(lzd%glr%d%n1+1)
              i1=i0+j1-j0
              tmp=0.d0
              tmp2=0.d0
              do i=i0,i1
                  tt=tt+weight_f(i,i2,i3)
                  tt2=tt2+weight_f(i,i2,i3)
                  ttt=ttt+sqrt(weight_f(i,i2,i3))
                  tmp=tmp+weight_f(i,i2,i3)
                  tmp2=tmp2+sqrt(weight_c(i,i2,i3))
                  iitot=iitot+1
                  if (jproc<nproc) then
                      if (tt>weights_f_startend(1,jproc)) then
                          if (jproc>0) then
                              if (iproc==jproc) then
                                  istartp_seg_f=iseg
                              end if
                              if (iproc==jproc-1) then
                                  iendp_seg_f=iseg
                                  weightp_f=tt2
                                  nvalp_f=nint(ttt)
                              end if
                              istartend_f(1,jproc)=iitot+1
                              tt2=0.d0
                              ttt=0.d0
                          end if
                          jproc=jproc+1
                      end if
                  end if
              end do
          end do loop_nseg_f
    
          do jproc=0,nproc-2
              istartend_f(2,jproc)=istartend_f(1,jproc+1)-1
          end do
          istartend_f(2,nproc-1)=lzd%glr%wfd%nvctr_f
          if(iproc==nproc-1) then
              nvalp_f=nint(ttt)
              weightp_f=tt2
              iendp_seg_f=iend
          end if
          !ii=iendp_seg_f-istartp_seg_f+1
          !call mpiallred(ii, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
          !if (ii/=lzd%glr%wfd%nseg_f) stop 'ii/=lzd%glr%wfd%nseg_f'
      end if
    
    
    
    
      iall = -product(shape(weights_c_startend))*kind(weights_c_startend)
      deallocate(weights_c_startend,stat=istat)
      call memocc(istat, iall, 'weights_c_startend', subname)
      iall = -product(shape(weights_f_startend))*kind(weights_f_startend)
      deallocate(weights_f_startend,stat=istat)
      call memocc(istat, iall, 'weights_f_startend', subname)
    
      nptsp_c=istartend_c(2,iproc)-istartend_c(1,iproc)+1
      nptsp_f=istartend_f(2,iproc)-istartend_f(1,iproc)+1
    
    
    
      ! some check
      ii_f=istartend_f(2,iproc)-istartend_f(1,iproc)+1
      if(nproc>1) call mpiallred(ii_f, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
      !if(ii_f/=lzd%glr%wfd%nvctr_f) stop 'assign_weight_to_process: ii_f/=lzd%glr%wfd%nvctr_f'
      if(ii_f/=lzd%glr%wfd%nvctr_f) then
         write(*,*) 'ii_f/=lzd%glr%wfd%nvctr_f',ii_f,lzd%glr%wfd%nvctr_f
         if (iproc==0) then
             do jproc=0,nproc-1
                 write(*,*) jproc, istartend_f(1,jproc), istartend_f(2,jproc)
             end do
         end if
         stop
      end if
     
      ii_c=istartend_c(2,iproc)-istartend_c(1,iproc)+1
      if(nproc>1) call mpiallred(ii_c, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
      if(ii_c/=lzd%glr%wfd%nvctr_c) then
         write(*,*) 'ii_c/=lzd%glr%wfd%nvctr_c',ii_c,lzd%glr%wfd%nvctr_c
         stop
      end if
    
    
      ! some checks
      if(nproc>1) then
          call mpi_allreduce(weightp_c, tt, 1, mpi_double_precision, mpi_sum, bigdft_mpi%mpi_comm, ierr)
      else
          tt=weightp_c
      end if
      if(tt/=weight_tot_c) stop 'wrong partition of coarse weights'
      if(nproc>1) then
          call mpi_allreduce(weightp_f, tt, 1, mpi_double_precision, mpi_sum, bigdft_mpi%mpi_comm, ierr)
      else
          tt=weightp_f
      end if     
      if(tt/=weight_tot_f) stop 'wrong partition of fine weights'
      if(nproc>1) then
          call mpi_allreduce(nptsp_c, ii, 1, mpi_integer, mpi_sum, bigdft_mpi%mpi_comm, ierr)
      else
          ii=nptsp_c
      end if
      if(ii/=lzd%glr%wfd%nvctr_c) stop 'wrong partition of coarse grid points'
      if(nproc>1) then
          call mpi_allreduce(nptsp_f, ii, 1, mpi_integer, mpi_sum, bigdft_mpi%mpi_comm, ierr)
      else
          ii=nptsp_f
      end if
      if(ii/=lzd%glr%wfd%nvctr_f) stop 'init_collective_comms: wrong partition of fine grid points'
      
     
    
    end subroutine assign_weight_to_process


    subroutine get_index_in_global2(lr, index_in_global_c, index_in_global_f)
    use module_base
    use module_types
    implicit none
    
    ! Calling arguments
    type(locreg_descriptors),intent(in) :: lr
    integer,dimension(0:lr%d%n1,0:lr%d%n2,0:lr%d%n3),intent(out) :: index_in_global_c, index_in_global_f
    
    ! Local variables
    integer :: iitot, iseg, j0, j1, ii, i1, i2, i3, i0, i, istart, iend
    
    
        iitot=0
        do iseg=1,lr%wfd%nseg_c
           j0=lr%wfd%keygloc(1,iseg)
           j1=lr%wfd%keygloc(2,iseg)
           ii=j0-1
           i3=ii/((lr%d%n1+1)*(lr%d%n2+1))
           ii=ii-i3*(lr%d%n1+1)*(lr%d%n2+1)
           i2=ii/(lr%d%n1+1)
           i0=ii-i2*(lr%d%n1+1)
           i1=i0+j1-j0
           do i=i0,i1
              iitot=iitot+1
              index_in_global_c(i,i2,i3)=iitot
           end do
        end do 
    
    
        iitot=0
        istart=lr%wfd%nseg_c+min(1,lr%wfd%nseg_f)
        iend=istart+lr%wfd%nseg_f-1
        do iseg=istart,iend
           j0=lr%wfd%keygloc(1,iseg)
           j1=lr%wfd%keygloc(2,iseg)
           ii=j0-1
           i3=ii/((lr%d%n1+1)*(lr%d%n2+1))
           ii=ii-i3*(lr%d%n1+1)*(lr%d%n2+1)
           i2=ii/(lr%d%n1+1)
           i0=ii-i2*(lr%d%n1+1)
           i1=i0+j1-j0
           do i=i0,i1
              iitot=iitot+1
              index_in_global_f(i,i2,i3)=iitot
           end do
        end do
    
    
    
    end subroutine get_index_in_global2




    subroutine determine_communication_arrays(iproc, nproc, npsidim_orbs, orbs, lzd, &
               istartend_c, istartend_f, index_in_global_c, index_in_global_f, &
               nvalp_c, nvalp_f,  nsendcounts_c, nsenddspls_c, nrecvcounts_c, nrecvdspls_c, &
               nsendcounts_f, nsenddspls_f, nrecvcounts_f, nrecvdspls_f)
      use module_base
      use module_types
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, npsidim_orbs
      type(orbitals_data),intent(in) :: orbs
      type(local_zone_descriptors),intent(in) :: lzd
      integer,dimension(2,0:nproc-1),intent(in) :: istartend_c, istartend_f
      integer,dimension(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3),intent(in) :: index_in_global_c, index_in_global_f
      integer,intent(in) :: nvalp_c, nvalp_f
      integer,dimension(0:nproc-1),intent(out) :: nsendcounts_c, nsenddspls_c, nrecvcounts_c, nrecvdspls_c
      integer,dimension(0:nproc-1),intent(out) :: nsendcounts_f, nsenddspls_f, nrecvcounts_f, nrecvdspls_f
      
      ! Local variables
      integer :: iorb, iiorb, i1, i2, i3, ii, jproc, jproctarget, ierr, jj, ilr, j0, j1, i0, i, ind
      integer :: istat, ii1, ii2, ii3, iseg, istart, iend, iall
      integer,dimension(:),allocatable :: nsendcounts_tmp, nsenddspls_tmp, nrecvcounts_tmp, nrecvdspls_tmp
      character(len=*),parameter :: subname='determine_communication_arrays'
    
      ! Determine values for mpi_alltoallv
      ! first nsendcounts
      nsendcounts_c=0
      nsendcounts_f=0
    
      !$omp parallel default(private) shared(ilr,nproc,orbs,lzd,index_in_global_c,istartend_c,nsendcounts_c,nsendcounts_f) &
      !$omp shared(istartend_f,index_in_global_f)
    
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          ilr=orbs%inwhichlocreg(iiorb)
          if (lzd%llr(ilr)%wfd%nseg_c>0) then
              !$omp do firstprivate(ilr) reduction(+:nsendcounts_c)
              do iseg=1,lzd%llr(ilr)%wfd%nseg_c
                  jj=lzd%llr(ilr)%wfd%keyvloc(iseg)
                  j0=lzd%llr(ilr)%wfd%keygloc(1,iseg)
                  j1=lzd%llr(ilr)%wfd%keygloc(2,iseg)
                  ii=j0-1
                  i3=ii/((lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1))
                  ii=ii-i3*(lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)
                  i2=ii/(lzd%llr(ilr)%d%n1+1)
                  i0=ii-i2*(lzd%llr(ilr)%d%n1+1)
                  i1=i0+j1-j0
                  ii2=i2+lzd%llr(ilr)%ns2
                  ii3=i3+lzd%llr(ilr)%ns3
                  do i=i0,i1
                      ii1=i+lzd%llr(ilr)%ns1
                      ind=index_in_global_c(ii1,ii2,ii3)
                      jproctarget=-1
                      do jproc=0,nproc-1
                          if(ind>=istartend_c(1,jproc) .and. ind<=istartend_c(2,jproc)) then
                              jproctarget=jproc
                              exit
                          end if
                      end do
                      if (jproctarget /= -1) &
                           nsendcounts_c(jproctarget)=nsendcounts_c(jproctarget)+1
                  end do
              end do
              !$omp end do
          end if
      end do
    
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          ilr=orbs%inwhichlocreg(iiorb)
          istart=lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)
          iend=istart+lzd%llr(ilr)%wfd%nseg_f-1
          if (istart<iend) then
              !$omp do firstprivate(ilr) reduction(+:nsendcounts_f)
              do iseg=istart,iend
                  jj=lzd%llr(ilr)%wfd%keyvloc(iseg)
                  j0=lzd%llr(ilr)%wfd%keygloc(1,iseg)
                  j1=lzd%llr(ilr)%wfd%keygloc(2,iseg)
                  ii=j0-1
                  i3=ii/((lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1))
                  ii=ii-i3*(lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)
                  i2=ii/(lzd%llr(ilr)%d%n1+1)
                  i0=ii-i2*(lzd%llr(ilr)%d%n1+1)
                  i1=i0+j1-j0
                  ii2=i2+lzd%llr(ilr)%ns2
                  ii3=i3+lzd%llr(ilr)%ns3
                  do i=i0,i1
                      ii1=i+lzd%llr(ilr)%ns1
                      !call get_index_in_global(lzd%glr, ii1, ii2, ii3, 'f', ind)
                      ind=index_in_global_f(ii1,ii2,ii3)
                      jproctarget=-1
                      do jproc=0,nproc-1
                          if(ind>=istartend_f(1,jproc) .and. ind<=istartend_f(2,jproc)) then
                              jproctarget=jproc
                              exit
                          end if
                      end do
                      if (jproctarget /= -1) &
                           nsendcounts_f(jproctarget)=nsendcounts_f(jproctarget)+1
                  end do
              end do
              !$omp end do
          end if
       end do
    
       !$omp end parallel
    
    
      ! The first check is to make sure that there is no stop in case this process has no orbitals (in which case
      ! npsidim_orbs is 1 and not 0 as assumed by the check)
      if(npsidim_orbs>1 .and. sum(nsendcounts_c)+7*sum(nsendcounts_f)/=npsidim_orbs) then
          write(*,'(a,2i10)') 'sum(nsendcounts_c)+sum(nsendcounts_f)/=npsidim_orbs', &
                              sum(nsendcounts_c)+sum(nsendcounts_f), npsidim_orbs
          stop
      end if
    
      
      ! now nsenddspls
      nsenddspls_c(0)=0
      do jproc=1,nproc-1
          nsenddspls_c(jproc)=nsenddspls_c(jproc-1)+nsendcounts_c(jproc-1)
      end do
      nsenddspls_f(0)=0
      do jproc=1,nproc-1
          nsenddspls_f(jproc)=nsenddspls_f(jproc-1)+nsendcounts_f(jproc-1)
      end do
    
    
    
      ! now nrecvcounts
      ! use an mpi_alltoallv to gather the data
      allocate(nsendcounts_tmp(0:nproc-1), stat=istat)
      call memocc(istat, nsendcounts_tmp, 'nsendcounts_tmp', subname)
      allocate(nsenddspls_tmp(0:nproc-1), stat=istat)
      call memocc(istat, nsenddspls_tmp, 'nsenddspls_tmp', subname)
      allocate(nrecvcounts_tmp(0:nproc-1), stat=istat)
      call memocc(istat, nrecvcounts_tmp, 'nrecvcounts_tmp', subname)
      allocate(nrecvdspls_tmp(0:nproc-1), stat=istat)
      call memocc(istat, nrecvdspls_tmp, 'nrecvdspls_tmp', subname)
      nsendcounts_tmp=1
      nrecvcounts_tmp=1
      do jproc=0,nproc-1
          nsenddspls_tmp(jproc)=jproc
          nrecvdspls_tmp(jproc)=jproc
      end do
      if(nproc>1) then
          call mpi_alltoallv(nsendcounts_c, nsendcounts_tmp, nsenddspls_tmp, mpi_integer, nrecvcounts_c, &
               nrecvcounts_tmp, nrecvdspls_tmp, mpi_integer, bigdft_mpi%mpi_comm, ierr)
          call mpi_alltoallv(nsendcounts_f, nsendcounts_tmp, nsenddspls_tmp, mpi_integer, nrecvcounts_f, &
               nrecvcounts_tmp, nrecvdspls_tmp, mpi_integer, bigdft_mpi%mpi_comm, ierr)
      else
          nrecvcounts_c=nsendcounts_c
          nrecvcounts_f=nsendcounts_f
      end if
      iall=-product(shape(nsendcounts_tmp))*kind(nsendcounts_tmp)
      deallocate(nsendcounts_tmp, stat=istat)
      call memocc(istat, iall, 'nsendcounts_tmp', subname)
      iall=-product(shape(nsenddspls_tmp))*kind(nsenddspls_tmp)
      deallocate(nsenddspls_tmp, stat=istat)
      call memocc(istat, iall, 'nsenddspls_tmp', subname)
      iall=-product(shape(nrecvcounts_tmp))*kind(nrecvcounts_tmp)
      deallocate(nrecvcounts_tmp, stat=istat)
      call memocc(istat, iall, 'nrecvcounts_tmp', subname)
      iall=-product(shape(nrecvdspls_tmp))*kind(nrecvdspls_tmp)
      deallocate(nrecvdspls_tmp, stat=istat)
      call memocc(istat, iall, 'nrecvdspls_tmp', subname)
    
      ! now recvdspls
      nrecvdspls_c(0)=0
      do jproc=1,nproc-1
          nrecvdspls_c(jproc)=nrecvdspls_c(jproc-1)+nrecvcounts_c(jproc-1)
      end do
      nrecvdspls_f(0)=0
      do jproc=1,nproc-1
          nrecvdspls_f(jproc)=nrecvdspls_f(jproc-1)+nrecvcounts_f(jproc-1)
      end do
    
      if(sum(nrecvcounts_c)/=nvalp_c) stop 'sum(nrecvcounts_c)/=nvalp_c'
      if(sum(nrecvcounts_f)/=nvalp_f) stop 'sum(nrecvcounts_f)/=nvalp_f'
    
    end subroutine determine_communication_arrays


    subroutine determine_num_orbs_per_gridpoint_new(iproc, nproc, lzd, istartend_c, istartend_f, &
               istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, &
               weightp_c, weightp_f, nptsp_c, nptsp_f, weight_c, weight_f, &
               norb_per_gridpoint_c, norb_per_gridpoint_f)
      use module_base
      use module_types
      implicit none
      
      ! Calling arguments
      integer,intent(in):: iproc, nproc, nptsp_c, nptsp_f, istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f
      type(local_zone_descriptors),intent(in):: lzd
      integer,dimension(2,0:nproc-1),intent(in):: istartend_c, istartend_f
      real(8),intent(in):: weightp_c, weightp_f
      real(8),dimension(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3),intent(in):: weight_c, weight_f
      integer,dimension(nptsp_c),intent(out):: norb_per_gridpoint_c
      integer,dimension(nptsp_f),intent(out):: norb_per_gridpoint_f
      
      ! Local variables
      integer:: ii, i1, i2, i3, iipt, iseg, jj, j0, j1, iitot, i, istart, iend, i0
      integer::icheck_c,icheck_f,iiorb_c,iiorb_f, npgp_c,npgp_f
      !!integer,dimension(:),allocatable:: iseg_start_c, iseg_start_f
    
    
      icheck_c = 0
      icheck_f = 0
      iiorb_f=0
      iiorb_c=0
      iipt=0
    
      !$omp parallel default(private) shared(lzd,iproc,istartend_c,istartend_f,istartp_seg_c,iendp_seg_c,istartp_seg_f,iendp_seg_f) &
      !$omp shared(nptsp_c, weight_c,norb_per_gridpoint_c,weightp_c,nptsp_f, weight_f,norb_per_gridpoint_f,weightp_f) &
      !$omp shared(icheck_f,iiorb_f,icheck_c,iiorb_c)
    
    
      if(istartp_seg_c<=iendp_seg_c) then
          !$omp do reduction(+:icheck_c) reduction(+:iiorb_c)
          do iseg=istartp_seg_c,iendp_seg_c
              jj=lzd%glr%wfd%keyvloc(iseg)
              j0=lzd%glr%wfd%keygloc(1,iseg)
              j1=lzd%glr%wfd%keygloc(2,iseg)
              ii=j0-1
              i3=ii/((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1))
              ii=ii-i3*(lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)
              i2=ii/(lzd%glr%d%n1+1)
              i0=ii-i2*(lzd%glr%d%n1+1)
              i1=i0+j1-j0
              do i=i0,i1
                  iitot=jj+i-i0
                  if(iitot>=istartend_c(1,iproc) .and. iitot<=istartend_c(2,iproc)) then
                      icheck_c = icheck_c + 1
                      iipt=jj-istartend_c(1,iproc)+i-i0+1
                      npgp_c = nint(sqrt(weight_c(i,i2,i3)))
                      iiorb_c=iiorb_c+nint(weight_c(i,i2,i3))
                      norb_per_gridpoint_c(iipt)=npgp_c
                  end if
              end do
          end do
          !$omp end do
      end if
    
      if(icheck_c/=nptsp_c) stop 'icheck_c/=nptsp_c'
      if(iiorb_c/=nint(weightp_c)) stop 'iiorb_c/=weightp_c'
    
    
      iitot=0
      iipt=0
      istart=lzd%glr%wfd%nseg_c+min(1,lzd%glr%wfd%nseg_f)
      iend=istart+lzd%glr%wfd%nseg_f-1
    
      if (istartp_seg_f<=iendp_seg_f) then
          !$omp do reduction(+:icheck_f) reduction(+:iiorb_f)
          do iseg=istartp_seg_f,iendp_seg_f
              jj=lzd%glr%wfd%keyvloc(iseg)
              j0=lzd%glr%wfd%keygloc(1,iseg)
              j1=lzd%glr%wfd%keygloc(2,iseg)
              ii=j0-1
              i3=ii/((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1))
              ii=ii-i3*(lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)
              i2=ii/(lzd%glr%d%n1+1)
              i0=ii-i2*(lzd%glr%d%n1+1)
              i1=i0+j1-j0
              do i=i0,i1
                  iitot=jj+i-i0
                  if(iitot>=istartend_f(1,iproc) .and. iitot<=istartend_f(2,iproc)) then
                      icheck_f = icheck_f +1
                      iipt=jj-istartend_f(1,iproc)+i-i0+1
                      npgp_f=0
                      npgp_f = nint(sqrt(weight_f(i,i2,i3)))
                      iiorb_f=iiorb_f+nint(weight_f(i,i2,i3))
                      norb_per_gridpoint_f(iipt)=npgp_f
                  end if
              end do
          end do
          !$omp end do
      end if
    
      !$omp end parallel
    
      if(icheck_f/=nptsp_f) stop 'icheck_f/=nptsp_f'
      if(iiorb_f/=nint(weightp_f)) stop 'iiorb_f/=weightp_f'
    
    end subroutine determine_num_orbs_per_gridpoint_new




    subroutine get_switch_indices(iproc, nproc, orbs, lzd, ndimpsi_c, ndimpsi_f, istartend_c, istartend_f, &
               nsendcounts_c, nsenddspls_c, ndimind_c, nrecvcounts_c, nrecvdspls_c, &
               nsendcounts_f, nsenddspls_f, ndimind_f, nrecvcounts_f, nrecvdspls_f, &
               index_in_global_c, index_in_global_f, &
               weightp_c, weightp_f,  isendbuf_c, irecvbuf_c, isendbuf_f, irecvbuf_f, &
               indexrecvorbital_c, iextract_c, iexpand_c, indexrecvorbital_f, iextract_f, iexpand_f)
      use module_base
      use module_types
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, ndimpsi_c, ndimpsi_f, ndimind_c,ndimind_f
      type(orbitals_data),intent(in) :: orbs
      type(local_zone_descriptors),intent(in) :: lzd
      integer,dimension(2,0:nproc-1),intent(in) :: istartend_c, istartend_f
      integer,dimension(0:nproc-1),intent(in) :: nsendcounts_c, nsenddspls_c, nrecvcounts_c, nrecvdspls_c
      integer,dimension(0:nproc-1),intent(in) :: nsendcounts_f, nsenddspls_f, nrecvcounts_f, nrecvdspls_f
      integer,dimension(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3),intent(in) :: index_in_global_c, index_in_global_f
      real(kind=8),intent(in) :: weightp_c, weightp_f
      integer,dimension(ndimpsi_c),intent(out) :: isendbuf_c, irecvbuf_c
      integer,dimension(ndimpsi_f),intent(out) :: isendbuf_f, irecvbuf_f
      integer,dimension(ndimind_c),intent(out) :: indexrecvorbital_c, iextract_c, iexpand_c
      integer,dimension(ndimind_f),intent(out) :: indexrecvorbital_f, iextract_f, iexpand_f
      
      ! Local variables
      integer :: i, iorb, iiorb, i1, i2, i3, ind, jproc, jproctarget, ii, ierr, jj, iseg, iitot, ilr
      integer :: istart, iend, indglob, ii1, ii2, ii3, j1, i0, j0, istat, iall
      integer,dimension(:),allocatable :: nsend_c,nsend_f, indexsendorbital2, indexrecvorbital2
      integer,dimension(:),allocatable :: gridpoint_start_c, gridpoint_start_f
      real(kind=8),dimension(:,:,:),allocatable :: weight_c, weight_f
      integer,dimension(:),allocatable :: indexsendorbital_c, indexsendbuf_c, indexrecvbuf_c
      integer,dimension(:),allocatable :: indexsendorbital_f, indexsendbuf_f, indexrecvbuf_f
      character(len=*),parameter :: subname='get_switch_indices'
    
    
      
      allocate(indexsendorbital_c(ndimpsi_c), stat=istat)
      call memocc(istat, indexsendorbital_c, 'indexsendorbital_c', subname)
      allocate(indexsendbuf_c(ndimpsi_c), stat=istat)
      call memocc(istat, indexsendbuf_c, 'indexsendbuf_c', subname)
      allocate(indexrecvbuf_c(sum(nrecvcounts_c)), stat=istat)
      call memocc(istat, indexrecvbuf_c, 'indexrecvbuf_c', subname)
      
      allocate(indexsendorbital_f(ndimpsi_f), stat=istat)
      call memocc(istat, indexsendorbital_f, 'indexsendorbital_f', subname)
      allocate(indexsendbuf_f(ndimpsi_f), stat=istat)
      call memocc(istat, indexsendbuf_f, 'indexsendbuf_f', subname)
      allocate(indexrecvbuf_f(sum(nrecvcounts_f)), stat=istat)
      call memocc(istat, indexrecvbuf_f, 'indexrecvbuf_f', subname)
      
      allocate(weight_c(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3), stat=istat)
      call memocc(istat, weight_c, 'weight_c', subname)
      allocate(weight_f(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3), stat=istat)
      call memocc(istat, weight_f, 'weight_f', subname)
      allocate(gridpoint_start_c((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(lzd%glr%d%n3+1)), stat=istat)
      call memocc(istat, gridpoint_start_c, 'gridpoint_start_c', subname)
      allocate(gridpoint_start_f((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(lzd%glr%d%n3+1)), stat=istat)
      call memocc(istat, gridpoint_start_f, 'gridpoint_start_f', subname)
      gridpoint_start_c=-1
      gridpoint_start_f=-1
    
    !write(*,*) 'ndimpsi_f, sum(nrecvcounts_f)', ndimpsi_f, sum(nrecvcounts_f)
    
      allocate(nsend_c(0:nproc-1), stat=istat)
      call memocc(istat, nsend_c, 'nsend_c', subname)
      allocate(nsend_f(0:nproc-1), stat=istat)
      call memocc(istat, nsend_f, 'nsend_f', subname)
    
      nsend_c=0
      nsend_f=0
    
      !$omp parallel default(private) shared(orbs,lzd,index_in_global_c,index_in_global_f,istartend_c,istartend_f)&
      !$omp shared(nsend_c,nsend_f,nsenddspls_c,nsenddspls_f,ndimpsi_c,ndimpsi_f,nsendcounts_c,nsendcounts_f,nproc) &
      !$omp shared(isendbuf_c,isendbuf_f,indexsendbuf_c,indexsendbuf_f,indexsendorbital_c,indexsendorbital_f)
    
      !$omp sections
      !$omp section
      iitot=0
     
      do iorb=1,orbs%norbp
        iiorb=orbs%isorb+iorb
        ilr=orbs%inwhichlocreg(iiorb)
        do iseg=1,lzd%llr(ilr)%wfd%nseg_c
           jj=lzd%llr(ilr)%wfd%keyvloc(iseg)
           j0=lzd%llr(ilr)%wfd%keygloc(1,iseg)
           j1=lzd%llr(ilr)%wfd%keygloc(2,iseg)
           ii=j0-1
           i3=ii/((lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1))
           ii=ii-i3*(lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)
           i2=ii/(lzd%llr(ilr)%d%n1+1)
           i0=ii-i2*(lzd%llr(ilr)%d%n1+1)
           i1=i0+j1-j0
           !write(*,'(a,8i8)') 'jj, ii, j0, j1, i0, i1, i2, i3',jj,ii,j0,j1,i0,i1,i2,i3
           do i=i0,i1
              ii1=i+lzd%llr(ilr)%ns1
              ii2=i2+lzd%llr(ilr)%ns2
              ii3=i3+lzd%llr(ilr)%ns3
              !call get_index_in_global(lzd%glr, ii1, ii2, ii3, 'c', indglob)
              indglob=index_in_global_c(ii1,ii2,ii3)
              iitot=iitot+1
                  jproctarget=-1
                  do jproc=0,nproc-1
                      if(indglob>=istartend_c(1,jproc) .and. indglob<=istartend_c(2,jproc)) then
                          jproctarget=jproc
                          exit
                      end if
                  end do
                  !write(600+iproc,'(a,2(i0,1x),i0,a,i0)') 'point ',ii1,ii2,ii3,' goes to process ',jproctarget
              
                  if (jproctarget/=-1) then
                     nsend_c(jproctarget)=nsend_c(jproctarget)+1
                     ind=nsenddspls_c(jproctarget)+nsend_c(jproctarget)
                     isendbuf_c(iitot)=ind
                     indexsendbuf_c(ind)=indglob
                     indexsendorbital_c(iitot)=iiorb
                  end if
                  !indexsendorbital(ind)=iiorb
              end do
          end do
          
      end do
     ! write(*,*) 'iitot,ndimpsi_c',iitot,ndimpsi_c
      if(iitot/=ndimpsi_c) stop 'iitot/=ndimpsi_c'
    
      !check
      do jproc=0,nproc-1
          if(nsend_c(jproc)/=nsendcounts_c(jproc)) stop 'nsend_c(jproc)/=nsendcounts_c(jproc)'
      end do
    
    
      !$omp section
      ! fine part
      iitot=0
     
      do iorb=1,orbs%norbp
        iiorb=orbs%isorb+iorb
        ilr=orbs%inwhichlocreg(iiorb)
        istart=lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)
        iend=istart+lzd%llr(ilr)%wfd%nseg_f-1
        do iseg=istart,iend
           jj=lzd%llr(ilr)%wfd%keyvloc(iseg)
           j0=lzd%llr(ilr)%wfd%keygloc(1,iseg)
           j1=lzd%llr(ilr)%wfd%keygloc(2,iseg)
           ii=j0-1
           i3=ii/((lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1))
           ii=ii-i3*(lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)
           i2=ii/(lzd%llr(ilr)%d%n1+1)
           i0=ii-i2*(lzd%llr(ilr)%d%n1+1)
           i1=i0+j1-j0
           !write(*,'(a,8i8)') 'jj, ii, j0, j1, i0, i1, i2, i3',jj,ii,j0,j1,i0,i1,i2,i3
           do i=i0,i1
              ii1=i+lzd%llr(ilr)%ns1
              ii2=i2+lzd%llr(ilr)%ns2
              ii3=i3+lzd%llr(ilr)%ns3
              !call get_index_in_global(lzd%glr, ii1, ii2, ii3, 'f', indglob)
              indglob=index_in_global_f(ii1,ii2,ii3)
                      iitot=iitot+1
                      jproctarget=-1
                      do jproc=0,nproc-1
                          if(indglob>=istartend_f(1,jproc) .and. indglob<=istartend_f(2,jproc)) then
                              jproctarget=jproc
                              exit
                          end if
                      end do
                      if (jproctarget/=-1) then
                         nsend_f(jproctarget)=nsend_f(jproctarget)+1
                         ind=nsenddspls_f(jproctarget)+nsend_f(jproctarget)
                         isendbuf_f(iitot)=ind
                         indexsendbuf_f(ind)=indglob
                         indexsendorbital_f(iitot)=iiorb
                      end if
                      !indexsendorbital(ind)=iiorb
              end do
          end do
     
      end do
      
      if(iitot/=ndimpsi_f) stop 'iitot/=ndimpsi_f'
    
      !$omp end sections
      !$omp end parallel
    
      !check
      do jproc=0,nproc-1
          !write(*,*) 'nsend(jproc), nsendcounts_f(jproc)', nsend(jproc), nsendcounts_f(jproc)
          if(nsend_f(jproc)/=nsendcounts_f(jproc)) stop 'nsend_f(jproc)/=nsendcounts_f(jproc)'
      end do
    
      allocate(indexsendorbital2(ndimpsi_c), stat=istat)
      call memocc(istat, indexsendorbital2, 'indexsendorbital2', subname)
      indexsendorbital2=indexsendorbital_c
      do i=1,ndimpsi_c
          ind=isendbuf_c(i)
          indexsendorbital_c(ind)=indexsendorbital2(i)
      end do
    
      ! Inverse of isendbuf
      call get_reverse_indices(ndimpsi_c, isendbuf_c, irecvbuf_c)
    
      iall=-product(shape(indexsendorbital2))*kind(indexsendorbital2)
      deallocate(indexsendorbital2, stat=istat)
      call memocc(istat, iall, 'indexsendorbital2', subname)
    
    
      allocate(indexsendorbital2(ndimpsi_f), stat=istat)
      call memocc(istat, indexsendorbital2, 'indexsendorbital2', subname)
      indexsendorbital2=indexsendorbital_f
      do i=1,ndimpsi_f
          ind=isendbuf_f(i)
          indexsendorbital_f(ind)=indexsendorbital2(i)
      end do
    
      ! Inverse of isendbuf
    
      call get_reverse_indices(ndimpsi_f, isendbuf_f, irecvbuf_f)
      iall=-product(shape(indexsendorbital2))*kind(indexsendorbital2)
      deallocate(indexsendorbital2, stat=istat)
      call memocc(istat, iall, 'indexsendorbital2', subname)
    
    
    
    
      if(nproc>1) then
          ! Communicate indexsendbuf
          call mpi_alltoallv(indexsendbuf_c, nsendcounts_c, nsenddspls_c, mpi_integer, indexrecvbuf_c, &
               nrecvcounts_c, nrecvdspls_c, mpi_integer, bigdft_mpi%mpi_comm, ierr)
          ! Communicate indexsendorbitals
          call mpi_alltoallv(indexsendorbital_c, nsendcounts_c, nsenddspls_c, mpi_integer, indexrecvorbital_c, &
               nrecvcounts_c, nrecvdspls_c, mpi_integer, bigdft_mpi%mpi_comm, ierr)
    
          ! Communicate indexsendbuf
          call mpi_alltoallv(indexsendbuf_f, nsendcounts_f, nsenddspls_f, mpi_integer, indexrecvbuf_f, &
               nrecvcounts_f, nrecvdspls_f, mpi_integer, bigdft_mpi%mpi_comm, ierr)
          ! Communicate indexsendorbitals
          call mpi_alltoallv(indexsendorbital_f, nsendcounts_f, nsenddspls_f, mpi_integer, indexrecvorbital_f, &
               nrecvcounts_f, nrecvdspls_f, mpi_integer, bigdft_mpi%mpi_comm, ierr)
       else
           indexrecvbuf_c=indexsendbuf_c
           indexrecvorbital_c=indexsendorbital_c
           indexrecvbuf_f=indexsendbuf_f
           indexrecvorbital_f=indexsendorbital_f
       end if
    
    
    
      !call get_gridpoint_start(iproc, nproc, norb, glr, llr, nrecvcounts, indexrecvbuf, weight, gridpoint_start)
      call get_gridpoint_start(iproc, nproc, lzd, sum(nrecvcounts_c), nrecvcounts_c, sum(nrecvcounts_f), &
                nrecvcounts_f, indexrecvbuf_c, indexrecvbuf_f, weight_c, weight_f, gridpoint_start_c, gridpoint_start_f)
    
    
    
      if(maxval(gridpoint_start_c)>sum(nrecvcounts_c)) stop '1: maxval(gridpoint_start_c)>sum(nrecvcounts_c)'
      if(maxval(gridpoint_start_f)>sum(nrecvcounts_f)) stop '1: maxval(gridpoint_start_f)>sum(nrecvcounts_f)'
      ! Rearrange the communicated data
      do i=1,sum(nrecvcounts_c)
          ii=indexrecvbuf_c(i)
          ind=gridpoint_start_c(ii)
          !if(ind==0) stop 'ind is zero!'
          iextract_c(i)=ind
          gridpoint_start_c(ii)=gridpoint_start_c(ii)+1  
      end do
      !write(*,'(a,2i12)') 'sum(iextract_c), nint(weightp_c*(weightp_c+1.d0)*.5d0)', sum(iextract_c), nint(weightp_c*(weightp_c+1.d0)*.5d0)
      !if(sum(iextract_c)/=nint(weightp_c*(weightp_c+1.d0)*.5d0)) stop 'sum(iextract_c)/=nint(weightp_c*(weightp_c+1.d0)*.5d0)'
      if(maxval(iextract_c)>sum(nrecvcounts_c)) stop 'maxval(iextract_c)>sum(nrecvcounts_c)'
      if(minval(iextract_c)<1) stop 'minval(iextract_c)<1'
    
      ! Rearrange the communicated data
      iextract_f = 0
      do i=1,sum(nrecvcounts_f)
          ii=indexrecvbuf_f(i)
          ind=gridpoint_start_f(ii)
          !if(ind==0) stop 'ind is zero!'
          iextract_f(i)=ind
          gridpoint_start_f(ii)=gridpoint_start_f(ii)+1  
      end do
      if(maxval(iextract_f)>sum(nrecvcounts_f)) stop 'maxval(iextract_f)>sum(nrecvcounts_f)'
      if(minval(iextract_f)<1) stop 'minval(iextract_f)<1'
    
    
    
    
      ! Get the array to transfrom back the data
      call get_reverse_indices(sum(nrecvcounts_c), iextract_c, iexpand_c)
      call get_reverse_indices(sum(nrecvcounts_f), iextract_f, iexpand_f)
      
    
    
    
      allocate(indexrecvorbital2(sum(nrecvcounts_c)), stat=istat)
      call memocc(istat, indexrecvorbital2, 'indexrecvorbital2', subname)
      indexrecvorbital2=indexrecvorbital_c
      do i=1,sum(nrecvcounts_c)
          ind=iextract_c(i)
          indexrecvorbital_c(ind)=indexrecvorbital2(i)
      end do
      iall=-product(shape(indexrecvorbital2))*kind(indexrecvorbital2)
      deallocate(indexrecvorbital2, stat=istat)
      call memocc(istat, iall, 'indexrecvorbital2', subname)
    
      allocate(indexrecvorbital2(sum(nrecvcounts_f)), stat=istat)
      call memocc(istat, indexrecvorbital2, 'indexrecvorbital2', subname)
      indexrecvorbital2=indexrecvorbital_f
      do i=1,sum(nrecvcounts_f)
          ind=iextract_f(i)
          indexrecvorbital_f(ind)=indexrecvorbital2(i)
      end do
      iall=-product(shape(indexrecvorbital2))*kind(indexrecvorbital2)
      deallocate(indexrecvorbital2, stat=istat)
      call memocc(istat, iall, 'indexrecvorbital2', subname)
    
    
      if(minval(indexrecvorbital_c)<1) stop 'minval(indexrecvorbital_c)<1'
      if(maxval(indexrecvorbital_c)>orbs%norb) stop 'maxval(indexrecvorbital_c)>orbs%norb'
      if(minval(indexrecvorbital_f)<1) stop 'minval(indexrecvorbital_f)<1'
      if(maxval(indexrecvorbital_f)>orbs%norb) stop 'maxval(indexrecvorbital_f)>orbs%norb'
    
    
    
      iall=-product(shape(indexsendorbital_c))*kind(indexsendorbital_c)
      deallocate(indexsendorbital_c, stat=istat)
      call memocc(istat, iall, 'indexsendorbital_c', subname)
      iall=-product(shape(indexsendbuf_c))*kind(indexsendbuf_c)
      deallocate(indexsendbuf_c, stat=istat)
      call memocc(istat, iall, 'indexsendbuf_c', subname)
      iall=-product(shape(indexrecvbuf_c))*kind(indexrecvbuf_c)
      deallocate(indexrecvbuf_c, stat=istat)
      call memocc(istat, iall, 'indexrecvbuf_c', subname)
    
      iall=-product(shape(indexsendorbital_f))*kind(indexsendorbital_f)
      deallocate(indexsendorbital_f, stat=istat)
      call memocc(istat, iall, 'indexsendorbital_f', subname)
      iall=-product(shape(indexsendbuf_f))*kind(indexsendbuf_f)
      deallocate(indexsendbuf_f, stat=istat)
      call memocc(istat, iall, 'indexsendbuf_f', subname)
      iall=-product(shape(indexrecvbuf_f))*kind(indexrecvbuf_f)
      deallocate(indexrecvbuf_f, stat=istat)
      call memocc(istat, iall, 'indexrecvbuf_f', subname)
    
      iall=-product(shape(weight_c))*kind(weight_c)
      deallocate(weight_c, stat=istat)
      call memocc(istat, iall, 'weight_c', subname)
      iall=-product(shape(weight_f))*kind(weight_f)
      deallocate(weight_f, stat=istat)
      call memocc(istat, iall, 'weight_f', subname)
    
      iall=-product(shape(gridpoint_start_c))*kind(gridpoint_start_c)
      deallocate(gridpoint_start_c, stat=istat)
      call memocc(istat, iall, 'gridpoint_start_c', subname)
      iall=-product(shape(gridpoint_start_f))*kind(gridpoint_start_f)
      deallocate(gridpoint_start_f, stat=istat)
      call memocc(istat, iall, 'gridpoint_start_f', subname)
    
      iall=-product(shape(nsend_c))*kind(nsend_c)
      deallocate(nsend_c, stat=istat)
      call memocc(istat, iall, 'nsend_c', subname)
    
       iall=-product(shape(nsend_f))*kind(nsend_f)
      deallocate(nsend_f, stat=istat)
      call memocc(istat, iall, 'nsend_f', subname)
    
    end subroutine get_switch_indices



    subroutine get_reverse_indices(n, indices, reverse_indices)
      use module_base
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: n
      integer,dimension(n),intent(in) :: indices
      integer,dimension(n),intent(out) :: reverse_indices
    
      ! Local variables
      integer :: i, j, m, j0, j1, j2, j3
    
      !$omp parallel default(private) &
      !$omp shared(n, m, indices, reverse_indices)
    
      m=mod(n,4)
      if (m/=0) then
          do i=1,m
              j=indices(i)
              reverse_indices(j)=i
          end do
      end if
    
      !$omp do
      do i=m+1,n,4
          j0=indices(i+0)
          reverse_indices(j0)=i+0
          j1=indices(i+1)
          reverse_indices(j1)=i+1
          j2=indices(i+2)
          reverse_indices(j2)=i+2
          j3=indices(i+3)
          reverse_indices(j3)=i+3
      end do
      !$omp end do
    
      !$omp end parallel
    
      !!do i=1,n
      !!    j=indices(i)
      !!    reverse_indices(j)=i
      !!end do
    
    end subroutine get_reverse_indices



    subroutine get_gridpoint_start(iproc, nproc, lzd, ndimind_c, nrecvcounts_c, ndimind_f, nrecvcounts_f, &
               indexrecvbuf_c, indexrecvbuf_f, weight_c, weight_f, gridpoint_start_c, gridpoint_start_f)
      use module_base
      use module_types
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc,ndimind_c,ndimind_f
      type(local_zone_descriptors),intent(in) :: lzd
      integer,dimension(0:nproc-1),intent(in) :: nrecvcounts_c, nrecvcounts_f
      integer,dimension(ndimind_c),intent(in) :: indexrecvbuf_c
      integer,dimension(ndimind_f),intent(in) :: indexrecvbuf_f
      real(kind=8),dimension(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3),intent(out) :: weight_c, weight_f
      integer,dimension((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(lzd%glr%d%n3+1)),intent(out) :: gridpoint_start_c, gridpoint_start_f
      
      ! Local variables
      integer :: i, ii, jj, i1, i2, i3
    
    
      !weight_c=0.d0
      call to_zero((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(lzd%glr%d%n3+1), weight_c(0,0,0))
      call to_zero((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(lzd%glr%d%n3+1), weight_f(0,0,0))
    
      !$omp parallel default(private) shared(lzd,nrecvcounts_c,indexrecvbuf_c,weight_c,gridpoint_start_c) &
      !$omp shared(nrecvcounts_f,indexrecvbuf_f,weight_f,gridpoint_start_f)
    
      !$omp sections
      !$omp section
      do i=1,sum(nrecvcounts_c)
          ii=indexrecvbuf_c(i)
          !write(650+iproc,*) i, ii
          jj=ii-1
          i3=jj/((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1))
          jj=jj-i3*(lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)
          i2=jj/(lzd%glr%d%n1+1)
          i1=jj-i2*(lzd%glr%d%n1+1)
          weight_c(i1,i2,i3)=weight_c(i1,i2,i3)+1.d0
      end do
    
      !write(*,*) 'in get_gridpoint_start: maxval(weight_c)', maxval(weight_c)
    
      ii=1
      i=0
      !gridpoint_start_c=0
      do i3=0,lzd%glr%d%n3
          do i2=0,lzd%glr%d%n2
              do i1=0,lzd%glr%d%n1
                  i=i+1
                  if(weight_c(i1,i2,i3)>0.d0) then
                      gridpoint_start_c(i)=ii
                      ii=ii+nint(weight_c(i1,i2,i3))
                  else
                      gridpoint_start_c(i) = 0
                  end if
              end do
          end do
      end do
    
      !$omp section
     
      do i=1,sum(nrecvcounts_f)
          ii=indexrecvbuf_f(i)
          jj=ii-1
          i3=jj/((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1))
          jj=jj-i3*(lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)
          i2=jj/(lzd%glr%d%n1+1)
          i1=jj-i2*(lzd%glr%d%n1+1)
          weight_f(i1,i2,i3)=weight_f(i1,i2,i3)+1.d0
      end do
    
    
      ii=1
      i=0
      !gridpoint_start_f=0
      do i3=0,lzd%glr%d%n3
          do i2=0,lzd%glr%d%n2
              do i1=0,lzd%glr%d%n1
                  i=i+1
                  if(weight_f(i1,i2,i3)>0.d0) then
                      gridpoint_start_f(i)=ii
                      ii=ii+nint(weight_f(i1,i2,i3))
                  else
                      gridpoint_start_f(i)=0
                  end if
              end do
          end do
      end do
    
      !$omp end sections
      !$omp end parallel
    
    
    end subroutine get_gridpoint_start




    subroutine transpose_switch_psi(npsidim_orbs, orbs, collcom, psi, psiwork_c, psiwork_f, lzd)
      use module_base
      use module_types
      implicit none
      
      ! Calling arguments
      integer, intent(in) :: npsidim_orbs
      type(orbitals_Data),intent(in) :: orbs
      type(collective_comms),intent(in) :: collcom
      real(kind=8),dimension(npsidim_orbs),intent(in) :: psi
      real(kind=8),dimension(collcom%ndimpsi_c),intent(out) :: psiwork_c
      real(kind=8),dimension(7*collcom%ndimpsi_f),intent(out) :: psiwork_f
      type(local_zone_descriptors),intent(in),optional :: lzd
      
      ! Local variables
      integer :: i_tot, i_c, i_f, iorb, iiorb, ilr, i, ind, istat, iall,m
      real(kind=8),dimension(:),allocatable :: psi_c, psi_f
      character(len=*),parameter :: subname='transpose_switch_psi'
    
      allocate(psi_c(collcom%ndimpsi_c), stat=istat)
      call memocc(istat, psi_c, 'psi_c', subname)
      allocate(psi_f(7*collcom%ndimpsi_f), stat=istat)
      call memocc(istat, psi_f, 'psi_f', subname)
    
    
      if(present(lzd)) then
      ! split up psi into coarse and fine part
    
      
          i_tot=0
          i_c=0
          i_f=0
    
          do iorb=1,orbs%norbp
             iiorb=orbs%isorb+iorb
             ilr=orbs%inwhichlocreg(iiorb)
    
             call vcopy(lzd%llr(ilr)%wfd%nvctr_c,psi(i_tot+1),1,psi_c(i_c+1),1)
    
             i_c = i_c + lzd%llr(ilr)%wfd%nvctr_c
             i_tot = i_tot + lzd%llr(ilr)%wfd%nvctr_c
    
             call vcopy(7*lzd%llr(ilr)%wfd%nvctr_f,psi(i_tot+1),1,psi_f(i_f+1),1)
    
             i_f = i_f + 7*lzd%llr(ilr)%wfd%nvctr_f
             i_tot = i_tot + 7*lzd%llr(ilr)%wfd%nvctr_f
    
          end do
        
    
      else
          ! only coarse part is used...
          call vcopy(collcom%ndimpsi_c, psi(1), 1, psi_c(1), 1)
      end if
    
      ! coarse part
    
      !$omp parallel default(private) &
      !$omp shared(collcom, psi, psiwork_c, psiwork_f, lzd, psi_c,psi_f,m)
    
      m = mod(collcom%ndimpsi_c,7)
      if(m/=0) then
         do i=1,m
            ind = collcom%isendbuf_c(i)
            psiwork_c(ind) = psi_c(i)
         end do
      end if
      !$omp do
      do i = m+1,collcom%ndimpsi_c,7
         psiwork_c(collcom%isendbuf_c(i+0)) = psi_c(i+0)
         psiwork_c(collcom%isendbuf_c(i+1)) = psi_c(i+1)
         psiwork_c(collcom%isendbuf_c(i+2)) = psi_c(i+2)
         psiwork_c(collcom%isendbuf_c(i+3)) = psi_c(i+3)
         psiwork_c(collcom%isendbuf_c(i+4)) = psi_c(i+4)
         psiwork_c(collcom%isendbuf_c(i+5)) = psi_c(i+5)
         psiwork_c(collcom%isendbuf_c(i+6)) = psi_c(i+6)
      end do
      !$omp end do
     
      ! fine part
    
      !$omp do
       do i=1,collcom%ndimpsi_f
           ind=collcom%isendbuf_f(i)
          psiwork_f(7*ind-6)=psi_f(7*i-6)
          psiwork_f(7*ind-5)=psi_f(7*i-5)
          psiwork_f(7*ind-4)=psi_f(7*i-4)
          psiwork_f(7*ind-3)=psi_f(7*i-3)
          psiwork_f(7*ind-2)=psi_f(7*i-2)
          psiwork_f(7*ind-1)=psi_f(7*i-1)
          psiwork_f(7*ind-0)=psi_f(7*i-0)
      end do
      !$omp end do
      !$omp end parallel
    
    
      iall=-product(shape(psi_c))*kind(psi_c)
      deallocate(psi_c, stat=istat)
      call memocc(istat, iall, 'psi_c', subname)
      iall=-product(shape(psi_f))*kind(psi_f)
      deallocate(psi_f, stat=istat)
      call memocc(istat, iall, 'psi_f', subname)
      
    end subroutine transpose_switch_psi



    subroutine transpose_communicate_psi(iproc, nproc, collcom, psiwork_c, psiwork_f, psitwork_c, psitwork_f)
      use module_base
      use module_types
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(collective_comms),intent(in) :: collcom
      real(kind=8),dimension(collcom%ndimpsi_c),intent(in) :: psiwork_c
      real(kind=8),dimension(7*collcom%ndimpsi_f),intent(in) :: psiwork_f
      !real(kind=8),dimension(sum(collcom%nrecvcounts_c)),intent(out) :: psitwork_c
      !real(kind=8),dimension(7*sum(collcom%nrecvcounts_f)),intent(out) :: psitwork_f
      real(kind=8),dimension(collcom%ndimind_c),intent(out) :: psitwork_c
      real(kind=8),dimension(7*collcom%ndimind_f),intent(out) :: psitwork_f
      
      ! Local variables
      integer :: ierr, istat, iall
      !!integer :: iisend, iirecv, ist, ist_c, ist_f, jproc
      real(kind=8),dimension(:),allocatable :: psiwork, psitwork
      integer,dimension(:),allocatable :: nsendcounts, nsenddspls, nrecvcounts, nrecvdspls
      character(len=*),parameter :: subname='transpose_communicate_psi'
    
      !call mpi_comm_size(bigdft_mpi%mpi_comm, nproc, ierr)
      !call mpi_comm_rank(bigdft_mpi%mpi_comm, iproc, ierr)
    
      allocate(psiwork(collcom%ndimpsi_c+7*collcom%ndimpsi_f), stat=istat)
      call memocc(istat, psiwork, 'psiwork', subname)
      allocate(psitwork(sum(collcom%nrecvcounts_c)+7*sum(collcom%nrecvcounts_f)), stat=istat)
      call memocc(istat, psitwork, 'psitwork', subname)
      allocate(nsendcounts(0:nproc-1), stat=istat)
      call memocc(istat, nsendcounts, 'nsendcounts', subname)
      allocate(nsenddspls(0:nproc-1), stat=istat)
      call memocc(istat, nsenddspls, 'nsenddspls', subname)
      allocate(nrecvcounts(0:nproc-1), stat=istat)
      call memocc(istat, nrecvcounts, 'nrecvcounts', subname)
      allocate(nrecvdspls(0:nproc-1), stat=istat)
      call memocc(istat, nrecvdspls, 'nrecvdspls', subname)
    
      !!ist=1
      !!ist_c=1
      !!ist_f=1
      !!iisend=0
      !!iirecv=0
      !!do jproc=0,nproc-1
      !!    if(collcom%nsendcounts_c(jproc)>0) call vcopy(collcom%nsendcounts_c(jproc), psiwork_c(ist_c), 1, psiwork(ist), 1)
      !!    ist_c=ist_c+collcom%nsendcounts_c(jproc)
      !!    ist=ist+collcom%nsendcounts_c(jproc)
      !!    if(collcom%nsendcounts_f(jproc)>0) call vcopy(7*collcom%nsendcounts_f(jproc), psiwork_f(ist_f), 1, psiwork(ist), 1)
      !!    ist_f=ist_f+7*collcom%nsendcounts_f(jproc)
      !!    ist=ist+7*collcom%nsendcounts_f(jproc)
      !!    nsendcounts(jproc)=collcom%nsendcounts_c(jproc)+7*collcom%nsendcounts_f(jproc)
      !!    nsenddspls(jproc)=iisend
      !!    nrecvcounts(jproc)=collcom%nrecvcounts_c(jproc)+7*collcom%nrecvcounts_f(jproc)
      !!    nrecvdspls(jproc)=iirecv
      !!    iisend=iisend+nsendcounts(jproc)
      !!    iirecv=iirecv+nrecvcounts(jproc)
      !!end do
    
      !write(*,'(a,i4,4x,100i8)') 'iproc, nsendcounts', iproc, nsendcounts
      !write(*,'(a,i4,4x,100i8)') 'iproc, nsenddspls', iproc, nsenddspls
      !write(*,'(a,i4,4x,100i8)') 'iproc, nrecvcounts', iproc, nrecvcounts
      !write(*,'(a,i4,4x,100i8)') 'iproc, nrecvdspls', iproc, nrecvdspls
      
      ! coarse part
      call mpi_alltoallv(psiwork_c, collcom%nsendcounts_c, collcom%nsenddspls_c, mpi_double_precision, psitwork_c, &
           collcom%nrecvcounts_c, collcom%nrecvdspls_c, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
      
      ! fine part
      call mpi_alltoallv(psiwork_f, 7*collcom%nsendcounts_f, 7*collcom%nsenddspls_f, mpi_double_precision, psitwork_f, &
           7*collcom%nrecvcounts_f, 7*collcom%nrecvdspls_f, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
      !!call mpi_alltoallv(psiwork, nsendcounts, nsenddspls, mpi_double_precision, psitwork, &
      !!     nrecvcounts, nrecvdspls, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
    
      !!ist=1
      !!ist_c=1
      !!ist_f=1
      !!do jproc=0,nproc-1
      !!    if(collcom%nrecvcounts_c(jproc)>0) call vcopy(collcom%nrecvcounts_c(jproc), psitwork(ist), 1, psitwork_c(ist_c), 1)
      !!    ist_c=ist_c+collcom%nrecvcounts_c(jproc)
      !!    ist=ist+collcom%nrecvcounts_c(jproc)
      !!    if(collcom%nrecvcounts_f(jproc)>0) call vcopy(7*collcom%nrecvcounts_f(jproc), psitwork(ist), 1, psitwork_f(ist_f), 1)
      !!    ist_f=ist_f+7*collcom%nrecvcounts_f(jproc)
      !!    ist=ist+7*collcom%nrecvcounts_f(jproc)
      !!end do
    
      iall=-product(shape(psiwork))*kind(psiwork)
      deallocate(psiwork, stat=istat)
      call memocc(istat, iall, 'psiwork', subname)
      iall=-product(shape(psitwork))*kind(psitwork)
      deallocate(psitwork, stat=istat)
      call memocc(istat, iall, 'psitwork', subname)
      iall=-product(shape(nsendcounts))*kind(nsendcounts)
      deallocate(nsendcounts, stat=istat)
      call memocc(istat, iall, 'nsendcounts', subname)
      iall=-product(shape(nsenddspls))*kind(nsenddspls)
      deallocate(nsenddspls, stat=istat)
      call memocc(istat, iall, 'nsenddspls', subname)
      iall=-product(shape(nrecvcounts))*kind(nrecvcounts)
      deallocate(nrecvcounts, stat=istat)
      call memocc(istat, iall, 'nrecvcounts', subname)
      iall=-product(shape(nrecvdspls))*kind(nrecvdspls)
      deallocate(nrecvdspls, stat=istat)
      call memocc(istat, iall, 'nrecvdspls', subname)
    
    
    end subroutine transpose_communicate_psi



    subroutine transpose_unswitch_psit(collcom, psitwork_c, psitwork_f, psit_c, psit_f)
      use module_base
      use module_types
      implicit none
      
      ! Calling arguments
      type(collective_comms),intent(in) :: collcom
      real(kind=8),dimension(collcom%ndimind_c),intent(in) :: psitwork_c
      real(kind=8),dimension(7*collcom%ndimind_f),intent(in) :: psitwork_f
      real(kind=8),dimension(collcom%ndimind_c),intent(out) :: psit_c
      real(kind=8),dimension(7*collcom%ndimind_f),intent(out) :: psit_f
      
      ! Local variables
      integer :: i, ind, sum_c,sum_f,m
    
      sum_c = sum(collcom%nrecvcounts_c)
      sum_f = sum(collcom%nrecvcounts_f)
    
      !$omp parallel private(i,ind) &
      !$omp shared(psit_c,psit_f, psitwork_c, psitwork_f,collcom,sum_c,sum_f,m)
    
    
      m = mod(sum_c,7)
    
      if(m/=0) then
        do i = 1,m
          ind=collcom%iextract_c(i)
          psit_c(ind)=psitwork_c(i)
        end do
      end if
    
      ! coarse part
    
      !$omp do
      do i=m+1, sum_c,7
          psit_c(collcom%iextract_c(i+0))=psitwork_c(i+0)
          psit_c(collcom%iextract_c(i+1))=psitwork_c(i+1)
          psit_c(collcom%iextract_c(i+2))=psitwork_c(i+2)
          psit_c(collcom%iextract_c(i+3))=psitwork_c(i+3)
          psit_c(collcom%iextract_c(i+4))=psitwork_c(i+4)
          psit_c(collcom%iextract_c(i+5))=psitwork_c(i+5)
          psit_c(collcom%iextract_c(i+6))=psitwork_c(i+6)
      end do
      !$omp end do
    
      ! fine part
    
      !$omp do
      do i=1,sum_f
          ind=collcom%iextract_f(i)
          psit_f(7*ind-6)=psitwork_f(7*i-6)
          psit_f(7*ind-5)=psitwork_f(7*i-5)
          psit_f(7*ind-4)=psitwork_f(7*i-4)
          psit_f(7*ind-3)=psitwork_f(7*i-3)
          psit_f(7*ind-2)=psitwork_f(7*i-2)
          psit_f(7*ind-1)=psitwork_f(7*i-1)
          psit_f(7*ind-0)=psitwork_f(7*i-0)
      end do
      !$omp end do
      
      !$omp end parallel
    
    end subroutine transpose_unswitch_psit



    subroutine transpose_switch_psit(collcom, psit_c, psit_f, psitwork_c, psitwork_f)
      use module_base
      use module_types
      implicit none
    
      ! Calling arguments
      type(collective_comms),intent(in) :: collcom
      real(kind=8),dimension(collcom%ndimind_c),intent(in) :: psit_c
      real(kind=8),dimension(7*collcom%ndimind_f),intent(in) :: psit_f
      real(kind=8),dimension(collcom%ndimind_c),intent(out) :: psitwork_c
      real(kind=8),dimension(7*collcom%ndimind_f),intent(out) :: psitwork_f
      
      ! Local variables
      integer :: i, ind, sum_c,sum_f,m
    
      sum_c = sum(collcom%nrecvcounts_c)
      sum_f = sum(collcom%nrecvcounts_f)
    
      !$omp parallel default(private) &
      !$omp shared(collcom, psit_c,psit_f, psitwork_c, psitwork_f,sum_c,sum_f,m)
    
      m = mod(sum_c,7)
    
      if(m/=0) then
        do i=1,m
           ind = collcom%iexpand_c(i)
           psitwork_c(ind) = psit_c(i)
        end do
      end if
    
      ! coarse part
    
      !$omp do
      do i=m+1,sum_c,7
          psitwork_c(collcom%iexpand_c(i+0))=psit_c(i+0)
          psitwork_c(collcom%iexpand_c(i+1))=psit_c(i+1)
          psitwork_c(collcom%iexpand_c(i+2))=psit_c(i+2)
          psitwork_c(collcom%iexpand_c(i+3))=psit_c(i+3)
          psitwork_c(collcom%iexpand_c(i+4))=psit_c(i+4)
          psitwork_c(collcom%iexpand_c(i+5))=psit_c(i+5)
          psitwork_c(collcom%iexpand_c(i+6))=psit_c(i+6)
      end do
      !$omp end do
    
      ! fine part
    
      !$omp do
      do i=1,sum_f
          ind=collcom%iexpand_f(i)
          psitwork_f(7*ind-6)=psit_f(7*i-6)
          psitwork_f(7*ind-5)=psit_f(7*i-5)
          psitwork_f(7*ind-4)=psit_f(7*i-4)
          psitwork_f(7*ind-3)=psit_f(7*i-3)
          psitwork_f(7*ind-2)=psit_f(7*i-2)
          psitwork_f(7*ind-1)=psit_f(7*i-1)
          psitwork_f(7*ind-0)=psit_f(7*i-0)
      end do
      !$omp end do
      !$omp end parallel
    
    end subroutine transpose_switch_psit



    subroutine transpose_communicate_psit(iproc, nproc, collcom, psitwork_c, psitwork_f, psiwork_c, psiwork_f)
      use module_base
      use module_types
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(collective_comms),intent(in) :: collcom
      real(kind=8),dimension(collcom%ndimind_c),intent(in) :: psitwork_c
      real(kind=8),dimension(7*collcom%ndimind_f),intent(in) :: psitwork_f
      real(kind=8),dimension(collcom%ndimpsi_c),intent(out) :: psiwork_c
      real(kind=8),dimension(7*collcom%ndimpsi_f),intent(out) :: psiwork_f
      
      ! Local variables
      integer :: ierr
      !!integer :: iall, ist, ist_c, ist_f, jproc, iisend, iirecv, istat
      !!real(kind=8),dimension(:),allocatable :: psiwork, psitwork
      !!integer,dimension(:),allocatable :: nsendcounts, nsenddspls, nrecvcounts, nrecvdspls
      !!character(len=*),parameter :: subname='transpose_communicate_psit'
    
      !call mpi_comm_size(bigdft_mpi%mpi_comm, nproc, ierr)
      !call mpi_comm_rank(bigdft_mpi%mpi_comm, iproc, ierr)
    
      !!allocate(psiwork(collcom%ndimpsi_c+7*collcom%ndimpsi_f), stat=istat)
      !!call memocc(istat, psiwork, 'psiwork', subname)
      !!allocate(psitwork(sum(collcom%nrecvcounts_c)+7*sum(collcom%nrecvcounts_f)), stat=istat)
      !!call memocc(istat, psitwork, 'psitwork', subname)
      !!allocate(nsendcounts(0:nproc-1), stat=istat)
      !!call memocc(istat, nsendcounts, 'nsendcounts', subname)
      !!allocate(nsenddspls(0:nproc-1), stat=istat)
      !!call memocc(istat, nsenddspls, 'nsenddspls', subname)
      !!allocate(nrecvcounts(0:nproc-1), stat=istat)
      !!call memocc(istat, nrecvcounts, 'nrecvcounts', subname)
      !!allocate(nrecvdspls(0:nproc-1), stat=istat)
      !!call memocc(istat, nrecvdspls, 'nrecvdspls', subname)
    
      !!ist=1
      !!ist_c=1
      !!ist_f=1
      !!iisend=0
      !!iirecv=0
      !!do jproc=0,nproc-1
      !!    if(collcom%nrecvcounts_c(jproc)>0) call vcopy(collcom%nrecvcounts_c(jproc), psitwork_c(ist_c), 1, psitwork(ist), 1)
      !!    ist_c=ist_c+collcom%nrecvcounts_c(jproc)
      !!    ist=ist+collcom%nrecvcounts_c(jproc)
      !!    if(collcom%nrecvcounts_f(jproc)>0) call vcopy(7*collcom%nrecvcounts_f(jproc), psitwork_f(ist_f), 1, psitwork(ist), 1)
      !!    ist_f=ist_f+7*collcom%nrecvcounts_f(jproc)
      !!    ist=ist+7*collcom%nrecvcounts_f(jproc)
      !!    nsendcounts(jproc)=collcom%nsendcounts_c(jproc)+7*collcom%nsendcounts_f(jproc)
      !!    nsenddspls(jproc)=iisend
      !!    nrecvcounts(jproc)=collcom%nrecvcounts_c(jproc)+7*collcom%nrecvcounts_f(jproc)
      !!    nrecvdspls(jproc)=iirecv
      !!    iisend=iisend+nsendcounts(jproc)
      !!    iirecv=iirecv+nrecvcounts(jproc)
      !!end do
    
    
      ! coarse part
       call mpi_alltoallv(psitwork_c, collcom%nrecvcounts_c, collcom%nrecvdspls_c, mpi_double_precision, psiwork_c, &
            collcom%nsendcounts_c, collcom%nsenddspls_c, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
    
      ! fine part
       call mpi_alltoallv(psitwork_f, 7*collcom%nrecvcounts_f, 7*collcom%nrecvdspls_f, mpi_double_precision, psiwork_f, &
            7*collcom%nsendcounts_f, 7*collcom%nsenddspls_f, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
      !!call mpi_alltoallv(psitwork, nrecvcounts, nrecvdspls, mpi_double_precision, psiwork, &
      !!     nsendcounts, nsenddspls, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
    
      !!ist=1
      !!ist_c=1
      !!ist_f=1
      !!do jproc=0,nproc-1
      !!    if(collcom%nsendcounts_c(jproc)>0) call vcopy(collcom%nsendcounts_c(jproc), psiwork(ist), 1, psiwork_c(ist_c), 1)
      !!    ist_c=ist_c+collcom%nsendcounts_c(jproc)
      !!    ist=ist+collcom%nsendcounts_c(jproc)
      !!    if(collcom%nsendcounts_f(jproc)>0) call vcopy(7*collcom%nsendcounts_f(jproc), psiwork(ist), 1, psiwork_f(ist_f), 1)
      !!    ist_f=ist_f+7*collcom%nsendcounts_f(jproc)
      !!    ist=ist+7*collcom%nsendcounts_f(jproc)
      !!end do
    
      !!iall=-product(shape(psiwork))*kind(psiwork)
      !!deallocate(psiwork, stat=istat)
      !!call memocc(istat, iall, 'psiwork', subname)
      !!iall=-product(shape(psitwork))*kind(psitwork)
      !!deallocate(psitwork, stat=istat)
      !!call memocc(istat, iall, 'psitwork', subname)
      !!iall=-product(shape(nsendcounts))*kind(nsendcounts)
      !!deallocate(nsendcounts, stat=istat)
      !!call memocc(istat, iall, 'nsendcounts', subname)
      !!iall=-product(shape(nsenddspls))*kind(nsenddspls)
      !!deallocate(nsenddspls, stat=istat)
      !!call memocc(istat, iall, 'nsenddspls', subname)
      !!iall=-product(shape(nrecvcounts))*kind(nrecvcounts)
      !!deallocate(nrecvcounts, stat=istat)
      !!call memocc(istat, iall, 'nrecvcounts', subname)
      !!iall=-product(shape(nrecvdspls))*kind(nrecvdspls)
      !!deallocate(nrecvdspls, stat=istat)
      !!call memocc(istat, iall, 'nrecvdspls', subname)
    
    end subroutine transpose_communicate_psit



    subroutine transpose_unswitch_psi(npsidim_orbs, orbs, collcom, psiwork_c, psiwork_f, psi, lzd)
      use module_base
      use module_types
      implicit none
      
      ! Caling arguments
      integer, intent(in) :: npsidim_orbs
      type(orbitals_data),intent(in) :: orbs
      type(collective_comms),intent(in) :: collcom
      real(kind=8),dimension(collcom%ndimpsi_c),intent(in) :: psiwork_c
      real(kind=8),dimension(7*collcom%ndimpsi_f),intent(in) :: psiwork_f
      real(kind=8),dimension(npsidim_orbs),intent(out) :: psi
      type(local_zone_descriptors),intent(in),optional :: lzd
      
      ! Local variables
      integer :: i, ind, iorb, iiorb, ilr, i_tot, i_c, i_f, istat, iall,m
      real(kind=8),dimension(:),allocatable :: psi_c, psi_f
      character(len=*),parameter :: subname='transpose_unswitch_psi'
      
      
      allocate(psi_c(collcom%ndimpsi_c), stat=istat)
      call memocc(istat, psi_c, 'psi_c', subname)
      allocate(psi_f(7*collcom%ndimpsi_f), stat=istat)
      call memocc(istat, psi_f, 'psi_f', subname)
      
      !$omp parallel default(private) &
      !$omp shared(collcom, psiwork_c, psi_c,psi_f,psiwork_f,m)
    
      m = mod(collcom%ndimpsi_c,7)
    
      if(m/=0) then
        do i = 1,m
         ind=collcom%irecvbuf_c(i)
         psi_c(ind)=psiwork_c(i) 
        end do
      end if
    
      ! coarse part
    
      !$omp do
        do i=m+1,collcom%ndimpsi_c,7
            psi_c(collcom%irecvbuf_c(i+0))=psiwork_c(i+0)
            psi_c(collcom%irecvbuf_c(i+1))=psiwork_c(i+1)
            psi_c(collcom%irecvbuf_c(i+2))=psiwork_c(i+2)
            psi_c(collcom%irecvbuf_c(i+3))=psiwork_c(i+3)
            psi_c(collcom%irecvbuf_c(i+4))=psiwork_c(i+4)
            psi_c(collcom%irecvbuf_c(i+5))=psiwork_c(i+5)
            psi_c(collcom%irecvbuf_c(i+6))=psiwork_c(i+6)
        end do
      !$omp end do
      
      ! fine part
     
      !$omp do
       do i=1,collcom%ndimpsi_f
            ind=collcom%irecvbuf_f(i)
            psi_f(7*ind-6)=psiwork_f(7*i-6)
            psi_f(7*ind-5)=psiwork_f(7*i-5)
            psi_f(7*ind-4)=psiwork_f(7*i-4)
            psi_f(7*ind-3)=psiwork_f(7*i-3)
            psi_f(7*ind-2)=psiwork_f(7*i-2)
            psi_f(7*ind-1)=psiwork_f(7*i-1)
            psi_f(7*ind-0)=psiwork_f(7*i-0)
        end do
      !$omp end do
      !$omp end parallel
    
        if(present(lzd)) then
            ! glue together coarse and fine part
    
            i_tot=0
            i_c=0
            i_f=0
            do iorb=1,orbs%norbp
                iiorb=orbs%isorb+iorb
                ilr=orbs%inwhichlocreg(iiorb)
    
                call vcopy(lzd%llr(ilr)%wfd%nvctr_c,psi_c(i_c+1),1,psi(i_tot+1),1)
    
                i_c = i_c + lzd%llr(ilr)%wfd%nvctr_c
                i_tot = i_tot + lzd%llr(ilr)%wfd%nvctr_c
                
                call vcopy(7*lzd%llr(ilr)%wfd%nvctr_f,psi_f(i_f+1),1,psi(i_tot+1),1)
    
    
                i_f = i_f + 7*lzd%llr(ilr)%wfd%nvctr_f
                i_tot = i_tot + 7*lzd%llr(ilr)%wfd%nvctr_f
    
    
            end do
        !!$omp end parallel 
    
        else
            call vcopy(collcom%ndimpsi_c, psi_c(1), 1, psi(1), 1)
        end if
      
      iall=-product(shape(psi_c))*kind(psi_c)
      deallocate(psi_c, stat=istat)
      call memocc(istat, iall, 'psi_c', subname)
      iall=-product(shape(psi_f))*kind(psi_f)
      deallocate(psi_f, stat=istat)
      call memocc(istat, iall, 'psi_f', subname)
    
    end subroutine transpose_unswitch_psi



    subroutine transpose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, psi, psit_c, psit_f, lzd)
      use module_base
      use module_types
      !use module_interfaces, except_this_one => transpose_localized
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, npsidim_orbs
      type(orbitals_data),intent(in) :: orbs
      type(collective_comms),intent(in) :: collcom
      real(kind=8),dimension(npsidim_orbs),intent(in) :: psi
      real(kind=8),dimension(collcom%ndimind_c),intent(out) :: psit_c
      real(kind=8),dimension(7*collcom%ndimind_f),intent(out) :: psit_f
      type(local_zone_descriptors),optional,intent(in) :: lzd
      
      ! Local variables
      real(kind=8),dimension(:),allocatable :: psiwork_c, psiwork_f, psitwork_c, psitwork_f
      integer :: istat, iall
      character(len=*),parameter :: subname='transpose_localized'
      
      allocate(psiwork_c(collcom%ndimpsi_c), stat=istat)
      call memocc(istat, psiwork_c, 'psiwork_c', subname)
      allocate(psiwork_f(7*collcom%ndimpsi_f), stat=istat)
      call memocc(istat, psiwork_f, 'psiwork_f', subname)
      allocate(psitwork_c(sum(collcom%nrecvcounts_c)), stat=istat)
      call memocc(istat, psitwork_c, 'psitwork_c', subname)
      allocate(psitwork_f(7*sum(collcom%nrecvcounts_f)), stat=istat)
      call memocc(istat, psitwork_f, 'psitwork_f', subname)
      
      call timing(iproc,'Un-TransSwitch','ON')
      if(present(lzd)) then
          call transpose_switch_psi(npsidim_orbs, orbs, collcom, psi, psiwork_c, psiwork_f, lzd)
      else
          call transpose_switch_psi(npsidim_orbs, orbs, collcom, psi, psiwork_c, psiwork_f)
      end if
      call timing(iproc,'Un-TransSwitch','OF')
    
      call timing(iproc,'Un-TransComm  ','ON')
      if(nproc>1) then
          call transpose_communicate_psi(iproc, nproc, collcom, psiwork_c, psiwork_f, psitwork_c, psitwork_f)
      else
          psitwork_c=psiwork_c
          psitwork_f=psiwork_f
      end if
      call timing(iproc,'Un-TransComm  ','OF')
    
      call timing(iproc,'Un-TransSwitch','ON')
      call transpose_unswitch_psit(collcom, psitwork_c, psitwork_f, psit_c, psit_f)
      call timing(iproc,'Un-TransSwitch','OF')
      
      iall=-product(shape(psiwork_c))*kind(psiwork_c)
      deallocate(psiwork_c, stat=istat)
      call memocc(istat, iall, 'psiwork_c', subname)
      iall=-product(shape(psiwork_f))*kind(psiwork_f)
      deallocate(psiwork_f, stat=istat)
      call memocc(istat, iall, 'psiwork_f', subname)
      iall=-product(shape(psitwork_c))*kind(psitwork_c)
      deallocate(psitwork_c, stat=istat)
      call memocc(istat, iall, 'psitwork_c', subname)
      iall=-product(shape(psitwork_f))*kind(psitwork_f)
      deallocate(psitwork_f, stat=istat)
      call memocc(istat, iall, 'psitwork_f', subname)
      
    end subroutine transpose_localized



    subroutine untranspose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, psit_c, psit_f, psi, lzd)
      use module_base
      use module_types
      !use module_interfaces, except_this_one => untranspose_localized
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, npsidim_orbs
      type(orbitals_data),intent(in) :: orbs
      type(collective_comms),intent(in) :: collcom
      real(kind=8),dimension(collcom%ndimind_c),intent(in) :: psit_c
      real(kind=8),dimension(7*collcom%ndimind_f),intent(in) :: psit_f
      real(kind=8),dimension(npsidim_orbs),intent(out) :: psi
      type(local_zone_descriptors),optional,intent(in) :: lzd
      
      ! Local variables
      real(kind=8),dimension(:),allocatable :: psiwork_c, psiwork_f, psitwork_c, psitwork_f
      integer :: istat, iall
      character(len=*),parameter :: subname='untranspose_localized'
      
      allocate(psiwork_c(collcom%ndimpsi_c), stat=istat)
      call memocc(istat, psiwork_c, 'psiwork_c', subname)
      allocate(psiwork_f(7*collcom%ndimpsi_f), stat=istat)
      call memocc(istat, psiwork_f, 'psiwork_f', subname)
      allocate(psitwork_c(sum(collcom%nrecvcounts_c)), stat=istat)
      call memocc(istat, psitwork_c, 'psitwork_c', subname)
      allocate(psitwork_f(7*sum(collcom%nrecvcounts_f)), stat=istat)
      call memocc(istat, psitwork_f, 'psitwork_f', subname)
    
      call timing(iproc,'Un-TransSwitch','ON')
      call transpose_switch_psit(collcom, psit_c, psit_f, psitwork_c, psitwork_f)
      call timing(iproc,'Un-TransSwitch','OF')
    
      call timing(iproc,'Un-TransComm  ','ON')
      if(nproc>1) then
          call transpose_communicate_psit(iproc, nproc, collcom, psitwork_c, psitwork_f, psiwork_c, psiwork_f)
      else
          psiwork_c=psitwork_c
          psiwork_f=psitwork_f
      end if
      call timing(iproc,'Un-TransComm  ','OF')
    
      call timing(iproc,'Un-TransSwitch','ON')
      if(present(lzd)) then
          call transpose_unswitch_psi(npsidim_orbs, orbs, collcom, psiwork_c, psiwork_f, psi, lzd)
      else
          call transpose_unswitch_psi(npsidim_orbs, orbs, collcom, psiwork_c, psiwork_f, psi)
      end if
      call timing(iproc,'Un-TransSwitch','OF')
      
      iall=-product(shape(psiwork_c))*kind(psiwork_c)
      deallocate(psiwork_c, stat=istat)
      call memocc(istat, iall, 'psiwork_c', subname)
      iall=-product(shape(psiwork_f))*kind(psiwork_f)
      deallocate(psiwork_f, stat=istat)
      call memocc(istat, iall, 'psiwork_f', subname)
      iall=-product(shape(psitwork_c))*kind(psitwork_c)
      deallocate(psitwork_c, stat=istat)
      call memocc(istat, iall, 'psitwork_c', subname)
      iall=-product(shape(psitwork_f))*kind(psitwork_f)
      deallocate(psitwork_f, stat=istat)
      call memocc(istat, iall, 'psitwork_f', subname)
      
    end subroutine untranspose_localized





    ! The sumrho routines

    subroutine init_collective_comms_sumrho(iproc, nproc, lzd, orbs, nscatterarr, collcom_sr)
      use module_base
      use module_types
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(local_zone_descriptors),intent(in) :: lzd
      type(orbitals_data),intent(in) :: orbs
      integer,dimension(0:nproc-1,4),intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
      type(collective_comms),intent(inout) :: collcom_sr
    
      ! Local variables
      integer :: ierr, istat, iall, ipt, ii
      real(kind=8) :: weight_tot, weight_ideal
      integer,dimension(:,:),allocatable :: istartend
      character(len=*),parameter :: subname='init_collective_comms_sumrho'
      real(kind=8),dimension(:),allocatable :: weights_per_slice, weights_per_zpoint
    
      ! Note: all weights are double precision to avoid integer overflow
      call timing(iproc,'init_collco_sr','ON')
    
      allocate(istartend(2,0:nproc-1), stat=istat)
      call memocc(istat, istartend, 'istartend', subname)
    
      allocate(weights_per_slice(0:nproc-1), stat=istat)
      call memocc(istat, weights_per_slice, 'weights_per_slice', subname)
    
      allocate(weights_per_zpoint(lzd%glr%d%n3i), stat=istat)
      call memocc(istat, weights_per_zpoint, 'weights_per_zpoint', subname)
    
      call get_weights_sumrho(iproc, nproc, orbs, lzd, nscatterarr, weight_tot, weight_ideal, &
           weights_per_slice, weights_per_zpoint)
    
      call assign_weight_to_process_sumrho(iproc, nproc, weight_tot, weight_ideal, weights_per_slice, &
           lzd, orbs, nscatterarr, istartend, collcom_sr%nptsp_c)
    
      iall = -product(shape(weights_per_slice))*kind(weights_per_slice)
      deallocate(weights_per_slice,stat=istat)
      call memocc(istat, iall, 'weights_per_slice', subname)
    

      call allocate_MPI_communication_arrays(nproc, collcom_sr, only_coarse=.true.)

      call determine_communication_arrays_sumrho(iproc, nproc, collcom_sr%nptsp_c, lzd, orbs, istartend, &
           collcom_sr%nsendcounts_c, collcom_sr%nsenddspls_c, &
           collcom_sr%nrecvcounts_c, collcom_sr%nrecvdspls_c, collcom_sr%ndimpsi_c)

      !Now set some integers in the collcomm structure
      collcom_sr%ndimind_c = sum(collcom_sr%nrecvcounts_c)

      call allocate_local_communications_arrays(collcom_sr, only_coarse=.true.)
    
      call determine_num_orbs_per_gridpoint_sumrho(iproc, nproc, collcom_sr%nptsp_c, lzd, orbs, &
           istartend, weight_tot, weights_per_zpoint, collcom_sr%norb_per_gridpoint_c)
    
      ! Some check
      ii=sum(collcom_sr%norb_per_gridpoint_c)
      if (ii/=collcom_sr%ndimind_c) stop 'ii/=collcom_sr%ndimind_c'
    
    
      collcom_sr%psit_c=f_malloc_ptr(collcom_sr%ndimind_c,id='collcom_sr%psit_c')
    
      call get_switch_indices_sumrho(iproc, nproc, collcom_sr%nptsp_c, collcom_sr%ndimpsi_c, collcom_sr%ndimind_c, lzd, &
           orbs, istartend, collcom_sr%norb_per_gridpoint_c, collcom_sr%nsendcounts_c, collcom_sr%nsenddspls_c, &
           collcom_sr%nrecvcounts_c, collcom_sr%nrecvdspls_c, collcom_sr%isendbuf_c, collcom_sr%irecvbuf_c, &
           collcom_sr%iextract_c, collcom_sr%iexpand_c, collcom_sr%indexrecvorbital_c)
    
      ! These variables are used in various subroutines to speed up the code
      collcom_sr%isptsp_c(1) = 0
      do ipt=2,collcom_sr%nptsp_c
            collcom_sr%isptsp_c(ipt) = collcom_sr%isptsp_c(ipt-1) + collcom_sr%norb_per_gridpoint_c(ipt-1)
      end do
    
      call allocate_MPI_communications_arrays_repartition(nproc, collcom_sr)
    
      call communication_arrays_repartitionrho(iproc, nproc, lzd, nscatterarr, istartend, &
           collcom_sr%nsendcounts_repartitionrho, collcom_sr%nsenddspls_repartitionrho, &
           collcom_sr%nrecvcounts_repartitionrho, collcom_sr%nrecvdspls_repartitionrho)
    
      call communication_arrays_repartitionrho_general(iproc, nproc, lzd, nscatterarr, istartend, & 
           collcom_sr%ncomms_repartitionrho, collcom_sr%commarr_repartitionrho)
    
      iall = -product(shape(weights_per_zpoint))*kind(weights_per_zpoint)
      deallocate(weights_per_zpoint,stat=istat)
      call memocc(istat, iall, 'weights_per_zpoint', subname)
    
      iall = -product(shape(istartend))*kind(istartend)
      deallocate(istartend,stat=istat)
      call memocc(istat, iall, 'istartend', subname)
    
      call timing(iproc,'init_collco_sr','OF')
    
    end subroutine init_collective_comms_sumrho



    subroutine get_weights_sumrho(iproc, nproc, orbs, lzd, nscatterarr, &
               weight_tot, weight_ideal, weights_per_slice, weights_per_zpoint)
      use module_base
      use module_types
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(orbitals_data),intent(in) :: orbs
      type(local_zone_descriptors),intent(in) :: lzd
      integer,dimension(0:nproc-1,4),intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
      real(kind=8),intent(out) :: weight_tot, weight_ideal
      real(kind=8),dimension(0:nproc-1),intent(out) :: weights_per_slice
      real(kind=8),dimension(lzd%glr%d%n3i),intent(out) :: weights_per_zpoint
    
      ! Local variables
      integer :: iorb, ilr, ierr, i3, i2, i1, is1, ie1, is2, ie2, is3, ie3, istat, iall
      real(kind=8) :: tt, zz
      real(kind=8),dimension(:,:),allocatable :: weight_xy
    
      call f_routine(id='get_weights_sumrho')
    
      call to_zero(lzd%glr%d%n3i, weights_per_zpoint(1))
    
      weight_xy=f_malloc((/lzd%glr%d%n1i,lzd%glr%d%n2i/),id='weight_xy')
    
      !write(*,*) 'iproc, nscatterarr', iproc, nscatterarr(iproc,:)
    
      tt=0.d0
      weights_per_slice(:) = 0.0d0
      do i3=nscatterarr(iproc,3)+1,nscatterarr(iproc,3)+nscatterarr(iproc,2)
          call to_zero(lzd%glr%d%n1i*lzd%glr%d%n2i, weight_xy(1,1))
          do iorb=1,orbs%norb
              ilr=orbs%inwhichlocreg(iorb)
              is3=1+lzd%Llr(ilr)%nsi3
              ie3=lzd%Llr(ilr)%nsi3+lzd%llr(ilr)%d%n3i
              if (is3>i3 .or. i3>ie3) cycle
              is1=1+lzd%Llr(ilr)%nsi1
              ie1=lzd%Llr(ilr)%nsi1+lzd%llr(ilr)%d%n1i
              is2=1+lzd%Llr(ilr)%nsi2
              ie2=lzd%Llr(ilr)%nsi2+lzd%llr(ilr)%d%n2i
              !$omp parallel default(none) shared(is2, ie2, is1, ie1, weight_xy) private(i2, i1)
              !$omp do
              do i2=is2,ie2
                  do i1=is1,ie1
                      weight_xy(i1,i2) = weight_xy(i1,i2)+1.d0
                  end do
              end do
              !$omp end do
              !$omp end parallel
          end do
          zz=0.d0
          !$omp parallel default(none) shared(lzd, weight_xy, zz, tt) private(i2, i1)
          !$omp do reduction(+: tt, zz)
          do i2=1,lzd%glr%d%n2i
              do i1=1,lzd%glr%d%n1i
                 tt = tt + .5d0*(weight_xy(i1,i2)*(weight_xy(i1,i2)+1.d0))
                 zz = zz + .5d0*(weight_xy(i1,i2)*(weight_xy(i1,i2)))
              end do
          end do
          !$omp end do
          !$omp end parallel
          weights_per_zpoint(i3)=zz
      end do
      weights_per_slice(iproc)=tt
      call mpiallred(weights_per_slice(0), nproc, mpi_sum, bigdft_mpi%mpi_comm, ierr)
      if (nproc>1) then
         call mpi_allreduce(tt, weight_tot, 1, mpi_double_precision, mpi_sum, bigdft_mpi%mpi_comm, ierr)
      else
         weight_tot=tt
      end if
      call mpiallred(weights_per_zpoint(1), lzd%glr%d%n3i, mpi_sum, bigdft_mpi%mpi_comm, ierr)
    
      call f_free(weight_xy)
    
      ! Ideal weight per process
      weight_ideal = weight_tot/dble(nproc)
    
      call f_release_routine()
    
    end subroutine get_weights_sumrho


    subroutine assign_weight_to_process_sumrho(iproc, nproc, weight_tot, weight_ideal, weights_per_slice, &
               lzd, orbs, nscatterarr, istartend, nptsp)
      use module_base
      use module_types
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      real(kind=8),intent(in) :: weight_tot, weight_ideal
      real(kind=8),dimension(0:nproc-1),intent(in) :: weights_per_slice
      type(local_zone_descriptors),intent(in) :: lzd
      type(orbitals_data),intent(in) :: orbs
      integer,dimension(0:nproc-1,4),intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
      integer,dimension(2,0:nproc-1),intent(out) :: istartend
      integer,intent(out) :: nptsp
    
      ! Local variables
      integer :: jproc, i1, i2, i3, ii, iorb, ilr, is1, ie1, is2, ie2, is3, ie3, ierr, istat, iall, jproc_out
      real(kind=8) :: tt, ttt
      real(kind=8),dimension(:,:),allocatable :: slicearr
      real(8),dimension(:,:),allocatable :: weights_startend
    
      call f_routine(id='assign_weight_to_process_sumrho')
    
      weights_startend=f_malloc((/1.to.2,0.to.nproc-1/),id='weights_startend')
    
      tt=0.d0
      weights_startend(1,0)=0.d0
      do jproc=0,nproc-2
          tt=tt+weight_ideal
          weights_startend(2,jproc)=dble(floor(tt,kind=8))
          weights_startend(1,jproc+1)=dble(floor(tt,kind=8))+1.d0
      end do
      weights_startend(2,nproc-1)=weight_tot
    
    
      ! Iterate through all grid points and assign them to processes such that the
      ! load balancing is optimal.
      if (nproc==1) then
          istartend(1,0)=1
          istartend(2,0)=lzd%glr%d%n1i*lzd%glr%d%n2i*lzd%glr%d%n3i
      else
          slicearr=f_malloc((/lzd%glr%d%n1i,lzd%glr%d%n2i/),id='slicearr')
          istartend(1,:)=0
          istartend(2,:)=0
          tt=0.d0
          jproc=0
          ii=0
          outer_loop: do jproc_out=0,nproc-1
              if (tt+weights_per_slice(jproc_out)<weights_startend(1,iproc)) then
                  tt=tt+weights_per_slice(jproc_out)
                  ii=ii+nscatterarr(jproc_out,2)*lzd%glr%d%n1i*lzd%glr%d%n2i
                  cycle outer_loop
              end if
              i3_loop: do i3=nscatterarr(jproc_out,3)+1,nscatterarr(jproc_out,3)+nscatterarr(jproc_out,2)
                  call to_zero(lzd%glr%d%n1i*lzd%glr%d%n2i, slicearr(1,1))
                  do iorb=1,orbs%norb
                      ilr=orbs%inwhichlocreg(iorb)
                      is1=1+lzd%Llr(ilr)%nsi1
                      ie1=lzd%Llr(ilr)%nsi1+lzd%llr(ilr)%d%n1i
                      is2=1+lzd%Llr(ilr)%nsi2
                      ie2=lzd%Llr(ilr)%nsi2+lzd%llr(ilr)%d%n2i
                      is3=1+lzd%Llr(ilr)%nsi3
                      ie3=lzd%Llr(ilr)%nsi3+lzd%llr(ilr)%d%n3i
                      if (is3>i3 .or. i3>ie3) cycle
                      !$omp parallel default(none) shared(lzd, slicearr, is1, ie1, is2, ie2) private(i1, i2)
                      !$omp do
                      do i2=1,lzd%glr%d%n2i
                          do i1=1,lzd%glr%d%n1i
                              if (is1<=i1 .and. i1<=ie1 .and. is2<=i2 .and. i2<=ie2) then
                                  slicearr(i1,i2)=slicearr(i1,i2)+1.d0
                              end if
                          end do
                      end do
                      !$omp end do
                      !$omp end parallel
                  end do
                  do i2=1,lzd%glr%d%n2i
                      do i1=1,lzd%glr%d%n1i
                          ii=ii+1
                          tt=tt+.5d0*slicearr(i1,i2)*(slicearr(i1,i2)+1.d0)
                          if (tt>=weights_startend(1,iproc)) then
                              istartend(1,iproc)=ii
                              exit outer_loop
                          end if
                      end do
                   end do
               end do i3_loop
            end do outer_loop
            call f_free(slicearr)
      end if
    
      call mpiallred(istartend(1,0), 2*nproc, mpi_sum, bigdft_mpi%mpi_comm, ierr)
    
      do jproc=0,nproc-2
          istartend(2,jproc)=istartend(1,jproc+1)-1
      end do
      istartend(2,nproc-1)=lzd%glr%d%n1i*lzd%glr%d%n2i*lzd%glr%d%n3i
    
      do jproc=0,nproc-1
          if (iproc==jproc) then
              nptsp=istartend(2,jproc)-istartend(1,jproc)+1
          end if
      end do
    
      call f_free(weights_startend)
    
    
      ! Some check
      ii=nptsp
      call mpiallred(ii, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
      if (ii/=lzd%glr%d%n1i*lzd%glr%d%n2i*lzd%glr%d%n3i) then
          stop 'ii/=lzd%glr%d%n1i*lzd%glr%d%n2i*lzd%glr%d%n3i'
      end if
    
      call f_release_routine()
    
    
    end subroutine assign_weight_to_process_sumrho



    subroutine determine_num_orbs_per_gridpoint_sumrho(iproc, nproc, nptsp, lzd, orbs, &
               istartend, weight_tot, weights_per_zpoint, norb_per_gridpoint)
      use module_base
      use module_types
      use yaml_output
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, nptsp
      type(local_zone_descriptors),intent(in) :: lzd
      type(orbitals_data),intent(in) :: orbs
      integer,dimension(2,0:nproc-1),intent(in) :: istartend
      real(kind=8),intent(in) :: weight_tot
      real(kind=8),dimension(lzd%glr%d%n3i),intent(in) :: weights_per_zpoint
      integer,dimension(nptsp),intent(out) :: norb_per_gridpoint
    
      ! Local variables
      integer :: i3, ii, i2, i1, ipt, ilr, is1, ie1, is2, ie2, is3, ie3, iorb, ierr, i
      real(8) :: tt, weight_check
    
    
      if (nptsp>0) then
          call to_zero(nptsp, norb_per_gridpoint(1))
      end if
      do i3=1,lzd%glr%d%n3i
          if (i3*lzd%glr%d%n1i*lzd%glr%d%n2i<istartend(1,iproc) .or. &
              (i3-1)*lzd%glr%d%n1i*lzd%glr%d%n2i+1>istartend(2,iproc)) then
              cycle
          end if
          if (weights_per_zpoint(i3)==0.d0) then
              cycle
          end if
          do iorb=1,orbs%norb
              ilr=orbs%inwhichlocreg(iorb)
              is3=1+lzd%Llr(ilr)%nsi3
              ie3=lzd%Llr(ilr)%nsi3+lzd%llr(ilr)%d%n3i
              if (is3>i3 .or. i3>ie3) cycle
              is2=1+lzd%Llr(ilr)%nsi2
              ie2=lzd%Llr(ilr)%nsi2+lzd%llr(ilr)%d%n2i
              is1=1+lzd%Llr(ilr)%nsi1
              ie1=lzd%Llr(ilr)%nsi1+lzd%llr(ilr)%d%n1i
              !$omp parallel default(none) &
              !$omp shared(i3, is2, ie2, is1, ie1, lzd, istartend, iproc, norb_per_gridpoint) private(i2, i1, ii, ipt)
              !$omp do
              do i2=is2,ie2
                  do i1=is1,ie1
                      ii=(i3-1)*lzd%glr%d%n1i*lzd%glr%d%n2i+(i2-1)*lzd%glr%d%n1i+i1
                      if (ii>=istartend(1,iproc) .and. ii<=istartend(2,iproc)) then
                          ipt=ii-istartend(1,iproc)+1
                          norb_per_gridpoint(ipt)=norb_per_gridpoint(ipt)+1
                      end if
                  end do
              end do
              !$omp end do
              !$omp end parallel
          end do
      end do
    
      tt=0.d0
      !$omp parallel default(none) shared(tt, nptsp, norb_per_gridpoint) private(i)
      !$omp do reduction(+:tt)
      do i=1,nptsp
          tt=tt+.5d0*dble(norb_per_gridpoint(i)*(norb_per_gridpoint(i)+1))
      end do
      !$omp end do
      !$omp end parallel
      weight_check=tt
    
    
      ! Some check
      call mpiallred(weight_check, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
      if (abs(weight_check-weight_tot) > 1.d-3) then
          stop '2: weight_check/=weight_tot'
      else if (abs(weight_check-weight_tot) > 0.d0) then
         call yaml_warning('The total weight for density seems inconsistent! Ref:'//&
               trim(yaml_toa(weight_tot,fmt='(1pe25.17)'))//', Check:'//&
               trim(yaml_toa(weight_check,fmt='(1pe25.17)')))
      end if
    
    end subroutine determine_num_orbs_per_gridpoint_sumrho


    subroutine determine_communication_arrays_sumrho(iproc, nproc, nptsp, lzd, orbs, &
               istartend, nsendcounts, nsenddspls, nrecvcounts, &
               nrecvdspls, ndimpsi)
      use module_base
      use module_types
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, nptsp
      type(local_zone_descriptors),intent(in) :: lzd
      type(orbitals_data),intent(in) :: orbs
      integer,dimension(2,0:nproc-1),intent(in) :: istartend
      integer,dimension(0:nproc-1),intent(out) :: nsendcounts, nsenddspls, nrecvcounts, nrecvdspls
      integer,intent(out) :: ndimpsi
    
      ! Local variables
      integer :: iorb, iiorb, ilr, is1, ie1, is2, ie2, is3, ie3, jproc, i3, i2, i1, ind, ii, istat, iall, ierr, ii0
      integer,dimension(:),allocatable :: nsendcounts_tmp, nsenddspls_tmp, nrecvcounts_tmp, nrecvdspls_tmp
      character(len=*),parameter :: subname='determine_communication_arrays_sumrho'
    
    
      call to_zero(nproc,nsendcounts(0))
    
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          ilr=orbs%inwhichlocreg(iiorb)
          is1=1+lzd%Llr(ilr)%nsi1
          ie1=lzd%Llr(ilr)%nsi1+lzd%llr(ilr)%d%n1i
          is2=1+lzd%Llr(ilr)%nsi2
          ie2=lzd%Llr(ilr)%nsi2+lzd%llr(ilr)%d%n2i
          is3=1+lzd%Llr(ilr)%nsi3
          ie3=lzd%Llr(ilr)%nsi3+lzd%llr(ilr)%d%n3i
          do jproc=0,nproc-1
              ii=0
              do i3=is3,ie3
                  if (i3*lzd%glr%d%n1i*lzd%glr%d%n2i<istartend(1,jproc) .or. &
                      (i3-1)*lzd%glr%d%n1i*lzd%glr%d%n2i+1>istartend(2,jproc)) then
                      cycle
                  end if
                  ii0=0
                  !$omp parallel default(none) &
                  !$omp shared(i3, is2, ie2, is1, ie1, lzd, istartend, jproc, ii0) private(i2, i1, ind)
                  !$omp do reduction(+:ii0)
                  do i2=is2,ie2
                      do i1=is1,ie1
                        ind = (i3-1)*lzd%glr%d%n1i*lzd%glr%d%n2i+(i2-1)*lzd%glr%d%n1i+i1
                        if (ind>=istartend(1,jproc) .and. ind<=istartend(2,jproc)) then
                            !nsendcounts(jproc)=nsendcounts(jproc)+1
                            ii0=ii0+1
                        end if
                      end do
                  end do
                  !$omp end do
                  !$omp end parallel
                  ii=ii+ii0
              end do
             nsendcounts(jproc)=nsendcounts(jproc)+ii
           end do
      end do
    
    
      ! Some check
      ii=0
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          ilr=orbs%inwhichlocreg(iiorb)
          ii = ii + lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i
      end do
      if (ii/=sum(nsendcounts)) then
          stop 'ii/=sum(nsendcounts)'
      end if
      ndimpsi=ii
    
    
      nsenddspls(0)=0
      do jproc=1,nproc-1
          nsenddspls(jproc)=nsenddspls(jproc-1)+nsendcounts(jproc-1)
      end do
    
      allocate(nsendcounts_tmp(0:nproc-1), stat=istat)
      call memocc(istat, nsendcounts_tmp, 'nsendcounts_tmp', subname)
      allocate(nsenddspls_tmp(0:nproc-1), stat=istat)
      call memocc(istat, nsenddspls_tmp, 'nsenddspls_tmp', subname)
      allocate(nrecvcounts_tmp(0:nproc-1), stat=istat)
      call memocc(istat, nrecvcounts_tmp, 'nrecvcounts_tmp', subname)
      allocate(nrecvdspls_tmp(0:nproc-1), stat=istat)
      call memocc(istat, nrecvdspls_tmp, 'nrecvdspls_tmp', subname)
      nsendcounts_tmp=1
      nrecvcounts_tmp=1
      do jproc=0,nproc-1
          nsenddspls_tmp(jproc)=jproc
          nrecvdspls_tmp(jproc)=jproc
      end do
      if(nproc>1) then
          call mpi_alltoallv(nsendcounts, nsendcounts_tmp, nsenddspls_tmp, mpi_integer, nrecvcounts, &
               nrecvcounts_tmp, nrecvdspls_tmp, mpi_integer, bigdft_mpi%mpi_comm, ierr)
      else
          nrecvcounts=nsendcounts
      end if
      iall=-product(shape(nsendcounts_tmp))*kind(nsendcounts_tmp)
      deallocate(nsendcounts_tmp, stat=istat)
      call memocc(istat, iall, 'nsendcounts_tmp', subname)
      iall=-product(shape(nsenddspls_tmp))*kind(nsenddspls_tmp)
      deallocate(nsenddspls_tmp, stat=istat)
      call memocc(istat, iall, 'nsenddspls_tmp', subname)
      iall=-product(shape(nrecvcounts_tmp))*kind(nrecvcounts_tmp)
      deallocate(nrecvcounts_tmp, stat=istat)
      call memocc(istat, iall, 'nrecvcounts_tmp', subname)
      iall=-product(shape(nrecvdspls_tmp))*kind(nrecvdspls_tmp)
      deallocate(nrecvdspls_tmp, stat=istat)
      call memocc(istat, iall, 'nrecvdspls_tmp', subname)
    
      !!ndimind = sum(nrecvcounts)
    
      !!! Some check
      !!ii=sum(norb_per_gridpoint)
      !!if (ii/=ndimind) stop 'ii/=sum(nrecvcounts)'
    
      nrecvdspls(0)=0
      do jproc=1,nproc-1
          nrecvdspls(jproc)=nrecvdspls(jproc-1)+nrecvcounts(jproc-1)
      end do
    
    end subroutine determine_communication_arrays_sumrho



    subroutine get_switch_indices_sumrho(iproc, nproc, nptsp, ndimpsi, ndimind, lzd, orbs, istartend, &
               norb_per_gridpoint, nsendcounts, nsenddspls, nrecvcounts, nrecvdspls, &
               isendbuf, irecvbuf, iextract, iexpand, indexrecvorbital)
      use module_base
      use module_types
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, nptsp, ndimpsi, ndimind
      type(local_zone_descriptors),intent(in) :: lzd
      type(orbitals_data),intent(in) :: orbs
      integer,dimension(2,0:nproc-1),intent(in) :: istartend
      integer,dimension(nptsp),intent(in) :: norb_per_gridpoint
      integer,dimension(0:nproc-1),intent(in) :: nsendcounts, nsenddspls, nrecvcounts, nrecvdspls
      integer,dimension(ndimpsi),intent(out) :: isendbuf, irecvbuf
      integer,dimension(ndimind),intent(out) :: iextract, iexpand, indexrecvorbital
    
      ! Local variables
      integer :: jproc, iitot, iiorb, ilr, is1, ie1, is2, ie2, is3, ie3, i3, i2, i1, ind, indglob, istat, iall, ierr, ii
      integer :: iorb, i, ipt
      integer,dimension(:),allocatable :: nsend, indexsendbuf, indexsendorbital, indexsendorbital2, indexrecvorbital2
      integer,dimension(:),allocatable :: gridpoint_start, indexrecvbuf
      character(len=*),parameter :: subname='get_switch_indices_sumrho'
    
    
      allocate(nsend(0:nproc-1), stat=istat)
      call memocc(istat, nsend, 'nsend', subname)
      nsend=0
      allocate(indexsendbuf(ndimpsi), stat=istat)
      call memocc(istat, indexsendbuf, 'indexsendbuf', subname)
      allocate(indexsendorbital(ndimpsi), stat=istat)
      call memocc(istat, indexsendorbital, 'indexsendorbital', subname)
      !!allocate(isendbuf(ndimpsi), stat=istat)
      !!call memocc(istat, isendbuf, 'isendbuf', subname)
    
      iitot=0
      !!$omp parallel default(shared) &
      !!$omp private(iorb, iiorb, ilr, is1, ie1, is2, ie2, is3, ie3, i3, i2, i1, indglob, ind)
      !!$omp do lastprivate(iitot)
      do jproc=0,nproc-1
          iitot=0
          do iorb=1,orbs%norbp
              iiorb=orbs%isorb+iorb
              ilr=orbs%inwhichlocreg(iiorb)
              is1=1+lzd%Llr(ilr)%nsi1
              ie1=lzd%Llr(ilr)%nsi1+lzd%llr(ilr)%d%n1i
              is2=1+lzd%Llr(ilr)%nsi2
              ie2=lzd%Llr(ilr)%nsi2+lzd%llr(ilr)%d%n2i
              is3=1+lzd%Llr(ilr)%nsi3
              ie3=lzd%Llr(ilr)%nsi3+lzd%llr(ilr)%d%n3i
              do i3=is3,ie3
                  if (i3*lzd%glr%d%n1i*lzd%glr%d%n2i<istartend(1,jproc) .or. &
                      (i3-1)*lzd%glr%d%n1i*lzd%glr%d%n2i+1>istartend(2,jproc)) then
                      iitot=iitot+(ie2-is2+1)*(ie1-is1+1)
                      cycle
                  end if
                  do i2=is2,ie2
                      do i1=is1,ie1
                          indglob = (i3-1)*lzd%glr%d%n1i*lzd%glr%d%n2i+(i2-1)*lzd%glr%d%n1i+i1
                          iitot=iitot+1
                          if (indglob>=istartend(1,jproc) .and. indglob<=istartend(2,jproc)) then
                              nsend(jproc)=nsend(jproc)+1
                              ind=nsenddspls(jproc)+nsend(jproc)
                              isendbuf(iitot)=ind
                              indexsendbuf(ind)=indglob
                              indexsendorbital(iitot)=iiorb
                              !exit
                          end if
                      end do
                  end do
              end do
          end do
      end do
      !!$omp end do
      !!$omp end parallel
    
    
      if(iitot/=ndimpsi) stop 'iitot/=ndimpsi'
    
      !check
      do jproc=0,nproc-1
          if(nsend(jproc)/=nsendcounts(jproc)) stop 'nsend(jproc)/=nsendcounts(jproc)'
      end do
    
    !!call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
    !!t2=mpi_wtime()
    !!tt=t2-t1
    !!if(iproc==0) write(*,*) 'time 5.1: iproc', iproc, tt
    
    
    
      !!allocate(irecvbuf(ndimpsi), stat=istat)
      !!call memocc(istat, irecvbuf, 'irecvbuf', subname)
    
      allocate(indexsendorbital2(ndimpsi), stat=istat)
      call memocc(istat, indexsendorbital2, 'indexsendorbital2', subname)
      indexsendorbital2=indexsendorbital
      do i=1,ndimpsi
          ind=isendbuf(i)
          indexsendorbital(ind)=indexsendorbital2(i)
      end do
    
      ! Inverse of isendbuf
      call get_reverse_indices(ndimpsi, isendbuf, irecvbuf)
    
      iall=-product(shape(indexsendorbital2))*kind(indexsendorbital2)
      deallocate(indexsendorbital2, stat=istat)
      call memocc(istat, iall, 'indexsendorbital2', subname)
    
    
      allocate(indexrecvbuf(ndimind), stat=istat)
      call memocc(istat, indexrecvbuf, 'indexrecvbuf', subname)
      !!allocate(indexrecvorbital(ndimind), stat=istat)
      !!call memocc(istat, indexrecvorbital, 'indexrecvorbital', subname)
    
      if(nproc>1) then
          ! Communicate indexsendbuf
          call mpi_alltoallv(indexsendbuf, nsendcounts, nsenddspls, mpi_integer, indexrecvbuf, &
               nrecvcounts, nrecvdspls, mpi_integer, bigdft_mpi%mpi_comm, ierr)
          ! Communicate indexsendorbitals
          call mpi_alltoallv(indexsendorbital, nsendcounts, nsenddspls, &
               mpi_integer, indexrecvorbital, &
               nrecvcounts, nrecvdspls, mpi_integer, bigdft_mpi%mpi_comm, ierr)
       else
           indexrecvbuf=indexsendbuf
           indexrecvorbital=indexsendorbital
       end if
    
      iall=-product(shape(indexsendbuf))*kind(indexsendbuf)
      deallocate(indexsendbuf, stat=istat)
      call memocc(istat, iall, 'indexsendbuf', subname)
    
      iall=-product(shape(indexsendorbital))*kind(indexsendorbital)
      deallocate(indexsendorbital, stat=istat)
      call memocc(istat, iall, 'indexsendorbital', subname)
    !!call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
    !!t2=mpi_wtime()
    !!tt=t2-t1
    !!if(iproc==0) write(*,*) 'time 5.2: iproc', iproc, tt
    
    
       allocate(gridpoint_start(istartend(1,iproc):istartend(2,iproc)), stat=istat)
       call memocc(istat, gridpoint_start, 'gridpoint_start', subname)
    
       ii=1
       do ipt=1,nptsp
           i=ipt+istartend(1,iproc)-1
           if (norb_per_gridpoint(ipt)>0) then
               gridpoint_start(i)=ii
           else
               gridpoint_start(i)=0
           end if
           ii=ii+norb_per_gridpoint(ipt)
       end do
    
       if (ii/=ndimind+1) stop '(ii/=ndimind+1)'
       if(maxval(gridpoint_start)>ndimind) stop '1: maxval(gridpoint_start)>sum(nrecvcountc)'
    
       !!allocate(iextract(ndimind), stat=istat)
       !!call memocc(istat, iextract, 'iextract', subname)
    
      ! Rearrange the communicated data
      do i=1,ndimind
          ii=indexrecvbuf(i)
          ind=gridpoint_start(ii)
          iextract(i)=ind
          gridpoint_start(ii)=gridpoint_start(ii)+1
      end do
    
      if(maxval(iextract)>ndimind) stop 'maxval(iextract)>ndimind'
      if(minval(iextract)<1) stop 'minval(iextract)<1'
    
      iall=-product(shape(indexrecvbuf))*kind(indexrecvbuf)
      deallocate(indexrecvbuf, stat=istat)
      call memocc(istat, iall, 'indexrecvbuf', subname)
    
    
      !! allocate(iexpand(ndimind), stat=istat)
      !! call memocc(istat, iexpand, 'iexpand', subname)
      ! Get the array to transfrom back the data
      call get_reverse_indices(ndimind, iextract, iexpand)
    
    !!call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
    !!t2=mpi_wtime()
    !!tt=t2-t1
    !!if(iproc==0) write(*,*) 'time 5.3: iproc', iproc, tt
    
      allocate(indexrecvorbital2(ndimind), stat=istat)
      call memocc(istat, indexrecvorbital2, 'indexrecvorbital2', subname)
    
      if (ndimind>0) then
          call vcopy(ndimind, indexrecvorbital(1), 1, indexrecvorbital2(1), 1)
      end if
    
      !$omp parallel default(none) &
      !$omp shared(ndimind, iextract, indexrecvorbital, indexrecvorbital2) private(i, ind)
      !$omp do
      do i=1,ndimind
          ind=iextract(i)
          indexrecvorbital(ind)=indexrecvorbital2(i)
      end do
      !$omp end do
      !$omp end parallel
    
      iall=-product(shape(indexrecvorbital2))*kind(indexrecvorbital2)
      deallocate(indexrecvorbital2, stat=istat)
      call memocc(istat, iall, 'indexrecvorbital2', subname)
    
      if(minval(indexrecvorbital)<1) stop 'minval(indexrecvorbital)<1'
      if(maxval(indexrecvorbital)>orbs%norb) stop 'maxval(indexrecvorbital)>orbs%norb'
    
    
      iall=-product(shape(gridpoint_start))*kind(gridpoint_start)
      deallocate(gridpoint_start, stat=istat)
      call memocc(istat, iall, 'gridpoint_start', subname)
    
      iall=-product(shape(nsend))*kind(nsend)
      deallocate(nsend, stat=istat)
      call memocc(istat, iall, 'nsend', subname)
    
    
    end subroutine get_switch_indices_sumrho



    subroutine communication_arrays_repartitionrho(iproc, nproc, lzd, nscatterarr, istartend, &
               nsendcounts_repartitionrho, nsenddspls_repartitionrho, &
               nrecvcounts_repartitionrho, nrecvdspls_repartitionrho)
      use module_base
      use module_types
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(local_zone_descriptors),intent(in) :: lzd
      integer,dimension(0:nproc-1,4),intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
      integer,dimension(2,0:nproc-1),intent(in) :: istartend
      integer,dimension(0:nproc-1),intent(out) :: nsendcounts_repartitionrho, nsenddspls_repartitionrho
      integer,dimension(0:nproc-1),intent(out) :: nrecvcounts_repartitionrho, nrecvdspls_repartitionrho
    
      ! Local variables
      integer :: jproc_send, jproc_recv, ii, i3, i2, i1, jproc
    
      jproc_send=0
      jproc_recv=0
      ii=0
      nsendcounts_repartitionrho=0
      nrecvcounts_repartitionrho=0
      do i3=1,lzd%glr%d%n3i
          do i2=1,lzd%glr%d%n2i
              do i1=1,lzd%glr%d%n1i
                  ii=ii+1
                  if (ii>istartend(2,jproc_send)) then
                      jproc_send=jproc_send+1
                  end if
                  if (i3>nscatterarr(jproc_recv,3)+nscatterarr(jproc_recv,2)) then
                      jproc_recv=jproc_recv+1
                  end if
                  if (iproc==jproc_send) then
                      nsendcounts_repartitionrho(jproc_recv)=nsendcounts_repartitionrho(jproc_recv)+1
                  end if
                  if (iproc==jproc_recv) then
                      nrecvcounts_repartitionrho(jproc_send)=nrecvcounts_repartitionrho(jproc_send)+1
                  end if
              end do
          end do
      end do
    
      nsenddspls_repartitionrho(0)=0
      nrecvdspls_repartitionrho(0)=0
      do jproc=1,nproc-1
          nsenddspls_repartitionrho(jproc)=nsenddspls_repartitionrho(jproc-1)+&
                                                      nsendcounts_repartitionrho(jproc-1)
          nrecvdspls_repartitionrho(jproc)=nrecvdspls_repartitionrho(jproc-1)+&
                                                      nrecvcounts_repartitionrho(jproc-1)
      end do
    
    
    end subroutine communication_arrays_repartitionrho
    
    
    
    subroutine communication_arrays_repartitionrho_general(iproc, nproc, lzd, nscatterarr, istartend, &
               ncomms_repartitionrho, commarr_repartitionrho)
      use module_base
      use module_types
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(local_zone_descriptors),intent(in) :: lzd
      integer,dimension(0:nproc-1,4),intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
      integer,dimension(2,0:nproc-1),intent(in) :: istartend
      integer,intent(out) :: ncomms_repartitionrho
      integer,dimension(:,:),pointer,intent(out) :: commarr_repartitionrho
      character(len=*),parameter :: subname='communication_arrays_repartitionrho_general'
    
      ! Local variables
      integer :: i1, i2, i3, ii, jproc, jproc_send, iidest, nel, ioverlaps, istat, ierr, iassign
      logical :: started
      integer,dimension(:),allocatable :: nel_array
    
      call f_routine(id='nel_array')
    
      ! only do this if task iproc has to receive a part of the potential
      if (nscatterarr(iproc,1)>0) then
        
        
          ! First process from which iproc has to receive data
          ncomms_repartitionrho=0
          i3=nscatterarr(iproc,3)-nscatterarr(iproc,4)
          ii=(i3)*(lzd%glr%d%n2i)*(lzd%glr%d%n1i)+1
          do jproc=nproc-1,0,-1
              if (ii>=istartend(1,jproc)) then
                  jproc_send=jproc
                  ncomms_repartitionrho=ncomms_repartitionrho+1
                  exit
              end if
          end do
        
        
        
          ! The remaining processes
          iidest=0
          nel=0
          started=.false.
          do i3=nscatterarr(iproc,3)-nscatterarr(iproc,4)+1,nscatterarr(iproc,3)-nscatterarr(iproc,4)+nscatterarr(iproc,1)
              ii=(i3-1)*(lzd%glr%d%n2i)*(lzd%glr%d%n1i)
              do i2=1,lzd%glr%d%n2i
                  do i1=1,lzd%glr%d%n1i
                      ii=ii+1
                      iidest=iidest+1
                      if (ii>=istartend(1,jproc_send) .and. ii<=istartend(2,jproc_send)) then
                          nel=nel+1
                      else
                          jproc_send=jproc_send+1
                          ncomms_repartitionrho=ncomms_repartitionrho+1
                      end if
                  end do
              end do
          end do
        
        
          call allocate_MPI_communications_arrays_repartitionp2p(ncomms_repartitionrho, commarr_repartitionrho)
        
        
          ! First process from which iproc has to receive data
          ioverlaps=0
          i3=nscatterarr(iproc,3)-nscatterarr(iproc,4)
          ii=(i3)*(lzd%glr%d%n2i)*(lzd%glr%d%n1i)+1
          do jproc=nproc-1,0,-1
              if (ii>=istartend(1,jproc)) then
                  jproc_send=jproc
                  ioverlaps=ioverlaps+1
                  exit
              end if
          end do
        
        
          ! The remaining processes
          iassign=0
          iidest=0
          nel=0
          started=.false.
          do i3=nscatterarr(iproc,3)-nscatterarr(iproc,4)+1,nscatterarr(iproc,3)-nscatterarr(iproc,4)+nscatterarr(iproc,1)
              ii=(i3-1)*(lzd%glr%d%n2i)*(lzd%glr%d%n1i)
              do i2=1,lzd%glr%d%n2i
                  do i1=1,lzd%glr%d%n1i
                      ii=ii+1
                      iidest=iidest+1
                      if (ii>=istartend(1,jproc_send) .and. ii<=istartend(2,jproc_send)) then
                          nel=nel+1
                      else
                          commarr_repartitionrho(4,ioverlaps)=nel
                          jproc_send=jproc_send+1
                          ioverlaps=ioverlaps+1
                          nel=1
                          started=.false.
                      end if
                      if (.not.started) then
                          if (jproc_send>=nproc) stop 'ERROR: jproc_send>=nproc'
                          commarr_repartitionrho(1,ioverlaps)=jproc_send
                          commarr_repartitionrho(2,ioverlaps)=ii-istartend(1,jproc_send)+1
                          commarr_repartitionrho(3,ioverlaps)=iidest
                          started=.true.
                          iassign=iassign+1
                      end if
                  end do
              end do
          end do
          commarr_repartitionrho(4,ioverlaps)=nel
          if (ioverlaps/=ncomms_repartitionrho) stop 'ERROR: ioverlaps/=ncomms_repartitionrho'
          if (iassign/=ncomms_repartitionrho) stop 'ERROR: iassign/=ncomms_repartitionrho'
        
          ! some checks
          nel=0
          !nel_array=f_malloc0(0.to.nproc-1,id='nel_array')
          do ioverlaps=1,ncomms_repartitionrho
              nel=nel+commarr_repartitionrho(4,ioverlaps)
              ii=commarr_repartitionrho(1,ioverlaps)
              !nel_array(ii)=nel_array(ii)+commarr_repartitionrho(4,ioverlaps)
          end do
          if (nel/=nscatterarr(iproc,1)*lzd%glr%d%n2i*lzd%glr%d%n1i) then
              stop 'nel/=nscatterarr(iproc,2)*lzd%glr%d%n2i*lzd%glr%d%n1i'
          end if
          !!call mpiallred(nel_array(0), nproc, mpi_sum, bigdft_mpi%mpi_comm, ierr)
          !!if (nel_array(iproc)/=istartend(2,iproc)-istartend(1,iproc)+1) then
          !!    !stop 'nel_array(iproc)/=istartend(2,iproc)-istartend(1,iproc)+1'
          !!end if
          !!call f_free(nel_array)
    
      else
          ncomms_repartitionrho=0
          call allocate_MPI_communications_arrays_repartitionp2p(1, commarr_repartitionrho)
    
      end if
    
    
    end subroutine communication_arrays_repartitionrho_general




    subroutine transpose_switch_psir(collcom_sr, psir, psirwork)
      use module_base
      use module_types
      implicit none
    
      ! Calling arguments
      type(collective_comms),intent(in) :: collcom_sr
      real(kind=8),dimension(collcom_sr%ndimpsi_c),intent(in) :: psir
      real(kind=8),dimension(collcom_sr%ndimpsi_c),intent(out) :: psirwork
    
      ! Local variables
      integer :: i, m, ind
    
      !$omp parallel default(private) &
      !$omp shared(collcom_sr, psir, psirwork, m)
    
      m = mod(collcom_sr%ndimpsi_c,7)
      if(m/=0) then
          do i=1,m
              ind = collcom_sr%isendbuf_c(i)
              psirwork(ind) = psir(i)
          end do
      end if
      !$omp do
      do i = m+1,collcom_sr%ndimpsi_c,7
         psirwork(collcom_sr%isendbuf_c(i+0)) = psir(i+0)
         psirwork(collcom_sr%isendbuf_c(i+1)) = psir(i+1)
         psirwork(collcom_sr%isendbuf_c(i+2)) = psir(i+2)
         psirwork(collcom_sr%isendbuf_c(i+3)) = psir(i+3)
         psirwork(collcom_sr%isendbuf_c(i+4)) = psir(i+4)
         psirwork(collcom_sr%isendbuf_c(i+5)) = psir(i+5)
         psirwork(collcom_sr%isendbuf_c(i+6)) = psir(i+6)
      end do
      !$omp end do
      !$omp end parallel
    
    
    end subroutine transpose_switch_psir
    
    subroutine transpose_communicate_psir(iproc, nproc, collcom_sr, psirwork, psirtwork)
      use module_base
      use module_types
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(collective_comms),intent(in) :: collcom_sr
      real(kind=8),dimension(collcom_sr%ndimpsi_c),intent(in) :: psirwork
      real(kind=8),dimension(collcom_sr%ndimind_c),intent(out) :: psirtwork
    
      ! Local variables
      integer :: ierr
    
    
      if (nproc>1) then
          call mpi_alltoallv(psirwork, collcom_sr%nsendcounts_c, collcom_sr%nsenddspls_c, mpi_double_precision, psirtwork, &
               collcom_sr%nrecvcounts_c, collcom_sr%nrecvdspls_c, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
      else
          call vcopy(collcom_sr%ndimpsi_c, psirwork(1), 1, psirtwork(1), 1)
      end if
    
    
    end subroutine transpose_communicate_psir
    
    subroutine transpose_unswitch_psirt(collcom_sr, psirtwork, psirt)
      use module_base
      use module_types
      implicit none
    
      ! Calling arguments
      type(collective_comms),intent(in) :: collcom_sr
      real(kind=8),dimension(collcom_sr%ndimind_c),intent(in) :: psirtwork
      real(kind=8),dimension(collcom_sr%ndimind_c),intent(out) :: psirt
    
      ! Local variables
      integer :: i, ind, sum_c, m
    
      sum_c = sum(collcom_sr%nrecvcounts_c)
    
      !$omp parallel private(i,ind) &
      !$omp shared(psirt, psirtwork, collcom_sr, sum_c, m)
    
      m = mod(sum_c,7)
    
      if(m/=0) then
        do i = 1,m
          ind=collcom_sr%iextract_c(i)
          psirt(ind)=psirtwork(i)
        end do
      end if
    
      !$omp do
      do i=m+1, sum_c,7
          psirt(collcom_sr%iextract_c(i+0))=psirtwork(i+0)
          psirt(collcom_sr%iextract_c(i+1))=psirtwork(i+1)
          psirt(collcom_sr%iextract_c(i+2))=psirtwork(i+2)
          psirt(collcom_sr%iextract_c(i+3))=psirtwork(i+3)
          psirt(collcom_sr%iextract_c(i+4))=psirtwork(i+4)
          psirt(collcom_sr%iextract_c(i+5))=psirtwork(i+5)
          psirt(collcom_sr%iextract_c(i+6))=psirtwork(i+6)
      end do
      !$omp end do
      !$omp end parallel
    
    end subroutine transpose_unswitch_psirt
    
    subroutine transpose_switch_psirt(collcom_sr, psirt, psirtwork)
      use module_base
      use module_types
      implicit none
    
      ! Calling arguments
      type(collective_comms),intent(in) :: collcom_sr
      real(kind=8),dimension(collcom_sr%ndimind_c),intent(in) :: psirt
      real(kind=8),dimension(collcom_sr%ndimind_c),intent(out) :: psirtwork
    
      ! Local variables
      integer :: i, ind, sum_c, m
    
      sum_c = sum(collcom_sr%nrecvcounts_c)
    
      !$omp parallel default(private) &
      !$omp shared(collcom_sr, psirt, psirtwork, sum_c, m)
    
      m = mod(sum_c,7)
    
      if(m/=0) then
        do i=1,m
           ind = collcom_sr%iexpand_c(i)
           psirtwork(ind) = psirt(i)
        end do
      end if
    
    
      !$omp do
      do i=m+1,sum_c,7
          psirtwork(collcom_sr%iexpand_c(i+0))=psirt(i+0)
          psirtwork(collcom_sr%iexpand_c(i+1))=psirt(i+1)
          psirtwork(collcom_sr%iexpand_c(i+2))=psirt(i+2)
          psirtwork(collcom_sr%iexpand_c(i+3))=psirt(i+3)
          psirtwork(collcom_sr%iexpand_c(i+4))=psirt(i+4)
          psirtwork(collcom_sr%iexpand_c(i+5))=psirt(i+5)
          psirtwork(collcom_sr%iexpand_c(i+6))=psirt(i+6)
      end do
      !$omp end do
      !$omp end parallel
    
    end subroutine transpose_switch_psirt
    
    subroutine transpose_communicate_psirt(iproc, nproc, collcom_sr, psirtwork, psirwork)
      use module_base
      use module_types
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(collective_comms),intent(in) :: collcom_sr
      real(kind=8),dimension(collcom_sr%ndimind_c),intent(in) :: psirtwork
      real(kind=8),dimension(collcom_sr%ndimpsi_c),intent(out) :: psirwork
    
      ! Local variables
      integer :: ierr
    
      if (nproc>1) then
      call mpi_alltoallv(psirtwork, collcom_sr%nrecvcounts_c, collcom_sr%nrecvdspls_c, mpi_double_precision, psirwork, &
           collcom_sr%nsendcounts_c, collcom_sr%nsenddspls_c, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
      else
          call vcopy(collcom_sr%ndimpsi_c, psirtwork(1), 1, psirwork(1), 1)
      end if
    
    end subroutine transpose_communicate_psirt
    
    subroutine transpose_unswitch_psir(collcom_sr, psirwork, psir)
      use module_base
      use module_types
      implicit none
    
      ! Caling arguments
      type(collective_comms),intent(in) :: collcom_sr
      real(kind=8),dimension(collcom_sr%ndimpsi_c),intent(in) :: psirwork
      real(kind=8),dimension(collcom_sr%ndimpsi_c),intent(out) :: psir
    
      ! Local variables
      integer :: i, ind, m
    
    
      !$omp parallel default(private) &
      !$omp shared(collcom_sr, psirwork, psir, m)
    
      m = mod(collcom_sr%ndimpsi_c,7)
    
      if(m/=0) then
        do i = 1,m
         ind=collcom_sr%irecvbuf_c(i)
         psir(ind)=psirwork(i)
        end do
      end if
    
      ! coarse part
    
      !$omp do
        do i=m+1,collcom_sr%ndimpsi_c,7
            psir(collcom_sr%irecvbuf_c(i+0))=psirwork(i+0)
            psir(collcom_sr%irecvbuf_c(i+1))=psirwork(i+1)
            psir(collcom_sr%irecvbuf_c(i+2))=psirwork(i+2)
            psir(collcom_sr%irecvbuf_c(i+3))=psirwork(i+3)
            psir(collcom_sr%irecvbuf_c(i+4))=psirwork(i+4)
            psir(collcom_sr%irecvbuf_c(i+5))=psirwork(i+5)
            psir(collcom_sr%irecvbuf_c(i+6))=psirwork(i+6)
        end do
      !$omp end do
      !$omp end parallel
    
    end subroutine transpose_unswitch_psir

end module communications
