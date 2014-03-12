module communications
  use module_base
  use module_types
  implicit none

  private


  public :: collective_comms_null
  public :: init_collective_comms
  public :: allocate_MPI_communications_arrays_repartition
  public :: get_reverse_indices

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

end module communications
