!> @file
! Intialization of the collective communications for the linear version
! @author
!    Copyright (C) 2011-2013 BigDFT group
!    This file is distributed under the terms of the
!    GNU General Public License, see ~/COPYING file
!    or http://www.gnu.org/copyleft/gpl.txt .
!    For the list of contributors, see ~/AUTHORS


subroutine init_collective_comms(iproc, nproc, npsidim_orbs, orbs, lzd, collcom, collcom_reference)
  use module_base
  use module_types
  use module_interfaces, except_this_one => init_collective_comms
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, npsidim_orbs
  type(orbitals_data),intent(in) :: orbs
  type(local_zone_descriptors),intent(in) :: lzd
  type(collective_comms),intent(inout) :: collcom
  type(collective_comms),optional,intent(in) :: collcom_reference
  
  ! Local variables
  integer :: ii, istat, iorb, iiorb, ilr, iall, istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, ierr
  integer :: ipt, nvalp_c, nvalp_f
  real(kind=8),dimension(:,:,:),allocatable :: weight_c, weight_f
  real(kind=8) :: weight_c_tot, weight_f_tot, weightp_c, weightp_f, tt, t1, t2
  integer,dimension(:,:),allocatable :: istartend_c, istartend_f
  integer,dimension(:,:,:),allocatable :: index_in_global_c, index_in_global_f
  integer,dimension(:),allocatable :: npts_par_c, npts_par_f
  character(len=*),parameter :: subname='init_collective_comms'
  
  call timing(iproc,'init_collcomm ','ON')
  
  allocate(weight_c(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3), stat=istat)
  call memocc(istat, weight_c, 'weight_c', subname)
  allocate(weight_f(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3), stat=istat)
  call memocc(istat, weight_f, 'weight_f', subname)
  allocate(index_in_global_c(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3), stat=istat)
  call memocc(istat, index_in_global_c, 'index_in_global_c', subname)
  allocate(index_in_global_f(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3), stat=istat)
  call memocc(istat, index_in_global_f, 'index_in_global_f', subname)


call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
t1=mpi_wtime()
call get_weights(iproc, nproc, orbs, lzd, weight_c, weight_f, weight_c_tot, weight_f_tot)
call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
t2=mpi_wtime()
!!if(iproc==0) write(*,'(a,es10.3)') 'time for part 1:',t2-t1
t1=mpi_wtime()
  ! Assign the grid points to the processes such that the work is equally dsitributed
  allocate(istartend_c(2,0:nproc-1), stat=istat)
  call memocc(istat, istartend_c, 'istartend_c', subname)
  allocate(istartend_f(2,0:nproc-1), stat=istat)
  call memocc(istat, istartend_f, 'istartend_f', subname)
  if(.not.present(collcom_reference)) then
      call assign_weight_to_process(iproc, nproc, lzd, weight_c, weight_f, weight_c_tot, weight_f_tot, &
           istartend_c, istartend_f, istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, &
           weightp_c, weightp_f, collcom%nptsp_c, collcom%nptsp_f, nvalp_c, nvalp_f)
  else
      stop 'THIS OPTION IS DEPRECTAED'
      !!allocate(npts_par_c(0:nproc-1), stat=istat)
      !! call memocc(istat, npts_par_c, 'npts_par_c', subname)
      !!allocate(npts_par_f(0:nproc-1), stat=istat)
      !! call memocc(istat, npts_par_f, 'npts_par_f', subname)
      !!npts_par_c=0
      !!npts_par_f=0
      !!npts_par_c(iproc)=collcom_reference%nptsp_c
      !!npts_par_f(iproc)=collcom_reference%nptsp_f
      !!call mpiallred(npts_par_c(0), nproc, mpi_sum, bigdft_mpi%mpi_comm, ierr)
      !!call mpiallred(npts_par_f(0), nproc, mpi_sum, bigdft_mpi%mpi_comm, ierr)
      !!call assign_weight_to_process2(iproc, nproc, lzd, weight_c, weight_f, weight_c_tot, weight_f_tot, &
      !!     npts_par_c, npts_par_f, &
      !!     istartend_c, istartend_f, istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, &
      !!     weightp_c, weightp_f, collcom%nptsp_c, collcom%nptsp_f)
      !!iall=-product(shape(npts_par_c))*kind(npts_par_c)
      !!deallocate(npts_par_c, stat=istat)
      !!call memocc(istat, iall, 'npts_par_c', subname)
      !!iall=-product(shape(npts_par_f))*kind(npts_par_f)
      !!deallocate(npts_par_f, stat=istat)
      !!call memocc(istat, iall, 'npts_par_f', subname)
  end if


call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
t2=mpi_wtime()
!!if(iproc==0) write(*,'(a,es10.3)') 'time for part 2:',t2-t1
t1=mpi_wtime()

  ! some checks
  if(nproc>1) then
      call mpi_allreduce(weightp_c, tt, 1, mpi_double_precision, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  else
      tt=weightp_c
  end if
  if(tt/=weight_c_tot) stop 'wrong partition of coarse weights'
  if(nproc>1) then
      call mpi_allreduce(weightp_f, tt, 1, mpi_double_precision, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  else
      tt=weightp_f
  end if     
  if(tt/=weight_f_tot) stop 'wrong partition of fine weights'
  if(nproc>1) then
      call mpi_allreduce(collcom%nptsp_c, ii, 1, mpi_integer, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  else
      ii=collcom%nptsp_c
  end if
  if(ii/=lzd%glr%wfd%nvctr_c) stop 'wrong partition of coarse grid points'
  if(nproc>1) then
      call mpi_allreduce(collcom%nptsp_f, ii, 1, mpi_integer, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  else
      ii=collcom%nptsp_f
  end if
  if(ii/=lzd%glr%wfd%nvctr_f) stop 'init_collective_comms: wrong partition of fine grid points'

  ! Allocate the keys
  allocate(collcom%norb_per_gridpoint_c(collcom%nptsp_c), stat=istat)
  call memocc(istat, collcom%norb_per_gridpoint_c, 'collcom%norb_per_gridpoint_c', subname)
  allocate(collcom%norb_per_gridpoint_f(collcom%nptsp_f), stat=istat)
  call memocc(istat, collcom%norb_per_gridpoint_f, 'collcom%norb_per_gridpoint_f', subname)
  call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
  call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
  t1=mpi_wtime()
  !call determine_num_orbs_per_gridpoint(iproc, nproc, orbs, lzd, istartend_c, istartend_f, &
  !     istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, &
  !     weightp_c, weightp_f, collcom%nptsp_c, collcom%nptsp_f, &
  !     collcom%norb_per_gridpoint_c, collcom%norb_per_gridpoint_f)
  call determine_num_orbs_per_gridpoint_new(iproc, nproc, lzd, istartend_c, istartend_f, &
       istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, &
       weightp_c, weightp_f, collcom%nptsp_c, collcom%nptsp_f, weight_c, weight_f, &
       collcom%norb_per_gridpoint_c, collcom%norb_per_gridpoint_f)
call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
t2=mpi_wtime()
!!if(iproc==0) write(*,'(a,es10.3)') 'time for part 3:',t2-t1
t1=mpi_wtime()

  ! Determine the index of a grid point i1,i2,i3 in the compressed array
  call get_index_in_global2(lzd%glr, index_in_global_c, index_in_global_f)
call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
t2=mpi_wtime()
!!if(iproc==0) write(*,'(a,es10.3)') 'time for part 4:',t2-t1
t1=mpi_wtime()




  iall=-product(shape(weight_c))*kind(weight_c)
  deallocate(weight_c, stat=istat)
  call memocc(istat, iall, 'weight_c', subname)
  iall=-product(shape(weight_f))*kind(weight_f)
  deallocate(weight_f, stat=istat)
  call memocc(istat, iall, 'weight_f', subname)

  ! Determine values for mpi_alltoallv
  allocate(collcom%nsendcounts_c(0:nproc-1), stat=istat)
  call memocc(istat, collcom%nsendcounts_c, 'collcom%nsendcounts_c', subname)
  allocate(collcom%nsenddspls_c(0:nproc-1), stat=istat)
  call memocc(istat, collcom%nsenddspls_c, 'collcom%nsenddspls_c', subname)
  allocate(collcom%nrecvcounts_c(0:nproc-1), stat=istat)
  call memocc(istat, collcom%nrecvcounts_c, 'collcom%nrecvcounts_c', subname)
  allocate(collcom%nrecvdspls_c(0:nproc-1), stat=istat)
  call memocc(istat, collcom%nrecvdspls_c, 'collcom%nrecvdspls_c', subname)
  allocate(collcom%nsendcounts_f(0:nproc-1), stat=istat)
  call memocc(istat, collcom%nsendcounts_f, 'collcom%nsendcounts_f', subname)
  allocate(collcom%nsenddspls_f(0:nproc-1), stat=istat)
  call memocc(istat, collcom%nsenddspls_f, 'collcom%nsenddspls_f', subname)
  allocate(collcom%nrecvcounts_f(0:nproc-1), stat=istat)
  call memocc(istat, collcom%nrecvcounts_f, 'collcom%nrecvcounts_f', subname)
  allocate(collcom%nrecvdspls_f(0:nproc-1), stat=istat)
  call memocc(istat, collcom%nrecvdspls_f, 'collcom%nrecvdspls_f', subname)
  call determine_communication_arrays(iproc, nproc, npsidim_orbs, orbs, lzd, istartend_c, istartend_f, &
       index_in_global_c, index_in_global_f, nvalp_c, nvalp_f, &
       collcom%nsendcounts_c, collcom%nsenddspls_c, collcom%nrecvcounts_c, collcom%nrecvdspls_c, &
       collcom%nsendcounts_f, collcom%nsenddspls_f, collcom%nrecvcounts_f, collcom%nrecvdspls_f)


  call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
  t2=mpi_wtime()
!!if(iproc==0) write(*,'(a,es10.3)') 'time for part 5:',t2-t1
t1=mpi_wtime()


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
  allocate(collcom%irecvbuf_c(collcom%ndimpsi_c), stat=istat)
  call memocc(istat, collcom%irecvbuf_c, 'collcom%irecvbuf_c', subname)
  allocate(collcom%indexrecvorbital_c(collcom%ndimind_c), stat=istat)
  call memocc(istat, collcom%indexrecvorbital_c, 'collcom%indexrecvorbital_c', subname)
  allocate(collcom%iextract_c(collcom%ndimind_c), stat=istat)
  call memocc(istat, collcom%iextract_c, 'collcom%iextract_c', subname)
  allocate(collcom%iexpand_c(collcom%ndimind_c), stat=istat)
  call memocc(istat, collcom%iexpand_c, 'collcom%iexpand_c', subname)
  allocate(collcom%isendbuf_c(collcom%ndimpsi_c), stat=istat)
  call memocc(istat, collcom%isendbuf_c, 'collcom%isendbuf_c', subname)

  collcom%ndimpsi_f=0
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      ilr=orbs%inwhichlocreg(iiorb)
      collcom%ndimpsi_f=collcom%ndimpsi_f+lzd%llr(ilr)%wfd%nvctr_f
  end do
  allocate(collcom%irecvbuf_f(collcom%ndimpsi_f), stat=istat)
  call memocc(istat, collcom%irecvbuf_f, 'collcom%irecvbuf_f', subname)
  allocate(collcom%indexrecvorbital_f(collcom%ndimind_f), stat=istat)
  call memocc(istat, collcom%indexrecvorbital_f, 'collcom%indexrecvorbital_f', subname)
  allocate(collcom%iextract_f(collcom%ndimind_f), stat=istat)
  call memocc(istat, collcom%iextract_f, 'collcom%iextract_f', subname)
  allocate(collcom%iexpand_f(collcom%ndimind_f), stat=istat)
  call memocc(istat, collcom%iexpand_f, 'collcom%iexpand_f', subname)
  allocate(collcom%isendbuf_f(collcom%ndimpsi_f), stat=istat)
  call memocc(istat, collcom%isendbuf_f, 'collcom%isendbuf_f', subname)


  call get_switch_indices(iproc, nproc, orbs, lzd, collcom%ndimpsi_c, collcom%ndimpsi_f, istartend_c, istartend_f, &
       collcom%nsendcounts_c, collcom%nsenddspls_c, collcom%ndimind_c, collcom%nrecvcounts_c, collcom%nrecvdspls_c, &
       collcom%nsendcounts_f, collcom%nsenddspls_f, collcom%ndimind_f, collcom%nrecvcounts_f, collcom%nrecvdspls_f, &
       index_in_global_c, index_in_global_f, &
       weightp_c, weightp_f, collcom%isendbuf_c, collcom%irecvbuf_c, collcom%isendbuf_f, collcom%irecvbuf_f, &
       collcom%indexrecvorbital_c, collcom%iextract_c, collcom%iexpand_c, &
       collcom%indexrecvorbital_f, collcom%iextract_f, collcom%iexpand_f)

  call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
  t2=mpi_wtime()
!!if(iproc==0) write(*,'(a,es10.3)') 'time for part 6:',t2-t1
t1=mpi_wtime()

  ! These variables are used in various subroutines to speed up the code
  allocate(collcom%isptsp_c(max(collcom%nptsp_c,1)), stat=istat)
  call memocc(istat, collcom%isptsp_c, 'collcom%isptsp_c', subname)
  allocate(collcom%isptsp_f(max(collcom%nptsp_f,1)), stat=istat)
  call memocc(istat, collcom%isptsp_f, 'collcom%isptsp_f', subname)
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
  iall=-product(shape(istartend_c))*kind(istartend_c)
  deallocate(istartend_c, stat=istat)
  call memocc(istat, iall, 'istartend_c', subname)

  iall=-product(shape(istartend_f))*kind(istartend_f)
  deallocate(istartend_f, stat=istat)
  call memocc(istat, iall, 'istartend_f', subname)

  iall=-product(shape(index_in_global_c))*kind(index_in_global_c)
  deallocate(index_in_global_c, stat=istat)
  call memocc(istat, iall, 'index_in_global_c', subname)

  iall=-product(shape(index_in_global_f))*kind(index_in_global_f)
  deallocate(index_in_global_f, stat=istat)
  call memocc(istat, iall, 'index_in_global_f', subname)

  
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



!!$$  ! First the coarse part...
!!$$ 
!!$$
!!$$  !$omp parallel default(private) shared(lzd,iproc,nproc)&
!!$$  !$omp shared(weight_f,weight_f_ideal,weight_tot_f,weight_c_ideal,weight_tot_c, weight_c,istartend_c,istartend_f)&
!!$$  !$omp shared(istartp_seg_c,iendp_seg_c,istartp_seg_f,iendp_seg_f,weightp_c,weightp_f,nptsp_c,nptsp_f)
!!$$
!!$$  !$omp sections 
!!$$
!!$$  !$omp section
!!$$
!!$$  jproc=0
!!$$  tt=0.d0
!!$$  tt2=0.d0
!!$$  iitot=0
!!$$  ii2=0
!!$$  iiseg=1
!!$$  jprocdone=-1
!!$$  weightp_c=0.d0
!!$$  loop_nseg_c: do iseg=1,lzd%glr%wfd%nseg_c
!!$$       jj=lzd%glr%wfd%keyvloc(iseg)
!!$$       j0=lzd%glr%wfd%keygloc(1,iseg)
!!$$       j1=lzd%glr%wfd%keygloc(2,iseg)
!!$$       ii=j0-1
!!$$       i3=ii/((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1))
!!$$       ii=ii-i3*(lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)
!!$$       i2=ii/(lzd%glr%d%n1+1)
!!$$       i0=ii-i2*(lzd%glr%d%n1+1)
!!$$       i1=i0+j1-j0
!!$$       do i=i0,i1
!!$$           tt=tt+weight_c(i,i2,i3)
!!$$           iitot=iitot+1
!!$$           if((tt>=weight_c_ideal .or. iseg==lzd%glr%wfd%nseg_c) .and. i==i1 .and. jproc<=nproc-2) then
!!$$               if(tt==weight_tot_c .and. jproc==nproc-1) then
!!$$                   ! this process also has to take the remaining points, even if they have no weight
!!$$                   iitot=lzd%glr%wfd%nvctr_c-ii2
!!$$               end if
!!$$               if(iproc==jproc) then
!!$$                   weightp_c=tt
!!$$                   nptsp_c=iitot
!!$$                   istartp_seg_c=iiseg
!!$$                   iendp_seg_c=iseg
!!$$                   if(tt==weight_tot_c .and. jproc==nproc-1) then
!!$$                       ! this process also has to take the remaining segments, even if they have no weight
!!$$                       iendp_seg_c=lzd%glr%wfd%nseg_c
!!$$                   end if
!!$$               end if
!!$$               istartend_c(1,jproc)=ii2+1
!!$$               istartend_c(2,jproc)=min(istartend_c(1,jproc)+iitot-1,lzd%glr%wfd%nvctr_c)
!!$$               tt2=tt2+tt
!!$$               tt=0.d0
!!$$               ii2=ii2+iitot
!!$$               iitot=0
!!$$               jproc=jproc+1
!!$$               iiseg=iseg
!!$$               if(ii2>=lzd%glr%wfd%nvctr_c) then
!!$$                   ! everything is distributed
!!$$                   jprocdone=jproc
!!$$                   exit loop_nseg_c
!!$$               end if
!!$$           end if
!!$$       end do
!!$$   end do loop_nseg_c
!!$$
!!$$  if(jprocdone>0) then
!!$$       do jproc=jprocdone,nproc-1
!!$$          ! these processes do nothing
!!$$          istartend_c(1,jproc)=lzd%glr%wfd%nvctr_c+1
!!$$          istartend_c(2,jproc)=lzd%glr%wfd%nvctr_c
!!$$          if(iproc==jproc) then
!!$$              weightp_c=0.d0
!!$$              nptsp_c=0
!!$$              istartp_seg_c=lzd%glr%wfd%nseg_c+1
!!$$              iendp_seg_c=lzd%glr%wfd%nseg_c
!!$$          end if
!!$$      end do
!!$$  else
!!$$      if(iproc==nproc-1) then
!!$$          ! Take the rest
!!$$          weightp_c=weight_tot_c-tt2
!!$$          nptsp_c=lzd%glr%wfd%nvctr_c-ii2
!!$$          istartp_seg_c=iiseg
!!$$          iendp_seg_c=lzd%glr%wfd%nseg_c
!!$$      end if
!!$$      istartend_c(1,nproc-1)=ii2+1
!!$$      istartend_c(2,nproc-1)=istartend_c(1,nproc-1)+iitot-1
!!$$  end if
!!$$
!!$$  ! some check
!!$$ 
!!$$  !write(*,*) 'subroutine', weightp_c
!!$$
!!$$  !$omp section
!!$$
!!$$  ! Now the fine part...
!!$$  jproc=0
!!$$  tt=0.d0
!!$$  tt2=0.d0
!!$$  iitot=0
!!$$  ii2=0
!!$$  weightp_f=0.d0
!!$$  istart=lzd%glr%wfd%nseg_c+min(1,lzd%glr%wfd%nseg_f)
!!$$  iend=istart+lzd%glr%wfd%nseg_f-1
!!$$  iiseg=istart
!!$$  jprocdone=-1
!!$$  loop_nseg_f: do iseg=istart,iend
!!$$       jj=lzd%glr%wfd%keyvloc(iseg)
!!$$       j0=lzd%glr%wfd%keygloc(1,iseg)
!!$$       j1=lzd%glr%wfd%keygloc(2,iseg)
!!$$       ii=j0-1
!!$$       i3=ii/((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1))
!!$$       ii=ii-i3*(lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)
!!$$       i2=ii/(lzd%glr%d%n1+1)
!!$$       i0=ii-i2*(lzd%glr%d%n1+1)
!!$$       i1=i0+j1-j0
!!$$       do i=i0,i1
!!$$           tt=tt+weight_f(i,i2,i3)
!!$$           iitot=iitot+1
!!$$           if((tt>=weight_f_ideal .or. iseg==iend) .and. i==i1 .and. jproc<=nproc-2) then
!!$$               if(tt==weight_tot_f .and. jproc==nproc-1) then
!!$$                   ! this process also has to take the remaining points, even if they have no weight
!!$$                   iitot=lzd%glr%wfd%nvctr_f-ii2
!!$$               end if
!!$$               if(iproc==jproc) then
!!$$                   weightp_f=tt
!!$$                   nptsp_f=iitot
!!$$                   istartp_seg_f=iiseg
!!$$                   iendp_seg_f=iseg
!!$$                   if(tt==weight_tot_f .and. jproc==nproc-1) then
!!$$                       ! this process also has to take the remaining segments, even if they have no weight
!!$$                       iendp_seg_f=iend
!!$$                   end if
!!$$               end if
!!$$               istartend_f(1,jproc)=ii2+1
!!$$               istartend_f(2,jproc)=min(istartend_f(1,jproc)+iitot-1,lzd%glr%wfd%nvctr_f)
!!$$               tt2=tt2+tt
!!$$               tt=0.d0
!!$$               ii2=ii2+iitot
!!$$               iitot=0
!!$$               jproc=jproc+1
!!$$               iiseg=iseg
!!$$               if(ii2>=lzd%glr%wfd%nvctr_f) then
!!$$                   ! everything is distributed
!!$$                   jprocdone=jproc
!!$$                   exit loop_nseg_f
!!$$               end if
!!$$           end if
!!$$       end do
!!$$   end do loop_nseg_f
!!$$  if(jprocdone>0) then
!!$$       do jproc=jprocdone,nproc-1
!!$$          ! these processes do nothing
!!$$          istartend_f(1,jproc)=lzd%glr%wfd%nvctr_f+1
!!$$          istartend_f(2,jproc)=lzd%glr%wfd%nvctr_f
!!$$          if(iproc==jproc) then
!!$$              weightp_f=0.d0
!!$$              nptsp_f=0
!!$$              istartp_seg_f=lzd%glr%wfd%nseg_f+1
!!$$              iendp_seg_f=lzd%glr%wfd%nseg_f
!!$$          end if
!!$$      end do
!!$$  else
!!$$      if(iproc==nproc-1) then
!!$$          ! Take the rest
!!$$          weightp_f=weight_tot_f-tt2
!!$$          nptsp_f=lzd%glr%wfd%nvctr_f-ii2
!!$$          istartp_seg_f=iiseg
!!$$          iendp_seg_f=iend
!!$$      end if
!!$$      istartend_f(1,nproc-1)=ii2+1
!!$$      istartend_f(2,nproc-1)=istartend_f(1,nproc-1)+iitot-1
!!$$  end if
!!$$
!!$$!$omp end sections
!!$$  !$omp end parallel


!!! NEW VERSION ####################################################
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
  
 

end subroutine assign_weight_to_process



! This subroutine distributes the components to the processes such that each process has the same grid points
! as indicated by npts_par_c, npts_par_f.
subroutine assign_weight_to_process2(iproc, nproc, lzd, weight_c, weight_f, weight_tot_c, weight_tot_f, &
           npts_par_c, npts_par_f, &
           istartend_c, istartend_f, istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, &
           weightp_c, weightp_f, nptsp_c, nptsp_f)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(local_zone_descriptors),intent(in) :: lzd
  real(kind=8),dimension(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3),intent(in) :: weight_c, weight_f
  real(kind=8),intent(in) :: weight_tot_c, weight_tot_f
  integer,dimension(0:nproc-1),intent(in) :: npts_par_c, npts_par_f
  integer,dimension(2,0:nproc-1),intent(out) :: istartend_c, istartend_f
  integer,intent(out) :: istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f
  real(kind=8),intent(out) :: weightp_c, weightp_f
  integer,intent(out) :: nptsp_c, nptsp_f
  
  ! Local variables
  integer :: jproc, i1, i2, i3, ii, istartp_c, iendp_c, ii2, istartp_f, iendp_f, istart, iend, jj, j0, j1
  integer :: i, iseg, i0, iitot, ierr, iiseg, jprocdone
  real(kind=8) :: tt, tt2, weight_c_ideal, weight_f_ideal

  ! Ideal weight per process
  weight_c_ideal=weight_tot_c/dble(nproc)
  weight_f_ideal=weight_tot_f/dble(nproc)

  jproc=0
  tt=0.d0
  tt2=0.d0
  iitot=0
  ii2=0
  iiseg=1
  weightp_c=0.d0
    do iseg=1,lzd%glr%wfd%nseg_c
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
           tt=tt+weight_c(i,i2,i3)
           iitot=iitot+1
           !if(tt>weight_c_ideal) then
           if(iitot==npts_par_c(jproc)) then
               if(iproc==jproc) then
                   weightp_c=tt
                   nptsp_c=iitot
                   istartp_c=ii2+1
                   iendp_c=istartp_c+iitot-1
                   istartp_seg_c=iiseg
                   iendp_seg_c=iseg
               end if
               istartend_c(1,jproc)=ii2+1
               istartend_c(2,jproc)=istartend_c(1,jproc)+iitot-1
               tt2=tt2+tt
               tt=0.d0
               ii2=ii2+iitot
               iitot=0
               jproc=jproc+1
               iiseg=iseg
           end if
       end do
   end do

   jprocdone=jproc
   do jproc=jprocdone,nproc-1
      ! these processes do nothing
      istartend_c(1,jproc)=lzd%glr%wfd%nvctr_c+1
      istartend_c(2,jproc)=lzd%glr%wfd%nvctr_c
      if(iproc==jproc) then
          weightp_c=0.d0
          nptsp_c=0
          istartp_seg_c=lzd%glr%wfd%nseg_c+1
          iendp_seg_c=lzd%glr%wfd%nseg_c
      end if
   end do
  !if(iproc==nproc-1) then
  !    ! Take the rest
  !    istartp_c=ii2+1
  !    iendp_c=istartp_c+iitot-1
  !    weightp_c=weight_tot_c-tt2
  !    nptsp_c=lzd%glr%wfd%nvctr_c-ii2
  !    istartp_seg_c=iiseg
  !    iendp_seg_c=lzd%glr%wfd%nseg_c
  !end if
  !istartend_c(1,nproc-1)=ii2+1
  !istartend_c(2,nproc-1)=istartend_c(1,nproc-1)+iitot-1

  ! some check
  ii=istartend_c(2,iproc)-istartend_c(1,iproc)+1
  if(nproc>1) call mpiallred(ii, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  if(ii/=lzd%glr%wfd%nvctr_c) stop 'assign_weight_to_process2: ii/=lzd%glr%wfd%nvctr_c'


  jproc=0
  tt=0.d0
  tt2=0.d0
  iitot=0
  ii2=0
  weightp_f=0.d0
  istart=lzd%glr%wfd%nseg_c+min(1,lzd%glr%wfd%nseg_f)
  iend=istart+lzd%glr%wfd%nseg_f-1
  iiseg=istart
    do iseg=istart,iend
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
           tt=tt+weight_f(i,i2,i3)
           iitot=iitot+1
           !if(tt>weight_f_ideal) then
           if(iitot==npts_par_f(jproc)) then
               if(iproc==jproc) then
                   weightp_f=tt
                   nptsp_f=iitot
                   istartp_f=ii2+1
                   iendp_f=istartp_f+iitot-1
                   istartp_seg_f=iiseg
                   iendp_seg_f=iseg
               end if
               istartend_f(1,jproc)=ii2+1
               istartend_f(2,jproc)=istartend_f(1,jproc)+iitot-1
               tt2=tt2+tt
               tt=0.d0
               ii2=ii2+iitot
               iitot=0
               jproc=jproc+1
               iiseg=iseg
           end if
       end do
   end do

   jprocdone=jproc
   do jproc=jprocdone,nproc-1
      ! these processes do nothing
      istartend_f(1,jproc)=lzd%glr%wfd%nvctr_f+1
      istartend_f(2,jproc)=lzd%glr%wfd%nvctr_f
      if(iproc==jproc) then
          weightp_f=0.d0
          nptsp_f=0
          istartp_seg_f=lzd%glr%wfd%nseg_f+1
          iendp_seg_f=lzd%glr%wfd%nseg_f
      end if
   end do
  !if(iproc==nproc-1) then
  !    ! Take the rest
  !    istartp_f=ii2+1
  !    iendp_f=istartp_f+iitot-1
  !    weightp_f=weight_tot_f-tt2
  !    nptsp_f=lzd%glr%wfd%nvctr_f-ii2
  !    istartp_seg_f=iiseg
  !    iendp_seg_f=iend
  !end if
  !istartend_f(1,nproc-1)=ii2+1
  !istartend_f(2,nproc-1)=istartend_f(1,nproc-1)+iitot-1

  ! some check
  ii=istartend_f(2,iproc)-istartend_f(1,iproc)+1
  if(nproc>1) call mpiallred(ii, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  if(ii/=lzd%glr%wfd%nvctr_f) stop 'assign_weight_to_process2: ii/=lzd%glr%wfd%nvctr_f'



end subroutine assign_weight_to_process2




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
                  do jproc=0,nproc-1
                      if(ind>=istartend_f(1,jproc) .and. ind<=istartend_f(2,jproc)) then
                          jproctarget=jproc
                          exit
                      end if
                  end do
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
              do jproc=0,nproc-1
                  if(indglob>=istartend_c(1,jproc) .and. indglob<=istartend_c(2,jproc)) then
                      jproctarget=jproc
                      exit
                  end if
              end do
              !write(600+iproc,'(a,2(i0,1x),i0,a,i0)') 'point ',ii1,ii2,ii3,' goes to process ',jproctarget
          
              nsend_c(jproctarget)=nsend_c(jproctarget)+1
              ind=nsenddspls_c(jproctarget)+nsend_c(jproctarget)
              isendbuf_c(iitot)=ind
              indexsendbuf_c(ind)=indglob
              indexsendorbital_c(iitot)=iiorb
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
                  do jproc=0,nproc-1
                      if(indglob>=istartend_f(1,jproc) .and. indglob<=istartend_f(2,jproc)) then
                          jproctarget=jproc
                          exit
                      end if
                  end do
                  nsend_f(jproctarget)=nsend_f(jproctarget)+1
                  ind=nsenddspls_f(jproctarget)+nsend_f(jproctarget)
                  isendbuf_f(iitot)=ind
                  indexsendbuf_f(ind)=indglob
                  indexsendorbital_f(iitot)=iiorb
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




subroutine check_gridpoint(nseg, n1, n2, noffset1, noffset2, noffset3, keyg, itarget1, itarget2, itarget3, iseg_start, found)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: nseg, n1, n2, noffset1, noffset2, noffset3, itarget1, itarget2, itarget3
  integer,dimension(2,nseg),intent(in) :: keyg
  integer,intent(inout) :: iseg_start
  logical,intent(out) :: found
  
  ! Local variables
  integer :: j0, j1, ii, i1, i2, i3, i0, ii1, ii2, ii3, iseg, i
  logical :: equal_possible, larger_possible, smaller_possible
  !integer :: iproc
  
  !call mpi_comm_rank(bigdft_mpi%mpi_comm, iproc, i)

  found=.false.
  loop_segments: do iseg=iseg_start,nseg
     j0=keyg(1,iseg)
     j1=keyg(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     ii2=i2+noffset2
     ii3=i3+noffset3
     equal_possible = (ii2==itarget2 .and. ii3==itarget3)
     larger_possible = (ii2>itarget2 .and. ii3>itarget3)
     smaller_possible = (ii2<itarget2 .and. ii3<itarget3)
     ! check whether quick exit is possible since there is no chance to find the point anymore...
     if(ii3>itarget3) then
            exit loop_segments
     end if
     if(ii3>=itarget3 .and. ii2>itarget2) then
         exit loop_segments
     end if
     larger_possible = (ii3>=itarget3 .and. ii2>=itarget2)

     do i=i0,i1
        ii1=i+noffset1
        if(equal_possible .and. ii1==itarget1) then
            found=.true.
            ! no need to search in smaller segments from now on, since the itargets will never decrease any more...
            iseg_start=iseg
            exit loop_segments
        end if
        if(larger_possible .and. ii1>itarget1) then
            ! there is no chance to find the point anymore...
            exit loop_segments
        end if
        if(smaller_possible .and. ii1<itarget1) then
            ! no need to search in these segments from now on, since the itargets will never decrease any more...
            iseg_start=iseg
        end if
     end do
  end do loop_segments


end subroutine check_gridpoint



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

         call dcopy(lzd%llr(ilr)%wfd%nvctr_c,psi(i_tot+1),1,psi_c(i_c+1),1)

         i_c = i_c + lzd%llr(ilr)%wfd%nvctr_c
         i_tot = i_tot + lzd%llr(ilr)%wfd%nvctr_c

         call dcopy(7*lzd%llr(ilr)%wfd%nvctr_f,psi(i_tot+1),1,psi_f(i_f+1),1)

         i_f = i_f + 7*lzd%llr(ilr)%wfd%nvctr_f
         i_tot = i_tot + 7*lzd%llr(ilr)%wfd%nvctr_f

      end do
    

  else
      ! only coarse part is used...
      call dcopy(collcom%ndimpsi_c, psi, 1, psi_c, 1)
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
  !!    if(collcom%nsendcounts_c(jproc)>0) call dcopy(collcom%nsendcounts_c(jproc), psiwork_c(ist_c), 1, psiwork(ist), 1)
  !!    ist_c=ist_c+collcom%nsendcounts_c(jproc)
  !!    ist=ist+collcom%nsendcounts_c(jproc)
  !!    if(collcom%nsendcounts_f(jproc)>0) call dcopy(7*collcom%nsendcounts_f(jproc), psiwork_f(ist_f), 1, psiwork(ist), 1)
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
  !!    if(collcom%nrecvcounts_c(jproc)>0) call dcopy(collcom%nrecvcounts_c(jproc), psitwork(ist), 1, psitwork_c(ist_c), 1)
  !!    ist_c=ist_c+collcom%nrecvcounts_c(jproc)
  !!    ist=ist+collcom%nrecvcounts_c(jproc)
  !!    if(collcom%nrecvcounts_f(jproc)>0) call dcopy(7*collcom%nrecvcounts_f(jproc), psitwork(ist), 1, psitwork_f(ist_f), 1)
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
  !!    if(collcom%nrecvcounts_c(jproc)>0) call dcopy(collcom%nrecvcounts_c(jproc), psitwork_c(ist_c), 1, psitwork(ist), 1)
  !!    ist_c=ist_c+collcom%nrecvcounts_c(jproc)
  !!    ist=ist+collcom%nrecvcounts_c(jproc)
  !!    if(collcom%nrecvcounts_f(jproc)>0) call dcopy(7*collcom%nrecvcounts_f(jproc), psitwork_f(ist_f), 1, psitwork(ist), 1)
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
  !!    if(collcom%nsendcounts_c(jproc)>0) call dcopy(collcom%nsendcounts_c(jproc), psiwork(ist), 1, psiwork_c(ist_c), 1)
  !!    ist_c=ist_c+collcom%nsendcounts_c(jproc)
  !!    ist=ist+collcom%nsendcounts_c(jproc)
  !!    if(collcom%nsendcounts_f(jproc)>0) call dcopy(7*collcom%nsendcounts_f(jproc), psiwork(ist), 1, psiwork_f(ist_f), 1)
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

            call dcopy(lzd%llr(ilr)%wfd%nvctr_c,psi_c(i_c+1),1,psi(i_tot+1),1)

            i_c = i_c + lzd%llr(ilr)%wfd%nvctr_c
            i_tot = i_tot + lzd%llr(ilr)%wfd%nvctr_c
            
            call dcopy(7*lzd%llr(ilr)%wfd%nvctr_f,psi_f(i_f+1),1,psi(i_tot+1),1)


            i_f = i_f + 7*lzd%llr(ilr)%wfd%nvctr_f
            i_tot = i_tot + 7*lzd%llr(ilr)%wfd%nvctr_f


        end do
    !!$omp end parallel 

    else
        call dcopy(collcom%ndimpsi_c, psi_c, 1, psi, 1)
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
  use module_interfaces, except_this_one => transpose_localized
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
  use module_interfaces, except_this_one => untranspose_localized
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

subroutine check_communications_locreg(iproc,nproc,orbs,Lzd,collcom,npsidim_orbs,npsidim_comp)
   use module_base
   use module_types
   use module_interfaces
   use yaml_output
   implicit none
   integer, intent(in) :: iproc,nproc
   type(orbitals_data), intent(in) :: orbs
   type(local_zone_descriptors), intent(in) :: lzd
   type(collective_comms), intent(in) :: collcom
   integer, intent(in) :: npsidim_orbs, npsidim_comp
   !local variables
   character(len=*), parameter :: subname='check_communications'
   integer, parameter :: ilog=6
   integer :: i,ispinor,iorb,indspin,i_stat,i_all,ikptsp
   integer :: ikpt,ierr,i0,ifine,ii,iiorb,ipt,jorb
   integer :: icomp
   !!$integer :: ipsi,ipsic,ipsif,ipsiworkc,ipsiworkf,jcomp,jkpt
   real(wp) :: psival,maxdiff,tt
   real(wp), dimension(:), allocatable :: psi,psit_c,psit_f
   real(wp), dimension(:,:), allocatable :: checksum
   real(wp) :: epsilon,tol
   logical :: abort

   !allocate the "wavefunction" and fill it, and also the workspace
   allocate(psi(max(npsidim_orbs,npsidim_comp)+ndebug),stat=i_stat)
   call memocc(i_stat,psi,'psi',subname)
   allocate(psit_c(sum(collcom%nrecvcounts_c)), stat=i_stat)
   call memocc(i_stat, psit_c, 'psit_c', subname)
   allocate(psit_f(7*sum(collcom%nrecvcounts_f)), stat=i_stat)
   call memocc(i_stat, psit_f, 'psit_f', subname)
   allocate(checksum(orbs%norb*orbs%nspinor,2), stat=i_stat)
   call memocc(i_stat, checksum, 'checksum', subname)
   if (orbs%norbp>0) then
      tol=1.e-10*real(npsidim_orbs,wp)/real(orbs%norbp,wp)
   else
      tol=0.0_wp
   end if

   checksum(:,:)=0.0_wp
   do iorb=1,orbs%norbp
      ikpt=(orbs%isorb+iorb-1)/orbs%norb+1
      do ispinor=1,orbs%nspinor
         indspin=(ispinor-1)*nvctr_orb(iorb)
         checksum(orbs%isorb+iorb+(ispinor-1)*orbs%nspinor,1)=0.0_wp
         do i=1,nvctr_orb(iorb)
            !vali=real(i,wp)/512.0_wp  ! *1.d-5
            call test_value_locreg(ikpt,orbs%isorb+iorb-(ikpt-1)*orbs%norb,ispinor,i,psival)
            psi(i+indspin+ind_orb(iorb))=psival!(valorb+vali)*(-1)**(ispinor-1)
            checksum(orbs%isorb+iorb+(ispinor-1)*orbs%nspinor,1)=&
                 checksum(orbs%isorb+iorb+(ispinor-1)*orbs%nspinor,1)+psival
         end do
      end do
   end do

   call transpose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, psi, psit_c, psit_f, lzd)
   
   !check the results of the transposed wavefunction
   maxdiff=0.0_wp
   if (iproc==0) call yaml_map('Number of coarse and fine DoF (MasterMPI task)',&
        (/collcom%nptsp_c,collcom%nptsp_f/),fmt='(i8)')

   do ikptsp=1,1!orbs%nkptsp !should be one for the moment
      ikpt=orbs%iskpts+ikptsp!orbs%ikptsp(ikptsp)
      ispinor=1 !for the (long?) moment
      icomp=1
      if (collcom%nptsp_c>0) then
         do ipt=1,collcom%nptsp_c 
            ii=collcom%norb_per_gridpoint_c(ipt)
            i0 = collcom%isptsp_c(ipt) 
            do i=1,ii
               iiorb=collcom%indexrecvorbital_c(i0+i)
!!$               !here a function which determin the address after mpi_alltoall
!!$               !procedure should be called
!!$               ipsitworkc=collcom%iexpand_c(icomp)
!!$               !ipsiglob=collcom%nrecvdspls_c(iproc)+1+(ipsitworkc-1)*sum(
!!$               ipsic=collcom%isendbuf_c(ipsiworkc)
!!$               ipsi=ipsic
!!$               do jorb=1,iiorb-1
!!$                  ipsi=ipsi-nvctr_c_orb(jorb)
!!$               end do
!!$               call test_value_locreg(ikpt,iiorb-(ikpt-1)*orbs%norb,ispinor,&
!!$                    ipsi,psival)
!!$               indspin=(ispinor-1)*nvctr_orb(iiorb)
!!$               maxdiff=max(abs(psit_c(i0+i)-psival),maxdiff)
               checksum(iiorb,2)=checksum(iiorb,2)+psit_c(i0+i)
               icomp=icomp+1
            end do
         end do
      end if
      icomp=1
      if (collcom%nptsp_f>0) then
         do ipt=1,collcom%nptsp_f 
            ii=collcom%norb_per_gridpoint_f(ipt) 
            i0 = collcom%isptsp_f(ipt) 
            do i=1,ii
               iiorb=collcom%indexrecvorbital_f(i0+i)
!!$               ipsitworkf=collcom%iexpand_f(icomp)
!!$               ipsif=collcom%isendbuf_f(ipsiworkf)
!!$               ipsi=ipsif
!!$               do jorb=1,iiorb-1
!!$                  ipsi=ipsi-nvctr_f_orb(jorb)
!!$               end do

               do ifine=1,7
!!$                  call test_value_locreg(ikpt,iiorb-(ikpt-1)*orbs%norb,ispinor,&
!!$                       nvctr_c_orb(iiorb)+7*(ipsi-1)+ifine,psival) 
!!$                  tt=abs(psit_f(7*(i0+i-1)+ifine)-psival)
!!$                  if (tt > maxdiff) then
!!$                     maxdiff=tt
!!$                     !call wrong_components(psival,jkpt,jorb,jcomp)
!!$                  end if
                  checksum(iiorb,2)=checksum(iiorb,2)+psit_f(7*(i0+i-1)+ifine)
               end do
               icomp=icomp+1
            end do
         end do
      end if
   end do
!!$
   if (iproc==0) call yaml_map('Tolerances for this check',&
        (/tol,real(orbs%norb,wp)*epsilon(1.0_wp)/),fmt='(1pe25.17)')

   call MPI_BARRIER(bigdft_mpi%mpi_comm, ierr)
   call mpiallred(checksum(1,1),2*orbs%norb*orbs%nspinor,MPI_SUM,&
        bigdft_mpi%mpi_comm,ierr)
   if (iproc==0) then
      maxdiff=0.0_wp
      do jorb=1,orbs%norb*orbs%nspinor
         tt=abs(checksum(jorb,1)-checksum(jorb,2))
         if (tt > maxdiff) then
            maxdiff=tt
            if (maxdiff > tol) then 
               call yaml_warning('ERROR of checksum for orbital'//trim(yaml_toa(jorb))//&
                    ': difference of '//trim(yaml_toa(tt,fmt='(1pe12.5)')))
            end if
         end if
      end do
   end if
   if (iproc==0) call yaml_map('Maxdiff for transpose (checksum)',&
        maxdiff,fmt='(1pe25.17)')

   abort = .false.
   if (abs(maxdiff) >tol) then
      call yaml_comment('ERROR (Transposition): process'//trim(yaml_toa(iproc))//&
           ' found an error of:'//trim(yaml_toa(maxdiff,fmt='(1pe15.7)')))
      !call yaml_map('Some wrong results in',(/jkpt,jorb,jcomp/),fmt='(i8)')
      abort=.true.
   end if

   if (abort) call MPI_ABORT(bigdft_mpi%mpi_comm,ierr)


   call untranspose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, psit_c, psit_f, psi, lzd)

   maxdiff=0.0_wp
   do iorb=1,orbs%norbp
      ikpt=(orbs%isorb+iorb-1)/orbs%norb+1
      do ispinor=1,orbs%nspinor
         indspin=(ispinor-1)*nvctr_orb(iorb)
         do i=1,nvctr_orb(iorb)
            call test_value_locreg(ikpt,orbs%isorb+iorb-(ikpt-1)*orbs%norb,ispinor,i,psival)
            maxdiff=max(abs(psi(i+indspin+ind_orb(iorb))-psival),maxdiff)
         end do
      end do
   end do


   abort = .false.
   if (abs(maxdiff) > real(orbs%norb,wp)*epsilon(1.0_wp)) then
      call yaml_comment('ERROR (Inverse Transposition): process'//trim(yaml_toa(iproc))//&
           ' found an error of:'//trim(yaml_toa(maxdiff,fmt='(1pe15.7)')))
      abort = .true.
   end if

   if (abort) call MPI_ABORT(bigdft_mpi%mpi_comm,ierr)
   call MPI_BARRIER(bigdft_mpi%mpi_comm, ierr)

   call mpiallred(maxdiff,1,MPI_MAX,bigdft_mpi%mpi_comm,ierr)
   if (iproc==0) call yaml_map('Maxdiff for untranspose',maxdiff,fmt='(1pe25.17)')

   i_all=-product(shape(psi))*kind(psi)
   deallocate(psi,stat=i_stat)
   call memocc(i_stat,i_all,'psi',subname)
   i_all=-product(shape(psit_c))*kind(psit_c)
   deallocate(psit_c, stat=i_stat)
   call memocc(i_stat, i_all, 'psit_c', subname)
   i_all=-product(shape(psit_f))*kind(psit_f)
   deallocate(psit_f, stat=i_stat)
   call memocc(i_stat, i_all, 'psit_f', subname)
   i_all=-product(shape(checksum))*kind(checksum)
   deallocate(checksum, stat=i_stat)
   call memocc(i_stat, i_all, 'checksum', subname)



 contains
   
   function ind_orb(iorb)
     implicit none
     integer, intent(in) :: iorb
     integer :: ind_orb
     !local variables
     integer :: jorb
     ind_orb=0
     do jorb=1,iorb-1
        ind_orb=ind_orb+nvctr_orb(jorb)
     end do
   end function ind_orb

   function nvctr_orb(iorb)
     implicit none
     integer, intent(in) :: iorb
     integer :: nvctr_orb
     !local variables
     integer :: jlr

     jlr = orbs%inwhichlocreg(iorb+orbs%isorb)
     nvctr_orb=(Lzd%Llr(jlr)%wfd%nvctr_c+7*Lzd%Llr(jlr)%wfd%nvctr_f)
     
   end function nvctr_orb

   function nvctr_c_orb(iorb)
     implicit none
     integer, intent(in) :: iorb
     integer :: nvctr_c_orb
     !local variables
     integer :: jlr

     jlr = orbs%inwhichlocreg(iorb+orbs%isorb)
     nvctr_c_orb=Lzd%Llr(jlr)%wfd%nvctr_c
     
   end function nvctr_c_orb

   function nvctr_f_orb(iorb)
     implicit none
     integer, intent(in) :: iorb
     integer :: nvctr_f_orb
     !local variables
     integer :: jlr

     jlr = orbs%inwhichlocreg(iorb+orbs%isorb)
     nvctr_f_orb=Lzd%Llr(jlr)%wfd%nvctr_f
     
   end function nvctr_f_orb


   !> define a value for the wavefunction which is dependent of the indices
   subroutine test_value_locreg(ikpt,iorb,ispinor,icomp,val)
     use module_base
     implicit none
     integer, intent(in) :: ikpt,icomp,iorb,ispinor
     real(wp), intent(out) :: val
     !local variables
     real(wp) :: valkpt,valorb,vali

     ! recognizable pattern, for debugging
     valkpt=real(10**ilog*(ikpt-1),wp)!real(512*ikpt,wp)
     valorb=real(iorb,wp)+valkpt
     vali=real(icomp,wp)*10.0_wp**(-ilog)  !real(icomp,wp)/512.0_wp  ! *1.d-5
     val=(valorb+vali)*(-1)**(ispinor-1)

   END SUBROUTINE test_value_locreg

   !>determine the components which were not communicated correctly
   !! works only with the recognizable pattern of test function
   subroutine wrong_components_locreg(psival,ikpt,iorb,icomp)
     use module_base
     implicit none
     real(wp), intent(in) :: psival
     integer, intent(out) :: ikpt,iorb,icomp

     icomp=nint((psival-real(floor(psival),wp))*10.0_wp**ilog)
     ikpt=floor(psival)/(10**ilog)
     iorb=floor(psival)-(ikpt-1)*(10**ilog)

   end subroutine wrong_components_locreg


 END SUBROUTINE check_communications_locreg





subroutine calculate_overlap_transposed(iproc, nproc, orbs, collcom, &
           psit_c1, psit_c2, psit_f1, psit_f2, ovrlp)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(orbitals_data),intent(in) :: orbs
  type(collective_comms),intent(in) :: collcom
  real(kind=8),dimension(collcom%ndimind_c),intent(in) :: psit_c1, psit_c2
  real(kind=8),dimension(7*collcom%ndimind_f),intent(in) :: psit_f1, psit_f2
  type(sparseMatrix),intent(inout) :: ovrlp

  ! Local variables
  integer :: i0, ipt, ii, iiorb, j, jjorb, i, ierr, istat, m, tid, norb, nthreads
  integer :: istart, iend, orb_rest, ind0, ind1, ind2, ind3
  integer,dimension(:),allocatable :: n
  !$ integer  :: omp_get_thread_num,omp_get_max_threads

  call timing(iproc,'ovrlptransComp','ON') !lr408t

  call to_zero(ovrlp%nvctr, ovrlp%matrix_compr(1))

  nthreads=1
  !$  nthreads = OMP_GET_max_threads()

  allocate(n(nthreads),stat=istat)

  norb = orbs%norb / nthreads
  
  orb_rest = orbs%norb-norb*nthreads

  n = norb

  do i = 1,orb_rest
     n(i) = n(i) + 1
  end do

  !$omp parallel default(private) &
  !$omp shared(collcom, ovrlp, psit_c1, psit_c2, psit_f1, psit_f2, n)
  tid=0
  !$ tid = OMP_GET_THREAD_NUM()
  istart = 1
  iend = 0

  do i = 0,tid-1
     istart = istart + n(i+1) 
  end do

  do i = 0,tid
     iend = iend + n(i+1)
  end do

 ! write(*,*) '#########################################',tid,istart,iend

  if (collcom%nptsp_c>0) then
  
      do ipt=1,collcom%nptsp_c 
          ii=collcom%norb_per_gridpoint_c(ipt) 
          i0 = collcom%isptsp_c(ipt)

          do i=1,ii
              iiorb=collcom%indexrecvorbital_c(i0+i)
              if(iiorb < istart .or. iiorb > iend) cycle

              m=mod(ii,4)
              if(m/=0) then
                  do j=1,m
                      jjorb=collcom%indexrecvorbital_c(i0+j)
                      ind0 = ovrlp%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
                      !ind0 = collcom%matrixindex_in_compressed(jjorb,iiorb)
                      ovrlp%matrix_compr(ind0) = ovrlp%matrix_compr(ind0) + psit_c1(i0+i)*psit_c2(i0+j)
                  end do
              end if
              do j=m+1,ii,4

                  jjorb=collcom%indexrecvorbital_c(i0+j+0)
                  ind0 = ovrlp%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
                  !ind0 = collcom%matrixindex_in_compressed(jjorb,iiorb)
                  ovrlp%matrix_compr(ind0) = ovrlp%matrix_compr(ind0) + psit_c1(i0+i)*psit_c2(i0+j+0)

                  jjorb=collcom%indexrecvorbital_c(i0+j+1)
                  ind1 = ovrlp%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
                  !ind1 = collcom%matrixindex_in_compressed(jjorb,iiorb)
                  ovrlp%matrix_compr(ind1) = ovrlp%matrix_compr(ind1) + psit_c1(i0+i)*psit_c2(i0+j+1)

                  jjorb=collcom%indexrecvorbital_c(i0+j+2)
                  ind2 = ovrlp%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
                  !ind2 = collcom%matrixindex_in_compressed(jjorb,iiorb)
                  ovrlp%matrix_compr(ind2) = ovrlp%matrix_compr(ind2) + psit_c1(i0+i)*psit_c2(i0+j+2)

                  jjorb=collcom%indexrecvorbital_c(i0+j+3)
                  ind3 = ovrlp%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
                  !ind3 = collcom%matrixindex_in_compressed(jjorb,iiorb)
                  ovrlp%matrix_compr(ind3) = ovrlp%matrix_compr(ind3) + psit_c1(i0+i)*psit_c2(i0+j+3)

              end do
          end do
      end do

  end if
  

  if (collcom%nptsp_f>0) then
   
      do ipt=1,collcom%nptsp_f 
          ii=collcom%norb_per_gridpoint_f(ipt) 
          i0 = collcom%isptsp_f(ipt)

          do i=1,ii

              iiorb=collcom%indexrecvorbital_f(i0+i)
              if(iiorb < istart .or. iiorb > iend) cycle

              do j=1,ii
                  jjorb=collcom%indexrecvorbital_f(i0+j)
                  ind0 = ovrlp%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
                  !ind0 = collcom%matrixindex_in_compressed(jjorb,iiorb)
                  ovrlp%matrix_compr(ind0) = ovrlp%matrix_compr(ind0) + psit_f1(7*(i0+i)-6)*psit_f2(7*(i0+j)-6)
                  ovrlp%matrix_compr(ind0) = ovrlp%matrix_compr(ind0) + psit_f1(7*(i0+i)-5)*psit_f2(7*(i0+j)-5)
                  ovrlp%matrix_compr(ind0) = ovrlp%matrix_compr(ind0) + psit_f1(7*(i0+i)-4)*psit_f2(7*(i0+j)-4)
                  ovrlp%matrix_compr(ind0) = ovrlp%matrix_compr(ind0) + psit_f1(7*(i0+i)-3)*psit_f2(7*(i0+j)-3)
                  ovrlp%matrix_compr(ind0) = ovrlp%matrix_compr(ind0) + psit_f1(7*(i0+i)-2)*psit_f2(7*(i0+j)-2)
                  ovrlp%matrix_compr(ind0) = ovrlp%matrix_compr(ind0) + psit_f1(7*(i0+i)-1)*psit_f2(7*(i0+j)-1)
                  ovrlp%matrix_compr(ind0) = ovrlp%matrix_compr(ind0) + psit_f1(7*(i0+i)-0)*psit_f2(7*(i0+j)-0)
              end do
          end do
      end do
  end if

  !$omp end parallel

  call timing(iproc,'ovrlptransComp','OF') !lr408t

  call timing(iproc,'ovrlptransComm','ON') !lr408t

  if(nproc>1) then
      call mpiallred(ovrlp%matrix_compr(1), ovrlp%nvctr, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  end if

  call timing(iproc,'ovrlptransComm','OF') !lr408t

end subroutine calculate_overlap_transposed

subroutine calculate_pulay_overlap(iproc, nproc, orbs1, orbs2, collcom1, collcom2, psit_c1, psit_c2, psit_f1, psit_f2, ovrlp)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(orbitals_data),intent(in) :: orbs1, orbs2
  type(collective_comms),intent(in) :: collcom1, collcom2
  real(kind=8),dimension(collcom1%ndimind_c),intent(in) :: psit_c1
  real(kind=8),dimension(collcom2%ndimind_c),intent(in) :: psit_c2
  real(kind=8),dimension(7*collcom1%ndimind_f),intent(in) :: psit_f1
  real(kind=8),dimension(7*collcom2%ndimind_f),intent(in) :: psit_f2
  real(kind=8),dimension(orbs1%norb,orbs2%norb),intent(out) :: ovrlp
  
  ! Local variables
  integer :: i0, j0, ipt, ii, iiorb, j, jj, jjorb, i, ierr  

  call timing(iproc,'ovrlptransComp','ON') !lr408t
  call to_zero(orbs1%norb*orbs2%norb, ovrlp(1,1))
  if(collcom1%nptsp_c/=collcom2%nptsp_c) then
      write(*,'(a,i0,a)') 'ERROR on process ',iproc,': collcom1%nptsp_c/=collcom2%nptsp_c'
      stop
  end if
  if(collcom1%nptsp_f/=collcom2%nptsp_f) then
      write(*,'(a,i0,a)') 'ERROR on process ',iproc,': collcom1%nptsp_f/=collcom2%nptsp_f'
      stop
  end if

  i0=0
  j0=0
  do ipt=1,collcom1%nptsp_c 
      ii=collcom1%norb_per_gridpoint_c(ipt)
      jj=collcom2%norb_per_gridpoint_c(ipt)
      do i=1,ii
          iiorb=collcom1%indexrecvorbital_c(i0+i)
          do j=1,jj
              jjorb=collcom2%indexrecvorbital_c(j0+j)
              ovrlp(iiorb,jjorb)=ovrlp(iiorb,jjorb)+psit_c1(i0+i)*psit_c2(j0+j)
          end do
      end do
      i0=i0+ii
      j0=j0+jj
  end do

  i0=0
  j0=0
  do ipt=1,collcom1%nptsp_f 
      ii=collcom1%norb_per_gridpoint_f(ipt)
      jj=collcom2%norb_per_gridpoint_f(ipt)
      do i=1,ii
          iiorb=collcom1%indexrecvorbital_f(i0+i)
          do j=1,jj
              jjorb=collcom2%indexrecvorbital_f(j0+j)
              ovrlp(iiorb,jjorb)=ovrlp(iiorb,jjorb)+psit_f1(7*(i0+i)-6)*psit_f2(7*(j0+j)-6)
              ovrlp(iiorb,jjorb)=ovrlp(iiorb,jjorb)+psit_f1(7*(i0+i)-5)*psit_f2(7*(j0+j)-5)
              ovrlp(iiorb,jjorb)=ovrlp(iiorb,jjorb)+psit_f1(7*(i0+i)-4)*psit_f2(7*(j0+j)-4)
              ovrlp(iiorb,jjorb)=ovrlp(iiorb,jjorb)+psit_f1(7*(i0+i)-3)*psit_f2(7*(j0+j)-3)
              ovrlp(iiorb,jjorb)=ovrlp(iiorb,jjorb)+psit_f1(7*(i0+i)-2)*psit_f2(7*(j0+j)-2)
              ovrlp(iiorb,jjorb)=ovrlp(iiorb,jjorb)+psit_f1(7*(i0+i)-1)*psit_f2(7*(j0+j)-1)
              ovrlp(iiorb,jjorb)=ovrlp(iiorb,jjorb)+psit_f1(7*(i0+i)-0)*psit_f2(7*(j0+j)-0)
          end do
      end do
      i0=i0+ii
      j0=j0+jj
  end do

  call timing(iproc,'ovrlptransComp','OF') !lr408t

  call timing(iproc,'ovrlptransComm','ON') !lr408t

  if(nproc>1) then
      call mpiallred(ovrlp(1,1), orbs1%norb*orbs2%norb, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  end if
  call timing(iproc,'ovrlptransComm','OF') !lr408t
end subroutine calculate_pulay_overlap

subroutine build_linear_combination_transposed(collcom, sparsemat, psitwork_c, psitwork_f, &
     reset, psit_c, psit_f, iproc)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  type(sparseMatrix),intent(in) :: sparsemat
  type(collective_comms),intent(in) :: collcom
  real(kind=8),dimension(collcom%ndimind_c),intent(in) :: psitwork_c
  real(kind=8),dimension(7*collcom%ndimind_f),intent(in) :: psitwork_f
  logical,intent(in) :: reset
  real(kind=8),dimension(collcom%ndimind_c),intent(inout) :: psit_c
  real(kind=8),dimension(7*collcom%ndimind_f),intent(inout) :: psit_f
  integer, intent(in) :: iproc
  ! Local variables
  integer :: i0, ipt, ii, j, iiorb, jjorb, i, m, ind0, ind1, ind2, ind3

  call timing(iproc,'lincombtrans  ','ON') !lr408t
  if(reset) then
      if(collcom%ndimind_c>0) call to_zero(collcom%ndimind_c, psit_c(1))
      if(collcom%ndimind_f>0) call to_zero(7*collcom%ndimind_f, psit_f(1))
  end if

 
  !$omp parallel default(private) &
  !$omp shared(collcom, psit_c, psitwork_c, psit_f, psitwork_f, sparsemat)

  !!write(*,'(a,i4,4i8)') 'iproc, lbound, ubound, minval, maxval',&
  !!iproc, lbound(sparsemat%matrixindex_in_compressed_fortransposed,2),&
  !!ubound(sparsemat%matrixindex_in_compressed_fortransposed,2),&
  !!minval(collcom%indexrecvorbital_c),maxval(collcom%indexrecvorbital_c)

  !$omp do
   do ipt=1,collcom%nptsp_c 
      ii=collcom%norb_per_gridpoint_c(ipt) 
      i0 = collcom%isptsp_c(ipt)
      do i=1,ii
          iiorb=collcom%indexrecvorbital_c(i0+i)
          m=mod(ii,4)
          if(m/=0) then
              do j=1,m
                  jjorb=collcom%indexrecvorbital_c(i0+j)
                  ind0 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
                  !ind0 = collcom%matrixindex_in_compressed(jjorb,iiorb)
                  !write(41,*) jjorb, iiorb, sparsemat%matrixindex_in_compressed(jjorb,iiorb)
                  !write(42,*) jjorb, iiorb, collcom%matrixindex_in_compressed(jjorb,iiorb)
                  psit_c(i0+i)=psit_c(i0+i)+sparsemat%matrix_compr(ind0)*psitwork_c(i0+j)
              end do
          end if
          do j=m+1,ii,4
              jjorb=collcom%indexrecvorbital_c(i0+j+0)
              ind0 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
              !ind0 = collcom%matrixindex_in_compressed(jjorb,iiorb)
              psit_c(i0+i)=psit_c(i0+i)+sparsemat%matrix_compr(ind0)*psitwork_c(i0+j+0)

              jjorb=collcom%indexrecvorbital_c(i0+j+1)
              ind1 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
              !ind1 = collcom%matrixindex_in_compressed(jjorb,iiorb)
              psit_c(i0+i)=psit_c(i0+i)+sparsemat%matrix_compr(ind1)*psitwork_c(i0+j+1)

              jjorb=collcom%indexrecvorbital_c(i0+j+2)
              ind2 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
              !ind2 = collcom%matrixindex_in_compressed(jjorb,iiorb)
              psit_c(i0+i)=psit_c(i0+i)+sparsemat%matrix_compr(ind2)*psitwork_c(i0+j+2)

              jjorb=collcom%indexrecvorbital_c(i0+j+3)
              ind3 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
              !ind3 = collcom%matrixindex_in_compressed(jjorb,iiorb)
              psit_c(i0+i)=psit_c(i0+i)+sparsemat%matrix_compr(ind3)*psitwork_c(i0+j+3)

          end do
      end do
  end do
  !$omp end do

  !$omp do
  do ipt=1,collcom%nptsp_f 
      ii=collcom%norb_per_gridpoint_f(ipt) 
       i0 = collcom%isptsp_f(ipt)
      do i=1,ii
          iiorb=collcom%indexrecvorbital_f(i0+i)
          do j=1,ii
              jjorb=collcom%indexrecvorbital_f(i0+j)
              ind0 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
              !ind0 = collcom%matrixindex_in_compressed(jjorb,iiorb)
              psit_f(7*(i0+i)-6) = psit_f(7*(i0+i)-6) + sparsemat%matrix_compr(ind0)*psitwork_f(7*(i0+j)-6)
              psit_f(7*(i0+i)-5) = psit_f(7*(i0+i)-5) + sparsemat%matrix_compr(ind0)*psitwork_f(7*(i0+j)-5)
              psit_f(7*(i0+i)-4) = psit_f(7*(i0+i)-4) + sparsemat%matrix_compr(ind0)*psitwork_f(7*(i0+j)-4)
              psit_f(7*(i0+i)-3) = psit_f(7*(i0+i)-3) + sparsemat%matrix_compr(ind0)*psitwork_f(7*(i0+j)-3)
              psit_f(7*(i0+i)-2) = psit_f(7*(i0+i)-2) + sparsemat%matrix_compr(ind0)*psitwork_f(7*(i0+j)-2)
              psit_f(7*(i0+i)-1) = psit_f(7*(i0+i)-1) + sparsemat%matrix_compr(ind0)*psitwork_f(7*(i0+j)-1)
              psit_f(7*(i0+i)-0) = psit_f(7*(i0+i)-0) + sparsemat%matrix_compr(ind0)*psitwork_f(7*(i0+j)-0)
          end do
      end do
     
  end do
  !$omp end do
  !$omp end parallel

  call timing(iproc,'lincombtrans  ','OF') !lr408t

end subroutine build_linear_combination_transposed




subroutine check_grid_point_from_boxes(i1, i2, i3, lr, overlap_possible)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: i1, i2, i3
  type(locreg_descriptors),intent(in) :: lr  
  logical,intent(out) :: overlap_possible

  ! Local variables
  logical :: ovrlpx, ovrlpy, ovrlpz
  
  ovrlpx = (i1>=lr%ns1 .and. i1<=lr%ns1+lr%d%n1)
  ovrlpy = (i2>=lr%ns2 .and. i2<=lr%ns2+lr%d%n2)
  ovrlpz = (i3>=lr%ns3 .and. i3<=lr%ns3+lr%d%n3)
  if(ovrlpx .and. ovrlpy .and. ovrlpz) then
      overlap_possible=.true.
  else
      overlap_possible=.true.
  end if

end subroutine check_grid_point_from_boxes


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


subroutine normalize_transposed(iproc, nproc, orbs, collcom, psit_c, psit_f, norm)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(orbitals_data),intent(in):: orbs
  type(collective_comms),intent(in):: collcom
  real(8),dimension(collcom%ndimind_c),intent(inout):: psit_c
  real(8),dimension(7*collcom%ndimind_f),intent(inout):: psit_f
  real(8),dimension(orbs%norb),intent(out):: norm
  
  ! Local variables
  integer:: i0, ipt, ii, iiorb, i, ierr, iorb

  call timing(iproc,'norm_trans','ON')

  call to_zero(orbs%norb, norm(1))

  !$omp parallel default(private) &
  !$omp shared(collcom, norm, psit_c,psit_f,orbs)

  if (collcom%nptsp_c>0) then
      !$omp do reduction(+:norm)
      do ipt=1,collcom%nptsp_c 
          ii=collcom%norb_per_gridpoint_c(ipt)
          i0 = collcom%isptsp_c(ipt) 
          do i=1,ii
              iiorb=collcom%indexrecvorbital_c(i0+i)
              norm(iiorb)=norm(iiorb)+psit_c(i0+i)**2
          end do
      end do
      !$omp end do
  end if

  if (collcom%nptsp_f>0) then
      !$omp do reduction(+:norm)
      do ipt=1,collcom%nptsp_f 
          ii=collcom%norb_per_gridpoint_f(ipt) 
          i0 = collcom%isptsp_f(ipt) 
          do i=1,ii
              iiorb=collcom%indexrecvorbital_f(i0+i)
              norm(iiorb)=norm(iiorb)+psit_f(7*(i0+i)-6)**2
              norm(iiorb)=norm(iiorb)+psit_f(7*(i0+i)-5)**2
              norm(iiorb)=norm(iiorb)+psit_f(7*(i0+i)-4)**2
              norm(iiorb)=norm(iiorb)+psit_f(7*(i0+i)-3)**2
              norm(iiorb)=norm(iiorb)+psit_f(7*(i0+i)-2)**2
              norm(iiorb)=norm(iiorb)+psit_f(7*(i0+i)-1)**2
              norm(iiorb)=norm(iiorb)+psit_f(7*(i0+i)-0)**2
          end do
      end do
      !$omp end do
  end if

  !$omp end parallel
  
  if(nproc>1) then
      call mpiallred(norm(1), orbs%norb, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  end if

  do iorb=1,orbs%norb
     norm(iorb)=1.d0/sqrt(norm(iorb))
  end do

  !$omp parallel default(private) shared(norm,orbs,collcom,psit_c,psit_f)  
  !$omp do
  do ipt=1,collcom%nptsp_c 
      ii=collcom%norb_per_gridpoint_c(ipt)
      i0=collcom%isptsp_c(ipt)
      do i=1,ii
          iiorb=collcom%indexrecvorbital_c(i0+i)
          psit_c(i0+i)=psit_c(i0+i)*norm(iiorb)
      end do
    
  end do
  !$omp end do
  !$omp do
  do ipt=1,collcom%nptsp_f 
      ii=collcom%norb_per_gridpoint_f(ipt)
      i0 = collcom%isptsp_f(ipt) 
      do i=1,ii
          iiorb=collcom%indexrecvorbital_f(i0+i)
          psit_f(7*(i0+i)-6)=psit_f(7*(i0+i)-6)*norm(iiorb)
          psit_f(7*(i0+i)-5)=psit_f(7*(i0+i)-5)*norm(iiorb)
          psit_f(7*(i0+i)-4)=psit_f(7*(i0+i)-4)*norm(iiorb)
          psit_f(7*(i0+i)-3)=psit_f(7*(i0+i)-3)*norm(iiorb)
          psit_f(7*(i0+i)-2)=psit_f(7*(i0+i)-2)*norm(iiorb)
          psit_f(7*(i0+i)-1)=psit_f(7*(i0+i)-1)*norm(iiorb)
          psit_f(7*(i0+i)-0)=psit_f(7*(i0+i)-0)*norm(iiorb)
      end do
  
  end do
!$omp end do
!$omp end parallel

  call timing(iproc,'norm_trans','OF')

end subroutine normalize_transposed



subroutine init_matrixindex_in_compressed_fortransposed(iproc, nproc, orbs, collcom, collcom_shamop, &
           collcom_sr, sparsemat)
  use module_base
  use module_types
  use module_interfaces, except_this_one => init_matrixindex_in_compressed_fortransposed
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(orbitals_data),intent(in) :: orbs
  type(collective_comms),intent(in) :: collcom, collcom_shamop, collcom_sr
  type(sparseMatrix), intent(inout) :: sparsemat
  
  ! Local variables
  integer :: iorb, jorb, istat, imin, imax
  !integer :: kproc,jproc,jjorbold,jjorb,isend,irecv,ilr,ijorb,iiorb,ind,ierr, irow, irowold, iseg
  integer :: compressed_index
!  integer,dimension(:,:),allocatable :: sendbuf, requests, iminmaxarr
  character(len=*),parameter :: subname='initSparseMatrix'


  ! for the calculation of overlaps and the charge density
  imin=minval(collcom%indexrecvorbital_c)
  imin=min(imin,minval(collcom%indexrecvorbital_f))
  imin=min(imin,minval(collcom_shamop%indexrecvorbital_c))
  imin=min(imin,minval(collcom_shamop%indexrecvorbital_f))
  imin=min(imin,minval(collcom_sr%indexrecvorbital_c))
  imax=maxval(collcom%indexrecvorbital_c)
  imax=max(imax,maxval(collcom%indexrecvorbital_f))
  imax=max(imax,maxval(collcom_shamop%indexrecvorbital_c))
  imax=max(imax,maxval(collcom_shamop%indexrecvorbital_f))
  imax=max(imax,maxval(collcom_sr%indexrecvorbital_c))

  allocate(sparsemat%matrixindex_in_compressed_fortransposed(imin:imax,imin:imax), stat=istat)
  call memocc(istat, sparsemat%matrixindex_in_compressed_fortransposed, &
      'sparsemat%matrixindex_in_compressed_fortransposed', subname)

!!  allocate(sendbuf(orbs%norbp,orbs%norb), stat=istat)
!!  call memocc(istat, sendbuf, 'sendbuf', subname)

  do iorb=imin,imax
      do jorb=imin,imax
          sparsemat%matrixindex_in_compressed_fortransposed(iorb,jorb)=compressed_index(iorb,jorb,orbs%norb,sparsemat)
          !sendbuf(jorb,iorb)=compressed_index(jorb,iiorb,orbs%norb,sparsemat)
          !sendbuf(iorb,jorb)=compressed_index(iiorb,jorb,orbs%norb,sparsemat)
      end do
  end do

!!  do iorb=1,orbs%norbp
!!      iiorb=orbs%isorb+iorb
!!      do jorb=1,orbs%norb
!!          !sendbuf(jorb,iorb)=compressed_index(iiorb,jorb,orbs%norb,sparsemat)
!!          !sendbuf(jorb,iorb)=compressed_index(jorb,iiorb,orbs%norb,sparsemat)
!!          sendbuf(iorb,jorb)=compressed_index(iiorb,jorb,orbs%norb,sparsemat)
!!      end do
!!  end do
!!
!!  allocate(iminmaxarr(2,0:nproc-1), stat=istat)
!!  call memocc(istat, iminmaxarr, 'iminmaxarr', subname)
!!  call to_zero(2*nproc, iminmaxarr(1,0))
!!  iminmaxarr(1,iproc)=imin
!!  iminmaxarr(2,iproc)=imax
!!  call mpiallred(iminmaxarr(1,0), 2*nproc, mpi_sum, bigdft_mpi%mpi_comm, ierr)
!!
!!  allocate(requests(maxval(orbs%norb_par(:,0))*nproc,2), stat=istat)
!!  call memocc(istat, requests, 'requests', subname)
!!
!!  if (nproc>1) then
!!
!!      isend=0
!!      irecv=0
!!      do jproc=0,nproc-1
!!          do jorb=1,orbs%norb_par(jproc,0)
!!              jjorb=jorb+orbs%isorb_par(jproc)
!!              do kproc=0,nproc-1
!!                  if (jjorb>=iminmaxarr(1,kproc) .and. jjorb<=iminmaxarr(2,kproc)) then
!!                      ! send from jproc to kproc
!!                      if (iproc==jproc) then
!!                          isend=isend+1
!!                          call mpi_isend(sendbuf(1,jorb), orbs%norb, &
!!                               mpi_integer, kproc, jjorb, bigdft_mpi%mpi_comm, requests(isend,1), ierr)
!!                      end if
!!                      if (iproc==kproc) then
!!                          irecv=irecv+1
!!                          call mpi_irecv(collcom%matrixindex_in_compressed(1,jjorb), orbs%norb, &
!!                               mpi_integer, jproc, jjorb, bigdft_mpi%mpi_comm, requests(irecv,2), ierr)
!!                      end if
!!                  end if
!!              end do
!!          end do
!!      end do
!!
!!      call mpi_waitall(isend, requests(1,1), mpi_statuses_ignore, ierr)
!!      call mpi_waitall(irecv, requests(1,2), mpi_statuses_ignore, ierr)
!!
!!  else
!!      call vcopy(orbs%norb*orbs%norbp, sendbuf(1,1), 1, collcom%matrixindex_in_compressed(1,1), 1)
!!  end if

  !!do iorb=imin,imax
  !!    do jorb=1,orbs%norb
  !!        write(200+iproc,*) iorb,jorb,collcom%matrixindex_in_compressed(jorb,iorb)
  !!    end do
  !!end do


!!  iall=-product(shape(iminmaxarr))*kind(iminmaxarr)
!!  deallocate(iminmaxarr, stat=istat)
!!  call memocc(istat, iall, 'iminmaxarr', subname)
!!
!!  iall=-product(shape(requests))*kind(requests)
!!  deallocate(requests, stat=istat)
!!  call memocc(istat, iall, 'requests', subname)
!!
!!  iall=-product(shape(sendbuf))*kind(sendbuf)
!!  deallocate(sendbuf, stat=istat)
!!  call memocc(istat, iall, 'sendbuf', subname)


end subroutine init_matrixindex_in_compressed_fortransposed
