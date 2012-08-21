!> @file
! Intialization of the collective communications for the linear version
! @author
!    Copyright (C) 2011-2012 BigDFT group
!    This file is distributed under the terms of the
!    GNU General Public License, see ~/COPYING file
!    or http://www.gnu.org/copyleft/gpl.txt .
!    For the list of contributors, see ~/AUTHORS


subroutine init_collective_comms(iproc, nproc, orbs, lzd, collcom, collcom_reference)
  use module_base
  use module_types
  use module_interfaces, except_this_one => init_collective_comms
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(orbitals_data),intent(in) :: orbs
  type(local_zone_descriptors),intent(in) :: lzd
  type(collective_comms),intent(out) :: collcom
  type(collective_comms),optional,intent(in) :: collcom_reference
  
  ! Local variables
  integer :: ii, istat, iorb, iiorb, ilr, iall, istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, ierr, ipt
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


  call mpi_barrier(mpi_comm_world, ierr)
  t1=mpi_wtime()
  call get_weights(iproc, nproc, orbs, lzd, weight_c, weight_f, weight_c_tot, weight_f_tot)
  call mpi_barrier(mpi_comm_world, ierr)
  t2=mpi_wtime()
  !if(iproc==0) write(*,'(a,es10.3)') 'time for part 1:',t2-t1

  ! Assign the grid points to the processes such that the work is equally dsitributed
  allocate(istartend_c(2,0:nproc-1), stat=istat)
  call memocc(istat, istartend_c, 'istartend_c', subname)
  allocate(istartend_f(2,0:nproc-1), stat=istat)
  call memocc(istat, istartend_f, 'istartend_f', subname)
  if(.not.present(collcom_reference)) then
      call mpi_barrier(mpi_comm_world, ierr)
      t1=mpi_wtime()
      call mpi_barrier(mpi_comm_world, ierr)
      t2=mpi_wtime()
      !if(iproc==0) write(*,'(a,es10.3)') 'time for part 2:',t2-t1
      call assign_weight_to_process(iproc, nproc, lzd, weight_c, weight_f, weight_c_tot, weight_f_tot, &
           istartend_c, istartend_f, istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, &
           weightp_c, weightp_f, collcom%nptsp_c, collcom%nptsp_f)
  else
      allocate(npts_par_c(0:nproc-1), stat=istat)
       call memocc(istat, npts_par_c, 'npts_par_c', subname)
      allocate(npts_par_f(0:nproc-1), stat=istat)
       call memocc(istat, npts_par_f, 'npts_par_f', subname)
      npts_par_c=0
      npts_par_f=0
      npts_par_c(iproc)=collcom_reference%nptsp_c
      npts_par_f(iproc)=collcom_reference%nptsp_f
      call mpiallred(npts_par_c(0), nproc, mpi_sum, mpi_comm_world, ierr)
      call mpiallred(npts_par_f(0), nproc, mpi_sum, mpi_comm_world, ierr)
      call assign_weight_to_process2(iproc, nproc, lzd, weight_c, weight_f, weight_c_tot, weight_f_tot, &
           npts_par_c, npts_par_f, &
           istartend_c, istartend_f, istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, &
           weightp_c, weightp_f, collcom%nptsp_c, collcom%nptsp_f)
      iall=-product(shape(npts_par_c))*kind(npts_par_c)
      deallocate(npts_par_c, stat=istat)
      call memocc(istat, iall, 'npts_par_c', subname)
      iall=-product(shape(npts_par_f))*kind(npts_par_f)
      deallocate(npts_par_f, stat=istat)
      call memocc(istat, iall, 'npts_par_f', subname)
  end if





  ! some checks
  if(nproc>1) then
      call mpi_allreduce(weightp_c, tt, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
  else
      tt=weightp_c
  end if
  if(tt/=weight_c_tot) stop 'wrong partition of coarse weights'
  if(nproc>1) then
      call mpi_allreduce(weightp_f, tt, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
  else
      tt=weightp_f
  end if     
  if(tt/=weight_f_tot) stop 'wrong partition of fine weights'
  if(nproc>1) then
      call mpi_allreduce(collcom%nptsp_c, ii, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr)
  else
      ii=collcom%nptsp_c
  end if
  if(ii/=lzd%glr%wfd%nvctr_c) stop 'wrong partition of coarse grid points'
  if(nproc>1) then
      call mpi_allreduce(collcom%nptsp_f, ii, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr)
  else
      ii=collcom%nptsp_f
  end if
  if(ii/=lzd%glr%wfd%nvctr_f) stop 'init_collective_comms: wrong partition of fine grid points'

  ! Allocate the keys
  allocate(collcom%norb_per_gridpoint_c(collcom%nptsp_c), stat=istat)
  call memocc(istat, collcom%norb_per_gridpoint_c, 'collcom%norb_per_gridpoint_c', subname)
  allocate(collcom%norb_per_gridpoint_f(collcom%nptsp_f), stat=istat)
  call memocc(istat, collcom%norb_per_gridpoint_f, 'collcom%norb_per_gridpoint_f', subname)
  call mpi_barrier(mpi_comm_world, ierr)
  call mpi_barrier(mpi_comm_world, ierr)
  t1=mpi_wtime()
  !call determine_num_orbs_per_gridpoint(iproc, nproc, orbs, lzd, istartend_c, istartend_f, &
  !     istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, &
  !     weightp_c, weightp_f, collcom%nptsp_c, collcom%nptsp_f, &
  !     collcom%norb_per_gridpoint_c, collcom%norb_per_gridpoint_f)
  call determine_num_orbs_per_gridpoint_new(iproc, nproc, orbs, lzd, istartend_c, istartend_f, &
       istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, &
       weightp_c, weightp_f, collcom%nptsp_c, collcom%nptsp_f, weight_c, weight_f, &
       collcom%norb_per_gridpoint_c, collcom%norb_per_gridpoint_f)
  call mpi_barrier(mpi_comm_world, ierr)
  t2=mpi_wtime()
  !if(iproc==0) write(*,'(a,es10.3)') 'time for part 3:',t2-t1

  ! Determine the index of a grid point i1,i2,i3 in the compressed array
  call mpi_barrier(mpi_comm_world, ierr)
  t1=mpi_wtime()
  call get_index_in_global2(lzd%glr, index_in_global_c, index_in_global_f)
  call mpi_barrier(mpi_comm_world, ierr)
  t2=mpi_wtime()
  !if(iproc==0) write(*,'(a,es10.3)') 'time for part 4:',t2-t1




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
call mpi_barrier(mpi_comm_world, ierr)
t1=mpi_wtime()
  call determine_communication_arrays(iproc, nproc, orbs, lzd, istartend_c, istartend_f, &
       index_in_global_c, index_in_global_f, weightp_c, weightp_f, &
       collcom%nsendcounts_c, collcom%nsenddspls_c, collcom%nrecvcounts_c, collcom%nrecvdspls_c, &
       collcom%nsendcounts_f, collcom%nsenddspls_f, collcom%nrecvcounts_f, collcom%nrecvdspls_f)
call mpi_barrier(mpi_comm_world, ierr)
t2=mpi_wtime()
!if(iproc==0) write(*,'(a,es10.3)') 'time for part 5:',t2-t1


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

call mpi_barrier(mpi_comm_world, ierr)
t1=mpi_wtime()
  call get_switch_indices(iproc, nproc, orbs, lzd, collcom%ndimpsi_c, collcom%ndimpsi_f, istartend_c, istartend_f, &
       collcom%nsendcounts_c, collcom%nsenddspls_c, collcom%ndimind_c, collcom%nrecvcounts_c, collcom%nrecvdspls_c, &
       collcom%nsendcounts_f, collcom%nsenddspls_f, collcom%ndimind_f, collcom%nrecvcounts_f, collcom%nrecvdspls_f, &
       index_in_global_c, index_in_global_f, &
       weightp_c, weightp_f, collcom%isendbuf_c, collcom%irecvbuf_c, collcom%isendbuf_f, collcom%irecvbuf_f, &
       collcom%indexrecvorbital_c, collcom%iextract_c, collcom%iexpand_c, &
       collcom%indexrecvorbital_f, collcom%iextract_f, collcom%iexpand_f)
call mpi_barrier(mpi_comm_world, ierr)
t2=mpi_wtime()
!if(iproc==0) write(*,'(a,es10.3)') 'time for part 6:',t2-t1

  ! These variables are used in various subroutines to speed up the code
  allocate(collcom%isptsp_c(collcom%nptsp_c), stat=istat)
  call memocc(istat, collcom%isptsp_c, 'collcom%isptsp_c', subname)
  allocate(collcom%isptsp_f(collcom%nptsp_f), stat=istat)
  call memocc(istat, collcom%isptsp_f, 'collcom%isptsp_f', subname)
  collcom%isptsp_c(1) = 0
  do ipt=2,collcom%nptsp_c
        collcom%isptsp_c(ipt) = collcom%isptsp_c(ipt-1) + collcom%norb_per_gridpoint_c(ipt-1)
  end do
  collcom%isptsp_f(1) = 0
  do ipt=2,collcom%nptsp_f
        collcom%isptsp_f(ipt) = collcom%isptsp_f(ipt-1) + collcom%norb_per_gridpoint_f(ipt-1)
  end do


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


  ! Calculate the weights for the coarse part.
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
              weight_c(ii1,ii2,ii3)=weight_c(ii1,ii2,ii3)+1.d0
              weight_c_tot=weight_c_tot+1.d0
          end do
      end do
  
      ! Calculate the weights for the fine part.
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
          do i=i0,i1
              ii1=i+lzd%llr(ilr)%ns1
              ii2=i2+lzd%llr(ilr)%ns2
              ii3=i3+lzd%llr(ilr)%ns3
              weight_f(ii1,ii2,ii3)=weight_f(ii1,ii2,ii3)+1.d0
              weight_f_tot=weight_f_tot+1.d0
          end do
      end do
  end do


  ! Sum up among all processes.
  if(nproc>1) then
      call mpiallred(weight_c_tot, 1, mpi_sum, mpi_comm_world, ierr)
      call mpiallred(weight_f_tot, 1, mpi_sum, mpi_comm_world, ierr)
      ii=(lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(lzd%glr%d%n3+1)
      call mpiallred(weight_c(0,0,0), ii,  mpi_sum, mpi_comm_world, ierr)
      call mpiallred(weight_f(0,0,0), ii,  mpi_sum, mpi_comm_world, ierr)
  end if


end subroutine get_weights



subroutine assign_weight_to_process(iproc, nproc, lzd, weight_c, weight_f, weight_tot_c, weight_tot_f, &
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
  integer,dimension(2,0:nproc-1),intent(out) :: istartend_c, istartend_f
  integer,intent(out) :: istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f
  real(kind=8),intent(out) :: weightp_c, weightp_f
  integer,intent(out) :: nptsp_c, nptsp_f
  
  ! Local variables
  integer :: jproc, i1, i2, i3, ii, ii2, istart, iend, jj, j0, j1, jprocdone
  integer :: i, iseg, i0, iitot, ierr, iiseg
  real(kind=8) :: tt, tt2, weight_c_ideal, weight_f_ideal

  ! Ideal weight per process.
  weight_c_ideal=weight_tot_c/dble(nproc)
  weight_f_ideal=weight_tot_f/dble(nproc)


  ! First the coarse part...
  jproc=0
  tt=0.d0
  tt2=0.d0
  iitot=0
  ii2=0
  iiseg=1
  jprocdone=-1
  weightp_c=0.d0
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
       do i=i0,i1
           tt=tt+weight_c(i,i2,i3)
           iitot=iitot+1
           if((tt>=weight_c_ideal .or. iseg==lzd%glr%wfd%nseg_c) .and. i==i1 .and. jproc<=nproc-2) then
               if(tt==weight_tot_c .and. jproc==nproc-1) then
                   ! this process also has to take the remaining points, even if they have no weight
                   iitot=lzd%glr%wfd%nvctr_c-ii2
               end if
               if(iproc==jproc) then
                   weightp_c=tt
                   nptsp_c=iitot
                   istartp_seg_c=iiseg
                   iendp_seg_c=iseg
                   if(tt==weight_tot_c .and. jproc==nproc-1) then
                       ! this process also has to take the remaining segments, even if they have no weight
                       iendp_seg_c=lzd%glr%wfd%nseg_c
                   end if
               end if
               istartend_c(1,jproc)=ii2+1
               istartend_c(2,jproc)=min(istartend_c(1,jproc)+iitot-1,lzd%glr%wfd%nvctr_c)
               tt2=tt2+tt
               tt=0.d0
               ii2=ii2+iitot
               iitot=0
               jproc=jproc+1
               iiseg=iseg
               if(ii2>=lzd%glr%wfd%nvctr_c) then
                   ! everything is distributed
                   jprocdone=jproc
                   exit loop_nseg_c
               end if
           end if
       end do
   end do loop_nseg_c

  if(jprocdone>0) then
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
  else
      if(iproc==nproc-1) then
          ! Take the rest
          weightp_c=weight_tot_c-tt2
          nptsp_c=lzd%glr%wfd%nvctr_c-ii2
          istartp_seg_c=iiseg
          iendp_seg_c=lzd%glr%wfd%nseg_c
      end if
      istartend_c(1,nproc-1)=ii2+1
      istartend_c(2,nproc-1)=istartend_c(1,nproc-1)+iitot-1
  end if

  ! some check
  ii=istartend_c(2,iproc)-istartend_c(1,iproc)+1
  if(nproc>1) call mpiallred(ii, 1, mpi_sum, mpi_comm_world, ierr)
  if(ii/=lzd%glr%wfd%nvctr_c) then
     write(*,*) 'ii/=lzd%glr%wfd%nvctr_c',ii,lzd%glr%wfd%nvctr_c
     stop
  end if


  ! Now the fine part...
  jproc=0
  tt=0.d0
  tt2=0.d0
  iitot=0
  ii2=0
  weightp_f=0.d0
  istart=lzd%glr%wfd%nseg_c+min(1,lzd%glr%wfd%nseg_f)
  iend=istart+lzd%glr%wfd%nseg_f-1
  iiseg=istart
  jprocdone=-1
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
       do i=i0,i1
           tt=tt+weight_f(i,i2,i3)
           iitot=iitot+1
           if((tt>=weight_f_ideal .or. iseg==iend) .and. i==i1 .and. jproc<=nproc-2) then
               if(tt==weight_tot_f .and. jproc==nproc-1) then
                   ! this process also has to take the remaining points, even if they have no weight
                   iitot=lzd%glr%wfd%nvctr_f-ii2
               end if
               if(iproc==jproc) then
                   weightp_f=tt
                   nptsp_f=iitot
                   istartp_seg_f=iiseg
                   iendp_seg_f=iseg
                   if(tt==weight_tot_f .and. jproc==nproc-1) then
                       ! this process also has to take the remaining segments, even if they have no weight
                       iendp_seg_f=iend
                   end if
               end if
               istartend_f(1,jproc)=ii2+1
               istartend_f(2,jproc)=min(istartend_f(1,jproc)+iitot-1,lzd%glr%wfd%nvctr_f)
               tt2=tt2+tt
               tt=0.d0
               ii2=ii2+iitot
               iitot=0
               jproc=jproc+1
               iiseg=iseg
               if(ii2>=lzd%glr%wfd%nvctr_f) then
                   ! everything is distributed
                   jprocdone=jproc
                   exit loop_nseg_f
               end if
           end if
       end do
   end do loop_nseg_f
  if(jprocdone>0) then
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
  else
      if(iproc==nproc-1) then
          ! Take the rest
          weightp_f=weight_tot_f-tt2
          nptsp_f=lzd%glr%wfd%nvctr_f-ii2
          istartp_seg_f=iiseg
          iendp_seg_f=iend
      end if
      istartend_f(1,nproc-1)=ii2+1
      istartend_f(2,nproc-1)=istartend_f(1,nproc-1)+iitot-1
  end if

  ! some check
  ii=istartend_f(2,iproc)-istartend_f(1,iproc)+1
  if(nproc>1) call mpiallred(ii, 1, mpi_sum, mpi_comm_world, ierr)
  if(ii/=lzd%glr%wfd%nvctr_f) stop 'assign_weight_to_process: ii/=lzd%glr%wfd%nvctr_f'



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
  if(nproc>1) call mpiallred(ii, 1, mpi_sum, mpi_comm_world, ierr)
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
  if(nproc>1) call mpiallred(ii, 1, mpi_sum, mpi_comm_world, ierr)
  if(ii/=lzd%glr%wfd%nvctr_f) stop 'assign_weight_to_process2: ii/=lzd%glr%wfd%nvctr_f'



end subroutine assign_weight_to_process2





subroutine determine_num_orbs_per_gridpoint_new(iproc, nproc, orbs, lzd, istartend_c, istartend_f, &
           istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, &
           weightp_c, weightp_f, nptsp_c, nptsp_f, weight_c, weight_f, &
           norb_per_gridpoint_c, norb_per_gridpoint_f)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc, nptsp_c, nptsp_f, istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f
  type(orbitals_data),intent(in):: orbs
  type(local_zone_descriptors),intent(in):: lzd
  integer,dimension(2,0:nproc-1),intent(in):: istartend_c, istartend_f
  real(8),intent(in):: weightp_c, weightp_f
  real(8),dimension(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3),intent(in):: weight_c, weight_f
  integer,dimension(nptsp_c),intent(out):: norb_per_gridpoint_c
  integer,dimension(nptsp_f),intent(out):: norb_per_gridpoint_f
  
  ! Local variables
  integer:: ii, iiorb, i1, i2, i3, iipt, iorb, iii, npgp, iseg, jj, j0, j1, iitot, ilr, i, istart, iend, i0, istat, iall
  logical:: found, overlap_possible
  integer,dimension(:),allocatable:: iseg_start_c, iseg_start_f
  character(len=*),parameter:: subname='determine_num_orbs_per_gridpoint'
  real(8):: t1, t2, t1tot, t2tot, t_check_gridpoint

  allocate(iseg_start_c(lzd%nlr), stat=istat)
  call memocc(istat, iseg_start_c, 'iseg_start_c', subname)
  allocate(iseg_start_f(lzd%nlr), stat=istat)
  call memocc(istat, iseg_start_f, 'iseg_start_f', subname)

  iseg_start_c=1
  iseg_start_f=1

  iitot=0
  iiorb=0
  iipt=0
t_check_gridpoint=0.d0
t1tot=mpi_wtime()
  !write(*,*) 'iproc, istartp_seg_c,iendp_seg_c', iproc, istartp_seg_c,iendp_seg_c
    !do iseg=1,lzd%glr%wfd%nseg_c
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
           !iitot=iitot+1
           iitot=jj+i-i0
           if(iitot>=istartend_c(1,iproc) .and. iitot<=istartend_c(2,iproc)) then
               !write(200+iproc,'(5i10)') iitot, iseg, iitot, jj, jj+i-i0
               iipt=iipt+1
               npgp=0
               !do iorb=1,orbs%norb
               !    ilr=orbs%inwhichlocreg(iorb)
               !    ! Check whether this orbitals extends here
               !    call check_grid_point_from_boxes(i, i2, i3, lzd%llr(ilr), overlap_possible)
               !    if(.not. overlap_possible) then
               !        found=.false.
               !    else
               !        t1=mpi_wtime()
               !        call check_gridpoint(lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, &
               !             lzd%llr(ilr)%ns1, lzd%llr(ilr)%ns2, lzd%llr(ilr)%ns3, lzd%llr(ilr)%wfd%keygloc, &
               !             i, i2, i3, iseg_start_c(ilr), found)
               !        t2=mpi_wtime()
               !        t_check_gridpoint=t_check_gridpoint+t2-t1
               !    end if
               !    if(found) then
               !        npgp=npgp+1
               !        iiorb=iiorb+1
               !    end if
               !end do
               npgp = weight_c(i,i2,i3)
               iiorb=iiorb+npgp
               norb_per_gridpoint_c(iipt)=npgp
           end if
      end do
  end do

  if(iipt/=nptsp_c) stop 'iipt/=nptsp_c'
  if(iiorb/=nint(weightp_c)) stop 'iiorb/=weightp_c'



  iitot=0
  iiorb=0
  iipt=0
    istart=lzd%glr%wfd%nseg_c+min(1,lzd%glr%wfd%nseg_f)
    iend=istart+lzd%glr%wfd%nseg_f-1
    !do iseg=istart,iend
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
           !iitot=iitot+1
           iitot=jj+i-i0
           if(iitot>=istartend_f(1,iproc) .and. iitot<=istartend_f(2,iproc)) then
               iipt=iipt+1
               npgp=0
               !do iorb=1,orbs%norb
               !    ilr=orbs%inwhichlocreg(iorb)
               !    ! Check whether this orbitals extends here
               !    call check_grid_point_from_boxes(i, i2, i3, lzd%llr(ilr), overlap_possible)
               !    if(.not. overlap_possible) then
               !        found=.false.
               !    else
               !        iii=lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)
               !        t1=mpi_wtime()
               !        call check_gridpoint(lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, &
               !             lzd%llr(ilr)%ns1, lzd%llr(ilr)%ns2, lzd%llr(ilr)%ns3, &
               !             lzd%llr(ilr)%wfd%keygloc(1,iii), &
               !             i, i2, i3, iseg_start_f(ilr), found)
               !        t2=mpi_wtime()
               !        t_check_gridpoint=t_check_gridpoint+t2-t1
               !    end if
               !    if(found) then
               !        npgp=npgp+1
               !        iiorb=iiorb+1
               !    end if
               !end do
               npgp = weight_f(i,i2,i3)
               iiorb=iiorb+npgp
               norb_per_gridpoint_f(iipt)=npgp
           end if
      end do
  end do

  if(iipt/=nptsp_f) stop 'iipt/=nptsp_f'
  !write(*,*) 'iiorb, weightp_f', iiorb, weightp_f
  if(iiorb/=nint(weightp_f)) stop 'iiorb/=weightp_f'


  iall=-product(shape(iseg_start_c))*kind(iseg_start_c)
  deallocate(iseg_start_c, stat=istat)
  call memocc(istat, iall, 'iseg_start_c', subname)
  iall=-product(shape(iseg_start_f))*kind(iseg_start_f)
  deallocate(iseg_start_f, stat=istat)
  call memocc(istat, iall, 'iseg_start_f', subname)

t2tot=mpi_wtime()
!if(iproc==0) write(*,'(a,es14.5)') 'in sub determine_num_orbs_per_gridpoint: iproc, total time', t2tot-t1tot
!if(iproc==0) write(*,'(a,es14.5)') 'in sub determine_num_orbs_per_gridpoint: iproc, time for check_gridpoint', t_check_gridpoint

end subroutine determine_num_orbs_per_gridpoint_new








subroutine determine_communication_arrays(iproc, nproc, orbs, lzd, istartend_c, istartend_f, &
           index_in_global_c, index_in_global_f, &
           weightp_c, weightp_f,  nsendcounts_c, nsenddspls_c, nrecvcounts_c, nrecvdspls_c, &
           nsendcounts_f, nsenddspls_f, nrecvcounts_f, nrecvdspls_f)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(orbitals_data),intent(in) :: orbs
  type(local_zone_descriptors),intent(in) :: lzd
  integer,dimension(2,0:nproc-1),intent(in) :: istartend_c, istartend_f
  integer,dimension(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3),intent(in) :: index_in_global_c, index_in_global_f
  real(kind=8),intent(in) :: weightp_c, weightp_f
  integer,dimension(0:nproc-1),intent(out) :: nsendcounts_c, nsenddspls_c, nrecvcounts_c, nrecvdspls_c
  integer,dimension(0:nproc-1),intent(out) :: nsendcounts_f, nsenddspls_f, nrecvcounts_f, nrecvdspls_f
  
  ! Local variables
  integer :: iorb, iiorb, i1, i2, i3, ii, jproc, jproctarget, ierr, jj, ilr, j0, j1, i0, i, ind
  integer :: istat, ii1, ii2, ii3, iseg, istart, iend, iall
  integer,dimension(:),allocatable :: nsendcounts_tmp, nsenddspls_tmp, nrecvcounts_tmp, nrecvdspls_tmp
  character(len=*),parameter :: subname='determine_communication_arrays'


  !if(iproc==0) then
  !    do jproc=0,nproc-1
  !        write(*,'(a,i6,3i10)') 'iproc, istartend_c(:,jproc), lzd%glr%wfd%nvctr_c', iproc, istartend_c(:,jproc), lzd%glr%wfd%nvctr_c
  !        write(*,'(a,i6,3i10)') 'iproc, istartend_f(:,jproc), lzd%glr%wfd%nvctr_f', iproc, istartend_f(:,jproc), lzd%glr%wfd%nvctr_f
  !    end do
  !end if

  ! Determine values for mpi_alltoallv
  ! first nsendcounts
  nsendcounts_c=0
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
       ii2=i2+lzd%llr(ilr)%ns2
       ii3=i3+lzd%llr(ilr)%ns3
       do i=i0,i1
          ii1=i+lzd%llr(ilr)%ns1
          !call get_index_in_global(lzd%glr, ii1, ii2, ii3, 'c', ind)
          ind=index_in_global_c(ii1,ii2,ii3)
          jproctarget=-1
          do jproc=0,nproc-1
              if(ind>=istartend_c(1,jproc) .and. ind<=istartend_c(2,jproc)) then
                  jproctarget=jproc
                  exit
              end if
          end do
          !if(jproctarget==-1) write(*,*) 'ind, lzd%glr%wfd%nvctr_c',ind, lzd%glr%wfd%nvctr_c
          nsendcounts_c(jproctarget)=nsendcounts_c(jproctarget)+1
        end do
     end do
   end do

   !write(*,'(a,i3,3x,100i8)') 'iproc, istartend_f(2,:)', iproc, istartend_f(2,:)

  nsendcounts_f=0
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
   end do



  ! The first check is to make sure that there is no stop in case this process has no orbitals (in which case
  ! orbs%npsidim_orbs is 1 and not 0 as assumed by the check)
  if(orbs%npsidim_orbs>1 .and. sum(nsendcounts_c)+7*sum(nsendcounts_f)/=orbs%npsidim_orbs) then
      write(*,'(a,2i10)') 'sum(nsendcounts_c)+sum(nsendcounts_f)/=orbs%npsidim_orbs', &
                          sum(nsendcounts_c)+sum(nsendcounts_f), orbs%npsidim_orbs
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
           nrecvcounts_tmp, nrecvdspls_tmp, mpi_integer, mpi_comm_world, ierr)
      call mpi_alltoallv(nsendcounts_f, nsendcounts_tmp, nsenddspls_tmp, mpi_integer, nrecvcounts_f, &
           nrecvcounts_tmp, nrecvdspls_tmp, mpi_integer, mpi_comm_world, ierr)
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

  !write(*,*) 'sum(nrecvcounts_c), nint(weightp_c)', sum(nrecvcounts_c), nint(weightp_c)
  !write(*,*) 'sum(nrecvcounts_f), nint(weightp_f)', sum(nrecvcounts_f), nint(weightp_f)
  if(sum(nrecvcounts_c)/=nint(weightp_c)) stop 'sum(nrecvcounts_c)/=nint(nweightp_c)'
  if(sum(nrecvcounts_f)/=nint(weightp_f)) stop 'sum(nrecvcounts_f)/=nint(nweightp_f)'

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
  integer,dimension(:),allocatable :: nsend, indexsendorbital2, gridpoint_start_c, gridpoint_start_f, indexrecvorbital2
  real(kind=8),dimension(:,:,:),allocatable :: weight_c, weight_f
  integer,dimension(:),allocatable :: indexsendorbital_c, indexsendbuf_c, indexrecvbuf_c
  integer,dimension(:),allocatable :: indexsendorbital_f, indexsendbuf_f, indexrecvbuf_f
  character(len=*),parameter :: subname='get_switch_indices'
  !real(kind=8) :: t1, t2, t1tot, t2tot, t_reverse
  
  !t_reverse=0.d0
  !t1tot=mpi_wtime()
  
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

  allocate(nsend(0:nproc-1), stat=istat)
  call memocc(istat, nsend, 'nsend', subname)

  iitot=0
  nsend=0
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
              nsend(jproctarget)=nsend(jproctarget)+1
              ind=nsenddspls_c(jproctarget)+nsend(jproctarget)
              isendbuf_c(iitot)=ind
              indexsendbuf_c(ind)=indglob
              indexsendorbital_c(iitot)=iiorb
              !indexsendorbital(ind)=iiorb
          end do
      end do
  end do

  if(iitot/=ndimpsi_c) stop 'iitot/=ndimpsi_c'

  !check
  do jproc=0,nproc-1
      if(nsend(jproc)/=nsendcounts_c(jproc)) stop 'nsend(jproc)/=nsendcounts_c(jproc)'
  end do


  ! fine part
  iitot=0
  nsend=0
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
                  nsend(jproctarget)=nsend(jproctarget)+1
                  ind=nsenddspls_f(jproctarget)+nsend(jproctarget)
                  isendbuf_f(iitot)=ind
                  indexsendbuf_f(ind)=indglob
                  indexsendorbital_f(iitot)=iiorb
                  !indexsendorbital(ind)=iiorb
          end do
      end do
  end do

  if(iitot/=ndimpsi_f) stop 'iitot/=ndimpsi_f'

  !check
  do jproc=0,nproc-1
      !write(*,*) 'nsend(jproc), nsendcounts_f(jproc)', nsend(jproc), nsendcounts_f(jproc)
      if(nsend(jproc)/=nsendcounts_f(jproc)) stop 'nsend(jproc)/=nsendcounts_f(jproc)'
  end do






  allocate(indexsendorbital2(ndimpsi_c), stat=istat)
  call memocc(istat, indexsendorbital2, 'indexsendorbital2', subname)
  indexsendorbital2=indexsendorbital_c
  do i=1,ndimpsi_c
      ind=isendbuf_c(i)
      indexsendorbital_c(ind)=indexsendorbital2(i)
  end do
  ! Inverse of isendbuf
!t1=mpi_wtime()
  call get_reverse_indices(ndimpsi_c, isendbuf_c, irecvbuf_c)
  !do i=1,ndimpsi_c
  !    do j=1,ndimpsi_c
  !        if(isendbuf_c(j)==i) then
  !            irecvbuf_c(i)=j
  !        end if
  !    end do
  !end do
!t2=mpi_wtime()
!t_reverse=t_reverse+t2-t1
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
!t1=mpi_wtime()
  call get_reverse_indices(ndimpsi_f, isendbuf_f, irecvbuf_f)
  !do i=1,ndimpsi_f
  !    do j=1,ndimpsi_f
  !        if(isendbuf_f(j)==i) then
  !            irecvbuf_f(i)=j
  !        end if
  !    end do
  !end do
!t2=mpi_wtime()
!t_reverse=t_reverse+t2-t1
  iall=-product(shape(indexsendorbital2))*kind(indexsendorbital2)
  deallocate(indexsendorbital2, stat=istat)
  call memocc(istat, iall, 'indexsendorbital2', subname)




  if(nproc>1) then
      ! Communicate indexsendbuf
      call mpi_alltoallv(indexsendbuf_c, nsendcounts_c, nsenddspls_c, mpi_integer, indexrecvbuf_c, &
           nrecvcounts_c, nrecvdspls_c, mpi_integer, mpi_comm_world, ierr)
      ! Communicate indexsendorbitals
      call mpi_alltoallv(indexsendorbital_c, nsendcounts_c, nsenddspls_c, mpi_integer, indexrecvorbital_c, &
           nrecvcounts_c, nrecvdspls_c, mpi_integer, mpi_comm_world, ierr)

      ! Communicate indexsendbuf
      call mpi_alltoallv(indexsendbuf_f, nsendcounts_f, nsenddspls_f, mpi_integer, indexrecvbuf_f, &
           nrecvcounts_f, nrecvdspls_f, mpi_integer, mpi_comm_world, ierr)
      ! Communicate indexsendorbitals
      call mpi_alltoallv(indexsendorbital_f, nsendcounts_f, nsenddspls_f, mpi_integer, indexrecvorbital_f, &
           nrecvcounts_f, nrecvdspls_f, mpi_integer, mpi_comm_world, ierr)
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
      jj=ii-1
      i3=jj/((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1))
      jj=jj-i3*(lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)
      i2=jj/(lzd%glr%d%n1+1)
      i1=jj-i2*(lzd%glr%d%n1+1)
      !if(weight_c(i1,i2,i3)==0.d0) stop 'coarse: weight is zero!'
      !if(weight_c(i1,i2,i3)>0.d0) then
      !    if(gridpoint_start_c(ii)==0) then
      !        write(*,'(a,5i8)') 'DEBUG: iproc, jj, i1, i2, i3', iproc, jj, i1, i2, i3
      !        stop 'coarse: weight>0, but gridpoint_start(ii)==0'
      !    end if
      !end if

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
      jj=ii-1
      i3=jj/((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1))
      jj=jj-i3*(lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)
      i2=jj/(lzd%glr%d%n1+1)
      i1=jj-i2*(lzd%glr%d%n1+1)
      !if(weight_f(i1,i2,i3)==0.d0) stop 'fine: weight is zero!'
      !if(weight_f(i1,i2,i3)>0.d0) then
      !    if(gridpoint_start_f(ii)==0) then
      !        write(*,'(a,5i8)') 'DEBUG: iproc, jj, i1, i2, i3', iproc, jj, i1, i2, i3
      !        stop 'fine: weight>0, but gridpoint_start(ii)==0'
      !    end if
      !end if

      ind=gridpoint_start_f(ii)
      !if(ind==0) stop 'ind is zero!'
      iextract_f(i)=ind
      gridpoint_start_f(ii)=gridpoint_start_f(ii)+1  
  end do
  !if(sum(iextract_f)/=nint(weightp_f*(weightp_f+1.d0)*.5d0,kind=8)) then
  !  print*,sum(real(iextract_f,dp)),nint(weightp_f*(weightp_f+1.d0)*.5d0,kind=8)
  !  stop 'sum(iextract_f)/=nint(weightp_f*(weightp_f+1.d0)*.5d0)'
  !end if
  if(maxval(iextract_f)>sum(nrecvcounts_f)) stop 'maxval(iextract_f)>sum(nrecvcounts_f)'
  if(minval(iextract_f)<1) stop 'minval(iextract_f)<1'




  ! Get the array to transfrom back the data
!t1=mpi_wtime()
  call get_reverse_indices(sum(nrecvcounts_c), iextract_c, iexpand_c)
  !do i=1,sum(nrecvcounts_c)
  !    do j=1,sum(nrecvcounts_c)
  !        if(iextract_c(j)==i) then
  !            iexpand_c(i)=j
  !        end if
  !    end do
  !end do
!t2=mpi_wtime()
!t_reverse=t_reverse+t2-t1

!t1=mpi_wtime()
  call get_reverse_indices(sum(nrecvcounts_f), iextract_f, iexpand_f)
  !do i=1,sum(nrecvcounts_f)
  !    do j=1,sum(nrecvcounts_f)
  !        if(iextract_f(j)==i) then
  !            iexpand_f(i)=j
  !        end if
  !    end do
  !end do
!t2=mpi_wtime()
!t_reverse=t_reverse+t2-t1
  



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

  iall=-product(shape(nsend))*kind(nsend)
  deallocate(nsend, stat=istat)
  call memocc(istat, iall, 'nsend', subname)

!t2tot=mpi_wtime()
!write(*,'(a,i6,es15.5)') 'in sub get_switch_indices: iproc, time reverse',iproc, t_reverse
!write(*,'(a,i6,es15.5)') 'in sub get_switch_indices: iproc, total time',iproc, t2tot-t1tot

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
  gridpoint_start_c=0
  do i3=0,lzd%glr%d%n3
      do i2=0,lzd%glr%d%n2
          do i1=0,lzd%glr%d%n1
              i=i+1
              if(weight_c(i1,i2,i3)>0.d0) then
                  gridpoint_start_c(i)=ii
                  ii=ii+nint(weight_c(i1,i2,i3))
              end if
          end do
      end do
  end do

  ! CHECK
  i=0
  do i3=0,lzd%glr%d%n3
      do i2=0,lzd%glr%d%n2
          do i1=0,lzd%glr%d%n1
              i=i+1
              if(weight_c(i1,i2,i3)>0.d0) then
                  if(gridpoint_start_c(i)==0) stop 'FIRST CHECK: ERROR'
              end if
          end do
      end do
  end do



  ! fine part
  !weight_f=0.d0
  call to_zero((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(lzd%glr%d%n3+1), weight_f(0,0,0))
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
  gridpoint_start_f=0
  do i3=0,lzd%glr%d%n3
      do i2=0,lzd%glr%d%n2
          do i1=0,lzd%glr%d%n1
              i=i+1
              if(weight_f(i1,i2,i3)>0.d0) then
                  gridpoint_start_f(i)=ii
                  ii=ii+nint(weight_f(i1,i2,i3))
              end if
          end do
      end do
  end do

  ! CHECK
  i=0
  do i3=0,lzd%glr%d%n3
      do i2=0,lzd%glr%d%n2
          do i1=0,lzd%glr%d%n1
              i=i+1
              if(weight_f(i1,i2,i3)>0.d0) then
                  if(gridpoint_start_f(i)==0) stop 'FIRST CHECK: ERROR'
              end if
          end do
      end do
  end do


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
  
  !call mpi_comm_rank(mpi_comm_world, iproc, i)

  found=.false.
  !write(300+iproc,*) '---start---'
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
  !write(300+iproc,*) 'new iseg_start:',iseg_start
  !write(300+iproc,*) '--- end ---'


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





subroutine transpose_switch_psi(orbs, collcom, psi, psiwork_c, psiwork_f, lzd)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  type(orbitals_Data),intent(in) :: orbs
  type(collective_comms),intent(in) :: collcom
  real(kind=8),dimension(orbs%npsidim_orbs),intent(in) :: psi
  real(kind=8),dimension(collcom%ndimpsi_c),intent(out) :: psiwork_c
  real(kind=8),dimension(7*collcom%ndimpsi_f),intent(out) :: psiwork_f
  type(local_zone_descriptors),intent(in),optional :: lzd
  
  ! Local variables
  integer :: i_tot, i_c, i_f, iorb, iiorb, ilr, i, ind, istat, iall
  real(kind=8),dimension(:),allocatable :: psi_c, psi_f
  character(len=*),parameter :: subname='transpose_switch_psi'
  
  allocate(psi_c(collcom%ndimpsi_c), stat=istat)
  call memocc(istat, psi_c, 'psi_c', subname)
  allocate(psi_f(7*collcom%ndimpsi_f), stat=istat)
  call memocc(istat, psi_f, 'psi_f', subname)
  
  if(present(lzd)) then
      ! split up psi into coarse and fine part
      

  !$omp parallel default(shared) &
  !$omp private(i,i_c,i_tot,i_f,iiorb,iorb,ilr) 
  
      i_tot=0
      i_c=0
      i_f=0

      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          ilr=orbs%inwhichlocreg(iiorb)
	  !$omp do
          do i=1,lzd%llr(ilr)%wfd%nvctr_c
              psi_c(i_c+i)=psi(i_tot+i)
          end do
	  !$omp end do
	  
	  i_c = i_c + lzd%llr(ilr)%wfd%nvctr_c
	  i_tot = i_tot + lzd%llr(ilr)%wfd%nvctr_c
	  
	  !$omp do
          do i=1,7*lzd%llr(ilr)%wfd%nvctr_f
              psi_f(i_f + i)=psi(i_tot+i)
          end do
	  !$omp end do

	  i_f = i_f + 7*lzd%llr(ilr)%wfd%nvctr_f
          i_tot = i_tot + 7*lzd%llr(ilr)%wfd%nvctr_f

      end do

  
  !$omp end parallel

  else
      ! only coarse part is used...
      call dcopy(collcom%ndimpsi_c, psi, 1, psi_c, 1)
  end if
  
  ! coarse part

  !$omp parallel default(private) &
  !$omp shared(orbs, collcom, psi, psiwork_c, psiwork_f, lzd, psi_c,psi_f)
  !$omp do

  do i=1,collcom%ndimpsi_c
      ind=collcom%isendbuf_c(i)
      psiwork_c(ind)=psi_c(i)
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
  integer :: ierr, istat, iall, ist, ist_c, ist_f, jproc, iisend, iirecv
  real(kind=8),dimension(:),allocatable :: psiwork, psitwork
  integer,dimension(:),allocatable :: nsendcounts, nsenddspls, nrecvcounts, nrecvdspls
  character(len=*),parameter :: subname='transpose_communicate_psi'

  !call mpi_comm_size(mpi_comm_world, nproc, ierr)
  !call mpi_comm_rank(mpi_comm_world, iproc, ierr)

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

  ist=1
  ist_c=1
  ist_f=1
  iisend=0
  iirecv=0
  do jproc=0,nproc-1
      if(collcom%nsendcounts_c(jproc)>0) call dcopy(collcom%nsendcounts_c(jproc), psiwork_c(ist_c), 1, psiwork(ist), 1)
      ist_c=ist_c+collcom%nsendcounts_c(jproc)
      ist=ist+collcom%nsendcounts_c(jproc)
      if(collcom%nsendcounts_f(jproc)>0) call dcopy(7*collcom%nsendcounts_f(jproc), psiwork_f(ist_f), 1, psiwork(ist), 1)
      ist_f=ist_f+7*collcom%nsendcounts_f(jproc)
      ist=ist+7*collcom%nsendcounts_f(jproc)
      nsendcounts(jproc)=collcom%nsendcounts_c(jproc)+7*collcom%nsendcounts_f(jproc)
      nsenddspls(jproc)=iisend
      nrecvcounts(jproc)=collcom%nrecvcounts_c(jproc)+7*collcom%nrecvcounts_f(jproc)
      nrecvdspls(jproc)=iirecv
      iisend=iisend+nsendcounts(jproc)
      iirecv=iirecv+nrecvcounts(jproc)
  end do

  !write(*,'(a,i4,4x,100i8)') 'iproc, nsendcounts', iproc, nsendcounts
  !write(*,'(a,i4,4x,100i8)') 'iproc, nsenddspls', iproc, nsenddspls
  !write(*,'(a,i4,4x,100i8)') 'iproc, nrecvcounts', iproc, nrecvcounts
  !write(*,'(a,i4,4x,100i8)') 'iproc, nrecvdspls', iproc, nrecvdspls
  
  !! coarse part
  !call mpi_alltoallv(psiwork_c, collcom%nsendcounts_c, collcom%nsenddspls_c, mpi_double_precision, psitwork_c, &
  !     collcom%nrecvcounts_c, collcom%nrecvdspls_c, mpi_double_precision, mpi_comm_world, ierr)
  !
  !! fine part
  !call mpi_alltoallv(psiwork_f, 7*collcom%nsendcounts_f, 7*collcom%nsenddspls_f, mpi_double_precision, psitwork_f, &
  !     7*collcom%nrecvcounts_f, 7*collcom%nrecvdspls_f, mpi_double_precision, mpi_comm_world, ierr)
  call mpi_alltoallv(psiwork, nsendcounts, nsenddspls, mpi_double_precision, psitwork, &
       nrecvcounts, nrecvdspls, mpi_double_precision, mpi_comm_world, ierr)

  ist=1
  ist_c=1
  ist_f=1
  do jproc=0,nproc-1
      if(collcom%nrecvcounts_c(jproc)>0) call dcopy(collcom%nrecvcounts_c(jproc), psitwork(ist), 1, psitwork_c(ist_c), 1)
      ist_c=ist_c+collcom%nrecvcounts_c(jproc)
      ist=ist+collcom%nrecvcounts_c(jproc)
      if(collcom%nrecvcounts_f(jproc)>0) call dcopy(7*collcom%nrecvcounts_f(jproc), psitwork(ist), 1, psitwork_f(ist_f), 1)
      ist_f=ist_f+7*collcom%nrecvcounts_f(jproc)
      ist=ist+7*collcom%nrecvcounts_f(jproc)
  end do

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
  integer :: i, ind, sum_c,sum_f

  sum_c = sum(collcom%nrecvcounts_c)
  sum_f = sum(collcom%nrecvcounts_f)

  !$omp parallel private(i,ind) &
  !$omp shared(psit_c,psit_f, psitwork_c, psitwork_f,collcom,sum_c,sum_f)

  ! coarse part

  !$omp do

  do i=1, sum_c
      ind=collcom%iextract_c(i)
      psit_c(ind)=psitwork_c(i)
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
  integer :: i, ind, sum_c,sum_f

  sum_c = sum(collcom%nrecvcounts_c)
  sum_f = sum(collcom%nrecvcounts_f)

  !$omp parallel default(private) &
  !$omp shared(collcom, psit_c,psit_f, psitwork_c, psitwork_f,sum_c,sum_f)

  ! coarse part

  !$omp do

  do i=1,sum_c
      ind=collcom%iexpand_c(i)
      psitwork_c(ind)=psit_c(i)
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
  integer :: ierr, istat, iall, ist, ist_c, ist_f, jproc, iisend, iirecv
  real(kind=8),dimension(:),allocatable :: psiwork, psitwork
  integer,dimension(:),allocatable :: nsendcounts, nsenddspls, nrecvcounts, nrecvdspls
  character(len=*),parameter :: subname='transpose_communicate_psit'

  !call mpi_comm_size(mpi_comm_world, nproc, ierr)
  !call mpi_comm_rank(mpi_comm_world, iproc, ierr)

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

  ist=1
  ist_c=1
  ist_f=1
  iisend=0
  iirecv=0
  do jproc=0,nproc-1
      if(collcom%nrecvcounts_c(jproc)>0) call dcopy(collcom%nrecvcounts_c(jproc), psitwork_c(ist_c), 1, psitwork(ist), 1)
      ist_c=ist_c+collcom%nrecvcounts_c(jproc)
      ist=ist+collcom%nrecvcounts_c(jproc)
      if(collcom%nrecvcounts_f(jproc)>0) call dcopy(7*collcom%nrecvcounts_f(jproc), psitwork_f(ist_f), 1, psitwork(ist), 1)
      ist_f=ist_f+7*collcom%nrecvcounts_f(jproc)
      ist=ist+7*collcom%nrecvcounts_f(jproc)
      nsendcounts(jproc)=collcom%nsendcounts_c(jproc)+7*collcom%nsendcounts_f(jproc)
      nsenddspls(jproc)=iisend
      nrecvcounts(jproc)=collcom%nrecvcounts_c(jproc)+7*collcom%nrecvcounts_f(jproc)
      nrecvdspls(jproc)=iirecv
      iisend=iisend+nsendcounts(jproc)
      iirecv=iirecv+nrecvcounts(jproc)
  end do


 !! coarse part
 ! call mpi_alltoallv(psitwork_c, collcom%nrecvcounts_c, collcom%nrecvdspls_c, mpi_double_precision, psiwork_c, &
 !      collcom%nsendcounts_c, collcom%nsenddspls_c, mpi_double_precision, mpi_comm_world, ierr)

 !! fine part
 ! call mpi_alltoallv(psitwork_f, 7*collcom%nrecvcounts_f, 7*collcom%nrecvdspls_f, mpi_double_precision, psiwork_f, &
 !      7*collcom%nsendcounts_f, 7*collcom%nsenddspls_f, mpi_double_precision, mpi_comm_world, ierr)
  call mpi_alltoallv(psitwork, nrecvcounts, nrecvdspls, mpi_double_precision, psiwork, &
       nsendcounts, nsenddspls, mpi_double_precision, mpi_comm_world, ierr)

  ist=1
  ist_c=1
  ist_f=1
  do jproc=0,nproc-1
      if(collcom%nsendcounts_c(jproc)>0) call dcopy(collcom%nsendcounts_c(jproc), psiwork(ist), 1, psiwork_c(ist_c), 1)
      ist_c=ist_c+collcom%nsendcounts_c(jproc)
      ist=ist+collcom%nsendcounts_c(jproc)
      if(collcom%nsendcounts_f(jproc)>0) call dcopy(7*collcom%nsendcounts_f(jproc), psiwork(ist), 1, psiwork_f(ist_f), 1)
      ist_f=ist_f+7*collcom%nsendcounts_f(jproc)
      ist=ist+7*collcom%nsendcounts_f(jproc)
  end do

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

end subroutine transpose_communicate_psit



subroutine transpose_unswitch_psi(orbs, collcom, psiwork_c, psiwork_f, psi, lzd)
  use module_base
  use module_types
  implicit none
  
  ! Caling arguments
  type(orbitals_data),intent(in) :: orbs
  type(collective_comms),intent(in) :: collcom
  real(kind=8),dimension(collcom%ndimpsi_c),intent(in) :: psiwork_c
  real(kind=8),dimension(7*collcom%ndimpsi_f),intent(in) :: psiwork_f
  real(kind=8),dimension(orbs%npsidim_orbs),intent(out) :: psi
  type(local_zone_descriptors),intent(in),optional :: lzd
  
  ! Local variables
  integer :: i, ind, iorb, iiorb, ilr, i_tot, i_c, i_f, istat, iall
  real(kind=8),dimension(:),allocatable :: psi_c, psi_f
  character(len=*),parameter :: subname='transpose_unswitch_psi'
  
  
  allocate(psi_c(collcom%ndimpsi_c), stat=istat)
  call memocc(istat, psi_c, 'psi_c', subname)
  allocate(psi_f(7*collcom%ndimpsi_f), stat=istat)
  call memocc(istat, psi_f, 'psi_f', subname)
  
  !$omp parallel default(private) &
  !$omp shared(collcom, psiwork_c, psi_c,psi_f,psiwork_f)

  ! coarse part

  !$omp do
    do i=1,collcom%ndimpsi_c
        ind=collcom%irecvbuf_c(i)
        psi_c(ind)=psiwork_c(i)
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

  !$omp parallel default(shared) &
  !$omp private(i,i_c,i_tot,i_f,iiorb,iorb,ilr) 

        i_tot=0
        i_c=0
        i_f=0
        do iorb=1,orbs%norbp
            iiorb=orbs%isorb+iorb
            ilr=orbs%inwhichlocreg(iiorb)
            !!$omp do
            do i=1,lzd%llr(ilr)%wfd%nvctr_c
                psi(i_tot+i)=psi_c(i_c+i)
            end do
            !!$omp end do
	    i_c = i_c + lzd%llr(ilr)%wfd%nvctr_c
            i_tot = i_tot + lzd%llr(ilr)%wfd%nvctr_c
            !!$omp do
            do i=1,7*lzd%llr(ilr)%wfd%nvctr_f
                psi(i_tot+i)=psi_f(i_f+i)
            end do
            !!$omp end do
   	
	    i_f = i_f + 7*lzd%llr(ilr)%wfd%nvctr_f
            i_tot = i_tot + 7*lzd%llr(ilr)%wfd%nvctr_f


        end do
    !$omp end parallel 

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




subroutine transpose_localized(iproc, nproc, orbs, collcom, psi, psit_c, psit_f, lzd)
  use module_base
  use module_types
  use module_interfaces, except_this_one => transpose_localized
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(orbitals_data),intent(in) :: orbs
  type(collective_comms),intent(in) :: collcom
  real(kind=8),dimension(orbs%npsidim_orbs),intent(in) :: psi
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
      call transpose_switch_psi(orbs, collcom, psi, psiwork_c, psiwork_f, lzd)
  else
      call transpose_switch_psi(orbs, collcom, psi, psiwork_c, psiwork_f)
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



subroutine untranspose_localized(iproc, nproc, orbs, collcom, psit_c, psit_f, psi, lzd)
  use module_base
  use module_types
  use module_interfaces, except_this_one => untranspose_localized
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(orbitals_data),intent(in) :: orbs
  type(collective_comms),intent(in) :: collcom
  real(kind=8),dimension(collcom%ndimind_c),intent(in) :: psit_c
  real(kind=8),dimension(7*collcom%ndimind_f),intent(in) :: psit_f
  real(kind=8),dimension(orbs%npsidim_orbs),intent(out) :: psi
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
      call transpose_unswitch_psi(orbs, collcom, psiwork_c, psiwork_f, psi, lzd)
  else
      call transpose_unswitch_psi(orbs, collcom, psiwork_c, psiwork_f, psi)
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



subroutine calculate_overlap_transposed(iproc, nproc, orbs, mad, collcom, psit_c1, psit_c2, psit_f1, psit_f2, ovrlp)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(orbitals_data),intent(in) :: orbs
  type(matrixDescriptors),intent(in) :: mad
  type(collective_comms),intent(in) :: collcom
  real(kind=8),dimension(collcom%ndimind_c),intent(in) :: psit_c1, psit_c2
  real(kind=8),dimension(7*collcom%ndimind_f),intent(in) :: psit_f1, psit_f2
  real(kind=8),dimension(orbs%norb,orbs%norb),intent(out) :: ovrlp
  
  real(8),dimension(collcom%nptsp_c) :: dummy_c,dummy_f  

  ! Local variables
  integer :: i0, ipt, ii, iiorb, j, jjorb, i, ierr, istat, iall, m
  real(kind=8),dimension(:),allocatable :: ovrlp_compr
  character(len=*),parameter :: subname='calculate_overlap_transposed'
  
  call timing(iproc,'ovrlptransComp','ON') !lr408t
  call to_zero(orbs%norb**2, ovrlp(1,1))

  !$omp parallel default(private) &
  !$omp shared(collcom, ovrlp, psit_c1, psit_c2, psit_f1, psit_f2)


  !$omp do reduction (+:ovrlp) 

  do ipt=1,collcom%nptsp_c 
      ii=collcom%norb_per_gridpoint_c(ipt) 
      i0 = collcom%isptsp_c(ipt)

      do i=1,ii
          iiorb=collcom%indexrecvorbital_c(i0+i)
          m=mod(ii,4)
          if(m/=0) then

              do j=1,m
                  jjorb=collcom%indexrecvorbital_c(i0+j)
                  ovrlp(jjorb,iiorb)=ovrlp(jjorb,iiorb)+psit_c1(i0+i)*psit_c2(i0+j)
              end do
    
          end if
       
          do j=m+1,ii,4
              jjorb=collcom%indexrecvorbital_c(i0+j+0)
              ovrlp(jjorb,iiorb)=ovrlp(jjorb,iiorb)+psit_c1(i0+i)*psit_c2(i0+j+0)
              jjorb=collcom%indexrecvorbital_c(i0+j+1)
              ovrlp(jjorb,iiorb)=ovrlp(jjorb,iiorb)+psit_c1(i0+i)*psit_c2(i0+j+1)
              jjorb=collcom%indexrecvorbital_c(i0+j+2)
              ovrlp(jjorb,iiorb)=ovrlp(jjorb,iiorb)+psit_c1(i0+i)*psit_c2(i0+j+2)
              jjorb=collcom%indexrecvorbital_c(i0+j+3)
              ovrlp(jjorb,iiorb)=ovrlp(jjorb,iiorb)+psit_c1(i0+i)*psit_c2(i0+j+3)
          end do
         
          !do j=1,ii
          !    jjorb=collcom%indexrecvorbital_c(i0+j)
          !    ovrlp(jjorb,iiorb)=ovrlp(jjorb,iiorb)+psit_c1(i0+i)*psit_c2(i0+j)
          !end do
      end do

      !!$omp end do

  end do
  !$omp end do
  
  !$omp do reduction(+:ovrlp) 

  do ipt=1,collcom%nptsp_f 
      ii=collcom%norb_per_gridpoint_f(ipt) 

      i0 = collcom%isptsp_f(ipt)

      do i=1,ii
          iiorb=collcom%indexrecvorbital_f(i0+i)
          do j=1,ii
              jjorb=collcom%indexrecvorbital_f(i0+j)
              ovrlp(jjorb,iiorb)=ovrlp(jjorb,iiorb)+psit_f1(7*(i0+i)-6)*psit_f2(7*(i0+j)-6)
              ovrlp(jjorb,iiorb)=ovrlp(jjorb,iiorb)+psit_f1(7*(i0+i)-5)*psit_f2(7*(i0+j)-5)
              ovrlp(jjorb,iiorb)=ovrlp(jjorb,iiorb)+psit_f1(7*(i0+i)-4)*psit_f2(7*(i0+j)-4)
              ovrlp(jjorb,iiorb)=ovrlp(jjorb,iiorb)+psit_f1(7*(i0+i)-3)*psit_f2(7*(i0+j)-3)
              ovrlp(jjorb,iiorb)=ovrlp(jjorb,iiorb)+psit_f1(7*(i0+i)-2)*psit_f2(7*(i0+j)-2)
              ovrlp(jjorb,iiorb)=ovrlp(jjorb,iiorb)+psit_f1(7*(i0+i)-1)*psit_f2(7*(i0+j)-1)
              ovrlp(jjorb,iiorb)=ovrlp(jjorb,iiorb)+psit_f1(7*(i0+i)-0)*psit_f2(7*(i0+j)-0)
          end do
      end do

  end do

  !$omp end do
  !$omp end parallel

  call timing(iproc,'ovrlptransComp','OF') !lr408t

  call timing(iproc,'ovrlptransComm','ON') !lr408t

  if(nproc>1) then
      allocate(ovrlp_compr(mad%nvctr), stat=istat)
      call memocc(istat, ovrlp_compr, 'ovrlp_compr', subname)
      call compress_matrix_for_allreduce(orbs%norb, mad, ovrlp, ovrlp_compr)
      call mpiallred(ovrlp_compr(1), mad%nvctr, mpi_sum, mpi_comm_world, ierr)
      call uncompressMatrix(orbs%norb, mad, ovrlp_compr,ovrlp)
      iall=-product(shape(ovrlp_compr))*kind(ovrlp_compr)
      deallocate(ovrlp_compr, stat=istat)
      call memocc(istat, iall, 'ovrlp_compr', subname)
  end if
  call timing(iproc,'ovrlptransComm','OF') !lr408t
end subroutine calculate_overlap_transposed


! This will work because the difference between collcom1 and collcom2 is only a factor 3 between the orbital numbers.
! Hence, nptsp_c and nptsp_f should be the same, only the norb_per_gridpoint will change.
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
      call mpiallred(ovrlp(1,1), orbs1%norb*orbs2%norb, mpi_sum, mpi_comm_world, ierr)
  end if
  call timing(iproc,'ovrlptransComm','OF') !lr408t
end subroutine calculate_pulay_overlap

subroutine build_linear_combination_transposed(norb, matrix, collcom, psitwork_c, psitwork_f, reset, psit_c, psit_f, &
     iproc)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: norb
  real(kind=8),dimension(norb,norb),intent(in) :: matrix
  type(collective_comms),intent(in) :: collcom
  real(kind=8),dimension(collcom%ndimind_c),intent(in) :: psitwork_c
  real(kind=8),dimension(7*collcom%ndimind_f),intent(in) :: psitwork_f
  logical,intent(in) :: reset
  real(kind=8),dimension(collcom%ndimind_c),intent(out) :: psit_c
  real(kind=8),dimension(7*collcom%ndimind_f),intent(out) :: psit_f
  integer, intent(in) :: iproc
  ! Local variables
  integer :: i0, ipt, ii, j, iiorb, jjorb, i, m

  call timing(iproc,'lincombtrans  ','ON') !lr408t
  if(reset) then
      if(collcom%ndimind_c>0) call to_zero(collcom%ndimind_c, psit_c(1))
      if(collcom%ndimind_f>0) call to_zero(7*collcom%ndimind_f, psit_f(1))
  end if

  i0=0
  do ipt=1,collcom%nptsp_c 
      ii=collcom%norb_per_gridpoint_c(ipt) 
      do i=1,ii
          iiorb=collcom%indexrecvorbital_c(i0+i)
          m=mod(ii,4)
          if(m/=0) then
              do j=1,m
                  jjorb=collcom%indexrecvorbital_c(i0+j)
                  psit_c(i0+i)=psit_c(i0+i)+matrix(jjorb,iiorb)*psitwork_c(i0+j)
              end do
          end if
          do j=m+1,ii,4
              jjorb=collcom%indexrecvorbital_c(i0+j+0)
              psit_c(i0+i)=psit_c(i0+i)+matrix(jjorb,iiorb)*psitwork_c(i0+j+0)
              jjorb=collcom%indexrecvorbital_c(i0+j+1)
              psit_c(i0+i)=psit_c(i0+i)+matrix(jjorb,iiorb)*psitwork_c(i0+j+1)
              jjorb=collcom%indexrecvorbital_c(i0+j+2)
              psit_c(i0+i)=psit_c(i0+i)+matrix(jjorb,iiorb)*psitwork_c(i0+j+2)
              jjorb=collcom%indexrecvorbital_c(i0+j+3)
              psit_c(i0+i)=psit_c(i0+i)+matrix(jjorb,iiorb)*psitwork_c(i0+j+3)
          end do
          !do j=1,ii
          !    jjorb=collcom%indexrecvorbital_c(i0+j)
          !    psit_c(i0+i)=psit_c(i0+i)+matrix(jjorb,iiorb)*psitwork_c(i0+j)
          !end do
      end do
      i0=i0+ii
  end do

  i0=0
  do ipt=1,collcom%nptsp_f 
      ii=collcom%norb_per_gridpoint_f(ipt) 
      do i=1,ii
          iiorb=collcom%indexrecvorbital_f(i0+i)
          do j=1,ii
              jjorb=collcom%indexrecvorbital_f(i0+j)
              psit_f(7*(i0+i)-6) = psit_f(7*(i0+i)-6) + matrix(jjorb,iiorb)*psitwork_f(7*(i0+j)-6)
              psit_f(7*(i0+i)-5) = psit_f(7*(i0+i)-5) + matrix(jjorb,iiorb)*psitwork_f(7*(i0+j)-5)
              psit_f(7*(i0+i)-4) = psit_f(7*(i0+i)-4) + matrix(jjorb,iiorb)*psitwork_f(7*(i0+j)-4)
              psit_f(7*(i0+i)-3) = psit_f(7*(i0+i)-3) + matrix(jjorb,iiorb)*psitwork_f(7*(i0+j)-3)
              psit_f(7*(i0+i)-2) = psit_f(7*(i0+i)-2) + matrix(jjorb,iiorb)*psitwork_f(7*(i0+j)-2)
              psit_f(7*(i0+i)-1) = psit_f(7*(i0+i)-1) + matrix(jjorb,iiorb)*psitwork_f(7*(i0+j)-1)
              psit_f(7*(i0+i)-0) = psit_f(7*(i0+i)-0) + matrix(jjorb,iiorb)*psitwork_f(7*(i0+j)-0)
          end do
      end do
      i0=i0+ii
  end do
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
  integer :: i, j

  do i=1,n
      j=indices(i)
      reverse_indices(j)=i
  end do

end subroutine get_reverse_indices


subroutine compress_matrix_for_allreduce(n, mad, mat, mat_compr)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: n
  type(matrixDescriptors),intent(in) :: mad
  real(kind=8),dimension(n**2),intent(in) :: mat
  real(kind=8),dimension(mad%nvctr),intent(out) :: mat_compr

  ! Local variables
  integer :: jj, iseg, jorb

  jj=0
  do iseg=1,mad%nseg
      do jorb=mad%keyg(1,iseg),mad%keyg(2,iseg)
          jj=jj+1
          mat_compr(jj)=mat(jorb)
      end do
  end do

end subroutine compress_matrix_for_allreduce



subroutine normalize_transposed(iproc, nproc, orbs, collcom, psit_c, psit_f)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(orbitals_data),intent(in):: orbs
  type(collective_comms),intent(in):: collcom
  real(8),dimension(collcom%ndimind_c),intent(inout):: psit_c
  real(8),dimension(7*collcom%ndimind_f),intent(inout):: psit_f
  
  ! Local variables
  integer:: i0, ipt, ii, iiorb, i, ierr, istat, iall, iorb
  real(8),dimension(:),allocatable:: norm
  character(len=*),parameter:: subname='normslize_transposed'

  allocate(norm(orbs%norb), stat=istat)
  call memocc(istat, norm, 'norm', subname)
  call to_zero(orbs%norb, norm(1))

  i0=0
  do ipt=1,collcom%nptsp_c 
      ii=collcom%norb_per_gridpoint_c(ipt) 
      do i=1,ii
          iiorb=collcom%indexrecvorbital_c(i0+i)
          norm(iiorb)=norm(iiorb)+psit_c(i0+i)**2
      end do
      i0=i0+ii
  end do

  i0=0
  do ipt=1,collcom%nptsp_f 
      ii=collcom%norb_per_gridpoint_f(ipt) 
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
      i0=i0+ii
  end do

  if(nproc>1) then
      call mpiallred(norm(1), orbs%norb, mpi_sum, mpi_comm_world, ierr)
  end if
  

  do iorb=1,orbs%norb
      norm(iorb)=1.d0/sqrt(norm(iorb))
  end do


  i0=0
  do ipt=1,collcom%nptsp_c 
      ii=collcom%norb_per_gridpoint_c(ipt) 
      do i=1,ii
          iiorb=collcom%indexrecvorbital_c(i0+i)
          psit_c(i0+i)=psit_c(i0+i)*norm(iiorb)
      end do
      i0=i0+ii
  end do

  i0=0
  do ipt=1,collcom%nptsp_f 
      ii=collcom%norb_per_gridpoint_f(ipt) 
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
      i0=i0+ii
  end do


  iall=-product(shape(norm))*kind(norm)
  deallocate(norm, stat=istat)
  call memocc(istat, iall, 'norm', subname)

end subroutine normalize_transposed
