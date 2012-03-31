subroutine init_collective_comms(iproc, nproc, orbs, lzd, collcom)
use module_base
use module_types
use module_interfaces, except_this_one => init_collective_comms
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data),intent(in):: orbs
type(local_zone_descriptors),intent(in):: lzd
type(collective_comms),intent(out):: collcom

! Local variables
integer:: ii, istat
real(8),dimension(:,:,:),allocatable:: weight_c, weight_c_temp, weight_f, weight_f_temp
real(8):: weight_c_tot, weight_f_tot, weightp_c, weightp_f, tt, ierr
integer,dimension(:,:),allocatable:: istartend_c, istartend_f

allocate(weight_c(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3), stat=istat)
allocate(weight_f(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3), stat=istat)

call get_weights(orbs, lzd, weight_c, weight_f, weight_c_tot, weight_f_tot)


!!
!!  tt=weight_tot
!!  call mpi_alzd%llreduce(tt, weight_tot, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
!!  weight_temp=weight
!!  ii=(lzd%glr%ie1-lzd%glr%is1+1)*(lzd%glr%ie2-lzd%glr%is2+1)*(lzd%glr%ie3-lzd%glr%is3+1)
!!  call mpi_alzd%llreduce(weight_temp, weight, ii, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
!!
  if(iproc==0) write(*,'(a,2es14.5)') 'total weights (coarse / fine):',weight_c_tot, weight_f_tot

!!  ! Assign the grid points to the processes such that the work is equally dsitributed
  allocate(istartend_c(2,0:nproc-1))
  allocate(istartend_f(2,0:nproc-1))
call assign_weight_to_process(iproc, nproc, lzd, weight_c, weight_f, weight_c_tot, weight_f_tot, &
           istartend_c, istartend_f, weightp_c, weightp_f, collcom%nptsp_c, collcom%nptsp_f)
!!  call assign_weight_to_process(iproc, nproc, lzd%glr, weight, weight_tot, istartend, weightp, collcom%nptsp)
  write(*,'(a,i0,a,i0,a,es14.5)') 'process ',iproc,' has ',collcom%nptsp_c,' coarse points and weight ',weightp_c
  write(*,'(a,i0,a,i0,a,es14.5)') 'process ',iproc,' has ',collcom%nptsp_f,' fine points and weight ',weightp_f
!!!!  if(istartend(2,nproc-1)/=(lzd%glr%ie1-lzd%glr%is1+1)*(lzd%glr%ie2-lzd%glr%is2+1)*(lzd%glr%ie3-lzd%glr%is3+1)) stop 'istartend(2,nproc-1)/=(lzd%glr%ie1-lzd%glr%is1+1)*(lzd%glr%ie2-lzd%glr%is2+1)*(lzd%glr%ie3-lzd%glr%is3+1)'
!!

  ! some checks
  call mpi_allreduce(weightp_c, tt, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
  if(tt/=weight_c_tot) stop 'wrong partition of coarse weights'
  call mpi_allreduce(weightp_f, tt, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
  if(tt/=weight_f_tot) stop 'wrong partition of fine weights'
  call mpi_allreduce(collcom%nptsp_c, ii, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr)
  if(ii/=lzd%glr%wfd%nvctr_c) stop 'wrong partition of coarse grid points'
  call mpi_allreduce(collcom%nptsp_f, ii, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr)
  if(ii/=lzd%glr%wfd%nvctr_f) stop 'wrong partition of fine grid points'

!!  ! Allocate the keys
  allocate(collcom%norb_per_gridpoint_c(collcom%nptsp_c), stat=istat)
  allocate(collcom%norb_per_gridpoint_f(collcom%nptsp_f), stat=istat)
!!  call determine_num_orbs_per_gridpoint(iproc, nproc, orbs%norb, collcom%nptsp, lzd%glr, lzd%llr, istartend, collcom%orbs%norb_per_gridpoint)
call determine_num_orbs_per_gridpoint(iproc, nproc, orbs, lzd, istartend_c, istartend_f, weightp_c, weightp_f, collcom%nptsp_c, collcom%nptsp_f, &
           collcom%norb_per_gridpoint_c, collcom%norb_per_gridpoint_f)
!!
!!
!!  ! Determine values for mpi_alltoallv
!!  allocate(collcom%nsendcounts(0:nproc-1))
!!  allocate(collcom%nsenddspls(0:nproc-1))
!!  allocate(collcom%nrecvcounts(0:nproc-1))
!!  allocate(collcom%nrecvdspls(0:nproc-1))
!!  call determine_communication_arrays(iproc, nproc, orbs%norb, orbs%orbs%norbp, orbs%isorb, orbs%npsidim, lzd%glr, lzd%llr, istartend, weightp, collcom%nsendcounts, collcom%nsenddspls, collcom%nrecvcounts, collcom%nrecvdspls)
!!
!!
!!  ! Now rearrange the data on the process to communicate them
!!  allocate(collcom%irecvbuf(orbs%npsidim))
!!  allocate(collcom%indexrecvorbital(sum(collcom%nrecvcounts)))
!!  allocate(collcom%iextract(sum(collcom%nrecvcounts)))
!!  allocate(collcom%iexpand(sum(collcom%nrecvcounts)))
!!allocate(collcom%isendbuf(orbs%npsidim_orbs))
!!  call get_switch_indices(iproc, nproc, orbs%orbs%norbp, orbs%norb, orbs%isorb, orbs%npsidim, lzd%glr, lzd%llr, istartend, collcom%nsendcounts, collcom%nsenddspls, collcom%nrecvcounts, &
!!           collcom%nrecvdspls, weightp,  collcom%isendbuf, collcom%irecvbuf, collcom%indexrecvorbital, collcom%iextract, collcom%iexpand)
!!
!!
end subroutine init_collective_comms

subroutine get_weights(orbs, lzd, weight_c, weight_f, weight_c_tot, weight_f_tot)
use module_base
use module_types
implicit none

! Calling arguments
type(orbitals_data),intent(in):: orbs
type(local_zone_descriptors),intent(in):: lzd
real(8),dimension(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3),intent(out):: weight_c, weight_f
real(8),intent(out):: weight_c_tot, weight_f_tot

! Local variables
integer:: iorb, iiorb, i0, i1, i2, i3, ii, jj, iseg, ierr, ilr, istart, iend, i, j0, j1, ii1, ii2, ii3


  weight_c=0.d0
  weight_f=0.d0
  weight_c_tot=0.d0
  weight_f_tot=0.d0

  ! coarse part
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
       enddo
    enddo
  
    ! fine part
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
          weight_f(ii1,ii2,ii3)=weight_f(ii1,ii2,ii3)+7.d0
          weight_f_tot=weight_f_tot+7.d0
       enddo
    enddo
  end do


  call mpiallred(weight_c_tot, 1, mpi_sum, mpi_comm_world, ierr)
  call mpiallred(weight_f_tot, 1, mpi_sum, mpi_comm_world, ierr)
  ii=(lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)*(lzd%glr%d%n3+1)
  call mpiallred(weight_c(0,0,0), ii,  mpi_sum, mpi_comm_world, ierr)
  call mpiallred(weight_f(0,0,0), ii,  mpi_sum, mpi_comm_world, ierr)


end subroutine get_weights



subroutine assign_weight_to_process(iproc, nproc, lzd, weight_c, weight_f, weight_tot_c, weight_tot_f, &
           istartend_c, istartend_f, weightp_c, weightp_f, nptsp_c, nptsp_f)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(local_zone_descriptors),intent(in):: lzd
real(8),dimension(0:lzd%glr%d%n1,0:lzd%glr%d%n2,0:lzd%glr%d%n3),intent(in):: weight_c, weight_f
real(8),intent(in):: weight_tot_c, weight_tot_f
integer,dimension(2,0:nproc-1),intent(out):: istartend_c, istartend_f
real(8),intent(out):: weightp_c, weightp_f
integer,intent(out):: nptsp_c, nptsp_f

! Local variables
integer:: jproc, i1, i2, i3, ii, istartp_c, iendp_c, ii2, istartp_f, iendp_f, istart, iend, jj, j0, j1, i, iseg, i0, iitot, ierr
real(8):: tt, tt2, weight_c_ideal, weight_f_ideal

  weight_c_ideal=weight_tot_c/dble(nproc)
  weight_f_ideal=weight_tot_f/dble(nproc)
  if(iproc==0) write(*,'(a,2es14.5)') 'ideal weight per process (coarse / fine):',weight_c_ideal,weight_f_ideal

  jproc=0
  tt=0.d0
  tt2=0.d0
  iitot=0
  ii2=0
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
           if(tt>weight_c_ideal) then
               if(iproc==jproc) then
                   weightp_c=tt
                   nptsp_c=iitot
                   istartp_c=ii2+1
                   iendp_c=istartp_c+iitot-1
               end if
               istartend_c(1,jproc)=ii2+1
               istartend_c(2,jproc)=istartend_c(1,jproc)+iitot-1
               tt2=tt2+tt
               tt=0.d0
               ii2=ii2+iitot
               iitot=0
               jproc=jproc+1
           end if
       end do
   end do
  if(iproc==nproc-1) then
      ! Take the rest
      istartp_c=ii2+1
      iendp_c=istartp_c+iitot-1
      weightp_c=weight_tot_c-tt2
      nptsp_c=lzd%glr%wfd%nvctr_c-ii2
  end if
  istartend_c(1,nproc-1)=ii2+1
  istartend_c(2,nproc-1)=istartend_c(1,nproc-1)+iitot-1

  ! some check
  ii=istartend_c(2,iproc)-istartend_c(1,iproc)+1
  call mpiallred(ii, 1, mpi_sum, mpi_comm_world, ierr)
  if(ii/=lzd%glr%wfd%nvctr_c) stop 'ii/=lzd%glr%wfd%nvctr_c'


  jproc=0
  tt=0.d0
  tt2=0.d0
  iitot=0
  ii2=0
  weightp_f=0.d0
    istart=lzd%glr%wfd%nseg_c+min(1,lzd%glr%wfd%nseg_f)
    iend=istart+lzd%glr%wfd%nseg_f-1
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
           if(tt>weight_f_ideal) then
               if(iproc==jproc) then
                   weightp_f=tt
                   nptsp_f=iitot
                   istartp_f=ii2+1
                   iendp_f=istartp_f+iitot-1
               end if
               istartend_f(1,jproc)=ii2+1
               istartend_f(2,jproc)=istartend_f(1,jproc)+iitot-1
               tt2=tt2+tt
               tt=0.d0
               ii2=ii2+iitot
               iitot=0
               jproc=jproc+1
           end if
       end do
   end do
  if(iproc==nproc-1) then
      ! Take the rest
      istartp_f=ii2+1
      iendp_f=istartp_f+iitot-1
      weightp_f=weight_tot_f-tt2
      nptsp_f=lzd%glr%wfd%nvctr_f-ii2
  end if
  istartend_f(1,nproc-1)=ii2+1
  istartend_f(2,nproc-1)=istartend_f(1,nproc-1)+iitot-1

  ! some check
  ii=istartend_f(2,iproc)-istartend_f(1,iproc)+1
  call mpiallred(ii, 1, mpi_sum, mpi_comm_world, ierr)
  if(ii/=lzd%glr%wfd%nvctr_f) stop 'ii/=lzd%glr%wfd%nvctr_f'



end subroutine assign_weight_to_process



subroutine determine_num_orbs_per_gridpoint(iproc, nproc, orbs, lzd, istartend_c, istartend_f, weightp_c, weightp_f, nptsp_c, nptsp_f, &
           norb_per_gridpoint_c, norb_per_gridpoint_f)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, nptsp_c, nptsp_f
type(orbitals_data),intent(in):: orbs
type(local_zone_descriptors),intent(in):: lzd
integer,dimension(2,0:nproc-1),intent(in):: istartend_c, istartend_f
real(8),intent(in):: weightp_c, weightp_f
integer,dimension(nptsp_c),intent(out):: norb_per_gridpoint_c
integer,dimension(nptsp_f),intent(out):: norb_per_gridpoint_f

! Local variables
integer:: ii, iiorb, i1, i2, i3, iipt, iorb, iii, npgp, iseg, jj, j0, j1, iitot, ilr, i, istart, iend, i0
logical:: found

  iitot=0
  iiorb=0
  iipt=0
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
              iitot=iitot+1
              if(iitot>=istartend_c(1,iproc) .and. iitot<=istartend_c(2,iproc)) then
                  iipt=iipt+1
                  npgp=0
                  do iorb=1,orbs%norb
                      ilr=orbs%inwhichlocreg(iorb)
                      ! Check whether this orbitals extends here
                      call check_gridpoint(lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, &
                                      lzd%llr(ilr)%ns1, lzd%llr(ilr)%ns2, lzd%llr(ilr)%ns3, lzd%llr(ilr)%wfd%keygloc, &
                                      i, i2, i3, found)
                      if(found) then
                          npgp=npgp+1
                          iiorb=iiorb+1
                      end if
                  end do
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
              iitot=iitot+1
              if(iitot>=istartend_f(1,iproc) .and. iitot<=istartend_f(2,iproc)) then
                  iipt=iipt+1
                  npgp=0
                  do iorb=1,orbs%norb
                      ilr=orbs%inwhichlocreg(iorb)
                      ! Check whether this orbitals extends here
                      iii=lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)
                      call check_gridpoint(lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, &
                                      lzd%llr(ilr)%ns1, lzd%llr(ilr)%ns2, lzd%llr(ilr)%ns3, lzd%llr(ilr)%wfd%keygloc(1,iii), &
                                      i, i2, i3, found)
                      if(found) then
                          npgp=npgp+1
                          iiorb=iiorb+7
                      end if
                  end do
                  norb_per_gridpoint_f(iipt)=npgp
              end if
      end do
  end do

  if(iipt/=nptsp_f) stop 'iipt/=nptsp_f'
  if(iiorb/=nint(weightp_f)) stop 'iiorb/=weightp_f'




end subroutine determine_num_orbs_per_gridpoint



!!subroutine determine_communication_arrays(iproc, nproc, norb, norbp, isorb, ndimpsi, glr, llr, istartend, weightp, nsendcounts, nsenddspls, nrecvcounts, nrecvdspls)
!!use types
!!implicit none
!!include 'mpif.h'
!!
!!! Calling arguments
!!integer,intent(in):: iproc, nproc, norb, norbp, isorb, ndimpsi
!!type(locreg_descriptors),intent(in):: glr
!!type(locreg_descriptors),dimension(norb),intent(in):: llr
!!integer,dimension(2,0:nproc-1),intent(in):: istartend
!!real(8),intent(in):: weightp
!!integer,dimension(0:nproc-1),intent(out):: nsendcounts, nsenddspls, nrecvcounts, nrecvdspls
!!
!!! Local variables
!!integer:: iorb, iiorb, i1, i2, i3, ii, jproc, jproctarget, ierr
!!integer,dimension(:),allocatable:: nsendcounts_tmp, nsenddspls_tmp, nrecvcounts_tmp, nrecvdspls_tmp
!!
!!  ! Determine values for mpi_alltoallv
!!  ! first nsendcounts
!!  !allocate(nsendcounts(0:nproc-1))
!!  nsendcounts=0
!!  do iorb=1,norbp
!!      iiorb=isorb+iorb
!!      do i3=llr(iiorb)%is3,llr(iiorb)%ie3
!!          do i2=llr(iiorb)%is2,llr(iiorb)%ie2
!!              do i1=llr(iiorb)%is1,llr(iiorb)%ie1
!!                  ii = (i3-1)*(glr%ie1-glr%is1+1)*(glr%ie2-glr%is2+1) + (i2-1)*(glr%ie1-glr%is1+1) + i1
!!                  do jproc=0,nproc-1
!!                      if(ii>=istartend(1,jproc) .and. ii<=istartend(2,jproc)) then
!!                          jproctarget=jproc
!!                          exit
!!                      end if
!!                  end do
!!                  nsendcounts(jproctarget)=nsendcounts(jproctarget)+1
!!              end do
!!          end do
!!      end do
!!  end do
!!
!!  if(sum(nsendcounts)/=ndimpsi) stop 'sum(nsendcounts)/=ndimpsi'
!!  
!!  ! now nsenddspls
!!  !allocate(nsenddspls(0:nproc-1))
!!  nsenddspls(0)=0
!!  do jproc=1,nproc-1
!!      nsenddspls(jproc)=nsenddspls(jproc-1)+nsendcounts(jproc-1)
!!  end do
!!
!!  ! now nrecvcounts
!!  ! use an mpi_alltoallv to gather the date
!!  !allocate(nrecvcounts(0:nproc-1))
!!  allocate(nsendcounts_tmp(0:nproc-1))
!!  allocate(nsenddspls_tmp(0:nproc-1))
!!  allocate(nrecvcounts_tmp(0:nproc-1))
!!  allocate(nrecvdspls_tmp(0:nproc-1))
!!  nsendcounts_tmp=1
!!  nrecvcounts_tmp=1
!!  do jproc=0,nproc-1
!!      nsenddspls_tmp(jproc)=jproc
!!      nrecvdspls_tmp(jproc)=jproc
!!  end do
!!  call mpi_alltoallv(nsendcounts, nsendcounts_tmp, nsenddspls_tmp, mpi_integer, nrecvcounts, &
!!       nrecvcounts_tmp, nrecvdspls_tmp, mpi_integer, mpi_comm_world, ierr)
!!  deallocate(nsendcounts_tmp)
!!  deallocate(nsenddspls_tmp)
!!  deallocate(nrecvcounts_tmp)
!!  deallocate(nrecvdspls_tmp)
!!
!!  ! now recvdspls
!!  !allocate(nrecvdspls(0:nproc-1))
!!  nrecvdspls(0)=0
!!  do jproc=1,nproc-1
!!      nrecvdspls(jproc)=nrecvdspls(jproc-1)+nrecvcounts(jproc-1)
!!  end do
!!
!!  if(sum(nrecvcounts)/=nint(weightp)) stop 'sum(nrecvcounts)/=nint(nweightp)'
!!
!!end subroutine determine_communication_arrays
!!
!!
!!subroutine get_switch_indices(iproc, nproc, norbp, norb, isorb, ndimpsi, glr, llr, istartend, nsendcounts, nsenddspls, nrecvcounts, &
!!           nrecvdspls, weightp, isendbuf, irecvbuf, indexrecvorbital, iextract, iexpand)
!!use types
!!implicit none
!!include 'mpif.h'
!!
!!! Calling arguments
!!integer,intent(in):: iproc, nproc, norbp, norb, isorb, ndimpsi
!!type(locreg_descriptors),intent(in):: glr
!!type(locreg_descriptors),dimension(norb),intent(in):: llr
!!integer,dimension(2,0:nproc-1),intent(in):: istartend
!!integer,dimension(0:nproc-1),intent(in):: nsendcounts, nsenddspls, nrecvcounts, nrecvdspls
!!real(8),intent(in):: weightp
!!integer,dimension(ndimpsi),intent(out):: isendbuf, irecvbuf
!!integer,dimension(sum(nrecvcounts)),intent(out):: indexrecvorbital, iextract, iexpand
!!
!!! Local variables
!!integer:: i, j, iorb, iiorb, i1, i2, i3, ind, jproc, jproctarget, ii, ierr, jj
!!integer,dimension(:),allocatable:: nsend, indexsendorbital2, gridpoint_start, indexrecvorbital2
!!real(8),dimension(:,:,:),allocatable:: weight
!!integer,dimension(:),allocatable:: indexsendorbital, indexsendbuf, indexrecvbuf
!!
!!allocate(indexsendorbital(ndimpsi))
!!allocate(indexsendbuf(ndimpsi))
!!allocate(indexrecvbuf(sum(nrecvcounts)))
!!
!!allocate(weight(glr%is1:glr%ie1,glr%is2:glr%ie2,glr%is3:glr%ie3))
!!allocate(gridpoint_start((glr%ie1-glr%is1+1)*(glr%ie2-glr%is2+1)*(glr%ie3-glr%is3+1)))
!!gridpoint_start=-1
!!
!!  allocate(nsend(0:nproc-1))
!!  i=0
!!  nsend=0
!!  do iorb=1,norbp
!!      iiorb=isorb+iorb
!!      do i3=llr(iiorb)%is3,llr(iiorb)%ie3
!!          do i2=llr(iiorb)%is2,llr(iiorb)%ie2
!!              do i1=llr(iiorb)%is1,llr(iiorb)%ie1
!!                  ii = (i3-1)*(glr%ie1-glr%is1+1)*(glr%ie2-glr%is2+1) + (i2-1)*(glr%ie1-glr%is1+1) + i1
!!                  i=i+1
!!                  do jproc=0,nproc-1
!!                      if(ii>=istartend(1,jproc) .and. ii<=istartend(2,jproc)) then
!!                          jproctarget=jproc
!!                          exit
!!                      end if
!!                  end do
!!                  nsend(jproctarget)=nsend(jproctarget)+1
!!                  ind=nsenddspls(jproctarget)+nsend(jproctarget)
!!                  isendbuf(i)=ind
!!                  indexsendbuf(ind)=ii
!!                  indexsendorbital(i)=iiorb
!!                  !indexsendorbital(ind)=iiorb
!!              end do
!!          end do
!!      end do
!!  end do
!!
!!
!!  allocate(indexsendorbital2(ndimpsi))
!!  indexsendorbital2=indexsendorbital
!!  do i=1,ndimpsi
!!      ind=isendbuf(i)
!!      indexsendorbital(ind)=indexsendorbital2(i)
!!  end do
!!
!!  ! Inverse of isendbuf
!!  do i=1,ndimpsi
!!      do j=1,ndimpsi
!!          if(isendbuf(j)==i) then
!!              irecvbuf(i)=j
!!          end if
!!      end do
!!  end do
!!
!!
!!  !check
!!  do jproc=0,nproc-1
!!      if(nsend(jproc)/=nsendcounts(jproc)) stop 'nsend(jproc)/=nsendcounts(jproc)'
!!  end do
!!
!!
!!
!!  ! Communicate indexsendbuf
!!  !allocate(indexrecvbuf(sum(nrecvcounts)))
!!  call mpi_alltoallv(indexsendbuf, nsendcounts, nsenddspls, mpi_integer, indexrecvbuf, &
!!       nrecvcounts, nrecvdspls, mpi_integer, mpi_comm_world, ierr)
!!
!!  ! Communicate indexsendorbitals
!!  !allocate(indexrecvorbital(sum(nrecvcounts)))
!!  call mpi_alltoallv(indexsendorbital, nsendcounts, nsenddspls, mpi_integer, indexrecvorbital, &
!!       nrecvcounts, nrecvdspls, mpi_integer, mpi_comm_world, ierr)
!!
!!
!!
!!  call get_gridpoint_start(iproc, nproc, norb, glr, llr, nrecvcounts, indexrecvbuf, weight, gridpoint_start)
!!
!!
!!
!!  if(maxval(gridpoint_start)>sum(nrecvcounts)) stop '1: maxval(gridpoint_start)>sum(nrecvcounts)'
!!  !write(*,*) 'maxval(gridpoint_start), sum(nrecvcounts)', maxval(gridpoint_start), sum(nrecvcounts)
!!  ! Rearrange the communicated data
!!  do i=1,sum(nrecvcounts)
!!      ii=indexrecvbuf(i)
!!      jj=ii-1
!!      i3=jj/((glr%ie1-glr%is1+1)*(glr%ie2-glr%is2+1))+1
!!      jj=jj-(i3-1)*(glr%ie1-glr%is1+1)*(glr%ie2-glr%is2+1)
!!      i2=jj/(glr%ie2-glr%is2+1)+1
!!      i1=jj-(i2-1)*(glr%ie2-glr%is2+1)+1
!!      if(weight(i1,i2,i3)==0.d0) stop 'weight is zero!'
!!      if(weight(i1,i2,i3)>0.d0) then
!!          if(gridpoint_start(ii)==0) then
!!              write(*,'(a,5i8)') 'DEBUG: iproc, jj, i1, i2, i3', iproc, jj, i1, i2, i3
!!              stop 'weight>0, but gridpoint_start(ii)==0'
!!          end if
!!      end if
!!
!!      !!if(gridpoint_start(ii)==0) then
!!      !!    if(weight(i1,i2,i3)>0.d0) stop 'weight>0, but gridpoint_start(ii)==0'
!!      !!end if
!!      ind=gridpoint_start(ii)
!!      if(ind==0) stop 'ind is zero!'
!!      iextract(i)=ind
!!      write(300+iproc,*) i, iextract(i)
!!      gridpoint_start(ii)=gridpoint_start(ii)+1  
!!  end do
!!  if(sum(iextract)/=nint(weightp*(weightp+1.d0)*.5d0)) stop 'sum(iextract)/=nint(weightp*(weightp+1.d0)*.5d0)'
!!  !write(*,*) 'AFTER: maxval(gridpoint_start), sum(nrecvcounts)', maxval(gridpoint_start), sum(nrecvcounts)
!!  if(maxval(iextract)>sum(nrecvcounts)) stop 'maxval(iextract)>sum(nrecvcounts)'
!!  if(minval(iextract)<1) stop 'minval(iextract)<1'
!!  write(*,*) 'minval(iextract)', minval(iextract)
!!  !if(maxval(gridpoint_start)>sum(nrecvcounts)) stop 'maxval(gridpoint_start)>sum(nrecvcounts)'
!!
!!
!!  ! Get the array to transfrom back the data
!!  !allocate(iexpand(sum(nrecvcounts)))
!!  do i=1,sum(nrecvcounts)
!!      do j=1,sum(nrecvcounts)
!!          if(iextract(j)==i) then
!!              iexpand(i)=j
!!          end if
!!      end do
!!  end do
!!  
!!
!!
!!  allocate(indexrecvorbital2(sum(nrecvcounts)))
!!  indexrecvorbital2=indexrecvorbital
!!  do i=1,sum(nrecvcounts)
!!      ind=iextract(i)
!!      indexrecvorbital(ind)=indexrecvorbital2(i)
!!  end do
!!
!!  if(minval(indexrecvorbital)<1) stop 'minval(indexrecvorbital)<1'
!!  if(maxval(indexrecvorbital)>norb) stop 'maxval(indexrecvorbital)>norb'
!!
!!
!!
!!
!!end subroutine get_switch_indices
!!
!!
!!
!!
!!subroutine get_gridpoint_start(iproc, nproc, norb, glr, llr, nrecvcounts, indexrecvbuf, weight, gridpoint_start)
!!use types
!!implicit none
!!
!!! Calling arguments
!!integer,intent(in):: iproc, nproc, norb
!!type(locreg_descriptors),intent(in):: glr
!!type(locreg_descriptors),dimension(norb),intent(in):: llr
!!integer,dimension(0:nproc-1),intent(in):: nrecvcounts
!!integer,dimension(sum(nrecvcounts)),intent(in):: indexrecvbuf
!!real(8),dimension(glr%is1:glr%ie1,glr%is2:glr%ie2,glr%is3:glr%ie3),intent(out):: weight
!!integer,dimension((glr%ie1-glr%is1+1)*(glr%ie2-glr%is2+1)*(glr%ie3-glr%is3+1)),intent(out):: gridpoint_start
!!
!!! Local variables
!!integer:: i, ii, jj, i1, i2, i3
!!
!!
!!  weight=0.d0
!!  do i=1,sum(nrecvcounts)
!!      ii=indexrecvbuf(i)
!!      jj=ii-1
!!      i3=jj/((glr%ie1-glr%is1+1)*(glr%ie2-glr%is2+1))+1
!!      jj=jj-(i3-1)*(glr%ie1-glr%is1+1)*(glr%ie2-glr%is2+1)
!!      i2=jj/(glr%ie2-glr%is2+1)+1
!!      i1=jj-(i2-1)*(glr%ie2-glr%is2+1)+1
!!      write(100+iproc,'(a,4i8)') 'ii, i1, i2, i3', ii, i1, i2, i3
!!      weight(i1,i2,i3)=weight(i1,i2,i3)+1.d0
!!  end do
!!  ii=1
!!  i=0
!!  gridpoint_start=0
!!  do i3=glr%is3,glr%ie3
!!      do i2=glr%is2,glr%ie2
!!          do i1=glr%is1,glr%ie1
!!              i=i+1
!!              if(weight(i1,i2,i3)>0.d0) then
!!                  gridpoint_start(i)=ii
!!                  ii=ii+nint(weight(i1,i2,i3))
!!              end if
!!          end do
!!      end do
!!  end do
!!  !! CHECK
!!  i=0
!!  do i3=glr%is3,glr%ie3
!!      do i2=glr%is2,glr%ie2
!!          do i1=glr%is1,glr%ie1
!!              i=i+1
!!              write(200+iproc,'(a,4i8,es13.5,i8)') 'i, i1, i2, i3, weight(i1,i2,i3), gridpoint_start(i)', i, i1, i2, i3, weight(i1,i2,i3), gridpoint_start(i)
!!              if(weight(i1,i2,i3)>0.d0) then
!!                  if(gridpoint_start(i)==0) stop 'FIRST CHECK: ERROR'
!!              end if
!!          end do
!!      end do
!!  end do
!!  write(*,*) 'maxval(gridpoint_start), sum(nrecvcounts), sum(weight)', maxval(gridpoint_start), sum(nrecvcounts), sum(weight)
!!
!!end subroutine get_gridpoint_start




subroutine check_gridpoint(nseg, n1, n2, noffset1, noffset2, noffset3, keyg, itarget1, itarget2, itarget3, found)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: nseg, n1, n2, noffset1, noffset2, noffset3, itarget1, itarget2, itarget3
integer,dimension(2,nseg),intent(in):: keyg
logical,intent(out):: found

! Local variables
integer:: j0, j1, ii, i1, i2, i3, i0, ii1, ii2, ii3, iseg, i


  !This could be optimized a lot...
  found=.false.
  loop_segments: do iseg=1,nseg
     j0=keyg(1,iseg)
     j1=keyg(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        ii1=i+noffset1
        ii2=i2+noffset2
        ii3=i3+noffset3
        if(ii1==itarget1 .and. ii2==itarget2 .and. ii3==itarget3) then
            found=.true.
            exit loop_segments
        end if
     enddo
  enddo loop_segments



end subroutine check_gridpoint
