subroutine init_collective_comms_vectors(iproc, nproc, nlr, orbs, orbsig, mlr, collcom)
  use module_base
  use module_types
  use module_interfaces, except_this_one => init_collective_comms_vectors
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc, nlr
  type(orbitals_data),intent(in):: orbs, orbsig
  type(matrixLocalizationRegion),dimension(nlr),intent(in):: mlr
  type(collective_comms),intent(out):: collcom

  ! Local variables
  integer:: iorb, iiorb, ilr, iall, istat, ierr, ii
  integer,dimension(:,:),allocatable:: istartend
  real(8),dimension(:),allocatable:: weight
  real(8):: weightp, tt, weight_tot
  character(len=*),parameter:: subname='init_collective_comms_vectors'

  collcom%ndimpsi_c=0
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      ilr=orbs%inwhichlocreg(iiorb)
      collcom%ndimpsi_c=collcom%ndimpsi_c+mlr(ilr)%norbinlr
  end do
  collcom%ndimpsi_f=0

  allocate(weight(orbsig%norb), stat=istat)
  call memocc(istat, weight, 'weight', subname)

  call get_weights_vectors(iproc, nproc, orbs, nlr, mlr, orbsig%norb, weight, weight_tot)

  allocate(istartend(2,0:nproc-1), stat=istat)
  call memocc(istat, istartend, 'istartend', subname)
  call assign_weight_to_process_vectors(iproc, nproc, orbsig%norb, weight, weight_tot, istartend, weightp, collcom%nptsp_c)
  collcom%nptsp_f=0


  iall=-product(shape(weight))*kind(weight)
  deallocate(weight, stat=istat)
  call memocc(istat, iall, 'weight', subname)


  ! some checks
  if(nproc>1) then
      call mpi_allreduce(weightp, tt, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
  else
      tt=weightp
  end if
  !!write(*,*) 'tt, weight_c_tot', tt, weight_c_tot
  if(tt/=weight_tot) stop 'wrong partition of coarse weights'
  if(nproc>1) then
      call mpi_allreduce(collcom%nptsp_c, ii, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr)
  else
      ii=collcom%nptsp_c
  end if
  if(ii/=orbsig%norb) stop 'wrong partition of orbitals'


  allocate(collcom%norb_per_gridpoint_c(collcom%nptsp_c), stat=istat)
  call memocc(istat, collcom%norb_per_gridpoint_c, 'collcom%norb_per_gridpoint_c', subname)
  call determine_num_orbs_per_gridpoint_vectors(iproc, nproc, orbsig%norb, nlr, collcom%nptsp_c, orbs, &
           istartend, mlr, weightp, collcom%norb_per_gridpoint_c)
  allocate(collcom%norb_per_gridpoint_f(collcom%nptsp_f), stat=istat)
  call memocc(istat, collcom%norb_per_gridpoint_f, 'collcom%norb_per_gridpoint_f', subname)
  collcom%norb_per_gridpoint_f=0

  allocate(collcom%nsendcounts_c(0:nproc-1), stat=istat)
  call memocc(istat, collcom%nsendcounts_c, 'collcom%nsendcounts_c', subname)
  allocate(collcom%nsenddspls_c(0:nproc-1), stat=istat)
  call memocc(istat, collcom%nsenddspls_c, 'collcom%nsenddspls_c', subname)
  allocate(collcom%nrecvcounts_c(0:nproc-1), stat=istat)
  call memocc(istat, collcom%nrecvcounts_c, 'collcom%nrecvcounts_c', subname)
  allocate(collcom%nrecvdspls_c(0:nproc-1), stat=istat)
  call memocc(istat, collcom%nrecvdspls_c, 'collcom%nrecvdspls_c', subname)

  call determine_communication_arrays_vectors(iproc, nproc, nlr, orbs, mlr, istartend, weightp, &
           collcom%nsendcounts_c, collcom%nsenddspls_c, collcom%nrecvcounts_c, collcom%nrecvdspls_c)

  allocate(collcom%nsendcounts_f(0:nproc-1), stat=istat)
  call memocc(istat, collcom%nsendcounts_f, 'collcom%nsendcounts_f', subname)
  allocate(collcom%nsenddspls_f(0:nproc-1), stat=istat)
  call memocc(istat, collcom%nsenddspls_f, 'collcom%nsenddspls_f', subname)
  allocate(collcom%nrecvcounts_f(0:nproc-1), stat=istat)
  call memocc(istat, collcom%nrecvcounts_f, 'collcom%nrecvcounts_f', subname)
  allocate(collcom%nrecvdspls_f(0:nproc-1), stat=istat)
  call memocc(istat, collcom%nrecvdspls_f, 'collcom%nrecvdspls_f', subname)
  collcom%nsendcounts_f=0
  collcom%nsenddspls_f=0
  collcom%nrecvcounts_f=0
  collcom%nrecvdspls_f=0

  allocate(collcom%irecvbuf_c(collcom%ndimpsi_c), stat=istat)
  call memocc(istat, collcom%irecvbuf_c, 'collcom%irecvbuf_c', subname)
  allocate(collcom%indexrecvorbital_c(sum(collcom%nrecvcounts_c)), stat=istat)
  call memocc(istat, collcom%indexrecvorbital_c, 'collcom%indexrecvorbital_c', subname)
  allocate(collcom%iextract_c(sum(collcom%nrecvcounts_c)), stat=istat)
  call memocc(istat, collcom%iextract_c, 'collcom%iextract_c', subname)
  allocate(collcom%iexpand_c(sum(collcom%nrecvcounts_c)), stat=istat)
  call memocc(istat, collcom%iexpand_c, 'collcom%iexpand_c', subname)
  allocate(collcom%isendbuf_c(collcom%ndimpsi_c), stat=istat)
  call memocc(istat, collcom%isendbuf_c, 'collcom%isendbuf_c', subname)


  call get_switch_indices_vectors(iproc, nproc, nlr, orbsig%norb, collcom%ndimpsi_c, orbs, mlr, &
           istartend, collcom%nsendcounts_c, collcom%nsenddspls_c, collcom%nrecvcounts_c, collcom%nrecvdspls_c, weightp, &
           collcom%isendbuf_c, collcom%irecvbuf_c, collcom%indexrecvorbital_c, collcom%iextract_c, collcom%iexpand_c)

  allocate(collcom%irecvbuf_f(collcom%ndimpsi_f), stat=istat)
  call memocc(istat, collcom%irecvbuf_f, 'collcom%irecvbuf_f', subname)
  allocate(collcom%indexrecvorbital_f(sum(collcom%nrecvcounts_f)), stat=istat)
  call memocc(istat, collcom%indexrecvorbital_f, 'collcom%indexrecvorbital_f', subname)
  allocate(collcom%iextract_f(sum(collcom%nrecvcounts_f)), stat=istat)
  call memocc(istat, collcom%iextract_f, 'collcom%iextract_f', subname)
  allocate(collcom%iexpand_f(sum(collcom%nrecvcounts_f)), stat=istat)
  call memocc(istat, collcom%iexpand_f, 'collcom%iexpand_f', subname)
  allocate(collcom%isendbuf_f(collcom%ndimpsi_f), stat=istat)
  call memocc(istat, collcom%isendbuf_f, 'collcom%isendbuf_f', subname)
  collcom%irecvbuf_f=0
  collcom%indexrecvorbital_f=0
  collcom%iextract_f=0
  collcom%iexpand_f=0
  collcom%isendbuf_f=0


  iall=-product(shape(istartend))*kind(istartend)
  deallocate(istartend, stat=istat)
  call memocc(istat, iall, 'istartend', subname)

end subroutine init_collective_comms_vectors


subroutine get_weights_vectors(iproc, nproc, orbs, nlr, mlr, norbig, weight, weight_tot)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc, norbig, nlr
  type(orbitals_data),intent(in):: orbs
  type(matrixLocalizationRegion),dimension(nlr),intent(in):: mlr
  real(8),dimension(norbig),intent(out):: weight
  real(8),intent(out):: weight_tot

  ! Local variables
  integer:: iorb, iiorb, ilr, jorb, ind, ierr

  weight=0.d0
  weight_tot=0.d0
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      ilr=orbs%inwhichlocreg(iiorb)
      do jorb=1,mlr(ilr)%norbinlr
          ind=mlr(ilr)%indexInGlobal(jorb)
          weight(ind)=weight(ind)+1.d0
          weight_tot=weight_tot+1.d0
      end do
  end do

  if(nproc>1) then
      call mpiallred(weight_tot, 1, mpi_sum, mpi_comm_world, ierr)
      call mpiallred(weight(1), norbig,  mpi_sum, mpi_comm_world, ierr)
  end if


end subroutine get_weights_vectors





subroutine assign_weight_to_process_vectors(iproc, nproc, norbig, weight, weight_tot, istartend, weightp, nptsp)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc, norbig
  real(8),dimension(norbig),intent(in):: weight
  real(8),intent(in):: weight_tot
  integer,dimension(2,0:nproc-1),intent(out):: istartend
  real(8),intent(out):: weightp
  integer,intent(out):: nptsp
  
  ! Local variables
  integer:: jproc, iitot, ii2, iorb, ii, ierr, jprocdone
  real(8):: tt, tt2, weight_ideal
   
  weight_ideal=weight_tot/dble(nproc)

  jproc=0
  tt=0.d0
  tt2=0.d0
  iitot=0
  ii2=0
  !!iiseg=1
  jprocdone=-1
  loop_norbig: do iorb=1,norbig
      tt=tt+weight(iorb) 
      iitot=iitot+1
      if(tt>=weight_ideal .or. iorb==norbig .and. jproc<=nproc-2) then
          if(tt==weight_tot .and. jproc==nproc-1) then
              ! this process also has to take the remaining points, even if they have no weight
              iitot=norbig-ii2
          end if
          if(iproc==jproc) then
              weightp=tt
              nptsp=iitot
              !istartp_seg_c=iiseg
              !iendp_seg_c=iseg
          end if
          istartend(1,jproc)=ii2+1
          istartend(2,jproc)=min(istartend(1,jproc)+iitot-1,norbig)
          tt2=tt2+tt
          tt=0.d0
          ii2=ii2+iitot
          iitot=0
          jproc=jproc+1
          if(ii2>=norbig) then
              ! everything is distributed
              jprocdone=jproc
              exit loop_norbig
          end if
          !!iiseg=iseg
      end if
  end do loop_norbig

  if(jprocdone>0) then
       do jproc=jprocdone,nproc-1
          ! these processes do nothing
          istartend(1,jproc)=norbig+1
          istartend(2,jproc)=norbig
          if(iproc==jproc) then
              weightp=0.d0
              nptsp=0
          end if
      end do
  else
      if(iproc==nproc-1) then
          ! Take the rest
          weightp=weight_tot-tt2
          nptsp=norbig-ii2
          !istartp_seg=iiseg
          !iendp_seg=lzd%glr%wfd%nseg_c
      end if
      istartend(1,nproc-1)=ii2+1
      istartend(2,nproc-1)=istartend(1,nproc-1)+iitot-1
  end if

  ! some check
  ii=istartend(2,iproc)-istartend(1,iproc)+1
  if(nproc>1) call mpiallred(ii, 1, mpi_sum, mpi_comm_world, ierr)
  if(ii/=norbig) stop 'ii/=norbig'


end subroutine assign_weight_to_process_vectors



subroutine determine_num_orbs_per_gridpoint_vectors(iproc, nproc, norbig, nlr, nptsp, orbs, &
           istartend, mlr, weightp, norb_per_gridpoint)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc, norbig, nlr, nptsp
  type(orbitals_data),intent(in):: orbs
  integer,dimension(2,0:nproc-1),intent(in):: istartend
  type(matrixLocalizationRegion),dimension(nlr),intent(in):: mlr
  real(8),intent(in):: weightp
  integer,dimension(nptsp),intent(out):: norb_per_gridpoint

  ! Local variables
  integer:: jjtot, jjorb, jjpt, iorb, npgp, jorb, jlr, korb
  logical:: found

  jjtot=0
  jjorb=0
  jjpt=0
  do iorb=1,norbig
      if(iorb>=istartend(1,iproc) .and. iorb<=istartend(2,iproc)) then
          jjpt=jjpt+1
          npgp=0
          do jorb=1,orbs%norb
              jlr=orbs%inwhichlocreg(jorb)
              ! Check whether this orbitals extends here
              !!call check_grid_point_from_boxes(i, i2, i3, lzd%llr(ilr), overlap_possible)
              !!if(.not. overlap_possible) then
              !!    found=.false.
              !!else
              !!    call check_gridpoint(lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, &
              !!         lzd%llr(ilr)%ns1, lzd%llr(ilr)%ns2, lzd%llr(ilr)%ns3, lzd%llr(ilr)%wfd%keygloc, &
              !!         i, i2, i3, iseg_start_c(ilr), found)
              !!end if
              found=.false.
              do korb=1,mlr(jlr)%norbinlr
                  if(mlr(jlr)%indexInGlobal(korb)==iorb) then
                      found=.true.
                      exit
                  end if
              end do
              if(found) then
                  npgp=npgp+1
                  jjorb=jjorb+1
              end if
          end do
          norb_per_gridpoint(jjpt)=npgp
      end if
  end do

  if(jjpt/=nptsp) stop 'jjpt/=nptsp'
  if(jjorb/=nint(weightp)) stop 'jjorb/=weightp'


end subroutine determine_num_orbs_per_gridpoint_vectors




subroutine determine_communication_arrays_vectors(iproc, nproc, nlr, orbs, mlr, istartend, weightp, &
           nsendcounts, nsenddspls, nrecvcounts, nrecvdspls)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc, nlr
  type(orbitals_data),intent(in):: orbs
  type(matrixLocalizationRegion),dimension(nlr),intent(in):: mlr
  integer,dimension(2,0:nproc-1),intent(in):: istartend
  real(8),intent(in):: weightp
  integer,dimension(0:nproc-1),intent(out):: nsendcounts, nsenddspls, nrecvcounts, nrecvdspls

  ! Local variables
  integer:: iorb, iiorb, ilr, jorb, ind, jproctarget, ii, istat, iall, ierr, jproc
  integer,dimension(:),allocatable:: nsendcounts_tmp, nsenddspls_tmp, nrecvcounts_tmp, nrecvdspls_tmp
  character(len=*),parameter:: subname='determine_communication_arrays_vectors'

  ! Determine values for mpi_alltoallv
  ! first nsendcounts
  nsendcounts=0
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      ilr=orbs%inwhichlocreg(iiorb)
      do jorb=1,mlr(ilr)%norbinlr
          ind=mlr(ilr)%indexInGlobal(jorb)
          jproctarget=-1
          do jproc=0,nproc-1
              if(ind>=istartend(1,jproc) .and. ind<=istartend(2,jproc)) then
                  jproctarget=jproc
                  exit
              end if
          end do
          nsendcounts(jproctarget)=nsendcounts(jproctarget)+1
      end do
  end do



  ! The first check is to make sure that there is no stop in case this process has no orbitals (in which case
  ! orbs%npsidim_orbs is 1 and not 0 as assumed by the check)
  ii=0
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      ilr=orbs%inwhichlocreg(iiorb)
      ii=ii+mlr(ilr)%norbinlr
  end do
  if(sum(nsendcounts)/=ii) stop 'sum(nsendcounts)/=ii'

  
  ! now nsenddspls
  nsenddspls(0)=0
  do jproc=1,nproc-1
      nsenddspls(jproc)=nsenddspls(jproc-1)+nsendcounts(jproc-1)
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
      call mpi_alltoallv(nsendcounts, nsendcounts_tmp, nsenddspls_tmp, mpi_integer, nrecvcounts, &
           nrecvcounts_tmp, nrecvdspls_tmp, mpi_integer, mpi_comm_world, ierr)
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

  ! now recvdspls
  nrecvdspls(0)=0
  do jproc=1,nproc-1
      nrecvdspls(jproc)=nrecvdspls(jproc-1)+nrecvcounts(jproc-1)
  end do

  !write(*,*) 'sum(nrecvcounts), nint(weightp)', sum(nrecvcounts), nint(weightp)
  !!write(*,*) 'sum(nrecvcounts_f), nint(weightp_f)', sum(nrecvcounts_f), nint(weightp_f)
  if(sum(nrecvcounts)/=nint(weightp)) stop 'sum(nrecvcounts)/=nint(nweightp)'

end subroutine determine_communication_arrays_vectors



subroutine get_switch_indices_vectors(iproc, nproc, nlr, norbig, ndimvec, orbs, mlr, &
           istartend, nsendcounts, nsenddspls, nrecvcounts, nrecvdspls, weightp, &
           isendbuf, irecvbuf, indexrecvorbital, iextract, iexpand)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc, nlr, norbig, ndimvec
  type(orbitals_data),intent(in):: orbs
  type(matrixLocalizationRegion),dimension(nlr),intent(in):: mlr
  integer,dimension(2,0:nproc-1),intent(in):: istartend
  integer,dimension(0:nproc-1),intent(in):: nsendcounts, nsenddspls, nrecvcounts, nrecvdspls
  real(8),intent(in):: weightp
  integer,dimension(ndimvec),intent(out):: isendbuf, irecvbuf
  integer,dimension(sum(nrecvcounts)),intent(out):: indexrecvorbital, iextract, iexpand

  ! Local variables
  integer:: iitot, iorb, iiorb, ilr, jorb, indglob, jproc, jproctarget, ind, iall, istat, ierr, i, ii
  integer,dimension(:),allocatable:: indexsendorbital, indexsendbuf, indexrecvbuf
  integer,dimension(:),allocatable:: gridpoint_start, nsend, indexsendorbital2, indexrecvorbital2
  real(8),dimension(:),allocatable:: weight 
  character(len=*),parameter:: subname='get_switch_indices_vectors'


allocate(indexsendorbital(ndimvec), stat=istat)
call memocc(istat, indexsendorbital, 'indexsendorbital', subname)
allocate(indexsendbuf(ndimvec), stat=istat)
call memocc(istat, indexsendbuf, 'indexsendbuf', subname)
allocate(indexrecvbuf(sum(nrecvcounts)), stat=istat)
call memocc(istat, indexrecvbuf, 'indexrecvbuf', subname)


allocate(weight(norbig), stat=istat)
call memocc(istat, weight, 'weight', subname)
allocate(gridpoint_start(norbig), stat=istat)
call memocc(istat, gridpoint_start, 'gridpoint_start', subname)
gridpoint_start=-1


  allocate(nsend(0:nproc-1), stat=istat)
  call memocc(istat, nsend, 'nsend', subname)

  iitot=0
  nsend=0
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      ilr=orbs%inwhichlocreg(iiorb)
      do jorb=1,mlr(ilr)%norbinlr
          indglob=mlr(ilr)%indexInGlobal(jorb)
          iitot=iitot+1
          do jproc=0,nproc-1
              if(indglob>=istartend(1,jproc) .and. indglob<=istartend(2,jproc)) then
                  jproctarget=jproc
                  exit
              end if
          end do
          nsend(jproctarget)=nsend(jproctarget)+1
          ind=nsenddspls(jproctarget)+nsend(jproctarget)
          isendbuf(iitot)=ind
          indexsendbuf(ind)=indglob
          indexsendorbital(iitot)=iiorb
      end do
  end do

  !!!ndimvec=0
  !!!do iorb=1,orbs%norbp
  !!!    iiorb=orbs%isorb+iorb
  !!!    ilr=orbs%inwhichlocreg(iiorb)
  !!!    ndimvec=ndimvec+mlr(ilr)%norbinlr
  !!!end do
  if(iitot/=ndimvec) stop 'iitot/=ndimvec'

  !check
  do jproc=0,nproc-1
      if(nsend(jproc)/=nsendcounts(jproc)) stop 'nsend(jproc)/=nsendcounts(jproc)'
  end do






  allocate(indexsendorbital2(ndimvec), stat=istat)
  call memocc(istat, indexsendorbital2, 'indexsendorbital2', subname)
  indexsendorbital2=indexsendorbital
  do i=1,ndimvec
      ind=isendbuf(i)
      indexsendorbital(ind)=indexsendorbital2(i)
  end do
  ! Inverse of isendbuf
  call get_reverse_indices(ndimvec, isendbuf, irecvbuf)
  iall=-product(shape(indexsendorbital2))*kind(indexsendorbital2)
  deallocate(indexsendorbital2, stat=istat)
  call memocc(istat, iall, 'indexsendorbital2', subname)





  if(nproc>1) then
      ! Communicate indexsendbuf
      call mpi_alltoallv(indexsendbuf, nsendcounts, nsenddspls, mpi_integer, indexrecvbuf, &
           nrecvcounts, nrecvdspls, mpi_integer, mpi_comm_world, ierr)
      ! Communicate indexsendorbitals
      call mpi_alltoallv(indexsendorbital, nsendcounts, nsenddspls, mpi_integer, indexrecvorbital, &
           nrecvcounts, nrecvdspls, mpi_integer, mpi_comm_world, ierr)

   else
       indexrecvbuf=indexsendbuf
       indexrecvorbital=indexsendorbital
   end if



  !!call get_gridpoint_start(iproc, nproc, norb, glr, llr, nrecvcounts, indexrecvbuf, weight, gridpoint_start)
  !!call get_gridpoint_start(iproc, nproc, lzd, nrecvcounts_c, nrecvcounts_f, indexrecvbuf_c, indexrecvbuf_f, &
  !!          weight_c, weight_f, gridpoint_start_c, gridpoint_start_f)
  call get_gridpoint_start_vectors(iproc, nproc, norbig, nrecvcounts, indexrecvbuf, weight, gridpoint_start)



  if(maxval(gridpoint_start)>sum(nrecvcounts)) stop '1: maxval(gridpoint_start)>sum(nrecvcounts)'
  ! Rearrange the communicated data
  do i=1,sum(nrecvcounts)
      ii=indexrecvbuf(i)
      ind=gridpoint_start(ii)
      !!if(ind==0) stop 'ind is zero!'
      iextract(i)=ind
      gridpoint_start(ii)=gridpoint_start(ii)+1  
  end do
  !!write(*,'(a,2i12)') 'sum(iextract_c), nint(weightp_c*(weightp_c+1.d0)*.5d0)', sum(iextract_c), nint(weightp_c*(weightp_c+1.d0)*.5d0)
  !!if(sum(iextract_c)/=nint(weightp_c*(weightp_c+1.d0)*.5d0)) stop 'sum(iextract_c)/=nint(weightp_c*(weightp_c+1.d0)*.5d0)'
  if(maxval(iextract)>sum(nrecvcounts)) stop 'maxval(iextract)>sum(nrecvcounts)'
  if(minval(iextract)<1) stop 'minval(iextract)<1'




  ! Get the array to transfrom back the data
  call get_reverse_indices(sum(nrecvcounts), iextract, iexpand)




  allocate(indexrecvorbital2(sum(nrecvcounts)), stat=istat)
  call memocc(istat, indexrecvorbital2, 'indexrecvorbital2', subname)
  indexrecvorbital2=indexrecvorbital
  do i=1,sum(nrecvcounts)
      ind=iextract(i)
      indexrecvorbital(ind)=indexrecvorbital2(i)
  end do
  iall=-product(shape(indexrecvorbital2))*kind(indexrecvorbital2)
  deallocate(indexrecvorbital2, stat=istat)
  call memocc(istat, iall, 'indexrecvorbital2', subname)



  if(minval(indexrecvorbital)<1) stop 'minval(indexrecvorbital)<1'
  if(maxval(indexrecvorbital)>orbs%norb) stop 'maxval(indexrecvorbital)>orbs%norb'



  iall=-product(shape(indexsendorbital))*kind(indexsendorbital)
  deallocate(indexsendorbital, stat=istat)
  call memocc(istat, iall, 'indexsendorbital', subname)
  iall=-product(shape(indexsendbuf))*kind(indexsendbuf)
  deallocate(indexsendbuf, stat=istat)
  call memocc(istat, iall, 'indexsendbuf', subname)
  iall=-product(shape(indexrecvbuf))*kind(indexrecvbuf)
  deallocate(indexrecvbuf, stat=istat)
  call memocc(istat, iall, 'indexrecvbuf', subname)

  iall=-product(shape(weight))*kind(weight)
  deallocate(weight, stat=istat)
  call memocc(istat, iall, 'weight', subname)

  iall=-product(shape(gridpoint_start))*kind(gridpoint_start)
  deallocate(gridpoint_start, stat=istat)
  call memocc(istat, iall, 'gridpoint_start', subname)

  iall=-product(shape(nsend))*kind(nsend)
  deallocate(nsend, stat=istat)
  call memocc(istat, iall, 'nsend', subname)


end subroutine get_switch_indices_vectors



subroutine get_gridpoint_start_vectors(iproc, nproc, norbig, nrecvcounts, indexrecvbuf, weight, gridpoint_start)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc, norbig
  integer,dimension(0:nproc-1),intent(in):: nrecvcounts
  integer,dimension(sum(nrecvcounts)),intent(in):: indexrecvbuf
  real(8),dimension(norbig),intent(out):: weight
  integer,dimension(norbig),intent(out):: gridpoint_start

  ! Local variables
  integer:: i, ii, iorb 

  weight=0.d0
  do i=1,sum(nrecvcounts)
      ii=indexrecvbuf(i)
      weight(ii)=weight(ii)+1.d0
  end do


  ii=1
  gridpoint_start=0
  do iorb=1,norbig
      if(weight(iorb)>0.d0) then
          gridpoint_start(iorb)=ii
          ii=ii+nint(weight(iorb))
      end if
  end do

  !! CHECK
  do iorb=1,norbig
      if(weight(iorb)>0.d0) then
          if(gridpoint_start(iorb)==0) stop 'FIRST CHECK: ERROR'
      end if
  end do

end subroutine get_gridpoint_start_vectors
