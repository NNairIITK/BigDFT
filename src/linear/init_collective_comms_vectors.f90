subroutine get_weights_vectors(orbs, nlr, mlr, norbig, weight)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: norbig, nlr
  type(orbitals_data),intent(in):: orbs
  type(matrixLocalizationRegion),dimension(nlr),intent(in):: mlr
  real(8),dimension(norbig),intent(out):: weight

  ! Local variables
  integer:: iorb, iiorb, ilr, jorb, ind

  weight=0.d0
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      ilr=orbs%inwhichlocreg(iiorb)
      do jorb=1,mlr(ilr)%norbinlr
          ind=mlr(ilr)%indexInGlobal(jorb)
          weight(ind)=weight(ind)+1.d0
      end do
  end do

end subroutine get_weights_vectors





subroutine assign_weight_to_process(iproc, nproc, norbig, weight, weight_tot, istartend, weightp, nptsp)
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
  integer:: jproc, iitot, ii2, iorb, istartp, iendp, ii, ierr
  real(8):: tt, tt2, weight_ideal
   
  weight_ideal=weight_tot/dble(nproc)

  jproc=0
  tt=0.d0
  tt2=0.d0
  iitot=0
  ii2=0
  !!iiseg=1
  do iorb=1,norbig
      tt=tt+weight(iorb) 
      iitot=iitot+1
      if(tt>weight_ideal) then
          if(iproc==jproc) then
              weightp=tt
              nptsp=iitot
              istartp=ii2+1
              iendp=istartp+iitot-1
              !istartp_seg_c=iiseg
              !iendp_seg_c=iseg
          end if
          istartend(1,jproc)=ii2+1
          istartend(2,jproc)=istartend(1,jproc)+iitot-1
          tt2=tt2+tt
          tt=0.d0
          ii2=ii2+iitot
          iitot=0
          jproc=jproc+1
          !!iiseg=iseg
      end if
  end do
  if(iproc==nproc-1) then
      ! Take the rest
      istartp=ii2+1
      iendp=istartp+iitot-1
      weightp=weight_tot-tt2
      nptsp=norbig-ii2
      !istartp_seg=iiseg
      !iendp_seg=lzd%glr%wfd%nseg_c
  end if
  istartend(1,nproc-1)=ii2+1
  istartend(2,nproc-1)=istartend(1,nproc-1)+iitot-1

  ! some check
  ii=istartend(2,iproc)-istartend(1,iproc)+1
  if(nproc>1) call mpiallred(ii, 1, mpi_sum, mpi_comm_world, ierr)
  if(ii/=norbig) stop 'ii/=norbig'


end subroutine assign_weight_to_process



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
