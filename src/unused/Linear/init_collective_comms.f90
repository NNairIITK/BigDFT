!!subroutine determine_num_orbs_per_gridpoint(iproc, nproc, orbs, lzd, istartend_c, istartend_f, &
!!           istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f, &
!!           weightp_c, weightp_f, nptsp_c, nptsp_f, &
!!           norb_per_gridpoint_c, norb_per_gridpoint_f)
!!  use module_base
!!  use module_types
!!  implicit none
!!  
!!  ! Calling arguments
!!  integer,intent(in):: iproc, nproc, nptsp_c, nptsp_f, istartp_seg_c, iendp_seg_c, istartp_seg_f, iendp_seg_f
!!  type(orbitals_data),intent(in):: orbs
!!  type(local_zone_descriptors),intent(in):: lzd
!!  integer,dimension(2,0:nproc-1),intent(in):: istartend_c, istartend_f
!!  real(8),intent(in):: weightp_c, weightp_f
!!  integer,dimension(nptsp_c),intent(out):: norb_per_gridpoint_c
!!  integer,dimension(nptsp_f),intent(out):: norb_per_gridpoint_f
!!  
!!  ! Local variables
!!  integer:: ii, iiorb, i1, i2, i3, iipt, iorb, iii, npgp, iseg, jj, j0, j1, iitot, ilr, i, istart, iend, i0, istat, iall
!!  logical:: found, overlap_possible
!!  integer,dimension(:),allocatable:: iseg_start_c, iseg_start_f
!!  character(len=*),parameter:: subname='determine_num_orbs_per_gridpoint'
!!  real(8):: t1, t2, t1tot, t2tot, t_check_gridpoint
!!
!!  allocate(iseg_start_c(lzd%nlr), stat=istat)
!!  call memocc(istat, iseg_start_c, 'iseg_start_c', subname)
!!  allocate(iseg_start_f(lzd%nlr), stat=istat)
!!  call memocc(istat, iseg_start_f, 'iseg_start_f', subname)
!!
!!  iseg_start_c=1
!!  iseg_start_f=1
!!
!!  iitot=0
!!  iiorb=0
!!  iipt=0
!!t_check_gridpoint=0.d0
!!t1tot=mpi_wtime()
!!  !write(*,*) 'iproc, istartp_seg_c,iendp_seg_c', iproc, istartp_seg_c,iendp_seg_c
!!    !do iseg=1,lzd%glr%wfd%nseg_c
!!    do iseg=istartp_seg_c,iendp_seg_c
!!       jj=lzd%glr%wfd%keyvloc(iseg)
!!       j0=lzd%glr%wfd%keygloc(1,iseg)
!!       j1=lzd%glr%wfd%keygloc(2,iseg)
!!       ii=j0-1
!!       i3=ii/((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1))
!!       ii=ii-i3*(lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)
!!       i2=ii/(lzd%glr%d%n1+1)
!!       i0=ii-i2*(lzd%glr%d%n1+1)
!!       i1=i0+j1-j0
!!       do i=i0,i1
!!           !iitot=iitot+1
!!           iitot=jj+i-i0
!!           if(iitot>=istartend_c(1,iproc) .and. iitot<=istartend_c(2,iproc)) then
!!               !write(200+iproc,'(5i10)') iitot, iseg, iitot, jj, jj+i-i0
!!               iipt=iipt+1
!!               npgp=0
!!               do iorb=1,orbs%norb
!!                   ilr=orbs%inwhichlocreg(iorb)
!!                   ! Check whether this orbitals extends here
!!                   call check_grid_point_from_boxes(i, i2, i3, lzd%llr(ilr), overlap_possible)
!!                   if(.not. overlap_possible) then
!!                       found=.false.
!!                   else
!!                       t1=mpi_wtime()
!!                       call check_gridpoint(lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, &
!!                            lzd%llr(ilr)%ns1, lzd%llr(ilr)%ns2, lzd%llr(ilr)%ns3, lzd%llr(ilr)%wfd%keygloc, &
!!                            i, i2, i3, iseg_start_c(ilr), found)
!!                       t2=mpi_wtime()
!!                       t_check_gridpoint=t_check_gridpoint+t2-t1
!!                   end if
!!                   if(found) then
!!                       npgp=npgp+1
!!                       iiorb=iiorb+1
!!                   end if
!!               end do
!!               norb_per_gridpoint_c(iipt)=npgp
!!           end if
!!      end do
!!  end do
!!
!!  if(iipt/=nptsp_c) stop 'iipt/=nptsp_c'
!!  if(iiorb/=nint(weightp_c)) stop 'iiorb/=weightp_c'
!!
!!
!!
!!  iitot=0
!!  iiorb=0
!!  iipt=0
!!    istart=lzd%glr%wfd%nseg_c+min(1,lzd%glr%wfd%nseg_f)
!!    iend=istart+lzd%glr%wfd%nseg_f-1
!!    !do iseg=istart,iend
!!    do iseg=istartp_seg_f,iendp_seg_f
!!       jj=lzd%glr%wfd%keyvloc(iseg)
!!       j0=lzd%glr%wfd%keygloc(1,iseg)
!!       j1=lzd%glr%wfd%keygloc(2,iseg)
!!       ii=j0-1
!!       i3=ii/((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1))
!!       ii=ii-i3*(lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)
!!       i2=ii/(lzd%glr%d%n1+1)
!!       i0=ii-i2*(lzd%glr%d%n1+1)
!!       i1=i0+j1-j0
!!       do i=i0,i1
!!           !iitot=iitot+1
!!           iitot=jj+i-i0
!!           if(iitot>=istartend_f(1,iproc) .and. iitot<=istartend_f(2,iproc)) then
!!               iipt=iipt+1
!!               npgp=0
!!               do iorb=1,orbs%norb
!!                   ilr=orbs%inwhichlocreg(iorb)
!!                   ! Check whether this orbitals extends here
!!                   call check_grid_point_from_boxes(i, i2, i3, lzd%llr(ilr), overlap_possible)
!!                   if(.not. overlap_possible) then
!!                       found=.false.
!!                   else
!!                       iii=lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)
!!                       t1=mpi_wtime()
!!                       call check_gridpoint(lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, &
!!                            lzd%llr(ilr)%ns1, lzd%llr(ilr)%ns2, lzd%llr(ilr)%ns3, &
!!                            lzd%llr(ilr)%wfd%keygloc(1,iii), &
!!                            i, i2, i3, iseg_start_f(ilr), found)
!!                       t2=mpi_wtime()
!!                       t_check_gridpoint=t_check_gridpoint+t2-t1
!!                   end if
!!                   if(found) then
!!                       npgp=npgp+1
!!                       iiorb=iiorb+1
!!                   end if
!!               end do
!!               norb_per_gridpoint_f(iipt)=npgp
!!           end if
!!      end do
!!  end do
!!
!!  if(iipt/=nptsp_f) stop 'iipt/=nptsp_f'
!!  !!write(*,*) 'iiorb, weightp_f', iiorb, weightp_f
!!  if(iiorb/=nint(weightp_f)) stop 'iiorb/=weightp_f'
!!
!!
!!  iall=-product(shape(iseg_start_c))*kind(iseg_start_c)
!!  deallocate(iseg_start_c, stat=istat)
!!  call memocc(istat, iall, 'iseg_start_c', subname)
!!  iall=-product(shape(iseg_start_f))*kind(iseg_start_f)
!!  deallocate(iseg_start_f, stat=istat)
!!  call memocc(istat, iall, 'iseg_start_f', subname)
!!
!!t2tot=mpi_wtime()
!!!if(iproc==0) write(*,'(a,es14.5)') 'in sub determine_num_orbs_per_gridpoint: iproc, total time', t2tot-t1tot
!!!if(iproc==0) write(*,'(a,es14.5)') 'in sub determine_num_orbs_per_gridpoint: iproc, time for check_gridpoint', t_check_gridpoint
!!
!!end subroutine determine_num_orbs_per_gridpoint




!!
!!subroutine get_index_in_global(lr, itarget1, itarget2, itarget3, region, ind)
!!use module_base
!!use module_types
!!implicit none
!!
!!! Calling arguments
!!type(locreg_descriptors),intent(in):: lr
!!integer,intent(in):: itarget1, itarget2, itarget3
!!character(len=1),intent(in):: region
!!integer,intent(out):: ind
!!
!!! Local variables
!!integer:: iitot, iseg, j0, j1, ii, i1, i2, i3, i0, i, istart, iend, ii1, ii2, ii3
!!
!!
!! if(region=='c') then
!!    iitot=0
!!    loop_segments_c: do iseg=1,lr%wfd%nseg_c
!!       j0=lr%wfd%keygloc(1,iseg)
!!       j1=lr%wfd%keygloc(2,iseg)
!!       ii=j0-1
!!       i3=ii/((lr%d%n1+1)*(lr%d%n2+1))
!!       ii=ii-i3*(lr%d%n1+1)*(lr%d%n2+1)
!!       i2=ii/(lr%d%n1+1)
!!       i0=ii-i2*(lr%d%n1+1)
!!       i1=i0+j1-j0
!!       do i=i0,i1
!!          iitot=iitot+1
!!          ii1=i+lr%ns1
!!          ii2=i2+lr%ns2
!!          ii3=i3+lr%ns3
!!          if(ii1==itarget1 .and. ii2==itarget2 .and. ii3==itarget3) then
!!              ind=iitot
!!              exit loop_segments_c
!!          end if
!!       end do
!!    end do loop_segments_c
!!
!!  else if(region=='f') then
!!
!!    iitot=0
!!    istart=lr%wfd%nseg_c+min(1,lr%wfd%nseg_f)
!!    iend=istart+lr%wfd%nseg_f-1
!!    loop_segments_f: do iseg=istart,iend
!!       j0=lr%wfd%keygloc(1,iseg)
!!       j1=lr%wfd%keygloc(2,iseg)
!!       ii=j0-1
!!       i3=ii/((lr%d%n1+1)*(lr%d%n2+1))
!!       ii=ii-i3*(lr%d%n1+1)*(lr%d%n2+1)
!!       i2=ii/(lr%d%n1+1)
!!       i0=ii-i2*(lr%d%n1+1)
!!       i1=i0+j1-j0
!!       do i=i0,i1
!!          ii1=i+lr%ns1
!!          ii2=i2+lr%ns2
!!          ii3=i3+lr%ns3
!!          iitot=iitot+1
!!          if(ii1==itarget1 .and. ii2==itarget2 .and. ii3==itarget3) then
!!              ind=iitot
!!              exit loop_segments_f
!!          end if
!!       end do
!!    end do loop_segments_f
!!
!!else
!!    stop 'wrong region'
!!end if
!!
!!
!!
!!end subroutine get_index_in_global



