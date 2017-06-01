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



!!subroutine calculate_pulay_overlap(iproc, nproc, orbs1, orbs2, collcom1, collcom2, psit_c1, psit_c2, psit_f1, psit_f2, ovrlp)
!!  use module_base
!!  use module_types
!!  implicit none
!!  
!!  ! Calling arguments
!!  integer,intent(in) :: iproc, nproc
!!  type(orbitals_data),intent(in) :: orbs1, orbs2
!!  type(comms_linear),intent(in) :: collcom1, collcom2
!!  real(kind=8),dimension(collcom1%ndimind_c),intent(in) :: psit_c1
!!  real(kind=8),dimension(collcom2%ndimind_c),intent(in) :: psit_c2
!!  real(kind=8),dimension(7*collcom1%ndimind_f),intent(in) :: psit_f1
!!  real(kind=8),dimension(7*collcom2%ndimind_f),intent(in) :: psit_f2
!!  real(kind=8),dimension(orbs1%norb,orbs2%norb),intent(out) :: ovrlp
!!  
!!  ! Local variables
!!  integer :: i0, j0, ipt, ii, iiorb, j, jj, jjorb, i, ierr  
!!
!!  call timing(iproc,'ovrlptransComp','ON') !lr408t
!!  call f_zero(ovrlp)
!!  if(collcom1%nptsp_c/=collcom2%nptsp_c) then
!!      write(*,'(a,i0,a)') 'ERROR on process ',iproc,': collcom1%nptsp_c/=collcom2%nptsp_c'
!!      stop
!!  end if
!!  if(collcom1%nptsp_f/=collcom2%nptsp_f) then
!!      write(*,'(a,i0,a)') 'ERROR on process ',iproc,': collcom1%nptsp_f/=collcom2%nptsp_f'
!!      stop
!!  end if
!!
!!  i0=0
!!  j0=0
!!  do ipt=1,collcom1%nptsp_c 
!!      ii=collcom1%norb_per_gridpoint_c(ipt)
!!      jj=collcom2%norb_per_gridpoint_c(ipt)
!!      do i=1,ii
!!          iiorb=collcom1%indexrecvorbital_c(i0+i)
!!          do j=1,jj
!!              jjorb=collcom2%indexrecvorbital_c(j0+j)
!!              ovrlp(iiorb,jjorb)=ovrlp(iiorb,jjorb)+psit_c1(i0+i)*psit_c2(j0+j)
!!          end do
!!      end do
!!      i0=i0+ii
!!      j0=j0+jj
!!  end do
!!
!!  i0=0
!!  j0=0
!!  do ipt=1,collcom1%nptsp_f 
!!      ii=collcom1%norb_per_gridpoint_f(ipt)
!!      jj=collcom2%norb_per_gridpoint_f(ipt)
!!      do i=1,ii
!!          iiorb=collcom1%indexrecvorbital_f(i0+i)
!!          do j=1,jj
!!              jjorb=collcom2%indexrecvorbital_f(j0+j)
!!              ovrlp(iiorb,jjorb)=ovrlp(iiorb,jjorb)+psit_f1(7*(i0+i)-6)*psit_f2(7*(j0+j)-6)
!!              ovrlp(iiorb,jjorb)=ovrlp(iiorb,jjorb)+psit_f1(7*(i0+i)-5)*psit_f2(7*(j0+j)-5)
!!              ovrlp(iiorb,jjorb)=ovrlp(iiorb,jjorb)+psit_f1(7*(i0+i)-4)*psit_f2(7*(j0+j)-4)
!!              ovrlp(iiorb,jjorb)=ovrlp(iiorb,jjorb)+psit_f1(7*(i0+i)-3)*psit_f2(7*(j0+j)-3)
!!              ovrlp(iiorb,jjorb)=ovrlp(iiorb,jjorb)+psit_f1(7*(i0+i)-2)*psit_f2(7*(j0+j)-2)
!!              ovrlp(iiorb,jjorb)=ovrlp(iiorb,jjorb)+psit_f1(7*(i0+i)-1)*psit_f2(7*(j0+j)-1)
!!              ovrlp(iiorb,jjorb)=ovrlp(iiorb,jjorb)+psit_f1(7*(i0+i)-0)*psit_f2(7*(j0+j)-0)
!!          end do
!!      end do
!!      i0=i0+ii
!!      j0=j0+jj
!!  end do
!!
!!  call timing(iproc,'ovrlptransComp','OF') !lr408t
!!
!!  call timing(iproc,'ovrlptransComm','ON') !lr408t
!!
!!  if(nproc>1) then
!!      call mpiallred(ovrlp, mpi_sum, comm=bigdft_mpi%mpi_comm)
!!  end if
!!  call timing(iproc,'ovrlptransComm','OF') !lr408t
!!end subroutine calculate_pulay_overlap





!!subroutine check_grid_point_from_boxes(i1, i2, i3, lr, overlap_possible)
!!  use module_base
!!  use module_types
!!  implicit none
!!  
!!  ! Calling arguments
!!  integer,intent(in) :: i1, i2, i3
!!  type(locreg_descriptors),intent(in) :: lr  
!!  logical,intent(out) :: overlap_possible
!!
!!  ! Local variables
!!  logical :: ovrlpx, ovrlpy, ovrlpz
!!  
!!  ovrlpx = (i1>=lr%ns1 .and. i1<=lr%ns1+lr%d%n1)
!!  ovrlpy = (i2>=lr%ns2 .and. i2<=lr%ns2+lr%d%n2)
!!  ovrlpz = (i3>=lr%ns3 .and. i3<=lr%ns3+lr%d%n3)
!!  if(ovrlpx .and. ovrlpy .and. ovrlpz) then
!!      overlap_possible=.true.
!!  else
!!      overlap_possible=.true.
!!  end if
!!
!!end subroutine check_grid_point_from_boxes

!!subroutine get_reverse_indices(n, indices, reverse_indices)
!!  use module_base
!!  implicit none
!!  
!!  ! Calling arguments
!!  integer,intent(in) :: n
!!  integer,dimension(n),intent(in) :: indices
!!  integer,dimension(n),intent(out) :: reverse_indices
!!
!!  ! Local variables
!!  integer :: i, j, m, j0, j1, j2, j3
!!
!!  !$omp parallel default(private) &
!!  !$omp shared(n, m, indices, reverse_indices)
!!
!!  m=mod(n,4)
!!  if (m/=0) then
!!      do i=1,m
!!          j=indices(i)
!!          reverse_indices(j)=i
!!      end do
!!  end if
!!
!!  !$omp do
!!  do i=m+1,n,4
!!      j0=indices(i+0)
!!      reverse_indices(j0)=i+0
!!      j1=indices(i+1)
!!      reverse_indices(j1)=i+1
!!      j2=indices(i+2)
!!      reverse_indices(j2)=i+2
!!      j3=indices(i+3)
!!      reverse_indices(j3)=i+3
!!  end do
!!  !$omp end do
!!
!!  !$omp end parallel
!!
!!  !!do i=1,n
!!  !!    j=indices(i)
!!  !!    reverse_indices(j)=i
!!  !!end do
!!
!!end subroutine get_reverse_indices
