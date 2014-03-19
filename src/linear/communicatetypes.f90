!> @file 
!!   Routines to communicate types
!! @author
!!   Copyright (C) 2011-2012 BigDFT group 
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS


subroutine communicate_locreg_descriptors_basics(iproc, nlr, rootarr, orbs, llr)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nlr
  integer,dimension(nlr),intent(in) :: rootarr
  type(orbitals_data),intent(in) :: orbs
  type(locreg_descriptors),dimension(nlr),intent(inout) :: llr

  ! Local variables
  integer:: ierr, istat, iall, ilr, iilr
  character(len=1),dimension(:),allocatable :: worksend_char, workrecv_char
  logical,dimension(:),allocatable :: worksend_log, workrecv_log
  integer,dimension(:,:),allocatable :: worksend_int, workrecv_int
  real(8),dimension(:,:),allocatable :: worksend_dbl, workrecv_dbl
  character(len=*),parameter :: subname='communicate_locreg_descriptors_basics'

  allocate(worksend_char(orbs%norbp), stat=istat)
  call memocc(istat, worksend_char, 'worksend_char', subname)
  allocate(worksend_log(orbs%norbp), stat=istat)
  call memocc(istat, worksend_log, 'worksend_log', subname)
  allocate(worksend_int(11,orbs%norbp), stat=istat)
  call memocc(istat, worksend_int, 'worksend_int', subname)
  allocate(worksend_dbl(5,orbs%norbp), stat=istat)
  call memocc(istat, worksend_dbl, 'worksend_dbl', subname)

  allocate(workrecv_char(orbs%norb), stat=istat)
  call memocc(istat, workrecv_char, 'workrecv_char', subname)
  allocate(workrecv_log(orbs%norb), stat=istat)
  call memocc(istat, workrecv_log, 'workrecv_log', subname)
  allocate(workrecv_int(11,orbs%norb), stat=istat)
  call memocc(istat, workrecv_int, 'workrecv_int', subname)
  allocate(workrecv_dbl(5,orbs%norb), stat=istat)
  call memocc(istat, workrecv_dbl, 'workrecv_dbl', subname)


  iilr=0
  do ilr=1,nlr
      if (iproc==rootarr(ilr)) then
          iilr=iilr+1
          worksend_char(iilr)=llr(ilr)%geocode
          worksend_log(iilr)=llr(ilr)%hybrid_on
          worksend_int(1,iilr)=llr(ilr)%ns1
          worksend_int(2,iilr)=llr(ilr)%ns2
          worksend_int(3,iilr)=llr(ilr)%ns3
          worksend_int(4,iilr)=llr(ilr)%nsi1
          worksend_int(5,iilr)=llr(ilr)%nsi2
          worksend_int(6,iilr)=llr(ilr)%nsi3
          worksend_int(7,iilr)=llr(ilr)%localnorb
          worksend_int(8:10,iilr)=llr(ilr)%outofzone(1:3)
          worksend_int(11,iilr)=ilr
          worksend_dbl(1:3,iilr)=llr(ilr)%locregCenter(1:3)
          worksend_dbl(4,iilr)=llr(ilr)%locrad
          worksend_dbl(5,iilr)=llr(ilr)%locrad_kernel
      end if
  end do

  call mpi_allgatherv(worksend_char, orbs%norbp, mpi_character, workrecv_char, orbs%norb_par(:,0), &
       orbs%isorb_par, mpi_character, bigdft_mpi%mpi_comm, ierr)
  call mpi_allgatherv(worksend_log, orbs%norbp, mpi_logical, workrecv_log, orbs%norb_par(:,0), &
       orbs%isorb_par, mpi_logical, bigdft_mpi%mpi_comm, ierr)
  call mpi_allgatherv(worksend_int, 11*orbs%norbp, mpi_integer, workrecv_int, 11*orbs%norb_par(:,0), &
       11*orbs%isorb_par, mpi_integer, bigdft_mpi%mpi_comm, ierr)
  call mpi_allgatherv(worksend_dbl, 5*orbs%norbp, mpi_double_precision, workrecv_dbl, 5*orbs%norb_par(:,0), &
       5*orbs%isorb_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)

  do ilr=1,nlr
      iilr=workrecv_int(11,ilr)
      llr(iilr)%geocode=workrecv_char(ilr)
      llr(iilr)%hybrid_on= workrecv_log(ilr)
      llr(iilr)%ns1=workrecv_int(1,ilr)
      llr(iilr)%ns2=workrecv_int(2,ilr)
      llr(iilr)%ns3=workrecv_int(3,ilr)
      llr(iilr)%nsi1=workrecv_int(4,ilr)
      llr(iilr)%nsi2=workrecv_int(5,ilr)
      llr(iilr)%nsi3=workrecv_int(6,ilr)
      llr(iilr)%localnorb=workrecv_int(7,ilr)
      llr(iilr)%outofzone(1:3)=workrecv_int(8:10,ilr)
      llr(iilr)%locregCenter(1:3)=workrecv_dbl(1:3,ilr)
      llr(iilr)%locrad=workrecv_dbl(4,ilr)
      llr(iilr)%locrad_kernel=workrecv_dbl(5,ilr)
  end do


  iall=-product(shape(worksend_int))*kind(worksend_int)
  deallocate(worksend_int,stat=istat)
  call memocc(istat, iall, 'worksend_int', subname)
  iall=-product(shape(workrecv_int))*kind(workrecv_int)
  deallocate(workrecv_int,stat=istat)
  call memocc(istat, iall, 'workrecv_int', subname)
  allocate(worksend_int(13,orbs%norbp), stat=istat)
  call memocc(istat, worksend_int, 'worksend_int', subname)
  allocate(workrecv_int(13,orbs%norb), stat=istat)
  call memocc(istat, workrecv_int, 'workrecv_int', subname)


  iilr=0
  do ilr=1,nlr
      if (iproc==rootarr(ilr)) then
          iilr=iilr+1
          worksend_int(1,iilr)=llr(ilr)%d%n1
          worksend_int(2,iilr)=llr(ilr)%d%n2
          worksend_int(3,iilr)=llr(ilr)%d%n3
          worksend_int(4,iilr)=llr(ilr)%d%nfl1
          worksend_int(5,iilr)=llr(ilr)%d%nfu1
          worksend_int(6,iilr)=llr(ilr)%d%nfl2
          worksend_int(7,iilr)=llr(ilr)%d%nfu2
          worksend_int(8,iilr)=llr(ilr)%d%nfl3
          worksend_int(9,iilr)=llr(ilr)%d%nfu3
          worksend_int(10,iilr)=llr(ilr)%d%n1i
          worksend_int(11,iilr)=llr(ilr)%d%n2i
          worksend_int(12,iilr)=llr(ilr)%d%n3i
          worksend_int(13,iilr)=ilr
      end if
  end do

  call mpi_allgatherv(worksend_int, 13*orbs%norbp, mpi_integer, workrecv_int, 13*orbs%norb_par(:,0), &
       13*orbs%isorb_par, mpi_integer, bigdft_mpi%mpi_comm, ierr)

  do ilr=1,nlr
      iilr=workrecv_int(13,ilr)
      llr(iilr)%d%n1=workrecv_int(1,ilr)
      llr(iilr)%d%n2=workrecv_int(2,ilr)
      llr(iilr)%d%n3=workrecv_int(3,ilr)
      llr(iilr)%d%nfl1=workrecv_int(4,ilr)
      llr(iilr)%d%nfu1=workrecv_int(5,ilr)
      llr(iilr)%d%nfl2=workrecv_int(6,ilr)
      llr(iilr)%d%nfu2=workrecv_int(7,ilr)
      llr(iilr)%d%nfl3=workrecv_int(8,ilr)
      llr(iilr)%d%nfu3=workrecv_int(9,ilr)
      llr(iilr)%d%n1i=workrecv_int(10,ilr)
      llr(iilr)%d%n2i=workrecv_int(11,ilr)
      llr(iilr)%d%n3i=workrecv_int(12,ilr)
  end do


  iall=-product(shape(worksend_char))*kind(worksend_char)
  deallocate(worksend_char,stat=istat)
  call memocc(istat, iall, 'worksend_char', subname)
  iall=-product(shape(worksend_log))*kind(worksend_log)
  deallocate(worksend_log,stat=istat)
  call memocc(istat, iall, 'worksend_log', subname)
  iall=-product(shape(worksend_int))*kind(worksend_int)
  deallocate(worksend_int,stat=istat)
  call memocc(istat, iall, 'worksend_int', subname)
  iall=-product(shape(worksend_dbl))*kind(worksend_dbl)
  deallocate(worksend_dbl,stat=istat)
  call memocc(istat, iall, 'worksend_dbl', subname)

  iall=-product(shape(workrecv_char))*kind(workrecv_char)
  deallocate(workrecv_char,stat=istat)
  call memocc(istat, iall, 'workrecv_char', subname)
  iall=-product(shape(workrecv_log))*kind(workrecv_log)
  deallocate(workrecv_log,stat=istat)
  call memocc(istat, iall, 'workrecv_log', subname)
  iall=-product(shape(workrecv_int))*kind(workrecv_int)
  deallocate(workrecv_int,stat=istat)
  call memocc(istat, iall, 'workrecv_int', subname)
  iall=-product(shape(workrecv_dbl))*kind(workrecv_dbl)
  deallocate(workrecv_dbl,stat=istat)
  call memocc(istat, iall, 'workrecv_dbl', subname)

end subroutine communicate_locreg_descriptors_basics


subroutine communicate_locreg_descriptors_keys(iproc, nproc, nlr, glr, llr, orbs, rootarr, onwhichmpi)
   use module_base
   use module_types
   use yaml_output
   implicit none

   ! Calling arguments
   integer,intent(in):: iproc, nproc, nlr
   type(locreg_descriptors),intent(in) :: glr
   type(locreg_descriptors),dimension(nlr),intent(inout) :: llr
   type(orbitals_data),intent(in) :: orbs
   integer,dimension(nlr),intent(in) :: rootarr
   integer,dimension(orbs%norb),intent(in) :: onwhichmpi

   ! Local variables
   integer:: ierr, istat, iall, jorb, ilr, jlr, itask, jtask, root, icomm, nrecv, nalloc, max_sim_comms
   integer :: maxrecvdim, maxsenddim, nsend
   logical :: isoverlap
   character(len=*),parameter:: subname='communicate_wavefunctions_descriptors2'
   integer,dimension(:),allocatable :: requests
   integer,dimension(:,:),allocatable :: worksend_int, workrecv_int
   logical,dimension(:,:),allocatable :: covered
   !integer :: total_sent, total_recv

   ! This maxval is put out of the allocate to avoid compiler crash with PathScale.
   jorb = maxval(orbs%norb_par(:,0))
   allocate(requests(8*nproc*jorb), stat=istat)
   call memocc(istat, requests, 'requests', subname)

   allocate(covered(nlr,0:nproc-1), stat=istat)
   call memocc(istat, covered, 'covered', subname)

   allocate(worksend_int(4,nlr), stat=istat)
   call memocc(istat, worksend_int, 'worksend_int', subname)

   allocate(workrecv_int(4,nlr), stat=istat)
   call memocc(istat, workrecv_int, 'workrecv_int', subname)

   nrecv=0
   !nsend=0
   icomm=0
   maxsenddim=0
   do ilr=1,nlr
       root=rootarr(ilr)
       covered(ilr,:)=.false.
       do jorb=1,orbs%norb
           jlr=orbs%inwhichlocreg(jorb)
           jtask=onwhichmpi(jorb)
           ! check we're on a sending or receiving proc
           if (iproc /= root .and. iproc /= jtask) cycle
           ! don't communicate to ourselves, or if we've already sent this locreg
           if (jtask == root .or. covered(ilr,jtask)) cycle
           call check_overlap_cubic_periodic(glr,llr(ilr),llr(jlr),isoverlap)
           if (isoverlap) then         
               covered(ilr,jtask)=.true.
               if (iproc == root) then
                  !write(*,'(5(a,i0))') 'process ',iproc,' sends locreg ',ilr,' to process ',&
                  !    jtask,' with tags ',4*ilr+0,'-',4*ilr+3
                  worksend_int(1,ilr)=llr(ilr)%wfd%nvctr_c
                  worksend_int(2,ilr)=llr(ilr)%wfd%nvctr_f
                  worksend_int(3,ilr)=llr(ilr)%wfd%nseg_c
                  worksend_int(4,ilr)=llr(ilr)%wfd%nseg_f
                  icomm=icomm+1
                  call mpi_isend(worksend_int(1,ilr), 4, mpi_integer, jtask,&
                       itag(ilr,jtask), bigdft_mpi%mpi_comm, requests(icomm), ierr)
                  maxsenddim=max(maxsenddim,llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f)
                  !nsend=nsend+1
               else if (iproc == jtask) then
                  !write(*,'(5(a,i0))') 'process ',iproc,' receives locreg ',ilr,' from process ',&
                  !    root,' with tags ',4*ilr+0,'-',4*ilr+3
                  icomm=icomm+1
                  call mpi_irecv(workrecv_int(1,ilr), 4, mpi_integer, root,&
                       itag(ilr,jtask), bigdft_mpi%mpi_comm, requests(icomm), ierr)
                  nrecv=nrecv+1
               end if
           end if
       end do
   end do
  
   call mpi_waitall(icomm, requests(1), mpi_statuses_ignore, ierr)
   call mpi_barrier(mpi_comm_world,ierr)

   iall=-product(shape(worksend_int))*kind(worksend_int)
   deallocate(worksend_int,stat=istat)
   call memocc(istat, iall, 'worksend_int', subname)

   nalloc=0
   maxrecvdim=0
   do jlr=1,nlr 
      if (covered(jlr,iproc)) then
         llr(jlr)%wfd%nvctr_c=workrecv_int(1,jlr)
         llr(jlr)%wfd%nvctr_f=workrecv_int(2,jlr)
         llr(jlr)%wfd%nseg_c=workrecv_int(3,jlr)
         llr(jlr)%wfd%nseg_f=workrecv_int(4,jlr)
!         call allocate_wfd(llr(jlr)%wfd,subname)
         nalloc=nalloc+1
         maxrecvdim=max(maxrecvdim,llr(jlr)%wfd%nseg_c+llr(jlr)%wfd%nseg_f)
      end if
   end do
   if (f_err_raise(nalloc /= nrecv,'problem in communicate locregs: mismatch in receives '//&
        trim(yaml_toa(nrecv))//' and allocates '//trim(yaml_toa(nalloc))//' for process '//trim(yaml_toa(iproc)),&
        err_name='BIGDFT_RUNTIME_ERROR')) return

   iall=-product(shape(workrecv_int))*kind(workrecv_int)
   deallocate(workrecv_int,stat=istat)
   call memocc(istat, iall, 'workrecv_int', subname)

   !should reduce memory by not allocating for all llr
   allocate(workrecv_int(6*maxrecvdim,nlr), stat=istat)
   call memocc(istat, workrecv_int, 'workrecv_int', subname)
   allocate(worksend_int(6*maxsenddim,nlr), stat=istat)
   call memocc(istat, worksend_int, 'worksend_int', subname)

   ! divide communications into chunks to avoid problems with memory (too many communications)
   ! set maximum number of simultaneous communications
   !total_sent=0
   !total_recv=0
   max_sim_comms=10000
   icomm=0
   do ilr=1,nlr
      root=rootarr(ilr)
      do jtask=0,nproc-1
         if (.not. covered(ilr,jtask)) cycle
         if (iproc == root) then
           !write(*,'(5(a,i0))') 'process ',iproc,' sends locreg ',ilr,' to process ',&
           !     jtask,' with tags ',4*ilr+0,'-',4*ilr+3
           call vcopy(2*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f),llr(ilr)%wfd%keyglob(1,1),1,worksend_int(1,ilr),1)
           call vcopy(2*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f),llr(ilr)%wfd%keygloc(1,1),1,&
                worksend_int(2*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f)+1,ilr),1)
           call vcopy((llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f),llr(ilr)%wfd%keyvloc(1),1,&
                worksend_int(4*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f)+1,ilr),1)
           call vcopy((llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f),llr(ilr)%wfd%keyvglob(1),1,&
                worksend_int(5*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f)+1,ilr),1)
           icomm=icomm+1
           call mpi_isend(worksend_int(1,ilr),6*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f), mpi_integer, &
                jtask, itag(ilr,jtask), bigdft_mpi%mpi_comm, requests(icomm), ierr)
         else if (iproc == jtask) then
            !write(*,'(5(a,i0))') 'process ',iproc,' receives locreg ',ilr,' from process ',&
            !    root,' with tags ',4*ilr+0,'-',4*ilr+3
            icomm=icomm+1
            call mpi_irecv(workrecv_int(1,ilr),6*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f), mpi_integer, &
                 root, itag(ilr,jtask), bigdft_mpi%mpi_comm, requests(icomm), ierr)
         end if
      end do
      if (mod(ilr,max_sim_comms)==0 .or. ilr==nlr) then
         do jlr=max(ilr-max_sim_comms+1,1),ilr
            if (covered(jlr,iproc))  call allocate_wfd(llr(jlr)%wfd)
         end do
         call mpi_waitall(icomm, requests(1), mpi_statuses_ignore, ierr)
         if (f_err_raise(ierr /= 0,'problem in communicate locregs: error in mpi_waitall '//&
              trim(yaml_toa(ierr))//' for process '//trim(yaml_toa(iproc)),&
              err_name='BIGDFT_RUNTIME_ERROR')) return
         call mpi_barrier(mpi_comm_world,ierr)
         icomm=0
      end if
   end do

   iall=-product(shape(worksend_int))*kind(worksend_int)
   deallocate(worksend_int,stat=istat)
   call memocc(istat, iall, 'worksend_int', subname)

   do ilr=1,nlr 
      if (covered(ilr,iproc)) then
         call vcopy(2*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f),workrecv_int(1,ilr),1,llr(ilr)%wfd%keyglob(1,1),1)
         call vcopy(2*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f),&
              workrecv_int(2*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f)+1,ilr),1,llr(ilr)%wfd%keygloc(1,1),1)
         call vcopy((llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f),&
              workrecv_int(4*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f)+1,ilr),1,llr(ilr)%wfd%keyvloc(1),1)
         call vcopy((llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f),&
              workrecv_int(5*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f)+1,ilr),1,llr(ilr)%wfd%keyvglob(1),1)
      end if
   end do

   iall=-product(shape(workrecv_int))*kind(workrecv_int)
   deallocate(workrecv_int,stat=istat)
   call memocc(istat, iall, 'workrecv_int', subname)

   !print*,'iproc,sent,received,num sent,num received',iproc,total_sent,total_recv,nsend,nrecv
   iall=-product(shape(requests))*kind(requests)
   deallocate(requests,stat=istat)
   call memocc(istat, iall, 'requests', subname)

   iall=-product(shape(covered))*kind(covered)
   deallocate(covered,stat=istat)
   call memocc(istat, iall, 'covered', subname)

contains

 pure function itag(ilr,recv)
 implicit none
 integer, intent(in) :: ilr,recv
 integer :: itag

 itag=ilr+recv*nlr

 end function itag

END SUBROUTINE communicate_locreg_descriptors_keys
