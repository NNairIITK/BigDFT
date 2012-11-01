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
  integer:: ierr, ncount, mpi_tmptype, istat, iall, iorb, iiorb
  type(locreg_descriptors) :: lr
  integer,dimension(12):: blocklengths, types
  integer(kind=mpi_address_kind),dimension(12):: dspls
  integer(kind=mpi_address_kind):: addr_geocode, addr_hybrid_on, addr_ns1, addr_ns2, addr_ns3, addr_nsi1, addr_nsi2
  integer(kind=mpi_address_kind):: addr_nsi3, addr_localnorb, addr_outofzone, addr_locregCenter, addr_locrad, addr_lr
  character(len=1),dimension(:),allocatable :: worksend_char, workrecv_char
  logical,dimension(:),allocatable :: worksend_log, workrecv_log
  integer,dimension(:,:),allocatable :: worksend_int, workrecv_int
  real(8),dimension(:,:),allocatable :: worksend_dbl, workrecv_dbl
  character(len=*),parameter :: subname='communicate_locreg_descriptors_basics'

  allocate(worksend_char(orbs%norbp), stat=istat)
  call memocc(istat, worksend_char, 'worksend_char', subname)
  allocate(worksend_log(orbs%norbp), stat=istat)
  call memocc(istat, worksend_log, 'worksend_log', subname)
  allocate(worksend_int(10,orbs%norbp), stat=istat)
  call memocc(istat, worksend_int, 'worksend_int', subname)
  allocate(worksend_dbl(4,orbs%norbp), stat=istat)
  call memocc(istat, worksend_dbl, 'worksend_dbl', subname)

  allocate(workrecv_char(orbs%norb), stat=istat)
  call memocc(istat, workrecv_char, 'workrecv_char', subname)
  allocate(workrecv_log(orbs%norb), stat=istat)
  call memocc(istat, workrecv_log, 'workrecv_log', subname)
  allocate(workrecv_int(10,orbs%norb), stat=istat)
  call memocc(istat, workrecv_int, 'workrecv_int', subname)
  allocate(workrecv_dbl(4,orbs%norb), stat=istat)
  call memocc(istat, workrecv_dbl, 'workrecv_dbl', subname)


  iiorb=0
  do iorb=1,orbs%norb
      if (iproc==rootarr(iorb)) then
          iiorb=iiorb+1
          worksend_char(iiorb)=llr(iorb)%geocode
          worksend_log(iiorb)=llr(iorb)%hybrid_on
          worksend_int(1,iiorb)=llr(iorb)%ns1
          worksend_int(2,iiorb)=llr(iorb)%ns2
          worksend_int(3,iiorb)=llr(iorb)%ns3
          worksend_int(4,iiorb)=llr(iorb)%nsi1
          worksend_int(5,iiorb)=llr(iorb)%nsi2
          worksend_int(6,iiorb)=llr(iorb)%nsi3
          worksend_int(7,iiorb)=llr(iorb)%localnorb
          worksend_int(8:10,iiorb)=llr(iorb)%outofzone(1:3)
          worksend_dbl(1:3,iiorb)=llr(iorb)%locregCenter(1:3)
          worksend_dbl(4,iiorb)=llr(iorb)%locrad
      end if
  end do

  call mpi_allgatherv(worksend_char, orbs%norbp, mpi_character, workrecv_char, orbs%norb_par(:,0), &
       orbs%isorb_par, mpi_character, bigdft_mpi%mpi_comm, ierr)
  call mpi_allgatherv(worksend_log, orbs%norbp, mpi_logical, workrecv_log, orbs%norb_par(:,0), &
       orbs%isorb_par, mpi_logical, bigdft_mpi%mpi_comm, ierr)
  call mpi_allgatherv(worksend_int, 10*orbs%norbp, mpi_integer, workrecv_int, 10*orbs%norb_par(:,0), &
       10*orbs%isorb_par, mpi_integer, bigdft_mpi%mpi_comm, ierr)
  call mpi_allgatherv(worksend_dbl, 4*orbs%norbp, mpi_double_precision, workrecv_dbl, 4*orbs%norb_par(:,0), &
       4*orbs%isorb_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)

  do iorb=1,orbs%norb
      iiorb=iiorb+1
      llr(iorb)%geocode=workrecv_char(iorb)
      llr(iorb)%hybrid_on= workrecv_log(iorb)
      llr(iorb)%ns1=workrecv_int(1,iorb)
      llr(iorb)%ns2=workrecv_int(2,iorb)
      llr(iorb)%ns3=workrecv_int(3,iorb)
      llr(iorb)%nsi1=workrecv_int(4,iorb)
      llr(iorb)%nsi2=workrecv_int(5,iorb)
      llr(iorb)%nsi3=workrecv_int(6,iorb)
      llr(iorb)%localnorb=workrecv_int(7,iorb)
      llr(iorb)%outofzone(1:3)=workrecv_int(8:10,iorb)
      llr(iorb)%locregCenter(1:3)=workrecv_dbl(1:3,iorb)
      llr(iorb)%locrad=workrecv_dbl(4,iorb)
  end do


  iall=-product(shape(worksend_int))*kind(worksend_int)
  deallocate(worksend_int,stat=istat)
  call memocc(istat, iall, 'worksend_int', subname)
  iall=-product(shape(workrecv_int))*kind(workrecv_int)
  deallocate(workrecv_int,stat=istat)
  call memocc(istat, iall, 'workrecv_int', subname)
  allocate(worksend_int(12,orbs%norbp), stat=istat)
  call memocc(istat, worksend_int, 'worksend_int', subname)
  allocate(workrecv_int(12,orbs%norb), stat=istat)
  call memocc(istat, workrecv_int, 'workrecv_int', subname)


  iiorb=0
  do iorb=1,orbs%norb
      if (iproc==rootarr(iorb)) then
          iiorb=iiorb+1
          worksend_int(1,iiorb)=llr(iorb)%d%n1
          worksend_int(2,iiorb)=llr(iorb)%d%n2
          worksend_int(3,iiorb)=llr(iorb)%d%n3
          worksend_int(4,iiorb)=llr(iorb)%d%nfl1
          worksend_int(5,iiorb)=llr(iorb)%d%nfu1
          worksend_int(6,iiorb)=llr(iorb)%d%nfl2
          worksend_int(7,iiorb)=llr(iorb)%d%nfu2
          worksend_int(8,iiorb)=llr(iorb)%d%nfl3
          worksend_int(9,iiorb)=llr(iorb)%d%nfu3
          worksend_int(10,iiorb)=llr(iorb)%d%n1i
          worksend_int(11,iiorb)=llr(iorb)%d%n2i
          worksend_int(12,iiorb)=llr(iorb)%d%n3i
      end if
  end do

  call mpi_allgatherv(worksend_int, 12*orbs%norbp, mpi_integer, workrecv_int, 12*orbs%norb_par(:,0), &
       12*orbs%isorb_par, mpi_integer, bigdft_mpi%mpi_comm, ierr)

  do iorb=1,orbs%norb
      llr(iorb)%d%n1=workrecv_int(1,iorb)
      llr(iorb)%d%n2=workrecv_int(2,iorb)
      llr(iorb)%d%n3=workrecv_int(3,iorb)
      llr(iorb)%d%nfl1=workrecv_int(4,iorb)
      llr(iorb)%d%nfu1=workrecv_int(5,iorb)
      llr(iorb)%d%nfl2=workrecv_int(6,iorb)
      llr(iorb)%d%nfu2=workrecv_int(7,iorb)
      llr(iorb)%d%nfl3=workrecv_int(8,iorb)
      llr(iorb)%d%nfu3=workrecv_int(9,iorb)
      llr(iorb)%d%n1i=workrecv_int(10,iorb)
      llr(iorb)%d%n2i=workrecv_int(11,iorb)
      llr(iorb)%d%n3i=workrecv_int(12,iorb)
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






subroutine communicate_locreg_descriptors_keys(iproc, nproc, nlr, glr, llr, orbs, orbsder, rootarr)
   use module_base
   use module_types
   implicit none

   ! Calling arguments
   integer,intent(in):: iproc, nproc, nlr
   type(locreg_descriptors),intent(in) :: glr
   type(locreg_descriptors),dimension(nlr),intent(inout) :: llr
   type(orbitals_data),intent(in) :: orbs, orbsder
   integer,dimension(orbs%norb),intent(in) :: rootarr

   ! Local variables
   integer:: ierr, ncount, commtype, istat, iall, iorb, jorb, ilr, jlr, itask, jtask, root, isend, irecv, jtaskder
   logical :: isoverlap
   character(len=*),parameter:: subname='communicate_wavefunctions_descriptors2'
   integer,dimension(4):: blocklengths,types
   integer(kind=mpi_address_kind):: addr_wfd, addr_nvctr_c, addr_nvctr_f, addr_nseg_c, addr_nseg_f
   integer(kind=mpi_address_kind),dimension(4):: dspls
   integer ,dimension(4):: itags
   integer, dimension(:), allocatable :: wrkarr
   integer,dimension(:,:),allocatable :: requests
   logical,dimension(:),allocatable :: covered

   allocate(requests(4*orbs%norb*orbs%norb,2), stat=istat)
   call memocc(istat, requests, 'requests', subname)

   allocate(covered(0:nproc-1), stat=istat)
   call memocc(istat, covered, 'covered', subname)


   isend=0
   irecv=0
   do iorb=1,orbs%norb
       ilr=orbs%inwhichlocreg(iorb)
       itask=orbs%onwhichmpi(iorb)
       root=rootarr(ilr)
       covered=.false.
       do jorb=1,orbs%norb
           jlr=orbs%inwhichlocreg(jorb)
           jtask=orbs%onwhichmpi(jorb)
           if (covered(jtask)) cycle
           !unambiguous mpi tags
           itags(1)=jtask+nproc*itask+(nproc**2)
           itags(2)=jtask+nproc*itask+(nproc**2)+1
           itags(3)=jtask+nproc*itask+(nproc**2)+2
           itags(4)=jtask+nproc*itask+(nproc**2)+3
           call check_overlap_cubic_periodic(glr,llr(ilr),llr(jlr),isoverlap)
           if (isoverlap) then
               covered(jtask)=.true.
               if (jtask /= root) then
                  if (iproc==root) then
                     isend=isend+1
                     call mpi_isend(llr(ilr)%wfd%nvctr_c, 1, mpi_integer, jtask,&
                          4*jtask+0, bigdft_mpi%mpi_comm, requests(isend,1), ierr)
                     isend=isend+1
                     call mpi_isend(llr(ilr)%wfd%nvctr_f, 1, mpi_integer, jtask,&
                          4*jtask+1, bigdft_mpi%mpi_comm, requests(isend,1), ierr)
                     isend=isend+1
                     call mpi_isend(llr(ilr)%wfd%nseg_c, 1, mpi_integer, jtask, &
                          4*jtask+2, bigdft_mpi%mpi_comm, requests(isend,1), ierr)
                     isend=isend+1
                     call mpi_isend(llr(ilr)%wfd%nseg_f, 1, mpi_integer, jtask, &
                          4*jtask+3, bigdft_mpi%mpi_comm, requests(isend,1), ierr)
                  else if (iproc==jtask) then
                     irecv=irecv+1
                     call mpi_irecv(llr(ilr)%wfd%nvctr_c, 1, mpi_integer, root,&
                          4*iproc+0, bigdft_mpi%mpi_comm, requests(irecv,2), ierr)
                     irecv=irecv+1
                     call mpi_irecv(llr(ilr)%wfd%nvctr_f, 1, mpi_integer, root,&
                          4*iproc+1, bigdft_mpi%mpi_comm, requests(irecv,2), ierr)
                     irecv=irecv+1
                     call mpi_irecv(llr(ilr)%wfd%nseg_c, 1, mpi_integer, root,&
                          4*iproc+2, bigdft_mpi%mpi_comm, requests(irecv,2), ierr)
                     irecv=irecv+1
                     call mpi_irecv(llr(ilr)%wfd%nseg_f, 1, mpi_integer, root,&
                          4*iproc+3, bigdft_mpi%mpi_comm, requests(irecv,2), ierr)
                  end if
               end if
           end if
       end do

       do jorb=1,orbsder%norb
           jlr=orbsder%inwhichlocreg(jorb)
           jtaskder=orbsder%onwhichmpi(jorb)
           if (covered(jtaskder)) cycle
           !unambiguous mpi tags
           itags(1)=jtaskder+nproc*itask+(nproc**2)
           itags(2)=jtaskder+nproc*itask+(nproc**2)+1
           itags(3)=jtaskder+nproc*itask+(nproc**2)+2
           itags(4)=jtaskder+nproc*itask+(nproc**2)+3

           call check_overlap_cubic_periodic(glr,llr(ilr),llr(jlr),isoverlap)
           if (isoverlap) then
               covered(jtaskder)=.true.
               if (jtaskder /= root) then
                  if (iproc==root) then
                     isend=isend+1
                     call mpi_isend(llr(ilr)%wfd%nvctr_c, 1, mpi_integer, jtaskder, 4*jtaskder+0, &
                          bigdft_mpi%mpi_comm, requests(isend,1), ierr)
                     isend=isend+1
                     call mpi_isend(llr(ilr)%wfd%nvctr_f, 1, mpi_integer, jtaskder, 4*jtaskder+1, &
                          bigdft_mpi%mpi_comm, requests(isend,1), ierr)
                     isend=isend+1
                     call mpi_isend(llr(ilr)%wfd%nseg_c, 1, mpi_integer, jtaskder, 4*jtaskder+2, &
                          bigdft_mpi%mpi_comm, requests(isend,1), ierr)
                     isend=isend+1
                     call mpi_isend(llr(ilr)%wfd%nseg_f, 1, mpi_integer, jtaskder, 4*jtaskder+3, &
                          bigdft_mpi%mpi_comm, requests(isend,1), ierr)
                  else if (iproc==jtaskder) then
                     irecv=irecv+1
                     call mpi_irecv(llr(ilr)%wfd%nvctr_c, 1, mpi_integer, root, 4*iproc+0, &
                          bigdft_mpi%mpi_comm, requests(irecv,2), ierr)
                     irecv=irecv+1
                     call mpi_irecv(llr(ilr)%wfd%nvctr_f, 1, mpi_integer, root, 4*iproc+1, &
                          bigdft_mpi%mpi_comm, requests(irecv,2), ierr)
                     irecv=irecv+1
                     call mpi_irecv(llr(ilr)%wfd%nseg_c, 1, mpi_integer, root, 4*iproc+2, &
                          bigdft_mpi%mpi_comm, requests(irecv,2), ierr)
                     irecv=irecv+1
                     call mpi_irecv(llr(ilr)%wfd%nseg_f, 1, mpi_integer, root, 4*iproc+3, &
                          bigdft_mpi%mpi_comm, requests(irecv,2), ierr)
                  end if
               end if
           end if
       end do
   end do

   call mpi_waitall(irecv, requests(1,2), mpi_statuses_ignore, ierr)
   call mpi_waitall(isend, requests(1,1), mpi_statuses_ignore, ierr)


   do iorb=1,orbs%norb
       ilr=orbs%inwhichlocreg(iorb)
       itask=orbs%onwhichmpi(iorb)
       root=rootarr(ilr)
       covered=.false.
       do jorb=1,orbs%norb
           jlr=orbs%inwhichlocreg(jorb)
           jtask=orbs%onwhichmpi(jorb)
           if (covered(jtask)) cycle
           call check_overlap_cubic_periodic(glr,llr(ilr),llr(jlr),isoverlap)
           if (isoverlap) then
               covered(jtask)=.true.
               if (iproc==root .and. iproc/=jtask) then
               else if (iproc==jtask .and. iproc/=root) then
                   call allocate_wfd(llr(ilr)%wfd,subname)
               end if
           end if
       end do
       do jorb=1,orbsder%norb
           jlr=orbsder%inwhichlocreg(jorb)
           jtaskder=orbsder%onwhichmpi(jorb)
           if (covered(jtaskder)) cycle
           call check_overlap_cubic_periodic(glr,llr(ilr),llr(jlr),isoverlap)
           if (isoverlap) then
               covered(jtaskder)=.true.
               if (iproc==root .and. iproc/=jtaskder) then
               else if (iproc==jtaskder .and. iproc/=root) then
                   call allocate_wfd(llr(ilr)%wfd,subname)
               end if
           end if
       end do
   end do



   isend=0
   irecv=0
   do iorb=1,orbs%norb
       ilr=orbs%inwhichlocreg(iorb)
       itask=orbs%onwhichmpi(iorb)
       root=rootarr(ilr)
       covered=.false.
       do jorb=1,orbs%norb
           jlr=orbs%inwhichlocreg(jorb)
           jtask=orbs%onwhichmpi(jorb)
           if (covered(jtask)) cycle
           !unambiguous mpi tags
           itags(1)=jtask+nproc*itask+(nproc**2)
           itags(2)=jtask+nproc*itask+(nproc**2)+1
           itags(3)=jtask+nproc*itask+(nproc**2)+2
           itags(4)=jtask+nproc*itask+(nproc**2)+3

           call check_overlap_cubic_periodic(glr,llr(ilr),llr(jlr),isoverlap)
           if (isoverlap) then
              covered(jtask)=.true.
              if (jtask /= root) then
                 if (iproc==root) then
                    isend=isend+1
                    call mpi_isend(llr(ilr)%wfd%keyglob, 2*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f), mpi_integer, &
                         jtask, 4*jtask+0, bigdft_mpi%mpi_comm, requests(isend,1), ierr)
                    isend=isend+1
                    call mpi_isend(llr(ilr)%wfd%keygloc, 2*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f), mpi_integer, &
                         jtask, 4*jtask+1, bigdft_mpi%mpi_comm, requests(isend,1), ierr)
                    isend=isend+1
                    call mpi_isend(llr(ilr)%wfd%keyvloc, llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f, mpi_integer, &
                         jtask, 4*jtask+2, bigdft_mpi%mpi_comm, requests(isend,1), ierr)
                    isend=isend+1
                    call mpi_isend(llr(ilr)%wfd%keyvglob, llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f, mpi_integer, &
                         jtask, 4*jtask+3, bigdft_mpi%mpi_comm, requests(isend,1), ierr)
                 else if (iproc==jtask) then
                    irecv=irecv+1
                    call mpi_irecv(llr(ilr)%wfd%keyglob, 2*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f), mpi_integer, &
                         root, 4*iproc+0, bigdft_mpi%mpi_comm, requests(irecv,2), ierr)
                    irecv=irecv+1
                    call mpi_irecv(llr(ilr)%wfd%keygloc, 2*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f), mpi_integer, &
                         root, 4*iproc+1, bigdft_mpi%mpi_comm, requests(irecv,2), ierr)
                    irecv=irecv+1
                    call mpi_irecv(llr(ilr)%wfd%keyvloc, llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f, mpi_integer, &
                         root, 4*iproc+2, bigdft_mpi%mpi_comm, requests(irecv,2), ierr)
                    irecv=irecv+1
                    call mpi_irecv(llr(ilr)%wfd%keyvglob, llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f, mpi_integer, &
                         root, 4*iproc+3, bigdft_mpi%mpi_comm, requests(irecv,2), ierr)
                 end if
              end if
           end if
        end do
        do jorb=1,orbsder%norb
           jlr=orbsder%inwhichlocreg(jorb)
           jtaskder=orbsder%onwhichmpi(jorb)
           if (covered(jtaskder)) cycle
           !unambiguous mpi tags
           itags(1)=jtaskder+nproc*itask+(nproc**2)
           itags(2)=jtaskder+nproc*itask+(nproc**2)+1
           itags(3)=jtaskder+nproc*itask+(nproc**2)+2
           itags(4)=jtaskder+nproc*itask+(nproc**2)+3

           call check_overlap_cubic_periodic(glr,llr(ilr),llr(jlr),isoverlap)
           if (isoverlap) then
               covered(jtaskder)=.true.
               if (jtaskder /= root) then
                  if (iproc==root) then
                     isend=isend+1
                     call mpi_isend(llr(ilr)%wfd%keyglob, 2*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f), mpi_integer, &
                          jtaskder, 4*jtaskder+0, bigdft_mpi%mpi_comm, requests(isend,1), ierr)
                     isend=isend+1
                     call mpi_isend(llr(ilr)%wfd%keygloc, 2*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f), mpi_integer, &
                          jtaskder, 4*jtaskder+1, bigdft_mpi%mpi_comm, requests(isend,1), ierr)
                     isend=isend+1
                     call mpi_isend(llr(ilr)%wfd%keyvloc, llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f, mpi_integer, &
                          jtaskder, 4*jtaskder+2, bigdft_mpi%mpi_comm, requests(isend,1), ierr)
                     isend=isend+1
                     call mpi_isend(llr(ilr)%wfd%keyvglob, llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f, mpi_integer, &
                          jtaskder, 4*jtaskder+3, bigdft_mpi%mpi_comm, requests(isend,1), ierr)
                  else if (iproc==jtaskder) then
                     irecv=irecv+1
                     call mpi_irecv(llr(ilr)%wfd%keyglob, 2*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f), mpi_integer, &
                          root, 4*iproc+0, bigdft_mpi%mpi_comm, requests(irecv,2), ierr)
                     irecv=irecv+1
                     call mpi_irecv(llr(ilr)%wfd%keygloc, 2*(llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f), mpi_integer, &
                          root, 4*iproc+1, bigdft_mpi%mpi_comm, requests(irecv,2), ierr)
                     irecv=irecv+1
                     call mpi_irecv(llr(ilr)%wfd%keyvloc, llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f, mpi_integer, &
                          root, 4*iproc+2, bigdft_mpi%mpi_comm, requests(irecv,2), ierr)
                     irecv=irecv+1
                     call mpi_irecv(llr(ilr)%wfd%keyvglob, llr(ilr)%wfd%nseg_c+llr(ilr)%wfd%nseg_f, mpi_integer, &
                          root, 4*iproc+3, bigdft_mpi%mpi_comm, requests(irecv,2), ierr)
                  end if
               end if
           end if
       end do
   end do

   call mpi_waitall(isend, requests(1,1), mpi_statuses_ignore, ierr)
   call mpi_waitall(irecv, requests(1,2), mpi_statuses_ignore, ierr)

   iall=-product(shape(requests))*kind(requests)
   deallocate(requests,stat=istat)
   call memocc(istat, iall, 'requests', subname)

   iall=-product(shape(covered))*kind(covered)
   deallocate(covered,stat=istat)
   call memocc(istat, iall, 'covered', subname)



END SUBROUTINE communicate_locreg_descriptors_keys
