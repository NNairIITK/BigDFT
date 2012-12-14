subroutine post_p2p_communication(iproc, nproc, nsendbuf, sendbuf, nrecvbuf, recvbuf, comm, lzd)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc, nsendbuf, nrecvbuf
  real(8),dimension(nsendbuf),intent(in):: sendbuf
  real(8),dimension(nrecvbuf),intent(out):: recvbuf
  type(p2pComms),intent(inout):: comm
  type(local_zone_descriptors),intent(in) :: lzd
  
  ! Local variables
  integer:: jproc, joverlap, nsends, nreceives, mpisource, istsource, ncount, mpidest, istdest, tag, ierr, it, nit
  integer :: ioffset_send, ioffset_recv, maxit, mpi_type, istat, iall
  integer :: ist, i2, i3, ist2, ist3, istrecv, istsend, info
  integer :: ncnt, nblocklength, nstride
  integer :: nsize, size_of_double
  integer,dimension(:),allocatable :: blocklengths, types
  integer(kind=mpi_address_kind),dimension(:),allocatable :: displacements
  character(len=*),parameter :: subname='post_p2p_communication'



  if(.not.comm%communication_complete) stop 'ERROR: there is already a p2p communication going on...'

! OLD VERSION
! #################################################################     
!!!  nproc_if: if (nproc>1) then
!!!
!!!      maxit = maxval(comm%comarr(7,:,:))
!!!
!!!      allocate(blocklengths(maxit), stat=istat)
!!!      call memocc(istat, blocklengths, 'blocklengths', subname)
!!!      allocate(types(maxit), stat=istat)
!!!      call memocc(istat, types, 'types', subname)
!!!      allocate(displacements(maxit), stat=istat)
!!!      !call memocc(istat, displacements, 'displacements', subname)
!!!
!!!      !!write(*,*) 'iproc, maxit', iproc, maxit
!!!
!!!      call mpi_type_size(mpi_double_precision, size_of_double, ierr)
!!!      
!!!      nreceives=0
!!!      nsends=0
!!!      do jproc=0,nproc-1
!!!          do joverlap=1,comm%noverlaps(jproc)
!!!              mpisource=comm%comarr(1,joverlap,jproc)
!!!              istsource=comm%comarr(2,joverlap,jproc)
!!!              ncount=comm%comarr(3,joverlap,jproc)
!!!              mpidest=comm%comarr(4,joverlap,jproc)
!!!              istdest=comm%comarr(5,joverlap,jproc)
!!!              tag=comm%comarr(6,joverlap,jproc)
!!!              nit=comm%comarr(7,joverlap,jproc)
!!!              ioffset_send=comm%comarr(8,joverlap,jproc)
!!!              ioffset_recv=comm%comarr(9,joverlap,jproc)
!!!              !!do it=1,nit
!!!              !!    blocklengths(it)=1
!!!              !!    displacements(it)=int(size_of_double*(it-1)*ioffset_send,kind=mpi_address_kind)
!!!              !!    types(it)=comm%mpi_datatypes(1,jproc)
!!!              !!end do
!!!              !!call mpi_type_create_struct(nit, blocklengths, displacements, &
!!!              !!     types,  mpi_type, ierr)
!!!              if (iproc==mpidest .or. iproc==mpisource) then
!!!                  call mpi_type_create_hvector(nit, 1, int(size_of_double*ioffset_send,kind=mpi_address_kind), &
!!!                       comm%mpi_datatypes(1,jproc), mpi_type, ierr)
!!!                  call mpi_type_commit(mpi_type, ierr)
!!!                  call mpi_type_size(mpi_type, nsize, ierr)
!!!                  nsize=nsize/size_of_double
!!!                  !if (iproc==0) write(*,'(a,4i12)') 'jproc, joverlap, mpi_type, nsize', jproc, joverlap, mpi_type, nsize
!!!                  if(nsize>0) then
!!!                      if(iproc==mpidest) then
!!!                          tag=mpidest
!!!                          nreceives=nreceives+1
!!!                          !!write(1200+iproc,'(5(a,i0))') 'process ',iproc,' receives ', nsize,&
!!!                          !!    ' elements from process ',mpisource,' with tag ',tag,' at position ',&
!!!                          !!    istdest
!!!                          call mpi_irecv(recvbuf(istdest), nsize, mpi_double_precision, mpisource, &
!!!                               tag, bigdft_mpi%mpi_comm, comm%requests(nreceives,2), ierr)
!!!                      end if
!!!                      if (iproc==mpisource) then
!!!                          !tag=mpidest*maxit
!!!                          tag=mpidest
!!!                          nsends=nsends+1
!!!                          !!write(1200+iproc,'(5(a,i0))') 'process ',mpisource,' sends ',ncount*nsize,&
!!!                          !!    ' elements from position ',istsource,' to process ',&
!!!                          !!    mpidest,' with tag ',tag
!!!                          call mpi_isend(sendbuf(istsource), ncount, mpi_type, mpidest, &
!!!                               tag, bigdft_mpi%mpi_comm, comm%requests(nsends,1), ierr)
!!!                      end if
!!!                  end if
!!!                  call mpi_type_free(mpi_type, ierr)
!!!              end if
!!!          end do
!!!      end do
!!!
!!!      iall=-product(shape(blocklengths))*kind(blocklengths)
!!!      deallocate(blocklengths,stat=istat)
!!!      call memocc(istat,iall,'blocklengths',subname)
!!!
!!!      iall=-product(shape(types))*kind(types)
!!!      deallocate(types,stat=istat)
!!!      call memocc(istat,iall,'types',subname)
!!!
!!!      !iall=-product(shape(displacements))*kind(displacements)
!!!      deallocate(displacements,stat=istat)
!!!      !call memocc(istat,iall,'displacements',subname)
!!!
!!!  else nproc_if
!!!
!!!      nreceives=0
!!!      nsends=0
!!!
!!!      ist=1
!!!      do i3=comm%ise(5,iproc),comm%ise(6,iproc)
!!!          ist3=(i3-1)*lzd%glr%d%n1i*lzd%glr%d%n2i
!!!          do i2=comm%ise(3,iproc),comm%ise(4,iproc)
!!!              ist2=(i2-1)*lzd%glr%d%n1i
!!!              !write(*,*) 'ist3+ist2+1, ist', ist3+ist2+1, ist
!!!              call dcopy(comm%ise(2,iproc)-comm%ise(1,iproc)+1, sendbuf(ist3+ist2+1), 1, recvbuf(ist), 1)
!!!              ist=ist+comm%ise(2,iproc)-comm%ise(1,iproc)+1
!!!          end do
!!!      end do
!!!
!!!
!!!  end if nproc_if
!!!  
!!!  
!!!  comm%nsend=nsends
!!!  comm%nrecv=nreceives
!!!  
!!!  ! Flag indicating whether the communication is complete or not
!!!  if(nproc>1) then
!!!      comm%communication_complete=.false.
!!!      comm%messages_posted=.true.
!!!  else
!!!      comm%communication_complete=.true.
!!!      comm%messages_posted=.false.
!!!  end if




! NEW VERSION
! #################################################################     
  nproc_if: if (nproc>1) then

      ! Allocate MPI memory window
      call mpi_type_size(mpi_double_precision, size_of_double, ierr)
      call mpi_info_create(info, ierr)
      call mpi_info_set(info, "no_locks", "true", ierr)
      !write(*,*) 'create window: iproc', iproc
      call mpi_win_create(sendbuf(1), int(nsendbuf*size_of_double,kind=mpi_address_kind), size_of_double, &
           info, bigdft_mpi%mpi_comm, comm%window, ierr)
      call mpi_info_free(info, ierr)

      call mpi_win_fence(mpi_mode_noprecede, comm%window, ierr)




      !!maxit = maxval(comm%comarr(7,:,:))

      !!allocate(blocklengths(maxit), stat=istat)
      !!call memocc(istat, blocklengths, 'blocklengths', subname)
      !!allocate(types(maxit), stat=istat)
      !!call memocc(istat, types, 'types', subname)
      !!allocate(displacements(maxit), stat=istat)
      !!!call memocc(istat, displacements, 'displacements', subname)

      !!write(*,*) 'iproc, maxit', iproc, maxit

      !!call mpi_type_size(mpi_double_precision, size_of_double, ierr)
      
      nreceives=0
      nsends=0
      do jproc=0,nproc-1
          do joverlap=1,comm%noverlaps(jproc)
              mpisource=comm%comarr(1,joverlap,jproc)
              istsource=comm%comarr(2,joverlap,jproc)
              ncount=comm%comarr(3,joverlap,jproc)
              mpidest=comm%comarr(4,joverlap,jproc)
              istdest=comm%comarr(5,joverlap,jproc)
              tag=comm%comarr(6,joverlap,jproc)
              nit=comm%comarr(7,joverlap,jproc)
              ioffset_send=comm%comarr(8,joverlap,jproc)
              ioffset_recv=comm%comarr(9,joverlap,jproc)
              !!do it=1,nit
              !!    blocklengths(it)=1
              !!    displacements(it)=int(size_of_double*(it-1)*ioffset_send,kind=mpi_address_kind)
              !!    types(it)=comm%mpi_datatypes(1,jproc)
              !!end do
              !!call mpi_type_create_struct(nit, blocklengths, displacements, &
              !!     types,  mpi_type, ierr)
              if (iproc==mpidest .or. iproc==mpisource) then
                  !!call mpi_type_create_hvector(nit, 1, int(size_of_double*ioffset_send,kind=mpi_address_kind), &
                  !!     comm%mpi_datatypes(1,jproc), mpi_type, ierr)
                  !!call mpi_type_commit(mpi_type, ierr)
                  !!call mpi_type_size(mpi_type, nsize, ierr)
                  !!nsize=nsize/size_of_double
                  !!!if (iproc==0) write(*,'(a,4i12)') 'jproc, joverlap, mpi_type, nsize', jproc, joverlap, mpi_type, nsize
                  !!if(nsize>0) then
                      call mpi_type_create_hvector(nit, 1, int(size_of_double*ioffset_send,kind=mpi_address_kind), &
                           comm%mpi_datatypes(1,jproc), mpi_type, ierr)
                      call mpi_type_commit(mpi_type, ierr)
                      call mpi_type_size(mpi_type, nsize, ierr)
                      nsize=nsize/size_of_double
                      if(iproc==mpidest .and. nsize>0) then
                          tag=mpidest
                          nreceives=nreceives+1
                          !!write(1200+iproc,'(5(a,i0))') 'process ',iproc,' receives ', nsize,&
                          !!    ' elements from process ',mpisource,' with tag ',tag,' at position ',&
                          !!    istdest
                          !!call mpi_irecv(recvbuf(istdest), nsize, mpi_double_precision, mpisource, &
                          !!     tag, bigdft_mpi%mpi_comm, comm%requests(nreceives,2), ierr)

                          !!istrecv=istdest
                          !!do i3=1,nit
                          !!    istsend=istsource+(i3-1)*lzd%glr%d%n1i*lzd%glr%d%n2i
                          !!    do i2=comm%ise(3,iproc),comm%ise(4,iproc)
                          !!        !!write(1000+iproc,'(a,4i9)') 'i3, i2, istrecv, nrecvbuf', i3, i2, istrecv, nrecvbuf
                          !!        call mpi_get(recvbuf(istrecv), comm%ise(2,iproc)-comm%ise(1,iproc)+1, &
                          !!             mpi_double_precision, mpisource, int((istsend-1),kind=mpi_address_kind), &
                          !!             comm%ise(2,iproc)-comm%ise(1,iproc)+1, mpi_double_precision, &
                          !!             comm%window, ierr)
                          !!        istsend=istsend+lzd%glr%d%n1i
                          !!        istrecv=istrecv+comm%ise(2,iproc)-comm%ise(1,iproc)+1
                          !!    end do
                          !!end do
                          call mpi_get(recvbuf(istdest), nsize, &
                               mpi_double_precision, mpisource, int((istsource-1),kind=mpi_address_kind), &
                               1, mpi_type, comm%window, ierr)

                      end if
                      call mpi_type_free(mpi_type, ierr)
                      !!if (iproc==mpisource) then
                      !!    !tag=mpidest*maxit
                      !!    tag=mpidest
                      !!    nsends=nsends+1
                      !!    !!write(1200+iproc,'(5(a,i0))') 'process ',mpisource,' sends ',ncount*nsize,&
                      !!    !!    ' elements from position ',istsource,' to process ',&
                      !!    !!    mpidest,' with tag ',tag
                      !!    call mpi_isend(sendbuf(istsource), ncount, mpi_type, mpidest, &
                      !!         tag, bigdft_mpi%mpi_comm, comm%requests(nsends,1), ierr)
                      !!end if
                  !!end if
                  !!call mpi_type_free(mpi_type, ierr)
              end if
          end do
      end do

      !!call mpi_win_fence(0, comm%window, ierr)
      !!call mpi_win_free(comm%window, ierr)
      !!do ist=1,nrecvbuf
      !!    write(200+iproc,*) ist, recvbuf(ist)
      !!end do

      !!iall=-product(shape(blocklengths))*kind(blocklengths)
      !!deallocate(blocklengths,stat=istat)
      !!call memocc(istat,iall,'blocklengths',subname)

      !!iall=-product(shape(types))*kind(types)
      !!deallocate(types,stat=istat)
      !!call memocc(istat,iall,'types',subname)

      !!!iall=-product(shape(displacements))*kind(displacements)
      !!deallocate(displacements,stat=istat)
      !!!call memocc(istat,iall,'displacements',subname)

  else nproc_if

      nreceives=0
      nsends=0

      ist=1
      do i3=comm%ise(5,iproc),comm%ise(6,iproc)
          ist3=(i3-1)*lzd%glr%d%n1i*lzd%glr%d%n2i
          do i2=comm%ise(3,iproc),comm%ise(4,iproc)
              ist2=(i2-1)*lzd%glr%d%n1i
              !write(*,*) 'ist3+ist2+1, ist', ist3+ist2+1, ist
              call dcopy(comm%ise(2,iproc)-comm%ise(1,iproc)+1, sendbuf(ist3+ist2+1), 1, recvbuf(ist), 1)
              ist=ist+comm%ise(2,iproc)-comm%ise(1,iproc)+1
          end do
      end do


  end if nproc_if
  
  
  comm%nsend=nsends
  comm%nrecv=nreceives
  









  !!if(nreceives/=comm%noverlaps(iproc)) then
  !!    write(*,'(1x,a,i0,a,i0,2x,i0)') 'ERROR on process ', iproc, ': nreceives/=comm%noverlaps(iproc)',&
  !!         nreceives, comm%noverlaps(iproc)
  !!  stop
  !!end if
  
  !!! Flag indicating whether the communication is complete or not
  if(nproc>1) then
      comm%communication_complete=.false.
      comm%messages_posted=.true.
  else
      comm%communication_complete=.true.
      comm%messages_posted=.false.
  end if




end subroutine post_p2p_communication


subroutine wait_p2p_communication(iproc, nproc, comm)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(p2pComms),intent(inout):: comm
  
  ! Local variables
  integer:: ierr
  
  
  if(.not.comm%communication_complete) then

      !!write(*,*) 'before fence window: iproc', iproc
      call mpi_win_fence(0, comm%window, ierr)
      !!write(*,*) 'after fence window: iproc', iproc
      call mpi_win_free(comm%window, ierr)

      !!if(.not.comm%messages_posted) stop 'ERROR: trying to wait for messages which have never been posted!'
      !!!!write(*,*) 'BEFORE WAIT SENDS: iproc', iproc

      !!! Wait for the sends to complete.
      !!if(comm%nsend>0) then
      !!    call mpi_waitall(comm%nsend, comm%requests(1,1), mpi_statuses_ignore, ierr)
      !!end if
      !!!!write(*,*) 'AFTER WAIT SENDS: iproc', iproc
 
      !!!!write(*,*) 'BEFORE WAIT RECEIVES: iproc', iproc
      !!! Wait for the receives to complete.
      !!if(comm%nrecv>0) then
      !!    call mpi_waitall(comm%nrecv, comm%requests(1,2), mpi_statuses_ignore, ierr)
      !!end if
      !!!!write(*,*) 'AFTER WAIT RECEIVES: iproc', iproc

  end if

  ! Flag indicating that the communication is complete
  comm%communication_complete=.true.

end subroutine wait_p2p_communication


!< Test whether the p2p communication hs completed. This is only called to make
!< sure that the communications really starts to execute.
subroutine test_p2p_communication(iproc, nproc, comm)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(p2pComms),intent(inout):: comm
  
  ! Local variables
  logical:: completed
  integer:: ierr
  

!!  !!write(10000+iproc,*) 'at beginning testing, iproc', iproc
!!  
!!  if(.not.comm%communication_complete) then
!!
!!      if(.not.comm%messages_posted) stop 'ERROR: trying to test for messages which have never been posted!'
!!
!!  !!write(20000+iproc,*) 'before testing sends, iproc', iproc
!!      ! Tests the sends.
!!      if(comm%nsend>0) then
!!          call mpi_testall(comm%nsend, comm%requests(1,1), completed, mpi_statuses_ignore, ierr)
!!      end if
!!  !!write(30000+iproc,*) 'after testing sends, iproc', iproc
!! 
!! 
!!  !!write(40000+iproc,*) 'before testing recvs, iproc', iproc
!!      ! Test the receives.
!!      if(comm%nrecv>0) then
!!          call mpi_testall(comm%nrecv, comm%requests(1,2), completed, mpi_statuses_ignore, ierr)
!!      end if
!!  !!write(50000+iproc,*) 'after testing recvs, iproc', iproc
!!
!!  end if


end subroutine test_p2p_communication



module p2p_tags_data
  use module_base
  implicit none
  logical,save:: initialized
  integer,dimension(:),allocatable,save:: tags
  integer,save:: tag_max
end module p2p_tags_data


subroutine init_p2p_tags(nproc)
  use module_base
  use module_types
  use p2p_tags_data
  implicit none

  ! Calling arguments
  integer,intent(in):: nproc
  character(len=*),parameter:: subname='init_p2p_tags'

  ! Local variables
  integer:: jproc, istat, ierr
  logical:: success

  if(initialized) stop 'trying to initialize the counter for the p2p tags which is already running!'

  allocate(tags(0:nproc-1),stat=istat)
  call memocc(istat,tags,'tags',subname)
  do jproc=0,nproc-1
      tags(jproc)=0
  end do
  initialized=.true.

  ! Determine the largest possible tag
  if (nproc > 1 ) then
     call mpi_attr_get(bigdft_mpi%mpi_comm, mpi_tag_ub, tag_max, success, ierr)
     if(.not.success) stop 'could not extract largest possible tag...'
  else
     tag_max=1
  end if
end subroutine init_p2p_tags


function p2p_tag(jproc)
  use module_base
  use p2p_tags_data
  implicit none

  ! Calling arguments
  integer,intent(in):: jproc
  integer:: p2p_tag

  if(.not.initialized) stop 'counter for tag was not properly initialized!'

  tags(jproc)=mod(tags(jproc)+1,tag_max)
  p2p_tag=tags(jproc)

end function p2p_tag


subroutine finalize_p2p_tags()
  use module_base
  use p2p_tags_data
  implicit none

  ! Local variables
  integer:: istat, iall
  character(len=*),parameter:: subname='finalize_p2p_tags'

  if(.not.initialized) stop 'trying to finalize the counter for the p2p tags which was not initialized!'

  iall=-product(shape(tags))*kind(tags)
  deallocate(tags,stat=istat)
  call memocc(istat,iall,'tags',subname)

  initialized=.false.

end subroutine finalize_p2p_tags
