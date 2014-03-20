!> @file
!!  Define routine for onesided communications (linear version)
!! @author
!!    Copyright (C) 2013-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

 
subroutine start_onesided_communication(iproc, nproc, nsendbuf, sendbuf, nrecvbuf, recvbuf, comm, lzd)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer, intent(in):: iproc, nproc, nsendbuf, nrecvbuf
  real(kind=8), dimension(nsendbuf), intent(in):: sendbuf
  real(kind=8), dimension(nrecvbuf), intent(out):: recvbuf
  type(p2pComms), intent(inout):: comm
  type(local_zone_descriptors), intent(in) :: lzd
  
  ! Local variables
  !character(len=*), parameter :: subname='start_onesided_communication'
  integer :: jproc, joverlap, mpisource, istsource, mpidest, istdest, ierr, nit
  integer :: ioffset_send, mpi_type, ist, i2, i3, ist2, ist3, info, nsize, size_of_double



  if(.not.comm%communication_complete) stop 'ERROR: there is already a p2p communication going on...'

  nproc_if: if (nproc>1) then

      ! Allocate MPI memory window
      call mpi_type_size(mpi_double_precision, size_of_double, ierr)
      call mpi_info_create(info, ierr)
      call mpi_info_set(info, "no_locks", "true", ierr)
      call mpi_win_create(sendbuf(1), int(nsendbuf*size_of_double,kind=mpi_address_kind), size_of_double, &
           info, bigdft_mpi%mpi_comm, comm%window, ierr)
      call mpi_info_free(info, ierr)

      call mpi_win_fence(mpi_mode_noprecede, comm%window, ierr)
      
      do jproc=0,nproc-1
          do joverlap=1,comm%noverlaps(jproc)
              mpisource=comm%comarr(1,joverlap,jproc)
              istsource=comm%comarr(2,joverlap,jproc)
              mpidest=comm%comarr(3,joverlap,jproc)
              istdest=comm%comarr(4,joverlap,jproc)
              nit=comm%comarr(5,joverlap,jproc)
              ioffset_send=comm%comarr(6,joverlap,jproc)
              call mpi_type_create_hvector(nit, 1, int(size_of_double*ioffset_send,kind=mpi_address_kind), &
                   comm%mpi_datatypes(1,jproc), comm%mpi_datatypes(2,jproc), ierr)
              call mpi_type_commit(comm%mpi_datatypes(2,jproc), ierr)
              if (iproc==mpidest) then
                  call mpi_type_size(comm%mpi_datatypes(2,jproc), nsize, ierr)
                  nsize=nsize/size_of_double
                  if(nsize>0) then
                      call mpi_get(recvbuf(istdest), nsize, &
                           mpi_double_precision, mpisource, int((istsource-1),kind=mpi_address_kind), &
                           1, comm%mpi_datatypes(2,jproc), comm%window, ierr)
                  end if
              end if
          end do
      end do

  else nproc_if

      ist=1
      do i3=comm%ise(5,iproc),comm%ise(6,iproc)
          ist3=(i3-1)*lzd%glr%d%n1i*lzd%glr%d%n2i
          do i2=comm%ise(3,iproc),comm%ise(4,iproc)
              ist2=(i2-1)*lzd%glr%d%n1i
              call dcopy(comm%ise(2,iproc)-comm%ise(1,iproc)+1, sendbuf(ist3+ist2+1), 1, recvbuf(ist), 1)
              ist=ist+comm%ise(2,iproc)-comm%ise(1,iproc)+1
          end do
      end do

  end if nproc_if
  
  
  ! Flag indicating whether the communication is complete or not
  if(nproc>1) then
      comm%communication_complete=.false.
  else
      comm%communication_complete=.true.
  end if


end subroutine start_onesided_communication


subroutine synchronize_onesided_communication(iproc, nproc, comm)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(p2pComms),intent(inout):: comm
  
  ! Local variables
  integer:: ierr, jproc
  
  
  if(.not.comm%communication_complete) then
      call mpi_win_fence(0, comm%window, ierr)
      do jproc=0,nproc-1
          call mpi_type_free(comm%mpi_datatypes(2,jproc), ierr)
      end do
      call mpi_win_free(comm%window, ierr)
  end if

  ! Flag indicating that the communication is complete
  comm%communication_complete=.true.

end subroutine synchronize_onesided_communication
