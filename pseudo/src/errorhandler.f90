subroutine errorhandler(ierr,iproc,nproc,errmsg)
! input ierr:
!             0  do nothing
!             1  write a NOTE    errmsg
!             2  write a WARNING errmsg
!             3  write a WARNING errmsg and STOP
implicit none
integer:: i,ierr,iproc,nproc
integer:: sendierr(nproc),getierr(nproc)
character(50):: errmsg
include 'mpif.h'

if(nproc>1)then
   !ierr: send to and receive from all processes
   sendierr=ierr
   call MPI_ALLTOALL(sendierr,1,MPI_INTEGER,  &
        getierr,1,MPI_INTEGER,MPI_COMM_WORLD,i)
   if(i/=0)   write(6,*)'Error in MPI_ALLTOALL occured-',i
else
   getierr=ierr 
end if


   if(any(getierr/=0)) write(6,*)
   if(any(getierr==1)) write(6,*)'             NOTE'
   if(any(getierr> 1)) write(6,*)'           WARNING'

if(any(getierr/=0))    write(6,*)             errmsg

if(nproc>1)then
   do i=1,nproc
      if(getierr(i)/=0)write(6,'(8x,a,i4)')'for process',i-1
   end do
   
   if(any(getierr==3))then
      write(6,*)'           EXITING'
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_FINALIZE(ierr)
      stop
   end if
else                        !serial case
   if(ierr==3)then
      write(6,*)
      write(6,*)'           EXITING'
      write(6,*)
      stop
   end if
end if

end subroutine

