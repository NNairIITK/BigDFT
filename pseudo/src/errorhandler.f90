subroutine errorhandler(ierr,iproc,nproc,errmsg)
   ! input ierr:
   !             0  do nothing
   !             1  write a note    errmsg
   !             2  write a warning errmsg
   !             3  write a warning errmsg and stop
   implicit none
   integer:: i,ierr,iproc,nproc
   integer:: sendierr(nproc),getierr(nproc)
   character(50):: errmsg
   include 'mpif.h'
   
   if(nproc>1)then
      !ierr: send to and receive from all processes
      sendierr=ierr
      call mpi_alltoall(sendierr,1,mpi_integer,  &
           getierr,1,mpi_integer,mpi_comm_world,i)
      if(i/=0)   write(6,*)'error in mpi_alltoall occured-',i
   else
      getierr=ierr 
   end if
   
   
   if(any(getierr/=0)) write(6,*)
   if(any(getierr==1)) write(6,*)'             note'
   if(any(getierr> 1)) write(6,*)'           warning'
   
   if(any(getierr/=0))    write(6,*)             errmsg
   
   if(nproc>1)then
      do i=1,nproc
         if(getierr(i)/=0)write(6,'(8x,a,i4)')'for process',i-1
      end do
      
      if(any(getierr==3))then
         write(6,*)'           exiting'
         call mpi_barrier(mpi_comm_world,ierr)
         call mpi_finalize(ierr)
         stop
      end if
   else                        !serial case
      if(ierr==3)then
         write(6,*)
         write(6,*)'           exiting'
         write(6,*)
         stop
      end if
   end if

end subroutine errorhandler
