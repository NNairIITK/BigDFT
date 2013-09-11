!> @file
!! atomic program for generating and optimizing HGH pseudo-potentials.
!! @author
!!    Alex Willand, under the supervision of Stefan Goedecker
!!    gpu accelerated routines by Raffael Widmer
!!    parts of this program were based on the fitting program by matthias krack
!!    http://cvs.berlios.de/cgi-bin/viewcvs.cgi/cp2k/potentials/goedecker/pseudo/v2.2/
!!
!!    Copyright (C) 2010-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Error handler for pseudo program
!! input ierr:
!!  0  do nothing
!!  1  write a NOTE    errmsg
!!  2  write a WARNING errmsg
!!  3  write a WARNING errmsg and STOP
subroutine errorhandler(ierr,iproc,nproc,errmsg)
   implicit none
   !Arguments
   integer, intent(in) :: ierr,iproc,nproc
   character(len=*), intent(in) :: errmsg
   !Local variables
   integer :: i
   integer :: sendierr(nproc),getierr(nproc)
   include 'mpif.h'

   if (nproc>1) then
      !ierr: send to and receive from all processes
      sendierr=ierr
      call MPI_ALLTOALL(sendierr,1,MPI_INTEGER,  &
           getierr,1,MPI_INTEGER,MPI_COMM_WORLD,i)
      if ( i /= 0)   write(6,'(1x,a,i0)') 'Error in MPI_ALLTOALL occured- ',i
   else
      getierr=ierr 
   end if


   if (any(getierr /= 0)) write(6,*)
   if (any(getierr == 1)) write(6,'(/,12x,a,/)') 'NOTE'
   if (any(getierr > 1))  write(6,'(/,12x,a,/)') 'WARNING'

   if (any(getierr /= 0)) write(6,'(1x,a)') errmsg

   if (nproc>1) then
      do i=1,nproc
         if (getierr(i) /= 0) write(6,'(8x,a,i4)')'for process',i-1
      end do
      
      if (any(getierr == 3)) then
         write(6,'(/,12x,a,/)') 'EXITING'
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         call MPI_FINALIZE(ierr)
         stop
      end if
   else                        !serial case
      if (ierr == 3) then
         write(6,'(/,12x,a,/)') 'EXITING'
         stop
      end if
   end if

end subroutine
