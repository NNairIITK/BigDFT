! This subroutine writes the atomic positions and others to a "refconfig" file
! which will be used a the reference point until a new events gets accepted
!
!  
! Copyright Normand Mousseau May 2001

subroutine write_refconfig()
   use defs
   implicit none
   integer :: i,ierror
   real(8) :: boxl

   boxl = box * scala  ! Update the box size 

   open(unit=FREFCONFIG,file=REFCONFIG,status='unknown',action='write',iostat=ierror) !switch replace for unknown
   write(FREFCONFIG,*) 'run_id: ', mincounter
   write(FREFCONFIG,*) 'total energy : '
   write(FREFCONFIG,*) total_energy
   write(FREFCONFIG,*) boxl
   do i=1, NATOMS
     write(FREFCONFIG,'(1x,i6, 3(2x,F16.8))') type(i),x(i),y(i), z(i)
   enddo
   close(FREFCONFIG)

   return
end subroutine
