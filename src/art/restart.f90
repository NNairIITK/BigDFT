!!****f* art/restart_states
!! FUNCTION
!!   This subroutine is used to restart an event during the activation. It uses the files:
!!     restart.dat
!!   It then continues the event where it stopped. 
!!
!! COPYRIGHT
!!    Copyright (C) 2001 Normand Mousseau
!!    Copyright (C) 2010 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! SOURCE
!! 
subroutine restart_states( istatus, ieventcurrent, iterations )

  use defs
  implicit none  

  !Arguments
  integer, intent(out) :: istatus
  integer, intent(out) :: ieventcurrent
  integer, intent(out) :: iterations 

  !Local variables
  integer :: i, j, ierror, ierr
  integer :: nat, id
  real(kind=8), dimension(:), pointer :: dx, dy, dz   ! Pointers for reference position

  logical :: exists_already
  real(kind=8),dimension(3) :: boxl
  character(len=100) :: fname
  character(len=100) :: f_xyz
  character(len=100) :: f_ascii
  character(len=150) :: commande
  character(len=9)   :: digit = "123456789"
  character(len=4), dimension(natoms) :: frzchain

  allocate(direction_restart(vecsize))
  dx => direction_restart(1:NATOMS)
  dy => direction_restart(NATOMS+1:2*NATOMS)
  dz => direction_restart(2*NATOMS+1:3*NATOMS)

  open(unit=FRESTART,file=RESTARTFILE,status='old',action='read',iostat=ierror)
  read(FRESTART,*) istatus, iterations, mincounter, refcounter, pas

  if ( istatus == 2 .or. istatus == 4 ) then 
     allocate(diis_forces_restart(DIIS_MEMORY,VECSIZE))
     allocate(diis_pos_restart(DIIS_MEMORY,VECSIZE))
     allocate(diis_norm_restart(DIIS_MEMORY))
     ! Initialization.
     diis_forces_restart = 0.0d0  
     diis_pos_restart    = 0.0d0
     diis_norm_restart   = 0.0d0
  end if 

  read(FRESTART,*) ieventcurrent, evalf_number, boundary
  read(FRESTART,*) box(1), box(2), box(3), scala
  read(FRESTART,*) boxref(1), boxref(2), boxref(3), scalaref
  read(FRESTART,*) total_energy, ref_energy
  read(FRESTART,*) ( typat(i), i = 1, NATOMS )
  read(FRESTART,*) ( xref(i), yref(i), zref(i), i = 1, NATOMS )
  read(FRESTART,*) ( x(i), y(i), z(i), i = 1, NATOMS )
  read(FRESTART,*) ( dx(i), dy(i), dz(i), i = 1, NATOMS )
  read(FRESTART,*) ( constr(i), i = 1, NATOMS )

  if ( istatus == 2 .or. istatus == 4 ) then 
     do j = 1, DIIS_MEMORY 
        dx => diis_forces_restart(j,1:NATOMS)
        dy => diis_forces_restart(j,NATOMS+1:2*NATOMS)
        dz => diis_forces_restart(j,2*NATOMS+1:3*NATOMS)
        read(FRESTART,*) ( dx(i), dy(i), dz(i), i = 1, NATOMS )
     end do
     do j = 1, DIIS_MEMORY 
        dx => diis_pos_restart(j,1:NATOMS)
        dy => diis_pos_restart(j,NATOMS+1:2*NATOMS)
        dz => diis_pos_restart(j,2*NATOMS+1:3*NATOMS)
        read(FRESTART,*) ( dx(i), dy(i), dz(i), i = 1, NATOMS )
     end do
    read(FRESTART,*) (diis_norm_restart(j), j = 1,DIIS_MEMORY)
    read(FRESTART,*) maxter_r,eigen_min_r,eigenvalue_r,nsteps_after_eigen_min_r 
  end if

  close(FRESTART)

  if ( iproc == 0 ) then
     ! Test posinp.xyz
     f_xyz = 'posinp.xyz' 
     fname = f_xyz 
     do i = 1, 9 
        inquire( file = fname, exist = exists_already )
        if ( exists_already ) then
           write(*,*) "BART: file ", trim(fname)," exists_already"
           fname = trim(f_xyz)//"."//digit(i:i) 
        else
           if ( i > 1 ) Then
              write(*,*) "BART: moving ",trim(f_xyz)," to ",trim(fname) 
              commande = "mv " //trim(f_xyz) // "  " //trim(fname)
              call system ( commande )
           end if
           exit
        end if
     end do
  end if

  if ( iproc == 0 ) then
     ! Test posinp.ascii
     f_ascii = 'posinp.ascii'
     fname = f_ascii
     do i = 1, 9 
        inquire( file = fname, exist = exists_already )
        if ( exists_already ) then
           write(*,*) "BART: file ", trim(fname)," exists_already"
           fname = trim(f_ascii)//"."//digit(i:i) 
        else
           if ( i > 1 ) Then
              write(*,*) "BART: moving ",trim(f_ascii)," to ",trim(fname) 
              commande = "mv " //trim(f_ascii)// "  " //trim(fname)
              call system ( commande )
           end if
           exit
        end if
     end do
  end if

  ! A security action
  if ( iproc == 0 ) then
     call system ( 'rm posinp.ascii' )
     call system ( 'rm posinp.xyz' ) 
  end if

  boxl = box * scala  ! Update the box size

  do i = 1, NATOMS
     Atom(i) = type_name(typat(i))
     if      ( constr(i)== 0) then
        frzchain(i)='    '
     else if (constr(i) == 1) then
        frzchain(i)='   f'
     else if (constr(i) == 2) then
        frzchain(i)='  fy'
     else if (constr(i) == 3) then
        frzchain(i)=' fxz'
     end if
  end do
 
  f_ascii = 'posinp.ascii'

  if (iproc == 0 ) then 
   write(*,*) "BART: writting new input file "
   open(unit=13,file=f_ascii,status='new',action='write',iostat=ierror)
   write(13,'(a,x,I5,x,a)') "# BART RESTART input file: ", NATOMS , 'angstroem' 
   write(13,*) boxl(1), 0.0, boxl(2)
   write(13,*)  0.0, 0.0, boxl(3)
   write(13,'(a)') "#keyword: angstroem"
   if (boundary == 'P') then
      write(13,'(a)') "#keyword: periodic"
   else if (boundary == 'S') then
      write(13,'(a)') "#keyword: surface"
   else
      write(13,'(a)') "#keyword: freeBC"
   end if
   do i=1, NATOMS
      write(13,'(1x,3(2x,f16.8),2x,A,2x,a4)')  x(i), y(i), z(i), Atom(i), frzchain(i)
   end do
   close(13)
  
                                      ! DEBUG 
   write(*,*) 'BART: restart file'
   write(*,*) 'BART: box: ', box
   write(*,'(a,(1p,e17.10,0p))') ' BART: total energy: ', total_energy
   write(*,'(a,3f20.8)') ' BART: xref(1): ', xref(1), yref(1), zref(1)
   write(*,'(a,3f20.8)') ' BART: x(1)   : ', x(1), y(1), z(1)
   write(*,'(a,3f20.8)') ' BART: x(NATOMS): ', x(NATOMS), y(NATOMS), z(NATOMS)
 
   if ( istatus == 2 .or. istatus == 4 ) then
   write(*,*) ' BART: previous_norm : ', diis_norm_restart
   write(*,*) ' BART: last line ',maxter_r,eigen_min_r,eigenvalue_r,nsteps_after_eigen_min_r 
   end if
  end if

END SUBROUTINE restart_states 
!!***


!!****f* restart/extra
!! FUNCTION
!!   saves the status for a restart
!!
!! SOURCE
!!
subroutine extra(ifrztyp,frzchain)
  implicit none
  integer, intent(in) :: ifrztyp
  character(len=4), intent(out) :: frzchain

  if (ifrztyp == 0) then
     frzchain='    '
  else if (ifrztyp == 1) then
     frzchain='   f'
  else if (ifrztyp == 2) then
     frzchain='  fy'
  else if (ifrztyp == 3) then
     frzchain=' fxz'
  end if
        
END SUBROUTINE extra 
!!***


!!****f* restart/save_state
!! FUNCTION
!!   saves the status for a restart
!!
!! SOURCE
!!
subroutine save_state( istatus, iter, direction )

  use defs
  implicit none  

  !Arguments
  integer, intent(in) :: istatus
  integer, intent(in) :: iter
  real(kind=8), dimension(VECSIZE), target :: direction

  !Local variables
  integer :: i, ierror
  real(kind=8), dimension(:), pointer :: dx, dy, dz   ! Pointers for reference position

  dx => direction(1:NATOMS)
  dy => direction(NATOMS+1:2*NATOMS)
  dz => direction(2*NATOMS+1:3*NATOMS)

  open(unit=FRESTART,file=RESTARTFILE,status='unknown',action='write',iostat=ierror)
  write(FRESTART,*) istatus, iter, mincounter, refcounter, pas+1
  write(FRESTART,*) ievent, evalf_number, boundary
  write(FRESTART,*) box(1), box(2), box(3), scala
  write(FRESTART,*) boxref(1), boxref(2), boxref(3), scalaref
  write(FRESTART,*) total_energy, ref_energy
  write(FRESTART,'(2x,10i6)') ( typat(i), i = 1, NATOMS )
  write(FRESTART,'(2x,3f20.8)') ( xref(i), yref(i), zref(i), i = 1, NATOMS )
  write(FRESTART,'(2x,3f20.8)') ( x(i), y(i), z(i), i = 1, NATOMS )
  write(FRESTART,'(2x,3f20.8)') ( dx(i), dy(i), dz(i), i = 1, NATOMS )
  write(FRESTART,'(2x,10i6)') ( constr(i), i = 1, NATOMS )
  close(FRESTART)

END SUBROUTINE save_state
!!***

!!****f* restart/save_state2
!! FUNCTION
!!   saves the status for a restart
!!
!! SOURCE
!!
subroutine save_state2( istatus, iter, direction, maxter, pf, px, pn, &
                        eigen_min, eigenvalue, nsteps_after_eigen_min   )

  use defs
  implicit none  

  !Arguments
  integer, intent(in) :: istatus
  integer, intent(in) :: iter
  real(kind=8), dimension(VECSIZE), target :: direction
  integer, intent(in) :: maxter 
  real(kind=8), dimension(DIIS_MEMORY,VECSIZE), target :: pf 
  real(kind=8), dimension(DIIS_MEMORY,VECSIZE), target :: px
  real(kind=8), dimension(DIIS_MEMORY), intent(in) :: pn
  real(kind=8), intent(in) :: eigen_min
  real(kind=8), intent(in) :: eigenvalue
  integer, intent(in) :: nsteps_after_eigen_min 

  !Local variables
  integer :: i, j, ierror
  real(kind=8), dimension(:), pointer :: dx, dy, dz   ! Pointers for reference position

  dx => direction(1:NATOMS)
  dy => direction(NATOMS+1:2*NATOMS)
  dz => direction(2*NATOMS+1:3*NATOMS)

  open(unit=FRESTART,file=RESTARTFILE,status='unknown',action='write',iostat=ierror)
  write(FRESTART,*) istatus, iter, mincounter, refcounter, pas+1
  write(FRESTART,*) ievent, evalf_number, boundary
  write(FRESTART,*) box(1), box(2), box(3), scala
  write(FRESTART,*) boxref(1), boxref(2), boxref(3), scalaref
  write(FRESTART,*) total_energy, ref_energy
  write(FRESTART,'(2x,10i6)') ( typat(i), i = 1, NATOMS )
  write(FRESTART,'(2x,3f20.8)') ( xref(i), yref(i), zref(i), i = 1, NATOMS )
  write(FRESTART,'(2x,3f20.8)') ( x(i), y(i), z(i), i = 1, NATOMS )
  write(FRESTART,'(2x,3f20.8)') ( dx(i), dy(i), dz(i), i = 1, NATOMS )
  write(FRESTART,'(2x,10i6)') ( constr(i), i = 1, NATOMS )
 
  do j = 1, DIIS_MEMORY 
     dx => pf(j,1:NATOMS)
     dy => pf(j,NATOMS+1:2*NATOMS)
     dz => pf(j,2*NATOMS+1:3*NATOMS)
     write(FRESTART,'(2x,3f20.8)') ( dx(i), dy(i), dz(i), i = 1, NATOMS )
  end do
  do j = 1, DIIS_MEMORY 
     dx => px(j,1:NATOMS)
     dy => px(j,NATOMS+1:2*NATOMS)
     dz => px(j,2*NATOMS+1:3*NATOMS)
     write(FRESTART,'(2x,3f20.8)') ( dx(i), dy(i), dz(i), i = 1, NATOMS )
  end do
  write(FRESTART,'(2x,3f20.8)') (pn(j), j = 1,DIIS_MEMORY)
  write(FRESTART,*) maxter, eigen_min, eigenvalue, nsteps_after_eigen_min 
  close(FRESTART)

END SUBROUTINE save_state2
!!***

!!****f* restart/save_input_file
!! FUNCTION
!!
!! SOURCE
!!
!subroutine save_input_file (file_p)
!
!  use defs
!  implicit none
!
!  !Arguments
!  character(len=*), intent(in) :: file_p
!
!  !Local variables
!  logical :: exists_already
!  character(len=100) :: fname
!  character(len=100) :: f_xyz
!  character(len=100) :: f_ascii
!  character(len=150) :: commande
!  character(len=9)   :: digit = "123456789"
!  integer :: i
!
!  integer :: ierror
!  real(kind=8),dimensio(3) :: boxl
!
!  if ( iproc == 0 ) then
!     ! Test posinp.xyz
!     f_xyz = trim(file_p)//".xyz"
!     fname = f_xyz 
!     do i = 1, 9 
!        inquire( file = fname, exist = exists_already )
!        write(*,*) "BART: exists_already", fname, exists_already
!        if ( exists_already ) then
!           fname = trim(f_xyz)//"."//digit(i:i) 
!        else
!           if ( i > 1 ) Then
!              commande = "mv " // f_xyz // "  " // fname
!              call system ( commande )
!           end if
!           exit
!        end if
!     end do
!  end if
!
!  !if ( iproc == 0 ) then
!  !   ! Test posinp.ascii
!  !   f_ascii = trim(file_p)//".ascii"
!  !   fname = f_ascii
!  !   do i = 1, 9 
!  !      inquire( file = fname, exist = exists_already )
!  !      write(*,*) "BART: ascii exists_already", fname, exists_already
!  !      if ( exists_already ) then
!  !         fname = trim(f_ascii)//"."//digit(i:i) 
!  !      write(*,*) "BART: fname", fname
!  !      else
!  !         if ( i > 1 ) Then
!  !            commande = "mv " //trim(f_ascii)// "  " //trim(fname)
!  !      write(*,*) "BART: commande", commande
!  !            call system ( commande )
!  !         end if
!  !         exit
!  !      end if
!  !   end do
!  !end if
!
!  boxl = box * scala  ! Update the box size
!
!  do i = 1, NATOMS
!     Atom(i) = type_name(typat(i))
!  end do
!
!  f_ascii = trim(file_p)//".ascii"
! 
! !if ( iproc == 0 ) then
! !inquire( file = f_ascii, exist = exists_already )
! !if ( exists_already ) then
! !commande = "rm " // f_ascii 
! !write(*,*) "BART: commande", commande
! !call system ( "rm posinp.ascii" )
! !end if
! !end if
!
!  if (iproc == 0 ) then 
!    write(*,*) "BART: new_input_file",  x(1), y(1), z(1), Atom(1)
!    write(*,*) "BART: new_input_file",  x(NATOMS), y(NATOMS), z(NATOMS), Atom(NATOMS)
!    open(unit=14,file=f_ascii,status='new',action='write',iostat=ierror)
!    write(*,*) "BART: f_ascii, ierror" ,  f_ascii, ierror
!    write(14,*) "# BART RESTART input file: ", NATOMS , 'angstroem' 
!    write(14,*) boxl, 0.0, boxl
!    write(14,*)  0.0, 0.0, boxl
!    write(14,*) "#keyword: angstroem"
!
!    do i=1, NATOMS
!       write(14,'(1x,3(2x,f16.8),2x,A)')  x(i), y(i), z(i), Atom(i)
!    end do
!
!    close(14)
!  end if
!
!end subroutine save_input_file 
!!***

!!****f* restart/new_input_file
!! FUNCTION
!!
!! SOURCE
!!
!subroutine new_input_file( file_p )
!
!  use defs
!  implicit none
!
!  !Arguments
!  character(len=*), intent(in) :: file_p
!
!  !Local variables
!  integer :: ierror
!  integer :: i
!  real(kind=8) :: boxl
!  character(len=100) :: f_ascii
!
!  boxl = box * scala  ! Update the box size
!
!  do i = 1, NATOMS
!     Atom(i) = type_name(typat(i))
!  end do
!
!  f_ascii = trim(file_p)//".ascii"
!
!  if (iproc == 0 ) then 
!    write(*,*) "BART: new_input_file",  x(1), y(1), z(1), Atom(1)
!    write(*,*) "BART: new_input_file",  x(NATOMS), y(NATOMS), z(NATOMS), Atom(NATOMS)
!    open(unit=13,file=f_ascii,status='new',action='write',iostat=ierror)
!    write(*,*) "BART  f_ascii, ierror" ,  f_ascii, ierror
!    write(13,*) "# BART RESTART input file: ", NATOMS , 'angstroem' 
!    write(13,*) boxl, 0.0, boxl
!    write(13,*)  0.0, 0.0, boxl
!    write(13,*) "#keyword: angstroem"
!
!    do i=1, NATOMS
!       write(13,'(1x,3(2x,f16.8),2x,A)')  x(i), y(i), z(i), Atom(i)
!    end do
!
!    close(13)
!  end if
!
!END SUBROUTINE new_input_file 
