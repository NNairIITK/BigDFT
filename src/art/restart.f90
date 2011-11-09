!> @file
!! @author
!!    Copyright (C) 2001 Normand Mousseau
!!    Copyright (C) 2010-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> ART restart_states
!!   This subroutine is used to restart an event during the activation. It uses the files:
!!   restartfile 
subroutine restart_states( istatus, ieventcurrent, iterations, iatp )

  use defs
  implicit none  

  !Arguments
  integer, intent(out) :: istatus
  integer, intent(out) :: ieventcurrent
  integer, intent(out) :: iterations 
  integer, intent(out) :: iatp 

  !Local variables
  integer :: i, j, ierror
  real(kind=8), dimension(:), pointer :: dx, dy, dz   ! Pointers for reference position

  allocate(direction_restart(vecsize))
  dx => direction_restart(1:NATOMS)
  dy => direction_restart(NATOMS+1:2*NATOMS)
  dz => direction_restart(2*NATOMS+1:3*NATOMS)

  open(unit=FRESTART,file=RESTARTFILE,status='old',action='read',iostat=ierror)
  read(FRESTART,*) istatus, iterations, mincounter, refcounter, pas, iatp

  if ( istatus == 2 .or. istatus == 4 ) then 
     allocate(diis_forces_restart(DIIS_MEMORY,VECSIZE))
     allocate(diis_pos_restart(DIIS_MEMORY,VECSIZE))
     allocate(diis_norm_restart(DIIS_MEMORY))
     ! Initialization.
     diis_forces_restart = 0.0d0  
     diis_pos_restart    = 0.0d0
     diis_norm_restart   = 0.0d0
  end if 

  read(FRESTART,*) ieventcurrent, evalf_number, boundary, central_atom 
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

  do i = 1, NATOMS
     Atom(i) = type_name(typat(i))
  end do
                                      ! DEBUG 
  if ( iproc == 0 ) then
     write(*,*) 'BART: restart file'
     write(*,*) 'BART: box: ', box
     write(*,*) 'BART: central_atom ', central_atom  
     write(*,'(a,(1p,e17.10,0p))') ' BART: total energy: ', total_energy
     write(*,'(a,(1p,e17.10,0p))') ' BART: ref   energy: ', ref_energy 
     write(*,'(a,3f20.8)') ' BART: xref(1): ', xref(1), yref(1), zref(1)
     write(*,'(a,3f20.8)') ' BART: x(1)   : ', x(1), y(1), z(1)
     write(*,'(a,3f20.8)') ' BART: x(NATOMS): ', x(NATOMS), y(NATOMS), z(NATOMS)
 
     if ( istatus == 2 .or. istatus == 4 ) then
     write(*,*) 'BART: previous_norm : ', diis_norm_restart
     write(*,*) 'BART: last line ',maxter_r,eigen_min_r,eigenvalue_r,nsteps_after_eigen_min_r 
     end if 
     write(*,*) 'BART: Reading atomic input positions from RESTART file,'
     write(*,*) 'BART: please, ignore the next message. '
  end if

END SUBROUTINE restart_states 


!> ART save_state
!!   Saves the status for a restart
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
  write(FRESTART,*) istatus, iter, mincounter, refcounter, pas+1, atp
  write(FRESTART,*) ievent, evalf_number, boundary, central_atom
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


!> ART save_state2
!!   Saves the status for a restart
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
  write(FRESTART,*) istatus, iter, mincounter, refcounter, pas+1, atp
  write(FRESTART,*) ievent, evalf_number, boundary, central_atom
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
