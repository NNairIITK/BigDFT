!!****m* art/saddles
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
module saddles

  implicit none
  
  integer :: natom_displaced    ! # of local atoms displaced
  integer :: KTER_MIN
  integer :: MAXITER, MAXKTER
  integer :: MAXIPERP, MAXKPERP
  integer, dimension(:), allocatable  :: atom_displaced  ! Id of local atoms displaced

  real(kind=8) :: INITSTEPSIZE 
  real(kind=8) :: LOCAL_CUTOFF
  real(kind=8) :: INCREMENT,BASIN_FACTOR 
  real(kind=8) :: FTHRESHOLD
  real(kind=8) :: EIGEN_THRESH 
  real(kind=8) :: EXITTHRESH

  character(len=20) :: TYPE_EVENTS

end module saddles
!!***


!!****f* find_saddle/find_saddle
!! FUNCTION
!!   This subroutine initiates the random displacement at the start
!!   of the ART algorithm. 
!!   After a random escape direction has been selected, the routine call 
!!   saddle_converge which will try to bring the configuration to a saddle point.
!!   If the convergence fails, find_saddle will restart; if it succeeds, it saves 
!!   the event and returns to the main loop.
!!   The displacement can be either LOCAL and involve a single atom
!!   and its neareast neighbours or...
!!
!!   GLOBAL and involve ALL the atoms.
!!
!!   For large cells, it is preferable to use a local initial
!!   displacement to prevent the creation of many trajectories
!!   at the same time in different sections of the cell.
!!
!! SOURCE
!! 
subroutine find_saddle( success, saddle_energy )

  use defs
  use saddles, only : TYPE_EVENTS 
  implicit none

  !Arguments
  logical, intent(out) :: success
  real(kind=8), intent(out) :: saddle_energy

  !Local variables
  integer :: ret, ierror, nat
  character(len=4)  :: scounter
  character(len=20) :: fname

  if ( ( .not. restart ) .and. new_event ) then 

     evalf_number = 0
                                     ! Broadcast pos to all nodes.
                                     ! IS IT NECCESARY ?
                                     
     nat = 3 * NATOMS
     call MPI_Bcast(posref,nat,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
     call MPI_Bcast(scalaref,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
     call MPI_Bcast(boxref,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
                                      ! We copy the reference half box 
                                      ! and scaling to the working ones.
     scala = scalaref
     box   = boxref
                                      ! And copy the reference position
                                      ! into the working vector.
     pos   = posref  
                                      ! These two subroutines modify the vector pos
                                      ! and generate a vector of length 1 indicating
                                      ! the direction of the random displacement.
     if ( TYPE_EVENTS == 'global' ) then
        call global_move( )
     else
        call symmetry_break( )        !Breaks the initial symmetry.
        call local_move( )
     end if
  end if

                                      ! Now, activate per se.
  call saddle_converge( ret, saddle_energy )

                                      ! If the activation did not converge, for
                                      ! whatever reason, we restart the routine
                                      ! and do not accept the new position.
  If_ret: if ( ( ret < 0 ) .or. ( ret > 30000 ) ) then

     success = .false.
                                      ! Write
     if ( iproc == 0 ) then
      open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
          & action = 'write', position = 'append', iostat = ierror )
      write(FLOG,*) ' Activation has not converged !!'
      write(FLOG,'(1X,A34,I17)'), ' - ret                          : ', ret
      close(FLOG)
     end if

  else                                ! Else If_ret

     success = .true.

                                      ! Write
     if ( iproc == 0 ) then
      open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
          & action = 'write', position = 'append', iostat = ierror )
      write(FLOG,'(1X,A34,(1p,e17.10,0p))') &
      &  ' - Total energy Saddle (eV)     : ', saddle_energy
      write(FLOG,'(1X,A34,I17)'), ' - ret                          : ', ret
      close(FLOG)
      write(*,'(1X,A,(1p,e17.10,0p))') 'BART: Total energy Saddle : ', saddle_energy
      write(*,'(1X,A,I12)'), 'BART: ret : ', ret
     end if

                                      ! We write the configuration in a sad.... file  
     if ( iproc == 0 ) call convert_to_chain( mincounter, scounter )
     fname = SADDLE // scounter
     if ( iproc == 0 ) call store( fname ) 
     conf_saddle = fname

                                      ! Write  
     if ( iproc == 0 ) then
      open( unit = FLOG, file = LOGFILE, status = 'unknown',&
          & action = 'write', position = 'append', iostat = ierror )
      write(*,*) 'BART: Configuration stored in file ',fname
      write(FLOG,'(1X,A34,A17)') ' - Configuration stored in file : ', trim(fname)
      write(FLOG,*) ' '
      close(FLOG)
     end if

  end if If_ret
 
END SUBROUTINE find_saddle
!!***


!!****f* find_saddle/local_move
!! FUNCTION
!!   The initial random direction is taken from a restricted space based on 
!!   the local bonding environment. For this, we need to know the list of neighbours
!!   and the cut-off of the potential. Other approaches could be used also.
!! SOURCE
!! 
subroutine local_move( )

  use defs
  use random
  use saddles
  implicit none

  !Local variables
  integer                               :: i, j, that, i_id, j_id, nat, ierror
  real(kind=8)                          :: lcutoff2    ! Cut-off for local moves, squared
  real(kind=8), dimension(VECSIZE), target :: dr
  real(kind=8), dimension(:), pointer   :: dx, dy, dz
  real(kind=8)                          :: dr2
  real(kind=8)                          :: xi, yi, zi, xij, yij, zij
  real(kind=8)                          :: boxl, invbox
  real(kind=8)                          :: xsum, ysum, zsum, xnorm, ynorm, znorm, norm
  real(kind=8)                          :: ran3

  allocate(atom_displaced(natoms))

                                      ! We assign a few pointers 
  dx => dr(1:NATOMS)
  dy => dr(NATOMS+1:2*NATOMS)
  dz => dr(2*NATOMS+1:3*NATOMS)

  boxl = box*scala
  invbox = 1.0d0 / boxl

                                      ! Square the cut-off
  lcutoff2 = LOCAL_CUTOFF * LOCAL_CUTOFF  

                                      ! Select an atom at random.
  if ( preferred_atom < 0 ) then
                                      ! Between 1 and NATOMS
     if ( iproc == 0 ) that = int( NATOMS * ran3() + 1 ) 
     call MPI_Bcast(that, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
  else 
     that = preferred_atom
  end if

                                      ! Write
  if ( iproc == 0 ) then
   open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
       & action = 'write', position = 'append', iostat = ierror )
   write(FLOG,'(1X,A34,I17)') ' - That atom                    : ', that
   write(FLOG,*) ' '
   close(FLOG)
   write(*,*) 'BART: That atom = ', that
  end if

  dr = 0.0d0  
  atom_displaced = 0  

  if ( iproc == 0 ) then              ! Work only on the master node.
     do  
        dx(that) = 0.5d0 - ran3()
        dy(that) = 0.5d0 - ran3()
        dz(that) = 0.5d0 - ran3()
        
        dr2 = dx(that)**2 + dy(that)**2 + dz(that)**2
                                      ! Ensures that the random displacement
                                      ! is isotropic
        if ( dr2 < 0.25d0 ) exit  
     end do

     natom_displaced = 1 
     atom_displaced(that) = 1
     i_id = type(that)

                                      ! Now we also displace all atoms within
                                      ! a cut-off distance, LOCAL_CUTOFF.
     xi = x(that)
     yi = y(that)
     zi = z(that)
     do j = 1, NATOMS
        j_id = type(j)
        xij = x(j) - xi - boxl * nint((x(j)-xi) * invbox)
        yij = y(j) - yi - boxl * nint((y(j)-yi) * invbox)
        zij = z(j) - zi - boxl * nint((z(j)-zi) * invbox)
        
        dr2 = xij*xij + yij*yij + zij*zij
        
        if ( dr2 < lcutoff2 ) then    ! Close enough, give a random displacement.
           do
              dx(j) = 0.5d0 - ran3()
              dy(j) = 0.5d0 - ran3()
              dz(j) = 0.5d0 - ran3()
              
              dr2 = dx(j)**2 + dy(j)**2 + dz(j)**2
                                      ! Ensures that the random displacement
                                      ! is isotropic
              if ( dr2 < 0.25d0 ) exit  
           end do
           natom_displaced = natom_displaced + 1
           atom_displaced(j) = 1
        end if
     end do
                                      ! We now center only on the atoms that have been
                                      ! randomly selected.
     xsum = sum(dx) / natom_displaced
     ysum = sum(dy) / natom_displaced
     zsum = sum(dz) / natom_displaced
     
     dx = dx - xsum * atom_displaced
     dy = dy - ysum * atom_displaced
     dz = dz - zsum * atom_displaced
                                      ! And normalize the total random displacement 
                                      ! effected to the value desired.
     norm = 0.0d0
     norm = dot_product(dr, dr)
                                      ! This renormalizes in angstroems to a
                                      ! displacement INITSTEPSIZE.     
     norm = INITSTEPSIZE / sqrt(norm)
     xnorm = norm 
     ynorm = norm
     znorm = norm

                                      ! The displacement is now in box units. 
     dx = dx * xnorm 
     dy = dy * ynorm
     dz = dz * znorm
     
  end if
  
                                      ! Redistribute the information.
  nat = 3*NATOMS
  call MPI_Bcast(dr,nat,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
  
                                      ! Update the position using this random displacement
  pos = pos +  dr  
                                      ! Now, we normalize dr to get the initial_direction
                                      ! (note that this had to be done after the transfer
                                      ! into box units.
  initial_direction = dr 
  norm = 0.0d0
  norm = dot_product( initial_direction, initial_direction ) 
  norm = 1.0 / sqrt (norm)
  initial_direction  = initial_direction * norm
  
  if ( iproc == 0 ) write(*,*) 'BART: Number of displaced atoms initially: ',natom_displaced

  deallocate(atom_displaced)

END SUBROUTINE local_move
!!***


!!****f* find_saddle/global_move
!! FUNCTION
!!   The initial random direction is taken from the full 3N-dimensional space
!! SOURCE
!!
subroutine global_move( )

  use defs
  use random
  use saddles
  implicit none

  !Local variables
  integer :: i, nat
  real(kind=8) :: norm, xnorm, ynorm, znorm, ierror
  real(kind=8), dimension(VECSIZE), target :: dr
  real(kind=8), dimension(:), pointer    :: dx, dy, dz
  real(kind=8) :: ran3

  allocate(atom_displaced(natoms))

                                      ! We assign a few pointers. 
  dx => dr(1:NATOMS)
  dy => dr(NATOMS+1:2*NATOMS)
  dz => dr(2*NATOMS+1:3*NATOMS)

                                      ! All atoms are displaced.
  atom_displaced = 1 
  natom_displaced = NATOMS
 
                                      ! Generate a random displacement.
  if ( iproc == 0 ) then              ! Only on the master node.
     do i = 1, VECSIZE
        dr(i) = 0.5d0 - ran3()
     end do

                                      ! Keep the center of mass fixed.
     call center(dr,VECSIZE)
  
                                      ! And renormalize the total displacement
                                      ! to the value desired.
     norm = 0.0d0
     norm = dot_product( dr, dr ) 
                                      ! This renormalizes in angstroems to a
                                      ! displacement INITSTEPSIZE.
     norm = INITSTEPSIZE / sqrt(norm)
                                      ! We need to go back to x, y, z in case
                                      ! the three lengths are different.
     xnorm = norm
     ynorm = norm
     znorm = norm

                                      ! The displacement is now in box units.
     dx = dx * xnorm  
     dy = dy * ynorm
     dz = dz * znorm

  end if

                                      ! Redistribute the information.
  nat = 3*NATOMS
  call MPI_Bcast(dr,nat,MPI_REAL8,0,MPI_COMM_WORLD,ierror)

                                      ! Update the position using this random displacement
  pos = pos +  dr  
                                      ! Now, we normalize dr to get the initial_direction
                                      ! (note that this had to be done after the transfer
                                      ! into box units.
  initial_direction = dr 
  norm = 0.0d0
  norm = dot_product( initial_direction, initial_direction ) 
  norm = 1.0 / sqrt (norm)
  initial_direction  = initial_direction * norm

  deallocate(atom_displaced)

END SUBROUTINE global_move
!!***


!!****f* find_sanddle/symmetry_break
!! SOURCE
!! 
subroutine symmetry_break( )

  use defs
  use random
  use saddles
  implicit none

  !Local variables
  integer                                  :: i, nat, ierror
  real(kind=8)                             :: lcutoff2 ! Cut-off for local moves, squared
  real(kind=8), dimension(VECSIZE), target :: dr
  real(kind=8), dimension(:), pointer      :: dx, dy, dz
  real(kind=8)                             :: dr2
  real(kind=8)                             :: boxl, invbox
  real(kind=8)                             :: xsum, ysum, zsum, xnorm, ynorm, znorm, norm
  real(kind=8)                             :: ran3

  allocate(atom_displaced(natoms))

                                      ! We assign a few pointers.
  dx => dr(1:NATOMS)
  dy => dr(NATOMS+1:2*NATOMS)
  dz => dr(2*NATOMS+1:3*NATOMS)

  boxl = box*scala
  invbox = 1.0d0 / boxl

                                      ! Square the cut-off.
  LOCAL_CUTOFF = 0.5
  lcutoff2 = LOCAL_CUTOFF * LOCAL_CUTOFF 

                                      ! Write
  if ( iproc == 0 ) then
   open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
       & action = 'write', position = 'append', iostat = ierror )
   write(FLOG,*) ' Starting symmetry break' 
   write(FLOG,'(1X,A34,F17.6)') ' - Push distance per atom       : ', sym_break_dist
   close(FLOG)
  end if

  if ( iproc == 0 ) then              ! Only on the master node.
        
     dr = 0.0d0 
     atom_displaced = 0 
     natom_displaced = 0 
     
     do i = 1, NATOMS
        do
           dx(i) = 0.5d0 - ran3()
           dy(i) = 0.5d0 - ran3()
           dz(i) = 0.5d0 - ran3()
           
           dr2 = dx(i)**2 + dy(i)**2 + dz(i)**2
           if ( dr2 < 0.25d0 ) exit   ! Ensure that the random displacement is isotropic.
        end do
        natom_displaced = natom_displaced + 1
        atom_displaced(i) = 1
     end do
                                      ! We now center only on the atoms that have been
                                      ! randomly selected.
     xsum = sum(dx) / natom_displaced
     ysum = sum(dy) / natom_displaced
     zsum = sum(dz) / natom_displaced
     
     dx = dx - xsum * atom_displaced
     dy = dy - ysum * atom_displaced
     dz = dz - zsum * atom_displaced
                                      ! And normalize the total random displacement 
                                      ! effected to the value desired.
     norm = dot_product(dr, dr)
                                      ! This renormalizes in angstroems to a 
                                      ! displacement sym_break_dist.
     norm = sym_break_dist / sqrt(norm)
     xnorm = norm 
     ynorm = norm
     znorm = norm

                                      ! The displacement is now in box units.
     dx = dx * xnorm 
     dy = dy * ynorm
     dz = dz * znorm
     
  endif

                                      ! Redistribute the information.
  nat = 3*NATOMS
  call MPI_Bcast(dr,nat,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
  
                                      ! Update the position using this random displacement
  pos = pos +  dr 
                                      ! Now, we normalize dr to get the initial_direction
                                      ! (note that this had to be done after the transfer
                                      ! into box units.
  initial_direction = dr 
  norm = 0.0d0
  norm = dot_product( initial_direction, initial_direction ) 
  norm = 1.0 / sqrt (norm)
  initial_direction  = initial_direction * norm

                                      ! Write  
  if (iproc == 0 ) then 
   open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
       & action = 'write', position = 'append', iostat = ierror )
   write(FLOG,*) ' Done breaking symmetry'
   close(FLOG)
  end if

  deallocate(atom_displaced)

END SUBROUTINE symmetry_break
!!***
