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
  integer :: MAXKTER
  integer :: MAXIPERP, MAXKPERP
  integer, dimension(:), allocatable  :: atom_displaced  ! Id of local atoms displaced
  real(kind=8), dimension(:), allocatable, target :: dr

  real(kind=8) :: INITSTEPSIZE 
  real(kind=8) :: LOCAL_CUTOFF
  real(kind=8) :: INCREMENT,BASIN_FACTOR 
  real(kind=8) :: FTHRESHOLD
  real(kind=8) :: EIGEN_THRESH 
  real(kind=8) :: EXITTHRESH

  character(len=20) :: TYPE_EVENTS

  real(kind=8) :: delta_e  ! current_energy - ref_energy
  integer :: m_perp   ! Accepted perpendicular iterations.

  integer :: try      ! Total # of iterationes in a given hyperplane
  real(kind=8) :: ftot     ! Norm of the total force...
  real(kind=8) :: fpar     ! Parallel projection of force
  real(kind=8) :: fperp    ! Norm of the perpendicular force
  real(kind=8) :: delr     ! Distance between two configurations.
  integer :: npart    ! number of particles having moved by more 
                      ! than a THRESHOLD
  logical :: switchDIIS
  logical :: end_activation
  integer :: nsteps_after_eigen_min 
  real(kind=8) :: eigen_min

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

  if ( ( .not. restart ) .and. new_event ) then 

     evalf_number = 0                 ! Initialization of counter.
                                      ! Broadcast pos to all nodes.
                                      ! IS IT NECCESARY ?
     nat = 3 * NATOMS
     call MPI_Bcast(posref,nat,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
     call MPI_Bcast(scalaref,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
     call MPI_Bcast(boxref,3,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
                                      ! We copy the reference half box 
                                      ! and scaling to the working ones.

                                      ! We start from the reference
                                      ! configuration.
     scala = scalaref
     box   = boxref
     pos   = posref  
                                      ! These  subroutines modify the vector pos
                                      ! and generate a vector of length 1 indicating
                                      ! the direction of the random displacement.
     if ( TYPE_EVENTS == 'global' ) then
        call global_move( )
     else if ( TYPE_EVENTS == 'local' ) then 
        call symmetry_break( )        !Breaks the initial symmetry.
        call local_move( )
     else if ( TYPE_EVENTS == 'list' ) then 
        call symmetry_break( )
        call list_of_atoms ( )
     else if ( TYPE_EVENTS == 'list_local' ) then 
        call symmetry_break( )
        call list_and_local ( )
     end if

  end if
                                      ! Now, activate per se.
  call saddle_converge( ret, saddle_energy )
                                      ! If the activation did not converge, for
                                      ! whatever reason, we restart the routine
                                      ! and do not accept the new position.
  call end_report( success, ret, saddle_energy )

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
  integer :: i, j, that, ierror
  real(kind=8) :: lcutoff2    ! Cut-off for local moves, squared
  real(kind=8) :: dr2
  real(kind=8) :: xi, yi, zi, xij, yij, zij
  real(kind=8) :: ran3
  real(kind=8), dimension(3) :: boxl, invbox
  real(kind=8), dimension(:), pointer  :: dx, dy, dz

                                      ! Select an atom at random.
  if ( preferred_atom < 0 ) then
                                      ! Only between totally free atoms 
     if ( iproc == 0 ) then
        do 
          that = int( NATOMS * ran3() + 1 ) 
          if ( constr(that) == 0 ) exit 
        end do
     end if
     call MPI_Bcast(that, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
  else 
     that = preferred_atom
     if ( constr(that) .ne. 0 ) then  ! die !
        write(*,*) ' ERROR: choosen atom is blocked in the geometry file '
        call end_art()                          
     end if 
  end if
                                      ! Write
  if ( iproc == 0 ) then
   open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
       & action = 'write', position = 'append', iostat = ierror )
   write(FLOG,*) ' '
   write(FLOG,'(1X,A34,I17)') ' - That atom                    : ', that
   close(FLOG)
   write(*,*) 'BART: That atom = ', that
  end if

  boxl = box*scala                    ! without periodic boundary conditions box
                                      ! is zero so we set invbox to 1 
  do i = 1, 3
     if ( boxl(i) .lt. 1.d-08 ) then 
        invbox(i) = 1.0d0
     else
        invbox(i) = 1.0d0 / boxl(i)
     end if
  end do
                                      ! Square the cut-off
  lcutoff2 = LOCAL_CUTOFF * LOCAL_CUTOFF  

  allocate(dr(3*natoms))
  allocate(atom_displaced(natoms))
                                      ! We assign a few pointers 
  dx => dr(1:NATOMS)
  dy => dr(NATOMS+1:2*NATOMS)
  dz => dr(2*NATOMS+1:3*NATOMS)

  dr = 0.0d0  
  natom_displaced = 0
  atom_displaced  = 0  

  if ( iproc == 0 ) then              ! Work only on the master node.
                                      ! Now we also displace all atoms within
                                      ! a cut-off distance, LOCAL_CUTOFF.
     xi = x(that)
     yi = y(that)
     zi = z(that)
     do j = 1, NATOMS
        if ( constr(j) == 0 ) then
           xij = x(j) - xi - boxl(1) * nint((x(j)-xi) * invbox(1))
           yij = y(j) - yi - boxl(2) * nint((y(j)-yi) * invbox(2))
           zij = z(j) - zi - boxl(3) * nint((z(j)-zi) * invbox(3))
           
           dr2 = xij*xij + yij*yij + zij*zij
           if ( dr2 < lcutoff2 ) then ! Close enough, give a random displacement.
              do
                dx(j) = 0.5d0 - ran3()
                dy(j) = 0.5d0 - ran3()
                dz(j) = 0.5d0 - ran3()
                                      ! Ensures that the random
                                      ! displacement is isotropic
                dr2 = dx(j)**2 + dy(j)**2 + dz(j)**2
                if ( dr2 < 0.25d0 ) exit  
              end do
              natom_displaced = natom_displaced + 1
              atom_displaced(j) = 1
           end if
        end if
     end do
   end if

   call center_and_norm ( INITSTEPSIZE )

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
  integer :: i
  real(kind=8) :: dr2
  real(kind=8) :: ran3
  real(kind=8), dimension(:), pointer :: dx, dy, dz

  allocate(dr(3*natoms))
  allocate(atom_displaced(natoms))
                                      ! We assign a few pointers. 
  dx => dr(1:NATOMS)
  dy => dr(NATOMS+1:2*NATOMS)
  dz => dr(2*NATOMS+1:3*NATOMS)

  atom_displaced = 0 
  natom_displaced = 0 
  dr = 0.0d0 

  if ( iproc == 0 ) then              ! Only on the master node.
                                      ! Generate a random displacement.
     do i = 1, natoms, 1
        if ( constr(i) == 0 ) then
           do
             dx(i) = 0.5d0 - ran3()
             dy(i) = 0.5d0 - ran3()
             dz(i) = 0.5d0 - ran3()
                                      ! Ensures that the random
                                      ! displacement is isotropic
             dr2 = dx(i)**2 + dy(i)**2 + dz(i)**2
             if ( dr2 < 0.25d0 ) exit 
           end do
           natom_displaced = natom_displaced + 1
           atom_displaced(i) = 1
        end if
     end do
  end if

  call center_and_norm ( INITSTEPSIZE )

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
  integer :: i
  real(kind=8) :: dr2
  real(kind=8) :: ran3
  real(kind=8), dimension(:), pointer :: dx, dy, dz

  allocate(dr(3*natoms))
  allocate(atom_displaced(natoms))
                                      ! We assign a few pointers.
  dx => dr(1:NATOMS)
  dy => dr(NATOMS+1:2*NATOMS)
  dz => dr(2*NATOMS+1:3*NATOMS)

  atom_displaced = 0 
  natom_displaced = 0 
  dr = 0.0d0                          ! Initialization of dr.

  if ( iproc == 0 ) then              ! Only on the master node.
                                      ! Generate a random displacement. 
     do i = 1, NATOMS
        if ( constr(i) == 0 ) then 
           do
             dx(i) = 0.5d0 - ran3()
             dy(i) = 0.5d0 - ran3()
             dz(i) = 0.5d0 - ran3()
                                      ! Ensures that the random
                                      ! displacement is isotropic
             dr2 = dx(i)**2 + dy(i)**2 + dz(i)**2
             if ( dr2 < 0.25d0 ) exit 
           end do
           natom_displaced = natom_displaced + 1
           atom_displaced(i) = 1
        end if
     end do
  end if

  call center_and_norm ( sym_break_dist )

END SUBROUTINE symmetry_break
!!***


!!****f* find_sanddle/list_of_atoms
!! SOURCE
!! 
subroutine list_of_atoms ( ) 

  use defs
  use random
  use saddles
  implicit none

  !Local variables  
  logical :: found
  integer :: i, j, i_stat, nlines, sizelist
  integer, dimension(:), allocatable :: list_atoms
  character(len = 128)                            :: filename
  character(len = 150)                            :: line
  character(len = 150), dimension(:), allocatable :: lines
  real(kind=8)                        :: ran3, dr2
  real(kind=8), dimension(:), pointer :: dx, dy, dz

  found    = .false.
  filename = 'list_atoms.dat'
  filename = adjustl( filename )

  ! Is there the file?  
  inquire( file = filename, exist = found )
  if ( .not. found ) then
     write(*,'(a)') 'ERROR: list_atoms.dat not found !! '
     call end_art() 
  end if
  
  open( unit= 99, file= filename, status= 'old' )
   
  ! First pass: to store of the file in a string buffer. 
  ! WARNING: Dimension of lines is 500.
  allocate (lines(500))
  nlines = 1
  do
     read( 99,'(a150)', iostat = i_stat) lines(nlines)
     if (i_stat /= 0) then
        exit
     end if
     nlines = nlines + 1
     if ( nlines > 500 ) then
        write(*,*) 'list_atoms.dat file too long (> 500 lines).'
        write(*,*) 'change this parameter in the source!!'
        call end_art() 
     end if
  end do
  nlines = nlines - 1
  close (99)
 
  if ( nlines < 2 ) then
     write(*,*) 'ERROR: list_atoms, file has less than 2 lines.'
     call end_art() 
  end if
 
  ! Second pass: to determine the correct number of atoms.
  sizelist = 0 
  do i = 1, nlines, 1
     write(line, "(a150)") adjustl(lines(i))
     if (len(trim(line)) /= 0 ) then 
        sizelist = sizelist + 1 
     end if
  end do

  allocate( list_atoms ( sizelist ) )

  ! Last pass. We store the atoms in list_atoms array.
  list_atoms = 0
  do i = 1, sizelist, 1
     write(line, "(a150)") adjustl(lines(i))
     read(line,*, iostat = i_stat)  list_atoms(i) 
  end do
  deallocate(lines)

  allocate(dr(3*natoms))
  allocate(atom_displaced(natoms))
                                      ! We assign a few pointers
  dx => dr(1:NATOMS)
  dy => dr(NATOMS+1:2*NATOMS)
  dz => dr(2*NATOMS+1:3*NATOMS)

  atom_displaced = 0 
  natom_displaced = 0 
  dr = 0.0d0  
                                      ! Now we apply a random displacement 
  if ( iproc == 0 ) then              ! Work only on the master node.
     do j = 1, NATOMS
        if ( constr(j) == 0 ) then
           do i = 1, sizelist 

              if ( j == list_atoms(i) ) then
                 do
                   dx(j) = 0.5d0 - ran3()
                   dy(j) = 0.5d0 - ran3()
                   dz(j) = 0.5d0 - ran3()

                   dr2 = dx(j)**2 + dy(j)**2 + dz(j)**2
                                                 ! Ensures that the random
                                                 ! displacement is isotropic
                   if ( dr2 < 0.25d0 ) exit
                 end do
                 natom_displaced = natom_displaced + 1
                 atom_displaced(j) = 1
                 exit
              end if     

           end do
        end if
     end do
  end if
  
  deallocate(list_atoms)

  call center_and_norm ( INITSTEPSIZE )

END SUBROUTINE list_of_atoms 
!!***


!!****f* find_sanddle/list_and_local
!! SOURCE
!! 
subroutine list_and_local () 

  use defs
  use random
  use saddles
  implicit none

  !Local variables  
  logical :: found
  integer :: i, j, i_stat, nlines, sizelist
  integer :: that, this, ierror
  integer, dimension(:), allocatable :: list_atoms
  character(len = 128)                            :: filename
  character(len = 150)                            :: line
  character(len = 150), dimension(:), allocatable :: lines
  real(kind=8) :: lcutoff2    ! Cut-off for local moves, squared
  real(kind=8) :: ran3, dr2
  real(kind=8) :: xi, yi, zi, xij, yij, zij
  real(kind=8), dimension(3) :: boxl, invbox
  real(kind=8), dimension(:), pointer :: dx, dy, dz

  found    = .false.
  filename = 'list_atoms.dat'
  filename = adjustl( filename )

  ! Is there the file?  
  inquire( file = filename, exist = found )
  if ( .not. found ) then
     write(*,'(a)') 'ERROR: list_atoms.dat not found !! '
     call end_art()  
  end if
  
  open( unit= 99, file= filename, status= 'old' )
   
  ! First pass: to store of the file in a string buffer. 
  ! WARNING: Dimension of lines is 500.
  allocate (lines(500))
  nlines = 1
  do
     read( 99,'(a150)', iostat = i_stat) lines(nlines)
     if (i_stat /= 0) then
        exit
     end if
     nlines = nlines + 1
     if ( nlines > 500 ) then
        write(*,*) 'list_atoms.dat file too long (> 500 lines).'
        write(*,*) 'change this parameter in the source!!'
        call end_art() 
     end if
  end do
  nlines = nlines - 1
  close (99)
 
  if ( nlines < 2 ) then
     write(*,*) 'ERROR: list_atoms, file has less than 2 lines.'
     call end_art() 
  end if
 
  ! Second pass: to determine the correct number of atoms.
  sizelist = 0 
  do i = 1, nlines, 1
     write(line, "(a150)") adjustl(lines(i))
     if (len(trim(line)) /= 0 ) then 
        sizelist = sizelist + 1 
     end if
  end do

  allocate( list_atoms ( sizelist ) )

  ! Last pass. We store the atoms in list_atoms array.
  list_atoms = 0
  do i = 1, sizelist, 1
     write(line, "(a150)") adjustl(lines(i))
     read(line,*, iostat = i_stat)  list_atoms(i) 
  end do
  deallocate(lines)

  if ( iproc == 0 ) then
     do 
       this = int( sizelist * ran3() + 1 ) 
       that = list_atoms(this)
       if ( constr(that) == 0 ) exit 
     end do
  end if
  call MPI_Bcast(that, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)

                                      ! Write
  if ( iproc == 0 ) then
   open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
       & action = 'write', position = 'append', iostat = ierror )
   write(FLOG,*) ' '
   write(FLOG,'(1X,A34,I17)') ' - That atom                    : ', that
   close(FLOG)
   write(*,*) 'BART: That atom = ', that
  end if

  boxl = box*scala                    ! without periodic boundary conditions box
                                      ! is zero so we set invbox to 1 
  do i = 1, 3
     if ( boxl(i) .lt. 1.d-08 ) then 
        invbox(i) = 1.0d0
     else
        invbox(i) = 1.0d0 / boxl(i)
     end if
  end do
                                      ! Square the cut-off
  lcutoff2 = LOCAL_CUTOFF * LOCAL_CUTOFF  

  allocate(dr(3*natoms))
  allocate(atom_displaced(natoms))
                                      ! We assign a few pointers 
  dx => dr(1:NATOMS)
  dy => dr(NATOMS+1:2*NATOMS)
  dz => dr(2*NATOMS+1:3*NATOMS)

  dr = 0.0d0  
  natom_displaced = 0
  atom_displaced  = 0  

  if ( iproc == 0 ) then              ! Work only on the master node.
                                      ! Now we also displace all atoms within
                                      ! a cut-off distance, LOCAL_CUTOFF.
     xi = x(that)
     yi = y(that)
     zi = z(that)
     do j = 1, NATOMS
        if ( constr(j) == 0 ) then
           xij = x(j) - xi - boxl(1) * nint((x(j)-xi) * invbox(1))
           yij = y(j) - yi - boxl(2) * nint((y(j)-yi) * invbox(2))
           zij = z(j) - zi - boxl(3) * nint((z(j)-zi) * invbox(3))
           
           dr2 = xij*xij + yij*yij + zij*zij
           if ( dr2 < lcutoff2 ) then ! Close enough, give a random displacement.
              do
                dx(j) = 0.5d0 - ran3()
                dy(j) = 0.5d0 - ran3()
                dz(j) = 0.5d0 - ran3()
                                      ! Ensures that the random
                                      ! displacement is isotropic
                dr2 = dx(j)**2 + dy(j)**2 + dz(j)**2
                if ( dr2 < 0.25d0 ) exit  
              end do
              natom_displaced = natom_displaced + 1
              atom_displaced(j) = 1
           end if
        end if
     end do
   end if

   call center_and_norm ( INITSTEPSIZE )

END SUBROUTINE list_and_local  
!!***


!!****f* find_sanddle/center_and_norm
!! SOURCE
!! 
subroutine center_and_norm ( step )

  use defs
  use saddles
  implicit none

  !Arguments
  real(kind=8), intent(in) :: step

  !Local variables 
  integer :: nat, i
  real(kind=8) :: ierror
  real(kind=8) :: xsum, ysum, zsum
  real(kind=8) :: xnorm, ynorm, znorm, norm
  real(kind=8), dimension(:), pointer :: dx, dy, dz

  dx => dr(1:NATOMS)
  dy => dr(NATOMS+1:2*NATOMS)
  dz => dr(2*NATOMS+1:3*NATOMS)

  if ( iproc == 0 ) then              ! Work only on the master node.   
                                      ! We now center only on the atoms that have
                                      ! been randomly selected.
     if ( natom_displaced > 1 ) then
        xsum = sum(dx) / natom_displaced
        ysum = sum(dy) / natom_displaced
        zsum = sum(dz) / natom_displaced
        
        dx = dx - xsum * atom_displaced
        dy = dy - ysum * atom_displaced
        dz = dz - zsum * atom_displaced
     end if
                                      ! And normalize the total random displacement 
                                      ! effected to the value desired.
     norm = 0.0d0
     norm = dot_product(dr, dr)
                                      ! This renormalizes in angstroems to a
                                      ! displacement INITSTEPSIZE.     
     norm = step / sqrt(norm)
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
  pos = pos + dr  
                                      ! Now, we normalize dr to get the initial_direction
                                      ! (note that this had to be done after the transfer
                                      ! into box units.
  initial_direction = dr 
  norm = dot_product( initial_direction, initial_direction ) 
  norm = 1.0d0 / sqrt(norm)
  initial_direction  = initial_direction * norm

  if ( iproc == 0 ) write(*,*) 'BART: Number of displaced atoms initially: ',natom_displaced

  deallocate(atom_displaced)
  deallocate(dr)

END SUBROUTINE center_and_norm 
