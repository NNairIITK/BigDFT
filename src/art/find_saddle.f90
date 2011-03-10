!> art/saddles
!!
!! @author
!!    Copyright (C) 2001 Normand Mousseau
!!    Copyright (C) 2010 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!! 
module saddles

  use defs
  implicit none

  character(len=20) :: TYPE_EVENTS

  real(kind=8) :: INITSTEPSIZE 
  real(kind=8) :: LOCAL_CUTOFF

  integer :: KTER_MIN, NVECTOR_LANCZOS
  integer :: MAXITER, MAXKTER
  integer :: MAXIPERP, MAXKPERP
  real(kind=8) :: INCREMENT,BASIN_FACTOR 
  real(kind=8) :: FTHRESHOLD, FTHRESH2
  real(kind=8) :: EIGEN_THRESH 
  real(kind=8) :: EXITTHRESH

  integer, dimension(:), allocatable  :: atom_displaced     ! Id of local atoms displaced
  integer                     :: natom_displaced    ! # of local atoms displaced
end module saddles



!>  This subroutine initiates the random displacement at the start
!!  of the ART algorithm. 
!!  After  random escape direction has been selected, the routine call 
!!  saddle_converge which will try to bring the configuration to a saddle point.
!!  If the convergence fails, find_saddle will restart; if it succeeds, it saves 
!!  the event and returns to the main loop.
!!  The displacement can be either LOCAL and involve a single atom
!!  and its neareast neighbours or
!!
!!  NON-LOCAL and involve ALL the atoms.
!!
!!  For large cells, it is preferable to use a local initial
!!  displacement to prevent the creation of many trajectories
!!  at the same time in different sections of the cell.
!! 
subroutine find_saddle(success)
  use random
  use defs
  use saddles
  use lanczos_defs
  implicit none

  logical, intent(out) :: success
  integer :: ret, npart, ierror, nat
  real(kind=8) :: boxl
  real(kind=8) :: fperp, fpar, del_r, saddle_energy
  character(len=4)  :: scounter
  character(len=20) :: fname


  allocate(atom_displaced(natoms))       

  if ( (.not. restart) .and. new_event ) then 

     evalf_number = 0

     ! Broadcast pos to all nodes
     nat = 3 * NATOMS
     call MPI_Bcast(posref,nat,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
     call MPI_Bcast(scalaref,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
     call MPI_Bcast(boxref,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)

     ! We first copy the reference position, half box and scaling to the working ones
     scala = scalaref
     box  = boxref


     ! And copy the reference position into the working vector
     pos = posref  !Vectorial operation

     ! We now define the correct box units
     boxl = box * scala

     ! These two subroutines modify the vector pos and generate a vector of length 1
     ! indicating the direction of the random displacement
     if (TYPE_EVENTS .eq. 'global') then
        call global_move()
     else
        call symmetry_break()   !Breaks the initial symmetry
        call local_move()
     endif
  endif

  ! Now, activate per se
  call saddle_converge(ret, saddle_energy, fpar, fperp)

  ! We compute the displacement (del_r) and number of involved particles (npart)
  call displacement(pos, posref, del_r,npart)
  
  ! We write out various information to both screen and file
  if (iproc .eq. 0 ) then
     write(*,"(' ','Total energy S: ',f16.4,'  npart: ', i4,'  delr: ',&
          & f12.6,'  fpar: ',f12.6,'  fperp: ',f12.6,'  ret: ',i6,' force eval: ',&
          & i6)") saddle_energy, npart,del_r,fpar,fperp,ret,evalf_number
     
     open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
     write(FLOG,"(' ','Total energy S: ',f16.4,'  npart: ', i4,'  delr: ',&
          & f12.6,'  fpar: ',f12.6,'  fperp: ',f12.6,'  ret: ',i6,' force eval: ',&
          & i6)") saddle_energy, npart,del_r,fpar,fperp,ret,evalf_number
  endif

  ! If the activation did not converge, for whatever reason, we restart the 
  ! routine and do not accept the new position

  if ( ( ret < 0) .or. (ret > 30000) ) then
    success = .false.
  else
    success = .true.
    ! We need to write the configuration in a sad.... file  
    call convert_to_chain(mincounter,scounter)
    if (iproc .eq. 0 ) write(*,*) ' Mincounter is : ', mincounter, ' and scounter is ', scounter
    fname =   SADDLE // scounter
    conf_saddle = fname
  
    ! We now store the configuration into fname
    if (iproc .eq. 0 ) then
       call store(fname)
  
       write(*,*) 'Configuration stored in file ',fname
       write(FLOG,*) 'Configuration stored in file ',fname
    endif
  endif
 
  if (iproc .eq. 0 ) close(FLOG)

  deallocate(atom_displaced)       

END SUBROUTINE find_saddle



!>   The initial random direction is taken from a restricted space based on 
!!   the local bonding environment. For this, we need to know the list of neighbours
!!   and the cut-off of the potential. Other approaches could be used also.
!!
!! 
subroutine local_move()
  use defs
  use random
  use saddles
  implicit none

  integer                               :: i, j, that, i_id, j_id, nat, ierror
  real(kind=8)                          :: lcutoff2    ! Cut-off for local moves, squared
  real(kind=8), dimension(VECSIZE), target :: dr
  real(kind=8), dimension(:), pointer   :: dx, dy, dz
  real(kind=8)                          :: boxl, invbox
  real(kind=8)                          :: xsum, ysum, zsum, xnorm, ynorm, znorm, norm
  real(kind=8) :: ran3, Numb,dr2,xi,xij,yi,yij,zi,zij

  ! We assign a few pointers 
  dx => dr(1:NATOMS)
  dy => dr(NATOMS+1:2*NATOMS)
  dz => dr(2*NATOMS+1:3*NATOMS)

  boxl = box*scala
  invbox = 1.0 / boxl
  if (iproc .eq. 0 ) Numb = ran3()
  call MPI_Bcast(Numb, 1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)

  ! Square the cut-off
  lcutoff2 = LOCAL_CUTOFF * LOCAL_CUTOFF   ! Vectorial operation

  ! Select an atom at random in this case it must be an atom having coordination number 3 
  if (preferred_atom .lt. 0 ) then
    if (iproc .eq. 0) that = int ( NATOMS * ran3() + 1)    ! Between 1 and NATOMS
    call MPI_Bcast(that, 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
  else 
    that = preferred_atom
  endif
  if (iproc .eq. 0) write(*,*) 'That = ', that

  if (iproc .eq. 0) then
     open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
     write(FLOG,*) 'That = ', that
     close(flog)
  endif

  dr = 0.0d0  ! Vectorial operation
  atom_displaced = 0  !Vectorial operation

  if (iproc .eq. 0 ) then     ! Work only on the master node
     do 
        dx(that) = 0.5d0 - ran3()
        dy(that) = 0.5d0 - ran3()
        dz(that) = 0.5d0 - ran3()
        
        dr2 = dx(that)**2 + dy(that)**2 + dz(that)**2
        if (dr2 < 0.25d0 ) exit  ! Ensures that the random displacement is isotropic
     end do

     natom_displaced = 1 
     atom_displaced(that) = 1
     
     i_id = type(that)
     ! Now we also displace all atoms within a cut-off distance, LOCAL_CUTOFF
     xi = x(that)
     yi = y(that)
     zi = z(that)
     do j=1, NATOMS
        j_id = type(j)
        xij = x(j) - xi - boxl * nint((x(j)-xi) * invbox)
        yij = y(j) - yi - boxl * nint((y(j)-yi) * invbox)
        zij = z(j) - zi - boxl * nint((z(j)-zi) * invbox)
        
        dr2 = xij*xij + yij*yij + zij*zij
        
        if(dr2 < lcutoff2 ) then  ! Close enough, give a random displacement
           do
              dx(j) = 0.5d0 - ran3()
              dy(j) = 0.5d0 - ran3()
              dz(j) = 0.5d0 - ran3()
              
              dr2 = dx(j)**2 + dy(j)**2 + dz(j)**2
              if (dr2 < 0.25d0 ) exit  ! Ensure that the random displacement is isotropic
           end do
           natom_displaced = natom_displaced + 1
           atom_displaced(j) = 1
        endif
     end do
     
     ! We now center only on the atoms that have been randomly selected
     xsum = sum(dx) / natom_displaced
     ysum = sum(dy) / natom_displaced
     zsum = sum(dz) / natom_displaced
     
     dx = dx - xsum * atom_displaced
     dy = dy - ysum * atom_displaced
     dz = dz - zsum * atom_displaced
     
     ! And  normalize the total random displacement effected to the value desired
     norm = dot_product(dr, dr)
     
     ! This renormalizes in angstroems to a displacement INITSTEPSIZE
     norm = INITSTEPSIZE / sqrt(norm)
     xnorm = norm 
     ynorm = norm
     znorm = norm
     
     ! The displacement is now in box units
     dx = dx * xnorm  ! Vectorial operation
     dy = dy * ynorm
     dz = dz * znorm
     
  endif
  
  ! Redistribute the information
  nat = 3*NATOMS
  call MPI_Bcast(dr,nat,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
  
  ! Update the position using this random displacement
  pos = pos +  dr  ! Vectorial operation
  
  ! Now, we normalize dr to get the initial_direction (note that this had to be
  ! done after the transfer into box units
  
  initial_direction = dr    ! Vectorial operation
  norm = 0.0d0
  do i=1, VECSIZE
     norm = norm + initial_direction(i) * initial_direction(i)
  end do
  norm = 1.0 / sqrt (norm)
  initial_direction  = initial_direction * norm
  
  if (iproc .eq. 0 ) write(*,*) 'Number of displaced atoms initially: ',natom_displaced
END SUBROUTINE local_move



!>   The initial random direction is taken from the full 3N-dimensional space
!!
!!
subroutine global_move()
  use defs
  use random
  use saddles
  implicit none

  integer :: i,ierror,nat
  real(kind=8) :: norm, xnorm, ynorm, znorm
  real(kind=8), dimension(VECSIZE), target :: dr
  real(kind=8), dimension(:), pointer    :: dx, dy, dz
  real(kind=8) :: ran3

  ! We assign a few pointers 
  dx => dr(1:NATOMS)
  dy => dr(NATOMS+1:2*NATOMS)
  dz => dr(2*NATOMS+1:3*NATOMS)

  ! All atoms are displaced
  atom_displaced = 1  ! Vectorial operation

  ! Generate a random displacement
  if (iproc .eq. 0) then  ! Only on the master node
     do i=1, VECSIZE
        dr(i) = 0.5d0 - ran3()
     end do

     ! Keep the center of mass fixed 
     call center(dr,VECSIZE)
  
     ! And renormalize the total displacement to the value desired
     norm = 0.0d0
     do i=1, VECSIZE
        norm = norm + dr(i) * dr(i)
     end do

     !  This renormalizes in angstroems to a displacement INITSTEPSIZE
     norm = INITSTEPSIZE / sqrt(norm)

     ! We need to go back to x, y, z in case the three lengths are different
     xnorm = norm
     ynorm = norm
     znorm = norm

     ! The displacement is now in box units
     dx = dx * xnorm  ! Vectorial operation
     dy = dy * ynorm
     dz = dz * znorm

  endif

  ! Redistribute the information
  nat = 3*NATOMS
  call MPI_Bcast(dr,nat,MPI_REAL8,0,MPI_COMM_WORLD,ierror)

  ! Update the position using this random displacement
  pos = pos +  dr  ! Vectorial operation

  ! Now, we normalize dr to get the initial_direction (note that this had to be
  ! done after the transfer into box units

  initial_direction = dr    ! Vectorial operation
  norm = 0.0d0
  do i=1, VECSIZE
     norm = norm + initial_direction(i) * initial_direction(i)
  end do
  norm = 1.0 / sqrt (norm)
  initial_direction  = initial_direction * norm

END SUBROUTINE global_move



!> art/symmetry_break
!!
!! 
subroutine symmetry_break()
  use defs
  use random
  use saddles
  implicit none

  integer                               :: i,ierror,nat
  real(kind=8)                          :: lcutoff2    ! Cut-off for local moves, squared
  real(kind=8), dimension(VECSIZE), target :: dr
  real(kind=8), dimension(:), pointer   :: dx, dy, dz
  real(kind=8)                          :: dr2
  real(kind=8)                          :: boxl, invbox
  real(kind=8)                          :: xsum, ysum, zsum, xnorm, ynorm, znorm, norm
  real(kind=8) :: ran3, Numb

  ! We assign a few pointers 
  dx => dr(1:NATOMS)
  dy => dr(NATOMS+1:2*NATOMS)
  dz => dr(2*NATOMS+1:3*NATOMS)

  boxl = box*scala
  invbox = 1.0 / boxl
  Numb = ran3()

  ! Square the cut-off
  LOCAL_CUTOFF = 0.5
  lcutoff2 = LOCAL_CUTOFF * LOCAL_CUTOFF   ! Vectorial operation

  if (iproc .eq. 0 ) then
     open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
     write(FLOG,*) 'Starting symmetry break' 
     write(FLOG,*) 'Push distance per atom:', sym_break_dist
     close(flog)
  endif


  if (iproc .eq. 0 ) then
        
     dr = 0.0d0  ! Vectorial operation
     atom_displaced = 0  !Vectorial operation
     natom_displaced = 0 
     
     do i=1, NATOMS
        do
           dx(i) = 0.5d0 - ran3()
           dy(i) = 0.5d0 - ran3()
           dz(i) = 0.5d0 - ran3()
           
           dr2 = dx(i)**2 + dy(i)**2 + dz(i)**2
           if (dr2 < 0.25d0 ) exit  ! Ensure that the random displacement is isotropic
        end do
        natom_displaced = natom_displaced + 1
        atom_displaced(i) = 1
     end do
  
     ! We now center only on the atoms that have been randomly selected
     xsum = sum(dx) / natom_displaced
     ysum = sum(dy) / natom_displaced
     zsum = sum(dz) / natom_displaced
     
     dx = dx - xsum * atom_displaced
     dy = dy - ysum * atom_displaced
     dz = dz - zsum * atom_displaced
     
     ! And  normalize the total random displacement effected to the value desired
     norm = dot_product(dr, dr)
     
     ! This renormalizes in angstroems to a displacement sym_break_dist
     norm = sym_break_dist / sqrt(norm)
     xnorm = norm 
     ynorm = norm
     znorm = norm
     
     ! The displacement is now in box units
     dx = dx * xnorm  ! Vectorial operation
     dy = dy * ynorm
     dz = dz * znorm
     
  endif

  ! Redistribute the information
  nat = 3*NATOMS
  call MPI_Bcast(dr,nat,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
  
  ! Update the position using this random displacement
  pos = pos +  dr  ! Vectorial operation

  ! Now, we normalize dr to get the initial_direction (note that this had to be
  ! done after the transfer into box units
  
  initial_direction = dr    ! Vectorial operation
  norm = 0.0d0
  do i=1, VECSIZE
     norm = norm + initial_direction(i) * initial_direction(i)
  end do
  norm = 1.0 / sqrt (norm)
  initial_direction  = initial_direction * norm
  
  if (iproc .eq. 0 ) then 
     open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
     write(FLOG,*) 'Done breaking symmetry'
     close(flog)
  endif

 END SUBROUTINE symmetry_break
 
