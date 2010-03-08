!!****m* art/defs
!! FUNCTION
!!    This module defines all variables used accross the program ART01
!!
!! COPYRIGHT
!!    Copyright N. Mousseau, May 2001
!!    Copyright (C) 2010 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! SOURCE
!!
module defs

  implicit none

  real(8), parameter :: VERSION_NUMBER  = 1.5     ! Version of the code

  integer :: iproc, nproc ! MPI proc identificators

  real(8) :: TEMPERATURE        ! Temperature in eV
  integer :: NATOMS             ! Number of atoms in the system
  integer :: MAXNEI             ! Maximum number of nearest neighbours
  integer :: VECSIZE     ! Length of the force and position vectors

  integer :: NUMBER_EVENTS      ! Total number of events in this run
  logical :: NEW_EVENT          ! Total number of events in this run

  integer, parameter :: FCONF        = 1         ! Units for printing/reading
  integer, parameter :: FCOUNTER     = 2        
  integer, parameter :: FLIST        = 3         
  integer, parameter :: FLOG         = 4         
  integer, parameter :: FREFCONFIG   = 11
  integer, parameter :: FSTARTCONF   = 12
  integer, parameter :: FRESTART     = 9        
  integer, parameter :: XYZ          = 13  
  ! Name of the file storing the current configurations
  character(len=20) :: conf_initial, conf_saddle, conf_final
  
  integer, dimension(:), allocatable          :: type       ! Atomic type
  real(8), dimension(:), allocatable, target  :: force      ! Working forces on the atoms
  real(8), dimension(:), allocatable, target  :: pos        ! Working positions of the atoms
  real(8), dimension(:), allocatable, target  :: posref     ! Reference position
  real(8), dimension(:), allocatable, target  :: direction_restart  
  real(8), dimension(:), allocatable :: initial_direction  ! Initial move for leaving harmonic  well
  character(len=5), dimension(:), allocatable :: Atom
  character(len=5), dimension(5) :: type_name

  real(8), dimension(:), pointer :: x, y, z   ! Pointers for working position
  real(8), dimension(:), pointer :: xref, yref, zref   ! Pointers for reference position
  real(8), dimension(:), pointer :: fx, fy, fz   ! Pointers for working force

  real(8) :: PUSH_OVER   ! Fraction of displacement for pushing of saddle point.

  real(8) :: boxref                         ! Reference boxsize
  real(8) :: box                            ! Working boxsize

  real(8) :: scalaref                       ! Reference volume scaling
  real(8) :: scala                          ! Working volume scaling
  real(8) :: fscala                         ! Working forces on volume

  real(8) :: sym_break_dist                 ! Distance atoms are pushed to
                                            ! break symmetry

  logical :: restart                        ! State of restart (true or false)
  integer :: state_restart                  ! start of restart (1 - harmonic well, 2 - 
                                            ! activation, 3 - relaxation)

  integer :: preferred_atom                 ! Atom at the center of the event
  real(8) :: radius_initial_deformation     ! Radius of the initial deformation

  integer :: ievent                         ! actual number of events
  integer :: iter_restart                   ! iteraction number of restart
  integer :: ievent_restart                 ! Event number at restart
  integer :: mincounter                     ! Counter for output files
  integer :: refcounter                     ! Id of reference file

  integer :: evalf_number                   ! Number of force evalutions

  real(8) :: total_energy, ref_energy       ! Energies

  logical :: USE_DIIS                       ! Use DIIS for final convergence to saddle
  logical :: DIIS_CHECK_EIGENVEC            ! Check whether the metastable point is a saddle
  integer :: DIIS_MEMORY                    ! Number of steps kept in memory
  integer :: DIIS_MAXITER                   ! Maximum number of iterations
  real(8) :: DIIS_FORCE_THRESHOLD           ! Force Threshold to leave the algorithm
  real(8) :: DIIS_STEP                      ! Step used to update positions in DIIS


  character(len=20) :: LOGFILE
  character(len=20) :: EVENTSLIST
  character(len=20) :: REFCONFIG
  character(len=3)  :: FINAL
  character(len=3)  :: SADDLE 
  character(len=11) :: COUNTER
  character(len=11) :: RESTARTFILE
  character(len=20) :: eventtype 

  include 'mpif.h'

end module defs
!!***
