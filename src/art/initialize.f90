!> @file
!! @author
!!    Copyright (C) 2001 Normand Mousseau
!!    Copyright (C) 2010 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!! Modified by Laurent Karim Beland, UdeM, 2011.
!!
!> ART initialize
!!   Initialization of art method
!!   It relaxes it into a local minimum without volume optimization
!!   The initial configuration is a "reference configuration", it will be the reference 
!!   configuration until a new event is accepted. 
!! 
!!   The reference configuration should contain the following information:
!! 
!!     first line:         refcounter 1232      
!!     second line:        total_energy: -867.33454
!!     third line:         8.34543 8.21345 8.67789
!!     next NATOMS lines:  0  0.0437844 0.96894 0.847555
!!     
!!     The first number is the configuration number associated with the reference 
!!     configuration. 
!! 
!!     The second line gives the total energy so that it can be compared with a 
!!     direct calculation.
!!     
!!     The third line indicates the full box size along the x, y and z directions.
!! 
!!     The NATOMS lines are the atomic species and coordinates in Angstroems
subroutine initialize()

  use defs
  use bigdft_forces, only : init_all_atoms
  use lanczos_defs, only: LANCZOS_MIN 
  implicit none

  !Local variables
  integer :: i, ierror
  character(len=20) :: dummy, fname
  character(len=4)  :: scounter
  logical           :: flag, success

  integer                   :: nat_test
  integer, pointer          :: typa(:)    ! Atomic type
  integer, pointer          :: const_(:)  ! Constraints
  real(kind=8), pointer     :: posa(:)    ! Working positions of the atoms
  real(kind=8),dimension(3) :: boxref_    ! Reference box from posinp file
  !_______________________
  !! WARNING EDUARDO: cuidado, aqui pierdo el refcounter de REFCONFIG

  !open(unit=FREFCONFIG,file=REFCONFIG,status='old',action='read',iostat=ierror)
  !read(FREFCONFIG, '(A10,i5)' ) dummy, refcounter
  !read(FREFCONFIG,*) dummy, ref_energy !! M-A Malouin
  !read(FREFCONFIG,*) boundary,boxref(1), boxref(2), boxref(3)
  !do i = 1,natoms
  !	read(FREFCONFIG,*) (typat(i),x(i),y(i),z(i),constr(i))
  !close(FREFCONFIG)
  
  ! Read the counter in order to continue the run where it stopped
  ! Format:
  ! Counter:     1000

  if (.not. restart) then 
     inquire( file = COUNTER, exist = flag )
     if ( flag .and. iproc == 0 ) then 
        open(unit=FCOUNTER,file=COUNTER,status='old',action='read',iostat=ierror)
        read(FCOUNTER,'(A12,I6)') dummy, mincounter 
        close(FCOUNTER)
     else
        mincounter = 1000
     end if
  end if 

  ! Read atomic file
  call init_all_atoms( nat_test, typa, posa, const_, boxref_, boundary, nproc, iproc )

  ! test nat_test and nat
  if ( nat_test /= NATOMS ) then
     if ( iproc == 0 ) write(*,*) "Different number of atoms"
     call end_art()
  end if
  !assign the data from the atomic file
  if ( .not. restart ) then
     typat(:)   = typa(:)
     pos(:)     = posa(:)
     constr(:)  = const_(:)
     boxref(:)  = boxref_(:)
     refcounter = mincounter
     box = boxref
  else if ( restart .and. dual_search) then
     call neighbours_local( )
  endif

  deallocate(posa)
  deallocate(typa)
  deallocate(const_)

  !atomic positions are now all read

  do i = 1, NATOMS
     Atom(i) = type_name(typat(i))
  end do
  !Atom(:) = type_name(typat(:))

  if ( iproc == 0 ) then 
   open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
   write(FLOG,'(1X,A34,I17)') ' - Mincounter                   : ', mincounter
   close(FLOG)
  end if

  ! We rescale the coordinates 
  scalaref = 1.0d0
  scala = scalaref

  call initialize_potential()         ! Initialize Potential 

  if ( iproc == 0 ) call convert_to_chain( refcounter, 4, scounter )
  fname = FINAL // scounter
  conf_initial = fname
                                      ! If this is a new event, 
  If_ne: if ( new_event .and. (.not. restart) ) then 
     call min_converge( success )     ! Converge the configuration to a local minimum

     posref = pos                     ! New reference configuration.
     mincounter = mincounter + 1
     ref_energy = total_energy
 
     if ( iproc == 0 ) then
                                      ! THIS iS NOT NECESSARY IN BIGDFT
        !call write_refconfig( )       ! Write reference in REFCONFIG. 
        call store( fname )           ! Store the configuration into fname.

        open( unit = FLOG, file = LOGFILE, status = 'unknown',&
        & action = 'write', position = 'append', iostat = ierror )
        write(*,*) 'BART: Configuration stored in file ',fname
        write(FLOG,'(1X,A34,A17)') ' - Configuration stored in file : ', trim(fname)
        if ( .not. success ) then
          write(FLOG,'(1X,A)') "ERROR: Initial configurations is not a minimum"
          call end_art()  
        end if
        close(FLOG)
     end if
                                      ! if dual_search we dont do this check at the
                                      ! beginning. It is not well defined.
     if ( LANCZOS_MIN .and. success .and. ( .not. dual_search ) ) call check_min( 'M' ) 
  else if ( (.not. new_event) .and. (.not. restart) ) then
     posref = pos
     ref_energy = total_energy
  end if If_ne

END SUBROUTINE initialize
