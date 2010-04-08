!!****f* art/initialize
!! FUNCTION
!!   Initialization of art method
!! 
!! DESCRIPTION
!!   then relaxes it into a  a local minimum with or without volume optimization
!!   depending on the compilation flags. 
!! 
!!   The initial configuration is a "reference configuration", it will be the reference 
!!   configuration until a new event is accepted. 
!! 
!!   The reference configuration should contain the following information:
!! 
!!     first line:         refcounter 1232      
!!     second line:        total energy: 
!!     third line:          -867.33454
!!     fourth line:        8.34543 8.21345 8.67789
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
subroutine initialize()

  use defs
  use lanczos_defs
  implicit none

  !Local variables
  integer :: i, ierror
  character(len=20) :: dummy,fname
  character(len=4) ::scounter
  logical :: flag,  new_projection
  
  ! Read the counter in order to continue the run where it stopped
  ! Format:
  ! Counter:     1000

  inquire( file = COUNTER, exist = flag )
  if (flag) then 
     open(unit=FCOUNTER,file=COUNTER,status='old',action='read',iostat=ierror)
     read(FCOUNTER,'(A12,I6)') dummy, mincounter 
     close(FCOUNTER)
  else
     mincounter = 1000
  end if

  do i = 1, NATOMS
     Atom(i) = type_name(type(i))
  end do

  if ( iproc == 0 ) then 
   open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
   write(FLOG,'(1X,A34,I17)') ' - Mincounter                   : ', mincounter
   close(FLOG)
  end if

  ! We rescale the coordinates and define the reference coordinates
  scalaref = 1.0d0
  scala = scalaref

  box = boxref

  xref = x ! Vectorial operation
  yref = y ! Vectorial operation
  zref = z ! Vectorial operation

  total_energy = ref_energy
                                      ! If this is a new event, 
                                      ! we minimize before starting.
  If_ne: if ( new_event ) then 
                                      ! Converge the configuration to a local minimum
     call min_converge( )
                                      ! Assign the relaxed configuration as the
                                      ! reference and write it to ref-file.
     posref = pos
     if ( iproc == 0 ) call write_refconfig()
 
     if ( iproc == 0 ) call convert_to_chain( mincounter, scounter )
     fname =   FINAL // scounter
     conf_initial = fname
                                      ! We now store the configuration into fname.
     if ( iproc == 0 ) then
        call store( fname )
                                      ! Write
        open( unit = FLOG, file = LOGFILE, status = 'unknown',&
        & action = 'write', position = 'append', iostat = ierror )
        write(*,*) 'BART: Configuration stored in file ',fname
        write(FLOG,'(1X,A34,A17)') ' - Configuration stored in file : ', trim(fname)
        close(FLOG)
     end if

     mincounter = mincounter + 1
     ref_energy = total_energy

  end if If_ne

END SUBROUTINE initialize
!!***

