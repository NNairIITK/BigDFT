!>  Main program to use BigDFT with art nouveau method
!!
!! @author
!!    Copyright Normand Mousseau, July 2001
!!    Copyright (C) 2010 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
program art90

  use DEFS
  use RANDOM
  use lanczos_defs

  implicit none

  integer :: ierror
  integer :: npart             ! Number of atoms participating to the eventt
  real(8) :: del_r
  real(8) :: ran3, random_number
  real(8) :: difpos
  real(8), dimension(:), allocatable :: del_pos

  logical :: success
  character(8) :: date
  character(10) :: time
  character(5) :: zone
  character(len=20) :: fname
  character(4) :: scounter
  integer, dimension(8) :: value

  call read_parameters()
  call initialize_potential()

  allocate(del_pos(VECSIZE))

  ! Only do for the main process
  if (iproc .eq. 0 ) then
     call date_and_time(date,time,zone,value)
     idum = -1 * mod( (1000 * value(7) + value(8)), 1024) 

     open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
     write(flog,*) 'The seed is : ', idum
     del_r = ran3()
  endif

  ! We now decide whether or not we restart if the restart exists, then we restart
  inquire(file=restartfile,exist=restart)
  write(flog,*) 'Restart : ', restart
  close(unit=flog)

  if (restart) then          ! Restarts the initial configuration
     open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
     write(flog,*) 'call restart'   
     call restart_states(state_restart,ievent_restart,iter_restart) 
     close(unit=flog)
  else 
     open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
     write(flog,*) 'Start with new event'
     close(unit=flog)
     call initialize()       ! Read in the initial configuration
     ievent_restart = 1
     
     if (iproc .eq. 0 ) call write_refconfig() 
  endif

  ! We define the name of the initial file
  mincounter = mincounter -1
  call convert_to_chain(mincounter,scounter)
  fname = FINAL // scounter
  conf_initial = fname
  mincounter = mincounter +1

  do ievent= ievent_restart, NUMBER_EVENTS         ! Main loop over the events

    if (iproc .eq. 0 )  call print_newevent(ievent,temperature)

    ! We look for a local saddle point b

    ! If it is a restart event for the activation, phase 1 or 2, or it is not a restart event
    ! then we call find_saddle
    if ( .not. (restart .and. (state_restart .eq. 3))  ) then
       do 
          call find_saddle( success )
          if ( success ) exit
       end do
    endif

    ! If not a new event, then a convergence to the saddle point, We do not go further
    if( .not.NEW_EVENT .and. (eventtype .eq. 'REFINE_SADDLE') ) stop

    ! If it is a restart event of type 3, we are starting at the right point
    if ( restart .and. (state_restart .eq. 3) ) then
       if (iproc .eq. 0 )  then
          open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
          write(flog,*) 'Restart event'
          write(flog,*) 'We restart at the saddle point '
          close(unit=flog)
       endif
       restart = .false.
    else ! We save the current state for a possible restart
       state_restart = 3
       iter_restart = 0
        if (iproc .eq. 0 ) call save_state(state_restart,iter_restart,projection)
    endif

    ! Push the configuration slightly over the saddle point in order to minimise 
    ! the odds that it falls back into its original state

    ! we first compute the displacement

    del_pos(:) =  pos(:) - posref(:)
    difpos = sqrt( dot_product(del_pos,del_pos) )

    ! Projection points away from the initial minimum
    pos(:) = pos(:) + PUSH_OVER * difpos * projection(:)

    ! And we converge to the new minimum.
    call min_converge()

    ! We need to write the configuration in a min.... file
    call convert_to_chain(mincounter,scounter)
     if (iproc .eq. 0 ) write(*,*) ' Mincounter is : ', mincounter, ' and scounter is ', scounter
    fname =   FINAL // scounter

    ! Compute the displacement and the number of atoms involved in the event     
    call displacement(pos,posref,del_r,npart)
    
    ! We now store the configuration into fname
    if (iproc .eq. 0 ) call store(fname)

    if (iproc .eq. 0 )  then
       open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
       write(*,*) 'Configuration stored in file ',fname
       write(FLOG,*) 'Configuration stored in file ',fname
    endif

    conf_final = fname
    mincounter = mincounter+1

    ! We write out various information to both screen and file
    if (iproc .eq. 0 )  then
       write(*,*) 'Total energy M: ', total_energy, '  npart: ',npart, &
            'del_r: ',del_r,' Number evaluations: ', evalf_number

       write(FLOG,*) 'Total energy M: ', total_energy, '  npart: ',npart, &
            'del_r: ',del_r,' Number evaluations: ', evalf_number

       ! Now, we accept or reject this move based on a Boltzmann weight

       open(unit=FLIST,file=EVENTSLIST,status='unknown',action='write',position='append',iostat=ierror)
    endif

    if (iproc .eq. 0 ) random_number = ran3()
    call MPI_Bcast(random_number,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)

    if( ( (total_energy - ref_energy) < -temperature * log(ran3()) ) .and. (temperature .ge. 0.0d0) ) then

       if (iproc .eq. 0 ) then        
          write(*,*) 'New configuration accepted, mincounter was : ', mincounter-1
          write(FLOG,*) 'New configuration accepted, mincounter was : ', mincounter-1
          write(FLIST,*) conf_initial, conf_saddle, conf_final,'    accepted'
       endif

      ! We now redefine the reference configuration
      scalaref = scala
      posref= pos          ! This is a vectorial copy
      conf_initial = conf_final

      ref_energy = total_energy

      if (eventtype .eq. "REFINE_AND_RELAX") then
         stop
      endif

      ! Update the reference configuration file, which serves as the initial
      ! configuration for events.
       if (iproc .eq. 0 ) call write_refconfig() 
    else

       if (iproc .eq. 0 ) then
        write(*,*) 'New configuration rejected, mincounter was : ', mincounter-1
        write(FLOG,*) 'New configuration rejected, mincounter was : ', mincounter-1

        ! The events is not accepted; we start from the previous refconfig
        if((total_energy - ref_energy) > 1.0d-5)  then
           write(FLIST,*) conf_initial, conf_saddle, conf_final,'    rejected'
        else  
           write(FLIST,*) conf_initial, conf_saddle, conf_final,'    exchanged'
        endif
     endif
    endif

    if (iproc .eq. 0 )  then
       close(FLIST)
       close(FLOG)

       open(unit=FCOUNTER,file=COUNTER,status='unknown',action='write',iostat=ierror)
       write(FCOUNTER,'(A12,I6)') 'Counter:    ', mincounter
       close(FCOUNTER)
    endif
  end do

  deallocate(del_pos)

  call finalise_potential()
end program art90



!>   This subroutine prints the initial details for a new events
!!
!!
subroutine print_newevent(ievent_current,temperat)
  use defs
  implicit none
  integer, intent(in) :: ievent_current
  real(8), intent(in) :: temperat
  integer :: ierror;

  write(*,*) 'Simulation : ', ievent_current
  write(*,*) 'Starting from minconf : ', mincounter
  write(*,*) 'Temperature : ', temperat

  open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
  write(FLOG,*) 'Simulation : ', ievent_current
  write(FLOG,*) 'Starting from minconf : ', mincounter
  write(FLOG,*) 'Initial energy : ', total_energy 
  write(FLOG,*) 'Temperature : ', temperat
  close(FLOG)

  return
END SUBROUTINE print_newevent

