!!****p* BigDFT/art90
!! FUNCTION
!!   Main program to use BigDFT with art nouveau method
!!
!! COPYRIGHT
!!    Copyright (C) Normand Mousseau, June 2001
!!    Copyright (C) 2010 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! SOURCE
!!
program art90

  use defs 
  use random 
  use lanczos_defs, Only: projection
  implicit None

  integer :: ierror, ierr
  integer :: npart             ! Number of atoms participating to the eventt
  real(8) :: delr
  real(8) :: ran3, random_number
  real(8) :: difpos
  real(8) :: saddle_energy
  real(8), dimension(:), allocatable :: del_pos

  logical       :: success
  character(8)  :: date
  character(10) :: time
  character(5)  :: zone
  character(20) :: fname
  character(4)  :: scounter
  integer, dimension(8) :: value

! _________
                                      ! Read options.
  call read_parameters( )

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

                                      ! Write options in LOGFILE 
  call write_parameters( )

                                      ! Initialize BigDFT
  call initialize_potential() 

  if ( iproc == 0 ) then              ! Work only on the master node. 

                                      ! Initialization of seed in random.
     call date_and_time(date,time,zone,value)
     idum = -1 * mod( (1000 * value(7) + value(8)), 1024) 

                                      ! We decide whether or not we restart,
                                      ! if the restart exists, then we restart.
     inquire ( file = restartfile, exist = restart )

                                      ! Write
     open( unit = FLOG, file = LOGFILE, status = 'unknown',&
         & action = 'write', position = 'append', iostat = ierror )
     write(FLOG,'(1X,A34,I17)') ' - The seed is                  : ', idum
     write(FLOG,'(1X,A34,L17)') ' - Restart                      : ', restart
     close(FLOG)

  end if

  if ( restart ) then                 ! Restarts the initial configuration.

                                      ! Write
     if ( iproc == 0 ) then 
      open( unit = FLOG, file = LOGFILE, status = 'unknown',&
          & action = 'write', position = 'append', iostat = ierror )
      write(FLOG,'(1X,A)') ' Call restart'   
      close(FLOG)
     end if 

     call restart_states( state_restart, ievent_restart, iter_restart )

  else 
                                      ! Write
     if ( iproc == 0 ) then
      open( unit = FLOG, file = LOGFILE, status = 'unknown',&
          & action = 'write', position = 'append', iostat = ierror )
      write(FLOG,'(1X,A)') ' Start with new event    '
      close(FLOG)
     end if 

     call initialize()                ! Read in the initial configuration
     ievent_restart = 1
    
     if ( iproc == 0 ) call write_refconfig( ) 

  end if

                                      ! We define the name of the initial file.
  mincounter = mincounter - 1
  if ( iproc == 0 ) call convert_to_chain( mincounter, scounter )
  fname = FINAL // scounter
  conf_initial = fname
  mincounter = mincounter + 1

                                      ! Main loop over the events.
  Do_ev: do ievent = ievent_restart, NUMBER_EVENTS  

     if ( iproc == 0 ) call print_newevent( ievent, temperature )

                                      ! If it is a restart event for the activation,
                                      ! phase 1 or 2, or it is not a restart event
                                      ! then we call find_saddle.
     if ( .not. ( restart .and. ( state_restart == 3 ) ) ) then
        do 
          call find_saddle( success, saddle_energy )
          if ( success ) exit
        end Do
     end if
                                      ! If not a new event, then a convergence to
                                      ! the saddle point, We do not go further.
     if ( .not. NEW_EVENT .and. ( eventtype == 'REFINE_SADDLE' ) ) stop
     
                                      ! If it is a restart event of type 3, we 
                                      ! are starting at the right point.
     if ( restart .and. ( state_restart == 3 ) ) then

                                      ! Write
        if ( iproc == 0 ) then
         open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
             & action = 'write', position = 'append', iostat = ierror )
         write(FLOG,'(1X,A)') ' - Restart event'
         write(FLOG,'(1X,A)') ' - We restart at the saddle point '
         close(FLOG)
        end if

        restart = .false.

     else                                  
                                      ! We save the current state for a possible 
                                      ! restart.
        state_restart = 3
        iter_restart  = 0
        if ( iproc == 0 ) call save_state( state_restart, iter_restart, projection )

     end if

                                      ! Push the configuration slightly over the
                                      ! saddle point in order to minimise  the odds 
                                      ! that it falls back into its original state.
     
                                      ! We first compute the displacement.
     allocate(del_pos(VECSIZE))
     del_pos =  pos - posref
     difpos  = sqrt( dot_product(del_pos,del_pos) )
     deallocate(del_pos)

                                      ! Projection points away from the initial 
                                      ! minimum.
     pos = pos + PUSH_OVER * difpos * projection

                                      ! And we converge to the new minimum.
     call min_converge( )

                                      ! We write the configuration in a min.... file
     if ( iproc == 0 ) call convert_to_chain( mincounter, scounter )
     if ( iproc == 0 ) write(*,*) 'BART: Mincounter is : ', mincounter,&
                       & ' and scounter is ', scounter
     fname = FINAL // scounter

                                      ! Magnitude of the displacement (utils.f90).
     call displacement( posref, pos, delr, npart )
     
     if ( iproc == 0 ) then
                                      ! We now store the configuration into fname
        Call store( fname ) 
                                      ! Open files.
        open( unit = FLOG, file = LOGFILE, status = 'unknown',&
            & action = 'write', position = 'append', iostat = ierror ) 
        open( unit = FLIST, file = EVENTSLIST, status = 'unknown',&
            & action = 'write', position = 'append', iostat = ierror )
                                      ! Write
        write(*,*) 'BART: Configuration stored in file ',fname
        write(FLOG,'(1X,A34,A17)') ' - Configuration stored in file : ', trim(fname)
     end if

     conf_final = fname
     mincounter = mincounter+1

                                      ! Write 
     if ( iproc == 0 )  then
      write(FLOG,'(1X,A34,(1p,e17.10,0p))')&
      &  ' - Total energy Minimum (eV)    : ', total_energy 
      write(FLOG,'(1X,A34,F17.6)')&
      &  ' - E( fin-ini )                 : ', total_energy - ref_energy
      write(FLOG,'(1X,A34,F17.6)')& 
      &  ' - E( fin-sad )                 : ', total_energy - saddle_energy  
      write(FLOG,'(1X,A34,I17)')&
      &  ' - npart                        : ', npart
      write(FLOG,'(1X,A34,F17.4)')&
      &  ' - delr( fin-ini )              : ', delr
      write(FLOG,'(1X,A34,I17)')&
      &  ' - evalf_number                 : ', evalf_number
      write(*,"(' ','BART: Total energy M: ',(1p,e17.10,0p),'  npart: ', i4,'  delr: ',&
      & f12.6,' Number evaluations: ',i6)") &
      & total_energy, npart, delr, evalf_number
     end if

                                      ! Now, we accept or reject this move based
                                      ! on a Boltzmann weight.
     if ( iproc == 0 ) random_number = ran3()

     If_bol: if ( ( (total_energy - ref_energy) < -temperature * log( random_number ) )& 
                & .and. ( temperature >= 0.0d0 ) ) then

                                      ! Write 
        if ( iproc == 0 ) Then        
         write(FLIST,*) conf_initial, conf_saddle, conf_final,'    accepted'
         close(FLIST)

         write(*,*) 'BART: New configuration accepted, mincounter was: ', mincounter-1
         write(FLOG,*) ' New configuration ACCEPTED'
         write(FLOG,'(1X,A34,F17.6)')&
         & ' - -T*log( random_number )      : ', -temperature*log( random_number ) 
         write(FLOG,'(1X,A34,I17)') ' - mincounter was               : ', mincounter-1
        end if

                                      ! We now redefine the reference configuration
        scalaref     = scala
        posref       = pos     
        conf_initial = conf_final
        ref_energy   = total_energy

        if ( eventtype == "REFINE_AND_RELAX" ) then
           close (FLIST) 
           close (FLOG)
           stop
        end if 
                                      ! Update the reference configuration file,
                                      ! which serves as the initial configuration
                                      ! for events.
        if ( iproc == 0 ) call write_refconfig( ) 

     else                             ! Else If_bol:
                                      ! The event is not accepted; we start
                                      ! from the previous refconfig.

                                      ! Write
        if ( iproc == 0 ) then
         if (( total_energy - ref_energy ) > 1.0d-5 )  then
             write(FLIST,*) conf_initial, conf_saddle, conf_final,'    rejected'
         else  
             write(FLIST,*) conf_initial, conf_saddle, conf_final,'    exchanged'
         end if

         write(*,*) 'BART: New configuration rejected, mincounter was : ', mincounter-1
         write(FLOG,*) ' New configuration REJECTED'
         write(FLOG,'(1X,A34,F17.6)')&
         & ' - log( random_number )         : ', log( random_number )
         write(FLOG,'(1X,A34,I17)') ' - mincounter was               : ', mincounter-1
        end if

     end if If_bol

     if ( iproc == 0 ) then

        open(unit=FCOUNTER,file=COUNTER,status='unknown',action='write',iostat=ierror)
        write(FCOUNTER,'(A12,I6)') 'Counter:    ', mincounter
        close(FCOUNTER)
        close(FLOG)
        close(FLIST)
       
     end if

  end do Do_ev

  call finalise_potential( )

end program art90
!!***


!!****f* art/print_newevent
!! FUNCTION
!!   This subroutine prints the initial details for a new events
!! SOURCE
!!
subroutine print_newevent( ievent_current, temperat )

  use defs
  implicit none

  !Arguments
  integer, intent(in) :: ievent_current
  real(8), intent(in) :: temperat

  !Local variables
  integer :: ierror

  write(*,*) 'BART: Simulation : ', ievent_current
  write(*,*) 'BART: Starting from minconf : ', mincounter
  write(*,*) 'BART: Reference Energy (eV) : ', ref_energy
  write(*,*) 'BART: Temperature : ', temperat

  open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
  write(FlOG,*) ''
  write(FLOG,*) ' _______________________________________'
  write(FLOG,'(1X,A34,I17)') ' - Simulation                   : ', ievent_current
  write(FLOG,'(1X,A34,I17)') ' - Starting from minconf        : ', mincounter
  write(FLOG,'(1X,A34,(1p,e17.10,0p))') ' - Reference Energy (eV)        : ', ref_energy 
  write(FLOG,'(1X,A34,F17.6)') ' - Temperature                  : ', temperat
  close(FLOG)

END SUBROUTINE print_newevent
!!***

