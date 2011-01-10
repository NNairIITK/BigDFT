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
  use lanczos_defs, only: projection, LANCZOS_MIN 
  implicit None

  integer :: ierror, ierr
  integer :: npart             ! Number of atoms participating to the eventt
  real(kind=8) :: delr
  real(kind=8) :: ran3, random_number
  real(kind=8) :: difpos
  real(kind=8) :: saddle_energy
  real(kind=8) :: delta_e  ! total_energy-ref_energy 
  real(kind=8), dimension(:), allocatable :: del_pos
  logical       :: success
  character(8)  :: accept 
  character(20) :: fname
  character(4)  :: scounter
  real(kind=8) :: a1, b1, c1
! _________
  call CPU_TIME( t1 )

  call read_parameters( )             ! Read options & initial allocation.
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
                                      ! If restartfile exists, then we restart.
  inquire ( file = restartfile, exist = restart )
  if ( restart ) & 
     & call restart_states( state_restart, ievent_restart, iter_restart, atp )
  call MPI_Barrier( MPI_COMM_WORLD, ierr )
  call write_parameters( )            ! Write options in LOGFILE. 
  call MPI_Barrier( MPI_COMM_WORLD, ierr )

  if ( iproc == 0 ) then              ! Report 
     open( unit = FLOG, file = LOGFILE, status = 'unknown',&
         & action = 'write', position = 'append', iostat = ierror )
     if ( restart ) then
        write(FLOG,'(1X,A)') ' Call restart'  
     else
        write(FLOG,'(1X,A)') ' Start with new event    '
     end if
     close(FLOG)
  end if

  call initialize_potential()         ! Initialize BigDFT

  if ( restart ) then        
     if ( iproc == 0 ) call convert_to_chain( refcounter, 4, scounter )
                             ! Information for FLIST 
     conf_initial = trim( FINAL // scounter )
  else 
                             ! Set up of mincounter & referecence configuration,
                             ! if new_event then relaxes it into a local minimum. 
     evalf_number = 0 
     call initialize( )                
     ievent_restart = 1
  end if
! _________
!                MAIN LOOP OVER THE EVENTS.

  Do_ev: do ievent = ievent_restart, NUMBER_EVENTS  
                                      ! If it is a restart event for the activation,
                                      ! phase 1, 2 or 4, or it is not a restart event
                                      ! then we call find_saddle.
     if ( .not. ( restart .and. ( state_restart == 3 ) ) ) then
                                      ! atp is the number of attempts for event.
        if (.not. restart) atp = 1 
        do 
          if ( iproc == 0 ) call print_event( ievent, temperature )
          call find_saddle( success, saddle_energy )
          if ( success ) then 
             exit
          else
             atp = atp + 1 
          end if 
        end do
     end if
                                      ! If not a new event, then a convergence to
                                      ! the saddle point, We do not go further.
     if ( .not. NEW_EVENT .and. ( eventtype == 'REFINE_SADDLE' ) ) call end_art () 
                                      ! If it is a restart event of type 3, we 
                                      ! are starting at the right point.
     if ( restart .and. ( state_restart == 3 ) ) then
        if ( iproc == 0 ) then        ! Report
           call print_event( ievent, temperature )
           open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
               & action = 'write', position = 'append', iostat = ierror )
           write(FLOG,'(1X,A)') ' - Restart event: We restart at the saddle point'
           close(FLOG)
           call convert_to_chain( mincounter, 4, scounter )
        end if
        saddle_energy = total_energy
        conf_saddle = trim( SADDLE // scounter )
        restart = .false.
     else if ( write_restart_file ) then 
                                      ! We save the current state for a possible 
                                      ! restart.
        state_restart = 3
        iter_restart  = 0
        if ( iproc == 0 ) call save_state( state_restart, iter_restart, projection )
     end if
                                      ! Push the configuration slightly over the
                                      ! saddle point in order to minimise  the odds 
                                      ! that it falls back into its original state.
     
     allocate(del_pos(VECSIZE))       ! We compute the displacement.
     del_pos =  pos - posref          
     a1 = dot_product(del_pos,projection) 
     if ( a1 < 0.0d0 ) then
        projection = -1.0d0 * projection
        a1 = -1.0d0 * a1
     end if
     b1 = dot_product(force,projection)
     c1 = dot_product(del_pos,force)
     difpos  = sqrt( dot_product(del_pos,del_pos) )
     !pos = pos + del_pos 
     deallocate(del_pos)
                                      
     pos = pos + PUSH_OVER * difpos * projection

     call min_converge( success )             ! And we converge to the new minimum.
     delta_e = total_energy - ref_energy 
     if ( iproc == 0 ) then
                                      ! We write the configuration in a min.... file.
        call convert_to_chain( mincounter, 4, scounter )
        write(*,*) 'BART: Mincounter is ', mincounter,', scounter is ', scounter
        fname = FINAL // scounter
        conf_final = fname
        call store( fname ) 
     end if
                                      ! Is a real minimum ?
     if ( LANCZOS_MIN .and. success ) call check_min( 'M' ) 
                                      ! Magnitude of the displacement (utils.f90).
     call displacement( posref, pos, delr, npart )
                                      ! Now, we accept or reject this move based
                                      ! on a Boltzmann weight.
     if ( iproc == 0 ) then           
        random_number = ran3( ) 
        open( unit = FLIST, file = EVENTSLIST, status = 'unknown', &
       & action = 'write', position = 'append', iostat = ierror )
     end if

     If_bol: if ( ( (total_energy - ref_energy) < -temperature * log( random_number ) )& 
                & .and. ( temperature >= 0.0d0 ) .and. success ) then
        accept = "ACCEPTED"
        if ( iproc == 0 )&
        &  write(FLIST,*) conf_initial, conf_saddle, conf_final,'    accepted'
                                      ! We now redefine the reference configuration
        scalaref     = scala
        posref       = pos     
        conf_initial = conf_final
        ref_energy   = total_energy
        refcounter   = mincounter
                                      ! THIS IS NOT NECESSARY IN BIGDFT
                                      ! Updating the reference configuration file,
                                      ! which serves as the initial configuration
                                      ! for events.
        !if ( iproc == 0 ) call write_refconfig( )
     else                             ! Else If_bol:
                                      ! If the event is not accepted we start
                                      ! from the previous refconfig.
        accept = "REJECTED"
        if ( iproc == 0 ) then        ! Write
           if (( total_energy - ref_energy ) > 1.0d-5 )  then
               write(FLIST,*) conf_initial, conf_saddle, conf_final,'    rejected'
           else  
               write(FLIST,*) conf_initial, conf_saddle, conf_final,'    exchanged'
           end if
        end if

     end if If_bol

     if ( iproc == 0 ) then           ! Report
        close(FLIST)
        open( unit = FLOG, file = LOGFILE, status = 'unknown',&
            & action = 'write', position = 'append', iostat = ierror ) 
        write(*,*) 'BART: Configuration stored in file ',fname
        write(FLOG,'(1X,A34,A17)') ' - Configuration stored in file : ', trim(fname)
        write(FLOG,'(1X,A34,(1p,e17.10,0p))')&
        &  ' - Total energy Minimum (eV)    : ', total_energy 
        write(FLOG,"(' ','MINIMUM',i5, a9,' |E(fin-ini)= ', f9.4,' |E(fin-sad)= ',&
          & f9.4,' |npart= ', i4,' |delr= ', f8.3,' |evalf=', i6,' |')")& 
          & mincounter, adjustr(accept), delta_e,                                 &
          & total_energy - saddle_energy, npart, delr, evalf_number
        write(*,"(' ','BART: MINIMUM',i5, a9,' |E(fin-ini)= ', f9.4,' |E(fin-sad)= ',&
          & f9.4,' |npart= ', i4,' |delr= ', f8.3,' |evalf=', i6,' |',f8.3,3f7.2)")       & 
          & mincounter, adjustr(accept), delta_e,                              &
          & total_energy - saddle_energy, npart, delr, evalf_number, difpos,   &
          & a1, b1, c1 

        ! Debug 
        !write(FLOG,'(1X,A34,F17.6)')&          
        !& ' - -T*log( random_number )      : ', -temperature*log( random_number )
        close(FLOG)
     end if

     if ( eventtype == "REFINE_AND_RELAX" ) call end_art() 

     mincounter = mincounter + 1

     if ( iproc == 0 ) then
        open(unit=FCOUNTER,file=COUNTER,status='unknown',action='write',iostat=ierror)
        write(FCOUNTER,'(A12,I6)') 'Counter:    ', mincounter
        close(FCOUNTER)
     end if

  end do Do_ev

  call end_art( )

END PROGRAM art90
!!***


!!****f* art/end_art( )
!! FUNCTION
!! SOURCE
!!
subroutine end_art( )

  use defs

  implicit none

  integer      :: ierror, ierr
  real(kind=8) :: t2  ! cputime

  call finalise_potential( )

  if ( iproc == 0 ) then 

    open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
        & action = 'write', position = 'append', iostat = ierror )
    write(FLOG,*) '********************** '
    write(FLOG,*)    '  A bientot !'
    write(FLOG,*) '********************** '
    call CPU_TIME( t2)
    write(FLOG,"(' CPU_TIME: ', f12.4, ' seg')") t2-t1
    call timestamp('End') 
    close(FLOG)
  end if 

  if (nproc > 1) call MPI_FINALIZE(ierr)
  stop

END SUBROUTINE end_art
