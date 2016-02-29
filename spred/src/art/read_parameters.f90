!> @file
!! Read the parameters for BigDFT+ART
!! @author
!!    Copyright (C) 2001 Normand Mousseau
!!    Copyright (C) 2010 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> ART read_parameters
!! Read the parameters defining the simulation 
subroutine read_parameters( )

  use defs
  use lanczos_defs
  use saddles
  implicit none  

  !Local variables
  Character(len=40)  :: temporary 

!SECTION____________________________ ATOMS  ( Obsolete ) 
 
  call get_environment_variable('NATOMS', temporary)
  if (temporary .eq. '') then
     write(*,*) 'Error: NATOMS is not defined'
     stop
  else
     read(temporary,*) natoms
  end if

  !!__________________
  ! We now get the types - define up to 5
  call get_environment_variable('type1',temporary)
  if (temporary .eq. '') then
     write(*,*) 'Error, must at least define 1 type of atoms -  use "setenv type1 Si", for example'
     stop
  else
     read(temporary,*) type_name(1)
  end if

  call get_environment_variable('type2',temporary)
  if (temporary .eq. '') then
     type_name(2) = ''
  else
     read(temporary,*) type_name(2)
  end if

  call get_environment_variable('type3',temporary)
  if (temporary .eq. '') then
     type_name(3) = ''
  else
     read(temporary,*) type_name(3)
  end if

  call get_environment_variable('type4',temporary)
  if (temporary .eq. '') then
     type_name(4) = ''
  else
     read(temporary,*) type_name(4)
  end if

  call get_environment_variable('type5',temporary)
  if (temporary .eq. '') then
     type_name(5) = ''
  else
     read(temporary,*) type_name(5)
  end if

  call get_environment_variable('type6',temporary)
  if (temporary .ne. '') then
     write(*,*) 'Error: The code can only handle 5 atomic types, change read_parameters.f90, to allow for more.'
     stop
  end if

!SECTION_____________________________ ART 

  ! Hide Option, by default .false.
  ! Only Lanczos analysis for the minimum and the inflection point.
  call get_environment_variable('Setup_Initial',temporary)
  if (temporary == '') then
     setup_initial = .false. 
  else
     read(temporary,*)  setup_initial
  end if

  !!__________________
  call get_environment_variable('EVENT_TYPE', temporary)
  selectcase( temporary )
     case ( 'NEW' )
          new_event = .true.
          eventtype = 'NEW'
     case ( 'REFINE_SADDLE' )
          new_event = .false.
          eventtype = 'REFINE_SADDLE'
     case ( 'REFINE_AND_RELAX' )
          new_event = .false.
          eventtype = 'REFINE_AND_RELAX'
     case ( 'GUESS_DIRECTION' )
          new_event = .true.
          eventtype = 'GUESS_DIRECTION'
     case default 
          write(*,*) 'Error: event_types permitted are:'
          write(*,*) 'NEW, REFINE_SADDLE, REFINE_AND_RELAX and GUESS_DIRECTION'
          stop
  end select

  !!__________________
  if ( eventtype == 'GUESS_DIRECTION' ) then

     !!__________________HIDE 
     call get_environment_variable('Noise',temporary)
     if (temporary == '') then
        guess_noise = 0.0d0
     else
        read(temporary,*) guess_noise 
     end if
     !!__________________NAME OF FILE
     call get_environment_variable('INITDIR', temporary)
     if (temporary .eq. '') then
        GUESSFILE = 'initdir'
     else
        read(temporary,*) GUESSFILE 
     end if
     GUESSFILE = adjustl(GUESSFILE)

  end if

  !!__________________
  ! type of energy and force calculation 
  call get_environment_variable('ENERGY_CALC', temporary)
  if (temporary .eq. '') then
     write(*,*) "Error: energy calculation type is not defined: ENERGY_CALC " 
     write(*,*) " choose: BIG or SWP, BSW, OTF (buggy), BAY (buggy)  "
     stop
  else
     read(temporary,*) energy_type   
  endif

  !!__________________
  ! Fictive temperature, if negative always reject the event
  call get_environment_variable('Temperature', temporary)
  if (temporary .eq. '') then
     write(*,*) 'Error: Metropolis temperature is not defined'
     stop
  else
     read(temporary,*) temperature
  end if

  !!__________________
  ! Maximum number of events
  call get_environment_variable('Max_Number_Events', temporary)
  if (temporary .eq. '') then
     number_events = 100
  else
     read(temporary,*) number_events
  end if

  !!__________________
  ! Read type of events 
  ! Activation: global, local, list_local, list, local_coord 
  call get_environment_variable('Type_of_Events', TYPE_EVENTS)
  if ( (TYPE_EVENTS .ne. 'global') .and. (TYPE_EVENTS .ne. 'local') .and. &
     & (TYPE_EVENTS .ne. 'list')   .and. (TYPE_EVENTS .ne. 'list_local') .and. &
     & (TYPE_EVENTS .ne. 'local_coord') ) then
     write(*,*) 'Error : only global, local, or list type of events are accepted - provided: ',&
     & TYPE_EVENTS
     stop
  end if

  !!__________________
  ! Cutoff for local and list_local (in angstroems)
  if (TYPE_EVENTS .eq. 'local' .or. TYPE_EVENTS .eq. 'list_local' ) then
     call get_environment_variable('Radius_Initial_Deformation', temporary)
     if (temporary .eq. '') then
        write(*,*) 'Error: Radius_Initial_Deformation must be defined when TYPE_EVENTS is local'
        stop
     else
        read(temporary,*) LOCAL_CUTOFF
     end if
   end if

  !!__________________
  ! Number of the atom around which the initial move takes place
  if (TYPE_EVENTS .eq. 'local' ) then 
     call get_environment_variable('Central_Atom',temporary)
     if (temporary .eq. '') then
        preferred_atom = -1
     else
        read(temporary,*) preferred_atom
     end if
  end if

  !!__________________
  ! Breaks the symmetry of the crystal by randomly displacing
  ! all atoms by this distance
  call get_environment_variable('sym_break_dist',temporary)
  if (temporary .eq. '') then
     sym_break_dist = 0.0
  else
     read(temporary,*) sym_break_dist
  end if

  !!__________________
  ! Maximum number of iteraction for reaching the saddle point
  call get_environment_variable('Activation_MaxIter',temporary)
  if (temporary .eq. '') then
     MAXPAS = 100 
  else
     read(temporary,*) MAXPAS
  end if

  !!__________________
  ! Info regarding initial displacement
  call get_environment_variable('Initial_Step_Size', temporary)
  if (temporary .eq. '') then
     INITSTEPSIZE = 0.001             ! Size of initial displacement in Ang.
  else
     read(temporary,*) INITSTEPSIZE
  end if

  !!__________________
  ! Increment size - overall scaling (in angstroems)
  call get_environment_variable('Increment_Size',temporary)
  if (temporary .eq. '') then
     INCREMENT = 0.01
  else
     read(temporary,*) INCREMENT
  end if

  !!__________________
  ! Force threshold for the perpendicular relaxation
  call get_environment_variable('Force_Threshold_Perp_Rel',temporary)
  if (temporary .eq. '') then
     FTHRESHOLD = 1.0
  else
     read(temporary,*) FTHRESHOLD
  end if

!SECTION_____________________________ HARMONIC WELL

  ! Size of the parallel displacement in for leavign the basin with respect to
  ! the overall increment
  call get_environment_variable('Basin_Factor',temporary)
  if (temporary .eq. '') then
     BASIN_FACTOR = 1.0d0
  else
     read(temporary,*) BASIN_FACTOR
  end if

  !!__________________
  ! Maximum number of perpendicular steps leaving basin
  call get_environment_variable('Max_Perp_Moves_Basin',temporary)
  if (temporary .eq. '') then
     MAXKPERP = 2
  else
     read(temporary,*) MAXKPERP
  end if

  !!__________________
  ! Minimum number of steps in kter-loop before we call lanczos
  call get_environment_variable('Min_Number_KSteps',temporary)
  if (temporary .eq. '') then
     KTER_MIN = 1
  else
     read(temporary,*) KTER_MIN
  end if

  !!__________________
  ! Eigenvalue threshold for leaving basin
  call get_environment_variable('Eigenvalue_Threshold',temporary)
  if (temporary .eq. '') then
     write(*,*) 'Error : No eigenvalue threshold provided  (Eigenvalue_Threshold)'
     stop
  else
     read(temporary,*) EIGEN_THRESH
  end if

  !!__________________
  ! Maximum number of iteration in the activation
  call get_environment_variable('Max_Iter_Basin',temporary)
  if (temporary .eq. '') then
     MAXKTER = 100
  else
     read(temporary,*) MAXKTER
  end if

!SECTION_____________________________ LANCZOS

  ! Convergence criterion for the wavefunction optimization in Lanczos
  ! If not defined, we use that one in input.dft 
  call get_environment_variable('gnrm',temporary)
  if (temporary == '') then
     my_gnrm = 1.0d0
  else
     read(temporary,*) my_gnrm 
  end if

  !!_________________HIDE
  ! inflection in the eigenvalue
  ! if it is true; it calculates the projection only at every two steps 
  ! but after 4 steps above of an inflection in the eigenvalue, and if the
  ! last a1 >0.9d0
  call get_environment_variable('calc_of_projection',temporary)
  if (temporary == '') then
     calc_proj = .false. 
  else
     read(temporary,*) calc_proj 
  end if

  !!__________________
  ! Calculation of the Hessian for each minimum
  call get_environment_variable('Lanczos_of_minimum',temporary)
  if (temporary .eq. '') then
     LANCZOS_MIN = .False. 
  else
     read(temporary,*)  LANCZOS_MIN
  end if

  !!__________________
  ! Number of iterations in the lanczos Self consistent loop
  call get_environment_variable('Lanczos_SCLoop',temporary)
  if (temporary .eq. '') then
     LANCZOS_SCL = 1 
  else
     read(temporary,*) LANCZOS_SCL
  end if

  !!__________________
  ! The convergence criteria in the Lanczos Self consistent loop
  call get_environment_variable('Lanczos_collinear',temporary)
  if (temporary .eq. '') then
     collinear_factor= 0.7d0
  else
     read(temporary,*) collinear_factor 
  end if

  !!__________________
  ! Number of vectors included in lanczos procedure in the Harmonic well
  call get_environment_variable('Number_Lanczos_Vectors_H',temporary)
  if (temporary .eq. '') then
     NVECTOR_LANCZOS_H = 16
  else
     read(temporary,*) NVECTOR_LANCZOS_H
  end if

  !!__________________
  ! Number of vectors included in lanczos procedure in convergence
  call get_environment_variable('Number_Lanczos_Vectors_C',temporary)
  if (temporary .eq. '') then
     NVECTOR_LANCZOS_C = NVECTOR_LANCZOS_H  
  else
     read(temporary,*) NVECTOR_LANCZOS_C
  end if

  !!__________________
  ! The step of the numerical derivative of forces for the Hessian (in Ang) 
  call get_environment_variable('delta_disp_Lanczos',temporary)
  if (temporary .eq. '') then
     DEL_LANCZOS =  0.01
  else
     read(temporary,*) DEL_LANCZOS 
  end if

!SECTION____________________________ CONVERGENCE

  ! Maximum number of perpendicular steps during activation
  call get_environment_variable('Max_Perp_Moves_Activ',temporary)
  if (temporary .eq. '') then
     MAXIPERP = 12
  else
     read(temporary,*) MAXIPERP
  end if

  !!__________________
  ! Force threshold for the convergence at the saddle point
  call get_environment_variable('Exit_Force_Threshold',temporary)
  if (temporary .eq. '') then
     EXITTHRESH = 0.1 
  else
     read(temporary,*) EXITTHRESH
  end if

  !!__________________
  ! The prefactor for pushing over the saddle point, fraction of distance from
  ! initial minimum to saddle point
  call get_environment_variable('Prefactor_Push_Over_Saddle',temporary)
  if (temporary .eq. '') then
     PUSH_OVER = 0.15
  else
     read(temporary,*) PUSH_OVER
  end if

  !!__________________
  ! if delta_e < delta_thr .and. delr < delr_thr => end_activation = .true. 
  ! Set up them to zero if you dont want to use these criteria.
  call get_environment_variable('delta_threshold',temporary)
  if (temporary == '') then
     delta_thr = 0.0d0
  else
     read(temporary,*) delta_thr  
  end if

  call get_environment_variable('delr_threshold',temporary)
  if (temporary == '') then
     delr_thr = 0.0d0
  else
     read(temporary,*) delr_thr  
  end if

!SECTION_____________________________ DIIS

  ! Do we use DIIS for refining and converging to saddle
  call get_environment_variable('Use_DIIS',temporary)
  if (temporary .eq. '') then
     USE_DIIS = .false.
  else
     read(temporary,*) USE_DIIS
  end if

  !!__________________
  if (USE_DIIS) then 
     ! Iterative use of Lanczos & DIIS
     call get_environment_variable('Iterative',temporary)
     if (temporary .eq. '') then
        ITERATIVE = .false.
     else
        read(temporary,*) ITERATIVE
     end if

     !!__________________
     ! Number of Lanczos steps after an inflection in the eigenvalue
     call get_environment_variable('Inflection', temporary)
     if (temporary == '') then
        INFLECTION = MAXPAS
     else
        read(temporary,*) INFLECTION
     end if

     !!__________________
     ! Force threshold for call DIIS
     call get_environment_variable('DIIS_Force_Threshold',temporary)
     if (temporary .eq. '') then
        DIIS_FORCE_THRESHOLD = 1.0d0
     else
        read(temporary,*) DIIS_FORCE_THRESHOLD
     end if

     !!__________________
     ! Number of vectors kepts in memory for algorithm
     call get_environment_variable('DIIS_Memory',temporary)
     if (temporary .eq. '') then
        DIIS_MEMORY = 12
     else
        read(temporary,*) DIIS_MEMORY
     end if

     !!__________________
     ! prefactor multiplying forces
     call get_environment_variable('DIIS_Step_size',temporary)
     if (temporary .eq. '') then
        DIIS_STEP = 0.01d0
     else
        read(temporary,*) DIIS_STEP
     end if

     !!__________________
     ! times Increment_Size, max allowed diis step size

     call get_environment_variable('FACTOR_DIIS',temporary)
     if (temporary .eq. '') then
        factor_diis = 5.0
     else
        read(temporary,*) factor_diis
     end if

     !!__________________
     ! max diis iterations per call
     call get_environment_variable('MAX_DIIS',temporary)
     if (temporary .eq. '') then
        maxdiis = 100 
     else
        read(temporary,*) maxdiis
     end if

     !!__________________
     ! Check that the final state is indeed a saddle
     call get_environment_variable('DIIS_Check_Eigenvector',temporary)
     if (temporary .eq. '') then
        DIIS_CHECK_EIGENVEC = .true.
     else
        read(temporary,*) DIIS_CHECK_EIGENVEC
     end if

  else  ! Is this necessary ??
     DIIS_FORCE_THRESHOLD = 1.0d0
     DIIS_STEP            = 1.0d0
     DIIS_MEMORY          = 1
     DIIS_CHECK_EIGENVEC  = .false.
  end if

!SECTION______________________ INPUT/OUTPUT 

  ! File tracking  the file (event) number - facultative
  call get_environment_variable('FILECOUNTER', temporary)
  if (temporary .eq. '') then
     COUNTER  = 'filecounter'
  else
     read(temporary,*) COUNTER
  end if

  !!__________________
  ! General output for message
  call get_environment_variable('LOGFILE', temporary)
  if (temporary .eq. '') then
     LOGFILE   = 'log.file'
  else
     read(temporary,*) LOGFILE
  end if

  !!__________________
  ! list of events with success or failure
  call get_environment_variable('EVENTSLIST', temporary)
  if (temporary .eq. '') then
     EVENTSLIST   = 'events.list'
  else
     read(temporary,*) EVENTSLIST
  end if

  !!__________________
  ! Save the configuration at every step?
  call get_environment_variable('Save_Conf_Int', temporary)
  if (temporary .eq. '') then
     SAVE_CONF_INT = .false. 
  else
     read(temporary,*) SAVE_CONF_INT
  end if

  !!__________________
  ! RESTART: It is useful only for ab-initio
  call get_environment_variable('Write_restart_file', temporary)
  if (temporary .eq. '') then
     write_restart_file = .false.
  else
     read(temporary,*) write_restart_file 
  endif

  !!__________________
  ! current data for restarting event
  call get_environment_variable('RESTART_FILE', temporary)
  if (temporary .eq. '') then
     RESTARTFILE = 'restart.dat'
  else
     read(temporary,*) RESTARTFILE
  end if

  !!__________________
  ! Writes min. and sad. configurations in .xyz format.
  call get_environment_variable('Write_xyz', temporary)
  if (temporary .eq. '') then
     write_xyz = .false.
  else
     read(temporary,*) write_xyz
  endif

  !!__________________
  ! Reference configuration for refine saddle. Without ext.
  call get_environment_variable('REFCONFIG', temporary)
  if (temporary .eq. '') then
     REFCONFIG   = 'refconfig.dat'
  else
     read(temporary,*) REFCONFIG
  end if

  !!__________________HIDE: name for minima files
  call get_environment_variable('FINAL', temporary)
  if (temporary .eq. '') then
     FINAL  = 'min'
  else
     read(temporary,*) FINAL
  end if
  FINAL = adjustl(FINAL)

  !!__________________HIDE: name for saddle file
  call get_environment_variable('SADDLE', temporary)
  if (temporary .eq. '') then
     SADDLE  = 'sad'
  else
     read(temporary,*) SADDLE
  end if

!_____________________________HIDE: Clean_wavefunct
  ! Hide Option, by default .false.
  ! we do not use the previously calculated wave function at some key
  ! points. This is for charged systems.
  call get_environment_variable('Clean_wavefunct',temporary)
  if (temporary == '') then
     clean_wf = .false. 
  else
     read(temporary,*) clean_wf  
     if (.not.( energy_type=="BSW" .or. energy_type=="OTF" .or. &
                energy_type=="BAY" .or. energy_type=="BIG" )     & 
          .and. clean_wf ) then
        write(*,*) "Error : Clean_wavefunct option is only for bigdft"
        stop
     end if
  end if
!_____________________________HIDE: QM/MM 
  call get_environment_variable('MAXNEI', temporary)
  if (temporary .eq. '') then
     maxnei = natoms
  else
     read(temporary,*) maxnei
  end if

   ! number of quantum atoms
  call get_environment_variable('NBR_QUNT', temporary)
  if (temporary .eq. '' .and. energy_type .ne. "BSW" &
       &  .and. energy_type .ne. "OTF" .and. energy_type .ne. "BAY") then
     nbr_quantum = natoms
  elseif ((energy_type == "BSW" .or. energy_type == "OTF" .or. energy_type == "BAY"&
           &) .and. temporary .eq. "" ) then
     write(*,*) "Error : you have not given the number of quantum atoms"
     write(*,*) "The program will stop"
     stop
  elseif (temporary .ne. "") then
     read(temporary,*) nbr_quantum
  endif

   ! number of quantum atoms to trash (buffer zone)
  call get_environment_variable('NBR_QUNT_BUF', temporary)
  if (temporary .ne. '' .and. energy_type .ne. "BSW" .and. energy_type .ne. "OTF" &
            & .and. energy_type .ne. "BAY") then
     write(*,*) "Error: number of quantum atoms only usefull for BSW energy_calc " 
     write(*,*) " we will stop  "
     stop
  elseif ( (energy_type == "BSW" .or. energy_type == "OTF" .or. energy_type == "BAY" &
              &) .and. temporary .eq. "" ) then
     write(*,*) "Error : you have not given the number of quantum atoms"
     write(*,*) "The program will stop"
     stop
  elseif (temporary .eq. '') then
     nbr_quantum_trash = 0
  else
     read(temporary,*) nbr_quantum_trash
  endif

  call get_environment_variable('PASSIVATE', temporary)
  if ( (temporary .eq. ".true.") .and. energy_type .ne. "BSW" .and. energy_type .ne. "OTF" &
              & .and. energy_type .ne. "BAY") then
     write(*,*) "Error: should only passivate if BSW energy_calc " 
     write(*,*) " we will stop  "
     write(*,*) temporary
     stop
  elseif (temporary .eq. "" .or. temporary .ne. ".true." ) then
     passivate = .false.
  else
     read(temporary,*) passivate
  endif

!_____________________________HIDE: local_coord
  if (TYPE_EVENTS == 'local_coord') then

     call get_environment_variable('Radius_Initial_Deformation', temporary)
     if (temporary .eq. '') then
        write(*,*) 'Error: Radius_Initial_Deformation must be defined when TYPE_EVENTS is local'
        stop
     else
        read(temporary,*) LOCAL_CUTOFF
     end if

     call get_environment_variable('Coord_radius', temporary)
     if (temporary .eq. '') then
        write(*,*) 'Error: Coord_radius must be defined when TYPE_EVENTS is local_coord'
        stop
     else
        read(temporary,*) coord_length
     end if

     call get_environment_variable('Coord_number',temporary)
     if (temporary .eq. '') then
        write(*,*) 'Error: Coord_number must be defined when TYPE_EVENTS is local_coord'
        stop
     else
        read(temporary,*) coord_number
     end if

     ! type of atom 
     call get_environment_variable('Type_selected',temporary)
     if (temporary .eq. '') then
        type_sel = 0
     else
        read(temporary,*) type_sel 
     end if
     
  end if
!_____________________________HIDE: Dual_system 

  call get_environment_variable('Dual_system',temporary)
  if (temporary == '') then
     dual_search = .false. 
  else
     read(temporary,*) dual_search 
     if ( .not. (TYPE_EVENTS=='local' .or. TYPE_EVENTS=='list_local') & 
          .and.  dual_search ) then
        write(*,*) 'Error: Dual_system is defined only for local or list_local'
        write(*,*)  dual_search, TYPE_EVENTS
        stop
     end if
  end if

  if ( dual_search ) then
     call get_environment_variable('Size_system', temporary)
     if (temporary .eq. '') then
        write(*,*) 'Error: size_system must be defined if dual_search'
        stop
     else    
        read(temporary,*) size_system
        if ( size_system <= LOCAL_CUTOFF ) then 
           write(*,*) 'Error: Radius_Initial_Deformation > Size_system'
           stop
        end if
     end if
  end if
!_____________________________

  ! We set up other related parameters
  vecsize = 3*natoms

  ! And we allocate the vectors
  allocate(typat(natoms))
  allocate(constr(natoms))
  allocate(force(vecsize))
  allocate(pos(vecsize))
  allocate(posref(vecsize))
  allocate(initial_direction(vecsize))
  allocate(Atom(vecsize))
  allocate(old_projection(VECSIZE))
  allocate(projection(VECSIZE))

  force(:) = 0.0d0
  pos(:) = 0.0d0
  posref(:) = 0.0d0
  initial_direction(:) = 0.0d0
  old_projection(:) = 0.0d0
  projection(:) = 0.0d0

  x => pos(1:NATOMS)
  y => pos(NATOMS+1:2*NATOMS)
  z => pos(2*NATOMS+1:3*NATOMS)

  xref => posref(1:NATOMS)
  yref => posref(NATOMS+1:2*NATOMS)
  zref => posref(2*NATOMS+1:3*NATOMS)

  fx => force(1:NATOMS)
  fy => force(NATOMS+1:2*NATOMS)
  fz => force(2*NATOMS+1:3*NATOMS)

  ! If LANCZOS_MIN, we check how it changes the energy of the system by applying the projection.
  ! This is done at the end of LANCZOS procedure, but this is done only for the minima.
  ! (see subroutine check_min )
  IN_MINIMUN = .False. 
  
END SUBROUTINE read_parameters


!>  ART write_parameters
!! * Define and open the log.file 
!! * Write in it the parameters defining the simulation
subroutine write_parameters( )

  use defs
  use lanczos_defs
  use saddles
  use random
  use module_defs
  implicit none  

  !Local variables
  integer :: ierror, i 
  integer, dimension(8) :: values
  character(len=20) :: fname  ! same len than LOGFILE in defs.f90
  character(len=4)  :: ext
  logical :: exists_already
  !_______________________
                                      ! Initialization of seed in random.
                                      ! WARNING: There will be an idum per
                                      ! iproc. Suggested by Laurent Karim Beland,
                                      ! UdeM 2011
  call date_and_time(values=values)

  if ( .not. setup_initial ) then
     idum = -1 * mod( (1000 * values(7) + values(8))+iproc, 1024)
  else
     idum = 0
  end if

  if ( iproc == 0 ) then

     ! We check what is the last 'LOGFILE.number' in the directory. LOGFILE will be LOGFILE.number+1
    
     i = 1
     do  

        write(ext,'(i4)')  i
        fname = trim(LOGFILE)// "." //adjustl(ext) 
        fname = adjustl(fname)
        inquire( file = fname, exist = exists_already )
         
        if ( .not. exists_already ) then
            logfile = fname
            exit
        end if 
        i = i + 1

     end do  

  ! We write down the various parameters for the simulation
  open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='rewind',iostat=ierror)  
  write(flog,*) '****************************** '
  write(flog,*) 'WELCOME TO BART : BigDFT + ART '
  write(flog,*) '****************************** '
  call timestamp ('Start')
  write(flog,'(1X,A39,f12.3)')  ' - Version number of  ART            : ', VERSION_NUMBER
  write(flog,'(1X,A39,A14)')    ' - Version number of  BIGDFT         : ', package_version 
  write(flog,*) ''
  write(flog,'(1X,A39,A16  )')  ' - Event type                        : ', trim(eventtype)
  write(flog,'(1X,A39,f12.4)')  ' - Temperature                       : ', temperature
  write(flog,'(1X,A39,I12  )')  ' - Number of atoms                   : ', natoms     
  write(flog,'(1X,A39,I12  )')  ' - Number of events                  : ', number_events
  write(flog,'(1X,A39,I12  )')  ' - Maximum number of neighbours      : ', maxnei    

  write(flog,'(1X,A39  )')  ' - Atomic types                      : '
  do i =1, 5
     if (type_name(i) .ne. '') then
        write(flog,'(1X,A31,I3,A5,a12  )')  ' Type ',i,': ', trim(type_name(i))
     end if
  end do

  write(flog,*) ' '
  write(flog,*) 'Selection of the event '
  write(flog,*) '********************** '
  write(flog,'(1X,A39,A15)')   ' - Type of events                    : ', trim(TYPE_EVENTS)
  if (TYPE_EVENTS .eq. 'local') then 
     write(flog,'(1X,A39,f12.4)')   ' - Radius of deformation (local ev.) : ', LOCAL_CUTOFF
     if (preferred_atom .gt. 0) then 
        write(flog,'(1X,A39,i12)')   ' - Central atom for events           : ', preferred_atom
     else
        write(flog,'(1X,A39,A12)')   ' - Central atom for events           :   none (all equal) '
     end if
  end if
  write(flog,*) ' '
  write(flog,*) 'Activation parameters '
  write(flog,*) '********************* '
  write(flog,'(1X,A51,F8.4)')  ' - Eigenvalue threshold                          : ', EIGEN_THRESH
  write(flog,'(1X,A51,F8.4)')  ' - Total force threshold (saddle)                : ', EXITTHRESH
  write(flog,'(1X,A51,F8.4)')  ' - Initial step size                             : ', INITSTEPSIZE
  write(flog,'(1X,A51,F8.4)')  ' - Increment size                                : ', INCREMENT
  write(flog,'(1X,A51,F8.4)')  ' - Atomic displacement for breaking the symmetry : ', sym_break_dist
  write(flog,'(1X,A51,2I8)')    ' - Number of vectors computed by Lanzcos         : ', NVECTOR_LANCZOS_H,&
  NVECTOR_LANCZOS_C
  write(flog,'(1X,A51,F8.4)')  ' - Delta displacement for derivative in Lanzcos  : ', DEL_LANCZOS
  write(flog,*) ' '
  write(flog,'(1X,A51,L8)')    ' - Calculation of the Hessian for each minimum   : ', LANCZOS_MIN
  write(flog,'(1X,A51,I8)')    ' - Min. number of ksteps before calling lanczos  : ', KTER_MIN
  write(flog,'(1X,A51,F8.4)')  ' - Factor mult. INCREMENT for leaving basin      : ', BASIN_FACTOR
  write(flog,'(1X,A51,I8)')    ' - Maximum number of iteractions (basin -kter)   : ', MAXKTER
  write(flog,'(1X,A51,I8)')    ' - Maximum number of perpendicular moves (basin) : ', MAXKPERP
  write(flog,*) ' '
  write(flog,'(1X,A51,I8)')    ' - Maximum number of iteractions (activation)    : ', MAXPAS
  write(flog,'(1X,A51,I8)')    ' - Maximum number of perpendicular moves (activ) : ', MAXIPERP
  write(flog,'(1X,A51,F8.4)')  ' - Force threshold for perpendicular relaxation  : ', FTHRESHOLD
  write(flog,'(1X,A51,F8.4)')  ' - Fraction of displacement over the saddle      : ', PUSH_OVER

  write(flog,*) ' '
  write(flog,'(1X,A39,L12)')    ' Use DIIS for convergence to saddle     : ', USE_DIIS
  if (USE_DIIS) then
     write(flog,'(1X,A40,F12.4)') ' - Total force threshold to call DIIS : ', DIIS_FORCE_THRESHOLD
     write(flog,'(1X,A40,F12.4)') ' - Step size to update positions      : ', DIIS_STEP
     write(flog,'(1X,A40,I12)')   ' - Memory (number of steps)           : ', DIIS_MEMORY
     write(flog,'(1X,A40,L12)')   ' - Check eigenvector at saddle        : ', DIIS_CHECK_EIGENVEC
  end if

  write(flog,*) ' '
  write(flog,*) 'Input / Output '
  write(flog,*) '********************* '
  write(flog,'(1X,A39,A15)')  ' - Name of log file                  : ', trim(LOGFILE) 
  write(flog,'(1X,A39,A15)')  ' - Liste of events                   : ', trim(eventslist) 
  write(flog,'(1X,A39,A15)')  ' - Reference configuration           : ', trim(refconfig) 
  write(flog,'(1X,A39,A15)')  ' - Restart file                      : ', trim(restartfile)
  write(flog,'(1X,A39,A15)')  ' - Prefix for minima (file)          : ', trim(FINAL)      
  write(flog,'(1X,A39,A15)')  ' - Prefix for saddle points          : ', trim(SADDLE)
  write(flog,'(1X,A39,A15)')  ' - File with filecounter             : ', trim(counter)
  write(flog,*) '********************* '
  write(flog,*) ' '
  write(flog,'(1X,A34,I17)') ' - The seed is                  : ', idum
  write(flog,'(1X,A34,L17)') ' - Restart                      : ', restart
  close(flog)

  end if

END SUBROUTINE write_parameters

!> ART timestamp
!! it stamps the time at the beggining and at the final of the simulation
subroutine timestamp (str)
   
  use defs, only : FLOG 
  implicit none
  !Arguments   
  character(len=*), intent(in) :: str
  !Local variables    
  integer         :: sec, min, hour, day, month, year
  character(len=1), parameter :: dash = "-"
  character(len=1), parameter :: colon = ":"
  character(len=3), parameter :: prefix = ">> "
  character(len=3)  :: month_str(12) =     &
    (/'JAN','FEB','MAR','APR','MAY','JUN', &
      'JUL','AUG','SEP','OCT','NOV','DEC'/)

  integer :: values(8)

  call date_and_time(values=values)
  year  = values(1)
  month = values(2)
  day   = values(3)
  hour  = values(5)
  min   = values(6)
  sec   = values(7)
  
  write(FLOG,1000)  prefix, trim(str), colon, day, dash, month_str(month),&
  &                 dash, year, hour, colon, min, colon, sec
  
1000 format(2a,a1,2x,i2,a1,a3,a1,i4,2x,i2,a1,i2.2,a1,i2.2)
    
END SUBROUTINE timestamp
