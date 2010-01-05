
subroutine read_parameters()

  use defs
  use lanczos_defs
  use saddles

  implicit none  
  integer :: ierror,i 
  logical :: exists_already
  character(len=40)  :: temporary 
  !character(len=20)  :: eventtype
  character(len=100) :: fname
  character(len=150) :: commande
  character(len=9)   :: digit = "123456789"

  ! We first read the parameters defining the run


  call getenv('EVENT_TYPE', temporary)
  if (temporary .eq. '') then
     new_event = .true.
     eventtype = 'NEW'
  else if (temporary .eq. 'NEW') then
     new_event = .true.
     eventtype = 'NEW'
  else if (temporary .eq. 'REFINE_SADDLE') then
     new_event = .false.
     eventtype = 'REFINE_SADDLE'
  else if (temporary .eq. 'REFINE_AND_RELAX') then
     new_event = .false.
     eventtype = 'REFINE_AND_RELAX'
  else
     write(*,*) 'Error: event_types permitted are NEW , REFINE_SADDLE and REFINE_AND_RELAX'
     stop
  endif

  call getenv('Temperature', temporary)
  if (temporary .eq. '') then
     write(*,*) 'Error: Metropolis temperature is not defined'
     stop
  else
     read(temporary,*) temperature
  endif

  call getenv('NATOMS', temporary)
  if (temporary .eq. '') then
     write(*,*) 'Error: NATOMS is not defined'
     stop
  else
     read(temporary,*) natoms
  endif

  call getenv('MAXNEI', temporary)
  if (temporary .eq. '') then
     maxnei = natoms
  else
     read(temporary,*) maxnei
  endif

  call getenv('Max_Number_Events', temporary)
  if (temporary .eq. '') then
     number_events = 100
  else
     read(temporary,*) number_events
  endif

  call getenv('MAXNEI', temporary)
  if (temporary .eq. '') then
     maxnei = natoms
  else
     read(temporary,*) maxnei
  endif


  ! File names
  call getenv('LOGFILE', temporary)
  if (temporary .eq. '') then
     LOGFILE   = 'log.file'
  else
     read(temporary,*) LOGFILE
  endif

  call getenv('EVENTSLIST', temporary)
  if (temporary .eq. '') then
     EVENTSLIST   = 'events.list'
  else
     read(temporary,*) EVENTSLIST
  endif

  call getenv('REFCONFIG', temporary)
  if (temporary .eq. '') then
     REFCONFIG   = 'refconfig.dat'
  else
     read(temporary,*) REFCONFIG
  endif

  call getenv('FINAL', temporary)
  if (temporary .eq. '') then
     FINAL  = 'min'
  else
     read(temporary,*) FINAL
  endif

  call getenv('SADDLE', temporary)
  if (temporary .eq. '') then
     SADDLE  = 'sad'
  else
     read(temporary,*) SADDLE
  endif

  call getenv('FILECOUNTER', temporary)
  if (temporary .eq. '') then
     COUNTER  = 'filecounter'
  else
     read(temporary,*) COUNTER
  endif

  call getenv('RESTART_FILE', temporary)
  if (temporary .eq. '') then
     RESTARTFILE = 'restart.dat'
  else
     read(temporary,*) RESTARTFILE
  endif


  ! Reading details for activation
  ! Read type of events - local or global
  call getenv('Type_of_Events', TYPE_EVENTS)
  if ( (TYPE_EVENTS .ne. 'global') .and. (TYPE_EVENTS .ne. 'local')) then
     write(*,*) 'Error : only global or local type of events are accepted - provided: ', TYPE_EVENTS
     stop
  endif

  if (TYPE_EVENTS .eq. 'local') then
     call getenv('Radius_Initial_Deformation', temporary)
     if (temporary .eq. '') then
        write(*,*) 'Error: Radius_Initial_Deformation must be defined when TYPE_EVENTS is local'
        stop
     else
        read(temporary,*) LOCAL_CUTOFF
     endif

     call getenv('Central_Atom',temporary)
     if (temporary .eq. '') then
       preferred_atom = -1
     else
       read(temporary,*) preferred_atom
     endif
  endif

  ! Info regarding initial displacement
  call getenv('Initial_Step_Size', temporary)
  if (temporary .eq. '') then
     INITSTEPSIZE = 0.001             ! Size of initial displacement in Ang.
  else
     read(temporary,*) INITSTEPSIZE
  endif

  ! Minimum number of steps in kter-loop before we call lanczos
  call getenv('Min_Number_KSteps',temporary)
  if (temporary .eq. '') then
     KTER_MIN = 1
  else
     read(temporary,*) KTER_MIN
  endif
  
  ! Size of the parallel displacement in for leavign the basin with respect to
  ! the overall increment
  call getenv('Basin_Factor',temporary)
  if (temporary .eq. '') then
    BASIN_FACTOR = 1.0d0
  else
    read(temporary,*) BASIN_FACTOR
  endif

  ! Maximum number of iteration in the activation
  call getenv('Max_Iter_Activation',temporary)
  if (temporary .eq. '') then
     MAXITER = 1000
  else
     read(temporary,*) MAXITER
  endif

  ! Maximum number of iteration in the activation
  call getenv('Max_Iter_Basin',temporary)
  if (temporary .eq. '') then
     MAXKTER = 100
  else
     read(temporary,*) MAXKTER
  endif


  ! Number of relaxation perpendicular moves - basin and activation
  call getenv('Max_Perp_Moves_Basin',temporary)
  if (temporary .eq. '') then
     MAXKPERP = 2
  else
     read(temporary,*) MAXKPERP
  endif

  call getenv('Max_Perp_Moves_Activ',temporary)
  if (temporary .eq. '') then
     MAXIPERP = 12
  else
     read(temporary,*) MAXIPERP
  endif

  ! Increment size - overall scaling (in angstroems)
  call getenv('Increment_Size',temporary)
  if (temporary .eq. '') then
     INCREMENT = 0.01
  else
     read(temporary,*) INCREMENT
  endif

  ! Eigenvalue threshold
  call getenv('Eigenvalue_Threshold',temporary)
  if (temporary .eq. '') then
     write(*,*) 'Error : No eigenvalue threshold provided  (Eigenvalue_Threshold)'
     stop
  else
     read(temporary,*) EIGEN_THRESH
  endif

  ! Force threshold for the perpendicular relaxation
  call getenv('Force_Threshold_Perp_Rel',temporary)
  if (temporary .eq. '') then
     FTHRESHOLD = 1.0
  else
     read(temporary,*) FTHRESHOLD
  endif
  FTHRESH2 = FTHRESHOLD * FTHRESHOLD

  ! Force threshold for the convergence at the saddle point
  call getenv('Exit_Force_Threshold',temporary)
  if (temporary .eq. '') then
     EXITTHRESH = 1.0
  else
     read(temporary,*) EXITTHRESH
  endif


  ! Force threshold for the convergence at the saddle point
  call getenv('Number_Lanczos_Vectors',temporary)
  if (temporary .eq. '') then
     NVECTOR_LANCZOS = 16
  else
     read(temporary,*) NVECTOR_LANCZOS
  endif

  ! The prefactor for pushing over the saddle point, fraction of distance from
  ! initial minimum to saddle point
  call getenv('Prefactor_Push_Over_Saddle',temporary)
  if (temporary .eq. '') then
    PUSH_OVER = 0.15
  else
    read(temporary,*) PUSH_OVER
  endif


  ! We now get the types - define up to 5
  call getenv('type1',temporary)
  if (temporary .eq. '') then
     write(*,*) 'Error, must at least define 1 type of atoms -  use "setenv type1 Si", for example'
     stop
  else
     read(temporary,*) type_name(1)
  endif

  call getenv('type2',temporary)
  if (temporary .eq. '') then
     type_name(2) = ''
  else
     read(temporary,*) type_name(2)
  endif

  call getenv('type3',temporary)
  if (temporary .eq. '') then
     type_name(3) = ''
  else
     read(temporary,*) type_name(3)
  endif

  call getenv('type4',temporary)
  if (temporary .eq. '') then
     type_name(4) = ''
  else
     read(temporary,*) type_name(4)
  endif

  call getenv('type5',temporary)
  if (temporary .eq. '') then
     type_name(5) = ''
  else
     read(temporary,*) type_name(5)
  endif

  call getenv('type6',temporary)
  if (temporary .ne. '') then
     write(*,*) 'Error: The code can only handle 5 atomic types, change read_parameters.f90, to allow for more.'
     stop
  endif

  call getenv('sym_break_dist',temporary)
  if (temporary .eq. '') then
     sym_break_dist = 0.0
  else
     read(temporary,*) sym_break_dist
  endif

  ! Do we use DIIS for refining and converging to saddle
  call getenv('Use_DIIS',temporary)
  if (temporary .eq. '') then
     USE_DIIS = .false.
  else
     read(temporary,*) USE_DIIS
  endif

  ! If diis is used, we define a number of parameters
  if (USE_DIIS) then
     call getenv('DIIS_Force_Threshold',temporary)
     if (temporary .eq. '') then
        DIIS_FORCE_THRESHOLD = 0.1d0
     else
        read(temporary,*) DIIS_FORCE_THRESHOLD
     endif

     call getenv('DIIS_Memory',temporary)
     if (temporary .eq. '') then
        DIIS_MEMORY = 12
     else
        read(temporary,*) DIIS_MEMORY
     endif

     call getenv('DIIS_MaxIter',temporary)
     if (temporary .eq. '') then
        DIIS_MAXITER = 50
     else
        read(temporary,*) DIIS_MAXITER
     endif

     call getenv('DISS_Check_Eigenvector',temporary)
     if (temporary .eq. '') then
        DIIS_CHECK_EIGENVEC = .true.
     else
        read(temporary,*) DIIS_CHECK_EIGENVEC
     endif

     call getenv('DIIS_Step_Size',temporary)
     if (temporary .eq. '') then
        DIIS_STEP = 0.01d0
     else
        read(temporary,*) DIIS_STEP
     endif
  else
     DIIS_FORCE_THRESHOLD = 1.0d0
     DIIS_STEP = 1.0d0
     DIIS_MEMORY          = 1
     DIIS_MAXITER         = 1
     DIIS_CHECK_EIGENVEC  = .false.
  endif
  
  ! We set up other related parameters
  vecsize = 3*natoms

  ! And we allocate the vectors
  allocate(type(natoms))
  allocate(force(vecsize))
  allocate(pos(vecsize))
  allocate(posref(vecsize))
  allocate(direction_restart(vecsize))
  allocate(initial_direction(vecsize))
  allocate(Atom(vecsize))

  ! Vectors for lanczos
  allocate(old_projection(VECSIZE))
  allocate(projection(VECSIZE))
  allocate(first_projection(VECSIZE))

  x => pos(1:NATOMS)
  y => pos(NATOMS+1:2*NATOMS)
  z => pos(2*NATOMS+1:3*NATOMS)

  xref => posref(1:NATOMS)
  yref => posref(NATOMS+1:2*NATOMS)
  zref => posref(2*NATOMS+1:3*NATOMS)

  fx => force(1:NATOMS)
  fy => force(NATOMS+1:2*NATOMS)
  fz => force(2*NATOMS+1:3*NATOMS)

  ! We first check whether the file "LOGFILE" exists. If so, we copy it before
  ! we start.
  fname = LOGFILE
  do i=1, 9 
    inquire(file=fname,exist=exists_already)
    if (exists_already)  then
       fname = trim(LOGFILE) // "." // digit(i:i)
    else 
       if (i .gt. 1 ) then
         commande = "mv " // LOGFILE // "  " // fname
         call system(commande)
       endif
       exit
    endif 
  end do

  ! We write down the various parameters for the simulation
  open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='rewind',iostat=ierror)  
  write(flog,'(A39,f16.3)')  ' Version number of siestart           : ', VERSION_NUMBER
  write(flog,*) ' '
  write(flog,'(A39,A16  )')  ' Event type                           : ', eventtype
  write(flog,'(A39,f16.4)')  ' Temperature                          : ', temperature
  write(flog,'(A39,I16  )')  ' Number of atoms                      : ', natoms     
  write(flog,'(A39,I16  )')  ' Number of events                     : ', number_events
  write(flog,'(A39,I16  )')  ' Maximum number of neighbours         : ', maxnei    

  write(flog,'(A39  )')  ' Atomic types                         : '
  do i =1, 5
     if (type_name(i) .ne. '') then
        write(flog,'(A31,i1,A5,a4  )')  '                          Type ',i,'   : ', type_name(i)
     endif
  end do

  write(flog,*) ' '
  write(flog,*) 'Selection of the event '
  write(flog,*) '********************** '
  write(flog,'(1X,A39,A12)')   ' - Type of events                    : ', TYPE_EVENTS
  if (TYPE_EVENTS .eq. 'local') then 
    write(flog,'(1X,A39,f12.4)')   ' - Radius of deformation (local ev.) : ', LOCAL_CUTOFF
    if (preferred_atom .gt. 0) then 
      write(flog,'(1X,A39,i12)')   ' - Central atom for events           : ', preferred_atom
    else
      write(flog,'(1X,A39,A12)')   ' - Central atom for events           :   none (all equal) '
    endif
  endif
  write(flog,*) ' '
  write(flog,*) 'Activation parameters '
  write(flog,*) '********************* '
  write(flog,'(1X,A39,F12.4)') ' - Eigenvalue threshold              : ', EIGEN_THRESH
  write(flog,'(1X,A39,F12.4)') ' - Total force threshold (saddle)    : ', EXITTHRESH

  write(flog,'(1X,A39,F12.4)') ' - Initial step size                 : ', INITSTEPSIZE
  write(flog,'(1X,A39,F12.4)') ' - Increment size                    : ', INCREMENT
  write(flog,'(1X,A51,I8   )') ' - Number of vectors computed by Lanzcos         : ', NVECTOR_LANCZOS
  write(flog,*) ' '
  write(flog,'(1X,A51,I8)')    ' - Min. number of ksteps before calling lanczos  : ', KTER_MIN
  write(flog,'(1X,A51,F8.4)')  ' - Factor mulp. INCREMENT for leaving basin      : ', BASIN_FACTOR
  write(flog,'(1X,A51,I8)')    ' - Maximum number of iteractions (basin -kter)   : ', MAXKTER
  write(flog,'(1X,A51,I8)')    ' - Maximum number of perpendicular moves (basin) : ', MAXKPERP
  write(flog,*) ' '
  write(flog,'(1X,A51,I8)')    ' - Maximum number of iteractions (activation)    : ', MAXITER
  write(flog,'(1X,A51,I8)')    ' - Maximum number of perpendicular moves (activ) : ', MAXIPERP
  write(flog,'(1X,A51,F8.4)')  ' - Force threshold for perpendicular relaxation  : ', FTHRESHOLD
  write(flog,'(1X,A51,F8.4)')  ' - Fraction of displacement over the saddle      : ', PUSH_OVER

  write(flog,*) ' '
  write(flog,'(1X,A36,L6)')    'Use DIIS for convergence to saddle : ', USE_DIIS
  if (USE_DIIS) then
     write(flog,'(1X,A36,F30.4)') ' - Total force threshold           : ', DIIS_FORCE_THRESHOLD
     write(flog,'(1X,A36,F30.4)') ' - Step size to update positions   : ', DIIS_STEP
     write(flog,'(1X,A36,I15)')   ' - Memory (number of steps)        : ', DIIS_MEMORY
     write(flog,'(1X,A36,I15)')   ' - Maximum number of iterations    : ', DIIS_MAXITER
     write(flog,'(1X,A36,L6)')    ' - Check eigenvector at saddle     : ', DIIS_CHECK_EIGENVEC
  endif


  write(flog,*) ' '
  write(flog,*) 'Input / Output '
  write(flog,*) '********************* '
  write(flog,'(A39,A16  )')  ' Name of log file                     : ', logfile
  write(flog,'(A39,A16  )')  ' Liste of events                      : ', eventslist 
  write(flog,'(A39,A16  )')  ' Reference configuration              : ', refconfig   
  write(flog,'(A39,A16  )')  ' Restart file                         : ', restartfile
  write(flog,'(A39,A16  )')  ' Prefix for minima (file)             : ', FINAL      
  write(flog,'(A39,A16  )')  ' Prefix for saddle points             : ', SADDLE
  write(flog,'(A39,A16  )')  ' File with filecounter                : ', counter

  close(flog)
end subroutine
