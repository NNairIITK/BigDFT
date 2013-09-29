!> @file
!! Define a module for BigDFT+ART
!! @author
!!    Eduardo Machado-Charry (EM) 2010
!!    Copyright (C) 2001 Normand Mousseau
!!    Copyright (C) 2010-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Module containing the DIIS global arrays for ART
module diis_defs

  implicit none

  real(kind=8), dimension(:,:), allocatable :: previous_forces
  real(kind=8), dimension(:,:), allocatable :: previous_pos
  real(kind=8), dimension(:),   allocatable :: previous_norm
  integer :: maxter

END MODULE diis_defs


!> ART allocate_activation
!! Allocation of ART-DIIS arrays
subroutine allocate_activation ()

  use defs
  use diis_defs
  implicit none

  allocate(previous_forces(DIIS_MEMORY,VECSIZE))
  allocate(previous_pos(DIIS_MEMORY,VECSIZE))
  allocate(previous_norm(DIIS_MEMORY))

  ! Initialization.
  previous_forces = 0.0d0  
  previous_pos    = 0.0d0
  previous_norm   = 0.0d0
  maxter = 0

END SUBROUTINE allocate_activation


!> ART deallocate_activation 
!! Deallocation of ART-DIIS arrays
subroutine deallocate_activation ()

  use diis_defs
  implicit none

  deallocate(previous_forces)
  deallocate(previous_pos)
  deallocate(previous_norm)

END SUBROUTINE deallocate_activation


!> ART apply_diis
subroutine apply_diis( diter, saddle_energy, ret )

  use defs
  use diis_defs
  use bigdft_forces
  use lanczos_defs,  only: projection, NVECTOR_LANCZOS_C, eigenvalue
  use saddles
  implicit none

  !Arguments
  integer,      intent(inout) :: diter
  real(kind=8), intent(out)   :: saddle_energy 
  integer,      intent(inout) :: ret

  !Local variables
  integer :: i
  real(kind=8) :: a1, smallest_error, invse
  real(kind=8) :: n_deltaGdiis
  real(kind=8), dimension(3) :: boxl
  real(kind=8), dimension(VECSIZE) :: deltaGdiis
  real(kind=8), dimension(:),   allocatable :: solution
  real(kind=8), dimension(:,:), allocatable :: error_vector
  logical :: new_projection 
  logical :: rejected_step
  !real(kind=8) :: costheta, n_deltaref,  degtheta, sumpos, sumneg ! For Farkas' Tests.
  !real(kind=8), dimension(VECSIZE) :: deltaref ! For Farkas' Tests.
  !_______________________
  boxl = box * scala

  allocate(error_vector(DIIS_MEMORY,VECSIZE))
  error_vector = 0.0d0

  if ( .not. restart ) then 
     ! We start with two configurations given by the last two lanczos steps :
     ! * They were saved as previous_*(1 & 2 ,:) in apply_lanczos. 
     maxter = 2                       ! Initialization of maxter
  else
     restart = .false.
  end if

  While_diis: do 

    ! Test of While_diis loop 
    if ( (ftot < EXITTHRESH) .or. (pas > MAXPAS) .or. &
       & (diter > maxdiis) .or.                       &
       & ( delta_e < delta_thr .and. delr < delr_thr .and. pas > (MAXKTER+5) ) ) then

       saddle_energy = total_energy
       end_activation = .true.
       pas = pas - 1                  ! the real number of steps in the event 

       if ( ftot < EXITTHRESH ) then
          ret = 20000 + pas
       else if  (delta_e < delta_thr .and. delr < delr_thr .and. pas > (MAXKTER+5) ) then
          ret = 80000 + pas 
       else if ( diter > maxdiis) then
          ret = 70000 + pas
       else if ( pas > maxpas ) then 
          ret = 50000 + pas 
       end if 

       exit
    end if

    ! We test a solution until either it is accepted or we do not have memory.
    do_rejected: do 
      ! "To avoid the effect of the size of the error vectors, we rescale them
      !  for the purpose of solving the DIIS equations so that the smallest error
      !  has unit length" Farkas,Schlegel. Phys.Chem.Chem.Phys.2002,4,11-15
      smallest_error = minval(previous_norm, MASK = previous_norm > 0.0)
      invse = DIIS_STEP/(smallest_error)
      do i = 1, maxter
         error_vector(i,:) = previous_forces(i,:)*invse
      end do

                                      ! Let's calculate the coefficients
      allocate(solution(maxter+1))
      call get_solution( maxter, error_vector, solution )

                                      ! The solution that interests us is made of two parts
      pos(:)   = 0.0d0
      force(:) = 0.0d0
      do i = 1, maxter
         pos(:)   = pos(:)   + solution(i) * previous_pos(i,:)
         force(:) = force(:) + solution(i) * previous_forces(i,:)*DIIS_STEP
      end do
                                      ! the DIIS_step 
      pos(:) = pos(:) + force(:)
                                      ! Let's be sure about frozen atoms
      do i = 1, natoms
         if ( constr(i) .ne. 0 ) then 
            pos(i)          = previous_pos(1,i)
            pos(i+natoms)   = previous_pos(1,i+natoms)
            pos(i+2*natoms) = previous_pos(1,i+2*natoms)
         end if
      end do  
      
      ! Can we accept this DIIS_step ?
      ! ----------------------

      ! Boundary conditions: Suggested by Laurent Karim Beland, UdeM 2011
      call boundary_cond( deltaGdiis, pos, previous_pos(maxter,:) )   
      n_deltaGdiis = sqrt(dot_product( deltaGdiis,deltaGdiis ))

      !deltaref(:) = DIIS_STEP* previous_forces(maxter,:)
      !n_deltaref   = sqrt(dot_product( deltaref,deltaref ))
      !costheta= dot_product(deltaGdiis,deltaref)/(n_deltaGdiis * n_deltaref)
      !degtheta= acos(costheta)*360/(8.d0*atan(1.d0))
      !sumpos = sum ( solution, MASK = solution > 0.0 ) 
      !sumneg = sum ( solution, MASK = solution < 0.0 )
      !if ( n_deltaGdiis/n_deltaref > 15.0d0 .or.       &
      !   & sumpos > 15.0d0 .or. solution(maxter+1) < 1.0E-08 ) then 
      !  rejected_step = .true.    

      if ( solution(maxter+1) < 1.0E-08 .or. n_deltaGdiis > factor_diis*INCREMENT) then
         rejected_step = .true.
      else
         rejected_step = .false.
      end if 
      if ( iproc == 0 ) then          ! REPORT DIIS 
        write (*,'(a,3I5,1x,2F11.4,2x,(1p,e14.5,0p),1x,L1)') 'BART DIIS:',  &
        &  pas, diter, maxter, n_deltaGdiis, factor_diis*INCREMENT, solution(maxter+1), &
        &  rejected_step 
      end if

      deallocate(solution)

      If_rejected: if ( rejected_step ) then 

         ! If the step is rejected we reduce the memory by one.
         do i = 1, maxter - 1
            previous_forces(i,:) = previous_forces(i+1,:)
            previous_pos(i,:)    = previous_pos(i+1,:)
            previous_norm(i)     = previous_norm(i+1)
         end do
         previous_forces(maxter,:) = 0.d0
         previous_pos(maxter,:)    = 0.d0
         previous_norm(maxter)     = 0.d0
         maxter = maxter - 1
  
         if ( maxter < 2 ) then   ! DIIS does not have memory to continue

            switchDIIS = .False.
            if ( .not. ITERATIVE ) then         

               ! If not ITERATIVE, we can not continue with this event. 
               end_activation = .true. 
               ret = 70000 + pas - 1
               if ( iproc == 0 ) write(*,*) 'BART: DIIS Failed, no memory'

            else if ( diter > 1 ) then

               ! If we have accepted at least one DIIS step, we check the
               ! projection vector up to four times if the resulting lowest 
               ! eigenvalue is positive at each iteration.
                                          ! clean_wf
               if ( clean_wf ) then
                  pas = pas - 1           ! only for the report 
                  call clean_wavefunction ( 0.0d0, .False. )
                  pas = pas + 1
               end if 

               new_projection = .false.   ! Let's start with the last projection.  

               Do_lanc: do i = 1, 4
                  call lanczos( NVECTOR_LANCZOS_C, new_projection, a1 )
                  if ( eigenvalue < 0.0d0 ) exit Do_lanc
               end do Do_lanc                         
                                      ! Orientation of the projection vector                               
               fpar = dot_product( force, projection )               
               if ( fpar > 0.0d0 ) projection = -1.0d0 * projection  
            else  
                                      ! This is the projection vector for the second 
                                      ! trial step when we called DIIS            
               new_projection = .false.       
               call lanczos( NVECTOR_LANCZOS_C, new_projection, a1 )
            end if
            
            return ! c'est la fin. If ITERATIVE Let's try again with Lanczos.

         end if ! if (maxter<2)

      else ! else of If_rejected. The solution is accepted. 
         
         ! Now the force is evaluted for our new pos configuration. 
         call calcforce( NATOMS, pos, boxl, force, total_energy, evalf_number, .false. )
         ftot = sqrt(dot_product(force,force))

         maxter = maxter + 1          ! update of history 

         ! If maxter is greater than DIIS_MEMORY, we move the previous
         ! solution up by one
         if ( maxter .gt. DIIS_MEMORY ) then
            maxter = DIIS_MEMORY 
            do i= 2, maxter
               previous_forces(i-1,:) = previous_forces(i,:)
               previous_pos(i-1,:)    = previous_pos(i,:)
               previous_norm(i-1)     = previous_norm(i)
            end do
         end if

         ! we add the last accepted configuration to our space
         previous_forces(maxter,:) = force(:)
         previous_pos(maxter,:)    = pos(:)
         previous_norm(maxter)     = ftot*DIIS_STEP

         call displacement( posref, pos, delr, npart )
         delta_e = total_energy - ref_energy

         ! write & save: Some values are meaningless in a DIIS step.
         m_perp = 0
         try = 0
         fpar = 0.0d0
         fperp= 0.0d0
         eigen_min = 0.0d0
         eigenvalue = 0.0d0
         nsteps_after_eigen_min = 0
         a1 = 0.0d0
         call write_step ( 'D', diter, a1, total_energy )

         if ( SAVE_CONF_INT ) call save_intermediate( 'D' ) 
                                      ! For restart ( restart.f90 )
         if ( write_restart_file ) then 
            state_restart = 4 
            if ( iproc == 0 ) then 
               call save_state2(state_restart,diter+1,projection,maxter,previous_forces,&
               & previous_pos,previous_norm,eigen_min,eigenvalue,nsteps_after_eigen_min )
            end if
         end if

         exit ! exit do_rejected

      end if If_rejected

    end do do_rejected

    saddle_energy = total_energy
                                                  ! clean_wf
    if ( clean_wf .and. ftot < EXITTHRESH ) then

       call clean_wavefunction ( 0.0d0, .False. ) 
                                                  ! The total force has changed.
       if ( ftot >= EXITTHRESH .and. cw_try_again )  then
                                                  ! We make Lanczos.
          cw_try_again = .False.                  ! But, this is done only once per event.
          switchDIIS = .False.
          previous_forces = 0.0d0                 ! Clean DIIS memory.
          previous_pos    = 0.0d0
          previous_norm   = 0.0d0
                                                  ! We need a projection vector
          new_projection = .false.                ! Let's start with the last projection.  
          Do_lanc_wf: do i = 1, 4
             call lanczos( NVECTOR_LANCZOS_C, new_projection, a1 )
             if ( eigenvalue < 0.0d0 ) exit Do_lanc_wf
          end do Do_lanc_wf                         
                                                  ! Orientation of the projection vector                               
          fpar = dot_product( force, projection )               
          if ( fpar > 0.0d0 ) projection = -1.0d0 * projection  
          return
       else
          ret = 40000 + pas
          saddle_energy = total_energy
          end_activation = .true.
          return 
       end if 
    end if

    pas = pas + 1
    diter = diter + 1

  end do While_diis

END SUBROUTINE apply_diis


!> ART get_solution 
!! of DIIS coefficients
subroutine get_solution( maxter, error_vector, solution )

  use defs, only: DIIS_MEMORY, VECSIZE, iproc
  implicit none

  !Arguments
  integer, intent(in) :: maxter
  real(kind=8), dimension(DIIS_MEMORY,VECSIZE), intent(in) :: error_vector
  real(kind=8), dimension(maxter+1), intent(out) :: solution

  ! Local variables for Lapack.
  integer :: i, j
  integer :: i_err, lwork, n, nrhs
  integer,      dimension(:),   allocatable :: interchanges
  real(kind=8), dimension(:),   allocatable :: work
  real(kind=8), dimension(:,:), allocatable :: matrice
  !_______________________
  ! Initialization.
  allocate(matrice(maxter+1,maxter+1))
  allocate(interchanges(maxter+1))
  solution = 0.d0
  matrice  = 0.d0

  ! We prepare the upper triangular matrix for lapack
  do i = 1, maxter
     do j = i, maxter
        matrice(i,j) = dot_product( error_vector(i,:), error_vector(j,:) )
        matrice(j,i) = matrice(i,j) 
     end do
  end do

  matrice(1:maxter,maxter+1) = -1.0d0
  matrice(maxter+1,1:maxter) = -1.0d0
  matrice(maxter+1,maxter+1) =  0.0d0

  solution(1:maxter) =  0.0d0
  solution(maxter+1) = -1.0d0

  ! We call the routine for diagonalizing a tridiagonal  matrix
  i_err = 0        ! Initialization of lapack infocode
  n = maxter + 1   ! size of arrays
  nrhs = 1         ! Number of right-hand arguments in B (Ax=B)
  
                   ! We find the best size for "work" array.
  allocate(work(100))
  call dsysv('U',n, nrhs, matrice, n, interchanges, solution,n,work,-1,i_err)
  lwork=int(work(1))
  deallocate(work)
  allocate(work(lwork))
  
                   ! We call the routine for diagonalizing a tridiagonal  matrix
  call dsysv('U', n, nrhs, matrice, n, interchanges, solution, n, work, lwork, i_err )
                   ! If something fails 
  if ( i_err /= 0 ) then
     if ( iproc == 0 ) write(*,*) 'BART WARNING DIIS: info calculation of solution', i_err
  end if

  deallocate(work)
  deallocate(interchanges)
  deallocate(matrice)

END SUBROUTINE get_solution


!> ART apply_lanczos
!! It is the loop over the hyperplanes. Once we decide to switch to DIIS we set up
!! the first two vectors of the memory.
subroutine  apply_lanczos ( liter, saddle_energy, ret )

  use defs
  use lanczos_defs
  use diis_defs
  use saddles
  implicit none

  !Arguments
  integer,      intent(inout) :: liter  ! # lanczos steps
  real(kind=8), intent(out)   :: saddle_energy 
  integer,      intent(inout) :: ret

  !Local variables
  logical      :: get_proj
  real(kind=8) :: a1
  !_______________________ 

  if ( .not. restart ) then 
     ! To be consistent with the restart file, each time we call
     ! this subroutine we make this initialization.
     previous_forces = 0.0d0  
     previous_pos    = 0.0d0
     previous_norm   = 0.0d0
     maxter = 0
     eigen_min = 0.0d0
     nsteps_after_eigen_min = 0
  else
     ! A warranty 
     eigen_min = min( eigen_min, eigenvalue )
     if ( eigenvalue == eigen_min ) then
        previous_forces = 0.0d0  
        previous_pos    = 0.0d0
        previous_norm   = 0.0d0
     end if
     restart = .false.
  end if

  a1 = 0.0d0

  While_lanczos: do
 
      ! Test of While_lanczos loop 
      if ( (ftot < EXITTHRESH) .or. (pas > MAXPAS) .or. (liter > 280) .or. & 
         & (eigenvalue > 0.0)  .or. &  
         & ( delta_e < delta_thr .and. delr < delr_thr .and. pas > (MAXKTER+5) ) ) then

         saddle_energy = total_energy 
         end_activation = .true.  
         pas = pas - 1                ! the real number of steps in the event 

         if ( ftot < EXITTHRESH ) then
            ret = 20000 + pas
         else if  (delta_e < delta_thr .and. delr < delr_thr .and. pas > (MAXKTER+5) ) then
            ret = 80000 + pas 
         else if (eigenvalue > 0.0) then
            ret = 60000 + pas
         else if ( pas > maxpas .or. liter > 280 ) then 
            ret = 50000 + pas 
         end if 

         exit
      end if 

      if ( USE_DIIS ) then 
                                      ! Criteria for calling DIIS:
         if ( (ftot < DIIS_FORCE_THRESHOLD .and. (.not. ITERATIVE) ) .or. &
              (ftot < DIIS_FORCE_THRESHOLD .and. liter > 2 .and. ITERATIVE) .or. &
              (nsteps_after_eigen_min >= INFLECTION .and. ITERATIVE) ) then

            ! If  previous_forces(1,:) is .ne. 0.0d0 is a restart
            ! from a lanczos_step in the next IF statement. 

            if ( all( previous_forces(1,:) .eq. 0.0d0 ) ) then
               ! We set the first trial vector
               previous_forces(1,:) = force(:)
               previous_pos(1,:)    = pos(:)
               previous_norm(1)     = ftot*DIIS_STEP
               ! The second trial step is the result of another
               ! call lanczos_step; but without the calculation
               ! of a new projection vector
               get_proj = .False.
               call lanczos_step ( saddle_energy, a1, liter, get_proj )
               pas = pas + 1
               liter = liter + 1
            end if
            
            ! if after this last lanczos step the criteria to end the activation are fulfilled
            ! we end here.          
            if ( (ftot < EXITTHRESH) .or. (pas > MAXPAS) .or. (liter > 280) .or. & 
                & (eigenvalue > 0.0) .or. &
                & ( delta_e < delta_thr .and. delr < delr_thr .and. pas > (MAXKTER+5) )) then 

               saddle_energy = total_energy 
               end_activation = .true.  
               pas = pas - 1                ! the real number of steps in the event 

               if ( ftot < EXITTHRESH ) then
                  ret = 20000 + pas
               else if  (delta_e < delta_thr .and. delr < delr_thr .and. pas > (MAXKTER+5) ) then
                  ret = 80000 + pas 
               else if (eigenvalue > 0.0) then
                  ret = 60000 + pas
               else if ( pas > maxpas .or. liter > 280 ) then 
                  ret = 50000 + pas 
               end if 

               exit
            else 
               ! We set the second trial vector 
               previous_forces(2,:) = force(:)
               previous_pos(2,:)    = pos(:)
               previous_norm(2)     = ftot*DIIS_STEP
               switchDIIS = .True.   ! let's make DIIS.
               exit                  ! We go back
            end if
         end if    
      end if


      ! Hide option: if the test is true; we calculate the projection only
      ! at every two steps but after 4 steps above of an inflection in the eigenvalue
      if ( nsteps_after_eigen_min>=4 .and. mod(pas,2)==0 .and. a1>0.9d0 .and. calc_proj ) then 
         get_proj = .False.
      else
         get_proj = .True.
      end if

      call lanczos_step ( saddle_energy, a1, liter, get_proj ) 
                                      ! clean_wf
      if_clean_wf: if ( clean_wf ) then                                     
          if ( eigenvalue == eigen_min .and. ftot > EXITTHRESH) then
             call clean_wavefunction ( a1, .True. )
          else if ( ftot < EXITTHRESH ) then

             call clean_wavefunction ( a1, .True. )
             if ( ftot >= EXITTHRESH .and. cw_try_again ) then
                                           ! we continue in lanczos,
                cw_try_again = .False.     ! but, this is done only once per event.  
             else
                ret = 40000 + pas
                saddle_energy = total_energy
                end_activation = .true.
                exit
             end if
          end if 
      end if if_clean_wf

      pas = pas + 1
      liter = liter + 1
        
  end do While_lanczos 

END SUBROUTINE apply_lanczos


!> ART lanczos_step
!! We move from one hyperplane to the next one. We relax the position in it.
!! Finally we get the eigenvector at the relaxed configuration.
!! Modified by:
!! -Laurent Karim Beland, UdeM 2011: Collinear
subroutine lanczos_step ( current_energy, a1, liter, get_proj ) 

  use defs 
  use bigdft_forces
  use saddles
  use lanczos_defs
  use diis_defs
  implicit none

  !Arguments
  real(kind=8), intent(out) :: current_energy      ! Accept energy.
  real(kind=8), intent(out) :: a1                  ! dot product between two eigenvectors
  integer, intent(in)  :: liter 
  logical, intent(in)  :: get_proj  ! If we need get the projection vector.

  !Local variables
  integer :: i, step_rejected

  real(kind=8), dimension(VECSIZE) :: pos_b        ! Position for evaluation.
  real(kind=8), dimension(VECSIZE) :: force_b      ! Total Force for evaluation.
  real(kind=8), dimension(VECSIZE) :: perp_force   ! Perpendicular force...
  real(kind=8), dimension(VECSIZE) :: perp_force_b ! ...& for evaluation.  

  real(kind=8), dimension(3) :: boxl
  real(kind=8) :: step                       ! This is the step in the hyperplane. 
  real(kind=8) :: ftot_b                     ! ...& for evaluation.
  real(kind=8) :: fpar_b                     ! ...& for evaluation. 
  real(kind=8) :: fperp_b                    ! ...& for evaluation.
  real(kind=8) :: current_fperp              ! fperp as a criteria of minimization.

  logical      :: new_projection 
  !_______________________
  boxl = box * scala                  ! We compute at constant volume.

  ! We move the configuration along the eigendirection corresponding to the lowest
  ! eigenvalue. With sign(1.0d0,fpar) we orient the direction of the eigenvector 
  ! corresponding to the negative eigendirection (called projection) such that it points
  ! away from minimum. 

  fpar = dot_product( force, projection )
  if ( abs(fpar) > EXITTHRESH*0.5 ) then 
     pos = pos - sign(1.0d0,fpar) * 2.0d0 * INCREMENT * projection / sqrt(1.0d0*liter) 
  else
     pos = pos - sign(1.0d0,fpar) * 1.0d0 * INCREMENT * projection / sqrt(1.0d0*liter) 
  end if
                                      ! Current energy and force.
  call calcforce( NATOMS, pos, boxl, force, current_energy, evalf_number, .false. )
                                      ! New force's components.
  call force_projection( fpar, perp_force, fperp, ftot, force, projection )
  current_fperp = fperp

  step = 0.04d0*INCREMENT          
  try = 0
  m_perp = 0
  step_rejected = 0

  ! We relax perpendicularly using a simple variable-step steepest descent
  While_perpi: do                     

    pos_b = pos + step * perp_force

    call calcforce( NATOMS, pos_b, boxl, force_b, total_energy, evalf_number, .false. )
                                      ! New force's components.
    call force_projection( fpar_b, perp_force_b, fperp_b, ftot_b, &
                  & force_b, projection )      
    
    if ( fperp_b < current_fperp ) then
       pos            = pos_b
       current_energy = total_energy
       fpar           = fpar_b
       perp_force     = perp_force_b
       fperp          = fperp_b
       ftot           = ftot_b
       force          = force_b

       current_fperp  = fperp_b
       m_perp = m_perp + 1
       step = 1.6 * step
       step_rejected = 0
    else
       step = 0.8 * step
       step_rejected = step_rejected + 1
    end if
    try = try + 1 
                                      ! exit criteria 
    if ( fperp < FTHRESHOLD .or. m_perp >= MAXIPERP .or. &
       & (try > 12 .or. (liter < 2 .and. try > 2 ) ) ) exit While_perpi
    
  end do  While_perpi 

  delta_e = current_energy - ref_energy
                                      ! Magnitude of the displacement (utils.f90).
  call displacement( posref, pos, delr, npart )
  if ( SAVE_CONF_INT ) call save_intermediate( 'L' )

  ! we get eigen direction for the minimum of this hyperplane.
  if ( get_proj ) then
     do i = 1, LANCZOS_SCL            ! Lanczos call, we start from the
        new_projection = .false.      ! previous direction each time.
        call lanczos( NVECTOR_LANCZOS_C, new_projection, a1 )
        if ( iproc == 0 ) write(*,'(a,3I5,f12.6,f7.2)') &
           & 'BART COLLINEAR:', pas, liter, i, eigenvalue, a1
        if ( a1 > collinear_factor ) exit 
     end do
  end if
                                      ! Debug
  !if ( iproc == 0 ) then                
  ! write(*,*) 'BART: eigenvalue : ', eigenvalue
  ! write(*,"(' ','BART: eigenvals: ',4f12.6)") (eigenvals(i),i=1,4)
  !end if
                                      ! Write 
  call write_step ( 'L', liter, a1, current_energy )
                                      ! We look for an inflection in the
                                      ! eigenvalue.
  eigen_min = min( eigen_min, eigenvalue )
  if ( eigenvalue == eigen_min ) nsteps_after_eigen_min = 0
  if ( eigenvalue > eigen_min  ) nsteps_after_eigen_min = nsteps_after_eigen_min + 1
                                      ! For restart ( restart.f90 
  if ( write_restart_file ) then
     state_restart = 2 
     total_energy = current_energy
     if ( iproc == 0 ) then
        call save_state2(state_restart,liter+1,projection,maxter,previous_forces,&
        & previous_pos,previous_norm,eigen_min,eigenvalue,nsteps_after_eigen_min )
     end if
  end if

END SUBROUTINE lanczos_step


!> ART apply_glisse
!! @author
!! Written by Laurent Karim Beland, UdeM 2011!!
!! IN DEVELOPMENT
!! This routine converges to the saddle point by using
!! the inertial properties of the solid
subroutine apply_glisse( giter, saddle_energy )

  use defs
  use lanczos_defs
  use diis_defs
  use saddles
  implicit none

  !Arguments
  integer, intent(inout) :: giter
  real(kind=8), intent(out)   :: saddle_energy

  !Local variables
  real(kind=8) :: current_energy              ! Accept energy.
  real(kind=8) :: a1
  logical :: should_reset_velo,new_projection

  real(kind=8), parameter :: dt = 0.015d0
  real(kind=8), dimension(3*natoms) :: fcur,poscur,velcur,fpred,pospred,velpred,perp_force
  integer :: liter,i,step_rejected
  logical :: conv

  real(kind=8), dimension(VECSIZE) :: pos_b        ! Position for evaluation.
  real(kind=8), dimension(VECSIZE) :: force_b      ! Total Force for evaluation.
  real(kind=8), dimension(VECSIZE) :: perp_force_b ! ...& for evaluation.  

  real(kind=8) :: step                        ! This is the step in the hyperplane. 
  real(kind=8) :: ftot_b                     ! ...& for evaluation.
  real(kind=8) :: fpar_b                     ! ...& for evaluation. 
  real(kind=8) :: fperp_b                    ! ...& for evaluation.
  real(kind=8) :: current_fperp              ! fperp as a criteria of minimization.
  !_______________________ 

  conv = .false.

  if ( .not. restart ) then
     ! To be consistent with the restart file, each time we call
     ! this subroutine we make this initialization.
     previous_forces = 0.0d0
     previous_pos    = 0.0d0
     previous_norm   = 0.0d0
     maxter = 0
     eigen_min = 0.0d0
  else
     ! A warranty 
     eigen_min = min( eigen_min, eigenvalue )
     if ( eigenvalue == eigen_min ) then
        previous_forces = 0.0d0
        previous_pos    = 0.0d0
        previous_norm   = 0.0d0
     end if
     restart = .false.
  end if

  ! We first get the current force and energy
  call calcforce( NATOMS, pos, box, force, current_energy, evalf_number, .false. )
  fcur = force
  poscur = pos

  !We then get the eigendirection

  liter = 1
  do i = 1,3
      new_projection = .false.        ! previous direction each time.
      call lanczos( NVECTOR_LANCZOS_C, new_projection, a1 )
       delta_e = current_energy - ref_energy
                                      ! Magnitude of the displacement (utils.f90).
       call displacement( posref, pos, delr, npart )
                                      ! Write 
       call write_step ( 'L', liter, a1, current_energy )

     if (a1 > 0.70d0) exit
  enddo

  giter = 0
  should_reset_velo = .true.
  do
     giter = giter + 1
     pas = pas + 1

     !we then use this eigendirection to establish the initial velocity
     if (should_reset_velo) then
        call force_projection( fpar, perp_force, fperp, ftot, force, projection )

        current_fperp = fperp
  
        step = 0.04d0*INCREMENT
        try = 0
        m_perp = 0
        step_rejected = 0

        do                     
           pos_b = pos + step * perp_force

           call calcforce( NATOMS, pos_b, box, force_b, total_energy, evalf_number, .false. )
                                      ! New force's components. 
           call force_projection( fpar_b, perp_force_b, fperp_b, ftot_b, &
                  & force_b, projection )      

           if ( fperp_b < current_fperp ) then
              pos            = pos_b
              current_energy = total_energy
              fpar           = fpar_b
              perp_force     = perp_force_b
              fperp          = fperp_b
              ftot           = ftot_b
              force          = force_b

              current_fperp  = fperp_b
              m_perp = m_perp + 1
              step = 1.6 * step
              step_rejected = 0
           else
              step = 0.8 * step
              step_rejected = step_rejected + 1
           end if
           try = try + 1
                                      ! 5 kills MAXIPERP 
           if ( fperp < FTHRESHOLD .or. &
              &  m_perp > MAXIPERP+3 .or. try > 8 ) exit 

        end do 

        velcur =   -3.0d0*sign(1.0d0,fpar) * abs(ftot) * projection  ! all mass eq 1

        should_reset_velo = .false.
     endif

     !we now start the verlet velocity integrator
     pospred=pos+dt*velcur + dt*dt*0.5d0*fcur  ! all mass = 1
     call calcforce(natoms,pospred,box,fpred,total_energy,evalf_number,conv)
     velpred = velcur + 0.5d0*dt*fpred + 0.5d0*dt*fcur

     !we then make sure the norm of the velocity goes up
     !otherwise, we stop and recalculate lanczos

     if ( abs(dot_product(velpred,projection)) >abs(dot_product(velcur,projection))) then
        should_reset_velo = .true.
     endif

     if (should_reset_velo ) then
 !       liter = liter + 1
 !       pos = poscur
 !       call calcforce( NATOMS, pos, box, force, current_energy, evalf_number, .false. )

 !       do j = 1,3

 !          new_projection = .false.        ! previous direction each time.
 !          call lanczos( NVECTOR_LANCZOS_C, new_projection, a1 )
 !          delta_e = current_energy - ref_energy
                                      ! Magnitude of the displacement (utils.f90).
 !          call displacement( posref, pos, delr, npart )
                                      ! Write 
 !          call write_step ( 'L', liter, a1, current_energy )
 !          if (a1 > 0.70d0) exit
 !       enddo
 !       cycle
     endif

     !we now look if we should stop this activation
     ftot = dsqrt(dot_product(fpred,fpred))
     if ( (ftot < EXITTHRESH) .or. (pas > MAXPAS) .or. eigenvalue > 0.0d0 .or. giter > 80) then
        pos = pospred
        force = fpred

        saddle_energy = total_energy
        end_activation = .true.
        exit
     endif

     !then we look if we should start diis
     if (ftot < DIIS_FORCE_THRESHOLD) then !
        pos = pospred
        force = fpred

        if ( all( previous_forces(1,:) .eq. 0.0d0 ) ) then
               ! We set the first trial vector
               previous_forces(1,:) = force(:)
               previous_pos(1,:)    = pos(:)
               previous_norm(1)     = ftot*DIIS_STEP

             !  get_proj = .true.
             !  liter = 1
             !  call lanczos_step ( saddle_energy, a1, liter, get_proj)

            if ( (ftot < EXITTHRESH) .or. (pas > MAXPAS) &
               & .or. (eigenvalue > 0.0) ) then
               end_activation = .true.
               exit
            else
               ! We set the second trial vector 
              ! previous_forces(2,:) = force(:)
              ! previous_pos(2,:)    = pos(:)
              ! previous_norm(2)     = ftot*DIIS_STEP
               switchDIIS = .True.   ! let's make DIIS.
               exit                  ! We go back
            end if
        endif
     endif

     fcur=fpred
     poscur=pospred
     pos = poscur
     force = fpred

     call force_projection( fpar, perp_force, fperp, ftot, force, projection )

     velcur = 0.60d0*velpred + 0.4d0*perp_force / fperp* dsqrt(dot_product(velpred,velpred))

     a1 = 0.0d0
     delta_e = total_energy - ref_energy
                             ! Magnitude of the displacement (utils.f90).
     call displacement( posref, pos, delr, npart )
                             ! Write 
     call write_step ( 'G', giter, a1, current_energy )
  enddo

END SUBROUTINE apply_glisse
