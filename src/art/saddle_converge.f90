!!****f* art/saddle_converge
!! FUNCTION
!!   This subroutine bring the configuration to a saddle point. It does that
!!   by first pushing the configuration outside of the harmonic well, using
!!   the initial direction selected in find_saddle. Once outside the harmonic
!!   well, as defined by the appearance of a negative eigenvalue (or
!!   reasonnable size) the configuration follows the direction corresponding
!!   to this eigenvalue until the force components parallel and perpdendicular
!!   to the eigendirection become close to zero.
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
subroutine saddle_converge( ret, saddle_energy )

  use defs
  use saddles
  use lanczos_defs
  use bigdft_forces
  implicit None

  !Arguments
  integer, intent(out) :: ret
  real(8), intent(out) :: saddle_energy  ! Energy at saddle point.

  !Local variables
  logical :: new_projection              ! For lanczos.
                                         ! Loop indeces :
  integer :: i, kter, kter_init, iter, iter_init, m_perp, m_rejected, try
  integer :: ierror                      ! File control.
  
  integer :: npart                       ! number of particles having moved by more 
                                         ! than a THRESHOLD
  real(8) :: delr                        ! Distance between two configurations.
  real(8) :: step                        ! This is the true step. 
  real(8) :: boxl

  real(8) :: current_energy                   ! Accept energy.
  real(8), dimension(VECSIZE) :: pos_b        ! Position for evaluation.
  real(8), dimension(VECSIZE) :: force_b      ! Total Force for evaluation.

  real(8), dimension(VECSIZE) :: perp_force   ! Perpendicular force...
  real(8), dimension(VECSIZE) :: perp_force_b ! ...& for evaluation.  

  real(8) :: ftot                       ! Norm of the total force...
  real(8) :: ftot_b                     ! ...& for evaluation.
  real(8) :: fpar                       ! Parallel projection of force
  real(8) :: fpar_b                     ! ...& for evaluation. 
  real(8) :: fperp                      ! Norm of the perpendicular force
  real(8) :: fperp_b                    ! ...& for evaluation.

  real(8) :: current_fperp              ! fperp as a criteria of minimization. 
  real(8) :: delta_e                    ! current_energy - ref_energy
                                        ! If we lost the eigendirection :
  real(8), dimension(VECSIZE) :: pos_p  ! previous position
  real(8) :: fpar_p                     ! previous fpar.
  
  real(8) :: a1

                                      ! We compute at constant volume.
  boxl = box * scala

                                      ! If restarts in harmonic well.
  if ( restart .and. ( state_restart == 1 ) ) then 

     initial_direction = direction_restart

     kter_init = iter_restart         ! kter loop, init value.

                                      ! Write
     if ( iproc == 0 ) Then 
      write(*,*) 'BART: Restart'
      write(*,*) 'BART: in harmonic well '
      write(*,*) 'BART: kter : ', kter_init
      write(*,*) 'BART: pos: ', pos(1), pos(2), pos(3)
     end if
     restart = .false.
  else
                                      ! kter loop, init value.
    kter_init = 0                   
  end if

                                      ! Write header in log.file.
  if ( iproc == 0 ) then
   open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
       & action = 'write', position = 'append', iostat = ierror )
   write(FLOG, '(a22,a8,a8,a12,a12,a12,a11,a7,a6,a5)')&
   &'E-Eref','m_perp','ftot','fpar','fperp','eigen','delr','npart','evalf','a1'
   write(FLOG, '(A22,A31,A27)' ) '( eV )', '( eV/Ang )', '( eV/Ang**2 )'
   close(FLOG)
  end if
! _________

  If_restart: if ( ( .not. restart ) .and. NEW_EVENT ) then 

     eigenvalue = 0.0d0
     a1         = 0.0d0

     call calcforce( NATOMS, pos, boxl, force, current_energy )
     evalf_number = evalf_number + 1

    ! We now project out the direction of the initial displacement from the
    ! minimum from the force vector so that we can minimize the energy in
    ! the direction perpendicular to the direction of escape from harmonic well

     call force_prj( fpar, perp_force, fperp, ftot, force, initial_direction )

                                      ! The true step is not the increment.
     step = 0.4d0*INCREMENT           ! ?? HERE OR INSIDE Do_kter ??
    
     new_projection = .true.          ! COSMIN, but is it necessary???

     Do_kter: do kter = kter_init, MAXKTER
      
        try = 0                       ! # iterationes of While loop.
        m_perp = 0                    ! Accepted perpendicular iterations.
        m_rejected = 0                ! Rejected movements.

                                      ! We relax perpendicularly using a simple
                                      ! variable-step steepest descent 
        While_perpk: do 

          pos_b = pos + step * perp_force

          call calcforce( NATOMS, pos_b, boxl, force_b, total_energy )
          evalf_number = evalf_number + 1

                                      ! Force's components.
          call force_prj( fpar_b, perp_force_b, fperp_b, ftot_b, &
                       & force_b, initial_direction ) 

          if ( total_energy < current_energy ) then
             pos            = pos_b
             current_energy = total_energy
             fpar           = fpar_b
             perp_force     = perp_force_b
             fperp          = fperp_b
             ftot           = ftot_b
             force          = force_b

             m_perp = m_perp + 1
             step = 1.2 * step
             m_rejected = 0
          else
             step = 0.6 * step
             m_rejected = m_rejected + 1
          end if
          try = try + 1
         
          if ( fperp < FTHRESHOLD .or. m_perp > MAXKPERP &
              &  .or. m_rejected > 5 ) exit While_perpk

        end do While_perpk

                                      ! We start checking of negative eigenvalues 
        if ( kter >= KTER_MIN ) then  ! eigenvalues only after a few steps.

                                      ! We do not use previously computed 
           new_projection = .true.    ! lowest direction.

                                      ! In lanczos.f90
           call lanczos( NVECTOR_LANCZOS, new_projection, a1 )
           new_projection = .false.   ! COSMIN

        end if

        current_energy = total_energy ! As computed in lanczos routine.
        delta_e = current_energy - ref_energy

                                      ! Magnitude of the displacement (utils.f90).
        call displacement( posref, pos, delr, npart )

        if ( SAVE_CONF_INT ) call save_intermediate( 'K', kter ) 

                                      ! Write 
        if ( iproc == 0 ) then
         open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
             & action = 'write', position = 'append', iostat = ierror )
         write( FLOG, '(a6,i4,x,F11.4,i3,i2,4f12.4,x,1f10.3,i4,i6,f7.2)') &
         & "kter=", kter, delta_e, m_perp, try, ftot, fpar, fperp,   & 
         &  eigenvalue, delr, npart, evalf_number, a1
         close( FLOG )
         write(*,'(a,i4,x,(1p,e17.10,0p),x,2i3,4f12.6,x,1f10.4,i4,i6)') &
         & " BART: kter: '", kter, current_energy, m_perp, try, ftot,   & 
         & fpar, fperp, eigenvalue, delr, npart, evalf_number
        end if

                                      ! For restart ( restart.f90 )
        state_restart = 1
        if ( iproc == 0 ) call save_state( state_restart, kter+1, initial_direction )

                                      ! is the configuration out of the harmonic basin?
        if ( eigenvalue < EIGEN_THRESH ) then 

           exit Do_kter
        else

                                      ! We move the configuration along the initial direction.
           pos = pos + BASIN_FACTOR * INCREMENT * initial_direction  

                                      ! Current energy and force.
           call calcforce( NATOMS, pos, boxl, force, current_energy )
           evalf_number = evalf_number + 1

                                      ! New force's components.
           call force_prj( fpar, perp_force, fperp, ftot, force, initial_direction )

        end if

     end do Do_kter
     
     ! The configuration is now out of the harmonic well, we can now bring
     ! it to the saddle point. Again, we split the move into a parallel and 
     ! perpendicular contribution.
      
     ! First, we must now orient the direction of the eigenvector corresponding to the
     ! negative eigendirection (called projection) such that it points away from minimum. 
     
     fpar = dot_product( force, projection )
     if ( fpar > 0.0d0 ) projection = -1.0d0 * projection

                                      ! We move the configuration along the projection.
     pos = pos +  1.0d0 * INCREMENT * projection 

                                      ! Current and energy and force.
     call calcforce( NATOMS, pos, boxl, force, current_energy )
     evalf_number = evalf_number + 1 

                                      ! New force's components.
     call force_prj( fpar, perp_force, fperp, ftot, force, projection )
     current_fperp = fperp 

                                      ! if we lost the eigendirection
     pos_p  = pos   
     fpar_p = fpar

     iter_init = 1
                                      ! Else If_restart
  else if ( restart .and. ( state_restart == 2 ) ) then 

                                      ! This is the restart from a previous ACTIVATION
                                      ! ACTIVATION process. 
     iter_init  = iter_restart
     projection = direction_restart

     call calcforce( NATOMS, pos, boxl, force, current_energy )
     evalf_number = evalf_number + 1

                                      ! Force's components.
     Call force_prj( fpar, perp_force, fperp, ftot, force, projection )
     current_fperp = fperp

                                      ! If we lost the eigendirection
     pos_p  = pos   
     fpar_p = fpar
      
     if ( iproc == 0 ) write(*,*) 'BART: Restart = 2'
     restart = .false.

  else if ( .not. new_event ) then    ! Else If_restart

                                      ! This is a convergence event, we must
                                      ! compute the eigenvector with precision.

     posref = pos                     ! ??? FOR WHAT ??

     new_projection = .true. 
     do i = 1, 5
        call lanczos( NVECTOR_LANCZOS, new_projection, a1 )
        new_projection = .false.
     end do

     call calcforce( NATOMS, pos, boxl, force, current_energy )
     evalf_number = evalf_number + 1

     fpar = dot_product(force,projection)
     if( fpar > 0.0d0 ) projection = -1.0d0 * projection

                                      ! Force's components.
     call force_prj( fpar, perp_force, fperp, ftot, force, projection )
     current_fperp = fperp

                                      ! If we lost the eigendirection.
     pos_p = pos   
     fpar_p = fpar
                                      ! For avoiding first 100 steps in next
     iter_init = 100                  ! loop ???; ummm.

  else                                ! Else If_restart

     write(*,*) 'BART: Problem with restart and state_restart : '
     write(*,*) 'BART: restart = ', restart
     write(*,*) 'BART: state_restart = ', state_restart
     stop

  end if If_restart
! _________ 
                                      ! ACTIVATION.
                                      ! If activation did not converged before
                                      ! MAXITER then ret = 0, and the the saddle
                                      ! is not accepted. 
  ret = 0 

  Do_iter: do iter = iter_init, MAXITER 

     step = 0.04d0*INCREMENT          ! ?? HERE OR OUTSIDE Do_iter ??
     reject = .false.                 ! If we lost the eigendirection

     try = 0
     m_perp = 0
     m_rejected = 0

                                      ! The next loop is on the perpendicular direction,
                                      ! again using a variable-step steepest descent
     While_perpi: do 

       pos_b = pos + step * perp_force

       call calcforce( NATOMS, pos_b, boxl, force_b, total_energy )
       evalf_number = evalf_number + 1

                                      ! New force's components.
       call force_prj( fpar_b, perp_force_b, fperp_b, ftot_b, &
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
          m_rejected = 0
       else
          step = 0.8 * step
          m_rejected = m_rejected + 1
       end if
       try = try + 1 

                                      ! 5 kills MAXIPERP 
       if ( fperp < FTHRESHOLD .or. m_perp > (iter - 10) &         
          & .or. m_perp > MAXIPERP .or. try > 5 ) exit While_perpi
       
     end do  While_perpi   

                                      ! Lanczos call, we start from the
     new_projection = .false.         ! previous direction each time.
     call lanczos( NVECTOR_LANCZOS, new_projection, a1 )

                                      ! Write     
     if ( iproc == 0 ) then 
      write(*,*) 'BART: eigenvalue : ', eigenvalue
      write(*,"(' ','BART: eigenvals: ',4f12.6)") (eigenvals(i),i=1,4)
     end if
    
     current_energy = total_energy    ! As computed in lanczos routine
     saddle_energy  = current_energy
     delta_e = current_energy - ref_energy

                                      ! Magnitude of the displacement (utils.f90).
     call displacement( posref, pos, delr, npart )

     if ( SAVE_CONF_INT ) call save_intermediate( 'I', iter )

                                      ! Write
     if ( iproc == 0 ) then
      open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
          & action = 'write', position = 'append', iostat = ierror )
      write(FLOG,'(a6,i4,x,F11.4,i3,i2,4f12.4,x,1f10.3,i4,i6,f7.2)') &
      & 'iter=', iter, delta_e, m_perp, try, ftot, fpar, fperp,   &
      &  eigenvalue, delr, npart, evalf_number, a1
      close(FLOG) 
      write(*,'(a,i4,x,(1p,e17.10,0p),i3,i2,4f12.6,x,1f10.4,i4,i6)') &
      & " BART: iter: ", iter, current_energy, m_perp, try, ftot,    &
      & fpar, fperp, eigenvalue, delr, npart, evalf_number
     end if

                                      ! If force is below some threshold
                                      ! of if we have a positive eigenvalue.
     If_ftot_eigen: if ( ftot < EXITTHRESH ) then

        If_diis1: if ( USE_DIIS ) then 

                                      ! Write 
           if ( iproc == 0 ) then
            open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
                 & action = 'write', position = 'append', iostat = ierror )
            write(FLOG, * ) "Call DIIS"
            close(FLOG)
           end if

                                       ! In diis.f90 
           call apply_diis( current_energy, ftot, delr, npart )
           saddle_energy = current_energy
           fpar  = 0.0d0
           fperp = 0.0d0

           If_check: if ( DIIS_CHECK_EIGENVEC ) then

                                      ! Lanczos several times.
              new_projection = .true.
              Do_lanc: do i = 1, 4

                 call lanczos( NVECTOR_LANCZOS, new_projection, a1 )
                 new_projection = .false.

                                      ! Write
                 if ( iproc == 0 ) then 
                  write(*,*) 'BART: Iter ', 1, ' : ', total_energy,  eigenvalue  
                 end if
                                      ! Exit of loop.  
                 if ( eigenvalue < 0.0d0 ) exit Do_lanc
              end do Do_lanc

              if ( eigenvalue >= 0.0 .or. ftot > DIIS_FORCE_THRESHOLD ) then
                                      ! failed.
                 ret = 60000 + iter
              else
                                      ! As computed in lanczos routine.
                 saddle_energy = total_energy 

                                      ! New fpar and fperp. 
                 call force_prj( fpar, perp_force, fperp, ftot, force, projection )
                 ret = 20000 + iter

              end if

           else                       ! Else of If_check 

                                      ! Write 
              if ( iproc == 0 ) then 
               open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
                   & action = 'write', position = 'append', iostat = ierror )
               write(FLOG,*) 'Warning : Not checking of saddle point' 
               close(FLOG)            
              end if
             
              eigenvalue = 0.0d0
              ret = 20000 + iter
              if ( ftot > DIIS_FORCE_THRESHOLD ) ret = 70000 + iter

           end if If_check

        else                          ! Else of If_diis1

           ret = 20000 + iter

        end if If_diis1

        exit Do_iter

                                      ! Else of If_ftot_eigen
     else if ( (abs(fpar) < 0.1*EXITTHRESH ) .and. ( fperp < EXITTHRESH ) ) then 

        If_diis2: if ( USE_DIIS ) then

                                      ! Write 
           if ( iproc == 0 ) then 
            open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
                & action = 'write', position = 'append', iostat = ierror )
            write(FLOG, * ) "Call DIIS (2)"
            close(FLOG)
           end if

                                      ! diis in diis.f90 
           call apply_diis( current_energy, ftot, delr, npart )
           saddle_energy = current_energy
           fpar  = 0.0d0
           fperp = 0.0d0

           If_check2: if ( DIIS_CHECK_EIGENVEC ) then

                                      ! Lanczos several times.
              new_projection = .true.
              Do_lanc2: do i = 1, 4

                 call lanczos( NVECTOR_LANCZOS, new_projection, a1 )
                 new_projection = .false.

                                      ! Write
                 if ( iproc == 0 ) then 
                  write(*,*) 'BART: Iter ', 1, ' : ', total_energy,  eigenvalue  
                 end if
                                      ! Exit of loop.
                 if ( eigenvalue < 0.0d0 ) exit  Do_lanc2
              end do Do_lanc2

              if ( eigenvalue >= 0.0 .or. ftot > DIIS_FORCE_THRESHOLD ) then
                                      ! failed.
                 ret = 60000 + iter
              else
                                      ! As computed in lanczos routine.
                 saddle_energy = total_energy 

                                      ! New fpar and fperp.  
                 call force_prj( fpar, perp_force, fperp, ftot, force, projection ) 
                 ret = 10000 + iter

              end if

           else                       ! Else of If_check2 

                                      ! Write 
              if ( iproc == 0 ) then 
               open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
                   & action = 'write', position = 'append', iostat = ierror )
               write(FLOG,*) 'Warning : Not checking of saddle point' 
               close(FLOG)
              end if

              eigenvalue = 0.0d0
              ret = 10000 + iter
              if ( ftot > DIIS_FORCE_THRESHOLD ) ret = 70000 + iter

           end if If_check2

        else                          ! Else of If_diis2 

           ret = 0000 + iter

        end if If_diis2    

        exit  Do_iter 

     else if ( eigenvalue > 0.0 ) then! Else of If_ftot_eigen
                                      ! Failed !!
        ret = 60000 + iter
        exit  Do_iter
           
     end if If_ftot_eigen

                                      ! If we lost the eigendirection ! 
     pos_p = pos
     fpar_p = fpar  
                                      ! We now move the configuration along the
                                      ! eigendirection corresponding to the lowest
                                      ! eigenvalue.
     if ( abs(fpar) > EXITTHRESH*0.5 ) then 
        pos = pos - sign(1.0d0,fpar) * 2.0d0 * INCREMENT * projection / sqrt(1.0d0*iter) 
     else
        pos = pos - sign(1.0d0,fpar) * 1.0d0 * INCREMENT * projection / sqrt(1.0d0*iter) 
     end if

                                      ! COSMIN 
     call calcforce( NATOMS, pos, boxl, force, current_energy ) 
     evalf_number = evalf_number + 1

                                      ! New force's components.
     call force_prj( fpar, perp_force, fperp, ftot, force, projection )
     current_fperp = fperp

                                      ! For restart ( restart.f90 )
     state_restart = 2
     if ( iproc == 0 ) Call save_state( state_restart, iter+1, projection )
    
  end do Do_iter
! _________ 

  if ( ( ret > 0 ) .and. ( ret < 30000 ) ) then 

     delta_e = saddle_energy - ref_energy

                                      ! Write
     if ( iproc == 0 ) then 
      open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
          & action = 'write', position = 'append', iostat = ierror )
      write(FlOG,*) ' '
      write(FLOG, '(a7,i5,F10.4,i3,i2,4f12.4,x,1f10.3,i4,i6,f7.2)') &
      & 'SADDLE', mincounter, delta_e, m_perp, try, ftot, fpar, fperp,   &
      &  eigenvalue, delr, npart, evalf_number, a1
      write(FLOG,*) ' '
      close(FLOG) 

      write( *, '(A,i5,F10.4,i3,i2,4f12.4,x,1f10.3,i4,i6,f7.2)') &
      & 'BART: SADDLE', mincounter, delta_e, m_perp, try, ftot, fpar, fperp,   &
      &  eigenvalue, delr, npart, evalf_number, a1
     end if

  end if

END SUBROUTINE saddle_converge
!!***


!!****f* saddle_converge/force_prj
!! FUNCTION
!!   It calculates: 
!!    the magnitude the force in the direction of some 'direction' (F_par),
!!    the force vector perpendicular to that 'direction' (F_perp_V) and its
!!    magnitude(F_perp), and the norm of the total force (norm_F).
!!
!! SOURCE
!! 
subroutine force_prj ( F_par, F_perp_V, F_perp, norm_F, ref_F, direction )

  use defs , only : VECSIZE
  implicit none

  ! Arguments
  real(8),                     intent(out) :: F_par     ! Norm of parallel force.
  real(8), dimension(VECSIZE), intent(out) :: F_perp_V  ! Perpendicular force Vector. 
  real(8),                     intent(out) :: F_perp    ! Norm of F_perp_V.
  real(8),                     intent(out) :: norm_F    ! Norm of the input force.
  real(8), dimension(VECSIZE), intent(in)  :: ref_F     ! Input force.
  real(8), dimension(VECSIZE), intent(in)  :: direction ! Eigen direction.

  ! Internal variables 
  real(8) :: F_perp2                  ! F_perp*F_perp.


  F_par  = dot_product( ref_F, direction )
  F_perp_V  = ref_F - F_par * direction  

  F_perp2 = dot_product( F_perp_V , F_perp_V )  
  F_perp  = dsqrt( F_perp2 )
                                      ! This is = dsqrt( dot_product( ref_F, ref_F ) ) 
                                      ! i.e, the norm of the force, but cheaper.
  norm_F = dsqrt( F_par*F_par + F_perp2 ) 

END SUBROUTINE force_prj
!!***
