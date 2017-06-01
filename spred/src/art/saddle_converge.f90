!> @file
!! Routines used by ART
!! @author
!!    Copyright (C) Normand Mousseau, June 2001
!!    Copyright (C) 2010-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> This subroutine brings the configuration to a saddle point. It does that
!! by first pushing the configuration outside of the harmonic well, using
!! the initial direction selected in find_saddle. Once outside the harmonic
!! well, as defined by the appearance of a negative eigenvalue (or
!! reasonnable size) the configuration follows the direction corresponding
!! to this eigenvalue until the force components parallel and perpendicular
!! to the eigendirection become close to zero.
subroutine saddle_converge( ret, saddle_energy )

   use defs
   use saddles
   use lanczos_defs
   use bigdft_forces
   use diis_defs
   use module_base
   implicit none

   !Arguments
   integer, intent(out) :: ret
   real(kind=8), intent(out) :: saddle_energy  ! Energy at saddle point.

   !Local variables
   logical :: new_projection              ! For lanczos.
   ! Loop indeces :
   integer :: i, kter, kter_init, liter, diter, step_rejected
   integer :: ierror, ierr                ! File and MPI control.

   real(kind=8) :: step                        ! This is the step in the hyperplane. 
   real(kind=8),dimension(3) :: boxl
   real(kind=8) :: a1
   real(kind=8) :: current_energy              ! Accepted energy.
   real(kind=8) :: ftot_b                      ! ftot  for evaluation.
   real(kind=8) :: fpar_b                      ! fpar  for evaluation. 
   real(kind=8) :: fperp_b                     ! fperp for evaluation.
   real(kind=8), dimension(VECSIZE) :: pos_b        ! Position for evaluation.
   real(kind=8), dimension(VECSIZE) :: force_b      ! Total Force for evaluation.
   real(kind=8), dimension(VECSIZE) :: perp_force   ! Perpendicular force...
   real(kind=8), dimension(VECSIZE) :: perp_force_b ! ...& for evaluation.  
   ! __________________
   ! Write header in log.file.
   if ( iproc == 0 ) then
      open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
         & action = 'write', position = 'append', iostat = ierror )
      write(FLOG, '(a23,a8,a8,a12,a12,a12,a11,a7,a6,a5)')&
         &   'E-Eref','m_perp','ftot','fpar','fperp','eigen','delr','npart','evalf','a1'
      write(FLOG, '(a23,a31,a27)' ) '( eV )', '( eV/Ang )', '( eV/Ang**2 )'
      close(FLOG)
   end if
   ! initialization
   saddle_energy = 0.0d0
   boxl = box * scala                  ! We compute at constant volume.
   cw_try_again = .True.               ! For clean_wf

   if ( .not. restart ) pas = 0
   ! If restarts in harmonic well.
   if ( restart .and. ( state_restart == 1 ) ) then 

      kter_init = iter_restart         ! init value of kter loop.
      initial_direction = direction_restart
      call displacement( posref, pos, delr, npart )
      deallocate(direction_restart)

      restart = .false.
      ! Write
      if ( iproc == 0 ) then
         write(*,*) 'BART: Restart  in harmonic well '
         write(*,*) 'BART: kter : ', kter_init
         write(*,*) 'BART: pos: ', pos(1), pos(2), pos(3)
      end if
   else
      kter_init = 0                    ! init value of kter loop.
      ! clean_wf
      if ( clean_wf ) then
         call displacement( posref, pos, delr, npart )
         call clean_wavefunction ( 0.0d0, .False. ) 
      end if 
   end if
   ! _________
   !                HARMONIC WELL 

   If_restart: if ( ( .not. restart ) .and. NEW_EVENT ) then 

      eigenvalue = 0.0d0               ! Only for the report.
      a1         = 0.0d0               ! Only for the report. 
      step = 0.4*INCREMENT             ! The step in the hyperplane.
      new_projection = .true.          ! We do not use previously computed 
      ! lowest direction.
      Do_kter: do kter = kter_init, MAXKTER
         ! Reference energy and force.
         call calcforce( NATOMS, pos, boxl, force, current_energy, evalf_number, .false. )

         ! We now project out the direction of the initial displacement from the
         ! minimum from the force vector so that we can minimize the energy in
         ! the direction perpendicular to the direction of escape from harmonic well

         call force_projection( fpar, perp_force, fperp, ftot, force, initial_direction )

         try = 0                       ! Total # of iterationes in a given hyperplane. 
         m_perp = 0                    ! Perpendicular iterations accepted.
         step_rejected = 0             ! Perpendicular steps rejected.

         ! We relax perpendicularly using a simple variable-step steepest descent
         While_perpk: do 

            pos_b = pos + step * perp_force

            call calcforce( NATOMS, pos_b, boxl, force_b, total_energy, evalf_number, .false. )
            ! New force's components.
            call force_projection( fpar_b, perp_force_b, fperp_b, ftot_b, &
               &   force_b, initial_direction ) 

            if ( total_energy < current_energy ) then
               pos            = pos_b
               current_energy = total_energy
               fpar           = fpar_b
               perp_force     = perp_force_b
               fperp          = fperp_b
               ftot           = ftot_b
               force          = force_b

               step = 1.2 * step
               m_perp = m_perp + 1
               step_rejected = 0
            else
               step = 0.6 * step
               step_rejected = step_rejected + 1
            end if
            try = try + 1

            if ( fperp < FTHRESHOLD .or. m_perp >= MAXKPERP &
               &   .or. step_rejected > 5 ) exit While_perpk

         end do While_perpk

         delta_e = current_energy - ref_energy
         ! Magnitude of the displacement (utils.f90).
         call displacement( posref, pos, delr, npart )

         if ( SAVE_CONF_INT ) call save_intermediate( 'K' ) 

         ! We start checking of negative eigenvalues only after a few steps.

         if ( kter == KTER_MIN ) then 
            ! First time, twice !!
            do i = 1, 1            
               call lanczos( NVECTOR_LANCZOS_H, new_projection, a1 )
               new_projection = .false. 
               if ( iproc == 0 ) write(*,'(a,3I5,f12.6,f7.2)') &
                  &   'BART COLLINEAR:', pas, kter, i, eigenvalue, a1
            end do

         else if ( setup_initial .and. kter > KTER_MIN + 1 ) then 

            call check_min( 'I' ) 
            call write_step ( 'K', kter, a1, current_energy )
            call MPI_Barrier( MPI_COMM_WORLD, ierr )
            call end_art( )

         else if ( kter > KTER_MIN  ) then
            ! we get eigen direction for the minimum of this hyperplane.
            call lanczos( NVECTOR_LANCZOS_H, new_projection, a1 )
            ! Lanczos call, we start from the
            new_projection = .false.      ! previous direction each time.
         end if
         ! Write 
         call write_step ( 'K', kter, a1, current_energy )
         ! For restart ( restart.f90 )
         if ( write_restart_file ) then 
            state_restart = 1
            total_energy = current_energy
            if ( iproc == 0 ) call save_state( state_restart, kter+1, initial_direction )
         end if
         pas = pas + 1
         ! clean_wf
         if (eigenvalue<EIGEN_THRESH .and. clean_wf) call clean_wavefunction(a1,.True.)
         ! Is the configuration out of the harmonic basin?
         if ( eigenvalue < EIGEN_THRESH .and. (.not. setup_initial) ) exit Do_kter 
         ! If not, we move the configuration along 
         ! the initial direction.
         pos = pos + BASIN_FACTOR * INCREMENT * initial_direction  

      end do Do_kter
      ! If after MAXKTER iterations we have not found
      ! an inflection the event is killed.
      if ( eigenvalue > EIGEN_THRESH ) then
         saddle_energy = current_energy ! For the report.
         ret = 90000 + kter 
         return
      end if

      ! The configuration is now out of the harmonic well, we can now bring
      ! it to the saddle point. 
      ! First, we must now orient the direction of the eigenvector corresponding to the
      ! negative eigendirection (called projection) such that it points away from minimum.

      fpar = dot_product( force, projection )
      if ( fpar > 0.0d0 ) projection = -1.0d0 * projection

      liter = 1                   ! init value of lanczos loop.
      diter = 1                   ! init value of diis loop. 

      switchDIIS= .False.
      ! _________
   else if ( restart .and. &
      &   ( state_restart == 2 .or. state_restart == 4 ) ) then ! Else If_restart

   !           RESTART FROM A PREVIOUS ACTIVATION PROCESS 

   call allocate_activation ()
   projection      = direction_restart    
   previous_forces = diis_forces_restart  
   previous_pos    = diis_pos_restart    
   previous_norm   = diis_norm_restart  
   maxter     = maxter_r                  
   eigen_min  = eigen_min_r               
   eigenvalue = eigenvalue_r              
   nsteps_after_eigen_min = nsteps_after_eigen_min_r 
   delta_e = total_energy - ref_energy

   deallocate(direction_restart) 
   deallocate(diis_forces_restart)
   deallocate(diis_pos_restart)
   deallocate(diis_norm_restart) 

   if ( state_restart == 2 ) then 
      liter = iter_restart 
      diter = 1
      switchDIIS= .False.           ! We go to apply_lanczos.
      if ( iproc == 0 ) write(*,*) "BART: Restart = 2"
      elseif ( state_restart == 4 ) then 
      diter = iter_restart 
      liter = 1
      switchDIIS= .True.            ! We go to apply_diis.
      if ( iproc == 0 ) write(*,*) "BART: Restart = 4"
   else                             !DEBUG
      if ( iproc == 0 ) write(*,*) "BART: HOUSTON, we've got a problem"
      call end_art () 
   end if 

   call displacement( posref, pos, delr, npart )
   if (iproc==0) write(*,*) "BART: delr npart", delr, npart
   call force_projection( fpar, perp_force, fperp, ftot, force, projection )
   ! _________
else if ( .not. new_event ) then    ! Else If_restart

   !                'REFINE' EVENT 
   ! we must compute the eigenvector
   ! with precision.
   new_projection = .true. 

   do i = 1, 5
      call lanczos( NVECTOR_LANCZOS_H, new_projection, a1 )
      new_projection = .false.
      if ( iproc == 0 ) write(*,'(a,2I5,f12.6,f7.2)') &
         &   'BART COLLINEAR:', pas, i, eigenvalue, a1
      if ( a1 > collinear_factor ) exit
   end do

   delta_e = total_energy - ref_energy
   call displacement( posref, pos, delr, npart )
   fpar = dot_product(force,projection)
   if( fpar > 0.0d0 ) projection = -1.0d0 * projection
   call force_projection( fpar, perp_force, fperp, ftot, force, projection )
   ! This will be just for divide the INCREMENT,
   liter = 10                       ! i.e, a small move.
   ! Write 
   call write_step ( 'L', liter, a1, total_energy )

   liter = liter + 1
   diter = 1                        ! init value of diis loop. 
   pas   = pas + 1
   switchDIIS= .False.
   ! _________
else                                ! Else If_restart
   write(*,*) 'BART: Problem with restart and state_restart : '
   write(*,*) 'BART: restart = ', restart, ' state_restart = ', state_restart
   stop
end if If_restart
! _________ 
!                ACTIVATION PART

if ( .not. restart ) call allocate_activation ()

ret = 0
end_activation = .false.

While_activation: do

   if ( .not. switchDIIS ) then 
      call apply_lanczos( liter, saddle_energy, ret )
   else 
      call apply_diis( diter, saddle_energy, ret )
      if ( diter > 1 ) liter = 1 
      diter = 1
   end if

   if ( end_activation ) exit 
end do While_activation

call deallocate_activation ()

END SUBROUTINE saddle_converge


!> ART force_projection
!!   It calculates: 
!!    the magnitude the force in the direction of some 'direction' (F_par),
!!    the force vector perpendicular to that 'direction' (F_perp_V) and its
!!    magnitude(F_perp), and the norm of the total force (norm_F).
subroutine force_projection ( F_par, F_perp_V, F_perp, norm_F, ref_F, direction )

   use defs , only : VECSIZE
   implicit none

   ! Arguments
   real(kind=8),                     intent(out) :: F_par     ! Norm of parallel force.
   real(kind=8), dimension(VECSIZE), intent(out) :: F_perp_V  ! Perpendicular force Vector. 
   real(kind=8),                     intent(out) :: F_perp    ! Norm of F_perp_V.
   real(kind=8),                     intent(out) :: norm_F    ! Norm of the input force.
   real(kind=8), dimension(VECSIZE), intent(in)  :: ref_F     ! Input force.
   real(kind=8), dimension(VECSIZE), intent(in)  :: direction ! Eigen direction.

   ! Internal variables 
   real(kind=8) :: F_perp2                  ! F_perp*F_perp.

   F_par  = dot_product( ref_F, direction )
   F_perp_V  = ref_F - F_par * direction  

   F_perp2 = dot_product( F_perp_V , F_perp_V )  
   F_perp  = sqrt( F_perp2 )
   ! This is = sqrt( dot_product( ref_F, ref_F ) ) 
   ! i.e, the norm of the force, but cheaper.
   norm_F = sqrt( F_par*F_par + F_perp2 ) 

END SUBROUTINE force_projection


!> ART write_step
subroutine write_step ( stage, it, a1, energy )

   use defs
   use saddles
   use lanczos_defs
   implicit None

   !Arguments
   character(len=1), intent(in) :: stage
   integer,          intent(in) :: it
   real(kind=8),     intent(in) :: a1
   real(kind=8),     intent(in) :: energy         

   !Local variables
   integer :: ierror
   character(len=2)  :: etape  

   etape = stage//"="

   if ( iproc == 0 ) then
      open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
         & action = 'write', position = 'append', iostat = ierror )

      write( FLOG,'(i4,2x,a2,i4,1x,F10.4,i3,i2,4f12.4,1x,1f10.3,i4,i6,f7.2)')  &
         &   pas, etape, it, delta_e, m_perp, try, ftot, fpar, fperp,  & 
         &  eigenvalue, delr, npart, evalf_number, a1
      close( FLOG )

      write(*,'(a,i4,2x,a2,i4,1x,(1p,e17.10,0p),1x,2i3,4f12.6,1x,1f10.4,i4,i6)') &
         &   " BART:", pas, etape, it, energy,  m_perp, try, ftot, fpar, fperp,   &
         &   eigenvalue, delr, npart, evalf_number
   end if

END SUBROUTINE write_step


!> ART end_report
subroutine end_report ( success, ret, saddle_energy )

   use defs
   use saddles
   use lanczos_defs
   implicit None

   !Arguments
   logical, intent(out)   :: success
   integer, intent(inout) :: ret
   real(kind=8), intent(in) :: saddle_energy 

   !Local variables
   integer :: i
   integer :: ierror                   ! File control.
   real(kind=8) :: a1
   real(kind=8), dimension(VECSIZE) :: perp_force 
   logical :: new_projection           ! For lanczos.
   character(len=10) :: converg        ! For report.
   character(len=4)  :: scounter
   character(len=20) :: fname
   ! __________________

   select case( ret ) 

   case( 20000 : 29999 )                ! ftot < EXITTHRESH  

      If_diis: if ( USE_DIIS ) then
         fpar  = 0.0d0
         fperp = 0.0d0

         If_check: if ( DIIS_CHECK_EIGENVEC ) then
            ! Lanczos several times.
            new_projection = .false.
            Do_lanc: do i = 1, 4
               call lanczos( NVECTOR_LANCZOS_C, new_projection, a1 )
               new_projection = .false.
               ! Exit of loop.  
               if ( eigenvalue < 0.0d0 ) exit Do_lanc
            end do Do_lanc

            if ( eigenvalue >= 0.0d0 ) then
               ret = 60000 + pas 
            else                            ! Else of eigenvalue.
               ! New fpar and fperp. 
               call force_projection( fpar, perp_force, fperp,&
                  &   ftot, force, projection )
               ret = 10000 + pas
            end if

         else                               ! Else of If_check 
            eigenvalue = 0.0d0
            ret = 10000 + pas 
         end if If_check

      else                                  ! Pure Lanczos 
         ! Else of If_diis
         if ( eigenvalue > 0.0d0 ) then
            ret = 60000 + pas
         else 
            ret = 20000 + pas
         end if
      end if If_diis

   case( 40000 : 49999 )                ! clean_wf
      ! We accept the event is the new total
      ! force changes only in max 0.1 eV
      if ((ftot -0.1d0) < EXITTHRESH) then
         ! we check always the eigenvalue
         new_projection = .false.
         Do_lanc_c: do i = 1, 4
            call lanczos( NVECTOR_LANCZOS_C, new_projection, a1 )
            new_projection = .false.
            ! Exit of loop.  
            if ( eigenvalue < 0.0d0 ) exit Do_lanc_c
         end do Do_lanc_c

         if ( eigenvalue >= 0.0d0 ) then
            ret = 60000 + pas 
         else                          ! Else of eigenvalue.
            ! New fpar and fperp. 
            call force_projection( fpar, perp_force, fperp,&
               &   ftot, force, projection )
            ret = 30000 + pas
         end if
      end if

   end select

   select case ( ret )

   case( 10000 : 39999 )

      converg = 'CONVERGED' 
      success = .true.
      ! We write the configuration in a sad.... file
      call convert_to_chain( mincounter, 4, scounter )
      fname = SADDLE // scounter
      if (iproc == 0 ) call store( fname ) 
      conf_saddle = fname

   case default  
      converg = 'FAILED'
      success = .false.
   end select

   delta_e = saddle_energy - ref_energy
   ! Final report of saddle point 
   if ( iproc == 0 ) then
      open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
         & action = 'write', position = 'append', iostat = ierror )
      write(FLOG,"(/' ','SADDLE',i5, a10,' |ret ',i6,' |delta energy= '," //  &
         &       "f9.4, ' |force_(tot,par,perp)= ', 3f10.4," // &
         &       "' |eigenval=',f9.4,' |npart= ', i4,' |delr= ', f8.3,' |evalf='," // &
         &       "i6, ' |')")                                          & 
         & mincounter, adjustr(converg), ret, delta_e, ftot,  &  
         & fpar, fperp, eigenvalue, npart, delr, evalf_number
      write(*,"(/' ','BART: SADDLE',i5, a10,' |ret ',i6,' |delta energy= '," //  &
         &    "f9.4, ' |force_(tot,par,perp)= ', 3f10.4," //       &
         &    "' |eigenval=',f9.4,' |npart= ', i4,' |delr= ', f8.3,' |evalf='," // &
         &    "i6, ' |')")                                          & 
         & mincounter, adjustr(converg), ret, delta_e, ftot,  &  
         & fpar, fperp, eigenvalue, npart, delr, evalf_number

      if ( success ) then  
         ! Write  
         write(*,*) 'BART: Configuration stored in file ',fname
         write(FLOG,'(1X,A34,A17)') ' - Configuration stored in file : ', trim(fname)
         write(FLOG,'(1X,A34,(1p,e17.10,0p))') &
            &   ' - Total energy Saddle (eV)     : ', saddle_energy
         write(FLOG,*) ' '
      end if
      close(FLOG) 
   end if

END SUBROUTINE end_report 


!> ART clean_wavefunction 
!! We only clean the wave function. No calculation of the eigenvector.
!! This is not going to be saved in the restart file.
subroutine clean_wavefunction( a1, call_fp )

   use defs
   use saddles
   use bigdft_forces
   use lanczos_defs
   implicit None

   !Arguments
   real(kind=8), intent(in) :: a1
   logical,      intent(in) :: call_fp

   !Local variables
   real(kind=8) :: ftot2
   real(kind=8), dimension(3) :: boxl
   real(kind=8), dimension(VECSIZE) :: perp_force   ! Perpendicular force...
   !_______________________ 


   ! we only clean the wave function. No calculation of the eigenvector. This is not
   ! going to be saved in the restart file.

   boxl = box * scala     ! We compute at constant volume.

   new_wf = .True.  
   call calcforce( NATOMS, pos, boxl, force, total_energy, evalf_number, .false. )
   new_wf = .False.

   delta_e = total_energy - ref_energy

   if ( .not. call_fp ) then 
      m_perp = 0
      try = 0
      fpar = 0.0d0
      fperp= 0.0d0
      eigenvalue = 0.0d0
      ftot2 = dot_product (force,force)
      ftot  = sqrt( ftot2 ) 
   else
      call force_projection( fpar, perp_force, fperp, ftot, force, projection ) 
   end if
   call write_step ( 'C', 0, a1, total_energy )

END SUBROUTINE  clean_wavefunction
