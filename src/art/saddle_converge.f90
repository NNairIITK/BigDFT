!> art/saddle_converge
!! :
!!    In the art method, converge to the saddle point
!!    This subroutine bring the configuration to a saddle point. It does that
!!    by first pushing the configuration outside of the harmonic well, using
!!    the initial direction selected in find_saddle. Once outside the harmonic
!!    well, as defined by the appearance of a negative eigenvalue (or
!!    reasonnable size) the configuration follows the direction corresponding
!!    to this eigenvalue until the force components parallel and perpdendicular
!!    to the eigendirection become close to zero.
!!
!!
!! Copyright:
!!
!!    Copyright (C) Normand Mousseau, June 2001
!!    Copyright (C) 2010 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!!
!!
subroutine saddle_converge(ret, saddle_energy, fpar, fperp)
  use random
  use defs
  use saddles
  use lanczos_defs
  use bigdft_forces
  implicit none

  !Arguments
  integer, intent(out) :: ret
  real(8), intent(out) :: saddle_energy, fpar, fperp
  logical :: new_projection
  !Local variables
  integer :: i, k, kter, iter, iperp, k_rejected, eigen_rejected, nat
  integer :: maxvec,npart,itry, iter_init, kter_init, ierror
  real(8) :: fdotinit, fperp2,current_fperp
  real(8) :: step,delr
  real(8), parameter :: one = 1.0d0
  real(8) :: boxl, current_energy, ftot
  real(8), dimension(VECSIZE) :: posb, perp_force, forceb,perp_forceb

  ! We compute at constant volume
  boxl = box * scala

  if (restart .and. (state_restart == 1) ) then ! Restarts in harmonic well
     initial_direction = direction_restart
     kter_init = iter_restart

     if (iproc.eq.0) then 
        write(*,*) 'Restart'
        write(*,*) 'in harmonic well '
        write(*,*) 'kter : ', kter_init
        write(*,*) 'pos: ', pos(1), pos(2), pos(3)
     endif
     restart = .false.
  else
     kter_init = 0
  endif

  if ( (.not. restart) .and. NEW_EVENT ) then 

     eigenvalue = 0.0d0
     call calcforce(NATOMS,pos,boxl,force,current_energy)
     
     ! We now project out the direction of the initial displacement from the
     ! minimum from the force vector so that we can minimize the energy in
     ! the direction perpendicular to the direction of escape from harmonic well
     fdotinit= 0.0d0
     do i=1, VECSIZE
        fdotinit = fdotinit + force(i) * initial_direction(i)
     end do
     perp_force  = force - fdotinit * initial_direction  ! Vectorial force
     
     step = 0.4d0*INCREMENT
     
     new_projection = .true.     !! COSMIN
     do kter = kter_init, MAXKTER
        
        k = 0
        k_rejected = 0
        
        eigen_rejected = 0
        ! We relax perpendicularly using a simple variable-step steepest descent 
        do 
           posb = pos + step * perp_force
           
           call calcforce(NATOMS,posb,boxl,forceb,total_energy)
           evalf_number = evalf_number + 1
           
           fdotinit= 0.0d0
           do i=1, VECSIZE
              fdotinit = fdotinit + forceb(i) * initial_direction(i)
           end do
           perp_forceb  = forceb - fdotinit * initial_direction  ! Vectorial force
           
           fperp2 = 0.0d0
           do i=1, VECSIZE
              fperp2 = fperp2 + perp_forceb(i) * perp_forceb(i)
           end do
           fperp = sqrt(fperp2)
           
           if(total_energy < current_energy ) then
              pos = posb
              force = forceb
              perp_force = perp_forceb
              step = 1.2 * step
              current_energy = total_energy
              k = k + 1
              k_rejected = 0
           else
              step = 0.6 * step
              k_rejected = k_rejected + 1
           endif
           
           if(fperp2 < FTHRESH2 .or. k > MAXKPERP .or. k_rejected > 5) exit
        end do
        
        ! We now move the configuration along the initial direction and check the lowest
        ! eigendirection 
        pos = pos + BASIN_FACTOR * INCREMENT * initial_direction    !Vectorial operation

        ! Broadcast pos to all nodes
        nat = 3 * NATOMS
        call MPI_Bcast(pos,nat,MPI_REAL8,0,MPI_COMM_WORLD,ierror)

        
        ! We start checking of negative eigenvalues only after a few steps
        if( kter>= KTER_MIN ) then 
           new_projection = .true.    ! We do not  use previously computed lowest direction as see
           maxvec = NVECTOR_LANCZOS
           call lanczos(maxvec,new_projection)
           new_projection = .false.     !! COSMIN
        endif
        current_energy = total_energy ! As computed in lanczos routine
        call displacement(posref, pos, delr,npart)
        if (iproc .eq. 0 ) then
           write(*,"(' ','kter: ',i4,'  Energy: ',f16.6,'  e-val: ',f12.6,'  delr: ',f12.6)") &
                &  kter,current_energy, eigenvalue, delr
           
           open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
           write(FLOG,"(' ','kter: ',i4,'  Energy: ',f16.6,'  e-val: ',f12.6,'  delr: ',f12.6)") &
                &  kter,current_energy, eigenvalue, delr
           close(FLOG)
        endif

        state_restart = 1
        if (iproc .eq. 0 ) call save_state(state_restart, kter+1, initial_direction)

        if(eigenvalue <  EIGEN_THRESH) exit
     end do
     
     
     ! The configuration is now out of the harmonic well, we can now bring
     ! it to the saddle point. Again, we split the move into a parallel and 
     ! perpendicular contribution.
     
     ! First, we must now orient the direction of the eigenvector corresponding to the
     ! negative eigendirection (called projection) such that it points away from minimum. 
     
     fpar= 0.0d0
     do i=1, VECSIZE
        fpar = fpar + forceb(i) * projection(i)
     end do
     if(fpar > 0.0d0 ) projection = -1.0d0 * projection
     pos = pos +  1.0d0 * INCREMENT * projection 
     
     fdotinit= 0.0d0
     do i=1, VECSIZE
        fdotinit = fdotinit + forceb(i) * projection(i)
     end do
     perp_force  = force - fdotinit * projection  ! Vectorial force
     
     fperp2 = 0.0d0
     do i=1, VECSIZE
        fperp2 = fperp2 + perp_force(i) * perp_force(i)
     end do
     current_fperp = sqrt(fperp2)

     iter_init = 1
  else if (restart .and. (state_restart == 2) ) then
        iter_init = iter_restart
        projection = direction_restart
        call calcforce(NATOMS,pos,boxl,force,current_energy)

        fpar= 0.0d0
        do i=1, VECSIZE
           fpar = fpar + force(i) * projection(i)
        end do

        perp_force  = force - fpar * projection  ! Vectorial force
        fperp2 = 0.0d0
        do i=1, VECSIZE
           fperp2 = fperp2 + perp_force(i) * perp_force(i)
        end do
        fperp = sqrt(fperp2)
        current_fperp = fperp
       
        if (iproc .eq. 0 ) write(*,*) 'Restart = 2'
        restart = .false.
  else if (.not.new_event) then       ! This is a convergence event
                                      ! We must compute the eigenvector with precision

     posref = pos

     maxvec = NVECTOR_LANCZOS
     new_projection = .true. 
     do i =1, 5
        call lanczos(maxvec,new_projection)
        new_projection = .false.
     end do

     call calcforce(NATOMS,pos,boxl,force,current_energy)
     fpar = dot_product(force,projection)
     if(fpar > 0.0d0 ) projection = -1.0d0 * projection

     perp_force  = force - fpar * projection  ! Vectorial force
     fperp = sqrt( dot_product(perp_force,perp_force) )
     current_fperp = fperp

     iter_init = 100
     
   else
     write(*,*) 'Problem with restart and state_restart : '
     write(*,*) 'restart = ', restart
     write(*,*) 'state_restart = ', state_restart
     stop

  endif

  do iter = iter_init, MAXITER 
    step = 0.04d0*INCREMENT
   
    reject = .false.
    itry = 0
    iperp = 0
    ! The next loop is on the perpendicular direction, again using a 
    ! variable-step steepest descent
    do  
       posb = pos + step * perp_force
       call displacement(posb, pos, delr,npart)

       call calcforce(NATOMS,posb,boxl,forceb,total_energy)
       evalf_number = evalf_number + 1
       
       fpar= 0.0d0
       do i=1, VECSIZE
          fpar = fpar + forceb(i) * projection(i)
       end do
       perp_force  = forceb - fpar * projection  ! Vectorial force
       
       fperp2 = 0.0d0
       do i=1, VECSIZE
          fperp2 = fperp2 + perp_force(i) * perp_force(i)
       end do
       fperp = sqrt(fperp2)
       
       if((fperp < current_fperp) ) then
          pos = posb
          force = forceb
          step = 1.2 * step
          current_energy = total_energy
          current_fperp = fperp
          iperp = iperp + 1
       else
          step = 0.6 * step
       endif
       itry = itry + 1
       
       if(fperp2 < FTHRESH2 .or. iperp > (iter - 10)  .or. iperp > MAXIPERP .or.  itry > 5) exit
       !       if(fperp2 < FTHRESH2 .or. iperp > 5 .or.  itry > 5) exit
       
    end do    
    new_projection = .false.    ! We start from the previuos direction each time
    maxvec = NVECTOR_LANCZOS
    call lanczos(maxvec,new_projection)
    
    if (iproc.eq.0) then 
       write(*,*) 'eigenvalue : ', eigenvalue
       write(*,*) 'eigenvals: ', (eigenvals(i),i=1,4)
    endif
    
    if (reject) eigen_rejected = eigen_rejected + 1  
    
    fpar = dot_product(force,projectioN)
    
    ! We now move the configuration along the eigendirection corresponding
    ! to the lowest eigenvalue
    
    if(abs(fpar) > EXITTHRESH*0.5 ) then 
       pos = pos - sign(one,fpar) * 2.0d0 * INCREMENT * projection / sqrt(1.0d0*iter) 
    else
       pos = pos - sign(one,fpar) * 1.0d0 * INCREMENT * projection / sqrt(1.0d0*iter) 
    endif

    ! Broadcast pos to all nodes
    nat = 3 * NATOMS
    call MPI_Bcast(pos,nat,MPI_REAL8,0,MPI_COMM_WORLD,ierror)

    
    current_energy = total_energy ! As computed in lanczos routine
    call displacement(posref, pos, delr,npart)

    if(iproc.eq.0) then
       open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
       write(FLOG,"(' ','iter: ',i4,'  iper: ',i4,'  ener: ',f12.4,'  fpar: ',f12.6,'  fperp: ',f12.6,&
            & '  e-val: ', f12.6,'  delr: ',f10.4,'  npart: ',i4,'  evalf: ',i6)")  &
            & iter, iperp, total_energy, fpar, fperp, eigenvalue, delr, npart, evalf_number
       close(FLOG) 
       write(*,"(' ','iter: ',i4,'  iper: ',i4,'  ener: ',f12.4,'  fpar: ',f12.6,'  fperp: ',f12.6,&
            & '  e-val: ', f12.6,'  delr: ',f10.4,'  npart: ',i4,'  evalf: ',i6)")  &
            & iter, iperp, total_energy, fpar, fperp, eigenvalue, delr, npart, evalf_number
    endif

    saddle_energy = current_energy
    if ( (abs(fpar)+fperp)< EXITTHRESH)  then

       if (USE_DIIS) then 
          if (iproc .eq. 0 ) then
             open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
             write(FLOG,*) "Call DIIS"
             close(FLOG)
          endif

          call apply_diis(current_energy,ftot)
          if (iproc.eq.0) print *, 'current energy ', current_energy
          saddle_energy = current_energy

          if (DIIS_CHECK_EIGENVEC) then
             new_projection = .true.
             maxvec = NVECTOR_LANCZOS
             do i=1, 4
                call lanczos(maxvec,new_projection)
                if (iproc.eq.0) then 
                   write(6,*) 'i: ', i, 'eigenvalue: ', eigenvalue," energy:",total_energy
                   open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',&
                        & iostat=ierror)
                   write(FLOG,*) 'i: ', i, 'eigenvalue: ', eigenvalue," energy:",total_energy
                   close(FLOG) 
                endif
                new_projection = .false.
                if (eigenvalue .lt. 0.0d0 ) exit
             end do
             if (eigenvalue .ge. 0.0 .or. ftot .gt. DIIS_FORCE_THRESHOLD) then
                ret = 60000 + iter
             else
                ret = 20000 + iter
             endif
          else
             ret = 20000 + iter
             if (ftot .gt. DIIS_FORCE_THRESHOLD) ret = 70000+iter
          endif
       else
          ret = 20000 + iter
       endif
       exit
    else if ( (abs(fpar)<0.1*EXITTHRESH) .and. (fperp<EXITTHRESH) ) then

       if (USE_DIIS) then
          if (iproc.eq. 0) then 
             open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
             write(FLOG,*) "Call DIIS (2)"
             close(FLOG)
          endif
          call apply_diis(current_energy,ftot)
          saddle_energy = current_energy

          if (DIIS_CHECK_EIGENVEC) then
             new_projection = .true.
             maxvec = NVECTOR_LANCZOS
             do i=1, 4
                call lanczos(maxvec,new_projection)
                if ( iproc.eq. 0) then
                   write(6,*) 'i: ', i, 'eigenvalue: ', eigenvalue," energy:",total_energy
                   open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',&
                        & iostat=ierror)
                   write(FLOG,*) 'i: ', i, 'eigenvalue: ', eigenvalue," energy:",total_energy
                   close(FLOG) 
                endif
                new_projection = .false.
                if (eigenvalue .lt. 0.0d0 ) exit
             end do
             if (eigenvalue .ge. 0.0 .or. ftot .gt. DIIS_FORCE_THRESHOLD) then
                ret = 60000 + iter
             else
                ret = 10000 + iter
             endif
          else
             ret = 10000 + iter
             if (ftot .gt. DIIS_FORCE_THRESHOLD) ret = 70000+iter
          endif
       else
          ret = 0000 + iter
       endif
       exit
    else if (eigenvalue>0.0) then
      ret = 60000 + iter
      exit
    endif

    call calcforce(NATOMS,pos,boxl,force,total_energy)  !! COSMIN
    perp_force  = forceb - fpar * projection  ! Vectorial force
    fperp2 = 0.0d0
    do i=1, VECSIZE
      fperp2 = fperp2 + perp_force(i) * perp_force(i)
    end do
    fperp = sqrt(fperp2)
    current_fperp = fperp

    state_restart = 2
    if (iproc .eq. 0 ) call save_state(state_restart, iter+1, projection)
    
  end do
END SUBROUTINE saddle_converge

