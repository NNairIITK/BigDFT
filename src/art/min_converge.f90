!> @file
!! @author
!!    Copyright (C) 2001 Normand Mousseau
!!    Copyright (C) 2010-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!! Modified by Laurent in 2011 for working with QM/MM and FIRE


!> ART min_converge
!!   Minimizes the energy at constant volume.
!!   This minimization is done with only a minimal knowledge of the 
!!   physics of the problem so that it is portable
subroutine min_converge ( success )

   use defs
   use bigdft_forces
   implicit none

   !Arguments
   logical, intent(out) :: success

   !Local variables
   integer :: i, ierror
   real(kind=8),dimension(3) :: boxl

   real(kind=8), dimension(3*natoms)         :: pos_temp
   integer :: nat
   integer, dimension(natoms) :: numnei
   integer, dimension(natoms,maxnei) :: nei
   real(kind=8), dimension(3) :: invbox
   integer :: j,k 
   logical, dimension(natoms) :: is_at_quantum
   real(kind=8) :: xij,yij,zij,rij2
   !_______________________
   if ( iproc == 0 ) then              ! Report
      open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
         & action = 'write', position = 'append', iostat = ierror )
      write(FLOG,'(1X,A)') ' RELAXATION'
      close(FLOG) 
   end if

   if (energy_type == "BIG" ) then
      ! Relaxation subroutine in bigdft_forces.f90                     
      boxl = box * scala          ! We compute at constant volume.
      ! Report.
      call mingeo( pos, force, boxl, evalf_number, total_energy, success )

      if ( clean_wf ) then 
         call check_force_clean_wf ( pos, boxl, evalf_number, total_energy, success )
      end if 

   else if (energy_type == "BSW") then
      pos_temp = 0.0d0
      invbox = 1.0d0/box
      is_at_quantum = .false. !vectorial operation

      !now passivate the quantum box, this will take a few steps. should make
      !a separate sub-routine

      !we should also make sure that the set we send to bigDFT forces did not
      !change in size since the last time


      call neighbours(natoms,pos,box,boundary,maxnei,numnei, nei)
      nat = 0
      do i = 1,nbr_quantum
         is_at_quantum(i) = .true.
         nat = nat + 1
         pos_temp(i) = pos(i)
         pos_temp(i+natoms) = pos(i+natoms)
         pos_temp(i+natoms+natoms) = pos(i+natoms+natoms)
      enddo

      do i = 1,nbr_quantum
         if (passivate) then 
            do j = 1,numnei(i)
               k = nei(i,j)
               if ( .not. is_at_quantum(k)) then 
                  xij = pos(k)-pos(i) - box(1) * nint((pos(k)-pos(i))*invbox(1))
                  yij = pos(k+natoms)-pos(i+natoms)
                  zij = pos(k+2*natoms)-pos(i+2*natoms) - box(3) * nint((pos(k+2*natoms)-pos(i+2*natoms))*invbox(3))
                  rij2 = xij*xij + yij*yij + zij*zij
                  if (rij2 .lt. 2.7d0*2.7d0) then
                     nat = nat + 1
                     pos_temp(nat) = pos(i) + 0.5d0*xij
                     pos_temp(nat+natoms) = pos(i+natoms) + 0.5d0*yij
                     pos_temp(nat+natoms+natoms) = pos(i+2*natoms) + 0.5d0*zij !we passivate with hydrogene at this distance
                  endif
               endif
            enddo
         endif
      enddo

      call  min_converge_fire(success)
      elseif (energy_type == "SWP") then
      call min_converge_fire(success)
      elseif (energy_type == "OTF") then
      call fit_SW_potential()
      call  min_converge_fire(success)
      elseif (energy_type == "BAY") then
      call  min_converge_fire(success)
   endif

   if ( iproc == 0 ) then 
      if ( .not. success ) then 
         open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
            & action = 'write', position = 'append', iostat = ierror )
         write(FLOG,'(1X,A)') "Minimization exited before the geometry optimization converged,"
         write(FLOG,'(1X,A)') "this minimum will be rejected."
         close(FLOG) 
      end if

      write(*,"('',' BART: Relaxed energy : ',(1p,e17.10,0p))") total_energy
   end if 


END SUBROUTINE min_converge


!> ART check_min
subroutine check_min( stage )
   use defs
   use lanczos_defs
   implicit none

   !Arguments
   character(len=1), intent(in) :: stage

   !Local variables
   integer :: i, ierror, repetition
   logical :: new_projection
   real(kind=8) :: a1
   real(kind=8) :: min_energy ! First reference energy in lanczos


   ! We check how it changes the energy of the system by applying the projection.
   IN_MINIMUN = .True.
   ! Report 

   if ( iproc == 0 ) then 
      open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
         & action = 'write', position = 'append', iostat = ierror )
      write(FLOG,*) ' Starting Lanczos'
      close(FLOG) 
   end if

   new_projection = .true.            ! We do not use any previously computed direction. 

   if ( .not. setup_initial ) then
      ! We call lanczos twice.
      repetition = 2
   else 
      ! if not, four times.
      repetition = 3 
   end if

   if ( iproc==0 ) write(*,*) "BART: INIT LANCZOS"  !debug
   do i = 1, repetition
      call lanczos( NVECTOR_LANCZOS_H, new_projection , a1 )
      ! Report
      if ( iproc == 0 ) then 
         open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
            & action = 'write', position = 'append', iostat = ierror )
         if ( i == 1 ) then           ! Our reference energy for the report.
            min_energy = lanc_energy
            write(FLOG,'(1X,A8,(1p,e17.10,0p),A12,1pe8.1,A3)') ' Em= ', min_energy, ' ( gnrm = ', my_gnrm, ' )'
            write(FLOG,'(A39)') '   Iter     Ep-Em (eV)   Eigenvalue  a1' 
         end if 
         write(FLOG,'(I6,3X,(1p,e10.2,0p),4X,F12.6,1X,F7.4)') i, proj_energy-min_energy, eigenvalue, a1
         close(FLOG) 
         write(*,*) 'BART: Iter ', i, ' : ', lanc_energy, proj_energy,  eigenvalue, a1  
      end if
      ! Now we start from the previous direction. 
      new_projection= .false.   
      ! let's see the projection
      if ( setup_initial ) call print_proj ( i, stage, projection, eigenvalue, DEL_LANCZOS )
   end do

   if (energy_type == "SWP") call reset_SW_potential()
   ! Report 
   if ( iproc == 0 ) then 
      write(*,*) "BART: END  LANCZOS"  !debug
      open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
         & action = 'write', position = 'append', iostat = ierror )
      write(FLOG,*) ' Done Lanczos'
      close(FLOG) 
   end if

   ! Default value in the activation part is false.
   IN_MINIMUN = .False.

END SUBROUTINE check_min


!>  This module defines a number of parameters used during the minimization
MODULE minimization_sd
   implicit none
   save
   integer, parameter :: MAX_ITER = 1000
   real(kind=8), parameter :: FTHRESHOLD = 3.5D-1
   real(kind=8), parameter :: STEPSIZE = 1.0D-4  ! Size in angstroems

   real(kind=8), parameter :: FTHRESH2 = FTHRESHOLD * FTHRESHOLD
END MODULE minimization_sd


subroutine min_converge_sd(minimized)
   use minimization_sd
   use defs
   implicit none

   logical, intent(inout)  :: minimized 

   real(kind=8),  dimension(VECSIZE):: forceb, posb, tmp_pos
   integer :: iter, i, npart
   real(kind=8) :: current_ftot2, ftot,ftot2, step, delr, current_energy
   logical :: conv

   conv = .false. !not sure about this variable

   ! We compute at constant volume
   call calcforce(NATOMS, pos, box, force, total_energy, evalf_number, conv)
   ftot2 = 0.0d0
   do i=1, VECSIZE
      ftot2 = ftot2 + force(i) * force(i)
   end do
   current_ftot2 = ftot2
   current_energy = total_energy ! if quantum, this energy is only quantum
   if (iproc == 0) write(*,*)  'initial energy : ',  total_energy


   step = STEPSIZE
   do iter = 1, MAX_ITER

      ! Apply PBC 	
      posb = pos + step * force
      tmp_pos = pos
      pos = posb
      call calcforce(NATOMS, pos, box, forceb, total_energy, evalf_number, conv)
      pos = tmp_pos
      ftot2 = 0.0d0
      do i=1, VECSIZE
         ftot2 = ftot2 + forceb(i) * forceb(i)
      end do

      ! if(ftot2 < current_ftot2 .or. total_energy < current_energy ) then
      if(ftot2 < current_ftot2  ) then    

         pos = posb  
         force = forceb
         step = 1.2 * step
         current_ftot2 = ftot2
         current_energy = total_energy

      else 
         step = 0.6 * step
      endif
      if(ftot2 < FTHRESH2) exit
      call displacement( posref, pos, delr, npart )
      ! Write 

      if (modulo(iter,5) == 0 ) then
         call write_step ( 'M', iter, 0.0d0, current_energy )
         if ( SAVE_CONF_INT ) call save_intermediate( 'M' )

      endif

   end do

   ftot = sqrt(ftot2)  
   if (ftot < FTHRESHOLD ) then
      if (iproc == 0) write(*,*) 'Minimization successful   ftot : ', ftot
      if (iproc == 0) write(*,*)  'final energy :', total_energy
      minimized = .true. 
   else
      if (iproc == 0) write(*,*) 'Minimization failed   ftot : ', ftot
      minimized = .false. 
   endif
END SUBROUTINE min_converge_sd


!> Implementation of the damped MD based geometry optimizer FIRE, PRL 97, 170201 (2006)
!! The MD-Integrator is the common velocity verlet, all masses are equal to 1.d0
!! Implemented in August 2010, Maximilian Amsler, Basel University 
!! Suggestion for maximal timestep as tmax=2*pi*sqrt(alphaVSSD)*1.2d-1
!! Choose the initial timestep as tinit=tmax*0.5d0
subroutine min_converge_fire(success)
   use defs
   use saddles
   use lanczos_defs

   implicit none

   !Arguments
   logical :: success
   !Local variables
   real(kind=8) :: fnrm
   real(kind=8) :: fmax,vmax
   real(kind=8), parameter :: pi = 3.14159265d0
   integer :: check
   integer :: iat

   logical :: conv

   !Fire parameters:
   real(kind=8) :: alpha,P,finc,fdec,falpha,alphastart,dt,dtmax,vnrm
   real(kind=8) :: velcur(3*natoms), velpred(3*natoms),poscur(3*natoms),pospred(3*natoms)
   real(kind=8) :: fcur(3*natoms),fpred(3*natoms),mass(3*natoms)
   real(kind=8) :: ecur,epred,anoise
   integer:: Nmin,nstep,it
   integer,parameter :: max_iter=100
   real(kind=8),parameter :: norm_criterium = 0.008d0
   real(kind=8),parameter :: fmax_criterium = 0.020d0
   integer :: miter

   miter = 0
   conv = .false.


   check=0
   !Set FIRE parameters
   Nmin=5
   finc=1.1d0
   fdec=0.5d0
   !  alphastart=0.25d0
   alphastart =0.1d0

   anoise=1.0d-8


   alpha=alphastart
   falpha=0.99d0
   nstep=1

   dtmax=2*pi*dsqrt(alpha)*1.2d-1
   dt = dtmax*0.5d0

   success=.false.
   fnrm=1.d10
   velcur=0.0d0
   poscur=pos

   call calcforce(natoms,pos,box,fpred,total_energy,evalf_number,conv)

   fcur=force
   mass=1.0d0
   ecur=total_energy
   epred=total_energy


   do it=1,max_iter
      miter = miter + 1
      pas = pas + 1
      do iat=1,3*natoms
         pospred(iat)=poscur(iat)+dt*velcur(iat)+dt*dt*0.5d0*fcur(iat)/mass(iat)
      enddo

      call calcforce(natoms,pospred,box,fpred,total_energy,evalf_number,conv)
      force=fpred
      call fnrmmax_fire(fpred,fnrm,fmax,natoms)
      !  call convcheck(fnrm,fmax,fluct*in%frac_fluct,in%forcemax,check)

      do iat=1,3*natoms
         velpred(iat)=velcur(iat)+0.5d0*dt*(fpred(iat))/mass(iat)+0.5d0*dt*fcur(iat)/mass(iat)
      enddo
      P=dot_product(fpred,velpred)
      call fnrmmax_fire(velpred,vnrm,vmax,natoms)

      force = fpred
      ftot = dsqrt(fnrm)
      delta_e = total_energy - ref_energy
      fpar = fmax
      fperp = 0.0d0
      eigenvalue = 0.0d0
      ! Magnitude of the displacement (utils.f90).
      call displacement( posref, poscur, delr, npart )
      ! Write 

      if (modulo(miter,5) == 0 ) then
         call write_step ( 'M', miter, 0.0d0, total_energy )
         if (iproc == 0) write(*,*) "fnrm",fnrm,"fmax",fmax
         pos = pospred
         if ( SAVE_CONF_INT ) call save_intermediate( 'M' )
      endif

      if( fnrm < norm_criterium .or. fmax < fmax_criterium) then
         pos = pospred
         success = .true.
         exit
      endif

      !Update variables
      fcur=fpred
      poscur=pospred
      !Normal verlet velocity update
      !  velcur=velpred

      !!FIRE Update
      call fnrmmax_fire(fpred,fnrm,fmax,natoms)
      fnrm=sqrt(fnrm)
      call fnrmmax_fire(velpred,vnrm,vmax,natoms)
      vnrm=sqrt(vnrm)
      !Modified velocity update, suggested by Alireza
      !velcur(:)=(1.0d0-alpha)*velpred(:)+fpred(:)*min(alpha*vnrm/fnrm,2.0d0*in%betax)!alpha*fpred(:)/fnrm*vnrm
      !Original FIRE velocitiy update
      velcur(:)=(1.0d0-alpha)*velpred(:)+alpha*fpred(:)/fnrm*vnrm
      if(P.gt.-anoise*vnrm .and. nstep.gt.Nmin) then
         dt=min(dt*finc,dtmax)
         !         alpha=max(alpha*falpha,0.1d0) !Limit the decrease of alpha
         alpha=alpha*falpha
         elseif(P.le.-anoise*vnrm) then
         nstep=0
         dt=dt*fdec
         velcur=0.d0
         alpha=alphastart
      endif
      nstep=nstep+1
   enddo

   return
END SUBROUTINE min_converge_fire


subroutine fnrmmax_fire(force,fnrm,fmax,natoms)
   implicit none
   integer,intent(in) :: natoms
   real(kind=8),dimension(3*natoms),intent(in) :: force
   real(kind=8),intent(out) :: fnrm,fmax

   integer :: i

   real(kind=8) :: nombre

   fnrm = 0.0d0
   fmax = 0.0d0
   !  force2 = dot_product(force)
   do i = 1,natoms
      nombre = force(i)**2 + force(i+natoms)**2 + force(i+natoms+natoms)**2
      if (nombre > fmax) fmax = nombre
      fnrm = fnrm+nombre
   enddo
   fmax = dsqrt(fmax)

END SUBROUTINE fnrmmax_fire
