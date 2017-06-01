!> @file
!!   Initialize the potential
!!
!! @author
!!    Copyright (C) 2001 Normand Mousseau
!!    Copyright (C) 2010-2016 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
!!  Modified by Laurent Karim Beland, UdeM, 2011. For working with QM/MM !!
subroutine initialize_potential( )
   use defs
   use bigdft_forces
   use wrapper_MPI
   implicit None

   !Local variables

   if (energy_type == "BIG") then
      call bigdft_init_art(nbr_quantum, iproc, nproc, my_gnrm,passivate,natoms )
      ! First force calculation.
      call calcforce( NATOMS, pos, boxref, force, total_energy, evalf_number, .false. )
      elseif (energy_type == "SWP") then
      call init_potential_SW()
      call calcforce( NATOMS, pos, boxref, force, total_energy, evalf_number, .false. )
      elseif (energy_type == "BSW") then
      call init_potential_SW()
      call bigdft_init_art(nbr_quantum, iproc, nproc, my_gnrm,passivate,natoms )
      ! First force calculation.
      call calcforce( NATOMS, pos, boxref, force, total_energy, evalf_number, .false. )
      elseif (energy_type == "OTF") then
      call init_potential_SW()
      call bigdft_init_art(nbr_quantum, iproc, nproc, my_gnrm,passivate,natoms )
      ! First force calculation.
      energy_type = "BSW"
      call calcforce( NATOMS, pos, boxref, force, total_energy, evalf_number, .false. )
      energy_type = "OTF"
      !now we define the fitting zone
      call initialize_fitting_zone()
      call fit_SW_potential()
      call mpibarrier()
      elseif (energy_type == "BAY") then
      call init_potential_SW()
      call bigdft_init_art(nbr_quantum, iproc, nproc, my_gnrm,passivate,natoms )
      ! First force calculation.
      energy_type = "BSW"
      call calcforce( NATOMS, pos, boxref, force, total_energy, evalf_number, .false. )
      energy_type = "BAY"
      if (iproc==0) write(*,*) "real forces below"
      if (iproc==0) write(*,*) force(1),force(1+natoms),force(1+2*natoms)
      if (iproc==0) write(*,*) ""

      !we prepare the first bayes data set
      call create_bayes_data()
      call mpibarrier()
   else
      write(*,*) "You have not chosen a proper energy type. Choose SWP or BIG in ENERGY_CALC"
      stop

   endif

END SUBROUTINE initialize_potential


!> Initialize the zone where we fit classical SW to the rest
subroutine initialize_fitting_zone()
   use defs

   implicit none
   integer, dimension(natoms) :: numnei
   integer, dimension(natoms,maxnei) :: nei
   integer :: i,j,k


   if (.not. allocated(should_fit)) allocate(should_fit(natoms))
   should_fit = .false. !vectorial operation

   call neighbours(natoms,pos,box,boundary,maxnei,numnei, nei)

   !we will fit all quantum atoms and their neighbours
   do i = 1,natoms
      should_fit(i) = .true.
      do j = 1,numnei(i)
         k = nei(i,j)
         should_fit(k) = .true.
      enddo
   enddo
   nbr_to_fit = 0
   do i = 1,natoms
      if (should_fit(i)) nbr_to_fit = nbr_to_fit + 1
   enddo

END SUBROUTINE initialize_fitting_zone


!> Finalize the potential
subroutine finalise_potential( )

   use bigdft_forces
   use defs, only : energy_type

   implicit none

   if (energy_type == "BIG" .or. energy_type == "BSW") call bigdft_finalise( )

END SUBROUTINE finalise_potential
