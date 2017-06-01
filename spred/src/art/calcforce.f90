!> @file
!!   Subroutine which calls the right for type inside art
!! @author
!!   Written by Laurent Karim Beland, UdeM 2011!!
!!   Copyright (C) 2010-2011 BigDFT group, Normand Mousseau
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 


!> Subroutine which calls the right for type inside art
subroutine calcforce(nat, posa, boxl, forca, energy, evalf_number, conv )
   use bigdft_forces
   use defs, only : energy_type
   implicit none

   !Arguments
   integer,      intent(in)                            :: nat
   real(kind=8), intent(in),  dimension(3*nat)         :: posa
   real(kind=8), dimension(3), intent(inout)           :: boxl
   real(kind=8), intent(out), dimension(3*nat)         :: forca
   real(kind=8), intent(out)                           :: energy
   integer,      intent(inout)                         :: evalf_number
   logical,      intent(in)                            :: conv


   if(energy_type == "SWP")  then
      call SWcalcforce(nat,posa,boxl,forca, energy)
      evalf_number = evalf_number +1
   endif

   if(energy_type == "BIG")                   call calcforce_bigdft( posa, forca, boxl, energy, evalf_number, conv )
   if(energy_type == "BSW")                   call calcforce_mix(posa,forca,energy,evalf_number, conv)
   if(energy_type == "OTF")                   call SWcalcforce(nat,posa,boxl,forca, energy)
   if(energy_type == "BAY")                   call calcforce_bayes(nat,posa,boxl,forca, energy)


END SUBROUTINE calcforce


! bart/calcfore
!! FUNCTION
!!   Subroutine which assigns QM and MM forces to each atom
!! COPYRIGHT
!!   Copyright (C) 2010 BigDFT group, Normand Mousseau
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 
!! SOURCE
!!

subroutine calcforce_mix(posa,forca,energy,evalf_numbera, conva) !be careful, the energy is bogus
   use defs
   use bigdft_forces
   implicit none
   real(kind=8), intent(in), dimension(3*natoms) :: posa
   real(kind=8), intent(out), dimension(3*natoms)         :: forca
   real(kind=8), intent(out)                           :: energy
   integer,      intent(inout)                         :: evalf_numbera
   logical,      intent(in)                            :: conva

   real(kind=8), dimension(3*natoms)  :: force_class
   real(kind=8), dimension(:),allocatable :: force_quant

   real(kind=8), dimension(:), allocatable ::  posquant
   real(kind=8), dimension(3*natoms)         :: pos_temp
   integer :: nat
   integer, dimension(natoms) :: numnei
   integer, dimension(natoms,maxnei) :: nei
   real(kind=8), dimension(3) :: invbox
   integer :: i,j,k 
   logical, dimension(natoms) :: is_at_quantum
   real(kind=8) :: xij,yij,zij,rij2

   invbox = 1.0d0/box
   is_at_quantum = .false. !vectorial operation


   !first do a classical potential force calc
   call SWcalcforce(natoms,posa,box,force_class, energy)

   !now passivate the quantum box, this will take a few steps. should make
   !a separate sub-routine

   !we should also make sure that the set we send to BigDFT forces did not
   !change in size since the last time


   call neighbours(natoms,pos,box,boundary,maxnei,numnei, nei)
   nat = 0
   do i = 1,nbr_quantum
      is_at_quantum(i) = .true.
      nat = nat + 1
      pos_temp(i) = posa(i)
      pos_temp(i+natoms) = posa(i+natoms)
      pos_temp(i+natoms+natoms) = posa(i+natoms+natoms)
   enddo

   do i = 1,nbr_quantum
      if (passivate) then 
         do j = 1,numnei(i)
            k = nei(i,j)
            if ( .not. is_at_quantum(k)) then 
               xij = posa(k)-posa(i) - box(1) * nint((posa(k)-posa(i))*invbox(1))
               yij = posa(k+natoms)-posa(i+natoms)
               zij = posa(k+2*natoms)-posa(i+2*natoms) - box(3) * nint((posa(k+2*natoms)-posa(i+2*natoms))*invbox(3))
               rij2 = xij*xij + yij*yij + zij*zij
               if (rij2 .lt. 2.7d0*2.7d0) then
                  nat = nat + 1
                  pos_temp(nat) = posa(i) + 0.5d0*xij
                  pos_temp(nat+natoms) = posa(i+natoms) + 0.5d0*yij
                  pos_temp(nat+natoms+natoms) = posa(i+2*natoms) + 0.5d0*zij !we passivate with hydrogene at this distance
               endif
            endif   
         enddo
      endif
   enddo

   allocate(posquant(3*nat))
   allocate(force_quant(3*nat))
   force_quant = 0.0d0
   posquant(1:nat) = pos_temp(1:nat)
   posquant(1+nat:nat+nat) = pos_temp(1+natoms:nat+natoms)
   posquant(1+nat+nat:nat+nat+nat) = pos_temp(1+natoms+natoms:nat+natoms+natoms)

   energy = 0.d0  !!!thus only the quantum region should contribute to the forces
   !now the box is passivated, we can send it to BigDFT
   call calcforce_bigdft(posquant,force_quant, box, energy, evalf_numbera, conva )

   !we assign the quantum forces
   do i = 1,nbr_quantum - nbr_quantum_trash
      forca(i) = force_quant(i)
      forca(i+natoms) = force_quant(i+nat)
      forca(i+natoms+natoms) = force_quant(i+nat+nat)   
   end do

   !we assign the classical forces
   do i = nbr_quantum - nbr_quantum_trash + 1,natoms
      forca(i) = force_class(i)
      forca(i+natoms) = force_class(i+natoms)
      forca(i+natoms+natoms) = force_class(i+natoms+natoms)
   enddo

   deallocate(posquant,force_quant)

END SUBROUTINE calcforce_mix
