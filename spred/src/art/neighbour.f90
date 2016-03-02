!> @file
!! This subroutine should depend only on the details of the physics and not
!! on the ART algorithm per se. It should be called only from the force or
!! energy routine
!! 
!! @author
!!   Written by Laurent Karim Beland, UdeM 2011!!
!!   Copyright (C) 2010-2011 BigDFT group, Normand Mousseau
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 


subroutine neighbours(natoms,pos,boxl,boundary,maxnei,numnei, nei)
   use cutoff
   implicit none

   !Arguments
   integer, intent(in)                     :: natoms
   real(8), dimension(natoms),intent(in)   :: pos
   real(8), dimension(idim), intent(in)                     :: boxl
   character(len=1),intent(in)              :: boundary
   integer, intent(in)                      :: maxnei
   integer, dimension(natoms), intent(out) :: numnei
   integer, dimension(natoms,maxnei), intent(out) :: nei
   !Local variables
   real(8), dimension(idim)    :: box
   integer i, i_ind, j, j_ind
   real(8) :: xi, yi, zi, xij, yij, zij, rij2
   real(8), dimension(idim) :: invbox

   box = boxl
   invbox = 1.0d0/box

   numnei = 0  ! Vectorial assignment
   nei = 0  

   do i=1, NATOMS

      !i_ind = types(i)
      i_ind = 1
      xi = pos(i)
      yi = pos(i+natoms)
      zi = pos(i+2*natoms)

      do j=i+1, NATOMS
         !        j_ind = types(j)
         j_ind = 1

         !be carefull with boundaries if surface
         if ( boundary .eq. 'P') xij = pos(j) - xi - box(1) * nint((pos(j)-xi)*invbox(1))
         if ( boundary .eq. 'P') yij = pos(j+natoms) - yi - box(2) * nint((pos(j+natoms)-yi)*invbox(2))
         if ( boundary .eq. 'P') zij = pos(j+2*natoms) - zi - box(3) * nint((pos(j+2*natoms)-zi)*invbox(3))

         if ( boundary .eq. 'S') xij = pos(j) - xi - box(1) * nint((pos(j)-xi)*invbox(1))
         if ( boundary .eq. 'S') zij = pos(j+2*natoms) - zi - box(3) * nint((pos(j+2*natoms)-zi)*invbox(3))
         if ( boundary .eq. 'S') yij = pos(j+natoms) - yi


         rij2 = xij*xij + yij*yij + zij*zij

         if (rij2 < rcut2(i_ind,j_ind)*1.1d0) then
            numnei(i) = numnei(i) + 1
            numnei(j) = numnei(j) + 1
            nei(i,numnei(i)) = j
            nei(j,numnei(j)) = i
         endif
      end do
   end do



END SUBROUTINE neighbours


subroutine neighbour_angles(numnei,nei)

   use defs
   use SWpotential

   implicit none

   integer,dimension(natoms),intent(in) :: numnei
   integer,dimension(natoms,maxnei),intent(in) :: nei

   real(8), dimension(3*natoms), target :: pos_normalised
   real(8), dimension(3*natoms) :: box_vec

   real(8), dimension(:), pointer :: xa, ya, za

   real(8), dimension(idim) :: invbox
   integer :: i, j, i_id, ind_j, k, ind_k

   real(8) :: xi, yi, zi, xij, yij, zij, rij2

   real(8) :: xik,yik,zik,rik2

   invbox = 1.0d0/box
   optimum_angle = 1.0d0/3.0d0
   angles_to_fit = 0
   jk_of_angle = 0

   box_vec(1:natoms) = box(1)
   box_vec(1+natoms:2*natoms) = box(2)
   box_vec(1+natoms+natoms:3*natoms) = box(3)

   pos_normalised = pos / box_vec

   xa => pos_normalised(1:NATOMS)
   ya => pos_normalised(NATOMS+1:2*NATOMS)
   za => pos_normalised(2*NATOMS+1:3*NATOMS)

   do i=1, NATOMS
      !     i_id = types(i)
      i_id = 1
      xi = xa(i)
      yi = ya(i)
      zi = za(i)
      do ind_j=1, numnei(i)
         j = nei(i,ind_j)

         ! Pair interactions of i and j
         ! Distance, with periodic boundary conditions
         if (boundary == "P") then
            xij = xa(j) - xi - 1.0d0 * nint( xa(j)-xi )
            yij = ya(j) - yi - 1.0d0 * nint( ya(j)-yi )
            zij = za(j) - zi - 1.0d0 * nint( za(j)-zi )
            elseif (boundary == "S") then
            xij = xa(j) - xi - 1.0d0 * nint( xa(j)-xi )
            yij = ya(j) - yi
            zij = za(j) - zi - 1.0d0 * nint( za(j)-zi )
         endif


         ! Rescale the lengths into Angstroems
         xij = xij * box(1)
         yij = yij * box(2)
         zij = zij * box(3)

         rij2 = xij*xij + yij*yij + zij*zij

         ! Check the cut-off before proceeding
         if( rij2 < 2.8d0*2.8d0 ) then

            do ind_k = ind_j+1, numnei(i)
               ! Triplet interaction with i in the middle; all interactions
               k = nei(i,ind_k)
               !k_id = types(k)

               ! Distance, with periodic boundary conditions
               if (boundary == "P") then
                  xik = xa(k) - xi - 1.0d0 * nint( xa(k)-xi )
                  yik = ya(k) - yi - 1.0d0 * nint( ya(k)-yi )
                  zik = za(k) - zi - 1.0d0 * nint( za(k)-zi )
                  elseif (boundary == "S") then
                  xik = xa(k) - xi - 1.0d0 * nint( xa(k)-xi )
                  yik = ya(k) - yi
                  zik = za(k) - zi - 1.0d0 * nint( za(k)-zi )
               endif

               ! Rescale the lengths into Angstroems
               xik = xik * box(1)
               yik = yik * box(2)
               zik = zik * box(3)

               rik2 = xik*xik + yik*yik + zik*zik

               ! Check whether the distance is too large 
               if (rik2<2.8d0*2.8d0)  then
                  angles_to_fit(i) = angles_to_fit(i)+1
                  jk_of_angle(i,angles_to_fit(i),1) = (j)
                  jk_of_angle(i,angles_to_fit(i),2) = (k)
               endif
            end do
         endif
      end do
   end do


END SUBROUTINE neighbour_angles

subroutine neighbour_bayes(natoms,pos,boxl,boundary,maxnei,numnei, nei)
   use cutoff
   implicit none

   integer, intent(in)                     :: natoms
   real(8), dimension(natoms),intent(in)   :: pos
   real(8), dimension(idim), intent(in)                     :: boxl
   character(len=1),intent(in)              :: boundary
   integer, intent(in)                      :: maxnei
   integer, dimension(natoms), intent(out) :: numnei
   integer, dimension(natoms,maxnei), intent(out) :: nei


   real(8), dimension(idim)    :: box
   integer i, i_ind, j, j_ind
   real(8) :: xi, yi, zi, xij, yij, zij, rij2
   real(8), dimension(idim) :: invbox

   box = boxl
   invbox = 1.0d0/box

   numnei = 0  ! Vectorial assignment
   nei = 0

   do i=1, NATOMS

      !i_ind = types(i)
      i_ind = 1
      xi = pos(i)
      yi = pos(i+natoms)
      zi = pos(i+2*natoms)

      do j=i+1, NATOMS
         !        j_ind = types(j)
         j_ind = 1

         !be carefull with boundaries if surface
         if ( boundary .eq. 'P') xij = pos(j) - xi - box(1) * nint((pos(j)-xi)*invbox(1))
         if ( boundary .eq. 'P') yij = pos(j+natoms) - yi - box(2) * nint((pos(j+natoms)-yi)*invbox(2))
         if ( boundary .eq. 'P') zij = pos(j+2*natoms) - zi - box(3) * nint((pos(j+2*natoms)-zi)*invbox(3))

         if ( boundary .eq. 'S') xij = pos(j) - xi - box(1) * nint((pos(j)-xi)*invbox(1))
         if ( boundary .eq. 'S') zij = pos(j+2*natoms) - zi - box(3) * nint((pos(j+2*natoms)-zi)*invbox(3))
         if ( boundary .eq. 'S') yij = pos(j+natoms) - yi


         rij2 = xij*xij + yij*yij + zij*zij

         if (rij2 < 2.8d0*2.8d0) then
            ! if (rij2 < 2.35d0*2.35d0*1.02d0) then
            numnei(i) = numnei(i) + 1
            numnei(j) = numnei(j) + 1
            nei(i,numnei(i)) = j
            nei(j,numnei(j)) = i
         endif
      end do
   end do

END SUBROUTINE neighbour_bayes
