!> Module to define DIIS inside ART
!!
!! @author
!!    Copyright (C) 2001 Normand Mousseau
!!    Copyright (C) 2010 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!! 
module diis_defs
  use defs
  implicit none
  real(8), dimension(:,:), allocatable  :: previous_forces
  real(8), dimension(:,:), allocatable  :: previous_pos
  real(8), dimension(:,:), allocatable :: product_matrix
  real(8), dimension(:), allocatable :: tildeforce
end module diis_defs



!> Apply DIIS for ART
!! 
subroutine apply_diis(current_energy, ftot)
  use defs
  use diis_defs
  use random
  use bigdft_forces
  implicit none

  real(8), intent(out) :: current_energy
  real(8), intent(out) :: ftot
  integer :: npart,lter, maxter, ierror
  real(8) :: boxl, delr
  real(8), dimension(VECSIZE) :: posb

  boxl = box * scala
  posb = pos

  allocate(previous_forces(DIIS_MEMORY,VECSIZE))
  allocate(previous_pos(DIIS_MEMORY,VECSIZE))
  allocate(product_matrix(DIIS_MEMORY,DIIS_MEMORY))
  allocate(tildeforce(VECSIZE))

  call calcforce(NATOMS,pos,boxl,tildeforce,total_energy)
  evalf_number = evalf_number + 1

  ! We set the first step and move to the second
  previous_forces(1,:) = tildeforce(:)
  previous_pos(1,:) = pos(:)
  pos(:) = pos(:) + DIIS_STEP * tildeforce(:)
  do lter = 2, DIIS_MAXITER
     maxter = min(lter,DIIS_MEMORY)
     call diis(lter,maxter,pos)
     ftot = dsqrt(dot_product(tildeforce,tildeforce))

     if (iproc .eq. 0 ) then 
        call displacement(posb, pos, delr,npart)
        
        open(unit=FLOG,file=LOGFILE,status='unknown',action='write',position='append',iostat=ierror)
        write(FLOG,"(' ','lter: ',i4,'  ener: ',f12.4,'  ftot: ',f12.6,&
             & '  delr: ',f10.4,'  npart: ',i4,' evalf: ',i7)")  &
             & lter, total_energy, ftot, delr, npart,evalf_number
        write(*,"(' ','lter: ',i4,'  ener: ',f12.4,'  ftot: ',f12.6,&
             & '  delr: ',f10.4,'  npart: ',i4,' evalf: ',i7)")  &
             & lter, total_energy, ftot, delr, npart,evalf_number
        close(unit=FLOG)
     endif
     if (ftot .lt. DIIS_FORCE_THRESHOLD) exit
  end do

  current_energy = total_energy
  deallocate(previous_forces)
  deallocate(previous_pos)
  deallocate(product_matrix)
  deallocate(tildeforce)

  return
END SUBROUTINE apply_diis



!>   This program implement the direct inversion in iterative subspace
!!   method which allows one to converge rapidly to a saddle point when
!!   we are in its vicinity.
!!
!!   maxvec is the iteration number. The matrix computed has one more dimension
!! 
subroutine diis(lter,maxvec, newpos)
  use defs
  use diis_defs
  use random
  use bigdft_forces
  implicit none

  integer, intent(in) :: lter, maxvec
  real(8), intent(out), dimension(VECSIZE) :: newpos
  real(8), dimension(VECSIZE) :: tildepos

  ! Variables for lapack routine
  integer, dimension(maxvec+1) :: interchanges
  real(8), dimension((maxvec+1)*(maxvec+1)) :: work
  real(8), dimension(maxvec+1) :: solution
  real(8), dimension(maxvec+1, maxvec+1) :: matrice
 
  integer :: i,j,i_err,lwork,n,nrhs
  real(8) :: boxl

  boxl = box * scala

  ! If lter is greater than maxvec, we move the previous solution up by one
  if (lter .gt. maxvec) then
     do i=2, maxvec
        previous_forces(i-1,:) = previous_forces(i,:)
        previous_pos(i-1,:)    = previous_pos(i,:)
        do j=2, maxvec
           product_matrix(i-1,j-1) = product_matrix(i,j)
        end do
     end do
  end if

  ! we first add the force to the previous_force vector and the
  ! position to the previous_pos vector
  previous_forces(maxvec,:) = tildeforce(:)
  previous_pos(maxvec,:) = pos(:)
  

  ! And we add the scalar products to the matrix
  do i = 1, maxvec
     product_matrix(i,maxvec) = dot_product(previous_forces(i,:),tildeforce(:))
     product_matrix(maxvec,i) = product_matrix(i,maxvec)
  end do

  matrice(1:maxvec,1:maxvec) = product_matrix(1:maxvec, 1:maxvec)
  matrice(1:maxvec,maxvec+1) = 1.0d0
  matrice(maxvec+1,1:maxvec) = 1.0d0
  matrice(maxvec+1,maxvec+1) = 0.0d0

  solution(1:maxvec) = 0.0d0
  solution(maxvec+1) = 1.0d0
  
  ! We now need the routines from Lapack. We define a few values
  i_err = 0

  ! We call the routine for diagonalizing a tridiagonal  matrix
  n = maxvec+1
  nrhs = 1    ! Number of right-hand arguments in B (Ax=B)
  lwork = n*n

  ! We prepare the upper triangular matrix for lapack
  call dsysv('U',n, nrhs, matrice, n, interchanges, solution,n,work,lwork,i_err)


  ! The solution that interests us is made of two parts

  tildepos = 0.0d0
  do i = 1, maxvec
     tildepos(:) = tildepos(:) + solution(i) * previous_pos(i,:)
  end do

  call calcforce(NATOMS,tildepos,boxl,tildeforce,total_energy)
  evalf_number = evalf_number + 1
  newpos(:) = tildepos(:) + DIIS_STEP * tildeforce(:)

END SUBROUTINE diis

