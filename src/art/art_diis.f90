!!****m* art/diis_def
!! COPYRIGHT
!!    Copyright (C) 2001 Normand Mousseau
!!    Copyright (C) 2010 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! SOURCE
!! 
module diis_defs

  implicit none

  real(8), dimension(:,:), allocatable :: previous_forces
  real(8), dimension(:,:), allocatable :: previous_pos
  real(8), dimension(:,:), allocatable :: product_matrix
  real(8), dimension(:),   allocatable :: tildeforce

end module diis_defs
!!***


!!****f* art/apply_diis
!! FUNCTION
!! SOURCE
!! 
subroutine apply_diis( current_energy, ftot, delr, npart )

  use defs
  use diis_defs
  use bigdft_forces
  use lanczos_defs,   Only: projection
  implicit none

  !Arguments
  real(8), intent(out) :: current_energy
  real(8), intent(out) :: ftot
  real(8), intent(out) :: delr
  integer, intent(out) :: npart

  !Local variables
  integer :: lter, maxter, ierror
  real(8) :: boxl
  real(8) :: delta_e                  ! current_energy - ref_energy


  boxl = box * scala

  allocate(previous_forces(DIIS_MEMORY,VECSIZE))
  allocate(previous_pos(DIIS_MEMORY,VECSIZE))
  allocate(product_matrix(DIIS_MEMORY,DIIS_MEMORY))
  allocate(tildeforce(VECSIZE))

  ! Initialization.
  previous_forces = 0.0d0  
  previous_pos    = 0.0d0
  product_matrix  = 0.0d0

  ! The forces of our input configuration will be the first
  ! residual vector.
  call calcforce(NATOMS,pos,boxl,tildeforce,total_energy)
  evalf_number = evalf_number + 1

  ! We set the first step.
  previous_forces(1,:) = tildeforce(:)
  previous_pos(1,:)    = pos(:)

  ! We set the first matrix element.
  product_matrix(1,1) = dot_product(previous_forces(1,:),previous_forces(1,:))

  ! We move to the second step, i.e, this is the second trial step.
  pos(:) = pos(:) + DIIS_STEP * projection(:)

  ! New force is evaluted. It will be stored as 
  ! previous_forces(2,:) at the beggining of next loop.

  call calcforce(NATOMS,pos,boxl,tildeforce,total_energy)
  evalf_number = evalf_number + 1

  do lter = 2, DIIS_MAXITER

     maxter = min(lter,DIIS_MEMORY)
     call diis( lter, maxter, pos )
     ftot = dsqrt(dot_product(tildeforce,tildeforce))
     call displacement( posref, pos, delr, npart )
     delta_e = total_energy - ref_energy

     if ( SAVE_CONF_INT ) call save_intermediate( 'L', lter ) 

                                      ! Write
     if (iproc == 0 ) then 
      open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
          & action = 'write', position = 'append', iostat = ierror )
      write(FLOG,'(a6,i4,x,F11.4,I3,I2,4f12.4,x,1f10.3,i4,i6,f7.2)') &
      & "lter=", lter, delta_e, 0, 0, ftot, 0.0, 0.0, 0.0,     &
      &  delr, npart, evalf_number, 0.0
      close(FLOG)
      write(*,'(a,i4,x,(1p,e17.10,0p),i3,i2,4f12.6,x,1f10.4,i4,i6)') &
      & " BART: lter: ", lter, total_energy, 0, 0, ftot,   & 
      & 0.0, 0.0, 0.0, delr, npart, evalf_number
     end if

     if ( ftot .lt. DIIS_FORCE_THRESHOLD ) exit

  end do

  current_energy = total_energy
  force          = tildeforce

  deallocate(previous_forces)
  deallocate(previous_pos)
  deallocate(product_matrix)
  deallocate(tildeforce)

END SUBROUTINE apply_diis
!!***


!!****f* art/diis
!! FUNCTION
!!   This program implement the direct inversion in iterative subspace
!!   method which allows one to converge rapidly to a saddle point when
!!   we are in its vicinity.
!!
!!   maxvec is the iteration number. The matrix computed has one more dimension
!!
!! SOURCE
!! 
subroutine diis( lter, maxvec, newpos )

  use defs
  use diis_defs
  use bigdft_forces
  implicit none

  !Arguments
  integer, intent(in) :: lter
  integer, intent(in) :: maxvec
  real(8), intent(out), dimension(VECSIZE) :: newpos

  !Local variables
  real(8), dimension(VECSIZE) :: tildepos

  ! Variables for lapack routine
  integer, dimension(maxvec+1) :: interchanges
  real(8), dimension(:), allocatable :: work
  real(8), dimension(maxvec+1) :: solution
  real(8), dimension(maxvec+1, maxvec+1) :: matrice
 
  integer :: i,j,i_err,lwork,n,nrhs
  real(8) :: boxl

  boxl = box * scala
 
  ! Initialization.
  solution = 0.0d0
  matrice  = 0.0d0

  ! If lter is greater than maxvec, we move the previous solution up by one
  if ( lter > maxvec ) then
     do i = 2, maxvec
        previous_forces(i-1,:) = previous_forces(i,:)
        previous_pos(i-1,:)    = previous_pos(i,:)
        do j = 2, maxvec
           product_matrix(i-1,j-1) = product_matrix(i,j)
        end do
     end do
  end if

  ! We first add the force to the previous_force vector and the
  ! position to the previous_pos vector
  previous_forces(maxvec,:) = tildeforce(:)
  previous_pos(maxvec,:)    = pos(:)

  ! And we add the scalar products to the matrix
  do i = 1, maxvec
     product_matrix(i,maxvec) = dot_product(previous_forces(i,:),tildeforce(:))
     product_matrix(maxvec,i) = product_matrix(i,maxvec)
  end do

  ! This is our matrix of overlaps.
  matrice(1:maxvec,1:maxvec) = product_matrix(1:maxvec,1:maxvec)
  matrice(1:maxvec,maxvec+1) = -1.0d0
  matrice(maxvec+1,1:maxvec) = -1.0d0
  matrice(maxvec+1,maxvec+1) =  0.0d0

  solution(1:maxvec) =  0.0d0
  solution(maxvec+1) = -1.0d0
  
  ! We now need the routines from Lapack. We define a few values
  i_err = 0

  ! We call the routine for diagonalizing a tridiagonal  matrix
  n = maxvec+1
  nrhs = 1    ! Number of right-hand arguments in B (Ax=B)

  ! We find the best size for "work" array.
  allocate(work(100))
  call dsysv('U',n, nrhs, matrice, n, interchanges, solution,n,work,-1,i_err)
  lwork=work(1)
  deallocate(work)
  allocate(work(lwork))

  ! We prepare the upper triangular matrix for lapack
  call dsysv('U',n, nrhs, matrice, n, interchanges, solution,n,work,lwork,i_err)

  if ( i_err /= 0 ) then
     if ( iproc == 0 ) write(*,*) 'WARNING DIIS: info calculation of solution', i_err
  end if

  ! The solution that interests us is made of two parts.
  ! This is the trial configuration. 
  tildepos = 0.0d0
  do i = 1, maxvec
     tildepos(:) = tildepos(:) + solution(i) * previous_pos(i,:)
  end do

  ! This is the residuum vector
  tildeforce = 0.0d0
  do i = 1, maxvec
     tildeforce(:) = tildeforce(:) + solution(i) * previous_forces(i,:)
  end do
  tildeforce(:) =  DIIS_STEP * tildeforce(:)

  ! Our new positions.
  newpos(:) = tildepos(:) +  tildeforce(:)

  ! Now the force is evaluted for our newpos configuration. 
  call calcforce(NATOMS,newpos,boxl,tildeforce,total_energy)
  evalf_number = evalf_number + 1

END SUBROUTINE diis
!!***
