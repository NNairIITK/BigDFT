!!****m* art/lanczos_defs
!! FUNCTION
!!   Module to use lanczos inside art
!!
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
module lanczos_defs

  implicit none

  logical :: first_time = .true., reject= .false.
  real(8) :: eigenvalue
  real(8), dimension(10) :: eigenvals

  ! Projection direction based on lanczos computations of lowest eigenvalues
  real(8), dimension(:), allocatable :: old_projection, projection 
  real(8) :: INC_LANCZOS
  logical :: LANCZOS_MIN
  integer :: NVECTOR_LANCZOS

end module lanczos_defs
!!***


!!****f* art_lanczos/lanczos
!! FUNCTION
!!   Lanczos routine to determine lowest frequencies
!! SOURCE
!!
subroutine lanczos( maxvec, new_projection, produit )

  use defs
  use lanczos_defs
  use random
  use bigdft_forces
  implicit none

  !Arguments
  integer, intent(in) :: maxvec
  logical, intent(in) ::  new_projection
  real(8), intent(out) :: produit 

  !Local variables
  real(8), dimension( 2 * maxvec -1 ) :: scratcha
  real(8), dimension(maxvec) :: diag
  real(8), dimension(maxvec-1) :: offdiag
  real(8), dimension(maxvec, maxvec) :: vector
  real(8), dimension(VECSIZE, maxvec), target :: lanc
  ! Vectors used to build the matrix for Lanzcos algorithm 
  real(8), dimension(:), pointer :: z0, z1, z2
 
  integer :: i, k, i_err, ivec, ierror, nat
 
  integer :: ivec2
  real(8), dimension(maxvec,maxvec):: test

  real(8) :: a1, a0,b2,b1
  real(8) :: boxl, excited_energy,c1,norm
  real(8) :: xsum, ysum, zsum, sum2, invsum
  real(8), dimension(VECSIZE) :: newpos,newforce,ref_force
  real(8), dimension(VECSIZE) :: forcea, forceb
  real(8) :: ran3

  boxl = box * scala

                                      ! We now take the current position as the 
                                      ! reference point and will make  a displacement
                                      ! in a random direction or using the previous 
                                      ! direction as the starting point.

  call calcforce( NATOMS, pos, boxl, ref_force, total_energy )
  evalf_number = evalf_number + 1

  z0 => lanc(:,1)

  if ( .not. new_projection ) then
     z0 = projection             
     old_projection = projection 
  else
                                      ! Initial movement.
     if ( iproc == 0 ) then 
       do i = 1, VECSIZE
          z0(i) = 0.5d0 - ran3()
       end do
                                      ! Center of mass at zero. 
       call center( z0, VECSIZE )
     end if
       
                                      ! Broadcast z0 to all notes
    nat = 3 * NATOMS
    call MPI_Bcast(z0,nat,MPI_REAL8,0,MPI_COMM_WORLD,ierror)

    if ( first_time ) then
       old_projection = z0 
       first_time = .false.
    else
       old_projection = projection
    end if 

  end if
                                      ! We normalize the displacement to  
                                      ! 1 total.
  sum2 = 0.0d0
  sum2 = dot_product ( z0, z0 )
  invsum = 1.0/sqrt(sum2)
  z0 = z0 * invsum
                                      ! New positions.
  newpos = pos + z0 * INC_LANCZOS

  call calcforce( NATOMS, newpos, boxl, newforce, excited_energy )
  evalf_number = evalf_number + 1
  
                                      ! We extract lanczos(1)
  newforce = newforce - ref_force  
  newforce = newforce*(-1.0d0/INC_LANCZOS)

                                      ! We get a0
  a0 = 0.0d0
  a0 = dot_product( z0, newforce )

  diag(1) = a0

  z1 => lanc(:,2)
  z1 = newforce - a0 * z0  

  b1 = 0.0d0
  b1 = dot_product( z1, z1 )

  offdiag(1) = sqrt(b1)

  invsum = 1.0d0 / sqrt(b1)
  z1 = z1 * invsum           

                                      ! We can now repeat this game for 
                                      ! the next vectors.
  do ivec = 2, maxvec-1
     z1 => lanc(:,ivec)

     newpos = pos + z1 * INC_LANCZOS 
     call calcforce( NATOMS, newpos, boxl, newforce, excited_energy )
     evalf_number = evalf_number + 1
     newforce = newforce - ref_force
     newforce = newforce*(-1.0d0/INC_LANCZOS)

     a1 = 0.0d0
     a1 = dot_product( z1, newforce )
     diag(ivec) = a1

     b1 = offdiag(ivec-1)
     z0 => lanc(:,ivec-1)
     z2 => lanc(:,ivec+1)
     z2 = newforce - a1*z1 - b1*z0

     b2=0.0d0
     b2 = dot_product( z2, z2 )

     offdiag(ivec) = sqrt(b2)
     
     invsum = 1.0/sqrt(b2)
     z2 = z2 * invsum

  end do
                                      ! We now consider the last line of
                                      ! our matrix

  ivec = maxvec
  z1 => lanc(:,maxvec)

  newpos = pos + z1 * INC_LANCZOS  

  call calcforce( NATOMS, newpos, boxl, newforce, excited_energy )
  evalf_number = evalf_number + 1

  newforce = newforce - ref_force
  newforce = newforce*(-1.0d0/INC_LANCZOS)
  a1 = 0.0d0
  a1 = dot_product( z1, newforce )
  diag(maxvec) = a1
                                      ! We now have everything we need in order to
                                      ! diagonalise and find the eigenvectors.

                                      ! We now need the routines from Lapack.
                                      ! We define a few values.
   i_err = 0
                                      ! We call the routine for diagonalizing a
                                      ! tridiagonal  matrix.
   call dstev( 'V', maxvec, diag, offdiag, vector, maxvec, scratcha, i_err )
                                      ! If unsuccessful dstev.
   if ( i_err /= 0 ) then
      if ( iproc == 0 ) then
       open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
           & action = 'write', position = 'append', iostat = ierror )
       write( FLOG, * ) " WARNING: i_err= ", i_err, " ->failed dstev in lanczos" 
       close( FLOG ) 
       write( *, * ) "WARNING: i_err= ", i_err, " ->failed dstev in lanczos"
      end if
   end if
                                      ! We now reconstruct the eigenvectors in the
                                      ! real space, of course, we need only the 
                                      ! first maxvec elements of vec.
   projection = 0.0d0 
   do k = 1, maxvec
      z1 => lanc(:,k)
      a1 = vector(k,1)
      projection = projection + a1 * z1  
   end do

   c1 = 0.0d0
   c1 = dot_product( projection, projection )
   norm = 1.0/sqrt(c1)
   projection = projection * norm 
   
   eigenvalue = diag(1)
   do i = 1, 4
      eigenvals(i) = diag(i) 
   end do
   
   a1 = 0.0d0
   a1 =  dot_product( old_projection, projection )
                                      ! Condition on the scalar product: we reject the
                                      ! point where we loose the eigenvalue and
                                      ! try to reduce the step size.

   if ( a1 < 0.0d0 ) then 
      projection = -1.0d0 * projection 
      produit = -1.0d0 * a1
   else
      produit = a1
   end if
                                      ! Write 
   if ( iproc == 0 ) then
    write(*,"(' ','BART: the scalar product between the two eigen directions: ',f12.6)") ,a1
    write(*,*) "BART: lanczos EIGENVALUE= ", eigenvalue
   end if

                                      ! Center of mass at zero. 
   call center ( projection, VECSIZE )
                                      ! Broadcast eigenvalue and projection to all nodes.
   nat = 3 * NATOMS
   call MPI_Bcast(eigenvalue,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
   call MPI_Bcast(projection,nat,MPI_REAL8,0,MPI_COMM_WORLD,ierror)

END SUBROUTINE lanczos
!!***
