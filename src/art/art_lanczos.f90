!> @file
!!    Contains routines for the lanczos procedure
!!
!! @details
!! Modified by:
!! - EM 2010, see ~/AUTHORS
!! - Laurent Karim Beland, UdeM 2011: Gramm-Schmidt orthogonalization
!! - ID and EM 2011, see ~/AUTHORS: Gramm-Schmidt, DSTEV->DGEEV
!!
!! @copyright
!!    Copyright (C) 2001 Normand Mousseau
!!    Copyright (C) 2010-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Module to use lanczos inside ART
module lanczos_defs

  implicit none
  save

  integer      :: NVECTOR_LANCZOS_H
  integer      :: NVECTOR_LANCZOS_C
  integer      :: LANCZOS_SCL

  real(kind=8) :: DEL_LANCZOS
  real(kind=8) :: eigenvalue
  !real(kind=8), dimension(10) :: eigenvals
  real(kind=8) :: lanc_energy
  real(kind=8) :: proj_energy
  real(kind=8) :: collinear_factor 

  ! Projection direction based on lanczos computations of lowest eigenvalues
  real(kind=8), dimension(:), allocatable :: old_projection, projection 
  logical :: LANCZOS_MIN, IN_MINIMUN

  ! for new projection vector but starting with a random displacement localized 
  ! only around a given atom. This is only for Type_of_events == local or list_local 

END MODULE lanczos_defs


!> ART lanczos
!! Determine lowest frequencies
subroutine lanczos( maxvec, new_projection, produit )

  use defs
  use lanczos_defs
  use random
  use bigdft_forces
  use module_interfaces
  implicit none

  !Arguments
  integer, intent(in)       :: maxvec
  logical, intent(in)       :: new_projection
  real(kind=8), intent(out) :: produit 

  !Local variables
  integer      :: i, j, k, i_err, ivec, ierror
  real(kind=8) :: a1, a0, b2, b1, c1, norm, ran3
  real(kind=8) :: excited_energy, sum2, invsum
  real(kind=8) :: dr2
  real(kind=8), dimension(3) :: boxl
  real(kind=8), dimension(VECSIZE)         :: newpos, newforce, ref_force
  real(kind=8), dimension(maxvec)          :: diag
  real(kind=8), dimension(maxvec-1)        :: offdiag
  real(kind=8), dimension(maxvec, maxvec)  :: vector
  real(kind=8), dimension(VECSIZE, maxvec), target :: lanc
  ! Vectors used to build the matrix for Lanzcos algorithm 
  real(kind=8), dimension(:), pointer :: z0, z1, z2
  real(kind=8), dimension(:), pointer :: dx, dy, dz

  real(kind=8), dimension(maxvec, maxvec)  :: Hmat
  real(kind=8), dimension(:), pointer      :: work
  real(kind=8), dimension(maxvec)          :: e_real
  real(kind=8), dimension(maxvec)          :: e_imag
  real(kind=8), dimension(1)               :: dummy_vl
  integer                                  :: lwork
  integer                                  :: i_min
  real(kind=8)                             :: e_min
  interface
     subroutine center( vector, vecsize )
       integer, intent(in) :: vecsize
       real(kind=8), dimension(vecsize), intent(inout), target :: vector
     end subroutine center
  end interface
  
  !_______________________
  newpos = 0.0d0
  diag = 0.0d0
  offdiag = 0.0d0
  lanc = 0.0d0
  
  boxl = box * scala
                                      ! We now take the current position as the 
                                      ! reference point and will make  a displacement
                                      ! in a random direction or using the previous 
                                      ! direction as the starting point.
  call calcforce( NATOMS, pos, boxl, ref_force, lanc_energy, evalf_number, .true. )

  z0 => lanc(:,1)

  dx => lanc(1:NATOMS,1)
  dy => lanc(NATOMS+1:2*NATOMS,1)
  dz => lanc(2*NATOMS+1:3*NATOMS,1)

  if ( .not. new_projection ) then
     z0 = projection             
     old_projection = projection 
  else
     dx = 0.0d0
     dy = 0.0d0
     dz = 0.0d0
     if ( iproc == 0 ) then           ! Initial movement.
        do i = 1, natoms, 1
           if ( constr(i) == 0 .and. in_system(i) == 0 ) then
              do
                dx(i) = 0.5d0 - ran3()
                dy(i) = 0.5d0 - ran3()
                dz(i) = 0.5d0 - ran3()
                                      ! displacement is isotropic
                dr2 = dx(i)**2 + dy(i)**2 + dz(i)**2
                if ( dr2 < 0.25d0 ) exit 
              end do
           end if
        end do

        call center( z0, VECSIZE )    ! Center of mass at zero.
     end if
                                      ! Broadcast z0 to all notes
     call MPI_Bcast(z0,3*NATOMS,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
     old_projection = 0.0d0
  end if
                                      ! We normalize the displacement
  sum2 = dot_product ( z0, z0 )
  invsum = 1.0d0/sqrt(sum2)
  z0 = z0 * invsum

  newpos = pos + z0 * DEL_LANCZOS     ! New positions.

  call calcforce( NATOMS, newpos, boxl, newforce, excited_energy, evalf_number, .true. )
                                      ! We extract lanczos(1)
  newforce = newforce - ref_force  
  newforce = newforce*(-1.0d0/DEL_LANCZOS)
                                      ! We get a0
  a0 = dot_product( z0, newforce )
  diag(1) = a0

  z1 => lanc(:,2)
  z1 = newforce - a0 * z0  
  b1 = dot_product( z1, z1 )
  offdiag(1) = sqrt(b1)
  invsum = 1.0d0/sqrt(b1)
  z1 = z1 * invsum           
                                      ! We can now repeat this game for 
                                      ! the next vectors.
  do ivec = 2, maxvec - 1

     z1 => lanc(:,ivec)
     newpos = pos + z1 * DEL_LANCZOS 
     call calcforce( NATOMS, newpos, boxl, newforce, excited_energy, evalf_number, .true. )
     newforce = newforce - ref_force
     newforce = newforce*(-1.0d0/DEL_LANCZOS)

     a1 = dot_product( z1, newforce )
     diag(ivec) = a1

     b1 = offdiag(ivec-1)
     z0 => lanc(:,ivec-1)
     z2 => lanc(:,ivec+1)
     z2 = newforce - a1*z1 - b1*z0
     
     do j = 1, 4                      !Gramm-Schmidt orthogonalization
          do i = 1, ivec
             ! write(*,*) i,dot_product(z2,lanc(:,i))
             z2 = z2 - dot_product( z2, lanc(:,i) )*lanc(:,i)
          end do
     
          b2 = dot_product( z2, z2 )
          invsum = 1.0d0/sqrt(b2)
          z2 = z2 * invsum
     end do

     b2 = dot_product( z2, newforce )
     offdiag(ivec) = b2
  end do

  ivec = maxvec                       ! We now consider the last line of
  z1 => lanc(:,maxvec)                ! our matrix  

  newpos = pos + z1 * DEL_LANCZOS  
  call calcforce( NATOMS, newpos, boxl, newforce, excited_energy, evalf_number, .true. )
  newforce = newforce - ref_force
  newforce = newforce*(-1.0d0/DEL_LANCZOS)

  a1 = dot_product( z1, newforce )
  diag(maxvec) = a1

  diag = 1.0d0 * diag
  offdiag = 1.0d0 * offdiag

  ! We now have everything we need in order to diagonalise and find the eigenvectors.
  ! We call the routine for diagonalizing a tridiagonal  matrix.

  i_err = 0                          ! Initialization of lapack infocode

                                     ! old lapack call
  !call dstev( 'V', maxvec, diag, offdiag, vector, maxvec, scratcha, i_err )
  !*******************************************************************
  ! fill the tridiagonal matrix
  Hmat = 0.0d0
  do i = 1, maxvec
     Hmat(i,i) = diag(i)
  end do
  do i= 1, maxvec - 1
     Hmat(i,i+1) = offdiag(i)
     Hmat(i+1,i) = offdiag(i)
  end do

  ! first call to get working dimensions
  lwork = -1
  allocate( work(maxvec) )
  call dgeev( 'N', 'V', maxvec, Hmat, maxvec, e_real, e_imag, dummy_vl, 1, vector, maxvec, work, lwork, i_err)
  lwork = int(work(1))
  deallocate( work )
  allocate( work(lwork) )
  ! second call to get eigen values
  call dgeev( 'N', 'V', maxvec, Hmat, maxvec, e_real, e_imag, dummy_vl, 1, vector, maxvec, work, lwork, i_err)
  deallocate( work ) 
  
  ! look for the smallest eigen value
  e_min = e_real(1)
  i_min = 1
  do i = 1, maxvec
     if ( e_real(i) < e_min ) then 
        e_min = e_real(i)
        i_min = i  
     end if
     if ( abs( e_imag(i) ) > 1.0d-8 ) then
        if ( iproc == 0 ) then
        open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
            & action = 'write', position = 'append', iostat = ierror )
        write( FLOG, * ) 'Error: not supposed to get imaginary value here'
        close( FLOG )
        end if
     end if
  end do  
  
  ! switch eigen vectors (scratch first vector)
  diag(1) = e_min
  do i = 1, maxvec
     vector(i,1) = vector(i,i_min)
  end do
  !*******************************************************************
  ! check error message
  if ( i_err /= 0 ) then             ! If unsuccessful dstev.
     if ( iproc == 0 ) then
        open( unit = FLOG, file = LOGFILE, status = 'unknown',& 
            & action = 'write', position = 'append', iostat = ierror )
        write( FLOG, * ) " WARNING: i_err= ", i_err, " ->failed dstev in lanczos" 
        close( FLOG ) 
        write( *, * ) "WARNING: i_err= ", i_err, " ->failed dstev in lanczos"
     end if
  end if

  ! We now reconstruct the eigenvectors in the real space, of course, we need only the 
  ! first maxvec elements of vec.
  projection = 0.0d0 
  do k = 1, maxvec
     z1 => lanc(:,k)
     a1 = vector(k,1)
     projection = projection + a1 * z1  
  end do

  c1 = dot_product( projection, projection )
  norm = 1.0d0/sqrt(c1)
  projection = projection * norm 
  eigenvalue = diag(1)
  
  a1 = dot_product( old_projection, projection )

  if ( a1 < 0.0d0 ) then 
     projection = -1.0d0 * projection 
     produit = -1.0d0 * a1
  else
     produit = a1
  end if
                                      ! Check: center of mass at zero. 
  call center ( projection, VECSIZE )
                                      ! This is for analyzing the minimum. 
  if ( IN_MINIMUN .and. LANCZOS_MIN ) then 
     newpos = pos + projection * DEL_LANCZOS 
     call calcforce( NATOMS, newpos, boxl, newforce, proj_energy, evalf_number, .true. )
  else
     proj_energy = lanc_energy 
  end if
                                      ! Broadcast eigenvalue and projection to all nodes.
  call MPI_Bcast(eigenvalue,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
  call MPI_Bcast(projection,3*NATOMS,MPI_REAL8,0,MPI_COMM_WORLD,ierror)

END SUBROUTINE lanczos


!> ART center
!! @author
!! Written by EM 2010, see ~/AUTHORS 
!! It places the center of mass of a 3D vector at (0,0,0). 
subroutine center( vector, vecsize )

  use defs, only : natoms, constr
  use bigdft_forces, only : in_system
  implicit none

  !Arguments
  integer, intent(in) :: vecsize
  real(kind=8), dimension(vecsize), intent(inout), target :: vector

  !Local variables
  integer      :: i, natoms_f  
  real(kind=8) :: xtotal, ytotal, ztotal
  real(kind=8), dimension(:), pointer :: x, y, z     ! Pointers for coordinates
  logical, dimension(natoms) :: mask
  !_______________________

  ! degrees of freedom 
  mask = constr .eq. 0 .and. in_system .eq. 0 
  natoms_f = count(mask)

  ! We first set-up pointers for the x, y, z components 
  x => vector(1:natoms)
  y => vector(natoms+1:2*natoms)
  z => vector(2*natoms+1:3*natoms)

  xtotal = 0.0d0
  ytotal = 0.0d0
  ztotal = 0.0d0

  do i = 1, natoms
     xtotal = xtotal + x(i)
     ytotal = ytotal + y(i)
     ztotal = ztotal + z(i)
  enddo 

  xtotal = xtotal / natoms_f 
  ytotal = ytotal / natoms_f
  ztotal = ztotal / natoms_f

  do i = 1, natoms
     if ( mask(i) ) then
        x(i) = x(i) - xtotal
        y(i) = y(i) - ytotal
        z(i) = z(i) - ztotal
     end if 
  end do

END SUBROUTINE center

