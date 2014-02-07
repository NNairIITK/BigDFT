!> @file
!! Module for NEB (Nudged-Elastic Band) calculation
!! @author
!!    Copyright (C) 2013-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Modules which contains minimizaton routines for NEB calculation
MODULE Minimization_routines
  use module_defs
    
  IMPLICIT NONE
  
  REAL (gp), PARAMETER, private :: epsi = 1.0D-16

  CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!! minimization algorithms !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! steepest descent with line minimization !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE steepest_descent( ndim, pos, grad, ds )

      IMPLICIT NONE

      integer, intent(in) :: ndim
      real(gp), intent(in) :: ds
      real(gp), dimension(ndim), intent(in) :: grad
      real(gp), dimension(ndim), intent(inout) :: pos

      IF ( nrm2(ndim,grad(1),1) >= epsi ) THEN
         call axpy(ndim, ds, grad(1), 1, pos(1), 1)
      END IF

    END SUBROUTINE steepest_descent 

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! conjugate gradient minimization !!
    !! Fletcher - Reeves               !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE fletcher_reeves( ndim, pos, grad, old_grad, delta, ds )
     
      IMPLICIT NONE

      integer, intent(in) :: ndim
      real(gp), intent(in) :: ds
      real(gp), dimension(ndim), intent(in) :: grad, old_grad
      real(gp), dimension(ndim), intent(inout) :: pos, delta

      REAL (gp) :: gamma, norm_grad
      REAL (gp) :: squared_old_grad_i, abs_conj_dir_i 
      REAL (gp), DIMENSION(:), ALLOCATABLE   :: conj_dir_i

      squared_old_grad_i = DOT_PRODUCT( old_grad, old_grad ) 
      IF ( squared_old_grad_i >= epsi ) THEN
         gamma = DOT_PRODUCT( grad, grad ) / squared_old_grad_i
      ELSE
         gamma = 0.D0
      END IF

      allocate(conj_dir_i(ndim))

      norm_grad = nrm2(ndim,grad(1),1)
      IF ( norm_grad >= epsi ) THEN
        conj_dir_i = - grad / norm_grad + gamma * delta
      ELSE
        conj_dir_i = 0.D0
      END IF

      abs_conj_dir_i = nrm2(ndim,conj_dir_i(1),1)

      IF ( abs_conj_dir_i >= epsi ) THEN
        delta(:) = conj_dir_i / abs_conj_dir_i
      ELSE
        delta(:) = 0.D0
      END IF

      call axpy(ndim, ds * norm_grad, conj_dir_i(1), 1, pos(1), 1)

      deallocate(conj_dir_i)
    END SUBROUTINE fletcher_reeves  

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! conjugate gradient minimization !!
    !! Polak - Ribiere                 !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE polak_ribiere( ndim, pos, grad, old_grad, delta, ds )
     
      IMPLICIT NONE

      integer, intent(in) :: ndim
      real(gp), intent(in) :: ds
      real(gp), dimension(ndim), intent(in) :: grad, old_grad
      real(gp), dimension(ndim), intent(inout) :: pos, delta

      REAL (gp) :: gamma, norm_grad
      REAL (gp) :: squared_old_grad_i, abs_conj_dir_i  
      REAL (gp), DIMENSION(:), ALLOCATABLE   :: conj_dir_i

      squared_old_grad_i = DOT_PRODUCT( old_grad , old_grad ) 
      IF ( squared_old_grad_i >= epsi ) THEN
        gamma = DOT_PRODUCT( ( grad - old_grad ) , grad ) / squared_old_grad_i
      ELSE
        gamma = 0.D0
      END IF

      allocate(conj_dir_i(ndim))

      norm_grad = nrm2(ndim,grad(1),1)
      IF ( norm_grad >= epsi ) THEN
        conj_dir_i = - grad / norm_grad + gamma * delta
      ELSE
        conj_dir_i = 0.D0
      END IF

      abs_conj_dir_i = nrm2(ndim,conj_dir_i(1),1)

      IF ( abs_conj_dir_i >= epsi ) THEN
        delta(:) = conj_dir_i / abs_conj_dir_i
      ELSE
        delta(:) = 0.D0
      END IF

      call axpy(ndim, ds * norm_grad, conj_dir_i(1), 1, pos(1), 1)

      deallocate(conj_dir_i)
    END SUBROUTINE  polak_ribiere

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Molecular Dynamics based algorithms !!
    !! velocity Verlet and quick min       !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE velocity_Verlet_first_step( ndim, pos, vel, grad, ds )

      IMPLICIT NONE

      integer, intent(in) :: ndim
      real(gp), intent(in) :: ds
      real(gp), dimension(ndim), intent(in) :: grad
      real(gp), dimension(ndim), intent(inout) :: pos, vel      
         
      call axpy(ndim, -ds / 2.D0, grad(1), 1, vel(1), 1)
      call axpy(ndim, ds, vel(1), 1, pos(1), 1)
      
    END SUBROUTINE velocity_Verlet_first_step
    
    
    SUBROUTINE velocity_Verlet_second_step( ndim, vel, grad, ds, damp )

      IMPLICIT NONE

      integer, intent(in) :: ndim
      real(gp), intent(in) :: ds
      real(gp), intent(in), optional :: damp
      real(gp), dimension(ndim), intent(in) :: grad
      real(gp), dimension(ndim), intent(inout) :: vel
    
      call axpy(ndim, -ds / 2.D0, grad(1), 1, vel(1), 1)
      IF ( present(damp) ) THEN
         vel = damp * vel
      END IF
      
    END SUBROUTINE velocity_Verlet_second_step
    
    
    SUBROUTINE quick_min_second_step( ndim, vel, grad, ds )

      IMPLICIT NONE

      integer, intent(in) :: ndim
      real(gp), intent(in) :: ds
      real(gp), dimension(ndim), intent(in) :: grad
      real(gp), dimension(ndim), intent(inout) :: vel

      REAL (gp), DIMENSION(ndim) :: force_versor
      REAL (gp)                  :: vel_component, nrm

      call axpy(ndim, -ds / 2.D0, grad(1), 1, vel(1), 1)
      
      nrm = nrm2(ndim,grad(1),1)
      IF ( nrm >= epsi ) THEN
        force_versor = - grad / nrm
      ELSE
        force_versor = 0.D0
      END IF
 
      vel_component = DOT_PRODUCT( vel , force_versor )
      IF ( vel_component > 0.D0 ) THEN
        vel(:) = vel_component * force_versor
      ELSE
        vel(:) = 0.D0
      END IF

    END SUBROUTINE quick_min_second_step

END MODULE Minimization_routines


module module_images
  use module_defs
  use module_types

  implicit none

  private

  CHARACTER (LEN=*), PARAMETER ::                                              &
  fmt1 = "(3(2X,F12.8),3(2X,I3),3(2X,F12.8))",                                 &
  fmt2 = "(3(2X,F12.8))",                                                      &
  fmt3 = "(2X,F16.8)"

  ! Calculation routines.
  public :: compute_local_tangent, compute_local_gradient

  ! Routines on one image.
  public :: image_init, image_deallocate, image_set_init_vel

  ! Routines for a list of images.
  public :: images_collect_results
  public :: compute_neb_pos
  public :: images_get_activation, images_get_energies, images_get_errors
  public :: images_output_step

  ! Misc.
  public :: write_restart, write_restart_vel, write_dat_files
  public :: termalization

  type, public :: NEB_data
     logical :: optimization, climbing
     integer :: max_iterations 
     real(gp) :: convergence 
     real(gp) :: ds, k_min, k_max
     real(gp) :: damp, temp_req
  end type NEB_data

  type, public :: run_image
     ! Objects to deals with call_bigdft.
     ! Contains the atomic structure definition,
     !  the total energy and the forces (modified to be gradients).
     type(run_objects) :: run
     type(DFT_global_output) :: outs
     character(len = 128) :: log_file
     ! Local convergence quantities
     real(gp) :: error, F
     ! Last running image.
     integer :: id
     ! Private work arrays.
     integer :: algorithm
     real(gp), dimension(:), pointer :: old_grad, delta_pos, vel
  end type run_image

contains

  FUNCTION norm(vect)
    IMPLICIT NONE

    REAL (gp), DIMENSION(:), INTENT(IN)  :: vect
    REAL (gp)                            :: norm

    norm = SQRT( DOT_PRODUCT( vect , vect ) )
  END FUNCTION norm

  FUNCTION cubic_pbc( vect, Lx, Ly, Lz )
    IMPLICIT NONE    

    REAL (gp), DIMENSION(:), INTENT(IN)  :: vect
    real (gp), intent(in) :: Lx, Ly, Lz
    REAL (gp), DIMENSION( SIZE( vect ) ) :: cubic_pbc
    REAL (gp)                            :: invLx, invLy, invLz
    INTEGER                                    :: i, dim

    invLx = 1.D0 / Lx
    invLy = 1.D0 / Ly
    invLz = 1.D0 / Lz

    dim = size(vect)
    DO i = 1, dim
       IF ( MOD(i,3) == 1 .and. Lx /= 0._gp) THEN
          cubic_pbc(i) = vect(i) - ANINT( vect(i) * invLx ) * Lx
       ELSE IF ( MOD(i,3) == 2 .and. Ly /= 0._gp) THEN
          cubic_pbc(i) = vect(i) - ANINT( vect(i) * invLy ) * Ly
       ELSE if ( MOD(i,3) == 0 .and. Lz /= 0._gp) then
          cubic_pbc(i) = vect(i) - ANINT( vect(i) * invLz ) * Lz
       else
          cubic_pbc(i) = vect(i)
       END IF
    END DO
  END FUNCTION cubic_pbc

  subroutine image_init(img, inputs, atoms, rst, algorithm)
    use module_interfaces, only: run_objects_associate
    implicit none
    type(run_image), intent(out) :: img
    type(input_variables), intent(in) :: inputs
    type(atoms_data), intent(in) :: atoms
    type(restart_objects), intent(in) :: rst
    integer, intent(in) :: algorithm

    integer :: ndim

    img%id = -1

    img%error = 999.d99
    img%F     = 999.d99

    nullify(img%old_grad)
    nullify(img%delta_pos)
    nullify(img%vel)

    call run_objects_nullify(img%run)
    call run_objects_associate(img%run, inputs, atoms, rst)
    call init_global_output(img%outs, atoms%astruct%nat)

    write(img%log_file, "(A,A)") trim(inputs%writing_directory), &
         & 'log-'//trim(inputs%run_name)//'.yaml'           

    img%algorithm = algorithm
    ndim = 3 * atoms%astruct%nat
    if (algorithm <= 3) then
       allocate(img%old_grad(ndim))
       call to_zero(ndim, img%old_grad(1))
       allocate(img%delta_pos(ndim))
       call to_zero(ndim, img%delta_pos(1))
    else
       allocate(img%vel(ndim))
       call to_zero(ndim, img%vel(1))
    end if
  end subroutine image_init

  subroutine image_set_init_vel(img, ndim, vel0)
    implicit none
    integer, intent(in) :: ndim
    type(run_image), intent(inout) :: img
    real(gp), dimension(ndim), intent(in) :: vel0

    call vcopy(ndim, vel0(1), 1, img%vel(1), 1)
  end subroutine image_set_init_vel

  subroutine image_deallocate(img, free_subs)
    implicit none
    type(run_image), intent(inout) :: img
    logical, intent(in) :: free_subs

    if (free_subs) then
       call run_objects_free_container(img%run)
       call deallocate_global_output(img%outs)
    end if

    if (associated(img%old_grad)) deallocate(img%old_grad)
    if (associated(img%delta_pos)) deallocate(img%delta_pos)
    if (associated(img%vel)) deallocate(img%vel)
  end subroutine image_deallocate

  function images_get_energies(imgs)
    implicit none
    type(run_image), dimension(:), intent(in) :: imgs

    real(gp), dimension(size(imgs)) :: images_get_energies
    integer :: i
    
    do i = 1, size(imgs)
       images_get_energies(i) = imgs(i)%outs%energy
    end do
  end function images_get_energies

  function images_get_activation(imgs)
    implicit none
    type(run_image), dimension(:), intent(in) :: imgs

    real(gp) :: images_get_activation
    integer :: i

    images_get_activation = -999.d99
    do i = 2, size(imgs) - 1
       images_get_activation = max(images_get_activation, imgs(i)%outs%energy - imgs(1)%outs%energy)
    end do
  end function images_get_activation

  function images_get_errors(imgs)
    implicit none
    type(run_image), dimension(:), intent(in) :: imgs

    real(gp), dimension(size(imgs)) :: images_get_errors
    integer :: i
    
    do i = 1, size(imgs)
       images_get_errors(i) = imgs(i)%error
    end do
  end function images_get_errors

  subroutine compute_local_tangent(tgt, ndim, V, posm1, pos0, posp1, Lx, Ly, Lz)
    implicit none

    integer, intent(in) :: ndim
    real(gp), intent(in) :: Lx, Ly, Lz
    real(gp), dimension(ndim), intent(out) :: tgt
    real(gp), dimension(-1:1), intent(in) :: V
    real(gp), dimension(ndim), intent(in) :: posm1, pos0, posp1

    REAL (gp) :: V_previous, V_actual, V_next
    REAL (gp) :: abs_next, abs_previous
    REAL (gp) :: delta_V_max, delta_V_min

    tgt = 0

    !! tangent to the path (normalized)
    V_previous = V( -1 )
    V_actual   = V( 0 )
    V_next     = V( +1 )

    IF ( ( V_next > V_actual ) .AND. ( V_actual > V_previous ) ) THEN
       tgt = cubic_pbc( posp1 - pos0, Lx, Ly, Lz )
    ELSE IF ( ( V_next < V_actual ) .AND. ( V_actual < V_previous ) ) THEN
       tgt = cubic_pbc( pos0 - posm1, Lx, Ly, Lz )
    ELSE

       abs_next     = ABS( V_next - V_actual ) 
       abs_previous = ABS( V_previous - V_actual ) 

       delta_V_max = MAX( abs_next , abs_previous ) 
       delta_V_min = MIN( abs_next , abs_previous )

       IF ( V_next > V_previous ) THEN
          tgt = &
               cubic_pbc( posp1 - pos0, Lx, Ly, Lz ) * delta_V_max + & 
               cubic_pbc( pos0 - posm1, Lx, Ly, Lz ) * delta_V_min
       ELSE IF ( V_next < V_previous ) THEN
          tgt = &
               cubic_pbc( posp1 - pos0, Lx, Ly, Lz ) * delta_V_min + & 
               cubic_pbc( pos0 - posm1, Lx, Ly, Lz ) * delta_V_max
       ELSE
          tgt = &
               cubic_pbc( posp1 - pos0, Lx, Ly, Lz ) 
       END IF

    END IF

    tgt = tgt / norm( tgt )
  END SUBROUTINE compute_local_tangent

  SUBROUTINE compute_local_gradient(ndim, grad, posm1, pos0, posp1, tgt, PES_forces, Lx, Ly, Lz, &
       & k, full, climbing)
    IMPLICIT NONE

    integer, intent(in) :: ndim
    real(gp), dimension(ndim), intent(out) :: grad
    real(gp), dimension(ndim), intent(in) :: tgt, PES_forces
    real(gp), dimension(ndim), intent(in) :: posm1, pos0, posp1
    real(gp), intent(in) :: Lx, Ly, Lz
    real(gp), dimension(2:3), intent(in) :: k
    logical, intent(in) :: full, climbing

    real(gp), dimension(:), allocatable :: elastic_gradient

    !! total gradient on each replica ( climbing_img image is used if needed )
    !! only the component of the PES gradient ortogonal to the path is taken into
    !! account
    if (climbing) then
       grad = - PES_forces + 2.D0 * dot_product( PES_forces , tgt ) * tgt
    else
       allocate(elastic_gradient(ndim))

       IF ( full ) THEN
          !! elastic gradient ( variable elastic constant is used )
          elastic_gradient = &
               ( k(2) * ( cubic_pbc( pos0 - posm1, Lx, Ly, Lz ) ) - &
                 k(3) * ( cubic_pbc( posp1 - pos0, Lx, Ly, Lz ) ) ) 
       ELSE
          !! elastic gradient only along the path ( variable elastic constant is used )  
          elastic_gradient = &
               ( k(2) * norm( cubic_pbc( pos0 - posm1, Lx, Ly, Lz ) ) - &
                 k(3) * norm( cubic_pbc( posp1 - pos0, Lx, Ly, Lz ) ) ) * tgt
       END IF

       grad = - PES_forces + elastic_gradient + dot_product( PES_forces, tgt ) * tgt

       deallocate(elastic_gradient)
    end if
  END SUBROUTINE compute_local_gradient

  subroutine compute_k(nimages, k, V, k_min, k_max)
    implicit none
    integer, intent(in) :: nimages
    real(gp), dimension(2:nimages), intent(out) :: k
    real(gp), dimension(nimages), intent(in) :: V
    real(gp), intent(in) :: k_min, k_max

    integer :: i
    real(gp) :: Ei, Eref, Emax

    Eref = MIN( V(1) , V(nimages) )
    Emax = maxval(V)

    elastic_const_loop: DO i = 2, nimages 
       IF ( i < nimages ) THEN
          Ei = MAX( MAX( V(i) , V(i-1) ) , MAX( V(i) , V(i+1) ) )
       ELSE
          Ei = MAX( V(i) , V(i-1) )  
       END IF

       IF ( Ei > Eref ) THEN
          k(i) = k_max - (k_max - k_min) * ( ( Emax - Ei ) / ( Emax - Eref ) )
       ELSE
          k(i) = k_min
       END IF
    end do elastic_const_loop
  END SUBROUTINE compute_k

  subroutine compute_neb_pos(imgs, iteration, neb)
    implicit none

    type(NEB_data), intent(in) :: neb
    type(run_image), dimension(:), intent(inout) :: imgs
    integer, intent(in) :: iteration

    integer :: i, n_in, n_fin, climbing_img
    real(gp), dimension(2:size(imgs)) :: k

    IF ( neb%optimization ) THEN
       N_in  = 1
       N_fin = size(imgs)
    ELSE
       N_in  = 2
       N_fin = size(imgs) - 1

       imgs(1)%error = 0.d0
       imgs(1)%F     = 0.d0
       imgs(size(imgs))%error = 0.d0
       imgs(size(imgs))%F     = 0.d0
    END IF

    ! Global line treatment.
    call compute_k(size(imgs), k, images_get_energies(imgs), neb%k_min, neb%k_max)
    climbing_img = 0
    if (neb%climbing) climbing_img = maxloc(images_get_energies(imgs), 1)

    ! Per image treatment.
    do i = N_in, N_fin
       call image_update_pos(imgs(i), iteration, &
            & imgs(max(1,i-1))%run%atoms%astruct%rxyz, &
            & imgs(min(i+1,size(imgs)))%run%atoms%astruct%rxyz, &
            & imgs(max(1,i-1))%outs%energy, imgs(min(i+1,size(imgs)))%outs%energy, &
            & k(i), k(i+1), (i == 1 .or. i == size(imgs)), (i == climbing_img), neb)
    end do

    ! Global line treatment.
    IF ( imgs(1)%algorithm == 6 ) CALL termalization(imgs, neb%temp_req)
  END SUBROUTINE compute_neb_pos
  
  SUBROUTINE termalization(imgs, temp_req)
    IMPLICIT NONE

    real(gp), intent(in) :: temp_req
    type(run_image), dimension(:), intent(inout) :: imgs

    REAL(gp) :: temp, fact
    INTEGER :: i, n

    n = 0
    temp = 0.D0
    DO i = 2, ( size(imgs) - 1 )
       temp = temp + DOT_PRODUCT( imgs(i)%vel , imgs(i)%vel )
       n = n + size(imgs(i)%vel)
    END DO

    temp = temp / DBLE( n )
    fact = SQRT( temp_req / temp )
    do i = 2, size(imgs) - 1
       imgs(i)%vel = imgs(i)%vel * fact
    end do
  END SUBROUTINE termalization

  SUBROUTINE write_restart(restart_file, imgs, fix_atom)
    IMPLICIT NONE

    character(len = *), intent(in) :: restart_file
    type(run_image), dimension(:), intent(in) :: imgs
    real(gp), dimension(:,:), intent(in) :: fix_atom

    INTEGER             :: i, j
    INTEGER, PARAMETER  :: unit = 10


    OPEN( UNIT = unit, FILE = trim(restart_file), STATUS = "UNKNOWN", ACTION = "WRITE" )
    DO i = 1, size(imgs)
       WRITE(unit,*) "Replica: ", i
       WRITE(unit,fmt3) imgs(i)%outs%energy

       DO j = 1, imgs(i)%run%atoms%astruct%nat
          WRITE(unit,fmt1) imgs(i)%run%atoms%astruct%rxyz(1,j),     & 
               imgs(i)%run%atoms%astruct%rxyz(2,j), &
               imgs(i)%run%atoms%astruct%rxyz(3,j), &
               int(fix_atom(1,j)),     &
               int(fix_atom(2,j)), &
               int(fix_atom(3,j)), &
               imgs(i)%outs%fxyz(1,j), &
               imgs(i)%outs%fxyz(2,j), &
               imgs(i)%outs%fxyz(3,j)
       END DO
    END DO
    CLOSE( UNIT = unit )
  END SUBROUTINE write_restart

  SUBROUTINE write_restart_vel(velocity_file, imgs)
    IMPLICIT NONE

    type(run_image), dimension(:), intent(in) :: imgs
    character(len = *), intent(in) :: velocity_file

    INTEGER             :: i, j
    INTEGER, PARAMETER  :: unit = 10

    OPEN( UNIT = unit, FILE = trim(velocity_file), STATUS = "UNKNOWN", ACTION = "WRITE" )
    DO i = 1, size(imgs)
       WRITE(unit,*) "Replica: ", i

       DO j = 1, size(imgs(i)%vel), 3 
          WRITE(unit,fmt2) imgs(i)%vel(j), imgs(i)%vel((j+1)), imgs(i)%vel((j+2))
       END DO
    END DO
    CLOSE( UNIT = unit )
  END SUBROUTINE write_restart_vel

  SUBROUTINE write_dat_files(job_name, imgs, iter)
    IMPLICIT NONE

    character(len = *), intent(in) :: job_name
    integer, intent(in) :: iter
    type(run_image), dimension(:), intent(in) :: imgs

    INTEGER :: i, j, ndim
    REAL (gp)                            :: R, delta_R, x
    REAL (gp), DIMENSION(:), ALLOCATABLE :: d_R
    REAL (gp), DIMENSION(:), ALLOCATABLE :: a, b, c, d
    REAL (gp), DIMENSION(:), ALLOCATABLE :: reaction_coordinate
    REAL (gp)                            :: E, E_0
    INTEGER, PARAMETER :: max_i = 1000 
    INTEGER, PARAMETER :: unit = 10
    character(len = 4) :: fn4

    write(fn4, "(I4.4)") iter

    ndim = imgs(1)%run%atoms%astruct%nat * 3
    ALLOCATE( d_R( ndim ) )

    ALLOCATE( a( size(imgs) - 1 ) )
    ALLOCATE( b( size(imgs) - 1 ) )
    ALLOCATE( c( size(imgs) - 1 ) )
    ALLOCATE( d( size(imgs) - 1 ) )
    ALLOCATE( reaction_coordinate( size(imgs) ) )

    reaction_coordinate(1) = 0.D0
    DO i = 1, ( size(imgs) - 1 )
       call vcopy(ndim, imgs(i + 1)%run%atoms%astruct%rxyz(1,1), 1, d_R(1), 1)
       call axpy(ndim, -1.d0, imgs(i)%run%atoms%astruct%rxyz(1,1), 1, d_R(1), 1)
       d_R = cubic_pbc( d_R, imgs(1)%run%atoms%astruct%cell_dim(1), &
            & imgs(1)%run%atoms%astruct%cell_dim(2), imgs(1)%run%atoms%astruct%cell_dim(3) ) 
       R = norm( d_R )

       reaction_coordinate(i+1) = reaction_coordinate(i) + R
       a(i) = 2.D0 * ( imgs(i)%outs%energy - imgs(i+1)%outs%energy ) / R**(3) - &
            ( imgs(i)%F + imgs(i+1)%F ) / R**(2)
       b(i) = 3.D0 * ( imgs(i+1)%outs%energy - imgs(i)%outs%energy ) / R**(2) + &
            ( 2.D0 * imgs(i)%F + imgs(i+1)%F ) / R
       c(i) = - imgs(i)%F
       d(i) = imgs(i)%outs%energy
    END DO

    OPEN( UNIT = unit, FILE = job_name // ".it" // fn4 // ".dat", STATUS = "UNKNOWN", &
         ACTION = "WRITE" )

    WRITE(unit,'(3(2X,F12.8))') 0.D0, 0.D0, 0.D0
    DO i = 2, size(imgs)
       WRITE(unit,'(3(2X,F12.8))') reaction_coordinate(i) * Bohr_Ang, &
            ( imgs(i)%outs%energy - imgs(1)%outs%energy ) * Ha_eV, imgs(i)%error * ( Ha_eV / Bohr_Ang )
    END DO
    CLOSE( UNIT = unit )  

    OPEN( UNIT = unit, FILE = job_name // ".it" // fn4 // ".int", STATUS = "UNKNOWN", &
         ACTION = "WRITE" )

    i = 1
    delta_R = reaction_coordinate(size(imgs)) / DBLE(max_i)
    DO j = 0, max_i
       R = DBLE(j) * delta_R 
       IF ( ( R > reaction_coordinate(i+1) ) .AND. &
            ( i < ( size(imgs) - 1 ) ) ) i = i + 1

       x = R - reaction_coordinate(i)
       E = a(i)*(x**3) + b(i)*(x**2) + c(i)*x + d(i) 
       IF ( j == 0 ) E_0 = E

       WRITE(unit,'(2(2X,F16.8))') ( R * Bohr_Ang ), ( E - E_0 ) * Ha_eV
    END DO
    CLOSE( UNIT = unit )

    DEALLOCATE( d_R )

    DEALLOCATE( a )
    DEALLOCATE( b )
    DEALLOCATE( c )
    DEALLOCATE( d )
    DEALLOCATE( reaction_coordinate )

  END SUBROUTINE write_dat_files

  subroutine images_output_step(imgs, full, iteration, tol)
    use yaml_output
    implicit none
    type(run_image), dimension(:), intent(in) :: imgs
    logical, intent(in), optional :: full
    integer, intent(in), optional :: iteration
    real(gp), intent(in), optional :: tol

    logical :: full_
    integer :: i
    character(len = size(imgs)) :: scheme
    logical, dimension(size(imgs)) :: update

    full_ = .false.
    if (present(full)) full_ = full

    if (full_) then
       call yaml_open_sequence("Energy and error per image")
       DO i = 1, size(imgs)
          call yaml_sequence(advance='no')
          call yaml_open_map(flow=.true.)
          call yaml_map("Energy (eV)", imgs(i)%outs%energy * Ha_eV,fmt='(F16.8)')
          call yaml_map("error (eV/ang)", imgs(i)%error * ( Ha_eV / Bohr_Ang ),fmt='(F8.5)')
          call yaml_close_map(advance='no')
!!$          call yaml_sequence( &
!!$               dict((/ "Energy (eV)"      .is. yaml_toa(V(i) * Ha_eV,fmt='(F16.8)') ,&
!!$                       "error (eV / ang)" .is. yaml_toa(error(i) * ( Ha_eV / Bohr_Ang ),fmt='(F8.5)') /), &
!!$                       advance = "no")
          call yaml_comment(trim(yaml_toa(i,fmt='(i2.2)')))
       END DO
       call yaml_close_sequence()
    else
       call yaml_open_map(flow=.true.)
       call yaml_map("Ea (eV)", images_get_activation(imgs) * Ha_eV,fmt='(F10.6)')
       if (present(tol)) then
          ! Print the update scheme.
          update = .true.
          where ( images_get_errors(imgs) * Ha_eV / Bohr_Ang <= tol ) update = .false.
          do i = 1, size(imgs)
             if (update(i)) then
                scheme(i:i) = '*'
             else
                scheme(i:i) = '.'
             end if
          end do
          call yaml_map("cv", '"' // trim(scheme) // '"')
       end if
       call yaml_map("max err (eV/ang)", maxval(images_get_errors(imgs)) * ( Ha_eV / Bohr_Ang ),fmt='(F10.6)')
       if (present(iteration)) then
          call yaml_close_map(advance="no")
          call yaml_comment(trim(yaml_toa(iteration, fmt='(i3.3)')))
       else
          call yaml_close_map()
       end if
    end if
  END SUBROUTINE images_output_step

  subroutine images_collect_results(imgs, igroup, nimages, mpi_env)
    implicit none
    integer, intent(in) :: nimages
    type(run_image), dimension(nimages), intent(inout) :: imgs
    integer, dimension(nimages), intent(in) :: igroup
    type(mpi_environment), intent(in) :: mpi_env

    integer :: i, ierr

    ! Reduce the results in case of taskgrouping.
    if (mpi_env%nproc > 1) then
       if (mpi_env%igroup == 0) then
          ! Broadcast between taskgroups.
          do i = 1, nimages
             if (igroup(i) > 0) then
                call mpi_bcast(imgs(i)%outs%energy, 1, MPI_DOUBLE_PRECISION, &
                     & igroup(i) - 1, mpi_env%mpi_comm, ierr)
                call mpi_bcast(imgs(i)%outs%fxyz(1,1), imgs(i)%outs%fdim * 3, MPI_DOUBLE_PRECISION, &
                     & igroup(i) - 1, mpi_env%mpi_comm, ierr)
             end if
          end do
       end if
       ! Broadcast inside taskgroups.
       do i = 1, nimages
          if (igroup(i) > 0) then
             call mpi_bcast(imgs(i)%outs%energy, 1, MPI_DOUBLE_PRECISION, &
                  & 0, bigdft_mpi%mpi_comm, ierr)
             call mpi_bcast(imgs(i)%outs%fxyz(1,1), imgs(i)%outs%fdim * 3, MPI_DOUBLE_PRECISION, &
                  & 0, bigdft_mpi%mpi_comm, ierr)
          end if
       end do
    end if
  END SUBROUTINE images_collect_results

END MODULE module_images


!> Public routines.
subroutine image_update_pos(img, iteration, posm1, posp1, Vm1, Vp1, &
     & km1, kp1, optimization, climbing, neb)
  use Minimization_routines
  use module_images
  implicit none
  type(run_image), intent(inout) :: img
  integer, intent(in) :: iteration
  real(gp), intent(in) :: km1, kp1
  real(gp), intent(in) :: Vm1, Vp1
  real(gp), dimension(3*img%run%atoms%astruct%nat), intent(in) :: posm1, posp1
  logical, intent(in) :: optimization, climbing
  type(NEB_data), intent(in) :: neb

  integer :: ndim
  real(gp) :: Lx, Ly, Lz
  real(gp), dimension(:), allocatable :: tangent, grad

  ! Aliasing
  ndim = 3 * img%run%atoms%astruct%nat
  Lx = img%run%atoms%astruct%cell_dim(1)
  Ly = img%run%atoms%astruct%cell_dim(2)
  Lz = img%run%atoms%astruct%cell_dim(3)

  allocate(grad(ndim))
  call to_zero(ndim, grad(1))

  if (.not. optimization) then
     allocate(tangent(ndim))
     call compute_local_tangent(tangent, ndim, (/ Vm1, img%outs%energy, Vp1 /), &
          & posm1, img%run%atoms%astruct%rxyz, posp1, Lx, Ly, Lz)
  end if

  if (iteration > 0 .and. img%algorithm >= 4 ) then
     if (optimization) then
        call to_zero(ndim, grad(1))
        call axpy(ndim, -1.d0, img%outs%fxyz(1,1), 1, grad(1), 1)
     else
        call compute_local_gradient(ndim, grad, posm1, img%run%atoms%astruct%rxyz, posp1, tangent, &
             & img%outs%fxyz, Lx,Ly,Lz, (/ km1, kp1 /), (img%algorithm == 6), climbing)
     end if
     IF ( img%algorithm == 4 ) THEN
        CALL quick_min_second_step( ndim, img%vel, grad, neb%ds)
     ELSE IF ( img%algorithm == 5 ) THEN
        CALL velocity_Verlet_second_step( ndim, img%vel, grad, neb%ds, neb%damp )
     ELSE IF ( img%algorithm == 6 ) THEN
        CALL velocity_Verlet_second_step( ndim, img%vel, grad, neb%ds )
     END IF

!!$       IF ( algorithm == 6 ) CALL termalization(neb%ndim, neb%nimages, neb%vel, neb%temp_req)
  end if

  if (iteration > 0 .or. neb%max_iterations == 1) then
     img%error = real(0, gp)
     img%F     = real(0, gp)
     if (optimization) then
        img%error = nrm2(ndim, img%outs%fxyz(1,1), 1)
     else
        img%F = - dot(ndim, img%outs%fxyz(1,1), 1, tangent(1), 1)
        call vcopy(ndim, img%outs%fxyz(1,1), 1, grad(1), 1)
        call axpy(ndim, img%F, tangent(1), 1, grad(1), 1)
        img%error = nrm2(ndim, grad(1), 1)
     end if
  else
     img%error = real(999, gp)
     img%F     = real(999, gp)
  end if

  ! Calculate new positions for next step.
  if (optimization) then
     call to_zero(ndim, grad(1))
     call axpy(ndim, -1.d0, img%outs%fxyz(1,1), 1, grad(1), 1)
  else
     call compute_local_gradient(ndim, grad, posm1, img%run%atoms%astruct%rxyz, posp1, tangent, &
          & img%outs%fxyz, Lx,Ly,Lz, (/ km1, kp1 /), (img%algorithm == 6), climbing)
  end if
  IF ( img%algorithm == 1 ) THEN
     CALL steepest_descent( ndim, img%run%atoms%astruct%rxyz, grad, neb%ds)
  ELSE IF ( img%algorithm == 2 ) THEN
     CALL fletcher_reeves( ndim, img%run%atoms%astruct%rxyz, grad, img%old_grad, img%delta_pos, neb%ds)
  ELSE IF ( img%algorithm == 3 ) THEN
     CALL polak_ribiere( ndim, img%run%atoms%astruct%rxyz, grad, img%old_grad, img%delta_pos, neb%ds)
  ELSE IF ( img%algorithm >= 4 ) THEN
     CALL velocity_Verlet_first_step( ndim, img%run%atoms%astruct%rxyz, img%vel, grad, neb%ds)
  END IF

  IF ( img%algorithm <= 3 ) call vcopy(ndim, grad(1), 1, img%old_grad(1), 1)

  if (.not. optimization) then
     deallocate(tangent)
  end if

  deallocate(grad)
END SUBROUTINE image_update_pos


subroutine image_update_pos_from_file(img, iteration, filem1, filep1, km1, kp1, climbing, neb)
  use Minimization_routines
  use module_types
  use module_images
  use module_interfaces, only: read_atomic_file
  implicit none
  type(run_image), intent(inout) :: img
  integer, intent(in) :: iteration
  character(len = *), intent(in) :: filem1, filep1
  real(gp), intent(in) :: km1, kp1
  logical, intent(in) :: climbing
  type(NEB_data), intent(in) :: neb

  character(len = *), parameter :: subname = "image_update_pos_from_file"
  real(gp), dimension(:,:), pointer :: rxyzm1, rxyzp1
  type(atomic_structure) :: astruct
  real(gp) :: Vm1, Vp1
  integer :: stat

  img%error = UNINITIALIZED(real(1, gp))
  nullify(rxyzm1)
  nullify(rxyzp1)
  call astruct_nullify(astruct)

  if (trim(filem1) /= "") then
     call read_atomic_file(trim(filem1), bigdft_mpi%iproc, astruct, &
          & status = stat, energy = Vm1)
     if (stat /= 0 .or. astruct%nat /= img%run%atoms%astruct%nat) then
        call free_me()
        return
     end if
     rxyzm1 => astruct%rxyz
     nullify(astruct%rxyz)
     call deallocate_atomic_structure(astruct, subname)
     call astruct_nullify(astruct)
  end if

  if (trim(filep1) /= "") then
     call read_atomic_file(trim(filep1), bigdft_mpi%iproc, astruct, &
          & status = stat, energy = Vp1)
     if (stat /= 0 .or. astruct%nat /= img%run%atoms%astruct%nat) then
        call free_me()
        return
     end if
     rxyzp1 => astruct%rxyz
     nullify(astruct%rxyz)
     call deallocate_atomic_structure(astruct, subname)
     call astruct_nullify(astruct)
  end if
  
  call image_update_pos(img, iteration, rxyzm1, rxyzp1, Vm1, Vp1, km1, kp1, &
       & .not. associated(rxyzm1) .or. .not. associated(rxyzp1), climbing, neb)
  call free_me()

contains

  subroutine free_me()
    implicit none
    integer :: i_all, i_stat
    character(len = *), parameter :: subname = "image_update_pos_from_file"

    if (associated(rxyzp1)) then
      i_all=-product(shape(rxyzp1))*kind(rxyzp1)
      deallocate(rxyzp1,stat=i_stat)
      call memocc(i_stat,i_all,'rxyzp1',subname)
    end if
    if (associated(rxyzm1)) then
      i_all=-product(shape(rxyzm1))*kind(rxyzm1)
      deallocate(rxyzm1,stat=i_stat)
      call memocc(i_stat,i_all,'rxyzm1',subname)
    end if
    call deallocate_atomic_structure(astruct, subname)
  end subroutine free_me
END SUBROUTINE image_update_pos_from_file

subroutine image_calculate(img, iteration, id)
  use yaml_output
  use module_types
  use module_images
  use module_interfaces, only: write_atomic_file
  implicit none
  type(run_image), intent(inout) :: img
  integer :: iteration
  integer, intent(in) :: id

  integer :: ierr, infocode
  character(len = 4) :: fn4

  !Why (TD) ??
  img%run%inputs%inputpsiid = 0
  if (iteration > 0 .and. abs(img%id - id) < 2) img%run%inputs%inputpsiid = 1

  img%id = id
  if (trim(img%log_file) /= "" .and. bigdft_mpi%iproc == 0) then
     call yaml_set_stream(unit = 9169 + id, filename = trim(img%log_file), istat = ierr)
     call yaml_comment("NEB iteration #" // trim(yaml_toa(iteration, fmt = "(I3.3)")), hfill="-")
  end if
  call call_bigdft(img%run, img%outs, bigdft_mpi%nproc, bigdft_mpi%iproc, infocode)
  if (trim(img%log_file) /= "" .and. bigdft_mpi%iproc == 0) call yaml_close_all_streams()

  ! Output the corresponding file.
  if (bigdft_mpi%iproc == 0) then
     write(fn4, "(I4.4)") iteration
     call write_atomic_file(trim(img%run%inputs%dir_output)//'posout_'//fn4, &
          & img%outs%energy, img%run%atoms%astruct%rxyz, img%run%atoms, "", forces = img%outs%fxyz)
  end if
end subroutine image_calculate

subroutine images_distribute_tasks(igroup, update, nimages, ngroup)
  implicit none
  integer, intent(in) :: nimages, ngroup
  logical, dimension(nimages), intent(in) :: update
  integer, dimension(nimages), intent(out) :: igroup

  integer :: alpha, beta
  integer :: i, l, m, n

  n = 0
  do i = 1, nimages
     if (update(i)) n = n + 1
  end do

  igroup(:) = -1
  alpha = 0
  beta = nimages + 1
  ! Taskgroups have l or (l+1) images to compute.
  l = n / ngroup
  ! There are m taskgroups with (l+1) images.
  m = n - l * ngroup
  do i = 1, m
     if (modulo(i, 2) == 1) then
        call span_group(alpha, l + 1, +1, i)
     else
        call span_group(beta,  l + 1, -1, i)
     end if
  end do
  ! There are n - m taskgroups with l images.
  do i = m + 1, ngroup * min(1, l)
     call span_group(alpha, l, +1, i)
  end do
contains
  subroutine span_group(it, n, dir, ig)
    integer, intent(inout) :: it
    integer, intent(in) :: n, dir, ig
    integer :: j

    j = 0
    do
       it = it + dir
       if (update(it)) then
          igroup(it) = ig
          j = j + 1
          if (j == n) exit
       end if
    end do
  end subroutine span_group
END SUBROUTINE images_distribute_tasks

! Routines for bindings.
subroutine image_new(img, run, outs, atoms, inputs, rst, algorithm)
  use module_types
  use module_images
  implicit none

  type(run_image), pointer :: img
  type(run_objects), pointer :: run
  type(DFT_global_output), pointer :: outs
  type(input_variables), intent(in) :: inputs
  type(atoms_data), intent(in) :: atoms
  type(restart_objects), intent(in) :: rst
  integer, intent(in) :: algorithm

  allocate(img)
  call image_init(img, inputs, atoms, rst, algorithm)
  run => img%run
  outs => img%outs
END SUBROUTINE image_new

subroutine image_free(img, run, outs)
  use module_types
  use module_images
  implicit none

  type(run_image), pointer :: img
  type(run_objects), pointer :: run
  type(DFT_global_output), pointer :: outs

  allocate(run)
  run = img%run

  allocate(outs)
  outs = img%outs

  call image_deallocate(img, .false.)
  deallocate(img)
END SUBROUTINE image_free

subroutine image_get_attributes(img, error, F, id)
  use module_images
  use module_types
  implicit none

  type(run_image), intent(in) :: img
  real(gp), intent(out) :: error, F
  integer, intent(out) :: id

  id = img%id
  error = img%error
  F = img%F
END SUBROUTINE image_get_attributes

subroutine neb_new(neb)
  use module_images
  implicit none
  type(NEB_data), pointer :: neb

  allocate(neb)
  neb%damp = 1.d0
  neb%ds = 0.8d0
  neb%k_min = 0.05d0
  neb%k_max = 0.1d0
  neb%max_iterations = 2
  neb%convergence = 0.5d0
end subroutine neb_new

subroutine neb_free(neb)
  use module_images
  implicit none
  type(NEB_data), pointer :: neb

  deallocate(neb)
end subroutine neb_free
