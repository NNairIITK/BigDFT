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

  implicit none

  private

  CHARACTER (LEN=*), PARAMETER ::                                              &
  fmt1 = "(3(2X,F12.8),3(2X,I1),3(2X,F12.8))",                                 &
  fmt2 = "(3(2X,F12.8))",                                                      &
  fmt3 = "(2X,F16.8)"

  public :: init_images, deallocate_images, set_init_vel
  public :: compute_tangent, compute_gradient, compute_error, compute_neb_pos
  public :: write_restart, write_restart_vel, write_dat_files
  public :: termalization

  type, public :: NEB_data
     logical :: optimization, climbing
     integer :: max_iterations 
     real(gp) :: convergence 
     real(gp) :: ds, k_min, k_max
     real(gp) :: damp, temp_req
     integer :: algorithm
     integer :: ndim, nimages
     real(gp), dimension(:,:), pointer :: old_grad, delta_pos, vel
  end type NEB_data

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
    REAL (gp)                            :: vect_i, invLx, invLy, invLz
    INTEGER                                    :: i, dim

    invLx = 1.D0 / Lx
    invLy = 1.D0 / Ly
    invLz = 1.D0 / Lz

    dim = size(vect)
    DO i = 1, dim
       vect_i = vect(i)
       IF ( MOD(i,3) == 1 ) THEN
          cubic_pbc(i) = vect_i - ANINT( vect_i * invLx ) * Lx
       ELSE IF ( MOD(i,3) == 2 ) THEN
          cubic_pbc(i) = vect_i - ANINT( vect_i * invLy ) * Ly
       ELSE
          cubic_pbc(i) = vect_i - ANINT( vect_i * invLz ) * Lz
       END IF
    END DO
  END FUNCTION cubic_pbc

  subroutine init_images(neb, ndim, nimages, algorithm)
    implicit none
    type(NEB_data), intent(out) :: neb
    integer, intent(in) :: ndim, nimages, algorithm

    neb%algorithm = algorithm
    neb%ndim = ndim
    neb%nimages = nimages
    if (algorithm <= 3) then
       allocate(neb%old_grad(ndim, nimages))
       call to_zero(ndim * nimages, neb%old_grad(1,1))
       allocate(neb%delta_pos(ndim, nimages))
       call to_zero(ndim * nimages, neb%delta_pos(1,1))
    else
       allocate(neb%vel(ndim, nimages))
       call to_zero(ndim * nimages, neb%vel(1,1))
    end if
  end subroutine init_images

  subroutine set_init_vel(neb, vel0)
    implicit none
    type(NEB_data), intent(inout) :: neb
    real(gp), dimension(neb%ndim, neb%nimages), intent(in) :: vel0

    call vcopy(neb%ndim * neb%nimages, vel0(1,1), 1, neb%vel(1,1), 1)
  end subroutine set_init_vel

  subroutine deallocate_images(neb)
    implicit none
    type(NEB_data), intent(inout) :: neb

    if (associated(neb%old_grad)) deallocate(neb%old_grad)
    if (associated(neb%delta_pos)) deallocate(neb%delta_pos)
    if (associated(neb%vel)) deallocate(neb%vel)
  end subroutine deallocate_images

  subroutine compute_tangent(tgt, ndim, nimages, V, pos, Lx, Ly, Lz)
    implicit none

    integer, intent(in) :: nimages, ndim
    real(gp), intent(in) :: Lx, Ly, Lz
    real(gp), dimension(ndim,nimages), intent(out) :: tgt
    real(gp), dimension(nimages), intent(in) :: V
    real(gp), dimension(ndim, nimages), intent(in) :: pos

    INTEGER :: i   
    REAL (gp) :: V_previous, V_actual, V_next
    REAL (gp) :: abs_next, abs_previous
    REAL (gp) :: delta_V_max, delta_V_min

    tgt = 0

    DO i = 2, ( nimages - 1 )

       !! tangent to the path (normalized)
       V_previous = V( i - 1 )
       V_actual   = V( i )
       V_next     = V( i + 1 )

       IF ( ( V_next > V_actual ) .AND. ( V_actual > V_previous ) ) THEN

          tgt(:,i) = cubic_pbc( pos(:,( i + 1 )) - pos(:,i), Lx, Ly, Lz )

       ELSE IF ( ( V_next < V_actual ) .AND. ( V_actual < V_previous ) ) THEN

          tgt(:,i) = cubic_pbc( pos(:,i) - pos(:,( i - 1 )), Lx, Ly, Lz )

       ELSE

          abs_next     = ABS( V_next - V_actual ) 
          abs_previous = ABS( V_previous - V_actual ) 

          delta_V_max = MAX( abs_next , abs_previous ) 
          delta_V_min = MIN( abs_next , abs_previous )

          IF ( V_next > V_previous ) THEN

             tgt(:,i) = &
                  cubic_pbc( pos(:,( i + 1 )) - pos(:,i), Lx, Ly, Lz ) * delta_V_max + & 
                  cubic_pbc( pos(:,i) - pos(:,( i - 1 )), Lx, Ly, Lz ) * delta_V_min

          ELSE IF ( V_next < V_previous ) THEN

             tgt(:,i) = &
                  cubic_pbc( pos(:,( i + 1 )) - pos(:,i), Lx, Ly, Lz ) * delta_V_min + & 
                  cubic_pbc( pos(:,i) - pos(:,( i - 1 )), Lx, Ly, Lz ) * delta_V_max

          ELSE

             tgt(:,i) = &
                  cubic_pbc( pos(:,( i + 1 )) - pos(:,( i - 1 )), Lx, Ly, Lz ) 

          END IF

       END IF

       tgt(:,i) = tgt(:,i)/ norm( tgt(:,i) )

    END DO

  END SUBROUTINE compute_tangent

  SUBROUTINE compute_gradient(ndim, nimages, grad, V, pos, tgt, PES_grad, Lx, Ly, Lz, &
       & k_min, k_max, optimization, algorithm, climbing)
    IMPLICIT NONE

    integer, intent(in) :: nimages, ndim, algorithm
    real(gp), dimension(ndim, nimages), intent(out) :: grad
    real(gp), dimension(nimages), intent(in) :: V
    real(gp), dimension(ndim, nimages), intent(in) :: pos, tgt, PES_grad
    real(gp), intent(in) :: k_max, k_min, Lx, Ly, Lz
    logical, intent(in) :: optimization, climbing

    INTEGER         :: i, Emax_index
    REAL (gp) :: Ei, Eref, Emax
    real(gp), dimension(nimages) :: k
    real(gp), dimension(:), allocatable :: elastic_gradient

    ALLOCATE( elastic_gradient( ndim ) )

    Eref = MIN( V(1) , V(nimages) )
    Emax = maxval(V)
    Emax_index = maxloc(V, 1)
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

    END DO elastic_const_loop

    call to_zero(ndim * nimages, grad(1,1))

    IF ( optimization ) THEN
       grad(:,1)       = PES_grad(:,1)
       grad(:,nimages) = PES_grad(:,nimages)
    END IF

    gradient_loop: DO i = 2, ( nimages - 1 )

       IF ( algorithm == 6 ) THEN

          !! elastic gradient ( variable elastic constant is used )
          elastic_gradient = &
               ( k(i) * ( cubic_pbc( pos(:,i) - pos(:,(i-1)), Lx, Ly, Lz ) ) - &
               k(i+1) * ( cubic_pbc( pos(:,(i+1)) - pos(:,i), Lx, Ly, Lz ) ) ) 

       ELSE

          !! elastic gradient only along the path ( variable elastic constant is used )  
          elastic_gradient = &
               ( k(i) * norm( cubic_pbc( pos(:,i) - pos(:,(i-1)), Lx, Ly, Lz ) ) - &
               k(i+1) * norm( cubic_pbc( pos(:,(i+1)) - pos(:,i), Lx, Ly, Lz ) ) ) * tgt(:,i)

       END IF

       !! total gradient on each replica ( climbing image is used if needed )
       !! only the component of the PES gradient ortogonal to the path is taken into
       !! account

       IF ( ( i == Emax_index ) .AND. ( climbing ) ) THEN

          grad(:,i) = PES_grad(:,i) - 2.D0 * &
               DOT_PRODUCT( PES_grad(:,i) , tgt(:,i) ) * tgt(:,i)

       ELSE

          grad(:,i) = PES_grad(:,i) + elastic_gradient - &
               DOT_PRODUCT( PES_grad(:,i) , tgt(:,i) ) * tgt(:,i)

       END IF

    END DO gradient_loop

    deallocate(elastic_gradient)

  END SUBROUTINE compute_gradient

  SUBROUTINE compute_error( error, ndim, nimages, grad, tgt, PES_grad )
    IMPLICIT NONE

    integer, intent(in) :: nimages, ndim
    real(gp), dimension(ndim, nimages), intent(in) :: tgt, PES_grad
    real(gp), dimension(ndim, nimages), intent(in) :: grad
    real(gp), dimension(nimages), intent(out) :: error

    INTEGER :: i

    error(1) = norm(grad(:,1))
    DO i = 2, ( nimages - 1 ) 

       !! the error is given by the norm of the component of the 
       !! total energy gradient ortogonal to the path
       error(i) = norm( PES_grad(:,i) - DOT_PRODUCT( PES_grad(:,i) , tgt(:,i) ) * tgt(:,i) )

    END DO
    error(nimages) = norm(grad(:,nimages ))
  END SUBROUTINE compute_error

  subroutine compute_neb_pos(loop, iteration, error, F, pos, V, PES_grad, Lx, Ly, Lz, neb)
    use Minimization_routines

    implicit none

    type(NEB_data), intent(inout) :: neb
    integer, intent(inout) :: iteration
    real(gp), intent(in) :: Lx, Ly, Lz
    real(gp), dimension(neb%nimages), intent(in) :: V
    real(gp), dimension(neb%ndim,neb%nimages), intent(in) :: PES_grad
    real(gp), dimension(neb%ndim,neb%nimages), intent(inout) :: pos
    real(gp), dimension(neb%nimages), intent(out) :: error, F
    logical, intent(out) :: loop

    integer :: i, n_in, n_fin, istat
    real(gp) :: err
    real(gp), dimension(:,:), allocatable :: tangent, grad
    logical :: file_exists
    CHARACTER (LEN=4), PARAMETER :: exit_file = "EXIT"  

    loop = .true.
    IF ( neb%optimization ) THEN
       N_in  = 1
       N_fin = neb%nimages
    ELSE
       N_in  = 2
       N_fin = neb%nimages - 1
    END IF
    allocate(tangent(neb%ndim, neb%nimages), grad(neb%ndim, neb%nimages))

!!$    do i = 2, neb%nimages - 1
!!$       CALL compute_tangent(tangent, neb%ndim, 3, V(i-1:i+1), pos(:,i-1:i+1), Lx, Ly, Lz)
!!$    end do
    CALL compute_tangent(tangent, neb%ndim, neb%nimages, V, pos, Lx, Ly, Lz)

    if (iteration > 0 .and. neb%algorithm >= 4 ) THEN

       CALL compute_gradient(neb%ndim, neb%nimages, grad, &
            & V, pos, tangent, PES_grad, Lx, Ly, Lz, &
            & neb%k_min, neb%k_max, neb%optimization, neb%algorithm, neb%climbing)
       second_minimization_loop: DO i = N_in, N_fin 
          IF ( neb%algorithm == 4 ) THEN
             CALL quick_min_second_step( neb%ndim, neb%vel(:,i), grad(:,i), neb%ds)
          ELSE IF ( neb%algorithm == 5 ) THEN
             CALL velocity_Verlet_second_step( neb%ndim, neb%vel(:,i), grad(:,i), neb%ds, neb%damp )
          ELSE IF ( neb%algorithm == 6 ) THEN
             CALL velocity_Verlet_second_step( neb%ndim, neb%vel(:,i), grad(:,i), neb%ds )
          END IF
       END DO second_minimization_loop
       IF ( neb%algorithm == 6 ) CALL termalization(neb%ndim, neb%nimages, neb%vel, neb%temp_req)

    end if

    if (iteration > 0 .or. neb%max_iterations == 1) then
       CALL compute_error( error, neb%ndim, neb%nimages, grad, tangent, PES_grad )
       err = maxval(error)

       F = 0.D0
       DO i = 2, ( neb%nimages - 1 )
          F(i) = DOT_PRODUCT( - PES_grad(:,i) , tangent(:,i) )
       END DO

       IF ( ( err * Ha_eV / Bohr_Ang ) <= neb%convergence .or. neb%max_iterations == 1)  THEN
          loop = .false.
          return
       END IF
    end if

    iteration = iteration + 1

    IF ( iteration > neb%max_iterations ) THEN
       loop = .false.
       return
    END IF

    inquire(FILE = exit_file, EXIST = file_exists)
    IF ( file_exists ) THEN
       call delete(trim(exit_file),len(trim(exit_file)),istat)

       WRITE(*,*) " WARNING :  soft exit required"
       WRITE(*,*) " STOPPING ...                 "

       loop = .false.
       return
    END IF

    IF ( neb%ALGORITHM <= 3 ) neb%old_grad = grad
    CALL compute_gradient(neb%ndim, neb%nimages, grad, V, pos, &
         & tangent, PES_grad, Lx, Ly, Lz, &
         & neb%k_min, neb%k_max, neb%optimization, neb%algorithm, neb%climbing)

    ! Calculate new positions for next step.
    first_minimization_loop: DO i = N_in, N_fin
       IF ( neb%algorithm == 1 ) THEN
          CALL steepest_descent( neb%ndim, pos(:,i), grad(:,i), neb%ds)
       ELSE IF ( neb%algorithm == 2 ) THEN
          CALL fletcher_reeves( neb%ndim, pos(:,i), grad(:,i), neb%old_grad(:,i), neb%delta_pos(:,i), neb%ds)
       ELSE IF ( neb%algorithm == 3 ) THEN
          CALL polak_ribiere( neb%ndim, pos(:,i), grad(:,i), neb%old_grad(:,i), neb%delta_pos(:,i), neb%ds)
       ELSE IF ( neb%algorithm >= 4 ) THEN
          CALL velocity_Verlet_first_step( neb%ndim, pos(:,i), neb%vel(:,i), grad(:,i), neb%ds)
       END IF
    END DO first_minimization_loop

    deallocate(tangent, grad)
  END SUBROUTINE compute_neb_pos
  
  SUBROUTINE termalization(ndim, nimages, vel, temp_req)
    IMPLICIT NONE

    integer, intent(in) :: ndim, nimages
    real(gp), intent(in) :: temp_req
    real(gp), dimension(ndim, nimages), intent(inout) :: vel

    REAL (gp)  :: temp
    INTEGER          :: i

    temp = 0.D0  
    DO i = 2, ( nimages - 1 )
       temp = temp + DOT_PRODUCT( vel(:,i) , vel(:,i) )
    END DO

    temp = temp / DBLE( ( nimages - 2 ) * ndim )
    vel = vel * SQRT( temp_req / temp )

  END SUBROUTINE termalization

  SUBROUTINE write_restart(restart_file, ndim, nimages, V, pos, fix_atom, PES_grad)
    IMPLICIT NONE

    integer, intent(in) :: ndim, nimages
    character(len = *), intent(in) :: restart_file
    real(gp), dimension(nimages), intent(in) :: V
    real(gp), dimension(ndim, nimages), intent(in) :: pos, PES_grad
    real(gp), dimension(ndim), intent(in) :: fix_atom

    INTEGER             :: i, j
    INTEGER, PARAMETER  :: unit = 10


    OPEN( UNIT = unit, FILE = trim(restart_file), STATUS = "UNKNOWN", ACTION = "WRITE" )
    DO i = 1, nimages
       WRITE(unit,*) "Replica: ", i
       WRITE(unit,fmt3) V(i)

       DO j = 1, ndim, 3 
          WRITE(unit,fmt1) pos(j,i),     & 
               pos((j+1),i), &
               pos((j+2),i), &
               int(fix_atom(j)),     &
               int(fix_atom((j+1))), &
               int(fix_atom((j+2))), &
               -PES_grad(j,i),     &
               -PES_grad((j+1),i), &
               -PES_grad((j+2),i)
       END DO
    END DO
    CLOSE( UNIT = unit )
  END SUBROUTINE write_restart

  SUBROUTINE write_restart_vel(neb, velocity_file)
    IMPLICIT NONE

    type(NEB_data), intent(in) :: neb
    character(len = *), intent(in) :: velocity_file

    INTEGER             :: i, j
    INTEGER, PARAMETER  :: unit = 10

    OPEN( UNIT = unit, FILE = trim(velocity_file), STATUS = "UNKNOWN", ACTION = "WRITE" )
    DO i = 1, neb%nimages
       WRITE(unit,*) "Replica: ", i

       DO j = 1, neb%ndim, 3 
          WRITE(unit,fmt2) neb%vel(j,i), neb%vel((j+1),i), neb%vel((j+2),i)
       END DO
    END DO
    CLOSE( UNIT = unit )
  END SUBROUTINE write_restart_vel

  SUBROUTINE write_dat_files(data_file, interpolation_file, ndim, nimages, pos, V, F, error, Lx, Ly, Lz)
    IMPLICIT NONE

    character(len = *), intent(in) :: data_file, interpolation_file
    integer, intent(in) :: ndim, nimages
    real(gp), intent(in) :: Lx, Ly, Lz
    real(gp), dimension(ndim, nimages), intent(in) :: pos
    real(gp), dimension(nimages), intent(in) :: F, V, error

    INTEGER                                    :: i, j
    REAL (gp)                            :: R, delta_R, x
    REAL (gp), DIMENSION(:), ALLOCATABLE :: d_R
    REAL (gp), DIMENSION(:), ALLOCATABLE :: a, b, c, d
    REAL (gp), DIMENSION(:), ALLOCATABLE :: reaction_coordinate
    REAL (gp)                            :: E, E_0
    INTEGER, PARAMETER                         :: max_i = 1000 
    INTEGER, PARAMETER                         :: unit = 10


    ALLOCATE( d_R( ndim ) )

    ALLOCATE( a( nimages - 1 ) )
    ALLOCATE( b( nimages - 1 ) )
    ALLOCATE( c( nimages - 1 ) )
    ALLOCATE( d( nimages - 1 ) )
    ALLOCATE( reaction_coordinate( nimages ) )

    reaction_coordinate(1) = 0.D0
    DO i = 1, ( nimages - 1 )
       d_R = cubic_pbc( pos(:,( i + 1 )) - pos(:,i), Lx, Ly, Lz ) 
       R = norm( d_R )

       reaction_coordinate(i+1) = reaction_coordinate(i) + R
       a(i) = 2.D0 * ( V(i) - V(i+1) ) / R**(3) - &
            ( F(i) + F(i+1) ) / R**(2)
       b(i) = 3.D0 * ( V(i+1) - V(i) ) / R**(2) + &
            ( 2.D0 * F(i) + F(i+1) ) / R
       c(i) = - F(i)
       d(i) = V(i)
    END DO

    OPEN( UNIT = unit, FILE = trim(data_file), STATUS = "UNKNOWN", &
         ACTION = "WRITE" )

    WRITE(unit,'(3(2X,F12.8))') 0.D0, 0.D0, 0.D0
    DO i = 2, nimages
       WRITE(unit,'(3(2X,F12.8))') reaction_coordinate(i) * Bohr_Ang, &
            ( V(i) - V(1) ) * Ha_eV, error(i) * ( Ha_eV / Bohr_Ang )
    END DO
    CLOSE( UNIT = unit )  

    OPEN( UNIT = unit, FILE = trim(interpolation_file), STATUS = "UNKNOWN", &
         ACTION = "WRITE" )

    i = 1
    delta_R = reaction_coordinate(nimages) / DBLE(max_i)
    DO j = 0, max_i
       R = DBLE(j) * delta_R 
       IF ( ( R > reaction_coordinate(i+1) ) .AND. &
            ( i < ( nimages - 1 ) ) ) i = i + 1

       x = R - reaction_coordinate(i)
       E = a(i)*(x**3) + b(i)*(x**2) + c(i)*x + d(i) 
       IF ( j == 0 ) E_0 = E

       WRITE(unit,'(2(2X,F16.8))') ( R * Bohr_Ang ), &
            ( E - E_0 ) * Ha_eV
    END DO
    CLOSE( UNIT = unit )

    DEALLOCATE( d_R )

    DEALLOCATE( a )
    DEALLOCATE( b )
    DEALLOCATE( c )
    DEALLOCATE( d )
    DEALLOCATE( reaction_coordinate )

  END SUBROUTINE write_dat_files

END MODULE module_images

! Routines for bindings.
subroutine neb_new(neb, nat, nimages, algorithm)
  use module_images
  implicit none

  type(NEB_data), pointer :: neb
  integer, intent(in) :: nat, nimages, algorithm

  allocate(neb)
  call init_images(neb, 3 * nat, nimages, algorithm)
END SUBROUTINE neb_new

subroutine neb_free(neb)
  use module_images
  implicit none

  type(NEB_data), pointer :: neb

  call deallocate_images(neb)
  deallocate(neb)
end subroutine neb_free
