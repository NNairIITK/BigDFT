!> @file 
!! NEB routines
!! The IO with the external program is performed using atomic units. 
!! The restart file is in atomic units too.
!! Both output files ( int and dat files ) are in angstrom and eV.
!!
!! PES energies and gradients are obtained calling the NEB_driver.sh and 
!! reading the gen_output_file.
!!
!! References :
!! - G. Henkelman, B.P. Uberuaga, H. Jonsson; J.Chem.Phys., 113, 9901, (2000)
!! - G. Henkelman and H. Jonsson; J.Chem.Phys., 113, 9978, (2000)
!! - H. Jonsson, G. Mills, K.W. Jacobsen, "Nudged elastic band method for finding
!!   minimum energy paths of transitions", in Classical and Quantum Dynamics in 
!!   Condensed Phase Simulations, edited by B.J.Berne, G.Ciccotti, D.F.Coker 
!!   (World Scientific, Singapore, 1998), pag. 385 .
!!
!! @author
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!          COPYRIGHT (C) 2003 Carlo Sbraccia.                                !!
!!                   modifications: 2009 Damien Caliste (DC)                  !!
!!          This file is distributed under the terms                          !!
!!          of the GNU General Public License.                                !!
!!          See http://www.gnu.org/copyleft/gpl.txt .                         !!
!!                                                                            !!
!!    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,         !!
!!    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF      !!
!!    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND                   !!
!!    NONINFRINGEMENT.  IN NO EVENT SHALL CARLO SBRACCIA BE LIABLE FOR ANY    !!
!!    CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,    !!
!!    TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE       !!
!!    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! @todo
!!  Group NEB modules


!> Module for NEB calculations (define kinds of real)
MODULE Numeric

  IMPLICIT NONE
  INTEGER, PARAMETER :: sgl = KIND(1.0)
  INTEGER, PARAMETER :: dbl = KIND(1.D0)

  REAL (KIND=dbl), PARAMETER   :: a_zero = 0.529177D0 ! from Bohr to angstrom
  REAL (KIND=dbl), PARAMETER   :: E_zero = 27.212D0   ! from a.u. to eV 

END MODULE Numeric


!> Module for NEB calculations (define formats)
MODULE Formats

  IMPLICIT NONE
  CHARACTER (LEN=*), PARAMETER ::                                              &
  fmt1 = "(3(2X,F12.8),3(2X,I1),3(2X,F12.8))",                                 &
  fmt2 = "(3(2X,F12.8))",                                                      &
  fmt3 = "(2X,F16.8)",                                                         &
  fmt4 = "(' iteration: ',I3,5X,'E activation =',F10.6,5X,'error =',F10.6)",   &
  fmt5 = "(' image: ',I2,'   Energy=  ',F16.8,'   Error=',F8.5)"
  
END MODULE Formats


!> Module for NEB calculations (variables)
MODULE NEB_variables

  USE numeric

  IMPLICIT NONE
 
  CHARACTER (LEN=80)                           :: first_config, last_config
  CHARACTER (LEN=80)                           :: scratch_dir
  CHARACTER (LEN=80)                           :: data_file, &
                                                  interpolation_file, &
                                                  barrier_file
  CHARACTER (LEN=80)                           :: restart_file
  CHARACTER (LEN=80)                           :: job_name
  LOGICAL                                      :: restart, climbing
  LOGICAL                                      :: optimization
  INTEGER                                      :: N, dim
  REAL (KIND=dbl)                              :: Lx, Ly, Lz
  REAL (KIND=dbl)                              :: k_max, k_min
  REAL (KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: pos
  REAL (KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: PES_gradient
  INTEGER, DIMENSION(:), ALLOCATABLE           :: fix_atom
  REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE   :: V, F, error
  REAL (KIND=dbl)                              :: Emax
  INTEGER(KIND=sgl)                            :: Emax_index

  INTEGER                                      :: algorithm
  INTEGER                                      :: num_of_images
  INTEGER                                      :: max_iterations 
  REAL (KIND=dbl)                              :: tolerance, convergence 
  REAL (KIND=dbl)                              :: ds 
  REAL (KIND=dbl)                              :: damp, temp_req

END MODULE NEB_variables


!> Module for NEB calculations (Calculate norm of a vector)
MODULE Miscellany

  USE Numeric
  
  IMPLICIT NONE
  
  CONTAINS
  
    FUNCTION norm(vect)

      IMPLICIT NONE

      REAL (KIND=dbl), DIMENSION(:), INTENT(IN)  :: vect
      REAL (KIND=dbl)                            :: norm


      norm = SQRT( DOT_PRODUCT( vect , vect ) )

    END FUNCTION norm
    
    FUNCTION file_exists( file )
    
      CHARACTER (LEN=*), INTENT(IN)   :: file
      LOGICAL                         :: file_exists
 
  
      inquire(FILE = file, EXIST = file_exists)
    
    END FUNCTION file_exists

    FUNCTION cubic_pbc( vect, Lx, Ly, Lz )

      IMPLICIT NONE    
      
      REAL (KIND=dbl), DIMENSION(:), INTENT(IN)  :: vect
      real (kind=dbl), intent(in) :: Lx, Ly, Lz
      REAL (KIND=dbl), DIMENSION( SIZE( vect ) ) :: cubic_pbc
      REAL (KIND=dbl)                            :: vect_i, invLx, invLy, invLz
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
  
END MODULE Miscellany


!> Modules which contains minimizaton routines for NEB calculation
MODULE Minimization_routines

  USE Numeric
  use module_defs
  USE Miscellany
    
  IMPLICIT NONE
  
  REAL (KIND=gp), PARAMETER  :: epsi = 1.0D-16

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

      IF ( norm(grad) >= epsi ) THEN
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

      REAL (KIND=dbl) :: gamma, norm_grad
      REAL (KIND=dbl) :: squared_old_grad_i, abs_conj_dir_i 
      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE   :: conj_dir_i

      squared_old_grad_i = DOT_PRODUCT( old_grad, old_grad ) 
      IF ( squared_old_grad_i >= epsi ) THEN
         gamma = DOT_PRODUCT( grad, grad ) / squared_old_grad_i
      ELSE
         gamma = 0.D0
      END IF

      allocate(conj_dir_i(ndim))

      norm_grad = norm(grad)
      IF ( norm_grad >= epsi ) THEN
        conj_dir_i = - grad / norm_grad + gamma * delta
      ELSE
        conj_dir_i = 0.D0
      END IF

      abs_conj_dir_i = norm( conj_dir_i ) 

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

      REAL (KIND=dbl) :: gamma, norm_grad
      REAL (KIND=dbl) :: squared_old_grad_i, abs_conj_dir_i  
      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE   :: conj_dir_i

      squared_old_grad_i = DOT_PRODUCT( old_grad , old_grad ) 
      IF ( squared_old_grad_i >= epsi ) THEN
        gamma = DOT_PRODUCT( ( grad - old_grad ) , grad ) / squared_old_grad_i
      ELSE
        gamma = 0.D0
      END IF

      allocate(conj_dir_i(ndim))

      norm_grad = norm(grad)
      IF ( norm_grad >= epsi ) THEN
        conj_dir_i = - grad / norm_grad + gamma * delta
      ELSE
        conj_dir_i = 0.D0
      END IF

      abs_conj_dir_i = norm( conj_dir_i ) 

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

      REAL (KIND=dbl), DIMENSION(ndim) :: force_versor
      REAL (KIND=dbl)                  :: vel_component

      call axpy(ndim, -ds / 2.D0, grad(1), 1, vel(1), 1)
      
      IF ( norm(grad) >= epsi ) THEN
        force_versor = - grad / norm(grad)
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
  use Miscellany

  implicit none

  private

  real(gp), dimension(:,:), allocatable :: old_grad, delta_pos, vel

  public :: init_images, deallocate_images, set_init_vel
  public :: compute_tangent, compute_gradient, compute_error, compute_neb_pos
  public :: write_restart, write_restart_vel, write_dat_files
  public :: termalization

contains

  subroutine init_images(ndim, nimages, algorithm)
    implicit none
    integer, intent(in) :: ndim, nimages, algorithm

    if (algorithm <= 3) then
       allocate(old_grad(ndim, nimages))
       old_grad = 0.
       allocate(delta_pos(ndim, nimages))
       delta_pos = 0.
    else
       allocate(vel(ndim, nimages))
       vel = 0.
    end if
  end subroutine init_images

  subroutine set_init_vel(ndim, nimages, vel0)
    implicit none
    integer, intent(in) :: ndim, nimages
    real(gp), dimension(ndim, nimages), intent(in) :: vel0

    call vcopy(ndim * nimages, vel0(1,1), 1, vel(1,1), 1)
  end subroutine set_init_vel

  subroutine deallocate_images()
    implicit none

    if (allocated(old_grad)) deallocate(old_grad)
    if (allocated(delta_pos)) deallocate(delta_pos)
    if (allocated(vel)) deallocate(vel)
  end subroutine deallocate_images

  SUBROUTINE compute_tangent(tgt, ndim, nimages, V, pos, Lx, Ly, Lz)
    IMPLICIT NONE

    integer, intent(in) :: nimages, ndim
    real(kind = gp), intent(in) :: Lx, Ly, Lz
    real(kind = gp), dimension(ndim,nimages), intent(out) :: tgt
    real(kind = gp), dimension(nimages), intent(in) :: V
    real(kind = gp), dimension(ndim, nimages), intent(in) :: pos

    INTEGER :: i   
    REAL (KIND=gp) :: V_previous, V_actual, V_next
    REAL (KIND=gp) :: abs_next, abs_previous
    REAL (KIND=gp) :: delta_V_max, delta_V_min

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
    real(kind = gp), dimension(ndim, nimages), intent(out) :: grad
    real(kind = gp), dimension(nimages), intent(in) :: V
    real(kind = gp), dimension(ndim, nimages), intent(in) :: pos, tgt, PES_grad
    real(kind = gp), intent(in) :: k_max, k_min, Lx, Ly, Lz
    logical, intent(in) :: optimization, climbing

    INTEGER         :: i, Emax_index
    REAL (KIND=gp) :: Ei, Eref, Emax
    real(kind = gp), dimension(nimages) :: k
    real(kind = gp), dimension(:), allocatable :: elastic_gradient

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
    real(kind = gp), dimension(ndim, nimages), intent(in) :: tgt, PES_grad
    real(kind = gp), dimension(ndim, nimages), intent(in) :: grad
    real(kind = gp), dimension(nimages), intent(out) :: error

    INTEGER                       :: i

    error(1) = norm(grad(:,1))
    DO i = 2, ( nimages - 1 ) 

       !! the error is given by the norm of the component of the 
       !! total energy gradient ortogonal to the path
       error(i) = norm( PES_grad(:,i) - DOT_PRODUCT( PES_grad(:,i) , tgt(:,i) ) * tgt(:,i) )

    END DO
    error(nimages) = norm(grad(:,nimages ))
  END SUBROUTINE compute_error

  subroutine compute_neb_pos(loop, ndim, nimages, iteration, error, F, V, pos, PES_grad, Lx, Ly, Lz, &
       & max_iterations, convergence, algorithm, optimization, climbing, k_min, k_max, ds, damp, temp_req)
    use module_defs
    use Minimization_routines
    use Miscellany

    implicit none

    integer, intent(in) :: ndim, nimages, max_iterations, algorithm
    integer, intent(inout) :: iteration
    logical, intent(in) :: optimization, climbing
    real(gp), intent(in) :: Lx, Ly, Lz, k_min, k_max, ds, damp, temp_req, convergence
    real(gp), dimension(nimages), intent(in) :: V
    real(gp), dimension(ndim,nimages), intent(in) :: PES_grad
    real(gp), dimension(ndim,nimages), intent(inout) :: pos
    real(gp), dimension(nimages), intent(out) :: error, F
    logical, intent(out) :: loop

    integer :: i, n_in, n_fin
    real(gp) :: err
    real(gp), dimension(:,:), allocatable :: tangent, grad
    CHARACTER (LEN=4), PARAMETER :: exit_file = "EXIT"  

    loop = .true.
    IF ( optimization ) THEN
       N_in  = 1
       N_fin = nimages
    ELSE
       N_in  = 2
       N_fin = nimages - 1
    END IF
    allocate(tangent(ndim, nimages), grad(ndim, nimages))

    CALL compute_tangent(tangent, ndim, nimages, V, pos, Lx, Ly, Lz)

    if (iteration > 0 .and. algorithm >= 4 ) THEN

       CALL compute_gradient(ndim, nimages, grad, &
            & V, pos, tangent, PES_grad, Lx, Ly, Lz, &
            & k_min, k_max, optimization, algorithm, climbing)
       second_minimization_loop: DO i = N_in, N_fin 
          IF ( algorithm == 4 ) THEN
             CALL quick_min_second_step( ndim, vel(:,i), grad(:,i), ds)
          ELSE IF ( algorithm == 5 ) THEN
             CALL velocity_Verlet_second_step( ndim, vel(:,i), grad(:,i), ds, damp )
          ELSE IF ( algorithm == 6 ) THEN
             CALL velocity_Verlet_second_step( ndim, vel(:,i), grad(:,i), ds )
          END IF
       END DO second_minimization_loop
       IF ( algorithm == 6 ) CALL termalization(ndim, nimages, vel, temp_req)

    end if

    if (iteration > 0 .or. max_iterations == 1) then
       CALL compute_error( error, ndim, nimages, grad, tangent, PES_grad )
       err = maxval(error)

       F = 0.D0
       DO i = 2, ( nimages - 1 )
          F(i) = DOT_PRODUCT( - PES_grad(:,i) , tangent(:,i) )
       END DO

       IF ( ( err * E_zero / a_zero ) <= convergence .or. max_iterations == 1)  THEN
          loop = .false.
          return
       END IF
    end if

    iteration = iteration + 1

    IF ( iteration > max_iterations ) THEN
       loop = .false.
       return
    END IF

    IF ( file_exists(exit_file) ) THEN
       CALL SYSTEM ("rm -f EXIT")

       WRITE(*,*) " WARNING :  soft exit required"
       WRITE(*,*) " STOPPING ...                 "

       loop = .false.
       return
    END IF

    IF ( ALGORITHM <= 3 ) old_grad = grad
    CALL compute_gradient(ndim, nimages, grad, V, pos, &
         & tangent, PES_grad, Lx, Ly, Lz, &
         & k_min, k_max, optimization, algorithm, climbing)

    ! Calculate new positions for next step.
    first_minimization_loop: DO i = N_in, N_fin
       IF ( algorithm == 1 ) THEN
          CALL steepest_descent( ndim, pos(:,i), grad(:,i), ds)
       ELSE IF ( algorithm == 2 ) THEN
          CALL fletcher_reeves( ndim, pos(:,i), grad(:,i), old_grad(:,i), delta_pos(:,i), ds)
       ELSE IF ( algorithm == 3 ) THEN
          CALL polak_ribiere( ndim, pos(:,i), grad(:,i), old_grad(:,i), delta_pos(:,i), ds)
       ELSE IF ( algorithm >= 4 ) THEN
          CALL velocity_Verlet_first_step( ndim, pos(:,i), vel(:,i), grad(:,i), ds)
       END IF
    END DO first_minimization_loop

    deallocate(tangent, grad)
  END SUBROUTINE compute_neb_pos
  
  SUBROUTINE termalization(ndim, nimages, vel, temp_req)
    IMPLICIT NONE

    integer, intent(in) :: ndim, nimages
    real(gp), intent(in) :: temp_req
    real(gp), dimension(ndim, nimages), intent(inout) :: vel

    REAL (KIND=dbl)  :: temp
    INTEGER          :: i

    temp = 0.D0  
    DO i = 2, ( nimages - 1 )
       temp = temp + DOT_PRODUCT( vel(:,i) , vel(:,i) )
    END DO

    temp = temp / DBLE( ( nimages - 2 ) * ndim )
    vel = vel * SQRT( temp_req / temp )

  END SUBROUTINE termalization

  SUBROUTINE write_restart(restart_file, ndim, nimages, V, pos, fix_atom, PES_grad)
    use Formats

    IMPLICIT NONE

    integer, intent(in) :: ndim, nimages
    character(len = *), intent(in) :: restart_file
    real(gp), dimension(nimages), intent(in) :: V
    real(gp), dimension(ndim, nimages), intent(in) :: pos, PES_grad
    integer, dimension(ndim), intent(in) :: fix_atom

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
               fix_atom(j),     &
               fix_atom((j+1)), &
               fix_atom((j+2)), &
               -PES_grad(j,i),     &
               -PES_grad((j+1),i), &
               -PES_grad((j+2),i)
       END DO
    END DO
    CLOSE( UNIT = unit )
  END SUBROUTINE write_restart

  SUBROUTINE write_restart_vel(velocity_file)
    use Formats

    IMPLICIT NONE

    character(len = *), intent(in) :: velocity_file

    INTEGER             :: i, j
    INTEGER, PARAMETER  :: unit = 10

    OPEN( UNIT = unit, FILE = trim(velocity_file), STATUS = "UNKNOWN", ACTION = "WRITE" )
    DO i = 1, size(vel, 2)
       WRITE(unit,*) "Replica: ", i

       DO j = 1, size(vel, 1), 3 
          WRITE(unit,fmt2) vel(j,i), vel((j+1),i), vel((j+2),i)
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
    REAL (KIND=dbl)                            :: R, delta_R, x
    REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: d_R
    REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: a, b, c, d
    REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: reaction_coordinate
    REAL (KIND=dbl)                            :: E, E_0
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
       WRITE(unit,'(3(2X,F12.8))') reaction_coordinate(i) * a_zero, &
            ( V(i) - V(1) ) * E_zero, error(i) * ( E_zero / a_zero )
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

       WRITE(unit,'(2(2X,F16.8))') ( R * a_zero ), &
            ( E - E_0 ) * E_zero
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

!> Module for NEB calculations
MODULE NEB_routines

  USE Numeric
  USE Formats
  USE NEB_variables
  USE Miscellany
  USE Minimization_routines
  use module_images
  
  IMPLICIT NONE

  CONTAINS

    SUBROUTINE read_input

      use module_types
      use module_interfaces

      IMPLICIT NONE

      INTEGER                                    :: i, j, n1, n2, ios, j0, k
      real(kind = dbl)                           :: r, a
      CHARACTER (LEN=20)                         :: minimization_scheme
      CHARACTER (LEN=95)                         :: vel_file
      INTEGER, PARAMETER                         :: unit = 10
      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: d_R
      REAL (KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: vel0
      type(atoms_data)                           :: at
      real(kind = dbl), dimension(:,:), pointer  :: rxyz
      real(kind = dbl), dimension(3, 1000)       :: xcart1, xcart2
      real(kind = dbl), dimension(3)             :: acell1, acell2

      NAMELIST /NEB/ first_config,        &
                     last_config,         &
         scratch_dir,         &
         job_name,            &
         restart,             &
         climbing,            &
         optimization,        &
         minimization_scheme, &
         damp,                &
         temp_req,            &
         k_max, k_min,        &
         ds,                  &
         max_iterations,      &
         tolerance,           &
         convergence,         &
         num_of_images


!! default values are assigned

      scratch_dir        = "./"

      job_name           = "neb"
      restart      = .FALSE.
      climbing     = .FALSE.
      optimization = .FALSE.
      
      minimization_scheme = "quick-min"
      damp                = 1.D0
      temp_req            = 0.D0
      
      k_max = 0.1D0
      k_min = 0.1D0
      
      ds = 0.5D0
      
      max_iterations = 1
      
      tolerance   = 1.0D-4
      convergence = 5.0D-2

      N = 0
      
      Lx = 0.D0
      Ly = 0.D0
      Lz = 0.D0

      READ( * , NML=NEB )

      barrier_file       = trim(job_name) // ".NEB.log"
      data_file          = trim(job_name) // ".NEB.dat"
      interpolation_file = trim(job_name) // ".NEB.int"
      restart_file       = trim(job_name) // ".NEB.restart"

!! initial and final configuration are read only if a new simulation
!! is started ( restart = .FALSE. ) 
      
      IF ( .NOT. restart ) THEN

         IF ( trim(first_config) == "" ) THEN

            WRITE(*,'(T2,"read_input: first_config not assigned ")') 

            STOP

         ELSE IF ( trim(last_config) == "" ) THEN       

            WRITE(*,'(T2,"read_input: last_config not assigned ")') 

            STOP

         END IF

      END IF
 
!!!      IF ( restart ) THEN
!!!         num_of_images = num_of_images - 2
!!!      END IF

      IF ( minimization_scheme == "steepest_descent" ) THEN

         algorithm = 1

      ELSE IF ( minimization_scheme == "fletcher-reeves" ) THEN

         algorithm = 2

      ELSE IF ( minimization_scheme == "polak-ribiere" ) THEN

         algorithm = 3

      ELSE IF ( minimization_scheme == "quick-min" ) THEN

         algorithm = 4

      ELSE IF ( minimization_scheme == "damped-verlet" ) THEN

         algorithm = 5

      ELSE IF ( minimization_scheme == "sim-annealing" ) THEN

         algorithm = 6

      ELSE

         WRITE(*,'(T2,"read_input: minimization_scheme ", A20)') &
              minimization_scheme
         WRITE(*,'(T2,"            does not exist")') 
         STOP 

      END IF

      call read_atomic_file(trim(first_config), 0, at, rxyz)
      n1 = at%nat
      acell1(1) = at%alat1
      acell1(2) = at%alat2
      acell1(3) = at%alat3
      xcart1(:,1:at%nat) = rxyz
      if (acell1(1) == 0.) acell1(1) = maxval(rxyz(1,:)) - minval(rxyz(1,:))
      if (acell1(2) == 0.) acell1(2) = maxval(rxyz(2,:)) - minval(rxyz(2,:))
      if (acell1(3) == 0.) acell1(3) = maxval(rxyz(3,:)) - minval(rxyz(3,:))
      deallocate(rxyz)

      call read_atomic_file(trim(last_config), 0, at, rxyz)
      n2 = at%nat
      acell2(1) = at%alat1
      acell2(2) = at%alat2
      acell2(3) = at%alat3
      xcart2(:,1:at%nat) = rxyz
      if (acell2(1) == 0.) acell2(1) = maxval(rxyz(1,:)) - minval(rxyz(1,:))
      if (acell2(2) == 0.) acell2(2) = maxval(rxyz(2,:)) - minval(rxyz(2,:))
      if (acell2(3) == 0.) acell2(3) = maxval(rxyz(3,:)) - minval(rxyz(3,:))
      deallocate(rxyz)

      if (at%geocode == 'F') then
        acell1 = max(acell1, acell2)
        acell2 = acell1
      end if

!! some consistency checks are done

      IF ( num_of_images <= 2 ) THEN

        WRITE(*,'(T1,"read_input: num_of_images must be larger than 2")')
        STOP

      END IF

      IF ( maxval(abs(acell2 - acell1)) > 1.d-6 ) THEN
      
         WRITE(*,'(T2,"read_input: box size is not constant")')
         WRITE(*,'(T2,"           dLx = ", F10.6 )') acell1(1) - acell2(1)
         WRITE(*,'(T2,"           dLy = ", F10.6 )') acell1(2) - acell2(2)
         WRITE(*,'(T2,"           dLz = ", F10.6 )') acell1(3) - acell2(3)
         STOP  

      END IF

      IF ( n1 /= n2 ) THEN
      
         WRITE(*,'(T2,"read_input: number of atoms is not constant")')
         WRITE(*,'(T2,"            N = ", I8, I8 )') n1, n2
         STOP  

      END IF

      N = n1
      Lx = acell1(1)
      Ly = acell1(2)
      Lz = acell1(3)

      dim = 3 * N
      
      CALL dyn_allocation

!! all the arrays are initialized

      V                = 0.D0
      PES_gradient     = 0.D0
      error            = 0.D0
      !TD k                = k_min

      IF ( restart ) THEN

        Emax = - 1.0D16

        OPEN( UNIT = unit, FILE = restart_file, STATUS = "OLD", &
              ACTION = "READ" )
    
        ! Read as many configurations as contained in the
        ! Restart file.
        DO i = 1, num_of_images

           READ(unit,*, iostat = ios)
           if (ios /= 0) exit
           READ(unit,fmt3) V(i)

           DO j = 1, dim, 3 

              READ(unit,*)    pos(j,i),     & 
                   pos((j+1),i), &
                   pos((j+2),i), &
                   fix_atom(j),     &
                   fix_atom((j+1)), &
                   fix_atom((j+2)), &
                   PES_gradient(j,i),     &
                   PES_gradient((j+1),i), &
                   PES_gradient((j+2),i)
           END DO
           PES_gradient(:, i) = PES_gradient(:, i) * (-1d0)

           IF ( V(i) >= Emax ) THEN

              Emax = V(i)         

              Emax_index = i

           END IF

        END DO

        CLOSE( UNIT = unit) 

        ! Check that we read at least two configurations...
        i = i - 1
        if (i < 2) then
           WRITE(*,'(T2,"read_input: number of replica in restart file is less than 2")')
           WRITE(*,'(T2,"            N = ", I8 )') i - 1
           STOP  
        end if

        ! Add some intermediate configurations if num_of_images > i
        if (num_of_images > i) then
           ALLOCATE( d_R(dim) )           
           r = real(i - 1, dbl) / real(num_of_images - 1, dbl)
           a = real(0, dbl)
           j0 = num_of_images
           do j = num_of_images, 1, -1
              ! j : the current position in pos array.
              ! i : the current position in read configurations.
              ! r : the ideal ratio (i_0 - 1) / (num_of_images - 1)
              ! a : the current ratio.
              ! We copy i to j if a < r, otherwise we create a new
              ! linear approximant.
              if (a <= r) then
                 pos(:, j) = pos(:, i)
                 i = i - 1
                 ! We create the new linear approx. replicas between j and j0.
                 if (j0 /= j) then
                    d_R = ( pos(:,j0) - pos(:,j) ) / &
                         DBLE( j0 - j )
                    do k = j0 - 1, j + 1, -1
                       write(*,*) k, ": new linear approx replica."
                       pos(:, k) = pos(:, k + 1) - d_R(:)
                    end do
                 end if
                 j0 = j
                 write(*,*) j, ": restart replica."
              end if
              a = (r * real(num_of_images - 1, dbl) - real(i, dbl) + real(1, dbl)) / &
                   & (real(num_of_images, dbl) - real(j, dbl) + real(1, dbl))
           end do
           deallocate(d_R)

           CALL write_restart(restart_file, dim, num_of_images, V, pos, fix_atom, PES_gradient)

        end if

        vel_file = TRIM( scratch_dir )//"/velocities_file"
        IF ( ( algorithm >= 4 ) .AND. ( file_exists( vel_file ) ) ) THEN
!!DEBUG          PRINT *, "reading ", vel_file
           allocate(vel0(dim, num_of_images))
           OPEN( UNIT = unit, FILE = vel_file, STATUS = "OLD", ACTION = "READ" )
           DO i = 1, num_of_images
              READ(unit,*)
              DO j = 1, dim, 3 
                 READ(unit,fmt2) vel0(j,i),     & 
                      vel0((j+1),i), &
                      vel0((j+2),i)
              END DO
           END DO
           CLOSE( UNIT = unit )
           call set_init_vel(dim, num_of_images, vel0)
           deallocate(vel0)
        END IF
      
      ELSE

         ! TODO, read XYZ here
      
        ALLOCATE( d_R(dim) )           

!!!        OPEN( UNIT = unit, FILE = first_config, STATUS = "OLD", &
!!!              ACTION = "READ" )
!!!    
!!!          DO i = 1, dim, 3 
!!!       
!!!            READ(unit,fmt2) pos(i,1),     & 
!!!                            pos((i+1),1), &
!!!                            pos((i+2),1)
!!!      
!!!          END DO
!!!
!!!        CLOSE( UNIT = unit )
        pos(:, 1) = reshape(xcart1, (/ dim /))

!!!        OPEN( UNIT = unit, FILE = last_config, STATUS = "OLD", &
!!!              ACTION = "READ" )
!!!
!!!          DO i = 1, dim, 3 
!!! 
!!!            READ(unit,fmt2) pos(i,num_of_images),     &
!!!                            pos((i+1),num_of_images), &
!!!                            pos((i+2),num_of_images)
!!!      
!!!          END DO
!!!
!!!        CLOSE( UNIT = unit )
        pos(:, num_of_images) = reshape(xcart2, (/ dim /))
  
        d_R = ( pos(:,num_of_images) - pos(:,1) ) / &
             DBLE( num_of_images - 1 )

        fix_atom = 1

        WHERE ( ABS( d_R ) <=  tolerance ) fix_atom = 0

        DO i = 2, ( num_of_images - 1 )

           pos(:,i) = pos(:,( i - 1 )) + d_R(:)

        END DO

        DEALLOCATE( d_R )

     END IF

   END SUBROUTINE read_input

    
    SUBROUTINE dyn_allocation

      IMPLICIT NONE
     
      ALLOCATE( pos( dim , num_of_images ) )

      ALLOCATE( PES_gradient( dim , num_of_images ) )          

      ALLOCATE( fix_atom( dim ) )

      ALLOCATE( F( num_of_images ) )
      ALLOCATE( V( num_of_images ) )
      ALLOCATE( error( num_of_images ) )

      call init_images(dim, num_of_images, algorithm)

    END SUBROUTINE dyn_allocation     


    SUBROUTINE search_MEP

      IMPLICIT NONE

      INTEGER         :: iteration
      LOGICAL         :: stat


      open(unit = 456, file = trim(barrier_file), action = "WRITE")
      write(456, "(A)") "# NEB barrier file"
      close(unit = 456)

      IF ( .NOT. restart) THEN
         CALL write_restart(restart_file, dim, num_of_images, V, pos, fix_atom, PES_gradient)
      END IF

      iteration = 0
      minimization: do
         CALL PES_IO(optimization .or. (.not. restart .and. iteration == 0),stat)
         if (.not. stat) exit minimization

         call compute_neb_pos(stat, dim, num_of_images, iteration, error, F, V, pos, PES_gradient, &
              & Lx, Ly, Lz, &
              & max_iterations, convergence, algorithm, optimization, climbing, &
              & k_min, k_max, ds, damp, temp_req)

         CALL write_restart(restart_file, dim, num_of_images, V, pos, fix_atom, PES_gradient)
         if (algorithm >= 4) then
            CALL write_restart_vel(trim(scratch_dir) // "velocities_file")
         end if

         if (iteration > 1) then
            CALL write_dat_files(data_file, interpolation_file, dim, num_of_images, pos, V, F, error, Lx, Ly, Lz)
            open(unit = 456, file = trim(barrier_file), action = "WRITE", position = "APPEND")
            WRITE(456, fmt4) iteration - 1, ( maxval(V(2:num_of_images - 1)) - V(1) ) * E_zero, &
                 & maxval(error) * ( E_zero / a_zero ) 
            WRITE(*, fmt4)   iteration - 1, ( maxval(V(2:num_of_images - 1)) - V(1) ) * E_zero, &
                 & maxval(error) * ( E_zero / a_zero ) 
            close(unit = 456)
         end if

         if (.not. stat) exit minimization
      end do minimization
    END SUBROUTINE search_MEP

    SUBROUTINE PES_IO( flag , stat )

      IMPLICIT NONE

      LOGICAL, INTENT(IN)        :: flag
      LOGICAL, INTENT(OUT)       :: stat
      INTEGER                    :: i, replica
      INTEGER                    :: N_in, N_fin
      REAL (KIND=dbl)            :: temp_V
      REAL (KIND=dbl), PARAMETER :: corruption_flag = 9999999.99999999
      INTEGER, PARAMETER         :: unit = 10     
      REAL (KIND=dbl), PARAMETER :: epsi = 1.0D-8

      IF ( flag ) THEN

         CALL SYSTEM( "./NEB_driver.sh all " // trim(job_name) // &
              & " " // trim(scratch_dir) // " " // trim(first_config))

        N_in  = 1
        N_fin = num_of_images

      ELSE
         
         CALL SYSTEM( "./NEB_driver.sh free_only " // trim(job_name) // &
              & " " // trim(scratch_dir) // " " // trim(first_config))

        N_in  = 2
        N_fin = ( num_of_images - 1 )

      END IF

      Emax = - 1.0D16
      stat = .TRUE.

      OPEN( UNIT = unit, FILE = "gen_output_file", STATUS = "OLD", &
            ACTION = "READ" )

        DO replica = N_in, N_fin

          READ(unit,*) temp_V

          IF ( ABS( temp_V - corruption_flag ) <= epsi ) THEN
                 
      stat = .FALSE. 
      
      RETURN
    
    END IF
    
    V(replica) = temp_V

          IF ( V(replica) >= Emax ) THEN
 
            Emax = V(replica)         

            Emax_index = replica

          END IF

          DO i = 1, dim, 3 

            READ(unit,*) PES_gradient(i,replica),     &
                 PES_gradient((i+1),replica), &
                 PES_gradient((i+2),replica)

          END DO

          PES_gradient(:,replica) = PES_gradient(:,replica) * &
                                    DBLE( fix_atom ) * (-1.d0)

        END DO

      CLOSE( UNIT = unit )   

    END SUBROUTINE PES_IO

    SUBROUTINE write_output

      IMPLICIT NONE

      INTEGER  :: i


      DO i = 1, num_of_images

        WRITE(*,fmt5) i, V(i) * E_zero, error(i) * ( E_zero / a_zero )

      END DO

    END SUBROUTINE write_output


    SUBROUTINE deallocation

      IMPLICIT NONE

      IF ( ALLOCATED( pos ) )              DEALLOCATE( pos )
      IF ( ALLOCATED( PES_gradient ) )     DEALLOCATE( PES_gradient )
      IF ( ALLOCATED( fix_atom ) )         DEALLOCATE( fix_atom )
      IF ( ALLOCATED( F ) )                DEALLOCATE( F )
      IF ( ALLOCATED( V ) )                DEALLOCATE( V )
      IF ( ALLOCATED( error ) )            DEALLOCATE( error )
      
      call deallocate_images()
    END SUBROUTINE deallocation

END MODULE NEB_routines

PROGRAM NEB

  USE NEB_routines

  IMPLICIT NONE

  CALL read_input

  CALL search_MEP

  CALL write_output

  CALL deallocation

  STOP

END PROGRAM NEB
