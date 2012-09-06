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
  REAL (KIND=dbl)                              :: Lx, Ly, Lz, &
                                                  invLx, invLy, invLz
  INTEGER                                      :: algorithm
  INTEGER                                      :: num_of_images
  REAL (KIND=dbl)                              :: k_max, k_min, delta_k
  REAL (KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: pos, delta_pos
  REAL (KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: vel
  REAL (KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: grad, old_grad
  REAL (KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: PES_gradient
  REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE   :: norm_grad 
  INTEGER, DIMENSION(:), ALLOCATABLE           :: fix_atom
  REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE   :: V, k, error
  REAL (KIND=dbl)                              :: Eref, Emax
  INTEGER(KIND=sgl)                            :: Emax_index
  REAL (KIND=dbl), DIMENSION(:,:), ALLOCATABLE :: tangent
  REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE   :: elastic_gradient
  REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE   :: conj_dir_i

  INTEGER                                      :: max_iterations 
  REAL (KIND=dbl)                              :: tolerance, convergence 
  REAL (KIND=dbl)                              :: ds 
  REAL (KIND=dbl)                              :: damp, temp_req

!! parameters

  CHARACTER (LEN=4), PARAMETER :: exit_file = "EXIT"  
  REAL (KIND=dbl), PARAMETER   :: a_zero = 0.529177D0 ! from Bohr to angstrom
  REAL (KIND=dbl), PARAMETER   :: E_zero = 27.212D0   ! from a.u. to eV 

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
  
END MODULE Miscellany


!> Modules which contains minimizaton routines for NEB calculation
MODULE Minimization_routines

  USE Numeric
  USE NEB_variables
  USE Miscellany
    
  IMPLICIT NONE
  
  CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!! minimization algorithms !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! steepest descent with line minimization !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE steepest_descent( index )

      IMPLICIT NONE

      INTEGER, INTENT(IN)         :: index
      REAL (KIND=dbl), PARAMETER  :: epsi = 1.0D-16 
      
      IF ( norm_grad(index) >= epsi ) THEN

        delta_pos(:,index) = grad(:,index) / norm_grad(index)
 
      ELSE

        delta_pos(:,index) = 0.D0

      END IF

      pos(:,index) = pos(:,index) - &
                     ds * norm_grad(index) * delta_pos(:,index)

    END SUBROUTINE steepest_descent 

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! conjugate gradient minimization !!
    !! Fletcher - Reeves               !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE fletcher_reeves( index )
     
      IMPLICIT NONE

      INTEGER (KIND=sgl), INTENT(IN)            :: index
      REAL (KIND=dbl)                           :: gamma
      REAL (KIND=dbl)                           :: squared_old_grad_i, &
                                                   abs_conj_dir_i 
      REAL (KIND=dbl), PARAMETER                :: epsi = 1.0D-16 


      squared_old_grad_i = DOT_PRODUCT( old_grad(:,index) , old_grad(:,index) ) 

      IF ( squared_old_grad_i >= epsi ) THEN

        gamma = DOT_PRODUCT( grad(:,index) , grad(:,index) ) / &
          squared_old_grad_i

      ELSE

        gamma = 0.D0

      END IF

      IF ( norm_grad(index) >= epsi ) THEN

        conj_dir_i = - grad(:,index) / norm_grad(index) + &
                       gamma * delta_pos(:,index)

      ELSE
      
        conj_dir_i = 0.D0
      
      END IF

      abs_conj_dir_i = norm( conj_dir_i ) 

      IF ( abs_conj_dir_i >= epsi ) THEN

        delta_pos(:,index) = conj_dir_i / abs_conj_dir_i

      ELSE

        delta_pos(:,index) = 0.D0

      END IF

      pos(:,index) = pos(:,index) + &
                     ds * norm_grad(index) * conj_dir_i 

    END SUBROUTINE fletcher_reeves  

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! conjugate gradient minimization !!
    !! Polak - Ribiere                 !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE polak_ribiere( index )

      IMPLICIT NONE

      INTEGER (KIND=sgl), INTENT(IN)            :: index
      REAL (KIND=dbl)                           :: gamma
      REAL (KIND=dbl)                           :: squared_old_grad_i, &
                                                   abs_conj_dir_i  
      REAL (KIND=dbl), PARAMETER                :: epsi = 1.0D-16 


      squared_old_grad_i = DOT_PRODUCT( old_grad(:,index) , old_grad(:,index) ) 

      IF ( squared_old_grad_i >= epsi ) THEN

        gamma = DOT_PRODUCT( ( grad(:,index) - old_grad(:,index) ) , &
          grad(:,index) ) / squared_old_grad_i

      ELSE

        gamma = 0.D0

      END IF

      IF ( norm_grad(index) >= epsi ) THEN

        conj_dir_i = - grad(:,index) / norm_grad(index) + &
                       gamma * delta_pos(:,index)

      ELSE
      
        conj_dir_i = 0.D0
      
      END IF

      abs_conj_dir_i = norm( conj_dir_i ) 

      IF ( abs_conj_dir_i >= epsi ) THEN

        delta_pos(:,index) = conj_dir_i / abs_conj_dir_i

      ELSE

        delta_pos(:,index) = 0.D0

      END IF

      pos(:,index) = pos(:,index) + &
                     ds * norm_grad(index) * conj_dir_i 

    END SUBROUTINE  polak_ribiere

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Molecular Dynamics based algorithms !!
    !! velocity Verlet and quick min       !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE velocity_Verlet_first_step( index )

      IMPLICIT NONE

      INTEGER (KIND=sgl), INTENT(IN) :: index
      
         
      vel(:,index) = vel(:,index) - ds / 2.D0 * grad(:,index)
      pos(:,index) = pos(:,index) + ds * vel(:,index)
      
    END SUBROUTINE velocity_Verlet_first_step
    
    
    SUBROUTINE velocity_Verlet_second_step( index )

      IMPLICIT NONE

      INTEGER (KIND=sgl), INTENT(IN) :: index
    
      IF ( algorithm == 5 ) THEN

        vel(:,index) = damp * ( vel(:,index) - ds / 2.D0 * grad(:,index) )
      
      ELSE IF ( algorithm == 6 ) THEN
              
  vel(:,index) = vel(:,index) - ds / 2.D0 * grad(:,index)

      END IF
      
    END SUBROUTINE velocity_Verlet_second_step
    
    
    SUBROUTINE quick_min_second_step( index )

      IMPLICIT NONE

      INTEGER (KIND=sgl), INTENT(IN)   :: index
      REAL (KIND=dbl), DIMENSION(dim)  :: force_versor
      REAL (KIND=dbl)                  :: vel_component
      REAL (KIND=dbl), PARAMETER       :: epsi = 1.0D-16


      vel(:,index) = vel(:,index) - ds / 2.D0 * grad(:,index)
      
      IF ( norm_grad(index) >= epsi ) THEN

        force_versor = - grad(:,index) / norm_grad(index)

      ELSE
      
        force_versor = 0.D0
      
      END IF
 
      vel_component = DOT_PRODUCT( vel(:,index) , force_versor )
      
      IF ( vel_component > 0.D0 ) THEN
      
        vel(:,index) = vel_component * force_versor
  
      ELSE
      
        vel(:,index) = 0.D0
  
      END IF    
      
    END SUBROUTINE quick_min_second_step
    
    
    SUBROUTINE termalization
     
      IMPLICIT NONE
      
      REAL (KIND=dbl)  :: temp
      INTEGER          :: i
      
      temp = 0.D0  

      DO i = 2, ( num_of_images - 1 )
       
        temp = temp + DOT_PRODUCT( vel(:,i) , vel(:,i) )

      END DO
  
      temp = temp / DBLE( ( num_of_images - 2 ) * dim ) 
   
      vel = vel * SQRT( temp_req / temp )  

    END SUBROUTINE termalization
    
END MODULE Minimization_routines


!> Module for NEB calculations
MODULE NEB_routines

  USE Numeric
  USE Formats
  USE NEB_variables
  USE Miscellany
  USE Minimization_routines
  
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
      
      invLx = 1.D0 / Lx
      invLy = 1.D0 / Ly
      invLz = 1.D0 / Lz

      CALL dyn_allocation

!! all the arrays are initialized

      V                = 0.D0
      PES_gradient     = 0.D0
      elastic_gradient = 0.D0
      tangent          = 0.D0
      grad             = 0.D0
      norm_grad        = 0.D0
      error            = 0.D0
      !TD k                = k_min

      delta_k = k_max - k_min

      IF ( algorithm <= 3 ) THEN

         delta_pos  = 0.D0
         old_grad   = 0.D0
         conj_dir_i = 0.D0

      ELSE

         vel = 0.D0
      
      END IF

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

           CALL write_restart     

        end if

        vel_file = TRIM( scratch_dir )//"/velocities_file"
  
        IF ( ( algorithm >= 4 ) .AND. &
       ( file_exists( vel_file ) ) ) THEN

!!DEBUG          PRINT *, "reading ", vel_file
      
          OPEN( UNIT = unit, FILE = vel_file, STATUS = "OLD", ACTION = "READ" )
    
            DO i = 1, num_of_images

              READ(unit,*)

              DO j = 1, dim, 3 
       
                READ(unit,fmt2) vel(j,i),     & 
                                vel((j+1),i), &
                                vel((j+2),i)
  
              END DO

            END DO

          CLOSE( UNIT = unit )
      
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
      
      IF ( algorithm <= 3 ) THEN
        
  ALLOCATE( delta_pos( dim , num_of_images ) )
      
      ELSE
      
        ALLOCATE( vel( dim , num_of_images ) )     
  
      END IF
      
      ALLOCATE( grad( dim , num_of_images ) )
      ALLOCATE( PES_gradient( dim , num_of_images ) )          
      
      IF ( algorithm <= 3 ) THEN
             
        ALLOCATE( old_grad( dim , num_of_images ) )    
  
      END IF    

      ALLOCATE( fix_atom( dim ) )
      
      ALLOCATE( norm_grad( num_of_images ) )

      ALLOCATE( V( num_of_images ) )
      ALLOCATE( k( num_of_images ) )
      ALLOCATE( error( num_of_images ) )

      ALLOCATE( tangent( dim , num_of_images ) ) 
      ALLOCATE( elastic_gradient( dim ) )

      IF ( algorithm <= 3 ) ALLOCATE( conj_dir_i( dim ) )           
      
            
    END SUBROUTINE dyn_allocation     


    SUBROUTINE search_MEP

      IMPLICIT NONE

      REAL (KIND=dbl) :: err 
      INTEGER         :: i, iteration
      INTEGER         :: N_in, N_fin
      LOGICAL         :: stat


      open(unit = 456, file = trim(barrier_file), action = "WRITE")
      write(456, "(A)") "# NEB barrier file"
      close(unit = 456)

      IF ( .NOT. restart) THEN

         CALL write_restart     

         CALL PES_IO(.TRUE.,stat)

         IF ( .NOT. stat ) THEN

            !!DEBUG    WRITE(*,*) " FATAL: corruption in the gen_output_file"
            !!DEBUG    WRITE(*,*) " STOPPING ...                            "

            RETURN

         END IF

         CALL compute_tangent

         IF ( max_iterations == 1 ) THEN

            CALL write_restart    

            CALL compute_error(err)

            CALL write_dat_files  

            open(unit = 456, file = trim(barrier_file), action = "WRITE", position = "APPEND")

            WRITE(456, fmt4) &
                 1, ( Emax - V(1) ) * E_zero, err * ( E_zero / a_zero ) 
            WRITE(*, fmt4) &
                 1, ( Emax - V(1) ) * E_zero, err * ( E_zero / a_zero ) 

            close(unit = 456)

            RETURN

         END IF

      ELSE

         if (maxval(PES_gradient(:, 1)) == 0.d0 .OR. &
              & maxval(PES_gradient(:, num_of_images)) == 0.d0) then
            CALL PES_IO(.TRUE.,stat)

            IF ( .NOT. stat ) THEN

               !!DEBUG    WRITE(*,*) " FATAL: corruption in the gen_output_file"
               !!DEBUG    WRITE(*,*) " STOPPING ...                            "

               RETURN

            END IF
         end if

         CALL compute_tangent   

         CALL write_restart     

      END IF

      IF ( optimization ) THEN

         N_in  = 1
         N_fin = num_of_images

      ELSE

         N_in  = 2
         N_fin = num_of_images - 1

      END IF

      iteration = 0

      minimization: DO

         iteration = iteration + 1

         IF ( iteration > max_iterations ) THEN

            CALL write_restart

            EXIT minimization

         END IF

         IF ( file_exists(exit_file) ) THEN

             CALL SYSTEM ("rm -f EXIT")

            WRITE(*,*) " WARNING :  soft exit required"
            WRITE(*,*) " STOPPING ...                 "

            CALL write_restart    

            EXIT minimization

         END IF

         CALL gradient

         first_minimization_loop: DO i = N_in, N_fin

            IF ( algorithm == 1 ) THEN

               CALL steepest_descent(i)

            ELSE IF ( algorithm == 2 ) THEN

               CALL fletcher_reeves(i)

            ELSE IF ( algorithm == 3 ) THEN

               CALL polak_ribiere(i)

            ELSE IF ( algorithm >= 4 ) THEN

               CALL velocity_Verlet_first_step(i)

            END IF

         END DO first_minimization_loop

         CALL write_restart

         CALL PES_IO(optimization,stat)

         IF ( .NOT. stat ) THEN

            !!DEBUG    WRITE(*,*) " FATAL: corruption in the gen_output_file"
            !!DEBUG    WRITE(*,*) " STOPPING ...                            "

            EXIT minimization

         END IF

         CALL compute_tangent

         IF ( algorithm >= 4 ) THEN

            CALL gradient

            second_minimization_loop: DO i = N_in, N_fin 

               IF ( algorithm == 4 ) THEN

                  CALL quick_min_second_step(i)

               ELSE IF ( algorithm >= 5 ) THEN

                  CALL velocity_Verlet_second_step(i)

               END IF

            END DO second_minimization_loop

            IF ( algorithm == 6 ) CALL termalization

         END IF

         CALL compute_error(err)

         CALL write_dat_files

         open(unit = 456, file = trim(barrier_file), action = "WRITE", position = "APPEND")

         WRITE(456, fmt4) &
              iteration, ( Emax - V(1) ) * E_zero, err * ( E_zero / a_zero ) 
         WRITE(*, fmt4) &
              iteration, ( Emax - V(1) ) * E_zero, err * ( E_zero / a_zero ) 

         close(unit = 456)

         IF ( ( err * E_zero / a_zero ) <= convergence )  THEN

            CALL write_restart

            EXIT minimization

         END IF

      END DO minimization

    END SUBROUTINE search_MEP

    
    SUBROUTINE compute_tangent
    
      IMPLICIT NONE
      
      INTEGER :: i   
      
    
      tangent = 0
    
      DO i = 2, ( num_of_images - 1 )

!! tangent to the path (normalized)

        tangent(:,i) = path_tangent(i)
        tangent(:,i) = tangent(:,i)/ norm( tangent(:,i) )
    
      END DO
    
    END SUBROUTINE compute_tangent
    
    
    SUBROUTINE gradient

      IMPLICIT NONE
      
      INTEGER         :: i
      REAL (KIND=dbl) :: Ei


      Eref = MIN( V(1) , V(num_of_images) )

      elastic_const_loop: DO i = 2, num_of_images 

        IF ( i < num_of_images ) THEN

          Ei = MAX( MAX( V(i) , V(i-1) ) , MAX( V(i) , V(i+1) ) )

        ELSE
  
          Ei = MAX( V(i) , V(i-1) )  
  
  END IF

        IF ( Ei > Eref ) THEN

          k(i) = k_max - delta_k * ( ( Emax - Ei ) / ( Emax - Eref ) )

        ELSE

          k(i) = k_min

        END IF
 
      END DO elastic_const_loop

      IF ( ALGORITHM <= 3 ) old_grad = grad

      IF ( optimization ) THEN
      
        grad(:,1)                = PES_gradient(:,1)
        grad(:,num_of_images)    = PES_gradient(:,num_of_images)
      
        norm_grad(1)             = norm( grad(:,1) )
  norm_grad(num_of_images) = norm( grad(:,num_of_images) )

      END IF

      gradient_loop: DO i = 2, ( num_of_images - 1 )

        IF ( algorithm == 6 ) THEN

!! elastic gradient ( variable elastic consatnt is used )

          elastic_gradient = &
      ( k(i) * ( cubic_pbc( pos(:,i) - pos(:,(i-1)) ) ) - &
              k(i+1) * ( cubic_pbc( pos(:,(i+1)) - pos(:,i) ) ) ) 

        ELSE

!! elastic gradient only along the path ( variable elastic consatnt is used )  

          elastic_gradient = &
      ( k(i) * norm( cubic_pbc( pos(:,i) - pos(:,(i-1)) ) ) - &
              k(i+1) * norm( cubic_pbc( pos(:,(i+1)) - pos(:,i) ) ) ) * &
      tangent(:,i)
  
  END IF

!! total gradient on each replica ( climbing image is used if needed )
!! only the component of the PES gradient ortogonal to the path is taken into
!! account

        IF ( ( i == Emax_index ) .AND. ( climbing ) ) THEN

          grad(:,i) = PES_gradient(:,i) - 2.D0 * &
                      DOT_PRODUCT( PES_gradient(:,i) , tangent(:,i) ) * &
          tangent(:,i)

        ELSE

          grad(:,i) = PES_gradient(:,i) + elastic_gradient - &
                      DOT_PRODUCT( PES_gradient(:,i) , tangent(:,i) ) * &
          tangent(:,i)

        END IF
  
        norm_grad(i) = norm( grad(:,i) )

      END DO gradient_loop
      
    END SUBROUTINE gradient


    SUBROUTINE compute_error( err  )

      IMPLICIT NONE

      REAL (KIND=dbl), INTENT(OUT)  :: err   
      INTEGER                       :: i


      error(1) = norm_grad(1)

      err = error(1)

      DO i = 2, ( num_of_images - 1 ) 

!! the error is given by the norm of the component of the 
!! total energy gradient ortogonal to the path

        error(i) = norm( PES_gradient(:,i) - &
                         DOT_PRODUCT( PES_gradient(:,i) , tangent(:,i) ) * &
             tangent(:,i) )

        IF ( error(i) > err ) err = error(i)
        
      END DO
      
      error(num_of_images) = norm_grad( num_of_images )

      IF ( error(num_of_images) > err ) err = error(num_of_images)

    END SUBROUTINE compute_error


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


    FUNCTION cubic_pbc( vect )

      IMPLICIT NONE    
      
      REAL (KIND=dbl), DIMENSION(:), INTENT(IN)  :: vect
      REAL (KIND=dbl), DIMENSION( SIZE( vect ) ) :: cubic_pbc
      REAL (KIND=dbl)                            :: vect_i
      INTEGER                                    :: i    
    

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


    FUNCTION path_tangent( index )

      IMPLICIT NONE

      INTEGER, INTENT(IN)             :: index
      REAL (KIND=dbl), DIMENSION(dim) :: path_tangent
      REAL (KIND=dbl)                 :: V_previous, &
                                         V_actual  , &
                                         V_next
      REAL (KIND=dbl)                 :: abs_next, abs_previous
      REAL (KIND=dbl)                 :: delta_V_max, delta_V_min


      V_previous = V( index - 1 )
      V_actual   = V( index )
      V_next     = V( index + 1 )

      IF ( ( V_next > V_actual ) .AND. ( V_actual > V_previous ) ) THEN

        path_tangent = cubic_pbc( pos(:,( index + 1 )) - pos(:,index) )

      ELSE IF ( ( V_next < V_actual ) .AND. ( V_actual < V_previous ) ) THEN

        path_tangent = cubic_pbc( pos(:,index) - pos(:,( index - 1 )) )

      ELSE

        abs_next     = ABS( V_next - V_actual ) 
        abs_previous = ABS( V_previous - V_actual ) 

        delta_V_max = MAX( abs_next , abs_previous ) 
        delta_V_min = MIN( abs_next , abs_previous )

        IF ( V_next > V_previous ) THEN

          path_tangent = &
            cubic_pbc( pos(:,( index + 1 )) - pos(:,index) ) * delta_V_max + & 
            cubic_pbc( pos(:,index) - pos(:,( index - 1 )) ) * delta_V_min

        ELSE IF ( V_next < V_previous ) THEN

          path_tangent = &
            cubic_pbc( pos(:,( index + 1 )) - pos(:,index) ) * delta_V_min + & 
            cubic_pbc( pos(:,index) - pos(:,( index - 1 )) ) * delta_V_max

        ELSE
  
          path_tangent = &
            cubic_pbc( pos(:,( index + 1 )) - pos(:,( index - 1 )) ) 

        END IF

      END IF 

    END FUNCTION path_tangent


    SUBROUTINE write_restart

      IMPLICIT NONE

      INTEGER             :: i, j
      CHARACTER (LEN=95)  :: vel_file
      INTEGER, PARAMETER  :: unit = 10


      OPEN( UNIT = unit, FILE = restart_file, STATUS = "UNKNOWN", &
            ACTION = "WRITE" )
    
        DO i = 1, num_of_images

          WRITE(unit,*) "Replica: ", i
          WRITE(unit,fmt3) V(i)

          DO j = 1, dim, 3 
       
            WRITE(unit,fmt1) pos(j,i),     & 
                             pos((j+1),i), &
                             pos((j+2),i), &
                             fix_atom(j),     &
                             fix_atom((j+1)), &
                             fix_atom((j+2)), &
                             -PES_gradient(j,i),     &
                             -PES_gradient((j+1),i), &
                             -PES_gradient((j+2),i)
  
          END DO

        END DO

      CLOSE( UNIT = unit )
      
      vel_file = TRIM( scratch_dir )//"/velocities_file"
      
!!DEBUG      PRINT *, "writing ", vel_file
      
      IF ( algorithm >= 4 ) THEN
      
        OPEN( UNIT = unit, FILE = vel_file, STATUS = "UNKNOWN", &
        ACTION = "WRITE" )
    
          DO i = 1, num_of_images

            WRITE(unit,*) "Replica: ", i

            DO j = 1, dim, 3 
       
              WRITE(unit,fmt2) vel(j,i),     & 
                               vel((j+1),i), &
                               vel((j+2),i)
  
            END DO

          END DO

        CLOSE( UNIT = unit )
      
      END IF

    END SUBROUTINE write_restart


    SUBROUTINE write_dat_files

      IMPLICIT NONE

      INTEGER                                    :: i, j
      REAL (KIND=dbl)                            :: R, delta_R, x
      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: d_R
      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: a, b, c, d, F
      REAL (KIND=dbl), DIMENSION(:), ALLOCATABLE :: reaction_coordinate
      REAL (KIND=dbl)                            :: E, E_0
      INTEGER, PARAMETER                         :: max_i = 1000 
      INTEGER, PARAMETER                         :: unit = 10


      ALLOCATE( d_R( dim ) )

      ALLOCATE( a( num_of_images - 1 ) )
      ALLOCATE( b( num_of_images - 1 ) )
      ALLOCATE( c( num_of_images - 1 ) )
      ALLOCATE( d( num_of_images - 1 ) )
      ALLOCATE( F( num_of_images ) )
      ALLOCATE( reaction_coordinate( num_of_images ) )

      F = 0.D0

      DO i = 2, ( num_of_images - 1 )

        F(i) = DOT_PRODUCT( - PES_gradient(:,i) , tangent(:,i) )

      END DO

      reaction_coordinate(1) = 0.D0

      DO i = 1, ( num_of_images - 1 )

        d_R = cubic_pbc( pos(:,( i + 1 )) - pos(:,i) ) 

        R = norm( d_R )

        reaction_coordinate(i+1) = reaction_coordinate(i) + R

        a(i) = 2.D0 * ( V(i) - V(i+1) ) / R**(3) - &
               ( F(i) + F(i+1) ) / R**(2)
            
        b(i) = 3.D0 * ( V(i+1) - V(i) ) / R**(2) + &
               ( 2.D0 * F(i) + F(i+1) ) / R

        c(i) = - F(i)

        d(i) = V(i)

      END DO
      
      OPEN( UNIT = unit, FILE = data_file, STATUS = "UNKNOWN", &
            ACTION = "WRITE" )

        WRITE(unit,'(3(2X,F12.8))') 0.D0, 0.D0, 0.D0

        DO i = 2, num_of_images
  
          WRITE(unit,'(3(2X,F12.8))') reaction_coordinate(i) * a_zero, &
                                      ( V(i) - V(1) ) * E_zero,        &
                                      error(i) * ( E_zero / a_zero )
      
        END DO

      CLOSE( UNIT = unit )  

      OPEN( UNIT = unit, FILE = interpolation_file, STATUS = "UNKNOWN", &
            ACTION = "WRITE" )

        i = 1

        delta_R = reaction_coordinate(num_of_images) / DBLE(max_i)

        DO j = 0, max_i

          R = DBLE(j) * delta_R 
 
          IF ( ( R > reaction_coordinate(i+1) ) .AND. &
               ( i < ( num_of_images - 1 ) ) ) i = i + 1

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
      DEALLOCATE( F )
      DEALLOCATE( reaction_coordinate )

    END SUBROUTINE write_dat_files


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
      IF ( ALLOCATED( delta_pos ) )        DEALLOCATE( delta_pos )
      IF ( ALLOCATED( pos ) )              DEALLOCATE( pos )
      IF ( ALLOCATED( grad ) )             DEALLOCATE( grad )
      IF ( ALLOCATED( old_grad ) )         DEALLOCATE( old_grad )
      IF ( ALLOCATED( PES_gradient ) )     DEALLOCATE( PES_gradient )
      IF ( ALLOCATED( fix_atom ) )         DEALLOCATE( fix_atom )
      IF ( ALLOCATED( V ) )                DEALLOCATE( V )
      IF ( ALLOCATED( k ) )                DEALLOCATE( k )
      IF ( ALLOCATED( norm_grad ) )        DEALLOCATE( norm_grad )
      IF ( ALLOCATED( error ) )            DEALLOCATE( error )
      
      IF ( ALLOCATED( tangent ) )          DEALLOCATE( tangent ) 
      IF ( ALLOCATED( elastic_gradient ) ) DEALLOCATE( elastic_gradient )

      IF ( ALLOCATED( conj_dir_i ) )       DEALLOCATE( conj_dir_i )

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
