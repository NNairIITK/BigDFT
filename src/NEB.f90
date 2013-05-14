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


!> Module for NEB calculations (variables)
MODULE NEB_variables
  use module_defs

  IMPLICIT NONE
 
  CHARACTER (LEN=80)                     :: first_config, last_config
  CHARACTER (LEN=80)                     :: scratch_dir
  CHARACTER (LEN=80)                     :: data_file, interpolation_file, barrier_file
  CHARACTER (LEN=80)                     :: restart_file
  CHARACTER (LEN=80)                     :: job_name
  LOGICAL                                :: restart
  INTEGER                                :: N
  REAL (gp)                              :: Lx, Ly, Lz
  REAL (gp), DIMENSION(:,:), ALLOCATABLE :: pos
  REAL (gp), DIMENSION(:,:), ALLOCATABLE :: PES_gradient
  INTEGER, DIMENSION(:), ALLOCATABLE     :: fix_atom
  REAL (gp), DIMENSION(:), ALLOCATABLE   :: V, F, error
  REAL (gp)                              :: tolerance

END MODULE NEB_variables

!> Module for NEB calculations
MODULE NEB_routines
  use module_defs
  use module_images
  USE NEB_variables
  
  IMPLICIT NONE

  CHARACTER (LEN=*), PARAMETER ::                                              &
  fmt1 = "(3(2X,F12.8),3(2X,I1),3(2X,F12.8))",                                 &
  fmt2 = "(3(2X,F12.8))",                                                      &
  fmt3 = "(2X,F16.8)",                                                         &
  fmt4 = "(' iteration: ',I3,5X,'E activation =',F10.6,5X,'error =',F10.6)",   &
  fmt5 = "(' image: ',I2,'   Energy=  ',F16.8,'   Error=',F8.5)"

  type(NEB_data), private :: neb_

  CONTAINS

    SUBROUTINE read_input

      use module_types
      use module_interfaces

      IMPLICIT NONE

      INTEGER                                    :: i, j, n1, n2, ios, j0, k, num_of_images
      real(gp)                           :: r, a
      CHARACTER (LEN=20)                         :: minimization_scheme
      CHARACTER (LEN=95)                         :: vel_file
      INTEGER, PARAMETER                         :: unit = 10
      REAL (gp), DIMENSION(:), ALLOCATABLE :: d_R
      REAL (gp), DIMENSION(:,:), ALLOCATABLE :: vel0
      type(atoms_data)                           :: at
      real(gp), dimension(:,:), pointer  :: rxyz
      real(gp), dimension(3, 1000)       :: xcart1, xcart2
      real(gp), dimension(3)             :: acell1, acell2
      logical :: file_exists, climbing, optimization
      integer :: max_iterations
      real(gp) :: convergence, damp, k_min, k_max, ds, temp_req

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

      scratch_dir       = "./"

      job_name          = "neb"
      restart           = .FALSE.
      neb_%climbing     = .FALSE.
      neb_%optimization = .FALSE.
      
      minimization_scheme = "quick-min"
      neb_%damp           = 1.D0
      neb_%temp_req       = 0.D0
      
      neb_%k_max = 0.1D0
      neb_%k_min = 0.1D0
      
      neb_%ds = 0.5D0
      
      neb_%max_iterations = 1
      
      tolerance   = 1.0D-4
      neb_%convergence = 5.0D-2

      N = 0
      
      Lx = 0.D0
      Ly = 0.D0
      Lz = 0.D0

      READ( * , NML=NEB )
      neb_%climbing = climbing
      neb_%optimization = optimization
      neb_%convergence = convergence
      neb_%damp = damp
      neb_%k_min = k_min
      neb_%k_max = k_max
      neb_%ds = ds
      neb_%temp_req = temp_req
      neb_%max_iterations = max_iterations
      neb_%nimages = num_of_images

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
!!!         neb_%nimages = neb_%nimages - 2
!!!      END IF

      IF ( minimization_scheme == "steepest_descent" ) THEN

         neb_%algorithm = 1

      ELSE IF ( minimization_scheme == "fletcher-reeves" ) THEN

         neb_%algorithm = 2

      ELSE IF ( minimization_scheme == "polak-ribiere" ) THEN

         neb_%algorithm = 3

      ELSE IF ( minimization_scheme == "quick-min" ) THEN

         neb_%algorithm = 4

      ELSE IF ( minimization_scheme == "damped-verlet" ) THEN

         neb_%algorithm = 5

      ELSE IF ( minimization_scheme == "sim-annealing" ) THEN

         neb_%algorithm = 6

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

      IF ( neb_%nimages <= 2 ) THEN

        WRITE(*,'(T1,"read_input: neb_%nimages must be larger than 2")')
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

      neb_%ndim = 3 * N
      
      CALL dyn_allocation

!! all the arrays are initialized

      V                = 0.D0
      PES_gradient     = 0.D0
      error            = 0.D0

      IF ( restart ) THEN

        OPEN( UNIT = unit, FILE = restart_file, STATUS = "OLD", &
              ACTION = "READ" )
    
        ! Read as many configurations as contained in the
        ! Restart file.
        DO i = 1, neb_%nimages

           READ(unit,*, iostat = ios)
           if (ios /= 0) exit
           READ(unit,fmt3) V(i)

           DO j = 1, neb_%ndim, 3 

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

        END DO

        CLOSE( UNIT = unit) 

        ! Check that we read at least two configurations...
        i = i - 1
        if (i < 2) then
           WRITE(*,'(T2,"read_input: number of replica in restart file is less than 2")')
           WRITE(*,'(T2,"            N = ", I8 )') i - 1
           STOP  
        end if

        ! Add some intermediate configurations if neb_%nimages > i
        if (neb_%nimages > i) then
           ALLOCATE( d_R(neb_%ndim) )           
           r = real(i - 1, gp) / real(neb_%nimages - 1, gp)
           a = real(0, gp)
           j0 = neb_%nimages
           do j = neb_%nimages, 1, -1
              ! j : the current position in pos array.
              ! i : the current position in read configurations.
              ! r : the ideal ratio (i_0 - 1) / (neb_%nimages - 1)
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
              a = (r * real(neb_%nimages - 1, gp) - real(i, gp) + real(1, gp)) / &
                   & (real(neb_%nimages, gp) - real(j, gp) + real(1, gp))
           end do
           deallocate(d_R)

           CALL write_restart(restart_file, neb_%ndim, neb_%nimages, V, pos, fix_atom, PES_gradient)

        end if

        vel_file = TRIM( scratch_dir )//"/velocities_file"
        inquire(FILE = vel_file, EXIST = file_exists)
        IF ( ( neb_%algorithm >= 4 ) .AND. file_exists ) THEN
!!DEBUG          PRINT *, "reading ", vel_file
           allocate(vel0(neb_%ndim, neb_%nimages))
           OPEN( UNIT = unit, FILE = vel_file, STATUS = "OLD", ACTION = "READ" )
           DO i = 1, neb_%nimages
              READ(unit,*)
              DO j = 1, neb_%ndim, 3 
                 READ(unit,fmt2) vel0(j,i),     & 
                      vel0((j+1),i), &
                      vel0((j+2),i)
              END DO
           END DO
           CLOSE( UNIT = unit )
           call set_init_vel(neb_, vel0)
           deallocate(vel0)
        END IF
      
      ELSE

        ALLOCATE( d_R(neb_%ndim) )           

        call vcopy(neb_%ndim, xcart1(1,1), 1, pos(1,1), 1)
        call vcopy(neb_%ndim, xcart2(1,1), 1, pos(1,neb_%nimages), 1)
  
        d_R = ( pos(:,neb_%nimages) - pos(:,1) ) / &
             DBLE( neb_%nimages - 1 )

        fix_atom = 1
        WHERE ( ABS( d_R ) <=  tolerance ) fix_atom = 0

        DO i = 2, ( neb_%nimages - 1 )
           pos(:,i) = pos(:,( i - 1 )) + d_R(:)
        END DO

        DEALLOCATE( d_R )

     END IF

   END SUBROUTINE read_input

    
    SUBROUTINE dyn_allocation

      IMPLICIT NONE
     
      ALLOCATE( pos( neb_%ndim , neb_%nimages ) )

      ALLOCATE( PES_gradient( neb_%ndim , neb_%nimages ) )          

      ALLOCATE( fix_atom( neb_%ndim ) )

      ALLOCATE( F( neb_%nimages ) )
      ALLOCATE( V( neb_%nimages ) )
      ALLOCATE( error( neb_%nimages ) )

      call init_images(neb_, neb_%ndim, neb_%nimages, neb_%algorithm)

    END SUBROUTINE dyn_allocation     


    SUBROUTINE search_MEP

      IMPLICIT NONE

      INTEGER         :: iteration
      LOGICAL         :: stat


      open(unit = 456, file = trim(barrier_file), action = "WRITE")
      write(456, "(A)") "# NEB barrier file"
      close(unit = 456)

      IF ( .NOT. restart) THEN
         CALL write_restart(restart_file, neb_%ndim, neb_%nimages, V, pos, fix_atom, PES_gradient)
      END IF

      iteration = 0
      minimization: do
         CALL PES_IO(neb_%optimization .or. (.not. restart .and. iteration == 0),stat)
         if (.not. stat) exit minimization

         call compute_neb_pos(stat, iteration, error, F, pos, V, PES_gradient, Lx, Ly, Lz, neb_)

         CALL write_restart(restart_file, neb_%ndim, neb_%nimages, V, pos, fix_atom, PES_gradient)
         if (neb_%algorithm >= 4) then
            CALL write_restart_vel(neb_, trim(scratch_dir) // "velocities_file")
         end if

         if (iteration > 1) then
            CALL write_dat_files(data_file, interpolation_file, neb_%ndim, neb_%nimages, pos, V, F, error, Lx, Ly, Lz)
            open(unit = 456, file = trim(barrier_file), action = "WRITE", position = "APPEND")
            WRITE(456, fmt4) iteration - 1, ( maxval(V(2:neb_%nimages - 1)) - V(1) ) * Ha_eV, &
                 & maxval(error) * ( Ha_eV / Bohr_Ang ) 
            WRITE(*, fmt4)   iteration - 1, ( maxval(V(2:neb_%nimages - 1)) - V(1) ) * Ha_eV, &
                 & maxval(error) * ( Ha_eV / Bohr_Ang ) 
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
      REAL (gp)            :: temp_V
      REAL (gp), PARAMETER :: corruption_flag = 9999999.99999999
      INTEGER, PARAMETER         :: unit = 10     
      REAL (gp), PARAMETER :: epsi = 1.0D-8

      IF ( flag ) THEN

         CALL SYSTEM( "./NEB_driver.sh all " // trim(job_name) // &
              & " " // trim(scratch_dir) // " " // trim(first_config))

        N_in  = 1
        N_fin = neb_%nimages

      ELSE
         
         CALL SYSTEM( "./NEB_driver.sh free_only " // trim(job_name) // &
              & " " // trim(scratch_dir) // " " // trim(first_config))

        N_in  = 2
        N_fin = ( neb_%nimages - 1 )

      END IF

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

          DO i = 1, neb_%ndim, 3 

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


      DO i = 1, neb_%nimages

        WRITE(*,fmt5) i, V(i) * Ha_eV, error(i) * ( Ha_eV / Bohr_Ang )

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
      
      call deallocate_images(neb_)
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
