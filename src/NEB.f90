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
  use module_types

  IMPLICIT NONE
 
  CHARACTER (LEN=80)                     :: first_config, last_config
  CHARACTER (LEN=80)                     :: scratch_dir
  CHARACTER (LEN=80)                     :: data_file, interpolation_file, barrier_file
  CHARACTER (LEN=80)                     :: restart_file
  CHARACTER (LEN=80)                     :: job_name
  LOGICAL                                :: restart, external_call
  INTEGER                                :: N
  REAL (gp)                              :: Lx, Ly, Lz
  REAL (gp), DIMENSION(:,:), ALLOCATABLE :: pos
  REAL (gp), DIMENSION(:,:), ALLOCATABLE :: PES_gradient
  real (gp), DIMENSION(:), ALLOCATABLE   :: fix_atom
  REAL (gp), DIMENSION(:), ALLOCATABLE   :: V, F, error
  REAL (gp)                              :: tolerance

  integer, dimension(4) :: mpi_info
  character(len=60) :: run_id
  type(input_variables) :: ins
  type(atoms_data) :: atoms
  type(restart_objects) :: rst

END MODULE NEB_variables

!> Module for NEB calculations
MODULE NEB_routines
  use module_defs
  use module_types
  use module_images
  use module_interfaces
  USE NEB_variables
  
  IMPLICIT NONE

  CHARACTER (LEN=*), PARAMETER ::                                              &
  fmt1 = "(3(2X,F12.8),3(2X,I1),3(2X,F12.8))",                                 &
  fmt2 = "(3(2X,F12.8))",                                                      &
  fmt3 = "(2X,F16.8)",                                                         &
  fmt4 = "(' iteration: ',I3,5X,'E activation =',F10.6,5X,'error =',F10.6)",   &
  fmt5 = "(' image: ',I2,'   Energy=  ',F16.8,'   Error=',F8.5)"

  type(NEB_data), private :: neb_
  type(mpi_environment), private :: neb_mpi

  CONTAINS

    SUBROUTINE read_input

      use module_interfaces

      IMPLICIT NONE

      INTEGER                                    :: i, j, n1, n2, ios, j0, k, num_of_images
      real(gp)                           :: r, a
      CHARACTER (LEN=20)                         :: minimization_scheme
      CHARACTER (LEN=95)                         :: vel_file
      INTEGER, PARAMETER                         :: unit = 10
      REAL (gp), DIMENSION(:), ALLOCATABLE :: d_R
      REAL (gp), DIMENSION(:,:), ALLOCATABLE :: vel0
      real(gp), dimension(:,:), pointer  :: rxyz
      real(gp), dimension(3, 1000)       :: xcart1, xcart2
      real(gp), dimension(3)             :: acell1, acell2
      logical :: file_exists, climbing, optimization
      integer :: max_iterations, ierr, nconfig
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


      call bigdft_init(mpi_info, nconfig, run_id, ierr)
      neb_mpi = mpi_environment_null()
      neb_mpi%igroup = 0
      neb_mpi%ngroup = 1
      neb_mpi%iproc  = mpi_info(3)
      neb_mpi%nproc  = mpi_info(4)
      neb_mpi%mpi_comm = MPI_COMM_NULL
      if (neb_mpi%nproc > 1 .and. mpi_info(1) == 0) then
         call create_rank_comm(bigdft_mpi%mpi_comm, neb_mpi%mpi_comm)
      end if

!! default values are assigned
      external_call     = (mpi_info(4) == 1) .and. (mpi_info(2) == 1)

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

      open(unit = 123, file = trim(run_id), action = "read")
      READ(123 , NML=NEB )
      close(123)

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

      call atoms_nullify(atoms)
      call read_atomic_file(trim(first_config), mpi_info(1), atoms, rxyz)
      n1 = atoms%nat
      acell1(1) = atoms%alat1
      acell1(2) = atoms%alat2
      acell1(3) = atoms%alat3
      xcart1(:,1:atoms%nat) = rxyz
      if (acell1(1) == 0.) acell1(1) = maxval(rxyz(1,:)) - minval(rxyz(1,:))
      if (acell1(2) == 0.) acell1(2) = maxval(rxyz(2,:)) - minval(rxyz(2,:))
      if (acell1(3) == 0.) acell1(3) = maxval(rxyz(3,:)) - minval(rxyz(3,:))
      deallocate(rxyz)
      call deallocate_atoms(atoms, "read_input")

      call atoms_nullify(atoms)
      call read_atomic_file(trim(last_config), mpi_info(1), atoms, rxyz)
      n2 = atoms%nat
      acell2(1) = atoms%alat1
      acell2(2) = atoms%alat2
      acell2(3) = atoms%alat3
      xcart2(:,1:atoms%nat) = rxyz
      if (acell2(1) == 0.) acell2(1) = maxval(rxyz(1,:)) - minval(rxyz(1,:))
      if (acell2(2) == 0.) acell2(2) = maxval(rxyz(2,:)) - minval(rxyz(2,:))
      if (acell2(3) == 0.) acell2(3) = maxval(rxyz(3,:)) - minval(rxyz(3,:))
      deallocate(rxyz)

      if (atoms%geocode == 'F') then
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

     if (.not. external_call) then
        call standard_inputfile_names(ins,trim(run_id), mpi_info(1))
        call default_input_variables(ins)
        call inputs_parse_params(ins, mpi_info(1), .true.)
        call inputs_parse_add(ins, atoms, mpi_info(1), .true.)

        call init_atomic_values((mpi_info(1) == 0), atoms, ins%ixc)
        call read_atomic_variables(atoms, trim(ins%file_igpop),ins%nspin)

        call restart_objects_new(rst)
        call restart_objects_set_mode(rst, ins%inputpsiid)
        call restart_objects_set_nat(rst, atoms%nat, "read_input")
        call restart_objects_set_mat_acc(rst, mpi_info(1), ins%matacc)
     end if

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
      use yaml_output

      IMPLICIT NONE

      INTEGER :: iteration, unt, ierr
      LOGICAL :: stat

      open(unit = 456, file = trim(barrier_file), action = "WRITE")
      write(456, "(A)") "# NEB barrier file"
      close(unit = 456)

      IF ( .NOT. restart) THEN
         CALL write_restart(restart_file, neb_%ndim, neb_%nimages, V, pos, fix_atom, PES_gradient)
      END IF

      error(:) = 999d99

      if (mpi_info(1) == 0 .and. mpi_info(3) == 0) &
           & call yaml_open_sequence("NEB minimization loop", unit = 6)
      iteration = 0
      minimization: do
         if (external_call) then
            CALL PES_IO(neb_%optimization .or. (.not. restart .and. iteration == 0),stat)
            if (.not. stat) exit minimization
         else
            call PES_internal(neb_%optimization .or. (.not. restart .and. iteration == 0), iteration)
         end if

         call compute_neb_pos(stat, iteration, error, F, pos, V, PES_gradient, Lx, Ly, Lz, neb_)

         if (mpi_info(1) == 0 .and. mpi_info(3) == 0) then
            CALL write_restart(restart_file, neb_%ndim, neb_%nimages, V, pos, fix_atom, PES_gradient)
            if (neb_%algorithm >= 4) then
               CALL write_restart_vel(neb_, trim(scratch_dir) // "velocities_file")
            end if

            if (iteration > 1) then
               CALL write_dat_files(data_file, interpolation_file, neb_%ndim, neb_%nimages, pos, V, F, error, Lx, Ly, Lz)
               open(unit = 456, file = trim(barrier_file), action = "WRITE", position = "APPEND")
               WRITE(456, fmt4) iteration - 1, ( maxval(V(2:neb_%nimages - 1)) - V(1) ) * Ha_eV, &
                    & maxval(error) * ( Ha_eV / Bohr_Ang ) 
               close(unit = 456)
               call yaml_swap_stream(6, unt, ierr)
               call yaml_sequence(advance='no')
               call images_write_step(V, error, neb_%nimages, iteration = iteration - 1)
               call yaml_set_default_stream(unt, ierr)
            end if
         end if

         if (.not. stat) exit minimization
      end do minimization
      if (mpi_info(1) == 0 .and. mpi_info(3) == 0) &
           & call yaml_close_sequence(unit = 6)

      if (mpi_info(1) == 0 .and. mpi_info(3) == 0) then
         call yaml_swap_stream(6, unt, ierr)
         call yaml_comment('Final results',hfill='-')
         call images_write_step(V, error, neb_%nimages, full = .true.)
         call yaml_set_default_stream(unt, ierr)
      end if
    END SUBROUTINE search_MEP

    subroutine PES_internal( flag, iteration )
      use yaml_output
      implicit none
      logical, intent(in) :: flag
      integer, intent(in) :: iteration

      logical, dimension(neb_%nimages) :: update
      integer, dimension(neb_%nimages) :: igroup
      integer :: i, infocode, ierr
      type(run_objects) :: runObj
      type(DFT_global_output) :: outs
      character(len = 80) :: file

      ! update() is a mask of images to compute.
      update = .true.
      where ( error * Ha_eV / Bohr_Ang <= neb_%convergence ) update = .false.
      if (.not. flag) then
         update(1) = .false.
         update(neb_%nimages) = .false.
      end if

      ! Put to zero V and PES_gradient for the updated images,
      ! or for all groups not 0 (for later reduction).
      do i = 1, neb_%nimages
         if (update(i) .or. neb_mpi%iproc > 0) then
            V(i) = 0.d0
            call to_zero(neb_%ndim, PES_gradient(1,i))
         end if
      end do
      
      ! Do the calculations, distributing along taskgroups.
      call images_distribute_tasks(igroup, update, neb_%nimages, neb_mpi%nproc)
      do i = 1, neb_%nimages
         if (igroup(i) - 1 == mpi_info(3)) then
            write(file, "(A,I2.2,A)") "log-img", i, ".yaml"
            call yaml_set_stream(unit = 9169 + i, filename = file, istat = ierr)
            call yaml_comment("NEB iteration #" // trim(yaml_toa(iteration, fmt = "(I3.3)")), hfill="-")
            call run_objects_init_container(runObj, ins, atoms, rst, pos(1,i))
            call init_global_output(outs, runObj%atoms%nat)
            call call_bigdft(runObj,outs,mpi_info(2),mpi_info(1),infocode)
            call run_objects_free_container(runObj)
            ins%inputpsiid = 1
            call dsbmv('L', neb_%ndim, 0, -1.d0, fix_atom(1), 1, outs%fxyz(1,1), 1, &
                 & 0.d0, PES_gradient(1, i), 1)
            V(i) = outs%energy
            call deallocate_global_output(outs)
            close(unit = 9169 + i)
         end if
      end do
      
      ! Reduce the result in case of taskgrouping.
      if (neb_mpi%nproc > 1) then
         if (neb_mpi%mpi_comm /= MPI_COMM_NULL) then
            call mpiallred(V(1), neb_%nimages, MPI_SUM, neb_mpi%mpi_comm, ierr)
            call mpiallred(PES_gradient(1,1), neb_%ndim * neb_%nimages, &
                 & MPI_SUM, neb_mpi%mpi_comm, ierr)
         end if
         call mpi_bcast(V(1), neb_%nimages, MPI_DOUBLE_PRECISION, 0, bigdft_mpi%mpi_comm, ierr)
         call mpi_bcast(PES_gradient(1,1), neb_%ndim * neb_%nimages, MPI_DOUBLE_PRECISION, &
              & 0, bigdft_mpi%mpi_comm, ierr)
      end if
    END SUBROUTINE PES_internal

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
                                    fix_atom * (-1.d0)

        END DO

      CLOSE( UNIT = unit )   

    END SUBROUTINE PES_IO

    SUBROUTINE deallocation
      use yaml_output
      use dynamic_memory
      IMPLICIT NONE

      integer :: ierr

      IF ( ALLOCATED( pos ) )              DEALLOCATE( pos )
      IF ( ALLOCATED( PES_gradient ) )     DEALLOCATE( PES_gradient )
      IF ( ALLOCATED( fix_atom ) )         DEALLOCATE( fix_atom )
      IF ( ALLOCATED( F ) )                DEALLOCATE( F )
      IF ( ALLOCATED( V ) )                DEALLOCATE( V )
      IF ( ALLOCATED( error ) )            DEALLOCATE( error )
      
      call deallocate_images(neb_)
      call deallocate_atoms(atoms, "deallocation")

      if (.not. external_call) then
         call free_input_variables(ins)
         call free_restart_objects(rst,"deallocation")
         call yaml_set_stream(unit = 6, istat = ierr)
         call f_finalize()
         call yaml_close_all_streams()
      end if

      call bigdft_finalize(ierr)
    END SUBROUTINE deallocation

END MODULE NEB_routines

PROGRAM NEB

  USE NEB_routines

  IMPLICIT NONE

  CALL read_input

  CALL search_MEP

  CALL deallocation

  STOP

END PROGRAM NEB
