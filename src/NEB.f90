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
 
  CHARACTER (LEN=80)                     :: scratch_dir
  CHARACTER (LEN=80)                     :: data_file, interpolation_file
  CHARACTER (LEN=80)                     :: restart_file
  CHARACTER (LEN=80)                     :: job_name
  LOGICAL                                :: external_call
  real (gp), DIMENSION(:,:), ALLOCATABLE :: fix_atom

  integer, dimension(4) :: mpi_info
  character(len=60), dimension(:), allocatable :: arr_posinp,arr_radical
  type(input_variables), dimension(:), allocatable :: ins
  type(atoms_data), dimension(:), allocatable :: atoms
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
  type(run_image), dimension(:), allocatable :: imgs
  type(mpi_environment), private :: neb_mpi

  CONTAINS

    SUBROUTINE read_input
      use yaml_output
      use dictionaries
      use module_interfaces
      use module_input_keys, only: input_keys_fill_all, PERF_VARIABLES, OUTDIR, &
           & GEOPT_VARIABLES
      use module_input_dicts
      use module_atoms, only: atomic_structure, &
           deallocate_atomic_structure, &
           read_atomic_file => set_astruct_from_file, &
           astruct_nullify => nullify_atomic_structure

      IMPLICIT NONE

      INTEGER :: i, j, num_of_images, istart, istop
      CHARACTER (LEN=max_field_length) :: minimization_scheme
      INTEGER, PARAMETER :: unit = 10
      REAL (gp), DIMENSION(:,:), ALLOCATABLE :: d_R
      real(gp), dimension(3) :: acell1, acell2
      integer :: ierr, nconfig, algorithm, unit_log
      type(mpi_environment) :: bigdft_mpi_svg
      character(len=60) :: run_id
      type(dictionary), pointer :: dict, dict_min
      REAL (gp) :: tolerance

      call f_lib_initialize()
      nullify(dict)
      call bigdft_init(mpi_info, nconfig, run_id, ierr)
      neb_mpi = mpi_environment_null()
      neb_mpi%igroup = mpi_info(1)
      neb_mpi%ngroup = mpi_info(2)
      neb_mpi%iproc  = mpi_info(3)
      neb_mpi%nproc  = mpi_info(4)
      neb_mpi%mpi_comm = MPI_COMM_NULL
      if (neb_mpi%nproc > 1) then
         call create_rank_comm(bigdft_mpi%mpi_comm, neb_mpi%mpi_comm)
      end if

!! default values are assigned
      external_call     = (mpi_info(4) == 1) .and. (mpi_info(2) == 1)

      call dict_init(dict)
      call read_input_dict_from_files(trim(run_id), bigdft_mpi, dict)
      !if (mpi_info(1) == 0 .and. mpi_info(3) == 0) call yaml_dict_dump(dict)
      call input_keys_fill_all(dict, dict_min)
      call neb_set_from_dict(dict, neb_%optimization, neb_%climbing, &
           & neb_%max_iterations, num_of_images, neb_%convergence, tolerance, &
           & neb_%ds, neb_%k_min, neb_%k_max, neb_%temp_req, neb_%damp, &
           & minimization_scheme)
      ! NEB is using cv criterion in ev per ang.
      neb_%convergence = neb_%convergence * Ha_eV / Bohr_Ang
      call dict_free(dict)
      call dict_free(dict_min)
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
              trim(minimization_scheme)
         WRITE(*,'(T2,"            does not exist")') 
         STOP 
      END IF

      scratch_dir       = "./"
      job_name          = "neb"
      if (trim(run_id) /= "input") write(job_name, "(A)") trim(run_id)

      allocate(arr_radical(abs(num_of_images)))
      allocate(arr_posinp(abs(num_of_images)))
      call bigdft_get_run_ids(num_of_images,trim(run_id),arr_radical,arr_posinp,ierr)

      allocate( ins(num_of_images), atoms(num_of_images) )
      allocate( imgs(num_of_images) )

      call dict_init(dict)
      ! trick to output the image logs where it should, on disk.
      call set(dict // PERF_VARIABLES // OUTDIR, "./")
      ! Trick here, only super master will read the input files...
      bigdft_mpi_svg = bigdft_mpi
      bigdft_mpi%mpi_comm = MPI_COMM_WORLD
      call mpi_comm_rank(MPI_COMM_WORLD, bigdft_mpi%iproc, ierr)
      call mpi_comm_size(MPI_COMM_WORLD, bigdft_mpi%nproc, ierr)
      bigdft_mpi%igroup = 0
      bigdft_mpi%ngroup = num_of_images
      do i = 1, num_of_images

         call user_dict_from_files(dict, trim(arr_radical(i)), &
              & trim(arr_posinp(i)), bigdft_mpi)
         ! Force no geometry relaxation
         call pop(dict, GEOPT_VARIABLES)
         call inputs_from_dict(ins(i), atoms(i), dict)

         if (.not. external_call .and. i == 1) then
            call restart_objects_new(rst)
            call restart_objects_set_mode(rst, ins(1)%inputpsiid)
            call restart_objects_set_nat(rst, atoms(1)%astruct%nat, "read_input")
            call restart_objects_set_mat_acc(rst, bigdft_mpi%iproc, ins(1)%matacc)
         end if

         ! Some consistency checks.
         IF ( atoms(1)%astruct%nat /= atoms(i)%astruct%nat ) THEN
            WRITE(*,'(T2,"read_input: number of atoms is not constant")')
            WRITE(*,'(T2,"            N = ", I8, I8 )') atoms(1)%astruct%nat, atoms(i)%astruct%nat
            STOP  
         END IF

         if (bigdft_mpi%iproc == 0) then
            ! Need to close streams here to avoid running out of available streams.
            call yaml_get_default_stream(unit_log)
            call yaml_close_stream(unit_log)
         end if

         call image_init(imgs(i), ins(i), atoms(i), rst, algorithm)
         ! Store forces if present (for restart).
         call global_output_set_from_dict(imgs(i)%outs, dict // "posinp")
      end do
      call dict_free(dict)

      data_file          = trim(job_name) // ".NEB.dat"
      interpolation_file = trim(job_name) // ".NEB.int"
      restart_file       = trim(job_name) // ".NEB.restart"

!! initial and final configuration are read only if a new simulation
!! is started ( restart = .FALSE. ) 
      
      acell1 = atoms(1)%astruct%cell_dim
      if (acell1(1) == 0.) acell1(1) = maxval(atoms(1)%astruct%rxyz(1,:)) - minval(atoms(1)%astruct%rxyz(1,:))
      if (acell1(2) == 0.) acell1(2) = maxval(atoms(1)%astruct%rxyz(2,:)) - minval(atoms(1)%astruct%rxyz(2,:))
      if (acell1(3) == 0.) acell1(3) = maxval(atoms(1)%astruct%rxyz(3,:)) - minval(atoms(1)%astruct%rxyz(3,:))

      acell2 = atoms(num_of_images)%astruct%cell_dim
      if (acell2(1) == 0.) acell2(1) = maxval(atoms(num_of_images)%astruct%rxyz(1,:)) - &
           & minval(atoms(num_of_images)%astruct%rxyz(1,:))
      if (acell2(2) == 0.) acell2(2) = maxval(atoms(num_of_images)%astruct%rxyz(2,:)) - &
           & minval(atoms(num_of_images)%astruct%rxyz(2,:))
      if (acell2(3) == 0.) acell2(3) = maxval(atoms(num_of_images)%astruct%rxyz(3,:)) - &
           & minval(atoms(num_of_images)%astruct%rxyz(3,:))

      if (atoms(1)%astruct%geocode == 'F') then
        acell1 = max(acell1, acell2)
        acell2 = acell1
      end if

!! some consistency checks are done
      IF ( maxval(abs(acell2 - acell1)) > 1.d-6 ) THEN
         WRITE(*,'(T2,"read_input: box size is not constant")')
         WRITE(*,'(T2,"           dLx = ", F10.6 )') acell1(1) - acell2(1)
         WRITE(*,'(T2,"           dLy = ", F10.6 )') acell1(2) - acell2(2)
         WRITE(*,'(T2,"           dLz = ", F10.6 )') acell1(3) - acell2(3)
         STOP  
      END IF

!!$      IF ( restart ) THEN
!!$        vel_file = TRIM( scratch_dir )//"/velocities_file"
!!$        inquire(FILE = vel_file, EXIST = file_exists)
!!$        IF ( ( neb_%algorithm >= 4 ) .AND. file_exists ) THEN
!!$!!DEBUG          PRINT *, "reading ", vel_file
!!$           allocate(vel0(ndim, neb_%nimages))
!!$           OPEN( UNIT = unit, FILE = vel_file, STATUS = "OLD", ACTION = "READ" )
!!$           DO i = 1, neb_%nimages
!!$              READ(unit,*)
!!$              DO j = 1, ndim, 3 
!!$                 READ(unit,fmt2) vel0(j,i),     & 
!!$                      vel0((j+1),i), &
!!$                      vel0((j+2),i)
!!$              END DO
!!$           END DO
!!$           CLOSE( UNIT = unit )
!!$           call set_init_vel(neb_, vel0)
!!$           deallocate(vel0)
!!$        END IF
!!$      
!!$      ELSE

      ALLOCATE( d_R(3, atoms(1)%astruct%nat) )           

      istart = 1
      ! We set the coordinates for all empty images.
      DO i = 2, num_of_images
         if (maxval(abs(atoms(i)%astruct%rxyz-atoms(i-1)%astruct%rxyz)) > 1d-6) then
            istop = i
            d_R = ( atoms(istop)%astruct%rxyz - atoms(istart)%astruct%rxyz ) / &
                 DBLE( istop - istart )
            do j = istart + 1, istop - 1, 1
               atoms(j)%astruct%rxyz = atoms(j - 1)%astruct%rxyz + d_R
               ! Dump generated image positions on disk.
               if (bigdft_mpi%iproc == 0) then
                  call write_atomic_file(trim(arr_posinp(j)) // ".in", UNINITIALIZED(1.d0), &
                       & atoms(j)%astruct%rxyz, atoms(j), "NEB generated")
               end if
               ! Erase forces.
               imgs(j)%outs%fxyz(:,:) = UNINITIALIZED(1.d0)
            end do
            istart = i
         end if         
      END DO
      
      d_R = ( atoms(num_of_images)%astruct%rxyz - atoms(1)%astruct%rxyz )
      ALLOCATE( fix_atom(3, atoms(1)%astruct%nat) )      
      fix_atom = 1
      WHERE ( ABS( d_R ) <=  tolerance ) fix_atom = 0

      DEALLOCATE( d_R )

      ! End of trick.
      bigdft_mpi = bigdft_mpi_svg

!!$     END IF
    END SUBROUTINE read_input

    
    SUBROUTINE search_MEP
      use yaml_output

      IMPLICIT NONE

      INTEGER :: iteration, unt, ierr, i
      real(gp) :: err
      LOGICAL :: stat, restart
      CHARACTER (LEN=4), PARAMETER :: exit_file = "EXIT"  
      character(len = 256) :: filename

!!$      IF ( .NOT. restart) THEN
!!$         CALL write_restart(restart_file, neb_%ndim, neb_%nimages, V, pos, fix_atom, PES_gradient)
!!$      END IF

      ! Test for restart
      restart = .false.
      if (external_call) then
         inquire(file = trim(restart_file), exist = restart)
         call yaml_map("NEB restart", restart, unit = 6)
      else if (mpi_info(1) == 0 .and. mpi_info(3) == 0) then
         call yaml_open_sequence("Restarting images", unit = 6)
         do i = 1, size(imgs), 1
            call yaml_sequence(trim(yaml_toa(all(imgs(i)%outs%fxyz /= UNINITIALIZED(1.d0)))), unit = 6, advance = "no")
            call yaml_comment(yaml_toa(i, fmt = "(I2.2)"), unit = 6)
            call yaml_newline(unit = 6)
         end do
         call yaml_close_sequence(unit = 6)
      end if

      if (mpi_info(1) == 0 .and. mpi_info(3) == 0) &
           & call yaml_open_sequence("NEB minimization loop", unit = 6)
      iteration = 0
      minimization: do
         if (external_call) then
            CALL PES_IO(imgs, (neb_%optimization .or. (.not. restart .and. iteration == 0)),stat)
            if (.not. stat) exit minimization
         else
            call PES_internal(imgs, neb_%optimization,  (iteration == 0), iteration)
         end if

         call compute_neb_pos(imgs, iteration, neb_)

         if (mpi_info(1) == 0 .and. mpi_info(3) == 0) then
!!$            CALL write_restart(restart_file, neb_%ndim, neb_%nimages, V, pos, fix_atom, PES_gradient)
            if (imgs(1)%algorithm >= 4) then
               CALL write_restart_vel(trim(scratch_dir) // "velocities_file", imgs)
            end if

            if (iteration > 0) then
               CALL write_dat_files(trim(job_name), imgs, iteration)
               call yaml_swap_stream(6, unt, ierr)
               call yaml_sequence(advance='no')
               call images_output_step(imgs, iteration = iteration, tol = neb_%convergence)
               call yaml_set_default_stream(unt, ierr)
            end if
         end if

         if (iteration > 0 .or. neb_%max_iterations == 1) then
            err = maxval(images_get_errors(imgs))

            IF ( ( err * Ha_eV / Bohr_Ang ) <= neb_%convergence .or. neb_%max_iterations == 1)  THEN
               exit minimization
            END IF
         end if

         iteration = iteration + 1

         IF ( iteration > neb_%max_iterations ) THEN
            exit minimization
         END IF

         inquire(FILE = exit_file, EXIST = stat)
         IF ( stat ) THEN
            call delete(trim(exit_file),len(trim(exit_file)), stat)

            WRITE(*,*) " WARNING :  soft exit required"
            WRITE(*,*) " STOPPING ...                 "

            exit minimization
         END IF
      end do minimization
      if (mpi_info(1) == 0 .and. mpi_info(3) == 0) &
           & call yaml_close_sequence(unit = 6)

      if (mpi_info(1) == 0 .and. mpi_info(3) == 0) then
         call yaml_swap_stream(6, unt, ierr)
         call yaml_comment('Final results',hfill='-')
         call images_output_step(imgs, full = .true.)
         call yaml_set_default_stream(unt, ierr)

         do i = 1, size(imgs), 1
            filename=trim('final_'//trim(arr_posinp(i)))
            call write_atomic_file(filename, imgs(i)%outs%energy,imgs(i)%run%atoms%astruct%rxyz, &
                 & imgs(i)%run%atoms,'FINAL CONFIGURATION',forces=imgs(i)%outs%fxyz)
         end do
      end if
    END SUBROUTINE search_MEP

    subroutine PES_internal( imgs, flag_optim, flag_restart, iteration )
      use yaml_output
      implicit none
      type(run_image), dimension(:), intent(inout) :: imgs
      integer, intent(in) :: iteration
      logical, intent(in) :: flag_optim, flag_restart

      integer :: i
      logical, dimension(size(imgs)) :: update
      integer, dimension(size(imgs)) :: igroup
      real(gp), dimension(size(imgs)) :: errors

      if (.not. flag_restart) then
         errors = images_get_errors(imgs) * Ha_eV / Bohr_Ang
         do i = 1, size(imgs)
            if (errors(i) > neb_%convergence .and. &
                 & (flag_optim .or. (i > 1 .and. i < size(imgs)))) &
                 & imgs(i)%outs%fxyz(1,1) = UNINITIALIZED(1.d0)
         end do
      end if

      ! update() is a mask of images to compute.
      update = .true.
      do i = 1, size(imgs)
         if (all(imgs(i)%outs%fxyz /= UNINITIALIZED(1.d0))) update(i) = .false.
      end do

      ! Do the calculations, distributing among taskgroups.
      call images_distribute_tasks(igroup, update, size(imgs), neb_mpi%nproc)
      do i = 1, size(imgs)
         if (igroup(i) - 1 == mpi_info(3)) then
            call image_calculate(imgs(i), iteration, i)
            imgs(i)%outs%fxyz = imgs(i)%outs%fxyz * fix_atom
        end if
      end do
      call images_collect_results(imgs, igroup, size(imgs), neb_mpi)
    END SUBROUTINE PES_internal

    SUBROUTINE PES_IO( imgs, flag , stat )

      IMPLICIT NONE

      type(run_image), dimension(:), intent(inout) :: imgs
      LOGICAL, INTENT(IN)        :: flag
      LOGICAL, INTENT(OUT)       :: stat
      INTEGER                    :: i, replica
      INTEGER                    :: N_in, N_fin
      REAL (gp)            :: temp_V
      REAL (gp), PARAMETER :: corruption_flag = 9999999.99999999
      INTEGER, PARAMETER         :: unit = 10     
      REAL (gp), PARAMETER :: epsi = 1.0D-8

      call write_restart(trim(restart_file), imgs, fix_atom)

      IF ( flag ) THEN

         CALL SYSTEM( "./NEB_driver.sh all " // trim(job_name) // &
              & " " // trim(scratch_dir) // " " // trim(arr_posinp(1)))

        N_in  = 1
        N_fin = size(imgs)

      ELSE
         
         CALL SYSTEM( "./NEB_driver.sh free_only " // trim(job_name) // &
              & " " // trim(scratch_dir) // " " // trim(arr_posinp(1)))

        N_in  = 2
        N_fin = ( size(imgs) - 1 )

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

          imgs(replica)%outs%energy = temp_V

          DO i = 1, imgs(replica)%outs%fdim, 1

            READ(unit,*) imgs(replica)%outs%fxyz(1,i), &
                 imgs(replica)%outs%fxyz(2,i), &
                 imgs(replica)%outs%fxyz(3,i)

          END DO

          imgs(replica)%outs%fxyz = imgs(replica)%outs%fxyz * fix_atom

        END DO

      CLOSE( UNIT = unit )   

    END SUBROUTINE PES_IO

    SUBROUTINE deallocation
      use yaml_output
      use dynamic_memory
      use module_atoms, only: deallocate_atoms_data
      IMPLICIT NONE

      integer :: i, ierr

      IF ( ALLOCATED( fix_atom ) )         DEALLOCATE( fix_atom )

      if (allocated(imgs)) then
         do i = 1, size(imgs)
            call image_deallocate(imgs(i), .true.)
         end do
         deallocate(imgs)
      end if

      if (allocated(atoms)) then
         do i = 1, size(atoms)
            call deallocate_atoms_data(atoms(i))
         end do
         deallocate(atoms)
      end if

      if (allocated(ins)) then
         do i = 1, size(ins)
            call free_input_variables(ins(i))
         end do
         deallocate( ins )
      end if

      deallocate(arr_posinp,arr_radical)

      if (.not. external_call) then
         call free_restart_objects(rst, "deallocation")
         call f_lib_finalize()
      end if

      call mpi_environment_free(neb_mpi)
      call bigdft_finalize(ierr)
      call f_lib_initialize()
    END SUBROUTINE deallocation

END MODULE NEB_routines

PROGRAM NEB

  USE NEB_routines

  IMPLICIT NONE

  CALL read_input

  CALL search_MEP

  CALL deallocation

END PROGRAM NEB
