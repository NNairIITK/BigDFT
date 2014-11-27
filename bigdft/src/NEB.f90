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


!> Module for NEB calculations
MODULE NEB_routines
  use module_defs
  use module_types
  use module_images
  use module_interfaces
  use module_atoms, only: astruct_dump_to_file
  use dictionaries
  
  IMPLICIT NONE

  integer, parameter :: PES_INTERNAL_DFT   = 0
  integer, parameter :: PES_EXTERNAL_SHELL = 1
  integer, parameter :: PES_PAIR_POTENTIAL = 2

  type(NEB_data), private :: neb_

  LOGICAL :: external_call
  type(run_image), dimension(:), allocatable :: imgs
  integer :: pes_algo
  ! Additional variables for external call
  character (len=80) :: scratch_dir
  character (len=80) :: restart_file, velocity_file
  character (len=80) :: job_name

  character(len=max_field_length), dimension(:), allocatable :: arr_posinp
  type(mpi_environment), private :: neb_mpi

  CONTAINS

    SUBROUTINE read_input()
       use module_defs, only: BIGDFT_INPUT_VARIABLES_ERROR
      use yaml_output
      use dictionaries
      use module_interfaces
      use module_input_keys, only: input_keys_fill_all
      use module_input_dicts
      use module_atoms, only: atomic_structure, &
           deallocate_atomic_structure, &
           read_atomic_file => set_astruct_from_file, &
           astruct_nullify => nullify_atomic_structure
      use bigdft_run

      IMPLICIT NONE

      INTEGER :: i, num_of_images
      CHARACTER (LEN=max_field_length) :: minimization_scheme
      INTEGER, PARAMETER :: unit = 10
      integer :: ierr
      type(mpi_environment) :: bigdft_mpi_svg
      character(len=60) :: run_id
      type(dictionary), pointer :: dict, dict_min,options
      REAL (gp) :: tolerance
      character(len=60), dimension(:), allocatable :: arr_radical

      call f_lib_initialize()
      nullify(dict,options)
      !no options fof BigDFT
      call bigdft_command_line_options(options)
      call bigdft_init(options)
      call bigdft_get_run_properties(options // 'BigDFT' // 0, run_id = run_id)
      neb_mpi = bigdft_mpi!mpi_environment_null()
      neb_mpi%mpi_comm = MPI_COMM_NULL
      !this is redundant
      if (neb_mpi%nproc > 1) then
         call create_rank_comm(bigdft_mpi%mpi_comm, neb_mpi%mpi_comm)
      end if

!! default values are assigned
      call dict_init(dict)
      call read_input_dict_from_files(trim(run_id), bigdft_mpi, dict)
      call input_keys_fill_all(dict, dict_min)
      call neb_set_from_dict(dict, neb_%optimization, neb_%climbing, &
           & neb_%max_iterations, num_of_images, neb_%convergence, tolerance, &
           & neb_%ds, neb_%k_min, neb_%k_max, neb_%temp_req, neb_%damp, &
           & minimization_scheme)
      ! NEB is using cv criterion in ev per ang.
      neb_%convergence = neb_%convergence * Ha_eV / Bohr_Ang
      call dict_free(dict,dict_min)
      !call dict_free(dict_min)

      allocate(arr_radical(abs(num_of_images)))
      allocate(arr_posinp(abs(num_of_images)))
      call bigdft_get_run_ids(num_of_images,trim(run_id),arr_radical,arr_posinp,ierr)

      allocate( imgs(num_of_images) )

      
      call dict_copy(dict, options // 'BigDFT' // 0)
      call dict_free(options)
      ! Trick here, only super master will read the input files...
      bigdft_mpi_svg = bigdft_mpi
      bigdft_mpi%mpi_comm = MPI_COMM_WORLD
      call mpi_comm_rank(MPI_COMM_WORLD, bigdft_mpi%iproc, ierr)
      call mpi_comm_size(MPI_COMM_WORLD, bigdft_mpi%nproc, ierr)
      bigdft_mpi%igroup = 0
      bigdft_mpi%ngroup = num_of_images
      
      external_call = (bigdft_mpi%nproc == 1)

      !Loop over the images (replica)
      do i = 1, num_of_images
         bigdft_mpi%igroup = i
         !!!!<<<to substitute
         ! trick to output the image logs where it should, on disk.
         call user_dict_from_files(dict, trim(arr_radical(i)), &
              & trim(arr_posinp(i)), bigdft_mpi)
         !!!!!>>>>>end to substitute

         if (.not. external_call .and. i > 1) then
            call image_init(imgs(i), dict, trim(minimization_scheme), img0 = imgs(1))
         else
            call image_init(imgs(i), dict, trim(minimization_scheme))
         end if
         ! Store forces if present (for restart).
         call global_output_set_from_dict(imgs(i)%outs, dict // "posinp")
      end do
      call dict_free(dict)
      deallocate(arr_radical)

      job_name      = "neb"
      if (trim(run_id) /= "input") write(job_name, "(A)") trim(run_id)
      if (external_call) then
         restart_file  = trim(job_name) // ".NEB.restart"
         velocity_file = trim(job_name) // ".NEB.velocity"
         scratch_dir   = "./"
      end if
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

      if (bigdft_mpi%iproc == 0) then
         call images_init_path(imgs, tolerance, arr_posinp)
      else
         call images_init_path(imgs, tolerance)
      end if

!!$      ! Optionally pre-optimize the path
!!$      pes_algo = PES_PAIR_POTENTIAL
!!$      if (pes_algo == PES_PAIR_POTENTIAL) then
!!$         call search_MEP
!!$      end if
      if (external_call) then
         pes_algo = PES_EXTERNAL_SHELL
      else
         pes_algo = PES_INTERNAL_DFT
      end if

      ! End of trick.
      bigdft_mpi = bigdft_mpi_svg

!!$     END IF
    END SUBROUTINE read_input

    
    SUBROUTINE search_MEP
      use yaml_output
      use bigdft_run

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
      else if (bigdft_mpi%iproc == 0 .and. bigdft_mpi%igroup == 0) then
         call yaml_sequence_open("Restarting images", unit = 6)
         do i = 1, size(imgs), 1
            call yaml_sequence(trim(yaml_toa(all(imgs(i)%outs%fxyz /= UNINITIALIZED(1.d0)))), unit = 6, advance = "no")
            call yaml_comment(yaml_toa(i, fmt = "(I2.2)"), unit = 6)
            call yaml_newline(unit = 6)
         end do
         call yaml_sequence_close(unit = 6)
      end if

      if (bigdft_mpi%iproc == 0 .and. bigdft_mpi%igroup == 0) &
           & call yaml_sequence_open("NEB minimization loop", unit = 6)
      iteration = 0
      minimization: do
         if (pes_algo == PES_EXTERNAL_SHELL) then
            CALL PES_IO(imgs, (neb_%optimization .or. (.not. restart .and. iteration == 0)),stat)
            if (.not. stat) exit minimization
         else if (pes_algo == PES_INTERNAL_DFT) then
            call PES_internal(imgs, neb_%optimization,  (iteration == 0), iteration)
!!$         else if (pes_algo == PES_PAIR_POTENTIAL) then
!!$            call PES_pair(imgs, iteration)
         end if

         call compute_neb_pos(imgs, iteration, neb_)

         if (external_call) then
            call write_restart_vel(trim(velocity_file), imgs)
         end if

         if (bigdft_mpi%iproc == 0 .and. bigdft_mpi%igroup == 0 .and. iteration > 0) then
            call write_dat_files(trim(job_name), imgs, iteration)
            call yaml_swap_stream(6, unt, ierr)
            call yaml_sequence(advance='no')
            call images_output_step(imgs, iteration = iteration, tol = neb_%convergence)
            call yaml_set_default_stream(unt, ierr)
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
            call yaml_warning("Soft exit required, stopping")
            exit minimization
         END IF
      end do minimization
      if (bigdft_mpi%iproc == 0 .and. bigdft_mpi%igroup == 0) &
           & call yaml_sequence_close(unit = 6)

      if (bigdft_mpi%iproc == 0 .and. bigdft_mpi%igroup == 0) then
         call yaml_swap_stream(6, unt, ierr)
         call yaml_comment('Final results',hfill='-')
         call images_output_step(imgs, full = .true.)
         call yaml_set_default_stream(unt, ierr)

         do i = 1, size(imgs), 1
            filename=trim('final_'//trim(arr_posinp(i)))
            call bigdft_write_atomic_file(imgs(i)%run,imgs(i)%outs,&
                 filename,'FINAL CONFIGURATION',cwd_path=.true.)
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
         if (igroup(i) - 1 == bigdft_mpi%igroup) then
            call image_calculate(imgs(i), iteration, i)
        end if
      end do
      call images_collect_results(imgs, igroup, size(imgs), neb_mpi)
    END SUBROUTINE PES_internal

    SUBROUTINE PES_IO( imgs, flag , stat )

      IMPLICIT NONE

      type(run_image), dimension(:), intent(inout) :: imgs
      LOGICAL, INTENT(IN)        :: flag
      LOGICAL, INTENT(OUT)       :: stat
      INTEGER                    :: N_in, N_fin

      call write_restart(trim(restart_file), imgs)

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

      call images_PES_from_file(stat, imgs, "gen_output_file", N_in, N_fin)

    END SUBROUTINE PES_IO

    SUBROUTINE deallocation
      use yaml_output
      use dynamic_memory
      IMPLICIT NONE

      integer :: i

      if (allocated(imgs)) then
         do i = 1, size(imgs)
            call image_deallocate(imgs(i), .true.)
         end do
         deallocate(imgs)
      end if

      call mpi_environment_free(neb_mpi)
    END SUBROUTINE deallocation

END MODULE NEB_routines


PROGRAM NEB

  USE NEB_routines

  IMPLICIT NONE

  integer :: ierr

  CALL read_input()

  CALL search_MEP()

  CALL deallocation()

  call bigdft_finalize(ierr)
  call f_lib_finalize()  !call f_lib_initialize() <LG: most likely it was an error?
END PROGRAM NEB
