!! @file
!! @author Bastian Schaefer
!! @section LICENCE
!!    Copyright (C) 2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
program globaltool
    use module_base
    use module_types
    use module_interfaces
    use yaml_output
    use bigdft_run
    use module_globaltool
    implicit none
    type(gt_data) :: gdat
    type(dictionary), pointer :: options
    type(dictionary), pointer :: run
    character(len=60)         :: run_id, naming_id

    call f_lib_initialize()

    call bigdft_command_line_options(options)
    call bigdft_init(options)!mpi_info,nconfig,run_id,ierr)
    if (bigdft_nruns(options) > 1) then
        call f_err_throw('runs-file not supported for globaltools')
    endif
    run => options // 'BigDFT' // 0

    call bigdft_get_run_properties(run,input_id=run_id,&
         naming_id=naming_id)
    call bigdft_set_run_properties(run,&
         posinp_id=trim(adjustl(filename))//trim(naming_id))

    call run_objects_init(runObj,run)

    !options and run are not needed
    call dict_free(options)
    nullify(run)

    call yaml_new_document()

    call read_globaltool_uinp(gdat)
    call write_globaltool_uinp(gdat)

    call init_gt_data(gdat)

    call read_and_merge_data(gdat)

    call write_merged(gdat)

    call write_transitionpairs(gdat)

    call finalize_gt_data(gdat)

    call yaml_release_document()
    call f_lib_finalize()
end program globaltool
