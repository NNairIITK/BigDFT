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
    use SPREDtypes
    use module_globaltool
    implicit none
    type(gt_data) :: gdat

    call f_lib_initialize()

    call yaml_new_document()

    call SPRED_read_uinp('',gdat%spredinputs,bigdft_mp)

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
