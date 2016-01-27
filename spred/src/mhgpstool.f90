!! @file
!! @author Bastian Schaefer
!! @section LICENCE
!!    Copyright (C) 2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
program mhgpstool
    use module_base
    use module_types
    use module_interfaces
    use yaml_output
    use bigdft_run
    use module_atoms, only: set_astruct_from_file
    use module_mhgpstool
    use module_userinput, read_mhgps_input => read_input
    use SPREDtypes
    implicit none
    type(atoms_data) :: atoms
    type(mhgpstool_data) :: mdat
    integer :: nat
    integer          :: nfolder
    character(len=500), allocatable :: folders(:)
    character(len=600) :: filename
    integer, allocatable :: nsad(:)
    real(gp) :: energy
    type(SPRED_inputs) :: spredinputs

    call f_lib_initialize()

    call yaml_new_document()

    call SPRED_read_uinp('',spredinputs,bigdft_mpi)

    call read_folders(nfolder,folders)
    call read_mhgps_input(mdat%mhgps_uinp)
    nsad = f_malloc((/ 1.to.nfolder/),id='nsad')
    write(filename,'(a,i5.5,a)')trim(adjustl(folders(1)))//'/sad',1,'_finalF'
    nat = bigdft_nat(filename=filename)
    call count_saddle_points(nfolder,folders,nsad)
    call init_mhgpstool_data(nat,nfolder,nsad,mdat)

    call set_astruct_from_file(trim(filename),0,mdat%astruct,energy=energy)
    call yaml_comment('Covalent radii ....',hfill='-')
    call give_rcov(0,mdat%astruct,mdat%astruct%nat,mdat%rcov)
    call read_and_merge_data(spredinputs,folders,nsad,mdat)

!write(fsaddle,'(a,i5.5,a)')trim(adjustl(folders(1)))//&
!                               '/sad',1,'_finalF'
!    call bigdft_get_rxyz(filename=trim(adjustl(fsaddle)),rxyz=rxyz,energy=energy)
!write(*,*)energy
    call write_data(spredinputs,mdat)

    call finalize_mhgpstool_data(mdat)
    call f_free_str(500,folders)
    call f_free(nsad)

    call yaml_release_document()
    call f_lib_finalize()
end program mhgpstool
