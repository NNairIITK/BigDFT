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
    use module_atoms, only: astruct_dump_to_file
    use module_mhgpstool
    use module_userinput, read_mhgps_input => read_input
    implicit none
    type(atoms_data) :: atoms
    type(userinput) :: mhgps_uinp
    type(mhgpstool_data) :: mdat
    integer :: nat
    integer          :: nfolder
    character(len=500), allocatable :: folders(:)
    character(len=600) :: filename
    integer, allocatable :: nsad(:)
!character(len=500) :: fsaddle,comment
!real(gp) :: energy
!real(gp) :: rxyz(3,38)


    call f_lib_initialize()

    call yaml_new_document()

    call read_folders(nfolder,folders)
    nsad = f_malloc((/ 1.to.nfolder/),id='nsad')
    write(filename,'(a,i5.5,a)')trim(adjustl(folders(1))),1,'_finalF'
    nat = bigdft_nat(filename=filename)
    call count_saddle_points(nfolder,folders,nsad)
    call init_mhgpstool_data(nat,nfolder,nsad,mdat)
    call read_data(folders,mdat)

    call read_mhgps_input(mhgps_uinp)
!write(fsaddle,'(a,i5.5,a)')trim(adjustl(folders(1)))//&
!                               '/sad',1,'_finalF'
!    call bigdft_get_rxyz(filename=trim(adjustl(fsaddle)),rxyz=rxyz,energy=energy)
!write(*,*)energy

    call finalize_mhgpstool_data(mdat)
    call f_free_str(500,folders)
    call f_free(nsad)

    call f_lib_finalize()
    call yaml_release_document()
end program mhgpstool
