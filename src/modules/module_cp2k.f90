!! @section LICENCE                                                    
!!    Copyright (C) 2015 BigDFT group                                  
!!    This file is distributed under the terms of the                  
!!    GNU General Public License, see ~/COPYING file                   
!!    or http://www.gnu.org/copyleft/gpl.txt .                         
!!    For the list of contributors, see ~/AUTHORS

!*************************************************************************
!
!  Subroutine cp2k_energy_forces calculates the energy and forces using the cp2k
!  code. 
!*************************************************************************
!
module module_cp2k
    use module_base
    implicit none

    !default private
    private

    public :: init_cp2k
    public :: finalize_cp2k
    public :: cp2k_energy_forces

    !parameters
    integer, save :: f_env_id
    logical,save :: initialized_cp2k=.false.
    character(len=1), save :: geocode='F'

    contains 
subroutine init_cp2k(paramfile,geocodeIn)
    use module_base
    use yaml_output
    implicit none
    !parameters
    character(len=*), intent(in) :: paramfile
    character(len=1), intent(in) :: geocodeIn
    !local
    logical :: exists
    integer :: ierr

    if(initialized_cp2k)stop'cp2k already initalized'
    if(bigdft_mpi%iproc==0)call yaml_comment('Initializing cp2k',hfill='-')

    initialized_cp2k=.false.

    if(trim(paramfile)/='none')then
        inquire(file=trim(adjustl(paramfile)),exist=exists)
        if(.not.exists)then
            call f_err_throw('Parameter file '//trim(adjustl(paramfile))//&
                 ' does not exist.')
        endif
        CALL cp_init_cp2k(0,ierr)
!!        CALL cp_init_cp2k(1,ierr)
        CALL cp_create_fenv(f_env_id,trim(adjustl(paramfile)),"cp2k.out",ierr)
        geocode=geocodeIn
        initialized_cp2k=.true.
    else
        call f_err_throw('CP2K input file is not specified in input.yaml.')
    endif
end subroutine init_cp2k
subroutine finalize_cp2k()
    use module_base
    use yaml_output
    implicit none
    !parameters
    !internal
    integer :: ierr
    character(len=5) :: cierr
    if(bigdft_mpi%iproc==0)call yaml_comment('Finalizing CP2K',hfill='-')
    CALL cp_finalize_cp2k(0,ierr)
    if (ierr/=0)then
        write(cierr,'(i5.5)')
        call f_err_throw('Error while finalizing cp2k, ierr: '//trim(adjustl(cierr)))
    endif
end subroutine
subroutine cp2k_energy_forces(nat,alat,rxyz, fxyz, epot,istat)
    !receives and returns atomic units
    use module_base
    implicit none 
    !parameter
    integer, intent(in) :: nat
    real(gp), intent(in) :: alat(3)
    real(gp), intent(in) :: rxyz(3*nat)
    real(gp), intent(out) :: fxyz(3*nat), epot
    integer, intent(out) :: istat
    !local
    real(gp) :: cell(3,3)
    integer:: ierr
    character(len=5) :: cierr
integer :: clck_counts_beg, clck_counts_end, clck_rate
call system_clock ( clck_counts_beg, clck_rate )
    if(.not.initialized_cp2k)then
        call f_err_throw('Method "cp2k" not initialized',&
             err_name='BIGDFT_RUNTIME_ERROR')
    endif
!write(199,*)
!do i=1,nat
!write(199,'(3(x,es24.17))')rxyz(3*i-2),rxyz(3*i-1),rxyz(3*i)
!enddo
    
    CALL cp_set_pos(f_env_id,rxyz,SIZE(rxyz),ierr)
    if (ierr/=0)then
        write(cierr,'(i5.5)')
        call f_err_throw('Error in CP2K while setting positions, ierr: '//trim(adjustl(cierr)))
    endif
    if(geocode/='F') then
        cell(1,1)=alat(1)
        cell(2,1)=0.0_gp
        cell(3,1)=0.0_gp
        cell(1,2)=0.0_gp
        cell(2,2)=alat(2)
        cell(3,2)=0.0_gp
        cell(1,3)=0.0_gp
        cell(2,3)=0.0_gp
        cell(3,3)=alat(3)
        CALL cp_set_cell(f_env_id,cell, ierr)
    endif
    if (ierr/=0)then
        write(cierr,'(i5.5)')
        call f_err_throw('Error in CP2K while setting cell, ierr: '//trim(adjustl(cierr)))
    endif
    CALL cp_calc_energy_force(f_env_id,1,istat)
    CALL cp_get_energy(f_env_id,epot,ierr)
    if (ierr/=0)then
        write(cierr,'(i5.5)')
        call f_err_throw('Error in CP2K while getting energies, ierr: '//trim(adjustl(cierr)))
    endif
    CALL cp_get_force(f_env_id,fxyz,SIZE(fxyz),ierr)
    if (ierr/=0)then
        write(cierr,'(i5.5)')
        call f_err_throw('Error in CP2K while getting forces, ierr: '//trim(adjustl(cierr)))
    endif
call system_clock ( clck_counts_end, clck_rate) 
write(333,*)(clck_counts_end - clck_counts_beg) / real (clck_rate), "seconds"
end subroutine cp2k_energy_forces

end module module_cp2k
