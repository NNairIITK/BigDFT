!> @file
!!    Energy and Forces for minima hopping guided path sampling
!! @author 
!!    Copyright (C) 2015-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Module returning energy and froces from minima hopping
module module_energyandforces
    implicit none

    private

    public :: mhgpsenergyandforces

contains


!> Returns energies in hartree and
!! forces in hartree/bohr
!! (except for LJ)
subroutine mhgpsenergyandforces(mhgpsst,runObj,outs,rxyz,fxyz,epot,infocode)
    !IMPORTANT:
    !receives distances in Bohr
    !use module_base
    !use yaml_output
    use module_mhgps_state
    use bigdft_run, only: bigdft_nat, bigdft_set_rxyz, bigdft_state
    use module_defs, only: gp
    use dynamic_memory
    use f_utils, only: char
    use yaml_output, only: yaml_warning
    implicit none
    !parameters
    type(mhgps_state), intent(inout) :: mhgpsst
    type(run_objects), intent(inout) :: runObj
    type(state_properties), intent(inout) :: outs
    !rxyz is modified if different MPI processes receive different coordinates:
    real(gp), intent(inout) :: rxyz(3,runObj%atoms%astruct%nat)
    real(gp), intent(out) :: fxyz(3,runObj%atoms%astruct%nat)
    real(gp), intent(out) :: epot
    integer, intent(out)  :: infocode
    !internal

    call bigdft_set_rxyz(runObj,rxyz=rxyz)
    mhgpsst%ef_counter=mhgpsst%ef_counter+1.0_gp
!debugging:
!if(mhgpsst%ef_counter <= 1.5d0)then
!    runObj%inputs%inputPsiId=0
!else
!    runObj%inputs%inputPsiId=1
!endif
    call bigdft_state(runObj,outs,infocode)
    if (bigdft_mpi%nproc >1) then
        call f_memcpy(src=runObj%atoms%astruct%rxyz,dest=rxyz)
    end if
    call f_memcpy(src=outs%fxyz,dest=fxyz)
    epot=outs%energy
end subroutine

end module
