!> @file
!!    Energy and Forces for minima hopping guided path sampling
!! @author 
!!    Copyright (C) 2014 BigDFT group
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
subroutine mhgpsenergyandforces(mhgpsst,runObj,outs,rxyz,fxyz,fnoise,epot,infocode)
    !IMPORTANT:
    !receives distances in Bohr
    use module_base
    use yaml_output
    use module_mhgps_state
    use bigdft_run
    implicit none
    !parameters
    type(mhgps_state), intent(inout) :: mhgpsst
    type(run_objects), intent(inout) :: runObj
    type(state_properties), intent(inout) :: outs
    real(gp), intent(in) :: rxyz(3,runObj%atoms%astruct%nat)
    real(gp), intent(out) :: fxyz(3,runObj%atoms%astruct%nat)
    real(gp), intent(out) :: fnoise
    real(gp), intent(out) :: epot
    integer, intent(out)  :: infocode
    !internal

    call bigdft_set_rxyz(runObj,rxyz=rxyz)
    mhgpsst%ef_counter=mhgpsst%ef_counter+1.0_gp
    call bigdft_state(runObj,outs,infocode)
    call f_memcpy(src=outs%fxyz,dest=fxyz)
    epot=outs%energy
    fnoise=outs%fnoise
end subroutine

end module
