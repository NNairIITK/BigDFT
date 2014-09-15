!! @file
!! @author Bastian Schaefer
!! @section LICENCE
!!    Copyright (C) 2014 UNIBAS
!!    This file is not freely distributed.
!!    A licence is necessary from UNIBAS

module module_energyandforces
    implicit none

    private

    public :: energyandforces

contains
!=====================================================================
subroutine energyandforces(nat,alat,rxyz,fxyz,fnoise,epot)
    !returns energies in hartree and
    !forces in hartree/bohr
    !(except for LJ)
    use module_base
    use module_lj
    use module_lenosky_si
    use module_types
    use module_interfaces
    use yaml_output
    use module_global_variables
    implicit none
    !parameters
    integer, intent(in) :: nat
    real(gp), intent(in) :: alat(3)
    real(gp), intent(in) :: rxyz(3,nat)
    real(gp), intent(out) :: fxyz(3,nat)
    real(gp), intent(out) :: fnoise
    real(gp), intent(out) :: epot
    !internal
    integer :: icc !for amber
    real(gp) :: rxyzint(3,nat)
    real(gp) :: alatint(3)
    if(nat/=fdim)stop 'nat /= fdim'
    ef_counter=ef_counter+1.0_gp 
    if(trim(adjustl(efmethod))=='LJ')then
        call lenjon(nat,rxyz(1,1),fxyz(1,1),epot)
        fnoise=0.0_gp
        return
    else if(trim(adjustl(efmethod))=='LENSIc')then!for clusters
        !convert from bohr to ansgtroem
        rxyzint=0.52917721092_gp*rxyz
        alatint=0.52917721092_gp*alat
        call lenosky_si_shift(nat,alatint,rxyzint(1,1),fxyz(1,1),epot)
        !convert energy from eV to Hartree
        epot=0.03674932379085202_gp * epot
        !convert forces from eV/Angstroem to hartree/bohr
        fxyz(1:3,1:nat)=fxyz(1:3,1:nat)*0.01944690466683907_gp
        fnoise=0.0_gp
        return
    else if(trim(adjustl(efmethod))=='LENSIb')then!for bulk
        !convert from bohr to ansgtroem
        rxyzint=0.52917721092_gp*rxyz
        alatint=0.52917721092_gp*alat
        call lenosky_si(nat,alatint,rxyzint,fxyz,epot)
        !convert energy from eV to Hartree
        epot=0.03674932379085202_gp * epot
        !convert forces from eV/Angstroem to hartree/bohr
        fxyz(1:3,1:nat)=fxyz(1:3,1:nat)*0.01944690466683907_gp
        fnoise=0.0_gp
        return
!    else if(trim(adjustl(efmethod))=='AMBER')then
!        icc=1
!        !convert from bohr to ansgtroem
!        rxyzint=0.52917721092_gp*rxyz
!        call call_nab_gradient(rxyzint(1,1),fxyz(1,1),epot,icc)
!        epot=epot*0.001593601437458137_gp !from kcal_th/mol to hartree
!                                          !(thermochemical calorie
!                                          !used: 1cal_th=4.184J)
!                                          !also see:
!                          !http://archive.ambermd.org/201009/0039.html
!        !convert from gradient in kcal_th/mol/angstrom to
!        !force in hartree/bohr
!        fxyz(1:3,1:nat)=-fxyz(1:3,1:nat)*0.0008432975639921999_gp
!        fnoise=0.0_gp
!        return
    else if(trim(adjustl(efmethod))=='BIGDFT')then
        if(nat/=runObj%atoms%astruct%nat)then
            call yaml_warning('nat /= runObj%atoms%astruct%nat in '//&
                              'energyandforces')
            stop
        endif
        call vcopy(3 * runObj%atoms%astruct%nat, rxyz(1,1),1,&
             runObj%atoms%astruct%rxyz(1,1), 1)
        runObj%inputs%inputPsiId=inputPsiId
        runObj%inputs%itermin=itermin
        call call_bigdft(runObj,outs,infocode)
        call vcopy(3 * outs%fdim, outs%fxyz(1,1), 1, fxyz(1,1), 1)
        call vcopy(3 * runObj%atoms%astruct%nat,&
             runObj%atoms%astruct%ixyz_int(1,1),1,&
             ixyz_int(1,1), 1)
        epot=outs%energy
        fnoise=outs%fnoise
        return
    endif
end subroutine
end module
