!! @file
!! @author Bastian Schaefer
!! @section LICENCE
!!    Copyright (C) 2014 UNIBAS
!!    This file is not freely distributed.
!!    A licence is necessary from UNIBAS

module module_energyandforces

contains

subroutine energyandforces(nat,alat,rxyz,fxyz,epot)
    use module_base
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
    real(gp), intent(out) :: epot
    !internal

    ef_counter=ef_counter+1.d0    
    if(trim(adjustl(efmethod))=='LJ')then
        call lenjon(nat,rxyz(1,1),fxyz(1,1),epot)
        return
    else if(trim(adjustl(efmethod))=='BIGDFT')then
        if(nat/=runObj%atoms%astruct%nat)then
            call yaml_warning('nat /= runObj%atoms%astruct%nat in energyandforces')
            stop
        endif
        call vcopy(3 * runObj%atoms%astruct%nat, rxyz(1,1),1,runObj%atoms%astruct%rxyz(1,1), 1)
        runObj%inputs%inputPsiId=inputPsiId
        runObj%inputs%itermin=itermin
        call call_bigdft(runObj,outs,bigdft_mpi%nproc,bigdft_mpi%iproc,infocode)
        call vcopy(3 * outs%fdim, outs%fxyz(1,1), 1, fxyz(1,1), 1)
        call vcopy(3 * runObj%atoms%astruct%nat, runObj%atoms%astruct%ixyz_int(1,1), 1, ixyz_int(1,1), 1)
        epot=outs%energy
        return
    endif
end subroutine

subroutine lenjon(nat,rxyz,fxyz,etot)
    use module_base
    !energy and forces for Lennard Jones potential
    !input: nat: number of atoms
    !       rxyz: positions of atoms
    !output: etot: energy
    !        fxyz: forces (negative derivative of energy with respect to
    !        positions
    implicit none
    !parameters
    integer, intent(in)   :: nat
    real(gp), intent(in)  :: rxyz(3,nat)
    real(gp), intent(out) :: fxyz(3,nat)
    real(gp), intent(out) :: etot
    !internal
    integer :: iat, jat
    real(g) :: dx,dy,dy,dd,dd2,dd6,dd12,tt,t1,t2,t3

    etot=0.d0
    do iat=1,nat
        fxyz(1,iat)=0.d0 ; fxyz(2,iat)=0.d0 ; fxyz(3,iat)=0.d0
    enddo
    do iat=1,nat
        do jat=1,iat-1
            dx=rxyz(1,iat)-rxyz(1,jat)
            dy=rxyz(2,iat)-rxyz(2,jat)
            dz=rxyz(3,iat)-rxyz(3,jat)
            dd=dx**2+dy**2+dz**2
            dd2=1.d0/dd
            dd6=dd2*dd2*dd2
            dd12=dd6*dd6
            etot=etot+4.d0*(dd12-dd6)
            tt=24.d0*dd2*(2.d0*dd12-dd6)
            t1=dx*tt ; t2=dy*tt ; t3=dz*tt
            fxyz(1,iat)=fxyz(1,iat)+t1 ; fxyz(1,jat)=fxyz(1,jat)-t1
            fxyz(2,iat)=fxyz(2,iat)+t2 ; fxyz(2,jat)=fxyz(2,jat)-t2
            fxyz(3,iat)=fxyz(3,iat)+t3 ; fxyz(3,jat)=fxyz(3,jat)-t3
        enddo
    enddo
end subroutine


end module
