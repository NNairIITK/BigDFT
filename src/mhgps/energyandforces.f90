module module_energyandforces

subroutine energyandforces(nat,alat,rxyz,fxyz,epot)
    use module_types
    implicit none
    
    if(trim(adjustl(efmethod))=='LJ')then
        call lenjon(runObj%atoms%astruct%nat,runObj%atoms%astruct%rxyz(1,1),outs%fxyz(1,1),outs%energy)
    else if(efmethod=='BIGDFT')then
        call call_bigdft(runObj,outs,bigdft_mpi%nproc,bigdft_mpi%iproc,infocode)
    endif
end subroutine

subroutine

subroutine lenjon(nat,rxyz,fxyz,etot)
    use module_base
    !energy and forces for Lennard Jones potential
    !input: nat: number of atoms
    !       rxyz: positions of atoms
    !output: etot: energy
    !        fxyz: forces (negative derivative of energy with respect to
    !        positions
    implicit real(gp) (a-h,o-z)
    dimension rxyz(3,nat),fxyz(3,nat)

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
    return
end


end module
