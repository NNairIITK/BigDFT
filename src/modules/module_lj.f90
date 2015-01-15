!! @section LICENCE
!!    Copyright (C) 2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


module module_lj
    implicit none
    private

    public lenjon
    public init_lj

contains
!=====================================================================
subroutine init_lj(paramset,paramfile,units)
    use module_base                                                    
    use yaml_output                                                    
    implicit none                                                      
    !parameters                                                        
    character(len=*), intent(in) :: paramset                           
    character(len=*), intent(in) :: paramfile 
    character(len=*), intent(in) :: units
    !local                                                             
    call yaml_comment('Initializing LJ',hfill='-')
    if(units/= 'atomicd0' .and. &
        units/= 'atomic' .and. &
        units/= 'bohrd0' .and. &
        units/= 'bohr')then
        call f_err_throw('For LJ, units must be atomic(d0) or bohr(d0)') 
    endif
    call yaml_scalar('Unit checks passed.')
    if(trim(paramfile)/='none')then
        call f_err_throw('Reading Parameters from file not '//&
             'implemented for LJ')
    else
        select case(trim(paramset))
        case('default')
            call yaml_scalar('Using normalized Units for LJ (sigma=1, epsilon=1)')
        case default
            call f_err_throw('Following parameter set for LJ force field '//&       
                'is unknown: '//trim(paramset))
        end select
    endif
end subroutine 
!=====================================================================
subroutine lenjon(nat,rxyz,fxyz,epot)
    use module_defs, only: gp
    !energy and forces for Lennard Jones potential
    !input: nat: number of atoms
    !       rxyz: positions of atoms
    !output: epot: energy
    !        fxyz: forces (negative derivative of energy with
    !              respect to positions
    implicit none
    !parameters
    integer, intent(in)   :: nat
    real(gp), intent(in)  :: rxyz(3,nat)
    real(gp), intent(out) :: fxyz(3,nat)
    real(gp), intent(out) :: epot
    !internal
    integer :: iat, jat
    real(gp) :: dx,dy,dz,dd,dd2,dd6,dd12,tt,t1,t2,t3
    real(gp) :: xiat, yiat, ziat

    epot=0.d0
    do iat=1,nat
        fxyz(1,iat)=0.d0 ; fxyz(2,iat)=0.d0 ; fxyz(3,iat)=0.d0
    enddo

    !$omp parallel default(private) shared(nat, rxyz, fxyz, epot)
    !$omp do schedule(dynamic) reduction(+:fxyz,epot)
    do iat=1,nat
        xiat=rxyz(1,iat);yiat=rxyz(2,iat) ;ziat=rxyz(3,iat)
        do jat=1,iat-1
            dx=xiat-rxyz(1,jat)
            dy=yiat-rxyz(2,jat)
            dz=ziat-rxyz(3,jat)
            dd=dx**2+dy**2+dz**2
            dd2=1.d0/dd
            dd6=dd2*dd2*dd2
            dd12=dd6*dd6
            epot=epot+4.d0*(dd12-dd6)
            tt=24.d0*dd2*(2.d0*dd12-dd6)
            t1=dx*tt ; t2=dy*tt ; t3=dz*tt
            fxyz(1,iat)=fxyz(1,iat)+t1 ; fxyz(1,jat)=fxyz(1,jat)-t1
            fxyz(2,iat)=fxyz(2,iat)+t2 ; fxyz(2,jat)=fxyz(2,jat)-t2
            fxyz(3,iat)=fxyz(3,iat)+t3 ; fxyz(3,jat)=fxyz(3,jat)-t3
        enddo
    enddo
    !$omp end do
    !$omp end parallel
end subroutine
end module
