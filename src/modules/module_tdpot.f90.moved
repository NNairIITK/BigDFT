!! @section LICENCE
!!    Copyright (C) 2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


module module_tdpot
    implicit none
    private

    public tdpot
    public init_tdpot

contains
!=====================================================================
subroutine init_tdpot(paramset,paramfile,units)
    use module_base
    use yaml_output
    implicit none
    !parameters
    character(len=*), intent(in) :: paramset
    character(len=*), intent(in) :: paramfile
    character(len=*), intent(in) :: units
    !local
    if (bigdft_mpi%iproc ==0) call yaml_comment('Initializing 2Dpot',hfill='-')
    if(units/= 'atomicd0' .and. &
        units/= 'atomic' .and. &
        units/= 'bohrd0' .and. &
        units/= 'bohr')then
        call f_err_throw('For 2Dpot, units must be atomic(d0) or bohr(d0)')
    endif
    !call yaml_scalar('Unit checks passed.')
    if(trim(paramfile)/='none')then
        call f_err_throw('Reading Parameters from file not '//&
             'implemented for 2Dpot')
    else
        select case(trim(paramset))
        case('default')
           if (bigdft_mpi%iproc ==0) call yaml_comment('Using normalized Units for LJ (sigma=1, epsilon=1)')
        case default
            call f_err_throw('Following parameter set for LJ force field '//&
                'is unknown: '//trim(paramset))
        end select
    endif
end subroutine
!=====================================================================
subroutine tdpot(nat,rxyz,fxyz,epot)
    use module_base
    use module_defs, only: gp
    use yaml_output
    !energy and forces for Lennard Jones potential
    !input: nat: number of atoms
    !       rxyz: positions of atoms
    !output: epot: energy
    !        fxyz: forces (negative derivative of energy with
    !              respect to positions
    implicit none
    !parameters
    integer, intent(in)   :: nat
    real(gp), intent(in)  :: rxyz(3)
    real(gp), intent(out) :: fxyz(3)
    real(gp), intent(out) :: epot
    !internal
    real(gp)  :: x,y,x2,y2,x3,y3,tt1,tt2,tt3

    if(nat/=1)then
            call f_err_throw('nat/=1 not allowed for 2Dpot')
    endif

    epot=0.d0
    fxyz(1)=0.d0 ; fxyz(2)=0.d0 ; fxyz(3)=0.d0

    x=rxyz(1)
    y=rxyz(2)
    x2=x**2
    y2=y**2
    x3=x2*x
    y3=y2*y
    tt1=x2+y2
    tt2=1.d0/(2.d0+tt1)
    tt3=tt2**2
    epot=x2*0.2d0+x3*0.1d0+y3*0.025d0+y2*tt2+(-2.d0+2*x2+y2)**2
    fxyz(1)=-0.1d0*(x*(-156.d0+x*(3.d0+160.d0*x)+y2*(80.d0-20.d0*tt3)))
    fxyz(2)=-(8.d0*(-1.d0+x2)*y+y2*0.075d0+4.d0*y3+2.d0*(2.d0+x2)*y*tt3)
end subroutine
end module
