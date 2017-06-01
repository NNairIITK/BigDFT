!! @section LICENCE                                                    
!!    Copyright (C) 2015 BigDFT group                                  
!!    This file is distributed under the terms of the                  
!!    GNU General Public License, see ~/COPYING file                   
!!    or http://www.gnu.org/copyleft/gpl.txt .                         
!!    For the list of contributors, see ~/AUTHORS

!   
!   Copyright (C) 1999-2006 David J. Wales
!   This file was part of GMIN.
!
!   GMIN is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   GMIN is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!*************************************************************************
!
!  Subroutine MORSE_bulk calculates the energy and forces analytically for 
!  the Morse potential with periodic boundary conditions. The potential has
!  a cutoff and is shifted to make the potential continuous.  It is not smooth
!  though.
!
!  
!  rho is the inverse width of the well. rcut is the cutoff.
!
!*************************************************************************
!
module module_morse_bulk
    use module_base
    implicit none

    !default private
    private

    !morse_bulk parameters
    real(gp),save :: rho
    real(gp),save :: rcut
    real(gp),save :: R0
    real(gp),save :: A
    logical, save :: periodic
    logical, save :: use_cutoff
    logical,save :: initialized_bulk=.false.
    logical,save :: initialized_slab=.false.

    real(gp), save :: alat_int(3)

    public morse_bulk_wrapper
    public init_morse_bulk
    public morse_slab_wrapper
    public init_morse_slab
    contains 
subroutine init_morse_slab(paramset,paramfile,geocode,units)
    use module_base
    use yaml_output
    implicit none
    !parameters
    character(len=*), intent(in) :: paramset
    character(len=*), intent(in) :: paramfile
    character(len=*), intent(in) :: geocode
    character(len=*), intent(in) :: units
    !local
    integer :: u
    logical :: exists
    call yaml_comment('Initializing Morse_Slab',hfill='-')


    initialized_slab=.false.

    if(trim(geocode)/='S')then
        initialized_slab=.false.
        call f_err_throw('Morse_Slab only works with surface '//&
             'boundary conditions. Specified boundary conditions are: '//&
             trim(adjustl(geocode)))
    endif

    select case(trim(paramset))
    case('Pt')
        !initialized_slab=false, because we still need spacing between slabs
        initialized_slab=.false.
        call yaml_mapping_open('Using Pt Parameters from'//&
             ' Bassett, D. W.; Webber, P. R. Surf. Sci. 1978, 70, 520.')
        rho = 1.6047_gp * Bohr_Ang !convert 1.6047 A^-1 to  1/Bohr
        rcut = 9.5_gp / Bohr_Ang !convert 9.5 Angstroem to Bohr
        R0 = 2.8970_gp / Bohr_Ang !convert 2.8960 Angstroem to Bohr
        A = 0.7102_gp * eV_Ha !convert 0.7102 eV to Hartree
        call yaml_map('rho (1/Bohr)', rho,  fmt='(1pe10.4)')
        call yaml_map('rcut (Bohr)',  rcut, fmt='(1pe10.4)')
        call yaml_map('R0 (Bohr)',    R0,   fmt='(1pe10.4)')
        call yaml_map('A (Hartree)',  A,    fmt='(1pe10.4)')
        call yaml_mapping_close()
    case('default')
        initialized_slab=.false.
        call f_err_throw('No "default" parameter set for morse_bulk defined.')
    case default
        initialized_slab=.false.
        call f_err_throw('Following parameter set for morse_bulk force field '//&       
            'is unknown: '//trim(paramset))
    end select
    if(trim(paramfile)/='none')then
        inquire(file=trim(adjustl(paramfile)),exist=exists)
        if(.not.exists)then
            call f_err_throw('Parameter file '//trim(adjustl(paramfile))//&
                 ' does not exist.')
        endif
        u=f_get_free_unit()
        open(unit=u,file=trim(adjustl(paramfile)))
            read(u,*) alat_int(2)
        close(u)
        if (units=='angstroem' .or. units=='angstroemd0') then         
           ! if Angstroem convert to Bohr                                    
           alat_int(1)=alat_int(1)/Bohr_Ang                            
           alat_int(2)=alat_int(2)/Bohr_Ang                            
           alat_int(3)=alat_int(3)/Bohr_Ang                            
        else if  (units=='atomic' .or. units=='bohr'  .or.&            
             units== 'atomicd0' .or. units== 'bohrd0') then            
             !do nothing                                               
        else if (units == 'reduced') then                              
           !assume that for reduced coordinates cell size is in bohr         
        else                                                           
           call f_err_throw('Length units in input file unrecognized.' // &
                'Recognized units are angstroem or atomic = bohr',err_id=BIGDFT_INPUT_VARIABLES_ERROR)
           return                                                      
           !write(*,*) 'length units in input file unrecognized'             
           !write(*,*) 'recognized units are angstroem or atomic = bohr'     
           !stop                                                             
        endif  
        call yaml_map('Read cell parameters from file',trim(adjustl(paramfile)))
        
        initialized_slab=.true.
    else
        call f_err_throw('Parameter file for slab spacing is missing.')
    endif
    
    periodic = .true.
    use_cutoff = .true.
    if (rcut / R0 < 1.e-1_gp .and. bigdft_mpi%iproc==0) then
       call yaml_warning('morse_bulk: warning the cutoff is very small')
    endif
end subroutine init_morse_slab
subroutine init_morse_bulk(paramset,paramfile,geocode)
    use module_base
    use yaml_output
    implicit none
    !parameters
    character(len=*), intent(in) :: paramset
    character(len=*), intent(in) :: paramfile
    character(len=*), intent(in) :: geocode
    !local
    call yaml_comment('Initializing Morse_Bulk',hfill='-')


    initialized_bulk=.true.

    if(trim(geocode)/='P')then
        initialized_bulk=.false.
        call f_err_throw('Morse_bulk only works with periodic '//&
             'boundary conditions. Specified boundary conditions are: '//&
             trim(adjustl(geocode)))
    endif

    if(trim(paramfile)/='none')then
        initialized_bulk=.false.
        call f_err_throw('Reading Parameters from file not '//&
             'implemented for morse_bulk')
    else
        select case(trim(paramset))
        case('Pt')
            call yaml_mapping_open('Using Pt Parameters from'//&
                 ' Bassett, D. W.; Webber, P. R. Surf. Sci. 1978, 70, 520.')
            rho = 1.6047_gp * Bohr_Ang !convert 1.6047 A^-1 to  1/Bohr
            rcut = 9.5_gp / Bohr_Ang !convert 9.5 Angstroem to Bohr
            R0 = 2.8970_gp / Bohr_Ang !convert 2.8960 Angstroem to Bohr
            A = 0.7102_gp * eV_Ha !convert 0.7102 eV to Hartree
            call yaml_map('rho (1/Bohr)', rho,  fmt='(1pe10.4)')
            call yaml_map('rcut (Bohr)',  rcut, fmt='(1pe10.4)')
            call yaml_map('R0 (Bohr)',    R0,   fmt='(1pe10.4)')
            call yaml_map('A (Hartree)',  A,    fmt='(1pe10.4)')
            call yaml_mapping_close()
        case('default')
            initialized_bulk=.false.
            call f_err_throw('No "default" parameter set for morse_bulk defined.')
        case default
            initialized_bulk=.false.
            call f_err_throw('Following parameter set for morse_bulk force field '//&       
                'is unknown: '//trim(paramset))
        end select
    endif
    periodic = .true.
    use_cutoff = .true.
    if (rcut / R0 < 1.e-1_gp .and. bigdft_mpi%iproc==0) then
       call yaml_warning('morse_bulk: warning the cutoff is very small')
    endif
end subroutine init_morse_bulk
subroutine morse_slab_wrapper(nat,alat,rxyz, fxyz, epot)
    use module_base
    implicit none 
    !parameter
    integer, intent(in) :: nat
    real(gp), intent(in) :: alat(3)
    real(gp), intent(in) :: rxyz(3*nat)
    real(gp), intent(out) :: fxyz(3*nat), epot
    !local
    if(.not.initialized_slab)then
        call f_err_throw('Potential "morse_slab" not initialized',&
             err_name='BIGDFT_RUNTIME_ERROR')
    endif
    
    alat_int(1)=alat(1)
    alat_int(3)=alat(3)

    call morse_bulk(rxyz(1),fxyz(1),epot, nat, rho, R0, A, periodic, & 
       alat_int, use_cutoff, rcut)
end subroutine morse_slab_wrapper
subroutine morse_bulk_wrapper(nat,alat,rxyz, fxyz, epot)
    use module_base
    implicit none 
    !parameter
    integer, intent(in) :: nat
    real(gp), intent(in) :: alat(3)
    real(gp), intent(in) :: rxyz(3*nat)
    real(gp), intent(out) :: fxyz(3*nat), epot
    !local
    if(.not.initialized_bulk)then
        call f_err_throw('Potential "morse_bulk" not initialized',&
             err_name='BIGDFT_RUNTIME_ERROR')
    endif
    
    call morse_bulk(rxyz(1),fxyz(1),epot, nat, rho, R0, A, periodic, & 
       alat, use_cutoff, rcut)
end subroutine morse_bulk_wrapper

SUBROUTINE MORSE_BULK(X,V,EMORSE, natoms, rho, R0, A, periodic, &
         boxvec, use_cutoff, rcut)
      ! R0 is the position of the bottom of the well
      ! rho is the width of the well and has units of inverse length
      ! A is the energy scale
      implicit none 
      logical, intent(in) :: periodic, use_cutoff
      integer, intent(in) :: NATOMS
      real(gp), intent(in) :: X(3*NATOMS), rho, R0, A, boxvec(3), rcut
      real(gp), intent(out) :: V(3*NATOMS), EMORSE
      integer ::  J1, J2, J3, J4
      real(gp) :: DIST, R, DUMMY, &
                       RR(NATOMS,NATOMS), &
                       XMUL2, iboxvec(3), dx(3), eshift
      if (periodic) iboxvec(:) = 1.0_gp / boxvec(:)

      if (use_cutoff) then
         Eshift = (1.0_gp - exp(rho * (r0 - rcut)))**2 - 1.0_gp
         !write(*,*) "Eshift", eshift, rcut
      else
        Eshift=0.0_gp
      endif

      V(:) = 0.0_gp
      EMORSE=0.0_gp
      !$omp parallel default(private) shared(natoms, x, V, emorse, &
      !$omp boxvec, iboxvec, RHO, R0, periodic, use_cutoff, rcut, &
      !$omp Eshift, A)
      !$omp do schedule(dynamic) reduction(+:V,emorse)
      DO J1=1,NATOMS
         J3=3*J1
         RR(J1,J1)=0.0_gp
         DO J2=J1+1,NATOMS
            J4=3*J2
            dx(:) = X(J3-2:j3)-X(J4-2:j4)
            if (periodic) then
               dx = dx - boxvec * nint(dx * iboxvec)
            endif
            dist = max(sqrt(sum(dx**2)), 1.0e-5_gp)

            if (use_cutoff .and. dist.ge.rcut) cycle

            R=exp(RHO*R0-RHO*DIST)
            DUMMY=R*(R-2.0_gp)
            EMORSE=EMORSE+DUMMY - Eshift
!            if (gtest) then
               xmul2 = 2.0_gp*R*(R-1.0_gp)/DIST * A
!               V(J3-2:j3) = V(j3-2:j3) - xmul2 * dx
!               V(J4-2:j4) = V(j4-2:j4) + xmul2 * dx
               V(J3-2:j3) = V(j3-2:j3) + xmul2 * dx
               V(J4-2:j4) = V(j4-2:j4) - xmul2 * dx
!            endif
         ENDDO
      ENDDO
      !$omp end do
      !$omp end parallel
      EMORSE = EMORSE * A

      RETURN
END subroutine morse_bulk
end module module_morse_bulk
