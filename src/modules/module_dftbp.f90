!! @section LICENCE                                                    
!!    Copyright (C) 2015 BigDFT group                                  
!!    This file is distributed under the terms of the                  
!!    GNU General Public License, see ~/COPYING file                   
!!    or http://www.gnu.org/copyleft/gpl.txt .                         
!!    For the list of contributors, see ~/AUTHORS

!*************************************************************************
!
!Interface to the DFTB+ code. !!OVER FILES!!
!
!*************************************************************************
!
module module_dftbp
    use module_base
    implicit none

    !default private
    private

    public :: dftbp_energy_forces

    contains 
subroutine dftbp_energy_forces(policy,nat,alat,astruct,geocode,rxyz, fxyz,strten,epot,istat)
    !receives and returns atomic units
    use module_base
    use module_types
    implicit none 
    !parameter
    integer, intent(in) :: policy
    integer, intent(in) :: nat
    real(gp), intent(in) :: alat(3)
    character(len=1), intent(in) :: geocode
    type(atomic_structure), intent(in) :: astruct
    real(gp), intent(in) :: rxyz(3*nat)
    real(gp), intent(out) :: fxyz(3*nat),strten(6), epot
    integer, intent(out) :: istat
    !local
    real(gp) :: cell(3,3)
    integer:: ierr
    character(len=5) :: cierr
!!integer :: clck_counts_beg, clck_counts_end, clck_rate
!!call system_clock ( clck_counts_beg, clck_rate )
    istat=0

    
    cell(1,1)=alat(1)
    cell(2,1)=0.0d0
    cell(3,1)=0.0d0
    cell(1,2)=0.0d0
    cell(2,2)=0.0d0
    cell(3,2)=alat(2)
    cell(1,3)=0.0d0
    cell(2,3)=0.0d0
    cell(3,3)=alat(3)
    call create_input_files(policy,nat,cell,astruct,geocode,rxyz)
    call run_dftb()
    call get_results(nat,geocode,fxyz,strten,epot)
!!call system_clock ( clck_counts_end, clck_rate) 
!!write(333,*)(clck_counts_end - clck_counts_beg) / real (clck_rate), "seconds"
end subroutine dftbp_energy_forces

subroutine get_results(nat,geocode,fxyz,strten,epot)
    use module_base
    implicit none
    !parameters
    integer, intent(in) :: nat
    character(len=1), intent(in) :: geocode
    real(gp), intent(out) :: fxyz(3,nat),strten(6)
    real(gp), intent(out) :: epot
    !internal
    integer :: u,n,k,m,l,iat
    real(gp) :: str_matrix(3,3)
    character(len=250) :: all_line

    epot=huge(1.0_gp)
    fxyz=huge(1.0_gp)
    strten=huge(1.0_gp)

    u=f_get_free_unit()
    open(unit=u,file="detailed.out")
    do
        read(u,'(a250)',end=99)all_line
        n = len_trim(all_line)
        k = index(all_line(1:n),"Total Mermin free energy")
        if(k.ne.0) then
            m = len_trim(all_line)
            l = scan(all_line(1:m),":",.true.)
            read(all_line(l+1:m),*) epot
            cycle
        endif
        k = index(all_line(1:n),"Total Forces")
        if(k.ne.0) then
          do iat=1,nat
              read(u,*) fxyz(:,iat)
          enddo
          cycle
        endif
        k = index(all_line(1:n),"Total stress tensor")
        if(k.ne.0) then
            do iat=1,3
              read(u,*) str_matrix(:,iat)
            enddo
            strten(1)=-str_matrix(1,1)
            strten(2)=-str_matrix(2,2)
            strten(3)=-str_matrix(3,3)
            strten(6)=-str_matrix(1,2)
            strten(5)=-str_matrix(1,3)
            strten(4)=-str_matrix(2,3)
            cycle
        endif
    enddo

    99 continue
    close(u)
    if(geocode=='F') strten=0.d0
    if(epot>=1.e10_gp.or.strten(1)>=1.e10_gp.or.fxyz(1,1)>=1.e10_gp)then
        call f_err_throw('DFTB+ interface: Could not find all requested variables in dftb+ output file')
    endif

end subroutine

subroutine run_dftb()
    implicit none
    call system('./rundftbp.sh')
end subroutine

subroutine create_input_files(policy,nat,cell,astruct,geocode,rxyz)
    use module_base
    use yaml_output
    use module_types
    !parameter
    integer, intent(in) :: policy
    integer, intent(in) :: nat
    real(gp), intent(in) :: cell(3,3)
    type(atomic_structure), intent(in) :: astruct
    character(len=1), intent(in) :: geocode
    real(gp), intent(in) :: rxyz(3,nat)
    !internal
    integer :: u
    integer :: iat,ityp

    u=f_get_free_unit()
    open(unit=u,file="input_geometry.gen")
    if(geocode=='F') then
        write(u,'(i5.5,a)') nat, " C"
    else if(geocode=='P') then
        write(u,'(i5.5,a)') nat, " S"
    else if(geocode=='S') then
        call f_err_throw('Surface boundary conditions not available for DFTB+')
    endif


    write(u,*) (trim(adjustl(astruct%atomnames(ityp)))//" ", ityp=1,astruct%ntypes)

    do iat = 1, nat
        write(u,'(i5,1x,i5,3(1x,es25.15))') iat, astruct%iatype(iat), rxyz(:,iat)* Bohr_Ang
    end do
    if(geocode/='F') then
        write(u,'(3(1x,es25.15))') 0.d0,0.d0,0.d0
        write(u,'(3(1x,es25.15))') cell(:, 1)*Bohr_Ang
        write(u,'(3(1x,es25.15))') cell(:, 2)*Bohr_Ang
        write(u,'(3(1x,es25.15))') cell(:, 3)*Bohr_Ang
    endif
    close(u)

    u=f_get_free_unit()
    open(unit=u,file="input_restart.gen")
    if (policy == 10001)then!INPUT_POLICY_MEMORY
        write(u,'(a)')'Yes'
    else
        write(u,'(a)')'No'
    endif
    close(u)

end subroutine

end module module_dftbp
