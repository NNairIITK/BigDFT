!>  @file
!!  forces module
!! @author
!!    Copyright (C) 2015-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


module module_forces

    implicit none

    private

    public :: clean_forces
    public :: clean_forces_base

contains

subroutine clean_forces(iproc,at,rxyz,fxyz,fnoise,run_mode)
    use module_base
    use f_utils                                                        
    use dynamic_memory                                                 
    use public_enums
    use module_atoms!types
    use yaml_output
    !parameters
    implicit none
    integer, intent(in) :: iproc
    type(atoms_data), intent(in) :: at
    real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
    real(gp), dimension(3,at%astruct%nat), intent(inout) :: fxyz
    real(gp), intent(out) :: fnoise
    type(f_enumerator), intent(in),optional :: run_mode
    !local
    logical :: QM_clean
    QM_clean=.true.
    if(present(run_mode))then
        if(trim(char(run_mode))/='QM_RUN_MODE')then
            QM_clean=.false.
        else
            return
        endif
    endif
    if(QM_clean)then
        call clean_forces_dft(iproc,at,rxyz,fxyz,fnoise)
    else
        call clean_forces_base(at,fxyz)
    endif
end subroutine clean_forces


subroutine clean_forces_base(at,fxyz)
  use module_base
  use module_atoms!types
  implicit none
  !parameters
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%astruct%nat), intent(inout) :: fxyz
  !internal
  integer :: iat,ixyz
  integer, dimension(3) :: ijk
  integer :: n_bloc1, n_bloc2                 !< Number of atoms allowed to move only as blocs.
  real(gp), dimension(3) :: f_bloc1, f_bloc2  !< Sum, then average of the forces in blocs.
  real(gp), dimension(3) :: u
  real(gp) :: scal
  !!Clean the forces for blocked atoms
  !!Modification by FL: atom possibly frozen in moving blocs.
  !!@todo Need a better handling of the given constraints
  f_bloc1 = 0.0_gp
  f_bloc2 = 0.0_gp
  n_bloc1 = 0
  n_bloc2 = 0
  do iat=1,at%astruct%nat
     if (at%astruct%ifrztyp(iat) < 1000) then
        if (at%astruct%ifrztyp(iat) < 200) then
           do ixyz=1,3
              if (.not. move_this_coordinate(at%astruct%ifrztyp(iat),ixyz)) fxyz(ixyz,iat)=0.0_gp
           end do
        else
           ! internal coordinates, will be handled separately
        end if
     else if (at%astruct%ifrztyp(iat) == 1001)   then   ! atom "iat" in bloc 1.
       f_bloc1 = f_bloc1 + fxyz(:,iat)
       n_bloc1 = n_bloc1 + 1                            ! could be done once, after reading the inputs.
     else if (at%astruct%ifrztyp(iat) == 1002)   then   ! atom "iat" in bloc 2. Can't be in 2 blocs.
       f_bloc2 = f_bloc2 + fxyz(:,iat)
       n_bloc2 = n_bloc2 + 1  ! could be done once, after reading the inputs.
     else
        ! Projection on a plane, defined by Miller indices stored in ifrztyp:
        !  ifrztyp(iat) = 9ijk
        ijk = (/ (at%astruct%ifrztyp(iat) - 9000) / 100, &
             & modulo(at%astruct%ifrztyp(iat) - 9000, 100) / 10, &
             & modulo(at%astruct%ifrztyp(iat) - 9000, 10) /)
        u = (/ at%astruct%cell_dim(1) / real(ijk(1), gp), &
             & at%astruct%cell_dim(2) / real(ijk(2), gp), &
             & at%astruct%cell_dim(3) / real(ijk(3), gp) /)
        u = u / nrm2(3, u(1), 1)
        scal = fxyz(1,iat) * u(1) + fxyz(2,iat) * u(2) + fxyz(3,iat) * u(3)
        fxyz(1,iat)=fxyz(1,iat) - scal * u(1)
        fxyz(2,iat)=fxyz(2,iat) - scal * u(2)
        fxyz(3,iat)=fxyz(3,iat) - scal * u(3)
     end if
  end do
  !--- We don't do the following in most of the cases ; only when blocs are defined:
  if ( n_bloc1 .ne. 0 )   f_bloc1 = f_bloc1 / n_bloc1
  if ( n_bloc2 .ne. 0 )   f_bloc2 = f_bloc2 / n_bloc2
  if_atoms_in_blocs: &
  if ( n_bloc1 .ne. 0  .or.  n_bloc2 .ne. 0 )   then
    !--- Forces of atoms in blocs are replaced by the average force in the bloc. Then
       ! - by action and reaction principle, internal forces are suppressed;
       ! - all atoms in a bloc have the same force => same displacments;
       ! - gradient of E relative to the bloc center of gravity is -n_bloc*f_bloc.
    do iat=1,at%astruct%nat
      if (at%astruct%ifrztyp(iat) == 1001)   then   ! atom "iat" in bloc 1.
         fxyz(:,iat) = f_bloc1
      else if (at%astruct%ifrztyp(iat) == 1002)   then   ! atom "iat" in bloc 2. Can't be in 2 blocs.
         fxyz(:,iat) = f_bloc2
      end if
    end do
  end if if_atoms_in_blocs
  !--- End of "Modification by FL: atom possibly frozen in moving blocs".
end subroutine clean_forces_base


subroutine clean_forces_dft(iproc,at,rxyz,fxyz,fnoise)
  use module_base
  use module_atoms!types
  use yaml_output
  implicit none
  !Arguments
  integer, intent(in) :: iproc
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  real(gp), dimension(3,at%astruct%nat), intent(inout) :: fxyz
  real(gp), intent(out) :: fnoise
  !local variables
  integer :: iat
  real(gp) :: sumx,sumy,sumz
  !my variables
  real(gp):: fmax1,t1,t2,t3,fnrm1
  real(gp):: fmax2,fnrm2
  !local variables for blocs (FL)


  !The maximum force and force norm is computed prior to modification of the forces
  fmax1=0._gp
  fnrm1=0._gp
  do iat=1,at%astruct%nat
     t1=fxyz(1,iat)**2
     t2=fxyz(2,iat)**2
     t3=fxyz(3,iat)**2
     fmax1=max(fmax1,sqrt(t1+t2+t3))
     fnrm1=fnrm1+t1+t2+t3
  enddo
  
  
  sumx=0.0_gp
  sumy=0.0_gp
  sumz=0.0_gp
  do iat=1,at%astruct%nat
     sumx=sumx+fxyz(1,iat)
     sumy=sumy+fxyz(2,iat)
     sumz=sumz+fxyz(3,iat)
  enddo
  if (at%astruct%nat /= 0) then 
     fnoise=sqrt((sumx**2+sumy**2+sumz**2)/real(at%astruct%nat,gp))
     sumx=sumx/real(at%astruct%nat,gp)
     sumy=sumy/real(at%astruct%nat,gp)
     sumz=sumz/real(at%astruct%nat,gp)
  else
     fnoise = 0.0_gp
  end if

  if (iproc==0) then 
     !write( *,'(1x,a,1x,3(1x,1pe9.2))') &
     !  'Subtracting center-mass shift of',sumx,sumy,sumz
!           write(*,'(1x,a)')'the sum of the forces is'

     call yaml_mapping_open('Average noise forces',flow=.true.)
     call yaml_map('x',sumx*sqrt(real(at%astruct%nat,gp)),fmt='(1pe16.8)')
     call yaml_map('y',sumy*sqrt(real(at%astruct%nat,gp)),fmt='(1pe16.8)')
     call yaml_map('z',sumz*sqrt(real(at%astruct%nat,gp)),fmt='(1pe16.8)')
     call yaml_map('total',sqrt(sumx**2+sumy**2+sumz**2)*sqrt(real(at%astruct%nat,gp)),fmt='(1pe16.8)')
     call yaml_mapping_close()
     !     write(*,'(a,1pe16.8)')' average noise along x direction: ',sumx*sqrt(real(at%astruct%nat,gp))
     !     write(*,'(a,1pe16.8)')' average noise along y direction: ',sumy*sqrt(real(at%astruct%nat,gp))
     !     write(*,'(a,1pe16.8)')' average noise along z direction: ',sumz*sqrt(real(at%astruct%nat,gp))
     !     write(*,'(a,1pe16.8)')' total average noise            : ',sqrt(sumx**2+sumy**2+sumz**2)*sqrt(real(at%astruct%nat,gp))
!!$
!!$     write(*,'(a,1x,1pe24.17)') 'translational force along x=', sumx  
!!$     write(*,'(a,1x,1pe24.17)') 'translational force along y=', sumy  
!!$     write(*,'(a,1x,1pe24.17)') 'translational force along z=', sumz  
  end if
  
  if (at%astruct%geocode == 'F') then
     do iat=1,at%astruct%nat
        fxyz(1,iat)=fxyz(1,iat)-sumx
        fxyz(2,iat)=fxyz(2,iat)-sumy
        fxyz(3,iat)=fxyz(3,iat)-sumz
     enddo
     
     call elim_torque_reza(at%astruct%nat,rxyz,fxyz)
     
  else if (at%astruct%geocode == 'S') then
     do iat=1,at%astruct%nat
        fxyz(2,iat)=fxyz(2,iat)-sumy
     enddo
  end if

  call clean_forces_base(at,fxyz)
  
  !the noise of the forces is the norm of the translational force
!  fnoise=real(at%astruct%nat,gp)**2*(sumx**2+sumy**2+sumz**2)

  !The maximum force and force norm is computed after modification of the forces
  fmax2=0._gp
  fnrm2=0._gp
  do iat=1,at%astruct%nat
     t1=fxyz(1,iat)**2
     t2=fxyz(2,iat)**2
     t3=fxyz(3,iat)**2
     fmax2=max(fmax2,sqrt(t1+t2+t3))
     fnrm2=fnrm2+t1+t2+t3
  enddo

  if (iproc==0) then
     call yaml_mapping_open('Clean forces norm (Ha/Bohr)',flow=.true.)
     call yaml_map('maxval', fmax2,fmt='(1pe20.12)')
     call yaml_map('fnrm2',  fnrm2,fmt='(1pe20.12)')
     call yaml_mapping_close()
     if (at%astruct%geocode /= 'P') then
        call yaml_mapping_open('Raw forces norm (Ha/Bohr)',flow=.true.)
        call yaml_map('maxval', fmax1,fmt='(1pe20.12)')
        call yaml_map('fnrm2',  fnrm1,fmt='(1pe20.12)')
        call yaml_mapping_close()
     end if
     !write(*,'(2(1x,a,1pe20.12))') 'clean forces norm (Ha/Bohr): maxval=', fmax2, ' fnrm2=', fnrm2
     !if (at%astruct%geocode /= 'P') &
     !&  write(*,'(2(1x,a,1pe20.12))') 'raw forces:                  maxval=', fmax1, ' fnrm2=', fnrm1
  end if
END SUBROUTINE clean_forces_dft

end module module_forces
