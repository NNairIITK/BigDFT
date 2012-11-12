!>  @file
!!  Routines to calculate the local part of atomic forces
!! @author
!!    Copyright (C) 2007-2012 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Calculate atomic forces via finite differences (test purpose)
subroutine forces_via_finite_differences(iproc,nproc,atoms,inputs,energy,fxyz,fnoise,rst,infocode)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,nproc
  integer, intent(inout) :: infocode
  real(gp), intent(inout) :: energy,fnoise
  type(input_variables), intent(inout) :: inputs
  type(atoms_data), intent(inout) :: atoms
  type(restart_objects), intent(inout) :: rst
  real(gp), dimension(3,atoms%nat), intent(inout) :: fxyz
  !local variables
  character(len=*), parameter :: subname='forces_via_finite_differences'
  character(len=4) :: cc
  integer :: ik,km,n_order,i_all,i_stat,iat,ii,i,k,order,iorb_ref
  real(gp) :: dd,alat,functional_ref,fd_alpha,energy_ref
  real(gp), dimension(3) :: fd_step
  real(gp), dimension(6) :: strten
  integer, dimension(:), allocatable :: kmoves
  real(gp), dimension(:), allocatable :: functional,dfunctional
  real(gp), dimension(:,:), allocatable :: rxyz_ref,fxyz_fake

  interface
     subroutine cluster(nproc,iproc,atoms,rxyz,energy,fxyz,strten,fnoise,&
          KSwfn,tmb,&!psi,Lzd,gaucoeffs,gbd,orbs,
          rxyz_old,hx_old,hy_old,hz_old,in,GPU,infocode)
       use module_base
       use module_types
       implicit none
       integer, intent(in) :: nproc,iproc
       integer, intent(out) :: infocode
       real(gp), intent(inout) :: hx_old,hy_old,hz_old
       type(input_variables), intent(in) :: in
       !type(local_zone_descriptors), intent(inout) :: Lzd
       type(atoms_data), intent(inout) :: atoms
       !type(gaussian_basis), intent(inout) :: gbd
       !type(orbitals_data), intent(inout) :: orbs
       type(GPU_pointers), intent(inout) :: GPU
       type(DFT_wavefunction), intent(inout) :: KSwfn,tmb
       real(gp), intent(out) :: energy,fnoise
       real(gp), dimension(3,atoms%nat), intent(inout) :: rxyz_old
       real(gp), dimension(3,atoms%nat), target, intent(inout) :: rxyz
       real(gp), dimension(6), intent(out) :: strten
       real(gp), dimension(3,atoms%nat), intent(out) :: fxyz
       !real(wp), dimension(:), pointer :: psi
       !real(wp), dimension(:,:), pointer :: gaucoeffs
     END SUBROUTINE cluster
  end interface

  if (iproc == 0) then
     write(*,*)
     write(*,'(1x,a,59("="))') '=Forces via finite Difference '
  end if

  !read the file (experimental version)
  open(unit=79,file='input.finite_difference_forces',status='unknown')
  read(79,*) order,fd_alpha
  read(79,*) iorb_ref
  close(unit=79)

  !read the step size
  ! Initialize freq_step (step to move atoms)
  fd_step(1) = fd_alpha*inputs%hx
  fd_step(2) = fd_alpha*inputs%hy
  fd_step(3) = fd_alpha*inputs%hz

  !first, mark the reference energy
  energy_ref=energy

  !assign the reference
  functional_ref=functional_definition(iorb_ref,energy)

  if (order == -1) then
     n_order = 1
     allocate(kmoves(n_order+ndebug),stat=i_stat)
     kmoves = (/ -1 /)
  else if (order == 1) then
     n_order = 1
     allocate(kmoves(n_order+ndebug),stat=i_stat)
     kmoves = (/ 1 /)
  else if (order == 2) then
     n_order = 2
     allocate(kmoves(n_order+ndebug),stat=i_stat)
     kmoves = (/ -1, 1 /)
  else if (order == 3) then
     n_order = 4
     allocate(kmoves(n_order+ndebug),stat=i_stat)
     kmoves = (/ -2, -1, 1, 2 /)
  else
     print *, "Finite Differences: This order",order," is not implemented!"
     stop
  end if
  call memocc(i_stat,kmoves,'kmoves',subname)

  allocate(functional(n_order+ndebug),stat=i_stat)
  call memocc(i_stat,functional,'functional',subname)
  allocate(dfunctional(3*atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,dfunctional,'dfunctional',subname)
  allocate(rxyz_ref(3,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,rxyz_ref,'rxyz_ref',subname)
  allocate(fxyz_fake(3,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,fxyz_fake,'fxyz_fake',subname)


  call razero(3*atoms%nat,dfunctional)

  !write reference in the array
  call dcopy(3*atoms%nat,rst%rxyz_new,1,rxyz_ref,1)

  do iat=1,atoms%nat

     if (atoms%ifrztyp(iat) == 1) then
        if (iproc == 0) write(*,"(1x,a,i0,a)") '=F:The atom ',iat,' is frozen.'
        cycle
     end if

     do i=1,3 !a step in each of the three directions
        ii = i+3*(iat-1)
        if (i==1) then
           alat=atoms%alat1
           cc(3:4)='*x'
        else if (i==2) then
           alat=atoms%alat2
           cc(3:4)='*y'
        else
           alat=atoms%alat3
           cc(3:4)='*z'
        end if
        km = 0
        functional=0.0_gp
        do ik=1,n_order
           k = kmoves(ik)
           !-1-> 1, 1 -> 2, y = ( x + 3 ) / 2
           km = km + 1
           write(cc(1:2),"(i2)") k
           !Displacement
           dd=real(k,gp)*fd_step(i)
           !We copy atomic positions (not necessary)
           call dcopy(3*atoms%nat,rxyz_ref,1,rst%rxyz_new,1)
           if (iproc == 0) then
              write(*,"(1x,a,i0,a,a,a,1pe20.10,a)") &
                   '=FD Move the atom ',iat,' in the direction ',cc,' by ',dd,' bohr'
           end if
           if (atoms%geocode == 'P') then
              rst%rxyz_new(i,iat)=modulo(rxyz_ref(i,iat)+dd,alat)
           else if (atoms%geocode == 'S') then
              rst%rxyz_new(i,iat)=modulo(rxyz_ref(i,iat)+dd,alat)
           else
              rst%rxyz_new(i,iat)=rxyz_ref(i,iat)+dd
           end if
           inputs%inputPsiId=1
           !here we should call cluster
           call cluster(nproc,iproc,atoms,rst%rxyz_new,energy,fxyz_fake,strten,fnoise,&
                rst%KSwfn,rst%tmb,&!psi,rst%Lzd,rst%gaucoeffs,rst%gbd,rst%orbs,&
                rst%rxyz_old,rst%hx_old,rst%hy_old,rst%hz_old,inputs,rst%GPU,infocode)

           !assign the quantity which should be differentiated
           functional(km)=functional_definition(iorb_ref,energy)
!!$           if (iorb_ref==0) then
!!$              functional(km)=energy
!!$           else if (iorb_ref==-1) then
!!$              functional(km)=rst%orbs%HLgap
!!$           else if(iorb_ref < -1) then      !definition which brings to the neutral fukui function (chemical potential)
!!$              !definition which brings Chemical potential
!!$              functional(km)=-abs(rst%orbs%eval(-iorb_ref)+ 0.5_gp*rst%orbs%HLgap)
!!$           else
!!$              functional(km)=rst%orbs%eval(iorb_ref)
!!$           end if
           
        end do
        ! Build the finite-difference quantity if the calculatio has converged properly
        if (infocode ==0) then
           !Force is -dE/dR
           if (order == -1) then
              dd = - (functional_ref - functional(1))/fd_step(i)
           else if (order == 1) then
              dd = - (functional(1) - functional_ref)/fd_step(i)
           else if (order == 2) then
              dd = - (functional(2) - functional(1))/(2.0_gp*fd_step(i))
           else if (order == 3) then
              dd = - (functional(4) + functional(3) - functional(2) - functional(1))/(6.d0*fd_step(i))
           else
              stop "BUG (FD_forces): this order is not defined"
           end if
           !if (abs(dd).gt.1.d-10) then
           dfunctional(ii) = dd
           !end if
        else
           if (iproc==0)&
                write(*,*)'ERROR: the wavefunctions have not converged properly, meaningless result. Exiting. Infocode:',infocode
           stop
        end if
        
     end do
  end do

  !copy the final value of the energy and of the dfunctional
  if (.not. experimental_modulebase_var_onlyfion) then !normal case
     call dcopy(3*atoms%nat,dfunctional,1,fxyz,1)
  else
     call axpy(3*atoms%nat,2.0_gp*rst%KSwfn%orbs%norb,dfunctional(1),1,fxyz(1,1),1)
  end if
  !clean the center mass shift and the torque in isolated directions
  call clean_forces(iproc,atoms,rxyz_ref,fxyz,fnoise)

  energy=functional_ref

  if (iproc == 0) then
     write(*,"(1x,2(a,1pe20.10))") &
          '=FD Step done, Internal Energy:',energy_ref,' functional value:', functional_ref
  end if


  i_all=-product(shape(kmoves))*kind(kmoves)
  deallocate(kmoves,stat=i_stat)
  call memocc(i_stat,i_all,'kmoves',subname)
  i_all=-product(shape(functional))*kind(functional)
  deallocate(functional,stat=i_stat)
  call memocc(i_stat,i_all,'functional',subname)
  i_all=-product(shape(dfunctional))*kind(dfunctional)
  deallocate(dfunctional,stat=i_stat)
  call memocc(i_stat,i_all,'dfunctional',subname)
  i_all=-product(shape(rxyz_ref))*kind(rxyz_ref)
  deallocate(rxyz_ref,stat=i_stat)
  call memocc(i_stat,i_all,'rxyz_ref',subname)
  i_all=-product(shape(fxyz_fake))*kind(fxyz_fake)
  deallocate(fxyz_fake,stat=i_stat)
  call memocc(i_stat,i_all,'fxyz_fake',subname)

contains
  
  function functional_definition(iorb_ref,energy)
    use module_base
    use module_types
    implicit none
    integer, intent(in) :: iorb_ref
    real(gp), intent(in) :: energy
    real(gp) :: functional_definition
    !local variables
    real(gp) :: mu

    !chemical potential =1/2(e_HOMO+e_LUMO)= e_HOMO + 1/2 GAP (the sign is to be decided - electronegativity?)
    !definition which brings to Chemical Potential
    if (rst%KSwfn%orbs%HLgap/=UNINITIALIZED(rst%KSwfn%orbs%HLgap) .and. iorb_ref< -1) then
       mu=-abs(rst%KSwfn%orbs%eval(-iorb_ref)+ 0.5_gp*rst%KSwfn%orbs%HLgap) 
    else
       mu=UNINITIALIZED(1.0_gp)
    end if

    !assign the reference
    if (iorb_ref==0) then
       functional_definition=energy
    else if (iorb_ref == -1) then
       if (rst%KSwfn%orbs%HLgap/=UNINITIALIZED(rst%KSwfn%orbs%HLgap)) then
          functional_definition=rst%KSwfn%orbs%HLgap !here we should add the definition which brings to Fukui function
       else
          stop ' ERROR (FDforces): gap not defined' 
       end if
    else if(iorb_ref < -1) then      !definition which brings to the neutral fukui function (chemical potential)
       if (rst%KSwfn%orbs%HLgap/=UNINITIALIZED(rst%KSwfn%orbs%HLgap)) then
          functional_definition=mu!-mu*real(2*orbs%norb,gp)+energy
       else
          stop ' ERROR (FDforces): gap not defined, chemical potential cannot be calculated' 
       end if
    else
       functional_definition=rst%KSwfn%orbs%eval(iorb_ref)
    end if
    
  end function functional_definition

end subroutine forces_via_finite_differences


subroutine calculate_forces(iproc,nproc,psolver_groupsize,Glr,atoms,orbs,nlpspd,rxyz,hx,hy,hz,proj,i3s,n3p,nspin,&
     refill_proj,ngatherarr,rho,pot,potxc,psi,fion,fdisp,fxyz,&
     ewaldstr,hstrten,xcstr,strten,fnoise,pressure,psoffset,imode,tmb,fpulay)
  use module_base
  use module_types
  use module_interfaces, except_this_one => calculate_forces
  use yaml_output
  implicit none
  logical, intent(in) :: refill_proj
  integer, intent(in) :: iproc,nproc,i3s,n3p,nspin,psolver_groupsize,imode
  real(gp), intent(in) :: hx,hy,hz,psoffset
  type(locreg_descriptors), intent(in) :: Glr
  type(atoms_data), intent(in) :: atoms
  type(orbitals_data), intent(in) :: orbs
  type(nonlocal_psp_descriptors), intent(in) :: nlpspd
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
  real(wp), dimension(nlpspd%nprojel), intent(inout) :: proj
  real(wp), dimension(Glr%d%n1i,Glr%d%n2i,n3p), intent(in) :: rho,pot,potxc
  real(wp), dimension(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(in) :: psi
  real(gp), dimension(6), intent(in) :: ewaldstr,hstrten,xcstr
  real(gp), dimension(3,atoms%nat), intent(in) :: rxyz,fion,fdisp
  real(gp), intent(out) :: fnoise,pressure
  real(gp), dimension(6), intent(out) :: strten
  real(gp), dimension(3,atoms%nat), intent(out) :: fxyz
  type(DFT_wavefunction),intent(in),optional :: tmb
  real(gp),dimension(3,atoms%nat),optional,intent(in) :: fpulay
  !local variables
  integer :: ierr,iat,i,j
  real(gp) :: charge,ucvol
  real(gp), dimension(6,4) :: strtens!local,nonlocal,kin,erf
  character(len=16), dimension(4) :: messages

  call to_zero(6,strten(1))

  call to_zero(6*4,strtens(1,1))

  call local_forces(iproc,atoms,rxyz,0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,&
       Glr%d%n1,Glr%d%n2,Glr%d%n3,n3p,i3s,Glr%d%n1i,Glr%d%n2i,rho,pot,fxyz,strtens(1,1),charge)

  !calculate forces originated by rhocore
  call rhocore_forces(iproc,atoms,nspin,Glr%d%n1,Glr%d%n2,Glr%d%n3,Glr%d%n1i,Glr%d%n2i,n3p,i3s,&
       0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,rxyz,potxc,fxyz)

  !for a taksgroup Poisson Solver, multiply by the ratio.
  !it is important that the forces are bitwise identical among the processors.
  if (psolver_groupsize < nproc) call vscal(3*atoms%nat,real(psolver_groupsize,gp)/real(nproc,gp),fxyz(1,1),1)
  
  !if (iproc == 0 .and. verbose > 1) write( *,'(1x,a)',advance='no')'Calculate nonlocal forces...'
 
  if (imode==0) then
      !cubic version of nonlocal forces
  call nonlocal_forces(iproc,Glr,hx,hy,hz,atoms,rxyz,&
       orbs,nlpspd,proj,Glr%wfd,psi,fxyz,refill_proj,strtens(1,2))
  else if (imode==1) then
      !linear version of nonlocal forces
      call nonlocal_forces_linear(iproc,nproc,tmb%lzd%glr,hx,hy,hz,atoms,rxyz,&
           tmb%orbs,nlpspd,proj,tmb%lzd,tmb%psi,tmb%wfnmd%density_kernel,fxyz,refill_proj,strtens(1,2))
  else
      stop 'wrong imode'
  end if

  !if (iproc == 0 .and. verbose > 1) write( *,'(1x,a)')'done.'
  
  if (iproc == 0 .and. verbose > 1) call yaml_map('Non Local forces calculated',.true.)
  
  if (atoms%geocode == 'P' .and. psolver_groupsize == nproc) then
     call local_hamiltonian_stress(orbs,Glr,hx,hy,hz,psi,strtens(1,3))

     call erf_stress(atoms,rxyz,0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,n3p,&
          iproc,nproc,ngatherarr,rho,strtens(1,4)) !shouud not be reduced for the moment
  end if

  ! Add up all the force contributions
  if (nproc > 1) then
     call mpiallred(fxyz(1,1),3*atoms%nat,MPI_SUM,bigdft_mpi%mpi_comm,ierr)
       if (atoms%geocode == 'P') &
            call mpiallred(strtens(1,1),6*3,MPI_SUM,bigdft_mpi%mpi_comm,ierr) !do not reduce erfstr
     call mpiallred(charge,1,MPI_SUM,bigdft_mpi%mpi_comm,ierr)
  end if

  !add to the forces the ionic and dispersion contribution 
  if (.not. experimental_modulebase_var_onlyfion) then !normal case
     do iat=1,atoms%nat
        fxyz(1,iat)=fxyz(1,iat)+fion(1,iat)+fdisp(1,iat)
        fxyz(2,iat)=fxyz(2,iat)+fion(2,iat)+fdisp(2,iat)
        fxyz(3,iat)=fxyz(3,iat)+fion(3,iat)+fdisp(3,iat)
     enddo
  else
     call vcopy(3*atoms%nat,fion(1,1),1,fxyz(1,1),1)
  end if

  ! Pulay correction if present
  if (present(fpulay)) then
      if (iproc==0) then
          !!do iat=1,atoms%nat
          !!    write(444,'(2es16.6)') fxyz(1,iat), fpulay(1,iat)
          !!    write(444,'(2es16.6)') fxyz(2,iat), fpulay(2,iat)
          !!    write(444,'(2es16.6)') fxyz(3,iat), fpulay(3,iat)
          !!end do
          !!write(444,'(a)') '============================================'
      end if
      !!write(*,*) 'WARNING: IGNORE PULAY FORCES'
      do iat=1,atoms%nat
          fxyz(1,iat) = fxyz(1,iat) - fpulay(1,iat)
          fxyz(2,iat) = fxyz(2,iat) - fpulay(2,iat)
          fxyz(3,iat) = fxyz(3,iat) - fpulay(3,iat)
      end do
  end if



  !clean the center mass shift and the torque in isolated directions
  call clean_forces(iproc,atoms,rxyz,fxyz,fnoise)

  !volume element for local stress
  strtens(:,1)=strtens(:,1)/real(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i,dp)
  strtens(1:3,1)=strtens(1:3,1)+charge*psoffset&
       /real(0.5_gp*hx*0.5_gp*hy*0.5_gp*hz,gp)**2.0_gp/real(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i,dp)**2.0_gp

  if (atoms%geocode == 'P') then
     ucvol=atoms%alat1*atoms%alat2*atoms%alat3 !orthorombic cell
     if (iproc==0) call yaml_open_map('Stress Tensor')
     !sum and symmetrize results
     if (iproc==0 .and. verbose > 2)call write_strten_info(.false.,ewaldstr,ucvol,pressure,'Ewald')
     if (iproc==0 .and. verbose > 2)call write_strten_info(.false.,hstrten,ucvol,pressure,'Hartree')
     if (iproc==0 .and. verbose > 2)call write_strten_info(.false.,xcstr,ucvol,pressure,'XC')
     do j=1,6
        strten(j)=ewaldstr(j)+hstrten(j)+xcstr(j)
     end do
     messages(1)='PSP Short Range'
     messages(2)='PSP Projectors'
     messages(3)='Kinetic'
     messages(4)='PSP Long Range'
     !here we should add the pretty printings
     do i=1,4
        if (atoms%sym%symObj >= 0) call symm_stress((iproc==0 .and. i==1),strtens(1,i),atoms%sym%symObj)
        if (iproc==0 .and. verbose>2)&
             call write_strten_info(.false.,strtens(1,i),ucvol,pressure,trim(messages(i)))
        do j=1,6
           strten(j)=strten(j)+strtens(j,i)
        end do
     end do
     !final result
     pressure=(strten(1)+strten(2)+strten(3))/3.0_gp
     if (iproc==0)call write_strten_info(.true.,strten,ucvol,pressure,'Total')
     if (iproc==0) call yaml_close_map()
  end if

!!$  if (iproc == 0) then
!!$     sumx=0.d0 ; sumy=0.d0 ; sumz=0.d0
!!$     fumx=0.d0 ; fumy=0.d0 ; fumz=0.d0
!!$     do iat=1,atoms%nat
!!$        sumx=sumx+fxyz(1,iat) ; sumy=sumy+fxyz(2,iat) ; sumz=sumz+fxyz(3,iat)
!!$        fumx=fumx+fion(1,iat) ; fumy=fumy+fion(2,iat) ; fumz=fumz+fion(3,iat)
!!$     enddo
!!$     write(77,'(a30,3(1x,e10.3))') 'translat. force total pot ',sumx,sumy,sumz
!!$     write(77,'(a30,3(1x,e10.3))') 'translat. force ionic pot ',fumx,fumy,fumz
!!$  endif

  ! Apply symmetries when needed
  if (atoms%sym%symObj >= 0) call symmetrise_forces(iproc,fxyz,atoms)
end subroutine calculate_forces

!> calculate the contribution to the forces given by the core density charge
subroutine rhocore_forces(iproc,atoms,nspin,n1,n2,n3,n1i,n2i,n3p,i3s,hxh,hyh,hzh,rxyz,potxc,fxyz)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,n1i,n2i,n3p,i3s,nspin,n1,n2,n3
  real(gp), intent(in) :: hxh,hyh,hzh
  type(atoms_data), intent(in) :: atoms
  real(wp), dimension(n1i*n2i*n3p*nspin), intent(in) :: potxc
  real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
  real(gp), dimension(3,atoms%nat), intent(inout) :: fxyz
  !local variables
  real(gp), parameter :: oneo4pi=.079577471545947_wp
  logical :: perx,pery,perz,gox,goy,goz
  integer :: ispin,ilcc,ityp,iat,jtyp,islcc,ngv,ngc,ig,ispinsh,ind
  integer :: nbl1,nbl2,nbl3,nbr1,nbr2,nbr3,isx,isy,isz,iex,iey,iez
  integer :: i1,i2,i3,j1,j2,j3
  real(gp) :: spinfac,rx,ry,rz,frcx,frcy,frcz,rloc,cutoff,x,y,z,r2
  real(gp) :: spherical_gaussian_value,drhoc,drhov,drhodr2

  if (atoms%donlcc) then
     if (iproc == 0) write(*,'(1x,a)',advance='no')'Calculate NLCC forces...'

     if (nspin==1) then
        spinfac=2.0_gp
     else if (nspin ==2) then
        spinfac=1.0_gp
     end if
     !perform the loop on any of the atoms which have this feature
     do iat=1,atoms%nat
        rx=rxyz(1,iat) 
        ry=rxyz(2,iat)
        rz=rxyz(3,iat)

        ityp=atoms%iatype(iat)
        frcx=0.0_gp
        frcy=0.0_gp
        frcz=0.0_gp
        if (atoms%nlcc_ngv(ityp)/=UNINITIALIZED(1) .or. atoms%nlcc_ngc(ityp)/=UNINITIALIZED(1) ) then

           !find the correct position of the nlcc parameters
           ilcc=0
           do jtyp=1,ityp-1
              ngv=atoms%nlcc_ngv(jtyp)
              if (ngv /= UNINITIALIZED(ngv)) ilcc=ilcc+(ngv*(ngv+1)/2)
              ngc=atoms%nlcc_ngc(jtyp)
              if (ngc /= UNINITIALIZED(ngc)) ilcc=ilcc+(ngc*(ngc+1))/2
           end do
           islcc=ilcc

           !find the maximum exponent of the core density
           ngv=atoms%nlcc_ngv(ityp)
           if (ngv==UNINITIALIZED(1)) ngv=0
           ngc=atoms%nlcc_ngc(ityp)
           if (ngc==UNINITIALIZED(1)) ngc=0
           rloc=0.0_gp
           do ig=1,(ngv*(ngv+1))/2+(ngc*(ngc+1))/2
              ilcc=ilcc+1
              rloc=max(rloc,atoms%nlccpar(0,ilcc))
           end do

           cutoff=10.d0*rloc

           !conditions for periodicity in the three directions
           perx=(atoms%geocode /= 'F')
           pery=(atoms%geocode == 'P')
           perz=(atoms%geocode /= 'F')

           call ext_buffers(perx,nbl1,nbr1)
           call ext_buffers(pery,nbl2,nbr2)
           call ext_buffers(perz,nbl3,nbr3)

           if (n3p >0) then

              isx=floor((rx-cutoff)/hxh)
              isy=floor((ry-cutoff)/hyh)
              isz=floor((rz-cutoff)/hzh)

              iex=ceiling((rx+cutoff)/hxh)
              iey=ceiling((ry+cutoff)/hyh)
              iez=ceiling((rz+cutoff)/hzh)
              do ispin=1,nspin
                 ispinsh=0
                 if (ispin==2) ispinsh=n1i*n2i*n3p
                 do i3=isz,iez
                    z=real(i3,kind=8)*hzh-rz
                    call ind_positions(perz,i3,n3,j3,goz)
                    j3=j3+nbl3+1
                    if (j3 >= i3s .and. j3 <= i3s+n3p-1) then
                       do i2=isy,iey
                          y=real(i2,kind=8)*hyh-ry
                          call ind_positions(pery,i2,n2,j2,goy)
                          if (goy) then
                             do i1=isx,iex
                                x=real(i1,kind=8)*hxh-rx
                                call ind_positions(perx,i1,n1,j1,gox)
                                if (gox) then
                                   r2=x**2+y**2+z**2
                                   ilcc=islcc
                                   drhov=0.0_dp
                                   do ig=1,(ngv*(ngv+1))/2
                                      ilcc=ilcc+1
                                      !derivative wrt r2
                                      drhov=drhov+&
                                           spherical_gaussian_value(r2,atoms%nlccpar(0,ilcc),atoms%nlccpar(1,ilcc),1)
                                   end do
                                   drhoc=0.0_dp
                                   do ig=1,(ngc*(ngc+1))/2
                                      ilcc=ilcc+1
                                      !derivative wrt r2
                                      drhoc=drhoc+&
                                           spherical_gaussian_value(r2,atoms%nlccpar(0,ilcc),atoms%nlccpar(1,ilcc),1)
                                   end do
                                   !forces in all the directions for the given atom
                                   ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s+1-1)*n1i*n2i+ispinsh
                                   drhodr2=drhoc-drhov
                                   frcx=frcx+potxc(ind)*x*drhodr2
                                   frcy=frcy+potxc(ind)*y*drhodr2
                                   frcz=frcz+potxc(ind)*z*drhodr2
                                endif
                             enddo
                          end if
                       enddo
                    end if
                 enddo
              end do
           end if
        end if
        !assign contribution per atom
        fxyz(1,iat)=fxyz(1,iat)+frcx*hxh*hyh*hzh*spinfac*oneo4pi
        fxyz(2,iat)=fxyz(2,iat)+frcy*hxh*hyh*hzh*spinfac*oneo4pi
        fxyz(3,iat)=fxyz(3,iat)+frcz*hxh*hyh*hzh*spinfac*oneo4pi
        !print *,'iat,iproc',iat,iproc,frcx*hxh*hyh*hzh*spinfac*oneo4pi
     end do

     if (iproc == 0 .and. verbose > 1) write( *,'(1x,a)')'done.'
  end if
end subroutine rhocore_forces


!> Calculates the local forces acting on the atoms belonging to iproc
subroutine local_forces(iproc,at,rxyz,hxh,hyh,hzh,&
     n1,n2,n3,n3pi,i3s,n1i,n2i,rho,pot,floc,locstrten,charge)
  use module_base
  use module_types
  implicit none
  !Arguments---------
  type(atoms_data), intent(in) :: at
  integer, intent(in) :: iproc,n1,n2,n3,n3pi,i3s,n1i,n2i
  real(gp), intent(in) :: hxh,hyh,hzh 
  real(gp),intent(out) :: charge
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(dp), dimension(*), intent(in) :: rho,pot
  real(gp), dimension(3,at%nat), intent(out) :: floc
  real(gp), dimension(6), intent(out) :: locstrten
  !Local variables---------
  logical :: perx,pery,perz,gox,goy,goz
  real(kind=8) :: pi,prefactor,cutoff,rloc,Vel,rhoel
  real(kind=8) :: fxerf,fyerf,fzerf,fxion,fyion,fzion,fxgau,fygau,fzgau,forceleaked,forceloc
  real(kind=8) :: rx,ry,rz,x,y,z,arg,r2,xp,tt,Txx,Tyy,Tzz,Txy,Txz,Tyz
  integer :: i1,i2,i3,ind,iat,ityp,nloc,iloc
  integer :: nbl1,nbr1,nbl2,nbr2,nbl3,nbr3,j1,j2,j3,isx,isy,isz,iex,iey,iez
  !array of coefficients of the derivative
  real(kind=8), dimension(4) :: cprime 
  
  pi=4.d0*atan(1.d0)

  locstrten=0.0_gp

charge=0.d0
do i3=1,n3pi
        do i2=1,n2i
           do i1=1,n1i
              ind=i1+(i2-1)*n1i+(i3-1)*n1i*n2i
charge=charge+rho(ind)
           enddo
        enddo
     enddo
charge=charge*hxh*hyh*hzh

  if (iproc == 0 .and. verbose > 1) write(*,'(1x,a)',advance='no')'Calculate local forces...'
  forceleaked=0.d0

  !conditions for periodicity in the three directions
  perx=(at%geocode /= 'F')
  pery=(at%geocode == 'P')
  perz=(at%geocode /= 'F')

  call ext_buffers(perx,nbl1,nbr1)
  call ext_buffers(pery,nbl2,nbr2)
  call ext_buffers(perz,nbl3,nbr3)

  do iat=1,at%nat
     ityp=at%iatype(iat)
     !coordinates of the center
     rx=rxyz(1,iat) 
     ry=rxyz(2,iat) 
     rz=rxyz(3,iat)
     !inizialization of the forces
     !ion-ion term
     fxion=0.d0
     fyion=0.d0
     fzion=0.d0
     !ion-electron term, error function part
     fxerf=0.d0
     fyerf=0.d0
     fzerf=0.d0
     !ion-electron term, gaussian part
     fxgau=0.d0
     fygau=0.d0
     fzgau=0.d0
     !local stress tensor component for this atom
     Txx=0.0_gp
     Tyy=0.0_gp
     Tzz=0.0_gp
     Txy=0.0_gp
     Txz=0.0_gp
     Tyz=0.0_gp

     !building array of coefficients of the derivative of the gaussian part
     cprime(1)=2.d0*at%psppar(0,2,ityp)-at%psppar(0,1,ityp)
     cprime(2)=4.d0*at%psppar(0,3,ityp)-at%psppar(0,2,ityp)
     cprime(3)=6.d0*at%psppar(0,4,ityp)-at%psppar(0,3,ityp)
     cprime(4)=-at%psppar(0,4,ityp)

     ! determine number of local terms
     nloc=0
     do iloc=1,4
        if (at%psppar(0,iloc,ityp) /= 0.d0) nloc=iloc
     enddo

     !local part
     rloc=at%psppar(0,0,ityp)
     prefactor=real(at%nelpsp(ityp),kind=8)/(2.d0*pi*sqrt(2.d0*pi)*rloc**5)
     !maximum extension of the gaussian
     cutoff=10.d0*rloc

     isx=floor((rx-cutoff)/hxh)
     isy=floor((ry-cutoff)/hyh)
     isz=floor((rz-cutoff)/hzh)
     
     iex=ceiling((rx+cutoff)/hxh)
     iey=ceiling((ry+cutoff)/hyh)
     iez=ceiling((rz+cutoff)/hzh)

     !calculate the forces near the atom due to the error function part of the potential
     !calculate forces for all atoms only in the distributed part of the simulation box
     if (n3pi >0 ) then
        do i3=isz,iez
           z=real(i3,kind=8)*hzh-rz
           call ind_positions(perz,i3,n3,j3,goz) 
           j3=j3+nbl3+1
           do i2=isy,iey
              y=real(i2,kind=8)*hyh-ry
              call ind_positions(pery,i2,n2,j2,goy)
              do i1=isx,iex
                 x=real(i1,kind=8)*hxh-rx
                 call ind_positions(perx,i1,n1,j1,gox)
                 r2=x**2+y**2+z**2
                 arg=r2/rloc**2
                 xp=exp(-.5d0*arg)
                 if (j3 >= i3s .and. j3 <= i3s+n3pi-1  .and. goy  .and. gox ) then
                    ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s+1-1)*n1i*n2i
                    !gaussian part
                    tt=0.d0
                    if (nloc /= 0) then
                       !derivative of the polynomial
                       tt=cprime(nloc)
                       do iloc=nloc-1,1,-1
                          tt=arg*tt+cprime(iloc)
                       enddo
                       rhoel=rho(ind)
                       forceloc=xp*tt*rhoel
                       fxgau=fxgau+forceloc*x
                       fygau=fygau+forceloc*y
                       fzgau=fzgau+forceloc*z
                       if (r2 /= 0.0_gp) then
                          Txx=Txx+forceloc*x*x
                          Tyy=Tyy+forceloc*y*y
                          Tzz=Tzz+forceloc*z*z
                          Txy=Txy+forceloc*x*y
                          Txz=Txz+forceloc*x*z
                          Tyz=Tyz+forceloc*y*z
                       end if
                    end if
                    !error function part
                    Vel=pot(ind)
                    fxerf=fxerf+xp*Vel*x
                    fyerf=fyerf+xp*Vel*y
                    fzerf=fzerf+xp*Vel*z
                 else if (.not. goz) then
                    !derivative of the polynomial
                    tt=cprime(nloc)
                    do iloc=nloc-1,1,-1
                       tt=arg*tt+cprime(iloc)
                    enddo
                    forceleaked=forceleaked+prefactor*xp*tt*rho(1) !(as a sample value)
                 endif
              end do
           end do
        end do
     end if

     !final result of the forces

     floc(1,iat)=fxion+(hxh*hyh*hzh*prefactor)*fxerf+(hxh*hyh*hzh/rloc**2)*fxgau
     floc(2,iat)=fyion+(hxh*hyh*hzh*prefactor)*fyerf+(hxh*hyh*hzh/rloc**2)*fygau
     floc(3,iat)=fzion+(hxh*hyh*hzh*prefactor)*fzerf+(hxh*hyh*hzh/rloc**2)*fzgau

     locstrten(1)=locstrten(1)+Txx/rloc/rloc
     locstrten(2)=locstrten(2)+Tyy/rloc/rloc
     locstrten(3)=locstrten(3)+Tzz/rloc/rloc
     locstrten(4)=locstrten(4)+Tyz/rloc/rloc
     locstrten(5)=locstrten(5)+Txz/rloc/rloc
     locstrten(6)=locstrten(6)+Txy/rloc/rloc
!!!     !only for testing purposes, printing the components of the forces for each atoms
!!!     write(10+iat,'(2(1x,3(1x,1pe12.5)))') &
!!!          (hxh*hyh*hzh*prefactor)*fxerf,(hxh*hyh*hzh*prefactor)*fyerf,&
!!!          (hxh*hyh*hzh*prefactor)*fzerf,(hxh*hyh*hzh/rloc**2)*fxgau,(hxh*hyh*hzh/rloc**2)*fygau,(hxh*hyh*hzh/rloc**2)*fzgau

  end do !iat

!write(*,*) 'iproc,charge:',iproc,charge

!locstrten(1:3)=locstrten(1:3)+charge*psoffset/(hxh*hyh*hzh)/real(n1i*n2i*n3pi,kind=8)

  forceleaked=forceleaked*hxh*hyh*hzh
  if (iproc == 0 .and. verbose > 1) write(*,'(a,1pe12.5)') 'done. Leaked force: ',forceleaked

END SUBROUTINE local_forces


!> Calculates the nonlocal forces on all atoms arising from the wavefunctions 
!! belonging to iproc and adds them to the force array
!! recalculate the projectors at the end if refill flag is .true.
subroutine nonlocal_forces(iproc,lr,hx,hy,hz,at,rxyz,&
     orbs,nlpspd,proj,wfd,psi,fsep,refill,strten)
  use module_base
  use module_types
  implicit none
  !Arguments-------------
  type(atoms_data), intent(in) :: at
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(nonlocal_psp_descriptors), intent(in) :: nlpspd
  logical, intent(in) :: refill
  integer, intent(in) :: iproc
  real(gp), intent(in) :: hx,hy,hz
  type(locreg_descriptors) :: lr
  type(orbitals_data), intent(in) :: orbs
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(wp), dimension((wfd%nvctr_c+7*wfd%nvctr_f)*orbs%norbp*orbs%nspinor), intent(inout) :: psi
  real(wp), dimension(nlpspd%nprojel), intent(inout) :: proj
  real(gp), dimension(3,at%nat), intent(inout) :: fsep
  real(gp), dimension(6), intent(out) :: strten
  !local variables--------------
  character(len=*), parameter :: subname='nonlocal_forces'
  integer :: istart_c,iproj,iat,ityp,i,j,l,m
  integer :: mbseg_c,mbseg_f,jseg_c,jseg_f
  integer :: mbvctr_c,mbvctr_f,iorb,nwarnings,nspinor,ispinor,jorbd
  real(gp) :: offdiagcoeff,hij,sp0,spi,sp0i,sp0j,spj,orbfac,strc,Enl,vol
  integer :: idir,i_all,i_stat,ncplx,icplx,isorb,ikpt,ieorb,istart_ck,ispsi_k,ispsi,jorb
  real(gp), dimension(2,2,3) :: offdiagarr
  real(gp), dimension(:,:), allocatable :: fxyz_orb
  real(dp), dimension(:,:,:,:,:,:,:), allocatable :: scalprod
  real(gp), dimension(6) :: sab

  call to_zero(6,strten(1)) 

  !quick return if no orbitals on this processor
  if (orbs%norbp == 0) return
     
  !always put complex scalprod
  !also nspinor for the moment is the biggest as possible

  !  allocate(scalprod(2,0:3,7,3,4,at%nat,orbs%norbp*orbs%nspinor+ndebug),stat=i_stat)
  ! need more components in scalprod to calculate terms like dp/dx*psi*x
  allocate(scalprod(2,0:9,7,3,4,at%nat,orbs%norbp*orbs%nspinor+ndebug),stat=i_stat)
  call memocc(i_stat,scalprod,'scalprod',subname)
  call razero(2*10*7*3*4*at%nat*orbs%norbp*orbs%nspinor,scalprod)


  Enl=0._gp
  !strten=0.d0
  vol=real(at%alat1*at%alat2*at%alat3,gp)
  sab=0.d0

  !calculate the coefficients for the off-diagonal terms
  do l=1,3
     do i=1,2
        do j=i+1,3
           offdiagcoeff=0.0_gp
           if (l==1) then
              if (i==1) then
                 if (j==2) offdiagcoeff=-0.5_gp*sqrt(3._gp/5._gp)
                 if (j==3) offdiagcoeff=0.5_gp*sqrt(5._gp/21._gp)
              else
                 offdiagcoeff=-0.5_gp*sqrt(100._gp/63._gp)
              end if
           else if (l==2) then
              if (i==1) then
                 if (j==2) offdiagcoeff=-0.5_gp*sqrt(5._gp/7._gp)
                 if (j==3) offdiagcoeff=1._gp/6._gp*sqrt(35._gp/11._gp)
              else
                 offdiagcoeff=-7._gp/3._gp*sqrt(1._gp/11._gp)
              end if
           else if (l==3) then
              if (i==1) then
                 if (j==2) offdiagcoeff=-0.5_gp*sqrt(7._gp/9._gp)
                 if (j==3) offdiagcoeff=0.5_gp*sqrt(63._gp/143._gp)
              else
                 offdiagcoeff=-9._gp*sqrt(1._gp/143._gp)
              end if
           end if
           offdiagarr(i,j-i,l)=offdiagcoeff
        end do
     end do
  end do

  !look for the strategy of projectors application
  if (DistProjApply) then
     !apply the projectors on the fly for each k-point of the processor
     !starting k-point
     ikpt=orbs%iokpt(1)
     ispsi_k=1
     jorb=0
     loop_kptD: do

        call orbs_in_kpt(ikpt,orbs,isorb,ieorb,nspinor)

        call ncplx_kpt(ikpt,orbs,ncplx)

        nwarnings=0 !not used, simply initialised 
        iproj=0 !should be equal to four times nproj at the end
        jorbd=jorb
        do iat=1,at%nat

           call plr_segs_and_vctrs(nlpspd%plr(iat),&
                mbseg_c,mbseg_f,mbvctr_c,mbvctr_f)
           jseg_c=1
           jseg_f=1

           do idir=0,9
!!$           mbseg_c=nlpspd%nseg_p(2*iat-1)-nlpspd%nseg_p(2*iat-2)
!!$           mbseg_f=nlpspd%nseg_p(2*iat  )-nlpspd%nseg_p(2*iat-1)
!!$           jseg_c=nlpspd%nseg_p(2*iat-2)+1
!!$           jseg_f=nlpspd%nseg_p(2*iat-1)+1
!!$           mbvctr_c=nlpspd%nvctr_p(2*iat-1)-nlpspd%nvctr_p(2*iat-2)
!!$           mbvctr_f=nlpspd%nvctr_p(2*iat  )-nlpspd%nvctr_p(2*iat-1)

           ityp=at%iatype(iat)
              !calculate projectors
              istart_c=1
              call atom_projector(ikpt,iat,idir,istart_c,iproj,nlpspd%nprojel,&
                   lr,hx,hy,hz,rxyz(1,iat),at,orbs,nlpspd%plr(iat),&
                   proj,nwarnings)
!              print *,'iat,ilr,idir,sum(proj)',iat,ilr,idir,sum(proj)
 
              !calculate the contribution for each orbital
              !here the nspinor contribution should be adjusted
              ! loop over all my orbitals
              ispsi=ispsi_k
              jorb=jorbd
              do iorb=isorb,ieorb
                 do ispinor=1,nspinor,ncplx
                    jorb=jorb+1
                    istart_c=1
                    do l=1,4
                       do i=1,3
                          if (at%psppar(l,i,ityp) /= 0.0_gp) then
                             do m=1,2*l-1
                                call wpdot_wrap(ncplx,&
                                     wfd%nvctr_c,wfd%nvctr_f,wfd%nseg_c,wfd%nseg_f,&
                                     wfd%keyvglob,wfd%keyglob,psi(ispsi),&
                                     mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
!!$                                     nlpspd%keyv_p(jseg_c),&
!!$                                     nlpspd%keyg_p(1,jseg_c),&
                                     nlpspd%plr(iat)%wfd%keyvglob(jseg_c),&
                                     nlpspd%plr(iat)%wfd%keyglob(1,jseg_c),&
                                     proj(istart_c),&
                                     scalprod(1,idir,m,i,l,iat,jorb))
                                istart_c=istart_c+(mbvctr_c+7*mbvctr_f)*ncplx
                             end do
                          end if
                       end do
                    end do
                    ispsi=ispsi+(wfd%nvctr_c+7*wfd%nvctr_f)*ncplx
                 end do
              end do
              if (istart_c-1  > nlpspd%nprojel) stop '2:applyprojectors'
           end do

        end do

        if (ieorb == orbs%norbp) exit loop_kptD
        ikpt=ikpt+1
        ispsi_k=ispsi
     end do loop_kptD

  else
     !calculate all the scalar products for each direction and each orbitals
     do idir=0,9

        if (idir /= 0) then !for the first run the projectors are already allocated
           call fill_projectors(iproc,lr,hx,hy,hz,at,orbs,rxyz,nlpspd,proj,idir)
        end if
        !apply the projectors  k-point of the processor
        !starting k-point
        ikpt=orbs%iokpt(1)
        istart_ck=1
        ispsi_k=1
        jorb=0
        loop_kpt: do

           call orbs_in_kpt(ikpt,orbs,isorb,ieorb,nspinor)

           call ncplx_kpt(ikpt,orbs,ncplx)

           ! calculate the scalar product for all the orbitals
           ispsi=ispsi_k
           do iorb=isorb,ieorb
              do ispinor=1,nspinor,ncplx
                 jorb=jorb+1
                 ! loop over all projectors of this k-point
                 iproj=0
                 istart_c=istart_ck
                 do iat=1,at%nat
                    call plr_segs_and_vctrs(nlpspd%plr(iat),&
                         mbseg_c,mbseg_f,mbvctr_c,mbvctr_f)
                    jseg_c=1
                    jseg_f=1

!!$                    mbseg_c=nlpspd%nseg_p(2*iat-1)-nlpspd%nseg_p(2*iat-2)
!!$                    mbseg_f=nlpspd%nseg_p(2*iat  )-nlpspd%nseg_p(2*iat-1)
!!$                    jseg_c=nlpspd%nseg_p(2*iat-2)+1
!!$                    jseg_f=nlpspd%nseg_p(2*iat-1)+1
!!$                    mbvctr_c=nlpspd%nvctr_p(2*iat-1)-nlpspd%nvctr_p(2*iat-2)
!!$                    mbvctr_f=nlpspd%nvctr_p(2*iat  )-nlpspd%nvctr_p(2*iat-1)
                    ityp=at%iatype(iat)
                    do l=1,4
                       do i=1,3
                          if (at%psppar(l,i,ityp) /= 0.0_gp) then
                             do m=1,2*l-1
                                iproj=iproj+1
                                call wpdot_wrap(ncplx,&
                                     wfd%nvctr_c,wfd%nvctr_f,wfd%nseg_c,wfd%nseg_f,&
                                     wfd%keyvglob,wfd%keyglob,psi(ispsi),  &
                                     mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
!!$                                     nlpspd%keyv_p(jseg_c),&
!!$                                     nlpspd%keyg_p(1,jseg_c),&
                                     nlpspd%plr(iat)%wfd%keyvglob(jseg_c),&
                                     nlpspd%plr(iat)%wfd%keyglob(1,jseg_c),&
                                     proj(istart_c),scalprod(1,idir,m,i,l,iat,jorb))
                                istart_c=istart_c+(mbvctr_c+7*mbvctr_f)*ncplx
                             end do
                          end if
                       end do
                    end do
                 end do
                 ispsi=ispsi+(wfd%nvctr_c+7*wfd%nvctr_f)*ncplx
              end do
              if (iproj /= nlpspd%nproj) stop '1:applyprojectors'
           end do
           istart_ck=istart_c
           if (ieorb == orbs%norbp) exit loop_kpt
           ikpt=ikpt+1
           ispsi_k=ispsi
        end do loop_kpt
        if (istart_ck-1  /= nlpspd%nprojel) stop '2:applyprojectors'

     end do

     !restore the projectors in the proj array (for on the run forces calc., tails or so)
     if (refill) then 
        call fill_projectors(iproc,lr,hx,hy,hz,at,orbs,rxyz,nlpspd,proj,0)
     end if

  end if

  allocate(fxyz_orb(3,at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,fxyz_orb,'fxyz_orb',subname)

  !apply the projectors  k-point of the processor
  !starting k-point
  ikpt=orbs%iokpt(1)
  jorb=0
  loop_kptF: do

     call orbs_in_kpt(ikpt,orbs,isorb,ieorb,nspinor)

     call ncplx_kpt(ikpt,orbs,ncplx)

     ! loop over all my orbitals for calculating forces
     do iorb=isorb,ieorb
        sab=0.0_gp
        ! loop over all projectors
        call to_zero(3*at%nat,fxyz_orb(1,1))
        do ispinor=1,nspinor,ncplx
           jorb=jorb+1
           do iat=1,at%nat
              ityp=at%iatype(iat)
              do l=1,4
                 do i=1,3
                    if (at%psppar(l,i,ityp) /= 0.0_gp) then
                       do m=1,2*l-1
                          do icplx=1,ncplx
                             ! scalar product with the derivatives in all the directions
                             sp0=real(scalprod(icplx,0,m,i,l,iat,jorb),gp)
                             !write(200+iproc,'(a,9i6,es18.8)') 'iorb,jorb,icplx,0,m,i,l,iat,iiat,sp0', &
                             !                                   iorb,jorb,icplx,0,m,i,l,iat,iat,sp0
                             do idir=1,3
                                spi=real(scalprod(icplx,idir,m,i,l,iat,jorb),gp)
                                !write(210+iproc,'(a,10i6,es18.8)') 'iorb,jorb,icplx,0,m,i,l,iat,iiat,&
                                !                                    &idir,fxyz_orb(idir,iat)', &
                                !                                    iorb,jorb,icplx,0,m,i,l,iat,iat,&
                                !                                    idir,fxyz_orb(idir,iat)
                                fxyz_orb(idir,iat)=fxyz_orb(idir,iat)+&
                                     at%psppar(l,i,ityp)*sp0*spi
                             end do

Enl=Enl+sp0*sp0*at%psppar(l,i,ityp)*&
orbs%occup(iorb+orbs%isorb)*orbs%kwgts(orbs%iokpt(iorb))
                            do idir=4,9 !for stress
strc=real(scalprod(icplx,idir,m,i,l,iat,jorb),gp)
sab(idir-3)=&
sab(idir-3)+&   
at%psppar(l,i,ityp)*sp0*2.0_gp*strc*&
orbs%occup(iorb+orbs%isorb)*orbs%kwgts(orbs%iokpt(iorb))
                            end do
                          end do
                       end do
                    end if
                 end do
              end do
              !HGH case, offdiagonal terms
              if (at%npspcode(ityp) == 3 .or. at%npspcode(ityp) == 10) then
                 do l=1,3 !no offdiagoanl terms for l=4 in HGH-K case
                    do i=1,2
                       if (at%psppar(l,i,ityp) /= 0.0_gp) then 
                          loop_j: do j=i+1,3
                             if (at%psppar(l,j,ityp) == 0.0_gp) exit loop_j
                             !offdiagonal HGH term
                             if (at%npspcode(ityp) == 3) then !traditional HGH convention
                                hij=offdiagarr(i,j-i,l)*at%psppar(l,j,ityp)
                             else !HGH-K convention
                                hij=at%psppar(l,i+j+1,ityp)
                             end if
                             do m=1,2*l-1
                                !F_t= 2.0*h_ij (<D_tp_i|psi><psi|p_j>+<p_i|psi><psi|D_tp_j>)
                                !(the two factor is below)
                                do icplx=1,ncplx
                                   sp0i=real(scalprod(icplx,0,m,i,l,iat,jorb),gp)
                                   sp0j=real(scalprod(icplx,0,m,j,l,iat,jorb),gp)
                                   do idir=1,3
                                      spi=real(scalprod(icplx,idir,m,i,l,iat,jorb),gp)
                                      spj=real(scalprod(icplx,idir,m,j,l,iat,jorb),gp)
                                      fxyz_orb(idir,iat)=fxyz_orb(idir,iat)+&
                                           hij*(sp0j*spi+spj*sp0i)
                                   end do

Enl=Enl+2.0_gp*sp0i*sp0j*hij&
*orbs%occup(iorb+orbs%isorb)*orbs%kwgts(orbs%iokpt(iorb))
                                  do idir=4,9
spi=real(scalprod(icplx,idir,m,i,l,iat,jorb),gp)
spj=real(scalprod(icplx,idir,m,j,l,iat,jorb),gp)
sab(idir-3)=&
sab(idir-3)+&   
2.0_gp*hij*(sp0j*spi+sp0i*spj)&
*orbs%occup(iorb+orbs%isorb)*orbs%kwgts(orbs%iokpt(iorb))
                                  end do
                                end do
                             end do
                          end do loop_j
                       end if
                    end do
                 end do
              end if
           end do
        end do

        !orbital-dependent factor for the forces
        orbfac=orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*2.0_gp

!seq: strten(1:6) =  11 22 33 23 13 12 
strten(1)=strten(1)+sab(1)/vol 
strten(2)=strten(2)+sab(2)/vol 
strten(3)=strten(3)+sab(3)/vol 
strten(4)=strten(4)+sab(5)/vol
strten(5)=strten(5)+sab(6)/vol
strten(6)=strten(6)+sab(4)/vol
        do iat=1,at%nat
           fsep(1,iat)=fsep(1,iat)+orbfac*fxyz_orb(1,iat)
           fsep(2,iat)=fsep(2,iat)+orbfac*fxyz_orb(2,iat)
           fsep(3,iat)=fsep(3,iat)+orbfac*fxyz_orb(3,iat)
        end do

     end do
     if (ieorb == orbs%norbp) exit loop_kptF
     ikpt=ikpt+1
     ispsi_k=ispsi
  end do loop_kptF


!Adding Enl to the diagonal components of strten after loop over kpts is finished...
do i=1,3
strten(i)=strten(i)+Enl/vol
end do

!!!  do iat=1,at%nat
!!!     write(20+iat,'(1x,i5,1x,3(1x,1pe12.5))') &
!!!          iat,fsep(1,iat),fsep(2,iat),fsep(3,iat)
!!!  end do

  i_all=-product(shape(fxyz_orb))*kind(fxyz_orb)
  deallocate(fxyz_orb,stat=i_stat)
  call memocc(i_stat,i_all,'fxyz_orb',subname)
  i_all=-product(shape(scalprod))*kind(scalprod)
  deallocate(scalprod,stat=i_stat)
  call memocc(i_stat,i_all,'scalprod',subname)

END SUBROUTINE nonlocal_forces


!> Calculates the coefficient of derivative of projectors
subroutine calc_coeff_derproj(l,i,m,nterm_max,rhol,nterm_arr,lxyz_arr,fac_arr)
  implicit none
  integer, intent(in) :: l,i,m,nterm_max
  integer, dimension(3), intent(out) :: nterm_arr
  real(kind=8), intent(in) :: rhol
  integer, dimension(3,nterm_max,3), intent(out) :: lxyz_arr
  real(kind=8), dimension(nterm_max,3), intent(out) :: fac_arr

if (l.eq.1 .and. i.eq.1 .and. m.eq.1) then
   nterm_arr(1)=1
   nterm_arr(2)=1
   nterm_arr(3)=1
   lxyz_arr(1,1,1)=1 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=0
   fac_arr(1,1)=-0.7071067811865475244008444d0/rhol**2d0
   lxyz_arr(1,1,2)=0 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=0
   fac_arr(1,2)=-0.7071067811865475244008444d0/rhol**2d0
   lxyz_arr(1,1,3)=0 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=1
   fac_arr(1,3)=-0.7071067811865475244008444d0/rhol**2d0
else if (l.eq.1 .and. i.eq.2 .and. m.eq.1) then
   nterm_arr(1)=4
   nterm_arr(2)=4
   nterm_arr(3)=4
   lxyz_arr(1,1,1)=1 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=3 ; lxyz_arr(2,2,1)=0 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=2 ; lxyz_arr(3,3,1)=0
   lxyz_arr(1,4,1)=1 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=2
   fac_arr(1,1)=0.730296743340221484609293d0
   fac_arr(2,1)=-0.3651483716701107423046465d0/rhol**2d0
   fac_arr(3,1)=-0.3651483716701107423046465d0/rhol**2d0
   fac_arr(4,1)=-0.3651483716701107423046465d0/rhol**2d0
   lxyz_arr(1,1,2)=0 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=2 ; lxyz_arr(2,2,2)=1 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=3 ; lxyz_arr(3,3,2)=0
   lxyz_arr(1,4,2)=0 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=2
   fac_arr(1,2)=0.730296743340221484609293d0
   fac_arr(2,2)=-0.3651483716701107423046465d0/rhol**2d0
   fac_arr(3,2)=-0.3651483716701107423046465d0/rhol**2d0
   fac_arr(4,2)=-0.3651483716701107423046465d0/rhol**2d0
   lxyz_arr(1,1,3)=0 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=2 ; lxyz_arr(2,2,3)=0 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=0 ; lxyz_arr(2,3,3)=2 ; lxyz_arr(3,3,3)=1
   lxyz_arr(1,4,3)=0 ; lxyz_arr(2,4,3)=0 ; lxyz_arr(3,4,3)=3
   fac_arr(1,3)=0.730296743340221484609293d0
   fac_arr(2,3)=-0.3651483716701107423046465d0/rhol**2d0
   fac_arr(3,3)=-0.3651483716701107423046465d0/rhol**2d0
   fac_arr(4,3)=-0.3651483716701107423046465d0/rhol**2d0
else if (l.eq.1 .and. i.eq.3 .and. m.eq.1) then
   nterm_arr(1)=9
   nterm_arr(2)=9
   nterm_arr(3)=9
   lxyz_arr(1,1,1)=3 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=1 ; lxyz_arr(2,2,1)=2 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=0 ; lxyz_arr(3,3,1)=2
   lxyz_arr(1,4,1)=5 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=0
   lxyz_arr(1,5,1)=3 ; lxyz_arr(2,5,1)=2 ; lxyz_arr(3,5,1)=0
   lxyz_arr(1,6,1)=1 ; lxyz_arr(2,6,1)=4 ; lxyz_arr(3,6,1)=0
   lxyz_arr(1,7,1)=3 ; lxyz_arr(2,7,1)=0 ; lxyz_arr(3,7,1)=2
   lxyz_arr(1,8,1)=1 ; lxyz_arr(2,8,1)=2 ; lxyz_arr(3,8,1)=2
   lxyz_arr(1,9,1)=1 ; lxyz_arr(2,9,1)=0 ; lxyz_arr(3,9,1)=4
   fac_arr(1,1)=0.3680349649825889161579343d0
   fac_arr(2,1)=0.3680349649825889161579343d0
   fac_arr(3,1)=0.3680349649825889161579343d0
   fac_arr(4,1)=-0.09200874124564722903948358d0/rhol**2d0
   fac_arr(5,1)=-0.1840174824912944580789672d0/rhol**2d0
   fac_arr(6,1)=-0.09200874124564722903948358d0/rhol**2d0
   fac_arr(7,1)=-0.1840174824912944580789672d0/rhol**2d0
   fac_arr(8,1)=-0.1840174824912944580789672d0/rhol**2d0
   fac_arr(9,1)=-0.09200874124564722903948358d0/rhol**2d0
   lxyz_arr(1,1,2)=2 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=0 ; lxyz_arr(2,2,2)=3 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=1 ; lxyz_arr(3,3,2)=2
   lxyz_arr(1,4,2)=4 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=0
   lxyz_arr(1,5,2)=2 ; lxyz_arr(2,5,2)=3 ; lxyz_arr(3,5,2)=0
   lxyz_arr(1,6,2)=0 ; lxyz_arr(2,6,2)=5 ; lxyz_arr(3,6,2)=0
   lxyz_arr(1,7,2)=2 ; lxyz_arr(2,7,2)=1 ; lxyz_arr(3,7,2)=2
   lxyz_arr(1,8,2)=0 ; lxyz_arr(2,8,2)=3 ; lxyz_arr(3,8,2)=2
   lxyz_arr(1,9,2)=0 ; lxyz_arr(2,9,2)=1 ; lxyz_arr(3,9,2)=4
   fac_arr(1,2)=0.3680349649825889161579343d0
   fac_arr(2,2)=0.3680349649825889161579343d0
   fac_arr(3,2)=0.3680349649825889161579343d0
   fac_arr(4,2)=-0.09200874124564722903948358d0/rhol**2d0
   fac_arr(5,2)=-0.1840174824912944580789672d0/rhol**2d0
   fac_arr(6,2)=-0.09200874124564722903948358d0/rhol**2d0
   fac_arr(7,2)=-0.1840174824912944580789672d0/rhol**2d0
   fac_arr(8,2)=-0.1840174824912944580789672d0/rhol**2d0
   fac_arr(9,2)=-0.09200874124564722903948358d0/rhol**2d0
   lxyz_arr(1,1,3)=2 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=0 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=0 ; lxyz_arr(2,3,3)=0 ; lxyz_arr(3,3,3)=3
   lxyz_arr(1,4,3)=4 ; lxyz_arr(2,4,3)=0 ; lxyz_arr(3,4,3)=1
   lxyz_arr(1,5,3)=2 ; lxyz_arr(2,5,3)=2 ; lxyz_arr(3,5,3)=1
   lxyz_arr(1,6,3)=0 ; lxyz_arr(2,6,3)=4 ; lxyz_arr(3,6,3)=1
   lxyz_arr(1,7,3)=2 ; lxyz_arr(2,7,3)=0 ; lxyz_arr(3,7,3)=3
   lxyz_arr(1,8,3)=0 ; lxyz_arr(2,8,3)=2 ; lxyz_arr(3,8,3)=3
   lxyz_arr(1,9,3)=0 ; lxyz_arr(2,9,3)=0 ; lxyz_arr(3,9,3)=5
   fac_arr(1,3)=0.3680349649825889161579343d0
   fac_arr(2,3)=0.3680349649825889161579343d0
   fac_arr(3,3)=0.3680349649825889161579343d0
   fac_arr(4,3)=-0.09200874124564722903948358d0/rhol**2d0
   fac_arr(5,3)=-0.1840174824912944580789672d0/rhol**2d0
   fac_arr(6,3)=-0.09200874124564722903948358d0/rhol**2d0
   fac_arr(7,3)=-0.1840174824912944580789672d0/rhol**2d0
   fac_arr(8,3)=-0.1840174824912944580789672d0/rhol**2d0
   fac_arr(9,3)=-0.09200874124564722903948358d0/rhol**2d0
else if (l.eq.2 .and. i.eq.1 .and. m.eq.1) then
   nterm_arr(1)=2
   nterm_arr(2)=1
   nterm_arr(3)=1
   lxyz_arr(1,1,1)=0 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=2 ; lxyz_arr(2,2,1)=0 ; lxyz_arr(3,2,1)=0
   fac_arr(1,1)=1.d0
   fac_arr(2,1)=-1.d0/rhol**2d0
   lxyz_arr(1,1,2)=1 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=0
   fac_arr(1,2)=-1.d0/rhol**2d0
   lxyz_arr(1,1,3)=1 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=1
   fac_arr(1,3)=-1.d0/rhol**2d0
else if (l.eq.2 .and. i.eq.1 .and. m.eq.2) then
   nterm_arr(1)=1
   nterm_arr(2)=2
   nterm_arr(3)=1
   lxyz_arr(1,1,1)=1 ; lxyz_arr(2,1,1)=1 ; lxyz_arr(3,1,1)=0
   fac_arr(1,1)=-1.d0/rhol**2d0
   lxyz_arr(1,1,2)=0 ; lxyz_arr(2,1,2)=0 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=0 ; lxyz_arr(2,2,2)=2 ; lxyz_arr(3,2,2)=0
   fac_arr(1,2)=1.d0
   fac_arr(2,2)=-1.d0/rhol**2d0
   lxyz_arr(1,1,3)=0 ; lxyz_arr(2,1,3)=1 ; lxyz_arr(3,1,3)=1
   fac_arr(1,3)=-1.d0/rhol**2d0
else if (l.eq.2 .and. i.eq.1 .and. m.eq.3) then
   nterm_arr(1)=1
   nterm_arr(2)=1
   nterm_arr(3)=2
   lxyz_arr(1,1,1)=1 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=1
   fac_arr(1,1)=-1.d0/rhol**2d0
   lxyz_arr(1,1,2)=0 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=1
   fac_arr(1,2)=-1.d0/rhol**2d0
   lxyz_arr(1,1,3)=0 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=0
   lxyz_arr(1,2,3)=0 ; lxyz_arr(2,2,3)=0 ; lxyz_arr(3,2,3)=2
   fac_arr(1,3)=1.d0
   fac_arr(2,3)=-1.d0/rhol**2d0
else if (l.eq.2 .and. i.eq.2 .and. m.eq.1) then
   nterm_arr(1)=6
   nterm_arr(2)=4
   nterm_arr(3)=4
   lxyz_arr(1,1,1)=2 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=0 ; lxyz_arr(2,2,1)=2 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=0 ; lxyz_arr(2,3,1)=0 ; lxyz_arr(3,3,1)=2
   lxyz_arr(1,4,1)=4 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=0
   lxyz_arr(1,5,1)=2 ; lxyz_arr(2,5,1)=2 ; lxyz_arr(3,5,1)=0
   lxyz_arr(1,6,1)=2 ; lxyz_arr(2,6,1)=0 ; lxyz_arr(3,6,1)=2
   fac_arr(1,1)=1.014185105674219893011542d0
   fac_arr(2,1)=0.3380617018914066310038473d0
   fac_arr(3,1)=0.3380617018914066310038473d0
   fac_arr(4,1)=-0.3380617018914066310038473d0/rhol**2d0
   fac_arr(5,1)=-0.3380617018914066310038473d0/rhol**2d0
   fac_arr(6,1)=-0.3380617018914066310038473d0/rhol**2d0
   lxyz_arr(1,1,2)=1 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=3 ; lxyz_arr(2,2,2)=1 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=1 ; lxyz_arr(2,3,2)=3 ; lxyz_arr(3,3,2)=0
   lxyz_arr(1,4,2)=1 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=2
   fac_arr(1,2)=0.6761234037828132620076947d0
   fac_arr(2,2)=-0.3380617018914066310038473d0/rhol**2d0
   fac_arr(3,2)=-0.3380617018914066310038473d0/rhol**2d0
   fac_arr(4,2)=-0.3380617018914066310038473d0/rhol**2d0
   lxyz_arr(1,1,3)=1 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=3 ; lxyz_arr(2,2,3)=0 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=1 ; lxyz_arr(2,3,3)=2 ; lxyz_arr(3,3,3)=1
   lxyz_arr(1,4,3)=1 ; lxyz_arr(2,4,3)=0 ; lxyz_arr(3,4,3)=3
   fac_arr(1,3)=0.6761234037828132620076947d0
   fac_arr(2,3)=-0.3380617018914066310038473d0/rhol**2d0
   fac_arr(3,3)=-0.3380617018914066310038473d0/rhol**2d0
   fac_arr(4,3)=-0.3380617018914066310038473d0/rhol**2d0
else if (l.eq.2 .and. i.eq.2 .and. m.eq.2) then
   nterm_arr(1)=4
   nterm_arr(2)=6
   nterm_arr(3)=4
   lxyz_arr(1,1,1)=1 ; lxyz_arr(2,1,1)=1 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=3 ; lxyz_arr(2,2,1)=1 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=3 ; lxyz_arr(3,3,1)=0
   lxyz_arr(1,4,1)=1 ; lxyz_arr(2,4,1)=1 ; lxyz_arr(3,4,1)=2
   fac_arr(1,1)=0.6761234037828132620076947d0
   fac_arr(2,1)=-0.3380617018914066310038473d0/rhol**2d0
   fac_arr(3,1)=-0.3380617018914066310038473d0/rhol**2d0
   fac_arr(4,1)=-0.3380617018914066310038473d0/rhol**2d0
   lxyz_arr(1,1,2)=2 ; lxyz_arr(2,1,2)=0 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=0 ; lxyz_arr(2,2,2)=2 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=0 ; lxyz_arr(3,3,2)=2
   lxyz_arr(1,4,2)=2 ; lxyz_arr(2,4,2)=2 ; lxyz_arr(3,4,2)=0
   lxyz_arr(1,5,2)=0 ; lxyz_arr(2,5,2)=4 ; lxyz_arr(3,5,2)=0
   lxyz_arr(1,6,2)=0 ; lxyz_arr(2,6,2)=2 ; lxyz_arr(3,6,2)=2
   fac_arr(1,2)=0.3380617018914066310038473d0
   fac_arr(2,2)=1.014185105674219893011542d0
   fac_arr(3,2)=0.3380617018914066310038473d0
   fac_arr(4,2)=-0.3380617018914066310038473d0/rhol**2d0
   fac_arr(5,2)=-0.3380617018914066310038473d0/rhol**2d0
   fac_arr(6,2)=-0.3380617018914066310038473d0/rhol**2d0
   lxyz_arr(1,1,3)=0 ; lxyz_arr(2,1,3)=1 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=2 ; lxyz_arr(2,2,3)=1 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=0 ; lxyz_arr(2,3,3)=3 ; lxyz_arr(3,3,3)=1
   lxyz_arr(1,4,3)=0 ; lxyz_arr(2,4,3)=1 ; lxyz_arr(3,4,3)=3
   fac_arr(1,3)=0.6761234037828132620076947d0
   fac_arr(2,3)=-0.3380617018914066310038473d0/rhol**2d0
   fac_arr(3,3)=-0.3380617018914066310038473d0/rhol**2d0
   fac_arr(4,3)=-0.3380617018914066310038473d0/rhol**2d0
else if (l.eq.2 .and. i.eq.2 .and. m.eq.3) then
   nterm_arr(1)=4
   nterm_arr(2)=4
   nterm_arr(3)=6
   lxyz_arr(1,1,1)=1 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=1
   lxyz_arr(1,2,1)=3 ; lxyz_arr(2,2,1)=0 ; lxyz_arr(3,2,1)=1
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=2 ; lxyz_arr(3,3,1)=1
   lxyz_arr(1,4,1)=1 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=3
   fac_arr(1,1)=0.6761234037828132620076947d0
   fac_arr(2,1)=-0.3380617018914066310038473d0/rhol**2d0
   fac_arr(3,1)=-0.3380617018914066310038473d0/rhol**2d0
   fac_arr(4,1)=-0.3380617018914066310038473d0/rhol**2d0
   lxyz_arr(1,1,2)=0 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=1
   lxyz_arr(1,2,2)=2 ; lxyz_arr(2,2,2)=1 ; lxyz_arr(3,2,2)=1
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=3 ; lxyz_arr(3,3,2)=1
   lxyz_arr(1,4,2)=0 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=3
   fac_arr(1,2)=0.6761234037828132620076947d0
   fac_arr(2,2)=-0.3380617018914066310038473d0/rhol**2d0
   fac_arr(3,2)=-0.3380617018914066310038473d0/rhol**2d0
   fac_arr(4,2)=-0.3380617018914066310038473d0/rhol**2d0
   lxyz_arr(1,1,3)=2 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=0
   lxyz_arr(1,2,3)=0 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=0
   lxyz_arr(1,3,3)=0 ; lxyz_arr(2,3,3)=0 ; lxyz_arr(3,3,3)=2
   lxyz_arr(1,4,3)=2 ; lxyz_arr(2,4,3)=0 ; lxyz_arr(3,4,3)=2
   lxyz_arr(1,5,3)=0 ; lxyz_arr(2,5,3)=2 ; lxyz_arr(3,5,3)=2
   lxyz_arr(1,6,3)=0 ; lxyz_arr(2,6,3)=0 ; lxyz_arr(3,6,3)=4
   fac_arr(1,3)=0.3380617018914066310038473d0
   fac_arr(2,3)=0.3380617018914066310038473d0
   fac_arr(3,3)=1.014185105674219893011542d0
   fac_arr(4,3)=-0.3380617018914066310038473d0/rhol**2d0
   fac_arr(5,3)=-0.3380617018914066310038473d0/rhol**2d0
   fac_arr(6,3)=-0.3380617018914066310038473d0/rhol**2d0
else if (l.eq.2 .and. i.eq.3 .and. m.eq.1) then
   nterm_arr(1)=12
   nterm_arr(2)=9
   nterm_arr(3)=9
   lxyz_arr(1,1,1)=4 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=2 ; lxyz_arr(2,2,1)=2 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=0 ; lxyz_arr(2,3,1)=4 ; lxyz_arr(3,3,1)=0
   lxyz_arr(1,4,1)=2 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=2
   lxyz_arr(1,5,1)=0 ; lxyz_arr(2,5,1)=2 ; lxyz_arr(3,5,1)=2
   lxyz_arr(1,6,1)=0 ; lxyz_arr(2,6,1)=0 ; lxyz_arr(3,6,1)=4
   lxyz_arr(1,7,1)=6 ; lxyz_arr(2,7,1)=0 ; lxyz_arr(3,7,1)=0
   lxyz_arr(1,8,1)=4 ; lxyz_arr(2,8,1)=2 ; lxyz_arr(3,8,1)=0
   lxyz_arr(1,9,1)=2 ; lxyz_arr(2,9,1)=4 ; lxyz_arr(3,9,1)=0
   lxyz_arr(1,10,1)=4 ; lxyz_arr(2,10,1)=0 ; lxyz_arr(3,10,1)=2
   lxyz_arr(1,11,1)=2 ; lxyz_arr(2,11,1)=2 ; lxyz_arr(3,11,1)=2
   lxyz_arr(1,12,1)=2 ; lxyz_arr(2,12,1)=0 ; lxyz_arr(3,12,1)=4
   fac_arr(1,1)=0.3397647942917503630913594d0
   fac_arr(2,1)=0.4077177531501004357096312d0
   fac_arr(3,1)=0.06795295885835007261827187d0
   fac_arr(4,1)=0.4077177531501004357096312d0
   fac_arr(5,1)=0.1359059177167001452365437d0
   fac_arr(6,1)=0.06795295885835007261827187d0
   fac_arr(7,1)=-0.06795295885835007261827187d0/rhol**2d0
   fac_arr(8,1)=-0.1359059177167001452365437d0/rhol**2d0
   fac_arr(9,1)=-0.06795295885835007261827187d0/rhol**2d0
   fac_arr(10,1)=-0.1359059177167001452365437d0/rhol**2d0
   fac_arr(11,1)=-0.1359059177167001452365437d0/rhol**2d0
   fac_arr(12,1)=-0.06795295885835007261827187d0/rhol**2d0
   lxyz_arr(1,1,2)=3 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=1 ; lxyz_arr(2,2,2)=3 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=1 ; lxyz_arr(2,3,2)=1 ; lxyz_arr(3,3,2)=2
   lxyz_arr(1,4,2)=5 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=0
   lxyz_arr(1,5,2)=3 ; lxyz_arr(2,5,2)=3 ; lxyz_arr(3,5,2)=0
   lxyz_arr(1,6,2)=1 ; lxyz_arr(2,6,2)=5 ; lxyz_arr(3,6,2)=0
   lxyz_arr(1,7,2)=3 ; lxyz_arr(2,7,2)=1 ; lxyz_arr(3,7,2)=2
   lxyz_arr(1,8,2)=1 ; lxyz_arr(2,8,2)=3 ; lxyz_arr(3,8,2)=2
   lxyz_arr(1,9,2)=1 ; lxyz_arr(2,9,2)=1 ; lxyz_arr(3,9,2)=4
   fac_arr(1,2)=0.2718118354334002904730875d0
   fac_arr(2,2)=0.2718118354334002904730875d0
   fac_arr(3,2)=0.2718118354334002904730875d0
   fac_arr(4,2)=-0.06795295885835007261827187d0/rhol**2d0
   fac_arr(5,2)=-0.1359059177167001452365437d0/rhol**2d0
   fac_arr(6,2)=-0.06795295885835007261827187d0/rhol**2d0
   fac_arr(7,2)=-0.1359059177167001452365437d0/rhol**2d0
   fac_arr(8,2)=-0.1359059177167001452365437d0/rhol**2d0
   fac_arr(9,2)=-0.06795295885835007261827187d0/rhol**2d0
   lxyz_arr(1,1,3)=3 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=1 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=1 ; lxyz_arr(2,3,3)=0 ; lxyz_arr(3,3,3)=3
   lxyz_arr(1,4,3)=5 ; lxyz_arr(2,4,3)=0 ; lxyz_arr(3,4,3)=1
   lxyz_arr(1,5,3)=3 ; lxyz_arr(2,5,3)=2 ; lxyz_arr(3,5,3)=1
   lxyz_arr(1,6,3)=1 ; lxyz_arr(2,6,3)=4 ; lxyz_arr(3,6,3)=1
   lxyz_arr(1,7,3)=3 ; lxyz_arr(2,7,3)=0 ; lxyz_arr(3,7,3)=3
   lxyz_arr(1,8,3)=1 ; lxyz_arr(2,8,3)=2 ; lxyz_arr(3,8,3)=3
   lxyz_arr(1,9,3)=1 ; lxyz_arr(2,9,3)=0 ; lxyz_arr(3,9,3)=5
   fac_arr(1,3)=0.2718118354334002904730875d0
   fac_arr(2,3)=0.2718118354334002904730875d0
   fac_arr(3,3)=0.2718118354334002904730875d0
   fac_arr(4,3)=-0.06795295885835007261827187d0/rhol**2d0
   fac_arr(5,3)=-0.1359059177167001452365437d0/rhol**2d0
   fac_arr(6,3)=-0.06795295885835007261827187d0/rhol**2d0
   fac_arr(7,3)=-0.1359059177167001452365437d0/rhol**2d0
   fac_arr(8,3)=-0.1359059177167001452365437d0/rhol**2d0
   fac_arr(9,3)=-0.06795295885835007261827187d0/rhol**2d0
else if (l.eq.2 .and. i.eq.3 .and. m.eq.2) then
   nterm_arr(1)=9
   nterm_arr(2)=12
   nterm_arr(3)=9
   lxyz_arr(1,1,1)=3 ; lxyz_arr(2,1,1)=1 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=1 ; lxyz_arr(2,2,1)=3 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=1 ; lxyz_arr(3,3,1)=2
   lxyz_arr(1,4,1)=5 ; lxyz_arr(2,4,1)=1 ; lxyz_arr(3,4,1)=0
   lxyz_arr(1,5,1)=3 ; lxyz_arr(2,5,1)=3 ; lxyz_arr(3,5,1)=0
   lxyz_arr(1,6,1)=1 ; lxyz_arr(2,6,1)=5 ; lxyz_arr(3,6,1)=0
   lxyz_arr(1,7,1)=3 ; lxyz_arr(2,7,1)=1 ; lxyz_arr(3,7,1)=2
   lxyz_arr(1,8,1)=1 ; lxyz_arr(2,8,1)=3 ; lxyz_arr(3,8,1)=2
   lxyz_arr(1,9,1)=1 ; lxyz_arr(2,9,1)=1 ; lxyz_arr(3,9,1)=4
   fac_arr(1,1)=0.2718118354334002904730875d0
   fac_arr(2,1)=0.2718118354334002904730875d0
   fac_arr(3,1)=0.2718118354334002904730875d0
   fac_arr(4,1)=-0.06795295885835007261827187d0/rhol**2d0
   fac_arr(5,1)=-0.1359059177167001452365437d0/rhol**2d0
   fac_arr(6,1)=-0.06795295885835007261827187d0/rhol**2d0
   fac_arr(7,1)=-0.1359059177167001452365437d0/rhol**2d0
   fac_arr(8,1)=-0.1359059177167001452365437d0/rhol**2d0
   fac_arr(9,1)=-0.06795295885835007261827187d0/rhol**2d0
   lxyz_arr(1,1,2)=4 ; lxyz_arr(2,1,2)=0 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=2 ; lxyz_arr(2,2,2)=2 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=4 ; lxyz_arr(3,3,2)=0
   lxyz_arr(1,4,2)=2 ; lxyz_arr(2,4,2)=0 ; lxyz_arr(3,4,2)=2
   lxyz_arr(1,5,2)=0 ; lxyz_arr(2,5,2)=2 ; lxyz_arr(3,5,2)=2
   lxyz_arr(1,6,2)=0 ; lxyz_arr(2,6,2)=0 ; lxyz_arr(3,6,2)=4
   lxyz_arr(1,7,2)=4 ; lxyz_arr(2,7,2)=2 ; lxyz_arr(3,7,2)=0
   lxyz_arr(1,8,2)=2 ; lxyz_arr(2,8,2)=4 ; lxyz_arr(3,8,2)=0
   lxyz_arr(1,9,2)=0 ; lxyz_arr(2,9,2)=6 ; lxyz_arr(3,9,2)=0
   lxyz_arr(1,10,2)=2 ; lxyz_arr(2,10,2)=2 ; lxyz_arr(3,10,2)=2
   lxyz_arr(1,11,2)=0 ; lxyz_arr(2,11,2)=4 ; lxyz_arr(3,11,2)=2
   lxyz_arr(1,12,2)=0 ; lxyz_arr(2,12,2)=2 ; lxyz_arr(3,12,2)=4
   fac_arr(1,2)=0.06795295885835007261827187d0
   fac_arr(2,2)=0.4077177531501004357096312d0
   fac_arr(3,2)=0.3397647942917503630913594d0
   fac_arr(4,2)=0.1359059177167001452365437d0
   fac_arr(5,2)=0.4077177531501004357096312d0
   fac_arr(6,2)=0.06795295885835007261827187d0
   fac_arr(7,2)=-0.06795295885835007261827187d0/rhol**2d0
   fac_arr(8,2)=-0.1359059177167001452365437d0/rhol**2d0
   fac_arr(9,2)=-0.06795295885835007261827187d0/rhol**2d0
   fac_arr(10,2)=-0.1359059177167001452365437d0/rhol**2d0
   fac_arr(11,2)=-0.1359059177167001452365437d0/rhol**2d0
   fac_arr(12,2)=-0.06795295885835007261827187d0/rhol**2d0
   lxyz_arr(1,1,3)=2 ; lxyz_arr(2,1,3)=1 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=0 ; lxyz_arr(2,2,3)=3 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=0 ; lxyz_arr(2,3,3)=1 ; lxyz_arr(3,3,3)=3
   lxyz_arr(1,4,3)=4 ; lxyz_arr(2,4,3)=1 ; lxyz_arr(3,4,3)=1
   lxyz_arr(1,5,3)=2 ; lxyz_arr(2,5,3)=3 ; lxyz_arr(3,5,3)=1
   lxyz_arr(1,6,3)=0 ; lxyz_arr(2,6,3)=5 ; lxyz_arr(3,6,3)=1
   lxyz_arr(1,7,3)=2 ; lxyz_arr(2,7,3)=1 ; lxyz_arr(3,7,3)=3
   lxyz_arr(1,8,3)=0 ; lxyz_arr(2,8,3)=3 ; lxyz_arr(3,8,3)=3
   lxyz_arr(1,9,3)=0 ; lxyz_arr(2,9,3)=1 ; lxyz_arr(3,9,3)=5
   fac_arr(1,3)=0.2718118354334002904730875d0
   fac_arr(2,3)=0.2718118354334002904730875d0
   fac_arr(3,3)=0.2718118354334002904730875d0
   fac_arr(4,3)=-0.06795295885835007261827187d0/rhol**2d0
   fac_arr(5,3)=-0.1359059177167001452365437d0/rhol**2d0
   fac_arr(6,3)=-0.06795295885835007261827187d0/rhol**2d0
   fac_arr(7,3)=-0.1359059177167001452365437d0/rhol**2d0
   fac_arr(8,3)=-0.1359059177167001452365437d0/rhol**2d0
   fac_arr(9,3)=-0.06795295885835007261827187d0/rhol**2d0
else if (l.eq.2 .and. i.eq.3 .and. m.eq.3) then
   nterm_arr(1)=9
   nterm_arr(2)=9
   nterm_arr(3)=12
   lxyz_arr(1,1,1)=3 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=1
   lxyz_arr(1,2,1)=1 ; lxyz_arr(2,2,1)=2 ; lxyz_arr(3,2,1)=1
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=0 ; lxyz_arr(3,3,1)=3
   lxyz_arr(1,4,1)=5 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=1
   lxyz_arr(1,5,1)=3 ; lxyz_arr(2,5,1)=2 ; lxyz_arr(3,5,1)=1
   lxyz_arr(1,6,1)=1 ; lxyz_arr(2,6,1)=4 ; lxyz_arr(3,6,1)=1
   lxyz_arr(1,7,1)=3 ; lxyz_arr(2,7,1)=0 ; lxyz_arr(3,7,1)=3
   lxyz_arr(1,8,1)=1 ; lxyz_arr(2,8,1)=2 ; lxyz_arr(3,8,1)=3
   lxyz_arr(1,9,1)=1 ; lxyz_arr(2,9,1)=0 ; lxyz_arr(3,9,1)=5
   fac_arr(1,1)=0.2718118354334002904730875d0
   fac_arr(2,1)=0.2718118354334002904730875d0
   fac_arr(3,1)=0.2718118354334002904730875d0
   fac_arr(4,1)=-0.06795295885835007261827187d0/rhol**2d0
   fac_arr(5,1)=-0.1359059177167001452365437d0/rhol**2d0
   fac_arr(6,1)=-0.06795295885835007261827187d0/rhol**2d0
   fac_arr(7,1)=-0.1359059177167001452365437d0/rhol**2d0
   fac_arr(8,1)=-0.1359059177167001452365437d0/rhol**2d0
   fac_arr(9,1)=-0.06795295885835007261827187d0/rhol**2d0
   lxyz_arr(1,1,2)=2 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=1
   lxyz_arr(1,2,2)=0 ; lxyz_arr(2,2,2)=3 ; lxyz_arr(3,2,2)=1
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=1 ; lxyz_arr(3,3,2)=3
   lxyz_arr(1,4,2)=4 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=1
   lxyz_arr(1,5,2)=2 ; lxyz_arr(2,5,2)=3 ; lxyz_arr(3,5,2)=1
   lxyz_arr(1,6,2)=0 ; lxyz_arr(2,6,2)=5 ; lxyz_arr(3,6,2)=1
   lxyz_arr(1,7,2)=2 ; lxyz_arr(2,7,2)=1 ; lxyz_arr(3,7,2)=3
   lxyz_arr(1,8,2)=0 ; lxyz_arr(2,8,2)=3 ; lxyz_arr(3,8,2)=3
   lxyz_arr(1,9,2)=0 ; lxyz_arr(2,9,2)=1 ; lxyz_arr(3,9,2)=5
   fac_arr(1,2)=0.2718118354334002904730875d0
   fac_arr(2,2)=0.2718118354334002904730875d0
   fac_arr(3,2)=0.2718118354334002904730875d0
   fac_arr(4,2)=-0.06795295885835007261827187d0/rhol**2d0
   fac_arr(5,2)=-0.1359059177167001452365437d0/rhol**2d0
   fac_arr(6,2)=-0.06795295885835007261827187d0/rhol**2d0
   fac_arr(7,2)=-0.1359059177167001452365437d0/rhol**2d0
   fac_arr(8,2)=-0.1359059177167001452365437d0/rhol**2d0
   fac_arr(9,2)=-0.06795295885835007261827187d0/rhol**2d0
   lxyz_arr(1,1,3)=4 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=0
   lxyz_arr(1,2,3)=2 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=0
   lxyz_arr(1,3,3)=0 ; lxyz_arr(2,3,3)=4 ; lxyz_arr(3,3,3)=0
   lxyz_arr(1,4,3)=2 ; lxyz_arr(2,4,3)=0 ; lxyz_arr(3,4,3)=2
   lxyz_arr(1,5,3)=0 ; lxyz_arr(2,5,3)=2 ; lxyz_arr(3,5,3)=2
   lxyz_arr(1,6,3)=0 ; lxyz_arr(2,6,3)=0 ; lxyz_arr(3,6,3)=4
   lxyz_arr(1,7,3)=4 ; lxyz_arr(2,7,3)=0 ; lxyz_arr(3,7,3)=2
   lxyz_arr(1,8,3)=2 ; lxyz_arr(2,8,3)=2 ; lxyz_arr(3,8,3)=2
   lxyz_arr(1,9,3)=0 ; lxyz_arr(2,9,3)=4 ; lxyz_arr(3,9,3)=2
   lxyz_arr(1,10,3)=2 ; lxyz_arr(2,10,3)=0 ; lxyz_arr(3,10,3)=4
   lxyz_arr(1,11,3)=0 ; lxyz_arr(2,11,3)=2 ; lxyz_arr(3,11,3)=4
   lxyz_arr(1,12,3)=0 ; lxyz_arr(2,12,3)=0 ; lxyz_arr(3,12,3)=6
   fac_arr(1,3)=0.06795295885835007261827187d0
   fac_arr(2,3)=0.1359059177167001452365437d0
   fac_arr(3,3)=0.06795295885835007261827187d0
   fac_arr(4,3)=0.4077177531501004357096312d0
   fac_arr(5,3)=0.4077177531501004357096312d0
   fac_arr(6,3)=0.3397647942917503630913594d0
   fac_arr(7,3)=-0.06795295885835007261827187d0/rhol**2d0
   fac_arr(8,3)=-0.1359059177167001452365437d0/rhol**2d0
   fac_arr(9,3)=-0.06795295885835007261827187d0/rhol**2d0
   fac_arr(10,3)=-0.1359059177167001452365437d0/rhol**2d0
   fac_arr(11,3)=-0.1359059177167001452365437d0/rhol**2d0
   fac_arr(12,3)=-0.06795295885835007261827187d0/rhol**2d0
else if (l.eq.3 .and. i.eq.1 .and. m.eq.1) then
   nterm_arr(1)=1
   nterm_arr(2)=2
   nterm_arr(3)=2
   lxyz_arr(1,1,1)=1 ; lxyz_arr(2,1,1)=1 ; lxyz_arr(3,1,1)=1
   fac_arr(1,1)=-1.414213562373095048801689d0/rhol**2d0
   lxyz_arr(1,1,2)=0 ; lxyz_arr(2,1,2)=0 ; lxyz_arr(3,1,2)=1
   lxyz_arr(1,2,2)=0 ; lxyz_arr(2,2,2)=2 ; lxyz_arr(3,2,2)=1
   fac_arr(1,2)=1.414213562373095048801689d0
   fac_arr(2,2)=-1.414213562373095048801689d0/rhol**2d0
   lxyz_arr(1,1,3)=0 ; lxyz_arr(2,1,3)=1 ; lxyz_arr(3,1,3)=0
   lxyz_arr(1,2,3)=0 ; lxyz_arr(2,2,3)=1 ; lxyz_arr(3,2,3)=2
   fac_arr(1,3)=1.414213562373095048801689d0
   fac_arr(2,3)=-1.414213562373095048801689d0/rhol**2d0
else if (l.eq.3 .and. i.eq.1 .and. m.eq.2) then
   nterm_arr(1)=2
   nterm_arr(2)=1
   nterm_arr(3)=2
   lxyz_arr(1,1,1)=0 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=1
   lxyz_arr(1,2,1)=2 ; lxyz_arr(2,2,1)=0 ; lxyz_arr(3,2,1)=1
   fac_arr(1,1)=1.414213562373095048801689d0
   fac_arr(2,1)=-1.414213562373095048801689d0/rhol**2d0
   lxyz_arr(1,1,2)=1 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=1
   fac_arr(1,2)=-1.414213562373095048801689d0/rhol**2d0
   lxyz_arr(1,1,3)=1 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=0
   lxyz_arr(1,2,3)=1 ; lxyz_arr(2,2,3)=0 ; lxyz_arr(3,2,3)=2
   fac_arr(1,3)=1.414213562373095048801689d0
   fac_arr(2,3)=-1.414213562373095048801689d0/rhol**2d0
else if (l.eq.3 .and. i.eq.1 .and. m.eq.3) then
   nterm_arr(1)=2
   nterm_arr(2)=2
   nterm_arr(3)=1
   lxyz_arr(1,1,1)=0 ; lxyz_arr(2,1,1)=1 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=2 ; lxyz_arr(2,2,1)=1 ; lxyz_arr(3,2,1)=0
   fac_arr(1,1)=1.414213562373095048801689d0
   fac_arr(2,1)=-1.414213562373095048801689d0/rhol**2d0
   lxyz_arr(1,1,2)=1 ; lxyz_arr(2,1,2)=0 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=1 ; lxyz_arr(2,2,2)=2 ; lxyz_arr(3,2,2)=0
   fac_arr(1,2)=1.414213562373095048801689d0
   fac_arr(2,2)=-1.414213562373095048801689d0/rhol**2d0
   lxyz_arr(1,1,3)=1 ; lxyz_arr(2,1,3)=1 ; lxyz_arr(3,1,3)=1
   fac_arr(1,3)=-1.414213562373095048801689d0/rhol**2d0
else if (l.eq.3 .and. i.eq.1 .and. m.eq.4) then
   nterm_arr(1)=3
   nterm_arr(2)=3
   nterm_arr(3)=2
   lxyz_arr(1,1,1)=1 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=3 ; lxyz_arr(2,2,1)=0 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=2 ; lxyz_arr(3,3,1)=0
   fac_arr(1,1)=1.414213562373095048801689d0
   fac_arr(2,1)=-0.7071067811865475244008444d0/rhol**2d0
   fac_arr(3,1)=0.7071067811865475244008444d0/rhol**2d0
   lxyz_arr(1,1,2)=0 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=2 ; lxyz_arr(2,2,2)=1 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=3 ; lxyz_arr(3,3,2)=0
   fac_arr(1,2)=-1.414213562373095048801689d0
   fac_arr(2,2)=-0.7071067811865475244008444d0/rhol**2d0
   fac_arr(3,2)=0.7071067811865475244008444d0/rhol**2d0
   lxyz_arr(1,1,3)=2 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=0 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=1
   fac_arr(1,3)=-0.7071067811865475244008444d0/rhol**2d0
   fac_arr(2,3)=0.7071067811865475244008444d0/rhol**2d0
else if (l.eq.3 .and. i.eq.1 .and. m.eq.5) then
   nterm_arr(1)=4
   nterm_arr(2)=4
   nterm_arr(3)=4
   lxyz_arr(1,1,1)=1 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=3 ; lxyz_arr(2,2,1)=0 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=2 ; lxyz_arr(3,3,1)=0
   lxyz_arr(1,4,1)=1 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=2
   fac_arr(1,1)=-0.816496580927726032732428d0
   fac_arr(2,1)=0.408248290463863016366214d0/rhol**2d0
   fac_arr(3,1)=0.408248290463863016366214d0/rhol**2d0
   fac_arr(4,1)=-0.816496580927726032732428d0/rhol**2d0
   lxyz_arr(1,1,2)=0 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=2 ; lxyz_arr(2,2,2)=1 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=3 ; lxyz_arr(3,3,2)=0
   lxyz_arr(1,4,2)=0 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=2
   fac_arr(1,2)=-0.816496580927726032732428d0
   fac_arr(2,2)=0.408248290463863016366214d0/rhol**2d0
   fac_arr(3,2)=0.408248290463863016366214d0/rhol**2d0
   fac_arr(4,2)=-0.816496580927726032732428d0/rhol**2d0
   lxyz_arr(1,1,3)=0 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=2 ; lxyz_arr(2,2,3)=0 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=0 ; lxyz_arr(2,3,3)=2 ; lxyz_arr(3,3,3)=1
   lxyz_arr(1,4,3)=0 ; lxyz_arr(2,4,3)=0 ; lxyz_arr(3,4,3)=3
   fac_arr(1,3)=1.632993161855452065464856d0
   fac_arr(2,3)=0.408248290463863016366214d0/rhol**2d0
   fac_arr(3,3)=0.408248290463863016366214d0/rhol**2d0
   fac_arr(4,3)=-0.816496580927726032732428d0/rhol**2d0
else if (l.eq.3 .and. i.eq.2 .and. m.eq.1) then
   nterm_arr(1)=4
   nterm_arr(2)=6
   nterm_arr(3)=6
   lxyz_arr(1,1,1)=1 ; lxyz_arr(2,1,1)=1 ; lxyz_arr(3,1,1)=1
   lxyz_arr(1,2,1)=3 ; lxyz_arr(2,2,1)=1 ; lxyz_arr(3,2,1)=1
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=3 ; lxyz_arr(3,3,1)=1
   lxyz_arr(1,4,1)=1 ; lxyz_arr(2,4,1)=1 ; lxyz_arr(3,4,1)=3
   fac_arr(1,1)=0.7126966450997983591588093d0
   fac_arr(2,1)=-0.3563483225498991795794046d0/rhol**2d0
   fac_arr(3,1)=-0.3563483225498991795794046d0/rhol**2d0
   fac_arr(4,1)=-0.3563483225498991795794046d0/rhol**2d0
   lxyz_arr(1,1,2)=2 ; lxyz_arr(2,1,2)=0 ; lxyz_arr(3,1,2)=1
   lxyz_arr(1,2,2)=0 ; lxyz_arr(2,2,2)=2 ; lxyz_arr(3,2,2)=1
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=0 ; lxyz_arr(3,3,2)=3
   lxyz_arr(1,4,2)=2 ; lxyz_arr(2,4,2)=2 ; lxyz_arr(3,4,2)=1
   lxyz_arr(1,5,2)=0 ; lxyz_arr(2,5,2)=4 ; lxyz_arr(3,5,2)=1
   lxyz_arr(1,6,2)=0 ; lxyz_arr(2,6,2)=2 ; lxyz_arr(3,6,2)=3
   fac_arr(1,2)=0.3563483225498991795794046d0
   fac_arr(2,2)=1.069044967649697538738214d0
   fac_arr(3,2)=0.3563483225498991795794046d0
   fac_arr(4,2)=-0.3563483225498991795794046d0/rhol**2d0
   fac_arr(5,2)=-0.3563483225498991795794046d0/rhol**2d0
   fac_arr(6,2)=-0.3563483225498991795794046d0/rhol**2d0
   lxyz_arr(1,1,3)=2 ; lxyz_arr(2,1,3)=1 ; lxyz_arr(3,1,3)=0
   lxyz_arr(1,2,3)=0 ; lxyz_arr(2,2,3)=3 ; lxyz_arr(3,2,3)=0
   lxyz_arr(1,3,3)=0 ; lxyz_arr(2,3,3)=1 ; lxyz_arr(3,3,3)=2
   lxyz_arr(1,4,3)=2 ; lxyz_arr(2,4,3)=1 ; lxyz_arr(3,4,3)=2
   lxyz_arr(1,5,3)=0 ; lxyz_arr(2,5,3)=3 ; lxyz_arr(3,5,3)=2
   lxyz_arr(1,6,3)=0 ; lxyz_arr(2,6,3)=1 ; lxyz_arr(3,6,3)=4
   fac_arr(1,3)=0.3563483225498991795794046d0
   fac_arr(2,3)=0.3563483225498991795794046d0
   fac_arr(3,3)=1.069044967649697538738214d0
   fac_arr(4,3)=-0.3563483225498991795794046d0/rhol**2d0
   fac_arr(5,3)=-0.3563483225498991795794046d0/rhol**2d0
   fac_arr(6,3)=-0.3563483225498991795794046d0/rhol**2d0
else if (l.eq.3 .and. i.eq.2 .and. m.eq.2) then
   nterm_arr(1)=6
   nterm_arr(2)=4
   nterm_arr(3)=6
   lxyz_arr(1,1,1)=2 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=1
   lxyz_arr(1,2,1)=0 ; lxyz_arr(2,2,1)=2 ; lxyz_arr(3,2,1)=1
   lxyz_arr(1,3,1)=0 ; lxyz_arr(2,3,1)=0 ; lxyz_arr(3,3,1)=3
   lxyz_arr(1,4,1)=4 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=1
   lxyz_arr(1,5,1)=2 ; lxyz_arr(2,5,1)=2 ; lxyz_arr(3,5,1)=1
   lxyz_arr(1,6,1)=2 ; lxyz_arr(2,6,1)=0 ; lxyz_arr(3,6,1)=3
   fac_arr(1,1)=1.069044967649697538738214d0
   fac_arr(2,1)=0.3563483225498991795794046d0
   fac_arr(3,1)=0.3563483225498991795794046d0
   fac_arr(4,1)=-0.3563483225498991795794046d0/rhol**2d0
   fac_arr(5,1)=-0.3563483225498991795794046d0/rhol**2d0
   fac_arr(6,1)=-0.3563483225498991795794046d0/rhol**2d0
   lxyz_arr(1,1,2)=1 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=1
   lxyz_arr(1,2,2)=3 ; lxyz_arr(2,2,2)=1 ; lxyz_arr(3,2,2)=1
   lxyz_arr(1,3,2)=1 ; lxyz_arr(2,3,2)=3 ; lxyz_arr(3,3,2)=1
   lxyz_arr(1,4,2)=1 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=3
   fac_arr(1,2)=0.7126966450997983591588093d0
   fac_arr(2,2)=-0.3563483225498991795794046d0/rhol**2d0
   fac_arr(3,2)=-0.3563483225498991795794046d0/rhol**2d0
   fac_arr(4,2)=-0.3563483225498991795794046d0/rhol**2d0
   lxyz_arr(1,1,3)=3 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=0
   lxyz_arr(1,2,3)=1 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=0
   lxyz_arr(1,3,3)=1 ; lxyz_arr(2,3,3)=0 ; lxyz_arr(3,3,3)=2
   lxyz_arr(1,4,3)=3 ; lxyz_arr(2,4,3)=0 ; lxyz_arr(3,4,3)=2
   lxyz_arr(1,5,3)=1 ; lxyz_arr(2,5,3)=2 ; lxyz_arr(3,5,3)=2
   lxyz_arr(1,6,3)=1 ; lxyz_arr(2,6,3)=0 ; lxyz_arr(3,6,3)=4
   fac_arr(1,3)=0.3563483225498991795794046d0
   fac_arr(2,3)=0.3563483225498991795794046d0
   fac_arr(3,3)=1.069044967649697538738214d0
   fac_arr(4,3)=-0.3563483225498991795794046d0/rhol**2d0
   fac_arr(5,3)=-0.3563483225498991795794046d0/rhol**2d0
   fac_arr(6,3)=-0.3563483225498991795794046d0/rhol**2d0
else if (l.eq.3 .and. i.eq.2 .and. m.eq.3) then
   nterm_arr(1)=6
   nterm_arr(2)=6
   nterm_arr(3)=4
   lxyz_arr(1,1,1)=2 ; lxyz_arr(2,1,1)=1 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=0 ; lxyz_arr(2,2,1)=3 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=0 ; lxyz_arr(2,3,1)=1 ; lxyz_arr(3,3,1)=2
   lxyz_arr(1,4,1)=4 ; lxyz_arr(2,4,1)=1 ; lxyz_arr(3,4,1)=0
   lxyz_arr(1,5,1)=2 ; lxyz_arr(2,5,1)=3 ; lxyz_arr(3,5,1)=0
   lxyz_arr(1,6,1)=2 ; lxyz_arr(2,6,1)=1 ; lxyz_arr(3,6,1)=2
   fac_arr(1,1)=1.069044967649697538738214d0
   fac_arr(2,1)=0.3563483225498991795794046d0
   fac_arr(3,1)=0.3563483225498991795794046d0
   fac_arr(4,1)=-0.3563483225498991795794046d0/rhol**2d0
   fac_arr(5,1)=-0.3563483225498991795794046d0/rhol**2d0
   fac_arr(6,1)=-0.3563483225498991795794046d0/rhol**2d0
   lxyz_arr(1,1,2)=3 ; lxyz_arr(2,1,2)=0 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=1 ; lxyz_arr(2,2,2)=2 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=1 ; lxyz_arr(2,3,2)=0 ; lxyz_arr(3,3,2)=2
   lxyz_arr(1,4,2)=3 ; lxyz_arr(2,4,2)=2 ; lxyz_arr(3,4,2)=0
   lxyz_arr(1,5,2)=1 ; lxyz_arr(2,5,2)=4 ; lxyz_arr(3,5,2)=0
   lxyz_arr(1,6,2)=1 ; lxyz_arr(2,6,2)=2 ; lxyz_arr(3,6,2)=2
   fac_arr(1,2)=0.3563483225498991795794046d0
   fac_arr(2,2)=1.069044967649697538738214d0
   fac_arr(3,2)=0.3563483225498991795794046d0
   fac_arr(4,2)=-0.3563483225498991795794046d0/rhol**2d0
   fac_arr(5,2)=-0.3563483225498991795794046d0/rhol**2d0
   fac_arr(6,2)=-0.3563483225498991795794046d0/rhol**2d0
   lxyz_arr(1,1,3)=1 ; lxyz_arr(2,1,3)=1 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=3 ; lxyz_arr(2,2,3)=1 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=1 ; lxyz_arr(2,3,3)=3 ; lxyz_arr(3,3,3)=1
   lxyz_arr(1,4,3)=1 ; lxyz_arr(2,4,3)=1 ; lxyz_arr(3,4,3)=3
   fac_arr(1,3)=0.7126966450997983591588093d0
   fac_arr(2,3)=-0.3563483225498991795794046d0/rhol**2d0
   fac_arr(3,3)=-0.3563483225498991795794046d0/rhol**2d0
   fac_arr(4,3)=-0.3563483225498991795794046d0/rhol**2d0
else if (l.eq.3 .and. i.eq.2 .and. m.eq.4) then
   nterm_arr(1)=6
   nterm_arr(2)=6
   nterm_arr(3)=6
   lxyz_arr(1,1,1)=3 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=1 ; lxyz_arr(2,2,1)=0 ; lxyz_arr(3,2,1)=2
   lxyz_arr(1,3,1)=5 ; lxyz_arr(2,3,1)=0 ; lxyz_arr(3,3,1)=0
   lxyz_arr(1,4,1)=1 ; lxyz_arr(2,4,1)=4 ; lxyz_arr(3,4,1)=0
   lxyz_arr(1,5,1)=3 ; lxyz_arr(2,5,1)=0 ; lxyz_arr(3,5,1)=2
   lxyz_arr(1,6,1)=1 ; lxyz_arr(2,6,1)=2 ; lxyz_arr(3,6,1)=2
   fac_arr(1,1)=0.7126966450997983591588093d0
   fac_arr(2,1)=0.3563483225498991795794046d0
   fac_arr(3,1)=-0.1781741612749495897897023d0/rhol**2d0
   fac_arr(4,1)=0.1781741612749495897897023d0/rhol**2d0
   fac_arr(5,1)=-0.1781741612749495897897023d0/rhol**2d0
   fac_arr(6,1)=0.1781741612749495897897023d0/rhol**2d0
   lxyz_arr(1,1,2)=0 ; lxyz_arr(2,1,2)=3 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=0 ; lxyz_arr(2,2,2)=1 ; lxyz_arr(3,2,2)=2
   lxyz_arr(1,3,2)=4 ; lxyz_arr(2,3,2)=1 ; lxyz_arr(3,3,2)=0
   lxyz_arr(1,4,2)=0 ; lxyz_arr(2,4,2)=5 ; lxyz_arr(3,4,2)=0
   lxyz_arr(1,5,2)=2 ; lxyz_arr(2,5,2)=1 ; lxyz_arr(3,5,2)=2
   lxyz_arr(1,6,2)=0 ; lxyz_arr(2,6,2)=3 ; lxyz_arr(3,6,2)=2
   fac_arr(1,2)=-0.7126966450997983591588093d0
   fac_arr(2,2)=-0.3563483225498991795794046d0
   fac_arr(3,2)=-0.1781741612749495897897023d0/rhol**2d0
   fac_arr(4,2)=0.1781741612749495897897023d0/rhol**2d0
   fac_arr(5,2)=-0.1781741612749495897897023d0/rhol**2d0
   fac_arr(6,2)=0.1781741612749495897897023d0/rhol**2d0
   lxyz_arr(1,1,3)=2 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=0 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=4 ; lxyz_arr(2,3,3)=0 ; lxyz_arr(3,3,3)=1
   lxyz_arr(1,4,3)=0 ; lxyz_arr(2,4,3)=4 ; lxyz_arr(3,4,3)=1
   lxyz_arr(1,5,3)=2 ; lxyz_arr(2,5,3)=0 ; lxyz_arr(3,5,3)=3
   lxyz_arr(1,6,3)=0 ; lxyz_arr(2,6,3)=2 ; lxyz_arr(3,6,3)=3
   fac_arr(1,3)=0.3563483225498991795794046d0
   fac_arr(2,3)=-0.3563483225498991795794046d0
   fac_arr(3,3)=-0.1781741612749495897897023d0/rhol**2d0
   fac_arr(4,3)=0.1781741612749495897897023d0/rhol**2d0
   fac_arr(5,3)=-0.1781741612749495897897023d0/rhol**2d0
   fac_arr(6,3)=0.1781741612749495897897023d0/rhol**2d0
else if (l.eq.3 .and. i.eq.2 .and. m.eq.5) then
   nterm_arr(1)=9
   nterm_arr(2)=9
   nterm_arr(3)=9
   lxyz_arr(1,1,1)=3 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=1 ; lxyz_arr(2,2,1)=2 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=0 ; lxyz_arr(3,3,1)=2
   lxyz_arr(1,4,1)=5 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=0
   lxyz_arr(1,5,1)=3 ; lxyz_arr(2,5,1)=2 ; lxyz_arr(3,5,1)=0
   lxyz_arr(1,6,1)=1 ; lxyz_arr(2,6,1)=4 ; lxyz_arr(3,6,1)=0
   lxyz_arr(1,7,1)=3 ; lxyz_arr(2,7,1)=0 ; lxyz_arr(3,7,1)=2
   lxyz_arr(1,8,1)=1 ; lxyz_arr(2,8,1)=2 ; lxyz_arr(3,8,1)=2
   lxyz_arr(1,9,1)=1 ; lxyz_arr(2,9,1)=0 ; lxyz_arr(3,9,1)=4
   fac_arr(1,1)=-0.4114755998989117606962519d0
   fac_arr(2,1)=-0.4114755998989117606962519d0
   fac_arr(3,1)=0.205737799949455880348126d0
   fac_arr(4,1)=0.102868899974727940174063d0/rhol**2d0
   fac_arr(5,1)=0.205737799949455880348126d0/rhol**2d0
   fac_arr(6,1)=0.102868899974727940174063d0/rhol**2d0
   fac_arr(7,1)=-0.102868899974727940174063d0/rhol**2d0
   fac_arr(8,1)=-0.102868899974727940174063d0/rhol**2d0
   fac_arr(9,1)=-0.205737799949455880348126d0/rhol**2d0
   lxyz_arr(1,1,2)=2 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=0 ; lxyz_arr(2,2,2)=3 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=1 ; lxyz_arr(3,3,2)=2
   lxyz_arr(1,4,2)=4 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=0
   lxyz_arr(1,5,2)=2 ; lxyz_arr(2,5,2)=3 ; lxyz_arr(3,5,2)=0
   lxyz_arr(1,6,2)=0 ; lxyz_arr(2,6,2)=5 ; lxyz_arr(3,6,2)=0
   lxyz_arr(1,7,2)=2 ; lxyz_arr(2,7,2)=1 ; lxyz_arr(3,7,2)=2
   lxyz_arr(1,8,2)=0 ; lxyz_arr(2,8,2)=3 ; lxyz_arr(3,8,2)=2
   lxyz_arr(1,9,2)=0 ; lxyz_arr(2,9,2)=1 ; lxyz_arr(3,9,2)=4
   fac_arr(1,2)=-0.4114755998989117606962519d0
   fac_arr(2,2)=-0.4114755998989117606962519d0
   fac_arr(3,2)=0.205737799949455880348126d0
   fac_arr(4,2)=0.102868899974727940174063d0/rhol**2d0
   fac_arr(5,2)=0.205737799949455880348126d0/rhol**2d0
   fac_arr(6,2)=0.102868899974727940174063d0/rhol**2d0
   fac_arr(7,2)=-0.102868899974727940174063d0/rhol**2d0
   fac_arr(8,2)=-0.102868899974727940174063d0/rhol**2d0
   fac_arr(9,2)=-0.205737799949455880348126d0/rhol**2d0
   lxyz_arr(1,1,3)=2 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=0 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=0 ; lxyz_arr(2,3,3)=0 ; lxyz_arr(3,3,3)=3
   lxyz_arr(1,4,3)=4 ; lxyz_arr(2,4,3)=0 ; lxyz_arr(3,4,3)=1
   lxyz_arr(1,5,3)=2 ; lxyz_arr(2,5,3)=2 ; lxyz_arr(3,5,3)=1
   lxyz_arr(1,6,3)=0 ; lxyz_arr(2,6,3)=4 ; lxyz_arr(3,6,3)=1
   lxyz_arr(1,7,3)=2 ; lxyz_arr(2,7,3)=0 ; lxyz_arr(3,7,3)=3
   lxyz_arr(1,8,3)=0 ; lxyz_arr(2,8,3)=2 ; lxyz_arr(3,8,3)=3
   lxyz_arr(1,9,3)=0 ; lxyz_arr(2,9,3)=0 ; lxyz_arr(3,9,3)=5
   fac_arr(1,3)=0.205737799949455880348126d0
   fac_arr(2,3)=0.205737799949455880348126d0
   fac_arr(3,3)=0.8229511997978235213925038d0
   fac_arr(4,3)=0.102868899974727940174063d0/rhol**2d0
   fac_arr(5,3)=0.205737799949455880348126d0/rhol**2d0
   fac_arr(6,3)=0.102868899974727940174063d0/rhol**2d0
   fac_arr(7,3)=-0.102868899974727940174063d0/rhol**2d0
   fac_arr(8,3)=-0.102868899974727940174063d0/rhol**2d0
   fac_arr(9,3)=-0.205737799949455880348126d0/rhol**2d0
else if (l.eq.3 .and. i.eq.3 .and. m.eq.1) then
   nterm_arr(1)=9
   nterm_arr(2)=12
   nterm_arr(3)=12
   lxyz_arr(1,1,1)=3 ; lxyz_arr(2,1,1)=1 ; lxyz_arr(3,1,1)=1
   lxyz_arr(1,2,1)=1 ; lxyz_arr(2,2,1)=3 ; lxyz_arr(3,2,1)=1
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=1 ; lxyz_arr(3,3,1)=3
   lxyz_arr(1,4,1)=5 ; lxyz_arr(2,4,1)=1 ; lxyz_arr(3,4,1)=1
   lxyz_arr(1,5,1)=3 ; lxyz_arr(2,5,1)=3 ; lxyz_arr(3,5,1)=1
   lxyz_arr(1,6,1)=1 ; lxyz_arr(2,6,1)=5 ; lxyz_arr(3,6,1)=1
   lxyz_arr(1,7,1)=3 ; lxyz_arr(2,7,1)=1 ; lxyz_arr(3,7,1)=3
   lxyz_arr(1,8,1)=1 ; lxyz_arr(2,8,1)=3 ; lxyz_arr(3,8,1)=3
   lxyz_arr(1,9,1)=1 ; lxyz_arr(2,9,1)=1 ; lxyz_arr(3,9,1)=5
   fac_arr(1,1)=0.2383947500094262395810797d0
   fac_arr(2,1)=0.2383947500094262395810797d0
   fac_arr(3,1)=0.2383947500094262395810797d0
   fac_arr(4,1)=-0.05959868750235655989526993d0/rhol**2d0
   fac_arr(5,1)=-0.1191973750047131197905399d0/rhol**2d0
   fac_arr(6,1)=-0.05959868750235655989526993d0/rhol**2d0
   fac_arr(7,1)=-0.1191973750047131197905399d0/rhol**2d0
   fac_arr(8,1)=-0.1191973750047131197905399d0/rhol**2d0
   fac_arr(9,1)=-0.05959868750235655989526993d0/rhol**2d0
   lxyz_arr(1,1,2)=4 ; lxyz_arr(2,1,2)=0 ; lxyz_arr(3,1,2)=1
   lxyz_arr(1,2,2)=2 ; lxyz_arr(2,2,2)=2 ; lxyz_arr(3,2,2)=1
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=4 ; lxyz_arr(3,3,2)=1
   lxyz_arr(1,4,2)=2 ; lxyz_arr(2,4,2)=0 ; lxyz_arr(3,4,2)=3
   lxyz_arr(1,5,2)=0 ; lxyz_arr(2,5,2)=2 ; lxyz_arr(3,5,2)=3
   lxyz_arr(1,6,2)=0 ; lxyz_arr(2,6,2)=0 ; lxyz_arr(3,6,2)=5
   lxyz_arr(1,7,2)=4 ; lxyz_arr(2,7,2)=2 ; lxyz_arr(3,7,2)=1
   lxyz_arr(1,8,2)=2 ; lxyz_arr(2,8,2)=4 ; lxyz_arr(3,8,2)=1
   lxyz_arr(1,9,2)=0 ; lxyz_arr(2,9,2)=6 ; lxyz_arr(3,9,2)=1
   lxyz_arr(1,10,2)=2 ; lxyz_arr(2,10,2)=2 ; lxyz_arr(3,10,2)=3
   lxyz_arr(1,11,2)=0 ; lxyz_arr(2,11,2)=4 ; lxyz_arr(3,11,2)=3
   lxyz_arr(1,12,2)=0 ; lxyz_arr(2,12,2)=2 ; lxyz_arr(3,12,2)=5
   fac_arr(1,2)=0.05959868750235655989526993d0
   fac_arr(2,2)=0.3575921250141393593716196d0
   fac_arr(3,2)=0.2979934375117827994763496d0
   fac_arr(4,2)=0.1191973750047131197905399d0
   fac_arr(5,2)=0.3575921250141393593716196d0
   fac_arr(6,2)=0.05959868750235655989526993d0
   fac_arr(7,2)=-0.05959868750235655989526993d0/rhol**2d0
   fac_arr(8,2)=-0.1191973750047131197905399d0/rhol**2d0
   fac_arr(9,2)=-0.05959868750235655989526993d0/rhol**2d0
   fac_arr(10,2)=-0.1191973750047131197905399d0/rhol**2d0
   fac_arr(11,2)=-0.1191973750047131197905399d0/rhol**2d0
   fac_arr(12,2)=-0.05959868750235655989526993d0/rhol**2d0
   lxyz_arr(1,1,3)=4 ; lxyz_arr(2,1,3)=1 ; lxyz_arr(3,1,3)=0
   lxyz_arr(1,2,3)=2 ; lxyz_arr(2,2,3)=3 ; lxyz_arr(3,2,3)=0
   lxyz_arr(1,3,3)=0 ; lxyz_arr(2,3,3)=5 ; lxyz_arr(3,3,3)=0
   lxyz_arr(1,4,3)=2 ; lxyz_arr(2,4,3)=1 ; lxyz_arr(3,4,3)=2
   lxyz_arr(1,5,3)=0 ; lxyz_arr(2,5,3)=3 ; lxyz_arr(3,5,3)=2
   lxyz_arr(1,6,3)=0 ; lxyz_arr(2,6,3)=1 ; lxyz_arr(3,6,3)=4
   lxyz_arr(1,7,3)=4 ; lxyz_arr(2,7,3)=1 ; lxyz_arr(3,7,3)=2
   lxyz_arr(1,8,3)=2 ; lxyz_arr(2,8,3)=3 ; lxyz_arr(3,8,3)=2
   lxyz_arr(1,9,3)=0 ; lxyz_arr(2,9,3)=5 ; lxyz_arr(3,9,3)=2
   lxyz_arr(1,10,3)=2 ; lxyz_arr(2,10,3)=1 ; lxyz_arr(3,10,3)=4
   lxyz_arr(1,11,3)=0 ; lxyz_arr(2,11,3)=3 ; lxyz_arr(3,11,3)=4
   lxyz_arr(1,12,3)=0 ; lxyz_arr(2,12,3)=1 ; lxyz_arr(3,12,3)=6
   fac_arr(1,3)=0.05959868750235655989526993d0
   fac_arr(2,3)=0.1191973750047131197905399d0
   fac_arr(3,3)=0.05959868750235655989526993d0
   fac_arr(4,3)=0.3575921250141393593716196d0
   fac_arr(5,3)=0.3575921250141393593716196d0
   fac_arr(6,3)=0.2979934375117827994763496d0
   fac_arr(7,3)=-0.05959868750235655989526993d0/rhol**2d0
   fac_arr(8,3)=-0.1191973750047131197905399d0/rhol**2d0
   fac_arr(9,3)=-0.05959868750235655989526993d0/rhol**2d0
   fac_arr(10,3)=-0.1191973750047131197905399d0/rhol**2d0
   fac_arr(11,3)=-0.1191973750047131197905399d0/rhol**2d0
   fac_arr(12,3)=-0.05959868750235655989526993d0/rhol**2d0
else if (l.eq.3 .and. i.eq.3 .and. m.eq.2) then
   nterm_arr(1)=12
   nterm_arr(2)=9
   nterm_arr(3)=12
   lxyz_arr(1,1,1)=4 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=1
   lxyz_arr(1,2,1)=2 ; lxyz_arr(2,2,1)=2 ; lxyz_arr(3,2,1)=1
   lxyz_arr(1,3,1)=0 ; lxyz_arr(2,3,1)=4 ; lxyz_arr(3,3,1)=1
   lxyz_arr(1,4,1)=2 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=3
   lxyz_arr(1,5,1)=0 ; lxyz_arr(2,5,1)=2 ; lxyz_arr(3,5,1)=3
   lxyz_arr(1,6,1)=0 ; lxyz_arr(2,6,1)=0 ; lxyz_arr(3,6,1)=5
   lxyz_arr(1,7,1)=6 ; lxyz_arr(2,7,1)=0 ; lxyz_arr(3,7,1)=1
   lxyz_arr(1,8,1)=4 ; lxyz_arr(2,8,1)=2 ; lxyz_arr(3,8,1)=1
   lxyz_arr(1,9,1)=2 ; lxyz_arr(2,9,1)=4 ; lxyz_arr(3,9,1)=1
   lxyz_arr(1,10,1)=4 ; lxyz_arr(2,10,1)=0 ; lxyz_arr(3,10,1)=3
   lxyz_arr(1,11,1)=2 ; lxyz_arr(2,11,1)=2 ; lxyz_arr(3,11,1)=3
   lxyz_arr(1,12,1)=2 ; lxyz_arr(2,12,1)=0 ; lxyz_arr(3,12,1)=5
   fac_arr(1,1)=0.2979934375117827994763496d0
   fac_arr(2,1)=0.3575921250141393593716196d0
   fac_arr(3,1)=0.05959868750235655989526993d0
   fac_arr(4,1)=0.3575921250141393593716196d0
   fac_arr(5,1)=0.1191973750047131197905399d0
   fac_arr(6,1)=0.05959868750235655989526993d0
   fac_arr(7,1)=-0.05959868750235655989526993d0/rhol**2d0
   fac_arr(8,1)=-0.1191973750047131197905399d0/rhol**2d0
   fac_arr(9,1)=-0.05959868750235655989526993d0/rhol**2d0
   fac_arr(10,1)=-0.1191973750047131197905399d0/rhol**2d0
   fac_arr(11,1)=-0.1191973750047131197905399d0/rhol**2d0
   fac_arr(12,1)=-0.05959868750235655989526993d0/rhol**2d0
   lxyz_arr(1,1,2)=3 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=1
   lxyz_arr(1,2,2)=1 ; lxyz_arr(2,2,2)=3 ; lxyz_arr(3,2,2)=1
   lxyz_arr(1,3,2)=1 ; lxyz_arr(2,3,2)=1 ; lxyz_arr(3,3,2)=3
   lxyz_arr(1,4,2)=5 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=1
   lxyz_arr(1,5,2)=3 ; lxyz_arr(2,5,2)=3 ; lxyz_arr(3,5,2)=1
   lxyz_arr(1,6,2)=1 ; lxyz_arr(2,6,2)=5 ; lxyz_arr(3,6,2)=1
   lxyz_arr(1,7,2)=3 ; lxyz_arr(2,7,2)=1 ; lxyz_arr(3,7,2)=3
   lxyz_arr(1,8,2)=1 ; lxyz_arr(2,8,2)=3 ; lxyz_arr(3,8,2)=3
   lxyz_arr(1,9,2)=1 ; lxyz_arr(2,9,2)=1 ; lxyz_arr(3,9,2)=5
   fac_arr(1,2)=0.2383947500094262395810797d0
   fac_arr(2,2)=0.2383947500094262395810797d0
   fac_arr(3,2)=0.2383947500094262395810797d0
   fac_arr(4,2)=-0.05959868750235655989526993d0/rhol**2d0
   fac_arr(5,2)=-0.1191973750047131197905399d0/rhol**2d0
   fac_arr(6,2)=-0.05959868750235655989526993d0/rhol**2d0
   fac_arr(7,2)=-0.1191973750047131197905399d0/rhol**2d0
   fac_arr(8,2)=-0.1191973750047131197905399d0/rhol**2d0
   fac_arr(9,2)=-0.05959868750235655989526993d0/rhol**2d0
   lxyz_arr(1,1,3)=5 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=0
   lxyz_arr(1,2,3)=3 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=0
   lxyz_arr(1,3,3)=1 ; lxyz_arr(2,3,3)=4 ; lxyz_arr(3,3,3)=0
   lxyz_arr(1,4,3)=3 ; lxyz_arr(2,4,3)=0 ; lxyz_arr(3,4,3)=2
   lxyz_arr(1,5,3)=1 ; lxyz_arr(2,5,3)=2 ; lxyz_arr(3,5,3)=2
   lxyz_arr(1,6,3)=1 ; lxyz_arr(2,6,3)=0 ; lxyz_arr(3,6,3)=4
   lxyz_arr(1,7,3)=5 ; lxyz_arr(2,7,3)=0 ; lxyz_arr(3,7,3)=2
   lxyz_arr(1,8,3)=3 ; lxyz_arr(2,8,3)=2 ; lxyz_arr(3,8,3)=2
   lxyz_arr(1,9,3)=1 ; lxyz_arr(2,9,3)=4 ; lxyz_arr(3,9,3)=2
   lxyz_arr(1,10,3)=3 ; lxyz_arr(2,10,3)=0 ; lxyz_arr(3,10,3)=4
   lxyz_arr(1,11,3)=1 ; lxyz_arr(2,11,3)=2 ; lxyz_arr(3,11,3)=4
   lxyz_arr(1,12,3)=1 ; lxyz_arr(2,12,3)=0 ; lxyz_arr(3,12,3)=6
   fac_arr(1,3)=0.05959868750235655989526993d0
   fac_arr(2,3)=0.1191973750047131197905399d0
   fac_arr(3,3)=0.05959868750235655989526993d0
   fac_arr(4,3)=0.3575921250141393593716196d0
   fac_arr(5,3)=0.3575921250141393593716196d0
   fac_arr(6,3)=0.2979934375117827994763496d0
   fac_arr(7,3)=-0.05959868750235655989526993d0/rhol**2d0
   fac_arr(8,3)=-0.1191973750047131197905399d0/rhol**2d0
   fac_arr(9,3)=-0.05959868750235655989526993d0/rhol**2d0
   fac_arr(10,3)=-0.1191973750047131197905399d0/rhol**2d0
   fac_arr(11,3)=-0.1191973750047131197905399d0/rhol**2d0
   fac_arr(12,3)=-0.05959868750235655989526993d0/rhol**2d0
else if (l.eq.3 .and. i.eq.3 .and. m.eq.3) then
   nterm_arr(1)=12
   nterm_arr(2)=12
   nterm_arr(3)=9
   lxyz_arr(1,1,1)=4 ; lxyz_arr(2,1,1)=1 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=2 ; lxyz_arr(2,2,1)=3 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=0 ; lxyz_arr(2,3,1)=5 ; lxyz_arr(3,3,1)=0
   lxyz_arr(1,4,1)=2 ; lxyz_arr(2,4,1)=1 ; lxyz_arr(3,4,1)=2
   lxyz_arr(1,5,1)=0 ; lxyz_arr(2,5,1)=3 ; lxyz_arr(3,5,1)=2
   lxyz_arr(1,6,1)=0 ; lxyz_arr(2,6,1)=1 ; lxyz_arr(3,6,1)=4
   lxyz_arr(1,7,1)=6 ; lxyz_arr(2,7,1)=1 ; lxyz_arr(3,7,1)=0
   lxyz_arr(1,8,1)=4 ; lxyz_arr(2,8,1)=3 ; lxyz_arr(3,8,1)=0
   lxyz_arr(1,9,1)=2 ; lxyz_arr(2,9,1)=5 ; lxyz_arr(3,9,1)=0
   lxyz_arr(1,10,1)=4 ; lxyz_arr(2,10,1)=1 ; lxyz_arr(3,10,1)=2
   lxyz_arr(1,11,1)=2 ; lxyz_arr(2,11,1)=3 ; lxyz_arr(3,11,1)=2
   lxyz_arr(1,12,1)=2 ; lxyz_arr(2,12,1)=1 ; lxyz_arr(3,12,1)=4
   fac_arr(1,1)=0.2979934375117827994763496d0
   fac_arr(2,1)=0.3575921250141393593716196d0
   fac_arr(3,1)=0.05959868750235655989526993d0
   fac_arr(4,1)=0.3575921250141393593716196d0
   fac_arr(5,1)=0.1191973750047131197905399d0
   fac_arr(6,1)=0.05959868750235655989526993d0
   fac_arr(7,1)=-0.05959868750235655989526993d0/rhol**2d0
   fac_arr(8,1)=-0.1191973750047131197905399d0/rhol**2d0
   fac_arr(9,1)=-0.05959868750235655989526993d0/rhol**2d0
   fac_arr(10,1)=-0.1191973750047131197905399d0/rhol**2d0
   fac_arr(11,1)=-0.1191973750047131197905399d0/rhol**2d0
   fac_arr(12,1)=-0.05959868750235655989526993d0/rhol**2d0
   lxyz_arr(1,1,2)=5 ; lxyz_arr(2,1,2)=0 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=3 ; lxyz_arr(2,2,2)=2 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=1 ; lxyz_arr(2,3,2)=4 ; lxyz_arr(3,3,2)=0
   lxyz_arr(1,4,2)=3 ; lxyz_arr(2,4,2)=0 ; lxyz_arr(3,4,2)=2
   lxyz_arr(1,5,2)=1 ; lxyz_arr(2,5,2)=2 ; lxyz_arr(3,5,2)=2
   lxyz_arr(1,6,2)=1 ; lxyz_arr(2,6,2)=0 ; lxyz_arr(3,6,2)=4
   lxyz_arr(1,7,2)=5 ; lxyz_arr(2,7,2)=2 ; lxyz_arr(3,7,2)=0
   lxyz_arr(1,8,2)=3 ; lxyz_arr(2,8,2)=4 ; lxyz_arr(3,8,2)=0
   lxyz_arr(1,9,2)=1 ; lxyz_arr(2,9,2)=6 ; lxyz_arr(3,9,2)=0
   lxyz_arr(1,10,2)=3 ; lxyz_arr(2,10,2)=2 ; lxyz_arr(3,10,2)=2
   lxyz_arr(1,11,2)=1 ; lxyz_arr(2,11,2)=4 ; lxyz_arr(3,11,2)=2
   lxyz_arr(1,12,2)=1 ; lxyz_arr(2,12,2)=2 ; lxyz_arr(3,12,2)=4
   fac_arr(1,2)=0.05959868750235655989526993d0
   fac_arr(2,2)=0.3575921250141393593716196d0
   fac_arr(3,2)=0.2979934375117827994763496d0
   fac_arr(4,2)=0.1191973750047131197905399d0
   fac_arr(5,2)=0.3575921250141393593716196d0
   fac_arr(6,2)=0.05959868750235655989526993d0
   fac_arr(7,2)=-0.05959868750235655989526993d0/rhol**2d0
   fac_arr(8,2)=-0.1191973750047131197905399d0/rhol**2d0
   fac_arr(9,2)=-0.05959868750235655989526993d0/rhol**2d0
   fac_arr(10,2)=-0.1191973750047131197905399d0/rhol**2d0
   fac_arr(11,2)=-0.1191973750047131197905399d0/rhol**2d0
   fac_arr(12,2)=-0.05959868750235655989526993d0/rhol**2d0
   lxyz_arr(1,1,3)=3 ; lxyz_arr(2,1,3)=1 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=1 ; lxyz_arr(2,2,3)=3 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=1 ; lxyz_arr(2,3,3)=1 ; lxyz_arr(3,3,3)=3
   lxyz_arr(1,4,3)=5 ; lxyz_arr(2,4,3)=1 ; lxyz_arr(3,4,3)=1
   lxyz_arr(1,5,3)=3 ; lxyz_arr(2,5,3)=3 ; lxyz_arr(3,5,3)=1
   lxyz_arr(1,6,3)=1 ; lxyz_arr(2,6,3)=5 ; lxyz_arr(3,6,3)=1
   lxyz_arr(1,7,3)=3 ; lxyz_arr(2,7,3)=1 ; lxyz_arr(3,7,3)=3
   lxyz_arr(1,8,3)=1 ; lxyz_arr(2,8,3)=3 ; lxyz_arr(3,8,3)=3
   lxyz_arr(1,9,3)=1 ; lxyz_arr(2,9,3)=1 ; lxyz_arr(3,9,3)=5
   fac_arr(1,3)=0.2383947500094262395810797d0
   fac_arr(2,3)=0.2383947500094262395810797d0
   fac_arr(3,3)=0.2383947500094262395810797d0
   fac_arr(4,3)=-0.05959868750235655989526993d0/rhol**2d0
   fac_arr(5,3)=-0.1191973750047131197905399d0/rhol**2d0
   fac_arr(6,3)=-0.05959868750235655989526993d0/rhol**2d0
   fac_arr(7,3)=-0.1191973750047131197905399d0/rhol**2d0
   fac_arr(8,3)=-0.1191973750047131197905399d0/rhol**2d0
   fac_arr(9,3)=-0.05959868750235655989526993d0/rhol**2d0
else if (l.eq.3 .and. i.eq.3 .and. m.eq.4) then
   nterm_arr(1)=13
   nterm_arr(2)=13
   nterm_arr(3)=12
   lxyz_arr(1,1,1)=5 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=3 ; lxyz_arr(2,2,1)=2 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=4 ; lxyz_arr(3,3,1)=0
   lxyz_arr(1,4,1)=3 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=2
   lxyz_arr(1,5,1)=1 ; lxyz_arr(2,5,1)=0 ; lxyz_arr(3,5,1)=4
   lxyz_arr(1,6,1)=7 ; lxyz_arr(2,6,1)=0 ; lxyz_arr(3,6,1)=0
   lxyz_arr(1,7,1)=5 ; lxyz_arr(2,7,1)=2 ; lxyz_arr(3,7,1)=0
   lxyz_arr(1,8,1)=3 ; lxyz_arr(2,8,1)=4 ; lxyz_arr(3,8,1)=0
   lxyz_arr(1,9,1)=1 ; lxyz_arr(2,9,1)=6 ; lxyz_arr(3,9,1)=0
   lxyz_arr(1,10,1)=5 ; lxyz_arr(2,10,1)=0 ; lxyz_arr(3,10,1)=2
   lxyz_arr(1,11,1)=1 ; lxyz_arr(2,11,1)=4 ; lxyz_arr(3,11,1)=2
   lxyz_arr(1,12,1)=3 ; lxyz_arr(2,12,1)=0 ; lxyz_arr(3,12,1)=4
   lxyz_arr(1,13,1)=1 ; lxyz_arr(2,13,1)=2 ; lxyz_arr(3,13,1)=4
   fac_arr(1,1)=0.1787960625070696796858098d0
   fac_arr(2,1)=0.1191973750047131197905399d0
   fac_arr(3,1)=-0.05959868750235655989526993d0
   fac_arr(4,1)=0.2383947500094262395810797d0
   fac_arr(5,1)=0.05959868750235655989526993d0
   fac_arr(6,1)=-0.02979934375117827994763496d0/rhol**2d0
   fac_arr(7,1)=-0.02979934375117827994763496d0/rhol**2d0
   fac_arr(8,1)=0.02979934375117827994763496d0/rhol**2d0
   fac_arr(9,1)=0.02979934375117827994763496d0/rhol**2d0
   fac_arr(10,1)=-0.05959868750235655989526993d0/rhol**2d0
   fac_arr(11,1)=0.05959868750235655989526993d0/rhol**2d0
   fac_arr(12,1)=-0.02979934375117827994763496d0/rhol**2d0
   fac_arr(13,1)=0.02979934375117827994763496d0/rhol**2d0
   lxyz_arr(1,1,2)=4 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=2 ; lxyz_arr(2,2,2)=3 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=5 ; lxyz_arr(3,3,2)=0
   lxyz_arr(1,4,2)=0 ; lxyz_arr(2,4,2)=3 ; lxyz_arr(3,4,2)=2
   lxyz_arr(1,5,2)=0 ; lxyz_arr(2,5,2)=1 ; lxyz_arr(3,5,2)=4
   lxyz_arr(1,6,2)=6 ; lxyz_arr(2,6,2)=1 ; lxyz_arr(3,6,2)=0
   lxyz_arr(1,7,2)=4 ; lxyz_arr(2,7,2)=3 ; lxyz_arr(3,7,2)=0
   lxyz_arr(1,8,2)=2 ; lxyz_arr(2,8,2)=5 ; lxyz_arr(3,8,2)=0
   lxyz_arr(1,9,2)=0 ; lxyz_arr(2,9,2)=7 ; lxyz_arr(3,9,2)=0
   lxyz_arr(1,10,2)=4 ; lxyz_arr(2,10,2)=1 ; lxyz_arr(3,10,2)=2
   lxyz_arr(1,11,2)=0 ; lxyz_arr(2,11,2)=5 ; lxyz_arr(3,11,2)=2
   lxyz_arr(1,12,2)=2 ; lxyz_arr(2,12,2)=1 ; lxyz_arr(3,12,2)=4
   lxyz_arr(1,13,2)=0 ; lxyz_arr(2,13,2)=3 ; lxyz_arr(3,13,2)=4
   fac_arr(1,2)=0.05959868750235655989526993d0
   fac_arr(2,2)=-0.1191973750047131197905399d0
   fac_arr(3,2)=-0.1787960625070696796858098d0
   fac_arr(4,2)=-0.2383947500094262395810797d0
   fac_arr(5,2)=-0.05959868750235655989526993d0
   fac_arr(6,2)=-0.02979934375117827994763496d0/rhol**2d0
   fac_arr(7,2)=-0.02979934375117827994763496d0/rhol**2d0
   fac_arr(8,2)=0.02979934375117827994763496d0/rhol**2d0
   fac_arr(9,2)=0.02979934375117827994763496d0/rhol**2d0
   fac_arr(10,2)=-0.05959868750235655989526993d0/rhol**2d0
   fac_arr(11,2)=0.05959868750235655989526993d0/rhol**2d0
   fac_arr(12,2)=-0.02979934375117827994763496d0/rhol**2d0
   fac_arr(13,2)=0.02979934375117827994763496d0/rhol**2d0
   lxyz_arr(1,1,3)=4 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=0 ; lxyz_arr(2,2,3)=4 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=2 ; lxyz_arr(2,3,3)=0 ; lxyz_arr(3,3,3)=3
   lxyz_arr(1,4,3)=0 ; lxyz_arr(2,4,3)=2 ; lxyz_arr(3,4,3)=3
   lxyz_arr(1,5,3)=6 ; lxyz_arr(2,5,3)=0 ; lxyz_arr(3,5,3)=1
   lxyz_arr(1,6,3)=4 ; lxyz_arr(2,6,3)=2 ; lxyz_arr(3,6,3)=1
   lxyz_arr(1,7,3)=2 ; lxyz_arr(2,7,3)=4 ; lxyz_arr(3,7,3)=1
   lxyz_arr(1,8,3)=0 ; lxyz_arr(2,8,3)=6 ; lxyz_arr(3,8,3)=1
   lxyz_arr(1,9,3)=4 ; lxyz_arr(2,9,3)=0 ; lxyz_arr(3,9,3)=3
   lxyz_arr(1,10,3)=0 ; lxyz_arr(2,10,3)=4 ; lxyz_arr(3,10,3)=3
   lxyz_arr(1,11,3)=2 ; lxyz_arr(2,11,3)=0 ; lxyz_arr(3,11,3)=5
   lxyz_arr(1,12,3)=0 ; lxyz_arr(2,12,3)=2 ; lxyz_arr(3,12,3)=5
   fac_arr(1,3)=0.1191973750047131197905399d0
   fac_arr(2,3)=-0.1191973750047131197905399d0
   fac_arr(3,3)=0.1191973750047131197905399d0
   fac_arr(4,3)=-0.1191973750047131197905399d0
   fac_arr(5,3)=-0.02979934375117827994763496d0/rhol**2d0
   fac_arr(6,3)=-0.02979934375117827994763496d0/rhol**2d0
   fac_arr(7,3)=0.02979934375117827994763496d0/rhol**2d0
   fac_arr(8,3)=0.02979934375117827994763496d0/rhol**2d0
   fac_arr(9,3)=-0.05959868750235655989526993d0/rhol**2d0
   fac_arr(10,3)=0.05959868750235655989526993d0/rhol**2d0
   fac_arr(11,3)=-0.02979934375117827994763496d0/rhol**2d0
   fac_arr(12,3)=0.02979934375117827994763496d0/rhol**2d0
else if (l.eq.3 .and. i.eq.3 .and. m.eq.5) then
   nterm_arr(1)=11
   nterm_arr(2)=11
   nterm_arr(3)=10
   lxyz_arr(1,1,1)=5 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=3 ; lxyz_arr(2,2,1)=2 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=4 ; lxyz_arr(3,3,1)=0
   lxyz_arr(1,4,1)=1 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=4
   lxyz_arr(1,5,1)=7 ; lxyz_arr(2,5,1)=0 ; lxyz_arr(3,5,1)=0
   lxyz_arr(1,6,1)=5 ; lxyz_arr(2,6,1)=2 ; lxyz_arr(3,6,1)=0
   lxyz_arr(1,7,1)=3 ; lxyz_arr(2,7,1)=4 ; lxyz_arr(3,7,1)=0
   lxyz_arr(1,8,1)=1 ; lxyz_arr(2,8,1)=6 ; lxyz_arr(3,8,1)=0
   lxyz_arr(1,9,1)=3 ; lxyz_arr(2,9,1)=0 ; lxyz_arr(3,9,1)=4
   lxyz_arr(1,10,1)=1 ; lxyz_arr(2,10,1)=2 ; lxyz_arr(3,10,1)=4
   lxyz_arr(1,11,1)=1 ; lxyz_arr(2,11,1)=0 ; lxyz_arr(3,11,1)=6
   fac_arr(1,1)=-0.1032279548185018340124748d0
   fac_arr(2,1)=-0.2064559096370036680249495d0
   fac_arr(3,1)=-0.1032279548185018340124748d0
   fac_arr(4,1)=0.1032279548185018340124748d0
   fac_arr(5,1)=0.01720465913641697233541246d0/rhol**2d0
   fac_arr(6,1)=0.05161397740925091700623738d0/rhol**2d0
   fac_arr(7,1)=0.05161397740925091700623738d0/rhol**2d0
   fac_arr(8,1)=0.01720465913641697233541246d0/rhol**2d0
   fac_arr(9,1)=-0.05161397740925091700623738d0/rhol**2d0
   fac_arr(10,1)=-0.05161397740925091700623738d0/rhol**2d0
   fac_arr(11,1)=-0.03440931827283394467082492d0/rhol**2d0
   lxyz_arr(1,1,2)=4 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=2 ; lxyz_arr(2,2,2)=3 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=5 ; lxyz_arr(3,3,2)=0
   lxyz_arr(1,4,2)=0 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=4
   lxyz_arr(1,5,2)=6 ; lxyz_arr(2,5,2)=1 ; lxyz_arr(3,5,2)=0
   lxyz_arr(1,6,2)=4 ; lxyz_arr(2,6,2)=3 ; lxyz_arr(3,6,2)=0
   lxyz_arr(1,7,2)=2 ; lxyz_arr(2,7,2)=5 ; lxyz_arr(3,7,2)=0
   lxyz_arr(1,8,2)=0 ; lxyz_arr(2,8,2)=7 ; lxyz_arr(3,8,2)=0
   lxyz_arr(1,9,2)=2 ; lxyz_arr(2,9,2)=1 ; lxyz_arr(3,9,2)=4
   lxyz_arr(1,10,2)=0 ; lxyz_arr(2,10,2)=3 ; lxyz_arr(3,10,2)=4
   lxyz_arr(1,11,2)=0 ; lxyz_arr(2,11,2)=1 ; lxyz_arr(3,11,2)=6
   fac_arr(1,2)=-0.1032279548185018340124748d0
   fac_arr(2,2)=-0.2064559096370036680249495d0
   fac_arr(3,2)=-0.1032279548185018340124748d0
   fac_arr(4,2)=0.1032279548185018340124748d0
   fac_arr(5,2)=0.01720465913641697233541246d0/rhol**2d0
   fac_arr(6,2)=0.05161397740925091700623738d0/rhol**2d0
   fac_arr(7,2)=0.05161397740925091700623738d0/rhol**2d0
   fac_arr(8,2)=0.01720465913641697233541246d0/rhol**2d0
   fac_arr(9,2)=-0.05161397740925091700623738d0/rhol**2d0
   fac_arr(10,2)=-0.05161397740925091700623738d0/rhol**2d0
   fac_arr(11,2)=-0.03440931827283394467082492d0/rhol**2d0
   lxyz_arr(1,1,3)=2 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=3
   lxyz_arr(1,2,3)=0 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=3
   lxyz_arr(1,3,3)=0 ; lxyz_arr(2,3,3)=0 ; lxyz_arr(3,3,3)=5
   lxyz_arr(1,4,3)=6 ; lxyz_arr(2,4,3)=0 ; lxyz_arr(3,4,3)=1
   lxyz_arr(1,5,3)=4 ; lxyz_arr(2,5,3)=2 ; lxyz_arr(3,5,3)=1
   lxyz_arr(1,6,3)=2 ; lxyz_arr(2,6,3)=4 ; lxyz_arr(3,6,3)=1
   lxyz_arr(1,7,3)=0 ; lxyz_arr(2,7,3)=6 ; lxyz_arr(3,7,3)=1
   lxyz_arr(1,8,3)=2 ; lxyz_arr(2,8,3)=0 ; lxyz_arr(3,8,3)=5
   lxyz_arr(1,9,3)=0 ; lxyz_arr(2,9,3)=2 ; lxyz_arr(3,9,3)=5
   lxyz_arr(1,10,3)=0 ; lxyz_arr(2,10,3)=0 ; lxyz_arr(3,10,3)=7
   fac_arr(1,3)=0.2064559096370036680249495d0
   fac_arr(2,3)=0.2064559096370036680249495d0
   fac_arr(3,3)=0.2064559096370036680249495d0
   fac_arr(4,3)=0.01720465913641697233541246d0/rhol**2d0
   fac_arr(5,3)=0.05161397740925091700623738d0/rhol**2d0
   fac_arr(6,3)=0.05161397740925091700623738d0/rhol**2d0
   fac_arr(7,3)=0.01720465913641697233541246d0/rhol**2d0
   fac_arr(8,3)=-0.05161397740925091700623738d0/rhol**2d0
   fac_arr(9,3)=-0.05161397740925091700623738d0/rhol**2d0
   fac_arr(10,3)=-0.03440931827283394467082492d0/rhol**2d0
else if (l.eq.4 .and. i.eq.1 .and. m.eq.1) then
   nterm_arr(1)=6
   nterm_arr(2)=4
   nterm_arr(3)=4
   lxyz_arr(1,1,1)=2 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=0 ; lxyz_arr(2,2,1)=2 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=0 ; lxyz_arr(2,3,1)=0 ; lxyz_arr(3,3,1)=2
   lxyz_arr(1,4,1)=4 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=0
   lxyz_arr(1,5,1)=2 ; lxyz_arr(2,5,1)=2 ; lxyz_arr(3,5,1)=0
   lxyz_arr(1,6,1)=2 ; lxyz_arr(2,6,1)=0 ; lxyz_arr(3,6,1)=2
   fac_arr(1,1)=0.9486832980505137995996681d0
   fac_arr(2,1)=0.3162277660168379331998894d0
   fac_arr(3,1)=-1.264911064067351732799557d0
   fac_arr(4,1)=-0.3162277660168379331998894d0/rhol**2d0
   fac_arr(5,1)=-0.3162277660168379331998894d0/rhol**2d0
   fac_arr(6,1)=1.264911064067351732799557d0/rhol**2d0
   lxyz_arr(1,1,2)=1 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=3 ; lxyz_arr(2,2,2)=1 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=1 ; lxyz_arr(2,3,2)=3 ; lxyz_arr(3,3,2)=0
   lxyz_arr(1,4,2)=1 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=2
   fac_arr(1,2)=0.6324555320336758663997787d0
   fac_arr(2,2)=-0.3162277660168379331998894d0/rhol**2d0
   fac_arr(3,2)=-0.3162277660168379331998894d0/rhol**2d0
   fac_arr(4,2)=1.264911064067351732799557d0/rhol**2d0
   lxyz_arr(1,1,3)=1 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=3 ; lxyz_arr(2,2,3)=0 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=1 ; lxyz_arr(2,3,3)=2 ; lxyz_arr(3,3,3)=1
   lxyz_arr(1,4,3)=1 ; lxyz_arr(2,4,3)=0 ; lxyz_arr(3,4,3)=3
   fac_arr(1,3)=-2.529822128134703465599115d0
   fac_arr(2,3)=-0.3162277660168379331998894d0/rhol**2d0
   fac_arr(3,3)=-0.3162277660168379331998894d0/rhol**2d0
   fac_arr(4,3)=1.264911064067351732799557d0/rhol**2d0
else if (l.eq.4 .and. i.eq.1 .and. m.eq.2) then
   nterm_arr(1)=4
   nterm_arr(2)=6
   nterm_arr(3)=4
   lxyz_arr(1,1,1)=1 ; lxyz_arr(2,1,1)=1 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=3 ; lxyz_arr(2,2,1)=1 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=3 ; lxyz_arr(3,3,1)=0
   lxyz_arr(1,4,1)=1 ; lxyz_arr(2,4,1)=1 ; lxyz_arr(3,4,1)=2
   fac_arr(1,1)=0.6324555320336758663997787d0
   fac_arr(2,1)=-0.3162277660168379331998894d0/rhol**2d0
   fac_arr(3,1)=-0.3162277660168379331998894d0/rhol**2d0
   fac_arr(4,1)=1.264911064067351732799557d0/rhol**2d0
   lxyz_arr(1,1,2)=2 ; lxyz_arr(2,1,2)=0 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=0 ; lxyz_arr(2,2,2)=2 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=0 ; lxyz_arr(3,3,2)=2
   lxyz_arr(1,4,2)=2 ; lxyz_arr(2,4,2)=2 ; lxyz_arr(3,4,2)=0
   lxyz_arr(1,5,2)=0 ; lxyz_arr(2,5,2)=4 ; lxyz_arr(3,5,2)=0
   lxyz_arr(1,6,2)=0 ; lxyz_arr(2,6,2)=2 ; lxyz_arr(3,6,2)=2
   fac_arr(1,2)=0.3162277660168379331998894d0
   fac_arr(2,2)=0.9486832980505137995996681d0
   fac_arr(3,2)=-1.264911064067351732799557d0
   fac_arr(4,2)=-0.3162277660168379331998894d0/rhol**2d0
   fac_arr(5,2)=-0.3162277660168379331998894d0/rhol**2d0
   fac_arr(6,2)=1.264911064067351732799557d0/rhol**2d0
   lxyz_arr(1,1,3)=0 ; lxyz_arr(2,1,3)=1 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=2 ; lxyz_arr(2,2,3)=1 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=0 ; lxyz_arr(2,3,3)=3 ; lxyz_arr(3,3,3)=1
   lxyz_arr(1,4,3)=0 ; lxyz_arr(2,4,3)=1 ; lxyz_arr(3,4,3)=3
   fac_arr(1,3)=-2.529822128134703465599115d0
   fac_arr(2,3)=-0.3162277660168379331998894d0/rhol**2d0
   fac_arr(3,3)=-0.3162277660168379331998894d0/rhol**2d0
   fac_arr(4,3)=1.264911064067351732799557d0/rhol**2d0
else if (l.eq.4 .and. i.eq.1 .and. m.eq.3) then
   nterm_arr(1)=4
   nterm_arr(2)=4
   nterm_arr(3)=6
   lxyz_arr(1,1,1)=1 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=1
   lxyz_arr(1,2,1)=3 ; lxyz_arr(2,2,1)=0 ; lxyz_arr(3,2,1)=1
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=2 ; lxyz_arr(3,3,1)=1
   lxyz_arr(1,4,1)=1 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=3
   fac_arr(1,1)=1.549193338482966754071706d0
   fac_arr(2,1)=-0.7745966692414833770358531d0/rhol**2d0
   fac_arr(3,1)=-0.7745966692414833770358531d0/rhol**2d0
   fac_arr(4,1)=0.5163977794943222513572354d0/rhol**2d0
   lxyz_arr(1,1,2)=0 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=1
   lxyz_arr(1,2,2)=2 ; lxyz_arr(2,2,2)=1 ; lxyz_arr(3,2,2)=1
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=3 ; lxyz_arr(3,3,2)=1
   lxyz_arr(1,4,2)=0 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=3
   fac_arr(1,2)=1.549193338482966754071706d0
   fac_arr(2,2)=-0.7745966692414833770358531d0/rhol**2d0
   fac_arr(3,2)=-0.7745966692414833770358531d0/rhol**2d0
   fac_arr(4,2)=0.5163977794943222513572354d0/rhol**2d0
   lxyz_arr(1,1,3)=2 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=0
   lxyz_arr(1,2,3)=0 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=0
   lxyz_arr(1,3,3)=0 ; lxyz_arr(2,3,3)=0 ; lxyz_arr(3,3,3)=2
   lxyz_arr(1,4,3)=2 ; lxyz_arr(2,4,3)=0 ; lxyz_arr(3,4,3)=2
   lxyz_arr(1,5,3)=0 ; lxyz_arr(2,5,3)=2 ; lxyz_arr(3,5,3)=2
   lxyz_arr(1,6,3)=0 ; lxyz_arr(2,6,3)=0 ; lxyz_arr(3,6,3)=4
   fac_arr(1,3)=0.7745966692414833770358531d0
   fac_arr(2,3)=0.7745966692414833770358531d0
   fac_arr(3,3)=-1.549193338482966754071706d0
   fac_arr(4,3)=-0.7745966692414833770358531d0/rhol**2d0
   fac_arr(5,3)=-0.7745966692414833770358531d0/rhol**2d0
   fac_arr(6,3)=0.5163977794943222513572354d0/rhol**2d0
else if (l.eq.4 .and. i.eq.1 .and. m.eq.4) then
   nterm_arr(1)=4
   nterm_arr(2)=3
   nterm_arr(3)=2
   lxyz_arr(1,1,1)=2 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=0 ; lxyz_arr(2,2,1)=2 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=4 ; lxyz_arr(2,3,1)=0 ; lxyz_arr(3,3,1)=0
   lxyz_arr(1,4,1)=2 ; lxyz_arr(2,4,1)=2 ; lxyz_arr(3,4,1)=0
   fac_arr(1,1)=1.224744871391589049098642d0
   fac_arr(2,1)=-1.224744871391589049098642d0
   fac_arr(3,1)=-0.408248290463863016366214d0/rhol**2d0
   fac_arr(4,1)=1.224744871391589049098642d0/rhol**2d0
   lxyz_arr(1,1,2)=1 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=3 ; lxyz_arr(2,2,2)=1 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=1 ; lxyz_arr(2,3,2)=3 ; lxyz_arr(3,3,2)=0
   fac_arr(1,2)=-2.449489742783178098197284d0
   fac_arr(2,2)=-0.408248290463863016366214d0/rhol**2d0
   fac_arr(3,2)=1.224744871391589049098642d0/rhol**2d0
   lxyz_arr(1,1,3)=3 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=1 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=1
   fac_arr(1,3)=-0.408248290463863016366214d0/rhol**2d0
   fac_arr(2,3)=1.224744871391589049098642d0/rhol**2d0
else if (l.eq.4 .and. i.eq.1 .and. m.eq.5) then
   nterm_arr(1)=3
   nterm_arr(2)=4
   nterm_arr(3)=2
   lxyz_arr(1,1,1)=1 ; lxyz_arr(2,1,1)=1 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=3 ; lxyz_arr(2,2,1)=1 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=3 ; lxyz_arr(3,3,1)=0
   fac_arr(1,1)=-2.449489742783178098197284d0
   fac_arr(2,1)=1.224744871391589049098642d0/rhol**2d0
   fac_arr(3,1)=-0.408248290463863016366214d0/rhol**2d0
   lxyz_arr(1,1,2)=2 ; lxyz_arr(2,1,2)=0 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=0 ; lxyz_arr(2,2,2)=2 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=2 ; lxyz_arr(2,3,2)=2 ; lxyz_arr(3,3,2)=0
   lxyz_arr(1,4,2)=0 ; lxyz_arr(2,4,2)=4 ; lxyz_arr(3,4,2)=0
   fac_arr(1,2)=-1.224744871391589049098642d0
   fac_arr(2,2)=1.224744871391589049098642d0
   fac_arr(3,2)=1.224744871391589049098642d0/rhol**2d0
   fac_arr(4,2)=-0.408248290463863016366214d0/rhol**2d0
   lxyz_arr(1,1,3)=2 ; lxyz_arr(2,1,3)=1 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=0 ; lxyz_arr(2,2,3)=3 ; lxyz_arr(3,2,3)=1
   fac_arr(1,3)=1.224744871391589049098642d0/rhol**2d0
   fac_arr(2,3)=-0.408248290463863016366214d0/rhol**2d0
else if (l.eq.4 .and. i.eq.1 .and. m.eq.6) then
   nterm_arr(1)=3
   nterm_arr(2)=3
   nterm_arr(3)=4
   lxyz_arr(1,1,1)=1 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=1
   lxyz_arr(1,2,1)=3 ; lxyz_arr(2,2,1)=0 ; lxyz_arr(3,2,1)=1
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=2 ; lxyz_arr(3,3,1)=1
   fac_arr(1,1)=2.d0
   fac_arr(2,1)=-1.d0/rhol**2d0
   fac_arr(3,1)=rhol**(-2)
   lxyz_arr(1,1,2)=0 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=1
   lxyz_arr(1,2,2)=2 ; lxyz_arr(2,2,2)=1 ; lxyz_arr(3,2,2)=1
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=3 ; lxyz_arr(3,3,2)=1
   fac_arr(1,2)=-2.d0
   fac_arr(2,2)=-1.d0/rhol**2d0
   fac_arr(3,2)=rhol**(-2)
   lxyz_arr(1,1,3)=2 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=0
   lxyz_arr(1,2,3)=0 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=0
   lxyz_arr(1,3,3)=2 ; lxyz_arr(2,3,3)=0 ; lxyz_arr(3,3,3)=2
   lxyz_arr(1,4,3)=0 ; lxyz_arr(2,4,3)=2 ; lxyz_arr(3,4,3)=2
   fac_arr(1,3)=1.d0
   fac_arr(2,3)=-1.d0
   fac_arr(3,3)=-1.d0/rhol**2d0
   fac_arr(4,3)=rhol**(-2)
else if (l.eq.4 .and. i.eq.1 .and. m.eq.7) then
   nterm_arr(1)=2
   nterm_arr(2)=2
   nterm_arr(3)=2
   lxyz_arr(1,1,1)=0 ; lxyz_arr(2,1,1)=1 ; lxyz_arr(3,1,1)=1
   lxyz_arr(1,2,1)=2 ; lxyz_arr(2,2,1)=1 ; lxyz_arr(3,2,1)=1
   fac_arr(1,1)=2.d0
   fac_arr(2,1)=-2.d0/rhol**2d0
   lxyz_arr(1,1,2)=1 ; lxyz_arr(2,1,2)=0 ; lxyz_arr(3,1,2)=1
   lxyz_arr(1,2,2)=1 ; lxyz_arr(2,2,2)=2 ; lxyz_arr(3,2,2)=1
   fac_arr(1,2)=2.d0
   fac_arr(2,2)=-2.d0/rhol**2d0
   lxyz_arr(1,1,3)=1 ; lxyz_arr(2,1,3)=1 ; lxyz_arr(3,1,3)=0
   lxyz_arr(1,2,3)=1 ; lxyz_arr(2,2,3)=1 ; lxyz_arr(3,2,3)=2
   fac_arr(1,3)=2.d0
   fac_arr(2,3)=-2.d0/rhol**2d0
else if (l.eq.4 .and. i.eq.2 .and. m.eq.1) then
   nterm_arr(1)=12
   nterm_arr(2)=9
   nterm_arr(3)=9
   lxyz_arr(1,1,1)=4 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=2 ; lxyz_arr(2,2,1)=2 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=0 ; lxyz_arr(2,3,1)=4 ; lxyz_arr(3,3,1)=0
   lxyz_arr(1,4,1)=2 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=2
   lxyz_arr(1,5,1)=0 ; lxyz_arr(2,5,1)=2 ; lxyz_arr(3,5,1)=2
   lxyz_arr(1,6,1)=0 ; lxyz_arr(2,6,1)=0 ; lxyz_arr(3,6,1)=4
   lxyz_arr(1,7,1)=6 ; lxyz_arr(2,7,1)=0 ; lxyz_arr(3,7,1)=0
   lxyz_arr(1,8,1)=4 ; lxyz_arr(2,8,1)=2 ; lxyz_arr(3,8,1)=0
   lxyz_arr(1,9,1)=2 ; lxyz_arr(2,9,1)=4 ; lxyz_arr(3,9,1)=0
   lxyz_arr(1,10,1)=4 ; lxyz_arr(2,10,1)=0 ; lxyz_arr(3,10,1)=2
   lxyz_arr(1,11,1)=2 ; lxyz_arr(2,11,1)=2 ; lxyz_arr(3,11,1)=2
   lxyz_arr(1,12,1)=2 ; lxyz_arr(2,12,1)=0 ; lxyz_arr(3,12,1)=4
   fac_arr(1,1)=0.3178208630818641051489253d0
   fac_arr(2,1)=0.3813850356982369261787104d0
   fac_arr(3,1)=0.06356417261637282102978506d0
   fac_arr(4,1)=-0.5720775535473553892680656d0
   fac_arr(5,1)=-0.1906925178491184630893552d0
   fac_arr(6,1)=-0.2542566904654912841191402d0
   fac_arr(7,1)=-0.06356417261637282102978506d0/rhol**2d0
   fac_arr(8,1)=-0.1271283452327456420595701d0/rhol**2d0
   fac_arr(9,1)=-0.06356417261637282102978506d0/rhol**2d0
   fac_arr(10,1)=0.1906925178491184630893552d0/rhol**2d0
   fac_arr(11,1)=0.1906925178491184630893552d0/rhol**2d0
   fac_arr(12,1)=0.2542566904654912841191402d0/rhol**2d0
   lxyz_arr(1,1,2)=3 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=1 ; lxyz_arr(2,2,2)=3 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=1 ; lxyz_arr(2,3,2)=1 ; lxyz_arr(3,3,2)=2
   lxyz_arr(1,4,2)=5 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=0
   lxyz_arr(1,5,2)=3 ; lxyz_arr(2,5,2)=3 ; lxyz_arr(3,5,2)=0
   lxyz_arr(1,6,2)=1 ; lxyz_arr(2,6,2)=5 ; lxyz_arr(3,6,2)=0
   lxyz_arr(1,7,2)=3 ; lxyz_arr(2,7,2)=1 ; lxyz_arr(3,7,2)=2
   lxyz_arr(1,8,2)=1 ; lxyz_arr(2,8,2)=3 ; lxyz_arr(3,8,2)=2
   lxyz_arr(1,9,2)=1 ; lxyz_arr(2,9,2)=1 ; lxyz_arr(3,9,2)=4
   fac_arr(1,2)=0.2542566904654912841191402d0
   fac_arr(2,2)=0.2542566904654912841191402d0
   fac_arr(3,2)=-0.3813850356982369261787104d0
   fac_arr(4,2)=-0.06356417261637282102978506d0/rhol**2d0
   fac_arr(5,2)=-0.1271283452327456420595701d0/rhol**2d0
   fac_arr(6,2)=-0.06356417261637282102978506d0/rhol**2d0
   fac_arr(7,2)=0.1906925178491184630893552d0/rhol**2d0
   fac_arr(8,2)=0.1906925178491184630893552d0/rhol**2d0
   fac_arr(9,2)=0.2542566904654912841191402d0/rhol**2d0
   lxyz_arr(1,1,3)=3 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=1 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=1 ; lxyz_arr(2,3,3)=0 ; lxyz_arr(3,3,3)=3
   lxyz_arr(1,4,3)=5 ; lxyz_arr(2,4,3)=0 ; lxyz_arr(3,4,3)=1
   lxyz_arr(1,5,3)=3 ; lxyz_arr(2,5,3)=2 ; lxyz_arr(3,5,3)=1
   lxyz_arr(1,6,3)=1 ; lxyz_arr(2,6,3)=4 ; lxyz_arr(3,6,3)=1
   lxyz_arr(1,7,3)=3 ; lxyz_arr(2,7,3)=0 ; lxyz_arr(3,7,3)=3
   lxyz_arr(1,8,3)=1 ; lxyz_arr(2,8,3)=2 ; lxyz_arr(3,8,3)=3
   lxyz_arr(1,9,3)=1 ; lxyz_arr(2,9,3)=0 ; lxyz_arr(3,9,3)=5
   fac_arr(1,3)=-0.3813850356982369261787104d0
   fac_arr(2,3)=-0.3813850356982369261787104d0
   fac_arr(3,3)=-1.017026761861965136476561d0
   fac_arr(4,3)=-0.06356417261637282102978506d0/rhol**2d0
   fac_arr(5,3)=-0.1271283452327456420595701d0/rhol**2d0
   fac_arr(6,3)=-0.06356417261637282102978506d0/rhol**2d0
   fac_arr(7,3)=0.1906925178491184630893552d0/rhol**2d0
   fac_arr(8,3)=0.1906925178491184630893552d0/rhol**2d0
   fac_arr(9,3)=0.2542566904654912841191402d0/rhol**2d0
else if (l.eq.4 .and. i.eq.2 .and. m.eq.2) then
   nterm_arr(1)=9
   nterm_arr(2)=12
   nterm_arr(3)=9
   lxyz_arr(1,1,1)=3 ; lxyz_arr(2,1,1)=1 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=1 ; lxyz_arr(2,2,1)=3 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=1 ; lxyz_arr(3,3,1)=2
   lxyz_arr(1,4,1)=5 ; lxyz_arr(2,4,1)=1 ; lxyz_arr(3,4,1)=0
   lxyz_arr(1,5,1)=3 ; lxyz_arr(2,5,1)=3 ; lxyz_arr(3,5,1)=0
   lxyz_arr(1,6,1)=1 ; lxyz_arr(2,6,1)=5 ; lxyz_arr(3,6,1)=0
   lxyz_arr(1,7,1)=3 ; lxyz_arr(2,7,1)=1 ; lxyz_arr(3,7,1)=2
   lxyz_arr(1,8,1)=1 ; lxyz_arr(2,8,1)=3 ; lxyz_arr(3,8,1)=2
   lxyz_arr(1,9,1)=1 ; lxyz_arr(2,9,1)=1 ; lxyz_arr(3,9,1)=4
   fac_arr(1,1)=0.2542566904654912841191402d0
   fac_arr(2,1)=0.2542566904654912841191402d0
   fac_arr(3,1)=-0.3813850356982369261787104d0
   fac_arr(4,1)=-0.06356417261637282102978506d0/rhol**2d0
   fac_arr(5,1)=-0.1271283452327456420595701d0/rhol**2d0
   fac_arr(6,1)=-0.06356417261637282102978506d0/rhol**2d0
   fac_arr(7,1)=0.1906925178491184630893552d0/rhol**2d0
   fac_arr(8,1)=0.1906925178491184630893552d0/rhol**2d0
   fac_arr(9,1)=0.2542566904654912841191402d0/rhol**2d0
   lxyz_arr(1,1,2)=4 ; lxyz_arr(2,1,2)=0 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=2 ; lxyz_arr(2,2,2)=2 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=4 ; lxyz_arr(3,3,2)=0
   lxyz_arr(1,4,2)=2 ; lxyz_arr(2,4,2)=0 ; lxyz_arr(3,4,2)=2
   lxyz_arr(1,5,2)=0 ; lxyz_arr(2,5,2)=2 ; lxyz_arr(3,5,2)=2
   lxyz_arr(1,6,2)=0 ; lxyz_arr(2,6,2)=0 ; lxyz_arr(3,6,2)=4
   lxyz_arr(1,7,2)=4 ; lxyz_arr(2,7,2)=2 ; lxyz_arr(3,7,2)=0
   lxyz_arr(1,8,2)=2 ; lxyz_arr(2,8,2)=4 ; lxyz_arr(3,8,2)=0
   lxyz_arr(1,9,2)=0 ; lxyz_arr(2,9,2)=6 ; lxyz_arr(3,9,2)=0
   lxyz_arr(1,10,2)=2 ; lxyz_arr(2,10,2)=2 ; lxyz_arr(3,10,2)=2
   lxyz_arr(1,11,2)=0 ; lxyz_arr(2,11,2)=4 ; lxyz_arr(3,11,2)=2
   lxyz_arr(1,12,2)=0 ; lxyz_arr(2,12,2)=2 ; lxyz_arr(3,12,2)=4
   fac_arr(1,2)=0.06356417261637282102978506d0
   fac_arr(2,2)=0.3813850356982369261787104d0
   fac_arr(3,2)=0.3178208630818641051489253d0
   fac_arr(4,2)=-0.1906925178491184630893552d0
   fac_arr(5,2)=-0.5720775535473553892680656d0
   fac_arr(6,2)=-0.2542566904654912841191402d0
   fac_arr(7,2)=-0.06356417261637282102978506d0/rhol**2d0
   fac_arr(8,2)=-0.1271283452327456420595701d0/rhol**2d0
   fac_arr(9,2)=-0.06356417261637282102978506d0/rhol**2d0
   fac_arr(10,2)=0.1906925178491184630893552d0/rhol**2d0
   fac_arr(11,2)=0.1906925178491184630893552d0/rhol**2d0
   fac_arr(12,2)=0.2542566904654912841191402d0/rhol**2d0
   lxyz_arr(1,1,3)=2 ; lxyz_arr(2,1,3)=1 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=0 ; lxyz_arr(2,2,3)=3 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=0 ; lxyz_arr(2,3,3)=1 ; lxyz_arr(3,3,3)=3
   lxyz_arr(1,4,3)=4 ; lxyz_arr(2,4,3)=1 ; lxyz_arr(3,4,3)=1
   lxyz_arr(1,5,3)=2 ; lxyz_arr(2,5,3)=3 ; lxyz_arr(3,5,3)=1
   lxyz_arr(1,6,3)=0 ; lxyz_arr(2,6,3)=5 ; lxyz_arr(3,6,3)=1
   lxyz_arr(1,7,3)=2 ; lxyz_arr(2,7,3)=1 ; lxyz_arr(3,7,3)=3
   lxyz_arr(1,8,3)=0 ; lxyz_arr(2,8,3)=3 ; lxyz_arr(3,8,3)=3
   lxyz_arr(1,9,3)=0 ; lxyz_arr(2,9,3)=1 ; lxyz_arr(3,9,3)=5
   fac_arr(1,3)=-0.3813850356982369261787104d0
   fac_arr(2,3)=-0.3813850356982369261787104d0
   fac_arr(3,3)=-1.017026761861965136476561d0
   fac_arr(4,3)=-0.06356417261637282102978506d0/rhol**2d0
   fac_arr(5,3)=-0.1271283452327456420595701d0/rhol**2d0
   fac_arr(6,3)=-0.06356417261637282102978506d0/rhol**2d0
   fac_arr(7,3)=0.1906925178491184630893552d0/rhol**2d0
   fac_arr(8,3)=0.1906925178491184630893552d0/rhol**2d0
   fac_arr(9,3)=0.2542566904654912841191402d0/rhol**2d0
else if (l.eq.4 .and. i.eq.2 .and. m.eq.3) then
   nterm_arr(1)=9
   nterm_arr(2)=9
   nterm_arr(3)=12
   lxyz_arr(1,1,1)=3 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=1
   lxyz_arr(1,2,1)=1 ; lxyz_arr(2,2,1)=2 ; lxyz_arr(3,2,1)=1
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=0 ; lxyz_arr(3,3,1)=3
   lxyz_arr(1,4,1)=5 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=1
   lxyz_arr(1,5,1)=3 ; lxyz_arr(2,5,1)=2 ; lxyz_arr(3,5,1)=1
   lxyz_arr(1,6,1)=1 ; lxyz_arr(2,6,1)=4 ; lxyz_arr(3,6,1)=1
   lxyz_arr(1,7,1)=3 ; lxyz_arr(2,7,1)=0 ; lxyz_arr(3,7,1)=3
   lxyz_arr(1,8,1)=1 ; lxyz_arr(2,8,1)=2 ; lxyz_arr(3,8,1)=3
   lxyz_arr(1,9,1)=1 ; lxyz_arr(2,9,1)=0 ; lxyz_arr(3,9,1)=5
   fac_arr(1,1)=0.6227991553292183767329405d0
   fac_arr(2,1)=0.6227991553292183767329405d0
   fac_arr(3,1)=0.1037998592215363961221568d0
   fac_arr(4,1)=-0.1556997888323045941832351d0/rhol**2d0
   fac_arr(5,1)=-0.3113995776646091883664703d0/rhol**2d0
   fac_arr(6,1)=-0.1556997888323045941832351d0/rhol**2d0
   fac_arr(7,1)=-0.05189992961076819806107838d0/rhol**2d0
   fac_arr(8,1)=-0.05189992961076819806107838d0/rhol**2d0
   fac_arr(9,1)=0.1037998592215363961221568d0/rhol**2d0
   lxyz_arr(1,1,2)=2 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=1
   lxyz_arr(1,2,2)=0 ; lxyz_arr(2,2,2)=3 ; lxyz_arr(3,2,2)=1
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=1 ; lxyz_arr(3,3,2)=3
   lxyz_arr(1,4,2)=4 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=1
   lxyz_arr(1,5,2)=2 ; lxyz_arr(2,5,2)=3 ; lxyz_arr(3,5,2)=1
   lxyz_arr(1,6,2)=0 ; lxyz_arr(2,6,2)=5 ; lxyz_arr(3,6,2)=1
   lxyz_arr(1,7,2)=2 ; lxyz_arr(2,7,2)=1 ; lxyz_arr(3,7,2)=3
   lxyz_arr(1,8,2)=0 ; lxyz_arr(2,8,2)=3 ; lxyz_arr(3,8,2)=3
   lxyz_arr(1,9,2)=0 ; lxyz_arr(2,9,2)=1 ; lxyz_arr(3,9,2)=5
   fac_arr(1,2)=0.6227991553292183767329405d0
   fac_arr(2,2)=0.6227991553292183767329405d0
   fac_arr(3,2)=0.1037998592215363961221568d0
   fac_arr(4,2)=-0.1556997888323045941832351d0/rhol**2d0
   fac_arr(5,2)=-0.3113995776646091883664703d0/rhol**2d0
   fac_arr(6,2)=-0.1556997888323045941832351d0/rhol**2d0
   fac_arr(7,2)=-0.05189992961076819806107838d0/rhol**2d0
   fac_arr(8,2)=-0.05189992961076819806107838d0/rhol**2d0
   fac_arr(9,2)=0.1037998592215363961221568d0/rhol**2d0
   lxyz_arr(1,1,3)=4 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=0
   lxyz_arr(1,2,3)=2 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=0
   lxyz_arr(1,3,3)=0 ; lxyz_arr(2,3,3)=4 ; lxyz_arr(3,3,3)=0
   lxyz_arr(1,4,3)=2 ; lxyz_arr(2,4,3)=0 ; lxyz_arr(3,4,3)=2
   lxyz_arr(1,5,3)=0 ; lxyz_arr(2,5,3)=2 ; lxyz_arr(3,5,3)=2
   lxyz_arr(1,6,3)=0 ; lxyz_arr(2,6,3)=0 ; lxyz_arr(3,6,3)=4
   lxyz_arr(1,7,3)=4 ; lxyz_arr(2,7,3)=0 ; lxyz_arr(3,7,3)=2
   lxyz_arr(1,8,3)=2 ; lxyz_arr(2,8,3)=2 ; lxyz_arr(3,8,3)=2
   lxyz_arr(1,9,3)=0 ; lxyz_arr(2,9,3)=4 ; lxyz_arr(3,9,3)=2
   lxyz_arr(1,10,3)=2 ; lxyz_arr(2,10,3)=0 ; lxyz_arr(3,10,3)=4
   lxyz_arr(1,11,3)=0 ; lxyz_arr(2,11,3)=2 ; lxyz_arr(3,11,3)=4
   lxyz_arr(1,12,3)=0 ; lxyz_arr(2,12,3)=0 ; lxyz_arr(3,12,3)=6
   fac_arr(1,3)=0.1556997888323045941832351d0
   fac_arr(2,3)=0.3113995776646091883664703d0
   fac_arr(3,3)=0.1556997888323045941832351d0
   fac_arr(4,3)=0.1556997888323045941832351d0
   fac_arr(5,3)=0.1556997888323045941832351d0
   fac_arr(6,3)=-0.5189992961076819806107838d0
   fac_arr(7,3)=-0.1556997888323045941832351d0/rhol**2d0
   fac_arr(8,3)=-0.3113995776646091883664703d0/rhol**2d0
   fac_arr(9,3)=-0.1556997888323045941832351d0/rhol**2d0
   fac_arr(10,3)=-0.05189992961076819806107838d0/rhol**2d0
   fac_arr(11,3)=-0.05189992961076819806107838d0/rhol**2d0
   fac_arr(12,3)=0.1037998592215363961221568d0/rhol**2d0
else if (l.eq.4 .and. i.eq.2 .and. m.eq.4) then
   nterm_arr(1)=10
   nterm_arr(2)=8
   nterm_arr(3)=7
   lxyz_arr(1,1,1)=4 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=2 ; lxyz_arr(2,2,1)=2 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=0 ; lxyz_arr(2,3,1)=4 ; lxyz_arr(3,3,1)=0
   lxyz_arr(1,4,1)=2 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=2
   lxyz_arr(1,5,1)=0 ; lxyz_arr(2,5,1)=2 ; lxyz_arr(3,5,1)=2
   lxyz_arr(1,6,1)=6 ; lxyz_arr(2,6,1)=0 ; lxyz_arr(3,6,1)=0
   lxyz_arr(1,7,1)=4 ; lxyz_arr(2,7,1)=2 ; lxyz_arr(3,7,1)=0
   lxyz_arr(1,8,1)=2 ; lxyz_arr(2,8,1)=4 ; lxyz_arr(3,8,1)=0
   lxyz_arr(1,9,1)=4 ; lxyz_arr(2,9,1)=0 ; lxyz_arr(3,9,1)=2
   lxyz_arr(1,10,1)=2 ; lxyz_arr(2,10,1)=2 ; lxyz_arr(3,10,1)=2
   fac_arr(1,1)=0.4103049699311091091141355d0
   fac_arr(2,1)=-0.4923659639173309309369626d0
   fac_arr(3,1)=-0.2461829819586654654684813d0
   fac_arr(4,1)=0.2461829819586654654684813d0
   fac_arr(5,1)=-0.2461829819586654654684813d0
   fac_arr(6,1)=-0.08206099398622182182282711d0/rhol**2d0
   fac_arr(7,1)=0.1641219879724436436456542d0/rhol**2d0
   fac_arr(8,1)=0.2461829819586654654684813d0/rhol**2d0
   fac_arr(9,1)=-0.08206099398622182182282711d0/rhol**2d0
   fac_arr(10,1)=0.2461829819586654654684813d0/rhol**2d0
   lxyz_arr(1,1,2)=3 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=1 ; lxyz_arr(2,2,2)=3 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=1 ; lxyz_arr(2,3,2)=1 ; lxyz_arr(3,3,2)=2
   lxyz_arr(1,4,2)=5 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=0
   lxyz_arr(1,5,2)=3 ; lxyz_arr(2,5,2)=3 ; lxyz_arr(3,5,2)=0
   lxyz_arr(1,6,2)=1 ; lxyz_arr(2,6,2)=5 ; lxyz_arr(3,6,2)=0
   lxyz_arr(1,7,2)=3 ; lxyz_arr(2,7,2)=1 ; lxyz_arr(3,7,2)=2
   lxyz_arr(1,8,2)=1 ; lxyz_arr(2,8,2)=3 ; lxyz_arr(3,8,2)=2
   fac_arr(1,2)=-0.3282439759448872872913084d0
   fac_arr(2,2)=-0.9847319278346618618739253d0
   fac_arr(3,2)=-0.4923659639173309309369626d0
   fac_arr(4,2)=-0.08206099398622182182282711d0/rhol**2d0
   fac_arr(5,2)=0.1641219879724436436456542d0/rhol**2d0
   fac_arr(6,2)=0.2461829819586654654684813d0/rhol**2d0
   fac_arr(7,2)=-0.08206099398622182182282711d0/rhol**2d0
   fac_arr(8,2)=0.2461829819586654654684813d0/rhol**2d0
   lxyz_arr(1,1,3)=3 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=1 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=5 ; lxyz_arr(2,3,3)=0 ; lxyz_arr(3,3,3)=1
   lxyz_arr(1,4,3)=3 ; lxyz_arr(2,4,3)=2 ; lxyz_arr(3,4,3)=1
   lxyz_arr(1,5,3)=1 ; lxyz_arr(2,5,3)=4 ; lxyz_arr(3,5,3)=1
   lxyz_arr(1,6,3)=3 ; lxyz_arr(2,6,3)=0 ; lxyz_arr(3,6,3)=3
   lxyz_arr(1,7,3)=1 ; lxyz_arr(2,7,3)=2 ; lxyz_arr(3,7,3)=3
   fac_arr(1,3)=0.1641219879724436436456542d0
   fac_arr(2,3)=-0.4923659639173309309369626d0
   fac_arr(3,3)=-0.08206099398622182182282711d0/rhol**2d0
   fac_arr(4,3)=0.1641219879724436436456542d0/rhol**2d0
   fac_arr(5,3)=0.2461829819586654654684813d0/rhol**2d0
   fac_arr(6,3)=-0.08206099398622182182282711d0/rhol**2d0
   fac_arr(7,3)=0.2461829819586654654684813d0/rhol**2d0
else if (l.eq.4 .and. i.eq.2 .and. m.eq.5) then
   nterm_arr(1)=8
   nterm_arr(2)=10
   nterm_arr(3)=7
   lxyz_arr(1,1,1)=3 ; lxyz_arr(2,1,1)=1 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=1 ; lxyz_arr(2,2,1)=3 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=1 ; lxyz_arr(3,3,1)=2
   lxyz_arr(1,4,1)=5 ; lxyz_arr(2,4,1)=1 ; lxyz_arr(3,4,1)=0
   lxyz_arr(1,5,1)=3 ; lxyz_arr(2,5,1)=3 ; lxyz_arr(3,5,1)=0
   lxyz_arr(1,6,1)=1 ; lxyz_arr(2,6,1)=5 ; lxyz_arr(3,6,1)=0
   lxyz_arr(1,7,1)=3 ; lxyz_arr(2,7,1)=1 ; lxyz_arr(3,7,1)=2
   lxyz_arr(1,8,1)=1 ; lxyz_arr(2,8,1)=3 ; lxyz_arr(3,8,1)=2
   fac_arr(1,1)=-0.9847319278346618618739253d0
   fac_arr(2,1)=-0.3282439759448872872913084d0
   fac_arr(3,1)=-0.4923659639173309309369626d0
   fac_arr(4,1)=0.2461829819586654654684813d0/rhol**2d0
   fac_arr(5,1)=0.1641219879724436436456542d0/rhol**2d0
   fac_arr(6,1)=-0.08206099398622182182282711d0/rhol**2d0
   fac_arr(7,1)=0.2461829819586654654684813d0/rhol**2d0
   fac_arr(8,1)=-0.08206099398622182182282711d0/rhol**2d0
   lxyz_arr(1,1,2)=4 ; lxyz_arr(2,1,2)=0 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=2 ; lxyz_arr(2,2,2)=2 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=4 ; lxyz_arr(3,3,2)=0
   lxyz_arr(1,4,2)=2 ; lxyz_arr(2,4,2)=0 ; lxyz_arr(3,4,2)=2
   lxyz_arr(1,5,2)=0 ; lxyz_arr(2,5,2)=2 ; lxyz_arr(3,5,2)=2
   lxyz_arr(1,6,2)=4 ; lxyz_arr(2,6,2)=2 ; lxyz_arr(3,6,2)=0
   lxyz_arr(1,7,2)=2 ; lxyz_arr(2,7,2)=4 ; lxyz_arr(3,7,2)=0
   lxyz_arr(1,8,2)=0 ; lxyz_arr(2,8,2)=6 ; lxyz_arr(3,8,2)=0
   lxyz_arr(1,9,2)=2 ; lxyz_arr(2,9,2)=2 ; lxyz_arr(3,9,2)=2
   lxyz_arr(1,10,2)=0 ; lxyz_arr(2,10,2)=4 ; lxyz_arr(3,10,2)=2
   fac_arr(1,2)=-0.2461829819586654654684813d0
   fac_arr(2,2)=-0.4923659639173309309369626d0
   fac_arr(3,2)=0.4103049699311091091141355d0
   fac_arr(4,2)=-0.2461829819586654654684813d0
   fac_arr(5,2)=0.2461829819586654654684813d0
   fac_arr(6,2)=0.2461829819586654654684813d0/rhol**2d0
   fac_arr(7,2)=0.1641219879724436436456542d0/rhol**2d0
   fac_arr(8,2)=-0.08206099398622182182282711d0/rhol**2d0
   fac_arr(9,2)=0.2461829819586654654684813d0/rhol**2d0
   fac_arr(10,2)=-0.08206099398622182182282711d0/rhol**2d0
   lxyz_arr(1,1,3)=2 ; lxyz_arr(2,1,3)=1 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=0 ; lxyz_arr(2,2,3)=3 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=4 ; lxyz_arr(2,3,3)=1 ; lxyz_arr(3,3,3)=1
   lxyz_arr(1,4,3)=2 ; lxyz_arr(2,4,3)=3 ; lxyz_arr(3,4,3)=1
   lxyz_arr(1,5,3)=0 ; lxyz_arr(2,5,3)=5 ; lxyz_arr(3,5,3)=1
   lxyz_arr(1,6,3)=2 ; lxyz_arr(2,6,3)=1 ; lxyz_arr(3,6,3)=3
   lxyz_arr(1,7,3)=0 ; lxyz_arr(2,7,3)=3 ; lxyz_arr(3,7,3)=3
   fac_arr(1,3)=-0.4923659639173309309369626d0
   fac_arr(2,3)=0.1641219879724436436456542d0
   fac_arr(3,3)=0.2461829819586654654684813d0/rhol**2d0
   fac_arr(4,3)=0.1641219879724436436456542d0/rhol**2d0
   fac_arr(5,3)=-0.08206099398622182182282711d0/rhol**2d0
   fac_arr(6,3)=0.2461829819586654654684813d0/rhol**2d0
   fac_arr(7,3)=-0.08206099398622182182282711d0/rhol**2d0
else if (l.eq.4 .and. i.eq.2 .and. m.eq.6) then
   nterm_arr(1)=6
   nterm_arr(2)=6
   nterm_arr(3)=8
   lxyz_arr(1,1,1)=3 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=1
   lxyz_arr(1,2,1)=1 ; lxyz_arr(2,2,1)=0 ; lxyz_arr(3,2,1)=3
   lxyz_arr(1,3,1)=5 ; lxyz_arr(2,3,1)=0 ; lxyz_arr(3,3,1)=1
   lxyz_arr(1,4,1)=1 ; lxyz_arr(2,4,1)=4 ; lxyz_arr(3,4,1)=1
   lxyz_arr(1,5,1)=3 ; lxyz_arr(2,5,1)=0 ; lxyz_arr(3,5,1)=3
   lxyz_arr(1,6,1)=1 ; lxyz_arr(2,6,1)=2 ; lxyz_arr(3,6,1)=3
   fac_arr(1,1)=0.8040302522073696603914988d0
   fac_arr(2,1)=0.4020151261036848301957494d0
   fac_arr(3,1)=-0.2010075630518424150978747d0/rhol**2d0
   fac_arr(4,1)=0.2010075630518424150978747d0/rhol**2d0
   fac_arr(5,1)=-0.2010075630518424150978747d0/rhol**2d0
   fac_arr(6,1)=0.2010075630518424150978747d0/rhol**2d0
   lxyz_arr(1,1,2)=0 ; lxyz_arr(2,1,2)=3 ; lxyz_arr(3,1,2)=1
   lxyz_arr(1,2,2)=0 ; lxyz_arr(2,2,2)=1 ; lxyz_arr(3,2,2)=3
   lxyz_arr(1,3,2)=4 ; lxyz_arr(2,3,2)=1 ; lxyz_arr(3,3,2)=1
   lxyz_arr(1,4,2)=0 ; lxyz_arr(2,4,2)=5 ; lxyz_arr(3,4,2)=1
   lxyz_arr(1,5,2)=2 ; lxyz_arr(2,5,2)=1 ; lxyz_arr(3,5,2)=3
   lxyz_arr(1,6,2)=0 ; lxyz_arr(2,6,2)=3 ; lxyz_arr(3,6,2)=3
   fac_arr(1,2)=-0.8040302522073696603914988d0
   fac_arr(2,2)=-0.4020151261036848301957494d0
   fac_arr(3,2)=-0.2010075630518424150978747d0/rhol**2d0
   fac_arr(4,2)=0.2010075630518424150978747d0/rhol**2d0
   fac_arr(5,2)=-0.2010075630518424150978747d0/rhol**2d0
   fac_arr(6,2)=0.2010075630518424150978747d0/rhol**2d0
   lxyz_arr(1,1,3)=4 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=0
   lxyz_arr(1,2,3)=0 ; lxyz_arr(2,2,3)=4 ; lxyz_arr(3,2,3)=0
   lxyz_arr(1,3,3)=2 ; lxyz_arr(2,3,3)=0 ; lxyz_arr(3,3,3)=2
   lxyz_arr(1,4,3)=0 ; lxyz_arr(2,4,3)=2 ; lxyz_arr(3,4,3)=2
   lxyz_arr(1,5,3)=4 ; lxyz_arr(2,5,3)=0 ; lxyz_arr(3,5,3)=2
   lxyz_arr(1,6,3)=0 ; lxyz_arr(2,6,3)=4 ; lxyz_arr(3,6,3)=2
   lxyz_arr(1,7,3)=2 ; lxyz_arr(2,7,3)=0 ; lxyz_arr(3,7,3)=4
   lxyz_arr(1,8,3)=0 ; lxyz_arr(2,8,3)=2 ; lxyz_arr(3,8,3)=4
   fac_arr(1,3)=0.2010075630518424150978747d0
   fac_arr(2,3)=-0.2010075630518424150978747d0
   fac_arr(3,3)=0.6030226891555272452936241d0
   fac_arr(4,3)=-0.6030226891555272452936241d0
   fac_arr(5,3)=-0.2010075630518424150978747d0/rhol**2d0
   fac_arr(6,3)=0.2010075630518424150978747d0/rhol**2d0
   fac_arr(7,3)=-0.2010075630518424150978747d0/rhol**2d0
   fac_arr(8,3)=0.2010075630518424150978747d0/rhol**2d0
else if (l.eq.4 .and. i.eq.2 .and. m.eq.7) then
   nterm_arr(1)=6
   nterm_arr(2)=6
   nterm_arr(3)=6
   lxyz_arr(1,1,1)=2 ; lxyz_arr(2,1,1)=1 ; lxyz_arr(3,1,1)=1
   lxyz_arr(1,2,1)=0 ; lxyz_arr(2,2,1)=3 ; lxyz_arr(3,2,1)=1
   lxyz_arr(1,3,1)=0 ; lxyz_arr(2,3,1)=1 ; lxyz_arr(3,3,1)=3
   lxyz_arr(1,4,1)=4 ; lxyz_arr(2,4,1)=1 ; lxyz_arr(3,4,1)=1
   lxyz_arr(1,5,1)=2 ; lxyz_arr(2,5,1)=3 ; lxyz_arr(3,5,1)=1
   lxyz_arr(1,6,1)=2 ; lxyz_arr(2,6,1)=1 ; lxyz_arr(3,6,1)=3
   fac_arr(1,1)=1.206045378311054490587248d0
   fac_arr(2,1)=0.4020151261036848301957494d0
   fac_arr(3,1)=0.4020151261036848301957494d0
   fac_arr(4,1)=-0.4020151261036848301957494d0/rhol**2d0
   fac_arr(5,1)=-0.4020151261036848301957494d0/rhol**2d0
   fac_arr(6,1)=-0.4020151261036848301957494d0/rhol**2d0
   lxyz_arr(1,1,2)=3 ; lxyz_arr(2,1,2)=0 ; lxyz_arr(3,1,2)=1
   lxyz_arr(1,2,2)=1 ; lxyz_arr(2,2,2)=2 ; lxyz_arr(3,2,2)=1
   lxyz_arr(1,3,2)=1 ; lxyz_arr(2,3,2)=0 ; lxyz_arr(3,3,2)=3
   lxyz_arr(1,4,2)=3 ; lxyz_arr(2,4,2)=2 ; lxyz_arr(3,4,2)=1
   lxyz_arr(1,5,2)=1 ; lxyz_arr(2,5,2)=4 ; lxyz_arr(3,5,2)=1
   lxyz_arr(1,6,2)=1 ; lxyz_arr(2,6,2)=2 ; lxyz_arr(3,6,2)=3
   fac_arr(1,2)=0.4020151261036848301957494d0
   fac_arr(2,2)=1.206045378311054490587248d0
   fac_arr(3,2)=0.4020151261036848301957494d0
   fac_arr(4,2)=-0.4020151261036848301957494d0/rhol**2d0
   fac_arr(5,2)=-0.4020151261036848301957494d0/rhol**2d0
   fac_arr(6,2)=-0.4020151261036848301957494d0/rhol**2d0
   lxyz_arr(1,1,3)=3 ; lxyz_arr(2,1,3)=1 ; lxyz_arr(3,1,3)=0
   lxyz_arr(1,2,3)=1 ; lxyz_arr(2,2,3)=3 ; lxyz_arr(3,2,3)=0
   lxyz_arr(1,3,3)=1 ; lxyz_arr(2,3,3)=1 ; lxyz_arr(3,3,3)=2
   lxyz_arr(1,4,3)=3 ; lxyz_arr(2,4,3)=1 ; lxyz_arr(3,4,3)=2
   lxyz_arr(1,5,3)=1 ; lxyz_arr(2,5,3)=3 ; lxyz_arr(3,5,3)=2
   lxyz_arr(1,6,3)=1 ; lxyz_arr(2,6,3)=1 ; lxyz_arr(3,6,3)=4
   fac_arr(1,3)=0.4020151261036848301957494d0
   fac_arr(2,3)=0.4020151261036848301957494d0
   fac_arr(3,3)=1.206045378311054490587248d0
   fac_arr(4,3)=-0.4020151261036848301957494d0/rhol**2d0
   fac_arr(5,3)=-0.4020151261036848301957494d0/rhol**2d0
   fac_arr(6,3)=-0.4020151261036848301957494d0/rhol**2d0
else if (l.eq.4 .and. i.eq.3 .and. m.eq.1) then
   nterm_arr(1)=20
   nterm_arr(2)=16
   nterm_arr(3)=16
   lxyz_arr(1,1,1)=6 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=4 ; lxyz_arr(2,2,1)=2 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=2 ; lxyz_arr(2,3,1)=4 ; lxyz_arr(3,3,1)=0
   lxyz_arr(1,4,1)=0 ; lxyz_arr(2,4,1)=6 ; lxyz_arr(3,4,1)=0
   lxyz_arr(1,5,1)=4 ; lxyz_arr(2,5,1)=0 ; lxyz_arr(3,5,1)=2
   lxyz_arr(1,6,1)=2 ; lxyz_arr(2,6,1)=2 ; lxyz_arr(3,6,1)=2
   lxyz_arr(1,7,1)=0 ; lxyz_arr(2,7,1)=4 ; lxyz_arr(3,7,1)=2
   lxyz_arr(1,8,1)=2 ; lxyz_arr(2,8,1)=0 ; lxyz_arr(3,8,1)=4
   lxyz_arr(1,9,1)=0 ; lxyz_arr(2,9,1)=2 ; lxyz_arr(3,9,1)=4
   lxyz_arr(1,10,1)=0 ; lxyz_arr(2,10,1)=0 ; lxyz_arr(3,10,1)=6
   lxyz_arr(1,11,1)=8 ; lxyz_arr(2,11,1)=0 ; lxyz_arr(3,11,1)=0
   lxyz_arr(1,12,1)=6 ; lxyz_arr(2,12,1)=2 ; lxyz_arr(3,12,1)=0
   lxyz_arr(1,13,1)=4 ; lxyz_arr(2,13,1)=4 ; lxyz_arr(3,13,1)=0
   lxyz_arr(1,14,1)=2 ; lxyz_arr(2,14,1)=6 ; lxyz_arr(3,14,1)=0
   lxyz_arr(1,15,1)=6 ; lxyz_arr(2,15,1)=0 ; lxyz_arr(3,15,1)=2
   lxyz_arr(1,16,1)=4 ; lxyz_arr(2,16,1)=2 ; lxyz_arr(3,16,1)=2
   lxyz_arr(1,17,1)=2 ; lxyz_arr(2,17,1)=4 ; lxyz_arr(3,17,1)=2
   lxyz_arr(1,18,1)=4 ; lxyz_arr(2,18,1)=0 ; lxyz_arr(3,18,1)=4
   lxyz_arr(1,19,1)=2 ; lxyz_arr(2,19,1)=2 ; lxyz_arr(3,19,1)=4
   lxyz_arr(1,20,1)=2 ; lxyz_arr(2,20,1)=0 ; lxyz_arr(3,20,1)=6
   fac_arr(1,1)=0.06372694925323242808889581d0
   fac_arr(2,1)=0.1365577483997837744762053d0
   fac_arr(3,1)=0.08193464903987026468572318d0
   fac_arr(4,1)=0.009103849893318918298413687d0
   fac_arr(5,1)=-0.09103849893318918298413687d0
   fac_arr(6,1)=-0.1092461987198270195809642d0
   fac_arr(7,1)=-0.01820769978663783659682737d0
   fac_arr(8,1)=-0.1911808477596972842666874d0
   fac_arr(9,1)=-0.06372694925323242808889581d0
   fac_arr(10,1)=-0.03641539957327567319365475d0
   fac_arr(11,1)=-0.009103849893318918298413687d0/rhol**2d0
   fac_arr(12,1)=-0.02731154967995675489524106d0/rhol**2d0
   fac_arr(13,1)=-0.02731154967995675489524106d0/rhol**2d0
   fac_arr(14,1)=-0.009103849893318918298413687d0/rhol**2d0
   fac_arr(15,1)=0.01820769978663783659682737d0/rhol**2d0
   fac_arr(16,1)=0.03641539957327567319365475d0/rhol**2d0
   fac_arr(17,1)=0.01820769978663783659682737d0/rhol**2d0
   fac_arr(18,1)=0.06372694925323242808889581d0/rhol**2d0
   fac_arr(19,1)=0.06372694925323242808889581d0/rhol**2d0
   fac_arr(20,1)=0.03641539957327567319365475d0/rhol**2d0
   lxyz_arr(1,1,2)=5 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=3 ; lxyz_arr(2,2,2)=3 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=1 ; lxyz_arr(2,3,2)=5 ; lxyz_arr(3,3,2)=0
   lxyz_arr(1,4,2)=3 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=2
   lxyz_arr(1,5,2)=1 ; lxyz_arr(2,5,2)=3 ; lxyz_arr(3,5,2)=2
   lxyz_arr(1,6,2)=1 ; lxyz_arr(2,6,2)=1 ; lxyz_arr(3,6,2)=4
   lxyz_arr(1,7,2)=7 ; lxyz_arr(2,7,2)=1 ; lxyz_arr(3,7,2)=0
   lxyz_arr(1,8,2)=5 ; lxyz_arr(2,8,2)=3 ; lxyz_arr(3,8,2)=0
   lxyz_arr(1,9,2)=3 ; lxyz_arr(2,9,2)=5 ; lxyz_arr(3,9,2)=0
   lxyz_arr(1,10,2)=1 ; lxyz_arr(2,10,2)=7 ; lxyz_arr(3,10,2)=0
   lxyz_arr(1,11,2)=5 ; lxyz_arr(2,11,2)=1 ; lxyz_arr(3,11,2)=2
   lxyz_arr(1,12,2)=3 ; lxyz_arr(2,12,2)=3 ; lxyz_arr(3,12,2)=2
   lxyz_arr(1,13,2)=1 ; lxyz_arr(2,13,2)=5 ; lxyz_arr(3,13,2)=2
   lxyz_arr(1,14,2)=3 ; lxyz_arr(2,14,2)=1 ; lxyz_arr(3,14,2)=4
   lxyz_arr(1,15,2)=1 ; lxyz_arr(2,15,2)=3 ; lxyz_arr(3,15,2)=4
   lxyz_arr(1,16,2)=1 ; lxyz_arr(2,16,2)=1 ; lxyz_arr(3,16,2)=6
   fac_arr(1,2)=0.05462309935991350979048212d0
   fac_arr(2,2)=0.1092461987198270195809642d0
   fac_arr(3,2)=0.05462309935991350979048212d0
   fac_arr(4,2)=-0.0728307991465513463873095d0
   fac_arr(5,2)=-0.0728307991465513463873095d0
   fac_arr(6,2)=-0.1274538985064648561777916d0
   fac_arr(7,2)=-0.009103849893318918298413687d0/rhol**2d0
   fac_arr(8,2)=-0.02731154967995675489524106d0/rhol**2d0
   fac_arr(9,2)=-0.02731154967995675489524106d0/rhol**2d0
   fac_arr(10,2)=-0.009103849893318918298413687d0/rhol**2d0
   fac_arr(11,2)=0.01820769978663783659682737d0/rhol**2d0
   fac_arr(12,2)=0.03641539957327567319365475d0/rhol**2d0
   fac_arr(13,2)=0.01820769978663783659682737d0/rhol**2d0
   fac_arr(14,2)=0.06372694925323242808889581d0/rhol**2d0
   fac_arr(15,2)=0.06372694925323242808889581d0/rhol**2d0
   fac_arr(16,2)=0.03641539957327567319365475d0/rhol**2d0
   lxyz_arr(1,1,3)=5 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=3 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=1 ; lxyz_arr(2,3,3)=4 ; lxyz_arr(3,3,3)=1
   lxyz_arr(1,4,3)=3 ; lxyz_arr(2,4,3)=0 ; lxyz_arr(3,4,3)=3
   lxyz_arr(1,5,3)=1 ; lxyz_arr(2,5,3)=2 ; lxyz_arr(3,5,3)=3
   lxyz_arr(1,6,3)=1 ; lxyz_arr(2,6,3)=0 ; lxyz_arr(3,6,3)=5
   lxyz_arr(1,7,3)=7 ; lxyz_arr(2,7,3)=0 ; lxyz_arr(3,7,3)=1
   lxyz_arr(1,8,3)=5 ; lxyz_arr(2,8,3)=2 ; lxyz_arr(3,8,3)=1
   lxyz_arr(1,9,3)=3 ; lxyz_arr(2,9,3)=4 ; lxyz_arr(3,9,3)=1
   lxyz_arr(1,10,3)=1 ; lxyz_arr(2,10,3)=6 ; lxyz_arr(3,10,3)=1
   lxyz_arr(1,11,3)=5 ; lxyz_arr(2,11,3)=0 ; lxyz_arr(3,11,3)=3
   lxyz_arr(1,12,3)=3 ; lxyz_arr(2,12,3)=2 ; lxyz_arr(3,12,3)=3
   lxyz_arr(1,13,3)=1 ; lxyz_arr(2,13,3)=4 ; lxyz_arr(3,13,3)=3
   lxyz_arr(1,14,3)=3 ; lxyz_arr(2,14,3)=0 ; lxyz_arr(3,14,3)=5
   lxyz_arr(1,15,3)=1 ; lxyz_arr(2,15,3)=2 ; lxyz_arr(3,15,3)=5
   lxyz_arr(1,16,3)=1 ; lxyz_arr(2,16,3)=0 ; lxyz_arr(3,16,3)=7
   fac_arr(1,3)=-0.03641539957327567319365475d0
   fac_arr(2,3)=-0.0728307991465513463873095d0
   fac_arr(3,3)=-0.03641539957327567319365475d0
   fac_arr(4,3)=-0.2549077970129297123555832d0
   fac_arr(5,3)=-0.2549077970129297123555832d0
   fac_arr(6,3)=-0.2184923974396540391619285d0
   fac_arr(7,3)=-0.009103849893318918298413687d0/rhol**2d0
   fac_arr(8,3)=-0.02731154967995675489524106d0/rhol**2d0
   fac_arr(9,3)=-0.02731154967995675489524106d0/rhol**2d0
   fac_arr(10,3)=-0.009103849893318918298413687d0/rhol**2d0
   fac_arr(11,3)=0.01820769978663783659682737d0/rhol**2d0
   fac_arr(12,3)=0.03641539957327567319365475d0/rhol**2d0
   fac_arr(13,3)=0.01820769978663783659682737d0/rhol**2d0
   fac_arr(14,3)=0.06372694925323242808889581d0/rhol**2d0
   fac_arr(15,3)=0.06372694925323242808889581d0/rhol**2d0
   fac_arr(16,3)=0.03641539957327567319365475d0/rhol**2d0
else if (l.eq.4 .and. i.eq.3 .and. m.eq.2) then
   nterm_arr(1)=16
   nterm_arr(2)=20
   nterm_arr(3)=16
   lxyz_arr(1,1,1)=5 ; lxyz_arr(2,1,1)=1 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=3 ; lxyz_arr(2,2,1)=3 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=5 ; lxyz_arr(3,3,1)=0
   lxyz_arr(1,4,1)=3 ; lxyz_arr(2,4,1)=1 ; lxyz_arr(3,4,1)=2
   lxyz_arr(1,5,1)=1 ; lxyz_arr(2,5,1)=3 ; lxyz_arr(3,5,1)=2
   lxyz_arr(1,6,1)=1 ; lxyz_arr(2,6,1)=1 ; lxyz_arr(3,6,1)=4
   lxyz_arr(1,7,1)=7 ; lxyz_arr(2,7,1)=1 ; lxyz_arr(3,7,1)=0
   lxyz_arr(1,8,1)=5 ; lxyz_arr(2,8,1)=3 ; lxyz_arr(3,8,1)=0
   lxyz_arr(1,9,1)=3 ; lxyz_arr(2,9,1)=5 ; lxyz_arr(3,9,1)=0
   lxyz_arr(1,10,1)=1 ; lxyz_arr(2,10,1)=7 ; lxyz_arr(3,10,1)=0
   lxyz_arr(1,11,1)=5 ; lxyz_arr(2,11,1)=1 ; lxyz_arr(3,11,1)=2
   lxyz_arr(1,12,1)=3 ; lxyz_arr(2,12,1)=3 ; lxyz_arr(3,12,1)=2
   lxyz_arr(1,13,1)=1 ; lxyz_arr(2,13,1)=5 ; lxyz_arr(3,13,1)=2
   lxyz_arr(1,14,1)=3 ; lxyz_arr(2,14,1)=1 ; lxyz_arr(3,14,1)=4
   lxyz_arr(1,15,1)=1 ; lxyz_arr(2,15,1)=3 ; lxyz_arr(3,15,1)=4
   lxyz_arr(1,16,1)=1 ; lxyz_arr(2,16,1)=1 ; lxyz_arr(3,16,1)=6
   fac_arr(1,1)=0.05462309935991350979048212d0
   fac_arr(2,1)=0.1092461987198270195809642d0
   fac_arr(3,1)=0.05462309935991350979048212d0
   fac_arr(4,1)=-0.0728307991465513463873095d0
   fac_arr(5,1)=-0.0728307991465513463873095d0
   fac_arr(6,1)=-0.1274538985064648561777916d0
   fac_arr(7,1)=-0.009103849893318918298413687d0/rhol**2d0
   fac_arr(8,1)=-0.02731154967995675489524106d0/rhol**2d0
   fac_arr(9,1)=-0.02731154967995675489524106d0/rhol**2d0
   fac_arr(10,1)=-0.009103849893318918298413687d0/rhol**2d0
   fac_arr(11,1)=0.01820769978663783659682737d0/rhol**2d0
   fac_arr(12,1)=0.03641539957327567319365475d0/rhol**2d0
   fac_arr(13,1)=0.01820769978663783659682737d0/rhol**2d0
   fac_arr(14,1)=0.06372694925323242808889581d0/rhol**2d0
   fac_arr(15,1)=0.06372694925323242808889581d0/rhol**2d0
   fac_arr(16,1)=0.03641539957327567319365475d0/rhol**2d0
   lxyz_arr(1,1,2)=6 ; lxyz_arr(2,1,2)=0 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=4 ; lxyz_arr(2,2,2)=2 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=2 ; lxyz_arr(2,3,2)=4 ; lxyz_arr(3,3,2)=0
   lxyz_arr(1,4,2)=0 ; lxyz_arr(2,4,2)=6 ; lxyz_arr(3,4,2)=0
   lxyz_arr(1,5,2)=4 ; lxyz_arr(2,5,2)=0 ; lxyz_arr(3,5,2)=2
   lxyz_arr(1,6,2)=2 ; lxyz_arr(2,6,2)=2 ; lxyz_arr(3,6,2)=2
   lxyz_arr(1,7,2)=0 ; lxyz_arr(2,7,2)=4 ; lxyz_arr(3,7,2)=2
   lxyz_arr(1,8,2)=2 ; lxyz_arr(2,8,2)=0 ; lxyz_arr(3,8,2)=4
   lxyz_arr(1,9,2)=0 ; lxyz_arr(2,9,2)=2 ; lxyz_arr(3,9,2)=4
   lxyz_arr(1,10,2)=0 ; lxyz_arr(2,10,2)=0 ; lxyz_arr(3,10,2)=6
   lxyz_arr(1,11,2)=6 ; lxyz_arr(2,11,2)=2 ; lxyz_arr(3,11,2)=0
   lxyz_arr(1,12,2)=4 ; lxyz_arr(2,12,2)=4 ; lxyz_arr(3,12,2)=0
   lxyz_arr(1,13,2)=2 ; lxyz_arr(2,13,2)=6 ; lxyz_arr(3,13,2)=0
   lxyz_arr(1,14,2)=0 ; lxyz_arr(2,14,2)=8 ; lxyz_arr(3,14,2)=0
   lxyz_arr(1,15,2)=4 ; lxyz_arr(2,15,2)=2 ; lxyz_arr(3,15,2)=2
   lxyz_arr(1,16,2)=2 ; lxyz_arr(2,16,2)=4 ; lxyz_arr(3,16,2)=2
   lxyz_arr(1,17,2)=0 ; lxyz_arr(2,17,2)=6 ; lxyz_arr(3,17,2)=2
   lxyz_arr(1,18,2)=2 ; lxyz_arr(2,18,2)=2 ; lxyz_arr(3,18,2)=4
   lxyz_arr(1,19,2)=0 ; lxyz_arr(2,19,2)=4 ; lxyz_arr(3,19,2)=4
   lxyz_arr(1,20,2)=0 ; lxyz_arr(2,20,2)=2 ; lxyz_arr(3,20,2)=6
   fac_arr(1,2)=0.009103849893318918298413687d0
   fac_arr(2,2)=0.08193464903987026468572318d0
   fac_arr(3,2)=0.1365577483997837744762053d0
   fac_arr(4,2)=0.06372694925323242808889581d0
   fac_arr(5,2)=-0.01820769978663783659682737d0
   fac_arr(6,2)=-0.1092461987198270195809642d0
   fac_arr(7,2)=-0.09103849893318918298413687d0
   fac_arr(8,2)=-0.06372694925323242808889581d0
   fac_arr(9,2)=-0.1911808477596972842666874d0
   fac_arr(10,2)=-0.03641539957327567319365475d0
   fac_arr(11,2)=-0.009103849893318918298413687d0/rhol**2d0
   fac_arr(12,2)=-0.02731154967995675489524106d0/rhol**2d0
   fac_arr(13,2)=-0.02731154967995675489524106d0/rhol**2d0
   fac_arr(14,2)=-0.009103849893318918298413687d0/rhol**2d0
   fac_arr(15,2)=0.01820769978663783659682737d0/rhol**2d0
   fac_arr(16,2)=0.03641539957327567319365475d0/rhol**2d0
   fac_arr(17,2)=0.01820769978663783659682737d0/rhol**2d0
   fac_arr(18,2)=0.06372694925323242808889581d0/rhol**2d0
   fac_arr(19,2)=0.06372694925323242808889581d0/rhol**2d0
   fac_arr(20,2)=0.03641539957327567319365475d0/rhol**2d0
   lxyz_arr(1,1,3)=4 ; lxyz_arr(2,1,3)=1 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=2 ; lxyz_arr(2,2,3)=3 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=0 ; lxyz_arr(2,3,3)=5 ; lxyz_arr(3,3,3)=1
   lxyz_arr(1,4,3)=2 ; lxyz_arr(2,4,3)=1 ; lxyz_arr(3,4,3)=3
   lxyz_arr(1,5,3)=0 ; lxyz_arr(2,5,3)=3 ; lxyz_arr(3,5,3)=3
   lxyz_arr(1,6,3)=0 ; lxyz_arr(2,6,3)=1 ; lxyz_arr(3,6,3)=5
   lxyz_arr(1,7,3)=6 ; lxyz_arr(2,7,3)=1 ; lxyz_arr(3,7,3)=1
   lxyz_arr(1,8,3)=4 ; lxyz_arr(2,8,3)=3 ; lxyz_arr(3,8,3)=1
   lxyz_arr(1,9,3)=2 ; lxyz_arr(2,9,3)=5 ; lxyz_arr(3,9,3)=1
   lxyz_arr(1,10,3)=0 ; lxyz_arr(2,10,3)=7 ; lxyz_arr(3,10,3)=1
   lxyz_arr(1,11,3)=4 ; lxyz_arr(2,11,3)=1 ; lxyz_arr(3,11,3)=3
   lxyz_arr(1,12,3)=2 ; lxyz_arr(2,12,3)=3 ; lxyz_arr(3,12,3)=3
   lxyz_arr(1,13,3)=0 ; lxyz_arr(2,13,3)=5 ; lxyz_arr(3,13,3)=3
   lxyz_arr(1,14,3)=2 ; lxyz_arr(2,14,3)=1 ; lxyz_arr(3,14,3)=5
   lxyz_arr(1,15,3)=0 ; lxyz_arr(2,15,3)=3 ; lxyz_arr(3,15,3)=5
   lxyz_arr(1,16,3)=0 ; lxyz_arr(2,16,3)=1 ; lxyz_arr(3,16,3)=7
   fac_arr(1,3)=-0.03641539957327567319365475d0
   fac_arr(2,3)=-0.0728307991465513463873095d0
   fac_arr(3,3)=-0.03641539957327567319365475d0
   fac_arr(4,3)=-0.2549077970129297123555832d0
   fac_arr(5,3)=-0.2549077970129297123555832d0
   fac_arr(6,3)=-0.2184923974396540391619285d0
   fac_arr(7,3)=-0.009103849893318918298413687d0/rhol**2d0
   fac_arr(8,3)=-0.02731154967995675489524106d0/rhol**2d0
   fac_arr(9,3)=-0.02731154967995675489524106d0/rhol**2d0
   fac_arr(10,3)=-0.009103849893318918298413687d0/rhol**2d0
   fac_arr(11,3)=0.01820769978663783659682737d0/rhol**2d0
   fac_arr(12,3)=0.03641539957327567319365475d0/rhol**2d0
   fac_arr(13,3)=0.01820769978663783659682737d0/rhol**2d0
   fac_arr(14,3)=0.06372694925323242808889581d0/rhol**2d0
   fac_arr(15,3)=0.06372694925323242808889581d0/rhol**2d0
   fac_arr(16,3)=0.03641539957327567319365475d0/rhol**2d0
else if (l.eq.4 .and. i.eq.3 .and. m.eq.3) then
   nterm_arr(1)=16
   nterm_arr(2)=16
   nterm_arr(3)=20
   lxyz_arr(1,1,1)=5 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=1
   lxyz_arr(1,2,1)=3 ; lxyz_arr(2,2,1)=2 ; lxyz_arr(3,2,1)=1
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=4 ; lxyz_arr(3,3,1)=1
   lxyz_arr(1,4,1)=3 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=3
   lxyz_arr(1,5,1)=1 ; lxyz_arr(2,5,1)=2 ; lxyz_arr(3,5,1)=3
   lxyz_arr(1,6,1)=1 ; lxyz_arr(2,6,1)=0 ; lxyz_arr(3,6,1)=5
   lxyz_arr(1,7,1)=7 ; lxyz_arr(2,7,1)=0 ; lxyz_arr(3,7,1)=1
   lxyz_arr(1,8,1)=5 ; lxyz_arr(2,8,1)=2 ; lxyz_arr(3,8,1)=1
   lxyz_arr(1,9,1)=3 ; lxyz_arr(2,9,1)=4 ; lxyz_arr(3,9,1)=1
   lxyz_arr(1,10,1)=1 ; lxyz_arr(2,10,1)=6 ; lxyz_arr(3,10,1)=1
   lxyz_arr(1,11,1)=5 ; lxyz_arr(2,11,1)=0 ; lxyz_arr(3,11,1)=3
   lxyz_arr(1,12,1)=3 ; lxyz_arr(2,12,1)=2 ; lxyz_arr(3,12,1)=3
   lxyz_arr(1,13,1)=1 ; lxyz_arr(2,13,1)=4 ; lxyz_arr(3,13,1)=3
   lxyz_arr(1,14,1)=3 ; lxyz_arr(2,14,1)=0 ; lxyz_arr(3,14,1)=5
   lxyz_arr(1,15,1)=1 ; lxyz_arr(2,15,1)=2 ; lxyz_arr(3,15,1)=5
   lxyz_arr(1,16,1)=1 ; lxyz_arr(2,16,1)=0 ; lxyz_arr(3,16,1)=7
   fac_arr(1,1)=0.1337987216011345233133409d0
   fac_arr(2,1)=0.2675974432022690466266818d0
   fac_arr(3,1)=0.1337987216011345233133409d0
   fac_arr(4,1)=0.1189321969787862429451919d0
   fac_arr(5,1)=0.1189321969787862429451919d0
   fac_arr(6,1)=-0.01486652462234828036814899d0
   fac_arr(7,1)=-0.02229978693352242055222348d0/rhol**2d0
   fac_arr(8,1)=-0.06689936080056726165667044d0/rhol**2d0
   fac_arr(9,1)=-0.06689936080056726165667044d0/rhol**2d0
   fac_arr(10,1)=-0.02229978693352242055222348d0/rhol**2d0
   fac_arr(11,1)=-0.02973304924469656073629797d0/rhol**2d0
   fac_arr(12,1)=-0.05946609848939312147259594d0/rhol**2d0
   fac_arr(13,1)=-0.02973304924469656073629797d0/rhol**2d0
   fac_arr(14,1)=0.007433262311174140184074493d0/rhol**2d0
   fac_arr(15,1)=0.007433262311174140184074493d0/rhol**2d0
   fac_arr(16,1)=0.01486652462234828036814899d0/rhol**2d0
   lxyz_arr(1,1,2)=4 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=1
   lxyz_arr(1,2,2)=2 ; lxyz_arr(2,2,2)=3 ; lxyz_arr(3,2,2)=1
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=5 ; lxyz_arr(3,3,2)=1
   lxyz_arr(1,4,2)=2 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=3
   lxyz_arr(1,5,2)=0 ; lxyz_arr(2,5,2)=3 ; lxyz_arr(3,5,2)=3
   lxyz_arr(1,6,2)=0 ; lxyz_arr(2,6,2)=1 ; lxyz_arr(3,6,2)=5
   lxyz_arr(1,7,2)=6 ; lxyz_arr(2,7,2)=1 ; lxyz_arr(3,7,2)=1
   lxyz_arr(1,8,2)=4 ; lxyz_arr(2,8,2)=3 ; lxyz_arr(3,8,2)=1
   lxyz_arr(1,9,2)=2 ; lxyz_arr(2,9,2)=5 ; lxyz_arr(3,9,2)=1
   lxyz_arr(1,10,2)=0 ; lxyz_arr(2,10,2)=7 ; lxyz_arr(3,10,2)=1
   lxyz_arr(1,11,2)=4 ; lxyz_arr(2,11,2)=1 ; lxyz_arr(3,11,2)=3
   lxyz_arr(1,12,2)=2 ; lxyz_arr(2,12,2)=3 ; lxyz_arr(3,12,2)=3
   lxyz_arr(1,13,2)=0 ; lxyz_arr(2,13,2)=5 ; lxyz_arr(3,13,2)=3
   lxyz_arr(1,14,2)=2 ; lxyz_arr(2,14,2)=1 ; lxyz_arr(3,14,2)=5
   lxyz_arr(1,15,2)=0 ; lxyz_arr(2,15,2)=3 ; lxyz_arr(3,15,2)=5
   lxyz_arr(1,16,2)=0 ; lxyz_arr(2,16,2)=1 ; lxyz_arr(3,16,2)=7
   fac_arr(1,2)=0.1337987216011345233133409d0
   fac_arr(2,2)=0.2675974432022690466266818d0
   fac_arr(3,2)=0.1337987216011345233133409d0
   fac_arr(4,2)=0.1189321969787862429451919d0
   fac_arr(5,2)=0.1189321969787862429451919d0
   fac_arr(6,2)=-0.01486652462234828036814899d0
   fac_arr(7,2)=-0.02229978693352242055222348d0/rhol**2d0
   fac_arr(8,2)=-0.06689936080056726165667044d0/rhol**2d0
   fac_arr(9,2)=-0.06689936080056726165667044d0/rhol**2d0
   fac_arr(10,2)=-0.02229978693352242055222348d0/rhol**2d0
   fac_arr(11,2)=-0.02973304924469656073629797d0/rhol**2d0
   fac_arr(12,2)=-0.05946609848939312147259594d0/rhol**2d0
   fac_arr(13,2)=-0.02973304924469656073629797d0/rhol**2d0
   fac_arr(14,2)=0.007433262311174140184074493d0/rhol**2d0
   fac_arr(15,2)=0.007433262311174140184074493d0/rhol**2d0
   fac_arr(16,2)=0.01486652462234828036814899d0/rhol**2d0
   lxyz_arr(1,1,3)=6 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=0
   lxyz_arr(1,2,3)=4 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=0
   lxyz_arr(1,3,3)=2 ; lxyz_arr(2,3,3)=4 ; lxyz_arr(3,3,3)=0
   lxyz_arr(1,4,3)=0 ; lxyz_arr(2,4,3)=6 ; lxyz_arr(3,4,3)=0
   lxyz_arr(1,5,3)=4 ; lxyz_arr(2,5,3)=0 ; lxyz_arr(3,5,3)=2
   lxyz_arr(1,6,3)=2 ; lxyz_arr(2,6,3)=2 ; lxyz_arr(3,6,3)=2
   lxyz_arr(1,7,3)=0 ; lxyz_arr(2,7,3)=4 ; lxyz_arr(3,7,3)=2
   lxyz_arr(1,8,3)=2 ; lxyz_arr(2,8,3)=0 ; lxyz_arr(3,8,3)=4
   lxyz_arr(1,9,3)=0 ; lxyz_arr(2,9,3)=2 ; lxyz_arr(3,9,3)=4
   lxyz_arr(1,10,3)=0 ; lxyz_arr(2,10,3)=0 ; lxyz_arr(3,10,3)=6
   lxyz_arr(1,11,3)=6 ; lxyz_arr(2,11,3)=0 ; lxyz_arr(3,11,3)=2
   lxyz_arr(1,12,3)=4 ; lxyz_arr(2,12,3)=2 ; lxyz_arr(3,12,3)=2
   lxyz_arr(1,13,3)=2 ; lxyz_arr(2,13,3)=4 ; lxyz_arr(3,13,3)=2
   lxyz_arr(1,14,3)=0 ; lxyz_arr(2,14,3)=6 ; lxyz_arr(3,14,3)=2
   lxyz_arr(1,15,3)=4 ; lxyz_arr(2,15,3)=0 ; lxyz_arr(3,15,3)=4
   lxyz_arr(1,16,3)=2 ; lxyz_arr(2,16,3)=2 ; lxyz_arr(3,16,3)=4
   lxyz_arr(1,17,3)=0 ; lxyz_arr(2,17,3)=4 ; lxyz_arr(3,17,3)=4
   lxyz_arr(1,18,3)=2 ; lxyz_arr(2,18,3)=0 ; lxyz_arr(3,18,3)=6
   lxyz_arr(1,19,3)=0 ; lxyz_arr(2,19,3)=2 ; lxyz_arr(3,19,3)=6
   lxyz_arr(1,20,3)=0 ; lxyz_arr(2,20,3)=0 ; lxyz_arr(3,20,3)=8
   fac_arr(1,3)=0.02229978693352242055222348d0
   fac_arr(2,3)=0.06689936080056726165667044d0
   fac_arr(3,3)=0.06689936080056726165667044d0
   fac_arr(4,3)=0.02229978693352242055222348d0
   fac_arr(5,3)=0.08919914773408968220889392d0
   fac_arr(6,3)=0.1783982954681793644177878d0
   fac_arr(7,3)=0.08919914773408968220889392d0
   fac_arr(8,3)=-0.03716631155587070092037247d0
   fac_arr(9,3)=-0.03716631155587070092037247d0
   fac_arr(10,3)=-0.1040656723564379625770429d0
   fac_arr(11,3)=-0.02229978693352242055222348d0/rhol**2d0
   fac_arr(12,3)=-0.06689936080056726165667044d0/rhol**2d0
   fac_arr(13,3)=-0.06689936080056726165667044d0/rhol**2d0
   fac_arr(14,3)=-0.02229978693352242055222348d0/rhol**2d0
   fac_arr(15,3)=-0.02973304924469656073629797d0/rhol**2d0
   fac_arr(16,3)=-0.05946609848939312147259594d0/rhol**2d0
   fac_arr(17,3)=-0.02973304924469656073629797d0/rhol**2d0
   fac_arr(18,3)=0.007433262311174140184074493d0/rhol**2d0
   fac_arr(19,3)=0.007433262311174140184074493d0/rhol**2d0
   fac_arr(20,3)=0.01486652462234828036814899d0/rhol**2d0
else if (l.eq.4 .and. i.eq.3 .and. m.eq.4) then
   nterm_arr(1)=18
   nterm_arr(2)=15
   nterm_arr(3)=14
   lxyz_arr(1,1,1)=6 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=4 ; lxyz_arr(2,2,1)=2 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=2 ; lxyz_arr(2,3,1)=4 ; lxyz_arr(3,3,1)=0
   lxyz_arr(1,4,1)=0 ; lxyz_arr(2,4,1)=6 ; lxyz_arr(3,4,1)=0
   lxyz_arr(1,5,1)=4 ; lxyz_arr(2,5,1)=0 ; lxyz_arr(3,5,1)=2
   lxyz_arr(1,6,1)=2 ; lxyz_arr(2,6,1)=2 ; lxyz_arr(3,6,1)=2
   lxyz_arr(1,7,1)=0 ; lxyz_arr(2,7,1)=4 ; lxyz_arr(3,7,1)=2
   lxyz_arr(1,8,1)=2 ; lxyz_arr(2,8,1)=0 ; lxyz_arr(3,8,1)=4
   lxyz_arr(1,9,1)=0 ; lxyz_arr(2,9,1)=2 ; lxyz_arr(3,9,1)=4
   lxyz_arr(1,10,1)=8 ; lxyz_arr(2,10,1)=0 ; lxyz_arr(3,10,1)=0
   lxyz_arr(1,11,1)=6 ; lxyz_arr(2,11,1)=2 ; lxyz_arr(3,11,1)=0
   lxyz_arr(1,12,1)=4 ; lxyz_arr(2,12,1)=4 ; lxyz_arr(3,12,1)=0
   lxyz_arr(1,13,1)=2 ; lxyz_arr(2,13,1)=6 ; lxyz_arr(3,13,1)=0
   lxyz_arr(1,14,1)=6 ; lxyz_arr(2,14,1)=0 ; lxyz_arr(3,14,1)=2
   lxyz_arr(1,15,1)=4 ; lxyz_arr(2,15,1)=2 ; lxyz_arr(3,15,1)=2
   lxyz_arr(1,16,1)=2 ; lxyz_arr(2,16,1)=4 ; lxyz_arr(3,16,1)=2
   lxyz_arr(1,17,1)=4 ; lxyz_arr(2,17,1)=0 ; lxyz_arr(3,17,1)=4
   lxyz_arr(1,18,1)=2 ; lxyz_arr(2,18,1)=2 ; lxyz_arr(3,18,1)=4
   fac_arr(1,1)=0.08227113772079145865717289d0
   fac_arr(2,1)=-0.05876509837199389904083778d0
   fac_arr(3,1)=-0.1762952951159816971225133d0
   fac_arr(4,1)=-0.03525905902319633942450267d0
   fac_arr(5,1)=0.1175301967439877980816756d0
   fac_arr(6,1)=-0.1410362360927853576980107d0
   fac_arr(7,1)=-0.07051811804639267884900533d0
   fac_arr(8,1)=0.03525905902319633942450267d0
   fac_arr(9,1)=-0.03525905902319633942450267d0
   fac_arr(10,1)=-0.01175301967439877980816756d0/rhol**2d0
   fac_arr(11,1)=0.01175301967439877980816756d0/rhol**2d0
   fac_arr(12,1)=0.05876509837199389904083778d0/rhol**2d0
   fac_arr(13,1)=0.03525905902319633942450267d0/rhol**2d0
   fac_arr(14,1)=-0.02350603934879755961633511d0/rhol**2d0
   fac_arr(15,1)=0.04701207869759511923267022d0/rhol**2d0
   fac_arr(16,1)=0.07051811804639267884900533d0/rhol**2d0
   fac_arr(17,1)=-0.01175301967439877980816756d0/rhol**2d0
   fac_arr(18,1)=0.03525905902319633942450267d0/rhol**2d0
   lxyz_arr(1,1,2)=5 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=3 ; lxyz_arr(2,2,2)=3 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=1 ; lxyz_arr(2,3,2)=5 ; lxyz_arr(3,3,2)=0
   lxyz_arr(1,4,2)=3 ; lxyz_arr(2,4,2)=1 ; lxyz_arr(3,4,2)=2
   lxyz_arr(1,5,2)=1 ; lxyz_arr(2,5,2)=3 ; lxyz_arr(3,5,2)=2
   lxyz_arr(1,6,2)=1 ; lxyz_arr(2,6,2)=1 ; lxyz_arr(3,6,2)=4
   lxyz_arr(1,7,2)=7 ; lxyz_arr(2,7,2)=1 ; lxyz_arr(3,7,2)=0
   lxyz_arr(1,8,2)=5 ; lxyz_arr(2,8,2)=3 ; lxyz_arr(3,8,2)=0
   lxyz_arr(1,9,2)=3 ; lxyz_arr(2,9,2)=5 ; lxyz_arr(3,9,2)=0
   lxyz_arr(1,10,2)=1 ; lxyz_arr(2,10,2)=7 ; lxyz_arr(3,10,2)=0
   lxyz_arr(1,11,2)=5 ; lxyz_arr(2,11,2)=1 ; lxyz_arr(3,11,2)=2
   lxyz_arr(1,12,2)=3 ; lxyz_arr(2,12,2)=3 ; lxyz_arr(3,12,2)=2
   lxyz_arr(1,13,2)=1 ; lxyz_arr(2,13,2)=5 ; lxyz_arr(3,13,2)=2
   lxyz_arr(1,14,2)=3 ; lxyz_arr(2,14,2)=1 ; lxyz_arr(3,14,2)=4
   lxyz_arr(1,15,2)=1 ; lxyz_arr(2,15,2)=3 ; lxyz_arr(3,15,2)=4
   fac_arr(1,2)=-0.02350603934879755961633511d0
   fac_arr(2,2)=-0.2350603934879755961633511d0
   fac_arr(3,2)=-0.211554354139178036547016d0
   fac_arr(4,2)=-0.09402415739519023846534044d0
   fac_arr(5,2)=-0.2820724721855707153960213d0
   fac_arr(6,2)=-0.07051811804639267884900533d0
   fac_arr(7,2)=-0.01175301967439877980816756d0/rhol**2d0
   fac_arr(8,2)=0.01175301967439877980816756d0/rhol**2d0
   fac_arr(9,2)=0.05876509837199389904083778d0/rhol**2d0
   fac_arr(10,2)=0.03525905902319633942450267d0/rhol**2d0
   fac_arr(11,2)=-0.02350603934879755961633511d0/rhol**2d0
   fac_arr(12,2)=0.04701207869759511923267022d0/rhol**2d0
   fac_arr(13,2)=0.07051811804639267884900533d0/rhol**2d0
   fac_arr(14,2)=-0.01175301967439877980816756d0/rhol**2d0
   fac_arr(15,2)=0.03525905902319633942450267d0/rhol**2d0
   lxyz_arr(1,1,3)=5 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=3 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=1 ; lxyz_arr(2,3,3)=4 ; lxyz_arr(3,3,3)=1
   lxyz_arr(1,4,3)=3 ; lxyz_arr(2,4,3)=0 ; lxyz_arr(3,4,3)=3
   lxyz_arr(1,5,3)=1 ; lxyz_arr(2,5,3)=2 ; lxyz_arr(3,5,3)=3
   lxyz_arr(1,6,3)=7 ; lxyz_arr(2,6,3)=0 ; lxyz_arr(3,6,3)=1
   lxyz_arr(1,7,3)=5 ; lxyz_arr(2,7,3)=2 ; lxyz_arr(3,7,3)=1
   lxyz_arr(1,8,3)=3 ; lxyz_arr(2,8,3)=4 ; lxyz_arr(3,8,3)=1
   lxyz_arr(1,9,3)=1 ; lxyz_arr(2,9,3)=6 ; lxyz_arr(3,9,3)=1
   lxyz_arr(1,10,3)=5 ; lxyz_arr(2,10,3)=0 ; lxyz_arr(3,10,3)=3
   lxyz_arr(1,11,3)=3 ; lxyz_arr(2,11,3)=2 ; lxyz_arr(3,11,3)=3
   lxyz_arr(1,12,3)=1 ; lxyz_arr(2,12,3)=4 ; lxyz_arr(3,12,3)=3
   lxyz_arr(1,13,3)=3 ; lxyz_arr(2,13,3)=0 ; lxyz_arr(3,13,3)=5
   lxyz_arr(1,14,3)=1 ; lxyz_arr(2,14,3)=2 ; lxyz_arr(3,14,3)=5
   fac_arr(1,3)=0.04701207869759511923267022d0
   fac_arr(2,3)=-0.09402415739519023846534044d0
   fac_arr(3,3)=-0.1410362360927853576980107d0
   fac_arr(4,3)=0.04701207869759511923267022d0
   fac_arr(5,3)=-0.1410362360927853576980107d0
   fac_arr(6,3)=-0.01175301967439877980816756d0/rhol**2d0
   fac_arr(7,3)=0.01175301967439877980816756d0/rhol**2d0
   fac_arr(8,3)=0.05876509837199389904083778d0/rhol**2d0
   fac_arr(9,3)=0.03525905902319633942450267d0/rhol**2d0
   fac_arr(10,3)=-0.02350603934879755961633511d0/rhol**2d0
   fac_arr(11,3)=0.04701207869759511923267022d0/rhol**2d0
   fac_arr(12,3)=0.07051811804639267884900533d0/rhol**2d0
   fac_arr(13,3)=-0.01175301967439877980816756d0/rhol**2d0
   fac_arr(14,3)=0.03525905902319633942450267d0/rhol**2d0
else if (l.eq.4 .and. i.eq.3 .and. m.eq.5) then
   nterm_arr(1)=15
   nterm_arr(2)=18
   nterm_arr(3)=14
   lxyz_arr(1,1,1)=5 ; lxyz_arr(2,1,1)=1 ; lxyz_arr(3,1,1)=0
   lxyz_arr(1,2,1)=3 ; lxyz_arr(2,2,1)=3 ; lxyz_arr(3,2,1)=0
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=5 ; lxyz_arr(3,3,1)=0
   lxyz_arr(1,4,1)=3 ; lxyz_arr(2,4,1)=1 ; lxyz_arr(3,4,1)=2
   lxyz_arr(1,5,1)=1 ; lxyz_arr(2,5,1)=3 ; lxyz_arr(3,5,1)=2
   lxyz_arr(1,6,1)=1 ; lxyz_arr(2,6,1)=1 ; lxyz_arr(3,6,1)=4
   lxyz_arr(1,7,1)=7 ; lxyz_arr(2,7,1)=1 ; lxyz_arr(3,7,1)=0
   lxyz_arr(1,8,1)=5 ; lxyz_arr(2,8,1)=3 ; lxyz_arr(3,8,1)=0
   lxyz_arr(1,9,1)=3 ; lxyz_arr(2,9,1)=5 ; lxyz_arr(3,9,1)=0
   lxyz_arr(1,10,1)=1 ; lxyz_arr(2,10,1)=7 ; lxyz_arr(3,10,1)=0
   lxyz_arr(1,11,1)=5 ; lxyz_arr(2,11,1)=1 ; lxyz_arr(3,11,1)=2
   lxyz_arr(1,12,1)=3 ; lxyz_arr(2,12,1)=3 ; lxyz_arr(3,12,1)=2
   lxyz_arr(1,13,1)=1 ; lxyz_arr(2,13,1)=5 ; lxyz_arr(3,13,1)=2
   lxyz_arr(1,14,1)=3 ; lxyz_arr(2,14,1)=1 ; lxyz_arr(3,14,1)=4
   lxyz_arr(1,15,1)=1 ; lxyz_arr(2,15,1)=3 ; lxyz_arr(3,15,1)=4
   fac_arr(1,1)=-0.211554354139178036547016d0
   fac_arr(2,1)=-0.2350603934879755961633511d0
   fac_arr(3,1)=-0.02350603934879755961633511d0
   fac_arr(4,1)=-0.2820724721855707153960213d0
   fac_arr(5,1)=-0.09402415739519023846534044d0
   fac_arr(6,1)=-0.07051811804639267884900533d0
   fac_arr(7,1)=0.03525905902319633942450267d0/rhol**2d0
   fac_arr(8,1)=0.05876509837199389904083778d0/rhol**2d0
   fac_arr(9,1)=0.01175301967439877980816756d0/rhol**2d0
   fac_arr(10,1)=-0.01175301967439877980816756d0/rhol**2d0
   fac_arr(11,1)=0.07051811804639267884900533d0/rhol**2d0
   fac_arr(12,1)=0.04701207869759511923267022d0/rhol**2d0
   fac_arr(13,1)=-0.02350603934879755961633511d0/rhol**2d0
   fac_arr(14,1)=0.03525905902319633942450267d0/rhol**2d0
   fac_arr(15,1)=-0.01175301967439877980816756d0/rhol**2d0
   lxyz_arr(1,1,2)=6 ; lxyz_arr(2,1,2)=0 ; lxyz_arr(3,1,2)=0
   lxyz_arr(1,2,2)=4 ; lxyz_arr(2,2,2)=2 ; lxyz_arr(3,2,2)=0
   lxyz_arr(1,3,2)=2 ; lxyz_arr(2,3,2)=4 ; lxyz_arr(3,3,2)=0
   lxyz_arr(1,4,2)=0 ; lxyz_arr(2,4,2)=6 ; lxyz_arr(3,4,2)=0
   lxyz_arr(1,5,2)=4 ; lxyz_arr(2,5,2)=0 ; lxyz_arr(3,5,2)=2
   lxyz_arr(1,6,2)=2 ; lxyz_arr(2,6,2)=2 ; lxyz_arr(3,6,2)=2
   lxyz_arr(1,7,2)=0 ; lxyz_arr(2,7,2)=4 ; lxyz_arr(3,7,2)=2
   lxyz_arr(1,8,2)=2 ; lxyz_arr(2,8,2)=0 ; lxyz_arr(3,8,2)=4
   lxyz_arr(1,9,2)=0 ; lxyz_arr(2,9,2)=2 ; lxyz_arr(3,9,2)=4
   lxyz_arr(1,10,2)=6 ; lxyz_arr(2,10,2)=2 ; lxyz_arr(3,10,2)=0
   lxyz_arr(1,11,2)=4 ; lxyz_arr(2,11,2)=4 ; lxyz_arr(3,11,2)=0
   lxyz_arr(1,12,2)=2 ; lxyz_arr(2,12,2)=6 ; lxyz_arr(3,12,2)=0
   lxyz_arr(1,13,2)=0 ; lxyz_arr(2,13,2)=8 ; lxyz_arr(3,13,2)=0
   lxyz_arr(1,14,2)=4 ; lxyz_arr(2,14,2)=2 ; lxyz_arr(3,14,2)=2
   lxyz_arr(1,15,2)=2 ; lxyz_arr(2,15,2)=4 ; lxyz_arr(3,15,2)=2
   lxyz_arr(1,16,2)=0 ; lxyz_arr(2,16,2)=6 ; lxyz_arr(3,16,2)=2
   lxyz_arr(1,17,2)=2 ; lxyz_arr(2,17,2)=2 ; lxyz_arr(3,17,2)=4
   lxyz_arr(1,18,2)=0 ; lxyz_arr(2,18,2)=4 ; lxyz_arr(3,18,2)=4
   fac_arr(1,2)=-0.03525905902319633942450267d0
   fac_arr(2,2)=-0.1762952951159816971225133d0
   fac_arr(3,2)=-0.05876509837199389904083778d0
   fac_arr(4,2)=0.08227113772079145865717289d0
   fac_arr(5,2)=-0.07051811804639267884900533d0
   fac_arr(6,2)=-0.1410362360927853576980107d0
   fac_arr(7,2)=0.1175301967439877980816756d0
   fac_arr(8,2)=-0.03525905902319633942450267d0
   fac_arr(9,2)=0.03525905902319633942450267d0
   fac_arr(10,2)=0.03525905902319633942450267d0/rhol**2d0
   fac_arr(11,2)=0.05876509837199389904083778d0/rhol**2d0
   fac_arr(12,2)=0.01175301967439877980816756d0/rhol**2d0
   fac_arr(13,2)=-0.01175301967439877980816756d0/rhol**2d0
   fac_arr(14,2)=0.07051811804639267884900533d0/rhol**2d0
   fac_arr(15,2)=0.04701207869759511923267022d0/rhol**2d0
   fac_arr(16,2)=-0.02350603934879755961633511d0/rhol**2d0
   fac_arr(17,2)=0.03525905902319633942450267d0/rhol**2d0
   fac_arr(18,2)=-0.01175301967439877980816756d0/rhol**2d0
   lxyz_arr(1,1,3)=4 ; lxyz_arr(2,1,3)=1 ; lxyz_arr(3,1,3)=1
   lxyz_arr(1,2,3)=2 ; lxyz_arr(2,2,3)=3 ; lxyz_arr(3,2,3)=1
   lxyz_arr(1,3,3)=0 ; lxyz_arr(2,3,3)=5 ; lxyz_arr(3,3,3)=1
   lxyz_arr(1,4,3)=2 ; lxyz_arr(2,4,3)=1 ; lxyz_arr(3,4,3)=3
   lxyz_arr(1,5,3)=0 ; lxyz_arr(2,5,3)=3 ; lxyz_arr(3,5,3)=3
   lxyz_arr(1,6,3)=6 ; lxyz_arr(2,6,3)=1 ; lxyz_arr(3,6,3)=1
   lxyz_arr(1,7,3)=4 ; lxyz_arr(2,7,3)=3 ; lxyz_arr(3,7,3)=1
   lxyz_arr(1,8,3)=2 ; lxyz_arr(2,8,3)=5 ; lxyz_arr(3,8,3)=1
   lxyz_arr(1,9,3)=0 ; lxyz_arr(2,9,3)=7 ; lxyz_arr(3,9,3)=1
   lxyz_arr(1,10,3)=4 ; lxyz_arr(2,10,3)=1 ; lxyz_arr(3,10,3)=3
   lxyz_arr(1,11,3)=2 ; lxyz_arr(2,11,3)=3 ; lxyz_arr(3,11,3)=3
   lxyz_arr(1,12,3)=0 ; lxyz_arr(2,12,3)=5 ; lxyz_arr(3,12,3)=3
   lxyz_arr(1,13,3)=2 ; lxyz_arr(2,13,3)=1 ; lxyz_arr(3,13,3)=5
   lxyz_arr(1,14,3)=0 ; lxyz_arr(2,14,3)=3 ; lxyz_arr(3,14,3)=5
   fac_arr(1,3)=-0.1410362360927853576980107d0
   fac_arr(2,3)=-0.09402415739519023846534044d0
   fac_arr(3,3)=0.04701207869759511923267022d0
   fac_arr(4,3)=-0.1410362360927853576980107d0
   fac_arr(5,3)=0.04701207869759511923267022d0
   fac_arr(6,3)=0.03525905902319633942450267d0/rhol**2d0
   fac_arr(7,3)=0.05876509837199389904083778d0/rhol**2d0
   fac_arr(8,3)=0.01175301967439877980816756d0/rhol**2d0
   fac_arr(9,3)=-0.01175301967439877980816756d0/rhol**2d0
   fac_arr(10,3)=0.07051811804639267884900533d0/rhol**2d0
   fac_arr(11,3)=0.04701207869759511923267022d0/rhol**2d0
   fac_arr(12,3)=-0.02350603934879755961633511d0/rhol**2d0
   fac_arr(13,3)=0.03525905902319633942450267d0/rhol**2d0
   fac_arr(14,3)=-0.01175301967439877980816756d0/rhol**2d0
else if (l.eq.4 .and. i.eq.3 .and. m.eq.6) then
   nterm_arr(1)=13
   nterm_arr(2)=13
   nterm_arr(3)=16
   lxyz_arr(1,1,1)=5 ; lxyz_arr(2,1,1)=0 ; lxyz_arr(3,1,1)=1
   lxyz_arr(1,2,1)=3 ; lxyz_arr(2,2,1)=2 ; lxyz_arr(3,2,1)=1
   lxyz_arr(1,3,1)=1 ; lxyz_arr(2,3,1)=4 ; lxyz_arr(3,3,1)=1
   lxyz_arr(1,4,1)=3 ; lxyz_arr(2,4,1)=0 ; lxyz_arr(3,4,1)=3
   lxyz_arr(1,5,1)=1 ; lxyz_arr(2,5,1)=0 ; lxyz_arr(3,5,1)=5
   lxyz_arr(1,6,1)=7 ; lxyz_arr(2,6,1)=0 ; lxyz_arr(3,6,1)=1
   lxyz_arr(1,7,1)=5 ; lxyz_arr(2,7,1)=2 ; lxyz_arr(3,7,1)=1
   lxyz_arr(1,8,1)=3 ; lxyz_arr(2,8,1)=4 ; lxyz_arr(3,8,1)=1
   lxyz_arr(1,9,1)=1 ; lxyz_arr(2,9,1)=6 ; lxyz_arr(3,9,1)=1
   lxyz_arr(1,10,1)=5 ; lxyz_arr(2,10,1)=0 ; lxyz_arr(3,10,1)=3
   lxyz_arr(1,11,1)=1 ; lxyz_arr(2,11,1)=4 ; lxyz_arr(3,11,1)=3
   lxyz_arr(1,12,1)=3 ; lxyz_arr(2,12,1)=0 ; lxyz_arr(3,12,1)=5
   lxyz_arr(1,13,1)=1 ; lxyz_arr(2,13,1)=2 ; lxyz_arr(3,13,1)=5
   fac_arr(1,1)=0.1727334068350121925245643d0
   fac_arr(2,1)=0.1151556045566747950163762d0
   fac_arr(3,1)=-0.05757780227833739750818811d0
   fac_arr(4,1)=0.2303112091133495900327524d0
   fac_arr(5,1)=0.05757780227833739750818811d0
   fac_arr(6,1)=-0.02878890113916869875409405d0/rhol**2d0
   fac_arr(7,1)=-0.02878890113916869875409405d0/rhol**2d0
   fac_arr(8,1)=0.02878890113916869875409405d0/rhol**2d0
   fac_arr(9,1)=0.02878890113916869875409405d0/rhol**2d0
   fac_arr(10,1)=-0.05757780227833739750818811d0/rhol**2d0
   fac_arr(11,1)=0.05757780227833739750818811d0/rhol**2d0
   fac_arr(12,1)=-0.02878890113916869875409405d0/rhol**2d0
   fac_arr(13,1)=0.02878890113916869875409405d0/rhol**2d0
   lxyz_arr(1,1,2)=4 ; lxyz_arr(2,1,2)=1 ; lxyz_arr(3,1,2)=1
   lxyz_arr(1,2,2)=2 ; lxyz_arr(2,2,2)=3 ; lxyz_arr(3,2,2)=1
   lxyz_arr(1,3,2)=0 ; lxyz_arr(2,3,2)=5 ; lxyz_arr(3,3,2)=1
   lxyz_arr(1,4,2)=0 ; lxyz_arr(2,4,2)=3 ; lxyz_arr(3,4,2)=3
   lxyz_arr(1,5,2)=0 ; lxyz_arr(2,5,2)=1 ; lxyz_arr(3,5,2)=5
   lxyz_arr(1,6,2)=6 ; lxyz_arr(2,6,2)=1 ; lxyz_arr(3,6,2)=1
   lxyz_arr(1,7,2)=4 ; lxyz_arr(2,7,2)=3 ; lxyz_arr(3,7,2)=1
   lxyz_arr(1,8,2)=2 ; lxyz_arr(2,8,2)=5 ; lxyz_arr(3,8,2)=1
   lxyz_arr(1,9,2)=0 ; lxyz_arr(2,9,2)=7 ; lxyz_arr(3,9,2)=1
   lxyz_arr(1,10,2)=4 ; lxyz_arr(2,10,2)=1 ; lxyz_arr(3,10,2)=3
   lxyz_arr(1,11,2)=0 ; lxyz_arr(2,11,2)=5 ; lxyz_arr(3,11,2)=3
   lxyz_arr(1,12,2)=2 ; lxyz_arr(2,12,2)=1 ; lxyz_arr(3,12,2)=5
   lxyz_arr(1,13,2)=0 ; lxyz_arr(2,13,2)=3 ; lxyz_arr(3,13,2)=5
   fac_arr(1,2)=0.05757780227833739750818811d0
   fac_arr(2,2)=-0.1151556045566747950163762d0
   fac_arr(3,2)=-0.1727334068350121925245643d0
   fac_arr(4,2)=-0.2303112091133495900327524d0
   fac_arr(5,2)=-0.05757780227833739750818811d0
   fac_arr(6,2)=-0.02878890113916869875409405d0/rhol**2d0
   fac_arr(7,2)=-0.02878890113916869875409405d0/rhol**2d0
   fac_arr(8,2)=0.02878890113916869875409405d0/rhol**2d0
   fac_arr(9,2)=0.02878890113916869875409405d0/rhol**2d0
   fac_arr(10,2)=-0.05757780227833739750818811d0/rhol**2d0
   fac_arr(11,2)=0.05757780227833739750818811d0/rhol**2d0
   fac_arr(12,2)=-0.02878890113916869875409405d0/rhol**2d0
   fac_arr(13,2)=0.02878890113916869875409405d0/rhol**2d0
   lxyz_arr(1,1,3)=6 ; lxyz_arr(2,1,3)=0 ; lxyz_arr(3,1,3)=0
   lxyz_arr(1,2,3)=4 ; lxyz_arr(2,2,3)=2 ; lxyz_arr(3,2,3)=0
   lxyz_arr(1,3,3)=2 ; lxyz_arr(2,3,3)=4 ; lxyz_arr(3,3,3)=0
   lxyz_arr(1,4,3)=0 ; lxyz_arr(2,4,3)=6 ; lxyz_arr(3,4,3)=0
   lxyz_arr(1,5,3)=4 ; lxyz_arr(2,5,3)=0 ; lxyz_arr(3,5,3)=2
   lxyz_arr(1,6,3)=0 ; lxyz_arr(2,6,3)=4 ; lxyz_arr(3,6,3)=2
   lxyz_arr(1,7,3)=2 ; lxyz_arr(2,7,3)=0 ; lxyz_arr(3,7,3)=4
   lxyz_arr(1,8,3)=0 ; lxyz_arr(2,8,3)=2 ; lxyz_arr(3,8,3)=4
   lxyz_arr(1,9,3)=6 ; lxyz_arr(2,9,3)=0 ; lxyz_arr(3,9,3)=2
   lxyz_arr(1,10,3)=4 ; lxyz_arr(2,10,3)=2 ; lxyz_arr(3,10,3)=2
   lxyz_arr(1,11,3)=2 ; lxyz_arr(2,11,3)=4 ; lxyz_arr(3,11,3)=2
   lxyz_arr(1,12,3)=0 ; lxyz_arr(2,12,3)=6 ; lxyz_arr(3,12,3)=2
   lxyz_arr(1,13,3)=4 ; lxyz_arr(2,13,3)=0 ; lxyz_arr(3,13,3)=4
   lxyz_arr(1,14,3)=0 ; lxyz_arr(2,14,3)=4 ; lxyz_arr(3,14,3)=4
   lxyz_arr(1,15,3)=2 ; lxyz_arr(2,15,3)=0 ; lxyz_arr(3,15,3)=6
   lxyz_arr(1,16,3)=0 ; lxyz_arr(2,16,3)=2 ; lxyz_arr(3,16,3)=6
   fac_arr(1,3)=0.02878890113916869875409405d0
   fac_arr(2,3)=0.02878890113916869875409405d0
   fac_arr(3,3)=-0.02878890113916869875409405d0
   fac_arr(4,3)=-0.02878890113916869875409405d0
   fac_arr(5,3)=0.1727334068350121925245643d0
   fac_arr(6,3)=-0.1727334068350121925245643d0
   fac_arr(7,3)=0.1439445056958434937704703d0
   fac_arr(8,3)=-0.1439445056958434937704703d0
   fac_arr(9,3)=-0.02878890113916869875409405d0/rhol**2d0
   fac_arr(10,3)=-0.02878890113916869875409405d0/rhol**2d0
   fac_arr(11,3)=0.02878890113916869875409405d0/rhol**2d0
   fac_arr(12,3)=0.02878890113916869875409405d0/rhol**2d0
   fac_arr(13,3)=-0.05757780227833739750818811d0/rhol**2d0
   fac_arr(14,3)=0.05757780227833739750818811d0/rhol**2d0
   fac_arr(15,3)=-0.02878890113916869875409405d0/rhol**2d0
   fac_arr(16,3)=0.02878890113916869875409405d0/rhol**2d0
else if (l.eq.4 .and. i.eq.3 .and. m.eq.7) then
   nterm_arr(1)=12
   nterm_arr(2)=12
   nterm_arr(3)=12
   lxyz_arr(1,1,1)=4 ; lxyz_arr(2,1,1)=1 ; lxyz_arr(3,1,1)=1
   lxyz_arr(1,2,1)=2 ; lxyz_arr(2,2,1)=3 ; lxyz_arr(3,2,1)=1
   lxyz_arr(1,3,1)=0 ; lxyz_arr(2,3,1)=5 ; lxyz_arr(3,3,1)=1
   lxyz_arr(1,4,1)=2 ; lxyz_arr(2,4,1)=1 ; lxyz_arr(3,4,1)=3
   lxyz_arr(1,5,1)=0 ; lxyz_arr(2,5,1)=3 ; lxyz_arr(3,5,1)=3
   lxyz_arr(1,6,1)=0 ; lxyz_arr(2,6,1)=1 ; lxyz_arr(3,6,1)=5
   lxyz_arr(1,7,1)=6 ; lxyz_arr(2,7,1)=1 ; lxyz_arr(3,7,1)=1
   lxyz_arr(1,8,1)=4 ; lxyz_arr(2,8,1)=3 ; lxyz_arr(3,8,1)=1
   lxyz_arr(1,9,1)=2 ; lxyz_arr(2,9,1)=5 ; lxyz_arr(3,9,1)=1
   lxyz_arr(1,10,1)=4 ; lxyz_arr(2,10,1)=1 ; lxyz_arr(3,10,1)=3
   lxyz_arr(1,11,1)=2 ; lxyz_arr(2,11,1)=3 ; lxyz_arr(3,11,1)=3
   lxyz_arr(1,12,1)=2 ; lxyz_arr(2,12,1)=1 ; lxyz_arr(3,12,1)=5
   fac_arr(1,1)=0.2878890113916869875409405d0
   fac_arr(2,1)=0.3454668136700243850491286d0
   fac_arr(3,1)=0.05757780227833739750818811d0
   fac_arr(4,1)=0.3454668136700243850491286d0
   fac_arr(5,1)=0.1151556045566747950163762d0
   fac_arr(6,1)=0.05757780227833739750818811d0
   fac_arr(7,1)=-0.05757780227833739750818811d0/rhol**2d0
   fac_arr(8,1)=-0.1151556045566747950163762d0/rhol**2d0
   fac_arr(9,1)=-0.05757780227833739750818811d0/rhol**2d0
   fac_arr(10,1)=-0.1151556045566747950163762d0/rhol**2d0
   fac_arr(11,1)=-0.1151556045566747950163762d0/rhol**2d0
   fac_arr(12,1)=-0.05757780227833739750818811d0/rhol**2d0
   lxyz_arr(1,1,2)=5 ; lxyz_arr(2,1,2)=0 ; lxyz_arr(3,1,2)=1
   lxyz_arr(1,2,2)=3 ; lxyz_arr(2,2,2)=2 ; lxyz_arr(3,2,2)=1
   lxyz_arr(1,3,2)=1 ; lxyz_arr(2,3,2)=4 ; lxyz_arr(3,3,2)=1
   lxyz_arr(1,4,2)=3 ; lxyz_arr(2,4,2)=0 ; lxyz_arr(3,4,2)=3
   lxyz_arr(1,5,2)=1 ; lxyz_arr(2,5,2)=2 ; lxyz_arr(3,5,2)=3
   lxyz_arr(1,6,2)=1 ; lxyz_arr(2,6,2)=0 ; lxyz_arr(3,6,2)=5
   lxyz_arr(1,7,2)=5 ; lxyz_arr(2,7,2)=2 ; lxyz_arr(3,7,2)=1
   lxyz_arr(1,8,2)=3 ; lxyz_arr(2,8,2)=4 ; lxyz_arr(3,8,2)=1
   lxyz_arr(1,9,2)=1 ; lxyz_arr(2,9,2)=6 ; lxyz_arr(3,9,2)=1
   lxyz_arr(1,10,2)=3 ; lxyz_arr(2,10,2)=2 ; lxyz_arr(3,10,2)=3
   lxyz_arr(1,11,2)=1 ; lxyz_arr(2,11,2)=4 ; lxyz_arr(3,11,2)=3
   lxyz_arr(1,12,2)=1 ; lxyz_arr(2,12,2)=2 ; lxyz_arr(3,12,2)=5
   fac_arr(1,2)=0.05757780227833739750818811d0
   fac_arr(2,2)=0.3454668136700243850491286d0
   fac_arr(3,2)=0.2878890113916869875409405d0
   fac_arr(4,2)=0.1151556045566747950163762d0
   fac_arr(5,2)=0.3454668136700243850491286d0
   fac_arr(6,2)=0.05757780227833739750818811d0
   fac_arr(7,2)=-0.05757780227833739750818811d0/rhol**2d0
   fac_arr(8,2)=-0.1151556045566747950163762d0/rhol**2d0
   fac_arr(9,2)=-0.05757780227833739750818811d0/rhol**2d0
   fac_arr(10,2)=-0.1151556045566747950163762d0/rhol**2d0
   fac_arr(11,2)=-0.1151556045566747950163762d0/rhol**2d0
   fac_arr(12,2)=-0.05757780227833739750818811d0/rhol**2d0
   lxyz_arr(1,1,3)=5 ; lxyz_arr(2,1,3)=1 ; lxyz_arr(3,1,3)=0
   lxyz_arr(1,2,3)=3 ; lxyz_arr(2,2,3)=3 ; lxyz_arr(3,2,3)=0
   lxyz_arr(1,3,3)=1 ; lxyz_arr(2,3,3)=5 ; lxyz_arr(3,3,3)=0
   lxyz_arr(1,4,3)=3 ; lxyz_arr(2,4,3)=1 ; lxyz_arr(3,4,3)=2
   lxyz_arr(1,5,3)=1 ; lxyz_arr(2,5,3)=3 ; lxyz_arr(3,5,3)=2
   lxyz_arr(1,6,3)=1 ; lxyz_arr(2,6,3)=1 ; lxyz_arr(3,6,3)=4
   lxyz_arr(1,7,3)=5 ; lxyz_arr(2,7,3)=1 ; lxyz_arr(3,7,3)=2
   lxyz_arr(1,8,3)=3 ; lxyz_arr(2,8,3)=3 ; lxyz_arr(3,8,3)=2
   lxyz_arr(1,9,3)=1 ; lxyz_arr(2,9,3)=5 ; lxyz_arr(3,9,3)=2
   lxyz_arr(1,10,3)=3 ; lxyz_arr(2,10,3)=1 ; lxyz_arr(3,10,3)=4
   lxyz_arr(1,11,3)=1 ; lxyz_arr(2,11,3)=3 ; lxyz_arr(3,11,3)=4
   lxyz_arr(1,12,3)=1 ; lxyz_arr(2,12,3)=1 ; lxyz_arr(3,12,3)=6
   fac_arr(1,3)=0.05757780227833739750818811d0
   fac_arr(2,3)=0.1151556045566747950163762d0
   fac_arr(3,3)=0.05757780227833739750818811d0
   fac_arr(4,3)=0.3454668136700243850491286d0
   fac_arr(5,3)=0.3454668136700243850491286d0
   fac_arr(6,3)=0.2878890113916869875409405d0
   fac_arr(7,3)=-0.05757780227833739750818811d0/rhol**2d0
   fac_arr(8,3)=-0.1151556045566747950163762d0/rhol**2d0
   fac_arr(9,3)=-0.05757780227833739750818811d0/rhol**2d0
   fac_arr(10,3)=-0.1151556045566747950163762d0/rhol**2d0
   fac_arr(11,3)=-0.1151556045566747950163762d0/rhol**2d0
   fac_arr(12,3)=-0.05757780227833739750818811d0/rhol**2d0
else
   stop 'PSP format error'
end if
END SUBROUTINE calc_coeff_derproj


!> Eliminate the translational forces before calling this subroutine!!!
!! Main subroutine: Input is nat (number of atoms), rat0 (atomic positions) and fat (forces on atoms)
!! The atomic positions will be returned untouched
!! In fat, the rotational forces will be eliminated with respect to the center of mass. 
!! All atoms are treated equally (same atomic mass) 
subroutine elim_torque_reza(nat,rat0,fat)
  use module_base
  implicit none
  integer, intent(in) :: nat
  real(gp), dimension(3*nat), intent(in) :: rat0
  real(gp), dimension(3*nat), intent(inout) :: fat
  !local variables
  character(len=*), parameter :: subname='elim_torque_reza'
  integer :: i,iat,i_all,i_stat
  real(gp) :: vrotnrm,cmx,cmy,cmz,alpha,totmass
  !this is an automatic array but it should be allocatable
  real(gp), dimension(3) :: evaleria
  real(gp), dimension(3,3) :: teneria
  real(gp), dimension(3*nat) :: rat
  real(gp), dimension(3*nat,3) :: vrot
  real(gp), dimension(:), allocatable :: amass
  
  allocate(amass(nat+ndebug),stat=i_stat)
  call memocc(i_stat,amass,'amass',subname)

  rat=rat0
  amass(1:nat)=1.0_gp
  !project out rotations
  totmass=0.0_gp
  cmx=0.0_gp 
  cmy=0.0_gp
  cmz=0.0_gp
  do i=1,3*nat-2,3
     iat=(i+2)/3
     cmx=cmx+amass(iat)*rat(i+0)
     cmy=cmy+amass(iat)*rat(i+1)
     cmz=cmz+amass(iat)*rat(i+2)
     totmass=totmass+amass(iat)
  enddo
  cmx=cmx/totmass 
  cmy=cmy/totmass 
  cmz=cmz/totmass
  do i=1,3*nat-2,3
     rat(i+0)=rat(i+0)-cmx
     rat(i+1)=rat(i+1)-cmy
     rat(i+2)=rat(i+2)-cmz
  enddo

  call moment_of_inertia(nat,rat,teneria,evaleria)
  do iat=1,nat
     i=iat*3-2
     call cross(teneria(1,1),rat(i),vrot(i,1))
     call cross(teneria(1,2),rat(i),vrot(i,2))
     call cross(teneria(1,3),rat(i),vrot(i,3))
  enddo
  call normalizevector(3*nat,vrot(1,1))
  call normalizevector(3*nat,vrot(1,2))
  call normalizevector(3*nat,vrot(1,3))
  
  do i=1,3*nat-2,3
     rat(i+0)=rat(i+0)+cmx
     rat(i+1)=rat(i+1)+cmy
     rat(i+2)=rat(i+2)+cmz
  enddo

  vrotnrm=nrm2(3*nat,vrot(1,1),1)
  if (vrotnrm /= 0.0_gp) vrot(1:3*nat,1)=vrot(1:3*nat,1)/vrotnrm
  vrotnrm=nrm2(3*nat,vrot(1,2),1)
  if (vrotnrm /= 0.0_gp) vrot(1:3*nat,2)=vrot(1:3*nat,2)/vrotnrm
  vrotnrm=nrm2(3*nat,vrot(1,3),1)
  if (vrotnrm /= 0.0_gp) vrot(1:3*nat,3)=vrot(1:3*nat,3)/vrotnrm
  
  do i=1,3
     alpha=0.0_gp  
     if(abs(evaleria(i)).gt.1.e-10_gp) then
        alpha=dot_product(vrot(:,i),fat(:))
        fat(:)=fat(:)-alpha*vrot(:,i) 
     endif
  enddo

  i_all=-product(shape(amass))*kind(amass)
  deallocate(amass,stat=i_stat)
  call memocc(i_stat,i_all,'amass',subname)

END SUBROUTINE elim_torque_reza


subroutine cross(a,b,c)
  use module_base
  implicit none
  real(gp), dimension(3), intent(in) :: a,b
  real(gp), dimension(3), intent(out) :: c

  c(1)=a(2)*b(3)-b(2)*a(3)
  c(2)=a(3)*b(1)-b(3)*a(1)
  c(3)=a(1)*b(2)-b(1)*a(2)
END SUBROUTINE cross


subroutine moment_of_inertia(nat,rat,teneria,evaleria)
  use module_base
  implicit none
  integer, intent(in) :: nat
  real(gp), dimension(3,nat), intent(in) :: rat
  real(gp), dimension(3), intent(out) :: evaleria
  real(gp), dimension(3,3), intent(out) :: teneria
  !local variables
  character(len=*), parameter :: subname='moment_of_inertia'
  integer, parameter::lwork=100
  integer :: iat,info,i_all,i_stat
  real(gp) :: tt
  real(gp), dimension(lwork) :: work
  real(gp), dimension(:), allocatable :: amass

  allocate(amass(nat+ndebug),stat=i_stat)
  call memocc(i_stat,amass,'amass',subname)
  
  !positions relative to center of geometry
  amass(1:nat)=1.0_gp
  !calculate inertia tensor
  teneria(1:3,1:3)=0.0_gp
  do iat=1,nat
     tt=amass(iat)
     teneria(1,1)=teneria(1,1)+tt*(rat(2,iat)*rat(2,iat)+rat(3,iat)*rat(3,iat))
     teneria(2,2)=teneria(2,2)+tt*(rat(1,iat)*rat(1,iat)+rat(3,iat)*rat(3,iat))
     teneria(3,3)=teneria(3,3)+tt*(rat(1,iat)*rat(1,iat)+rat(2,iat)*rat(2,iat))
     teneria(1,2)=teneria(1,2)-tt*(rat(1,iat)*rat(2,iat))
     teneria(1,3)=teneria(1,3)-tt*(rat(1,iat)*rat(3,iat))
     teneria(2,3)=teneria(2,3)-tt*(rat(2,iat)*rat(3,iat))
     teneria(2,1)=teneria(1,2)
     teneria(3,1)=teneria(1,3)
     teneria(3,2)=teneria(2,3)
  enddo
  !diagonalize inertia tensor
  call DSYEV('V','L',3,teneria,3,evaleria,work,lwork,info)
  i_all=-product(shape(amass))*kind(amass)
  deallocate(amass,stat=i_stat)
  call memocc(i_stat,i_all,'amass',subname)
  
END SUBROUTINE moment_of_inertia


subroutine normalizevector(n,v)
  use module_base
  implicit none
  integer, intent(in) :: n
  real(gp), dimension(n), intent(inout) :: v
  !local variables
  integer :: i
  real(gp) :: vnrm

  vnrm=0.0_gp
  do i=1,n
     vnrm=vnrm+v(i)**2
  enddo
  vnrm=sqrt(vnrm)
  if (vnrm /= 0.0_gp) v(1:n)=v(1:n)/vnrm

END SUBROUTINE normalizevector


subroutine clean_forces(iproc,at,rxyz,fxyz,fnoise)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(gp), dimension(3,at%nat), intent(inout) :: fxyz
  real(gp), intent(out) :: fnoise
  !local variables
  logical :: move_this_coordinate
  integer :: iat,ixyz
  real(gp) :: sumx,sumy,sumz
  !my variables
  real(gp):: fmax1,t1,t2,t3,fnrm1
  real(gp):: fmax2,fnrm2

  !The maximum force and force norm is computed prior to modification of the forces
  fmax1=0._gp
  fnrm1=0._gp
  do iat=1,at%nat
     t1=fxyz(1,iat)**2
     t2=fxyz(2,iat)**2
     t3=fxyz(3,iat)**2
     fmax1=max(fmax1,sqrt(t1+t2+t3))
     fnrm1=fnrm1+t1+t2+t3
  enddo
  
  
  sumx=0.0_gp
  sumy=0.0_gp
  sumz=0.0_gp
  do iat=1,at%nat
     sumx=sumx+fxyz(1,iat)
     sumy=sumy+fxyz(2,iat)
     sumz=sumz+fxyz(3,iat)
  enddo
  if (at%nat /= 0) then 
     fnoise=sqrt((sumx**2+sumy**2+sumz**2)/real(at%nat,gp))
     sumx=sumx/real(at%nat,gp)
     sumy=sumy/real(at%nat,gp)
     sumz=sumz/real(at%nat,gp)
  else
     fnoise = 0.0_gp
  end if

  if (iproc==0) then 
     !write( *,'(1x,a,1x,3(1x,1pe9.2))') &
     !  'Subtracting center-mass shift of',sumx,sumy,sumz
!           write(*,'(1x,a)')'the sum of the forces is'
           write(*,'(a,1pe16.8)')' average noise along x direction: ',sumx*sqrt(real(at%nat,gp))
           write(*,'(a,1pe16.8)')' average noise along y direction: ',sumy*sqrt(real(at%nat,gp))
           write(*,'(a,1pe16.8)')' average noise along z direction: ',sumz*sqrt(real(at%nat,gp))
           write(*,'(a,1pe16.8)')' total average noise            : ',sqrt(sumx**2+sumy**2+sumz**2)*sqrt(real(at%nat,gp))
!!$
!!$     write(*,'(a,1x,1pe24.17)') 'translational force along x=', sumx  
!!$     write(*,'(a,1x,1pe24.17)') 'translational force along y=', sumy  
!!$     write(*,'(a,1x,1pe24.17)') 'translational force along z=', sumz  
  end if
  
  if (at%geocode == 'F') then
     do iat=1,at%nat
        fxyz(1,iat)=fxyz(1,iat)-sumx
        fxyz(2,iat)=fxyz(2,iat)-sumy
        fxyz(3,iat)=fxyz(3,iat)-sumz
     enddo
     
     call elim_torque_reza(at%nat,rxyz,fxyz)
     
  else if (at%geocode == 'S') then
     do iat=1,at%nat
        fxyz(2,iat)=fxyz(2,iat)-sumy
     enddo
  end if
  
  !clean the forces for blocked atoms
  do iat=1,at%nat
     do ixyz=1,3
        if (.not. move_this_coordinate(at%ifrztyp(iat),ixyz)) fxyz(ixyz,iat)=0.0_gp
     end do
  end do
  
  !the noise of the forces is the norm of the translational force
!  fnoise=real(at%nat,gp)**2*(sumx**2+sumy**2+sumz**2)

  !The maximum force and force norm is computed after modification of the forces
  fmax2=0._gp
  fnrm2=0._gp
  do iat=1,at%nat
     t1=fxyz(1,iat)**2
     t2=fxyz(2,iat)**2
     t3=fxyz(3,iat)**2
     fmax2=max(fmax2,sqrt(t1+t2+t3))
     fnrm2=fnrm2+t1+t2+t3
  enddo

  if (iproc==0) then
     write(*,'(2(1x,a,1pe20.12))') 'clean forces norm (Ha/Bohr): maxval=', fmax2, ' fnrm2=', fnrm2
     if (at%geocode /= 'P') &
  &  write(*,'(2(1x,a,1pe20.12))') 'raw forces:                  maxval=', fmax1, ' fnrm2=', fnrm1
  end if
END SUBROUTINE clean_forces


!> Symmetrize stress (important with special k points)
!@todo: modifiy the arguments of this routine
subroutine symm_stress(dump,tens,symobj)
  use defs_basis
  use module_base, only: verbose,gp
  use m_ab6_symmetry
  use module_types
  use yaml_output
  implicit none
  !Arguments
  logical, intent(in) :: dump
  integer, intent(in) :: symobj
  real(gp), dimension(6), intent(inout) :: tens
  !Local variables
  integer, pointer  :: sym(:,:,:)
  integer, pointer  :: symAfm(:)
  real(gp), pointer :: transNon(:,:)
  integer :: isym, errno, nsym,k,l
  integer, allocatable :: symrec(:,:,:)
  real(gp),dimension(3,3) :: symtens

  call symmetry_get_matrices_p(symObj, nsym, sym, transNon, symAfm, errno)
  if (errno /= AB6_NO_ERROR) stop
  if (nsym < 2) return

  if (dump) call yaml_map('Number of Symmetries',nsym,fmt='(i0)')
  !write(*,"(1x,A,I0,A)") "Symmetrize stress tensor with ", nsym, "symmetries."

  !Get the symmetry matrices in terms of reciprocal basis
  allocate(symrec(3, 3, nsym))
  do isym = 1, nsym, 1
     call mati3inv(sym(:,:,isym), symrec(:,:,isym))
  end do

  symtens=0.0_gp
  do isym = 1,nsym
     do k=1,3
        do l=1,3
           symtens(k,l)=&
                symtens(k,l)+&
                sym(1,k,isym)*tens(1)*sym(1,l,isym)+&
                sym(1,k,isym)*tens(6)*sym(2,l,isym)+&
                sym(1,k,isym)*tens(5)*sym(3,l,isym)+&
                sym(2,k,isym)*tens(6)*sym(1,l,isym)+&
                sym(2,k,isym)*tens(2)*sym(2,l,isym)+&
                sym(2,k,isym)*tens(4)*sym(3,l,isym)+&
                sym(3,k,isym)*tens(5)*sym(1,l,isym)+&
                sym(3,k,isym)*tens(4)*sym(2,l,isym)+&
                sym(3,k,isym)*tens(3)*sym(3,l,isym)
        end do
     end do
  end do
  symtens=symtens / real(nsym,gp)

  tens(1)=symtens(1,1)
  tens(2)=symtens(2,2)
  tens(3)=symtens(3,3)
  tens(4)=symtens(2,3)
  tens(5)=symtens(1,3)
  tens(6)=symtens(1,2)

!  if (iproc == 0 .and. verbose > 2) then
!     write(*,*) '=== SYMMETRISED ==='
!     write(*,*) tens(:)
!  end if

end subroutine symm_stress

!> Symmetrise the atomic forces (needed with special k points)
subroutine symmetrise_forces(iproc, fxyz, at)
  use defs_basis
  use m_ab6_symmetry
  use module_types

  implicit none

  integer, intent(in) :: iproc
  type(atoms_data), intent(in) :: at
  real(gp), intent(inout) :: fxyz(3, at%nat)
  integer :: ia, mu, isym, errno, ind, nsym
  integer :: indsym(4, AB6_MAX_SYMMETRIES)
  real(gp) :: summ
  real(gp) :: alat(3)
  real(gp), allocatable :: dedt(:,:)
  integer, allocatable :: symrec(:,:,:)
  integer, pointer  :: sym(:,:,:)
  integer, pointer  :: symAfm(:)
  real(gp), pointer :: transNon(:,:)

  call symmetry_get_matrices_p(at%sym%symObj, nsym, sym, transNon, symAfm, errno)
  if (errno /= AB6_NO_ERROR) stop
  if (nsym < 2) return

  if (iproc == 0) write(*,"(1x,A,I0,A)") "Symmetrise forces with ", nsym, " symmetries."

  !Get the symmetry matrices in terms of reciprocal basis
  allocate(symrec(3, 3, nsym))
  do isym = 1, nsym, 1
     call mati3inv(sym(:,:,isym), symrec(:,:,isym))
  end do

  alat = (/ at%alat1, at%alat2, at%alat3 /)
  if (at%geocode == 'S') alat(2) = real(1, gp)

  !Save fxyz into dedt.
  allocate(dedt(3,at%nat))
  do ia = 1, at%nat
     dedt(:, ia) = fxyz(:, ia) / alat
  end do

  ! actually conduct symmetrization
  do ia = 1, at%nat
     call symmetry_get_equivalent_atom(at%sym%symObj, indsym, ia, errno)
     if (errno /= AB6_NO_ERROR) stop
     do mu = 1, 3
        summ = real(0, gp)
        do isym = 1, nsym
           ind = indsym(4, isym)
           summ = summ + real(symrec(mu,1,isym), gp) * dedt(1, ind) + &
                & real(symrec(mu,2,isym), gp) * dedt(2, ind) + &
                & real(symrec(mu,3,isym), gp) * dedt(3, ind)
        end do
        fxyz(mu, ia) = summ / real(nsym, gp)
        ! if (abs(fred(mu, ia))<tol)fred(mu,ia)=0.0_dp
     end do
  end do

  deallocate(dedt)
  deallocate(symrec)
  
  ! fxyz is in reduced coordinates, we expand here.
  do ia = 1, at%nat
     fxyz(:, ia) = fxyz(:, ia) * alat
  end do
end subroutine symmetrise_forces


subroutine local_hamiltonian_stress(orbs,lr,hx,hy,hz,psi,tens)
  use module_base
  use module_types
  use module_interfaces
  use libxc_functionals
  implicit none
  real(gp), intent(in) :: hx,hy,hz
  type(orbitals_data), intent(in) :: orbs
  type(locreg_descriptors), intent(in) :: lr
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*orbs%norbp), intent(in) :: psi
  real(gp), intent(inout) :: tens(6)
   real(gp) :: ekin_sum,epot_sum
  !local variables
  character(len=*), parameter :: subname='local_hamiltonian_stress'
  integer :: i_all,i_stat,iorb,npot,oidx
  real(wp) :: exctXcoeff,kinstr(6)
  real(gp) :: ekin,kx,ky,kz,etest
  type(workarr_locham) :: wrk_lh
  real(wp), dimension(:,:), allocatable :: psir,hpsi

  exctXcoeff=libxc_functionals_exctXfac()

  !initialise the work arrays
  call initialize_work_arrays_locham(lr,orbs%nspinor,wrk_lh)  

  tens=0.d0

  !components of the potential
  npot=orbs%nspinor
  if (orbs%nspinor == 2) npot=1

  allocate(hpsi(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*orbs%norbp),stat=i_stat)
  call memocc(i_stat,hpsi,'hpsi',subname)
  hpsi=0.0_wp
  ! Wavefunction in real space
  allocate(psir(lr%d%n1i*lr%d%n2i*lr%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
  call memocc(i_stat,psir,'psir',subname)
  call razero(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspinor,psir)



  ekin_sum=0.0_gp
  epot_sum=0.0_gp

  etest=0.0_gp


  do iorb=1,orbs%norbp
     kinstr=0._wp
     oidx=(iorb-1)*orbs%nspinor+1

     call daub_to_isf_locham(orbs%nspinor,lr,wrk_lh,psi(1,oidx),psir)

     kx=orbs%kpts(1,orbs%iokpt(iorb))
     ky=orbs%kpts(2,orbs%iokpt(iorb))
     kz=orbs%kpts(3,orbs%iokpt(iorb))

     call isf_to_daub_kinetic(hx,hy,hz,kx,ky,kz,orbs%nspinor,lr,wrk_lh,&
          psir,hpsi(1,oidx),ekin,kinstr)

     kinstr = -kinstr*8.0_gp/(hx*hy*hz)/real(lr%d%n1i*lr%d%n2i*lr%d%n3i,gp)
     tens=tens+kinstr*2.0_gp*&
     orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)

  end do !loop over orbitals: finished

  !deallocations of work arrays
  i_all=-product(shape(psir))*kind(psir)
  deallocate(psir,stat=i_stat)
  call memocc(i_stat,i_all,'psir',subname)

  i_all=-product(shape(hpsi))*kind(hpsi)
  deallocate(hpsi,stat=i_stat)
  call memocc(i_stat,i_all,'hpsi',subname)

  call deallocate_work_arrays_locham(lr,wrk_lh)

END SUBROUTINE local_hamiltonian_stress


subroutine erf_stress(at,rxyz,hxh,hyh,hzh,n1i,n2i,n3i,n3p,iproc,nproc,ngatherarr,rho,tens)
  use module_base
  use module_types
  use module_fft_sg
  implicit none
  !passed var
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%nat), target, intent(in) :: rxyz
  real(gp), intent(in) :: hxh,hyh,hzh
  integer,intent(in) :: n1i,n2i,n3i,n3p,iproc,nproc
  real(kind=8), dimension(n1i*n2i*max(n3p,1)), intent(in), target :: rho
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
  real(dp),dimension(6), intent(out) :: tens
  !local var
  character(len=*), parameter :: subname='erf_stress'
  real(kind=8),allocatable :: rhog(:,:,:,:,:)
  real(kind=8),dimension(:),pointer :: rhor
  integer :: ierr,i_stat,i_all
  real(kind=8) :: pi,p(3),g2,rloc,set,fac
  real(kind=8) :: rx,ry,rz,sfr,sfi,rhore,rhoim
  real(kind=8) :: potg,potg2
  real(kind=8) :: Zion
  integer :: iat,ityp
  integer :: j1,j2,j3,i1,i2,i3,inzee,ind

  !write(*,*) 'iproc,n3i,n3p',iproc,n3i,n3p
  !write(*,*) 'iproc',iproc, ngatherarr(iproc-1,1),ngatherarr(iproc-1,2)

  if (nproc > 1) then
     allocate(rhor(n1i*n2i*n3i),stat=i_stat)
     call memocc(i_stat,rhor,'rhor',subname)
     call MPI_ALLGATHERV(rho(1),ngatherarr(iproc,1),&
          &   mpidtypw,rhor(1),ngatherarr(0,1),&
          ngatherarr(0,2),mpidtypw,bigdft_mpi%mpi_comm,ierr)
  else
     rhor => rho        
  end if

  pi = 4.0_gp*atan(1.0_gp)
  allocate(rhog(2,n1i+1,n2i+1,n3i+1,2))
  tens(:)=0.0_dp ; p(:)=0.0_dp

  ! calculate total rho(G)
  rhog=0.0_dp
  do i3=1,n3i
     do i2=1,n2i
        do i1=1,n1i
           ind=i1+(i2-1)*n1i+(i3-1)*n1i*n2i
           rhog(1,i1,i2,i3,1)=rhor(ind)
           rhog(2,i1,i2,i3,1)=0.d0
        end do
     end do
  end do

  ! DO FFT OF DENSITY: Rho(r) -FFT-> Rho(G)
  inzee=1
  call FFT(n1i,n2i,n3i,n1i+1,n2i+1,n3i+1,rhog,-1,inzee)   
  rhog=rhog/real(n1i*n2i*n3i,kind=8)

  tens=0.0_dp

  do iat=1,at%nat                          ! SUM OVER ATOMS

     ityp=at%iatype(iat)                ! ityp
     rloc=at%psppar(0,0,ityp)           ! take corresp. r_loc
     Zion = real(at%nelpsp(ityp),kind=8)

     !coordinates of new center (atom)
     rx=real(rxyz(1,iat),kind=8) 
     ry=real(rxyz(2,iat),kind=8) 
     rz=real(rxyz(3,iat),kind=8)

     potg=0.0_dp

     do i3=1,n3i
        j3=i3-(i3/(n3i/2+2))*n3i-1
        p(3)=real(j3,dp)/(n3i*hzh)         

        do i2=1,n2i
           j2=i2-(i2/(n2i/2+2))*n2i-1
           p(2)=real(j2,dp)/(n2i*hyh)    

           do i1=1,n1i
              j1=i1-(i1/(n1i/2+2))*n1i-1
              p(1)=real(j1,dp)/(n1i*hxh)

              ! calculate structural factor exp(-2*pi*i*R*G)
              sfr=cos(-2.0_gp*pi*(rx*p(1)+ry*p(2)+rz*p(3)))
              sfi=sin(-2.0_gp*pi*(rx*p(1)+ry*p(2)+rz*p(3)))

              ! multiply density by shift /structural/ factor 
              rhore=rhog(1,i1,i2,i3,inzee)*sfr
              rhoim=rhog(2,i1,i2,i3,inzee)*sfi

              ! use analytical expression for rhog^el
              ! multiply by  Kernel 1/pi/G^2 in reciprocal space

              g2=real(p(1)**2+p(2)**2+p(3)**2,kind=8)

              !set = rhog^el (analytic)
              fac=(Zion/rloc**3.0_gp)/sqrt(2.0_gp*pi)/(2.0_gp*pi)
              fac=fac/real(n1i*hxh*n2i*hyh*n3i*hzh,kind=8)                    !Division by Volume
              set=((sqrt(pi*2.0_gp*rloc**2.0_gp))**3)*fac*exp(-pi*pi*g2*2.0_gp*rloc**2.0_gp)

              if (g2 /= 0) then

                 potg = -set/(pi*g2)  ! V^el(G)
                 potg2 = (set/pi)*((real(1.d0,kind=8)/g2**2.d0)&
                      +real(pi*pi*2.d0*rloc**2.d0/g2,kind=8))

                 !STRESS TENSOR
                 tens(1)=tens(1)-(rhore+rhoim)*(potg2*2.0_gp*p(1)*p(1)+potg)
                 tens(2)=tens(2)-(rhore+rhoim)*(potg2*2.0_gp*p(2)*p(2)+potg)
                 tens(3)=tens(3)-(rhore+rhoim)*(potg2*2.0_gp*p(3)*p(3)+potg)
                 tens(6)=tens(6)-(rhore+rhoim)*(potg2*2.0_gp*p(1)*p(2))
                 tens(5)=tens(5)-(rhore+rhoim)*(potg2*2.0_gp*p(1)*p(3))
                 tens(4)=tens(4)-(rhore+rhoim)*(potg2*2.0_gp*p(2)*p(3))

              end if  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! g2 /=0
           end do !i1
        end do !i2
     end do !i3

  end do !iat -atoms

!!$  if (iproc ==0) then
!!$     write(*,*)
!!$     write(*,*)'--------------------------------------------------------------------'
!!$     write(*,*) 'STRESS TENSOR: ERF PART'
!!$     write(*,*) tens(1),tens(2),tens(3)
!!$     write(*,*) tens(4),tens(5),tens(6)
!!$     write(*,*)
!!$  end if

  if (nproc>1) then
     i_all=-product(shape(rhor))*kind(rhor)
     deallocate(rhor,stat=i_stat)
     call memocc(i_stat,i_all,'rhor',subname)
  else
     nullify(rhor)
  end if
  deallocate(rhog)

END SUBROUTINE erf_stress







!> Calculates the nonlocal forces on all atoms arising from the wavefunctions 
!! belonging to iproc and adds them to the force array
!! recalculate the projectors at the end if refill flag is .true.
subroutine nonlocal_forces_linear(iproc,nproc,lr,hx,hy,hz,at,rxyz,&
     orbs,nlpspd,proj,lzd,phi,kernel,fsep,refill,strten)
  use module_base
  use module_types
  implicit none
  !Arguments-------------
  type(atoms_data), intent(in) :: at
  type(local_zone_descriptors), intent(in) :: lzd
  type(nonlocal_psp_descriptors), intent(in) :: nlpspd
  logical, intent(in) :: refill
  integer, intent(in) :: iproc, nproc
  real(gp), intent(in) :: hx,hy,hz
  type(locreg_descriptors) :: lr
  type(orbitals_data), intent(in) :: orbs
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(wp), dimension(orbs%npsidim_orbs), intent(in) :: phi
  real(gp), dimension(orbs%norb,orbs%norb),intent(in) :: kernel
  real(wp), dimension(nlpspd%nprojel), intent(inout) :: proj
  real(gp), dimension(3,at%nat), intent(inout) :: fsep
  real(gp), dimension(6), intent(out) :: strten
  !local variables--------------
  character(len=*), parameter :: subname='nonlocal_forces'
  integer :: istart_c,iproj,iat,ityp,i,j,l,m,iorbout,iiorb,ilr
  integer :: mbseg_c,mbseg_f,jseg_c,jseg_f
  integer :: mbvctr_c,mbvctr_f,iorb,nwarnings,nspinor,ispinor,jorbd
  real(gp) :: offdiagcoeff,hij,sp0,spi,sp0i,sp0j,spj,orbfac,strc,Enl,vol
  integer :: idir,i_all,i_stat,ncplx,icplx,isorb,ikpt,ieorb,istart_ck,ispsi_k,ispsi,jorb,jproc,ii,ist,ierr,iiat
  real(gp), dimension(2,2,3) :: offdiagarr
  real(gp), dimension(:,:), allocatable :: fxyz_orb
  real(dp), dimension(:,:,:,:,:,:,:), allocatable :: scalprod
  real(gp), dimension(6) :: sab
  integer,dimension(:),allocatable :: nat_par, isat_par, sendcounts, recvcounts, senddspls, recvdspls
  real(dp),dimension(:,:,:,:,:,:,:),allocatable :: scalprod_sendbuf
  real(dp),dimension(:),allocatable :: scalprod_recvbuf


  ! Determine how many atoms each MPI task will handle
  allocate(nat_par(0:nproc-1),stat=i_stat)
  call memocc(i_stat,nat_par,'nat_par',subname)
  allocate(isat_par(0:nproc-1),stat=i_stat)
  call memocc(i_stat,isat_par,'isat_par',subname)
  ii=at%nat/nproc
  nat_par(0:nproc-1)=ii
  ii=at%nat-ii*nproc
  do i=0,ii-1
      nat_par(i)=nat_par(i)+1
  end do
  isat_par(0)=0
  do jproc=1,nproc-1
      isat_par(jproc)=isat_par(jproc-1)+nat_par(jproc-1)
  end do

  allocate(sendcounts(0:nproc-1),stat=i_stat)
  call memocc(i_stat,sendcounts,'sendcounts',subname)
  allocate(recvcounts(0:nproc-1),stat=i_stat)
  call memocc(i_stat,recvcounts,'recvcounts',subname)
  allocate(senddspls(0:nproc-1),stat=i_stat)
  call memocc(i_stat,senddspls,'senddspls',subname)
  allocate(recvdspls(0:nproc-1),stat=i_stat)
  call memocc(i_stat,recvdspls,'recvdspls',subname)

  do jproc=0,nproc-1
      sendcounts(jproc)=2*10*7*3*4*orbs%norbp*nat_par(jproc)
      recvcounts(jproc)=2*10*7*3*4*orbs%norb_par(jproc,0)*nat_par(iproc)
  end do
  senddspls(0)=0
  recvdspls(0)=0
  do jproc=1,nproc-1
      senddspls(jproc)=senddspls(jproc-1)+sendcounts(jproc-1)
      recvdspls(jproc)=recvdspls(jproc-1)+recvcounts(jproc-1)
  end do

  call to_zero(6,strten(1)) 

     
  !always put complex scalprod
  !also nspinor for the moment is the biggest as possible

  !  allocate(scalprod(2,0:3,7,3,4,at%nat,orbs%norbp*orbs%nspinor+ndebug),stat=i_stat)
  ! need more components in scalprod to calculate terms like dp/dx*psi*x
  allocate(scalprod(2,0:9,7,3,4,at%nat,orbs%norbp*orbs%nspinor+ndebug),stat=i_stat)
  call memocc(i_stat,scalprod,'scalprod',subname)
  call razero(2*10*7*3*4*at%nat*orbs%norbp*orbs%nspinor,scalprod(1,0,1,1,1,1,1))

  allocate(scalprod_sendbuf(2,0:9,7,3,4,orbs%norbp*orbs%nspinor+ndebug,at%nat),stat=i_stat)
  call memocc(i_stat,scalprod_sendbuf,'scalprod_sendbuf',subname)
  call razero(2*10*7*3*4*at%nat*orbs%norbp*orbs%nspinor,scalprod_sendbuf(1,0,1,1,1,1,1))

  allocate(scalprod_recvbuf(2*10*7*3*4*nat_par(iproc)*orbs%norb*orbs%nspinor+ndebug),stat=i_stat)
  call memocc(i_stat,scalprod_recvbuf,'scalprod_recvbuf',subname)
  call razero(2*10*7*3*4*nat_par(iproc)*orbs%norb*orbs%nspinor,scalprod_recvbuf(1))



  Enl=0._gp
  !strten=0.d0
  vol=real(at%alat1*at%alat2*at%alat3,gp)
  sab=0.d0

  !calculate the coefficients for the off-diagonal terms
  do l=1,3
     do i=1,2
        do j=i+1,3
           offdiagcoeff=0.0_gp
           if (l==1) then
              if (i==1) then
                 if (j==2) offdiagcoeff=-0.5_gp*sqrt(3._gp/5._gp)
                 if (j==3) offdiagcoeff=0.5_gp*sqrt(5._gp/21._gp)
              else
                 offdiagcoeff=-0.5_gp*sqrt(100._gp/63._gp)
              end if
           else if (l==2) then
              if (i==1) then
                 if (j==2) offdiagcoeff=-0.5_gp*sqrt(5._gp/7._gp)
                 if (j==3) offdiagcoeff=1._gp/6._gp*sqrt(35._gp/11._gp)
              else
                 offdiagcoeff=-7._gp/3._gp*sqrt(1._gp/11._gp)
              end if
           else if (l==3) then
              if (i==1) then
                 if (j==2) offdiagcoeff=-0.5_gp*sqrt(7._gp/9._gp)
                 if (j==3) offdiagcoeff=0.5_gp*sqrt(63._gp/143._gp)
              else
                 offdiagcoeff=-9._gp*sqrt(1._gp/143._gp)
              end if
           end if
           offdiagarr(i,j-i,l)=offdiagcoeff
        end do
     end do
  end do

  norbp_if: if (orbs%norbp>0) then

      !look for the strategy of projectors application
      if (DistProjApply) then
         !apply the projectors on the fly for each k-point of the processor
         !starting k-point
         ikpt=orbs%iokpt(1)
         ispsi_k=1
         jorb=0
         loop_kptD: do
    
            call orbs_in_kpt(ikpt,orbs,isorb,ieorb,nspinor)
    
            call ncplx_kpt(ikpt,orbs,ncplx)
    
            nwarnings=0 !not used, simply initialised 
            iproj=0 !should be equal to four times nproj at the end
            jorbd=jorb
            do iat=1,at%nat
    
               call plr_segs_and_vctrs(nlpspd%plr(iat),&
                    mbseg_c,mbseg_f,mbvctr_c,mbvctr_f)
               jseg_c=1
               jseg_f=1
    
               do idir=0,9
    !!$           mbseg_c=nlpspd%nseg_p(2*iat-1)-nlpspd%nseg_p(2*iat-2)
    !!$           mbseg_f=nlpspd%nseg_p(2*iat  )-nlpspd%nseg_p(2*iat-1)
    !!$           jseg_c=nlpspd%nseg_p(2*iat-2)+1
    !!$           jseg_f=nlpspd%nseg_p(2*iat-1)+1
    !!$           mbvctr_c=nlpspd%nvctr_p(2*iat-1)-nlpspd%nvctr_p(2*iat-2)
    !!$           mbvctr_f=nlpspd%nvctr_p(2*iat  )-nlpspd%nvctr_p(2*iat-1)
    
               ityp=at%iatype(iat)
                  !calculate projectors
                  istart_c=1
                  call atom_projector(ikpt,iat,idir,istart_c,iproj,nlpspd%nprojel,&
                       lr,hx,hy,hz,rxyz(1,iat),at,orbs,nlpspd%plr(iat),&
                       proj,nwarnings)
    !              print *,'iat,ilr,idir,sum(proj)',iat,ilr,idir,sum(proj)
     
                  !calculate the contribution for each orbital
                  !here the nspinor contribution should be adjusted
                  ! loop over all my orbitals
                  ispsi=ispsi_k
                  jorb=jorbd
                  do iorb=isorb,ieorb
                     iiorb=orbs%isorb+iorb
                     ilr=orbs%inwhichlocreg(iiorb)
                     do ispinor=1,nspinor,ncplx
                        jorb=jorb+1
                        istart_c=1
                        do l=1,4
                           do i=1,3
                              if (at%psppar(l,i,ityp) /= 0.0_gp) then
                                 do m=1,2*l-1
                                    !!do i_stat=ispsi,ispsi+lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f-1
                                    !!    write(200+iproc,*) phi(i_stat)
                                    !!end do
                                    !!do i_stat=istart_c,istart_c+mbvctr_c+7*mbvctr_f-1
                                    !!    write(250+iproc,*) proj(i_stat)
                                    !!end do
                                    call wpdot_wrap(ncplx,&
                                         lzd%llr(ilr)%wfd%nvctr_c,lzd%llr(ilr)%wfd%nvctr_f,&
                                         lzd%llr(ilr)%wfd%nseg_c,lzd%llr(ilr)%wfd%nseg_f,&
                                         lzd%llr(ilr)%wfd%keyvglob,lzd%llr(ilr)%wfd%keyglob,phi(ispsi),&
                                         mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
    !!$                                     nlpspd%keyv_p(jseg_c),&
    !!$                                     nlpspd%keyg_p(1,jseg_c),&
                                         nlpspd%plr(iat)%wfd%keyvglob(jseg_c),&
                                         nlpspd%plr(iat)%wfd%keyglob(1,jseg_c),&
                                         proj(istart_c),&
                                         scalprod(1,idir,m,i,l,iat,jorb))
                                    istart_c=istart_c+(mbvctr_c+7*mbvctr_f)*ncplx
                                 end do
                              end if
                           end do
                        end do
                        ispsi=ispsi+(lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f)*ncplx
                     end do
                  end do
                  if (istart_c-1  > nlpspd%nprojel) stop '2:applyprojectors'
               end do
    
            end do
    
            if (ieorb == orbs%norbp) exit loop_kptD
            ikpt=ikpt+1
            ispsi_k=ispsi
         end do loop_kptD
    
      else
         !calculate all the scalar products for each direction and each orbitals
         do idir=0,9
    
            if (idir /= 0) then !for the first run the projectors are already allocated
               call fill_projectors(iproc,lr,hx,hy,hz,at,orbs,rxyz,nlpspd,proj,idir)
            end if
            !apply the projectors  k-point of the processor
            !starting k-point
            ikpt=orbs%iokpt(1)
            istart_ck=1
            ispsi_k=1
            jorb=0
            loop_kpt: do
    
               call orbs_in_kpt(ikpt,orbs,isorb,ieorb,nspinor)
    
               call ncplx_kpt(ikpt,orbs,ncplx)
    
               ! calculate the scalar product for all the orbitals
               ispsi=ispsi_k
               do iorb=isorb,ieorb
                  iiorb=orbs%isorb+iorb
                  ilr=orbs%inwhichlocreg(iiorb)
                  do ispinor=1,nspinor,ncplx
                     jorb=jorb+1
                     ! loop over all projectors of this k-point
                     iproj=0
                     istart_c=istart_ck
                     do iat=1,at%nat
                        call plr_segs_and_vctrs(nlpspd%plr(iat),&
                             mbseg_c,mbseg_f,mbvctr_c,mbvctr_f)
                        jseg_c=1
                        jseg_f=1
    
    !!$                    mbseg_c=nlpspd%nseg_p(2*iat-1)-nlpspd%nseg_p(2*iat-2)
    !!$                    mbseg_f=nlpspd%nseg_p(2*iat  )-nlpspd%nseg_p(2*iat-1)
    !!$                    jseg_c=nlpspd%nseg_p(2*iat-2)+1
    !!$                    jseg_f=nlpspd%nseg_p(2*iat-1)+1
    !!$                    mbvctr_c=nlpspd%nvctr_p(2*iat-1)-nlpspd%nvctr_p(2*iat-2)
    !!$                    mbvctr_f=nlpspd%nvctr_p(2*iat  )-nlpspd%nvctr_p(2*iat-1)
                        ityp=at%iatype(iat)
                        do l=1,4
                           do i=1,3
                              if (at%psppar(l,i,ityp) /= 0.0_gp) then
                                 do m=1,2*l-1
                                    iproj=iproj+1
                                    call wpdot_wrap(ncplx,&
                                         lzd%llr(ilr)%wfd%nvctr_c,lzd%llr(ilr)%wfd%nvctr_f,&
                                         lzd%llr(ilr)%wfd%nseg_c,lzd%llr(ilr)%wfd%nseg_f,&
                                         lzd%llr(ilr)%wfd%keyvglob,lzd%llr(ilr)%wfd%keyglob,phi(ispsi),  &
                                         mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
    !!$                                     nlpspd%keyv_p(jseg_c),&
    !!$                                     nlpspd%keyg_p(1,jseg_c),&
                                         nlpspd%plr(iat)%wfd%keyvglob(jseg_c),&
                                         nlpspd%plr(iat)%wfd%keyglob(1,jseg_c),&
                                         proj(istart_c),scalprod(1,idir,m,i,l,iat,jorb))
                                    istart_c=istart_c+(mbvctr_c+7*mbvctr_f)*ncplx
                                 end do
                              end if
                           end do
                        end do
                     end do
                     ispsi=ispsi+(lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f)*ncplx
                  end do
                  if (iproj /= nlpspd%nproj) stop '1:applyprojectors'
               end do
               istart_ck=istart_c
               if (ieorb == orbs%norbp) exit loop_kpt
               ikpt=ikpt+1
               ispsi_k=ispsi
            end do loop_kpt
            if (istart_ck-1  /= nlpspd%nprojel) stop '2:applyprojectors'
    
         end do
    
         !restore the projectors in the proj array (for on the run forces calc., tails or so)
         if (refill) then 
            call fill_projectors(iproc,lr,hx,hy,hz,at,orbs,rxyz,nlpspd,proj,0)
         end if
    
      end if

  end if norbp_if

  ! Copy scalprod to auxiliary array for communication
  do iorb=1,orbs%norbp
      do iat=1,at%nat
          call dcopy(2*10*7*3*4, scalprod(1,0,1,1,1,iat,iorb), 1, scalprod_sendbuf(1,0,1,1,1,iorb,iat), 1)
          !write(*,'(a,3i7,es18.8)') 'FIRST: iproc, iorb, iat, scalprod(1,0,1,1,1,iat,iorb)', &
          !                           iproc, iorb, iat, scalprod(1,0,1,1,1,iat,iorb)
          !write(*,'(a,3i7,es18.8)') 'TEMP: iproc, iorb, iat, scalprod_sendbuf(1,0,1,1,1,iorb,iat)', &
          !                           iproc, iorb, iat, scalprod_sendbuf(1,0,1,1,1,iorb,iat)
      end do
  end do

  if (nproc>1) then
      call mpi_alltoallv(scalprod_sendbuf, sendcounts, senddspls, mpi_double_precision, &
                         scalprod_recvbuf, recvcounts, recvdspls, mpi_double_precision, &
                         bigdft_mpi%mpi_comm, ierr)
  else
      call dcopy(2*10*7*3*4*at%nat*orbs%norb*orbs%nspinor, scalprod_sendbuf(1,0,1,1,1,1,1), 1, scalprod_recvbuf(1), 1)
  end if
  !write(*,'(a,i7,es18.8)') 'iproc, scalprod_recvbuf(1)', iproc, scalprod_recvbuf(1)
  
  !allocate(scalprod(2,0:9,7,3,4,at%nat,orbs%norbp*orbs%nspinor+ndebug),stat=i_stat)
  !allocate(scalprod_sendbuf(2,0:9,7,3,4,orbs%norbp*orbs%nspinor+ndebug,at%nat),stat=i_stat)
  !allocate(scalprod_recvbuf(2*10*7*3*4*nat_par(iproc)*orbs%norb*orbs%nspinor+ndebug),stat=i_stat)

  i_all=-product(shape(scalprod))*kind(scalprod)
  deallocate(scalprod,stat=i_stat)
  call memocc(i_stat,i_all,'scalprod',subname)
  allocate(scalprod(2,0:9,7,3,4,nat_par(iproc),orbs%norb*orbs%nspinor+ndebug),stat=i_stat)
  call memocc(i_stat,scalprod,'scalprod',subname)
  call to_zero(2*10*7*3*4*nat_par(iproc)*orbs%norb*orbs%nspinor,scalprod(1,0,1,1,1,1,1))

  ist=1
  do jproc=0,nproc-1
      do iat=1,nat_par(iproc)
          iiorb=orbs%isorb_par(jproc)
          !!write(*,'(a,4i8)') 'iproc, jproc, orbs%isorb_par(jproc), iiorb', iproc, jproc, orbs%isorb_par(jproc), iiorb
          do iorb=1,orbs%norb_par(jproc,0)
              iiorb=iiorb+1
              call dcopy(2*10*7*3*4, scalprod_recvbuf(ist), 1, scalprod(1,0,1,1,1,iat,iiorb), 1)
              !write(*,'(a,5i7,2es18.8)') 'SECOND: iproc, jproc, iiorb, iat, ist, scalprod(1,0,1,1,1,iat,iiorb), &
              !                           &scalprod_recvbuf(ist)', &
              !                           iproc, jproc, iiorb, iat, ist, scalprod(1,0,1,1,1,iat,iiorb), scalprod_recvbuf(ist)
              ist=ist+2*10*7*3*4
          end do
      end do
  end do



  allocate(fxyz_orb(3,at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,fxyz_orb,'fxyz_orb',subname)

  natp_if: if (nat_par(iproc)>0) then

      !apply the projectors  k-point of the processor
      !starting k-point
      ikpt=orbs%iokpt(1)
      jorb=0
      loop_kptF: do

         call orbs_in_kpt(ikpt,orbs,isorb,ieorb,nspinor)

         call ncplx_kpt(ikpt,orbs,ncplx)

         call to_zero(3*at%nat,fxyz_orb(1,1))

         ! loop over all my orbitals for calculating forces
         !do iorbout=isorb,ieorb
         do iorbout=1,orbs%norb
            jorb=0 !THIS WILL CREATE PROBLEMS FOR K-POINTS!!
            sab=0.0_gp
            ! loop over all projectors
            !do iorb=isorb,ieorb
            do iorb=1,orbs%norb
               do ispinor=1,nspinor,ncplx
                  jorb=jorb+1
                  !do iat=1,at%nat
                  do iat=1,nat_par(iproc)
                     iiat=isat_par(iproc)+iat
                     ityp=at%iatype(iiat)
                     do l=1,4
                        do i=1,3
                           if (at%psppar(l,i,ityp) /= 0.0_gp) then
                              do m=1,2*l-1
                                 do icplx=1,ncplx
                                    ! scalar product with the derivatives in all the directions
                                    sp0=real(scalprod(icplx,0,m,i,l,iat,iorbout),gp)
                                    !if (kernel(jorb,iorbout)/=0.d0) then
                                        !!write(100+iproc,'(a,9i6,es18.8)') 'iorbout,jorb,icplx,0,m,i,l,iat,iiat,sp0', &
                                        !!                                   iorbout,jorb,icplx,0,m,i,l,iat,iiat,sp0
                                    !end if
                                    do idir=1,3
                                       spi=real(scalprod(icplx,idir,m,i,l,iat,jorb),gp)
                                       fxyz_orb(idir,iiat)=fxyz_orb(idir,iiat)+&
                                            kernel(jorb,iorbout)*at%psppar(l,i,ityp)*sp0*spi
                                       !if (kernel(jorb,iorbout)/=0.d0) then
                                           !!write(110+iproc,'(a,10i6,es18.8)') 'iorbout,jorb,icplx,0,m,i,l,iat,iiat,&
                                           !!                                    &idir,fxyz_orb(idir,iat)', &
                                           !!                                    iorbout,jorb,icplx,0,m,i,l,iat,iiat,&
                                           !!                                    idir,fxyz_orb(idir,iat)
                                       !end if
                                    end do
                                    spi=real(scalprod(icplx,0,m,i,l,iat,jorb),gp)
                                    !!Enl=Enl+sp0*spi*at%psppar(l,i,ityp)*&
                                    !!orbs%occup(iorb+orbs%isorb)*orbs%kwgts(orbs%iokpt(iorb))
                                    !!do idir=4,9 !for stress
                                    !!    strc=real(scalprod(icplx,idir,m,i,l,iat,jorb),gp)
                                    !!    sab(idir-3)=&
                                    !!    sab(idir-3)+&   
                                    !!    at%psppar(l,i,ityp)*sp0*2.0_gp*strc*&
                                    !!    orbs%occup(iorb+orbs%isorb)*orbs%kwgts(orbs%iokpt(iorb))
                                    !!end do
                                 end do
                              end do
                           end if
                        end do
                     end do
                     !HGH case, offdiagonal terms
                     if (at%npspcode(ityp) == 3 .or. at%npspcode(ityp) == 10) then
                        do l=1,3 !no offdiagoanl terms for l=4 in HGH-K case
                           do i=1,2
                              if (at%psppar(l,i,ityp) /= 0.0_gp) then 
                                 loop_j: do j=i+1,3
                                    if (at%psppar(l,j,ityp) == 0.0_gp) exit loop_j
                                    !offdiagonal HGH term
                                    if (at%npspcode(ityp) == 3) then !traditional HGH convention
                                       hij=offdiagarr(i,j-i,l)*at%psppar(l,j,ityp)
                                    else !HGH-K convention
                                       hij=at%psppar(l,i+j+1,ityp)
                                    end if
                                    do m=1,2*l-1
                                       !F_t= 2.0*h_ij (<D_tp_i|psi><psi|p_j>+<p_i|psi><psi|D_tp_j>)
                                       !(the two factor is below)
                                       do icplx=1,ncplx
                                          sp0i=real(scalprod(icplx,0,m,i,l,iat,iorbout),gp)
                                          sp0j=real(scalprod(icplx,0,m,j,l,iat,iorbout),gp)
                                          do idir=1,3
                                             spi=real(scalprod(icplx,idir,m,i,l,iat,jorb),gp)
                                             spj=real(scalprod(icplx,idir,m,j,l,iat,jorb),gp)
                                             fxyz_orb(idir,iiat)=fxyz_orb(idir,iiat)+&
                                                  kernel(jorb,iorbout)*hij*(sp0j*spi+spj*sp0i)
                                          end do
                                          sp0i=real(scalprod(icplx,0,m,i,l,iat,jorb),gp)
                                          !!Enl=Enl+2.0_gp*sp0i*sp0j*hij&
                                          !!*orbs%occup(iorb+orbs%isorb)*orbs%kwgts(orbs%iokpt(iorb))
                                          !!do idir=4,9
                                          !!    spi=real(scalprod(icplx,idir,m,i,l,iat,jorb),gp)
                                          !!    spj=real(scalprod(icplx,idir,m,j,l,iat,jorb),gp)
                                          !!    sab(idir-3)=&
                                          !!    sab(idir-3)+&   
                                          !!    2.0_gp*hij*(sp0j*spi+sp0i*spj)&
                                          !!    *orbs%occup(iorb+orbs%isorb)*orbs%kwgts(orbs%iokpt(iorb))
                                          !!end do
                                       end do
                                    end do
                                 end do loop_j
                              end if
                           end do
                        end do
                     end if
                  end do
               end do
       
               !!!orbital-dependent factor for the forces
               !!orbfac=orbs%kwgts(orbs%iokpt(iorbout))*orbs%occup(iorbout+orbs%isorb)*2.0_gp
       
               !seq: strten(1:6) =  11 22 33 23 13 12 
               strten(1)=strten(1)+sab(1)/vol 
               strten(2)=strten(2)+sab(2)/vol 
               strten(3)=strten(3)+sab(3)/vol 
               strten(4)=strten(4)+sab(5)/vol
               strten(5)=strten(5)+sab(6)/vol
               strten(6)=strten(6)+sab(4)/vol
            end do
         end do
         !do iat=1,at%nat
         do iat=1,nat_par(iproc)
            iiat=isat_par(iproc)+iat
            !write(120+iproc,'(a,2i7,2es18.8)') 'iorbout, iat, fsep(1,iiat), fxyz_orb(1,iiat)', &
            !                                    iorbout, iat, fsep(1,iiat), fxyz_orb(1,iiat)
            !!fsep(1,iat)=fsep(1,iat)+orbfac*fxyz_orb(1,iat)
            !!fsep(2,iat)=fsep(2,iat)+orbfac*fxyz_orb(2,iat)
            !!fsep(3,iat)=fsep(3,iat)+orbfac*fxyz_orb(3,iat)
            fsep(1,iiat)=fsep(1,iiat)+2.d0*fxyz_orb(1,iiat)
            fsep(2,iiat)=fsep(2,iiat)+2.d0*fxyz_orb(2,iiat)
            fsep(3,iiat)=fsep(3,iiat)+2.d0*fxyz_orb(3,iiat)
         end do
         if (ieorb == orbs%norbp) exit loop_kptF
         ikpt=ikpt+1
         ispsi_k=ispsi
      end do loop_kptF

  end if natp_if


!!!Adding Enl to the diagonal components of strten after loop over kpts is finished...
!!do i=1,3
!!strten(i)=strten(i)+Enl/vol
!!end do

!!!  do iat=1,at%nat
!!!     write(20+iat,'(1x,i5,1x,3(1x,1pe12.5))') &
!!!          iat,fsep(1,iat),fsep(2,iat),fsep(3,iat)
!!!  end do

  i_all=-product(shape(fxyz_orb))*kind(fxyz_orb)
  deallocate(fxyz_orb,stat=i_stat)
  call memocc(i_stat,i_all,'fxyz_orb',subname)
  i_all=-product(shape(scalprod))*kind(scalprod)
  deallocate(scalprod,stat=i_stat)
  call memocc(i_stat,i_all,'scalprod',subname)
  i_all=-product(shape(nat_par))*kind(nat_par)
  deallocate(nat_par,stat=i_stat)
  call memocc(i_stat,i_all,'nat_par',subname)
  i_all=-product(shape(isat_par))*kind(isat_par)
  deallocate(isat_par,stat=i_stat)
  call memocc(i_stat,i_all,'isat_par',subname)
  i_all=-product(shape(sendcounts))*kind(sendcounts)
  deallocate(sendcounts,stat=i_stat)
  call memocc(i_stat,i_all,'sendcounts',subname)
  i_all=-product(shape(recvcounts))*kind(recvcounts)
  deallocate(recvcounts,stat=i_stat)
  call memocc(i_stat,i_all,'recvcounts',subname)
  i_all=-product(shape(senddspls))*kind(senddspls)
  deallocate(senddspls,stat=i_stat)
  call memocc(i_stat,i_all,'senddspls',subname)
  i_all=-product(shape(recvdspls))*kind(recvdspls)
  deallocate(recvdspls,stat=i_stat)
  call memocc(i_stat,i_all,'recvdspls',subname)

  i_all=-product(shape(scalprod_sendbuf))*kind(scalprod_sendbuf)
  deallocate(scalprod_sendbuf,stat=i_stat)
  call memocc(i_stat,i_all,'scalprod_sendbuf',subname)
  i_all=-product(shape(scalprod_recvbuf))*kind(scalprod_recvbuf)
  deallocate(scalprod_recvbuf,stat=i_stat)
  call memocc(i_stat,i_all,'scalprod_recvbuf',subname)


END SUBROUTINE nonlocal_forces_linear
