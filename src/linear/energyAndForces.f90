subroutine potentialAndEnergySub(iproc, nproc, n3d, n3p, Glr, orbs, atoms, in, lin, phi, psi, rxyz, rxyzParab, &
    rhopot, nscatterarr, ngatherarr, GPU, irrzon, phnons, pkernel, pot_ion, rhocore, potxc, PSquiet, &
    proj, nlpspd, pkernelseq, eion, edisp, eexctX, scpot, coeff, ebsMod, energy)
!
! Purpose:
! ========
!   Calculates the potential and energy and writes them. This is subroutine is copied
!   from cluster.
!
! Calling arguments:
! ==================
!   Input arguments:
!   -----------------
!     iproc       process ID
!     nproc       total number of processes
!     n3d         ??
!     n3p         ??
!     Glr         type describing the localization region
!     orbs        type describing the physical orbitals psi
!     atoms       type containing the parameters for the atoms
!     in          type  containing some very general parameters
!     lin         type containing parameters for the linear version
!     psi         the physical orbitals
!     rxyz        atomic positions
!     rhopot      the charge density
!     nscatterarr ??
!     ngatherarr  ??
!     nlpspd      ??
!     proj        ??
!     GPU         parameters for GPUs?
!     pkernelseq  ??
!     radii_cf    coarse and fine radii around the atoms
!     irrzon      ??
!     phnons      ??
!     pkernel     ??
!     pot_ion     the ionic potential
!     rhocore     ??
!     potxc       ??
!     PSquiet     flag to control the output from the Poisson solver
!     eion        ionic energy
!     edisp       dispersion energy
!     eexctX      ??
!     scpot       flag indicating whether we have a self consistent calculation
!     fion        ionic forces
!     fdisp       dispersion forces
!   Input / Output arguments
!   ------------------------
!     GPU         parameters for GPUs?
!     rhopot      the charge density
!   Output arguments:
!   -----------------
!     energy      the total energy


use module_base
use module_types
use module_interfaces, exceptThisOne => potentialAndEnergySub
use Poisson_Solver
implicit none

! Calling arguments
integer:: iproc, nproc, n3d, n3p
type(locreg_descriptors) :: Glr
type(orbitals_data):: orbs
type(atoms_data):: atoms
type(input_variables):: in
type(linearParameters):: lin
real(8),dimension(lin%orbs%npsidim):: phi
real(8),dimension(orbs%npsidim):: psi
real(dp), dimension(lin%as%size_rhopot) :: rhopot
integer,dimension(0:nproc-1,4) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
integer,dimension(0:nproc-1,2),intent(in) :: ngatherarr
type(GPU_pointers),intent(in out):: GPU
integer, dimension(lin%as%size_irrzon(1),lin%as%size_irrzon(2),lin%as%size_irrzon(3)) :: irrzon
real(dp), dimension(lin%as%size_phnons(1),lin%as%size_phnons(2),lin%as%size_phnons(3)) :: phnons
real(dp), dimension(lin%as%size_pkernel):: pkernel
real(wp), dimension(lin%as%size_pot_ion):: pot_ion
!real(wp), dimension(lin%as%size_rhocore):: rhocore 
real(wp), dimension(:),pointer:: rhocore 
real(wp), dimension(lin%as%size_potxc(1),lin%as%size_potxc(2),lin%as%size_potxc(3),lin%as%size_potxc(4)):: potxc
character(len=3):: PSquiet
type(nonlocal_psp_descriptors),intent(in) :: nlpspd
real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
real(dp),dimension(lin%as%size_pkernelseq),intent(in):: pkernelseq
real(8),dimension(3,atoms%nat),intent(in):: rxyz
real(8),dimension(3,atoms%nat),intent(in):: rxyzParab
real(gp):: eion, edisp, eexctX, energy
real(8),dimension(lin%orbs%norb,orbs%norb):: coeff
real(8):: ebsMod
logical:: scpot

! Local variables
real(8):: hxh, hyh, hzh, ehart, eexcu, vexcu, ekin_sum, epot_sum, eproj_sum, energybs, energyMod
real(8):: energyMod2, ehartMod
real(wp), dimension(:), pointer :: potential
real(8),dimension(:),allocatable:: hpsi, hphi
real(8),dimension(:,:,:),allocatable:: matrixElements
integer:: istat, iall
character(len=*),parameter:: subname='potentialAndEnergy'


hxh=0.5d0*in%hx
hyh=0.5d0*in%hy
hzh=0.5d0*in%hz

allocate(hpsi(orbs%npsidim), stat=istat)
call memocc(istat, hpsi, 'hpsi', subname)
allocate(hphi(lin%orbs%npsidim), stat=istat)
call memocc(istat, hphi, 'hphi', subname)
allocate(matrixElements(lin%orbs%norb,lin%orbs%norb,2), stat=istat)
call memocc(istat, matrixElements, 'matrixElements', subname)


if(iproc==0) write(*,'(x,a)') '-------------------------------------------------- Calculation of energy and forces.'

  !calculate the self-consistent potential
  if (scpot) then
     ! Potential from electronic charge density
     call sumrho(iproc,nproc,orbs,Glr,in%ixc,hxh,hyh,hzh,psi,rhopot,&
          Glr%d%n1i*Glr%d%n2i*n3d,nscatterarr,in%nspin,GPU,atoms%symObj,irrzon,phnons)

     if(orbs%nspinor==4) then
        !this wrapper can be inserted inside the poisson solver 
        call PSolverNC(atoms%geocode,'D',iproc,nproc,Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,n3d,&
             in%ixc,hxh,hyh,hzh,&
             rhopot,pkernel,pot_ion,ehart,eexcu,vexcu,0.d0,.true.,4)
     else
        call XC_potential(atoms%geocode,'D',iproc,nproc,&
             Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,in%ixc,hxh,hyh,hzh,&
             rhopot,eexcu,vexcu,in%nspin,rhocore,potxc)

        call H_potential(atoms%geocode,'D',iproc,nproc,&
             Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,hxh,hyh,hzh,&
             rhopot,pkernel,pot_ion,ehart,0.0_dp,.true.,&
             quiet=PSquiet) !optional argument

        !sum the two potentials in rhopot array
        !fill the other part, for spin, polarised
        if (in%nspin == 2) then
           call dcopy(Glr%d%n1i*Glr%d%n2i*n3p,rhopot(1),1,&
                rhopot(1+Glr%d%n1i*Glr%d%n2i*n3p),1)
        end if
        !spin up and down together with the XC part
        call axpy(Glr%d%n1i*Glr%d%n2i*n3p*in%nspin,1.0_dp,potxc(1,1,1,1),1,&
             rhopot(1),1)

     end if
  end if

  !allocate the potential in the full box
  call full_local_potential(iproc,nproc,Glr%d%n1i*Glr%d%n2i*n3p,Glr%d%n1i*Glr%d%n2i*Glr%d%n3i,in%nspin,&
       orbs%norb,orbs%norbp,ngatherarr,rhopot,potential)


  call HamiltonianApplication(iproc,nproc,atoms,orbs,in%hx,in%hy,in%hz,rxyz,&
       nlpspd,proj,Glr,ngatherarr,potential,psi,hpsi,ekin_sum,epot_sum,eexctX,eproj_sum,&
       in%nspin,GPU,pkernel=pkernelseq)
  if(iproc==0) write(*,'(x,a)') 'done.'

  !deallocate potential
  call free_full_potential(nproc,potential,subname)

  energybs=ekin_sum+epot_sum+eproj_sum !the potential energy contains also exctX
  energy=energybs-ehart+eexcu-vexcu-eexctX+eion+edisp


  ! Now calculate the modified band structure energy.
  call HamiltonianApplicationConfinement(iproc,nproc,atoms,lin%orbs,lin,in%hx,in%hy,in%hz,rxyz,&
       nlpspd,proj,Glr,ngatherarr,Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2),&
       rhopot(1),&
       phi(1),hphi(1),ekin_sum,epot_sum,eexctX,eproj_sum,in%nspin,GPU, rxyzParab, pkernel=pkernelseq)
  call getMatrixElements(iproc, nproc, Glr, lin, phi, hphi, matrixElements)
  if(trim(lin%getCoeff)=='diag') then
      call modifiedBSEnergyModified(in%nspin, orbs, lin, coeff(1,1), matrixElements(1,1,1), ebsMod)
  else if(trim(lin%getCoeff)=='min') then
  
      call optimizeCoefficients(iproc, orbs, lin, matrixElements, coeff)
      call modifiedBSEnergyModified(in%nspin, orbs, lin, coeff(1,1), matrixElements(1,1,1), ebsMod)
  else
      if(iproc==0) write(*,'(a,a,a)') "ERROR: lin%getCoeff can have the values 'diag' or 'min' , &
          & but we found '", lin%getCoeff, "'."
      stop
  end if


  ! IS THIS CORRECT??
  energyMod=ebsMod-ehart+eexcu-vexcu-eexctX+eion+edisp

  if(iproc==0) write(*,'(x,a,2es26.17)') 'energybs, ebsMod', energybs, ebsMod

  if(iproc==0) write( *,'(1x,a,3(1x,1pe18.11))') 'ekin_sum,epot_sum,eproj_sum',  &
                  ekin_sum,epot_sum,eproj_sum
  if(iproc==0) write( *,'(1x,a,3(1x,1pe18.11))') '   ehart,   eexcu,    vexcu',ehart,eexcu,vexcu
  if(iproc==0) write(*,'(x,a,2es26.17)') 'total energy, modified energy', energy, energyMod

  !!if (iproc == 0) then
  !!   if (verbose > 0 .and. in%itrpmax==1) then
  !!      write( *,'(1x,a,3(1x,1pe18.11))') 'ekin_sum,epot_sum,eproj_sum',  &
  !!           ekin_sum,epot_sum,eproj_sum
  !!      write( *,'(1x,a,3(1x,1pe18.11))') '   ehart,   eexcu,    vexcu',ehart,eexcu,vexcu
  !!   end if
  !!   if (.not. scpot) then
  !!      if (gnrm_zero == 0.0_gp) then
  !!         write( *,'(1x,a,i6,2x,1pe24.17,1x,1pe9.2)') 'iter, tr(H),gnrm',iter,trH,gnrm
  !!      else
  !!         write( *,'(1x,a,i6,2x,1pe24.17,2(1x,1pe9.2))') 'iter, tr(H),gnrm,gnrm_zero',&
  !!             iter,trH,gnrm,gnrm_zero
  !!      end if
  !!   else
  !!      if (gnrm_zero == 0.0_gp) then
  !!         write( *,'(1x,a,i6,2x,1pe24.17,1x,1pe9.2)') 'iter,total energy,gnrm',iter,energy,gnrm
  !!      else
  !!         write( *,'(1x,a,i6,2x,1pe24.17,2(1x,1pe9.2))') 'iter,total energy,gnrm,gnrm_zero',&
  !!             iter,energy,gnrm,gnrm_zero
  !!      end if
  !!   end if
  !!endif

  iall=-product(shape(hpsi))*kind(hpsi)
  deallocate(hpsi, stat=istat)
  call memocc(istat, iall, 'hpsi', subname)

  iall=-product(shape(hphi))*kind(hphi)
  deallocate(hphi, stat=istat)
  call memocc(istat, iall, 'hphi', subname)

  iall=-product(shape(matrixElements))*kind(matrixElements)
  deallocate(matrixElements, stat=istat)
  call memocc(istat, iall, 'matrixElements', subname)



end subroutine potentialAndEnergySub




subroutine calculateForcesSub(iproc, nproc, n3p, i3s, i3xcsh, Glr, orbs, atoms, in, lin, nlpspd, proj, &
    ngatherarr, nscatterarr, GPU, irrzon, phnons, pkernel, rxyz, fion, fdisp, psi, fxyz, fnoise)
! Purpose:
! ========
!   Calculates the forces we get with psi. It is copied from cluster, with an additional
!   write statement to print the forces.
!
! Calling arguments:
! ==================
!   Input arguments:
!   -----------------
!     iproc       process ID
!     nproc       total number of processes
!     n3p         ??
!     i3s         ??
!     i3xcsh      ??
!     Glr         type describing the localization region
!     orbs        type describing the physical orbitals psi
!     atoms       type containing the parameters for the atoms
!     in          type  containing some very general parameters
!     lin         type containing parameters for the linear version
!     nlpspd      ??
!     ngatherarr  ??
!     GPU         parameters for GPUs?
!     irrzon      ??
!     phnons      ??
!     pkernel     ??
!     rxyz        atomic positions
!     fion        ionic forces
!     fdisp       dispersion forces
!     psi         the physical orbitals
!   Input / Output arguments
!   ------------------------
!     proj        ??
!     nscatterarr ??
!     GPU         parameters for GPUs?
!   Output arguments:
!   -----------------
!     fxyz        the forces
!     fnoise      noise of the forces
!
use module_base
use module_types
use Poisson_Solver
use module_interfaces, exceptThisOne => calculateForcesSub
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, n3p, i3s, i3xcsh
type(locreg_descriptors),intent(in):: Glr
type(orbitals_data),intent(in):: orbs
type(atoms_data),intent(in):: atoms
type(input_variables),intent(in):: in
type(linearParameters),intent(in):: lin
type(nonlocal_psp_descriptors),intent(in) :: nlpspd
real(wp),dimension(nlpspd%nprojel),intent(inout) :: proj
integer,dimension(0:nproc-1,2),intent(in) :: ngatherarr   !!! NOT NEEDED
integer,dimension(0:nproc-1,4),intent(inout) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
type(GPU_pointers),intent(inout):: GPU
integer,dimension(lin%as%size_irrzon(1),lin%as%size_irrzon(2),lin%as%size_irrzon(3)),intent(in) :: irrzon
real(dp),dimension(lin%as%size_phnons(1),lin%as%size_phnons(2),lin%as%size_phnons(3)),intent(in) :: phnons
real(dp),dimension(lin%as%size_pkernel),intent(in):: pkernel
real(8),dimension(3,atoms%nat),intent(in):: rxyz, fion, fdisp
real(8),dimension(3,atoms%nat),intent(out):: fxyz
real(8),intent(out):: fnoise
real(8),dimension(orbs%npsidim),intent(in):: psi

! Local variables
integer:: jproc, i_stat, i_all, iat, ierr, j
real(8):: hxh, hyh, hzh, ehart_fake
real(kind=8), dimension(:), allocatable :: rho
real(gp), dimension(:,:), allocatable :: gxyz
real(kind=8), dimension(:,:,:,:), allocatable :: pot
character(len=*),parameter:: subname='calculateForcesSub'
logical:: refill_proj

  hxh=0.5d0*in%hx
  hyh=0.5d0*in%hy
  hzh=0.5d0*in%hz

  if (iproc==0) then
     write( *,'(1x,a)')&
          '----------------------------------------------------------------- Forces Calculation'
  end if

  ! Selfconsistent potential is saved in rhopot, 
  ! new arrays rho,pot for calculation of forces ground state electronic density

  ! Potential from electronic charge density

  !manipulate scatter array for avoiding the GGA shift
  do jproc=0,nproc-1
     !n3d=n3p
     nscatterarr(jproc,1)=nscatterarr(jproc,2)
     !i3xcsh=0
     nscatterarr(jproc,4)=0
  end do

  if (n3p>0) then
     allocate(rho(Glr%d%n1i*Glr%d%n2i*n3p*in%nspin+ndebug),stat=i_stat)
     call memocc(i_stat,rho,'rho',subname)
  else
     allocate(rho(1+ndebug),stat=i_stat)
     call memocc(i_stat,rho,'rho',subname)
  end if
  call sumrho(iproc,nproc,orbs,Glr,0,hxh,hyh,hzh,psi,rho,Glr%d%n1i*Glr%d%n2i*n3p,&
          nscatterarr,in%nspin,GPU,atoms%symObj,irrzon,phnons)

  !calculate the total density in the case of nspin==2
  if (in%nspin==2) then
     call axpy(Glr%d%n1i*Glr%d%n2i*n3p,1.0_dp,rho(1+Glr%d%n1i*Glr%d%n2i*n3p),1,rho(1),1)
  end if
  if (n3p>0) then
     allocate(pot(Glr%d%n1i,Glr%d%n2i,n3p,1+ndebug),stat=i_stat)
     call memocc(i_stat,pot,'pot',subname)
  else
     allocate(pot(1,1,1,1+ndebug),stat=i_stat)
     call memocc(i_stat,pot,'pot',subname)
  end if

  !calculate electrostatic potential
  call dcopy(Glr%d%n1i*Glr%d%n2i*n3p,rho,1,pot,1)
  call H_potential(atoms%geocode,'D',iproc,nproc,&
       Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,hxh,hyh,hzh,pot,pkernel,pot,ehart_fake,0.0_dp,.false.)


  allocate(gxyz(3,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,gxyz,'gxyz',subname)

  call timing(iproc,'Forces        ','ON')
  ! calculate local part of the forces gxyz
  call local_forces(iproc,atoms,rxyz,hxh,hyh,hzh,&
       Glr%d%n1,Glr%d%n2,Glr%d%n3,n3p,i3s+i3xcsh,Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,rho,pot,gxyz)

  i_all=-product(shape(rho))*kind(rho)
  deallocate(rho,stat=i_stat)
  call memocc(i_stat,i_all,'rho',subname)
  i_all=-product(shape(pot))*kind(pot)
  deallocate(pot,stat=i_stat)
  call memocc(i_stat,i_all,'pot',subname)

  if (iproc == 0 .and. verbose > 1) write( *,'(1x,a)',advance='no')'Calculate nonlocal forces...'

  !refill projectors for tails, davidson
  !refill_proj=(in%calc_tail .or. DoDavidson) .and. DoLastRunThings
  refill_proj=.false.  !! IS THIS CORRECT??

  call nonlocal_forces(iproc,Glr%d%n1,Glr%d%n2,Glr%d%n3,in%hx,in%hy,in%hz,atoms,rxyz,&
       orbs,nlpspd,proj,Glr%wfd,psi,gxyz,refill_proj)

  if (iproc == 0 .and. verbose > 1) write( *,'(1x,a)')'done.'

  ! Add up all the force contributions
  if (nproc > 1) then
     call MPI_ALLREDUCE(gxyz,fxyz,3*atoms%nat,mpidtypg,MPI_SUM,MPI_COMM_WORLD,ierr)
  else
     do iat=1,atoms%nat
        fxyz(1,iat)=gxyz(1,iat)
        fxyz(2,iat)=gxyz(2,iat)
        fxyz(3,iat)=gxyz(3,iat)
     enddo
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

    !add to the forces the ionic and dispersion contribution 
    do iat=1,atoms%nat
       fxyz(1,iat)=fxyz(1,iat)+fion(1,iat)+fdisp(1,iat)
       fxyz(2,iat)=fxyz(2,iat)+fion(2,iat)+fdisp(2,iat)
       fxyz(3,iat)=fxyz(3,iat)+fion(3,iat)+fdisp(3,iat)
    enddo

    !i_all=-product(shape(fion))*kind(fion)
    !deallocate(fion,stat=i_stat)
    !call memocc(i_stat,i_all,'fion',subname)
    !i_all=-product(shape(fdisp))*kind(fdisp)
    !deallocate(fdisp,stat=i_stat)
    !call memocc(i_stat,i_all,'fdisp',subname)
    i_all=-product(shape(gxyz))*kind(gxyz)
    deallocate(gxyz,stat=i_stat)
    call memocc(i_stat,i_all,'gxyz',subname)


  !!do iat=1,atoms%nat
  !!   if(iproc==0) write(*,'(a,i0,3es14.5)') 'forces for atom ',iat, fxyz(1,iat), fxyz(2,iat), fxyz(3,iat)
  !!end do

  !subtraction of zero of the forces, disabled for the moment
  !the zero of the forces depends on the atomic positions
  !if (in%gaussian_help .and. .false.) then
  call clean_forces(iproc,atoms,rxyz,fxyz,fnoise)
  !end if

  if(iproc==0) then
      write(*,'(x,a)') 'Force values for all atoms in x, y, z direction.'
      do iat=1,atoms%nat
         write(*,'(3x,i0,x,a6,x,3(x,es12.5))') &
              iat,trim(atoms%atomnames(atoms%iatype(iat))),(fxyz(j,iat),j=1,3)
      end do
  end if

  call timing(iproc,'Forces        ','OF')

end subroutine calculateForcesSub

