subroutine potentialAndEnergySub(iproc, nproc, Glr, orbs, atoms, in, lin, psi, rhopot, &
    nscatterarr, ngatherarr, GPU, irrzon, phnons, pkernel, pot_ion, rhocore, potxc, PSquiet, &
    proj, nlpspd, pkernelseq, rxyz, eion, edisp, eexctX, scpot, n3d, n3p)
  !
  ! Purpose:
  ! ========
  !   Calculates the potential and energy and writes them. It is just copy&paste from above.
  !

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
real(8),dimension(orbs%npsidim):: psi
real(dp), dimension(lin%as%size_rhopot) :: rhopot
integer,dimension(0:nproc-1,4) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
integer,dimension(0:nproc-1,2),intent(in) :: ngatherarr
type(GPU_pointers):: GPU
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
real(gp):: eion, edisp, eexctX
logical:: scpot

! Local variables
real(8):: hxh, hyh, hzh, ehart, eexcu, vexcu, ekin_sum, epot_sum, eproj_sum, energybs, energy
real(wp), dimension(:), pointer :: potential
real(8),dimension(:),allocatable:: hpsi
integer:: istat
character(len=*),parameter:: subname='potentialAndEnergy'


hxh=0.5d0*in%hx
hyh=0.5d0*in%hy
hzh=0.5d0*in%hz

allocate(hpsi(orbs%npsidim), stat=istat)
hpsi=0.d0



  !calculate the self-consistent potential
  if (scpot) then
      !write(*,*) 'SCPOT'
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

  !deallocate potential
  call free_full_potential(nproc,potential,subname)

  !call HamiltonianApplication(iproc,nproc,atoms,orbs,hx,hy,hz,rxyz,&
  !     nlpspd,proj,Glr,ngatherarr,n1i*n2i*n3p,&
  !     rhopot,psi,hpsi,ekin_sum,epot_sum,eexctX,eproj_sum,&
  !     in%nspin,GPU,pkernel=pkernelseq)

  energybs=ekin_sum+epot_sum+eproj_sum !the potential energy contains also exctX
  energy=energybs-ehart+eexcu-vexcu-eexctX+eion+edisp

  if(iproc==0) write( *,'(1x,a,3(1x,1pe18.11))') 'ekin_sum,epot_sum,eproj_sum',  &
                  ekin_sum,epot_sum,eproj_sum
  if(iproc==0) write( *,'(1x,a,3(1x,1pe18.11))') '   ehart,   eexcu,    vexcu',ehart,eexcu,vexcu
  if(iproc==0) write(*,*) 'energy', energy

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
end subroutine potentialAndEnergySub

