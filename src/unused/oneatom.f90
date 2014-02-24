!> @file
!!  Program to do one atom calculation
!! @author
!!    Copyright (C) 2010-2011 ESRF, PoliTo
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Compute one atom system
!! @deprecated
program oneatom
  use BigDFT_API
  use Poisson_Solver
  use gaussians, only: gaussian_basis
  implicit none
  character(len=*), parameter :: subname='oneatom'
  logical :: dokernel=.false.
  logical :: endloop,endlooprp
  integer :: n1i,n2i,n3i,iproc,nproc,i_stat,i_all,nelec
  integer :: n3d,n3p,n3pi,i3xcsh,i3s,n1,n2,n3,ndegree_ip
  integer :: idsx_actual,ndiis_sd_sw,idsx_actual_before,iter,istat
  integer :: iatyp
  real(gp) :: hxh,hyh,hzh
  real(gp) :: tt,gnrm !n(c) gnrm_zero,alpha
  real(gp) :: energy,energy_min,scprsum !n(c) energy_old
  type(atoms_data) :: atoms
  type(input_variables) :: in
  type(orbitals_data) :: orbs
  type(locreg_descriptors) :: Glr
  type(local_zone_descriptors) :: Lzd
  type(nonlocal_psp_descriptors) :: nlpspd
  type(communications_arrays) :: comms
  type(denspot_distribution) :: denspotd
  type(GPU_pointers) :: GPU
  type(diis_objects) :: diis
  type(rho_descriptors)  :: rhodsc
  type(energy_terms) :: energs
  type(gaussian_basis),dimension(:),allocatable::proj_G
  type(paw_objects)::paw
  character(len=4) :: itername
  real(gp), dimension(3) :: shift
  integer, dimension(:,:), allocatable :: nscatterarr,ngatherarr
  real(gp), dimension(:,:), allocatable :: radii_cf
  real(wp), dimension(:), pointer :: hpsi,psit,psi,proj,pot
  real(dp), dimension(:), pointer :: pkernel,pot_ion
  real(gp), dimension(:,:), pointer :: rxyz
  type(confpot_data), dimension(:), allocatable :: confdatarr
  character(len=60) :: radical

  !for the moment no need to have parallelism
  iproc=0
  nproc=1

  !Initilization
  idsx_actual = huge(1)
  energy_min = huge(1.0_gp)

  call memocc_set_memory_limit(memorylimit)

  ! Read a possible radical format argument.
  call get_command_argument(1, value = radical, status = istat)
  if (istat > 0) then
     write(radical, "(A)") "input"
  end if


  !initalise the variables for the calculation
  call standard_inputfile_names(in,radical,nproc)
  call read_input_variables(iproc,'posinp',in, atoms, rxyz)

  if (iproc == 0) then
     call print_general_parameters(in,atoms)
  end if
       
  allocate(radii_cf(atoms%ntypes,3+ndebug),stat=i_stat)
  call memocc(i_stat,radii_cf,'radii_cf',subname)

! Nullify paw objects:
  allocate(proj_G(atoms%ntypes),stat=i_stat)
  !call memocc(i_stat,proj_G,'proj_G',subname)
  do iatyp=1,atoms%ntypes
     call nullify_gaussian_basis(proj_G(iatyp))
  end do
  paw%usepaw=0 !Not using PAW
  call nullify_paw_objects(paw)

  call system_properties(iproc,nproc,in,atoms,orbs,radii_cf,nelec)

  ! Determine size alat of overall simulation cell and shift atom positions
  ! then calculate the size in units of the grid space
  call system_size(iproc,atoms,rxyz,radii_cf,in%crmult,in%frmult,in%hx,in%hy,in%hz,&
       Glr,shift)

  !variables substitution for the PSolver part
  hxh=0.5d0*in%hx
  hyh=0.5d0*in%hy
  hzh=0.5d0*in%hz
  n1i=Glr%d%n1i
  n2i=Glr%d%n2i
  n3i=Glr%d%n3i

  n1=Glr%d%n1
  n2=Glr%d%n2
  n3=Glr%d%n3

  call timing(iproc,'CrtDescriptors','ON')
  call createWavefunctionsDescriptors(iproc,in%hx,in%hy,in%hz,&
       atoms,rxyz,radii_cf,in%crmult,in%frmult,Glr)
  call timing(iproc,'CrtDescriptors','OF')
  ! Calculate all projectors, or allocate array for on-the-fly calculation
  call timing(iproc,'CrtProjectors ','ON')
  call createProjectorsArrays(iproc,Glr,rxyz,atoms,orbs,&
       radii_cf,in%frmult,in%frmult,in%hx,in%hy,in%hz,nlpspd,proj_G,proj)
  call timing(iproc,'CrtProjectors ','OF')

  !allocate communications arrays
  call orbitals_communicators(iproc,nproc,Glr,orbs,comms)  

  call check_linear_and_create_Lzd(iproc,nproc,in%linear,Lzd,atoms,orbs,in%nspin,rxyz)

  allocate(nscatterarr(0:nproc-1,4+ndebug),stat=i_stat)
  call memocc(i_stat,nscatterarr,'nscatterarr',subname)
  allocate(ngatherarr(0:nproc-1,2+ndebug),stat=i_stat)
  call memocc(i_stat,ngatherarr,'ngatherarr',subname)
!oneatom is to be rewritten almost completely
!!$  call createDensPotDescriptors(iproc,nproc,atoms,Glr%d,hxh,hyh,hzh,&
!!$       rxyz,in%crmult,in%frmult,radii_cf,in%nspin,'D',0,in%rho_commun,&
!!$       n3d,n3p,n3pi,i3xcsh,i3s,nscatterarr,ngatherarr,rhodsc)

  call local_potential_dimensions(Lzd,orbs,ngatherarr(0,1))
  !commented out, to be used in the future
  if (dokernel) then
     ndegree_ip=16 !default value 
     call createKernel(iproc,nproc,atoms%geocode,n1i,n2i,n3i,hxh,hyh,hzh,ndegree_ip,&
          pkernel,.true.,0)
  else
     nullify(pkernel)
  end if

  !allocate ionic potential
  if (n3pi > 0) then
     allocate(pot_ion(n1i*n2i*n3pi+ndebug),stat=i_stat)
     call memocc(i_stat,pot_ion,'pot_ion',subname)
  else
     allocate(pot_ion(1+ndebug),stat=i_stat)
     call memocc(i_stat,pot_ion,'pot_ion',subname)
  end if

  allocate(psi(orbs%npsidim_orbs+ndebug),stat=i_stat)
  call memocc(i_stat,psi,'psi',subname)
  allocate(hpsi(orbs%npsidim_orbs+ndebug),stat=i_stat)
  call memocc(i_stat,hpsi,'hpsi',subname)
  if (nproc == 1) then
     nullify(psit)
  else
     stop 'for the moment only sequential runs'
  end if

  ! allocate arrays necessary for DIIS convergence acceleration
  !the allocation with npsidim is not necessary here since DIIS arrays
  !are always calculated in the transposed form
  call allocate_diis_objects(in%idsx,in%alphadiis,sum(comms%ncntt(0:nproc-1)),&
       orbs%nkptsp,orbs%nspinor,diis,subname)  

  !write the local potential in pot_ion array
  call createPotential(atoms%geocode,iproc,nproc,atoms,rxyz,hxh,hyh,hzh,&
       n1,n2,n3,n3pi,i3s+i3xcsh,n1i,n2i,n3i,pkernel,pot_ion,0.0_dp) !n(m)

  !pot_ion=0.0d0

  !create input guess wavefunction
  call psi_from_gaussians(iproc,nproc,atoms,orbs,Lzd,rxyz,in%hx,in%hy,in%hz,in%nspin,psi)

!!$  psi=0.0d0
!!$  ttsum=0.0d0
!!$  do i=1,orbs%npsidim
!!$     do j=0,iproc-1
!!$        call random_number(ttr)
!!$     end do
!!$     call random_number(ttr)
!!$     psi(i)=real(ttr,wp)*0.01_wp
!!$     ttsum=ttsum+psi(i)
!!$     do j=iproc+1,nproc
!!$        call random_number(ttr)
!!$     end do
!!$  end do
!!$
!!$  psi=1.d0


  call plot_wf_oneatom('iter0',1,atoms,Glr,hxh,hyh,hzh,rxyz,psi)


  !othogonalise them
  !transpose the psi wavefunction
  call transpose_v(iproc,nproc,orbs,Glr%wfd,comms,&
       psi,work=hpsi)
  call orthogonalize(iproc,nproc,orbs,comms,psi,in%orthpar)
  !untranspose psi
  call untranspose_v(iproc,nproc,orbs,Glr%wfd,comms,psi,work=hpsi)

  allocate(orbs%eval(orbs%norb*orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,orbs%eval,'eval',subname)

  orbs%eval(1:orbs%norb*orbs%nkpts)=-0.5d0

  !n(c) alpha=2.d0
  energy=1.d10
  gnrm=1.d10
  !n(c) gnrm_zero=0.0_gp

  !number of switching betweed DIIS and SD during self-consistent loop
  ndiis_sd_sw=0
  !previous value of idsx_actual to control if switching has appeared
  idsx_actual_before=diis%idsx

  !allocate the potential in the full box
  call full_local_potential(iproc,nproc,orbs,Lzd,0,denspotd,pot_ion,pot)
!!$  call full_local_potential(iproc,nproc,Glr%d%n1i*Glr%d%n2i*n3p,Glr%d%n1i*Glr%d%n2i*Glr%d%n3i,&
!!$       in%nspin,&
!!$       Glr%d%n1i*Glr%d%n2i*n3d*in%nspin,0,&
!!$       orbs,Lzd,0,ngatherarr,pot_ion,pot)
!!$  
  allocate(confdatarr(orbs%norbp))
  call default_confinement_data(confdatarr,orbs%norbp)

  wfn_loop: do iter=1,in%itermax

     if (iproc == 0 .and. verbose > 0) then 
        write( *,'(1x,a,i0)') &
             & repeat('~',76 - int(log(real(iter))/log(10.))) // ' iter= ', iter
     endif
     !control whether the minimisation iterations ended
     endloop= gnrm <= in%gnrm_cv .or. iter == in%itermax
     
     !control how many times the DIIS has switched into SD
     if (idsx_actual /= idsx_actual_before) ndiis_sd_sw=ndiis_sd_sw+1

     !terminate SCF loop if forced to switch more than once from DIIS to SD
     endloop=endloop .or. ndiis_sd_sw > 2

     call FullHamiltonianApplication(iproc,nproc,atoms,orbs,rxyz,&
          proj,Lzd,nlpspd,confdatarr,ngatherarr,pot_ion,psi,hpsi,&
          energs,in%SIC,GPU)

     call total_energies(energs, iter, iproc)
     energy=energs%eKS

     !check for convergence or whether max. numb. of iterations exceeded
     if (endloop) then 
        if (iproc == 0) then 
           if (verbose > 1) write( *,'(1x,a,i0,a)')'done. ',iter,' minimization iterations required'
           write( *,'(1x,a)') &
                '------------------------------------------- End of Virtual Wavefunction Optimisation'
              if ((in%itrpmax >1 .and. endlooprp) .or. in%itrpmax == 1) then

                 call write_energies(iter,0,energs,gnrm,0.0_gp,"FINAL")
                 
                 !write(61,*)hx,hy,hz,energy,ekin_sum,epot_sum,eproj_sum,ehart,eexcu,vexcu
                 if (energy > energy_min) write( *,'(1x,a,1pe9.2)')&
                      'WARNING: Found an "energy" value lower than the FINAL "energy", delta:',energy-energy_min
              end if
           end if
           exit wfn_loop 
        endif
     !control the previous value of idsx_actual
     idsx_actual_before=idsx_actual
     stop
     !call hpsitopsi(iproc,nproc,orbs,Glr,comms,iter,diis,in%idsx,psi,psit,hpsi,in%orthpar) 

     write(itername,'(i4.4)')iter
     call plot_wf_oneatom('iter'//itername,1,atoms,Glr,hxh,hyh,hzh,rxyz,psi)

     tt=(energs%ebs-scprsum)/scprsum
     if (((abs(tt) > 1.d-10 .and. .not. GPUconv) .or.&
          (abs(tt) > 1.d-8 .and. GPUconv)) .and. iproc==0) then 
        write( *,'(1x,a,1pe9.2,2(1pe22.14))') &
             'ERROR: inconsistency between gradient and energy',tt,energs%ebs,scprsum
     endif
     if (iproc.eq.0) then
        if (verbose > 0) then
           call write_energies(iter,0,energs,gnrm,0.0_gp,"FINAL")
           !write( *,'(1x,a,3(1x,1pe18.11))') 'ekin_sum,epot_sum,eproj_sum',  & 
           !     ekin_sum,epot_sum,eproj_sum
        end if
        write( *,'(1x,a,i6,2x,1pe24.17,1x,1pe9.2)') 'iter,total "energy",gnrm',iter,energy,gnrm
     endif
  end do wfn_loop
  if (iter == in%itermax .and. iproc == 0 ) &
       write( *,'(1x,a)')'No convergence within the allowed number of minimization steps'

  !this deallocates also hpsivirt and psitvirt
  stop
  !call last_orthon(iproc,nproc,orbs,Glr%wfd,in%nspin,&
  !     comms,psi,hpsi,psit,evsum)
  
  call deallocate_diis_objects(diis,subname)

  !i_all=-product(shape(proj_G))*kind(proj_G)
  deallocate(proj_G,stat=i_stat)
  !call memocc(i_stat,i_all,'proj_G',subname)

  if (nproc > 1) then
     i_all=-product(shape(psit))*kind(psit)
     deallocate(psit,stat=i_stat)
     call memocc(i_stat,i_all,'psit',subname)
  end if

  i_all=-product(shape(psi))*kind(psi)
  deallocate(psi,stat=i_stat)
  call memocc(i_stat,i_all,'psi',subname)

  i_all=-product(shape(proj))*kind(proj)
  deallocate(proj,stat=i_stat)
  call memocc(i_stat,i_all,'proj',subname)

  i_all=-product(shape(orbs%eval))*kind(orbs%eval)
  deallocate(orbs%eval,stat=i_stat)
  call memocc(i_stat,i_all,'eval',subname)


  i_all=-product(shape(nscatterarr))*kind(nscatterarr)
  deallocate(nscatterarr,stat=i_stat)
  call memocc(i_stat,i_all,'nscatterarr',subname)
  i_all=-product(shape(ngatherarr))*kind(ngatherarr)
  deallocate(ngatherarr,stat=i_stat)
  call memocc(i_stat,i_all,'ngatherarr',subname)

  i_all=-product(shape(radii_cf))*kind(radii_cf)
  deallocate(radii_cf,stat=i_stat)
  call memocc(i_stat,i_all,'radii_cf',subname)


  call deallocate_lr(Glr,subname)
  call deallocate_comms(comms,subname)
  call deallocate_orbs(orbs,subname)
  call deallocate_proj_descr(nlpspd,subname)

  if (dokernel) then
     i_all=-product(shape(pkernel))*kind(pkernel)
     deallocate(pkernel,stat=i_stat)
     call memocc(i_stat,i_all,'kernel',subname)
  end if

  !deallocate potential
  call free_full_potential(nproc,0,pot,subname)

  i_all=-product(shape(pot_ion))*kind(pot_ion)
  deallocate(pot_ion,stat=i_stat)
  call memocc(i_stat,i_all,'pot_ion',subname)
  
  call deallocate_atoms(atoms,subname) 
  i_all=-product(shape(rxyz))*kind(rxyz)
  deallocate(rxyz,stat=i_stat)
  call memocc(i_stat,i_all,'rxyz',subname)
  call free_input_variables(in)

  call deallocate_rho_descriptors(rhodsc,subname)
  deallocate(confdatarr)
  !finalize memory counting
  call memocc(0,0,'count','stop')


end program oneatom



!>
!!
!!
subroutine createPotential(geocode,iproc,nproc,at,rxyz,& !n(c) elecfield (arg:9)
     hxh,hyh,hzh,n1,n2,n3,n3pi,i3s,n1i,n2i,n3i,pkernel,pot_ion,psoffset)
  use module_base
  use module_types
  !  use module_interfaces, except_this_one => createIonicPotential
  use Poisson_Solver
  implicit none
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  integer, intent(in) :: iproc,nproc,n1,n2,n3,n3pi,i3s,n1i,n2i,n3i
  real(gp), intent(in) :: hxh,hyh,hzh,psoffset
  type(atoms_data), intent(in) :: at
  !n(c) real(gp), intent(in) :: elecfield(3)
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(dp), dimension(*), intent(in) :: pkernel
  real(wp), dimension(*), intent(inout) :: pot_ion
  !local variables
  character(len=*), parameter :: subname='createPotential'
  logical :: perx,pery,perz,gox,goy,goz
  logical :: htoobig=.false.
  logical :: check_potion=.false.
  logical :: harmonic=.true.
  integer :: iat,i1,i2,i3,j1,j2,j3,isx,isy,isz,iex,iey,iez,ierr,ityp,nspin
  integer :: nl1,nl2,nl3,nu1,nu2,nu3
  integer :: ind,i_all,i_stat,nbl1,nbr1,nbl2,nbr2,nbl3,nbr3,nloc,iloc
  real(kind=8) :: pi,rholeaked,rloc,charge,cutoff,x,y,z,r2,arg,xp,tt,rx,ry,rz
  real(kind=8) :: tt_tot,rholeaked_tot,potxyz,potcoeff
  real(wp) :: maxdiff
  real(gp) :: ehart
  real(dp), dimension(2) :: charges_mpi
  real(dp), dimension(:), allocatable :: potion_corr

  read(1,*)potcoeff

  call timing(iproc,'CrtLocPot     ','ON')

  if (iproc.eq.0) then
     write(*,'(1x,a)')&
          '----------------------------------------------------------------- Potential Creation'
  end if

  pi=4.d0*atan(1.d0)
  ! Ionic charge (must be calculated for the PS active processes)
  rholeaked=0.d0
  ! Ionic energy (can be calculated for all the processors)

  !Creates charge density arising from the ionic PSP cores
  call razero(n1i*n2i*n3pi,pot_ion)

  !conditions for periodicity in the three directions
  perx=(geocode /= 'F')
  pery=(geocode == 'P')
  perz=(geocode /= 'F')

  call ext_buffers(perx,nbl1,nbr1)
  call ext_buffers(pery,nbl2,nbr2)
  call ext_buffers(perz,nbl3,nbr3)

  if (harmonic) then
     if (n3pi >0) then

        do iat=1,at%nat
           ityp=at%iatype(iat)
           rx=rxyz(1,iat) 
           ry=rxyz(2,iat)
           rz=rxyz(3,iat)

           rloc=at%psppar(0,0,ityp)
           charge=real(at%nelpsp(ityp),kind=8)/(2.d0*pi*sqrt(2.d0*pi)*rloc**3)
           cutoff=max(at%alat1,at%alat2,at%alat3)!10.d0*rloc

           isx=floor((rx-cutoff)/hxh)
           isy=floor((ry-cutoff)/hyh)
           isz=floor((rz-cutoff)/hzh)

           iex=ceiling((rx+cutoff)/hxh)
           iey=ceiling((ry+cutoff)/hyh)
           iez=ceiling((rz+cutoff)/hzh)

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
                    arg=r2*potcoeff
                    xp=exp(-.5d0*arg)
                    if (j3 >= i3s .and. j3 <= i3s+n3pi-1  .and. goy  .and. gox ) then
                       ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s+1-1)*n1i*n2i
                       pot_ion(ind)=pot_ion(ind)+arg
                    else if (.not. goz ) then
                       rholeaked=rholeaked+arg
                    endif
                 enddo
              enddo
           enddo

        enddo
     end if
  else
     if (n3pi >0 .and. .not. htoobig) then

        do iat=1,at%nat
           ityp=at%iatype(iat)
           rx=rxyz(1,iat) 
           ry=rxyz(2,iat)
           rz=rxyz(3,iat)

           rloc=at%psppar(0,0,ityp)
           charge=real(at%nelpsp(ityp),kind=8)/(2.d0*pi*sqrt(2.d0*pi)*rloc**3)
           cutoff=10.d0*rloc

           isx=floor((rx-cutoff)/hxh)
           isy=floor((ry-cutoff)/hyh)
           isz=floor((rz-cutoff)/hzh)

           iex=ceiling((rx+cutoff)/hxh)
           iey=ceiling((ry+cutoff)/hyh)
           iez=ceiling((rz+cutoff)/hzh)

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
                       pot_ion(ind)=pot_ion(ind)-xp*charge
                    else if (.not. goz ) then
                       rholeaked=rholeaked+xp*charge
                    endif
                 enddo
              enddo
           enddo

        enddo

     end if

     ! Check
     tt=0.d0
     do j3=1,n3pi
        do i2= -nbl2,2*n2+1+nbr2
           do i1= -nbl1,2*n1+1+nbr1
              ind=i1+1+nbl1+(i2+nbl2)*n1i+(j3-1)*n1i*n2i
              tt=tt+pot_ion(ind)
           enddo
        enddo
     enddo

     tt=tt*hxh*hyh*hzh
     rholeaked=rholeaked*hxh*hyh*hzh

     !print *,'test case input_rho_ion',iproc,i3start,i3end,n3pi,2*n3+16,tt

     if (nproc > 1) then
        charges_mpi(1)=tt
        charges_mpi(2)=rholeaked

        call mpiallred(charges_mpi(1),2,MPI_SUM,MPI_COMM_WORLD,ierr)

        tt_tot=charges_mpi(1)
        rholeaked_tot=charges_mpi(2)
     else
        tt_tot=tt
        rholeaked_tot=rholeaked
     end if

     if (iproc == 0) write(*,'(1x,a,f26.12,2x,1pe10.3)') &
          'total ionic charge, leaked charge ',tt_tot,rholeaked_tot

     if (.not. htoobig) then
        call timing(iproc,'CrtLocPot     ','OF')
        !here the value of the datacode must be kept fixed
        nspin=1

        call H_potential(geocode,'D',iproc,nproc,&
             n1i,n2i,n3i,hxh,hyh,hzh,&
             pot_ion,pkernel,pot_ion,ehart,-psoffset,.false.)

        call timing(iproc,'CrtLocPot     ','ON')

        if (check_potion) then
           if (iproc == 0) write(*,'(1x,a)',advance='no') &
                'Check the ionic potential...'

           allocate(potion_corr(n1i*n2i*n3pi+ndebug),stat=i_stat)
           call memocc(i_stat,potion_corr,'potion_corr',subname)

           call razero(n1i*n2i*n3pi,potion_corr)

           !calculate pot_ion with an explicit error function to correct in the case of big grid spacings
           !for the moment works only in the isolated BC case
           do i3=1,n3pi
              z=real(i3+i3s-1-nbl3-1,gp)*hzh
              do i2=1,n2i
                 y=real(i2-nbl2-1,gp)*hyh
                 do i1=1,n1i
                    x=real(i1-nbl1-1,gp)*hxh
                    ind=i1+(i2-1)*n1i+(i3-1)*n1i*n2i
                    !if (i1==49 .and. i2==46 .and. i3==44) then
                    call sum_erfcr(at%nat,at%ntypes,x,y,z,at%iatype,at%nelpsp,at%psppar,rxyz,potxyz)
                    !   stop
                    !end if
                    potion_corr(ind)=potion_corr(ind)+potxyz
                    !write(18,'(3(i6),i12,3(1x,1pe24.17))')i1,i2,i3,ind,potion_corr(ind),pot_ion(ind)
                 end do
              end do
           end do

           !then calculate the maximum difference in the sup norm
           maxdiff=0.0_wp
           do i3=1,n3pi
              do i2=1,n2i
                 do i1=1,n1i
                    ind=i1+(i2-1)*n1i+(i3-1)*n1i*n2i
                    maxdiff=max(maxdiff,abs(potion_corr(ind)-pot_ion(ind)))
                    !write(17,'(3(i6),i12,3(1x,1pe24.17))')i1,i2,i3,ind,potion_corr(ind),pot_ion(ind),maxdiff
                 end do
              end do
           end do

           !call mpiallred(maxdiff,1,MPI_MAX,MPI_COMM_WORLD,ierr)

           if (iproc == 0) write(*,'(1x,a,1pe24.17)')'...done. MaxDiff=',maxdiff

           stop

           i_all=-product(shape(potion_corr))*kind(potion_corr)
           deallocate(potion_corr,stat=i_stat)
           call memocc(i_stat,i_all,'potion_corr',subname)

        end if

     end if


!!!  !calculate the value of the offset to be put
!!!  tt_tot=0.d0
!!!  do ind=1,n1i*n2i*n3i
!!!     tt_tot=tt_tot+pot_ion(ind)
!!!  end do
!!!  print *,'previous offset',tt_tot*hxh*hyh*hzh

     if (n3pi > 0) then
        do iat=1,at%nat
           ityp=at%iatype(iat)

           rx=rxyz(1,iat)
           ry=rxyz(2,iat)
           rz=rxyz(3,iat)

           ! determine number of local terms
           nloc=0
           do iloc=1,4
              if (at%psppar(0,iloc,ityp) /= 0.d0) nloc=iloc
           enddo
           rloc=at%psppar(0,0,ityp)
           cutoff=10.d0*rloc

           isx=floor((rx-cutoff)/hxh)
           isy=floor((ry-cutoff)/hyh)
           isz=floor((rz-cutoff)/hzh)

           iex=ceiling((rx+cutoff)/hxh)
           iey=ceiling((ry+cutoff)/hyh)
           iez=ceiling((rz+cutoff)/hzh)

           !do not add the local part for the vacancy
           if (nloc /= 0) then

              do i3=isz,iez
                 z=real(i3,kind=8)*hzh-rz
                 call ind_positions(perz,i3,n3,j3,goz) 
                 j3=j3+nbl3+1
                 if (goz .and. j3 >= i3s .and. j3 <=  i3s+n3pi-1) then
                    do i2=isy,iey
                       y=real(i2,kind=8)*hyh-ry
                       call ind_positions(pery,i2,n2,j2,goy)
                       if (goy) then
                          do i1=isx,iex
                             x=real(i1,kind=8)*hxh-rx
                             call ind_positions(perx,i1,n1,j1,gox)
                             if (gox) then
                                r2=x**2+y**2+z**2
                                arg=r2/rloc**2
                                xp=exp(-.5d0*arg)
                                tt=at%psppar(0,nloc,ityp)
                                do iloc=nloc-1,1,-1
                                   tt=arg*tt+at%psppar(0,iloc,ityp)
                                enddo
                                ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s+1-1)*n1i*n2i
                                pot_ion(ind)=pot_ion(ind)+xp*tt
                             end if
                          enddo
                       end if
                    enddo
                 end if
              end do

           end if

        enddo

        if (htoobig) then
           !add to pot_ion an explicit error function to correct in the case of big grid spacing
           !for the moment works only in the isolated BC case
           do i3=1,n3pi
              z=real(i3+i3s-1-nbl3-1,gp)*hzh
              do i2=1,n2i
                 y=real(i2-nbl2-1,gp)*hyh
                 do i1=1,n1i
                    x=real(i1-nbl1-1,gp)*hxh
                    ind=i1+(i2-1)*n1i+(i3-1)*n1i*n2i
                    call sum_erfcr(at%nat,at%ntypes,x,y,z,at%iatype,at%nelpsp,at%psppar,rxyz,potxyz)
                    pot_ion(ind)=pot_ion(ind)+potxyz
                 end do
              end do
           end do
        end if

     end if

  end if
  call timing(iproc,'CrtLocPot     ','OF')


  !plot the created potential (valid only with nproc =1 and geocode == 'F')
     nl1=14
     nu1=15
     nl2=14
     nu2=15
     nl3=14
     nu3=15

  do iat=1,at%nat

     open(unit=22,file='potential',status='unknown')

     do i3=-nl3,2*n3+1+nu3
        z=hzh*real(i3,gp)
        do i2=-nl2,2*n2+1+nu2
           y=hyh*real(i2,gp)
           do i1=-nl1,2*n1+1+nu1
              x=hxh*real(i1,gp)
              if (z == rxyz(3,iat) .and. y ==  rxyz(2,iat)) then
                 !print *,'value of the center',rxyz(1,iat)
                 ind=i1+nl1+1+(i2+nl2)*n1i+(i3+nl3)*n1i*n2i
                 write(22,'(1x,f9.5,1pe18.10)')x,pot_ion(ind)
              end if
           end do
        end do
     end do
     close(22)
  end do


END SUBROUTINE createPotential



!>
!!
!!
subroutine psi_from_gaussians(iproc,nproc,at,orbs,Lzd,rxyz,hx,hy,hz,nspin,psi)
  use module_base
  use module_types
  use module_interfaces
  implicit none
  integer, intent(in) :: iproc,nproc,nspin
  real(gp), intent(in) :: hx,hy,hz
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  type(local_zone_descriptors), intent(in) :: Lzd
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(wp), dimension(orbs%npsidim_orbs), intent(out) :: psi
  !local variables
  character(len=*), parameter :: subname='psi_from_gaussians'
  logical ::  randinp
  integer :: iorb,icoeff,i_all,i_stat,jproc,nwork,info,jorb
  real(kind=4) :: tt
  real(wp), dimension(:,:), allocatable :: gaucoeffs
  real(gp), dimension(:), allocatable :: work,ev
  real(gp), dimension(:,:), allocatable :: ovrlp
  type(gaussian_basis) :: G
  real(wp), dimension(:), pointer :: gbd_occ


  !initialise some coefficients in the gaussian basis
  !nullify the G%rxyz pointer
  nullify(G%rxyz)
  !extract the gaussian basis from the pseudowavefunctions
  !use a better basis than the input guess
  call gaussian_pswf_basis(31,.false.,iproc,nspin,at,rxyz,G,gbd_occ)

  allocate(gaucoeffs(G%ncoeff,orbs%norbp*orbs%nspinor+ndebug),stat=i_stat)
  call memocc(i_stat,gaucoeffs,'gaucoeffs',subname)

  !the kinetic overlap is correctly calculated only with Free BC
  randinp = .true.!.false.!lr%geocode /= 'F'

  if (randinp) then
     !fill randomly the gaussian coefficients for the orbitals considered
     do iorb=1,orbs%norbp*orbs%nspinor
        do icoeff=1,G%ncoeff
           !be sure to call always a different random number
           do jproc=0,iproc-1
              call random_number(tt)
           end do
           call random_number(tt)
           gaucoeffs(icoeff,iorb)=real(tt,wp)
           do jproc=iproc+1,nproc-1
              call random_number(tt)
           end do
        end do
     end do

     !othogonalise the gaussian basis (wrong with k-points)
     !call gaussian_orthogonality(iproc,nproc,norb,norbp,G,coeffs)

  else
     !as an alternative strategy we may take the eigenvectors of the kinetic+k hamiltonian

     !in view of complete gaussian calculation
     allocate(ovrlp(G%ncoeff,G%ncoeff),stat=i_stat)
     call memocc(i_stat,ovrlp,'ovrlp',subname)

     !overlap calculation of the kinetic operator, upper triangular part
     !call kinetic_overlap(G,G,ovrlp)
     call gaussian_overlap(G,G,ovrlp)
     nwork=3*G%ncoeff+1
     allocate(work(nwork+ndebug),stat=i_stat)
     call memocc(i_stat,work,'work',subname)
     allocate(ev(G%ncoeff+ndebug),stat=i_stat)
     call memocc(i_stat,ev,'ev',subname)

!!$  if (iproc == 0) then
!!$     do iat=1,G%ncoeff
!!$        write(*,'(a,i0,10(1pe15.7))')'T',iat,ovrlp(1:iat,iat)
!!$     end do
!!$  end if

     !print *,'nwork',nwork,3*nbasis-1
     call dsyev('V','U',G%ncoeff,ovrlp(1,1),G%ncoeff,ev(1),work(1),nwork,info)
     if (info /= 0) then
        if (iproc == 0) then
           write(*,*)'DSyev Error',info
        end if
        stop
     end if

!!$  if (iproc == 0) then
!!$     do iat=1,G%ncoeff
!!$        write(*,'(a,i0,10(1pe15.7))')'Ev',iat,ovrlp(:,iat)
!!$     end do
!!$     do iat=1,G%ncoeff
!!$        write(*,'(a,i0,10(1pe15.7))')'K',iat,ev(iat)
!!$     end do
!!$  end if

     !copy the eigenvectors to the matrix
     call razero(G%ncoeff*orbs%norbp*orbs%nspinor,gaucoeffs)
     if (orbs%norb > G%ncoeff) stop 'wrong gaussian basis'
     jorb=mod(orbs%isorb,orbs%norb)
     do iorb=1,orbs%norbp
        jorb=jorb+1
        if (jorb == orbs%norb+1) jorb=1 !for k-points calculation
        call vcopy(G%ncoeff,ovrlp(1,jorb),1,gaucoeffs(1,orbs%nspinor*(iorb-1)+1),orbs%nspinor)
     end do


     i_all=-product(shape(ovrlp))*kind(ovrlp)
     deallocate(ovrlp,stat=i_stat)
     call memocc(i_stat,i_all,'ovrlp',subname)
     i_all=-product(shape(work))*kind(work)
     deallocate(work,stat=i_stat)
     call memocc(i_stat,i_all,'work',subname)
     i_all=-product(shape(ev))*kind(ev)
     deallocate(ev,stat=i_stat)
     call memocc(i_stat,i_all,'ev',subname)

     !call MPI_BARRIER(MPI_COMM_WORLD,info)
     !stop

  end if

  call gaussians_to_wavelets_new(iproc,nproc,Lzd,orbs,G,gaucoeffs,psi)

  !deallocate the gaussian basis descriptors
  call deallocate_gwf(G,subname)

  !deallocate gaussian array
  i_all=-product(shape(gaucoeffs))*kind(gaucoeffs)
  deallocate(gaucoeffs,stat=i_stat)
  call memocc(i_stat,i_all,'gaucoeffs',subname)
  i_all=-product(shape(gbd_occ))*kind(gbd_occ)
  deallocate(gbd_occ,stat=i_stat)
  call memocc(i_stat,i_all,'gbd_occ',subname)

  
END SUBROUTINE psi_from_gaussians


subroutine plot_wf_oneatom(orbname,nexpo,at,lr,hxh,hyh,hzh,rxyz,psi)
  use module_base
  use module_types
  implicit none
  character(len=*) :: orbname
  integer, intent(in) :: nexpo
  real(gp), intent(in) :: hxh,hyh,hzh
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  type(locreg_descriptors), intent(in) :: lr
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f), intent(in) :: psi
  !local variables
  character(len=*), parameter :: subname='plot_wf_oneatom'
  integer :: i_stat,i_all
  integer :: nl1,nl2,nl3,n1i,n2i,n3i,n1,n2,n3,i1,i2,i3,nu1,nu2,nu3,iat
  real(gp) :: x,y,z,maxval
  type(workarr_sumrho) :: w
  real(wp), dimension(:), allocatable :: psi2
  real(wp), dimension(:,:,:), allocatable :: psir,psir2


  n1=lr%d%n1
  n2=lr%d%n2
  n3=lr%d%n3
  n1i=lr%d%n1i
  n2i=lr%d%n2i
  n3i=lr%d%n3i

  if (at%geocode == 'F') then
     nl1=14
     nu1=15
     nl2=14
     nu2=15
     nl3=14
     nu3=15
  else if (at%geocode == 'S') then
     nl1=0
     nu1=0
     nl2=14
     nu2=15
     nl3=0
     nu3=0
  else if (at%geocode == 'P') then
     nl1=0
     nu1=0
     nl2=0
     nu2=0
     nl3=0
     nu3=0
  end if

  call initialize_work_arrays_sumrho(lr,w)
 
  allocate(psir(-nl1:2*n1+1+nu1,-nl2:2*n2+1+nu2,-nl3:2*n3+1+nu3+ndebug),stat=i_stat)
  call memocc(i_stat,psir,'psir',subname)
  !initialisation
  if (lr%geocode == 'F') then
     call razero(lr%d%n1i*lr%d%n2i*lr%d%n3i,psir)
  end if

  allocate(psi2(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f+ndebug),stat=i_stat)
  call memocc(i_stat,psi2,'psi2',subname)
  allocate(psir2(-nl1:2*n1+1+nu1,-nl2:2*n2+1+nu2,-nl3:2*n3+1+nu3+ndebug),stat=i_stat)
  call memocc(i_stat,psir2,'psir2',subname)

  print *,'normDaub',dot(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,psi(1),1,psi(1),1)
  
  call daub_to_isf(lr,w,psi,psir)

  print *,'normISF',dot(n1i*n2i*n3i,psir(-nl1,-nl2,-nl3),1,psir(-nl1,-nl2,-nl3),1)

  call vcopy(n1i*n2i*n3i,psir(-nl1,-nl2,-nl3),1,psir2(-nl1,-nl2,-nl3),1)
  call isf_to_daub(lr,w,psir2,psi2)
  !psir2 is destroyed
  call vcopy(n1i*n2i*n3i,psir(-nl1,-nl2,-nl3),1,psir2(-nl1,-nl2,-nl3),1)

  print *,'normDaub2',dot(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,psi2(1),1,psi2(1),1)

  call daub_to_isf(lr,w,psi2,psir)

  print *,'normISF2',dot(n1i*n2i*n3i,psir(-nl1,-nl2,-nl3),1,psir(-nl1,-nl2,-nl3),1)

  maxval=0.0_wp
  do i3=-nl3,2*n3+1+nu3
     do i2=-nl2,2*n2+1+nu2
        do i1=-nl1,2*n1+1+nu1
           maxval=max(maxval,abs(psir(i1,i2,i3)-psir2(i1,i2,i3)))
         end do
     end do
  end do
  print *,'maxvalISF',maxval

  maxval=0.0_wp
  do i3=1,lr%wfd%nvctr_c+7*lr%wfd%nvctr_f
     maxval=max(maxval,abs(psi(i3)-psi2(i3)))
  end do
  print *,'maxvalDaub',maxval

  i_all=-product(shape(psir2))*kind(psir2)
  deallocate(psir2,stat=i_stat)
  call memocc(i_stat,i_all,'psir2',subname)

  i_all=-product(shape(psi2))*kind(psi2)
  deallocate(psi2,stat=i_stat)
  call memocc(i_stat,i_all,'psi2',subname)

  stop

  do iat=1,at%nat

     open(unit=22,file=trim(orbname),status='unknown')

     do i3=-nl3,2*n3+1+nu3
        z=hzh*real(i3,gp)
        do i2=-nl2,2*n2+1+nu2
           y=hyh*real(i2,gp)
           do i1=-nl1,2*n1+1+nu1
              x=hxh*real(i1,gp)
              if (z == rxyz(3,iat) .and. y ==  rxyz(2,iat)) then
                 !print *,'value of the center',rxyz(1,iat)
                 write(22,'(1x,f9.5,1pe18.10)')x,psir(i1,i2,i3)**nexpo
              end if
           end do
        end do
     end do
     close(22)
  end do

  i_all=-product(shape(psir))*kind(psir)
  deallocate(psir,stat=i_stat)
  call memocc(i_stat,i_all,'psir',subname)

  call deallocate_work_arrays_sumrho(w)

END SUBROUTINE plot_wf_oneatom

