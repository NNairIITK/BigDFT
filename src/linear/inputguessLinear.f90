!>   input guess wavefunction diagonalization
subroutine inputguessConfinement(iproc, nproc, at, &
     comms, Glr, input, lin, rxyz, n3p, rhopot, rhocore, pot_ion,&
     nlpspd, proj, pkernel, pkernelseq, &
     nscatterarr, ngatherarr, potshortcut, irrzon, phnons, GPU, radii_cf,  &
     phi, ehart, eexcu, vexcu)
  ! Input wavefunctions are found by a diagonalization in a minimal basis set
  ! Each processors write its initial wavefunctions into the wavefunction file
  ! The files are then read by readwave
  use module_base
  use module_interfaces, exceptThisOne => inputguessConfinement
  use module_types
  use Poisson_Solver
  implicit none
  !Arguments
  integer, intent(in) :: iproc,nproc,n3p
  type(atoms_data), intent(inout) :: at
  type(nonlocal_psp_descriptors), intent(in) :: nlpspd
  type(locreg_descriptors), intent(in) :: Glr
  type(communications_arrays), intent(in) :: comms
  type(GPU_pointers), intent(inout) :: GPU
  type(input_variables):: input
  type(linearParameters),intent(inout):: lin
  integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
  real(dp),dimension(max(Glr%d%n1i*Glr%d%n2i*n3p,1)*input%nspin),intent(inout) :: rhopot
  real(wp), dimension(lin%as%size_pot_ion),intent(inout):: pot_ion
  real(wp), dimension(:), pointer :: rhocore
  real(dp), dimension(lin%as%size_pkernel),intent(in):: pkernel
  real(dp), dimension(:), pointer :: pkernelseq
  integer, intent(in) ::potshortcut
  integer, dimension(lin%as%size_irrzon(1),lin%as%size_irrzon(2),lin%as%size_irrzon(3)),intent(in) :: irrzon
  real(dp), dimension(lin%as%size_phnons(1),lin%as%size_phnons(2),lin%as%size_phnons(3)),intent(in) :: phnons
  real(8),dimension(at%ntypes,3),intent(in):: radii_cf
  real(8),dimension(lin%orbs%npsidim),intent(out):: phi
  real(8),intent(out):: ehart, eexcu, vexcu
  !local variables
  type(gaussian_basis):: G !basis for davidson IG
  character(len=*), parameter :: subname='input_wf_diag'
  logical :: switchGPUconv,switchOCLconv
  integer :: i_stat,i_all,iat,nspin_ig,iorb,idum=0, nvirt, norbat
  real(kind=4) :: tt,builtin_rand
  real(gp) :: hxh,hyh,hzh,eks,epot_sum,ekin_sum,eexctX,eproj_sum,etol,accurex
  type(orbitals_data) :: orbsig
  type(communications_arrays) :: commsig
  integer, dimension(:,:), allocatable :: norbsc_arr
  real(wp), dimension(:), allocatable :: potxc
  real(gp), dimension(:), allocatable :: locrad
  type(locreg_descriptors), dimension(:), allocatable :: Llr
  real(wp), dimension(:), pointer :: pot
  real(wp), dimension(:,:,:), pointer :: psigau
type(linear_zone_descriptors):: lzdig
type(p2pCommsGatherPot):: comgp
type(orbitals_data):: orbsLIN
type(communications_arrays):: commsLIN
real(8),dimension(:),allocatable:: eval
integer:: istat
real(8),dimension(:),allocatable:: chi, lchi
real(8),dimension(:,:),allocatable:: hchi, lhchi
integer,dimension(:),allocatable:: onWhichAtom, onWhichAtomp, norbsPerAt, onWhichAtomTemp, onWhichAtomPhi
integer, parameter :: nmax=6,lmax=3,noccmax=2,nelecmax=32
logical:: withConfinement

integer, dimension(lmax+1) :: nl
real(gp), dimension(noccmax,lmax+1) :: occup
real(8):: dnrm2, ddot, dasum, t1, t2, time
integer:: ist, jst, jorb, iiAt, i, iadd, ii, jj, ndimpot, ilr, ind1, ind2, ldim, gdim, ierr



  allocate(norbsc_arr(at%natsc+1,input%nspin+ndebug),stat=i_stat)
  call memocc(i_stat,norbsc_arr,'norbsc_arr',subname)
  allocate(locrad(at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,locrad,'locrad',subname)
  allocate(norbsPerAt(at%nat), stat=istat)
  call memocc(istat, norbsPerAt, 'norbsPerAt', subname)

  ! number of localization regions
  lzdig%nlr=at%nat


  if (iproc == 0) then
     write(*,'(1x,a)')&
          '------------------------------------------------------- Input Wavefunctions Creation'
  end if

  !spin for inputguess orbitals
  if (input%nspin == 4) then
     nspin_ig=1
  else
     nspin_ig=input%nspin
  end if

  ! Determine how many atomic orbitals we have. Maybe we have to increase this number to more than
  ! its 'natural' value.
  norbat=0
  ist=0
  do iat=1,at%nat
      ii=lin%norbsPerType(at%iatype(iat))
      iadd=0
      do 
          ! Count the number of atomic orbitals and increase the number if necessary until we have more
          ! (or equal) atomic orbitals than basis functions per atom.
          jj=1*nint(at%aocc(1,iat))+3*nint(at%aocc(3,iat))+5*nint(at%aocc(7,iat))+7*nint(at%aocc(13,iat))
          if(jj>=ii) then
              ! we have enough atomic orbitals
              exit
          else
              ! add additional orbitals
              iadd=iadd+1
              select case(iadd)
                  case(1) 
                      at%aocc(1,iat)=1.d0
                  case(2) 
                      at%aocc(3,iat)=1.d0
                  case(3) 
                      at%aocc(7,iat)=1.d0
                  case(4) 
                      at%aocc(13,iat)=1.d0
                  case default 
                      write(*,'(x,a)') 'ERROR: more than 16 basis functions per atom are not possible!'
                      stop
              end select
          end if
      end do
      norbsPerAt(iat)=jj
      norbat=norbat+norbsPerAt(iat)
  end do

  ! Assign the orbitals to the atoms. onWhichAtom(i)=j means that orbital i belongs to atom j.
  ! onWhichAtom is the 'global' distribution and onWhichAtomp is the one for each MPI process.
  allocate(onWhichAtom(norbat),stat=i_stat)
  call memocc(i_stat, onWhichAtom, 'onWhichAtom', subname)

  ! This allocate is garbage (orbsig%norbp is not defined yet!) Anyway this array is not needed...
  !!allocate(onWhichAtomp(orbsig%norbp),stat=i_stat)
  !!call memocc(i_stat, onWhichAtomp, 'onWhichAtomp', subname)
  ist=0
  do iat=1,at%nat
      do i=1,norbsPerAt(iat)
          onWhichAtom(ist+i)=iat
      end do
      ist=ist+norbsPerAt(iat)
  end do
  !write(*,'(a,i3,3x,100i4)') 'iproc, owa', iproc, onWhichAtom(:)


  ! Create the atomic orbitals in a Gaussian basis.
  nvirt=0
  call inputguess_gaussian_orbitals(iproc,nproc,at,rxyz,Glr,nvirt,nspin_ig,&
       lin%orbs,orbsig,norbsc_arr,locrad,G,psigau,eks)
  call inputguess_gaussian_orbitals(iproc,nproc,at,rxyz,Glr,nvirt,nspin_ig,&
       lin%orbs,lzdig%orbs,norbsc_arr,locrad,G,psigau,eks)

  ! Allocate communications arrays for inputguess orbitals.
  call orbitals_communicators(iproc,nproc,Glr,orbsig,commsig)  
  call orbitals_communicators(iproc,nproc,Glr,lzdig%orbs,lzdig%comms)  

  allocate(onWhichAtomp(lzdig%orbs%norbp), stat=istat)
  call memocc(i_stat, onWhichAtomp, 'onWhichAtomp', subname)
  !call assignOrbitalsToAtoms(iproc, lzdig%orbs, at%nat, norbsPerAt, onWhichAtomp)
  call assignOrbitalsToAtoms(iproc, lzdig%orbs, at%nat, norbsPerAt, onWhichAtomp, onWhichAtom)
  ! This is the same as above, but with orbs%inWhichLocreg instead of lin%onWhichAtom
  call assignToLocreg2(iproc, at%nat, lzdig%nlr, input%nspin, norbsPerAt, lzdig%orbs)
  !!write(*,'(a,i3,3x,100i4)') 'iproc, owa', iproc, onWhichAtom(:)
  !!write(*,'(a,i3,3x,100i4)') 'iproc, iwi', iproc, lzdig%orbs%inwhichlocreg(:)
  !!write(*,'(a,i3,3x,100i4)') 'iproc, owap', iproc, onWhichAtomp(:)

  !call initLocregs2(iproc, at%nat, rxyz, lzdig, input, Glr, locrad)
  call initLocregs2(iproc, at%nat, rxyz, lzdig, input, Glr, lin%locrad)
  allocate(lchi(lzdig%orbs%npsidim+ndebug),stat=i_stat)
  call memocc(i_stat,chi,'lchi',subname)
  allocate(lhchi(lzdig%orbs%npsidim,at%nat),stat=i_stat)
  call memocc(i_stat,lhchi,'lhchi',subname)
  lchi=0.d0
  lhchi=0.d0


  hxh=.5_gp*input%hx
  hyh=.5_gp*input%hy
  hzh=.5_gp*input%hz


  ! Allocate the atomic orbitals chi and hchi (i.e. the Hamiltonian applied to them).
  allocate(chi(orbsig%npsidim+ndebug),stat=i_stat)
  call memocc(i_stat,chi,'chi',subname)
  allocate(hchi(orbsig%npsidim,at%nat),stat=i_stat)
  call memocc(i_stat,hchi,'hchi',subname)
  chi=0.d0
  hchi=0.d0

  !allocate arrays for the GPU if a card is present
  switchGPUconv=.false.
  switchOCLconv=.false.
  if (GPUconv .and. potshortcut ==0 ) then
     call prepare_gpu_for_locham(Glr%d%n1,Glr%d%n2,Glr%d%n3,nspin_ig,&
          input%hx,input%hy,input%hz,Glr%wfd,orbsig,GPU)
  else if (OCLconv .and. potshortcut ==0) then
     call allocate_data_OCL(Glr%d%n1,Glr%d%n2,Glr%d%n3,at%geocode,&
          nspin_ig,input%hx,input%hy,input%hz,Glr%wfd,orbsig,GPU)
     if (iproc == 0) write(*,*)&
          'GPU data allocated'
  else if (GPUconv .and. potshortcut >0 ) then
     switchGPUconv=.true.
     GPUconv=.false.
  else if (OCLconv .and. potshortcut >0 ) then
     switchOCLconv=.true.
     OCLconv=.false.
  end if


  ! Transform the Gaussian based orbitals to wavelets.
  !call gaussians_to_wavelets_new(iproc,nproc,Glr,orbsig,input%hx,input%hy,input%hz,G,&
  !     psigau(1,1,min(orbsig%isorb+1,orbsig%norb)),chi)
  call gaussians_to_wavelets_new(iproc,nproc,Glr,lzdig%orbs,input%hx,input%hy,input%hz,G,&
       psigau(1,1,min(lzdig%orbs%isorb+1,lzdig%orbs%norb)),chi)


  i_all=-product(shape(locrad))*kind(locrad)
  deallocate(locrad,stat=i_stat)
  call memocc(i_stat,i_all,'locrad',subname)

  
  ! Create the potential.
  call sumrho(iproc,nproc,orbsig,Glr,input%ixc,hxh,hyh,hzh,chi,rhopot,&
       & Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,1),nscatterarr,input%nspin,GPU, &
       & at%symObj, irrzon, phnons)
     
  !-- if spectra calculation uses a energy dependent potential
  !    input_wf_diag will write (to be used in abscalc)
  !    the density to the file electronic_density.cube
  !  The writing is activated if  5th bit of  in%potshortcut is on.
  if( iand( potshortcut,16)==0 .and. potshortcut /= 0) then
     call plot_density_cube_old(at%geocode,'electronic_density',&
          iproc,nproc,Glr%d%n1,Glr%d%n2,Glr%d%n3,Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,nscatterarr(iproc,2),  & 
          input%nspin,hxh,hyh,hzh,at,rxyz,ngatherarr,rhopot(1+nscatterarr(iproc,4)*Glr%d%n1i*Glr%d%n2i))
  endif
  !---
  
  if(lin%orbs%nspinor==4) then
     !this wrapper can be inserted inside the poisson solver 
     call PSolverNC(at%geocode,'D',iproc,nproc,Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,&
          nscatterarr(iproc,1),& !this is n3d
          input%ixc,hxh,hyh,hzh,&
          rhopot,pkernel,pot_ion,ehart,eexcu,vexcu,0.d0,.true.,4)
  else
     !Allocate XC potential
     if (nscatterarr(iproc,2) >0) then
        allocate(potxc(Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2)*input%nspin+ndebug),stat=i_stat)
        call memocc(i_stat,potxc,'potxc',subname)
     else
        allocate(potxc(1+ndebug),stat=i_stat)
        call memocc(i_stat,potxc,'potxc',subname)
     end if

     call XC_potential(at%geocode,'D',iproc,nproc,&
          Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,input%ixc,hxh,hyh,hzh,&
          rhopot,eexcu,vexcu,input%nspin,rhocore,potxc)


     if( iand(potshortcut,4)==0) then
        call H_potential(at%geocode,'D',iproc,nproc,&
             Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,hxh,hyh,hzh,&
             rhopot,pkernel,pot_ion,ehart,0.0_dp,.true.)
     endif


     !sum the two potentials in rhopot array
     !fill the other part, for spin, polarised
     if (input%nspin == 2) then
        call dcopy(Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2),rhopot(1),1,&
             rhopot(Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2)+1),1)
     end if
     !spin up and down together with the XC part
     call axpy(Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2)*input%nspin,1.0_dp,potxc(1),1,&
          rhopot(1),1)


     i_all=-product(shape(potxc))*kind(potxc)
     deallocate(potxc,stat=i_stat)
     call memocc(i_stat,i_all,'potxc',subname)

  end if


  if(potshortcut>0) then
!!$    if (GPUconv) then
!!$       call free_gpu(GPU,orbs%norbp)
!!$    end if
     if (switchGPUconv) then
        GPUconv=.true.
     end if
     if (switchOCLconv) then
        OCLconv=.true.
     end if

     call deallocate_orbs(orbsig,subname)
     
     !deallocate the gaussian basis descriptors
     call deallocate_gwf(G,subname)
    
     i_all=-product(shape(psigau))*kind(psigau)
     deallocate(psigau,stat=i_stat)
     call memocc(i_stat,i_all,'psigau',subname)
     call deallocate_comms(commsig,subname)
     i_all=-product(shape(norbsc_arr))*kind(norbsc_arr)
     deallocate(norbsc_arr,stat=i_stat)
     call memocc(i_stat,i_all,'norbsc_arr',subname)
    return 
  end if


  !call dcopy(orbsig%npsidim,psi,1,hpsi,1)
  if (input%exctxpar == 'OP2P') eexctX = -99.0_gp


  ! Copy Glr to lzd
  lzdig%Glr = Glr
  
  ! Copy nlpspd to lin%lzd
  lzdig%Gnlpspd = nlpspd
  
  ! Set localnorb
  do ilr=1,lzdig%nlr
      lzdig%Llr(ilr)%localnorb=0
      do iorb=1,lzdig%orbs%norbp
          if(onWhichAtomp(iorb)==ilr) then
              lzdig%Llr(ilr)%localnorb = lzdig%Llr(ilr)%localnorb+1
          end if
      end do
  end do


  ! Initialize the parameters for the communications of the potential.
  call initializeCommunicationPotential(iproc, nproc, nscatterarr, lzdig%orbs, lzdig, comgp, onWhichAtom)

  ! Post the messages for the communication of the potential.
  ndimpot = lin%lzd%Glr%d%n1i*lin%lzd%Glr%d%n2i*nscatterarr(iproc,2)
  call postCommunicationsPotential(iproc, nproc, ndimpot, rhopot, comgp)


  ! Orthogonalize the atomic basis functions (Loewdin).
  call orthonormalizeAtomicOrbitals(iproc, nproc, orbsig, commsig, Glr, chi)
  
  ! Gather the potential
  call gatherPotential(iproc, nproc, comgp)

  ! Build the potential.
  ! This is not required at the moment since HamiltonianApplicationConfinement has a different (old) structure.
  call full_local_potential(iproc,nproc,Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2),Glr%d%n1i*Glr%d%n2i*Glr%d%n3i,input%nspin,&
       orbsig%norb,orbsig%norbp,ngatherarr,rhopot,pot)

  ! Apply the Hamiltonian for each atom.
  ! onWhichAtomTemp indicates indicating that all orbitals feel the potential from atom iat.
  allocate(onWhichAtomTemp(lzdig%orbs%norbp), stat=istat)
  call memocc(i_stat,onWhichAtomTemp,'onWhichAtomTemp',subname)
  if(iproc==0) write(*,'(x,a)') 'Hamiltonian application for all atoms. This may take some time.'

  ! Transform chi to the localization region. This is not needed if we really habe O(N).
  ind1=1
  ind2=1
  do iorb=1,lzdig%orbs%norbp
      ilr = onWhichAtomp(iorb)
      ldim=lzdig%Llr(ilr)%wfd%nvctr_c+7*lzdig%Llr(ilr)%wfd%nvctr_f
      gdim=Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
      call psi_to_locreg2(iproc, nproc, ldim, gdim, lin%Llr(ilr), Glr, chi(ind1), lchi(ind2))
      ind1=ind1+Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
      ind2=ind2+lzdig%Llr(ilr)%wfd%nvctr_c+7*lzdig%Llr(ilr)%wfd%nvctr_f
  end do

  hchi=0.d0
  call cpu_time(t1)
  withConfinement=.true.
  do iat=1,at%nat
      do iorb=1,orbsig%norbp
          onWhichAtomTemp(iorb)=iat
      end do
      if(iproc==0) write(*,'(3x,a,i0,a)', advance='no') 'Hamiltonian application for atom ', iat, '... '
      !!call HamiltonianApplicationConfinement(iproc, nproc, at, orbsig, lin, input%hx, input%hy, input%hz, rxyz,&
      !!     nlpspd, proj, Glr, ngatherarr, Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2), &
      !!     rhopot(1), &
      !!     chi(1), hchi(1,iat), ekin_sum, epot_sum, eexctX, eproj_sum, input%nspin, GPU, rxyz, onWhichAtomTemp, pkernel=pkernelseq)
      call HamiltonianApplicationConfinement2(input, iproc, nproc, at, lzdig, lin, input%hx, input%hy, input%hz, rxyz,&
           proj, ngatherarr, comgp%nrecvBuf, comgp%recvBuf, lchi, lhchi(1,iat), &
           ekin_sum, epot_sum, eexctX, eproj_sum, input%nspin, GPU, radii_cf, comgp, onWhichAtomTemp, withConfinement, pkernel=pkernelseq)
      ind1=1
      ind2=1
      do iorb=1,lzdig%orbs%norbp
          ilr = onWhichAtomp(iorb)
          ldim=lzdig%Llr(ilr)%wfd%nvctr_c+7*lzdig%Llr(ilr)%wfd%nvctr_f
          gdim=Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
          call Lpsi_to_global2(iproc, nproc, ldim, gdim, lzdig%orbs%norb, lin%orbs%nspinor, input%nspin, Glr, lzdig%Llr(ilr), lhchi(ind2,iat), hchi(ind1,iat))
          ind1=ind1+Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
          ind2=ind2+lzdig%Llr(ilr)%wfd%nvctr_c+7*lzdig%Llr(ilr)%wfd%nvctr_f
      end do

      if(iproc==0) write(*,'(a)') 'done.'
  end do
  call cpu_time(t2)
  time=t2-t1
  call mpiallred(time, 1, mpi_sum, mpi_comm_world, ierr)
  if(iproc==0) write(*,'(x,a,es10.3)') 'time for applying potential:', time/dble(nproc)
  

  ! Deallocate potential.
  call free_full_potential(nproc,pot,subname)


  !free GPU if it is the case
  if (GPUconv) then
     call free_gpu(GPU,orbsig%norbp)
  else if (OCLconv) then
     call free_gpu_OCL(GPU,orbsig,nspin_ig)
  end if


   ! onWhichAtomPhi indicates on which atom the basis functions phi are centered.
   allocate(onWhichAtomPhi(lin%orbs%norb))
   call memocc(i_stat,onWhichAtomPhi,'onWhichAtomPhi',subname)
   onWhichAtomPhi=0
   do iorb=1,lin%orbs%norbp
       onWhichAtomPhi(lin%orbs%isorb+iorb)=lin%onWhichAtom(iorb)
   end do
   call mpiallred(onWhichAtomPhi(1), lin%orbs%norb, mpi_sum, mpi_comm_world, iorb)

   ! Build a linear combination of atomic orbitals to get a goor input guess for phi.
   !! THIS IS THE ORIGINAL
   call cpu_time(t1)
   !call buildLinearCombinations(iproc, nproc, orbsig, lin%orbs, commsig, lin%comms, at, Glr, lin%norbsPerType, &
   !     onWhichAtom, chi, hchi, phi, rxyz, onWhichAtomPhi, lin)
   call buildLinearCombinations2(iproc, nproc, orbsig, lin%orbs, commsig, lin%comms, at, Glr, lin%norbsPerType, &
        onWhichAtom, chi, hchi, phi, rxyz, onWhichAtomPhi, lin)
   call cpu_time(t2)
   time=t2-t1
   call mpiallred(time, 1, mpi_sum, mpi_comm_world, ierr)
   if(iproc==0) write(*,'(x,a,es10.3)') 'time for "buildLinearCombinations":', time/dble(nproc)

  if(iproc==0) write(*,'(x,a)') '------------------------------------------------------------- Input guess generated.'



  ! Deallocate all local arrays.
  call deallocate_gwf(G,subname)

  i_all=-product(shape(psigau))*kind(psigau)
  deallocate(psigau,stat=i_stat)
  call memocc(i_stat,i_all,'psigau',subname)

  i_all=-product(shape(norbsc_arr))*kind(norbsc_arr)
  deallocate(norbsc_arr,stat=i_stat)
  call memocc(i_stat,i_all,'norbsc_arr',subname)

  call deallocate_comms(commsig,subname)

  call deallocate_orbs(orbsig,subname)

  i_all=-product(shape(onWhichAtom))*kind(onWhichAtom)
  deallocate(onWhichAtom, stat=i_stat)
  call memocc(i_stat, i_all, 'onWhichAtom',subname)

  i_all=-product(shape(onWhichAtomPhi))*kind(onWhichAtomPhi)
  deallocate(onWhichAtomPhi, stat=i_stat)
  call memocc(i_stat, i_all, 'onWhichAtomPhi',subname)

  !!i_all=-product(shape(onWhichAtomp))*kind(onWhichAtomp)
  !!deallocate(onWhichAtomp, stat=i_stat)
  !!call memocc(i_stat, i_all, 'onWhichAtomp',subname)

  i_all=-product(shape(onWhichAtomTemp))*kind(onWhichAtomTemp)
  deallocate(onWhichAtomTemp, stat=i_stat)
  call memocc(i_stat, i_all, 'onWhichAtomTemp',subname)
  
  i_all=-product(shape(norbsPerAt))*kind(norbsPerAt)
  deallocate(norbsPerAt, stat=i_stat)
  call memocc(i_stat, i_all, 'norbsPerAt',subname)

  i_all=-product(shape(chi))*kind(chi)
  deallocate(chi, stat=i_stat)
  call memocc(i_stat, i_all, 'chi',subname)

  i_all=-product(shape(hchi))*kind(hchi)
  deallocate(hchi, stat=i_stat)
  call memocc(i_stat, i_all, 'hchi',subname)


END SUBROUTINE inputguessConfinement





subroutine buildLinearCombinations2(iproc, nproc, orbsig, orbs, commsig, comms, at, Glr, norbsPerType, &
           onWhichAtom, chi, hchi, phi, rxyz, onWhichAtomPhi, lin)
!
use module_base
use module_types
use module_interfaces
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data),intent(in):: orbsig, orbs
type(communications_arrays),intent(in):: commsig, comms
type(atoms_data),intent(in):: at
type(locreg_descriptors),intent(in):: Glr
integer,dimension(at%ntypes):: norbsPerType
integer,dimension(orbsig%norb),intent(in):: onWhichAtom
real(8),dimension(orbsig%npsidim):: chi
real(8),dimension(orbsig%npsidim,at%nat):: hchi
real(8),dimension(orbs%npsidim):: phi
real(8),dimension(3,at%nat):: rxyz
integer,dimension(orbs%norb):: onWhichAtomPhi
type(linearParameters),intent(in):: lin

! Local variables
integer:: iorb, jorb, korb, iat, ist, jst, nvctrp, iall, istat, ierr, infoCoeff, k, l,it, iiAt, jjAt
real(8),dimension(:),allocatable:: alpha, coeffPad, coeff2, gradTemp, gradOld, fnrmArr, fnrmOvrlpArr, fnrmOldArr, grad
real(8),dimension(:,:),allocatable:: ovrlp, ovrlpTemp
real(8),dimension(:,:),allocatable:: coeff, lagMat, coeffOld
real(8),dimension(:,:,:),allocatable:: Ham, HamPad
real(8),dimension(:),pointer:: chiw
integer,dimension(:),allocatable:: recvcounts, displs, norb_par
real(8):: ddot, cosangle, tt, dnrm2, fnrm, meanAlpha, cut, trace, traceOld, fnrmMax
logical:: converged
character(len=*),parameter:: subname='buildLinearCombinations'
real(4):: ttreal
integer:: wholeGroup, newGroup, newComm, norbtot
integer,dimension(:),allocatable:: newID
  
! new
real(8),dimension(:),allocatable:: work, eval, evals
real(8),dimension(:,:),allocatable:: tempMat
integer:: lwork, ii, info, iiAtprev, i, jproc, norbTarget, sendcount
type(inguessParameters):: ip

  if(iproc==0) write(*,'(x,a)') '------------------------------- Minimizing trace in the basis of the atomic orbitals'

  ! Allocate the local arrays that are hold by all processes.
  allocate(Ham(orbsig%norb,orbsig%norb,at%nat), stat=istat)
  call memocc(istat, Ham, 'Ham', subname)
  allocate(coeff(orbsig%norb,orbs%norb), stat=istat)
  call memocc(istat, coeff, 'coeff', subname)

  ! Transpose all the wavefunctions for having a piece of all the orbitals 
  ! for each processor
  allocate(chiw(orbsig%npsidim+ndebug),stat=istat)
  call memocc(istat,chiw, 'chiw', subname)
  call transpose_v(iproc, nproc, orbsig, Glr%wfd, commsig, chi(1), work=chiw)
  do iat=1,at%nat
      call transpose_v(iproc, nproc, orbsig, Glr%wfd, commsig, hchi(1,iat), work=chiw)
  end do
  iall=-product(shape(chiw))*kind(chiw)
  deallocate(chiw,stat=istat)
  call memocc(istat, iall, 'chiw', subname)

  ! Calculate Hamiltonian matrix.
  Ham=0.d0
  nvctrp=commsig%nvctr_par(iproc,1) ! 1 for k-points
  do iat=1,at%nat
      ist=1
      do iorb=1,orbsig%norb
          jst=1
          do jorb=1,orbsig%norb
              Ham(jorb,iorb,iat)=ddot(nvctrp, chi(ist), 1, hchi(jst,iat), 1)
              jst=jst+nvctrp
          end do
          ist=ist+nvctrp
      end do
  end do
  call mpiallred(Ham(1,1,1), orbsig%norb**2*at%nat, mpi_sum, mpi_comm_world, ierr)



  ! Determine the number of processes we need, which will be stored in ip%nproc.
  ! This is the only variable in ip that all processes (also those which do not
  ! participate in the minimization of the trace) will know. The other variables
  ! are known only by the active processes and will be determined by initializeInguessParameters.
  if(lin%norbsPerProcIG>orbs%norb) then
      norbTarget=orbs%norb
  else
     norbTarget=lin%norbsperProcIG
  end if
  ip%nproc=ceiling(dble(orbs%norb)/dble(norbTarget))
  ip%nproc=min(ip%nproc,nproc)
  if(iproc==0) write(*,'(a,i0,a)') 'The minimization is performed using ', ip%nproc, ' processes.'

  ! Create the new communicator newComm.
  allocate(newID(0:ip%nproc-1), stat=istat)
  call memocc(istat, newID, 'newID', subname)
  do jproc=0,ip%nproc-1
     newID(jproc)=jproc
  end do
  call mpi_comm_group(mpi_comm_world, wholeGroup, ierr)
  call mpi_group_incl(wholeGroup, ip%nproc, newID, newGroup, ierr)
  call mpi_comm_create(mpi_comm_world, newGroup, newComm, ierr)

  ! Everything inside this if statements is only executed by the processes in newComm.
  processIf: if(iproc<ip%nproc) then

      ! Initialize the parameters for performing tha calculations in parallel.
      call initializeInguessParameters(iproc, orbs, orbsig, newComm, ip)
    
      ! Allocate the local arrays.
      call allocateArrays()
    
      ! Initialize the coefficient vectors. Put random number to palces where it is
      ! reasonable (i.e. close to the atom where the basis function is centered).
      call initRandomSeed(iproc, 1)
      do ii=1,orbsig%isorb*ip%norbtotPad
          call random_number(ttreal)
      end do
      coeffPad=0.d0
    
      ii=0
      do iorb=1,ip%norb_par(iproc)
          iiAt=onWhichAtom(ip%isorb+iorb)
          ! Do not fill up to the boundary of the localization region, but only up to one fourth of it.
          cut=0.0625d0*lin%locrad(at%iatype(iiAt))**2
          do jorb=1,ip%norbtot
              jjAt=onWhichAtom(jorb)
              tt = (rxyz(1,iiat)-rxyz(1,jjAt))**2 + (rxyz(2,iiat)-rxyz(2,jjAt))**2 + (rxyz(3,iiat)-rxyz(3,jjAt))**2
              if(tt>cut) then
                   coeffPad((iorb-1)*ip%norbtotPad+jorb)=0.d0
              else
                  call random_number(ttreal)
                  coeffPad((iorb-1)*ip%norbtotPad+jorb)=dble(ttreal)
              end if
          end do
      end do
    
    
      ! Pad the Hamiltonian with zeros.
      do iat=1,at%nat
          do iorb=1,ip%norbtot
              call dcopy(ip%norbtot, Ham(1,iorb,iat), 1, HamPad(1,iorb,iat), 1)
              do i=ip%norbtot+1,ip%norbtotPad
                  HamPad(i,iorb,iat)=0.d0
              end do
          end do
          do iorb=ip%norbtot+1,ip%norbtotPad
              do i=1,ip%norbtotPad
                  HamPad(i,iorb,iat)=0.d0
              end do
          end do
      end do
    
      
      ! Initial step size for the optimization
      alpha=1.d-3
    
      ! Flag which checks convergence.
      converged=.false.
    
      if(iproc==0) write(*,'(x,a)') '============================== optmizing coefficients =============================='
    
      ! The optimization loop.
    
      ! Transpose the coefficients for the first orthonormalization.
      call transposeInguess(iproc, ip, newComm, coeffPad)
    
      iterLoop: do it=1,lin%nItInguess
    
          if (iproc==0 .and. mod(it,1)==0) then
              write( *,'(1x,a,i0)') repeat('-',77 - int(log(real(it))/log(10.))) // ' iter=', it
          endif
    
    
          ! Othonormalize the coefficients.
          !call orthonormalizeCoefficients(orbs, orbsig, coeff)
          call orthonormalizeCoefficients_parallel(iproc, ip, newComm, coeffPad)
          call untransposeInguess(iproc, ip, newComm, coeffPad)
    
    
          ! Calculate the gradient grad. At the same time we determine whether the step size shall be increased
          ! or decreased (depending on gradient feedback).
          meanAlpha=0.d0
          grad=0.d0
          do iorb=1,ip%norb_par(iproc)
              iiAt=onWhichAtom(ip%isorb+iorb)
              call dgemv('n', ip%norbtotPad, ip%norbtotPad, 1.d0, HamPad(1,1,iiAt), ip%norbtotPad, &
                   coeffPad((iorb-1)*ip%norbtotPad+1), 1, 0.d0, grad((iorb-1)*ip%norbtotPad+1), 1)
          end do
      
          ! Apply the orthoconstraint to the gradient. To do so first calculate the Lagrange
          ! multiplier matrix.
          lagMat=0.d0
          if(it>1) then
              traceOld=trace
          else
              traceOld=1.d10
          end if
          trace=0.d0
          call transposeInguess(iproc, ip, newComm, coeffPad)
          call transposeInguess(iproc, ip, newComm, grad)
          call gemm('t', 'n', ip%norb, ip%norb, ip%nvctrp, 1.0d0, coeffPad(1), &
               ip%nvctrp, grad(1), ip%nvctrp, 0.d0, lagMat(1,1), ip%norb)
          call mpiallred(lagMat(1,1), ip%norb**2, mpi_sum, newComm, ierr)
          do iorb=1,orbs%norb
              trace=trace+lagMat(iorb,iorb)
          end do
    
    
          ! Now apply the orthoconstraint.
          call dgemm('n', 'n', ip%nvctrp, ip%norb, ip%norb, -.5d0, coeffPad(1), ip%nvctrp, &
               lagMat(1,1), ip%norb, 1.d0, grad(1), ip%nvctrp)
          call dgemm('n', 't', ip%nvctrp, ip%norb, ip%norb, -.5d0, coeffPad(1), ip%nvctrp, &
               lagMat(1,1), ip%norb, 1.d0, grad(1), ip%nvctrp)
    
    
          ! Calculate the gradient norm.
          fnrm=0.d0
          do iorb=1,ip%norb
              fnrmArr(iorb)=ddot(ip%nvctrp, grad((iorb-1)*ip%nvctrp+1), 1, grad((iorb-1)*ip%nvctrp+1), 1)
              if(it>1) fnrmOvrlpArr(iorb)=ddot(ip%nvctrp, grad((iorb-1)*ip%nvctrp+1), 1, gradOld((iorb-1)*ip%nvctrp+1), 1)
          end do
          call mpiallred(fnrmArr(1), ip%norb, mpi_sum,newComm, ierr)
          call mpiallred(fnrmOvrlpArr(1), ip%norb, mpi_sum, newComm, ierr)
          call dcopy(ip%nvctrp*ip%norb, grad(1), 1, gradOld(1), 1)
    
          ! Keep the gradient for the next iteration.
          if(it>1) then
              call dcopy(ip%norb, fnrmArr(1), 1, fnrmOldArr(1), 1)
          end if
    
          fnrmMax=0.d0
          do iorb=1,ip%norb
              fnrm=fnrm+fnrmArr(iorb)
              if(fnrmArr(iorb)>fnrmMax) fnrmMax=fnrmArr(iorb)
              if(it>1) then
              ! Adapt step size for the steepest descent minimization.
                  tt=fnrmOvrlpArr(iorb)/sqrt(fnrmArr(iorb)*fnrmOldArr(iorb))
                  if(tt>.9d0) then
                      alpha(iorb)=alpha(iorb)*1.05d0
                  else
                      alpha(iorb)=alpha(iorb)*.5d0
                  end if
              end if
          end do
    
         fnrm=sqrt(fnrm)
         fnrmMax=sqrt(fnrmMax)
    
         ! Determine the mean step size for steepest descent iterations.
         tt=sum(alpha)
         meanAlpha=tt/dble(ip%norb)
    
          ! Precondition the gradient.
          !!if(fnrm<1.d0) call preconditionGradient(iproc, nproc, orbsig, orbs, at, Ham, lagMat, onWhichAtomPhi, grad, it, evals)
      
    
          ! Write some informations to the screen, but only every 1000th iteration.
          if(iproc==0 .and. mod(it,1)==0) write(*,'(x,a,es11.2,es22.13,es10.2)') 'fnrm, trace, mean alpha', &
              fnrm, trace, meanAlpha
          
          ! Check for convergence.
          if(fnrm<1.d-3) then
              if(iproc==0) write(*,'(x,a,i0,a)') 'converged in ', it, ' iterations.'
              if(iproc==0) write(*,'(3x,a,2es14.5)') 'Final values for fnrm, trace:', fnrm, trace
              converged=.true.
              infoCoeff=it
              exit
          end if
      
          ! Quit if the maximal number of iterations is reached.
          if(it==lin%nItInguess) then
              if(iproc==0) write(*,'(x,a,i0,a)') 'WARNING: not converged within ', it, &
                  ' iterations! Exiting loop due to limitations of iterations.'
              if(iproc==0) write(*,'(x,a,2es15.7,f12.7)') 'Final values for fnrm, trace: ', fnrm, trace
              infoCoeff=-1
              exit
          end if
    
          ! Improve the coefficients (by steepet descent).
          do iorb=1,orbs%norb
              call daxpy(ip%nvctrp, -alpha(iorb), grad((iorb-1)*ip%nvctrp+1), 1, coeffPad((iorb-1)*ip%nvctrp+1), 1)
          end do
    
    
      end do iterLoop
    
    
      if(iproc==0) write(*,'(x,a)') '===================================================================================='
    
    
      ! Cut out the zeros
      call untransposeInguess(iproc, ip, newComm, coeffPad)
      allocate(coeff2(ip%norbtot*ip%norb_par(iproc)), stat=istat)
      call memocc(istat, coeff2, 'coeff2', subname)
      do iorb=1,ip%norb_par(iproc)
          call dcopy(ip%norbtot, coeffPad((iorb-1)*ip%norbtotPad+1), 1, coeff2((iorb-1)*ip%norbtot+1), 1)
      end do


  end if processIf
  
  call mpi_barrier(mpi_comm_world, ierr)
  
  ! Allocate coeff2 for those processes which did not allocate it
  ! during the previous if statement.
  if(iproc>=ip%nproc) then
      allocate(coeff2(1), stat=istat)
      call memocc(istat, coeff2, 'coeff2', subname)
  end if
  
  
  ! Now collect all coefficients on all processes.
  allocate(recvcounts(0:nproc-1), stat=istat)
  call memocc(istat, recvcounts, 'recvcounts', subname)
  allocate(displs(0:nproc-1), stat=istat)
  call memocc(istat, displs, 'displs', subname)
  
  ! Send ip%norb_par and ip%norbtot to all processes.
  allocate(norb_par(0:ip%nproc-1), stat=istat)
  call memocc(istat, norb_par, 'norb_par', subname)
  if(iproc==0) then
      do i=0,ip%nproc-1
          norb_par(i)=ip%norb_par(i)
      end do
      norbtot=ip%norbtot
  end if
  call mpi_bcast(norb_par(0), ip%nproc, mpi_integer, 0, mpi_comm_world, ierr)
  call mpi_bcast(norbtot, 1, mpi_integer, 0, mpi_comm_world, ierr)
  
  ! Define the parameters, for the mpi_allgatherv.
  ii=0
  do jproc=0,ip%nproc-1
      recvcounts(jproc)=norbtot*norb_par(jproc)
      displs(jproc)=ii
      ii=ii+recvcounts(jproc)
  end do
  do jproc=ip%nproc,nproc-1
      recvcounts(jproc)=0
      displs(jproc)=ii
      ii=ii+recvcounts(jproc)
  end do
  if(iproc<ip%nproc) then
      sendcount=ip%norbtot*ip%norb_par(iproc)
  else
      sendcount=0
  end if

  ! Gather the coefficients.
  call mpi_allgatherv(coeff2(1), sendcount, mpi_double_precision, coeff(1,1), recvcounts, &
       displs, mpi_double_precision, mpi_comm_world, ierr)


  ! Now every process has all coefficients, so we can build the linear combinations.
  phi=0.d0
  ist=1
  do iorb=1,orbs%norb
      jst=1
      do jorb=1,orbsig%norb
          call daxpy(nvctrp, coeff(jorb,iorb), chi(jst), 1, phi(ist), 1)
          jst=jst+nvctrp
      end do
      ist=ist+nvctrp
  end do

  

  ! Untranpose the orbitals.
  allocate(chiw(orbs%npsidim+ndebug),stat=istat)
  call memocc(istat, chiw, 'chiw', subname)
  call untranspose_v(iproc, nproc, orbs, Glr%wfd, comms, phi, work=chiw)
  iall=-product(shape(chiw))*kind(chiw)
  deallocate(chiw, stat=istat)
  call memocc(istat, iall, 'chiw', subname)


  ! Deallocate the local arrays.
  iall=-product(shape(Ham))*kind(Ham)
  deallocate(Ham, stat=istat)
  call memocc(istat, iall, 'Ham', subname)

  iall=-product(shape(coeff))*kind(coeff)
  deallocate(coeff, stat=istat)
  call memocc(istat, iall, 'coeff', subname)


  
  contains

    subroutine allocateArrays()
      allocate(coeffPad(max(ip%norbtotPad*ip%norb_par(iproc), ip%nvctrp*ip%norb)), stat=istat)
      call memocc(istat, coeffPad, 'coeffPad', subname)
      allocate(HamPad(ip%norbtotPad,ip%norbtotPad,at%nat), stat=istat)
      call memocc(istat, HamPad, 'HamPad', subname)
      allocate(grad(max(ip%norbtotPad*ip%norb_par(iproc), ip%nvctrp*ip%norb)), stat=istat)
      call memocc(istat, grad, 'grad', subname)
      allocate(gradOld(max(ip%norbtotPad*ip%norb_par(iproc), ip%nvctrp*ip%norb)), stat=istat)
      call memocc(istat, gradOld, 'gradOld', subname)
      allocate(fnrmArr(ip%norb), stat=istat)
      call memocc(istat, fnrmArr, 'fnrmArr', subname)
      allocate(fnrmOvrlpArr(ip%norb), stat=istat)
      call memocc(istat, fnrmOvrlpArr, 'fnrmOvrlpArr', subname)
      allocate(fnrmOldArr(ip%norb), stat=istat)
      call memocc(istat, fnrmOldArr, 'fnrmOldArr', subname)
      allocate(alpha(orbs%norb), stat=istat)
      call memocc(istat, alpha, 'alpha', subname)
      allocate(lagMat(orbs%norb,orbs%norb), stat=istat)
      call memocc(istat, lagMat, 'lagMat', subname)
    end subroutine allocateArrays


    subroutine deallocateArrays()
      iall=-product(shape(grad))*kind(grad)
      deallocate(grad, stat=istat)
      call memocc(istat, iall, 'grad', subname)

      iall=-product(shape(gradOld))*kind(gradOld)
      deallocate(gradOld, stat=istat)
      call memocc(istat, iall, 'gradOld', subname)

      iall=-product(shape(alpha))*kind(alpha)
      deallocate(alpha, stat=istat)
      call memocc(istat, iall, 'alpha', subname)

      iall=-product(shape(lagMat))*kind(lagMat)
      deallocate(lagMat, stat=istat)
      call memocc(istat, iall, 'lagMat', subname)
     
      iall=-product(shape(coeffPad))*kind(coeffPad)
      deallocate(coeffPad, stat=istat)
      call memocc(istat, iall, 'coeffPad', subname)

      iall=-product(shape(HamPad))*kind(HamPad)
      deallocate(HamPad, stat=istat)
      call memocc(istat, iall, 'HamPad', subname)

      iall=-product(shape(fnrmArr))*kind(fnrmArr)
      deallocate(fnrmArr, stat=istat)
      call memocc(istat, iall, 'fnrmArr', subname)

      iall=-product(shape(fnrmOvrlpArr))*kind(fnrmOvrlpArr)
      deallocate(fnrmOvrlpArr, stat=istat)
      call memocc(istat, iall, 'fnrmOvrlpArr', subname)

      iall=-product(shape(fnrmOldArr))*kind(fnrmOldArr)
      deallocate(fnrmOldArr, stat=istat)
      call memocc(istat, iall, 'fnrmOldArr', subname)


    end subroutine deallocateArrays


end subroutine buildLinearCombinations2



!!!subroutine buildLinearCombinations(iproc, nproc, orbsig, orbs, commsig, comms, at, Glr, norbsPerType, &
!!!           onWhichAtom, chi, hchi, phi, rxyz, onWhichAtomPhi, lin)
!!!!
!!!use module_base
!!!use module_types
!!!use module_interfaces
!!!implicit none
!!!
!!!! Calling arguments
!!!integer,intent(in):: iproc, nproc
!!!type(orbitals_data),intent(in):: orbsig, orbs
!!!type(communications_arrays),intent(in):: commsig, comms
!!!type(atoms_data),intent(in):: at
!!!type(locreg_descriptors),intent(in):: Glr
!!!integer,dimension(at%ntypes):: norbsPerType
!!!integer,dimension(orbsig%norb),intent(in):: onWhichAtom
!!!real(8),dimension(orbsig%npsidim):: chi
!!!real(8),dimension(orbsig%npsidim,at%nat):: hchi
!!!real(8),dimension(orbs%npsidim):: phi
!!!real(8),dimension(3,at%nat):: rxyz
!!!integer,dimension(orbs%norb):: onWhichAtomPhi
!!!type(linearParameters),intent(in):: lin
!!!
!!!! Local variables
!!!integer:: iorb, jorb, korb, iat, ist, jst, nvctrp, iall, istat, ierr, infoCoeff, k, l,it, iiAt, jjAt
!!!real(8),dimension(:),allocatable:: alpha, tmparr
!!!real(8),dimension(:,:),allocatable:: ovrlp, ovrlpTemp, coeffTemp
!!!real(8),dimension(:,:),allocatable:: coeff, grad, gradOld, lagMat, ovrlpLarge, coeffOld
!!!real(8),dimension(:,:,:),allocatable:: Ham, tempArr
!!!real(8),dimension(:),pointer:: chiw
!!!real(8):: ddot, cosangle, tt, dnrm2, fnrm, meanAlpha, cut, trace, traceOld
!!!logical:: converged
!!!character(len=*),parameter:: subname='buildLinearCombinations'
!!!real(4):: ttreal
!!!  
!!!! new
!!!real(8),dimension(:),allocatable:: work, eval, evals
!!!real(8),dimension(:,:),allocatable:: tempMat
!!!integer:: lwork, ii, info, iiAtprev
!!!
!!!  if(iproc==0) write(*,'(x,a)') '------------------------------- Minimizing trace in the basis of the atomic orbitals'
!!!
!!!  ! Allocate the local arrays
!!!  allocate(coeff(orbsig%norb,orbs%norb), stat=istat)
!!!  call memocc(istat, coeff, 'coeff', subname)
!!!  allocate(alpha(orbs%norb), stat=istat)
!!!  call memocc(istat, alpha, 'alpha', subname)
!!!  allocate(grad(orbsig%norb,orbs%norb), stat=istat)
!!!  call memocc(istat, grad, 'grad', subname)
!!!  allocate(gradOld(orbsig%norb,orbs%norb), stat=istat)
!!!  call memocc(istat, gradOld, 'gradOld', subname)
!!!  allocate(lagMat(orbs%norb,orbs%norb), stat=istat)
!!!  call memocc(istat, lagMat, 'lagMat', subname)
!!!  allocate(Ham(orbsig%norb,orbsig%norb,at%nat), stat=istat)
!!!  call memocc(istat, Ham, 'Ham', subname)
!!!
!!!  ! Transpose all the wavefunctions for having a piece of all the orbitals 
!!!  ! for each processor
!!!  allocate(chiw(orbsig%npsidim+ndebug),stat=istat)
!!!  call memocc(istat,chiw, 'chiw', subname)
!!!  call transpose_v(iproc, nproc, orbsig, Glr%wfd, commsig, chi(1), work=chiw)
!!!  do iat=1,at%nat
!!!      call transpose_v(iproc, nproc, orbsig, Glr%wfd, commsig, hchi(1,iat), work=chiw)
!!!  end do
!!!  iall=-product(shape(chiw))*kind(chiw)
!!!  deallocate(chiw,stat=istat)
!!!  call memocc(istat, iall, 'chiw', subname)
!!!
!!!  ! Calculate Hamiltonian matrix.
!!!  Ham=0.d0
!!!  nvctrp=commsig%nvctr_par(iproc,1) ! 1 for k-points
!!!  do iat=1,at%nat
!!!      ist=1
!!!      do iorb=1,orbsig%norb
!!!          jst=1
!!!          do jorb=1,orbsig%norb
!!!              Ham(jorb,iorb,iat)=ddot(nvctrp, chi(ist), 1, hchi(jst,iat), 1)
!!!              jst=jst+nvctrp
!!!          end do
!!!          ist=ist+nvctrp
!!!      end do
!!!  end do
!!!  call mpiallred(Ham(1,1,1), orbsig%norb**2*at%nat, mpi_sum, mpi_comm_world, ierr)
!!!
!!!
!!!!!! Calculate proconditioning matrix
!!!!!ncplx=1
!!!!!cprecr=.5d0
!!!!!do iat=1,at%nat
!!!!!    do iorb=1,orbsig%norbp
!!!!!        call allocate_work_arrays(Glr%geocode, Glr%hybrid_on, ncplx, Glr%d, w)
!!!!!        call differentiateBetweenBoundaryConditions(ncplx, Glr, inut%hx, inout%hy, input%hz, kx, ky, kz, cprecr, chi, pchi, ,w, scal, rxyzParab, orbs, potentialPrefac, confPotOrder, it)
!!!!!        call deallocate_work_arrays(lr%geocode,lr%hybrid_on,ncplx,w)
!!!!!    end do
!!!!!end do
!!!
!!!
!!!
!!!
!!!  processIf: if(iproc==0) then
!!!
!!!      ! Initialize the coefficient vectors. Put random number to palces where it is
!!!      ! reasonable (i.e. close to the atom where the basis function is centered).
!!!allocate(tmparr(196))
!!!open(unit=123, file='coeffs')
!!!do iorb=1,196
!!!    read(123,*) tmparr(iorb)
!!!end do
!!!      ii=0
!!!      do iorb=1,orbs%norb
!!!          iiAt=onWhichAtomPhi(iorb)
!!!          !cut=1.d0/lin%potentialPrefac(at%iatype(iiAt))
!!!          !cut=cut**(1.d0/dble(lin%confPotOrder))
!!!          cut=lin%locrad(at%iatype(iiAt))**2
!!!          do jorb=1,orbsig%norb
!!!              jjAt=onWhichAtom(jorb)
!!!              tt = (rxyz(1,iiat)-rxyz(1,jjAt))**2 + (rxyz(2,iiat)-rxyz(2,jjAt))**2 + (rxyz(3,iiat)-rxyz(3,jjAt))**2
!!!              if(tt>cut) then
!!!                  coeff(jorb,iorb)=0.d0
!!!              else
!!!                  call random_number(ttreal)
!!!                  coeff(jorb,iorb)=dble(ttreal)
!!!                  ii=ii+1
!!!                  coeff(jorb,iorb)=tmparr(ii)
!!!              end if
!!!          end do
!!!      end do
!!!
!!!    !!! Diagonalize all Hamiltonian matrices
!!!    !!allocate(tempMat(orbsig%norb,orbsig%norb))
!!!    !!allocate(eval(orbsig%norb))
!!!    !!lwork=1000*orbsig%norb
!!!    !!allocate(work(lwork))
!!!    !!allocate(evals(orbs%norb))
!!!    !!iiAtprev=0
!!!    !!ii=0
!!!    !!do iorb=1,orbs%norb
!!!    !!    iiAt=onWhichAtomPhi(iorb)
!!!    !!    ii=ii+1
!!!    !!    if(iiAt>iiAtprev) then
!!!    !!        ii=1
!!!    !!        call dcopy(orbsig%norb**2, Ham(1,1,iiAt), 1, tempMat(1,1), 1)
!!!    !!        call dsyev('n', 'l', orbsig%norb, tempMat, orbsig%norb, eval, work, lwork, info)
!!!    !!    end if
!!!    !!    evals(iorb)=eval(ii)
!!!    !!    iiAtprev=iiAt
!!!    !!end do
!!!    !!deallocate(tempMat)
!!!    !!deallocate(eval)
!!!    !!deallocate(work)
!!!
!!!
!!!    
!!!    ! Initial step size for the optimization
!!!    alpha=5.d-1
!!!
!!!    ! Flag which checks convergence.
!!!    converged=.false.
!!!
!!!    if(iproc==0) write(*,'(x,a)') '============================== optmizing coefficients =============================='
!!!
!!!    ! The optimization loop.
!!!    iterLoop: do it=1,lin%nItInguess
!!!
!!!        if (iproc==0 .and. mod(it,1)==0) then
!!!            write( *,'(1x,a,i0)') repeat('-',77 - int(log(real(it))/log(10.))) // ' iter=', it
!!!        endif
!!!
!!!
!!!        ! Othonormalize the coefficients.
!!!        call orthonormalizeCoefficients(orbs, orbsig, coeff)
!!!
!!!
!!!
!!!        ! Calculate the gradient grad. At the same time we determine whether the step size shall be increased
!!!        ! or decreased (depending on gradient feedback).
!!!        meanAlpha=0.d0
!!!        grad=0.d0
!!!        do iorb=1,orbs%norb
!!!            iiAt=onWhichAtomPhi(iorb)
!!!            call dgemv('n', orbsig%norb, orbsig%norb, 1.d0, Ham(1,1,iiAt), orbsig%norb, coeff(1,iorb), 1, 0.d0, grad(1,iorb), 1)
!!!            if(it>1) then
!!!                cosangle=ddot(orbsig%norb, grad(1,iorb), 1, gradOld(1,iorb), 1)
!!!                cosangle=cosangle/dnrm2(orbsig%norb, grad(1,iorb), 1)
!!!                cosangle=cosangle/dnrm2(orbsig%norb, gradOld(1,iorb), 1)
!!!                if(cosangle>.8d0 .and. trace<traceOld+1.d-6*abs(traceOld)) then
!!!                    alpha(iorb)=max(alpha(iorb)*1.05d0,1.d-6)
!!!                else
!!!                    alpha(iorb)=max(alpha(iorb)*.5d0,1.d-6)
!!!                end if
!!!            end if
!!!            call dcopy(orbsig%norb, grad(1,iorb), 1, gradOld(1,iorb), 1)
!!!            meanAlpha=meanAlpha+alpha(iorb)
!!!        end do
!!!        meanAlpha=meanAlpha/orbs%norb
!!!do iall=1,orbsig%norb
!!!    write(23000+iproc,'(100f12.5)') (coeff(iall,iorb), iorb=1,orbs%norb)
!!!    write(24000+iproc,'(100f12.5)') (grad(iall,iorb), iorb=1,orbs%norb)
!!!end do
!!!write(23000+iproc,*) '----------------------'
!!!write(24000+iproc,*) '----------------------'
!!!    
!!!    
!!!        ! Apply the orthoconstraint to the gradient. To do so first calculate the Lagrange
!!!        ! multiplier matrix.
!!!        lagMat=0.d0
!!!        if(it>1) then
!!!            traceOld=trace
!!!        else
!!!            traceOld=1.d10
!!!        end if
!!!        trace=0.d0
!!!        call gemm('t', 'n', orbs%norb, orbs%norb, orbsig%norb, 1.0d0, coeff(1,1), &
!!!             orbsig%norb, grad(1,1), orbsig%norb, 0.d0, lagMat(1,1), orbs%norb)
!!!        do iorb=1,orbs%norb
!!!            trace=trace+lagMat(iorb,iorb)
!!!        end do
!!!
!!!
!!!        ! Now apply the orthoconstraint.
!!!        call dgemm('n', 'n', orbsig%norb, orbs%norb, orbs%norb, -.5d0, coeff(1,1), orbsig%norb, &
!!!             lagMat(1,1), orbs%norb, 1.d0, grad(1,1), orbsig%norb)
!!!        call dgemm('n', 't', orbsig%norb, orbs%norb, orbs%norb, -.5d0, coeff(1,1), orbsig%norb, &
!!!             lagMat(1,1), orbs%norb, 1.d0, grad(1,1), orbsig%norb)
!!!do iall=1,orbsig%norb
!!!    write(25000+iproc,'(100f12.5)') (grad(iall,iorb), iorb=1,orbs%norb)
!!!end do
!!!write(25000+iproc,*) '----------------------'
!!!
!!!
!!!        ! Calculate the gradient norm.
!!!        fnrm=0.d0
!!!        do iorb=1,orbs%norb
!!!            fnrm=fnrm+dnrm2(orbsig%norb, grad(1,iorb), 1)
!!!        end do
!!!
!!!        ! Precondition the gradient.
!!!        !!if(fnrm<1.d0) call preconditionGradient(iproc, nproc, orbsig, orbs, at, Ham, lagMat, onWhichAtomPhi, grad, it, evals)
!!!    
!!!
!!!        ! Write some informations to the screen, but only every 1000th iteration.
!!!        if(iproc==0 .and. mod(it,1)==0) write(*,'(x,a,es11.2,es22.13,es10.2)') 'fnrm, trace, mean alpha', &
!!!            fnrm, trace, meanAlpha
!!!        
!!!        ! Check for convergence.
!!!        if(fnrm<1.d-2) then
!!!            if(iproc==0) write(*,'(x,a,i0,a)') 'converged in ', it, ' iterations.'
!!!            if(iproc==0) write(*,'(3x,a,2es14.5)') 'Final values for fnrm, trace:', fnrm, trace
!!!            converged=.true.
!!!            infoCoeff=it
!!!            exit
!!!        end if
!!!  
!!!        ! Quit if the maximal number of iterations is reached.
!!!        if(it==lin%nItInguess) then
!!!            if(iproc==0) write(*,'(x,a,i0,a)') 'WARNING: not converged within ', it, &
!!!                ' iterations! Exiting loop due to limitations of iterations.'
!!!            if(iproc==0) write(*,'(x,a,2es15.7,f12.7)') 'Final values for fnrm, trace: ', fnrm, trace
!!!            infoCoeff=-1
!!!            exit
!!!        end if
!!!
!!!        ! Improve the coefficients (by steepet descent).
!!!        do iorb=1,orbs%norb
!!!            call daxpy(orbsig%norb, -alpha(iorb), grad(1,iorb), 1, coeff(1,iorb), 1)
!!!        end do
!!!    
!!!
!!!    end do iterLoop
!!!
!!!
!!!    if(iproc==0) write(*,'(x,a)') '===================================================================================='
!!!    !do iorb=1,orbs%norb
!!!    !    do jorb=1,orbsig%norb
!!!    !        write(999,'(2i8,es20.12)') iorb, jorb, coeff(jorb,iorb)
!!!    !    end do
!!!    !end do
!!!
!!!end if processIf
!!!
!!!
!!!! Now broadcast the result to all processes
!!!call mpi_bcast(coeff(1,1), orbsig%norb*orbs%norb, mpi_double_precision, 0, mpi_comm_world, ierr)
!!!call mpi_bcast(infoCoeff, 1, mpi_integer, 0, mpi_comm_world, ierr)
!!!call mpi_barrier(mpi_comm_world, ierr)
!!!
!!!
!!!
!!!  ! Build new linear combination
!!!  phi=0.d0
!!!  ist=1
!!!  do iorb=1,orbs%norb
!!!      jst=1
!!!      do jorb=1,orbsig%norb
!!!          call daxpy(nvctrp, coeff(jorb,iorb), chi(jst), 1, phi(ist), 1)
!!!          jst=jst+nvctrp
!!!      end do
!!!      ist=ist+nvctrp
!!!  end do
!!!
!!!  
!!!
!!!  ! Untranpose the orbitals.
!!!  allocate(chiw(orbs%npsidim+ndebug),stat=istat)
!!!  call memocc(istat, chiw, 'chiw', subname)
!!!  call untranspose_v(iproc, nproc, orbs, Glr%wfd, comms, phi, work=chiw)
!!!  iall=-product(shape(chiw))*kind(chiw)
!!!  deallocate(chiw, stat=istat)
!!!  call memocc(istat, iall, 'chiw', subname)
!!!
!!!
!!!  ! Deallocate the local arrays.
!!!  iall=-product(shape(Ham))*kind(Ham)
!!!  deallocate(Ham, stat=istat)
!!!  call memocc(istat, iall, 'Ham', subname)
!!!
!!!  iall=-product(shape(coeff))*kind(coeff)
!!!  deallocate(coeff, stat=istat)
!!!  call memocc(istat, iall, 'coeff', subname)
!!!
!!!  iall=-product(shape(alpha))*kind(alpha)
!!!  deallocate(alpha, stat=istat)
!!!  call memocc(istat, iall, 'alpha', subname)
!!!
!!!  iall=-product(shape(grad))*kind(grad)
!!!  deallocate(grad, stat=istat)
!!!  call memocc(istat, iall, 'grad', subname)
!!!
!!!  iall=-product(shape(gradOld))*kind(gradOld)
!!!  deallocate(gradOld, stat=istat)
!!!  call memocc(istat, iall, 'gradOld', subname)
!!!
!!!  iall=-product(shape(lagMat))*kind(lagMat)
!!!  deallocate(lagMat, stat=istat)
!!!  call memocc(istat, iall, 'lagMat', subname)
!!!
!!!  
!!!
!!!
!!!end subroutine buildLinearCombinations




subroutine orthonormalizeAtomicOrbitals(iproc, nproc, orbsig, commsig, Glr, chi)
!
! Purpose:
! ========
!  Orthonormalizes the atomic orbitals chi using a Lowedin orthonormalization.
!
! Calling arguments:
!    

use module_base
use module_types
use module_interfaces
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data),intent(in):: orbsig
type(communications_arrays),intent(in):: commsig
type(locreg_descriptors),intent(in):: Glr
real(8),dimension(orbsig%npsidim),intent(inout):: chi

! Local variables
integer:: iorb, jorb, istat, iall, lwork, info, nvctrp, ierr
real(8),dimension(:),allocatable:: eval, work
real(8),dimension(:,:),allocatable:: chiTemp, ovrlp
real(8),dimension(:,:,:),allocatable:: tempArr
real(8),dimension(:),pointer:: chiw
character(len=*),parameter:: subname='orthonormalizeAtomicOrbitals'


allocate(chiw(orbsig%npsidim),stat=istat)
call memocc(istat, chiw, 'chiw', subname)

call transpose_v(iproc, nproc, orbsig, Glr%wfd, commsig, chi, work=chiw)

allocate(ovrlp(orbsig%norb,orbsig%norb), stat=istat)
call memocc(istat, ovrlp, 'ovrlp', subname)


nvctrp=commsig%nvctr_par(iproc,1) ! 1 for k-points

! Calculate the overlap matrix
call dsyrk('l', 't', orbsig%norb, nvctrp, 1.d0, chi(1), nvctrp, &
     0.d0, ovrlp(1,1), orbsig%norb)
! MPI
call mpiallred(ovrlp(1,1), orbsig%norb**2, mpi_sum, mpi_comm_world, ierr)


allocate(eval(orbsig%norb), stat=istat)
call memocc(istat, eval, 'eval', subname)

! Diagonalize overlap matrix.
allocate(work(1), stat=istat)
call memocc(istat, work, 'work', subname)
call dsyev('v', 'l', orbsig%norb, ovrlp(1,1), orbsig%norb, eval, &
     work, -1, info)
lwork=work(1)
iall=-product(shape(work))*kind(work)
deallocate(work, stat=istat)
call memocc(istat, iall, 'work', subname)
allocate(work(lwork), stat=istat)
call memocc(istat, work, 'work', subname)
call dsyev('v', 'l', orbsig%norb, ovrlp(1,1), orbsig%norb, eval, &
     work, lwork, info)

! Calculate S^{-1/2}. 
! First calulate ovrlp*diag(1/sqrt(eval)) (ovrlp is the diagonalized overlap
! matrix and diag(1/sqrt(eval)) the diagonal matrix consisting of the inverse square roots of the eigenvalues...
allocate(tempArr(orbsig%norb,orbsig%norb,2), stat=istat)
call memocc(istat, tempArr, 'tempArr', subname)
do iorb=1,orbsig%norb
    do jorb=1,orbsig%norb
        tempArr(jorb,iorb,1)=ovrlp(jorb,iorb)*1.d0/sqrt(eval(iorb))
    end do
end do

! ...and now apply the diagonalized overlap matrix to the matrix constructed above.
! This will give S^{-1/2}.
call dgemm('n', 't', orbsig%norb, orbsig%norb, orbsig%norb, 1.d0, ovrlp(1,1), &
     orbsig%norb, tempArr(1,1,1), orbsig%norb, 0.d0, &
     tempArr(1,1,2), orbsig%norb)

! Now calculate the orthonormal orbitals by applying S^{-1/2} to the orbitals.
! This requires the use of a temporary variable phidTemp.
allocate(chiTemp(nvctrp,1:orbsig%norb), stat=istat)
call memocc(istat, chiTemp, 'chiTemp', subname)
call dgemm('n', 'n', nvctrp, orbsig%norb, orbsig%norb, 1.d0, chi(1), &
     nvctrp, tempArr(1,1,2),  orbsig%norb, 0.d0, &
     chiTemp(1,1), nvctrp)

! Now copy the orbitals from the temporary variable to phid.
call dcopy(orbsig%norb*nvctrp, chiTemp(1,1), 1, chi(1), 1)


call untranspose_v(iproc, nproc, orbsig, Glr%wfd, commsig, chi, work=chiw)


iall=-product(shape(ovrlp))*kind(ovrlp)
deallocate(ovrlp, stat=istat)
call memocc(istat, iall, 'ovrlp', subname)

iall=-product(shape(work))*kind(work)
deallocate(work, stat=istat)
call memocc(istat, iall, 'work', subname)

iall=-product(shape(eval))*kind(eval)
deallocate(eval, stat=istat)
call memocc(istat, iall, 'eval', subname)

iall=-product(shape(chiTemp))*kind(chiTemp)
deallocate(chiTemp, stat=istat)
call memocc(istat, iall, 'chiTemp', subname)

iall=-product(shape(tempArr))*kind(tempArr)
deallocate(tempArr, stat=istat)
call memocc(istat, iall, 'tempArr', subname)

iall=-product(shape(chiw))*kind(chiw)
deallocate(chiw, stat=istat)
call memocc(istat, iall, 'chiw', subname)


end subroutine orthonormalizeAtomicOrbitals





subroutine orthonormalizeCoefficients(orbs, orbsig, coeff)
use module_base
use module_types
implicit none

! Calling arguments
type(orbitals_data),intent(in):: orbs, orbsig
real(8),dimension(orbsig%norb,orbs%norb),intent(inout):: coeff

! Local variables
integer:: iorb, jorb, istat, iall, lwork, info
real(8),dimension(:),allocatable:: work, eval
real(8),dimension(:,:),allocatable:: ovrlp, coeffTemp
real(8),dimension(:,:,:),allocatable:: tempArr
character(len=*),parameter:: subname='orthonormalizeCoefficients'
real(8):: ddot

        allocate(ovrlp(orbs%norb,orbs%norb), stat=istat)
        call memocc(istat, ovrlp, 'ovrlp', subname)
        allocate(eval(orbs%norb), stat=istat)
        call memocc(istat, eval, 'eval', subname)
        allocate(tempArr(orbs%norb,orbs%norb,2), stat=istat)
        call memocc(istat, tempArr, 'tempArr', subname)
        allocate(coeffTemp(orbsig%norb,orbs%norb), stat=istat)
        call memocc(istat, coeffTemp, 'coeffTemp', subname)


        !!! Orthonormalize the coefficient vectors (Gram-Schmidt).
        !!do iorb=1,orbs%norb
        !!    do jorb=1,iorb-1
        !!        tt=ddot(orbsig%norb, coeff(1,iorb), 1, coeff(1,jorb), 1)
        !!        call daxpy(orbsig%norb, -tt, coeff(1,jorb), 1, coeff(1,iorb), 1)
        !!    end do
        !!    tt=dnrm2(orbsig%norb, coeff(1,iorb), 1)
        !!    call dscal(orbsig%norb, 1/tt, coeff(1,iorb), 1)
        !!end do

        !!! Orthonormalize the coefficient vectors (Loewdin).
        !!do iorb=1,orbs%norb
        !!    do jorb=1,orbs%norb
        !!        ovrlp(iorb,jorb)=ddot(orbsig%norb, coeff(1,iorb), 1, coeff(1,jorb), 1)
        !!    end do
        !!end do

        allocate(work(1), stat=istat)
        call memocc(istat, work, 'work', subname)
        call dsyev('v', 'l', orbs%norb, ovrlp(1,1), orbs%norb, eval, work, -1, info)
        lwork=work(1)
        iall=-product(shape(work))*kind(work)
        deallocate(work, stat=istat)
        call memocc(istat, iall, 'work', subname)
        allocate(work(lwork), stat=istat)
        call memocc(istat, work, 'work', subname)
        call dsyev('v', 'l', orbs%norb, ovrlp(1,1), orbs%norb, eval, work, lwork, info)
        iall=-product(shape(work))*kind(work)
        deallocate(work, stat=istat)
        call memocc(istat, iall, 'work', subname)

        ! Calculate S^{-1/2}. 
        ! First calulate ovrlp*diag(1/sqrt(evall)) (ovrlp is the diagonalized overlap
        ! matrix and diag(1/sqrt(evall)) the diagonal matrix consisting of the inverse square roots of the eigenvalues...
        do iorb=1,orbs%norb
            do jorb=1,orbs%norb
                tempArr(jorb,iorb,1)=ovrlp(jorb,iorb)*1.d0/sqrt(eval(iorb))
            end do
        end do

        ! ...and now apply the diagonalized overlap matrix to the matrix constructed above.
        ! This will give S^{-1/2}.
        call dgemm('n', 't', orbs%norb, orbs%norb, orbs%norb, 1.d0, ovrlp(1,1), &
        orbs%norb, tempArr(1,1,1), orbs%norb, 0.d0, &
        tempArr(1,1,2), orbs%norb)

        ! Now calculate the orthonormal orbitals by applying S^{-1/2} to the orbitals.
        ! This requires the use of a temporary variable phidTemp.
        call dgemm('n', 'n', orbsig%norb, orbs%norb, orbs%norb, 1.d0, coeff(1,1), &
             orbsig%norb, tempArr(1,1,2), orbs%norb, 0.d0, &
             coeffTemp(1,1), orbsig%norb)
        
        ! Now copy the orbitals from the temporary variable to phid.
        call dcopy(orbs%norb*orbsig%norb, coeffTemp(1,1), 1, coeff(1,1), 1)

        iall=-product(shape(ovrlp))*kind(ovrlp)
        deallocate(ovrlp, stat=istat)
        call memocc(istat, iall, 'ovrlp', subname)

        iall=-product(shape(eval))*kind(eval)
        deallocate(eval, stat=istat)
        call memocc(istat, iall, 'eval', subname)

        iall=-product(shape(tempArr))*kind(tempArr)
        deallocate(tempArr, stat=istat)
        call memocc(istat, iall, 'tempArr', subname)

        iall=-product(shape(coeffTemp))*kind(coeffTemp)
        deallocate(coeffTemp, stat=istat)
        call memocc(istat, iall, 'coeffTemp', subname)

end subroutine orthonormalizeCoefficients




subroutine preconditionGradient(iproc, nproc, orbsig, orbs, at, Ham, lagMat, onWhichAtomPhi, grad, it, evals)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, it
type(orbitals_data),intent(in):: orbsig, orbs
type(atoms_data),intent(in):: at
real(8),dimension(orbsig%norb,orbsig%norb,at%nat),intent(in):: Ham
real(8),dimension(orbs%norb,orbs%norb),intent(in):: lagMat
integer,dimension(orbs%norb),intent(in):: onWhichAtomPhi
real(8),dimension(orbsig%norb,orbs%norb),intent(inout):: grad
real(8),dimension(orbs%norb),intent(in):: evals

! Local variables
integer:: iorb, jorb, korb, iiAt, info, istat, iall
integer,dimension(:),allocatable:: ipiv
real(8),dimension(:,:,:),allocatable:: matc
real(8),dimension(:,:),allocatable:: solc
character(len=*),parameter:: subname='preconditionGradient'
real(8):: ddot


allocate(ipiv(orbsig%norb), stat=istat)
call memocc(istat, ipiv, 'ipiv', subname)
allocate(matc(2,orbsig%norb,orbsig%norb), stat=istat)
call memocc(istat, matc, 'matc', subname)
allocate(solc(2,orbsig%norb), stat=istat)
call memocc(istat, solc, 'solc', subname)

do iorb=1,orbs%norb
    iiAt=onWhichAtomPhi(iorb)
    ! Build matrix that has to be inverted
!write(*,*) 'iorb, lagMat(iorb,iorb)', iorb, lagMat(iorb,iorb)
write(*,'(a,i5,4x,2es15.7)') 'iorb, evals(iorb), lagMat(iorb,iorb)', iorb, evals(iorb), lagMat(iorb,iorb)
    do jorb=1,orbsig%norb
        do korb=1,orbsig%norb
            matc(1,korb,jorb)=Ham(korb,jorb,iiAt)
            matc(2,korb,jorb)=0.d0
        end do
        !matc(1,jorb,jorb)=matc(1,jorb,jorb)+lagMat(iorb,iorb)**2
        !matc(1,jorb,jorb)=matc(1,jorb,jorb)-.5d0
        matc(1,jorb,jorb)=matc(1,jorb,jorb)-evals(iorb)
        !matc(1,jorb,jorb)=1.d0
        matc(2,jorb,jorb)=-1.d-2
        solc(1,jorb)=grad(jorb,iiAt)
        solc(2,jorb)=0.d0
    end do
    !do jorb=1,orbsig%norb
    !    mat(jorb,jorb)=mat(jorb,jorb)-lagMat(iorb,iorb)
    !end do
    !call dcopy(orbsig%norb, grad(1,iorb), 1, sol(1), 1)
write(100+iorb,'(a,i0,es15.7)') 'grad ', it, ddot(orbsig%norb, grad(1,iorb), 1, grad(1,iorb), 1)
do iall=1,orbsig%norb
  write(100+iorb,*) iall, grad(iall,iorb)
end do
    !call dgesv(orbsig%norb, 1, mat(1,1), orbsig%norb, ipiv, grad(1,iorb), orbsig%norb, info)
    call zgesv(orbsig%norb, 1, matc(1,1,1), orbsig%norb, ipiv, solc(1,1), orbsig%norb, info)
    if(info/=0) then
        write(*,'(x,a,i0)') 'ERROR in zgesv (subroutine preconditionGradient), infocode=', info
    end if
    do jorb=1,orbsig%norb
        grad(jorb,iiAt)=solc(1,jorb)
    end do
    !call dscal(orbsig%norb, -1.d0, grad(1,iorb), 1)
 write(200+iorb,'(a,i0,es15.7)') 'grad ', it, ddot(orbsig%norb, grad(1,iorb), 1, grad(1,iorb), 1)
 do iall=1,orbsig%norb
   write(200+iorb,*) iall, grad(iall,iorb)
 end do
end do

iall=-product(shape(ipiv))*kind(ipiv)
deallocate(ipiv, stat=istat)
call memocc(istat, iall, 'ipiv', subname)

iall=-product(shape(matc))*kind(matc)
deallocate(matc, stat=istat)
call memocc(istat, iall, 'matc', subname)

iall=-product(shape(solc))*kind(solc)
deallocate(solc, stat=istat)
call memocc(istat, iall, 'solc', subname)

end subroutine preconditionGradient




subroutine transposeInguess(iproc, ip, newComm, chi)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, newComm
type(inguessParameters),intent(in):: ip
real(8),dimension(ip%norbtotPad*ip%norb_par(iproc)),intent(inout):: chi

! Local variables
integer:: ii, jj, jproc, iorb, j, ierr, istat, iall
real(8),dimension(:),allocatable:: chiw
character(len=*),parameter:: subname='transposeInguess'


  ! Allocate the work array.
  allocate(chiw(ip%sizeWork), stat=istat)
  call memocc(istat, chiw, 'chiw', subname)
  chiw=0.d0
  
  ! Rearranges the elements on the processor to enable the transposition with only one call
  ! to mpi_alltoallv.
  jj=1
  ii=0
  do jproc=0,ip%nproc-1
     do iorb=0,ip%norb_par(iproc)-1
        ! Copy the non-zero parts
        !write(*,'(a,6i8)') 'iproc, jproc, iorb, ii, ip%norbtotPad, ii+iorb*ip%norbtotPad+1', iproc, jproc, iorb, ii, ip%norbtotPad, ii+iorb*ip%norbtotPad+1
        call dcopy(ip%nvctrp_nz(jproc), chi(ii+iorb*ip%norbtotPad+1), 1, chiw(jj), 1)
        jj=jj+ip%nvctrp_nz(jproc)
        do j=ip%nvctrp_nz(jproc)+1,ip%nvctrp
           ! "Copy" the zeros. This happens only if ip%nvctrp_nz(jproc) < ip%nvctrp
           chiw(jj)=0.d0
           jj=jj+1
        end do
     end do
     ii=ii+ip%nvctrp_nz(jproc)
  end do

  ! Communicate the vectors.
  !call mpi_alltoallv(chiw(1), ip%sendcounts, ip%senddispls, mpi_double_precision, chi(1), &
  !     ip%recvcounts, ip%recvdispls, mpi_double_precision, mpi_comm_world, ierr)
  call mpi_alltoallv(chiw(1), ip%sendcounts, ip%senddispls, mpi_double_precision, chi(1), &
       ip%recvcounts, ip%recvdispls, mpi_double_precision, newComm, ierr)

  ! Dellocate the work array.
  iall=-product(shape(chiw))*kind(chiw)
  deallocate(chiw, stat=istat)
  call memocc(istat, iall, 'chiw', subname)

end subroutine transposeInguess





subroutine untransposeInguess(iproc, ip, newComm, chi)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, newComm
type(inguessParameters),intent(in):: ip
real(8),dimension(ip%norbtotPad*ip%norb_par(iproc)),intent(inout):: chi

! Local variables
integer:: ii, jj, jproc, iorb, istat, iall, ierr
real(8),dimension(:),allocatable:: chiw
character(len=*),parameter:: subname='untransposeInguess'


  ! Allocate the work array.
  allocate(chiw(ip%sizeWork), stat=istat)
  call memocc(istat, chiw, 'chiw', subname)

  ! Communicate the data.
  !call mpi_alltoallv(chi(1), ip%recvcounts, ip%recvdispls, mpi_double_precision, chiw(1), &
  !     ip%sendcounts, ip%senddispls, mpi_double_precision, mpi_comm_world, ierr)
  call mpi_alltoallv(chi(1), ip%recvcounts, ip%recvdispls, mpi_double_precision, chiw(1), &
       ip%sendcounts, ip%senddispls, mpi_double_precision, newComm, ierr)

  ! Rearrange back the elements.
  chi=0.d0
  ii=0
  jj=1
  do jproc=0,ip%nproc-1
     do iorb=0,ip%norb_par(iproc)-1
        call dcopy(ip%nvctrp, chiw(jj), 1, chi(ii+iorb*ip%norbtotPad+1), 1)
        jj=jj+ip%nvctrp
     end do
     ii=ii+ip%nvctrp_nz(jproc)
  end do

     !!ii=0
     !!jj=1
     !!do i=nprocSt,nprocSt+nproc-1
     !!   do iorb=0,norbp-1
     !!      call dcopy(norbtotp*nspinor, psiW(jj), 1, psi(ii+iorb*norbtot*nspinor+1), 1)
     !!      jj=jj+norbtotp*nspinor
     !!   end do
     !!   ii=ii+norbtotpArr(i)*nspinor
     !!end do


  ! Dellocate the work array.
  iall=-product(shape(chiw))*kind(chiw)
  deallocate(chiw, stat=istat)
  call memocc(istat, iall, 'chiw', subname)

end subroutine untransposeInguess




!subroutine initializeInguessParameters(iproc, nproc, orbs, orbsig, lin, ip)
subroutine initializeInguessParameters(iproc, orbs, orbsig, newComm, ip)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc
type(orbitals_data),intent(in):: orbs, orbsig
integer,intent(in):: newComm
type(inguessParameters),intent(inout):: ip

! Local variables
integer:: ii, kk, jproc, istat, ierr, norbTarget
real(8):: tt
character(len=*),parameter:: subname='initializeInguessParameters'


  


  ip%norb=orbs%norb
  ip%norbtot=orbsig%norb

  ! In order to symplify the transposing/untransposing, the orbitals are padded with zeros such that 
  ! they can be distributed evenly over all processes when being transposed. The new length of the 
  ! orbitals after this padding is then given by ip%norbtotPad.
  ip%norbtotPad=ip%norbtot
  do
      if(mod(ip%norbtotPad, ip%nproc)==0) exit
      ip%norbtotPad=ip%norbtotPad+1
  end do

  ! Distribute the orbitals among the processes.
  allocate(ip%norb_par(0:ip%nproc-1), stat=istat)
  call memocc(istat, ip%norb_par, 'ip%norb_par', subname)
  ip%norb_par=0
  tt=dble(ip%norb)/dble(ip%nproc)
  ii=floor(tt)
  ! ii is now the number of orbitals that every process has. Distribute the remaining ones.
  ip%norb_par(0:ip%nproc-1)=ii
  kk=ip%norb-ip%nproc*ii
  ip%norb_par(0:kk-1)=ii+1

  ! Starting orbital for each process
  ii=0
  do jproc=0,iproc-1
     ii=ii+ip%norb_par(jproc)
  end do
  !reference orbital for process
  ip%isorb=ii


  ! Calculate the number of elements that each process has when the vectors are transposed.
  ! nvctrp is the total number, nvctrp_nz is the nonzero numbers.
  allocate(ip%nvctrp_nz(0:ip%nproc-1), stat=istat)
  call memocc(istat, ip%nvctrp_nz, 'ip%nvctrp_nz', subname)
  tt=ip%norbtot/dble(ip%nproc)
  ii=floor(tt)
  ! ii is now the number of elements that every process has. Distribute the remaining ones.
  ip%nvctrp_nz=ii
  kk=ip%norbtot-ip%nproc*ii
  ip%nvctrp_nz(0:kk-1)=ii+1
  ! Check wheter this distribution is correct
  ii=0
  do jproc=0,ip%nproc-1
     ii=ii+ip%nvctrp_nz(jproc)
  end do
  if(ii/=ip%norbtot) then
     if(iproc==0) write(*,'(3x,a)') 'ERROR: wrong partition of ip%norbtot!'
     call mpi_barrier(newComm, ierr)
     stop
  end if

  ! With the padded zeros, the elemts can be distributed evenly.
  ip%nvctrp=ip%norbtotPad/ip%nproc

  ! Define the values for the mpi_alltoallv.
  ! sendcounts: number of elements that a given  process sends to another process.
  ! senddispls: offset of the starting index on a given process for the send operation to another process.
  allocate(ip%sendcounts(0:ip%nproc-1), stat=istat)
  call memocc(istat, ip%sendcounts, 'ip%sendcounts', subname)
  allocate(ip%senddispls(0:ip%nproc-1), stat=istat)
  call memocc(istat, ip%senddispls, 'ip%senddispls', subname)
  ii=0
  do jproc=0,ip%nproc-1
      ip%sendcounts(jproc)=ip%nvctrp*ip%norb_par(iproc)
      ip%senddispls(jproc)=ii
      ii=ii+ip%sendcounts(jproc)
  end do
  ! recvcounts: number of elements that a given process receives from another process.
  ! recvdispls: offset of the starting index on a given process for the receive operation from another process.
  allocate(ip%recvcounts(0:ip%nproc-1), stat=istat)
  call memocc(istat, ip%recvcounts, 'ip%recvcounts', subname)
  allocate(ip%recvdispls(0:ip%nproc-1), stat=istat)
  call memocc(istat, ip%recvdispls, 'ip%recvdispls', subname)
  ii=0
  do jproc=0,ip%nproc-1
      ip%recvcounts(jproc)=ip%nvctrp*ip%norb_par(jproc)
      ip%recvdispls(jproc)=ii
      ii=ii+ip%recvcounts(jproc)
  end do

  ! Determine the size of the work array needed for the transposition.
  ip%sizeWork=max(ip%norbtotPad*ip%norb_par(iproc),sum(ip%recvcounts(:)))

end subroutine initializeInguessParameters





subroutine orthonormalizeCoefficients_parallel(iproc, ip, newComm, coeff)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, newComm
type(inguessParameters),intent(in):: ip
real(8),dimension(ip%nvctrp,ip%norb),intent(inout):: coeff

! Local variables
integer:: iorb, jorb, istat, iall, lwork, info, ierr
real(8),dimension(:),allocatable:: work, eval
real(8),dimension(:,:),allocatable:: ovrlp, coeffTemp
real(8),dimension(:,:,:),allocatable:: tempArr
character(len=*),parameter:: subname='orthonormalizeCoefficients'
real(8):: ddot


        allocate(ovrlp(ip%norb,ip%norb), stat=istat)
        call memocc(istat, ovrlp, 'ovrlp', subname)
        allocate(eval(ip%norb), stat=istat)
        call memocc(istat, eval, 'eval', subname)
        allocate(tempArr(ip%norb,ip%norb,2), stat=istat)
        call memocc(istat, tempArr, 'tempArr', subname)
        allocate(coeffTemp(ip%nvctrp,ip%norb), stat=istat)
        call memocc(istat, coeffTemp, 'coeffTemp', subname)


        !!! Orthonormalize the coefficient vectors (Gram-Schmidt).
        !!do iorb=1,orbs%norb
        !!    do jorb=1,iorb-1
        !!        tt=ddot(orbsig%norb, coeff(1,iorb), 1, coeff(1,jorb), 1)
        !!        call daxpy(orbsig%norb, -tt, coeff(1,jorb), 1, coeff(1,iorb), 1)
        !!    end do
        !!    tt=dnrm2(orbsig%norb, coeff(1,iorb), 1)
        !!    call dscal(orbsig%norb, 1/tt, coeff(1,iorb), 1)
        !!end do


        ! Orthonormalize the coefficient vectors (Loewdin).
        do iorb=1,ip%norb
            do jorb=1,ip%norb
                !ovrlp(iorb,jorb)=ddot(orbsig%norb, coeff(1,iorb), 1, coeff(1,jorb), 1)
                ovrlp(iorb,jorb)=ddot(ip%nvctrp_nz(iproc), coeff(1,iorb), 1, coeff(1,jorb), 1)
            end do
        end do

        ! Sum up over all processes.
        !if(ip%nproc>1) call mpiallred(ovrlp(1,1), ip%norb**2, mpi_sum, mpi_comm_world, ierr)
        if(ip%nproc>1) call mpiallred(ovrlp(1,1), ip%norb**2, mpi_sum, newComm, ierr)

        !!do iorb=1,ip%norb
        !!    do jorb=1,ip%norb
        !!        if(abs(ovrlp(iorb,jorb)-ovrlp(jorb,iorb))>1.d-10) stop 'not symmetric'
        !!        write(3000+iproc,*) iorb, jorb, ovrlp(jorb,iorb)
        !!    end do
        !!end do
        !!write(3000+iproc,*) '=================================='

        !!allocate(work(1), stat=istat)
        !!call memocc(istat, work, 'work', subname)
        !!call dsyev('v', 'l', ip%norb, ovrlp(1,1), ip%norb, eval, work, -1, info)
        !!lwork=work(1)
        !!iall=-product(shape(work))*kind(work)
        !!deallocate(work, stat=istat)
        !!call memocc(istat, iall, 'work', subname)
        lwork=100*ip%norb
        allocate(work(lwork), stat=istat)
        call memocc(istat, work, 'work', subname)
        call dsyev('v', 'l', ip%norb, ovrlp(1,1), ip%norb, eval, work, lwork, info)
        if(info/=0) then
            write(*,'(a,i0)') 'ERROR in dsyev, info=', info
        end if
        iall=-product(shape(work))*kind(work)
        deallocate(work, stat=istat)
        call memocc(istat, iall, 'work', subname)

        ! Calculate S^{-1/2}. 
        ! First calulate ovrlp*diag(1/sqrt(evall)) (ovrlp is the diagonalized overlap
        ! matrix and diag(1/sqrt(evall)) the diagonal matrix consisting of the inverse square roots of the eigenvalues...
        do iorb=1,ip%norb
            do jorb=1,ip%norb
                tempArr(jorb,iorb,1)=ovrlp(jorb,iorb)*1.d0/sqrt(eval(iorb))
            end do
        end do

        ! ...and now apply the diagonalized overlap matrix to the matrix constructed above.
        ! This will give S^{-1/2}.
        call dgemm('n', 't', ip%norb, ip%norb, ip%norb, 1.d0, ovrlp(1,1), &
             ip%norb, tempArr(1,1,1), ip%norb, 0.d0, &
             tempArr(1,1,2), ip%norb)

        ! Now calculate the orthonormal orbitals by applying S^{-1/2} to the orbitals.
        ! This requires the use of a temporary variable phidTemp.
        coeffTemp=0.d0
        call dgemm('n', 'n', ip%nvctrp_nz(iproc), ip%norb, ip%norb, 1.d0, coeff(1,1), &
             ip%nvctrp, tempArr(1,1,2), ip%norb, 0.d0, &
             coeffTemp(1,1), ip%nvctrp)
        
        ! Now copy the orbitals from the temporary variable to phid.
        call dcopy(ip%norb*ip%nvctrp, coeffTemp(1,1), 1, coeff(1,1), 1)

        iall=-product(shape(ovrlp))*kind(ovrlp)
        deallocate(ovrlp, stat=istat)
        call memocc(istat, iall, 'ovrlp', subname)

        iall=-product(shape(eval))*kind(eval)
        deallocate(eval, stat=istat)
        call memocc(istat, iall, 'eval', subname)

        iall=-product(shape(tempArr))*kind(tempArr)
        deallocate(tempArr, stat=istat)
        call memocc(istat, iall, 'tempArr', subname)

        iall=-product(shape(coeffTemp))*kind(coeffTemp)
        deallocate(coeffTemp, stat=istat)
        call memocc(istat, iall, 'coeffTemp', subname)

end subroutine orthonormalizeCoefficients_parallel




subroutine getHamiltonianMatrix(iproc, nproc, lzdig, Glr, onWhichAtom, onWhichAtomp, chi, hchi, ham)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(linear_zone_descriptors),intent(in):: lzdig
type(locreg_descriptors),intent(in):: Glr
integer,dimension(lzdig%orbs%norb),intent(in):: onWhichAtom
integer,dimension(lzdig%orbs%norbp),intent(in):: onWhichAtomp
real(8),dimension(lzdig%orbs%npsidim),intent(in):: chi, hchi
real(8),dimension(lzdig%orbs%norb,lzdig%orbs%norb),intent(out):: ham

! Local variables
integer:: sizeChi, istat, iorb, ilr, iall, ind1, ind2, ldim, gdim
type(overlapParameters):: op
type(p2pCommsOrthonormality):: comon
real(8),dimension(:),allocatable:: lchi, lhchi, lphiovrlp
character(len=*),parameter:: subname='getHamiltonianMatrix'



! Initialize the parameters for calculating the matrix.
call initCommsOrtho(iproc, nproc, lzdig, onWhichAtom, op, comon)

allocate(lphiovrlp(op%ndim_lphiovrlp), stat=istat)
call memocc(istat, lphiovrlp, 'lphiovrlp',subname)



! Calculate the dimension of the wave function for each process.
sizeChi=0
do iorb=1,lzdig%orbs%norbp
    ilr=onWhichAtomp(iorb)
    sizeChi = sizeChi + (lzdig%Llr(ilr)%wfd%nvctr_c+7*lzdig%Llr(ilr)%wfd%nvctr_f)
end do

allocate(lchi(sizeChi), stat=istat)
call memocc(istat, lchi, 'lchi', subname)
allocate(lhchi(sizeChi), stat=istat)
call memocc(istat, lhchi, 'lhchi', subname)

! Transform chi to the localization region. This is not needed if we really habe O(N).
ind1=1
ind2=1
do iorb=1,lzdig%orbs%norbp
    ilr = onWhichAtomp(iorb)
    ldim=lzdig%Llr(ilr)%wfd%nvctr_c+7*lzdig%Llr(ilr)%wfd%nvctr_f
    gdim=Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
    call psi_to_locreg2(iproc, nproc, ldim, gdim, lzdig%Llr(ilr), Glr, chi(ind1), lchi(ind2))
    ind1=ind1+Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
    ind2=ind2+lzdig%Llr(ilr)%wfd%nvctr_c+7*lzdig%Llr(ilr)%wfd%nvctr_f
end do

! Put lphi in the sendbuffer, i.e. lphi will be sent to other processes' receive buffer.
call extractOrbital(iproc, nproc, lzdig%orbs, sizeChi, onWhichAtom, lzdig, op, lchi, comon)
call postCommsOverlap(iproc, nproc, comon)
call gatherOrbitals(iproc, nproc, comon)
! Put lhphi to the sendbuffer, so we can the calculate <lphi|lhphi>
call extractOrbital(iproc, nproc, lzdig%orbs, sizeChi, onWhichAtom, lzdig, op, lhchi, comon)
call calculateOverlapMatrix2(iproc, nproc, lzdig%orbs, op, comon, onWhichAtom, ham)


iall=-product(shape(lchi))*kind(lchi)
deallocate(lchi, stat=istat)
call memocc(istat, iall, 'lchi', subname)

iall=-product(shape(lhchi))*kind(lhchi)
deallocate(lhchi, stat=istat)
call memocc(istat, iall, 'lhchi', subname)


end subroutine getHamiltonianMatrix

