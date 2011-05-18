!>   input guess wavefunction diagonalization
subroutine inputguessConfinement(iproc, nproc, at, &
     comms, Glr, input, lin, rxyz, n3p, rhopot, rhocore, pot_ion,&
     nlpspd, proj, pkernel, pkernelseq, &
     nscatterarr, ngatherarr, potshortcut, irrzon, phnons, GPU, &
     phi)
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
  real(8),dimension(lin%orbs%npsidim),intent(out):: phi
  !local variables
  type(gaussian_basis):: G !basis for davidson IG
  character(len=*), parameter :: subname='input_wf_diag'
  logical :: switchGPUconv,switchOCLconv
  integer :: i_stat,i_all,iat,nspin_ig,iorb,idum=0, nvirt, norbat
  real(kind=4) :: tt,builtin_rand
  real(gp) :: hxh,hyh,hzh,eks,eexcu,vexcu,epot_sum,ekin_sum,ehart,eexctX,eproj_sum,etol,accurex
  type(orbitals_data) :: orbsig
  type(communications_arrays) :: commsig
  integer, dimension(:,:), allocatable :: norbsc_arr
  real(wp), dimension(:), allocatable :: potxc
  real(gp), dimension(:), allocatable :: locrad
  type(locreg_descriptors), dimension(:), allocatable :: Llr
  real(wp), dimension(:), pointer :: pot
  real(wp), dimension(:,:,:), pointer :: psigau
type(orbitals_data):: orbsLIN
type(communications_arrays):: commsLIN
real(8),dimension(:,:),allocatable:: HamSmall
real(8),dimension(:,:),pointer:: hpsiForAllLocregs
real(8),dimension(:),allocatable:: eval
integer:: istat
real(8),dimension(:),pointer:: phiWorkPointer
real(8),dimension(:),allocatable:: chi
real(8),dimension(:,:),allocatable:: hchi
integer,dimension(:),allocatable:: onWhichAtom
integer, parameter :: nmax=6,lmax=3,noccmax=2,nelecmax=32
integer, dimension(lmax+1) :: nl
real(gp), dimension(noccmax,lmax+1) :: occup



  allocate(norbsc_arr(at%natsc+1,input%nspin+ndebug),stat=i_stat)
  call memocc(i_stat,norbsc_arr,'norbsc_arr',subname)
  allocate(locrad(at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,locrad,'locrad',subname)


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

  ! Determine how many atomic orbitals we have
  norbat=0
  do iat=1,at%nat
      call count_atomic_shells(lmax+1, noccmax, nelecmax, input%nspin, lin%orbs%nspinor, at%aocc(1,iat), occup, nl)
      !write(*,'(a,i4,2x,10i4)') 'iat, nl', iat, nl
      norbat=norbat+(nl(1)+3*nl(2)+5*nl(3)+7*nl(4))
  end do
  write(*,*) 'norbat',norbat
  orbsig%norb=norbat ! This will be assigned later anyway, but here it is needed
                     ! to pass the array onWhichAtom to the subroutine inputguess_gaussian_orbitals
  allocate(onWhichAtom(orbsig%norb),stat=i_stat)
  call memocc(i_stat, onWhichAtom, 'onWhichAtom', subname)


  nvirt=0
  call inputguess_gaussian_orbitals_withOnWhichAtom(iproc,nproc,at,rxyz,Glr,nvirt,nspin_ig,&
       lin%orbs,orbsig,norbsc_arr,locrad,G,psigau,eks,onWhichAtom)
  write(*,*) 'orbsig%norb',orbsig%norb

  !allocate communications arrays for inputguess orbitals
  !call allocate_comms(nproc,orbsig,commsig,subname)
  call orbitals_communicators(iproc,nproc,Glr,orbsig,commsig)  

  hxh=.5_gp*input%hx
  hyh=.5_gp*input%hy
  hzh=.5_gp*input%hz

  !check the communication distribution
  !call check_communications(iproc,nproc,orbsig,Glr,commsig)

  !once the wavefunction coefficients are known perform a set 
  !of nonblocking send-receive operations to calculate overlap matrices

!!!  !create mpirequests array for controlling the success of the send-receive operation
!!!  allocate(mpirequests(nproc-1+ndebug),stat=i_stat)
!!!  call memocc(i_stat,mpirequests,'mpirequests',subname)
!!!
!!!  call nonblocking_transposition(iproc,nproc,G%ncoeff,orbsig%isorb+orbsig%norbp,&
!!!       orbsig%nspinor,psigau,orbsig%norb_par,mpirequests)

  !experimental part for building the localisation regions
  if (at%geocode == 'F') then
     !allocate the array of localisation regions
     allocate(Llr(at%nat+ndebug),stat=i_stat)
     !call memocc(i_stat,Llr,'Llr',subname)

     !print *,'locrad',locrad

     call determine_locreg(at%nat,rxyz,locrad,input%hx,input%hy,input%hz,Glr,Llr)

     do iat=1,at%nat
        call deallocate_lr(Llr(iat),subname)
!!$        call deallocate_wfd(Llr(iat)%wfd,subname)
!!$        if (Llr(iat)%geocode=='F') then
!!$           call deallocate_bounds(Llr(iat)%bounds,subname)
!!$        end if
     end do

     !i_all=-product(shape(Llr))*kind(Llr)
     deallocate(Llr,stat=i_stat) !these allocation are special
     !call memocc(i_stat,i_all,'Llr',subname)
  end if

  !allocate the wavefunction in the transposed way to avoid allocations/deallocations
  !allocate(psi(orbsig%npsidim+ndebug),stat=i_stat)
  !call memocc(i_stat,psi,'psi',subname)
  allocate(chi(orbsig%npsidim+ndebug),stat=i_stat)
  call memocc(i_stat,chi,'chi',subname)

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


  !use only the part of the arrays for building the hamiltonian matrix
  call gaussians_to_wavelets_new(iproc,nproc,Glr,orbsig,input%hx,input%hy,input%hz,G,&
       psigau(1,1,min(orbsig%isorb+1,orbsig%norb)),chi)


  i_all=-product(shape(locrad))*kind(locrad)
  deallocate(locrad,stat=i_stat)
  call memocc(i_stat,i_all,'locrad',subname)

  !application of the hamiltonian for gaussian based treatment
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

!!!  if (nproc == 1) then
!!!     !calculate the overlap matrix as well as the kinetic overlap
!!!     !in view of complete gaussian calculation
!!!     allocate(ovrlp(G%ncoeff*G%ncoeff),stat=i_stat)
!!!     call memocc(i_stat,ovrlp,'ovrlp',subname)
!!!     allocate(tmp(G%ncoeff,orbsig%norb),stat=i_stat)
!!!     call memocc(i_stat,tmp,'tmp',subname)
!!!     allocate(smat(orbsig%norb,orbsig%norb),stat=i_stat)
!!!     call memocc(i_stat,smat,'smat',subname)
!!!
!!!     !overlap calculation of the gaussian matrix
!!!     call gaussian_overlap(G,G,ovrlp)
!!!     call dsymm('L','U',G%ncoeff,orbsig%norb,1.0_gp,ovrlp(1),G%ncoeff,&
!!!          gaucoeff(1,1),G%ncoeff,0.d0,tmp(1,1),G%ncoeff)
!!!
!!!     call gemm('T','N',orbsig%norb,orbsig%norb,G%ncoeff,1.0_gp,&
!!!          gaucoeff(1,1),G%ncoeff,tmp(1,1),G%ncoeff,0.0_wp,smat(1,1),orbsig%norb)
!!!
!!!     !print overlap matrices
!!!     do i=1,orbsig%norb
!!!        write(*,'(i5,30(1pe15.8))')i,(smat(i,iorb),iorb=1,orbsig%norb)
!!!     end do
!!!
!!!     !overlap calculation of the kinetic operator
!!!     call kinetic_overlap(G,G,ovrlp)
!!!     call dsymm('L','U',G%ncoeff,orbsig%norb,1.0_gp,ovrlp(1),G%ncoeff,&
!!!          gaucoeff(1,1),G%ncoeff,0.d0,tmp(1,1),G%ncoeff)
!!!
!!!     call gemm('T','N',orbsig%norb,orbsig%norb,G%ncoeff,1.0_gp,&
!!!          gaucoeff(1,1),G%ncoeff,tmp(1,1),G%ncoeff,0.0_wp,smat(1,1),orbsig%norb)
!!!
!!!     !print overlap matrices
!!!     tt=0.0_wp
!!!     do i=1,orbsig%norb
!!!        write(*,'(i5,30(1pe15.8))')i,(smat(i,iorb),iorb=1,orbsig%norb)
!!!        !write(12,'(i5,30(1pe15.8))')i,(smat(i,iorb),iorb=1,orbsig%norb)
!!!        tt=tt+smat(i,i)
!!!     end do
!!!     print *,'trace',tt
!!!
!!!     !overlap calculation of the kinetic operator
!!!     call cpu_time(t0)
!!!     call potential_overlap(G,G,rhopot,Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,hxh,hyh,hzh,&
!!!          ovrlp)
!!!     call cpu_time(t1)
!!!     call dsymm('L','U',G%ncoeff,orbsig%norb,1.0_gp,ovrlp(1),G%ncoeff,&
!!!          gaucoeff(1,1),G%ncoeff,0.d0,tmp(1,1),G%ncoeff)
!!!
!!!     call gemm('T','N',orbsig%norb,orbsig%norb,G%ncoeff,1.0_gp,&
!!!          gaucoeff(1,1),G%ncoeff,tmp(1,1),G%ncoeff,0.0_wp,smat(1,1),orbsig%norb)
!!!
!!!     !print overlap matrices
!!!     tt=0.0_wp
!!!     do i=1,orbsig%norb
!!!        write(*,'(i5,30(1pe15.8))')i,(smat(i,iorb),iorb=1,orbsig%norb)
!!!        !write(12,'(i5,30(1pe15.8))')i,(smat(i,iorb),iorb=1,orbsig%norb)
!!!        tt=tt+smat(i,i)
!!!     end do
!!!     print *,'trace',tt
!!!     print *, 'time',t1-t0
!!!
!!!     i_all=-product(shape(ovrlp))*kind(ovrlp)
!!!     deallocate(ovrlp,stat=i_stat)
!!!     call memocc(i_stat,i_all,'ovrlp',subname)
!!!     i_all=-product(shape(tmp))*kind(tmp)
!!!     deallocate(tmp,stat=i_stat)
!!!     call memocc(i_stat,i_all,'tmp',subname)
!!!     i_all=-product(shape(smat))*kind(smat)
!!!     deallocate(smat,stat=i_stat)
!!!     call memocc(i_stat,i_all,'smat',subname)
!!!  end if

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

  !allocate the wavefunction in the transposed way to avoid allocations/deallocations
  allocate(hchi(orbsig%npsidim,at%nat),stat=i_stat)
  call memocc(i_stat,hchi,'hchi',subname)

  !call dcopy(orbsig%npsidim,psi,1,hpsi,1)
  if (input%exctxpar == 'OP2P') eexctX = -99.0_gp

  call full_local_potential(iproc,nproc,Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2),Glr%d%n1i*Glr%d%n2i*Glr%d%n3i,input%nspin,&
       orbsig%norb,orbsig%norbp,ngatherarr,rhopot,pot)




  !call HamiltonianApplication(iproc,nproc,at,orbsig,hx,hy,hz,rxyz,&
  !     nlpspd,proj,Glr,ngatherarr,pot,&
  !     psi,hpsi,ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU,pkernel=pkernelseq)

  !! ATTENTION:
  !! this works ONLY if the number of atomic orbitals used here is identical to the number of
  !! basis functions phi used later!
  !! HAS TO BE IMPROVED (adapt lin%onWhichAtom)

  !allocate(hpsiForAllLocregs(orbsig%npsidim,at%nat), stat=i_stat)
  !call memocc(i_stat,hpsiForAllLocregs,'hpsiForAllLocregs',subname)
  call HamiltonianApplicationConfinementForAllLocregs(iproc,nproc,at,orbsig,lin,input%hx,input%hy,input%hz,rxyz,&
       nlpspd,proj,Glr,ngatherarr,Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2),&
       rhopot(1),&
       chi,hchi,ekin_sum,epot_sum,eexctX,eproj_sum,input%nspin,GPU, rxyz, onWhichAtom, pkernel=pkernelseq)



  !deallocate potential
  call free_full_potential(nproc,pot,subname)

!!!  !calculate the overlap matrix knowing that the original functions are gaussian-based
!!!  allocate(thetaphi(2,G%nat+ndebug),stat=i_stat)
!!!  call memocc(i_stat,thetaphi,'thetaphi',subname)
!!!  thetaphi=0.0_gp
!!!
!!!  !calculate the scalar product between the hamiltonian and the gaussian basis
!!!  allocate(hpsigau(G%ncoeff,orbsig%norbp+ndebug),stat=i_stat)
!!!  call memocc(i_stat,hpsigau,'hpsigau',subname)
!!!
!!!
!!!  call wavelets_to_gaussians(at%geocode,orbsig%norbp,Glr%d%n1,Glr%d%n2,Glr%d%n3,G,&
!!!       thetaphi,hx,hy,hz,Glr%wfd,hpsi,hpsigau)
!!!
!!!  i_all=-product(shape(thetaphi))*kind(thetaphi)
!!!  deallocate(thetaphi,stat=i_stat)
!!!  call memocc(i_stat,i_all,'thetaphi',subname)

  accurex=abs(eks-ekin_sum)
  !tolerance for comparing the eigenvalues in the case of degeneracies
  etol=accurex/real(orbsig%norbu,gp)
  if (iproc == 0 .and. verbose > 1) write(*,'(1x,a,2(f19.10))') 'done. ekin_sum,eks:',ekin_sum,eks
  if (iproc == 0) then
     write(*,'(1x,a,3(1x,1pe18.11))') 'ekin_sum,epot_sum,eproj_sum',  & 
          ekin_sum,epot_sum,eproj_sum
     write(*,'(1x,a,3(1x,1pe18.11))') '   ehart,   eexcu,    vexcu',ehart,eexcu,vexcu
  endif

!!!  call Gaussian_DiagHam(iproc,nproc,at%natsc,nspin,orbs,G,mpirequests,&
!!!       psigau,hpsigau,orbsig,etol,norbsc_arr)


!!!  i_all=-product(shape(mpirequests))*kind(mpirequests)
!!!  deallocate(mpirequests,stat=i_stat)
!!!  call memocc(i_stat,i_all,'mpirequests',subname)

!!!  i_all=-product(shape(hpsigau))*kind(hpsigau)
!!!  deallocate(hpsigau,stat=i_stat)
!!!  call memocc(i_stat,i_all,'hpsigau',subname)

  !free GPU if it is the case
  if (GPUconv) then
     call free_gpu(GPU,orbsig%norbp)
  else if (OCLconv) then
     call free_gpu_OCL(GPU,orbsig,nspin_ig)
  end if

  if (iproc == 0 .and. verbose > 1) write(*,'(1x,a)')&
       'Input Wavefunctions Orthogonalization:'

  !psivirt can be eliminated here, since it will be allocated before davidson
  !with a gaussian basis
!!$  call DiagHam(iproc,nproc,at%natsc,nspin_ig,orbs,Glr%wfd,comms,&
!!$       psi,hpsi,psit,orbsig,commsig,etol,norbsc_arr,orbsv,psivirt)

  !call DiagHamForAllLocregs(iproc,nproc,at%natsc,nspin_ig,orbs,Glr%wfd,comms,&
  !     psi,hpsi,psit,input,at,orbsig,commsig,etol,norbsc_arr)
   call buildLinearCombinations(iproc, nproc, orbsig, lin%orbs, commsig, comms, at, Glr, lin%norbsPerType, &
              onWhichAtom, chi, hchi, phi)

  if (input%itrpmax > 1 .or. input%Tel > 0.0_gp) then
     !use the eval array of orbsig structure to save the original values
     allocate(orbsig%eval(lin%orbs%norb*lin%orbs%nkpts+ndebug),stat=i_stat)
     call memocc(i_stat,orbsig%eval,'orbsig%eval',subname)
     
     call dcopy(lin%orbs%norb*lin%orbs%nkpts,lin%orbs%eval(1),1,orbsig%eval(1),1)

     !add a small displacement in the eigenvalues
     do iorb=1,lin%orbs%norb*lin%orbs%nkpts
        tt=builtin_rand(idum)
        lin%orbs%eval(iorb)=lin%orbs%eval(iorb)*(1.0_gp+max(input%Tel,1.0e-3_gp)*real(tt,gp))
     end do

     !correct the occupation numbers wrt fermi level
     call evaltoocc(iproc,nproc,.false.,input%Tel,lin%orbs)

     !restore the occupation numbers
     call dcopy(lin%orbs%norb*lin%orbs%nkpts,orbsig%eval(1),1,lin%orbs%eval(1),1)

     i_all=-product(shape(orbsig%eval))*kind(orbsig%eval)
     deallocate(orbsig%eval,stat=i_stat)
     call memocc(i_stat,i_all,'orbsig%eval',subname)
  end if

  call deallocate_comms(commsig,subname)

  i_all=-product(shape(norbsc_arr))*kind(norbsc_arr)
  deallocate(norbsc_arr,stat=i_stat)
  call memocc(i_stat,i_all,'norbsc_arr',subname)

  if (iproc == 0) then
     !gaussian estimation valid only for Free BC
     if (at%geocode == 'F') then
        write(*,'(1x,a,1pe9.2)') 'expected accuracy in energy ',accurex
        write(*,'(1x,a,1pe9.2)') &
          'expected accuracy in energy per orbital ',accurex/real(lin%orbs%norb,kind=8)
        !write(*,'(1x,a,1pe9.2)') &
        !     'suggested value for gnrm_cv ',accurex/real(orbs%norb,kind=8)
     end if
  endif

  !here we can define the subroutine which generates the coefficients for the virtual orbitals
  call deallocate_gwf(G,subname)

  i_all=-product(shape(psigau))*kind(psigau)
  deallocate(psigau,stat=i_stat)
  call memocc(i_stat,i_all,'psigau',subname)

  call deallocate_orbs(orbsig,subname)


  i_all=-product(shape(onWhichAtom))*kind(onWhichAtom)
  deallocate(onWhichAtom, stat=i_stat)
  call memocc(i_stat, i_all, 'onWhichAtom',subname)
  



END SUBROUTINE inputguessConfinement





subroutine buildLinearCombinations(iproc, nproc, orbsig, orbs, commsig, comms, at, Glr, norbsPerType, &
           onWhichAtom, chi, hchi, phi)
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

! Local variables
integer:: iorb, jorb, iat, ist, jst, nvctrp, iall, istat, lwork, info, ierr
real(8),dimension(:),allocatable:: eval, work
real(8),dimension(:,:),allocatable:: ovrlp, ovrlpTemp
real(8),dimension(:,:,:),allocatable:: Ham
real(8),dimension(:),pointer:: psiw
real(8):: ddot
character(len=*),parameter:: subname='buildLinearCombinations'
     

  ! Transpose all the wavefunctions for having a piece of all the orbitals 
  ! for each processor
  allocate(psiw(orbsig%npsidim+ndebug),stat=istat)
  call memocc(istat, psiw, 'psiw', subname)
  call transpose_v(iproc, nproc, orbsig, Glr%wfd, commsig, chi, work=psiw)
  do iat=1,at%nat
      call transpose_v(iproc, nproc, orbsig, Glr%wfd, commsig, hchi(1,iat), work=psiw)
  end do
  iall=-product(shape(psiw))*kind(psiw)
  deallocate(psiw,stat=istat)
  call memocc(istat, iall, 'psiw', subname)

  ! Calculate Hamiltonian matrix and overlap matrix.
  allocate(Ham(orbsig%norb,orbsig%norb,at%nat), stat=istat)
  call memocc(istat, Ham, 'Ham', subname)
  allocate(ovrlp(orbsig%norb,orbsig%norb), stat=istat)
  call memocc(istat, ovrlp, 'ovrlp', subname)
  Ham=0.d0
  ovrlp=0.d0
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
  ist=1
  do iorb=1,orbsig%norb
      jst=1
      do jorb=1,orbsig%norb
          ovrlp(jorb,iorb)=ddot(nvctrp, chi(ist), 1, chi(jst), 1)
          jst=jst+nvctrp
      end do
      ist=ist+nvctrp
  end do
  call mpiallred(Ham(1,1,1), orbsig%norb**2*at%nat, mpi_sum, mpi_comm_world, ierr)
  call mpiallred(ovrlp(1,1), orbsig%norb**2, mpi_sum, mpi_comm_world, ierr)

  ! Diagonalize Hamiltonian matrix
  lwork=100*orbsig%norb
  allocate(work(lwork), stat=istat)
  call memocc(istat, work, 'work', subname)
  allocate(eval(orbsig%norb), stat=istat)
  call memocc(istat, eval, 'eval', subname)
  allocate(ovrlpTemp(orbsig%norb,orbsig%norb), stat=istat)
  call memocc(istat, ovrlpTemp, 'ovrlpTemp', subname)
  do iat=1,at%nat
      ! Copy ovrlp, since dsygv overwrites it.
      call dcopy(orbsig%norb**2, ovrlp(1,1), 1, ovrlpTemp(1,1), 1)
      call dsygv(1, 'v', 'u', orbsig%norb, Ham(1,1,iat), orbsig%norb, ovrlpTemp(1,1), orbsig%norb, &
           eval(1), work(1), lwork, info)
      if(info/=0) stop 'ERROR in dsygev: buildLinearCombinations'
  end do

  ! Build new linear combination
  phi=0.d0
  ist=1
  do iat=1,at%nat
      do iorb=1,norbsPerType(at%iatype(iat))
          jst=1
          do jorb=1,orbsig%norb
              call daxpy(nvctrp, Ham(jorb,iorb,iat), chi(jst), 1, phi(ist), 1)
              jst=jst+nvctrp
          end do
          ist=ist+nvctrp
      end do
  end do
  
  ! Untranpose
  allocate(psiw(orbs%npsidim+ndebug),stat=istat)
  call memocc(istat, psiw, 'psiw', subname)
  call transpose_v(iproc, nproc, orbs, Glr%wfd, comms, phi, work=psiw)
  iall=-product(shape(psiw))*kind(psiw)
  deallocate(psiw, stat=istat)
  call memocc(istat, iall, 'psiw', subname)

  iall=-product(shape(Ham))*kind(Ham)
  deallocate(Ham, stat=istat)
  call memocc(istat, iall, 'Ham', subname)
  
  iall=-product(shape(ovrlp))*kind(ovrlp)
  deallocate(ovrlp, stat=istat)
  call memocc(istat, iall, 'ovrlp', subname)

  iall=-product(shape(ovrlpTemp))*kind(ovrlpTemp)
  deallocate(ovrlpTemp, stat=istat)
  call memocc(istat, iall, 'ovrlpTemp', subname)

  iall=-product(shape(work))*kind(work)
  deallocate(work, stat=istat)
  call memocc(istat, iall, 'work', subname)

  iall=-product(shape(eval))*kind(eval)
  deallocate(eval, stat=istat)
  call memocc(istat, iall, 'eval', subname)



end subroutine buildLinearCombinations
