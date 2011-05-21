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
integer,dimension(:),allocatable:: onWhichAtom, onWhichAtomp, norbsPerAt, onWhichAtomTemp, onWhichAtomPhi
integer, parameter :: nmax=6,lmax=3,noccmax=2,nelecmax=32
integer, dimension(lmax+1) :: nl
real(gp), dimension(noccmax,lmax+1) :: occup
real(8):: dnrm2, ddot, dasum
integer:: ist, jst, jorb, iiAt, i
real(8),dimension(:,:),allocatable:: ttarr1, ttarr2



  allocate(norbsc_arr(at%natsc+1,input%nspin+ndebug),stat=i_stat)
  call memocc(i_stat,norbsc_arr,'norbsc_arr',subname)
  allocate(locrad(at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,locrad,'locrad',subname)
  allocate(norbsPerAt(at%nat), stat=istat)
  call memocc(istat, norbsPerAt, 'norbsPerAt', subname)


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
      !! WARNING: add more atomic orbitals
      if(at%iatype(iat)==1) then
          !at%aocc(7,iat)=1.d0
      else if(at%iatype(iat)==2) then
          !at%aocc(3,iat)=1.d0
      else
          stop 'only carbon and hydrogen ATM.'
      end if
      write(*,'(a,i6,100es10.2)') 'iat, at%aocc(:,iat)',iat, at%aocc(:,iat)
      call count_atomic_shells(lmax+1, noccmax, nelecmax, input%nspin, lin%orbs%nspinor, at%aocc(1,iat), occup, nl)
      !write(*,'(a,i4,2x,10i4)') 'iat, nl', iat, nl
      norbsPerAt(iat)=(nl(1)+3*nl(2)+5*nl(3)+7*nl(4))
      norbat=norbat+norbsPerAt(iat)
  end do
  write(*,*) 'norbat',norbat
  orbsig%norb=norbat ! This will be assigned later anyway, but here it is needed
                     ! to pass the array onWhichAtom to the subroutine inputguess_gaussian_orbitals

  allocate(onWhichAtom(norbat),stat=i_stat)
  call memocc(i_stat, onWhichAtom, 'onWhichAtom', subname)

  nvirt=0
  call inputguess_gaussian_orbitals_withOnWhichAtom(iproc,nproc,at,rxyz,Glr,nvirt,nspin_ig,&
       lin%orbs,orbsig,norbsc_arr,locrad,G,psigau,eks,onWhichAtom)
  !call inputguess_gaussian_orbitals(iproc,nproc,at,rxyz,Glr,nvirt,nspin_ig,&
  !     lin%orbs,orbsig,norbsc_arr,locrad,G,psigau,eks)
  !!write(*,*) 'orbsig%norb',orbsig%norb
  !!do iorb=1,orbsig%norb
  !!    write(*,'(a,3i8)') 'iproc, iorb, onWhichAtom(iorb)', iproc, iorb, onWhichAtom(iorb)
  !!end do

  !allocate communications arrays for inputguess orbitals
  !call allocate_comms(nproc,orbsig,commsig,subname)
  call orbitals_communicators(iproc,nproc,Glr,orbsig,commsig)  

  !write(*,*) '>>>>> iproc, orbsig%norbp', iproc, orbsig%norbp
  allocate(onWhichAtomp(orbsig%norbp),stat=i_stat)
  call memocc(i_stat, onWhichAtomp, 'onWhichAtomp', subname)

  call assignOrbitalsToAtoms(iproc, orbsig, at%nat, norbsPerAt, onWhichAtomp)
  do iorb=1,orbsig%norbp
      write(*,'(a,3i8)') 'iproc, iorb, onWhichAtom(iorb)', iproc, iorb, onWhichAtomp(iorb)
  end do


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
  chi=0.d0

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
  hchi=0.d0

  !call dcopy(orbsig%npsidim,psi,1,hpsi,1)
  if (input%exctxpar == 'OP2P') eexctX = -99.0_gp


  ! Orthogonalize the atomic basis functions (Loewdin)
  call orthonormalizeAtomicOrbitals(iproc, nproc, orbsig, commsig, Glr, chi)
  


  call full_local_potential(iproc,nproc,Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2),Glr%d%n1i*Glr%d%n2i*Glr%d%n3i,input%nspin,&
       orbsig%norb,orbsig%norbp,ngatherarr,rhopot,pot)




  !call HamiltonianApplication(iproc,nproc,at,orbsig,hx,hy,hz,rxyz,&
  !     nlpspd,proj,Glr,ngatherarr,pot,&
  !     psi,hpsi,ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU,pkernel=pkernelseq)

  !!!! THIS IS STUPID LIKE THIS, CHANGE LATER
  !!do iat=1,at%nat
  !!call HamiltonianApplication(iproc,nproc,at,orbsig,input%hx,input%hy,input%hz,rxyz,&
  !!     nlpspd,proj,Glr,ngatherarr,pot,&
  !!     chi,hchi(1,iat),ekin_sum,epot_sum,eexctX,eproj_sum,input%nspin,GPU,pkernel=pkernelseq)
  !!end do

  !! ATTENTION:
  !! this works ONLY if the number of atomic orbitals used here is identical to the number of
  !! basis functions phi used later!
  !! HAS TO BE IMPROVED (adapt lin%onWhichAtom)

  !allocate(hpsiForAllLocregs(orbsig%npsidim,at%nat), stat=i_stat)
  !call memocc(i_stat,hpsiForAllLocregs,'hpsiForAllLocregs',subname)
  !!call HamiltonianApplicationConfinementForAllLocregs(iproc,nproc,at,orbsig,lin,input%hx,input%hy,input%hz,rxyz,&
  !!     nlpspd,proj,Glr,ngatherarr,Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2),&
  !!     rhopot(1),&
  !!     chi,hchi,ekin_sum,epot_sum,eexctX,eproj_sum,input%nspin,GPU, rxyz, onWhichAtom, pkernel=pkernelseq)

  allocate(onWhichAtomTemp(orbsig%norbp), stat=istat)
  do iat=1,at%nat
      do iorb=1,orbsig%norbp
          iiAt=onWhichAtomp(iorb)
          tt = (rxyz(1,iiAt)-rxyz(1,iat))**2 + (rxyz(2,iiAt)-rxyz(2,iat))**2 + (rxyz(3,iiAt)-rxyz(3,iat))**2
          if(tt>50) then
              !onWhichAtomTemp(iorb)=-1
              onWhichAtomTemp(iorb)=iat
          else
              onWhichAtomTemp(iorb)=iat
              write(*,*) 'iat, iiAt', iat, iiAt
          end if
      end do
      if(iproc==0) write(*,'(a,i0)') 'iat=',iat
      call HamiltonianApplicationConfinement(iproc, nproc, at, orbsig, lin, input%hx, input%hy, input%hz, rxyz,&
           nlpspd, proj, Glr, ngatherarr, Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2), &
           rhopot(1), &
           chi(1), hchi(1,iat), ekin_sum, epot_sum, eexctX, eproj_sum, input%nspin, GPU, rxyz, onWhichAtomTemp, pkernel=pkernelseq)
      write(*,'(a,2i8,es15.6)') 'iproc, iat, dsum', iproc, iat, dasum((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*orbsig%norbp, hchi(1,iat), 1)
      do iorb=1,orbsig%norbp
          ist=0
          !if(onWhichAtomTemp(iorb)==-1) then
          !    write(*,'(a,3i7)') 'in main: set to zero. iat, iproc, iorb', iat, iproc, iorb
          !    call razero(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f, hchi(ist+1,iat))
          !end if
          !do i=1,Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
          !    write(10000+1000*iat+100*iorb+10*iproc,*) hchi(i+ist,iat)
          !end do
          ist=ist+Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
      end do
      call mpi_barrier(mpi_comm_world, iorb)
  end do



  !!if(iproc==4) then
  !!    allocate(ttarr1(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i,8))
  !!    allocate(ttarr2(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i,8))
  !!    open(unit=1,file='fort.2001')
  !!        do iorb=1,Glr%d%n1i*Glr%d%n2i*Glr%d%n3i
  !!            read(1,*) jorb, ttarr1(iorb,1)
  !!        end do
  !!    close(unit=1)
  !!    open(unit=1,file='fort.2002')
  !!        do iorb=1,Glr%d%n1i*Glr%d%n2i*Glr%d%n3i
  !!            read(1,*) jorb, ttarr1(iorb,2)
  !!        end do
  !!    close(unit=1)
  !!    open(unit=1,file='fort.2003')
  !!        do iorb=1,Glr%d%n1i*Glr%d%n2i*Glr%d%n3i
  !!            read(1,*) jorb, ttarr1(iorb,3)
  !!        end do
  !!    close(unit=1)
  !!    open(unit=1,file='fort.2004')
  !!        do iorb=1,Glr%d%n1i*Glr%d%n2i*Glr%d%n3i
  !!            read(1,*) jorb, ttarr1(iorb,4)
  !!        end do
  !!    close(unit=1)
  !!    open(unit=1,file='fort.2005')
  !!        do iorb=1,Glr%d%n1i*Glr%d%n2i*Glr%d%n3i
  !!            read(1,*) jorb, ttarr1(iorb,5)
  !!        end do
  !!    close(unit=1)
  !!    open(unit=1,file='fort.2006')
  !!        do iorb=1,Glr%d%n1i*Glr%d%n2i*Glr%d%n3i
  !!            read(1,*) jorb, ttarr1(iorb,6)
  !!        end do
  !!    close(unit=1)
  !!    open(unit=1,file='fort.2007')
  !!        do iorb=1,Glr%d%n1i*Glr%d%n2i*Glr%d%n3i
  !!            read(1,*) jorb, ttarr1(iorb,7)
  !!        end do
  !!    close(unit=1)
  !!    open(unit=1,file='fort.2008')
  !!        do iorb=1,Glr%d%n1i*Glr%d%n2i*Glr%d%n3i
  !!            read(1,*) jorb, ttarr1(iorb,8)
  !!        end do
  !!    close(unit=1)
  !!    open(unit=1,file='fort.2101')
  !!        do iorb=1,Glr%d%n1i*Glr%d%n2i*Glr%d%n3i
  !!            read(1,*) jorb, ttarr2(iorb,1)
  !!        end do
  !!    close(unit=1)
  !!    open(unit=1,file='fort.2102')
  !!        do iorb=1,Glr%d%n1i*Glr%d%n2i*Glr%d%n3i
  !!            read(1,*) jorb, ttarr2(iorb,2)
  !!        end do
  !!    close(unit=1)
  !!    open(unit=1,file='fort.2103')
  !!        do iorb=1,Glr%d%n1i*Glr%d%n2i*Glr%d%n3i
  !!            read(1,*) jorb, ttarr2(iorb,3)
  !!        end do
  !!    close(unit=1)
  !!    open(unit=1,file='fort.2104')
  !!        do iorb=1,Glr%d%n1i*Glr%d%n2i*Glr%d%n3i
  !!            read(1,*) jorb, ttarr2(iorb,4)
  !!        end do
  !!    close(unit=1)
  !!    open(unit=1,file='fort.2105')
  !!        do iorb=1,Glr%d%n1i*Glr%d%n2i*Glr%d%n3i
  !!            read(1,*) jorb, ttarr2(iorb,5)
  !!        end do
  !!    close(unit=1)
  !!    open(unit=1,file='fort.2106')
  !!        do iorb=1,Glr%d%n1i*Glr%d%n2i*Glr%d%n3i
  !!            read(1,*) jorb, ttarr2(iorb,6)
  !!        end do
  !!    close(unit=1)
  !!    open(unit=1,file='fort.2107')
  !!        do iorb=1,Glr%d%n1i*Glr%d%n2i*Glr%d%n3i
  !!            read(1,*) jorb, ttarr2(iorb,7)
  !!        end do
  !!    close(unit=1)
  !!    open(unit=1,file='fort.2108')
  !!        do iorb=1,Glr%d%n1i*Glr%d%n2i*Glr%d%n3i
  !!            read(1,*) jorb, ttarr2(iorb,8)
  !!        end do
  !!    close(unit=1)
  !!    do iorb=1,orbsig%norbp
  !!        do jorb=1,orbsig%norbp
  !!            write(*,'(a,2i7,3es15.5)') 'debug: iorb, jorb, <chi|chi>, <chi|hchi>, <hchi|hchi>', iorb, jorb, ddot(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i, ttarr1(1,iorb), 1, ttarr1(1,jorb), 1), &
  !!                                      ddot(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i, ttarr1(1,iorb), 1, ttarr2(1,jorb), 1), ddot(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i, ttarr2(1,iorb), 1, ttarr2(1,jorb), 1)
  !!        end do  
  !!    end do  
  !!end if

  
  !write(*,*) 'iproc, orbsig%norbp', iproc, orbsig%norbp
  !do iat=1,at%nat
  !    ist=1
  !    do iorb=1,orbsig%norbp
  !        jst=1
  !        do jorb=1,orbsig%norbp
  !            write(*,'(a,5i8,3es13.5)') 'iproc, iat, iorb, jorb, onWhichAtom(orbsig%isorb+iorb), <hchi|hchi>, <chi|hchi>, <chi|chi>', iproc, iat, iorb, jorb, onWhichAtom(orbsig%isorb+iorb), &
  !                                       ddot(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f, hchi(ist,iat), 1, hchi(jst,iat), 1), ddot(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f, chi(ist), 1, hchi(jst,iat), 1), &
  !                                       ddot(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f, chi(ist), 1, chi(jst), 1)
  !            jst=jst+Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
  !        end do
  !        ist=ist+Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
  !    end do
  !end do



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

   allocate(onWhichAtomPhi(lin%orbs%norb))
   onWhichAtomPhi=0
   do iorb=1,lin%orbs%norbp
       onWhichAtomPhi(lin%orbs%isorb+iorb)=lin%onWhichAtom(iorb)
   end do
   call mpiallred(onWhichAtomPhi(1), lin%orbs%norb, mpi_sum, mpi_comm_world, iorb)
   do iorb=1,lin%orbs%norb
       if(iproc==0) write(*,*) 'iorb, owap', iorb, onWhichAtomPhi(iorb)
   end do

   call buildLinearCombinations(iproc, nproc, orbsig, lin%orbs, commsig, lin%comms, at, Glr, lin%norbsPerType, &
              onWhichAtom, chi, hchi, phi, rxyz, onWhichAtomPhi, lin)

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
  
  i_all=-product(shape(norbsPerAt))*kind(norbsPerAt)
  deallocate(norbsPerAt, stat=i_stat)
  call memocc(i_stat, i_all, 'norbsPerAt',subname)



END SUBROUTINE inputguessConfinement





subroutine buildLinearCombinations(iproc, nproc, orbsig, orbs, commsig, comms, at, Glr, norbsPerType, &
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
integer:: iorb, jorb, korb, iat, ist, jst, nvctrp, iall, istat, lwork, info, ierr, infoCoeff, k, l,it, iiAt, jjAt
real(8),dimension(:),allocatable:: eval, work, alpha
real(8),dimension(:,:),allocatable:: ovrlp, ovrlpTemp, coeffTemp
real(8),dimension(:,:),allocatable:: coeff, grad, gradOld, lagMat, ovrlpLarge, coeffOld
real(8),dimension(:,:,:),allocatable:: Ham, tempArr
real(8),dimension(:),pointer:: chiw
real(8):: ddot, cosangle, ebsMod, ebsModOld, tt, dnrm2, fnrm, meanAlpha
logical:: converged
character(len=*),parameter:: subname='buildLinearCombinations'
real(4):: ttreal
     

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

  ! Calculate Hamiltonian matrix and overlap matrix.
  allocate(Ham(orbsig%norb,orbsig%norb,at%nat), stat=istat)
  call memocc(istat, Ham, 'Ham', subname)
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


  !! ATTENTION: modify Hamiltonian matrix
  do iat=1,at%nat
      do iorb=1,orbsig%norb
          do jorb=1,orbsig%norb
              !!!!!iiAt=onWhichAtom(iorb)
              !!!!jjAt=onWhichAtom(jorb)
              !!!!tt = (rxyz(1,iat)-rxyz(1,jjAt))**2 + (rxyz(2,iat)-rxyz(2,jjAt))**2 + (rxyz(3,iat)-rxyz(3,jjAt))**2
              !!!!tt=lin%potentialPrefac(at%iatype(iat))*tt**(lin%confPotOrder/2)
              !!!!!if(tt>50.d0) then
              !!!!!    if(iproc==0) write(*,'(a,2i7,es14.5)') 'set matrix to zero: iorb, jorb, tt',  iorb, jorb, tt
              !!!!!    Ham(iorb,jorb,iat)=0.d0
              !!!!    Ham(iorb,jorb,iat)=Ham(iorb,jorb,iat)+tt
              !!!!!end if
              !! NEW  !!!!!!!!!!
              !if(iorb/=jorb) Ham(iorb,jorb,iat)=0.d0
          end do
      end do
  end do   
  !!


  allocate(coeff(orbsig%norb,orbs%norb), stat=istat)
  call memocc(istat, coeff, 'coeff', subname)
  allocate(alpha(orbs%norb), stat=istat)
  call memocc(istat, alpha, 'alpha', subname)
  allocate(grad(orbsig%norb,orbs%norb), stat=istat)
  call memocc(istat, grad, 'grad', subname)
  allocate(gradOld(orbsig%norb,orbs%norb), stat=istat)
  call memocc(istat, gradOld, 'gradOld', subname)
  allocate(lagMat(orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, lagMat, 'lagMat', subname)

  allocate(coeffOld(orbsig%norb,orbs%norb), stat=istat)

        allocate(ovrlp(orbs%norb,orbs%norb))
        allocate(eval(orbs%norb))
        allocate(tempArr(orbs%norb,orbs%norb,2), stat=istat)
        call memocc(istat, tempArr, 'tempArr', subname)
        allocate(coeffTemp(orbsig%norb,orbs%norb), stat=istat)
        call memocc(istat, coeffTemp, 'coeffTemp', subname)

  processIf: if(iproc==0) then

      !call random_number(coeff)
      do iorb=1,orbs%norb
          iiAt=onWhichAtomPhi(iorb)
          do jorb=1,orbsig%norb
              jjAt=onWhichAtom(jorb)
              write(*,'(a,4i8)') 'iorb, jorb, onWhichAtomPhi(iorb), onWhichAtom(jorb)', iorb, jorb, onWhichAtomPhi(iorb), onWhichAtom(jorb)
              tt = (rxyz(1,iiat)-rxyz(1,jjAt))**2 + (rxyz(2,iiat)-rxyz(2,jjAt))**2 + (rxyz(3,iiat)-rxyz(3,jjAt))**2
              if(tt>50.d0) then
                  coeff(jorb,iorb)=0.d0
              else
                  call random_number(ttreal)
                  coeff(jorb,iorb)=dble(ttreal)
              end if
          end do
      end do

      do iat=1,at%nat
          do iorb=1,orbsig%norb
              do jorb=1,orbsig%norb
                   if(iproc==0) write(200+iat,*) iorb, jorb, Ham(iorb,jorb,iat)
              end do
          end do
      end do
      do iorb=1,orbs%norb
          do jorb=1,orbsig%norb
              write(998,*) iorb, jorb, coeff(jorb,iorb)
          end do
      end do

      !tt=0.d0
      !do iorb=1,orbsig%norb
      !    tt=tt+Ham(iorb,iorb)
      !end do
      !if(iproc==0) write(*,*) 'tr(Ham)',tt


    
    
    ! Initial step size for the optimization
    alpha=5.d-3

    ! Flag which checks convergence.
    converged=.false.

    if(iproc==0) write(*,'(x,a)') '============================== optmizing coefficients =============================='

    ! The optimization loop.
    iterLoop: do it=1,50000

        if (iproc==0) then
            write( *,'(1x,a,i0)') repeat('-',77 - int(log(real(it))/log(10.))) // ' iter=', it
        endif

      !! Put to zero the bad elements
      !do iorb=1,orbs%norb
      !    iiAt=onWhichAtomPhi(iorb)
      !    do jorb=1,orbsig%norb
      !        jjAt=onWhichAtom(jorb)
      !        tt = (rxyz(1,iiat)-rxyz(1,jjAt))**2 + (rxyz(2,iiat)-rxyz(2,jjAt))**2 + (rxyz(3,iiat)-rxyz(3,jjAt))**2
      !        if(tt>50.d0) then
      !            coeff(jorb,iorb)=0.d0
      !        end if
      !    end do
      !end do

 
        ! copy coeff to coeffOld (only for debugging)
        coeffOld=coeff


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
        do iorb=1,orbs%norb
            do jorb=1,orbs%norb
                ovrlp(iorb,jorb)=ddot(orbsig%norb, coeff(1,iorb), 1, coeff(1,jorb), 1)
            end do
        end do

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

        !!!! CHECK
        !!do iorb=1,orbs%norb
        !!    do jorb=1,orbs%norb
        !!        write(777,'(a,2i8,es15.5)') 'iorb, jorb, ddot', iorb, jorb, ddot(orbsig%norb, coeff(1,iorb), 1, coeff(1,jorb), 1)
        !!    end do
        !!end do




        !!write(888,*) '>>>>>>>>>>>>>>>>>it=',it
        !!do iorb=1,orbs%norb
        !!    do jorb=1,orbs%norb
        !!        tt=ddot(orbsig%norb, coeff(1,iorb), 1, coeff(1,jorb), 1)
        !!        write(888,*) iorb, jorb, tt
        !!    end do
        !!    do jorb=1,orbsig%norb
        !!        write(999,*) iorb, jorb, coeff(jorb,iorb)
        !!    end do
        !!end do


        ! Calculate the gradient grad. At the same time we determine whether the step size shall be increased
        ! or decreased (depending on gradient feedback).
        meanAlpha=0.d0
        grad=0.d0
        !if(ebsMod<ebsModOld+1.d-6*abs(ebsModOld)) then
        !    write(*,*) 'trace decreasing by', ebsMod-ebsModOld
        !else
        !    write(*,*) 'trace increasing by', ebsMod-ebsModOld
        !end if
        do iorb=1,orbs%norb
            !iiAt=onWhichAtom(iorb)
            iiAt=onWhichAtomPhi(iorb)
            !!do l=1,orbsig%norb
            !!    do k=1,orbsig%norb
            !!        !grad(l,iorb)=grad(l,iorb)+coeff(k,iorb)*(Ham(k,l,iiAt)+Ham(l,k,iiAt))
            !!        grad(l,iorb)=grad(l,iorb)+Ham(l,k,iiAt)*coeff(k,iorb)
            !!    end do
            !!end do
            call dgemv('n', orbsig%norb, orbsig%norb, 1.d0, Ham(1,1,iiAt), orbsig%norb, coeff(1,iorb), 1, 0.d0, grad(1,iorb), 1)
            if(it>1) then
                cosangle=ddot(orbsig%norb, grad(1,iorb), 1, gradOld(1,iorb), 1)
                cosangle=cosangle/dnrm2(orbsig%norb, grad(1,iorb), 1)
                cosangle=cosangle/dnrm2(orbsig%norb, gradOld(1,iorb), 1)
                !write(*,'(a,es14.5,2es17.8)') 'cosangle, ebsMod, ebsModOld', cosangle, ebsMod, ebsModOld
                if(cosangle>.8d0 .and. ebsMod<ebsModOld+1.d-6*abs(ebsModOld)) then
                    alpha(iorb)=max(alpha(iorb)*1.05d0,1.d-6)
                else
                    alpha(iorb)=max(alpha(iorb)*.5d0,1.d-6)
                end if
            end if
            call dcopy(orbsig%norb, grad(1,iorb), 1, gradOld(1,iorb), 1)
            meanAlpha=meanAlpha+alpha(iorb)
        end do
        meanAlpha=meanAlpha/orbs%norb
        !!do iorb=1,orbs%norb
        !!    do jorb=1,orbsig%norb
        !!        write(1111,'(2i7,es16.6,6x,2es16.6,4x,es12.3)') iorb, jorb, grad(jorb,iorb), coeff(jorb,iorb), coeffOld(jorb,iorb), alpha(iorb)
        !!    end do
        !!end do

        !do iorb=1,orbs%norb
        !    do k=1,orbsig%norb
        !        write(400,*) iorb, k, grad(k,iorb)
        !    end do
        !end do
    
    
        ! Apply the orthoconstraint to the gradient. To do so first calculate the Lagrange
        ! multiplier matrix.
        lagMat=0.d0
        if(it>1) then
            ebsModOld=ebsMod
        else
            ebsModOld=1.d10
        end if
        ebsMod=0.d0
        do iorb=1,orbs%norb
            do jorb=1,orbs%norb
                lagMat(jorb,iorb)=ddot(orbsig%norb, coeff(1,jorb), 1, grad(1,iorb), 1)
                !do k=1,orbsig%norb
                !    lagMat(iorb,jorb)=lagMat(iorb,jorb)+coeff(k,iorb)*grad(k,jorb)
                !end do
                !write(500,*) iorb, jorb, lagMat(jorb,iorb)
            end do
            ebsMod=ebsMod+lagMat(iorb,iorb)
        end do

        ! Now apply the orthoconstraint.
        do iorb=1,orbs%norb
            do k=1,orbsig%norb
                do jorb=1,orbs%norb
                    grad(k,iorb)=grad(k,iorb)-.5d0*(lagMat(iorb,jorb)*coeff(k,jorb)+lagMat(jorb,iorb)*coeff(k,jorb))
                    !if(iorb==jorb) grad(k,iorb)=grad(k,iorb)-.5d0*(lagMat(iorb,jorb)*coeff(k,jorb)+lagMat(jorb,iorb)*coeff(k,jorb))
                end do
            end do
        end do
        !!do iorb=1,orbs%norb
        !!    do jorb=1,orbsig%norb
        !!        write(1112,'(2i7,es16.6)') iorb, jorb, grad(jorb,iorb)
        !!    end do
        !!end do

        !do iorb=1,orbs%norb
        !    do k=1,orbsig%norb
        !        write(600,*) iorb, k, grad(k,iorb)
        !    end do
        !end do
    
        
        !!$! Calculate the modified band structure energy and the gradient norm.
        !!$if(it>1) then
        !!$    ebsModOld=ebsMod
        !!$else
        !!$    ebsModOld=1.d10
        !!$end if
        !ebsMod=0.d0
        fnrm=0.d0
        do iorb=1,orbs%norb
            !iiAt=onWhichAtom(iorb)
            fnrm=fnrm+dnrm2(orbsig%norb, grad(1,iorb), 1)
            !do jorb=1,orbsig%norb
            !    do korb=1,orbsig%norb
            !        ebsMod=ebsMod+coeff(korb,iorb)*coeff(jorb,iorb)*Ham(korb,jorb,iiAt)
            !    end do
            !end do
        end do
    
        ! Multiply the energy with a factor of 2 if we have a closed-shell system.
        !if(nspin==1) then
        !    ebsMod=2.d0*ebsMod
        !end if

        !if(iproc==0) write(*,'(x,a,4x,i0,es12.4,3x,es10.3, es19.9)') 'iter, fnrm, meanAlpha, Energy', &
        if(iproc==0) write(*,'(x,a,es11.2,es22.13,es10.2)') 'fnrm, band structure energy, mean alpha', &
            fnrm, ebsMod, meanAlpha
        
        ! Check for convergence.
        if(fnrm<1.d-5) then
            if(iproc==0) write(*,'(x,a,i0,a)') 'converged in ', it, ' iterations.'
            if(iproc==0) write(*,'(3x,a,2es14.5)') 'Final values for fnrm, Energy:', fnrm, ebsMod
            converged=.true.
            infoCoeff=it
            exit
        end if
  
        if(it==50000) then
            if(iproc==0) write(*,'(x,a,i0,a)') 'WARNING: not converged within ', it, &
                ' iterations! Exiting loop due to limitations of iterations.'
            if(iproc==0) write(*,'(x,a,2es15.7,f12.7)') 'Final values for fnrm, Energy: ', fnrm, ebsMod
            infoCoeff=-1
            exit
        end if

        ! Improve the coefficients (by steepet descent).
        do iorb=1,orbs%norb
            do l=1,orbsig%norb
                coeff(l,iorb)=coeff(l,iorb)-alpha(iorb)*grad(l,iorb)
            end do
        end do
    

    end do iterLoop


    if(iproc==0) write(*,'(x,a)') '===================================================================================='

end if processIf


! Now broadcast the result to all processes
call mpi_bcast(coeff(1,1), orbsig%norb*orbs%norb, mpi_double_precision, 0, mpi_comm_world, ierr)
call mpi_bcast(infoCoeff, 1, mpi_integer, 0, mpi_comm_world, ierr)
call mpi_barrier(mpi_comm_world, ierr)
write(*,*) 'after MPI stuff'







!#############################################

  !!ist=1
  !!do iorb=1,orbsig%norb
  !!    jst=1
  !!    do jorb=1,orbsig%norb
  !!        Ham(jorb,iorb)=ddot(nvctrp, chi(ist), 1, chi(jst), 1)
  !!        jst=jst+nvctrp
  !!    end do
  !!    ist=ist+nvctrp
  !!end do
  !!call mpiallred(Ham(1,1), orbsig%norb**2, mpi_sum, mpi_comm_world, ierr)
  !!do iorb=1,orbsig%norb
  !!    do jorb=1,orbsig%norb
  !!        write(300,*) iorb, jorb, Ham(iorb,jorb)
  !!    end do
  !!end do

  ! Build new linear combination
  phi=0.d0
  ist=1
  do iorb=1,orbs%norb
      jst=1
      do jorb=1,orbsig%norb
          !write(*,'(a,3i7,es15.5)') 'iproc, iorb, jorb, Ham(jorb,iorb,iat)', iproc, iorb, jorb, Ham(jorb,iorb,iat)
          call daxpy(nvctrp, coeff(jorb,iorb), chi(jst), 1, phi(ist), 1)
          jst=jst+nvctrp
      end do
      ist=ist+nvctrp
  end do


  ! Check orthogonality
  ist=1
  do iorb=1,orbs%norb
      jst=1
      do jorb=1,orbs%norb
          ovrlp(jorb,iorb)=ddot(nvctrp, phi(ist), 1, phi(jst), 1)
          jst=jst+nvctrp
      end do
      ist=ist+nvctrp
  end do
  call mpiallred(ovrlp(1,1), orbs%norb**2, mpi_sum, mpi_comm_world, ierr)
  do iorb=1,orbs%norb
      do jorb=1,orbs%norb
          if(iproc==0) write(*,'(a,2i7,es15.5)') 'iorb, jorb, ovrlp(iorb,jorb)', iorb, jorb, ovrlp(iorb,jorb)
      end do
  end do
  

  ! Untranpose
  allocate(chiw(orbs%npsidim+ndebug),stat=istat)
  call memocc(istat, chiw, 'chiw', subname)
  call untranspose_v(iproc, nproc, orbs, Glr%wfd, comms, phi, work=chiw)
  iall=-product(shape(chiw))*kind(chiw)
  deallocate(chiw, stat=istat)
  call memocc(istat, iall, 'chiw', subname)

  iall=-product(shape(Ham))*kind(Ham)
  deallocate(Ham, stat=istat)
  call memocc(istat, iall, 'Ham', subname)

  iall=-product(shape(coeff))*kind(coeff)
  deallocate(coeff, stat=istat)
  call memocc(istat, iall, 'coeff', subname)
  
  iall=-product(shape(ovrlp))*kind(ovrlp)
  deallocate(ovrlp, stat=istat)
  call memocc(istat, iall, 'ovrlp', subname)



end subroutine buildLinearCombinations


subroutine orthonormalizeAtomicOrbitals(iproc, nproc, orbsig, commsig, Glr, chi)

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
