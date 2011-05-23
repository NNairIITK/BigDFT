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
real(8),dimension(:),allocatable:: eval
integer:: istat
real(8),dimension(:),allocatable:: chi
real(8),dimension(:,:),allocatable:: hchi
integer,dimension(:),allocatable:: onWhichAtom, onWhichAtomp, norbsPerAt, onWhichAtomTemp, onWhichAtomPhi
integer, parameter :: nmax=6,lmax=3,noccmax=2,nelecmax=32
integer, dimension(lmax+1) :: nl
real(gp), dimension(noccmax,lmax+1) :: occup
real(8):: dnrm2, ddot, dasum
integer:: ist, jst, jorb, iiAt, i, iadd, ii, jj



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

  ! Determine how many atomic orbitals we have.
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

  allocate(onWhichAtom(norbat),stat=i_stat)
  call memocc(i_stat, onWhichAtom, 'onWhichAtom', subname)
  ist=0
  do iat=1,at%nat
      ! Assign the orbitals to the atoms, i.e. onWhichAtom(i)=j means that orbital i belongs
      ! to atom j.
      do i=1,norbsPerAt(iat)
          onWhichAtom(ist+i)=iat
      end do
      ist=ist+norbsPerAt(iat)
  end do


  nvirt=0
  call inputguess_gaussian_orbitals(iproc,nproc,at,rxyz,Glr,nvirt,nspin_ig,&
       lin%orbs,orbsig,norbsc_arr,locrad,G,psigau,eks)

  !allocate communications arrays for inputguess orbitals
  !call allocate_comms(nproc,orbsig,commsig,subname)
  call orbitals_communicators(iproc,nproc,Glr,orbsig,commsig)  

  allocate(onWhichAtomp(orbsig%norbp),stat=i_stat)
  call memocc(i_stat, onWhichAtomp, 'onWhichAtomp', subname)

  call assignOrbitalsToAtoms(iproc, orbsig, at%nat, norbsPerAt, onWhichAtomp)
  !do iorb=1,orbsig%norbp
  !    write(*,'(a,3i8)') 'iproc, iorb, onWhichAtom(iorb)', iproc, iorb, onWhichAtomp(iorb)
  !end do


  hxh=.5_gp*input%hx
  hyh=.5_gp*input%hy
  hzh=.5_gp*input%hz

  !check the communication distribution
  !call check_communications(iproc,nproc,orbsig,Glr,commsig)

  !once the wavefunction coefficients are known perform a set 
  !of nonblocking send-receive operations to calculate overlap matrices


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

  ! Allocate the atomic orbitals.
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

  
  ! Create the potential

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
  

  ! This is not required at teh moment since HamiltonianApplicationConfinement has a different structure.
  call full_local_potential(iproc,nproc,Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2),Glr%d%n1i*Glr%d%n2i*Glr%d%n3i,input%nspin,&
       orbsig%norb,orbsig%norbp,ngatherarr,rhopot,pot)

  allocate(onWhichAtomTemp(orbsig%norbp), stat=istat)
  call memocc(i_stat,onWhichAtomTemp,'onWhichAtomTemp',subname)
  if(iproc==0) write(*,'(x,a)') 'Hamiltonian application for all atoms. This may take some time.'
  do iat=1,at%nat
      do iorb=1,orbsig%norbp
          onWhichAtomTemp(iorb)=iat
      end do
      if(iproc==0) write(*,'(3x,a,i0,a)', advance='no') 'Hamiltonian application for atom ', iat, '... '
      call HamiltonianApplicationConfinement(iproc, nproc, at, orbsig, lin, input%hx, input%hy, input%hz, rxyz,&
           nlpspd, proj, Glr, ngatherarr, Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2), &
           rhopot(1), &
           chi(1), hchi(1,iat), ekin_sum, epot_sum, eexctX, eproj_sum, input%nspin, GPU, rxyz, onWhichAtomTemp, pkernel=pkernelseq)
      if(iproc==0) write(*,'(a)') 'done.'
  end do
  !!call HamiltonianApplicationConfinementForAllLocregs(iproc, nproc, at, orbsig, lin, input%hx, input%hy, input%hz, rxyz,&
  !!     nlpspd, proj, Glr, ngatherarr, Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2), rhopot, chi, hchi, &
  !!     ekin_sum, epot_sum, eexctX, eproj_sum, input%nspin, GPU, rxyz, onWhichAtom, pkernel=pkernelseq)

  !deallocate potential
  call free_full_potential(nproc,pot,subname)


  ! This does not make sense here...
  !!accurex=abs(eks-ekin_sum)
  !!!tolerance for comparing the eigenvalues in the case of degeneracies
  !!etol=accurex/real(orbsig%norbu,gp)
  !!if (iproc == 0 .and. verbose > 1) write(*,'(1x,a,2(f19.10))') 'done. ekin_sum,eks:',ekin_sum,eks
  !!if (iproc == 0) then
  !!   write(*,'(1x,a,3(1x,1pe18.11))') 'ekin_sum,epot_sum,eproj_sum',  & 
  !!        ekin_sum,epot_sum,eproj_sum
  !!   write(*,'(1x,a,3(1x,1pe18.11))') '   ehart,   eexcu,    vexcu',ehart,eexcu,vexcu
  !!endif


  !free GPU if it is the case
  if (GPUconv) then
     call free_gpu(GPU,orbsig%norbp)
  else if (OCLconv) then
     call free_gpu_OCL(GPU,orbsig,nspin_ig)
  end if



   allocate(onWhichAtomPhi(lin%orbs%norb))
   call memocc(i_stat,onWhichAtomPhi,'onWhichAtomPhi',subname)
   onWhichAtomPhi=0
   do iorb=1,lin%orbs%norbp
       onWhichAtomPhi(lin%orbs%isorb+iorb)=lin%onWhichAtom(iorb)
   end do
   call mpiallred(onWhichAtomPhi(1), lin%orbs%norb, mpi_sum, mpi_comm_world, iorb)
   !do iorb=1,lin%orbs%norb
   !    if(iproc==0) write(*,*) 'iorb, owap', iorb, onWhichAtomPhi(iorb)
   !end do

   call buildLinearCombinations(iproc, nproc, orbsig, lin%orbs, commsig, lin%comms, at, Glr, lin%norbsPerType, &
        onWhichAtom, chi, hchi, phi, rxyz, onWhichAtomPhi, lin)

  if(iproc==0) write(*,'(x,a)') '------------------------------------------------------------- Input guess generated.'

  ! This is probably not needed..
  !!if (input%itrpmax > 1 .or. input%Tel > 0.0_gp) then
  !!   !use the eval array of orbsig structure to save the original values
  !!   allocate(orbsig%eval(lin%orbs%norb*lin%orbs%nkpts+ndebug),stat=i_stat)
  !!   call memocc(i_stat,orbsig%eval,'orbsig%eval',subname)
  !!   
  !!   call dcopy(lin%orbs%norb*lin%orbs%nkpts,lin%orbs%eval(1),1,orbsig%eval(1),1)

  !!   !add a small displacement in the eigenvalues
  !!   do iorb=1,lin%orbs%norb*lin%orbs%nkpts
  !!      tt=builtin_rand(idum)
  !!      lin%orbs%eval(iorb)=lin%orbs%eval(iorb)*(1.0_gp+max(input%Tel,1.0e-3_gp)*real(tt,gp))
  !!   end do

  !!   !correct the occupation numbers wrt fermi level
  !!   call evaltoocc(iproc,nproc,.false.,input%Tel,lin%orbs)

  !!   !restore the occupation numbers
  !!   call dcopy(lin%orbs%norb*lin%orbs%nkpts,orbsig%eval(1),1,lin%orbs%eval(1),1)

  !!   i_all=-product(shape(orbsig%eval))*kind(orbsig%eval)
  !!   deallocate(orbsig%eval,stat=i_stat)
  !!   call memocc(i_stat,i_all,'orbsig%eval',subname)
  !!end if

  call deallocate_comms(commsig,subname)

  i_all=-product(shape(norbsc_arr))*kind(norbsc_arr)
  deallocate(norbsc_arr,stat=i_stat)
  call memocc(i_stat,i_all,'norbsc_arr',subname)

  ! This does not make sense here...
  !!if (iproc == 0) then
  !!   !gaussian estimation valid only for Free BC
  !!   if (at%geocode == 'F') then
  !!      write(*,'(1x,a,1pe9.2)') 'expected accuracy in energy ',accurex
  !!      write(*,'(1x,a,1pe9.2)') &
  !!        'expected accuracy in energy per orbital ',accurex/real(lin%orbs%norb,kind=8)
  !!      !write(*,'(1x,a,1pe9.2)') &
  !!      !     'suggested value for gnrm_cv ',accurex/real(orbs%norb,kind=8)
  !!   end if
  !!endif

  !here we can define the subroutine which generates the coefficients for the virtual orbitals
  call deallocate_gwf(G,subname)

  i_all=-product(shape(psigau))*kind(psigau)
  deallocate(psigau,stat=i_stat)
  call memocc(i_stat,i_all,'psigau',subname)

  call deallocate_orbs(orbsig,subname)


  i_all=-product(shape(onWhichAtom))*kind(onWhichAtom)
  deallocate(onWhichAtom, stat=i_stat)
  call memocc(i_stat, i_all, 'onWhichAtom',subname)

  i_all=-product(shape(onWhichAtomPhi))*kind(onWhichAtomPhi)
  deallocate(onWhichAtomPhi, stat=i_stat)
  call memocc(i_stat, i_all, 'onWhichAtomPhi',subname)

  i_all=-product(shape(onWhichAtomp))*kind(onWhichAtomp)
  deallocate(onWhichAtomp, stat=i_stat)
  call memocc(i_stat, i_all, 'onWhichAtomp',subname)

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
real(8):: ddot, cosangle, tt, dnrm2, fnrm, meanAlpha, cut, trace, traceOld
logical:: converged
character(len=*),parameter:: subname='buildLinearCombinations'
real(4):: ttreal
     

  if(iproc==0) write(*,'(x,a)') '------------------------------- Minimizing trace in the basis of the atomic orbitals'

  ! Allocate the local arrays
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
  allocate(Ham(orbsig%norb,orbsig%norb,at%nat), stat=istat)
  call memocc(istat, Ham, 'Ham', subname)

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





  processIf: if(iproc==0) then

      ! Initialize the coefficient vectors. Put random number to palces where it is
      ! reasonable (i.e. close to the atom where the basis function is centered).
      do iorb=1,orbs%norb
          iiAt=onWhichAtomPhi(iorb)
          cut=1.d0/lin%potentialPrefac(at%iatype(iiAt))
          cut=cut**(1.d0/dble(lin%confPotOrder))
          do jorb=1,orbsig%norb
              jjAt=onWhichAtom(jorb)
              tt = (rxyz(1,iiat)-rxyz(1,jjAt))**2 + (rxyz(2,iiat)-rxyz(2,jjAt))**2 + (rxyz(3,iiat)-rxyz(3,jjAt))**2
              if(tt>cut) then
                  coeff(jorb,iorb)=0.d0
              else
                  call random_number(ttreal)
                  coeff(jorb,iorb)=dble(ttreal)
              end if
          end do
      end do


    
    ! Initial step size for the optimization
    alpha=5.d-3

    ! Flag which checks convergence.
    converged=.false.

    if(iproc==0) write(*,'(x,a)') '============================== optmizing coefficients =============================='

    ! The optimization loop.
    iterLoop: do it=1,lin%nItInguess

        if (iproc==0 .and. mod(it,1000)==1) then
            write( *,'(1x,a,i0)') repeat('-',77 - int(log(real(it))/log(10.))) // ' iter=', it
        endif


        ! Othonormalize the coefficients.
        call orthonormalizeCoefficients(orbs, orbsig, coeff)



        ! Calculate the gradient grad. At the same time we determine whether the step size shall be increased
        ! or decreased (depending on gradient feedback).
        meanAlpha=0.d0
        grad=0.d0
        do iorb=1,orbs%norb
            iiAt=onWhichAtomPhi(iorb)
            call dgemv('n', orbsig%norb, orbsig%norb, 1.d0, Ham(1,1,iiAt), orbsig%norb, coeff(1,iorb), 1, 0.d0, grad(1,iorb), 1)
            if(it>1) then
                cosangle=ddot(orbsig%norb, grad(1,iorb), 1, gradOld(1,iorb), 1)
                cosangle=cosangle/dnrm2(orbsig%norb, grad(1,iorb), 1)
                cosangle=cosangle/dnrm2(orbsig%norb, gradOld(1,iorb), 1)
                if(cosangle>.8d0 .and. trace<traceOld+1.d-6*abs(traceOld)) then
                    alpha(iorb)=max(alpha(iorb)*1.05d0,1.d-6)
                else
                    alpha(iorb)=max(alpha(iorb)*.5d0,1.d-6)
                end if
            end if
            call dcopy(orbsig%norb, grad(1,iorb), 1, gradOld(1,iorb), 1)
            meanAlpha=meanAlpha+alpha(iorb)
        end do
        meanAlpha=meanAlpha/orbs%norb
    
    
        ! Apply the orthoconstraint to the gradient. To do so first calculate the Lagrange
        ! multiplier matrix.
        lagMat=0.d0
        if(it>1) then
            traceOld=trace
        else
            traceOld=1.d10
        end if
        trace=0.d0
        call gemm('t', 'n', orbs%norb, orbs%norb, orbsig%norb, 1.0d0, coeff(1,1), &
             orbsig%norb, grad(1,1), orbsig%norb, 0.d0, lagMat(1,1), orbs%norb)
        do iorb=1,orbs%norb
            trace=trace+lagMat(iorb,iorb)
        end do


        ! Now apply the orthoconstraint.
        call dgemm('n', 'n', orbsig%norb, orbs%norb, orbs%norb, -.5d0, coeff(1,1), orbsig%norb, &
             lagMat(1,1), orbs%norb, 1.d0, grad(1,1), orbsig%norb)
        call dgemm('n', 't', orbsig%norb, orbs%norb, orbs%norb, -.5d0, coeff(1,1), orbsig%norb, &
             lagMat(1,1), orbs%norb, 1.d0, grad(1,1), orbsig%norb)


        ! Calculate the gradient norm.
        fnrm=0.d0
        do iorb=1,orbs%norb
            fnrm=fnrm+dnrm2(orbsig%norb, grad(1,iorb), 1)
        end do
    

        ! Write some informations to the screen, but only every 1000th iteration.
        if(iproc==0 .and. mod(it,1000)==1) write(*,'(x,a,es11.2,es22.13,es10.2)') 'fnrm, trace, mean alpha', &
            fnrm, trace, meanAlpha
        
        ! Check for convergence.
        if(fnrm<1.d-5) then
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
            call daxpy(orbsig%norb, -alpha(iorb), grad(1,iorb), 1, coeff(1,iorb), 1)
        end do
    

    end do iterLoop


    if(iproc==0) write(*,'(x,a)') '===================================================================================='

end if processIf


! Now broadcast the result to all processes
call mpi_bcast(coeff(1,1), orbsig%norb*orbs%norb, mpi_double_precision, 0, mpi_comm_world, ierr)
call mpi_bcast(infoCoeff, 1, mpi_integer, 0, mpi_comm_world, ierr)
call mpi_barrier(mpi_comm_world, ierr)



  ! Build new linear combination
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

  iall=-product(shape(alpha))*kind(alpha)
  deallocate(alpha, stat=istat)
  call memocc(istat, iall, 'alpha', subname)

  iall=-product(shape(grad))*kind(grad)
  deallocate(grad, stat=istat)
  call memocc(istat, iall, 'grad', subname)

  iall=-product(shape(gradOld))*kind(gradOld)
  deallocate(gradOld, stat=istat)
  call memocc(istat, iall, 'gradOld', subname)

  iall=-product(shape(lagMat))*kind(lagMat)
  deallocate(lagMat, stat=istat)
  call memocc(istat, iall, 'lagMat', subname)

  


end subroutine buildLinearCombinations


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
