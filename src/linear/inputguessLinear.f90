!>   input guess wavefunction diagonalization
subroutine inputguessConfinement(iproc, nproc, at, &
     comms, Glr, input, lin, rxyz, n3p, rhopot, rhocore, pot_ion,&
     nlpspd, proj, pkernel, pkernelseq, &
     nscatterarr, ngatherarr, potshortcut, irrzon, phnons, GPU, radii_cf,  &
     lphi, ehart, eexcu, vexcu)
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
  real(8),dimension(lin%orbs%npsidim),intent(out):: lphi
  real(8),intent(out):: ehart, eexcu, vexcu
  !local variables
  type(gaussian_basis):: G !basis for davidson IG
  character(len=*), parameter :: subname='inputguessConfinement'
  logical :: switchGPUconv,switchOCLconv
  integer :: i_stat,i_all,iat,nspin_ig,iorb,idum=0, nvirt, norbat
  real(kind=4) :: tt,builtin_rand
  real(gp) :: hxh,hyh,hzh,eks,epot_sum,ekin_sum,eexctX,eproj_sum,etol,accurex
  integer, dimension(:,:), allocatable :: norbsc_arr
  real(wp), dimension(:), allocatable :: potxc
  real(gp), dimension(:), allocatable :: locrad
  real(wp), dimension(:), pointer :: pot
  real(wp), dimension(:,:,:), pointer :: psigau
type(linear_zone_descriptors):: lzdig, lzdGauss
type(p2pCommsGatherPot):: comgp
integer:: istat, tag
real(8),dimension(:),allocatable:: lchi, lchi2
real(8),dimension(:,:),allocatable::  lhchi
real(8),dimensioN(:,:,:),allocatable:: ham
integer,dimension(:),allocatable:: norbsPerAt, onWhichAtomTemp
integer, parameter :: nmax=6,lmax=3,noccmax=2,nelecmax=32
logical:: withConfinement, ovrlpx, ovrlpy, ovrlpz
logical,dimension(:),allocatable:: doNotCalculate
integer, dimension(lmax+1) :: nl
real(gp), dimension(noccmax,lmax+1) :: occup
real(8):: dnrm2, ddot, dasum, t1, t2, time
integer:: ist, jst, jorb, iiAt, i, iadd, ii, jj, ndimpot, ilr, ind1, ind2, ldim, gdim, ierr, jlr
integer:: is1, ie1, is2, ie2, is3, ie3, js1, je1, js2, je2, js3, je3


  if (iproc == 0) then
     write(*,'(1x,a)')&
          '------------------------------------------------------- Input Wavefunctions Creation'
  end if

  call nullify_linear_zone_descriptors(lzdig)
  call nullify_linear_zone_descriptors(lzdGauss)


  ! Allocate some arrays we need for the input guess.
  allocate(norbsc_arr(at%natsc+1,input%nspin+ndebug),stat=i_stat)
  call memocc(i_stat,norbsc_arr,'norbsc_arr',subname)
  allocate(locrad(at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,locrad,'locrad',subname)
  allocate(norbsPerAt(at%nat), stat=istat)
  call memocc(istat, norbsPerAt, 'norbsPerAt', subname)

  ! Number of localization regions
  lzdig%nlr=at%nat
  lzdGauss%nlr=at%nat



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

  ! Create the atomic orbitals in a Gaussian basis.
  nvirt=0
  call inputguess_gaussian_orbitals(iproc,nproc,at,rxyz,Glr,nvirt,nspin_ig,&
       lin%orbs,lzdig%orbs,norbsc_arr,locrad,G,psigau,eks)

  ! lzdig%orbs%inWhichLocreg has been allocated in inputguess_gaussian_orbitals. SInce it will again be allcoated
  ! in assignToLocreg2, deallocate it first.
  i_all=-product(shape(lzdig%orbs%inWhichLocreg))*kind(lzdig%orbs%inWhichLocreg)
  deallocate(lzdig%orbs%inWhichLocreg,stat=i_stat)
  call memocc(i_stat,i_all,'lzdig%orbs%inWhichLocreg',subname)

  call assignToLocreg2(iproc, at%nat, lzdig%nlr, input%nspin, norbsPerAt, lzdig%orbs)
  !lzdGauss%orbs=lzdig%orbs
  call copy_orbitals_data(lzdig%orbs, lzdGauss%orbs, subname)
  !lzdGauss%Glr=Glr
  call copy_locreg_descriptors(Glr, lzdGauss%Glr, subname)


  locrad=20.d0
  call initLocregs2(iproc, at%nat, rxyz, lzdGauss, input, Glr, locrad)
  call initLocregs2(iproc, at%nat, rxyz, lzdig, input, Glr, lin%locrad)
allocate(lchi(lzdig%orbs%npsidim+ndebug),stat=i_stat)
call memocc(i_stat,lchi,'lchi',subname)
allocate(lhchi(lzdig%orbs%npsidim,at%nat),stat=i_stat)
call memocc(i_stat,lhchi,'lhchi',subname)
allocate(lchi2(lzdGauss%orbs%npsidim),stat=i_stat)
call memocc(i_stat,lchi2,'lchi2',subname)
  lchi=0.d0
  lhchi=0.d0
  lchi2=0.d0


  hxh=.5_gp*input%hx
  hyh=.5_gp*input%hy
  hzh=.5_gp*input%hz


  lzdig%lpsidimtot=lzdig%orbs%npsidim
  lzdGauss%lpsidimtot=lzdGauss%orbs%npsidim

  
  lchi2=0.d0
  call gaussians_to_wavelets_new2(iproc, nproc, lzdGauss, input%hx, input%hy, input%hz, G, &
       psigau(1,1,min(lzdGauss%orbs%isorb+1, lzdGauss%orbs%norb)), lchi2(1))
  call orthonormalizeAtomicOrbitalsLocalized(iproc, nproc, lzdGauss, input, lchi2)
  ! Transform chi to the localization region - should be used if locrad for the input guess is larger than our localization radius.
  ind1=1
  ind2=1
  do iorb=1,lzdGauss%orbs%norbp
      ilr = lzdig%orbs%inWhichLocregp(iorb)
      ldim=lzdig%Llr(ilr)%wfd%nvctr_c+7*lzdig%Llr(ilr)%wfd%nvctr_f
      gdim=lzdGauss%llr(ilr)%wfd%nvctr_c+7*lzdGauss%llr(ilr)%wfd%nvctr_f
      call psi_to_locreg2(iproc, nproc, ldim, gdim, lzdig%llr(ilr), lzdGauss%llr(ilr), lchi2(ind1), lchi(ind2))
      ind1=ind1+lzdGauss%llr(ilr)%wfd%nvctr_c+7*lzdGauss%llr(ilr)%wfd%nvctr_f
      ind2=ind2+lzdig%Llr(ilr)%wfd%nvctr_c+7*lzdig%Llr(ilr)%wfd%nvctr_f
  end do
  if(ind1/=lzdGauss%orbs%npsidim+1) then
      write(*,'(2(a,i0))') 'ERROR on process ',iproc,': ind1/=lzdGauss%orbs%npsidim+1',ind1,lzdGauss%orbs%npsidim+1
      stop
  end if
  if(ind2/=lzdig%orbs%npsidim+1) then
      write(*,'(2(a,i0))') 'ERROR on process ',iproc,': ind2/=lzdig%orbs%npsidim+1',ind2,lzdig%orbs%npsidim+1
      stop
  end if


i_all=-product(shape(locrad))*kind(locrad)
deallocate(locrad,stat=i_stat)
call memocc(i_stat,i_all,'locrad',subname)

  
  !!! Create the potential.
  if(iproc==0) write(*,'(x,a)',advance='no') 'Calculating charge density...'
  call sumrhoLinear(iproc, nproc, lzdGauss, input%ixc, hxh, hyh, hzh, lchi2, rhopot,&
    & lzdGauss%Glr%d%n1i*lzdGauss%Glr%d%n2i*nscatterarr(iproc,1), nscatterarr, input%nspin, GPU, &
    & at%symObj, irrzon, phnons)
  if(iproc==0) write(*,'(a)') 'done.'

     
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


!!!! WHAT DOES THIS MEAN??
!!  if(potshortcut>0) then
!!!!$    if (GPUconv) then
!!!!$       call free_gpu(GPU,orbs%norbp)
!!!!$    end if
!!     if (switchGPUconv) then
!!        GPUconv=.true.
!!     end if
!!     if (switchOCLconv) then
!!        OCLconv=.true.
!!     end if
!!
!!     
!!     !deallocate the gaussian basis descriptors
!!     call deallocate_gwf(G,subname)
!!    
!!     i_all=-product(shape(psigau))*kind(psigau)
!!     deallocate(psigau,stat=i_stat)
!!     call memocc(i_stat,i_all,'psigau',subname)
!!     call deallocate_comms(commsig,subname)
!!     i_all=-product(shape(norbsc_arr))*kind(norbsc_arr)
!!     deallocate(norbsc_arr,stat=i_stat)
!!     call memocc(i_stat,i_all,'norbsc_arr',subname)
!!    return 
!!  end if


  !call dcopy(orbsig%npsidim,psi,1,hpsi,1)
  if (input%exctxpar == 'OP2P') eexctX = -99.0_gp


  ! Copy Glr to lzd
  !lzdig%Glr = Glr
  call copy_locreg_descriptors(Glr, lzdig%Glr, subname)
  
  ! Copy nlpspd to lin%lzd
  !lzdig%Gnlpspd = nlpspd
  call copy_nonlocal_psp_descriptors(nlpspd, lzdig%Gnlpspd, subname)
  
  ! Set localnorb
  do ilr=1,lzdig%nlr
      lzdig%Llr(ilr)%localnorb=0
      do iorb=1,lzdig%orbs%norbp
          !if(onWhichAtomp(iorb)==ilr) then
          if(lzdig%orbs%inWhichLocregp(iorb)==ilr) then
              lzdig%Llr(ilr)%localnorb = lzdig%Llr(ilr)%localnorb+1
          end if
      end do
  end do


  ! Initialize the parameters for the communications of the potential.
  tag=20000 !! CHANGE LATER!!
  call initializeCommunicationPotential(iproc, nproc, nscatterarr, lzdig%orbs, lzdig, comgp, lzdig%orbs%inWhichLocreg, tag)

  ! Post the messages for the communication of the potential.
  ndimpot = lin%lzd%Glr%d%n1i*lin%lzd%Glr%d%n2i*nscatterarr(iproc,2)
  call allocateCommunicationsBuffersPotential(comgp, subname)
  call postCommunicationsPotential(iproc, nproc, ndimpot, rhopot, comgp)

  
  ! Gather the potential
  call gatherPotential(iproc, nproc, comgp)

  ! Apply the Hamiltonian for each atom.
  ! onWhichAtomTemp indicates indicating that all orbitals feel the potential from atom iat.
allocate(onWhichAtomTemp(lzdig%orbs%norbp), stat=istat)
call memocc(i_stat,onWhichAtomTemp,'onWhichAtomTemp',subname)
allocate(doNotCalculate(lzdig%nlr), stat=istat)
call memocc(i_stat, doNotCalculate, 'doNotCalculate', subname)
  if(iproc==0) write(*,'(x,a)') 'Hamiltonian application for all atoms. This may take some time.'
  lhchi=0.d0
  call cpu_time(t1)
  withConfinement=.true.
  do iat=1,at%nat
      doNotCalculate=.false.
      call mpi_barrier(mpi_comm_world, ierr)
      call getIndices(lzdig%llr(iat), is1, ie1, is2, ie2, is3, ie3)
      do jorb=1,lzdig%orbs%norbp
          onWhichAtomTemp(jorb)=iat
          !jlr=onWhichAtomp(jorb)
          jlr=lzdig%orbs%inWhichLocreg(jorb)
          call getIndices(lzdig%llr(jlr), js1, je1, js2, je2, js3, je3)
          ovrlpx = ( is1<=je1 .and. ie1>=js1 )
          ovrlpy = ( is2<=je2 .and. ie2>=js2 )
          ovrlpz = ( is3<=je3 .and. ie3>=js3 )
          if(ovrlpx .and. ovrlpy .and. ovrlpz) then
          else
              !doNotCalculate(ilr)=.true.
              doNotCalculate(jlr)=.true.
          end if
      end do
      !write(*,'(a,2i4,4x,100l4)') 'iat, iproc, doNotCalculate', iat, iproc, doNotCalculate
      if(iproc==0) write(*,'(3x,a,i0,a)', advance='no') 'Hamiltonian application for atom ', iat, '... '
      call HamiltonianApplicationConfinement2(input, iproc, nproc, at, lzdig, lzdig%orbs, lin, input%hx, input%hy, input%hz, rxyz,&
           ngatherarr, comgp%nrecvBuf, comgp%recvBuf, lchi, lhchi(1,iat), &
           ekin_sum, epot_sum, eexctX, eproj_sum, input%nspin, GPU, radii_cf, comgp, onWhichAtomTemp,&
           withConfinement, pkernel=pkernelseq)

      if(iproc==0) write(*,'(a)') 'done.'
  end do
  call cpu_time(t2)
  time=t2-t1
  call mpiallred(time, 1, mpi_sum, mpi_comm_world, ierr)
  if(iproc==0) write(*,'(x,a,es10.3)') 'time for applying potential:', time/dble(nproc)

   call deallocateCommunicationsBuffersPotential(comgp, subname)
   call deallocate_p2pCommsGatherPot(comgp, subname)


   call cpu_time(t1)
 allocate(ham(lzdig%orbs%norb,lzdig%orbs%norb,at%nat), stat=istat)
 call memocc(i_stat,ham,'ham',subname)
   call getHamiltonianMatrix2(iproc, nproc, lzdig, Glr, input, lzdig%orbs%inWhichLocreg, lzdig%orbs%inWhichLocregp, at%nat, lchi, lhchi, ham)
   call buildLinearCombinationsLocalized(iproc, nproc, lzdig%orbs, lin%orbs, lin%comms, at, Glr, input, lin%norbsPerType, &
        lzdig%orbs%inWhichLocreg, lchi, lphi, rxyz, lin%orbs%inWhichLocreg, lin, lzdig, ham)
   !!do istat=1,size(lphi)
   call cpu_time(t2)
   time=t2-t1
   call mpiallred(time, 1, mpi_sum, mpi_comm_world, ierr)
   if(iproc==0) write(*,'(x,a,es10.3)') 'time for "buildLinearCombinations":', time/dble(nproc)

  if(iproc==0) write(*,'(x,a)') '------------------------------------------------------------- Input guess generated.'


  ! Deallocate all local arrays.
  call deallocate_gwf(G,subname)

  call deallocate_linear_zone_descriptors(lzdig, subname)
  call deallocate_linear_zone_descriptors(lzdGauss, subname)

  i_all=-product(shape(psigau))*kind(psigau)
  deallocate(psigau,stat=i_stat)
  call memocc(i_stat,i_all,'psigau',subname)

  i_all=-product(shape(norbsc_arr))*kind(norbsc_arr)
  deallocate(norbsc_arr,stat=i_stat)
  call memocc(i_stat,i_all,'norbsc_arr',subname)

  i_all=-product(shape(onWhichAtomTemp))*kind(onWhichAtomTemp)
  deallocate(onWhichAtomTemp, stat=i_stat)
  call memocc(i_stat, i_all, 'onWhichAtomTemp',subname)
  
  i_all=-product(shape(norbsPerAt))*kind(norbsPerAt)
  deallocate(norbsPerAt, stat=i_stat)
  call memocc(i_stat, i_all, 'norbsPerAt',subname)

  i_all=-product(shape(lchi))*kind(lchi)
  deallocate(lchi, stat=i_stat)
  call memocc(i_stat, i_all, 'lchi',subname)

  i_all=-product(shape(lhchi))*kind(lhchi)
  deallocate(lhchi, stat=i_stat)
  call memocc(i_stat, i_all, 'lhchi',subname)

  i_all=-product(shape(lchi2))*kind(lchi2)
  deallocate(lchi2, stat=i_stat)
  call memocc(i_stat, i_all, 'lchi2',subname)

  i_all=-product(shape(doNotCalculate))*kind(doNotCalculate)
  deallocate(doNotCalculate, stat=i_stat)
  call memocc(i_stat, i_all, 'doNotCalculate',subname)

  i_all=-product(shape(ham))*kind(ham)
  deallocate(ham, stat=i_stat)
  call memocc(i_stat, i_all, 'ham',subname)


END SUBROUTINE inputguessConfinement





subroutine buildLinearCombinations2(iproc, nproc, orbsig, orbs, commsig, comms, at, Glr, norbsPerType, &
           onWhichAtom, chi, hchi, phi, rxyz, onWhichAtomPhi, lin, ham)
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
real(8),dimension(orbsig%norb,orbsig%norb,at%nat),intent(inout):: ham

! Local variables
integer:: iorb, jorb, korb, iat, ist, jst, nvctrp, iall, istat, ierr, infoCoeff, k, l,it, iiAt, jjAt
real(8),dimension(:),allocatable:: alpha, coeffPad, coeff2, gradTemp, gradOld, fnrmArr, fnrmOvrlpArr, fnrmOldArr, grad
real(8),dimension(:,:),allocatable:: ovrlp, ovrlpTemp
real(8),dimension(:,:),allocatable:: coeff, lagMat, coeffOld
real(8),dimension(:,:,:),allocatable:: HamPad
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
          iiAt=onWhichAtomPhi(ip%isorb+iorb)
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
  nvctrp=commsig%nvctr_par(iproc,1) ! 1 for k-points
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





subroutine orthonormalizeAtomicOrbitalsLocalized(iproc, nproc, lzd, input, lchi)
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
type(linear_zone_descriptors),intent(in):: lzd
type(input_variables),intent(in):: input
real(8),dimension(lzd%orbs%npsidim),intent(inout):: lchi

! Local variables
integer:: iorb, jorb, istat, iall, lwork, info, nvctrp, ierr, tag, ilr
real(8),dimension(:,:),allocatable:: ovrlp
character(len=*),parameter:: subname='orthonormalizeAtomicOrbitalsLocalized'
type(overlapParameters):: op
type(p2pCommsOrthonormality):: comon


! Initialize the communication parameters.
tag=5000
call initCommsOrtho(iproc, nproc, lzd, lzd%orbs%inWhichLocreg, input, op, comon, tag)
allocate(ovrlp(lzd%orbs%norb,lzd%orbs%norb), stat=istat)
call memocc(istat, ovrlp, 'ovrlp', subname)

call orthonormalizeLocalized(iproc, nproc, 2, lzd%orbs, op, comon, lzd, lzd%orbs%inWhichLocreg, 1.d-6, input, lchi, ovrlp)

iall=-product(shape(ovrlp))*kind(ovrlp)
deallocate(ovrlp, stat=istat)
call memocc(istat, iall, 'ovrlp', subname)

call deallocate_overlapParameters(op, subname)
call deallocate_p2pCommsOrthonormality(comon, subname)

end subroutine orthonormalizeAtomicOrbitalsLocalized



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
integer:: ii, kk, jproc, istat, ierr, norbTarget, iorb, iiorb
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

  ! ip%onWhichMPI indicates on which MPI process a given orbital is located.
  allocate(ip%onWhichMPI(ip%norbtot), stat=istat)
  call memocc(istat, ip%onWhichMPI, 'ip%onWhichMPI', subname)
  iiorb=0
  do jproc=0,ip%nproc-1
      do iorb=1,ip%norb_par(jproc)
          iiorb=iiorb+1
          ip%onWhichMPI(iiorb)=jproc
      end do
  end do
  if(iproc==0) write(*,'(a,100i5)') 'owmpi', ip%onWhichMPI


  ! Starting orbital for each process
  allocate(ip%isorb_par(0:ip%nproc-1), stat=istat)
  call memocc(istat, ip%isorb_par, 'ip%isorb_par', subname)
  ip%isorb_par=0
  ii=0
  do jproc=0,iproc-1
     ii=ii+ip%norb_par(jproc)
  end do
  !reference orbital for process
  ip%isorb=ii
  ip%isorb_par(iproc)=ip%isorb
  call mpiallred(ip%isorb_par(0), ip%nproc, mpi_sum, newComm, ierr)
  if(iproc==0) write(*,'(a,100i5)') 'isorb_par', ip%isorb_par
  


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

  ! With the padded zeros, the elements can be distributed evenly.
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

        do iorb=1,ip%norb
            do jorb=1,ip%norb
                if(abs(ovrlp(iorb,jorb)-ovrlp(jorb,iorb))>1.d-10) stop 'not symmetric'
                write(3000+iproc,'(2i6,es26.17)') iorb, jorb, ovrlp(iorb,jorb)
            end do
        end do
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




subroutine getHamiltonianMatrix(iproc, nproc, lzdig, Glr, input, onWhichAtom, onWhichAtomp, nat, chi, hchi, ham, orbsig)
use module_base
use module_types
use module_interfaces, exceptThisOne => getHamiltonianMatrix
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, nat
type(linear_zone_descriptors),intent(in):: lzdig
type(locreg_descriptors),intent(in):: Glr
type(input_variables),intent(in):: input
integer,dimension(lzdig%orbs%norb),intent(in):: onWhichAtom
integer,dimension(lzdig%orbs%norbp),intent(in):: onWhichAtomp
type(orbitals_data),intent(in):: orbsig
!real(8),dimension(lzdig%orbs%npsidim),intent(in):: chi
!real(8),dimension(lzdig%orbs%npsidim,nat),intent(in):: hchi
real(8),dimension(orbsig%npsidim),intent(in):: chi
real(8),dimension(orbsig%npsidim,nat),intent(in):: hchi
real(8),dimension(lzdig%orbs%norb,lzdig%orbs%norb,nat),intent(out):: ham

! Local variables
integer:: sizeChi, istat, iorb, ilr, iall, ind1, ind2, ldim, gdim, iat, tag
type(overlapParameters):: op
type(p2pCommsOrthonormality):: comon
real(8),dimension(:),allocatable:: lchi, lhchi, lphiovrlp
character(len=*),parameter:: subname='getHamiltonianMatrix'


!! CHANGE THIS LATER?
tag=10000
! Initialize the parameters for calculating the matrix.
call initCommsOrtho(iproc, nproc, lzdig, onWhichAtom, input, op, comon, tag)

!allocate(lphiovrlp(op%ndim_lphiovrlp), stat=istat)
!call memocc(istat, lphiovrlp, 'lphiovrlp',subname)



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


call allocateCommuncationBuffersOrtho(comon, subname)
if(iproc==0) write(*,'(x,a)') 'Calculating Hamiltonian matrix for all atoms. This may take some time.'
do iat=1,nat
    if(iproc==0) write(*,'(3x,a,i0,a)', advance='no') 'Calculating matrix for atom ', iat, '... '
    ! Transform chi to the localization region. This is not needed if we really habe O(N).
    ind1=1
    ind2=1
    do iorb=1,lzdig%orbs%norbp
        ilr = onWhichAtomp(iorb)
        ldim=lzdig%Llr(ilr)%wfd%nvctr_c+7*lzdig%Llr(ilr)%wfd%nvctr_f
        gdim=Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
        call psi_to_locreg2(iproc, nproc, ldim, gdim, lzdig%Llr(ilr), Glr, chi(ind1), lchi(ind2))
        call psi_to_locreg2(iproc, nproc, ldim, gdim, lzdig%Llr(ilr), Glr, hchi(ind1,iat), lhchi(ind2))
        ind1=ind1+Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
        ind2=ind2+lzdig%Llr(ilr)%wfd%nvctr_c+7*lzdig%Llr(ilr)%wfd%nvctr_f
    end do

    ! Put lphi in the sendbuffer, i.e. lphi will be sent to other processes' receive buffer.
    !call extractOrbital(iproc, nproc, lzdig%orbs, sizeChi, onWhichAtom, lzdig, op, lchi(1), comon)
    call extractOrbital2(iproc, nproc, lzdig%orbs, sizeChi, onWhichAtom, lzdig, op, lchi(1), comon)
    call postCommsOverlap(iproc, nproc, comon)
    !call gatherOrbitals(iproc, nproc, comon)
    call gatherOrbitals2(iproc, nproc, comon)
    ! Put lhphi to the sendbuffer, so we can the calculate <lphi|lhphi>
    !call extractOrbital(iproc, nproc, lzdig%orbs, sizeChi, onWhichAtom, lzdig, op, lhchi(1), comon)
    call extractOrbital2(iproc, nproc, lzdig%orbs, sizeChi, onWhichAtom, lzdig, op, lhchi(1), comon)
    call calculateOverlapMatrix2(iproc, nproc, lzdig%orbs, op, comon, onWhichAtom, ham(1,1,iat))
    if(iproc==0) write(*,'(a)') 'done.'
end do
call deallocateCommuncationBuffersOrtho(comon, subname)
call deallocate_overlapParameters(op, subname)
call deallocate_p2pCommsOrthonormality(comon, subname)


iall=-product(shape(lchi))*kind(lchi)
deallocate(lchi, stat=istat)
call memocc(istat, iall, 'lchi', subname)

iall=-product(shape(lhchi))*kind(lhchi)
deallocate(lhchi, stat=istat)
call memocc(istat, iall, 'lhchi', subname)


end subroutine getHamiltonianMatrix


subroutine getHamiltonianMatrix2(iproc, nproc, lzdig, Glr, input, onWhichAtom, onWhichAtomp, nat, lchi, lhchi, ham)
use module_base
use module_types
use module_interfaces, exceptThisOne => getHamiltonianMatrix2
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, nat
type(linear_zone_descriptors),intent(in):: lzdig
type(locreg_descriptors),intent(in):: Glr
type(input_variables),intent(in):: input
integer,dimension(lzdig%orbs%norb),intent(in):: onWhichAtom
integer,dimension(lzdig%orbs%norbp),intent(in):: onWhichAtomp
!real(8),dimension(lzdig%orbs%npsidim),intent(in):: chi
!real(8),dimension(lzdig%orbs%npsidim,nat),intent(in):: hchi
real(8),dimension(lzdig%orbs%npsidim),intent(in):: lchi
real(8),dimension(lzdig%orbs%npsidim,nat),intent(in):: lhchi
real(8),dimension(lzdig%orbs%norb,lzdig%orbs%norb,nat),intent(out):: ham

! Local variables
integer:: sizeChi, istat, iorb, ilr, iall, ind1, ind2, ldim, gdim, iat, tag
type(overlapParameters):: op
type(p2pCommsOrthonormality):: comon
!real(8),dimension(:),allocatable:: lchi, lhchi, lphiovrlp
real(8),dimension(:),allocatable:: lphiovrlp
character(len=*),parameter:: subname='getHamiltonianMatrix'


!! CHANGE THIS LATER?
tag=10000
! Initialize the parameters for calculating the matrix.
call initCommsOrtho(iproc, nproc, lzdig, onWhichAtom, input, op, comon, tag)

!allocate(lphiovrlp(op%ndim_lphiovrlp), stat=istat)
!call memocc(istat, lphiovrlp, 'lphiovrlp',subname)



!!! Calculate the dimension of the wave function for each process.
!!sizeChi=0
!!do iorb=1,lzdig%orbs%norbp
!!    ilr=onWhichAtomp(iorb)
!!    sizeChi = sizeChi + (lzdig%Llr(ilr)%wfd%nvctr_c+7*lzdig%Llr(ilr)%wfd%nvctr_f)
!!end do
!!
!!allocate(lchi(sizeChi), stat=istat)
!!call memocc(istat, lchi, 'lchi', subname)
!!allocate(lhchi(sizeChi), stat=istat)
!!call memocc(istat, lhchi, 'lhchi', subname)


call allocateCommuncationBuffersOrtho(comon, subname)
if(iproc==0) write(*,'(x,a)') 'Calculating Hamiltonian matrix for all atoms. This may take some time.'
do iat=1,nat
    if(iproc==0) write(*,'(3x,a,i0,a)', advance='no') 'Calculating matrix for atom ', iat, '... '
    !!! Transform chi to the localization region. This is not needed if we really habe O(N).
    !!ind1=1
    !!ind2=1
    !!do iorb=1,lzdig%orbs%norbp
    !!    ilr = onWhichAtomp(iorb)
    !!    ldim=lzdig%Llr(ilr)%wfd%nvctr_c+7*lzdig%Llr(ilr)%wfd%nvctr_f
    !!    gdim=Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
    !!    call psi_to_locreg2(iproc, nproc, ldim, gdim, lzdig%Llr(ilr), Glr, chi(ind1), lchi(ind2))
    !!    call psi_to_locreg2(iproc, nproc, ldim, gdim, lzdig%Llr(ilr), Glr, hchi(ind1,iat), lhchi(ind2))
    !!    ind1=ind1+Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
    !!    ind2=ind2+lzdig%Llr(ilr)%wfd%nvctr_c+7*lzdig%Llr(ilr)%wfd%nvctr_f
    !!end do

    ! Put lphi in the sendbuffer, i.e. lphi will be sent to other processes' receive buffer.
    !call extractOrbital(iproc, nproc, lzdig%orbs, sizeChi, onWhichAtom, lzdig, op, lchi(1), comon)
    call extractOrbital2(iproc, nproc, lzdig%orbs, lzdig%orbs%npsidim, onWhichAtom, lzdig, op, lchi(1), comon)
    call postCommsOverlap(iproc, nproc, comon)
    !call gatherOrbitals(iproc, nproc, comon)
    call gatherOrbitals2(iproc, nproc, comon)
    ! Put lhphi to the sendbuffer, so we can the calculate <lphi|lhphi>
    !call extractOrbital(iproc, nproc, lzdig%orbs, sizeChi, onWhichAtom, lzdig, op, lhchi(1), comon)
    call extractOrbital2(iproc, nproc, lzdig%orbs, lzdig%orbs%npsidim, onWhichAtom, lzdig, op, lhchi(1,iat), comon)
    call calculateOverlapMatrix2(iproc, nproc, lzdig%orbs, op, comon, onWhichAtom, ham(1,1,iat))
    if(iproc==0) write(*,'(a)') 'done.'
end do
call deallocateCommuncationBuffersOrtho(comon, subname)

call deallocate_overlapParameters(op, subname)
call deallocate_p2pCommsOrthonormality(comon, subname)

!!iall=-product(shape(lchi))*kind(lchi)
!!deallocate(lchi, stat=istat)
!!call memocc(istat, iall, 'lchi', subname)
!!
!!iall=-product(shape(lhchi))*kind(lhchi)
!!deallocate(lhchi, stat=istat)
!!call memocc(istat, iall, 'lhchi', subname)


end subroutine getHamiltonianMatrix2




subroutine buildLinearCombinationsLocalized(iproc, nproc, orbsig, orbs, comms, at, Glr, input, norbsPerType, &
           onWhichAtom, lchi, lphi, rxyz, onWhichAtomPhi, lin, lzdig, ham)
!
use module_base
use module_types
use module_interfaces, exceptThisOne => buildLinearCombinationsLocalized
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data),intent(in):: orbsig, orbs
type(communications_arrays),intent(in):: comms
type(atoms_data),intent(in):: at
type(locreg_descriptors),intent(in):: Glr
type(input_variables),intent(in):: input
type(linearParameters),intent(in):: lin
type(linear_zone_descriptors),intent(inout):: lzdig
integer,dimension(at%ntypes):: norbsPerType
integer,dimension(orbsig%norb),intent(in):: onWhichAtom
real(8),dimension(lzdig%orbs%npsidim):: lchi
real(8),dimension(lin%orbs%npsidim):: lphi
real(8),dimension(3,at%nat):: rxyz
integer,dimension(orbs%norb):: onWhichAtomPhi
real(8),dimension(orbsig%norb,orbsig%norb,at%nat),intent(inout):: ham

! Local variables
integer:: iorb, jorb, korb, iat, ist, jst, nvctrp, iall, istat, ierr, infoCoeff, k, l,it, iiAt, jjAt
real(8),dimension(:),allocatable:: alpha, coeffPad, coeff2, gradTemp, gradOld, fnrmArr, fnrmOvrlpArr, fnrmOldArr, grad
real(8),dimension(:,:),allocatable:: ovrlp, ovrlpTemp
real(8),dimension(:,:),allocatable:: coeff, lagMat, coeffOld, lcoeff, lgrad, lgradold
real(8),dimension(:,:,:),allocatable:: HamPad
real(8),dimension(:),pointer:: chiw
integer,dimension(:),allocatable:: recvcounts, displs, norb_par
real(8):: ddot, cosangle, tt, dnrm2, fnrm, meanAlpha, cut, trace, traceOld, fnrmMax
logical:: converged
character(len=*),parameter:: subname='buildLinearCombinationsLocalized'
real(4):: ttreal
integer:: wholeGroup, newGroup, newComm, norbtot
integer,dimension(:),allocatable:: newID
  
! new
real(8),dimension(:),allocatable:: work, eval, evals
real(8),dimension(:,:),allocatable:: tempMat
integer:: lwork, ii, info, iiAtprev, i, jproc, norbTarget, sendcount, ilr, iilr, tag
type(inguessParameters):: ip
real(8),dimension(:,:,:),pointer:: hamextract
type(p2pCommsOrthonormalityMatrix):: comom
type(matrixMinimization):: matmin

  if(iproc==0) write(*,'(x,a)') '------------------------------- Minimizing trace in the basis of the atomic orbitals'

  ! Allocate the local arrays that are hold by all processes.
  allocate(coeff(orbsig%norb,orbs%norb), stat=istat)
  call memocc(istat, coeff, 'coeff', subname)

  call nullify_matrixMinimization(matmin)


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

      call determineLocalizationRegions(iproc, ip%nproc, lzdig%nlr, lzdig%orbs%norb, at, onWhichAtom, lin%locrad, rxyz, lin%lzd, matmin%mlr)
      call extractMatrix(iproc, ip%nproc, lin%orbs%norb, ip%norb_par(iproc), lzdig%orbs, onWhichAtomPhi, ip%onWhichMPI, at%nat, ham, matmin, hamextract)
      call determineOverlapRegionMatrix(iproc, ip%nproc, lin%lzd, matmin%mlr, lin%orbs, lzdig%orbs, onWhichAtom, onWhichAtomPhi, comom)
      tag=1
      call initCommsMatrixOrtho(iproc, ip%nproc, lin%orbs%norb, ip%norb_par, ip%isorb_par, onWhichAtomPhi, ip%onWhichMPI, tag, comom)
      allocate(lcoeff(matmin%norbmax,ip%norb_par(iproc)), stat=istat)
      call memocc(istat, lcoeff, 'lcoeff', subname)
      allocate(lgrad(matmin%norbmax,ip%norb_par(iproc)), stat=istat)
      call memocc(istat, lgrad, 'lgrad', subname)
      allocate(lgradold(matmin%norbmax,ip%norb_par(iproc)), stat=istat)
      call memocc(istat, lgradold, 'lgradold', subname)
    
      ! Initialize the coefficient vectors. Put random number to palces where it is
      ! reasonable (i.e. close to the atom where the basis function is centered).
      call initRandomSeed(iproc, 1)
      do ii=1,orbsig%isorb*ip%norbtotPad
          call random_number(ttreal)
      end do
      coeffPad=0.d0
    
      ii=0
      do iorb=1,ip%norb_par(iproc)
          iiAt=onWhichAtomPhi(ip%isorb+iorb)
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
    
      ! Transform to localization regions.
      do iorb=1,ip%norb_par(iproc)
          ilr=matmin%inWhichLocregExtracted(iorb)
          call vectorGlobalToLocal(ip%norbtotPad, matmin%mlr(ilr), coeffPad((iorb-1)*ip%norbtotPad+1), lcoeff(1,iorb))
      end do
    
      iterLoop: do it=1,lin%nItInguess
    
          if (iproc==0 .and. mod(it,1)==0) then
              write( *,'(1x,a,i0)') repeat('-',77 - int(log(real(it))/log(10.))) // ' iter=', it
          endif
    
    
          ! Othonormalize the coefficients.
          call orthonormalizeVectors(iproc, ip%nproc, lin%orbs, onWhichAtom, ip%onWhichMPI, ip%isorb_par, &
               matmin%norbmax, ip%norb_par(iproc), ip%isorb_par(iproc), lin%lzd%nlr, newComm, &
               matmin%mlr, lcoeff, comom)
    
          ! Calculate the gradient grad.
          do iorb=1,ip%norb_par(iproc)
              ilr=onWhichAtom(ip%isorb+iorb)
              iilr=matmin%inWhichLocregOnMPI(iorb)
              call dgemv('n', matmin%mlr(ilr)%norbinlr, matmin%mlr(ilr)%norbinlr, 1.d0, hamextract(1,1,iilr), matmin%norbmax, &
                   lcoeff(1,iorb), 1, 0.d0, lgrad(1,iorb), 1)
          end do
      
          !!do jorb=1,matmin%norbmax
          !!    write(650+iproc,'(100f15.5)') (lgrad(jorb,iorb), iorb=1,ip%norb_par(iproc))
          !!end do
          if(it>1) then
              traceOld=trace
          else
              traceOld=1.d10
          end if
          ! Apply the orthoconstraint to the gradient. To do so first calculate the Lagrange
          ! multiplier matrix.
          call orthoconstraintVectors(iproc, ip%nproc, lin%orbs, onWhichAtom, ip%onWhichMPI, ip%isorb_par, &
               matmin%norbmax, ip%norb_par(iproc), ip%isorb_par(iproc), lin%lzd%nlr, newComm, &
               matmin%mlr, lcoeff, lgrad, comom, trace)
    
          ! Calculate the gradient norm.
          fnrm=0.d0
          do iorb=1,ip%norb_par(iproc)
              ilr=onWhichAtom(ip%isorb+iorb)
              iilr=matmin%inWhichLocregOnMPI(iorb)
              fnrmArr(iorb)=ddot(matmin%mlr(ilr), lgrad(1,iorb), 1, lgrad(1,iorb), 1)
              if(it>1) fnrmOvrlpArr(iorb)=ddot(matmin%mlr(ilr), lgrad(1,iorb), 1, lgradold(1,iorb), 1)
          end do
          call mpi_barrier(newComm, ierr)
          call dcopy(ip%norb_par(iproc)*matmin%norbmax, lgrad(1,1), 1, lgradold(1,1), 1)
    
          ! Keep the gradient for the next iteration.
          if(it>1) then
              call dcopy(ip%norb_par(iproc), fnrmArr(1), 1, fnrmOldArr(1), 1)
          end if
    
          fnrmMax=0.d0
          meanAlpha=0.d0
          do iorb=1,ip%norb_par(iproc)
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
              meanAlpha=meanAlpha+alpha(iorb)
          end do
         call mpiallred(fnrm, 1, mpi_sum, newComm, ierr)
         call mpiallred(fnrmMax, 1, mpi_max, newComm, ierr)
         fnrm=sqrt(fnrm)
         fnrmMax=sqrt(fnrmMax)
    
         ! Determine the mean step size for steepest descent iterations.
         call mpiallred(meanAlpha, 1, mpi_sum, newComm, ierr)
         meanAlpha=meanAlpha/dble(ip%norb)

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
              ! Transform back to global ragion.
              do iorb=1,ip%norb_par(iproc)
                  ilr=matmin%inWhichLocregExtracted(iorb)
                  call vectorLocalToGlobal(ip%norbtotPad, matmin%mlr(ilr), lcoeff(1,iorb), coeffPad((iorb-1)*ip%norbtotPad+1))
              end do
              exit
          end if
      
          ! Quit if the maximal number of iterations is reached.
          if(it==lin%nItInguess) then
              if(iproc==0) write(*,'(x,a,i0,a)') 'WARNING: not converged within ', it, &
                  ' iterations! Exiting loop due to limitations of iterations.'
              if(iproc==0) write(*,'(x,a,2es15.7,f12.7)') 'Final values for fnrm, trace: ', fnrm, trace
              infoCoeff=-1
              ! Transform back to global ragion.
              do iorb=1,ip%norb_par(iproc)
                  ilr=matmin%inWhichLocregExtracted(iorb)
                  call vectorLocalToGlobal(ip%norbtotPad, matmin%mlr(ilr), lcoeff(1,iorb), coeffPad((iorb-1)*ip%norbtotPad+1))
              end do
              exit
          end if
    
          ! Improve the coefficients (by steepet descent).
          do iorb=1,ip%norb_par(iproc)
              ilr=onWhichAtom(ip%isorb+iorb)
              call daxpy(matmin%mlr(ilr)%norbinlr, -alpha(iorb), lgrad(1,iorb), 1, lcoeff(1,iorb), 1)
          end do

    
    
      end do iterLoop
    
    
      if(iproc==0) write(*,'(x,a)') '===================================================================================='
    
    
      ! Cut out the zeros
      allocate(coeff2(ip%norbtot*ip%norb_par(iproc)), stat=istat)
      call memocc(istat, coeff2, 'coeff2', subname)
      do iorb=1,ip%norb_par(iproc)
          call dcopy(ip%norbtot, coeffPad((iorb-1)*ip%norbtotPad+1), 1, coeff2((iorb-1)*ip%norbtot+1), 1)
      end do

      call deallocateArrays()

      iall=-product(shape(hamextract))*kind(hamextract)
      deallocate(hamextract, stat=istat)
      call memocc(istat, iall, 'hamextract', subname)


  end if processIf
  
  call mpi_barrier(mpi_comm_world, ierr)
  
  !! Allocate coeff2 for those processes which did not allocate it
  !! during the previous if statement.
  if(iproc>=ip%nproc) then
      allocate(coeff2(1), stat=istat)
      call memocc(istat, coeff2, 'coeff2', subname)
  end if
  
  
  ! Now collect all coefficients on all processes.
  allocate(recvcounts(0:nproc-1), stat=istat)
  call memocc(istat, recvcounts, 'recvcounts', subname)
  allocate(displs(0:nproc-1), stat=istat)
  call memocc(istat, displs, 'displs', subname)
  
  !!! Send ip%norb_par and ip%norbtot to all processes.
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
  ! Do this in a localized way -- TEST
  call buildLinearCombinations(iproc, nproc, lzdig, lin%lzd, input, coeff, lchi, lphi)

  if(iproc<ip%nproc) then
      call deallocate_inguessParameters(ip, subname)
      call deallocate_p2pCommsOrthonormalityMatrix(comom, subname)
      call deallocate_matrixMinimization(matmin,subname)

      iall=-product(shape(lcoeff))*kind(lcoeff)
      deallocate(lcoeff, stat=istat)
      call memocc(istat, iall, 'lcoeff', subname)

      iall=-product(shape(lgrad))*kind(lgrad)
      deallocate(lgrad, stat=istat)
      call memocc(istat, iall, 'lgrad', subname)

      iall=-product(shape(lgradold))*kind(lgradold)
      deallocate(lgradold, stat=istat)
      call memocc(istat, iall, 'lgradold', subname)
  end if

  ! Deallocate the local arrays.
  iall=-product(shape(newID))*kind(newID)
  deallocate(newID, stat=istat)
  call memocc(istat, iall, 'newID', subname)

  iall=-product(shape(coeff))*kind(coeff)
  deallocate(coeff, stat=istat)
  call memocc(istat, iall, 'coeff', subname)

  iall=-product(shape(coeff2))*kind(coeff2)
  deallocate(coeff2, stat=istat)
  call memocc(istat, iall, 'coeff2', subname)

  iall=-product(shape(recvcounts))*kind(recvcounts)
  deallocate(recvcounts, stat=istat)
  call memocc(istat, iall, 'recvcounts', subname)

  iall=-product(shape(displs))*kind(displs)
  deallocate(displs, stat=istat)
  call memocc(istat, iall, 'displs', subname)

  iall=-product(shape(norb_par))*kind(norb_par)
  deallocate(norb_par, stat=istat)
  call memocc(istat, iall, 'norb_par', subname)
  
  
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


end subroutine buildLinearCombinationsLocalized





subroutine determineLocalizationRegions(iproc, nproc, nlr, norb, at, onWhichAtomAll, locrad, rxyz, lzd, mlr)
use module_base
use module_types
use module_interfaces, exceptThisOne => determineLocalizationRegions
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, nlr, norb
type(atoms_data),intent(in):: at
integer,dimension(norb),intent(in):: onWhichAtomAll
real(8),dimension(at%nat),intent(in):: locrad
real(8),dimension(3,at%nat),intent(in):: rxyz
type(linear_zone_descriptors),intent(in):: lzd
type(matrixLocalizationRegion),dimension(:),pointer,intent(out):: mlr

! Local variables
integer:: ilr, jlr, jorb, ii, istat, is1, ie1, is2, ie2, is3, ie3, js1, je1, js2, je2, js3, je3
real(8):: cut, tt
logical:: ovrlpx, ovrlpy, ovrlpz
character(len=*),parameter:: subname='determineLocalizationRegions'


allocate(mlr(nlr), stat=istat)

! Count for each localization region the number of matrix elements within the cutoff.
do ilr=1,nlr
    mlr(ilr)%norbinlr=0
    call getIndices(lzd%llr(ilr), is1, ie1, is2, ie2, is3, ie3)
    do jorb=1,norb
        jlr=onWhichAtomAll(jorb)
        call getIndices(lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
        ovrlpx = ( is1<=je1 .and. ie1>=js1 )
        ovrlpy = ( is2<=je2 .and. ie2>=js2 )
        ovrlpz = ( is3<=je3 .and. ie3>=js3 )
        if(ovrlpx .and. ovrlpy .and. ovrlpz) then
            mlr(ilr)%norbinlr=mlr(ilr)%norbinlr+1
        end if
    end do
    allocate(mlr(ilr)%indexInGlobal(mlr(ilr)%norbinlr), stat=istat)
    call memocc(istat, mlr(ilr)%indexInGlobal, 'mlr(ilr)%indexInGlobal', subname)
    !if(iproc==0) write(*,'(a,i4,i7)') 'ilr, mlr(ilr)%norbinlr', ilr, mlr(ilr)%norbinlr
end do


! Now determine the indices of the elements with an overlap.
do ilr=1,nlr
    ii=0
    call getIndices(lzd%llr(ilr), is1, ie1, is2, ie2, is3, ie3)
    do jorb=1,norb
        jlr=onWhichAtomAll(jorb)
        call getIndices(lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
        ovrlpx = ( is1<=je1 .and. ie1>=js1 )
        ovrlpy = ( is2<=je2 .and. ie2>=js2 )
        ovrlpz = ( is3<=je3 .and. ie3>=js3 )
        if(ovrlpx .and. ovrlpy .and. ovrlpz) then
            ii=ii+1
            mlr(ilr)%indexInGlobal(ii)=jorb
        end if
    end do
    if(ii/=mlr(ilr)%norbinlr) then
        write(*,'(a,i0,a,2(2x,i0))') 'ERROR on process ', iproc, ': ii/=mlr(ilr)%norbinlr', ii, mlr(ilr)%norbinlr
    end if
    !if(iproc==0) write(*,'(a,i6,200i5)') 'ilr, mlr(ilr)%indexInGlobal(ii)', ilr, mlr(ilr)%indexInGlobal(:)
end do


end subroutine determineLocalizationRegions




subroutine extractMatrix(iproc, nproc, norb, norbp, orbstot, onWhichAtomPhi, onWhichMPI, nmat, ham, matmin, hamextract)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, nmat, norb, norbp
type(orbitals_data),intent(in):: orbstot
integer,dimension(norb),intent(in):: onWhichAtomPhi, onWhichMPI
real(8),dimension(orbstot%norb,orbstot%norb,nmat),intent(in):: ham
type(matrixMinimization),intent(inout):: matmin
real(8),dimension(:,:,:),pointer,intent(out):: hamextract

! Local variables
integer:: jorb, jlr, jproc, jlrold, jjorb, ind, indlarge, jnd, jndlarge, ii, istat, jjlr
character(len=*),parameter:: subname='extractMatrix'

allocate(matmin%inWhichLocregExtracted(norbp), stat=istat)
call memocc(istat, matmin%inWhichLocregExtracted, 'matmin%inWhichLocregExtracted', subname)
allocate(matmin%inWhichLocregOnMPI(norbp), stat=istat)
call memocc(istat, matmin%inWhichLocregOnMPI, 'matmin%inWhichLocregOnMPI', subname)

! Allocate the matrices holding the extracted quantities. In principle this matrix
! has a difefrent size for each localization region. To simplify the program, allocate them
! with the same size for all localization regions on a given MPI process.
matmin%norbmax=0
jlrold=-1
matmin%nlrp=0
jjorb=0
do jorb=1,norb
    jlr=onWhichAtomPhi(jorb)
    jproc=onWhichMPI(jorb)
    if(iproc==jproc) then
        jjorb=jjorb+1
        if(jlr>jlrold) then
           matmin%nlrp=matmin%nlrp+1
        end if
        if(matmin%mlr(jlr)%norbinlr>matmin%norbmax) then
            matmin%norbmax=matmin%mlr(jlr)%norbinlr
        end if
        matmin%inWhichLocregExtracted(jjorb)=jlr
        matmin%inWhichLocregOnMPI(jjorb)=matmin%nlrp
        jlrold=jlr
    end if
end do

allocate(matmin%indexInLocreg(matmin%nlrp), stat=istat)
call memocc(istat, matmin%indexInLocreg, 'matmin%indexInLocreg', subname)

! Allocate the matrix
allocate(hamextract(matmin%norbmax,matmin%norbmax,matmin%nlrp), stat=istat)
call memocc(istat, hamextract, 'hamextract', subname)
hamextract=0.d0

! Exctract the data from the large Hamiltonian.
jlrold=-1
jjlr=0
do jorb=1,norb
    jlr=onWhichAtomPhi(jorb)
    jproc=onWhichMPI(jorb)
    if(iproc==jproc) then
        if(jlr>jlrold) then
            jjlr=jjlr+1
            matmin%indexInLocreg(jjlr)=jlr
            ! To make it work for both input guess (where we have nmat>1 different matrices) and
            ! for the iterative diagonalization (where we have only nmat=1 matrix).
            ii=min(jlr,nmat)
            do ind=1,matmin%mlr(jlr)%norbinlr
                indlarge=matmin%mlr(jlr)%indexInGlobal(ind)
                do jnd=1,matmin%mlr(jlr)%norbinlr
                    jndlarge=matmin%mlr(jlr)%indexInGlobal(jnd)
                    hamextract(jnd,ind,jjlr)=ham(jndlarge,indlarge,ii)
                end do
            end do
            jlrold=jlr
        end if
    end if
end do

!!do jlr=1,nmat
!!    do ind=1,orbstot%norb
!!        write(200+10*iproc+jlr,'(100es9.1)') (ham(ind,jnd,jlr), jnd=1,orbstot%norb)
!!    end do
!!end do
!!
!!do jlr=1,matmin%nlrp
!!    jjlr=matmin%indexInLocreg(jlr)
!!    do ind=1,matmin%mlr(jjlr)%norbinlr
!!        write(100+10*iproc+jlr,'(100es9.1)') (hamextract(ind,jnd,jlr), jnd=1,matmin%mlr(jjlr)%norbinlr)
!!    end do
!!end do




end subroutine extractMatrix




subroutine vectorGlobalToLocal(norbtot, mlr, vglobal, vlocal)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: norbtot
type(matrixLocalizationRegion),intent(in):: mlr
real(8),dimension(norbtot),intent(in):: vglobal
real(8),dimension(mlr%norbinlr),intent(out):: vlocal

! Local variables
integer:: ilocal, iglobal

do ilocal=1,mlr%norbinlr
    iglobal=mlr%indexInGlobal(ilocal)
    vlocal(ilocal)=vglobal(iglobal)
end do


end subroutine vectorGlobalToLocal




subroutine vectorLocalToGlobal(norbtot, mlr, vlocal, vglobal)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: norbtot
type(matrixLocalizationRegion),intent(in):: mlr
real(8),dimension(mlr%norbinlr),intent(in):: vlocal
real(8),dimension(norbtot),intent(out):: vglobal

! Local variables
integer:: ilocal, iglobal

vglobal=0.d0
do ilocal=1,mlr%norbinlr
    iglobal=mlr%indexInGlobal(ilocal)
    vglobal(iglobal)=vlocal(ilocal)
end do


end subroutine vectorLocalToGlobal





subroutine determineOverlapRegionMatrix(iproc, nproc, lzd, mlr, orbs, orbstot, onWhichAtom, onWhichAtomPhi, comom)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(linear_zone_descriptors),intent(in):: lzd
type(orbitals_data),intent(in):: orbs, orbstot
integer,dimension(orbstot%norb),intent(in):: onWhichAtom
integer,dimension(orbs%norb),intent(in):: onWhichAtomPhi
type(matrixLocalizationRegion),dimension(lzd%nlr),intent(in):: mlr
type(p2pCommsOrthonormalityMatrix),intent(out):: comom

! Local variables
integer:: ilr, jlr, klr, novrlp, korb, istat, jlrold, jjlr, jjorb, jorb, kkorb, lorb, iorb
integer:: is1, ie1, is2, ie2, is3, ie3, js1, je1, js2, je2, js3, je3, ks1, ke1, ks2, ke2, ks3, ke3
logical:: ovrlpx_ki, ovrlpy_ki, ovrlpz_ki, ovrlpx_kj, ovrlpy_kj, ovrlpz_kj, ovrlpx, ovrlpy, ovrlpz
logical:: overlapFound
character(len=*),parameter:: subname='determineOverlapRegionMatrix'


allocate(comom%noverlap(lzd%nlr), stat=istat)
call memocc(istat, comom%noverlap, 'comom%noverlap', subname)

! First count the number of overlapping localization regions for each localization region.
do ilr=1,lzd%nlr
    call getIndices(lzd%llr(ilr), is1, ie1, is2, ie2, is3, ie3)
    novrlp=0
    do jorb=1,orbs%norb
        jlr=onWhichAtomPhi(jorb)
        call getIndices(lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
        ovrlpx = ( is1<=je1 .and. ie1>=js1 )
        ovrlpy = ( is2<=je2 .and. ie2>=js2 )
        ovrlpz = ( is3<=je3 .and. ie3>=js3 )
        if(ovrlpx .and. ovrlpy .and. ovrlpz) then
            novrlp=novrlp+1
        end if
    end do
    comom%noverlap(ilr)=novrlp
    !!if(iproc==0) write(*,*) 'ilr, comom%noverlap(ilr)', ilr, comom%noverlap(ilr) 
end do


allocate(comom%overlaps(maxval(comom%noverlap(:)),lzd%nlr), stat=istat)
call memocc(istat, comom%overlaps, 'comom%overlaps', subname)
do ilr=1,lzd%nlr
    comom%overlaps(:,ilr)=0
    call getIndices(lzd%llr(ilr), is1, ie1, is2, ie2, is3, ie3)
    novrlp=0
    do jorb=1,orbs%norb
        jlr=onWhichAtomPhi(jorb)
        call getIndices(lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
        ovrlpx = ( is1<=je1 .and. ie1>=js1 )
        ovrlpy = ( is2<=je2 .and. ie2>=js2 )
        ovrlpz = ( is3<=je3 .and. ie3>=js3 )
        if(ovrlpx .and. ovrlpy .and. ovrlpz) then
            novrlp=novrlp+1
            comom%overlaps(novrlp,ilr)=jorb
        end if
    end do
    !!if(iproc==0) write(*,'(a,i4,3x,100i5)') 'ilr, comom%overlaps(,ilr)', ilr, comom%overlaps(:,ilr) 
end do

allocate(comom%olr(maxval(comom%noverlap(:)),lzd%nlr), stat=istat)



! Now determine which orbitals (corresponding to basis functions) will be in the overlap localization region.
do ilr=1,lzd%nlr
    call getIndices(lzd%llr(ilr), is1, ie1, is2, ie2, is3, ie3)
    comom%olr(:,ilr)%norbinlr=0
    do jorb=1,comom%noverlap(ilr)
        jjorb=comom%overlaps(jorb,ilr)
        jlr=onWhichAtomPhi(jjorb)
        call getIndices(lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
        do korb=1,mlr(jlr)%norbinlr
            lorb=mlr(jlr)%indexInGlobal(korb)
            klr=onWhichAtom(lorb)
            call getIndices(lzd%llr(klr), ks1, ke1, ks2, ke2, ks3, ke3)
            ovrlpx_ki = ( ks1<=ie1 .and. ke1>=is1 )
            ovrlpy_ki = ( ks2<=ie2 .and. ke2>=is2 )
            ovrlpz_ki = ( ks3<=ie3 .and. ke3>=is3 )
            ovrlpx_kj = ( ks1<=je1 .and. ke1>=js1 )
            ovrlpy_kj = ( ks2<=je2 .and. ke2>=js2 )
            ovrlpz_kj = ( ks3<=je3 .and. ke3>=js3 )
            ovrlpx = ( ovrlpx_ki .and. ovrlpx_kj )
            ovrlpy = ( ovrlpy_ki .and. ovrlpy_kj )
            ovrlpz = ( ovrlpz_ki .and. ovrlpz_kj )
            if(ovrlpx .and. ovrlpy .and. ovrlpz) then
                comom%olr(jorb,ilr)%norbinlr=comom%olr(jorb,ilr)%norbinlr+1
            end if
        end do
        allocate(comom%olr(jorb,ilr)%indexInGlobal(comom%olr(jorb,ilr)%norbinlr), stat=istat)
        call memocc(istat, comom%olr(jorb,ilr)%indexInGlobal, 'comom%olr(jorb,ilr)%indexInGlobal', subname)
    end do
end do

!!do ilr=1,lzd%nlr
!!    do jorb=1,comom%noverlap(ilr)
!!        if(iproc==0) write(*,'(a,2i5,2i8)') 'ilr, jjorb, comom%overlaps(jorb,ilr), comom%olr(jorb,ilr)%norbinlr', ilr, jorb, comom%overlaps(jorb,ilr), comom%olr(jorb,ilr)%norbinlr
!!    end do
!!end do


! Determine the indices to switch from global region to localization region.
do ilr=1,lzd%nlr
    call getIndices(lzd%llr(ilr), is1, ie1, is2, ie2, is3, ie3)
    do jorb=1,comom%noverlap(ilr)
        jjorb=comom%overlaps(jorb,ilr)
        jlr=onWhichAtomPhi(jjorb)
        call getIndices(lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
        kkorb=0
        comom%olr(jorb,ilr)%indexInGlobal(:)=0
        do korb=1,mlr(jlr)%norbinlr
            lorb=mlr(jlr)%indexInGlobal(korb)
            klr=onWhichAtom(lorb)
            call getIndices(lzd%llr(klr), ks1, ke1, ks2, ke2, ks3, ke3)
            ovrlpx_ki = ( ks1<=ie1 .and. ke1>=is1 )
            ovrlpy_ki = ( ks2<=ie2 .and. ke2>=is2 )
            ovrlpz_ki = ( ks3<=ie3 .and. ke3>=is3 )
            ovrlpx_kj = ( ks1<=je1 .and. ke1>=js1 )
            ovrlpy_kj = ( ks2<=je2 .and. ke2>=js2 )
            ovrlpz_kj = ( ks3<=je3 .and. ke3>=js3 )
            ovrlpx = ( ovrlpx_ki .and. ovrlpx_kj )
            ovrlpy = ( ovrlpy_ki .and. ovrlpy_kj )
            ovrlpz = ( ovrlpz_ki .and. ovrlpz_kj )
            if(ovrlpx .and. ovrlpy .and. ovrlpz) then
                kkorb=kkorb+1
                comom%olr(jorb,ilr)%indexInGlobal(kkorb)=korb
            end if
        end do
    end do
end do

! With these indices it is possible to extract data from the global region to the
! overlap region. For example: comom%olr(jorb,ilr) allows to extract data from orbital
! jorb (the jorb-th orbital overlapping with region ilr) to the overlap region. To expand
! this overlap region to the whole region ilr, we need comom%(iorb,jlr), where jlr is the 
! localization region of jorb and iorb the iorb-th overlapping orbital of region jlr.
! This information is stored in comom%olrForExpansion:
! comom%olrForExpansion(1,jorb,ilr)=jlr
! comom%olrForExpansion(2,jorb,ilr)=iorb
allocate(comom%olrForExpansion(2,maxval(comom%noverlap(:)),lzd%nlr), stat=istat)
call memocc(istat, comom%olrForExpansion, 'comom%olrForExpansion', subname)
do ilr=1,lzd%nlr
    do iorb=1,comom%noverlap(ilr)
        jorb=comom%overlaps(iorb,ilr)
        jlr=onWhichAtom(jorb)
        comom%olrForExpansion(1,iorb,ilr)=jlr
        do korb=1,comom%noverlap(jlr)
             kkorb=comom%overlaps(korb,jlr)
             klr=onWhichAtom(kkorb)
             if(klr==ilr) then
                 comom%olrForExpansion(2,iorb,ilr)=korb
             end if
        end do
    end do
end do




!!do ilr=1,lzd%nlr
!!    do jorb=1,comom%noverlap(ilr)
!!        if(iproc==0) write(*,'(a,2i5,100i5)') 'ilr, jorb, comom%olr(jorb,ilr)%indexInGlobal', ilr, jorb, comom%olr(jorb,ilr)%indexInGlobal
!!    end do
!!end do


end subroutine determineOverlapRegionMatrix





subroutine initCommsMatrixOrtho(iproc, nproc, norb, norb_par, isorb_par, onWhichAtomPhi, onWhichMPI, tag, comom)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, norb
integer,dimension(norb),intent(in):: onWhichAtomPhi, onWhichMPI
integer,dimension(0:nproc-1),intent(in):: norb_par, isorb_par
integer,intent(inout):: tag
type(p2pCommsOrthonormalityMatrix),intent(inout):: comom

! Local variables
integer:: jlrold, jproc, jj, jorb, jjorb, jlr, jjmax, istat, jkorb, mpisource, mpidest, istsource, istdest, ncount, korb, iall, kkorb
integer,dimension(:),allocatable:: istsourcearr, istdestarr
character(len=*),parameter:: subname='initCommsMatrixOrtho'


allocate(istsourcearr(0:nproc-1), stat=istat)
call memocc(istat, istsourcearr, 'istsourcearr', subname)
istsourcearr=1
allocate(istdestarr(0:nproc-1), stat=istat)
call memocc(istat, istdestarr, 'istdestarr', subname)
istdestarr=1

comom%nrecvbuf=0
allocate(comom%indexInRecvBuf(maxval(comom%noverlap(:)),0:nproc-1), stat=istat)
call memocc(istat, comom%indexInRecvBuf, 'comom%indexInRecvBuf', subname)

allocate(comom%noverlapProc(0:nproc-1), stat=istat)
call memocc(istat, comom%noverlapProc, 'comom%noverlapProc', subname)

! Count how many orbitals each process will receive
do jproc=0,nproc-1
    jlrold=0
    comom%noverlapProc(jproc)=0
    do jorb=1,norb_par(jproc)
        jjorb=isorb_par(jproc)+jorb
        jlr=onWhichAtomPhi(jjorb)
        if(jlr==jlrold) cycle
        do korb=1,comom%noverlap(jlr)
           comom%noverlapProc(jproc)=comom%noverlapProc(jproc)+1
        end do
        jlrold=jlr
    end do 
end do

allocate(comom%comarr(8,maxval(comom%noverlapProc(:)),0:nproc-1), stat=istat)
call memocc(istat, comom%comarr, 'comom%comarr', subname)
allocate(comom%overlapsProc(maxval(comom%noverlapProc(:)),0:nproc-1), stat=istat)
call memocc(istat, comom%overlapsProc, 'comom%overlapsProc', subname)
allocate(comom%communComplete(maxval(comom%noverlapProc(:)),0:nproc-1), stat=istat)
call memocc(istat, comom%communComplete, 'comom%communComplete', subname)

comom%nsendBuf=0
comom%nrecvBuf=0
do jproc=0,nproc-1
    jkorb=0
    jlrold=0
    do jorb=1,norb_par(jproc)
        jjorb=isorb_par(jproc)+jorb
        jlr=onWhichAtomPhi(jjorb)
        if(jlr==jlrold) cycle
        do korb=1,comom%noverlap(jlr)
            jkorb=jkorb+1
            kkorb=comom%overlaps(korb,jlr)
            mpidest=jproc
            mpisource=onWhichMPI(kkorb)
            istsource=istsourcearr(mpisource)
            ncount=comom%olr(korb,jlr)%norbinlr 
            istdest=istdestarr(mpidest)
            tag=tag+1
            comom%overlapsProc(jkorb,jproc)=kkorb
            call setCommsParameters(mpisource, mpidest, istsource, istdest, ncount, tag, comom%comarr(1,jkorb,jproc))
            if(iproc==mpisource) then
                comom%nsendBuf=comom%nsendBuf+ncount
                !write(*,'(3(a,i0),6x,a,4i8)') 'process ',iproc,' adds ',ncount,' elements to comom%nsendBuf. Value after that: ',comom%nsendBuf,'. jjorb, jlr, korb, kkorb', jjorb, jlr, korb, kkorb
            end if
            if(iproc==mpidest) then
                !write(*,'(2(a,i0),a)') 'process ',iproc,' adds ',ncount,' elements to comom%nrecvBuf.'
                comom%nrecvbuf=comom%nrecvbuf+ncount
                comom%indexInRecvBuf(korb,jproc)=istdest
            end if
            istdestarr(mpidest)=istdestarr(mpidest)+ncount
            istsourcearr(mpisource)=istsourcearr(mpisource)+ncount
        end do
        jlrold=jlr
    end do 
end do

allocate(comom%recvBuf(comom%nrecvbuf), stat=istat)
call memocc(istat, comom%recvBuf, 'comom%recvBuf', subname)
allocate(comom%sendBuf(comom%nsendbuf), stat=istat)
call memocc(istat, comom%sendBuf, 'comom%sendBuf', subname)

iall=-product(shape(istsourcearr))*kind(istsourcearr)
deallocate(istsourcearr, stat=istat)
call memocc(istat, iall, 'istsourcearr', subname)
iall=-product(shape(istdestarr))*kind(istdestarr)
deallocate(istdestarr, stat=istat)
call memocc(istat, iall, 'istdestarr', subname)

end subroutine initCommsMatrixOrtho




subroutine orthonormalizeVectors(iproc, nproc, orbs, onWhichAtom, onWhichMPI, isorb_par, norbmax, norbp, isorb, nlr, newComm, mlr, vec, comom)
use module_base
use module_types
use module_interfaces, exceptThisOne => orthonormalizeVectors
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, norbmax, norbp, isorb, nlr, newComm
type(orbitals_data),intent(in):: orbs
integer,dimension(orbs%norb),intent(in):: onWhichAtom, onWhichMPI
integer,dimension(0:nproc-1),intent(in):: isorb_par
type(matrixLocalizationRegion),dimension(nlr),intent(in):: mlr
real(8),dimension(norbmax,norbp),intent(inout):: vec
type(p2pCommsOrthonormalityMatrix),intent(inout):: comom

! Local variables
integer:: noverlaps, iorb, iiorb, ilr, istat, ilrold, jorb, iall
real(8),dimension(:,:),allocatable:: vecOvrlp, ovrlp
character(len=*),parameter:: subname='orthonormalizeVectors'

noverlaps=0
ilrold=0
do iorb=1,norbp
    iiorb=isorb+iorb
    ilr=onWhichAtom(iiorb)
    if(ilr>ilrold) then
        noverlaps=noverlaps+comom%noverlap(ilr)
    end if
    ilrold=ilr
end do
allocate(vecOvrlp(norbmax,noverlaps), stat=istat)
call memocc(istat, vecOvrlp, 'vecOvrlp', subname)
allocate(ovrlp(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, ovrlp, 'ovrlp', subname)

call extractToOverlapregion(iproc, nproc, orbs%norb, onWhichAtom, onWhichMPI, isorb_par, norbmax, norbp, vec, comom)
call postCommsVectorOrthonormalization(iproc, nproc, newComm, comom)
call gatherVectors(iproc, nproc, newComm, comom)

call expandFromOverlapregion(iproc, nproc, isorb, norbp, orbs, onWhichAtom, comom, norbmax, noverlaps, vecOvrlp)

! Calculate the overlap matrix.
call calculateOverlap(iproc, nproc, nlr, norbmax, norbp, noverlaps, isorb, orbs%norb, comom, mlr, onWhichAtom, vec, vecOvrlp, newComm, ovrlp)
call transformOverlapMatrix(iproc, nproc, orbs%norb, ovrlp)
call orthonormalLinearCombinations(iproc, nproc, nlr, norbmax, norbp, noverlaps, isorb, orbs%norb, comom, mlr, onWhichAtom, vecOvrlp, ovrlp, vec)

iall=-product(shape(ovrlp))*kind(ovrlp)
deallocate(ovrlp, stat=istat)
call memocc(istat, iall, 'ovrlp', subname)
iall=-product(shape(vecOvrlp))*kind(vecOvrlp)
deallocate(vecOvrlp, stat=istat)
call memocc(istat, iall, 'vecOvrlp', subname)

end subroutine orthonormalizeVectors




subroutine orthoconstraintVectors(iproc, nproc, orbs, onWhichAtom, onWhichMPI, isorb_par, norbmax, norbp, isorb, nlr, newComm, mlr, vec, grad, comom, trace)
use module_base
use module_types
use module_interfaces, exceptThisOne => orthoconstraintVectors
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, norbmax, norbp, isorb, nlr, newComm
type(orbitals_data),intent(in):: orbs
integer,dimension(orbs%norb),intent(in):: onWhichAtom, onWhichMPI
integer,dimension(0:nproc-1),intent(in):: isorb_par
type(matrixLocalizationRegion),dimension(nlr),intent(in):: mlr
real(8),dimension(norbmax,norbp),intent(inout):: vec, grad
type(p2pCommsOrthonormalityMatrix),intent(inout):: comom
real(8),intent(out):: trace

! Local variables
integer:: noverlaps, iorb, iiorb, ilr, istat, ilrold, jorb, iall
real(8),dimension(:,:),allocatable:: gradOvrlp, vecOvrlp, lagmat, ovrlp
character(len=*),parameter:: subname='orthonormalizeVectors'


noverlaps=0
ilrold=0
do iorb=1,norbp
    iiorb=isorb+iorb
    ilr=onWhichAtom(iiorb)
    if(ilr>ilrold) then
        noverlaps=noverlaps+comom%noverlap(ilr)
    end if
    ilrold=ilr
end do
allocate(gradOvrlp(norbmax,noverlaps), stat=istat)
call memocc(istat, gradOvrlp, 'gradOvrlp', subname)
allocate(vecOvrlp(norbmax,noverlaps), stat=istat)
call memocc(istat, vecOvrlp, 'vecOvrlp', subname)
allocate(lagmat(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, lagmat, 'lagmat', subname)
allocate(ovrlp(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, ovrlp, 'ovrlp', subname)

call extractToOverlapregion(iproc, nproc, orbs%norb, onWhichAtom, onWhichMPI, isorb_par, norbmax, norbp, grad, comom)
call postCommsVectorOrthonormalization(iproc, nproc, newComm, comom)
call gatherVectors(iproc, nproc, newComm, comom)

call expandFromOverlapregion(iproc, nproc, isorb, norbp, orbs, onWhichAtom, comom, norbmax, noverlaps, gradOvrlp)

! Calculate the Lagrange multiplier matrix <vec|grad>.
call calculateOverlap(iproc, nproc, nlr, norbmax, norbp, noverlaps, isorb, orbs%norb, comom, mlr, onWhichAtom, vec, gradOvrlp, newComm, lagmat)
!!do iorb=1,orbs%norb
!!    do jorb=1,orbs%norb
!!        write(500+iproc,*) iorb, jorb, lagmat(iorb,jorb)
!!    end do
!!end do
trace=0.d0
do iorb=1,orbs%norb
    trace=trace+lagmat(iorb,iorb)
end do

! Now we also have to calculate the overlap matrix.
call extractToOverlapregion(iproc, nproc, orbs%norb, onWhichAtom, onWhichMPI, isorb_par, norbmax, norbp, vec, comom)
call postCommsVectorOrthonormalization(iproc, nproc, newComm, comom)
call gatherVectors(iproc, nproc, newComm, comom)
call expandFromOverlapregion(iproc, nproc, isorb, norbp, orbs, onWhichAtom, comom, norbmax, noverlaps, vecOvrlp)
call calculateOverlap(iproc, nproc, nlr, norbmax, norbp, noverlaps, isorb, orbs%norb, comom, mlr, onWhichAtom, vec, vecOvrlp, newComm, ovrlp)

! Now apply the orthoconstraint.
call applyOrthoconstraintVectors(iproc, nproc, orbs%norb, norbmax, norbp, isorb, nlr, noverlaps, onWhichAtom, vecOvrlp, ovrlp, lagmat, comom, mlr, grad)

!call transformOverlapMatrix(iproc, nproc, orbs%norb, ovrlp)
!call orthonormalLinearCombinations(iproc, nproc, nlr, norbmax, norbp, noverlaps, isorb, orbs%norb, comom, mlr, onWhichAtom, vecOvrlp, ovrlp, vec)

iall=-product(shape(lagmat))*kind(lagmat)
deallocate(lagmat, stat=istat)
call memocc(istat, iall, 'lagmat', subname)

iall=-product(shape(gradOvrlp))*kind(gradOvrlp)
deallocate(gradOvrlp, stat=istat)
call memocc(istat, iall, 'gradOvrlp', subname)

iall=-product(shape(ovrlp))*kind(ovrlp)
deallocate(ovrlp, stat=istat)
call memocc(istat, iall, 'ovrlp', subname)

iall=-product(shape(vecOvrlp))*kind(vecOvrlp)
deallocate(vecOvrlp, stat=istat)
call memocc(istat, iall, 'vecOvrlp', subname)

end subroutine orthoconstraintVectors





subroutine postCommsVectorOrthonormalization(iproc, nproc, newComm, comom)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, newComm
type(p2pCommsOrthonormalityMatrix),intent(inout):: comom

! Local variables
integer:: nsends, nreceives, jproc, iorb, mpisource, istsource, ncount, mpidest, istdest, tag, ierr

nsends=0
nreceives=0
comom%communComplete=.false.
do jproc=0,nproc-1
    do iorb=1,comom%noverlapProc(jproc)
        mpisource=comom%comarr(1,iorb,jproc)
        istsource=comom%comarr(2,iorb,jproc)
        ncount=comom%comarr(3,iorb,jproc)
        mpidest=comom%comarr(4,iorb,jproc)
        istdest=comom%comarr(5,iorb,jproc)
        tag=comom%comarr(6,iorb,jproc)
        if(mpisource/=mpidest) then
            ! The orbitals are on different processes, so we need a point to point communication.
            if(iproc==mpisource) then
                !write(*,'(6(a,i0))') 'process ', mpisource, ' sends ', ncount, ' elements from position ', istsource, ' to position ', istdest, ' on process ', mpidest, ', tag=',tag
                call mpi_isend(comom%sendBuf(istsource), ncount, mpi_double_precision, mpidest, tag, newComm, comom%comarr(7,iorb,jproc), ierr)
                !call mpi_isend(sendBuf(istsource), ncount, mpi_double_precision, mpidest, tag, mpi_comm_world, lin%comsr%comarr(8,iorb,jproc), ierr)
                comom%comarr(8,iorb,jproc)=mpi_request_null !is this correct?
                nsends=nsends+1
            else if(iproc==mpidest) then
                !write(*,'(6(a,i0))') 'process ', mpidest, ' receives ', ncount, ' elements at position ', istdest, ' from position ', istsource, ' on process ', mpisource, ', tag=',tag
                call mpi_irecv(comom%recvBuf(istdest), ncount, mpi_double_precision, mpisource, tag, newComm, comom%comarr(8,iorb,jproc), ierr)
                comom%comarr(7,iorb,jproc)=mpi_request_null !is this correct?
                nreceives=nreceives+1
            else
                comom%comarr(7,iorb,jproc)=mpi_request_null
                comom%comarr(8,iorb,jproc)=mpi_request_null
            end if
        else
            ! The orbitals are on the same process, so simply copy them.
            if(iproc==mpisource) then
                call dcopy(ncount, comom%sendBuf(istsource), 1, comom%recvBuf(istdest), 1)
                !write(*,'(6(a,i0))') 'process ', iproc, ' copies ', ncount, ' elements from position ', istsource, ' to position ', istdest, ' on process ', iproc, ', tag=',tag
                comom%comarr(7,iorb,jproc)=mpi_request_null
                comom%comarr(8,iorb,jproc)=mpi_request_null
                nsends=nsends+1
                nreceives=nreceives+1
                comom%communComplete(iorb,iproc)=.true.
            else
                comom%comarr(7,iorb,jproc)=mpi_request_null
                comom%comarr(8,iorb,jproc)=mpi_request_null
            end if

        end if
    end do
end do



end subroutine postCommsVectorOrthonormalization




subroutine gatherVectors(iproc, nproc, newComm, comom)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, newComm
type(p2pCommsOrthonormalityMatrix),intent(inout):: comom

! Local variables
integer:: jorb, mpisource, mpidest, nfast, nslow, nsameproc, ierr, jproc
integer,dimension(mpi_status_size):: stat
logical:: sendComplete, receiveComplete

! Check whether the communications have completed.
nfast=0
nsameproc=0
testLoop: do
    do jproc=0,nproc-1
        do jorb=1,comom%noverlapProc(jproc)
            if(comom%communComplete(jorb,jproc)) cycle
            call mpi_test(comom%comarr(7,jorb,jproc), sendComplete, stat, ierr)
            call mpi_test(comom%comarr(8,jorb,jproc), receiveComplete, stat, ierr)
            if(sendComplete .and. receiveComplete) comom%communComplete(jorb,jproc)=.true.
            if(comom%communComplete(jorb,jproc)) then
                !write(*,'(2(a,i0))') 'fast communication; process ', iproc, ' has received orbital ', jorb
                mpisource=comom%comarr(1,jorb,jproc)
                mpidest=comom%comarr(4,jorb,jproc)
                if(mpisource/=mpidest) then
                    nfast=nfast+1
                else
                    nsameproc=nsameproc+1
                end if
            end if
        end do
    end do
    ! If we made it until here, either all all the communication is
    ! complete or we better wait for each single orbital.
    exit testLoop
end do testLoop


! Wait for the communications that have not completed yet
nslow=0
do jproc=0,nproc-1
    do jorb=1,comom%noverlapProc(jproc)
        if(comom%communComplete(jorb,jproc)) cycle
        !write(*,'(2(a,i0))') 'process ', iproc, ' is waiting for orbital ', korb
        nslow=nslow+1
        call mpi_wait(comom%comarr(7,jorb,jproc), stat, ierr)   !COMMENTED BY PB
        call mpi_wait(comom%comarr(8,jorb,jproc), stat, ierr)   !COMMENTED BY PB
        comom%communComplete(jorb,jproc)=.true.
    end do
end do

!call mpiallred(nreceives, 1, mpi_sum, mpi_comm_world, ierr)
call mpiallred(nfast, 1, mpi_sum, newComm, ierr)
call mpiallred(nslow, 1, mpi_sum, newComm, ierr)
call mpiallred(nsameproc, 1, mpi_sum, newComm, ierr)
!if(iproc==0) write(*,'(x,2(a,i0),a)') 'statistics: - ', nfast+nslow, ' point to point communications, of which ', &
!                       nfast, ' could be overlapped with computation.'
!if(iproc==0) write(*,'(x,a,i0,a)') '            - ', nsameproc, ' copies on the same processor.'


end subroutine gatherVectors




subroutine extractToOverlapregion(iproc, nproc, norb, onWhichAtom, onWhichMPI, isorb_par, norbmax, norbp, vec, comom)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, norbmax, norb, norbp
integer,dimension(norb),intent(in):: onWhichAtom, onWhichMPI
integer,dimension(0:nproc-1),intent(in):: isorb_par
real(8),dimension(norbmax,norbp),intent(in):: vec
type(p2pCommsOrthonormalityMatrix),intent(inout):: comom

! Local variables
integer:: ilrold, iiprocold, ist, ilr, iiproc, jjorb, jorb, iorb, jjproc, korb, ind, i, jjlr

!write(*,*) 'iproc, norbmax, norbp', iproc, norbmax, norbp

ilrold=-1
iiprocold=-1
ist=0
do iorb=1,norb
    ilr=onWhichAtom(iorb)
    iiproc=onWhichMPI(iorb)
    if(ilr==ilrold .and.  iiproc==iiprocold) cycle ! otherwise we would extract the same again
    do jorb=1,comom%noverlap(ilr)
        jjorb=comom%overlaps(jorb,ilr)
        jjlr=onWhichAtom(jjorb)
        jjproc=onWhichMPI(jjorb)
        korb=jjorb-isorb_par(jjproc)
        !if(iproc==0) write(*,'(a,8i8)') 'iorb, ilr, iiproc, jorb, jjorb, jjlr, jjproc, korb', iorb, ilr, iiproc, jorb, jjorb, jjlr, jjproc, korb
        if(iproc==jjproc) then
            !write(*,'(a,6i9)') 'iproc, jorb, jjorb, ilr, comom%overlaps(jorb,ilr), comom%olr(jorb,ilr)%norbinlr', iproc, jorb, jjorb, ilr, comom%overlaps(jorb,ilr), comom%olr(jorb,ilr)%norbinlr
            !write(*,'(3(a,i0),6x,a,4i8)') 'process ',iproc,' adds ',comom%olr(jorb,ilr)%norbinlr,' elements to position ',ist,'. iorb, ilr, jorb, jjorb', iorb, ilr, jorb, jjorb
            do i=1,comom%olr(jorb,ilr)%norbinlr
                ind=comom%olr(jorb,ilr)%indexInGlobal(i)
                comom%sendBuf(ist+i)=vec(ind,korb)
            end do
            ist=ist+comom%olr(jorb,ilr)%norbinlr
        end if
    end do
    ilrold=ilr
    iiprocold=iiproc
end do

if(ist/=comom%nsendBuf) then
    write(*,'(x,a,i0,a,2(2x,i0))') 'ERROR on process ',iproc,': ist/=comom%nsendBuf',ist,comom%nsendBuf
end if


end subroutine extractToOverlapregion





subroutine expandFromOverlapregion(iproc, nproc, isorb, norbp, orbs, onWhichAtom, comom, norbmax, noverlaps, vecOvrlp)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, isorb, norbp, norbmax, noverlaps
type(orbitals_data),intent(in):: orbs
integer,dimension(orbs%norb),intent(in):: onWhichAtom
type(p2pCommsOrthonormalityMatrix),intent(in):: comom
real(8),dimension(norbmax,noverlaps),intent(out):: vecOvrlp

! Local variables
integer:: ilrold, ist, iorb, iiorb, ilr, jorb, klr, korb, i, ind, ijorb


vecOvrlp=0.d0
ilrold=0
ist=0
ijorb=0
do iorb=1,norbp
    iiorb=isorb+iorb
    ilr=onWhichAtom(iiorb)
    if(ilr==ilrold) cycle
    do jorb=1,comom%noverlap(ilr)
        ijorb=ijorb+1
        klr=comom%olrForExpansion(1,jorb,ilr)
        korb=comom%olrForExpansion(2,jorb,ilr)
        do i=1,comom%olr(korb,klr)%norbinlr
            ind=comom%olr(korb,klr)%indexInGlobal(i)
            vecOvrlp(ind,ijorb)=comom%recvBuf(ist+i)
            !if(iproc==4) write(*,'(a,9i9)') 'iorb, iiorb, ilr, jorb, klr, korb, ist, i, ind', iorb, iiorb, ilr, jorb, klr, korb, ist, i, ind
        end do
        ist=ist+comom%olr(korb,klr)%norbinlr
    end do
    ilrold=ilr
end do

if(ist/=comom%nrecvBuf) then
    write(*,'(a,i0,a,2(2x,i0))') 'ERROR on process ',iproc,': ist/=comom%nrecvBuf',ist,comom%nrecvBuf
end if


end subroutine expandFromOverlapregion



subroutine calculateOverlap(iproc, nproc, nlr, norbmax, norbp, noverlaps, isorb, norb, comom, mlr, onWhichAtom, vec, vecOvrlp, newComm, ovrlp)
use module_base
use module_types
implicit none

! Calling arguments 
integer,intent(in):: iproc, nproc, nlr, norbmax, norbp, noverlaps, isorb, norb, newComm
type(p2pCommsOrthonormalityMatrix),intent(in):: comom
type(matrixLocalizationRegion),dimension(nlr),intent(in):: mlr
integer,dimension(norb),intent(in):: onWhichAtom
real(8),dimension(norbmax,norbp),intent(in):: vec
real(8),dimension(norbmax,noverlaps),intent(in):: vecOvrlp
real(8),dimension(norb,norb),intent(out):: ovrlp

! Local variables
integer:: ijorb, ilrold, ilr, iorb, iiorb, ncount, jjorb, jorb, ierr
real(8):: ddot

ovrlp=0.d0

ijorb=0
ilrold=0
do iorb=1,norbp
    iiorb=isorb+iorb
    ilr=onWhichAtom(iiorb)
    if(ilr==ilrold) then
         ! Put back index if we are in the same localization region, since then we can use the same vecOvrlp again.
         ijorb=ijorb-comom%noverlap(ilr) 
    end if
    ncount=mlr(ilr)%norbinlr
    do jorb=1,comom%noverlap(ilr)
        ijorb=ijorb+1
        jjorb=comom%overlaps(jorb,ilr)
        ovrlp(iiorb,jjorb)=ddot(ncount, vec(1,iorb), 1, vecOvrlp(1,ijorb), 1)
    end do
    ilrold=ilr
end do

call mpiallred(ovrlp(1,1), norb**2, mpi_sum, newComm, ierr)


end subroutine calculateOverlap




subroutine orthonormalLinearCombinations(iproc, nproc, nlr, norbmax, norbp, noverlaps, isorb, norb, comom, mlr, onWhichAtom, vecOvrlp, ovrlp, vec)
use module_base
use module_types
implicit none

! Calling arguments 
integer,intent(in):: iproc, nproc, nlr, norbmax, norbp, noverlaps, isorb, norb
type(p2pCommsOrthonormalityMatrix),intent(in):: comom
type(matrixLocalizationRegion),dimension(nlr),intent(in):: mlr
integer,dimension(norb),intent(in):: onWhichAtom
real(8),dimension(norbmax,noverlaps),intent(in):: vecOvrlp
real(8),dimension(norb,norb),intent(in):: ovrlp
real(8),dimension(norbmax,norbp),intent(inout):: vec

! Local variables
integer:: ijorb, ilrold, ilr, iorb, iiorb, ncount, jjorb, jorb, ierr, istat, iall
real(8):: ddot
real(8),dimension(:,:),allocatable:: vecTemp
character(len=*),parameter:: subname='orthonormalLinearCombinations'

allocate(vecTemp(norbmax,norbp), stat=istat)
call memocc(istat, vecTemp, 'vecTemp',subname)

call dcopy(norbmax*norbp, vec(1,1), 1, vecTemp(1,1), 1)

vec=0.d0

ijorb=0
ilrold=0
do iorb=1,norbp
    iiorb=isorb+iorb
    ilr=onWhichAtom(iiorb)
    if(ilr==ilrold) then
         ! Put back index if we are in the same localization region, since then we can use the same vecOvrlp again.
         ijorb=ijorb-comom%noverlap(ilr) 
    end if
    ncount=mlr(ilr)%norbinlr
    do jorb=1,comom%noverlap(ilr)
        ijorb=ijorb+1
        jjorb=comom%overlaps(jorb,ilr)
        call daxpy(ncount, ovrlp(jjorb,iiorb), vecOvrlp(1,ijorb), 1, vec(1,iorb), 1)
    end do
    ilrold=ilr
end do

iall=-product(shape(vecTemp))*kind(vecTemp)
deallocate(vecTemp, stat=istat)
call memocc(istat, iall, 'vecTemp', subname)

end subroutine orthonormalLinearCombinations




subroutine applyOrthoconstraintVectors(iproc, nproc, norb, norbmax, norbp, isorb, nlr, noverlaps, onWhichAtom, vecOvrlp, ovrlp, lagmat, comom, mlr, grad)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, norb, norbmax, norbp, isorb, nlr, noverlaps
integer,dimension(norb),intent(in):: onWhichAtom
real(8),dimension(norbmax,noverlaps),intent(in):: vecOvrlp
real(8),dimension(norb,norb),intent(in):: ovrlp, lagmat
type(p2pCommsOrthonormalityMatrix),intent(in):: comom
type(matrixLocalizationRegion),dimension(nlr),intent(in):: mlr
real(8),dimension(norbmax,norbp),intent(inout):: grad

! Local variables
integer:: info, iorb, ilrold, iiorb, jjorb, ilr, ncount, jorb, ijorb, istat, iall
real(8),dimension(:,:),allocatable:: ovrlp2, ovrlp_minus_one_lagmat, ovrlp_minus_one_lagmat_trans
character(len=*),parameter:: subname='applyOrthoconstarintVectors'


allocate(ovrlp_minus_one_lagmat(norb,norb), stat=istat)
call memocc(istat, ovrlp_minus_one_lagmat, 'ovrlp_minus_one_lagmat', subname)
allocate(ovrlp_minus_one_lagmat_trans(norb,norb), stat=istat)
call memocc(istat, ovrlp_minus_one_lagmat_trans, 'ovrlp_minus_one_lagmat_trans', subname)
allocate(ovrlp2(norb,norb), stat=istat)
call memocc(istat, ovrlp2, 'ovrlp2', subname)

call dcopy(norb**2, ovrlp(1,1), 1, ovrlp2(1,1), 1)
! Invert the overlap matrix
call dpotrf('l', norb, ovrlp2(1,1), norb, info)
if(info/=0) then
    write(*,'(x,a,i0)') 'ERROR in dpotrf, info=',info
    stop
end if
call dpotri('l', norb, ovrlp2(1,1), norb, info)
if(info/=0) then
    write(*,'(x,a,i0)') 'ERROR in dpotri, info=',info
    stop
end if


! Multiply the Lagrange multiplier matrix with S^-1/2.
! First fill the upper triangle.
do iorb=1,norb
    do jorb=1,iorb-1
        ovrlp2(jorb,iorb)=ovrlp2(iorb,jorb)
    end do
end do
call dgemm('n', 'n', norb, norb, norb, 1.d0, ovrlp2(1,1), norb, lagmat(1,1), norb, &
     0.d0, ovrlp_minus_one_lagmat(1,1), norb)
call dgemm('n', 't', norb, norb, norb, 1.d0, ovrlp2(1,1), norb, lagmat(1,1), norb, &
     0.d0, ovrlp_minus_one_lagmat_trans(1,1), norb)


ilrold=-1
ijorb=0
do iorb=1,norbp
    iiorb=isorb+iorb
    ilr=onWhichAtom(iiorb)
    if(ilr==ilrold) then
        ! Set back the index of lphiovrlp, since we again need the same orbitals.
        ijorb=ijorb-comom%noverlap(ilr)
    end if
    ncount=mlr(ilr)%norbinlr
    do jorb=1,comom%noverlap(ilr)
        ijorb=ijorb+1
        jjorb=comom%overlaps(jorb,ilr)
        call daxpy(ncount, -.5d0*ovrlp_minus_one_lagmat(jjorb,iiorb), vecOvrlp(1,ijorb), 1, grad(1,iorb), 1)
        call daxpy(ncount, -.5d0*ovrlp_minus_one_lagmat_trans(jjorb,iiorb), vecOvrlp(1,ijorb), 1, grad(1,iorb), 1)
    end do
    ilrold=ilr
end do


iall=-product(shape(ovrlp_minus_one_lagmat))*kind(ovrlp_minus_one_lagmat)
deallocate(ovrlp_minus_one_lagmat, stat=istat)
call memocc(istat, iall, 'ovrlp_minus_one_lagmat', subname)
iall=-product(shape(ovrlp_minus_one_lagmat_trans))*kind(ovrlp_minus_one_lagmat_trans)
deallocate(ovrlp_minus_one_lagmat_trans, stat=istat)
call memocc(istat, iall, 'ovrlp_minus_one_lagmat_trans', subname)
iall=-product(shape(ovrlp2))*kind(ovrlp2)
deallocate(ovrlp2, stat=istat)
call memocc(istat, iall, 'ovrlp2', subname)


end subroutine applyOrthoconstraintVectors




subroutine buildLinearCombinations(iproc, nproc, lzdig, lzd, input, coeff, lchi, lphi)
use module_base
use module_types
use module_interfaces, exceptThisOne => buildLinearCombinations
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(linear_zone_descriptors),intent(in):: lzdig, lzd
type(input_variables),intent(in):: input
real(8),dimension(lzdig%orbs%norb,lzd%orbs%norb),intent(in):: coeff
real(8),dimension(lzdig%orbs%npsidim),intent(in):: lchi
real(8),dimension(lzd%orbs%npsidim),intent(out):: lphi

! Local variables
integer:: tag, istat, iall, ist, jst, ilr, ilrold, iorb, iiorb, ncount, jorb, jjorb
type(overlapParameters):: op
type(p2pCommsOrthonormality):: comon
real(8),dimension(:),allocatable:: lchiovrlp
character(len=*),parameter:: subname='buildLinearCombinations'

tag=10000
call initCommsOrtho(iproc, nproc, lzdig, lzdig%orbs%inWhichLocreg, input, op, comon, tag)
allocate(lchiovrlp(op%ndim_lphiovrlp), stat=istat)
call memocc(istat, lchiovrlp, 'lchiovrlp',subname)

call allocateCommuncationBuffersOrtho(comon, subname)
call extractOrbital2(iproc, nproc, lzdig%orbs, lzdig%orbs%npsidim, lzdig%orbs%inWhichLocreg, lzdig, op, lchi, comon)
call postCommsOverlap(iproc, nproc, comon)
call gatherOrbitals2(iproc, nproc, comon)
call expandOrbital2(iproc, nproc, lzdig%orbs, input, lzdig%orbs%inWhichLocreg, lzdig, op, comon, lchiovrlp)
call deallocateCommuncationBuffersOrtho(comon, subname)



lphi=0.d0

ist=1
jst=1
ilrold=-1
do iorb=1,lzd%orbs%norbp
    iiorb=lzd%orbs%isorb+iorb
    ilr=lzd%orbs%inWhichLocreg(iiorb)
    if(ilr==ilrold) then
        ! Set back the index of lphiovrlp, since we again need the same orbitals.
        jst=jst-op%noverlaps(iiorb)*ncount
    end if
    ncount=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
    do jorb=1,op%noverlaps(iiorb)
        jjorb=op%overlaps(jorb,iiorb)
        !call daxpy(ncount, ovrlp(jjorb,iiorb), lphiovrlp(jst), 1, lphi(ist), 1)
        call daxpy(ncount, coeff(jjorb,iiorb), lchiovrlp(jst), 1, lphi(ist), 1)
        jst=jst+ncount
    end do

    ist=ist+ncount
    ilrold=ilr

end do





call deallocate_overlapParameters(op, subname)
call deallocate_p2pCommsOrthonormality(comon, subname)


iall=-product(shape(lchiovrlp))*kind(lchiovrlp)
deallocate(lchiovrlp, stat=istat)
call memocc(istat, iall, 'lchiovrlp', subname)



end subroutine buildLinearCombinations
