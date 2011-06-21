subroutine getLinearPsi(iproc, nproc, nspin, Glr, orbs, comms, at, lin, rxyz, rxyzParab, &
    nscatterarr, ngatherarr, nlpspd, proj, rhopot, GPU, input, pkernelseq, phi, psi, psit, updatePhi, &
    infoBasisFunctions, infoCoeff, itSCC, n3p, n3pi, n3d, irrzon, phnons, pkernel, pot_ion, rhocore, potxc, PSquiet, &
    i3s, i3xcsh, fion, fdisp, fxyz, eion, edisp, fnoise, ebs, coeff, lphi, radii_cf)
!
! Purpose:
! ========
!   This subroutine creates the orbitals psi out of a linear combination of localized basis functions
!   phi. To do so, it proceeds as follows:
!    1. Create the basis functions (with subroutine 'getLocalizedBasis')
!    2. Write the Hamiltonian in this new basis.
!    3. Diagonalize this Hamiltonian matrix.
!    4. Build the new linear combinations. 
!   The basis functions are localized by adding a confining quartic potential to the ordinary DFT 
!   Hamiltonian. There is no self consistency cycle for the potential, i.e. the basis functionsi
!   are optimized with a fixed potential.
!
! Calling arguments:
! ==================
!   Input arguments:
!   ----------------
!     iproc           process ID
!     nproc           total number of processes
!     nspin           npsin==1 -> closed shell; npsin==2 -> spin polarized
!     Glr             type describing the localization region
!     orbs            type describing the physical orbitals psi
!     comms           type containing the communication parameters for the physical orbitals psi
!     at              type containing the paraneters for the atoms
!     lin             type containing parameters for the linear version
!     rxyz            the atomic positions
!     rxyzParab       the center of the confinement potential (at the moment identical rxyz)
!     nscatterarr     ???
!     ngatherarr      ???
!     nlpsp           ???
!     proj            ???
!     rhopot          the charge density
!     GPU             parameters for GPUs
!     input           type containing some very general parameters
!     pkernelseq      ???
!     n3p             ???
!     itSCC           iteration in the self consistency cycle
!  Input/Output arguments
!  ---------------------
!     phi             the localized basis functions. It is assumed that they have been initialized
!                     somewhere else
!   Output arguments
!   ----------------
!     psi             the physical orbitals, which will be a linear combinations of the localized
!                     basis functions phi
!     psit            psi transposed
!     infoBasisFunctions  indicated wheter the basis functions converged to the specified limit (value is the
!                         number of iterations it took to converge) or whether the iteration stopped due to 
!                         the iteration limit (value is -1). This info is returned by 'getLocalizedBasis'
!     infoCoeff           the same as infoBasisFunctions, just for the coefficients. This value is returned
!                         by 'optimizeCoefficients'
!
use module_base
use module_types
use module_interfaces, exceptThisOne => getLinearPsi
use Poisson_Solver
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, nspin, n3p, n3pi, n3d, i3s, i3xcsh, itSCC
type(locreg_descriptors),intent(in):: Glr
type(orbitals_data),intent(in) :: orbs
type(communications_arrays),intent(in) :: comms
type(atoms_data),intent(in):: at
type(linearParameters),intent(inout):: lin
type(input_variables),intent(in):: input
real(8),dimension(3,at%nat),intent(in):: rxyz, fion, fdisp
real(8),dimension(3,at%nat),intent(inout):: rxyzParab
integer,dimension(0:nproc-1,4),intent(inout):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
integer,dimension(0:nproc-1,2),intent(inout):: ngatherarr
type(nonlocal_psp_descriptors),intent(in):: nlpspd
real(wp),dimension(nlpspd%nprojel),intent(inout):: proj
real(dp),dimension(max(Glr%d%n1i*Glr%d%n2i*n3p,1)*input%nspin),intent(inout) :: rhopot
type(GPU_pointers),intent(inout):: GPU
integer, dimension(lin%as%size_irrzon(1),lin%as%size_irrzon(2),lin%as%size_irrzon(3)),intent(in) :: irrzon 
real(dp), dimension(lin%as%size_phnons(1),lin%as%size_phnons(2),lin%as%size_phnons(3)),intent(in) :: phnons 
real(dp), dimension(lin%as%size_pkernel),intent(in):: pkernel
logical,intent(in):: updatePhi
real(wp), dimension(lin%as%size_pot_ion),intent(inout):: pot_ion
!real(wp), dimension(lin%as%size_rhocore):: rhocore 
real(wp), dimension(:),pointer,intent(in):: rhocore                  
real(wp), dimension(lin%as%size_potxc(1),lin%as%size_potxc(2),lin%as%size_potxc(3),lin%as%size_potxc(4)),intent(inout):: potxc
real(dp),dimension(:),pointer,intent(in):: pkernelseq
real(8),dimension(lin%lb%orbs%npsidim),intent(inout):: phi
real(8),dimension(orbs%npsidim),intent(out):: psi, psit
integer,intent(out):: infoBasisFunctions, infoCoeff
character(len=3),intent(in):: PSquiet
real(8),intent(out):: ebs
real(8),dimension(lin%lb%orbs%norb,orbs%norb),intent(in out):: coeff
real(8),dimension(3,at%nat),intent(out):: fxyz
real(8):: eion, edisp, fnoise
real(8),dimension(lin%lb%Lorbs%npsidim),intent(inout):: lphi
real(8),dimension(at%ntypes,3),intent(in):: radii_cf

! Local variables 
integer:: istat, iall, ind1, ind2, ldim, gdim, ilr, istr, nphibuff
real(8),dimension(:),allocatable:: hphi, eval, lhphi
real(8),dimension(:,:),allocatable:: HamSmall
real(8),dimension(:,:,:),allocatable:: matrixElements
real(8),dimension(:),pointer:: phiWork
real(8):: epot_sum, ekin_sum, eexctX, eproj_sum, ddot, trace, lastAlpha
real(wp),dimension(:),pointer:: potential 
character(len=*),parameter:: subname='getLinearPsi' 
logical:: withConfinement
type(workarr_sumrho):: w


real(8):: hxh, hyh, hzh, ehart, eexcu, vexcu, tt, energybs
integer:: iorb, jorb, it, istart, korb
character(len=11):: procName, orbNumber, orbName
character(len=30):: filename


integer:: ist, ierr

  ! Allocate the local arrays.  
  allocate(hphi(lin%lb%orbs%npsidim), stat=istat) 
  call memocc(istat, hphi, 'hphi', subname)
  allocate(matrixElements(lin%lb%orbs%norb,lin%lb%orbs%norb,2), stat=istat)
  call memocc(istat, matrixElements, 'matrixElements', subname)
  allocate(HamSmall(lin%lb%orbs%norb,lin%lb%orbs%norb), stat=istat)
  call memocc(istat, HamSmall, 'HamSmall', subname)
  allocate(phiWork(max(size(phi),size(psi))), stat=istat)
  call memocc(istat, phiWork, 'phiWork', subname)
  allocate(eval(lin%orbs%norb), stat=istat)
  call memocc(istat, eval, 'eval', subname)
  

  ! This is a flag whether the basis functions shall be updated.
  if(updatePhi) then
      if(lin%useDerivativeBasisFunctions) then
          call dcopy(lin%orbs%npsidim, lin%phiRestart(1), 1, phi(1), 1)
      end if
      ! Optimize the localized basis functions by minimizing the trace of <phi|H|phi>.
      call getLocalizedBasis(iproc, nproc, at, orbs, Glr, input, lin, rxyz, nspin, nlpspd, proj, &
          nscatterarr, ngatherarr, rhopot, GPU, pkernelseq, phi, trace, rxyzParab, &
          itSCC, lastAlpha, infoBasisFunctions, radii_cf)
  end if

!!! 3D plot of the basis functions
!!write(procName,'(i0)') iproc
!!istart=1
!!do iorb=1,lin%orbs%norbp
!!  write(orbNumber,'(i0)') iorb
!!  orbName='orb_'//trim(procName)//'_'//trim(orbNumber)
!!  call plot_wfSquare_cube(orbName, at, Glr, input%hx, input%hy, input%hz, rxyz, phi(istart), 'comment   ')
!!  istart=istart+Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
!!end do


  if(lin%useDerivativeBasisFunctions) then
      ! Create the derivative basis functions.
      nphibuff=0
      do iorb=1,lin%orbs%norbp
          nphibuff = nphibuff + Glr%wfd%nvctr_c + 7*Glr%wfd%nvctr_f
      end do
      call dcopy(nphibuff, phi(1), 1, lin%phiRestart(1), 1)
      call getDerivativeBasisFunctions(iproc, nproc, input%hx, Glr, lin, nphibuff, lin%phiRestart, phi)

      ! Orthonormalize
      call transpose_v(iproc, nproc, lin%lb%orbs, Glr%wfd, lin%lb%comms, phi, work=phiWork)
      call orthonormalizeOnlyDerivatives(iproc, nproc, lin, phi)
      call untranspose_v(iproc, nproc, lin%lb%orbs, Glr%wfd, lin%lb%comms, phi, work=phiWork)
  end if

  ! Transform the global phi to the local phi
  ! This part will not be needed if we really have O(N)
  ind1=1
  ind2=1
  do iorb=1,lin%lb%orbs%norbp
      ilr = lin%lb%onWhichAtom(iorb)
      ldim=lin%Llr(ilr)%wfd%nvctr_c+7*lin%Llr(ilr)%wfd%nvctr_f
      gdim=Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
      call psi_to_locreg2(iproc, nproc, ldim, gdim, lin%Llr(ilr), Glr, phi(ind1), lphi(ind2))
      ind1=ind1+Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
      ind2=ind2+lin%Llr(ilr)%wfd%nvctr_c+7*lin%Llr(ilr)%wfd%nvctr_f
  end do
  ! Transfom all orbitals to real space
  ist=1
  istr=1
  do iorb=1,lin%lb%orbs%norbp
      ilr=lin%lb%onWhichAtom(iorb)
      call initialize_work_arrays_sumrho(lin%Llr(ilr), w)
      !call daub_to_isf(lin%Llr(ilr), w, lphi(ist), lphir(istr))
      call daub_to_isf(lin%Llr(ilr), w, lphi(ist), lin%comsr%sendBuf(istr))
      call deallocate_work_arrays_sumrho(w)
      ist = ist + lin%Llr(ilr)%wfd%nvctr_c + 7*lin%Llr(ilr)%wfd%nvctr_f
      istr = istr + lin%Llr(ilr)%d%n1i*lin%Llr(ilr)%d%n2i*lin%Llr(ilr)%d%n3i
  end do
  
  ! Post the MPI messages for the communication of sumrho. Since we use non blocking point
  ! to point communication, the program will continue immediately.
  call postCommunicationSumrho2(iproc, nproc, lin, lin%comsr%sendBuf, lin%comsr%recvBuf)
  

  if(iproc==0) write(*,'(x,a)') '----------------------------------- Determination of the orbitals in this new basis.'


  if(.not.updatePhi) then
      ! Otherwise the potential is gathered in getLocalizedBasis.
      call gatherPotential(iproc, nproc, lin%comgp)
  end if

  if(trim(lin%getCoeff)=='min') then

      !!!! THIS IS THE ORIGINAL
      !!call full_local_potential(iproc,nproc,Glr%d%n1i*Glr%d%n2i*n3p,Glr%d%n1i*Glr%d%n2i*Glr%d%n3i,input%nspin,&
      !!     lin%orbs%norb,lin%orbs%norbp,ngatherarr,rhopot,potential)
      !!call HamiltonianApplication(iproc,nproc,at,lin%lb%orbs,input%hx,input%hy,input%hz,rxyz,&
      !!     nlpspd,proj,Glr,ngatherarr,potential,&
      !!     phi(1),hphi(1),ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU,pkernel=pkernelseq)
      !!!deallocate potential
      !!call free_full_potential(nproc,potential,subname)

      ! Transform the global phi to the local phi
      ! This part will not be needed if we really have O(N)
      allocate(lhphi(lin%lb%Lorbs%npsidim), stat=istat)
      call memocc(istat, lhphi, 'lhphi', subname)

      ind1=1
      ind2=1
      !lphi=0.d0
      do iorb=1,lin%lb%orbs%norbp
          ilr = lin%lb%onWhichAtom(iorb)
          ldim=lin%Llr(ilr)%wfd%nvctr_c+7*lin%Llr(ilr)%wfd%nvctr_f
          gdim=Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
          call psi_to_locreg2(iproc, nproc, ldim, gdim, lin%Llr(ilr), Glr, phi(ind1), lphi(ind2))
          ind1=ind1+Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
          ind2=ind2+lin%Llr(ilr)%wfd%nvctr_c+7*lin%Llr(ilr)%wfd%nvctr_f
      end do

      withConfinement=.false.
      call HamiltonianApplicationConfinement2(input, iproc, nproc, at, lin%lzd, lin, input%hx, input%hy, input%hz, rxyz,&
           proj, ngatherarr, lin%comgp%nrecvBuf, lin%comgp%recvBuf, lphi, lhphi, &
           ekin_sum, epot_sum, eexctX, eproj_sum, nspin, GPU, radii_cf, lin%comgp, lin%onWhichAtom, withConfinement, &
           pkernel=pkernelseq)
      ind1=1
      ind2=1
      hphi=0.d0
      do iorb=1,lin%lb%orbs%norbp
          ilr = lin%lb%onWhichAtom(iorb)
          ldim=lin%Llr(ilr)%wfd%nvctr_c+7*lin%Llr(ilr)%wfd%nvctr_f
          gdim=Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
          call Lpsi_to_global2(iproc, nproc, ldim, gdim, lin%orbs%norb, lin%orbs%nspinor, input%nspin, Glr, lin%Llr(ilr), lhphi(ind2), hphi(ind1))
          ind1=ind1+Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
          ind2=ind2+lin%Llr(ilr)%wfd%nvctr_c+7*lin%Llr(ilr)%wfd%nvctr_f
      end do

      iall=-product(shape(lhphi))*kind(lhphi)
      deallocate(lhphi, stat=istat)
      call memocc(istat, iall, 'lhphi', subname)

      if(iproc==0) write(*,'(x,a)', advance='no') 'done.'
      

      ! Calculate the matrix elements <phi|H|phi>.
      call getMatrixElements(iproc, nproc, Glr, lin%lb%orbs, lin%lb%comms, phi, hphi, matrixElements)

          ! Calculate the coefficients which minimize the band structure energy
          ! ebs = \sum_i \sum_{k,l} c_{ik}*c_{il}*<phi_k|H|phi_l>
          ! for the given basis functions.
          !if(.not. updatePhi) call optimizeCoefficients(iproc, orbs, lin, nspin, matrixElements, coeff, infoCoeff)
          call optimizeCoefficients(iproc, orbs, lin, nspin, matrixElements, coeff, infoCoeff)
      !!call diagonalizeHamiltonian(iproc, nproc, lin%orbs, matrixElements(1,1,1), eval)
      !!call dcopy(lin%lb%orbs%norb*orbs%norb, matrixElements(1,1,1), 1, coeff(1,1), 1)

      ! Calculate the band structure energy with matrixElements.
      ebs=0.d0
      do iorb=1,orbs%norb
          do jorb=1,lin%lb%orbs%norb
              do korb=1,lin%lb%orbs%norb
                  ebs = ebs + coeff(jorb,iorb)*coeff(korb,iorb)*matrixElements(korb,jorb,1)
              end do
          end do
      end do
      ! If closed shell multiply by two
      if(input%nspin==1) ebs=2.d0*ebs

  else if(trim(lin%getCoeff)=='diag') then
  
      !allocate the potential in the full box
      call full_local_potential(iproc,nproc,Glr%d%n1i*Glr%d%n2i*n3p,Glr%d%n1i*Glr%d%n2i*Glr%d%n3i,input%nspin,&
           lin%orbs%norb,lin%orbs%norbp,ngatherarr,rhopot,potential)
      
      call HamiltonianApplication(iproc,nproc,at,lin%orbs,input%hx,input%hy,input%hz,rxyz,&
           nlpspd,proj,Glr,ngatherarr,potential,&
           phi(1),hphi(1),ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU,pkernel=pkernelseq)
      if(iproc==0) write(*,'(x,a)', advance='no') 'done.'
      
      !deallocate potential
      call free_full_potential(nproc,potential,subname)

  end if
  
  
  call transpose_v(iproc, nproc, lin%lb%orbs, Glr%wfd, lin%lb%comms, phi, work=phiWork)
  call transpose_v(iproc, nproc, lin%lb%orbs, Glr%wfd, lin%lb%comms, hphi, work=phiWork)
      
  
  if(trim(lin%getCoeff)=='diag') then
      call transformHam(iproc, nproc, lin%orbs, lin%comms, phi, hphi, HamSmall)
      if(iproc==0) write(*,'(a)', advance='no') ' Diagonalization... '
      call diagonalizeHamiltonian(iproc, nproc, lin%orbs, HamSmall, eval)
      call dcopy(lin%orbs%norb*orbs%norb, HamSmall(1,1), 1, coeff(1,1), 1)
      if(iproc==0) write(*,'(a)') 'done.'
  end if

  
  if(iproc==0) then
      write(*,'(x,a)', advance='no') '------------------------------------- Building linear combinations... '
      
  end if
  if(trim(lin%getCoeff)=='diag') then
      call buildWavefunction(iproc, nproc, orbs, lin%orbs, comms, lin%comms, phi, psi, HamSmall)
  else if(trim(lin%getCoeff)=='min') then
      call buildWavefunctionModified(iproc, nproc, orbs, lin%lb%orbs, comms, lin%lb%comms, phi, psi, coeff)
  else
      if(iproc==0) write(*,'(a,a,a)') "ERROR: lin%getCoeff can have the values 'diag' or 'min' , &
          & but we found '", lin%getCoeff, "'."
      stop
  end if

  
  call dcopy(orbs%npsidim, psi, 1, psit, 1)
  if(iproc==0) write(*,'(a)') 'done.'
  
  
  if(.not.lin%useDerivativeBasisFunctions) then
      call untranspose_v(iproc, nproc, lin%lb%orbs, Glr%wfd, lin%lb%comms, phi, work=phiWork)
      call untranspose_v(iproc, nproc, orbs, Glr%wfd, comms, psi, work=phiWork)
  else
      call untranspose_v(iproc, nproc, orbs, Glr%wfd, comms, psi, work=phiWork)
  end if


  !!!!! Improve the coefficients -- EXPERIMENTAL
  !!!if(updatePhi) then
  !!!    do iorb=1,orbs%norb
  !!!        do jorb=1,lin%lb%orbs%norb
  !!!            tt=0.d0
  !!!            do korb=1,lin%lb%orbs%norb
  !!!                tt=tt+coeff(korb,iorb)*matrixElements(jorb,korb,1)
  !!!            end do
  !!!            coeff(jorb,iorb) = coeff(jorb,iorb)-1.d-1*tt
  !!!        end do 
  !!!    end do
  !!!end if
  !!!!!


  
  iall=-product(shape(HamSmall))*kind(HamSmall)
  deallocate(HamSmall, stat=istat)
  call memocc(istat, iall, 'HamSmall', subname)

  iall=-product(shape(hphi))*kind(hphi)
  deallocate(hphi, stat=istat)
  call memocc(istat, iall, 'hphi', subname)

  iall=-product(shape(phiWork))*kind(phiWork)
  deallocate(phiWork, stat=istat)
  call memocc(istat, iall, 'phiWork', subname)

  iall=-product(shape(matrixElements))*kind(matrixElements)
  deallocate(matrixElements, stat=istat)
  call memocc(istat, iall, 'matrixElements', subname)
  
  iall=-product(shape(eval))*kind(eval)
  deallocate(eval, stat=istat)
  call memocc(istat, iall, 'eval', subname)




end subroutine getLinearPsi








subroutine getLocalizedBasis(iproc, nproc, at, orbs, Glr, input, lin, rxyz, nspin, nlpspd, &
    proj, nscatterarr, ngatherarr, rhopot, GPU, pkernelseq, phi, trH, rxyzParabola, &
    itScc, lastAlpha, infoBasisFunctions, radii_cf)
!
! Purpose:
! ========
!   Calculates the localized basis functions phi. These basis functions are obtained by adding a
!   quartic potential centered on the atoms to the ordinary Hamiltonian. The eigenfunctions are then
!   determined by minimizing the trace until the gradient norm is below the convergence criterion.
!
! Calling arguments:
! ==================
!   Input arguments:
!   ----------------
!     iproc           process ID
!     nproc           total number of processes
!     at              type containing the paraneters for the atoms
!     orbs            type describing the physical orbitals psi
!     Glr             type describing the localization region
!     input           type containing some very general parameters
!     lin             type containing parameters for the linear version
!     rxyz            the atomic positions
!     nspin           npsin==1 -> closed shell; npsin==2 -> spin polarized
!     nlpsp           ???
!     proj            ???
!     nscatterarr     ???
!     ngatherarr      ???
!     rhopot          the charge density
!     GPU             parameters for GPUs
!     pkernelseq      ???
!     rxyzParab       the center of the confinement potential (at the moment identical rxyz)
!     n3p             ???
!     itSCC           iteration in the self consistency cycle
!  Input/Output arguments
!  ---------------------
!     phi             the localized basis functions. It is assumed that they have been initialized
!                     somewhere else
!   Output arguments
!   ----------------
!     hphi            the modified Hamiltonian applied to phi
!     trH             the trace of the Hamiltonian
!     infoBasisFunctions  indicates wheter the basis functions converged to the specified limit (value is 0)
!                         or whether the iteration stopped due to the iteration limit (value is -1). This info
!                         is returned by 'getLocalizedBasis'
!
! Calling arguments:
!   Input arguments
!   Output arguments
!    phi   the localized basis functions
!
use module_base
use module_types
use module_interfaces, except_this_one => getLocalizedBasis
!  use Poisson_Solver
!use allocModule
implicit none

! Calling arguments
integer:: iproc, nproc, infoBasisFunctions, itSCC
type(atoms_data), intent(in) :: at
type(orbitals_data):: orbs
type(locreg_descriptors), intent(in) :: Glr
type(input_variables):: input
type(linearParameters):: lin
real(8),dimension(3,at%nat):: rxyz, rxyzParabola
integer:: nspin
type(nonlocal_psp_descriptors), intent(in) :: nlpspd
real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
real(dp), dimension(*), intent(inout) :: rhopot
type(GPU_pointers), intent(inout) :: GPU
real(dp), dimension(:), pointer :: pkernelseq
real(8),dimension(lin%orbs%npsidim):: phi
real(8):: trH, lastAlpha
real(8),dimension(at%ntypes,3),intent(in):: radii_cf

! Local variables
real(8) ::epot_sum, ekin_sum, eexctX, eproj_sum, evalmax, eval_zero
real(8):: tt, ddot, fnrm, fnrmMax, meanAlpha, gnrm, gnrm_zero, gnrmMax, t1, t2, time
integer:: iorb, icountSDSatur, icountSwitch, idsx, icountDIISFailureTot, icountDIISFailureCons, itBest
integer:: istat, istart, ierr, ii, it, nbasisPerAtForDebug, ncong, iall, nvctrp, nit, ind1, ind2
integer:: ldim, gdim, ilr
real(8),dimension(:),allocatable:: hphi, hphiold, alpha, fnrmOldArr, lagMatDiag, alphaDIIS, lphi, lhphi
real(8),dimension(:,:),allocatable:: HamSmall, fnrmArr, fnrmOvrlpArr
real(8),dimension(:),pointer:: phiWork
logical:: quiet, allowDIIS, startWithSD, adapt, withConfinement
character(len=*),parameter:: subname='getLocalizedBasis'
character(len=1):: message
type(diis_objects):: diisLIN

real(8),dimension(:,:),allocatable:: lagMult
real(8),dimension(:,:,:),allocatable:: ovrlp
integer:: jstart, jorb

allocate(lagMult(lin%orbs%norb,lin%orbs%norb), stat=istat)
lagMult=1.d-1
allocate(ovrlp(lin%orbs%norb,lin%orbs%norb,2), stat=istat)

allocate(lagMatDiag(lin%orbs%norb), stat=istat)

  ! Allocate all local arrays
  call allocateLocalArrays()
  
  ! Initialize the DIIS parameters 
  icountSDSatur=0
  icountSwitch=0
  icountDIISFailureTot=0
  icountDIISFailureCons=0
  call initializeDIISParameters(lin%DIISHistMax)
  if(lin%startWithSD) then
      allowDIIS=.false.
      diisLIN%switchSD=.false.
      startWithSD=.true.
  else
      allowDIIS=.true.
      startWithSD=.false.
  end if
  
  if(iproc==0) write(*,'(x,a)') '======================== Creation of the basis functions... ========================'

  ! Assign the step size for SD iterations.
  alpha=lin%alphaSD
  alphaDIIS=lin%alphaDIIS
  adapt=.false.

  ! Cut off outside localization region -- experimental
  call cutoffOutsideLocreg(iproc, nproc, Glr, at, input, lin, rxyz, phi)
  ! Transpose phi
  call transpose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, phi, work=phiWork)

  if(itSCC==1) then
      nit=lin%nItBasisFirst
  else
      nit=lin%nItBasis
  end if

  ! Gather the potential
  call gatherPotential(iproc, nproc, lin%comgp)

  iterLoop: do it=1,nit
      fnrmMax=0.d0
      fnrm=0.d0
  
      if (iproc==0) then
          write( *,'(1x,a,i0)') repeat('-',77 - int(log(real(it))/log(10.))) // ' iter=', it
      endif
  
      ! Orthonormalize the orbitals.
      if(iproc==0) then
          write(*,'(x,a)', advance='no') 'Orthonormalization... '
      end if
      call orthogonalize(iproc, nproc, lin%orbs, lin%comms, Glr%wfd, phi, input)

      ! Untranspose phi
      call untranspose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, phi, work=phiWork)
  
  
      ! Calculate the unconstrained gradient.
      if(iproc==0) then
          write(*,'(a)', advance='no') 'Hamiltonian application... '
      end if

      !!!! THIS IS THE ORIGINAL
      !!call HamiltonianApplicationConfinement(iproc,nproc,at,lin%orbs,lin,input%hx,input%hy,input%hz,rxyz,&
      !!     nlpspd,proj,Glr,ngatherarr,Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2),&
      !!     rhopot(1),&
      !!     phi(1),hphi(1),ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU, rxyzParabola, lin%onWhichAtom, pkernel=pkernelseq)

      ! Transform the global phi to the local phi
      ! This part will not be needed if we really have O(N)
      ind1=1
      ind2=1
      do iorb=1,lin%orbs%norbp
          ilr = lin%onWhichAtom(iorb)
          ldim=lin%Llr(ilr)%wfd%nvctr_c+7*lin%Llr(ilr)%wfd%nvctr_f
          gdim=Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
          call psi_to_locreg2(iproc, nproc, ldim, gdim, lin%Llr(ilr), Glr, phi(ind1), lphi(ind2))
          ind1=ind1+Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
          ind2=ind2+lin%Llr(ilr)%wfd%nvctr_c+7*lin%Llr(ilr)%wfd%nvctr_f
      end do
      withConfinement=.true.
if(it==1) then
    do iall=1,size(lin%comgp%recvBuf)
        !write(8000+iproc,*) iall, lin%comgp%recvBuf(iall)
    end do
end if
      call HamiltonianApplicationConfinement2(input, iproc, nproc, at, lin%lzd, lin, input%hx, input%hy, input%hz, rxyz,&
           proj, ngatherarr, lin%comgp%nrecvBuf, lin%comgp%recvBuf, lphi, lhphi, &
           ekin_sum, epot_sum, eexctX, eproj_sum, nspin, GPU, radii_cf, lin%comgp, lin%onWhichAtom, withConfinement, &
           pkernel=pkernelseq)
      ind1=1
      ind2=1
      hphi=0.d0
      do iorb=1,lin%orbs%norbp
          ilr = lin%onWhichAtom(iorb)
          ldim=lin%Llr(ilr)%wfd%nvctr_c+7*lin%Llr(ilr)%wfd%nvctr_f
          gdim=Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
          call Lpsi_to_global2(iproc, nproc, ldim, gdim, lin%orbs%norb, lin%orbs%nspinor, input%nspin, Glr, lin%Llr(ilr), lhphi(ind2), hphi(ind1))
          ind1=ind1+Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
          ind2=ind2+lin%Llr(ilr)%wfd%nvctr_c+7*lin%Llr(ilr)%wfd%nvctr_f
      end do
!!write(*,'(a,i5,2es16.7)') 'iproc, ekin_sum, epot_sum', iproc, ekin_sum, epot_sum
!!do iorb=1,lin%orbs%npsidim
!!  write(30000+iproc*1000,*) iorb, hphi(iorb)
!!end do
!!call mpi_barrier(mpi_comm_world, ierr)
!!stop
!!! Plot the gradients
!!call plotOrbitals(iproc, lin%orbs, Glr, hphi, at%nat, rxyz, lin%onWhichAtom, .5d0*input%hx, &
!!    .5d0*input%hy, .5d0*input%hz, 500+it)
  
  
      ! Apply the orthoconstraint to the gradient. This subroutine also calculates the trace trH.
      if(iproc==0) then
          write(*,'(a)', advance='no') 'orthoconstraint... '
      end if
      call transpose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, hphi, work=phiWork)
      call transpose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, phi, work=phiWork)
      call orthoconstraintNotSymmetric(iproc, nproc, lin%orbs, lin%comms, Glr%wfd, phi, hphi, trH, lagMatDiag)
  
  
      ! Calculate the norm of the gradient (fnrmArr) and determine the angle between the current gradient and that
      ! of the previous iteration (fnrmOvrlpArr).
      nvctrp=lin%comms%nvctr_par(iproc,1) ! 1 for k-point
      istart=1
      do iorb=1,lin%orbs%norb
          if(it>1) fnrmOvrlpArr(iorb,2)=ddot(nvctrp*orbs%nspinor, hphi(istart), 1, hphiold(istart), 1)
          fnrmArr(iorb,2)=ddot(nvctrp*orbs%nspinor, hphi(istart), 1, hphi(istart), 1)
          istart=istart+nvctrp*orbs%nspinor
      end do
      !in case mpiallred is used
      !call mpiallred(fnrmArr(1,1),lin%orbs%norb,mpi_sum,mpi_comm_world, ierr) 
      call mpi_allreduce(fnrmArr(1,2), fnrmArr(1,1), lin%orbs%norb, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
      call mpi_allreduce(fnrmOvrlpArr(1,2), fnrmOvrlpArr(1,1), lin%orbs%norb, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
  
      ! Keep the gradient for the next iteration.
      if(it>1) then
          call dcopy(lin%orbs%norb, fnrmArr(1,1), 1, fnrmOldArr(1), 1)
      end if
  
      ! Determine the gradient norm and its maximal component. In addition, adapt the
      ! step size for the steepest descent minimization (depending on the angle 
      ! between the current gradient and the one from the previous iteration).
      ! This is of course only necessary if we are using steepest descent and not DIIS.
      do iorb=1,lin%orbs%norb
          fnrm=fnrm+fnrmArr(iorb,1)
          if(fnrmArr(iorb,1)>fnrmMax) fnrmMax=fnrmArr(iorb,1)
          if(it>1 .and. diisLIN%idsx==0 .and. .not.diisLIN%switchSD) then
          ! Adapt step size for the steepest descent minimization.
              tt=fnrmOvrlpArr(iorb,1)/sqrt(fnrmArr(iorb,1)*fnrmOldArr(iorb))
              !if(tt>.7d0) then
              if(tt>.9d0) then
                  alpha(iorb)=alpha(iorb)*1.05d0
              else
                  alpha(iorb)=alpha(iorb)*.5d0
              end if
          end if
      end do
      fnrm=sqrt(fnrm)
      fnrmMax=sqrt(fnrmMax)
      ! Copy the gradient (will be used in the next iteration to adapt the step size).
      call dcopy(lin%orbs%norb*nvctrp*orbs%nspinor, hphi(1), 1, hphiold(1), 1)
  
      ! Untranspose hphi.
      call untranspose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, hphi, work=phiWork)

  
      ! Adapt the preconditioning constant
      if(fnrmMax<1.d-9) then
          if(.not.adapt .and. iproc==0) then
              write(*,'(x,a)') 'Adapting the preconditioning constant from now on'
              adapt=.true.
          end if
          do iorb=1,lin%orbs%norb
              lin%orbs%eval(iorb)=lagMatDiag(iorb)
          end do
      end if
  

      ! Precondition the gradient
      if(iproc==0) then
          write(*,'(a)') 'preconditioning. '
      end if
      gnrm=1.d3 ; gnrm_zero=1.d3
      call cpu_time(t1)
      !!!! THIS IS THE ORIGINAL
      !!call choosePreconditioner(iproc, nproc, lin%orbs, lin, Glr, input%hx, input%hy, input%hz, &
      !!    lin%nItPrecond, hphi, at%nat, rxyz, at, it)

      evalmax=lin%orbs%eval(lin%orbs%isorb+1)
      do iorb=1,lin%orbs%norbp
        evalmax=max(lin%orbs%eval(lin%orbs%isorb+iorb),evalmax)
      enddo
      call MPI_ALLREDUCE(evalmax,eval_zero,1,mpidtypd,&
           MPI_MAX,MPI_COMM_WORLD,ierr)

      ! Transform the global phi to the local phi
      ! This part will not be needed if we really have O(N)
      ind1=1
      ind2=1
      do iorb=1,lin%orbs%norbp
          ilr = lin%onWhichAtom(iorb)
          ldim=lin%Llr(ilr)%wfd%nvctr_c+7*lin%Llr(ilr)%wfd%nvctr_f
          gdim=Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
          call psi_to_locreg2(iproc, nproc, ldim, gdim, lin%Llr(ilr), Glr, hphi(ind1), lhphi(ind2))
          call choosePreconditioner2(iproc, nproc, lin%orbs, lin, lin%Llr(ilr), input%hx, input%hy, input%hz, &
              lin%nItPrecond, lhphi(ind2), at%nat, rxyz, at, it, iorb, eval_zero)
          call razero(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f, hphi(ind1))
          call Lpsi_to_global2(iproc, nproc, ldim, gdim, lin%orbs%norb, lin%orbs%nspinor, input%nspin, Glr, lin%Llr(ilr), lhphi(ind2), hphi(ind1))
          ind1=ind1+Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
          ind2=ind2+lin%Llr(ilr)%wfd%nvctr_c+7*lin%Llr(ilr)%wfd%nvctr_f
      end do
      call cpu_time(t2)
      time=t2-t1
      call mpiallred(time, 1, mpi_sum, mpi_comm_world, ierr)
      if(iproc==0) write(*,'(x,a,es10.3)') 'time for preconditioning:', time/dble(nproc)

!!do iorb=1,lin%orbs%npsidim
!!  write(40000+iproc*1000,*) iorb, hphi(iorb)
!!end do
!!call mpi_barrier(mpi_comm_world, ierr)
!!stop

      ! Determine the mean step size for steepest descent iterations.
      tt=sum(alpha)
      meanAlpha=tt/dble(lin%orbs%norb)
  
      ! Write some informations to the screen.
      if(iproc==0) write(*,'(x,a,i6,2es15.7,f17.10)') 'iter, fnrm, fnrmMax, trace', it, fnrm, fnrmMax, trH
      if(iproc==0) write(1000,'(i6,2es15.7,f18.10,es12.4)') it, fnrm, fnrmMax, trH, meanAlpha
      if(fnrmMax<lin%convCrit .or. it>=nit) then
          if(it>=nit) then
              if(iproc==0) write(*,'(x,a,i0,a)') 'WARNING: not converged within ', it, &
                  ' iterations! Exiting loop due to limitations of iterations.'
              if(iproc==0) write(*,'(x,a,2es15.7,f12.7)') 'Final values for fnrm, fnrmMax, trace: ', fnrm, fnrmMax, trH
              infoBasisFunctions=-1
          else
              if(iproc==0) then
                  write(*,'(x,a,i0,a,2es15.7,f12.7)') 'converged in ', it, ' iterations.'
                  write (*,'(x,a,2es15.7,f12.7)') 'Final values for fnrm, fnrmMax, trace: ', fnrm, fnrmMax, trH
              end if
              infoBasisFunctions=it
          end if
          if(iproc==0) write(*,'(x,a)') '============================= Basis functions created. ============================='
          call untranspose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, phi, work=phiWork)
          if(lin%plotBasisFunctions) then
              call plotOrbitals(iproc, lin%orbs, Glr, phi, at%nat, rxyz, lin%onWhichAtom, .5d0*input%hx, &
                  .5d0*input%hy, .5d0*input%hz, 1)
          end if
          exit iterLoop
      end if
  
  
      call DIISorSD()
      if(iproc==0) then
          if(diisLIN%idsx>0) then
              write(*,'(x,3(a,i0))') 'DIIS informations: history length=',diisLIN%idsx, ', consecutive failures=', &
                  icountDIISFailureCons, ', total failures=', icountDIISFailureTot
          else
              if(allowDIIS) then
                  message='y'
              else
                  message='n'
              end if
              write(*,'(x,a,es9.3,a,i0,a,a)') 'steepest descent informations: mean alpha=', meanAlpha, &
              ', consecutive successes=', icountSDSatur, ', DIIS=', message
          end if
      end if
      if(.not. diisLIN%switchSD) call improveOrbitals()
  
     ! Flush the standard output
      call flush(6) 
  end do iterLoop

  ! Store the mean alpha.
  lastAlpha=meanAlpha

  call deallocateLocalArrays()

contains

    subroutine initializeDIISParameters(idsxHere)
    ! Purpose:
    ! ========
    !   Initializes all parameters needed for the DIIS procedure.
    !
    ! Calling arguments
    !   idsx    DIIS history length
    !
    implicit none
    
    ! Calling arguments
    integer:: idsxHere

      diisLIN%switchSD=.false.
      diisLIN%idiistol=0
      diisLIN%mids=1
      diisLIN%ids=0
      diisLIN%idsx=idsxHere
      diisLIN%energy_min=1.d10
      diisLIN%energy_old=1.d10
      diisLIN%energy=1.d10
      diisLIN%alpha=2.d0
      call allocate_diis_objects(diisLIN%idsx, lin%orbs%npsidim, 1, diisLIN, subname) ! 1 for k-points

    end subroutine initializeDIISParameters


    subroutine DIISorSD()
    !
    ! Purpose:
    ! ========
    !   This subroutine decides whether one should use DIIS or variable step size
    !   steepest descent to improve the orbitals. In the beginning we start with DIIS
    !   with history length lin%DIISHistMax. If DIIS becomes unstable, we switch to
    !   steepest descent. If the steepest descent iterations are successful, we switch
    !   back to DIIS, but decrease the DIIS history length by one. However the DIIS
    !   history length is limited to be larger or equal than lin%DIISHistMin.
    !


      ! First there are some checks whether the force is small enough to allow DIIS.

      ! Decide whether the force is small eneough to allow DIIS
      if(fnrmMax<lin%startDIIS .and. .not.allowDIIS) then
          allowDIIS=.true.
          if(iproc==0) write(*,'(x,a)') 'The force is small enough to allow DIIS.'
          ! This is to get the correct DIIS history 
          ! (it is chosen as max(lin%DIISHistMin,lin%DIISHistMax-icountSwitch).
          icountSwitch=icountSwitch-1
      else if(fnrmMax>lin%startDIIS .and. allowDIIS) then
          allowDIIS=.false.
          if(iproc==0) write(*,'(x,a)') 'The force is too large to allow DIIS.'
      end if    

      ! Switch to SD if the flag indicating that we should start with SD is true.
      ! If this is the case, this flag is set to false, since this flag concerns only the beginning.
      if(startWithSD .and. diisLIN%idsx>0) then
          call deallocate_diis_objects(diisLIN, subname)
          diisLIN%idsx=0
          diisLIN%switchSD=.false.
          startWithSD=.false.
      end if

      ! Decide whether we should switch from DIIS to SD in case we are using DIIS and it 
      ! is not allowed.
      if(.not.startWithSD .and. .not.allowDIIS .and. diisLIN%idsx>0) then
          if(iproc==0) write(*,'(x,a,es10.3)') 'The force is too large, switch to SD with stepsize', alpha(1)
          call deallocate_diis_objects(diisLIN, subname)
          diisLIN%idsx=0
          diisLIN%switchSD=.true.
      end if

      ! If we swicthed to SD in the previous iteration, reset this flag.
      if(diisLIN%switchSD) diisLIN%switchSD=.false.

      ! Now come some checks whether the trace is descreasing or not. This further decides
      ! whether we should use DIIS or SD.

      ! Determine wheter the trace is decreasing (as it should) or increasing.
      ! This is done by comparing the current value with diisLIN%energy_min, which is
      ! the minimal value of the trace so far.
      if(trH<=diisLIN%energy_min) then
          ! Everything ok
          diisLIN%energy_min=trH
          diisLIN%switchSD=.false.
          itBest=it
          icountSDSatur=icountSDSatur+1
          icountDIISFailureCons=0

          ! If we are using SD (i.e. diisLIN%idsx==0) and the trace has been decreasing
          ! for at least 10 iterations, switch to DIIS. However the history length is decreased.
          if(icountSDSatur>=10 .and. diisLIN%idsx==0 .and. allowDIIS) then
              icountSwitch=icountSwitch+1
              idsx=max(lin%DIISHistMin,lin%DIISHistMax-icountSwitch)
              if(idsx>0) then
                  if(iproc==0) write(*,'(x,a,i0)') 'switch to DIIS with new history length ', idsx
                  call initializeDIISParameters(idsx)
                  icountDIISFailureTot=0
                  icountDIISFailureCons=0
              end if
          end if
      else
          ! The trace is growing.
          ! Count how many times this occurs and (if we are using DIIS) switch to SD after 3 
          ! total failures or after 2 consecutive failures.
          icountDIISFailureCons=icountDIISFailureCons+1
          icountDIISFailureTot=icountDIISFailureTot+1
          icountSDSatur=0
          if((icountDIISFailureCons>=2 .or. icountDIISFailureTot>=3) .and. diisLIN%idsx>0) then
              ! Switch back to SD. The initial step size is 1.d0.
              alpha=lin%alphaSD
              if(iproc==0) then
                  if(icountDIISFailureCons>=2) write(*,'(x,a,i0,a,es10.3)') 'DIIS failed ', &
                      icountDIISFailureCons, ' times consecutievly. Switch to SD with stepsize', alpha(1)
                  if(icountDIISFailureTot>=3) write(*,'(x,a,i0,a,es10.3)') 'DIIS failed ', &
                      icountDIISFailureTot, ' times in total. Switch to SD with stepsize', alpha(1)
              end if
              ! Try to get back the orbitals of the best iteration. This is possible if
              ! these orbitals are still present in the DIIS history.
              if(it-itBest<diisLIN%idsx) then
                 if(iproc==0) then
                     if(iproc==0) write(*,'(x,a,i0,a)')  'Recover the orbitals from iteration ', &
                         itBest, ' which are the best so far.'
                 end if
                 ii=modulo(diisLIN%mids-(it-itBest),diisLIN%mids)
                 nvctrp=lin%comms%nvctr_par(iproc,1) ! 1 for k-point
                 call dcopy(lin%orbs%norb*nvctrp, diisLIN%psidst(ii*nvctrp*lin%orbs%norb+1), 1, phi(1), 1)
              end if
              call deallocate_diis_objects(diisLIN, subname)
              diisLIN%idsx=0
              diisLIN%switchSD=.true.
          end if
      end if

    end subroutine DIISorSD


    subroutine improveOrbitals()
    !
    ! Purpose:
    ! ========
    !   This subroutine improves the basis functions by following the gradient 
    ! For DIIS 
    if (diisLIN%idsx > 0) then
       diisLIN%mids=mod(diisLIN%ids,diisLIN%idsx)+1
       diisLIN%ids=diisLIN%ids+1
    end if

    ! Follow the gradient using steepest descent.
    ! The same, but transposed
    call transpose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, hphi, work=phiWork)
    
    ! steepest descent
    if(diisLIN%idsx==0) then
        istart=1
        nvctrp=lin%comms%nvctr_par(iproc,1) ! 1 for k-point
        do iorb=1,lin%orbs%norb
            call daxpy(nvctrp*orbs%nspinor, -alpha(iorb), hphi(istart), 1, phi(istart), 1)
            istart=istart+nvctrp*orbs%nspinor
        end do
    else
        ! DIIS
        quiet=.true. ! less output
        istart=1
        nvctrp=lin%comms%nvctr_par(iproc,1) ! 1 for k-point
        do iorb=1,lin%orbs%norb
            call dscal(nvctrp, alphaDIIS(iorb), hphi(istart), 1)
            istart=istart+nvctrp*orbs%nspinor
        end do
        call psimix(iproc, nproc, lin%orbs, lin%comms, diisLIN, hphi, phi, quiet)
    end if
    end subroutine improveOrbitals



    subroutine allocateLocalArrays()
    !
    ! Purpose:
    ! ========
    !   This subroutine allocates all local arrays.
    !
      allocate(hphi(lin%orbs%npsidim), stat=istat)
      call memocc(istat, hphi, 'hphi', subname)

      allocate(hphiold(lin%orbs%npsidim), stat=istat)
      call memocc(istat, hphiold, 'hphiold', subname)

      allocate(alpha(lin%orbs%norb), stat=istat)
      call memocc(istat, alpha, 'alpha', subname)

      allocate(alphaDIIS(lin%orbs%norb), stat=istat)
      call memocc(istat, alphaDIIS, 'alphaDIIS', subname)

      allocate(fnrmArr(lin%orbs%norb,2), stat=istat)
      call memocc(istat, fnrmArr, 'fnrmArr', subname)

      allocate(fnrmOldArr(lin%orbs%norb), stat=istat)
      call memocc(istat, fnrmOldArr, 'fnrmOldArr', subname)

      allocate(fnrmOvrlpArr(lin%orbs%norb,2), stat=istat)
      call memocc(istat, fnrmOvrlpArr, 'fnrmOvrlpArr', subname)

      allocate(phiWork(size(phi)), stat=istat)
      call memocc(istat, phiWork, 'phiWork', subname)
      
      allocate(lphi(lin%Lorbs%npsidim), stat=istat)
      call memocc(istat, lphi, 'lphi', subname)

      allocate(lhphi(lin%Lorbs%npsidim), stat=istat)
      call memocc(istat, lhphi, 'lhphi', subname)
    

    end subroutine allocateLocalArrays


    subroutine deallocateLocalArrays()
    !
    ! Purpose:
    ! ========
    !   This subroutine deallocates all local arrays.
    !

      iall=-product(shape(hphiold))*kind(hphiold)
      deallocate(hphiold, stat=istat)
      call memocc(istat, iall, 'hphiold', subname)

      iall=-product(shape(hphi))*kind(hphi)
      deallocate(hphi, stat=istat)
      call memocc(istat, iall, 'hphi', subname)
      
      iall=-product(shape(alpha))*kind(alpha)
      deallocate(alpha, stat=istat)
      call memocc(istat, iall, 'alpha', subname)

      iall=-product(shape(alphaDIIS))*kind(alphaDIIS)
      deallocate(alphaDIIS, stat=istat)
      call memocc(istat, iall, 'alphaDIIS', subname)

      iall=-product(shape(fnrmArr))*kind(fnrmArr)
      deallocate(fnrmArr, stat=istat)
      call memocc(istat, iall, 'fnrmArr', subname)

      iall=-product(shape(fnrmOldArr))*kind(fnrmOldArr)
      deallocate(fnrmOldArr, stat=istat)
      call memocc(istat, iall, 'fnrmOldArr', subname)

      iall=-product(shape(fnrmOvrlpArr))*kind(fnrmOvrlpArr)
      deallocate(fnrmOvrlpArr, stat=istat)
      call memocc(istat, iall, 'fnrmOvrlpArr', subname)

      iall=-product(shape(phiWork))*kind(phiWork)
      deallocate(phiWork, stat=istat)
      call memocc(istat, iall, 'phiWork', subname)

      iall=-product(shape(lphi))*kind(lphi)
      deallocate(lphi, stat=istat)
      call memocc(istat, iall, 'lphi', subname)

      iall=-product(shape(lhphi))*kind(lhphi)
      deallocate(lhphi, stat=istat)
      call memocc(istat, iall, 'lhphi', subname)
      
      ! if diisLIN%idsx==0, these arrays have already been deallocated
      if(diisLIN%idsx>0 .and. lin%DIISHistMax>0) call deallocate_diis_objects(diisLIN,subname)

    end subroutine deallocateLocalArrays


end subroutine getLocalizedBasis






subroutine transformHam(iproc, nproc, orbs, comms, phi, hphi, HamSmall)
!
! Purpose:
! =======
!   Builds the Hamiltonian in the basis of the localized basis functions phi. To do so, it gets all basis
!   functions |phi_i> and H|phi_i> and then calculates H_{ij}=<phi_i|H|phi_j>. The basis functions phi are
!   provided in the transposed form.
!
! Calling arguments:
! ==================
!   Input arguments:
!   ----------------
!     iproc      process ID
!     nproc      total number of processes
!     orbs       type describing the basis functions psi
!     comms      type containing the communication parameters for the physical orbitals phi
!     phi        basis functions 
!     hphi       the Hamiltonian applied to the basis functions 
!   Output arguments:
!   -----------------
!     HamSmall   Hamiltonian in small basis
!
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data), intent(in) :: orbs
type(communications_arrays), intent(in) :: comms
real(8),dimension(sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor,orbs%norb), intent(in) :: phi, hphi
real(8),dimension(orbs%norb,orbs%norb),intent(out):: HamSmall

! Local variables
integer:: istat, ierr, nvctrp, iall
real(8),dimension(:,:),allocatable:: HamTemp
character(len=*),parameter:: subname='transformHam'



  ! Allocate a temporary array if there are several MPI processes
  if(nproc>1) then
      allocate(HamTemp(orbs%norb,orbs%norb), stat=istat)
      call memocc(istat, HamTemp, 'HamTemp', subname)
  end if
  
  ! nvctrp is the amount of each phi hold by the current process
  nvctrp=sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor
  
  ! Build the Hamiltonian. In the parallel case, each process writes its Hamiltonian in HamTemp
  ! and a mpi_allreduce sums up the contribution from all processes.
  if(nproc==1) then
      call dgemm('t', 'n', orbs%norb, orbs%norb, nvctrp, 1.d0, phi(1,1), nvctrp, &
                 hphi(1,1), nvctrp, 0.d0, HamSmall(1,1), orbs%norb)
  else
      call dgemm('t', 'n', orbs%norb, orbs%norb, nvctrp, 1.d0, phi(1,1), nvctrp, &
                 hphi(1,1), nvctrp, 0.d0, HamTemp(1,1), orbs%norb)
  end if
  if(nproc>1) then
      call mpi_allreduce(HamTemp(1,1), HamSmall(1,1), orbs%norb**2, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
  end if
  
  if(nproc>1) then
     iall=-product(shape(HamTemp))*kind(HamTemp)
     deallocate(HamTemp,stat=istat)
     call memocc(istat, iall, 'HamTemp', subname)
  end if

end subroutine transformHam




subroutine diagonalizeHamiltonian(iproc, nproc, orbs, HamSmall, eval)
!
! Purpose:
! ========
!   Diagonalizes the Hamiltonian HamSmall and makes sure that all MPI processes give
!   the same result. This is done by requiring that the first entry of each vector
!   is positive.
!
! Calling arguments:
! ==================
!   Input arguments:
!   ----------------
!     iproc     process ID
!     nproc     number of MPI processes
!     orbs      type describing the physical orbitals psi
!   Input / Putput arguments
!     HamSmall  on input: the Hamiltonian
!               on exit: the eigenvectors
!   Output arguments
!     eval      the associated eigenvalues 
!
use module_base
use module_types
implicit none

! Calling arguments
integer:: iproc, nproc
type(orbitals_data), intent(inout) :: orbs
real(8),dimension(orbs%norb, orbs%norb):: HamSmall
real(8),dimension(orbs%norb):: eval

! Local variables
integer:: lwork, info, istat, iall, i, iorb, jorb
real(8),dimension(:),allocatable:: work
character(len=*),parameter:: subname='diagonalizeHamiltonian'

  ! Get the optimal work array size
  lwork=-1 
  allocate(work(1), stat=istat)
  call memocc(istat, work, 'work', subname)
  call dsyev('v', 'l', orbs%norb, HamSmall(1,1), orbs%norb, eval(1), work(1), lwork, info) 
  lwork=work(1) 

  ! Deallocate the work array ane reallocate it with the optimal size
  iall=-product(shape(work))*kind(work)
  deallocate(work, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating work' 
  call memocc(istat, iall, 'work', subname)
  allocate(work(lwork), stat=istat) ; if(istat/=0) stop 'ERROR in allocating work' 
  call memocc(istat, work, 'work', subname)

  ! Diagonalize the Hamiltonian
  call dsyev('v', 'l', orbs%norb, HamSmall(1,1), orbs%norb, eval(1), work(1), lwork, info) 

  ! Deallocate the work array.
  iall=-product(shape(work))*kind(work)
  deallocate(work, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating work' 
  call memocc(istat, iall, 'work', subname)
  
  ! Make sure that the eigenvectors are the same for all MPI processes. To do so, require that 
  ! the first entry of each vector is positive.
  do iorb=1,orbs%norb
      if(HamSmall(1,iorb)<0.d0) then
          do jorb=1,orbs%norb
              HamSmall(jorb,iorb)=-HamSmall(jorb,iorb)
          end do
      end if
  end do


end subroutine diagonalizeHamiltonian





subroutine buildWavefunction(iproc, nproc, orbs, orbsLIN, comms, commsLIN, phi, psi, HamSmall)
!
! Purpose:
! =======
!   Builds the physical orbitals psi as a linear combination of the basis functions phi. The coefficients
!   for this linear combination are obtained by diagonalizing the Hamiltonian matrix HamSmall.
!
! Calling arguments:
! ==================
!   Input arguments:
!   ----------------
!     iproc      process ID
!     nproc      total number of processes
!     orbs       type describing the physical orbitals psi
!     orbsLIN    type describing the basis functions phi
!     comms      type containing the communication parameters for the physical orbitals psi
!     commsLIN   type containing the communication parameters for the basis functions phi
!     phi        the basis functions 
!     HamSmall   the  Hamiltonian matrix
!   Output arguments:
!   -----------------
!     psi        the physical orbitals 
!

use module_base
use module_types
implicit none

! Calling arguments
integer:: iproc, nproc
type(orbitals_data), intent(in) :: orbs
type(orbitals_data), intent(in) :: orbsLIN
type(communications_arrays), intent(in) :: comms
type(communications_arrays), intent(in) :: commsLIN
real(8),dimension(sum(commsLIN%nvctr_par(iproc,1:orbsLIN%nkptsp))*orbsLIN%nspinor,orbsLIN%norb) :: phi
real(8),dimension(sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor,orbs%norb) :: psi
real(8),dimension(orbsLIN%norb,orbsLIN%norb):: HamSmall

! Local variables
integer:: nvctrp


  nvctrp=sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor
  call dgemm('n', 'n', nvctrp, orbs%norb, orbsLIN%norb, 1.d0, phi(1,1), nvctrp, HamSmall(1,1), &
             orbsLIN%norb, 0.d0, psi(1,1), nvctrp)
  

end subroutine buildWavefunction





subroutine buildWavefunctionModified(iproc, nproc, orbs, orbsLIN, comms, commsLIN, phi, psi, coeff)
!
! Purpose:
! =======
!   Builds the physical orbitals psi as a linear combination of the basis functions phi. The coefficients
!   for this linear combination are obtained by diagonalizing the Hamiltonian matrix HamSmall.
!
! Calling arguments:
! ==================
!   Input arguments:
!   ----------------
!     iproc      process ID
!     nproc      total number of processes
!     orbs       type describing the physical orbitals psi
!     orbsLIN    type describing the basis functions phi
!     comms      type containing the communication parameters for the physical orbitals psi
!     commsLIN   type containing the communication parameters for the basis functions phi
!     phi        the basis functions 
!     coeff      the coefficients for the linear combination
!   Output arguments:
!   -----------------
!     psi        the physical orbitals 
!

use module_base
use module_types
implicit none

! Calling arguments
integer:: iproc, nproc
type(orbitals_data), intent(in) :: orbs
type(orbitals_data), intent(in) :: orbsLIN
type(communications_arrays), intent(in) :: comms
type(communications_arrays), intent(in) :: commsLIN
real(8),dimension(sum(commsLIN%nvctr_par(iproc,1:orbsLIN%nkptsp))*orbsLIN%nspinor,orbsLIN%norb) :: phi
real(8),dimension(sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor,orbs%norb) :: psi
real(8),dimension(orbsLIN%norb,orbs%norb):: coeff

! Local variables
integer:: nvctrp


  nvctrp=sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor
  call dgemm('n', 'n', nvctrp, orbs%norb, orbsLIN%norb, 1.d0, phi(1,1), nvctrp, coeff(1,1), &
             orbsLIN%norb, 0.d0, psi(1,1), nvctrp)
  

end subroutine buildWavefunctionModified


subroutine getMatrixElements(iproc, nproc, Glr, orbs, comms, phi, hphi, matrixElements)
!
! Purpose:
! ========
!
! Calling arguments:
! ==================
!
use module_base
use module_types
use module_interfaces, exceptThisOne => getMatrixElements
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(locreg_descriptors),intent(in):: Glr
type(orbitals_data),intent(in):: orbs
type(communications_arrays),intent(in):: comms
real(8),dimension(orbs%npsidim),intent(inout):: phi, hphi
real(8),dimension(orbs%norb,orbs%norb,2),intent(out):: matrixElements

! Local variables
integer:: istart, jstart, nvctrp, iorb, jorb, istat, iall, ierr
real(8):: ddot
real(8),dimension(:),pointer:: phiWork
character(len=*),parameter:: subname='getMatrixELements'


  allocate(phiWork(orbs%npsidim), stat=istat)
  call memocc(istat, phiWork, 'phiWork', subname)


  call transpose_v(iproc, nproc, orbs, Glr%wfd, comms, phi, work=phiWork)
  call transpose_v(iproc, nproc, orbs, Glr%wfd, comms, hphi, work=phiWork)

  matrixElements=0.d0

  ! Calculate <phi_i|H_j|phi_j>
  nvctrp=sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor
  jstart=1
  do jorb=1,orbs%norb
      istart=1
      do iorb=1,orbs%norb
          matrixElements(iorb,jorb,2)=ddot(nvctrp, phi(istart), 1, hphi(jstart), 1)
          istart=istart+nvctrp
      end do
      jstart=jstart+nvctrp
  end do
  call mpi_allreduce(matrixElements(1,1,2), matrixElements(1,1,1), orbs%norb**2, &
      mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
!!if(iproc==0) then
!!    write(*,*) 'matrix Elements'
!!    do iorb=1,lin%orbs%norb
!!        write(*,'(80es9.2)') (matrixElements(iorb,jorb,1), jorb=1,lin%orbs%norb)
!!    end do
!!end if


  call untranspose_v(iproc, nproc, orbs, Glr%wfd, comms, phi, work=phiWork)
  call untranspose_v(iproc, nproc, orbs, Glr%wfd, comms, hphi, work=phiWork)

  iall=-product(shape(phiWork))*kind(phiWork)
  deallocate(phiWork)
  call memocc(istat, iall, 'phiWork', subname)


  !!! Calculate the modified band structure energy
  !!tt=0.d0
  !!do iorb=1,orbs%norb
  !!    do jorb=1,orbsLIN%norb
  !!        do korb=1,orbsLIN%norb
  !!            tt=tt+HamSmall(korb,iorb)*HamSmall(jorb,iorb)*matrixElements(korb,jorb,1)
  !!        end do
  !!    end do
  !!end do
  !!if(present(ebs_mod)) then
  !!    if(nspin==1) ebs_mod=2.d0*tt ! 2 for closed shell
  !!end if



end subroutine getMatrixElements





subroutine modifiedBSEnergy(nspin, orbs, lin, HamSmall, matrixElements, ebsMod)
!
! Purpose:
! ========
!
! Calling arguments:
! ==================
!
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: nspin
type(orbitals_data),intent(in) :: orbs
type(linearParameters),intent(in):: lin
real(8),dimension(lin%orbs%norb,lin%orbs%norb),intent(in):: HamSmall, matrixElements
real(8),intent(out):: ebsMod

! Local variables
integer:: iorb, jorb, korb
real(8):: tt

  ! Calculate the modified band structure energy
  tt=0.d0
  do iorb=1,orbs%norb
      do jorb=1,lin%orbs%norb
          do korb=1,lin%orbs%norb
              tt=tt+HamSmall(korb,iorb)*HamSmall(jorb,iorb)*matrixElements(korb,jorb)
          end do
      end do
  end do
  if(nspin==1) then
      ebsMod=2.d0*tt ! 2 for closed shell
  else
      ebsMod=tt
  end if



end subroutine modifiedBSEnergy





subroutine modifiedBSEnergyModified(nspin, orbs, lin, coeff, matrixElements, ebsMod)
!
! Purpose:
! ========
!
! Calling arguments:
! ==================
!
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: nspin
type(orbitals_data),intent(in) :: orbs
type(linearParameters),intent(in):: lin
real(8),dimension(lin%orbs%norb,orbs%norb),intent(in):: coeff
real(8),dimension(lin%orbs%norb,lin%orbs%norb),intent(in):: matrixElements
real(8),intent(out):: ebsMod

! Local variables
integer:: iorb, jorb, korb
real(8):: tt

  ! Calculate the modified band structure energy
  tt=0.d0
  do iorb=1,orbs%norb
      do jorb=1,lin%orbs%norb
          do korb=1,lin%orbs%norb
              tt=tt+coeff(korb,iorb)*coeff(jorb,iorb)*matrixElements(korb,jorb)
          end do
      end do
  end do
  if(nspin==1) then
      ebsMod=2.d0*tt ! 2 for closed shell
  else
      ebsMod=tt
  end if



end subroutine modifiedBSEnergyModified












subroutine optimizeCoefficients(iproc, orbs, lin, nspin, matrixElements, coeff, infoCoeff)
!
! Purpose:
! ========
!   Determines the optimal coefficients which minimize the modified band structure energy, i.e.
!   E = sum_{i}sum_{k,l}c_{ik}c_{il}<phi_k|H_l|phi_l>.
!   This is done by a steepest descen minimization using the gradient of the above expression with
!   respect to the coefficients c_{ik}.
!
! Calling arguments:
! ==================
!   Input arguments:
!   ----------------
!     iproc            process ID
!     orbs             type describing the physical orbitals psi
!     lin              type containing parameters for the linear version
!     nspin            nspin==1 -> closed shell, npsin==2 -> open shell
!     matrixElements   contains the matrix elements <phi_k|H_l|phi_l>
!   Output arguments:
!   -----------------
!     coeff            the optimized coefficients 
!     infoCoeff        if infoCoeff=0, the optimization converged
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nspin
type(orbitals_data),intent(in):: orbs
type(linearParameters),intent(in):: lin
real(8),dimension(lin%lb%orbs%norb,lin%lb%orbs%norb),intent(in):: matrixElements
real(8),dimension(lin%lb%orbs%norb,orbs%norb),intent(inout):: coeff
integer,intent(out):: infoCoeff

! Local variables
integer:: it, iorb, jorb, k, l, istat, iall, korb, ierr
real(8):: tt, fnrm, ddot, dnrm2, meanAlpha, cosangle, ebsMod, ebsModOld
real(8),dimension(:,:),allocatable:: grad, gradOld, lagMat
real(8),dimension(:),allocatable:: alpha
character(len=*),parameter:: subname='optimizeCoefficients'
logical:: converged


! Allocate all local arrays.
allocate(grad(lin%lb%orbs%norb,orbs%norb), stat=istat)
call memocc(istat, grad, 'grad', subname)
allocate(gradOld(lin%lb%orbs%norb,orbs%norb), stat=istat)
call memocc(istat, gradOld, 'gradOld', subname)
allocate(lagMat(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, lagMat, 'lagMat', subname)
allocate(alpha(orbs%norb), stat=istat)
call memocc(istat, alpha, 'alpha', subname)

! trace of matrixElements
if(iproc==0) then
    tt=0.d0
    do iorb=1,lin%lb%orbs%norb
        do jorb=1,lin%lb%orbs%norb
            if(iorb==jorb) tt=tt+matrixElements(iorb,jorb)
            !write(777,*) iorb,jorb,matrixElements(jorb,iorb)
        end do
    end do
    !write(*,*) 'trace',tt
end if

! Do everything only on the root process and then broadcast to all processes.
! Maybe this part can be parallelized later, but at the moment it is not necessary since
! it is fast enough.
processIf: if(iproc==0) then
    

    !!! Orthonormalize the coefficient vectors (Gram-Schmidt).
    !!do iorb=1,orbs%norb
    !!    do jorb=1,iorb-1
    !!        tt=ddot(lin%lb%orbs%norb, coeff(1,iorb), 1, coeff(1,jorb), 1)
    !!        call daxpy(lin%lb%orbs%norb, -tt, coeff(1,jorb), 1, coeff(1,iorb), 1)
    !!    end do
    !!    tt=dnrm2(lin%lb%orbs%norb, coeff(1,iorb), 1)
    !!    call dscal(lin%lb%orbs%norb, 1/tt, coeff(1,iorb), 1)
    !!end do
    
    ! Initial step size for the optimization
    alpha=5.d-3

    ! Flag which checks convergence.
    converged=.false.

    if(iproc==0) write(*,'(x,a)') '============================== optmizing coefficients =============================='

    ! The optimization loop.
    iterLoop: do it=1,lin%nItCoeff

        if (iproc==0) then
            write( *,'(1x,a,i0)') repeat('-',77 - int(log(real(it))/log(10.))) // ' iter=', it
        endif


        ! Orthonormalize the coefficient vectors (Gram-Schmidt).
        do iorb=1,orbs%norb
            do jorb=1,iorb-1
                tt=ddot(lin%lb%orbs%norb, coeff(1,iorb), 1, coeff(1,jorb), 1)
                call daxpy(lin%lb%orbs%norb, -tt, coeff(1,jorb), 1, coeff(1,iorb), 1)
            end do
            tt=dnrm2(lin%lb%orbs%norb, coeff(1,iorb), 1)
            call dscal(lin%lb%orbs%norb, 1/tt, coeff(1,iorb), 1)
        end do


        ! Calculate the gradient grad. At the same time we determine whether the step size shall be increased
        ! or decreased (depending on gradient feedback).
        meanAlpha=0.d0
        grad=0.d0
        do iorb=1,orbs%norb
            do l=1,lin%lb%orbs%norb
                do k=1,lin%lb%orbs%norb
                    grad(l,iorb)=grad(l,iorb)+coeff(k,iorb)*(matrixElements(k,l)+matrixElements(l,k))
                end do
            end do
            if(it>1) then
                cosangle=ddot(lin%lb%orbs%norb, grad(1,iorb), 1, gradOld(1,iorb), 1)
                cosangle=cosangle/dnrm2(lin%lb%orbs%norb, grad(1,iorb), 1)
                cosangle=cosangle/dnrm2(lin%lb%orbs%norb, gradOld(1,iorb), 1)
                !write(*,*) 'cosangle, ebsMod, ebsModOld', cosangle, ebsMod, ebsModOld
                if(cosangle>.8d0 .and. ebsMod<ebsModOld+1.d-6*abs(ebsModOld)) then
                    alpha(iorb)=max(alpha(iorb)*1.05d0,1.d-6)
                else
                    alpha(iorb)=max(alpha(iorb)*.5d0,1.d-6)
                end if
            end if
            call dcopy(lin%lb%orbs%norb, grad(1,iorb), 1, gradOld(1,iorb), 1)
            meanAlpha=meanAlpha+alpha(iorb)
        end do
        meanAlpha=meanAlpha/orbs%norb
    
    
        ! Apply the orthoconstraint to the gradient. To do so first calculate the Lagrange
        ! multiplier matrix.
        lagMat=0.d0
        do iorb=1,orbs%norb
            do jorb=1,orbs%norb
                do k=1,lin%lb%orbs%norb
                    lagMat(iorb,jorb)=lagMat(iorb,jorb)+coeff(k,iorb)*grad(k,jorb)
                end do
            end do
        end do

        ! Now apply the orthoconstraint.
        do iorb=1,orbs%norb
            do k=1,lin%lb%orbs%norb
                do jorb=1,orbs%norb
                    grad(k,iorb)=grad(k,iorb)-.5d0*(lagMat(iorb,jorb)*coeff(k,jorb)+lagMat(jorb,iorb)*coeff(k,jorb))
                end do
            end do
        end do
    
        
        ! Calculate the modified band structure energy and the gradient norm.
        if(it>1) then
            ebsModOld=ebsMod
        else
            ebsModOld=1.d10
        end if
        ebsMod=0.d0
        fnrm=0.d0
        do iorb=1,orbs%norb
            fnrm=fnrm+dnrm2(lin%lb%orbs%norb, grad(1,iorb), 1)
            do jorb=1,lin%lb%orbs%norb
                do korb=1,lin%lb%orbs%norb
                    ebsMod=ebsMod+coeff(korb,iorb)*coeff(jorb,iorb)*matrixElements(korb,jorb)
                end do
            end do
        end do
    
        ! Multiply the energy with a factor of 2 if we have a closed-shell system.
        if(nspin==1) then
            ebsMod=2.d0*ebsMod
        end if

        !if(iproc==0) write(*,'(x,a,4x,i0,es12.4,3x,es10.3, es19.9)') 'iter, fnrm, meanAlpha, Energy', &
        if(iproc==0) write(*,'(x,a,es11.2,es22.13,es10.2)') 'fnrm, band structure energy, mean alpha', &
            fnrm, ebsMod, meanAlpha
        
        ! Check for convergence.
        if(fnrm<lin%convCritCoeff) then
            if(iproc==0) write(*,'(x,a,i0,a)') 'converged in ', it, ' iterations.'
            if(iproc==0) write(*,'(3x,a,2es14.5)') 'Final values for fnrm, Energy:', fnrm, ebsMod
            converged=.true.
            infoCoeff=it
            exit
        end if
  
        if(it==lin%nItCoeff) then
            if(iproc==0) write(*,'(x,a,i0,a)') 'WARNING: not converged within ', it, &
                ' iterations! Exiting loop due to limitations of iterations.'
            if(iproc==0) write(*,'(x,a,2es15.7,f12.7)') 'Final values for fnrm, Energy: ', fnrm, ebsMod
            infoCoeff=-1
            exit
        end if

        ! Improve the coefficients (by steepet descent).
        do iorb=1,orbs%norb
            do l=1,lin%lb%orbs%norb
                coeff(l,iorb)=coeff(l,iorb)-alpha(iorb)*grad(l,iorb)
            end do
        end do
    

    end do iterLoop

    !!if(.not.converged) then
    !!    if(iproc==0) write(*,'(x,a,i0,a)') 'WARNING: not converged within ', it, &
    !!        ' iterations! Exiting loop due to limitations of iterations.'
    !!    if(iproc==0) write(*,'(x,a,2es15.7,f12.7)') 'Final values for fnrm, Energy: ', fnrm, ebsMod
    !!    infoCoeff=-1
    !!    ! Orthonormalize the coefficient vectors (Gram-Schmidt).
    !!    do iorb=1,orbs%norb
    !!        do jorb=1,iorb-1
    !!            tt=ddot(lin%lb%orbs%norb, coeff(1,iorb), 1, coeff(1,jorb), 1)
    !!            call daxpy(lin%lb%orbs%norb, -tt, coeff(1,jorb), 1, coeff(1,iorb), 1)
    !!        end do
    !!        tt=dnrm2(lin%lb%orbs%norb, coeff(1,iorb), 1)
    !!        call dscal(lin%lb%orbs%norb, 1/tt, coeff(1,iorb), 1)
    !!    end do
    !!end if

    if(iproc==0) write(*,'(x,a)') '===================================================================================='

end if processIf


! Now broadcast the result to all processes
call mpi_bcast(coeff(1,1), lin%lb%orbs%norb*orbs%norb, mpi_double_precision, 0, mpi_comm_world, ierr)
call mpi_bcast(infoCoeff, 1, mpi_integer, 0, mpi_comm_world, ierr)


! Deallocate all local arrays.
iall=-product(shape(grad))*kind(grad)
deallocate(grad, stat=istat)
call memocc(istat, iall, 'grad', subname)

iall=-product(shape(gradOld))*kind(gradOld)
deallocate(gradOld, stat=istat)
call memocc(istat, iall, 'gradOld', subname)

iall=-product(shape(lagMat))*kind(lagMat)
deallocate(lagMat, stat=istat)
call memocc(istat, iall, 'lagMat', subname)

iall=-product(shape(alpha))*kind(alpha)
deallocate(alpha, stat=istat)
call memocc(istat, iall, 'alpha', subname)

end subroutine optimizeCoefficients








subroutine getDerivativeBasisFunctions(iproc, nproc, hgrid, Glr, lin, nphi, phi, phid)
use module_base
use module_types
use module_interfaces, exceptThisOne => getDerivativeBasisFunctions
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, nphi
real(8),intent(in):: hgrid
type(locreg_descriptors),intent(in):: Glr
type(linearParameters),intent(in):: lin
real(8),dimension(nphi),intent(in):: phi
real(8),dimension(lin%lb%orbs%npsidim),intent(out):: phid

! Local variables
integer:: ist1_c, ist1_f, ist2_c, ist2_f, nf, istat, iall, iorb, jproc, ierr
integer:: ist0_c, istx_c, isty_c, istz_c, ist0_f, istx_f, isty_f, istz_f, istLoc, istRoot
real(8),dimension(0:3),parameter:: scal=1.d0
real(8),dimension(:),allocatable:: w_f1, w_f2, w_f3, phiLoc, phiRoot
real(8),dimension(:,:,:),allocatable:: w_c, phix_c, phiy_c, phiz_c
real(8),dimension(:,:,:,:),allocatable:: w_f, phix_f, phiy_f, phiz_f
character(len=*),parameter:: subname='getDerivativeBasisFunctions'
integer,dimension(:),allocatable:: recvcounts, sendcounts, displs



  ! THIS IS COPIED FROM allocate_work_arrays. Works only for free boundary.
  nf=(Glr%d%nfu1-Glr%d%nfl1+1)*(Glr%d%nfu2-Glr%d%nfl2+1)*(Glr%d%nfu3-Glr%d%nfl3+1)
  ! Allocate work arrays
  allocate(w_c(0:Glr%d%n1,0:Glr%d%n2,0:Glr%d%n3+ndebug), stat=istat)
  call memocc(istat, w_c, 'w_c', subname)
  allocate(w_f(7,Glr%d%nfl1:Glr%d%nfu1,Glr%d%nfl2:Glr%d%nfu2,Glr%d%nfl3:Glr%d%nfu3+ndebug), stat=istat)
  call memocc(istat, w_f, 'w_f', subname)

  allocate(w_f1(nf+ndebug), stat=istat)
  call memocc(istat, w_f1, 'w_f1', subname)
  allocate(w_f2(nf+ndebug), stat=istat)
  call memocc(istat, w_f2, 'w_f2', subname)
  allocate(w_f3(nf+ndebug), stat=istat)
  call memocc(istat, w_f3, 'w_f3', subname)


  allocate(phix_f(7,Glr%d%nfl1:Glr%d%nfu1,Glr%d%nfl2:Glr%d%nfu2,Glr%d%nfl3:Glr%d%nfu3), stat=istat)
  call memocc(istat, phix_f, 'phix_f', subname)
  allocate(phix_c(0:Glr%d%n1,0:Glr%d%n2,0:Glr%d%n3), stat=istat)
  call memocc(istat, phix_c, 'phix_c', subname)
  allocate(phiy_f(7,Glr%d%nfl1:Glr%d%nfu1,Glr%d%nfl2:Glr%d%nfu2,Glr%d%nfl3:Glr%d%nfu3), stat=istat)
  call memocc(istat, phiy_f, 'phiy_f', subname)
  allocate(phiy_c(0:Glr%d%n1,0:Glr%d%n2,0:Glr%d%n3), stat=istat)
  call memocc(istat, phiy_c, 'phiy_c', subname)
  allocate(phiz_f(7,Glr%d%nfl1:Glr%d%nfu1,Glr%d%nfl2:Glr%d%nfu2,Glr%d%nfl3:Glr%d%nfu3), stat=istat)
  call memocc(istat, phiz_f, 'phiz_f', subname)
  allocate(phiz_c(0:Glr%d%n1,0:Glr%d%n2,0:Glr%d%n3), stat=istat)
  call memocc(istat, phiz_c, 'phiz_c', subname)

  allocate(phiLoc(4*lin%orbs%norbp*(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)), stat=istat)
  call memocc(istat, phiLoc, 'phiLoc', subname)
 
  ! Initialize to zero
  w_c=0.d0
  w_f=0.d0
  w_f1=0.d0
  w_f2=0.d0
  w_f3=0.d0
  phix_c=0.d0
  phix_f=0.d0
  phiy_c=0.d0
  phiy_f=0.d0
  phiz_c=0.d0
  phiz_f=0.d0



!write(*,'(a,i12,i5,i12)') 'Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f, lind%orbs%norbp, lind%orbs%npsidim', Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f, lind%orbs%norbp, lind%orbs%npsidim
!write(*,*) 'size(phid)', size(phid)

  ist1_c=1
  ist1_f=1+Glr%wfd%nvctr_c
  ist2_c=1
  ist2_f=1+Glr%wfd%nvctr_c

  ! These are the starting indices in phiLoc:
  ! ist0: start index of the orginal phi
  ist0_c=1
  ist0_f=ist0_c+Glr%wfd%nvctr_c
  ! istx: start index of the derivative with respect to x
  istx_c=1+lin%orbs%norbp*(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)
  istx_f=istx_c+Glr%wfd%nvctr_c
  ! isty: start index of the derivative with respect to y
  isty_c=1+lin%orbs%norbp*2*(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)
  isty_f=isty_c+Glr%wfd%nvctr_c
  ! istz: start index of the derivative with respect to z
  istz_c=1+lin%orbs%norbp*3*(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)
  istz_f=istz_c+Glr%wfd%nvctr_c
  do iorb=1,lin%orbs%norbp
      ! Uncompress the wavefunction.
      call uncompress_forstandard(Glr%d%n1, Glr%d%n2, Glr%d%n3, Glr%d%nfl1, Glr%d%nfu1, & 
           Glr%d%nfl2, Glr%d%nfu2, Glr%d%nfl3, Glr%d%nfu3,  &
           Glr%wfd%nseg_c, Glr%wfd%nvctr_c, Glr%wfd%keyg, Glr%wfd%keyv,  &
           Glr%wfd%nseg_f, Glr%wfd%nvctr_f, Glr%wfd%keyg(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)), &
           Glr%wfd%keyv(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),  &
           scal, phi(ist1_c), phi(ist1_f), w_c, w_f, w_f1, w_f2, w_f3)

      call createDerivativeBasis(Glr%d%n1, Glr%d%n2, Glr%d%n3, &
           Glr%d%nfl1, Glr%d%nfu1, Glr%d%nfl2, Glr%d%nfu2, Glr%d%nfl3, Glr%d%nfu3,  &
           hgrid, Glr%bounds%kb%ibyz_c, Glr%bounds%kb%ibxz_c, Glr%bounds%kb%ibxy_c, &
           Glr%bounds%kb%ibyz_f, Glr%bounds%kb%ibxz_f, Glr%bounds%kb%ibxy_f, &
           w_c, w_f, w_f1, w_f2, w_f3, phix_c, phix_f, phiy_c, phiy_f, phiz_c, phiz_f)

      ! Copy phi to phiLoc
      call dcopy(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f, phi(ist1_c), 1, phiLoc(ist0_c), 1)
      ist0_c = ist0_c + Glr%wfd%nvctr_c + 7*Glr%wfd%nvctr_f
      ist0_f = ist0_f + Glr%wfd%nvctr_c + 7*Glr%wfd%nvctr_f

      ! Compress the x wavefunction.
      call compress_forstandard(Glr%d%n1, Glr%d%n2, Glr%d%n3, Glr%d%nfl1, Glr%d%nfu1, &
           Glr%d%nfl2, Glr%d%nfu2, Glr%d%nfl3, Glr%d%nfu3, &
           Glr%wfd%nseg_c, Glr%wfd%nvctr_c, Glr%wfd%keyg, Glr%wfd%keyv, &
           Glr%wfd%nseg_f, Glr%wfd%nvctr_f, Glr%wfd%keyg(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)), &
           Glr%wfd%keyv(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),  &
           scal, phix_c, phix_f, phiLoc(istx_c), phiLoc(istx_f))
      istx_c = istx_c + Glr%wfd%nvctr_c + 7*Glr%wfd%nvctr_f
      istx_f = istx_f + Glr%wfd%nvctr_c + 7*Glr%wfd%nvctr_f

      ! Compress the y wavefunction.
      call compress_forstandard(Glr%d%n1, Glr%d%n2, Glr%d%n3, Glr%d%nfl1, Glr%d%nfu1, &
           Glr%d%nfl2, Glr%d%nfu2, Glr%d%nfl3, Glr%d%nfu3, &
           Glr%wfd%nseg_c, Glr%wfd%nvctr_c, Glr%wfd%keyg, Glr%wfd%keyv, &
           Glr%wfd%nseg_f, Glr%wfd%nvctr_f, Glr%wfd%keyg(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)), &
           Glr%wfd%keyv(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),  &
           scal, phiy_c, phiy_f, phiLoc(isty_c), phiLoc(isty_f))
      isty_c = isty_c + Glr%wfd%nvctr_c + 7*Glr%wfd%nvctr_f
      isty_f = isty_f + Glr%wfd%nvctr_c + 7*Glr%wfd%nvctr_f

      ! Compress the z wavefunction.
      call compress_forstandard(Glr%d%n1, Glr%d%n2, Glr%d%n3, Glr%d%nfl1, Glr%d%nfu1, &
           Glr%d%nfl2, Glr%d%nfu2, Glr%d%nfl3, Glr%d%nfu3, &
           Glr%wfd%nseg_c, Glr%wfd%nvctr_c, Glr%wfd%keyg, Glr%wfd%keyv, &
           Glr%wfd%nseg_f, Glr%wfd%nvctr_f, Glr%wfd%keyg(1,Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)), &
           Glr%wfd%keyv(Glr%wfd%nseg_c+min(1,Glr%wfd%nseg_f)),  &
           scal, phiz_c, phiz_f, phiLoc(istz_c), phiLoc(istz_f))
      istz_c = istz_c + Glr%wfd%nvctr_c + 7*Glr%wfd%nvctr_f
      istz_f = istz_f + Glr%wfd%nvctr_c + 7*Glr%wfd%nvctr_f

      ist1_c = ist1_c + Glr%wfd%nvctr_c + 7*Glr%wfd%nvctr_f
      ist1_f = ist1_f + Glr%wfd%nvctr_c + 7*Glr%wfd%nvctr_f

  end do

  ! Now copy phiLoc to phid. This requires some communication, since the partition of phiLoc
  ! is not identical to the partition of phid.
  ! Collect on root and then redistribute.
  ! THIS MUST BE IMPROVED.
  allocate(recvcounts(0:nproc-1), stat=istat)
  call memocc(istat, recvcounts, 'recvcounts', subname)
  allocate(sendcounts(0:nproc-1), stat=istat)
  call memocc(istat, sendcounts, 'sendcounts', subname)
  allocate(displs(0:nproc-1), stat=istat)
  call memocc(istat, displs, 'displs', subname)
  !if(iproc==0) then
      allocate(phiRoot(lin%lb%orbs%norb*(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)), stat=istat)
      call memocc(istat, phiRoot, 'phiRoot', subname)
  !end if

  ! Gather the original phi
  istLoc=1
  istRoot=1
  displs(0)=1
  do jproc=0,nproc-1
      if(jproc>0) displs(jproc)=displs(jproc-1)+recvcounts(jproc-1)
      recvcounts(jproc)=lin%orbs%norb_par(jproc)*(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)
      !if(iproc==0) write(*,'(a,i4,2i6)') '0: jproc, recvcounts(jproc)/size, displs(jproc)/size', jproc, recvcounts(jproc)/(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f), displs(jproc)/(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)
  end do
  call mpi_gatherv(phiLoc(istLoc), lin%orbs%norbp*(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f), mpi_double_precision, &
       phiRoot(istRoot), recvcounts, displs, mpi_double_precision, 0, mpi_comm_world, ierr)

  ! Gather the derivatives with respect to x
  istLoc=lin%orbs%norbp*(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)+1
  istRoot=lin%orbs%norb*(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)+1
  displs(0)=1
  do jproc=0,nproc-1
      if(jproc>0) displs(jproc)=displs(jproc-1)+recvcounts(jproc-1)
      recvcounts(jproc)=lin%orbs%norb_par(jproc)*(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)
      !if(iproc==0) write(*,'(a,i4,2i6)') 'x: jproc, recvcounts(jproc)/size, displs(jproc)/size', jproc, recvcounts(jproc)/(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f), displs(jproc)/(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)
  end do
  call mpi_gatherv(phiLoc(istLoc), lin%orbs%norbp*(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f), mpi_double_precision, &
       phiRoot(istRoot), recvcounts, displs, mpi_double_precision, 0, mpi_comm_world, ierr)

  ! Gather the derivatives with respect to y
  istLoc=2*lin%orbs%norbp*(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)+1
  istRoot=2*lin%orbs%norb*(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)+1
  displs(0)=1
  do jproc=0,nproc-1
      if(jproc>0) displs(jproc)=displs(jproc-1)+recvcounts(jproc-1)
      recvcounts(jproc)=lin%orbs%norb_par(jproc)*(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)
      !if(iproc==0) write(*,'(a,i4,2i6)') 'y: jproc, recvcounts(jproc)/size, displs(jproc)/size', jproc, recvcounts(jproc)/(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f), displs(jproc)/(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)
  end do
  call mpi_gatherv(phiLoc(istLoc), lin%orbs%norbp*(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f), mpi_double_precision, &
       phiRoot(istRoot), recvcounts, displs, mpi_double_precision, 0, mpi_comm_world, ierr)

  ! Gather the derivatives with respect to z
  istLoc=3*lin%orbs%norbp*(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)+1
  istRoot=3*lin%orbs%norb*(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)+1
  displs(0)=1
  do jproc=0,nproc-1
      if(jproc>0) displs(jproc)=displs(jproc-1)+recvcounts(jproc-1)
      recvcounts(jproc)=lin%orbs%norb_par(jproc)*(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)
      !if(iproc==0) write(*,'(a,i4,2i6)') 'z: jproc, recvcounts(jproc)/size, displs(jproc)/size', jproc, recvcounts(jproc)/(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f), displs(jproc)/(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)
  end do
  call mpi_gatherv(phiLoc(istLoc), lin%orbs%norbp*(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f), mpi_double_precision, &
       phiRoot(istRoot), recvcounts, displs, mpi_double_precision, 0, mpi_comm_world, ierr)


  displs(0)=1
  do jproc=0,nproc-1
      if(jproc>0) displs(jproc)=displs(jproc-1)+sendcounts(jproc-1)
      sendcounts(jproc)=lin%lb%orbs%norb_par(jproc)*(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)
       !if(iproc==0) write(*,'(a,i4,2i6)') 'jproc, sendcounts(jproc)/size, displs(jproc)/size', jproc, sendcounts(jproc)/(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f), displs(jproc)/(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)
  end do
  call mpi_scatterv(phiRoot(1), sendcounts, displs, mpi_double_precision, phid(1), &
       lin%lb%orbs%norbp*(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f), mpi_double_precision, 0, mpi_comm_world, ierr)

  ! Deallocate all local arrays
  iall=-product(shape(w_c))*kind(w_c)
  deallocate(w_c, stat=istat)
  call memocc(istat, iall, 'w_c', subname)

  iall=-product(shape(w_f))*kind(w_f)
  deallocate(w_f, stat=istat)
  call memocc(istat, iall, 'w_f', subname)

  iall=-product(shape(w_f1))*kind(w_f1)
  deallocate(w_f1, stat=istat)
  call memocc(istat, iall, 'w_f1', subname)

  iall=-product(shape(w_f2))*kind(w_f2)
  deallocate(w_f2, stat=istat)
  call memocc(istat, iall, 'w_f2', subname)

  iall=-product(shape(w_f3))*kind(w_f3)
  deallocate(w_f3, stat=istat)
  call memocc(istat, iall, 'w_f3', subname)

  iall=-product(shape(phix_f))*kind(phix_f)
  deallocate(phix_f, stat=istat)
  call memocc(istat, iall, 'phix_f', subname)

  iall=-product(shape(phix_c))*kind(phix_c)
  deallocate(phix_c, stat=istat)
  call memocc(istat, iall, 'phix_c', subname)

  iall=-product(shape(phiy_f))*kind(phiy_f)
  deallocate(phiy_f, stat=istat)
  call memocc(istat, iall, 'phiy_f', subname)

  iall=-product(shape(phiy_c))*kind(phiy_c)
  deallocate(phiy_c, stat=istat)
  call memocc(istat, iall, 'phiy_c', subname)

  iall=-product(shape(phiz_f))*kind(phiz_f)
  deallocate(phiz_f, stat=istat)
  call memocc(istat, iall, 'phiz_f', subname)

  iall=-product(shape(phiz_c))*kind(phiz_c)
  deallocate(phiz_c, stat=istat)
  call memocc(istat, iall, 'phiz_c', subname)

  iall=-product(shape(phiLoc))*kind(phiLoc)
  deallocate(phiLoc, stat=istat)
  call memocc(istat, iall, 'phiLoc', subname)

  iall=-product(shape(phiRoot))*kind(phiRoot)
  deallocate(phiRoot, stat=istat)
  call memocc(istat, iall, 'phiRoot', subname)

  iall=-product(shape(recvcounts))*kind(recvcounts)
  deallocate(recvcounts, stat=istat)
  call memocc(istat, iall, 'recvcounts', subname)

  iall=-product(shape(sendcounts))*kind(sendcounts)
  deallocate(sendcounts, stat=istat)
  call memocc(istat, iall, 'sendcounts', subname)

  iall=-product(shape(displs))*kind(displs)
  deallocate(displs, stat=istat)
  call memocc(istat, iall, 'displs', subname)


end subroutine getDerivativeBasisFunctions





subroutine orthonormalizeOnlyDerivatives(iproc, nproc, lin, phid)
!
! Purpose:
! ========
!   Firts orthogonalizes the derivative basis functions to the original basis functions
!   using Gram Schmidt. Then orthonormalizes the derivative basis functions with Loewdin.
!
use module_base
use module_defs
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(linearParameters),intent(in):: lin
real(8),dimension(lin%lb%orbs%npsidim),intent(inout):: phid

! Local varibles
integer:: iorb, jorb, ist, jst, nvctrp, ierr, istat, iall, lwork, info

real(8):: tt, ddot, dnrm2
real(8),dimension(:),allocatable:: work, eval
real(8),dimension(:,:),allocatable:: ovrlp, A, phidTemp
real(8),dimension(:,:,:),allocatable:: tempArr
character(len=*),parameter:: subname='orthonormalizeOnlyDerivatives'


allocate(ovrlp(lin%orbs%norb+1:4*lin%orbs%norb,1:lin%orbs%norb), stat=istat)
call memocc(istat, ovrlp, 'ovrlp', subname)

! Orthogonalize the derivative basis functions to the original ones.
nvctrp=lin%lb%comms%nvctr_par(iproc,1) ! 1 for k-point

allocate(A(nvctrp,lin%orbs%norb+1:4*lin%orbs%norb), stat=istat)
call memocc(istat, A, 'A', subname)


! Overlap matrix
call dgemm('t', 'n', 3*lin%orbs%norb, lin%orbs%norb, nvctrp, 1.d0, phid(lin%orbs%norb*nvctrp+1), &
     nvctrp, phid(1), nvctrp, 0.d0, ovrlp(lin%orbs%norb+1,1), 3*lin%orbs%norb) 

! MPI
call mpiallred(ovrlp(lin%orbs%norb+1,1), 3*lin%orbs%norb**2, mpi_sum, mpi_comm_world, ierr)

! The components to be projected out
call dgemm('n', 't', nvctrp, 3*lin%orbs%norb, lin%orbs%norb, 1.d0, phid(1), &
     nvctrp, ovrlp(lin%orbs%norb+1,1), 3*lin%orbs%norb, 0.d0, A(1,lin%orbs%norb+1), nvctrp)

! Project out
call daxpy(3*lin%orbs%norb*nvctrp, -1.d0, A(1,lin%orbs%norb+1), 1, phid(lin%orbs%norb*nvctrp+1), 1)



iall=-product(shape(ovrlp))*kind(ovrlp)
deallocate(ovrlp, stat=istat)
call memocc(istat, iall, 'ovrlp', subname)

allocate(ovrlp(lin%orbs%norb+1:4*lin%orbs%norb,lin%orbs%norb+1:4*lin%orbs%norb), stat=istat)
call memocc(istat, ovrlp, 'ovrlp', subname)

! Calculate the overlap matrix
call dsyrk('l', 't', 3*lin%orbs%norb, nvctrp, 1.d0, phid(lin%orbs%norb*nvctrp+1), nvctrp, &
     0.d0, ovrlp(lin%orbs%norb+1,lin%orbs%norb+1), 3*lin%orbs%norb)
! MPI
call mpiallred(ovrlp(lin%orbs%norb+1,lin%orbs%norb+1), 9*lin%orbs%norb**2, mpi_sum, mpi_comm_world, ierr)


allocate(eval(lin%orbs%norb+1:4*lin%orbs%norb), stat=istat)
call memocc(istat, eval, 'eval', subname)

! Diagonalize overlap matrix.
allocate(work(1), stat=istat)
call memocc(istat, work, 'work', subname)
call dsyev('v', 'l', 3*lin%orbs%norb, ovrlp(lin%orbs%norb+1,lin%orbs%norb+1), 3*lin%orbs%norb, eval, &
     work, -1, info)
lwork=work(1)
iall=-product(shape(work))*kind(work)
deallocate(work, stat=istat)
call memocc(istat, iall, 'work', subname)
allocate(work(lwork), stat=istat)
call memocc(istat, work, 'work', subname)
call dsyev('v', 'l', 3*lin%orbs%norb, ovrlp(lin%orbs%norb+1,lin%orbs%norb+1), 3*lin%orbs%norb, eval, &
     work, lwork, info)

! Calculate S^{-1/2}. 
! First calulate ovrlp*diag(1/sqrt(evall)) (ovrlp is the diagonalized overlap
! matrix and diag(1/sqrt(evall)) the diagonal matrix consisting of the inverse square roots of the eigenvalues...
allocate(tempArr(lin%orbs%norb+1:4*lin%orbs%norb,lin%orbs%norb+1:4*lin%orbs%norb,2), stat=istat)
call memocc(istat, tempArr, 'tempArr', subname)
do iorb=lin%orbs%norb+1,4*lin%orbs%norb
    do jorb=lin%orbs%norb+1,4*lin%orbs%norb
        tempArr(jorb,iorb,1)=ovrlp(jorb,iorb)*1.d0/sqrt(eval(iorb))
    end do
end do

! ...and now apply the diagonalized overlap matrix to the matrix constructed above.
! This will give S^{-1/2}.
call dgemm('n', 't', 3*lin%orbs%norb, 3*lin%orbs%norb, 3*lin%orbs%norb, 1.d0, ovrlp(lin%orbs%norb+1,lin%orbs%norb+1), &
     3*lin%orbs%norb, tempArr(lin%orbs%norb+1,lin%orbs%norb+1,1), 3*lin%orbs%norb, 0.d0, &
     tempArr(lin%orbs%norb+1,lin%orbs%norb+1,2), 3*lin%orbs%norb)

! Now calculate the orthonormal orbitals by applying S^{-1/2} to the orbitals.
! This requires the use of a temporary variable phidTemp.
allocate(phidTemp(nvctrp,lin%orbs%norb+1:4*lin%orbs%norb), stat=istat)
call memocc(istat, phidTemp, 'phidTemp', subname)
call dgemm('n', 'n', nvctrp, 3*lin%orbs%norb, 3*lin%orbs%norb, 1.d0, phid(lin%orbs%norb*nvctrp+1), &
     nvctrp, tempArr(lin%orbs%norb+1,lin%orbs%norb+1,2),  3*lin%orbs%norb, 0.d0, &
     phidTemp(1,lin%orbs%norb+1), nvctrp)

! Now copy the orbitals from the temporary variable to phid.
call dcopy(3*lin%orbs%norb*nvctrp, phidTemp(1,lin%orbs%norb+1), 1, phid(lin%orbs%norb*nvctrp+1), 1)

iall=-product(shape(A))*kind(A)
deallocate(A, stat=istat)
call memocc(istat, iall, 'A', subname)

iall=-product(shape(ovrlp))*kind(ovrlp)
deallocate(ovrlp, stat=istat)
call memocc(istat, iall, 'ovrlp', subname)

iall=-product(shape(work))*kind(work)
deallocate(work, stat=istat)
call memocc(istat, iall, 'work', subname)

iall=-product(shape(eval))*kind(eval)
deallocate(eval, stat=istat)
call memocc(istat, iall, 'eval', subname)

iall=-product(shape(phidTemp))*kind(phidTemp)
deallocate(phidTemp, stat=istat)
call memocc(istat, iall, 'phidTemp', subname)

iall=-product(shape(tempArr))*kind(tempArr)
deallocate(tempArr, stat=istat)
call memocc(istat, iall, 'tempArr', subname)


end subroutine orthonormalizeOnlyDerivatives








subroutine postCommunicationSumrho2(iproc, nproc, lin, sendBuf, recvBuf)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(linearParameters),intent(inout):: lin
real(8),dimension(lin%comsr%nsendBuf),intent(inout):: sendBuf
real(8),dimension(lin%comsr%nrecvBuf),intent(out):: recvBuf

! Local variables
integer:: jproc, nreceives, nsends, iorb, mpisource, istsource, ncount, lrsource, mpidest, istdest, tag, ierr
integer:: ist, istr, ilr


! Communicate the orbitals for the calculation of the charge density.
! Since we use non blocking point to point communication, only post the message
! and continues with other calculations.
! Be aware that you must not modify the send buffer without checking whether
! the communications has completed.
if(iproc==0) write(*,'(x,a)', advance='no') 'Posting sends / receives for the calculation of the charge density... '
nreceives=0
nsends=0
procLoop1: do jproc=0,nproc-1
    orbsLoop1: do iorb=1,lin%comsr%noverlaps(jproc)
        mpisource=lin%comsr%comarr(1,iorb,jproc)
        istsource=lin%comsr%comarr(2,iorb,jproc)
        ncount=lin%comsr%comarr(3,iorb,jproc)
        lrsource=lin%comsr%comarr(4,iorb,jproc)
        mpidest=lin%comsr%comarr(5,iorb,jproc)
        istdest=lin%comsr%comarr(6,iorb,jproc)
        tag=lin%comsr%comarr(7,iorb,jproc)
        if(mpisource/=mpidest) then
            ! The orbitals are on different processes, so we need a point to point communication.
            if(iproc==mpisource) then
                !write(*,'(6(a,i0))') 'process ', mpisource, ' sends ', ncount, ' elements from position ', istsource, ' to position ', istdest, ' on process ', mpidest, ', tag=',tag
                !call mpi_isend(lphi(istsource), ncount, mpi_double_precision, mpidest, tag, mpi_comm_world, lin%comsr%comarr(8,iorb,jproc), ierr)
                call mpi_isend(sendBuf(istsource), ncount, mpi_double_precision, mpidest, tag, mpi_comm_world, lin%comsr%comarr(8,iorb,jproc), ierr)
                lin%comsr%comarr(9,iorb,jproc)=mpi_request_null !is this correct?
                nsends=nsends+1
            else if(iproc==mpidest) then
                !write(*,'(6(a,i0))') 'process ', mpidest, ' receives ', ncount, ' elements at position ', istdest, ' from position ', istsource, ' on process ', mpisource, ', tag=',tag
                call mpi_irecv(recvBuf(istdest), ncount, mpi_double_precision, mpisource, tag, mpi_comm_world, lin%comsr%comarr(9,iorb,jproc), ierr)
                lin%comsr%comarr(8,iorb,jproc)=mpi_request_null !is this correct?
                nreceives=nreceives+1
            else
                lin%comsr%comarr(8,iorb,jproc)=mpi_request_null
                lin%comsr%comarr(9,iorb,jproc)=mpi_request_null
            end if
        else
            ! The orbitals are on the same process, so simply copy them.
            if(iproc==mpisource) then
                !write(*,'(6(a,i0))') 'process ', iproc, ' copies ', ncount, ' elements from position ', istsource, ' to position ', istdest, ' on process ', iproc, ', tag=',tag
                call dcopy(ncount, sendBuf(istsource), 1, recvBuf(istdest), 1)
                lin%comsr%comarr(8,iorb,jproc)=mpi_request_null
                lin%comsr%comarr(9,iorb,jproc)=mpi_request_null
                nsends=nsends+1
                nreceives=nreceives+1
                lin%comsr%communComplete(iorb,iproc)=.true.
            end if
        end if
    end do orbsLoop1
end do procLoop1
if(iproc==0) write(*,'(a)') 'done.'

if(nreceives/=lin%comsr%noverlaps(iproc)) then
    write(*,'(x,a,i0,a,i0,2x,i0)') 'ERROR on process ', iproc, ': nreceives/=lin%comsr%noverlaps(iproc)', nreceives, lin%comsr%noverlaps(iproc)
    stop
end if
call mpi_barrier(mpi_comm_world, ierr)

end subroutine postCommunicationSumrho2





!subroutine initializeCommunicationPotential(iproc, nproc, nscatterarr, lin)
subroutine initializeCommunicationPotential(iproc, nproc, nscatterarr, orbs, lzd, comgp, onWhichAtomAll)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
integer,dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
!type(linearParameters),intent(inout):: lin
type(orbitals_data),intent(in):: orbs
type(linear_zone_descriptors),intent(in):: lzd
type(p2pCommsGatherPot),intent(out):: comgp
integer,dimension(orbs%norb),intent(in):: onWhichAtomAll

! Local variables
integer:: is1, ie1, is2, ie2, is3, ie3, ilr, ii, iorb, iiorb, jproc, kproc, istat
integer:: ioverlap, is3j, ie3j, is3k, ie3k, mpidest, tag, istdest, ioffset, is3min, ie3max
integer,dimension(:,:),allocatable:: iStartEnd
character(len=*),parameter:: subname='setCommunicationPotential'


! Determine the bounds of the potential that we need for
! the orbitals on this process.
allocate(iStartEnd(6,0:nproc-1), stat=istat)
call memocc(istat, iStartEnd, 'iStartEnd', subname)
is1=0
ie1=0
is2=0
ie2=0
is3=0
ie3=0
iiorb=0
do jproc=0,nproc-1
    do iorb=1,lzd%orbs%norb_par(jproc)
        
        iiorb=iiorb+1 
        ilr=onWhichAtomAll(iiorb)
    
        ii=lzd%Llr(ilr)%nsi1
        if(ii < is1 .or. iorb==1) then
            is1=ii
        end if
        ii=lzd%Llr(ilr)%nsi1+lzd%Llr(ilr)%d%n1i
        if(ii > ie1 .or. iorb==1) then
            ie1=ii
        end if
    
        ii=lzd%Llr(ilr)%nsi2
        if(ii < is2 .or. iorb==1) then
            is2=ii
        end if
        ii=lzd%Llr(ilr)%nsi2+lzd%Llr(ilr)%d%n2i
        if(ii > ie2 .or. iorb==1) then
            ie2=ii
        end if
    
        ii=lzd%Llr(ilr)%nsi3
        if(ii < is3 .or. iorb==1) then
            is3=ii
        end if
        ii=lzd%Llr(ilr)%nsi3+lzd%Llr(ilr)%d%n3i
        if(ii > ie3 .or. iorb==1) then
            ie3=ii
        end if
    
    end do
    iStartEnd(1,jproc)=is1
    iStartEnd(2,jproc)=ie1
    iStartEnd(3,jproc)=is2
    iStartEnd(4,jproc)=ie2
    iStartEnd(5,jproc)=is3
    iStartEnd(6,jproc)=ie3
end do

! Determine how many slices each process receives.
allocate(comgp%noverlaps(0:nproc-1), stat=istat)
call memocc(istat, comgp%noverlaps, 'comgp%noverlaps', subname)
do jproc=0,nproc-1
    is3j=istartEnd(5,jproc)
    ie3j=istartEnd(6,jproc)
    mpidest=jproc
    ioverlap=0
    do kproc=0,nproc-1
        is3k=nscatterarr(kproc,3)+1
        ie3k=is3k+nscatterarr(kproc,2)-1
        if(is3j<=ie3k .and. ie3j>=is3k) then
            ioverlap=ioverlap+1
            !if(iproc==0) write(*,'(2(a,i0),a)') 'process ',jproc,' gets potential from process ',kproc,'.' 
        end if
    end do
    comgp%noverlaps(jproc)=ioverlap
    !if(iproc==0) write(*,'(2(a,i0),a)') 'Process ',jproc,' gets ',ioverlap,' potential slices.'
end do

! Determine the parameters for the communications.
allocate(comgp%overlaps(comgp%noverlaps(iproc)), stat=istat)
call memocc(istat, comgp%overlaps, 'comgp%overlaps', subname)
allocate(comgp%comarr(8,maxval(comgp%noverlaps),0:nproc-1))
call memocc(istat, comgp%comarr, 'comgp%comarr', subname)
allocate(comgp%ise3(2,0:nproc-1), stat=istat)
call memocc(istat, comgp%ise3, 'comgp%ise3', subname)
tag=1000
comgp%nrecvBuf = 0
is3min=0
ie3max=0
do jproc=0,nproc-1
    is3j=istartEnd(5,jproc)
    ie3j=istartEnd(6,jproc)
    mpidest=jproc
    ioverlap=0
    istdest=1
    do kproc=0,nproc-1
        is3k=nscatterarr(kproc,3)+1
        ie3k=is3k+nscatterarr(kproc,2)-1
        if(is3j<=ie3k .and. ie3j>=is3k) then
            ioverlap=ioverlap+1
            tag=tag+1
            is3=max(is3j,is3k) ! starting index in z dimension for data to be sent
            ie3=min(ie3j,ie3k) ! ending index in z dimension for data to be sent
            ioffset=is3-is3k ! starting index (in z direction) of data to be sent (actually it is the index -1)
            if(is3<is3min .or. ioverlap==1) then
                is3min=is3
            end if
            if(ie3>ie3max .or. ioverlap==1) then
                ie3max=ie3
            end if
            !write(*,'(a,8i8)') 'jproc, kproc, is3j, ie3j, is3k, ie3k, is3, ie3', jproc, kproc, is3j, ie3j, is3k, ie3k, is3, ie3
            call setCommunicationPotential(kproc, is3, ie3, ioffset, lzd%Glr%d%n1i, lzd%Glr%d%n2i, jproc, istdest, tag, comgp%comarr(1,ioverlap,jproc))
            !if(iproc==0) write(*,'(6(a,i0))') 'process ',comgp%comarr(1,ioverlap,jproc),' sends ',comgp%comarr(3,ioverlap,jproc),' elements from position ',&
            !                        comgp%comarr(2,ioverlap,jproc),' to position ',comgp%comarr(5,ioverlap,jproc),' on process ',&
            !                        comgp%comarr(4,ioverlap,jproc),'; tag=',comgp%comarr(6,ioverlap,jproc)
            istdest = istdest + (ie3-is3+1)*lzd%Glr%d%n1i*lzd%Glr%d%n2i
            !write(*,'(a,4i8)') 'jproc, kproc, (ie3-is3+1),lzd%Glr%d%n1i*lzd%Glr%d%n2i', jproc, kproc, (ie3-is3+1),lzd%Glr%d%n1i*lzd%Glr%d%n2i
            if(iproc==jproc) then
                comgp%nrecvBuf = comgp%nrecvBuf + (ie3-is3+1)*lzd%Glr%d%n1i*lzd%Glr%d%n2i
            end if
        end if
    end do
    comgp%ise3(1,jproc)=is3min
    comgp%ise3(2,jproc)=ie3max
    !if(iproc==0) write(*,'(a,3i8)') 'jproc, comgp%ise3(1,jproc), comgp%ise3(2,jproc)', jproc, comgp%ise3(1,jproc), comgp%ise3(2,jproc)
end do

!write(*,'(a,i4,i12)') 'iproc, comgp%nrecvBuf', iproc, comgp%nrecvBuf
allocate(comgp%recvBuf(comgp%nrecvBuf), stat=istat)
call memocc(istat, comgp%recvBuf, 'comgp%recvBuf', subname)

allocate(comgp%communComplete(maxval(comgp%noverlaps),0:nproc-1), stat=istat)
call memocc(istat, comgp%communComplete, 'comgp%communComplete', subname)


end subroutine initializeCommunicationPotential




subroutine setCommunicationPotential(mpisource, is3, ie3, ioffset, n1i, n2i, mpidest, istdest, tag, comarr)
use module_base

! Calling arguments
integer,intent(in):: mpisource, is3, ie3, ioffset, n1i, n2i, mpidest, istdest, tag
integer,dimension(8),intent(out):: comarr

! Local variables
integer:: istsource, ncount

! From which MPI process shall the slice be sent
comarr(1)=mpisource

! Starting index on the sending process
istsource=ioffset*n1i*n2i+1
comarr(2)=istsource

! Amount of data to be sent
ncount=(ie3-is3+1)*n1i*n2i
comarr(3)=ncount

! To which MPI process shall the slice be sent
comarr(4)=mpidest

! Starting index on the receiving index
comarr(5)=istdest

! Tag for the communication
comarr(6)=tag

! comarr(7): this entry is used as request for the mpi_isend.

! comarr(8): this entry is used as request for the mpi_irecv.


end subroutine setCommunicationPotential




subroutine postCommunicationsPotential(iproc, nproc, ndimpot, pot, comgp)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, ndimpot
real(8),dimension(ndimpot),intent(in):: pot
type(p2pCommsGatherPot),intent(inout):: comgp

! Local variables
integer:: jproc, kproc, nsends, nreceives, istat, mpisource, istsource, ncount, mpidest, istdest, tag, ierr


! Post the messages
if(iproc==0) write(*,'(x,a)', advance='no') 'Posting sends / receives for communicating the potential... '
nreceives=0
nsends=0
destLoop: do jproc=0,nproc-1
    sourceLoop: do kproc=1,comgp%noverlaps(jproc)
        mpisource=comgp%comarr(1,kproc,jproc)
        istsource=comgp%comarr(2,kproc,jproc)
        ncount=comgp%comarr(3,kproc,jproc)
        mpidest=comgp%comarr(4,kproc,jproc)
        istdest=comgp%comarr(5,kproc,jproc)
        tag=comgp%comarr(6,kproc,jproc)
        if(mpisource/=mpidest) then
            if(iproc==mpisource) then
                !write(*,'(6(a,i0))') 'process ', mpisource, ' sends ', ncount, ' elements from position ', istsource, ' to position ', istdest, ' on process ', mpidest, ', tag=',tag
                call mpi_isend(pot(istsource), ncount, mpi_double_precision, mpidest, tag, mpi_comm_world, comgp%comarr(7,kproc,jproc), ierr)
                comgp%comarr(8,kproc,jproc)=mpi_request_null !is this correct?
                nsends=nsends+1
            else if(iproc==mpidest) then
                !write(*,'(6(a,i0))') 'process ', mpidest, ' receives ', ncount, ' elements at position ', istdest, ' from position ', istsource, ' on process ', mpisource, ', tag=',tag
                call mpi_irecv(comgp%recvBuf(istdest), ncount, mpi_double_precision, mpisource, tag, mpi_comm_world, comgp%comarr(8,kproc,jproc), ierr)
                comgp%comarr(7,kproc,jproc)=mpi_request_null !is this correct?
                nreceives=nreceives+1
            else
                comgp%comarr(7,kproc,jproc)=mpi_request_null
                comgp%comarr(8,kproc,jproc)=mpi_request_null
            end if
        else
            ! The orbitals are on the same process, so simply copy them.
            if(iproc==mpisource) then
                !write(*,'(6(a,i0))') 'process ', iproc, ' copies ', ncount, ' elements from position ', istsource, ' to position ', istdest, ' on process ', iproc, ', tag=',tag
                call dcopy(ncount, pot(istsource), 1, comgp%recvBuf(istdest), 1)
                comgp%comarr(7,kproc,jproc)=mpi_request_null
                comgp%comarr(8,kproc,jproc)=mpi_request_null
                nsends=nsends+1
                nreceives=nreceives+1
                comgp%communComplete(kproc,iproc)=.true.
            end if
        end if
    end do sourceLoop
end do destLoop
if(iproc==0) write(*,'(a)') 'done.'

if(nreceives/=comgp%noverlaps(iproc)) then
    write(*,'(x,a,i0,a,i0,2x,i0)') 'ERROR on process ', iproc, ': nreceives/=comgp%noverlaps(iproc)', nreceives, comgp%noverlaps(iproc)
    stop
end if


end subroutine postCommunicationsPotential





subroutine gatherPotential(iproc, nproc, comgp)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(p2pCommsGatherPot),intent(inout):: comgp

! Local variables
integer:: kproc, mpisource, mpidest, nfast, nslow, nsameproc, ierr, nreceives
integer,dimension(mpi_status_size):: stat
logical:: sendComplete, receiveComplete


! Check whether the communications have completed.
comgp%communComplete=.false.
nfast=0
nsameproc=0
testLoop: do 
    do kproc=1,comgp%noverlaps(iproc)
       if(comgp%communComplete(kproc,iproc)) cycle
        call mpi_test(comgp%comarr(7,kproc,iproc), sendComplete, stat, ierr)      
        call mpi_test(comgp%comarr(8,kproc,iproc), receiveComplete, stat, ierr)   
        if(sendComplete .and. receiveComplete) comgp%communComplete(kproc,iproc)=.true.
        if(comgp%communComplete(kproc,iproc)) then
            !write(*,'(2(a,i0))') 'fast communication; process ', iproc, ' has received orbital ', korb
            mpisource=comgp%comarr(1,kproc,iproc)
            mpidest=comgp%comarr(5,kproc,iproc)
            if(mpisource/=mpidest) then
                nfast=nfast+1
            else
                nsameproc=nsameproc+1
            end if
        end if
    end do
    ! If we made it until here, either all all the communication is
    ! complete or we better wait for each single orbital.
    exit testLoop
end do testLoop


! Wait for the communications that have not completed yet
!do korb=1,nreceives
nslow=0
do kproc=1,comgp%noverlaps(iproc)
    if(comgp%communComplete(kproc,iproc)) cycle
    !write(*,'(2(a,i0))') 'process ', iproc, ' is waiting for orbital ', korb
    nslow=nslow+1
    call mpi_wait(comgp%comarr(7,kproc,iproc), stat, ierr)   !COMMENTED BY PB
    call mpi_wait(comgp%comarr(8,kproc,iproc), stat, ierr)   !COMMENTED BY PB
    comgp%communComplete(kproc,iproc)=.true.
end do

call mpiallred(nreceives, 1, mpi_sum, mpi_comm_world, ierr)
call mpiallred(nfast, 1, mpi_sum, mpi_comm_world, ierr)
call mpiallred(nslow, 1, mpi_sum, mpi_comm_world, ierr)
call mpiallred(nsameproc, 1, mpi_sum, mpi_comm_world, ierr)
if(iproc==0) write(*,'(x,2(a,i0),a)') 'statistics: - ', nfast+nslow, ' point to point communications, of which ', &
                       nfast, ' could be overlapped with computation.'
if(iproc==0) write(*,'(x,a,i0,a)') '            - ', nsameproc, ' copies on the same processor.'





end subroutine gatherPotential
