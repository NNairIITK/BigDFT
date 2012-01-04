subroutine getLinearPsi(iproc, nproc, nspin, Glr, orbs, comms, at, lin, rxyz, rxyzParab, &
    nscatterarr, ngatherarr, rhopot, GPU, input, pkernelseq, phi, psi, psit, updatePhi, &
    infoBasisFunctions, infoCoeff, itSCC, n3p, n3pi, n3d, pkernel, &
    i3s, i3xcsh, ebs, coeff, lphi, radii_cf, nlpspd, proj, communicate_lphi)
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
!use deallocatePointers
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, nspin, n3p, n3pi, n3d, i3s, i3xcsh, itSCC
type(locreg_descriptors),intent(in):: Glr
type(orbitals_data),intent(in) :: orbs
type(communications_arrays),intent(in) :: comms
type(atoms_data),intent(in):: at
type(linearParameters),intent(inout):: lin
type(input_variables),intent(in):: input
real(8),dimension(3,at%nat),intent(in):: rxyz
real(8),dimension(3,at%nat),intent(inout):: rxyzParab
integer,dimension(0:nproc-1,4),intent(inout):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
integer,dimension(0:nproc-1,2),intent(inout):: ngatherarr
real(dp),dimension(max(Glr%d%n1i*Glr%d%n2i*n3p,1)*input%nspin),intent(inout) :: rhopot
type(GPU_pointers),intent(inout):: GPU
real(dp), dimension(lin%as%size_pkernel),intent(in):: pkernel
logical,intent(in):: updatePhi
real(dp),dimension(:),pointer,intent(in):: pkernelseq
real(8),dimension(lin%lb%gorbs%npsidim),intent(inout):: phi
real(8),dimension(orbs%npsidim),intent(out):: psi, psit
integer,intent(out):: infoBasisFunctions, infoCoeff
real(8),intent(out):: ebs
real(8),dimension(lin%lb%orbs%norb,orbs%norb),intent(in out):: coeff
real(8),dimension(lin%lb%orbs%npsidim),intent(inout):: lphi
real(8),dimension(at%ntypes,3),intent(in):: radii_cf
type(nonlocal_psp_descriptors),intent(in):: nlpspd
real(wp),dimension(nlpspd%nprojel),intent(inout):: proj
logical,intent(in):: communicate_lphi

! Local variables 
integer:: istat, iall, ind1, ind2, ldim, gdim, ilr, istr, nphibuff, iorb, jorb, istart, korb, jst, nvctrp, ncount, jlr, ii
real(8),dimension(:),allocatable:: hphi, eval, lhphi, lphiold, phiold, lhphiold, hphiold, eps, temparr, work
real(8),dimension(:,:),allocatable:: HamSmall, ovrlp, ovrlpold, hamold, overlapmatrix, lambda, coeffold
real(8),dimension(:,:,:),allocatable:: matrixElements
real(8),dimension(:),pointer:: phiWork
real(8):: epot_sum, ekin_sum, eexctX, eproj_sum, trace, tt, ddot, tt2, dnrm2, t1, t2, time
character(len=*),parameter:: subname='getLinearPsi' 
logical:: withConfinement
type(workarr_sumrho):: w
integer:: ist, ierr, iiorb, info, lorb, lwork, norbtot, k, l
real(8),dimension(:),pointer:: lpot

  ! Allocate the local arrays.  
  allocate(matrixElements(lin%lb%orbs%norb,lin%lb%orbs%norb,2), stat=istat)
  call memocc(istat, matrixElements, 'matrixElements', subname)
  allocate(HamSmall(lin%lb%orbs%norb,lin%lb%orbs%norb), stat=istat)
  call memocc(istat, HamSmall, 'HamSmall', subname)
  allocate(phiWork(max(size(phi),size(psi))), stat=istat)
  call memocc(istat, phiWork, 'phiWork', subname)
  allocate(eval(lin%lb%orbs%norb), stat=istat)
  call memocc(istat, eval, 'eval', subname)
  allocate(ovrlp(lin%lb%orbs%norb,lin%lb%orbs%norb), stat=istat)
  call memocc(istat, ovrlp, 'ovrlp', subname)


  ! This is a flag whether the basis functions shall be updated.
  if(updatePhi) then
      ! If we use the derivative basis functions, the trace minimizing orbitals of the last iteration are
      ! stored in lin%lphiRestart.
      ! WARNING: Will probably not work if we use the random input guess
      if(lin%useDerivativeBasisFunctions) then
          call dcopy(lin%orbs%npsidim, lin%lphiRestart(1), 1, lphi(1), 1)
      end if

      ! Improve the trace minimizing orbitals.
      call getLocalizedBasis(iproc, nproc, at, orbs, Glr, input, lin, rxyz, nspin, &
          nscatterarr, ngatherarr, rhopot, GPU, pkernelseq, lphi, trace, rxyzParab, &
          itSCC, infoBasisFunctions, radii_cf, ovrlp, nlpspd, proj)
  end if

  ! Calculate the derivative basis functions. Copy the trace minimizing orbitals to lin%lphiRestart.
  if(lin%useDerivativeBasisFunctions .and. (updatePhi .or. itSCC==0)) then
      call dcopy(lin%orbs%npsidim, lphi(1), 1, lin%lphiRestart(1), 1)
      if(iproc==0) write(*,'(1x,a)',advance='no') 'calculating derivative basis functions...'
      call getDerivativeBasisFunctions2(iproc, nproc, input%hx, Glr, lin, lin%orbs%npsidim, lin%lphiRestart, lphi)
      if(iproc==0) write(*,'(a)') 'done.'
  end if

  ! Get the overlap matrix.
  if(.not.lin%useDerivativeBasisFunctions) then
      call getOverlapMatrix2(iproc, nproc, lin%lzd, lin%orbs, lin%comon, lin%op, lphi, lin%mad, ovrlp)
  else
      call getOverlapMatrix2(iproc, nproc, lin%lzd, lin%lb%orbs, lin%lb%comon, lin%lb%op, lphi, lin%lb%mad, ovrlp)
  end if


  if(communicate_lphi) then
      ! Allocate the communication buffers for the calculation of the charge density.
      !call allocateCommunicationbufferSumrho(iproc, lin%comsr, subname)
      ! Transform all orbitals to real space.
      ist=1
      istr=1
      do iorb=1,lin%lb%orbs%norbp
          ilr=lin%lb%orbs%inWhichLocregp(iorb)
          call initialize_work_arrays_sumrho(lin%lzd%Llr(ilr), w)
          call daub_to_isf(lin%lzd%Llr(ilr), w, lphi(ist), lin%comsr%sendBuf(istr))
          call deallocate_work_arrays_sumrho(w)
          ist = ist + lin%lzd%Llr(ilr)%wfd%nvctr_c + 7*lin%lzd%Llr(ilr)%wfd%nvctr_f
          istr = istr + lin%lzd%Llr(ilr)%d%n1i*lin%lzd%Llr(ilr)%d%n2i*lin%lzd%Llr(ilr)%d%n3i
      end do
      if(istr/=lin%comsr%nsendBuf+1) then
          write(*,'(a,i0,a)') 'ERROR on process ',iproc,' : istr/=lin%comsr%nsendBuf+1'
          stop
      end if
      
      ! Post the MPI messages for the communication of sumrho. Since we use non blocking point
      ! to point communication, the program will continue immediately. The messages will be gathered
      ! in the subroutine sumrhoForLocalizedBasis2.
      call postCommunicationSumrho2(iproc, nproc, lin, lin%comsr%sendBuf, lin%comsr%recvBuf)
  end if
  

  if(iproc==0) write(*,'(1x,a)') '----------------------------------- Determination of the orbitals in this new basis.'

  ! Gather the potential (it has been posted in the subroutine linearScaling) if the basis functions
  ! have not been updated (in that case it was gathered there).
  if(.not.updatePhi) then
      call gatherPotential(iproc, nproc, lin%comgp)
  end if
  ! If we use the derivative basis functions the potential has to be gathered anyway.
  if(lin%useDerivativeBasisFunctions) call gatherPotential(iproc, nproc, lin%lb%comgp)

  if(.not.lin%useDerivativeBasisFunctions) then
      call full_local_potential2(iproc, nproc, lin%lzd%glr%d%n1i*lin%lzd%glr%d%n2i*nscatterarr(iproc,2), &
           lin%lzd%glr%d%n1i*lin%lzd%glr%d%n2i*lin%lzd%glr%d%n3i,&
           lin%lzd%glr%d%n1i*lin%lzd%glr%d%n2i*nscatterarr(iproc,1)*input%nspin,&
           input%nspin, lin%orbs,lin%lzd, ngatherarr, rhopot, lpot, 2, lin%comgp)
  else
      call full_local_potential2(iproc, nproc, lin%lzd%glr%d%n1i*lin%lzd%glr%d%n2i*nscatterarr(iproc,2), &
           lin%lzd%glr%d%n1i*lin%lzd%glr%d%n2i*lin%lzd%glr%d%n3i,&
           lin%lzd%glr%d%n1i*lin%lzd%glr%d%n2i*nscatterarr(iproc,1)*input%nspin,&
           input%nspin, lin%lb%orbs,lin%lzd, ngatherarr, rhopot, lpot, 2, lin%lb%comgp)
  end if

  ! Apply the Hamitonian to the orbitals. The flag withConfinement=.false. indicates that there is no
  ! confining potential added to the Hamiltonian.
  allocate(lhphi(lin%lb%orbs%npsidim), stat=istat)
  call memocc(istat, lhphi, 'lhphi', subname)
  withConfinement=.false.
  if(iproc==0) write(*,'(1x,a)',advance='no') 'Hamiltonian application...'
  allocate(lin%lzd%doHamAppl(lin%lzd%nlr), stat=istat)
  call memocc(istat, lin%lzd%doHamAppl, 'lin%lzd%doHamAppl', subname)
  lin%lzd%doHamAppl=.true.
  if(.not.lin%useDerivativeBasisFunctions) then
      call HamiltonianApplication3(iproc, nproc, at, lin%orbs, input%hx, input%hy, input%hz, rxyz, &
           proj, lin%lzd, ngatherarr, lpot, lphi, lhphi, &
           ekin_sum, epot_sum, eexctX, eproj_sum, nspin, GPU, withConfinement, .true., pkernel=pkernelseq)
  else
      !! ATTENTION NEW!!
      ! Modify the value of lzd%lnpsidimtot to take into account the derivatives
      ii=lin%lzd%lpsidimtot
      lin%lzd%lpsidimtot=lin%lzd%lpsidimtot_der

      ! Deallocate old PSP structures and rebuild them for the derivatives.
      do ilr=1,lin%lzd%nlr
          call deallocate_nonlocal_psp_descriptors(lin%lzd%lnlpspd(ilr), subname)
          if(associated(lin%lzd%llr(ilr)%projflg)) then
              if(size(lin%lzd%llr(ilr)%projflg)>0) then
                  iall=-product(shape(lin%lzd%llr(ilr)%projflg))*kind(lin%lzd%llr(ilr)%projflg)
                  deallocate(lin%lzd%llr(ilr)%projflg, stat=istat)
                  call memocc(istat, iall, 'lin%lzd%llr(ilr)%projflg', subname)
              else
                  nullify(lin%lzd%llr(ilr)%projflg)
              end if
          end if
      end do
      call prepare_lnlpspd(iproc, at, input, lin%lb%orbs, rxyz, radii_cf, lin%locregShape, lin%lzd)

      call HamiltonianApplication3(iproc, nproc, at, lin%lb%orbs, input%hx, input%hy, input%hz, rxyz, &
           proj, lin%lzd, ngatherarr, lpot, lphi, lhphi, &
           ekin_sum, epot_sum, eexctX, eproj_sum, nspin, GPU, withConfinement, .true., pkernel=pkernelseq)
      lin%lzd%lpsidimtot=ii
      
      ! Deallocate the derivative PSP structures and rebuild them for the standard case (i.e. without derivatives)
      do ilr=1,lin%lzd%nlr
          call deallocate_nonlocal_psp_descriptors(lin%lzd%lnlpspd(ilr), subname)
          if(associated(lin%lzd%llr(ilr)%projflg)) then
              if(size(lin%lzd%llr(ilr)%projflg)>0) then
                  iall=-product(shape(lin%lzd%llr(ilr)%projflg))*kind(lin%lzd%llr(ilr)%projflg)
                  deallocate(lin%lzd%llr(ilr)%projflg, stat=istat)
                  call memocc(istat, iall, 'lin%lzd%llr(ilr)%projflg', subname)
              else
                  nullify(lin%lzd%llr(ilr)%projflg)
              end if
          end if
      end do
      call prepare_lnlpspd(iproc, at, input, lin%orbs, rxyz, radii_cf, lin%locregShape, lin%lzd)
  end if
  iall=-product(shape(lin%lzd%doHamAppl))*kind(lin%lzd%doHamAppl)
  deallocate(lin%lzd%doHamAppl, stat=istat)
  call memocc(istat, iall, 'lin%lzd%doHamAppl', subname)

  if(.not.lin%useDerivativeBasisFunctions)then
     iall=-product(shape(lin%orbs%ispot))*kind(lin%orbs%ispot)
     deallocate(lin%orbs%ispot, stat=istat)
     call memocc(istat, iall, 'lin%orbs%ispot', subname)
  else
     iall=-product(shape(lin%lb%orbs%ispot))*kind(lin%lb%orbs%ispot) 
     deallocate(lin%lb%orbs%ispot, stat=istat)
     call memocc(istat, iall, 'lin%lb%orbs%ispot', subname)
  end if


  iall=-product(shape(lpot))*kind(lpot)
  deallocate(lpot, stat=istat)
  call memocc(istat, iall, 'lpot', subname)

  if(iproc==0) write(*,'(1x,a)') 'done.'

  ! Deallocate the buffers needed for the communication of the potential.
  call deallocateCommunicationsBuffersPotential(lin%comgp, subname)
  if(lin%useDerivativeBasisFunctions) call deallocateCommunicationsBuffersPotential(lin%lb%comgp, subname)

  ! Calculate the matrix elements <phi|H|phi>.
  call allocateCommuncationBuffersOrtho(lin%lb%comon, subname)
  if(.not. lin%useDerivativeBasisFunctions) then
      call getMatrixElements2(iproc, nproc, lin%lzd, lin%lb%orbs, lin%lb%op, lin%lb%comon, lphi, lhphi, lin%mad, matrixElements)
  else
      call getMatrixElements2(iproc, nproc, lin%lzd, lin%lb%orbs, lin%lb%op, lin%lb%comon, lphi, lhphi, lin%lb%mad, matrixElements)
  end if
  call deallocateCommuncationBuffersOrtho(lin%lb%comon, subname)

  tt=0.d0
  do iall=1,lin%lb%orbs%norb
       tt=tt+matrixElements(iall,iall,1)
  end do
  if(iproc==0) write(*,*) 'trace of H without confinement:',tt


  allocate(overlapmatrix(lin%lb%orbs%norb,lin%lb%orbs%norb), stat=istat)
  call memocc(istat, overlapmatrix, 'overlapmatrix', subname)
  overlapmatrix=ovrlp
  

  ! Diagonalize the Hamiltonian, either iteratively or with lapack.
  call mpi_barrier(mpi_comm_world, ierr) !To measure the time correctly.
  call cpu_time(t1)
  if(trim(lin%getCoeff)=='min') then
      call optimizeCoefficients(iproc, orbs, lin, nspin, matrixElements, coeff, infoCoeff)
  else if(trim(lin%getCoeff)=='diag') then
      ! Make a copy of the matrix elements since dsyev overwrites the matrix and the matrix elements
      ! are still needed later.
      call dcopy(lin%lb%orbs%norb**2, matrixElements(1,1,1), 1, matrixElements(1,1,2), 1)
      if(lin%blocksize_pdsyev<0) then
          if(iproc==0) write(*,'(1x,a)',advance='no') 'Diagonalizing the Hamiltonian, sequential version... '
          call diagonalizeHamiltonian2(iproc, nproc, lin%lb%orbs, matrixElements(1,1,2), ovrlp, eval)
      else
          if(iproc==0) write(*,'(1x,a)',advance='no') 'Diagonalizing the Hamiltonian, parallel version... '
          call dsygv_parallel(iproc, nproc, lin%blocksize_pdsyev, lin%nproc_pdsyev, mpi_comm_world, 1, 'v', 'l', lin%lb%orbs%norb,&
               matrixElements(1,1,2), lin%lb%orbs%norb, ovrlp, lin%lb%orbs%norb, eval, info)
      end if
      if(iproc==0) write(*,'(a)') 'done.'
      call dcopy(lin%lb%orbs%norb*orbs%norb, matrixElements(1,1,2), 1, coeff(1,1), 1)
      if(.not.updatePhi) call dcopy(lin%lb%orbs%norb*(orbs%norb+lin%norbvirt), matrixElements(1,1,2), 1, lin%coeffall(1,1), 1)
      infoCoeff=0

      ! Write some eigenvalues. Don't write all, but only a few around the last occupied orbital.
      if(iproc==0) then
          do iorb=max(orbs%norb-8,1),min(orbs%norb+8,lin%orbs%norb)
              if(iorb==orbs%norb) then
                  write(*,'(1x,a,i0,a,es12.5,a)') 'eval(',iorb,')=',eval(iorb),'  <-- last occupied orbital'
              else if(iorb==orbs%norb+1) then
                  write(*,'(1x,a,i0,a,es12.5,a)') 'eval(',iorb,')=',eval(iorb),'  <-- first virtual orbital'
              else
                  write(*,'(1x,a,i0,a,es12.5)') 'eval(',iorb,')=',eval(iorb)
              end if
          end do
      end if
  else if(trim(lin%getCoeff)=='new') then
      !stop 'not yet ready'
      !! THIS IS THE NEW PART ###########################################################################
      call getCoefficients_new(iproc, nproc, lin, orbs, lin%hamold, lphi, ovrlp, coeff)
      !! ################################################################################################
  end if
  call cpu_time(t2)
  time=t2-t1
  if(iproc==0) write(*,'(1x,a,es10.3)') 'time for diagonalizing the Hamiltonian:',time


  ! Calculate the band structure energy with matrixElements.
  ebs=0.d0
  do iorb=1,orbs%norb
      do jorb=1,lin%lb%orbs%norb
          do korb=1,lin%lb%orbs%norb
              ebs = ebs + coeff(jorb,iorb)*coeff(korb,iorb)*matrixElements(korb,jorb,1)
          end do
      end do
  end do
  ! If closed shell multiply by two.
  if(input%nspin==1) ebs=2.d0*ebs
  

  ! Copy the basis functions for the next iterations
  call dcopy(lin%orbs%npsidim, lphi(1), 1, lin%lphiold(1), 1)

  ! Copy the Hamiltonian matrix for the next iteration
  call dcopy(lin%orbs%norb**2, matrixElements(1,1,1), 1, lin%hamold(1,1), 1)

  ! Deallocate all local arrays.
  iall=-product(shape(HamSmall))*kind(HamSmall)
  deallocate(HamSmall, stat=istat)
  call memocc(istat, iall, 'HamSmall', subname)

  iall=-product(shape(lhphi))*kind(lhphi)
  deallocate(lhphi, stat=istat)
  call memocc(istat, iall, 'lhphi', subname)

  iall=-product(shape(phiWork))*kind(phiWork)
  deallocate(phiWork, stat=istat)
  call memocc(istat, iall, 'phiWork', subname)

  iall=-product(shape(matrixElements))*kind(matrixElements)
  deallocate(matrixElements, stat=istat)
  call memocc(istat, iall, 'matrixElements', subname)
  
  iall=-product(shape(eval))*kind(eval)
  deallocate(eval, stat=istat)
  call memocc(istat, iall, 'eval', subname)

  iall=-product(shape(ovrlp))*kind(ovrlp)
  deallocate(ovrlp, stat=istat)
  call memocc(istat, iall, 'ovrlp', subname)

  iall=-product(shape(overlapmatrix))*kind(overlapmatrix)
  deallocate(overlapmatrix, stat=istat)
  call memocc(istat, iall, 'overlapmatrix', subname)



end subroutine getLinearPsi








subroutine getLocalizedBasis(iproc, nproc, at, orbs, Glr, input, lin, rxyz, nspin, &
    nscatterarr, ngatherarr, rhopot, GPU, pkernelseq, lphi, trH, rxyzParabola, &
    itScc, infoBasisFunctions, radii_cf, ovrlp, nlpspd, proj)
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
integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
real(dp), dimension(*), intent(inout) :: rhopot
type(GPU_pointers), intent(inout) :: GPU
real(dp), dimension(:), pointer :: pkernelseq
real(8),dimension(lin%orbs%npsidim):: lphi
real(8):: trH
real(8),dimension(at%ntypes,3),intent(in):: radii_cf
real(8),dimension(lin%orbs%norb,lin%orbs%norb),intent(out):: ovrlp
type(nonlocal_psp_descriptors),intent(in):: nlpspd
real(wp),dimension(nlpspd%nprojel),intent(inout):: proj

! Local variables
real(8) ::epot_sum, ekin_sum, eexctX, eproj_sum, evalmax, eval_zero, t1tot
real(8) :: t2tot, timetot, tt1, tt2, tt3, tt4, tt5, lstep, dfactorial
real(8):: tt, ddot, fnrm, fnrmMax, meanAlpha, gnrm, gnrm_zero, gnrmMax, t1, t2
real(8) :: timecommunp2p, timeextract, timecommuncoll, timeoverlap, timecompress, energyconf_0, energyconf_trial
real(8):: trHold
integer:: iorb, icountSDSatur, icountSwitch, idsx, icountDIISFailureTot, itinner, consecutive_rejections
integer :: icountDIISFailureCons, itBest, info, lwork, ndim_lchi, ndim_lhchi
integer:: istat, istart, ierr, ii, it, iall, nit, ind1, ind2, jorb, i, ist, jst, iiorb, jjorb, ilrold, k
integer:: ldim, gdim, ilr, ncount, offset, istsource, istdest
real(8),dimension(:),allocatable:: alpha, fnrmOldArr, alphaDIIS, lhphi, lhphiold
real(8),dimension(:),allocatable:: eval, lvphi, lvphiovrlp, alpha2, lhpsiold, work, rwork, lphiold, lphiovrlporig, lphi2
real(8),dimension(:),allocatable:: lchi, lphidebug
real(8),dimension(:,:),allocatable:: HamSmall, fnrmArr, fnrmOvrlpArr, W, ttmat, Kmat, Gmat, Umat, lhchi, omatr, omat2
real(8),dimension(:,:,:),allocatable:: Gmatc, tempmat, Omat, tempmat2, ham3, tempmat3, omatrtot
logical:: quiet, allowDIIS, startWithSD, withConfinement, calc
character(len=*),parameter:: subname='getLocalizedBasis'
character(len=1):: message
type(localizedDIISParameters):: ldiis
real(8),dimension(5):: time
real(8),dimension(:),pointer:: potential
real(8),dimension(:),pointer:: phiWork
real(8),dimension(:),pointer:: lpot
real(8), external :: mpi_wtime1
complex(8):: ttc, ttc2, ttc3, zdotu
complex(8),dimension(:),allocatable:: Gmat_c, expD_cmplx
character(len=3):: orbname, comment
integer,dimension(:),allocatable:: onwhichatomtemp, norb_parTemp, onWhichMPITemp
logical,dimension(:),allocatable:: doNotCalculate, skip
integer:: iat, is1, ie1, is2, ie2, is3, ie3, jlr, js1, je1, js2, je2, js3, je3
integer :: norbTarget, nprocTemp, kk, jlrold, nlocregPerMPI, tag, jproc
logical:: ovrlpx, ovrlpy, ovrlpz, check_whether_locregs_overlap, resetDIIS, immediateSwitchToSD




  ! Allocate all local arrays.
  call allocateLocalArrays()
  
  
  if(iproc==0) write(*,'(1x,a)') '======================== Creation of the basis functions... ========================'

  ! Initialize the arrays and variable needed for DIIS.
  call initializeDIIS(lin%DIISHistMax, lin%lzd, lin%orbs, lin%orbs%inWhichLocregp, lin%startWithSD, lin%alphaSD, lin%alphaDIIS, &
       lin%orbs%norb, icountSDSatur, icountSwitch, icountDIISFailureTot, icountDIISFailureCons, allowDIIS, &
       startWithSD, ldiis, alpha, alphaDIIS)

  ! Set the maximal number of iterations.
  if(itSCC==1) then
      nit=lin%nItBasisFirst
  else
      nit=lin%nItBasis
  end if

  ! Gather the potential that each process needs for the Hamiltonian application for all its orbitals.
  ! The messages for this point to point communication have been posted in the subroutine linearScaling.
  call gatherPotential(iproc, nproc, lin%comgp)

  ! Build the required potential
  call full_local_potential2(iproc, nproc, lin%lzd%glr%d%n1i*lin%lzd%glr%d%n2i*nscatterarr(iproc,2), &
       lin%lzd%glr%d%n1i*lin%lzd%glr%d%n2i*lin%lzd%glr%d%n3i,&
       lin%lzd%glr%d%n1i*lin%lzd%glr%d%n2i*nscatterarr(iproc,1)*input%nspin,&
       input%nspin, lin%orbs,lin%lzd, ngatherarr, rhopot, lpot, 2, lin%comgp)

  allocate(lphiold(size(lphi)), stat=istat)
  call memocc(istat, lphiold, 'lphiold', subname)
  allocate(W(lin%orbs%norb,lin%orbs%norb), stat=istat)
  call memocc(istat, W, 'W', subname)
  allocate(eval(lin%orbs%norb), stat=istat)
  call memocc(istat, eval, 'eval', subname)

  time=0.d0
  resetDIIS=.false.
  immediateSwitchToSD=.false.
  t1tot=mpi_wtime()
  iterLoop: do it=1,nit
      fnrmMax=0.d0
      fnrm=0.d0
  
      if (iproc==0) then
          write( *,'(1x,a,i0)') repeat('-',77 - int(log(real(it))/log(10.))) // ' iter=', it
      endif


      !! Flatten at the edges -  EXPERIMENTAL
      !call flatten_at_edges(iproc, nproc, lin, at, input, lin%orbs, lin%lzd, rxyz, lphi)

  
      ! Orthonormalize the orbitals. If the localization regions are smaller that the global box (which
      ! will be the usual case), the orthogonalization can not be done exactly, but only approximately.
      if(iproc==0) then
          write(*,'(1x,a)') 'Orthonormalization... '
      end if
      t1=mpi_wtime()
      if(.not.ldiis%switchSD) call orthonormalizeLocalized(iproc, nproc, lin%methTransformOverlap, lin%nItOrtho, &
           lin%blocksize_pdsyev, lin%blocksize_pdgemm, lin%orbs, lin%op, lin%comon, lin%lzd, lin%orbs%inWhichLocreg, &
           lin%convCritOrtho, input, lin%mad, lphi, ovrlp, 'new')
      t2=mpi_wtime()
      time(1)=time(1)+t2-t1
      !!! Flatten at the edges -  EXPERIMENTAL
      !!call flatten_at_edges(iproc, nproc, lin, at, input, lin%orbs, lin%lzd, rxyz, lphi)


      !!! NEW #########################################################
      t1=mpi_wtime()
      if(.not.ldiis%switchSD) then
          call unitary_optimization(iproc, nproc, lin, lin%lzd, lin%orbs, at, input, lin%op, &
                                    lin%comon, rxyz, lin%nItInnerLoop, lphi)
          !!! Flatten at the edges -  EXPERIMENTAL
          !!call flatten_at_edges(iproc, nproc, lin, at, input, lin%orbs, lin%lzd, rxyz, lphi)
      end if
      t2=mpi_wtime()
      time(5)=time(5)+t2-t1
      !!!call unitary_optimization2(iproc, nproc, lin, lin%lzd, lin%orbs, at, input, lin%op, lin%comon, rxyz, lin%nItInnerLoop, lphi)
      !!! #############################################################

      ! Post the sends again to calculate the overlap matrix (will be needed for the orthoconstraint).
      call allocateSendBufferOrtho(lin%comon, subname)
      call allocateRecvBufferOrtho(lin%comon, subname)
      ! Extract the overlap region from the orbitals phi and store them in comon%sendBuf.
      call extractOrbital3(iproc, nproc, lin%orbs, lin%orbs%npsidim, lin%orbs%inWhichLocreg, &
           lin%lzd, lin%op, lphi, lin%comon%nsendBuf, lin%comon%sendBuf)
      ! Post the send messages.
      call postCommsOverlapNew(iproc, nproc, lin%orbs, lin%op, lin%lzd, lphi, lin%comon, timecommunp2p, timeextract)

  
      ! Calculate the unconstrained gradient by applying the Hamiltonian.
      if(iproc==0) then
          write(*,'(1x,a)', advance='no') 'Hamiltonian application... '
      end if
      t1=mpi_wtime()
      !withConfinement=.false.
      withConfinement=.true.
      allocate(lin%lzd%doHamAppl(lin%orbs%norb), stat=istat)
      call memocc(istat, lin%lzd%doHamAppl, 'lin%lzd%doHamAppl', subname)
      lin%lzd%doHamAppl=.true.


      call HamiltonianApplication3(iproc, nproc, at, lin%orbs, input%hx, input%hy, input%hz, rxyz, &
           proj, lin%lzd, ngatherarr, lpot, lphi, lhphi, &
           ekin_sum, epot_sum, eexctX, eproj_sum, nspin, GPU, withConfinement, .true., &
           pkernel=pkernelseq, lin=lin, confinementCenter=lin%orbs%inWhichLocregp)

      iall=-product(shape(lin%lzd%doHamAppl))*kind(lin%lzd%doHamAppl)
      deallocate(lin%lzd%doHamAppl, stat=istat)
      call memocc(istat, iall, 'lin%lzd%doHamAppl', subname)

      t2=mpi_wtime()


      time(2)=time(2)+t2-t1
      if(iproc==0) then
          write(*,'(a)', advance='no') 'done. '
      end if


  
      ! Apply the orthoconstraint to the gradient. This subroutine also calculates the trace trH.
      if(iproc==0) then
          write(*,'(a)', advance='no') 'Orthoconstraint... '
      end if
      !!!!!!!!!!!!!!!call orthoconstraintLocalized(iproc, nproc, lin, input, lphi, lhphi, trH)

      ! Gather the messages and calculate the overlap matrix.
      call collectnew(iproc, nproc, lin%comon, lin%mad, lin%op, lin%orbs, input, lin%lzd, lin%comon%nsendbuf, &
           lin%comon%sendbuf, lin%comon%nrecvbuf, lin%comon%recvbuf, ovrlp, timecommunp2p, timecommuncoll, timecompress)
      call calculateOverlapMatrix3(iproc, nproc, lin%orbs, lin%op, lin%orbs%inWhichLocreg, lin%comon%nsendBuf, &
           lin%comon%sendBuf, lin%comon%nrecvBuf, lin%comon%recvBuf, lin%mad, ovrlp)
      call deallocateRecvBufferOrtho(lin%comon, subname)
      call deallocateSendBufferOrtho(lin%comon, subname)

      t1=mpi_wtime()
      !!call flatten_at_edges(iproc, nproc, lin, at, input, lin%orbs, lin%lzd, rxyz, lhphi)
      call orthoconstraintNonorthogonal(iproc, nproc, lin, input, ovrlp, lphi, lhphi, lin%mad, trH, W, eval)


      ! Flatten at the edges -  EXPERIMENTAL
      !!call flatten_at_edges(iproc, nproc, lin, at, input, lin%orbs, lin%lzd, rxyz, lhphi)

      ! Cycle if the trace increased (steepest descent only)
      if(iproc==0) write(*,*) 'ldiis%switchSD, ldiis%isx', ldiis%switchSD, ldiis%isx
      if(.not. ldiis%switchSD .and. ldiis%isx==0) then
           if(iproc==0) write(*,*) 'trH, trHold', trH, trHold
           if(trH>trHold) then
               consecutive_rejections=consecutive_rejections+1
               if(consecutive_rejections==300000) then
                   if(fnrmMax<lin%convCrit .or. it>=nit) then
                       if(it>=nit) then
                           if(iproc==0) write(*,'(1x,a,i0,a)') 'WARNING: not converged within ', it, &
                               ' iterations! Exiting loop due to consective failures of SD.'
                           if(iproc==0) write(*,'(1x,a,2es15.7,f12.7)') 'CHECK THIS Final values for fnrm, fnrmMax, trace: ', &
                                        fnrm, fnrmMax, trHold
                           infoBasisFunctions=-1
                       else
                           if(iproc==0) then
                               write(*,'(1x,a,i0,a,2es15.7,f12.7)') 'converged in ', it, ' iterations.'
                               write (*,'(1x,a,2es15.7,f12.7)') 'CHECK THIS Final values for fnrm, fnrmMax, trace: ', &
                               fnrm, fnrmMax, trHold
                           end if
                           infoBasisFunctions=it
                       end if
                       if(iproc==0) write(*,'(1x,a)') '============================= Basis functions created. &
                                    &============================='
                       !!if(lin%plotBasisFunctions) then
                       !!    call plotOrbitals(iproc, lin%orbs, Glr, phi, at%nat, rxyz, lin%onWhichAtom, .5d0*input%hx, &
                       !!        .5d0*input%hy, .5d0*input%hz, 1)
                       !!end if
                       exit iterLoop
                   end if
               end if
               alpha=alpha*.5d0
               call dcopy(size(lphi), lphiold, 1, lphi, 1)
               if(iproc==0) write(*,'(x,a)') 'trace increased; reject orbitals and cycle'
               cycle iterLoop
           end if
      end if
      consecutive_rejections=0


      t2=mpi_wtime()
      time(3)=time(3)+t2-t1
      if(iproc==0) then
          write(*,'(a)', advance='no') 'done. '
      end if
      if(iproc==0) write(*,*) 'trH',trH


  
      ! Calculate the norm of the gradient (fnrmArr) and determine the angle between the current gradient and that
      ! of the previous iteration (fnrmOvrlpArr).
      istart=1
      do iorb=1,lin%orbs%norbp
          ilr=lin%orbs%inWhichLocregp(iorb)
          ncount=lin%lzd%llr(ilr)%wfd%nvctr_c+7*lin%lzd%llr(ilr)%wfd%nvctr_f
          if(it>1) fnrmOvrlpArr(iorb,1)=ddot(ncount, lhphi(istart), 1, lhphiold(istart), 1)
          fnrmArr(iorb,1)=ddot(ncount, lhphi(istart), 1, lhphi(istart), 1)
          istart=istart+ncount
      end do

      ! Keep the gradient for the next iteration.
      if(it>1) then
          call dcopy(lin%orbs%norbp, fnrmArr(1,1), 1, fnrmOldArr(1), 1)
      end if
  
      ! Determine the gradient norm and its maximal component. In addition, adapt the
      ! step size for the steepest descent minimization (depending on the angle 
      ! between the current gradient and the one from the previous iteration).
      ! This is of course only necessary if we are using steepest descent and not DIIS.
      do iorb=1,lin%orbs%norbp
          fnrm=fnrm+fnrmArr(iorb,1)
          if(fnrmArr(iorb,1)>fnrmMax) fnrmMax=fnrmArr(iorb,1)
          if(it>1 .and. ldiis%isx==0 .and. .not.ldiis%switchSD) then
          ! Adapt step size for the steepest descent minimization.
              tt=fnrmOvrlpArr(iorb,1)/sqrt(fnrmArr(iorb,1)*fnrmOldArr(iorb))
              !if(tt>.7d0) then
              if(tt>.9d0 .and. trH<trHold) then
                  alpha(iorb)=alpha(iorb)*1.05d0
              else
                  alpha(iorb)=alpha(iorb)*.5d0
              end if
          end if
      end do
      call mpiallred(fnrm, 1, mpi_sum, mpi_comm_world, ierr)
      call mpiallred(fnrmMax, 1, mpi_max, mpi_comm_world, ierr)
      !fnrm=sqrt(fnrm)
      fnrm=sqrt(fnrm/dble(lin%orbs%norb))
      fnrmMax=sqrt(fnrmMax)
      ! Copy the gradient (will be used in the next iteration to adapt the step size).
      call dcopy(lin%orbs%npsidim, lhphi, 1, lhphiold, 1)
      trHold=trH
  

      ! Precondition the gradient.
      if(iproc==0) then
          write(*,'(a)', advance='no') 'Preconditioning... '
      end if
      gnrm=1.d3 ; gnrm_zero=1.d3
      t1=mpi_wtime()

      ind2=1
      do iorb=1,lin%orbs%norbp
          ilr = lin%orbs%inWhichLocregp(iorb)
          call choosePreconditioner2(iproc, nproc, lin%orbs, lin, lin%lzd%Llr(ilr), input%hx, input%hy, input%hz, &
              lin%nItPrecond, lhphi(ind2), at%nat, rxyz, at, it, iorb, eval_zero)
          ind2=ind2+lin%lzd%Llr(ilr)%wfd%nvctr_c+7*lin%lzd%Llr(ilr)%wfd%nvctr_f
      end do
      !!!end if
      !call preconditionall(iproc, nproc, lin%orbs, lin%lzd%glr, input%hx, input%hy, input%hz, lin%nItPrecond, lhphi, tt, tt2)
      t2=mpi_wtime()
      time(4)=time(4)+t2-t1
      if(iproc==0) then
          write(*,'(a)') 'done. '
      end if


      !!!! plot the orbitals -- EXPERIMENTAL ##################################################
      !!!allocate(lvphiovrlp(lin%lzd%glr%wfd%nvctr_c+7*lin%lzd%glr%wfd%nvctr_f))
      !!!ist=1
      !!!write(comment,'(i3.3)') it
      !!!do iorb=1,lin%orbs%norbp
      !!!    iiorb=iorb+lin%orbs%isorb
      !!!    ilr=lin%orbs%inwhichlocreg(iiorb)
      !!!    write(orbname,'(i3.3)') iiorb
      !!!    write(*,'(a,i0)') 'plotting orbital ',iiorb
      !!!    lvphiovrlp=0.d0
      !!!    call Lpsi_to_global2(iproc, nproc, lin%lzd%llr(ilr)%wfd%nvctr_c+7*lin%lzd%llr(ilr)%wfd%nvctr_f, &
      !!!         lin%lzd%glr%wfd%nvctr_c+7*lin%lzd%glr%wfd%nvctr_f, lin%orbs%norb, lin%orbs%nspinor, input%nspin, &
      !!!         lin%lzd%Glr, lin%lzd%Llr(ilr), lphi(ist), lvphiovrlp(1))
      !!!    call plot_wf(orbname//'_'//comment, 2, at, 1.d0, lin%lzd%glr, input%hx, input%hy, input%hz, rxyz, lvphiovrlp(1))
      !!!    ist=ist+lin%lzd%llr(ilr)%wfd%nvctr_c+7*lin%lzd%llr(ilr)%wfd%nvctr_f
      !!!end do
      !!!deallocate(lvphiovrlp)
      !!!! ####################################################################################





      ! Determine the mean step size for steepest descent iterations.
      tt=sum(alpha)
      meanAlpha=tt/dble(lin%orbs%norb)
  
      ! Write some informations to the screen.
      if(iproc==0) write(*,'(1x,a,i6,2es15.7,f17.10)') 'iter, fnrm, fnrmMax, trace', it, fnrm, fnrmMax, trH
      if(fnrmMax<lin%convCrit .or. it>=nit) then
          if(it>=nit) then
              if(iproc==0) write(*,'(1x,a,i0,a)') 'WARNING: not converged within ', it, &
                  ' iterations! Exiting loop due to limitations of iterations.'
              if(iproc==0) write(*,'(1x,a,2es15.7,f12.7)') 'Final values for fnrm, fnrmMax, trace: ', fnrm, fnrmMax, trH
              infoBasisFunctions=-1
          else
              if(iproc==0) then
                  write(*,'(1x,a,i0,a,2es15.7,f12.7)') 'converged in ', it, ' iterations.'
                  write (*,'(1x,a,2es15.7,f12.7)') 'Final values for fnrm, fnrmMax, trace: ', fnrm, fnrmMax, trH
              end if
              infoBasisFunctions=it
          end if
          if(iproc==0) write(*,'(1x,a)') '============================= Basis functions created. ============================='
          !!if(lin%plotBasisFunctions) then
          !!    call plotOrbitals(iproc, lin%orbs, Glr, phi, at%nat, rxyz, lin%onWhichAtom, .5d0*input%hx, &
          !!        .5d0*input%hy, .5d0*input%hz, 1)
          !!end if
          exit iterLoop
      end if
  
  
      ! Determine whether the basis functions shall be further optimized using DIIS or steepest descent.
      call DIISorSD()
      if(iproc==0) then
          if(ldiis%isx>0) then
              write(*,'(1x,3(a,i0))') 'DIIS informations: history length=',ldiis%isx, ', consecutive failures=', &
                  icountDIISFailureCons, ', total failures=', icountDIISFailureTot
          else
              if(allowDIIS) then
                  message='y'
              else
                  message='n'
              end if
              write(*,'(1x,a,es9.3,a,i0,a,a)') 'steepest descent informations: mean alpha=', meanAlpha, &
              ', consecutive successes=', icountSDSatur, ', DIIS=', message
          end if
      end if


      ! Improve the orbitals, depending on the choice made above.
      if(.not.ldiis%switchSD) then
          call improveOrbitals()
      else
          if(iproc==0) write(*,'(x,a)') 'no improvement of the orbitals, recalculate gradient'
      end if




      !!! plot the orbitals -- EXPERIMENTAL ##################################################
      !!allocate(lvphiovrlp(lin%lzd%glr%wfd%nvctr_c+7*lin%lzd%glr%wfd%nvctr_f))
      !!ist=1
      !!write(comment,'(i3.3)') it
      !!do iorb=1,lin%orbs%norbp
      !!    iiorb=iorb+lin%orbs%isorb
      !!    ilr=lin%orbs%inwhichlocreg(iiorb)
      !!    lvphiovrlp=0.d0
      !!    call Lpsi_to_global2(iproc, nproc, lin%lzd%llr(ilr)%wfd%nvctr_c+7*lin%lzd%llr(ilr)%wfd%nvctr_f, &
      !!         lin%lzd%glr%wfd%nvctr_c+7*lin%lzd%glr%wfd%nvctr_f, lin%orbs%norb, lin%orbs%nspinor, input%nspin, &
      !!         lin%lzd%Glr, lin%lzd%Llr(ilr), lphi(ist), lvphiovrlp(1))
      !!    call plotOrbitals(iproc, orbs, lin%lzd%Glr, lvphiovrlp, at%nat, rxyz, lin%orbs%inwhichlocreg, &
      !!         .5d0*input%hx, .5d0*input%hy, .5d0*input%hz, it)
      !!    ist=ist+lin%lzd%llr(ilr)%wfd%nvctr_c+7*lin%lzd%llr(ilr)%wfd%nvctr_f
      !!end do
      !!deallocate(lvphiovrlp)
      !!! ####################################################################################




     ! Flush the standard output
     !flush(unit=6) 

  end do iterLoop


  iall=-product(shape(lphiold))*kind(lphiold)
  deallocate(lphiold, stat=istat)
  call memocc(istat, iall, 'lphiold', subname)
  iall=-product(shape(W))*kind(W)
  deallocate(W, stat=istat)
  call memocc(istat, iall, 'W', subname)
  iall=-product(shape(eval))*kind(eval)
  deallocate(eval, stat=istat)
  call memocc(istat, iall, 'eval', subname)


  iall=-product(shape(lin%orbs%ispot))*kind(lin%orbs%ispot)
  deallocate(lin%orbs%ispot, stat=istat)
  call memocc(istat, iall, 'lin%orbs%ispot', subname)
  ! Deallocate potential
  iall=-product(shape(lpot))*kind(lpot)
  deallocate(lpot, stat=istat)
  call memocc(istat, iall, 'lpot', subname)

  ! Deallocate PSP stuff
  !call free_lnlpspd(lin%orbs, lin%lzd)

  t2tot=mpi_wtime()
  timetot=t2tot-t1tot

  ! Sum up the timings.
  call mpiallred(time(1), 5, mpi_sum, mpi_comm_world, ierr)
  call mpiallred(timetot, 1, mpi_sum, mpi_comm_world, ierr)
  time=time/dble(nproc)
  timetot=timetot/dble(nproc)
  if(iproc==0) then
      write(*,'(1x,a)') 'timings:'
      write(*,'(3x,a,es10.3)') '-total time:', timetot
      write(*,'(5x,a,es10.3,a,f4.1,a)') '- orthonormalization:', time(1), '=', 100.d0*time(1)/timetot, '%'
      write(*,'(5x,a,es10.3,a,f4.1,a)') '- Hamiltonian application:', time(2),  '=', 100.d0*time(2)/timetot, '%'
      write(*,'(5x,a,es10.3,a,f4.1,a)') '- orthoconstraint:', time(3),  '=', 100.d0*time(3)/timetot, '%'
      write(*,'(5x,a,es10.3,a,f4.1,a)') '- preconditioning:', time(4),  '=', 100.d0*time(4)/timetot, '%'
      write(*,'(5x,a,es10.3,a,f4.1,a)') '- unitary optimization:', time(5),  '=', 100.d0*time(5)/timetot, '%'
      tt=time(1)+time(2)+time(3)+time(4)+time(5)
      tt=timetot-tt
      write(*,'(5x,a,es10.3,a,f4.1,a)') '- other:', tt,  '=', 100.d0*tt/timetot, '%'
  end if


  ! Deallocate all quantities related to DIIS,
  if(ldiis%isx>0) call deallocateDIIS(ldiis)

  ! Deallocate all local arrays.
  call deallocateLocalArrays()

contains



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
          if(iproc==0) write(*,'(1x,a)') 'The force is small enough to allow DIIS.'
          ! This is to get the correct DIIS history 
          ! (it is chosen as max(lin%DIISHistMin,lin%DIISHistMax-icountSwitch).
          icountSwitch=icountSwitch-1
      else if(fnrmMax>lin%startDIIS .and. allowDIIS) then
          allowDIIS=.false.
          if(iproc==0) write(*,'(1x,a)') 'The force is too large to allow DIIS.'
      end if    

      ! Switch to SD if the flag indicating that we should start with SD is true.
      ! If this is the case, this flag is set to false, since this flag concerns only the beginning.
      if(startWithSD .and. ldiis%isx>0) then
          call deallocateDIIS(ldiis)
          ldiis%isx=0
          ldiis%switchSD=.false.
          startWithSD=.false.
      end if

      ! Decide whether we should switch from DIIS to SD in case we are using DIIS and it 
      ! is not allowed.
      if(.not.startWithSD .and. .not.allowDIIS .and. ldiis%isx>0) then
          if(iproc==0) write(*,'(1x,a,es10.3)') 'The force is too large, switch to SD with stepsize', alpha(1)
          call deallocateDIIS(ldiis)
          ldiis%isx=0
          ldiis%switchSD=.true.
      end if

      ! If we swicthed to SD in the previous iteration, reset this flag.
      if(ldiis%switchSD) ldiis%switchSD=.false.

      ! Now come some checks whether the trace is descreasing or not. This further decides
      ! whether we should use DIIS or SD.

      ! Determine wheter the trace is decreasing (as it should) or increasing.
      ! This is done by comparing the current value with diisLIN%energy_min, which is
      ! the minimal value of the trace so far.
      if(trH<=ldiis%trmin .and. .not.resetDIIS) then
          ! Everything ok
          ldiis%trmin=trH
          ldiis%switchSD=.false.
          itBest=it
          icountSDSatur=icountSDSatur+1
          icountDIISFailureCons=0

          ! If we are using SD (i.e. diisLIN%idsx==0) and the trace has been decreasing
          ! for at least 10 iterations, switch to DIIS. However the history length is decreased.
          if(icountSDSatur>=10 .and. ldiis%isx==0 .and. allowDIIS .or. immediateSwitchToSD) then
              icountSwitch=icountSwitch+1
              idsx=max(lin%DIISHistMin,lin%DIISHistMax-icountSwitch)
              if(idsx>0) then
                  if(iproc==0) write(*,'(1x,a,i0)') 'switch to DIIS with new history length ', idsx
                  call initializeDIIS(lin%DIISHistMax, lin%lzd, lin%orbs, lin%orbs%inWhichLocregp, lin%startWithSD, lin%alphaSD, &
                       lin%alphaDIIS, lin%orbs%norb, icountSDSatur, icountSwitch, icountDIISFailureTot, &
                       icountDIISFailureCons, allowDIIS, startWithSD, ldiis, alpha, alphaDIIS)
                  icountDIISFailureTot=0
                  icountDIISFailureCons=0
                  immediateSwitchToSD=.false.
              end if
          end if
      else
          ! The trace is growing.
          ! Count how many times this occurs and (if we are using DIIS) switch to SD after 3 
          ! total failures or after 2 consecutive failures.
          icountDIISFailureCons=icountDIISFailureCons+1
          icountDIISFailureTot=icountDIISFailureTot+1
          icountSDSatur=0
          if((icountDIISFailureCons>=2 .or. icountDIISFailureTot>=3 .or. resetDIIS) .and. ldiis%isx>0) then
          !if((icountDIISFailureCons>=400 .or. icountDIISFailureTot>=600 .or. resetDIIS) .and. ldiis%isx>0) then
              ! Switch back to SD. The initial step size is 1.d0.
              alpha=lin%alphaSD
              if(iproc==0) then
                  if(icountDIISFailureCons>=2) write(*,'(1x,a,i0,a,es10.3)') 'DIIS failed ', &
                      icountDIISFailureCons, ' times consecutively. Switch to SD with stepsize', alpha(1)
                  if(icountDIISFailureTot>=3) write(*,'(1x,a,i0,a,es10.3)') 'DIIS failed ', &
                      icountDIISFailureTot, ' times in total. Switch to SD with stepsize', alpha(1)
                  if(resetDIIS) write(*,'(1x,a)') 'reset DIIS due to flag'
              end if
              if(resetDIIS) then
                  resetDIIS=.false.
                  immediateSwitchToSD=.true.
                  ldiis%trmin=1.d100
              end if
              ! Try to get back the orbitals of the best iteration. This is possible if
              ! these orbitals are still present in the DIIS history.
              if(it-itBest<ldiis%isx) then
                 if(iproc==0) then
                     if(iproc==0) write(*,'(1x,a,i0,a)')  'Recover the orbitals from iteration ', &
                         itBest, ' which are the best so far.'
                 end if
                 ii=modulo(ldiis%mis-(it-itBest),ldiis%mis)
                 offset=0
                 istdest=1
                 do iorb=1,lin%orbs%norbp
                     ilr=lin%orbs%inWhichLocregp(iorb)
                     ncount=lin%lzd%llr(ilr)%wfd%nvctr_c+7*lin%lzd%llr(ilr)%wfd%nvctr_f
                     istsource=offset+ii*ncount+1
                     !write(*,'(a,4i9)') 'iproc, ncount, istsource, istdest', iproc, ncount, istsource, istdest
                     call dcopy(ncount, ldiis%phiHist(istsource), 1, lphi(istdest), 1)
                     call dcopy(ncount, ldiis%phiHist(istsource), 1, lphiold(istdest), 1)
                     offset=offset+ldiis%isx*ncount
                     istdest=istdest+ncount
                 end do
             else
                 ! else copy the orbitals of the last iteration to lphiold
                 call dcopy(size(lphi), lphi(1), 1, lphiold(1), 1)
              end if
              call deallocateDIIS(ldiis)
              ldiis%isx=0
              ldiis%switchSD=.true.
          end if
      end if

    end subroutine DIISorSD


    subroutine improveOrbitals()
    !
    ! Purpose:
    ! ========
    !   This subroutine improves the basis functions by following the gradient 
    ! For DIIS 
    !!if (diisLIN%idsx > 0) then
    !!   diisLIN%mids=mod(diisLIN%ids,diisLIN%idsx)+1
    !!   diisLIN%ids=diisLIN%ids+1
    !!end if
    if (ldiis%isx > 0) then
        ldiis%mis=mod(ldiis%is,ldiis%isx)+1
        ldiis%is=ldiis%is+1
    end if

    ! Follow the gradient using steepest descent.
    ! The same, but transposed
    
    ! steepest descent
    if(ldiis%isx==0) then
        istart=1
        do iorb=1,lin%orbs%norbp
            ilr=lin%orbs%inWhichLocregp(iorb)
            ncount=lin%lzd%llr(ilr)%wfd%nvctr_c+7*lin%lzd%llr(ilr)%wfd%nvctr_f
            call daxpy(ncount, -alpha(iorb), lhphi(istart), 1, lphi(istart), 1)
            istart=istart+ncount
        end do
    else
        ! DIIS
        if(lin%alphaDIIS/=0.d0) then
            !call dscal(lin%lorbs%npsidim, lin%alphaDIIS, lphi, 1)
            call dscal(lin%orbs%npsidim, lin%alphaDIIS, lhphi, 1)
        end if
        call optimizeDIIS(iproc, nproc, lin%orbs, lin%orbs, lin%lzd, lin%orbs%inWhichLocregp, lhphi, lphi, ldiis, it)
    end if
    end subroutine improveOrbitals



    subroutine allocateLocalArrays()
    !
    ! Purpose:
    ! ========
    !   This subroutine allocates all local arrays.
    !
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

      !allocate(lhphi(lin%Lorbs%npsidim), stat=istat)
      allocate(lhphi(lin%orbs%npsidim), stat=istat)
      call memocc(istat, lhphi, 'lhphi', subname)
    
      !allocate(lhphiold(lin%Lorbs%npsidim), stat=istat)
      allocate(lhphiold(lin%orbs%npsidim), stat=istat)
      call memocc(istat, lhphiold, 'lhphiold', subname)

    end subroutine allocateLocalArrays


    subroutine deallocateLocalArrays()
    !
    ! Purpose:
    ! ========
    !   This subroutine deallocates all local arrays.
    !
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

      iall=-product(shape(lhphi))*kind(lhphi)
      deallocate(lhphi, stat=istat)
      call memocc(istat, iall, 'lhphi', subname)

      iall=-product(shape(lhphiold))*kind(lhphiold)
      deallocate(lhphiold, stat=istat)
      call memocc(istat, iall, 'lhphiold', subname)

      ! if diisLIN%idsx==0, these arrays have already been deallocated
      !if(diisLIN%idsx>0 .and. lin%DIISHistMax>0) call deallocate_diis_objects(diisLIN,subname)

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



subroutine diagonalizeHamiltonian2(iproc, nproc, orbs, HamSmall, ovrlp, eval)
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
real(8),dimension(orbs%norb, orbs%norb):: HamSmall, ovrlp
real(8),dimension(orbs%norb):: eval

! Local variables
integer:: lwork, info, istat, iall, i, iorb, jorb
real(8),dimension(:),allocatable:: work
character(len=*),parameter:: subname='diagonalizeHamiltonian'

  ! Get the optimal work array size
  lwork=-1 
  allocate(work(1), stat=istat)
  call memocc(istat, work, 'work', subname)
  call dsygv(1, 'v', 'l', orbs%norb, HamSmall(1,1), orbs%norb, ovrlp(1,1), orbs%norb, eval(1), work(1), lwork, info) 
  lwork=work(1) 

  ! Deallocate the work array ane reallocate it with the optimal size
  iall=-product(shape(work))*kind(work)
  deallocate(work, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating work' 
  call memocc(istat, iall, 'work', subname)
  allocate(work(lwork), stat=istat) ; if(istat/=0) stop 'ERROR in allocating work' 
  call memocc(istat, work, 'work', subname)

  ! Diagonalize the Hamiltonian
  call dsygv(1, 'v', 'l', orbs%norb, HamSmall(1,1), orbs%norb, ovrlp(1,1), orbs%norb, eval(1), work(1), lwork, info) 

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


end subroutine diagonalizeHamiltonian2




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
real(8),dimension(:),allocatable:: phi2, hphi2



  allocate(phi2(orbs%npsidim), stat=istat)
  call memocc(istat, phi2, 'phi2', subname)
  allocate(hphi2(orbs%npsidim), stat=istat)
  call memocc(istat, hphi2, 'hphi2', subname)
  allocate(phiWork(orbs%npsidim), stat=istat)
  call memocc(istat, phiWork, 'phiWork', subname)
  phiWork=0.d0
  call dcopy(orbs%npsidim, phi(1), 1, phi2(1), 1)
  call dcopy(orbs%npsidim, hphi(1), 1, hphi2(1), 1)

  call transpose_v(iproc, nproc, orbs, Glr%wfd, comms, phi2, work=phiWork)
  call transpose_v(iproc, nproc, orbs, Glr%wfd, comms, hphi2, work=phiWork)

  matrixElements=0.d0

  ! Calculate <phi_i|H_j|phi_j>
  nvctrp=sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor
  jstart=1
  do jorb=1,orbs%norb
      istart=1
      do iorb=1,orbs%norb
          matrixElements(iorb,jorb,2)=ddot(nvctrp, phi2(istart), 1, hphi2(jstart), 1)
          !write(*,'(a,3i7,3es16.7)') 'iorb, jorb, iproc, ddot', iorb, jorb, iproc,  matrixElements(iorb,jorb,2), phi2(istart), hphi2(jstart)
          istart=istart+nvctrp
      end do
      jstart=jstart+nvctrp
  end do
  call mpi_allreduce(matrixElements(1,1,2), matrixElements(1,1,1), orbs%norb**2, &
      mpi_double_precision, mpi_sum, mpi_comm_world, ierr)

  call untranspose_v(iproc, nproc, orbs, Glr%wfd, comms, phi2, work=phiWork)
  call untranspose_v(iproc, nproc, orbs, Glr%wfd, comms, hphi2, work=phiWork)



  iall=-product(shape(phiWork))*kind(phiWork)
  deallocate(phiWork, stat=istat)
  call memocc(istat, iall, 'phiWork', subname)

  iall=-product(shape(phi2))*kind(phi2)
  deallocate(phi2, stat=istat)
  call memocc(istat, iall, 'phi2', subname)

  iall=-product(shape(hphi2))*kind(hphi2)
  deallocate(hphi2, stat=istat)
  call memocc(istat, iall, 'hphi2', subname)


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

    if(iproc==0) write(*,'(1x,a)') '============================== optmizing coefficients =============================='

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

        !if(iproc==0) write(*,'(1x,a,4x,i0,es12.4,3x,es10.3, es19.9)') 'iter, fnrm, meanAlpha, Energy', &
        if(iproc==0) write(*,'(1x,a,es11.2,es22.13,es10.2)') 'fnrm, band structure energy, mean alpha', &
            fnrm, ebsMod, meanAlpha
        
        ! Check for convergence.
        if(fnrm<lin%convCritCoeff) then
            if(iproc==0) write(*,'(1x,a,i0,a)') 'converged in ', it, ' iterations.'
            if(iproc==0) write(*,'(3x,a,2es14.5)') 'Final values for fnrm, Energy:', fnrm, ebsMod
            converged=.true.
            infoCoeff=it
            exit
        end if
  
        if(it==lin%nItCoeff) then
            if(iproc==0) write(*,'(1x,a,i0,a)') 'WARNING: not converged within ', it, &
                ' iterations! Exiting loop due to limitations of iterations.'
            if(iproc==0) write(*,'(1x,a,2es15.7,f12.7)') 'Final values for fnrm, Energy: ', fnrm, ebsMod
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
    !!    if(iproc==0) write(*,'(1x,a,i0,a)') 'WARNING: not converged within ', it, &
    !!        ' iterations! Exiting loop due to limitations of iterations.'
    !!    if(iproc==0) write(*,'(1x,a,2es15.7,f12.7)') 'Final values for fnrm, Energy: ', fnrm, ebsMod
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

    if(iproc==0) write(*,'(1x,a)') '===================================================================================='

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
if(iproc==0) write(*,'(1x,a)', advance='no') 'Posting sends / receives for the calculation of the charge density... '
nreceives=0
nsends=0
lin%comsr%communComplete=.false.
procLoop1: do jproc=0,nproc-1
    orbsLoop1: do iorb=1,lin%comsr%noverlaps(jproc)
        mpisource=lin%comsr%comarr(1,iorb,jproc)
        istsource=lin%comsr%comarr(2,iorb,jproc)
        ncount=lin%comsr%comarr(3,iorb,jproc)
        lrsource=lin%comsr%comarr(4,iorb,jproc)
        mpidest=lin%comsr%comarr(5,iorb,jproc)
        istdest=lin%comsr%comarr(6,iorb,jproc)
        tag=lin%comsr%comarr(7,iorb,jproc)
        if(ncount==0) then
            ! No communication is needed. This should be improved in the initialization, i.e. this communication
            ! with 0 elements should be removed from comgp%noverlaps etc.
            lin%comsr%comarr(8,iorb,jproc)=mpi_request_null
            lin%comsr%comarr(9,iorb,jproc)=mpi_request_null
            lin%comsr%communComplete(iorb,jproc)=.true.
            if(iproc==mpidest) then
                ! This is just to make the check at the end happy.
                nreceives=nreceives+1
            end if
        else
            if(mpisource/=mpidest) then
                ! The orbitals are on different processes, so we need a point to point communication.
                if(iproc==mpisource) then
                    !write(*,'(6(a,i0))') 'sumrho: process ', mpisource, ' sends ', ncount, ' elements from position ', istsource, ' to position ', istdest, ' on process ', mpidest, ', tag=',tag
                    !call mpi_isend(lphi(istsource), ncount, mpi_double_precision, mpidest, tag, mpi_comm_world, lin%comsr%comarr(8,iorb,jproc), ierr)
                    call mpi_isend(sendBuf(istsource), ncount, mpi_double_precision, mpidest, tag, mpi_comm_world,&
                         lin%comsr%comarr(8,iorb,jproc), ierr)
                    lin%comsr%comarr(9,iorb,jproc)=mpi_request_null !is this correct?
                    nsends=nsends+1
                else if(iproc==mpidest) then
                   !write(*,'(6(a,i0))') 'sumrho: process ', mpidest, ' receives ', ncount, &
                   !     ' elements at position ', istdest, ' from position ', istsource, ' on process ', mpisource, ', tag=',tag
                    call mpi_irecv(recvBuf(istdest), ncount, mpi_double_precision, mpisource, tag, mpi_comm_world,&
                         lin%comsr%comarr(9,iorb,jproc), ierr)
                    lin%comsr%comarr(8,iorb,jproc)=mpi_request_null !is this correct?
                    nreceives=nreceives+1
                else
                    lin%comsr%comarr(8,iorb,jproc)=mpi_request_null
                    lin%comsr%comarr(9,iorb,jproc)=mpi_request_null
                end if
            else
                ! The orbitals are on the same process, so simply copy them.
                if(iproc==mpisource) then
                    !write(*,'(6(a,i0))') 'sumrho: process ', iproc, ' copies ', ncount, ' elements from position ', istsource, ' to position ', istdest, ' on process ', iproc, ', tag=',tag
                    call dcopy(ncount, sendBuf(istsource), 1, recvBuf(istdest), 1)
                    lin%comsr%comarr(8,iorb,jproc)=mpi_request_null
                    lin%comsr%comarr(9,iorb,jproc)=mpi_request_null
                    nsends=nsends+1
                    nreceives=nreceives+1
                    lin%comsr%communComplete(iorb,mpisource)=.true.
                else
                    lin%comsr%comarr(8,iorb,jproc)=mpi_request_null
                    lin%comsr%comarr(9,iorb,jproc)=mpi_request_null
                    lin%comsr%communComplete(iorb,mpisource)=.true.
                end if
            end if
        end if
    end do orbsLoop1
end do procLoop1
if(iproc==0) write(*,'(a)') 'done.'

if(nreceives/=lin%comsr%noverlaps(iproc)) then
    write(*,'(1x,a,i0,a,i0,2x,i0)') 'ERROR on process ', iproc, ': nreceives/=lin%comsr%noverlaps(iproc)', nreceives,&
         lin%comsr%noverlaps(iproc)
    stop
end if
call mpi_barrier(mpi_comm_world, ierr)

end subroutine postCommunicationSumrho2





!subroutine initializeCommunicationPotential(iproc, nproc, nscatterarr, lin)
subroutine initializeCommunicationPotential(iproc, nproc, nscatterarr, orbs, lzd, comgp, onWhichAtomAll, tag)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
integer,dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
!type(linearParameters),intent(inout):: lin
type(orbitals_data),intent(in):: orbs
type(local_zone_descriptors),intent(in):: lzd
type(p2pCommsGatherPot),intent(out):: comgp
integer,dimension(orbs%norb),intent(in):: onWhichAtomAll
integer,intent(inout):: tag

! Local variables
integer:: is1, ie1, is2, ie2, is3, ie3, ilr, ii, iorb, iiorb, jproc, kproc, istat, iall
integer:: ioverlap, is3j, ie3j, is3k, ie3k, mpidest, istdest, ioffset, is3min, ie3max
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
    do iorb=1,orbs%norb_par(jproc,0)
        
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
            is3=max(is3j,is3k) ! starting index in z dimension for data to be sent
            ie3=min(ie3j,ie3k) ! ending index in z dimension for data to be sent
            ioffset=is3-is3k ! starting index (in z direction) of data to be sent (actually it is the index -1)
            ioverlap=ioverlap+1
            tag=tag+1
            if(is3<is3min .or. ioverlap==1) then
                is3min=is3
            end if
            if(ie3>ie3max .or. ioverlap==1) then
                ie3max=ie3
            end if
            !write(*,'(a,8i8)') 'jproc, kproc, is3j, ie3j, is3k, ie3k, is3, ie3', jproc, kproc, is3j, ie3j, is3k, ie3k, is3, ie3
            call setCommunicationPotential(kproc, is3, ie3, ioffset, lzd%Glr%d%n1i, lzd%Glr%d%n2i, jproc,&
                 istdest, tag, comgp%comarr(1,ioverlap,jproc))
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

allocate(comgp%communComplete(maxval(comgp%noverlaps),0:nproc-1), stat=istat)
call memocc(istat, comgp%communComplete, 'comgp%communComplete', subname)

iall=-product(shape(iStartEnd))*kind(iStartEnd)
deallocate(iStartEnd, stat=istat)
call memocc(istat, iall, 'iStartEnd', subname)

end subroutine initializeCommunicationPotential


subroutine allocateCommunicationsBuffersPotential(comgp, subname)
use module_base
use module_types
implicit none

! Calling arguments
type(p2pCommsGatherPot),intent(inout):: comgp
character(len=*),intent(in):: subname

! Local variables
integer:: istat

allocate(comgp%recvBuf(comgp%nrecvBuf), stat=istat)
call memocc(istat, comgp%recvBuf, 'comgp%recvBuf', subname)

end subroutine allocateCommunicationsBuffersPotential



subroutine deallocateCommunicationsBuffersPotential(comgp, subname)
use module_base
use module_types
implicit none

! Calling arguments
type(p2pCommsGatherPot),intent(inout):: comgp
character(len=*),intent(in):: subname

! Local variables
integer:: istat, iall

iall=-product(shape(comgp%recvBuf))*kind(comgp%recvBuf)
deallocate(comgp%recvBuf, stat=istat)
call memocc(istat, iall, 'comgp%recvBuf', subname)

end subroutine deallocateCommunicationsBuffersPotential



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
if(iproc==0) write(*,'(1x,a)', advance='no') 'Posting sends / receives for communicating the potential... '
nreceives=0
nsends=0
comgp%communComplete=.false.
destLoop: do jproc=0,nproc-1
    sourceLoop: do kproc=1,comgp%noverlaps(jproc)
        mpisource=comgp%comarr(1,kproc,jproc)
        istsource=comgp%comarr(2,kproc,jproc)
        ncount=comgp%comarr(3,kproc,jproc)
        mpidest=comgp%comarr(4,kproc,jproc)
        istdest=comgp%comarr(5,kproc,jproc)
        tag=comgp%comarr(6,kproc,jproc)
        if(ncount==0) then
            ! No communication is needed. This should be improved in the initialization, i.e. this communication
            ! with 0 elements should be removed from comgp%noverlaps etc.
            comgp%comarr(7,kproc,jproc)=mpi_request_null
            comgp%comarr(8,kproc,jproc)=mpi_request_null
            comgp%communComplete(kproc,jproc)=.true.
            if(iproc==mpidest) then
                ! This is just to make the check at the end happy.
                nreceives=nreceives+1
            end if
        else
            if(mpisource/=mpidest) then
                if(iproc==mpisource) then
                    !write(*,'(6(a,i0))') 'process ', mpisource, ' sends ', ncount, ' elements from position ', istsource, ' to position ', istdest, ' on process ', mpidest, ', tag=',tag
                    call mpi_isend(pot(istsource), ncount, mpi_double_precision, mpidest, tag, mpi_comm_world,&
                         comgp%comarr(7,kproc,jproc), ierr)
                    comgp%comarr(8,kproc,jproc)=mpi_request_null !is this correct?
                    nsends=nsends+1
                else if(iproc==mpidest) then
                   !write(*,'(6(a,i0))') 'process ', mpidest, ' receives ', ncount, &
                   !    ' elements at position ', istdest, ' from position ', istsource, ' on process ', mpisource, ', tag=',tag
                    call mpi_irecv(comgp%recvBuf(istdest), ncount, mpi_double_precision, mpisource, tag, mpi_comm_world,&
                        comgp%comarr(8,kproc,jproc), ierr)
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
                else
                    comgp%comarr(7,kproc,jproc)=mpi_request_null
                    comgp%comarr(8,kproc,jproc)=mpi_request_null
                    comgp%communComplete(kproc,jproc)=.true.
                end if
            end if
        end if
    end do sourceLoop
end do destLoop
if(iproc==0) write(*,'(a)') 'done.'

if(nreceives/=comgp%noverlaps(iproc)) then
    write(*,'(1x,a,i0,a,i0,2x,i0)') 'ERROR on process ', iproc, ': nreceives/=comgp%noverlaps(iproc)',&
         nreceives, comgp%noverlaps(iproc)
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
integer:: kproc, mpisource, mpidest, nfast, nslow, nsameproc, ierr, jproc
integer,dimension(mpi_status_size):: stat
logical:: sendComplete, receiveComplete


! Check whether the communications have completed.
nfast=0
nsameproc=0
testLoop: do 
    do jproc=0,nproc-1
        do kproc=1,comgp%noverlaps(jproc)
           if(comgp%communComplete(kproc,jproc)) cycle
            call mpi_test(comgp%comarr(7,kproc,jproc), sendComplete, stat, ierr)      
            call mpi_test(comgp%comarr(8,kproc,jproc), receiveComplete, stat, ierr)   
            ! Attention: mpi_test is a local function.
            if(sendComplete .and. receiveComplete) comgp%communComplete(kproc,jproc)=.true.
            !!if(comgp%communComplete(kproc,jproc)) then
            !!    !write(*,'(2(a,i0))') 'fast communication; process ', iproc, ' has received orbital ', korb
            !!    mpisource=comgp%comarr(1,kproc,jproc)
            !!    mpidest=comgp%comarr(4,kproc,jproc)
            !!    if(mpisource/=mpidest) then
            !!        nfast=nfast+1
            !!    else
            !!        nsameproc=nsameproc+1
            !!    end if
            !!end if
        end do
    end do
    ! If we made it until here, either all all the communication is
    ! complete or we better wait for each single orbital.
    exit testLoop
end do testLoop

! Since mpi_test is a local function, check whether the communication has completed on all processes.
call mpiallred(comgp%communComplete(1,0), nproc*maxval(comgp%noverlaps), mpi_land, mpi_comm_world, ierr)

! Wait for the communications that have not completed yet
nslow=0
do jproc=0,nproc-1
    do kproc=1,comgp%noverlaps(jproc)
        if(comgp%communComplete(kproc,jproc)) then
            mpisource=comgp%comarr(1,kproc,jproc)
            mpidest=comgp%comarr(4,kproc,jproc)
            if(mpisource==mpidest) then
                nsameproc=nsameproc+1
            else
                nfast=nfast+1
            end if
            cycle
        end if
        !write(*,'(2(a,i0))') 'process ', iproc, ' is waiting for orbital ', korb
        nslow=nslow+1
        call mpi_wait(comgp%comarr(7,kproc,jproc), stat, ierr)   !COMMENTED BY PB
        call mpi_wait(comgp%comarr(8,kproc,jproc), stat, ierr)   !COMMENTED BY PB
        comgp%communComplete(kproc,jproc)=.true.
    end do
end do

call mpiallred(nfast, 1, mpi_sum, mpi_comm_world, ierr)
call mpiallred(nslow, 1, mpi_sum, mpi_comm_world, ierr)
call mpiallred(nsameproc, 1, mpi_sum, mpi_comm_world, ierr)
if(iproc==0) write(*,'(1x,2(a,i0),a)') 'statistics: - ', nfast+nslow, ' point to point communications, of which ', &
                       nfast, ' could be overlapped with computation.'
if(iproc==0) write(*,'(1x,a,i0,a)') '            - ', nsameproc, ' copies on the same processor.'


end subroutine gatherPotential




!!!subroutine diagonalizeHamiltonianParallel(iproc, nproc, norb, ham, ovrlp, eval)
!!!use module_base
!!!use module_types
!!!implicit none
!!!
!!!! Calling arguments
!!!integer,intent(in):: iproc, nproc, norb
!!!real(8),dimension(norb,norb),intent(inout):: ham, ovrlp
!!!real(8),dimension(norb),intent(out):: eval
!!!
!!!! Local variables
!!!integer:: ierr, mbrow, mbcol, i, j, istat, lwork, info, ii1, ii2, nproc_scalapack, iall
!!!integer:: nprocrow, nproccol, context, irow, icol, lnrow, lncol, numroc, jproc, liwork, neval_found, neval_computed
!!!real(8):: tt1, tt2
!!!real(8),dimension(:,:),allocatable:: lmat, loverlap, levec
!!!real(8),dimension(:),allocatable:: work, gap
!!!integer,dimension(9):: desc_levec, desc_lmat, desc_loverlap
!!!integer,dimension(:),allocatable:: iwork, ifail, icluster
!!!character(len=*),parameter:: subname='diagonalizeHamiltonianParallel'
!!!
!!!
!!!
!!!
!!!! Block size for scalapack
!!!mbrow=64
!!!mbcol=64
!!!
!!!! Number of processes that will be involved in the calculation
!!!tt1=dble(norb)/dble(mbrow)
!!!tt2=dble(norb)/dble(mbcol)
!!!ii1=ceiling(tt1)
!!!ii2=ceiling(tt2)
!!!nproc_scalapack = min(ii1*ii2,nproc)
!!!!nproc_scalapack = nproc
!!!if(iproc==0) write(*,'(a,i0,a)') 'scalapack will use ',nproc_scalapack,' processes.'
!!!
!!!! process grid: number of processes per row and column
!!!tt1=sqrt(dble(nproc_scalapack))
!!!ii1=ceiling(tt1)
!!!do i=ii1,nproc_scalapack
!!!    if(mod(nproc_scalapack,i)==0) then
!!!        nprocrow=i
!!!        exit
!!!    end if
!!!end do
!!!nproccol=nproc_scalapack/nprocrow
!!!if(iproc==0) write(*,'(a,i0,a,i0,a)') 'calculation is done on process grid with dimension ',nprocrow,' x ',nproccol,'.'
!!!
!!!
!!!! Initialize blacs context
!!!call blacs_get(-1, 0, context)
!!!call blacs_gridinit(context, 'r', nprocrow, nproccol )
!!!call blacs_gridinfo(context,nprocrow, nproccol, irow, icol)
!!!!write(*,*) 'iproc, irow, icol', iproc, irow, icol
!!!
!!!! Initialize the matrix mat to zero for processes that don't do the calculation.
!!!! For processes participating in the diagonalization, 
!!!! it will be partially (only at the position that process was working on) overwritten with the result. 
!!!! At the end we can the make an allreduce to get the correct result on all processes.
!!!if(irow==-1) ham=0.d0
!!!
!!!! Everything that follows is only done if the current process is part of the grid.
!!!processIf: if(irow/=-1) then
!!!    ! Determine the size of the matrix (lnrow x lncol):
!!!    lnrow = numroc(norb, mbrow, irow, 0, nprocrow)
!!!    lncol = numroc(norb, mbcol, icol, 0, nproccol)
!!!    write(*,'(a,i0,a,i0,a,i0)') 'iproc ',iproc,' will have a local matrix of size ',lnrow,' x ',lncol
!!!
!!!    ! Initialize descriptor arrays.
!!!    call descinit(desc_lmat, norb, norb, mbrow, mbcol, 0, 0, context, lnrow, info)
!!!    call descinit(desc_loverlap, norb, norb, mbrow, mbcol, 0, 0, context, lnrow, info)
!!!    call descinit(desc_levec, norb, norb, mbrow, mbcol, 0, 0, context, lnrow, info)
!!!
!!!    ! Allocate the local array lmat
!!!    allocate(lmat(lnrow,lncol), stat=istat)
!!!    call memocc(istat, lmat, 'lmat', subname)
!!!    allocate(loverlap(lnrow,lncol), stat=istat)
!!!    call memocc(istat, loverlap, 'loverlap', subname)
!!!
!!!    ! Copy the global array mat to the local array lmat.
!!!    ! The same for loverlap and overlap, respectively.
!!!    !call dcopy(norb**2, ham(1,1), 1, mat(1,1), 1)
!!!    !call dcopy(norb**2, ovrlp(1,1), 1, overlap(1,1), 1)
!!!    do i=1,norb
!!!        do j=1,norb
!!!            call pdelset(lmat(1,1), j, i, desc_lmat, ham(j,i))
!!!            call pdelset(loverlap(1,1), j, i, desc_loverlap, ovrlp(j,i))
!!!        end do
!!!    end do
!!!
!!!
!!!    ! Solve the generalized eigenvalue problem.
!!!    allocate(levec(lnrow,lncol), stat=istat)
!!!    call memocc(istat, levec, 'levec', subname)
!!!    allocate(ifail(norb), stat=istat)
!!!    call memocc(istat, ifail, 'ifail', subname)
!!!    allocate(icluster(2*nprocrow*nproccol), stat=istat)
!!!    call memocc(istat, icluster, 'icluster', subname)
!!!    allocate(gap(nprocrow*nproccol), stat=istat)
!!!    call memocc(istat, gap, 'gap', subname)
!!!
!!!    ! workspace query
!!!    lwork=-1
!!!    liwork=-1
!!!    allocate(work(1), stat=istat)
!!!    call memocc(istat, work, 'work', subname)
!!!    allocate(iwork(1), stat=istat) ; if(istat/=0) stop 'ERROR in allocating'
!!!    call memocc(istat, iwork, 'iwork', subname)
!!!    call pdsygvx(1, 'v', 'a', 'l', norb, lmat(1,1), 1, 1, desc_lmat, loverlap(1,1), 1, 1, &
!!!                 desc_loverlap, 0.d0, 1.d0, 0, 1, -1.d0, neval_found, neval_computed, eval(1), &
!!!                 -1.d0, levec(1,1), 1, 1, desc_levec, work, lwork, iwork, liwork, &
!!!                 ifail, icluster, gap, info)
!!!    lwork=ceiling(work(1))
!!!    liwork=iwork(1)
!!!    !write(*,*) 'iproc, lwork, liwork', iproc, lwork, liwork
!!!    iall=-product(shape(work))*kind(work)
!!!    deallocate(work, stat=istat)
!!!    call memocc(istat, iall, 'work', subname)
!!!    iall=-product(shape(iwork))*kind(iwork)
!!!    deallocate(iwork, stat=istat)
!!!    call memocc(istat, iall, 'iwork', subname)
!!!
!!!    allocate(work(lwork), stat=istat)
!!!    call memocc(istat, work, 'work', subname)
!!!    allocate(iwork(liwork), stat=istat)
!!!    call memocc(istat, iwork, 'iwork', subname)
!!!
!!!    call pdsygvx(1, 'v', 'a', 'l', norb, lmat(1,1), 1, 1, desc_lmat, loverlap(1,1), 1, 1, &
!!!                 desc_loverlap, 0.d0, 1.d0, 0, 1, -1.d0, neval_found, neval_computed, eval(1), &
!!!                 -1.d0, levec(1,1), 1, 1, desc_levec, work, lwork, iwork, liwork, &
!!!                 ifail, icluster, gap, info)
!!!
!!!    ! Gather together the eigenvectors from all processes and store them in mat.
!!!    do i=1,norb
!!!        do j=1,norb
!!!            call pdelset2(ham(j,i), levec(1,1), j, i, desc_lmat, 0.d0)
!!!        end do
!!!    end do
!!!
!!!
!!!    iall=-product(shape(lmat))*kind(lmat)
!!!    deallocate(lmat, stat=istat)
!!!    call memocc(istat, iall, 'lmat', subname)
!!!
!!!    iall=-product(shape(levec))*kind(levec)
!!!    deallocate(levec, stat=istat)
!!!    call memocc(istat, iall, 'levec', subname)
!!!
!!!    iall=-product(shape(loverlap))*kind(loverlap)
!!!    deallocate(loverlap, stat=istat)
!!!    call memocc(istat, iall, 'loverlap', subname)
!!!
!!!    iall=-product(shape(work))*kind(work)
!!!    deallocate(work, stat=istat)
!!!    call memocc(istat, iall, 'work', subname)
!!!
!!!    iall=-product(shape(iwork))*kind(iwork)
!!!    deallocate(iwork, stat=istat)
!!!    call memocc(istat, iall, 'iwork', subname)
!!!
!!!    iall=-product(shape(ifail))*kind(ifail)
!!!    deallocate(ifail, stat=istat)
!!!    call memocc(istat, iall, 'ifail', subname)
!!!
!!!    iall=-product(shape(icluster))*kind(icluster)
!!!    deallocate(icluster, stat=istat)
!!!    call memocc(istat, iall, 'icluster', subname)
!!!
!!!    iall=-product(shape(gap))*kind(gap)
!!!    deallocate(gap, stat=istat)
!!!    call memocc(istat, iall, 'gap', subname)
!!!
!!!end if processIF
!!!
!!!! Gather the eigenvectors on all processes.
!!!call mpiallred(ham(1,1), norb**2, mpi_sum, mpi_comm_world, ierr)
!!!
!!!! Broadcast the eigenvalues if required. If nproc_scalapack==nproc, then all processes
!!!! diagonalized the matrix and therefore have the eigenvalues.
!!!if(nproc_scalapack/=nproc) then
!!!    call mpi_bcast(eval(1), norb, mpi_double_precision, 0, mpi_comm_world, ierr)
!!!end if
!!!
!!!
!!!
!!!
!!!end subroutine diagonalizeHamiltonianParallel





subroutine prepare_lnlpspd(iproc, at, input, orbs, rxyz, radii_cf, locregShape, lzd)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => prepare_lnlpspd
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc
  type(atoms_data),intent(in):: at
  type(input_variables),intent(in):: input
  type(orbitals_data),intent(in):: orbs
  real(8),dimension(3,at%nat),intent(in):: rxyz
  real(8),dimension(at%ntypes,3),intent(in):: radii_cf
  character(len=1),intent(in):: locregShape
  type(local_zone_descriptors),intent(inout):: lzd
  
  ! Local variables
  integer:: ilr, istat, iorb
  logical:: calc
  character(len=*),parameter:: subname='prepare_lnlpspd'


  allocate(Lzd%Lnlpspd(Lzd%nlr), stat=istat)
  do ilr=1,Lzd%nlr
      call nullify_nonlocal_psp_descriptors(Lzd%Lnlpspd(ilr))
  end do

  do ilr=1,Lzd%nlr

      nullify(Lzd%Llr(ilr)%projflg) !to avoid problems when deallocating
      calc=.false.
      do iorb=1,orbs%norbp
          if(ilr == orbs%inwhichLocreg(iorb+orbs%isorb)) calc=.true.
      end do
      if (.not. calc) cycle !calculate only for the locreg on this processor, without repeating for same locreg.
      ! allocate projflg
      allocate(Lzd%Llr(ilr)%projflg(at%nat), stat=istat)
      call memocc(istat, Lzd%Llr(ilr)%projflg, 'Lzd%Llr(ilr)%projflg', subname)

      call nlpspd_to_locreg(input, iproc, Lzd%Glr, Lzd%Llr(ilr), rxyz, at, orbs, &
           radii_cf, input%frmult, input%frmult, input%hx, input%hy, input%hz, locregShape, lzd%Gnlpspd, &
           Lzd%Lnlpspd(ilr), Lzd%Llr(ilr)%projflg)

  end do

end subroutine prepare_lnlpspd




subroutine free_lnlpspd(orbs, lzd)
  use module_base
  use module_types
  !use deallocatePointers
  use module_interfaces, exceptThisOne => free_lnlpspd
  implicit none
  
  ! Calling arguments
  type(orbitals_data),intent(in):: orbs
  type(local_zone_descriptors),intent(inout):: lzd

  ! Local variables
  integer:: ilr, iorb, istat, iall
  logical:: go
  character(len=*),parameter:: subname='free_lnlpspd'

  do ilr=1,lzd%nlr

      go=.false.
      do iorb=1,orbs%norbp
         if(ilr == orbs%inwhichLocreg(iorb+orbs%isorb)) go=.true.
      end do
      if (.not. go) cycle !deallocate only for the locreg on this processor, without repeating for same locreg.

      ! Deallocate projflg.
      !call checkAndDeallocatePointer(lzd%llr(ilr)%projflg, 'lzd%llr(ilr)%projflg', subname)
      iall=-product(shape(lzd%llr(ilr)%projflg))*kind(lzd%llr(ilr)%projflg)
      deallocate(lzd%llr(ilr)%projflg, stat=istat)
      call memocc(istat, iall, 'lzd%llr(ilr)%projflg', subname)

      call deallocate_nonlocal_psp_descriptors(lzd%lnlpspd(ilr), subname)
  end do

end subroutine free_lnlpspd




subroutine getCoefficients_new(iproc, nproc, lin, orbs, hamold, lphi, ovrlp, coeff)
use module_base
use module_types
use module_interfaces, exceptThisOne => getCoefficients_new
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(linearParameters),intent(inout):: lin
type(orbitals_data),intent(in):: orbs
real(8),dimension(lin%orbs%norb,lin%orbs%norb),intent(in):: hamold
real(8),dimension(lin%lzd%lpsidimtot),intent(in):: lphi
real(8),dimension(lin%orbs%norb,lin%orbs%norb),intent(inout):: ovrlp
real(8),dimension(lin%orbs%norb,orbs%norb),intent(inout):: coeff

! Local variables
integer:: iorb, jorb, j, k, l, info, istat, iall, jproc, ierr
real(8),dimension(:,:),allocatable:: Q, lambda
real(8):: tt, tt2, tt3, alpha
integer,dimension(:),allocatable:: sendcounts, displs
character(len=*),parameter:: subname='getCoefficients_new'

! this is the "step size"
alpha=1.d-1

allocate(Q(lin%orbs%norb,lin%orbs%norb), stat=istat)
call memocc(istat, Q, 'Q', subname)
allocate(lambda(lin%orbs%norb,orbs%norbp), stat=istat)
call memocc(istat, lambda, 'lambda', subname)
allocate(sendcounts(0:nproc-1), stat=istat)
call memocc(istat, sendcounts, 'sendcounts', subname)
allocate(displs(0:nproc-1), stat=istat)
call memocc(istat, displs, 'displs', subname)


! Calculate the matrices Q=<phi|phiold>
call allocateCommuncationBuffersOrtho(lin%lb%comon, subname)
call getMatrixElements2(iproc, nproc, lin%lzd, lin%lb%orbs, lin%lb%op, lin%lb%comon, lphi, lin%lphiold, lin%mad, Q)
call deallocateCommuncationBuffersOrtho(lin%lb%comon, subname)

! Calculate the right hand sides for all physical orbitals handled by this process.
do iorb=1,orbs%norbp
    do jorb=1,lin%orbs%norb
        tt=0.d0
        ! First part. Check indices of Q.
        do j=1,lin%orbs%norb
            !tt=tt+coeff(j,iorb)*Q(jorb,j)
            tt=tt+coeff(j,iorb)*Q(j,jorb)
        end do
        ! Second part. Keep the value of tt.
        tt2=0.d0
        do j=1,lin%orbs%norb
            do k=1,lin%orbs%norb
                tt2=tt2+coeff(j,iorb)*coeff(k,iorb)*hamold(k,j)
            end do
        end do
        ! Check signs of Q.
        tt3=0.d0
        do l=1,lin%orbs%norb
            !tt3=tt3+Q(jorb,l)
            tt3=tt3+Q(l,jorb)
        end do
        lambda(jorb,iorb)=tt-alpha*(tt-tt2*tt3)
    end do
end do

! Solve the system of linear equations.
! Copy the overlap matrix.
call dposv('l', lin%orbs%norb, orbs%norbp, ovrlp(1,1), lin%orbs%norb, lambda(1,1), lin%orbs%norb, info)

! Communicate the coefficients to all processes.
displs(0)=0
do jproc=0,nproc-1
    sendcounts(jproc)=lin%orbs%norb*orbs%norb_par(jproc,0)
    if(jproc>0) displs(jproc)=displs(jproc-1)+sendcounts(jproc-1)
end do
if (nproc > 1) then
   call mpi_allgatherv(lambda, sendcounts(iproc), mpi_double_precision, coeff, sendcounts, displs, &
        mpi_double_precision, mpi_comm_world, ierr)
else
   call vcopy(sendcounts(iproc),lambda(1,1),1,coeff(1,1),1)
end if

iall=-product(shape(Q))*kind(Q)
deallocate(Q, stat=istat)
call memocc(istat, iall, 'Q', subname)
iall=-product(shape(lambda))*kind(lambda)
deallocate(lambda, stat=istat)
call memocc(istat, iall, 'lambda', subname)
iall=-product(shape(sendcounts))*kind(sendcounts)
deallocate(sendcounts, stat=istat)
call memocc(istat, iall, 'sendcounts', subname)
iall=-product(shape(displs))*kind(displs)
deallocate(displs, stat=istat)
call memocc(istat, iall, 'displs', subname)

end subroutine getCoefficients_new





subroutine apply_orbitaldependent_potential(iproc, nproc, lin, at, input, orbs, lzd, rxyz, psi, vpsi)
use module_base
use module_types
use module_interfaces
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(linearParameters),intent(in):: lin
type(atoms_data),intent(in):: at
type(input_variables),intent(in):: input
type(orbitals_data),intent(in):: orbs
type(local_zone_descriptors),intent(in):: lzd
real(8),dimension(3,at%nat),intent(in):: rxyz
!real(8),dimension(lzd%lpsidimtot),intent(in):: psi  !!!! ATENTION, intent should be in !
real(8),dimension(lzd%lpsidimtot),intent(inout):: psi
real(8),dimension(lzd%lpsidimtot),intent(out):: vpsi

! Local variables
integer:: oidx, iorb, ilr, npot, icenter, i_stat, i_all, ist_c, ist_f, ist, iiorb, iall, ierr
real(8):: hxh, hyh, hzh, ddot, tt, t1, t2, time
type(workarr_sumrho):: work_sr
real(8),dimension(:,:),allocatable:: psir, vpsir
type(workarr_precond) :: work
type(workarrays_quartic_convolutions):: work_conv
character(len=*),parameter:: subname='apply_orbitaldependent_potential'
real(8),dimension(0:3),parameter:: scal=1.d0
real(8),dimension(:,:,:),allocatable:: ypsitemp_c
real(8),dimension(:,:,:,:),allocatable:: ypsitemp_f



  !!psi=1.d0
  vpsi=0.d0
  oidx = 0
  do iorb=1,orbs%norbp
     ilr = orbs%inwhichlocreg(iorb+orbs%isorb)
  
     !initialise the work arrays
     call initialize_work_arrays_sumrho(lzd%llr(ilr), work_sr)

     ! Wavefunction in real space
     allocate(psir(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
     call memocc(i_stat,psir,'psir',subname)
     call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psir)

     !transform the wavefunction in Daubechies basis to the wavefunction in ISF basis
     !the psir wavefunction is given in the spinorial form

     !psi(1+oidx+lzd%llr(ilr)%wfd%nvctr_c:1+oidx+lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f-1)=0.d0

     call daub_to_isf(lzd%llr(ilr), work_sr, psi(1+oidx), psir)
     !apply the potential to the psir wavefunction and calculate potential energy
     hxh=.5d0*input%hx
     hyh=.5d0*input%hy
     hzh=.5d0*input%hz
     !icenter=confinementCenter(iorb)
     icenter=lin%orbs%inWhichLocregp(iorb)
     !components of the potential
     npot=orbs%nspinor
     if (orbs%nspinor == 2) npot=1

     call apply_confinement(iproc, lzd%llr(ilr)%d%n1,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,1,1,1,0,orbs%nspinor, psir, &
          rxyz(1,icenter), hxh, hyh, hzh, lin%potentialprefac(at%iatype(icenter)), lin%confpotorder, &
          lzd%llr(ilr)%nsi1, lzd%llr(ilr)%nsi2, lzd%llr(ilr)%nsi3,  &
          lzd%llr(ilr)%bounds%ibyyzz_r) !optional
     !!call apply_confinement(iproc, lzd%llr(ilr)%d%n1,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,1,1,1,0,orbs%nspinor, psir, &
     !!     rxyz(1,icenter), hxh, hyh, hzh, .01d0*lin%potentialprefac(at%iatype(icenter)), 6, &
     !!     lzd%llr(ilr)%nsi1, lzd%llr(ilr)%nsi2, lzd%llr(ilr)%nsi3,  &
     !!     lzd%llr(ilr)%bounds%ibyyzz_r) !optional


     call isf_to_daub(lzd%llr(ilr), work_sr, psir, vpsi(1+oidx))

     i_all=-product(shape(psir))*kind(psir)
     deallocate(psir,stat=i_stat)
     call memocc(i_stat,i_all,'psir',subname)

     call deallocate_work_arrays_sumrho(work_sr)

     oidx = oidx + (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor

  enddo



  !!vpsi=0.d0
  !!ist_c=1
  !!ist_f=1
  !!do iorb=1,orbs%norbp
  !!   iiorb=iorb+orbs%isorb
  !!   ilr = orbs%inwhichlocreg(iorb+orbs%isorb)
  !!
  !!   icenter=lin%orbs%inWhichLocregp(iorb)
  !!   !components of the potential
  !!   npot=orbs%nspinor
  !!   if (orbs%nspinor == 2) npot=1


  !!   ist_f=ist_f+lzd%llr(ilr)%wfd%nvctr_c

  !!   call allocate_workarrays_quartic_convolutions(lzd%llr(ilr), subname, work_conv)

  !!       call uncompress_for_quartic_convolutions(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
  !!            lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
  !!            lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, &
  !!            lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, &
  !!            lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, &
  !!            lzd%llr(ilr)%wfd%keyg, lzd%llr(ilr)%wfd%keyv,  & 
  !!            lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
  !!            lzd%llr(ilr)%wfd%keyg(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
  !!            lzd%llr(ilr)%wfd%keyv(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)),  & 
  !!            scal, psi(ist_c), psi(ist_f), &
  !!            work_conv)

  !!   if(lin%confpotorder==4) then
  !!    call ConvolQuartic4(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
  !!         lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
  !!         lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, &
  !!         lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, & 
  !!         input%hx, lzd%llr(ilr)%ns1, lzd%llr(ilr)%ns2, lzd%llr(ilr)%ns3, &
  !!         lzd%llr(ilr)%bounds%kb%ibyz_c, lzd%llr(ilr)%bounds%kb%ibxz_c, lzd%llr(ilr)%bounds%kb%ibxy_c, &
  !!         lzd%llr(ilr)%bounds%kb%ibyz_f, lzd%llr(ilr)%bounds%kb%ibxz_f, lzd%llr(ilr)%bounds%kb%ibxy_f, &
  !!         rxyz(1,ilr), lin%potentialprefac(at%iatype(icenter)), .false., 0.d0, &
  !!         work_conv%xx_c, work_conv%xx_f1, work_conv%xx_f, &
  !!         work_conv%xy_c, work_conv%xy_f2, work_conv%xy_f, &
  !!         work_conv%xz_c, work_conv%xz_f4, work_conv%xz_f, &
  !!         work_conv%y_c, work_conv%y_f)
  !!    !!call ConvolQuarticShifted(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
  !!    !!     lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
  !!    !!     lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, &
  !!    !!     lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, & 
  !!    !!     input%hx, lzd%llr(ilr)%ns1, lzd%llr(ilr)%ns2, lzd%llr(ilr)%ns3, &
  !!    !!     lzd%llr(ilr)%bounds%kb%ibyz_c, lzd%llr(ilr)%bounds%kb%ibxz_c, lzd%llr(ilr)%bounds%kb%ibxy_c, &
  !!    !!     lzd%llr(ilr)%bounds%kb%ibyz_f, lzd%llr(ilr)%bounds%kb%ibxz_f, lzd%llr(ilr)%bounds%kb%ibxy_f, &
  !!    !!     rxyz(1,ilr), lin%potentialprefac(at%iatype(icenter)), .false., 0.d0, 2.d0,  &
  !!    !!     work_conv%xx_c, work_conv%xx_f1, work_conv%xx_f, &
  !!    !!     work_conv%xy_c, work_conv%xy_f2, work_conv%xy_f, &
  !!    !!     work_conv%xz_c, work_conv%xz_f4, work_conv%xz_f, &
  !!    !!     work_conv%y_c, work_conv%y_f)
  !!    else if(lin%confpotorder==6) then
  !!    !!allocate(ypsitemp_c(0:lzd%llr(ilr)%d%n1, 0:lzd%llr(ilr)%d%n2, 0:lzd%llr(ilr)%d%n3))
  !!    !!allocate(ypsitemp_f(7,lzd%llr(ilr)%d%nfl1:lzd%llr(ilr)%d%nfu1,lzd%llr(ilr)%d%nfl2:lzd%llr(ilr)%d%nfu2,lzd%llr(ilr)%d%nfl3:lzd%llr(ilr)%d%nfu3))

  !!    call ConvolSextic(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
  !!         lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
  !!         lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, &
  !!         lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, & 
  !!         input%hx, lzd%llr(ilr)%ns1, lzd%llr(ilr)%ns2, lzd%llr(ilr)%ns3, &
  !!         lzd%llr(ilr)%bounds%kb%ibyz_c, lzd%llr(ilr)%bounds%kb%ibxz_c, lzd%llr(ilr)%bounds%kb%ibxy_c, &
  !!         lzd%llr(ilr)%bounds%kb%ibyz_f, lzd%llr(ilr)%bounds%kb%ibxz_f, lzd%llr(ilr)%bounds%kb%ibxy_f, &
  !!         rxyz(1,ilr), lin%potentialprefac(at%iatype(icenter)), .false., 0.d0, &
  !!         work_conv%xx_c, work_conv%xx_f1, work_conv%xx_f, &
  !!         work_conv%xy_c, work_conv%xy_f2, work_conv%xy_f, &
  !!         work_conv%xz_c, work_conv%xz_f4, work_conv%xz_f, &
  !!         work_conv%y_c, work_conv%y_f)

  !!     else
  !!         stop 'wronf conf pot'

  !!     end if
  !!   !!ypsitemp_c=work_conv%y_c
  !!   !!ypsitemp_f=work_conv%y_f


  !!   !! call ConvolQuartic4(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
  !!   !!      lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
  !!   !!      lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, &
  !!   !!      lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, & 
  !!   !!      input%hx, lzd%llr(ilr)%ns1, lzd%llr(ilr)%ns2, lzd%llr(ilr)%ns3, &
  !!   !!      lzd%llr(ilr)%bounds%kb%ibyz_c, lzd%llr(ilr)%bounds%kb%ibxz_c, lzd%llr(ilr)%bounds%kb%ibxy_c, &
  !!   !!      lzd%llr(ilr)%bounds%kb%ibyz_f, lzd%llr(ilr)%bounds%kb%ibxz_f, lzd%llr(ilr)%bounds%kb%ibxy_f, &
  !!   !!      rxyz(1,ilr), lin%potentialprefac(at%iatype(icenter)), .false., 0.d0, &
  !!   !!      work_conv%xx_c, work_conv%xx_f1, work_conv%xx_f, &
  !!   !!      work_conv%xy_c, work_conv%xy_f2, work_conv%xy_f, &
  !!   !!      work_conv%xz_c, work_conv%xz_f4, work_conv%xz_f, &
  !!   !!      work_conv%y_c, work_conv%y_f)

  !!   !!work_conv%y_c=.5d0*work_conv%y_c+.5d0*ypsitemp_c
  !!   !!work_conv%y_f=.5d0*work_conv%y_f+.5d0*ypsitemp_f

  !!    !!deallocate(ypsitemp_c)
  !!    !!deallocate(ypsitemp_f)

  !!   call compress_forstandard(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
  !!        lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
  !!        lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, &
  !!        lzd%llr(Ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, &
  !!        lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, &
  !!        lzd%llr(ilr)%wfd%keyg, lzd%llr(ilr)%wfd%keyv,  & 
  !!        lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
  !!        lzd%llr(ilr)%wfd%keyg(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
  !!        lzd%llr(ilr)%wfd%keyv(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)),  & 
  !!        scal, work_conv%y_c, work_conv%y_f, vpsi(ist_c), vpsi(ist_f))

  !!   call deallocate_workarrays_quartic_convolutions(lzd%llr(ilr), subname, work_conv)

  !!   ist_f = ist_f + 7*lzd%llr(ilr)%wfd%nvctr_f
  !!   ist_c = ist_c + lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f

  !!enddo



end subroutine apply_orbitaldependent_potential




!> Expands the compressed wavefunction in vector form (psi_c,psi_f) into the psig format
subroutine uncompress_for_quartic_convolutions(n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3,  & 
     mseg_c, mvctr_c, keyg_c, keyv_c,  & 
     mseg_f, mvctr_f, keyg_f, keyv_f,  & 
     scal, psi_c, psi_f, &
     work)
  use module_base
  use module_types
  implicit none
  integer,intent(in):: n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3, mseg_c, mvctr_c, mseg_f, mvctr_f
  integer,dimension(mseg_c),intent(in):: keyv_c
  integer,dimension(mseg_f),intent(in):: keyv_f
  integer,dimension(2,mseg_c),intent(in):: keyg_c
  integer,dimension(2,mseg_f),intent(in):: keyg_f
  real(wp),dimension(0:3),intent(in):: scal
  real(wp),dimension(mvctr_c),intent(in):: psi_c
  real(wp),dimension(7,mvctr_f),intent(in):: psi_f
  type(workarrays_quartic_convolutions),intent(out):: work
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

  !!!$omp parallel default(private) &
  !!!$omp shared(scal,psig_c,psig_f,x_f1,x_f2,x_f3) &
  !!!$omp shared(psi_c,psi_f,keyv_c,keyg_c,keyv_f,keyg_f,n1,n2,n3,mseg_c,mseg_f)
  !!! coarse part
  !!!$omp do
  do iseg=1,mseg_c
     jj=keyv_c(iseg)
     j0=keyg_c(1,iseg)
     j1=keyg_c(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        work%xx_c(i,i2,i3)=psi_c(i-i0+jj)*scal(0)
        work%xy_c(i2,i,i3)=psi_c(i-i0+jj)*scal(0)
        work%xz_c(i3,i,i2)=psi_c(i-i0+jj)*scal(0)
     enddo
  enddo
  !!!$omp enddo
  !!! fine part
  !!!$omp do
  do iseg=1,mseg_f
     jj=keyv_f(iseg)
     j0=keyg_f(1,iseg)
     j1=keyg_f(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        work%xx_f1(i,i2,i3)=psi_f(1,i-i0+jj)*scal(1)
        work%xx_f(1,i,i2,i3)=psi_f(1,i-i0+jj)*scal(1)
        work%xy_f(1,i2,i,i3)=psi_f(1,i-i0+jj)*scal(1)
        work%xz_f(1,i3,i,i2)=psi_f(1,i-i0+jj)*scal(1)

        work%xy_f2(i2,i,i3)=psi_f(2,i-i0+jj)*scal(1)
        work%xx_f(2,i,i2,i3)=psi_f(2,i-i0+jj)*scal(1)
        work%xy_f(2,i2,i,i3)=psi_f(2,i-i0+jj)*scal(1)
        work%xz_f(2,i3,i,i2)=psi_f(2,i-i0+jj)*scal(1)

        work%xx_f(3,i,i2,i3)=psi_f(3,i-i0+jj)*scal(2)
        work%xy_f(3,i2,i,i3)=psi_f(3,i-i0+jj)*scal(2)
        work%xz_f(3,i3,i,i2)=psi_f(3,i-i0+jj)*scal(2)

        work%xz_f4(i3,i,i2)=psi_f(4,i-i0+jj)*scal(1)
        work%xx_f(4,i,i2,i3)=psi_f(4,i-i0+jj)*scal(1)
        work%xy_f(4,i2,i,i3)=psi_f(4,i-i0+jj)*scal(1)
        work%xz_f(4,i3,i,i2)=psi_f(4,i-i0+jj)*scal(1)

        work%xx_f(5,i,i2,i3)=psi_f(5,i-i0+jj)*scal(2)
        work%xy_f(5,i2,i,i3)=psi_f(5,i-i0+jj)*scal(2)
        work%xz_f(5,i3,i,i2)=psi_f(5,i-i0+jj)*scal(2)

        work%xx_f(6,i,i2,i3)=psi_f(6,i-i0+jj)*scal(2)
        work%xy_f(6,i2,i,i3)=psi_f(6,i-i0+jj)*scal(2)
        work%xz_f(6,i3,i,i2)=psi_f(6,i-i0+jj)*scal(2)

        work%xx_f(7,i,i2,i3)=psi_f(7,i-i0+jj)*scal(3)
        work%xy_f(7,i2,i,i3)=psi_f(7,i-i0+jj)*scal(3)
        work%xz_f(7,i3,i,i2)=psi_f(7,i-i0+jj)*scal(3)
     enddo
  enddo
 !!!$omp enddo
 !!!$omp end parallel

END SUBROUTINE uncompress_for_quartic_convolutions







function dfactorial(n)
implicit none

! Calling arguments
integer,intent(in):: n
real(8):: dfactorial

! Local variables
integer:: i

  dfactorial=1.d0
  do i=1,n
      dfactorial=dfactorial*dble(i)
  end do

end function dfactorial








subroutine minimize_in_subspace(iproc, nproc, lin, at, input, lpot, GPU, ngatherarr, proj, rxyz, pkernelseq, nlpspd, lphi)
use module_base
use module_types
use module_interfaces, exceptThisOne => minimize_in_subspace
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(linearParameters),intent(inout):: lin
type(atoms_data),intent(in):: at
type(input_variables),intent(in):: input
real(8),dimension(lin%lzd%ndimpotisf),intent(in):: lpot
type(GPU_pointers),intent(inout):: GPU
integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
type(nonlocal_psp_descriptors),intent(in):: nlpspd
real(wp),dimension(nlpspd%nprojel),intent(inout):: proj
real(8),dimension(3,at%nat),intent(in):: rxyz
real(dp), dimension(:), pointer :: pkernelseq
real(8),dimension(lin%orbs%npsidim),intent(inout):: lphi

! Local variables
integer:: ndim_lhchi, jorb, jjorb, jlr, ii, istat, iall, jproc, iiorb, kk, iorb, norbTarget, nprocTemp
integer:: is1, ie1, is2, ie2, is3, ie3, js1, je1, js2, je2, js3, je3, iat, ilr, ierr, tag, jlrold, nlocregPerMPI
integer,dimension(:),allocatable:: onWhichAtomTemp, norb_parTemp, onWhichMPITemp
logical,dimension(:),allocatable:: skip, doNotCalculate
logical:: ovrlpx, ovrlpy, ovrlpz, check_whether_locregs_overlap, withConfinement
real(8),dimension(:),allocatable:: lchi
real(8),dimension(:,:),allocatable:: lhchi
real(8):: ekin_sum, epot_sum, eexctX, eproj_sum, t1, t2, time, tt, tt1, tt2, tt3
real(8),dimension(:,:,:),allocatable:: ham3
character(len=*),parameter:: subname='minimize_in_subspace'




           !resetDIIS=.true.
           ! Apply the Hamiltonian for each atom.
           ! onWhichAtomTemp indicates that all orbitals feel the confining potential
           ! centered on atom iat.
           allocate(onWhichAtomTemp(lin%orbs%norbp), stat=istat)
           call memocc(istat,onWhichAtomTemp,'onWhichAtomTemp',subname)
           allocate(doNotCalculate(lin%lzd%nlr), stat=istat)
           call memocc(istat, doNotCalculate, 'doNotCalculate', subname)
           allocate(skip(lin%lzd%nlr), stat=istat)
           call memocc(istat, skip, 'skip', subname)
           !allocate(skipGlobal(lin%lig%lzdig%nlr,0:nproc-1), stat=istat)
           !call memocc(istat, skipGlobal, 'skipGlobal', subname)


           ! Determine for how many localization regions we need a Hamiltonian application.
           ndim_lhchi=0
           do iat=1,at%nat
               call getIndices(lin%lzd%llr(iat), is1, ie1, is2, ie2, is3, ie3)
               skip(iat)=.true.
               do jorb=1,lin%orbs%norbp
                  jjorb=jorb+lin%orbs%isorb
                   onWhichAtomTemp(jorb)=iat
                   jlr=lin%orbs%inWhichLocreg(jjorb)
                   if(lin%orbs%inWhichlocreg(jorb+lin%orbs%isorb)/=jlr) stop 'this should not happen'
                   call getIndices(lin%lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
                   ovrlpx = ( is1<=je1 .and. ie1>=js1 )
                   ovrlpy = ( is2<=je2 .and. ie2>=js2 )
                   ovrlpz = ( is3<=je3 .and. ie3>=js3 )
                   if(ovrlpx .and. ovrlpy .and. ovrlpz) then
                       if(check_whether_locregs_overlap(lin%lzd%llr(iat), lin%lzd%llr(jlr), lin%lzd%glr)) then
                           skip(iat)=.false.
                       end if
                   end if
               end do
               if(.not.skip(iat)) then
                   ndim_lhchi=ndim_lhchi+1
               end if
           end do


           allocate(lhchi(lin%orbs%npsidim,ndim_lhchi),stat=istat)
           call memocc(istat, lhchi, 'lhchi', subname)
           lhchi=0.d0


           if(iproc==0) write(*,'(1x,a)') 'Hamiltonian application for all atoms. This may take some time.'
           call mpi_barrier(mpi_comm_world, ierr)
           call cpu_time(t1)
           !!call prepare_lnlpspd(iproc, at, input, lin%lig%orbsig, rxyz, radii_cf, lin%locregShape, lin%lig%lzdig)
           !!call full_local_potential2(iproc, nproc, lin%lig%lzdig%glr%d%n1i*lin%lig%lzdig%glr%d%n2i*nscatterarr(iproc,2), &
           !!     lin%lig%lzdig%glr%d%n1i*lin%lig%lzdig%glr%d%n2i*lin%lig%lzdig%glr%d%n3i,&
           !!     lin%lig%lzdig%glr%d%n1i*lin%lig%lzdig%glr%d%n2i*nscatterarr(iproc,1)*input%nspin, input%nspin, lin%lig%orbsig,&
           !!     lin%lig%lzdig, ngatherarr, rhopot, lpot, 2, lin%lig%comgp)

           allocate(lin%lzd%doHamAppl(lin%lzd%nlr), stat=istat)
           call memocc(istat, lin%lzd%doHamAppl, 'lin%lzd%doHamAppl', subname)
           withConfinement=.true.
           ii=0
           allocate(lchi(size(lphi)), stat=istat)
           call memocc(istat, lchi, 'lchi', subname)
           lchi=lphi
           do iat=1,at%nat
               doNotCalculate=.true.
               lin%lzd%doHamAppl=.false.
               !!call mpi_barrier(mpi_comm_world, ierr)
               call getIndices(lin%lzd%llr(iat), is1, ie1, is2, ie2, is3, ie3)
               skip(iat)=.true.
               do jorb=1,lin%orbs%norbp
                   onWhichAtomTemp(jorb)=iat
                   !jlr=onWhichAtomp(jorb)
                   jlr=lin%orbs%inWhichLocregp(jorb)
                   call getIndices(lin%lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
                   ovrlpx = ( is1<=je1 .and. ie1>=js1 )
                   ovrlpy = ( is2<=je2 .and. ie2>=js2 )
                   ovrlpz = ( is3<=je3 .and. ie3>=js3 )
                   if(ovrlpx .and. ovrlpy .and. ovrlpz) then
                       doNotCalculate(jlr)=.false.
                       lin%lzd%doHamAppl(jlr)=.true.
                       skip(iat)=.false.
                   else
                       doNotCalculate(jlr)=.true.
                       lin%lzd%doHamAppl(jlr)=.false.
                   end if
               end do
               !write(*,'(a,2i4,4x,100l4)') 'iat, iproc, doNotCalculate', iat, iproc, doNotCalculate
               if(iproc==0) write(*,'(3x,a,i0,a)', advance='no') 'Hamiltonian application for atom ', iat, '... '
               if(.not.skip(iat)) then
                   ii=ii+1
                   if(lin%nItInguess>0) then
                       call HamiltonianApplication3(iproc, nproc, at, lin%orbs, input%hx, input%hy, input%hz, rxyz, &
                            proj, lin%lzd, ngatherarr, lpot, lchi, lhchi(1,ii), &
                            ekin_sum, epot_sum, eexctX, eproj_sum, input%nspin, GPU, withConfinement, .false., &
                            pkernel=pkernelseq, lin=lin, confinementCenter=onWhichAtomTemp)
                   end if

               else
               end if


               if(iproc==0) write(*,'(a)') 'done.'
           end do


           !!iall=-product(shape(lpot))*kind(lpot)
           !!deallocate(lpot, stat=istat)
           !!call memocc(istat, iall, 'lpot', subname)
            if(ii/=ndim_lhchi) then
               write(*,'(a,i0,a,2(a2,i0))') 'ERROR on process ',iproc,': ii/=ndim_lhchi',ii,ndim_lhchi
               stop
           end if
           call mpi_barrier(mpi_comm_world, ierr)
           call cpu_time(t2)
           time=t2-t1
           if(iproc==0) write(*,'(1x,a,es10.3)') 'time for applying potential:', time



           ! The input guess is possibly performed only with a subset of all processes.
           if(lin%norbsPerProcIG>lin%orbs%norb) then
               norbTarget=lin%orbs%norb
           else
              norbTarget=lin%norbsperProcIG
           end if
           nprocTemp=ceiling(dble(lin%orbs%norb)/dble(norbTarget))
           nprocTemp=min(nprocTemp,nproc)
           if(iproc==0) write(*,'(a,i0,a)') 'The minimization is performed using ', nprocTemp, ' processes.'

           ! Create temporary norb_parTemp, onWhichMPITemp
           allocate(norb_parTemp(0:nprocTemp-1), stat=istat)
           call memocc(istat, norb_parTemp, 'norb_parTemp', subname)
           norb_parTemp=0
           tt=dble(lin%orbs%norb)/dble(nprocTemp)
           ii=floor(tt)
           ! ii is now the number of orbitals that every process has. Distribute the remaining ones.
           norb_parTemp(0:nprocTemp-1)=ii
           kk=lin%orbs%norb-nprocTemp*ii
           norb_parTemp(0:kk-1)=ii+1

           allocate(onWhichMPITemp(lin%orbs%norb), stat=istat)
           call memocc(istat, onWhichMPITemp, 'onWhichMPITemp', subname)
           iiorb=0
           do jproc=0,nprocTemp-1
               do iorb=1,norb_parTemp(jproc)
                   iiorb=iiorb+1
                   onWhichMPITemp(iiorb)=jproc
               end do
           end do

           ! Calculate the number of different matrices that have to be stored on a given MPI process.
           jlrold=0
           nlocregPerMPI=0
           do jorb=1,lin%orbs%norb
               jlr=lin%orbs%inWhichLocreg(jorb)
               !jproc=lin%orbs%onWhichMPI(jorb)
               jproc=onWhichMPITemp(jorb)
               !if(iproc==0) write(*,'(a,5i7)') 'jorb, jlr, jlrold, jproc, nlocregPerMPI', jorb, jlr, jlrold, jproc, nlocregPerMPI
               if(iproc==jproc) then
                   if(jlr/=jlrold) then
                       nlocregPerMPI=nlocregPerMPI+1
                       jlrold=jlr
                   end if
               end if
           end do



           ! Calculate the Hamiltonian matrix.
           call cpu_time(t1)
           allocate(ham3(lin%orbs%norb,lin%orbs%norb,nlocregPerMPI), stat=istat)
           call memocc(istat,ham3,'ham3',subname)
           if(lin%nItInguess>0) then
               if(iproc==0) write(*,*) 'calling getHamiltonianMatrix6'
               call getHamiltonianMatrix6(iproc, nproc, nprocTemp, lin%lzd, lin%orbs, lin%orbs, &
                    onWhichMPITemp, input, lin%orbs%inWhichLocreg, ndim_lhchi, &
                    nlocregPerMPI, lchi, lhchi, skip, lin%mad, &
                    lin%memoryForCommunOverlapIG, lin%locregShape, tag, ham3)
           end if

           iall=-product(shape(lhchi))*kind(lhchi)
           deallocate(lhchi, stat=istat)
           call memocc(istat, iall, 'lhchi',subname)


           ! Build the orbitals phi as linear combinations of the atomic orbitals.
           if(iproc==0) write(*,*) 'calling buildLinearCombinationsLocalized3'
           call buildLinearCombinationsLocalized3(iproc, nproc, lin%orbs, lin%orbs, lin%comms,&
                at, lin%lzd%Glr, input, lin%norbsPerType, &
                lin%orbs%inWhichLocreg, lchi, lphi, rxyz, lin%orbs%inWhichLocreg, &
                lin, lin%lzd, nlocregPerMPI, tag, ham3)

           iall=-product(shape(lchi))*kind(lchi)
           deallocate(lchi, stat=istat)
           call memocc(istat, iall, 'lchi',subname)

           iall=-product(shape(lin%lzd%doHamAppl))*kind(lin%lzd%doHamAppl)
           deallocate(lin%lzd%doHamAppl, stat=istat)
           call memocc(istat, iall, 'lin%lzd%doHamAppl',subname)

           iall=-product(shape(norb_parTemp))*kind(norb_parTemp)
           deallocate(norb_parTemp, stat=istat)
           call memocc(istat, iall, 'norb_parTemp',subname)

           iall=-product(shape(ham3))*kind(ham3)
           deallocate(ham3, stat=istat)
           call memocc(istat, iall, 'ham3',subname)

           ! Deallocate all remaining local arrays.
           iall=-product(shape(onWhichAtomTemp))*kind(onWhichAtomTemp)
           deallocate(onWhichAtomTemp, stat=istat)
           call memocc(istat, iall, 'onWhichAtomTemp',subname)

           iall=-product(shape(doNotCalculate))*kind(doNotCalculate)
           deallocate(doNotCalculate, stat=istat)
           call memocc(istat, iall, 'doNotCalculate',subname)

           iall=-product(shape(skip))*kind(skip)
           deallocate(skip, stat=istat)
           call memocc(istat, iall, 'skip',subname)

           iall=-product(shape(onWhichMPITemp))*kind(onWhichMPITemp)
           deallocate(onWhichMPITemp, stat=istat)
           call memocc(istat, iall, 'onWhichMPITemp',subname)


end subroutine minimize_in_subspace









subroutine unitary_optimization(iproc, nproc, lin, lzd, orbs, at, input, op, comon, rxyz, nit, lphi)
use module_base
use module_types
use module_interfaces, exceptThisOne => unitary_optimization
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, nit
type(linearParameters),intent(inout):: lin
type(local_zone_descriptors),intent(in):: lzd
type(orbitals_data),intent(in):: orbs
type(atoms_data),intent(in):: at
type(input_variables),intent(in):: input
type(overlapParameters),intent(inout):: op
type(p2pCommsOrthonormality),intent(inout):: comon
real(8),dimension(3,at%nat),intent(in):: rxyz
real(8),dimension(orbs%npsidim),intent(inout):: lphi

! Local variables
integer:: it, info, lwork, k, istat, iorb, jorb, iall, ierr, ist, jst, ilrold, ncount, jjorb, iiorb, ilr
real(8):: trace, lstep, dfactorial, energyconf_trial, energyconf_0, energyconf_der0, lstep_optimal, ddot
real(8):: tt1, tt2, tt3, tt4, tt5
real(8):: t1, t2, t1_tot, t2_tot
real(8):: time_convol, time_commun, time_lincomb, time_linalg, time_matrixmodification, time_exponential, time_tot
real(8):: time_matrixelements
complex(8):: ttc
real(8),dimension(:,:),allocatable:: gmat, hamtrans, ttmat, Kmat
complex(8),dimension(:,:),allocatable:: gmatc, omatc
complex(8),dimension(:,:,:),allocatable:: tempmatc
complex(8),dimension(:),allocatable:: work, expD_cmplx
real(8),dimension(:),allocatable:: rwork, eval, lphiovrlp, lvphi
real(8),dimension(:,:,:),allocatable:: tempmat3
character(len=*),parameter:: subname='unitary_optimization'
type(p2pCommsOrthonormality):: comon_local

allocate(gmat(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, gmat, 'gmat', subname)
allocate(hamtrans(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, hamtrans, 'hamtrans', subname)
allocate(gmatc(orbs%norb,orbs%norb), stat=istat)
!call memocc(istat, gmatc, 'gmatc', subname)
allocate(omatc(orbs%norb,orbs%norb), stat=istat)
!call memocc(istat, omatc, 'omatc', subname)
allocate(tempmat3(orbs%norb,orbs%norb,3), stat=istat)
call memocc(istat, tempmat3, 'tempmat3', subname)
allocate(eval(orbs%norb), stat=istat)
call memocc(istat, eval, 'eval', subname)
allocate(expD_cmplx(orbs%norb), stat=istat)
call memocc(istat, expD_cmplx, 'expD_cmplx', subname)
allocate(tempmatc(orbs%norb,orbs%norb,2), stat=istat)
!call memocc(istat, tempmatc, 'tempmatc', subname)
allocate(lphiovrlp(op%ndim_lphiovrlp), stat=istat)
call memocc(istat, lphiovrlp, 'lphiovrlp',subname)
allocate(lvphi(orbs%npsidim), stat=istat)
call memocc(istat, lvphi, 'lvphi', subname)
allocate(Kmat(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, Kmat, 'Kmat', subname)



  call allocateSendBufferOrtho(comon, subname)
  call allocateRecvBufferOrtho(comon, subname)
  ! Extract the overlap region from the orbitals phi and store them in comon%sendBuf.
  if(nit>0) then
      call extractOrbital3(iproc, nproc, orbs, orbs%npsidim, orbs%inWhichLocreg, lzd, op, lphi, comon%nsendBuf, comon%sendBuf)
      call postCommsOverlapNew(iproc, nproc, orbs, op, lzd, lphi, comon, tt1, tt2)
  end if

  energyconf_trial=0.d0 !just to initialize this variable and make the compiler happy
  lstep=0.d0 !just to initialize this variable and make the compiler happy
  lstep_optimal=0.d0 !just to initialize this variable and make the compiler happy

  time_convol=0.d0
  time_lincomb=0.d0
  time_commun=0.d0
  time_linalg=0.d0
  time_exponential=0.d0
  time_matrixmodification=0.d0
  time_matrixElements=0.d0
  t1_tot=mpi_wtime()
  innerLoop: do it=1,nit

  !write(*,*) '1: iproc, associated(comon%recvbuf)', iproc, associated(comon%recvbuf)

      t1=mpi_wtime()
      call apply_orbitaldependent_potential(iproc, nproc, lin, at, input, lin%orbs, lin%lzd, rxyz, lphi, lvphi)
      t2=mpi_wtime()
      time_convol=time_convol+t2-t1
      energyconf_0=ddot(orbs%npsidim, lphi(1), 1, lvphi(1), 1)
      call mpiallred(energyconf_0, 1, mpi_sum, mpi_comm_world, ierr)
      if(iproc==0) write(*,'(a,i6,2es20.10,2es17.7)') 'it, energyconf_0, energyvonf_trial, lstep, lstep_optimal', &
                   it, energyconf_0, energyconf_trial, lstep, lstep_optimal

      t1=mpi_wtime()
      allocate(ttmat(lin%orbs%norb,lin%orbs%norb))
      call collectnew(iproc, nproc, comon, lin%mad,lin%op, lin%orbs, input, lin%lzd, comon%nsendbuf, &
           comon%sendbuf, comon%nrecvbuf, comon%recvbuf, ttmat, tt3, tt4, tt5)
      deallocate(ttmat)
      !write(*,*) '2: iproc, associated(comon%recvbuf)', iproc, associated(comon%recvbuf)
      t2=mpi_wtime()
      time_commun=time_commun+t2-t1

      t1=mpi_wtime()
      !call getMatrixElements2(iproc, nproc, lin%lzd, lin%lb%orbs, lin%lb%op, lin%lb%comon, lphi, lvphi, lin%mad, Kmat)
      !call deallocateRecvBufferOrtho(comon, subname)
      !call deallocateSendBufferOrtho(comon, subname)
      call getMatrixElements2(iproc, nproc, lin%lzd, lin%orbs, lin%op, comon, lphi, lvphi, lin%mad, Kmat)
      !call allocateSendBufferOrtho(comon, subname)
      !call allocateRecvBufferOrtho(comon, subname)
      !write(*,*) '3: iproc, associated(comon%recvbuf)', iproc, associated(comon%recvbuf)
      t2=mpi_wtime()
      time_matrixElements=time_matrixElements+t2-t1


      t1=mpi_wtime()
      ! Construct antisymmtric matrix Gmat
      do iorb=1,orbs%norb
          do jorb=1,orbs%norb
              gmat(jorb,iorb)=2.d0*(Kmat(jorb,iorb)-Kmat(iorb,jorb))
          end do
      end do 
      t2=mpi_wtime()
      time_matrixmodification=time_matrixmodification+t2-t1


      t1=mpi_wtime()
      !Build the complex matrix -iGmat
      do iorb=1,orbs%norb
          do jorb=1,orbs%norb
              gmatc(jorb,iorb)=cmplx(0.d0,-gmat(jorb,iorb),kind=8)
          end do
      end do 
      t2=mpi_wtime()
      time_matrixmodification=time_matrixmodification+t2-t1



      ! Diagonalize Gmatc
      t1=mpi_wtime()
      lwork=10*orbs%norb
      allocate(work(lwork), stat=istat) ! factor of 2 since it is assumed to be complex
      allocate(rwork(lwork), stat=istat)
      call zheev('v', 'l', orbs%norb, gmatc(1,1), orbs%norb, eval(1), work, lwork, rwork, info)
      if(info/=0) stop 'ERROR in zheev'
      deallocate(work)
      deallocate(rwork)
      t2=mpi_wtime()
      time_linalg=time_linalg+t2-t1


      ! Calculate step size
      if(it==1) then
          lstep=5.d-2/(maxval(eval))
      else
          lstep=2.d0*lstep_optimal
      end if

      t1=mpi_wtime()
      ! Calculate exp(-i*l*D) (with D diagonal matrix of eigenvalues).
      ! This is also a diagonal matrix, so only calculate the diagonal part.
      do iorb=1,orbs%norb
         ttc=cmplx(0.d0,-lstep*eval(iorb),kind=8)
         expD_cmplx(iorb)=(0.d0,0.d0)
          do k=0,50
              expD_cmplx(iorb)=expD_cmplx(iorb)+ttc**k/dfactorial(k)
          end do
      end do
      t2=mpi_wtime()
      time_exponential=time_exponential+t2-t1

      t1=mpi_wtime()
      do iorb=1,orbs%norb
          do jorb=1,orbs%norb
              if(iorb==jorb) then
                  tempmatc(jorb,iorb,1)=expD_cmplx(iorb)
              else
                  tempmatc(jorb,iorb,1)=cmplx(0.d0,0.d0,kind=8)
              end if
          end do
      end do
      t2=mpi_wtime()
      time_matrixmodification=time_matrixmodification+t2-t1

      t1=mpi_wtime()
      call zgemm('n', 'c', orbs%norb, orbs%norb, orbs%norb, (1.d0,0.d0), tempmatc(1,1,1), orbs%norb, &
           gmatc(1,1), orbs%norb, (0.d0,0.d0), tempmatc(1,1,2), orbs%norb)
      call zgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, (1.d0,0.d0), gmatc(1,1), orbs%norb, &
           tempmatc(1,1,2), orbs%norb, (0.d0,0.d0), omatc(1,1), orbs%norb)
      t2=mpi_wtime()
      time_linalg=time_linalg+t2-t1

      t1=mpi_wtime()
      ! Build new lphi
      do iorb=1,orbs%norb
          do jorb=1,orbs%norb
              tempmat3(jorb,iorb,1)=real(omatc(jorb,iorb))
          end do
      end do
      t2=mpi_wtime()
      time_matrixmodification=time_matrixmodification+t2-t1

      t1=mpi_wtime()
      !write(*,*) '5: iproc, associated(comon%recvbuf)', iproc, associated(comon%recvbuf)
      call build_new_linear_combinations(lin%lzd, lin%orbs, lin%op, comon%nrecvbuf, &
           comon%recvbuf, tempmat3(1,1,1), .true., lphi)
      t2=mpi_wtime()
      time_lincomb=time_lincomb+t2-t1

      t1=mpi_wtime()
      call apply_orbitaldependent_potential(iproc, nproc, lin, at, input, lin%orbs, lin%lzd, rxyz, lphi, lvphi)
      t2=mpi_wtime()
      time_convol=time_convol+t2-t1

      energyconf_trial=ddot(orbs%npsidim, lphi(1), 1, lvphi(1), 1)
      call mpiallred(energyconf_trial, 1, mpi_sum, mpi_comm_world, ierr)

      ! Calculate the gradient of the confinement
      energyconf_der0=0.d0
      do iorb=1,orbs%norb
          do jorb=1,orbs%norb
              energyconf_der0=energyconf_der0+gmat(jorb,iorb)**2
          end do
      end do
      energyconf_der0=-.5d0*energyconf_der0

      ! Calculate optimal step size
      lstep_optimal = -energyconf_der0*lstep**2/(2.d0*(energyconf_trial-energyconf_0-lstep*energyconf_der0))
      lstep_optimal=min(lstep_optimal,lstep)

      t1=mpi_wtime()
      ! Calculate exp(-i*l*D) (with D diagonal matrix of eigenvalues).
      ! This is also a diagonal matrix, so only calculate the diagonal part.
      do iorb=1,orbs%norb
         ttc=cmplx(0.d0,-lstep_optimal*eval(iorb),kind=8)
         expD_cmplx(iorb)=(0.d0,0.d0)
          do k=0,50
              expD_cmplx(iorb)=expD_cmplx(iorb)+ttc**k/dfactorial(k)
          end do
      end do
      t2=mpi_wtime()
      time_exponential=time_exponential+t2-t1

      t1=mpi_wtime()
      do iorb=1,orbs%norb
          do jorb=1,orbs%norb
              if(iorb==jorb) then
                  tempmatc(jorb,iorb,1)=expD_cmplx(iorb)
              else
                  tempmatc(jorb,iorb,1)=cmplx(0.d0,0.d0,kind=8)
              end if
          end do
      end do
      t2=mpi_wtime()
      time_matrixmodification=time_matrixmodification+t2-t1

      t1=mpi_wtime()
      call zgemm('n', 'c', orbs%norb, orbs%norb, orbs%norb, (1.d0,0.d0), tempmatc(1,1,1), orbs%norb, &
           gmatc(1,1), orbs%norb, (0.d0,0.d0), tempmatc(1,1,2), orbs%norb)
      call zgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, (1.d0,0.d0), gmatc(1,1), orbs%norb, &
           tempmatc(1,1,2), orbs%norb, (0.d0,0.d0), omatc(1,1), orbs%norb)
      t2=mpi_wtime()
      time_linalg=time_linalg+t2-t1


      ! Build new lphi
      do iorb=1,orbs%norb
          do jorb=1,orbs%norb
              tempmat3(jorb,iorb,1)=real(omatc(jorb,iorb))
          end do
      end do
      t1=mpi_wtime()
      call build_new_linear_combinations(lin%lzd, lin%orbs, lin%op, comon%nrecvbuf, &
           comon%recvbuf, tempmat3(1,1,1), .true., lphi)
      t2=mpi_wtime()
      time_lincomb=time_lincomb+t2-t1


      if(it<nit) then
          call extractOrbital3(iproc, nproc, orbs, orbs%npsidim, orbs%inWhichLocreg, lzd, op, lphi, comon%nsendBuf, comon%sendBuf)
          call postCommsOverlapNew(iproc, nproc, orbs, op, lzd, lphi, comon, tt1, tt2)
      end if


  end do innerLoop

  t2_tot=mpi_wtime()
  time_tot=t2_tot-t1_tot
  call mpiallred(time_convol, 1, mpi_max, mpi_comm_world, ierr)
  call mpiallred(time_commun, 1, mpi_max, mpi_comm_world, ierr)
  call mpiallred(time_lincomb, 1, mpi_max, mpi_comm_world, ierr)
  call mpiallred(time_linalg, 1, mpi_max, mpi_comm_world, ierr)
  call mpiallred(time_matrixmodification, 1, mpi_max, mpi_comm_world, ierr)
  call mpiallred(time_exponential, 1, mpi_max, mpi_comm_world, ierr)
  call mpiallred(time_matrixelements, 1, mpi_max, mpi_comm_world, ierr)
  call mpiallred(time_tot, 1, mpi_max, mpi_comm_world, ierr)
  !time_convol=time_convol/dble(nproc)
  !time_commun=time_commun/dble(nproc)
  !time_lincomb=time_lincomb/dble(nproc)
  !time_linalg=time_linalg/dble(nproc)
  !time_matrixmodification=time_matrixmodification/dble(nproc)
  !time_exponential_=time_exponential/dble(nproc)
  !time_tot=time_tot/dble(nproc)
  if(iproc==0) then
      write(*,'(a,es16.6)') 'total time: ',time_tot
      write(*,'(a,es15.7,a,f5.2,a)') 'convolutions: ',time_convol,' (',time_convol/time_tot*100.d0,'%)'
      write(*,'(a,es15.7,a,f5.2,a)') 'linear combinations: ',time_lincomb,' (',time_lincomb/time_tot*100.d0,'%)'
      write(*,'(a,es15.7,a,f5.2,a)') 'communication: ',time_commun,' (',time_commun/time_tot*100.d0,'%)'
      write(*,'(a,es15.7,a,f5.2,a)') 'linear algebra: ',time_linalg,' (',time_linalg/time_tot*100.d0,'%)'
      write(*,'(a,es15.7,a,f5.2,a)') 'matrix modification: ',time_matrixmodification, &
                                     ' (',time_matrixmodification/time_tot*100.d0,'%)'
      write(*,'(a,es15.7,a,f5.2,a)') 'building exponential: ',time_exponential,' (',time_exponential/time_tot*100.d0,'%)'
      write(*,'(a,es15.7,a,f5.2,a)') 'matrix elements ',time_matrixelements,' (',time_matrixelements/time_tot*100.d0,'%)'
  end if


  iall=-product(shape(gmat))*kind(gmat)
  deallocate(gmat, stat=istat)
  call memocc(istat, iall, 'gmat', subname)
  iall=-product(shape(gmatc))*kind(gmatc)
  deallocate(gmatc, stat=istat)
  !call memocc(istat, iall, 'gmatc', subname)
  iall=-product(shape(omatc))*kind(omatc)
  deallocate(omatc, stat=istat)
  !call memocc(istat, iall, 'omatc', subname)
  iall=-product(shape(tempmat3))*kind(tempmat3)
  deallocate(tempmat3, stat=istat)
  call memocc(istat, iall, 'tempmat3', subname)
  iall=-product(shape(eval))*kind(eval)
  deallocate(eval, stat=istat)
  call memocc(istat, iall, 'eval', subname)
  iall=-product(shape(expD_cmplx))*kind(expD_cmplx)
  deallocate(expD_cmplx, stat=istat)
  call memocc(istat, iall, 'expD_cmplx', subname)
  iall=-product(shape(tempmatc))*kind(tempmatc)
  deallocate(tempmatc, stat=istat)
  !call memocc(istat, iall, 'tempmatc', subname)
  iall=-product(shape(hamtrans))*kind(hamtrans)
  deallocate(hamtrans, stat=istat)
  call memocc(istat, iall, 'hamtrans', subname)
  iall=-product(shape(lphiovrlp))*kind(lphiovrlp)
  deallocate(lphiovrlp, stat=istat)
  call memocc(istat, iall, 'lphiovrlp', subname)
  iall=-product(shape(lvphi))*kind(lvphi)
  deallocate(lvphi, stat=istat)
  call memocc(istat, iall, 'lvphi', subname)
  iall=-product(shape(Kmat))*kind(Kmat)
  deallocate(Kmat, stat=istat)
  call memocc(istat, iall, 'Kmat', subname)

  call deallocateRecvBufferOrtho(comon, subname)
  call deallocateSendBufferOrtho(comon, subname)

end subroutine unitary_optimization






subroutine build_new_linear_combinations(lzd, orbs, op, nrecvbuf, recvbuf, omat, reset, lphi)
use module_base
use module_types
implicit none

!Calling arguments
type(local_zone_descriptors),intent(in):: lzd
type(orbitals_data),intent(in):: orbs
type(overlapParameters),intent(in):: op
integer,intent(in):: nrecvbuf
real(8),dimension(nrecvbuf),intent(in):: recvbuf
real(8),dimension(orbs%norb,orbs%norb),intent(in):: omat
logical,intent(in):: reset
real(8),dimension(orbs%npsidim),intent(out):: lphi

! Local variables
integer:: ist, jst, ilrold, iorb, iiorb, ilr, ncount, jorb, jjorb, i, ldim, ind, indout, gdim, iorbref, m, ii
integer:: istart, iend, iseg
real(8):: tt


      ! Build new lphi
      if(reset) then
          lphi=0.d0
      end if

      indout=1
      ilrold=-1
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          ilr=orbs%inwhichlocreg(iiorb)
          if(ilr>ilrold) then
              iorbref=iorb
          end if
          gdim=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
          do jorb=1,op%noverlaps(iiorb)
              jjorb=op%overlaps(jorb,iiorb)
              jst=op%indexInRecvBuf(iorbref,jjorb)
              ldim=op%olr(jorb,iorbref)%wfd%nvctr_c+7*op%olr(jorb,iorbref)%wfd%nvctr_f
              tt=omat(jjorb,iiorb)
              do iseg=1,op%expseg(jorb,iorbref)%nseg
                  istart=op%expseg(jorb,iorbref)%segborders(1,iseg)
                  iend=op%expseg(jorb,iorbref)%segborders(2,iseg)
                  ncount=iend-istart+1
                  call daxpy(ncount, tt, recvBuf(jst), 1, lphi(indout+istart-1), 1)
                  jst=jst+ncount
              end do
          end do
          indout=indout+gdim
          ilrold=ilr

      end do

          

end subroutine build_new_linear_combinations




subroutine flatten_at_edges(iproc, nproc, lin, at, input, orbs, lzd, rxyz, psi)
use module_base
use module_types
use module_interfaces
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(linearParameters),intent(in):: lin
type(atoms_data),intent(in):: at
type(input_variables),intent(in):: input
type(orbitals_data),intent(in):: orbs
type(local_zone_descriptors),intent(in):: lzd
real(8),dimension(3,at%nat),intent(in):: rxyz
real(8),dimension(lzd%lpsidimtot),intent(inout):: psi

! Local variables
integer:: oidx, iorb, ilr, npot, icenter, i_stat, i_all
real(8):: hxh, hyh, hzh, alpha
type(workarr_sumrho):: work_sr
real(8),dimension(:,:),allocatable:: psir
character(len=*),parameter:: subname='flatten_at_edges'



  oidx = 0
  do iorb=1,orbs%norbp
     ilr = orbs%inwhichlocreg(iorb+orbs%isorb)
  
     !initialise the work arrays
     call initialize_work_arrays_sumrho(lzd%llr(ilr), work_sr)

     ! Wavefunction in real space
     allocate(psir(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
     call memocc(i_stat,psir,'psir',subname)
     call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psir)

     !transform the wavefunction in Daubechies basis to the wavefunction in ISF basis
     !the psir wavefunction is given in the spinorial form

     call daub_to_isf(lzd%llr(ilr), work_sr, psi(1+oidx), psir)
     !apply the potential to the psir wavefunction and calculate potential energy
     hxh=.5d0*input%hx
     hyh=.5d0*input%hy
     hzh=.5d0*input%hz
     !icenter=confinementCenter(iorb)
     icenter=lin%orbs%inWhichLocregp(iorb)
     !components of the potential
     npot=orbs%nspinor
     if (orbs%nspinor == 2) npot=1

     alpha=-log(1.d-5)/(1.d0-.7d0*lin%locrad(icenter))**2
     !write(*,*) 'iproc, iorb, alpha', iproc, iorb, alpha
     call flatten(iproc, lzd%llr(ilr)%d%n1,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,1,1,1,0,orbs%nspinor, psir, &
          rxyz(1,icenter), hxh, hyh, hzh, lin%potentialprefac(at%iatype(icenter)), lin%confpotorder, &
          lzd%llr(ilr)%nsi1, lzd%llr(ilr)%nsi2, lzd%llr(ilr)%nsi3, .7d0*lin%locrad(icenter), alpha, &
          lzd%llr(ilr)%bounds%ibyyzz_r) !optional



     call isf_to_daub(lzd%llr(ilr), work_sr, psir, psi(1+oidx))

     i_all=-product(shape(psir))*kind(psir)
     deallocate(psir,stat=i_stat)
     call memocc(i_stat,i_all,'psir',subname)

     call deallocate_work_arrays_sumrho(work_sr)

     oidx = oidx + (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor

  enddo


end subroutine flatten_at_edges


subroutine flatten(iproc, n1, n2, n3, nl1, nl2, nl3, nbuf, nspinor, psir, &
     rxyzConfinement, hxh, hyh, hzh, potentialPrefac, confPotOrder, offsetx, offsety, offsetz, cut, alpha, &
     ibyyzz_r) !optional
use module_base
implicit none
integer, intent(in) :: iproc, n1,n2,n3,nl1,nl2,nl3,nbuf,nspinor, confPotOrder, offsetx, offsety, offsetz
real(wp), dimension(-14*nl1:2*n1+1+15*nl1,-14*nl2:2*n2+1+15*nl2,-14*nl3:2*n3+1+15*nl3,nspinor), intent(inout) :: psir
integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r
real(8),dimension(3),intent(in):: rxyzConfinement
real(8),intent(in):: hxh, hyh, hzh, potentialPrefac, cut, alpha
!local variables
integer :: i1,i2,i3,i1s,i1e,ispinor, order
real(wp) :: tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42,tt
real(wp) :: psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4
real(gp) :: epot_p, epot, values, valuesold


  !the Tail treatment is allowed only in the Free BC case
  if (nbuf /= 0 .and. nl1*nl2*nl3 == 0) stop 'NONSENSE: nbuf/=0 only for Free BC'

  ! The order of the cofinement potential (we take order divided by two, 
  ! since later we calculate (r**2)**order.
  if(confPotOrder==2) then
      ! parabolic potential
      order=1
  else if(confPotOrder==4) then
      ! quartic potential
      order=2
  else if(confPotOrder==6) then
      ! sextic potential
      order=3
  end if
  
  values=0.d0
  valuesold=0.d0

!!!$omp parallel default(private)&
!!!$omp shared(psir,n1,n2,n3,epot,ibyyzz_r,nl1,nl2,nl3,nbuf,nspinor)
  !case without bounds
  i1s=-14*nl1
  i1e=2*n1+1+15*nl1
  !epot_p=0._gp
!!!$omp do
  do i3=-14*nl3,2*n3+1+15*nl3
     if (i3 >= -14+2*nbuf .and. i3 <= 2*n3+16-2*nbuf) then !check for the nbuf case
        do i2=-14*nl2,2*n2+1+15*nl2
           if (i2 >= -14+2*nbuf .and. i2 <= 2*n2+16-2*nbuf) then !check for the nbuf case
              !this if statement is inserted here for avoiding code duplication
              !it is to be seen whether the code results to be too much unoptimised
              if (present(ibyyzz_r)) then
                 !in this case we are surely in Free BC
                 !the min is to avoid to calculate for no bounds
                 do i1=-14+2*nbuf,min(ibyyzz_r(1,i2,i3),ibyyzz_r(2,i2,i3))-14-1
                    psir(i1,i2,i3,:)=0.0_wp
                 enddo
                 i1s=max(ibyyzz_r(1,i2,i3)-14,-14+2*nbuf)
                 i1e=min(ibyyzz_r(2,i2,i3)-14,2*n1+16-2*nbuf)
              end if
              !write(*,'(a,5i8)') 'iproc, i1, i2, i1s, i1e', iproc, i1, i2, i1s, i1e
              
              !here we put the branchments wrt to the spin
              if (nspinor == 4) then
                 stop 'this part is not yet implemented'
              else
                 do ispinor=1,nspinor
                    do i1=i1s,i1e
                       tt=(hxh*dble(i1+offsetx)-rxyzConfinement(1))**2 + (hyh*dble(i2+offsety)-rxyzConfinement(2))**2 + &
                          (hzh*dble(i3+offsetz)-rxyzConfinement(3))**2
                       tt=sqrt(tt)
                       tt=tt-cut

                       if(tt>0.d0) then
                           tt=exp(-alpha*tt**2)
                       else
                           tt=1.d0
                       end if
                       values=values+tt*abs(psir(i1,i2,i3,ispinor))
                       valuesold=valuesold+abs(psir(i1,i2,i3,ispinor))
     
                       tt=tt*psir(i1,i2,i3,ispinor)
                       psir(i1,i2,i3,ispinor)=tt
                    end do
                 end do
              end if
              
              if (present(ibyyzz_r)) then
                 !the max is to avoid the calculation for no bounds
                 do i1=max(ibyyzz_r(1,i2,i3),ibyyzz_r(2,i2,i3))-14+1,2*n1+16-2*nbuf
                    psir(i1,i2,i3,:)=0.0_wp
                 enddo
              end if

           else
              do i1=-14,2*n1+16
                 psir(i1,i2,i3,:)=0.0_wp
              enddo
           endif
        enddo
     else
        do i2=-14,2*n2+16
           do i1=-14,2*n1+16
              psir(i1,i2,i3,:)=0.0_wp
           enddo
        enddo
     endif
  enddo
!!!$omp end do

!!!$omp end parallel

write(*,*) 'values',values/valuesold

END SUBROUTINE flatten
!!***
