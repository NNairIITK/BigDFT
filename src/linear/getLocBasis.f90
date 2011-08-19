subroutine getLinearPsi(iproc, nproc, nspin, Glr, orbs, comms, at, lin, rxyz, rxyzParab, &
    nscatterarr, ngatherarr, rhopot, GPU, input, pkernelseq, phi, psi, psit, updatePhi, &
    infoBasisFunctions, infoCoeff, itSCC, n3p, n3pi, n3d, pkernel, &
    i3s, i3xcsh, ebs, coeff, lphi, radii_cf, nlpspd, proj)
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

! Local variables 
integer:: istat, iall, ind1, ind2, ldim, gdim, ilr, istr, nphibuff, iorb, jorb, istart, korb, jst, nvctrp, ncount, jlr
real(8),dimension(:),allocatable:: hphi, eval, lhphi, lphiold, phiold, lhphiold, hphiold, eps, temparr
real(8),dimension(:,:),allocatable:: HamSmall, ovrlp, ovrlpold, hamold
real(8),dimension(:,:,:),allocatable:: matrixElements
real(8),dimension(:),pointer:: phiWork
real(8):: epot_sum, ekin_sum, eexctX, eproj_sum, trace, tt, ddot, tt2, dnrm2, t1, t2, time
character(len=*),parameter:: subname='getLinearPsi' 
logical:: withConfinement
type(workarr_sumrho):: w
integer:: ist, ierr, iiorb, info

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
      if(lin%useDerivativeBasisFunctions) then
          call dcopy(lin%orbs%npsidim, lin%lphiRestart(1), 1, lphi(1), 1)
      end if

      ! Improve the trace minimizing orbitals.
      call getLocalizedBasis(iproc, nproc, at, orbs, Glr, input, lin, rxyz, nspin, &
          nscatterarr, ngatherarr, rhopot, GPU, pkernelseq, lphi, trace, rxyzParab, &
          itSCC, infoBasisFunctions, radii_cf, ovrlp, nlpspd, proj)
  end if

  ! Calculate the derivative basis functions. Copy the trace minimizing orbitals to lin%lphiRestart.
  if(lin%useDerivativeBasisFunctions) then
      call dcopy(lin%orbs%npsidim, lphi(1), 1, lin%lphiRestart(1), 1)
      if(iproc==0) write(*,'(x,a)',advance='no') 'calculating derivative basis functions...'
      call getDerivativeBasisFunctions2(iproc, nproc, input%hx, Glr, lin, lin%orbs%npsidim, lin%lphiRestart, lphi)
      if(iproc==0) write(*,'(a)') 'done.'

      ! Normalize the derivative basis functions. Normalize all of them (i.e. also the trace minimizing
      ! orbitals) to keep it easy.
      ! Do not orthogonalize them, since the 'normal' phis are not exactly orthogonal either.
      ist=1
      do iorb=1,lin%lb%orbs%norbp
          ilr=lin%lb%orbs%inWhichLocregp(iorb)
          ncount=lin%lzd%llr(ilr)%wfd%nvctr_c+7*lin%lzd%llr(ilr)%wfd%nvctr_f
          tt=dnrm2(ncount, lphi(ist), 1)
          call dscal(ncount, 1/tt, lphi(ist), 1)
          ist=ist+ncount
      end do
  end if

  ! Get the overlap matrix.
  !if(.not.updatePhi .and. .not.lin%useDerivativeBasisFunctions) then
  if(.not.lin%useDerivativeBasisFunctions) then
      !call getOverlapMatrix(iproc, nproc, lin, input, lphi, lin%mad, ovrlp)
      call getOverlapMatrix2(iproc, nproc, lin%lzd, lin%orbs, lin%comon, lin%op, lphi, lin%mad, ovrlp)
  end if
  if(lin%useDerivativeBasisFunctions) then
      !call getOverlapMatrix2(iproc, nproc, lin%lb%lzd, lin%lb%orbs, lin%lb%comon, lin%lb%op, lphi, ovrlp)
      call getOverlapMatrix2(iproc, nproc, lin%lzd, lin%lb%orbs, lin%lb%comon, lin%lb%op, lphi, lin%mad, ovrlp)
  end if

  !ierr=0
  !do iorb=1,lin%orbs%norb
  !    do jorb=1,lin%orbs%norb
  !        if(ovrlp(jorb,iorb)==0.d0) then
  !            ierr=ierr+1
  !        else
  !            if(iproc==0) write(350,*) iorb,jorb
  !        end if
  !    end do
  !end do
  !if(iproc==0) write(*,*) 'zero comp, total', ierr, lin%orbs%norb**2

  ! Allocate the communication buffers for the calculation of the charge density.
  call allocateCommunicationbufferSumrho(lin%comsr, subname)
  ! Transform all orbitals to real space.
  ist=1
  istr=1
  do iorb=1,lin%lb%orbs%norbp
      ilr=lin%lb%orbs%inWhichLocregp(iorb)
      call initialize_work_arrays_sumrho(lin%lzd%Llr(ilr), w)
      !call daub_to_isf(lin%lzd%llr(ilr), w, lphi(ist), lphir(istr))
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
  

  if(iproc==0) write(*,'(x,a)') '----------------------------------- Determination of the orbitals in this new basis.'

  ! Gather the potential (it has been posted in the subroutine linearScaling) if the basis functions
  ! have not been updated (in that case it was gathered there).
  if(.not.updatePhi) then
      call gatherPotential(iproc, nproc, lin%comgp)
  end if
  ! If we use the derivative basis functions the potential has to be gathered anyway.
  if(lin%useDerivativeBasisFunctions) call gatherPotential(iproc, nproc, lin%lb%comgp)


  ! Apply the Hamitonian to the orbitals. The flag withConfinement=.false. indicates that there is no
  ! confining potential added to the Hamiltonian.
  allocate(lhphi(lin%lb%orbs%npsidim), stat=istat)
  call memocc(istat, lhphi, 'lhphi', subname)
  withConfinement=.false.
  if(iproc==0) write(*,'(x,a)',advance='no') 'Hamiltonian application...'
  if(.not.lin%useDerivativeBasisFunctions) then
      call HamiltonianApplicationConfinement2(input, iproc, nproc, at, lin%lzd, lin%orbs, lin, input%hx, input%hy, input%hz, rxyz,&
           ngatherarr, lin%comgp%nrecvBuf, lin%comgp%recvBuf, lphi, lhphi, &
           ekin_sum, epot_sum, eexctX, eproj_sum, nspin, GPU, radii_cf, lin%comgp, lin%orbs%inWhichLocregp, withConfinement, .true., &
           pkernel=pkernelseq)
  else
      !!call HamiltonianApplicationConfinement2(input, iproc, nproc, at, lin%lb%lzd, lin%lb%orbs, lin, input%hx, input%hy, input%hz, rxyz,&
      !!     ngatherarr, lin%lb%comgp%nrecvBuf, lin%lb%comgp%recvBuf, lphi, lhphi, &
      !!     ekin_sum, epot_sum, eexctX, eproj_sum, nspin, GPU, radii_cf, lin%lb%comgp, lin%orbs%inWhichLocregp, withConfinement, .true., &
      !!     pkernel=pkernelseq)
      call HamiltonianApplicationConfinement2(input, iproc, nproc, at, lin%lzd, lin%lb%orbs, lin, input%hx, input%hy, input%hz, rxyz,&
           ngatherarr, lin%lb%comgp%nrecvBuf, lin%lb%comgp%recvBuf, lphi, lhphi, &
           ekin_sum, epot_sum, eexctX, eproj_sum, nspin, GPU, radii_cf, lin%lb%comgp, lin%lb%orbs%inWhichLocregp, withConfinement, .true., &
           pkernel=pkernelseq)
  end if

  if(iproc==0) write(*,'(x,a)') 'done.'

  ! Deallocate the buffers needed for the communication of the potential.
  call deallocateCommunicationsBuffersPotential(lin%comgp, subname)
  if(lin%useDerivativeBasisFunctions) call deallocateCommunicationsBuffersPotential(lin%lb%comgp, subname)

  ! Calculate the matrix elements <phi|H|phi>.
  !!do iall=1,size(lphi)
  !!    write(24000+iproc,*) iall, lphi(iall)
  !!end do
  !!do iall=1,size(lhphi)
  !!    write(25000+iproc,*) iall, lhphi(iall)
  !!end do
  !call getMatrixElements2(iproc, nproc, lin%lb%lzd, lin%lb%orbs, lin%lb%op, lin%lb%comon, lphi, lhphi, matrixElements)
  call getMatrixElements2(iproc, nproc, lin%lzd, lin%lb%orbs, lin%lb%op, lin%lb%comon, lphi, lhphi, lin%mad, matrixElements)

  !!!if(iproc==0) then
  !!    ierr=0
  !!    do iall=1,lin%lb%orbs%norb
  !!        do istat=1,lin%lb%orbs%norb
  !!            ierr=ierr+1
  !!            write(23000+iproc,*) iall, istat, matrixElements(istat,iall,1)
  !!        end do
  !!    end do
  !!    write(23000+iproc,*) '=============================='
  !!!end if

  

  ! Diagonalize the Hamiltonian, either iteratively or with lapack.
  call mpi_barrier(mpi_comm_world, ierr) !To measure the time correctly.
  call cpu_time(t1)
  if(trim(lin%getCoeff)=='min') then
      call optimizeCoefficients(iproc, orbs, lin, nspin, matrixElements, coeff, infoCoeff)
  else
      ! Make a copy of the matrix elements since dsyev overwrites the matrix and the matrix elements
      ! are still needed later.
      call dcopy(lin%lb%orbs%norb**2, matrixElements(1,1,1), 1, matrixElements(1,1,2), 1)
      !if(trim(lin%diagMethod)=='seq') then
      if(lin%blocksize_pdsyev<0) then
          if(iproc==0) write(*,'(x,a)',advance='no') 'Diagonalizing the Hamiltonian, sequential version... '
          call diagonalizeHamiltonian2(iproc, nproc, lin%lb%orbs, matrixElements(1,1,2), ovrlp, eval)
      !else if(trim(lin%diagMethod)=='par') then
      else
          if(iproc==0) write(*,'(x,a)',advance='no') 'Diagonalizing the Hamiltonian, parallel version... '
          !call diagonalizeHamiltonianParallel(iproc, nproc, lin%lb%orbs%norb, matrixElements(1,1,2), ovrlp, eval)
          call dsygv_parallel(iproc, nproc, lin%blocksize_pdsyev, lin%nproc_pdsyev, mpi_comm_world, 1, 'v', 'l', lin%lb%orbs%norb, &
               matrixElements(1,1,2), lin%lb%orbs%norb, ovrlp, lin%lb%orbs%norb, eval, info)
      end if
      if(iproc==0) write(*,'(a)') 'done.'
      call dcopy(lin%lb%orbs%norb*orbs%norb, matrixElements(1,1,2), 1, coeff(1,1), 1)
      infoCoeff=0

      ! Write some eigenvalues. Don't write all, but only a few around the last occupied orbital.
      if(iproc==0) then
          do iorb=max(orbs%norb-8,1),min(orbs%norb+8,lin%orbs%norb)
              if(iorb==orbs%norb) then
                  write(*,'(x,a,i0,a,es12.5,a)') 'eval(',iorb,')=',eval(iorb),'  <-- last occupied orbital'
              else if(iorb==orbs%norb+1) then
                  write(*,'(x,a,i0,a,es12.5,a)') 'eval(',iorb,')=',eval(iorb),'  <-- first virtual orbital'
              else
                  write(*,'(x,a,i0,a,es12.5)') 'eval(',iorb,')=',eval(iorb)
              end if
          end do
      end if
  end if
  call cpu_time(t2)
  time=t2-t1
  if(iproc==0) write(*,'(x,a,es10.3)') 'time for diagonalizing the Hamiltonian:',time
  !!do iorb=1,lin%lb%orbs%norb
  !!    write(2000+iproc,'(100es9.2)') (coeff(iorb,jorb), jorb=1,orbs%norb)
  !!end do

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
  

  ind1=1
  ind2=1
  phi=0.d0
  do iorb=1,lin%lb%orbs%norbp
      ilr = lin%lb%orbs%inWhichLocregp(iorb)
      ldim=lin%lzd%Llr(ilr)%wfd%nvctr_c+7*lin%lzd%Llr(ilr)%wfd%nvctr_f
      gdim=Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
      call Lpsi_to_global2(iproc, nproc, ldim, gdim, lin%lb%orbs%norb, lin%lb%orbs%nspinor, input%nspin, Glr, lin%lzd%Llr(ilr), lphi(ind2), phi(ind1))
      ind1=ind1+Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
      ind2=ind2+lin%lzd%Llr(ilr)%wfd%nvctr_c+7*lin%lzd%Llr(ilr)%wfd%nvctr_f
  end do
  call transpose_v(iproc, nproc, lin%lb%orbs, Glr%wfd, lin%lb%comms, phi, work=phiWork)

  
  if(iproc==0) then
      write(*,'(x,a)', advance='no') '------------------------------------- Building linear combinations... '
      
  end if
  ! Build the extended orbital psi as a linear combination of localized basis functions phi. for real O(N)
  ! this has to replaced, but at the moment it is still needed.
  call buildWavefunctionModified(iproc, nproc, orbs, lin%lb%orbs, comms, lin%lb%comms, phi, psi, coeff)

  !!!! ATTENTION -- DEBUG
  !!call dcopy(lin%orbs%npsidim, lphi(1), 1, phi(1), 1)
  !!call transpose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, phi, work=phiWork)
  !!nvctrp=sum(lin%comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor
  !!!!call dgemm('n', 'n', nvctrp, lin%orbs%norb, lin%orbs%norb, 1.d0, phi(1), nvctrp, matrixElements(1,1,2), &
  !!!!           lin%orbs%norb, 0.d0, phiWork(1), nvctrp)
  !!phiWork=0.d0
  !!ist=1
  !!do iorb=1,lin%orbs%norb
  !!    ilr=lin%orbs%inWhichLocreg(iorb)
  !!    jst=1
  !!    do jorb=1,lin%orbs%norb
  !!        jlr=lin%orbs%inWhichLocreg(jorb)
  !!        if(jlr==ilr) then
  !!            call daxpy(nvctrp, matrixElements(jorb,iorb,2), phi(jst), 1, phiWork(ist), 1)
  !!        end if
  !!        jst=jst+nvctrp
  !!    end do
  !!    ist=ist+nvctrp
  !!end do
  !!call dcopy(lin%orbs%npsidim, phiWork(1), 1, phi(1), 1)
  !!call untranspose_v(iproc, nproc, lin%orbs, Glr%wfd, lin%comms, phi, work=phiWork)
  !!call dcopy(lin%orbs%npsidim, phi(1), 1, lphi(1), 1)


  call dcopy(orbs%npsidim, psi, 1, psit, 1)
  if(iproc==0) write(*,'(a)') 'done.'
  
  
  if(.not.lin%useDerivativeBasisFunctions) then
      call untranspose_v(iproc, nproc, lin%lb%orbs, Glr%wfd, lin%lb%comms, phi, work=phiWork)
      call untranspose_v(iproc, nproc, orbs, Glr%wfd, comms, psi, work=phiWork)
  else
      ! The first one was commented..?
      call untranspose_v(iproc, nproc, lin%lb%orbs, Glr%wfd, lin%lb%comms, phi, work=phiWork)
      call untranspose_v(iproc, nproc, orbs, Glr%wfd, comms, psi, work=phiWork)
  end if




  ! Copy phi.
  call dcopy(lin%lb%orbs%npsidim, lphi, 1, lin%lphiold, 1)

  ! Copy lhphi.
  call dcopy(lin%lb%orbs%npsidim, lhphi, 1, lin%lhphiold, 1)

  ! Copy the Hamiltonian
  call dcopy(lin%lb%orbs%norb**2, matrixElements(1,1,1), 1, lin%hamold(1,1), 1)


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
real(8) ::epot_sum, ekin_sum, eexctX, eproj_sum, evalmax, eval_zero, t1tot, t2tot, timetot
real(8):: tt, tt2, ddot, fnrm, fnrmMax, meanAlpha, gnrm, gnrm_zero, gnrmMax, t1, t2
integer:: iorb, icountSDSatur, icountSwitch, idsx, icountDIISFailureTot, icountDIISFailureCons, itBest
integer:: istat, istart, ierr, ii, it, iall, nit, ind1, ind2, jorb, i
integer:: ldim, gdim, ilr, ncount, offset, istsource, istdest
real(8),dimension(:),allocatable:: alpha, fnrmOldArr, alphaDIIS, lhphi, lhphiold
real(8),dimension(:,:),allocatable:: HamSmall, fnrmArr, fnrmOvrlpArr
logical:: quiet, allowDIIS, startWithSD, withConfinement
character(len=*),parameter:: subname='getLocalizedBasis'
character(len=1):: message
type(localizedDIISParameters):: ldiis
real(8),dimension(4):: time
real(8),dimension(:),pointer:: potential
real(8),dimension(:),pointer:: phiWork

  !!do iall=1,Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2)
  !!    write(20000+iproc,*) rhopot(iall)
  !!end do
  !!do iall=1,size(lphi)
  !!    write(21000+iproc,*) lphi(iall)
  !!end do



  ! Allocate all local arrays.
  call allocateLocalArrays()
  
  
  if(iproc==0) write(*,'(x,a)') '======================== Creation of the basis functions... ========================'

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

  time=0.d0
  call cpu_time(t1tot)
  iterLoop: do it=1,nit
      fnrmMax=0.d0
      fnrm=0.d0
  
      if (iproc==0) then
          write( *,'(1x,a,i0)') repeat('-',77 - int(log(real(it))/log(10.))) // ' iter=', it
      endif

  
      ! Orthonormalize the orbitals. If the localization regions are smaller that the global box (which
      ! will be the usual case), the orthogonalization can not be done exactly, but only approximately.
      if(iproc==0) then
          write(*,'(x,a)') 'Orthonormalization... '
      end if
      call cpu_time(t1)
      call orthonormalizeLocalized(iproc, nproc, lin%methTransformOverlap, lin%nItOrtho, lin%blocksize_pdsyev, &
           lin%blocksize_pdgemm, lin%orbs, lin%op, lin%comon, lin%lzd, lin%orbs%inWhichLocreg, lin%convCritOrtho, &
           input, lin%mad, lphi, ovrlp)
      !!do iorb=1,lin%orbs%norb
      !!    do jorb=1,lin%orbs%norb
      !!        if(iproc==0) write(5000,*) iorb, jorb, ovrlp(jorb,iorb)
      !!    end do
      !!end do
      call cpu_time(t2)
      time(1)=time(1)+t2-t1
  
      ! Calculate the unconstrained gradient by applying the Hamiltonian.
      if(iproc==0) then
          write(*,'(x,a)', advance='no') 'Hamiltonian application... '
      end if
      call cpu_time(t1)
      withConfinement=.true.
      call HamiltonianApplicationConfinement2(input, iproc, nproc, at, lin%lzd, lin%orbs, lin, input%hx, input%hy, input%hz, rxyz,&
           ngatherarr, lin%comgp%nrecvBuf, lin%comgp%recvBuf, lphi, lhphi, &
           ekin_sum, epot_sum, eexctX, eproj_sum, nspin, GPU, radii_cf, lin%comgp, lin%orbs%inWhichLocregp, withConfinement, .true., &
           pkernel=pkernelseq)
      call cpu_time(t2)
      time(2)=time(2)+t2-t1
      if(iproc==0) then
          write(*,'(a)', advance='no') 'done. '
      end if


      !!! OLD VERSION FOR DEBUG !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!call full_local_potential(iproc, nproc, Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2), Glr%d%n1i*Glr%d%n2i*Glr%d%n3i, input%nspin, &
      !!     lin%orbs%norb, lin%orbs%norbp, ngatherarr, rhopot, potential)

      !!call HamiltonianApplication(iproc, nproc, at, lin%orbs, input%hx, input%hy, input%hz, rxyz,&
      !!     nlpspd, proj, Glr, ngatherarr, potential, lphi, lhphi, ekin_sum, epot_sum, eexctX, eproj_sum,&
      !!     input%nspin, GPU, pkernel=pkernelseq)

      !!!deallocate potential
      !!call free_full_potential(nproc, potential, subname)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
      ! Apply the orthoconstraint to the gradient. This subroutine also calculates the trace trH.
      if(iproc==0) then
          write(*,'(a)', advance='no') 'Orthoconstraint... '
      end if
      call cpu_time(t1)
      !!!!!!!!!!!!!!!call orthoconstraintLocalized(iproc, nproc, lin, input, lphi, lhphi, trH)

      call orthoconstraintNonorthogonal(iproc, nproc, lin, input, ovrlp, lphi, lhphi, lin%mad, trH)
      call cpu_time(t2)
      time(3)=time(3)+t2-t1
      if(iproc==0) then
          write(*,'(a)', advance='no') 'done. '
      end if

      !!! OLD VERSION FOR DEBUG !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!allocate(phiWork(size(lphi)))
      !!call transpose_v(iproc, nproc, lin%orbs, glr%wfd, lin%comms, lphi, work=phiWork)
      !!call transpose_v(iproc, nproc, lin%orbs, glr%wfd, lin%comms, lhphi, work=phiWork)
      !!call orthoconstraint(iproc, nproc, lin%orbs, lin%comms, glr%wfd, lphi, lhphi, tt)
      !!call untranspose_v(iproc, nproc, lin%orbs, glr%wfd, lin%comms, lphi, work=phiWork)
      !!call untranspose_v(iproc, nproc, lin%orbs, glr%wfd, lin%comms, lhphi, work=phiWork)
      !!deallocate(phiWork)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
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
              if(tt>.9d0) then
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
  

      ! Precondition the gradient.
      if(iproc==0) then
          write(*,'(a)', advance='no') 'Preconditioning... '
      end if
      gnrm=1.d3 ; gnrm_zero=1.d3
      call cpu_time(t1)

      !!evalmax=lin%orbs%eval(lin%orbs%isorb+1)
      !!do iorb=1,lin%orbs%norbp
      !!  evalmax=max(lin%orbs%eval(lin%orbs%isorb+iorb),evalmax)
      !!enddo
      !!call MPI_ALLREDUCE(evalmax,eval_zero,1,mpidtypd,&
      !!     MPI_MAX,MPI_COMM_WORLD,ierr)

      ind2=1
      do iorb=1,lin%orbs%norbp
          ilr = lin%orbs%inWhichLocregp(iorb)
          call choosePreconditioner2(iproc, nproc, lin%orbs, lin, lin%lzd%Llr(ilr), input%hx, input%hy, input%hz, &
              lin%nItPrecond, lhphi(ind2), at%nat, rxyz, at, it, iorb, eval_zero)
          ind2=ind2+lin%lzd%Llr(ilr)%wfd%nvctr_c+7*lin%lzd%Llr(ilr)%wfd%nvctr_f
      end do
      !call preconditionall(iproc, nproc, lin%orbs, lin%lzd%glr, input%hx, input%hy, input%hz, lin%nItPrecond, lhphi, tt, tt2)
      call cpu_time(t2)
      time(4)=time(4)+t2-t1
      if(iproc==0) then
          write(*,'(a)') 'done. '
      end if

  !!do iall=1,size(lphi)
  !!    write(22000+iproc,*) lphi(iall)
  !!end do

      ! Determine the mean step size for steepest descent iterations.
      tt=sum(alpha)
      meanAlpha=tt/dble(lin%orbs%norb)
  
      ! Write some informations to the screen.
      if(iproc==0) write(*,'(x,a,i6,2es15.7,f17.10)') 'iter, fnrm, fnrmMax, trace', it, fnrm, fnrmMax, trH
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
          !!if(lin%plotBasisFunctions) then
          !!    call plotOrbitals(iproc, lin%orbs, Glr, phi, at%nat, rxyz, lin%onWhichAtom, .5d0*input%hx, &
          !!        .5d0*input%hy, .5d0*input%hz, 1)
          !!end if
          exit iterLoop
      end if
  
  
      ! Determine whether tha basis functions shall be further optimized using DIIS or steepest descent.
      call DIISorSD()
      if(iproc==0) then
          if(ldiis%isx>0) then
              write(*,'(x,3(a,i0))') 'DIIS informations: history length=',ldiis%isx, ', consecutive failures=', &
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


      ! Improve the orbitals, depending on the choice made above.
      call improveOrbitals()


     ! Flush the standard output
      call flush(6) 

  end do iterLoop

  call cpu_time(t2tot)
  timetot=t2tot-t1tot

  ! Sum up the timings.
  call mpiallred(time(1), 4, mpi_sum, mpi_comm_world, ierr)
  call mpiallred(timetot, 1, mpi_sum, mpi_comm_world, ierr)
  time=time/dble(nproc)
  timetot=timetot/dble(nproc)
  if(iproc==0) then
      write(*,'(x,a)') 'timings:'
      write(*,'(3x,a,es10.3)') '-total time:', timetot
      write(*,'(5x,a,es10.3,a,f4.1,a)') '- orthonormalization:', time(1), '=', 100.d0*time(1)/timetot, '%'
      write(*,'(5x,a,es10.3,a,f4.1,a)') '- Hamiltonian application:', time(2),  '=', 100.d0*time(2)/timetot, '%'
      write(*,'(5x,a,es10.3,a,f4.1,a)') '- orthoconstraint:', time(3),  '=', 100.d0*time(3)/timetot, '%'
      write(*,'(5x,a,es10.3,a,f4.1,a)') '- preconditioning:', time(4),  '=', 100.d0*time(4)/timetot, '%'
      tt=time(1)+time(2)+time(3)+time(4)
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
      if(startWithSD .and. ldiis%isx>0) then
          call deallocateDIIS(ldiis)
          ldiis%isx=0
          ldiis%switchSD=.false.
          startWithSD=.false.
      end if

      ! Decide whether we should switch from DIIS to SD in case we are using DIIS and it 
      ! is not allowed.
      if(.not.startWithSD .and. .not.allowDIIS .and. ldiis%isx>0) then
          if(iproc==0) write(*,'(x,a,es10.3)') 'The force is too large, switch to SD with stepsize', alpha(1)
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
      if(trH<=ldiis%trmin) then
          ! Everything ok
          ldiis%trmin=trH
          ldiis%switchSD=.false.
          itBest=it
          icountSDSatur=icountSDSatur+1
          icountDIISFailureCons=0

          ! If we are using SD (i.e. diisLIN%idsx==0) and the trace has been decreasing
          ! for at least 10 iterations, switch to DIIS. However the history length is decreased.
          if(icountSDSatur>=10 .and. ldiis%isx==0 .and. allowDIIS) then
              icountSwitch=icountSwitch+1
              idsx=max(lin%DIISHistMin,lin%DIISHistMax-icountSwitch)
              if(idsx>0) then
                  if(iproc==0) write(*,'(x,a,i0)') 'switch to DIIS with new history length ', idsx
                  call initializeDIIS(lin%DIISHistMax, lin%lzd, lin%orbs, lin%orbs%inWhichLocregp, lin%startWithSD, lin%alphaSD, &
                       lin%alphaDIIS, lin%orbs%norb, icountSDSatur, icountSwitch, icountDIISFailureTot, &
                       icountDIISFailureCons, allowDIIS, startWithSD, ldiis, alpha, alphaDIIS)
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
          if((icountDIISFailureCons>=2 .or. icountDIISFailureTot>=3) .and. ldiis%isx>0) then
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
              if(it-itBest<ldiis%isx) then
                 if(iproc==0) then
                     if(iproc==0) write(*,'(x,a,i0,a)')  'Recover the orbitals from iteration ', &
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
                     offset=offset+ldiis%isx*ncount
                     istdest=istdest+ncount
                 end do
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
                    !write(*,'(6(a,i0))') 'sumrho: process ', mpidest, ' receives ', ncount, ' elements at position ', istdest, ' from position ', istsource, ' on process ', mpisource, ', tag=',tag
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
    write(*,'(x,a,i0,a,i0,2x,i0)') 'ERROR on process ', iproc, ': nreceives/=lin%comsr%noverlaps(iproc)', nreceives,&
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
type(linear_zone_descriptors),intent(in):: lzd
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
    do iorb=1,orbs%norb_par(jproc)
        
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
if(iproc==0) write(*,'(x,a)', advance='no') 'Posting sends / receives for communicating the potential... '
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
                    !write(*,'(6(a,i0))') 'process ', mpidest, ' receives ', ncount, ' elements at position ', istdest, ' from position ', istsource, ' on process ', mpisource, ', tag=',tag
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
    write(*,'(x,a,i0,a,i0,2x,i0)') 'ERROR on process ', iproc, ': nreceives/=comgp%noverlaps(iproc)',&
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
if(iproc==0) write(*,'(x,2(a,i0),a)') 'statistics: - ', nfast+nslow, ' point to point communications, of which ', &
                       nfast, ' could be overlapped with computation.'
if(iproc==0) write(*,'(x,a,i0,a)') '            - ', nsameproc, ' copies on the same processor.'


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
