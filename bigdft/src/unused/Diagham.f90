!>    Diagonalise the hamiltonian in a basis set of norbe orbitals and select the first
!!    norb eigenvectors. Works also with the spin-polarisation case and perform also the 
!!    treatment of semicore atoms. 
!!    In the absence of norbe parameters, it simply diagonalizes the hamiltonian in the given
!!    orbital basis set.
!! @author
!! INPUT VARIABLES
!!    @param iproc  process id
!!    @param nproc  number of mpi processes
!!    @param natsc  number of semicore atoms for the orthogonalisation treatment
!!                  used as a dimension for the array of semicore atoms
!!    @param nspin  spin polarised id; 1 => non spin-polarised; 2 => spin-polarised (collinear)
!!    @param norbu  number of up orbitals in the spin-polarised case; for non spin-pol equal to norb
!!    @param norbd  number of down orbitals in the spin-polarised case; for non spin-pol equal to 0
!!    @param norb   total number of orbitals of the resulting eigenfunctions
!!    @param norbp  number of orbitals in parallel. For nproc=1 norbp=norb
!!    @param nvirte number of virtual orbitals to be saved as input guess 
!!                  for the Davidson method (for both spins)
!!    @param nvctrp number of points of the wavefunctions for each orbital in the transposed sense
!!    @param wfd    data structure of the wavefunction descriptors
!!    @param norbe  (optional) number of orbitals of the initial set of wavefunction, to be reduced
!!    @param etol    tolerance for which a degeneracy should be printed. Set to zero if absent
!! INPUT-OUTPUT VARIABLES
!!    @param psi    wavefunctions. 
!!                  - If norbe is absent: on input, set of norb wavefunctions, 
!!                                      on output eigenfunctions
!!                  - If norbe is present: on input, set of norbe wavefunctions, 
!!                                      on output the first norb eigenfunctions
!!    @param hpsi   hamiltonian on the wavefunctions
!!                  - If norbe is absent: on input, set of norb arrays, 
!!                                      destroyed on output
!!                  - If norbe is present: on input, set of norbe wavefunctions, 
!!                                      destroyed on output
!! OUTPUT VARIABLES
!!    @param psit   wavefunctions in the transposed form.
!!           On input: nullified
!!           on Output: transposed wavefunction but only if nproc>1, nullified otherwise
!!    @param psivirt wavefunctions for input guess of the Davidson method in gaussian form
!!           On input, if present: coefficients of the orbitals in the gaussian basis set 
!!           if nvirte >0: on Output, eigenvectors after input guess
!!           if nvirte=0: unchanged on output
!!    @param eval   array of the first norb eigenvalues       
!! Author:
!!    Luigi Genovese (2007-2009-2011)
!!    Stephan Mohr (2010-2011)
subroutine DiagHam(iproc,nproc,natsc,nspin,orbs,wfd,comms,&
      &   psi,hpsi,psit,orthpar,passmat,& !mandatory
      orbse,commse,etol,norbsc_arr,orbsv,psivirt) !optional
   use module_base
   use module_types
   use yaml_output
   use module_interfaces, except_this_one => DiagHam
   implicit none
   integer, intent(in) :: iproc,nproc,natsc,nspin
   type(wavefunctions_descriptors), intent(in) :: wfd
   type(comms_cubic), target, intent(in) :: comms
   type(orbitals_data), target, intent(inout) :: orbs
   type(orthon_data), intent(inout) :: orthpar
   real(wp), dimension(*), intent(out) :: passmat !< passage matrix for building the eigenvectors (the size depends of the optional arguments)
   real(wp), dimension(:), pointer :: psi,hpsi,psit
   !optional arguments
   real(gp), optional, intent(in) :: etol
   type(orbitals_data), optional, intent(in) :: orbsv
   type(orbitals_data), optional, target, intent(in) :: orbse
   type(comms_cubic), optional, target, intent(in) :: commse
   integer, optional, dimension(natsc+1,nspin), intent(in) :: norbsc_arr
   real(wp), dimension(:), pointer, optional :: psivirt
   !local variables
   character(len=*), parameter :: subname='DiagHam'
   !n(c) real(kind=8), parameter :: eps_mach=1.d-12
  logical :: semicore,minimal
   integer :: ikptp,ikpt,nvctrp
   integer :: i,ndim_hamovr,i_all,i_stat,ierr,norbi_max,j,noncoll,ispm,ncplx
   integer :: norbtot,natsceff,norbsc,ndh1,ispin,npsidim,nspinor,ispsi,ispsie,ispsiv !n(c) nvctr
   real(gp) :: tolerance
   type(orbitals_data), pointer :: orbsu
   type(comms_cubic), pointer :: commu
   integer, dimension(:,:), allocatable :: norbgrp
   real(wp), dimension(:,:,:), allocatable :: hamovr
   real(wp), dimension(:), pointer :: psiw

   !performs some check of the arguments
   if (present(orbse) .neqv. present(commse)) then
      !if (iproc ==0) 
      call yaml_warning('(DiagHam) The variables orbse and commse must be present at the same time')
      !write(*,'(1x,a)') 'ERROR (DiagHam): the variables orbse and commse must be present at the same time'
      stop
   else
      minimal=present(orbse)
   end if

   if (present(etol)) then
      tolerance=etol
   else
      tolerance=0.0_gp
   end if

   semicore=present(norbsc_arr)

   !assign total orbital number for calculating the overlap matrix and diagonalise the system

   if(minimal) then
      norbtot=orbse%norb !beware that norbe is equal both for spin up and down
      commu => commse
      orbsu => orbse
      npsidim=max(orbse%npsidim_orbs,orbse%npsidim_comp)
      nspinor=orbse%nspinor
   else
      norbtot=orbs%norb
      commu => comms
      orbsu => orbs
      npsidim=max(orbs%npsidim_orbs,orbs%npsidim_comp)
      nspinor=orbs%nspinor
   end if
   if (nproc > 1) then
      allocate(psiw(npsidim+ndebug),stat=i_stat)
      call memocc(i_stat,psiw,'psiw',subname)
   else
      nullify(psiw)
   end if

   !transpose all the wavefunctions for having a piece of all the orbitals 
   !for each processor
   call transpose_v(iproc,nproc,orbsu,wfd,commu,psi,work=psiw)
   call transpose_v(iproc,nproc,orbsu,wfd,commu,hpsi,work=psiw)

   if (nproc > 1) then
      i_all=-product(shape(psiw))*kind(psiw)
      deallocate(psiw,stat=i_stat)
      call memocc(i_stat,i_all,'psiw',subname)
   end if

   !define the grouping of the orbitals: for the semicore case, follow the semicore atoms,
   !otherwise use the number of orbitals, separated in the spin-polarised case
   !for the spin-polarised case it is supposed that the semicore orbitals are disposed equally
   !fon non-collinear spins, the orbitals atoms are doubled
   !calculate the maximum of the dimension for each k-point
   if (orbs%nspinor == 4) then
      noncoll=2
   else
      noncoll=1
   end if

   if (semicore) then
      if (present(orbsv)) then
         norbi_max=max(noncoll*maxval(norbsc_arr),orbsv%norb)
      else
         norbi_max=noncoll*maxval(norbsc_arr)
      end if

      !calculate the dimension of the overlap matrix
      !take the maximum as the two spin dimensions
      ndim_hamovr=0
      do ispin=1,nspin
         ndh1=0
         norbsc=0
         do i=1,natsc+1
            ndh1=ndh1+(noncoll*norbsc_arr(i,ispin))**2
         end do
         ndim_hamovr=max(ndim_hamovr,ndh1)
      end do
      if (natsc > 0) then
         if (nspin == 2) then
            if (sum(norbsc_arr(1:natsc,1)) /= sum(norbsc_arr(1:natsc,2))) then
               call yaml_warning('(DiagHam) The number of semicore orbitals must be the same for both spins')
               !write(*,'(1x,a)') 'ERROR (DiagHam): The number of semicore orbitals must be the same for both spins'
               stop
            end if
         end if
         norbsc=noncoll*sum(norbsc_arr(1:natsc,1))
      else
         norbsc=0
      end if

      natsceff=natsc
      allocate(norbgrp(natsceff+1,nspin+ndebug),stat=i_stat)
      call memocc(i_stat,norbgrp,'norbgrp',subname)

      !assign the grouping of the orbitals
      do j=1,nspin
         do i=1,natsceff+1
            norbgrp(i,j)=noncoll*norbsc_arr(i,j)
         end do
      end do
   else
      !this works also for non spin-polarised since there norbu=norb
      norbi_max=max(orbs%norbu,orbs%norbd) 
      ndim_hamovr=norbi_max**2
      natsceff=0
      allocate(norbgrp(1,nspin+ndebug),stat=i_stat)
      call memocc(i_stat,norbgrp,'norbgrp',subname)

      norbsc=0
      norbgrp(1,1)=orbs%norbu
      if (nspin == 2) norbgrp(1,2)=orbs%norbd

   end if

   !for complex matrices the dimension is doubled
   if (nspinor /=1) then
      ndim_hamovr=2*ndim_hamovr
   end if

   allocate(hamovr(nspin*ndim_hamovr,2,orbsu%nkpts+ndebug),stat=i_stat)
   call memocc(i_stat,hamovr,'hamovr',subname)

   !initialise hamovr
   call to_zero(nspin*ndim_hamovr*2*orbsu%nkpts,hamovr)

   if (iproc == 0 .and. verbose > 1) call yaml_comment('Overlap Matrix...')
   !if (iproc == 0 .and. verbose > 1) write(*,'(1x,a)',advance='no') 'Overlap Matrix...'

   !after having applied the hamiltonian to all the atomic orbitals
   !we split the semicore orbitals from the valence ones
   !this is possible since the semicore orbitals are the first in the 
   !order, so the linear algebra on the transposed wavefunctions 
   !may be splitted
   ispsi=1
   do ikptp=1,orbsu%nkptsp
      ikpt=orbsu%iskpts+ikptp!orbsu%ikptsp(ikptp)

      nvctrp=commu%nvctr_par(iproc,ikpt )
      if (nvctrp == 0) cycle

      !print *,'iproc,nvctrp,nspin,norb,ispsi,ndimovrlp',iproc,nvctrp,nspin,norb,ispsi,ndimovrlp(ispin,ikpt-1)
      call overlap_matrices(norbtot,nvctrp,natsceff,nspin,nspinor,&
         &   ndim_hamovr,norbgrp,hamovr(1,1,ikpt),psi(ispsi),hpsi(ispsi))

      ispsi=ispsi+nvctrp*norbtot*orbsu%nspinor
   end do

   !  if(iproc==0 .and. verbose>1) write(*,'(a)') ' done.'
   !if (iproc == 0) print *,'hamovr,iproc:',iproc,hamovr

   if (minimal) then
      !deallocate hpsi in the case of a minimal basis
      i_all=-product(shape(hpsi))*kind(hpsi)
      deallocate(hpsi,stat=i_stat)
      call memocc(i_stat,i_all,'hpsi',subname)
   end if

   if (nproc > 1) then
      !reduce the overlap matrix between all the processors
      call mpiallred(hamovr(1,1,1),2*nspin*ndim_hamovr*orbsu%nkpts,&
         &   MPI_SUM,bigdft_mpi%mpi_comm,ierr)
   end if

! DEBUG
!  if(iproc == 0) then
!     print *,'size(hamovr)',size(hamovr,1),size(hamovr,2),size(hamovr,3)
!     do i_all=1,size(hamovr,1)
!        print *,'iel, ham, ovr:',i_all,hamovr(i_all,1,:),hamovr(i_all,2,:)
!     end do
!  end if
! END DEBUG

   !in the case of minimal basis allocate now the transposed wavefunction
   !otherwise do it only in parallel
   if (.not. associated(psit)) then
      if (minimal .or. nproc > 1) then
         allocate(psit(max(orbs%npsidim_orbs,orbs%npsidim_comp)+ndebug),stat=i_stat)
         call memocc(i_stat,psit,'psit',subname)
      else
         psit => hpsi
      end if
   end if

   ! There are two possibilities to generate the input guess
   differentInputGuess: if(.not. orthpar%directDiag) then
      if(iproc==0) call yaml_comment('Iterative diagonalization...')
      !if(iproc==0) write(*,'(1x,a)') 'Iterative diagonalization...'

      if(present(orbsv)) then
         write(*,'(a)') 'ERROR: Virtual orbitals cannot be handled with the iterative input guess at the moment.'
         write(*,'(a)') "Change the value of input%directDiag in 'input.perf' to 'T'."
         stop
      end if
      call inputguessParallel(iproc, nproc, orbs, norbsc_arr, hamovr, &
         &   psi, psit, orthpar, nspin, nspinor, npsidim, comms, natsc, ndim_hamovr, norbsc)

   else

      if(iproc==0) write(*,'(1x,a)') 'Direct diagonalization...'

      call timing(iproc, 'Input_comput', 'ON')

      ispsi=1
      !it is important that the k-points repartition of the inputguess orbitals
      !coincides with the one of the SCF orbitals
      do ikptp=1,orbsu%nkptsp
         ikpt=orbsu%iskpts+ikptp!orbs%ikptsp(ikptp)
         call solve_eigensystem(norbi_max,&
            &   ndim_hamovr,sum(norbgrp),natsceff,nspin,nspinor,norbgrp,hamovr(1,1,ikpt),&
            &   orbsu%eval((ikpt-1)*orbsu%norb+1)) !changed from orbs

         !assign the value for the orbital
         call vcopy(orbs%norbu,orbsu%eval((ikpt-1)*orbsu%norb+1),1,&
              orbs%eval((ikpt-1)*orbs%norb+1),1)
         if (orbs%norbd >0) then
            call vcopy(orbs%norbd,orbsu%eval((ikpt-1)*orbsu%norb+orbsu%norbu+1),1,orbs%eval((ikpt-1)*orbs%norb+orbs%norbu+1),1)
         end if
         !do iorb=1,orbs%norbu
         !   orbs%eval((ikpt-1)*orbs%norb+iorb)=orbsu%eval((ikpt-1)*orbsu%norb+iorb)
         !end do
         !case for collinear spin
         !do iorb=1,orbs%norbd
         !   orbs%eval((ikpt-1)*orbs%norb+iorb+orbs%norbu)=orbsu%eval((ikpt-1)*orbsu%norb+iorb+orbsu%norbu)
         !end do
      end do

      !broadcast values for k-points 
      call broadcast_kpt_objects(nproc, orbsu%nkpts, orbsu%norb, &
           orbsu%eval(1), orbsu%ikptproc)

      if (iproc ==0) then 
         call write_ig_eigenvectors(tolerance,orbsu,nspin,orbs%norb,orbs%norbu,orbs%norbd)
      end if
      !!$  !not necessary anymore since psivirt is gaussian
      !allocate the pointer for virtual orbitals
      if(present(orbsv) .and. present(psivirt)) then
         if (orbsv%norb > 0) then
            allocate(psivirt(max(orbsv%npsidim_orbs,orbsv%npsidim_comp)+ndebug),stat=i_stat)
            call memocc(i_stat,psivirt,'psivirt',subname)
         end if
      else if(present(psivirt)) then
         if (orbsv%norb == 0) then
            allocate(psivirt(1+ndebug),stat=i_stat)
            call memocc(i_stat,psivirt,'psivirt',subname)
         end if
      end if

      if (iproc == 0 .and. verbose > 1) write(*,'(1x,a)',advance='no')'Building orthogonal Wavefunctions...'
      !n(c) nvctr=wfd%nvctr_c+7*wfd%nvctr_f

      ispsi=1
      ispsie=1
      ispsiv=1
      ispm=1
      !number of complex components
      ncplx=1
      if (nspinor > 1) ncplx=2
      do ikptp=1,orbsu%nkptsp
         ikpt=orbsu%iskpts+ikptp!orbsu%ikptsp(ikptp)

         nvctrp=commu%nvctr_par(iproc,ikpt )
         if (nvctrp == 0) cycle

         if (.not. present(orbsv)) then
            call build_eigenvectors(orbs%norbu,orbs%norbd,orbs%norb,norbtot,nvctrp,&
               &   natsceff,nspin,nspinor,orbs%nspinor,ndim_hamovr,norbgrp,hamovr(1,1,ikpt),&
               &   psi(ispsie:),psit(ispsi:),passmat(ispm))
         else
            call build_eigenvectors(orbs%norbu,orbs%norbd,orbs%norb,norbtot,nvctrp,&
               &   natsceff,nspin,nspinor,orbs%nspinor,ndim_hamovr,norbgrp,hamovr(1,1,ikpt),&
               &   psi(ispsie:),psit(ispsi:),passmat(ispm),&
               &   (/orbsv%norbu,orbsv%norbd/),psivirt(ispsiv:))
         end if
         ispsi=ispsi+nvctrp*orbs%norb*orbs%nspinor
         ispsie=ispsie+nvctrp*norbtot*orbs%nspinor
         ispm=ispm+ncplx*(orbsu%norbu*orbs%norbu+orbsu%norbd*orbs%norbd)
         !print *,'iproc,nkptsp,',iproc,orbsu%nkptsp,ispm
         if (present(orbsv)) ispsiv=ispsiv+nvctrp*orbsv%norb*orbs%nspinor
      end do

      !if(nproc==1.and.nspinor==4) call psitransspi(nvctrp,norbu+norbd,psit,.false.)
      !if (iproc == 0 .and. verbose > 1) write(*,'(1x,a)') 'done.'

      if(present(psivirt)) then
         if (orbsv%norb == 0) then
            i_all=-product(shape(psivirt))*kind(psivirt)
            deallocate(psivirt,stat=i_stat)
            call memocc(i_stat,i_all,'psivirt',subname)
         end if
      end if
      call timing(iproc, 'Input_comput', 'OF')

   end if differentInputGuess


   i_all=-product(shape(hamovr))*kind(hamovr)
   deallocate(hamovr,stat=i_stat)
   call memocc(i_stat,i_all,'hamovr',subname)
   i_all=-product(shape(norbgrp))*kind(norbgrp)
   deallocate(norbgrp,stat=i_stat)
   call memocc(i_stat,i_all,'norbgrp',subname)

   !print * ,'debug2,iproc',iproc,orbsv%norb,orbsv%norbp,orbsv%norbu,orbsv%norbd,orbsv%npsidim
   if (minimal) then
      !deallocate the old psi
      i_all=-product(shape(psi))*kind(psi)
      deallocate(psi,stat=i_stat)
      call memocc(i_stat,i_all,'psi',subname)
   else if (nproc == 1) then
      !!$     !reverse objects for the normal diagonalisation in serial
      !!$     !at this stage hpsi is the eigenvectors and psi is the old wavefunction
      !!$     !this will restore the correct identification
      !!$     nullify(hpsi)
      !!$     hpsi => psi
      !!$     !     if(nspinor==4) call psitransspi(nvctrp,norb,psit,.false.) 
      !!$     nullify(psi)
      !!$     psi => psit
   end if

   !orthogonalise the orbitals in the case of semi-core atoms
   if (norbsc > 0) then
      call orthogonalize(iproc,nproc,orbs,comms,psit,orthpar)
   end if
   if (minimal) then
      allocate(hpsi(max(orbs%npsidim_orbs,orbs%npsidim_comp)+ndebug),stat=i_stat)
      call memocc(i_stat,hpsi,'hpsi',subname)
      !     hpsi=0.0d0
      if (nproc > 1) then
         !allocate the direct wavefunction
         allocate(psi(max(orbs%npsidim_orbs,orbs%npsidim_comp)+ndebug),stat=i_stat)
         call memocc(i_stat,psi,'psi',subname)
      else
         psi => psit
      end if
   end if
   !this untranspose also the wavefunctions 
   call untranspose_v(iproc,nproc,orbs,wfd,comms,&
      &   psit,work=hpsi,outadd=psi(1))

   if (nproc == 1 .and. minimal) then
      nullify(psit)
   end if

END SUBROUTINE DiagHam
