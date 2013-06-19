! Needs to be cleaned up, changed the data structure of much Lzd !!!
!!subroutine LinearDiagHam(iproc,at,etol,Lzd,orbs,nspin,natsc,Lhpsi,Lpsi,psit,orbsv,norbsc_arr,orbse,commse)
!!  use module_base
!!  use module_types
!!  use module_interfaces, exceptThisOne => LinearDiagHam
!!  implicit none
!!  integer, intent(in) :: iproc                                          !> current processor
!!  integer, intent(in) :: nspin                                          !> number of spin in spin polarized calculation
!!  integer, intent(in) :: natsc                                          !> number of atoms having semicore states
!!  real(gp),intent(in) :: etol                                           !> Tolerance on energy for solve_eigensystem
!!  type(atoms_data),intent(in) :: at                                     !> Atoms data
!!  type(local_zone_descriptors) :: Lzd                                  !> Information about the locregs
!!  type(orbitals_data), intent(in) :: orbs                               !> description of orbitals after the input guess (less orbitals)
!!  type(orbitals_data), optional, intent(in) :: orbsv                    !> description of virtual orbitals (for davidson?)
!!  real(wp),dimension(Lzd%Lpsidimtot),intent(in):: Lhpsi                   !> All the local H|Psi>
!!  real(wp),dimension(Lzd%Lpsidimtot),intent(in):: Lpsi                    !> All the local |Psi>
!!  real(wp),dimension(orbs%npsidim),intent(inout):: psit                 !> Eigenvectors
!!  integer, optional, dimension(natsc+1,nspin), intent(in) :: norbsc_arr !> semicore description
!!  type(orbitals_data),intent(in),optional,target :: orbse                      !> description of orbitals for input guess
!!  type(communications_arrays), optional, target, intent(in) :: commse   !> communicators
!!  !Local variables
!!  integer :: ilr,ilr2,ikptp,ikpt,ii,jj                       !> loop integers
!!  integer :: i_stat,i_all                                    !> error handling for allocation/deallocation and memocc
!!  integer :: ndim_hamovr                                     !> dimension of Global Hamiltonian/overlap matrix
!!  integer :: psishift1,psishift2                             !> shifting index of wavefunctions for locreg(1) and locreg(2)
!!  integer :: psidim1,psidim2                                 !> dimension of wavefunctions in locreg(1) and locreg(2)
!!  integer :: firstrow,lastrow                                !> index of first and last row of hamovr calculted locally
!!  integer :: firstcol,lastcol                                !> index of first and last col of hamovr calculted locally
!!  integer :: isovrlp                                         !> number of overlap between locreg(1) and locreg(2) [more then one because of periodic BC]
!!  integer :: dim_Lhamovr                                     !> dimension of the local, to locreg(1), hamiltonian/overlap matrix
!!  integer :: norbi_max                                       !> Maximum number of orbitals
!!  integer :: natsceff,ispsi,norbsc,ldim                      !> 
!!  integer :: nvctrp,norbtot,ispsie,ispsiv                    !>
!!  integer :: totshift                                        !> Total shift for the wavefunction when passing from local_to_Global
!!  integer :: Gpsidim                                         !> Global wavefunction dimension
!!  integer :: noncoll                                         !> =2 for non collinear spins, =1 otherwise
!!  integer :: norb1,norb2
!!  integer :: i,ispin,ndh1,iat,iel,scshift,scstr              !> used for loops with semicore states
!!  integer :: npsidim,nspinor
!!  type(orbitals_data),pointer :: orbsu                               !> used to define the orbitals
!!  logical :: minimal
!!  character(len=*),parameter :: subname='Linear_DiagHam'     ! name of subroutine
!!  integer, dimension(:,:), allocatable :: norbgrp            !>
!!  real(wp), dimension(:), pointer :: psivirt                 !> 
!!  real(wp),dimension(:,:,:),allocatable :: hamovr            !> Global Hamiltonian/overlap matrix
!!  integer,dimension(5) :: sizes                              !> array with the sizes for the reshaping of hamovr
!!  integer,dimension(Lzd%nlr,nspin) :: orbscToAtom            !> mapping the number of semicore orbitals on the atoms
!!  real(wp), dimension(:), pointer :: psi                     !> global quantities for now   SHOULD BE REPLACED
!!  real(wp),dimension(:,:,:,:,:),allocatable :: work1  !> working arrays for hamovr
!!  real(wp),dimension(:,:,:),allocatable :: Lhamovr           !> Local,to locreg(1), hamiltonian/overlap matrix
!!  logical :: semicore
!!
!!  semicore=present(norbsc_arr)
!!
!!  if (present(orbse) .neqv. present(commse)) then
!!     !if (iproc ==0) 
!!           write(*,'(1x,a)')&
!!          'ERROR (DiagHam): the variables orbse and commse must be present at the same time'
!!     stop
!!  else
!!     minimal=present(orbse)
!!  end if
!!
!!  if(minimal) then
!!     norbtot=orbse%norb !beware that norbe is equal both for spin up and down
!!     orbsu => orbse
!!     npsidim=orbse%npsidim
!!     nspinor=orbse%nspinor
!!  else
!!   write(*,*) 'Linear InputGuess only with minimal.'
!!   stop
!!  end if
!!
!!
!!  if (iproc == 0 .and. verbose > 1) write(*,'(1x,a)',advance='no')&
!!        'Overlap Matrix...'
!!
!!  if (orbs%nspinor == 4) then
!!     noncoll=2
!!  else
!!     noncoll=1
!!  end if
!!
!!  if (semicore) then
!!     if (present(orbsv)) then
!!        norbi_max=max(noncoll*maxval(norbsc_arr),orbsv%norb)
!!     else
!!        norbi_max=noncoll*maxval(norbsc_arr)
!!     end if
!!
!!     scshift = 0
!!     do ii=1,natsc
!!        scshift = scshift + (noncoll*norbsc_arr(ii,1))**2
!!     end do
!!
!!     !calculate the dimension of the overlap matrix
!!     !take the maximum as the two spin dimensions
!!     ndim_hamovr=0
!!     do ispin=1,nspin
!!        ndh1=0
!!        norbsc=0
!!        do i=1,natsc+1
!!           ndh1=ndh1+(noncoll*norbsc_arr(i,ispin))**2
!!        end do
!!        ndim_hamovr=max(ndim_hamovr,ndh1)
!!     end do
!!     if (natsc > 0) then
!!        if (nspin == 2) then
!!           if (sum(norbsc_arr(1:natsc,1)) /= sum(norbsc_arr(1:natsc,2))) then
!!              write(*,'(1x,a)')&
!!                'ERROR (DiagHam): The number of semicore orbitals must be the same for both spins'
!!              stop
!!           end if
!!        end if
!!        norbsc=noncoll*sum(norbsc_arr(1:natsc,1))
!!     else
!!        norbsc=0
!!     end if
!!  else
!!     norbsc = 0
!!     norbi_max=max(orbse%norbu,orbse%norbd)
!!     ndim_hamovr = norbi_max**2
!!  end if
!!
!!  !allocation of Global hamiltonian/overlap matrix
!!  allocate(hamovr(nspin*ndim_hamovr,2,orbse%nkpts+ndebug),stat=i_stat)
!!  call memocc(i_stat,hamovr,'hamovr',subname)
!!
!!  ! put it to zero
!!  call razero(nspin*ndim_hamovr*2*orbse%nkpts,hamovr)
!!
!!  ! First build the number of semicore states for each atom (hence each locreg)
!!  ! ONLY WORKS FOR INPUT GUESS where each atom has it's locreg
!!  orbscToAtom = 0
!!  do ispin=1,nspin
!!     iat = 0
!!     do ilr=1,Lzd%nlr
!!        if(at%iasctype(ilr) == 0)cycle
!!         iat = iat + 1
!!        orbscToAtom(ilr,ispin) = norbsc_arr(iat,ispin)
!!     end do
!!  end do
!!
!!  psishift1 = 1
!!  firstrow  = 1 + norbsc  !This assumes that number of semicore orbitals are the same for each spin. CHECK IF NONCOLL is a problem
!!  scstr = 1
!!  ! The loop on ilr gives the row indexes, the loop on ilr2 gives the column indexes
!!  do ilr = 1, Lzd%nlr
!!     norb1 = Lzd%Llr(ilr)%Localnorb/orbsu%nspin
!!     if (orbscToAtom(ilr,1) > 0) then  !This assumes that number of semicore orbitals are the same for each spin.
!!        ! Calculate the hamiltonian/overlap matrix for the semicore states of locreg ilr
!!        ! These states only intersect with themselves, so only need to consider this locreg
!!        call semicore_overlap_matrices(ilr,nspin,nspinor,orbse%norb,Lzd,orbsu,orbscToAtom,&
!!             ndim_hamovr,hamovr,scstr,Lpsi,Lhpsi)
!!        scstr= scstr + (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbscToAtom(ilr,1)*nspinor
!!        psishift1 = psishift1 + (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbscToAtom(ilr,1)*nspinor
!!     end if
!!
!!! DEBUG
!!!  print *,'size(hamovr)',size(hamovr,1),size(hamovr,2),size(hamovr,3)
!!!  do i_all=1,size(hamovr,1)
!!!     print *,'iel, ham, ovr:',i_all,hamovr(i_all,1,:),hamovr(i_all,2,:)
!!!  end do
!!! END DEBUG
!!
!!     !Now calculate the hamiltonian/overlap matrix for the non semicore states (FOR SPINS THIS DOES NOT WORK)
!!     firstcol = 1 + norbsc   !This assumes that number of semicore orbitals are the same for each spin. CHECK IF NONCOLL is a problem
!!     lastrow  = firstrow-1  + norb1 - orbscToAtom(ilr,1) !This assumes that number of semicore orbitals are the same for each spin.
!!     psidim1 = (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*&
!!               (norb1-orbscToAtom(ilr,1))*nspinor
!!     psishift2 = 1
!!     do ilr2 = 1,Lzd%nlr
!!        norb2 = Lzd%Llr(ilr2)%Localnorb/orbsu%nspin
!!        ! don't use the semicore states
!!        psishift2 = psishift2 + (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*&
!!                    orbscToAtom(ilr2,1)*nspinor
!!
!!        call get_number_of_overlap_region(ilr,ilr2,Lzd%Glr,isovrlp,Lzd%Llr,Lzd%nlr)!,outofzone)
!!
!!        psidim2=(Lzd%Llr(ilr2)%wfd%nvctr_c+7*Lzd%Llr(ilr2)%wfd%nvctr_f)*norb2*nspinor
!!
!!        ! If no overlap, increment the index of Lpsi and overlap matrix then cycle
!!        if(isovrlp == 0)then
!!           psishift2 = psishift2 + psidim2*nspin
!!           lastcol  = firstcol-1  + norb2 - orbscToAtom(ilr2,1)
!!           firstcol = firstcol + norb2 - orbscToAtom(ilr2,1)
!!           cycle
!!        end if
!!
!!        !redefine psidim2 to remove the semicore states
!!        psidim2=(Lzd%Llr(ilr2)%wfd%nvctr_c+7*Lzd%Llr(ilr2)%wfd%nvctr_f)*(norb2-orbscToAtom(ilr2,1))*orbse%nspinor
!!
!!        ! dimensions and allocation of Local hamiltonian/overlap matrix
!!        dim_Lhamovr=(norb1-orbscToAtom(ilr,1))*(norb2-orbscToAtom(ilr2,1))
!!        allocate(Lhamovr(nspin*dim_Lhamovr,2,orbse%nkpts+ndebug),stat=i_stat)
!!        call memocc(i_stat,Lhamovr,'Lhamovr',subname)
!!
!!     ! In this routine, we begin by calculating the hamiltonian/overlap matrix between two locregs.
!!        call overlap_matrix_between_locreg(ilr,ilr2,isovrlp,nspin,orbscToAtom,psidim1,psidim2,psishift1,&
!!           psishift2,Lzd,orbsu,Lpsi,Lhpsi,dim_Lhamovr,Lhamovr)
!!
!!! DEBUG
!!!  print *,'size(Lhamovr)',size(Lhamovr,1),size(Lhamovr,2),size(Lhamovr,3)
!!!  do i_all=1,size(Lhamovr,1)
!!!     print *,'iel, ham, ovr:',i_all,Lhamovr(i_all,1,:),Lhamovr(i_all,2,:)
!!!  end do
!!! END DEBUG
!!
!!     ! update the shift for second wavefunction
!!        psishift2 = psishift2 + psidim2*nspin
!!
!!     ! reshape the hamiltonian/overlap matrix for easy assignations  
!!       allocate(work1(norb1-orbscToAtom(ilr,1),norb2-orbscToAtom(ilr2,1),nspin,2,orbse%nkpts+ndebug),stat=i_stat)
!!       call memocc(i_stat,work1,'work1',subname)
!!       sizes = (/ norb1-orbscToAtom(ilr,1),norb2-orbscToAtom(ilr2,1),nspin, 2, orbse%nkpts+ndebug /)
!!       work1 = reshape(Lhamovr,sizes)
!!
!!     ! Assign the calculated values inside global matrix (for truly O(N) this should be replaced) 
!!       lastcol  = firstcol-1  + norb2 - orbscToAtom(ilr2,1)
!!
!!
!!       do ispin=1,nspin
!!          iel = (firstrow-norbsc-1)*norbsc_arr(natsc+1,1)+(firstcol-norbsc-1)+ispin*scshift+(ispin-1)*norbsc_arr(natsc+1,1)**2
!!          do ii=1,size(work1,1)
!!             do jj=1,size(work1,2)
!!                hamovr(iel+jj,:,:) = work1(ii,jj,ispin,:,:)
!!             end do
!!             iel = iel + norbsc_arr(natsc+1,1)
!!          end do
!!       end do
!!
!!     ! deallocate this instance of Lhamovr
!!        i_all=-product(shape(work1))*kind(work1)
!!        deallocate(work1,stat=i_stat)
!!        call memocc(i_stat,i_all,'work1',subname)
!!
!!        i_all=-product(shape(Lhamovr))*kind(Lhamovr)
!!        deallocate(Lhamovr,stat=i_stat)
!!        call memocc(i_stat,i_all,'Lhamovr',subname)
!!
!!     ! update indexes
!!       firstcol = firstcol + norb2 - orbscToAtom(ilr2,1)
!!     end do
!!     ! increment the shift of wavefunctions
!!     psishift1 = psishift1 + psidim1*nspin
!!     firstrow = firstrow + norb1 - orbscToAtom(ilr,1)
!!  end do
!!
!!
!!! DEBUG
!!  print *,'size(hamovr)',size(hamovr,1),size(hamovr,2),size(hamovr,3)
!!  do i_all=1,size(hamovr,1)
!!     print *,'iel, ham, ovr:',i_all,hamovr(i_all,1,:),hamovr(i_all,2,:)
!!  end do
!!! END DEBUG
!!
!!  ! Don't need Lhpsi anymore
!!!  i_all=-product(shape(Lhpsi))*kind(Lhpsi)
!!!  deallocate(Lhpsi,stat=i_stat)
!!!  call memocc(i_stat,i_all,'Lhpsi',subname)
!!
!!  ! Now solve the eigensystem: H |Lpsi> = epsilon S |Lpsi>  
!!  if(iproc==0) write(*,'(1x,a)') 'Direct diagonalization...'
!!
!!!  call timing(iproc, 'Input_comput', 'ON')
!!
!!  ! SET SOME VARIABLE FOR NOW (NO SEMICORE)
!!  ispsi=1
!!  if (semicore) then
!!     natsceff=natsc
!!     allocate(norbgrp(natsceff+1,nspin+ndebug),stat=i_stat)
!!     call memocc(i_stat,norbgrp,'norbgrp',subname)
!!
!!     !assign the grouping of the orbitals
!!     do jj=1,nspin
!!        do ii=1,natsceff+1
!!           norbgrp(ii,jj)=noncoll*norbsc_arr(ii,jj)
!!        end do
!!     end do
!!  else
!!     natsceff = 0
!!     allocate(norbgrp(1,nspin+ndebug),stat=i_stat)
!!     call memocc(i_stat,norbgrp,'norbgrp',subname)
!!     norbsc=0
!!     norbgrp(1,1)=orbse%norbu
!!     if (nspin == 2) norbgrp(1,2)=orbse%norbd
!!  end if
!!  !it is important that the k-points repartition of the inputguess orbitals
!!  !coincides with the one of the SCF orbitals
!!  ! NOTE: USE ORBS OR ORBSE??
!!  do ikptp=1,orbse%nkptsp
!!     ikpt=orbse%iskpts+ikptp!orbs%ikptsp(ikptp)
!!     call solve_eigensystem(iproc,orbs%norb,orbs%norbu,orbs%norbd,norbi_max,&
!!          ndim_hamovr,natsceff,nspin,orbse%nspinor,etol,norbgrp,hamovr(1,1,ikpt),&
!!          orbs%eval((ikpt-1)*orbs%norb+1))
!!  end do
!!
!!  if (iproc == 0 .and. verbose > 1) write(*,'(1x,a)',advance='no')'Building orthogonal Wavefunctions...'
!!
!!
!!! FOR NOW, just transform the Lpsi to psi in global region.
!!  Gpsidim = (Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*orbse%norb*orbse%nspinor
!!  allocate(psi(Gpsidim+ndebug),stat=i_stat)
!!  call memocc(i_stat,psi,'psi',subname)
!!  call razero(Gpsidim+ndebug,psi)
!!
!!! WATCH OUT, does not work for nspinor > 1
!!  psishift1 = 1
!!  totshift = 0
!!  do ilr = 1,Lzd%nlr
!!     norb1 = Lzd%Llr(ilr)%Localnorb/orbsu%nspin
!!     ldim = (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*Lzd%Llr(ilr)%Localnorb*nspinor
!!     call Lpsi_to_global(Lzd%Glr,Gpsidim,Lzd%Llr(ilr),Lpsi(psishift1),&
!!          ldim,norb1,nspinor,nspin,totshift,psi)
!!     psishift1 = psishift1 + ldim
!!     totshift = totshift + (Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*norb1*orbse%nspinor
!!  end do
!!
!!!  allocate(psit(orbs%npsidim+ndebug),stat=i_stat)
!!!  call memocc(i_stat,psit,'psit',subname)
!!
!!  ispsi=1
!!  ispsie=1
!!  ispsiv=1
!!  norbtot = orbse%norb
!!  do ikptp=1,orbse%nkptsp
!!     ikpt=orbse%iskpts+ikptp!orbsu%ikptsp(ikptp)\
!!
!!!    nvctrp is not a simple quantity anymore has it depends on the locregs (can be different for every locreg)
!!!    for an O(N) code, should change these routines.
!!!    FOR NOW, just transform the Lpsi to psi in global region.
!!     nvctrp=Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f
!!     if (nvctrp == 0) cycle
!!     call build_eigenvectors(orbs%norbu,orbs%norbd,orbs%norb,norbtot,nvctrp,&
!!          natsceff,nspin,orbse%nspinor,orbse%nspinor,ndim_hamovr,norbgrp,hamovr(1,1,ikpt),&
!!          psi(ispsie:),psit(ispsi:))
!!     ispsi=ispsi+nvctrp*orbs%norb*orbse%nspinor
!!     ispsie=ispsie+nvctrp*norbtot*orbse%nspinor
!!  end do
!!
!!
!!  !if(nproc==1.and.nspinor==4) call psitransspi(nvctrp,norbu+norbd,psit,.false.)
!!  if (iproc == 0 .and. verbose > 1) write(*,'(1x,a)') 'done.'
!!
!!  !deallocate psi
!!  i_all=-product(shape(psi))*kind(psi)
!!  deallocate(psi,stat=i_stat)
!!  call memocc(i_stat,i_all,'psi',subname)
!!
!!end subroutine LinearDiagHam

