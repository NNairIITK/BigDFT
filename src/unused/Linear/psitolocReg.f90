!> Tranform wavefunction between Global region and localisation region
!! @warning 
!!     Only coded for sequential, not parallel cases !! For parallel should change increment and loc_psi dimensions
!!subroutine psi_to_locreg(Glr,ilr,ldim,Olr,lpsi,nlr,orbs,psi)
!!
!!  use module_base
!!  use module_types
!! 
!! implicit none
!!
!!  ! Subroutine Scalar Arguments
!!  integer, intent(in) :: nlr                  ! number of localization regions
!!  integer :: ilr           ! index of the localization region we are considering
!!  integer :: ldim          ! dimension of lpsi 
!!  type(orbitals_data),intent(in) :: orbs      ! orbital descriptor
!!  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
!!  
!!  !Subroutine Array Arguments
!!  type(locreg_descriptors), dimension(nlr), intent(in) :: Olr  ! Localization grid descriptors 
!!  real(wp),dimension((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*orbs%norbp*orbs%nspinor),intent(in) :: psi       !Wavefunction (compressed format)
!!  real(wp),dimension(ldim),intent(inout) :: lpsi !Wavefunction in localization region
!!  
!!  !local variables
!!  integer :: igrid,isegloc,isegG,ix,iorbs
!!  integer :: lmin,lmax,Gmin,Gmax
!!  integer :: icheck      ! check to make sure the dimension of loc_psi does not overflow 
!!  integer :: offset      ! gives the difference between the starting point of Lseg and Gseg
!!  integer :: length      ! Length of the overlap between Lseg and Gseg
!!  integer :: lincrement  ! Increment for writing orbitals in loc_psi
!!  integer :: Gincrement  ! Increment for reading orbitals in psi
!!  integer :: nseg        ! total number of segments in Llr
!!  integer, allocatable :: keymask(:,:)  ! shift for every segment of Llr (with respect to Glr)
!!  character(len=*), parameter :: subname='psi_to_locreg'
!!  integer :: i_stat,i_all
!!  integer :: start,Gstart
!!  integer :: lfinc,Gfinc
!!
!!! Define integers
!!  nseg = Olr(ilr)%wfd%nseg_c + Olr(ilr)%wfd%nseg_f
!!  lincrement = Olr(ilr)%wfd%nvctr_c + 7*Olr(ilr)%wfd%nvctr_f
!!  Gincrement = Glr%wfd%nvctr_c + 7*Glr%wfd%nvctr_f
!!  icheck = 0
!!
!!! Initialize loc_psi
!!  if (ldim > 0) call to_zero(lincrement*orbs%norbp*orbs%nspinor,lpsi(1))
!!
!!! Get the keymask: shift for every segment of Llr (with respect to Glr)
!!  allocate(keymask(2,nseg),stat=i_stat)
!!  call memocc(i_stat,keymask,'keymask',subname)
!!
!!  call shift_locreg_indexes(Glr,Olr(ilr),keymask,nseg)
!!
!!!####################################################
!!! Do coarse region
!!!####################################################
!!  do isegloc = 1,Olr(ilr)%wfd%nseg_c
!!     lmin = keymask(1,isegloc)
!!     lmax = keymask(2,isegloc)
!! 
!!! Could optimize the routine by looping only on Gsegs not looped on before (TO DO)
!!     do isegG = 1,Glr%wfd%nseg_c
!!        Gmin = Glr%wfd%keygloc(1,isegG)
!!        Gmax = Glr%wfd%keygloc(2,isegG)
!!
!!        ! For each segment in Llr check if there is a collision with the segment in Glr
!!        ! if not, cycle
!!        if((lmin > Gmax) .or. (lmax < Gmin)) cycle
!!        
!!        ! Define the offset between the two segments
!!        offset = lmin - Gmin
!!        if(offset < 0) then
!!           offset = 0
!!        end if
!!    
!!        ! Define the length of the two segments
!!        length = min(lmax,Gmax)-max(lmin,Gmin)
!! 
!!        !Find the common elements and write them to the new localized wavefunction
!!        ! WARNING: index goes from 0 to length because it is the offset of the element
!!        do ix = 0,length
!!           icheck = icheck + 1
!!           ! loop over the orbitals
!!           do iorbs=1,orbs%norbp*orbs%nspinor
!!              lpsi(icheck+lincrement*(iorbs-1))=psi(Glr%wfd%keyvloc(isegG)+offset+ix+Gincrement*(iorbs-1))
!!           end do
!!        end do
!!     end do
!!  end do
!!
!!! Check if the number of elements in loc_psi is valid
!!  if(icheck .ne. Olr(ilr)%wfd%nvctr_c) then
!!    write(*,*)'There is an error in psi_to_locreg: number of coarse points used',icheck
!!    write(*,*)'is not equal to the number of coarse points in the region',Olr(ilr)%wfd%nvctr_c
!!  end if
!!
!!!##############################################################
!!! Now do fine region
!!!##############################################################
!!
!!  icheck = 0
!!  start = Olr(ilr)%wfd%nvctr_c
!!  Gstart = Glr%wfd%nvctr_c
!!  lfinc  = Olr(ilr)%wfd%nvctr_f
!!  Gfinc = Glr%wfd%nvctr_f
!!
!!  do isegloc = Olr(ilr)%wfd%nseg_c+1,nseg
!!     lmin = keymask(1,isegloc)
!!     lmax = keymask(2,isegloc)
!! 
!!! Could optimize the routine by looping only on Gsegs not looped on before (TO DO)
!!     do isegG = Glr%wfd%nseg_c+1,Glr%wfd%nseg_c+Glr%wfd%nseg_f
!!
!!        Gmin = Glr%wfd%keygloc(1,isegG)
!!        Gmax = Glr%wfd%keygloc(2,isegG)
!!
!!        ! For each segment in Llr check if there is a collision with the segment in Glr
!!        ! if not, cycle
!!        if((lmin > Gmax) .or. (lmax < Gmin)) cycle
!!
!!        offset = lmin - Gmin
!!        if(offset < 0) offset = 0
!!
!!        length = min(lmax,Gmax)-max(lmin,Gmin)
!!
!!        !Find the common elements and write them to the new localized wavefunction
!!        ! WARNING: index goes from 0 to length because it is the offset of the element
!!        do ix = 0,length
!!           icheck = icheck + 1
!!           do igrid=0,6
!!              do iorbs=1,orbs%norbp*orbs%nspinor
!!                 lpsi(start+icheck+lincrement*(iorbs-1)+igrid*lfinc)=&
!!&                psi(Gstart+Glr%wfd%keyvloc(isegG)+offset+ix+Gincrement*(iorbs-1)+igrid*Gfinc)
!!              end do
!!           end do
!!        end do
!!     end do
!!  end do
!!  
!! ! Check if the number of elements in loc_psi is valid
!!  if(icheck .ne. Olr(ilr)%wfd%nvctr_f) then
!!    write(*,*)'There is an error in psi_to_locreg: number of fine points used',icheck
!!    write(*,*)'is not equal to the number of fine points in the region',Olr(ilr)%wfd%nvctr_f
!!  end if
!!
!!  i_all=-product(shape(keymask))*kind(keymask)
!!  deallocate(keymask,stat=i_stat)
!!  call memocc(i_stat,i_all,'keymask',subname)
!!
!!END SUBROUTINE psi_to_locreg

!!$subroutine overlap_matrix_for_locreg(ilr,nlr,nspin,psidim1,psishift1,npsidim,at,orbsi,Glr,Llr,Lpsi,&
!!$           Lhpsi,Localnorb,Lnorbovr,outofzone,dim_Lhamovr,Lhamovr)
!!$
!!$  use module_base
!!$  use module_types
!!$
!!$ implicit none
!!$
!!$  ! Subroutine Scalar Arguments
!!$  integer,intent(in) :: ilr          ! index of current locreg
!!$  integer,intent(in) :: nlr          ! number of localization regions 
!!$  integer,intent(in) :: nspin        ! number of spins
!!$  integer,intent(in) :: psidim1      ! dimension of the wavefunctions in locreg(ilr)
!!$  integer,intent(in) :: psishift1    ! starting index of the first orbital of locreg(ilr)
!!$  integer,intent(in) :: npsidim      ! total dimension of the wavefunction
!!$  integer,intent(in) :: dim_Lhamovr  ! dimension of the Local Hamiltonian/Overlap Matrix
!!$  integer,intent(in) :: Lnorbovr     ! number of orbitals in Local overlap/hamiltonian matrix
!!$  type(atoms_data), intent(in) :: at ! atoms data
!!$  type(orbitals_data),intent(in) :: orbsi      ! orbital descriptor   
!!$  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
!!$  type(locreg_descriptors),dimension(nlr), intent(in) :: Llr  ! Localization grid descriptors 
  
!!$  !Subroutine Array Arguments
!!$  real(wp),dimension(npsidim),intent(in) :: Lpsi       ! Wavefunction (compressed format)
!!$  real(wp),dimension(npsidim),intent(in) :: Lhpsi      ! Wavefunction in localization region
!!$  integer,dimension(nlr),intent(in) :: Localnorb       ! Number of orbitals in each locreg
!!$  integer,dimension(3,nlr),intent(in) :: outofzone     ! Periodicity of the localization regions
!!$  real(wp),dimension(dim_Lhamovr,2,orbsi%nkpts),intent(out) :: Lhamovr           ! Local Hamiltonian/Overlap matrix
  
!!$  ! Local Variables
!!$  integer :: ilr2,iolr,i_all,i_stat,psidim2,psishift2
!!$  integer :: isovrlp           ! number of overlap regions (periodicity)
!!$  integer :: ndim_hamovr       ! dimension of hamiltonian and overlap matrix
!!$  integer :: ldim,ldim1,ldim2
!!$  integer :: shiftdiag,shiftnd,shift1,shift2,ind
!!$  integer :: i1,i2,i3,lspin,indspin,ispin
!!$  integer :: fmin
!!$  type(locreg_descriptors),dimension(:), allocatable :: Olr  ! Localization grid descriptors 
!!$  integer, dimension(1,nspin) :: norbgrp    !assign orbtial to groups TO DO : check for semi-core states
!!$  real(wp), dimension(:,:,:), allocatable :: Lohamovr,Ahamovr,OnSitehamovr  !hamiltonian/overlap matrices
!!$  real(wp),dimension(:),pointer:: Lopsi,Lohpsi
!!$  real(wp) :: scpr
!!$  character(len=*), parameter :: subname='overlap_matrix_for_locreg'
!!$  type(orbitals_data) :: orbs      ! orbital descriptor   
!!$
!!$  orbs = orbsi 
!!$
!!$  psishift2 = 1
!!$  do ilr2 = 1, nlr
!!$     print *,'Looking for overlap of regions:',ilr,ilr2
!!$     psidim2 = (Llr(ilr2)%wfd%nvctr_c+7*Llr(ilr2)%wfd%nvctr_f)*Localnorb(ilr2)*orbs%nspinor
!!$     !Calculate the number of overlap regions between two logregs (can be more then one because
!!$     ! of periodicity). The number of overlap regions is stored in the isvorlp integer.
!!$     call get_number_of_overlap_region(ilr,ilr2,Glr,isovrlp,Llr,nlr,outofzone)
!!$     
!!$     if (isovrlp > 0 .and. (ilr .ne. ilr2)) then
!!$          !allocate the overlap regions (memocc not used because it is a type) 
!!$          allocate(Olr(isovrlp),stat=i_stat)
!!$
!!$          ! assign orbtial to groups (no semi-core)         TO-DO: Semi-core states
!!$          norbgrp(1,1)= (Localnorb(ilr)+Localnorb(ilr2))                            ! up orbitals in first group
!!$          if (nspin == 2) norbgrp(1,2)= (Localnorb(ilr)+Localnorb(ilr2))            ! down orbitals in second group
!!$
!!$          ! assign dimension of Overlap matrix     
!!$          ndim_hamovr=(Localnorb(ilr)+Localnorb(ilr2))**2  !TO-DO : Check spins !!
!!$          ! Now allocate the total overlap matrix for the overlap region [Lohamovr]
!!$          ! and a temporary overlap matrix [Ahamovr] for each zone (Olr fractured because of periodicity)
!!$          ! It has dimensions: Lohamovr(nspin*norb**2,2,nkpt). The second dimension is:
!!$          ! 1: Hamiltonian Matrix =  <Obr(iorb1)|H|Orb(iorb2)>
!!$          ! 2: Overlap Matrix = <Obr(iorb1)|Orb(iorb2)>
!!$          allocate(Lohamovr(nspin*ndim_hamovr,2,orbs%nkpts+ndebug),stat=i_stat)
!!$          call memocc(i_stat,Lohamovr,'Lohamovr',subname)
!!$          allocate(Ahamovr(nspin*ndim_hamovr,2,orbs%nkpts+ndebug),stat=i_stat)    !temporary
!!$          call memocc(i_stat,Ahamovr,'Ahamovr',subname)
!!$
!!$          !since Lohamovr is an accumulator
!!$          call razero(nspin*ndim_hamovr*2*(orbs%nkpts+ndebug),Lohamovr)
!!$
!!$          ! Construct the overlap region descriptors
!!$          call get_overlap_region_periodic(ilr,ilr2,Glr,isovrlp,Llr,nlr,Olr)
!!$
!!$          ! Third, transform the wavefunction to overlap regions
!!$          do iolr=1,isovrlp
!!$            print *,'Treating overlap region :',iolr
!!$            ldim  =  Olr(iolr)%wfd%nvctr_c+7*Olr(iolr)%wfd%nvctr_f
!!$            ldim1 = ldim * Localnorb(ilr) * orbs%nspinor 
!!$            ldim2 = ldim * Localnorb(ilr2)* orbs%nspinor 
!!$
!!$            ! Allocate the local wavefunction (in one overlap region)
!!$            allocate(Lopsi((ldim1+ldim2)*nspin+ndebug), stat=i_stat)
!!$            call memocc(i_stat,lopsi,'lopsi',subname)
!!$            allocate(Lohpsi((ldim1+ldim2)*nspin+ndebug), stat=i_stat)
!!$            call memocc(i_stat,lohpsi,'lohpsi',subname)
!!$
!!$            ! Project the wavefunctions inside the overlap region (first for Llr(ilr)and second for Llr(ilr2))
!!$            ! Spin pas bien ordonnee (doit mettre tout les spin up, ensuite tout les spins down)
!!$            orbs%norbp = Localnorb(ilr) 
!!$            !orbs%npsidim = psidim1 !always from Glr
!!$            do ispin = 1, nspin
!!$               call psi_to_locreg(Llr(ilr),iolr,ldim1,Olr(iolr),&
!!$                         Lopsi(1+(ispin-1)*(ldim1+ldim2):ispin*ldim1+(ispin-1)*ldim2),isovrlp,orbs,&
!!$                         Lpsi(psishift1+(ispin-1)*psidim1:psishift1+ispin*psidim1-1))
!!$               call psi_to_locreg(Llr(ilr),iolr,ldim1,Olr(iolr),&
!!$                        Lohpsi(1+(ispin-1)*(ldim1+ldim2):ispin*ldim1+(ispin-1)*ldim2),isovrlp,orbs,&
!!$                        Lhpsi(psishift1+(ispin-1)*psidim1:psishift1+ispin*psidim1-1))
!!$            end do
!!$            !second region
!!$            orbs%norbp = Localnorb(ilr2)
!!$            !orbs%npsidim = psidim2 !always from Glr
!!$            do ispin = 1, nspin
!!$               call psi_to_locreg(Llr(ilr2),iolr,ldim2,Olr(iolr),&
!!$                        Lopsi(1+ldim1+(ispin-1)*(ldim1+ldim2):ispin*(ldim1+ldim2)),&
!!$                        isovrlp,orbs,Lpsi(psishift2+(ispin-1)*psidim2:psishift2+ispin*psidim2-1))
!!$               call psi_to_locreg(Llr(ilr2),iolr,ldim2,Olr(iolr),&
!!$                        Lohpsi(1+ldim1+(ispin-1)*(ldim1+ldim2):ispin*(ldim1+ldim2)),&
!!$                        isovrlp,orbs,Lhpsi(psishift2+(ispin-1)*psidim2:psishift2+ispin*psidim2-1))
!!$            end do
!!$
!!$            ! Now we can build the overlap matrix of the intersection region
!!$!            ispsi=1                                                     !TO DO: K-points and parallelization
!!$!            do ikptp=1,orbs%nkptsp
!!$!               ikpt=orbs%iskpts+ikptp!orbs%ikptsp(ikptp)
!!$!          
!!$!               nvctrp=commu%nvctr_par(iproc,ikptp)
!!$!               if (nvctrp == 0) cycle
!!$
!!$               !print *,'iproc,nvctrp,nspin,norb,ispsi,ndimovrlp',iproc,nvctrp,nspin,norb,ispsi,ndimovrlp(ispin,ikpt-1)
!!$               call overlap_matrices((Localnorb(ilr)+Localnorb(ilr2))*nspin,ldim,at%natsc,nspin,orbs%nspinor,&
!!$                    & ndim_hamovr,norbgrp,Ahamovr(1,1,1),Lopsi(1),Lohpsi(1))  !Lohamovr(1,1,ikpt),Lopsi(ispsi),Lohpsi(ispsi))
!!$!               ispsi=ispsi+nvctrp*norbtot*orbs%nspinor
!!$!            end do
!!$
!!$            ! The overlap matrix should be accumulated for each intersection region 
!!$              do i1=1,size(Ahamovr,1)
!!$                 do i2=1,size(Ahamovr,2)
!!$                    do i3=1,size(Ahamovr,3)
!!$                       Lohamovr(i1,i2,i3) = Lohamovr(i1,i2,i3) + Ahamovr(i1,i2,i3) !Careful how we accumulate (orbital ordering)                
!!$                     end do
!!$                 end do
!!$              end do
!!$
!!$            ! deallocate the arrays depending on Olr
!!$              i_all=-product(shape(Lopsi))*kind(Lopsi)
!!$              deallocate(Lopsi,stat=i_stat)
!!$              call memocc(i_stat,i_all,'Lopsi',subname)
!!$            
!!$              i_all=-product(shape(Lohpsi))*kind(Lohpsi)
!!$              deallocate(Lohpsi,stat=i_stat)
!!$              call memocc(i_stat,i_all,'Lohpsi',subname)
!!$          end do
!!$          
!!$          ! deallocate Olr
!!$          deallocate(Olr,stat=i_stat)
!!$
!!$     else if (ilr == ilr2) then
!!$          ! allocate overlap matrix inside the locreg (OnSitehamovr)
!!$          ndim_hamovr=Localnorb(ilr)**2
!!$          allocate(OnSitehamovr(nspin*ndim_hamovr,2,orbs%nkpts+ndebug),stat=i_stat)
!!$          call memocc(i_stat,OnSitehamovr,'OnSitehamovr',subname)
!!$
!!$          ! assign orbtial to groups (no semi-core)         TO-DO: Semi-core states
!!$          norbgrp(1,1)= Localnorb(ilr)                      ! up orbitals in first group
!!$          if (nspin == 2) norbgrp(1,2)= Localnorb(ilr)      ! down orbitals in second group
!!$
!!$          ldim1 = (Llr(ilr)%wfd%nvctr_c+7*Llr(ilr)%wfd%nvctr_f) * orbs%nspinor
!!$          call overlap_matrices(Localnorb(ilr)*nspin,ldim1,at%natsc,nspin,orbs%nspinor,ndim_hamovr,norbgrp,&
!!$            OnSitehamovr(1,1,1),Lpsi(psishift1:psishift1+psidim1-1),Lhpsi(psishift1:psishift1+psidim1-1))
!!$           
!!$     end if
!!$
!!$     psishift2 = psishift2 + psidim2*nspin
!!$
!!$     shift1 = 0
!!$     do i1=1,ilr-1
!!$        shift1 = shift1 + Localnorb(i1)
!!$     end do
!!$
!!$     ! Now we must add all the components to form the whole overlap matrix for locreg(ilr)
!!$     ! This only includes the orbitals from intersecting locregs(ilr2)
!!$     ! First add the overlap of wavefunctions inside one locreg(ilr) [OnSitehamovr]
!!$     indspin = 0
!!$     lspin = 0
!!$     do ispin = 1, nspin
!!$        shift2 = 0
!!$        do i1=1,ilr2-1
!!$           shift2 = shift2 + Localnorb(i1)
!!$        end do
!!$        ind = 0
!!$        shiftdiag = 0
!!$        if(ilr > 1) then
!!$           shiftdiag = shiftdiag + (ilr-1) * Lnorbovr * Localnorb(ilr-1) + shift1   ! shift to the correct spot in matrix 
!!$           shift2 = shift2 + (ilr-1) * Lnorbovr * Localnorb(ilr-1) 
!!$        end if
!!$
!!$        shiftnd = 0
!!$        if(ilr2 .ne. ilr ) then
!!$           shiftnd = shiftnd + Localnorb(ilr)
!!$        end if
!!$
!!$        do i1 = 1, Localnorb(ilr)
!!$           do i2 = 1, Localnorb(ilr2)
!!$              ind = ind + 1
!!$              if (isovrlp > 0 .and. ilr .ne. ilr2) then
!!$                 Lhamovr(shift2+i2+indspin,:,:) = Lohamovr(ind+shiftnd,:,:)
!!$              else if (ilr == ilr2) then
!!$                 print *,'shiftdiag+i2+indspin',shiftdiag+i2+indspin,shiftdiag,i2,indspin
!!$                 Lhamovr(shiftdiag+i2+indspin,:,:) = OnSitehamovr(ind+lspin,:,:)
!!$              end if
!!$           end do
!!$           shiftdiag = shiftdiag + Lnorbovr
!!$           shiftnd = shiftnd + Lnorbovr
!!$           shift2 = shift2 + Lnorbovr
!!$        end do
!!$        indspin = indspin + Lnorbovr * Lnorbovr
!!$        lspin = lspin + Localnorb(ilr) * Localnorb(ilr2)
!!$     end do
!!$
!!$     !deallocate the arrays depending on ilr2 (if they are allocated)
!!$     if(allocated(OnSitehamovr)) then
!!$        i_all=-product(shape(OnSitehamovr))*kind(OnSitehamovr)
!!$        deallocate(OnSitehamovr,stat=i_stat)
!!$        call memocc(i_stat,i_all,'OnSitehamovr',subname)
!!$     end if
!!$     if(allocated(Ahamovr)) then
!!$        i_all=-product(shape(Ahamovr))*kind(Ahamovr)
!!$        deallocate(Ahamovr,stat=i_stat)
!!$        call memocc(i_stat,i_all,'Ahamovr',subname)
!!$        i_all=-product(shape(Lohamovr))*kind(Lohamovr)
!!$        deallocate(Lohamovr,stat=i_stat)
!!$        call memocc(i_stat,i_all,'Lohamovr',subname)
!!$     end if
!!$  end do
!!$
!!$END SUBROUTINE overlap_matrix_for_locreg

!!$subroutine overlap_matrix_between_locreg(ilr,ilr2,isovrlp,nspin,orbscToAtom,psidim1,psidim2,psishift1,psishift2,&
!!$           Lzd,orbs,Lpsi,Lhpsi,dim_Lhamovr,Lhamovr)
!!$
!!$  use module_base
!!$  use module_types
!!$
!!$ implicit none
!!$
!!$  ! Subroutine Scalar Arguments
!!$  integer,intent(in) :: ilr,ilr2                     ! index of current locregs
!!$  integer,intent(in) :: isovrlp                      ! number of overlaps(because of periodicity)
!!$  integer,intent(in) :: nspin                        ! number of spins
!!$  integer,intent(in) :: psidim1                      ! dimension of the wavefunctions in locreg(ilr)
!!$  integer,intent(in) :: psidim2                      ! dimension of the wavefunctions in locreg(ilr2)
!!$  integer,intent(in) :: psishift1                    ! starting index of the first orbital of locreg(ilr)
!!$  integer,intent(in) :: psishift2                    ! starting index of the first orbital of locreg(ilr2)
!!$  integer,intent(in) :: dim_Lhamovr                  ! dimension of the Local Hamiltonian/Overlap Matrix
!!$  type(local_zone_descriptors), intent(in) :: Lzd   ! Descriptors of regions for linear scaling
!!$  type(orbitals_data),intent(in) :: orbs
!!$
!!$  !Subroutine Array Arguments
!!$  integer , dimension(Lzd%nlr,nspin),intent(in) :: orbscToAtom
!!$  real(wp),dimension(Lzd%Lpsidimtot),intent(in) :: Lpsi       ! Wavefunction (compressed format)
!!$  real(wp),dimension(Lzd%Lpsidimtot),intent(in) :: Lhpsi      ! Wavefunction in localization region
!!$  real(wp),dimension(nspin*dim_Lhamovr,2,orbs%nkpts),intent(out) :: Lhamovr           ! Local Hamiltonian/Overlap matrix
!!$
!!$  ! Local Variables
!!$  integer :: iolr,i_all,i_stat
!!$  integer :: ldim,ldim1,ldim2
!!$  integer :: i1,i2,i3,ispin
!!$  integer :: norb1,norb2
!!$  integer :: spinshift,iorbst,iorbst2,imatrst
!!$  type(locreg_descriptors),dimension(:), allocatable :: Olr  ! Localization grid descriptors 
!!$  real(wp), dimension(:,:,:), allocatable :: Ahamovr  !hamiltonian/overlap matrices for intersection region
!!$  real(wp),dimension(:),pointer:: Lopsi1,Lopsi2,Lohpsi
!!$  character(len=*), parameter :: subname='overlap_matrix_for_locreg'
!!$  type(orbitals_data) :: orbsu      ! orbital descriptor   
!!$
!!$  orbsu = orbs
!!$  norb1 = Lzd%Llr(ilr)%Localnorb / orbs%nspin
!!$  norb2 = Lzd%LLr(ilr2)%Localnorb / orbs%nspin
!!$ 
!!$  !Calculate the number of overlap regions between two logregs (can be more then one because
!!$  ! of periodicity). The number of overlap regions is stored in the isvorlp integer.
!!$  if (ilr .ne. ilr2) then
!!$     !allocate the overlap regions (memocc not used because it is a type) 
!!$    allocate(Olr(isovrlp),stat=i_stat)
!!$ 
!!$    ! Now allocate the total overlap matrix for the overlap region [Lohamovr]
!!$    ! and a temporary overlap matrix [Ahamovr] for each zone (Olr fractured because of periodicity)
!!$    ! It has dimensions: Lohamovr(nspin*norb**2,2,nkpt). The second dimension is:
!!$    ! 1: Hamiltonian Matrix =  <Obr(iorb1)|H|Orb(iorb2)>
!!$    ! 2: Overlap Matrix = <Obr(iorb1)|Orb(iorb2)>
!!$    allocate(Ahamovr(nspin*dim_Lhamovr,2,orbs%nkpts+ndebug),stat=i_stat)    !temporary
!!$    call memocc(i_stat,Ahamovr,'Ahamovr',subname)
!!$ 
!!$    !since Lhamovr is an accumulator
!!$    call razero(nspin*dim_Lhamovr*2*(orbs%nkpts+ndebug),Lhamovr)
!!$
!!$    ! Construct the overlap region descriptors
!!$    call get_overlap_region_periodic(ilr,ilr2,Lzd%Glr,isovrlp,Lzd%Llr,Lzd%nlr,Olr)
!!$ 
!!$    ! Third, transform the wavefunction to overlap regions (assuming same for all spin)
!!$    do iolr=1,isovrlp
!!$       ldim  =  Olr(iolr)%wfd%nvctr_c+7*Olr(iolr)%wfd%nvctr_f
!!$       ldim1 = ldim * (norb1-orbscToAtom(ilr,1)) * orbs%nspinor
!!$       ldim2 = ldim * (norb2-orbscToAtom(ilr2,1))* orbs%nspinor
!!$ 
!!$       ! Allocate the local wavefunction (in one overlap region)
!!$       allocate(Lopsi1(ldim1*nspin+ndebug), stat=i_stat)
!!$       call memocc(i_stat,Lopsi1,'Lopsi1',subname)
!!$       allocate(Lopsi2(ldim2*nspin+ndebug), stat=i_stat)
!!$       call memocc(i_stat,Lopsi2,'Lopsi2',subname)
!!$       allocate(Lohpsi(ldim2*nspin+ndebug), stat=i_stat)
!!$       call memocc(i_stat,Lohpsi,'Lohpsi',subname)
!!$
!!$       ! Project the wavefunctions inside the overlap region (first for Llr(ilr)and second for Llr(ilr2))
!!$       orbsu%npsidim = psidim1
!!$       do ispin = 1, nspin
!!$          orbsu%norbp = norb1-orbscToAtom(ilr,ispin)
!!$          call psi_to_locreg(Lzd%Llr(ilr),iolr,ldim1,Olr(iolr),&
!!$                    Lopsi1(1+(ispin-1)*ldim1:ispin*ldim1),isovrlp,orbs,&
!!$                    Lpsi(psishift1+(ispin-1)*psidim1:psishift1+ispin*psidim1-1))
!!$!          call psi_to_locreg(Llr(ilr),iolr,ldim1,Olr(iolr),&
!!$!                   Lohpsi(1+(ispin-1)*(ldim1+ldim2):ispin*ldim1+(ispin-1)*ldim2),isovrlp,orbs,&
!!$!                   Lhpsi(psishift1+(ispin-1)*psidim1:psishift1+ispin*psidim1-1))
!!$       end do
!!$       !second region
!!$       orbsu%npsidim = psidim2
!!$       do ispin = 1, nspin
!!$          orbsu%norbp = norb2-orbscToAtom(ilr2,ispin)
!!$          call psi_to_locreg(Lzd%Llr(ilr2),iolr,ldim2,Olr(iolr),&
!!$                   Lopsi2(1+(ispin-1)*ldim2:ispin*ldim2),&
!!$                   isovrlp,orbs,Lpsi(psishift2+(ispin-1)*psidim2:psishift2+ispin*psidim2-1))
!!$          call psi_to_locreg(Lzd%Llr(ilr2),iolr,ldim2,Olr(iolr),&
!!$                   Lohpsi(1+(ispin-1)*ldim2:ispin*ldim2),&
!!$                   isovrlp,orbs,Lhpsi(psishift2+(ispin-1)*psidim2:psishift2+ispin*psidim2-1))
!!$       end do
!!$
!!$      ! Now we can build the overlap matrix of the intersection region
!!$!      ispsi=1                                                     !TO DO: K-points and parallelization
!!$!      do ikptp=1,orbs%nkptsp
!!$!         ikpt=orbs%iskpts+ikptp!orbs%ikptsp(ikptp)
!!$!      
!!$!         nvctrp=commu%nvctr_par(iproc,ikptp)
!!$!         if (nvctrp == 0) cycle
!!$       iorbst =1
!!$       iorbst2=1
!!$       imatrst=1
!!$       do ispin=1,nspin !this construct assumes that the semicore is identical for both the spins
!!$  
!!$          call local_overlap_matrices((norb1+norb2-orbscToAtom(ilr,1)-orbscToAtom(ilr2,1))*nspin,&
!!$             norb1-orbscToAtom(ilr,1),norb2-orbscToAtom(ilr2,1),ldim,nspin,orbs%nspinor,dim_Lhamovr,&
!!$             Ahamovr(1,1,1),Lopsi1(1),Lopsi2,Lohpsi(1),iorbst,iorbst2,imatrst)  
!!$          iorbst =iorbst + norb1-orbscToAtom(ilr,1)
!!$          iorbst2=iorbst2+ norb2-orbscToAtom(ilr2,1)
!!$       end do
!!$
!!$!      ispsi=ispsi+nvctrp*norbtot*orbs%nspinor
!!$!      end do
!!$
!!$       ! The overlap matrix should be accumulated for each intersection region 
!!$       do i1=1,size(Ahamovr,1)
!!$          do i2=1,size(Ahamovr,2)
!!$             do i3=1,size(Ahamovr,3)
!!$                Lhamovr(i1,i2,i3) = Lhamovr(i1,i2,i3) + Ahamovr(i1,i2,i3) !Careful how we accumulate (orbital ordering)                
!!$              end do
!!$          end do
!!$       end do
!!$
!!$       ! deallocate the arrays depending on Olr
!!$       i_all=-product(shape(Lopsi1))*kind(Lopsi1)
!!$       deallocate(Lopsi1,stat=i_stat)
!!$       call memocc(i_stat,i_all,'Lopsi1',subname)
!!$
!!$       i_all=-product(shape(Lopsi2))*kind(Lopsi2)
!!$       deallocate(Lopsi2,stat=i_stat)
!!$       call memocc(i_stat,i_all,'Lopsi2',subname)
!!$
!!$       i_all=-product(shape(Lohpsi))*kind(Lohpsi)
!!$       deallocate(Lohpsi,stat=i_stat)
!!$       call memocc(i_stat,i_all,'Lohpsi',subname)
!!$    end do
!!$
!!$    ! deallocate Olr
!!$    deallocate(Olr,stat=i_stat)
!!$
!!$  else if (ilr == ilr2) then
!!$     ldim1 = (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f) * orbs%nspinor
!!$
!!$     imatrst=1
!!$     spinshift = 0
!!$     do ispin=1,nspin !this construct assumes that the semicore is identical for both the spins
!!$        call local_overlap_matrices((norb1-orbscToAtom(ilr,1)),&
!!$          norb1-orbscToAtom(ilr,1),norb2-orbscToAtom(ilr2,1),ldim1,&
!!$          nspin,orbs%nspinor,dim_Lhamovr,Lhamovr(1,1,1),Lpsi(psishift1+spinshift:psishift1+spinshift+psidim1-1),&
!!$          Lpsi(psishift2+spinshift:psishift2+spinshift+psidim2-1),&
!!$          Lhpsi(psishift2+spinshift:psishift2+spinshift+psidim2-1),1,1,imatrst)
!!$        spinshift = spinshift + (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f) * orbs%nspinor * norb1 
!!$     end do
!!$  end if
!!$
!!$  !deallocate the arrays depending on ilr2 (if they are allocated)
!!$  if(allocated(Ahamovr)) then
!!$     i_all=-product(shape(Ahamovr))*kind(Ahamovr)
!!$     deallocate(Ahamovr,stat=i_stat)
!!$     call memocc(i_stat,i_all,'Ahamovr',subname)
!!$  end if
!!$
!!$END SUBROUTINE overlap_matrix_between_locreg

!!$subroutine semicore_overlap_matrices(ilr,nspin,nspinor,norbtot,Lzd,orbs,orbscToAtom,ndim_hamovr,hamovr,scstr,psi,hpsi)
!!$  use module_base
!!$  use module_types
!!$  implicit none
!!$  integer, intent(in) :: ilr         ! localization region
!!$  integer, intent(in) :: ndim_hamovr,nspinor,nspin,scstr,norbtot
!!$  type(local_zone_descriptors),intent(in) :: Lzd
!!$  type(orbitals_data) :: orbs
!!$  integer, dimension(Lzd%nlr,nspin),intent(in) :: orbscToAtom 
!!$  real(wp), dimension(nspin*ndim_hamovr,2,orbs%nkpts+ndebug), intent(inout) :: hamovr
!!$  real(wp), dimension(orbs%npsidim), intent(in) :: psi,hpsi
!!$  !local variables
!!$  integer :: ncomp,ncplx,pos,ii,ispsi,ikptp,ikpt,i_all
!!$  integer :: ispin,norbe,nvctrp,ndim_hamsc,start,iel
!!$  integer :: i_stat,jj,kk
!!$  real(wp),allocatable :: psitmp(:,:),hpsitmp(:,:)
!!$  character(len=*), parameter :: subname='semicore_overlap_matrices'
!!$  real(wp), allocatable :: hamsc(:,:)
!!$  !WARNING: here nspin=1 for nspinor=4
!!$  if(nspinor == 1) then
!!$     ncplx=1
!!$  elseif(nspinor == 2) then
!!$     ncplx=2
!!$     ncomp=1
!!$  else if (nspinor == 4) then
!!$     ncplx=2
!!$     ncomp=2
!!$  end if
!!$
!!$!WARNING: What about the kpoints? Does this work? Not parallel, because nvctrp not distributed
!!$  ispsi=0
!!$  do ikptp=1,orbs%nkptsp
!!$     ikpt = orbs%iskpts+ikptp!orbsu%ikptsp(ikptp)
!!$
!!$     do ispin=1,nspin
!!$        norbe = orbscToAtom(ilr,ispin)
!!$        nvctrp = Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f
!!$        ndim_hamsc = norbe**2
!!$        start  = scstr + (ispin-1)*nvctrp*Lzd%Llr(ilr)%Localnorb/nspin + ispsi
!!$        iel = start-1 + nvctrp
!!$
!!$        allocate(hamsc(ndim_hamsc,2),stat=i_stat)
!!$        call memocc(i_stat,hamsc,'hamsc',subname)
!!$        allocate(psitmp(nvctrp*nspinor,norbe))
!!$        call memocc(i_stat,psitmp,'psitmp',subname)
!!$        allocate(hpsitmp(nvctrp*nspinor,norbe))
!!$        call memocc(i_stat,hpsitmp,'hpsitmp',subname)
!!$
!!$        kk=0 
!!$        do jj=1,norbe
!!$           do ii=1,nvctrp*nspinor
!!$              psitmp(ii,jj)=psi(start+kk)
!!$              hpsitmp(ii,jj)=hpsi(start+kk)
!!$              kk=kk+1
!!$           end do
!!$        end do
!!$
!!$        !calculate the overlap matrix for each group of the semicore atoms
!!$        !       hamovr(jorb,iorb,1)=+psit(k,jorb)*hpsit(k,iorb)
!!$        !       hamovr(jorb,iorb,2)=+psit(k,jorb)* psit(k,iorb)
!!$        if (nspinor ==1) then
!!$           call gemm('T','N',norbe,norbe,nvctrp,1.0_wp,psitmp(1,1),max(1,nvctrp),&
!!$                hpsitmp(1,1),max(1,nvctrp),0.0_wp,hamsc(1,1),norbe)
!!$           !here probably dsyrk can be used
!!$           call gemm('T','N',norbe,norbe,nvctrp,1.0_wp,psitmp(1,1),max(1,nvctrp),&
!!$                psitmp(1,1),max(1,nvctrp),0.0_wp,hamsc(1,2),norbe)
!!$        else
!!$           call c_gemm('C','N',norbe,norbe,ncomp*nvctrp,(1.0_wp,0.0_wp),psitmp(1,1),&
!!$                max(1,ncomp*nvctrp),hpsitmp(1,1),max(1,ncomp*nvctrp),&
!!$                (0.0_wp,0.0_wp),hamsc(1,1),norbe)
!!$           !here probably zherk can be used
!!$           call c_gemm('C','N',norbe,norbe,ncomp*nvctrp,(1.0_wp,0.0_wp),psitmp(1,1),&
!!$                max(1,ncomp*nvctrp),psitmp(1,1),max(1,ncomp*nvctrp),&
!!$                (0.0_wp,0.0_wp),hamsc(1,2),norbe)
!!$        end if
!!$
!!$        ! Calculate the starting position of hamsc in hamovr   
!!$        pos = 0
!!$        do ii=1,ilr-1
!!$           pos = pos + orbscToAtom(ii,ispin)**2
!!$        end do
!!$
!!$        ! Put into the correct order in hamovr
!!$        do ii=1,ndim_hamsc
!!$           hamovr(ii+pos+(ispin-1)*ndim_hamovr,:,ikpt) = hamsc(ii,:)
!!$        end do
!!$
!!$        ! deallocations
!!$        i_all = -product(shape(hamsc))*kind(hamsc)
!!$        deallocate(hamsc,stat=i_stat)
!!$        call memocc(i_stat,i_all,'hamsc',subname)
!!$        i_all = -product(shape(psitmp))*kind(psitmp)
!!$        deallocate(psitmp,stat=i_stat)
!!$        call memocc(i_stat,i_all,'psitmp',subname)
!!$        i_all = -product(shape(hpsitmp))*kind(hpsitmp)
!!$        deallocate(hpsitmp,stat=i_stat)
!!$        call memocc(i_stat,i_all,'hpsitmp',subname)
!!$
!!$
!!$     end do
!!$     ispsi=ispsi+nvctrp*norbtot*nspinor
!!$  end do
!!$
!!$END SUBROUTINE semicore_overlap_matrices

!> Tranform wavefunction between localisation region and the global region. 
!! Only assigns the indices for this transformation:
!! indexLpsi(i)=j means that element i of the localization region muss be copied to position j in the global region.
!!
!! @warning 
!!     Psi must be initialized to zero before entering this routine. 
!!     Each Lpsi is added to the corresponding place in Global.
!!     Only coded for sequential, not parallel cases !! For parallel should change increment and loc_psi dimensions
!subroutine Lpsi_to_global2(Glr,Gdim,Llr,lpsi,Ldim,norb,nspinor,nspin,shift,psi)
!!subroutine index_of_Lpsi_to_global2(iproc, nproc, ldim, gdim, norb, nspinor, nspin, Glr, Llr, indexLpsi)
!!
!!  use module_base
!!  use module_types
!!
!! implicit none
!!
!!  ! Subroutine Scalar Arguments
!!  integer,intent(in):: iproc, nproc
!!  integer :: Gdim          ! dimension of psi 
!!  integer :: Ldim          ! dimension of lpsi
!!  integer :: norb          ! number of orbitals
!!  integer :: nspinor       ! number of spinors
!!  integer :: nspin         ! number of spins 
!!  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
!!  type(locreg_descriptors), intent(in) :: Llr  ! Localization grid descriptors 
!!  
!!  !Subroutine Array Arguments
!!  integer,dimension(Ldim),intent(out) :: indexLpsi         !Wavefunction in localization region
!!  
!!  !local variables
!!  integer :: igrid,isegloc,isegG,ix,iorbs
!!  integer :: lmin,lmax,Gmin,Gmax
!!  integer :: icheck      ! check to make sure the dimension of loc_psi does not overflow 
!!  integer :: offset      ! gives the difference between the starting point of Lseg and Gseg
!!  integer :: length      ! Length of the overlap between Lseg and Gseg
!!  integer :: lincrement  ! Increment for writing orbitals in loc_psi
!!  integer :: Gincrement  ! Increment for reading orbitals in psi
!!  integer :: nseg        ! total number of segments in Llr
!!  integer, allocatable :: keymask(:,:)  ! shift for every segment of Llr (with respect to Glr)
!!  character(len=*), parameter :: subname='index_of_Lpsi_to_global'
!!  integer :: i_stat,i_all
!!  integer :: start,Gstart,Lindex
!!  integer :: lfinc,Gfinc,spinshift,ispin,Gindex,isegstart
!!  ! debug
!!  integer:: lxs, lys, lzs, lxe, lye, lze, gxe, gye, gze, locallength, actuallength
!!
!!! Define integers
!!  nseg = Llr%wfd%nseg_c + Llr%wfd%nseg_f
!!  lincrement = Llr%wfd%nvctr_c + 7*Llr%wfd%nvctr_f
!!  Gincrement = Glr%wfd%nvctr_c + 7*Glr%wfd%nvctr_f
!!  icheck = 0
!!  spinshift = Gdim / nspin
!! 
!!! Get the keymask: shift for every segment of Llr (with respect to Glr)
!!  allocate(keymask(2,nseg),stat=i_stat)
!!  call memocc(i_stat,keymask,'keymask',subname)
!!
!!  call shift_locreg_indexes(Glr,Llr,keymask,nseg)
!!!####################################################
!!! Do coarse region
!!!####################################################
!!  isegstart=1
!!  local_loop_c: do isegloc = 1,Llr%wfd%nseg_c
!!     lmin = keymask(1,isegloc)
!!     lmax = keymask(2,isegloc)
!!
!!     locallength=lmax-lmin+1
!!     actuallength=0
!!     ! Could optimize the routine by looping only on Gsegs not looped on before (TO DO)... DONE!
!!     !do isegG = 1,Glr%wfd%nseg_c
!!     global_loop_c: do isegG = isegstart,Glr%wfd%nseg_c
!!        Gmin = Glr%wfd%keygloc(1,isegG)
!!        Gmax = Glr%wfd%keygloc(2,isegG)
!!
!!        ! For each segment in Llr check if there is a collision with the segment in Glr
!!        ! if not, cycle
!!        if(lmin > Gmax) then
!!            isegstart=isegG
!!        end if
!!        if(Gmin > lmax) cycle local_loop_c
!!        if((lmin > Gmax) .or. (lmax < Gmin)) cycle global_loop_c
!!        !!if(iproc==0) write(*,'(a,4i9)') 'iproc, isegloc, isegstart, isegG', iproc, isegloc, isegstart, isegG
!!        !!if(isegG<isegstart) write(*,*) 'ERROR: isegG, isegstart', isegG, isegstart
!!
!!        ! Define the offset between the two segments
!!        offset = lmin - Gmin
!!        if(offset < 0) then
!!           offset = 0
!!        end if
!!
!!        ! Define the length of the two segments
!!        length = min(lmax,Gmax)-max(lmin,Gmin)
!!        actuallength=actuallength+length+1
!!
!!
!!        !Find the common elements and write them to the new global wavefunction
!!        ! WARNING: index goes from 0 to length because it is the offset of the element
!!        do ix = 0,length
!!           icheck = icheck + 1
!!           ! loop over the orbitals
!!           do ispin=1,nspin
!!              Gindex = Glr%wfd%keyvloc(isegG)+offset+ix+spinshift*(ispin-1)
!!              Lindex = icheck+lincrement*norb*(ispin-1)
!!              !psi(Gindex) = psi(Gindex) + lpsi(Lindex)
!!              indexLpsi(Lindex)=Gindex
!!           end do
!!           !!do iorbs=1,norb*nspinor
!!           !!   do ispin=1,nspin
!!           !!      Gindex = Glr%wfd%keyv(isegG)+offset+ix+Gincrement*(iorbs-1)+shift+spinshift*(ispin-1)
!!           !!      Lindex = icheck+lincrement*(iorbs-1)+lincrement*norb*(ispin-1)
!!           !!      psi(Gindex) = psi(Gindex) + lpsi(Lindex)
!!           !!   end do
!!           !!end do
!!        end do
!!     end do global_loop_c
!!     if(locallength/=actuallength) then
!!         write(*,'(2(a,i0),a)') 'ERROR: locallength = ',locallength, ' /= ', actuallength, ' = actuallength' 
!!     end if
!!  end do local_loop_c
!!
!!! Check if the number of elements in loc_psi is valid
!!  if(icheck .ne. Llr%wfd%nvctr_c) then
!!    write(*,*)'There is an error in index_of_Lpsi_to_global2: number of coarse points used',icheck
!!    write(*,*)'is not equal to the number of coarse points in the region',Llr%wfd%nvctr_c
!!  end if
!!
!!!##############################################################
!!! Now do fine region
!!!##############################################################
!!
!!  icheck = 0
!!  start = Llr%wfd%nvctr_c
!!  Gstart = Glr%wfd%nvctr_c
!!  lfinc  = Llr%wfd%nvctr_f
!!  Gfinc = Glr%wfd%nvctr_f
!!
!!  isegstart=Glr%wfd%nseg_c+1
!!  local_loop_f: do isegloc = Llr%wfd%nseg_c+1,nseg
!!     lmin = keymask(1,isegloc)
!!     lmax = keymask(2,isegloc)
!!
!!     ! Could optimize the routine by looping only on Gsegs not looped on before (TO DO)... DONE!
!!     !do isegG = Glr%wfd%nseg_c+1,Glr%wfd%nseg_c+Glr%wfd%nseg_f
!!     global_loop_f: do isegG = isegstart,Glr%wfd%nseg_c+Glr%wfd%nseg_f
!!
!!        Gmin = Glr%wfd%keygloc(1,isegG)
!!        Gmax = Glr%wfd%keygloc(2,isegG)
!!
!!        ! For each segment in Llr check if there is a collision with the segment in Glr
!!        ! if not, cycle
!!        if(lmin > Gmax) then
!!            isegstart=isegG
!!        end if
!!        if(Gmin > lmax) cycle local_loop_f
!!        if((lmin > Gmax) .or. (lmax < Gmin)) cycle global_loop_f
!!
!!        offset = lmin - Gmin
!!        if(offset < 0) offset = 0
!!
!!        length = min(lmax,Gmax)-max(lmin,Gmin)
!!
!!        !Find the common elements and write them to the new global wavefunction
!!        ! WARNING: index goes from 0 to length because it is the offset of the element
!!        do ix = 0,length
!!           icheck = icheck + 1
!!           !!do igrid=0,6
!!           !!   do iorbs=1,norb*nspinor
!!           !!     do ispin = 1, nspin
!!           !!        Gindex = Gstart+Glr%wfd%keyv(isegG)+offset+ix+Gincrement*(iorbs-1)+igrid*Gfinc+&
!!           !!                 shift + spinshift*(ispin-1)
!!           !!        Lindex = start+icheck+lincrement*(iorbs-1)+igrid*lfinc + lincrement*norb*(ispin-1) 
!!           !!        psi(Gindex) = psi(Gindex) + lpsi(Lindex)
!!           !!     end do
!!           !!   end do
!!           !!end do
!!           do igrid=1,7
!!              do ispin = 1, nspin
!!                 Gindex = Gstart + (Glr%wfd%keyvloc(isegG)+offset+ix-1)*7+igrid + spinshift*(ispin-1)
!!                 Lindex = start+(icheck-1)*7+igrid + lincrement*norb*(ispin-1) 
!!                 !psi(Gindex) = psi(Gindex) + lpsi(Lindex)
!!                 indexLpsi(Lindex)=Gindex
!!              end do
!!           end do
!!        end do
!!     end do global_loop_f
!!  end do local_loop_f
!!
!! ! Check if the number of elements in loc_psi is valid
!!  if(icheck .ne. Llr%wfd%nvctr_f) then
!!    write(*,*)'There is an error in index_of_Lpsi_to_global: number of fine points used',icheck
!!    write(*,*)'is not equal to the number of fine points in the region',Llr%wfd%nvctr_f
!!  end if
!!
!!  i_all=-product(shape(keymask))*kind(keymask)
!!  deallocate(keymask,stat=i_stat)
!!  call memocc(i_stat,i_all,'keymask',subname)
!!
!!END SUBROUTINE index_of_Lpsi_to_global2


!> Tranform wavefunction between global region and the localization region. 
!! Only assigns the indices for this transformation:
!! indexLpsi(i)=j means that element i of the localization region must be copied 
!! from position j in the global region.
!!
!! @warning 
!!     Only coded for sequential, not parallel cases 
!!     For parallel should change increment and loc_psi dimensions
!!$subroutine index_of_psi_to_locreg2(iproc, nproc, ldim, gdim, Llr, Glr, indexLpsi)
!!$
!!$  use module_base
!!$  use module_types
!!$ 
!!$ implicit none
!!$
!!$  ! Subroutine Scalar Arguments
!!$  integer,intent(in) :: iproc                  ! process ID
!!$  integer,intent(in) :: nproc                  ! number of processes
!!$  integer,intent(in) :: ldim          ! dimension of lpsi 
!!$  integer,intent(in) :: gdim          ! dimension of gpsi 
!!$  type(locreg_descriptors),intent(in) :: Llr  ! Local grid descriptor
!!$  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
!!$  
!!$  !Subroutine Array Arguments
!!$  integer,dimension(ldim),intent(out) :: indexLpsi   !Wavefunction in localization region
!!$  
!!$  !local variables
!!$  integer :: igrid,isegloc,isegG,ix,iorbs
!!$  integer :: lmin,lmax,Gmin,Gmax
!!$  integer :: icheck      ! check to make sure the dimension of loc_psi does not overflow 
!!$  integer :: offset      ! gives the difference between the starting point of Lseg and Gseg
!!$  integer :: length      ! Length of the overlap between Lseg and Gseg
!!$  integer :: lincrement  ! Increment for writing orbitals in loc_psi
!!$  integer :: Gincrement  ! Increment for reading orbitals in psi
!!$  integer :: nseg        ! total number of segments in Llr
!!$  integer, allocatable :: keymask(:,:)  ! shift for every segment of Llr (with respect to Glr)
!!$  character(len=*), parameter :: subname='psi_to_locreg'
!!$  integer :: i_stat,i_all
!!$  integer :: start,Gstart
!!$  integer :: lfinc,Gfinc,isegstart
!!$
!!$! Define integers
!!$  nseg = Llr%wfd%nseg_c + Llr%wfd%nseg_f
!!$  lincrement = Llr%wfd%nvctr_c + 7*Llr%wfd%nvctr_f
!!$  Gincrement = Glr%wfd%nvctr_c + 7*Glr%wfd%nvctr_f
!!$  icheck = 0
!!$
!!$! Initialize loc_psi
!!$  !call razero(ldim, lpsi)
!!$
!!$! Get the keymask: shift for every segment of Llr (with respect to Glr)
!!$  allocate(keymask(2,nseg),stat=i_stat)
!!$  call memocc(i_stat,keymask,'keymask',subname)
!!$
!!$  call shift_locreg_indexes(Glr,Llr,keymask,nseg)
!!$
!!$!####################################################
!!$! Do coarse region
!!$!####################################################
!!$  isegstart=1
!!$  local_loop_c: do isegloc = 1,Llr%wfd%nseg_c
!!$     lmin = keymask(1,isegloc)
!!$     lmax = keymask(2,isegloc)
!!$ 
!!$     ! Could optimize the routine by looping only on Gsegs not looped on before (TO DO).. DONE
!!$     global_loop_c: do isegG = isegstart,Glr%wfd%nseg_c
!!$        Gmin = Glr%wfd%keygloc(1,isegG)
!!$        Gmax = Glr%wfd%keygloc(2,isegG)
!!$
!!$        ! For each segment in Llr check if there is a collision with the segment in Glr
!!$        ! if not, cycle
!!$        if(lmin > Gmax) then
!!$            isegstart=isegG
!!$        end if
!!$        if(Gmin > lmax) cycle local_loop_c
!!$        if((lmin > Gmax) .or. (lmax < Gmin)) cycle global_loop_c
!!$        
!!$        ! Define the offset between the two segments
!!$        offset = lmin - Gmin
!!$        if(offset < 0) then
!!$           offset = 0
!!$        end if
!!$    
!!$        ! Define the length of the two segments
!!$        length = min(lmax,Gmax)-max(lmin,Gmin)
!!$ 
!!$        !Find the common elements and write them to the new localized wavefunction
!!$        ! WARNING: index goes from 0 to length because it is the offset of the element
!!$        do ix = 0,length
!!$           icheck = icheck + 1
!!$           !lpsi(icheck) = gpsi(Glr%wfd%keyv(isegG)+offset+ix)
!!$           indexLpsi(icheck)=Glr%wfd%keyvloc(isegG)+offset+ix
!!$           !!! loop over the orbitals
!!$           !!do iorbs=1,orbs%norbp*orbs%nspinor
!!$           !!   lpsi(icheck+lincrement*(iorbs-1))=psi(Glr%wfd%keyv(isegG)+offset+ix+Gincrement*(iorbs-1))
!!$           !!end do
!!$        end do
!!$     end do global_loop_c
!!$  end do local_loop_c
!!$
!!$! Check if the number of elements in loc_psi is valid
!!$  if(icheck .ne. Llr%wfd%nvctr_c) then
!!$    write(*,'(2(a,i0))')'There is an error in index_of_psi_to_locreg2: number of coarse points used, ',icheck, &
!!$              ', is not equal to the number of coarse points in the region, ',Llr%wfd%nvctr_c
!!$  end if
!!$
!!$!##############################################################
!!$! Now do fine region
!!$!##############################################################
!!$
!!$  icheck = 0
!!$  start = Llr%wfd%nvctr_c
!!$  Gstart = Glr%wfd%nvctr_c
!!$  lfinc  = Llr%wfd%nvctr_f
!!$  Gfinc = Glr%wfd%nvctr_f
!!$
!!$  isegstart=Glr%wfd%nseg_c+1
!!$  local_loop_f: do isegloc = Llr%wfd%nseg_c+1,nseg
!!$     lmin = keymask(1,isegloc)
!!$     lmax = keymask(2,isegloc)
!!$ 
!!$! Could optimize the routine by looping only on Gsegs not looped on before (TO DO)
!!$     global_loop_f: do isegG = isegstart,Glr%wfd%nseg_c+Glr%wfd%nseg_f
!!$
!!$        Gmin = Glr%wfd%keygloc(1,isegG)
!!$        Gmax = Glr%wfd%keygloc(2,isegG)
!!$
!!$        ! For each segment in Llr check if there is a collision with the segment in Glr
!!$        ! if not, cycle
!!$        if(lmin > Gmax) then
!!$            isegstart=isegG
!!$        end if
!!$        if(Gmin > lmax) cycle local_loop_f
!!$        if((lmin > Gmax) .or. (lmax < Gmin)) cycle global_loop_f
!!$
!!$        offset = lmin - Gmin
!!$        if(offset < 0) offset = 0
!!$
!!$        length = min(lmax,Gmax)-max(lmin,Gmin)
!!$
!!$        !Find the common elements and write them to the new localized wavefunction
!!$        ! WARNING: index goes from 0 to length because it is the offset of the element
!!$        do ix = 0,length
!!$           icheck = icheck + 1
!!$           do igrid=1,7
!!$              !lpsi(start+(icheck-1)*7+igrid) = gpsi(Gstart+(Glr%wfd%keyv(isegG)+offset+ix-1)*7+igrid)
!!$              indexLpsi(start+(icheck-1)*7+igrid) = Gstart+(Glr%wfd%keyvloc(isegG)+offset+ix-1)*7+igrid
!!$              !lpsi(start+(icheck-1)*7+igrid) = gpsi(Gstart+(Glr%wfd%keyv(isegG)+ix-1)*7+offset+igrid)
!!$              !!do iorbs=1,orbs%norbp*orbs%nspinor
!!$              !!   lpsi(start+icheck+lincrement*(iorbs-1)+igrid*lfinc)=&
!!$              !!   psi(Gstart+Glr%wfd%keyv(isegG)+offset+ix+Gincrement*(iorbs-1)+igrid*Gfinc)
!!$              !!end do
!!$           end do
!!$        end do
!!$     end do global_loop_f
!!$  end do local_loop_f
!!$  
!!$ ! Check if the number of elements in loc_psi is valid
!!$  if(icheck .ne. Llr%wfd%nvctr_f) then
!!$    write(*,*)'There is an error in psi_to_locreg: number of fine points used',icheck
!!$    write(*,*)'is not equal to the number of fine points in the region',Llr%wfd%nvctr_f
!!$  end if
!!$
!!$  i_all=-product(shape(keymask))*kind(keymask)
!!$  deallocate(keymask,stat=i_stat)
!!$  call memocc(i_stat,i_all,'keymask',subname)
!!$
!!$END SUBROUTINE index_of_psi_to_locreg2

!!subroutine local_overlap_matrices(norbe,norb1,norb2,nvctrp,nspin,nspinor,ndim_hamovr,hamovr,psi,psi2,hpsi,&
!!                                  iorbst,iorbst2,imatrst)
!!  use module_base
!!  implicit none
!!  integer, intent(in) :: norbe          ! total number of orbitals for overlap region
!!  integer, intent(in) :: norb1          ! number of orbitals in first locreg
!!  integer, intent(in) :: norb2          ! number of orbitals in second locreg
!!  integer, intent(in) :: nvctrp,ndim_hamovr,nspin,nspinor
!!  integer, intent(in) :: iorbst,iorbst2
!!  integer, intent(inout) :: imatrst
!!  real(wp), dimension(nspin*ndim_hamovr,2), intent(out) :: hamovr
!!  real(wp), dimension(nvctrp*nspinor,norbe), intent(in) :: psi,psi2,hpsi
!!  !local variables
!!  integer :: ispin,ncomp,ncplx
!!  !WARNING: here nspin=1 for nspinor=4
!!  if(nspinor == 1) then
!!     ncplx=1
!!  elseif(nspinor == 2) then
!!     ncplx=2
!!     ncomp=1
!!  else if (nspinor == 4) then
!!     ncplx=2
!!     ncomp=2
!!  end if
!!
!!     if (nspinor ==1) then
!!        call gemm('T','N',norb1,norb2,nvctrp,1.0_wp,psi(1,iorbst),max(1,nvctrp),&
!!             hpsi(1,iorbst2),max(1,nvctrp),&
!!             0.0_wp,hamovr(imatrst,1),norb1)
!!        !here probably dsyrk can be used
!!        call gemm('T','N',norb1,norb2,nvctrp,1.0_wp,psi(1,iorbst),max(1,nvctrp),&
!!             psi2(1,iorbst2),max(1,nvctrp),0.0_wp,hamovr(imatrst,2),norb1)
!!     else
!!        call c_gemm('C','N',norb1,norb2,ncomp*nvctrp,(1.0_wp,0.0_wp),psi(1,iorbst),&
!!             max(1,ncomp*nvctrp),hpsi(1,iorbst2),max(1,ncomp*nvctrp),&
!!             (0.0_wp,0.0_wp),hamovr(imatrst,1),norb1)
!!        !here probably zherk can be used
!!        call c_gemm('C','N',norb2,norb1,ncomp*nvctrp,(1.0_wp,0.0_wp),psi(1,iorbst),&
!!             max(1,ncomp*nvctrp),psi2(1,iorbst2),max(1,ncomp*nvctrp),&
!!             (0.0_wp,0.0_wp),hamovr(imatrst,2),norb1)
!!     end if
!!
!!     imatrst =imatrst+ncplx*norb1*norb2
!!
!!END SUBROUTINE local_overlap_matrices
