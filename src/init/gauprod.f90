!!****f* BigDFT/restart_from_gaussians
!! FUNCTION
!!  Restart from gaussian functions
!!
!! COPYRIGHT
!!    Copyright (C) 2007-2011 CEA (LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! SOURCE
!!
subroutine restart_from_gaussians(iproc,nproc,orbs,lr,hx,hy,hz,psi,G,coeffs)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,nproc
  real(gp), intent(in) :: hx,hy,hz
  type(orbitals_data), intent(in) :: orbs
  type(locreg_descriptors), intent(in) :: lr
  type(gaussian_basis), intent(inout) :: G
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%norbp*orbs%nspinor), intent(out) :: psi
  real(wp), dimension(:,:), pointer :: coeffs
  !local variables
  character(len=*), parameter :: subname='restart_from_gaussians'
  integer :: i_stat,i_all

  !the atomic positions corresponds to the new values
  !calculate the dual coefficients with the new positions
  
  !call gaussian_orthogonality(iproc,nproc,norb,norbp,G,coeffs)

  call dual_gaussian_coefficients(orbs%norbp,G,coeffs)

  !call gaussians_to_wavelets(iproc,nproc,lr%geocode,orbs,lr%d,hx,hy,hz,lr%wfd,G,coeffs,psi)
  call gaussians_to_wavelets_new(iproc,nproc,lr,orbs,hx,hy,hz,G,coeffs,psi)

  !deallocate gaussian structure and coefficients
  call deallocate_gwf(G,subname)
  i_all=-product(shape(coeffs))*kind(coeffs)
  deallocate(coeffs,stat=i_stat)
  call memocc(i_stat,i_all,'coeffs',subname)

  nullify(G%rxyz)

END SUBROUTINE restart_from_gaussians
!!***


!!****f* BigDFT/read_gaussian_information
!! FUNCTION
!!   Read information for gaussian basis set (from CP2K) or for restarting
!!
!! SOURCE
!!
subroutine read_gaussian_information(orbs,G,coeffs,filename, opt_fillrxyz)
  use module_base
  use module_types
  implicit none
  character(len=*), intent(in) :: filename
  type(orbitals_data), intent(inout) :: orbs
  type(gaussian_basis), intent(out) :: G
  real(wp), dimension(:,:),   pointer :: coeffs
  logical , optional :: opt_fillrxyz
  !local variables
  character(len=*), parameter :: subname='read_gaussian_information'
  logical :: exists
  integer :: i_stat,jexpo,iexpo,iat,iorb,jat,icoeff,jcoeff,jorb,j
  real(gp) :: rx,ry
  real(gp), dimension(4) :: coeff
  logical fillrxyz

  if (present(opt_fillrxyz)) then
     fillrxyz=opt_fillrxyz
  else
     fillrxyz=.false.
  endif

  !read the information from a file
  inquire(file=filename,exist=exists)
  if (.not. exists) then
     !if (iproc == 0) 
           write(*,'(1x,3a)')&
          'ERROR: The gaussian wavefunctions file "',trim(filename),'" is lacking, exiting...'
     stop
  end if

  open(unit=99,file=filename,status='unknown')
  read(99,*)G%nat,G%nshltot,G%nexpo,G%ncoeff

  allocate(G%nshell(G%nat+ndebug),stat=i_stat)
  call memocc(i_stat,G%nshell,'G%nshell',subname)
  allocate(G%nam(G%nshltot+ndebug),stat=i_stat)
  call memocc(i_stat,G%nam,'G%nam',subname)
  allocate(G%ndoc(G%nshltot+ndebug),stat=i_stat)
  call memocc(i_stat,G%ndoc,'G%ndoc',subname)
  allocate(G%xp(G%nexpo+ndebug),stat=i_stat)
  call memocc(i_stat,G%xp,'G%xp',subname)
  allocate(G%psiat(G%nexpo+ndebug),stat=i_stat)
  call memocc(i_stat,G%psiat,'G%psiat',subname)

  allocate(coeffs(G%ncoeff,orbs%norbp*orbs%nspinor+ndebug),stat=i_stat)
  call memocc(i_stat,coeffs,'coeffs',subname)
  

  if(fillrxyz) then

     allocate(G%rxyz (3,G%nat+ndebug),stat=i_stat)
     call memocc(i_stat,G%rxyz,'G%rxyz',subname)
     do iat=1,G%nat
        read(99,*)jat,G%rxyz(1, iat),G%rxyz(2, iat),G%rxyz(3, iat)  ,G%nshell(iat)
     end do
  else
     do iat=1,G%nat
        read(99,*)jat,rx,ry ,ry  ,G%nshell(iat)
     end do
  endif



  read(99,*)G%ndoc(1:G%nshltot),G%nam(1:G%nshltot)
  do iexpo=1,G%nexpo
     read(99,*)jexpo,G%xp(jexpo),G%psiat(jexpo)
  end do
  do iorb=1,orbs%norb
     read(99,*)jorb,orbs%eval(jorb)
     do icoeff=1,G%ncoeff
        read(99,*)jorb,jcoeff,(coeff(j),j=1,orbs%nspinor)
        if (orbs%isorb < iorb .and. iorb <= orbs%isorb+orbs%norbp) then
           do j=1,orbs%nspinor
              coeffs(jcoeff,(jorb-1-orbs%isorb)*orbs%nspinor+j)=coeff(j)
           end do
        end if
     end do
  end do
  close(99)
 
END SUBROUTINE read_gaussian_information
!!***


!!****f* BigDFT/write_gaussian_information
!! FUNCTION
!!   Write gaussian informatio for another program or for restarting
!!
!! SOURCE
!!
subroutine write_gaussian_information(iproc,nproc,orbs,G,coeffs,filename)
  use module_base
  use module_types
  implicit none
  character(len=*), intent(in) :: filename
  integer, intent(in) :: iproc,nproc
  type(gaussian_basis), intent(in) :: G
  type(orbitals_data), intent(in) :: orbs
  real(wp), dimension(G%ncoeff,orbs%norbp*orbs%nspinor), intent(in) :: coeffs
  !local variables
  character(len=*), parameter :: subname='write_gaussian_information'
  integer :: jproc,i_stat,i_all,ierr,iexpo,iat,iorb,icoeff,jorb,j,norb_tot
  integer, dimension(:,:), allocatable :: gatherarr
  real(gp), dimension(:,:), allocatable :: gaupsi

  allocate(gaupsi(G%ncoeff,orbs%norb*orbs%nspinor+ndebug),stat=i_stat)
  call memocc(i_stat,gaupsi,'gaupsi',subname)


  if (nproc > 1) then
     allocate(gatherarr(0:nproc-1,2+ndebug),stat=i_stat)
     call memocc(i_stat,gatherarr,'gatherarr',subname)
     
     norb_tot=0
     gatherarr(0,1)=G%ncoeff*orbs%norb_par(0)*orbs%nspinor
     gatherarr(0,2)=G%ncoeff*norb_tot*orbs%nspinor
     !gather the coefficients in a unique array
     do jproc=1,nproc-1
        norb_tot=norb_tot+orbs%norb_par(jproc-1)
        gatherarr(jproc,1)=G%ncoeff*orbs%norb_par(jproc)
        gatherarr(jproc,2)=G%ncoeff*norb_tot*orbs%nspinor
     end do

     call MPI_GATHERV(coeffs,gatherarr(iproc,1),mpidtypw,gaupsi,gatherarr(0,1),gatherarr(0,2),&
          mpidtypw,0,MPI_COMM_WORLD,ierr)

     i_all=-product(shape(gatherarr))*kind(gatherarr)
     deallocate(gatherarr,stat=i_stat)
     call memocc(i_stat,i_all,'gatherarr',subname)
  else
     gaupsi(1:G%ncoeff,1:orbs%norb*orbs%nspinor)=&
          coeffs(1:G%ncoeff,1:orbs%norb*orbs%nspinor)
  end if

  !write the information on a file
  if (iproc == 0) then
     open(unit=99,file=filename,status='unknown')

     write(99,'(4(i6))')G%nat,G%nshltot,G%nexpo,G%ncoeff
     do iat=1,G%nat
        write(99,'(i6,3(1x,1pe21.14),i6)')iat,(G%rxyz(j,iat),j=1,3),G%nshell(iat)
     end do
     write(99,*)G%ndoc,G%nam
     do iexpo=1,G%nexpo
        write(99,'(i6,2(1x,1pe21.14))')iexpo,G%xp(iexpo),G%psiat(iexpo)
     end do
     do iorb=1,orbs%norb
        write(99,'(i6,1x,1pe21.14)')iorb,orbs%eval(iorb)
        do icoeff=1,G%ncoeff
           write(99,'(2(i6),4(1x,1pe21.14))')iorb,icoeff,&
                (gaupsi(icoeff,orbs%nspinor*(iorb-1)+jorb),jorb=1,orbs%nspinor)
        end do
     end do
     close(99)
  end if

  i_all=-product(shape(gaupsi))*kind(gaupsi)
  deallocate(gaupsi,stat=i_stat)
  call memocc(i_stat,i_all,'gaupsi',subname)
  
END SUBROUTINE write_gaussian_information
!!***


!!****f* BigDFT/gaussian_pswf_basis
!! FUNCTION
!!   gaussian section
!!   Create gaussian structure from input guess pseudo wavefunctions
!!
!! SOURCE
!!
subroutine gaussian_pswf_basis(ng,enlargerprb,iproc,nspin,at,rxyz,G,Gocc)
  use module_base
  use module_types
  implicit none
  logical, intent(in) :: enlargerprb
  integer, intent(in) :: iproc,nspin,ng
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%nat), target, intent(in) :: rxyz
  type(gaussian_basis), intent(out) :: G
  real(wp), dimension(:), pointer :: Gocc
  !local variables
  character(len=*), parameter :: subname='gaussian_pswf_basis'
  integer, parameter :: noccmax=2,lmax=4,nmax=6,nelecmax=32
  logical :: occeq
  integer :: i_stat,i_all,iat,ityp,ishell,iexpo,l,i,ig,ictotpsi,norbe,norbsc,ishltmp
  integer :: ityx,ntypesx,nspinor,jat,noncoll,icoeff,iocc,nlo,ispin,m,icoll
  real(gp) :: ek
  integer, dimension(lmax) :: nl
  real(gp), dimension(noccmax,lmax) :: occup
  logical, dimension(:,:,:), allocatable :: scorb
  integer, dimension(:), allocatable :: iatypex
  integer, dimension(:,:), allocatable :: norbsc_arr
  real(gp), dimension(:), allocatable :: psiatn,locrad
  real(gp), dimension(:,:), allocatable :: xpt
  real(gp), dimension(:,:,:), allocatable :: psiat  

  !quick return if possible
  !if the positions are already associated it means that the basis is generated
  if (associated(G%rxyz)) then
     return
  end if

  allocate(scorb(4,2,at%natsc+ndebug),stat=i_stat)
  call memocc(i_stat,scorb,'scorb',subname)
  allocate(norbsc_arr(at%natsc+1,1+ndebug),stat=i_stat)
  call memocc(i_stat,norbsc_arr,'norbsc_arr',subname)
  allocate(locrad(at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,locrad,'locrad',subname)

  !for the moment, only collinear
  nspinor=1
  !if non-collinear it is like nspin=1 but with the double of orbitals
  if (nspinor == 4) then
     noncoll=2
  else
     noncoll=1
  end if

  call readAtomicOrbitals(at,norbe,norbsc,nspin,nspinor,scorb,norbsc_arr,locrad)

  i_all=-product(shape(locrad))*kind(locrad)
  deallocate(locrad,stat=i_stat)
  call memocc(i_stat,i_all,'locrad',subname)

  !Generate the input guess via the inguess_generator
  !take also into account the IG polarisations

  !the number of gaussian centers are thus nat
  G%nat=at%nat
  !this pointing creates problems if at the next call the positions are given by a different array.
  !presumably the best if to copy the values and allocate the pointer
  G%rxyz => rxyz

  !copy the parsed values in the gaussian structure
  !count also the total number of shells
  allocate(G%nshell(at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,G%nshell,'G%nshell',subname)

  !calculate the number of atom types by taking into account the occupation
  allocate(iatypex(at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,iatypex,'iatypex',subname)
  
  ntypesx=0
  G%nshltot=0
  count_shells: do iat=1,at%nat
     ityp=at%iatype(iat)
     call count_atomic_shells(lmax,noccmax,nelecmax,nspin,nspinor,at%aocc(1,iat),occup,nl)
     G%nshell(iat)=(nl(1)+nl(2)+nl(3)+nl(4))
     G%nshltot=G%nshltot+G%nshell(iat)
     !check the occupation numbers and the atoms type
     !once you find something equal exit the procedure
     do jat=1,iat-1
        if (at%iatype(jat) == ityp) then
           occeq=.true.
           do i=1,nelecmax
              occeq = occeq .and. (at%aocc(i,jat) == at%aocc(i,iat))
           end do
           !have found another similar atoms
           if (occeq) then
              iatypex(iat)=iatypex(jat)
              cycle count_shells
           end if
        end if
     end do
     ntypesx=ntypesx+1
     iatypex(iat)=ntypesx
  end do count_shells

  allocate(G%ndoc(G%nshltot+ndebug),stat=i_stat)
  call memocc(i_stat,G%ndoc,'G%ndoc',subname)
  allocate(G%nam(G%nshltot+ndebug),stat=i_stat)
  call memocc(i_stat,G%nam,'G%nam',subname)


  !the default value for the gaussians is chosen to be 21
  allocate(xpt(ng,ntypesx+ndebug),stat=i_stat)
  call memocc(i_stat,xpt,'xpt',subname)
  allocate(psiat(ng,5,ntypesx+ndebug),stat=i_stat)
  call memocc(i_stat,psiat,'psiat',subname)
  allocate(psiatn(ng+ndebug),stat=i_stat)
  call memocc(i_stat,psiatn,'psiatn',subname)


  !assign shell IDs and count the number of exponents and coefficients
  G%nexpo=0
  G%ncoeff=0
  ishell=0
  ntypesx=0
  do iat=1,at%nat
     ityp=at%iatype(iat)
     ityx=iatypex(iat)
     ishltmp=0
     call count_atomic_shells(lmax,noccmax,nelecmax,nspin,nspinor,at%aocc(1,iat),occup,nl)
     if (ityx > ntypesx) then
        if (iproc == 0 .and. verbose > 1) then
           write(*,'(1x,a,a6,a)')&
                'Generation of input wavefunction data for atom ',&
                trim(at%atomnames(ityp)),&
                ':'
           call print_eleconf(nspin,nspinor,noccmax,nelecmax,lmax,&
                at%aocc(1,iat),at%iasctype(iat))
        end if

        call iguess_generator(at%nzatom(ityp),at%nelpsp(ityp),&
             real(at%nelpsp(ityp),gp),at%psppar(0,0,ityp),&
             at%npspcode(ityp),&
             ng-1,nl,5,noccmax,lmax,occup,xpt(1,ityx),&
             psiat(1,1,ityx),enlargerprb)
        ntypesx=ntypesx+1
        if (iproc == 0 .and. verbose > 1) write(*,'(1x,a)')'done.'
     end if

     do l=1,4
        do i=1,nl(l)
           ishell=ishell+1
           ishltmp=ishltmp+1
           G%ndoc(ishell)=ng!(ity)
           G%nam(ishell)=l
           G%nexpo=G%nexpo+ng!(ity)
           G%ncoeff=G%ncoeff+2*l-1
           !print *,'iat,i,l',iat,i,l,norbe,G%ncoeff
        end do
     end do
     if (ishltmp /= G%nshell(iat)) then
        write(*,*)'ERROR: ishelltmp <> nshell',ishell,G%nshell(iat)
        stop 
     end if
  end do

  !up to here the code of this section is identical to the 
  !atomic orbitals part in inputguess.f90

  !now we have to allocate the array of the "occupation numbers"
  !of the molecular orbitals
  allocate(Gocc(G%ncoeff+ndebug),stat=i_stat)
  call memocc(i_stat,Gocc,'Gocc',subname)
  call razero(G%ncoeff,Gocc)

  !allocate and assign the exponents and the coefficients
  allocate(G%psiat(G%nexpo+ndebug),stat=i_stat)
  call memocc(i_stat,G%psiat,'G%psiat',subname)
  allocate(G%xp(G%nexpo+ndebug),stat=i_stat)
  call memocc(i_stat,G%xp,'G%xp',subname)

  ishell=0
  iexpo=0
  icoeff=1
  do iat=1,at%nat
     ityp=at%iatype(iat)
     ityx=iatypex(iat)
     call count_atomic_shells(lmax,noccmax,nelecmax,nspin,nspinor,at%aocc(1,iat),occup,nl)
     ictotpsi=0
     iocc=0
     do l=1,4
        iocc=iocc+1
        nlo=nint(at%aocc(iocc,iat)) !just to increase the counting 
        do i=1,nl(l)
           ishell=ishell+1
           ictotpsi=ictotpsi+1
           call atomkin(l-1,ng,xpt(1,ityx),psiat(1,ictotpsi,ityx),psiatn,ek)
           do ig=1,G%ndoc(ishell)
              iexpo=iexpo+1
              G%psiat(iexpo)=psiatn(ig)
              G%xp(iexpo)=xpt(ig,ityp)
           end do

           do ispin=1,nspin
              do m=1,2*l-1
                 !each orbital has two electrons in the case of the 
                 !non-collinear case
                 do icoll=1,noncoll !non-trivial only for nspinor=4
                    iocc=iocc+1
                    Gocc(icoeff)=Gocc(icoeff)+at%aocc(iocc,iat)
                 end do
                 icoeff=icoeff+1
              end do
              icoeff=icoeff-(2*l-1)
           end do
           icoeff=icoeff+(2*l-1)
        end do
     end do
  end do
  if (iexpo /= G%nexpo) then
     write(*,*)'ERROR: iexpo <> nexpo',iexpo,G%nexpo
     stop 
  end if

  i_all=-product(shape(scorb))*kind(scorb)
  deallocate(scorb,stat=i_stat)
  call memocc(i_stat,i_all,'scorb',subname)
  i_all=-product(shape(xpt))*kind(xpt)
  deallocate(xpt,stat=i_stat)
  call memocc(i_stat,i_all,'xpt',subname)
  i_all=-product(shape(psiat))*kind(psiat)
  deallocate(psiat,stat=i_stat)
  call memocc(i_stat,i_all,'psiat',subname)
  i_all=-product(shape(norbsc_arr))*kind(norbsc_arr)
  deallocate(norbsc_arr,stat=i_stat)
  call memocc(i_stat,i_all,'norbsc_arr',subname)
  i_all=-product(shape(psiatn))*kind(psiatn)
  deallocate(psiatn,stat=i_stat)
  call memocc(i_stat,i_all,'psiatn',subname)
  i_all=-product(shape(iatypex))*kind(iatypex)
  deallocate(iatypex,stat=i_stat)
  call memocc(i_stat,i_all,'iatypex',subname)

END SUBROUTINE gaussian_pswf_basis
!!***


!!****f* BigDFT/gaussian_psp_basis
!! FUNCTION
!!   Extract the pseudopotential basis
!! WARNING
!!   This is not the complete PSP basis set. 
!!   The radial power term is lacking in the gaussian descriptors should be added if needed
!! SOURCE
!!
subroutine gaussian_psp_basis(at,rxyz,G)
  use module_base
  use module_types
  implicit none
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%nat), target, intent(in) :: rxyz
  type(gaussian_basis), intent(out) :: G  
  !local variables
  character(len=*), parameter :: subname='gaussian_psp_basis'
  integer :: iat,nshell,ityp,iexpo,l,ishell,i_stat

  G%nat=at%nat
  G%rxyz => rxyz
  allocate(G%nshell(at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,G%nshell,'G%nshell',subname)
 
  G%nshltot=0
  do iat=1,G%nat
     ityp=at%iatype(iat) 
     nshell=0
     do l=1,4 
        if (at%psppar(l,0,ityp) /= 0.0_gp) nshell=nshell+1
     enddo
     G%nshell(iat)=nshell
     G%nshltot=G%nshltot+nshell
  end do

  allocate(G%ndoc(G%nshltot+ndebug),stat=i_stat)
  call memocc(i_stat,G%ndoc,'G%ndoc',subname)
  allocate(G%nam(G%nshltot+ndebug),stat=i_stat)
  call memocc(i_stat,G%nam,'G%nam',subname)

  !assign shell IDs and count the number of exponents and coefficients
  G%nexpo=0
  G%ncoeff=0
  ishell=0
  do iat=1,G%nat
     ityp=at%iatype(iat)
     do l=1,4 
        if (at%psppar(l,0,ityp) /= 0.0_gp) then
           ishell=ishell+1
           G%ndoc(ishell)=1
           G%nam(ishell)=l
           G%nexpo=G%nexpo+1
           G%ncoeff=G%ncoeff+2*l-1
        end if
     enddo
  end do

  !allocate and assign the exponents and the coefficients
  allocate(G%xp(G%nexpo+ndebug),stat=i_stat)
  call memocc(i_stat,G%xp,'G%xp',subname)
  allocate(G%psiat(G%nexpo+ndebug),stat=i_stat)
  call memocc(i_stat,G%psiat,'G%psiat',subname)

  ishell=0
  iexpo=0
  do iat=1,G%nat
     ityp=at%iatype(iat)
     do l=1,4 
        if (at%psppar(l,0,ityp) /= 0.0_gp) then
           ishell=ishell+1
           iexpo=iexpo+1
           G%psiat(iexpo)=1.0_gp
           G%xp(iexpo)=at%psppar(l,0,ityp)
        end if
     end do
  end do

END SUBROUTINE gaussian_psp_basis
!!***


!!****f* BigDFT/gaussian_orthogonality
!! FUNCTION
!!
!! SOURCE
!!
subroutine gaussian_orthogonality(iproc,nproc,norb,norbp,G,coeffs)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,nproc,norb,norbp
  type(gaussian_basis), intent(in) :: G
  real(wp), dimension(G%ncoeff,norbp), intent(inout) :: coeffs
  !local variables
  character(len=*), parameter :: subname='gaussian_orthogonality' 
  integer :: iorb,i,jproc,i_stat,i_all,info,ierr
  integer, dimension(:,:), allocatable :: gatherarr
  real(gp), dimension(:,:), allocatable :: ovrlp,gaupsi,tmp,smat

  allocate(ovrlp(G%ncoeff,G%ncoeff+ndebug),stat=i_stat)
  call memocc(i_stat,ovrlp,'ovrlp',subname)
  allocate(gaupsi(G%ncoeff,norb+ndebug),stat=i_stat)
  call memocc(i_stat,gaupsi,'gaupsi',subname)
  allocate(tmp(G%ncoeff,norb+ndebug),stat=i_stat)
  call memocc(i_stat,tmp,'tmp',subname)
  allocate(smat(norb,norb+ndebug),stat=i_stat)
  call memocc(i_stat,smat,'smat',subname)


  if (nproc > 1) then
     allocate(gatherarr(0:nproc-1,2+ndebug),stat=i_stat)
     call memocc(i_stat,gatherarr,'gatherarr',subname)

     !gather the coefficients in a unique array
     do jproc=0,nproc-1
        gatherarr(jproc,1)=G%ncoeff*min(max(norb-jproc*norbp,0),norbp)
        gatherarr(jproc,2)=G%ncoeff*min(jproc*norbp,norb)
     end do

     call MPI_ALLGATHERV(coeffs,gatherarr(iproc,1),mpidtypw,gaupsi,gatherarr(0,1),gatherarr(0,2),&
          mpidtypw,MPI_COMM_WORLD,ierr)

     i_all=-product(shape(gatherarr))*kind(gatherarr)
     deallocate(gatherarr,stat=i_stat)
     call memocc(i_stat,i_all,'gatherarr',subname)
  else
     gaupsi(1:G%ncoeff,1:norb)=coeffs(1:G%ncoeff,1:norb)
  end if

  !overlap of the basis
  call gaussian_overlap(G,G,ovrlp)

  call dsymm('L','U',G%ncoeff,norb,1.0_gp,ovrlp(1,1),G%ncoeff,gaupsi(1,1),G%ncoeff,&
       0.0_gp,tmp(1,1),G%ncoeff)

  call gemm('T','N',norb,norb,G%ncoeff,1.0_gp,gaupsi(1,1),G%ncoeff,tmp(1,1),G%ncoeff,&
       0.0_wp,smat(1,1),norb)
  !print overlap matrices
  if (iproc==0) then
     do i=1,norb
        write(*,'(a,i3,i3,30(1pe19.12))')'ovrlp',iproc,i,(smat(i,iorb),iorb=1,norb)
     end do
  end if

  !orthogonalise the overlap matrix (from orthon_p)
  ! Cholesky factorization
  call potrf( 'L',norb,smat(1,1),norb,info)
  if (info.ne.0) write(6,*) 'info Cholesky factorization',info

  ! calculate L^{-1}
  call trtri( 'L','N',norb,smat(1,1),norb,info)
  if (info.ne.0) write(6,*) 'info L^-1',info

  ! new vectors   
  call trmm('R','L','T','N',G%ncoeff,norb,1.0_wp,smat(1,1),norb,gaupsi(1,1),G%ncoeff)

  !copy the result in the portion of the array
  do iorb=1,norbp
     if (iorb+iproc*norbp <= norb) then
        do i=1,G%ncoeff
           coeffs(i,iorb)=gaupsi(i,iorb+iproc*norbp)
        end do
     end if
  end do
  
  i_all=-product(shape(ovrlp))*kind(ovrlp)
  deallocate(ovrlp,stat=i_stat)
  call memocc(i_stat,i_all,'ovrlp',subname)
  i_all=-product(shape(tmp))*kind(tmp)
  deallocate(tmp,stat=i_stat)
  call memocc(i_stat,i_all,'tmp',subname)
  i_all=-product(shape(smat))*kind(smat)
  deallocate(smat,stat=i_stat)
  call memocc(i_stat,i_all,'smat',subname)
  i_all=-product(shape(gaupsi))*kind(gaupsi)
  deallocate(gaupsi,stat=i_stat)
  call memocc(i_stat,i_all,'gaupsi',subname)

END SUBROUTINE gaussian_orthogonality
!!***


!!****f* BigDFT/dual_gaussian_coefficients
!! FUNCTION
!!
!! SOURCE
!!
subroutine dual_gaussian_coefficients(norbp,G,coeffs)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: norbp
  type(gaussian_basis), intent(in) :: G
  real(gp), dimension(G%ncoeff,norbp), intent(inout) :: coeffs !warning: the precision here should be wp
  !local variables
  character(len=*), parameter :: subname='dual_gaussian_coefficients'
  integer :: nwork,info,i_stat,i_all
  integer, dimension(:), allocatable :: iwork
  real(gp), dimension(:), allocatable :: ovrlp,work
  
  allocate(iwork(G%ncoeff+ndebug),stat=i_stat)
  call memocc(i_stat,iwork,'iwork',subname)
  allocate(ovrlp(G%ncoeff*G%ncoeff+ndebug),stat=i_stat)
  call memocc(i_stat,ovrlp,'ovrlp',subname)

  !temporary allocation of the work array, workspace query in dsysv
  allocate(work(100+ndebug),stat=i_stat)
  call memocc(i_stat,work,'work',subname)

  if (norbp > 0) then
     call dsysv('U',G%ncoeff,norbp,ovrlp(1),G%ncoeff,iwork(1),coeffs(1,1),&
          G%ncoeff,work(1),-1,info)
  end if
  nwork=work(1)

  i_all=-product(shape(work))*kind(work)
  deallocate(work,stat=i_stat)
  call memocc(i_stat,i_all,'work',subname)
  allocate(work(nwork+ndebug),stat=i_stat)
  call memocc(i_stat,work,'work',subname)
  
  call gaussian_overlap(G,G,ovrlp)

!!  !overlap matrix, print it for convenience
!!  do icoeff=1,G%ncoeff
!!     write(30,'(1x,200(1pe12.3))')&
!!          (ovrlp(jcoeff+(icoeff-1)*G%ncoeff),jcoeff=1,G%ncoeff)
!!  end do

  if (norbp > 0) then
     call dsysv('U',G%ncoeff,norbp,ovrlp(1),G%ncoeff,iwork(1),coeffs(1,1),&
          G%ncoeff,work,nwork,info)
  end if

  i_all=-product(shape(iwork))*kind(iwork)
  deallocate(iwork,stat=i_stat)
  call memocc(i_stat,i_all,'iwork',subname)
  i_all=-product(shape(work))*kind(work)
  deallocate(work,stat=i_stat)
  call memocc(i_stat,i_all,'work',subname)
  i_all=-product(shape(ovrlp))*kind(ovrlp)
  deallocate(ovrlp,stat=i_stat)
  call memocc(i_stat,i_all,'ovrlp',subname)

!!!  do iorb=1,norbp
!!!     print *,'iorb, dual,coeffs',iorb,coeffs(:,iorb)
!!!  end do
  
END SUBROUTINE dual_gaussian_coefficients
!!***


!!****f* BigDFT/normalize_shell
!! FUNCTION
!!   Normalize a given atomic shell following the angular momentum
!! SOURCE
!!
subroutine normalize_shell(ng,l,expo,coeff)
  use module_base
  implicit none
  integer, intent(in) :: ng,l
  real(gp), dimension(ng), intent(in) :: expo
  real(gp), dimension(ng), intent(inout) :: coeff
  !local variables
  integer :: i,j
  real(gp) :: norm,tt,e1,ex,c1,c2,gauint0

  norm=0.d0
  do i=1,ng
     e1=expo(i)
     c1=coeff(i)
     do j=1,ng
        ex=expo(j)+e1
        c2=coeff(j)
        tt=gauint0(ex,2*l+2)
        norm=norm+c1*tt*c2
     end do
  end do
  norm=sqrt(0.5_gp*norm)
  norm=1.0_gp/norm
  do i=1,ng
     coeff(i)=coeff(i)*norm
  end do

  !print *,'l=',l,'norm=',norm

END SUBROUTINE normalize_shell
!!***


!!****f* BigDFT/gauinth
!! FUNCTION
!!   Calculates \int_0^\infty \exp^{-a*x^2} x^l dx
!! SOURCE
!!
function gauinth(a,l)
  use module_base
  implicit none
  integer, intent(in) :: l
  real(gp), intent(in) :: a
  real(gp) :: gauinth
  !local variables
  real(gp), parameter :: gammaonehalf=1.772453850905516027298d0
  integer :: p
  real(gp) :: xfac,prefac,tt,sh
  !build the prefactor
  prefac=sqrt(a)
  prefac=1.d0/prefac
  prefac=prefac**(l+1)
  p=1-l+2*(l/2)
  tt=0.5d0*gammaonehalf**p
  prefac=prefac*tt
  sh=-0.5d0*real(p,gp)
  p=l/2
  tt=xfac(1,p,sh)

  !final result
  gauinth=prefac*tt
  
END FUNCTION gauinth
!!***


!!****f* BigDFT/gprod
!! FUNCTION
!!   Calculates a dot product between two differents gaussians times spherical harmonics
!!   valid only for shell which belongs to different atoms, and with also dy/=0/=dx dz/=0
!!   to be rearranged when only some of them is zero
!! SOURCE
!!
subroutine gprod(a1,a2,dx,dy,dz,l1,m1,l2,m2,niw,nrw,iw,rw,ovrlp)
  use module_base
  implicit none
  integer, intent(in) :: l1,l2,m1,m2,niw,nrw 
  real(gp), intent(in) :: a1,a2,dx,dy,dz
  integer, dimension(niw) :: iw !work array of the exponents of the two polynomials
  real(gp), dimension(nrw) :: rw !work array of the polynomials coefficients 
  real(gp), intent(out) :: ovrlp
  !local variables
  integer, parameter :: nx=3
  integer :: n1,n2,i1,i2,px,py,pz,qx,qy,qz
  real(gp) :: fx,fy,fz,fa,fb,govrlp

  !calculates the number of different couples
  call calc_coeff_inguess(l1,m1,nx,n1,&
       iw(1),iw(nx+1),iw(2*nx+1),rw(1))
  call calc_coeff_inguess(l2,m2,nx,n2,&
       iw(3*nx+1),iw(4*nx+1),iw(5*nx+1),rw(n1+1))
  ovrlp=0.d0
  do i2=1,n2
     qx=iw(3*nx+i2)
     qy=iw(4*nx+i2)
     qz=iw(5*nx+i2)
     fb=rw(n1+i2)
     do i1=1,n1
        px=iw(i1)
        py=iw(nx+i1)
        pz=iw(2*nx+i1)
        fa=rw(i1)

        fx=govrlp(a1,a2,dx,px,qx)
        fy=govrlp(a1,a2,dy,py,qy)
        fz=govrlp(a1,a2,dz,pz,qz)

        ovrlp=ovrlp+fa*fb*fx*fy*fz
        !print *,i1,i2,fx,fy,fz,fa,fb
     end do
  end do
 
END SUBROUTINE gprod
!!***


!!****f* BigDFT/kinprod
!! FUNCTION
!!   Kinetic overlap between gaussians, based on cartesian coordinates
!!   calculates a dot product between two differents gaussians times spherical harmonics
!!   valid only for shell which belongs to different atoms, and with also dy/=0/=dx dz/=0
!!   to be rearranged when only some of them is zero
!!
!! SOURCE
!!
subroutine kinprod(a1,a2,dx,dy,dz,l1,m1,l2,m2,niw,nrw,iw,rw,ovrlp)
  use module_base
  implicit none
  integer, intent(in) :: l1,l2,m1,m2,niw,nrw 
  real(gp), intent(in) :: a1,a2,dx,dy,dz
  integer, dimension(niw) :: iw !work array of the exponents of the two polynomials
  real(gp), dimension(nrw) :: rw !work array of the polynomials coefficients 
  real(gp), intent(out) :: ovrlp
  !local variables
  integer, parameter :: nx=3
  integer :: n1,n2,i1,i2,px,py,pz,qx,qy,qz
  real(gp) :: fx,fy,fz,fa,fb,govrlp,kinovrlp,d2fx,d2fy,d2fz

  !calculates the number of different couples
  call calc_coeff_inguess(l1,m1,nx,n1,&
       iw(1),iw(nx+1),iw(2*nx+1),rw(1))
  call calc_coeff_inguess(l2,m2,nx,n2,&
       iw(3*nx+1),iw(4*nx+1),iw(5*nx+1),rw(n1+1))
  ovrlp=0.d0
  do i2=1,n2
     qx=iw(3*nx+i2)
     qy=iw(4*nx+i2)
     qz=iw(5*nx+i2)
     fb=rw(n1+i2)
     do i1=1,n1
        px=iw(i1)
        py=iw(nx+i1)
        pz=iw(2*nx+i1)
        fa=rw(i1)

        fx=govrlp(a1,a2,dx,px,qx)
        fy=govrlp(a1,a2,dy,py,qy)
        fz=govrlp(a1,a2,dz,pz,qz)

        d2fx=kinovrlp(a1,a2,dx,px,qx)
        d2fy=kinovrlp(a1,a2,dy,py,qy)
        d2fz=kinovrlp(a1,a2,dz,pz,qz)

        ovrlp=ovrlp-0.5_gp*fa*fb*(d2fx*fy*fz+fx*d2fy*fz+fx*fy*d2fz)
        !print *,i1,i2,fx,fy,fz,fa,fb
     end do
  end do
 
END SUBROUTINE kinprod
!!***


!!****f* BigDFT/kinovrlp
!! FUNCTION
!!   Calculates \int d^2/dx^2(\exp^{-a1*x^2} x^l1) \exp^{-a2*(x-d)^2} (x-d)^l2 dx
!!   in terms of the govrlp function below
!! SOURCE
!!
function kinovrlp(a1,a2,d,l1,l2)
  use module_base
  implicit none
  integer, intent(in) :: l1,l2
  real(gp), intent(in) :: a1,a2,d
  real(gp) :: kinovrlp
  !local variables
  real(gp) :: govrlp,fac,ovrlp

  !case l1+2
  fac=4._gp*a1**2
  ovrlp=govrlp(a1,a2,d,l1+2,l2)
  kinovrlp=fac*ovrlp
  !case l1
  fac=2._gp*a1*real(2*l1+1,gp)
  ovrlp=govrlp(a1,a2,d,l1,l2)
  kinovrlp=kinovrlp-fac*ovrlp
  !case l1-2 (if applicable)
  if (l1 >=2) then
     fac=real(l1*(l1-1),gp)
     ovrlp=govrlp(a1,a2,d,l1-2,l2)
     kinovrlp=kinovrlp+fac*ovrlp
  end if
END FUNCTION kinovrlp
!!***


!!****f* BigDFT/govrlp
!! FUNCTION
!!   Calculates \int \exp^{-a1*x^2} x^l1 \exp^{-a2*(x-d)^2} (x-d)^l2 dx
!! SOURCE
!!
function govrlp(a1,a2,d,l1,l2)
  use module_base
  implicit none
  integer, intent(in) :: l1,l2
  real(gp), intent(in) :: a1,a2,d
  real(gp) :: govrlp
  !local variables
  integer :: p
  real(gp) :: prefac,rfac,stot,aeff,ceff,tt,fsum,gauint

  !build the prefactor
  prefac=a1+a2
  prefac=a2/prefac
  prefac=a1*prefac
  prefac=-d**2*prefac
  prefac=dexp(prefac)

  !build the effective exponent and coefficients
  aeff=a1+a2
  ceff=a2*d
  ceff=ceff/aeff

  !build the first term in the sum
  stot=(-d)**l2
  stot=gauint(aeff,ceff,l1)*stot

  !perform the sum
  do p=1,l2/2
     tt=rfac(1,p)
     fsum=rfac(l2-p+1,l2)
     fsum=fsum/tt
     tt=(-d)**(l2-p)
     fsum=fsum*tt
     tt=gauint(aeff,ceff,l1+p)
     fsum=fsum*tt
     stot=stot+fsum
  end do
  do p=l2/2+1,l2
     tt=rfac(1,l2-p)
     fsum=rfac(p+1,l2)
     fsum=fsum/tt
     tt=(-d)**(l2-p)
     fsum=fsum*tt
     tt=gauint(aeff,ceff,l1+p)
     fsum=fsum*tt
     stot=stot+fsum
  end do

  !final result
  govrlp=prefac*stot
END FUNCTION govrlp
!!***


!!****f* BigDFT/gauint
!! FUNCTION
!!   Calculates \int \exp^{-a*(x-c)^2} x^l dx
!!   this works ALSO when c/=0.d0
!! SOURCE
!!
function gauint(a,c,l)
  use module_base
  implicit none
  integer, intent(in) :: l
  real(gp), intent(in) :: a,c
  real(gp) :: gauint
  !local variables
  real(gp), parameter :: gammaonehalf=1.772453850905516027298d0
  integer :: p
  real(gp) :: rfac,prefac,stot,fsum,tt,firstprod
  !build the prefactor
  prefac=sqrt(a)
  prefac=1.d0/prefac
  prefac=gammaonehalf*prefac

  !the first term of the sum is one
  !but we have to multiply for the prefactor
  stot=c**l

  !calculate the sum
  do p=1,l/4
     tt=rfac(p+1,2*p)
     fsum=rfac(l-2*p+1,l)
     fsum=fsum/tt
     tt=firstprod(p)
     fsum=fsum*tt
     tt=c**(l-2*p)
     tt=tt/a**p
     fsum=fsum*tt
     stot=stot+fsum
  end do
  do p=l/4+1,l/3
     tt=rfac(p+1,l-2*p)
     fsum=rfac(2*p+1,l)
     fsum=fsum/tt
     tt=firstprod(p)
     fsum=fsum*tt
     tt=c**(l-2*p)
     tt=tt/a**p
     fsum=fsum*tt
     stot=stot+fsum
  end do
  do p=l/3+1,l/2
     tt=rfac(l-2*p+1,p)
     fsum=rfac(2*p+1,l)
     fsum=fsum*tt
     tt=firstprod(p)
     fsum=fsum*tt
     tt=c**(l-2*p)
     tt=tt/a**p
     fsum=fsum*tt
     stot=stot+fsum
  end do

  !final result
  gauint=stot*prefac
  
END FUNCTION gauint
!!***


!!****f* BigDFT/gauint0
!! FUNCTION
!!   Calculates \int \exp^{-a*x^2} x^l dx
!!   this works only when l is even (if not equal to zero)
!! SOURCE
!!
function gauint0(a,l)
  use module_base
  implicit none
  integer, intent(in) :: l
  real(gp), intent(in) :: a
  real(gp) :: gauint0
  !local variables
  real(gp), parameter :: gammaonehalf=1.772453850905516027298d0
  integer :: p
  real(gp) :: xfac,prefac,tt
  !build the prefactor
  prefac=sqrt(a)
  prefac=1.d0/prefac
  prefac=gammaonehalf*prefac**(l+1)

  p=l/2
  tt=xfac(1,p,-0.5d0)
  !final result
  gauint0=prefac*tt
  
END FUNCTION gauint0
!!***


!!****f* BigDFT/firstprod
!! FUNCTION
!!
!! SOURCE
!!
function firstprod(p)
  use module_base
  implicit none
  integer, intent(in) :: p
  real(gp) :: firstprod
  !local variables
  integer :: i
  real(gp) :: tt
  firstprod=1.0_gp
  do i=1,p
     tt=real(2*i,gp)
     tt=1.0_gp/tt
     tt=1.0_gp-tt
     firstprod=firstprod*tt
  end do
END FUNCTION firstprod
!!***


!!****f* BigDFT/rfac
!! FUNCTION
!!
!! SOURCE
!!
function rfac(is,ie)
  use module_base
  implicit none
  integer, intent(in) :: is,ie
  real(gp) :: rfac
  !local variables
  integer :: i
  real(gp) :: tt
  rfac=1.d0
  do i=is,ie
     tt=real(i,gp)
     rfac=rfac*tt
  end do
END FUNCTION rfac
!!***


!!****f* BigDFT/xfac
!! FUNCTION
!!   With this function n!=xfac(1,n,0.d0)
!! SOURCE
!!
function xfac(is,ie,sh)
  use module_base
  implicit none
  integer, intent(in) :: is,ie
  real(gp), intent(in) :: sh
  real(gp) :: xfac
  !local variables
  integer :: i
  real(gp) :: tt
  xfac=1.d0
  do i=is,ie
     tt=real(i,gp)+sh
     xfac=xfac*tt
  end do
END FUNCTION xfac
!!***


!End of the interesting part


!!****f* BigDFT/gauinti
!! FUNCTION
!!   The same function but with integer factorials (valid ONLY if l<=18)
!!   not a visible improvement in speed with respect to the analogous real
!! SOURCE
!!
function gauinti(a,c,l)
  use module_base
  implicit none
  integer, intent(in) :: l
  real(gp), intent(in) :: a,c
  real(gp) :: gauinti
  !local variables
  real(gp), parameter :: gammaonehalf=1.772453850905516027298d0
  integer :: p,ifac
  real(gp) :: prefac,xsum,stot,fsum,tt,firstprod
  !build the prefactor
  prefac=sqrt(a)
  prefac=c**l/prefac
  prefac=gammaonehalf*prefac

  !object in the sum
  xsum=a*c**2
  xsum=1.d0/xsum

  !the first term of the sum is one
  stot=1.d0

  !calculate the sum
  do p=1,l/4
     tt=real(ifac(p+1,2*p),gp)
     fsum=real(ifac(l-2*p+1,l),gp)
     fsum=fsum/tt
     tt=firstprod(p)
     fsum=fsum*tt
     fsum=fsum*xsum**p
     stot=stot+fsum
  end do
  do p=l/4+1,l/3
     tt=real(ifac(p+1,l-2*p),gp)
     fsum=real(ifac(2*p+1,l),gp)
     fsum=fsum/tt
     tt=firstprod(p)
     fsum=fsum*tt
     fsum=fsum*xsum**p
     stot=stot+fsum
  end do
  do p=l/3+1,l/2
     tt=real(ifac(l-2*p+1,p),gp)
     fsum=real(ifac(2*p+1,l),gp)
     fsum=fsum*tt
     tt=firstprod(p)
     fsum=fsum*tt
     fsum=fsum*xsum**p
     stot=stot+fsum
  end do

  !final result
  gauinti=stot*prefac
  
END FUNCTION gauinti
!!***


!!****f* BigDFT/secondprod1
!! FUNCTION
!!   Valid if p<l/4 AND p/=0
!! SOURCE
!!
function secondprod1(p,l)
  use module_base
  implicit none
  integer, intent(in) :: p,l
  real(gp) :: secondprod1
  !local variables
  !integer :: i
  real(gp) :: tt,part1,rfac
  part1=rfac(p+1,2*p)
  !divide by the last value
  part1=real(l,gp)/part1
  tt=rfac(l-2*p+1,l-1)
!!!  part1=1.d0
!!!  do i=p+1,2*p !in the second case the bound must be changed here
!!!     tt=real(i,gp)
!!!     part1=part1*tt
!!!  end do
  secondprod1=tt*part1
END FUNCTION secondprod1
!!***


!!****f* BigDFT/secondprod2
!! FUNCTION
!!   Valid if p>=l/4 AND p<l/3
!! SOURCE
!!
function secondprod2(p,l)
  use module_base
  implicit none
  integer, intent(in) :: p,l
  real(gp) :: secondprod2
  !local variables
  !integer :: i
  real(gp) :: tt,part1,rfac
  part1=rfac(p+1,l-2*p)
  !divide by the last value
  part1=real(l,gp)/part1
  tt=rfac(2*p+1,l-1)
!!!  part1=1.d0
!!!  do i=p+1,2*p !in the second case the bound must be changed here
!!!     tt=real(i,gp)
!!!     part1=part1*tt
!!!  end do
  secondprod2=tt*part1
END FUNCTION secondprod2
!!***


!!****f* BigDFT/ifac
!! FUNCTION
!!   Integer version of factorial
!! SOURCE
!!
function ifac(is,ie)
  implicit none
  integer, intent(in) :: is,ie
  integer :: ifac
  !local variables
  integer :: i
  ifac=1
  do i=is,ie
     ifac=ifac*i
  end do
END FUNCTION ifac
!!***


!!****f* BigDFT/wavelets_to_gaussians
!! FUNCTION
!!   Calculate the projection of norb wavefunctions on a gaussian basis set
!! SOURCE
!!
subroutine wavelets_to_gaussians(geocode,norbp,nspinor,n1,n2,n3,G,thetaphi,hx,hy,hz,wfd,psi,coeffs)
  use module_base
  use module_types
  implicit none
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: norbp,n1,n2,n3,nspinor
  real(gp), intent(in) :: hx,hy,hz
  type(gaussian_basis), intent(in) :: G
  type(wavefunctions_descriptors), intent(in) :: wfd
  real(gp), dimension(2,G%nat), intent(in) :: thetaphi
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,nspinor,norbp), intent(in) :: psi
  real(wp), dimension(G%ncoeff,nspinor,norbp), intent(out) :: coeffs
  !local variables
  integer :: iorb,ispinor
  
  do iorb=1,norbp
     do ispinor=1,nspinor
        call orbital_projection(geocode,n1,n2,n3,G%nat,G%rxyz,thetaphi,&
             G%nshell,G%ndoc,G%nam,G%xp,G%psiat,G%nshltot,G%nexpo,G%ncoeff,&
             hx,hy,hz,wfd,psi(1,ispinor,iorb),coeffs(1,ispinor,iorb))
        !print *,'iorb, coeffs',iorb,coeffs(:,1,iorb)
     end do
  end do
  
END SUBROUTINE wavelets_to_gaussians
!!***


!!****f* BigDFT/orbital_projection
!! FUNCTION
!!   Calculate the projection of a given orbital on a gaussian basis centered on 
!!   a set of points
!! SOURCE
!!
subroutine orbital_projection(geocode,n1,n2,n3,nat,rxyz,thetaphi,nshell,ndoc,nam,xp,psiat,&
     nshltot,nexpo,ncoeff,hx,hy,hz,wfd,psi,coeffs)
  use module_base
  use module_types
  implicit none
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: n1,n2,n3,nat,nshltot,nexpo,ncoeff
  real(gp), intent(in) :: hx,hy,hz
  type(wavefunctions_descriptors), intent(in) :: wfd
  integer, dimension(nat), intent(in) :: nshell
  integer, dimension(nshltot), intent(in) :: ndoc,nam
  real(gp), dimension(nexpo), intent(in) :: xp,psiat
  real(gp), dimension(2,nat), intent(in) :: thetaphi
  real(gp), dimension(3,nat), intent(in) :: rxyz
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f), intent(in) :: psi
  real(wp), dimension(ncoeff), intent(out) :: coeffs
  !local variables
  integer :: ishell,iexpo,icoeff,iat,isat,ng,l
  
  ishell=0
  iexpo=1
  icoeff=1
  do iat=1,nat
     do isat=1,nshell(iat)
        ishell=ishell+1
        ng=ndoc(ishell)
        l=nam(ishell)
        call lsh_projection(geocode,l,ng,xp(iexpo),psiat(iexpo),n1,n2,n3,rxyz(1,iat),&
             thetaphi(1,iat),hx,hy,hz,wfd,psi,coeffs(icoeff))
        iexpo=iexpo+ng
        icoeff=icoeff+2*l-1
     end do
  end do

  call gaudim_check(iexpo,icoeff,ishell,nexpo,ncoeff,nshltot)
  
END SUBROUTINE orbital_projection
!!***


!!****f* BigDFT/gaudim_check
!! FUNCTION
!!
!! SOURCE
!!
subroutine gaudim_check(iexpo,icoeff,ishell,nexpo,ncoeff,nshltot)
  implicit none
  integer, intent(in) :: iexpo,icoeff,ishell,nexpo,ncoeff,nshltot
  !check of the dimensions
  if (iexpo /= nexpo+1) then
     write(*,*)' ERROR: nexpo+1 <> iexpo',nexpo,iexpo
     stop
  else if (icoeff /= ncoeff+1) then
     write(*,*)' ERROR: ncoeff+1 <> icoeff',ncoeff,icoeff
     stop
  else if (ishell /= nshltot) then
     write(*,*)' ERROR: nshltot <> ishell',nshltot,ishell
     stop
  end if
END SUBROUTINE gaudim_check
!!***


!!****f* BigDFT/lsh_projection
!! FUNCTION
!!   Calculate the projection of a gaussian for a given eigenspace of spherical harmonics
!!   centered on a given point and rotated by theta(along z) and phi(along x)
!! SOURCE
!!
subroutine lsh_projection(geocode,l,ng,xp,psiat,n1,n2,n3,rxyz,thetaphi,hx,hy,hz,&
     wfd,psi,coeffs)
  use module_base
  use module_types
  implicit none
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: l,ng,n1,n2,n3
  real(gp), intent(in) :: hx,hy,hz
  type(wavefunctions_descriptors), intent(in) :: wfd
  real(gp), dimension(2), intent(in) :: thetaphi
  real(gp), dimension(3), intent(in) :: rxyz
  real(gp), dimension(ng), intent(in) :: xp,psiat
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f), intent(in) :: psi
  real(wp), dimension(2*l-1), intent(out) :: coeffs
  !local variables
  integer, parameter :: nterm_max=3 !to be enlarged for generic theta and phi
  integer :: nterm,m
  integer, dimension(nterm_max) :: lx,ly,lz
  real(gp), dimension(nterm_max) :: fac_arr
  
  do m=1,2*l-1
     !for the moment theta and phi are ignored but a new routine should be made
     call calc_coeff_inguess(l,m,nterm_max,nterm,lx,ly,lz,fac_arr)

     call wavetogau(geocode,n1,n2,n3,ng,nterm,lx,ly,lz,fac_arr,xp,psiat,&
          rxyz(1),rxyz(2),rxyz(3),hx,hy,hz,wfd%nseg_c,wfd%nvctr_c,wfd%keyg(1,1),wfd%keyv(1),&
          wfd%nseg_f,wfd%nvctr_f,wfd%keyg(1,wfd%nseg_c+1),wfd%keyv(wfd%nseg_c+1),&
          psi(1),psi(wfd%nvctr_c+1),coeffs(m))

     !print '(a,2(i4),5(1pe12.5))','l,m,rxyz,coeffs(m)',l,m,rxyz(:),coeffs(m)
  end do

  !here we should modify the coefficients following the direction of the axis
  call lsh_rotation(l-1,thetaphi(1),thetaphi(2),coeffs)

END SUBROUTINE lsh_projection
!!***


!!****f* BigDFT/lsh_rotation
!! FUNCTION
!!
!! SOURCE
!!
subroutine lsh_rotation(l,theta,phi,coeffs)
  use module_base
  implicit none
  integer, intent(in) :: l !beware the change in notation
  real(gp), intent(in) :: theta,phi
  real(wp), dimension(2*l+1), intent(inout) :: coeffs
  !local variables
  real(gp), parameter :: degrad=0.0174532925199432957692369076849_gp
  integer :: m,m1
  real(gp) :: t,p,res
  real(gp), dimension(7) :: incoef ! calculated projection
  real(gp), dimension(7,7) :: hrot ! rotation coefficients

  !quick return if possible
  if (theta == 0._gp .and. phi == 0._gp) then
     return
  end if

  !angles in radiants
  t=theta*degrad
  p=phi*degrad

  !extract coefficients for the rotation
  call rotation_matrix(l+1,t,p,hrot)
  
  !copy input variables
  do m=1,2*l+1
     incoef(m)=real(coeffs(m),gp)
  end do
  
  !apply rotation matrix
  do m=1,2*l+1
     res=0._gp
     do m1=1,2*l+1
        res=res+incoef(m1)*hrot(m,m1) !sense to be verified
     end do
     coeffs(m)=real(res,wp)
  end do
  
END SUBROUTINE lsh_rotation
!!***


!!****f* BigDFT/wavetogau
!! FUNCTION
!!   Calculate the scalar product between a sum of gaussians times polynomials and a wavefunction
!!   \int dx dy dz 
!!             \sum_i=1..ntp fac_arr(i) {
!!                  \sum_j=1..nterm psiat(j) [exp(-r^2/(2*(xp(j)^2)))] 
!!                      *((x-rx)^lx(i) *(y-ry)^ly(i) * (z-rz)^lz(i) ))} psi(x,y,z)
!!   Expressed in Daubechies Basis in the arrays psi_c, psi_f
!! SOURCE
!!
subroutine wavetogau(geocode,n1,n2,n3,nterm,ntp,lx,ly,lz,fac_arr,xp,psiat,rx,ry,rz,hx,hy,hz, & 
     nseg_c,mvctr_c,keyg_c,keyv_c,nseg_f,mvctr_f,keyg_f,keyv_f,psi_c,psi_f,overlap)
  use module_base
  implicit none
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: n1,n2,n3,nterm,ntp,nseg_c,nseg_f,mvctr_c,mvctr_f
  real(gp), intent(in) :: rx,ry,rz,hx,hy,hz
  integer, dimension(ntp), intent(in) :: lx,ly,lz
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(nseg_f), intent(in) :: keyv_f
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  integer, dimension(2,nseg_f), intent(in) :: keyg_f
  real(gp), dimension(ntp), intent(in) :: fac_arr
  real(gp), dimension(nterm), intent(in) :: xp,psiat
  real(wp), dimension(mvctr_c), intent(in) :: psi_c
  real(wp), dimension(7,mvctr_f), intent(in) :: psi_f
  real(wp), intent(out) :: overlap
  !local variables
  character(len=*), parameter :: subname='wavetogau'
  integer, parameter ::nw=32000
  logical :: perx,pery,perz
  integer:: iterm,itp,n_gau,ml1,mu1,ml2,mu2,ml3,mu3,i1,i2,i3,i_all,i_stat,iseg,ii,jj,j0,j1,i0,i
  real(wp) :: ovlp_c,ovlp_f1,ovlp_f2,ovlp_f3,ovlp_f4,ovlp_f5,ovlp_f6,ovlp_f7,ovlp
  real(gp) :: gau_a,te
  real(wp), dimension(0:nw,2) :: work
  real(wp), dimension(:,:), allocatable :: wprojx,wprojy,wprojz

  !conditions for periodicity in the three directions
  perx=(geocode /= 'F')
  pery=(geocode == 'P')
  perz=(geocode /= 'F')

  allocate(wprojx(0:n1,2+ndebug),stat=i_stat)
  call memocc(i_stat,wprojx,'wprojx',subname)
  allocate(wprojy(0:n2,2+ndebug),stat=i_stat)
  call memocc(i_stat,wprojy,'wprojy',subname)
  allocate(wprojz(0:n3,2+ndebug),stat=i_stat)
  call memocc(i_stat,wprojz,'wprojz',subname)

  overlap=0.0_wp

  do itp=1,ntp

     do iterm=1,nterm
        gau_a=xp(iterm)
        n_gau=lx(itp)
        call gauss_to_daub(hx,fac_arr(itp),rx,gau_a,n_gau,n1,ml1,mu1,wprojx(0,1),te,work,nw,&
             perx)
        n_gau=ly(itp)
        call gauss_to_daub(hy,1.0_gp,ry,gau_a,n_gau,n2,ml2,mu2,wprojy(0,1),te,work,nw,pery)
        n_gau=lz(itp)
        call gauss_to_daub(hz,psiat(iterm),rz,gau_a,n_gau,n3,ml3,mu3,wprojz(0,1),te,work,nw,&
             perz)

        !scalar product (in wavefunction precision)

        ! coarse part
        ovlp_c=0.0_wp
        do iseg=1,nseg_c
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
              ovlp_c=ovlp_c+psi_c(i-i0+jj)*wprojx(i,1)*wprojy(i2,1)*wprojz(i3,1)
           enddo
        enddo

        ! fine part
        ovlp_f1=0.0_wp
        ovlp_f2=0.0_wp
        ovlp_f3=0.0_wp
        ovlp_f4=0.0_wp
        ovlp_f5=0.0_wp
        ovlp_f6=0.0_wp
        ovlp_f7=0.0_wp
        do iseg=1,nseg_f
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
              !itp=itp+1
              ovlp_f1=ovlp_f1+psi_f(1,i-i0+jj)*wprojx(i,2)*wprojy(i2,1)*wprojz(i3,1)
              ovlp_f2=ovlp_f2+psi_f(2,i-i0+jj)*wprojx(i,1)*wprojy(i2,2)*wprojz(i3,1)
              ovlp_f3=ovlp_f3+psi_f(3,i-i0+jj)*wprojx(i,2)*wprojy(i2,2)*wprojz(i3,1)
              ovlp_f4=ovlp_f4+psi_f(4,i-i0+jj)*wprojx(i,1)*wprojy(i2,1)*wprojz(i3,2)
              ovlp_f5=ovlp_f5+psi_f(5,i-i0+jj)*wprojx(i,2)*wprojy(i2,1)*wprojz(i3,2)
              ovlp_f6=ovlp_f6+psi_f(6,i-i0+jj)*wprojx(i,1)*wprojy(i2,2)*wprojz(i3,2)
              ovlp_f7=ovlp_f7+psi_f(7,i-i0+jj)*wprojx(i,2)*wprojy(i2,2)*wprojz(i3,2)
           enddo
        enddo

        ovlp=ovlp_c+ovlp_f1+ovlp_f2+ovlp_f3+ovlp_f4+ovlp_f5+ovlp_f6+ovlp_f7

        overlap=overlap+ovlp

     end do

  end do

  i_all=-product(shape(wprojx))*kind(wprojx)
  deallocate(wprojx,stat=i_stat)
  call memocc(i_stat,i_all,'wprojx',subname)
  i_all=-product(shape(wprojy))*kind(wprojy)
  deallocate(wprojy,stat=i_stat)
  call memocc(i_stat,i_all,'wprojy',subname)
  i_all=-product(shape(wprojz))*kind(wprojz)
  deallocate(wprojz,stat=i_stat)
  call memocc(i_stat,i_all,'wprojz',subname)

END SUBROUTINE wavetogau
!!***


!!****f* BigDFT/rotation_matrix
!! FUNCTION
!!   Coefficients of the rotation matrix
!! SOURCE
!!
subroutine rotation_matrix(l,t,p,hrot)
  use module_base
  implicit none
  integer, intent(in) :: l
  real(gp), intent(in) :: t,p
  real(gp), dimension(2*l+1,2*l+1), intent(out) :: hrot
  !local variables
  real(gp) :: cp,c2p,c3p,sp,s2p,s3p,ct,c2t,c3t,st,s2t,s3t

  !trigonometrical functions
  cp=cos(p)
  sp=sin(p)
  c2p=cos(2._gp*p)
  s2p=sin(2._gp*p)
  c3p=cos(3._gp*p)
  s3p=sin(3._gp*p)
  ct=cos(t)
  st=sin(t)
  c2t=cos(2._gp*t)
  s2t=sin(2._gp*t)
  c3t=cos(3._gp*t)
  s3t=sin(3._gp*t)
  
  if (l == 1) then
     hrot(1,1)=cp*ct
     hrot(2,1)=ct*sp
     hrot(3,1)=-1._gp*st
     hrot(1,2)=-1._gp*sp
     hrot(2,2)=cp
     hrot(3,2)=0_gp
     hrot(1,3)=cp*st
     hrot(2,3)=sp*st
     hrot(3,3)=ct
  else if (l == 2) then
     hrot(1,1)=cp*ct
     hrot(2,1)=-1._gp*ct*sp
     hrot(3,1)=c2p*st
     hrot(4,1)=-2._gp*cp*sp*st
     hrot(5,1)=0_gp
     hrot(1,2)=c2t*sp
     hrot(2,2)=c2t*cp
     hrot(3,2)=2._gp*cp*ct*sp*st
     hrot(4,2)=0.5_gp*c2p*s2t
     hrot(5,2)=-1.7320508075688772935_gp*ct*st
     hrot(1,3)=-1._gp*cp*st
     hrot(2,3)=sp*st
     hrot(3,3)=c2p*ct
     hrot(4,3)=-2._gp*cp*ct*sp
     hrot(5,3)=0_gp
     hrot(1,4)=-1._gp*ct*sp*st
     hrot(2,4)=-1._gp*cp*ct*st
     hrot(3,4)=cp*(1._gp + ct**2)*sp
     hrot(4,4)=0.25_gp*c2p*(3._gp + c2t)
     hrot(5,4)=0.86602540378443864676_gp*st**2
     hrot(1,5)=1.7320508075688772935_gp*ct*sp*st
     hrot(2,5)=1.7320508075688772935_gp*cp*ct*st
     hrot(3,5)=1.7320508075688772935_gp*cp*sp*st**2
     hrot(4,5)=0.86602540378443864676_gp*c2p*st**2
     hrot(5,5)=0.25_gp*(1._gp + 3._gp*c2t)
  else if (l == 3) then
     hrot(1,1)=0.0625_gp*cp*(15._gp*c3t + ct)
     hrot(2,1)=0.0625_gp*(15._gp*c3t + ct)*sp
     hrot(3,1)=-0.15309310892394863114_gp*(5._gp*s3t + st)
     hrot(4,1)=-0.96824583655185422129_gp*c3p*ct*st**2
     hrot(5,1)=0.96824583655185422129_gp*ct*s3p*st**2
     hrot(6,1)=0.19764235376052370825_gp*c2p*(-3._gp*s3t + st)
     hrot(7,1)=-0.3952847075210474165_gp*(1._gp + 3._gp*c2t)*s2p*st
     hrot(1,2)=-0.125_gp*(3._gp + 5._gp*c2t)*sp
     hrot(2,2)=0.125_gp*(3._gp + 5._gp*c2t)*cp
     hrot(3,2)=0_gp
     hrot(4,2)=0.96824583655185422129_gp*s3p*st**2
     hrot(5,2)=0.96824583655185422129_gp*c3p*st**2
     hrot(6,2)=3.162277660168379332_gp*cp*ct*sp*st
     hrot(7,2)=-0.790569415042094833_gp*c2p*s2t
     hrot(1,3)=0.15309310892394863114_gp*cp*(5._gp*s3t + st)
     hrot(2,3)=0.15309310892394863114_gp*sp*(5._gp*s3t + st)
     hrot(3,3)=0.125_gp*(5._gp*c3t + 3._gp*ct)
     hrot(4,3)=0.790569415042094833_gp*(1._gp - 2._gp*c2p)*cp*st**3
     hrot(5,3)=0.790569415042094833_gp*(1._gp + 2._gp*c2p)*sp*st**3
     hrot(6,3)=-1.9364916731037084426_gp*c2p*ct*st**2
     hrot(7,3)=-3.8729833462074168852_gp*cp*ct*sp*st**2
     hrot(1,4)=-0.96824583655185422129_gp*cp*ct*st**2
     hrot(2,4)=-0.96824583655185422129_gp*ct*sp*st**2
     hrot(3,4)=0.790569415042094833_gp*st**3
     hrot(4,4)=0.0625_gp*c3p*(c3t + 15._gp*ct)
     hrot(5,4)=-0.0625_gp*(c3t + 15._gp*ct)*s3p
     hrot(6,4)=-0.15309310892394863114_gp*c2p*(s3t + 5._gp*st)
     hrot(7,4)=-0.15309310892394863114_gp*s2p*(s3t + 5._gp*st)
     hrot(1,5)=-0.96824583655185422129_gp*sp*st**2
     hrot(2,5)=0.96824583655185422129_gp*cp*st**2
     hrot(3,5)=0_gp
     hrot(4,5)=0.125_gp*(5._gp + 3._gp*c2t)*s3p
     hrot(5,5)=0.125_gp*(5._gp + 3._gp*c2t)*c3p
     hrot(6,5)=-2.4494897427831780982_gp*cp*ct*sp*st
     hrot(7,5)=0.61237243569579452455_gp*c2p*s2t
     hrot(1,6)=-0.19764235376052370825_gp*cp*(-3._gp*s3t + st)
     hrot(2,6)=-0.19764235376052370825_gp*sp*(-3._gp*s3t + st)
     hrot(3,6)=-1.9364916731037084426_gp*ct*st**2
     hrot(4,6)=0.15309310892394863114_gp*c3p*(s3t + 5._gp*st)
     hrot(5,6)=-0.15309310892394863114_gp*s3p*(s3t + 5._gp*st)
     hrot(6,6)=0.25_gp*c2p*(1._gp + 3._gp*c2t)*ct
     hrot(7,6)=0.25_gp*(1._gp + 3._gp*c2t)*ct*s2p
     hrot(1,7)=-1.581138830084189666_gp*ct*sp*st
     hrot(2,7)=1.581138830084189666_gp*cp*ct*st
     hrot(3,7)=0_gp
     hrot(4,7)=-0.61237243569579452455_gp*s2t*s3p
     hrot(5,7)=-0.61237243569579452455_gp*c3p*s2t
     hrot(6,7)=-1._gp*c2t*s2p
     hrot(7,7)=c2p*c2t
  else
     write(*,*)'ERROR: rotation matrix for angular momentum l=',l,'not implemented'
     stop
  end if
END SUBROUTINE rotation_matrix
!!***
