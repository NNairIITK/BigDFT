!> @file
!!   Routines to handle Gaussian basis set
!! @author
!!    Copyright (C) 2007-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Restart from gaussian functions
subroutine restart_from_gaussians(iproc,nproc,orbs,Lzd,hx,hy,hz,psi,G,coeffs)
  use module_base
  use module_types
  use gaussians, only: gaussian_basis, deallocate_gwf
  implicit none
  integer, intent(in) :: iproc,nproc
  real(gp), intent(in) :: hx,hy,hz
  type(orbitals_data), intent(in) :: orbs
  type(local_zone_descriptors), intent(in) :: Lzd
  type(gaussian_basis), intent(inout) :: G
  real(wp), dimension(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f,orbs%norbp*orbs%nspinor), intent(out) :: psi
  real(wp), dimension(:,:), pointer :: coeffs
  !local variables
  character(len=*), parameter :: subname='restart_from_gaussians'

  !the atomic positions corresponds to the new values
  !calculate the dual coefficients with the new positions
  
  !call gaussian_orthogonality(iproc,nproc,norb,norbp,G,coeffs)

  call dual_gaussian_coefficients(orbs%norbp,G,coeffs)

  !call gaussians_to_wavelets(iproc,nproc,lr%geocode,orbs,lr%d,hx,hy,hz,lr%wfd,G,coeffs,psi)
  call gaussians_to_wavelets_new(iproc,nproc,Lzd,orbs,G,coeffs,psi)

  !deallocate gaussian structure and coefficients
  call deallocate_gwf(G)
  call f_free_ptr(coeffs)

  nullify(G%rxyz)

END SUBROUTINE restart_from_gaussians


!> Read information for gaussian basis set (from CP2K) or for restarting
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
  integer :: jexpo,iexpo,iat,iorb,jat,icoeff,jcoeff,jorb,j
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
  G%ncplx=1 !2 only for PAW or XANES

  G%nshell = f_malloc_ptr(G%nat,id='G%nshell')
  G%nam = f_malloc_ptr(G%nshltot,id='G%nam')
  G%ndoc = f_malloc_ptr(G%nshltot,id='G%ndoc')
  G%xp = f_malloc_ptr((/ G%ncplx, G%nexpo /),id='G%xp')
  G%psiat = f_malloc_ptr((/ G%ncplx, G%nexpo /),id='G%psiat')

  coeffs = f_malloc_ptr((/ G%ncoeff, orbs%norbp*orbs%nspinor /),id='coeffs')
  
  if(fillrxyz) then

     G%rxyz  = f_malloc_ptr((/ 3, G%nat /),id='G%rxyz ')
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
     read(99,*)jexpo,G%xp(1,jexpo),G%psiat(1,jexpo)
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



!>   Write gaussian informatio for another program or for restarting
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
  integer :: jproc,ierr,iexpo,iat,iorb,icoeff,jorb,j,norb_tot
  integer, dimension(:,:), allocatable :: gatherarr
  real(gp), dimension(:,:), allocatable :: gaupsi

  gaupsi = f_malloc((/ G%ncoeff, orbs%norb*orbs%nspinor /),id='gaupsi')


  if (nproc > 1) then
     gatherarr = f_malloc((/ 0.to.nproc-1, 1.to.2 /),id='gatherarr')
     
     norb_tot=0
     gatherarr(0,1)=G%ncoeff*orbs%norb_par(0,0)*orbs%nspinor
     gatherarr(0,2)=G%ncoeff*norb_tot*orbs%nspinor
     !gather the coefficients in a unique array
     do jproc=1,nproc-1
        norb_tot=norb_tot+orbs%norb_par(jproc-1,0)
        gatherarr(jproc,1)=G%ncoeff*orbs%norb_par(jproc,0)
        gatherarr(jproc,2)=G%ncoeff*norb_tot*orbs%nspinor
     end do

     call MPI_GATHERV(coeffs,gatherarr(iproc,1),mpidtypw,gaupsi,gatherarr(0,1),gatherarr(0,2),&
          mpidtypw,0,bigdft_mpi%mpi_comm,ierr)

     call f_free(gatherarr)
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
        write(99,'(i6,2(1x,1pe21.14))')iexpo,G%xp(1,iexpo),G%psiat(1,iexpo)
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

  call f_free(gaupsi)
  
END SUBROUTINE write_gaussian_information


!> Create gaussian structure from input guess pseudo wavefunctions
subroutine gaussian_pswf_basis(ng,enlargerprb,iproc,nspin,at,rxyz,G,Gocc, gaenes, &
     iorbtolr,iorbto_l, iorbto_m,  iorbto_ishell,iorbto_iexpobeg)
  use module_base
  use ao_inguess, only: iguess_generator,print_eleconf,ao_nspin_ig,nmax_occ_ao
  use module_types
  use yaml_output
  implicit none
  logical, intent(in) :: enlargerprb
  integer, intent(in) :: iproc,nspin,ng
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%astruct%nat), target, intent(in) :: rxyz
  type(gaussian_basis), intent(inout) :: G
  real(wp), dimension(:), pointer :: Gocc

  !! the following arguments are used when building PPD : the preconditioner for CG spectra
  real(gp), dimension(:), pointer, optional :: gaenes
  integer, dimension(:), pointer, optional :: iorbtolr
  integer, dimension(:), pointer, optional :: iorbto_l
  integer, dimension(:), pointer, optional :: iorbto_m
  integer, dimension(:), pointer, optional :: iorbto_ishell
  integer, dimension(:), pointer, optional :: iorbto_iexpobeg

  !local variables
  character(len=*), parameter :: subname='gaussian_pswf_basis'
  !integer, parameter :: noccmax=2,lmax=4,nelecmax=32 !n(c) nmax=6
  logical :: occeq
  integer :: iat,ityp,ishell,iexpo,l,i,ig,ictotpsi,norbe,norbsc,ishltmp
  integer :: ityx,ntypesx,nspinor,jat,noncoll,icoeff,iocc,nlo,ispin,m,icoll,ngv,ngc,islcc
  real(gp) :: ek
  !integer, dimension(lmax) :: nl
  !real(gp), dimension(noccmax,lmax) :: occup
  integer, dimension(:), allocatable :: iatypex
  integer, dimension(:,:), allocatable :: norbsc_arr
  real(gp), dimension(:), allocatable :: psiatn,locrad
  real(gp), dimension(:,:), allocatable :: xpt
  real(gp), dimension(:,:,:), allocatable :: psiat  

  !! auxiliary variables used when creating optional arrays for PPD
  real(gp)  :: gaenes_aux(nmax_occ_ao*at%astruct%nat)
  integer :: last_aux, firstperityx(at%astruct%nat)
  integer :: nspin_print !< to be removed, shouldpass in input variables


  !quick return if possible
  !if the positions are already associated it means that the basis is generated
  if (associated(G%rxyz)) then
     nullify(Gocc) !to avoid problem with the initialization
     return
  end if

  norbsc_arr = f_malloc((/ at%natsc+1, 1 /),id='norbsc_arr')
  locrad = f_malloc(at%astruct%nat,id='locrad')

  !for the moment, only collinear
  nspinor=1
  !if non-collinear it is like nspin=1 but with the double of orbitals
  if (nspinor == 4) then
     noncoll=2
  else
     noncoll=1
  end if

  nspin_print=ao_nspin_ig(nspin,nspinor=nspinor)

  call readAtomicOrbitals(at,norbe,norbsc,nspin,nspinor,norbsc_arr,locrad)

  call f_free(locrad)

  !Generate the input guess via the inguess_generator
  !take also into account the IG polarisations

  !the number of gaussian centers are thus nat
  G%nat=at%astruct%nat
  !this pointing creates problems if at the next call the positions are given by a different array.
  !presumably the best if to copy the values and allocate the pointer
  G%rxyz => rxyz

  !copy the parsed values in the gaussian structure
  !count also the total number of shells
  G%nshell = f_malloc_ptr(at%astruct%nat,id='G%nshell')

  !calculate the number of atom types by taking into account the occupation
  iatypex = f_malloc(at%astruct%nat,id='iatypex')
  
  ntypesx=0
  G%nshltot=0
  count_shells: do iat=1,at%astruct%nat
     ityp=at%astruct%iatype(iat)
     !call count_atomic_shells(nspin_print,at%aoig(iat)%aocc,occup,nl)
     G%nshell(iat)=sum(at%aoig(iat)%nl)!(nl(1)+nl(2)+nl(3)+nl(4))
     G%nshltot=G%nshltot+G%nshell(iat)
     !check the occupation numbers and the atoms type
     !once you find something equal exit the procedure
     do jat=1,iat-1
        if (at%astruct%iatype(jat) == ityp) then
           occeq= all(at%aoig(jat)%aocc == at%aoig(iat)%aocc)!.true.
           !do i=1,nelecmax
           !   occeq = occeq .and. &
           !        (at%aoig(jat)%aocc(i) == at%aoig(iat)%aocc(i))
           !end do
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

  G%ndoc = f_malloc_ptr(G%nshltot,id='G%ndoc')
  G%nam = f_malloc_ptr(G%nshltot,id='G%nam')

  !the default value for the gaussians is chosen to be 21
  xpt = f_malloc((/ ng, ntypesx /),id='xpt')
  psiat = f_malloc((/ ng , nmax_occ_ao , ntypesx /),id='psiat')
  psiatn = f_malloc(ng,id='psiatn')


  !assign shell IDs and count the number of exponents and coefficients
  G%ncplx=1
  G%nexpo=0
  G%ncoeff=0
  ishell=0
  ntypesx=0
  do iat=1,at%astruct%nat
     ityp=at%astruct%iatype(iat)
     ityx=iatypex(iat)
     ishltmp=0
     !call count_atomic_shells(nspin_print,at%aoig(iat)%aocc,occup,nl)
     if (ityx > ntypesx) then
        if (iproc == 0 .and. verbose > 1) then
           call yaml_map('Generation of input wavefunction data for atom ', trim(at%astruct%atomnames(ityp)))
           call print_eleconf(nspin_print,&
                at%aoig(iat)%aocc,at%aoig(iat)%nl_sc)
        end if

        firstperityx( ityx)=iat
        !positions for the nlcc arrays
        call nlcc_start_position(ityp,at,ngv,ngc,islcc)
         !eliminate the nlcc parameters from the IG, since XC is always LDA
         ngv=0
         ngc=0


        if( present(gaenes)) then
           call iguess_generator(at%nzatom(ityp),at%nelpsp(ityp),& !_modified
                real(at%nelpsp(ityp),gp),nspin_print,at%aoig(iat)%aocc,at%psppar(0:,0:,ityp),&
                at%npspcode(ityp),ngv,ngc,at%nlccpar(0:,max(islcc,1)),&
                ng-1,xpt(1,ityx),psiat(1,1,ityx),enlargerprb, &
                gaenes_aux=gaenes_aux(1+nmax_occ_ao*( firstperityx( ityx)-1))  )
        else
           call iguess_generator(at%nzatom(ityp),at%nelpsp(ityp),&
                real(at%nelpsp(ityp),gp),nspin_print,at%aoig(iat)%aocc,at%psppar(0:,0:,ityp),&
                at%npspcode(ityp),ngv,ngc,at%nlccpar(0:,max(islcc,1)),&
                ng-1,xpt(1,ityx),psiat(1,1,ityx),enlargerprb)
        endif

        ntypesx=ntypesx+1
        !if (iproc == 0 .and. verbose > 1) write(*,'(1x,a)')'done.'
     end if

     do l=1,size(at%aoig(iat)%nl)
        do i=1,at%aoig(iat)%nl(l-1)!nl(l)
           ishell=ishell+1
           ishltmp=ishltmp+1
           G%ndoc(ishell)=ng!(ity)
           G%nam(ishell)=l
           G%nexpo=G%nexpo+ng!(ity)
           G%ncoeff=G%ncoeff+2*l-1
           !print *,'iat,i,l',iat,i,l,norbe,G%ncoeff
           if( present(gaenes)) then
              gaenes_aux(ishltmp+nmax_occ_ao*(iat-1))=gaenes_aux(ishltmp+nmax_occ_ao*(firstperityx(ityx)-1))
           endif
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
  Gocc = f_malloc0_ptr(G%ncoeff,id='Gocc')
  !call to_zero(G%ncoeff,Gocc)

  if( present(gaenes)) then

     gaenes = f_malloc0_ptr(G%ncoeff,id='gaenes')
     !call to_zero(G%ncoeff,gaenes)

     iorbtolr = f_malloc_ptr(G%ncoeff,id='iorbtolr')
     iorbto_l = f_malloc_ptr(G%ncoeff,id='iorbto_l')
     iorbto_m = f_malloc_ptr(G%ncoeff,id='iorbto_m')
     iorbto_ishell = f_malloc_ptr(G%ncoeff,id='iorbto_ishell')
     iorbto_iexpobeg = f_malloc_ptr(G%ncoeff,id='iorbto_iexpobeg')
     
  endif

  !allocate and assign the exponents and the coefficients
  G%psiat = f_malloc_ptr((/ G%ncplx, G%nexpo /),id='G%psiat')
  G%xp = f_malloc_ptr((/ G%ncplx, G%nexpo /),id='G%xp')

  ishell=0
  iexpo=0
  icoeff=1
  do iat=1,at%astruct%nat
     if( present(gaenes))  last_aux=ishell
     !print *, 'debug',iat,present(gaenes),nspin,noncoll
     ityp=at%astruct%iatype(iat)
     ityx=iatypex(iat)
     !call count_atomic_shells(ao_nspin_ig(nspin,nspinor=nspinor),&
     !     at%aoig(iat)%aocc,occup,nl)
     ictotpsi=0
     iocc=0
     do l=1,size(at%aoig(iat)%nl)
        iocc=iocc+1
        nlo=nint(at%aoig(iat)%aocc(iocc)) !just to increase the counting 
        do i=1,at%aoig(iat)%nl(l-1)
           ishell=ishell+1
           ictotpsi=ictotpsi+1
           call atomkin(l-1,ng,xpt(1,ityx),psiat(1,ictotpsi,ityx),psiatn,ek)
           do ig=1,G%ndoc(ishell)
              iexpo=iexpo+1
              G%psiat(1,iexpo)=psiatn(ig)*sign(1.0_gp,psiatn(1))
              G%xp(1,iexpo)=xpt(ig,ityp)
           end do

           do ispin=1,nspin
              do m=1,2*l-1
                 !each orbital has two electrons in the case of the 
                 !non-collinear case
                 do icoll=1,noncoll !non-trivial only for nspinor=4
                    iocc=iocc+1
                    Gocc(icoeff)=Gocc(icoeff)+at%aoig(iat)%aocc(iocc)
                    !print *,'test',iocc,icoeff,shape(at%aocc),'test2',shape(Gocc)
                    if( present(gaenes)) then
                        gaenes(icoeff)=gaenes_aux( ishell-last_aux+  nmax_occ_ao*(iat-1) )
                        iorbtolr       (icoeff)=iat
                        iorbto_l       (icoeff)=l        
                        iorbto_m       (icoeff)=m
                        iorbto_ishell  (icoeff)=ishell
                        iorbto_iexpobeg(icoeff)=iexpo - G%ndoc(ishell)  +1
                        
                    endif
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

  call f_free(xpt)
  call f_free(psiat)
  call f_free(norbsc_arr)
  call f_free(psiatn)
  call f_free(iatypex)

END SUBROUTINE gaussian_pswf_basis


!> Create gaussian structure from ptildes   _for_paw
subroutine gaussian_pswf_basis_for_paw(at,rxyz,G,  &
     iorbtolr,iorbto_l, iorbto_m,  iorbto_ishell,iorbto_iexpobeg, iorbto_paw_nchannels, iorbto_imatrixbeg )
  use module_base
  use module_types
  implicit none
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%astruct%nat), target, intent(in) :: rxyz
  type(gaussian_basis_c), intent(inout) :: G

  integer, pointer :: iorbtolr(:)
  integer, pointer :: iorbto_l(:)
  integer, pointer :: iorbto_paw_nchannels(:)
  integer, pointer :: iorbto_m(:)
  integer, pointer :: iorbto_ishell(:)
  integer, pointer :: iorbto_iexpobeg(:)
  integer, pointer :: iorbto_imatrixbeg(:)

  !local variables
  character(len=*), parameter :: subname='gaussian_pswf_basis_for_paw'
  !integer, parameter :: noccmax=2,lmax=4,nmax=6,nelecmax=32
  integer :: il, j
  integer :: iat,ityp,jtyp,ishell,iexpo,l,i
  integer :: ig,iexpoat_qs , iexpoat_coeffs
  integer :: natpaw, imatrix
  integer :: nspinor,noncoll,icoeff,m

  !quick return if possible
  !if the positions are already associated it means that the basis is generated
  if (associated(G%rxyz)) then
     return
  end if


  !for the moment, only collinear
  nspinor=1
  !if non-collinear it is like nspin=1 but with the double of orbitals
  if (nspinor == 4) then
     noncoll=2
  else
     noncoll=1
  end if

  natpaw=0
  do iat=1, at%astruct%nat
     if(  at%paw_NofL(at%astruct%iatype(iat)).gt.0) then
        natpaw=natpaw+1
     end if
  end do

  !the number of gaussian centers are thus natpaw
  G%nat= natpaw   
  G%rxyz  = f_malloc_ptr((/ 3, natpaw /),id='G%rxyz ')

  natpaw=0
  do iat=1, at%astruct%nat
     if(  at%paw_NofL(at%astruct%iatype(iat)).gt.0) then
        natpaw=natpaw+1
        G%rxyz (:, natpaw) = rxyz(:,iat)
     end if
  end do

  G%nshell = f_malloc_ptr(natpaw  ,id='G%nshell')

  G%nshltot=0
  natpaw=0
  count_shells: do iat=1,at%astruct%nat
     ityp=at%astruct%iatype(iat)
     if(  at%paw_NofL(ityp).gt.0) then
        natpaw=natpaw+1
        G%nshell(natpaw)=0
        il=0
        do jtyp=1,ityp-1
           il=il+at%paw_NofL(jtyp)
        enddo
        do i=1, at%paw_NofL(ityp)
           il=il+1
           G%nshell(natpaw)=G%nshell(natpaw)+at%paw_nofchannels(il)
        end do
        G%nshltot=G%nshltot+G%nshell(natpaw)
     end if
  end do count_shells

  G%ndoc = f_malloc_ptr(G%nshltot,id='G%ndoc')
  G%nam = f_malloc_ptr(G%nshltot,id='G%nam')

  !assign shell IDs and count the number of exponents and coefficients
  G%nexpo=0
  G%ncoeff=0
  ishell=0
  il=0
  do iat=1,at%astruct%nat
     ityp=at%astruct%iatype(iat) 
     if(  at%paw_NofL(ityp).gt.0) then
        il=0
        do jtyp=1,ityp-1
           il=il+at%paw_NofL(jtyp)
        enddo
        do i=1, at%paw_NofL(ityp)
           il=il+1
           do j=1, at%paw_nofchannels(il)
              ishell=ishell+1
              G%ndoc(ishell)=at%paw_nofgaussians(il)

              if(at%paw_l(il).ge.0) then
                 G%nam(ishell)= at%paw_l(il) + 1 
              else
                 G%nam(ishell)=  at%paw_l(il)  
              endif
              G%nexpo=G%nexpo+ G%ndoc(ishell)
              G%ncoeff=G%ncoeff + abs(2*G%nam(ishell))   -1
           enddo
        end do
     end if
  end do


  iorbtolr = f_malloc_ptr(G%ncoeff,id='iorbtolr')
  iorbto_l = f_malloc_ptr(G%ncoeff,id='iorbto_l')
  iorbto_paw_nchannels = f_malloc_ptr(G%ncoeff,id='iorbto_paw_nchannels')
  iorbto_m = f_malloc_ptr(G%ncoeff,id='iorbto_m')
  iorbto_ishell = f_malloc_ptr(G%ncoeff,id='iorbto_ishell')
  iorbto_iexpobeg = f_malloc_ptr(G%ncoeff,id='iorbto_iexpobeg')
  iorbto_imatrixbeg = f_malloc_ptr(G%ncoeff,id='iorbto_imatrixbeg')


  !allocate and assign the exponents and the coefficients
  G%psiat = f_malloc_ptr(G%nexpo,id='G%psiat')

  G%expof = f_malloc_ptr(G%nexpo,id='G%expof')

  iexpo=0
  icoeff=0

  natpaw=0
  
  do iat=1,at%astruct%nat

     il=0
     ishell=0
     iexpoat_qs=0
     iexpoat_coeffs =0
     imatrix=1

     ityp=at%astruct%iatype(iat)

     if(  at%paw_NofL(ityp).gt.0) then
        natpaw=natpaw+1
        do jtyp=1,ityp-1
           do i=l, at%paw_NofL(jtyp)
              il=il+1
              ishell=ishell+at%paw_nofchannels(il)
              iexpoat_coeffs =iexpoat_coeffs+G%ndoc(il)*at%paw_nofchannels(il)
              iexpoat_qs=iexpoat_qs+G%ndoc(il)
              imatrix=imatrix+at%paw_nofchannels(il)*at%paw_nofchannels(il)
           end do
        enddo

        do i=1, at%paw_NofL(ityp)
           il=il+1
           do j=1, at%paw_nofchannels(il)
              ishell=ishell+1

              do ig=1,G%ndoc(ishell)
                 iexpo=iexpo+1
                 iexpoat_coeffs =iexpoat_coeffs +1

                 G%psiat(iexpo)=CMPLX(at%paw_Gcoeffs(2*iexpoat_coeffs -1) , at%paw_Gcoeffs(2*iexpoat_coeffs ) )   

                 G%expof (iexpo)   =CMPLX(at%paw_Greal(il), at%paw_Gimag( iexpoat_qs+ ig ) )
              enddo
              

              l=abs(G%nam(ishell))
              do m=1,2*l-1
                 icoeff=icoeff+1
                 iorbtolr       (icoeff)=natpaw
                 iorbto_l       (icoeff)=  G%nam(ishell)      
                 iorbto_paw_nchannels(icoeff)=  at%paw_nofchannels(il)      
                 iorbto_m       (icoeff)=m
                 iorbto_ishell  (icoeff)=ishell         
                 iorbto_iexpobeg(icoeff)=iexpo - G%ndoc(ishell)  +1
                 iorbto_imatrixbeg(icoeff) = imatrix
              end do
           enddo
           iexpoat_qs=iexpoat_qs+G%ndoc(il)
           imatrix=imatrix+at%paw_nofchannels(il)*at%paw_nofchannels(il)
        end do
     end if
  end do

!!  gaudim_check(iexpo,icoeff,ishell,nexpo,ncoeff,nshltot)

END SUBROUTINE gaussian_pswf_basis_for_paw


!>   Extract the pseudopotential basis
!! @warning
!!   This is not the complete PSP basis set. 
!!   The radial power term is lacking in the gaussian descriptors should be added if needed
subroutine gaussian_psp_basis(at,rxyz,G)
  use module_base
  use module_types
  implicit none
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%astruct%nat), target, intent(in) :: rxyz
  type(gaussian_basis), intent(out) :: G  
  !local variables
  character(len=*), parameter :: subname='gaussian_psp_basis'
  integer :: iat,nshell,ityp,iexpo,l,ishell

  G%nat=at%astruct%nat
  G%rxyz => rxyz
  G%nshell = f_malloc_ptr(at%astruct%nat,id='G%nshell')
 
  G%nshltot=0
  do iat=1,G%nat
     ityp=at%astruct%iatype(iat) 
     nshell=0
     do l=1,4 
        if (at%psppar(l,0,ityp) /= 0.0_gp) nshell=nshell+1
     enddo
     G%nshell(iat)=nshell
     G%nshltot=G%nshltot+nshell
  end do

  G%ndoc = f_malloc_ptr(G%nshltot,id='G%ndoc')
  G%nam = f_malloc_ptr(G%nshltot,id='G%nam')

  !assign shell IDs and count the number of exponents and coefficients
  G%ncplx=1
  G%nexpo=0
  G%ncoeff=0
  ishell=0
  do iat=1,G%nat
     ityp=at%astruct%iatype(iat)
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
  G%xp = f_malloc_ptr((/ G%ncplx, G%nexpo /),id='G%xp')
  G%psiat = f_malloc_ptr((/ G%ncplx, G%nexpo /),id='G%psiat')

  ishell=0
  iexpo=0
  do iat=1,G%nat
     ityp=at%astruct%iatype(iat)
     do l=1,4 
        if (at%psppar(l,0,ityp) /= 0.0_gp) then
           ishell=ishell+1
           iexpo=iexpo+1
           G%psiat(1,iexpo)=1.0_gp
           G%xp(1,iexpo)=at%psppar(l,0,ityp)
        end if
     end do
  end do

END SUBROUTINE gaussian_psp_basis



!>
!!
!!
subroutine gaussian_orthogonality(iproc,nproc,norb,norbp,G,coeffs)
  use module_base
  use module_types
  use gaussians
  implicit none
  integer, intent(in) :: iproc,nproc,norb,norbp
  type(gaussian_basis), intent(in) :: G
  real(wp), dimension(G%ncoeff,norbp), intent(inout) :: coeffs
  !local variables
  character(len=*), parameter :: subname='gaussian_orthogonality' 
  integer :: iorb,i,jproc,info,ierr
  integer, dimension(:,:), allocatable :: gatherarr
  real(gp), dimension(:,:), allocatable :: ovrlp,gaupsi,tmp,smat

  ovrlp = f_malloc((/ G%ncoeff, G%ncoeff /),id='ovrlp')
  gaupsi = f_malloc((/ G%ncoeff, norb /),id='gaupsi')
  tmp = f_malloc((/ G%ncoeff, norb /),id='tmp')
  smat = f_malloc((/ norb, norb /),id='smat')


  if (nproc > 1) then
     gatherarr = f_malloc((/ 0.to.nproc-1, 1.to.2 /),id='gatherarr')

     !gather the coefficients in a unique array
     do jproc=0,nproc-1
        gatherarr(jproc,1)=G%ncoeff*min(max(norb-jproc*norbp,0),norbp)
        gatherarr(jproc,2)=G%ncoeff*min(jproc*norbp,norb)
     end do

     call MPI_ALLGATHERV(coeffs,gatherarr(iproc,1),mpidtypw,gaupsi,gatherarr(0,1),gatherarr(0,2),&
          mpidtypw,bigdft_mpi%mpi_comm,ierr)

     call f_free(gatherarr)
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
  
  call f_free(ovrlp)
  call f_free(tmp)
  call f_free(smat)
  call f_free(gaupsi)

END SUBROUTINE gaussian_orthogonality


!>   Calculates @f$\int_0^\infty \exp^{-a*x^2} x^l dx@f$
!!
!!
function gauinth(a,l)
  use module_base
  implicit none
  integer, intent(in) :: l
  real(gp), intent(in) :: a
  real(gp) :: gauinth
  !local variables
  real(gp), parameter :: gammaonehalf=1.772453850905516027298_gp
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


!>
!!
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



!>   With this function n!=xfac(1,n,0.d0)
!!
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



!End of the interesting part
!!$
!!$
!!$!>   The same function but with integer factorials (valid ONLY if l<=18)
!!$!!   not a visible improvement in speed with respect to the analogous real
!!$!!
!!$!!
!!$function gauinti(a,c,l)
!!$  use module_base
!!$  implicit none
!!$  integer, intent(in) :: l
!!$  real(gp), intent(in) :: a,c
!!$  real(gp) :: gauinti
!!$  !local variables
!!$  real(gp), parameter :: gammaonehalf=1.772453850905516027298d0
!!$  integer :: p,ifac
!!$  real(gp) :: prefac,xsum,stot,fsum,tt,firstprod
!!$  !build the prefactor
!!$  prefac=sqrt(a)
!!$  prefac=c**l/prefac
!!$  prefac=gammaonehalf*prefac
!!$
!!$  !object in the sum
!!$  xsum=a*c**2
!!$  xsum=1.d0/xsum
!!$
!!$  !the first term of the sum is one
!!$  stot=1.d0
!!$
!!$  !calculate the sum
!!$  do p=1,l/4
!!$     tt=real(ifac(p+1,2*p),gp)
!!$     fsum=real(ifac(l-2*p+1,l),gp)
!!$     fsum=fsum/tt
!!$     tt=firstprod(p)
!!$     fsum=fsum*tt
!!$     fsum=fsum*xsum**p
!!$     stot=stot+fsum
!!$  end do
!!$  do p=l/4+1,l/3
!!$     tt=real(ifac(p+1,l-2*p),gp)
!!$     fsum=real(ifac(2*p+1,l),gp)
!!$     fsum=fsum/tt
!!$     tt=firstprod(p)
!!$     fsum=fsum*tt
!!$     fsum=fsum*xsum**p
!!$     stot=stot+fsum
!!$  end do
!!$  do p=l/3+1,l/2
!!$     tt=real(ifac(l-2*p+1,p),gp)
!!$     fsum=real(ifac(2*p+1,l),gp)
!!$     fsum=fsum*tt
!!$     tt=firstprod(p)
!!$     fsum=fsum*tt
!!$     fsum=fsum*xsum**p
!!$     stot=stot+fsum
!!$  end do
!!$
!!$  !final result
!!$  gauinti=stot*prefac
!!$  
!!$END FUNCTION gauinti
!!$
!!$
!!$
!!$!>   Valid if p<l/4 AND p/=0
!!$!!
!!$!!
!!$function secondprod1(p,l)
!!$  use module_base
!!$  implicit none
!!$  integer, intent(in) :: p,l
!!$  real(gp) :: secondprod1
!!$  !local variables
!!$  !integer :: i
!!$  real(gp) :: tt,part1,rfac
!!$  part1=rfac(p+1,2*p)
!!$  !divide by the last value
!!$  part1=real(l,gp)/part1
!!$  tt=rfac(l-2*p+1,l-1)
!!$!!!  part1=1.d0
!!$!!!  do i=p+1,2*p !in the second case the bound must be changed here
!!$!!!     tt=real(i,gp)
!!$!!!     part1=part1*tt
!!$!!!  end do
!!$  secondprod1=tt*part1
!!$END FUNCTION secondprod1
!!$
!!$
!!$
!!$!>   Valid if p>=l/4 AND p<l/3
!!$!!
!!$!!
!!$function secondprod2(p,l)
!!$  use module_base
!!$  implicit none
!!$  integer, intent(in) :: p,l
!!$  real(gp) :: secondprod2
!!$  !local variables
!!$  !integer :: i
!!$  real(gp) :: tt,part1,rfac
!!$  part1=rfac(p+1,l-2*p)
!!$  !divide by the last value
!!$  part1=real(l,gp)/part1
!!$  tt=rfac(2*p+1,l-1)
!!$!!!  part1=1.d0
!!$!!!  do i=p+1,2*p !in the second case the bound must be changed here
!!$!!!     tt=real(i,gp)
!!$!!!     part1=part1*tt
!!$!!!  end do
!!$  secondprod2=tt*part1
!!$END FUNCTION secondprod2
!!$
!!$
!!$
!!$!>   Integer version of factorial
!!$!!
!!$!!
!!$function ifac(is,ie)
!!$  implicit none
!!$  integer, intent(in) :: is,ie
!!$  integer :: ifac
!!$  !local variables
!!$  integer :: i
!!$  ifac=1
!!$  do i=is,ie
!!$     ifac=ifac*i
!!$  end do
!!$END FUNCTION ifac
!!$


!>   Calculate the projection of norb wavefunctions on a gaussian basis set
!!
!!
subroutine wavelets_to_gaussians(geocode,norbp,nspinor,n1,n2,n3,G,thetaphi,hx,hy,hz,wfd,psi,coeffs)
  use module_base
  use module_types
  implicit none
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
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
             G%nshell,G%ndoc,G%nam,G%xp(1,:),G%psiat(1,:),G%nshltot,G%nexpo,G%ncoeff,&
             hx,hy,hz,wfd,psi(1,ispinor,iorb),coeffs(1,ispinor,iorb))
        !print *,'iorb, coeffs',iorb,coeffs(:,1,iorb)
     end do
  end do
  
END SUBROUTINE wavelets_to_gaussians



!>   Calculate the projection of a given orbital on a gaussian basis centered on 
!!   a set of points
!!
!!
subroutine orbital_projection(geocode,n1,n2,n3,nat,rxyz,thetaphi,nshell,ndoc,nam,xp,psiat,&
     nshltot,nexpo,ncoeff,hx,hy,hz,wfd,psi,coeffs)
  use module_base
  use module_types
  use gaussians
  implicit none
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
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


subroutine dual_gaussian_coefficients(norbp,G,coeffs)
  use module_base
  use gaussians
  implicit none
  integer, intent(in) :: norbp
  type(gaussian_basis), intent(in) :: G
  real(gp), dimension(G%ncoeff,norbp), intent(inout) :: coeffs !warning: the precision here should be wp
  !local variables
  character(len=*), parameter :: subname='dual_gaussian_coefficients'
  integer :: nwork,info
  integer, dimension(:), allocatable :: iwork
  real(gp), dimension(:), allocatable :: ovrlp,work

  iwork = f_malloc(G%ncoeff,id='iwork')
  ovrlp = f_malloc(G%ncoeff*G%ncoeff,id='ovrlp')

  !temporary allocation of the work array, workspace query in dsysv
  work = f_malloc(100,id='work')

  if (norbp > 0) then
     call dsysv('U',G%ncoeff,norbp,ovrlp(1),G%ncoeff,iwork(1),coeffs(1,1),&
          G%ncoeff,work(1),-1,info)
  end if
  nwork=int(work(1))

  call f_free(work)
  work = f_malloc(nwork,id='work')

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

  call f_free(iwork)
  call f_free(work)
  call f_free(ovrlp)

!!!  do iorb=1,norbp
!!!     print *,'iorb, dual,coeffs',iorb,coeffs(:,iorb)
!!!  end do

END SUBROUTINE dual_gaussian_coefficients



!>   Calculate the projection of a gaussian for a given eigenspace of spherical harmonics
!!   centered on a given point and rotated by theta(along z) and phi(along x)
!!
!!
subroutine lsh_projection(geocode,l,ng,xp,psiat,n1,n2,n3,rxyz,thetaphi,hx,hy,hz,&
     wfd,psi,coeffs)
  use module_base
  use module_types
  implicit none
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
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
        & rxyz(1),rxyz(2),rxyz(3),hx,hy,hz,wfd%nseg_c,wfd%nvctr_c,wfd%keygloc,wfd%keyvloc,&
        & wfd%nseg_f,wfd%nvctr_f,&
        & wfd%keygloc(1,wfd%nseg_c+min(1,wfd%nseg_f)),&
        & wfd%keyvloc(wfd%nseg_c+min(1,wfd%nseg_f)),&
        & psi(1),psi(wfd%nvctr_c+min(1,wfd%nvctr_f)),coeffs(m))

     !print '(a,2(i4),5(1pe12.5))','l,m,rxyz,coeffs(m)',l,m,rxyz(:),coeffs(m)
  end do

  !here we should modify the coefficients following the direction of the axis
  call lsh_rotation(l-1,thetaphi(1),thetaphi(2),coeffs)

END SUBROUTINE lsh_projection


!>
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


!>   Calculate the scalar product between a sum of gaussians times polynomials and a wavefunction
!!   @f$\int dx dy dz 
!!             \sum_i=1..ntp fac_arr(i) {
!!                  \sum_j=1..nterm psiat(j) [exp(-r^2/(2*(xp(j)^2)))] 
!!                      *((x-rx)^lx(i) *(y-ry)^ly(i) * (z-rz)^lz(i) ))} psi(x,y,z)@f$
!!   Expressed in Daubechies Basis in the arrays psi_c, psi_f
subroutine wavetogau(geocode,n1,n2,n3,nterm,ntp,lx,ly,lz,fac_arr,xp,psiat,rx,ry,rz,hx,hy,hz, & 
     nseg_c,mvctr_c,keyg_c,keyv_c,nseg_f,mvctr_f,keyg_f,keyv_f,psi_c,psi_f,overlap)
  use module_base
  implicit none
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
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
  integer:: iterm,itp,n_gau,ml1,mu1,ml2,mu2,ml3,mu3,i1,i2,i3,iseg,ii,jj,j0,j1,i0,i
  real(wp) :: ovlp_c,ovlp_f1,ovlp_f2,ovlp_f3,ovlp_f4,ovlp_f5,ovlp_f6,ovlp_f7,ovlp
  real(gp) :: gau_a,te
  real(wp), dimension(0:nw,2) :: work
  real(wp), dimension(:,:), allocatable :: wprojx,wprojy,wprojz

  !conditions for periodicity in the three directions
  perx=(geocode /= 'F')
  pery=(geocode == 'P')
  perz=(geocode /= 'F')

  wprojx = f_malloc((/ 0.to.n1, 1.to.2 /),id='wprojx')
  wprojy = f_malloc((/ 0.to.n2, 1.to.2 /),id='wprojy')
  wprojz = f_malloc((/ 0.to.n3, 1.to.2 /),id='wprojz')

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

  call f_free(wprojx)
  call f_free(wprojy)
  call f_free(wprojz)

END SUBROUTINE wavetogau


!>   Coefficients of the rotation matrix
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

