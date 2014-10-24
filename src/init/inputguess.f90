!> @file
!!  Routines to generate the input guess
!! @author
!!    Copyright (C) 2007-2013 (LG) BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!> Generate the input guess via the inguess_generator
subroutine inputguess_gaussian_orbitals(iproc,nproc,at,rxyz,nvirt,nspin,&
      &   orbs,orbse,norbsc_arr,locrad,G,psigau,eks,iversion,mapping,quartic_prefactor)
   use module_base
   use module_types
   use module_interfaces, except_this_one => inputguess_gaussian_orbitals
   use yaml_output
   implicit none
   integer, intent(in) :: iproc,nproc,nspin
   integer, intent(inout) :: nvirt
   type(atoms_data), intent(in) :: at
   type(orbitals_data), intent(in) :: orbs
   real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
   real(gp), intent(out) :: eks
   integer, dimension(at%natsc+1,nspin), intent(out) :: norbsc_arr
   real(gp), dimension(at%astruct%nat), intent(out) :: locrad
   type(orbitals_data), intent(out) :: orbse
   type(gaussian_basis), intent(out) :: G
   real(wp), dimension(:,:,:), pointer :: psigau
   integer,intent(in) :: iversion !< 1:cubic, 2:linear
   integer,dimension(orbs%norb),intent(in),optional:: mapping
   real(gp),dimension(at%astruct%ntypes),intent(in),optional:: quartic_prefactor
   !local variables
   character(len=*), parameter :: subname='inputguess_gaussian_orbitals'
   !n(c) integer, parameter :: ngx=31
   integer :: norbe,norbme,norbyou,norbsc,nvirte,ikpt
   integer :: ispin,jproc,ist,jpst,nspinorfororbse,noncoll
   integer, dimension(:), allocatable :: iorbtolr

   !Generate the input guess via the inguess_generator
   !here we should allocate the gaussian basis descriptors 
   !the prescriptions can be found in the creation of psp basis
   call readAtomicOrbitals(at,norbe,norbsc,nspin,orbs%nspinor,&
        norbsc_arr,locrad)

   !in the non-collinear case the number of orbitals double
   if (orbs%nspinor == 4) then
      noncoll=2
   else
      noncoll=1
   end if

   if (iproc ==0) then
      call yaml_map('Total No. of Atomic Input Orbitals',nspin*noncoll*norbe,fmt='(i6)')
      !write(*,'(1x,a,i0,a)')'Generating ',nspin*noncoll*norbe,' Atomic Input Orbitals'
      if (norbsc /=0) then
         !write(*,'(1x,a,i0,a)')'  of which ',nspin*noncoll*norbsc,&
         !&   ' are semicore orbitals'
         call yaml_map('No. of Semicore Orbitals',nspin*noncoll*norbsc,fmt='(i6)')
      end if
   end if

   if (nvirt /= 0) then
      do ispin=1,nspin
         !Check for max number of virtual orbitals
         !the unoccupied orbitals available as a LCAO
         !this is well defined only for closed-shell systems
         !if (ispin == 1) nvirte=min(noncoll*norbe-orbs%norbu,nvirt)
         !if (ispin == 2) nvirte=min(noncoll*norbe-orbs%norbd,nvirt)
         !alternative test, put always the limit to the number of elements of the input guess
         nvirte=noncoll*norbe
         if(nvirt < nvirte .and. iproc==0) then
            call yaml_warning('A bigger number of virtual orbitals may be needed for better convergence')
            call yaml_map('Suggested nvirt',nvirte,fmt='(i0)')
!!$            write(*,'(1x,a)')&
!!$               &   "WARNING: A bigger number of virtual orbitals may be needed for better convergence."
!!$            write(*,'(1x,a,i0)')'         Put nvirt= ',nvirte
         end if
      end do
   end if

   !!!orbitals descriptor for inguess orbitals
   nspinorfororbse=orbs%nspinor

   !the number of orbitals to be considered is doubled 
   !in the case of a spin-polarised calculation
   !also for non-collinear case
   !nspin*noncoll is always <= 2
   select case(iversion)
   case(1)
       ! SM: use for basedistd the same basedist as for the up orbitals. CHECK THIS!!!
       call orbitals_descriptors(iproc,nproc,nspin*noncoll*norbe,noncoll*norbe,(nspin-1)*norbe, &
            nspin,nspinorfororbse,orbs%nkpts,orbs%kpts,orbs%kwgts,orbse,LINEAR_PARTITION_NONE,&
            basedist=orbs%norb_par(0:,1:),basedistu=orbs%norbu_par(0:,1:),basedistd=orbs%norbu_par(0:,1:))
   case(2)
       call orbitals_descriptors(iproc,nproc,nspin*noncoll*norbe,noncoll*norbe,(nspin-1)*norbe, &
            nspin,nspinorfororbse,orbs%nkpts,orbs%kpts,orbs%kwgts,orbse,LINEAR_PARTITION_SIMPLE)
   case default
       call f_err_throw('wrong value of iversion',err_id=BIGDFT_RUNTIME_ERROR)
   end select
   do ikpt = 1, orbse%nkpts
      ist=1 + (ikpt - 1 ) * nspin*noncoll*norbe
      do ispin=1,nspin
         orbse%spinsgn(ist:ist+norbe-1)=real(1-2*(ispin-1),gp)
         ist=ist+norbe
      end do
   end do

   !this is the distribution procedure for cubic code
   !should be referred to another routine
   if (iproc == 0 .and. nproc > 1) then
      call yaml_newline()
      call yaml_mapping_open('Inputguess Orbitals Repartition')
      jpst=0
      do jproc=0,nproc-1
         norbme=orbse%norb_par(jproc,0)
         norbyou=orbse%norb_par(min(jproc+1,nproc-1),0)
         if (norbme /= norbyou .or. jproc == nproc-1) then
            call yaml_map('MPI tasks '//trim(yaml_toa(jpst,fmt='(i0)'))//'-'//trim(yaml_toa(jproc,fmt='(i0)')),norbme,fmt='(i0)')
            !!this is a screen output that must be modified
            !write(*,'(3(a,i0),a)')&
            !   &   ' Processes from ',jpst,' to ',jproc,' treat ',norbme,' inguess orbitals '
            jpst=jproc+1
         end if
      end do
      call yaml_mapping_close()
      !write(*,'(3(a,i0),a)')&
         !     ' Processes from ',jpst,' to ',nproc-1,' treat ',norbyou,' inguess orbitals '
   end if

  !write(*,'(a,3i6)') 'iproc, orbse%isorb, orbse%norbp', iproc, orbse%isorb,orbse%norbp
  !write(*,'(a,3i6)') 'norbe, orbse%nspinor, orbse%isorb+orbse%norbp+ndebug', norbe, orbse%nspinor, orbse%isorb+orbse%norbp+ndebug
   !allocate the gaussian coefficients for the number of orbitals which is needed
   psigau = f_malloc_ptr((/ norbe , orbse%nspinor , max(orbse%norbp,1) /),id='psigau')
   iorbtolr = f_malloc(orbse%norbp,id='iorbtolr')

   !fill just the interesting part of the orbital
   if (present(mapping)) then
       ! this will be use for the linear scaling part
       if(present(quartic_prefactor)) then
           call AtomicOrbitals(iproc,at,rxyz,norbe,orbse,norbsc,nspin,eks,G,&
                psigau(1,1,1),&
                iorbtolr,mapping,quartic_prefactor)
       else
           call AtomicOrbitals(iproc,at,rxyz,norbe,orbse,norbsc,nspin,eks,G,&
                psigau(1,1,1),&
                iorbtolr,mapping)
       end if
   else
       call AtomicOrbitals(iproc,at,rxyz,norbe,orbse,norbsc,nspin,eks,G,&
            psigau(1,1,1),iorbtolr)
   end if

   call f_free(iorbtolr)


END SUBROUTINE inputguess_gaussian_orbitals

!> Read atomic orbitals
subroutine readAtomicOrbitals(at,norbe,norbsc,nspin,nspinor,norbsc_arr,locrad)
   use module_base, only: gp
   use ao_inguess, only: atomic_info,ao_nspin_ig,count_atomic_shells
   use module_types
   implicit none
   !Arguments
   integer, intent(in) :: nspin,nspinor
   integer, intent(out) :: norbe,norbsc
   type(atoms_data), intent(in) :: at
   !logical, dimension(4,2,at%natsc), intent(out) :: scorb
   integer, dimension(at%natsc+1,nspin), intent(out) :: norbsc_arr
   real(gp), dimension(at%astruct%nat), intent(out) :: locrad
   !local variables
   !n(c) character(len=*), parameter :: subname='readAtomicOrbitals'
   !integer, parameter :: lmax=3,noccmax=2,nelecmax=32
   !character(len=2) :: symbol
   integer :: ity,iatsc,iat !,i,lsc
   !integer :: nsccode,mxpl,mxchg
   integer :: norbat,iorbsc_count !,niasc,nlsc
   real(gp) :: ehomo!rcov,rprb,ehomo,amu
   !integer, dimension(nmax,lmax+1) :: neleconf
   !real(kind=8), dimension(nmax,lmax+1) :: neleconf
   !integer, dimension(lmax+1) :: nl
   !real(gp), dimension(noccmax,lmax+1) :: occup

   ! number of orbitals, total and semicore
   norbe=0
   norbsc=0
   iatsc=0
   !scorb(:,:,:)=.false.
   do iat=1,at%astruct%nat
      ity=at%astruct%iatype(iat)
      !call count_atomic_shells(ao_nspin_ig(nspin,nspinor=nspinor),&
      !     at%aoig(iat)%aocc,occup,nl)

      norbat=at%aoig(iat)%nao!(nl(1)+3*nl(2)+5*nl(3)+7*nl(4))

      norbe=norbe+norbat
      !print *,'iat',iat,l,norbe,norbat,nl(:)
      !calculate the localisation radius for the input orbitals 
      !call eleconf(at%nzatom(ity),at%nelpsp(ity),symbol,rcov,rprb,ehomo,&
      !   &   neleconf,nsccode,mxpl,mxchg,amu)
      call atomic_info(at%nzatom(ity),at%nelpsp(ity),ehomo=ehomo)

      locrad(iat)=5._gp/sqrt(abs(2._gp*ehomo))
      !count the number of semicore orbitals we have for this atom
      iorbsc_count=at%aoig(iat)%nao_sc
      if (iorbsc_count /= 0) then
         iatsc=iatsc+1
         norbsc_arr(iatsc,1)=iorbsc_count
         norbsc=norbsc+iorbsc_count
      end if

!!$      nsccode=at%aoig(iat)%iasctype
!!$      if (nsccode/=0) then !the atom has some semicore orbitals
!!$         iatsc=iatsc+1
!!$         niasc=nsccode
!!$         !count the semicore orbitals for this atom
!!$         iorbsc_count=0
!!$         do lsc=4,1,-1
!!$            nlsc=niasc/4**(lsc-1)
!!$            iorbsc_count=iorbsc_count+nlsc*(2*lsc-1)
!!$            if (nlsc > 2) then
!!$               write(*,*)'ERROR, atom:',iat,&
!!$                  &   ': cannot admit more than two semicore shells per channel',nlsc
!!$               stop
!!$            end if
!!$            do i=1,nlsc
!!$               scorb(lsc,i,iatsc)=.true.
!!$            end do
!!$            niasc=niasc-nlsc*4**(lsc-1)
!!$         end do
!!$         norbsc_arr(iatsc,1)=iorbsc_count
!!$         norbsc=norbsc+iorbsc_count
!!$         !if (iproc == 0) write(*,*) iat,nsccode,iorbsc_count,norbsc,scorb(:,:,iatsc)
!!$      end if

   end do

   !print *,'NL',nl,norbe

   !orbitals which are non semicore
   norbsc_arr(at%natsc+1,1)=norbe-norbsc

   !duplicate the values in the case of spin-polarization
   if (nspin == 2) norbsc_arr(:,2)=norbsc_arr(:,1)

END SUBROUTINE readAtomicOrbitals

!> Generate atomic orbitals
subroutine AtomicOrbitals(iproc,at,rxyz,norbe,orbse,norbsc,&
      &   nspin,eks,G,gaucoeff,iorbtolr,mapping,quartic_prefactor)
   use module_base
   use ao_inguess, only: iguess_generator,print_eleconf,ao_nspin_ig,count_atomic_shells,&
        nmax_occ => nmax_occ_ao
   use module_types
   use module_interfaces, except_this_one => AtomicOrbitals
   use yaml_output
   implicit none
   integer, intent(in) :: norbe,iproc
   integer, intent(in) :: norbsc,nspin
   type(atoms_data), intent(in) :: at
   real(gp), dimension(3,at%astruct%nat), intent(in), target :: rxyz
   type(orbitals_data), intent(inout) :: orbse
   type(gaussian_basis), intent(out) :: G
   real(gp), intent(out) :: eks
   integer, dimension(orbse%norbp), intent(out) :: iorbtolr !assign the localisation region
   real(wp), dimension(norbe,orbse%nspinor,orbse%norbp), intent(out) :: gaucoeff !norbe=G%ncoeff
   integer,dimension(orbse%norb), optional, intent(in):: mapping
   real(gp),dimension(at%astruct%ntypes),intent(in),optional:: quartic_prefactor
   !local variables
   character(len=*), parameter :: subname= 'AtomicOrbitals'
   !integer, parameter :: noccmax=2,lmax=4,nelecmax=32,nmax_occ=10!actually is 24
   !integer, parameter :: nterm_max=3,nmax=7
   logical :: orbpol_nc,occeq
   integer :: iatsc,i_stat,ispin,iexpo,ishltmp,ngv,ngc,islcc,iiorb,jjorb
   integer :: iorb,jorb,iat,ity,i,ictot,inl,l,m,nctot,iocc,ictotpsi,ishell,icoeff
   integer :: noncoll,ig,ispinor,icoll,ikpts,ikorb,nlo,ntypesx,ityx,jat,ng,nspin_print
   !integer :: nsccode
   real(gp) :: ek,mx,my,mz,ma,mb,mc,md
   real(gp) :: mnorm,fac
   !logical, dimension(lmax,noccmax) :: semicore
   integer, dimension(2) :: iorbsc,iorbv
   !integer, dimension(lmax) :: nl
   !real(gp), dimension(noccmax,lmax) :: occup
   integer, dimension(:), allocatable :: iatypex
   real(gp), dimension(:), allocatable :: psiatn
   real(gp), dimension(:,:), allocatable :: atmoments,xp
   real(gp), dimension(:,:,:), allocatable :: psiat

   !if (iproc == 0 .and. verbose > 1) then
      !write(*,'(1x,a)')'Calculating AIO wavefunctions: '
   !end if

   !gaussian basis structure informations
   !insert these things in the loops above
   !here we can create a gaussian basis structure for the input guess functions
   !the number of gaussian centers are thus nat
   G%nat=at%astruct%nat
   G%rxyz => rxyz
   !copy the parsed values in the gaussian structure
   !count also the total number of shells
   G%nshell = f_malloc_ptr(at%astruct%nat,id='G%nshell')

   !if non-collinear it is like nspin=1 but with the double of orbitals
   if (orbse%nspinor == 4) then
      noncoll=2
      nspin_print=4
   else
      noncoll=1
      nspin_print=nspin
   end if

   !calculate the number of atom types by taking into account the occupation
   iatypex = f_malloc(at%astruct%nat,id='iatypex')

   ntypesx=0
   G%nshltot=0
   count_shells: do iat=1,at%astruct%nat
      !print *,'atom,aocc',iat,at%aocc(1:nelecmax,iat)
      ity=at%astruct%iatype(iat)
      !call count_atomic_shells(nspin_print,at%aoig(iat)%aocc,occup,nl)
      G%nshell(iat)=sum(at%aoig(iat)%nl)!(nl(1)+nl(2)+nl(3)+nl(4))
      G%nshltot=G%nshltot+G%nshell(iat)
      !check the occupation numbers and the atoms type
      !once you find something equal exit the procedure
      do jat=1,iat-1
         if (at%astruct%iatype(jat) == ity) then
            occeq=all(at%aoig(jat)%aocc == at%aoig(iat)%aocc)
!!$            do i=1,nelecmax
!!$               occeq = occeq .and. &
!!$                    (at%aoig(jat)%aocc(i) == at%aoig(iat)%aocc(i))
!!$            end do
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

   !print *,'ntypesx',ntypesx,iatypex

   G%ndoc = f_malloc_ptr(G%nshltot,id='G%ndoc')
   G%nam = f_malloc_ptr(G%nshltot,id='G%nam')

   !the default value for the gaussians is chosen to be 21
   ng=21
   !allocate arrays for the inequivalent wavefunctions
   xp = f_malloc((/ ng, ntypesx /),id='xp')
   psiat = f_malloc((/ ng, nmax_occ, ntypesx /),id='psiat')

   !print *,'atomx types',ntypesx

   if (iproc == 0 .and. verbose > 1) then
      call yaml_newline()
      call yaml_sequence_open('Atomic Input Orbital Generation')
      call yaml_newline()
   end if


   !assign shell IDs and count the number of exponents and coefficients
   !also calculate the wavefunctions 
   G%ncplx=1 !2 only for PAW and projectors
   G%nexpo=0
   G%ncoeff=0
   ishell=0
   ntypesx=0
   do iat=1,at%astruct%nat
      ity=at%astruct%iatype(iat)
      ityx=iatypex(iat)
      ishltmp=0
      !call count_atomic_shells(nspin_print,at%aoig(iat)%aocc,occup,nl)
      if (ityx > ntypesx) then
         if (iproc == 0 .and. verbose > 1) then
            call yaml_sequence(advance='no')
            call yaml_mapping_open(flow=.true.)
            call yaml_map('Atom Type',trim(at%astruct%atomnames(ity)))
            call print_eleconf(nspin_print,&
                 at%aoig(iat)%aocc,at%aoig(iat)%nl_sc)
         end if

         !positions for the nlcc arrays
         call nlcc_start_position(ity,at,ngv,ngc,islcc)
         !print *,'debug',ity,ngv,ngc,islcc,at%nlccpar(:,:),'acc',shape(at%nlccpar),'end'
         !eliminate the nlcc parameters from the IG, since XC is always LDA
         ngv=0
         ngc=0
         if(present(quartic_prefactor)) then
             call iguess_generator(at%nzatom(ity),at%nelpsp(ity),&
                  real(at%nelpsp(ity),gp),nspin_print,at%aoig(iat)%aocc,at%psppar(0:,0:,ity),&
                  at%npspcode(ity),ngv,ngc,at%nlccpar(0:,max(islcc,1)),&
                  ng-1,xp(1,ityx),psiat(1:,1:,ityx),.false.,&
                  quartic_prefactor=quartic_prefactor(ity))
        else
             call iguess_generator(at%nzatom(ity),at%nelpsp(ity),&
                  real(at%nelpsp(ity),gp),nspin_print,at%aoig(iat)%aocc,at%psppar(0:,0:,ity),&
                  at%npspcode(ity),ngv,ngc,at%nlccpar(0:,max(islcc,1)),&
                  ng-1,xp(1,ityx),psiat(1:,1:,ityx),.false.)
         end if
         ntypesx=ntypesx+1
         if (iproc == 0 .and. verbose > 1) then
            !write(*,'(1x,a)')'done.'
            call yaml_mapping_close()
         end if
      end if

      do l=1,4
         do i=1,at%aoig(iat)%nl(l-1)
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
         call f_err_throw('Input Guess ishelltmp (' // trim(yaml_toa(ishell)) // ') /= nshell (' &
             & // trim(yaml_toa(G%nshell(iat))) // ')', err_id=BIGDFT_INPUT_VARIABLES_ERROR)
         !write(*,*)'ERROR: ishelltmp <> nshell',ishell,G%nshell(iat)
         !stop 
      end if
   end do
   if (iproc == 0 .and. verbose > 1) then
      call yaml_sequence_close()
      call yaml_newline()
   end if

   !print *,'nl',nl,norbe,G%ncoeff
   if (norbe /= G%ncoeff) then
         call f_err_throw('Input Guess norbe(' // trim(yaml_toa(norbe)) // ') /= G%ncoeff (' &
             & // trim(yaml_toa(G%ncoeff)) // ')', err_id=BIGDFT_INPUT_VARIABLES_ERROR)
      !write(*,*)'ERROR: norbe /= G%ncoeff',norbe,G%ncoeff
      !stop 
   end if

   call to_zero(orbse%norbp*orbse%nspinor*G%ncoeff,gaucoeff)

   !allocate and assign the exponents and the coefficients
   G%psiat = f_malloc_ptr((/ G%ncplx, G%nexpo /),id='G%psiat')
   G%xp = f_malloc_ptr((/ G%ncplx, G%nexpo /),id='G%xp')
   psiatn = f_malloc(ng,id='psiatn')

   !read the atomic moments if non-collinear
   !WARNING: units are not good for the moment
   !the moments can be inserted in the atoms_data structure
   if (orbse%nspinor == 4) then
      atmoments = f_malloc((/ 3, at%astruct%nat /),id='atmoments')

      open(unit=22,file='moments',form='formatted',action='read',status='old')
      !this part can be transferred on the atomic orbitals section
      do iat=1,at%astruct%nat
         read(unit=22,fmt=*,iostat=i_stat) mx,my,mz
         if (i_stat > 0) then
            call f_err_throw('The file "moments" is not correct!' // &
               & 'The file "moments" has the line ' // trim(yaml_toa(iat)) // &
               & ' which have not 3 numbers for the atom ' // trim(yaml_toa(iat)) // '.', &
               & err_id=BIGDFT_INPUT_VARIABLES_ERROR)
            !write(unit=*,fmt='(a,i0,a,i0,a)') 'The file "moments" has the line ',iat,&
            !   &   ' which have not 3 numbers for the atom ',iat,'.'
            !stop 'The file "moments" is not correct!'
         end if
         atmoments(1,iat)=mx
         atmoments(2,iat)=my
         atmoments(3,iat)=mz
      end do
   end if

   eks=0.0_gp
   iorb=0
   iatsc=0

   !initialise the orbital counters
   iorbsc(1)=0
   iorbv(1)=norbsc
   !used in case of spin-polarisation, ignored otherwise
   iorbsc(2)=norbe
   iorbv(2)=norbsc+norbe

   ishell=0
   iexpo=0
   icoeff=1
   do iat=1,at%astruct%nat

      ity=at%astruct%iatype(iat)
      ityx=iatypex(iat)
      !call count_atomic_shells(nspin_print,at%aoig(iat)%aocc,occup,nl)

      !the scorb array was already corrected in readAtomicOrbitals routine
      !nsccode=at%aoig(iat)%iasctype
      !if (nsccode/=0) then !the atom has some semicore orbitals
      if (at%aoig(iat)%nao_sc/=0) then !the atom has some semicore orbitals
         iatsc=iatsc+1
!!$         semicore(:,:)=scorb(:,:,iatsc)
      else
!!$         semicore(:,:)=.false.
      end if

      !calculate the atomic input orbitals
      ictot=0
      ictotpsi=0
      nctot=sum(at%aoig(iat)%nl)!(nl(1)+nl(2)+nl(3)+nl(4))
      if (iorbsc(1)+nctot > norbe .and. iorbv(1)+nctot > norbe) then
         call f_err_throw('transgpw occupe ' // trim(yaml_toa(at%aoig(iat)%nl)) // ' ' // &
           & trim(yaml_toa(norbe)), err_id=BIGDFT_RUNTIME_ERROR)
         stop
         !print *,'transgpw occupe',at%aoig(iat)%nl,norbe
         !stop
      end if
      iocc=0
      do l=1,4
         iocc=iocc+1
         nlo=ceiling(at%aoig(iat)%aocc(iocc))!nint(at%aocc(iocc,iat))
         do inl=1,nlo
            ictotpsi=ictotpsi+1
            ictot=ictot+1
            ishell=ishell+1

            !contribution to the kinetic energy given by the electrons in this shell
            call atomkin(l-1,ng,xp(1,ityx),psiat(1,ictotpsi,ityx),psiatn,ek)
            do ig=1,G%ndoc(ishell)
               iexpo=iexpo+1
               G%psiat(1,iexpo)=psiatn(ig)
               G%xp(1,iexpo)=xp(ig,ityx)
            end do

            do ispin=1,nspin
               !in the new language the semicore orbitals 
               !should be  if (inl <= nsc(l))
               !the order of the orbitals (iorb,jorb) must put in the beginning
               !the semicore orbitals
               !if (semicore(l,inl)) then
               if (inl <= at%aoig(iat)%nl_sc(l-1)) then
                  !the orbital is semi-core
                  iorb=iorbsc(ispin)
               else
                  !normal case, the orbital is a valence orbital
                  iorb=iorbv(ispin)
               end if

               do m=1,2*l-1
                  !each orbital has two electrons in the case of the 
                  !non-collinear case
                  do icoll=1,noncoll !non-trivial only for nspinor=4
                     iocc=iocc+1
                     iorb=iorb+1
                     !in noncollinear case see if the orbital is polarised
                     if (noncoll==2) then
                        orbpol_nc=.true.
                        !full orbital, non polarised
                        !if icoll=1 occ(iocc)+occ(iocc+1)
                        !if icoll=2 occ(iocc-1)+occ(iocc)
                        if (at%aoig(iat)%aocc(iocc+1-icoll)+&
                             at%aoig(iat)%aocc(iocc+2-icoll) == 2.0_gp) then
                           orbpol_nc=.false.
                        end if
                     end if
                     do ikpts=1,orbse%nkpts
                        ikorb=iorb+(ikpts-1)*orbse%norb
                        jorb=ikorb-orbse%isorb
                        orbse%occup(ikorb)=at%aoig(iat)%aocc(iocc)

                        !!write(*,'(a,4i6,2es12.4)') 'iat, l, m, iocc, at%aocc(iocc,iat), ek', iat, l, m, iocc, at%aocc(iocc,iat), ek
                        eks=eks+ek*at%aoig(iat)%aocc(iocc)*orbse%kwgts(ikpts)
                        !!write(*,*) 'iat, at%aocc(iocc,iat)', iat, at%aocc(iocc,iat)
                        if (present(mapping)) then
                           iiorb=mapping(iorb)
                           jjorb=iiorb-orbse%isorb
                        else
                           iiorb=ikorb
                           jjorb=jorb
                        end if
                        !print *,'iat',iat,l,m,icoll,orbse%occup(ikorb),orbpol_nc
                        !if (orbse%isorb < ikorb .and. ikorb <= orbse%isorb+orbse%norbp) then
                        if (orbse%isorb < iiorb .and. iiorb <= orbse%isorb+orbse%norbp) then
                           if (orbse%nspinor == 1) then
                              do ispinor=1,orbse%nspinor
                                 !here we put only the case nspinor==1
                                 !we can put a phase for check with the complex wavefunction
                                 gaucoeff(icoeff,ispinor,jjorb)=1.0_wp
                              end do
                           else if (orbse%nspinor == 2) then
                              gaucoeff(icoeff,1,jjorb)=1.0_wp
                              gaucoeff(icoeff,2,jjorb)=0.0_wp
                              !!$                             !we can put a phase for check with the complex wavefunction
                              !!$                             gaucoeff(icoeff,1,jorb)=0.5_wp*sqrt(3.0_wp)
                              !!$                             gaucoeff(icoeff,2,jorb)=0.5_wp
                           else if (orbse%nspinor == 4) then
                              !assign the input orbitals according to the atomic moments
                              fac=0.5_gp

                              !if the shell is not polarised
                              !put one electron up and the other down
                              !otherwise, put a small unbalancing on the orbitals
                              !such that the momentum of the
                              if (orbpol_nc) then
                                 !in the case of unoccupied orbital, 
                                 !choose the orthogonal direction
                                 mx=atmoments(1,iat)
                                 my=atmoments(2,iat)
                                 mz=atmoments(3,iat)

                                 if (orbse%occup(ikorb) == 0.0_gp) then
                                    mx=-mx
                                    my=-my
                                    mz=-mz
                                 end if
                              else
                                 mx=0.0_gp
                                 my=0.0_gp
                                 mz=1.0_gp-2.0_gp*real(icoll-1,gp)
                              end if

                              mnorm=sqrt(mx**2+my**2+mz**2)
                              if (mnorm /= 0.0_gp) then
                                 mx=mx/mnorm
                                 my=my/mnorm
                                 mz=mz/mnorm
                              end if

                              ma=0.0_gp
                              mb=0.0_gp
                              mc=0.0_gp
                              md=0.0_gp

                              if(mz > 0.0_gp) then 
                                 ma=ma+mz
                              else
                                 mc=mc+abs(mz)
                              end if
                              if(mx > 0.0_gp) then 
                                 ma=ma+fac*mx
                                 mb=mb+fac*mx
                                 mc=mc+fac*mx
                                 md=md+fac*mx
                              else
                                 ma=ma-fac*abs(mx)
                                 mb=mb-fac*abs(mx)
                                 mc=mc+fac*abs(mx)
                                 md=md+fac*abs(mx)
                              end if
                              if(my > 0.0_gp) then 
                                 ma=ma+fac*my
                                 mb=mb-fac*my
                                 mc=mc+fac*my
                                 md=md+fac*my
                              else
                                 ma=ma-fac*abs(my)
                                 mb=mb+fac*abs(my)
                                 mc=mc+fac*abs(my)
                                 md=md+fac*abs(my)
                              end if
                              if(mx==0.0_gp .and. my==0.0_gp .and. mz==0.0_gp) then
                                 ma=1.0_gp/sqrt(2.0_gp)
                                 mb=0.0_gp
                                 mc=1.0_gp/sqrt(2.0_gp)
                                 md=0.0_gp
                              end if

                              !assign the gaussian coefficients for each
                              !spinorial component
                              gaucoeff(icoeff,1,jorb)=real(ma,wp)
                              gaucoeff(icoeff,2,jorb)=real(mb,wp)
                              gaucoeff(icoeff,3,jorb)=real(mc,wp)
                              gaucoeff(icoeff,4,jorb)=real(md,wp)
                           end if

                           !associate to each orbital the reference localisation region
                           !here identified by the atom
                           iorbtolr(jjorb)=iat 
                        endif
                     end do
                  end do
                  icoeff=icoeff+1
               end do
               !if (semicore(l,inl)) then
               if (inl <= at%aoig(iat)%nl_sc(l-1)) then
                  !increase semicore orbitals
                  iorbsc(ispin)=iorb
               else
                  !increase valence orbitals
                  iorbv(ispin)=iorb
               end if
               icoeff=icoeff-(2*l-1)
            end do
            icoeff=icoeff+(2*l-1)

         end do
      end do

      !if (ictotpsi /= nctot) stop 'Atomic orbitals: error (nctot)'
      if (ictotpsi /= nctot) &
         & call f_err_throw('Atomic orbitals: error (nctot)', err_id=BIGDFT_INPUT_VARIABLES_ERROR)
   end do
   if (iexpo /= G%nexpo) then
      call f_err_throw('iexpo <> nexpo ' // trim(yaml_toa(iexpo)) // trim(yaml_toa(G%nexpo)), &
           & err_id=BIGDFT_INPUT_VARIABLES_ERROR)
      !write(*,*)'ERROR: iexpo <> nexpo',iexpo,G%nexpo
      !stop 
   end if

   !print *,'icoeff,ncoeff',icoeff,G%ncoeff

   if (iorbsc(1) /= norbsc) then
      call f_err_throw('Atomic orbitals error (iorbsc) ' // trim(yaml_toa(iorbsc(1))) // trim(yaml_toa(norbsc)), &
           & err_id=BIGDFT_INPUT_VARIABLES_ERROR)
      !print *,iorbsc(1),norbsc
      !stop 'Atomic orbitals: error (iorbsc)'
   end if
   if (iorbv(1) /= noncoll*norbe) call f_err_throw('Atomic orbitals: error (iorbv)',err_id=BIGDFT_INPUT_VARIABLES_ERROR)
   if (iatsc /= at%natsc) call f_err_throw('Atomic orbitals: error (iatsc)',err_id=BIGDFT_INPUT_VARIABLES_ERROR)
   !if (iorbv(1) /= noncoll*norbe) stop 'Atomic orbitals: error (iorbv)'
   !if (iatsc /= at%natsc) stop 'Atomic orbitals: error (iatsc)'

   if (nspin==2) then
      if (iorbsc(2) /= norbsc+norbe) call f_err_throw('createAtomic orbitals: error (iorbsc) nspin=2',&
                                     & err_id=BIGDFT_INPUT_VARIABLES_ERROR)
      if (iorbv(2) /= 2*norbe) call f_err_throw('createAtomic orbitals: error (iorbv) nspin=2', &
                                     & err_id=BIGDFT_INPUT_VARIABLES_ERROR)
      !if (iorbsc(2) /= norbsc+norbe) stop 'createAtomic orbitals: error (iorbsc) nspin=2'
      !if (iorbv(2) /= 2*norbe) stop 'createAtomic orbitals: error (iorbv) nspin=2'
   end if

   call f_free(xp)
   call f_free(psiat)
   call f_free(psiatn)
   call f_free(iatypex)


   if (orbse%nspinor == 4) then
      call f_free(atmoments)
   end if

   !  if (iproc ==0 .and. verbose > 1) then
   !     write(*,'(1x,a)')'done.'
   !  end if

END SUBROUTINE AtomicOrbitals


!> Calculates the kinetic energy of an atomic wavefunction expressed in Gaussians
!! the output psiatn is a normalized version of psiat
subroutine atomkin(l,ng,xp,psiat,psiatn,ek)
   use module_base
   !use module_types, only: f_err_throw, BIGDFT_RUNTIME_ERROR
   use yaml_output, only: yaml_toa
   implicit none
   integer, intent(in) :: l,ng
   real(gp), dimension(ng), intent(in) :: xp,psiat
   real(gp), intent(out) :: ek
   real(gp), dimension(ng), intent(out) :: psiatn
   !local variables
   integer :: i,j
   real(gp) :: gml,tt,xpi,xpj,d,sxp,const,hij,sij

   !        gml=.5d0*gamma_restricted(.5d0+l)
   gml = 0.0_gp

   select case(l)
   case(0) 
      gml=0.88622692545275801365_gp
   case(1)
      gml=0.44311346272637900682_gp
   case(2)
      gml=0.66467019408956851024_gp
   case(3) 
      gml=1.6616754852239212756_gp
   case default
      call f_err_throw('atomkin, l too big (' // trim(yaml_toa(l)) // ')',err_id=BIGDFT_RUNTIME_ERROR)
      !stop 'atomkin'
   end select

   ek=0.0_gp
   tt=0.0_gp
   do i=1,ng
      xpi=.5_gp/xp(i)**2
      do j=1,ng
         xpj=.5_gp/xp(j)**2
         d=xpi+xpj
         sxp=1.0_gp/d
         const=gml*sqrt(sxp)**(2*l+1)
         ! kinetic energy  matrix element hij
         hij=.5_gp*const*sxp**2* ( 3._gp*xpi*xpj +                  &
            &   real(l,gp)*(6._gp*xpi*xpj-xpi**2-xpj**2) -        &
            &   real(l**2,gp)*(xpi-xpj)**2  ) + .5_gp*real(l,gp)*(real(l,gp)+1._gp)*const
         sij=const*sxp*(real(l,gp)+.5_gp)
         ek=ek+hij*psiat(i)*psiat(j)
         tt=tt+sij*psiat(i)*psiat(j)
      enddo
   enddo

   if (abs(tt-1._gp).gt.1.e-2_gp) write(*,*) 'presumably wrong inguess data',l,tt
   ! energy expectation value
   ek=ek/tt
   !write(*,*) 'ek=',ek,tt,l,ng
   ! scale atomic wavefunction
   tt=sqrt(1._gp/tt)
   !!!        if (l.eq.0) then  ! multiply with 1/sqrt(4*pi)
   !!!        tt=tt*0.28209479177387814347_gp
   !!!        else if (l.eq.1) then  ! multiply with sqrt(3/(4*pi))
   !!!        tt=tt*0.48860251190291992159_gp
   !!!        !decide the value of the normalization to be used
   !!!        endif
   do i=1,ng
      psiatn(i)=psiat(i)*tt
   enddo

END SUBROUTINE atomkin


subroutine calc_coeff_inguess(l,m,nterm_max,nterm,lx,ly,lz,fac_arr)
   use module_base
   use yaml_output, only: yaml_toa
   implicit none
   integer, intent(in) :: l,m,nterm_max
   integer, intent(out) :: nterm
   integer, dimension(nterm_max), intent(out) :: lx,ly,lz
   real(gp), dimension(nterm_max), intent(out) :: fac_arr

   if (l.eq.1 .and. m.eq.1) then
      nterm=1
      lx(1)=0 ; ly(1)=0 ; lz(1)=0
      fac_arr(1)=0.28209479177387814347_gp

   else if (l.eq.2  .and. m.eq.1) then
      nterm=1
      lx(1)=1 ; ly(1)=0 ; lz(1)=0
      fac_arr(1)=0.48860251190291992159_gp
   else if (l.eq.2  .and. m.eq.2) then
      nterm=1
      lx(1)=0 ; ly(1)=1 ; lz(1)=0
      fac_arr(1)=0.48860251190291992159_gp
   else if (l.eq.2  .and. m.eq.3) then
      nterm=1
      lx(1)=0 ; ly(1)=0 ; lz(1)=1
      fac_arr(1)=0.48860251190291992159_gp

   else if (l.eq.3  .and. m.eq.1) then
      nterm=1
      lx(1)=0 ; ly(1)=1 ; lz(1)=1
      fac_arr(1)=1.092548430592079_gp
   else if (l.eq.3  .and. m.eq.2) then
      nterm=1
      lx(1)=1 ; ly(1)=0 ; lz(1)=1
      fac_arr(1)=1.092548430592079_gp
   else if (l.eq.3  .and. m.eq.3) then
      nterm=1
      lx(1)=1 ; ly(1)=1 ; lz(1)=0
      fac_arr(1)=1.092548430592079_gp
   else if (l.eq.3  .and. m.eq.4) then
      nterm=2
      lx(1)=2 ; ly(1)=0 ; lz(1)=0
      lx(2)=0 ; ly(2)=2 ; lz(2)=0
      fac_arr(1)=0.5462742152960396_gp
      fac_arr(2)=-0.5462742152960396_gp
   else if (l.eq.3  .and. m.eq.5) then 
      nterm=3
      lx(1)=2 ; ly(1)=0 ; lz(1)=0
      lx(2)=0 ; ly(2)=2 ; lz(2)=0
      lx(3)=0 ; ly(3)=0 ; lz(3)=2
      fac_arr(1)=-0.3153915652525201_gp
      fac_arr(2)=-0.3153915652525201_gp
      fac_arr(3)=2._gp*0.3153915652525201_gp

   else if (l.eq.4  .and. m.eq.1) then
      nterm=3
      lx(1)=3 ; ly(1)=0 ; lz(1)=0
      lx(2)=1 ; ly(2)=2 ; lz(2)=0
      lx(3)=1 ; ly(3)=0 ; lz(3)=2
      fac_arr(1)=0.4570457994644658_gp
      fac_arr(2)=0.4570457994644658_gp
      fac_arr(3)=-4._gp*0.4570457994644658_gp
   else if (l.eq.4  .and. m.eq.2) then
      nterm=3
      lx(1)=2 ; ly(1)=1 ; lz(1)=0
      lx(2)=0 ; ly(2)=3 ; lz(2)=0
      lx(3)=0 ; ly(3)=1 ; lz(3)=2
      fac_arr(1)=0.4570457994644658_gp
      fac_arr(2)=0.4570457994644658_gp
      fac_arr(3)=-4._gp*0.4570457994644658_gp
   else if (l.eq.4  .and. m.eq.3) then
      nterm=3
      lx(1)=2 ; ly(1)=0 ; lz(1)=1
      lx(2)=0 ; ly(2)=2 ; lz(2)=1
      lx(3)=0 ; ly(3)=0 ; lz(3)=3
      fac_arr(1)=3._gp*0.3731763325901154_gp
      fac_arr(2)=3._gp*0.3731763325901154_gp
      fac_arr(3)=-2._gp*0.3731763325901154_gp
   else if (l.eq.4  .and. m.eq.4) then
      nterm=2
      lx(1)=3 ; ly(1)=0 ; lz(1)=0
      lx(2)=1 ; ly(2)=2 ; lz(2)=0
      fac_arr(1)=0.5900435899266436_gp
      fac_arr(2)=-3._gp*0.5900435899266436_gp
   else if (l.eq.4  .and. m.eq.5) then
      nterm=2
      lx(1)=2 ; ly(1)=1 ; lz(1)=0
      lx(2)=0 ; ly(2)=3 ; lz(2)=0
      fac_arr(1)=-3._gp*0.5900435899266436_gp
      fac_arr(2)=0.5900435899266436_gp
   else if (l.eq.4  .and. m.eq.6) then
      nterm=2
      lx(1)=2 ; ly(1)=0 ; lz(1)=1
      lx(2)=0 ; ly(2)=2 ; lz(2)=1
      fac_arr(1)=1.445305721320277_gp
      fac_arr(2)=-1.445305721320277_gp
   else if (l.eq.4  .and. m.eq.7) then
      nterm=1
      lx(1)=1 ; ly(1)=1 ; lz(1)=1
      fac_arr(1)=2.890611442640554_gp
   else
      write(*,*) 'l,m',l,m
      call f_err_throw('input guess format error (l=' // trim(yaml_toa(l)) // ', m=' &
         & // trim(yaml_toa(m)) // ')',err_id=BIGDFT_RUNTIME_ERROR)
      !write(*,*) 'l,m',l,m
      !stop 'input guess format error'
   endif

END SUBROUTINE calc_coeff_inguess


!!$subroutine iguess_generator_modified(izatom,ielpsp,zion,psppar,npspcode,ngv,ngc,nlccpar,ng,nl,&
!!$      &   nmax_occ,noccmax,lmax,occup,expo,psiat,enlargerprb,gaenes_aux)
!!$   use module_base
!!$   implicit none
!!$   logical, intent(in) :: enlargerprb
!!$   integer, intent(in) :: ng,npspcode,nmax_occ,lmax,noccmax,ielpsp,izatom,ngv,ngc
!!$   real(gp), intent(in) :: zion
!!$   integer, dimension(lmax+1), intent(in) :: nl
!!$   real(gp), dimension(0:4,0:6), intent(in) :: psppar
!!$   real(gp), dimension(0:4,max((ngv*(ngv+1)/2)+(ngc*(ngc+1)/2),1)), intent(in) :: nlccpar
!!$   real(gp), dimension(noccmax,lmax+1), intent(in) :: occup
!!$   real(gp), dimension(ng+1), intent(out) :: expo
!!$   real(gp), dimension(ng+1,nmax_occ), intent(out) :: psiat
!!$   real(gp), dimension(nmax_occ),intent (out) :: gaenes_aux
!!$
!!$
!!$   !local variables
!!$   character(len=*), parameter :: subname='iguess_generator'
!!$   integer, parameter :: n_int=100
!!$   real(gp), parameter :: fact=4.0_gp
!!$   character(len=2) :: symbol
!!$   integer :: lpx,nsccode,mxpl,mxchg
!!$   integer :: l,i,j,iocc,i_all,i_stat
!!$   real(gp) :: alpz,alpl,amu,rprb,rij,a,a0,a0in,tt,ehomo,rcov
!!$   real(kind=8), dimension(6,4) :: neleconf
!!$   real(gp), dimension(4) :: gpot
!!$   real(gp), dimension(noccmax,lmax+1) :: aeval,chrg,res
!!$   real(gp), dimension(:), allocatable :: xp,alps
!!$   real(gp), dimension(:,:), allocatable :: vh,hsep,ofdcoef
!!$   real(gp), dimension(:,:,:), allocatable :: psi
!!$   real(gp), dimension(:,:,:,:), allocatable :: rmt
!!$
!!$   !filename = 'psppar.'//trim(atomname)
!!$   call to_zero(nmax_occ,gaenes_aux(1))
!!$
!!$   lpx=0
!!$   lpx_determination: do i=1,4
!!$      if (psppar(i,0) == 0.0_gp) then
!!$         exit lpx_determination
!!$      else
!!$         lpx=i-1
!!$      end if
!!$   end do lpx_determination
!!$
!!$   allocate(alps(lpx+1+ndebug),stat=i_stat)
!!$   call memocc(i_stat,alps,'alps',subname)
!!$   allocate(hsep(6,lpx+1+ndebug),stat=i_stat)
!!$   call memocc(i_stat,hsep,'hsep',subname)
!!$
!!$   !assignation of radii and coefficients of the local part
!!$   alpz=psppar(0,0)
!!$   alpl=psppar(0,0)
!!$   alps(1:lpx+1)=psppar(1:lpx+1,0)
!!$   gpot(1:4)=psppar(0,1:4)
!!$
!!$   !assignation of the coefficents for the nondiagonal terms
!!$   if (npspcode == 2) then !GTH case
!!$      do l=1,lpx+1
!!$         hsep(1,l)=psppar(l,1)
!!$         hsep(2,l)=0.0_gp
!!$         hsep(3,l)=psppar(l,2)
!!$         hsep(4,l)=0.0_gp
!!$         hsep(5,l)=0.0_gp
!!$         hsep(6,l)=psppar(l,3)
!!$      end do
!!$   else if (npspcode == 3) then !HGH case
!!$      allocate(ofdcoef(3,4+ndebug),stat=i_stat)
!!$      call memocc(i_stat,ofdcoef,'ofdcoef',subname)
!!$
!!$      ofdcoef(1,1)=-0.5_gp*sqrt(3._gp/5._gp) !h2
!!$      ofdcoef(2,1)=0.5_gp*sqrt(5._gp/21._gp) !h4
!!$      ofdcoef(3,1)=-0.5_gp*sqrt(100.0_gp/63._gp) !h5
!!$
!!$      ofdcoef(1,2)=-0.5_gp*sqrt(5._gp/7._gp) !h2
!!$      ofdcoef(2,2)=1._gp/6._gp*sqrt(35._gp/11._gp) !h4
!!$      ofdcoef(3,2)=-7._gp/3._gp*sqrt(1._gp/11._gp) !h5
!!$
!!$      ofdcoef(1,3)=-0.5_gp*sqrt(7._gp/9._gp) !h2
!!$      ofdcoef(2,3)=0.5_gp*sqrt(63._gp/143._gp) !h4
!!$      ofdcoef(3,3)=-9._gp*sqrt(1._gp/143._gp) !h5
!!$
!!$      ofdcoef(1,4)=0.0_gp !h2
!!$      ofdcoef(2,4)=0.0_gp !h4
!!$      ofdcoef(3,4)=0.0_gp !h5
!!$
!!$      !define the values of hsep starting from the pseudopotential file
!!$      do l=1,lpx+1
!!$         hsep(1,l)=psppar(l,1)
!!$         hsep(2,l)=psppar(l,2)*ofdcoef(1,l)
!!$         hsep(3,l)=psppar(l,2)
!!$         hsep(4,l)=psppar(l,3)*ofdcoef(2,l)
!!$         hsep(5,l)=psppar(l,3)*ofdcoef(3,l)
!!$         hsep(6,l)=psppar(l,3)
!!$      end do
!!$      i_all=-product(shape(ofdcoef))*kind(ofdcoef)
!!$      deallocate(ofdcoef,stat=i_stat)
!!$      call memocc(i_stat,i_all,'ofdcoef',subname)
!!$   else if (npspcode == 10 .or. npspcode == 12) then !HGH-K case
!!$      do l=1,lpx+1
!!$         hsep(1,l)=psppar(l,1) !h11
!!$         hsep(2,l)=psppar(l,4) !h12
!!$         hsep(3,l)=psppar(l,2) !h22
!!$         hsep(4,l)=psppar(l,5) !h13
!!$         hsep(5,l)=psppar(l,6) !h23
!!$         hsep(6,l)=psppar(l,3) !h33
!!$      end do
!!$   end if
!!$
!!$   !!Just for extracting the covalent radius and rprb
!!$   call eleconf(izatom,ielpsp,symbol,rcov,rprb,ehomo,neleconf,nsccode,mxpl,mxchg,amu)
!!$
!!$   if (enlargerprb) then
!!$      !experimental
!!$      rprb=100.0_gp
!!$   end if
!!$
!!$   !  occup(:,:)=0.0_gp
!!$   !   do l=0,lmax-1
!!$   !     iocc=0
!!$   !     do i=1,6
!!$   !        if (elecorbs(i,l+1) > 0.0_gp) then
!!$   !           iocc=iocc+1
!!$   !           !print *,'elecorbs',i,l,elecorbs(i,l+1),noccmax
!!$   !            if (iocc > noccmax) stop 'iguess_generator: noccmax too small'
!!$   !           occup(iocc,l+1)=elecorbs(i,l+1)
!!$   !        endif
!!$   !     end do
!!$   !     nl(l+1)=iocc
!!$   !  end do
!!$
!!$   !allocate arrays for the gatom routine
!!$   allocate(vh(4*(ng+1)**2,4*(ng+1)**2+ndebug),stat=i_stat)
!!$   call memocc(i_stat,vh,'vh',subname)
!!$   allocate(psi(0:ng,noccmax,lmax+ndebug),stat=i_stat)
!!$   call memocc(i_stat,psi,'psi',subname)
!!$   allocate(xp(0:ng+ndebug),stat=i_stat)
!!$   call memocc(i_stat,xp,'xp',subname)
!!$   allocate(rmt(n_int,0:ng,0:ng,lmax+ndebug),stat=i_stat)
!!$   call memocc(i_stat,rmt,'rmt',subname)
!!$
!!$   !can be switched on for debugging
!!$   !if (iproc.eq.0) write(*,'(1x,a,a7,a9,i3,i3,a9,i3,f5.2)')&
!!$   !     'Input Guess Generation for atom',trim(atomname),&
!!$   !     'Z,Zion=',izatom,ielpsp,'ng,rprb=',ng+1,rprb
!!$
!!$   rij=3._gp
!!$   ! exponents of gaussians
!!$   a0in=alpz
!!$   a0=a0in/rij
!!$   !       tt=sqrt(sqrt(2._gp))
!!$   tt=2._gp**.3_gp
!!$   do i=0,ng
!!$      a=a0*tt**i
!!$      xp(i)=.5_gp/a**2
!!$   end do
!!$
!!$   ! initial guess
!!$   do l=0,lmax-1
!!$      do iocc=1,noccmax
!!$         do i=0,ng
!!$            psi(i,iocc,l+1)=0.0_gp
!!$         end do
!!$      end do
!!$   end do
!!$
!!$   call crtvh(ng,lmax-1,xp,vh,rprb,fact,n_int,rmt)
!!$   call gatom(rcov,rprb,lmax-1,lpx,noccmax,occup,&
!!$      &   zion,alpz,gpot,alpl,hsep,alps,ngv,ngc,nlccpar,vh,xp,rmt,fact,n_int,&
!!$      &   aeval,ng,psi,res,chrg,2)
!!$
!!$   !post-treatment of the inguess data
!!$   do i=1,ng+1
!!$      expo(i)=sqrt(0.5_gp/xp(i-1))
!!$   end do
!!$
!!$   i=0
!!$   do l=1,4
!!$      do iocc=1,nl(l)
!!$         i=i+1
!!$         !occupat(i)=occup(iocc,l)
!!$         do j=1,ng+1
!!$            psiat(j,i)=psi(j-1,iocc,l)
!!$         end do
!!$         gaenes_aux(i) = aeval(iocc,l)
!!$      end do
!!$   end do
!!$
!!$   i_all=-product(shape(vh))*kind(vh)
!!$   deallocate(vh,stat=i_stat)
!!$   call memocc(i_stat,i_all,'vh',subname)
!!$   i_all=-product(shape(psi))*kind(psi)
!!$   deallocate(psi,stat=i_stat)
!!$   call memocc(i_stat,i_all,'psi',subname)
!!$   i_all=-product(shape(xp))*kind(xp)
!!$   deallocate(xp,stat=i_stat)
!!$   call memocc(i_stat,i_all,'xp',subname)
!!$   i_all=-product(shape(rmt))*kind(rmt)
!!$   deallocate(rmt,stat=i_stat)
!!$   call memocc(i_stat,i_all,'rmt',subname)
!!$   i_all=-product(shape(hsep))*kind(hsep)
!!$   deallocate(hsep,stat=i_stat)
!!$   call memocc(i_stat,i_all,'hsep',subname)
!!$   i_all=-product(shape(alps))*kind(alps)
!!$   deallocate(alps,stat=i_stat)
!!$   call memocc(i_stat,i_all,'alps',subname)
!!$
!!$END SUBROUTINE iguess_generator_modified


!>  Calculates the solution of the radial Schroedinger equation for a given
!!  pseudoptential.
subroutine gatom(rcov,rprb,lmax,lpx,noccmax,occup,&
      &   zion,alpz,gpot,alpl,hsep,alps,ngv,ngc,nlccpar,vh,xp,rmt,fact,nintp,&
      &   aeval,ng,psi,res,chrg,iorder)
   use module_base, only: gp,f_err_throw,BIGDFT_RUNTIME_ERROR,safe_exp
   use yaml_output, only: yaml_toa
   implicit none
   integer, parameter :: n_int=100
   !Arguments
   integer, intent(in) :: lmax,lpx,noccmax,ngv,ngc,nintp,ng,iorder
   real(gp), intent(in) :: rcov,rprb,zion,alpz,alpl
   real(gp), dimension(0:4,max((ngv*(ngv+1)/2)+(ngc*(ngc+1)/2),1)), intent(in) :: nlccpar
   real(gp), dimension(noccmax,lmax+1), intent(in) :: occup !<occupation numbers of the atomic orbitals, ordered by angular momentum
                                                            ! there is no need to specify the principal quantum number
                                                            ! but there cannot be more than noccmax orbitals per shell
   real(gp), dimension(noccmax,lmax+1), intent(out) :: chrg !<charge of the corresponding shell
   real(gp), dimension(0:ng,0:ng,4,0:ng,0:ng,4), intent(in) :: vh !<hartree potential matrix in the auxiliary basis set
   real(gp) :: psi(0:ng,noccmax,lmax+1),aeval(noccmax,lmax+1),&
   !   &   hh(0:ng,0:ng),
   ss(0:ng,0:ng),eval(0:ng),evec(0:ng,0:ng),&
      &   gpot(4),hsep(6,lpx+1),rmt(n_int,0:ng,0:ng,lmax+1),&
   !   &   pp1(0:ng,lpx+1),pp2(0:ng,lpx+1),pp3(0:ng,lpx+1),
      alps(lpx+1),&
      &   potgrd(n_int),&
      &   rho(0:ng,0:ng,lmax+1),rhoold(0:ng,0:ng,lmax+1),xcgrd(n_int),&
   !   &   occup(noccmax,lmax+1),
   !   chrg(noccmax,lmax+1),&
   !   &   vh(0:ng,0:ng,4,0:ng,0:ng,4),&
      &   res(noccmax,lmax+1),xp(0:ng)!,work(8*(ng+1)),ar(0:ng),ai(0:ng),beta(0:ng)
   !Local variables
   logical :: noproj
   integer :: i,l,k,j,it,iocc,ilcc,ig,lcx,info
   real(gp) :: sxp,rmix,rk,r,ttt,gml,gml1,gml2,gml3,tt,texp,sd,terf,evsum,evsumold
   real(gp) :: emuxc,dr,d,fact,const
   real(gp), dimension(0:ng,lpx+1) :: pp1,pp2,pp3 !<scalar products with the coefficients
   real(gp), dimension(0:ng,0:ng) :: hh !<hamiltonian matrix in the gaussian basis
   !Functions
   real(kind=8), external :: ddot,gamma_restricted
   real(gp), external :: spherical_gaussian_value,dlamch

   if(iorder /= 2 .and. iorder /= 4) then
       call f_err_throw('Can only use quadratic or quartic potential', err_id=BIGDFT_RUNTIME_ERROR)
       stop 'ERROR: can only use quadratic or quartic potential'
   end if

   if (nintp.ne.n_int) then
      call f_err_throw('n_int /= nintp',err_id=BIGDFT_RUNTIME_ERROR)
      !stop 'n_int/=nintp'
   end if

   do l=0,lmax
      if (occup(1,l+1) > 0._gp) lcx=l
   end do
   !write(6,*) 'lcx',lcx

   noproj=.true.
   do l=1,lpx+1
      noproj = noproj .and. (alps(l) .eq. 0._gp)
   end do


   ! projectors, just in case
   if (.not. noproj) then
      do l=0,lpx
         gml1=sqrt( gamma_restricted(real(l,gp)+1.5_gp) / (2._gp*alps(l+1)**(2*l+3)) )
         gml2=sqrt( gamma_restricted(real(l,gp)+3.5_gp) / (2._gp*alps(l+1)**(2*l+7)) )&
            &   /(real(l,gp)+2.5_gp)
         gml3=sqrt( gamma_restricted(real(l,gp)+5.5_gp) / (2._gp*alps(l+1)**(2*l+11)) )&
            &   /((real(l,gp)+3.5_gp)*(real(l,gp)+4.5_gp))
         tt=1._gp/(2._gp*alps(l+1)**2)
         do i=0,ng
            ttt=1._gp/(xp(i)+tt)
            pp1(i,l+1)=gml1*(sqrt(ttt)**(2*l+3))
            pp2(i,l+1)=gml2*ttt*(sqrt(ttt)**(2*l+3))
            pp3(i,l+1)=gml3*ttt**2*(sqrt(ttt)**(2*l+3))
         end do
      end do
   else
      pp1(:,:)=0._gp
      pp2(:,:)=0._gp
      pp3(:,:)=0._gp
   end if

   do l=0,lmax
      do j=0,ng
         do i=0,ng
            rho(i,j,l+1)=0._gp
         end do
      end do
   end do

   evsum=1.d30
   big_loop: do it=1,50
      evsumold=evsum
      evsum=0._gp

      ! coefficients of charge density
      do l=0,lmax
         do j=0,ng
            do i=0,ng
               rhoold(i,j,l+1)=rho(i,j,l+1)
               rho(i,j,l+1)=0._gp        
            end do
         end do
      end do

      do l=0,lmax
         do iocc=1,noccmax
            if (occup(iocc,l+1).gt.0._gp) then
               do j=0,ng
                  do i=0,ng
                     rho(i,j,l+1)=rho(i,j,l+1) + &
                        &   psi(i,iocc,l+1)*psi(j,iocc,l+1)*occup(iocc,l+1)
                  end do
               end do
            end if
         end do
      end do


      rmix=.5_gp
      if (it.eq.1) rmix=1._gp
      do l=0,lmax
         do j=0,ng
            do i=0,ng
               tt=rmix*rho(i,j,l+1) + (1._gp-rmix)*rhoold(i,j,l+1)
               rho(i,j,l+1)=tt
            end do
         end do
      end do

      ! XC potential on grid
      !        do k=1,n_int
      !           xcgrd(k)=0._gp
      !        end do
      !        do l=0,lmax
      !           do j=0,ng
      !              do i=0,ng
      !                 do k=1,n_int
      !                    xcgrd(k)=xcgrd(k)+rmt(k,i,j,l+1)*rho(i,j,l+1)
      !                 end do
      !              end do
      !           end do
      !        end do
      call DGEMV('N',n_int,(lcx+1)*(ng+1)**2,1._gp,&
         &   rmt,n_int,rho,1,0._gp,xcgrd,1)

      dr=fact*rprb/real(n_int,gp)
      do k=1,n_int
         r=(real(k,gp)-.5_gp)*dr
         !terms for nlcc, if present
         ilcc=0
         do ig=1,(ngv*(ngv+1))/2
            ilcc=ilcc+1
            xcgrd(k)=xcgrd(k)-&
               &   spherical_gaussian_value(r*r,nlccpar(0,ilcc),nlccpar(1,ilcc),0)
         end do
         do ig=1,(ngc*(ngc+1))/2
            ilcc=ilcc+1
            xcgrd(k)=xcgrd(k)+&
               &   spherical_gaussian_value(r*r,nlccpar(0,ilcc),nlccpar(1,ilcc),0)
         end do
         ! divide by 4 pi
         tt=xcgrd(k)*0.07957747154594768_gp
         ! multiply with r^2 to speed up calculation of matrix elements
         xcgrd(k)=emuxc(tt)*r**2
      end do

      loop_l: do l=0,lmax
         gml=.5_gp*gamma_restricted(.5_gp+real(l,gp))

         !  lower triangles only
         loop_i: do i=0,ng
            loop_j: do j=0,i
               d=xp(i)+xp(j)
               sxp=1._gp/d
               const=gml*sqrt(sxp)**(2*l+1)
               ! overlap
               ss(i,j)=const*sxp*(real(l,gp)+.5_gp)
               ! kinetic energy
               hh(i,j)=.5_gp*const*sxp**2* ( 3._gp*xp(i)*xp(j) +&
                  &   real(l,gp)*(6._gp*xp(i)*xp(j)-xp(i)**2-xp(j)**2) -&
                  &   real(l,gp)**2*(xp(i)-xp(j))**2  ) + .5_gp*real(l,gp)*(real(l,gp)+1._gp)*const
               ! potential energy from parabolic or quartic potential
               if (iorder==2) then
                   ! parabolic potential
                   hh(i,j)=hh(i,j) +&
                      &   .5_gp*const*sxp**2*(real(l,gp)+.5_gp)*(real(l,gp)+1.5_gp)/rprb**4 
               else if (iorder==4) then
                   ! quartic potential
                   hh(i,j)=hh(i,j) +&
                      &   .5_gp*const*sxp**2*(real(l,gp)+.5_gp)*(real(l,gp)+1.5_gp)/rprb**4 *sxp*(l+5/2)
               end if
               ! Hartree potential from ionic core charge
               tt=sqrt(1._gp+2._gp*alpz**2*d)
               select case(l)
               case(0)
                  hh(i,j)=hh(i,j) -zion/(2._gp*d*tt)
               case(1)
                  hh(i,j)=hh(i,j) -zion* &
                     &   (1._gp + 3._gp*alpz**2*d)/(2._gp*d**2*tt**3)
               case(2)
                  hh(i,j)=hh(i,j) -zion* &
                     &   (2._gp + 10._gp*alpz**2*d + 15._gp*alpz**4*d**2)/(2._gp*d**3*tt**5)
               case(3)
                  hh(i,j)=hh(i,j) -zion*3._gp* &
                     &   (2._gp +14._gp*alpz**2*d +35._gp*alpz**4*d**2 +35._gp*alpz**6*d**3)/&
                     &   (2._gp*d**4*tt**7)
               case default
                  call f_err_throw('l too big',err_id=BIGDFT_RUNTIME_ERROR)
                  !stop 'l too big'
               end select
               ! potential from repulsive gauss potential
               tt=alpl**2/(.5_gp+d*alpl**2)
               hh(i,j)=hh(i,j)+ gpot(1)*.5_gp*gamma_restricted(1.5_gp+real(l,gp))*&
                  &   tt**(1.5_gp+real(l,gp))&
                  &   + (gpot(2)/alpl**2)*.5_gp*gamma_restricted(2.5_gp+real(l,gp))*&
                  &   tt**(2.5_gp+real(l,gp))&
                  &   + (gpot(3)/alpl**4)*.5_gp*gamma_restricted(3.5_gp+real(l,gp))*&
                  &   tt**(3.5_gp+real(l,gp))&
                  &   + (gpot(4)/alpl**6)*.5_gp*gamma_restricted(4.5_gp+real(l,gp))*&
                  &   tt**(4.5_gp+real(l,gp))
               ! separable terms
               if (l.le.lpx) then
                  hh(i,j)=hh(i,j) + pp1(i,l+1)*hsep(1,l+1)*pp1(j,l+1)&
                     &   + pp1(i,l+1)*hsep(2,l+1)*pp2(j,l+1)&
                     &   + pp2(i,l+1)*hsep(2,l+1)*pp1(j,l+1)&
                     &   + pp2(i,l+1)*hsep(3,l+1)*pp2(j,l+1)&
                     &   + pp1(i,l+1)*hsep(4,l+1)*pp3(j,l+1)&
                     &   + pp3(i,l+1)*hsep(4,l+1)*pp1(j,l+1)&
                     &   + pp2(i,l+1)*hsep(5,l+1)*pp3(j,l+1)&
                     &   + pp3(i,l+1)*hsep(5,l+1)*pp2(j,l+1)&
                     &   + pp3(i,l+1)*hsep(6,l+1)*pp3(j,l+1)
               end if
               ! hartree potential from valence charge distribution
               !              tt=0._gp
               !              do lp=0,lcx
               !                 do jp=0,ng
               !                    do ip=0,ng
               !                       tt=tt + vh(ip,jp,lp+1,i,j,l+1)*rho(ip,jp,lp+1)
               !                    end do
               !                 end do
               !              end do
               tt=DDOT((lcx+1)*(ng+1)**2,vh(0,0,1,i,j,l+1),1,rho(0,0,1),1)
               hh(i,j)=hh(i,j) + tt
               ! potential from XC potential
               dr=fact*rprb/real(n_int,gp)
               !              tt=0._gp
               !              do k=1,n_int
               !                 tt=tt+xcgrd(k)*rmt(k,i,j,l+1)
               !              end do
               tt=DDOT(n_int,rmt(1,i,j,l+1),1,xcgrd(1),1)
               hh(i,j)=hh(i,j)+tt*dr
            end do loop_j
         end do loop_i

         !symmetrize matrices
         do i=0,ng
            do j=0,i-1
               ss(j,i)=ss(i,j)
               hh(j,i)=hh(i,j)
               !print *,'i,j,ss',i,j,ss(i,j)
               !print *,'i,j,ss',i,j,hh(i,j)
            end do
         end do
         ! ESSL
         !        call DSYGV(1,hh,ng+1,ss,ng+1,eval,evec,ng+1,ng+1,aux,2*ng+2)
         ! LAPACK
!!$         call dggev('N','V',ng+1,hh,ng+1,ss,ng+1,ar,ai,beta,&
!!$              ai,1,evec,ng+1,work,8*(ng+1),info)
!!$         if (info.ne.0) write(6,*) 'LAPACK',info
!!$         !naive division
!!$         do i=0,ng
!!$            call yaml_newline()
!!$            call yaml_map('beta,ai,ar',(/beta(i),ai(i),ar(i)/),&
!!$                 fmt='(1pe12.5)')
!!$            call yaml_map('e',ar(i)/beta(i),fmt='(1pe12.5)')
!!$         end do

         call DSYGV(1,'V','L',ng+1,hh,ng+1,ss,ng+1,eval,evec,(ng+1)**2,info)
!!$         do i=0,ng
!!$            call yaml_newline()
!!$            call yaml_map('e',eval(i),fmt='(1pe12.5)')
!!$         end do

         do iocc=0,noccmax-1
            do i=0,ng
               evec(i,iocc)=hh(i,iocc)
            end do
         end do
         ! end LAPACK

         do iocc=1,noccmax
            evsum=evsum+eval(iocc-1)
            aeval(iocc,l+1)=eval(iocc-1)
            do i=0,ng
               psi(i,iocc,l+1)=evec(i,iocc-1)
            end do
         end do
         !        write(6,*) 'eval',l
         !55      format(5(e14.7))
         !        write(6,55) eval 
         !        write(6,*) 'diff eval'
         !        write(6,55) (eval(i)-eval(i-1),i=1,ng)
         !        write(6,*) 'evec',l
         !33      format(10(e9.2))
         !        do i=0,ng
         !           write(6,33) (evec(i,iocc),iocc=0,noccmax-1)
         !        end do

      end do loop_l

      tt=abs(evsum-evsumold)
      !     write(6,*) 'evdiff',it,tt
      if (tt.lt.1.e-12_gp) then
         exit big_loop
      end if
   end do big_loop
   ! End of the big loopq

   call resid(lmax,lpx,noccmax,rprb,xp,aeval,psi,rho,ng,res,&
      &   zion,alpz,alpl,gpot,pp1,pp2,pp3,alps,hsep,fact,n_int,&
      &   potgrd,xcgrd,noproj)

   ! charge up to radius rcov
   if (lmax > 3) call f_err_throw('cannot calculate chrg', err_id=BIGDFT_RUNTIME_ERROR)
   !if (lmax > 3) stop 'cannot calculate chrg'
   do l=0,lmax
      do iocc=1,noccmax
         chrg(iocc,l+1)=0._gp
      end do
   end do

   do iocc=1,noccmax
      do j=0,ng
         do i=0,ng
            d=xp(i)+xp(j)
            sd=sqrt(d)
            call derf_ab(terf, sd*rcov) 
            texp=exp(-d*rcov**2)

            tt=0.4431134627263791_gp*terf/sd**3 - 0.5_gp*rcov*texp/d
            chrg(iocc,1)=chrg(iocc,1) + psi(i,iocc,1)*psi(j,iocc,1)*tt
            if (lmax == 0) cycle
            tt=0.6646701940895686_gp*terf/sd**5 + &
               &   (-0.75_gp*rcov*texp - 0.5_gp*d*rcov**3*texp)/d**2
            chrg(iocc,2)=chrg(iocc,2) + psi(i,iocc,2)*psi(j,iocc,2)*tt
            if (lmax == 1) cycle
            tt=1.661675485223921_gp*terf/sd**7 + &
               &   (-1.875_gp*rcov*texp-1.25_gp*d*rcov**3*texp-.5_gp*d**2*rcov**5*texp) &
               &   /d**3
            chrg(iocc,3)=chrg(iocc,3) + psi(i,iocc,3)*psi(j,iocc,3)*tt
            if (lmax == 2) cycle
            tt=5.815864198283725_gp*terf/sd**9 + &
               &   (-6.5625_gp*rcov*texp - 4.375_gp*d*rcov**3*texp - &
               &   1.75_gp*d**2*rcov**5*texp - .5_gp*d**3*rcov**7*texp)/d**4
            chrg(iocc,4)=chrg(iocc,4) + psi(i,iocc,4)*psi(j,iocc,4)*tt
         end do
      end do
   end do

   ! writing lines suppressed
   !!!        write(66,*)  lmax+1
   !!!        write(66,*) ' #LINETYPE{1324}' 
   !!!        write(66,*) ' $' 
   !!!  do l=0,lmax
   !!!           write(66,*) ' 161'
   !!!     r=0._gp
   !!!     do
   !!!        tt= wave(ng,l,xp,psi(0,1,l+1),r)
   !!!              write(66,*) r,tt
   !!!        r=r+.025_gp
   !!!        if(r > 4.00001_gp) exit
   !!!     end do
   !!!  end do
   ! writing lines suppressed
   !!!        write(67,*) min(lmax+1,3)
   !!!        write(67,*) ' #LINETYPE{132}'
   !!!        write(67,*) ' #TITLE{FOURIER}' 
   !!!        write(67,*) ' $'
   dr=6.28_gp/rprb/200._gp
   !!!        write(67,*) ' 200'
   rk=0._gp
   loop_rk1: do 
      tt=0._gp
      do i=0,ng
         texp=safe_exp(-.25_gp*rk**2/xp(i))
         !        texp=exp(-.5_gp*energy/xp(i))
         sd=sqrt(xp(i))
         tt=tt+psi(i,1,1)*0.4431134627263791_gp*texp/sd**3
      end do
      !!!           write(67,*) rk,tt
      rk=rk+dr
      if(rk > 6.28_gp/rprb-.5_gp*dr) exit loop_rk1
   end do loop_rk1
   if (lmax.ge.1) then
      !!!           write(67,*) ' 200'
      rk=0._gp
      loop_rk2: do 
         tt=0._gp
         do i=0,ng
            texp=safe_exp(-.25_gp*rk**2/xp(i))
            sd=sqrt(xp(i))
            tt=tt+psi(i,1,2)*0.2215567313631895_gp*rk*texp/sd**5
         end do
         !!!              write(67,*) rk,tt
         rk=rk+dr
         if (rk > 6.28_gp/rprb-.5_gp*dr) exit loop_rk2
      end do loop_rk2
   end if
   if (lmax.ge.2) then
      !!!           write(67,*) ' 200'
      rk=0._gp
      do 
         tt=0._gp
         do i=0,ng
            texp=safe_exp(-.25_gp*rk**2/xp(i))
            sd=sqrt(xp(i))
            tt=tt+psi(i,1,3)*0.1107783656815948_gp*rk**2*texp/sd**7
         end do
         !!!              write(67,*) rk,tt
         rk=rk+dr
         if (rk > 6.28_gp/rprb-.5_gp*dr) exit
      end do
   end if

END SUBROUTINE gatom


subroutine resid(lmax,lpx,noccmax,rprb,xp,aeval,psi,rho,&
      &   ng,res,zion,alpz,alpl,gpot,pp1,pp2,pp3,alps,hsep,fact,n_int,&
      &   potgrd,xcgrd,noproj)
   use module_base, only: gp,safe_exp
   implicit real(gp) (a-h,o-z)
   logical :: noproj
   dimension psi(0:ng,noccmax,lmax+1),rho(0:ng,0:ng,lmax+1),&
      &   gpot(4),pp1(0:ng,lmax+1),pp2(0:ng,lmax+1),pp3(0:ng,lmax+1),&
      &   alps(lmax+1),hsep(6,lmax+1),res(noccmax,lmax+1),xp(0:ng),&
      &   xcgrd(n_int),aeval(noccmax,lmax+1),potgrd(n_int)
   real(gp) :: derf_val

   ! potential on grid 
   dr=fact*rprb/real(n_int,gp)
   do k=1,n_int
      r=(real(k,gp)-.5_gp)*dr
      call derf_ab(derf_val, r/(sqrt(2._gp)*alpz))
      potgrd(k)= .5_gp*(r/rprb**2)**2 - &
         &   zion*derf_val/r &
         &   + safe_exp(-.5_gp*(r/alpl)**2)*&
         &   ( gpot(1) + gpot(2)*(r/alpl)**2 + gpot(3)*(r/alpl)**4 + gpot(4)*(r/alpl)**6  )&
         &   + xcgrd(k)/r**2
      do j=0,ng
         do i=0,ng
            spi=1.772453850905516_gp
            d=xp(i)+xp(j)
            sd=sqrt(d)
            tx=safe_exp(-d*r**2)
            call derf_ab(tt, sd*r)
            tt=spi*tt
            u_gp=tt/(4._gp*sd**3*r)
            potgrd(k)=potgrd(k)+u_gp*rho(i,j,1)
            ud1=-tx/(4._gp*d**2) + 3._gp*tt/(8._gp*sd**5*r)
            if (lmax.ge.1) potgrd(k)=potgrd(k)+ud1*rho(i,j,2)
            ud2=-tx*(7._gp + 2._gp*d*r**2)/(8._gp*d**3) +&
               &   15._gp*tt/(16._gp*sd**7*r)
            if (lmax.ge.2) potgrd(k)=potgrd(k)+ud2*rho(i,j,3)
            ud3=-tx*(57._gp+22._gp*d*r**2+4._gp*d**2*r**4)/(16._gp*d**4) + &
               &   105._gp*tt/(32._gp*sd**9*r)
            if (lmax.ge.3) potgrd(k)=potgrd(k)+ud3*rho(i,j,4)
         end do
      end do
   end do

   loop_ll: do ll=0,lmax
      if (ll.le.lpx .and. .not. noproj) then
         rnrm1=1._gp/sqrt(.5_gp*gamma_restricted(real(ll,gp)+1.5_gp)*alps(ll+1)**(2*ll+3))
         rnrm2=1._gp/sqrt(.5_gp*gamma_restricted(real(ll,gp)+3.5_gp)*alps(ll+1)**(2*ll+7))
         rnrm3=1._gp/sqrt(.5_gp*gamma_restricted(real(ll,gp)+5.5_gp)*alps(ll+1)**(2*ll+11))
      end if
      loop_iocc: do iocc=1,noccmax
         ! separable part
         if (ll.le.lpx) then
            scpr1=DDOT(ng+1,psi(0,iocc,ll+1),1,pp1(0,ll+1),1)
            scpr2=DDOT(ng+1,psi(0,iocc,ll+1),1,pp2(0,ll+1),1)
            scpr3=DDOT(ng+1,psi(0,iocc,ll+1),1,pp3(0,ll+1),1)
         end if
         res(iocc,ll+1)=0._gp
         loop_j: do j=1,n_int
            ! wavefunction on grid
            r=(real(j,gp)-.5_gp)*dr
            psigrd = wave(ng,ll,xp,psi(0,iocc,ll+1),r)
            ! kinetic energy        
            rkin=0._gp
            do i=0,ng
               rkin=rkin + psi(i,iocc,ll+1) *  (&
                  &   xp(i)*(3._gp+2._gp*real(ll,gp)-2._gp*xp(i)*r**2)*&
                  safe_exp(-xp(i)*r**2) )
            end do
            rkin=rkin*r**ll
            ! separable part
            if (ll.le.lpx .and. .not. noproj) then
               sep =& 
               (scpr1*hsep(1,ll+1) + scpr2*hsep(2,ll+1) + scpr3*hsep(4,ll+1))&
                  &   *rnrm1*r**ll*safe_exp(-.5_gp*(r/alps(ll+1))**2)   +&
                  &   (scpr1*hsep(2,ll+1) + scpr2*hsep(3,ll+1) + scpr3*hsep(5,ll+1))&
                  &   *rnrm2*r**(ll+2)*safe_exp(-.5_gp*(r/alps(ll+1))**2)   +&
                  &   (scpr1*hsep(4,ll+1) + scpr2*hsep(5,ll+1) + scpr3*hsep(6,ll+1))&
                  &   *rnrm3*r**(ll+4)*safe_exp(-.5_gp*(r/alps(ll+1))**2)
            else
               sep=0._gp
            end if
            ! residue
            tt=rkin+sep+(potgrd(j)-aeval(iocc,ll+1))*psigrd
            !384        format(6(e12.5))
            !12        format(i2,i2,e9.2,3(e12.5),e10.3)
            res(iocc,ll+1)=res(iocc,ll+1) + tt**2*dr
         end do loop_j
      end do loop_iocc
   end do loop_ll
   !  do l=0,lmax
   !     do iocc=1,noccmax
   !        write(6,*) 'res',l,iocc,res(iocc,l+1)
   !     end do
   !  end do

END SUBROUTINE resid


subroutine crtvh(ng,lmax,xp,vh,rprb,fact,n_int,rmt)
   use module_base, only: gp,f_err_throw,BIGDFT_RUNTIME_ERROR,safe_exp
   use yaml_output, only: yaml_toa
   implicit none
   !implicit real(gp) (a-h,o-z)
   !Arguments
   integer, intent(in) :: ng,lmax,n_int
   real(gp), intent(in) ::rprb,fact
   real(gp), dimension(0:ng), intent(in) :: xp
   real(gp), dimension(0:ng,0:ng,0:3,0:ng,0:ng,0:3), intent(out) :: vh
   real(gp), dimension(n_int,0:ng,0:ng,lmax+1), intent(out) :: rmt
   !Local variables
   integer :: i,j,ip,jp,k,l
   real(gp) :: c,d,scpd,dr,r

   if (lmax > 3) call f_err_throw('crtvh', err_id=BIGDFT_RUNTIME_ERROR)
   !if (lmax > 3) stop 'crtvh'

   dr=fact*rprb/real(n_int,gp)
   do l=0,lmax
      do k=1,n_int
         r=(real(k,gp)-.5_gp)*dr
         do j=0,ng
            do i=0,ng
               !rmt(k,i,j,l+1)=(r**2)**l*exp(-(xp(i)+xp(j))*r**2)
               rmt(k,i,j,l+1)=(r**2)**l*safe_exp(-(xp(i)+xp(j))*r**2)
            end do
         end do
      end do
   end do

   loop_j: do j=0,ng
      loop_i: do i=0,ng
         c=xp(i)+xp(j)
         loop_jp: do jp=0,ng
            loop_ip: do ip=0,ng
               d=xp(ip)+xp(jp)
               scpd=sqrt(c+d)
               vh(ip,jp,0,i,j,0)=0.2215567313631895_gp/(c*d*scpd)
               vh(ip,jp,1,i,j,0)=&
                  &   .1107783656815948_gp*(2._gp*c+3._gp*d)/(c*d**2*scpd**3)
               vh(ip,jp,2,i,j,0)=.05538918284079739_gp*&
                  &   (8._gp*c**2+20._gp*c*d+15._gp*d**2)/(c*d**3*scpd**5)
               vh(ip,jp,3,i,j,0)=.0830837742611961_gp*&
                  &   (16._gp*c**3+56._gp*c**2*d+70._gp*c*d**2+35._gp*d**3)/&
                  &   (c*d**4*scpd**7)

               vh(ip,jp,0,i,j,1)=&
                  &   .1107783656815948_gp*(3._gp*c+2._gp*d)/(c**2*d*scpd**3)
               vh(ip,jp,1,i,j,1)=&
                  &   .05538918284079739_gp*(6._gp*c**2+15._gp*c*d+6._gp*d**2)/&
                  &   (c**2*d**2*scpd**5)
               vh(ip,jp,2,i,j,1)=.02769459142039869_gp*&
                  &   (24._gp*c**3+84._gp*c**2*d+105._gp*c*d**2+30._gp*d**3)/&
                  &   (c**2*d**3*scpd**7)
               vh(ip,jp,3,i,j,1)=0.04154188713059803_gp*&
                  &   (48._gp*c**4+216._gp*c**3*d+378._gp*c**2*d**2+&
                  &   315._gp*c*d**3+70._gp*d**4)/(c**2*d**4*scpd**9)

               vh(ip,jp,0,i,j,2)=&
                  &   .05538918284079739_gp*(15._gp*c**2+20._gp*c*d+8._gp*d**2)/&
                  &   (c**3*d*scpd**5)
               vh(ip,jp,1,i,j,2)=.02769459142039869_gp*&
                  &   (30._gp*c**3+105._gp*c**2*d+84._gp*c*d**2+24._gp*d**3)/&
                  &   (c**3*d**2*scpd**7)
               vh(ip,jp,2,i,j,2)=&
                  &   .2077094356529901_gp*(8._gp*c**4+36._gp*c**3*d+63._gp*c**2*d**2+&
                  &   36._gp*c*d**3+8._gp*d**4)/(c**3*d**3*scpd**9)
               vh(ip,jp,3,i,j,2)=&
                  &   .1038547178264951_gp*(48._gp*c**5+264._gp*c**4*d+594._gp*c**3*d**2+&
                  &   693._gp*c**2*d**3+308._gp*c*d**4+56._gp*d**5)/&
                  &   (c**3*d**4*scpd**11)

               vh(ip,jp,0,i,j,3)=.0830837742611961_gp*&
                  &   (35._gp*c**3+70._gp*c**2*d+56._gp*c*d**2+16._gp*d**3)/&
                  &   (c**4*d*scpd**7)
               vh(ip,jp,1,i,j,3)=&
                  &   .04154188713059803_gp*(70._gp*c**4+315._gp*c**3*d+378._gp*c**2*d**2+&
                  &   216._gp*c*d**3+48._gp*d**4)/(c**4*d**2*scpd**9)
               vh(ip,jp,2,i,j,3)=&
                  &   .1038547178264951_gp*(56._gp*c**5+308._gp*c**4*d+693._gp*c**3*d**2+&
                  &   594._gp*c**2*d**3+264._gp*c*d**4+48._gp*d**5)/&
                  &   (c**4*d**3*scpd**11)
               vh(ip,jp,3,i,j,3)=&
                  &   1.090474537178198_gp*(16._gp*c**6+104._gp*c**5*d+286._gp*c**4*d**2+&
                  &   429._gp*c**3*d**3+286._gp*c**2*d**4+104._gp*c*d**5+16._gp*d**6)/&
                  &   (c**4*d**4*scpd**13)
            end do loop_ip
         end do loop_jp
      end do loop_i
   end do loop_j

END SUBROUTINE crtvh


function wave(ng,ll,xp,psi,r)
   use module_base, only: gp,safe_exp
   implicit none
   !Arguments
   integer, intent(in) :: ll,ng
   real(gp) :: r,wave
   real(gp) :: psi(0:ng),xp(0:ng)
   !Local variables
   integer :: i

   wave=0._gp
   do i=0,ng
      wave=wave + psi(i)*safe_exp(-xp(i)*r**2)
   end do
   if(ll>0)then
      wave=wave*r**ll
   endif
END FUNCTION wave



!>
!!
!!
function emuxc(rho)
   use module_base, only: gp
   implicit none
   real(gp), intent(in) :: rho
   real(gp) :: emuxc
   real(gp), parameter :: &
      &   a0p=.4581652932831429_gp,&
      &   a1p=2.217058676663745_gp,&
      &   a2p=0.7405551735357053_gp,&
      &   a3p=0.01968227878617998_gp
   real(gp), parameter :: &
      &   b1p=1.0_gp,&
      &   b2p=4.504130959426697_gp,&
      &   b3p=1.110667363742916_gp,&
      &   b4p=0.02359291751427506_gp
   real(gp), parameter :: rsfac=.6203504908994000_gp,ot=1._gp/3._gp
   real(gp), parameter :: &
      &   c1=4._gp*a0p*b1p/3.0_gp,  &
      &   c2=5.0_gp*a0p*b2p/3.0_gp+a1p*b1p,&
      &   c3=2.0_gp*a0p*b3p+4.0_gp*a1p*b2p/3.0_gp+2.0_gp*a2p*b1p/3.0_gp,&
      &   c4=7.0_gp*a0p*b4p/3.0_gp+5.0_gp*a1p*b3p/3.0_gp+a2p*b2p+a3p*b1p/3.0_gp,&
      &   c5=2.0_gp*a1p*b4p+4.0_gp*a2p*b3p/3.0_gp+2.0_gp*a3p*b2p/3.0_gp,&
      &   c6=5.0_gp*a2p*b4p/3.0_gp+a3p*b3p,c7=4.0_gp*a3p*b4p/3.0_gp
   real(gp) :: bot,rs,top

   if(rho.lt.1.e-24_gp) then
      emuxc=0._gp
   else
      if(rho.lt.0._gp) write(6,*) ' rho less than zero',rho
      rs=rsfac*rho**(-ot)
      top=-rs*(c1+rs*(c2+rs*(c3+rs*(c4+rs*(c5+rs*(c6+rs*c7))))))
      bot=rs*(b1p+rs*(b2p+rs*(b3p+rs*b4p)))
      emuxc=top/(bot*bot)
   end if
END FUNCTION emuxc


!> Restricted version of the Gamma function
function gamma_restricted(x)
   use module_base, only: gp,f_err_throw,BIGDFT_RUNTIME_ERROR
   use yaml_output, only: yaml_toa
   implicit none
   !Arguments
   real(gp), intent(in) :: x
   real(gp) :: gamma_restricted
   !Local variables
   integer :: ii,i

   if (x <= 0._gp) call f_err_throw('wrong argument for gamma_restricted', err_id=BIGDFT_RUNTIME_ERROR)
   !if (x <= 0._gp) stop 'wrong argument for gamma_restricted'
   if (mod(x,1._gp) == 0._gp) then
      ii=int(x)
      gamma_restricted=1.0_gp
      do i=2,ii
         gamma_restricted=gamma_restricted*real(i-1,gp)
      end do
   else if (mod(x,.5_gp) == 0._gp) then
      ii=int(x-.5_gp)
      !     gamma_restricted=sqrt(3.14159265358979_gp)
      gamma_restricted=1.772453850905516027_gp
      do i=1,ii
         gamma_restricted=gamma_restricted*(real(i,gp)-.5_gp)
      end do
   else
      call f_err_throw('wrong argument for gamma_restricted', err_id=BIGDFT_RUNTIME_ERROR)
      !stop 'wrong argument for gamma_restricted'
   end if
END FUNCTION gamma_restricted


!  call psitospi0(iproc,nproc,norbe,norbep,&
!       wfd%nvctr_c,wfd%nvctr_f,nspin,spinsgne,psi)
subroutine psitospi0(iproc,nproc,norbe,norbep,&
      &   nvctr_c,nvctr_f,nspin,spinsgne,psi)
   use module_base
  use yaml_output
   implicit none
   !Arguments
   integer, intent(in) :: norbe,norbep,iproc,nproc
   integer, intent(in) :: nvctr_c,nvctr_f
   integer, intent(in) :: nspin
   integer, dimension(norbe*nspin), intent(in) :: spinsgne
   real(kind=8), dimension(nvctr_c+7*nvctr_f,norbep*nspin), intent(inout) :: psi
   !Local variables
   character(len=*), parameter :: subname='psitospi0'
   logical :: myorbital
   integer :: nvctr
   integer :: iorb,jorb,i
   real(kind=8) :: facu,facd
   real(kind=8), dimension(:,:), allocatable :: psi_o
   !n(c) integer, dimension(2) :: iorbsc,iorbv

   !initialise the orbital counters
   !n(c) iorbsc(1)=0
   !n(c) iorbv(1)=norbsc
   !used in case of spin-polarisation, ignored otherwise
   !n(c) iorbsc(2)=norbe
   !n(c) iorbv(2)=norbsc+norbe

   !if (iproc ==0) write(*,'(1x,a)',advance='no')'Transforming AIO to spinors...'
   if (iproc ==0) call yaml_map('Transforming AIO to spinors',.true.)

   nvctr=nvctr_c+7*nvctr_f
   psi_o = f_malloc((/ nvctr, norbep /),id='psi_o')

   do iorb=1,norbep
      do i=1,nvctr
         psi_o(i,iorb)=psi(i,iorb)
      end do
   end do

   call to_zero(nvctr*nspin*norbep,psi)

   do iorb=1,norbe
      jorb=iorb-iproc*norbep
      if (myorbital(iorb,nspin*norbe,iproc,nproc)) then
         if(spinsgne(jorb)>0.0d0) then
            facu=1.0d0
            facd=0.0d0
         else
            facu=0.0d0
            facd=1.0d0
         end if
         do i=1,nvctr
            psi(i,iorb*4-3) = facu*psi_o(i,iorb)
            psi(i,iorb*4-2) = .0d0*psi_o(i,iorb)
            psi(i,iorb*4-1) = facd*psi_o(i,iorb)
            psi(i,iorb*4)   = .0d0*psi_o(i,iorb)
         end do
      end if
   end do
   call f_free(psi_o)

   !if (iproc ==0) write(*,'(1x,a)')'done.'

END SUBROUTINE psitospi0

!> Control whether the occupation number can be rounded by a shell-dependent fraction 
!! denominator
subroutine write_fraction_string(l,occ,string,nstring)
   use module_base
   implicit none
   integer, intent(in) :: l
   real(gp), intent(in) :: occ
   character(len=10), intent(out) :: string
   integer, intent(out) :: nstring
   !local variables
   real(gp), parameter :: occtol=1e-10_gp
   integer :: num,den
   real(gp) :: check

   den=(2*l-1)
   check=occ*real(den,gp)
   !find the numerator
   num=nint(check)
   if (abs(check - real(num,gp)) < occtol .and. &
      &   (occ /= 0.0_gp .and. occ /= 1.0_gp .and. occ /= 2.0_gp)) then
   !the length of nstring depends of the l value
   if (l >3) then
      nstring=6
      write(string,'(1x,i2,a,i2)') num,'/',den
   else
      nstring=4
      write(string,'(1x,i1,a,i1)') num,'/',den
   end if
else
   nstring=5
   write(string,'(1x,f4.2)') occ
end if

END SUBROUTINE write_fraction_string


!!$!> Print the electronic configuration, with the semicore orbitals
!!$subroutine aocc_to_dict(dict, nzatom, nelpsp, nspin, nspinor, aocc, nsccode)
!!$   use module_defs, only: gp
!!$   use dictionaries
!!$   implicit none
!!$   integer, parameter :: nelecmax=32
!!$   type(dictionary), pointer :: dict
!!$   integer, intent(in) :: nzatom,nelpsp,nspinor,nspin,nsccode
!!$   real(gp), dimension(nelecmax), intent(in) :: aocc
!!$   !local variables
!!$   character(len=10) :: tmp
!!$   character(len=500) :: string
!!$   integer :: i,m,iocc,icoll,inl,noncoll,l,ispin,is,nl,niasc,lsc,nlsc,ntmp,iss
!!$   logical, dimension(4,2) :: scorb
!!$
!!$   !if non-collinear it is like nspin=1 but with the double of orbitals
!!$   if (nspinor == 4) then
!!$      noncoll=2
!!$   else
!!$      noncoll=1
!!$   end if
!!$   scorb=.false.
!!$   if (nsccode/=0) then !the atom has some semicore orbitals
!!$      niasc=nsccode
!!$      do lsc=4,1,-1
!!$         nlsc=niasc/4**(lsc-1)
!!$         do i=1,nlsc
!!$            scorb(lsc,i)=.true.
!!$         end do
!!$         niasc=niasc-nlsc*4**(lsc-1)
!!$      end do
!!$   end if
!!$
!!$   call dict_init(dict)
!!$
!!$   !initalise string
!!$   string=repeat(' ',len(string))
!!$
!!$   is=1
!!$   do i=1,noccmax
!!$      iocc=0
!!$      do l=1,lmax
!!$         iocc=iocc+1
!!$         nl=nint(aocc(iocc))
!!$         do inl=1,nl
!!$            !write to the string the angular momentum
!!$            if (inl == i) then
!!$               iss=is
!!$               if (scorb(l,inl)) then
!!$                  string(is:is)='('
!!$                  is=is+1
!!$               end if
!!$               select case(l)
!!$               case(1)
!!$                  string(is:is)='s'
!!$               case(2)
!!$                  string(is:is)='p'
!!$               case(3)
!!$                  string(is:is)='d'
!!$               case(4)
!!$                  string(is:is)='f'
!!$               case default
!!$                  stop 'l not admitted'
!!$               end select
!!$               is=is+1
!!$               if (scorb(l,inl)) then
!!$                  string(is:is)=')'
!!$                  is=is+1
!!$               end if
!!$               call yaml_sequence_open(string(iss:is))
!!$            end if
!!$            do ispin=1,nspin
!!$               do m=1,2*l-1
!!$                  do icoll=1,noncoll !non-trivial only for nspinor=4
!!$                     iocc=iocc+1
!!$                     !write to the string the value of the occupation numbers
!!$                     if (inl == i) then
!!$                        call write_fraction_string(l,aocc(iocc),tmp,ntmp)
!!$                        string(is:is+ntmp-1)=tmp(1:ntmp)
!!$                        call yaml_sequence(tmp(1:ntmp))
!!$                        is=is+ntmp
!!$                     end if
!!$                  end do
!!$               end do
!!$            end do
!!$            if (inl == i) then
!!$               string(is:is+2)=' , '
!!$               is=is+3
!!$               call yaml_sequence_close()
!!$            end if
!!$         end do
!!$      end do
!!$   end do
!!$
!!$   !write(*,'(2x,a,1x,a,1x,a)',advance='no')' Elec. Configuration:',trim(string),'...'
!!$
!!$ END SUBROUTINE aocc_to_dict
