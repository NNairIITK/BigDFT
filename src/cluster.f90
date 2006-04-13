

        implicit real*8 (a-h,o-z)
! atomic coordinates, forces
        real*8, allocatable, dimension(:,:) :: rxyz, fxyz
        logical parallel
        character*20 tatonam
! atomic types
        integer, allocatable, dimension(:) :: iatype
        character*20 :: atomnames(100), units
        include 'mpif.h'

        parallel=.true.

! For parallel MPI execution set parallel=.true., for serial parallel=.false.
! Start MPI in parallel version
        if (parallel) then
        call MPI_INIT(ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
        write(6,*) 'mpi started',iproc,nproc
        else
        nproc=1
        iproc=0
        endif

! read atomic positions
        open(unit=9,file='posinp',status='old')
        read(9,*) nat,units
        if (iproc.eq.0) write(6,*) 'nat=',nat
        allocate(rxyz(3,nat),iatype(nat),fxyz(3,nat))
        ntypes=0
        do iat=1,nat
        read(9,*) rxyz(1,iat),rxyz(2,iat),rxyz(3,iat),tatonam
         do ityp=1,ntypes
           if (tatonam.eq.atomnames(ityp)) then
              iatype(iat)=ityp
              goto 200
           endif
         enddo
         ntypes=ntypes+1
         if (ntypes.gt.100) stop 'more than 100 atomnames not permitted'
         atomnames(ityp)=tatonam
         iatype(iat)=ntypes
200        continue
        if (units.eq.'angstroem') then
! if Angstroem convert to Bohr
        do i=1,3 ;  rxyz(i,iat)=rxyz(i,iat)/.529177d0  ; enddo
        else if  (units.eq.'atomic' .or. units.eq.'bohr') then
        else
        write(*,*) 'length units in input file unrecognized'
        write(*,*) 'recognized units are angstroem or atomic = bohr'
        stop 
        endif
        enddo
        close(9)
        do ityp=1,ntypes
        if (iproc.eq.0) write(*,*) 'atoms of type ',ityp,' are ',atomnames(ityp)
        enddo

        call cluster(parallel,nproc,iproc,nat,ntypes,iatype,atomnames,rxyz,energy,fxyz)

        deallocate(rxyz,iatype,fxyz)

        if (parallel) call MPI_FINALIZE(ierr)

	end


        subroutine cluster(parallel,nproc,iproc,nat,ntypes,iatype,atomnames,rxyz,energy,fxyz)
! does an electronic structure calculation. Output is the total energy 
! and the forces (not yet impleneted)
        implicit real*8 (a-h,o-z)
        integer count1,count2,count_rate,count_max
        character*30 label ; character*27 filename ; character*20 atomnames ; character*4 f4
        logical logrid,calc_inp_wf,parallel,loregion
! work array for ALLREDUCE
        dimension wrkallred(3,2) 
! work array for crtproj
        dimension lx(3),ly(3),lz(3)
! atomic coordinates
        dimension rxyz(3,nat),fxyz(3,nat),iatype(nat),atomnames(100)
! active grid points, segments of real space grid
        allocatable :: logrid(:,:,:)
! occupation numbers
        allocatable :: occup(:)

! wavefunction segments
        allocatable :: keyv(:), nseg(:)
! wavefunction segments on real space grid
        allocatable :: keyg(:,:),nvctr(:)
! wavefunction 
        allocatable :: psi(:)
! wavefunction gradients
        allocatable :: hpsi(:)

! Charge density/potential,ionic potential, pkernel
        allocatable :: rhopot(:),pot_ion(:),pkernel(:)

! projector segments
        allocatable :: keyv_p(:), nseg_p(:)
! projector segments on real space grid
        allocatable :: keyg_p(:,:),nvctr_p(:)
! projectors 
        allocatable :: proj(:)
! Parameters for the boxes containing the projectors
        allocatable :: nboxp_c(:,:,:),nboxp_f(:,:,:)

! pseudopotential parameters
        allocatable :: psppar(:,:,:),nelpsp(:),radii_cf(:,:),rad_cov(:)
! Parameters determining the localization region of the Wannier functions
        allocatable :: wan_par(:,:),loregion(:,:)
! Parameters for the corresponding boxes
        allocatable :: nbox_c(:,:,:),nbox_f(:,:,:)
        include 'mpif.h'

	nloewe=1
        if (iproc.eq.0) write(*,*) 'loewdin orthogonalization performed ',nloewe,' times'

        open(unit=78,file='timings')
        open(unit=79,file='malloc')
        call cpu_time(tcpu1)

        open(unit=1,file='input.dat',status='old')
	read(1,*) hgrid
	read(1,*) crmult
	read(1,*) frmult
	read(1,*) cpmult
	read(1,*) fpmult
                  if (fpmult.gt.frmult) write(*,*) 'NONSENSE: fpmult > frmult'
        read(1,*) radlocmult
	read(1,*) gnrm_cv
	read(1,*) itermax
	read(1,*) calc_inp_wf
        close(1)

        if (iproc.eq.0) then 
          write(*,*) 'hgrid=',hgrid
          write(*,*) 'crmult=',crmult
          write(*,*) 'frmult=',frmult
          write(*,*) 'cpmult=',cpmult
          write(*,*) 'fpmult=',fpmult
          write(*,*) 'gnrm_cv=',gnrm_cv
          write(*,*) 'itermax=',itermax
          write(*,*) 'calc_inp_wf=',calc_inp_wf
        endif


! grid spacing (same in x,y and z direction)
        hgridh=.5d0*hgrid

! store PSP parameters
        allocate(psppar(0:2,0:4,ntypes),nelpsp(ntypes),radii_cf(ntypes,2),rad_cov(ntypes))
      do ityp=1,ntypes
        filename = 'psppar.'//atomnames(ityp)
!        if (iproc.eq.0) write(*,*) 'opening PSP file ',filename
        open(unit=11,file=filename,status='old')
	read(11,'(a30)') label
	read(11,*) radii_cf(ityp,1),radii_cf(ityp,2),rad_cov(ityp)
	read(11,*) nelpsp(ityp)
        if (iproc.eq.0) write(*,'(a,1x,a,a,i3,a,a)') 'atom type ',atomnames(ityp), & 
                        ' is described by a ',nelpsp(ityp),' electron',label
        do i=0,2
	read(11,*) (psppar(i,j,ityp),j=0,4)
	enddo
        close(11)
      enddo


! Number of orbitals and their occupation number
         norb_vir=0

         nelec=0
         do iat=1,nat
         ityp=iatype(iat)
         nelec=nelec+nelpsp(ityp)
         enddo
         if (iproc.eq.0) write(*,*) 'number of electrons',nelec
         if (mod(nelec,2).ne.0) write(*,*) 'WARNING: odd number of electrons, no closed shell system'
         norb=(nelec+1)/2+norb_vir

         allocate(occup(norb))

         nt=0
         do iorb=1,norb
         it=min(2,nelec-nt)
         occup(iorb)=it
         nt=nt+it
         enddo

         if (iproc.eq.0) then 
            write(*,*) 'number of orbitals',norb
         do iorb=1,norb
           write(*,*) 'occup(',iorb,')=',occup(iorb)
         enddo
         endif

! determine size alat of overall simulation cell
        call system_size(nat,rxyz,radii_cf(1,1),crmult,iatype,ntypes, &
                   cxmin,cxmax,cymin,cymax,czmin,czmax)
        alat1=(cxmax-cxmin)
        alat2=(cymax-cymin)
        alat3=(czmax-czmin)

! shift atomic positions such that molecule is inside cell
        if (iproc.eq.0) write(*,'(a,3(1x,e12.5))') 'atomic positions shifted',-cxmin,-cymin,-czmin
!        if (iproc.eq.0) write(*,'(a,3(1x,e24.17))') 'atomic positions shifted',-cxmin,-cymin,-czmin
	do iat=1,nat
        rxyz(1,iat)=rxyz(1,iat)-cxmin
        rxyz(2,iat)=rxyz(2,iat)-cymin
        rxyz(3,iat)=rxyz(3,iat)-czmin
        enddo

!    Find grid sizes n1,n2,n3
         n1=int(alat1/hgrid)
         n2=int(alat2/hgrid)
         n3=int(alat3/hgrid)
        alat1=n1*hgrid ; alat2=n2*hgrid ; alat3=n3*hgrid
         if (iproc.eq.0) then 
           write(*,*) 'n1,n2,n3',n1,n2,n3
           write(*,*) 'total number of grid points',(n1+1)*(n2+1)*(n3+1)
           write(*,'(a,3(1x,e12.5))') 'simulation cell',alat1,alat2,alat3
!           write(*,'(a,3(1x,e24.17))') 'simulation cell',alat1,alat2,alat3
         endif

! determine localization regions in terms of atomic sphere centers
        allocate(wan_par(0:3,norb))
	call wannier_par(iproc,nat,norb,rxyz,ntypes,iatype,rad_cov,wan_par)
        allocate(loregion(nat,norb))
        call localizationregion(iproc,nat,norb,rxyz,wan_par,radlocmult,loregion)

! determine localization region for all Wanier orbitals, but do not yet fill the descriptor arrays
    allocate(logrid(0:n1,0:n2,0:n3))
    allocate(nbox_c(2,3,norb),nbox_f(2,3,norb))
    allocate(nseg(0:2*norb))
    allocate(nvctr(0:2*norb))
    nseg(0)=0 
    nvctr(0)=0 

    do iorb=1,norb
    if (iproc.eq.0) write(*,*) 'iorb=',iorb

! coarse grid quantities
        call  loregion_size(nat,rxyz,radii_cf(1,1),crmult,iatype,ntypes,loregion(1,iorb), &
                   hgrid,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)
        if (iproc.eq.0) write(*,'(a,6(i4))') 'coarse grid',nl1,nu1,nl2,nu2,nl3,nu3
        nbox_c(1,1,iorb)=nl1 ; nbox_c(2,1,iorb)=nu1
        nbox_c(1,2,iorb)=nl2 ; nbox_c(2,2,iorb)=nu2
        nbox_c(1,3,iorb)=nl3 ; nbox_c(2,3,iorb)=nu3
        call fill_logrid(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,nat,  &
                         ntypes,iatype,rxyz,radii_cf(1,1),crmult,hgrid,loregion(1,iorb),logrid)
!         if (iproc.eq.0) then
!          write(f4,'(i4.4)') iorb
!          filename = 'grid'//f4//'.ascii'
!          open(unit=22,file=filename,status='unknown')
!          write(22,*) nat
!          write(22,*) alat1,' 0. ',alat2
!          write(22,*) ' 0. ',' 0. ',alat3
!          do iat=1,nat
!           write(22,'(3(1x,e12.5),3x,a20)') rxyz(1,iat),rxyz(2,iat),rxyz(3,iat),atomnames(iatype(iat))
!          enddo
!           write(22,'(3(1x,e12.5),1x,a4)') wan_par(1,iorb),wan_par(2,iorb),wan_par(3,iorb),' W '
! 	  do i3=nl3,nu3 ; do i2=nl2,nu2 ; do i1=nl1,nu1
!           if (logrid(i1,i2,i3)) write(22,'(3(1x,e10.3),1x,a4)') i1*hgrid,i2*hgrid,i3*hgrid,'  g '
!          enddo ; enddo ; enddo 
!         endif
	 call num_segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,mvctr)
         nseg(2*iorb-1)=nseg(2*iorb-2) + mseg
         nvctr(2*iorb-1)=nvctr(2*iorb-2) + mvctr
        if (iproc.eq.0) write(*,*) iorb,'-th orbital has coarse segment, elements',mseg,mvctr

! fine grid quantities
        call  loregion_size(nat,rxyz,radii_cf(1,2),frmult,iatype,ntypes,loregion(1,iorb), &
                   hgrid,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)
        if (iproc.eq.0) write(*,'(a,6(i4))') 'fine   grid',nl1,nu1,nl2,nu2,nl3,nu3
        nbox_f(1,1,iorb)=nl1 ; nbox_f(2,1,iorb)=nu1
        nbox_f(1,2,iorb)=nl2 ; nbox_f(2,2,iorb)=nu2
        nbox_f(1,3,iorb)=nl3 ; nbox_f(2,3,iorb)=nu3
        call fill_logrid(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,nat,  &
                         ntypes,iatype,rxyz,radii_cf(1,2),frmult,hgrid,loregion(1,iorb),logrid)
!         if (iproc.eq.0) then
! 	  do i3=nl3,nu3 ; do i2=nl2,nu2 ; do i1=nl1,nu1
!           if (logrid(i1,i2,i3)) write(22,'(3(1x,e10.3),1x,a4)') i1*hgrid,i2*hgrid,i3*hgrid,'  G '
!          enddo ; enddo ; enddo 
!         endif
	 call num_segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,mvctr)
         nseg(2*iorb)=nseg(2*iorb-1) + mseg
         nvctr(2*iorb)=nvctr(2*iorb-1) + 7*mvctr
        if (iproc.eq.0) write(*,*) iorb,'-th orbital has  fine  segment, elements',mseg,7*mvctr

!        if (iproc.eq.0) close(22)
    enddo

! allocations for arrays holding the wavefunctions and their data descriptors
        allocate(keyg(2,nseg(2*norb)),keyv(nseg(2*norb)))
        write(79,*) 'words for psi and hpsi ',2*nvctr(2*norb)
        allocate(psi(nvctr(2*norb)),hpsi(nvctr(2*norb)))

! now fill the wavefunction descriptor arrays
    do iorb=1,norb

! coarse grid quantities
        nl1=nbox_c(1,1,iorb) ; nu1=nbox_c(2,1,iorb)
        nl2=nbox_c(1,2,iorb) ; nu2=nbox_c(2,2,iorb)
        nl3=nbox_c(1,3,iorb) ; nu3=nbox_c(2,3,iorb)
        call fill_logrid(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,nat,  &
                         ntypes,iatype,rxyz,radii_cf(1,1),crmult,hgrid,loregion(1,iorb),logrid)
         iseg=nseg(2*iorb-2)+1
         mseg=nseg(2*iorb-1)-nseg(2*iorb-2)
	 call segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,keyg(1,iseg),keyv(iseg))

! fine grid quantities
        nl1=nbox_f(1,1,iorb) ; nu1=nbox_f(2,1,iorb)
        nl2=nbox_f(1,2,iorb) ; nu2=nbox_f(2,2,iorb)
        nl3=nbox_f(1,3,iorb) ; nu3=nbox_f(2,3,iorb)
        call fill_logrid(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,nat,  &
                         ntypes,iatype,rxyz,radii_cf(1,2),frmult,hgrid,loregion(1,iorb),logrid)
         iseg=nseg(2*iorb-1)+1
         mseg=nseg(2*iorb)-nseg(2*iorb-1)
         mvctr_f=(nvctr(2*iorb  )-nvctr(2*iorb-1))
         if (abs(dble(mvctr_f/7)-mvctr_f/7.d0).gt.1.d-14) stop 'mvctr_f no multiple of 7'
	 call segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,keyg(1,iseg),keyv(iseg))

      enddo


!       call compresstest(iproc,nproc,n1,n2,n3,norb,nbox_c,nseg,nvctr,keyg,keyv,psi,hpsi)


! determine localization region for all projectors, but do not yet fill the descriptor arrays
    allocate(nboxp_c(2,3,nat),nboxp_f(2,3,nat))
    allocate(nseg_p(0:2*nat))
    allocate(nvctr_p(0:2*nat))
    nseg_p(0)=0 
    nvctr_p(0)=0 

    istart=1
    nproj=0
    do iat=1,nat

    call numb_proj(iatype(iat),ntypes,psppar,mproj)
    if (mproj.ne.0) then 
    
    if (iproc.eq.0) write(*,*) 'projector descriptors for atom with mproj ',iat,mproj
    nproj=nproj+mproj

! coarse grid quantities
        call  loregion_size(1,rxyz(1,iat),radii_cf(1,2),cpmult,iatype(iat),ntypes,.true., &
                   hgrid,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)
        if (iproc.eq.0) write(*,'(a,6(i4))') 'coarse grid',nl1,nu1,nl2,nu2,nl3,nu3
        nboxp_c(1,1,iat)=nl1 ; nboxp_c(2,1,iat)=nu1
        nboxp_c(1,2,iat)=nl2 ; nboxp_c(2,2,iat)=nu2
        nboxp_c(1,3,iat)=nl3 ; nboxp_c(2,3,iat)=nu3
        call fill_logrid(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,1,  &
                         ntypes,iatype(iat),rxyz(1,iat),radii_cf(1,2),cpmult,hgrid,.true.,logrid)
	 call num_segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,mvctr)
         if (iproc.eq.0) write(*,*) 'mseg,mvctr,coarse projectors ',mseg,mvctr
         nseg_p(2*iat-1)=nseg_p(2*iat-2) + mseg
         nvctr_p(2*iat-1)=nvctr_p(2*iat-2) + mvctr
         istart=istart+mvctr*mproj

! fine grid quantities
        call  loregion_size(1,rxyz(1,iat),radii_cf(1,2),fpmult,iatype(iat),ntypes,.true., &
                   hgrid,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)
        if (iproc.eq.0) write(*,'(a,6(i4))') 'fine   grid',nl1,nu1,nl2,nu2,nl3,nu3
        nboxp_f(1,1,iat)=nl1 ; nboxp_f(2,1,iat)=nu1
        nboxp_f(1,2,iat)=nl2 ; nboxp_f(2,2,iat)=nu2
        nboxp_f(1,3,iat)=nl3 ; nboxp_f(2,3,iat)=nu3
        call fill_logrid(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,1,  &
                         ntypes,iatype(iat),rxyz(1,iat),radii_cf(1,2),fpmult,hgrid,.true.,logrid)
	 call num_segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,mvctr)
         if (iproc.eq.0) write(*,*) 'mseg,mvctr, fine  projectors ',mseg,mvctr
         nseg_p(2*iat)=nseg_p(2*iat-1) + mseg
         nvctr_p(2*iat)=nvctr_p(2*iat-1) + mvctr
         istart=istart+7*mvctr*mproj

    else  !(atom has no nonlocal PSP, e.g. H)
         nseg_p(2*iat-1)=nseg_p(2*iat-2) 
         nvctr_p(2*iat-1)=nvctr_p(2*iat-2) 
         nseg_p(2*iat)=nseg_p(2*iat-1) 
         nvctr_p(2*iat)=nvctr_p(2*iat-1) 
    endif
    enddo

        if (iproc.eq.0) write(*,*) 'total number of projectors',nproj
! allocations for arrays holding the projectors and their data descriptors
        allocate(keyg_p(2,nseg_p(2*nat)),keyv_p(nseg_p(2*nat)))
        nprojel=istart-1
        write(79,*) 'words for proj ',nprojel
        allocate(proj(nprojel))
        

! After having determined the size of the projector descriptor arrays fill them
    istart_c=1
    do iat=1,nat
    call numb_proj(iatype(iat),ntypes,psppar,mproj)
    if (mproj.ne.0) then 

! coarse grid quantities
        nl1=nboxp_c(1,1,iat) ; nu1=nboxp_c(2,1,iat)
        nl2=nboxp_c(1,2,iat) ; nu2=nboxp_c(2,2,iat)
        nl3=nboxp_c(1,3,iat) ; nu3=nboxp_c(2,3,iat)
        call fill_logrid(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,1,  &
                         ntypes,iatype(iat),rxyz(1,iat),radii_cf(1,2),cpmult,hgrid,.true.,logrid)

         iseg=nseg_p(2*iat-2)+1
         mseg=nseg_p(2*iat-1)-nseg_p(2*iat-2)
	 call segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
                           logrid,mseg,keyg_p(1,iseg),keyv_p(iseg))

! fine grid quantities
        nl1=nboxp_f(1,1,iat) ; nu1=nboxp_f(2,1,iat)
        nl2=nboxp_f(1,2,iat) ; nu2=nboxp_f(2,2,iat)
        nl3=nboxp_f(1,3,iat) ; nu3=nboxp_f(2,3,iat)
        call fill_logrid(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,1,  &
                         ntypes,iatype(iat),rxyz(1,iat),radii_cf(1,2),fpmult,hgrid,.true.,logrid)
         iseg=nseg_p(2*iat-1)+1
         mseg=nseg_p(2*iat)-nseg_p(2*iat-1)
	 call segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
                           logrid,mseg,keyg_p(1,iseg),keyv_p(iseg))

    endif
    enddo

! Calculate all projectors
    iproj=0
    fpi=(4.d0*atan(1.d0))**(-.75d0)
    do iat=1,nat
     rx=rxyz(1,iat) ; ry=rxyz(2,iat) ; rz=rxyz(3,iat)
     ityp=iatype(iat)

! ONLY GTH PSP form (not HGH)
     do i=1,2
     do j=1,2
!     if (iproc.eq.0) write(*,'(a,3(1x,i3),2(1x,e10.3))') 'iat,i,j,gau_a,ampl',iat,i,j,psppar(i,0,ityp),psppar(i,j,ityp)
     if (psppar(i,j,ityp).ne.0.d0) then
     gau_a=psppar(i,0,ityp)
     do l=1,2*i-1
        if (i.eq.1 .and. j.eq.1) then    ! first s type projector
              nterm=1 ; lx(1)=0 ; ly(1)=0 ; lz(1)=0 
              factor=fpi/(sqrt(gau_a)**3)
        else if (i.eq.1 .and. j.eq.2) then   ! second s type projector
              nterm=3 ; lx(1)=2 ; ly(1)=0 ; lz(1)=0 
                        lx(2)=0 ; ly(2)=2 ; lz(2)=0 
                        lx(3)=0 ; ly(3)=0 ; lz(3)=2 
              factor=sqrt(4.d0/15.d0)*fpi/(sqrt(gau_a)**7)
        else if (i.eq.2 .and. l.eq.1) then  ! px type projector
              nterm=1 ; lx(1)=1 ; ly(1)=0 ; lz(1)=0
              factor=sqrt(2.d0)*fpi/(sqrt(gau_a)**5)
        else if (i.eq.2 .and. l.eq.2) then  ! py type projector
              nterm=1 ; lx(1)=0 ; ly(1)=1 ; lz(1)=0 
              factor=sqrt(2.d0)*fpi/(sqrt(gau_a)**5)
        else if (i.eq.2 .and. l.eq.3) then  ! pz type projector
              nterm=1 ; lx(1)=0 ; ly(1)=0 ; lz(1)=1 
              factor=sqrt(2.d0)*fpi/(sqrt(gau_a)**5)
        else
          stop 'PSP format error'
        endif

        mvctr_c=nvctr_p(2*iat-1)-nvctr_p(2*iat-2)
        mvctr_f=nvctr_p(2*iat  )-nvctr_p(2*iat-1)
        istart_f=istart_c+mvctr_c
        nl1_c=nboxp_c(1,1,iat) ; nu1_c=nboxp_c(2,1,iat)
        nl2_c=nboxp_c(1,2,iat) ; nu2_c=nboxp_c(2,2,iat)
        nl3_c=nboxp_c(1,3,iat) ; nu3_c=nboxp_c(2,3,iat)
        nl1_f=nboxp_f(1,1,iat) ; nu1_f=nboxp_f(2,1,iat)
        nl2_f=nboxp_f(1,2,iat) ; nu2_f=nboxp_f(2,2,iat)
        nl3_f=nboxp_f(1,3,iat) ; nu3_f=nboxp_f(2,3,iat)

        call crtproj(iproc,nterm,n1,n2,n3, & 
                     nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,  & 
                     radii_cf(iatype(iat),2),cpmult,fpmult,hgrid,gau_a,factor,rx,ry,rz,lx,ly,lz, & 
                     mvctr_c,mvctr_f,proj(istart_c),proj(istart_f))
          iproj=iproj+1
! testing
	  call wnrm(mvctr_c,mvctr_f,proj(istart_c),proj(istart_f),scpr)
          if (abs(1.d0-scpr).gt.1.d-1) stop 'norm projector'
! testing end
         istart_c=istart_f+7*mvctr_f
         if (istart_c.gt.istart) stop 'istart_c > istart'

        do iterm=1,nterm
        if (iproc.eq.0) write(*,'(a,i3,1x,a,e10.3,3(i2))') 'projector: iat,atomname,gau_a,lx,ly,lz,err', & 
                             iat,atomnames(iatype(iat)),gau_a,lx(iterm),ly(iterm),lz(iterm)
        enddo


     enddo
     endif
     enddo
     enddo
   enddo
     if (iproj.ne.nproj) stop 'incorrect number of projectors created'
! projector part finished

      deallocate(logrid)


       if (iproc.eq.0) write(*,*) 'Size of real space grids',(2*n1+31),(2*n2+31),(2*n3+31)
! Charge density, Potential in real space
        write(79,*) 'words for rhopot and pot_ion ',2*(2*n1+31)*(2*n2+31)*(2*n3+31)
       allocate(rhopot((2*n1+31)*(2*n2+31)*(2*n3+31)),pot_ion((2*n1+31)*(2*n2+31)*(2*n3+31)))
                 call zero((2*n1+31)*(2*n2+31)*(2*n3+31),pot_ion)
! Allocate and calculate the 1/|r-r'| kernel for the solution of Poisson's equation and test it
       ndegree_ip=8
       if (parallel) then
          call calculate_pardimensions(2*n1+31,2*n2+31,2*n3+31,m1,m2,m3,nf1,nf2,nf3,md1,md2,md3,nfft1,nfft2,nfft3,nproc)
          !call Dimensions_FFT(2*n1+31,2*n2+31,2*n3+31,nfft1,nfft2,nfft3)
          write(*,*) 'dimension of FFT grid',nf1,nf2,nf3
          write(*,*) 'dimension of kernel',nfft1,nfft2,nfft3/nproc
          allocate(pkernel(nfft1*nfft2*nfft3/nproc))
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
          call ParBuild_Kernel(2*n1+31,2*n2+31,2*n3+31,nf1,nf2,nf3,nfft1,nfft2,nfft3, &
               hgridh,ndegree_ip,iproc,nproc,pkernel)
          print *,"kernel built! iproc=",iproc
          call PARtest_kernel(2*n1+31,2*n2+31,2*n3+31,nfft1,nfft2,nfft3,hgridh,pkernel,pot_ion,rhopot,iproc,nproc) 
         
       else
          call Dimensions_FFT(2*n1+31,2*n2+31,2*n3+31,nfft1,nfft2,nfft3)
          write(*,*) 'dimension of FFT grid',nfft1,nfft2,nfft3
          write(*,*) 'dimension of kernel',nfft1/2+1,nfft2/2+1,nfft3/2+1
          allocate(pkernel((nfft1/2+1)*(nfft2/2+1)*(nfft3/2+1)))
          call Build_Kernel(2*n1+31,2*n2+31,2*n3+31,nfft1,nfft2,nfft3, &
               hgridh,ndegree_ip,pkernel)
          
          call test_kernel(2*n1+31,2*n2+31,2*n3+31,nfft1,nfft2,nfft3,hgridh,pkernel,pot_ion,rhopot)
       end if
       
!  Precalculate ionic potential from PSP charge densities and local Gaussian terms
       call input_rho_ion(iproc,ntypes,nat,iatype,atomnames,rxyz,psppar,nelpsp,n1,n2,n3,hgrid,pot_ion,eion)
       if (iproc.eq.0) write(*,*) 'ion-ion interaction energy',eion
       if (parallel) then
          call ParPSolver_Kernel(2*n1+31,2*n2+31,2*n3+31,nfft1,nfft2,nfft3,hgridh,pkernel,.false.,  &
               rhopot,pot_ion,ehart,eexcu,vexcu,iproc,nproc)
          else
       call PSolver_Kernel(2*n1+31,2*n2+31,2*n3+31,nfft1,nfft2,nfft3,hgridh,pkernel,.false.,  &
                          rhopot,pot_ion,ehart,eexcu,vexcu)
       end if

       call addlocgauspsp(iproc,ntypes,nat,iatype,atomnames,rxyz,psppar,n1,n2,n3,hgrid,pot_ion)

! INPUT WAVEFUNCTIONS
       if (calc_inp_wf) then
!	call input_wf_wan(iproc,n1,n2,n3,hgrid,nat,wan_par,ntypes,iatype,psppar, & 
!                          norb,nbox_c,keyg,keyv,nseg,nvctr,psi)
        call crtinpwave(iproc,nproc,n1,n2,n3,nbox_c,nbox_f,norb,  & 
                              wan_par,nseg,nvctr,keyg,keyv,hgrid,psi)
        if (parallel) call commallwaves(iproc,nproc,norb,nvctr,psi)
        do iloewe=1,nloewe
	call loewe(iproc,norb,nvctr,nseg,keyg,keyv,psi)
        enddo
        if (iproc.eq.0) write(*,*) 'input wavefunction has been calculated'
       else
        call readmywaves(iproc,nproc,norb,n1,n2,n3,nbox_c,hgrid,rxyz(1,1),nseg,nvctr,keyg,keyv,psi)
        if (iproc.eq.0) write(*,*) 'Read input wavefunctions from file'
        if (parallel) call commallwaves(iproc,nproc,norb,nvctr,psi)
       endif

        toler=1.d-6
	call checkortho(iproc,toler,norb,nvctr,nseg,keyg,keyv,psi)

       alpha=1.d0
       energy=1.d100
       gnrm=1.d100
! loop for wavefunction minimization
  loop_wavemin: do iter=1,itermax
        if (iproc.eq.0) write(*,*) '-------------------------------------- iter= ',iter

! Potential from electronic charge density
       call sumrho(parallel,iproc,nproc,norb,n1,n2,n3,nbox_c,hgrid,occup,  & 
                       nseg,nvctr,keyg,keyv,psi,rhopot)
       if (parallel) then
          call ParPSolver_Kernel(2*n1+31,2*n2+31,2*n3+31,nfft1,nfft2,nfft3,hgridh,pkernel,.true., & 
               pot_ion,rhopot,ehart,eexcu,vexcu,iproc,nproc)
       else
          call PSolver_Kernel(2*n1+31,2*n2+31,2*n3+31,nfft1,nfft2,nfft3,hgridh,pkernel,.true., & 
               pot_ion,rhopot,ehart,eexcu,vexcu)
       end if


! local potential and kinetic energy for all orbitals belonging to iproc
        call applylocpotkin(iproc,nproc,norb,n1,n2,n3,nbox_c,hgrid,occup,  & 
                   nseg,nvctr,keyg,keyv,psi,rhopot,hpsi,epot_sum,ekin_sum)
 
! apply all PSP projectors for all orbitals belonging to iproc
        call applyprojectors(iproc,nproc,ntypes,nat,iatype,psppar,occup, &
                    nprojel,nproj,nseg_p,keyg_p,keyv_p,nvctr_p,proj,  &
                    norb,nseg,keyg,keyv,nvctr,psi,hpsi,eproj_sum)

       energy_p=ekin_sum+epot_sum+eproj_sum

! check whether CPU time exceeded
         open(unit=55,file='CPUlimit',status='unknown')
         read(55,*,end=555) cpulimit
         close(55)
         call cpu_time(tcpu2)
           if (tcpu2-tcpu1.ge.cpulimit) then
           if (iproc.eq.0) write(6,*) 'CPU time exceeded'
!            call system("rm CPUlimit")
           exit
         endif
555      continue
         close(55)

       if (gnrm.le.gnrm_cv .or. iter.eq.itermax) exit

! Apply  orthogonality constraints to all orbitals belonging to iproc
        call  orthoconstraint(iproc,nproc,norb,occup,  & 
                   nseg,nvctr,keyg,keyv,psi,hpsi,gnrm_p,scprsum_p)

! Preconditions all orbitals belonging to iproc
        call preconditionall(iproc,nproc,norb,n1,n2,n3,nbox_c,hgrid,  & 
                   nseg,nvctr,keyg,keyv,hpsi)

       energyold=energy
       if (parallel) then
       call cpu_time(tr0)
       call system_clock(count1,count_rate,count_max)
       wrkallred(1,2)=energy_p ; wrkallred(2,2)=gnrm_p ; wrkallred(3,2)=scprsum_p
       call MPI_ALLREDUCE(wrkallred(1,2),wrkallred(1,1),3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
       energy=wrkallred(1,1) ; gnrm=wrkallred(2,1) ; scprsum=wrkallred(3,1)
       call cpu_time(tr1)
       call system_clock(count2,count_rate,count_max)
       tel=dble(count2-count1)/dble(count_rate)
       write(78,*) 'ENERGY: ALLREDUCE TIME',iproc,tr1-tr0,tel
       write(78,*) '----------------------------------------------'
       else
       energy=energy_p
       gnrm=gnrm_p
       scprsum=scprsum_p
       endif
       energybs=energy
       energy=energy-ehart+eexcu-vexcu+eion
       gnrm=sqrt(gnrm)


! update wavefunction
       if (energy.gt.energyold) then
       alpha=max(.125d0,.5d0*alpha)
       if (alpha.eq..125d0) write(*,*) 'Convergence problem or limit'
       else
       alpha=min(1.05d0*alpha,2.d0)
       endif
       if (iproc.eq.0) write(*,*) 'alpha=',alpha

! update all wavefunctions (belonging to iproc)  with the preconditioned gradient
       call update(iproc,nproc,norb,alpha,nseg,nvctr,keyg,keyv,hpsi,psi)

! All processors need now the updated psi
       if (parallel) call commallwaves(iproc,nproc,norb,nvctr,psi)

       do iloewe=1,nloewe
       call loewe(iproc,norb,nvctr,nseg,keyg,keyv,psi)
       enddo
       toler=1.d-6
       call checkortho(iproc,toler,norb,nvctr,nseg,keyg,keyv,psi)

       tt=energybs-scprsum
       if (abs(tt).gt.1.d-10) then 
          write(*,*) 'ERROR: inconsistency between gradient and energy',tt
       endif
       if (iproc.eq.0) then
       write(*,'(a,3(1x,e18.11))') 'ekin_sum,epot_sum,eproj_sum',  & 
                                    ekin_sum,epot_sum,eproj_sum
       write(*,'(a,3(1x,e18.11))') 'ehart,eexcu,vexcu',ehart,eexcu,vexcu
       write(*,'(a,i6,2x,e19.12,1x,e9.2))') 'iter,total energy,gnrm',iter,energy,gnrm
       endif


   enddo loop_wavemin

       write(*,'(a,i3,3(1x,e18.11))') 'iproc,ehart,eexcu,vexcu',iproc,ehart,eexcu,vexcu
     if (iproc.eq.0) then
       write(*,*) 'FINAL RESULT'
       write(*,'(a,3(1x,e18.11))') 'ekin_sum,epot_sum,eproj_sum',ekin_sum,epot_sum,eproj_sum
       write(*,'(a,3(1x,e18.11))') 'ehart,eexcu,vexcu',ehart,eexcu,vexcu
       write(*,'(a,i6,2x,e19.12,1x,e9.2))') 'iter,total energy,gnrm',iter,energy,gnrm
     endif
       call cpu_time(tcpu2)
       write(6,*) 'total CPU time', tcpu2-tcpu1
       write(78,*) 'total CPU time', tcpu2-tcpu1

! transform to KS orbitals
!       call KStrans(parallel,iproc,nproc,norb,nvctr_c,nvctr_f,occup,  & 
!                          hpsi_c,hpsi_f,psi_c,psi_f,evsum)
!       if (abs(evsum-energybs).gt.1.d-8) write(*,*) 'Difference:evsum,energybs',evsum,energybs
!
!
!  write all the wavefunctions into files
      write(*,*) iproc,' start writing waves'
        call  writemywaves(iproc,nproc,norb,n1,n2,n3,nbox_c,hgrid,  & 
                   rxyz(1,1),nseg,nvctr,keyg,keyv,psi)
      write(*,*) iproc,' finished writing waves'
        deallocate(wan_par)
        deallocate(pkernel)
        deallocate(keyg_p,keyv_p,proj)
        deallocate(rhopot,pot_ion)
        deallocate(occup,loregion)
        deallocate(nvctr,nseg)
        deallocate(nvctr_p,nseg_p)
        deallocate(psi,hpsi)
        deallocate(keyg,keyv,nbox_c,nbox_f)
        deallocate(psppar,nelpsp,radii_cf,rad_cov)

	end


       subroutine commallwaves(iproc,nproc,norb,nvctr,psi)
! Before calling this routine each processor has only the latest version of the 
! orbitals on which it works, afterwards the updated version of all orbitals
       implicit real*8 (a-h,o-z)
        integer count1,count2,count_rate,count_max
       integer recvcounts,displs,sendcount
! automatic arrays
       dimension recvcounts(0:nproc-1),displs(0:nproc-1)
        dimension nvctr(0:2*norb),psi(nvctr(2*norb))
        real*8, allocatable :: tpsi(:)
       parameter(eps_mach=1.d-12)
       include 'mpif.h'

       allocate(tpsi(nvctr(2*norb)))

        onem=1.d0-eps_mach
        norb_p=int(onem+dble(norb)/dble(nproc))
      do jproc=0,nproc-1
         iaorb=min(jproc*norb_p+1,norb)
         iborb=min((jproc+1)*norb_p+1,norb+1)
!         write(*,*) 'jproc,iaorb,iborb',jproc,iaorb,iborb,nvctr(2*iborb-2),nvctr(2*iaorb-2)
!         displs(jproc)=nvctr(2*iaorb-2)+1
         displs(jproc)=nvctr(2*iaorb-2)
         recvcounts(jproc)=nvctr(2*iborb-2)-nvctr(2*iaorb-2)
      enddo
!      ipsi=displs(iproc)
      ipsi=displs(iproc)+1
      sendcount=recvcounts(iproc)

!     write(*,*) 'iproc',iproc
!     write(*,*) 'displs',displs,ipsi
!     write(*,*) 'recvcounts',recvcounts,sendcount
!     write(*,*) 'values',psi(ipsi),psi(ipsi+sendcount-1)

     call cpu_time(tb0)
     call system_clock(count1,count_rate,count_max)

      call MPI_ALLGATHERV(psi(ipsi),sendcount,MPI_DOUBLE_PRECISION,  &
                          tpsi,recvcounts,displs, MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

!      write(*,*) iproc,' ALLGATHERV done'
      do i=1,nvctr(2*norb)
      psi(i)=tpsi(i)
      enddo

     call cpu_time(tb1)
     call system_clock(count2,count_rate,count_max)
     tel=dble(count2-count1)/dble(count_rate)
     write(78,*) 'ALLGATHERV TIME for all wavefunctions',iproc,tb1-tb0,tel
     write(78,*) '-------------------------------------------------------'

     deallocate(tpsi)

     return
     end




       subroutine commallwaves_in_place(iproc,nproc,norb,nvctr,psi)
! Requires full MPI2 implementation
! Before calling this routine each processor has only the latest version of the 
! orbitals on which it works, afterwards the updated version of all orbitals
       implicit real*8 (a-h,o-z)
        integer count1,count2,count_rate,count_max
       integer recvcounts,displs,sendcount
! automatic arrays
       dimension recvcounts(0:nproc-1),displs(0:nproc-1)
        dimension nvctr(0:2*norb),psi(nvctr(2*norb))
       parameter(eps_mach=1.d-12)
       include 'mpif.h'


        onem=1.d0-eps_mach
        norb_p=int(onem+dble(norb)/dble(nproc))
      do jproc=0,nproc-1
         iaorb=min(jproc*norb_p+1,norb)
         iborb=min((jproc+1)*norb_p+1,norb+1)
!         write(*,*) 'jproc,iaorb,iborb',jproc,iaorb,iborb,nvctr(2*iborb-2),nvctr(2*iaorb-2)
!         displs(jproc)=nvctr(2*iaorb-2)+1
         displs(jproc)=nvctr(2*iaorb-2)
         recvcounts(jproc)=nvctr(2*iborb-2)-nvctr(2*iaorb-2)
      enddo
!      ipsi=displs(iproc)
      ipsi=displs(iproc)+1
      sendcount=recvcounts(iproc)

!     write(*,*) 'iproc',iproc
!     write(*,*) 'displs',displs,ipsi
!     write(*,*) 'recvcounts',recvcounts,sendcount
!     write(*,*) 'values',psi(ipsi),psi(ipsi+sendcount-1)

     call cpu_time(tb0)
     call system_clock(count1,count_rate,count_max)

!      call MPI_ALLGATHERV(psi(ipsi),sendcount,MPI_DOUBLE_PRECISION,  &
!                          tpsi,recvcounts,displs, MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

!      write(*,*) 'MPI_IN_PLACE',MPI_IN_PLACE
      call MPI_ALLGATHERV(MPI_IN_PLACE,sendcount,MPI_DOUBLE_PRECISION,  &
                          psi,recvcounts,displs, MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

     call cpu_time(tb1)
     call system_clock(count2,count_rate,count_max)
     tel=dble(count2-count1)/dble(count_rate)
     write(78,*) 'ALLGATHERV TIME for all wavefunctions',iproc,tb1-tb0,tel
     write(78,*) '-------------------------------------------------------'

     return
     end


