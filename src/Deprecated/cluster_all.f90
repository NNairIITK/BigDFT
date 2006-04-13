

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
! does an electronic structure calculation. Output is the total energy and the forces
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


       call compresstest(iproc,nproc,n1,n2,n3,norb,nbox_c,nseg,nvctr,keyg,keyv,psi,hpsi)

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
       call chargedens(parallel,iproc,nproc,norb,n1,n2,n3,nbox_c,hgrid,occup,  & 
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


        subroutine writemywaves(iproc,nproc,norb,n1,n2,n3,nbox_c,hgrid,  & 
                   rxyz,nseg,nvctr,keyg,keyv,psi)
        implicit real*8 (a-h,o-z)
        parameter(eps_mach=1.d-12)
        dimension nbox_c(2,3,norb),rxyz(3)
        dimension nseg(0:2*norb),nvctr(0:2*norb)
        dimension keyg(2,nseg(2*norb)),keyv(nseg(2*norb))
        dimension psi(nvctr(2*norb))

       onem=1.d0-eps_mach
       norb_p=int(onem+dble(norb)/dble(nproc))
       do iorb=iproc*norb_p+1,min((iproc+1)*norb_p,norb)
         mseg_c=nseg(2*iorb-1)-nseg(2*iorb-2) 
         mseg_f=nseg(2*iorb  )-nseg(2*iorb-1)
         iseg_c=nseg(2*iorb-2)+1
         iseg_f=nseg(2*iorb-1)+1
         mvctr_c= nvctr(2*iorb-1)-nvctr(2*iorb-2) 
         mvctr_f=(nvctr(2*iorb  )-nvctr(2*iorb-1))/7
         ipsi_c=nvctr(2*iorb-2)+1
         ipsi_f=nvctr(2*iorb-1)+1
         nl1=nbox_c(1,1,iorb) ; nu1=nbox_c(2,1,iorb)
         nl2=nbox_c(1,2,iorb) ; nu2=nbox_c(2,2,iorb)
         nl3=nbox_c(1,3,iorb) ; nu3=nbox_c(2,3,iorb)

       call writeonewave(iorb,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,hgrid,rxyz,  & 
                         mseg_c,mvctr_c,keyg(1,iseg_c),keyv(iseg_c)  & 
                        ,mseg_f,mvctr_f,keyg(1,iseg_f),keyv(iseg_f), & 
                             psi(ipsi_c),psi(ipsi_f))
       enddo
       return
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


        subroutine orthoconstraint(iproc,nproc,norb,occup,  & 
                   nseg,nvctr,keyg,keyv,psi,hpsi,gnrm_p,scprsum_p)
!Effect of orthogonality constraints on gradient and then preconditioning
        implicit real*8 (a-h,o-z)
        parameter(eps_mach=1.d-12)
        dimension nseg(0:2*norb),nvctr(0:2*norb)
        dimension keyg(2,nseg(2*norb)),keyv(nseg(2*norb))
        dimension psi(nvctr(2*norb)),hpsi(nvctr(2*norb)),occup(norb)

        onem=1.d0-eps_mach
        norb_p=int(onem+dble(norb)/dble(nproc))

     gnrm_p=0.d0
     scprsum_p=0.d0
     do iorb=iproc*norb_p+1,min((iproc+1)*norb_p,norb)
       maseg_c=nseg(2*iorb-1)-nseg(2*iorb-2)
       maseg_f=nseg(2*iorb  )-nseg(2*iorb-1)
       iseg_c=nseg(2*iorb-2)+1
       iseg_f=nseg(2*iorb-1)+1
       mavctr_c= nvctr(2*iorb-1)-nvctr(2*iorb-2)
       mavctr_f=(nvctr(2*iorb  )-nvctr(2*iorb-1))/7
       ipsi_c=nvctr(2*iorb-2)+1
       ipsi_f=nvctr(2*iorb-1)+1

         do jorb=1,norb
           mbseg_c=nseg(2*jorb-1)-nseg(2*jorb-2)
           mbseg_f=nseg(2*jorb  )-nseg(2*jorb-1)
           jseg_c=nseg(2*jorb-2)+1
           jseg_f=nseg(2*jorb-1)+1
           mbvctr_c= nvctr(2*jorb-1)-nvctr(2*jorb-2)
           mbvctr_f=(nvctr(2*jorb  )-nvctr(2*jorb-1))/7
           jpsi_c=nvctr(2*jorb-2)+1
           jpsi_f=nvctr(2*jorb-1)+1

	call wdot(  & 
                   mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv(jseg_c),keyv(jseg_f),  & 
                   keyg(1,jseg_c),keyg(1,jseg_f),psi(jpsi_c),psi(jpsi_f),  &
                   mavctr_c,mavctr_f,maseg_c,maseg_f,keyv(iseg_c),keyv(iseg_f),  & 
                   keyg(1,iseg_c),keyg(1,iseg_f),hpsi(ipsi_c),hpsi(ipsi_f),scpr)

         if (jorb.eq.iorb) scprsum_p=scprsum_p+occup(iorb)*scpr
	call waxpy(  & 
                  -scpr,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv(jseg_c),keyv(jseg_f),  & 
                   keyg(1,jseg_c),keyg(1,jseg_f),psi(jpsi_c),psi(jpsi_f), &
                   mavctr_c,mavctr_f,maseg_c,maseg_f,keyv(iseg_c),keyv(iseg_f),  & 
                   keyg(1,iseg_c),keyg(1,iseg_f),hpsi(ipsi_c),hpsi(ipsi_f))
       enddo
       call wnrm(mavctr_c,mavctr_f,hpsi(ipsi_c),hpsi(ipsi_f),scpr) 
       gnrm_p=gnrm_p+scpr
     enddo

     return
     end


        subroutine update(iproc,nproc,norb,alpha,nseg,nvctr,keyg,keyv,hpsi,psi)
!Effect of orthogonality constraints on gradient and then preconditioning
        implicit real*8 (a-h,o-z)
        parameter(eps_mach=1.d-12)
        dimension nseg(0:2*norb),nvctr(0:2*norb)
        dimension keyg(2,nseg(2*norb)),keyv(nseg(2*norb))
        dimension hpsi(nvctr(2*norb)),psi(nvctr(2*norb))

        onem=1.d0-eps_mach
        norb_p=int(onem+dble(norb)/dble(nproc))

     do iorb=iproc*norb_p+1,min((iproc+1)*norb_p,norb)
         mseg_c=nseg(2*iorb-1)-nseg(2*iorb-2) 
         mseg_f=nseg(2*iorb  )-nseg(2*iorb-1)
         iseg_c=nseg(2*iorb-2)+1
         iseg_f=nseg(2*iorb-1)+1
         mvctr_c= nvctr(2*iorb-1)-nvctr(2*iorb-2) 
         mvctr_f=(nvctr(2*iorb  )-nvctr(2*iorb-1))/7
         ipsi_c=nvctr(2*iorb-2)+1
         ipsi_f=nvctr(2*iorb-1)+1

!	call waxpy(  & 
!                  -alpha,mvctr_c,mvctr_f,mseg_c,mseg_f,keyv(iseg_c),keyv(iseg_f),  & 
!                   keyg(1,iseg_c),keyg(1,iseg_f),hpsi(ipsi_c),hpsi(ipsi_f), &
!                   mvctr_c,mvctr_f,mseg_c,mseg_f,keyv(iseg_c),keyv(iseg_f),  & 
!                   keyg(1,iseg_c),keyg(1,iseg_f),psi(ipsi_c),psi(ipsi_f))

        call DAXPY(mvctr_c,-alpha,hpsi(ipsi_c),1,psi(ipsi_c),1)
        call DAXPY(7*mvctr_f,-alpha,hpsi(ipsi_f),1,psi(ipsi_f),1)
     enddo

     return
     end


        subroutine preconditionall(iproc,nproc,norb,n1,n2,n3,nbox_c,hgrid,  & 
                   nseg,nvctr,keyg,keyv,hpsi)
!Effect of orthogonality constraints on gradient and then preconditioning
        implicit real*8 (a-h,o-z)
        parameter(eps_mach=1.d-12)
        dimension nseg(0:2*norb),nvctr(0:2*norb)
        dimension keyg(2,nseg(2*norb)),keyv(nseg(2*norb))
        dimension hpsi(nvctr(2*norb)),nbox_c(2,3,norb)

        onem=1.d0-eps_mach
        norb_p=int(onem+dble(norb)/dble(nproc))

     do iorb=iproc*norb_p+1,min((iproc+1)*norb_p,norb)
         mseg_c=nseg(2*iorb-1)-nseg(2*iorb-2) 
         mseg_f=nseg(2*iorb  )-nseg(2*iorb-1)
         iseg_c=nseg(2*iorb-2)+1
         iseg_f=nseg(2*iorb-1)+1
         mvctr_c= nvctr(2*iorb-1)-nvctr(2*iorb-2) 
         mvctr_f=(nvctr(2*iorb  )-nvctr(2*iorb-1))/7
         ipsi_c=nvctr(2*iorb-2)+1
         ipsi_f=nvctr(2*iorb-1)+1
         nl1=nbox_c(1,1,iorb) ; nu1=nbox_c(2,1,iorb)
         nl2=nbox_c(1,2,iorb) ; nu2=nbox_c(2,2,iorb)
         nl3=nbox_c(1,3,iorb) ; nu3=nbox_c(2,3,iorb)

!       call precondition_simple(mvctr_c,mvctr_f,hpsi(ipsi_c),hpsi(ipsi_f))
!The whole bo is used instead of subbox
       call precondition(n1,n2,n3,hgrid,mseg_c,mvctr_c,mvctr_f,  & 
            keyg(1,iseg_c),keyv(iseg_c),hpsi(ipsi_c),hpsi(ipsi_f))

     enddo

     return
     end


     subroutine precondition_simple(mvctr_c,mvctr_f,hpsi_c,hpsi_f)
     implicit real*8 (a-h,o-z)
     dimension hpsi_c(mvctr_c),hpsi_f(7,mvctr_f)

      small=5.d-3
      do i=1,mvctr_c
        hpsi_c(i)=small*hpsi_c(i)
      enddo

      do i=1,mvctr_f
      do l=1,7
        hpsi_f(l,i)=small*hpsi_f(l,i)
      enddo
      enddo

     return
     end


	subroutine wannier_par(iproc,nat,norb,rxyz,ntypes,iatype,rad_cov,wan_par)
        implicit real*8 (a-h,o-z)
        dimension rxyz(3,nat),iatype(nat),rad_cov(ntypes),wan_par(0:3,norb)
! Wannier centers are in between any pair of atoms whose distance is less than the sum of the covalent radii

        iorb=0
	do iat=1,nat
        ityp=iatype(iat)
	do jat=1,iat-1
        jtyp=iatype(jat)
        dist2=(rxyz(1,iat)-rxyz(1,jat))**2 + (rxyz(2,iat)-rxyz(2,jat))**2+ & 
              (rxyz(3,iat)-rxyz(3,jat))**2
! add a safety margin of 20 percent to the covalent radii
        rcut2=(1.2d0*(rad_cov(ityp)+rad_cov(jtyp)))**2
 
        if (dist2.lt.rcut2) then 
           iorb=iorb+1
           if (iorb.gt.norb) stop 'wannier_par, iorb > norb'
           wan_par(1,iorb)=.5d0*(rxyz(1,iat)+rxyz(1,jat))
           wan_par(2,iorb)=.5d0*(rxyz(2,iat)+rxyz(2,jat))
           wan_par(3,iorb)=.5d0*(rxyz(3,iat)+rxyz(3,jat))
           wan_par(0,iorb)=max(rad_cov(ityp),rad_cov(jtyp))
!           if (iorb.eq.norb) goto 1232
        endif
        enddo
        enddo

!1232    continue
        if (iproc.eq.0) then 
        write(*,*) ' Wannier centers and their radii'
        do jorb=1,iorb
        write(*,'(i4,4(1x,e10.3))')  & 
             jorb,wan_par(1,jorb),wan_par(2,jorb),wan_par(3,jorb),wan_par(0,jorb)
        enddo
        endif
        if (iorb.ne.norb) stop 'wannier_par, iorb <> norb'

        return
	end


        subroutine localizationregion(iproc,nat,norb,rxyz,wan_par,radlocmult,loregion)
        implicit real*8 (a-h,o-z)
        logical loregion
        dimension rxyz(3,nat),wan_par(0:3,norb),loregion(nat,norb)

        do iorb=1,norb
          ic=0
          do iat=1,nat
            dd=(rxyz(1,iat)-wan_par(1,iorb))**2+(rxyz(2,iat)-wan_par(2,iorb))**2  &
                                               +(rxyz(3,iat)-wan_par(3,iorb))**2 
            if (dd.le.(radlocmult*wan_par(0,iorb))**2) then
              ic=ic+1
              loregion(iat,iorb)=.true.
            else
              loregion(iat,iorb)=.false.
            endif
          enddo
          if (iproc.eq.0) write(*,'(a,i4,a,i3,a)')   & 
             'localization region of orbital ',iorb,' consists of ',ic,' atom centered spheres'
          if (ic.lt.1) stop 'vanishing localization region'
        enddo

        return
        end

        subroutine system_size(nat,rxyz,radii,rmult,iatype,ntypes, &
                   cxmin,cxmax,cymin,cymax,czmin,czmax)
        implicit real*8 (a-h,o-z)
        parameter(eps_mach=1.d-12)
        dimension rxyz(3,nat),radii(ntypes),iatype(nat)

        cxmax=-1.d100 ; cxmin=1.d100
        cymax=-1.d100 ; cymin=1.d100
        czmax=-1.d100 ; czmin=1.d100
        do iat=1,nat
            rad=radii(iatype(iat))*rmult
            cxmax=max(cxmax,rxyz(1,iat)+rad) ; cxmin=min(cxmin,rxyz(1,iat)-rad)
            cymax=max(cymax,rxyz(2,iat)+rad) ; cymin=min(cymin,rxyz(2,iat)-rad)
            czmax=max(czmax,rxyz(3,iat)+rad) ; czmin=min(czmin,rxyz(3,iat)-rad)
        enddo
  
      cxmax=cxmax-eps_mach ; cxmin=cxmin+eps_mach
      cymax=cymax-eps_mach ; cymin=cymin+eps_mach
      czmax=czmax-eps_mach ; czmin=czmin+eps_mach

        return
        end


        subroutine loregion_size(nat,rxyz,radii,rmult,iatype,ntypes,loregion, &
                   hgrid,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)
        implicit real*8 (a-h,o-z)
        logical loregion
        parameter(eps_mach=1.d-12)
        dimension rxyz(3,nat),radii(ntypes),iatype(nat),loregion(nat)

        cxmax=-1.d100 ; cxmin=1.d100
        cymax=-1.d100 ; cymin=1.d100
        czmax=-1.d100 ; czmin=1.d100
        do iat=1,nat
        if (loregion(iat)) then
            rad=radii(iatype(iat))*rmult
            cxmax=max(cxmax,rxyz(1,iat)+rad) ; cxmin=min(cxmin,rxyz(1,iat)-rad)
            cymax=max(cymax,rxyz(2,iat)+rad) ; cymin=min(cymin,rxyz(2,iat)-rad)
            czmax=max(czmax,rxyz(3,iat)+rad) ; czmin=min(czmin,rxyz(3,iat)-rad)
!        write(*,*) radii(iatype(iat)),rmult
!        write(*,*) rxyz(1,iat),rxyz(2,iat),rxyz(3,iat)
!        write(*,*) 'loregion_size',cxmin,cxmax
!        write(*,*) '             ',cymin,cymax
!        write(*,*) '             ',czmin,czmax
        endif
        enddo
  
      cxmax=cxmax-eps_mach ; cxmin=cxmin+eps_mach
      cymax=cymax-eps_mach ; cymin=cymin+eps_mach
      czmax=czmax-eps_mach ; czmin=czmin+eps_mach
      onem=1.d0-eps_mach
      nl1=int(onem+cxmin/hgrid)   
      nl2=int(onem+cymin/hgrid)   
      nl3=int(onem+czmin/hgrid)   
      nu1=int(cxmax/hgrid)  
      nu2=int(cymax/hgrid)  
      nu3=int(czmax/hgrid)  
!        write(*,'(a,6(i4))') 'loregion_size',nl1,nu1,nl2,nu2,nl3,nu3
!        write(*,*) 'loregion_size',cxmin,cxmax
!        write(*,*) '             ',cymin,cymax
!        write(*,*) '             ',czmin,czmax
      if (nl1.lt.0)   stop 'nl1: localization region outside cell'
      if (nl2.lt.0)   stop 'nl2: localization region outside cell'
      if (nl3.lt.0)   stop 'nl3: localization region outside cell'
      if (nu1.gt.n1)   stop 'nu1: localization region outside cell'
      if (nu2.gt.n2)   stop 'nu2: localization region outside cell'
      if (nu3.gt.n3)   stop 'nu3: localization region outside cell'

        return
        end





        subroutine readmywaves(iproc,nproc,norb,n1,n2,n3,nbox_c,hgrid,rxyz,nseg,nvctr,keyg,keyv,psi)
! reads wavefunction from file and transforms it properly if hgrid or size of simulation cell have changed
        implicit real*8 (a-h,o-z)
        character*30 filename
        character*4 f4
        logical cif1,cif2,cif3
        parameter(eps_mach=1.d-12)
        dimension nseg(0:2*norb),nvctr(0:2*norb)
        dimension keyg(2,nseg(2*norb)),keyv(nseg(2*norb))
        dimension psi(nvctr(2*norb))
        dimension xya(-1:1,-1:1),xa(-1:1)
        dimension rxyz(3),rxyz_old(3),nbox_c(2,3,norb)
        allocatable :: psifscf(:,:,:)
        allocatable :: psigold(:,:,:,:,:,:),psifscfold(:,:,:),psifscfoex(:,:,:),wwold(:)

        onem=1.d0-eps_mach
        norb_p=int(onem+dble(norb)/dble(nproc))

     do 100,iorb=iproc*norb_p+1,min((iproc+1)*norb_p,norb)
         mseg_c=nseg(2*iorb-1)-nseg(2*iorb-2) 
         mseg_f=nseg(2*iorb  )-nseg(2*iorb-1)
         iseg_c=nseg(2*iorb-2)+1
         iseg_f=nseg(2*iorb-1)+1
         mvctr_c= nvctr(2*iorb-1)-nvctr(2*iorb-2) 
         mvctr_f=(nvctr(2*iorb  )-nvctr(2*iorb-1))/7
         ipsi_c=nvctr(2*iorb-2)+1
         ipsi_f=nvctr(2*iorb-1)+1
         nl1=nbox_c(1,1,iorb) ; nu1=nbox_c(2,1,iorb)
         nl2=nbox_c(1,2,iorb) ; nu2=nbox_c(2,2,iorb)
         nl3=nbox_c(1,3,iorb) ; nu3=nbox_c(2,3,iorb)


        write(f4,'(i4.4)') iorb
        filename = 'wavefunction.'//f4
        open(unit=99,file=filename,status='unknown')

         read(99,*) iorbold
         if (iorbold.ne.iorb) stop 'readallwaves'
         read(99,*) hgridold
         read(99,*) nold1,nold2,nold3
         read(99,*) nlold1,nuold1,nlold2,nuold2,nlold3,nuold3
         read(99,*) (rxyz_old(j),j=1,3)
         read(99,*) mvctrold_c, mvctrold_f

      if (hgridold.eq. hgrid .and. mvctrold_c.eq.mvctr_c .and. mvctrold_f.eq.mvctr_f  & 
          .and. nold1.eq.n1  .and. nold2.eq.n2 .and. nold3.eq.n3  & 
          .and. nlold1.eq.nl1 .and. nuold1.eq.nu1  & 
          .and. nlold2.eq.nl2 .and. nuold2.eq.nu2  &  
          .and. nlold3.eq.nl3 .and. nuold3.eq.nu3) then

      if (iproc.eq.0) write(*,*) 'wavefunction ',iorb,' needs NO transformation'
	do j=1,mvctrold_c
            read(99,*) i1,i2,i3,tt
            psi(ipsi_c+j-1)=tt
         enddo
	do j=1,7*mvctrold_f-6,7
            read(99,*) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
            psi(ipsi_f+j-1+0)=t1
            psi(ipsi_f+j-1+1)=t2
            psi(ipsi_f+j-1+2)=t3
            psi(ipsi_f+j-1+3)=t4
            psi(ipsi_f+j-1+4)=t5
            psi(ipsi_f+j-1+5)=t6
            psi(ipsi_f+j-1+6)=t7
         enddo

       else
       if (iproc.eq.0) write(*,*) 'wavefunction ',iorb,' is transformed'

         allocate(psifscf(-7+2*nl1:2*nu1+8,-7+2*nl2:2*nu2+8,-7+2*nl3:2*nu3+8))

         allocate(psigold(nlold1:nuold1,2,nlold2:nuold2,2,nlold3:nuold3,2),  & 
               psifscfold(-7+2*nlold1:2*nuold1+8,-7+2*nlold2:2*nuold2+8,-7+2*nlold3:2*nuold3+8), &
               psifscfoex(-8+2*nlold1:2*nuold1+9,-8+2*nlold2:2*nuold2+9,-8+2*nlold3:2*nuold3+9), &
                   wwold((2*(nuold1-nlold1)+16)*(2*(nold2-nlold2)+16)*(2*(nold3-nlold3)+16)))

         call zero(8*(nuold1-nlold1+1)*(nuold2-nlold2+1)*(nuold3-nlold3+1),psigold)

	do iel=1,mvctrold_c
            read(99,*) i1,i2,i3,tt
            psigold(i1,1,i2,1,i3,1)=tt
         enddo
	do iel=1,mvctrold_f
            read(99,*) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
            psigold(i1,2,i2,1,i3,1)=t1
            psigold(i1,1,i2,2,i3,1)=t2
            psigold(i1,2,i2,2,i3,1)=t3
            psigold(i1,1,i2,1,i3,2)=t4
            psigold(i1,2,i2,1,i3,2)=t5
            psigold(i1,1,i2,2,i3,2)=t6
            psigold(i1,2,i2,2,i3,2)=t7
         enddo

! calculate fine scaling functions
	call synthese_grow(nuold1-nlold1,nuold2-nlold2,nuold3-nlold3,wwold,psigold,psifscfold)

         do i3=-7+2*nlold3,2*nuold3+8
         do i2=-7+2*nlold2,2*nuold2+8
           i1=-8+2*nlold1
           psifscfoex(i1,i2,i3)=0.d0
           do i1=-7+2*nlold1,2*nuold1+8
             psifscfoex(i1,i2,i3)=psifscfold(i1,i2,i3)
           enddo
           i1=2*nuold1+9
           psifscfoex(i1,i2,i3)=0.d0
         enddo
         enddo

         i3=-8+2*nlold3
         do i2=-8+2*nlold2,2*nuold2+9
         do i1=-8+2*nlold1,2*nuold1+9
           psifscfoex(i1,i2,i3)=0.d0
         enddo
         enddo
         i3=2*nuold3+9
         do i2=-8+2*nlold2,2*nuold2+9
         do i1=-8+2*nlold1,2*nuold1+9
           psifscfoex(i1,i2,i3)=0.d0
         enddo
         enddo

         i2=-8+2*nlold2
         do i3=-8+2*nlold3,2*nuold3+9
         do i1=-8+2*nlold1,2*nuold1+9
           psifscfoex(i1,i2,i3)=0.d0
         enddo
         enddo
         i2=2*nuold2+9
         do i3=-8+2*nlold3,2*nuold3+9
         do i1=-8+2*nlold1,2*nuold1+9
           psifscfoex(i1,i2,i3)=0.d0
         enddo
         enddo

! transform to new structure	
        dx=rxyz(1)-rxyz_old(1)
        dy=rxyz(2)-rxyz_old(2)
        dz=rxyz(3)-rxyz_old(3)
        write(*,*) 'dxyz',dx,dy,dz
        hgridh=.5d0*hgrid
        hgridhold=.5d0*hgridold
        call zero((2*(nu1-nl1)+16)*(2*(nu2-nl2)+16)*(2*(nu3-nl3)+16),psifscf)
        do i3=-7+2*nl3,2*nu3+8
        z=i3*hgridh
        j3=nint((z-dz)/hgridhold)
        cif3=(j3.ge.-7+2*nlold3 .and. j3.le.2*nuold3+8)
        do i2=-7+2*nl2,2*nu2+8
        y=i2*hgridh
        j2=nint((y-dy)/hgridhold)
        cif2=(j2.ge.-7+2*nlold2 .and. j2.le.2*nuold2+8)
        do i1=-7+2*nl1,2*nu1+8
        x=i1*hgridh
        j1=nint((x-dx)/hgridhold)
        cif1=(j1.ge.-72*nlold1 .and. j1.le.2*nuold1+8)

!        if (cif1 .and. cif2 .and. cif3) psifscf(i1,i2,i3)=psifscfold(j1,j2,j3)
!        if (cif1 .and. cif2 .and. cif3) psifscf(i1,i2,i3)=psifscfoex(j1,j2,j3)

        if (cif1 .and. cif2 .and. cif3) then 
        zr = ((z-dz)-j3*hgridhold)/hgridhold
        do l2=-1,1
        do l1=-1,1
        ym1=psifscfoex(j1+l1,j2+l2,j3-1)
        y00=psifscfoex(j1+l1,j2+l2,j3  )
        yp1=psifscfoex(j1+l1,j2+l2,j3+1)
        xya(l1,l2)=ym1 + (1.d0 + zr)*(y00 - ym1 + zr*(.5d0*ym1 - y00  + .5d0*yp1))
        enddo
        enddo

        yr = ((y-dy)-j2*hgridhold)/hgridhold
        do l1=-1,1
        ym1=xya(l1,-1)
        y00=xya(l1,0)
        yp1=xya(l1,1)
        xa(l1)=ym1 + (1.d0 + yr)*(y00 - ym1 + yr*(.5d0*ym1 - y00  + .5d0*yp1))
        enddo

        xr = ((x-dx)-j1*hgridhold)/hgridhold
        ym1=xa(-1)
        y00=xa(0)
        yp1=xa(1)
        psifscf(i1,i2,i3)=ym1 + (1.d0 + xr)*(y00 - ym1 + xr*(.5d0*ym1 - y00  + .5d0*yp1))

        endif

        enddo
        enddo
        enddo


        call   compress(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
                    mseg_c,mvctr_c,keyg(1,iseg_c),keyv(iseg_c),   & 
                    mseg_f,mvctr_f,keyg(1,iseg_f),keyv(iseg_f),   & 
                    psifscf,psi(ipsi_c),psi(ipsi_f))

         deallocate(psigold,psifscfold,psifscfoex,wwold)
         deallocate(psifscf)

      endif

        close(99)
100  continue


	end



        subroutine writeonewave(iorb,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,hgrid,rxyz,  & 
                           mseg_c,mvctr_c,keyg_c,keyv_c,  & 
                           mseg_f,mvctr_f,keyg_f,keyv_f, & 
                              psi_c,psi_f)
        implicit real*8 (a-h,o-z)
        character*30 filename
        character*4 f4
        dimension keyg_c(2,mseg_c),keyv_c(mseg_c),keyg_f(2,mseg_f),keyv_f(mseg_f)
        dimension psi_c(mvctr_c),psi_f(7,mvctr_f),rxyz(3)


        write(f4,'(i4.4)') iorb
        filename = 'wavefunction.'//f4
        write(*,*) 'opening ',filename
        open(unit=99,file=filename,status='unknown')
         write(99,*) iorb
         write(99,*) hgrid
         write(99,*) n1,n2,n3
         write(99,*) nl1,nu1,nl2,nu2,nl3,nu3
         write(99,'(3(1x,e24.17))') (rxyz(j),j=1,3)
         write(99,*) mvctr_c, mvctr_f

! coarse part
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
            tt=psi_c(i-i0+jj) 
            write(99,'(3(i4),1x,e19.12)') i,i2,i3,tt
          enddo
         enddo
                                                                                                                             
! fine part
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
            t1=psi_f(1,i-i0+jj)
            t2=psi_f(2,i-i0+jj)
            t3=psi_f(3,i-i0+jj)
            t4=psi_f(4,i-i0+jj)
            t5=psi_f(5,i-i0+jj)
            t6=psi_f(6,i-i0+jj)
            t7=psi_f(7,i-i0+jj)
            write(99,'(3(i4),7(1x,e17.10))') i,i2,i3,t1,t2,t3,t4,t5,t6,t7
          enddo
         enddo

          close(99)

	write(*,*) iorb,'th wavefunction written'


	end




	subroutine compresstest(iproc,nproc,n1,n2,n3,norb,nbox_c,nseg,nvctr,keyg,keyv,psi,hpsi)
        implicit real*8 (a-h,o-z)
        parameter(eps_mach=1.d-12)
        dimension nbox_c(2,3,norb),nseg(0:2*norb),nvctr(0:2*norb)
        dimension psi(nvctr(2*norb)),hpsi(nvctr(2*norb))
        dimension keyg(2,nseg(2*norb)),keyv(nseg(2*norb))
        allocatable :: psifscf(:)


! begin test compress-uncompress
       onem=1.d0-eps_mach
       norb_p=int(onem+dble(norb)/dble(nproc))
       do iorb=iproc*norb_p+1,min((iproc+1)*norb_p,norb)
         mseg_c=nseg(2*iorb-1)-nseg(2*iorb-2) 
         mseg_f=nseg(2*iorb  )-nseg(2*iorb-1)
         iseg_c=nseg(2*iorb-2)+1
         iseg_f=nseg(2*iorb-1)+1
         mvctr_c= nvctr(2*iorb-1)-nvctr(2*iorb-2) 
         mvctr_f=(nvctr(2*iorb  )-nvctr(2*iorb-1))/7
         ipsi_c=nvctr(2*iorb-2)+1
         ipsi_f=nvctr(2*iorb-1)+1
         nl1=nbox_c(1,1,iorb) ; nu1=nbox_c(2,1,iorb)
         nl2=nbox_c(1,2,iorb) ; nu2=nbox_c(2,2,iorb)
         nl3=nbox_c(1,3,iorb) ; nu3=nbox_c(2,3,iorb)

        allocate(psifscf((2*(nu1-nl1)+16)*(2*(nu2-nl2)+16)*(2*(nu3-nl3)+16)) )

	do i=1,mvctr_c
        call random_number(tt)
        psi(ipsi_c+i-1)=tt
        hpsi(ipsi_c+i-1)=tt
        enddo

	do i=1,7*mvctr_f
        call random_number(tt)
        psi(ipsi_f+i-1)=tt
        hpsi(ipsi_f+i-1)=tt
        enddo

        call uncompress(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
                    mseg_c,mvctr_c,keyg(1,iseg_c),keyv(iseg_c),   & 
                    mseg_f,mvctr_f,keyg(1,iseg_f),keyv(iseg_f),   & 
                    psi(ipsi_c),psi(ipsi_f),psifscf)
        call   compress(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
                    mseg_c,mvctr_c,keyg(1,iseg_c),keyv(iseg_c),   & 
                    mseg_f,mvctr_f,keyg(1,iseg_f),keyv(iseg_f),   & 
                    psifscf,psi(ipsi_c),psi(ipsi_f))

        tc=0.d0
	do i=1,mvctr_c
        tc=max(tc,abs(psi(ipsi_c+i-1)-hpsi(ipsi_c+i-1)))
        enddo
        if (tc.gt.1.d-10) stop 'coarse compress error'

        tf=0.d0
	do i=1,7*mvctr_f
        tf=max(tf,abs(psi(ipsi_f+i-1)-hpsi(ipsi_f+i-1))) 
        enddo
        write(*,'(a,i4,2(1x,e9.2))') 'COMPRESSION TEST: iorb, coarse, fine error:',iorb,tc,tf
        if (tf.gt.1.d-10) stop 'fine compress error'

        deallocate(psifscf)

     enddo

	return
	end


        subroutine input_rho_ion(iproc,ntypes,nat,iatype,atomnames,rxyz,psppar,nelpsp,n1,n2,n3,hgrid,rho,eion)
! Creates Initial charge density 
        implicit real*8 (a-h,o-z)
        character*20 :: atomnames(100)
        dimension psppar(0:2,0:4,ntypes),rxyz(3,nat),iatype(nat),nelpsp(ntypes)
        dimension rho(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16)

        hgridh=hgrid*.5d0 
        pi=4.d0*atan(1.d0)
	call zero((2*n1+31)*(2*n2+31)*(2*n3+31),rho)

! Ionic charge (Simple implementation, 
!               Should  be calculated on a finer grid and transformed to coarse grid)
   rholeaked=0.d0
   eion=0.d0
   do iat=1,nat
   ityp=iatype(iat)
   rx=rxyz(1,iat) ; ry=rxyz(2,iat) ; rz=rxyz(3,iat)
   ix=nint(rx/hgridh) ; iy=nint(ry/hgridh) ; iz=nint(rz/hgridh)
!    ion-ion interaction
     do jat=1,iat-1
     dist=sqrt( (rx-rxyz(1,jat))**2+(ry-rxyz(2,jat))**2+(rz-rxyz(3,jat))**2 )
     jtyp=iatype(jat)
     eion=eion+nelpsp(jtyp)*nelpsp(ityp)/dist
     enddo

     rloc=psppar(0,0,ityp)
     if (iproc.eq.0) write(*,'(a,i3,a,a,a,e10.3)') 'atom ',iat,' of type ',atomnames(ityp),' has an ionic charge with rloc',rloc
     charge=nelpsp(ityp)/(2.d0*pi*sqrt(2.d0*pi)*rloc**3)
     cutoff=10.d0*rloc
     ii=nint(cutoff/hgridh)

      do i3=iz-ii,iz+ii
      do i2=iy-ii,iy+ii
      do i1=ix-ii,ix+ii
         x=i1*hgridh-rx
         y=i2*hgridh-ry
         z=i3*hgridh-rz
         r2=x**2+y**2+z**2
         arg=r2/rloc**2
        xp=exp(-.5d0*arg)
        if (i3.ge.-14 .and. i3.le.2*n3+16  .and.  & 
            i2.ge.-14 .and. i2.le.2*n2+16  .and.  & 
            i1.ge.-14 .and. i1.le.2*n1+16 ) then
        rho(i1,i2,i3)=rho(i1,i2,i3)-xp*charge
        else
        rholeaked=rholeaked+xp*charge
        endif
      enddo
      enddo
      enddo

    enddo

! Check
	tt=0.d0
        do i3= -14,2*n3+16
        do i2= -14,2*n2+16
        do i1= -14,2*n1+16
        tt=tt+rho(i1,i2,i3)
        enddo
        enddo
        enddo
        tt=tt*hgridh**3
        rholeaked=rholeaked*hgridh**3
	if (iproc.eq.0) write(*,'(a,e21.14,1x,e10.3)') 'total ionic charge,leaked charge: ',tt,rholeaked

        return
	end


!subroutine input_wf_wan(iproc,n1,n2,n3,hgrid,nat,wan_par,ntypes,iatype,psppar, & 
!                        norb,nbox_c,keyg,keyv,nseg,nvctr,psi)
!an wavefunctions in centered in the middle of the bonds
!
!implicit real*8 (a-h,o-z)
!dimension psppar(0:2,0:4,ntypes),iatype(nat)
!dimension nbox_c(2,3,norb),wan_par(0:3,norb)
!dimension nseg(0:2*norb),nvctr(0:2*norb)
!dimension keyg(2,nseg(2*norb)),keyv(nseg(2*norb))
!dimension psi(nvctr(2*norb))
!allocatable :: psifscf(:,:,:)
!
!write(*,*) 'STARTING CALCULATING WANNIER INPUT WAVEFUNCTIONS'
!
!hgridh=.5d0*hgrid
!40
!iorb=1,norb
!iunit+1
!nl1=nbox_c(1,1,iorb) ; nu1=nbox_c(2,1,iorb)
!nl2=nbox_c(1,2,iorb) ; nu2=nbox_c(2,2,iorb)
!nl3=nbox_c(1,3,iorb) ; nu3=nbox_c(2,3,iorb)
!
!allocate(psifscf(-7+2*nl1:2*nu1+8,-7+2*nl2:2*nu2+8,-7+2*nl3:2*nu3+8))
!call zero((2*(nu1-nl1)+16)*(2*(nu2-nl2)+16)*(2*(nu3-nl3)+16),psifscf)
!
!rx=wan_par(1,iorb)
!ry=wan_par(2,iorb)
!rz=wan_par(3,iorb)
!rbond=.5d0*wan_par(0,iorb)
!ii=nint(5.d0*rbond/hgridh)
!
!ix=nint(rx/hgridh) ; iy=nint(ry/hgridh) ; iz=nint(rz/hgridh)
! write(*,*) 'W, center, width',rx,ry,rz,rbond
! write(*,*) 'W, ix,iy,iz',ix,iy,iz
! write(*,*) 'W, boundsx',max(-7+2*nl1,ix-ii),min(2*nu1+8,ix+ii)
! write(*,*) 'W, boundsy',max(-7+2*nl2,iy-ii),min(2*nu2+8,iy+ii)
! write(*,*) 'W, boundsz',max(-7+2*nl3,iz-ii),min(2*nu3+8,iz+ii)
!
!t=0.d0
! do i3=max(-7+2*nl3,iz-ii),min(2*nu3+8,iz+ii)
! do i2=max(-7+2*nl2,iy-ii),min(2*nu2+8,iy+ii)
! do i1=max(-7+2*nl1,ix-ii),min(2*nu1+8,ix+ii)
!    x=i1*hgridh-rx
!    y=i2*hgridh-ry
!    z=i3*hgridh-rz
!    r2=x**2+y**2+z**2
!    arg=r2/rbond**2
!    xp=arg*exp(-.5d0*arg)
!    psifscf(i1,i2,i3)=xp
!+xp**2
! enddo ; enddo ; enddo
!
!(*,*) 'W, norm',tt
!
!g
!agonal
!do i=max(-7+2*nl1,-7+2*nl2,-7+2*nl3),min(2*nu1+8,2*nu2+8,2*nu3+8)
!    i1=i ; i2=i ; i3=i
!    write(iunit,'(3(1x,e10.3),1x,e12.5)') i1*hgridh,i2*hgridh,i3*hgridh,psifscf(i1,i2,i3)
!enddo
!
!nc1=2*(nu1-nl1)+1
!nc2=2*(nu2-nl2)+1
!nc3=2*(nu3-nl3)+1
!diagonal
!do i=max(-7+2*nl1,-7+2*nl2,-7+2*nl3),min(2*nu1+8,2*nu2+8,2*nu3+8)
!    i1=i ; i2=nc2-i ; i3=nc3-i
!    write(iunit,'(3(1x,e10.3),1x,e12.5))') i1*hgridh,i2*hgridh,i3*hgridh,psifscf(i1,i2,i3)
!enddo
!
!diagonal
!do i=max(-7+2*nl1,-7+2*nl2,-7+2*nl3),min(2*nu1+8,2*nu2+8,2*nu3+8)
!    i1=nc1-i ; i2=i ; i3=nc3-i
!    write(iunit,'(3(1x,e10.3),1x,e12.5)') i1*hgridh,i2*hgridh,i3*hgridh,psifscf(i1,i2,i3)
!enddo
!
!diagonal
!do i=max(-7+2*nl1,-7+2*nl2,-7+2*nl3),min(2*nu1+8,2*nu2+8,2*nu3+8)
!    i1=nc1-i ; i2=nc2-i ; i3=i
!    write(iunit,'(3(1x,e10.3),1x,e12.5)') i1*hgridh,i2*hgridh,i3*hgridh,psifscf(i1,i2,i3)
!enddo
!sting
!
!
! mseg_c=nseg(2*iorb-1)-nseg(2*iorb-2) 
! mseg_f=nseg(2*iorb  )-nseg(2*iorb-1)
! iseg_c=nseg(2*iorb-2)+1
! iseg_f=nseg(2*iorb-1)+1
! mvctr_c= nvctr(2*iorb-1)-nvctr(2*iorb-2) 
! mvctr_f=(nvctr(2*iorb  )-nvctr(2*iorb-1))/7
! ipsi_c=nvctr(2*iorb-2)+1
! ipsi_f=nvctr(2*iorb-1)+1
!
! call   compress(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
!            mseg_c,mvctr_c,keyg(1,iseg_c),keyv(iseg_c),   & 
!            mseg_f,mvctr_f,keyg(1,iseg_f),keyv(iseg_f),   & 
!            psifscf,psi(ipsi_c),psi(ipsi_f))
!
!deallocate(psifscf)
!do
!
!g
!call loewe(iproc,norb,nvctr,nseg,keyg,keyv,psi)
! toler=1.d-6
!call checkortho(iproc,toler,norb,nvctr,nseg,keyg,keyv,psi)
!
!write(*,*) 'FINISHED CALCULATING INPUT WAVEFUNCTIONS'
!
!return
!end


	subroutine checkortho(iproc,toler,norb,nvctr,nseg,keyg,keyv,psi)
        implicit real*8 (a-h,o-z)
        dimension nseg(0:2*norb),nvctr(0:2*norb)
        dimension keyg(2,nseg(2*norb)),keyv(nseg(2*norb))
        dimension psi(nvctr(2*norb))

        write(*,*) 'Start checking orthogonality',iproc
        dev=0.d0
     do 100,iorb=1,norb
       maseg_c=nseg(2*iorb-1)-nseg(2*iorb-2)
       maseg_f=nseg(2*iorb  )-nseg(2*iorb-1)
       iseg_c=nseg(2*iorb-2)+1
       iseg_f=nseg(2*iorb-1)+1
       mavctr_c= nvctr(2*iorb-1)-nvctr(2*iorb-2)
       mavctr_f=(nvctr(2*iorb  )-nvctr(2*iorb-1))/7
       ipsi_c=nvctr(2*iorb-2)+1
       ipsi_f=nvctr(2*iorb-1)+1

         do 100,jorb=1,norb
           mbseg_c=nseg(2*jorb-1)-nseg(2*jorb-2)
           mbseg_f=nseg(2*jorb  )-nseg(2*jorb-1)
           jseg_c=nseg(2*jorb-2)+1
           jseg_f=nseg(2*jorb-1)+1
           mbvctr_c= nvctr(2*jorb-1)-nvctr(2*jorb-2)
           mbvctr_f=(nvctr(2*jorb  )-nvctr(2*jorb-1))/7
           jpsi_c=nvctr(2*jorb-2)+1
           jpsi_f=nvctr(2*jorb-1)+1

	call wdot(  & 
                   mavctr_c,mavctr_f,maseg_c,maseg_f,keyv(iseg_c),keyv(iseg_f),  & 
                   keyg(1,iseg_c),keyg(1,iseg_f),psi(ipsi_c),psi(ipsi_f),  & 
                   mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv(jseg_c),keyv(jseg_f),  & 
                   keyg(1,jseg_c),keyg(1,jseg_f),psi(jpsi_c),psi(jpsi_f),scpr)
        if (iorb.eq.jorb) then
        dev=dev+(scpr-1.d0)**2
        else
        dev=dev+scpr**2
        endif
!        if (iorb.eq.jorb .and. abs(scpr-1.d0).gt.toler)  write(*,*) 'ERROR ORTHO',iorb,jorb,scpr
!        if (iorb.ne.jorb .and. abs(scpr).gt.toler)  write(*,*) 'ERROR ORTHO',iorb,jorb,scpr
100    continue
        write(*,*) 'Deviation from orthogonality ',iproc,dev

        return
	end


	subroutine loewe(iproc,norb,nvctr,nseg,keyg,keyv,psi)
! loewdin orthogonalisation
        implicit real*8 (a-h,o-z)
        dimension nseg(0:2*norb),nvctr(0:2*norb)
        dimension keyg(2,nseg(2*norb)),keyv(nseg(2*norb))
        dimension psi(nvctr(2*norb))
        real*8, allocatable :: ovrlp(:,:,:),evall(:),tpsi(:)

 if (norb.eq.1) then
      iorb=1
      mvctr_c= nvctr(2*iorb-1)-nvctr(2*iorb-2)
      mvctr_f=(nvctr(2*iorb  )-nvctr(2*iorb-1))/7
      ipsi_c=nvctr(2*iorb-2)+1
      ipsi_f=nvctr(2*iorb-1)+1
      call  wnrm(mvctr_c,mvctr_f,psi(ipsi_c),psi(ipsi_f),scpr) ; scpr=1.d0/sqrt(scpr)
      call wscal(mvctr_c,mvctr_f,scpr,psi(ipsi_c),psi(ipsi_f)) 
 else

        allocate(ovrlp(norb,norb,3),evall(norb),tpsi(nvctr(2*norb)))

        offdiag=0.d0
     do 100,iorb=1,norb
       maseg_c=nseg(2*iorb-1)-nseg(2*iorb-2)
       maseg_f=nseg(2*iorb  )-nseg(2*iorb-1)
       iseg_c=nseg(2*iorb-2)+1
       iseg_f=nseg(2*iorb-1)+1
       mavctr_c= nvctr(2*iorb-1)-nvctr(2*iorb-2)
       mavctr_f=(nvctr(2*iorb  )-nvctr(2*iorb-1))/7
       ipsi_c=nvctr(2*iorb-2)+1
       ipsi_f=nvctr(2*iorb-1)+1
       call wcopy(mavctr_c,mavctr_f,psi(ipsi_c),psi(ipsi_f),tpsi(ipsi_c),tpsi(ipsi_f))

! Full matrix for testing
!         do 100,jorb=1,norb
! Lower triangle
         do 100,jorb=1,iorb
           mbseg_c=nseg(2*jorb-1)-nseg(2*jorb-2)
           mbseg_f=nseg(2*jorb  )-nseg(2*jorb-1)
           jseg_c=nseg(2*jorb-2)+1
           jseg_f=nseg(2*jorb-1)+1
           mbvctr_c= nvctr(2*jorb-1)-nvctr(2*jorb-2)
           mbvctr_f=(nvctr(2*jorb  )-nvctr(2*jorb-1))/7
           jpsi_c=nvctr(2*jorb-2)+1
           jpsi_f=nvctr(2*jorb-1)+1

	call wdot(  & 
                   mavctr_c,mavctr_f,maseg_c,maseg_f,keyv(iseg_c),keyv(iseg_f),  & 
                   keyg(1,iseg_c),keyg(1,iseg_f),psi(ipsi_c),psi(ipsi_f),  & 
                   mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv(jseg_c),keyv(jseg_f),  & 
                   keyg(1,jseg_c),keyg(1,jseg_f),psi(jpsi_c),psi(jpsi_f),scpr)

        if (iorb.ne.jorb) offdiag=offdiag+abs(scpr)
       ovrlp(iorb,jorb,1)=scpr
100    continue
       if (iproc.eq.0) write(*,*) 'offdiagonal overlap',offdiag
      
! testing
!       tt=0.d0
!       do jorb=1,norb
!       do iorb=1,jorb-1
!       tt=max(tt,abs(ovrlp(jorb,iorb,1)-ovrlp(iorb,jorb,1)))
!!       write(*,'(i3,i3,1x,e21.14,1x,e21.14,1x,e9.2)') jorb,iorb,ovrlp(jorb,iorb,1),ovrlp(iorb,jorb,1),ovrlp(jorb,iorb,1)-ovrlp(iorb,jorb,1)
!       enddo ; enddo
!       write(*,*) 'iproc,overlap: max violation of symmetry:',iproc,tt

!       if (iproc.eq.0) then 
!        write(6,*) 'loewe S'
!        do i=1,norb
!!        write(6,'(14(1x,e6.1))') (abs(ovrlp(i,j,1)),j=1,min(14,norb))
!        write(6,'(i2,i3,14(1x,e21.14))') iproc,i,(abs(ovrlp(i,j,1)),j=1,min(14,norb))
!        enddo
!       endif

! LAPACK
        call DSYEV('V','L',norb,ovrlp(1,1,1),norb,evall,ovrlp(1,1,3),norb**2,info)
        if (info.ne.0) write(6,*) 'info loewe', info
!        if (iproc.eq.0) then 
!          write(6,*) 'overlap eigenvalues'
!77        format(8(1x,e10.3))
!          if (norb.le.16) then
!          write(6,77) evall
!          else
!          write(6,77) (evall(i),i=1,4), (evall(i),i=norb-3,norb)
!          endif
!        endif

! calculate S^{-1/2} ovrlp(*,*,3)
        do 2935,lorb=1,norb
        do 2935,jorb=1,norb
2935    ovrlp(jorb,lorb,2)=ovrlp(jorb,lorb,1)*sqrt(1.d0/evall(lorb))
!        do 3985,j=1,norb
!        do 3985,i=1,norb
!        ovrlp(i,j,3)=0.d0
!        do 3985,l=1,norb
!3985    ovrlp(i,j,3)=ovrlp(i,j,3)+ovrlp(i,l,1)*ovrlp(j,l,2)
! BLAS:
        call DGEMM('N','T',norb,norb,norb,1.d0,ovrlp(1,1,1),norb,ovrlp(1,1,2),norb,0.d0,ovrlp(1,1,3),norb)

! new eigenvectors
     do 200,iorb=1,norb
       maseg_c=nseg(2*iorb-1)-nseg(2*iorb-2)
       maseg_f=nseg(2*iorb  )-nseg(2*iorb-1)
       iseg_c=nseg(2*iorb-2)+1
       iseg_f=nseg(2*iorb-1)+1
       mavctr_c= nvctr(2*iorb-1)-nvctr(2*iorb-2)
       mavctr_f=(nvctr(2*iorb  )-nvctr(2*iorb-1))/7
       ipsi_c=nvctr(2*iorb-2)+1
       ipsi_f=nvctr(2*iorb-1)+1
       call wzero(mavctr_c,mavctr_f,psi(ipsi_c),psi(ipsi_f))

         do 200,jorb=1,norb
           mbseg_c=nseg(2*jorb-1)-nseg(2*jorb-2)
           mbseg_f=nseg(2*jorb  )-nseg(2*jorb-1)
           jseg_c=nseg(2*jorb-2)+1
           jseg_f=nseg(2*jorb-1)+1
           mbvctr_c= nvctr(2*jorb-1)-nvctr(2*jorb-2)
           mbvctr_f=(nvctr(2*jorb  )-nvctr(2*jorb-1))/7
           jpsi_c=nvctr(2*jorb-2)+1
           jpsi_f=nvctr(2*jorb-1)+1

        scpr=ovrlp(jorb,iorb,3)
	call waxpy(  & 
                   scpr,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv(jseg_c),keyv(jseg_f),  & 
                   keyg(1,jseg_c),keyg(1,jseg_f),tpsi(jpsi_c),tpsi(jpsi_f), &
                   mavctr_c,mavctr_f,maseg_c,maseg_f,keyv(iseg_c),keyv(iseg_f),  & 
                   keyg(1,iseg_c),keyg(1,iseg_f),psi(ipsi_c),psi(ipsi_f))
200     continue

        deallocate(ovrlp,evall,tpsi)

 endif

        return
        end




        subroutine crtparabolicpot(n1,n2,n3,hgrid,pot)
!  parabolic potential centered in the middle of the box
        implicit real*8 (a-h,o-z)
        dimension pot(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16)

        hgridh=hgrid*.5d0 

        rx=n1*hgridh ; ry=n2*hgridh ; rz=n3*hgridh
      do i3=-14,2*n3+16
         z=i3*hgridh-rz
      do i2=-14,2*n2+16
         y=i2*hgridh-ry
      do i1=-14,2*n1+16
         x=i1*hgridh-rx
         r2=x**2+y**2+z**2
        pot(i1,i2,i3)=.5d0*r2
      enddo
      enddo
      enddo

        return
	end


        subroutine crtplanepot(n1,n2,n3,hgrid,pot)
!  linear potential centered in the middle of the box
        implicit real*8 (a-h,o-z)
        dimension pot(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16)

        hgridh=hgrid*.5d0 

         rz=n3*hgridh
      do i3=-14,2*n3+16
         z=i3*hgridh-rz
      do i2=-14,2*n2+16
      do i1=-14,2*n1+16
        pot(i1,i2,i3)=z
        pot(i1,i2,i3)=1.d0
      enddo
      enddo
      enddo

        return
	end

        subroutine addlocgauspsp(iproc,ntypes,nat,iatype,atomnames,rxyz,psppar,n1,n2,n3,hgrid,pot)
! Add local Gaussian terms of the PSP to pot 
        implicit real*8 (a-h,o-z)
        character*20 :: atomnames(100)
        dimension psppar(0:2,0:4,ntypes),rxyz(3,nat),iatype(nat)
        dimension pot(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16)

        hgridh=hgrid*.5d0 

   do iat=1,nat
   ityp=iatype(iat)
   rx=rxyz(1,iat) ; ry=rxyz(2,iat) ; rz=rxyz(3,iat)
   ix=nint(rx/hgridh) ; iy=nint(ry/hgridh) ; iz=nint(rz/hgridh)
! determine number of local terms
     nloc=0
     do iloc=1,4
     if (psppar(0,iloc,ityp).ne.0.d0) nloc=iloc
     enddo

     rloc=psppar(0,0,ityp)
     if (iproc.eq.0) write(*,'(a,i4,a,a,a,i3,a,1x,e9.2)')  & 
     'atom ',iat,' is of type ',atomnames(ityp),' and has ',nloc,' local terms with rloc',rloc
     cutoff=10.d0*rloc
     ii=nint(cutoff/hgridh)

      do i3=max(-14,iz-ii),min(2*n3+16,iz+ii)
      do i2=max(-14,iy-ii),min(2*n2+16,iy+ii)
      do i1=max(-14,ix-ii),min(2*n1+16,ix+ii)
         x=i1*hgridh-rx
         y=i2*hgridh-ry
         z=i3*hgridh-rz
         r2=x**2+y**2+z**2
         arg=r2/rloc**2
        xp=exp(-.5d0*arg)
        tt=psppar(0,nloc,ityp)
	do iloc=nloc-1,1,-1
        tt=arg*tt+psppar(0,iloc,ityp)
        enddo
        pot(i1,i2,i3)=pot(i1,i2,i3)+xp*tt
      enddo
      enddo
      enddo

!! For testing only: Add erf part (in the final version that should be part of Hartree pot)
!      do i3=-14,2*n3+16
!      do i2=-14,2*n2+16
!      do i1=-14,2*n1+16
!         x=(i1-ix)*hgridh
!         y=(i2-iy)*hgridh
!         z=(i3-iz)*hgridh
!         r2=x**2+y**2+z**2
!         r=sqrt(r2)
!         arg=r*(sqrt(.5d0)/rloc)
!         if (arg.lt.1.d-7) then 
!! Taylor expansion
!         x=arg**2
!         tt=   -0.37612638903183752463d0*x + 1.1283791670955125739d0
!         tt=tt*(sqrt(.5d0)/rloc)
!         else
!          tt=derf(arg)/r
!         endif
!        pot(i1,i2,i3)=pot(i1,i2,i3)+nelpsp(ityp)*tt
!      enddo
!      enddo
!      enddo


    enddo

        return
	end


        subroutine chargedens(parallel,iproc,nproc,norb,n1,n2,n3,nbox_c,hgrid,occup,  & 
                              nseg,nvctr,keyg,keyv,psi,rho)
! Calculates the charge density 
! Input: psi
! Output: rho
        implicit real*8 (a-h,o-z)
        logical parallel
        integer count1,count2,count_rate,count_max
        parameter(eps_mach=1.d-12)
	dimension rho((2*n1+31)*(2*n2+31)*(2*n3+31)),occup(norb)
        dimension nbox_c(2,3,norb)
        dimension nseg(0:2*norb),nvctr(0:2*norb)
        dimension keyg(2,nseg(2*norb)),keyv(nseg(2*norb))
        dimension psi(nvctr(2*norb))
        real*8, allocatable, dimension(:) :: psifscf,psir,rho_p
        include 'mpif.h'


        onem=1.d0-eps_mach
        norb_p=int(onem+dble(norb)/dble(nproc))

        hgridh=hgrid*.5d0 

! Determine aximal size of work arrays
      nl1=10000 ; nu1=0
      nl2=10000 ; nu2=0
      nl3=10000 ; nu3=0
      do iorb=iproc*norb_p+1,min((iproc+1)*norb_p,norb)
        nl1=min(nl1,nbox_c(1,1,iorb)) ; nu1=max(nu1,nbox_c(2,1,iorb))
        nl2=min(nl2,nbox_c(1,2,iorb)) ; nu2=max(nu2,nbox_c(2,2,iorb))
        nl3=min(nl3,nbox_c(1,3,iorb)) ; nu3=max(nu3,nbox_c(2,3,iorb))
      enddo
! Wavefunction expressed everywhere in fine scaling functions (for potential and kinetic energy)
        allocate(psifscf((2*(nu1-nl1)+16)*(2*(nu2-nl2)+16)*(2*(nu3-nl3)+16)) )
! Wavefunction in real space
        allocate(psir((2*(nu1-nl1)+31)*(2*(nu2-nl2)+31)*(2*(nu3-nl3)+31)))

 if (parallel) then
        allocate(rho_p((2*n1+31)*(2*n2+31)*(2*n3+31)))
	call zero((2*n1+31)*(2*n2+31)*(2*n3+31),rho_p)

      do iorb=iproc*norb_p+1,min((iproc+1)*norb_p,norb)
         mseg_c=nseg(2*iorb-1)-nseg(2*iorb-2) 
         mseg_f=nseg(2*iorb  )-nseg(2*iorb-1)
         iseg_c=nseg(2*iorb-2)+1
         iseg_f=nseg(2*iorb-1)+1
         mvctr_c= nvctr(2*iorb-1)-nvctr(2*iorb-2) 
         mvctr_f=(nvctr(2*iorb  )-nvctr(2*iorb-1))/7
         ipsi_c=nvctr(2*iorb-2)+1
         ipsi_f=nvctr(2*iorb-1)+1
         nl1=nbox_c(1,1,iorb) ; nu1=nbox_c(2,1,iorb)
         nl2=nbox_c(1,2,iorb) ; nu2=nbox_c(2,2,iorb)
         nl3=nbox_c(1,3,iorb) ; nu3=nbox_c(2,3,iorb)

        call uncompress(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
                    mseg_c,mvctr_c,keyg(1,iseg_c),keyv(iseg_c),   & 
                    mseg_f,mvctr_f,keyg(1,iseg_f),keyv(iseg_f),   & 
                    psi(ipsi_c),psi(ipsi_f),psifscf)

        call convolut_magic_n(2*(nu1-nl1)+15,2*(nu2-nl2)+15,2*(nu3-nl3)+15,psifscf,psir) 

        const=occup(iorb)/hgridh**3
	call addpartrho(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,const,psir,rho_p)

      enddo

        call cpu_time(tr0)
        call system_clock(count1,count_rate,count_max)
        call MPI_ALLREDUCE(rho_p,rho,(2*n1+31)*(2*n2+31)*(2*n3+31),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call cpu_time(tr1)
        call system_clock(count2,count_rate,count_max)
        tel=dble(count2-count1)/dble(count_rate)
        write(78,*) 'RHO: ALLREDUCE TIME',iproc,tr1-tr0,tel
        write(78,*) '---------------------------------------------'

        deallocate(rho_p)
 else

	call zero((2*n1+31)*(2*n2+31)*(2*n3+31),rho)

     do iorb=1,norb
         mseg_c=nseg(2*iorb-1)-nseg(2*iorb-2) 
         mseg_f=nseg(2*iorb  )-nseg(2*iorb-1)
         iseg_c=nseg(2*iorb-2)+1
         iseg_f=nseg(2*iorb-1)+1
         mvctr_c= nvctr(2*iorb-1)-nvctr(2*iorb-2) 
         mvctr_f=(nvctr(2*iorb  )-nvctr(2*iorb-1))/7
         ipsi_c=nvctr(2*iorb-2)+1
         ipsi_f=nvctr(2*iorb-1)+1
         nl1=nbox_c(1,1,iorb) ; nu1=nbox_c(2,1,iorb)
         nl2=nbox_c(1,2,iorb) ; nu2=nbox_c(2,2,iorb)
         nl3=nbox_c(1,3,iorb) ; nu3=nbox_c(2,3,iorb)

        call uncompress(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
                    mseg_c,mvctr_c,keyg(1,iseg_c),keyv(iseg_c),   & 
                    mseg_f,mvctr_f,keyg(1,iseg_f),keyv(iseg_f),   & 
                    psi(ipsi_c),psi(ipsi_f),psifscf)

        call convolut_magic_n(2*(nu1-nl1)+15,2*(nu2-nl2)+15,2*(nu3-nl3)+15,psifscf,psir) 

        const=occup(iorb)/hgridh**3
	call addpartrho(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,const,psir,rho)

     enddo
 endif

! Check
        tt=0.d0
        do i=1,(2*n1+31)*(2*n2+31)*(2*n3+31)
         tt=tt+rho(i)
        enddo
        tt=tt*hgridh**3
	if (iproc.eq.0) write(*,*) 'Total charge from routine chargedens',tt,iproc


        deallocate(psifscf,psir)

        return
	end


	subroutine addpartrho(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,const,psir,rho)
        implicit real*8 (a-h,o-z)
        dimension rho(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16)
        dimension psir(-14+2*nl1:2*nu1+16,-14+2*nl2:2*nu2+16,-14+2*nl3:2*nu3+16)

        do i3=-14+2*nl3,2*nu3+16
        do i2=-14+2*nl2,2*nu2+16
        do i1=-14+2*nl1,2*nu1+16
         rho(i1,i2,i3)=rho(i1,i2,i3)+const*psir(i1,i2,i3)**2
        enddo ; enddo ; enddo

        return
        end


        subroutine applylocpotkin(iproc,nproc,norb,n1,n2,n3,nbox_c,hgrid,occup,  & 
                   nseg,nvctr,keyg,keyv,psi,pot,hpsi,epot_sum,ekin_sum)
!  Applies the local potential and kinetic energy operator to a wavefunction
! Input: pot,psi
! Output: hpsi,epot,ekin
        implicit real*8 (a-h,o-z)
        parameter(eps_mach=1.d-12)
        dimension nbox_c(2,3,norb),occup(norb),pot((2*n1+31)*(2*n2+31)*(2*n3+31))
        dimension nseg(0:2*norb),nvctr(0:2*norb)
        dimension keyg(2,nseg(2*norb)),keyv(nseg(2*norb))
        dimension psi(nvctr(2*norb)),hpsi(nvctr(2*norb))
        real*8, allocatable, dimension(:) :: psifscf,psifscfk,psir
        include 'mpif.h'

        onem=1.d0-eps_mach
        norb_p=int(onem+dble(norb)/dble(nproc))

        hgridh=hgrid*.5d0

! Determine aximal size of work arrays
      nl1=10000 ; nu1=0
      nl2=10000 ; nu2=0
      nl3=10000 ; nu3=0
      do iorb=iproc*norb_p+1,min((iproc+1)*norb_p,norb)
        nl1=min(nl1,nbox_c(1,1,iorb)) ; nu1=max(nu1,nbox_c(2,1,iorb))
        nl2=min(nl2,nbox_c(1,2,iorb)) ; nu2=max(nu2,nbox_c(2,2,iorb))
        nl3=min(nl3,nbox_c(1,3,iorb)) ; nu3=max(nu3,nbox_c(2,3,iorb))
      enddo
! Wavefunction expressed everywhere in fine scaling functions (for potential and kinetic energy)
        allocate(psifscf((2*(nu1-nl1)+16)*(2*(nu2-nl2)+16)*(2*(nu3-nl3)+16)) )
        allocate(psifscfk((2*(nu1-nl1)+16)*(2*(nu2-nl2)+16)*(2*(nu3-nl3)+16)) )
! Wavefunction in real space
        allocate(psir((2*(nu1-nl1)+31)*(2*(nu2-nl2)+31)*(2*(nu3-nl3)+31)))


        ekin_sum=0.d0
        epot_sum=0.d0
     do iorb=iproc*norb_p+1,min((iproc+1)*norb_p,norb)
         mseg_c=nseg(2*iorb-1)-nseg(2*iorb-2) 
         mseg_f=nseg(2*iorb  )-nseg(2*iorb-1)
         iseg_c=nseg(2*iorb-2)+1
         iseg_f=nseg(2*iorb-1)+1
         mvctr_c= nvctr(2*iorb-1)-nvctr(2*iorb-2) 
         mvctr_f=(nvctr(2*iorb  )-nvctr(2*iorb-1))/7
         ipsi_c=nvctr(2*iorb-2)+1
         ipsi_f=nvctr(2*iorb-1)+1
         nl1=nbox_c(1,1,iorb) ; nu1=nbox_c(2,1,iorb)
         nl2=nbox_c(1,2,iorb) ; nu2=nbox_c(2,2,iorb)
         nl3=nbox_c(1,3,iorb) ; nu3=nbox_c(2,3,iorb)

        call uncompress(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  &
                    mseg_c,mvctr_c,keyg(1,iseg_c),keyv(iseg_c),   &
                    mseg_f,mvctr_f,keyg(1,iseg_f),keyv(iseg_f),   &
                    psi(ipsi_c),psi(ipsi_f),psifscf)

        call convolut_kinetic(2*(nu1-nl1)+15,2*(nu2-nl2)+15,2*(nu3-nl3)+15,hgridh,psifscf,psifscfk)

        ekin=0.d0 
        do i=1,(2*(nu1-nl1)+16)*(2*(nu2-nl2)+16)*(2*(nu3-nl3)+16)
          ekin=ekin+psifscf(i)*psifscfk(i)
        enddo
        ekin_sum=ekin_sum+occup(iorb)*ekin 

        call convolut_magic_n(2*(nu1-nl1)+15,2*(nu2-nl2)+15,2*(nu3-nl3)+15,psifscf,psir) 

        call  potloc(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,pot,psir,epot)
        epot_sum=epot_sum+occup(iorb)*epot
        write(*,'(a,i5,2(1x,e17.10))') 'iorb,ekin,epot',iorb,ekin,epot

        call convolut_magic_t(2*(nu1-nl1)+15,2*(nu2-nl2)+15,2*(nu3-nl3)+15,psir,psifscf)

        do i=1,(2*(nu1-nl1)+16)*(2*(nu2-nl2)+16)*(2*(nu3-nl3)+16)
          psifscf(i)=psifscf(i)+psifscfk(i)
        enddo

        call   compress(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
                    mseg_c,mvctr_c,keyg(1,iseg_c),keyv(iseg_c),   & 
                    mseg_f,mvctr_f,keyg(1,iseg_f),keyv(iseg_f),   & 
                    psifscf,hpsi(ipsi_c),hpsi(ipsi_f))
     enddo

        deallocate(psifscf,psifscfk,psir)
   
        return
	end


        subroutine potloc(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,pot,psir,epot)
        implicit real*8 (a-h,o-z)
        dimension pot(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16)
        dimension psir(-14+2*nl1:2*nu1+16,-14+2*nl2:2*nu2+16,-14+2*nl3:2*nu3+16)

        epot=0.d0 
        do i3=-14+2*nl3,2*nu3+16
        do i2=-14+2*nl2,2*nu2+16
        do i1=-14+2*nl1,2*nu1+16
          tt=pot(i1,i2,i3)*psir(i1,i2,i3)
          epot=epot+tt*psir(i1,i2,i3)
          psir(i1,i2,i3)=tt
        enddo ; enddo ; enddo

!        do i3=-14+2*nl3,2*nu3+16
!        i1=i3 ; i2=i3 
!        write(77,*) i3,pot(i1,i2,i3),psir(i1,i2,i3)
!        enddo

	return
        end


	

        subroutine numb_proj(ityp,ntypes,psppar,mproj)
! Determines the number of projectors
        implicit real*8 (a-h,o-z)
        dimension psppar(0:2,0:4,ntypes)

        mproj=0
! ONLY GTH PSP form (not HGH)
          do i=1,2
          do j=1,2
            if (psppar(i,j,ityp).ne.0.d0) mproj=mproj+2*i-1
          enddo
          enddo

        return
        end


        subroutine crtinpwave(iproc,nproc,n1,n2,n3,nbox_c,nbox_f,norb,  & 
                              wan_par,nseg,nvctr,keyg,keyv,hgrid,psi)

        implicit real*8 (a-h,o-z)
        parameter(eps_mach=1.d-12)
        dimension nseg(0:2*norb),nvctr(0:2*norb)
        dimension keyg(2,nseg(2*norb)),keyv(nseg(2*norb))
        dimension psi(nvctr(2*norb))
        dimension nbox_c(2,3,norb),nbox_f(2,3,norb),wan_par(0:3,norb)

       onem=1.d0-eps_mach
       norb_p=int(onem+dble(norb)/dble(nproc))
    do iorb=iproc*norb_p+1,min((iproc+1)*norb_p,norb)

         rx=wan_par(1,iorb)
         ry=wan_par(2,iorb)
         rz=wan_par(3,iorb)
         gau_a=wan_par(0,iorb)

         mseg_c=nseg(2*iorb-1)-nseg(2*iorb-2)
         mseg_f=nseg(2*iorb  )-nseg(2*iorb-1)
         iseg_c=nseg(2*iorb-2)+1
         iseg_f=nseg(2*iorb-1)+1
         mvctr_c= nvctr(2*iorb-1)-nvctr(2*iorb-2)
         mvctr_f=(nvctr(2*iorb  )-nvctr(2*iorb-1))/7
         ipsi_c=nvctr(2*iorb-2)+1
         ipsi_f=nvctr(2*iorb-1)+1
         nl1_c=nbox_c(1,1,iorb) ; nu1_c=nbox_c(2,1,iorb)
         nl2_c=nbox_c(1,2,iorb) ; nu2_c=nbox_c(2,2,iorb)
         nl3_c=nbox_c(1,3,iorb) ; nu3_c=nbox_c(2,3,iorb)
         nl1_f=nbox_f(1,1,iorb) ; nu1_f=nbox_f(2,1,iorb)
         nl2_f=nbox_f(1,2,iorb) ; nu2_f=nbox_f(2,2,iorb)
         nl3_f=nbox_f(1,3,iorb) ; nu3_f=nbox_f(2,3,iorb)

          call crtonewave(n1,n2,n3, & 
          nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,  & 
          mseg_c,mvctr_c,keyg(1,iseg_c),keyv(iseg_c),mseg_f,mvctr_f,keyg(1,iseg_f),keyv(iseg_f),  &
          hgrid,gau_a,rx,ry,rz,psi(ipsi_c),psi(ipsi_f))
     enddo

        return
        end


        subroutine crtonewave(n1,n2,n3, & 
                   nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,  & 
                   mseg_c,mvctr_c,keyg_c,keyv_c,mseg_f,mvctr_f,keyg_f,keyv_f,  &
                   hgrid,gau_a,rx,ry,rz,psi_c,psi_f)
! returns the compressed form of a Gaussian 
! exp (-1/(2*gau_a^2) *((x-cntrx)^2 + (y-cntry)^2 + (z-cntrz)^2 ))
! in the arrays psi_c, psi_f
        implicit real*8 (a-h,o-z)
        parameter(nw=16000)
        dimension keyg_c(2,mseg_c),keyv_c(mseg_c),keyg_f(2,mseg_f),keyv_f(mseg_f)
        dimension psi_c(mvctr_c),psi_f(7,mvctr_f)
        real*8, allocatable, dimension(:,:) :: wprojx, wprojy, wprojz
        real*8, allocatable, dimension(:,:) :: work
        real*8, allocatable :: psig_c(:,:,:), psig_f(:,:,:,:)

        allocate(wprojx(0:n1,2),wprojy(0:n2,2),wprojz(0:n3,2),work(0:nw,2))
        allocate(psig_c(nl1_c:nu1_c,nl2_c:nu2_c,nl3_c:nu3_c))
        allocate(psig_f(7,nl1_f:nu1_f,nl2_f:nu2_f,nl3_f:nu3_f))

        n_gau=0
        CALL GAUSS_TO_DAUB(hgrid,1.d0,rx,gau_a,n_gau,n1,ml1,mu1,wprojx(0,1),te,work,nw)
        CALL GAUSS_TO_DAUB(hgrid,1.d0,ry,gau_a,n_gau,n2,ml2,mu2,wprojy(0,1),te,work,nw)
        CALL GAUSS_TO_DAUB(hgrid,1.d0,rz,gau_a,n_gau,n3,ml3,mu3,wprojz(0,1),te,work,nw)
       if (ml1.gt.min(nl1_c,nl1_f)) write(*,*) 'ERRORi: ml1',ml1
       if (ml2.gt.min(nl2_c,nl2_f)) write(*,*) 'ERRORi: ml2',ml2
       if (ml3.gt.min(nl3_c,nl3_f)) write(*,*) 'ERRORi: ml3',ml3
       if (mu1.lt.max(nu1_c,nu1_f)) write(*,*) 'ERRORi: mu1',mu1
       if (mu2.lt.max(nu2_c,nu2_f)) write(*,*) 'ERRORi: mu2',mu2
       if (mu3.lt.max(nu3_c,nu3_f)) write(*,*) 'ERRORi: mu3',mu3

! First term: coarse projector components
          do i3=nl3_c,nu3_c
          do i2=nl2_c,nu2_c
          do i1=nl1_c,nu1_c
            psig_c(i1,i2,i3)=wprojx(i1,1)*wprojy(i2,1)*wprojz(i3,1)
          enddo ; enddo ; enddo

! First term: fine projector components
          do i3=nl3_f,nu3_f
          do i2=nl2_f,nu2_f
          do i1=nl1_f,nu1_f
            psig_f(1,i1,i2,i3)=wprojx(i1,2)*wprojy(i2,1)*wprojz(i3,1)
            psig_f(2,i1,i2,i3)=wprojx(i1,1)*wprojy(i2,2)*wprojz(i3,1)
            psig_f(3,i1,i2,i3)=wprojx(i1,2)*wprojy(i2,2)*wprojz(i3,1)
            psig_f(4,i1,i2,i3)=wprojx(i1,1)*wprojy(i2,1)*wprojz(i3,2)
            psig_f(5,i1,i2,i3)=wprojx(i1,2)*wprojy(i2,1)*wprojz(i3,2)
            psig_f(6,i1,i2,i3)=wprojx(i1,1)*wprojy(i2,2)*wprojz(i3,2)
            psig_f(7,i1,i2,i3)=wprojx(i1,2)*wprojy(i2,2)*wprojz(i3,2)
          enddo ; enddo ; enddo

! coarse part
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
            psi_c(i-i0+jj)=psig_c(i,i2,i3)
          enddo
        enddo

! fine part
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
            psi_f(1,i-i0+jj)=psig_f(1,i,i2,i3)
            psi_f(2,i-i0+jj)=psig_f(2,i,i2,i3)
            psi_f(3,i-i0+jj)=psig_f(3,i,i2,i3)
            psi_f(4,i-i0+jj)=psig_f(4,i,i2,i3)
            psi_f(5,i-i0+jj)=psig_f(5,i,i2,i3)
            psi_f(6,i-i0+jj)=psig_f(6,i,i2,i3)
            psi_f(7,i-i0+jj)=psig_f(7,i,i2,i3)
          enddo
        enddo
  
          deallocate(wprojx,wprojy,wprojz,work,psig_c,psig_f)

	return
	end



        subroutine crtproj(iproc,nterm,n1,n2,n3, & 
                     nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,  & 
                     radius_f,cpmult,fpmult,hgrid,gau_a,factor,rx,ry,rz,lx,ly,lz, & 
                     mvctr_c,mvctr_f,proj_c,proj_f)
! returns the compressed form of a Gaussian projector 
! x^lx * y^ly * z^lz * exp (-1/(2*gau_a^2) *((x-cntrx)^2 + (y-cntry)^2 + (z-cntrz)^2 ))
! in the arrays proj_c, proj_f
        implicit real*8 (a-h,o-z)
        parameter(ntermx=3,nw=16000)
        parameter(eps_mach=1.d-12)
        dimension lx(3),ly(3),lz(3)
        dimension proj_c(mvctr_c),proj_f(7,mvctr_f)
        real*8, allocatable, dimension(:,:,:) :: wprojx, wprojy, wprojz
        real*8, allocatable, dimension(:,:) :: work

        allocate(wprojx(0:n1,2,ntermx),wprojy(0:n2,2,ntermx),wprojz(0:n3,2,ntermx),work(0:nw,2))


        onem=1.d0-eps_mach
        rad_c=radius_f*cpmult
        rad_f=radius_f*fpmult

! make sure that the coefficients returned by CALL GAUSS_TO_DAUB are zero outside [ml:mr] 
        err_norm=0.d0 
      do 100,iterm=1,nterm
        n_gau=lx(iterm) 
        CALL GAUSS_TO_DAUB(hgrid,factor,rx,gau_a,n_gau,n1,ml1,mu1,wprojx(0,1,iterm),te,work,nw)
        err_norm=max(err_norm,te) 
        n_gau=ly(iterm) 
        CALL GAUSS_TO_DAUB(hgrid,1.d0,ry,gau_a,n_gau,n2,ml2,mu2,wprojy(0,1,iterm),te,work,nw)
        err_norm=max(err_norm,te) 
        n_gau=lz(iterm) 
        CALL GAUSS_TO_DAUB(hgrid,1.d0,rz,gau_a,n_gau,n3,ml3,mu3,wprojz(0,1,iterm),te,work,nw)
        err_norm=max(err_norm,te) 
       if (ml1.gt.min(nl1_c,nl1_f)) write(*,*) 'ERROR: ml1'
       if (ml2.gt.min(nl2_c,nl2_f)) write(*,*) 'ERROR: ml2'
       if (ml3.gt.min(nl3_c,nl3_f)) write(*,*) 'ERROR: ml3'
       if (mu1.lt.max(nu1_c,nu1_f)) write(*,*) 'ERROR: mu1'
       if (mu2.lt.max(nu2_c,nu2_f)) write(*,*) 'ERROR: mu2'
       if (mu3.lt.max(nu3_c,nu3_f)) write(*,*) 'ERROR: mu3'
100   continue
        if (iproc.eq.0) write(*,*) 'max err_norm ',err_norm

! First term: coarse projector components
!          write(*,*) 'rad_c=',rad_c
          mvctr=0
          do i3=nl3_c,nu3_c
          dz2=(i3*hgrid-rz)**2
          do i2=nl2_c,nu2_c
          dy2=(i2*hgrid-ry)**2
          do i1=nl1_c,nu1_c
          dx=i1*hgrid-rx
          if (dx**2+(dy2+dz2).lt.rad_c**2) then
            mvctr=mvctr+1
            proj_c(mvctr)=wprojx(i1,1,1)*wprojy(i2,1,1)*wprojz(i3,1,1)
          endif
          enddo ; enddo ; enddo
          if (mvctr.ne.mvctr_c) then 
            write(*,*)  'mvctr,mvctr_c',mvctr,mvctr_c
            stop 'mvctr >< mvctr_c'
          endif

! First term: fine projector components
          mvctr=0
          do i3=nl3_f,nu3_f
          dz2=(i3*hgrid-rz)**2
          do i2=nl2_f,nu2_f
          dy2=(i2*hgrid-ry)**2
          do i1=nl1_f,nu1_f
          dx=i1*hgrid-rx
          if (dx**2+(dy2+dz2).lt.rad_f**2) then
            mvctr=mvctr+1
            proj_f(1,mvctr)=wprojx(i1,2,1)*wprojy(i2,1,1)*wprojz(i3,1,1)
            proj_f(2,mvctr)=wprojx(i1,1,1)*wprojy(i2,2,1)*wprojz(i3,1,1)
            proj_f(3,mvctr)=wprojx(i1,2,1)*wprojy(i2,2,1)*wprojz(i3,1,1)
            proj_f(4,mvctr)=wprojx(i1,1,1)*wprojy(i2,1,1)*wprojz(i3,2,1)
            proj_f(5,mvctr)=wprojx(i1,2,1)*wprojy(i2,1,1)*wprojz(i3,2,1)
            proj_f(6,mvctr)=wprojx(i1,1,1)*wprojy(i2,2,1)*wprojz(i3,2,1)
            proj_f(7,mvctr)=wprojx(i1,2,1)*wprojy(i2,2,1)*wprojz(i3,2,1)
          endif
          enddo ; enddo ; enddo
          if (mvctr.ne.mvctr_f) stop 'mvctr >< mvctr_f'
  

         do iterm=2,nterm

! Other terms: coarse projector components
         mvctr=0
         do i3=nl3_c,nu3_c
         dz2=(i3*hgrid-rz)**2
         do i2=nl2_c,nu2_c
         dy2=(i2*hgrid-ry)**2
         do i1=nl1_c,nu1_c
         dx=i1*hgrid-rx
         if (dx**2+(dy2+dz2).lt.rad_c**2) then
           mvctr=mvctr+1
           proj_c(mvctr)=proj_c(mvctr)+wprojx(i1,1,iterm)*wprojy(i2,1,iterm)*wprojz(i3,1,iterm)
         endif
         enddo ; enddo ; enddo

! Other terms: fine projector components
         mvctr=0
         do i3=nl3_f,nu3_f
         dz2=(i3*hgrid-rz)**2
         do i2=nl2_f,nu2_f
         dy2=(i2*hgrid-ry)**2
         do i1=nl1_f,nu1_f
         dx=i1*hgrid-rx
         if (dx**2+(dy2+dz2).lt.rad_f**2) then
           mvctr=mvctr+1
           proj_f(1,mvctr)=proj_f(1,mvctr)+wprojx(i1,2,iterm)*wprojy(i2,1,iterm)*wprojz(i3,1,iterm)
           proj_f(2,mvctr)=proj_f(2,mvctr)+wprojx(i1,1,iterm)*wprojy(i2,2,iterm)*wprojz(i3,1,iterm)
           proj_f(3,mvctr)=proj_f(3,mvctr)+wprojx(i1,2,iterm)*wprojy(i2,2,iterm)*wprojz(i3,1,iterm)
           proj_f(4,mvctr)=proj_f(4,mvctr)+wprojx(i1,1,iterm)*wprojy(i2,1,iterm)*wprojz(i3,2,iterm)
           proj_f(5,mvctr)=proj_f(5,mvctr)+wprojx(i1,2,iterm)*wprojy(i2,1,iterm)*wprojz(i3,2,iterm)
           proj_f(6,mvctr)=proj_f(6,mvctr)+wprojx(i1,1,iterm)*wprojy(i2,2,iterm)*wprojz(i3,2,iterm)
           proj_f(7,mvctr)=proj_f(7,mvctr)+wprojx(i1,2,iterm)*wprojy(i2,2,iterm)*wprojz(i3,2,iterm)
         endif
         enddo ; enddo ; enddo
          

          enddo
  
          deallocate(wprojx,wprojy,wprojz,work)

	return
	end


        subroutine applyprojectors(iproc,nproc,ntypes,nat,iatype,psppar,occup, &
                    nprojel,nproj,nseg_p,keyg_p,keyv_p,nvctr_p,proj,  &
                    norb,nseg,keyg,keyv,nvctr,psi,hpsi,eproj_sum)
! Applies all the projectors onto a wavefunction
! Input: psi_c,psi_f
! In/Output: hpsi_c,hpsi_f (both are updated, i.e. not initilized to zero at the beginning)
        implicit real*8 (a-h,o-z)
        parameter(eps_mach=1.d-12)
        dimension psppar(0:2,0:4,ntypes),iatype(nat)
        dimension nseg(0:2*norb),nvctr(0:2*norb)
        dimension keyg(2,nseg(2*norb)),keyv(nseg(2*norb))
        dimension psi(nvctr(2*norb)),hpsi(nvctr(2*norb))
        dimension nseg_p(0:2*nat),nvctr_p(0:2*nat)
        dimension keyg_p(2,nseg_p(2*nat)),keyv_p(nseg_p(2*nat))
        dimension proj(nprojel),occup(norb)

  eproj_sum=0.d0
  onem=1.d0-eps_mach
  norb_p=int(onem+dble(norb)/dble(nproc))
! loop over all my orbitals
  do 100,iorb=iproc*norb_p+1,min((iproc+1)*norb_p,norb)
    maseg_c=nseg(2*iorb-1)-nseg(2*iorb-2)
    maseg_f=nseg(2*iorb  )-nseg(2*iorb-1)
    iseg_c=nseg(2*iorb-2)+1
    iseg_f=nseg(2*iorb-1)+1
    mavctr_c= nvctr(2*iorb-1)-nvctr(2*iorb-2)
    mavctr_f=(nvctr(2*iorb  )-nvctr(2*iorb-1))/7
    ipsi_c=nvctr(2*iorb-2)+1
    ipsi_f=nvctr(2*iorb-1)+1

! loop over all projectors
    iproj=0
    eproj=0.d0
    istart_c=1
    do iat=1,nat
        mbseg_c=nseg_p(2*iat-1)-nseg_p(2*iat-2)
        mbseg_f=nseg_p(2*iat  )-nseg_p(2*iat-1)
        jseg_c=nseg_p(2*iat-2)+1
        jseg_f=nseg_p(2*iat-1)+1
        mbvctr_c=nvctr_p(2*iat-1)-nvctr_p(2*iat-2)
        mbvctr_f=nvctr_p(2*iat  )-nvctr_p(2*iat-1)
     ityp=iatype(iat)
! ONLY GTH PSP form (not HGH)
     do i=1,2
     do j=1,2
     if (psppar(i,j,ityp).ne.0.d0) then
     do l=1,2*i-1
        iproj=iproj+1
        istart_f=istart_c+mbvctr_c
        call wdot(  &
                   mavctr_c,mavctr_f,maseg_c,maseg_f,keyv(iseg_c),keyv(iseg_f),  &
                   keyg(1,iseg_c),keyg(1,iseg_f),psi(ipsi_c),psi(ipsi_f),  &
                   mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p(jseg_c),keyv_p(jseg_f),  &
                   keyg_p(1,jseg_c),keyg_p(1,jseg_f),proj(istart_c),proj(istart_f),scpr)
! test
        call wdot(  &
                   mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p(jseg_c),keyv_p(jseg_f),  &
                   keyg_p(1,jseg_c),keyg_p(1,jseg_f),proj(istart_c),proj(istart_f),  &
                   mavctr_c,mavctr_f,maseg_c,maseg_f,keyv(iseg_c),keyv(iseg_f),  &
                   keyg(1,iseg_c),keyg(1,iseg_f),psi(ipsi_c),psi(ipsi_f),tcpr)
        if (scpr.ne.tcpr) stop 'projectors: scpr.ne.tcpr'
! testend

        scprp=scpr*psppar(i,j,ityp)
!        write(*,*) 'pdot scpr',scpr,psppar(i,j,ityp)
        eproj=eproj+scprp*scpr

        call waxpy(  &
             scprp,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p(jseg_c),keyv_p(jseg_f),  &
                   keyg_p(1,jseg_c),keyg_p(1,jseg_f),proj(istart_c),proj(istart_f),  &
                   mavctr_c,mavctr_f,maseg_c,maseg_f,keyv(iseg_c),keyv(iseg_f),  &
                   keyg(1,iseg_c),keyg(1,iseg_f),hpsi(ipsi_c),hpsi(ipsi_f))

        istart_c=istart_f+7*mbvctr_f
     enddo
     endif
     enddo
     enddo
   enddo
     write(*,*) 'iorb,eproj',iorb,eproj
     eproj_sum=eproj_sum+occup(iorb)*eproj
     if (iproj.ne.nproj) stop '1:applyprojectors'
     if (istart_c-1.ne.nprojel) stop '2:applyprojectors'
100 continue

         return
         end



	subroutine wcopy(mvctr_c,mvctr_f,bpsi_c,bpsi_f,apsi_c,apsi_f)
! multiplies a wavefunction psi_c,psi_f (in vector form) with a scalar (scal)
        implicit real*8 (a-h,o-z)
        dimension apsi_c(mvctr_c),apsi_f(7,mvctr_f)
        dimension bpsi_c(mvctr_c),bpsi_f(7,mvctr_f)

	do i=1,mvctr_c
           apsi_c(i)=bpsi_c(i)
        enddo
	do i=1,mvctr_f
           apsi_f(1,i)=bpsi_f(1,i)
           apsi_f(2,i)=bpsi_f(2,i)
           apsi_f(3,i)=bpsi_f(3,i)
           apsi_f(4,i)=bpsi_f(4,i)
           apsi_f(5,i)=bpsi_f(5,i)
           apsi_f(6,i)=bpsi_f(6,i)
           apsi_f(7,i)=bpsi_f(7,i)
        enddo

	return
	end


	subroutine wzero(mvctr_c,mvctr_f,psi_c,psi_f)
! multiplies a wavefunction psi_c,psi_f (in vector form) with a scalar (scal)
        implicit real*8 (a-h,o-z)
        dimension psi_c(mvctr_c),psi_f(7,mvctr_f)

	do i=1,mvctr_c
           psi_c(i)=0.d0
        enddo
	do i=1,mvctr_f
           psi_f(1,i)=0.d0
           psi_f(2,i)=0.d0
           psi_f(3,i)=0.d0
           psi_f(4,i)=0.d0
           psi_f(5,i)=0.d0
           psi_f(6,i)=0.d0
           psi_f(7,i)=0.d0
        enddo

	return
	end



	subroutine wscal(mvctr_c,mvctr_f,scal,psi_c,psi_f)
! multiplies a wavefunction psi_c,psi_f (in vector form) with a scalar (scal)
        implicit real*8 (a-h,o-z)
        dimension psi_c(mvctr_c),psi_f(7,mvctr_f)

	do i=1,mvctr_c
           psi_c(i)=psi_c(i)*scal
        enddo
	do i=1,mvctr_f
           psi_f(1,i)=psi_f(1,i)*scal
           psi_f(2,i)=psi_f(2,i)*scal
           psi_f(3,i)=psi_f(3,i)*scal
           psi_f(4,i)=psi_f(4,i)*scal
           psi_f(5,i)=psi_f(5,i)*scal
           psi_f(6,i)=psi_f(6,i)*scal
           psi_f(7,i)=psi_f(7,i)*scal
        enddo

	return
	end



	subroutine wnrm(mvctr_c,mvctr_f,psi_c,psi_f,scpr)
! calculates the norm SQUARED (scpr) of a wavefunction (in vector form)
        implicit real*8 (a-h,o-z)
        dimension psi_c(mvctr_c),psi_f(7,mvctr_f)

        scpr=0.d0
	do i=1,mvctr_c
           scpr=scpr+psi_c(i)**2
        enddo
        scpr1=0.d0
        scpr2=0.d0
        scpr3=0.d0
        scpr4=0.d0
        scpr5=0.d0
        scpr6=0.d0
        scpr7=0.d0
	do i=1,mvctr_f
           scpr1=scpr1+psi_f(1,i)**2
           scpr2=scpr2+psi_f(2,i)**2
           scpr3=scpr3+psi_f(3,i)**2
           scpr4=scpr4+psi_f(4,i)**2
           scpr5=scpr5+psi_f(5,i)**2
           scpr6=scpr6+psi_f(6,i)**2
           scpr7=scpr7+psi_f(7,i)**2
        enddo
        scpr=scpr+scpr1+scpr2+scpr3+scpr4+scpr5+scpr6+scpr7

	return
	end



	subroutine wdot(  & 
        mavctr_c,mavctr_f,maseg_c,maseg_f,keyav_c,keyav_f,keyag_c,keyag_f,apsi_c,apsi_f,  & 
        mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keybv_c,keybv_f,keybg_c,keybg_f,bpsi_c,bpsi_f,scpr)
! calculates the dot product between two wavefunctions in compressed form
        implicit real*8 (a-h,o-z)
        dimension keyav_c(maseg_c),keyag_c(2,maseg_c),keyav_f(maseg_f),keyag_f(2,maseg_f)
        dimension keybv_c(mbseg_c),keybg_c(2,mbseg_c),keybv_f(mbseg_f),keybg_f(2,mbseg_f)
        dimension apsi_c(mavctr_c),apsi_f(7,mavctr_f),bpsi_c(mbvctr_c),bpsi_f(7,mbvctr_f)

!        llc=0
        scpr=0.d0
! coarse part
        ibseg=1
        do iaseg=1,maseg_c
          jaj=keyav_c(iaseg)
          ja0=keyag_c(1,iaseg)
          ja1=keyag_c(2,iaseg)

100       jb1=keybg_c(2,ibseg)
          if (jb1.lt.ja0) then
             ibseg=ibseg+1
             if (ibseg.gt.mbseg_c) goto 111
             goto 100
          endif
          jb0=keybg_c(1,ibseg)
          jbj=keybv_c(ibseg)
          if (ja0 .gt. jb0) then 
             iaoff=0
             iboff=ja0-jb0
             length=min(ja1,jb1)-ja0
          else
             iaoff=jb0-ja0
             iboff=0
             length=min(ja1,jb1)-jb0
          endif
!           write(*,*) 'ja0,ja1,jb0,jb1',ja0,ja1,jb0,jb1,length
!          write(*,'(5(a,i5))') 'C:from ',jaj+iaoff,' to ',jaj+iaoff+length,' and from ',jbj+iboff,' to ',jbj+iboff+length
          do i=0,length
!          llc=llc+1
          scpr=scpr+apsi_c(jaj+iaoff+i)*bpsi_c(jbj+iboff+i) 
          enddo
        enddo
111     continue


!        llf=0
        scpr1=0.d0
        scpr2=0.d0
        scpr3=0.d0
        scpr4=0.d0
        scpr5=0.d0
        scpr6=0.d0
        scpr7=0.d0
! fine part
        ibseg=1
        do iaseg=1,maseg_f
          jaj=keyav_f(iaseg)
          ja0=keyag_f(1,iaseg)
          ja1=keyag_f(2,iaseg)

200       jb1=keybg_f(2,ibseg)
          if (jb1.lt.ja0) then
             ibseg=ibseg+1
             if (ibseg.gt.mbseg_f) goto 222
             goto 200
          endif
          jb0=keybg_f(1,ibseg)
          jbj=keybv_f(ibseg)
          if (ja0 .gt. jb0) then 
             iaoff=0
             iboff=ja0-jb0
             length=min(ja1,jb1)-ja0
          else
             iaoff=jb0-ja0
             iboff=0
             length=min(ja1,jb1)-jb0
          endif
          do i=0,length
!          llf=llf+1
          scpr1=scpr1+apsi_f(1,jaj+iaoff+i)*bpsi_f(1,jbj+iboff+i) 
          scpr2=scpr2+apsi_f(2,jaj+iaoff+i)*bpsi_f(2,jbj+iboff+i) 
          scpr3=scpr3+apsi_f(3,jaj+iaoff+i)*bpsi_f(3,jbj+iboff+i) 
          scpr4=scpr4+apsi_f(4,jaj+iaoff+i)*bpsi_f(4,jbj+iboff+i) 
          scpr5=scpr5+apsi_f(5,jaj+iaoff+i)*bpsi_f(5,jbj+iboff+i) 
          scpr6=scpr6+apsi_f(6,jaj+iaoff+i)*bpsi_f(6,jbj+iboff+i) 
          scpr7=scpr7+apsi_f(7,jaj+iaoff+i)*bpsi_f(7,jbj+iboff+i) 
          enddo
        enddo
222     continue

        scpr=scpr+scpr1+scpr2+scpr3+scpr4+scpr5+scpr6+scpr7
!        write(*,*) 'llc,llf',llc,llf

	return
	end



	subroutine waxpy(  & 
        scpr,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keybv_c,keybv_f,keybg_c,keybg_f,bpsi_c,bpsi_f, & 
        mavctr_c,mavctr_f,maseg_c,maseg_f,keyav_c,keyav_f,keyag_c,keyag_f,apsi_c,apsi_f)
! rank 1 update of wavefunction a with wavefunction b: apsi=apsi+scpr*bpsi
! The update is only done in the localization region of apsi
        implicit real*8 (a-h,o-z)
        dimension keyav_c(maseg_c),keyag_c(2,maseg_c),keyav_f(maseg_f),keyag_f(2,maseg_f)
        dimension keybv_c(mbseg_c),keybg_c(2,mbseg_c),keybv_f(mbseg_f),keybg_f(2,mbseg_f)
        dimension apsi_c(mavctr_c),apsi_f(7,mavctr_f),bpsi_c(mbvctr_c),bpsi_f(7,mbvctr_f)

!        llc=0
! coarse part
        ibseg=1
        do iaseg=1,maseg_c
          jaj=keyav_c(iaseg)
          ja0=keyag_c(1,iaseg)
          ja1=keyag_c(2,iaseg)

100       jb1=keybg_c(2,ibseg)
          if (jb1.lt.ja0) then
             ibseg=ibseg+1
             if (ibseg.gt.mbseg_c) goto 111
             goto 100
          endif
          jb0=keybg_c(1,ibseg)
          jbj=keybv_c(ibseg)
          if (ja0 .gt. jb0) then 
             iaoff=0
             iboff=ja0-jb0
             length=min(ja1,jb1)-ja0
          else
             iaoff=jb0-ja0
             iboff=0
             length=min(ja1,jb1)-jb0
          endif
          do i=0,length
!          llc=llc+1
          apsi_c(jaj+iaoff+i)=apsi_c(jaj+iaoff+i)+scpr*bpsi_c(jbj+iboff+i) 
          enddo
        enddo
111     continue

!        llf=0
! fine part
        ibseg=1
        do iaseg=1,maseg_f
          jaj=keyav_f(iaseg)
          ja0=keyag_f(1,iaseg)
          ja1=keyag_f(2,iaseg)

200       jb1=keybg_f(2,ibseg)
          if (jb1.lt.ja0) then
             ibseg=ibseg+1
             if (ibseg.gt.mbseg_f) goto 222
             goto 200
          endif
          jb0=keybg_f(1,ibseg)
          jbj=keybv_f(ibseg)
          if (ja0 .gt. jb0) then 
             iaoff=0
             iboff=ja0-jb0
             length=min(ja1,jb1)-ja0
          else
             iaoff=jb0-ja0
             iboff=0
             length=min(ja1,jb1)-jb0
          endif
          do i=0,length
!          llf=llf+1
          apsi_f(1,jaj+iaoff+i)=apsi_f(1,jaj+iaoff+i)+scpr*bpsi_f(1,jbj+iboff+i) 
          apsi_f(2,jaj+iaoff+i)=apsi_f(2,jaj+iaoff+i)+scpr*bpsi_f(2,jbj+iboff+i) 
          apsi_f(3,jaj+iaoff+i)=apsi_f(3,jaj+iaoff+i)+scpr*bpsi_f(3,jbj+iboff+i) 
          apsi_f(4,jaj+iaoff+i)=apsi_f(4,jaj+iaoff+i)+scpr*bpsi_f(4,jbj+iboff+i) 
          apsi_f(5,jaj+iaoff+i)=apsi_f(5,jaj+iaoff+i)+scpr*bpsi_f(5,jbj+iboff+i) 
          apsi_f(6,jaj+iaoff+i)=apsi_f(6,jaj+iaoff+i)+scpr*bpsi_f(6,jbj+iboff+i) 
          apsi_f(7,jaj+iaoff+i)=apsi_f(7,jaj+iaoff+i)+scpr*bpsi_f(7,jbj+iboff+i) 
          enddo
        enddo
222     continue
!        write(*,*) 'waxpy,llc,llf',llc,llf

	return
	end




       subroutine fill_logrid(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,nat,  &
                               ntypes,iatype,rxyz,radii,rmult,hgrid,loregion,logrid)
! set up an array logrid(i1,i2,i3) that specifies whether the grid point
! i1,i2,i3 is the center of a scaling function/wavelet
        implicit real*8 (a-h,o-z)
        logical logrid, loregion
        parameter(eps_mach=1.d-12)
        dimension rxyz(3,nat),iatype(nat),radii(ntypes)
        dimension logrid(0:n1,0:n2,0:n3),loregion(nat)

        do i3=nl3,nu3 ; do i2=nl2,nu2 ; do i1=nl1,nu1
         logrid(i1,i2,i3)=.false.
        enddo ; enddo ; enddo

      do iat=1,nat
      if (loregion(iat)) then
        rad=radii(iatype(iat))*rmult
!        write(*,*) 'iat,nat,rad',iat,nat,rad
        onem=1.d0-eps_mach
        ml1=int(onem+(rxyz(1,iat)-rad)/hgrid)  ; mu1=int((rxyz(1,iat)+rad)/hgrid)
        ml2=int(onem+(rxyz(2,iat)-rad)/hgrid)  ; mu2=int((rxyz(2,iat)+rad)/hgrid)
        ml3=int(onem+(rxyz(3,iat)-rad)/hgrid)  ; mu3=int((rxyz(3,iat)+rad)/hgrid)
!        write(*,'(a,6(i4))') 'fill grid',ml1,mu1,ml2,mu2,ml3,mu3
        if (ml1.lt.nl1) stop 'ml1 < nl1' ; if (mu1.gt.nu1) stop 'mu1 > nu1'
        if (ml2.lt.nl2) stop 'ml2 < nl2' ; if (mu2.gt.nu2) stop 'mu2 > nu2'
        if (ml3.lt.nl3) stop 'ml3 < nl3' ; if (mu3.gt.nu3) stop 'mu3 > nu3'
        do i3=ml3,mu3
        dz2=(i3*hgrid-rxyz(3,iat))**2
        do i2=ml2,mu2
        dy2=(i2*hgrid-rxyz(2,iat))**2
        do i1=ml1,mu1
        dx=i1*hgrid-rxyz(1,iat)
        if (dx**2+(dy2+dz2).lt.rad**2) logrid(i1,i2,i3)=.true.
        enddo ; enddo ; enddo
      endif
      enddo

        return
        end



	subroutine num_segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,mvctr)
! Calculates the keys describing a wavefunction data structure
        implicit real*8 (a-h,o-z)
        logical logrid,plogrid
        dimension logrid(0:n1,0:n2,0:n3)

        mvctr=0
        nsrt=0
        nend=0
        do i3=nl3,nu3 ; do i2=nl2,nu2

        plogrid=.false.
        do i1=nl1,nu1
         if (logrid(i1,i2,i3)) then
           mvctr=mvctr+1
           if (plogrid .eqv. .false.) then
             nsrt=nsrt+1
           endif
         else
           if (plogrid .eqv. .true.) then
             nend=nend+1
           endif
         endif
         plogrid=logrid(i1,i2,i3)
        enddo 
           if (plogrid .eqv. .true.) then
             nend=nend+1
           endif
        enddo ; enddo
        if (nend.ne.nsrt) then 
           write(*,*) 'nend , nsrt',nend,nsrt
           stop 'nend <> nsrt'
        endif
        mseg=nend

	return
        end


	subroutine segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
                           logrid,mseg,keyg,keyv)
! Calculates the keys describing a wavefunction data structure
        implicit real*8 (a-h,o-z)
        logical logrid,plogrid
        dimension logrid(0:n1,0:n2,0:n3),keyg(2,mseg),keyv(mseg)

        mvctr=0
        nsrt=0
        nend=0
        do i3=nl3,nu3 ; do i2=nl2,nu2

        plogrid=.false.
        do i1=nl1,nu1
         ngridp=i3*((n1+1)*(n2+1)) + i2*(n1+1) + i1+1
         if (logrid(i1,i2,i3)) then
           mvctr=mvctr+1
           if (plogrid .eqv. .false.) then
             nsrt=nsrt+1
             keyg(1,nsrt)=ngridp
             keyv(nsrt)=mvctr
           endif
         else
           if (plogrid .eqv. .true.) then
             nend=nend+1
             keyg(2,nend)=ngridp-1
           endif
         endif
         plogrid=logrid(i1,i2,i3)
        enddo 
           if (plogrid .eqv. .true.) then
             nend=nend+1
             keyg(2,nend)=ngridp
           endif
        enddo ; enddo
        if (nend.ne.nsrt) then 
           write(*,*) 'nend , nsrt',nend,nsrt
           stop 'nend <> nsrt'
        endif
        mseg=nend

	return
        end


        subroutine compress(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
                            mseg_c,mvctr_c,keyg_c,keyv_c,  & 
                            mseg_f,mvctr_f,keyg_f,keyv_f,  & 
                            psifscf,psi_c,psi_f)
! Compresses a wavefunction that is given in terms of fine scaling functions (psifscf) into 
! the retained coarse scaling functions and wavelet coefficients (psi_c,psi_f)
        implicit real*8 (a-h,o-z)
        dimension keyg_c(2,mseg_c),keyv_c(mseg_c),keyg_f(2,mseg_f),keyv_f(mseg_f)
        dimension psi_c(mvctr_c),psi_f(7,mvctr_f)
        dimension psifscf((2*(nu1-nl1)+16)*(2*(nu2-nl2)+16)*(2*(nu3-nl3)+16))
        real*8, allocatable :: psig(:,:,:,:,:,:),ww(:)

        allocate(psig(nl1:nu1,2,nl2:nu2,2,nl3:nu3,2),  & 
                 ww((2*(nu1-nl1)+16)*(2*(nu2-nl2)+16)*(2*(nu3-nl3)+16)))
! decompose wavelets into coarse scaling functions and wavelets
	call analyse_shrink(nu1-nl1,nu2-nl2,nu3-nl3,ww,psifscf,psig)

! coarse part
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
            psi_c(i-i0+jj)=psig(i,1,i2,1,i3,1)
          enddo
        enddo

! fine part
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
            psi_f(1,i-i0+jj)=psig(i,2,i2,1,i3,1)
            psi_f(2,i-i0+jj)=psig(i,1,i2,2,i3,1)
            psi_f(3,i-i0+jj)=psig(i,2,i2,2,i3,1)
            psi_f(4,i-i0+jj)=psig(i,1,i2,1,i3,2)
            psi_f(5,i-i0+jj)=psig(i,2,i2,1,i3,2)
            psi_f(6,i-i0+jj)=psig(i,1,i2,2,i3,2)
            psi_f(7,i-i0+jj)=psig(i,2,i2,2,i3,2)
          enddo
        enddo

        deallocate(psig,ww)

	end


        subroutine uncompress(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
                              mseg_c,mvctr_c,keyg_c,keyv_c,  & 
                              mseg_f,mvctr_f,keyg_f,keyv_f,  & 
                              psi_c,psi_f,psifscf)
! Expands the compressed wavefunction in vector form (psi_c,psi_f) 
! into fine scaling functions (psifscf)
        implicit real*8 (a-h,o-z)
        dimension keyg_c(2,mseg_c),keyv_c(mseg_c),keyg_f(2,mseg_f),keyv_f(mseg_f)
        dimension psi_c(mvctr_c),psi_f(7,mvctr_f)
        dimension psifscf((2*(nu1-nl1)+16)*(2*(nu2-nl2)+16)*(2*(nu3-nl3)+16))
        real*8, allocatable :: psig(:,:,:,:,:,:),ww(:)

        allocate(psig(nl1:nu1,2,nl2:nu2,2,nl3:nu3,2),  & 
                 ww((2*(nu1-nl1)+16)*(2*(nu2-nl2)+16)*(2*(nu3-nl3)+16)))

        call zero(8*(nu1-nl1+1)*(nu2-nl2+1)*(nu3-nl3+1),psig)

! coarse part
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
            psig(i,1,i2,1,i3,1)=psi_c(i-i0+jj)
          enddo
         enddo

! fine part
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
            psig(i,2,i2,1,i3,1)=psi_f(1,i-i0+jj)
            psig(i,1,i2,2,i3,1)=psi_f(2,i-i0+jj)
            psig(i,2,i2,2,i3,1)=psi_f(3,i-i0+jj)
            psig(i,1,i2,1,i3,2)=psi_f(4,i-i0+jj)
            psig(i,2,i2,1,i3,2)=psi_f(5,i-i0+jj)
            psig(i,1,i2,2,i3,2)=psi_f(6,i-i0+jj)
            psig(i,2,i2,2,i3,2)=psi_f(7,i-i0+jj)
          enddo
         enddo

! calculate fine scaling functions
	call synthese_grow(nu1-nl1,nu2-nl2,nu3-nl3,ww,psig,psifscf)

        deallocate(psig,ww)

	end




	subroutine synthese_grow(n1,n2,n3,ww,x,y)
! A synthesis wavelet transformation where the size of the data is allowed to grow
! The input array x is not overwritten
        implicit real*8 (a-h,o-z)
        dimension x(0:n1,2,0:n2,2,0:n3,2)
        dimension ww(-7:2*n2+8,-7:2*n3+8,-7:2*n1+8)
        dimension  y(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8)

! i1,i2,i3 -> i2,i3,I1
        nt=(2*n2+2)*(2*n3+2)
        call  syn_rot_grow(n1,nt,x,y)
! i2,i3,I1 -> i3,I1,I2
        nt=(2*n3+2)*(2*n1+16)
        call  syn_rot_grow(n2,nt,y,ww)
! i3,I1,I2  -> I1,I2,I3
        nt=(2*n1+16)*(2*n2+16)
        call  syn_rot_grow(n3,nt,ww,y)

	return
        end


	subroutine analyse_shrink(n1,n2,n3,ww,y,x)
! A analysis wavelet transformation where the size of the data is forced to shrink
! The input array y is overwritten
        implicit real*8 (a-h,o-z)
        dimension ww(-7:2*n2+8,-7:2*n3+8,-7:2*n1+8)
        dimension  y(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8)
        dimension x(0:n1,2,0:n2,2,0:n3,2)

! I1,I2,I3 -> I2,I3,i1
        nt=(2*n2+16)*(2*n3+16)
        call  ana_rot_shrink(n1,nt,y,ww)
! I2,I3,i1 -> I3,i1,i2
        nt=(2*n3+16)*(2*n1+2)
        call  ana_rot_shrink(n2,nt,ww,y)
! I3,i1,i2 -> i1,i2,i3
        nt=(2*n1+2)*(2*n2+2)
        call  ana_rot_shrink(n3,nt,y,x)

	return
        end


        subroutine syn_rot_grow(n,nt,x,y)
        implicit real*8 (a-h,o-z)
!  Daubechies_16 symetric 
!        REAL*8 CH(-M:M)/0.d0,&
!        -0.0033824159510050025955D0,-0.00054213233180001068935D0,&
!        0.031695087811525991431D0,0.0076074873249766081919D0,&
!        -0.14329423835127266284D0,-0.061273359067811077843D0,&
!        0.48135965125905339159D0,0.77718575169962802862D0,0.36444189483617893676D0,&
!        -0.051945838107881800736D0,-0.027219029917103486322D0,&
!        0.049137179673730286787D0,0.0038087520138944894631D0,&
!        -0.014952258337062199118D0,-0.00030292051472413308126D0,&
!        0.0018899503327676891843D0&

        dimension x(0:2*n+1,nt),y(nt,-7:2*n+8)
        include 'sym_16.inc'

        do 200,l=1,nt

       DO I=-7,2*n+8
         y(l,I)=0.D0
       enddo

       DO I=0,n
         I2=2*I
         DO J=-M+1,M
           y(l,J+I2)=y(l,J+I2)+CH(J)*x(I,l)+CG(J)*x(n+1+I,l)
         ENDDO
       ENDDO

200     continue

        return
        end



        subroutine ana_rot_shrink(n,nt,x,y)
        implicit real*8 (a-h,o-z)
!  Daubechies_8 symetric 
!        parameter(chm3=.230377813308896501d0, chm2=.714846570552915647d0, & 
!                  chm1=.630880767929858908d0, ch00=-.279837694168598542d-1, & 
!                  chp1=-.187034811719093084d0, chp2=.308413818355607636d-1, & 
!                  chp3=.328830116668851997d-1, chp4=-.105974017850690321d-1)


       dimension x(-7:2*n+8,nt),y(nt,0:2*n+1)
       INCLUDE 'sym_16.inc'

       do 200,l=1,nt

       DO I=0,n
         I2=2*I
         CI=0.D0
         DI=0.D0
         DO J=-M+1,M
           CI=CI+CHT(J)*x(J+I2,l)
           DI=DI+CGT(J)*x(J+I2,l)
         ENDDO
         y(l,I)=CI
         y(l,n+1+I)=DI
       ENDDO
       
200     continue

        return
        end



        subroutine convolut_magic_n(n1,n2,n3,x,y)
! Applies the magic filter matrix ( no transposition) ; data set grows
! The input array x is not overwritten
        implicit real*8 (a-h,o-z)
        parameter(lowfil=-8,lupfil=7) ! has to be consistent with values in convrot
        dimension x(0:n1,0:n2,0:n3),y(-lupfil:n1-lowfil,-lupfil:n2-lowfil,-lupfil:n3-lowfil)
        real*8, allocatable :: ww(:,:,:)

        allocate(ww(-lupfil:n1-lowfil,-lupfil:n2-lowfil,0:n3))
 
!  (i1,i2*i3) -> (i2*i3,I1)
        ndat=(n2+1)*(n3+1)
        call convrot_grow(n1,ndat,x,y)
!  (i2,i3*I1) -> (i3*i1,I2)
        ndat=(n3+1)*(n1+1+lupfil-lowfil)
        call convrot_grow(n2,ndat,y,ww)
!  (i3,I1*I2) -> (iI*I2,I3)
        ndat=(n1+1+lupfil-lowfil)*(n2+1+lupfil-lowfil)
        call convrot_grow(n3,ndat,ww,y)

        deallocate(ww)

	return
	end


        subroutine convolut_magic_t(n1,n2,n3,x,y)
! Applies the magic filter matrix transposed ; data set shrinks
! The input array x is overwritten
        implicit real*8 (a-h,o-z)
        parameter(lowfil=-8,lupfil=7) ! has to be consistent with values in convrot
        dimension x(-lupfil:n1-lowfil,-lupfil:n2-lowfil,-lupfil:n3-lowfil),y(0:n1,0:n2,0:n3)
        real*8, allocatable :: ww(:,:,:)

        allocate(ww(0:n1,-lupfil:n2-lowfil,-lupfil:n3-lowfil))

!  (I1,I2*I3) -> (I2*I3,i1)
        ndat=(n2+1+lupfil-lowfil)*(n3+1+lupfil-lowfil)
        call convrot_shrink(n1,ndat,x,ww)
!  (I2,I3*i1) -> (I3*i1,i2)
        ndat=(n3+1+lupfil-lowfil)*(n1+1)
        call convrot_shrink(n2,ndat,ww,x)
!  (I3,i1*i2) -> (i1*i2,i3)
        ndat=(n1+1)*(n2+1)
        call convrot_shrink(n3,ndat,x,y)

        deallocate(ww)

	return
	end


!   subroutine zero(n,x)
!     implicit real*8 (a-h,o-z)
!     dimension x(n)
!     do i=1,n
!        x(i)=0.d0
!     end do
!   end subroutine zero

        subroutine plot_psifscf(iunit,hgrid,nl1,nu1,nl2,nu2,nl3,nu3,psifscf)
        implicit real*8 (a-h,o-z)
        dimension psifscf(-7+2*nl1:2*nu1+8,-7+2*nl2:2*nu2+8,-7+2*nl3:2*nu3+8)

        hgridh=.5d0*hgrid

! 111 diagonal
        do i=max(-7+2*nl1,-7+2*nl2,-7+2*nl3),min(2*nu1+8,2*nu2+8,2*nu3+8)
            i1=i ; i2=i ; i3=i
            write(iunit,'(3(1x,e10.3),1x,e12.5)') i1*hgridh,i2*hgridh,i3*hgridh,psifscf(i1,i2,i3)
        enddo

        nc1=2*(nu1-nl1)+1
        nc2=2*(nu2-nl2)+1
        nc3=2*(nu3-nl3)+1
! 1-1-1 diagonal
        do i=max(-7+2*nl1,-7+2*nl2,-7+2*nl3),min(2*nu1+8,2*nu2+8,2*nu3+8)
            i1=i ; i2=nc2-i ; i3=nc3-i
            write(iunit,'(3(1x,e10.3),1x,e12.5))') i1*hgridh,i2*hgridh,i3*hgridh,psifscf(i1,i2,i3)
        enddo

! -11-1 diagonal
        do i=max(-7+2*nl1,-7+2*nl2,-7+2*nl3),min(2*nu1+8,2*nu2+8,2*nu3+8)
            i1=nc1-i ; i2=i ; i3=nc3-i
            write(iunit,'(3(1x,e10.3),1x,e12.5)') i1*hgridh,i2*hgridh,i3*hgridh,psifscf(i1,i2,i3)
        enddo

! -1-11 diagonal
        do i=max(-7+2*nl1,-7+2*nl2,-7+2*nl3),min(2*nu1+8,2*nu2+8,2*nu3+8)
            i1=nc1-i ; i2=nc2-i ; i3=i
            write(iunit,'(3(1x,e10.3),1x,e12.5)') i1*hgridh,i2*hgridh,i3*hgridh,psifscf(i1,i2,i3)
        enddo

        return
        end

