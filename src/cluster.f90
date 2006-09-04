

        implicit real*8 (a-h,o-z)
! atomic coordinates, forces
        real*8, allocatable, dimension(:,:) :: rxyz, fxyz
        logical parallel
        character*20 tatonam
! atomic types
        integer, allocatable, dimension(:) :: iatype
        character*20 :: atomnames(100), units
!$      interface
!$        integer ( kind=4 ) function omp_get_num_threads ( )
!$        end function omp_get_num_threads
!$      end interface
!$      interface
!$        integer ( kind=4 ) function omp_get_thread_num ( )
!$        end function omp_get_thread_num
!$      end interface
        include 'mpif.h'

        parallel=.false.


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

!$omp parallel private(iam)  shared (npr)
!$       iam=omp_get_thread_num()
!$       if (iam.eq.0) npr=omp_get_num_threads()
!$       write(*,*) 'iproc,iam,npr',iproc,iam,npr
!$omp end parallel


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

           ampl=2.d-1  ! amplitude for random displacement away from input file geometry (usually equilibrium geom.)
           do iat=1,nat
              call random_number(tt)
              rxyz(1,iat)=rxyz(1,iat)+ampl*tt
              call random_number(tt)
              rxyz(2,iat)=rxyz(2,iat)+ampl*tt
              call random_number(tt)
              rxyz(3,iat)=rxyz(3,iat)+ampl*tt
           enddo
! geometry optimization
        betax=5.d0
        beta=.75d0*betax
        energyold=1.d100
       fluctsum=0.d0
       ngeostep=100
       do 500, igeostep=1,ngeostep

        if (parallel) call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call cluster(parallel,nproc,iproc,nat,ntypes,iatype,atomnames,rxyz,energy,fxyz)

        
        if (energy.gt.energyold) then 
           beta=.5d0*beta
        else
           beta=min(1.05d0*beta,betax)
        endif


           sumx=0.d0
           sumy=0.d0
           sumz=0.d0
           sum=0.d0
           do iat=1,nat
              rxyz(1,iat)=rxyz(1,iat)+beta*fxyz(1,iat)
              rxyz(2,iat)=rxyz(2,iat)+beta*fxyz(2,iat)
              rxyz(3,iat)=rxyz(3,iat)+beta*fxyz(3,iat)
              sumx=sumx+fxyz(1,iat)
              sumy=sumy+fxyz(2,iat)
              sumz=sumz+fxyz(3,iat)
              sum=sum+fxyz(1,iat)**2+fxyz(2,iat)**2+fxyz(3,iat)**2
              if (iproc.eq.0) write(*,'(a,i3,3(x,e14.7))') 'fxyz ',iat,(fxyz(j,iat),j=1,3)
           end do
           fluctsum=fluctsum+sqrt(sumx**2+sumy**2+sumz**2)
           fluct=fluctsum/igeostep
        if (iproc.eq.0) then
           write(*,'(a,x,e21.14,x,e10.3)')'ANALYSIS OF FORCES energy, beta',energy,beta
           write(*,'(a,3(x,e11.4))')'the norm of the forces is', sqrt(sum),sqrt(sum/nat),sqrt(sum/(3*nat))
           write(*,*) 'fluct',fluct
           write(*,*)'the sum of the forces is'
           write(*,'(a16,3x,e16.8)')'x direction',sumx
           write(*,'(a16,3x,e16.8)')'y direction',sumy
           write(*,'(a16,3x,e16.8)')'z direction',sumz
        endif

        if (sqrt(sum).lt.fluct) then 
        if (iproc.eq.0) then
           write(*,*) 'Final positions'
           do iat=1,nat
              write(*,'(i3,3(x,e14.7))') iat,(rxyz(j,iat),j=1,3)
           enddo
        endif
           write(*,*) 'No better convergence possible'
           goto 501
        endif
        energyold=energy
500 continue
501 continue

        deallocate(rxyz,iatype,fxyz)

        if (parallel) call MPI_FINALIZE(ierr)

	end


        subroutine cluster(parallel,nproc,iproc,nat,ntypes,iatype,atomnames,rxyz,energy,fxyz)
! does an electronic structure calculation. Output is the total energy and the forces 
        implicit real*8 (a-h,o-z)
        character*30 label ; character*27 filename ; character*20 atomnames 
        logical logrid,calc_inp_wf,parallel
        logical :: firstrun=.true.
        parameter(eps_mach=1.d-12,onem=1.d0-eps_mach)
! work array for ALLREDUCE
        dimension wrkallred(5,2) 
! work array for crtproj
        dimension lx(3),ly(3),lz(3)
! atomic coordinates
        dimension rxyz(3,nat),fxyz(3,nat),iatype(nat),atomnames(100)
        allocatable :: gxyz(:,:)
! active grid points, segments of real space grid
        allocatable :: logrid(:,:,:)
! occupation numbers
        allocatable :: occup(:)

! wavefunction segments
        allocatable :: keyv(:)
! wavefunction segments on real space grid
        allocatable :: keyg(:,:)
! wavefunction 
        allocatable :: psi(:,:),psit(:,:)
! wavefunction gradients
        allocatable :: hpsi(:,:),hpsit(:,:)

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
        allocatable :: psppar(:,:,:),nelpsp(:),radii_cf(:,:)
        include 'mpif.h'

        write(*,*) 'CLUSTER CLUSTER CLUSTER CLUSTER CLUSTER CLUSTER CLUSTER CLUSTER CLUSTER'

        open(unit=77,file='timings_comput')
        open(unit=78,file='timings_commun')
        open(unit=79,file='malloc')
        call cpu_time(tcpu1)

        open(unit=1,file='input.dat',status='old')
	read(1,*) hgrid
	read(1,*) crmult
	read(1,*) frmult
	read(1,*) cpmult
	read(1,*) fpmult
                  if (fpmult.gt.frmult) write(*,*) 'NONSENSE: fpmult > frmult'
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
        allocate(psppar(0:2,0:4,ntypes),nelpsp(ntypes),radii_cf(ntypes,2))
      do ityp=1,ntypes
        filename = 'psppar.'//atomnames(ityp)
!        if (iproc.eq.0) write(*,*) 'opening PSP file ',filename
        open(unit=11,file=filename,status='old')
	read(11,'(a30)') label
	read(11,*) radii_cf(ityp,1),radii_cf(ityp,2)
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
        if (iproc.eq.0) write(*,'(a,3(1x,e14.7))') 'atomic positions shifted by',-cxmin,-cymin,-czmin
	do iat=1,nat
        rxyz(1,iat)=rxyz(1,iat)-cxmin
        rxyz(2,iat)=rxyz(2,iat)-cymin
        rxyz(3,iat)=rxyz(3,iat)-czmin
        enddo

        if (iproc.eq.0) then
           write(*,*) 'Shifted atomic positions'
           do iat=1,nat
              write(*,'(i3,3(x,e14.7))') iat,(rxyz(j,iat),j=1,3)
           enddo
        endif

!    grid sizes n1,n2,n3
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

! fine grid size (needed for creation of input wavefunction, preconditioning)
      nfl1=n1 ; nfl2=n2 ; nfl3=n3
      nfu1=0 ; nfu2=0 ; nfu3=0
      do iat=1,nat
      rad=radii_cf(iatype(iat),2)*frmult
      nfl1=min(nfl1,int(onem+(rxyz(1,iat)-rad)/hgrid)) ; nfu1=max(nfu1,int((rxyz(1,iat)+rad)/hgrid))
      nfl2=min(nfl2,int(onem+(rxyz(2,iat)-rad)/hgrid)) ; nfu2=max(nfu2,int((rxyz(2,iat)+rad)/hgrid))
      nfl3=min(nfl3,int(onem+(rxyz(3,iat)-rad)/hgrid)) ; nfu3=max(nfu3,int((rxyz(3,iat)+rad)/hgrid))
      enddo
         if (iproc.eq.0) then
           write(*,*) 'nfl1,nfu1 ',nfl1,nfu1
           write(*,*) 'nfl2,nfu2 ',nfl2,nfu2
           write(*,*) 'nfl3,nfu3 ',nfl3,nfu3
         endif

! determine localization region for all Wanier orbitals, but do not yet fill the descriptor arrays
    allocate(logrid(0:n1,0:n2,0:n3))

! coarse grid quantities
        call fill_logrid(n1,n2,n3,0,n1,0,n2,0,n3,nat,ntypes,iatype,rxyz, & 
                         radii_cf(1,1),crmult,hgrid,logrid)
         if (iproc.eq.0) then
          open(unit=22,file='grid.ascii',status='unknown')
          write(22,*) nat
          write(22,*) alat1,' 0. ',alat2
          write(22,*) ' 0. ',' 0. ',alat3
          do iat=1,nat
           write(22,'(3(1x,e12.5),3x,a20)') rxyz(1,iat),rxyz(2,iat),rxyz(3,iat),atomnames(iatype(iat))
          enddo
 	  do i3=0,n3 ; do i2=0,n2 ; do i1=0,n1
           if (logrid(i1,i2,i3)) write(22,'(3(1x,e10.3),1x,a4)') i1*hgrid,i2*hgrid,i3*hgrid,'  g '
          enddo ; enddo ; enddo 
         endif
	 call num_segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid,nseg_c,nvctr_c)
        if (iproc.eq.0) write(*,*) 'orbitals have coarse segment, elements',nseg_c,nvctr_c

! fine grid quantities
        call fill_logrid(n1,n2,n3,0,n1,0,n2,0,n3,nat,ntypes,iatype,rxyz, & 
                         radii_cf(1,2),frmult,hgrid,logrid)
         if (iproc.eq.0) then
 	  do i3=0,n3 ; do i2=0,n2 ; do i1=0,n1
           if (logrid(i1,i2,i3)) write(22,'(3(1x,e10.3),1x,a4)') i1*hgrid,i2*hgrid,i3*hgrid,'  G '
          enddo ; enddo ; enddo 
         endif
	 call num_segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid,nseg_f,nvctr_f)
        if (iproc.eq.0) write(*,*) 'orbitals have fine   segment, elements',nseg_f,7*nvctr_f

        if (iproc.eq.0) close(22)

! allocations for arrays holding the wavefunctions and their data descriptors
        allocate(keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f))
        tt=dble(norb)/dble(nproc)
        norbp=int((1.d0-eps_mach*tt) + tt)
        write(79,'(a40,i10)') 'words for psi and hpsi ',2*(nvctr_c+7*nvctr_f)*norbp
        allocate(psi(nvctr_c+7*nvctr_f,norbp),hpsi(nvctr_c+7*nvctr_f,norbp))
        write(79,*) 'allocation done'
        norbme=max(min((iproc+1)*norbp,norb)-iproc*norbp,0)
        write(*,*) 'iproc ',iproc,' treats ',norbme,' orbitals '

        tt=dble(nvctr_c+7*nvctr_f)/dble(nproc)
        nvctrp=int((1.d0-eps_mach*tt) + tt)
        if (parallel) then
          write(79,'(a40,i10)') 'words for psit',nvctrp*norbp*nproc
          allocate(psit(nvctrp,norbp*nproc))
        endif

        if (iproc.eq.0) write(*,*) 'norbp,nvctrp=',norbp,nvctrp

! now fill the wavefunction descriptor arrays
! coarse grid quantities
        call fill_logrid(n1,n2,n3,0,n1,0,n2,0,n3,nat,ntypes,iatype,rxyz, & 
                         radii_cf(1,1),crmult,hgrid,logrid)
        call segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid,nseg_c,keyg(1,1),keyv(1))

! fine grid quantities
        call fill_logrid(n1,n2,n3,0,n1,0,n2,0,n3,nat,ntypes,iatype,rxyz, & 
                         radii_cf(1,2),frmult,hgrid,logrid)
        call segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid,nseg_f,keyg(1,nseg_c+1),keyv(nseg_c+1))

       call compresstest(iproc,n1,n2,n3,norb,norbp,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,hpsi)


! determine localization region for all projectors, but do not yet fill the descriptor arrays
    allocate(nboxp_c(2,3,nat),nboxp_f(2,3,nat))
    allocate(nseg_p(0:2*nat))
    allocate(nvctr_p(0:2*nat))
    nseg_p(0)=0 
    nvctr_p(0)=0 

    istart=1
    nproj=0
    do iat=1,nat

     if (iproc.eq.0) write(*,*) '+++++++++++++++++++++++++++++++++++++++++++ iat=',iat

    call numb_proj(iatype(iat),ntypes,psppar,mproj)
    if (mproj.ne.0) then 
    
    if (iproc.eq.0) write(*,*) 'projector descriptors for atom with mproj ',iat,mproj
    nproj=nproj+mproj

! coarse grid quantities
        call  pregion_size(rxyz(1,iat),radii_cf(1,2),cpmult,iatype(iat),ntypes, &
                   hgrid,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)
        if (iproc.eq.0) write(*,'(a,6(i4))') 'coarse grid',nl1,nu1,nl2,nu2,nl3,nu3
        nboxp_c(1,1,iat)=nl1 ; nboxp_c(2,1,iat)=nu1
        nboxp_c(1,2,iat)=nl2 ; nboxp_c(2,2,iat)=nu2
        nboxp_c(1,3,iat)=nl3 ; nboxp_c(2,3,iat)=nu3
        call fill_logrid(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,1,  &
                         ntypes,iatype(iat),rxyz(1,iat),radii_cf(1,2),cpmult,hgrid,logrid)
	 call num_segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,mvctr)
         if (iproc.eq.0) write(*,*) 'mseg,mvctr,coarse projectors ',mseg,mvctr
         nseg_p(2*iat-1)=nseg_p(2*iat-2) + mseg
         nvctr_p(2*iat-1)=nvctr_p(2*iat-2) + mvctr
         istart=istart+mvctr*mproj

! fine grid quantities
        call  pregion_size(rxyz(1,iat),radii_cf(1,2),fpmult,iatype(iat),ntypes, &
                   hgrid,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)
        if (iproc.eq.0) write(*,'(a,6(i4))') 'fine   grid',nl1,nu1,nl2,nu2,nl3,nu3
        nboxp_f(1,1,iat)=nl1 ; nboxp_f(2,1,iat)=nu1
        nboxp_f(1,2,iat)=nl2 ; nboxp_f(2,2,iat)=nu2
        nboxp_f(1,3,iat)=nl3 ; nboxp_f(2,3,iat)=nu3
        call fill_logrid(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,1,  &
                         ntypes,iatype(iat),rxyz(1,iat),radii_cf(1,2),fpmult,hgrid,logrid)
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
        write(79,'(a40,i10)') 'words for proj ',nprojel
        allocate(proj(nprojel))
        write(79,*) 'allocation done'
        

     if (iproc.eq.0) write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
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
                         ntypes,iatype(iat),rxyz(1,iat),radii_cf(1,2),cpmult,hgrid,logrid)

         iseg=nseg_p(2*iat-2)+1
         mseg=nseg_p(2*iat-1)-nseg_p(2*iat-2)
	 call segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
                           logrid,mseg,keyg_p(1,iseg),keyv_p(iseg))

! fine grid quantities
        nl1=nboxp_f(1,1,iat) ; nu1=nboxp_f(2,1,iat)
        nl2=nboxp_f(1,2,iat) ; nu2=nboxp_f(2,2,iat)
        nl3=nboxp_f(1,3,iat) ; nu3=nboxp_f(2,3,iat)
        call fill_logrid(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,1,  &
                         ntypes,iatype(iat),rxyz(1,iat),radii_cf(1,2),fpmult,hgrid,logrid)
         iseg=nseg_p(2*iat-1)+1
         mseg=nseg_p(2*iat)-nseg_p(2*iat-1)
	 call segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
                           logrid,mseg,keyg_p(1,iseg),keyv_p(iseg))

    endif
    enddo

     if (iproc.eq.0) write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

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
        write(79,'(a40,i10)') 'words for rhopot and pot_ion ',2*(2*n1+31)*(2*n2+31)*(2*n3+31)
       allocate(rhopot((2*n1+31)*(2*n2+31)*(2*n3+31)),pot_ion((2*n1+31)*(2*n2+31)*(2*n3+31)))
                 call zero((2*n1+31)*(2*n2+31)*(2*n3+31),pot_ion)
        write(79,*) 'allocation done'
! Allocate and calculate the 1/|r-r'| kernel for the solution of Poisson's equation and test it
       ndegree_ip=14
       if (parallel) then
          call calculate_pardimensions(2*n1+31,2*n2+31,2*n3+31,m1,m2,m3,nf1,nf2,nf3,md1,md2,md3,nfft1,nfft2,nfft3,nproc)
          !call Dimensions_FFT(2*n1+31,2*n2+31,2*n3+31,nfft1,nfft2,nfft3)
          write(*,*) 'dimension of FFT grid',nf1,nf2,nf3
          write(*,*) 'dimension of kernel',nfft1,nfft2,nfft3/nproc
          write(79,'(a40,i10)') 'words for kernel ',nfft1*nfft2*nfft3/nproc
          allocate(pkernel(nfft1*nfft2*nfft3/nproc))
           write(79,*) 'allocation done'
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
          call ParBuild_Kernel(2*n1+31,2*n2+31,2*n3+31,nf1,nf2,nf3,nfft1,nfft2,nfft3, &
               hgridh,ndegree_ip,iproc,nproc,pkernel)
          print *,"kernel built! iproc=",iproc
          call PARtest_kernel(2*n1+31,2*n2+31,2*n3+31,nfft1,nfft2,nfft3,hgridh,pkernel,pot_ion,rhopot,iproc,nproc) 
         
       else
          call Dimensions_FFT(2*n1+31,2*n2+31,2*n3+31,nfft1,nfft2,nfft3)
          write(*,*) 'dimension of FFT grid',nfft1,nfft2,nfft3
          write(*,*) 'dimension of kernel',nfft1/2+1,nfft2/2+1,nfft3/2+1
          write(79,'(a40,i10)') 'words for kernel ',(nfft1/2+1)*(nfft2/2+1)*(nfft3/2+1)
          allocate(pkernel((nfft1/2+1)*(nfft2/2+1)*(nfft3/2+1)))
           write(79,*) 'allocation done'
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
       if (calc_inp_wf .and. firstrun ) then
        firstrun=.false.
	call input_wf_diag(parallel,iproc,nproc,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
                    nat,norb,norbp,n1,n2,n3,nfft1,nfft2,nfft3,nvctr_c,nvctr_f,nvctrp,hgrid,rxyz, & 
                    rhopot,pot_ion,nseg_c,nseg_f,keyg,keyv, &
                    nprojel,nproj,nseg_p,keyg_p,keyv_p,nvctr_p,proj,  &
                    atomnames,ntypes,iatype,pkernel,psppar,psi,cprec,accurex)
        if (iproc.eq.0) then
        write(*,*) 'expected accuracy in total energy due to grid size',accurex
        write(*,*) 'suggested value for gnrm_cv ',accurex
        endif
        if (iproc.eq.0) write(*,*) 'input wavefunction has been calculated'
       else
        call readmywaves(iproc,norb,norbp,n1,n2,n3,hgrid,rxyz(1,1),nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,cprec)
        write(*,*) iproc,' readmywaves finished'
        write(*,*) iproc,'Read input wavefunctions from file'
        if (iproc.eq.0) write(*,*) 'Read input wavefunctions from file'
       endif


       if (parallel) then
       call transallwaves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,psi,psit)
       call loewe_p(iproc,nproc,norb,norbp,nvctrp,psit)
       call checkortho_p(iproc,nproc,norb,norbp,nvctrp,psit)
       call untransallwaves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,psit,psi)
       else
       call loewe(norb,norbp,nvctrp,psi)
       call checkortho(norb,norbp,nvctrp,psi)
       endif

       alpha=1.d0
       energy=1.d100
       gnrm=1.d100
! loop for wavefunction minimization
  loop_wavemin: do iter=1,itermax
        if (iproc.eq.0) then 
          write(*,*) '-------------------------------------- iter= ',iter
          if (gnrm.le.gnrm_cv) then
           write(*,'(a,i3,3(1x,e18.11))') 'iproc,ehart,eexcu,vexcu',iproc,ehart,eexcu,vexcu
           write(*,'(a,3(1x,e18.11))') 'final ekin_sum,epot_sum,eproj_sum',ekin_sum,epot_sum,eproj_sum
           write(*,'(a,3(1x,e18.11))') 'final ehart,eexcu,vexcu',ehart,eexcu,vexcu
           write(*,'(a,i6,2x,e19.12,1x,e9.2))') 'FINAL iter,total energy,gnrm',iter,energy,gnrm
         endif
       endif

! Potential from electronic charge density
       call sumrho(parallel,iproc,norb,norbp,n1,n2,n3,hgrid,occup,  & 
                   nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,rhopot)

       call cpu_time(tr0)
       call system_clock(ncount1,ncount_rate,ncount_max)

       if (parallel) then
          call ParPSolver_Kernel(2*n1+31,2*n2+31,2*n3+31,nfft1,nfft2,nfft3,hgridh,pkernel,.true., & 
               pot_ion,rhopot,ehart,eexcu,vexcu,iproc,nproc)
       else
          call PSolver_Kernel(2*n1+31,2*n2+31,2*n3+31,nfft1,nfft2,nfft3,hgridh,pkernel,.true., & 
               pot_ion,rhopot,ehart,eexcu,vexcu)
       end if

       call cpu_time(tr1)
       call system_clock(ncount2,ncount_rate,ncount_max)
       tel=dble(ncount2-ncount1)/dble(ncount_rate)
       write(77,'(a40,i4,2(x,e10.3))') 'PSOLVER TIME',iproc,tr1-tr0,tel


! local potential and kinetic energy for all orbitals belonging to iproc
        call applylocpotkin(iproc,norb,norbp,n1,n2,n3,hgrid,occup,  & 
                   nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,rhopot,hpsi,epot_sum,ekin_sum)

 
! apply all PSP projectors for all orbitals belonging to iproc
        call applyprojectors(iproc,ntypes,nat,iatype,psppar,occup, &
                    nprojel,nproj,nseg_p,keyg_p,keyv_p,nvctr_p,proj,  &
                    norb,norbp,nseg_c,nseg_f,keyg,keyv,nvctr_c,nvctr_f,psi,hpsi,eproj_sum)


       if (parallel) then
       wrkallred(1,2)=ekin_sum ; wrkallred(2,2)=epot_sum ; wrkallred(3,2)=eproj_sum
       call MPI_ALLREDUCE(wrkallred(1,2),wrkallred(1,1),3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
       ekin_sum=wrkallred(1,1) ; epot_sum=wrkallred(2,1) ; eproj_sum=wrkallred(3,1) 
       endif
       energybs=ekin_sum+epot_sum+eproj_sum
       energyold=energy
       energy=energybs-ehart+eexcu-vexcu+eion

!check for convergence or wheher max. numb. of iterations exceeded
       if (gnrm.le.gnrm_cv .or. iter.eq.itermax) then 
       if (iproc.eq.0) write(*,*) iter,' minimization iterations required'
          exit
       endif

! check whether CPU time exceeded
         open(unit=55,file='CPUlimit',status='unknown')
         read(55,*,end=555) cpulimit
         close(55)
         call cpu_time(tcpu2)
           if (tcpu2-tcpu1.ge.cpulimit) then
           write(6,*) iproc, ' CPU time exceeded'
           exit
         endif
555      continue
         close(55)



! Apply  orthogonality constraints to all orbitals belonging to iproc
      if (parallel) then
        allocate(hpsit(nvctrp,norbp*nproc))
        call transallwaves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,hpsi,hpsit)
        call  orthoconstraint_p(iproc,nproc,norb,norbp,occup,nvctrp,psit,hpsit,scprsum)
        call untransallwaves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,hpsit,hpsi)
        deallocate(hpsit)
      else
        call  orthoconstraint(norb,norbp,occup,nvctrp,psi,hpsi,scprsum)
      endif

! norm of residue
     gnrm=0.d0
     do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
       scpr=dnrm2(nvctr_c+7*nvctr_f,hpsi(1,iorb-iproc*norbp),1) 
       gnrm=gnrm+scpr**2
     enddo


       if (parallel) then
       tt=gnrm
       call MPI_ALLREDUCE(tt,gnrm,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
       endif
       gnrm=sqrt(gnrm)

! Preconditions all orbitals belonging to iproc
        call preconditionall(iproc,nproc,norb,norbp,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,hgrid,  & 
                   nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,hpsi,cprec)
      if (parallel) then
        allocate(hpsit(nvctrp,norbp*nproc))
        call transallwaves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,hpsi,hpsit)
      endif

! update wavefunction
       if (energy.gt.energyold) then
       alpha=max(.125d0,.5d0*alpha)
       if (alpha.eq..125d0) write(*,*) 'Convergence problem or limit'
       else
       alpha=min(1.05d0*alpha,2.d0)
       endif
       if (iproc.eq.0) write(*,*) 'alpha=',alpha

! update all wavefunctions (belonging to iproc)  with the preconditioned gradient
     if (parallel) then
        do iorb=1,norb
           call DAXPY(nvctrp,-alpha,hpsit(1,iorb),1,psit(1,iorb),1)
        enddo
        deallocate(hpsit)
     else
        do iorb=1,norb
           call DAXPY(nvctrp,-alpha,hpsi(1,iorb),1,psi(1,iorb),1)
        enddo
     endif

     if (parallel) then
       call loewe_p(iproc,nproc,norb,norbp,nvctrp,psit)
!       call checkortho_p(iproc,nproc,norb,norbp,nvctrp,psit)
    else
       call loewe(norb,norbp,nvctrp,psi)
!       call checkortho(norb,norbp,nvctrp,psi)
    endif


       tt=energybs-scprsum
       if (abs(tt).gt.1.d-10) then 
          write(*,*) 'ERROR: inconsistency between gradient and energy',tt,energybs,scprsum
       endif
       if (iproc.eq.0) then
       write(*,'(a,3(1x,e18.11))') 'ekin_sum,epot_sum,eproj_sum',  & 
                                    ekin_sum,epot_sum,eproj_sum
       write(*,'(a,3(1x,e18.11))') 'ehart,eexcu,vexcu',ehart,eexcu,vexcu
       write(*,'(a,i6,2x,e19.12,1x,e9.2))') 'iter,total energy,gnrm',iter,energy,gnrm
       endif

       write(77,*) '----------------------------------------------'
       write(78,*) '----------------------------------------------'

     if (parallel) then
        call untransallwaves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,psit,psi)
     endif

   enddo loop_wavemin


! transform to KS orbitals
     if (parallel) then
        allocate(hpsit(nvctrp,norbp*nproc))
        call transallwaves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,hpsi,hpsit)
        call KStrans_p(iproc,nproc,norb,norbp,nvctrp,occup,hpsit,psit,evsum,cprec)
        deallocate(hpsit)
        call untransallwaves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,psit,psi)
     else
	call KStrans(norb,norbp,nvctrp,occup,hpsi,psi,evsum,cprec)
     endif
       if (abs(evsum-energybs).gt.1.d-8) write(*,*) 'Difference:evsum,energybs',evsum,energybs

       call cpu_time(tcpu2)
       write(6,*) 'total CPU time without forces', tcpu2-tcpu1
       write(78,*) 'total CPU time without forces', tcpu2-tcpu1

! here we start the calculation of the forces
  if (iproc.eq.0) write(*,*)'calculation of forces'

! ground state electronic density
       call sumrho(parallel,iproc,norb,norbp,n1,n2,n3,hgrid,occup,  &
                   nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,rhopot)

! electronic potential
! for the moment let us use the pot_ion array for the electronic potential  
! and save rho in rhopot for calculation of fsep
    call DCOPY((2*n1+31)*(2*n2+31)*(2*n3+31),rhopot,1,pot_ion,1)
      if (parallel) then
          call ParPSolver_Kernel(2*n1+31,2*n2+31,2*n3+31,nfft1,nfft2,nfft3,hgridh,pkernel,.false., &
               rhopot,pot_ion,ehart_fake,eexcu_fake,vexcu_fake,iproc,nproc)
       else
          call PSolver_Kernel(2*n1+31,2*n2+31,2*n3+31,nfft1,nfft2,nfft3,hgridh,pkernel,.false., &
               rhopot,pot_ion,ehart_fake,eexcu_fake,vexcu_fake)
       end if

  if (iproc.eq.0) write(*,*)'electronic potential calculated'
  allocate(gxyz(3,nat))

! calculate local part of the forces gxyz
   call local_forces(iproc,nproc,ntypes,nat,iatype,atomnames,rxyz,psppar,nelpsp,hgrid,&
                     n1,n2,n3,rhopot,pot_ion,gxyz)

! Add the nonlocal part of the forces to gxyz
! calculating derivatives of the projectors (for the moment recalculate projectors)
  call nonlocal_forces(iproc,nproc,n1,n2,n3,nboxp_c,nboxp_f, &
     ntypes,nat,norb,norbp,istart,nprojel,nproj,&
     iatype,psppar,occup,nseg_c,nseg_f,nvctr_c,nvctr_f,nseg_p,nvctr_p,proj,  &
     keyg,keyv,keyg_p,keyv_p,psi,rxyz,radii_cf,cpmult,fpmult,hgrid,gxyz)

! Add up all the force contributions
  if (parallel) then
     call MPI_ALLREDUCE(gxyz,fxyz,3*nat,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  else
    do iat=1,nat
     fxyz(1,iat)=gxyz(1,iat)
     fxyz(2,iat)=gxyz(2,iat)
     fxyz(3,iat)=gxyz(3,iat)
    enddo
  end if


  deallocate(gxyz)

       call cpu_time(tcpu2)
       write(6,*) 'total CPU time with forces', tcpu2-tcpu1
       write(78,*) 'total CPU time with forces', tcpu2-tcpu1


!  write all the wavefunctions into files
      write(*,*) iproc,' start writing waves,cprec=',cprec
        call  writemywaves(iproc,norb,norbp,n1,n2,n3,hgrid,  & 
                   rxyz(1,1),nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,cprec)
      write(*,*) iproc,' finished writing waves'
        deallocate(pkernel)
        deallocate(keyg_p,keyv_p,proj)
        deallocate(rhopot,pot_ion)
        deallocate(occup)
        deallocate(nvctr_p,nseg_p)
        deallocate(psi,hpsi)
        deallocate(keyg,keyv)
        deallocate(psppar,nelpsp,radii_cf)

        if (parallel) call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        if (iproc.eq.0) call system("rm CPUlimit")

	end


