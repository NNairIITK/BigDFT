subroutine solveKS(parallel, iproc, nproc, norb, norbp, norbe, norbep, nvctr_c, &
  & nvctr_f, nvctrp, psi, hpsi, ppsi, cprec)

implicit real*8 (a-h,o-z)

logical, intent(in) :: parallel
integer, intent(in) :: nvctrp, norbe, norbep, nproc, iproc
integer, intent(in) :: nvctr_c, nvctr_f, norb, norbp
dimension psi(nvctr_c+7*nvctr_f,norbep), hpsi(nvctr_c+7*nvctr_f,norbep)
real*8, intent(out) :: cprec
real*8, intent(out) :: ppsi(nvctr_c+7*nvctr_f,norbp)

allocatable :: ppsit(:,:), psit(:,:), hpsit(:,:), hamovr(:,:,:),evale(:),work_lp(:)

        allocate(hamovr(norbe,norbe,4))
 if (parallel) then
        write(79,'(a40,i10)') 'words for psit inguess',nvctrp*norbep*nproc
        allocate(psit(nvctrp,norbep*nproc))
        write(79,*) 'allocation done'

        call  transallwaves(iproc,nproc,norbe,norbep,nvctr_c,nvctr_f,nvctrp,psi,psit)

        write(79,'(a40,i10)') 'words for hpsit inguess',2*nvctrp*norbep*nproc
        allocate(hpsit(nvctrp,norbep*nproc))
        write(79,*) 'allocation done'

        call  transallwaves(iproc,nproc,norbe,norbep,nvctr_c,nvctr_f,nvctrp,hpsi,hpsit)

!       hamovr(jorb,iorb,3)=+psit(k,jorb)*hpsit(k,iorb)
!       hamovr(jorb,iorb,4)=+psit(k,jorb)* psit(k,iorb)
      call DGEMM('T','N',norbe,norbe,nvctrp,1.d0,psit,nvctrp,hpsit,nvctrp,0.d0,hamovr(1,1,3),norbe)
      call DGEMM('T','N',norbe,norbe,nvctrp,1.d0,psit,nvctrp, psit,nvctrp,0.d0,hamovr(1,1,4),norbe)
        deallocate(hpsit)

        call MPI_ALLREDUCE (hamovr(1,1,3),hamovr(1,1,1),2*norbe**2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

! calculate  KS orbitals
!      if (iproc.eq.0) then
!        write(*,*) 'KS Hamiltonian',iproc
!        do iorb=1,norbe
!        write(*,'(10(1x,e10.3))') (hamovr(iorb,jorb,1),jorb=1,norbe)
!        enddo
!        write(*,*) 'Overlap',iproc
!        do iorb=1,norbe
!        write(*,'(10(1x,e10.3))') (hamovr(iorb,jorb,2),jorb=1,norbe)
!        enddo
!     endif

        n_lp=5000
        allocate(work_lp(n_lp),evale(norbe))
        call  DSYGV(1,'V','U',norbe,hamovr(1,1,1),norbe,hamovr(1,1,2),norbe,evale, work_lp, n_lp, info )
        if (info.ne.0) write(*,*) 'DSYGV ERROR',info
        if (iproc.eq.0) then
        do iorb=1,norbe
        write(*,*) 'evale(',iorb,')=',evale(iorb)
        enddo
        endif
        cprec=evale(1)
        deallocate(work_lp,evale)

        write(79,'(a40,i10)') 'words for ppsit ',nvctrp*norbep*nproc
        allocate(ppsit(nvctrp,norbp*nproc))
        write(79,*) 'allocation done'

! ppsit(k,iorb)=+psit(k,jorb)*hamovr(jorb,iorb,1)
      call DGEMM('N','N',nvctrp,norb,norbe,1.d0,psit,nvctrp,hamovr,norbe,0.d0,ppsit,nvctrp)

       call  untransallwaves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,ppsit,ppsi)

        deallocate(psit,ppsit)

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  else !serial case
!       hamovr(jorb,iorb,3)=+psi(k,jorb)*hpsi(k,iorb)
      call DGEMM('T','N',norbe,norbe,nvctrp,1.d0,psi,nvctrp,hpsi,nvctrp,0.d0,hamovr(1,1,1),norbe)
      call DGEMM('T','N',norbe,norbe,nvctrp,1.d0,psi,nvctrp, psi,nvctrp,0.d0,hamovr(1,1,2),norbe)

! calculate  KS orbitals
        write(*,*) 'KS Hamiltonian'
        do iorb=1,norbe
        write(*,'(10(1x,e10.3))') (hamovr(iorb,jorb,1),jorb=1,norbe)
        enddo
        write(*,*) 'Overlap'
        do iorb=1,norbe
        write(*,'(10(1x,e10.3))') (hamovr(iorb,jorb,2),jorb=1,norbe)
        enddo

        n_lp=5000
        allocate(work_lp(n_lp),evale(norbe))
        call  DSYGV(1,'V','U',norbe,hamovr(1,1,1),norbe,hamovr(1,1,2),norbe,evale, work_lp, n_lp, info )
        if (info.ne.0) write(*,*) 'DSYGV ERROR',info
        if (iproc.eq.0) then
        do iorb=1,norbe
        write(*,*) 'evale(',iorb,')=',evale(iorb)
        enddo
        endif
        cprec=evale(1)
        deallocate(work_lp,evale)

! ppsi(k,iorb)=+psi(k,jorb)*hamovr(jorb,iorb,1)
        call DGEMM('N','N',nvctrp,norb,norbe,1.d0,psi,nvctrp,hamovr,norbe,0.d0,ppsi,nvctrp)

  endif
  deallocate(hamovr)
end subroutine solveKS

subroutine readAtomicOrbitals(iproc, filename, ngx, npsp, xp, psiat, occupat, ng, &
                            & ns, np, pspatomnames, norbe, atomnames, iatype, nat)
  character(len = *), intent(in) :: filename
	integer, intent(in) :: ngx, npsp, iproc
	real*8, intent(out) :: xp(ngx, npsp), psiat(ngx, 2, npsp), occupat(2, npsp)
	integer, intent(out) :: ng(npsp), ns(npsp), np(npsp)
	character(len = 20), intent(out) :: pspatomnames(npsp)
	character(len = 20), intent(in) :: atomnames(100)
	integer, intent(out) :: norbe
	integer, intent(in) :: nat
	integer, intent(in) :: iatype(nat)
	
	integer :: ity, iat, i, j, ipsp
	
	! Read the data file.
	open(unit=24,file=trim(filename),form='formatted',status='unknown')
	do ity=1,15
33	format(30(e12.5))
19	format(a)
	read(24,19) pspatomnames(ity)
	read(24,*) ns(ity),(occupat(i,ity),i=1,ns(ity)),  &
                   np(ity),(occupat(i,ity),i=ns(ity)+1,ns(ity)+np(ity))
	if (ns(ity)+np(ity).gt.2) stop 'error ns+np'
	read(24,*) ng(ity)
        if (ng(ity).gt.ngx) stop 'enlarge ngx'
	read(24,33) (xp(i,ity)  ,i=1,ng(ity))
	do i=1,ng(ity) 
	read(24,*) (psiat(i,j,ity),j=1,ns(ity)+np(ity))
	enddo
	enddo
	close(unit=24)

  ! number of orbitals 
  norbe=0
	do iat=1,nat
	  ity=iatype(iat)
          do i=1,15
          if (pspatomnames(i).eq.atomnames(ity)) then
             ipsp=i
             goto 333
          endif
          enddo
          if (iproc.eq.0) write(*,*) 'no PSP for ',atomnames(ity)
          stop 
333       continue
          if (iproc.eq.0) write(*,*) 'found input wavefunction data for ',iat,atomnames(ity)
          norbe=norbe+ns(ipsp)+3*np(ipsp)
!        if (iproc.eq.0) write(*,*) 'ns(ipsp),np(ipsp)',ns(ipsp),np(ipsp)
        enddo
        if (iproc.eq.0) write(*,*) 'number of orbitals used in the construction of input guess ',norbe

end subroutine readAtomicOrbitals


subroutine createAtomicOrbitals(iproc, npsp, pspatomnames, atomnames, nat, rxyz, &
  & norbe, norbep, occupe, occupat, ngx, xp, psiat, psiatn, ng, ns, np, nvctr_c, nvctr_f, &
  & n1, n2, n3, hgrid, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3, nseg_c, nseg_f, &
  & keyg, keyv, iatype, psi, eks)

  integer, intent(in) :: npsp, nat, norbe, norbep, ngx
  integer, intent(in) :: nvctr_c, nvctr_f, n1, n2, n3, nseg_c
  integer, intent(in) :: nfl1, nfu1, nfl2, nfu2, nfl3, nfu3
  integer, intent(in) :: keyg(2, nseg_c + nseg_f), keyv(nseg_c + nseg_f)
  integer, intent(in) :: iatype(nat)
  real*8, intent(in) :: hgrid
  real*8, intent(out) :: eks
	character(len = 20), intent(in) :: pspatomnames(npsp)
	character(len = 20), intent(in) :: atomnames(100)
	integer, intent(in) :: ng(npsp), ns(npsp), np(npsp)
	real*8, intent(in) :: rxyz(3, nat)
	real*8, intent(in) :: xp(ngx, npsp), psiat(ngx, 2, npsp)
	real*8, intent(in) :: occupat(2, npsp), psiatn(ngx)
	real*8, intent(out) :: psi(nvctr_c + 7 * nvctr_f, norbep), occupe(norbe)

  integer :: iorb, jorb, iat, ity, ipsp, i
  real*8 :: rx, ry, rz, ek, scpr
  logical :: myorb


  eks=0.d0
	iorb=0
	do 100,iat=1,nat

    rx=rxyz(1,iat) ; ry=rxyz(2,iat) ; rz=rxyz(3,iat)

	  ity=iatype(iat)
    do i=1,15
      if (pspatomnames(i).eq.atomnames(ity)) then
        ipsp=i
        goto 444
      endif
    enddo
    if (iproc.eq.0) write(*,*) 'no PSP for ',atomnames(ity)
    stop 
444     continue

! H:
  	if (ns(ipsp).eq.1 .and. np(ipsp).eq.0) then
      iorb=iorb+1
      jorb=iorb-iproc*norbep
    	if (iorb.gt.norbe) stop 'transgpw occupe'
      occupe(iorb)=occupat(1,ipsp)
      call atomkin(0,ng(ipsp),xp(1,ipsp),psiat(1,1,ipsp),psiatn,ek) ; eks=eks+ek*occupe(iorb)
      call myorbital(myorb,iorb,norbe,iproc,nproc)
      if (myorb) then
        call crtonewave(n1,n2,n3,ng(ipsp),0,0,0,xp(1,ipsp),psiatn,rx,ry,rz,hgrid, & 
              0,n1,0,n2,0,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
              nseg_c,nvctr_c,keyg(1,1),keyv(1),nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), & 
              psi(1,jorb),psi(nvctr_c+1,jorb))
        call wnrm(nvctr_c,nvctr_f,psi(1,jorb),psi(nvctr_c+1,jorb),scpr) 
        write(*,*) 'ATOMIC INPUT ORBITAL iorb,norm',iorb,scpr ;  scpr=1.d0/sqrt(scpr)
        call wscal(nvctr_c,nvctr_f,scpr,psi(1,jorb),psi(nvctr_c+1,jorb))
      endif

! Li,Be , Na,Mg
  	else if (ns(ipsp).eq.2 .and. np(ipsp).eq.0) then
      iorb=iorb+1
      jorb=iorb-iproc*norbep
	    if (iorb+1.gt.norbe) stop 'transgpw occupe'
!    semicore s
	    occupe(iorb)=occupat(1,ipsp)
      call atomkin(0,ng(ipsp),xp(1,ipsp),psiat(1,1,ipsp),psiatn,ek) ; eks=eks+ek*occupe(iorb)
      call myorbital(myorb,iorb,norbe,iproc,nproc)
      if (myorb) then
        call crtonewave(n1,n2,n3,ng(ipsp),0,0,0,xp(1,ipsp),psiatn,rx,ry,rz,hgrid, & 
             0,n1,0,n2,0,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
             nseg_c,nvctr_c,keyg(1,1),keyv(1),nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), & 
             psi(1,jorb),psi(nvctr_c+1,jorb))
	      call wnrm(nvctr_c,nvctr_f,psi(1,jorb),psi(nvctr_c+1,jorb),scpr) 
        write(*,*) 'ATOMIC INPUT ORBITAL iorb,norm',iorb,scpr ;  scpr=1.d0/sqrt(scpr)
	      call wscal(nvctr_c,nvctr_f,scpr,psi(1,jorb),psi(nvctr_c+1,jorb))
      endif
!    valence s
      iorb=iorb+1
      jorb=iorb-iproc*norbep
	    occupe(iorb)=occupat(2,ipsp)
      call atomkin(0,ng(ipsp),xp(1,ipsp),psiat(1,1,ipsp),psiatn,ek) ; eks=eks+ek*occupe(iorb)
      call myorbital(myorb,iorb,norbe,iproc,nproc)
      if (myorb) then
        call crtonewave(n1,n2,n3,ng(ipsp),0,0,0,xp(1,ipsp),psiatn,rx,ry,rz,hgrid, & 
             0,n1,0,n2,0,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
             nseg_c,nvctr_c,keyg(1,1),keyv(1),nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), & 
             psi(1,jorb),psi(nvctr_c+1,jorb))
	      call wnrm(nvctr_c,nvctr_f,psi(1,jorb),psi(nvctr_c+1,jorb),scpr) ; scpr=1.d0/sqrt(scpr)
        write(*,*) 'ATOMIC INPUT ORBITAL iorb,norm',iorb,scpr 
	      call wscal(nvctr_c,nvctr_f,scpr,psi(1,jorb),psi(nvctr_c+1,jorb))
      endif


! B,C,N,O,F , Al,Si,P,S,Cl
	  else if (ns(ipsp).eq.1 .and. np(ipsp).eq.1) then
      iorb=iorb+1
      jorb=iorb-iproc*norbep
	    if (iorb+3.gt.norbe) stop 'transgpw occupe'
!   valence s
	    occupe(iorb)=occupat(1,ipsp)
      call atomkin(0,ng(ipsp),xp(1,ipsp),psiat(1,1,ipsp),psiatn,ek) ; eks=eks+ek*occupe(iorb)
      call myorbital(myorb,iorb,norbe,iproc,nproc)
      if (myorb) then
        call crtonewave(n1,n2,n3,ng(ipsp),0,0,0,xp(1,ipsp),psiatn,rx,ry,rz,hgrid, & 
             0,n1,0,n2,0,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
             nseg_c,nvctr_c,keyg(1,1),keyv(1),nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), & 
             psi(1,jorb),psi(nvctr_c+1,jorb))
	      call wnrm(nvctr_c,nvctr_f,psi(1,jorb),psi(nvctr_c+1,jorb),scpr) 
        write(*,*) 'ATOMIC INPUT ORBITAL iorb,norm',iorb,scpr ;  scpr=1.d0/sqrt(scpr)
	      call wscal(nvctr_c,nvctr_f,scpr,psi(1,jorb),psi(nvctr_c+1,jorb))
      endif
!   valence px
      iorb=iorb+1
      jorb=iorb-iproc*norbep
	    occupe(iorb)=occupat(2,ipsp)*(1.d0/3.d0)
      call atomkin(1,ng(ipsp),xp(1,ipsp),psiat(1,2,ipsp),psiatn,ek) ; eks=eks+3.d0*ek*occupe(iorb)
      call myorbital(myorb,iorb,norbe,iproc,nproc)
      if (myorb) then
        call crtonewave(n1,n2,n3,ng(ipsp),1,0,0,xp(1,ipsp),psiatn,rx,ry,rz,hgrid, & 
             0,n1,0,n2,0,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
             nseg_c,nvctr_c,keyg(1,1),keyv(1),nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), & 
             psi(1,jorb),psi(nvctr_c+1,jorb))
	      call wnrm(nvctr_c,nvctr_f,psi(1,jorb),psi(nvctr_c+1,jorb),scpr) 
        write(*,*) 'ATOMIC INPUT ORBITAL iorb,norm',iorb,scpr ;  scpr=1.d0/sqrt(scpr)
	      call wscal(nvctr_c,nvctr_f,scpr,psi(1,jorb),psi(nvctr_c+1,jorb))
      endif
!   valence py
      iorb=iorb+1
      jorb=iorb-iproc*norbep
	    occupe(iorb)=occupat(2,ipsp)*(1.d0/3.d0)
	    call myorbital(myorb,iorb,norbe,iproc,nproc)
      if (myorb) then
        call crtonewave(n1,n2,n3,ng(ipsp),0,1,0,xp(1,ipsp),psiatn,rx,ry,rz,hgrid, & 
             0,n1,0,n2,0,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
             nseg_c,nvctr_c,keyg(1,1),keyv(1),nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), & 
             psi(1,jorb),psi(nvctr_c+1,jorb))
	      call wnrm(nvctr_c,nvctr_f,psi(1,jorb),psi(nvctr_c+1,jorb),scpr) 
        write(*,*) 'ATOMIC INPUT ORBITAL iorb,norm',iorb,scpr ;  scpr=1.d0/sqrt(scpr)
	      call wscal(nvctr_c,nvctr_f,scpr,psi(1,jorb),psi(nvctr_c+1,jorb))
      endif
!   valence pz
      iorb=iorb+1
      jorb=iorb-iproc*norbep
	    occupe(iorb)=occupat(2,ipsp)*(1.d0/3.d0)
	    call myorbital(myorb,iorb,norbe,iproc,nproc)
      if (myorb) then
        call crtonewave(n1,n2,n3,ng(ipsp),0,0,1,xp(1,ipsp),psiatn,rx,ry,rz,hgrid, & 
             0,n1,0,n2,0,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
             nseg_c,nvctr_c,keyg(1,1),keyv(1),nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), & 
             psi(1,jorb),psi(nvctr_c+1,jorb))
	      call wnrm(nvctr_c,nvctr_f,psi(1,jorb),psi(nvctr_c+1,jorb),scpr) 
        write(*,*) 'ATOMIC INPUT ORBITAL iorb,norm',iorb,scpr ;  scpr=1.d0/sqrt(scpr)
	      call wscal(nvctr_c,nvctr_f,scpr,psi(1,jorb),psi(nvctr_c+1,jorb))
      endif

! ERROR
	  else 
	    stop 'error in inguess.dat'
  	endif
100	continue

end subroutine createAtomicOrbitals


        subroutine transallwaves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,psi,psit)
        implicit real*8 (a-h,o-z)
        integer recvcount,sendcount
        dimension psi(nvctr_c+7*nvctr_f,norbp),psit(nvctrp,norbp*nproc)
        real*8, allocatable :: psiw(:,:,:)
        include 'mpif.h'
 
         allocate(psiw(nvctrp,norbp,nproc))
 
          sendcount=nvctrp*norbp
          recvcount=nvctrp*norbp
 
        call cpu_time(tr0) ; call system_clock(ncount0,ncount_rate,ncount_max)

! reformatting: psiw(i,iorb,j,jorb) <- psi(ij,iorb,jorb)
      do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
         ij=1
         do j=1,nproc
         do i=1,nvctrp
         if (ij .le. nvctr_c+7*nvctr_f) then
         psiw(i,iorb-iproc*norbp,j)=psi(ij,iorb-iproc*norbp)
         else
         psiw(i,iorb-iproc*norbp,j)=0.d0
         endif
         ij=ij+1
         enddo
         enddo
      enddo

        call cpu_time(tr1) ; call system_clock(ncount1,ncount_rate,ncount_max)
 
! transposition: psit(i,iorb,jorb,j) <- psiw(i,iorb,j,jorb) 
       call MPI_ALLTOALL(psiw,sendcount,MPI_DOUBLE_PRECISION,  &
                         psit,recvcount,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
 
        call cpu_time(tr2) ; call system_clock(ncount2,ncount_rate,ncount_max)
        tel=dble(ncount1-ncount0)/dble(ncount_rate)
        write(77,'(a40,i4,2(x,e10.3))') 'TRANSALL TIME',iproc,tr1-tr0,tel
        tel=dble(ncount2-ncount1)/dble(ncount_rate)
        write(78,'(a40,i4,2(x,e10.3))') 'TRANSALL TIME',iproc,tr2-tr1,tel
 
      deallocate(psiw)
 
        return
        end subroutine transallwaves


	
        subroutine untransallwaves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,psit,psi)
        implicit real*8 (a-h,o-z)
        integer recvcount,sendcount
        dimension psi(nvctr_c+7*nvctr_f,norbp),psit(nvctrp,norbp*nproc)
        real*8, allocatable :: psiw(:,:,:)
        include 'mpif.h'
 
         allocate(psiw(nvctrp,norbp,nproc))
 
          sendcount=nvctrp*norbp
          recvcount=nvctrp*norbp
 
! transposition: psiw(i,iorb,j,jorb) <- psit(i,iorb,jorb,j) 
       call MPI_ALLTOALL(psit,sendcount,MPI_DOUBLE_PRECISION,  &
                         psiw,recvcount,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
 
! reformatting: psi(ij,iorb,jorb) <- psiw(i,iorb,j,jorb)
      do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
         ij=1
         do j=1,nproc
         do i=1,nvctrp
         psi(ij,iorb-iproc*norbp)=psiw(i,iorb-iproc*norbp,j)
         ij=ij+1
         if (ij.gt. nvctr_c+7*nvctr_f) goto 333
         enddo
         enddo
333     continue
      enddo
 
        deallocate(psiw)
 
        return
        end subroutine untransallwaves
 


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


        subroutine input_rho_ion(iproc,ntypes,nat,iatype,atomnames,rxyz,psppar,nelpsp,n1,n2,n3,hgrid,rho,eion)
! Creates charge density arising from the ionoc PSP cores
        implicit real*8 (a-h,o-z)
        character*20 :: atomnames(100)
        dimension psppar(0:2,0:4,ntypes),rxyz(3,nat),iatype(nat),nelpsp(ntypes)
        dimension rho(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16)

        hgridh=hgrid*.5d0 
        pi=4.d0*atan(1.d0)
	call zero((2*n1+31)*(2*n2+31)*(2*n3+31),rho)

! Ionic charge 
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
      end subroutine input_rho_ion


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
      end subroutine addlocgauspsp
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


        subroutine compresstest(iproc,n1,n2,n3,norb,norbp,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,hpsi)
! test compress-uncompress
        implicit real*8 (a-h,o-z)
        dimension psi(nvctr_c+7*nvctr_f,norbp),hpsi(nvctr_c+7*nvctr_f,norbp)
        dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
        allocatable :: psifscf(:)


! begin test compress-uncompress
       do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)

        allocate(psifscf((2*n1+16)*(2*n2+16)*(2*n3+16)) )

	do i=1,nvctr_c+7*nvctr_f
        call random_number(tt)
        psi(i,iorb-iproc*norbp)=tt
        hpsi(i,iorb-iproc*norbp)=tt
        enddo

        call uncompress(n1,n2,n3,nseg_c,nvctr_c,keyg(1,1),keyv(1),   & 
                    nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   & 
                    psi(1,iorb-iproc*norbp),psi(nvctr_c+1,iorb-iproc*norbp),psifscf)
        call   compress(n1,n2,n3,nseg_c,nvctr_c,keyg(1,1),keyv(1),   & 
                    nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   & 
                    psifscf,psi(1,iorb-iproc*norbp),psi(nvctr_c+1,iorb-iproc*norbp))

        tc=0.d0
	do i=1,nvctr_c
        tc=max(tc,abs(psi(i,iorb-iproc*norbp)-hpsi(i,iorb-iproc*norbp)))
        enddo
        if (tc.gt.1.d-10) stop 'coarse compress error'

        tf=0.d0
	do i=1,7*nvctr_f
        tf=max(tf,abs(psi(nvctr_c+i,iorb-iproc*norbp)-hpsi(nvctr_c+i,iorb-iproc*norbp))) 
        enddo
        write(*,'(a,i4,2(1x,e9.2))') 'COMPRESSION TEST: iorb, coarse, fine error:',iorb,tc,tf
        if (tf.gt.1.d-10) stop 'fine compress error'

        deallocate(psifscf)

     enddo

	return
	end subroutine compresstest


        subroutine compress(n1,n2,n3,nseg_c,nvctr_c,keyg_c,keyv_c,  & 
                            nseg_f,nvctr_f,keyg_f,keyv_f,  & 
                            psifscf,psi_c,psi_f)
! Compresses a wavefunction that is given in terms of fine scaling functions (psifscf) into 
! the retained coarse scaling functions and wavelet coefficients (psi_c,psi_f)
        implicit real*8 (a-h,o-z)
        dimension keyg_c(2,nseg_c),keyv_c(nseg_c),keyg_f(2,nseg_f),keyv_f(nseg_f)
        dimension psi_c(nvctr_c),psi_f(7,nvctr_f)
        dimension psifscf((2*n1+16)*(2*n2+16)*(2*n3+16))
        real*8, allocatable :: psig(:,:,:,:,:,:),ww(:)

        allocate(psig(0:n1,2,0:n2,2,0:n3,2),  & 
                 ww((2*n1+16)*(2*n2+16)*(2*n3+16)))
! decompose wavelets into coarse scaling functions and wavelets
	call analyse_shrink(n1,n2,n3,ww,psifscf,psig)

! coarse part
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
            psi_c(i-i0+jj)=psig(i,1,i2,1,i3,1)
          enddo
        enddo

! fine part
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

      end subroutine compress


        subroutine uncompress(n1,n2,n3,nseg_c,nvctr_c,keyg_c,keyv_c,  & 
                              nseg_f,nvctr_f,keyg_f,keyv_f,  & 
                              psi_c,psi_f,psifscf)
! Expands the compressed wavefunction in vector form (psi_c,psi_f) 
! into fine scaling functions (psifscf)
        implicit real*8 (a-h,o-z)
        dimension keyg_c(2,nseg_c),keyv_c(nseg_c),keyg_f(2,nseg_f),keyv_f(nseg_f)
        dimension psi_c(nvctr_c),psi_f(7,nvctr_f)
        dimension psifscf((2*n1+16)*(2*n2+16)*(2*n3+16))
        real*8, allocatable :: psig(:,:,:,:,:,:),ww(:)

        allocate(psig(0:n1,2,0:n2,2,0:n3,2),  & 
                 ww((2*n1+16)*(2*n2+16)*(2*n3+16)))

        call zero(8*(n1+1)*(n2+1)*(n3+1),psig)

! coarse part
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
            psig(i,1,i2,1,i3,1)=psi_c(i-i0+jj)
          enddo
         enddo

! fine part
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
	call synthese_grow(n1,n2,n3,ww,psig,psifscf)

        deallocate(psig,ww)

      end subroutine uncompress


        subroutine applylocpotkin(iproc,norb,norbp,n1,n2,n3,hgrid,occup,  & 
                   nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,pot,hpsi,epot_sum,ekin_sum)
!  Applies the local potential and kinetic energy operator to all wavefunctions
! Input: pot,psi
! Output: hpsi,epot,ekin
        implicit real*8 (a-h,o-z)
        dimension occup(norb),pot((2*n1+31)*(2*n2+31)*(2*n3+31))
        dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
        dimension psi(nvctr_c+7*nvctr_f,norbp),hpsi(nvctr_c+7*nvctr_f,norbp)
        real*8, allocatable, dimension(:) :: psifscf,psifscfk,psir
        include 'mpif.h'

       call cpu_time(tr0)
       call system_clock(ncount1,ncount_rate,ncount_max)


        hgridh=hgrid*.5d0

! Wavefunction expressed everywhere in fine scaling functions (for potential and kinetic energy)
        allocate(psifscf((2*n1+16)*(2*n2+16)*(2*n3+16)) )
        allocate(psifscfk((2*n1+16)*(2*n2+16)*(2*n3+16)) )
! Wavefunction in real space
        allocate(psir((2*n1+31)*(2*n2+31)*(2*n3+31)))


        ekin_sum=0.d0
        epot_sum=0.d0
     do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)

        call uncompress(n1,n2,n3,nseg_c,nvctr_c,keyg(1,1),keyv(1),   &
                    nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
                    psi(1,iorb-iproc*norbp),psi(nvctr_c+1,iorb-iproc*norbp),psifscf)

        call convolut_kinetic(2*n1+15,2*n2+15,2*n3+15,hgridh,psifscf,psifscfk)

        ekin=0.d0 
        do i=1,(2*n1+16)*(2*n2+16)*(2*n3+16)
          ekin=ekin+psifscf(i)*psifscfk(i)
        enddo
        ekin_sum=ekin_sum+occup(iorb)*ekin 

        call convolut_magic_n(2*n1+15,2*n2+15,2*n3+15,psifscf,psir) 
        epot=0.d0
        do i=1,(2*n1+31)*(2*n2+31)*(2*n3+31)
          tt=pot(i)*psir(i)
          epot=epot+tt*psir(i)
          psir(i)=tt
        enddo
        epot_sum=epot_sum+occup(iorb)*epot
        write(*,'(a,i5,2(1x,e17.10))') 'iorb,ekin,epot',iorb,ekin,epot

        call convolut_magic_t(2*n1+15,2*n2+15,2*n3+15,psir,psifscf)

        do i=1,(2*n1+16)*(2*n2+16)*(2*n3+16)
          psifscf(i)=psifscf(i)+psifscfk(i)
        enddo

        call   compress(n1,n2,n3,nseg_c,nvctr_c,keyg(1,1),keyv(1),   & 
                    nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   & 
                    psifscf,hpsi(1,iorb-iproc*norbp),hpsi(nvctr_c+1,iorb-iproc*norbp))
     enddo

        deallocate(psifscf,psifscfk,psir)

       call cpu_time(tr1)
       call system_clock(ncount2,ncount_rate,ncount_max)
       tel=dble(ncount2-ncount1)/dble(ncount_rate)
       write(77,'(a40,i4,2(x,e10.3))') 'APPLYLOCPOTKIN TIME',iproc,tr1-tr0,tel

   
        return
      end subroutine applylocpotkin


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
      end subroutine convolut_magic_n


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
      end subroutine convolut_magic_t


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
      end subroutine synthese_grow



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
      end subroutine analyse_shrink




!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


        subroutine sumrho(parallel,iproc,norb,norbp,n1,n2,n3,hgrid,occup,  & 
                              nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,rho)
! Calculates the charge density by summing the square of all orbitals
! Input: psi
! Output: rho
        implicit real*8 (a-h,o-z)
        logical parallel
	dimension rho((2*n1+31)*(2*n2+31)*(2*n3+31)),occup(norb)
        dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
        dimension psi(nvctr_c+7*nvctr_f,norbp)
        real*8, allocatable, dimension(:) :: psifscf,psir,rho_p
        include 'mpif.h'

       call cpu_time(tr0)
       call system_clock(ncount1,ncount_rate,ncount_max)


        hgridh=hgrid*.5d0 

! Wavefunction expressed everywhere in fine scaling functions (for potential and kinetic energy)
        allocate(psifscf((2*n1+16)*(2*n2+16)*(2*n3+16)) )
! Wavefunction in real space
        allocate(psir((2*n1+31)*(2*n2+31)*(2*n3+31)))

 if (parallel) then
        allocate(rho_p((2*n1+31)*(2*n2+31)*(2*n3+31)))
	call zero((2*n1+31)*(2*n2+31)*(2*n3+31),rho_p)

      do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)

        call uncompress(n1,n2,n3,nseg_c,nvctr_c,keyg(1,1),keyv(1),   & 
                    nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   & 
                    psi(1,iorb-iproc*norbp),psi(nvctr_c+1,iorb-iproc*norbp),psifscf)

        call convolut_magic_n(2*n1+15,2*n2+15,2*n3+15,psifscf,psir) 

       do i=1,(2*n1+31)*(2*n2+31)*(2*n3+31)
         rho_p(i)=rho_p(i)+(occup(iorb)/hgridh**3)*psir(i)**2
        enddo

      enddo

        call cpu_time(tr0)
        call system_clock(ncount1,ncount_rate,ncount_max)
        call MPI_ALLREDUCE(rho_p,rho,(2*n1+31)*(2*n2+31)*(2*n3+31),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call cpu_time(tr1)
        call system_clock(ncount2,ncount_rate,ncount_max)
        tel=dble(ncount2-ncount1)/dble(ncount_rate)
        write(78,'(a40,i4,2(x,e10.3))') 'RHO: ALLREDUCE TIME',iproc,tr1-tr0,tel

        deallocate(rho_p)
 else

	call zero((2*n1+31)*(2*n2+31)*(2*n3+31),rho)

     do iorb=1,norb

        call uncompress(n1,n2,n3,nseg_c,nvctr_c,keyg(1,1),keyv(1),   & 
                    nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   & 
                    psi(1,iorb),psi(nvctr_c+1,iorb),psifscf)

        call convolut_magic_n(2*n1+15,2*n2+15,2*n3+15,psifscf,psir) 

       do i=1,(2*n1+31)*(2*n2+31)*(2*n3+31)
         rho(i)=rho(i)+(occup(iorb)/hgridh**3)*psir(i)**2
        enddo

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

       call cpu_time(tr1)
       call system_clock(ncount2,ncount_rate,ncount_max)
       tel=dble(ncount2-ncount1)/dble(ncount_rate)
       write(77,'(a40,i4,2(x,e10.3))') 'SUMRHO TIME',iproc,tr1-tr0,tel

        return
      end subroutine sumrho

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
	

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
      end subroutine numb_proj


        subroutine applyprojectors(iproc,ntypes,nat,iatype,psppar,occup, &
                    nprojel,nproj,nseg_p,keyg_p,keyv_p,nvctr_p,proj,  &
                    norb,norbp,nseg_c,nseg_f,keyg,keyv,nvctr_c,nvctr_f,psi,hpsi,eproj_sum)
! Applies all the projectors onto a wavefunction
! Input: psi_c,psi_f
! In/Output: hpsi_c,hpsi_f (both are updated, i.e. not initilized to zero at the beginning)
        implicit real*8 (a-h,o-z)
        dimension psppar(0:2,0:4,ntypes),iatype(nat)
        dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
        dimension psi(nvctr_c+7*nvctr_f,norbp),hpsi(nvctr_c+7*nvctr_f,norbp)
        dimension nseg_p(0:2*nat),nvctr_p(0:2*nat)
        dimension keyg_p(2,nseg_p(2*nat)),keyv_p(nseg_p(2*nat))
        dimension proj(nprojel),occup(norb)

       call cpu_time(tr0)
       call system_clock(ncount1,ncount_rate,ncount_max)


  eproj_sum=0.d0
! loop over all my orbitals
  do 100,iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)

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
                   nvctr_c,nvctr_f,nseg_c,nseg_f,keyv(1),keyv(nseg_c+1),  &
                   keyg(1,1),keyg(1,nseg_c+1),psi(1,iorb-iproc*norbp),psi(nvctr_c+1,iorb-iproc*norbp),  &
                   mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p(jseg_c),keyv_p(jseg_f),  &
                   keyg_p(1,jseg_c),keyg_p(1,jseg_f),proj(istart_c),proj(istart_f),scpr)
! test (will sometimes give wrong result)
        call wdot(  &
                   mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p(jseg_c),keyv_p(jseg_f),  &
                   keyg_p(1,jseg_c),keyg_p(1,jseg_f),proj(istart_c),proj(istart_f),  &
                   nvctr_c,nvctr_f,nseg_c,nseg_f,keyv(1),keyv(nseg_c+1),  &
                   keyg(1,1),keyg(1,nseg_c+1),psi(1,iorb-iproc*norbp),psi(nvctr_c+1,iorb-iproc*norbp),tcpr)
        if (scpr.ne.tcpr) stop 'projectors: scpr.ne.tcpr'
! testend

        scprp=scpr*psppar(i,j,ityp)
!        write(*,*) 'pdot scpr',scpr,psppar(i,j,ityp)
        eproj=eproj+scprp*scpr

        call waxpy(  &
             scprp,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p(jseg_c),keyv_p(jseg_f),  &
                   keyg_p(1,jseg_c),keyg_p(1,jseg_f),proj(istart_c),proj(istart_f),  &
                   nvctr_c,nvctr_f,nseg_c,nseg_f,keyv(1),keyv(nseg_c+1),  &
                   keyg(1,1),keyg(1,nseg_c+1),hpsi(1,iorb-iproc*norbp),hpsi(nvctr_c+1,iorb-iproc*norbp))

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

       call cpu_time(tr1)
       call system_clock(ncount2,ncount_rate,ncount_max)
       tel=dble(ncount2-ncount1)/dble(ncount_rate)
       write(77,'(a40,i4,2(x,e10.3))') 'APPLYPROJECTOR TIME',iproc,tr1-tr0,tel


         return
      end subroutine applyprojectors


        subroutine crtproj(iproc,nterm,n1,n2,n3, & 
                     nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,  & 
                     radius_f,cpmult,fpmult,hgrid,gau_a,factor,rx,ry,rz,lx,ly,lz, & 
                     mvctr_c,mvctr_f,proj_c,proj_f)
! returns the compressed form of a Gaussian projector 
! x^lx * y^ly * z^lz * exp (-1/(2*gau_a^2) *((x-cntrx)^2 + (y-cntry)^2 + (z-cntrz)^2 ))
! in the arrays proj_c, proj_f
        implicit real*8 (a-h,o-z)
        parameter(ntermx=3,nw=16000)
        dimension lx(3),ly(3),lz(3)
        dimension proj_c(mvctr_c),proj_f(7,mvctr_f)
        real*8, allocatable, dimension(:,:,:) :: wprojx, wprojy, wprojz
        real*8, allocatable, dimension(:,:) :: work

        allocate(wprojx(0:n1,2,ntermx),wprojy(0:n2,2,ntermx),wprojz(0:n3,2,ntermx),work(0:nw,2))


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
!        if (iproc.eq.0) write(*,*) 'max err_norm ',err_norm

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
          if (mvctr.ne.mvctr_c) stop 'mvctr >< mvctr_c'

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
     end subroutine crtproj
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

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
      end subroutine wdot



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
      end subroutine waxpy



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
      end subroutine wnrm



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
      end subroutine wscal


	subroutine wzero(mvctr_c,mvctr_f,psi_c,psi_f)
! initializes a wavefunction to zero
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
      end subroutine wzero


        subroutine orthoconstraint_p(iproc,nproc,norb,norbp,occup,nvctrp,psit,hpsit,scprsum)
!Effect of orthogonality constraints on gradient 
        implicit real*8 (a-h,o-z)
        dimension psit(nvctrp,norbp*nproc),hpsit(nvctrp,norbp*nproc),occup(norb)
        allocatable :: alag(:,:,:)
        include 'mpif.h'

       call cpu_time(tr0)
       call system_clock(ncount0,ncount_rate,ncount_max)

      allocate(alag(norb,norb,2))
!     alag(jorb,iorb,2)=+psit(k,jorb)*hpsit(k,iorb)
      call DGEMM('T','N',norb,norb,nvctrp,1.d0,psit,nvctrp,hpsit,nvctrp,0.d0,alag(1,1,2),norb)

       call cpu_time(tr1)
       call system_clock(ncount1,ncount_rate,ncount_max)
     call MPI_ALLREDUCE(alag(1,1,2),alag(1,1,1),norb**2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
       call cpu_time(tr2)
       call system_clock(ncount2,ncount_rate,ncount_max)

!        if (iproc.eq.0) then
!        write(*,*) 'ALAG',iproc
!        do iorb=1,norb
!        write(*,'(10(1x,e10.3))') (alag(iorb,jorb,1),jorb=1,norb)
!        enddo
!        endif

     scprsum=0.d0
     do iorb=1,norb
         scprsum=scprsum+occup(iorb)*alag(iorb,iorb,1)
     enddo

! hpsit(k,iorb)=-psit(k,jorb)*alag(jorb,iorb,1)
      call DGEMM('N','N',nvctrp,norb,norb,-1.d0,psit,nvctrp,alag,norb,1.d0,hpsit,nvctrp)
     deallocate(alag)
 
       call cpu_time(tr3)
       call system_clock(ncount3,ncount_rate,ncount_max)

       tel=dble(ncount1-ncount0+ncount3-ncount2)/dble(ncount_rate)
       write(77,'(a40,i4,2(x,e10.3))') 'ORTHOCONSTRAINT TIME',iproc,tr1-tr0+tr3-tr2,tel
       tel=dble(ncount2-ncount1)/dble(ncount_rate)
       write(78,'(a40,i4,2(x,e10.3))') 'ORTHOCONSTRAINT TIME',iproc,tr2-tr1,tel
 

       call cpu_time(tr1)
       call system_clock(ncount2,ncount_rate,ncount_max)


     return
   end subroutine orthoconstraint_p
     


        subroutine orthoconstraint(norb,norbp,occup,nvctrp,psi,hpsi,scprsum)
!Effect of orthogonality constraints on gradient 
        implicit real*8 (a-h,o-z)
        dimension psi(nvctrp,norbp),hpsi(nvctrp,norbp),occup(norb)
        allocatable :: alag(:,:,:)

       call cpu_time(tr0)
       call system_clock(ncount0,ncount_rate,ncount_max)

     allocate(alag(norb,norb,2))
 
!     alag(jorb,iorb,2)=+psi(k,jorb)*hpsi(k,iorb)
      call DGEMM('T','N',norb,norb,nvctrp,1.d0,psi,nvctrp,hpsi,nvctrp,0.d0,alag(1,1,1),norb)

     scprsum=0.d0
     do iorb=1,norb
         scprsum=scprsum+occup(iorb)*alag(iorb,iorb,1)
     enddo

! hpsit(k,iorb)=-psit(k,jorb)*alag(jorb,iorb,1)
      call DGEMM('N','N',nvctrp,norb,norb,-1.d0,psi,nvctrp,alag,norb,1.d0,hpsi,nvctrp)

     deallocate(alag)
       call cpu_time(tr1)
       call system_clock(ncount1,ncount_rate,ncount_max)
       tel=dble(ncount1-ncount0)/dble(ncount_rate)
       write(77,'(a40,i4,2(x,e10.3))') 'ORTHOCONSTRAINT TIME',0,tr1-tr0,tel

     return
   end subroutine orthoconstraint



!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

        subroutine system_size(nat,rxyz,radii,rmult,iatype,ntypes, &
                   cxmin,cxmax,cymin,cymax,czmin,czmax)
! calculates the overall size of the simulation cell (cxmin,cxmax,cymin,cymax,czmin,czmax)
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
      end subroutine system_size


        subroutine pregion_size(rxyz,radii,rmult,iatype,ntypes, &
                   hgrid,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)
! finds the size of the smallest subbox that contains a localization region made 
! out of atom centered spheres
        implicit real*8 (a-h,o-z)
        parameter(eps_mach=1.d-12)
        dimension rxyz(3),radii(ntypes)

            rad=radii(iatype)*rmult
            cxmax=rxyz(1)+rad ; cxmin=rxyz(1)-rad
            cymax=rxyz(2)+rad ; cymin=rxyz(2)-rad
            czmax=rxyz(3)+rad ; czmin=rxyz(3)-rad
!        write(*,*) radii(iatype(iat)),rmult
!        write(*,*) rxyz(1),rxyz(2),rxyz(3)
  
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
!        write(*,'(a,6(i4))') 'projector region ',nl1,nu1,nl2,nu2,nl3,nu3
!        write(*,*) ' projector region size ',cxmin,cxmax
!        write(*,*) '             ',cymin,cymax
!        write(*,*) '             ',czmin,czmax
      if (nl1.lt.0)   stop 'nl1: projector region outside cell'
      if (nl2.lt.0)   stop 'nl2: projector region outside cell'
      if (nl3.lt.0)   stop 'nl3: projector region outside cell'
      if (nu1.gt.n1)   stop 'nu1: projector region outside cell'
      if (nu2.gt.n2)   stop 'nu2: projector region outside cell'
      if (nu3.gt.n3)   stop 'nu3: projector region outside cell'

        return
      end subroutine pregion_size


	subroutine num_segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,mvctr)
! Calculates the length of the keys describing a wavefunction data structure
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
      end subroutine num_segkeys


	subroutine segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,keyg,keyv)
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
      end subroutine segkeys




       subroutine fill_logrid(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,nat,  &
                               ntypes,iatype,rxyz,radii,rmult,hgrid,logrid)
! set up an array logrid(i1,i2,i3) that specifies whether the grid point
! i1,i2,i3 is the center of a scaling function/wavelet
        implicit real*8 (a-h,o-z)
        logical logrid
        parameter(eps_mach=1.d-12,onem=1.d0-eps_mach)
        dimension rxyz(3,nat),iatype(nat),radii(ntypes)
        dimension logrid(0:n1,0:n2,0:n3)

        do i3=nl3,nu3 ; do i2=nl2,nu2 ; do i1=nl1,nu1
         logrid(i1,i2,i3)=.false.
        enddo ; enddo ; enddo

      do iat=1,nat
        rad=radii(iatype(iat))*rmult
!        write(*,*) 'iat,nat,rad',iat,nat,rad
        ml1=int(onem+(rxyz(1,iat)-rad)/hgrid)  ; mu1=int((rxyz(1,iat)+rad)/hgrid)
        ml2=int(onem+(rxyz(2,iat)-rad)/hgrid)  ; mu2=int((rxyz(2,iat)+rad)/hgrid)
        ml3=int(onem+(rxyz(3,iat)-rad)/hgrid)  ; mu3=int((rxyz(3,iat)+rad)/hgrid)
        if (ml1.lt.nl1) stop 'ml1 < nl1' ; if (mu1.gt.nu1) stop 'mu1 > nu1'
        if (ml2.lt.nl2) stop 'ml2 < nl2' ; if (mu2.gt.nu2) stop 'mu2 > nu2'
        if (ml3.lt.nl3) stop 'ml3 < nl3' ; if (mu3.gt.nu3) stop 'mu3 > nu3'
        do i3=ml3,mu3
        dz2=(i3*hgrid-rxyz(3,iat))**2
        do i2=ml2,mu2
        dy2=(i2*hgrid-rxyz(2,iat))**2
        do i1=ml1,mu1
        dx=i1*hgrid-rxyz(1,iat)
        if (dx**2+(dy2+dz2).lt.rad**2) then 
              logrid(i1,i2,i3)=.true.
        endif
        enddo ; enddo ; enddo
      enddo

        return
     end subroutine fill_logrid
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

	subroutine loewe_p(iproc,nproc,norb,norbp,nvctrp,psit)
! loewdin orthogonalisation
        implicit real*8 (a-h,o-z)
        dimension psit(nvctrp,norbp*nproc)
        real*8, allocatable :: ovrlp(:,:,:),evall(:),psitt(:,:)
        include 'mpif.h'

 if (norb.eq.1) stop 'more than one orbital needed for a parallel run'

        allocate(ovrlp(norb,norb,3),evall(norb))

       call cpu_time(tr0)
       call system_clock(ncount0,ncount_rate,ncount_max)

! Upper triangle of overlap matrix using BLAS
!     ovrlp(iorb,jorb)=psit(k,iorb)*psit(k,jorb) ; upper triangle
        call DSYRK('U','T',norb,nvctrp,1.d0,psit,nvctrp,0.d0,ovrlp(1,1,2),norb)

! Full overlap matrix using  BLAS
!     ovrlap(jorb,iorb,2)=+psit(k,jorb)*psit(k,iorb)
!      call DGEMM('T','N',norb,norb,nvctrp,1.d0,psit,nvctrp,psit,nvctrp,0.d0,ovrlp(1,1,2),norb)

       call cpu_time(tr1)
       call system_clock(ncount1,ncount_rate,ncount_max)
        call MPI_ALLREDUCE (ovrlp(1,1,2),ovrlp(1,1,1),norb**2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
       call cpu_time(tr2)
       call system_clock(ncount2,ncount_rate,ncount_max)

! LAPACK
        call DSYEV('V','U',norb,ovrlp(1,1,1),norb,evall,ovrlp(1,1,3),norb**2,info)
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

        allocate(psitt(nvctrp,norbp*nproc))
! new eigenvectors
!   psitt(i,iorb)=psit(i,jorb)*ovrlp(jorb,iorb,3)
      call DGEMM('N','N',nvctrp,norb,norb,1.d0,psit,nvctrp,ovrlp(1,1,3),norb,0.d0,psitt,nvctrp)
      call DCOPY(nvctrp*norbp*nproc,psitt,1,psit,1)
      deallocate(psitt)


       call cpu_time(tr3)
       call system_clock(ncount3,ncount_rate,ncount_max)

        deallocate(ovrlp,evall)

       tel=dble(ncount1-ncount0+ncount3-ncount2)/dble(ncount_rate)
       write(77,'(a40,i4,2(x,e10.3))') 'LOEWE TIME',iproc,tr1-tr0+tr3-tr2,tel
       tel=dble(ncount2-ncount1)/dble(ncount_rate)
       write(78,'(a40,i4,2(x,e10.3))') 'LOEWE TIME',iproc,tr2-tr1,tel

        return
     end subroutine loewe_p



	subroutine loewe(norb,norbp,nvctrp,psi)
! loewdin orthogonalisation
        implicit real*8 (a-h,o-z)
        dimension psi(nvctrp,norbp)
        real*8, allocatable :: ovrlp(:,:,:),evall(:),tpsi(:,:)

       call cpu_time(tr0)
       call system_clock(ncount0,ncount_rate,ncount_max)

 if (norb.eq.1) then
        tt=0.d0
        do i=1,nvctrp
        tt=tt+psi(i,1)**2
        enddo
        tt=1.d0/sqrt(tt)
        do i=1,nvctrp
        psi(i,1)=psi(i,1)*tt
        enddo

 else

        allocate(ovrlp(norb,norb,3),evall(norb))

! Overlap matrix using BLAS
!     ovrlp(iorb,jorb)=psi(k,iorb)*psi(k,jorb) ; upper triangle
        call DSYRK('U','T',norb,nvctrp,1.d0,psi,nvctrp,0.d0,ovrlp(1,1,1),norb)

! LAPACK
        call DSYEV('V','U',norb,ovrlp(1,1,1),norb,evall,ovrlp(1,1,3),norb**2,info)
        if (info.ne.0) write(6,*) 'info loewe', info
!          write(6,*) 'overlap eigenvalues'
!77        format(8(1x,e10.3))
!          if (norb.le.16) then
!          write(6,77) evall
!          else
!          write(6,77) (evall(i),i=1,4), (evall(i),i=norb-3,norb)
!          endif

! calculate S^{-1/2} ovrlp(*,*,3)
        do 3935,lorb=1,norb
        do 3935,jorb=1,norb
3935    ovrlp(jorb,lorb,2)=ovrlp(jorb,lorb,1)*sqrt(1.d0/evall(lorb))
!        do 3985,j=1,norb
!        do 3985,i=1,norb
!        ovrlp(i,j,3)=0.d0
!        do 3985,l=1,norb
!3985    ovrlp(i,j,3)=ovrlp(i,j,3)+ovrlp(i,l,1)*ovrlp(j,l,2)
! BLAS:
        call DGEMM('N','T',norb,norb,norb,1.d0,ovrlp(1,1,1),norb,ovrlp(1,1,2),norb,0.d0,ovrlp(1,1,3),norb)

! new eigenvectors
      allocate(tpsi(nvctrp,norb))
!   tpsi(i,iorb)=psi(i,jorb)*ovrlp(jorb,iorb,3)
      call DGEMM('N','N',nvctrp,norb,norb,1.d0,psi,nvctrp,ovrlp(1,1,3),norb,0.d0,tpsi,nvctrp)
      call DCOPY(nvctrp*norb,tpsi,1,psi,1)
      deallocate(tpsi)

        deallocate(ovrlp,evall)

 endif

       call cpu_time(tr1)
       call system_clock(ncount1,ncount_rate,ncount_max)
       tel=dble(ncount1-ncount0)/dble(ncount_rate)
       write(77,'(a40,i4,2(x,e10.3))') 'LOEWE TIME',iproc,tr1-tr0,tel

        return
     end subroutine loewe


	subroutine checkortho_p(iproc,nproc,norb,norbp,nvctrp,psit)
        implicit real*8 (a-h,o-z)
        dimension psit(nvctrp,norbp*nproc)
        real*8, allocatable :: ovrlp(:,:,:)
        include 'mpif.h'

        allocate(ovrlp(norb,norb,2))

     do 100,iorb=1,norb
        do 100,jorb=1,norb
        ovrlp(iorb,jorb,2)=ddot(nvctrp,psit(1,iorb),1,psit(1,jorb),1)
100    continue

     call MPI_ALLREDUCE(ovrlp(1,1,2),ovrlp(1,1,1),norb**2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

        toler=1.d-10
        dev=0.d0
     do 110,iorb=1,norb
        do 110,jorb=1,norb
        scpr=ovrlp(iorb,jorb,1)
        if (iorb.eq.jorb) then
        dev=dev+(scpr-1.d0)**2
        else
        dev=dev+scpr**2
        endif
        if (iorb.eq.jorb .and. abs(scpr-1.d0).gt.toler)  write(*,*) 'ERROR ORTHO',iorb,jorb,scpr
        if (iorb.ne.jorb .and. abs(scpr).gt.toler)  write(*,*) 'ERROR ORTHO',iorb,jorb,scpr
110     continue

        write(*,*) 'Deviation from orthogonality ',iproc,dev

        deallocate(ovrlp)

        return
     end subroutine checkortho_p


	subroutine checkortho(norb,norbp,nvctrp,psi)
        implicit real*8 (a-h,o-z)
        dimension psi(nvctrp,norbp)
        real*8, allocatable :: ovrlp(:,:,:)

        allocate(ovrlp(norb,norb,1))

     do iorb=1,norb
     do jorb=1,norb
        ovrlp(iorb,jorb,1)=ddot(nvctrp,psi(1,iorb),1,psi(1,jorb),1)
     enddo ; enddo

        toler=1.d-10
        dev=0.d0
     do iorb=1,norb
     do jorb=1,norb
        scpr=ovrlp(iorb,jorb,1)
        if (iorb.eq.jorb) then
        dev=dev+(scpr-1.d0)**2
        else
        dev=dev+scpr**2
        endif
        if (iorb.eq.jorb .and. abs(scpr-1.d0).gt.toler)  write(*,*) 'ERROR ORTHO',iorb,jorb,scpr
        if (iorb.ne.jorb .and. abs(scpr).gt.toler)  write(*,*) 'ERROR ORTHO',iorb,jorb,scpr
     enddo  ; enddo

        write(*,*) 'Deviation from orthogonality ',0,dev

        deallocate(ovrlp)


        return
	end subroutine checkortho

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

        subroutine myorbital(myorb, iorb,norbe,iproc,nproc)
        implicit real*8 (a-h,o-z)
        logical myorb
        parameter(eps_mach=1.d-12)

        tt=dble(norbe)/dble(nproc)
        norbep=int((1.d0-eps_mach*tt) + tt)
        if (iorb .ge. iproc*norbep+1 .and. iorb .le. min((iproc+1)*norbep,norbe)) then
        myorb=.true.
        else
        myorb=.false.
        endif

        return
        end subroutine myorbital


	subroutine KStrans_p(iproc,nproc,norb,norbp,nvctrp,occup,  & 
                           hpsit,psit,evsum,cprec)
! at the start each processor has all the Psi's but only its part of the HPsi's
! at the end each processor has only its part of the Psi's
        implicit real*8 (a-h,o-z)
        dimension occup(norb)
        dimension psit(nvctrp,norbp*nproc),hpsit(nvctrp,norbp*nproc)
! arrays for KS orbitals
        allocatable :: hamks(:,:,:),eval(:),work_lp(:),psitt(:,:)
        include 'mpif.h'

! set up Hamiltonian matrix
        allocate(hamks(norb,norb,2))
        do jorb=1,norb
        do iorb=1,norb
        hamks(iorb,jorb,2)=0.d0
	enddo
	enddo
        do iorb=1,norb
        do jorb=1,norb
        scpr=ddot(nvctrp,psit(1,jorb),1,hpsit(1,iorb),1)
        hamks(iorb,jorb,2)=scpr
        enddo
        enddo


        call MPI_ALLREDUCE (hamks(1,1,2),hamks(1,1,1),norb**2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        
!        write(*,*) 'KS Hamiltonian',iproc
!        do iorb=1,norb
!        write(*,'(10(1x,e10.3))') (hamks(iorb,jorb,1),jorb=1,norb)
!        enddo

        n_lp=1000
        allocate(work_lp(n_lp),eval(norb))
        call  DSYEV('V','U',norb,hamks,norb,eval, work_lp, n_lp, info )
        cprec=eval(1)
        evsum=0.d0
        do iorb=1,norb
        evsum=evsum+eval(iorb)*occup(iorb)
        if (iproc.eq.0) write(*,*) 'eval(',iorb,')=',eval(iorb)
        enddo
        deallocate(work_lp,eval)
        if (info.ne.0) write(*,*) 'DSYEV ERROR',info

        allocate(psitt(nvctrp,norbp*nproc))
! Transform to KS orbitals
      do iorb=1,norb
        call zero(nvctrp,psitt(1,iorb))
        do jorb=1,norb
        alpha=hamks(jorb,iorb,1)
        call daxpy(nvctrp,alpha,psit(1,jorb),1,psitt(1,iorb),1)
        enddo
      enddo
        deallocate(hamks)

        call DCOPY(nvctrp*norbp*nproc,psitt,1,psit,1)
        deallocate(psitt)

	return
        end subroutine KStrans_p



	subroutine KStrans(norb,norbp,nvctrp,occup,hpsi,psi,evsum,cprec)
! at the start each processor has all the Psi's but only its part of the HPsi's
! at the end each processor has only its part of the Psi's
        implicit real*8 (a-h,o-z)
        dimension occup(norb)
        dimension psi(nvctrp,norbp),hpsi(nvctrp,norbp)
! arrays for KS orbitals
        allocatable :: hamks(:,:,:),eval(:),work_lp(:),psitt(:,:)

! set up Hamiltonian matrix
        allocate(hamks(norb,norb,2))
        do jorb=1,norb
        do iorb=1,norb
        hamks(iorb,jorb,2)=0.d0
	enddo
	enddo
        do iorb=1,norb
        do jorb=1,norb
        scpr=ddot(nvctrp,psi(1,jorb),1,hpsi(1,iorb),1)
        hamks(iorb,jorb,1)=scpr
        enddo
        enddo


!        write(*,*) 'KS Hamiltonian',0
!        do iorb=1,norb
!        write(*,'(10(1x,e10.3))') (hamks(iorb,jorb,1),jorb=1,norb)
!        enddo

        n_lp=1000
        allocate(work_lp(n_lp),eval(norb))
        call  DSYEV('V','U',norb,hamks,norb,eval, work_lp, n_lp, info )
        cprec=eval(1)
        evsum=0.d0
        do iorb=1,norb
        evsum=evsum+eval(iorb)*occup(iorb)
        write(*,*) 'eval(',iorb,')=',eval(iorb)
        enddo
        deallocate(work_lp,eval)
        if (info.ne.0) write(*,*) 'DSYEV ERROR',info

        allocate(psitt(nvctrp,norbp))
! Transform to KS orbitals
      do iorb=1,norb
        call zero(nvctrp,psitt(1,iorb))
        do jorb=1,norb
        alpha=hamks(jorb,iorb,1)
        call daxpy(nvctrp,alpha,psi(1,jorb),1,psitt(1,iorb),1)
        enddo
      enddo
        deallocate(hamks)

        call DCOPY(nvctrp*norbp,psitt,1,psi,1)
        deallocate(psitt)

	return
      end subroutine KStrans



        subroutine crtonewave(n1,n2,n3,nterm,lx,ly,lz,xp,psiat,rx,ry,rz,hgrid, & 
                   nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,  & 
                   nseg_c,mvctr_c,keyg_c,keyv_c,nseg_f,mvctr_f,keyg_f,keyv_f,psi_c,psi_f)
! returns an input guess orbital that is a Gaussian centered at a Wannier center
! exp (-1/(2*gau_a^2) *((x-cntrx)^2 + (y-cntry)^2 + (z-cntrz)^2 ))
! in the arrays psi_c, psi_f
        implicit real*8 (a-h,o-z)
        parameter(nw=16000)
        dimension xp(nterm),psiat(nterm)
        dimension keyg_c(2,nseg_c),keyv_c(nseg_c),keyg_f(2,nseg_f),keyv_f(nseg_f)
        dimension psi_c(mvctr_c),psi_f(7,mvctr_f)
        real*8, allocatable, dimension(:,:) :: wprojx, wprojy, wprojz
        real*8, allocatable, dimension(:,:) :: work
        real*8, allocatable :: psig_c(:,:,:), psig_f(:,:,:,:)

        allocate(wprojx(0:n1,2),wprojy(0:n2,2),wprojz(0:n3,2),work(0:nw,2))
        allocate(psig_c(nl1_c:nu1_c,nl2_c:nu2_c,nl3_c:nu3_c))
        allocate(psig_f(7,nl1_f:nu1_f,nl2_f:nu2_f,nl3_f:nu3_f))


      iterm=1
        gau_a=xp(iterm)
        n_gau=lx
        CALL GAUSS_TO_DAUB(hgrid,1.d0,rx,gau_a,n_gau,n1,ml1,mu1,wprojx(0,1),te,work,nw)
        n_gau=ly
        CALL GAUSS_TO_DAUB(hgrid,1.d0,ry,gau_a,n_gau,n2,ml2,mu2,wprojy(0,1),te,work,nw)
        n_gau=lz
        CALL GAUSS_TO_DAUB(hgrid,psiat(iterm),rz,gau_a,n_gau,n3,ml3,mu3,wprojz(0,1),te,work,nw)

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

      do 100,iterm=2,nterm
        gau_a=xp(iterm)
        n_gau=lx
        CALL GAUSS_TO_DAUB(hgrid,1.d0,rx,gau_a,n_gau,n1,ml1,mu1,wprojx(0,1),te,work,nw)
        n_gau=ly
        CALL GAUSS_TO_DAUB(hgrid,1.d0,ry,gau_a,n_gau,n2,ml2,mu2,wprojy(0,1),te,work,nw)
        n_gau=lz
        CALL GAUSS_TO_DAUB(hgrid,psiat(iterm),rz,gau_a,n_gau,n3,ml3,mu3,wprojz(0,1),te,work,nw)

! First term: coarse projector components
          do i3=nl3_c,nu3_c
          do i2=nl2_c,nu2_c
          do i1=nl1_c,nu1_c
            psig_c(i1,i2,i3)=psig_c(i1,i2,i3)+wprojx(i1,1)*wprojy(i2,1)*wprojz(i3,1)
          enddo ; enddo ; enddo

! First term: fine projector components
          do i3=nl3_f,nu3_f
          do i2=nl2_f,nu2_f
          do i1=nl1_f,nu1_f
            psig_f(1,i1,i2,i3)=psig_f(1,i1,i2,i3)+wprojx(i1,2)*wprojy(i2,1)*wprojz(i3,1)
            psig_f(2,i1,i2,i3)=psig_f(2,i1,i2,i3)+wprojx(i1,1)*wprojy(i2,2)*wprojz(i3,1)
            psig_f(3,i1,i2,i3)=psig_f(3,i1,i2,i3)+wprojx(i1,2)*wprojy(i2,2)*wprojz(i3,1)
            psig_f(4,i1,i2,i3)=psig_f(4,i1,i2,i3)+wprojx(i1,1)*wprojy(i2,1)*wprojz(i3,2)
            psig_f(5,i1,i2,i3)=psig_f(5,i1,i2,i3)+wprojx(i1,2)*wprojy(i2,1)*wprojz(i3,2)
            psig_f(6,i1,i2,i3)=psig_f(6,i1,i2,i3)+wprojx(i1,1)*wprojy(i2,2)*wprojz(i3,2)
            psig_f(7,i1,i2,i3)=psig_f(7,i1,i2,i3)+wprojx(i1,2)*wprojy(i2,2)*wprojz(i3,2)
          enddo ; enddo ; enddo

100	continue

! coarse part
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
            psi_c(i-i0+jj)=psig_c(i,i2,i3)
          enddo
        enddo

! fine part
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
     end subroutine crtonewave
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


        subroutine readmywaves(iproc,norb,norbp,n1,n2,n3,hgrid,rxyz,  & 
                               nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,cprec)
! reads wavefunction from file and transforms it properly if hgrid or size of simulation cell have changed
        implicit real*8 (a-h,o-z)
        character*50 filename
        character*4 f4
        logical cif1,cif2,cif3
        dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
        dimension psi(nvctr_c+7*nvctr_f,norbp)
        dimension xya(-1:1,-1:1),xa(-1:1)
        dimension rxyz(3),rxyz_old(3)
        allocatable :: psifscf(:,:,:)
        allocatable :: psigold(:,:,:,:,:,:),psifscfold(:,:,:),psifscfoex(:,:,:),wwold(:)

       call cpu_time(tr0)
       call system_clock(ncount1,ncount_rate,ncount_max)

        allocate(psifscf(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8))

  do 100,iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)

        write(f4,'(i4.4)') iorb
!        filename = '/scratch/tmp/stefan/wavefunction.'//f4
        filename = 'wavefunction.'//f4
        open(unit=99,file=filename,status='unknown')

         read(99,*) iorbold,cprec
         if (iorbold.ne.iorb) stop 'readallwaves'
         read(99,*) hgridold
         read(99,*) nold1,nold2,nold3
         read(99,*) (rxyz_old(j),j=1,3)
         write(*,'(i2,6(x,e14.7))') iproc,(rxyz(j),j=1,3),(rxyz_old(j),j=1,3)
         read(99,*) nvctrold_c, nvctrold_f

!        write(*,*) iorb,' hgridold,hgrid ',hgridold,hgrid
!        write(*,*) iorb,' nvctrold_c,nvctr_c ',nvctrold_c,nvctr_c
!        write(*,*) iorb,' nvctrold_f,nvctr_f ',nvctrold_f,nvctr_f
!        write(*,*) iorb,' nold1,n1 ',nold1,n1
!        write(*,*) iorb,' nold2,n2 ',nold2,n2
!        write(*,*) iorb,' nold3,n3 ',nold3,n3

      if (hgridold.eq. hgrid .and. nvctrold_c.eq.nvctr_c .and. nvctrold_f.eq.nvctr_f  & 
          .and. nold1.eq.n1  .and. nold2.eq.n2 .and. nold3.eq.n3 ) then

       write(*,*) 'wavefunction ',iorb,' needs NO transformation on processor',iproc
	do j=1,nvctrold_c
            read(99,*) i1,i2,i3,tt
            psi(j,iorb-iproc*norbp)=tt
         enddo
	do j=1,7*nvctrold_f-6,7
            read(99,*) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
            psi(nvctr_c+j+0,iorb-iproc*norbp)=t1
            psi(nvctr_c+j+1,iorb-iproc*norbp)=t2
            psi(nvctr_c+j+2,iorb-iproc*norbp)=t3
            psi(nvctr_c+j+3,iorb-iproc*norbp)=t4
            psi(nvctr_c+j+4,iorb-iproc*norbp)=t5
            psi(nvctr_c+j+5,iorb-iproc*norbp)=t6
            psi(nvctr_c+j+6,iorb-iproc*norbp)=t7
         enddo

       else
       write(*,*) 'wavefunction ',iorb,' needs transformation on processor',iproc
       if (hgridold.ne.hgrid) write(*,*) 'because hgridold >< hgrid'
       if (nvctrold_c.eq.nvctr_c) write(*,*) 'because nvctrold_c >< nvctr_c'
       if (nvctrold_f.eq.nvctr_f) write(*,*) 'because nvctrold_f >< nvctr_f'
       if (nold1.ne.n1  .or. nold2.ne.n2 .or. nold3.ne.n3 ) write(*,*) 'because cell size has changed'

         allocate(psigold(0:nold1,2,0:nold2,2,0:nold3,2),  & 
                  psifscfold(-7:2*nold1+8,-7:2*nold2+8,-7:2*nold3+8), &
                  psifscfoex(-8:2*nold1+9,-8:2*nold2+9,-8:2*nold3+9), &
                   wwold((2*nold1+16)*(2*nold2+16)*(2*nold3+16)))

         call zero(8*(nold1+1)*(nold2+1)*(nold3+1),psigold)
	do iel=1,nvctrold_c
            read(99,*) i1,i2,i3,tt
            psigold(i1,1,i2,1,i3,1)=tt
         enddo
	do iel=1,nvctrold_f
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
	call synthese_grow(nold1,nold2,nold3,wwold,psigold,psifscfold)
!        call plot_psifscf(50+iorb,hgridold,nold1,nold2,nold3,psifscfold,0.d0,0.d0,0.d0)

         do i3=-7,2*nold3+8
         do i2=-7,2*nold2+8
           i1=-8
           psifscfoex(i1,i2,i3)=0.d0
           do i1=-7,2*nold1+8
             psifscfoex(i1,i2,i3)=psifscfold(i1,i2,i3)
           enddo
           i1=2*nold1+9
           psifscfoex(i1,i2,i3)=0.d0
         enddo
         enddo

         i3=-8
         do i2=-8,2*nold2+9
         do i1=-8,2*nold1+9
           psifscfoex(i1,i2,i3)=0.d0
         enddo
         enddo
         i3=2*nold3+9
         do i2=-8,2*nold2+9
         do i1=-8,2*nold1+9
           psifscfoex(i1,i2,i3)=0.d0
         enddo
         enddo

         i2=-8
         do i3=-8,2*nold3+9
         do i1=-8,2*nold1+9
           psifscfoex(i1,i2,i3)=0.d0
         enddo
         enddo
         i2=2*nold2+9
         do i3=-8,2*nold3+9
         do i1=-8,2*nold1+9
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
        call zero((2*n1+16)*(2*n2+16)*(2*n3+16),psifscf)
        do i3=-7,2*n3+8
        z=i3*hgridh
        j3=nint((z-dz)/hgridhold)
        cif3=(j3.ge.-7 .and. j3.le.2*nold3+8)
        do i2=-7,2*n2+8
        y=i2*hgridh
        j2=nint((y-dy)/hgridhold)
        cif2=(j2.ge.-7 .and. j2.le.2*nold2+8)
        do i1=-7,2*n1+8
        x=i1*hgridh
        j1=nint((x-dx)/hgridhold)
        cif1=(j1.ge.-7 .and. j1.le.2*nold1+8)

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

        call   compress(n1,n2,n3,nseg_c,nvctr_c,keyg(1,1),keyv(1),   & 
                    nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   & 
                    psifscf,psi(1,iorb-iproc*norbp),psi(nvctr_c+1,iorb-iproc*norbp))

         deallocate(psigold,psifscfold,psifscfoex,wwold)

      endif

        close(99)
  100  continue

         deallocate(psifscf)

       call cpu_time(tr1)
       call system_clock(ncount2,ncount_rate,ncount_max)
       tel=dble(ncount2-ncount1)/dble(ncount_rate)
       write(77,'(a40,i4,2(x,e10.3))') 'READING WAVES TIME',iproc,tr1-tr0,tel



	end subroutine readmywaves



        subroutine writemywaves(iproc,norb,norbp,n1,n2,n3,hgrid,  & 
                   rxyz,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,cprec)
! write all my wavefunctions in files by calling writeonewave
        implicit real*8 (a-h,o-z)
        dimension rxyz(3)
        dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
        dimension psi(nvctr_c+7*nvctr_f,norbp)

       call cpu_time(tr0)
       call system_clock(ncount1,ncount_rate,ncount_max)

       do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)

       call writeonewave(iorb,n1,n2,n3,hgrid,rxyz,  & 
                         nseg_c,nvctr_c,keyg(1,1),keyv(1)  & 
                        ,nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), & 
                        psi(1,iorb-iproc*norbp),psi(nvctr_c+1,iorb-iproc*norbp),cprec)
       enddo

       call cpu_time(tr1)
       call system_clock(ncount2,ncount_rate,ncount_max)
       tel=dble(ncount2-ncount1)/dble(ncount_rate)
       write(77,'(a40,i4,2(x,e10.3))') 'WRITE WAVES TIME',iproc,tr1-tr0,tel


       return
     end subroutine writemywaves



        subroutine writeonewave(iorb,n1,n2,n3,hgrid,rxyz,  & 
                           nseg_c,nvctr_c,keyg_c,keyv_c,  & 
                           nseg_f,nvctr_f,keyg_f,keyv_f, & 
                              psi_c,psi_f,cprec)
        implicit real*8 (a-h,o-z)
        character*50 filename
        character*4 f4
        dimension keyg_c(2,nseg_c),keyv_c(nseg_c),keyg_f(2,nseg_f),keyv_f(nseg_f)
        dimension psi_c(nvctr_c),psi_f(7,nvctr_f),rxyz(3)


        write(f4,'(i4.4)') iorb
!        filename = '/scratch/tmp/stefan/wavefunction.'//f4
        filename = 'wavefunction.'//f4
        write(*,*) 'opening ',filename
        open(unit=99,file=filename,status='unknown')
         write(99,*) iorb,cprec
         write(99,*) hgrid
         write(99,*) n1,n2,n3
         write(99,'(3(1x,e24.17))') (rxyz(j),j=1,3)
         write(99,*) nvctr_c, nvctr_f

! coarse part
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
            tt=psi_c(i-i0+jj) 
            write(99,'(3(i4),1x,e19.12)') i,i2,i3,tt
          enddo
         enddo
                                                                                                                             
! fine part
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


      end subroutine writeonewave
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

        subroutine plot_wf(n1,n2,n3,hgrid,nseg_c,nvctr_c,keyg_c,keyv_c,nseg_f,nvctr_f,keyg_f,keyv_f, & 
                                        psi_c,psi_f)
        implicit real*8 (a-h,o-z)
        dimension keyg_c(2,nseg_c),keyv_c(nseg_c),keyg_f(2,nseg_f),keyv_f(nseg_f)
        dimension psi_c(nvctr_c),psi_f(7,nvctr_f)
        real*8, allocatable, dimension(:) :: psifscf,psir

        allocate(psifscf((2*n1+16)*(2*n2+16)*(2*n3+16)))
        allocate(psir((2*n1+31)*(2*n2+31)*(2*n3+31)))

        call uncompress(n1,n2,n3,nseg_c,nvctr_c,keyg_c,keyv_c,nseg_f,nvctr_f,keyg_f,keyv_f, & 
                        psi_c,psi_f,psifscf)
        call convolut_magic_n(2*n1+15,2*n2+15,2*n3+15,psifscf,psir) 

        call plot_pot(hgrid,n1,n2,n3,10,psir)

        deallocate(psifscf,psir)
        return
      end subroutine plot_wf




        subroutine plot_pot(hgrid,n1,n2,n3,iounit,pot)
        implicit real*8 (a-h,o-z)
        dimension pot(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16)

        hgridh=.5d0*hgrid
         rx=8.d0
         ry=8.d0
         rz=8.d0
        open(iounit) 
        open(iounit+1) 
        open(iounit+2) 

        i3=nint(rz/hgridh)
        i2=nint(ry/hgridh)
        write(*,*) 'plot_p, i2,i3 ',i2,i3
        do i1=-14,2*n1+16
        write(iounit,*) i1*hgridh,pot(i1,i2,i3)
        enddo

        i1=nint(rx/hgridh)
        i2=nint(ry/hgridh)
        write(*,*) 'plot_p, i1,i2 ',i1,i2
        do i3=-14,2*n3+16
        write(iounit+1,*) i3*hgridh,pot(i1,i2,i3)
        enddo

        i1=nint(rx/hgridh)
        i3=nint(rz/hgridh)
        write(*,*) 'plot_p, i1,i3 ',i1,i3
        do i2=-14,2*n2+16
        write(iounit+2,*) i2*hgridh,pot(i1,i2,i3)
        enddo

        close(iounit) 
        close(iounit+1) 
        close(iounit+2) 

        return
      end subroutine plot_pot



        subroutine plot_psifscf(iunit,hgrid,n1,n2,n3,psifscf)
        implicit real*8 (a-h,o-z)
        dimension psifscf(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8)

	hgridh=.5d0*hgrid

! along x-axis
	i3=n3
	i2=n2
	do i1=-7,2*n1+8
            write(iunit,'(3(1x,e10.3),1x,e12.5)') i1*hgridh,i2*hgridh,i3*hgridh,psifscf(i1,i2,i3)
	enddo 

! 111 diagonal
	do i=-7,min(2*n1+8,2*n2+8,2*n3+8)
        i1=i ; i2=i ; i3=i
            write(iunit,'(3(1x,e10.3),1x,e12.5)') i1*hgridh,i2*hgridh,i3*hgridh,psifscf(i1,i2,i3)
	enddo 

! 1-1-1 diagonal
	do i=-7,min(2*n1+8,2*n2+8,2*n3+8)
        i1=i ; i2=-i ; i3=-i
            write(iunit,'(3(1x,e10.3),1x,e12.5)') i1*hgridh,i2*hgridh,i3*hgridh,psifscf(i1,i2,i3)
	enddo 

! -11-1 diagonal
	do i=-7,min(2*n1+8,2*n2+8,2*n3+8)
        i1=-i ; i2=i ; i3=-i
            write(iunit,'(3(1x,e10.3),1x,e12.5)') i1*hgridh,i2*hgridh,i3*hgridh,psifscf(i1,i2,i3)
	enddo 

! -1-11 diagonal
	do i=-7,min(2*n1+8,2*n2+8,2*n3+8)
        i1=-i ; i2=-i ; i3=i
            write(iunit,'(3(1x,e10.3),1x,e12.5)') i1*hgridh,i2*hgridh,i3*hgridh,psifscf(i1,i2,i3)
	enddo 

        return
      end subroutine plot_psifscf

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

      subroutine atomkin(l,ng,xp,psiat,psiatn,ek)
! calculates the kinetic energy of an atomic wavefunction expressed in Gaussians
        implicit real*8 (a-h,o-z)
        dimension xp(ng),psiat(ng),psiatn(ng)

!        gml=.5d0*gamma(.5d0+l)
        if (l.eq.0) then 
            gml=0.88622692545275801365d0
        else if (l.eq.1) then 
            gml=0.44311346272637900682d0
!        else if (l.eq.2) then 
!            gml=0.66467019408956851024d0
!        else if (l.eq.3) then 
!            gml=1.6616754852239212756d0
        else
          stop 'atomkin'
        endif

        ek=0.d0
        tt=0.d0
        do i=1,ng
        xpi=.5d0/xp(i)**2
        do j=1,ng
        xpj=.5d0/xp(j)**2
        d=xpi+xpj
        sxp=1.d0/d
        const=gml*sqrt(sxp)**(2*l+1)
! kinetic energy  matrix element hij
        hij=.5d0*const*sxp**2* ( 3.d0*xpi*xpj +                  &
                     l*(6.d0*xpi*xpj-xpi**2-xpj**2) -        &
                     l**2*(xpi-xpj)**2  ) + .5d0*l*(l+1.d0)*const
        sij=const*sxp*(l+.5d0)
        ek=ek+hij*psiat(i)*psiat(j)
        tt=tt+sij*psiat(i)*psiat(j)
        enddo
        enddo

        if (abs(tt-1.d0).gt.1.d-2) write(*,*) 'presumably wrong inguess data',tt
! energy expectation value
        ek=ek/tt
!        write(*,*) 'ek=',ek,tt,l,ng
! scale atomic wavefunction
        tt=sqrt(1.d0/tt)
        if (l.eq.0) then  ! multiply with 1/sqrt(4*pi)
        tt=tt*0.28209479177387814347d0
        else if (l.eq.1) then  ! multiply with sqrt(3/(4*pi))
        tt=tt*0.48860251190291992159d0
        endif
        do i=1,ng
        psiatn(i)=psiat(i)*tt
        enddo

        return
      end subroutine atomkin

