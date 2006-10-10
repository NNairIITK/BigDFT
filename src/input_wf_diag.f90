	subroutine input_wf_diag(parallel,iproc,nproc,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
                    nat,norb,norbp,n1,n2,n3,nfft1,nfft2,nfft3,nvctr_c,nvctr_f,nvctrp,hgrid,rxyz, & 
                   rhopot,pot_ion,nseg_c,nseg_f,keyg,keyv, &
                    nprojel,nproj,nseg_p,keyg_p,keyv_p,nvctr_p,proj,  &
                    atomnames,ntypes,iatype,pkernel,psppar,ppsi,cprec,accurex)
! Input wavefunctions are found by a diagonalization in a minimal basis set
! Each processors writes its initial wavefunctions into the wavefunction file
! The files are then read by readwave
        implicit real*8 (a-h,o-z)
        logical parallel, myorb
	parameter (ngx=31)
	parameter (npsp=15)
	character*20 atomnames(100)
	character*20 pspatomnames(npsp)
        parameter(eps_mach=1.d-12)
	dimension  xp(ngx,npsp),psiat(ngx,2,npsp),occupat(2,npsp),ng(npsp),ns(npsp),np(npsp),psiatn(ngx)
        dimension rxyz(3,nat),iatype(nat)
        dimension rhopot((2*n1+31)*(2*n2+31)*(2*n3+31)),pot_ion((2*n1+31)*(2*n2+31)*(2*n3+31))
	dimension pkernel(*)
        dimension psppar(0:2,0:4,ntypes)
        dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
        dimension nseg_p(0:2*nat),nvctr_p(0:2*nat)
        dimension keyg_p(2,nseg_p(2*nat)),keyv_p(nseg_p(2*nat))
        dimension proj(nprojel)
        dimension ppsi(nvctr_c+7*nvctr_f,norbp)

        allocatable :: psi(:,:),hpsi(:,:),psit(:,:),hpsit(:,:),ppsit(:,:),occupe(:)
        allocatable :: hamovr(:,:,:),evale(:),work_lp(:)
        include 'mpif.h'


  call readAtomicOrbitals(iproc, 'inguess.dat', ngx, npsp, xp, psiat, occupat, ng, &
                        & ns, np, pspatomnames, norbe, atomnames, iatype, nat)
  
!  allocate wavefunctions and their occupation numbers
        allocate(occupe(norbe))
        tt=dble(norbe)/dble(nproc)
        norbep=int((1.d0-eps_mach*tt) + tt)
        write(79,'(a40,i10)') 'words for (h)psi inguess',2*(nvctr_c+7*nvctr_f)*norbep
        allocate(psi(nvctr_c+7*nvctr_f,norbep))
        allocate(hpsi(nvctr_c+7*nvctr_f,norbep))
        write(79,*) 'allocation done'
        norbeme=max(min((iproc+1)*norbep,norbe)-iproc*norbep,0)
        write(*,*) 'iproc ',iproc,' treats ',norbeme,' orbitals '

        hgridh=.5d0*hgrid

    call createAtomicOrbitals(iproc, npsp, pspatomnames, atomnames, nat, rxyz, &
  & norbe, norbep, occupe, ngx, xp, psiat, psiatn, ng, ns, np, nvctr_c, nvctr_f, &
  & n1, n2, n3, hgrid, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3, nseg_c, nseg_f, &
  & keyg, keyv, iatype, psi, eks)
  

! resulting charge density and potential
       call sumrho(parallel,iproc,norbe,norbeme,n1,n2,n3,hgrid,occupe,  & 
                   nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,rhopot)
      if (parallel) then
          call ParPSolver_Kernel(2*n1+31,2*n2+31,2*n3+31,nfft1,nfft2,nfft3,hgridh,pkernel,.true., &
               pot_ion,rhopot,ehart,eexcu,vexcu,iproc,nproc)
       else
          call PSolver_Kernel(2*n1+31,2*n2+31,2*n3+31,nfft1,nfft2,nfft3,hgridh,pkernel,.true., &
               pot_ion,rhopot,ehart,eexcu,vexcu)
       end if

! set up subspace Hamiltonian 
!        allocate(hamovr(norbe,norbe,4))

        call applylocpotkin(iproc,norbe,norbeme,n1,n2,n3,hgrid,occupe,  & 
                   nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,psi,rhopot,hpsi,epot_sum,ekin_sum)

       if (parallel) then
       tt=ekin_sum
       call MPI_ALLREDUCE(tt,ekin_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
       endif

       accurex=abs(eks-ekin_sum)
       write(*,*) 'ekin_sum,eks',ekin_sum,eks

        call applyprojectors(iproc,ntypes,nat,iatype,psppar,occupe, &
                    nprojel,nproj,nseg_p,keyg_p,keyv_p,nvctr_p,proj,  &
                    norbe,norbeme,nseg_c,nseg_f,keyg,keyv,nvctr_c,nvctr_f,psi,hpsi,eproj_sum)

      call solveKS(parallel, iproc, nproc, norb, norbp, norbe, norbep, nvctr_c, nvctr_f, &
        & nvctrp, psi, hpsi, ppsi, cprec)
        
      deallocate(occupe, hpsi, psi)

! if (parallel) then
!        write(79,'(a40,i10)') 'words for psit inguess',nvctrp*norbep*nproc
!        allocate(psit(nvctrp,norbep*nproc))
!        write(79,*) 'allocation done'
!
!        call  transallwaves(iproc,nproc,norbe,norbep,nvctr_c,nvctr_f,nvctrp,psi,psit)
!
!        deallocate(psi)
!
!        write(79,'(a40,i10)') 'words for hpsit inguess',2*nvctrp*norbep*nproc
!        allocate(hpsit(nvctrp,norbep*nproc))
!        write(79,*) 'allocation done'
!
!        call  transallwaves(iproc,nproc,norbe,norbep,nvctr_c,nvctr_f,nvctrp,hpsi,hpsit)
!
!        deallocate(hpsi)
!
!       hamovr(jorb,iorb,3)=+psit(k,jorb)*hpsit(k,iorb)
!       hamovr(jorb,iorb,4)=+psit(k,jorb)* psit(k,iorb)
!      call DGEMM('T','N',norbe,norbe,nvctrp,1.d0,psit,nvctrp,hpsit,nvctrp,0.d0,hamovr(1,1,3),norbe)
!      call DGEMM('T','N',norbe,norbe,nvctrp,1.d0,psit,nvctrp, psit,nvctrp,0.d0,hamovr(1,1,4),norbe)
!        deallocate(hpsit)
!
!        call MPI_ALLREDUCE (hamovr(1,1,3),hamovr(1,1,1),2*norbe**2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!
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
!
!        n_lp=5000
!        allocate(work_lp(n_lp),evale(norbe))
!        call  DSYGV(1,'V','U',norbe,hamovr(1,1,1),norbe,hamovr(1,1,2),norbe,evale, work_lp, n_lp, info )
!        if (info.ne.0) write(*,*) 'DSYGV ERROR',info
!        if (iproc.eq.0) then
!        do iorb=1,norbe
!        write(*,*) 'evale(',iorb,')=',evale(iorb)
!        enddo
!        endif
!        cprec=evale(1)
!        deallocate(work_lp,evale)
!
!        write(79,'(a40,i10)') 'words for ppsit ',nvctrp*norbp*nproc
!        allocate(ppsit(nvctrp,norbp*nproc))
!        write(79,*) 'allocation done'
!
!! ppsit(k,iorb)=+psit(k,jorb)*hamovr(jorb,iorb,1)
!      call DGEMM('N','N',nvctrp,norb,norbe,1.d0,psit,nvctrp,hamovr,norbe,0.d0,ppsit,nvctrp)
!
!       call  untransallwaves(iproc,nproc,norb,norbp,nvctr_c,nvctr_f,nvctrp,ppsit,ppsi)
!
!        deallocate(psit,ppsit)
!
!        if (parallel) call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!  else !serial case
!!       hamovr(jorb,iorb,3)=+psi(k,jorb)*hpsi(k,iorb)
!      call DGEMM('T','N',norbe,norbe,nvctrp,1.d0,psi,nvctrp,hpsi,nvctrp,0.d0,hamovr(1,1,1),norbe)
!      call DGEMM('T','N',norbe,norbe,nvctrp,1.d0,psi,nvctrp, psi,nvctrp,0.d0,hamovr(1,1,2),norbe)
!        deallocate(hpsi)
!
!! calculate  KS orbitals
!        write(*,*) 'KS Hamiltonian'
!        do iorb=1,norbe
!        write(*,'(10(1x,e10.3))') (hamovr(iorb,jorb,1),jorb=1,norbe)
!        enddo
!        write(*,*) 'Overlap'
!        do iorb=1,norbe
!        write(*,'(10(1x,e10.3))') (hamovr(iorb,jorb,2),jorb=1,norbe)
!        enddo
!
!        n_lp=5000
!        allocate(work_lp(n_lp),evale(norbe))
!        call  DSYGV(1,'V','U',norbe,hamovr(1,1,1),norbe,hamovr(1,1,2),norbe,evale, work_lp, !n_lp, info )
!        if (info.ne.0) write(*,*) 'DSYGV ERROR',info
!        if (iproc.eq.0) then
!        do iorb=1,norbe
!        write(*,*) 'evale(',iorb,')=',evale(iorb)
!        enddo
!        endif
!        cprec=evale(1)
!        deallocate(work_lp,evale)
!
!! ppsi(k,iorb)=+psi(k,jorb)*hamovr(jorb,iorb,1)
!        call DGEMM('N','N',nvctrp,norb,norbe,1.d0,psi,nvctrp,hamovr,norbe,0.d0,ppsi,nvctrp)
!        deallocate(psi)
!
!  endif

!        deallocate(hamovr,occupe)

	return
	end subroutine input_wf_diag
