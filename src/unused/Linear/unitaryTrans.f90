!!!!subroutine MLWF(iproc, nproc, lzd, orbs, at, op, comon, mad, rxyz, nit, kernel, &
!!!!           newgradient, confdatarr, hx, lphi, Umat)
!!!!use module_base
!!!!use module_types
!!!!use module_interfaces, exceptThisOne => MLWF
!!!!implicit none
!!!!
!!!!! Calling arguments
!!!!integer,intent(in):: iproc, nproc, nit
!!!!type(local_zone_descriptors),intent(in):: lzd
!!!!type(orbitals_data),intent(in):: orbs
!!!!type(atoms_data),intent(in):: at
!!!!type(overlapParameters),intent(inout):: op
!!!!type(p2pComms),intent(inout):: comon
!!!!type(matrixDescriptors),intent(in):: mad
!!!!real(8),dimension(3,at%nat),intent(in):: rxyz
!!!!real(8),dimension(orbs%norb,orbs%norb),intent(in):: kernel
!!!!logical,intent(in):: newgradient
!!!!real(8),intent(in):: hx
!!!!type(confpot_data),dimension(orbs%norbp),intent(in):: confdatarr
!!!!real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(inout):: lphi
!!!!real(8),dimension(orbs%norb,orbs%norb),intent(out):: Umat
!!!!
!!!!! Local variables
!!!!integer:: it, info, lwork, k, istat, iorb, jorb, iall, ierr, ist, jst, ilrold, ncount, jjorb, iiorb, ilr, lorb, jlr
!!!!integer:: nlocregOnMPI, jlrold, jj, ii
!!!!real(8):: trace, lstep, dfactorial, energyconf_trial, energyconf_0, energyconf_der0, lstep_optimal, ddot
!!!!real(8):: tt1, tt2, tt3, tt4, tt5, tt
!!!!real(8):: t1, t2, t1_tot, t2_tot, omega
!!!!real(8):: time_convol, time_commun, time_lincomb, time_linalg, time_matrixmodification, time_exponential, time_tot
!!!!real(8):: time_matrixelements
!!!!complex(8):: ttc
!!!!real(8),dimension(:,:),allocatable:: gmat, hamtrans, ttmat, Kmat, X, Y, Z, Xd, Yd, Zd, commutX, commutY, commutZ
!!!!real(8),dimension(:,:),allocatable:: Xprime, Yprime, Zprime, Xprimesquare, Yprimesquare, Zprimesquare
!!!!real(8),dimension(:,:,:),allocatable:: potmat, potmatsmall
!!!!complex(8),dimension(:,:),allocatable:: gmatc, omatc
!!!!complex(8),dimension(:,:,:),allocatable:: tempmatc
!!!!complex(8),dimension(:),allocatable:: work, expD_cmplx
!!!!real(8),dimension(:),allocatable:: rwork, eval, lphiovrlp, lvphi, recvbuf, lxphi, lyphi, lzphi
!!!!real(8),dimension(:,:,:),allocatable:: tempmat3
!!!!character(len=*),parameter:: subname='unitary_optimization'
!!!!type(p2pComms):: comon_local
!!!!
!!!!! Quick return if possible
!!!!if(nit==0) return
!!!!
!!!!allocate(gmat(orbs%norb,orbs%norb), stat=istat)
!!!!call memocc(istat, gmat, 'gmat', subname)
!!!!allocate(hamtrans(orbs%norb,orbs%norb), stat=istat)
!!!!call memocc(istat, hamtrans, 'hamtrans', subname)
!!!!allocate(gmatc(orbs%norb,orbs%norb), stat=istat)
!!!!!call memocc(istat, gmatc, 'gmatc', subname)
!!!!allocate(omatc(orbs%norb,orbs%norb), stat=istat)
!!!!!call memocc(istat, omatc, 'omatc', subname)
!!!!allocate(tempmat3(orbs%norb,orbs%norb,3), stat=istat)
!!!!call memocc(istat, tempmat3, 'tempmat3', subname)
!!!!allocate(eval(orbs%norb), stat=istat)
!!!!call memocc(istat, eval, 'eval', subname)
!!!!allocate(expD_cmplx(orbs%norb), stat=istat)
!!!!call memocc(istat, expD_cmplx, 'expD_cmplx', subname)
!!!!allocate(tempmatc(orbs%norb,orbs%norb,2), stat=istat)
!!!!!call memocc(istat, tempmatc, 'tempmatc', subname)
!!!!allocate(lphiovrlp(op%ndim_lphiovrlp), stat=istat)
!!!!call memocc(istat, lphiovrlp, 'lphiovrlp',subname)
!!!!allocate(lvphi(max(orbs%npsidim_orbs,orbs%npsidim_comp)), stat=istat)
!!!!call memocc(istat, lvphi, 'lvphi', subname)
!!!!allocate(Kmat(orbs%norb,orbs%norb), stat=istat)
!!!!call memocc(istat, Kmat, 'Kmat', subname)
!!!!allocate(recvbuf(comon%nrecvbuf), stat=istat)
!!!!call memocc(istat, recvbuf, 'recvbuf', subname)
!!!!
!!!!allocate(potmat(orbs%norb,orbs%norb,at%nat), stat=istat)
!!!!call memocc(istat, potmat, 'potmat', subname)
!!!!
!!!!allocate(lxphi(max(orbs%npsidim_orbs,orbs%npsidim_comp)), stat=istat)
!!!!call memocc(istat, lxphi, 'lxphi', subname)
!!!!allocate(lyphi(max(orbs%npsidim_orbs,orbs%npsidim_comp)), stat=istat)
!!!!call memocc(istat, lyphi, 'lyphi', subname)
!!!!allocate(lzphi(max(orbs%npsidim_orbs,orbs%npsidim_comp)), stat=istat)
!!!!call memocc(istat, lzphi, 'lzphi', subname)
!!!!
!!!!allocate(X(orbs%norb,orbs%norb), stat=istat)
!!!!call memocc(istat, X, 'X', subname)
!!!!allocate(Y(orbs%norb,orbs%norb), stat=istat)
!!!!call memocc(istat, Y, 'Y', subname)
!!!!allocate(Z(orbs%norb,orbs%norb), stat=istat)
!!!!call memocc(istat, Z, 'Z', subname)
!!!!allocate(Xd(orbs%norb,orbs%norb), stat=istat)
!!!!call memocc(istat, Xd, 'Xd', subname)
!!!!allocate(Yd(orbs%norb,orbs%norb), stat=istat)
!!!!call memocc(istat, Yd, 'Yd', subname)
!!!!allocate(Zd(orbs%norb,orbs%norb), stat=istat)
!!!!call memocc(istat, Zd, 'Zd', subname)
!!!!allocate(Xprime(orbs%norb,orbs%norb), stat=istat)
!!!!call memocc(istat, Xprime, 'Xprime', subname)
!!!!allocate(Yprime(orbs%norb,orbs%norb), stat=istat)
!!!!call memocc(istat, Yprime, 'Yprime', subname)
!!!!allocate(Zprime(orbs%norb,orbs%norb), stat=istat)
!!!!call memocc(istat, Zprime, 'Zprime', subname)
!!!!allocate(Xprimesquare(orbs%norb,orbs%norb), stat=istat)
!!!!call memocc(istat, Xprimesquare, 'Xprimesquare', subname)
!!!!allocate(Yprimesquare(orbs%norb,orbs%norb), stat=istat)
!!!!call memocc(istat, Yprimesquare, 'Yprimesquare', subname)
!!!!allocate(Zprimesquare(orbs%norb,orbs%norb), stat=istat)
!!!!call memocc(istat, Zprimesquare, 'Zprimesquare', subname)
!!!!allocate(commutX(orbs%norb,orbs%norb), stat=istat)
!!!!call memocc(istat, commutX, 'commutX', subname)
!!!!allocate(commutY(orbs%norb,orbs%norb), stat=istat)
!!!!call memocc(istat, commutY, 'commutY', subname)
!!!!allocate(commutZ(orbs%norb,orbs%norb), stat=istat)
!!!!call memocc(istat, commutZ, 'commutZ', subname)
!!!!
!!!!
!!!!! Count how many locregs each process handles
!!!!ilrold=-1
!!!!nlocregOnMPI=0
!!!!do iorb=1,orbs%norbp
!!!!    iiorb=orbs%isorb+iorb
!!!!    ilr=orbs%inwhichlocreg(iiorb)
!!!!    !if(ilr>ilrold) then
!!!!    if(ilr/=ilrold) then
!!!!        nlocregOnMPI=nlocregOnMPI+1
!!!!    end if
!!!!    ilrold=ilr
!!!!end do
!!!!allocate(potmatsmall(orbs%norb,orbs%norb,nlocregOnMPI), stat=istat)
!!!!call memocc(istat, potmatsmall, 'potmatsmall', subname)
!!!!
!!!!
!!!!
!!!!  call allocateSendBufferOrtho(comon, subname)
!!!!  call allocateRecvBufferOrtho(comon, subname)
!!!!  ! Extract the overlap region from the orbitals phi and store them in comon%sendBuf.
!!!!  !if(nit>0) then
!!!!  !    call extractOrbital3(iproc, nproc, orbs, orbs%npsidim, orbs%inWhichLocreg, lzd, op, lphi, comon%nsendBuf, comon%sendBuf)
!!!!  !    call postCommsOverlapNew(iproc, nproc, orbs, op, lzd, lphi, comon, tt1, tt2)
!!!!  !end if
!!!!
!!!!  energyconf_trial=0.d0 !just to initialize this variable and make the compiler happy
!!!!  lstep=0.d0 !just to initialize this variable and make the compiler happy
!!!!  lstep_optimal=0.d0 !just to initialize this variable and make the compiler happy
!!!!  energyconf_der0=0.d0 !just to initialize this variable and make the compiler happy
!!!!
!!!!  time_convol=0.d0
!!!!  time_lincomb=0.d0
!!!!  time_commun=0.d0
!!!!  time_linalg=0.d0
!!!!  time_exponential=0.d0
!!!!  time_matrixmodification=0.d0
!!!!  time_matrixElements=0.d0
!!!!  t1_tot=mpi_wtime()
!!!!
!!!!  ! Initialize Umat
!!!!  do iorb=1,orbs%norb
!!!!      do jorb=1,orbs%norb
!!!!          if(jorb==iorb) then
!!!!              Umat(jorb,iorb)=1.d0
!!!!          else
!!!!              Umat(jorb,iorb)=0.d0
!!!!          end if
!!!!      end do
!!!!  end do
!!!!
!!!!  innerLoop: do it=1,nit
!!!!
!!!!  !write(*,*) '1: iproc, associated(comon%recvbuf)', iproc, associated(comon%recvbuf)
!!!!
!!!!      t1=mpi_wtime()
!!!!      !!call apply_orbitaldependent_potential(iproc, nproc, at, orbs, lzd, rxyz, &
!!!!      !!     confdatarr, hx, lphi, -1, lvphi)
!!!!      call apply_position_operators(iproc, nproc, orbs, lzd, hx, hx, hx, lphi, lxphi, lyphi, lzphi)
!!!!
!!!!
!!!!      t2=mpi_wtime()
!!!!      time_convol=time_convol+t2-t1
!!!!
!!!!      t1=mpi_wtime()
!!!!      !allocate(ttmat(lin%orbs%norb,lin%orbs%norb))
!!!!      !call collectnew(iproc, nproc, comon, lin%mad,lin%op, lin%orbs, input, lin%lzd, comon%nsendbuf, &
!!!!      !     comon%sendbuf, comon%nrecvbuf, comon%recvbuf, ttmat, tt3, tt4, tt5)
!!!!      !deallocate(ttmat)
!!!!      !write(*,*) '2: iproc, associated(comon%recvbuf)', iproc, associated(comon%recvbuf)
!!!!      t2=mpi_wtime()
!!!!      time_commun=time_commun+t2-t1
!!!!
!!!!      t1=mpi_wtime()
!!!!      !call getMatrixElements2(iproc, nproc, lin%lzd, lin%lb%orbs, lin%lb%op, lin%lb%comon, lphi, lvphi, lin%mad, Kmat)
!!!!      !call deallocateRecvBufferOrtho(comon, subname)
!!!!      !call deallocateSendBufferOrtho(comon, subname)
!!!!      !!call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lvphi, mad, Kmat)
!!!!
!!!!      ! Build the matrices X, Y, Z
!!!!      call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lxphi, mad, X)
!!!!      call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lyphi, mad, Y)
!!!!      call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lzphi, mad, Z)
!!!!
!!!!      ! Build the matrices Xd, Yd, Zd
!!!!      do iorb=1,orbs%norb
!!!!          do jorb=1,orbs%norb
!!!!              if(jorb==iorb) then
!!!!                  Xd(jorb,iorb)=X(jorb,iorb)
!!!!                  Yd(jorb,iorb)=Y(jorb,iorb)
!!!!                  Zd(jorb,iorb)=Z(jorb,iorb)
!!!!              else
!!!!                  Xd(jorb,iorb)=0.d0
!!!!                  Yd(jorb,iorb)=0.d0
!!!!                  Zd(jorb,iorb)=0.d0
!!!!              end if
!!!!          end do
!!!!      end do
!!!!
!!!!      ! Build the matrices Xprime, Yprime, Zprime
!!!!      Xprime=X-Xd
!!!!      Yprime=Y-Yd
!!!!      Zprime=Z-Zd
!!!!
!!!!      ! Calculate value of Omega
!!!!      call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, Xprime(1,1), orbs%norb, &
!!!!           Xprime(1,1), orbs%norb, 0.d0, Xprimesquare, orbs%norb)
!!!!      call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, Yprime(1,1), orbs%norb, &
!!!!           Yprime(1,1), orbs%norb, 0.d0, Yprimesquare, orbs%norb)
!!!!      call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, Zprime(1,1), orbs%norb, &
!!!!           Zprime(1,1), orbs%norb, 0.d0, Zprimesquare, orbs%norb)
!!!!      omega=0.d0
!!!!      do iorb=1,orbs%norb
!!!!          omega=omega+Xprimesquare(iorb,iorb)+Yprimesquare(iorb,iorb)+Zprimesquare(iorb,iorb)
!!!!      end do
!!!!      if(iproc==0) write(*,'(a,i7,2es16.7)') 'it, omega, lstep', it, omega, lstep
!!!!
!!!!      call commutator(orbs%norb, Xprime, Xd, commutX)
!!!!      call commutator(orbs%norb, Yprime, Yd, commutY)
!!!!      call commutator(orbs%norb, Zprime, Zd, commutZ)
!!!!
!!!!      gmat=-2.d0*(commutX+commutY+commutZ)
!!!!
!!!!
!!!!
!!!!      !do iorb=1,orbs%norb
!!!!      !    do jorb=1,orbs%norb
!!!!      !        if(iproc==0) write(66,*) iorb,jorb,Kmat(jorb,iorb)
!!!!      !    end do
!!!!      !end do
!!!!      !call allocateSendBufferOrtho(comon, subname)
!!!!      !call allocateRecvBufferOrtho(comon, subname)
!!!!      !write(*,*) '3: iproc, associated(comon%recvbuf)', iproc, associated(comon%recvbuf)
!!!!      t2=mpi_wtime()
!!!!      time_matrixElements=time_matrixElements+t2-t1
!!!!
!!!!
!!!!      !!if(newgradient) then
!!!!      !!    call get_potential_matrices(iproc, nproc, at, orbs, lzd, op, comon, mad, rxyz, &
!!!!      !!         confdatarr, hx, lphi, potmat)
!!!!      !!    !call get_potential_matrices_new(iproc, nproc, lin, at, input, orbs, lzd, op, comon, rxyz, lphi, &
!!!!      !!    !     nlocregOnMPI, potmatsmall)
!!!!      !!end if
!!!!
!!!!
!!!!      !!if(.not.newgradient) then
!!!!          !energyconf_0=ddot(orbs%npsidim, lphi(1), 1, lvphi(1), 1)
!!!!          !call mpiallred(energyconf_0, 1, mpi_sum, mpi_comm_world, ierr)
!!!!          !!$$energyconf_0=0.d0
!!!!          !!$$do iorb=1,orbs%norb
!!!!          !!$$    energyconf_0 = energyconf_0 + Kmat(iorb,iorb)
!!!!          !!$$end do
!!!!      !!else
!!!!      !!    energyconf_0=0.d0
!!!!      !!    do iorb=1,orbs%norb
!!!!      !!        do jorb=1,orbs%norb
!!!!      !!            energyconf_0 = energyconf_0 + kernel(jorb,iorb)*Kmat(jorb,iorb)
!!!!      !!            !energyconf_0 = energyconf_0 + kernel(jorb,iorb)*Kmat(iorb,jorb)
!!!!      !!        end do
!!!!      !!    end do
!!!!      !!end if
!!!!      !!$$if(iproc==0) write(*,'(a,i6,3es20.10,2es17.7)') &
!!!!      !!$$             'it, energyconf_0, energyvonf_trial, energyconf_der0, lstep, lstep_optimal', &
!!!!      !!$$             it, energyconf_0, energyconf_trial, energyconf_der0, lstep, lstep_optimal
!!!!
!!!!      t1=mpi_wtime()
!!!!      !!if(.not.newgradient) then
!!!!          ! Construct antisymmtric matrix Gmat
!!!!          !!$$do iorb=1,orbs%norb
!!!!          !!$$    do jorb=1,orbs%norb
!!!!          !!$$        gmat(jorb,iorb)=2.d0*(Kmat(jorb,iorb)-Kmat(iorb,jorb))
!!!!          !!$$    end do
!!!!          !!$$end do 
!!!!      !!else
!!!!      !!    !!! THIS IS THE OLD VERSION #############################################################################
!!!!      !!    do iorb=1,orbs%norb
!!!!      !!        ilr=orbs%inwhichlocreg(iorb)
!!!!      !!        do jorb=1,orbs%norb
!!!!      !!            jlr=orbs%inwhichlocreg(jorb)
!!!!      !!            tt=0.d0
!!!!      !!            do lorb=1,orbs%norb
!!!!      !!                tt = tt + kernel(jorb,lorb)*Kmat(lorb,iorb) - kernel(iorb,lorb)*Kmat(lorb,jorb) + &
!!!!      !!                          kernel(jorb,lorb)*potmat(lorb,iorb,jlr) - kernel(iorb,lorb)*potmat(lorb,jorb,ilr)
!!!!      !!            end do
!!!!      !!            gmat(jorb,iorb)=-tt
!!!!      !!            !if(iproc==0) then
!!!!      !!            !    write(77,*) iorb, jorb, gmat(jorb,iorb)
!!!!      !!            !end if
!!!!      !!        end do
!!!!      !!    end do 
!!!!      !!    ! ########################################################################################################
!!!!      !!    !!! THIS IS THE NEW VERSION
!!!!      !!    !!gmat=0.d0
!!!!      !!    !!ii=0
!!!!      !!    !!ilrold=-1
!!!!      !!    !!do iorb=1,orbs%norbp
!!!!      !!    !!    iiorb=orbs%isorb+iorb
!!!!      !!    !!    ilr=orbs%inwhichlocreg(iiorb)
!!!!      !!    !!    if(ilr>ilrold) then
!!!!      !!    !!        ii=ii+1
!!!!      !!    !!    end if
!!!!      !!    !!    do jorb=1,orbs%norb
!!!!      !!    !!        jlr=orbs%inwhichlocreg(jorb)
!!!!      !!    !!        tt=0.d0
!!!!      !!    !!        do lorb=1,orbs%norb
!!!!      !!    !!            !tt = tt + kernel(jorb,lorb)*Kmat(lorb,iiorb) - kernel(iiorb,lorb)*Kmat(lorb,jorb) + &
!!!!      !!    !!            !          - kernel(iiorb,lorb)*potmat(lorb,jorb,ilr)
!!!!      !!    !!            tt = tt + kernel(jorb,lorb)*Kmat(lorb,iiorb) - kernel(iiorb,lorb)*Kmat(lorb,jorb) + &
!!!!      !!    !!                      - kernel(iiorb,lorb)*potmatsmall(lorb,jorb,ii)
!!!!      !!    !!        end do
!!!!      !!    !!        gmat(jorb,iiorb)=-tt
!!!!      !!    !!    end do
!!!!      !!    !!    ilrold=ilr
!!!!      !!    !!end do 
!!!!      !!    !!do iorb=1,orbs%norb
!!!!      !!    !!    ilr=orbs%inwhichlocreg(iorb)
!!!!      !!    !!    jlrold=-1
!!!!      !!    !!    jj=0
!!!!      !!    !!    do jorb=1,orbs%norbp
!!!!      !!    !!        jjorb=orbs%isorb+jorb
!!!!      !!    !!        jlr=orbs%inwhichlocreg(jjorb)
!!!!      !!    !!        if(jlr>jlrold) then
!!!!      !!    !!            jj=jj+1
!!!!      !!    !!        end if
!!!!      !!    !!        tt=0.d0
!!!!      !!    !!        do lorb=1,orbs%norb
!!!!      !!    !!            !tt = tt + kernel(jjorb,lorb)*potmat(lorb,iorb,jlr)
!!!!      !!    !!            tt = tt + kernel(jjorb,lorb)*potmatsmall(lorb,iorb,jj)
!!!!      !!    !!        end do
!!!!      !!    !!        gmat(jjorb,iorb)=gmat(jjorb,iorb)-tt
!!!!      !!    !!        jlrold=jlr
!!!!      !!    !!    end do
!!!!      !!    !!end do 
!!!!      !!    !!call mpiallred(gmat(1,1), orbs%norb**2, mpi_sum, mpi_comm_world, ierr)
!!!!      !!    !!do iorb=1,orbs%norb
!!!!      !!    !!    do jorb=1,orbs%norb
!!!!      !!    !!        if(iproc==0) then
!!!!      !!    !!            write(77,*) iorb, jorb, gmat(jorb,iorb)
!!!!      !!    !!        end if
!!!!      !!    !!    end do
!!!!      !!    !!end do
!!!!      !!end if
!!!!      t2=mpi_wtime()
!!!!      time_matrixmodification=time_matrixmodification+t2-t1
!!!!
!!!!      t1=mpi_wtime()
!!!!      !Build the complex matrix -iGmat
!!!!      do iorb=1,orbs%norb
!!!!          do jorb=1,orbs%norb
!!!!              gmatc(jorb,iorb)=cmplx(0.d0,-gmat(jorb,iorb),kind=8)
!!!!              if(iproc==0) write(999,'(a,2i8,2es16.7)') 'iorb, jorb, gmatc(jorb,iorb)', iorb, jorb, gmatc(jorb,iorb)
!!!!          end do
!!!!      end do 
!!!!      t2=mpi_wtime()
!!!!      time_matrixmodification=time_matrixmodification+t2-t1
!!!!
!!!!
!!!!
!!!!      ! Diagonalize Gmatc
!!!!      t1=mpi_wtime()
!!!!      lwork=10*orbs%norb
!!!!      allocate(work(lwork), stat=istat) ! factor of 2 since it is assumed to be complex
!!!!      allocate(rwork(lwork), stat=istat)
!!!!      call zheev('v', 'l', orbs%norb, gmatc(1,1), orbs%norb, eval(1), work, lwork, rwork, info)
!!!!      if(info/=0) stop 'ERROR in zheev'
!!!!      deallocate(work)
!!!!      deallocate(rwork)
!!!!      t2=mpi_wtime()
!!!!      time_linalg=time_linalg+t2-t1
!!!!
!!!!
!!!!      !!$$! Calculate step size
!!!!      !!$$if(it==1) then
!!!!      !!$$    !!if(.not.newgradient) then
!!!!      !!$$        lstep=5.d-2/(maxval(eval))
!!!!      !!$$    !!else
!!!!      !!$$    !!    lstep=1.d-4/(maxval(eval))
!!!!      !!$$    !!end if
!!!!      !!$$else
!!!!      !!$$    lstep=2.d0*lstep_optimal
!!!!      !!$$    !lstep=1.d-3/(maxval(eval))
!!!!      !!$$end if
!!!!      lstep=1.d-2/(maxval(eval))
!!!!
!!!!      t1=mpi_wtime()
!!!!      ! Calculate exp(-i*l*D) (with D diagonal matrix of eigenvalues).
!!!!      ! This is also a diagonal matrix, so only calculate the diagonal part.
!!!!      do iorb=1,orbs%norb
!!!!         ttc=cmplx(0.d0,-lstep*eval(iorb),kind=8)
!!!!         expD_cmplx(iorb)=(0.d0,0.d0)
!!!!          do k=0,50
!!!!              expD_cmplx(iorb)=expD_cmplx(iorb)+ttc**k/dfactorial(k)
!!!!          end do
!!!!      end do
!!!!      t2=mpi_wtime()
!!!!      time_exponential=time_exponential+t2-t1
!!!!
!!!!      t1=mpi_wtime()
!!!!      do iorb=1,orbs%norb
!!!!          do jorb=1,orbs%norb
!!!!              if(iorb==jorb) then
!!!!                  tempmatc(jorb,iorb,1)=expD_cmplx(iorb)
!!!!              else
!!!!                  tempmatc(jorb,iorb,1)=cmplx(0.d0,0.d0,kind=8)
!!!!              end if
!!!!          end do
!!!!      end do
!!!!      t2=mpi_wtime()
!!!!      time_matrixmodification=time_matrixmodification+t2-t1
!!!!
!!!!      t1=mpi_wtime()
!!!!      call zgemm('n', 'c', orbs%norb, orbs%norb, orbs%norb, (1.d0,0.d0), tempmatc(1,1,1), orbs%norb, &
!!!!           gmatc(1,1), orbs%norb, (0.d0,0.d0), tempmatc(1,1,2), orbs%norb)
!!!!      call zgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, (1.d0,0.d0), gmatc(1,1), orbs%norb, &
!!!!           tempmatc(1,1,2), orbs%norb, (0.d0,0.d0), omatc(1,1), orbs%norb)
!!!!      t2=mpi_wtime()
!!!!      time_linalg=time_linalg+t2-t1
!!!!
!!!!      t1=mpi_wtime()
!!!!      ! Build new lphi
!!!!      do iorb=1,orbs%norb
!!!!          do jorb=1,orbs%norb
!!!!              tempmat3(jorb,iorb,1)=real(omatc(jorb,iorb))
!!!!          end do
!!!!      end do
!!!!      t2=mpi_wtime()
!!!!      time_matrixmodification=time_matrixmodification+t2-t1
!!!!
!!!!      ! Update Umat
!!!!      call dcopy(orbs%norb**2, Umat(1,1), 1, tempmat3(1,1,2), 1)
!!!!      call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,2), orbs%norb, &
!!!!           tempmat3(1,1,1), orbs%norb, 0.d0, Umat(1,1), orbs%norb)
!!!!
!!!!      t1=mpi_wtime()
!!!!      !write(*,*) '5: iproc, associated(comon%recvbuf)', iproc, associated(comon%recvbuf)
!!!!      call build_new_linear_combinations(iproc, nproc, lzd, orbs, op, comon%nrecvbuf, &
!!!!           comon%recvbuf, tempmat3(1,1,1), .true., lphi)
!!!!      t2=mpi_wtime()
!!!!      time_lincomb=time_lincomb+t2-t1
!!!!
!!!!      !!$$t1=mpi_wtime()
!!!!      !!$$call apply_orbitaldependent_potential(iproc, nproc, at, orbs, lzd, rxyz, &
!!!!      !!$$     confdatarr, hx, lphi, -1, lvphi)
!!!!      !!$$t2=mpi_wtime()
!!!!      !!$$time_convol=time_convol+t2-t1
!!!!
!!!!      !!$$!!if(.not.newgradient) then
!!!!      !!$$    energyconf_trial=ddot(max(orbs%npsidim_orbs,orbs%npsidim_comp), lphi(1), 1, lvphi(1), 1)
!!!!      !!$$    call mpiallred(energyconf_trial, 1, mpi_sum, mpi_comm_world, ierr)
!!!!      !!$$!!else
!!!!
!!!!      !!$$!!    call dcopy(comon%nrecvbuf, comon%recvbuf, 1, recvbuf, 1)
!!!!      !!$$!!    !call extractOrbital3(iproc, nproc, orbs, orbs%npsidim, orbs%inWhichLocreg, lzd, op, lphi, comon%nsendBuf, comon%sendBuf)
!!!!      !!$$!!    !call postCommsOverlapNew(iproc, nproc, orbs, op, lzd, lphi, comon, tt1, tt2)
!!!!      !!$$!!    !allocate(ttmat(lin%orbs%norb,lin%orbs%norb))
!!!!      !!$$!!    !call collectnew(iproc, nproc, comon, lin%mad,lin%op, lin%orbs, input, lin%lzd, comon%nsendbuf, &
!!!!      !!$$!!    !     comon%sendbuf, comon%nrecvbuf, comon%recvbuf, ttmat, tt3, tt4, tt5)
!!!!      !!$$!!    !deallocate(ttmat)
!!!!      !!$$!!    call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lvphi, mad, Kmat)
!!!!      !!$$!!    call dcopy(comon%nrecvbuf, recvbuf, 1, comon%recvbuf, 1)
!!!!
!!!!      !!$$!!    energyconf_trial=0.d0
!!!!      !!$$!!    do iorb=1,orbs%norb
!!!!      !!$$!!        do jorb=1,orbs%norb
!!!!      !!$$!!            energyconf_trial = energyconf_trial + kernel(jorb,iorb)*Kmat(jorb,iorb)
!!!!      !!$$!!            !energyconf_trial = energyconf_trial + kernel(jorb,iorb)*Kmat(iorb,jorb)
!!!!      !!$$!!        end do
!!!!      !!$$!!    end do
!!!!      !!$$!!end if
!!!!
!!!!      !!$$! Calculate the gradient of the confinement
!!!!      !!$$energyconf_der0=0.d0
!!!!      !!$$do iorb=1,orbs%norb
!!!!      !!$$    do jorb=1,orbs%norb
!!!!      !!$$        energyconf_der0=energyconf_der0+gmat(jorb,iorb)**2
!!!!      !!$$    end do
!!!!      !!$$end do
!!!!      !!$$energyconf_der0=-.5d0*energyconf_der0
!!!!
!!!!      !!$$! Calculate optimal step size
!!!!      !!$$lstep_optimal = -energyconf_der0*lstep**2/(2.d0*(energyconf_trial-energyconf_0-lstep*energyconf_der0))
!!!!      !!$$!!if(.not.newgradient) then
!!!!      !!$$    lstep_optimal=min(lstep_optimal,lstep)
!!!!      !!$$!!else
!!!!      !!$$!!    if(lstep_optimal<0) then
!!!!      !!$$!!        lstep_optimal=lstep
!!!!      !!$$!!    else
!!!!      !!$$!!        lstep_optimal=min(lstep_optimal,lstep)
!!!!      !!$$!!    end if
!!!!      !!$$!!end if
!!!!
!!!!      !!$$t1=mpi_wtime()
!!!!      !!$$! Calculate exp(-i*l*D) (with D diagonal matrix of eigenvalues).
!!!!      !!$$! This is also a diagonal matrix, so only calculate the diagonal part.
!!!!      !!$$do iorb=1,orbs%norb
!!!!      !!$$   ttc=cmplx(0.d0,-lstep_optimal*eval(iorb),kind=8)
!!!!      !!$$   expD_cmplx(iorb)=(0.d0,0.d0)
!!!!      !!$$    do k=0,50
!!!!      !!$$        expD_cmplx(iorb)=expD_cmplx(iorb)+ttc**k/dfactorial(k)
!!!!      !!$$    end do
!!!!      !!$$end do
!!!!      !!$$t2=mpi_wtime()
!!!!      !!$$time_exponential=time_exponential+t2-t1
!!!!
!!!!      !!$$t1=mpi_wtime()
!!!!      !!$$do iorb=1,orbs%norb
!!!!      !!$$    do jorb=1,orbs%norb
!!!!      !!$$        if(iorb==jorb) then
!!!!      !!$$            tempmatc(jorb,iorb,1)=expD_cmplx(iorb)
!!!!      !!$$        else
!!!!      !!$$            tempmatc(jorb,iorb,1)=cmplx(0.d0,0.d0,kind=8)
!!!!      !!$$        end if
!!!!      !!$$    end do
!!!!      !!$$end do
!!!!      !!$$t2=mpi_wtime()
!!!!      !!$$time_matrixmodification=time_matrixmodification+t2-t1
!!!!
!!!!      !!$$t1=mpi_wtime()
!!!!      !!$$call zgemm('n', 'c', orbs%norb, orbs%norb, orbs%norb, (1.d0,0.d0), tempmatc(1,1,1), orbs%norb, &
!!!!      !!$$     gmatc(1,1), orbs%norb, (0.d0,0.d0), tempmatc(1,1,2), orbs%norb)
!!!!      !!$$call zgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, (1.d0,0.d0), gmatc(1,1), orbs%norb, &
!!!!      !!$$     tempmatc(1,1,2), orbs%norb, (0.d0,0.d0), omatc(1,1), orbs%norb)
!!!!      !!$$t2=mpi_wtime()
!!!!      !!$$time_linalg=time_linalg+t2-t1
!!!!
!!!!
!!!!      !!$$! Build new lphi
!!!!      !!$$do iorb=1,orbs%norb
!!!!      !!$$    do jorb=1,orbs%norb
!!!!      !!$$        tempmat3(jorb,iorb,1)=real(omatc(jorb,iorb))
!!!!      !!$$    end do
!!!!      !!$$end do
!!!!      !!$$t1=mpi_wtime()
!!!!      !!$$call build_new_linear_combinations(iproc, nproc, lzd, orbs, op, comon%nrecvbuf, &
!!!!      !!$$     comon%recvbuf, tempmat3(1,1,1), .true., lphi)
!!!!      !!$$t2=mpi_wtime()
!!!!      !!$$time_lincomb=time_lincomb+t2-t1
!!!!
!!!!
!!!!      !!$$!if(it<nit) then
!!!!      !!$$!    call extractOrbital3(iproc, nproc, orbs, orbs%npsidim, orbs%inWhichLocreg, lzd, op, lphi, comon%nsendBuf, comon%sendBuf)
!!!!      !!$$!    call postCommsOverlapNew(iproc, nproc, orbs, op, lzd, lphi, comon, tt1, tt2)
!!!!      !!$$!end if
!!!!
!!!!
!!!!  end do innerLoop
!!!!
!!!!  t2_tot=mpi_wtime()
!!!!  time_tot=t2_tot-t1_tot
!!!!  call mpiallred(time_convol, 1, mpi_max, mpi_comm_world, ierr)
!!!!  call mpiallred(time_commun, 1, mpi_max, mpi_comm_world, ierr)
!!!!  call mpiallred(time_lincomb, 1, mpi_max, mpi_comm_world, ierr)
!!!!  call mpiallred(time_linalg, 1, mpi_max, mpi_comm_world, ierr)
!!!!  call mpiallred(time_matrixmodification, 1, mpi_max, mpi_comm_world, ierr)
!!!!  call mpiallred(time_exponential, 1, mpi_max, mpi_comm_world, ierr)
!!!!  call mpiallred(time_matrixelements, 1, mpi_max, mpi_comm_world, ierr)
!!!!  call mpiallred(time_tot, 1, mpi_max, mpi_comm_world, ierr)
!!!!  !time_convol=time_convol/dble(nproc)
!!!!  !time_commun=time_commun/dble(nproc)
!!!!  !time_lincomb=time_lincomb/dble(nproc)
!!!!  !time_linalg=time_linalg/dble(nproc)
!!!!  !time_matrixmodification=time_matrixmodification/dble(nproc)
!!!!  !time_exponential_=time_exponential/dble(nproc)
!!!!  !time_tot=time_tot/dble(nproc)
!!!!  !!if(iproc==0) then
!!!!  !!    write(*,'(a,es16.6)') 'total time: ',time_tot
!!!!  !!    write(*,'(a,es15.7,a,f5.2,a)') 'convolutions: ',time_convol,' (',time_convol/time_tot*100.d0,'%)'
!!!!  !!    write(*,'(a,es15.7,a,f5.2,a)') 'linear combinations: ',time_lincomb,' (',time_lincomb/time_tot*100.d0,'%)'
!!!!  !!    write(*,'(a,es15.7,a,f5.2,a)') 'communication: ',time_commun,' (',time_commun/time_tot*100.d0,'%)'
!!!!  !!    write(*,'(a,es15.7,a,f5.2,a)') 'linear algebra: ',time_linalg,' (',time_linalg/time_tot*100.d0,'%)'
!!!!  !!    write(*,'(a,es15.7,a,f5.2,a)') 'matrix modification: ',time_matrixmodification, &
!!!!  !!                                   ' (',time_matrixmodification/time_tot*100.d0,'%)'
!!!!  !!    write(*,'(a,es15.7,a,f5.2,a)') 'building exponential: ',time_exponential,' (',time_exponential/time_tot*100.d0,'%)'
!!!!  !!    write(*,'(a,es15.7,a,f5.2,a)') 'matrix elements ',time_matrixelements,' (',time_matrixelements/time_tot*100.d0,'%)'
!!!!  !!end if
!!!!
!!!!
!!!!  iall=-product(shape(gmat))*kind(gmat)
!!!!  deallocate(gmat, stat=istat)
!!!!  call memocc(istat, iall, 'gmat', subname)
!!!!  iall=-product(shape(gmatc))*kind(gmatc)
!!!!  deallocate(gmatc, stat=istat)
!!!!  !call memocc(istat, iall, 'gmatc', subname)
!!!!  iall=-product(shape(omatc))*kind(omatc)
!!!!  deallocate(omatc, stat=istat)
!!!!  !call memocc(istat, iall, 'omatc', subname)
!!!!  iall=-product(shape(tempmat3))*kind(tempmat3)
!!!!  deallocate(tempmat3, stat=istat)
!!!!  call memocc(istat, iall, 'tempmat3', subname)
!!!!  iall=-product(shape(eval))*kind(eval)
!!!!  deallocate(eval, stat=istat)
!!!!  call memocc(istat, iall, 'eval', subname)
!!!!  iall=-product(shape(expD_cmplx))*kind(expD_cmplx)
!!!!  deallocate(expD_cmplx, stat=istat)
!!!!  call memocc(istat, iall, 'expD_cmplx', subname)
!!!!  iall=-product(shape(tempmatc))*kind(tempmatc)
!!!!  deallocate(tempmatc, stat=istat)
!!!!  !call memocc(istat, iall, 'tempmatc', subname)
!!!!  iall=-product(shape(hamtrans))*kind(hamtrans)
!!!!  deallocate(hamtrans, stat=istat)
!!!!  call memocc(istat, iall, 'hamtrans', subname)
!!!!  iall=-product(shape(lphiovrlp))*kind(lphiovrlp)
!!!!  deallocate(lphiovrlp, stat=istat)
!!!!  call memocc(istat, iall, 'lphiovrlp', subname)
!!!!  iall=-product(shape(lvphi))*kind(lvphi)
!!!!  deallocate(lvphi, stat=istat)
!!!!  call memocc(istat, iall, 'lvphi', subname)
!!!!  iall=-product(shape(Kmat))*kind(Kmat)
!!!!  deallocate(Kmat, stat=istat)
!!!!  call memocc(istat, iall, 'Kmat', subname)
!!!!  iall=-product(shape(recvbuf))*kind(recvbuf)
!!!!  deallocate(recvbuf, stat=istat)
!!!!  call memocc(istat, iall, 'recvbuf', subname)
!!!!
!!!!  iall=-product(shape(potmat))*kind(potmat)
!!!!  deallocate(potmat, stat=istat)
!!!!  call memocc(istat, iall, 'potmat', subname)
!!!!  iall=-product(shape(potmatsmall))*kind(potmatsmall)
!!!!  deallocate(potmatsmall, stat=istat)
!!!!  call memocc(istat, iall, 'potmatsmall', subname)
!!!!
!!!!
!!!!  iall=-product(shape(lxphi))*kind(lxphi)
!!!!  deallocate(lxphi, stat=istat)
!!!!  call memocc(istat, iall, 'lxphi', subname)
!!!!  iall=-product(shape(lyphi))*kind(lyphi)
!!!!  deallocate(lyphi, stat=istat)
!!!!  call memocc(istat, iall, 'lyphi', subname)
!!!!  iall=-product(shape(lzphi))*kind(lzphi)
!!!!  deallocate(lzphi, stat=istat)
!!!!  call memocc(istat, iall, 'lzphi', subname)
!!!!
!!!!
!!!!
!!!!  iall=-product(shape(X))*kind(X)
!!!!  deallocate(X, stat=istat)
!!!!  call memocc(istat, iall, 'X', subname)
!!!!  iall=-product(shape(Y))*kind(Y)
!!!!  deallocate(Y, stat=istat)
!!!!  call memocc(istat, iall, 'Y', subname)
!!!!  iall=-product(shape(Z))*kind(Z)
!!!!  deallocate(Z, stat=istat)
!!!!  call memocc(istat, iall, 'Z', subname)
!!!!  iall=-product(shape(Xd))*kind(Xd)
!!!!  deallocate(Xd, stat=istat)
!!!!  call memocc(istat, iall, 'Xd', subname)
!!!!  iall=-product(shape(Yd))*kind(Yd)
!!!!  deallocate(Yd, stat=istat)
!!!!  call memocc(istat, iall, 'Yd', subname)
!!!!  iall=-product(shape(Zd))*kind(Zd)
!!!!  deallocate(Zd, stat=istat)
!!!!  call memocc(istat, iall, 'Zd', subname)
!!!!  iall=-product(shape(Xprime))*kind(Xprime)
!!!!  deallocate(Xprime, stat=istat)
!!!!  call memocc(istat, iall, 'Xprime', subname)
!!!!  iall=-product(shape(Yprime))*kind(Yprime)
!!!!  deallocate(Yprime, stat=istat)
!!!!  call memocc(istat, iall, 'Yprime', subname)
!!!!  iall=-product(shape(Zprime))*kind(Zprime)
!!!!  deallocate(Zprime, stat=istat)
!!!!  call memocc(istat, iall, 'Zprime', subname)
!!!!  iall=-product(shape(Xprimesquare))*kind(Xprimesquare)
!!!!  deallocate(Xprimesquare, stat=istat)
!!!!  call memocc(istat, iall, 'Xprimesquare', subname)
!!!!  iall=-product(shape(Yprimesquare))*kind(Yprimesquare)
!!!!  deallocate(Yprimesquare, stat=istat)
!!!!  call memocc(istat, iall, 'Yprimesquare', subname)
!!!!  iall=-product(shape(Zprimesquare))*kind(Zprimesquare)
!!!!  deallocate(Zprimesquare, stat=istat)
!!!!  call memocc(istat, iall, 'Zprimesquare', subname)
!!!!  iall=-product(shape(commutX))*kind(commutX)
!!!!  deallocate(commutX, stat=istat)
!!!!  call memocc(istat, iall, 'commutX', subname)
!!!!  iall=-product(shape(commutY))*kind(commutY)
!!!!  deallocate(commutY, stat=istat)
!!!!  call memocc(istat, iall, 'commutY', subname)
!!!!  iall=-product(shape(commutZ))*kind(commutZ)
!!!!  deallocate(commutZ, stat=istat)
!!!!  call memocc(istat, iall, 'commutZ', subname)
!!!!
!!!!
!!!!
!!!!  call deallocateRecvBufferOrtho(comon, subname)
!!!!  call deallocateSendBufferOrtho(comon, subname)
!!!!
!!!!end subroutine MLWF
