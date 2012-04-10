subroutine MLWFnew(iproc, nproc, lzd, orbs, at, op, comon, mad, rxyz, nit, kernel, &
           confdatarr, hx, locregCenters, maxDispl, lphi, Umat, centers)
use module_base
use module_types
use module_interfaces, exceptThisOne => MLWFnew
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, nit
type(local_zone_descriptors),intent(in):: lzd
type(orbitals_data),intent(in):: orbs
type(atoms_data),intent(in):: at
type(overlapParameters),intent(inout):: op
type(p2pComms),intent(inout):: comon
type(matrixDescriptors),intent(in):: mad
real(8),dimension(3,at%nat),intent(in):: rxyz
real(8),dimension(orbs%norb,orbs%norb),intent(in):: kernel
!logical,intent(in):: newgradient
real(8),intent(in):: hx, maxDispl
type(confpot_data),dimension(orbs%norbp),intent(in):: confdatarr
real(8),dimension(3,lzd%nlr),intent(in):: locregCenters
real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(inout):: lphi
real(8),dimension(orbs%norb,orbs%norb),intent(out):: Umat
real(8),dimension(3,lzd%nlr),intent(out):: centers

! Local variables
integer:: it, info, lwork, k, istat, iorb, jorb, iall, ierr, ist, jst, ilrold, ncount, jjorb, iiorb, ilr, lorb, jlr
integer:: nlocregOnMPI, jlrold, jj, ii, ncnt, order
real(8):: trace, lstep, dfactorial, energyconf_trial, energyconf_0, energyconf_der0, lstep_optimal, ddot
real(8):: tt1, tt2, tt3, tt4, tt5, tt, var
real(8):: t1, t2, t1_tot, t2_tot, omega
real(8):: time_convol, time_commun, time_lincomb, time_linalg, time_matrixmodification, time_exponential, time_tot
real(8):: time_matrixelements, rspread, locdiff
complex(8):: ttc
real(8),dimension(:,:),allocatable:: gmat, hamtrans, ttmat, Kmat, X, Y, Z, Xd, Yd, Zd, X2, Y2, Z2, R, R2
real(8),dimension(:,:),allocatable:: commutX1, commutY1, commutZ1, commutX2, commutY2, commutZ2, commutX3, commutY3, commutZ3
real(8),dimension(:,:),allocatable:: Xprime, Yprime, Zprime, Xprimesquare, Yprimesquare, Zprimesquare
real(8),dimension(:,:),allocatable:: X0, H, gHmat, Y0, M, gMmat, Z0, N, gNmat, centers_start, centers_end
real(8),dimension(:,:),allocatable:: HXd, HXdsquare, HXdX0, HX0, HX0square, Xdsquare, XdX0, X0square
real(8),dimension(:,:),allocatable:: MYd, MYdsquare, MYdY0, MY0, MY0square, Ydsquare, YdY0, Y0square
real(8),dimension(:,:),allocatable:: NZd, NZdsquare, NZdZ0, NZ0, NZ0square, Zdsquare, ZdZ0, Z0square
real(8),dimension(:,:,:),allocatable:: potmat, potmatsmall
complex(8),dimension(:,:),allocatable:: gmatc, omatc
complex(8),dimension(:,:,:),allocatable:: tempmatc
complex(8),dimension(:),allocatable:: work, expD_cmplx
real(8),dimension(:),allocatable:: rwork, eval, lphiovrlp, lvphi, lxphi, lyphi, lzphi, normarr
real(8),dimension(:,:,:),allocatable:: tempmat3
character(len=*),parameter:: subname='MLWFnew'
type(p2pComms):: comon_local

! Quick return if possible. In this way the localization regions will remain unchanged.
if(nit==-1) return

allocate(gmat(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, gmat, 'gmat', subname)
allocate(hamtrans(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, hamtrans, 'hamtrans', subname)
allocate(gmatc(orbs%norb,orbs%norb), stat=istat)
!call memocc(istat, gmatc, 'gmatc', subname)
allocate(omatc(orbs%norb,orbs%norb), stat=istat)
!call memocc(istat, omatc, 'omatc', subname)
allocate(tempmat3(orbs%norb,orbs%norb,3), stat=istat)
call memocc(istat, tempmat3, 'tempmat3', subname)
allocate(eval(orbs%norb), stat=istat)
call memocc(istat, eval, 'eval', subname)
allocate(expD_cmplx(orbs%norb), stat=istat)
call memocc(istat, expD_cmplx, 'expD_cmplx', subname)
allocate(tempmatc(orbs%norb,orbs%norb,2), stat=istat)
!call memocc(istat, tempmatc, 'tempmatc', subname)
allocate(lphiovrlp(op%ndim_lphiovrlp), stat=istat)
call memocc(istat, lphiovrlp, 'lphiovrlp',subname)
allocate(lvphi(max(orbs%npsidim_orbs,orbs%npsidim_comp)), stat=istat)
call memocc(istat, lvphi, 'lvphi', subname)
allocate(Kmat(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, Kmat, 'Kmat', subname)
!!allocate(recvbuf(comon%nrecvbuf), stat=istat)
!!call memocc(istat, recvbuf, 'recvbuf', subname)

allocate(potmat(orbs%norb,orbs%norb,at%nat), stat=istat)
call memocc(istat, potmat, 'potmat', subname)

allocate(lxphi(max(orbs%npsidim_orbs,orbs%npsidim_comp)), stat=istat)
call memocc(istat, lxphi, 'lxphi', subname)
allocate(lyphi(max(orbs%npsidim_orbs,orbs%npsidim_comp)), stat=istat)
call memocc(istat, lyphi, 'lyphi', subname)
allocate(lzphi(max(orbs%npsidim_orbs,orbs%npsidim_comp)), stat=istat)
call memocc(istat, lzphi, 'lzphi', subname)

allocate(X(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, X, 'X', subname)
allocate(Y(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, Y, 'Y', subname)
allocate(Z(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, Z, 'Z', subname)
allocate(X2(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, X2, 'X2', subname)
allocate(Y2(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, Y2, 'Y2', subname)
allocate(Z2(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, Z2, 'Z2', subname)
allocate(R(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, R, 'R', subname)
allocate(R2(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, R2, 'R2', subname)
allocate(Xd(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, Xd, 'Xd', subname)
allocate(Yd(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, Yd, 'Yd', subname)
allocate(Zd(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, Zd, 'Zd', subname)
allocate(Xprime(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, Xprime, 'Xprime', subname)
allocate(Yprime(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, Yprime, 'Yprime', subname)
allocate(Zprime(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, Zprime, 'Zprime', subname)
allocate(Xprimesquare(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, Xprimesquare, 'Xprimesquare', subname)
allocate(Yprimesquare(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, Yprimesquare, 'Yprimesquare', subname)
allocate(Zprimesquare(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, Zprimesquare, 'Zprimesquare', subname)

allocate(commutX1(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, commutX1, 'commutX1', subname)
allocate(commutY1(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, commutY1, 'commutY1', subname)
allocate(commutZ1(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, commutZ1, 'commutZ1', subname)
allocate(commutX2(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, commutX2, 'commutX2', subname)
allocate(commutY2(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, commutY2, 'commutY2', subname)
allocate(commutZ2(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, commutZ2, 'commutZ2', subname)
allocate(commutX3(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, commutX3, 'commutX3', subname)
allocate(commutY3(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, commutY3, 'commutY3', subname)
allocate(commutZ3(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, commutZ3, 'commutZ3', subname)

!!allocate(X0(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, X0, 'X0', subname)
!!allocate(H(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, H, 'H', subname)
!!allocate(gHmat(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, gHmat, 'gHmat', subname)
!!
!!allocate(Y0(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, Y0, 'Y0', subname)
!!allocate(M(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, M, 'M', subname)
!!allocate(gMmat(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, gMmat, 'gMmat', subname)
!!
!!allocate(Z0(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, Z0, 'Z0', subname)
!!allocate(N(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, N, 'N', subname)
!!allocate(gNmat(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, gNmat, 'gNmat', subname)
!!
!!allocate(HXd(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, HXd, 'HXd', subname)
!!allocate(HXdsquare(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, HXdsquare, 'HXdsquare', subname)
!!allocate(HXdX0(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, HXdX0, 'HXdX0', subname)
!!allocate(HX0(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, HX0, 'HX0', subname)
!!allocate(HX0square(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, HX0square, 'HX0square', subname)
!!allocate(Xdsquare(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, Xdsquare, 'Xdsquare', subname)
!!allocate(XdX0(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, XdX0, 'XdX0', subname)
!!allocate(X0square(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, X0square, 'X0square', subname)
!!
!!allocate(MYd(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, MYd, 'MYd', subname)
!!allocate(MYdsquare(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, MYdsquare, 'MYdsquare', subname)
!!allocate(MYdY0(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, MYdY0, 'MYdY0', subname)
!!allocate(MY0(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, MY0, 'MY0', subname)
!!allocate(MY0square(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, MY0square, 'MY0square', subname)
!!allocate(Ydsquare(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, Ydsquare, 'Ydsquare', subname)
!!allocate(YdY0(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, YdY0, 'YdY0', subname)
!!allocate(Y0square(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, Y0square, 'Y0square', subname)
!!
!!allocate(NZd(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, NZd, 'NZd', subname)
!!allocate(NZdsquare(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, NZdsquare, 'NZdsquare', subname)
!!allocate(NZdZ0(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, NZdZ0, 'NZdZ0', subname)
!!allocate(NZ0(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, NZ0, 'NZ0', subname)
!!allocate(NZ0square(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, NZ0square, 'NZ0square', subname)
!!allocate(Zdsquare(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, Zdsquare, 'Zdsquare', subname)
!!allocate(ZdZ0(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, ZdZ0, 'ZdZ0', subname)
!!allocate(Z0square(orbs%norb,orbs%norb), stat=istat)
!!call memocc(istat, Z0square, 'Z0square', subname)

allocate(centers_start(3,lzd%nlr), stat=istat)
call memocc(istat, centers_start, 'centers_start', subname)
allocate(centers_end(3,lzd%nlr), stat=istat)
call memocc(istat, centers_end, 'centers_end', subname)
allocate(normarr(orbs%norb), stat=istat)
call memocc(istat, normarr, 'normarr', subname)


! Count how many locregs each process handles
ilrold=-1
nlocregOnMPI=0
do iorb=1,orbs%norbp
    iiorb=orbs%isorb+iorb
    ilr=orbs%inwhichlocreg(iiorb)
    !if(ilr>ilrold) then
    if(ilr/=ilrold) then
        nlocregOnMPI=nlocregOnMPI+1
    end if
    ilrold=ilr
end do
allocate(potmatsmall(orbs%norb,orbs%norb,nlocregOnMPI), stat=istat)
call memocc(istat, potmatsmall, 'potmatsmall', subname)



  call allocateSendBufferOrtho(comon, subname)
  call allocateRecvBufferOrtho(comon, subname)
  ! Extract the overlap region from the orbitals phi and store them in comon%sendBuf.
  !if(nit>0) then
  !    call extractOrbital3(iproc, nproc, orbs, orbs%npsidim, orbs%inWhichLocreg, lzd, op, lphi, comon%nsendBuf, comon%sendBuf)
  !    call postCommsOverlapNew(iproc, nproc, orbs, op, lzd, lphi, comon, tt1, tt2)
  !end if

  energyconf_trial=0.d0 !just to initialize this variable and make the compiler happy
  lstep=0.d0 !just to initialize this variable and make the compiler happy
  lstep_optimal=0.d0 !just to initialize this variable and make the compiler happy
  energyconf_der0=0.d0 !just to initialize this variable and make the compiler happy

  time_convol=0.d0
  time_lincomb=0.d0
  time_commun=0.d0
  time_linalg=0.d0
  time_exponential=0.d0
  time_matrixmodification=0.d0
  time_matrixElements=0.d0
  t1_tot=mpi_wtime()

  ! Initialize Umat
  do iorb=1,orbs%norb
      do jorb=1,orbs%norb
          if(jorb==iorb) then
              Umat(jorb,iorb)=1.d0
          else
              Umat(jorb,iorb)=0.d0
          end if
      end do
  end do

  !!! Initialize Lagragme multipliers to zero
  !!H=0.d0
  !!M=0.d0
  !!N=0.d0

  !!! Initialize X0 etc.
  !!X0=0.d0
  !!Y0=0.d0
  !!Z0=0.d0
  !!do iorb=1,orbs%norbp
  !!    iiorb=orbs%isorb+iorb
  !!    X0(iiorb,iiorb)=confdatarr(iorb)%rxyzConf(1)
  !!    Y0(iiorb,iiorb)=confdatarr(iorb)%rxyzConf(2)
  !!    Z0(iiorb,iiorb)=confdatarr(iorb)%rxyzConf(3)
  !!end do
  !!call mpiallred(X0(1,1), orbs%norb**2, mpi_sum, mpi_comm_world, ierr)
  !!call mpiallred(Y0(1,1), orbs%norb**2, mpi_sum, mpi_comm_world, ierr)
  !!call mpiallred(Z0(1,1), orbs%norb**2, mpi_sum, mpi_comm_world, ierr)


  ! Calculate the norm of all basis functions
  normarr=0.d0
  ist=1
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      ilr=orbs%inwhichlocreg(iiorb)
      ncnt=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
      normarr(iiorb)=ddot(ncnt, lphi(ist), 1, lphi(ist), 1)
      ist=ist+ncnt
  end do
  call mpiallred(normarr(1), orbs%norb, mpi_sum, mpi_comm_world, ierr)
  !!if(iproc==0) then
  !!    do iorb=1,orbs%norb
  !!        write(*,'(a,i6,es16.6)') 'new norm before loop: iorb, normarr(iorb)', iorb, normarr(iorb)
  !!    end do
  !!end if


  !!! OTHER CHECK
  !!call apply_rminusmu_operator(iproc, nproc, orbs, lzd, hx, hx, hx, confdatarr, lphi, locregCenters, lvphi)
  !!call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lvphi, mad, R2)
  !!! New calculation of the spread
  !!var=0.d0
  !!do iorb=1,orbs%norb
  !!    tt = R2(iorb,iorb)/normarr(iorb)
  !!    var=var+tt
  !!end do
  !!!!if(iproc==0) write(*,'(a,es18.8)') 'FIRST: total variance', var

  innerLoop: do it=1,nit

  !write(*,*) '1: iproc, associated(comon%recvbuf)', iproc, associated(comon%recvbuf)

      t1=mpi_wtime()
      !!call apply_orbitaldependent_potential(iproc, nproc, at, orbs, lzd, rxyz, &
      !!     confdatarr, hx, lphi, -1, lvphi)

      if(it==1) then
          order=1
          call apply_position_operators(iproc, nproc, orbs, lzd, hx, hx, hx, confdatarr, lphi, order, lxphi, lyphi, lzphi)
          ! Build the matrices X, Y, Z
          call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lxphi, mad, X)
          call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lyphi, mad, Y)
          call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lzphi, mad, Z)

          order=2
          call apply_position_operators(iproc, nproc, orbs, lzd, hx, hx, hx, confdatarr, lphi, order, lxphi, lyphi, lzphi)
          ! Build the matrices X, Y, Z
          call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lxphi, mad, X2)
          call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lyphi, mad, Y2)
          call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lzphi, mad, Z2)

          order=1
          call apply_r_operators(iproc, nproc, orbs, lzd, hx, hx, hx, confdatarr, lphi, order, lvphi)
          call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lvphi, mad, R)

          order=2
          call apply_r_operators(iproc, nproc, orbs, lzd, hx, hx, hx, confdatarr, lphi, order, lvphi)
          call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lvphi, mad, R2)


          do iorb=1,orbs%norb
              ilr=orbs%inwhichlocreg(iorb)
              !!if(iproc==0) write(*,'(a,2i5,3f10.4,4x,3f10.4)') 'START: iorb, ilr, centers: ', &
              !!  iorb, ilr, X(iorb,iorb)/normarr(iorb), Y(iorb,iorb)/normarr(iorb), Z(iorb,iorb)/normarr(iorb), &
              !!  locregCenters(1,ilr),  locregCenters(2,ilr), locregCenters(3,ilr)
              centers_start(1,iorb)=X(iorb,iorb)/normarr(iorb)
              centers_start(2,iorb)=Y(iorb,iorb)/normarr(iorb)
              centers_start(3,iorb)=Z(iorb,iorb)/normarr(iorb)
          end do
      else
        call dgemm('t', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,1), orbs%norb, &
             X(1,1), orbs%norb, 0.d0, tempmat3(1,1,2), orbs%norb)
        call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,2), orbs%norb, &
             tempmat3(1,1,1), orbs%norb, 0.d0, X(1,1), orbs%norb)

        call dgemm('t', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,1), orbs%norb, &
             Y(1,1), orbs%norb, 0.d0, tempmat3(1,1,2), orbs%norb)
        call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,2), orbs%norb, &
             tempmat3(1,1,1), orbs%norb, 0.d0, Y(1,1), orbs%norb)

        call dgemm('t', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,1), orbs%norb, &
             Z(1,1), orbs%norb, 0.d0, tempmat3(1,1,2), orbs%norb)
        call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,2), orbs%norb, &
             tempmat3(1,1,1), orbs%norb, 0.d0, Z(1,1), orbs%norb)


        call dgemm('t', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,1), orbs%norb, &
             X2(1,1), orbs%norb, 0.d0, tempmat3(1,1,2), orbs%norb)
        call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,2), orbs%norb, &
             tempmat3(1,1,1), orbs%norb, 0.d0, X2(1,1), orbs%norb)

        call dgemm('t', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,1), orbs%norb, &
             Y2(1,1), orbs%norb, 0.d0, tempmat3(1,1,2), orbs%norb)
        call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,2), orbs%norb, &
             tempmat3(1,1,1), orbs%norb, 0.d0, Y2(1,1), orbs%norb)

        call dgemm('t', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,1), orbs%norb, &
             Z2(1,1), orbs%norb, 0.d0, tempmat3(1,1,2), orbs%norb)
        call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,2), orbs%norb, &
             tempmat3(1,1,1), orbs%norb, 0.d0, Z2(1,1), orbs%norb)


        call dgemm('t', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,1), orbs%norb, &
             R(1,1), orbs%norb, 0.d0, tempmat3(1,1,2), orbs%norb)
        call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,2), orbs%norb, &
             tempmat3(1,1,1), orbs%norb, 0.d0, R(1,1), orbs%norb)

        call dgemm('t', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,1), orbs%norb, &
             R2(1,1), orbs%norb, 0.d0, tempmat3(1,1,2), orbs%norb)
        call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,2), orbs%norb, &
             tempmat3(1,1,1), orbs%norb, 0.d0, R2(1,1), orbs%norb)
      end if

      ! New calculation of the spread
      !tt=0.d0
      var=0.d0
      do iorb=1,orbs%norb
          !!if(iproc==0) write(*,'(a,i7,9es15.6)') 'iorb, diagonal part of X,Y,Z,X2,Y2,Z2,R,R2,normarr', &
          !!    iorb, X(iorb,iorb), Y(iorb,iorb), Z(iorb,iorb), X2(iorb,iorb), Y2(iorb,iorb), Z2(iorb,iorb), &
          !!    R(iorb,iorb), R2(iorb,iorb), normarr(iorb)
          !!tt = X2(iorb,iorb)/normarr(iorb) + Y2(iorb,iorb)/normarr(iorb) + Z2(iorb,iorb)/normarr(iorb) &
          !!    -X(iorb,iorb)**2/normarr(iorb) - Y(iorb,iorb)**2/normarr(iorb) - Z(iorb,iorb)**2/normarr(iorb)
          !tt = R2(iorb,iorb)/normarr(iorb)-R(iorb,iorb)**2/normarr(iorb)
          tt = R2(iorb,iorb)/normarr(iorb)-(R(iorb,iorb)/normarr(iorb))**2
          var=var+tt
          !if(iproc==0) write(*,'(a,i8,es15.6)') 'iorb, tt', iorb, tt
      end do
      !if(iproc==0) write(*,'(a,i6,es18.8)') 'it, total variance', it, var
      !if(iproc==0) write(*,'(a,es18.8)') 'NEW SPREAD: ',tt

      ! Check the deviation from the center of the original localization region
      !do ilr=1,lzd%nlr
      do iorb=1,orbs%norb
          ilr=orbs%inwhichlocreg(iorb)
          tt = (locregCenters(1,ilr)-X(iorb,iorb)/normarr(iorb))**2 &
              +(locregCenters(2,ilr)-Y(iorb,iorb)/normarr(iorb))**2 &
              +(locregCenters(3,ilr)-Z(iorb,iorb)/normarr(iorb))**2 
          tt=sqrt(tt)
          if(tt>maxDispl) then
              if(iproc==0) write(*,'(a,i0,2x,es11.2)') 'WARNING: too large displacement for locreg ',ilr,tt
          end if
      end do

      !!h=0.d0
      !!m=0.d0
      !!n=0.d0
      !!y=0.d0
      !!z=0.d0
      !!y0=0.d0
      !!z0=0.d0

      !!tt=0.d0
      !!do iorb=1,orbs%norb
      !!    tt = tt + (X(iorb,iorb)-X0(iorb,iorb))**2 + (Y(iorb,iorb)-Y0(iorb,iorb))**2 + (Z(iorb,iorb)-Z0(iorb,iorb))**2
      !!end do
      !!if(iproc==0) write(*,'(a,es16.7)') 'difference to locreg center', tt

      ! Build the matrices Xd, Yd, Zd
      do iorb=1,orbs%norb
          do jorb=1,orbs%norb
              if(jorb==iorb) then
                  Xd(jorb,iorb)=X(jorb,iorb)
                  Yd(jorb,iorb)=Y(jorb,iorb)
                  Zd(jorb,iorb)=Z(jorb,iorb)
              else
                  Xd(jorb,iorb)=0.d0
                  Yd(jorb,iorb)=0.d0
                  Zd(jorb,iorb)=0.d0
              end if
          end do
      end do

      ! Build the matrices Xprime, Yprime, Zprime
      Xprime=X-Xd
      Yprime=Y-Yd
      Zprime=Z-Zd

      ! Calculate value of Omega
      call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, Xprime(1,1), orbs%norb, &
           Xprime(1,1), orbs%norb, 0.d0, Xprimesquare, orbs%norb)
      call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, Yprime(1,1), orbs%norb, &
           Yprime(1,1), orbs%norb, 0.d0, Yprimesquare, orbs%norb)
      call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, Zprime(1,1), orbs%norb, &
           Zprime(1,1), orbs%norb, 0.d0, Zprimesquare, orbs%norb)

      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, H(1,1), orbs%norb, &
      !!     Xd(1,1), orbs%norb, 0.d0, HXd, orbs%norb)
      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, M(1,1), orbs%norb, &
      !!     Yd(1,1), orbs%norb, 0.d0, MYd, orbs%norb)
      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, N(1,1), orbs%norb, &
      !!     Zd(1,1), orbs%norb, 0.d0, NZd, orbs%norb)

      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, HXd(1,1), orbs%norb, &
      !!     Xd(1,1), orbs%norb, 0.d0, HXdsquare, orbs%norb)
      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, MYd(1,1), orbs%norb, &
      !!     Yd(1,1), orbs%norb, 0.d0, MYdsquare, orbs%norb)
      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, NZd(1,1), orbs%norb, &
      !!     Zd(1,1), orbs%norb, 0.d0, NZdsquare, orbs%norb)

      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, HXd(1,1), orbs%norb, &
      !!     X0(1,1), orbs%norb, 0.d0, HXdX0, orbs%norb)
      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, MYd(1,1), orbs%norb, &
      !!     Y0(1,1), orbs%norb, 0.d0, MYdY0, orbs%norb)
      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, NZd(1,1), orbs%norb, &
      !!     Z0(1,1), orbs%norb, 0.d0, NZdZ0, orbs%norb)

      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, H(1,1), orbs%norb, &
      !!     X0(1,1), orbs%norb, 0.d0, HX0, orbs%norb)
      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, M(1,1), orbs%norb, &
      !!     Y0(1,1), orbs%norb, 0.d0, MY0, orbs%norb)
      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, N(1,1), orbs%norb, &
      !!     Z0(1,1), orbs%norb, 0.d0, NZ0, orbs%norb)

      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, HX0(1,1), orbs%norb, &
      !!     X0(1,1), orbs%norb, 0.d0, HX0square, orbs%norb)
      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, MY0(1,1), orbs%norb, &
      !!     Y0(1,1), orbs%norb, 0.d0, MY0square, orbs%norb)
      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, NZ0(1,1), orbs%norb, &
      !!     Z0(1,1), orbs%norb, 0.d0, NZ0square, orbs%norb)
      omega=0.d0
      rspread=0.d0
      do iorb=1,orbs%norb
          !!omega = omega + Xprimesquare(iorb,iorb)+Yprimesquare(iorb,iorb)+Zprimesquare(iorb,iorb) &
          !!              - (HXdsquare(iorb,iorb)-2.d0*HXdX0(iorb,iorb)+HX0square(iorb,iorb)) &
          !!              - (MYdsquare(iorb,iorb)-2.d0*MYdY0(iorb,iorb)+MY0square(iorb,iorb)) &
          !!              - (NZdsquare(iorb,iorb)-2.d0*NZdZ0(iorb,iorb)+NZ0square(iorb,iorb))
          omega = omega + Xprimesquare(iorb,iorb)+Yprimesquare(iorb,iorb)+Zprimesquare(iorb,iorb) 
          rspread = rspread + Xprimesquare(iorb,iorb)+Yprimesquare(iorb,iorb)+Zprimesquare(iorb,iorb)
      end do
      !if(iproc==0) write(*,'(a,i7,2es16.7)') 'it, omega, lstep', it, omega, lstep

      call commutator(orbs%norb, Xprime, Xd, commutX1)
      call commutator(orbs%norb, Yprime, Yd, commutY1)
      call commutator(orbs%norb, Zprime, Zd, commutZ1)

      !!call commutator(orbs%norb, HXd, X, commutX2)
      !!call commutator(orbs%norb, MYd, Y, commutY2)
      !!call commutator(orbs%norb, NZd, Z, commutZ2)

      !!call commutator(orbs%norb, HX0, X, commutX3)
      !!call commutator(orbs%norb, MY0, Y, commutY3)
      !!call commutator(orbs%norb, NZ0, Z, commutZ3)

      !!gmat=2.d0*( (commutX1+commutY1+commutZ1) &
      !!             -commutX2+commutX3 &
      !!             -commutY2+commutY3 &
      !!             -commutZ2+commutZ3 )
      gmat=2.d0*( commutX1+commutY1+commutZ1 )


      !!! Calculate gradient for Lagrange multipliers
      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, Xd(1,1), orbs%norb, &
      !!     Xd(1,1), orbs%norb, 0.d0, Xdsquare(1,1), orbs%norb)
      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, Xd(1,1), orbs%norb, &
      !!     X0(1,1), orbs%norb, 0.d0, XdX0(1,1), orbs%norb)
      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, X0(1,1), orbs%norb, &
      !!     X0(1,1), orbs%norb, 0.d0, X0square(1,1), orbs%norb)
      !!gHmat=Xdsquare-2.d0*XdX0+X0square
      !!H=H-1.d-3*gHmat

      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, Yd(1,1), orbs%norb, &
      !!     Yd(1,1), orbs%norb, 0.d0, Ydsquare(1,1), orbs%norb)
      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, Yd(1,1), orbs%norb, &
      !!     Y0(1,1), orbs%norb, 0.d0, YdY0(1,1), orbs%norb)
      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, Y0(1,1), orbs%norb, &
      !!     Y0(1,1), orbs%norb, 0.d0, Y0square(1,1), orbs%norb)
      !!gMmat=Ydsquare-2.d0*YdY0+Y0square
      !!M=M-1.d-3*gMmat

      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, Zd(1,1), orbs%norb, &
      !!     Zd(1,1), orbs%norb, 0.d0, Zdsquare(1,1), orbs%norb)
      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, Zd(1,1), orbs%norb, &
      !!     Z0(1,1), orbs%norb, 0.d0, ZdZ0(1,1), orbs%norb)
      !!call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, Z0(1,1), orbs%norb, &
      !!     Z0(1,1), orbs%norb, 0.d0, Z0square(1,1), orbs%norb)
      !!gNmat=Zdsquare-2.d0*ZdZ0+Z0square
      !!N=N-1.d-3*gNmat

      locdiff=0.d0
      !!do iorb=1,orbs%norb
      !!    locdiff = locdiff + gHmat(iorb,iorb) + gMmat(iorb,iorb) + gNmat(iorb,iorb)
      !!end do
      if(iproc==0 .and. verbose >2) write(*,'(a,i8,3es16.7,es12.4)')&
           'it, rspread, locdiff, omega, lstep', it, rspread, locdiff, omega, lstep



      !do iorb=1,orbs%norb
      !    do jorb=1,orbs%norb
      !        if(iproc==0) write(66,*) iorb,jorb,Kmat(jorb,iorb)
      !    end do
      !end do
      !call allocateSendBufferOrtho(comon, subname)
      !call allocateRecvBufferOrtho(comon, subname)
      !write(*,*) '3: iproc, associated(comon%recvbuf)', iproc, associated(comon%recvbuf)
      t2=mpi_wtime()
      time_matrixElements=time_matrixElements+t2-t1


      !!if(newgradient) then
      !!    call get_potential_matrices(iproc, nproc, at, orbs, lzd, op, comon, mad, rxyz, &
      !!         confdatarr, hx, lphi, potmat)
      !!    !call get_potential_matrices_new(iproc, nproc, lin, at, input, orbs, lzd, op, comon, rxyz, lphi, &
      !!    !     nlocregOnMPI, potmatsmall)
      !!end if


      !!if(.not.newgradient) then
          !energyconf_0=ddot(orbs%npsidim, lphi(1), 1, lvphi(1), 1)
          !call mpiallred(energyconf_0, 1, mpi_sum, mpi_comm_world, ierr)
          !!$$energyconf_0=0.d0
          !!$$do iorb=1,orbs%norb
          !!$$    energyconf_0 = energyconf_0 + Kmat(iorb,iorb)
          !!$$end do
      !!else
      !!    energyconf_0=0.d0
      !!    do iorb=1,orbs%norb
      !!        do jorb=1,orbs%norb
      !!            energyconf_0 = energyconf_0 + kernel(jorb,iorb)*Kmat(jorb,iorb)
      !!            !energyconf_0 = energyconf_0 + kernel(jorb,iorb)*Kmat(iorb,jorb)
      !!        end do
      !!    end do
      !!end if
      !!$$if(iproc==0) write(*,'(a,i6,3es20.10,2es17.7)') &
      !!$$             'it, energyconf_0, energyvonf_trial, energyconf_der0, lstep, lstep_optimal', &
      !!$$             it, energyconf_0, energyconf_trial, energyconf_der0, lstep, lstep_optimal

      t1=mpi_wtime()
      !!if(.not.newgradient) then
          ! Construct antisymmtric matrix Gmat
          !!$$do iorb=1,orbs%norb
          !!$$    do jorb=1,orbs%norb
          !!$$        gmat(jorb,iorb)=2.d0*(Kmat(jorb,iorb)-Kmat(iorb,jorb))
          !!$$    end do
          !!$$end do 
      !!else
      !!    !!! THIS IS THE OLD VERSION #############################################################################
      !!    do iorb=1,orbs%norb
      !!        ilr=orbs%inwhichlocreg(iorb)
      !!        do jorb=1,orbs%norb
      !!            jlr=orbs%inwhichlocreg(jorb)
      !!            tt=0.d0
      !!            do lorb=1,orbs%norb
      !!                tt = tt + kernel(jorb,lorb)*Kmat(lorb,iorb) - kernel(iorb,lorb)*Kmat(lorb,jorb) + &
      !!                          kernel(jorb,lorb)*potmat(lorb,iorb,jlr) - kernel(iorb,lorb)*potmat(lorb,jorb,ilr)
      !!            end do
      !!            gmat(jorb,iorb)=-tt
      !!            !if(iproc==0) then
      !!            !    write(77,*) iorb, jorb, gmat(jorb,iorb)
      !!            !end if
      !!        end do
      !!    end do 
      !!    ! ########################################################################################################
      !!    !!! THIS IS THE NEW VERSION
      !!    !!gmat=0.d0
      !!    !!ii=0
      !!    !!ilrold=-1
      !!    !!do iorb=1,orbs%norbp
      !!    !!    iiorb=orbs%isorb+iorb
      !!    !!    ilr=orbs%inwhichlocreg(iiorb)
      !!    !!    if(ilr>ilrold) then
      !!    !!        ii=ii+1
      !!    !!    end if
      !!    !!    do jorb=1,orbs%norb
      !!    !!        jlr=orbs%inwhichlocreg(jorb)
      !!    !!        tt=0.d0
      !!    !!        do lorb=1,orbs%norb
      !!    !!            !tt = tt + kernel(jorb,lorb)*Kmat(lorb,iiorb) - kernel(iiorb,lorb)*Kmat(lorb,jorb) + &
      !!    !!            !          - kernel(iiorb,lorb)*potmat(lorb,jorb,ilr)
      !!    !!            tt = tt + kernel(jorb,lorb)*Kmat(lorb,iiorb) - kernel(iiorb,lorb)*Kmat(lorb,jorb) + &
      !!    !!                      - kernel(iiorb,lorb)*potmatsmall(lorb,jorb,ii)
      !!    !!        end do
      !!    !!        gmat(jorb,iiorb)=-tt
      !!    !!    end do
      !!    !!    ilrold=ilr
      !!    !!end do 
      !!    !!do iorb=1,orbs%norb
      !!    !!    ilr=orbs%inwhichlocreg(iorb)
      !!    !!    jlrold=-1
      !!    !!    jj=0
      !!    !!    do jorb=1,orbs%norbp
      !!    !!        jjorb=orbs%isorb+jorb
      !!    !!        jlr=orbs%inwhichlocreg(jjorb)
      !!    !!        if(jlr>jlrold) then
      !!    !!            jj=jj+1
      !!    !!        end if
      !!    !!        tt=0.d0
      !!    !!        do lorb=1,orbs%norb
      !!    !!            !tt = tt + kernel(jjorb,lorb)*potmat(lorb,iorb,jlr)
      !!    !!            tt = tt + kernel(jjorb,lorb)*potmatsmall(lorb,iorb,jj)
      !!    !!        end do
      !!    !!        gmat(jjorb,iorb)=gmat(jjorb,iorb)-tt
      !!    !!        jlrold=jlr
      !!    !!    end do
      !!    !!end do 
      !!    !!call mpiallred(gmat(1,1), orbs%norb**2, mpi_sum, mpi_comm_world, ierr)
      !!    !!do iorb=1,orbs%norb
      !!    !!    do jorb=1,orbs%norb
      !!    !!        if(iproc==0) then
      !!    !!            write(77,*) iorb, jorb, gmat(jorb,iorb)
      !!    !!        end if
      !!    !!    end do
      !!    !!end do
      !!end if
      t2=mpi_wtime()
      time_matrixmodification=time_matrixmodification+t2-t1

      t1=mpi_wtime()
      !Build the complex matrix -iGmat
      do iorb=1,orbs%norb
          do jorb=1,orbs%norb
              gmatc(jorb,iorb)=cmplx(0.d0,-gmat(jorb,iorb),kind=8)
              !!if(iproc==0) write(999,'(a,2i8,2es16.7)') 'iorb, jorb, gmatc(jorb,iorb)', iorb, jorb, gmatc(jorb,iorb)
          end do
      end do 
      t2=mpi_wtime()
      time_matrixmodification=time_matrixmodification+t2-t1



      ! Diagonalize Gmatc
      t1=mpi_wtime()
      lwork=10*orbs%norb
      allocate(work(lwork), stat=istat) ! factor of 2 since it is assumed to be complex
      allocate(rwork(lwork), stat=istat)
      call zheev('v', 'l', orbs%norb, gmatc(1,1), orbs%norb, eval(1), work, lwork, rwork, info)
      if(info/=0) stop 'ERROR in zheev 3'
      deallocate(work)
      deallocate(rwork)
      t2=mpi_wtime()
      time_linalg=time_linalg+t2-t1


      !!$$! Calculate step size
      !!$$if(it==1) then
      !!$$    !!if(.not.newgradient) then
      !!$$        lstep=5.d-2/(maxval(eval))
      !!$$    !!else
      !!$$    !!    lstep=1.d-4/(maxval(eval))
      !!$$    !!end if
      !!$$else
      !!$$    lstep=2.d0*lstep_optimal
      !!$$    !lstep=1.d-3/(maxval(eval))
      !!$$end if
      lstep=-2.d-3/(maxval(eval))

      t1=mpi_wtime()
      ! Calculate exp(-i*l*D) (with D diagonal matrix of eigenvalues).
      ! This is also a diagonal matrix, so only calculate the diagonal part.
      do iorb=1,orbs%norb
         ttc=cmplx(0.d0,-lstep*eval(iorb),kind=8)
         expD_cmplx(iorb)=(0.d0,0.d0)
          do k=0,50
              expD_cmplx(iorb)=expD_cmplx(iorb)+ttc**k/dfactorial(k)
          end do
      end do
      t2=mpi_wtime()
      time_exponential=time_exponential+t2-t1

      t1=mpi_wtime()
      do iorb=1,orbs%norb
          do jorb=1,orbs%norb
              if(iorb==jorb) then
                  tempmatc(jorb,iorb,1)=expD_cmplx(iorb)
              else
                  tempmatc(jorb,iorb,1)=cmplx(0.d0,0.d0,kind=8)
              end if
          end do
      end do
      t2=mpi_wtime()
      time_matrixmodification=time_matrixmodification+t2-t1

      t1=mpi_wtime()
      call zgemm('n', 'c', orbs%norb, orbs%norb, orbs%norb, (1.d0,0.d0), tempmatc(1,1,1), orbs%norb, &
           gmatc(1,1), orbs%norb, (0.d0,0.d0), tempmatc(1,1,2), orbs%norb)
      call zgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, (1.d0,0.d0), gmatc(1,1), orbs%norb, &
           tempmatc(1,1,2), orbs%norb, (0.d0,0.d0), omatc(1,1), orbs%norb)
      t2=mpi_wtime()
      time_linalg=time_linalg+t2-t1

      t1=mpi_wtime()
      ! Build new lphi
      do iorb=1,orbs%norb
          do jorb=1,orbs%norb
              tempmat3(jorb,iorb,1)=real(omatc(jorb,iorb))
          end do
      end do
      t2=mpi_wtime()
      time_matrixmodification=time_matrixmodification+t2-t1

      ! Update Umat
      call dcopy(orbs%norb**2, Umat(1,1), 1, tempmat3(1,1,2), 1)
      call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,2), orbs%norb, &
           tempmat3(1,1,1), orbs%norb, 0.d0, Umat(1,1), orbs%norb)

      !!t1=mpi_wtime()
      !!!write(*,*) '5: iproc, associated(comon%recvbuf)', iproc, associated(comon%recvbuf)
      !!call build_new_linear_combinations(iproc, nproc, lzd, orbs, op, comon%nrecvbuf, &
      !!     comon%recvbuf, tempmat3(1,1,1), .true., lphi)
      !!t2=mpi_wtime()
      !!time_lincomb=time_lincomb+t2-t1

      !!$$t1=mpi_wtime()
      !!$$call apply_orbitaldependent_potential(iproc, nproc, at, orbs, lzd, rxyz, &
      !!$$     confdatarr, hx, lphi, -1, lvphi)
      !!$$t2=mpi_wtime()
      !!$$time_convol=time_convol+t2-t1

      !!$$!!if(.not.newgradient) then
      !!$$    energyconf_trial=ddot(max(orbs%npsidim_orbs,orbs%npsidim_comp), lphi(1), 1, lvphi(1), 1)
      !!$$    call mpiallred(energyconf_trial, 1, mpi_sum, mpi_comm_world, ierr)
      !!$$!!else

      !!$$!!    call dcopy(comon%nrecvbuf, comon%recvbuf, 1, recvbuf, 1)
      !!$$!!    !call extractOrbital3(iproc, nproc, orbs, orbs%npsidim, orbs%inWhichLocreg, lzd, op, lphi, comon%nsendBuf, comon%sendBuf)
      !!$$!!    !call postCommsOverlapNew(iproc, nproc, orbs, op, lzd, lphi, comon, tt1, tt2)
      !!$$!!    !allocate(ttmat(lin%orbs%norb,lin%orbs%norb))
      !!$$!!    !call collectnew(iproc, nproc, comon, lin%mad,lin%op, lin%orbs, input, lin%lzd, comon%nsendbuf, &
      !!$$!!    !     comon%sendbuf, comon%nrecvbuf, comon%recvbuf, ttmat, tt3, tt4, tt5)
      !!$$!!    !deallocate(ttmat)
      !!$$!!    call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lvphi, mad, Kmat)
      !!$$!!    call dcopy(comon%nrecvbuf, recvbuf, 1, comon%recvbuf, 1)

      !!$$!!    energyconf_trial=0.d0
      !!$$!!    do iorb=1,orbs%norb
      !!$$!!        do jorb=1,orbs%norb
      !!$$!!            energyconf_trial = energyconf_trial + kernel(jorb,iorb)*Kmat(jorb,iorb)
      !!$$!!            !energyconf_trial = energyconf_trial + kernel(jorb,iorb)*Kmat(iorb,jorb)
      !!$$!!        end do
      !!$$!!    end do
      !!$$!!end if

      !!$$! Calculate the gradient of the confinement
      !!$$energyconf_der0=0.d0
      !!$$do iorb=1,orbs%norb
      !!$$    do jorb=1,orbs%norb
      !!$$        energyconf_der0=energyconf_der0+gmat(jorb,iorb)**2
      !!$$    end do
      !!$$end do
      !!$$energyconf_der0=-.5d0*energyconf_der0

      !!$$! Calculate optimal step size
      !!$$lstep_optimal = -energyconf_der0*lstep**2/(2.d0*(energyconf_trial-energyconf_0-lstep*energyconf_der0))
      !!$$!!if(.not.newgradient) then
      !!$$    lstep_optimal=min(lstep_optimal,lstep)
      !!$$!!else
      !!$$!!    if(lstep_optimal<0) then
      !!$$!!        lstep_optimal=lstep
      !!$$!!    else
      !!$$!!        lstep_optimal=min(lstep_optimal,lstep)
      !!$$!!    end if
      !!$$!!end if

      !!$$t1=mpi_wtime()
      !!$$! Calculate exp(-i*l*D) (with D diagonal matrix of eigenvalues).
      !!$$! This is also a diagonal matrix, so only calculate the diagonal part.
      !!$$do iorb=1,orbs%norb
      !!$$   ttc=cmplx(0.d0,-lstep_optimal*eval(iorb),kind=8)
      !!$$   expD_cmplx(iorb)=(0.d0,0.d0)
      !!$$    do k=0,50
      !!$$        expD_cmplx(iorb)=expD_cmplx(iorb)+ttc**k/dfactorial(k)
      !!$$    end do
      !!$$end do
      !!$$t2=mpi_wtime()
      !!$$time_exponential=time_exponential+t2-t1

      !!$$t1=mpi_wtime()
      !!$$do iorb=1,orbs%norb
      !!$$    do jorb=1,orbs%norb
      !!$$        if(iorb==jorb) then
      !!$$            tempmatc(jorb,iorb,1)=expD_cmplx(iorb)
      !!$$        else
      !!$$            tempmatc(jorb,iorb,1)=cmplx(0.d0,0.d0,kind=8)
      !!$$        end if
      !!$$    end do
      !!$$end do
      !!$$t2=mpi_wtime()
      !!$$time_matrixmodification=time_matrixmodification+t2-t1

      !!$$t1=mpi_wtime()
      !!$$call zgemm('n', 'c', orbs%norb, orbs%norb, orbs%norb, (1.d0,0.d0), tempmatc(1,1,1), orbs%norb, &
      !!$$     gmatc(1,1), orbs%norb, (0.d0,0.d0), tempmatc(1,1,2), orbs%norb)
      !!$$call zgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, (1.d0,0.d0), gmatc(1,1), orbs%norb, &
      !!$$     tempmatc(1,1,2), orbs%norb, (0.d0,0.d0), omatc(1,1), orbs%norb)
      !!$$t2=mpi_wtime()
      !!$$time_linalg=time_linalg+t2-t1


      !!$$! Build new lphi
      !!$$do iorb=1,orbs%norb
      !!$$    do jorb=1,orbs%norb
      !!$$        tempmat3(jorb,iorb,1)=real(omatc(jorb,iorb))
      !!$$    end do
      !!$$end do
      !!$$t1=mpi_wtime()
      !!$$call build_new_linear_combinations(iproc, nproc, lzd, orbs, op, comon%nrecvbuf, &
      !!$$     comon%recvbuf, tempmat3(1,1,1), .true., lphi)
      !!$$t2=mpi_wtime()
      !!$$time_lincomb=time_lincomb+t2-t1


      !!$$!if(it<nit) then
      !!$$!    call extractOrbital3(iproc, nproc, orbs, orbs%npsidim, orbs%inWhichLocreg, lzd, op, lphi, comon%nsendBuf, comon%sendBuf)
      !!$$!    call postCommsOverlapNew(iproc, nproc, orbs, op, lzd, lphi, comon, tt1, tt2)
      !!$$!end if


  end do innerLoop


  !!!!!! EXPERIMENTAL
  !!!!write(*,*) 'WARNING: TRANSPOSE UMAT!!!!'
  !!!!call dcopy(orbs%norb**2, Umat(1,1), 1, tempmat3(1,1,1), 1)
  !!!!do iorb=1,orbs%norb
  !!!!    do jorb=1,orbs%norb
  !!!!        Umat(jorb,iorb)=tempmat3(iorb,jorb,1)
  !!!!    end do
  !!!!end do


  !!if(iproc==0) then
  !!    do iorb=1,orbs%norb
  !!        do jorb=1,orbs%norb
  !!            write(999,*) iorb,jorb,Umat(jorb,iorb)
  !!        end do
  !!    end do
  !!end if


  t1=mpi_wtime()
  !write(*,*) '5: iproc, associated(comon%recvbuf)', iproc, associated(comon%recvbuf)
  call extractOrbital3(iproc, nproc, orbs, orbs, orbs%npsidim_orbs, orbs%inWhichLocreg, lzd, lzd, &
       op, op, lphi, comon%nsendBuf, comon%sendBuf)
  call postCommsOverlapNew(iproc, nproc, orbs, op, lzd, lphi, comon, tt1, tt2)
  call collectnew(iproc, nproc, comon, mad, op, orbs, lzd, comon%nsendbuf, &
       comon%sendbuf, comon%nrecvbuf, comon%recvbuf, tt3, tt4, tt5)
  call build_new_linear_combinations(iproc, nproc, lzd, orbs, op, comon%nrecvbuf, &
       comon%recvbuf, Umat(1,1), .true., lphi)
  t2=mpi_wtime()
  time_lincomb=time_lincomb+t2-t1


  ! Recalculate the norm of all basis functions
  normarr=0.d0
  ist=1
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      ilr=orbs%inwhichlocreg(iiorb)
      ncnt=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
      normarr(iiorb)=ddot(ncnt, lphi(ist), 1, lphi(ist), 1)
      ist=ist+ncnt
  end do
  call mpiallred(normarr(1), orbs%norb, mpi_sum, mpi_comm_world, ierr)
  !!if(iproc==0) then
  !!    do iorb=1,orbs%norb
  !!        write(*,'(a,i6,es16.6)') 'new norm after loop: iorb, normarr(iorb)', iorb, normarr(iorb)
  !!    end do
  !!end if

  ! Recalculate X,Y,Z
  order=1
  call apply_position_operators(iproc, nproc, orbs, lzd, hx, hx, hx, confdatarr, lphi, order, lxphi, lyphi, lzphi)
  ! Build the matrices X, Y, Z
  call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lxphi, mad, X)
  call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lyphi, mad, Y)
  call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lzphi, mad, Z)


!!!  !!! CHECK #####################################
!!!  order=1
!!!  call apply_r_operators(iproc, nproc, orbs, lzd, hx, hx, hx, confdatarr, lphi, order, lvphi)
!!!  call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lvphi, mad, R)
!!!
!!!  order=2
!!!  call apply_r_operators(iproc, nproc, orbs, lzd, hx, hx, hx, confdatarr, lphi, order, lvphi)
!!!  call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lvphi, mad, R2)
!!!
!!!  call dgemm('t', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,1), orbs%norb, &
!!!       R(1,1), orbs%norb, 0.d0, tempmat3(1,1,2), orbs%norb)
!!!  call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,2), orbs%norb, &
!!!       tempmat3(1,1,1), orbs%norb, 0.d0, R(1,1), orbs%norb)
!!!
!!!  call dgemm('t', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,1), orbs%norb, &
!!!       R2(1,1), orbs%norb, 0.d0, tempmat3(1,1,2), orbs%norb)
!!!  call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, tempmat3(1,1,2), orbs%norb, &
!!!       tempmat3(1,1,1), orbs%norb, 0.d0, R2(1,1), orbs%norb)
!!!
!!!  ! New calculation of the spread
!!!  var=0.d0
!!!  do iorb=1,orbs%norb
!!!      tt = R2(iorb,iorb)/normarr(iorb)-(R(iorb,iorb)/normarr(iorb))**2
!!!      var=var+tt
!!!  end do
!!!  if(iproc==0) write(*,'(a,es18.8)') 'FINAL: total variance', var



  do iorb=1,orbs%norb
      ilr=orbs%inwhichlocreg(iorb)
      !!if(iproc==0) write(*,'(a,2i5,3f10.4)') 'END: iorb, ilr, centers: ', &
      !!    iorb, ilr, X(iorb,iorb)/normarr(iorb), Y(iorb,iorb)/normarr(iorb), Z(iorb,iorb)/normarr(iorb)
      !!centers(1,iorb)=X(iorb,iorb)/normarr(iorb)
      !!centers(2,iorb)=Y(iorb,iorb)/normarr(iorb)
      !!centers(3,iorb)=Z(iorb,iorb)/normarr(iorb)
      centers(1,ilr)=X(iorb,iorb)/normarr(iorb)
      centers(2,ilr)=Y(iorb,iorb)/normarr(iorb)
      centers(3,ilr)=Z(iorb,iorb)/normarr(iorb)
      centers_end(1,iorb)=X(iorb,iorb)/normarr(iorb)
      centers_end(2,iorb)=Y(iorb,iorb)/normarr(iorb)
      centers_end(3,iorb)=Z(iorb,iorb)/normarr(iorb)
  end do

  !!write(*,*) 'ATTENTION DEBUG'
  !!centers=locregCenters

  if(nit>0) then
      do iorb=1,lzd%nlr
          !!write(1000+iproc,'(6es12.4)') centers_start(1,iorb), centers_start(2,iorb), centers_start(3,iorb), &
          !!                    centers_end(1,iorb), centers_end(2,iorb), centers_end(3,iorb)
          tt = (centers_start(1,iorb)-centers_end(1,iorb))**2 &
              +(centers_start(2,iorb)-centers_end(2,iorb))**2 &
              +(centers_start(3,iorb)-centers_end(3,iorb))**2 
          tt=sqrt(tt)
          !!if(iproc==0) write(*,'(a,i5,es9.2)') 'iorb, shift of the center: ', iorb, tt
      end do
  end if


!!  ! OTHER CHECK
!!  call apply_rminusmu_operator(iproc, nproc, orbs, lzd, hx, hx, hx, confdatarr, lphi, centers, lvphi)
!!  call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lvphi, mad, R2)
!!  ! New calculation of the spread
!!  var=0.d0
!!  do iorb=1,orbs%norb
!!      tt = R2(iorb,iorb)/normarr(iorb)
!!      var=var+tt
!!  end do
!!  if(iproc==0) write(*,'(a,es18.8)') 'FINAL2: total variance', var
!!
!!  call apply_rminusmu_operator(iproc, nproc, orbs, lzd, hx, hx, hx, confdatarr, lphi, locregCenters, lvphi)
!!  call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lvphi, mad, R2)
!!  ! New calculation of the spread
!!  var=0.d0
!!  do iorb=1,orbs%norb
!!      tt = R2(iorb,iorb)/normarr(iorb)
!!      var=var+tt
!!  end do
!!  if(iproc==0) write(*,'(a,es18.8)') 'FINAL3: total variance', var



  t2_tot=mpi_wtime()
  time_tot=t2_tot-t1_tot
  call mpiallred(time_convol, 1, mpi_max, mpi_comm_world, ierr)
  call mpiallred(time_commun, 1, mpi_max, mpi_comm_world, ierr)
  call mpiallred(time_lincomb, 1, mpi_max, mpi_comm_world, ierr)
  call mpiallred(time_linalg, 1, mpi_max, mpi_comm_world, ierr)
  call mpiallred(time_matrixmodification, 1, mpi_max, mpi_comm_world, ierr)
  call mpiallred(time_exponential, 1, mpi_max, mpi_comm_world, ierr)
  call mpiallred(time_matrixelements, 1, mpi_max, mpi_comm_world, ierr)
  call mpiallred(time_tot, 1, mpi_max, mpi_comm_world, ierr)
  !time_convol=time_convol/dble(nproc)
  !time_commun=time_commun/dble(nproc)
  !time_lincomb=time_lincomb/dble(nproc)
  !time_linalg=time_linalg/dble(nproc)
  !time_matrixmodification=time_matrixmodification/dble(nproc)
  !time_exponential_=time_exponential/dble(nproc)
  !time_tot=time_tot/dble(nproc)
  !!if(iproc==0) then
  !!    write(*,'(a,es16.6)') 'total time: ',time_tot
  !!    write(*,'(a,es15.7,a,f5.2,a)') 'convolutions: ',time_convol,' (',time_convol/time_tot*100.d0,'%)'
  !!    write(*,'(a,es15.7,a,f5.2,a)') 'linear combinations: ',time_lincomb,' (',time_lincomb/time_tot*100.d0,'%)'
  !!    write(*,'(a,es15.7,a,f5.2,a)') 'communication: ',time_commun,' (',time_commun/time_tot*100.d0,'%)'
  !!    write(*,'(a,es15.7,a,f5.2,a)') 'linear algebra: ',time_linalg,' (',time_linalg/time_tot*100.d0,'%)'
  !!    write(*,'(a,es15.7,a,f5.2,a)') 'matrix modification: ',time_matrixmodification, &
  !!                                   ' (',time_matrixmodification/time_tot*100.d0,'%)'
  !!    write(*,'(a,es15.7,a,f5.2,a)') 'building exponential: ',time_exponential,' (',time_exponential/time_tot*100.d0,'%)'
  !!    write(*,'(a,es15.7,a,f5.2,a)') 'matrix elements ',time_matrixelements,' (',time_matrixelements/time_tot*100.d0,'%)'
  !!end if


  iall=-product(shape(gmat))*kind(gmat)
  deallocate(gmat, stat=istat)
  call memocc(istat, iall, 'gmat', subname)
  iall=-product(shape(gmatc))*kind(gmatc)
  deallocate(gmatc, stat=istat)
  !call memocc(istat, iall, 'gmatc', subname)
  iall=-product(shape(omatc))*kind(omatc)
  deallocate(omatc, stat=istat)
  !call memocc(istat, iall, 'omatc', subname)
  iall=-product(shape(tempmat3))*kind(tempmat3)
  deallocate(tempmat3, stat=istat)
  call memocc(istat, iall, 'tempmat3', subname)
  iall=-product(shape(eval))*kind(eval)
  deallocate(eval, stat=istat)
  call memocc(istat, iall, 'eval', subname)
  iall=-product(shape(expD_cmplx))*kind(expD_cmplx)
  deallocate(expD_cmplx, stat=istat)
  call memocc(istat, iall, 'expD_cmplx', subname)
  iall=-product(shape(tempmatc))*kind(tempmatc)
  deallocate(tempmatc, stat=istat)
  !call memocc(istat, iall, 'tempmatc', subname)
  iall=-product(shape(hamtrans))*kind(hamtrans)
  deallocate(hamtrans, stat=istat)
  call memocc(istat, iall, 'hamtrans', subname)
  iall=-product(shape(lphiovrlp))*kind(lphiovrlp)
  deallocate(lphiovrlp, stat=istat)
  call memocc(istat, iall, 'lphiovrlp', subname)
  iall=-product(shape(lvphi))*kind(lvphi)
  deallocate(lvphi, stat=istat)
  call memocc(istat, iall, 'lvphi', subname)
  iall=-product(shape(Kmat))*kind(Kmat)
  deallocate(Kmat, stat=istat)
  call memocc(istat, iall, 'Kmat', subname)
  !!iall=-product(shape(recvbuf))*kind(recvbuf)
  !!deallocate(recvbuf, stat=istat)
  !!call memocc(istat, iall, 'recvbuf', subname)

  iall=-product(shape(potmat))*kind(potmat)
  deallocate(potmat, stat=istat)
  call memocc(istat, iall, 'potmat', subname)
  iall=-product(shape(potmatsmall))*kind(potmatsmall)
  deallocate(potmatsmall, stat=istat)
  call memocc(istat, iall, 'potmatsmall', subname)


  iall=-product(shape(lxphi))*kind(lxphi)
  deallocate(lxphi, stat=istat)
  call memocc(istat, iall, 'lxphi', subname)
  iall=-product(shape(lyphi))*kind(lyphi)
  deallocate(lyphi, stat=istat)
  call memocc(istat, iall, 'lyphi', subname)
  iall=-product(shape(lzphi))*kind(lzphi)
  deallocate(lzphi, stat=istat)
  call memocc(istat, iall, 'lzphi', subname)



  iall=-product(shape(X))*kind(X)
  deallocate(X, stat=istat)
  call memocc(istat, iall, 'X', subname)
  iall=-product(shape(Y))*kind(Y)
  deallocate(Y, stat=istat)
  call memocc(istat, iall, 'Y', subname)
  iall=-product(shape(Z))*kind(Z)
  deallocate(Z, stat=istat)
  call memocc(istat, iall, 'Z', subname)
  iall=-product(shape(X2))*kind(X2)
  deallocate(X2, stat=istat)
  call memocc(istat, iall, 'X2', subname)
  iall=-product(shape(R))*kind(R)
  deallocate(R, stat=istat)
  call memocc(istat, iall, 'R', subname)
  iall=-product(shape(R2))*kind(R2)
  deallocate(R2, stat=istat)
  call memocc(istat, iall, 'R2', subname)
  iall=-product(shape(Y2))*kind(Y2)
  deallocate(Y2, stat=istat)
  call memocc(istat, iall, 'Y2', subname)
  iall=-product(shape(Z2))*kind(Z2)
  deallocate(Z2, stat=istat)
  call memocc(istat, iall, 'Z2', subname)
  iall=-product(shape(Xd))*kind(Xd)
  deallocate(Xd, stat=istat)
  call memocc(istat, iall, 'Xd', subname)
  iall=-product(shape(Yd))*kind(Yd)
  deallocate(Yd, stat=istat)
  call memocc(istat, iall, 'Yd', subname)
  iall=-product(shape(Zd))*kind(Zd)
  deallocate(Zd, stat=istat)
  call memocc(istat, iall, 'Zd', subname)
  iall=-product(shape(Xprime))*kind(Xprime)
  deallocate(Xprime, stat=istat)
  call memocc(istat, iall, 'Xprime', subname)
  iall=-product(shape(Yprime))*kind(Yprime)
  deallocate(Yprime, stat=istat)
  call memocc(istat, iall, 'Yprime', subname)
  iall=-product(shape(Zprime))*kind(Zprime)
  deallocate(Zprime, stat=istat)
  call memocc(istat, iall, 'Zprime', subname)
  iall=-product(shape(Xprimesquare))*kind(Xprimesquare)
  deallocate(Xprimesquare, stat=istat)
  call memocc(istat, iall, 'Xprimesquare', subname)
  iall=-product(shape(Yprimesquare))*kind(Yprimesquare)
  deallocate(Yprimesquare, stat=istat)
  call memocc(istat, iall, 'Yprimesquare', subname)
  iall=-product(shape(Zprimesquare))*kind(Zprimesquare)
  deallocate(Zprimesquare, stat=istat)
  call memocc(istat, iall, 'Zprimesquare', subname)

  iall=-product(shape(commutX1))*kind(commutX1)
  deallocate(commutX1, stat=istat)
  call memocc(istat, iall, 'commutX1', subname)
  iall=-product(shape(commutY1))*kind(commutY1)
  deallocate(commutY1, stat=istat)
  call memocc(istat, iall, 'commutY1', subname)
  iall=-product(shape(commutZ1))*kind(commutZ1)
  deallocate(commutZ1, stat=istat)
  call memocc(istat, iall, 'commutZ1', subname)
  iall=-product(shape(commutX2))*kind(commutX2)
  deallocate(commutX2, stat=istat)
  call memocc(istat, iall, 'commutX2', subname)
  iall=-product(shape(commutY2))*kind(commutY2)
  deallocate(commutY2, stat=istat)
  call memocc(istat, iall, 'commutY2', subname)
  iall=-product(shape(commutZ2))*kind(commutZ2)
  deallocate(commutZ2, stat=istat)
  call memocc(istat, iall, 'commutZ2', subname)
  iall=-product(shape(commutX3))*kind(commutX3)
  deallocate(commutX3, stat=istat)
  call memocc(istat, iall, 'commutX3', subname)
  iall=-product(shape(commutY3))*kind(commutY3)
  deallocate(commutY3, stat=istat)
  call memocc(istat, iall, 'commutY3', subname)
  iall=-product(shape(commutZ3))*kind(commutZ3)
  deallocate(commutZ3, stat=istat)
  call memocc(istat, iall, 'commutZ3', subname)

  !!iall=-product(shape(X0))*kind(X0)
  !!deallocate(X0, stat=istat)
  !!call memocc(istat, iall, 'X0', subname)
  !!iall=-product(shape(H))*kind(H)
  !!deallocate(H, stat=istat)
  !!call memocc(istat, iall, 'H', subname)
  !!iall=-product(shape(gHmat))*kind(gHmat)
  !!deallocate(gHmat, stat=istat)
  !!call memocc(istat, iall, 'gHmat', subname)

  !!iall=-product(shape(Y0))*kind(Y0)
  !!deallocate(Y0, stat=istat)
  !!call memocc(istat, iall, 'Y0', subname)
  !!iall=-product(shape(M))*kind(M)
  !!deallocate(M, stat=istat)
  !!call memocc(istat, iall, 'M', subname)
  !!iall=-product(shape(gMmat))*kind(gMmat)
  !!deallocate(gMmat, stat=istat)
  !!call memocc(istat, iall, 'gMmat', subname)

  !!iall=-product(shape(Z0))*kind(Z0)
  !!deallocate(Z0, stat=istat)
  !!call memocc(istat, iall, 'Z0', subname)
  !!iall=-product(shape(N))*kind(N)
  !!deallocate(N, stat=istat)
  !!call memocc(istat, iall, 'N', subname)
  !!iall=-product(shape(gNmat))*kind(gNmat)
  !!deallocate(gNmat, stat=istat)
  !!call memocc(istat, iall, 'gNmat', subname)


  !!iall=-product(shape(HXd))*kind(HXd)
  !!deallocate(HXd, stat=istat)
  !!call memocc(istat, iall, 'HXd', subname)
  !!iall=-product(shape(HXdsquare))*kind(HXdsquare)
  !!deallocate(HXdsquare, stat=istat)
  !!call memocc(istat, iall, 'HXdsquare', subname)
  !!iall=-product(shape(HXdX0))*kind(HXdX0)
  !!deallocate(HXdX0, stat=istat)
  !!call memocc(istat, iall, 'HXdX0', subname)
  !!iall=-product(shape(HX0))*kind(HX0)
  !!deallocate(HX0, stat=istat)
  !!call memocc(istat, iall, 'HX0', subname)
  !!iall=-product(shape(HX0square))*kind(HX0square)
  !!deallocate(HX0square, stat=istat)
  !!call memocc(istat, iall, 'HX0square', subname)
  !!iall=-product(shape(Xdsquare))*kind(Xdsquare)
  !!deallocate(Xdsquare, stat=istat)
  !!call memocc(istat, iall, 'Xdsquare', subname)
  !!iall=-product(shape(XdX0))*kind(XdX0)
  !!deallocate(XdX0, stat=istat)
  !!call memocc(istat, iall, 'XdX0', subname)
  !!iall=-product(shape(X0square))*kind(X0square)
  !!deallocate(X0square, stat=istat)
  !!call memocc(istat, iall, 'X0square', subname)

  !!iall=-product(shape(MYd))*kind(MYd)
  !!deallocate(MYd, stat=istat)
  !!call memocc(istat, iall, 'MYd', subname)
  !!iall=-product(shape(MYdsquare))*kind(MYdsquare)
  !!deallocate(MYdsquare, stat=istat)
  !!call memocc(istat, iall, 'MYdsquare', subname)
  !!iall=-product(shape(MYdY0))*kind(MYdY0)
  !!deallocate(MYdY0, stat=istat)
  !!call memocc(istat, iall, 'MYdY0', subname)
  !!iall=-product(shape(MY0))*kind(MY0)
  !!deallocate(MY0, stat=istat)
  !!call memocc(istat, iall, 'MY0', subname)
  !!iall=-product(shape(MY0square))*kind(MY0square)
  !!deallocate(MY0square, stat=istat)
  !!call memocc(istat, iall, 'MY0square', subname)
  !!iall=-product(shape(Ydsquare))*kind(Ydsquare)
  !!deallocate(Ydsquare, stat=istat)
  !!call memocc(istat, iall, 'Ydsquare', subname)
  !!iall=-product(shape(YdY0))*kind(YdY0)
  !!deallocate(YdY0, stat=istat)
  !!call memocc(istat, iall, 'YdY0', subname)
  !!iall=-product(shape(Y0square))*kind(Y0square)
  !!deallocate(Y0square, stat=istat)
  !!call memocc(istat, iall, 'Y0square', subname)

  !!iall=-product(shape(NZd))*kind(NZd)
  !!deallocate(NZd, stat=istat)
  !!call memocc(istat, iall, 'NZd', subname)
  !!iall=-product(shape(NZdsquare))*kind(NZdsquare)
  !!deallocate(NZdsquare, stat=istat)
  !!call memocc(istat, iall, 'NZdsquare', subname)
  !!iall=-product(shape(NZdZ0))*kind(NZdZ0)
  !!deallocate(NZdZ0, stat=istat)
  !!call memocc(istat, iall, 'NZdZ0', subname)
  !!iall=-product(shape(NZ0))*kind(NZ0)
  !!deallocate(NZ0, stat=istat)
  !!call memocc(istat, iall, 'NZ0', subname)
  !!iall=-product(shape(NZ0square))*kind(NZ0square)
  !!deallocate(NZ0square, stat=istat)
  !!call memocc(istat, iall, 'NZ0square', subname)
  !!iall=-product(shape(Zdsquare))*kind(Zdsquare)
  !!deallocate(Zdsquare, stat=istat)
  !!call memocc(istat, iall, 'Zdsquare', subname)
  !!iall=-product(shape(ZdZ0))*kind(ZdZ0)
  !!deallocate(ZdZ0, stat=istat)
  !!call memocc(istat, iall, 'ZdZ0', subname)
  !!iall=-product(shape(Z0square))*kind(Z0square)
  !!deallocate(Z0square, stat=istat)
  !!call memocc(istat, iall, 'Z0square', subname)

  iall=-product(shape(centers_start))*kind(centers_start)
  deallocate(centers_start, stat=istat)
  call memocc(istat, iall, 'centers_start', subname)
  iall=-product(shape(centers_end))*kind(centers_end)
  deallocate(centers_end, stat=istat)
  call memocc(istat, iall, 'centers_end', subname)
  iall=-product(shape(normarr))*kind(normarr)
  deallocate(normarr, stat=istat)
  call memocc(istat, iall, 'normarr', subname)


  call deallocateRecvBufferOrtho(comon, subname)
  call deallocateSendBufferOrtho(comon, subname)

end subroutine MLWFnew









!> Some description of the routine goes here
subroutine unitary_optimization(iproc, nproc, lzd, orbs, at, op, comon, mad, rxyz, nit, kernel, &
           newgradient, confdatarr, hx, lphi)
use module_base
use module_types
use module_interfaces, exceptThisOne => unitary_optimization
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, nit
type(local_zone_descriptors),intent(in):: lzd
type(orbitals_data),intent(in):: orbs
type(atoms_data),intent(in):: at
type(overlapParameters),intent(inout):: op
type(p2pComms),intent(inout):: comon
type(matrixDescriptors),intent(in):: mad
real(8),dimension(3,at%nat),intent(in):: rxyz
real(8),dimension(orbs%norb,orbs%norb),intent(in):: kernel
logical,intent(in):: newgradient
real(8),intent(in):: hx
type(confpot_data),dimension(orbs%norbp),intent(in):: confdatarr
real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(inout):: lphi

! Local variables
integer:: it, info, lwork, k, istat, iorb, jorb, iall, ierr, ist, jst, ilrold, ncount, jjorb, iiorb, ilr, lorb, jlr
integer:: nlocregOnMPI, jlrold, jj, ii
real(8):: trace, lstep, dfactorial, energyconf_trial, energyconf_0, energyconf_der0, lstep_optimal, ddot
real(8):: tt1, tt2, tt3, tt4, tt5, tt
real(8):: t1, t2, t1_tot, t2_tot
real(8):: time_convol, time_commun, time_lincomb, time_linalg, time_matrixmodification, time_exponential, time_tot
real(8):: time_matrixelements
complex(8):: ttc
real(8),dimension(:,:),allocatable:: gmat, hamtrans, ttmat, Kmat
real(8),dimension(:,:,:),allocatable:: potmat, potmatsmall
complex(8),dimension(:,:),allocatable:: gmatc, omatc
complex(8),dimension(:,:,:),allocatable:: tempmatc
complex(8),dimension(:),allocatable:: work, expD_cmplx
real(8),dimension(:),allocatable:: rwork, eval, lphiovrlp, lvphi, recvbuf
real(8),dimension(:,:,:),allocatable:: tempmat3
character(len=*),parameter:: subname='unitary_optimization'
type(p2pComms):: comon_local

! Quick return if possible
if(nit==0) return

allocate(gmat(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, gmat, 'gmat', subname)
allocate(hamtrans(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, hamtrans, 'hamtrans', subname)
allocate(gmatc(orbs%norb,orbs%norb), stat=istat)
!call memocc(istat, gmatc, 'gmatc', subname)
allocate(omatc(orbs%norb,orbs%norb), stat=istat)
!call memocc(istat, omatc, 'omatc', subname)
allocate(tempmat3(orbs%norb,orbs%norb,3), stat=istat)
call memocc(istat, tempmat3, 'tempmat3', subname)
allocate(eval(orbs%norb), stat=istat)
call memocc(istat, eval, 'eval', subname)
allocate(expD_cmplx(orbs%norb), stat=istat)
call memocc(istat, expD_cmplx, 'expD_cmplx', subname)
allocate(tempmatc(orbs%norb,orbs%norb,2), stat=istat)
!call memocc(istat, tempmatc, 'tempmatc', subname)
allocate(lphiovrlp(op%ndim_lphiovrlp), stat=istat)
call memocc(istat, lphiovrlp, 'lphiovrlp',subname)
allocate(lvphi(max(orbs%npsidim_orbs,orbs%npsidim_comp)), stat=istat)
call memocc(istat, lvphi, 'lvphi', subname)
allocate(Kmat(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, Kmat, 'Kmat', subname)
allocate(recvbuf(comon%nrecvbuf), stat=istat)
call memocc(istat, recvbuf, 'recvbuf', subname)

allocate(potmat(orbs%norb,orbs%norb,at%nat), stat=istat)
call memocc(istat, potmat, 'potmat', subname)


! Count how many locregs each process handles
ilrold=-1
nlocregOnMPI=0
do iorb=1,orbs%norbp
    iiorb=orbs%isorb+iorb
    ilr=orbs%inwhichlocreg(iiorb)
    !if(ilr>ilrold) then
    if(ilr/=ilrold) then
        nlocregOnMPI=nlocregOnMPI+1
    end if
    ilrold=ilr
end do
allocate(potmatsmall(orbs%norb,orbs%norb,nlocregOnMPI), stat=istat)
call memocc(istat, potmatsmall, 'potmatsmall', subname)



  call allocateSendBufferOrtho(comon, subname)
  call allocateRecvBufferOrtho(comon, subname)
  ! Extract the overlap region from the orbitals phi and store them in comon%sendBuf.
  !if(nit>0) then
  !    call extractOrbital3(iproc, nproc, orbs, orbs%npsidim, orbs%inWhichLocreg, lzd, op, lphi, comon%nsendBuf, comon%sendBuf)
  !    call postCommsOverlapNew(iproc, nproc, orbs, op, lzd, lphi, comon, tt1, tt2)
  !end if

  energyconf_trial=0.d0 !just to initialize this variable and make the compiler happy
  lstep=0.d0 !just to initialize this variable and make the compiler happy
  lstep_optimal=0.d0 !just to initialize this variable and make the compiler happy
  energyconf_der0=0.d0 !just to initialize this variable and make the compiler happy

  time_convol=0.d0
  time_lincomb=0.d0
  time_commun=0.d0
  time_linalg=0.d0
  time_exponential=0.d0
  time_matrixmodification=0.d0
  time_matrixElements=0.d0
  t1_tot=mpi_wtime()
  innerLoop: do it=1,nit

  !write(*,*) '1: iproc, associated(comon%recvbuf)', iproc, associated(comon%recvbuf)

      t1=mpi_wtime()
      call apply_orbitaldependent_potential(iproc, nproc, at, orbs, lzd, rxyz, &
           confdatarr, hx, lphi, -1, lvphi)
      t2=mpi_wtime()
      time_convol=time_convol+t2-t1

      t1=mpi_wtime()
      !allocate(ttmat(lin%orbs%norb,lin%orbs%norb))
      !call collectnew(iproc, nproc, comon, lin%mad,lin%op, lin%orbs, input, lin%lzd, comon%nsendbuf, &
      !     comon%sendbuf, comon%nrecvbuf, comon%recvbuf, ttmat, tt3, tt4, tt5)
      !deallocate(ttmat)
      !write(*,*) '2: iproc, associated(comon%recvbuf)', iproc, associated(comon%recvbuf)
      t2=mpi_wtime()
      time_commun=time_commun+t2-t1

      t1=mpi_wtime()
      !call getMatrixElements2(iproc, nproc, lin%lzd, lin%lb%orbs, lin%lb%op, lin%lb%comon, lphi, lvphi, lin%mad, Kmat)
      !call deallocateRecvBufferOrtho(comon, subname)
      !call deallocateSendBufferOrtho(comon, subname)
      call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lvphi, mad, Kmat)
      !do iorb=1,orbs%norb
      !    do jorb=1,orbs%norb
      !        if(iproc==0) write(66,*) iorb,jorb,Kmat(jorb,iorb)
      !    end do
      !end do
      !call allocateSendBufferOrtho(comon, subname)
      !call allocateRecvBufferOrtho(comon, subname)
      !write(*,*) '3: iproc, associated(comon%recvbuf)', iproc, associated(comon%recvbuf)
      t2=mpi_wtime()
      time_matrixElements=time_matrixElements+t2-t1


      !!if(newgradient) then
      !!    call get_potential_matrices(iproc, nproc, at, orbs, lzd, op, comon, mad, rxyz, &
      !!         confdatarr, hx, lphi, potmat)
      !!    !call get_potential_matrices_new(iproc, nproc, lin, at, input, orbs, lzd, op, comon, rxyz, lphi, &
      !!    !     nlocregOnMPI, potmatsmall)
      !!end if


      !!if(.not.newgradient) then
          !energyconf_0=ddot(orbs%npsidim, lphi(1), 1, lvphi(1), 1)
          !call mpiallred(energyconf_0, 1, mpi_sum, mpi_comm_world, ierr)
          energyconf_0=0.d0
          do iorb=1,orbs%norb
              energyconf_0 = energyconf_0 + Kmat(iorb,iorb)
          end do
      !!else
      !!    energyconf_0=0.d0
      !!    do iorb=1,orbs%norb
      !!        do jorb=1,orbs%norb
      !!            energyconf_0 = energyconf_0 + kernel(jorb,iorb)*Kmat(jorb,iorb)
      !!            !energyconf_0 = energyconf_0 + kernel(jorb,iorb)*Kmat(iorb,jorb)
      !!        end do
      !!    end do
      !!end if
      if(iproc==0) write(*,'(a,i6,3es20.10,2es17.7)') &
                   'it, energyconf_0, energyvonf_trial, energyconf_der0, lstep, lstep_optimal', &
                   it, energyconf_0, energyconf_trial, energyconf_der0, lstep, lstep_optimal

      t1=mpi_wtime()
      !!if(.not.newgradient) then
          ! Construct antisymmtric matrix Gmat
          do iorb=1,orbs%norb
              do jorb=1,orbs%norb
                  gmat(jorb,iorb)=2.d0*(Kmat(jorb,iorb)-Kmat(iorb,jorb))
              end do
          end do 
      !!else
      !!    !!! THIS IS THE OLD VERSION #############################################################################
      !!    do iorb=1,orbs%norb
      !!        ilr=orbs%inwhichlocreg(iorb)
      !!        do jorb=1,orbs%norb
      !!            jlr=orbs%inwhichlocreg(jorb)
      !!            tt=0.d0
      !!            do lorb=1,orbs%norb
      !!                tt = tt + kernel(jorb,lorb)*Kmat(lorb,iorb) - kernel(iorb,lorb)*Kmat(lorb,jorb) + &
      !!                          kernel(jorb,lorb)*potmat(lorb,iorb,jlr) - kernel(iorb,lorb)*potmat(lorb,jorb,ilr)
      !!            end do
      !!            gmat(jorb,iorb)=-tt
      !!            !if(iproc==0) then
      !!            !    write(77,*) iorb, jorb, gmat(jorb,iorb)
      !!            !end if
      !!        end do
      !!    end do 
      !!    ! ########################################################################################################
      !!    !!! THIS IS THE NEW VERSION
      !!    !!gmat=0.d0
      !!    !!ii=0
      !!    !!ilrold=-1
      !!    !!do iorb=1,orbs%norbp
      !!    !!    iiorb=orbs%isorb+iorb
      !!    !!    ilr=orbs%inwhichlocreg(iiorb)
      !!    !!    if(ilr>ilrold) then
      !!    !!        ii=ii+1
      !!    !!    end if
      !!    !!    do jorb=1,orbs%norb
      !!    !!        jlr=orbs%inwhichlocreg(jorb)
      !!    !!        tt=0.d0
      !!    !!        do lorb=1,orbs%norb
      !!    !!            !tt = tt + kernel(jorb,lorb)*Kmat(lorb,iiorb) - kernel(iiorb,lorb)*Kmat(lorb,jorb) + &
      !!    !!            !          - kernel(iiorb,lorb)*potmat(lorb,jorb,ilr)
      !!    !!            tt = tt + kernel(jorb,lorb)*Kmat(lorb,iiorb) - kernel(iiorb,lorb)*Kmat(lorb,jorb) + &
      !!    !!                      - kernel(iiorb,lorb)*potmatsmall(lorb,jorb,ii)
      !!    !!        end do
      !!    !!        gmat(jorb,iiorb)=-tt
      !!    !!    end do
      !!    !!    ilrold=ilr
      !!    !!end do 
      !!    !!do iorb=1,orbs%norb
      !!    !!    ilr=orbs%inwhichlocreg(iorb)
      !!    !!    jlrold=-1
      !!    !!    jj=0
      !!    !!    do jorb=1,orbs%norbp
      !!    !!        jjorb=orbs%isorb+jorb
      !!    !!        jlr=orbs%inwhichlocreg(jjorb)
      !!    !!        if(jlr>jlrold) then
      !!    !!            jj=jj+1
      !!    !!        end if
      !!    !!        tt=0.d0
      !!    !!        do lorb=1,orbs%norb
      !!    !!            !tt = tt + kernel(jjorb,lorb)*potmat(lorb,iorb,jlr)
      !!    !!            tt = tt + kernel(jjorb,lorb)*potmatsmall(lorb,iorb,jj)
      !!    !!        end do
      !!    !!        gmat(jjorb,iorb)=gmat(jjorb,iorb)-tt
      !!    !!        jlrold=jlr
      !!    !!    end do
      !!    !!end do 
      !!    !!call mpiallred(gmat(1,1), orbs%norb**2, mpi_sum, mpi_comm_world, ierr)
      !!    !!do iorb=1,orbs%norb
      !!    !!    do jorb=1,orbs%norb
      !!    !!        if(iproc==0) then
      !!    !!            write(77,*) iorb, jorb, gmat(jorb,iorb)
      !!    !!        end if
      !!    !!    end do
      !!    !!end do
      !!end if
      t2=mpi_wtime()
      time_matrixmodification=time_matrixmodification+t2-t1


      t1=mpi_wtime()
      !Build the complex matrix -iGmat
      do iorb=1,orbs%norb
          do jorb=1,orbs%norb
              gmatc(jorb,iorb)=cmplx(0.d0,-gmat(jorb,iorb),kind=8)
          end do
      end do 
      t2=mpi_wtime()
      time_matrixmodification=time_matrixmodification+t2-t1



      ! Diagonalize Gmatc
      t1=mpi_wtime()
      lwork=10*orbs%norb
      allocate(work(lwork), stat=istat) ! factor of 2 since it is assumed to be complex
      allocate(rwork(lwork), stat=istat)
      call zheev('v', 'l', orbs%norb, gmatc(1,1), orbs%norb, eval(1), work, lwork, rwork, info)
      if(info/=0) stop 'ERROR in zheev 1'
      deallocate(work)
      deallocate(rwork)
      t2=mpi_wtime()
      time_linalg=time_linalg+t2-t1


      ! Calculate step size
      if(it==1) then
          !!if(.not.newgradient) then
              lstep=5.d-2/(maxval(eval))
          !!else
          !!    lstep=1.d-4/(maxval(eval))
          !!end if
      else
          lstep=2.d0*lstep_optimal
          !lstep=1.d-3/(maxval(eval))
      end if

      t1=mpi_wtime()
      ! Calculate exp(-i*l*D) (with D diagonal matrix of eigenvalues).
      ! This is also a diagonal matrix, so only calculate the diagonal part.
      do iorb=1,orbs%norb
         ttc=cmplx(0.d0,-lstep*eval(iorb),kind=8)
         expD_cmplx(iorb)=(0.d0,0.d0)
          do k=0,50
              expD_cmplx(iorb)=expD_cmplx(iorb)+ttc**k/dfactorial(k)
          end do
      end do
      t2=mpi_wtime()
      time_exponential=time_exponential+t2-t1

      t1=mpi_wtime()
      do iorb=1,orbs%norb
          do jorb=1,orbs%norb
              if(iorb==jorb) then
                  tempmatc(jorb,iorb,1)=expD_cmplx(iorb)
              else
                  tempmatc(jorb,iorb,1)=cmplx(0.d0,0.d0,kind=8)
              end if
          end do
      end do
      t2=mpi_wtime()
      time_matrixmodification=time_matrixmodification+t2-t1

      t1=mpi_wtime()
      call zgemm('n', 'c', orbs%norb, orbs%norb, orbs%norb, (1.d0,0.d0), tempmatc(1,1,1), orbs%norb, &
           gmatc(1,1), orbs%norb, (0.d0,0.d0), tempmatc(1,1,2), orbs%norb)
      call zgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, (1.d0,0.d0), gmatc(1,1), orbs%norb, &
           tempmatc(1,1,2), orbs%norb, (0.d0,0.d0), omatc(1,1), orbs%norb)
      t2=mpi_wtime()
      time_linalg=time_linalg+t2-t1

      t1=mpi_wtime()
      ! Build new lphi
      do iorb=1,orbs%norb
          do jorb=1,orbs%norb
              tempmat3(jorb,iorb,1)=real(omatc(jorb,iorb))
          end do
      end do
      t2=mpi_wtime()
      time_matrixmodification=time_matrixmodification+t2-t1

      t1=mpi_wtime()
      !write(*,*) '5: iproc, associated(comon%recvbuf)', iproc, associated(comon%recvbuf)
      call build_new_linear_combinations(iproc, nproc, lzd, orbs, op, comon%nrecvbuf, &
           comon%recvbuf, tempmat3(1,1,1), .true., lphi)
      t2=mpi_wtime()
      time_lincomb=time_lincomb+t2-t1

      t1=mpi_wtime()
      call apply_orbitaldependent_potential(iproc, nproc, at, orbs, lzd, rxyz, &
           confdatarr, hx, lphi, -1, lvphi)
      t2=mpi_wtime()
      time_convol=time_convol+t2-t1

      !!if(.not.newgradient) then
          energyconf_trial=ddot(max(orbs%npsidim_orbs,orbs%npsidim_comp), lphi(1), 1, lvphi(1), 1)
          call mpiallred(energyconf_trial, 1, mpi_sum, mpi_comm_world, ierr)
      !!else

      !!    call dcopy(comon%nrecvbuf, comon%recvbuf, 1, recvbuf, 1)
      !!    !call extractOrbital3(iproc, nproc, orbs, orbs%npsidim, orbs%inWhichLocreg, lzd, op, lphi, comon%nsendBuf, comon%sendBuf)
      !!    !call postCommsOverlapNew(iproc, nproc, orbs, op, lzd, lphi, comon, tt1, tt2)
      !!    !allocate(ttmat(lin%orbs%norb,lin%orbs%norb))
      !!    !call collectnew(iproc, nproc, comon, lin%mad,lin%op, lin%orbs, input, lin%lzd, comon%nsendbuf, &
      !!    !     comon%sendbuf, comon%nrecvbuf, comon%recvbuf, ttmat, tt3, tt4, tt5)
      !!    !deallocate(ttmat)
      !!    call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, lphi, lvphi, mad, Kmat)
      !!    call dcopy(comon%nrecvbuf, recvbuf, 1, comon%recvbuf, 1)

      !!    energyconf_trial=0.d0
      !!    do iorb=1,orbs%norb
      !!        do jorb=1,orbs%norb
      !!            energyconf_trial = energyconf_trial + kernel(jorb,iorb)*Kmat(jorb,iorb)
      !!            !energyconf_trial = energyconf_trial + kernel(jorb,iorb)*Kmat(iorb,jorb)
      !!        end do
      !!    end do
      !!end if

      ! Calculate the gradient of the confinement
      energyconf_der0=0.d0
      do iorb=1,orbs%norb
          do jorb=1,orbs%norb
              energyconf_der0=energyconf_der0+gmat(jorb,iorb)**2
          end do
      end do
      energyconf_der0=-.5d0*energyconf_der0

      ! Calculate optimal step size
      lstep_optimal = -energyconf_der0*lstep**2/(2.d0*(energyconf_trial-energyconf_0-lstep*energyconf_der0))
      !!if(.not.newgradient) then
          lstep_optimal=min(lstep_optimal,lstep)
      !!else
      !!    if(lstep_optimal<0) then
      !!        lstep_optimal=lstep
      !!    else
      !!        lstep_optimal=min(lstep_optimal,lstep)
      !!    end if
      !!end if

      t1=mpi_wtime()
      ! Calculate exp(-i*l*D) (with D diagonal matrix of eigenvalues).
      ! This is also a diagonal matrix, so only calculate the diagonal part.
      do iorb=1,orbs%norb
         ttc=cmplx(0.d0,-lstep_optimal*eval(iorb),kind=8)
         expD_cmplx(iorb)=(0.d0,0.d0)
          do k=0,50
              expD_cmplx(iorb)=expD_cmplx(iorb)+ttc**k/dfactorial(k)
          end do
      end do
      t2=mpi_wtime()
      time_exponential=time_exponential+t2-t1

      t1=mpi_wtime()
      do iorb=1,orbs%norb
          do jorb=1,orbs%norb
              if(iorb==jorb) then
                  tempmatc(jorb,iorb,1)=expD_cmplx(iorb)
              else
                  tempmatc(jorb,iorb,1)=cmplx(0.d0,0.d0,kind=8)
              end if
          end do
      end do
      t2=mpi_wtime()
      time_matrixmodification=time_matrixmodification+t2-t1

      t1=mpi_wtime()
      call zgemm('n', 'c', orbs%norb, orbs%norb, orbs%norb, (1.d0,0.d0), tempmatc(1,1,1), orbs%norb, &
           gmatc(1,1), orbs%norb, (0.d0,0.d0), tempmatc(1,1,2), orbs%norb)
      call zgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, (1.d0,0.d0), gmatc(1,1), orbs%norb, &
           tempmatc(1,1,2), orbs%norb, (0.d0,0.d0), omatc(1,1), orbs%norb)
      t2=mpi_wtime()
      time_linalg=time_linalg+t2-t1


      ! Build new lphi
      do iorb=1,orbs%norb
          do jorb=1,orbs%norb
              tempmat3(jorb,iorb,1)=real(omatc(jorb,iorb))
          end do
      end do
      t1=mpi_wtime()
      call build_new_linear_combinations(iproc, nproc, lzd, orbs, op, comon%nrecvbuf, &
           comon%recvbuf, tempmat3(1,1,1), .true., lphi)
      t2=mpi_wtime()
      time_lincomb=time_lincomb+t2-t1


      !if(it<nit) then
      !    call extractOrbital3(iproc, nproc, orbs, orbs%npsidim, orbs%inWhichLocreg, lzd, op, lphi, comon%nsendBuf, comon%sendBuf)
      !    call postCommsOverlapNew(iproc, nproc, orbs, op, lzd, lphi, comon, tt1, tt2)
      !end if


  end do innerLoop

  t2_tot=mpi_wtime()
  time_tot=t2_tot-t1_tot
  call mpiallred(time_convol, 1, mpi_max, mpi_comm_world, ierr)
  call mpiallred(time_commun, 1, mpi_max, mpi_comm_world, ierr)
  call mpiallred(time_lincomb, 1, mpi_max, mpi_comm_world, ierr)
  call mpiallred(time_linalg, 1, mpi_max, mpi_comm_world, ierr)
  call mpiallred(time_matrixmodification, 1, mpi_max, mpi_comm_world, ierr)
  call mpiallred(time_exponential, 1, mpi_max, mpi_comm_world, ierr)
  call mpiallred(time_matrixelements, 1, mpi_max, mpi_comm_world, ierr)
  call mpiallred(time_tot, 1, mpi_max, mpi_comm_world, ierr)
  !time_convol=time_convol/dble(nproc)
  !time_commun=time_commun/dble(nproc)
  !time_lincomb=time_lincomb/dble(nproc)
  !time_linalg=time_linalg/dble(nproc)
  !time_matrixmodification=time_matrixmodification/dble(nproc)
  !time_exponential_=time_exponential/dble(nproc)
  !time_tot=time_tot/dble(nproc)
  !!if(iproc==0) then
  !!    write(*,'(a,es16.6)') 'total time: ',time_tot
  !!    write(*,'(a,es15.7,a,f5.2,a)') 'convolutions: ',time_convol,' (',time_convol/time_tot*100.d0,'%)'
  !!    write(*,'(a,es15.7,a,f5.2,a)') 'linear combinations: ',time_lincomb,' (',time_lincomb/time_tot*100.d0,'%)'
  !!    write(*,'(a,es15.7,a,f5.2,a)') 'communication: ',time_commun,' (',time_commun/time_tot*100.d0,'%)'
  !!    write(*,'(a,es15.7,a,f5.2,a)') 'linear algebra: ',time_linalg,' (',time_linalg/time_tot*100.d0,'%)'
  !!    write(*,'(a,es15.7,a,f5.2,a)') 'matrix modification: ',time_matrixmodification, &
  !!                                   ' (',time_matrixmodification/time_tot*100.d0,'%)'
  !!    write(*,'(a,es15.7,a,f5.2,a)') 'building exponential: ',time_exponential,' (',time_exponential/time_tot*100.d0,'%)'
  !!    write(*,'(a,es15.7,a,f5.2,a)') 'matrix elements ',time_matrixelements,' (',time_matrixelements/time_tot*100.d0,'%)'
  !!end if


  iall=-product(shape(gmat))*kind(gmat)
  deallocate(gmat, stat=istat)
  call memocc(istat, iall, 'gmat', subname)
  iall=-product(shape(gmatc))*kind(gmatc)
  deallocate(gmatc, stat=istat)
  !call memocc(istat, iall, 'gmatc', subname)
  iall=-product(shape(omatc))*kind(omatc)
  deallocate(omatc, stat=istat)
  !call memocc(istat, iall, 'omatc', subname)
  iall=-product(shape(tempmat3))*kind(tempmat3)
  deallocate(tempmat3, stat=istat)
  call memocc(istat, iall, 'tempmat3', subname)
  iall=-product(shape(eval))*kind(eval)
  deallocate(eval, stat=istat)
  call memocc(istat, iall, 'eval', subname)
  iall=-product(shape(expD_cmplx))*kind(expD_cmplx)
  deallocate(expD_cmplx, stat=istat)
  call memocc(istat, iall, 'expD_cmplx', subname)
  iall=-product(shape(tempmatc))*kind(tempmatc)
  deallocate(tempmatc, stat=istat)
  !call memocc(istat, iall, 'tempmatc', subname)
  iall=-product(shape(hamtrans))*kind(hamtrans)
  deallocate(hamtrans, stat=istat)
  call memocc(istat, iall, 'hamtrans', subname)
  iall=-product(shape(lphiovrlp))*kind(lphiovrlp)
  deallocate(lphiovrlp, stat=istat)
  call memocc(istat, iall, 'lphiovrlp', subname)
  iall=-product(shape(lvphi))*kind(lvphi)
  deallocate(lvphi, stat=istat)
  call memocc(istat, iall, 'lvphi', subname)
  iall=-product(shape(Kmat))*kind(Kmat)
  deallocate(Kmat, stat=istat)
  call memocc(istat, iall, 'Kmat', subname)
  iall=-product(shape(recvbuf))*kind(recvbuf)
  deallocate(recvbuf, stat=istat)
  call memocc(istat, iall, 'recvbuf', subname)

  iall=-product(shape(potmat))*kind(potmat)
  deallocate(potmat, stat=istat)
  call memocc(istat, iall, 'potmat', subname)
  iall=-product(shape(potmatsmall))*kind(potmatsmall)
  deallocate(potmatsmall, stat=istat)
  call memocc(istat, iall, 'potmatsmall', subname)

  call deallocateRecvBufferOrtho(comon, subname)
  call deallocateSendBufferOrtho(comon, subname)

end subroutine unitary_optimization




