!> @file
!! Derivatives of the basis set (linear version)
!! PB: Note that now they contain ONLY the derivatives, the original phi is not in phid
!! @author
!!    Copyright (C) 2011-2012 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


subroutine getDerivativeBasisFunctions(iproc, nproc, hgrid, lzd, lorbs, lborbs, comrp, nphi, phi, phid)
use module_base
use module_types
use module_interfaces, exceptThisOne => getDerivativeBasisFunctions
implicit none

! Calling arguments
integer,intent(in) :: iproc, nproc, nphi
real(kind=8),intent(in) :: hgrid
type(local_zone_descriptors),intent(in) :: lzd
type(orbitals_data),intent(in) :: lorbs, lborbs
type(p2pComms),intent(inout) :: comrp
real(kind=8),dimension(nphi),intent(in) :: phi !< Basis functions
real(kind=8),dimension(3*max(lorbs%npsidim_orbs,lorbs%npsidim_comp)),target,intent(inout) :: phid  !< Derivative basis functions

! Local variables
integer :: ist1_c, ist1_f, nf, istat, iall, iorb, jproc
integer :: ist0_c, istx_c, isty_c, istz_c, ist0_f, istx_f, isty_f, istz_f
integer :: jlr, offset, ilr, iiorb, i0, i1, i2, i3
real(kind=8),dimension(0:3),parameter :: scal=1.d0
real(kind=8),dimension(:),allocatable :: w_f1, w_f2, w_f3
real(kind=8),dimension(:),pointer :: phiLoc
real(kind=8),dimension(:,:,:),allocatable :: w_c, phix_c, phiy_c, phiz_c
real(kind=8),dimension(:,:,:,:),allocatable :: w_f, phix_f, phiy_f, phiz_f
logical :: repartition
character(len=*),parameter :: subname='getDerivativeBasisFunctions'


  ! Determine whether the orbitals must be redistributed after the calculation of the derivatives.
  ! If each processor has the same number of orbitals, this is never required.
  repartition=.false.
  !!do jproc=1,nproc-1
  !!   if(lorbs%norb_par(jproc,0)/=lorbs%norb_par(jproc-1,0)) then 
  !!       repartition=.true.
  !!       exit
  !!   end if
  !!end do

  if(repartition) then
      allocate(phiLoc(3*max(lorbs%npsidim_orbs,lorbs%npsidim_comp)), stat=istat)
      call memocc(istat, phiLoc, 'phiLoc', subname)
      call to_zero(3*max(lorbs%npsidim_orbs,lorbs%npsidim_comp), phiLoc(1))
  else
      phiLoc => phid
      call to_zero(3*max(lorbs%npsidim_orbs,lorbs%npsidim_comp), phid(1))
  end if


  ist1_c=1
  ist0_c=1
  ! Dimension of the first orbital on each process
  if(lorbs%norbp>0) then !otherwise this make no sense
      ilr=lorbs%inWhichLocreg(lorbs%isorb+1)
      offset=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f

      istx_c=1
      isty_c=offset+1
      istz_c=2*offset+1
  end if

  do iorb=1,lorbs%norbp

      !ilr=lin%orbs%inWhichLocregp(iorb)
      iiorb=lorbs%isorb+iorb
      ilr=lorbs%inWhichLocreg(iiorb)
      call allocateWorkarrays()

      ist1_f=ist1_c+lzd%llr(ilr)%wfd%nvctr_c
      ! ist0: start index of the orginal phi
      ist0_f=ist0_c+lzd%llr(ilr)%wfd%nvctr_c
      ! istx: start index of the derivative with respect to x
      istx_f=istx_c+lzd%llr(ilr)%wfd%nvctr_c
      ! isty: start index of the derivative with respect to y
      isty_f=isty_c+lzd%llr(ilr)%wfd%nvctr_c
      ! istz: start index of the derivative with respect to z
      istz_f=istz_c+lzd%llr(ilr)%wfd%nvctr_c

      ! Uncompress the wavefunction.
      !phi(ist1_f:ist1_f+7*lzd%llr(ilr)%wfd%nvctr_f-1)=0.d0
      !phi(ist1_c:ist1_c+lzd%llr(ilr)%wfd%nvctr_c-1)=0.d0
      call uncompress_forstandard(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
           lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, & 
           lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3,  &
           lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc,  &
           lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
           lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
           lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)),  &
           scal, phi(ist1_c), phi(ist1_f), w_c, w_f, w_f1, w_f2, w_f3)


      call createDerivativeBasis(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
           lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, &
           lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3,  &
           hgrid, lzd%llr(ilr)%bounds%kb%ibyz_c, lzd%llr(ilr)%bounds%kb%ibxz_c, lzd%llr(ilr)%bounds%kb%ibxy_c, &
           lzd%llr(ilr)%bounds%kb%ibyz_f, lzd%llr(ilr)%bounds%kb%ibxz_f, lzd%llr(ilr)%bounds%kb%ibxy_f, &
           w_c, w_f, w_f1, w_f2, w_f3, phix_c, phix_f, phiy_c, phiy_f, phiz_c, phiz_f)


      ! Copy phi to phiLoc
      !call dcopy(lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f, phi(ist1_c), 1, phiLoc(ist0_c), 1)
      !ist0_c = ist0_c + 4*(lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f)
      ist1_c = ist1_c + lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f

      ! Compress the x wavefunction.
      call compress_forstandard(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
           lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
           lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, &
           lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc, &
           lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
           lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
           lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)),  &
           scal, phix_c, phix_f, phiLoc(istx_c), phiLoc(istx_f))
      if(iorb<lorbs%norbp) then
          !jlr=lin%orbs%inWhichLocregp(iorb+1)
          jlr=lorbs%inWhichLocreg(iiorb+1)
          istx_c = istx_c + 3*(lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f)
      end if

      ! Compress the y wavefunction.
      call compress_forstandard(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
           lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
           lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, &
           lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc, &
           lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
           lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
           lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)),  &
           scal, phiy_c, phiy_f, phiLoc(isty_c), phiLoc(isty_f))
      if(iorb<lorbs%norbp) then
          !jlr=lin%orbs%inWhichLocregp(iorb+1)
          jlr=lorbs%inWhichLocreg(iiorb+1)
          isty_c = isty_c + 2*(lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f) + &
              (lzd%llr(jlr)%wfd%nvctr_c + 7*lzd%llr(jlr)%wfd%nvctr_f)
      end if

      ! Compress the z wavefunction.
      call compress_forstandard(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
           lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
           lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, &
           lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc, &
           lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
           lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
           lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)),  &
           scal, phiz_c, phiz_f, phiLoc(istz_c), phiLoc(istz_f))
      if(iorb<lorbs%norbp) then
          !jlr=lin%orbs%inWhichLocregp(iorb+1)
          jlr=lorbs%inWhichLocreg(iiorb+1)
          istz_c = istz_c + lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f + &
              2*(lzd%llr(jlr)%wfd%nvctr_c + 7*lzd%llr(jlr)%wfd%nvctr_f)
      end if

      call deallocateWorkarrays()

  end do


  if(repartition) then
      ! Communicate the orbitals to meet the partition.
      !!call postCommsRepartition(iproc, nproc, lorbs, comrp, size(phiLoc), phiLoc, size(phid), phid)
      !!call gatherDerivativeOrbitals(iproc, nproc, lorbs, comrp)
      call post_p2p_communication(iproc, nproc, size(phiLoc), phiLoc, size(phid), phid, comrp)
      call wait_p2p_communication(iproc, nproc, comrp)

      iall=-product(shape(phiLoc))*kind(phiLoc)
      deallocate(phiLoc, stat=istat)
      call memocc(istat, iall, 'phiLoc', subname)
  end if
  


contains

  subroutine allocateWorkarrays()
  
    ! THIS IS COPIED FROM allocate_work_arrays. Works only for free boundary.
    nf=(lzd%llr(ilr)%d%nfu1-lzd%llr(ilr)%d%nfl1+1)*(lzd%llr(ilr)%d%nfu2-lzd%llr(ilr)%d%nfl2+1)* &
       (lzd%llr(ilr)%d%nfu3-lzd%llr(ilr)%d%nfl3+1)

    ! Allocate work arrays
    allocate(w_c(0:lzd%llr(ilr)%d%n1,0:lzd%llr(ilr)%d%n2,0:lzd%llr(ilr)%d%n3+ndebug), stat=istat)
    call memocc(istat, w_c, 'w_c', subname)
    !!w_c=0.d0
    call to_zero((lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)*(lzd%llr(ilr)%d%n3+1), w_c(0,0,0))

    allocate(w_f(7,lzd%llr(ilr)%d%nfl1:lzd%llr(ilr)%d%nfu1,lzd%llr(ilr)%d%nfl2:lzd%llr(ilr)%d%nfu2, &
                 lzd%llr(ilr)%d%nfl3:lzd%llr(ilr)%d%nfu3+ndebug), stat=istat)
    call memocc(istat, w_f, 'w_f', subname)
    !!w_f=0.d0
    call to_zero(7*nf, w_f(1,lzd%llr(ilr)%d%nfl1,lzd%llr(ilr)%d%nfl2,lzd%llr(ilr)%d%nfl3))

  
    allocate(w_f1(nf+ndebug), stat=istat)
    call memocc(istat, w_f1, 'w_f1', subname)
    !!w_f1=0.d0
    call to_zero(nf, w_f1(1))
    
    allocate(w_f2(nf+ndebug), stat=istat)
    call memocc(istat, w_f2, 'w_f2', subname)
    !!w_f2=0.d0
    call to_zero(nf, w_f2(1))

    allocate(w_f3(nf+ndebug), stat=istat)
    call memocc(istat, w_f3, 'w_f3', subname)
    !!w_f3=0.d0
    call to_zero(nf, w_f3(1))
  
  
    allocate(phix_f(7,lzd%llr(ilr)%d%nfl1:lzd%llr(ilr)%d%nfu1,lzd%llr(ilr)%d%nfl2:lzd%llr(ilr)%d%nfu2, &
                    lzd%llr(ilr)%d%nfl3:lzd%llr(ilr)%d%nfu3), stat=istat)
    call memocc(istat, phix_f, 'phix_f', subname)
    !!phix_f=0.d0
    call to_zero(7*nf, phix_f(1,lzd%llr(ilr)%d%nfl1,lzd%llr(ilr)%d%nfl2,lzd%llr(ilr)%d%nfl3))

    allocate(phix_c(0:lzd%llr(ilr)%d%n1,0:lzd%llr(ilr)%d%n2,0:lzd%llr(ilr)%d%n3), stat=istat)
    call memocc(istat, phix_c, 'phix_c', subname)
    !!phix_c=0.d0
    call to_zero((lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)*(lzd%llr(ilr)%d%n3+1), phix_c(0,0,0))

    allocate(phiy_f(7,lzd%llr(ilr)%d%nfl1:lzd%llr(ilr)%d%nfu1,lzd%llr(ilr)%d%nfl2:lzd%llr(ilr)%d%nfu2, &
                    lzd%llr(ilr)%d%nfl3:lzd%llr(ilr)%d%nfu3), stat=istat)
    call memocc(istat, phiy_f, 'phiy_f', subname)
    !!phiy_f=0.d0
    call to_zero(7*nf, phiy_f(1,lzd%llr(ilr)%d%nfl1,lzd%llr(ilr)%d%nfl2,lzd%llr(ilr)%d%nfl3))

    allocate(phiy_c(0:lzd%llr(ilr)%d%n1,0:lzd%llr(ilr)%d%n2,0:lzd%llr(ilr)%d%n3), stat=istat)
    call memocc(istat, phiy_c, 'phiy_c', subname)
    !!phiy_c=0.d0
    call to_zero((lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)*(lzd%llr(ilr)%d%n3+1), phiy_c(0,0,0))

    allocate(phiz_f(7,lzd%llr(ilr)%d%nfl1:lzd%llr(ilr)%d%nfu1,lzd%llr(ilr)%d%nfl2:lzd%llr(ilr)%d%nfu2, &
                    lzd%llr(ilr)%d%nfl3:lzd%llr(ilr)%d%nfu3), stat=istat)
    call memocc(istat, phiz_f, 'phiz_f', subname)
    !!phiz_f=0.d0
    call to_zero(7*nf, phiz_f(1,lzd%llr(ilr)%d%nfl1,lzd%llr(ilr)%d%nfl2,lzd%llr(ilr)%d%nfl3))

    allocate(phiz_c(0:lzd%llr(ilr)%d%n1,0:lzd%llr(ilr)%d%n2,0:lzd%llr(ilr)%d%n3), stat=istat)
    call memocc(istat, phiz_c, 'phiz_c', subname)
    !!phiz_c=0.d0
    call to_zero((lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)*(lzd%llr(ilr)%d%n3+1), phiz_c(0,0,0))
  
  end subroutine allocateWorkarrays


  subroutine deallocateWorkarrays

    iall=-product(shape(w_c))*kind(w_c)
    deallocate(w_c, stat=istat)
    call memocc(istat, iall, 'w_c', subname)

    iall=-product(shape(w_f))*kind(w_f)
    deallocate(w_f, stat=istat)
    call memocc(istat, iall, 'w_f', subname)

    iall=-product(shape(w_f1))*kind(w_f1)
    deallocate(w_f1, stat=istat)
    call memocc(istat, iall, 'w_f1', subname)

    iall=-product(shape(w_f2))*kind(w_f2)
    deallocate(w_f2, stat=istat)
    call memocc(istat, iall, 'w_f2', subname)

    iall=-product(shape(w_f3))*kind(w_f3)
    deallocate(w_f3, stat=istat)
    call memocc(istat, iall, 'w_f3', subname)

    iall=-product(shape(phix_f))*kind(phix_f)
    deallocate(phix_f, stat=istat)
    call memocc(istat, iall, 'phix_f', subname)

    iall=-product(shape(phix_c))*kind(phix_c)
    deallocate(phix_c, stat=istat)
    call memocc(istat, iall, 'phix_c', subname)

    iall=-product(shape(phiy_f))*kind(phiy_f)
    deallocate(phiy_f, stat=istat)
    call memocc(istat, iall, 'phiy_f', subname)

    iall=-product(shape(phiy_c))*kind(phiy_c)
    deallocate(phiy_c, stat=istat)
    call memocc(istat, iall, 'phiy_c', subname)

    iall=-product(shape(phiz_f))*kind(phiz_f)
    deallocate(phiz_f, stat=istat)
    call memocc(istat, iall, 'phiz_f', subname)

    iall=-product(shape(phiz_c))*kind(phiz_c)
    deallocate(phiz_c, stat=istat)
    call memocc(istat, iall, 'phiz_c', subname)

  end subroutine deallocateWorkarrays

end subroutine getDerivativeBasisFunctions






subroutine initializeRepartitionOrbitals(iproc, nproc, tag, lorbs, llborbs, lzd, comrp)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in) :: iproc, nproc
integer,intent(inout) :: tag
type(orbitals_data),intent(in) :: lorbs, llborbs
type(local_zone_descriptors),intent(in) :: lzd
type(p2pComms),intent(out) :: comrp

! Local variables
integer :: jproc, jorb, kproc, korb, istat, mpidest, mpisource, istdest, istsource, ncount, jlr, klr, norbdest, jjorb, iall
integer :: isend, irecv, p2p_tag
integer,dimension(:,:,:),allocatable :: move
character(len=*),parameter :: subname='initializeRepartitionOrbitals'

call timing(iproc,'init_repart   ','ON')

! Nullify at the beginning
call nullify_p2pComms(comrp)

! To which position has the orbital to be moved:
!  - move(1,i,j)=k -> orbital i on process j has to be sent to process k
!  - move(2,i,j)=l -> orbital i on process j has to be sent to position l (orbital number l)

allocate(move(2,3*maxval(lorbs%norb_par(:,0)),0:nproc-1), stat=istat)
call memocc(istat, move, 'move', subname)

kproc=0
korb=0
do jproc=0,nproc-1
    do jorb=1,3*lorbs%norb_par(jproc,0)
        korb=korb+1
        if(korb>llborbs%norb_par(kproc,0)) then
            kproc=kproc+1
            korb=1
        end if
        move(1,jorb,jproc)=kproc
        move(2,jorb,jproc)=korb
    end do
end do


allocate(comrp%noverlaps(0:nproc-1), stat=istat)
call memocc(istat, comrp%noverlaps, 'comrp%noverlaps', subname)

allocate(comrp%comarr(6,3*maxval(lorbs%norb_par(:,0)),0:nproc-1), stat=istat)
call memocc(istat, comrp%comarr, 'comrp%comarr', subname)

! Determine the indices of starting and receive buffer.
isend=0
irecv=0
comrp%noverlaps=0
do jproc=0,nproc-1
    istsource=1
    do jorb=1,3*lorbs%norb_par(jproc,0)
        jjorb=ceiling(dble(jorb)/3.d0)
        jlr=lorbs%inWhichLocreg(jjorb+lorbs%isorb_par(jproc))
        mpisource=jproc
        ncount=lzd%llr(jlr)%wfd%nvctr_c+7*lzd%llr(jlr)%wfd%nvctr_f
        mpidest=move(1,jorb,jproc)
        comrp%noverlaps(mpidest)=comrp%noverlaps(mpidest)+1
        norbdest=move(2,jorb,jproc)
        istdest=1
        do korb=1,norbdest-1
            klr=llborbs%inWhichLocreg(korb+llborbs%isorb_par(mpidest))
            istdest=istdest+lzd%llr(klr)%wfd%nvctr_c+7*lzd%llr(klr)%wfd%nvctr_f
        end do
        !tag=tag+1
        tag=p2p_tag(mpidest)
        !write(*,'(7(a,i0))') 'init on iproc=',iproc,': process ',mpisource,' sends ',ncount,' elements from position ',istsource,' to position ',istdest,' on process ',mpidest,'; tag=',tag
        if(iproc==mpisource) then
            isend=isend+1
        end if
        if(iproc==mpidest) then
            irecv=irecv+1
            !write(*,'(a,2(i0,a))') 'INIT: process ',iproc,' receives data from process ',mpisource,'.'
        end if
        !call setCommsParameters(mpisource, mpidest, istsource, istdest, ncount, tag, comrp%comarr(1,jorb,jproc))
        call setCommsParameters(mpisource, mpidest, istsource, istdest, ncount, &
             tag, comrp%comarr(1,comrp%noverlaps(mpidest),mpidest))
        istsource=istsource+ncount
    end do
end do
comrp%nsend=isend
comrp%nrecv=irecv
if(irecv/=comrp%noverlaps(iproc)) stop 'irecv/=comrp%noverlaps(iproc)'

! Allocate the requests for the point to point communication.
allocate(comrp%requests(max(comrp%nsend,comrp%nrecv),2), stat=istat)
call memocc(istat, comrp%requests, 'comrp%requests', subname)
      

iall=-product(shape(move))*kind(move)
deallocate(move, stat=istat)
call memocc(istat, iall, 'move', subname)

! To indicate that no communication is going on.
comrp%communication_complete=.true.
comrp%messages_posted=.false.


call timing(iproc,'init_repart   ','OF')

end subroutine initializeRepartitionOrbitals





subroutine get_derivative_supportfunctions(ndim, hgrid, lzd, lorbs, phi, phid)
use module_base
use module_types
use module_interfaces, except_this_one => get_derivative_supportfunctions
implicit none

! Calling arguments
integer,intent(in):: ndim
real(kind=8),intent(in) :: hgrid
type(local_zone_descriptors),intent(in) :: lzd
type(orbitals_data),intent(in) :: lorbs
real(kind=8),dimension(ndim),intent(in) :: phi !< Basis functions
real(kind=8),dimension(3*ndim),intent(inout) :: phid  !< Derivative basis functions

! Local variables
integer :: ist1_c, ist1_f, nf, istat, iall, iorb, jproc
integer :: ist0_c, istx_c, isty_c, istz_c, ist0_f, istx_f, isty_f, istz_f
integer :: jlr, offset, ilr, iiorb
real(kind=8),dimension(0:3),parameter :: scal=1.d0
real(kind=8),dimension(:),allocatable :: w_f1, w_f2, w_f3
real(kind=8),dimension(:),pointer :: phiLoc
real(kind=8),dimension(:,:,:),allocatable :: w_c, phix_c, phiy_c, phiz_c
real(kind=8),dimension(:,:,:,:),allocatable :: w_f, phix_f, phiy_f, phiz_f
logical :: repartition
character(len=*),parameter :: subname='get_derivative_supportfunctions'


  ! Determine whether the orbitals must be redistributed after the calculation of the derivatives.
  ! If each processor has the same number of orbitals, this is never required.
  repartition=.false.
  !!do jproc=1,nproc-1
  !!   if(lorbs%norb_par(jproc,0)/=lorbs%norb_par(jproc-1,0)) then 
  !!       repartition=.true.
  !!       exit
  !!   end if
  !!end do

  !!if(repartition) then
  !!    allocate(phiLoc(3*max(lorbs%npsidim_orbs,lorbs%npsidim_comp)), stat=istat)
  !!    call memocc(istat, phiLoc, 'phiLoc', subname)
  !!else
  !!    phiLoc => phid
  !!end if
 

  ist1_c=1
  ist0_c=1
  ! Dimension of the first orbital on each process
  if(lorbs%norbp>0) then !otherwise this make no sense
      ilr=lorbs%inWhichLocreg(lorbs%isorb+1)
      offset=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f

      istx_c=1
      isty_c=offset+1
      istz_c=2*offset+1
  end if

  do iorb=1,lorbs%norbp

      !ilr=lin%orbs%inWhichLocregp(iorb)
      iiorb=lorbs%isorb+iorb
      ilr=lorbs%inWhichLocreg(iiorb)
      call allocateWorkarrays()

      ist1_f=ist1_c+lzd%llr(ilr)%wfd%nvctr_c
      ! ist0: start index of the orginal phi
      ist0_f=ist0_c+lzd%llr(ilr)%wfd%nvctr_c
      ! istx: start index of the derivative with respect to x
      istx_f=istx_c+lzd%llr(ilr)%wfd%nvctr_c
      ! isty: start index of the derivative with respect to y
      isty_f=isty_c+lzd%llr(ilr)%wfd%nvctr_c
      ! istz: start index of the derivative with respect to z
      istz_f=istz_c+lzd%llr(ilr)%wfd%nvctr_c

      ! Uncompress the wavefunction.
      !phi(ist1_f:ist1_f+7*lzd%llr(ilr)%wfd%nvctr_f-1)=0.d0
      !phi(ist1_c:ist1_c+lzd%llr(ilr)%wfd%nvctr_c-1)=0.d0
      call uncompress_forstandard(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
           lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, & 
           lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3,  &
           lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc,  &
           lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
           lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
           lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)),  &
           scal, phi(ist1_c), phi(ist1_f), w_c, w_f, w_f1, w_f2, w_f3)


      call createDerivativeBasis(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
           lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, &
           lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3,  &
           hgrid, lzd%llr(ilr)%bounds%kb%ibyz_c, lzd%llr(ilr)%bounds%kb%ibxz_c, lzd%llr(ilr)%bounds%kb%ibxy_c, &
           lzd%llr(ilr)%bounds%kb%ibyz_f, lzd%llr(ilr)%bounds%kb%ibxz_f, lzd%llr(ilr)%bounds%kb%ibxy_f, &
           w_c, w_f, w_f1, w_f2, w_f3, phix_c, phix_f, phiy_c, phiy_f, phiz_c, phiz_f)

      ! Copy phi to phid
      !call dcopy(lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f, phi(ist1_c), 1, phid(ist0_c), 1)
      !ist0_c = ist0_c + 4*(lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f)
      ist1_c = ist1_c + lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f

      ! Compress the x wavefunction.
      call compress_forstandard(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
           lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
           lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, &
           lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc, &
           lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
           lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
           lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)),  &
           scal, phix_c, phix_f, phid(istx_c), phid(istx_f))
      if(iorb<lorbs%norbp) then
          !jlr=lin%orbs%inWhichLocregp(iorb+1)
          jlr=lorbs%inWhichLocreg(iiorb+1)
          istx_c = istx_c + 3*(lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f)
      end if

      ! Compress the y wavefunction.
      call compress_forstandard(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
           lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
           lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, &
           lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc, &
           lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
           lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
           lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)),  &
           scal, phiy_c, phiy_f, phid(isty_c), phid(isty_f))
      if(iorb<lorbs%norbp) then
          !jlr=lin%orbs%inWhichLocregp(iorb+1)
          jlr=lorbs%inWhichLocreg(iiorb+1)
          isty_c = isty_c + 2*(lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f) + &
              (lzd%llr(jlr)%wfd%nvctr_c + 7*lzd%llr(jlr)%wfd%nvctr_f)
      end if

      ! Compress the z wavefunction.
      call compress_forstandard(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
           lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
           lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, &
           lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc, &
           lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
           lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
           lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)),  &
           scal, phiz_c, phiz_f, phid(istz_c), phid(istz_f))
      if(iorb<lorbs%norbp) then
          !jlr=lin%orbs%inWhichLocregp(iorb+1)
          jlr=lorbs%inWhichLocreg(iiorb+1)
          istz_c = istz_c + lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f + &
              2*(lzd%llr(jlr)%wfd%nvctr_c + 7*lzd%llr(jlr)%wfd%nvctr_f)
      end if

      call deallocateWorkarrays()

  end do


  !!if(repartition) then
  !!    ! Communicate the orbitals to meet the partition.
  !!    !!call postCommsRepartition(iproc, nproc, lorbs, comrp, size(phid), phid, size(phid), phid)
  !!    !!call gatherDerivativeOrbitals(iproc, nproc, lorbs, comrp)
  !!    call post_p2p_communication(iproc, nproc, size(phid), phid, size(phid), phid, comrp)
  !!    call wait_p2p_communication(iproc, nproc, comrp)

  !!    iall=-product(shape(phid))*kind(phid)
  !!    deallocate(phid, stat=istat)
  !!    call memocc(istat, iall, 'phid', subname)
  !!end if
  


contains

  subroutine allocateWorkarrays()
  
    ! THIS IS COPIED FROM allocate_work_arrays. Works only for free boundary.
    nf=(lzd%llr(ilr)%d%nfu1-lzd%llr(ilr)%d%nfl1+1)*(lzd%llr(ilr)%d%nfu2-lzd%llr(ilr)%d%nfl2+1)* &
       (lzd%llr(ilr)%d%nfu3-lzd%llr(ilr)%d%nfl3+1)

    ! Allocate work arrays
    allocate(w_c(0:lzd%llr(ilr)%d%n1,0:lzd%llr(ilr)%d%n2,0:lzd%llr(ilr)%d%n3+ndebug), stat=istat)
    call memocc(istat, w_c, 'w_c', subname)
    !!w_c=0.d0
    call to_zero((lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)*(lzd%llr(ilr)%d%n3+1), w_c(0,0,0))

    allocate(w_f(7,lzd%llr(ilr)%d%nfl1:lzd%llr(ilr)%d%nfu1,lzd%llr(ilr)%d%nfl2:lzd%llr(ilr)%d%nfu2, &
                 lzd%llr(ilr)%d%nfl3:lzd%llr(ilr)%d%nfu3+ndebug), stat=istat)
    call memocc(istat, w_f, 'w_f', subname)
    !!w_f=0.d0
    call to_zero(7*nf, w_f(1,lzd%llr(ilr)%d%nfl1,lzd%llr(ilr)%d%nfl2,lzd%llr(ilr)%d%nfl3))

  
    allocate(w_f1(nf+ndebug), stat=istat)
    call memocc(istat, w_f1, 'w_f1', subname)
    !!w_f1=0.d0
    call to_zero(nf, w_f1(1))
    
    allocate(w_f2(nf+ndebug), stat=istat)
    call memocc(istat, w_f2, 'w_f2', subname)
    !!w_f2=0.d0
    call to_zero(nf, w_f2(1))

    allocate(w_f3(nf+ndebug), stat=istat)
    call memocc(istat, w_f3, 'w_f3', subname)
    !!w_f3=0.d0
    call to_zero(nf, w_f3(1))
  
  
    allocate(phix_f(7,lzd%llr(ilr)%d%nfl1:lzd%llr(ilr)%d%nfu1,lzd%llr(ilr)%d%nfl2:lzd%llr(ilr)%d%nfu2, &
                    lzd%llr(ilr)%d%nfl3:lzd%llr(ilr)%d%nfu3), stat=istat)
    call memocc(istat, phix_f, 'phix_f', subname)
    !!phix_f=0.d0
    call to_zero(7*nf, phix_f(1,lzd%llr(ilr)%d%nfl1,lzd%llr(ilr)%d%nfl2,lzd%llr(ilr)%d%nfl3))

    allocate(phix_c(0:lzd%llr(ilr)%d%n1,0:lzd%llr(ilr)%d%n2,0:lzd%llr(ilr)%d%n3), stat=istat)
    call memocc(istat, phix_c, 'phix_c', subname)
    !!phix_c=0.d0
    call to_zero((lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)*(lzd%llr(ilr)%d%n3+1), phix_c(0,0,0))

    allocate(phiy_f(7,lzd%llr(ilr)%d%nfl1:lzd%llr(ilr)%d%nfu1,lzd%llr(ilr)%d%nfl2:lzd%llr(ilr)%d%nfu2, &
                    lzd%llr(ilr)%d%nfl3:lzd%llr(ilr)%d%nfu3), stat=istat)
    call memocc(istat, phiy_f, 'phiy_f', subname)
    !!phiy_f=0.d0
    call to_zero(7*nf, phiy_f(1,lzd%llr(ilr)%d%nfl1,lzd%llr(ilr)%d%nfl2,lzd%llr(ilr)%d%nfl3))

    allocate(phiy_c(0:lzd%llr(ilr)%d%n1,0:lzd%llr(ilr)%d%n2,0:lzd%llr(ilr)%d%n3), stat=istat)
    call memocc(istat, phiy_c, 'phiy_c', subname)
    !!phiy_c=0.d0
    call to_zero((lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)*(lzd%llr(ilr)%d%n3+1), phiy_c(0,0,0))

    allocate(phiz_f(7,lzd%llr(ilr)%d%nfl1:lzd%llr(ilr)%d%nfu1,lzd%llr(ilr)%d%nfl2:lzd%llr(ilr)%d%nfu2, &
                    lzd%llr(ilr)%d%nfl3:lzd%llr(ilr)%d%nfu3), stat=istat)
    call memocc(istat, phiz_f, 'phiz_f', subname)
    !!phiz_f=0.d0
    call to_zero(7*nf, phiz_f(1,lzd%llr(ilr)%d%nfl1,lzd%llr(ilr)%d%nfl2,lzd%llr(ilr)%d%nfl3))

    allocate(phiz_c(0:lzd%llr(ilr)%d%n1,0:lzd%llr(ilr)%d%n2,0:lzd%llr(ilr)%d%n3), stat=istat)
    call memocc(istat, phiz_c, 'phiz_c', subname)
    !!phiz_c=0.d0
    call to_zero((lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)*(lzd%llr(ilr)%d%n3+1), phiz_c(0,0,0))
  
  end subroutine allocateWorkarrays


  subroutine deallocateWorkarrays

    iall=-product(shape(w_c))*kind(w_c)
    deallocate(w_c, stat=istat)
    call memocc(istat, iall, 'w_c', subname)

    iall=-product(shape(w_f))*kind(w_f)
    deallocate(w_f, stat=istat)
    call memocc(istat, iall, 'w_f', subname)

    iall=-product(shape(w_f1))*kind(w_f1)
    deallocate(w_f1, stat=istat)
    call memocc(istat, iall, 'w_f1', subname)

    iall=-product(shape(w_f2))*kind(w_f2)
    deallocate(w_f2, stat=istat)
    call memocc(istat, iall, 'w_f2', subname)

    iall=-product(shape(w_f3))*kind(w_f3)
    deallocate(w_f3, stat=istat)
    call memocc(istat, iall, 'w_f3', subname)

    iall=-product(shape(phix_f))*kind(phix_f)
    deallocate(phix_f, stat=istat)
    call memocc(istat, iall, 'phix_f', subname)

    iall=-product(shape(phix_c))*kind(phix_c)
    deallocate(phix_c, stat=istat)
    call memocc(istat, iall, 'phix_c', subname)

    iall=-product(shape(phiy_f))*kind(phiy_f)
    deallocate(phiy_f, stat=istat)
    call memocc(istat, iall, 'phiy_f', subname)

    iall=-product(shape(phiy_c))*kind(phiy_c)
    deallocate(phiy_c, stat=istat)
    call memocc(istat, iall, 'phiy_c', subname)

    iall=-product(shape(phiz_f))*kind(phiz_f)
    deallocate(phiz_f, stat=istat)
    call memocc(istat, iall, 'phiz_f', subname)

    iall=-product(shape(phiz_c))*kind(phiz_c)
    deallocate(phiz_c, stat=istat)
    call memocc(istat, iall, 'phiz_c', subname)

  end subroutine deallocateWorkarrays

end subroutine get_derivative_supportfunctions
