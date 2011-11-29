subroutine getDerivativeBasisFunctions2(iproc, nproc, hgrid, Glr, lin, nphi, phi, phid)
use module_base
use module_types
use module_interfaces, exceptThisOne => getDerivativeBasisFunctions2
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, nphi
real(8),intent(in):: hgrid
type(locreg_descriptors),intent(in):: Glr
type(linearParameters),intent(inout):: lin
real(8),dimension(nphi),intent(in):: phi
real(8),dimension(lin%lb%orbs%npsidim),target,intent(out):: phid

! Local variables
integer:: ist1_c, ist1_f, ist2_c, ist2_f, nf, istat, iall, iorb, jproc, ierr
integer:: ist0_c, istx_c, isty_c, istz_c, ist0_f, istx_f, isty_f, istz_f, istLoc, istRoot
integer:: jjorb, jlr, jj, offset, ilr, jorb
real(8),dimension(0:3),parameter:: scal=1.d0
real(8),dimension(:),allocatable:: w_f1, w_f2, w_f3
real(8),dimension(:),pointer:: phiLoc
real(8),dimension(:,:,:),allocatable:: w_c, phix_c, phiy_c, phiz_c
real(8),dimension(:,:,:,:),allocatable:: w_f, phix_f, phiy_f, phiz_f
character(len=*),parameter:: subname='getDerivativeBasisFunctions'
integer,dimension(:),allocatable:: recvcounts, sendcounts, displs
logical:: repartition

integer:: i1, i2, i3
logical,dimension(:,:,:),allocatable:: logrid_c, logrid_f

  !!do jj=1,nphi
  !!    write(700+iproc,*) jj, phi(jj)
  !!end do

  ! Determine whether the orbitals must be redistributed after the calculation of the derivatives.
  ! If each orbital has the same number of orbitals, this is never required.
  repartition=.false.
  do jproc=1,nproc-1
     if(lin%orbs%norb_par(jproc,0)/=lin%orbs%norb_par(jproc-1,0)) then 
         repartition=.true.
         exit
     end if
  end do


  if(repartition) then
      allocate(phiLoc(4*lin%orbs%npsidim), stat=istat)
      call memocc(istat, phiLoc, 'phiLoc', subname)
  else
      phiLoc => phid
  end if
 

  ist1_c=1
  ist2_c=1
  ist0_c=1
  ! Dimension of the first orbital on each process
  ilr=lin%orbs%inWhichLocregp(1)
  offset=lin%lzd%llr(ilr)%wfd%nvctr_c+7*lin%lzd%llr(ilr)%wfd%nvctr_f

  istx_c=offset+1
  isty_c=2*offset+1
  istz_c=3*offset+1

  do iorb=1,lin%orbs%norbp

      ilr=lin%orbs%inWhichLocregp(iorb)
      call allocateWorkarrays()

      ist1_f=ist1_c+lin%lzd%llr(ilr)%wfd%nvctr_c
      ist2_f=ist2_c+lin%lzd%llr(ilr)%wfd%nvctr_c
      ! ist0: start index of the orginal phi
      ist0_f=ist0_c+lin%lzd%llr(ilr)%wfd%nvctr_c
      ! istx: start index of the derivative with respect to x
      istx_f=istx_c+lin%lzd%llr(ilr)%wfd%nvctr_c
      ! isty: start index of the derivative with respect to y
      isty_f=isty_c+lin%lzd%llr(ilr)%wfd%nvctr_c
      ! istz: start index of the derivative with respect to z
      istz_f=istz_c+lin%lzd%llr(ilr)%wfd%nvctr_c

      ! Uncompress the wavefunction.
      !phi(ist1_f:ist1_f+7*lin%lzd%llr(ilr)%wfd%nvctr_f-1)=0.d0
      !phi(ist1_c:ist1_c+lin%lzd%llr(ilr)%wfd%nvctr_c-1)=0.d0
      call uncompress_forstandard(lin%lzd%llr(ilr)%d%n1, lin%lzd%llr(ilr)%d%n2, lin%lzd%llr(ilr)%d%n3, &
           lin%lzd%llr(ilr)%d%nfl1, lin%lzd%llr(ilr)%d%nfu1, & 
           lin%lzd%llr(ilr)%d%nfl2, lin%lzd%llr(ilr)%d%nfu2, lin%lzd%llr(ilr)%d%nfl3, lin%lzd%llr(ilr)%d%nfu3,  &
           lin%lzd%llr(ilr)%wfd%nseg_c, lin%lzd%llr(ilr)%wfd%nvctr_c, lin%lzd%llr(ilr)%wfd%keyg, lin%lzd%llr(ilr)%wfd%keyv,  &
           lin%lzd%llr(ilr)%wfd%nseg_f, lin%lzd%llr(ilr)%wfd%nvctr_f, &
           lin%lzd%llr(ilr)%wfd%keyg(1,lin%lzd%llr(ilr)%wfd%nseg_c+min(1,lin%lzd%llr(ilr)%wfd%nseg_f)), &
           lin%lzd%llr(ilr)%wfd%keyv(lin%lzd%llr(ilr)%wfd%nseg_c+min(1,lin%lzd%llr(ilr)%wfd%nseg_f)),  &
           scal, phi(ist1_c), phi(ist1_f), w_c, w_f, w_f1, w_f2, w_f3)


      call createDerivativeBasis(lin%lzd%llr(ilr)%d%n1, lin%lzd%llr(ilr)%d%n2, lin%lzd%llr(ilr)%d%n3, &
           lin%lzd%llr(ilr)%d%nfl1, lin%lzd%llr(ilr)%d%nfu1, lin%lzd%llr(ilr)%d%nfl2, lin%lzd%llr(ilr)%d%nfu2, &
           lin%lzd%llr(ilr)%d%nfl3, lin%lzd%llr(ilr)%d%nfu3,  &
           hgrid, lin%lzd%llr(ilr)%bounds%kb%ibyz_c, lin%lzd%llr(ilr)%bounds%kb%ibxz_c, lin%lzd%llr(ilr)%bounds%kb%ibxy_c, &
           lin%lzd%llr(ilr)%bounds%kb%ibyz_f, lin%lzd%llr(ilr)%bounds%kb%ibxz_f, lin%lzd%llr(ilr)%bounds%kb%ibxy_f, &
           w_c, w_f, w_f1, w_f2, w_f3, phix_c, phix_f, phiy_c, phiy_f, phiz_c, phiz_f)

      ! Copy phi to phiLoc
      call dcopy(lin%lzd%llr(ilr)%wfd%nvctr_c+7*lin%lzd%llr(ilr)%wfd%nvctr_f, phi(ist1_c), 1, phiLoc(ist0_c), 1)
      ist0_c = ist0_c + 4*(lin%lzd%llr(ilr)%wfd%nvctr_c + 7*lin%lzd%llr(ilr)%wfd%nvctr_f)
      ist1_c = ist1_c + lin%lzd%llr(ilr)%wfd%nvctr_c + 7*lin%lzd%llr(ilr)%wfd%nvctr_f

      ! Compress the x wavefunction.
      call compress_forstandard(lin%lzd%llr(ilr)%d%n1, lin%lzd%llr(ilr)%d%n2, lin%lzd%llr(ilr)%d%n3, &
           lin%lzd%llr(ilr)%d%nfl1, lin%lzd%llr(ilr)%d%nfu1, &
           lin%lzd%llr(ilr)%d%nfl2, lin%lzd%llr(ilr)%d%nfu2, lin%lzd%llr(ilr)%d%nfl3, lin%lzd%llr(ilr)%d%nfu3, &
           lin%lzd%llr(ilr)%wfd%nseg_c, lin%lzd%llr(ilr)%wfd%nvctr_c, lin%lzd%llr(ilr)%wfd%keyg, lin%lzd%llr(ilr)%wfd%keyv, &
           lin%lzd%llr(ilr)%wfd%nseg_f, lin%lzd%llr(ilr)%wfd%nvctr_f, &
           lin%lzd%llr(ilr)%wfd%keyg(1,lin%lzd%llr(ilr)%wfd%nseg_c+min(1,lin%lzd%llr(ilr)%wfd%nseg_f)), &
           lin%lzd%llr(ilr)%wfd%keyv(lin%lzd%llr(ilr)%wfd%nseg_c+min(1,lin%lzd%llr(ilr)%wfd%nseg_f)),  &
           scal, phix_c, phix_f, phiLoc(istx_c), phiLoc(istx_f))
      if(iorb<lin%orbs%norbp) then
          jlr=lin%orbs%inWhichLocregp(iorb+1)
          istx_c = istx_c + 3*(lin%lzd%llr(ilr)%wfd%nvctr_c + 7*lin%lzd%llr(ilr)%wfd%nvctr_f) + &
              lin%lzd%llr(jlr)%wfd%nvctr_c + 7*lin%lzd%llr(jlr)%wfd%nvctr_f
      end if

      ! Compress the y wavefunction.
      call compress_forstandard(lin%lzd%llr(ilr)%d%n1, lin%lzd%llr(ilr)%d%n2, lin%lzd%llr(ilr)%d%n3, &
           lin%lzd%llr(ilr)%d%nfl1, lin%lzd%llr(ilr)%d%nfu1, &
           lin%lzd%llr(ilr)%d%nfl2, lin%lzd%llr(ilr)%d%nfu2, lin%lzd%llr(ilr)%d%nfl3, lin%lzd%llr(ilr)%d%nfu3, &
           lin%lzd%llr(ilr)%wfd%nseg_c, lin%lzd%llr(ilr)%wfd%nvctr_c, lin%lzd%llr(ilr)%wfd%keyg, lin%lzd%llr(ilr)%wfd%keyv, &
           lin%lzd%llr(ilr)%wfd%nseg_f, lin%lzd%llr(ilr)%wfd%nvctr_f, &
           lin%lzd%llr(ilr)%wfd%keyg(1,lin%lzd%llr(ilr)%wfd%nseg_c+min(1,lin%lzd%llr(ilr)%wfd%nseg_f)), &
           lin%lzd%llr(ilr)%wfd%keyv(lin%lzd%llr(ilr)%wfd%nseg_c+min(1,lin%lzd%llr(ilr)%wfd%nseg_f)),  &
           scal, phiy_c, phiy_f, phiLoc(isty_c), phiLoc(isty_f))
      if(iorb<lin%orbs%norbp) then
          jlr=lin%orbs%inWhichLocregp(iorb+1)
          isty_c = isty_c + 2*(lin%lzd%llr(ilr)%wfd%nvctr_c + 7*lin%lzd%llr(ilr)%wfd%nvctr_f) + &
              2*(lin%lzd%llr(jlr)%wfd%nvctr_c + 7*lin%lzd%llr(jlr)%wfd%nvctr_f)
      end if

      ! Compress the z wavefunction.
      call compress_forstandard(lin%lzd%llr(ilr)%d%n1, lin%lzd%llr(ilr)%d%n2, lin%lzd%llr(ilr)%d%n3, &
           lin%lzd%llr(ilr)%d%nfl1, lin%lzd%llr(ilr)%d%nfu1, &
           lin%lzd%llr(ilr)%d%nfl2, lin%lzd%llr(ilr)%d%nfu2, lin%lzd%llr(ilr)%d%nfl3, lin%lzd%llr(ilr)%d%nfu3, &
           lin%lzd%llr(ilr)%wfd%nseg_c, lin%lzd%llr(ilr)%wfd%nvctr_c, lin%lzd%llr(ilr)%wfd%keyg, lin%lzd%llr(ilr)%wfd%keyv, &
           lin%lzd%llr(ilr)%wfd%nseg_f, lin%lzd%llr(ilr)%wfd%nvctr_f, &
           lin%lzd%llr(ilr)%wfd%keyg(1,lin%lzd%llr(ilr)%wfd%nseg_c+min(1,lin%lzd%llr(ilr)%wfd%nseg_f)), &
           lin%lzd%llr(ilr)%wfd%keyv(lin%lzd%llr(ilr)%wfd%nseg_c+min(1,lin%lzd%llr(ilr)%wfd%nseg_f)),  &
           scal, phiz_c, phiz_f, phiLoc(istz_c), phiLoc(istz_f))
      if(iorb<lin%orbs%norbp) then
          jlr=lin%orbs%inWhichLocregp(iorb+1)
          istz_c = istz_c + lin%lzd%llr(ilr)%wfd%nvctr_c + 7*lin%lzd%llr(ilr)%wfd%nvctr_f + &
              3*(lin%lzd%llr(jlr)%wfd%nvctr_c + 7*lin%lzd%llr(jlr)%wfd%nvctr_f)
      end if

      call deallocateWorkarrays()

  end do


  if(repartition) then
      ! Communicate the orbitals to meet the partition.
      call postCommsRepartition(iproc, nproc, lin%orbs, lin%lb%comrp, size(phiLoc), phiLoc, size(phid), phid)
      call gatherDerivativeOrbitals(iproc, nproc, lin%orbs, lin%lb%comrp)

      iall=-product(shape(phiLoc))*kind(phiLoc)
      deallocate(phiLoc, stat=istat)
      call memocc(istat, iall, 'phiLoc', subname)
  end if
  


contains

  subroutine allocateWorkarrays()
  
    ! THIS IS COPIED FROM allocate_work_arrays. Works only for free boundary.
    nf=(lin%lzd%llr(ilr)%d%nfu1-lin%lzd%llr(ilr)%d%nfl1+1)*(lin%lzd%llr(ilr)%d%nfu2-lin%lzd%llr(ilr)%d%nfl2+1)* &
       (lin%lzd%llr(ilr)%d%nfu3-lin%lzd%llr(ilr)%d%nfl3+1)

    ! Allocate work arrays
    allocate(w_c(0:lin%lzd%llr(ilr)%d%n1,0:lin%lzd%llr(ilr)%d%n2,0:lin%lzd%llr(ilr)%d%n3+ndebug), stat=istat)
    call memocc(istat, w_c, 'w_c', subname)
    w_c=0.d0

    allocate(w_f(7,lin%lzd%llr(ilr)%d%nfl1:lin%lzd%llr(ilr)%d%nfu1,lin%lzd%llr(ilr)%d%nfl2:lin%lzd%llr(ilr)%d%nfu2, &
                 lin%lzd%llr(ilr)%d%nfl3:lin%lzd%llr(ilr)%d%nfu3+ndebug), stat=istat)
    call memocc(istat, w_f, 'w_f', subname)
    w_f=0.d0
  
    allocate(w_f1(nf+ndebug), stat=istat)
    call memocc(istat, w_f1, 'w_f1', subname)
    w_f1=0.d0
    
    allocate(w_f2(nf+ndebug), stat=istat)
    call memocc(istat, w_f2, 'w_f2', subname)
    w_f2=0.d0

    allocate(w_f3(nf+ndebug), stat=istat)
    call memocc(istat, w_f3, 'w_f3', subname)
    w_f3=0.d0
  
  
    allocate(phix_f(7,lin%lzd%llr(ilr)%d%nfl1:lin%lzd%llr(ilr)%d%nfu1,lin%lzd%llr(ilr)%d%nfl2:lin%lzd%llr(ilr)%d%nfu2, &
                    lin%lzd%llr(ilr)%d%nfl3:lin%lzd%llr(ilr)%d%nfu3), stat=istat)
    call memocc(istat, phix_f, 'phix_f', subname)
    phix_f=0.d0

    allocate(phix_c(0:lin%lzd%llr(ilr)%d%n1,0:lin%lzd%llr(ilr)%d%n2,0:lin%lzd%llr(ilr)%d%n3), stat=istat)
    call memocc(istat, phix_c, 'phix_c', subname)
    phix_c=0.d0

    allocate(phiy_f(7,lin%lzd%llr(ilr)%d%nfl1:lin%lzd%llr(ilr)%d%nfu1,lin%lzd%llr(ilr)%d%nfl2:lin%lzd%llr(ilr)%d%nfu2, &
                    lin%lzd%llr(ilr)%d%nfl3:lin%lzd%llr(ilr)%d%nfu3), stat=istat)
    call memocc(istat, phiy_f, 'phiy_f', subname)
    phiy_f=0.d0

    allocate(phiy_c(0:lin%lzd%llr(ilr)%d%n1,0:lin%lzd%llr(ilr)%d%n2,0:lin%lzd%llr(ilr)%d%n3), stat=istat)
    call memocc(istat, phiy_c, 'phiy_c', subname)
    phiy_c=0.d0

    allocate(phiz_f(7,lin%lzd%llr(ilr)%d%nfl1:lin%lzd%llr(ilr)%d%nfu1,lin%lzd%llr(ilr)%d%nfl2:lin%lzd%llr(ilr)%d%nfu2, &
                    lin%lzd%llr(ilr)%d%nfl3:lin%lzd%llr(ilr)%d%nfu3), stat=istat)
    call memocc(istat, phiz_f, 'phiz_f', subname)
    phiz_f=0.d0

    allocate(phiz_c(0:lin%lzd%llr(ilr)%d%n1,0:lin%lzd%llr(ilr)%d%n2,0:lin%lzd%llr(ilr)%d%n3), stat=istat)
    call memocc(istat, phiz_c, 'phiz_c', subname)
    phiz_c=0.d0
  
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

end subroutine getDerivativeBasisFunctions2






subroutine initializeRepartitionOrbitals(iproc, nproc, tag, lin)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
integer,intent(inout):: tag
type(linearParameters),intent(inout):: lin

! Local variables
integer:: jproc, jorb, kproc, korb, istat, mpidest, mpisource, istdest, istsource, ncount, jlr, klr, norbdest, jjorb, iall
integer,dimension(:,:,:),allocatable:: move
character(len=*),parameter:: subname='initializeRepartitionOrbitals'


! To which position has the orbital to be moved:
!  - move(1,i,j)=k -> orbital i on process j has to be sent to process k
!  - move(2,i,j)=l -> orbital i on process j has to be sent to position l (orbital number l)

allocate(move(2,4*maxval(lin%orbs%norb_par(:,0)),0:nproc-1), stat=istat)
call memocc(istat, move, 'move', subname)

kproc=0
korb=0
do jproc=0,nproc-1
    do jorb=1,4*lin%orbs%norb_par(jproc,0)
        korb=korb+1
        if(korb>lin%lb%orbs%norb_par(kproc,0)) then
            kproc=kproc+1
            korb=1
        end if
        move(1,jorb,jproc)=kproc
        move(2,jorb,jproc)=korb
    end do
end do


allocate(lin%lb%comrp%communComplete(4*maxval(lin%orbs%norb_par(:,0)),0:nproc-1), stat=istat)
call memocc(istat, lin%lb%comrp%communComplete, 'lin%lb%comrp%communComplete', subname)



allocate(lin%lb%comrp%comarr(8,4*maxval(lin%orbs%norb_par(:,0)),0:nproc-1), stat=istat)
call memocc(istat, lin%lb%comrp%comarr, 'lin%lb%comrp%comarr', subname)

! Determine the indices of starting and receive buffer.
do jproc=0,nproc-1
    istsource=1
    do jorb=1,4*lin%orbs%norb_par(jproc,0)
        jjorb=ceiling(dble(jorb)/4.d0)
        jlr=lin%orbs%inWhichLocreg(jjorb+lin%orbs%isorb_par(jproc))
        mpisource=jproc
        ncount=lin%lzd%llr(jlr)%wfd%nvctr_c+7*lin%lzd%llr(jlr)%wfd%nvctr_f
        mpidest=move(1,jorb,jproc)
        norbdest=move(2,jorb,jproc)
        istdest=1
        do korb=1,norbdest-1
            klr=lin%lb%orbs%inWhichLocreg(korb+lin%lb%orbs%isorb_par(mpidest))
            istdest=istdest+lin%lzd%llr(klr)%wfd%nvctr_c+7*lin%lzd%llr(klr)%wfd%nvctr_f
        end do
        tag=tag+1
        !write(*,'(7(a,i0))') 'init on iproc=',iproc,': process ',mpisource,' sends ',ncount,' elements from position ',istsource,' to position ',istdest,' on process ',mpidest,'; tag=',tag
        call setCommsParameters(mpisource, mpidest, istsource, istdest, ncount, tag, lin%lb%comrp%comarr(1,jorb,jproc))
        istsource=istsource+ncount
    end do
end do

iall=-product(shape(move))*kind(move)
deallocate(move, stat=istat)
call memocc(istat, iall, 'move', subname)


end subroutine initializeRepartitionOrbitals




subroutine postCommsRepartition(iproc, nproc, orbs, comrp, nsendBuf, sendBuf, nrecvBuf, recvBuf)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, nsendBuf, nrecvBuf
type(orbitals_data),intent(in):: orbs
type(p2pCommsRepartition),intent(inout):: comrp
real(8),dimension(nsendBuf),intent(in):: sendBuf
real(8),dimension(nrecvBuf),intent(out):: recvBuf

! Local variables
integer:: nsends, nreceives, jproc, jorb, mpisource, mpidest, istsource, istdest, ncount, tag, ierr



nsends=0
nreceives=0
comrp%communComplete=.false.
do jproc=0,nproc-1
    do jorb=1,4*orbs%norb_par(jproc,0)
        mpisource=comrp%comarr(1,jorb,jproc)
        istsource=comrp%comarr(2,jorb,jproc)
        ncount=comrp%comarr(3,jorb,jproc)
        mpidest=comrp%comarr(4,jorb,jproc)
        istdest=comrp%comarr(5,jorb,jproc)
        tag=comrp%comarr(6,jorb,jproc)
        !write(*,'(10(a,i0))') 'post on iproc=',iproc,': process ',mpisource,' sends ',ncount,' elements from position ',istsource,' to position ',istdest,' on process ',mpidest,'; tag=',tag,'  jproc=',jproc,' jorb=',jorb,' ncnt=',comrp%comarr(3,jorb,jproc)
        if(mpisource/=mpidest) then
            ! The orbitals are on different processes, so we need a point to point communication.
            if(iproc==mpisource) then
                !write(*,'(9(a,i0))') 'process ', mpisource, ' sends ', ncount, ' elements from position ', istsource, ' to position ', istdest, ' on process ', mpidest, ', tag=',tag,'  jproc=',jproc,' jorb=',jorb,' ncnt=',comrp%comarr(3,jorb,jproc)
                call mpi_isend(sendBuf(istsource), ncount, mpi_double_precision, mpidest, tag, mpi_comm_world, &
                     comrp%comarr(7,jorb,jproc), ierr)
                !!do ierr=istsource,istsource+ncount-1
                !!    write(tag,*) ierr, sendBuf(ierr)
                !!end do
                !call mpi_isend(sendBuf(istsource), ncount, mpi_double_precision, mpidest, tag, mpi_comm_world, lin%comsr%comarr(8,iorb,jproc), ierr)
                comrp%comarr(8,jorb,jproc)=mpi_request_null !is this correct?
                nsends=nsends+1
            else if(iproc==mpidest) then
                !write(*,'(9(a,i0))') 'process ', mpidest, ' receives ', ncount, ' elements at position ', istdest, ' from position ', istsource, ' on process ', mpisource, ', tag=',tag,'  jproc=',jproc,' jorb=',jorb,' ncnt=',comrp%comarr(3,jorb,jproc)
                call mpi_irecv(recvBuf(istdest), ncount, mpi_double_precision, mpisource, tag, mpi_comm_world, &
                     comrp%comarr(8,jorb,jproc), ierr)
                comrp%comarr(7,jorb,jproc)=mpi_request_null !is this correct?
                nreceives=nreceives+1
            else
                comrp%comarr(7,jorb,jproc)=mpi_request_null
                comrp%comarr(8,jorb,jproc)=mpi_request_null
            end if
        else
            ! The orbitals are on the same process, so simply copy them.
            if(iproc==mpisource) then
                call dcopy(ncount, sendBuf(istsource), 1, recvBuf(istdest), 1)
                !!do ierr=istsource,istsource+ncount-1
                !!    write(tag,*) ierr, sendBuf(ierr)
                !!end do
                !write(*,'(6(a,i0))') 'process ', iproc, ' copies ', ncount, ' elements from position ', istsource, ' to position ', istdest, ' on process ', iproc, ', tag=',tag
                comrp%comarr(7,jorb,jproc)=mpi_request_null
                comrp%comarr(8,jorb,jproc)=mpi_request_null
                nsends=nsends+1
                nreceives=nreceives+1
                comrp%communComplete(jorb,iproc)=.true.
            else
                comrp%comarr(7,jorb,jproc)=mpi_request_null
                comrp%comarr(8,jorb,jproc)=mpi_request_null
                comrp%communComplete(jorb,iproc)=.true.
            end if

        end if
    end do
end do


end subroutine postCommsRepartition




subroutine gatherDerivativeOrbitals(iproc, nproc, orbs, comrp)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data),intent(in):: orbs
type(p2pCommsRepartition),intent(inout):: comrp

! Local variables
integer:: jorb, mpisource, mpidest, nfast, nslow, nsameproc, ierr, jproc
integer,dimension(mpi_status_size):: stat
logical:: sendComplete, receiveComplete


! Check whether the communications have completed.
nfast=0
nsameproc=0
testLoop: do
    do jproc=0,nproc-1
        do jorb=1,4*orbs%norb_par(jproc,0)
            if(comrp%communComplete(jorb,jproc)) cycle
            call mpi_test(comrp%comarr(7,jorb,jproc), sendComplete, stat, ierr)
            call mpi_test(comrp%comarr(8,jorb,jproc), receiveComplete, stat, ierr)
            ! Attention: mpi_test is a local function.
            if(sendComplete .and. receiveComplete) comrp%communComplete(jorb,jproc)=.true.
            !!if(comrp%communComplete(jorb,jproc)) then
            !!    !write(*,'(2(a,i0))') 'fast communication; process ', iproc, ' has received orbital ', jorb
            !!    mpisource=comrp%comarr(1,jorb,jproc)
            !!    mpidest=comrp%comarr(4,jorb,jproc)
            !!    if(mpisource/=mpidest) then
            !!        nfast=nfast+1
            !!    else
            !!        nsameproc=nsameproc+1
            !!    end if
            !!end if
        end do
    end do
    ! If we made it until here, either all all the communication is
    ! complete or we better wait for each single orbital.
    exit testLoop
end do testLoop

! Since mpi_test is a local function, check whether the communication has completed on all processes.
call mpiallred(comrp%communComplete(1,0), nproc*4*maxval(orbs%norb_par(:,0)), mpi_land, mpi_comm_world, ierr)

! Wait for the communications that have not completed yet
nslow=0
do jproc=0,nproc-1
    do jorb=1,4*orbs%norb_par(jproc,0)
        if(comrp%communComplete(jorb,jproc)) then
            mpisource=comrp%comarr(1,jorb,jproc)
            mpidest=comrp%comarr(4,jorb,jproc)
            if(mpisource==mpidest) then
                nsameproc=nsameproc+1
            else
                nfast=nfast+1
            end if
            cycle
        end if
        !write(*,'(2(a,i0))') 'process ', iproc, ' is waiting for orbital ', korb
        nslow=nslow+1
        call mpi_wait(comrp%comarr(7,jorb,jproc), stat, ierr)   !COMMENTED BY PB
        call mpi_wait(comrp%comarr(8,jorb,jproc), stat, ierr)   !COMMENTED BY PB
        comrp%communComplete(jorb,jproc)=.true.
    end do
end do

!call mpiallred(nreceives, 1, mpi_sum, mpi_comm_world, ierr)
!call mpiallred(nfast, 1, mpi_sum, mpi_comm_world, ierr)
!call mpiallred(nslow, 1, mpi_sum, mpi_comm_world, ierr)
!call mpiallred(nsameproc, 1, mpi_sum, mpi_comm_world, ierr)
if(iproc==0) write(*,'(1x,2(a,i0),a)') 'statistics: - ', nfast+nslow, ' point to point communications, of which ', &
                       nfast, ' could be overlapped with computation.'
if(iproc==0) write(*,'(1x,a,i0,a)') '            - ', nsameproc, ' copies on the same processor.'





end subroutine gatherDerivativeOrbitals
