subroutine optimizeDIIS(iproc, nproc, orbs, lorbs, lzd, hphi, phi, ldiis, it)
use module_base
use module_types
use module_interfaces, exceptThisOne => optimizeDIIS
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, it
type(orbitals_data),intent(in):: orbs, lorbs
type(local_zone_descriptors),intent(in):: lzd
real(8),dimension(max(lorbs%npsidim_orbs,lorbs%npsidim_comp)),intent(in):: hphi
real(8),dimension(max(lorbs%npsidim_orbs,lorbs%npsidim_comp)),intent(inout):: phi
type(localizedDIISParameters),intent(inout):: ldiis

! Local variables
integer:: iorb, jorb, ist, ilr, ncount, jst, i, j, mi, ist1, ist2, jlr, istat, lwork, info
integer:: mj, jj, k, jjst, isthist, ierr, iall
real(8):: ddot
real(8),dimension(:,:),allocatable:: mat
real(8),dimension(:),allocatable:: rhs, work
integer,dimension(:),allocatable:: ipiv
character(len=*),parameter:: subname='optimizeDIIS'

call timing(iproc,'optimize_DIIS ','ON')

! Allocate the local arrays.
allocate(mat(ldiis%isx+1,ldiis%isx+1), stat=istat)
call memocc(istat, mat, 'mat', subname)
allocate(rhs(ldiis%isx+1), stat=istat)
call memocc(istat, rhs, 'rhs', subname)
!lwork=100*ldiis%isx
!allocate(work(lwork), stat=istat)
!call memocc(istat, work, 'work', subname)
allocate(ipiv(ldiis%isx+1), stat=istat)
call memocc(istat, ipiv, 'ipiv', subname)

!!mat=0.d0
!!rhs=0.d0
call to_zero((ldiis%isx+1)**2, mat(1,1))
call to_zero(ldiis%isx+1, rhs(1))

! Copy phi and hphi to history.
ist=1
do iorb=1,orbs%norbp
    jst=1
    do jorb=1,iorb-1
        !jlr=onWhichAtom(jorb)
        jlr=orbs%inwhichlocreg(orbs%isorb+jorb)
        ncount=lzd%llr(jlr)%wfd%nvctr_c+7*lzd%llr(jlr)%wfd%nvctr_f
        jst=jst+ncount*ldiis%isx
    end do
    !ilr=onWhichAtom(iorb)
    ilr=orbs%inwhichlocreg(orbs%isorb+iorb)
    ncount=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
    jst=jst+(ldiis%mis-1)*ncount
    call dcopy(ncount, phi(ist), 1, ldiis%phiHist(jst), 1)
    call dcopy(ncount, hphi(ist), 1, ldiis%hphiHist(jst), 1)


    !ilr=onWhichAtom(iorb)
    ilr=orbs%inwhichlocreg(orbs%isorb+iorb)
    ncount=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
    ist=ist+ncount
end do

do iorb=1,orbs%norbp
    ! Shift the DIIS matrix left up if we reached the maximal history length.
    if(ldiis%is>ldiis%isx) then
       do i=1,ldiis%isx-1
          do j=1,i
             ldiis%mat(j,i,iorb)=ldiis%mat(j+1,i+1,iorb)
             !!write(3100+iproc,*) ldiis%mat(j,i,iorb)
          end do
       end do
    end if
end do



do iorb=1,orbs%norbp

    ! Calculate a new line for the matrix.
    i=max(1,ldiis%is-ldiis%isx+1)
    jst=1
    ist1=1
    do jorb=1,iorb-1
        !jlr=onWhichAtom(jorb)
        jlr=orbs%inwhichlocreg(orbs%isorb+jorb)
        ncount=lzd%llr(jlr)%wfd%nvctr_c+7*lzd%llr(jlr)%wfd%nvctr_f
        jst=jst+ncount*ldiis%isx
        ist1=ist1+ncount
    end do
    !ilr=onWhichAtom(iorb)
    ilr=orbs%inwhichlocreg(orbs%isorb+iorb)
    ncount=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
    do j=i,ldiis%is
       mi=mod(j-1,ldiis%isx)+1
       ist2=jst+(mi-1)*ncount
       if(ist2>size(ldiis%hphiHist)) then
           write(*,'(a,7i8)') 'ERROR ist2: iproc, iorb, ldiis%is, mi, ncount, ist2, size(ldiis%hphiHist)', iproc, iorb, ldiis%is,&
                               mi, ncount, ist2, size(ldiis%hphiHist)
       end if
       ldiis%mat(j-i+1,min(ldiis%isx,ldiis%is),iorb)=ddot(ncount, hphi(ist1), 1, ldiis%hphiHist(ist2), 1)
       !!write(3000+iproc,*) ldiis%mat(j-i+1,min(ldiis%isx,ldiis%is),iorb)
       !!write(3200+iproc,'(4i8,2es20.12)') mi, ldiis%is, ist1, ist2, hphi(ist1), ldiis%hphiHist(ist2)
       ist2=ist2+ncount
    end do
end do


ist=1
do iorb=1,orbs%norbp
    
    ! Copy the matrix to an auxiliary array and fill with the zeros and ones.
    do i=1,min(ldiis%isx,ldiis%is)
        mat(i,min(ldiis%isx,ldiis%is)+1)=1.d0
        rhs(i)=0.d0
        do j=i,min(ldiis%isx,ldiis%is)
            mat(i,j)=ldiis%mat(i,j,iorb)
            !if(iproc==0) write(*,'(a,2i8,es14.3)') 'i, j, mat(i,j)', i, j, mat(i,j)
            !!write(*,'(a,3i8,es14.3)') 'proc, i, j, mat(i,j)', iproc, i, j, mat(i,j)
        end do
    end do
    mat(min(ldiis%isx,ldiis%is)+1,min(ldiis%isx,ldiis%is)+1)=0.d0
    rhs(min(ldiis%isx,ldiis%is)+1)=1.d0

    !make the matrix symmetric (hermitian) to use DGESV (ZGESV) (no work array, more stable)
    do i=1,min(ldiis%isx,ldiis%is)+1
       do j=1,min(ldiis%isx,ldiis%is)+1
          mat(j,i) = mat(i,j)
       end do
    end do
    ! solve linear system, supposing it is general. More stable, no need of work array
    if(ldiis%is>1) then
     call dgesv(min(ldiis%isx,ldiis%is)+1,1,mat(1,1),ldiis%isx+1,  & 
                   ipiv(1),rhs(1),ldiis%isx+1,info)
     if (info /= 0) then
        write(*,'(a,i0)') 'ERROR in dgesv (subroutine optimizeDIIS), info=', info
        stop
     end if
    else
       rhs(1)=1.d0
    endif


    ! Solve the linear system
    !!if(ldiis%is>1) then
    !!   call dsysv('u', min(ldiis%isx,ldiis%is)+1, 1, mat, ldiis%isx+1,  & 
    !!        ipiv, rhs(1), ldiis%isx+1, work, lwork, info)
    !!   
    !!   if (info /= 0) then
    !!      write(*,'(a,i0)') 'ERROR in dsysv (subroutine optimizeDIIS), info=', info
    !!      stop
    !!   end if
    !!else
    !!   rhs(1)=1.d0
    !!endif

    ! Make a new guess for the orbital.
    !ilr=onWhichAtom(iorb)
    ilr=orbs%inwhichlocreg(orbs%isorb+iorb)
    ncount=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
    call razero(ncount, phi(ist))
    isthist=max(1,ldiis%is-ldiis%isx+1)
    jj=0
    jst=0
    do jorb=1,iorb-1
        !jlr=onWhichAtom(jorb)
        jlr=orbs%inwhichlocreg(orbs%isorb+jorb)
        ncount=lzd%llr(jlr)%wfd%nvctr_c+7*lzd%llr(jlr)%wfd%nvctr_f
        jst=jst+ncount*ldiis%isx
    end do
    !!write(2000+iproc,'(a,i8,100es9.2)') 'iproc, rhs',iproc, rhs(1:(ldiis%is-isthist+1))
    do j=isthist,ldiis%is
        jj=jj+1
        mj=mod(j-1,ldiis%isx)+1
        !ilr=onWhichAtom(iorb)
        ilr=orbs%inwhichlocreg(orbs%isorb+iorb)
        ncount=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
        jjst=jst+(mj-1)*ncount
        do k=1,ncount
            phi(ist+k-1) = phi(ist+k-1) + rhs(jj)*(ldiis%phiHist(jjst+k)-ldiis%hphiHist(jjst+k))
        end do
    end do


    !ilr=onWhichAtom(iorb)
    ilr=orbs%inwhichlocreg(orbs%isorb+iorb)
    ncount=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
    ist=ist+ncount
end do


iall=-product(shape(mat))*kind(mat)
deallocate(mat, stat=istat)
call memocc(istat, iall, 'mat', subname)

iall=-product(shape(rhs))*kind(rhs)
deallocate(rhs, stat=istat)
call memocc(istat, iall, 'rhs', subname)

!iall=-product(shape(work))*kind(work)
!deallocate(work, stat=istat)
!call memocc(istat, iall, 'work', subname)

iall=-product(shape(ipiv))*kind(ipiv)
deallocate(ipiv, stat=istat)
call memocc(istat, iall, 'ipiv', subname)

call timing(iproc,'optimize_DIIS ','OF')


end subroutine optimizeDIIS




subroutine initializeDIIS(isx, lzd, orbs, norb, ldiis)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: isx, norb
type(local_zone_descriptors),intent(in):: lzd
type(orbitals_data),intent(in):: orbs
type(localizedDIISParameters),intent(out):: ldiis

! Local variables
integer:: iorb, ii, istat, ilr
character(len=*),parameter:: subname='initializeDIIS'


ldiis%isx=isx
ldiis%is=0
ldiis%switchSD=.false.
ldiis%trmin=1.d100
ldiis%trold=1.d100
allocate(ldiis%mat(ldiis%isx,ldiis%isx,orbs%norbp), stat=istat)
call memocc(istat, ldiis%mat, 'ldiis%mat', subname)
ii=0
do iorb=1,orbs%norbp
    !ilr=onWhichAtom(iorb)
    ilr=orbs%inwhichlocreg(orbs%isorb+iorb)
    ii=ii+ldiis%isx*(lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f)
end do
allocate(ldiis%phiHist(ii), stat=istat)
call memocc(istat, ldiis%phiHist, 'ldiis%phiHist', subname)
allocate(ldiis%hphiHist(ii), stat=istat)
call memocc(istat, ldiis%hphiHist, 'ldiis%hphiHist', subname)
ldiis%phiHist=555.d0
ldiis%hphiHist=555.d0


end subroutine initializeDIIS



subroutine deallocateDIIS(ldiis)
use module_base
use module_types
implicit none

! Calling arguments
type(localizedDIISParameters),intent(inout):: ldiis

! Local variables
integer:: istat, iall
character(len=*),parameter:: subname='deallocateDIIS'


iall=-product(shape(ldiis%mat))*kind(ldiis%mat)
deallocate(ldiis%mat, stat=istat)
call memocc(istat, iall, 'ldiis%mat', subname)

iall=-product(shape(ldiis%phiHist))*kind(ldiis%phiHist)
deallocate(ldiis%phiHist, stat=istat)
call memocc(istat, iall, 'ldiis%phiHist', subname)

iall=-product(shape(ldiis%hphiHist))*kind(ldiis%hphiHist)
deallocate(ldiis%hphiHist, stat=istat)
call memocc(istat, iall, 'ldiis%hphiHist', subname)

end subroutine deallocateDIIS





!!!subroutine initializeDIIS_inguess(isx, norbp, matmin, onWhichAtomp, ldiis)
!!!use module_base
!!!use module_types
!!!implicit none
!!!
!!!! Calling arguments
!!!integer,intent(in):: isx, norbp
!!!type(matrixMinimization),intent(in):: matmin
!!!integer,dimension(norbp):: onWhichAtomp
!!!type(localizedDIISParameters),intent(out):: ldiis
!!!
!!!! Local variables
!!!integer:: iorb, ii, istat, ilr
!!!character(len=*),parameter:: subname='initializeDIIS_inguess'
!!!
!!!
!!!ldiis%isx=isx
!!!ldiis%is=0
!!!ldiis%switchSD=.false.
!!!ldiis%trmin=1.d100
!!!ldiis%trold=1.d100
!!!allocate(ldiis%mat(ldiis%isx,ldiis%isx,norbp), stat=istat)
!!!call memocc(istat, ldiis%mat, 'ldiis%mat', subname)
!!!ii=0
!!!do iorb=1,norbp
!!!    ilr=onWhichAtomp(iorb)
!!!    !ilr=orbs%inwhichlocreg(orbs%isorb+iorb)
!!!    !ii=ii+ldiis%isx*(lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f)
!!!    ii=ii+ldiis%isx*matmin%mlr(ilr)%norbinlr
!!!end do
!!!allocate(ldiis%phiHist(ii), stat=istat)
!!!call memocc(istat, ldiis%phiHist, 'ldiis%phiHist', subname)
!!!allocate(ldiis%hphiHist(ii), stat=istat)
!!!call memocc(istat, ldiis%hphiHist, 'ldiis%hphiHist', subname)
!!!
!!!!! Initialize the DIIS parameters 
!!!!icountSDSatur=0
!!!!icountSwitch=0
!!!!icountDIISFailureTot=0
!!!!icountDIISFailureCons=0
!!!!
!!!!! Assign the step size for SD iterations.
!!!!alpha=alphaSDx
!!!!alphaDIIS=alphaDIISx
!!!
!!!
!!!end subroutine initializeDIIS_inguess



!!subroutine optimizeDIIS_inguess(iproc, nproc, norbp, onWhichAtomp, matmin, lgrad, lcoeff, ldiis)
!!use module_base
!!use module_types
!!use module_interfaces, exceptThisOne => optimizeDIIS
!!implicit none
!!
!!! Calling arguments
!!integer,intent(in):: iproc, nproc, norbp
!!integer,dimension(norbp),intent(in):: onWhichAtomp
!!type(matrixMinimization),intent(in):: matmin
!!real(8),dimension(matmin%norbmax,norbp),intent(in):: lgrad
!!real(8),dimension(matmin%norbmax,norbp),intent(inout):: lcoeff
!!type(localizedDIISParameters),intent(inout):: ldiis
!!
!!! Local variables
!!integer:: iorb, jorb, ist, ilr, ncount, jst, i, j, mi, ist1, ist2, jlr, istat, lwork, info
!!integer:: mj, jj, jst2, k, jjst, isthist, ierr, iall
!!real(8):: ddot
!!real(8),dimension(:,:),allocatable:: mat
!!real(8),dimension(:),allocatable:: rhs, work
!!integer,dimension(:),allocatable:: ipiv
!!character(len=*),parameter:: subname='optimizeDIIS'
!!
!!! Allocate the local arrays.
!!allocate(mat(ldiis%isx+1,ldiis%isx+1), stat=istat)
!!call memocc(istat, mat, 'mat', subname)
!!allocate(rhs(ldiis%isx+1), stat=istat)
!!call memocc(istat, rhs, 'rhs', subname)
!!lwork=100*ldiis%isx
!!allocate(work(lwork), stat=istat)
!!call memocc(istat, work, 'work', subname)
!!allocate(ipiv(ldiis%isx+1), stat=istat)
!!call memocc(istat, ipiv, 'ipiv', subname)
!!
!!mat=0.d0
!!rhs=0.d0
!!
!!! Copy phi and hphi to history.
!!!ist=1
!!do iorb=1,norbp
!!    jst=1
!!    do jorb=1,iorb-1
!!        jlr=onWhichAtomp(jorb)
!!        !jlr=orbs%inwhichlocreg(orbs%isorb+jorb)
!!        ncount=matmin%mlr(jlr)%norbinlr
!!        jst=jst+ncount*ldiis%isx
!!    end do
!!    ilr=onWhichAtomp(iorb)
!!    !ilr=orbs%inwhichlocreg(orbs%isorb+iorb)
!!    ncount=matmin%mlr(ilr)%norbinlr
!!    jst=jst+(ldiis%mis-1)*ncount
!!    call dcopy(ncount, lcoeff(1,iorb), 1, ldiis%phiHist(jst), 1)
!!    call dcopy(ncount, lgrad(1,iorb), 1, ldiis%hphiHist(jst), 1)
!!
!!
!!    !ilr=onWhichAtom(iorb)
!!    !ilr=orbs%inwhichlocreg(orbs%isorb+iorb)
!!    !ncount=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
!!    !ist=ist+ncount
!!end do
!!
!!do iorb=1,norbp
!!    ! Shift the DIIS matrix left up if we reached the maximal history length.
!!    if(ldiis%is>ldiis%isx) then
!!       do i=1,ldiis%isx-1
!!          do j=1,i
!!             ldiis%mat(j,i,iorb)=ldiis%mat(j+1,i+1,iorb)
!!          end do
!!       end do
!!    end if
!!end do
!!
!!
!!
!!do iorb=1,norbp
!!
!!    ! Calculate a new line for the matrix.
!!    i=max(1,ldiis%is-ldiis%isx+1)
!!    jst=1
!!    ist1=1
!!    do jorb=1,iorb-1
!!        jlr=onWhichAtomp(jorb)
!!        !jlr=orbs%inwhichlocreg(orbs%isorb+jorb)
!!        ncount=matmin%mlr(jlr)%norbinlr
!!        jst=jst+ncount*ldiis%isx
!!        ist1=ist1+ncount
!!    end do
!!    ilr=onWhichAtomp(iorb)
!!    !ilr=orbs%inwhichlocreg(orbs%isorb+iorb)
!!    ncount=matmin%mlr(ilr)%norbinlr
!!    do j=i,ldiis%is
!!       mi=mod(j-1,ldiis%isx)+1
!!       ist2=jst+(mi-1)*ncount
!!       if(ist2>size(ldiis%hphiHist)) then
!!           write(*,'(a,7i8)') 'ERROR ist2: iproc, iorb, ldiis%is, mi, ncount, ist2, size(ldiis%hphiHist)', iproc, iorb, ldiis%is,&
!!                               mi, ncount, ist2, size(ldiis%hphiHist)
!!       end if
!!       ldiis%mat(j-i+1,min(ldiis%isx,ldiis%is),iorb)=ddot(ncount, lgrad(1,iorb), 1, ldiis%hphiHist(ist2), 1)
!!       ist2=ist2+ncount
!!    end do
!!end do
!!
!!
!!ist=1
!!do iorb=1,norbp
!!    
!!    ! Copy the matrix to an auxiliary array and fill with the zeros and ones.
!!    do i=1,min(ldiis%isx,ldiis%is)
!!        mat(i,min(ldiis%isx,ldiis%is)+1)=1.d0
!!        rhs(i)=0.d0
!!        do j=i,min(ldiis%isx,ldiis%is)
!!            mat(i,j)=ldiis%mat(i,j,iorb)
!!        end do
!!    end do
!!    mat(min(ldiis%isx,ldiis%is)+1,min(ldiis%isx,ldiis%is)+1)=0.d0
!!    rhs(min(ldiis%isx,ldiis%is)+1)=1.d0
!!
!!
!!    ! Solve the linear system
!!    if(ldiis%is>1) then
!!       call dsysv('u', min(ldiis%isx,ldiis%is)+1, 1, mat, ldiis%isx+1,  & 
!!            ipiv, rhs(1), ldiis%isx+1, work, lwork, info)
!!       
!!       if (info /= 0) then
!!          write(*,'(a,i0)') 'ERROR in dsysv (subroutine optimizeDIIS), info=', info
!!          stop
!!       end if
!!    else
!!       rhs(1)=1.d0
!!    endif
!!
!!
!!    ! Make a new guess for the orbital.
!!    ilr=onWhichAtomp(iorb)
!!    !ilr=orbs%inwhichlocreg(orbs%isorb+iorb)
!!    ncount=matmin%mlr(ilr)%norbinlr
!!    call razero(ncount, lcoeff(1,iorb))
!!    isthist=max(1,ldiis%is-ldiis%isx+1)
!!    jj=0
!!    jst=0
!!    do jorb=1,iorb-1
!!        jlr=onWhichAtomp(jorb)
!!        !jlr=orbs%inwhichlocreg(orbs%isorb+jorb)
!!        ncount=matmin%mlr(jlr)%norbinlr
!!        jst=jst+ncount*ldiis%isx
!!    end do
!!    do j=isthist,ldiis%is
!!        jj=jj+1
!!        mj=mod(j-1,ldiis%isx)+1
!!        ilr=onWhichAtomp(iorb)
!!        !ilr=orbs%inwhichlocreg(orbs%isorb+iorb)
!!        ncount=matmin%mlr(ilr)%norbinlr
!!        jjst=jst+(mj-1)*ncount
!!        do k=1,ncount
!!            !phi(ist+k-1) = phi(ist+k-1) + rhs(jj)*(ldiis%phiHist(jjst+k)-ldiis%hphiHist(jjst+k))
!!            lcoeff(k,iorb) = lcoeff(k,iorb) + rhs(jj)*(ldiis%phiHist(jjst+k)-ldiis%hphiHist(jjst+k))
!!        end do
!!    end do
!!
!!
!!    ilr=onWhichAtomp(iorb)
!!    !ilr=orbs%inwhichlocreg(orbs%isorb+iorb)
!!    ncount=matmin%mlr(ilr)%norbinlr
!!    ist=ist+ncount
!!end do
!!
!!
!!iall=-product(shape(mat))*kind(mat)
!!deallocate(mat, stat=istat)
!!call memocc(istat, iall, 'mat', subname)
!!
!!iall=-product(shape(rhs))*kind(rhs)
!!deallocate(rhs, stat=istat)
!!call memocc(istat, iall, 'rhs', subname)
!!
!!iall=-product(shape(work))*kind(work)
!!deallocate(work, stat=istat)
!!call memocc(istat, iall, 'work', subname)
!!
!!iall=-product(shape(ipiv))*kind(ipiv)
!!deallocate(ipiv, stat=istat)
!!call memocc(istat, iall, 'ipiv', subname)
!!
!!
!!
!!end subroutine optimizeDIIS_inguess
