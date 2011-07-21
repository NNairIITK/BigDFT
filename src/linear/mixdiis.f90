subroutine mixrhopotDIIS(iproc, nproc, ndimpot, rhopot, rhopotold, mixdiis, ndimtot, alphaMix, mixMeth, pnrm)
use module_base
use module_types
use libxc_functionals
use module_interfaces, exceptThisOne => mixrhopotDIIS
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, ndimpot, ndimtot, mixMeth
real(8),dimension(ndimpot),intent(in):: rhopotold
real(8),dimension(ndimpot),intent(out):: rhopot
type(mixrhopotDIISParameters),intent(inout):: mixdiis
real(8),intent(in):: alphaMix
real(8),intent(out):: pnrm

! Local variables
integer:: ist, jst, i, j, mi, istat, lwork, info
integer:: mj, jj, k, jjst, isthist, ierr, iall
real(8):: ddot, tt
real(8),dimension(:,:),allocatable:: mat
real(8),dimension(:),allocatable:: rhs, work, rhopotres
integer,dimension(:),allocatable:: ipiv
character(len=*),parameter:: subname='optimizeDIIS'

! Allocate the local arrays.
allocate(mat(mixdiis%isx+1,mixdiis%isx+1), stat=istat)
call memocc(istat, mat, 'mat', subname)
allocate(rhs(mixdiis%isx+1), stat=istat)
call memocc(istat, rhs, 'rhs', subname)
lwork=100*mixdiis%isx
allocate(work(lwork), stat=istat)
call memocc(istat, work, 'work', subname)
allocate(ipiv(mixdiis%isx+1), stat=istat)
call memocc(istat, ipiv, 'ipiv', subname)
allocate(rhopotres(ndimpot), stat=istat)
call memocc(istat, rhopotres, 'rhopotres', subname)

! Calculate the residue.
pnrm=0.d0
!tt=1.d0-alphaMix
!do i=1,max(Glr%d%n1i*Glr%d%n2i*n3p,1)*input%nspin
do i=1,ndimpot
    !rhopotres(i) = alphaMix*rhopot(i) + tt*rhopotold(i)
    tt = rhopotold(i) - rhopot(i)
    rhopotres(i) = alphaMix*tt
    pnrm=pnrm+tt**2
end do
call mpiallred(pnrm, 1, mpi_sum, mpi_comm_world, ierr)
pnrm=sqrt(pnrm)/ndimtot


mat=0.d0
rhs=0.d0

! Copy rhopot and rhopotres to the DIIS history.
jst=(mixdiis%mis-1)*ndimpot+1
call dcopy(ndimpot, rhopotold(1), 1, mixdiis%rhopotHist(jst), 1)
call dcopy(ndimpot, rhopotres(1), 1, mixdiis%rhopotresHist(jst), 1)

! Shift the DIIS matrix left up if we reached the maximal history length.
if(mixdiis%is>mixdiis%isx) then
   do i=1,mixdiis%isx-1
      do j=1,i
         mixdiis%mat(j,i)=mixdiis%mat(j+1,i+1)
      end do
   end do
end if




! Calculate a new line for the matrix.
i=max(1,mixdiis%is-mixdiis%isx+1)
do j=i,mixdiis%is
   mi=mod(j-1,mixdiis%isx)+1
   ist=(mi-1)*ndimpot+1
   mixdiis%mat(j-i+1,min(mixdiis%isx,mixdiis%is))=ddot(ndimpot, rhopotres(1), 1, mixdiis%rhopotresHist(ist), 1)
end do

! Sum up over all processes.
call mpiallred(mixdiis%mat(1,min(mixdiis%isx,mixdiis%is)), min(mixdiis%is,mixdiis%isx), mpi_sum, mpi_comm_world, ierr)



! Copy the matrix to an auxiliary array and fill with the zeros and ones.
do i=1,min(mixdiis%isx,mixdiis%is)
    mat(i,min(mixdiis%isx,mixdiis%is)+1)=1.d0
    rhs(i)=0.d0
    do j=i,min(mixdiis%isx,mixdiis%is)
        mat(i,j)=mixdiis%mat(i,j)
    end do
end do
mat(min(mixdiis%isx,mixdiis%is)+1,min(mixdiis%isx,mixdiis%is)+1)=0.d0
rhs(min(mixdiis%isx,mixdiis%is)+1)=1.d0


! Solve the linear system
if(mixdiis%is>1) then
   call dsysv('u', min(mixdiis%isx,mixdiis%is)+1, 1, mat, mixdiis%isx+1,  & 
        ipiv, rhs(1), mixdiis%isx+1, work, lwork, info)
   if (info /= 0) then
      write(*,'(a,i0)') 'ERROR in dsysv (subroutine mixrhopotDIIS), info=', info
      stop
   end if
else
   rhs(1)=1.d0
endif

!!do i=1,ndimpot
!!    write(100+iproc,*) i, rhopot(i)
!!    write(110+iproc,*) i, rhopotres(i)
!!end do


! Make a new guess for the density/potential.
! If we are mixing the density (mixMeth==1) it is initialized to 0 or 10^-20, depending on the functional.
! If we are mixing the potential (mixMeth==2) it is always initialized to 0.
if (libxc_functionals_isgga() .or. mixMeth==2) then
    call razero(ndimpot, rhopot)
else
    ! There is no mpi_allreduce, therefore directly initialize to
    ! 10^-20 and not 10^-20/nproc.
    rhopot=1.d-20
    !call tenminustwenty(nrho, rho, nproc)
end if

isthist=max(1,mixdiis%is-mixdiis%isx+1)
jj=0
do j=isthist,mixdiis%is
    jj=jj+1
    mj=mod(j-1,mixdiis%isx)+1
    jjst=(mj-1)*ndimpot
    do k=1,ndimpot
        rhopot(k) = rhopot(k) + rhs(jj)*(mixdiis%rhopotHist(jjst+k)-mixdiis%rhopotresHist(jjst+k))
    end do
end do
!!do i=1,ndimpot
!!    write(120+iproc,*) i, rhopot(i)
!!end do


iall=-product(shape(mat))*kind(mat)
deallocate(mat, stat=istat)
call memocc(istat, iall, 'mat', subname)

iall=-product(shape(rhs))*kind(rhs)
deallocate(rhs, stat=istat)
call memocc(istat, iall, 'rhs', subname)

iall=-product(shape(work))*kind(work)
deallocate(work, stat=istat)
call memocc(istat, iall, 'work', subname)

iall=-product(shape(ipiv))*kind(ipiv)
deallocate(ipiv, stat=istat)
call memocc(istat, iall, 'ipiv', subname)

iall=-product(shape(rhopotres))*kind(rhopotres)
deallocate(rhopotres, stat=istat)
call memocc(istat, iall, 'rhopotres', subname)


end subroutine mixrhopotDIIS




subroutine initializeMixrhopotDIIS(isx, ndimpot, mixdiis)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: isx, ndimpot
type(mixrhopotDIISParameters),intent(out):: mixdiis

! Local variables
integer:: istat
character(len=*),parameter:: subname='initializeMixrhopotDIIS'


mixdiis%isx=isx
mixdiis%is=0
allocate(mixdiis%mat(mixdiis%isx,mixdiis%isx), stat=istat)
call memocc(istat, mixdiis%mat, 'mixdiis%mat', subname)
allocate(mixdiis%rhopotHist(mixdiis%isx*ndimpot), stat=istat)
call memocc(istat, mixdiis%rhopotHist, 'mixdiis%rhopotHist', subname)
allocate(mixdiis%rhopotresHist(mixdiis%isx*ndimpot), stat=istat)
call memocc(istat, mixdiis%rhopotresHist, 'mixdiis%rhopotresHist', subname)


end subroutine initializeMixrhopotDIIS



subroutine deallocateMixrhopotDIIS(mixdiis)
use module_base
use module_types
implicit none

! Calling arguments
type(mixrhopotDIISParameters),intent(inout):: mixdiis

! Local variables
integer:: istat, iall
character(len=*),parameter:: subname='deallocateMixrhopotDIIS'


iall=-product(shape(mixdiis%mat))*kind(mixdiis%mat)
deallocate(mixdiis%mat, stat=istat)
call memocc(istat, iall, 'mixdiis%mat', subname)

iall=-product(shape(mixdiis%rhopotHist))*kind(mixdiis%rhopotHist)
deallocate(mixdiis%rhopotHist, stat=istat)
call memocc(istat, iall, 'mixdiis%rhopotHist', subname)

iall=-product(shape(mixdiis%rhopotresHist))*kind(mixdiis%rhopotresHist)
deallocate(mixdiis%rhopotresHist, stat=istat)
call memocc(istat, iall, 'mixdiis%rhopotresHist', subname)

end subroutine deallocateMixrhopotDIIS
