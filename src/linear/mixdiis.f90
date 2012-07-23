subroutine mixPotential(iproc, n3p, Glr, input, alphaMix, rhopotOld, rhopot, pnrm)
!
! Purpose:
! ========
!   Mixes the potential in order to get a self consistent potential.
!
! Calling arguments:
! ==================
!   Input arguments
!   ---------------
!
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, n3p
type(locreg_descriptors),intent(in) :: Glr
type(input_variables),intent(in):: input
real(8),intent(in):: alphaMix
real(dp),dimension(max(Glr%d%n1i*Glr%d%n2i*n3p,1)*input%nspin),intent(in):: rhopotOld
real(dp),dimension(max(Glr%d%n1i*Glr%d%n2i*n3p,1)*input%nspin),intent(in out):: rhopot
real(8),intent(out):: pnrm

! Local variables
integer:: i, ierr
real(8):: tt

  call timing(iproc,'mix_linear    ','ON')

  pnrm=0.d0
  tt=1.d0-alphaMix
  !do i=1,max(Glr%d%n1i*Glr%d%n2i*n3p,1)*input%nspin
  !do i=1,max(Glr%d%n1i*Glr%d%n2i*n3p,1)
  do i=1,Glr%d%n1i*Glr%d%n2i*n3p
      pnrm=pnrm+(rhopot(i)-rhopotOld(i))**2
      rhopot(i)=tt*rhopotOld(i)+alphaMix*rhopot(i)
  end do
  call mpiallred(pnrm, 1, mpi_sum, mpi_comm_world, ierr)
  pnrm=sqrt(pnrm)/(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i*input%nspin)
  pnrm=pnrm/alphaMix

  if(pnrm<input%lin%convCritMix) then
      if(iproc==0) write(*,*) 'keep old density / potential'
      call dcopy(Glr%d%n1i*Glr%d%n2i*n3p, rhopotOld(1), 1, rhopot(1), 1)
  end if

  call timing(iproc,'mix_linear    ','OF')

end subroutine mixPotential



subroutine mix_main(iproc, nproc, mixHist, compare_outer_loop, input, glr, alpha_mix, &
           denspot, mixdiis, rhopotold, rhopotold_out, pnrm, pnrm_out)
  use module_base
  use module_types
  use module_interfaces, except_this_one => mix_main
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc, mixHist
  logical,intent(in):: compare_outer_loop
  type(input_variables),intent(in):: input
  type(locreg_descriptors),intent(in):: glr
  real(8),intent(in):: alpha_mix
  type(DFT_local_fields),intent(inout):: denspot
  type(mixrhopotDIISParameters),intent(inout):: mixdiis
  real(8),dimension(max(glr%d%n1i*glr%d%n2i*denspot%dpbox%n3p,1)*input%nspin),intent(inout):: rhopotold, rhopotold_out
  real(8),intent(out):: pnrm, pnrm_out
  
  ! Local variables
  integer:: ndimtot, i, ierr

  ! Mix the density.
  if(mixHist==0) then
      call mixPotential(iproc, denspot%dpbox%n3p, glr, input, alpha_mix, rhopotOld, denspot%rhov, pnrm)
  else 
      ndimtot=glr%d%n1i*glr%d%n2i*glr%d%n3i
      mixdiis%mis=mod(mixdiis%is,mixdiis%isx)+1
      mixdiis%is=mixdiis%is+1
      call mixrhopotDIIS(iproc, nproc, denspot%dpbox%ndimpot,&
           denspot%rhov, rhopotold, mixdiis, ndimtot, alpha_mix, 1, pnrm)
  end if
  ! Determine the change in the density between this iteration and the last iteration in the outer loop.
  if(compare_outer_loop) then
      pnrm_out=0.d0
      do i=1,glr%d%n1i*glr%d%n2i*denspot%dpbox%n3p
          pnrm_out=pnrm_out+(denspot%rhov(i)-rhopotOld_out(i))**2
      end do
      call mpiallred(pnrm_out, 1, mpi_sum, mpi_comm_world, ierr)
      pnrm_out=sqrt(pnrm_out)/(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i*input%nspin)
      ! Do not divide by alpha_mix here since it is the difference in the outer loop.
      call dcopy(max(Glr%d%n1i*Glr%d%n2i*denspot%dpbox%n3p,1)*input%nspin, denspot%rhov(1), 1, rhopotOld_out(1), 1)
  end if

  ! Copy the current charge density.
  call dcopy(max(Glr%d%n1i*Glr%d%n2i*denspot%dpbox%n3p,1)*input%nspin, denspot%rhov(1), 1, rhopotOld(1), 1)

end subroutine mix_main


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

call timing(iproc,'mix_DIIS      ','ON')

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


!!mat=0.d0
call to_zero((mixdiis%isx+1)**2, mat(1,1))
!!rhs=0.d0
call to_zero(mixdiis%isx+1, rhs(1))

! Copy rhopot and rhopotres to the DIIS history.
jst=(mixdiis%mis-1)*ndimpot+1
!call dcopy(ndimpot, rhopotold, 1, mixdiis%rhopotHist, 1)
!call dcopy(ndimpot, rhopotres, 1, mixdiis%rhopotresHist, 1)
call dcopy(ndimpot, rhopotold, 1, mixdiis%rhopotHist(jst), 1)
call dcopy(ndimpot, rhopotres, 1, mixdiis%rhopotresHist(jst), 1)

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
   if(ndimpot>0) then
       mixdiis%mat(j-i+1,min(mixdiis%isx,mixdiis%is))=ddot(ndimpot, rhopotres(1), 1, mixdiis%rhopotresHist(ist), 1)
   else
       mixdiis%mat(j-i+1,min(mixdiis%isx,mixdiis%is))=0.d0
   end if
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

call timing(iproc,'mix_DIIS      ','OF')

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
