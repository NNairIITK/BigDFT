!> @file
!! Mixes the potential in order to get a self consistent potential.
!! @author
!!    Copyright (C) 2013-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


subroutine mixPotential(iproc, nproc, n3d, n3p, Glr, input, alphaMix, rhopotOld, ioffset, rhopot, pnrm)

  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  integer,intent(in):: n3d, n3p, ioffset
  type(locreg_descriptors),intent(in) :: Glr
  type(input_variables),intent(in) :: input
  real(8),intent(in) :: alphaMix
  real(dp),dimension(max(Glr%d%n1i*Glr%d%n2i*n3d,1)*input%nspin),intent(in) :: rhopotOld
  real(dp),dimension(max(Glr%d%n1i*Glr%d%n2i*n3d,1)*input%nspin),intent(in out) :: rhopot
  real(8),intent(out) :: pnrm

  ! Local variables
  integer :: i, ierr
  real(8) :: tt

  call timing(iproc,'mix_linear    ','ON')

  pnrm=0.d0
  tt=1.d0-alphaMix
  !do i=1,max(Glr%d%n1i*Glr%d%n2i*n3p,1)*input%nspin
  !do i=1,max(Glr%d%n1i*Glr%d%n2i*n3p,1)
  do i=1,Glr%d%n1i*Glr%d%n2i*n3p
      pnrm=pnrm+(rhopot(ioffset+i)-rhopotOld(ioffset+i))**2
  end do
  do i=1,Glr%d%n1i*Glr%d%n2i*n3d
      rhopot(i)=tt*rhopotOld(i)+alphaMix*rhopot(i)
  end do

  if (nproc > 1) then
    call mpiallred(pnrm, 1, mpi_sum, bigdft_mpi%mpi_comm)
  end if

  pnrm=sqrt(pnrm)/(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i*input%nspin)
  pnrm=pnrm/alphaMix

  !!if(pnrm<input%lin%convCritMix) then
  !!    if(iproc==0) write(*,*) 'keep old density / potential'
  !!    call vcopy(Glr%d%n1i*Glr%d%n2i*n3p, rhopotOld(1), 1, rhopot(1), 1)
  !!end if

  call timing(iproc,'mix_linear    ','OF')

end subroutine mixPotential


subroutine mix_main(iproc, nproc, mix_mode, mixHist, input, glr, alpha_mix, &
           denspot, mixdiis, rhopotold, pnrm)
  use module_base
  use module_types
  use module_interfaces, except_this_one => mix_main
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc, mix_mode, mixHist
  type(input_variables),intent(in):: input
  type(locreg_descriptors),intent(in):: glr
  real(8),intent(in):: alpha_mix
  type(DFT_local_fields),intent(inout):: denspot
  type(mixrhopotDIISParameters),intent(inout):: mixdiis
  real(8),dimension(max(glr%d%n1i*glr%d%n2i*denspot%dpbox%n3d,1)*input%nspin),intent(inout):: rhopotold
  real(8),intent(out):: pnrm
  
  ! Local variables
  integer:: ndimtot, ioffset
  
  ! Offset to calculate the change in the density / potential. Since the density
  ! contains some buffer in the case of GGA, a non-zero shift is required. For the
  ! potential, however, no buffers are present, so no shift is needed.
  if (mix_mode==LINEAR_MIXPOT_SIMPLE) then
      ioffset=0
  else 
      ioffset=glr%d%n1i*glr%d%n2i*denspot%dpbox%i3xcsh
  end if


  ! Mix the density.
  if(mixHist==0) then
      call mixPotential(iproc, nproc, denspot%dpbox%n3d, denspot%dpbox%n3p, &
           glr, input, alpha_mix, rhopotOld, ioffset, denspot%rhov, pnrm)
  else 
      ndimtot=glr%d%n1i*glr%d%n2i*glr%d%n3i
      mixdiis%mis=mod(mixdiis%is,mixdiis%isx)+1
      mixdiis%is=mixdiis%is+1
      call mixrhopotDIIS(iproc, nproc, denspot%dpbox%n3d, denspot%dpbox%n3p, glr, input, &
           denspot%rhov, rhopotold, mixdiis, alpha_mix, ioffset, 1, pnrm, denspot%xc)
  end if

  ! Copy the current charge density.
  call vcopy(max(Glr%d%n1i*Glr%d%n2i*denspot%dpbox%n3d,1)*input%nspin, denspot%rhov(1), 1, rhopotOld(1), 1)

end subroutine mix_main


subroutine mixrhopotDIIS(iproc, nproc, n3d, n3p, glr, input, rhopot, rhopotold, mixdiis, alphaMix, ioffset, mixMeth, pnrm, xc)

  use module_base
  use module_types
  use module_xc
  use module_interfaces, exceptThisOne => mixrhopotDIIS
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc, n3d, n3p, mixMeth, ioffset
  type(locreg_descriptors),intent(in) :: glr
  type(input_variables),intent(in):: input
  real(8),dimension(max(glr%d%n1i*glr%d%n2i*n3d,1)*input%nspin),intent(in):: rhopotold
  real(8),dimension(max(glr%d%n1i*glr%d%n2i*n3d,1)*input%nspin),intent(out):: rhopot
  type(mixrhopotDIISParameters),intent(inout):: mixdiis
  real(8),intent(in):: alphaMix
  real(8),intent(out):: pnrm
  type(xc_info), intent(in) :: xc

  ! Local variables
  integer:: ist, jst, i, j, mi, istat, lwork, info
  integer:: mj, jj, k, jjst, isthist, ierr, iall, ndim, ndimtot
  real(8):: ddot, tt
  real(8),dimension(:,:),allocatable:: mat
  real(8),dimension(:),allocatable:: rhs, work, rhopotres
  integer,dimension(:),allocatable:: ipiv
  character(len=*),parameter:: subname='optimizeDIIS'

  call timing(iproc,'mix_DIIS      ','ON')

  ndimtot=max(glr%d%n1i*glr%d%n2i*n3d,1)*input%nspin
  ndim=max(glr%d%n1i*glr%d%n2i*n3p,1)*input%nspin

  ! Allocate the local arrays.
  mat = f_malloc0((/ mixdiis%isx+1, mixdiis%isx+1 /),id='mat')
  rhs = f_malloc0(mixdiis%isx+1,id='rhs')
  lwork=100*mixdiis%isx
  work = f_malloc(lwork,id='work')
  ipiv = f_malloc(mixdiis%isx+1,id='ipiv')
  rhopotres = f_malloc(ndimtot,id='rhopotres')

  ! Calculate the residue.
  pnrm=0.d0
  !tt=1.d0-alphaMix
  !do i=1,max(Glr%d%n1i*Glr%d%n2i*n3p,1)*input%nspin
  !do i=1,ndimtot
  do i=1,Glr%d%n1i*Glr%d%n2i*n3d
      tt = rhopotold(i) - rhopot(i)
      rhopotres(i) = alphaMix*tt
  end do
  do i=1,Glr%d%n1i*Glr%d%n2i*n3p
      tt = rhopotold(ioffset+i) - rhopot(ioffset+i)
      pnrm=pnrm+tt**2
  end do

  if (nproc > 1) then
    call mpiallred(pnrm, 1, mpi_sum, bigdft_mpi%mpi_comm)
  end if

  pnrm=sqrt(pnrm)/(glr%d%n1i*glr%d%n2i*glr%d%n3i)


  !!mat=0.d0
  !call to_zero((mixdiis%isx+1)**2, mat(1,1))
  !!rhs=0.d0
  !call to_zero(mixdiis%isx+1, rhs(1))

  ! Copy rhopot and rhopotres to the DIIS history.
  jst=(mixdiis%mis-1)*ndimtot+1
  !call vcopy(ndimtot, rhopotold, 1, mixdiis%rhopotHist, 1)
  !call vcopy(ndimtot, rhopotres, 1, mixdiis%rhopotresHist, 1)
  call vcopy(ndimtot, rhopotold(1), 1, mixdiis%rhopotHist(jst), 1)
  call vcopy(ndimtot, rhopotres(1), 1, mixdiis%rhopotresHist(jst), 1)

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
     ist=(mi-1)*ndimtot+1
     if(ndim>0) then
         !mixdiis%mat(j-i+1,min(mixdiis%isx,mixdiis%is))=ddot(ndimtot, rhopotres(1), 1, mixdiis%rhopotresHist(ist), 1)
         ! only use the non-overlapping part
         mixdiis%mat(j-i+1,min(mixdiis%isx,mixdiis%is))=ddot(ndim, rhopotres(ioffset+1), 1, mixdiis%rhopotresHist(ioffset+ist), 1)
     else
         mixdiis%mat(j-i+1,min(mixdiis%isx,mixdiis%is))=0.d0
     end if
  end do

  ! Sum up over all processes.
  if (nproc > 1) then
    call mpiallred(mixdiis%mat(1,min(mixdiis%isx,mixdiis%is)), min(mixdiis%is,mixdiis%isx), mpi_sum, bigdft_mpi%mpi_comm)
  end if

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

  !!if (iproc==0) then
  !!    do i=1,min(mixdiis%isx,mixdiis%is)+1
  !!        do j=1,min(mixdiis%isx,mixdiis%is)+1
  !!            write(*,*) 'i, j, mat(j,i)', i, j, mat(j,i)
  !!        end do
  !!    end do 
  !!end if


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

  !!do i=1,ndimtot
  !!    write(100+iproc,*) i, rhopot(i)
  !!    write(110+iproc,*) i, rhopotres(i)
  !!end do


  ! Make a new guess for the density/potential.
  ! If we are mixing the density (mixMeth==1) it is initialized to 0 or 10^-20, depending on the functional.
  ! If we are mixing the potential (mixMeth==2) it is always initialized to 0.
  if (xc_isgga(xc) .or. mixMeth==2) then
      call f_zero(rhopot)
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
      jjst=(mj-1)*ndimtot
      do k=1,ndimtot
          rhopot(k) = rhopot(k) + rhs(jj)*(mixdiis%rhopotHist(jjst+k)-mixdiis%rhopotresHist(jjst+k))
      end do
  end do



  call f_free(mat)
  call f_free(rhs)
  call f_free(work)
  call f_free(ipiv)
  call f_free(rhopotres)

  call timing(iproc,'mix_DIIS      ','OF')

end subroutine mixrhopotDIIS




subroutine initializeMixrhopotDIIS(isx, ndim, mixdiis)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: isx, ndim
type(mixrhopotDIISParameters),intent(out):: mixdiis

! Local variables
integer:: istat
character(len=*),parameter:: subname='initializeMixrhopotDIIS'

call f_routine(id='initializeMixrhopotDIIS')


mixdiis%isx=isx
mixdiis%is=0
mixdiis%mat = f_malloc_ptr((/mixdiis%isx,mixdiis%isx/),id='mixdiis%mat')
mixdiis%rhopotHist = f_malloc_ptr(mixdiis%isx*max(1,ndim),id='mixdiis%rhopotHist')
mixdiis%rhopotresHist = f_malloc_ptr(mixdiis%isx*max(1,ndim),id='mixdiis%rhopotresHist')

call f_release_routine()

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

call f_routine(id='deallocateMixrhopotDIIS')

call f_free_ptr(mixdiis%mat)
call f_free_ptr(mixdiis%rhopotHist)
call f_free_ptr(mixdiis%rhopotresHist)

call f_release_routine()

end subroutine deallocateMixrhopotDIIS
