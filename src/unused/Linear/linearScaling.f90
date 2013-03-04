subroutine derivatives_with_orthoconstraint(iproc, nproc, tmb, tmbder)
  use module_base
  use module_types
  use module_interfaces, except_this_one => derivatives_with_orthoconstraint
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(DFT_wavefunction),intent(in) :: tmb
  type(DFT_wavefunction),intent(inout) :: tmbder

  ! Local variables
  integer :: i0, j0, ii, jj, ipt, i, iiorb, jjorb, istat, iall, j
  real(8),dimension(:),allocatable :: psit_c, psit_f, psidert_c, psidert_f
  real(8),dimension(:,:),allocatable :: matrix
  character(len=*),parameter :: subname='derivatives_with_orthoconstraint'


  write(*,*) 'WARNING: in derivatives_with_orthoconstraint'

  allocate(psit_c(tmb%collcom%ndimind_c), stat=istat)
  call memocc(istat, psit_c, 'psit_c', subname)
  allocate(psit_f(7*tmb%collcom%ndimind_f), stat=istat)
  call memocc(istat, psit_f, 'psit_f', subname)

  allocate(psidert_c(tmbder%collcom%ndimind_c), stat=istat)
  call memocc(istat, psidert_c, 'psidert_c', subname)
  allocate(psidert_f(7*tmbder%collcom%ndimind_f), stat=istat)
  call memocc(istat, psidert_f, 'psidert_f', subname)

  do istat=1,size(tmbder%psi)
      write(200,*) istat, tmbder%psi(istat)
  end do

  ! Transpose the support functions
  call transpose_localized(iproc, nproc, tmb%orbs%npsidim_orbs, tmb%orbs, tmb%collcom, &
       tmb%psi, psit_c, psit_f, tmb%lzd)
  do istat=1,size(psit_c)
      write(201,*) psit_c(istat)
  end do
  do istat=1,size(psit_f)
      write(201,*) psit_f(istat)
  end do

  ! Transpose the derivatives
  call transpose_localized(iproc, nproc, tmbder%orbs%npsidim_orbs, tmb%orbs, tmbder%collcom, &
       tmbder%psi, psidert_c, psidert_f, tmb%lzd)


  allocate(matrix(tmbder%orbs%norb,tmb%orbs%norb), stat=istat)
  call memocc(istat, matrix, 'matrix', subname)

  ! Calculate the matrix <dphi_i|phi_j>
  call calculate_pulay_overlap(iproc, nproc, tmbder%orbs, tmb%orbs, tmbder%collcom, &
       tmb%collcom, psidert_c, psit_c, psidert_f, psit_f, matrix)
  do i=1,tmb%orbs%norb
      do j=1,tmbder%orbs%norb
          if(iproc==0) write(400,*) i,j,matrix(j,i)
      end do
  end do


  i0=0
  j0=0
  do ipt=1,tmb%collcom%nptsp_c 
      ii=tmb%collcom%norb_per_gridpoint_c(ipt) 
      jj=tmbder%collcom%norb_per_gridpoint_c(ipt) 
      do i=1,jj
          jjorb=tmbder%collcom%indexrecvorbital_c(j0+i)
          do j=1,ii
              iiorb=tmb%collcom%indexrecvorbital_c(i0+j)
              write(333,'(3i8,3es15.5)') jjorb, iiorb, ceiling(dble(jjorb)/3.d0), &
                                         5d0*matrix(jjorb,iiorb), psidert_c(j0+i), psit_c(i0+j)
              psidert_c(j0+i)=psidert_c(j0+i)-.5d0*matrix(jjorb,iiorb)*psit_c(i0+j)
              if (iiorb==ceiling(dble(jjorb)/3.d0)) then
                  psidert_c(j0+i)=psidert_c(j0+i)-.5d0*matrix(iiorb,jjorb)*psit_c(i0+j)
              end if
          end do
      end do
      i0=i0+ii
      j0=j0+jj
  end do

  i0=0
  j0=0
  do ipt=1,tmb%collcom%nptsp_f 
      ii=tmb%collcom%norb_per_gridpoint_f(ipt) 
      jj=tmbder%collcom%norb_per_gridpoint_f(ipt) 
      do i=1,jj
          jjorb=tmbder%collcom%indexrecvorbital_f(j0+i)
          do j=1,ii
              iiorb=tmb%collcom%indexrecvorbital_f(i0+j)
              psidert_f(7*(j0+i)-6)=psidert_f(7*(j0+i)-6)-.5d0*matrix(jjorb,iiorb)*psit_f(7*(i0+j)-6)
              psidert_f(7*(j0+i)-5)=psidert_f(7*(j0+i)-5)-.5d0*matrix(jjorb,iiorb)*psit_f(7*(i0+j)-5)
              psidert_f(7*(j0+i)-4)=psidert_f(7*(j0+i)-4)-.5d0*matrix(jjorb,iiorb)*psit_f(7*(i0+j)-4)
              psidert_f(7*(j0+i)-3)=psidert_f(7*(j0+i)-3)-.5d0*matrix(jjorb,iiorb)*psit_f(7*(i0+j)-3)
              psidert_f(7*(j0+i)-2)=psidert_f(7*(j0+i)-2)-.5d0*matrix(jjorb,iiorb)*psit_f(7*(i0+j)-2)
              psidert_f(7*(j0+i)-1)=psidert_f(7*(j0+i)-1)-.5d0*matrix(jjorb,iiorb)*psit_f(7*(i0+j)-1)
              psidert_f(7*(j0+i)-0)=psidert_f(7*(j0+i)-0)-.5d0*matrix(jjorb,iiorb)*psit_f(7*(i0+j)-0)
              if (iiorb==ceiling(dble(jjorb)/3.d0)) then
                  psidert_f(7*(j0+i)-6)=psidert_f(7*(j0+i)-6)-.5d0*matrix(iiorb,jjorb)*psit_f(7*(i0+j)-6)
                  psidert_f(7*(j0+i)-5)=psidert_f(7*(j0+i)-5)-.5d0*matrix(iiorb,jjorb)*psit_f(7*(i0+j)-5)
                  psidert_f(7*(j0+i)-4)=psidert_f(7*(j0+i)-4)-.5d0*matrix(iiorb,jjorb)*psit_f(7*(i0+j)-4)
                  psidert_f(7*(j0+i)-3)=psidert_f(7*(j0+i)-3)-.5d0*matrix(iiorb,jjorb)*psit_f(7*(i0+j)-3)
                  psidert_f(7*(j0+i)-2)=psidert_f(7*(j0+i)-2)-.5d0*matrix(iiorb,jjorb)*psit_f(7*(i0+j)-2)
                  psidert_f(7*(j0+i)-1)=psidert_f(7*(j0+i)-1)-.5d0*matrix(iiorb,jjorb)*psit_f(7*(i0+j)-1)
                  psidert_f(7*(j0+i)-0)=psidert_f(7*(j0+i)-0)-.5d0*matrix(iiorb,jjorb)*psit_f(7*(i0+j)-0)
              end if
          end do
      end do
      i0=i0+ii
      j0=j0+jj
  end do

  !! TEST ONLY
  call calculate_pulay_overlap(iproc, nproc, tmbder%orbs, tmb%orbs, tmbder%collcom, &
       tmb%collcom, psidert_c, psit_c, psidert_f, psit_f, matrix)
  do i=1,tmb%orbs%norb
      do j=1,tmbder%orbs%norb
          if(iproc==0) write(450,*) i,j,matrix(j,i)
      end do
  end do

  !!do istat=1,size(tmbder%psi)
  !!    write(200+iproc,*) istat, tmbder%psi(istat)
  !!end do

  ! Untranpose the derivatives
  call untranspose_localized(iproc, nproc, tmbder%orbs%npsidim_orbs, tmbder%orbs, tmbder%collcom, &
       psidert_c, psidert_f, tmbder%psi, tmb%lzd)
  !!do istat=1,size(tmbder%psi)
  !!    write(300+iproc,*) istat, tmbder%psi(istat)
  !!end do
  do istat=1,size(tmbder%psi)
      write(250,*) istat, tmbder%psi(istat)
  end do

  iall=-product(shape(matrix))*kind(matrix)
  deallocate(matrix,stat=istat)
  call memocc(istat,iall,'matrix',subname)

  iall=-product(shape(psit_c))*kind(psit_c)
  deallocate(psit_c,stat=istat)
  call memocc(istat,iall,'psit_c',subname)

  iall=-product(shape(psit_f))*kind(psit_f)
  deallocate(psit_f,stat=istat)
  call memocc(istat,iall,'psit_f',subname)

  iall=-product(shape(psidert_c))*kind(psidert_c)
  deallocate(psidert_c,stat=istat)
  call memocc(istat,iall,'psidert_c',subname)

  iall=-product(shape(psidert_f))*kind(psidert_f)
  deallocate(psidert_f,stat=istat)
  call memocc(istat,iall,'psidert_f',subname)

end subroutine derivatives_with_orthoconstraint


subroutine transformToGlobal(iproc,nproc,lzd,lorbs,orbs,comms,input,coeff,lphi,psi,psit)
use module_base
use module_types
use module_interfaces, exceptThisOne => transformToGlobal
implicit none

! Calling arguments
integer,intent(in) :: iproc, nproc
type(local_zone_descriptors),intent(in) :: lzd
type(orbitals_data),intent(in) :: lorbs, orbs
type(communications_arrays) :: comms
type(input_variables),intent(in) :: input
real(8),dimension(lorbs%norb,orbs%norb),intent(in) :: coeff
real(8),dimension(max(lorbs%npsidim_orbs,lorbs%npsidim_comp)),intent(inout) :: lphi
real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),target,intent(out) :: psi
real(8),dimension(:),pointer,intent(inout) :: psit

! Local variables
integer :: ind1, ind2, istat, iall, iorb, ilr, ldim, gdim, nvctrp
real(8),dimension(:),pointer :: phiWork
real(8),dimension(:),allocatable :: phi
character(len=*),parameter :: subname='transformToGlobal'
type(orbitals_data) :: gorbs
type(communications_arrays) :: gcomms

  call nullify_orbitals_data(gorbs)
  call copy_orbitals_data(lorbs, gorbs, subname)
  call orbitals_communicators(iproc,nproc,lzd%glr,gorbs,gcomms)

  allocate(phi(max(gorbs%npsidim_orbs,gorbs%npsidim_comp)+ndebug), stat=istat)
  call memocc(istat, phi, 'phi', subname)
  allocate(phiWork(max(size(phi),size(psi))), stat=istat)
  call memocc(istat, phiWork, 'phiWork', subname)

  ind1=1
  ind2=1
  if (max(gorbs%npsidim_orbs,gorbs%npsidim_comp) > 0) &
       call to_zero(max(gorbs%npsidim_orbs,gorbs%npsidim_comp),phi(1))

  do iorb=1,lorbs%norbp
      ilr = lorbs%inWhichLocreg(lorbs%isorb+iorb)
      ldim=lzd%Llr(ilr)%wfd%nvctr_c+7*lzd%Llr(ilr)%wfd%nvctr_f
      gdim=lzd%Glr%wfd%nvctr_c+7*lzd%Glr%wfd%nvctr_f

      call Lpsi_to_global2(iproc,ldim,gdim,lorbs%norb,lorbs%nspinor,input%nspin,lzd%Glr,&
           lzd%Llr(ilr),lphi(ind2),phi(ind1))
      ind1=ind1+lzd%Glr%wfd%nvctr_c+7*lzd%Glr%wfd%nvctr_f
      ind2=ind2+lzd%Llr(ilr)%wfd%nvctr_c+7*lzd%Llr(ilr)%wfd%nvctr_f

  end do
  call transpose_v(iproc, nproc, gorbs, lzd%Glr%wfd, gcomms, phi, work=phiWork)

  if(iproc==0) then
      write(*,'(1x,a)', advance='no') '------------------------------------- Building linear combinations... '
  end if
  ! Build the extended orbital psi as a linear combination of localized basis functions phi. for real O(N)
  ! this has to replaced, but at the moment it is still needed.
  !call buildWavefunctionModified(iproc, nproc, orbs, gorbs, comms, gcomms, phi, psi, coeff)

  nvctrp=sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor
  call dgemm('n', 'n', nvctrp, orbs%norb, lorbs%norb, 1.d0, phi(1), nvctrp, coeff(1,1), &
       lorbs%norb, 0.d0, psi(1), nvctrp)

  ! not used in linearscaling
  !if(nproc>1) then
  !    call dcopy(orbs%npsidim_comp, psi, 1, psit, 1)
  !else
  !    psit => psi
  !end if

  call untranspose_v(iproc, nproc, orbs, lzd%Glr%wfd, comms, psi, work=phiWork)

  if(iproc==0) write(*,'(a)') 'done.'


  iall=-product(shape(phi))*kind(phi)
  deallocate(phi, stat=istat)
  call memocc(istat, iall, 'phi', subname)
  iall=-product(shape(phiWork))*kind(phiWork)
  deallocate(phiWork, stat=istat)
  call memocc(istat, iall, 'phiWork', subname)

  call deallocate_orbitals_data(gorbs, subname)
  call deallocate_communications_arrays(gcomms, subname)


end subroutine transformToGlobal



