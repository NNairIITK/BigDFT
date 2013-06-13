!> @file
!! Optimize the coefficients
!! @author
!!    Copyright (C) 2011-2012 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
!! NOTES: Coefficients are defined for Ntmb KS orbitals so as to maximize the number
!!        of orthonormality constraints. This should speedup the convergence by
!!        reducing the effective number of degrees of freedom.

subroutine optimize_coeffs(iproc, nproc, orbs, tmb, ldiis_coeff, fnrm)
  use module_base
  use module_types
  use module_interfaces
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(orbitals_data),intent(in):: orbs
  type(DFT_wavefunction),intent(inout):: tmb
  type(localizedDIISParameters),intent(inout):: ldiis_coeff
  real(8),intent(out):: fnrm

  ! Local variables
  integer:: iorb, jorb, istat, iall, info, iiorb, ierr
  real(8),dimension(:,:),allocatable:: rhs, coeffp, grad_cov
  real(8) :: tt, ddot
  character(len=*),parameter:: subname='optimize_coeffs'

  allocate(rhs(tmb%orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, rhs, 'rhs', subname)

  allocate(grad_cov(tmb%orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, grad_cov, 'grad_cov', subname)

  call calculate_coeff_gradient(iproc,nproc,tmb,orbs,grad_cov,rhs)

  ! Precondition the gradient (only making things worse...)
  !call precondition_gradient_coeff(tmb%orbs%norb, tmb%orbs%norbp, tmb%linmat%ham%matrix, tmb%linmat%ovrlp%matrix, rhs(1,orbs%isorb+1))

  call timing(iproc,'dirmin_sddiis','ON')
  ! Improve the coefficients
  if (ldiis_coeff%isx > 0) then
      ldiis_coeff%mis=mod(ldiis_coeff%is,ldiis_coeff%isx)+1
      ldiis_coeff%is=ldiis_coeff%is+1
  end if  

  if (ldiis_coeff%isx > 0) then !do DIIS
     !TO DO: make sure DIIS works
     print *,'in DIIS'
     call DIIS_coeff(iproc, orbs, tmb, rhs(1,orbs%isorb+1), tmb%coeff, ldiis_coeff)
  else  !steepest descent
     allocate(coeffp(tmb%orbs%norb,orbs%norbp),stat=istat)
     call memocc(istat, coeffp, 'coeffp', subname)
     do iorb=1,orbs%norbp
        iiorb = orbs%isorb + iorb
        do jorb=1,tmb%orbs%norb
           coeffp(jorb,iorb)=tmb%coeff(jorb,iiorb)-ldiis_coeff%alpha_coeff*rhs(jorb,iiorb)
        end do
     end do

     if(nproc > 1) then 
        call mpi_allgatherv(coeffp, tmb%orbs%norb*orbs%norbp, mpi_double_precision, tmb%coeff, &
           tmb%orbs%norb*orbs%norb_par(:,0), tmb%orbs%norb*orbs%isorb_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
     else
        call dcopy(tmb%orbs%norb*orbs%norb,coeffp(1,1),1,tmb%coeff(1,1),1)
     end if
     iall=-product(shape(coeffp))*kind(coeffp)
     deallocate(coeffp, stat=istat)
     call memocc(istat, iall, 'coeffp', subname)
  end if

  !For fnrm, we only sum on the occupied KS orbitals
  tt=0.d0
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      tt=tt+ddot(tmb%orbs%norb, grad_cov(1,iiorb), 1, rhs(1,iiorb), 1)
  end do
  call mpiallred(tt, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  fnrm=sqrt(2.d0*tt/dble(orbs%norb))

  iall=-product(shape(grad_cov))*kind(grad_cov)
  deallocate(grad_cov, stat=istat)
  call memocc(istat, iall, 'grad_cov', subname)

  iall=-product(shape(rhs))*kind(rhs)
  deallocate(rhs, stat=istat)
  call memocc(istat, iall, 'rhs', subname)

  call timing(iproc,'dirmin_sddiis','OF')

  ! do twice with approx S^_1/2, as not quite good enough at preserving charge if only once, but exact too expensive
  ! instead of twice could add some criterion to check accuracy?
  call reorthonormalize_coeff(iproc, nproc, orbs, -8, -8, 1, tmb%orbs, tmb%linmat%ovrlp, tmb%coeff)
  call reorthonormalize_coeff(iproc, nproc, orbs, -8, -8, 1, tmb%orbs, tmb%linmat%ovrlp, tmb%coeff)


end subroutine optimize_coeffs


! calculate grad_cov_i^a = f_i (I_ab - S_ag K^gb) H_bg c_i^d
! then grad=S^-1grad_cov
subroutine calculate_coeff_gradient(iproc,nproc,tmb,KSorbs,grad_cov,grad)
  use module_base
  use module_types
  implicit none

  integer, intent(in) :: iproc, nproc
  type(DFT_wavefunction), intent(inout) :: tmb
  type(orbitals_data), intent(in) :: KSorbs
  real(gp), dimension(tmb%orbs%norb,KSorbs%norb), intent(out) :: grad_cov, grad  ! could make grad_cov KSorbs%norbp

  integer :: iorb, iiorb, info, ierr
  real(gp),dimension(:,:),allocatable :: sk, skh, skhp
  integer,dimension(:),allocatable:: ipiv
  character(len=*),parameter:: subname='calculate_coeff_gradient'

  call f_routine(id='calculate_coeff_gradient')
  call timing(iproc,'dirmin_lagmat1','ON')

  ! we have the kernel already, but need it to not contain occupations so recalculate here
  tmb%linmat%denskern%matrix=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='denskern')

  call calculate_density_kernel(iproc, nproc, .false., KSorbs, tmb%orbs, tmb%coeff, tmb%linmat%denskern%matrix)

  sk=f_malloc0((/tmb%orbs%norbp,tmb%orbs%norb/), id='sk')

  ! calculate I-S*K - first set sk to identity
  do iorb=1,tmb%orbs%norbp
     iiorb=tmb%orbs%isorb+iorb
     sk(iorb,iiorb) = 1.d0
  end do 

  if (tmb%orbs%norbp>0) then
     call dgemm('n', 'n', tmb%orbs%norbp, tmb%orbs%norb, tmb%orbs%norb, -1.d0, &
          tmb%linmat%ovrlp%matrix(tmb%orbs%isorb+1:tmb%orbs%isorb+tmb%orbs%norbp,:),&
          tmb%orbs%norbp, tmb%linmat%denskern%matrix(1,1), tmb%orbs%norb, 1.d0, sk, tmb%orbs%norbp)
  end if

  ! coeffs and therefore kernel will change, so no need to keep it
  call f_free_ptr(tmb%linmat%denskern%matrix)

  skhp=f_malloc((/tmb%orbs%norb,tmb%orbs%norbp/), id='skhp')

  ! multiply by H to get (I_ab - S_ag K^gb) H_bg, or in this case the transpose of the above
  if (tmb%orbs%norbp>0) then
     call dgemm('t', 't', tmb%orbs%norb, tmb%orbs%norbp, tmb%orbs%norb, 1.d0, tmb%linmat%ham%matrix(1,1), &
          tmb%orbs%norb, sk(1,1), tmb%orbs%norbp, 0.d0, skhp(1,1), tmb%orbs%norb)
  end if

  call f_free(sk)

  skh=f_malloc((/tmb%orbs%norb,tmb%orbs%norb/), id='skh')

  ! gather together
  if(nproc > 1) then
     call mpi_allgatherv(skhp(1,1), tmb%orbs%norb*tmb%orbs%norbp, mpi_double_precision, skh(1,1), &
        tmb%orbs%norb*tmb%orbs%norb_par(:,0), tmb%orbs%norb*tmb%orbs%isorb_par, &
        mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
  else
     call dcopy(tmb%orbs%norbp*tmb%orbs%norb,skhp(1,1),1,skh(1,1),1)
  end if

  call f_free(skhp)

  ! calc for i on this proc: (I_ab - S_ag K^gb) H_bg c_i^d
  if (KSorbs%norbp>0) then
     call dgemm('t', 'n', tmb%orbs%norb, KSorbs%norbp, tmb%orbs%norb, 1.d0, skh(1,1), &
          tmb%orbs%norb, tmb%coeff(1,KSorbs%isorb+1), tmb%orbs%norb, 0.d0, grad_cov(1,KSorbs%isorb+1), tmb%orbs%norb)
  end if

  call f_free(skh)

  ! multiply by f_i to get grad_i^a
  do iorb=1,KSorbs%norbp
     iiorb=KSorbs%isorb+iorb
     grad_cov(:,iiorb)=grad_cov(:,iiorb)*KSorbs%occup(iiorb)
  end do

  ! keep the covariant gradient to calculate fnrm correctly
  call dcopy(tmb%orbs%norb*KSorbs%norb,grad_cov,1,grad,1)

  call timing(iproc,'dirmin_lagmat1','OF')
  call timing(iproc,'dirmin_dgesv','ON') !lr408t

  info = 0 ! needed for when some processors have orbs%norbp=0
  ipiv=f_malloc(tmb%orbs%norb,id='ipiv')

  ! Solve the linear system ovrlp*grad=grad_cov
  if(tmb%orthpar%blocksize_pdsyev<0) then
      if (KSorbs%norbp>0) then
          call dgesv(tmb%orbs%norb, KSorbs%norbp, tmb%linmat%ovrlp%matrix(1,1), tmb%orbs%norb, ipiv(1), &
               grad(1,KSorbs%isorb+1), tmb%orbs%norb, info)
      end if
  else
      call mpiallred(grad(1,1), tmb%orbs%norb*KSorbs%norb, mpi_sum, bigdft_mpi%mpi_comm, ierr)
      call dgesv_parallel(iproc, tmb%orthpar%nproc_pdsyev, tmb%orthpar%blocksize_pdsyev, bigdft_mpi%mpi_comm, &
           tmb%orbs%norb, KSorbs%norb, tmb%linmat%ovrlp%matrix, tmb%orbs%norb, grad, tmb%orbs%norb, info)
  end if

  if(info/=0) then
      write(*,'(a,i0)') 'ERROR in dgesv: info=',info
      stop
  end if

  call f_free(ipiv)
  call timing(iproc,'dirmin_dgesv','OF') !lr408t
  call f_release_routine()

end subroutine calculate_coeff_gradient



subroutine precondition_gradient_coeff(ntmb, norb, ham, ovrlp, grad)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: ntmb, norb
  real(8),dimension(ntmb,ntmb),intent(in):: ham, ovrlp
  real(8),dimension(ntmb,norb),intent(inout):: grad
  
  ! Local variables
  integer:: iorb, itmb, jtmb, info, istat, iall
  complex(8),dimension(:,:),allocatable:: mat
  complex(8),dimension(:,:),allocatable:: rhs
  integer,dimension(:),allocatable:: ipiv
  character(len=*),parameter:: subname='precondition_gradient_coeff'
  
  allocate(mat(ntmb,ntmb), stat=istat)
  !call memocc(istat, mat, 'mat', subname)
  allocate(rhs(ntmb,norb), stat=istat)
  !call memocc(istat, mat, 'mat', subname)
  
  ! Build the matrix to be inverted
  do itmb=1,ntmb
      do jtmb=1,ntmb
          mat(jtmb,itmb) = cmplx(ham(jtmb,itmb)+.5d0*ovrlp(jtmb,itmb),0.d0,kind=8)
      end do
      mat(itmb,itmb)=mat(itmb,itmb)+cmplx(0.d0,-1.d-1,kind=8)
      !mat(itmb,itmb)=mat(itmb,itmb)-cprec
  end do
  do iorb=1,norb
      do itmb=1,ntmb
          rhs(itmb,iorb)=cmplx(grad(itmb,iorb),0.d0,kind=8)
      end do
  end do
  
  
  allocate(ipiv(ntmb), stat=istat)
  call memocc(istat, ipiv, 'ipiv', subname)
  
  call zgesv(ntmb, norb, mat(1,1), ntmb, ipiv, rhs(1,1), ntmb, info)
  if(info/=0) then
      stop 'ERROR in dgesv'
  end if
  !call dcopy(nel, rhs(1), 1, grad(1), 1)
  do iorb=1,norb
      do itmb=1,ntmb
          grad(itmb,iorb)=real(rhs(itmb,iorb))
      end do
  end do
  
  iall=-product(shape(ipiv))*kind(ipiv)
  deallocate(ipiv, stat=istat)
  call memocc(istat, iall, 'ipiv', subname)
  
  iall=-product(shape(mat))*kind(mat)
  deallocate(mat, stat=istat)
  !call memocc(istat, iall, 'mat', subname)
  
  iall=-product(shape(rhs))*kind(rhs)
  deallocate(rhs, stat=istat)
  !call memocc(istat, iall, 'rhs', subname)

end subroutine precondition_gradient_coeff



subroutine DIIS_coeff(iproc, orbs, tmb, grad, coeff, ldiis)
  use module_base
  use module_types
  use module_interfaces, except_this_one => DIIS_coeff
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc
  type(orbitals_data),intent(in):: orbs
  type(DFT_wavefunction),intent(in):: tmb
  real(8),dimension(tmb%orbs%norb*tmb%orbs%norbp),intent(in):: grad
  real(8),dimension(tmb%orbs%norb*tmb%orbs%norb),intent(inout):: coeff
  type(localizedDIISParameters),intent(inout):: ldiis
  
  ! Local variables
  integer:: iorb, jorb, ist, ncount, jst, i, j, mi, ist1, ist2, istat, lwork, info
  integer:: mj, jj, k, jjst, isthist, iall
  real(8):: ddot
  real(8),dimension(:,:),allocatable:: mat
  real(8),dimension(:),allocatable:: rhs, work
  integer,dimension(:),allocatable:: ipiv
  character(len=*),parameter:: subname='DIIS_coeff'
  
  !!call timing(iproc,'optimize_DIIS ','ON')
  
  ! Allocate the local arrays.
  allocate(mat(ldiis%isx+1,ldiis%isx+1), stat=istat)
  call memocc(istat, mat, 'mat', subname)
  allocate(rhs(ldiis%isx+1), stat=istat)
  call memocc(istat, rhs, 'rhs', subname)
  allocate(ipiv(ldiis%isx+1), stat=istat)
  call memocc(istat, ipiv, 'ipiv', subname)
  
  mat=0.d0
  rhs=0.d0
  call to_zero((ldiis%isx+1)**2, mat(1,1))
  call to_zero(ldiis%isx+1, rhs(1))
  
  ncount=tmb%orbs%norb

  ! Copy coeff and grad to history.
  ist=1
  do iorb=1,tmb%orbs%norbp
      jst=1
      do jorb=1,iorb-1
          jst=jst+ncount*ldiis%isx
      end do
      jst=jst+(ldiis%mis-1)*ncount
      call dcopy(ncount, coeff(ist+tmb%orbs%isorb*tmb%orbs%norb), 1, ldiis%phiHist(jst), 1)
      call dcopy(ncount, grad(ist), 1, ldiis%hphiHist(jst), 1)
      ist=ist+ncount
  end do
  
  do iorb=1,tmb%orbs%norbp
      ! Shift the DIIS matrix left up if we reached the maximal history length.
      if(ldiis%is>ldiis%isx) then
         do i=1,ldiis%isx-1
            do j=1,i
               ldiis%mat(j,i,iorb)=ldiis%mat(j+1,i+1,iorb)
            end do
         end do
      end if
  end do
  
  do iorb=1,tmb%orbs%norbp
      ! Calculate a new line for the matrix.
      i=max(1,ldiis%is-ldiis%isx+1)
      jst=1
      ist1=1
      do jorb=1,iorb-1
          jst=jst+ncount*ldiis%isx
          ist1=ist1+ncount
      end do
      do j=i,ldiis%is
         mi=mod(j-1,ldiis%isx)+1
         ist2=jst+(mi-1)*ncount
         if(ist2>size(ldiis%hphiHist)) then
             write(*,'(a,7i8)') 'ERROR ist2: iproc, iorb, ldiis%is, mi, ncount, ist2, size(ldiis%hphiHist)', iproc, iorb, ldiis%is,&
                                 mi, ncount, ist2, size(ldiis%hphiHist)
         end if
         ldiis%mat(j-i+1,min(ldiis%isx,ldiis%is),iorb)=ddot(ncount, grad(ist1), 1, ldiis%hphiHist(ist2), 1)
         ist2=ist2+ncount
      end do
  end do
  
  ist=1+tmb%orbs%isorb*tmb%orbs%norb
  do iorb=1,tmb%orbs%norbp
      ! Copy the matrix to an auxiliary array and fill with the zeros and ones.
      do i=1,min(ldiis%isx,ldiis%is)
          mat(i,min(ldiis%isx,ldiis%is)+1)=1.d0
          rhs(i)=0.d0
          do j=i,min(ldiis%isx,ldiis%is)
              mat(i,j)=ldiis%mat(i,j,iorb)
          end do
      end do
      mat(min(ldiis%isx,ldiis%is)+1,min(ldiis%isx,ldiis%is)+1)=0.d0
      rhs(min(ldiis%isx,ldiis%is)+1)=1.d0
   
      ! Solve the linear system
      !!do istat=1,ldiis%isx+1
          !!do iall=1,ldiis%isx+1
              !!if(iproc==0) write(500,*) istat, iall, mat(iall,istat)
          !!end do
      !!end do

      if(ldiis%is>1) then
         lwork=-1   !100*ldiis%isx
         allocate(work(1000), stat=istat)
         call memocc(istat, work, 'work', subname)
         call dsysv('u', min(ldiis%isx,ldiis%is)+1, 1, mat, ldiis%isx+1,  & 
              ipiv, rhs(1), ldiis%isx+1, work, lwork, info)
         lwork=nint(work(1))
         iall=-product(shape(work))*kind(work)
         deallocate(work,stat=istat)
         call memocc(istat,iall,'work',subname)
         allocate(work(lwork), stat=istat)
         call memocc(istat, work, 'work', subname)
         call dsysv('u', min(ldiis%isx,ldiis%is)+1, 1, mat, ldiis%isx+1,  & 
              ipiv, rhs(1), ldiis%isx+1, work, lwork, info)
         iall=-product(shape(work))*kind(work)
         deallocate(work, stat=istat)
         call memocc(istat, iall, 'work', subname)
         
         if (info /= 0) then
            write(*,'(a,i0)') 'ERROR in dsysv (DIIS_coeff), info=', info
            stop
         end if
      else
         rhs(1)=1.d0
      endif
    
      ! Make a new guess for the orbital.
      call razero(ncount, coeff(ist))
      isthist=max(1,ldiis%is-ldiis%isx+1)
      jj=0
      jst=0
      do jorb=1,iorb-1
          jst=jst+ncount*ldiis%isx
      end do
      do j=isthist,ldiis%is
          jj=jj+1
          mj=mod(j-1,ldiis%isx)+1
          jjst=jst+(mj-1)*ncount
          do k=1,ncount
              coeff(ist+k-1) = coeff(ist+k-1) + rhs(jj)*(ldiis%phiHist(jjst+k)-ldiis%hphiHist(jjst+k))
          end do
      end do
      ist=ist+ncount
  end do
    
  iall=-product(shape(mat))*kind(mat)
  deallocate(mat, stat=istat)
  call memocc(istat, iall, 'mat', subname)
  
  iall=-product(shape(rhs))*kind(rhs)
  deallocate(rhs, stat=istat)
  call memocc(istat, iall, 'rhs', subname)

  iall=-product(shape(ipiv))*kind(ipiv)
  deallocate(ipiv, stat=istat)
  call memocc(istat, iall, 'ipiv', subname)
  
  !!call timing(iproc,'optimize_DIIS ','OF')

end subroutine DIIS_coeff


subroutine initialize_DIIS_coeff(isx, ldiis)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: isx
  type(localizedDIISParameters),intent(inout):: ldiis
  
  ! Local variables
  character(len=*),parameter:: subname='initialize_DIIS_coeff'
    
  ldiis%isx=isx
  ldiis%is=0
  ldiis%switchSD=.false.
  ldiis%trmin=1.d100
  ldiis%trold=1.d100
  ldiis%alpha_coeff=0.1d0

end subroutine initialize_DIIS_coeff


subroutine allocate_DIIS_coeff(tmb, ldiis)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  type(DFT_wavefunction),intent(in):: tmb
  type(localizedDIISParameters),intent(inout):: ldiis
  
  ! Local variables
  integer:: ii, istat
  character(len=*),parameter:: subname='allocate_DIIS_coeff'

  allocate(ldiis%mat(ldiis%isx,ldiis%isx,tmb%orbs%norbp),stat=istat)
  call memocc(istat, ldiis%mat, 'ldiis%mat', subname)

  ii=ldiis%isx*tmb%orbs%norb*tmb%orbs%norbp
  allocate(ldiis%phiHist(ii), stat=istat)
  call memocc(istat, ldiis%phiHist, 'ldiis%phiHist', subname)
  allocate(ldiis%hphiHist(ii), stat=istat)
  call memocc(istat, ldiis%hphiHist, 'ldiis%hphiHist', subname)

end subroutine allocate_DIIS_coeff

