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
subroutine optimize_coeffs(iproc, nproc, orbs, tmb, ldiis_coeff, fnrm, fnrm_crit, itmax, energy, sd_fit_curve)
  use module_base
  use module_types
  use module_interfaces
  use diis_sd_optimization
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc, itmax
  type(orbitals_data),intent(in):: orbs
  type(DFT_wavefunction),intent(inout):: tmb
  type(DIIS_obj), intent(inout) :: ldiis_coeff
  real(kind=gp),intent(in):: fnrm_crit
  real(kind=gp),intent(out):: fnrm
  real(kind=gp), intent(inout) :: energy
  logical, intent(in) :: sd_fit_curve

  ! Local variables
  integer:: iorb, jorb, istat, iall, info, iiorb, ierr, it
  real(kind=gp),dimension(:,:),allocatable:: coeffp, grad, grad_cov
  real(kind=gp) :: tt, ddot, energy0, pred_e

  call f_routine(id='optimize_coeffs')

  grad=f_malloc((/tmb%orbs%norb,orbs%norbp/), id='grad')

  if (ldiis_coeff%idsx == 0 .and. sd_fit_curve) then
     ! calculate initial energy for SD line fitting and printing (maybe don't need to (re)calculate kernel here?)
     call calculate_kernel_and_energy(iproc,nproc,tmb%linmat%denskern,tmb%linmat%ham,energy0,&
          tmb%coeff,orbs,tmb%orbs,.true.)
  else
     energy0=energy
  end if


  do it=1,itmax

     grad_cov=f_malloc((/tmb%orbs%norb,orbs%norbp/), id='grad_cov')

     call calculate_coeff_gradient(iproc,nproc,tmb,orbs,grad_cov,grad)

     ! Precondition the gradient (only making things worse...)
     !call precondition_gradient_coeff(tmb%orbs%norb, orbs%norbp, tmb%linmat%ham%matrix, tmb%linmat%ovrlp%matrix, grad)

     !For fnrm, we only sum on the occupied KS orbitals
     tt=0.d0
     do iorb=1,orbs%norbp
         tt=tt+ddot(tmb%orbs%norb, grad_cov(1,iorb), 1, grad(1,iorb), 1)
     end do
     call mpiallred(tt, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
     call f_free(grad_cov)
     fnrm=2.0_gp*tt

     call timing(iproc,'dirmin_sddiis','ON')
     coeffp=f_malloc((/tmb%orbs%norb,orbs%norbp/),id='coeffp')

     if (ldiis_coeff%idsx > 0) then !do DIIS
        !TO DO: make sure DIIS works
        ldiis_coeff%mids=mod(ldiis_coeff%ids,ldiis_coeff%idsx)+1
        ldiis_coeff%ids=ldiis_coeff%ids+1

        call dcopy(tmb%orbs%norb*orbs%norbp,tmb%coeff(1,orbs%isorb+1),1,coeffp(1,1),1)

        call diis_opt(iproc,nproc,1,0,1,(/iproc/),(/tmb%orbs%norb*orbs%norbp/),tmb%orbs%norb*orbs%norbp,&
             coeffp,grad(1,1),ldiis_coeff) 
     else  !steepest descent with curve fitting for line minimization
        if (sd_fit_curve) call find_alpha_sd(iproc,nproc,ldiis_coeff%alpha_coeff,tmb,orbs,coeffp,grad,energy0,fnrm,pred_e)

        do iorb=1,orbs%norbp
           iiorb = orbs%isorb + iorb
           do jorb=1,tmb%orbs%norb
              coeffp(jorb,iorb)=tmb%coeff(jorb,iiorb)-ldiis_coeff%alpha_coeff*grad(jorb,iorb)
           end do
        end do
     end if

     if(nproc > 1) then 
        call mpi_allgatherv(coeffp, tmb%orbs%norb*orbs%norbp, mpi_double_precision, tmb%coeff, &
           tmb%orbs%norb*orbs%norb_par(:,0), tmb%orbs%norb*orbs%isorb_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
     else
        call dcopy(tmb%orbs%norb*orbs%norb,coeffp(1,1),1,tmb%coeff(1,1),1)
     end if

     call f_free(coeffp)

     fnrm=sqrt(fnrm/dble(orbs%norb))

     ! do twice with approx S^_1/2, as not quite good enough at preserving charge if only once, but exact too expensive
     ! instead of twice could add some criterion to check accuracy?
     call reorthonormalize_coeff(iproc, nproc, orbs%norb, -8, -8, 1, tmb%orbs, tmb%linmat%ovrlp, tmb%coeff, orbs)
     call reorthonormalize_coeff(iproc, nproc, orbs%norb, -8, -8, 1, tmb%orbs, tmb%linmat%ovrlp, tmb%coeff, orbs)

     call calculate_kernel_and_energy(iproc,nproc,tmb%linmat%denskern,tmb%linmat%ham,energy,&
          tmb%coeff,orbs,tmb%orbs,.true.)
     !write(127,*) ldiis_coeff%alpha_coeff,energy
     !close(127)

     ! can't check ebs only, need to check Etot in linearScaling, but if we do it>1 without curve fitting will still need to update here
     !if (ldiis_coeff%idsx == 0 .and. (.not. sd_fit_curve) .and. energy0/=0.0_gp) then ! only update alpha after first iteration
     !   ! apply a cap so that alpha_coeff never goes below around 1.d-2 or above 2
     !   if ((energy-energy0)<0.d0 .and. ldiis_coeff%alpha_coeff < 1.8d0) then
     !      ldiis_coeff%alpha_coeff=1.1d0*ldiis_coeff%alpha_coeff
     !   else if (ldiis_coeff%alpha_coeff > 1.7d-3) then
     !      ldiis_coeff%alpha_coeff=0.5d0*ldiis_coeff%alpha_coeff
     !   end if
     !   if (iproc==0) print*,'EBSdiff,alpha',energy-energy0,ldiis_coeff%alpha_coeff,energy,energy0
     !end if

     if (iproc==0) write(*,*) ''
     if (sd_fit_curve .and. ldiis_coeff%idsx == 0) then
        if (iproc==0) write(*,'(a,I4,2x,6(ES16.6e3,2x))')'DminSD: it, fnrm, ebs, ebsdiff, alpha, pred E, diff',&
             it,fnrm,energy0,energy-energy0,ldiis_coeff%alpha_coeff,pred_e,pred_e-energy
     else if (ldiis_coeff%idsx == 0) then
        if (iproc==0) write(*,'(a,I4,2x,4(ES16.6e3,2x))')'DminSD: it, fnrm, ebs, ebsdiff, alpha',&
             it,fnrm,energy0,energy-energy0,ldiis_coeff%alpha_coeff
     else
        if (iproc==0) write(*,'(a,I4,2x,3(ES16.6e3,2x))')'DminDIIS: it, fnrm, ebs, ebsdiff',&
             it,fnrm,energy0,energy-energy0
     end if

     energy0=energy

     if (fnrm<fnrm_crit) exit

  end do

  call f_free(grad)

  call timing(iproc,'dirmin_sddiis','OF')

  call f_release_routine()

end subroutine optimize_coeffs

subroutine find_alpha_sd(iproc,nproc,alpha,tmb,orbs,coeffp,grad,energy0,fnrm,pred_e)
  use module_base
  use module_types
  use module_interfaces
  implicit none
  integer, intent(in) :: iproc, nproc
  real(kind=gp), intent(inout) :: alpha
  type(DFT_wavefunction) :: tmb
  type(orbitals_data), intent(in) :: orbs
  real(kind=gp), dimension(tmb%orbs%norb,orbs%norbp), intent(inout) :: coeffp
  real(kind=gp), dimension(tmb%orbs%norb,orbs%norbp), intent(in) :: grad
  real(kind=gp), intent(in) :: energy0, fnrm
  real(kind=gp), intent(out) :: pred_e
  integer :: iorb, iiorb, jorb, ierr
  real(kind=gp) :: tt, ddot, energy1, a, b, c, alpha_old
  real(kind=gp),dimension(:,:),allocatable :: coeff_tmp

  ! take an initial step to get 2nd point
  coeff_tmp=f_malloc((/tmb%orbs%norb,tmb%orbs%norb/),id='coeff_tmp')
  do iorb=1,orbs%norbp
     iiorb = orbs%isorb + iorb
     do jorb=1,tmb%orbs%norb
        coeffp(jorb,iorb)=tmb%coeff(jorb,iiorb)-alpha*grad(jorb,iorb)
     end do
  end do

  if(nproc > 1) then 
     call mpi_allgatherv(coeffp, tmb%orbs%norb*orbs%norbp, mpi_double_precision, coeff_tmp, &
        tmb%orbs%norb*orbs%norb_par(:,0), tmb%orbs%norb*orbs%isorb_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
  else
     call dcopy(tmb%orbs%norb*orbs%norb,coeffp(1,1),1,coeff_tmp(1,1),1)
  end if

  ! do twice with approx S^_1/2, as not quite good enough at preserving charge if only once, but exact too expensive
  ! instead of twice could add some criterion to check accuracy?
  call reorthonormalize_coeff(iproc, nproc, orbs%norb, -8, -8, 1, tmb%orbs, tmb%linmat%ovrlp, coeff_tmp, orbs)
  call reorthonormalize_coeff(iproc, nproc, orbs%norb, -8, -8, 1, tmb%orbs, tmb%linmat%ovrlp, coeff_tmp, orbs)
  call calculate_kernel_and_energy(iproc,nproc,tmb%linmat%denskern,tmb%linmat%ham,energy1,&
       coeff_tmp,orbs,tmb%orbs,.true.)
  call f_free(coeff_tmp)

  ! find ideal alpha using both points
  alpha_old=alpha
  a=fnrm/alpha_old+(energy1-energy0)/alpha_old**2
  b=-fnrm
  c=energy0
  alpha=-0.5_gp*b/a
  ! don't update if we have found a maximum, or negative alpha is predicted
  ! do something better here - don't just want to do the same thing twice, so at least check if energy has decreased
  if (alpha<0.0_gp .or. a<0.0_gp) alpha=alpha_old
  pred_e=a*alpha**2+b*alpha+c

  !open(127)
  !write(127,*) '#',a,b,c,(energy1-energy0)/alpha_old,b-(energy1-energy0)/alpha_old
  !write(127,*) 0.0_gp,energy0
  !write(127,*) alpha_old,energy1


end subroutine find_alpha_sd


subroutine calculate_kernel_and_energy(iproc,nproc,denskern,ham,energy,coeff,orbs,tmb_orbs,calculate_kernel)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc, nproc
  type(sparseMatrix), intent(in) :: ham
  type(sparseMatrix), intent(inout) :: denskern
  logical, intent(in) :: calculate_kernel
  real(kind=gp), intent(out) :: energy
  type(orbitals_data), intent(in) :: orbs, tmb_orbs
  real(kind=gp), dimension(tmb_orbs%norb,tmb_orbs%norb), intent(in) :: coeff

  integer :: iorb, jorb, ind_ham, ind_denskern, ierr, iorbp
  integer :: matrixindex_in_compressed

  if (calculate_kernel) then 
     call calculate_density_kernel(iproc, nproc, .true., orbs, tmb_orbs, coeff, denskern)
  end if

  energy=0.0_gp
  do iorbp=1,tmb_orbs%norbp
     iorb=iorbp+tmb_orbs%isorb
     do jorb=1,tmb_orbs%norb
        ind_ham = matrixindex_in_compressed(ham,iorb,jorb)
        ind_denskern = matrixindex_in_compressed(denskern,jorb,iorb)
        if (ind_ham==0.or.ind_denskern==0) cycle
        energy = energy + denskern%matrix_compr(ind_denskern)*ham%matrix_compr(ind_ham)
     end do
  end do
  if (nproc>1) then
     call mpiallred(energy, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  end if


end subroutine calculate_kernel_and_energy


! calculate grad_cov_i^a = f_i (I_ab - S_ag K^gb) H_bg c_i^d
! then grad=S^-1grad_cov
subroutine calculate_coeff_gradient(iproc,nproc,tmb,KSorbs,grad_cov,grad)
  use module_base
  use module_types
  implicit none

  integer, intent(in) :: iproc, nproc
  type(DFT_wavefunction), intent(inout) :: tmb
  type(orbitals_data), intent(in) :: KSorbs
  real(gp), dimension(tmb%orbs%norb,KSorbs%norbp), intent(out) :: grad_cov, grad  ! could make grad_cov KSorbs%norbp

  integer :: iorb, iiorb, info, ierr
  real(gp),dimension(:,:),allocatable :: sk, skh, skhp, inv_ovrlp
  integer :: matrixindex_in_compressed
  integer,dimension(:),allocatable:: ipiv
  real(kind=gp), dimension(:,:), allocatable:: grad_full
  character(len=*),parameter:: subname='calculate_coeff_gradient'

  call f_routine(id='calculate_coeff_gradient')
  call timing(iproc,'dirmin_lagmat1','ON')

  ! we have the kernel already, but need it to not contain occupations so recalculate here
  call calculate_density_kernel(iproc, nproc, .false., KSorbs, tmb%orbs, tmb%coeff, tmb%linmat%denskern)
  tmb%linmat%denskern%matrix=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='denskern')
  call uncompressMatrix(iproc,tmb%linmat%denskern)

  sk=f_malloc0((/tmb%orbs%norbp,tmb%orbs%norb/), id='sk')

  ! calculate I-S*K - first set sk to identity
  do iorb=1,tmb%orbs%norbp
     iiorb=tmb%orbs%isorb+iorb
     sk(iorb,iiorb) = 1.d0
  end do 

  if (tmb%orbs%norbp>0) then
     call dgemm('t', 'n', tmb%orbs%norbp, tmb%orbs%norb, tmb%orbs%norb, -1.d0, &
          tmb%linmat%ovrlp%matrix(1,tmb%orbs%isorb+1), tmb%orbs%norb, &
          tmb%linmat%denskern%matrix(1,1), tmb%orbs%norb, 1.d0, sk, tmb%orbs%norbp)
  end if

  ! coeffs and therefore kernel will change, so no need to keep it
  call f_free_ptr(tmb%linmat%denskern%matrix)

  skhp=f_malloc((/tmb%orbs%norb,tmb%orbs%norbp/), id='skhp')

  ! multiply by H to get (I_ab - S_ag K^gb) H_bd, or in this case the transpose of the above
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
          tmb%orbs%norb, tmb%coeff(1,KSorbs%isorb+1), tmb%orbs%norb, 0.d0, grad_cov(1,1), tmb%orbs%norb)
  end if

  call f_free(skh)

  ! multiply by f_i to get grad_i^a
  do iorb=1,KSorbs%norbp
     iiorb=KSorbs%isorb+iorb
     grad_cov(:,iorb)=grad_cov(:,iorb)*KSorbs%occup(iiorb)
  end do

  call timing(iproc,'dirmin_lagmat1','OF')
  call timing(iproc,'dirmin_dgesv','ON') !lr408t

  info = 0 ! needed for when some processors have orbs%norbp=0
  ! Solve the linear system ovrlp*grad=grad_cov
  if(tmb%orthpar%blocksize_pdsyev<0) then
     !! keep the covariant gradient to calculate fnrm correctly
     !call dcopy(tmb%orbs%norb*KSorbs%norbp,grad_cov,1,grad,1)
     !if (KSorbs%norbp>0) then
     !   ipiv=f_malloc(tmb%orbs%norb,id='ipiv')
     !   call dgesv(tmb%orbs%norb, KSorbs%norbp, tmb%linmat%ovrlp%matrix(1,1), tmb%orbs%norb, ipiv(1), &
     !        grad(1,1), tmb%orbs%norb, info)
     !   call f_free(ipiv)
     !end if
     inv_ovrlp=f_malloc((/tmb%orbs%norb,tmb%orbs%norb/),id='inv_ovrlp')
     call overlapPowerMinusOne(iproc, nproc, 1, -8, tmb%orbs%norb, tmb%linmat%ovrlp%matrix, inv_ovrlp)
     if (KSorbs%norbp>0) then
        call dgemm('n', 'n', tmb%orbs%norb, KSorbs%norbp, tmb%orbs%norb, 1.d0, inv_ovrlp(1,1), &
             tmb%orbs%norb, grad_cov(1,1), tmb%orbs%norb, 0.d0, grad(1,1), tmb%orbs%norb)
     else
        call dcopy(tmb%orbs%norb*KSorbs%norbp,grad_cov,1,grad,1)
     end if
     call f_free(inv_ovrlp)
  else
      grad_full=f_malloc((/tmb%orbs%norb,KSorbs%norb/),id='grad_full')
      ! do allgather instead of allred so we can keep grad as per proc
      if(nproc > 1) then 
         call mpi_allgatherv(grad_cov, tmb%orbs%norb*KSorbs%norbp, mpi_double_precision, grad_full, &
            tmb%orbs%norb*KSorbs%norb_par(:,0), tmb%orbs%norb*KSorbs%isorb_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
      else
         call dcopy(tmb%orbs%norb*KSorbs%norb,grad_cov(1,1),1,grad_full(1,1),1)
      end if
      !call mpiallred(grad(1,1), tmb%orbs%norb*KSorbs%norb, mpi_sum, bigdft_mpi%mpi_comm, ierr)

      call dgesv_parallel(iproc, tmb%orthpar%nproc_pdsyev, tmb%orthpar%blocksize_pdsyev, bigdft_mpi%mpi_comm, &
           tmb%orbs%norb, KSorbs%norb, tmb%linmat%ovrlp%matrix, tmb%orbs%norb, grad_full, tmb%orbs%norb, info)

      call dcopy(tmb%orbs%norb*KSorbs%norbp,grad_full(1,KSorbs%isorb+1),1,grad(1,1),1)

      call f_free(grad_full)
  end if

  if(info/=0) then
      write(*,'(a,i0)') 'ERROR in dgesv: info=',info
      stop
  end if


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

