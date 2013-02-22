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

subroutine optimize_coeffs_sparse(iproc, nproc, orbs, tmb, ldiis_coeff, fnrm)
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
  integer:: iorb, jorb, korb, lorb, istat, iall, info, iiorb, ierr, ind, indh, indo, kkorb
  integer :: npts_per_proc, ind_start, ind_end, indc, iseg, segn
  real(8),dimension(:,:),allocatable:: rhs, coeffp, sk, skh, skhp !gradp, lagmat, ovrlp_tmp, ovrlp_coeff
  integer,dimension(:),allocatable:: ipiv
  real(8) :: tt, ddot, tt2
  logical :: dense
  character(len=*),parameter:: subname='optimize_coeffs'

  ! we have the kernel already, but need it to not contain occupations so recalculate here
  allocate(tmb%linmat%denskern%matrix(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
  call memocc(istat, tmb%linmat%denskern%matrix, 'tmb%linmat%denskern%matrix', subname)
  call calculate_density_kernel(iproc, nproc, .false., orbs, tmb%orbs, tmb%coeff, tmb%linmat%denskern%matrix)

  call timing(iproc,'dirmin_lagmat1','ON')

  dense=.false.

  ! calculate rhs_i^a = f_i (I_ab - S_ag K^gb) H_bg c_i^d
  if (dense) then ! non-parallelized for now wrt tmb%orbs%norb
    
     allocate(sk(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
     call memocc(istat, sk, 'sk', subname)

     ! calculate I-S*K - first set sk to identity
     call to_zero(tmb%orbs%norb*tmb%orbs%norb, sk(1,1))
     do iorb=1,tmb%orbs%norb
        sk(iorb,iorb) = 1.d0
     end do 
 
     call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norb, tmb%orbs%norb, -1.d0, tmb%linmat%ovrlp%matrix(1,1), &
          tmb%orbs%norb, tmb%linmat%denskern%matrix(1,1), tmb%orbs%norb, 1.d0, sk, tmb%orbs%norb)

     ! coeffs and therefore kernel will change, so no need to keep it
     iall=-product(shape(tmb%linmat%denskern%matrix))*kind(tmb%linmat%denskern%matrix)
     deallocate(tmb%linmat%denskern%matrix, stat=istat)
     call memocc(istat, iall, 'tmb%linmat%denskern%matrix', subname)

     allocate(skh(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
     call memocc(istat, skh, 'skh', subname)

     ! multiply by H to get (I_ab - S_ag K^gb) H_bg
     call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norb, tmb%orbs%norb, 1.d0, sk(1,1), &
          tmb%orbs%norb, tmb%linmat%ham%matrix(1,1), tmb%orbs%norb, 0.d0, skh, tmb%orbs%norb)

     iall=-product(shape(sk))*kind(sk)
     deallocate(sk, stat=istat)
     call memocc(istat, iall, 'sk', subname)

     allocate(rhs(tmb%orbs%norb,orbs%norb), stat=istat)
     call memocc(istat, rhs, 'rhs', subname)

     ! calc for i on this proc: (I_ab - S_ag K^gb) H_bg c_i^d
     call dgemm('n', 'n', tmb%orbs%norb, orbs%norbp, tmb%orbs%norb, 1.d0, skh(1,1), &
          tmb%orbs%norb, tmb%coeff(1,orbs%isorb+1), tmb%orbs%norb, 0.d0, rhs(1,orbs%isorb+1), tmb%orbs%norb)

     iall=-product(shape(skh))*kind(skh)
     deallocate(skh, stat=istat)
     call memocc(istat, iall, 'skh', subname)

     ! multiply by f_i to get rhs_i^a
     do iorb=1,orbs%norbp
        iiorb=orbs%isorb+iorb
        rhs(:,iiorb)=rhs(:,iiorb)*orbs%occup(iiorb)
     end do

  else ! parallel version wrt tmb%orbs%norb as well as orbs%norb

     allocate(sk(tmb%orbs%norbp,tmb%orbs%norb), stat=istat)
     call memocc(istat, sk, 'sk', subname)

     ! calculate I-S*K - first set sk to identity
     call to_zero(tmb%orbs%norbp*tmb%orbs%norb, sk(1,1))
     do iorb=1,tmb%orbs%norbp
        iiorb=tmb%orbs%isorb+iorb
        sk(iorb,iiorb) = 1.d0
     end do 
 
     call dgemm('n', 'n', tmb%orbs%norbp, tmb%orbs%norb, tmb%orbs%norb, -1.d0, &
          tmb%linmat%ovrlp%matrix(tmb%orbs%isorb+1:tmb%orbs%isorb+tmb%orbs%norbp,:),&
          tmb%orbs%norbp, tmb%linmat%denskern%matrix(1,1), tmb%orbs%norb, 1.d0, sk, tmb%orbs%norbp)

     ! coeffs and therefore kernel will change, so no need to keep it
     iall=-product(shape(tmb%linmat%denskern%matrix))*kind(tmb%linmat%denskern%matrix)
     deallocate(tmb%linmat%denskern%matrix, stat=istat)
     call memocc(istat, iall, 'tmb%linmat%denskern%matrix', subname)

     allocate(skhp(tmb%orbs%norb,tmb%orbs%norbp), stat=istat)
     call memocc(istat, skhp, 'skhp', subname)

     ! multiply by H to get (I_ab - S_ag K^gb) H_bg, or in this case the transpose of the above
     call dgemm('t', 't', tmb%orbs%norb, tmb%orbs%norbp, tmb%orbs%norb, 1.d0, tmb%linmat%ham%matrix(1,1), &
          tmb%orbs%norb, sk(1,1), tmb%orbs%norbp, 0.d0, skhp(1,1), tmb%orbs%norb)

     iall=-product(shape(sk))*kind(sk)
     deallocate(sk, stat=istat)
     call memocc(istat, iall, 'sk', subname)

     allocate(skh(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
     call memocc(istat, skh, 'skh', subname)

     ! gather together
     if(nproc > 1) then
        call mpi_allgatherv(skhp(1,1), tmb%orbs%norb*tmb%orbs%norbp, mpi_double_precision, skh(1,1), &
           tmb%orbs%norb*tmb%orbs%norb_par(:,0), tmb%orbs%norb*tmb%orbs%isorb_par, &
           mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
     else
        call dcopy(tmb%orbs%norbp*tmb%orbs%norb,skhp(1,1),1,skh(1,1),1)
     end if

     iall=-product(shape(skhp))*kind(skhp)
     deallocate(skhp, stat=istat)
     call memocc(istat, iall, 'skhp', subname)

     allocate(rhs(tmb%orbs%norb,orbs%norb), stat=istat)
     call memocc(istat, rhs, 'rhs', subname)

     ! calc for i on this proc: (I_ab - S_ag K^gb) H_bg c_i^d
     if (orbs%norbp>0) then
        call dgemm('t', 'n', tmb%orbs%norb, orbs%norbp, tmb%orbs%norb, 1.d0, skh(1,1), &
             tmb%orbs%norb, tmb%coeff(1,orbs%isorb+1), tmb%orbs%norb, 0.d0, rhs(1,orbs%isorb+1), tmb%orbs%norb)
     end if

     iall=-product(shape(skh))*kind(skh)
     deallocate(skh, stat=istat)
     call memocc(istat, iall, 'skh', subname)

     ! multiply by f_i to get rhs_i^a
     do iorb=1,orbs%norbp
        iiorb=orbs%isorb+iorb
        rhs(:,iiorb)=rhs(:,iiorb)*orbs%occup(iiorb)
     end do

  !else ! sparse
  !   ! need matrices in compressed form
  !   call compress_matrix_for_allreduce(iproc,tmb%linmat%denskern)
  !   iall=-product(shape(tmb%linmat%denskern%matrix))*kind(tmb%linmat%denskern%matrix)
  !   deallocate(tmb%linmat%denskern%matrix, stat=istat)
  !   call memocc(istat, iall, 'tmb%linmat%denskern%matrix', subname)
  !   call compress_matrix_for_allreduce(iproc,tmb%linmat%ham)
  !   call compress_matrix_for_allreduce(iproc,tmb%linmat%ovrlp)

     !.....

  end if ! sparse/dense

  call timing(iproc,'dirmin_lagmat1','OF')

  ! Solve the linear system ovrlp*grad=rhs
  call timing(iproc,'dirmin_dgesv','ON') !lr408t

  info = 0 ! needed for when some processors have orbs%norbp=0

  allocate(ipiv(tmb%orbs%norb), stat=istat)
  call memocc(istat, ipiv, 'ipiv', subname)

  if(tmb%orthpar%blocksize_pdsyev<0) then
      if (orbs%norbp>0) then
          call dgesv(tmb%orbs%norb, orbs%norbp, tmb%linmat%ovrlp%matrix(1,1), tmb%orbs%norb, ipiv(1), &
               rhs(1,orbs%isorb+1), tmb%orbs%norb, info)
      end if
  else
      call mpiallred(rhs(1,1), tmb%orbs%norb*orbs%norb, mpi_sum, bigdft_mpi%mpi_comm, ierr)
      call dgesv_parallel(iproc, tmb%orthpar%nproc_pdsyev, tmb%orthpar%blocksize_pdsyev, bigdft_mpi%mpi_comm, &
           tmb%orbs%norb, orbs%norb, tmb%linmat%ovrlp%matrix, tmb%orbs%norb, rhs, tmb%orbs%norb, info)
  end if

  if(info/=0) then
      write(*,'(a,i0)') 'ERROR in dgesv: info=',info
      stop
  end if

  iall=-product(shape(ipiv))*kind(ipiv)
  deallocate(ipiv, stat=istat)
  call memocc(istat, iall, 'ipiv', subname)

  call timing(iproc,'dirmin_dgesv','OF') !lr408t

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
      tt=tt+ddot(tmb%orbs%norb, rhs(1,iiorb), 1, rhs(1,iiorb), 1)
  end do
  call mpiallred(tt, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  fnrm=sqrt(tt/dble(orbs%norb))

  iall=-product(shape(rhs))*kind(rhs)
  deallocate(rhs, stat=istat)
  call memocc(istat, iall, 'rhs', subname)

  call timing(iproc,'dirmin_sddiis','OF')

  ! ovrlp%matrix will have been destroyed by dgesv, so need to uncompress again - done in reorthonorm_coeff if needed
  !call uncompressMatrix(iproc,tmb%linmat%ovrlp)
  ! do twice with approx S^_1/2, as not quite good enough at preserving charge if only once, but exact too expensive
  ! instead of twice could add some criterion to check accuracy?
  call reorthonormalize_coeff(iproc, nproc, orbs, -8, -8, 1, tmb%orbs, tmb%linmat%ovrlp, tmb%coeff)
  call reorthonormalize_coeff(iproc, nproc, orbs, -8, -8, 1, tmb%orbs, tmb%linmat%ovrlp, tmb%coeff)

end subroutine optimize_coeffs_sparse





subroutine optimize_coeffs(iproc, nproc, orbs, tmb, ldiis_coeff, fnrm)
  use module_base
  use module_types
  use module_interfaces, except_this_one => optimize_coeffs
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(orbitals_data),intent(in):: orbs
  type(DFT_wavefunction),intent(inout):: tmb
  type(localizedDIISParameters),intent(inout):: ldiis_coeff
  real(8),intent(out):: fnrm

  ! Local variables
  integer:: iorb, jorb, korb, lorb, istat, iall, info, iiorb, ierr, ind, indh, indo, kkorb
  integer :: npts_per_proc, ind_start, ind_end, indc, iseg, segn
  real(8),dimension(:,:),allocatable:: lagmat, rhs, ovrlp_coeff, gradp, coeffp !, ovrlp_tmp
  integer,dimension(:),allocatable:: ipiv
  real(8):: tt, ddot, tt2
  character(len=*),parameter:: subname='optimize_coeffs'

  allocate(lagmat(orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, lagmat, 'lagmat', subname)

  allocate(rhs(tmb%orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, rhs, 'rhs', subname)

  allocate(gradp(tmb%orbs%norb,orbs%norbp), stat=istat)
  call memocc(istat, gradp, 'gradp', subname)

  allocate(ipiv(tmb%orbs%norb), stat=istat)
  call memocc(istat, ipiv, 'ipiv', subname)

  !allocate(ovrlp_tmp(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
  !call memocc(istat, ovrlp_tmp, 'ovrlp_tmp', subname)

  !allocate(ovrlp_coeff(orbs%norb,orbs%norb), stat=istat)
  !call memocc(istat, ovrlp_coeff, 'ovrlp_coeff', subname)

  call compress_matrix_for_allreduce(iproc,tmb%linmat%ham)
  call compress_matrix_for_allreduce(iproc,tmb%linmat%ovrlp)

  call timing(iproc,'dirmin_lagmat1','ON')

  ! Calculate the Lagrange multiplier matrix. Use ovrlp_coeff as temporary array.
  !do iorb=1,orbs%norbp
  !    iiorb=orbs%isorb+iorb
  !    do jorb=1,orbs%norb
  !        tt=0.d0
  !        do korb=1,tmb%orbs%norb
  !            do lorb=1,tmb%orbs%norb
  !                tt=tt+tmb%coeff(korb,jorb)*tmb%coeff(lorb,iiorb)*tmb%linmat%ham%matrix(lorb,korb)
  !            end do
  !        end do
  !        ovrlp_coeff(jorb,iorb)=tt
  !    end do
  !end do

  ! Calculate the Lagrange multiplier matrix. Use ovrlp_coeff as temporary array.
  ! Gather together the complete matrix
  !if (nproc > 1) then
  !   call mpi_allgatherv(ovrlp_coeff(1,1), orbs%norb*orbs%norbp, mpi_double_precision, lagmat(1,1), &
  !        orbs%norb*orbs%norb_par(:,0), orbs%norb*orbs%isorb_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
  !else
  !   call vcopy(orbs%norb*orbs%norb,ovrlp_coeff(1,1),1,lagmat(1,1),1)
  !end if

  call to_zero(orbs%norb*orbs%norb, lagmat(1,1))

  !inelegant way of dividing up
  npts_per_proc = nint(real(tmb%linmat%ham%nvctr + tmb%linmat%ham%full_dim1,dp) / real(nproc*2,dp))
  ind_start = 1+iproc*npts_per_proc
  ind_end = (iproc+1)*npts_per_proc
  if (iproc==nproc-1) ind_end = tmb%linmat%ham%nvctr!ceiling(0.5d0*real(basis_overlap%nvctr + basis_overlap%full_dim1,dp))

  indc=0
  do ind = 1, tmb%linmat%ham%nvctr
     lorb = tmb%linmat%ham%orb_from_index(ind,1)
     kkorb = tmb%linmat%ham%orb_from_index(ind,2)

     if (lorb<kkorb) cycle ! so still only doing half
     indc = indc + 1
     if (indc < ind_start .or. indc > ind_end) cycle

     do iorb=1,orbs%norb
        if (kkorb==lorb) then
           tt=tmb%linmat%ham%matrix_compr(ind)*tmb%coeff(lorb,iorb)
           do jorb=iorb,orbs%norb
              lagmat(jorb,iorb)=lagmat(jorb,iorb)+tmb%coeff(kkorb,jorb)*tt
           end do
        else
           do jorb=iorb,orbs%norb
              lagmat(jorb,iorb)=lagmat(jorb,iorb)&
                   +(tmb%coeff(kkorb,jorb)*tmb%coeff(lorb,iorb)+tmb%coeff(kkorb,iorb)*tmb%coeff(lorb,jorb))&
                   *tmb%linmat%ham%matrix_compr(ind)
           end do
        end if
     end do
  end do

  ! use symmetry to calculate other half
  do iorb=1,orbs%norb
     do jorb=iorb+1,orbs%norb
        lagmat(iorb,jorb) = lagmat(jorb,iorb)
     end do
  end do

  ! All reduce the complete matrix
  if (nproc > 1) then
      call mpiallred(lagmat(1,1), orbs%norb**2, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  end if

  call timing(iproc,'dirmin_lagmat1','OF')
  call timing(iproc,'dirmin_lagmat2','ON')

  ! Calculate the right hand side
  !rhs=0.d0
  !do iorb=1,orbs%norbp
  !    iiorb=orbs%isorb+iorb
  !    do lorb=1,tmb%orbs%norb
  !        tt=0.d0
  !        do korb=1,tmb%orbs%norb
  !            tt=tt+tmb%coeff(korb,iiorb)*tmb%linmat%ham%matrix(korb,lorb)
  !        end do
  !        do jorb=1,orbs%norb
  !            do korb=1,tmb%orbs%norb
  !                tt=tt-lagmat(jorb,iiorb)*tmb%coeff(korb,jorb)*tmb%linmat%ovrlp%matrix(korb,lorb)
  !            end do
  !        end do
  !        rhs(lorb,iiorb)=tt*orbs%occup(iiorb)
  !    end do
  !end do

  ! Calculate the right hand side
  !!rhs=0.d0
  !call to_zero(tmb%orbs%norb*orbs%norbp, rhs(1,orbs%isorb+1))
  !do iorb=1,orbs%norbp
  !    iiorb=orbs%isorb+iorb
  !    do lorb=1,tmb%orbs%norb
  !        tt=0.d0
  !        do korb=1,tmb%orbs%norb
  !            ind=tmb%linmat%ham%matrixindex_in_compressed(korb,lorb)
  !            if (ind==0) cycle
  !            tt=tt+tmb%coeff(korb,iiorb)*tmb%linmat%ham%matrix_compr(ind)
  !        end do
  !        do korb=1,tmb%orbs%norb
  !            ind=tmb%linmat%ovrlp%matrixindex_in_compressed(korb,lorb)
  !            if (ind==0) cycle
  !            do jorb=1,orbs%norb
  !                tt=tt-lagmat(jorb,iiorb)*tmb%coeff(korb,jorb)*tmb%linmat%ovrlp%matrix_compr(ind)
  !            end do
  !        end do
  !        rhs(lorb,iiorb)=tt*orbs%occup(iiorb)
  !    end do
  !end do

  ! option 1 - calculation divided by orbs%norbp, wasteful if some procs have no orbs
  if(tmb%orthpar%blocksize_pdsyev<0 .and. orbs%norbp>0) then
     call to_zero(tmb%orbs%norb*orbs%norbp, rhs(1,orbs%isorb+1))
  else if (tmb%orthpar%blocksize_pdsyev>0) then
     call to_zero(tmb%orbs%norb*orbs%norb, rhs(1,1))
  end if

  if (orbs%norbp>0) then ! don't need to bother if we have no orbs on this proc
     do lorb=1,tmb%orbs%norb
        do korb=lorb,tmb%orbs%norb
           indh=tmb%linmat%ham%matrixindex_in_compressed(korb,lorb)
           if (indh==0) cycle ! H should always be less sparse than S

           do iorb=1,orbs%norbp
              iiorb=orbs%isorb+iorb
              rhs(lorb,iiorb)=rhs(lorb,iiorb)+tmb%coeff(korb,iiorb)*tmb%linmat%ham%matrix_compr(indh)
              if (korb/=lorb) rhs(korb,iiorb)=rhs(korb,iiorb)+tmb%coeff(lorb,iiorb)*tmb%linmat%ham%matrix_compr(indh)
           end do

           indo=tmb%linmat%ovrlp%matrixindex_in_compressed(korb,lorb)
           if (indo==0) cycle

           do jorb=1,orbs%norb
              tt=tmb%coeff(korb,jorb)*tmb%linmat%ovrlp%matrix_compr(indo)
              tt2=tmb%coeff(lorb,jorb)*tmb%linmat%ovrlp%matrix_compr(indo)
              do iorb=1,orbs%norbp
                 iiorb=orbs%isorb+iorb
                 rhs(lorb,iiorb)=rhs(lorb,iiorb)-lagmat(jorb,iiorb)*tt
                 if (lorb/=korb) rhs(korb,iiorb)=rhs(korb,iiorb)-lagmat(jorb,iiorb)*tt2
              end do
           end do

         end do
     end do

     do iorb=1,orbs%norbp
        iiorb=orbs%isorb+iorb
        rhs(:,iiorb)=rhs(:,iiorb)*orbs%occup(iiorb)
     end do
  end if

  ! Solve the linear system ovrlp*grad=rhs
  !call dcopy(tmb%orbs%norb**2, tmb%linmat%ovrlp%matrix(1,1), 1, ovrlp_tmp(1,1), 1)
  call timing(iproc,'dirmin_lagmat2','OF') !lr408t
  call timing(iproc,'dirmin_dgesv','ON') !lr408t

  info = 0 ! needed for when some processors have orbs%orbp=0

  if(tmb%orthpar%blocksize_pdsyev<0) then
      if (orbs%norbp>0) then
          call dgesv(tmb%orbs%norb, orbs%norbp, tmb%linmat%ovrlp%matrix(1,1), tmb%orbs%norb, ipiv(1), &
               rhs(1,orbs%isorb+1), tmb%orbs%norb, info)
      end if
  else
      call mpiallred(rhs(1,1), tmb%orbs%norb*orbs%norb, mpi_sum, bigdft_mpi%mpi_comm, ierr)
      call dgesv_parallel(iproc, tmb%orthpar%nproc_pdsyev, tmb%orthpar%blocksize_pdsyev, bigdft_mpi%mpi_comm, &
           tmb%orbs%norb, orbs%norb, tmb%linmat%ovrlp%matrix, tmb%orbs%norb, rhs, tmb%orbs%norb, info)
  end if

  if(info/=0) then
      write(*,'(a,i0)') 'ERROR in dgesv: info=',info
      stop
  end if

  call dcopy(tmb%orbs%norb*orbs%norbp, rhs(1,orbs%isorb+1), 1, gradp(1,1), 1)
  !if(iproc==0) write(200,*) gradp
  call timing(iproc,'dirmin_dgesv','OF') !lr408t

  ! Precondition the gradient (only making things worse...)
  !call precondition_gradient_coeff(tmb%orbs%norb, tmb%orbs%norbp, tmb%linmat%ham%matrix, tmb%linmat%ovrlp%matrix, gradp)

  call timing(iproc,'dirmin_sddiis','ON')

  ! Improve the coefficients
  if (ldiis_coeff%isx > 0) then
      ldiis_coeff%mis=mod(ldiis_coeff%is,ldiis_coeff%isx)+1
      ldiis_coeff%is=ldiis_coeff%is+1
  end if  

  if (ldiis_coeff%isx > 0) then !do DIIS
     !TO DO: make sure DIIS works
     print *,'in DIIS'
     call DIIS_coeff(iproc, orbs, tmb, gradp, tmb%coeff, ldiis_coeff)
  else  !steepest descent
     allocate(coeffp(tmb%orbs%norb,orbs%norbp),stat=istat)
     call memocc(istat, coeffp, 'coeffp', subname)
     do iorb=1,orbs%norbp
        iiorb = orbs%isorb + iorb
        do jorb=1,tmb%orbs%norb
           coeffp(jorb,iorb)=tmb%coeff(jorb,iiorb)-ldiis_coeff%alpha_coeff*gradp(jorb,iorb)
        end do
     end do

     if(nproc > 1) then 
        call mpi_allgatherv(coeffp(1,1), tmb%orbs%norb*orbs%norbp, mpi_double_precision, tmb%coeff(1,1), &
           tmb%orbs%norb*orbs%norb_par(:,0), tmb%orbs%norb*orbs%isorb_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
     else
        call dcopy(tmb%orbs%norb*orbs%norb,coeffp(1,1),1,tmb%coeff(1,1),1)
     end if
     iall=-product(shape(coeffp))*kind(coeffp)
     deallocate(coeffp, stat=istat)
     call memocc(istat, iall, 'coeffp', subname)
     !if(iproc==0) write(100,*) tmb%coeff
  end if

  !For fnrm, we only sum on the occupied KS orbitals
  tt=0.d0
  do iorb=1,orbs%norbp
      tt=tt+ddot(tmb%orbs%norb, gradp(1,iorb), 1, gradp(1,iorb), 1)
  end do
  call mpiallred(tt, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  fnrm=sqrt(tt/dble(orbs%norb))

  call timing(iproc,'dirmin_sddiis','OF')

  ! ovrlp%matrix will have been destroyed by dgesv, so need to uncompress again - done in reorthonorm_coeff if needed
  !call uncompressMatrix(iproc,tmb%linmat%ovrlp)
  call reorthonormalize_coeff(iproc, nproc, orbs, -8, -8, 0, tmb%orbs, tmb%linmat%ovrlp, tmb%coeff)

  iall=-product(shape(lagmat))*kind(lagmat)
  deallocate(lagmat, stat=istat)
  call memocc(istat, iall, 'lagmat', subname)

  iall=-product(shape(rhs))*kind(rhs)
  deallocate(rhs, stat=istat)
  call memocc(istat, iall, 'rhs', subname)

  iall=-product(shape(gradp))*kind(gradp)
  deallocate(gradp, stat=istat)
  call memocc(istat, iall, 'gradp', subname)

  iall=-product(shape(ipiv))*kind(ipiv)
  deallocate(ipiv, stat=istat)
  call memocc(istat, iall, 'ipiv', subname)

  !iall=-product(shape(ovrlp_tmp))*kind(ovrlp_tmp)
  !deallocate(ovrlp_tmp, stat=istat)
  !call memocc(istat, iall, 'ovrlp_tmp', subname)

  !iall=-product(shape(ovrlp_coeff))*kind(ovrlp_coeff)
  !deallocate(ovrlp_coeff, stat=istat)
  !call memocc(istat, iall, 'ovrlp_coeff', subname)

end subroutine optimize_coeffs

!!Just to test without MPI
!subroutine optimize_coeffs2(iproc, nproc, orbs, tmb, ldiis_coeff, fnrm)
!  use module_base
!  use module_types
!  use module_interfaces, except_this_one => optimize_coeffs
!  implicit none
!
!  ! Calling arguments
!  integer,intent(in):: iproc, nproc
!  type(orbitals_data),intent(in):: orbs
!  type(DFT_wavefunction),intent(inout):: tmb
!  type(localizedDIISParameters),intent(inout):: ldiis_coeff
!  real(8),intent(out):: fnrm
!
!  ! Local variables
!  integer:: iorb, jorb, istat, iall, info
!  integer:: ialpha, ibeta
!  real(8),dimension(:,:),allocatable:: lagmat, rhs, ovrlp_tmp, ovrlp_coeff, grad
!  integer,dimension(:),allocatable:: ipiv
!  real(8):: tt, ddot
!  character(len=*),parameter:: subname='optimize_coeffs2'
!
!  allocate(lagmat(orbs%norb,orbs%norb), stat=istat)
!  call memocc(istat, lagmat, 'lagmat', subname)
!
!  allocate(rhs(tmb%orbs%norb,orbs%norb), stat=istat)
!  call memocc(istat, rhs, 'rhs', subname)
!
!  allocate(grad(tmb%orbs%norb,orbs%norb), stat=istat)
!  call memocc(istat, grad, 'grad', subname)
!
!
!  allocate(ovrlp_tmp(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
!  call memocc(istat, ovrlp_tmp, 'ovrlp_tmp', subname)
!
!  allocate(ovrlp_coeff(orbs%norb,orbs%norb), stat=istat)
!  call memocc(istat, ovrlp_coeff, 'ovrlp_coeff', subname)
!
!
!  call timing(iproc,'dirmin_lagmat1','ON') !lr408t
!
!  ! Calculate the Lagrange multiplier matrix. Use ovrlp_coeff as temporary array.
!  do iorb=1,orbs%norb
!      do jorb=1,orbs%norb
!          tt=0.d0
!          do ialpha=1,tmb%orbs%norb
!              do ibeta=1,tmb%orbs%norb
!                  tt=tt+tmb%coeff(ialpha,jorb)*tmb%coeff(ibeta,iorb)*tmb%linmat%ham%matrix(ialpha,ibeta)
!              end do
!          end do
!          ovrlp_coeff(jorb,iorb)=tt
!      end do
!  end do
!
!  call vcopy(orbs%norb*orbs%norb,ovrlp_coeff(1,1),1,lagmat(1,1),1)
!
!  call timing(iproc,'dirmin_lagmat1','OF')
!  call timing(iproc,'dirmin_lagmat2','ON')
!
!  ! Calculate the right hand side
!  rhs=0.d0
!  do iorb=1,orbs%norb
!      do ialpha=1,tmb%orbs%norb
!          tt=0.d0
!          do ibeta=1,tmb%orbs%norb
!              tt=tt+tmb%coeff(ibeta,iorb)*tmb%linmat%ham%matrix(ibeta,ialpha)
!          end do
!          do jorb=1,orbs%norb
!              do ibeta=1,tmb%orbs%norb
!                  tt=tt-lagmat(jorb,iorb)*tmb%coeff(ibeta,jorb)*tmb%linmat%ovrlp%matrix(ibeta,ialpha)
!              end do
!          end do
!          rhs(ialpha,iorb)=tt*orbs%occup(iorb)
!      end do
!  end do
!
!  ! Solve the linear system ovrlp*grad=rhs
!  call dcopy(tmb%orbs%norb**2, tmb%linmat%ovrlp%matrix(1,1), 1, ovrlp_tmp(1,1), 1)
!
!  call timing(iproc,'dirmin_lagmat2','OF')
!  call timing(iproc,'dirmin_dgesv','ON')
!
!  info = 0 ! needed for when some processors have orbs%orbp=0
!  if(tmb%orthpar%blocksize_pdsyev<0) then
!     allocate(ipiv(tmb%orbs%norb), stat=istat)
!     call memocc(istat, ipiv, 'ipiv', subname)
!     call dgesv(tmb%orbs%norb, orbs%norb, ovrlp_tmp(1,1), tmb%orbs%norb, ipiv(1), &
!          rhs(1,1), tmb%orbs%norb, info)
!     iall=-product(shape(ipiv))*kind(ipiv)
!     deallocate(ipiv, stat=istat)
!     call memocc(istat, iall, 'ipiv', subname)
!  else
!      call dgesv_parallel(iproc, tmb%orthpar%nproc_pdsyev, tmb%orthpar%blocksize_pdsyev, bigdft_mpi%mpi_comm, &
!          tmb%orbs%norb, orbs%norb, ovrlp_tmp, tmb%orbs%norb, rhs, tmb%orbs%norb, info)
!  end if
!
!  if(info/=0) then
!      write(*,'(a,i0)') 'ERROR in dgesv: info=',info
!      stop
!  end if
!
!  call dcopy(tmb%orbs%norb*orbs%norb, rhs(1,1), 1, grad(1,1), 1)
!  call timing(iproc,'dirmin_dgesv','OF') 
!
!  ! Precondition the gradient (only making things worse...)
!  !call precondition_gradient_coeff(tmb%orbs%norb, orbs%norbp, tmb%linmat%ham%matrix, ovrlp, gradp)
!
!  call timing(iproc,'dirmin_sddiis','ON') !lr408t
!
!  ! Improve the coefficients
!  if (ldiis_coeff%isx > 0) then
!      ldiis_coeff%mis=mod(ldiis_coeff%is,ldiis_coeff%isx)+1
!      ldiis_coeff%is=ldiis_coeff%is+1
!  end if  
!
!  if (.false. .and. ldiis_coeff%isx > 0) then !do DIIS, must change this for non parallel
!     call DIIS_coeff(iproc, orbs, tmb, grad, tmb%coeff, ldiis_coeff)
!  else  !steepest descent
!     do iorb=1,orbs%norb
!        do ialpha=1,tmb%orbs%norb
!           tmb%coeff(ialpha,iorb)=tmb%coeff(ialpha,iorb)-ldiis_coeff%alpha_coeff*grad(ialpha,iorb)
!        end do
!     end do
!  end if
!
!  tt=0.d0
!  do iorb=1,orbs%norb
!      print *,'norm gradp',iorb+tmb%orbs%isorb, ddot(tmb%orbs%norb, grad(1,iorb), 1, grad(1,iorb), 1)
!      tt=tt+ddot(tmb%orbs%norb, grad(1,iorb), 1, grad(1,iorb), 1)
!  end do
!  fnrm=sqrt(tt/dble(orbs%norb))
!
!  call timing(iproc,'dirmin_sddiis','OF') !lr408t
!  
!  ! Normalize the coefficients (Loewdin)
!  call reorthonormalize_coeff(iproc, nproc, orbs%norb, -8, -8, 0, tmb%orbs, tmb%linmat%ovrlp%matrix, tmb%coeff)
!
!  iall=-product(shape(lagmat))*kind(lagmat)
!  deallocate(lagmat, stat=istat)
!  call memocc(istat, iall, 'lagmat', subname)
!
!  iall=-product(shape(rhs))*kind(rhs)
!  deallocate(rhs, stat=istat)
!  call memocc(istat, iall, 'rhs', subname)
!
!  iall=-product(shape(grad))*kind(grad)
!  deallocate(grad, stat=istat)
!  call memocc(istat, iall, 'grad', subname)
!
!  iall=-product(shape(ovrlp_tmp))*kind(ovrlp_tmp)
!  deallocate(ovrlp_tmp, stat=istat)
!  call memocc(istat, iall, 'ovrlp_tmp', subname)
!
!  iall=-product(shape(ovrlp_coeff))*kind(ovrlp_coeff)
!  deallocate(ovrlp_coeff, stat=istat)
!  call memocc(istat, iall, 'ovrlp_coeff', subname)
!
!end subroutine optimize_coeffs2

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
  integer:: iorb, ii, istat
  character(len=*),parameter:: subname='allocate_DIIS_coeff'

  allocate(ldiis%mat(ldiis%isx,ldiis%isx,tmb%orbs%norbp),stat=istat)
  call memocc(istat, ldiis%mat, 'ldiis%mat', subname)

  ii=ldiis%isx*tmb%orbs%norb*tmb%orbs%norbp
  allocate(ldiis%phiHist(ii), stat=istat)
  call memocc(istat, ldiis%phiHist, 'ldiis%phiHist', subname)
  allocate(ldiis%hphiHist(ii), stat=istat)
  call memocc(istat, ldiis%hphiHist, 'ldiis%hphiHist', subname)

end subroutine allocate_DIIS_coeff

