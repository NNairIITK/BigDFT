!> @file
!! Optimize the coefficients
!! @author
!!    Copyright (C) 2011-2012 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


subroutine optimize_coeffs(iproc, nproc, orbs, ham, ovrlp, tmb, ldiis_coeff, fnrm)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(orbitals_data),intent(in):: orbs
  type(DFT_wavefunction),intent(inout):: tmb
  real(8),dimension(tmb%orbs%norb,tmb%orbs%norb),intent(in):: ham, ovrlp
  type(localizedDIISParameters),intent(inout):: ldiis_coeff
  real(8),intent(out):: fnrm

  ! Local variables
  integer:: iorb, jorb, korb, lorb, istat, iall, info, iiorb, ierr
  real(8),dimension(:,:),allocatable:: lagmat, rhs, ovrlp_tmp, coeff_tmp, ovrlp_coeff, gradp
  integer,dimension(:),allocatable:: ipiv
  real(8):: tt, ddot, mean_alpha, dnrm2
  character(len=*),parameter:: subname='optimize_coeffs'

  allocate(lagmat(orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, lagmat, 'lagmat', subname)

  allocate(rhs(tmb%orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, rhs, 'rhs', subname)

  !!allocate(grad(tmb%orbs%norb,orbs%norb), stat=istat)
  !!call memocc(istat, grad, 'grad', subname)

  allocate(gradp(tmb%orbs%norb,orbs%norbp), stat=istat)
  call memocc(istat, gradp, 'gradp', subname)

  allocate(ipiv(tmb%orbs%norb), stat=istat)
  call memocc(istat, ipiv, 'ipiv', subname)

  allocate(ovrlp_tmp(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
  call memocc(istat, ovrlp_tmp, 'ovrlp_tmp', subname)

  allocate(coeff_tmp(tmb%orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, coeff_tmp, 'coeff_tmp', subname)

  allocate(ovrlp_coeff(orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, ovrlp_coeff, 'ovrlp_coeff', subname)


  !!$! Check normalization
  !!$call dgemm('n', 'n', tmb%orbs%norb, orbs%norb, tmb%orbs%norb, 1.d0, ovrlp(1,1), tmb%orbs%norb, &
  !!$     tmb%wfnmd%coeff(1,1), tmb%orbs%norb, 0.d0, coeff_tmp(1,1), tmb%orbs%norb)
  !!$do iorb=1,orbs%norb
  !!$    do jorb=1,orbs%norb
  !!$        tt=ddot(tmb%orbs%norb, tmb%wfnmd%coeff(1,iorb), 1, coeff_tmp(1,jorb), 1)
  !!$        tt2=ddot(tmb%orbs%norb, coeff_tmp(1,iorb), 1, tmb%wfnmd%coeff(1,jorb), 1)
  !!$        tt3=ddot(tmb%orbs%norb, tmb%wfnmd%coeff(1,iorb), 1, tmb%wfnmd%coeff(1,jorb), 1)
  !!$        !!if(iproc==0) write(100,'(2i6,3es15.5)') iorb, jorb, tt, tt2, tt3
  !!$    end do
  !!$end do

  call distribute_coefficients(orbs, tmb)

  ! Calculate the Lagrange multiplier matrix. Use ovrlp_coeff as temporary array.
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      do jorb=1,orbs%norb
          tt=0.d0
          do korb=1,tmb%orbs%norb
              do lorb=1,tmb%orbs%norb
                  tt=tt+tmb%wfnmd%coeff(korb,jorb)*tmb%wfnmd%coeffp(lorb,iorb)*ham(lorb,korb)
              end do
          end do
          !lagmat(jorb,iorb)=tt
          ovrlp_coeff(jorb,iorb)=tt
          !!if(iproc==0) write(510,*) iorb, jorb, lagmat(jorb,iorb)
      end do
  end do

  ! Gather together the complete matrix
  call mpi_allgatherv(ovrlp_coeff(1,1), orbs%norb*orbs%norbp, mpi_double_precision, lagmat(1,1), &
       orbs%norb*orbs%norb_par(:,0), orbs%norb*orbs%isorb_par, mpi_double_precision, mpi_comm_world, ierr)


  ! ##############################################################################
  ! ################################ OLD #########################################
  ! Calculate the right hand side
  !!do iorb=1,orbs%norb
  rhs=0.d0
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      do lorb=1,tmb%orbs%norb
          tt=0.d0
          do korb=1,tmb%orbs%norb
              tt=tt+tmb%wfnmd%coeffp(korb,iorb)*ham(korb,lorb)
          end do
          do jorb=1,orbs%norb
              do korb=1,tmb%orbs%norb
                  tt=tt-lagmat(jorb,iiorb)*tmb%wfnmd%coeff(korb,jorb)*ovrlp(korb,lorb)
              end do
          end do
          !rhs(lorb,iorb)=tt
          rhs(lorb,iiorb)=tt
      end do
  end do
  call mpiallred(rhs(1,1), orbs%norb*tmb%orbs%norb, mpi_sum, mpi_comm_world, ierr)

  ! Solve the linear system ovrlp*grad=rhs
  call dcopy(tmb%orbs%norb**2, ovrlp(1,1), 1, ovrlp_tmp(1,1), 1)

  if(tmb%wfnmd%bpo%blocksize_pdsyev<0) then
      if (orbs%norbp>0) then
          call dgesv(tmb%orbs%norb, orbs%norbp, ovrlp_tmp(1,1), tmb%orbs%norb, ipiv(1), &
               rhs(1,orbs%isorb+1), tmb%orbs%norb, info)
      end if
  else
      call dgesv_parallel(iproc, tmb%wfnmd%bpo%nproc_pdsyev, tmb%wfnmd%bpo%blocksize_pdsyev, mpi_comm_world, &
           tmb%orbs%norb, orbs%norb, ovrlp_tmp, tmb%orbs%norb, rhs, tmb%orbs%norb, info)
  end if

  if(info/=0) then
      write(*,'(a,i0)') 'ERROR in dgesv: info=',info
      stop
  end if
  !call dcopy(tmb%orbs%norb*orbs%norb, rhs(1,1), 1, grad(1,1), 1)

  if(tmb%wfnmd%bpo%blocksize_pdsyev<0) then
      call dcopy(tmb%orbs%norb*orbs%norbp, rhs(1,orbs%isorb+1), 1, gradp(1,1), 1)
  else
      call dcopy(tmb%orbs%norb*orbs%norbp, rhs(1,orbs%isorb+1), 1, gradp(1,1), 1)
  end if

  ! ##############################################################################
  ! ############################ END OLD #########################################


  !! NEW VERSION - TEST ######################################################
  !do iorb=1,orbs%norb
  !    do jorb=1,tmb%orbs%norb
  !         tt=0.d0
  !         do korb=1,tmb%orbs%norb
  !             tt=tt+ham(korb,jorb)*tmb%wfnmd%coeff(korb,iorb)
  !         end do
  !         do korb=1,orbs%norb
  !             do lorb=1,tmb%orbs%norb
  !                 tt=tt-lagmat(korb,iorb)*ovrlp(lorb,jorb)*tmb%wfnmd%coeff(lorb,korb)
  !             end do
  !         end do
  !         grad(jorb,iorb)=tt
  !    end do
  !end do

  !!! #############################################################################
  !!! ########################## NEW ##############################################
  !!do iorb=1,orbs%norbp
  !!    iiorb=orbs%isorb+iorb
  !!    do jorb=1,tmb%orbs%norb
  !!         tt=0.d0
  !!         do korb=1,tmb%orbs%norb
  !!             tt=tt+ham(korb,jorb)*tmb%wfnmd%coeffp(korb,iorb)
  !!             !write(100+iproc,'(4i8,es14.5)') iorb, iiorb, jorb, korb, tmb%wfnmd%coeffp(korb,iorb)
  !!         end do
  !!         do korb=1,orbs%norb
  !!             do lorb=1,tmb%orbs%norb
  !!                 tt=tt-lagmat(korb,iiorb)*ovrlp(lorb,jorb)*tmb%wfnmd%coeff(lorb,korb)
  !!                 !write(100+iproc,'(6i8,3es14.5)') iorb, iiorb, jorb, lorb, korb, korb, tmb%wfnmd%coeff(lorb,korb), ovrlp(lorb,jorb), lagmat(korb,iiorb)
  !!             end do
  !!         end do
  !!         gradp(jorb,iorb)=tt
  !!    end do
  !!end do
  !!! #############################################################################
  !!! ###################### END NEW ##############################################

  !! #########################################################################, 
  !!do iorb=1,orbs%norbp
  !!    iiorb=orbs%isorb+iorb
  !!    do jorb=1,tmb%orbs%norb
  !!        write(300+iproc,'(3i8,es14.6)') iorb,iiorb,jorb,gradp(jorb,iorb)
  !!    end do
  !!end do
  !!do iorb=1,orbs%norb
  !!    do jorb=1,tmb%orbs%norb
  !!        write(200+iproc,*) iorb,jorb,grad(jorb,iorb)
  !!    end do
  !!end do

  ! Precondition the gradient
  call precondition_gradient_coeff(tmb%orbs%norb, orbs%norbp, ham, ovrlp, gradp)


  ! Improve the coefficients
  if (ldiis_coeff%isx > 0) then
      ldiis_coeff%mis=mod(ldiis_coeff%is,ldiis_coeff%isx)+1
      ldiis_coeff%is=ldiis_coeff%is+1
  end if  
  !!do iorb=1,orbs%norb
  !!    call dscal(tmb%orbs%norb, tmb%wfnmd%alpha_coeff(iorb), grad(1,iorb), 1)
  !!end do

  call DIIS_coeff(iproc, nproc, orbs, tmb, gradp, tmb%wfnmd%coeffp, ldiis_coeff)


  tt=0.d0
  do iorb=1,orbs%norbp
      do jorb=1,tmb%orbs%norb
          !if(iproc==0) write(500,'(a,2i8,2es14.6)') 'iorb, jorb, tmb%wfnmd%coeff(jorb,iorb), grad(jorb,iorb)', iorb, jorb, tmb%wfnmd%coeff(jorb,iorb), grad(jorb,iorb)
          !tmb%wfnmd%coeff(jorb,iorb)=tmb%wfnmd%coeff(jorb,iorb)-tmb%wfnmd%alpha_coeff(iorb)*grad(jorb,iorb)
      end do
      tt=tt+ddot(tmb%orbs%norb, gradp(1,iorb), 1, gradp(1,iorb), 1)
  end do
  call mpiallred(tt, 1, mpi_sum, mpi_comm_world, ierr)
  fnrm=sqrt(tt/dble(orbs%norb))
  !if(iproc==0) write(*,'(a,es13.5)') 'coeff gradient: ',tt
  tmb%wfnmd%it_coeff_opt=tmb%wfnmd%it_coeff_opt+1
  if(tmb%wfnmd%it_coeff_opt>1) then
      mean_alpha=0.d0
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          tt=ddot(tmb%orbs%norb, gradp(1,iorb), 1, tmb%wfnmd%grad_coeff_old(1,iorb), 1)
          tt=tt/(dnrm2(tmb%orbs%norb, gradp(1,iorb), 1)*dnrm2(tmb%orbs%norb, tmb%wfnmd%grad_coeff_old(1,iorb), 1))
          !if(iproc==0) write(*,*) 'iorb, tt', iorb, tt
          if(tt>.85d0) then
              tmb%wfnmd%alpha_coeff(iiorb)=1.1d0*tmb%wfnmd%alpha_coeff(iiorb)
          else
              tmb%wfnmd%alpha_coeff(iiorb)=0.5d0*tmb%wfnmd%alpha_coeff(iiorb)
          end if
          mean_alpha=mean_alpha+tmb%wfnmd%alpha_coeff(iiorb)
      end do
      mean_alpha=mean_alpha/dble(orbs%norb)
      call mpiallred(mean_alpha, 1, mpi_sum, mpi_comm_world, ierr)
      if(iproc==0) write(*,*) 'mean_alpha',mean_alpha
  end if
  call collect_coefficients(orbs, tmb, tmb%wfnmd%coeffp, tmb%wfnmd%coeff)
  !call collect_coefficients(orbs, tmb, gradp, grad)

  call dcopy(tmb%orbs%norb*orbs%norbp, gradp(1,1), 1, tmb%wfnmd%grad_coeff_old(1,1), 1)


  ! Normalize the coeffiecients (Loewdin)

  ! Calculate the overlap matrix among the coefficients with resct to ovrlp. Use lagmat as temporary array.
  call dgemm('n', 'n', tmb%orbs%norb, orbs%norbp, tmb%orbs%norb, 1.d0, ovrlp(1,1), tmb%orbs%norb, &
       tmb%wfnmd%coeffp(1,1), tmb%orbs%norb, 0.d0, coeff_tmp(1,1), tmb%orbs%norb)
  do iorb=1,orbs%norbp
      do jorb=1,orbs%norb
          lagmat(jorb,iorb)=ddot(tmb%orbs%norb, tmb%wfnmd%coeff(1,jorb), 1, coeff_tmp(1,iorb), 1)
      end do
  end do
  ! Gather together the complete matrix
  call mpi_allgatherv(lagmat(1,1), orbs%norb*orbs%norbp, mpi_double_precision, ovrlp_coeff(1,1), &
       orbs%norb*orbs%norb_par(:,0), orbs%norb*orbs%isorb_par, mpi_double_precision, mpi_comm_world, ierr)

  ! WARNING: this is the wrong mad, but it does not matter for iorder=0
  call overlapPowerMinusOneHalf(iproc, nproc, mpi_comm_world, 0, -8, -8, orbs%norb, ovrlp_coeff)

  ! Build the new linear combinations
  call dgemm('n', 'n', tmb%orbs%norb, orbs%norbp, orbs%norb, 1.d0, tmb%wfnmd%coeff(1,1), tmb%orbs%norb, &
       ovrlp_coeff(1,orbs%isorb+1), orbs%norb, 0.d0, coeff_tmp(1,1), tmb%orbs%norb)
  ! Gather together the results partial results.
  call collect_coefficients(orbs, tmb, coeff_tmp(1,1), tmb%wfnmd%coeff)

  !!! Gram schmidt
  !!do iorb=1,orbs%norb
  !!    do jorb=1,iorb-1
  !!        call dgemv('n', tmb%orbs%norb, tmb%orbs%norb, 1.d0, ovrlp(1,1), &
  !!             tmb%orbs%norb, tmb%wfnmd%coeff(1,jorb), 1, 0.d0, coeff_tmp(1,jorb), 1)
  !!        tt=ddot(tmb%orbs%norb, tmb%wfnmd%coeff(1,iorb), 1, coeff_tmp(1,jorb), 1)
  !!        call daxpy(tmb%orbs%norb, -tt, tmb%wfnmd%coeff(1,jorb), 1, tmb%wfnmd%coeff(1,iorb), 1)
  !!    end do
  !!    call dgemv('n', tmb%orbs%norb, tmb%orbs%norb, 1.d0, ovrlp(1,1), &
  !!         tmb%orbs%norb, tmb%wfnmd%coeff(1,iorb), 1, 0.d0, coeff_tmp(1,iorb), 1)
  !!    tt=ddot(tmb%orbs%norb, tmb%wfnmd%coeff(1,iorb), 1, coeff_tmp(1,iorb), 1)
  !!    call dscal(tmb%orbs%norb, 1/sqrt(tt), tmb%wfnmd%coeff(1,iorb), 1)
  !!end do

  !!! Check normalization
  !!call dgemm('n', 'n', tmb%orbs%norb, orbs%norb, tmb%orbs%norb, 1.d0, ovrlp(1,1), tmb%orbs%norb, &
  !!     tmb%wfnmd%coeff(1,1), tmb%orbs%norb, 0.d0, coeff_tmp(1,1), tmb%orbs%norb)
  !!do iorb=1,orbs%norb
  !!    do jorb=1,orbs%norb
  !!        tt=ddot(tmb%orbs%norb, tmb%wfnmd%coeff(1,iorb), 1, coeff_tmp(1,jorb), 1)
  !!        tt2=ddot(tmb%orbs%norb, coeff_tmp(1,iorb), 1, tmb%wfnmd%coeff(1,jorb), 1)
  !!        tt3=ddot(tmb%orbs%norb, tmb%wfnmd%coeff(1,iorb), 1, tmb%wfnmd%coeff(1,jorb), 1)
  !!        if(iproc==0) write(200,'(2i6,3es15.5)') iorb, jorb, tt, tt2, tt3
  !!    end do
  !!end do


  iall=-product(shape(lagmat))*kind(lagmat)
  deallocate(lagmat, stat=istat)
  call memocc(istat, iall, 'lagmat', subname)

  iall=-product(shape(rhs))*kind(rhs)
  deallocate(rhs, stat=istat)
  call memocc(istat, iall, 'rhs', subname)

  !!iall=-product(shape(grad))*kind(grad)
  !!deallocate(grad, stat=istat)
  !!call memocc(istat, iall, 'grad', subname)

  iall=-product(shape(gradp))*kind(gradp)
  deallocate(gradp, stat=istat)
  call memocc(istat, iall, 'gradp', subname)

  iall=-product(shape(ipiv))*kind(ipiv)
  deallocate(ipiv, stat=istat)
  call memocc(istat, iall, 'ipiv', subname)

  iall=-product(shape(ovrlp_tmp))*kind(ovrlp_tmp)
  deallocate(ovrlp_tmp, stat=istat)
  call memocc(istat, iall, 'ovrlp_tmp', subname)

  iall=-product(shape(coeff_tmp))*kind(coeff_tmp)
  deallocate(coeff_tmp, stat=istat)
  call memocc(istat, iall, 'coeff_tmp', subname)

  iall=-product(shape(ovrlp_coeff))*kind(ovrlp_coeff)
  deallocate(ovrlp_coeff, stat=istat)
  call memocc(istat, iall, 'ovrlp_coeff', subname)

end subroutine optimize_coeffs


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



subroutine DIIS_coeff(iproc, nproc, orbs, tmb, grad, coeff, ldiis)
  use module_base
  use module_types
  use module_interfaces, except_this_one => DIIS_coeff
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(orbitals_data),intent(in):: orbs
  type(DFT_wavefunction),intent(in):: tmb
  real(8),dimension(tmb%orbs%norb*orbs%norbp),intent(in):: grad
  real(8),dimension(tmb%orbs%norb*orbs%norbp),intent(inout):: coeff
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
  lwork=100*ldiis%isx
  allocate(work(lwork), stat=istat)
  call memocc(istat, work, 'work', subname)
  allocate(ipiv(ldiis%isx+1), stat=istat)
  call memocc(istat, ipiv, 'ipiv', subname)
  
  mat=0.d0
  rhs=0.d0
  call to_zero((ldiis%isx+1)**2, mat(1,1))
  call to_zero(ldiis%isx+1, rhs(1))
  
  ! Copy coeff and grad to history.
  ist=1
  do iorb=1,orbs%norbp
      jst=1
      do jorb=1,iorb-1
          ncount=tmb%orbs%norb
          jst=jst+ncount*ldiis%isx
      end do
      ncount=tmb%orbs%norb
      jst=jst+(ldiis%mis-1)*ncount
      call dcopy(ncount, coeff(ist), 1, ldiis%phiHist(jst), 1)
      call dcopy(ncount, grad(ist), 1, ldiis%hphiHist(jst), 1)
      ist=ist+ncount
  end do
  
  do iorb=1,orbs%norbp
      ! Shift the DIIS matrix left up if we reached the maximal history length.
      if(ldiis%is>ldiis%isx) then
         do i=1,ldiis%isx-1
            do j=1,i
               ldiis%mat(j,i,iorb)=ldiis%mat(j+1,i+1,iorb)
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
          ncount=tmb%orbs%norb
          jst=jst+ncount*ldiis%isx
          ist1=ist1+ncount
      end do
      ncount=tmb%orbs%norb
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
  
  
  ist=1
  do iorb=1,orbs%norbp
      
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
      do istat=1,ldiis%isx+1
          do iall=1,ldiis%isx+1
              !!if(iproc==0) write(500,*) istat, iall, mat(iall,istat)
          end do
      end do
      if(ldiis%is>1) then
         call dsysv('u', min(ldiis%isx,ldiis%is)+1, 1, mat, ldiis%isx+1,  & 
              ipiv, rhs(1), ldiis%isx+1, work, lwork, info)
         
         if (info /= 0) then
            write(*,'(a,i0)') 'ERROR in dsysv (DIIS_coeff), info=', info
            stop
         end if
      else
         rhs(1)=1.d0
      endif
  
  
      ! Make a new guess for the orbital.
      ncount=tmb%orbs%norb
      call razero(ncount, coeff(ist))
      isthist=max(1,ldiis%is-ldiis%isx+1)
      jj=0
      jst=0
      do jorb=1,iorb-1
          ncount=tmb%orbs%norb
          jst=jst+ncount*ldiis%isx
      end do
      do j=isthist,ldiis%is
          jj=jj+1
          mj=mod(j-1,ldiis%isx)+1
          ncount=tmb%orbs%norb
          jjst=jst+(mj-1)*ncount
          do k=1,ncount
              coeff(ist+k-1) = coeff(ist+k-1) + rhs(jj)*(ldiis%phiHist(jjst+k)-ldiis%hphiHist(jjst+k))
          end do
      end do
  
      ncount=tmb%orbs%norb
      ist=ist+ncount
  end do
  
  
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
  
  !!call timing(iproc,'optimize_DIIS ','OF')


end subroutine DIIS_coeff



subroutine initialize_DIIS_coeff(isx, ldiis)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: isx
  type(localizedDIISParameters),intent(out):: ldiis
  
  ! Local variables
  character(len=*),parameter:: subname='initialize_DIIS_coeff'
  
  
  ldiis%isx=isx
  ldiis%is=0
  ldiis%switchSD=.false.
  ldiis%trmin=1.d100
  ldiis%trold=1.d100

end subroutine initialize_DIIS_coeff


subroutine allocate_DIIS_coeff(tmb, orbs, ldiis)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  type(DFT_wavefunction),intent(in):: tmb
  type(orbitals_data),intent(in):: orbs
  type(localizedDIISParameters),intent(out):: ldiis
  
  ! Local variables
  integer:: iorb, ii, istat
  character(len=*),parameter:: subname='allocate_DIIS_coeff'

  allocate(ldiis%mat(ldiis%isx,ldiis%isx,orbs%norb),stat=istat)
  call memocc(istat, ldiis%mat, 'ldiis%mat', subname)
  ii=0
  do iorb=1,orbs%norb
      ii=ii+ldiis%isx*tmb%orbs%norb
  end do
  allocate(ldiis%phiHist(ii), stat=istat)
  call memocc(istat, ldiis%phiHist, 'ldiis%phiHist', subname)
  allocate(ldiis%hphiHist(ii), stat=istat)
  call memocc(istat, ldiis%hphiHist, 'ldiis%hphiHist', subname)

end subroutine allocate_DIIS_coeff





subroutine distribute_coefficients(orbs, tmb)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(orbitals_data),intent(in):: orbs
  type(DFT_wavefunction),intent(inout):: tmb

  ! Local variables
  integer:: iorb, iiorb

  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      call dcopy(tmb%orbs%norb, tmb%wfnmd%coeff(1,iiorb), 1, tmb%wfnmd%coeffp(1,iorb), 1)
  end do

end subroutine distribute_coefficients



subroutine collect_coefficients(orbs, tmb, coeffp, coeff)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(orbitals_data),intent(in):: orbs
  type(DFT_wavefunction),intent(in):: tmb
  real(8),dimension(tmb%orbs%norb,orbs%norbp),intent(in):: coeffp
  real(8),dimension(tmb%orbs%norb,orbs%norb),intent(out):: coeff

  ! Local variables
  integer:: ierr

  call mpi_allgatherv(coeffp(1,1), tmb%orbs%norb*orbs%norbp, mpi_double_precision, coeff(1,1), &
       tmb%orbs%norb*orbs%norb_par(:,0), tmb%orbs%norb*orbs%isorb_par, mpi_double_precision, mpi_comm_world, ierr)

end subroutine collect_coefficients
