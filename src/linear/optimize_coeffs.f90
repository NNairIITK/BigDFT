subroutine optimize_coeffs(iproc, nproc, orbs, ham, ovrlp, tmb, fnrm)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(orbitals_data),intent(in):: orbs
  type(DFT_wavefunction),intent(inout):: tmb
  real(8),dimension(tmb%orbs%norb,tmb%orbs%norb),intent(in):: ham, ovrlp
  real(8),intent(out):: fnrm

  ! Local variables
  integer:: iorb, jorb, korb, lorb, istat, iall, info
  real(8),dimension(:,:),allocatable:: lagmat, rhs, grad, ovrlp_tmp, coeff_tmp, ovrlp_coeff
  integer,dimension(:),allocatable:: ipiv
  real(8):: tt, ddot, tt2, tt3, mean_alpha, dnrm2
  character(len=*),parameter:: subname='optimize_coeffs'

  allocate(lagmat(orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, lagmat, 'lagmat', subname)

  allocate(rhs(tmb%orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, rhs, 'rhs', subname)

  allocate(grad(tmb%orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, grad, 'grad', subname)

  allocate(ipiv(tmb%orbs%norb), stat=istat)
  call memocc(istat, ipiv, 'ipiv', subname)

  allocate(ovrlp_tmp(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
  call memocc(istat, ovrlp_tmp, 'ovrlp_tmp', subname)

  allocate(coeff_tmp(tmb%orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, coeff_tmp, 'coeff_tmp', subname)

  allocate(ovrlp_coeff(orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, ovrlp_coeff, 'ovrlp_coeff', subname)


  ! Check normalization
  call dgemm('n', 'n', tmb%orbs%norb, orbs%norb, tmb%orbs%norb, 1.d0, ovrlp(1,1), tmb%orbs%norb, &
       tmb%wfnmd%coeff(1,1), tmb%orbs%norb, 0.d0, coeff_tmp(1,1), tmb%orbs%norb)
  do iorb=1,orbs%norb
      do jorb=1,orbs%norb
          tt=ddot(tmb%orbs%norb, tmb%wfnmd%coeff(1,iorb), 1, coeff_tmp(1,jorb), 1)
          tt2=ddot(tmb%orbs%norb, coeff_tmp(1,iorb), 1, tmb%wfnmd%coeff(1,jorb), 1)
          tt3=ddot(tmb%orbs%norb, tmb%wfnmd%coeff(1,iorb), 1, tmb%wfnmd%coeff(1,jorb), 1)
          if(iproc==0) write(100,'(2i6,3es15.5)') iorb, jorb, tt, tt2, tt3
      end do
  end do

  ! Calculate the Lagrange multiplier matrix
  do iorb=1,orbs%norb
      do jorb=1,orbs%norb
          tt=0.d0
          do korb=1,tmb%orbs%norb
              do lorb=1,tmb%orbs%norb
                  tt=tt+tmb%wfnmd%coeff(korb,jorb)*tmb%wfnmd%coeff(lorb,iorb)*ham(lorb,korb)
              end do
          end do
          lagmat(jorb,iorb)=tt
          if(iproc==0) write(510,*) iorb, jorb, lagmat(jorb,iorb)
      end do
  end do

  ! Calculate the right hand side
  do iorb=1,orbs%norb
      do lorb=1,tmb%orbs%norb
          tt=0.d0
          do korb=1,tmb%orbs%norb
              tt=tt+tmb%wfnmd%coeff(korb,iorb)*ham(korb,lorb)
          end do
          do jorb=1,orbs%norb
              do korb=1,tmb%orbs%norb
                  tt=tt-lagmat(jorb,iorb)*tmb%wfnmd%coeff(korb,jorb)*ovrlp(korb,lorb)
              end do
          end do
          rhs(lorb,iorb)=tt
          if(iproc==0) write(520,*) iorb, lorb, rhs(lorb,iorb)
      end do
  end do

  ! Solve the linear system ovrlp*grad=rhs
  call dcopy(tmb%orbs%norb**2, ovrlp(1,1), 1, ovrlp_tmp(1,1), 1)
  call dgesv(tmb%orbs%norb, orbs%norb, ovrlp_tmp(1,1), tmb%orbs%norb, ipiv(1), rhs(1,1), tmb%orbs%norb, info)
  if(info/=0) then
      write(*,'(a,i0)') 'ERROR in dgesv: info=',info
      stop
  end if
  call dcopy(tmb%orbs%norb*orbs%norb, rhs(1,1), 1, grad(1,1), 1)


  ! Improve the coefficients
  tt=0.d0
  do iorb=1,orbs%norb
      do jorb=1,tmb%orbs%norb
          if(iproc==0) write(500,'(a,2i8,2es14.6)') 'iorb, jorb, tmb%wfnmd%coeff(jorb,iorb), grad(jorb,iorb)', iorb, jorb, tmb%wfnmd%coeff(jorb,iorb), grad(jorb,iorb)
          tmb%wfnmd%coeff(jorb,iorb)=tmb%wfnmd%coeff(jorb,iorb)-tmb%wfnmd%alpha_coeff(iorb)*grad(jorb,iorb)
      end do
      tt=tt+ddot(tmb%orbs%norb, grad(1,iorb), 1, grad(1,iorb), 1)
  end do
  tt=sqrt(tt)
  fnrm=tt
  if(iproc==0) write(*,'(a,es13.5)') 'coeff gradient: ',tt
  tmb%wfnmd%it_coeff_opt=tmb%wfnmd%it_coeff_opt+1
  if(tmb%wfnmd%it_coeff_opt>1) then
      mean_alpha=0.d0
      do iorb=1,orbs%norb
          tt=ddot(tmb%orbs%norb, grad(1,iorb), 1, tmb%wfnmd%grad_coeff_old(1,iorb), 1)
          tt=tt/(dnrm2(tmb%orbs%norb, grad(1,iorb), 1)*dnrm2(tmb%orbs%norb, tmb%wfnmd%grad_coeff_old(1,iorb), 1))
          if(iproc==0) write(*,*) 'iorb, tt', iorb, tt
          if(tt>.8d0) then
              tmb%wfnmd%alpha_coeff(iorb)=1.2d0*tmb%wfnmd%alpha_coeff(iorb)
          else
              tmb%wfnmd%alpha_coeff(iorb)=0.5d0*tmb%wfnmd%alpha_coeff(iorb)
          end if
          mean_alpha=mean_alpha+tmb%wfnmd%alpha_coeff(iorb)
      end do
      mean_alpha=mean_alpha/dble(orbs%norb)
      if(iproc==0) write(*,*) 'mean_alpha',mean_alpha
  end if
  call dcopy(tmb%orbs%norb*orbs%norb, grad(1,1), 1, tmb%wfnmd%grad_coeff_old(1,1), 1)

  ! Normalize the coeffiecients.
  ! Loewdin
  !call random_number(tmb%wfnmd%coeff)
  call dgemm('n', 'n', tmb%orbs%norb, orbs%norb, tmb%orbs%norb, 1.d0, ovrlp(1,1), tmb%orbs%norb, &
       tmb%wfnmd%coeff(1,1), tmb%orbs%norb, 0.d0, coeff_tmp(1,1), tmb%orbs%norb)
  do iorb=1,orbs%norb
      do jorb=1,orbs%norb
          ovrlp_coeff(jorb,iorb)=ddot(tmb%orbs%norb, tmb%wfnmd%coeff(1,jorb), 1, coeff_tmp(1,iorb), 1)
          if(iproc==0) write(400,'(a,2i8,es15.6)') 'iorb, jorb, ovrlp_coeff(jorb,iorb)', iorb, jorb, ovrlp_coeff(jorb,iorb)
      end do
  end do
  ! WARNING: this is the wrong mad, but it does not matter for iorder=0
  call overlapPowerMinusOneHalf(iproc, nproc, mpi_comm_world, 0, -8, -8, orbs%norb, tmb%mad, ovrlp_coeff)

  call dgemm('n', 'n', tmb%orbs%norb, orbs%norb, orbs%norb, 1.d0, tmb%wfnmd%coeff(1,1), tmb%orbs%norb, &
       ovrlp_coeff(1,1), orbs%norb, 0.d0, coeff_tmp(1,1), tmb%orbs%norb)
  call dcopy(tmb%orbs%norb*orbs%norb, coeff_tmp(1,1), 1, tmb%wfnmd%coeff(1,1), 1)

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

  ! Check normalization
  call dgemm('n', 'n', tmb%orbs%norb, orbs%norb, tmb%orbs%norb, 1.d0, ovrlp(1,1), tmb%orbs%norb, &
       tmb%wfnmd%coeff(1,1), tmb%orbs%norb, 0.d0, coeff_tmp(1,1), tmb%orbs%norb)
  do iorb=1,orbs%norb
      do jorb=1,orbs%norb
          tt=ddot(tmb%orbs%norb, tmb%wfnmd%coeff(1,iorb), 1, coeff_tmp(1,jorb), 1)
          tt2=ddot(tmb%orbs%norb, coeff_tmp(1,iorb), 1, tmb%wfnmd%coeff(1,jorb), 1)
          tt3=ddot(tmb%orbs%norb, tmb%wfnmd%coeff(1,iorb), 1, tmb%wfnmd%coeff(1,jorb), 1)
          if(iproc==0) write(200,'(2i6,3es15.5)') iorb, jorb, tt, tt2, tt3
      end do
  end do


  iall=-product(shape(lagmat))*kind(lagmat)
  deallocate(lagmat, stat=istat)
  call memocc(istat, iall, 'lagmat', subname)

  iall=-product(shape(rhs))*kind(rhs)
  deallocate(rhs, stat=istat)
  call memocc(istat, iall, 'rhs', subname)

  iall=-product(shape(grad))*kind(grad)
  deallocate(grad, stat=istat)
  call memocc(istat, iall, 'grad', subname)

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
