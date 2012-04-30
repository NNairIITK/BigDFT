subroutine optimize_coeffs(iproc, nproc, orbs, ham, ovrlp, tmb)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(orbitals_data),intent(in):: orbs
  type(DFT_wavefunction),intent(inout):: tmb
  real(8),dimension(tmb%orbs%norb,tmb%orbs%norb),intent(in):: ham, ovrlp

  ! Local variables
  integer:: iorb, jorb, korb, lorb, istat, iall, info
  real(8),dimension(:,:),allocatable:: lagmat, rhs, grad, ovrlp_tmp, coeff_tmp, ovrlp_coeff
  integer,dimension(:),allocatable:: ipiv
  real(8):: tt, alpha, ddot, tt2, tt3
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
  call memocc(istat, ovrlp_tmp, 'ovrlp_tmp', subname)

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
      end do
  end do
  write(*,*) '1'

  ! Calculate the right hand side
  do iorb=1,orbs%norb
      do lorb=1,tmb%orbs%norb
          tt=0.d0
          do korb=1,tmb%orbs%norb
              tt=tt+tmb%wfnmd%coeff(korb,iorb)*ham(korb,lorb)
          end do
          do jorb=1,orbs%norb
              do korb=1,tmb%orbs%norb
                  tt=tt+lagmat(jorb,iorb)*tmb%wfnmd%coeff(korb,jorb)*ovrlp(korb,lorb)
              end do
          end do
          rhs(lorb,iorb)=tt
      end do
  end do
  write(*,*) '2'

  ! Solve the linear system ovrlp*grad=rhs
  call dcopy(tmb%orbs%norb**2, ovrlp(1,1), 1, ovrlp_tmp(1,1), 1)
  call dgesv(tmb%orbs%norb, orbs%norb, ovrlp_tmp(1,1), tmb%orbs%norb, ipiv(1), rhs(1,1), tmb%orbs%norb, info)
  if(info/=0) then
      write(*,'(a,i0)') 'ERROR in dgesv: info=',info
      stop
  end if
  write(*,*) '3'


  ! Improve the coefficients
  alpha=0.d-1
  do iorb=1,orbs%norb
      do jorb=1,tmb%orbs%norb
          tmb%wfnmd%coeff(jorb,iorb)=tmb%wfnmd%coeff(jorb,iorb)-alpha*grad(jorb,iorb)
      end do
  end do
  write(*,*) '4'

  !!!!! Normalize the coeffiecients.
  !!!!!!call dgemm('n', 'n', tmb%orbs%norb, orbs%norb, tmb%orbs%norb, 1.d0, ovrlp(1,1), tmb%orbs%norb, &
  !!!!!!     tmb%wfnmd%coeff(1,1), tmb%orbs%norb, 0.d0, coeff_tmp(1,1), tmb%orbs%norb)
  !!!!do iorb=1,orbs%norb
  !!!!    do jorb=1,orbs%norb
  !!!!        ovrlp_coeff(jorb,iorb)=ddot(tmb%orbs%norb, tmb%wfnmd%coeff(1,jorb),1 , tmb%wfnmd%coeff(1,iorb), 1)
  !!!!        !!ovrlp_coeff(jorb,iorb)=ddot(tmb%orbs%norb, tmb%wfnmd%coeff(1,jorb),1 , coeff_tmp(1,iorb), 1)
  !!!!    end do
  !!!!end do
  !!!!! WARNING: this is the wrong mad, but it does not matter for iorder=0
  !!!!call overlapPowerMinusOneHalf(iproc, nproc, mpi_comm_world, 0, -8, -8, orbs%norb, tmb%mad, ovrlp_coeff)

  !!!!call dcopy(tmb%orbs%norb**2, ovrlp(1,1), 1, ovrlp_tmp(1,1), 1)
  !!!!call overlapPowerMinusOneHalf(iproc, nproc, mpi_comm_world, 0, -8, -8, tmb%orbs%norb, tmb%mad, ovrlp_tmp)
  !!!!! WARNING: this is the wrong mad, but it does not matter for iorder=0
  !!!!!call overlapPowerMinusOne(iproc, nproc, 0, tmb%orbs%norb, tmb%mad, tmb%orbs, ovrlp_tmp)
  !!!!call dgemm('n', 'n', tmb%orbs%norb, orbs%norb, orbs%norb, 1.d0, ovrlp_tmp(1,1), tmb%orbs%norb, &
  !!!!     tmb%wfnmd%coeff(1,1), tmb%orbs%norb, 0.d0, coeff_tmp(1,1), tmb%orbs%norb)

  !!!!write(*,*) '4.2'
  !!!!!!call dgemm('n', 'n', tmb%orbs%norb, orbs%norb, orbs%norb, 1.d0, tmb%wfnmd%coeff(1,1), tmb%orbs%norb, &
  !!!!!!     ovrlp_coeff(1,1), orbs%norb, 0.d0, coeff_tmp(1,1), tmb%orbs%norb)
  !!!!call dgemm('n', 'n', tmb%orbs%norb, orbs%norb, orbs%norb, 1.d0, coeff_tmp(1,1), tmb%orbs%norb, &
  !!!!     ovrlp_coeff(1,1), orbs%norb, 0.d0, tmb%wfnmd%coeff(1,1), tmb%orbs%norb)
  !!!!!!call dgemm('n', 'n', tmb%orbs%norb, orbs%norb, tmb%orbs%norb, 1.d0, ovrlp_tmp(1,1), tmb%orbs%norb, &
  !!!!!!     coeff_tmp(1,1), tmb%orbs%norb, 0.d0, tmb%wfnmd%coeff(1,1), tmb%orbs%norb)
  !!!!!!call dcopy(tmb%orbs%norb*orbs%norb, coeff_tmp(1,1), 1, tmb%wfnmd%coeff(1,1), 1)
  !!!!write(*,*) '4.3'
  !!!!write(*,*) '5'


  ! gram schmidt
  if(iproc==0) write(300,*) coeff_tmp
  if(iproc==0) write(400,*) tmb%wfnmd%coeff
  do iorb=1,orbs%norb
      do jorb=1,iorb-1
          call dgemv('n', tmb%orbs%norb, tmb%orbs%norb, 1.d0, ovrlp(1,1), &
               tmb%orbs%norb, tmb%wfnmd%coeff(1,jorb), 1, 0.d0, coeff_tmp(1,jorb), 1)
          if(iproc==0) write(310,*) coeff_tmp(:,jorb)
          if(iproc==0) write(310,'(a,2i8,es14.5)') 'iorb, jorb, tt', iorb, jorb, tt
          tt=ddot(tmb%orbs%norb, tmb%wfnmd%coeff(1,iorb), 1, coeff_tmp(1,jorb), 1)
          if(iproc==0) write(*,'(a,2i8,es14.5)') 'iorb, jorb, tt', iorb, jorb, tt
          call daxpy(tmb%orbs%norb, -tt, tmb%wfnmd%coeff(1,jorb), 1, tmb%wfnmd%coeff(1,iorb), 1)
      end do
      call dgemv('n', tmb%orbs%norb, tmb%orbs%norb, 1.d0, ovrlp(1,1), &
           tmb%orbs%norb, tmb%wfnmd%coeff(1,iorb), 1, 0.d0, coeff_tmp(1,iorb), 1)
      if(iproc==0) write(320,*) tmb%wfnmd%coeff(:,iorb),'//', coeff_tmp(:,iorb)
      tt=ddot(tmb%orbs%norb, tmb%wfnmd%coeff(1,iorb), 1, coeff_tmp(1,iorb), 1)
      if(iproc==0) write(*,'(a,i8,es14.5)') 'iorb, tt', iorb, tt
      call dscal(tmb%orbs%norb, 1/sqrt(tt), tmb%wfnmd%coeff(1,iorb), 1)
  end do

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
  write(*,*) '6'


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
