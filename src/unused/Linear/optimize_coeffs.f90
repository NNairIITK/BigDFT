subroutine transform_coeffs_to_derivatives(iproc, nproc, orbs, lzd, tmb, tmbder)
  use module_base
  use module_types
  use module_interfaces, except_this_one => transform_coeffs_to_derivatives
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(orbitals_data),intent(in):: orbs
  type(local_zone_descriptors),intent(in):: lzd
  type(DFT_wavefunction),intent(in):: tmb
  type(DFT_wavefunction),intent(inout):: tmbder

  ! Local variables
  integer:: iorb, jorb, korb, istat, iall, info, kkorb
  real(8):: tt, ddot
  real(8),dimension(:),allocatable:: psit_c, psit_f
  real(8),dimension(:,:,:),allocatable:: ovrlp
  real(8),dimension(:,:),allocatable:: coeff_tmp, ovrlp_coeff
  integer,dimension(:),allocatable:: ipiv
  character(len=*),parameter:: subname='transform_coeffs_to_derivatives'

  allocate(ovrlp(tmbder%orbs%norb,tmbder%orbs%norb,2), stat=istat)
  call memocc(istat, ovrlp, 'ovrlp', subname)

  allocate(ipiv(tmbder%orbs%norb), stat=istat)
  call memocc(istat, ipiv, 'ipiv', subname)

  allocate(coeff_tmp(tmbder%orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, coeff_tmp, 'coeff_tmp', subname)

  allocate(ovrlp_coeff(orbs%norb,orbs%norb), stat=istat)
  call memocc(istat, ovrlp_coeff, 'ovrlp_coeff', subname)

  ! Calculate overlap matrix for derivatives
  allocate(psit_c(tmbder%collcom%ndimind_c))
  call memocc(istat, psit_c, 'psit_c', subname)
  allocate(psit_f(7*tmbder%collcom%ndimind_f))
  call memocc(istat, psit_f, 'psit_f', subname)
  call transpose_localized(iproc, nproc, tmbder%orbs, tmbder%collcom, tmbder%psi, psit_c, psit_f, lzd)
  call calculate_overlap_transposed(iproc, nproc, tmbder%orbs, tmbder%mad, tmbder%collcom, psit_c, &
       psit_c, psit_f, psit_f, ovrlp(1,1,1))
  call untranspose_localized(iproc, nproc, tmbder%orbs, tmbder%collcom, psit_c, psit_f, tmbder%psi, lzd)
  iall=-product(shape(psit_c))*kind(psit_c)
  deallocate(psit_c, stat=istat)
  call memocc(istat, iall, 'psit_c', subname)
  iall=-product(shape(psit_f))*kind(psit_f)
  deallocate(psit_f, stat=istat)
  call memocc(istat, iall, 'psit_f', subname)

  ! Calculate right hand side
  do iorb=1,orbs%norb
      do jorb=1,tmbder%orbs%norb
          tt=0.d0
          kkorb=0
          do korb=1,tmbder%orbs%norb,4
              kkorb=kkorb+1
              tt=tt+ovrlp(korb,jorb,1)*tmb%wfnmd%coeff(kkorb,iorb)
          end do
          tmbder%wfnmd%coeff(jorb,iorb)=tt
      end do
  end do


  call dcopy(tmbder%orbs%norb**2, ovrlp(1,1,1), 1, ovrlp(1,1,2), 1)
  call dgesv(tmbder%orbs%norb, orbs%norb, ovrlp(1,1,2), tmbder%orbs%norb, ipiv(1), tmbder%wfnmd%coeff(1,1), tmbder%orbs%norb, info)
  if(info/=0) then
      write(*,'(a,i0)') 'ERROR in dgsesv (transform_coeffs_to_derivatives): info=',info
      stop
  end if

  !!if(iproc==0) then
  !!    do iorb=1,orbs%norb
  !!        do jorb=1,tmbder%orbs%norb
  !!            write(200,'(2i8,es14.6)') iorb, jorb, tmbder%wfnmd%coeff(jorb,iorb)
  !!        end do
  !!    end do
  !!end if

  ! Normalize the coeffiecients.
  ! Loewdin
  call dgemm('n', 'n', tmbder%orbs%norb, orbs%norb, tmbder%orbs%norb, 1.d0, ovrlp(1,1,1), tmbder%orbs%norb, &
       tmbder%wfnmd%coeff(1,1), tmbder%orbs%norb, 0.d0, coeff_tmp(1,1), tmbder%orbs%norb)
  !!if(iproc==0) then
  !!    do iorb=1,orbs%norb
  !!        do jorb=1,tmbder%orbs%norb
  !!            write(200,'(2i8,2es14.6)') iorb, jorb, tmbder%wfnmd%coeff(jorb,iorb), coeff_tmp(jorb,iorb)
  !!        end do
  !!    end do
  !!end if
  do iorb=1,orbs%norb
      do jorb=1,orbs%norb
          ovrlp_coeff(jorb,iorb)=ddot(tmbder%orbs%norb, tmbder%wfnmd%coeff(1,jorb), 1, coeff_tmp(1,iorb), 1)
          !!if(iproc==0) write(400,'(a,2i8,es15.6)') 'iorb, jorb, ovrlp_coeff(jorb,iorb)', iorb, jorb, ovrlp_coeff(jorb,iorb)
      end do
  end do
  ! WARNING: this is the wrong mad, but it does not matter for iorder=0
  call overlapPowerMinusOneHalf(iproc, nproc, mpi_comm_world, 0, -8, -8, orbs%norb, ovrlp_coeff)

  call dgemm('n', 'n', tmbder%orbs%norb, orbs%norb, orbs%norb, 1.d0, tmbder%wfnmd%coeff(1,1), tmbder%orbs%norb, &
       ovrlp_coeff(1,1), orbs%norb, 0.d0, coeff_tmp(1,1), tmbder%orbs%norb)
  call dcopy(tmbder%orbs%norb*orbs%norb, coeff_tmp(1,1), 1, tmbder%wfnmd%coeff(1,1), 1)

  !!if(iproc==0) then
  !!    do iorb=1,orbs%norb
  !!        do jorb=1,tmbder%orbs%norb
  !!            write(210,'(2i8,es14.6)') iorb, jorb, tmbder%wfnmd%coeff(jorb,iorb)
  !!        end do
  !!    end do
  !!end if
  !!call mpi_barrier(mpi_comm_world,istat)


  iall=-product(shape(ovrlp))*kind(ovrlp)
  deallocate(ovrlp, stat=istat)
  call memocc(istat, iall, 'ovrlp', subname)

  iall=-product(shape(ipiv))*kind(ipiv)
  deallocate(ipiv, stat=istat)
  call memocc(istat, iall, 'ipiv', subname)

  iall=-product(shape(coeff_tmp))*kind(coeff_tmp)
  deallocate(coeff_tmp, stat=istat)
  call memocc(istat, iall, 'coeff_tmp', subname)

  iall=-product(shape(ovrlp_coeff))*kind(ovrlp_coeff)
  deallocate(ovrlp_coeff, stat=istat)
  call memocc(istat, iall, 'ovrlp_coeff', subname)

end subroutine transform_coeffs_to_derivatives


