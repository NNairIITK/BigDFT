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


!!subroutine distribute_coefficients(orbs, tmb)
!!  use module_base
!!  use module_types
!!  implicit none
!!
!!  ! Calling arguments
!!  type(orbitals_data),intent(in):: orbs
!!  type(DFT_wavefunction),intent(inout):: tmb
!!
!!  ! Local variables
!!  integer:: iorb, iiorb
!!
!!  do iorb=1,orbs%norbp
!!      iiorb=orbs%isorb+iorb
!!      call dcopy(tmb%orbs%norb, tmb%wfnmd%coeff(1,iiorb), 1, tmb%wfnmd%coeffp(1,iorb), 1)
!!  end do
!!
!!end subroutine distribute_coefficients
!!
!!
!!
!!subroutine collect_coefficients(nproc, orbs, tmb, coeffp, coeff)
!!  use module_base
!!  use module_types
!!  implicit none
!!
!!  ! Calling arguments
!!  integer, intent(in) :: nproc
!!  type(orbitals_data),intent(in):: orbs
!!  type(DFT_wavefunction),intent(inout):: tmb
!!  real(8),dimension(tmb%orbs%norb,orbs%norbp),intent(in):: coeffp
!!  real(8),dimension(tmb%orbs%norb,orbs%norb),intent(out):: coeff
!!
!!  ! Local variables
!!  integer:: ierr
!!
!!  if (nproc > 1) then
!!     call mpi_allgatherv(coeffp(1,1), tmb%orbs%norb*orbs%norbp, mpi_double_precision, coeff(1,1), &
!!          tmb%orbs%norb*orbs%norb_par(:,0), tmb%orbs%norb*orbs%isorb_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
!!  else
!!     call vcopy(tmb%orbs%norb*orbs%norb,coeffp(1,1),1,coeff(1,1),1)
!!  end if
!!
!!end subroutine collect_coefficients



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
      call vcopy(ncount, coeff(ist+tmb%orbs%isorb*tmb%orbs%norb), 1, ldiis%phiHist(jst), 1)
      call vcopy(ncount, grad(ist), 1, ldiis%hphiHist(jst), 1)
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
      call to_zero(ncount, coeff(ist))
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

