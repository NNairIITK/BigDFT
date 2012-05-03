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
          !!if(iproc==0) write(510,*) iorb, jorb, lagmat(jorb,iorb)
      end do
  end do

  !!! Calculate the right hand side
  !!do iorb=1,orbs%norb
  !!    do lorb=1,tmb%orbs%norb
  !!        tt=0.d0
  !!        do korb=1,tmb%orbs%norb
  !!            tt=tt+tmb%wfnmd%coeff(korb,iorb)*ham(korb,lorb)
  !!        end do
  !!        do jorb=1,orbs%norb
  !!            do korb=1,tmb%orbs%norb
  !!                tt=tt-lagmat(jorb,iorb)*tmb%wfnmd%coeff(korb,jorb)*ovrlp(korb,lorb)
  !!            end do
  !!        end do
  !!        rhs(lorb,iorb)=tt
  !!        !!if(iproc==0) write(520,*) iorb, lorb, rhs(lorb,iorb)
  !!    end do
  !!end do

  !!! Solve the linear system ovrlp*grad=rhs
  !!call dcopy(tmb%orbs%norb**2, ovrlp(1,1), 1, ovrlp_tmp(1,1), 1)
  !!call dgesv(tmb%orbs%norb, orbs%norb, ovrlp_tmp(1,1), tmb%orbs%norb, ipiv(1), rhs(1,1), tmb%orbs%norb, info)
  !!if(info/=0) then
  !!    write(*,'(a,i0)') 'ERROR in dgesv: info=',info
  !!    stop
  !!end if
  !!call dcopy(tmb%orbs%norb*orbs%norb, rhs(1,1), 1, grad(1,1), 1)

  !! NEW VERSION - TEST ######################################################
  do iorb=1,orbs%norb
      do jorb=1,tmb%orbs%norb
           tt=0.d0
           do korb=1,tmb%orbs%norb
               tt=tt+ham(korb,jorb)*tmb%wfnmd%coeff(korb,iorb)
           end do
           do korb=1,orbs%norb
               do lorb=1,tmb%orbs%norb
                   tt=tt-lagmat(korb,iorb)*ovrlp(lorb,jorb)*tmb%wfnmd%coeff(lorb,korb)
               end do
           end do
           grad(jorb,iorb)=tt
      end do
  end do
  !! #########################################################################

  ! Precondition the gradient
  call precondition_gradient_coeff(tmb%orbs%norb, orbs%norb, ham, ovrlp, grad)

  ! Improve the coefficients
  if (ldiis_coeff%isx > 0) then
      ldiis_coeff%mis=mod(ldiis_coeff%is,ldiis_coeff%isx)+1                                                                               
      ldiis_coeff%is=ldiis_coeff%is+1                                                                                               
  end if  
  !!do iorb=1,orbs%norb
  !!    call dscal(tmb%orbs%norb, tmb%wfnmd%alpha_coeff(iorb), grad(1,iorb), 1)
  !!end do
  call DIIS_coeff(iproc, nproc, orbs, tmb, grad, tmb%wfnmd%coeff, ldiis_coeff)
  tt=0.d0
  do iorb=1,orbs%norb
      do jorb=1,tmb%orbs%norb
          !if(iproc==0) write(500,'(a,2i8,2es14.6)') 'iorb, jorb, tmb%wfnmd%coeff(jorb,iorb), grad(jorb,iorb)', iorb, jorb, tmb%wfnmd%coeff(jorb,iorb), grad(jorb,iorb)
          !tmb%wfnmd%coeff(jorb,iorb)=tmb%wfnmd%coeff(jorb,iorb)-tmb%wfnmd%alpha_coeff(iorb)*grad(jorb,iorb)
      end do
      tt=tt+ddot(tmb%orbs%norb, grad(1,iorb), 1, grad(1,iorb), 1)
  end do
  tt=sqrt(tt)
  fnrm=tt
  !if(iproc==0) write(*,'(a,es13.5)') 'coeff gradient: ',tt
  tmb%wfnmd%it_coeff_opt=tmb%wfnmd%it_coeff_opt+1
  if(tmb%wfnmd%it_coeff_opt>1) then
      mean_alpha=0.d0
      do iorb=1,orbs%norb
          tt=ddot(tmb%orbs%norb, grad(1,iorb), 1, tmb%wfnmd%grad_coeff_old(1,iorb), 1)
          tt=tt/(dnrm2(tmb%orbs%norb, grad(1,iorb), 1)*dnrm2(tmb%orbs%norb, tmb%wfnmd%grad_coeff_old(1,iorb), 1))
          !if(iproc==0) write(*,*) 'iorb, tt', iorb, tt
          if(tt>.85d0) then
              tmb%wfnmd%alpha_coeff(iorb)=1.1d0*tmb%wfnmd%alpha_coeff(iorb)
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
          !if(iproc==0) write(400,'(a,2i8,es15.6)') 'iorb, jorb, ovrlp_coeff(jorb,iorb)', iorb, jorb, ovrlp_coeff(jorb,iorb)
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
real(8),dimension(tmb%orbs%norb*orbs%norb),intent(in):: grad
real(8),dimension(tmb%orbs%norb*orbs%norb),intent(inout):: coeff
type(localizedDIISParameters),intent(inout):: ldiis

! Local variables
integer:: iorb, jorb, ist, ilr, ncount, jst, i, j, mi, ist1, ist2, jlr, istat, lwork, info
integer:: mj, jj, k, jjst, isthist, ierr, iall
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

! Copy coeff and grad to history.
ist=1
do iorb=1,orbs%norb
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

do iorb=1,orbs%norb
    ! Shift the DIIS matrix left up if we reached the maximal history length.
    if(ldiis%is>ldiis%isx) then
       do i=1,ldiis%isx-1
          do j=1,i
             ldiis%mat(j,i,iorb)=ldiis%mat(j+1,i+1,iorb)
          end do
       end do
    end if
end do



do iorb=1,orbs%norb

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
do iorb=1,orbs%norb
    
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
            if(iproc==0) write(500,*) istat, iall, mat(iall,istat)
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



subroutine initialize_DIIS_coeff(isx, tmb, orbs, ldiis)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: isx
type(DFT_wavefunction),intent(in):: tmb
type(orbitals_data),intent(in):: orbs
type(localizedDIISParameters),intent(out):: ldiis

! Local variables
integer:: iorb, ii, istat, ilr
character(len=*),parameter:: subname='initialize_DIIS_coeff'


ldiis%isx=isx
ldiis%is=0
ldiis%switchSD=.false.
ldiis%trmin=1.d100
ldiis%trold=1.d100
allocate(ldiis%mat(ldiis%isx,ldiis%isx,orbs%norb), stat=istat)
call memocc(istat, ldiis%mat, 'ldiis%mat', subname)
ii=0
do iorb=1,orbs%norb
    ii=ii+ldiis%isx*tmb%orbs%norb
end do
allocate(ldiis%phiHist(ii), stat=istat)
call memocc(istat, ldiis%phiHist, 'ldiis%phiHist', subname)
allocate(ldiis%hphiHist(ii), stat=istat)
call memocc(istat, ldiis%hphiHist, 'ldiis%hphiHist', subname)


end subroutine initialize_DIIS_coeff




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

  if(iproc==0) then
      do iorb=1,orbs%norb
          do jorb=1,tmbder%orbs%norb
              write(200,'(2i8,es14.6)') iorb, jorb, tmbder%wfnmd%coeff(jorb,iorb)
          end do
      end do
  end if

  ! Normalize the coeffiecients.
  ! Loewdin
  call dgemm('n', 'n', tmbder%orbs%norb, orbs%norb, tmbder%orbs%norb, 1.d0, ovrlp(1,1,1), tmbder%orbs%norb, &
       tmbder%wfnmd%coeff(1,1), tmbder%orbs%norb, 0.d0, coeff_tmp(1,1), tmbder%orbs%norb)
  if(iproc==0) then
      do iorb=1,orbs%norb
          do jorb=1,tmbder%orbs%norb
              write(200,'(2i8,2es14.6)') iorb, jorb, tmbder%wfnmd%coeff(jorb,iorb), coeff_tmp(jorb,iorb)
          end do
      end do
  end if
  do iorb=1,orbs%norb
      do jorb=1,orbs%norb
          ovrlp_coeff(jorb,iorb)=ddot(tmbder%orbs%norb, tmbder%wfnmd%coeff(1,jorb), 1, coeff_tmp(1,iorb), 1)
          if(iproc==0) write(400,'(a,2i8,es15.6)') 'iorb, jorb, ovrlp_coeff(jorb,iorb)', iorb, jorb, ovrlp_coeff(jorb,iorb)
      end do
  end do
  ! WARNING: this is the wrong mad, but it does not matter for iorder=0
  call overlapPowerMinusOneHalf(iproc, nproc, mpi_comm_world, 0, -8, -8, orbs%norb, tmbder%mad, ovrlp_coeff)

  call dgemm('n', 'n', tmbder%orbs%norb, orbs%norb, orbs%norb, 1.d0, tmbder%wfnmd%coeff(1,1), tmbder%orbs%norb, &
       ovrlp_coeff(1,1), orbs%norb, 0.d0, coeff_tmp(1,1), tmbder%orbs%norb)
  call dcopy(tmbder%orbs%norb*orbs%norb, coeff_tmp(1,1), 1, tmbder%wfnmd%coeff(1,1), 1)

  if(iproc==0) then
      do iorb=1,orbs%norb
          do jorb=1,tmbder%orbs%norb
              write(210,'(2i8,es14.6)') iorb, jorb, tmbder%wfnmd%coeff(jorb,iorb)
          end do
      end do
  end if
  call mpi_barrier(mpi_comm_world,istat)


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
