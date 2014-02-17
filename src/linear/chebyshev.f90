!> @file
!!  Linear version: Define Chebyshev polynomials
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

 
!> Again assuming all matrices have same sparsity, still some tidying to be done
subroutine chebyshev_clean(iproc, nproc, npl, cc, orbs, foe_obj, sparsemat, kernel, ham_compr, &
           ovrlp_compr, calculate_SHS, nsize_polynomial, SHS, fermi, penalty_ev, chebyshev_polynomials)
  use module_base
  use module_types
  use module_interfaces, except_this_one => chebyshev_clean
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, npl, nsize_polynomial
  real(8),dimension(npl,3),intent(in) :: cc
  type(orbitals_data),intent(in) :: orbs
  type(foe_data),intent(in) :: foe_obj
  type(sparseMatrix), intent(in) :: sparsemat, kernel
  real(kind=8),dimension(sparsemat%nvctr),intent(in) :: ham_compr, ovrlp_compr
  logical,intent(in) :: calculate_SHS
  real(kind=8),dimension(sparsemat%nvctr),intent(inout) :: SHS
  real(kind=8),dimension(orbs%norb,orbs%norbp),intent(out) :: fermi
  real(kind=8),dimension(orbs%norb,orbs%norbp,2),intent(out) :: penalty_ev
  real(kind=8),dimension(nsize_polynomial,npl),intent(out) :: chebyshev_polynomials
  ! Local variables
  integer :: istat, iorb,iiorb, jorb, iall,ipl,norb,norbp,isorb, ierr, nseq, nmaxsegk, nmaxvalk
  integer :: isegstart, isegend, iseg, ii, jjorb, nout
  character(len=*),parameter :: subname='chebyshev_clean'
  real(8), dimension(:,:,:), allocatable :: vectors
  real(kind=8),dimension(:),allocatable :: ham_compr_seq, ovrlp_compr_seq, SHS_seq
  real(kind=8),dimension(:,:),allocatable :: matrix
  real(kind=8) :: tt
  integer,dimension(:,:,:),allocatable :: istindexarr
  integer,dimension(:),allocatable :: ivectorindex
  integer,parameter :: one=1, three=3
  integer,parameter :: number_of_matmuls=one
  integer,dimension(:,:),pointer :: onedimindices

  call timing(iproc, 'chebyshev_comp', 'ON')

  norb = orbs%norb
  norbp = orbs%norbp
  isorb = orbs%isorb

  if (norbp>0) then

      call init_onedimindices(norbp, isorb, foe_obj, sparsemat, nout, onedimindices)
    
      call determine_sequential_length(norbp, isorb, norb, foe_obj, sparsemat, nseq, nmaxsegk, nmaxvalk)
    
    
      ham_compr_seq = f_malloc(nseq,id='ham_compr_seq')
      ovrlp_compr_seq = f_malloc(nseq,id='ovrlp_compr_seq')
      istindexarr = f_malloc((/ nmaxvalk, nmaxsegk, norbp /),id='istindexarr')
      ivectorindex = f_malloc(nseq,id='ivectorindex')
    
      call get_arrays_for_sequential_acces(norbp, isorb, norb, foe_obj, sparsemat, nseq, nmaxsegk, nmaxvalk, &
           istindexarr, ivectorindex)
    
    
      if (number_of_matmuls==one) then
          matrix = f_malloc((/ orbs%norb, orbs%norbp /),id='matrix')
          SHS_seq = f_malloc(nseq,id='SHS_seq')
    
          if (norbp>0) then
              call to_zero(norb*norbp, matrix(1,1))
          end if
          !write(*,*) 'WARNING CHEBYSHEV: MODIFYING MATRIX MULTIPLICATION'
          if (orbs%norbp>0) then
              isegstart=sparsemat%istsegline(orbs%isorb_par(iproc)+1)
              if (orbs%isorb+orbs%norbp<orbs%norb) then
                  isegend=sparsemat%istsegline(orbs%isorb_par(iproc+1)+1)-1
              else
                  isegend=sparsemat%nseg
              end if
              do iseg=isegstart,isegend
                  ii=sparsemat%keyv(iseg)-1
                  do jorb=sparsemat%keyg(1,iseg),sparsemat%keyg(2,iseg)
                      ii=ii+1
                      iiorb = (jorb-1)/orbs%norb + 1
                      jjorb = jorb - (iiorb-1)*orbs%norb
                      matrix(jjorb,iiorb-orbs%isorb)=ovrlp_compr(ii)
                      !if (jjorb==iiorb) then
                      !    matrix(jjorb,iiorb-orbs%isorb)=1.d0
                      !else
                      !    matrix(jjorb,iiorb-orbs%isorb)=0.d0
                      !end if
                  end do
              end do
          end if
      end if
    
      call sequential_acces_matrix(norbp, isorb, norb, foe_obj, sparsemat, ham_compr, nseq, nmaxsegk, nmaxvalk, &
           ham_compr_seq)
    
    
      call sequential_acces_matrix(norbp, isorb, norb, foe_obj, sparsemat, ovrlp_compr, nseq, nmaxsegk, nmaxvalk, &
           ovrlp_compr_seq)
    
      vectors = f_malloc((/ norb, norbp, 4 /),id='vectors')
      if (norbp>0) then
          call to_zero(norb*norbp, vectors(1,1,1))
      end if
    
  end if
    
  if (number_of_matmuls==one) then
  
      if (calculate_SHS) then
  
          if (norbp>0) then
              call sparsemm(nseq, ham_compr_seq, matrix(1,1), vectors(1,1,1), &
                   norb, norbp, ivectorindex, nout, onedimindices)
              call to_zero(norbp*norb, matrix(1,1))
              call sparsemm(nseq, ovrlp_compr_seq, vectors(1,1,1), matrix(1,1), &
                   norb, norbp, ivectorindex, nout, onedimindices)
              !call to_zero(sparsemat%nvctr, SHS(1))
          end if
          call to_zero(sparsemat%nvctr, SHS(1))
          
          if (orbs%norbp>0) then
              isegstart=sparsemat%istsegline(orbs%isorb_par(iproc)+1)
              if (orbs%isorb+orbs%norbp<orbs%norb) then
                  isegend=sparsemat%istsegline(orbs%isorb_par(iproc+1)+1)-1
              else
                  isegend=sparsemat%nseg
              end if
              do iseg=isegstart,isegend
                  ii=sparsemat%keyv(iseg)-1
                  do jorb=sparsemat%keyg(1,iseg),sparsemat%keyg(2,iseg)
                      ii=ii+1
                      iiorb = (jorb-1)/orbs%norb + 1
                      jjorb = jorb - (iiorb-1)*orbs%norb
                      SHS(ii)=matrix(jjorb,iiorb-orbs%isorb)
                  end do
              end do
          end if
  
          call mpiallred(SHS(1), sparsemat%nvctr, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  
      end if
  
      if (orbs%norbp>0) then
          call sequential_acces_matrix(norbp, isorb, norb, foe_obj, sparsemat, SHS, nseq, nmaxsegk, &
               nmaxvalk, SHS_seq)
      end if
  
  end if

  !!if (iproc==0) then
  !!    do istat=1,sparsemat%nvctr
  !!        write(300,*) ham_compr(istat), SHS(istat)
  !!    end do
  !!end if
    
  if (norbp>0) then
    
      ! No need to set to zero the 3rd and 4th entry since they will be overwritten
      ! by copies of the 1st entry.
      if (norbp>0) then
          call to_zero(2*norb*norbp, vectors(1,1,1))
      end if
      do iorb=1,norbp
          iiorb=isorb+iorb
          vectors(iiorb,iorb,1)=1.d0
      end do
    
      if (norbp>0) then
    
          call vcopy(norb*norbp, vectors(1,1,1), 1, vectors(1,1,3), 1)
          call vcopy(norb*norbp, vectors(1,1,1), 1, vectors(1,1,4), 1)
        
          ! apply(3/2 - 1/2 S) H (3/2 - 1/2 S)
          if (number_of_matmuls==three) then
              call sparsemm(nseq, ovrlp_compr_seq, vectors(1,1,3), vectors(1,1,1), &
                   norb, norbp, ivectorindex, nout, onedimindices)
              call sparsemm(nseq, ham_compr_seq, vectors(1,1,1), vectors(1,1,3), &
                   norb, norbp, ivectorindex, nout, onedimindices)
              call sparsemm(nseq, ovrlp_compr_seq, vectors(1,1,3), vectors(1,1,1), &
                   norb, norbp, ivectorindex, nout, onedimindices)
          else if (number_of_matmuls==one) then
              call sparsemm(nseq, SHS_seq, vectors(1,1,3), vectors(1,1,1), &
                   norb, norbp, ivectorindex, nout, onedimindices)
          end if
        
        
          call vcopy(norb*norbp, vectors(1,1,1), 1, vectors(1,1,2), 1)
        
          !initialize fermi
          call to_zero(norbp*norb, fermi(1,1))
          call to_zero(2*norb*norbp, penalty_ev(1,1,1))
          call compress_polynomial_vector(iproc, nsize_polynomial, orbs, kernel, vectors(1,1,4), chebyshev_polynomials(1,1))
          call axpy_kernel_vectors(norbp, norb, nout, onedimindices, 0.5d0*cc(1,1), vectors(1,1,4), fermi(:,1))
          call axpy_kernel_vectors(norbp, norb, nout, onedimindices, 0.5d0*cc(1,3), vectors(1,1,4), penalty_ev(:,1,1))
          call axpy_kernel_vectors(norbp, norb, nout, onedimindices, 0.5d0*cc(1,3), vectors(1,1,4), penalty_ev(:,1,2))
          call compress_polynomial_vector(iproc, nsize_polynomial, orbs, kernel, vectors(1,1,2), chebyshev_polynomials(1,2))
          call axpy_kernel_vectors(norbp, norb, nout, onedimindices, cc(2,1), vectors(1,1,2), fermi(:,1))
          call axpy_kernel_vectors(norbp, norb, nout, onedimindices, cc(2,3), vectors(1,1,2), penalty_ev(:,1,1))
          call axpy_kernel_vectors(norbp, norb, nout, onedimindices, -cc(2,3), vectors(1,1,2), penalty_ev(:,1,2))
        
        
          do ipl=3,npl
              ! apply (3/2 - 1/2 S) H (3/2 - 1/2 S)
              if (number_of_matmuls==three) then
                  call sparsemm(nseq, ovrlp_compr_seq, vectors(1,1,1), vectors(1,1,2), &
                       norb, norbp, ivectorindex, nout, onedimindices)
                  call sparsemm(nseq, ham_compr_seq, vectors(1,1,2), vectors(1,1,3), &
                       norb, norbp, ivectorindex, nout, onedimindices)
                  call sparsemm(nseq, ovrlp_compr_seq, vectors(1,1,3), vectors(1,1,2), &
                       norb, norbp, ivectorindex, nout, onedimindices)
              else if (number_of_matmuls==one) then
                  call sparsemm(nseq, SHS_seq, vectors(1,1,1), vectors(1,1,2), &
                       norb, norbp, ivectorindex, nout, onedimindices)
              end if
              call axbyz_kernel_vectors(norbp, norb, nout, onedimindices, 2.d0, vectors(1,1,2), &
                   -1.d0, vectors(1,1,4), vectors(1,1,3))
              call compress_polynomial_vector(iproc, nsize_polynomial, orbs, kernel, vectors(1,1,3), chebyshev_polynomials(1,ipl))
              call axpy_kernel_vectors(norbp, norb, nout, onedimindices, cc(ipl,1), vectors(1,1,3), fermi(:,1))
              call axpy_kernel_vectors(norbp, norb, nout, onedimindices, cc(ipl,3), vectors(1,1,3), penalty_ev(:,1,1))
         
              if (mod(ipl,2)==1) then
                  tt=cc(ipl,3)
              else
                  tt=-cc(ipl,3)
              end if
              call axpy_kernel_vectors(norbp, norb, nout, onedimindices, tt, vectors(1,1,3), penalty_ev(:,1,2))
         
              call copy_kernel_vectors(norbp, norb, nout, onedimindices, vectors(1,1,1), vectors(1,1,4))
              call copy_kernel_vectors(norbp, norb, nout, onedimindices, vectors(1,1,3), vectors(1,1,1))
          end do
    
      end if

 
    
      call f_free(vectors)
      call f_free(ham_compr_seq)
      call f_free(ovrlp_compr_seq)
      call f_free(istindexarr)
      call f_free(ivectorindex)
      call f_free_ptr(onedimindices)
    
      if (number_of_matmuls==one) then
          call f_free(matrix)
          call f_free(SHS_seq)
      end if

  end if

  call timing(iproc, 'chebyshev_comp', 'OF')

end subroutine chebyshev_clean



! Performs z = a*x + b*y
subroutine axbyz_kernel_vectors(norbp, norb, nout, onedimindices, a, x, b, y, z)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: norbp, norb, nout
  integer,dimension(4,nout),intent(in) :: onedimindices
  real(8),intent(in) :: a, b
  real(kind=8),dimension(norb,norbp),intent(in) :: x, y
  real(kind=8),dimension(norb,norbp),intent(out) :: z

  ! Local variables
  integer :: i, jorb, iorb

  !$omp parallel default(private) shared(nout, onedimindices,a, b, x, y, z)
  !$omp do
  do i=1,nout
      iorb=onedimindices(1,i)
      jorb=onedimindices(2,i)
      z(jorb,iorb)=a*x(jorb,iorb)+b*y(jorb,iorb)
  end do
  !$omp end do
  !$omp end parallel

end subroutine axbyz_kernel_vectors



subroutine sparsemm(nseq, a_seq, b, c, norb, norbp, ivectorindex, nout, onedimindices)
  use module_base
  use module_types

  implicit none

  !Calling Arguments
  integer, intent(in) :: norb,norbp,nseq
  real(kind=8), dimension(norb,norbp),intent(in) :: b
  real(kind=8), dimension(nseq),intent(in) :: a_seq
  real(kind=8), dimension(norb,norbp), intent(out) :: c
  integer,dimension(nseq),intent(in) :: ivectorindex
  integer,intent(in) :: nout
  integer,dimension(4,nout) :: onedimindices

  !Local variables
  !character(len=*), parameter :: subname='sparsemm'
  integer :: i,jorb,jjorb,m,mp1
  integer :: iorb, ii0, ii2, ilen, jjorb0, jjorb1, jjorb2, jjorb3, jjorb4, jjorb5, jjorb6, iout
  real(kind=8) :: tt


  !$omp parallel default(private) shared(ivectorindex, a_seq, b, c, onedimindices, nout)
  !$omp do
  do iout=1,nout
      i=onedimindices(1,iout)
      iorb=onedimindices(2,iout)
      ilen=onedimindices(3,iout)
      ii0=onedimindices(4,iout)
      ii2=0
      tt=0.d0

      m=mod(ilen,7)
      if (m/=0) then
          do jorb=1,m
             jjorb=ivectorindex(ii0+ii2)
             tt = tt + b(jjorb,i)*a_seq(ii0+ii2)
             ii2=ii2+1
          end do
      end if
      mp1=m+1
      do jorb=mp1,ilen,7

         jjorb0=ivectorindex(ii0+ii2+0)
         tt = tt + b(jjorb0,i)*a_seq(ii0+ii2+0)

         jjorb1=ivectorindex(ii0+ii2+1)
         tt = tt + b(jjorb1,i)*a_seq(ii0+ii2+1)

         jjorb2=ivectorindex(ii0+ii2+2)
         tt = tt + b(jjorb2,i)*a_seq(ii0+ii2+2)

         jjorb3=ivectorindex(ii0+ii2+3)
         tt = tt + b(jjorb3,i)*a_seq(ii0+ii2+3)

         jjorb4=ivectorindex(ii0+ii2+4)
         tt = tt + b(jjorb4,i)*a_seq(ii0+ii2+4)

         jjorb5=ivectorindex(ii0+ii2+5)
         tt = tt + b(jjorb5,i)*a_seq(ii0+ii2+5)

         jjorb6=ivectorindex(ii0+ii2+6)
         tt = tt + b(jjorb6,i)*a_seq(ii0+ii2+6)

         ii2=ii2+7
      end do
      c(iorb,i)=tt
  end do 
  !$omp end do
  !$omp end parallel

    
end subroutine sparsemm



subroutine copy_kernel_vectors(norbp, norb, nout, onedimindices, a, b)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: norbp, norb, nout
  integer,dimension(4,nout),intent(in) :: onedimindices
  real(kind=8),dimension(norb,norbp),intent(in) :: a
  real(kind=8),dimension(norb,norbp),intent(out) :: b

  ! Local variables
  integer :: i, jorb, iorb


  !$omp parallel default(private) shared(nout, onedimindices,a, b)
  !$omp do
  do i=1,nout
      iorb=onedimindices(1,i)
      jorb=onedimindices(2,i)
      b(jorb,iorb)=a(jorb,iorb)
  end do
  !$omp end do
  !$omp end parallel


end subroutine copy_kernel_vectors




subroutine axpy_kernel_vectors(norbp, norb, nout, onedimindices, a, x, y)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: norbp, norb, nout
  integer,dimension(4,nout),intent(in) :: onedimindices
  real(kind=8),intent(in) :: a
  real(kind=8),dimension(norb,norbp),intent(in) :: x
  real(kind=8),dimension(norb,norbp),intent(inout) :: y

  ! Local variables
  integer :: i, jorb, iorb

  !$omp parallel default(private) shared(nout, onedimindices, y, x, a)
  !$omp do
  do i=1,nout
      iorb=onedimindices(1,i)
      jorb=onedimindices(2,i)
      y(jorb,iorb)=y(jorb,iorb)+a*x(jorb,iorb)
  end do
  !$omp end do
  !$omp end parallel


end subroutine axpy_kernel_vectors




subroutine determine_sequential_length(norbp, isorb, norb, foe_obj, sparsemat, nseq, nmaxsegk, nmaxvalk)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: norbp, isorb, norb
  type(foe_data),intent(in) :: foe_obj
  type(sparseMatrix),intent(in) :: sparsemat
  integer,intent(out) :: nseq, nmaxsegk, nmaxvalk

  ! Local variables
  integer :: i,iseg,jorb,iorb,jseg,ii

  nseq=0
  nmaxsegk=0
  nmaxvalk=0
  do i = 1,norbp
     ii=isorb+i
     nmaxsegk=max(nmaxsegk,foe_obj%kernel_nseg(ii))
     do iseg=1,foe_obj%kernel_nseg(ii)
          nmaxvalk=max(nmaxvalk,foe_obj%kernel_segkeyg(2,iseg,ii)-foe_obj%kernel_segkeyg(1,iseg,ii)+1)
          do iorb=foe_obj%kernel_segkeyg(1,iseg,ii),foe_obj%kernel_segkeyg(2,iseg,ii)
              do jseg=sparsemat%istsegline(iorb),sparsemat%istsegline(iorb)+sparsemat%nsegline(iorb)-1
                  do jorb = sparsemat%keyg(1,jseg),sparsemat%keyg(2,jseg)
                      nseq=nseq+1
                  end do
              end do
          end do
     end do
  end do 

end subroutine determine_sequential_length




subroutine init_onedimindices(norbp, isorb, foe_obj, sparsemat, nout, onedimindices)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: norbp, isorb
  type(foe_data),intent(in) :: foe_obj
  type(sparseMatrix),intent(in) :: sparsemat
  integer,intent(out) :: nout
  integer,dimension(:,:),pointer :: onedimindices

  ! Local variables
  integer :: i, iii, iseg, iorb, istat, ii, jseg, ilen, itot
  character(len=*),parameter :: subname='init_onedimindices'


  nout=0
  do i = 1,norbp
     iii=isorb+i
     do iseg=1,foe_obj%kernel_nseg(iii)
          do iorb=foe_obj%kernel_segkeyg(1,iseg,iii),foe_obj%kernel_segkeyg(2,iseg,iii)
              nout=nout+1
          end do
      end do
  end do

! allocate(onedimindices(4,nout), stat=istat)
  onedimindices = f_malloc_ptr((/ 4, nout /),id='onedimindices')

  ii=0
  itot=1
  do i = 1,norbp
     iii=isorb+i
     do iseg=1,foe_obj%kernel_nseg(iii)
          do iorb=foe_obj%kernel_segkeyg(1,iseg,iii),foe_obj%kernel_segkeyg(2,iseg,iii)
              ii=ii+1
              onedimindices(1,ii)=i
              onedimindices(2,ii)=iorb
              ilen=0
              do jseg=sparsemat%istsegline(iorb),sparsemat%istsegline(iorb)+sparsemat%nsegline(iorb)-1
                  ilen=ilen+sparsemat%keyg(2,jseg)-sparsemat%keyg(1,jseg)+1
              end do
              onedimindices(3,ii)=ilen
              onedimindices(4,ii)=itot
              itot=itot+ilen
          end do
      end do
  end do

end subroutine init_onedimindices



subroutine get_arrays_for_sequential_acces(norbp, isorb, norb, foe_obj, sparsemat, nseq, nmaxsegk, nmaxvalk, &
           istindexarr, ivectorindex)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: norbp, isorb, norb, nseq, nmaxsegk, nmaxvalk
  type(foe_data),intent(in) :: foe_obj
  type(sparseMatrix),intent(in) :: sparsemat
  integer,dimension(nmaxvalk,nmaxsegk,norbp),intent(out) :: istindexarr
  integer,dimension(nseq),intent(out) :: ivectorindex

  ! Local variables
  integer :: i,iseg,jorb,jjorb,iorb,jseg,ii,iii


  ii=1
  do i = 1,norbp
     iii=isorb+i
     do iseg=1,foe_obj%kernel_nseg(iii)
          do iorb=foe_obj%kernel_segkeyg(1,iseg,iii),foe_obj%kernel_segkeyg(2,iseg,iii)
              istindexarr(iorb-foe_obj%kernel_segkeyg(1,iseg,iii)+1,iseg,i)=ii
              do jseg=sparsemat%istsegline(iorb),sparsemat%istsegline(iorb)+sparsemat%nsegline(iorb)-1
                  do jorb = sparsemat%keyg(1,jseg),sparsemat%keyg(2,jseg)
                      jjorb = jorb - (iorb-1)*norb
                      ivectorindex(ii)=jjorb
                      ii = ii+1
                  end do
              end do
          end do
     end do
  end do 

end subroutine get_arrays_for_sequential_acces




subroutine sequential_acces_matrix(norbp, isorb, norb, foe_obj, sparsemat, a, nseq, nmaxsegk, nmaxvalk, a_seq)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: norbp, isorb, norb, nseq, nmaxsegk, nmaxvalk
  type(foe_data),intent(in) :: foe_obj
  type(sparseMatrix),intent(in) :: sparsemat
  real(kind=8),dimension(sparsemat%nvctr),intent(in) :: a
  real(kind=8),dimension(nseq),intent(out) :: a_seq

  ! Local variables
  integer :: i,iseg,jorb,jj,iorb,jseg,ii,iii


  ii=1
  do i = 1,norbp
     iii=isorb+i
     do iseg=1,foe_obj%kernel_nseg(iii)
          do iorb=foe_obj%kernel_segkeyg(1,iseg,iii),foe_obj%kernel_segkeyg(2,iseg,iii)
              do jseg=sparsemat%istsegline(iorb),sparsemat%istsegline(iorb)+sparsemat%nsegline(iorb)-1
                  jj=1
                  do jorb = sparsemat%keyg(1,jseg),sparsemat%keyg(2,jseg)
                      a_seq(ii)=a(sparsemat%keyv(jseg)+jj-1)
                      jj = jj+1
                      ii = ii+1
                  end do
              end do
          end do
     end do
  end do 

end subroutine sequential_acces_matrix



subroutine chebyshev_fast(iproc, nsize_polynomial, npl, orbs, fermi, chebyshev_polynomials, cc, kernelp)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nsize_polynomial, npl
  type(orbitals_data),intent(in) :: orbs
  type(sparseMatrix),intent(in) :: fermi
  real(kind=8),dimension(nsize_polynomial,npl),intent(in) :: chebyshev_polynomials
  real(kind=8),dimension(npl),intent(in) :: cc
  real(kind=8),dimension(orbs%norb,orbs%norbp),intent(out) :: kernelp

  ! Local variables
  integer :: ipl, istat, iall
  real(kind=8),dimension(:),allocatable :: kernel_compressed
  character(len=*),parameter :: subname='chebyshev_fast'


  if (nsize_polynomial>0) then
!     allocate(kernel_compressed(nsize_polynomial),stat=istat)
      kernel_compressed = f_malloc(nsize_polynomial,id='kernel_compressed')

      call to_zero(nsize_polynomial,kernel_compressed(1))
      !write(*,*) 'ipl, first element', 1, chebyshev_polynomials(1,1)
      call daxpy(nsize_polynomial, 0.5d0*cc(1), chebyshev_polynomials(1,1), 1, kernel_compressed(1), 1)
      do ipl=2,npl
      !write(*,*) 'ipl, first element', ipl, chebyshev_polynomials(1,ipl)
          call daxpy(nsize_polynomial, cc(ipl), chebyshev_polynomials(1,ipl), 1, kernel_compressed(1), 1)
      end do

      call uncompress_polynomial_vector(iproc, nsize_polynomial, orbs, fermi, kernel_compressed, kernelp)

      iall=-product(shape(kernel_compressed))*kind(kernel_compressed)
!     deallocate(kernel_compressed, stat=istat)
      call f_free(kernel_compressed)
  end if

end subroutine chebyshev_fast
