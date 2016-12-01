!> @file
!!   File containing the routines to calculate the coefficients (eigenvectors) of the Hamiltonian matrix
!! @author
!!   Copyright (C) 2016 CheSS developers
!!
!!   This file is part of CheSS.
!!   
!!   CheSS is free software: you can redistribute it and/or modify
!!   it under the terms of the GNU Lesser General Public License as published by
!!   the Free Software Foundation, either version 3 of the License, or
!!   (at your option) any later version.
!!   
!!   CheSS is distributed in the hope that it will be useful,
!!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!   GNU Lesser General Public License for more details.
!!   
!!   You should have received a copy of the GNU Lesser General Public License
!!   along with CheSS.  If not, see <http://www.gnu.org/licenses/>.


!> Module coefficients
module coeffs
  use sparsematrix_base
  !use dynamic_memory
  use wrapper_mpi
  !use yaml_output
  use wrapper_linalg
  use futile
  implicit none

  private

  !> Public routines
  public :: get_coeffs_diagonalization
  public :: calculate_kernel_and_energy
  ! SM: Should try to make the following routine private
  public :: calculate_density_kernel 

  contains

    subroutine get_coeffs_diagonalization(iproc, nproc, comm, nfvctr, norbu, norbd, norb, blocksize_pdsyev, &
               smats, smatm, ovrlp, ham, coeff, eval_all, eval_occup, info_coeff)
      use futile
      use wrapper_mpi
      use wrapper_linalg, only: vcopy
      use sparsematrix, only: uncompress_matrix2, diagonalizehamiltonian2
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm, nfvctr, norbu, norbd, norb, blocksize_pdsyev
      type(sparse_matrix),intent(in) :: smats, smatm
      type(matrices),intent(in) :: ovrlp, ham
      real(mp),dimension(nfvctr,norb),intent(out) :: coeff
      real(mp),dimension(nfvctr*smats%nspin),intent(out) :: eval_all
      real(mp),dimension(norb),intent(out) :: eval_occup
      integer,intent(out) :: info_coeff
    
      ! Local variables
      integer :: ispin, iorb
      real(mp) :: maxdiff
      real(mp),dimension(:),allocatable :: eval
      real(mp),dimension(:,:,:),allocatable :: ovrlp_full, ham_full
    
      ovrlp_full = f_malloc((/nfvctr,nfvctr,smats%nspin/),id='ovrlp_full')
      ham_full = f_malloc((/nfvctr,nfvctr,smatm%nspin/),id='ham_full')
    
      call uncompress_matrix2(iproc, nproc, comm, smats, ovrlp%matrix_compr, ovrlp_full)
      call uncompress_matrix2(iproc, nproc, comm, smatm, ham%matrix_compr, ham_full)
    
      ! Keep the Hamiltonian and the overlap since they will be overwritten by the diagonalization.
      eval = f_malloc(nfvctr,id='eval')
    
      do ispin=1,smats%nspin
          if (iproc==0) call yaml_map('method','diagonalization')
          call diagonalizeHamiltonian2(iproc, nproc, comm, &
               blocksize_pdsyev, nfvctr, &
               ham_full(1,1,ispin), ovrlp_full(1,1,ispin), eval)
    
          ! Broadcast the results (eigenvectors and eigenvalues) from task 0 to
          ! all other tasks (in this way avoiding that different MPI tasks have different values)
          if (nproc>1) then
              if (iproc==0) call yaml_mapping_open('Cross-check among MPI tasks')
              call mpibcast(ham_full(:,:,ispin), comm=comm, maxdiff=maxdiff)
              if (iproc==0) call yaml_map('max diff of eigenvectors',maxdiff,fmt='(es8.2)')
              call mpibcast(eval, comm=comm, maxdiff=maxdiff)
              if (iproc==0) call yaml_map('max diff of eigenvalues',maxdiff,fmt='(es8.2)')
              if (iproc==0) call yaml_mapping_close()
          end if
    
          ! copy all the eigenvalues
          call vcopy(nfvctr, eval(1), 1, eval_all((ispin-1)*nfvctr+1), 1)
          ! copy the eigenvalues of the occupied states
          if (ispin==1) then
              call vcopy(norbu, eval(1), 1, eval_occup(1), 1)
          else
              call vcopy(norbd, eval(1), 1, eval_occup(norbu+1), 1)
          end if
    
          ! Make sure that the eigenvectors have the same sign on all MPI tasks.
          ! To do so, ensure that the first entry is always positive.
          do iorb=1,nfvctr
              if (ham_full(1,iorb,1)<0.d0) then
                  call dscal(nfvctr, -1.d0, ham_full(1,iorb,1), 1)
              end if
          end do
    
          ! Copy the diagonalized matrix to the coeff array.
          if (smats%nspin/=1) then
              ! Only copy the occupied states
              if (ispin==1) then
                  call vcopy(norbu*nfvctr, ham_full(1,1,ispin), 1, coeff(1,1), 1)
              else if (ispin==2) then
                  call vcopy(norbd*nfvctr, ham_full(1,1,ispin), 1, coeff(1,norbu+1), 1)
              end if
          else
              ! Copy all states
              call vcopy(norb*nfvctr, ham_full(1,1,1), 1, coeff(1,1), 1)
          end if
          info_coeff = 0
      end do
    
      call f_free(eval)
      call f_free(ovrlp_full)
      call f_free(ham_full)
    
    end subroutine get_coeffs_diagonalization



    subroutine calculate_kernel_and_energy(iproc,nproc,comm,denskern,ham,denskern_mat,ham_mat, &
               energy,coeff,norbp,isorb,norbu,norb,occup,calculate_kernel)
      use sparsematrix_highlevel, only: trace_AB
      implicit none
      integer, intent(in) :: iproc, nproc, comm, norbp, isorb, norbu, norb
      real(kind=8),dimension(norb),intent(in) :: occup
      type(sparse_matrix), intent(in) :: ham
      type(sparse_matrix), intent(in) :: denskern
      type(matrices),intent(in) :: ham_mat
      type(matrices),intent(out) :: denskern_mat
      logical, intent(in) :: calculate_kernel
      real(kind=mp), intent(out) :: energy
      real(kind=mp), dimension(denskern%nfvctr,norb), intent(in) :: coeff
    
      integer :: iorb, jorb, ind_ham, ind_denskern, ierr, iorbp, is, ie, ispin
    
      if (calculate_kernel) then 
         !!call extract_taskgroup_inplace(denskern, denskern_mat)
         call calculate_density_kernel(iproc, nproc, comm, .true., norbp, isorb, norbu, norb, occup, &
              coeff, denskern, denskern_mat)
         !call gather_matrix_from_taskgroups_inplace(iproc, nproc, denskern, denskern_mat)
         !denskern%matrix_compr = denskern_mat%matrix_compr
      end if
    
      call timing(iproc,'calc_energy','ON')
      energy=0.0_mp
      do ispin=1,denskern%nspin
          energy = energy + trace_AB(iproc, nproc, comm, ham, denskern, ham_mat, denskern_mat, ispin)
      end do
    
      !!do iorbp=1,tmb_orbs%norbp
      !!   iorb=iorbp+tmb_orbs%isorb
      !!   if (tmb_orbs%spinsgn(iorb)>0.d0) then
      !!       ! spin up support function or non-polarized case
      !!       is=1
      !!       ie=tmb_orbs%norbu
      !!       !ispin=1
      !!   else
      !!       ! spin down support function
      !!       is=tmb_orbs%norbu+1
      !!       ie=tmb_orbs%norb
      !!       !ispin=2
      !!   end if
      !!   !$omp parallel default(private) shared(is,ie,iorb,denskern,ham,denskern_mat,ham_mat,tmb_orbs,energy,ispin)
      !!   !$omp do reduction(+:energy)
      !!   !do jorb=1,tmb_orbs%norb
      !!   do jorb=is,ie
      !!      ind_ham = matrixindex_in_compressed(ham,iorb,jorb)
      !!      ind_denskern = matrixindex_in_compressed(denskern,jorb,iorb)
      !!      if (ind_ham==0.or.ind_denskern==0) cycle
      !!      energy = energy + &
      !!          denskern_mat%matrix_compr(ind_denskern-denskern%isvctrp_tg)*ham_mat%matrix_compr(ind_ham-ham%isvctrp_tg)
      !!          !!write(*,'(a,5i8,2es16.7)') 'iorb, jorb, ispin, ind_denskern, ind_ham, val_denskern, val_ham', &
      !!          !!    iorb, jorb, ispin, mod(ind_denskern-denskern%isvctrp_tg-1,denskern%nvctr)+1, mod(ind_ham-ham%isvctrp_tg-1,ham%nvctr)+1, denskern_mat%matrix_compr(ind_denskern-denskern%isvctrp_tg), ham_mat%matrix_compr(ind_ham-ham%isvctrp_tg)
      !!   end do
      !!   !$omp end do
      !!   !$omp end parallel
      !!end do
      !!if (nproc>1) then
      !!   call mpiallred(energy, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
      !!end if
      call timing(iproc,'calc_energy','OF')
    
    end subroutine calculate_kernel_and_energy


    subroutine calculate_density_kernel(iproc, nproc, comm, isKernel, norbp, isorb, norbu, norb, occup, &
               coeff, denskern, denskern_, keep_uncompressed_)
      use sparsematrix, only: compress_matrix, extract_taskgroup
      implicit none
    
      ! Calling arguments
      integer,intent(in):: iproc, nproc, comm, norbp, isorb, norbu, norb
      real(kind=8),dimension(norb),intent(in) :: occup
      logical, intent(in) :: isKernel
      type(sparse_matrix), intent(in) :: denskern
      real(kind=8),dimension(denskern%nfvctr,norb),intent(in):: coeff   !only use the first (occupied) orbitals
      type(matrices), intent(out) :: denskern_
      logical,intent(in),optional :: keep_uncompressed_ !< keep the uncompressed kernel in denskern_%matrix (requires that this array is already allocated outside of the routine)
    
      ! Local variables
      integer :: ierr, sendcount, jproc, iorb, itmb, iiorb, ispin, jorb
      real(kind=8),dimension(:,:),allocatable :: density_kernel_partial, fcoeff
      real(kind=8),dimension(:),allocatable :: tmparr
    ! real(kind=8), dimension(:,:,), allocatable :: ks,ksk,ksksk
      character(len=*),parameter :: subname='calculate_density_kernel'
      integer,dimension(:),allocatable :: recvcounts, dspls
      integer,parameter :: ALLGATHERV=1, ALLREDUCE=2
      integer,parameter :: communication_strategy=ALLREDUCE
      logical :: keep_uncompressed
    
      call f_routine(id='calculate_density_kernel')
    
      if (present(keep_uncompressed_)) then
          keep_uncompressed = keep_uncompressed_
      else
          keep_uncompressed = .false.
      end if
    
      if (keep_uncompressed) then
          if (.not.associated(denskern_%matrix)) stop 'ERROR: denskern_%matrix must be associated if keep_uncompressed is true'
      end if
    
      !!write(*,*) 'iproc, orbs_tmb%norbp, orbs_tmb%isorb, orbs_tmb%norb', &
      !!            iproc, orbs_tmb%norbp, orbs_tmb%isorb, orbs_tmb%norb
      !!write(*,*) 'iproc, denskern%nfvctrp, denskern%isfvctr, denskern%nfvctr', &
      !!            iproc, denskern%nfvctrp, denskern%isfvctr, denskern%nfvctr
    
      if (communication_strategy==ALLGATHERV) then
          if (iproc==0) call yaml_map('communication strategy kernel','ALLGATHERV')
          stop 'calculate_density_kernel: ALLGATHERV option needs reworking due to the spin'
          call timing(iproc,'calc_kernel','ON')
          !if(iproc==0) write(*,'(1x,a)',advance='no') 'calculate density kernel... '
          density_kernel_partial=f_malloc((/denskern%nfvctr,max(denskern%nfvctrp,1)/), id='density_kernel_partial')
          fcoeff=f_malloc0((/denskern%nfvctrp,norb/), id='fcoeff')
          if(denskern%nfvctrp>0) then
              !decide whether we calculate the density kernel or just transformation matrix
              if(isKernel) then
                 do iorb=1,norb
                    !call daxpy(denskern%nfvctrp,orbs%occup(iorb),coeff(1+denskern%isfvctr,iorb),1,fcoeff(1+denskern%isfvctr,iorb),1)
                    do itmb=1,denskern%nfvctrp
                         fcoeff(itmb,iorb) = occup(iorb)*coeff(denskern%isfvctr+itmb,iorb)
                    end do
                 end do
              else
                 do iorb=1,norb
                    do itmb=1,denskern%nfvctrp
                         fcoeff(itmb,iorb) = coeff(denskern%isfvctr+itmb,iorb)
                    end do
                 end do
              end if
    
              call dgemm('n', 't', denskern%nfvctr, denskern%nfvctrp, norb, 1.d0, coeff(1,1), denskern%nfvctr, &
                   fcoeff(1,1), denskern%nfvctrp, 0.d0, density_kernel_partial(1,1), denskern%nfvctr)
          end if
          call f_free(fcoeff)
          call timing(iproc,'calc_kernel','OF')
    
          call timing(iproc,'waitAllgatKern','ON')
          call mpibarrier(comm)
          call timing(iproc,'waitAllgatKern','OF')
    
          !denskern_%matrix=f_malloc_ptr((/denskern%nfvctr,denskern%nfvctr/), id='denskern_%matrix')
    
          if (.not.keep_uncompressed) then
              denskern_%matrix=sparsematrix_malloc_ptr(denskern,iaction=DENSE_FULL,id='denskern_%matrix')
          end if
    
          if (nproc > 1) then
             call timing(iproc,'commun_kernel','ON')
             recvcounts=f_malloc((/0.to.nproc-1/),id='recvcounts')
             dspls=f_malloc((/0.to.nproc-1/),id='dspls')
             do jproc=0,nproc-1
                 recvcounts(jproc)=denskern%nfvctr*denskern%nfvctr_par(jproc)
                 dspls(jproc)=denskern%nfvctr*denskern%isfvctr_par(jproc)
             end do
             sendcount=denskern%nfvctr*denskern%nfvctrp
             call mpi_allgatherv(density_kernel_partial(1,1), sendcount, mpi_double_precision, &
                  denskern_%matrix(1,1,1), recvcounts, dspls, mpi_double_precision, &
                  comm, ierr)
             call f_free(recvcounts)
             call f_free(dspls)
             call timing(iproc,'commun_kernel','OF')
          else
             call vcopy(denskern%nfvctr*denskern%nfvctrp,density_kernel_partial(1,1),1,denskern_%matrix(1,1,1),1)
          end if
    
          call f_free(density_kernel_partial)
    
          call compress_matrix(iproc,nproc,denskern,inmat=denskern_%matrix,outmat=denskern_%matrix_compr)
          if (.not.keep_uncompressed) then
              call f_free_ptr(denskern_%matrix)
          end if
      else if (communication_strategy==ALLREDUCE) then
          if (iproc==0) call yaml_map('communication strategy kernel','ALLREDUCE')
          call timing(iproc,'calc_kernel','ON')
          !!if(iproc==0) write(*,'(1x,a)',advance='no') 'calculate density kernel... '
          !denskern_%matrix=f_malloc_ptr((/denskern%nfvctr,denskern%nfvctr/), id='denskern_%matrix_compr')
          if (.not.keep_uncompressed) then
              denskern_%matrix=sparsematrix_malloc_ptr(denskern,iaction=DENSE_FULL,id='denskern_%matrix')
          end if
          if(norbp>0) then
              fcoeff=f_malloc((/denskern%nfvctr,norbp/), id='fcoeff')
              !decide wether we calculate the density kernel or just transformation matrix
              if(isKernel)then
                 do iorb=1,norbp
                    !call f_zero(denskern%nfvctr,f_coeff(1,iorb))
                    !call daxpy(denskern%nfvctr,orbs%occup(orbs%isorb+iorb),coeff(1,orbs%isorb+iorb),1,fcoeff(1,iorb),1)
                    !write(*,*) 'iorb, occup', iorb, orbs%occup(orbs%isorb+iorb)
                    do itmb=1,denskern%nfvctr
                        fcoeff(itmb,iorb) = occup(isorb+iorb)*coeff(itmb,isorb+iorb)
                    end do
                 end do
              else
                 do iorb=1,norbp
                    call vcopy(denskern%nfvctr,coeff(1,isorb+iorb),1,fcoeff(1,iorb),1)
                 end do
              end if
          !!if (iproc==0) then
          !!    do iorb=1,orbs%norbp
          !!        do jorb=1,denskern%nfvctr
          !!            write(970,'(a,2i9,f14.7)') 'iorb, jorb, fcoeff(jorb,iorb)', iorb, jorb, fcoeff(jorb,iorb)
          !!        end do
          !!    end do
          !!end if
              !call dgemm('n', 't', denskern%nfvctr, denskern%nfvctr, orbs%norbp, 1.d0, coeff(1,orbs%isorb+1), denskern%nfvctr, &
              !     fcoeff(1,1), denskern%nfvctr, 0.d0, denskern_%matrix(1,1,1), denskern%nfvctr)
              call f_zero(denskern%nspin*denskern%nfvctr**2, denskern_%matrix(1,1,1))
              !!write(*,*) 'iproc, orbs%spinsgn',iproc, orbs%spinsgn
              !!write(*,*) 'iproc, orbs%norbu, orbs%norbd', iproc, orbs%norbu, orbs%norbd 
              do iorb=1,norbp
                  iiorb=isorb+iorb
                  !if (orbs%spinsgn(iiorb)>0.d0) then
                  if (iiorb<=norbu) then
                      ispin=1
                  else
                      ispin=2
                  end if
                  call dgemm('n', 't', denskern%nfvctr, denskern%nfvctr, 1, 1.d0, coeff(1,isorb+iorb), denskern%nfvctr, &
                       fcoeff(1,iorb), denskern%nfvctr, 1.d0, denskern_%matrix(1,1,ispin), denskern%nfvctr)
              end do
              call f_free(fcoeff)
          else
              call f_zero(denskern%nspin*denskern%nfvctr**2, denskern_%matrix(1,1,1))
          end if
          call timing(iproc,'calc_kernel','OF')
    
          !!if (iproc==0) then
          !!    do ispin=1,denskern%nspin
          !!        do iorb=1,denskern%nfvctr
          !!            do jorb=1,denskern%nfvctr
          !!                write(940+ispin,'(a,3i9,f14.7)') 'ispin, iorb, jorb, denskern_%matrix(jorb,iorb,ispin)', ispin, iorb, jorb, denskern_%matrix(jorb,iorb,ispin)
          !!            end do
          !!        end do
          !!    end do
          !!end if
    
          call timing(iproc,'waitAllgatKern','ON')
          call mpi_barrier(comm,ierr)
          call timing(iproc,'waitAllgatKern','OF')
    
          tmparr = sparsematrix_malloc(denskern,iaction=SPARSE_FULL,id='tmparr')
          call compress_matrix(iproc,nproc,denskern,inmat=denskern_%matrix,outmat=tmparr)
          if (keep_uncompressed) then
              if (nproc > 1) then
                  call timing(iproc,'commun_kernel','ON')
                  call mpiallred(denskern_%matrix(1,1,1), denskern%nspin*denskern%nfvctr**2, &
                       mpi_sum, comm=comm)
                  call timing(iproc,'commun_kernel','OF')
              end if
          end if
          if (.not.keep_uncompressed) then
              call f_free_ptr(denskern_%matrix)
          end if
          if (nproc > 1) then
              call timing(iproc,'commun_kernel','ON')
              call mpiallred(tmparr(1), denskern%nspin*denskern%nvctr, mpi_sum, comm=comm)
              call timing(iproc,'commun_kernel','OF')
          end if
          call extract_taskgroup(denskern, tmparr, denskern_%matrix_compr)
          call f_free(tmparr)
    
          !call compress_matrix(iproc,denskern)
          !call f_free_ptr(denskern%matrix)
      end if
      call f_release_routine()
    
     ! Purify Kernel
     !call dgemm('n', 't', orbs_tmb%norb, orbs_tmb%norb, orbs_tmb%norbp, 1.d0, kernel(1,orbs_tmb%isorb+1), orbs_tmb%norb, &
     !           overlap(1,orbs_tmb%isorb+1), orbs_tmb%norb, 0.d0, ks(1,1), orbs_tmb%norb) 
     !call dgemm('n', 't', orbs_tmb%norb, orbs_tmb%norb, orbs_tmb%norbp, 1.d0, ks(1,orbs_tmb%isorb+1), orbs_tmb%norb, &
     !           kernel(1,orbs_tmb%isorb+1), orbs_tmb%norb, 0.d0, ksk(1,1), orbs_tmb%norb)
     !call dgemm('n', 't', orbs_tmb%norb, orbs_tmb%norb, orbs_tmb%norbp, 1.d0, ks(1,orbs_tmb%isorb+1), orbs_tmb%norb, &
     !           ksk(1,orbs_tmb%isorb+1), orbs_tmb%norb, 0.d0, ksksk(1,1), orbs_tmb%norb)
    
     !!if(present(overlap)) then
       !!allocate(ks(orbs_tmb%norb,orbs_tmb%norb),stat=istat)
       !!call memocc(istat, ks, 'ks', subname) 
       !!allocate(ksk(orbs_tmb%norb,orbs_tmb%norb),stat=istat)
       !!call memocc(istat, ksk, 'ksk', subname) 
       !!allocate(ksksk(orbs_tmb%norb,orbs_tmb%norb),stat=istat)
       !!call memocc(istat, ksksk, 'ksksk', subname) 
    
       !!call dgemm('n', 't', orbs_tmb%norb, orbs_tmb%norb, orbs_tmb%norb, 1.d0, kernel(1,1), orbs_tmb%norb, &
       !!           overlap(1,1), orbs_tmb%norb, 0.d0, ks(1,1), orbs_tmb%norb) 
       !!call dgemm('n', 't', orbs_tmb%norb, orbs_tmb%norb, orbs_tmb%norb, 1.d0, ks(1,1), orbs_tmb%norb, &
       !!           kernel(1,1), orbs_tmb%norb, 0.d0, ksk(1,1), orbs_tmb%norb)
       !!call dgemm('n', 't', orbs_tmb%norb, orbs_tmb%norb, orbs_tmb%norb, 1.d0, ks(1,1), orbs_tmb%norb, &
       !!           ksk(1,1), orbs_tmb%norb, 0.d0, ksksk(1,1), orbs_tmb%norb)
       !!print *,'PURIFYING THE KERNEL'
       !!kernel = 3*ksk-2*ksksk
       !!
       !!iall = -product(shape(ks))*kind(ks)
       !!deallocate(ks,stat=istat)
       !!call memocc(istat, iall, 'ks', subname)
       !!iall = -product(shape(ksk))*kind(ksk)
       !!deallocate(ksk,stat=istat)
       !!call memocc(istat, iall, 'ksk', subname)
       !!iall = -product(shape(ksksk))*kind(ksksk)
       !!deallocate(ksksk,stat=istat)
       !!call memocc(istat, iall, 'ksksk', subname)
     !!end if
    
    
    end subroutine calculate_density_kernel

end module coeffs
