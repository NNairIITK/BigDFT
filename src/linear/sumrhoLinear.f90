!> @file 
!!   sumrho: linear version
!! @author
!!   Copyright (C) 2013-2014 BigDFT group 
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 


!> Here starts the routine for building partial density inside the localisation region
!! This routine should be treated as a building-block for the linear scaling code
subroutine local_partial_densityLinear(nproc,rsflag,nscatterarr,&
     nrhotot,Lzd,hxh,hyh,hzh,xc,nspin,orbs,mapping,psi,rho)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => local_partial_densityLinear
  use module_xc
  use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
  implicit none
  logical, intent(in) :: rsflag
  integer, intent(in) :: nproc
  integer,intent(in) :: nrhotot
  integer, intent(in) :: nspin
  real(gp), intent(in) :: hxh,hyh,hzh
  type(xc_info), intent(in) :: xc
  type(local_zone_descriptors), intent(in) :: Lzd
  type(orbitals_data),intent(in) :: orbs
  integer,dimension(orbs%norb),intent(in) :: mapping
  integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
  real(dp),dimension(max(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*nrhotot,1),max(nspin,orbs%nspinor)),intent(out) :: rho
  !local variables
  character(len=*), parameter :: subname='local_partial_densityLinear'
  integer :: iorb,i_stat,i_all,ii, ind, indSmall, indLarge
  integer :: oidx,sidx,nspinn,npsir,ncomplex, i1, i2, i3, ilr, ispin
  integer :: nspincomp,ii1,ii2,ii3
  real(gp) :: hfac,spinval
  type(workarr_sumrho) :: w
  real(wp), dimension(:,:), allocatable :: psir
  real(dp), dimension(:),allocatable :: rho_p
  integer, dimension(:,:), allocatable :: Lnscatterarr
 !components of wavefunction in real space which must be considered simultaneously
  !and components of the charge density
  if (orbs%nspinor ==4) then
     npsir=4
     nspinn=4
     ncomplex=0
  else
     npsir=1
     nspinn=nspin
     ncomplex=orbs%nspinor-1
  end if
  nspincomp = 1
  if (nspin > 1) nspincomp = 2

 !allocate and define Lnscatterarr which is just a fake
  Lnscatterarr = f_malloc((/ 0.to.nproc-1, 1.to.4 /),id='Lnscatterarr')
  Lnscatterarr(:,3) = 0
  Lnscatterarr(:,4) = 0

  !initialize the rho array at 10^-20 instead of zero, due to the invcb ABINIT routine
  !otherwise use libXC routine
  call xc_init_rho(xc, max(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*nrhotot,1)*max(nspin,orbs%nspinor),rho,nproc)

  ind=1
  orbitalsLoop: do ii=1,orbs%norbp

     iorb = ii + orbs%isorb
     ilr = orbs%inwhichLocreg(iorb)

     Lnscatterarr(:,1) = Lzd%Llr(ilr)%d%n3i 
     Lnscatterarr(:,2) = Lzd%Llr(ilr)%d%n3i 


     call initialize_work_arrays_sumrho(1,Lzd%Llr(ilr),.true.,w)
     rho_p = f_malloc0(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*nspinn,id='rho_p')
     psir = f_malloc((/ Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i, npsir /),id='psir')
  
     if (Lzd%Llr(ilr)%geocode == 'F') then
        call f_zero(psir)
     end if
 
     !Need to zero rho_p
     !call f_zero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*nspinn, rho_p)

     !print *,'norbp',orbs%norbp,orbs%norb,orbs%nkpts,orbs%kwgts,orbs%iokpt,orbs%occup
     !hfac=orbs%kwgts(orbs%iokpt(ii))*(orbs%occup(iorb)/(hxh*hyh*hzh))
     hfac=orbs%kwgts(orbs%iokpt(ii))*(orbs%occup(mapping(iorb))/(hxh*hyh*hzh))
     spinval=orbs%spinsgn(iorb)

     !this shoudl have aleady been defined
     if (Lzd%Llr(ilr)%hybrid_on) stop 'ERROR, ilr not initialized'
     !Lzd%Llr(ilr)%hybrid_on=.false.

     if (hfac /= 0.d0) then

        !sum for complex function case, npsir=1 in that case
        do oidx=0,ncomplex

           do sidx=1,npsir
              call daub_to_isf(Lzd%Llr(ilr),w,psi(ind),psir(1,sidx))
              ind=ind+Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f
           end do
           

           select case(Lzd%Llr(ilr)%geocode)
           case('F')
              !write(*,*) 'WARNING: MODIFIED CALLING SEQUENCE OF partial_density_free!!!!'
              call partial_density_free((rsflag .and. .not. Lzd%linear),nproc,Lzd%Llr(ilr)%d%n1i,&
                   Lzd%Llr(ilr)%d%n2i,Lzd%Llr(ilr)%d%n3i,npsir,nspinn,Lzd%Llr(ilr)%d%n3i,&!nrhotot,&
                   hfac,Lnscatterarr,spinval,psir,rho_p,Lzd%Llr(ilr)%bounds%ibyyzz_r)
           case('P')

              call partial_density(rsflag,nproc,Lzd%Llr(ilr)%d%n1i,Lzd%Llr(ilr)%d%n2i,Lzd%Llr(ilr)%d%n3i,&
                   npsir,nspinn,Lzd%Llr(ilr)%d%n3i,&!nrhotot,&
                   hfac,nscatterarr,spinval,psir,rho_p)

           case('S')

              call partial_density(rsflag,nproc,Lzd%Llr(ilr)%d%n1i,Lzd%Llr(ilr)%d%n2i,Lzd%Llr(ilr)%d%n3i,&
                   npsir,nspinn,Lzd%Llr(ilr)%d%n3i,&!nrhotot,&
                   hfac,nscatterarr,spinval,psir,rho_p)

           end select

           ! Copy rho_p to the correct place in rho
           indSmall=0
           do ispin=1,nspinn
               do i3=1,Lzd%Llr(ilr)%d%n3i !min(Lzd%Llr(ilr)%d%n3i,nscatterarr(iproc,1)) 
                   ii3 = i3 + Lzd%Llr(ilr)%nsi3 - 1
                   if(ii3 < 0 .and. Lzd%Glr%geocode /='F') ii3=ii3+Lzd%Glr%d%n3i
                   if(ii3+1 > Lzd%Glr%d%n3i .and. Lzd%Glr%geocode /='F') &
                        ii3 = modulo(ii3+1,Lzd%Glr%d%n3i+1)
                   do i2=1,Lzd%Llr(ilr)%d%n2i
                       ii2 = i2 + Lzd%Llr(ilr)%nsi2 - 1
                       if(ii2 < 0 .and. Lzd%Glr%geocode =='P') ii2=ii2+Lzd%Glr%d%n2i
                       if(ii2+1 > Lzd%Glr%d%n2i .and. Lzd%Glr%geocode =='P') &
                            ii2 = modulo(ii2+1,Lzd%Glr%d%n2i+1)
                       do i1=1,Lzd%Llr(ilr)%d%n1i
                           ii1=i1 + Lzd%Llr(ilr)%nsi1-1
                           if(ii1<0 .and. Lzd%Glr%geocode /= 'F') ii1=ii1+Lzd%Glr%d%n1i
                           if(ii1+1 > Lzd%Glr%d%n1i.and.Lzd%Glr%geocode/='F') &
                                ii1 = modulo(ii1+1,Lzd%Glr%d%n1i+1)
                           ! indSmall is the index in the currect localization region
                           indSmall=indSmall+1
                           ! indLarge is the index in the whole box. 
                           indLarge=ii3*Lzd%Glr%d%n2i*Lzd%Glr%d%n1i +&
                               ii2*Lzd%Glr%d%n1i + ii1 + 1
                           rho(indLarge,ispin)=rho(indLarge,ispin)+rho_p(indSmall)
                       end do
                   end do
               end do
           end do
        end do
     else
        ind=ind+(Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*max(ncomplex,1)*npsir
     end if

     call f_free(rho_p)
     call f_free(psir)

     call deallocate_work_arrays_sumrho(w)
  end do orbitalsLoop
 
  call f_free(Lnscatterarr)
 

END SUBROUTINE local_partial_densityLinear


subroutine calculate_density_kernel(iproc, nproc, isKernel, orbs, orbs_tmb, &
           coeff, denskern, denskern_, keep_uncompressed_)
  use module_base
  use module_types
  use yaml_output
  use sparsematrix_base, only: sparse_matrix, sparsematrix_malloc_ptr, DENSE_FULL, assignment(=), &
                               sparsematrix_malloc, SPARSE_FULL
  use sparsematrix, only: compress_matrix, extract_taskgroup
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(orbitals_data),intent(in) :: orbs, orbs_tmb
  logical, intent(in) :: isKernel
  type(sparse_matrix), intent(inout) :: denskern
  real(kind=8),dimension(denskern%nfvctr,orbs%norb),intent(in):: coeff   !only use the first (occupied) orbitals
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

  if (communication_strategy==ALLGATHERV) then
      if (iproc==0) call yaml_map('communication strategy kernel','ALLGATHERV')
      stop 'calculate_density_kernel: ALLGATHERV option needs reworking due to the spin'
      call timing(iproc,'calc_kernel','ON')
      !if(iproc==0) write(*,'(1x,a)',advance='no') 'calculate density kernel... '
      density_kernel_partial=f_malloc((/orbs_tmb%norb,max(orbs_tmb%norbp,1)/), id='density_kernel_partial')
      fcoeff=f_malloc0((/orbs_tmb%norbp,orbs%norb/), id='fcoeff')
      if(orbs_tmb%norbp>0) then
          !decide whether we calculate the density kernel or just transformation matrix
          if(isKernel) then
             do iorb=1,orbs%norb
                !call daxpy(orbs_tmb%norbp,orbs%occup(iorb),coeff(1+orbs_tmb%isorb,iorb),1,fcoeff(1+orbs_tmb%isorb,iorb),1)
                do itmb=1,orbs_tmb%norbp
                     fcoeff(itmb,iorb) = orbs%occup(iorb)*coeff(orbs_tmb%isorb+itmb,iorb)
                end do
             end do
          else
             do iorb=1,orbs%norb
                do itmb=1,orbs_tmb%norbp
                     fcoeff(itmb,iorb) = coeff(orbs_tmb%isorb+itmb,iorb)
                end do
             end do
          end if

          call dgemm('n', 't', orbs_tmb%norb, orbs_tmb%norbp, orbs%norb, 1.d0, coeff(1,1), orbs_tmb%norb, &
               fcoeff(1,1), orbs_tmb%norbp, 0.d0, density_kernel_partial(1,1), orbs_tmb%norb)
      end if
      call f_free(fcoeff)
      call timing(iproc,'calc_kernel','OF')

      call timing(iproc,'waitAllgatKern','ON')
      call mpi_barrier(bigdft_mpi%mpi_comm,ierr)
      call timing(iproc,'waitAllgatKern','OF')

      !denskern_%matrix=f_malloc_ptr((/orbs_tmb%norb,orbs_tmb%norb/), id='denskern_%matrix')

      if (.not.keep_uncompressed) then
          denskern_%matrix=sparsematrix_malloc_ptr(denskern,iaction=DENSE_FULL,id='denskern_%matrix')
      end if

      if (nproc > 1) then
         call timing(iproc,'commun_kernel','ON')
         recvcounts=f_malloc((/0.to.nproc-1/),id='recvcounts')
         dspls=f_malloc((/0.to.nproc-1/),id='dspls')
         do jproc=0,nproc-1
             recvcounts(jproc)=orbs_tmb%norb*orbs_tmb%norb_par(jproc,0)
             dspls(jproc)=orbs_tmb%norb*orbs_tmb%isorb_par(jproc)
         end do
         sendcount=orbs_tmb%norb*orbs_tmb%norbp
         call mpi_allgatherv(density_kernel_partial(1,1), sendcount, mpi_double_precision, &
              denskern_%matrix(1,1,1), recvcounts, dspls, mpi_double_precision, &
              bigdft_mpi%mpi_comm, ierr)
         call f_free(recvcounts)
         call f_free(dspls)
         call timing(iproc,'commun_kernel','OF')
      else
         call vcopy(orbs_tmb%norb*orbs_tmb%norbp,density_kernel_partial(1,1),1,denskern_%matrix(1,1,1),1)
      end if

      call f_free(density_kernel_partial)

      call compress_matrix(iproc,denskern,inmat=denskern_%matrix,outmat=denskern_%matrix_compr)
      if (.not.keep_uncompressed) then
          call f_free_ptr(denskern_%matrix)
      end if
  else if (communication_strategy==ALLREDUCE) then
      if (iproc==0) call yaml_map('communication strategy kernel','ALLREDUCE')
      call timing(iproc,'calc_kernel','ON')
      !!if(iproc==0) write(*,'(1x,a)',advance='no') 'calculate density kernel... '
      !denskern_%matrix=f_malloc_ptr((/orbs_tmb%norb,orbs_tmb%norb/), id='denskern_%matrix_compr')
      if (.not.keep_uncompressed) then
          denskern_%matrix=sparsematrix_malloc_ptr(denskern,iaction=DENSE_FULL,id='denskern_%matrix')
      end if
      if(orbs%norbp>0) then
          fcoeff=f_malloc((/denskern%nfvctr,orbs%norbp/), id='fcoeff')
          !decide wether we calculate the density kernel or just transformation matrix
          if(isKernel)then
             do iorb=1,orbs%norbp
                !call f_zero(orbs_tmb%norb,f_coeff(1,iorb))
                !call daxpy(orbs_tmb%norb,orbs%occup(orbs%isorb+iorb),coeff(1,orbs%isorb+iorb),1,fcoeff(1,iorb),1)
                do itmb=1,denskern%nfvctr
                    fcoeff(itmb,iorb) = orbs%occup(orbs%isorb+iorb)*coeff(itmb,orbs%isorb+iorb)
                end do
             end do
          else
             do iorb=1,orbs%norbp
                call vcopy(denskern%nfvctr,coeff(1,orbs%isorb+iorb),1,fcoeff(1,iorb),1)
             end do
          end if
      !!if (iproc==0) then
      !!    do iorb=1,orbs%norbp
      !!        do jorb=1,denskern%nfvctr
      !!            write(970,'(a,2i9,f14.7)') 'iorb, jorb, fcoeff(jorb,iorb)', iorb, jorb, fcoeff(jorb,iorb)
      !!        end do
      !!    end do
      !!end if
          !call dgemm('n', 't', orbs_tmb%norb, orbs_tmb%norb, orbs%norbp, 1.d0, coeff(1,orbs%isorb+1), orbs_tmb%norb, &
          !     fcoeff(1,1), orbs_tmb%norb, 0.d0, denskern_%matrix(1,1,1), orbs_tmb%norb)
          call f_zero(denskern%nspin*denskern%nfvctr**2, denskern_%matrix(1,1,1))
          do iorb=1,orbs%norbp
              iiorb=orbs%isorb+iorb
              if (orbs%spinsgn(iiorb)>0.d0) then
                  ispin=1
              else
                  ispin=2
              end if
              call dgemm('n', 't', denskern%nfvctr, denskern%nfvctr, 1, 1.d0, coeff(1,orbs%isorb+iorb), denskern%nfvctr, &
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
      call mpi_barrier(bigdft_mpi%mpi_comm,ierr)
      call timing(iproc,'waitAllgatKern','OF')

      tmparr = sparsematrix_malloc(denskern,iaction=SPARSE_FULL,id='tmparr')
      call compress_matrix(iproc,denskern,inmat=denskern_%matrix,outmat=tmparr)
      if (keep_uncompressed) then
          if (nproc > 1) then
              call timing(iproc,'commun_kernel','ON')
              call mpiallred(denskern_%matrix(1,1,1), denskern%nspin*denskern%nfvctr**2, mpi_sum, bigdft_mpi%mpi_comm)
              call timing(iproc,'commun_kernel','OF')
          end if
      end if
      if (.not.keep_uncompressed) then
          call f_free_ptr(denskern_%matrix)
      end if
      if (nproc > 1) then
          call timing(iproc,'commun_kernel','ON')
          call mpiallred(tmparr(1), denskern%nspin*denskern%nvctr, mpi_sum, bigdft_mpi%mpi_comm)
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



!!subroutine calculate_density_kernel_uncompressed(iproc, nproc, isKernel, orbs, orbs_tmb, coeff, kernel)
!!  use module_base
!!  use module_types
!!  use yaml_output
!!  implicit none
!!
!!  ! Calling arguments
!!  integer,intent(in):: iproc, nproc
!!  type(orbitals_data),intent(in) :: orbs, orbs_tmb
!!  logical, intent(in) :: isKernel
!!  real(kind=8),dimension(orbs_tmb%norb,orbs%norb),intent(in):: coeff   !only use the first (occupied) orbitals
!!  real(kind=8),dimension(orbs_tmb%norb,orbs_tmb%norb),intent(out) :: kernel
!!
!!  ! Local variables
!!  integer :: istat, iall, ierr, sendcount, jproc, iorb, itmb
!!  real(kind=8),dimension(:,:),allocatable :: density_kernel_partial, fcoeff
!!! real(kind=8), dimension(:,:,), allocatable :: ks,ksk,ksksk
!!  character(len=*),parameter :: subname='calculate_density_kernel_uncompresses'
!!  integer,dimension(:),allocatable :: recvcounts, dspls
!!  integer,parameter :: ALLGATHERV=1, ALLREDUCE=2
!!  integer,parameter :: communication_strategy=ALLREDUCE
!!
!!  if (communication_strategy==ALLGATHERV) then
!!      if (iproc==0) call yaml_map('calculate density kernel, communication strategy','ALLGATHERV')
!!      call timing(iproc,'calc_kernel','ON')
!!      !if(iproc==0) write(*,'(1x,a)',advance='no') 'calculate density kernel... '
!!      density_kernel_partial = f_malloc((/ orbs_tmb%norb, max(orbs_tmb%norbp, 1) /),id='density_kernel_partial')
!!      fcoeff = f_malloc0((/ orbs_tmb%norb, orbs%norb /),id='fcoeff')
!!      if(orbs_tmb%norbp>0) then
!!          !decide whether we calculate the density kernel or just transformation matrix
!!          if(isKernel) then
!!             do iorb=1,orbs%norb
!!                !call daxpy(orbs_tmb%norbp,orbs%occup(iorb),coeff(1+orbs_tmb%isorb,iorb),1,fcoeff(1+orbs_tmb%isorb,iorb),1)
!!                do itmb=1,orbs_tmb%norbp
!!                     fcoeff(orbs_tmb%isorb+itmb,iorb) = orbs%occup(iorb)*coeff(orbs_tmb%isorb+itmb,iorb)
!!                end do
!!             end do
!!          else
!!             do iorb=1,orbs%norb
!!                do itmb=1,orbs_tmb%norbp
!!                     fcoeff(orbs_tmb%isorb+itmb,iorb) = coeff(orbs_tmb%isorb+itmb,iorb)
!!                end do
!!             end do
!!          end if
!!
!!          call dgemm('n', 't', orbs_tmb%norb, orbs_tmb%norbp, orbs%norb, 1.d0, coeff(1,1), orbs_tmb%norb, &
!!               fcoeff(orbs_tmb%isorb+1,1), orbs_tmb%norb, 0.d0, density_kernel_partial(1,1), orbs_tmb%norb)
!!      end if
!!      call f_free(fcoeff)
!!      call timing(iproc,'calc_kernel','OF')
!!
!!      call timing(iproc,'waitAllgatKern','ON')
!!      call mpi_barrier(bigdft_mpi%mpi_comm,ierr)
!!      call timing(iproc,'waitAllgatKern','OF')
!!
!!      if (nproc > 1) then
!!         call timing(iproc,'commun_kernel','ON')
!!         recvcounts = f_malloc(0.to.nproc-1,id='recvcounts')
!!         dspls = f_malloc(0.to.nproc-1,id='dspls')
!!         do jproc=0,nproc-1
!!             recvcounts(jproc)=orbs_tmb%norb*orbs_tmb%norb_par(jproc,0)
!!             dspls(jproc)=orbs_tmb%norb*orbs_tmb%isorb_par(jproc)
!!         end do
!!         sendcount=orbs_tmb%norb*orbs_tmb%norbp
!!         call mpi_allgatherv(density_kernel_partial(1,1), sendcount, mpi_double_precision, &
!!              kernel(1,1), recvcounts, dspls, mpi_double_precision, &
!!              bigdft_mpi%mpi_comm, ierr)
!!         call f_free(recvcounts)
!!         call f_free(dspls)
!!         call timing(iproc,'commun_kernel','OF')
!!      else
!!         call vcopy(orbs_tmb%norb*orbs_tmb%norbp,density_kernel_partial(1,1),1,kernel(1,1),1)
!!      end if
!!
!!      call f_free(density_kernel_partial)
!!  else if (communication_strategy==ALLREDUCE) then
!!      if (iproc==0) call yaml_map('calculate density kernel, communication strategy','ALLREDUCE')
!!      call timing(iproc,'calc_kernel','ON')
!!      !!if(iproc==0) write(*,'(1x,a)',advance='no') 'calculate density kernel... '
!!      if(orbs%norbp>0) then
!!          fcoeff = f_malloc0((/ orbs_tmb%norb, orbs%norb /),id='fcoeff')
!!
!!          !decide wether we calculate the density kernel or just transformation matrix
!!          if(isKernel)then
!!             do iorb=1,orbs%norbp
!!                !call daxpy(orbs_tmb%norb,orbs%occup(orbs%isorb+iorb),coeff(1,orbs%isorb+iorb),1,fcoeff(1,orbs%isorb+iorb),1)
!!                do itmb=1,orbs_tmb%norb
!!                     fcoeff(itmb,orbs%isorb+iorb) = orbs%occup(orbs%isorb+iorb)*coeff(itmb,orbs%isorb+iorb)
!!                end do
!!             end do
!!          else
!!             do iorb=1,orbs%norbp
!!                !call daxpy(orbs_tmb%norb,orbs%occup(orbs%isorb+iorb),coeff(1,orbs%isorb+iorb),1,fcoeff(1,orbs%isorb+iorb),1)
!!                do itmb=1,orbs_tmb%norb
!!                     fcoeff(itmb,orbs%isorb+iorb) = coeff(itmb,orbs%isorb+iorb)
!!                end do
!!          end do
!!
!!          end if
!!          call dgemm('n', 't', orbs_tmb%norb, orbs_tmb%norb, orbs%norbp, 1.d0, coeff(1,orbs%isorb+1), orbs_tmb%norb, &
!!               fcoeff(1,orbs%isorb+1), orbs_tmb%norb, 0.d0, kernel(1,1), orbs_tmb%norb)
!!          call f_free(fcoeff)
!!      else
!!          call f_zero(orbs_tmb%norb**2, kernel(1,1))
!!      end if
!!      call timing(iproc,'calc_kernel','OF')
!!
!!      call timing(iproc,'waitAllgatKern','ON')
!!      call mpi_barrier(bigdft_mpi%mpi_comm,ierr)
!!      call timing(iproc,'waitAllgatKern','OF')
!!      if (nproc > 1) then
!!          call timing(iproc,'commun_kernel','ON')
!!          call mpiallred(kernel(1,1),orbs_tmb%norb**2, mpi_sum, bigdft_mpi%mpi_comm)
!!          call timing(iproc,'commun_kernel','OF')
!!      end if
!!  end if
!!
!! ! Purify Kernel
!! !call dgemm('n', 't', orbs_tmb%norb, orbs_tmb%norb, orbs_tmb%norbp, 1.d0, kernel(1,orbs_tmb%isorb+1), orbs_tmb%norb, &
!! !           overlap(1,orbs_tmb%isorb+1), orbs_tmb%norb, 0.d0, ks(1,1), orbs_tmb%norb) 
!! !call dgemm('n', 't', orbs_tmb%norb, orbs_tmb%norb, orbs_tmb%norbp, 1.d0, ks(1,orbs_tmb%isorb+1), orbs_tmb%norb, &
!! !           kernel(1,orbs_tmb%isorb+1), orbs_tmb%norb, 0.d0, ksk(1,1), orbs_tmb%norb)
!! !call dgemm('n', 't', orbs_tmb%norb, orbs_tmb%norb, orbs_tmb%norbp, 1.d0, ks(1,orbs_tmb%isorb+1), orbs_tmb%norb, &
!! !           ksk(1,orbs_tmb%isorb+1), orbs_tmb%norb, 0.d0, ksksk(1,1), orbs_tmb%norb)
!!
!!
!! !!if(present(overlap)) then
!!   !!allocate(ks(orbs_tmb%norb,orbs_tmb%norb),stat=istat)
!!   !!call memocc(istat, ks, 'ks', subname) 
!!   !!allocate(ksk(orbs_tmb%norb,orbs_tmb%norb),stat=istat)
!!   !!call memocc(istat, ksk, 'ksk', subname) 
!!   !!allocate(ksksk(orbs_tmb%norb,orbs_tmb%norb),stat=istat)
!!   !!call memocc(istat, ksksk, 'ksksk', subname) 
!!
!!   !!call dgemm('n', 't', orbs_tmb%norb, orbs_tmb%norb, orbs_tmb%norb, 1.d0, kernel(1,1), orbs_tmb%norb, &
!!   !!           overlap(1,1), orbs_tmb%norb, 0.d0, ks(1,1), orbs_tmb%norb) 
!!   !!call dgemm('n', 't', orbs_tmb%norb, orbs_tmb%norb, orbs_tmb%norb, 1.d0, ks(1,1), orbs_tmb%norb, &
!!   !!           kernel(1,1), orbs_tmb%norb, 0.d0, ksk(1,1), orbs_tmb%norb)
!!   !!call dgemm('n', 't', orbs_tmb%norb, orbs_tmb%norb, orbs_tmb%norb, 1.d0, ks(1,1), orbs_tmb%norb, &
!!   !!           ksk(1,1), orbs_tmb%norb, 0.d0, ksksk(1,1), orbs_tmb%norb)
!!   !!print *,'PURIFYING THE KERNEL'
!!   !!kernel = 3*ksk-2*ksksk
!!   !!
!!   !!iall = -product(shape(ks))*kind(ks)
!!   !!deallocate(ks,stat=istat)
!!   !!call memocc(istat, iall, 'ks', subname)
!!   !!iall = -product(shape(ksk))*kind(ksk)
!!   !!deallocate(ksk,stat=istat)
!!   !!call memocc(istat, iall, 'ksk', subname)
!!   !!iall = -product(shape(ksksk))*kind(ksksk)
!!   !!deallocate(ksksk,stat=istat)
!!   !!call memocc(istat, iall, 'ksksk', subname)
!! !!end if
!!
!!end subroutine calculate_density_kernel_uncompressed





subroutine sumrho_for_TMBs(iproc, nproc, hx, hy, hz, collcom_sr, denskern, denskern_, ndimrho, rho, rho_negative, &
        print_results)
  use module_base
  use module_types
  use yaml_output
  use sparsematrix_base, only: sparse_matrix
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, ndimrho
  real(kind=8),intent(in) :: hx, hy, hz
  type(comms_linear),intent(in) :: collcom_sr
  type(sparse_matrix),intent(in) :: denskern
  type(matrices),intent(in) :: denskern_
  real(kind=8),dimension(ndimrho),intent(out) :: rho
  logical,intent(out) :: rho_negative
  logical,intent(in),optional :: print_results

  ! Local variables
  integer :: ipt, ii, i0, iiorb, jjorb, istat, iall, i, j, ierr, ind, ispin, ishift, ishift_mat, iorb_shift
  real(8) :: tt, total_charge, hxh, hyh, hzh, factor, tt1, rho_neg
  integer,dimension(:),allocatable :: isend_total
  real(kind=8),dimension(:),allocatable :: rho_local
  character(len=*),parameter :: subname='sumrho_for_TMBs'
  logical :: print_local
  integer :: size_of_double, info, mpisource, istsource, istdest, nsize, jproc, ishift_dest, ishift_source

  call f_routine('sumrho_for_TMBs')

  ! check whether all entries of the charge density are positive
  rho_negative=.false.

  if (present(print_results)) then
      if (print_results) then
          print_local=.true.
      else
          print_local=.false.
      end if
  else
      print_local=.true.
  end if


  rho_local = f_malloc(collcom_sr%nptsp_c*denskern%nspin,id='rho_local')

  ! Define some constant factors.
  hxh=.5d0*hx
  hyh=.5d0*hy
  hzh=.5d0*hz
  factor=1.d0/(hxh*hyh*hzh)

  call timing(iproc,'sumrho_TMB    ','ON')
  
  ! Initialize rho. (not necessary for the moment)
  !if (xc_isgga()) then
  !    call f_zero(collcom_sr%nptsp_c, rho_local)
  !else
   !   ! There is no mpi_allreduce, therefore directly initialize to
   !   ! 10^-20 and not 10^-20/nproc.
  !    rho_local=1.d-20
  !end if


  !!if (print_local .and. iproc==0) write(*,'(a)', advance='no') 'Calculating charge density... '

  total_charge=0.d0
  rho_neg=0.d0

!ispin=1
  !SM: check if the modulo operations take a lot of time. If so, try to use an
  !auxiliary array with shifted bounds in order to access smat%matrixindex_in_compressed_fortransposed
  do ispin=1,denskern%nspin
      if (ispin==1) then
          ishift=0
      else
          ishift=collcom_sr%ndimind_c/2
      end if
      iorb_shift=(ispin-1)*denskern%nfvctr
      ishift_mat=(ispin-1)*denskern%nvctr
      !$omp parallel default(private) &
      !$omp shared(total_charge, collcom_sr, factor, denskern, denskern_, rho_local, rho_neg, ispin, ishift, ishift_mat, iorb_shift)
      !$omp do schedule(static,200) reduction(+:total_charge, rho_neg)
      do ipt=1,collcom_sr%nptsp_c
          ii=collcom_sr%norb_per_gridpoint_c(ipt)
    
          i0=collcom_sr%isptsp_c(ipt)+ishift
          tt=1.e-20_dp
          do i=1,ii
              iiorb=collcom_sr%indexrecvorbital_c(i0+i) - iorb_shift
              !iiorb=mod(iiorb-1,denskern%nfvctr)+1
    !ispin=spinsgn(iiorb) 
              tt1=collcom_sr%psit_c(i0+i)
              ind=denskern%matrixindex_in_compressed_fortransposed(iiorb,iiorb)
              ind=ind+ishift_mat-denskern%isvctrp_tg
              tt=tt+denskern_%matrix_compr(ind)*tt1*tt1
    !tt(ispin)=tt(ispin)+denskern_%matrix_compr(ind)*tt1*tt1
              do j=i+1,ii
                  jjorb=collcom_sr%indexrecvorbital_c(i0+j) - iorb_shift
                  !jjorb=mod(jjorb-1,denskern%nfvctr)+1
                  ind=denskern%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
                  if (ind==0) cycle
                  ind=ind+ishift_mat-denskern%isvctrp_tg
                  tt=tt+2.0_dp*denskern_%matrix_compr(ind)*tt1*collcom_sr%psit_c(i0+j)
              end do
          end do
          tt=factor*tt
          total_charge=total_charge+tt
          rho_local(ipt+(ispin-1)*collcom_sr%nptsp_c)=tt
    !rho_local(ipt,ispin)=tt(ispin)
          if (tt<0.d0) rho_neg=rho_neg+1.d0
      end do
      !$omp end do
      !$omp end parallel
  end do


  !if (print_local .and. iproc==0) write(*,'(a)') 'done.'

  call timing(iproc,'sumrho_TMB    ','OF')


  call communicate_density()


  !!if (nproc > 1) then
  !!   call mpiallred(irho, 1, mpi_sum, bigdft_mpi%mpi_comm)
  !!end if

  if (rho_neg>0.d0) then
      rho_negative=.true.
  end if

  call f_free(rho_local)

  call f_release_routine()

  contains

    subroutine communicate_density()
      implicit none
      real(kind=8),dimension(2) :: reducearr

      call f_routine(id='communicate_density')
      call timing(iproc,'sumrho_allred','ON')
    
      ! Communicate the density to meet the shape required by the Poisson solver.
      !!if (nproc>1) then
      !!    call mpi_alltoallv(rho_local, collcom_sr%nsendcounts_repartitionrho, collcom_sr%nsenddspls_repartitionrho, &
      !!                       mpi_double_precision, rho, collcom_sr%nrecvcounts_repartitionrho, &
      !!                       collcom_sr%nrecvdspls_repartitionrho, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
      !!else
      !!    call vcopy(ndimrho, rho_local(1), 1, rho(1), 1)
      !!end if
    
      !!!!do ierr=1,size(rho)
      !!!!    write(200+iproc,*) ierr, rho(ierr)
      !!!!end do
    
    
    
      if (nproc>1) then
          call mpi_type_size(mpi_double_precision, size_of_double, ierr)
          call mpi_info_create(info, ierr)
          call mpi_info_set(info, "no_locks", "true", ierr)
          call mpi_win_create(rho_local(1), int(collcom_sr%nptsp_c*denskern%nspin*size_of_double,kind=mpi_address_kind), &
               size_of_double, info, bigdft_mpi%mpi_comm, collcom_sr%window, ierr)
          call mpi_info_free(info, ierr)
    
          call mpi_win_fence(mpi_mode_noprecede, collcom_sr%window, ierr)
    
          ! This is a bit quick and dirty. Could be done in a better way, but
          ! would probably required to pass additional arguments to the subroutine
          isend_total = f_malloc0(0.to.nproc-1,id='isend_total')
          isend_total(iproc)=collcom_sr%nptsp_c
          call mpiallred(isend_total(0), nproc, mpi_sum, bigdft_mpi%mpi_comm)
    
    
          do ispin=1,denskern%nspin
              !ishift_dest=(ispin-1)*sum(collcom_sr%commarr_repartitionrho(4,:)) !spin shift for the receive buffer
              ishift_dest=(ispin-1)*ndimrho/denskern%nspin
              do jproc=1,collcom_sr%ncomms_repartitionrho
                  mpisource=collcom_sr%commarr_repartitionrho(1,jproc)
                  istsource=collcom_sr%commarr_repartitionrho(2,jproc)
                  istdest=collcom_sr%commarr_repartitionrho(3,jproc)
                  nsize=collcom_sr%commarr_repartitionrho(4,jproc)
                  ishift_source=(ispin-1)*isend_total(mpisource) !spin shift for the send buffer
                  if (nsize>0) then
                      !!write(*,'(6(a,i0))') 'process ',iproc, ' gets ',nsize,' elements at position ',istdest+ishift_dest, &
                      !!                     ' from position ',istsource+ishift_source,' on process ',mpisource, &
                      !!                     '; error code=',ierr
                      call mpi_get(rho(istdest+ishift_dest), nsize, mpi_double_precision, mpisource, &
                           int((istsource-1+ishift_source),kind=mpi_address_kind), &
                           nsize, mpi_double_precision, collcom_sr%window, ierr)
                  end if
              end do
          end do
          call mpi_win_fence(0, collcom_sr%window, ierr)
          !!write(*,'(a,i0)') 'mpi_win_fence error code: ',ierr
          call mpi_win_free(collcom_sr%window, ierr)
          !!write(*,'(a,i0)') 'mpi_win_free error code: ',ierr
    
          call f_free(isend_total)
      else
          call vcopy(ndimrho, rho_local(1), 1, rho(1), 1)
      end if
    
      !do ierr=1,size(rho)
      !    write(300+iproc,*) ierr, rho(ierr)
      !end do
      !call mpi_finalize(ierr)
      !stop
    
      if (nproc > 1) then
          reducearr(1) = total_charge
          reducearr(2) = rho_neg
          call mpiallred(reducearr(1), 2, mpi_sum, bigdft_mpi%mpi_comm)
          total_charge = reducearr(1)
          rho_neg = reducearr(2)
         !call mpiallred(total_charge, 1, mpi_sum, bigdft_mpi%mpi_comm)
      end if
    
      !!if(print_local .and. iproc==0) write(*,'(3x,a,es20.12)') 'Calculation finished. TOTAL CHARGE = ', total_charge*hxh*hyh*hzh
      if (iproc==0 .and. print_local) then
          call yaml_map('Total charge',total_charge*hxh*hyh*hzh,fmt='(es20.12)')
      end if
      
      call timing(iproc,'sumrho_allred','OF')
      call f_release_routine()

    end subroutine communicate_density

  !!write(*,*) 'after deallocate'
  !!call mpi_finalize(ierr)
  !!stop


end subroutine sumrho_for_TMBs

!!$!> this routine is important to create the density needed for GGA XC functionals calculation
!!$subroutine fill_global_density(rho_local)
!!$  use module_base
!!$  implicit none
!!$  
!!$  real(dp), dimension(:,:,:,:), allocatable :: rho_global !>density in the global box to be reduced 
!!$
!!$  call f_malloc_routine_id('fill_global_density')
!!$
!!$  !malloc and initialization to zero to be defined
!!$  rho_global=f_malloc0((/dpbox%ndims(1),dpbox%ndims(2),dpbox%ndims(3),1/),id='rho_global')
!!$
!!$  !rho is the density in the LDA_like poisson solver distribution, therefore starting from i3s+i3xcsh until i3s+i3xcsh+n3p-1
!!$  do jproc=0,nproc-1
!!$     !quantity of data to be copied
!!$     nrho=collcom_sr%nsendcounts_repartitionrho(jproc)
!!$     !starting point of the local array to be put at the place of 
!!$     irls=collcom_sr%nsenddspls_repartitionrho(jproc)+1
!!$     !starting point of the global density in the processor jproc in z direction
!!$     isz=dpbox%nscatterarr(jproc,3)
!!$     irgs=
!!$     call vcopy(nrho,rho_local(irhos),1,rho_global(),1)
!!$  end do
!!$
!!$  call f_free(rho_global)
!!$  !this call might free all the allocatable arrays in the routine if the address of the original object is known
!!$  call f_malloc_free_routine()
!!$
!!$end subroutine fill_global_density




subroutine check_negative_rho(ndimrho, rho, rho_negative)
  use module_base
  use dynamic_memory
  implicit none

  ! Calling arguments
  integer,intent(in) :: ndimrho
  real(kind=8),dimension(ndimrho),intent(in) :: rho
  logical,intent(out) :: rho_negative

  ! Local variables
  integer :: i, irho, ierr

  call f_routine(id='check_negative_rho')

  irho=0
  do i=1,ndimrho
      if (rho(i)<0.d0) then
          irho=1
          exit
      end if
  end do

  if (bigdft_mpi%nproc > 1) then
     call mpiallred(irho, 1, mpi_sum, bigdft_mpi%mpi_comm)
  end if

  if (irho>0) then
      rho_negative=.true.
  else
      rho_negative=.false.
  end if

  call f_release_routine()

end subroutine check_negative_rho
