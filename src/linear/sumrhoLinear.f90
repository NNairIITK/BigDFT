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


     call initialize_work_arrays_sumrho(Lzd%Llr(ilr),w)
     rho_p = f_malloc(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*nspinn,id='rho_p')
     psir = f_malloc((/ Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i, npsir /),id='psir')
  
     if (Lzd%Llr(ilr)%geocode == 'F') then
        call to_zero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*npsir,psir)
     end if
 
     !Need to zero rho_p
     call to_zero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*nspinn, rho_p)

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


subroutine calculate_density_kernel(iproc, nproc, isKernel, orbs, orbs_tmb, coeff, denskern, denskern_)
  use module_base
  use module_types
  use yaml_output
  use sparsematrix_base, only: sparse_matrix
  use sparsematrix, only: compress_matrix
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(orbitals_data),intent(in) :: orbs, orbs_tmb
  logical, intent(in) :: isKernel
  real(kind=8),dimension(orbs_tmb%norb,orbs%norb),intent(in):: coeff   !only use the first (occupied) orbitals
  type(sparse_matrix), intent(inout) :: denskern
  type(matrices), intent(out) :: denskern_

  ! Local variables
  integer :: ierr, sendcount, jproc, iorb, itmb
  real(kind=8),dimension(:,:),allocatable :: density_kernel_partial, fcoeff
! real(kind=8), dimension(:,:,), allocatable :: ks,ksk,ksksk
  character(len=*),parameter :: subname='calculate_density_kernel'
  integer,dimension(:),allocatable :: recvcounts, dspls
  integer,parameter :: ALLGATHERV=1, ALLREDUCE=2
  integer,parameter :: communication_strategy=ALLREDUCE

  call f_routine(id='calculate_density_kernel')

  if (communication_strategy==ALLGATHERV) then
      if (iproc==0) call yaml_map('communication strategy kernel','ALLGATHERV')
      call timing(iproc,'calc_kernel','ON') !lr408t
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
      call timing(iproc,'calc_kernel','OF') !lr408t

      call timing(iproc,'waitAllgatKern','ON')
      call mpi_barrier(bigdft_mpi%mpi_comm,ierr)
      call timing(iproc,'waitAllgatKern','OF')

      denskern_%matrix=f_malloc_ptr((/orbs_tmb%norb,orbs_tmb%norb/), id='denskern_%matrix')

      if (nproc > 1) then
         call timing(iproc,'commun_kernel','ON') !lr408t
         recvcounts=f_malloc((/0.to.nproc-1/),id='recvcounts')
         dspls=f_malloc((/0.to.nproc-1/),id='dspls')
         do jproc=0,nproc-1
             recvcounts(jproc)=orbs_tmb%norb*orbs_tmb%norb_par(jproc,0)
             dspls(jproc)=orbs_tmb%norb*orbs_tmb%isorb_par(jproc)
         end do
         sendcount=orbs_tmb%norb*orbs_tmb%norbp
         call mpi_allgatherv(density_kernel_partial(1,1), sendcount, mpi_double_precision, &
              denskern_%matrix(1,1), recvcounts, dspls, mpi_double_precision, &
              bigdft_mpi%mpi_comm, ierr)
         call f_free(recvcounts)
         call f_free(dspls)
         call timing(iproc,'commun_kernel','OF') !lr408t
      else
         call vcopy(orbs_tmb%norb*orbs_tmb%norbp,density_kernel_partial(1,1),1,denskern_%matrix(1,1),1)
      end if

      call f_free(density_kernel_partial)

      call compress_matrix(iproc,denskern,inmat=denskern_%matrix,outmat=denskern_%matrix_compr)
      call f_free_ptr(denskern_%matrix)
  else if (communication_strategy==ALLREDUCE) then
      if (iproc==0) call yaml_map('communication strategy kernel','ALLREDUCE')
      call timing(iproc,'calc_kernel','ON') !lr408t
      !!if(iproc==0) write(*,'(1x,a)',advance='no') 'calculate density kernel... '
      denskern_%matrix=f_malloc_ptr((/orbs_tmb%norb,orbs_tmb%norb/), id='denskern_%matrix_compr')
      if(orbs%norbp>0) then
          fcoeff=f_malloc((/orbs_tmb%norb,orbs%norbp/), id='fcoeff')
          !decide wether we calculate the density kernel or just transformation matrix
          if(isKernel)then
             do iorb=1,orbs%norbp
                !call to_zero(orbs_tmb%norb,f_coeff(1,iorb))
                !call daxpy(orbs_tmb%norb,orbs%occup(orbs%isorb+iorb),coeff(1,orbs%isorb+iorb),1,fcoeff(1,iorb),1)
                do itmb=1,orbs_tmb%norb
                    fcoeff(itmb,iorb) = orbs%occup(orbs%isorb+iorb)*coeff(itmb,orbs%isorb+iorb)
                end do
             end do
          else
             do iorb=1,orbs%norbp
                call vcopy(orbs_tmb%norb,coeff(1,orbs%isorb+iorb),1,fcoeff(1,iorb),1)
             end do
          end if
          call dgemm('n', 't', orbs_tmb%norb, orbs_tmb%norb, orbs%norbp, 1.d0, coeff(1,orbs%isorb+1), orbs_tmb%norb, &
               fcoeff(1,1), orbs_tmb%norb, 0.d0, denskern_%matrix(1,1), orbs_tmb%norb)
          call f_free(fcoeff)
      else
          call to_zero(orbs_tmb%norb**2, denskern_%matrix(1,1))
      end if
      call timing(iproc,'calc_kernel','OF') !lr408t

      call timing(iproc,'waitAllgatKern','ON')
      call mpi_barrier(bigdft_mpi%mpi_comm,ierr)
      call timing(iproc,'waitAllgatKern','OF')

      call compress_matrix(iproc,denskern,inmat=denskern_%matrix,outmat=denskern_%matrix_compr)
      call f_free_ptr(denskern_%matrix)
      if (nproc > 1) then
          call timing(iproc,'commun_kernel','ON') !lr408t
          call mpiallred(denskern_%matrix_compr(1), denskern%nvctr, mpi_sum, bigdft_mpi%mpi_comm)
          call timing(iproc,'commun_kernel','OF') !lr408t
      end if

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



subroutine calculate_density_kernel_uncompressed(iproc, nproc, isKernel, orbs, orbs_tmb, coeff, kernel)
  use module_base
  use module_types
  use yaml_output
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(orbitals_data),intent(in) :: orbs, orbs_tmb
  logical, intent(in) :: isKernel
  real(kind=8),dimension(orbs_tmb%norb,orbs%norb),intent(in):: coeff   !only use the first (occupied) orbitals
  real(kind=8),dimension(orbs_tmb%norb,orbs_tmb%norb),intent(out) :: kernel

  ! Local variables
  integer :: istat, iall, ierr, sendcount, jproc, iorb, itmb
  real(kind=8),dimension(:,:),allocatable :: density_kernel_partial, fcoeff
! real(kind=8), dimension(:,:,), allocatable :: ks,ksk,ksksk
  character(len=*),parameter :: subname='calculate_density_kernel'
  integer,dimension(:),allocatable :: recvcounts, dspls
  integer,parameter :: ALLGATHERV=1, ALLREDUCE=2
  integer,parameter :: communication_strategy=ALLREDUCE

  if (communication_strategy==ALLGATHERV) then
      if (iproc==0) call yaml_map('calculate density kernel, communication strategy','ALLGATHERV')
      call timing(iproc,'calc_kernel','ON') !lr408t
      !if(iproc==0) write(*,'(1x,a)',advance='no') 'calculate density kernel... '
      density_kernel_partial = f_malloc((/ orbs_tmb%norb, max(orbs_tmb%norbp, 1) /),id='density_kernel_partial')
      fcoeff = f_malloc0((/ orbs_tmb%norb, orbs%norb /),id='fcoeff')
      if(orbs_tmb%norbp>0) then
          !decide whether we calculate the density kernel or just transformation matrix
          if(isKernel) then
             do iorb=1,orbs%norb
                !call daxpy(orbs_tmb%norbp,orbs%occup(iorb),coeff(1+orbs_tmb%isorb,iorb),1,fcoeff(1+orbs_tmb%isorb,iorb),1)
                do itmb=1,orbs_tmb%norbp
                     fcoeff(orbs_tmb%isorb+itmb,iorb) = orbs%occup(iorb)*coeff(orbs_tmb%isorb+itmb,iorb)
                end do
             end do
          else
             do iorb=1,orbs%norb
                do itmb=1,orbs_tmb%norbp
                     fcoeff(orbs_tmb%isorb+itmb,iorb) = coeff(orbs_tmb%isorb+itmb,iorb)
                end do
             end do
          end if

          call dgemm('n', 't', orbs_tmb%norb, orbs_tmb%norbp, orbs%norb, 1.d0, coeff(1,1), orbs_tmb%norb, &
               fcoeff(orbs_tmb%isorb+1,1), orbs_tmb%norb, 0.d0, density_kernel_partial(1,1), orbs_tmb%norb)
      end if
      call f_free(fcoeff)
      call timing(iproc,'calc_kernel','OF') !lr408t

      call timing(iproc,'waitAllgatKern','ON')
      call mpi_barrier(bigdft_mpi%mpi_comm,ierr)
      call timing(iproc,'waitAllgatKern','OF')

      if (nproc > 1) then
         call timing(iproc,'commun_kernel','ON') !lr408t
         recvcounts = f_malloc(0.to.nproc-1,id='recvcounts')
         dspls = f_malloc(0.to.nproc-1,id='dspls')
         do jproc=0,nproc-1
             recvcounts(jproc)=orbs_tmb%norb*orbs_tmb%norb_par(jproc,0)
             dspls(jproc)=orbs_tmb%norb*orbs_tmb%isorb_par(jproc)
         end do
         sendcount=orbs_tmb%norb*orbs_tmb%norbp
         call mpi_allgatherv(density_kernel_partial(1,1), sendcount, mpi_double_precision, &
              kernel(1,1), recvcounts, dspls, mpi_double_precision, &
              bigdft_mpi%mpi_comm, ierr)
         call f_free(recvcounts)
         call f_free(dspls)
         call timing(iproc,'commun_kernel','OF') !lr408t
      else
         call vcopy(orbs_tmb%norb*orbs_tmb%norbp,density_kernel_partial(1,1),1,kernel(1,1),1)
      end if

      call f_free(density_kernel_partial)
  else if (communication_strategy==ALLREDUCE) then
      if (iproc==0) call yaml_map('calculate density kernel, communication strategy','ALLREDUCE')
      call timing(iproc,'calc_kernel','ON') !lr408t
      !!if(iproc==0) write(*,'(1x,a)',advance='no') 'calculate density kernel... '
      if(orbs%norbp>0) then
          fcoeff = f_malloc0((/ orbs_tmb%norb, orbs%norb /),id='fcoeff')

          !decide wether we calculate the density kernel or just transformation matrix
          if(isKernel)then
             do iorb=1,orbs%norbp
                !call daxpy(orbs_tmb%norb,orbs%occup(orbs%isorb+iorb),coeff(1,orbs%isorb+iorb),1,fcoeff(1,orbs%isorb+iorb),1)
                do itmb=1,orbs_tmb%norb
                     fcoeff(itmb,orbs%isorb+iorb) = orbs%occup(orbs%isorb+iorb)*coeff(itmb,orbs%isorb+iorb)
                end do
             end do
          else
             do iorb=1,orbs%norbp
                !call daxpy(orbs_tmb%norb,orbs%occup(orbs%isorb+iorb),coeff(1,orbs%isorb+iorb),1,fcoeff(1,orbs%isorb+iorb),1)
                do itmb=1,orbs_tmb%norb
                     fcoeff(itmb,orbs%isorb+iorb) = coeff(itmb,orbs%isorb+iorb)
                end do
          end do

          end if
          call dgemm('n', 't', orbs_tmb%norb, orbs_tmb%norb, orbs%norbp, 1.d0, coeff(1,orbs%isorb+1), orbs_tmb%norb, &
               fcoeff(1,orbs%isorb+1), orbs_tmb%norb, 0.d0, kernel(1,1), orbs_tmb%norb)
          call f_free(fcoeff)
      else
          call to_zero(orbs_tmb%norb**2, kernel(1,1))
      end if
      call timing(iproc,'calc_kernel','OF') !lr408t

      call timing(iproc,'waitAllgatKern','ON')
      call mpi_barrier(bigdft_mpi%mpi_comm,ierr)
      call timing(iproc,'waitAllgatKern','OF')
      if (nproc > 1) then
          call timing(iproc,'commun_kernel','ON') !lr408t
          call mpiallred(kernel(1,1),orbs_tmb%norb**2, mpi_sum, bigdft_mpi%mpi_comm)
          call timing(iproc,'commun_kernel','OF') !lr408t
      end if
  end if

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

end subroutine calculate_density_kernel_uncompressed





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
  integer :: ipt, ii, i0, iiorb, jjorb, istat, iall, i, j, ierr, ind
  real(8) :: tt, total_charge, hxh, hyh, hzh, factor, tt1
  real(kind=8),dimension(:),allocatable :: rho_local
  character(len=*),parameter :: subname='sumrho_for_TMBs'
  logical :: print_local
  integer :: size_of_double, info, mpisource, istsource, istdest, nsize, jproc, irho

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


  rho_local = f_malloc(collcom_sr%nptsp_c,id='rho_local')

  ! Define some constant factors.
  hxh=.5d0*hx
  hyh=.5d0*hy
  hzh=.5d0*hz
  factor=1.d0/(hxh*hyh*hzh)

  call timing(iproc,'sumrho_TMB    ','ON')
  
  ! Initialize rho. (not necessary for the moment)
  !if (xc_isgga()) then
  !    call to_zero(collcom_sr%nptsp_c, rho_local)
  !else
   !   ! There is no mpi_allreduce, therefore directly initialize to
   !   ! 10^-20 and not 10^-20/nproc.
  !    rho_local=1.d-20
  !end if


  !!if (print_local .and. iproc==0) write(*,'(a)', advance='no') 'Calculating charge density... '

  total_charge=0.d0
  irho=0

!ispin=1
  !$omp parallel default(private) &
  !$omp shared(total_charge, collcom_sr, factor, denskern, denskern_, rho_local, irho)
  !$omp do schedule(static,50) reduction(+:total_charge, irho)
  do ipt=1,collcom_sr%nptsp_c
      ii=collcom_sr%norb_per_gridpoint_c(ipt)

      i0=collcom_sr%isptsp_c(ipt)
      tt=1.e-20_dp
      do i=1,ii
          iiorb=collcom_sr%indexrecvorbital_c(i0+i)
!ispin=spinsgn(iiorb) 
          tt1=collcom_sr%psit_c(i0+i)
          ind=denskern%matrixindex_in_compressed_fortransposed(iiorb,iiorb)
          tt=tt+denskern_%matrix_compr(ind)*tt1*tt1
!tt(ispin)=tt(ispin)+denskern_%matrix_compr(ind)*tt1*tt1
          do j=i+1,ii
              jjorb=collcom_sr%indexrecvorbital_c(i0+j)
              ind=denskern%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
              if (ind==0) cycle
              tt=tt+2.0_dp*denskern_%matrix_compr(ind)*tt1*collcom_sr%psit_c(i0+j)
          end do
      end do
      tt=factor*tt
      total_charge=total_charge+tt
      rho_local(ipt)=tt
!rho_local(ipt,ispin)=tt(ispin)
      if (tt<0.d0) irho=irho+1
  end do
  !$omp end do
  !$omp end parallel

  if (nproc > 1) then
     call mpiallred(irho, 1, mpi_sum, bigdft_mpi%mpi_comm)
  end if

  if (irho>0) then
      rho_negative=.true.
  end if

  !if (print_local .and. iproc==0) write(*,'(a)') 'done.'

  call timing(iproc,'sumrho_TMB    ','OF')

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
      call mpi_win_create(rho_local(1), int(collcom_sr%nptsp_c*size_of_double,kind=mpi_address_kind), size_of_double, &
           info, bigdft_mpi%mpi_comm, collcom_sr%window, ierr)
      call mpi_info_free(info, ierr)

      call mpi_win_fence(mpi_mode_noprecede, collcom_sr%window, ierr)

      do jproc=1,collcom_sr%ncomms_repartitionrho
          mpisource=collcom_sr%commarr_repartitionrho(1,jproc)
          istsource=collcom_sr%commarr_repartitionrho(2,jproc)
          istdest=collcom_sr%commarr_repartitionrho(3,jproc)
          nsize=collcom_sr%commarr_repartitionrho(4,jproc)
          if (nsize>0) then
              call mpi_get(rho(istdest), nsize, mpi_double_precision, mpisource, &
                   int((istsource-1),kind=mpi_address_kind), &
                   nsize, mpi_double_precision, collcom_sr%window, ierr)
              !!write(*,'(6(a,i0))') 'process ',iproc, ' gets ',nsize,' elements at position ',istdest, &
              !!                     ' from position ',istsource,' on process ',mpisource, &
              !!                     '; error code=',ierr
          end if
      end do
      call mpi_win_fence(0, collcom_sr%window, ierr)
      !!write(*,'(a,i0)') 'mpi_win_fence error code: ',ierr
      call mpi_win_free(collcom_sr%window, ierr)
      !!write(*,'(a,i0)') 'mpi_win_free error code: ',ierr
  else
      call vcopy(ndimrho, rho_local(1), 1, rho(1), 1)
  end if

  !do ierr=1,size(rho)
  !    write(300+iproc,*) ierr, rho(ierr)
  !end do
  !call mpi_finalize(ierr)
  !stop

  if (nproc > 1) then
     call mpiallred(total_charge, 1, mpi_sum, bigdft_mpi%mpi_comm)
  end if

  !!if(print_local .and. iproc==0) write(*,'(3x,a,es20.12)') 'Calculation finished. TOTAL CHARGE = ', total_charge*hxh*hyh*hzh
  if (iproc==0 .and. print_local) then
      call yaml_map('Total charge',total_charge*hxh*hyh*hzh,fmt='(es20.12)')
  end if
  
  call timing(iproc,'sumrho_allred','OF')

  call f_free(rho_local)

  !!write(*,*) 'after deallocate'
  !!call mpi_finalize(ierr)
  !!stop

  call f_release_routine()

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

!> perform the communication needed for the potential and verify that the results is as expected
subroutine check_communication_potential(iproc,denspot,tmb)
  use module_base, only:dp,bigdft_mpi,mpi_sum,mpi_max,mpiallred
  use module_types
  use module_interfaces
  use yaml_output
  use dictionaries, only: f_err_throw
  use communications, only: start_onesided_communication
  implicit none
  integer,intent(in) :: iproc
  type(DFT_wavefunction), intent(inout) :: tmb
  type(DFT_local_fields), intent(inout) :: denspot
  !local variables
  logical :: dosome
  integer :: i1,i2,i3,ind,i3s,n3p,ilr,iorb,ilr_orb,n2i,n1i,ierr,numtot,i_stat,i_all
  real(dp) :: maxdiff,sumdiff,testval
  real(dp),parameter :: tol_calculation_mean=1.d-12
  real(dp),parameter :: tol_calculation_max=1.d-10
  character(len=200), parameter :: subname='check_communication_potential'

  call timing(bigdft_mpi%iproc,'check_pot','ON')

  !assign constants
  i3s=denspot%dpbox%nscatterarr(bigdft_mpi%iproc,3)+1 !< starting point of the planes in the z direction
  n3p=denspot%dpbox%nscatterarr(bigdft_mpi%iproc,2) !< number of planes for the potential
  n2i=denspot%dpbox%ndims(2) !< size of the global domain in y direction
  n1i=denspot%dpbox%ndims(1) !< size of the global domain in x direction

  !fill the values of the rhov array
  ind=0
  do i3=i3s,i3s+n3p-1
     do i2=1,n2i
        do i1=1,n1i
           ind=ind+1
           denspot%rhov(ind)=real(i1+(i2-1)*n1i+(i3-1)*n1i*n2i,dp)
        end do
     end do
  end do

  !calculate the dimensions and communication of the potential element with mpi_get
  call local_potential_dimensions(iproc,tmb%ham_descr%lzd,tmb%orbs,denspot%xc,denspot%dpbox%ngatherarr(0,1))
  call start_onesided_communication(bigdft_mpi%iproc, bigdft_mpi%nproc, max(denspot%dpbox%ndimpot,1), denspot%rhov, &
       tmb%ham_descr%comgp%nrecvbuf, tmb%ham_descr%comgp%recvbuf, tmb%ham_descr%comgp, tmb%ham_descr%lzd)

  !check the fetching of the potential element, destroy the MPI window, results in pot_work
  call full_local_potential(bigdft_mpi%iproc,bigdft_mpi%nproc,tmb%orbs,tmb%ham_descr%lzd,&
       2,denspot%dpbox,denspot%xc,denspot%rhov,denspot%pot_work,tmb%ham_descr%comgp)

  maxdiff=0.0_dp
  sumdiff=0.0_dp
  numtot=0
  loop_lr: do ilr=1,tmb%ham_descr%Lzd%nlr
     !check if this localisation region is used by one of the orbitals
     dosome=.false.
     do iorb=1,tmb%orbs%norbp
        dosome = (tmb%orbs%inwhichlocreg(iorb+tmb%orbs%isorb) == ilr)
        if (dosome) exit
     end do
     if (.not. dosome) cycle loop_lr

     loop_orbs: do iorb=1,tmb%orbs%norbp
        ilr_orb=tmb%orbs%inwhichlocreg(iorb+tmb%orbs%isorb)
        if (ilr_orb /= ilr) cycle loop_orbs

        ind=tmb%orbs%ispot(iorb)-1
        do i3=1,tmb%ham_descr%Lzd%Llr(ilr)%d%n3i
           do i2=1,tmb%ham_descr%Lzd%Llr(ilr)%d%n2i
              do i1=1,tmb%ham_descr%Lzd%Llr(ilr)%d%n1i
                 ind=ind+1
                 testval=real(i1+tmb%ham_descr%Lzd%Llr(ilr)%nsi1+&
                      (i2+tmb%ham_descr%Lzd%Llr(ilr)%nsi2-1)*n1i+&
                      (i3+tmb%ham_descr%Lzd%Llr(ilr)%nsi3-1)*n1i*n2i,dp)
                 testval=abs(denspot%pot_work(ind)-testval)
                 maxdiff=max(maxdiff,testval)
                 sumdiff=sumdiff+testval
                 numtot=numtot+1
              end do
           end do
        end do

     enddo loop_orbs

  end do loop_lr

  if (numtot>0) sumdiff = sumdiff/numtot

  ! Reduce the results
  if (bigdft_mpi%nproc>1) then
      call mpiallred(sumdiff, 1, mpi_sum, bigdft_mpi%mpi_comm)
      call mpiallred(maxdiff, 1, mpi_max, bigdft_mpi%mpi_comm)
  end if
    
  ! Get mean value for the sum
  sumdiff=sqrt(sumdiff)

  if (bigdft_mpi%iproc==0) call yaml_open_map('Checking operations for potential communication')    
  ! Print the results
  if (bigdft_mpi%iproc==0) then
      call yaml_map('Tolerance for the following test',tol_calculation_mean,fmt='(1es25.18)')
      if (sumdiff>tol_calculation_mean) then
         call yaml_warning('CALCULATION ERROR: total difference of '//trim(yaml_toa(sumdiff,fmt='(1es25.18)')))
      else
         call yaml_map('calculation check, error sum', sumdiff,fmt='(1es25.18)')
      end if
      call yaml_map('Tolerance for the following test',tol_calculation_max,fmt='(1es25.18)')
      if (sumdiff>tol_calculation_max) then
         call yaml_warning('CALCULATION ERROR: max difference of '//trim(yaml_toa(maxdiff,fmt='(1es25.18)')))
         call f_err_throw('The communication of the potential is not correct for this setup, check communication routines',&
                 err_name='BIGDFT_MPI_ERROR')
      else
         call yaml_map('calculation check, error max', maxdiff,fmt='(1es25.18)')
      end if
  end if
  if (bigdft_mpi%iproc==0) call yaml_close_map()

  call f_free_ptr(denspot%pot_work)

  nullify(denspot%pot_work)

  call timing(bigdft_mpi%iproc,'check_pot','OF')

end subroutine check_communication_potential


subroutine check_communication_sumrho(iproc, nproc, orbs, lzd, collcom_sr, denspot, denskern, denskern_, check_sumrho)
  use module_base
  use module_types
  use module_interfaces, except_this_one => check_communication_sumrho
  use yaml_output
  use communications, only: transpose_switch_psir, transpose_communicate_psir, transpose_unswitch_psirt
  use sparsematrix_base, only: sparse_matrix, matrices
  use sparsematrix_init, only: matrixindex_in_compressed
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(local_zone_descriptors),intent(in) :: lzd
  type(orbitals_data),intent(in) :: orbs
  type(comms_linear),intent(inout) :: collcom_sr
  type(DFT_local_fields),intent(in) :: denspot
  type(sparse_matrix),intent(inout) :: denskern
  type(matrices),intent(inout) :: denskern_
  integer,intent(in) :: check_sumrho

  ! Local variables
  integer :: ist, iorb, iiorb, ilr, i, iz, ii, iy, ix, iix, iiy, iiz, iixyz, nxyz, ipt, i0, ierr, jproc
  integer :: i1, i2, i3, is1, is2, is3, ie1, ie2, ie3, ii3s, ii3e, nmax, jj, j, ind, ikernel
  integer :: iorbmin, iorbmax, jorb, iall, istat
  real(kind=8) :: maxdiff, sumdiff, tt, tti, ttj, hxh, hyh, hzh, factor, ref_value
  real(kind=8) :: diff
  real(kind=8),dimension(:),allocatable :: psir, psirwork, psirtwork, rho, rho_check
  integer,dimension(:,:,:),allocatable :: weight
  integer,dimension(:,:,:,:),allocatable :: orbital_id
  integer,dimension(:),allocatable :: istarr
  integer,dimension(:,:),allocatable :: matrixindex_in_compressed_auxilliary
  real(kind=8),parameter :: tol_transpose=1.d-14
  real(kind=8),parameter :: tol_calculation_mean=1.d-12
  real(kind=8),parameter :: tol_calculation_max=1.d-10
  character(len=*), parameter :: subname='check_sumrho'
  logical :: rho_negative

  call timing(iproc,'check_sumrho','ON')

  if (iproc==0) call yaml_open_map('Checking operations for sumrho')

  call f_routine(id='check_communication_sumrho')

  ! Allocate all the main arrays arrays
  psir=f_malloc(collcom_sr%ndimpsi_c,id='psir') !direct array

  psirwork=f_malloc(collcom_sr%ndimpsi_c,id='psirwork') !direct workarray

  psirtwork=f_malloc(collcom_sr%ndimind_c,id='psirtwork') !transposed workarray

  ! Size of global box
  nxyz=lzd%glr%d%n1i*lzd%glr%d%n2i*lzd%glr%d%n3i

  call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
  call mpi_barrier(bigdft_mpi%mpi_comm, ierr)

  ! Fill the direct array with a recognizable pattern
  ist=0
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      ilr=orbs%inWhichLocreg(iiorb)
      !$omp parallel default(none) &
      !$omp shared(orbs, lzd, psir, iorb, iiorb, ilr, ist, nxyz) &
      !$omp private(i, ii, iz, iy, ix, iix, iiy, iiz, iixyz)
      !$omp do
      do i=1,lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i
          ! coordinates within locreg
          ii=i-1
          iz=ii/(lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i)+1
          ii=ii-(iz-1)*lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i
          iy=ii/lzd%llr(ilr)%d%n1i+1
          ix=ii-(iy-1)*lzd%llr(ilr)%d%n1i+1
          ! coordinates within global region
          iix=ix+lzd%llr(ilr)%nsi1
          iiy=iy+lzd%llr(ilr)%nsi2
          iiz=iz+lzd%llr(ilr)%nsi3
          iixyz=(iiz-1)*lzd%glr%d%n1i*lzd%glr%d%n2i+(iiy-1)*lzd%glr%d%n1i+iix
          ! assign unique value
          psir(ist+i)=test_value_sumrho(iiorb,iixyz,nxyz)
      end do
      !$omp end do
      !$omp end parallel
      ist = ist + lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i
  end do
  if(ist/=collcom_sr%ndimpsi_c) then
      write(*,'(a,i0,a)') 'ERROR on process ',iproc,' : ist/=collcom_sr%ndimpsi_c'
      stop
  end if

  call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
  call mpi_barrier(bigdft_mpi%mpi_comm, ierr)

  ! Rearrange data
  call transpose_switch_psir(collcom_sr, psir, psirwork)

  ! direct array not needed anymore
  call f_free(psir)

  ! Communicate the data
  call transpose_communicate_psir(iproc, nproc, collcom_sr, psirwork, psirtwork)

  ! Direct workarray not needed anymore
  call f_free(psirwork)

  ! Rearrange array
  call transpose_unswitch_psirt(collcom_sr, psirtwork, collcom_sr%psit_c)

  ! Transposed workarray not needed anymore
  call f_free(psirtwork)

  call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
  call mpi_barrier(bigdft_mpi%mpi_comm, ierr)

  ! Check the layout of the transposed data
  maxdiff=0.d0
  sumdiff=0.d0

  ! Get the starting point of each MPI task
  istarr=f_malloc((/0.to.nproc-1/),id='istarr')
  istarr=0
  istarr(iproc)=collcom_sr%nptsp_c

  if (nproc > 1) then
     call mpiallred(istarr(0), nproc, mpi_sum, bigdft_mpi%mpi_comm)
  end if

  ist=0
  do jproc=0,iproc-1
      ist=ist+istarr(jproc)
  end do
  call f_free(istarr)
  
  ! Iterate through all the transposed values and check whether they are correct
  !$omp parallel default(none) &
  !$omp shared(collcom_sr, ist, nxyz, maxdiff, sumdiff) &
  !$omp private(ipt, ii, i0, iixyz, i, iiorb, tt, ref_value, diff)
  !$omp do reduction(+:sumdiff) reduction(max:maxdiff)
  do ipt=1,collcom_sr%nptsp_c
      ii=collcom_sr%norb_per_gridpoint_c(ipt)
      i0=collcom_sr%isptsp_c(ipt)
      iixyz=ist+ipt
      do i=1,ii
          iiorb=collcom_sr%indexrecvorbital_c(i0+i)
          tt=collcom_sr%psit_c(i0+i)
          ref_value=test_value_sumrho(iiorb,iixyz,nxyz)
          diff=abs(tt-ref_value)
          if (diff>maxdiff) maxdiff=diff
          sumdiff=sumdiff+diff**2
      end do
  end do
  !$omp end do
  !$omp end parallel


  ! Reduce the results
  if (nproc>1) then
      call mpiallred(sumdiff, 1, mpi_sum, bigdft_mpi%mpi_comm)
      call mpiallred(maxdiff, 1, mpi_max, bigdft_mpi%mpi_comm)
  end if
  call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
  call mpi_barrier(bigdft_mpi%mpi_comm, ierr)

  ! Get mean value for the sum
  sumdiff = sumdiff/(lzd%glr%d%n1i*lzd%glr%d%n2i*lzd%glr%d%n3i)
  sumdiff=sqrt(sumdiff)

  ! Print the results
  if (iproc==0) then
      call yaml_map('Tolerance for the following test',tol_transpose,fmt='(1es25.18)')
      if (sumdiff>tol_transpose) then
         call yaml_warning('TRANSPOSITION ERROR: mean difference of '//trim(yaml_toa(sumdiff,fmt='(1es25.18)')))
      else
         call yaml_map('transposition check, mean error ', sumdiff,fmt='(1es25.18)')
      end if
      if (maxdiff>tol_transpose) then
         call yaml_warning('TRANSPOSITION ERROR: max difference of '//trim(yaml_toa(maxdiff,fmt='(1es25.18)')))
      else
         call yaml_map('transposition check, max error ', maxdiff,fmt='(1es25.18)')
      end if
  end if


  ! Now comes the full check.. Do it depending on the value of check_sumrho
  if (check_sumrho==2) then
  
      ! Now simulate the calculation of the charge density. Take the same reproducable
      ! values as above. In this way the charge density can be calculated without
      ! the communication.
      
      ! First determine how many orbitals one has for each grid point in the current slice
      ii3s=denspot%dpbox%nscatterarr(iproc,3)-denspot%dpbox%nscatterarr(iproc,4)+1
      ii3e=denspot%dpbox%nscatterarr(iproc,3)-denspot%dpbox%nscatterarr(iproc,4)+denspot%dpbox%nscatterarr(iproc,1)
      weight=f_malloc0((/lzd%glr%d%n1i,lzd%glr%d%n2i,ii3e-ii3s+1/),lbounds=(/1,1,ii3s/),id='weight')

      if (denspot%dpbox%nscatterarr(iproc,1)>0) then
          call to_zero(lzd%glr%d%n1i*lzd%glr%d%n2i*denspot%dpbox%nscatterarr(iproc,1), weight(1,1,ii3s))
      end if

      do i3=ii3s,ii3e
          do iorb=1,orbs%norb
              ilr=orbs%inwhichlocreg(iorb)
              is3=1+lzd%Llr(ilr)%nsi3
              ie3=lzd%Llr(ilr)%nsi3+lzd%llr(ilr)%d%n3i
              if (is3>i3 .or. i3>ie3) cycle
              is1=1+lzd%Llr(ilr)%nsi1
              ie1=lzd%Llr(ilr)%nsi1+lzd%llr(ilr)%d%n1i
              is2=1+lzd%Llr(ilr)%nsi2
              ie2=lzd%Llr(ilr)%nsi2+lzd%llr(ilr)%d%n2i
              !$omp parallel default(none) &
              !$omp shared(is2, ie2, is1, ie1, weight, i3) &
              !$omp private(i2, i1) 
              !$omp do
              do i2=is2,ie2
                  do i1=is1,ie1
                      weight(i1,i2,i3) = weight(i1,i2,i3)+1
                  end do
              end do
              !$omp end do
              !$omp end parallel
          end do
      end do
    
      ! The array orbital_id contains the IDs of the orbitals touching a given gridpoint
      nmax=maxval(weight)

      orbital_id=f_malloc((/nmax,lzd%glr%d%n1i,lzd%glr%d%n2i,ii3e-ii3s+1/),lbounds=(/1,1,1,ii3s/),id='orbital_id')

      if (denspot%dpbox%nscatterarr(iproc,1)>0) then
          call to_zero(lzd%glr%d%n1i*lzd%glr%d%n2i*denspot%dpbox%nscatterarr(iproc,1), weight(1,1,ii3s))
      end if
      iorbmin=1000000000
      iorbmax=-1000000000
      do i3=ii3s,ii3e
          do iorb=1,orbs%norb
              ilr=orbs%inwhichlocreg(iorb)
              is3=1+lzd%Llr(ilr)%nsi3
              ie3=lzd%Llr(ilr)%nsi3+lzd%llr(ilr)%d%n3i
              if (is3>i3 .or. i3>ie3) cycle
              is1=1+lzd%Llr(ilr)%nsi1
              ie1=lzd%Llr(ilr)%nsi1+lzd%llr(ilr)%d%n1i
              is2=1+lzd%Llr(ilr)%nsi2
              ie2=lzd%Llr(ilr)%nsi2+lzd%llr(ilr)%d%n2i
              !$omp parallel default(none) &
              !$omp shared(is2, ie2, is1, ie1, weight, orbital_id, i3, iorb, iorbmin, iorbmax) &
              !$omp private(i2, i1, jj)
              !$omp do reduction(min:iorbmin) reduction(max:iorbmax)
              do i2=is2,ie2
                  do i1=is1,ie1
                      jj=weight(i1,i2,i3)+1
                      !weight(i1,i2,i3) = weight(i1,i2,i3)+1
                      !orbital_id(weight(i1,i2,i3),i1,i2,i3) = iorb
                      orbital_id(jj,i1,i2,i3) = iorb
                      if (iorb<iorbmin) iorbmin=iorb
                      if (iorb>iorbmax) iorbmax=iorb
                      weight(i1,i2,i3)=jj
                  end do
              end do
              !$omp end do
              !$omp end parallel
          end do
      end do
    
      ! Make sure that the bounds are okay for all processes
      if (iorbmin>iorbmax) then
          iorbmin=1
          iorbmax=1
      end if
    
    
      ! Now calculate the charge density. Of course this is only possible since the
      ! value of each gridpoint is given by the special pattern and therefore always known.
    
      ! First fill the kernel with some numbers.
      do i=1,denskern%nvctr
          denskern%matrix_compr(i)=sine_taylor(real(denskern%nvctr-i+1,kind=8))
      end do
    
      hxh=.5d0*lzd%hgrids(1)
      hyh=.5d0*lzd%hgrids(2)
      hzh=.5d0*lzd%hgrids(3)
      factor=1.d0/(hxh*hyh*hzh)
    
      ! Use an auxilliary array to store the indices of the kernel in the compressed
      ! format. The usual denskern%matrixindex_in_compressed_fortransposed can not be used 
      ! since we are not in the transposed distribution.
      matrixindex_in_compressed_auxilliary=f_malloc((/iorbmin.to.iorbmax,iorbmin.to.iorbmax /), &
          id='matrixindex_in_compressed_auxilliary')

      !$omp parallel default(none) &
      !$omp shared(iorbmin, iorbmax, matrixindex_in_compressed_auxilliary, denskern) &
      !$omp private(iorb, jorb)
      !$omp do
      do iorb=iorbmin,iorbmax
          do jorb=iorbmin,iorbmax
              matrixindex_in_compressed_auxilliary(jorb,iorb)=matrixindex_in_compressed(denskern, jorb, iorb)
          end do
      end do
      !$omp end do
      !$omp end parallel
    
      ! Now calculate the charge density and store the result in rho_check
      rho_check=f_malloc(max(lzd%glr%d%n1i*lzd%glr%d%n2i*(ii3e-ii3s+1),1),id='rho_check')
      !$omp parallel default (none) &
      !$omp private (i3, i2, i1, iixyz, ind, tt, i,j, ii, tti, ikernel, jj, ttj) &
      !$omp shared (ii3s, ii3e, lzd, weight, orbital_id, denskern, rho_check) &
      !$omp shared (nxyz, factor, matrixindex_in_compressed_auxilliary)
      do i3=ii3s,ii3e
          !$omp do
          do i2=1,lzd%glr%d%n2i
              do i1=1,lzd%glr%d%n1i
                  iixyz=(i3-1)*lzd%glr%d%n1i*lzd%glr%d%n2i+(i2-1)*lzd%glr%d%n1i+i1
                  ind=(i3-ii3s)*lzd%glr%d%n1i*lzd%glr%d%n2i+(i2-1)*lzd%glr%d%n1i+i1
                  tt=1.d-20
                  do i=1,weight(i1,i2,i3) !the number of orbitals touching this grid point
                      ii=orbital_id(i,i1,i2,i3)
                      tti=test_value_sumrho(ii,iixyz,nxyz)
                      ikernel=matrixindex_in_compressed_auxilliary(ii,ii)
                      tt=tt+denskern%matrix_compr(ikernel)*tti*tti
                      do j=i+1,weight(i1,i2,i3)
                          jj=orbital_id(j,i1,i2,i3)
                          ikernel=matrixindex_in_compressed_auxilliary(jj,ii)
                          if (ikernel==0) cycle
                          ttj=test_value_sumrho(jj,iixyz,nxyz)
                          tt=tt+2.d0*denskern%matrix_compr(ikernel)*tti*ttj
                      end do
                  end do
                  tt=tt*factor
                  rho_check(ind)=tt
              end do
          end do
          !$omp end do
      end do
      !$omp end parallel
    
      call f_free(matrixindex_in_compressed_auxilliary)
    
      ! Now calculate the charge density in the transposed way using the standard routine
      rho=f_malloc(max(lzd%glr%d%n1i*lzd%glr%d%n2i*(ii3e-ii3s+1),1),id='rho')
      denskern_%matrix_compr = denskern%matrix_compr
      call sumrho_for_TMBs(iproc, nproc, lzd%hgrids(1), lzd%hgrids(2), lzd%hgrids(3), collcom_sr, denskern, denskern_, &
           lzd%glr%d%n1i*lzd%glr%d%n2i*denspot%dpbox%n3d, rho, rho_negative, .false.)
    
      ! Determine the difference between the two versions
      sumdiff=0.d0
      maxdiff=0.d0
      !$omp parallel default(none) shared(lzd,ii3e,ii3s,rho,rho_check,sumdiff,maxdiff) private(i,tt)
      !$omp do reduction(+:sumdiff) reduction(max:maxdiff) 
      do i=1,lzd%glr%d%n1i*lzd%glr%d%n2i*(ii3e-ii3s+1)
          tt=abs(rho(i)-rho_check(i))
          sumdiff = sumdiff + tt**2
          if (tt>maxdiff) maxdiff=tt
      end do
      !$omp end do
      !$omp end parallel
    
      ! Reduce the results
      if (nproc>1) then
          call mpiallred(sumdiff, 1, mpi_sum, bigdft_mpi%mpi_comm)
          call mpiallred(maxdiff, 1, mpi_max, bigdft_mpi%mpi_comm)
      end if
    
      ! Get mean value for the sum
      sumdiff = sumdiff/(lzd%glr%d%n1i*lzd%glr%d%n2i*lzd%glr%d%n3i)
      sumdiff=sqrt(sumdiff)
    
      ! Print the results
      if (iproc==0) then
          call yaml_map('Tolerance for the following test',tol_calculation_mean,fmt='(1es25.18)')
          if (sumdiff>tol_calculation_mean) then
             call yaml_warning('CALCULATION ERROR: total difference of '//trim(yaml_toa(sumdiff,fmt='(1es25.18)')))
          else
             call yaml_map('calculation check, error sum', sumdiff,fmt='(1es25.18)')
          end if
          call yaml_map('Tolerance for the following test',tol_calculation_max,fmt='(1es25.18)')
          if (sumdiff>tol_calculation_max) then
             call yaml_warning('CALCULATION ERROR: max difference of '//trim(yaml_toa(maxdiff,fmt='(1es25.18)')))
          else
             call yaml_map('calculation check, error max', maxdiff,fmt='(1es25.18)')
          end if
      end if
    
    
      call f_free(weight)
      call f_free(orbital_id)
      call f_free(rho_check)
      call f_free(rho)

  end if

  if (iproc==0) call yaml_close_map()

  call timing(iproc,'check_sumrho','OF')

  call f_release_routine()

  contains

    function test_value_sumrho(i, j, n)
      implicit none

      ! Calling arguments
      integer,intent(in) :: i, j, n
      real(kind=8) :: test_value_sumrho

      ! Local variables
      real(kind=8) :: ri, rj, rn
      real(kind=8),parameter :: fac=1.d-8

      ri=real(i,kind=8)
      rj=real(j,kind=8)
      rn=real(n,kind=8)
      !test_value=fac*real((i-1)*n+j,dp)
      !test_value_sumrho=fac*(ri-1.d0)*rn+rj
      test_value_sumrho=sine_taylor((ri-1.d0)*rn)*cosine_taylor(rj)
      !test_value_sumrho=0.d0

    end function test_value_sumrho

    function sine_taylor(xx)
      implicit none

      ! Calling arguments
      real(kind=8),intent(in) :: xx
      real(kind=8) :: sine_taylor

      ! Local variables
      real(kind=8) :: xxtmp, x, x2, x3, x5, x7, x9, x11, x13, x15
      real(kind=8),parameter :: pi=3.14159265358979323846d0
      real(kind=8),parameter :: pi2=6.28318530717958647693d0
      real(kind=8),parameter :: inv6=1.66666666666666666667d-1
      real(kind=8),parameter :: inv120=8.33333333333333333333d-3
      real(kind=8),parameter :: inv5040=1.98412698412698412698d-4
      real(kind=8),parameter :: inv362880=2.75573192239858906526d-6
      real(kind=8),parameter :: inv39916800=2.50521083854417187751d-8
      real(kind=8),parameter :: inv6227020800=1.60590438368216145994d-10
      real(kind=8),parameter :: inv1307674368000=7.6471637318198164759d-13

      ! The Taylor approximation is most accurate around 0, so shift by pi to be centered around this point.
      ! This first part is equivalent to x=mod(xx,pi2)-pi
      x=xx/pi2
      x=real(int(x,kind=8),kind=8)*pi2
      x=xx-x-pi

      x2=x*x
      x3=x2*x
      x5=x3*x2
      x7=x5*x2
      x9=x7*x2
      x11=x9*x2
      x13=x11*x2
      x15=x13*x2

      ! Calculate the value
      sine_taylor = x - x3*inv6 + x5*inv120 - x7*inv5040 + x9*inv362880 &
                    - x11*inv39916800 + x13*inv6227020800 - x15*inv1307674368000

      ! Undo the shift of pi, which corresponds to a multiplication with -1
      sine_taylor=-1.d0*sine_taylor

    end function sine_taylor

    function cosine_taylor(xx)
      implicit none

      ! Calling arguments
      real(kind=8),intent(in) :: xx
      real(kind=8) :: cosine_taylor

      ! Local variables
      real(kind=8) :: x, x2, x4, x6, x8, x10, x12, x14
      real(kind=8),parameter :: pi=3.14159265358979323846d0
      real(kind=8),parameter :: pi2=6.28318530717958647693d0
      real(kind=8),parameter :: inv2=5.d-1
      real(kind=8),parameter :: inv24=4.16666666666666666667d-2
      real(kind=8),parameter :: inv720=1.38888888888888888889d-3
      real(kind=8),parameter :: inv40320=2.48015873015873015873d-5
      real(kind=8),parameter :: inv3628800=2.75573192239858906526d-7
      real(kind=8),parameter :: inv479001600=2.08767569878680989792d-9
      real(kind=8),parameter :: inv87178291200=1.14707455977297247139d-11

      ! The Taylor approximation is most accurate around 0, so shift by pi to be centered around this point.
      ! This first part is equivalent to x=mod(xx,pi2)-pi
      x=xx/pi2
      x=real(int(x,kind=8),kind=8)*pi2
      x=xx-x-pi

      x2=x*x
      x4=x2*x2
      x6=x4*x2
      x8=x6*x2
      x10=x8*x2
      x12=x10*x2
      x14=x12*x2

      ! Calculate the value
      cosine_taylor = 1.d0 - x2*inv2 + x4*inv24 - x6*inv720 + x8*inv40320 &
                      - x10*inv3628800 + x12*inv479001600 - x14*inv87178291200

      ! Undo the shift of pi, which corresponds to a multiplication with -1
      cosine_taylor=-1.d0*cosine_taylor

    end function cosine_taylor

end subroutine check_communication_sumrho



subroutine check_negative_rho(ndimrho, rho, rho_negative)
  use module_base
  implicit none

  ! Calling arguments
  integer,intent(in) :: ndimrho
  real(kind=8),dimension(ndimrho),intent(in) :: rho
  logical,intent(out) :: rho_negative

  ! Local variables
  integer :: i, irho, ierr

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

end subroutine check_negative_rho
