!> @file 
!!   sumrho: linear version
!! @author
!!   Copyright (C) 2011-2012 BigDFT group 
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 

!> Here starts the routine for building partial density inside the localisation region
!! This routine should be treated as a building-block for the linear scaling code
subroutine local_partial_densityLinear(nproc,rsflag,nscatterarr,&
     nrhotot,Lzd,hxh,hyh,hzh,nspin,orbs,mapping,psi,rho)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => local_partial_densityLinear
  use module_xc
  use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
  implicit none
  logical, intent(in) :: rsflag
  integer, intent(in) :: nproc
  integer,intent(inout) :: nrhotot
  integer, intent(in) :: nspin
  real(gp), intent(in) :: hxh,hyh,hzh
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
  allocate(Lnscatterarr(0:nproc-1,4+ndebug),stat=i_stat)
  call memocc(i_stat,Lnscatterarr,'Lnscatterarr',subname)
  Lnscatterarr(:,3) = 0
  Lnscatterarr(:,4) = 0

  !initialize the rho array at 10^-20 instead of zero, due to the invcb ABINIT routine
  !otherwise use libXC routine
  call xc_init_rho(max(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*nrhotot,1)*max(nspin,orbs%nspinor),rho,nproc)

  ind=1
  orbitalsLoop: do ii=1,orbs%norbp

     iorb = ii + orbs%isorb
     ilr = orbs%inwhichLocreg(iorb)

     Lnscatterarr(:,1) = Lzd%Llr(ilr)%d%n3i 
     Lnscatterarr(:,2) = Lzd%Llr(ilr)%d%n3i 


     call initialize_work_arrays_sumrho(Lzd%Llr(ilr),w)
     allocate(rho_p(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*nspinn), stat=i_stat) !must redefine the size of rho_p?
     call memocc(i_stat,rho_p,'rho_p',subname)
     allocate(psir(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,npsir+ndebug),stat=i_stat)
     call memocc(i_stat,psir,'psir',subname)
  
     if (Lzd%Llr(ilr)%geocode == 'F') then
        call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*npsir,psir)
     end if
 
     !Need to zero rho_p
     call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*nspinn, rho_p)

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

     i_all=-product(shape(rho_p))*kind(rho_p)
     deallocate(rho_p,stat=i_stat)
     call memocc(i_stat,i_all,'rho_p',subname)
     i_all=-product(shape(psir))*kind(psir)
     deallocate(psir,stat=i_stat)
     call memocc(i_stat,i_all,'psir',subname)

     call deallocate_work_arrays_sumrho(w)
  end do orbitalsLoop
 
  i_all=-product(shape(Lnscatterarr))*kind(Lnscatterarr)
  deallocate(Lnscatterarr,stat=i_stat)
  call memocc(i_stat,i_all,'Lnscatterarr',subname)
 

END SUBROUTINE local_partial_densityLinear
!
!!!
!!!
subroutine partial_density_linear(rsflag,nproc,n1i,n2i,n3i,npsir,nspinn,nrhotot,&
     hfac,nscatterarr,spinsgn,psir,rho_p,ibyyzz_r) 
  use module_base
  use module_types
  implicit none
  logical, intent(in) :: rsflag
  integer, intent(in) :: nproc,n1i,n2i,n3i,nrhotot,nspinn,npsir
  real(gp), intent(in) :: hfac,spinsgn
  integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr
  real(wp), dimension(n1i,n2i,n3i,npsir), intent(in) :: psir
  real(dp), dimension(n1i,n2i,nrhotot,nspinn), intent(inout) :: rho_p
  integer, dimension(:,:,:),pointer :: ibyyzz_r 
  !local variables
  integer :: i3s,jproc,i3off,n3d,isjmp,i1,i2,i3,i1s,i1e,j3,i3sg
  real(gp) :: hfac2
  real(dp) :: psisq,p1,p2,p3,p4,r1,r2,r3,r4
!  integer :: ncount0,ncount1,ncount_rate,ncount_max
!!!  integer :: ithread,nthread,omp_get_thread_num,omp_get_num_threads
  !sum different slices by taking into account the overlap
  i3sg=0
!$omp parallel default(private) shared(n1i,nproc,rsflag,nspinn,nscatterarr,spinsgn) &
!$omp shared(n2i,npsir,hfac,psir,rho_p,n3i,i3sg,ibyyzz_r)
  i3s=0
!!!   ithread=omp_get_thread_num()
!!!   nthread=omp_get_num_threads()
  hfac2=2.0_gp*hfac

!  call system_clock(ncount0,ncount_rate,ncount_max)

  !case without bounds
  i1s=1
  i1e=n1i
  loop_xc_overlap: do jproc=0,nproc-1
     !case for REDUCE_SCATTER approach, not used for GGA since it enlarges the 
     !communication buffer
     if (rsflag) then
        i3off=nscatterarr(jproc,3)-nscatterarr(jproc,4)
        n3d=nscatterarr(jproc,1)
        if (n3d==0) exit loop_xc_overlap
     else
        i3off=0
        n3d=n3i
     end if
     !here the condition for the MPI_ALLREDUCE should be entered
     if(spinsgn > 0.0_gp) then
        isjmp=1
     else
        isjmp=2
     end if
     do i3=i3off+1,i3off+n3d
        !this allows the presence of GGA with non-isolated BC. If i3 is between 1 and n3i
        !j3=i3. This is useful only when dealing with rsflags and GGA, so we can comment it out
        !j3=modulo(i3-1,n3i)+1 
        j3=i3
        i3s=i3s+1
!!!    if(mod(i3s,nthread) .eq. ithread) then
     !$omp do
        do i2=1,n2i
              i1s=ibyyzz_r(1,i2-15,j3-15)+1
              i1e=ibyyzz_r(2,i2-15,j3-15)+1
           if (npsir == 1) then
              do i1=i1s,i1e
                 !conversion between the different types
                 psisq=real(psir(i1,i2,j3,1),dp)
                 psisq=psisq*psisq
                 rho_p(i1,i2,i3s,isjmp)=rho_p(i1,i2,i3s,isjmp)+real(hfac,dp)*psisq
              end do
           else !similar loop for npsir=4
              do i1=i1s,i1e
                 !conversion between the different types
                 p1=real(psir(i1,i2,j3,1),dp)
                 p2=real(psir(i1,i2,j3,2),dp)
                 p3=real(psir(i1,i2,j3,3),dp)
                 p4=real(psir(i1,i2,j3,4),dp)

                 !density values
                 r1=p1*p1+p2*p2+p3*p3+p4*p4
                 r2=p1*p3+p2*p4
                 r3=p1*p4-p2*p3
                 r4=p1*p1+p2*p2-p3*p3-p4*p4

                 rho_p(i1,i2,i3s,1)=rho_p(i1,i2,i3s,1)+real(hfac,dp)*r1
                 rho_p(i1,i2,i3s,2)=rho_p(i1,i2,i3s,2)+real(hfac2,dp)*r2
                 rho_p(i1,i2,i3s,3)=rho_p(i1,i2,i3s,3)+real(hfac2,dp)*r3
                 rho_p(i1,i2,i3s,4)=rho_p(i1,i2,i3s,4)+real(hfac,dp)*r4
              end do
           end if
        end do
     !$omp enddo
!!!    end if

!$omp critical
        i3sg=max(i3sg,i3s)
!$omp end critical

     end do
     if (.not. rsflag) exit loop_xc_overlap !the whole range is already done
  end do loop_xc_overlap
!$omp end parallel

  if (i3sg /= nrhotot) then
     write(*,'(1x,a,i0,1x,i0)')'ERROR: problem with rho_p: i3s,nrhotot,',i3sg,nrhotot
     stop
  end if

!  call system_clock(ncount1,ncount_rate,ncount_max)
!  write(*,*) 'TIMING:PDF',real(ncount1-ncount0)/real(ncount_rate)
END SUBROUTINE partial_density_linear



subroutine calculate_density_kernel(iproc, nproc, isKernel, orbs, orbs_tmb, coeff, kernel)
  use module_base
  use module_types
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
      call timing(iproc,'calc_kernel','ON') !lr408t
      if(iproc==0) write(*,'(1x,a)',advance='no') 'calculate density kernel... '
      allocate(density_kernel_partial(orbs_tmb%norb,max(orbs_tmb%norbp,1)), stat=istat)
      call memocc(istat, density_kernel_partial, 'density_kernel_partial', subname)
      allocate(fcoeff(orbs_tmb%norb,orbs%norb), stat=istat)
      call memocc(istat, fcoeff, 'fcoeff', subname)
      call to_zero(orbs_tmb%norb*orbs%norb,fcoeff(1,1))
      if(orbs_tmb%norbp>0) then
          !decide wether we calculate the density kernel or just transformation matrix
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
      iall = -product(shape(fcoeff))*kind(fcoeff)
      deallocate(fcoeff,stat=istat)
      call memocc(istat, iall, 'fcoeff', subname)
      call timing(iproc,'calc_kernel','OF') !lr408t

      call timing(iproc,'waitAllgatKern','ON')
      call mpi_barrier(bigdft_mpi%mpi_comm,ierr)
      call timing(iproc,'waitAllgatKern','OF')

      if (nproc > 1) then
         call timing(iproc,'commun_kernel','ON') !lr408t
         allocate(recvcounts(0:nproc-1),stat=istat)
         call memocc(istat,recvcounts,'recvcounts',subname)
         allocate(dspls(0:nproc-1),stat=istat)
         call memocc(istat,recvcounts,'recvcounts',subname)
         do jproc=0,nproc-1
             recvcounts(jproc)=orbs_tmb%norb*orbs_tmb%norb_par(jproc,0)
             dspls(jproc)=orbs_tmb%norb*orbs_tmb%isorb_par(jproc)
         end do
         sendcount=orbs_tmb%norb*orbs_tmb%norbp
         call mpi_allgatherv(density_kernel_partial(1,1), sendcount, mpi_double_precision, &
              kernel(1,1), recvcounts, dspls, mpi_double_precision, &
              bigdft_mpi%mpi_comm, ierr)
         iall=-product(shape(recvcounts))*kind(recvcounts)
         deallocate(recvcounts,stat=istat)
         call memocc(istat,iall,'recvcounts',subname)
         iall=-product(shape(dspls))*kind(dspls)
         deallocate(dspls,stat=istat)
         call memocc(istat,iall,'dspls',subname)
         call timing(iproc,'commun_kernel','OF') !lr408t
      else
         call vcopy(orbs_tmb%norb*orbs_tmb%norbp,density_kernel_partial(1,1),1,kernel(1,1),1)
      end if

      iall=-product(shape(density_kernel_partial))*kind(density_kernel_partial)
      deallocate(density_kernel_partial,stat=istat)
      call memocc(istat,iall,'density_kernel_partial',subname)
  else if (communication_strategy==ALLREDUCE) then
      call timing(iproc,'calc_kernel','ON') !lr408t
      if(iproc==0) write(*,'(1x,a)',advance='no') 'calculate density kernel... '
      if(orbs%norbp>0) then
          allocate(fcoeff(orbs_tmb%norb,orbs%norb), stat=istat)
          call memocc(istat, fcoeff, 'fcoeff', subname)
          call to_zero(orbs_tmb%norb*orbs%norb,fcoeff(1,1))

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
          iall = -product(shape(fcoeff))*kind(fcoeff)
          deallocate(fcoeff,stat=istat)
          call memocc(istat, iall, 'fcoeff', subname)
      else
          call to_zero(orbs_tmb%norb**2, kernel(1,1))
      end if
      call timing(iproc,'calc_kernel','OF') !lr408t

      call timing(iproc,'waitAllgatKern','ON')
      call mpi_barrier(bigdft_mpi%mpi_comm,ierr)
      call timing(iproc,'waitAllgatKern','OF')
      if (nproc > 1) then
          call timing(iproc,'commun_kernel','ON') !lr408t
          call mpiallred(kernel(1,1),orbs_tmb%norb**2, mpi_sum, bigdft_mpi%mpi_comm, ierr)
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

end subroutine calculate_density_kernel

subroutine init_collective_comms_sumro(iproc, nproc, lzd, orbs, nscatterarr, collcom_sr)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(local_zone_descriptors),intent(in) :: lzd
  type(orbitals_data),intent(in) :: orbs
  integer,dimension(0:nproc-1,4),intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  type(collective_comms),intent(inout) :: collcom_sr

  ! Local variables
  integer :: ierr, istat, iall, ipt
  real(kind=8) :: weight_tot, weight_ideal
  integer,dimension(:,:),allocatable :: istartend
  character(len=*),parameter :: subname='determine_weights_sumrho'
  real(kind=8),dimension(:),allocatable :: weights_per_slice, weights_per_zpoint

  ! Note: all weights are double precision to avoid integer overflow
  call timing(iproc,'init_collco_sr','ON')

  allocate(istartend(2,0:nproc-1), stat=istat)
  call memocc(istat, istartend, 'istartend', subname)

  allocate(weights_per_slice(0:nproc-1), stat=istat)
  call memocc(istat, weights_per_slice, 'weights_per_slice', subname)

  allocate(weights_per_zpoint(lzd%glr%d%n3i), stat=istat)
  call memocc(istat, weights_per_zpoint, 'weights_per_zpoint', subname)

  call get_weights_sumrho(iproc, nproc, orbs, lzd, nscatterarr, weight_tot, weight_ideal, &
       weights_per_slice, weights_per_zpoint)
  call mpi_barrier(mpi_comm_world, ierr)

  call assign_weight_to_process_sumrho(iproc, nproc, weight_tot, weight_ideal, weights_per_slice, &
       lzd, orbs, nscatterarr, istartend, collcom_sr%nptsp_c)
  call mpi_barrier(mpi_comm_world, ierr)

  iall = -product(shape(weights_per_slice))*kind(weights_per_slice)
  deallocate(weights_per_slice,stat=istat)
  call memocc(istat, iall, 'weights_per_slice', subname)

  allocate(collcom_sr%norb_per_gridpoint_c(collcom_sr%nptsp_c), stat=istat)
  call memocc(istat, collcom_sr%norb_per_gridpoint_c, 'collcom_sr%norb_per_gridpoint_c', subname)

  call determine_num_orbs_per_gridpoint_sumrho(iproc, nproc, collcom_sr%nptsp_c, lzd, orbs, &
       istartend, weight_tot, weights_per_zpoint, collcom_sr%norb_per_gridpoint_c)
  call mpi_barrier(mpi_comm_world, ierr)

  allocate(collcom_sr%nsendcounts_c(0:nproc-1), stat=istat)
  call memocc(istat, collcom_sr%nsendcounts_c, 'collcom_sr%nsendcounts_c', subname)
  allocate(collcom_sr%nsenddspls_c(0:nproc-1), stat=istat)
  call memocc(istat, collcom_sr%nsenddspls_c, 'collcom_sr%nsenddspls_c', subname)
  allocate(collcom_sr%nrecvcounts_c(0:nproc-1), stat=istat)
  call memocc(istat, collcom_sr%nrecvcounts_c, 'collcom_sr%nrecvcounts_c', subname)
  allocate(collcom_sr%nrecvdspls_c(0:nproc-1), stat=istat)
  call memocc(istat, collcom_sr%nrecvdspls_c, 'collcom_sr%nrecvdspls_c', subname)

  call determine_communication_arrays_sumrho(iproc, nproc, collcom_sr%nptsp_c, lzd, orbs, istartend, &
       collcom_sr%norb_per_gridpoint_c, collcom_sr%nsendcounts_c, collcom_sr%nsenddspls_c, &
       collcom_sr%nrecvcounts_c, collcom_sr%nrecvdspls_c, collcom_sr%ndimpsi_c, collcom_sr%ndimind_c)
  call mpi_barrier(mpi_comm_world, ierr)

  allocate(collcom_sr%psit_c(collcom_sr%ndimind_c), stat=istat)
  call memocc(istat, collcom_sr%psit_c, 'collcom_sr%psit_c', subname)

  allocate(collcom_sr%isendbuf_c(collcom_sr%ndimpsi_c), stat=istat)
  call memocc(istat, collcom_sr%isendbuf_c, 'collcom_sr%isendbuf_c', subname)
  allocate(collcom_sr%irecvbuf_c(collcom_sr%ndimpsi_c), stat=istat)
  call memocc(istat, collcom_sr%irecvbuf_c, 'collcom_sr%irecvbuf_c', subname)
  allocate(collcom_sr%indexrecvorbital_c(collcom_sr%ndimind_c), stat=istat)
  call memocc(istat, collcom_sr%indexrecvorbital_c, 'collcom_sr%indexrecvorbital_c', subname)
  allocate(collcom_sr%iextract_c(collcom_sr%ndimind_c), stat=istat)
  call memocc(istat, collcom_sr%iextract_c, 'collcom_sr%iextract_c', subname)
  allocate(collcom_sr%iexpand_c(collcom_sr%ndimind_c), stat=istat)
  call memocc(istat, collcom_sr%iexpand_c, 'collcom_sr%iexpand_c', subname)

  call get_switch_indices_sumrho(iproc, nproc, collcom_sr%nptsp_c, collcom_sr%ndimpsi_c, collcom_sr%ndimind_c, lzd, &
       orbs, istartend, collcom_sr%norb_per_gridpoint_c, collcom_sr%nsendcounts_c, collcom_sr%nsenddspls_c, &
       collcom_sr%nrecvcounts_c, collcom_sr%nrecvdspls_c, collcom_sr%isendbuf_c, collcom_sr%irecvbuf_c, &
       collcom_sr%iextract_c, collcom_sr%iexpand_c, collcom_sr%indexrecvorbital_c)
  call mpi_barrier(bigdft_mpi%mpi_comm, ierr)

  ! These variables are used in various subroutines to speed up the code
  allocate(collcom_sr%isptsp_c(max(collcom_sr%nptsp_c,1)), stat=istat)
  call memocc(istat, collcom_sr%isptsp_c, 'collcom_sr%isptsp_c', subname)
  collcom_sr%isptsp_c(1) = 0
  do ipt=2,collcom_sr%nptsp_c
        collcom_sr%isptsp_c(ipt) = collcom_sr%isptsp_c(ipt-1) + collcom_sr%norb_per_gridpoint_c(ipt-1)
  end do

  allocate(collcom_sr%nsendcounts_repartitionrho(0:nproc-1), stat=istat)
  call memocc(istat, collcom_sr%nsendcounts_repartitionrho, 'collcom_sr%nsendcounts_repartitionrho', subname)
  allocate(collcom_sr%nrecvcounts_repartitionrho(0:nproc-1), stat=istat)
  call memocc(istat, collcom_sr%nrecvcounts_repartitionrho, 'collcom_sr%nrecvcounts_repartitionrho', subname)
  allocate(collcom_sr%nsenddspls_repartitionrho(0:nproc-1), stat=istat)
  call memocc(istat, collcom_sr%nsenddspls_repartitionrho, 'collcom_sr%nsenddspls_repartitionrho', subname)
  allocate(collcom_sr%nrecvdspls_repartitionrho(0:nproc-1), stat=istat)
  call memocc(istat, collcom_sr%nrecvdspls_repartitionrho, 'collcom_sr%nrecvdspls_repartitionrho', subname)

  call communication_arrays_repartitionrho(iproc, nproc, lzd, nscatterarr, istartend, &
       collcom_sr%nsendcounts_repartitionrho, collcom_sr%nsenddspls_repartitionrho, &
       collcom_sr%nrecvcounts_repartitionrho, collcom_sr%nrecvdspls_repartitionrho)

  iall = -product(shape(weights_per_zpoint))*kind(weights_per_zpoint)
  deallocate(weights_per_zpoint,stat=istat)
  call memocc(istat, iall, 'weights_per_zpoint', subname)

  iall = -product(shape(istartend))*kind(istartend)
  deallocate(istartend,stat=istat)
  call memocc(istat, iall, 'istartend', subname)

  call timing(iproc,'init_collco_sr','OF')

end subroutine init_collective_comms_sumro

subroutine get_weights_sumrho(iproc, nproc, orbs, lzd, nscatterarr, &
           weight_tot, weight_ideal, weights_per_slice, weights_per_zpoint)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(orbitals_data),intent(in) :: orbs
  type(local_zone_descriptors),intent(in) :: lzd
  integer,dimension(0:nproc-1,4),intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  real(kind=8),intent(out) :: weight_tot, weight_ideal
  real(kind=8),dimension(0:nproc-1),intent(out) :: weights_per_slice
  real(kind=8),dimension(lzd%glr%d%n3i),intent(out) :: weights_per_zpoint

  ! Local variables
  integer :: iorb, ilr, ierr, i3, i2, i1, is1, ie1, is2, ie2, is3, ie3, istat, iall
  real(kind=8) :: tt
  real(kind=8),dimension(:,:),allocatable :: weight_xy
  character(len=*),parameter :: subname='get_weights_sumrho'

  !!! Determine the total weight.
  !!weight_tot=0.d0
  !!do iorb=1,orbs%norbp
  !!    iiorb=orbs%isorb+iorb
  !!    ilr=orbs%inwhichlocreg(iiorb)
  !!    ncount = lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i
  !!    weight_tot = weight_tot + dble(ncount)
  !!end do
  !!call mpiallred(weight_tot, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)

  !!! Ideal weight per process
  !!weight_ideal = weight_tot/dble(nproc)


  !!weight_tot=0.d0
  !!!!$omp parallel default(shared) &
  !!!!$omp private(i2, i1, iorb, ilr, is1, ie1, is2, ie2, is3, ie3)
  !!do i3=nscatterarr(iproc,3)+1,nscatterarr(iproc,3)+nscatterarr(iproc,1)
  !!    !!$omp do reduction(+:tt)
  !!    do i2=1,lzd%glr%d%n2i
  !!        do i1=1,lzd%glr%d%n1i
  !!            tt=0.d0
  !!            do iorb=1,orbs%norb
  !!                ilr=orbs%inwhichlocreg(iorb)
  !!                is1=1+lzd%Llr(ilr)%nsi1
  !!                ie1=lzd%Llr(ilr)%nsi1+lzd%llr(ilr)%d%n1i
  !!                is2=1+lzd%Llr(ilr)%nsi2
  !!                ie2=lzd%Llr(ilr)%nsi2+lzd%llr(ilr)%d%n2i
  !!                is3=1+lzd%Llr(ilr)%nsi3
  !!                ie3=lzd%Llr(ilr)%nsi3+lzd%llr(ilr)%d%n3i
  !!                if (is1<=i1 .and. i1<=ie1 .and. is2<=i2 .and. i2<=ie2 .and. is3<=i3 .and. i3<=ie3) then
  !!                    tt=tt+1.d0
  !!                end if
  !!            end do
  !!            weight_tot=weight_tot+tt**2
  !!        end do
  !!    end do
  !!    !!$omp end do
  !!end do
  !!!!$omp end parallel

  call to_zero(lzd%glr%d%n3i, weights_per_zpoint(1))

  allocate(weight_xy(lzd%glr%d%n1i,lzd%glr%d%n2i), stat=istat)
  call memocc(istat, weight_xy, 'weight_xy', subname)

  tt=0.d0
  weights_per_slice(:) = 0.0d0
  do i3=nscatterarr(iproc,3)+1,nscatterarr(iproc,3)+nscatterarr(iproc,1)
      !tmp=0.d0
      call to_zero(lzd%glr%d%n1i*lzd%glr%d%n2i, weight_xy(1,1))
      !!$omp parallel default(shared) &
      !!$omp private(i2, i1, iorb, ilr, is1, ie1, is2, ie2, is3, ie3, ttt)
      !!$omp do reduction(+:tmp)
      do iorb=1,orbs%norb
          ilr=orbs%inwhichlocreg(iorb)
          is3=1+lzd%Llr(ilr)%nsi3
          ie3=lzd%Llr(ilr)%nsi3+lzd%llr(ilr)%d%n3i
          if (is3>i3 .or. i3>ie3) cycle
          is1=1+lzd%Llr(ilr)%nsi1
          ie1=lzd%Llr(ilr)%nsi1+lzd%llr(ilr)%d%n1i
          is2=1+lzd%Llr(ilr)%nsi2
          ie2=lzd%Llr(ilr)%nsi2+lzd%llr(ilr)%d%n2i
          do i2=is2,ie2
              do i1=is1,ie1
                  weight_xy(i1,i2) = weight_xy(i1,i2)+1.d0
              end do
          end do
      end do
      !!$omp end do
      !!$omp end parallel
      weights_per_zpoint(i3)=0.d0
      do i2=1,lzd%glr%d%n2i
          do i1=1,lzd%glr%d%n1i
             tt = tt + .5d0*(weight_xy(i1,i2)*(weight_xy(i1,i2)+1.d0))
             weights_per_zpoint(i3) = weights_per_zpoint(i3) + .5d0*(weight_xy(i1,i2)*(weight_xy(i1,i2)))
          end do
      end do
  end do
  weights_per_slice(iproc)=tt
  call mpiallred(weights_per_slice(0), nproc, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  if (nproc>1) then
     call mpi_allreduce(tt, weight_tot, 1, mpi_double_precision, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  else
     weight_tot=tt
  end if
  call mpiallred(weights_per_zpoint(1), lzd%glr%d%n3i, mpi_sum, bigdft_mpi%mpi_comm, ierr)


  iall = -product(shape(weight_xy))*kind(weight_xy)
  deallocate(weight_xy,stat=istat)
  call memocc(istat, iall, 'weight_xy', subname)
  !!call mpiallred(weight_tot, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)

  ! Ideal weight per process
  weight_ideal = weight_tot/dble(nproc)

end subroutine get_weights_sumrho

subroutine assign_weight_to_process_sumrho(iproc, nproc, weight_tot, weight_ideal, weights_per_slice, &
           lzd, orbs, nscatterarr, istartend, nptsp)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  real(kind=8),intent(in) :: weight_tot, weight_ideal
  real(kind=8),dimension(0:nproc-1),intent(in) :: weights_per_slice
  type(local_zone_descriptors),intent(in) :: lzd
  type(orbitals_data),intent(in) :: orbs
  integer,dimension(0:nproc-1,4),intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  integer,dimension(2,0:nproc-1),intent(out) :: istartend
  integer,intent(out) :: nptsp

  ! Local variables
  integer :: jproc, i1, i2, i3, ii, iorb, ilr, is1, ie1, is2, ie2, is3, ie3, ierr, istat, iall, jproc_out
  real(kind=8) :: tt, ttt
  real(kind=8),dimension(:,:),allocatable :: slicearr
  real(8),dimension(:,:),allocatable :: weights_startend
  character(len=*),parameter :: subname='assign_weight_to_process_sumrho'

  !!allocate(weights_per_slice(0:nproc-1), stat=istat)
  !!call memocc(istat, weights_per_slice, 'weights_per_slice', subname)

  allocate(weights_startend(2,0:nproc-1), stat=istat)
  call memocc(istat, weights_startend, 'weights_startend', subname)

  tt=0.d0
  weights_startend(1,0)=0.d0
  do jproc=0,nproc-2
      tt=tt+weight_ideal
      weights_startend(2,jproc)=dble(floor(tt,kind=8))
      weights_startend(1,jproc+1)=dble(floor(tt,kind=8))+1.d0
  end do
  weights_startend(2,nproc-1)=weight_tot

  !!call to_zero(nproc, weights_per_slice(0))
  ! Iterate through all grid points and assign them to processes such that the
  ! load balancing is optimal.

  !!do jproc=0,nproc-1
  !!    if (iproc==0) write(*,'(a,i7,2f16.1)') 'jproc, start, end', iproc, weights_startend(1,jproc), weights_startend(2,jproc)
  !!end do

  !!if (nproc>1) then
  !!    tt=0.d0
  !!    jproc=0
  !!    istartend(1,jproc)=1
  !!    !$omp parallel default(shared) &
  !!    !$omp private(i2, i1, iorb, ilr, is1, ie1, is2, ie2, is3, ie3)
  !!    do i3=nscatterarr(iproc,3)+1,nscatterarr(iproc,3)+nscatterarr(iproc,1)
  !!        !$omp do reduction(+:tt)
  !!        do i2=1,lzd%glr%d%n2i
  !!            do i1=1,lzd%glr%d%n1i
  !!                ttt=0.d0
  !!                do iorb=1,orbs%norb
  !!                    ilr=orbs%inwhichlocreg(iorb)
  !!                    is1=1+lzd%Llr(ilr)%nsi1
  !!                    ie1=lzd%Llr(ilr)%nsi1+lzd%llr(ilr)%d%n1i
  !!                    is2=1+lzd%Llr(ilr)%nsi2
  !!                    ie2=lzd%Llr(ilr)%nsi2+lzd%llr(ilr)%d%n2i
  !!                    is3=1+lzd%Llr(ilr)%nsi3
  !!                    ie3=lzd%Llr(ilr)%nsi3+lzd%llr(ilr)%d%n3i
  !!                    if (is1<=i1 .and. i1<=ie1 .and. is2<=i2 .and. i2<=ie2 .and. is3<=i3 .and. i3<=ie3) then
  !!                        !tt=tt+1.d0
  !!                        ttt=ttt+1.d0
  !!                    end if
  !!                end do
  !!                tt=tt+ttt**2
  !!            end do
  !!        end do
  !!        !$omp end do
  !!    end do
  !!    !$omp end parallel
  !!    weights_per_slice(iproc)=tt
  !!    call mpiallred(weights_per_slice(0), nproc, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  !!end if


  ! Iterate through all grid points and assign them to processes such that the
  ! load balancing is optimal.
  if (nproc==1) then
      istartend(1,0)=1
      istartend(2,0)=lzd%glr%d%n1i*lzd%glr%d%n2i*lzd%glr%d%n3i
  else
      allocate(slicearr(lzd%glr%d%n1i,lzd%glr%d%n2i), stat=istat)
      call memocc(istat, slicearr, 'slicearr', subname)
      istartend(1,:)=0
      istartend(2,:)=0
      tt=0.d0
      jproc=0
      ii=0
      outer_loop: do jproc_out=0,nproc-1
          if (tt+weights_per_slice(jproc_out)<weights_startend(1,iproc)) then
              tt=tt+weights_per_slice(jproc_out)
              ii=ii+nscatterarr(jproc_out,1)*lzd%glr%d%n1i*lzd%glr%d%n2i
              cycle outer_loop
          end if
          !!i3_loop: do i3=nscatterarr(jproc_out,3)+1,nscatterarr(jproc_out,3)+nscatterarr(jproc_out,1)
          !!    do i2=1,lzd%glr%d%n2i
          !!        do i1=1,lzd%glr%d%n1i
          !!            ii=ii+1
          !!            ttt=0.d0
          !!            !!$omp parallel if (orbs%norb>512) &
          !!            !!$omp default(shared) &
          !!            !!$omp private(iorb, ilr, is1, ie1, is2, ie2, is3, ie3)
          !!            !!$omp do reduction(+:ttt)
          !!            do iorb=1,orbs%norb
          !!                ilr=orbs%inwhichlocreg(iorb)
          !!                is1=1+lzd%Llr(ilr)%nsi1
          !!                ie1=lzd%Llr(ilr)%nsi1+lzd%llr(ilr)%d%n1i
          !!                is2=1+lzd%Llr(ilr)%nsi2
          !!                ie2=lzd%Llr(ilr)%nsi2+lzd%llr(ilr)%d%n2i
          !!                is3=1+lzd%Llr(ilr)%nsi3
          !!                ie3=lzd%Llr(ilr)%nsi3+lzd%llr(ilr)%d%n3i
          !!                if (is1<=i1 .and. i1<=ie1 .and. is2<=i2 .and. i2<=ie2 .and. is3<=i3 .and. i3<=ie3) then
          !!                    ttt=ttt+1.d0
          !!                end if
          !!            end do
          !!            !!$omp end do
          !!            !!$omp end parallel
          !!            !tt=tt+ttt
          !!            !tt=tt+ttt**2
          !!            tt=tt+.5d0*ttt*(ttt+1.d0)
          !!            if (tt>=weights_startend(1,iproc)) then
          !!                istartend(1,iproc)=ii
          !!                exit outer_loop
          !!            end if
          !!        end do
          !!    end do
          !!end do i3_loop
          i3_loop: do i3=nscatterarr(jproc_out,3)+1,nscatterarr(jproc_out,3)+nscatterarr(jproc_out,1)
              call to_zero(lzd%glr%d%n1i*lzd%glr%d%n2i, slicearr(1,1))
              do iorb=1,orbs%norb
                  ilr=orbs%inwhichlocreg(iorb)
                  is1=1+lzd%Llr(ilr)%nsi1
                  ie1=lzd%Llr(ilr)%nsi1+lzd%llr(ilr)%d%n1i
                  is2=1+lzd%Llr(ilr)%nsi2
                  ie2=lzd%Llr(ilr)%nsi2+lzd%llr(ilr)%d%n2i
                  is3=1+lzd%Llr(ilr)%nsi3
                  ie3=lzd%Llr(ilr)%nsi3+lzd%llr(ilr)%d%n3i
                  if (is3>i3 .or. i3>ie3) cycle
                  do i2=1,lzd%glr%d%n2i
                      do i1=1,lzd%glr%d%n1i
                          ttt=0.d0
                          !!$omp parallel if (orbs%norb>512) &
                          !!$omp default(shared) &
                          !!$omp private(iorb, ilr, is1, ie1, is2, ie2, is3, ie3)
                          !!$omp do reduction(+:ttt)
                          if (is1<=i1 .and. i1<=ie1 .and. is2<=i2 .and. i2<=ie2) then
                              slicearr(i1,i2)=slicearr(i1,i2)+1.d0
                          end if
                          !!$omp end do
                          !!$omp end parallel
                          !tt=tt+ttt
                          !tt=tt+ttt**2
                      end do
                  end do
               end do
              do i2=1,lzd%glr%d%n2i
                  do i1=1,lzd%glr%d%n1i
                      ii=ii+1
                      tt=tt+.5d0*slicearr(i1,i2)*(slicearr(i1,i2)+1.d0)
                      if (tt>=weights_startend(1,iproc)) then
                          istartend(1,iproc)=ii
                          exit outer_loop
                      end if
                  end do
               end do
           end do i3_loop
        end do outer_loop
      iall = -product(shape(slicearr))*kind(slicearr)
      deallocate(slicearr,stat=istat)
      call memocc(istat, iall, 'slicearr', subname)
  end if

  call mpiallred(istartend(1,0), 2*nproc, mpi_sum, bigdft_mpi%mpi_comm, ierr)

  do jproc=0,nproc-2
      istartend(2,jproc)=istartend(1,jproc+1)-1
  end do
  istartend(2,nproc-1)=lzd%glr%d%n1i*lzd%glr%d%n2i*lzd%glr%d%n3i

  !weightp=istartend(2,iproc)-istartend(1,iproc)+1
  !!if (iproc==0) then
  !!    do jproc=0,nproc-1
  !!        write(*,'(a,2i9)') 'istartend(1,jproc), istartend(2,jproc)', istartend(1,jproc), istartend(2,jproc)
  !!    end do
  !!end if


  do jproc=0,nproc-1
      if (iproc==jproc) then
          nptsp=istartend(2,jproc)-istartend(1,jproc)+1
      end if
  end do


  !!iall = -product(shape(weights_per_slice))*kind(weights_per_slice)
  !!deallocate(weights_per_slice,stat=istat)
  !!call memocc(istat, iall, 'weights_per_slice', subname)
  iall = -product(shape(weights_startend))*kind(weights_startend)
  deallocate(weights_startend,stat=istat)
  call memocc(istat, iall, 'weights_startend', subname)


  ! Some check
  ii=nptsp
  call mpiallred(ii, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  if (ii/=lzd%glr%d%n1i*lzd%glr%d%n2i*lzd%glr%d%n3i) then
      stop 'ii/=lzd%glr%d%n1i*lzd%glr%d%n2i*lzd%glr%d%n3i'
  end if



end subroutine assign_weight_to_process_sumrho

subroutine determine_num_orbs_per_gridpoint_sumrho(iproc, nproc, nptsp, lzd, orbs, &
           istartend, weight_tot, weights_per_zpoint, norb_per_gridpoint)
  use module_base
  use module_types
  use yaml_output
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, nptsp
  type(local_zone_descriptors),intent(in) :: lzd
  type(orbitals_data),intent(in) :: orbs
  integer,dimension(2,0:nproc-1),intent(in) :: istartend
  real(kind=8),intent(in) :: weight_tot
  real(kind=8),dimension(lzd%glr%d%n3i),intent(in) :: weights_per_zpoint
  integer,dimension(nptsp),intent(out) :: norb_per_gridpoint

  ! Local variables
  integer :: i3, ii, i2, i1, ipt, ilr, is1, ie1, is2, ie2, is3, ie3, iorb, ierr, i
  real(8) :: tt, weight_check
  !logical :: fast


!!t1=mpi_wtime()
!!  weight_check=0.d0
!!  do i3=1,lzd%glr%d%n3i
!!      if (i3*lzd%glr%d%n1i*lzd%glr%d%n2i<istartend(1,iproc) .or. &
!!          (i3-1)*lzd%glr%d%n1i*lzd%glr%d%n2i+1>istartend(2,iproc)) then
!!          cycle
!!      end if
!!      if (weights_per_zpoint(i3)==0.d0) then
!!          fast=.true.
!!      else
!!          fast=.false.
!!      end if
!!      write(*,*) 'iproc, fast', iproc, fast
!!      tt=0.d0
!!      !$omp parallel default(shared) &
!!      !$omp private(i2, i1, ii, ipt, norb, iorb, ilr, is1, ie1, is2, ie2, is3, ie3)
!!      !$omp do reduction(+:tt)
!!      do i2=1,lzd%glr%d%n2i
!!          do i1=1,lzd%glr%d%n1i
!!              ii=(i3-1)*lzd%glr%d%n1i*lzd%glr%d%n2i+(i2-1)*lzd%glr%d%n1i+i1
!!              if (ii>=istartend(1,iproc) .and. ii<=istartend(2,iproc)) then
!!                  ipt=ii-istartend(1,iproc)+1
!!                  norb=0
!!                  if (.not.fast) then
!!                      do iorb=1,orbs%norb
!!                          ilr=orbs%inwhichlocreg(iorb)
!!                          is1=1+lzd%Llr(ilr)%nsi1
!!                          ie1=lzd%Llr(ilr)%nsi1+lzd%llr(ilr)%d%n1i
!!                          is2=1+lzd%Llr(ilr)%nsi2
!!                          ie2=lzd%Llr(ilr)%nsi2+lzd%llr(ilr)%d%n2i
!!                          is3=1+lzd%Llr(ilr)%nsi3
!!                          ie3=lzd%Llr(ilr)%nsi3+lzd%llr(ilr)%d%n3i
!!                          if (is1<=i1 .and. i1<=ie1 .and. is2<=i2 .and. i2<=ie2 .and. is3<=i3 .and. i3<=ie3) then
!!                          norb=norb+1
!!                          end if
!!                      end do
!!                  end if
!!                  norb_per_gridpoint(ipt)=norb
!!                  !tt=tt+dble(norb**2)
!!                  tt=tt+.5d0*dble(norb*(norb+1))
!!              end if
!!          end do
!!      end do
!!      !$omp end do
!!      !$omp end parallel
!!      weight_check=weight_check+tt
!!  end do
!!t2=mpi_wtime()
!!write(*,*) 'iproc, individual time', iproc, t2-t1
!!
!!  ! Some check
!!  call mpiallred(weight_check, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
!!  if (weight_check/=weight_tot) then
!!      stop '2: tt/=weight_tot'
!!  end if






!t1=mpi_wtime()
  call to_zero(nptsp, norb_per_gridpoint(1))
  do i3=1,lzd%glr%d%n3i
      if (i3*lzd%glr%d%n1i*lzd%glr%d%n2i<istartend(1,iproc) .or. &
          (i3-1)*lzd%glr%d%n1i*lzd%glr%d%n2i+1>istartend(2,iproc)) then
          cycle
      end if
      if (weights_per_zpoint(i3)==0.d0) then
          cycle
      end if
      !!$omp parallel default(shared) &
      !!$omp private(i2, i1, ii, ipt, iorb, ilr, is1, ie1, is2, ie2, is3, ie3)
      !!$omp do
      do iorb=1,orbs%norb
          ilr=orbs%inwhichlocreg(iorb)
          is3=1+lzd%Llr(ilr)%nsi3
          ie3=lzd%Llr(ilr)%nsi3+lzd%llr(ilr)%d%n3i
          if (is3>i3 .or. i3>ie3) cycle
          is2=1+lzd%Llr(ilr)%nsi2
          ie2=lzd%Llr(ilr)%nsi2+lzd%llr(ilr)%d%n2i
          is1=1+lzd%Llr(ilr)%nsi1
          ie1=lzd%Llr(ilr)%nsi1+lzd%llr(ilr)%d%n1i
          do i2=is2,ie2
              do i1=is1,ie1
                  ii=(i3-1)*lzd%glr%d%n1i*lzd%glr%d%n2i+(i2-1)*lzd%glr%d%n1i+i1
                  if (ii>=istartend(1,iproc) .and. ii<=istartend(2,iproc)) then
                      ipt=ii-istartend(1,iproc)+1
                      norb_per_gridpoint(ipt)=norb_per_gridpoint(ipt)+1
                  end if
              end do
          end do
      end do
      !!$omp end do
      !!$omp end parallel
  end do

  tt=0.d0
  do i=1,nptsp
      tt=tt+.5d0*dble(norb_per_gridpoint(i)*(norb_per_gridpoint(i)+1))
  end do
  weight_check=tt
!t2=mpi_wtime()
!write(*,*) 'iproc, individual time', iproc, t2-t1


  ! Some check
  call mpiallred(weight_check, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  if (abs(weight_check-weight_tot) > 1.d-3) then
      stop '2: tt/=weight_tot'
  else if (abs(weight_check-weight_tot) > 0.d0) then
     call yaml_warning('The total weight for density seems inconsistent! Ref:'//&
           trim(yaml_toa(weight_tot,fmt='(1pe25.17)'))//', Check:'//&
           trim(yaml_toa(weight_check,fmt='(1pe25.17)')))
  end if

end subroutine determine_num_orbs_per_gridpoint_sumrho

subroutine determine_communication_arrays_sumrho(iproc, nproc, nptsp, lzd, orbs, &
           istartend, norb_per_gridpoint, nsendcounts, nsenddspls, nrecvcounts, &
           nrecvdspls, ndimpsi, ndimind)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, nptsp
  type(local_zone_descriptors),intent(in) :: lzd
  type(orbitals_data),intent(in) :: orbs
  integer,dimension(2,0:nproc-1),intent(in) :: istartend
  integer,dimension(nptsp),intent(in) :: norb_per_gridpoint
  integer,dimension(0:nproc-1),intent(out) :: nsendcounts, nsenddspls, nrecvcounts, nrecvdspls
  integer,intent(out) :: ndimpsi, ndimind

  ! Local variables
  integer :: iorb, iiorb, ilr, is1, ie1, is2, ie2, is3, ie3, jproc, i3, i2, i1, ind, ii, istat, iall, ierr
  integer,dimension(:),allocatable :: nsendcounts_tmp, nsenddspls_tmp, nrecvcounts_tmp, nrecvdspls_tmp
  character(len=*),parameter :: subname='determine_communication_arrays_sumrho'


  nsendcounts=0
  !!$omp parallel default(shared) &
  !!$omp private(jproc, i3, i2, i1, ind)
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      ilr=orbs%inwhichlocreg(iiorb)
      is1=1+lzd%Llr(ilr)%nsi1
      ie1=lzd%Llr(ilr)%nsi1+lzd%llr(ilr)%d%n1i
      is2=1+lzd%Llr(ilr)%nsi2
      ie2=lzd%Llr(ilr)%nsi2+lzd%llr(ilr)%d%n2i
      is3=1+lzd%Llr(ilr)%nsi3
      ie3=lzd%Llr(ilr)%nsi3+lzd%llr(ilr)%d%n3i
      !!$omp do
      do jproc=0,nproc-1
          do i3=is3,ie3
              if (i3*lzd%glr%d%n1i*lzd%glr%d%n2i<istartend(1,jproc) .or. &
                  (i3-1)*lzd%glr%d%n1i*lzd%glr%d%n2i+1>istartend(2,jproc)) then
                  cycle
              end if
              do i2=is2,ie2
                  do i1=is1,ie1
                    ind = (i3-1)*lzd%glr%d%n1i*lzd%glr%d%n2i+(i2-1)*lzd%glr%d%n1i+i1
                    if (ind>=istartend(1,jproc) .and. ind<=istartend(2,jproc)) then
                        nsendcounts(jproc)=nsendcounts(jproc)+1
                    end if
                  end do
              end do
          end do
       end do
       !!$omp end do
  end do
  !!$omp end parallel


  ! Some check
  ii=0
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      ilr=orbs%inwhichlocreg(iiorb)
      ii = ii + lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i
  end do
  if (ii/=sum(nsendcounts)) then
      stop 'ii/=sum(nsendcounts)'
  end if
  ndimpsi=ii


  nsenddspls(0)=0
  do jproc=1,nproc-1
      nsenddspls(jproc)=nsenddspls(jproc-1)+nsendcounts(jproc-1)
  end do

  allocate(nsendcounts_tmp(0:nproc-1), stat=istat)
  call memocc(istat, nsendcounts_tmp, 'nsendcounts_tmp', subname)
  allocate(nsenddspls_tmp(0:nproc-1), stat=istat)
  call memocc(istat, nsenddspls_tmp, 'nsenddspls_tmp', subname)
  allocate(nrecvcounts_tmp(0:nproc-1), stat=istat)
  call memocc(istat, nrecvcounts_tmp, 'nrecvcounts_tmp', subname)
  allocate(nrecvdspls_tmp(0:nproc-1), stat=istat)
  call memocc(istat, nrecvdspls_tmp, 'nrecvdspls_tmp', subname)
  nsendcounts_tmp=1
  nrecvcounts_tmp=1
  do jproc=0,nproc-1
      nsenddspls_tmp(jproc)=jproc
      nrecvdspls_tmp(jproc)=jproc
  end do
  if(nproc>1) then
      call mpi_alltoallv(nsendcounts, nsendcounts_tmp, nsenddspls_tmp, mpi_integer, nrecvcounts, &
           nrecvcounts_tmp, nrecvdspls_tmp, mpi_integer, bigdft_mpi%mpi_comm, ierr)
  else
      nrecvcounts=nsendcounts
  end if
  iall=-product(shape(nsendcounts_tmp))*kind(nsendcounts_tmp)
  deallocate(nsendcounts_tmp, stat=istat)
  call memocc(istat, iall, 'nsendcounts_tmp', subname)
  iall=-product(shape(nsenddspls_tmp))*kind(nsenddspls_tmp)
  deallocate(nsenddspls_tmp, stat=istat)
  call memocc(istat, iall, 'nsenddspls_tmp', subname)
  iall=-product(shape(nrecvcounts_tmp))*kind(nrecvcounts_tmp)
  deallocate(nrecvcounts_tmp, stat=istat)
  call memocc(istat, iall, 'nrecvcounts_tmp', subname)
  iall=-product(shape(nrecvdspls_tmp))*kind(nrecvdspls_tmp)
  deallocate(nrecvdspls_tmp, stat=istat)
  call memocc(istat, iall, 'nrecvdspls_tmp', subname)

  ndimind = sum(nrecvcounts)

  ! Some check
  ii=sum(norb_per_gridpoint)
  if (ii/=ndimind) stop 'ii/=sum(nrecvcounts)'

  nrecvdspls(0)=0
  do jproc=1,nproc-1
      nrecvdspls(jproc)=nrecvdspls(jproc-1)+nrecvcounts(jproc-1)
  end do

end subroutine determine_communication_arrays_sumrho

subroutine get_switch_indices_sumrho(iproc, nproc, nptsp, ndimpsi, ndimind, lzd, orbs, istartend, &
           norb_per_gridpoint, nsendcounts, nsenddspls, nrecvcounts, nrecvdspls, &
           isendbuf, irecvbuf, iextract, iexpand, indexrecvorbital)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, nptsp, ndimpsi, ndimind
  type(local_zone_descriptors),intent(in) :: lzd
  type(orbitals_data),intent(in) :: orbs
  integer,dimension(2,0:nproc-1),intent(in) :: istartend
  integer,dimension(nptsp),intent(in) :: norb_per_gridpoint
  integer,dimension(0:nproc-1),intent(in) :: nsendcounts, nsenddspls, nrecvcounts, nrecvdspls
  integer,dimension(ndimpsi),intent(out) :: isendbuf, irecvbuf
  integer,dimension(ndimind),intent(out) :: iextract, iexpand, indexrecvorbital

  ! Local variables
  integer :: jproc, iitot, iiorb, ilr, is1, ie1, is2, ie2, is3, ie3, i3, i2, i1, ind, indglob, istat, iall, ierr, ii
  integer :: iorb, i, ipt
  integer,dimension(:),allocatable :: nsend, indexsendbuf, indexsendorbital, indexsendorbital2, indexrecvorbital2
  integer,dimension(:),allocatable :: gridpoint_start, indexrecvbuf
  character(len=*),parameter :: subname='get_switch_indices_sumrho'


  allocate(nsend(0:nproc-1), stat=istat)
  call memocc(istat, nsend, 'nsend', subname)
  nsend=0
  allocate(indexsendbuf(ndimpsi), stat=istat)
  call memocc(istat, indexsendbuf, 'indexsendbuf', subname)
  allocate(indexsendorbital(ndimpsi), stat=istat)
  call memocc(istat, indexsendorbital, 'indexsendorbital', subname)
  !!allocate(isendbuf(ndimpsi), stat=istat)
  !!call memocc(istat, isendbuf, 'isendbuf', subname)

  
  !!$omp parallel default(shared) &
  !!$omp private(iorb, iiorb, ilr, is1, ie1, is2, ie2, is3, ie3, i3, i2, i1, indglob, ind)
  !!$omp do lastprivate(iitot)
  do jproc=0,nproc-1
      iitot=0
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          ilr=orbs%inwhichlocreg(iiorb)
          is1=1+lzd%Llr(ilr)%nsi1
          ie1=lzd%Llr(ilr)%nsi1+lzd%llr(ilr)%d%n1i
          is2=1+lzd%Llr(ilr)%nsi2
          ie2=lzd%Llr(ilr)%nsi2+lzd%llr(ilr)%d%n2i
          is3=1+lzd%Llr(ilr)%nsi3
          ie3=lzd%Llr(ilr)%nsi3+lzd%llr(ilr)%d%n3i
          do i3=is3,ie3
              if (i3*lzd%glr%d%n1i*lzd%glr%d%n2i<istartend(1,jproc) .or. &
                  (i3-1)*lzd%glr%d%n1i*lzd%glr%d%n2i+1>istartend(2,jproc)) then
                  iitot=iitot+(ie2-is2+1)*(ie1-is1+1)
                  cycle
              end if
              do i2=is2,ie2
                  do i1=is1,ie1
                      indglob = (i3-1)*lzd%glr%d%n1i*lzd%glr%d%n2i+(i2-1)*lzd%glr%d%n1i+i1
                      iitot=iitot+1
                      if (indglob>=istartend(1,jproc) .and. indglob<=istartend(2,jproc)) then
                          nsend(jproc)=nsend(jproc)+1
                          ind=nsenddspls(jproc)+nsend(jproc)
                          isendbuf(iitot)=ind
                          indexsendbuf(ind)=indglob
                          indexsendorbital(iitot)=iiorb
                          !exit
                      end if
                  end do
              end do
          end do
      end do
  end do
  !!$omp end do
  !!$omp end parallel


  if(iitot/=ndimpsi) stop 'iitot/=ndimpsi'

  !check
  do jproc=0,nproc-1
      if(nsend(jproc)/=nsendcounts(jproc)) stop 'nsend(jproc)/=nsendcounts(jproc)'
  end do

!!call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
!!t2=mpi_wtime()
!!tt=t2-t1
!!if(iproc==0) write(*,*) 'time 5.1: iproc', iproc, tt



  !!allocate(irecvbuf(ndimpsi), stat=istat)
  !!call memocc(istat, irecvbuf, 'irecvbuf', subname)

  allocate(indexsendorbital2(ndimpsi), stat=istat)
  call memocc(istat, indexsendorbital2, 'indexsendorbital2', subname)
  indexsendorbital2=indexsendorbital
  do i=1,ndimpsi
      ind=isendbuf(i)
      indexsendorbital(ind)=indexsendorbital2(i)
  end do

  ! Inverse of isendbuf
  call get_reverse_indices(ndimpsi, isendbuf, irecvbuf)

  iall=-product(shape(indexsendorbital2))*kind(indexsendorbital2)
  deallocate(indexsendorbital2, stat=istat)
  call memocc(istat, iall, 'indexsendorbital2', subname)


  allocate(indexrecvbuf(ndimind), stat=istat)
  call memocc(istat, indexrecvbuf, 'indexrecvbuf', subname)
  !!allocate(indexrecvorbital(ndimind), stat=istat)
  !!call memocc(istat, indexrecvorbital, 'indexrecvorbital', subname)

  if(nproc>1) then
      ! Communicate indexsendbuf
      call mpi_alltoallv(indexsendbuf, nsendcounts, nsenddspls, mpi_integer, indexrecvbuf, &
           nrecvcounts, nrecvdspls, mpi_integer, bigdft_mpi%mpi_comm, ierr)
      ! Communicate indexsendorbitals
      call mpi_alltoallv(indexsendorbital, nsendcounts, nsenddspls, &
           mpi_integer, indexrecvorbital, &
           nrecvcounts, nrecvdspls, mpi_integer, bigdft_mpi%mpi_comm, ierr)
   else
       indexrecvbuf=indexsendbuf
       indexrecvorbital=indexsendorbital
   end if

  iall=-product(shape(indexsendbuf))*kind(indexsendbuf)
  deallocate(indexsendbuf, stat=istat)
  call memocc(istat, iall, 'indexsendbuf', subname)

  iall=-product(shape(indexsendorbital))*kind(indexsendorbital)
  deallocate(indexsendorbital, stat=istat)
  call memocc(istat, iall, 'indexsendorbital', subname)
!!call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
!!t2=mpi_wtime()
!!tt=t2-t1
!!if(iproc==0) write(*,*) 'time 5.2: iproc', iproc, tt


   allocate(gridpoint_start(istartend(1,iproc):istartend(2,iproc)), stat=istat)
   call memocc(istat, gridpoint_start, 'gridpoint_start', subname)

   ii=1
   do ipt=1,nptsp
       i=ipt+istartend(1,iproc)-1
       if (norb_per_gridpoint(ipt)>0) then
           gridpoint_start(i)=ii
       else
           gridpoint_start(i)=0
       end if
       ii=ii+norb_per_gridpoint(ipt)
   end do

   if (ii/=ndimind+1) stop '(ii/=ndimind+1)'
   if(maxval(gridpoint_start)>ndimind) stop '1: maxval(gridpoint_start)>sum(nrecvcountc)'

   !!allocate(iextract(ndimind), stat=istat)
   !!call memocc(istat, iextract, 'iextract', subname)

  ! Rearrange the communicated data
  do i=1,ndimind
      ii=indexrecvbuf(i)
      ind=gridpoint_start(ii)
      iextract(i)=ind
      gridpoint_start(ii)=gridpoint_start(ii)+1
  end do

  if(maxval(iextract)>ndimind) stop 'maxval(iextract)>ndimind'
  if(minval(iextract)<1) stop 'minval(iextract)<1'

  iall=-product(shape(indexrecvbuf))*kind(indexrecvbuf)
  deallocate(indexrecvbuf, stat=istat)
  call memocc(istat, iall, 'indexrecvbuf', subname)


  !! allocate(iexpand(ndimind), stat=istat)
  !! call memocc(istat, iexpand, 'iexpand', subname)
  ! Get the array to transfrom back the data
  call get_reverse_indices(ndimind, iextract, iexpand)

!!call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
!!t2=mpi_wtime()
!!tt=t2-t1
!!if(iproc==0) write(*,*) 'time 5.3: iproc', iproc, tt

  allocate(indexrecvorbital2(ndimind), stat=istat)
  call memocc(istat, indexrecvorbital2, 'indexrecvorbital2', subname)

  call vcopy(ndimind, indexrecvorbital(1), 1, indexrecvorbital2(1), 1)

  !!$omp parallel default(shared) private(i, ind)
  !!$omp do
  do i=1,ndimind
      ind=iextract(i)
      indexrecvorbital(ind)=indexrecvorbital2(i)
  end do
  !!$omp end do
  !!$omp end parallel

  iall=-product(shape(indexrecvorbital2))*kind(indexrecvorbital2)
  deallocate(indexrecvorbital2, stat=istat)
  call memocc(istat, iall, 'indexrecvorbital2', subname)

  if(minval(indexrecvorbital)<1) stop 'minval(indexrecvorbital)<1'
  if(maxval(indexrecvorbital)>orbs%norb) stop 'maxval(indexrecvorbital)>orbs%norb'


  iall=-product(shape(gridpoint_start))*kind(gridpoint_start)
  deallocate(gridpoint_start, stat=istat)
  call memocc(istat, iall, 'gridpoint_start', subname)

  iall=-product(shape(nsend))*kind(nsend)
  deallocate(nsend, stat=istat)
  call memocc(istat, iall, 'nsend', subname)


end subroutine get_switch_indices_sumrho

subroutine communication_arrays_repartitionrho(iproc, nproc, lzd, nscatterarr, istartend, &
           nsendcounts_repartitionrho, nsenddspls_repartitionrho, &
           nrecvcounts_repartitionrho, nrecvdspls_repartitionrho)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(local_zone_descriptors),intent(in) :: lzd
  integer,dimension(0:nproc-1,4),intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  integer,dimension(2,0:nproc-1),intent(in) :: istartend
  integer,dimension(0:nproc-1),intent(out) :: nsendcounts_repartitionrho, nsenddspls_repartitionrho
  integer,dimension(0:nproc-1),intent(out) :: nrecvcounts_repartitionrho, nrecvdspls_repartitionrho

  ! Local variables
  integer :: jproc_send, jproc_recv, ii, i3, i2, i1, jproc

  jproc_send=0
  jproc_recv=0
  ii=0
  nsendcounts_repartitionrho=0
  nrecvcounts_repartitionrho=0
  do i3=1,lzd%glr%d%n3i
      do i2=1,lzd%glr%d%n2i
          do i1=1,lzd%glr%d%n1i
              ii=ii+1
              if (ii>istartend(2,jproc_send)) then
                  jproc_send=jproc_send+1
              end if
              if (i3>nscatterarr(jproc_recv,3)+nscatterarr(jproc_recv,1)) then
                  jproc_recv=jproc_recv+1
              end if
              if (iproc==jproc_send) then
                  nsendcounts_repartitionrho(jproc_recv)=nsendcounts_repartitionrho(jproc_recv)+1
              end if
              if (iproc==jproc_recv) then
                  nrecvcounts_repartitionrho(jproc_send)=nrecvcounts_repartitionrho(jproc_send)+1
              end if
          end do
      end do
  end do

  nsenddspls_repartitionrho(0)=0
  nrecvdspls_repartitionrho(0)=0
  do jproc=1,nproc-1
      nsenddspls_repartitionrho(jproc)=nsenddspls_repartitionrho(jproc-1)+&
                                                  nsendcounts_repartitionrho(jproc-1)
      nrecvdspls_repartitionrho(jproc)=nrecvdspls_repartitionrho(jproc-1)+&
                                                  nrecvcounts_repartitionrho(jproc-1)
  end do


end subroutine communication_arrays_repartitionrho

subroutine transpose_switch_psir(collcom_sr, psir, psirwork)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(collective_comms),intent(in) :: collcom_sr
  real(kind=8),dimension(collcom_sr%ndimpsi_c),intent(in) :: psir
  real(kind=8),dimension(collcom_sr%ndimpsi_c),intent(out) :: psirwork

  ! Local variables
  integer :: i, m, ind

  !$omp parallel default(private) &
  !$omp shared(collcom_sr, psir, psirwork, m)

  m = mod(collcom_sr%ndimpsi_c,7)
  if(m/=0) then
      do i=1,m
          ind = collcom_sr%isendbuf_c(i)
          psirwork(ind) = psir(i)
      end do
  end if
  !$omp do
  do i = m+1,collcom_sr%ndimpsi_c,7
     psirwork(collcom_sr%isendbuf_c(i+0)) = psir(i+0)
     psirwork(collcom_sr%isendbuf_c(i+1)) = psir(i+1)
     psirwork(collcom_sr%isendbuf_c(i+2)) = psir(i+2)
     psirwork(collcom_sr%isendbuf_c(i+3)) = psir(i+3)
     psirwork(collcom_sr%isendbuf_c(i+4)) = psir(i+4)
     psirwork(collcom_sr%isendbuf_c(i+5)) = psir(i+5)
     psirwork(collcom_sr%isendbuf_c(i+6)) = psir(i+6)
  end do
  !$omp end do
  !$omp end parallel


end subroutine transpose_switch_psir

subroutine transpose_communicate_psir(iproc, nproc, collcom_sr, psirwork, psirtwork)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(collective_comms),intent(in) :: collcom_sr
  real(kind=8),dimension(collcom_sr%ndimpsi_c),intent(in) :: psirwork
  real(kind=8),dimension(collcom_sr%ndimind_c),intent(out) :: psirtwork

  ! Local variables
  integer :: ierr


  if (nproc>1) then
      call mpi_alltoallv(psirwork, collcom_sr%nsendcounts_c, collcom_sr%nsenddspls_c, mpi_double_precision, psirtwork, &
           collcom_sr%nrecvcounts_c, collcom_sr%nrecvdspls_c, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
  else
      call vcopy(collcom_sr%ndimpsi_c, psirwork(1), 1, psirtwork(1), 1)
  end if


end subroutine transpose_communicate_psir

subroutine transpose_unswitch_psirt(collcom_sr, psirtwork, psirt)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(collective_comms),intent(in) :: collcom_sr
  real(kind=8),dimension(collcom_sr%ndimind_c),intent(in) :: psirtwork
  real(kind=8),dimension(collcom_sr%ndimind_c),intent(out) :: psirt

  ! Local variables
  integer :: i, ind, sum_c, m

  sum_c = sum(collcom_sr%nrecvcounts_c)

  !$omp parallel private(i,ind) &
  !$omp shared(psirt, psirtwork, collcom_sr, sum_c, m)

  m = mod(sum_c,7)

  if(m/=0) then
    do i = 1,m
      ind=collcom_sr%iextract_c(i)
      psirt(ind)=psirtwork(i)
    end do
  end if

  !$omp do
  do i=m+1, sum_c,7
      psirt(collcom_sr%iextract_c(i+0))=psirtwork(i+0)
      psirt(collcom_sr%iextract_c(i+1))=psirtwork(i+1)
      psirt(collcom_sr%iextract_c(i+2))=psirtwork(i+2)
      psirt(collcom_sr%iextract_c(i+3))=psirtwork(i+3)
      psirt(collcom_sr%iextract_c(i+4))=psirtwork(i+4)
      psirt(collcom_sr%iextract_c(i+5))=psirtwork(i+5)
      psirt(collcom_sr%iextract_c(i+6))=psirtwork(i+6)
  end do
  !$omp end do
  !$omp end parallel

end subroutine transpose_unswitch_psirt

subroutine transpose_switch_psirt(collcom_sr, psirt, psirtwork)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(collective_comms),intent(in) :: collcom_sr
  real(kind=8),dimension(collcom_sr%ndimind_c),intent(in) :: psirt
  real(kind=8),dimension(collcom_sr%ndimind_c),intent(out) :: psirtwork

  ! Local variables
  integer :: i, ind, sum_c, m

  sum_c = sum(collcom_sr%nrecvcounts_c)

  !$omp parallel default(private) &
  !$omp shared(collcom_sr, psirt, psirtwork, sum_c, m)

  m = mod(sum_c,7)

  if(m/=0) then
    do i=1,m
       ind = collcom_sr%iexpand_c(i)
       psirtwork(ind) = psirt(i)
    end do
  end if


  !$omp do
  do i=m+1,sum_c,7
      psirtwork(collcom_sr%iexpand_c(i+0))=psirt(i+0)
      psirtwork(collcom_sr%iexpand_c(i+1))=psirt(i+1)
      psirtwork(collcom_sr%iexpand_c(i+2))=psirt(i+2)
      psirtwork(collcom_sr%iexpand_c(i+3))=psirt(i+3)
      psirtwork(collcom_sr%iexpand_c(i+4))=psirt(i+4)
      psirtwork(collcom_sr%iexpand_c(i+5))=psirt(i+5)
      psirtwork(collcom_sr%iexpand_c(i+6))=psirt(i+6)
  end do
  !$omp end do
  !$omp end parallel

end subroutine transpose_switch_psirt

subroutine transpose_communicate_psirt(iproc, nproc, collcom_sr, psirtwork, psirwork)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(collective_comms),intent(in) :: collcom_sr
  real(kind=8),dimension(collcom_sr%ndimind_c),intent(in) :: psirtwork
  real(kind=8),dimension(collcom_sr%ndimpsi_c),intent(out) :: psirwork

  ! Local variables
  integer :: ierr

  if (nproc>1) then
  call mpi_alltoallv(psirtwork, collcom_sr%nrecvcounts_c, collcom_sr%nrecvdspls_c, mpi_double_precision, psirwork, &
       collcom_sr%nsendcounts_c, collcom_sr%nsenddspls_c, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
  else
      call vcopy(collcom_sr%ndimpsi_c, psirtwork(1), 1, psirwork(1), 1)
  end if

end subroutine transpose_communicate_psirt

subroutine transpose_unswitch_psir(collcom_sr, psirwork, psir)
  use module_base
  use module_types
  implicit none

  ! Caling arguments
  type(collective_comms),intent(in) :: collcom_sr
  real(kind=8),dimension(collcom_sr%ndimpsi_c),intent(in) :: psirwork
  real(kind=8),dimension(collcom_sr%ndimpsi_c),intent(out) :: psir

  ! Local variables
  integer :: i, ind, m


  !$omp parallel default(private) &
  !$omp shared(collcom_sr, psirwork, psir, m)

  m = mod(collcom_sr%ndimpsi_c,7)

  if(m/=0) then
    do i = 1,m
     ind=collcom_sr%irecvbuf_c(i)
     psir(ind)=psirwork(i)
    end do
  end if

  ! coarse part

  !$omp do
    do i=m+1,collcom_sr%ndimpsi_c,7
        psir(collcom_sr%irecvbuf_c(i+0))=psirwork(i+0)
        psir(collcom_sr%irecvbuf_c(i+1))=psirwork(i+1)
        psir(collcom_sr%irecvbuf_c(i+2))=psirwork(i+2)
        psir(collcom_sr%irecvbuf_c(i+3))=psirwork(i+3)
        psir(collcom_sr%irecvbuf_c(i+4))=psirwork(i+4)
        psir(collcom_sr%irecvbuf_c(i+5))=psirwork(i+5)
        psir(collcom_sr%irecvbuf_c(i+6))=psirwork(i+6)
    end do
  !$omp end do
  !$omp end parallel

end subroutine transpose_unswitch_psir


subroutine sumrho_for_TMBs(iproc, nproc, hx, hy, hz, collcom_sr, denskern, ndimrho, rho)
  use module_base
  use module_types
  use libxc_functionals
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, ndimrho
  real(kind=8),intent(in) :: hx, hy, hz
  type(collective_comms),intent(in) :: collcom_sr
  type(sparseMatrix),intent(in) :: denskern
  real(kind=8),dimension(ndimrho),intent(out) :: rho

  ! Local variables
  integer :: ipt, ii, i0, iiorb, jjorb, istat, iall, i, j, ierr, ind
  real(8) :: tt, total_charge, hxh, hyh, hzh, factor, tt1
  real(kind=8),dimension(:),allocatable :: rho_local
  character(len=*),parameter :: subname='sumrho_for_TMBs'


  allocate(rho_local(collcom_sr%nptsp_c), stat=istat)
  call memocc(istat, rho_local, 'rho_local', subname)

  ! Define some constant factors.
  hxh=.5d0*hx
  hyh=.5d0*hy
  hzh=.5d0*hz
  factor=1.d0/(hxh*hyh*hzh)

  call timing(iproc,'sumrho_TMB    ','ON')
  
  ! Initialize rho. (not necessary for the moment)
  !if (libxc_functionals_isgga()) then
  !    call razero(collcom_sr%nptsp_c, rho_local)
  !else
   !   ! There is no mpi_allreduce, therefore directly initialize to
   !   ! 10^-20 and not 10^-20/nproc.
  !    rho_local=1.d-20
  !end if

  if (iproc==0) write(*,'(a)', advance='no') 'Calculating charge density... '

  total_charge=0.d0
  !$omp parallel default(private) &
  !$omp shared(total_charge, collcom_sr, factor, denskern, rho_local)
  !$omp do schedule(static,50) reduction(+:total_charge)
  do ipt=1,collcom_sr%nptsp_c
      ii=collcom_sr%norb_per_gridpoint_c(ipt)
      i0=collcom_sr%isptsp_c(ipt)
      tt=1.e-20_dp
      do i=1,ii
          iiorb=collcom_sr%indexrecvorbital_c(i0+i)
          tt1=collcom_sr%psit_c(i0+i)
          ind=denskern%matrixindex_in_compressed_fortransposed(iiorb,iiorb)
          tt=tt+denskern%matrix_compr(ind)*tt1*tt1
          do j=i+1,ii
              jjorb=collcom_sr%indexrecvorbital_c(i0+j)
              ind=denskern%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
              if (ind==0) cycle
              tt=tt+2.0_dp*denskern%matrix_compr(ind)*tt1*collcom_sr%psit_c(i0+j)
          end do
      end do
      tt=factor*tt
      total_charge=total_charge+tt
      rho_local(ipt)=tt
  end do
  !$omp end do
  !$omp end parallel

  if (iproc==0) write(*,'(a)') 'done.'

  call timing(iproc,'sumrho_TMB    ','OF')

  call timing(iproc,'sumrho_allred','ON')

  ! Communicate the density to meet the shape required by the Poisson solver.
  if (nproc>1) then
      call mpi_alltoallv(rho_local, collcom_sr%nsendcounts_repartitionrho, collcom_sr%nsenddspls_repartitionrho, &
                         mpi_double_precision, rho, collcom_sr%nrecvcounts_repartitionrho, &
                         collcom_sr%nrecvdspls_repartitionrho, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
  else
      call vcopy(ndimrho, rho_local(1), 1, rho(1), 1)
  end if

  call mpiallred(total_charge, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)

  if(iproc==0) write(*,'(3x,a,es20.12)') 'Calculation finished. TOTAL CHARGE = ', total_charge*hxh*hyh*hzh
  
  call timing(iproc,'sumrho_allred','OF')

  iall=-product(shape(rho_local))*kind(rho_local)
  deallocate(rho_local, stat=istat)
  call memocc(istat, iall, 'rho_local', subname)

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




subroutine check_communication_sumrho(iproc, nproc, orbs, lzd, collcom_sr)
  use module_base
  use module_types
  use module_interfaces
  use yaml_output
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(local_zone_descriptors),intent(in) :: lzd
  type(orbitals_data),intent(in) :: orbs
  type(collective_comms),intent(in) :: collcom_sr

  ! Local variables
  integer :: ist, iorb, iiorb, ilr, i, iz, ii, iy, ix, iix, iiy, iiz, iixyz, nxyz, ipt, i0, ierr, jproc
  real(kind=8) :: maxdiff, sumdiff, tt
  real(kind=8),dimension(:),allocatable :: psir, psirwork, psirt, psirtwork
  integer,dimension(:),allocatable :: istarr
  real(kind=8),parameter :: tol=1.d-10

  if (iproc==0) call yaml_open_map('Checking transposition for sumrho')

  call f_routine(id='check_communication_sumrho')

  ! Allocate dummy array
  psir=f_malloc(collcom_sr%ndimpsi_c,id='collcom_sr%ndimpsi_c')

  ! Size of global box
  nxyz=lzd%glr%d%n1i*lzd%glr%d%n2i*lzd%glr%d%n3i

  ! Fill with a recognizable pattern
  ist=0
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      ilr=orbs%inWhichLocreg(iiorb)
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
          psir(ist+i)=real((iiorb-1)*nxyz+iixyz,dp)
      end do
      ist = ist + lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i
  end do
  if(ist/=collcom_sr%ndimpsi_c) then
      write(*,'(a,i0,a)') 'ERROR on process ',iproc,' : ist/=collcom_sr%ndimpsi_c'
      stop
  end if

  ! Direct workarray
  psirwork=f_malloc(collcom_sr%ndimpsi_c,id='collcom_sr%ndimpsi_c')

  ! Rearrange data
  call transpose_switch_psir(collcom_sr, psir, psirwork)

  ! Original array not needed anymore
  call f_free(psir)

  ! Allocate transposed workarray
  psirtwork=f_malloc(collcom_sr%ndimind_c,id='psirtwork')

  ! Communicate the data
  call transpose_communicate_psir(iproc, nproc, collcom_sr, psirwork, psirtwork)

  ! Direct workarray not needed anymore
  call f_free(psirwork)

  ! Allocate transposed array
  psirt=f_malloc(collcom_sr%ndimind_c,id='psirt')
  call transpose_unswitch_psirt(collcom_sr, psirtwork, psirt)

  ! Transposed workarray not needed anymore
  call f_free(psirtwork)


  ! Check the layout of the transposed data
  maxdiff=0.d0
  sumdiff=0.d0

  ! Get the starting point of each MPI task
  !istarr=f_malloc(0:nproc-1,id='istarr')
  allocate(istarr(0:nproc-1))
  istarr=0
  istarr(iproc)=collcom_sr%nptsp_c
  call mpiallred(istarr(0), nproc, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  ist=0
  do jproc=0,iproc-1
      ist=ist+istarr(jproc)
  end do
  deallocate(istarr)
  !call f_free(istarr)
  
  ! Iterate through all the transposed values and check whether they are correct
  do ipt=1,collcom_sr%nptsp_c
      ii=collcom_sr%norb_per_gridpoint_c(ipt)
      i0=collcom_sr%isptsp_c(ipt)
      iixyz=ist+ipt
      do i=1,ii
          iiorb=collcom_sr%indexrecvorbital_c(i0+i)
          tt=psirt(i0+i)
          maxdiff=max(maxdiff,abs(tt-real((iiorb-1)*nxyz+iixyz,dp)))
          sumdiff=sumdiff+abs(tt-real((iiorb-1)*nxyz+iixyz,dp))
      end do
  end do

  call f_free(psirt)

  ! Reduce the results
  if (nproc>1) then
      call mpiallred(maxdiff, 1, mpi_max, bigdft_mpi%mpi_comm, ierr)
      call mpiallred(maxdiff, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  end if
  if (iproc==0) then
      call yaml_map('Tolerances for this check',tol,fmt='(1es23.15)')
      if (sumdiff>tol) then
         call yaml_warning('TRANSPOSITION ERROR: total difference of '//trim(yaml_toa(sumdiff,fmt='(1es23.15)')))
         call yaml_warning('TRANSPOSITION ERROR: max difference of '//trim(yaml_toa(maxdiff,fmt='(1es23.15)')))
      else
         call yaml_map('transposition check, sum', sumdiff,fmt='(1es23.15)')
         call yaml_map('transposition check, max', maxdiff,fmt='(1es23.15)')
      end if
  end if

  if (iproc==0) call yaml_close_map()


end subroutine check_communication_sumrho
