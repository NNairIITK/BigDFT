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
  use Poisson_Solver
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

     Lzd%Llr(ilr)%hybrid_on=.false.

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


subroutine sumrhoForLocalizedBasis2(iproc,nproc,lzd,input,hx,hy,hz,orbs,&
     comsr,densKern,nrho,rho,at,nscatterarr)
!
use module_base
use module_types
use libxc_functionals
use module_interfaces, exceptThisOne => sumrhoForLocalizedBasis2
implicit none

! Calling arguments
integer,intent(in) :: iproc, nproc, nrho
real(gp),intent(in) :: hx, hy, hz
type(local_zone_descriptors),intent(in) :: lzd
type(input_variables),intent(in) :: input
type(orbitals_data),intent(in) :: orbs
!type(p2pCommsSumrho),intent(inout) :: comsr
type(p2pComms),intent(inout) :: comsr
!real(kind=8),dimension(orbs%norb,norb),intent(in) :: coeff
!real(kind=8),dimension(ld_coeff,norb),intent(in) :: coeff
real(kind=8),dimension(orbs%norb,orbs%norb),intent(in) :: densKern
real(8),dimension(nrho),intent(out),target:: rho
type(atoms_data),intent(in) :: at
integer, dimension(0:nproc-1,4),intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh

! Local variables
integer :: iorb, jorb, indLarge, i1, i2, i3, ilr, jlr
integer :: i1s, i1e, i2s, i2e, i3s, i3e, i1d, j1d, i2d, j2d, i3d, j3d, indri, indrj, istri, istrj
integer :: indi2, indi3, indj2, indj3, indl2, indl3, iiorb, jjorb
integer :: ierr, is, ie
integer :: m, i1d0, j1d0, indri0, indrj0, indLarge0
integer :: azones,bzones,ii,izones,jzones,x,y,z,ishift1,ishift2,ishift3,jshift1,jshift2,jshift3
integer,dimension(3,4) :: astart, aend, bstart,bend
real(kind=8) :: tt, hxh, hyh, hzh, factor, totalCharge, tt0, tt1, tt2, tt3, factorTimesDensKern
!real(kind=8),dimension(:,:),allocatable :: densKern
!character(len=*),parameter :: subname='sumrhoForLocalizedBasis2'

if(iproc==0) write(*,'(a)',advance='no') 'Calculating charge density... '


!call mpi_barrier(mpi_comm_world, ierr)
!call cpu_time(t1)
! Calculate the density kernel.
!!if(iproc==0) write(*,'(3x,a)',advance='no') 'calculating the density kernel... '
!!call timing(iproc,'sumrho_TMB    ','ON')
!!if(norbp>0) then
!!    call dgemm('n', 't', orbs%norb, orbs%norb, norbp, 1.d0, coeff(1,isorb+1), ld_coeff, &
!!         coeff(1,isorb+1), ld_coeff, 0.d0, densKern(1,1), orbs%norb)
!!else
!!    call to_zero(orbs%norb**2, densKern(1,1))
!!end if
!!call mpiallred(densKern(1,1), orbs%norb**2, mpi_sum, mpi_comm_world, ierr)
!!call timing(iproc,'sumrho_TMB    ','OF')


!call mpi_barrier(mpi_comm_world, ierr)
!call cpu_time(t2)
!time=t2-t1
!if(iproc==0) write(*,'(a,es12.4)') 'time for kernel:',time

!DEBUG test idempotency of density matrix
!!allocate(rhoprime(orbs%norb,orbs%norb))
!!call dgemm('n','t', orbs%norb,orbs%norb,orbs%norb,1.d0,densKern(1,1),orbs%norb,&
!!     densKern(1,1),orbs%norb,0.d0,rhoprime(1,1),orbs%norb)
!!do i = 1, orbs%norb
!! do j = 1, orbs%norb
!!    if(abs(rhoprime(i,j) - densKern(i,j))>1.0d-5) then
!!      write(*,*) 'Not indempotent',i,j,rhoprime(i,j), densKern(i,j)
!!    end if
!! end do
!!end do
!!deallocate(rhoprime)
!END DEBUG


! Define some constant factors.
hxh=.5d0*hx
hyh=.5d0*hy
hzh=.5d0*hz
if(input%nspin==1) then
    factor=2.d0/(hxh*hyh*hzh)
else
    factor=1.d0/(hxh*hyh*hzh)
end if

! Initialize rho.
if (libxc_functionals_isgga()) then
    call razero(nrho, rho)
else
    ! There is no mpi_allreduce, therefore directly initialize to
    ! 10^-20 and not 10^-20/nproc.
    rho=1.d-20
    !call tenminustwenty(nrho, rho, nproc)
end if
call timing(iproc,'p2pSumrho_wait','ON')


call wait_p2p_communication(iproc, nproc, comsr)



! Now calculate the charge density. Each process calculates only one slice of the total charge density.
! Such a slice has the full extent in the x and y direction, but is limited in the z direction.
! The bounds of the slice are given by nscatterarr. To do so, each process has received all orbitals that
! extend into this slice. The number of these orbitals is given by lin%comsr%noverlaps(iproc).
!call mpi_barrier(mpi_comm_world, ierr)
!call cpu_time(t1)

call timing(iproc,'p2pSumrho_wait','OF')


!!call calculate_charge_density(iproc, nproc, lzd, hxh, hyh, hzh, orbs, densKern, factor, &
!!           nscatterarr, maxval(comsr%noverlaps), comsr%noverlaps, comsr%overlaps, comsr%comarr, comsr%ise3, comsr%nrecvBuf, comsr%recvBuf, nrho, rho)

call timing(iproc,'sumrho_TMB    ','ON')

! Bounds of the slice in global coordinates.
is=nscatterarr(iproc,3) 
ie=is+nscatterarr(iproc,1)-1

! This sum is "symmetric", so only do the second loop (jorb) only up to iorb and multiply by two if iorb/=jorb.
totalCharge=0.d0
!print*,'comsr%noverlaps(iproc)',iproc,comsr%noverlaps(iproc)
do iorb=1,comsr%noverlaps(iproc)
    iiorb=comsr%overlaps(iorb) !global index of orbital iorb
    ilr=orbs%inwhichlocreg(iiorb) !localization region of orbital iorb
    istri=comsr%comarr(5,iorb,iproc)-1 !starting index of orbital iorb in the receive buffer
    do jorb=1,iorb
        jjorb=comsr%overlaps(jorb) !global indes of orbital jorb
        jlr=orbs%inwhichlocreg(jjorb) !localization region of orbital jorb
        istrj=comsr%comarr(5,jorb,iproc)-1 !starting index of orbital jorb in the receive buffer
        !!tt = (lzd%llr(ilr)%locregCenter(1)-lzd%llr(jlr)%locregCenter(1))**2 &
        !!    +(lzd%llr(ilr)%locregCenter(2)-lzd%llr(jlr)%locregCenter(2))**2 &
        !!    +(lzd%llr(ilr)%locregCenter(3)-lzd%llr(jlr)%locregCenter(3))**2
        !!tt=sqrt(tt)
        !!write(200,*) tt, densKern(iiorb,jjorb)
        !!if(tt>6.d0) cycle

        azones = 1
        bzones = 1
        !Calculate the number of regions to cut alr and blr
        do ii=1,2
           if(lzd%llr(ilr)%outofzone(ii) > 0) azones = azones * 2
           if(lzd%llr(jlr)%outofzone(ii) > 0) bzones = bzones * 2
        end do

        !FRACTURE THE FIRST LOCALIZATION REGION
        call fracture_periodic_zone_ISF(azones,lzd%Glr,lzd%Llr(ilr),lzd%Llr(ilr)%outofzone(:),astart,aend)
       
        !FRACTURE SECOND LOCREG
        call fracture_periodic_zone_ISF(bzones,lzd%Glr,lzd%Llr(jlr),lzd%Llr(jlr)%outofzone(:),bstart,bend)
        do izones=1,azones
           do jzones=1,bzones
              ! Bounds of the overlap of orbital iorb and jorb in global coordinates.
              i1s=max(astart(1,izones),bstart(1,jzones))
              i1e=min(aend(1,izones)-1,bend(1,jzones)-1)
              i2s=max(astart(2,izones),bstart(2,jzones))
              i2e=min(aend(2,izones)-1,bend(2,jzones)-1)
              i3s=max(comsr%ise3(iorb,1),comsr%ise3(jorb,1))
              i3e=min(comsr%ise3(iorb,2),comsr%ise3(jorb,2))
              call transform_ISFcoordinates(1,i1s,i2s,i3s,lzd%Glr,lzd%Llr(ilr),x,y,z,ishift1, ishift2, ishift3)
              call transform_ISFcoordinates(1,i1s,i2s,i3s,lzd%Glr,lzd%Llr(jlr),x,y,z,jshift1, jshift2, jshift3)
              factorTimesDensKern = factor*densKern(iiorb,jjorb)
              if(iorb/=jorb) then
                  ! Multiply by two since these elements appear twice (but are calculated only once).
                  factorTimesDensKern = 2.d0*factorTimesDensKern
              end if
              ! Now loop over all points in the box in which the orbitals overlap.
              !if(i3s>i3e) write(*,*) 'no calculation done'

              do i3=i3s,i3e !bounds in z direction
                  i3d=i3 -max(is,-ishift3) !z coordinate of orbital iorb with respect to the overlap box
                  j3d=i3 -max(is,-jshift3) !z coordinate of orbital jorb with respect to the overlap box
                  indi3=i3d*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n1i !z-part of the index of orbital iorb in the 1-dim receive buffer
                  indj3=j3d*lzd%llr(jlr)%d%n2i*lzd%llr(jlr)%d%n1i !z-part of the index of orbital jorb in the 1-dim receive buffer
                  !indl3=(i3-is)*lzd%Glr%d%n2i*lzd%Glr%d%n1i !z-part of the index for which the charge density is beeing calculated
                  !indl3=(modulo(i3-1,Lzd%Glr%d%n3i)-is+1)*lzd%Glr%d%n2i*lzd%Glr%d%n1i !z-part of the index for which the charge density is beeing calculated
                  indl3=(modulo(i3,Lzd%Glr%d%n3i+1)-is)*lzd%Glr%d%n2i*lzd%Glr%d%n1i !z-part of the index for which the charge density is beeing calculated
                  do i2=i2s,i2e !bounds in y direction
                      i2d=i2 + ishift2 !y coordinate of orbital iorb with respect to the overlap box
                      j2d=i2 + jshift2 !y coordinate of orbital jorb with respect to the overlap box
                      indi2=i2d*lzd%llr(ilr)%d%n1i !y-part of the index of orbital iorb in the 1-dim receive buffer
                      indj2=j2d*lzd%llr(jlr)%d%n1i !y-part of the index of orbital jorb in the 1-dim receive buffer
                      !indl2=i2*lzd%Glr%d%n1i !y-part of the index for which the charge density is beeing calculated
                      !indl2=(modulo(i2-1,Lzd%Glr%d%n2i)+1)*lzd%Glr%d%n1i !y-part of the index for which the charge density is beeing calculated
                      indl2=(modulo(i2,Lzd%Glr%d%n2i+1))*lzd%Glr%d%n1i !y-part of the index for which the charge density is beeing calculated
                      ! For all other than free BC, choose m such that the unrolled part is never used.
                      if(Lzd%Glr%geocode=='F') then
                          m=mod(i1e-i1s+1,4)
                      else
                          m=i1e-i1s+1
                      end if
                      if(m/=0) then
                          ! The following five variables hold some intermediate results to speed up the code.
                          i1d0= ishift1 
                          j1d0= jshift1
                          indri0 = indi3 + indi2 + istri + 1
                          indrj0 = indj3 + indj2 + istrj + 1
                          indLarge0 = indl3 + indl2 + 1   
                          do i1=i1s,i1s+m-1
                              i1d=i1d0+i1 !x coordinate of orbital iorb with respect to the overlap box
                              j1d=j1d0+i1 !x coordinate of orbital jorb with respect to the overlap box
                              indri = indri0 + i1d !index of orbital iorb in the 1-dim receive buffer
                              indrj = indrj0 + j1d !index of orbital jorb in the 1-dim receive buffer
                              !indLarge = indLarge0 + i1 !index for which the charge density is beeing calculated
                              !indLarge = indLarge0 + modulo(i1-1,Lzd%Glr%d%n1i)+1 !index for which the charge density is beeing calculated
                              indLarge = indLarge0 + modulo(i1,Lzd%Glr%d%n1i+1) !index for which the charge density is beeing calculated
                              tt = factorTimesDensKern*comsr%recvBuf(indri)*comsr%recvBuf(indrj)
                              rho(indLarge) = rho(indLarge) + tt !update the charge density at point indLarge
                              totalCharge = totalCharge + tt !add the contribution to the total charge
                          end do
                      end if
                      ! This is the same again, this time with unrolled loops.
                      if(i1e-i1s+1>4) then
                          i1d0= ishift1 
                          j1d0= jshift1
                          indri0 = indi3 + indi2 + istri + 1
                          indrj0 = indj3 + indj2 + istrj + 1
                          indLarge0 = indl3 + indl2 + 1
                          do i1=i1s+m,i1e,4
                              i1d=i1d0+i1
                              j1d=j1d0+i1
                              indri = indri0 + i1d
                              indrj = indrj0 + j1d
                              !indLarge = indLarge0 + i1
                              tt0 = factorTimesDensKern*comsr%recvBuf(indri  )*comsr%recvBuf(indrj  )
                              tt1 = factorTimesDensKern*comsr%recvBuf(indri+1)*comsr%recvBuf(indrj+1)
                              tt2 = factorTimesDensKern*comsr%recvBuf(indri+2)*comsr%recvBuf(indrj+2)
                              tt3 = factorTimesDensKern*comsr%recvBuf(indri+3)*comsr%recvBuf(indrj+3)
                              !indLarge = indLarge0 + modulo(i1-1,Lzd%Glr%d%n1i)+1
                              indLarge = indLarge0 + i1
                              rho(indLarge  ) = rho(indLarge  ) + tt0
                              !indLarge = indLarge0 +modulo(i1,Lzd%Glr%d%n1i)+1
                              indLarge = indLarge0 + i1+1
                              rho(indLarge) = rho(indLarge) + tt1
                              !rho(indLarge+1) = rho(indLarge+1) + tt1
                              !indLarge = indLarge0 + modulo(i1+1,Lzd%Glr%d%n1i)+1
                              indLarge = indLarge0 + i1+2
                              rho(indLarge) = rho(indLarge) + tt2
                              !rho(indLarge+2) = rho(indLarge+2) + tt1
                              !indLarge = indLarge0 + modulo(i1+2,Lzd%Glr%d%n1i)+1
                              indLarge = indLarge0 + i1+3
                              rho(indLarge) = rho(indLarge) + tt3
                              !rho(indLarge+3) = rho(indLarge+3) + tt1
                              totalCharge = totalCharge + tt0 + tt1 + tt2 + tt3
                          end do
                      end if
                  end do
              end do
          end do !jzones
       end do !izones
    end do
end do

if(iproc==0) write(*,'(a)') 'done.'

call timing(iproc,'sumrho_TMB    ','OF')

call timing(iproc,'sumrho_allred','ON')

call mpiallred(totalCharge, 1, mpi_sum, mpi_comm_world, ierr)
if(iproc==0) write(*,'(3x,a,es20.12)') 'Calculation finished. TOTAL CHARGE = ', totalCharge*hxh*hyh*hzh

call timing(iproc,'sumrho_allred','OF')



end subroutine sumrhoForLocalizedBasis2






subroutine calculate_density_kernel(iproc, nproc, ld_coeff, orbs, orbs_tmb, coeff, kernel,  ovrlp)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc, ld_coeff
  type(orbitals_data),intent(in):: orbs, orbs_tmb
  real(8),dimension(ld_coeff,orbs%norb),intent(in):: coeff
  real(8),dimension(orbs_tmb%norb,orbs_tmb%norb),intent(out):: kernel
  real(8),dimension(orbs_tmb%norb,orbs_tmb%norb),optional,intent(in):: ovrlp

  ! Local variables
  integer:: istat, iall, ierr
  real(8),dimension(:,:),allocatable:: density_kernel_partial
  character(len=*),parameter:: subname='calculate_density_kernel'
  call timing(iproc,'calc_kernel','ON') !lr408t

  if(iproc==0) write(*,'(1x,a)',advance='no') 'calculate density kernel... '
  allocate(density_kernel_partial(orbs_tmb%norb,max(orbs_tmb%norbp,1)), stat=istat)
  call memocc(istat, density_kernel_partial, 'density_kernel_partial', subname)
  if(orbs_tmb%norbp>0) then
      call dgemm('n', 't', orbs_tmb%norb, orbs_tmb%norbp, orbs%norb, 1.d0, coeff(1,1), ld_coeff, &
           coeff(orbs_tmb%isorb+1,1), ld_coeff, 0.d0, density_kernel_partial(1,1), orbs_tmb%norb)
  end if
  call timing(iproc,'calc_kernel','OF') !lr408t

  if (nproc > 1) then
     call timing(iproc,'commun_kernel','ON') !lr408t
     call mpi_allgatherv(density_kernel_partial(1,1), orbs_tmb%norb*orbs_tmb%norbp, mpi_double_precision, &
          kernel(1,1), orbs_tmb%norb*orbs_tmb%norb_par(:,0), orbs_tmb%norb*orbs_tmb%isorb_par, mpi_double_precision, &
          mpi_comm_world, ierr)
     call timing(iproc,'commun_kernel','OF') !lr408t
  else
     call vcopy(orbs_tmb%norb*orbs_tmb%norbp,density_kernel_partial(1,1),1,kernel(1,1),1)
  end if

  iall=-product(shape(density_kernel_partial))*kind(density_kernel_partial)
  deallocate(density_kernel_partial,stat=istat)
  call memocc(istat,iall,'density_kernel_partial',subname)


end subroutine calculate_density_kernel

subroutine calculate_energy_kernel(iproc, nproc, ld_coeff, orbs, orbs_tmb, coeff, kernel,  ovrlp)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, ld_coeff
  type(orbitals_data),intent(in) :: orbs, orbs_tmb
  real(kind=8),dimension(ld_coeff,orbs%norb),intent(in) :: coeff
  real(kind=8),dimension(orbs_tmb%norb,orbs_tmb%norb),intent(out) :: kernel
  real(kind=8),dimension(orbs_tmb%norb,orbs_tmb%norb),optional,intent(in) :: ovrlp

  ! Local variables
  integer:: istat, iall, ierr, iorb
  real(kind=8),dimension(:,:),allocatable :: density_kernel_partial, coeff_en
  character(len=*),parameter :: subname='calculate_density_kernel'

  call timing(iproc,'calc_kernel','ON') !lr408t

  if(iproc==0) write(*,'(3x,a)') 'calculating the density kernel... '

  allocate(density_kernel_partial(orbs_tmb%norb,max(orbs_tmb%norbp,1)), stat=istat)
  call memocc(istat, density_kernel_partial, 'density_kernel_partial', subname)

  if(orbs_tmb%norbp>0) then
      allocate(coeff_en(ld_coeff,orbs%norb), stat=istat)
      call memocc(istat, coeff_en, 'coeff_en', subname)
      call dcopy(ld_coeff*orbs%norb, coeff,1,coeff_en,1)
      do iorb=1,orbs%norb
         call dscal(ld_coeff,orbs_tmb%eval(iorb),coeff_en(1,iorb),1)
      end do
      call dgemm('n', 't', orbs_tmb%norb, orbs_tmb%norbp, orbs%norb, 1.d0, coeff_en(1,1), ld_coeff, &
           coeff(orbs_tmb%isorb+1,1), ld_coeff, 0.d0, density_kernel_partial(1,1), orbs_tmb%norb)
      iall=-product(shape(coeff_en))*kind(coeff_en)
      deallocate(coeff_en,stat=istat)
      call memocc(istat,iall,'coeff_en',subname)
      
  end if

  call timing(iproc,'calc_kernel','OF') !lr408t

  call timing(iproc,'commun_kernel','ON') !lr408t
  call mpi_allgatherv(density_kernel_partial(1,1), orbs_tmb%norb*orbs_tmb%norbp, mpi_double_precision, &
       kernel(1,1), orbs_tmb%norb*orbs_tmb%norb_par(:,0), orbs_tmb%norb*orbs_tmb%isorb_par, mpi_double_precision, &
       mpi_comm_world, ierr)
  call timing(iproc,'commun_kernel','OF') !lr408t

  iall=-product(shape(density_kernel_partial))*kind(density_kernel_partial)
  deallocate(density_kernel_partial,stat=istat)
  call memocc(istat,iall,'density_kernel_partial',subname)

end subroutine calculate_energy_kernel
