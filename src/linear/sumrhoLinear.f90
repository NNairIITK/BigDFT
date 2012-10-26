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




! Define some constant factors.
hxh=.5d0*hx
hyh=.5d0*hy
hzh=.5d0*hz
!if(input%nspin==1) then
    factor=1.d0/(hxh*hyh*hzh)
!else
!    factor=1.d0/(hxh*hyh*hzh)
!end if

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
!call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
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

call mpiallred(totalCharge, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
if(iproc==0) write(*,'(3x,a,es20.12)') 'Calculation finished. TOTAL CHARGE = ', totalCharge*hxh*hyh*hzh

call timing(iproc,'sumrho_allred','OF')



end subroutine sumrhoForLocalizedBasis2

subroutine calculate_density_kernel(iproc, nproc, isKernel, ld_coeff, orbs, orbs_tmb, coeff, kernel,overlap)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc, ld_coeff
  type(orbitals_data),intent(in):: orbs, orbs_tmb
  logical, intent(in) :: isKernel
  real(8),dimension(ld_coeff,orbs%norb),intent(in):: coeff
  real(8),dimension(orbs_tmb%norb,orbs_tmb%norb),intent(out):: kernel
  real(8),dimension(orbs_tmb%norb,orbs_tmb%norb),intent(in), optional:: overlap

  ! Local variables
  integer:: istat, iall, ierr, sendcount, jproc, iorb, itmb
  real(8),dimension(:,:),allocatable:: density_kernel_partial, fcoeff,ks,ksk,ksksk
  character(len=*),parameter:: subname='calculate_density_kernel'
  integer,dimension(:),allocatable:: recvcounts, dspls
  integer,parameter:: ALLGATHERV=1, ALLREDUCE=2
  integer,parameter:: communication_strategy=ALLREDUCE

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

          call dgemm('n', 't', orbs_tmb%norb, orbs_tmb%norbp, orbs%norb, 1.d0, coeff(1,1), ld_coeff, &
               fcoeff(orbs_tmb%isorb+1,1), ld_coeff, 0.d0, density_kernel_partial(1,1), orbs_tmb%norb)
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
  end if


  if (communication_strategy==ALLREDUCE) then
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
          call dgemm('n', 't', orbs_tmb%norb, orbs_tmb%norb, orbs%norbp, 1.d0, coeff(1,orbs%isorb+1), ld_coeff, &
               fcoeff(1,orbs%isorb+1), ld_coeff, 0.d0, kernel(1,1), orbs_tmb%norb)
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
 !call dgemm('n', 't', orbs_tmb%norb, orbs_tmb%norb, orbs_tmb%norbp, 1.d0, kernel(1,orbs_tmb%isorb+1), ld_coeff, &
 !           overlap(1,orbs_tmb%isorb+1), ld_coeff, 0.d0, ks(1,1), orbs_tmb%norb) 
 !call dgemm('n', 't', orbs_tmb%norb, orbs_tmb%norb, orbs_tmb%norbp, 1.d0, ks(1,orbs_tmb%isorb+1), ld_coeff, &
 !           kernel(1,orbs_tmb%isorb+1), ld_coeff, 0.d0, ksk(1,1), orbs_tmb%norb)
 !call dgemm('n', 't', orbs_tmb%norb, orbs_tmb%norb, orbs_tmb%norbp, 1.d0, ks(1,orbs_tmb%isorb+1), ld_coeff, &
 !           ksk(1,orbs_tmb%isorb+1), ld_coeff, 0.d0, ksksk(1,1), orbs_tmb%norb)


 !!if(present(overlap)) then
   !!allocate(ks(orbs_tmb%norb,orbs_tmb%norb),stat=istat)
   !!call memocc(istat, ks, 'ks', subname) 
   !!allocate(ksk(orbs_tmb%norb,orbs_tmb%norb),stat=istat)
   !!call memocc(istat, ksk, 'ksk', subname) 
   !!allocate(ksksk(orbs_tmb%norb,orbs_tmb%norb),stat=istat)
   !!call memocc(istat, ksksk, 'ksksk', subname) 

   !!call dgemm('n', 't', orbs_tmb%norb, orbs_tmb%norb, orbs_tmb%norb, 1.d0, kernel(1,1), ld_coeff, &
   !!           overlap(1,1), ld_coeff, 0.d0, ks(1,1), orbs_tmb%norb) 
   !!call dgemm('n', 't', orbs_tmb%norb, orbs_tmb%norb, orbs_tmb%norb, 1.d0, ks(1,1), ld_coeff, &
   !!           kernel(1,1), ld_coeff, 0.d0, ksk(1,1), orbs_tmb%norb)
   !!call dgemm('n', 't', orbs_tmb%norb, orbs_tmb%norb, orbs_tmb%norb, 1.d0, ks(1,1), ld_coeff, &
   !!           ksk(1,1), ld_coeff, 0.d0, ksksk(1,1), orbs_tmb%norb)
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
  type(collective_comms),intent(out) :: collcom_sr

  ! Local variables
  integer :: iorb, iiorb, ilr, ncount, is1, ie1, is2, ie2, is3, ie3, ii, i1, i2, i3, ierr, istat, iall, jproc, norb, ipt
  integer ::  ind, indglob, iitot, i, jproc_send, jproc_recv, i3e, jproc_out, p2p_tag
  integer,dimension(:),allocatable :: nsendcounts_tmp, nsenddspls_tmp, nrecvcounts_tmp, nrecvdspls_tmp, nsend, indexsendbuf
  integer,dimension(:),allocatable :: indexsendorbital
  integer,dimension(:),allocatable :: indexsendorbital2, indexrecvbuf
  integer,dimension(:),allocatable :: gridpoint_start, indexrecvorbital2
  real(kind=8) :: weight_tot, weight_ideal, tt, weightp, ttp
  integer,dimension(:,:),allocatable :: istartend
  character(len=*),parameter :: subname='determine_weights_sumrho'
  real(8) :: t1, t2, weight_start, weight_end, ttt
  real(kind=8),dimension(:),allocatable :: weights_per_slice
  real(kind=8),dimension(:,:),allocatable :: weights_startend

  ! Note: all weights are double precision to avoid integer overflow

  call timing(iproc,'init_collco_sr','ON')

  allocate(istartend(2,0:nproc-1), stat=istat)
  call memocc(istat, istartend, 'istartend', subname)
 
!!call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
!!t1=mpi_wtime()

  call get_weights_sumrho(iproc, nproc, orbs, lzd, nscatterarr, weight_tot, weight_ideal)

  call assign_weight_to_process_sumrho(iproc, nproc, weight_tot, weight_ideal, lzd, orbs, &
       nscatterarr, istartend, collcom_sr%nptsp_c)

!!call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
!!t2=mpi_wtime()
!!tt=t2-t1
!!if(iproc==0) write(*,*) 'time 1: iproc', iproc, tt
!!t1=mpi_wtime()

  allocate(collcom_sr%norb_per_gridpoint_c(collcom_sr%nptsp_c), stat=istat)
  call memocc(istat, collcom_sr%norb_per_gridpoint_c, 'collcom_sr%norb_per_gridpoint_c', subname)

  call determine_num_orbs_per_gridpoint_sumrho(iproc, nproc, collcom_sr%nptsp_c, lzd, orbs, &
       istartend, weight_tot, collcom_sr%norb_per_gridpoint_c)


!!call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
!!t2=mpi_wtime()
!!tt=t2-t1
!!if(iproc==0) write(*,*) 'time 2: iproc', iproc, tt
!!t1=mpi_wtime()

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

  allocate(collcom_sr%psit_c(collcom_sr%ndimind_c), stat=istat)
  call memocc(istat, collcom_sr%psit_c, 'collcom_sr%psit_c', subname)


!!call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
!!t2=mpi_wtime()
!!tt=t2-t1
!!if(iproc==0) write(*,*) 'time 3: iproc', iproc, tt
!!t1=mpi_wtime()


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


!!t2=mpi_wtime()
!!tt=t2-t1
!!if(iproc==0) write(*,*) 'time 4: iproc', iproc, tt
!!t1=mpi_wtime()

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


  !!! The tags for the self-made non blocking version of the mpi_alltoallv
  !!allocate(collcom_sr%tags(0:nproc-1), stat=istat)
  !!call memocc(istat, collcom_sr%tags, 'collcom_sr%tags', subname)
  !!do jproc=0,nproc-1
  !!    collcom_sr%tags(jproc)=p2p_tag(jproc)
  !!end do
  !!collcom_sr%messages_posted=.false.
  !!collcom_sr%communication_complete=.false.



  iall = -product(shape(istartend))*kind(istartend)
  deallocate(istartend,stat=istat)
  call memocc(istat, iall, 'istartend', subname)

!!call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
!!t2=mpi_wtime()
!!tt=t2-t1
!!if(iproc==0) write(*,*) 'time 5: iproc', iproc, tt

  call timing(iproc,'init_collco_sr','OF')

end subroutine init_collective_comms_sumro


subroutine get_weights_sumrho(iproc, nproc, orbs, lzd, nscatterarr, weight_tot, weight_ideal)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(orbitals_data),intent(in) :: orbs
  type(local_zone_descriptors),intent(in) :: lzd
  integer,dimension(0:nproc-1,4),intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  real(kind=8),intent(out) :: weight_tot, weight_ideal

  ! Local variables
  integer :: iorb, iiorb, ilr, ncount, ierr, i3, i2, i1, is1, ie1, is2, ie2, is3, ie3
  real(kind=8) :: tt

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


  weight_tot=0.d0
  !!$omp parallel default(shared) &
  !!$omp private(i2, i1, iorb, ilr, is1, ie1, is2, ie2, is3, ie3)
  do i3=nscatterarr(iproc,3)+1,nscatterarr(iproc,3)+nscatterarr(iproc,1)
      !!$omp do reduction(+:tt)
      do i2=1,lzd%glr%d%n2i
          do i1=1,lzd%glr%d%n1i
              tt=0.d0
              do iorb=1,orbs%norb
                  ilr=orbs%inwhichlocreg(iorb)
                  is1=1+lzd%Llr(ilr)%nsi1
                  ie1=lzd%Llr(ilr)%nsi1+lzd%llr(ilr)%d%n1i
                  is2=1+lzd%Llr(ilr)%nsi2
                  ie2=lzd%Llr(ilr)%nsi2+lzd%llr(ilr)%d%n2i
                  is3=1+lzd%Llr(ilr)%nsi3
                  ie3=lzd%Llr(ilr)%nsi3+lzd%llr(ilr)%d%n3i
                  if (is1<=i1 .and. i1<=ie1 .and. is2<=i2 .and. i2<=ie2 .and. is3<=i3 .and. i3<=ie3) then
                      tt=tt+1.d0
                  end if
              end do
              weight_tot=weight_tot+tt**2
          end do
      end do
      !!$omp end do
  end do
  !!$omp end parallel

  call mpiallred(weight_tot, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)

  ! Ideal weight per process
  weight_ideal = weight_tot/dble(nproc)

end subroutine get_weights_sumrho


subroutine assign_weight_to_process_sumrho(iproc, nproc, weight_tot, weight_ideal, lzd, orbs, &
           nscatterarr, istartend, nptsp)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  real(kind=8),intent(in) :: weight_tot, weight_ideal
  type(local_zone_descriptors),intent(in) :: lzd
  type(orbitals_data),intent(in) :: orbs
  integer,dimension(0:nproc-1,4),intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  integer,dimension(2,0:nproc-1),intent(out) :: istartend
  integer,intent(out) :: nptsp

  ! Local variables
  integer :: jproc, i1, i2, i3, ii, iorb, ilr, is1, ie1, is2, ie2, is3, ie3, ierr, istat, iall, jproc_out
  real(kind=8) :: tt, ttt
  real(8),dimension(:),allocatable :: weights_per_slice
  real(8),dimension(:,:),allocatable :: weights_startend
  character(len=*),parameter :: subname='assign_weight_to_process_sumrho'

  allocate(weights_per_slice(0:nproc-1), stat=istat)
  call memocc(istat, weights_per_slice, 'weights_per_slice', subname)

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

  call to_zero(nproc, weights_per_slice(0))
  ! Iterate through all grid points and assign them to processes such that the
  ! load balancing is optimal.

  !!do jproc=0,nproc-1
  !!    if (iproc==0) write(*,'(a,i7,2f16.1)') 'jproc, start, end', iproc, weights_startend(1,jproc), weights_startend(2,jproc)
  !!end do


  if (nproc>1) then
      tt=0.d0
      jproc=0
      istartend(1,jproc)=1
      !$omp parallel default(shared) &
      !$omp private(i2, i1, iorb, ilr, is1, ie1, is2, ie2, is3, ie3)
      do i3=nscatterarr(iproc,3)+1,nscatterarr(iproc,3)+nscatterarr(iproc,1)
          !$omp do reduction(+:tt)
          do i2=1,lzd%glr%d%n2i
              do i1=1,lzd%glr%d%n1i
                  ttt=0.d0
                  do iorb=1,orbs%norb
                      ilr=orbs%inwhichlocreg(iorb)
                      is1=1+lzd%Llr(ilr)%nsi1
                      ie1=lzd%Llr(ilr)%nsi1+lzd%llr(ilr)%d%n1i
                      is2=1+lzd%Llr(ilr)%nsi2
                      ie2=lzd%Llr(ilr)%nsi2+lzd%llr(ilr)%d%n2i
                      is3=1+lzd%Llr(ilr)%nsi3
                      ie3=lzd%Llr(ilr)%nsi3+lzd%llr(ilr)%d%n3i
                      if (is1<=i1 .and. i1<=ie1 .and. is2<=i2 .and. i2<=ie2 .and. is3<=i3 .and. i3<=ie3) then
                          !tt=tt+1.d0
                          ttt=ttt+1.d0
                      end if
                  end do
                  tt=tt+ttt**2
              end do
          end do
          !$omp end do
      end do
      !$omp end parallel
      weights_per_slice(iproc)=tt
      call mpiallred(weights_per_slice(0), nproc, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  end if



  ! Iterate through all grid points and assign them to processes such that the
  ! load balancing is optimal.
  if (nproc==1) then
      istartend(1,0)=1
      istartend(2,0)=lzd%glr%d%n1i*lzd%glr%d%n2i*lzd%glr%d%n3i
  else
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
          i3_loop: do i3=nscatterarr(jproc_out,3)+1,nscatterarr(jproc_out,3)+nscatterarr(jproc_out,1)
              do i2=1,lzd%glr%d%n2i
                  do i1=1,lzd%glr%d%n1i
                      ii=ii+1
                      ttt=0.d0
                      !$omp parallel if (orbs%norb>512) &
                      !$omp default(shared) &
                      !$omp private(iorb, ilr, is1, ie1, is2, ie2, is3, ie3)
                      !$omp do reduction(+:ttt)
                      do iorb=1,orbs%norb
                          ilr=orbs%inwhichlocreg(iorb)
                          is1=1+lzd%Llr(ilr)%nsi1
                          ie1=lzd%Llr(ilr)%nsi1+lzd%llr(ilr)%d%n1i
                          is2=1+lzd%Llr(ilr)%nsi2
                          ie2=lzd%Llr(ilr)%nsi2+lzd%llr(ilr)%d%n2i
                          is3=1+lzd%Llr(ilr)%nsi3
                          ie3=lzd%Llr(ilr)%nsi3+lzd%llr(ilr)%d%n3i
                          if (is1<=i1 .and. i1<=ie1 .and. is2<=i2 .and. i2<=ie2 .and. is3<=i3 .and. i3<=ie3) then
                              ttt=ttt+1.d0
                          end if
                      end do
                      !$omp end do
                      !$omp end parallel
                      !tt=tt+ttt
                      tt=tt+ttt**2
                      if (tt>=weights_startend(1,iproc)) then
                          istartend(1,iproc)=ii
                          exit outer_loop
                      end if
                  end do
              end do
          end do i3_loop
      end do outer_loop
  end if






  call mpiallred(istartend(1,0), 2*nproc, mpi_sum, bigdft_mpi%mpi_comm, ierr)


  do jproc=0,nproc-2
      istartend(2,jproc)=istartend(1,jproc+1)-1
  end do
  istartend(2,nproc-1)=lzd%glr%d%n1i*lzd%glr%d%n2i*lzd%glr%d%n3i

  !weightp=istartend(2,iproc)-istartend(1,iproc)+1


  do jproc=0,nproc-1
      if (iproc==jproc) then
          nptsp=istartend(2,jproc)-istartend(1,jproc)+1
      end if
  end do


  iall = -product(shape(weights_per_slice))*kind(weights_per_slice)
  deallocate(weights_per_slice,stat=istat)
  call memocc(istat, iall, 'weights_per_slice', subname)
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
           istartend, weight_tot, norb_per_gridpoint)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, nptsp
  type(local_zone_descriptors),intent(in) :: lzd
  type(orbitals_data),intent(in) :: orbs
  integer,dimension(2,0:nproc-1),intent(in) :: istartend
  real(kind=8),intent(in) :: weight_tot
  integer,dimension(nptsp),intent(out) :: norb_per_gridpoint

  ! Local variables
  integer :: i3, ii, i2, i1, ipt, norb, ilr, is1, ie1, is2, ie2, is3, ie3, iorb, ierr
  real(8) :: tt

  !!$omp parallel default(shared) &
  !!$omp private(i2, i1, ii, ipt, norb, iorb, ilr, is1, ie1, is2, ie2, is3, ie3)
  do i3=1,lzd%glr%d%n3i
      if (i3*lzd%glr%d%n1i*lzd%glr%d%n2i<istartend(1,iproc) .or. &
          (i3-1)*lzd%glr%d%n1i*lzd%glr%d%n2i+1>istartend(2,iproc)) then
          !!ii=ii+lzd%glr%d%n2i*lzd%glr%d%n1i
          cycle
      end if
      !!$omp do
      do i2=1,lzd%glr%d%n2i
          do i1=1,lzd%glr%d%n1i
              ii=(i3-1)*lzd%glr%d%n1i*lzd%glr%d%n2i+(i2-1)*lzd%glr%d%n1i+i1
              if (ii>=istartend(1,iproc) .and. ii<=istartend(2,iproc)) then
                  ipt=ii-istartend(1,iproc)+1
                  norb=0
                  do iorb=1,orbs%norb
                      ilr=orbs%inwhichlocreg(iorb)
                      is1=1+lzd%Llr(ilr)%nsi1
                      ie1=lzd%Llr(ilr)%nsi1+lzd%llr(ilr)%d%n1i
                      is2=1+lzd%Llr(ilr)%nsi2
                      ie2=lzd%Llr(ilr)%nsi2+lzd%llr(ilr)%d%n2i
                      is3=1+lzd%Llr(ilr)%nsi3
                      ie3=lzd%Llr(ilr)%nsi3+lzd%llr(ilr)%d%n3i
                      if (is1<=i1 .and. i1<=ie1 .and. is2<=i2 .and. i2<=ie2 .and. is3<=i3 .and. i3<=ie3) then
                          norb=norb+1.d0
                      end if
                  end do
                  norb_per_gridpoint(ipt)=norb
              end if
          end do
      end do
      !!$omp end do
  end do
  !!$omp end parallel
  !write(*,*) 'after loop', iproc

  ! Some check
  tt=0.d0
  do ipt=1,nptsp
      tt=tt+dble(norb_per_gridpoint(ipt))**2
  end do
  !tt=dble(sum(norb_per_gridpoint))
  call mpiallred(tt, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  write(*,*) 'tt, weight_tot', tt, weight_tot
  if (tt/=weight_tot) then
      stop '2: tt/=weight_tot'
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
  real(kind=8) :: tt
  character(len=*),parameter :: subname='determine_communication_arrays_sumrho'


  nsendcounts=0
  !$omp parallel default(shared) &
  !$omp private(jproc, i3, i2, i1, ind)
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      ilr=orbs%inwhichlocreg(iiorb)
      is1=1+lzd%Llr(ilr)%nsi1
      ie1=lzd%Llr(ilr)%nsi1+lzd%llr(ilr)%d%n1i
      is2=1+lzd%Llr(ilr)%nsi2
      ie2=lzd%Llr(ilr)%nsi2+lzd%llr(ilr)%d%n2i
      is3=1+lzd%Llr(ilr)%nsi3
      ie3=lzd%Llr(ilr)%nsi3+lzd%llr(ilr)%d%n3i
      !$omp do
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
       !$omp end do
  end do
  !$omp end parallel


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

  
  !$omp parallel default(shared) &
  !$omp private(iorb, iiorb, ilr, is1, ie1, is2, ie2, is3, ie3, i3, i2, i1, indglob, ind)
  !$omp do lastprivate(iitot)
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
  !$omp end do
  !$omp end parallel


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

  !$omp parallel default(shared) private(i, ind)
  !$omp do
  do i=1,ndimind
      ind=iextract(i)
      indexrecvorbital(ind)=indexrecvorbital2(i)
  end do
  !$omp end do
  !$omp end parallel

  iall=-product(shape(indexrecvorbital2))*kind(indexrecvorbital2)
  deallocate(indexrecvorbital2, stat=istat)
  call memocc(istat, iall, 'indexrecvorbital2', subname)

  if(minval(indexrecvorbital)<1) stop 'minval(indexrecvorbital)<1'
  if(maxval(indexrecvorbital)>orbs%norb) stop 'maxval(indexrecvorbital)>orbs%norb'




  iall=-product(shape(indexsendorbital))*kind(indexsendorbital)
  deallocate(indexsendorbital, stat=istat)
  call memocc(istat, iall, 'indexsendorbital', subname)
  iall=-product(shape(indexsendbuf))*kind(indexsendbuf)
  deallocate(indexsendbuf, stat=istat)
  call memocc(istat, iall, 'indexsendbuf', subname)
  iall=-product(shape(indexrecvbuf))*kind(indexrecvbuf)
  deallocate(indexrecvbuf, stat=istat)
  call memocc(istat, iall, 'indexrecvbuf', subname)

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



subroutine transpose_switch_psir(orbs, collcom_sr, psir, psirwork)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(orbitals_data),intent(in) :: orbs
  type(collective_comms),intent(in) :: collcom_sr
  real(kind=8),dimension(collcom_sr%ndimpsi_c),intent(in) :: psir
  real(kind=8),dimension(collcom_sr%ndimpsi_c),intent(out) :: psirwork

  ! Local variables
  integer :: i, m, ind


  !$omp parallel default(private) &
  !$omp shared(orbs, collcom_sr, psir, psirwork, m)

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
  integer :: ierr, istat, iall, ist, ist_c, ist_f, jproc, iisend, iirecv
  real(kind=8),dimension(:),allocatable :: psiwork, psitwork
  integer,dimension(:),allocatable :: nsendcounts, nsenddspls, nrecvcounts, nrecvdspls
  character(len=*),parameter :: subname='transpose_communicate_psit'

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



subroutine sumrho_for_TMBs(iproc, nproc, hx, hy, hz, orbs, collcom_sr, kernel, ndimrho, rho)
  use module_base
  use module_types
  use libxc_functionals
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, ndimrho
  real(kind=8),intent(in) :: hx, hy, hz
  type(orbitals_data),intent(in) :: orbs
  type(collective_comms),intent(in) :: collcom_sr
  real(kind=8),dimension(orbs%norb,orbs%norb),intent(in) :: kernel
  real(kind=8),dimension(ndimrho),intent(out) :: rho

  ! Local variables
  integer :: ipt, ii, i0, iiorb, jjorb, istat, iall, i, j, ierr
  real(8) :: tt, total_charge, hxh, hyh, hzh, factor, ddot, op
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
  
  ! Initialize rho.
  if (libxc_functionals_isgga()) then
      call razero(collcom_sr%nptsp_c, rho_local)
  else
      ! There is no mpi_allreduce, therefore directly initialize to
      ! 10^-20 and not 10^-20/nproc.
      rho_local=1.d-20
  end if

  if (iproc==0) write(*,'(a)', advance='no') 'Calculating charge density... '

  !!$omp parallel default(private) &
  !!$omp shared(total_charge, collcom_sr, factor, kernel, rho_local)

  total_charge=0.d0
  !!$omp do reduction(+:total_charge)
  op=0.d0
  do ipt=1,collcom_sr%nptsp_c
      ii=collcom_sr%norb_per_gridpoint_c(ipt)
      i0 = collcom_sr%isptsp_c(ipt)
      do i=1,ii
          iiorb=collcom_sr%indexrecvorbital_c(i0+i)
          do j=1,ii
              jjorb=collcom_sr%indexrecvorbital_c(i0+j)
              tt=factor*kernel(iiorb,jjorb)*collcom_sr%psit_c(i0+i)*collcom_sr%psit_c(i0+j)
              rho_local(ipt)=rho_local(ipt)+tt
              total_charge=total_charge+tt
              op=op+1.d0
          end do
      end do
  end do
  !!$omp end do
  !!$omp end parallel

  call mpi_allreduce(op, tt, 1, mpi_double_precision, mpi_min, bigdft_mpi%mpi_comm, ierr)
  if (iproc==0) write(*,'(a,es18.8)') 'minimal value', tt
  call mpi_allreduce(op, tt, 1, mpi_double_precision, mpi_max, bigdft_mpi%mpi_comm, ierr)
  if (iproc==0) write(*,'(a,es18.8)') 'maximal value', tt

  if (iproc==0) write(*,'(a)') 'done.'

  call timing(iproc,'sumrho_TMB    ','OF')


  call timing(iproc,'sumrho_allred','ON')

  call mpiallred(total_charge, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)

  if(iproc==0) write(*,'(3x,a,es20.12)') 'Calculation finished. TOTAL CHARGE = ', total_charge*hxh*hyh*hzh
  
  ! Communicate the density to meet the shape required by the Poisson solver.
  if (nproc>1) then
      call mpi_alltoallv(rho_local, collcom_sr%nsendcounts_repartitionrho, collcom_sr%nsenddspls_repartitionrho, &
                         mpi_double_precision, rho, collcom_sr%nrecvcounts_repartitionrho, &
                         collcom_sr%nrecvdspls_repartitionrho, mpi_double_precision, &
                         bigdft_mpi%mpi_comm, ierr)
  else
      call vcopy(ndimrho, rho_local(1), 1, rho(1), 1)
  end if

  call timing(iproc,'sumrho_allred','OF')


  iall=-product(shape(rho_local))*kind(rho_local)
  deallocate(rho_local, stat=istat)
  call memocc(istat, iall, 'rho_local', subname)


end subroutine sumrho_for_TMBs
