!> @file
!!  Routines to calculate electronic density from wavefunctions
!! @author
!!    Copyright (C) 2007-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Calculates the charge density by summing the square of all orbitals
!! Input: psi
!! Output: rho
!subroutine sumrhoLinear(iproc,nproc,Lzd,orbs,hxh,hyh,hzh,psi,rho,&
!     & nscatterarr,nspin,GPU,symObj,irrzon,phnons,rhodsc)
!  use module_base!, only: gp,dp,wp,ndebug,memocc
!  use module_types
!  use module_xc
!
!  implicit none
!  !Arguments
!  integer, intent(in) :: iproc,nproc,nspin,symObj
!  real(gp), intent(in) :: hxh,hyh,hzh
!  type(rho_descriptors),intent(in) :: rhodsc
!  type(orbitals_data),intent(in) :: orbs
!  type(local_zone_descriptors), intent(in) :: Lzd
!  integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
!  real(wp), dimension(Lzd%Lpsidimtot), intent(in) :: psi
!  real(dp), dimension(max(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*nscatterarr(iproc,1),1),nspin), intent(out), target :: rho
!  type(GPU_pointers), intent(inout) :: GPU
!  integer, dimension(*), intent(in) :: irrzon
!  real(dp), dimension(*), intent(in) :: phnons
!  !Local variables
!  character(len=*), parameter :: subname='sumrhoLinear'
!  integer :: nrhotot,n3d,itmred
!  integer :: nspinn
!  integer :: i1,i2,i3,i3off,i3s,i,ispin,jproc,i_all,i_stat,ierr,j3,j3p,j
!  real(dp) :: charge,tt, dasum
!  real(dp), dimension(:,:), allocatable :: tmred
!  real(dp), dimension(:,:), pointer :: rho_p
!!!  real(dp), dimension(:,:), allocatable :: rho_p_OCL
!!!  real(dp), dimension(:,:), allocatable :: psi_OCL
!
!!  integer :: ncount0,ncount1,ncount2,ncount3,ncountmpi0,ncountmpi1,ncount_max,ncount_rate
!  real(gp),dimension(:,:),allocatable :: dprho_comp
!  real(4) ,dimension(:,:),allocatable :: sprho_comp
!
!  call timing(iproc,'Rho_comput    ','ON')
!
!  if (iproc==0 .and. verbose >= 1) then
!     write(*,'(1x,a)',advance='no')&
!          'Calculation of charge density...'
!  end if
!
!  !components of the charge density
!  if (orbs%nspinor ==4) then
!     nspinn=4
!  else
!     nspinn=nspin
!  end if
!
!  !flag for toggling the REDUCE_SCATTER stategy (deprecated, icomm used instead of ixc value)
!  !rsflag=.not. ((ixc >= 11 .and. ixc <= 16) .or. &
!  !     & (ixc < 0 .and. module_xc_isgga()))
!
!!  write(*,*) 'RSFLAG stuffs ',(ixc >= 11 .and. ixc <= 16),&
!!             (ixc < 0 .and. module_xc_isgga()), have_mpi2,rsflag
!
!  !calculate dimensions of the complete array to be allocated before the reduction procedure
!  if (rhodsc%icomm==1) then
!     nrhotot=0
!     do jproc=0,nproc-1
!        nrhotot=nrhotot+nscatterarr(jproc,1)
!     end do
!  else
!     nrhotot=Lzd%Glr%d%n3i
!  end if
!
!  if (nproc > 1) then
!     allocate(rho_p(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*nrhotot,nspinn+ndebug),stat=i_stat)
!     call memocc(i_stat,rho_p,'rho_p',subname)
!  else
!     rho_p => rho
!  end if
!
!  !call system_clock(ncount1,ncount_rate,ncount_max)
!!!$  if (OCLconv) then
!!!$     allocate(rho_p_OCL(max(nrho,1),nspin),stat=i_stat)
!!!$     call memocc(i_stat,rho_p_OCL,'rho_p_OCL',subname)
!!!$     allocate(psi_OCL((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f),orbs%nspinor*orbs%norbp),stat=i_stat)
!!!$     call memocc(i_stat,psi_OCL,'psi_OCL',subname)
!!!$     call dcopy((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp, psi, 1, psi_OCL, 1)
!!!$  end if
!
!
!
!  !switch between GPU/CPU treatment of the density
!  if (GPUconv) then
!     if(Lzd%linear) stop 'local_partial_density_GPU not implemented with Linear scaling!'
!     call local_partial_density_GPU(orbs,nrhotot,Lzd%Glr,hxh,hyh,hzh,nspin,psi,rho_p,GPU)
!  else if (OCLconv) then
!     if(Lzd%linear) stop 'local_partial_density_OCL not implemented with Linear scaling!'
!     call local_partial_density_OCL(iproc,nproc,orbs,nrhotot,Lzd%Glr,hxh,hyh,hzh,nspin,psi,rho_p,GPU)
!  else if(Lzd%linear) then
!     call local_partial_densityLinear(iproc,nproc,(rhodsc%icomm==1),nscatterarr,nrhotot,&
!          Lzd,hxh,hyh,hzh,nspin,orbs,psi,rho_p)
!  else
!     !initialize the rho array at 10^-20 instead of zero, due to the invcb ABINIT routine
!     !otherwise use libXC routine
!     call xc_init_rho(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*nrhotot*nspinn,rho_p,nproc)
!
!     !for each of the orbitals treated by the processor build the partial densities
!     call local_partial_density(nproc,(rhodsc%icomm==1),nscatterarr,&
!          nrhotot,Lzd%Glr,hxh,hyh,hzh,nspin,orbs,psi,rho_p)
!  end if
!
!!!  if (OCLconv) then
!!!     call local_partial_density_OCL(iproc,nproc,orbs,nrhotot,Glr,hxh,hyh,hzh,nspin,psi_OCL,rho_p_OCL,GPU)
!!!     maxdiff=0.0_wp
!!!     do i=1,max(nrho,1)
!!!       do j=1,nspin
!!!        maxdiff=max(maxdiff,abs(rho_p(i,j)-rho_p_OCL(i,j)))
!!!       end do
!!!     end do
!!!     print *,''
!!!     print *,'maxdiff',maxdiff
!!!     i_all=-product(shape(rho_p_OCL))*kind(rho_p_OCL)
!!!     deallocate(rho_p_OCL,stat=i_stat)
!!!     call memocc(i_stat,i_all,'rho_p_OCL',subname)
!!!     i_all=-product(shape(psi_OCL))*kind(psi_OCL)
!!!     deallocate(psi_OCL,stat=i_stat)
!!!     call memocc(i_stat,i_all,'psi_OCL',subname)
!!!  end if
!
!  ! Symmetrise density, TODO...
!  !after validation this point can be deplaced after the allreduce such as to reduce the number of operations
!  if (symObj >= 0) then
!     call symmetrise_density(0,1,Lzd%Glr%d%n1i,Lzd%Glr%d%n2i,Lzd%Glr%d%n3i,nspin,&
!          Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i,&
!          rho_p,symObj,irrzon,phnons)
!  end if
!
!
!  !write(*,*) 'iproc,TIMING:SR1',iproc,real(ncount1-ncount0)/real(ncount_rate)
!  !the density must be communicated to meet the shape of the poisson solver
!  if (nproc > 1) then
!     call timing(iproc,'Rho_comput    ','OF')
!     call timing(iproc,'Rho_commun    ','ON')
!     !communication strategy for the density
!     !LDA case (icomm==1)
!     if (rhodsc%icomm==1) then
!        do ispin=1,nspin
!          call MPI_REDUCE_SCATTER(rho_p(1,ispin),rho(1,ispin),&
!               Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*nscatterarr(:,1),&
!               MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!        end do
!     ! splitted single-double precision communication (icomm=2)
!     else if (rhodsc%icomm==2) then
!        !        if (rho_compress .and. rhodsc%geocode.eq.'F') then
!        !if (rhodsc%geocode == 'F') then
!        !call system_clock(ncount0,ncount_rate,ncount_max)
!        !write(*,*) 'geocode=F, compress rho called'
!        allocate(sprho_comp(rhodsc%sp_size,nspin),stat=i_stat)
!        call memocc(i_stat,sprho_comp,'sprho_comp',subname)
!        allocate(dprho_comp(rhodsc%dp_size,nspin),stat=i_stat)
!        call memocc(i_stat,dprho_comp,'dprho_comp',subname)
!        call compress_rho(iproc,nproc,rho_p,Lzd%Glr,nspin,rhodsc,sprho_comp,dprho_comp)
!        !call system_clock(ncount1,ncount_rate,ncount_max)
!        !write(*,*) 'TIMING:ARED1',real(ncount1-ncount0)/real(ncount_rate)
!        call mpiallred(sprho_comp(1,1),rhodsc%sp_size*nspin,MPI_SUM,MPI_COMM_WORLD,ierr)
!        call mpiallred(dprho_comp(1,1),rhodsc%dp_size*nspin,MPI_SUM,MPI_COMM_WORLD,ierr)
!        !call system_clock(ncount2,ncount_rate,ncount_max)
!        !write(*,*) 'TIMING:ARED2',real(ncount2-ncount1)/real(ncount_rate)
!        i3s=nscatterarr(iproc,3)-nscatterarr(iproc,4)
!        n3d=nscatterarr(iproc,1)
!        call uncompress_rho(iproc,nproc,sprho_comp,dprho_comp,&
!             Lzd%Glr,nspin,rhodsc,rho_p,i3s,n3d)
!        !call system_clock(ncount3,ncount_rate,ncount_max)
!        !write(*,*) 'TIMING:ARED3',real(ncount3-ncount2)/real(ncount_rate)
!        i_all=-product(shape(sprho_comp))*kind(sprho_comp)
!        deallocate(sprho_comp,stat=i_stat)
!        call memocc(i_stat,i_all,'sprho_comp',subname)
!        i_all=-product(shape(dprho_comp))*kind(dprho_comp)
!        deallocate(dprho_comp,stat=i_stat)
!        call memocc(i_stat,i_all,'dprho_comp',subname)
!
!     !naive communication (unsplitted GGA case) (icomm=0)
!     else if (rhodsc%icomm==0) then
!        call mpiallred(rho_p(1,1),Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*nspin,&
!             MPI_SUM,MPI_COMM_WORLD,ierr)
!     else
!        STOP 'DENSITY COMMUNICATION KEY UNVALID' 
!     endif
!
!     call timing(iproc,'Rho_commun    ','OF')
!     call timing(iproc,'Rho_comput    ','ON')
!
!     if (rhodsc%icomm /= 1) then
!        !treatment which includes the periodic GGA
!        !the density should meet the poisson solver distribution
!        i3s=nscatterarr(iproc,3)-nscatterarr(iproc,4)
!        n3d=nscatterarr(iproc,1)
!        do ispin=1,nspin
!           do i3=1,n3d
!              j3=i3+i3s
!              j3p=modulo(j3-1,Lzd%Glr%d%n3i)+1
!              do i2=1,Lzd%Glr%d%n2i
!                 do i1=1,Lzd%Glr%d%n1i
!                    i=i1+(i2-1)*Lzd%Glr%d%n1i+Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*(i3-1)
!                    j=i1+(i2-1)*Lzd%Glr%d%n1i+Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*(j3p-1)
!                    rho(i,ispin)=rho_p(j,ispin)
!                 end do
!              end do
!           end do
!        end do
!     end if
!  end if
!  ! Check
!  tt=0.d0
!  i3off=Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*nscatterarr(iproc,4)
!
!  !allocation of the magnetic density orientation array
!  if (nproc > 1) then
!     itmred=2
!  else
!     itmred=1
!  end if
!
!  !use this check also for the magnetic density orientation
!  allocate(tmred(nspin+1,itmred+ndebug),stat=i_stat)
!  call memocc(i_stat,tmred,'tmred',subname)
!
!  tmred(nspin+1,itmred)=0.0_dp
!  do ispin=1,nspin!n
!     tmred(ispin,itmred)=0.0_dp
!     do i=1,Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*nscatterarr(iproc,2)
!!!        tt=tt+rho(i+i3off,ispin)
!        tmred(ispin,itmred)=tmred(ispin,itmred)+rho(i+i3off,ispin)
!        !temporary check for debugging purposes
!!!        if (rho(i+i3off,ispin)/rho(i+i3off,ispin) /= 1.d0) then
!!!           print *,iproc,'error in density construction',rho(i+i3off,ispin)
!!!        end if
!     enddo
!     tmred(nspin+1,itmred)=tmred(nspin+1,itmred)+tmred(ispin,itmred)
!  end do
!
!  if (nproc > 1) then
!     i_all=-product(shape(rho_p))*kind(rho_p)
!     deallocate(rho_p,stat=i_stat)
!     call memocc(i_stat,i_all,'rho_p',subname)
!
!     call timing(iproc,'Rho_comput    ','OF')
!     call timing(iproc,'Rho_commun    ','ON')
!!!     call MPI_REDUCE(tt,charge,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!     call MPI_REDUCE(tmred(1,2),tmred(1,1),nspin+1,mpidtypd,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!     call timing(iproc,'Rho_commun    ','OF')
!     call timing(iproc,'Rho_comput    ','ON')
!  else
!     !useless, only for completeness
!     nullify(rho_p)
!
!     !charge=tt
!  endif
!
!  !write the results
!  if (iproc == 0 .and. verbose >= 1) then
!     if(nspin==4) then
!        charge=tmred(1,1)
!        tt=sqrt(tmred(2,1)**2+tmred(3,1)**2+tmred(4,1)**2)
!     else
!        charge=0._dp
!        do ispin=1,nspin
!           charge=charge+tmred(ispin,1)
!        end do
!     end if
!     write(*,'(1x,a,f21.12)')&
!          'done. Total electronic charge=',real(charge,gp)*hxh*hyh*hzh
!     if(nspin == 4 .and. tt > 0._dp)&
!          write(*,'(a,5f10.4)')'  Magnetic density orientation:',&
!          (tmred(ispin,1)/tmred(1,1),ispin=2,nspin)
!  end if
!  
!  i_all=-product(shape(tmred))*kind(tmred)
!  deallocate(tmred,stat=i_stat)
!  call memocc(i_stat,i_all,'tmred',subname)
!
!  call timing(iproc,'Rho_comput    ','OF')
!
!END SUBROUTINE sumrhoLinear
!

!!>   Here starts the routine for building partial density inside the localisation region
!!!   This routine should be treated as a building-block for the linear scaling code
subroutine local_partial_densityLinear(iproc,nproc,rsflag,nscatterarr,&
     nrhotot,Lzd,hxh,hyh,hzh,nspin,orbs,psi,rho)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => local_partial_densityLinear
  use module_xc
  implicit none
  logical, intent(in) :: rsflag
  integer, intent(in) :: iproc,nproc
  integer,intent(inout):: nrhotot
  integer, intent(in) :: nspin
  real(gp), intent(in) :: hxh,hyh,hzh
  type(local_zone_descriptors), intent(in) :: Lzd
  type(orbitals_data),intent(in) :: orbs
  integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
  real(dp),dimension(max(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*nrhotot,1),max(nspin,orbs%nspinor)),intent(out):: rho
  !local variables
  character(len=*), parameter :: subname='local_partial_densityLinear'
  integer :: iorb,i_stat,i_all,ii, ind, indSmall, indLarge, orbtot
  integer :: oidx,sidx,nspinn,npsir,ncomplex, i1, i2, i3, ilr, ispin
  integer :: nspincomp,i3s,i3e
  real(gp) :: hfac,spinval
  type(workarr_sumrho) :: w
  real(wp), dimension(:,:), allocatable :: psir
  real(dp), dimension(:),allocatable :: rho_p
  real(8):: dnrm2
  integer, dimension(:,:), allocatable :: Lnscatterarr
  integer :: n3d,n3p,n3pi,i3xcsh,i3tmp

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
     hfac=orbs%kwgts(orbs%iokpt(ii))*(orbs%occup(iorb)/(hxh*hyh*hzh))
     spinval=orbs%spinsgn(iorb)

     if (hfac /= 0.d0) then

        !sum for complex function case, npsir=1 in that case
        do oidx=0,ncomplex

           do sidx=1,npsir
              call daub_to_isf(Lzd%Llr(ilr),w,psi(ind),psir(1,sidx))
              ind=ind+Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f
           end do

           select case(Lzd%Llr(ilr)%geocode)
           case('F')
              call partial_density_free((rsflag .and. .not. Lzd%linear),nproc,Lzd%Llr(ilr)%d%n1i,&
                   Lzd%Llr(ilr)%d%n2i,Lzd%Llr(ilr)%d%n3i,npsir,nspinn,Lzd%Llr(ilr)%d%n3i,&!nrhotot,&
                   hfac,Lnscatterarr,spinval,psir,rho_p,Lzd%Llr(ilr)%bounds%ibyyzz_r)

           case('P')

              call partial_density(rsflag,nproc,Lzd%Glr%d%n1i,Lzd%Glr%d%n2i,Lzd%Glr%d%n3i,&
                   npsir,nspinn,Lzd%Llr(ilr)%d%n3i,&!nrhotot,&
                   hfac,Lnscatterarr,spinval,psir,rho_p)

           case('S')

              call partial_density(rsflag,nproc,Lzd%Glr%d%n1i,Lzd%Glr%d%n2i,Lzd%Glr%d%n3i,&
                   npsir,nspinn,Lzd%Llr(ilr)%d%n3i,&!nrhotot,&
                   hfac,Lnscatterarr,spinval,psir,rho_p)

           end select

           ! Copy rho_p to the correct place in rho
           indSmall=0
           do ispin=1,nspinn
               do i3=1,Lzd%Llr(ilr)%d%n3i !min(Lzd%Llr(ilr)%d%n3i,nscatterarr(iproc,1)) 
                   if(Lzd%Llr(ilr)%nsi3 + i3 - 1 < 0) cycle   !throwing away the extra buffer of the locreg, related to the change of boundary conditions
                   if(Lzd%Llr(ilr)%nsi3 + i3  > Lzd%Glr%d%n3i) cycle
                   do i2=1,Lzd%Llr(ilr)%d%n2i
                       if(Lzd%Llr(ilr)%nsi2+ i2 - 1 < 0) cycle !same
                       if(Lzd%Llr(ilr)%nsi2 + i2  > Lzd%Glr%d%n2i) cycle
                       do i1=1,Lzd%Llr(ilr)%d%n1i
                           if(Lzd%Llr(ilr)%nsi1+ i1 - 1 < 0) cycle ! same
                           if(Lzd%Llr(ilr)%nsi1 + i1  > Lzd%Glr%d%n1i) cycle
                           ! indSmall is the index in the currect localization region
                           indSmall=indSmall+1
                           ! indLarge is the index in the whole box. 
                           indLarge=(Lzd%Llr(ilr)%nsi3+i3-1)*Lzd%Glr%d%n2i*Lzd%Glr%d%n1i +&
                               (Lzd%Llr(ilr)%nsi2+i2-1)*Lzd%Glr%d%n1i + Lzd%Llr(ilr)%nsi1+i1
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




subroutine sumrhoForLocalizedBasis2(iproc, nproc, orbs, Glr, input, lin, coeff, phi, nrho, rho, at, nscatterarr)
!
use module_base
use module_types
use libxc_functionals
use module_interfaces, exceptThisOne => sumrhoForLocalizedBasis2
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, nrho
type(orbitals_data),intent(in):: orbs
type(locreg_descriptors),intent(in):: Glr
type(input_variables),intent(in):: input
type(linearParameters),intent(inout):: lin
real(8),dimension(lin%lb%orbs%norb,orbs%norb),intent(in):: coeff
real(8),dimension(lin%lb%orbs%npsidim_orbs),intent(in):: phi
real(8),dimension(nrho),intent(out),target:: rho
type(atoms_data),intent(in):: at
integer, dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh

! Local variables
integer:: iorb, jorb, korb, istat, indLarge, i1, i2, i3, ilr, jlr
integer:: i1s, i1e, i2s, i2e, i3s, i3e, i1d, j1d, i2d, j2d, i3d, j3d, indri, indrj, ldim, iall, istr, istri, istrj
integer:: indi2, indi3, indj2, indj3, indl2, indl3, mpisource, mpidest, iiorb, jjorb
integer:: ierr, jproc, is, ie, nreceives
integer:: nfast, nslow, nsameproc, m, i1d0, j1d0, indri0, indrj0, indLarge0
real(8):: tt, hxh, hyh, hzh, factor, totalCharge, tt0, tt1, tt2, tt3, factorTimesDensKern, t1, t2, time
real(8),dimension(:,:),allocatable:: densKern
character(len=*),parameter:: subname='sumrhoForLocalizedBasis2'
integer,dimension(mpi_status_size):: stat
logical:: sendComplete, receiveComplete


if(iproc==0) write(*,'(1x,a)',advance='no') 'Calculating charge density...'

!lin%comsr%communComplete=.false.
!lin%comsr%computComplete=.false.


! Allocate the density kernel.
allocate(densKern(lin%lb%orbs%norb,lin%lb%orbs%norb), stat=istat)
call memocc(istat, densKern, 'densKern', subname)

call mpi_barrier(mpi_comm_world, ierr)
call cpu_time(t1)
! Calculate the density kernel.
call dgemm('n', 't', lin%lb%orbs%norb, lin%lb%orbs%norb, orbs%norb, 1.d0, coeff(1,1), lin%lb%orbs%norb, &
     coeff(1,1), lin%lb%orbs%norb, 0.d0, densKern(1,1), lin%lb%orbs%norb)
call mpi_barrier(mpi_comm_world, ierr)
call cpu_time(t2)
time=t2-t1
if(iproc==0) write(*,'(a,es12.4)') 'time for kernel:',time


! Define some constant factors.
hxh=.5d0*input%hx
hyh=.5d0*input%hy
hzh=.5d0*input%hz
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


! Check whether the communication has completed.
call mpi_barrier(mpi_comm_world, ierr)
call cpu_time(t1)
nfast=0
nsameproc=0
testLoop: do
    do jproc=0,nproc-1
        do korb=1,lin%comsr%noverlaps(jproc)
            if(lin%comsr%communComplete(korb,jproc)) cycle
            call mpi_test(lin%comsr%comarr(8,korb,jproc), sendComplete, stat, ierr)
            call mpi_test(lin%comsr%comarr(9,korb,jproc), receiveComplete, stat, ierr)
            ! Attention: mpi_test is a local function.
            if(sendComplete .and. receiveComplete) lin%comsr%communComplete(korb,jproc)=.true.
        end do
    end do
    ! If we made it until here, either all all the communication is
    ! complete or we better wait for each single orbital.
    exit testLoop
end do testLoop
call mpi_barrier(mpi_comm_world, ierr)
call cpu_time(t2)
time=t2-t1
if(iproc==0) write(*,'(a,es12.4)') 'time for test:',time

! Since mpi_test is a local function, check whether the communication has completed on all processes.
call mpi_barrier(mpi_comm_world, ierr)
call cpu_time(t1)
call mpiallred(lin%comsr%communComplete(1,0), nproc*maxval(lin%comsr%noverlaps), mpi_land, mpi_comm_world, ierr)
call mpi_barrier(mpi_comm_world, ierr)
call cpu_time(t2)
time=t2-t1
if(iproc==0) write(*,'(a,es12.4)') 'time for allreduce:',time



call mpi_barrier(mpi_comm_world, ierr)
call cpu_time(t1)
! Wait for the communications that have not completed yet
nslow=0
do jproc=0,nproc-1
    do korb=1,lin%comsr%noverlaps(jproc)
        if(lin%comsr%communComplete(korb,jproc)) then
            mpisource=lin%comsr%comarr(1,korb,jproc)
            mpidest=lin%comsr%comarr(5,korb,jproc)
            if(mpisource==mpidest) then
                nsameproc=nsameproc+1
            else
                nfast=nfast+1
            end if
            cycle
        end if
        !write(*,'(2(a,i0))') 'process ', iproc, ' is waiting for orbital ', korb
        nslow=nslow+1
        call mpi_wait(lin%comsr%comarr(8,korb,jproc), stat, ierr)
        call mpi_wait(lin%comsr%comarr(9,korb,jproc), stat, ierr)
        lin%comsr%communComplete(korb,jproc)=.true.
        lin%comsr%computComplete(korb,jproc)=.true.
    end do
end do
call mpi_barrier(mpi_comm_world, ierr)
call cpu_time(t2)
time=t2-t1
if(iproc==0) write(*,'(a,es12.4)') 'time for wait:',time
!call mpiallred(nreceives, 1, mpi_sum, mpi_comm_world, ierr)
!call mpiallred(nfast, 1, mpi_sum, mpi_comm_world, ierr)
!call mpiallred(nslow, 1, mpi_sum, mpi_comm_world, ierr)
!call mpiallred(nsameproc, 1, mpi_sum, mpi_comm_world, ierr)
if(iproc==0) write(*,'(1x,2(a,i0),a)') 'statistics: - ', nfast+nslow, ' point to point communications, of which ', &
                       nfast, ' could be overlapped with computation.'
if(iproc==0) write(*,'(1x,a,i0,a)') '            - ', nsameproc, ' copies on the same processor.'


do iorb=1,lin%comsr%noverlaps(iproc)
    if(.not. lin%comsr%communComplete(iorb,iproc)) then
        write(*,'(a,i0,a,i0,a)') 'ERROR: communication of orbital ', iorb, ' to process ', iproc, ' failed!'
        stop
    end if
    !!if(.not. lin%comsr%computComplete(iorb,iproc)) then
    !!    write(*,'(a,i0,a,i0,a)') 'ERROR: computation of orbital ', iorb, ' on process ', iproc, ' failed!'
    !!    stop
    !!end if
end do



! Now calculate the charge density. Each process calculates only one slice of the total charge density.
! Such a slice has the full extent in the x and y direction, but is limited in the z direction.
! The bounds of the slice are given by nscatterarr. To do so, each process has received all orbitals that
! extend into this slice. The number of these orbitals is given by lin%comsr%noverlaps(iproc).
call mpi_barrier(mpi_comm_world, ierr)
call cpu_time(t1)


! Bounds of the slice in global coordinates.
is=nscatterarr(iproc,3)-14
ie=is+nscatterarr(iproc,1)-1

totalCharge=0.d0
do iorb=1,lin%comsr%noverlaps(iproc)
    iiorb=lin%comsr%overlaps(iorb) !global index of orbital iorb
    ilr=lin%comsr%comarr(4,iorb,iproc) !localization region of orbital iorb
    istri=lin%comsr%comarr(6,iorb,iproc)-1 !starting index of orbital iorb in the receive buffer
    do jorb=1,lin%comsr%noverlaps(iproc)
        jjorb=lin%comsr%overlaps(jorb) !global indes of orbital jorb
        jlr=lin%comsr%comarr(4,jorb,iproc) !localization region of orbital jorb
        istrj=lin%comsr%comarr(6,jorb,iproc)-1 !starting index of orbital jorb in the receive buffer
        ! Bounds of the overlap of orbital iorb and jorb in global coordinates.
        i1s=max(2*lin%lzd%llr(ilr)%ns1-14,2*lin%lzd%llr(jlr)%ns1-14)
        i1e=min(2*lin%lzd%llr(ilr)%ns1-14+lin%lzd%llr(ilr)%d%n1i-1,2*lin%lzd%llr(jlr)%ns1-14+lin%lzd%llr(jlr)%d%n1i-1)
        i2s=max(2*lin%lzd%llr(ilr)%ns2-14,2*lin%lzd%llr(jlr)%ns2-14)
        i2e=min(2*lin%lzd%llr(ilr)%ns2-14+lin%lzd%llr(ilr)%d%n2i-1,2*lin%lzd%llr(jlr)%ns2-14+lin%lzd%llr(jlr)%d%n2i-1)
        i3s=max(2*lin%lzd%llr(ilr)%ns3-14,2*lin%lzd%llr(jlr)%ns3-14,is)
        i3e=min(2*lin%lzd%llr(ilr)%ns3-14+lin%lzd%llr(ilr)%d%n3i-1,2*lin%lzd%llr(jlr)%ns3-14+lin%lzd%llr(jlr)%d%n3i-1,ie)
        factorTimesDensKern = factor*densKern(iiorb,jjorb)
        ! Now loop over all points in the box in which the orbitals overlap.
        do i3=i3s,i3e !bounds in z direction
            !!i3d=i3-i3s+1 !z coordinate of orbital iorb with respect to the overlap box
            !!j3d=i3-i3s+1 !z coordinate of orbital jorb with respect to the overlap box
            i3d=i3-max(is,2*lin%lzd%llr(ilr)%ns3-14)+1 !z coordinate of orbital iorb with respect to the overlap box
            j3d=i3-max(is,2*lin%lzd%llr(jlr)%ns3-14)+1 !z coordinate of orbital jorb with respect to the overlap box
            indi3=(i3d-1)*lin%lzd%llr(ilr)%d%n2i*lin%lzd%llr(ilr)%d%n1i !z-part of the index of orbital iorb in the 1-dim receive buffer
            indj3=(j3d-1)*lin%lzd%llr(jlr)%d%n2i*lin%lzd%llr(jlr)%d%n1i !z-part of the index of orbital jorb in the 1-dim receive buffer
            indl3=(i3-is)*Glr%d%n2i*Glr%d%n1i !z-part of the index for which the charge density is beeing calculated
            do i2=i2s,i2e !bounds in y direction
                i2d=i2-2*lin%lzd%llr(ilr)%ns2 !y coordinate of orbital iorb with respect to the overlap box
                j2d=i2-2*lin%lzd%llr(jlr)%ns2 !y coordinate of orbital jorb with respect to the overlap box
                indi2=(i2d+15-1)*lin%lzd%llr(ilr)%d%n1i !y-part of the index of orbital iorb in the 1-dim receive buffer
                indj2=(j2d+15-1)*lin%lzd%llr(jlr)%d%n1i !y-part of the index of orbital jorb in the 1-dim receive buffer
                indl2=(i2+15-1)*Glr%d%n1i !y-part of the index for which the charge density is beeing calculated
                !!!! This is the old version.
                !!do i1=i1s,i1e
                !!    i1d=i1-2*lin%lzd%llr(ilr)%ns1
                !!    j1d=i1-2*lin%lzd%llr(jlr)%ns1
                !!    ! Now calculate the index in the boxes.
                !!    indri = indi3 + indi2 + i1d+15 + istri
                !!    indrj = indj3 + indj2 + j1d+15 + istrj
                !!    indLarge = indl3 + indl2 + i1+15
                !!    tt = factor*densKern(iiorb,jjorb)*lin%comsr%recvBuf(indri)*lin%comsr%recvBuf(indrj)
                !!    rho(indLarge) = rho(indLarge) + tt
                !!    totalCharge = totalCharge + tt
                !!end do
                ! #####################################################################
                ! This is the new version.
                m=mod(i1e-i1s+1,4)
                if(m/=0) then
                    ! The following five variables hold some intermediate results to speed up the code.
                    i1d0=-2*lin%lzd%llr(ilr)%ns1 
                    j1d0=-2*lin%lzd%llr(jlr)%ns1
                    indri0 = indi3 + indi2 + 15 + istri
                    indrj0 = indj3 + indj2 + 15 + istrj
                    indLarge0 = indl3 + indl2 + 15
                    do i1=i1s,i1s+m-1
                        i1d=i1d0+i1 !x coordinate of orbital iorb with respect to the overlap box
                        j1d=j1d0+i1 !x coordinate of orbital jorb with respect to the overlap box
                        indri = indri0 + i1d !index of orbital iorb in the 1-dim receive buffer
                        indrj = indrj0 + j1d !index of orbital jorb in the 1-dim receive buffer
                        indLarge = indLarge0 + i1 !index for which the charge density is beeing calculated
                        tt = factorTimesDensKern*lin%comsr%recvBuf(indri)*lin%comsr%recvBuf(indrj)
                        rho(indLarge) = rho(indLarge) + tt !update the charge density at point indLarge
                        totalCharge = totalCharge + tt !add the contribution to the total charge
                    end do
                end if
                ! This is the same again, this time with unrolled loops.
                if(i1e-i1s+1>4) then
                    i1d0=-2*lin%lzd%llr(ilr)%ns1
                    j1d0=-2*lin%lzd%llr(jlr)%ns1
                    indri0 = indi3 + indi2 + 15 + istri
                    indrj0 = indj3 + indj2 + 15 + istrj
                    indLarge0 = indl3 + indl2 + 15
                    do i1=i1s+m,i1e,4
                        i1d=i1d0+i1
                        j1d=j1d0+i1
                        indri = indri0 + i1d
                        indrj = indrj0 + j1d
                        indLarge = indLarge0 + i1
                        tt0 = factorTimesDensKern*lin%comsr%recvBuf(indri  )*lin%comsr%recvBuf(indrj  )
                        tt1 = factorTimesDensKern*lin%comsr%recvBuf(indri+1)*lin%comsr%recvBuf(indrj+1)
                        tt2 = factorTimesDensKern*lin%comsr%recvBuf(indri+2)*lin%comsr%recvBuf(indrj+2)
                        tt3 = factorTimesDensKern*lin%comsr%recvBuf(indri+3)*lin%comsr%recvBuf(indrj+3)
                        rho(indLarge  ) = rho(indLarge  ) + tt0
                        rho(indLarge+1) = rho(indLarge+1) + tt1
                        rho(indLarge+2) = rho(indLarge+2) + tt2
                        rho(indLarge+3) = rho(indLarge+3) + tt3
                        totalCharge = totalCharge + tt0 + tt1 + tt2 + tt3
                    end do
                end if
            end do
        end do
    end do
end do
call mpi_barrier(mpi_comm_world, ierr)
call cpu_time(t2)
time=t2-t1
if(iproc==0) write(*,'(a,es12.4)') 'time for large loop:',time

call mpiallred(totalCharge, 1, mpi_sum, mpi_comm_world, ierr)
if(iproc==0) write(*,'(1x,a,es20.12)') 'done. TOTAL CHARGE = ', totalCharge*hxh*hyh*hzh

iall=-product(shape(densKern))*kind(densKern)
deallocate(densKern, stat=istat)
call memocc(istat, iall, 'densKern', subname)


end subroutine sumrhoForLocalizedBasis2





subroutine sumrholinear_auxiliary(iproc, nproc, orbs, Glr, input, lin, coeff, phi, at, nscatterarr)
!
use module_base
use module_types
use libxc_functionals
use module_interfaces, exceptThisOne => sumrholinear_auxiliary
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data),intent(in):: orbs
type(locreg_descriptors),intent(in):: Glr
type(input_variables),intent(in):: input
type(linearParameters),intent(inout):: lin
real(8),dimension(lin%lb%orbs%norb,orbs%norb),intent(in):: coeff
real(8),dimension(max(lin%lb%orbs%npsidim_orbs,lin%lb%orbs%npsidim_comp)),intent(in):: phi
type(atoms_data),intent(in):: at
integer, dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh

! Local variables
integer:: iorb, jorb, korb, istat, indLarge, i1, i2, i3, ilr, jlr, ind
integer:: i1s, i1e, i2s, i2e, i3s, i3e, i1d, j1d, i2d, j2d, i3d, j3d, indri, indrj, ldim, iall, istr, istri, istrj
integer:: indi2, indi3, indj2, indj3, indl2, indl3, mpisource, mpidest, iiorb, jjorb
integer:: ierr, jproc, is, ie, nreceives
integer:: nfast, nslow, nsameproc, m, i1d0, j1d0, indri0, indrj0, indLarge0
real(8):: tt, hxh, hyh, hzh, factor, totalCharge, tt0, tt1, tt2, tt3, factorTimesDensKern, t1, t2, time
real(8),dimension(:,:),allocatable:: densKern
character(len=*),parameter:: subname='sumrhoForLocalizedBasis2'
integer,dimension(mpi_status_size):: stat
logical:: sendComplete, receiveComplete


if(iproc==0) write(*,'(1x,a)',advance='no') 'Calculating auxiliary arrays for charge density...'

!lin%comsr%communComplete=.false.
!lin%comsr%computComplete=.false.


!!! Allocate the density kernel.
!!allocate(densKern(lin%lb%orbs%norb,lin%lb%orbs%norb), stat=istat)
!!call memocc(istat, densKern, 'densKern', subname)
!!
!!call mpi_barrier(mpi_comm_world, ierr)
!!call cpu_time(t1)
!!! Calculate the density kernel.
!!call dgemm('n', 't', lin%lb%orbs%norb, lin%lb%orbs%norb, orbs%norb, 1.d0, coeff(1,1), lin%lb%orbs%norb, &
!!     coeff(1,1), lin%lb%orbs%norb, 0.d0, densKern(1,1), lin%lb%orbs%norb)
!!call mpi_barrier(mpi_comm_world, ierr)
!!call cpu_time(t2)
!!time=t2-t1
!!if(iproc==0) write(*,'(a,es12.4)') 'time for kernel:',time
!!
!!
!!! Define some constant factors.
!!hxh=.5d0*input%hx
!!hyh=.5d0*input%hy
!!hzh=.5d0*input%hz
!!if(input%nspin==1) then
!!    factor=2.d0/(hxh*hyh*hzh)
!!else
!!    factor=1.d0/(hxh*hyh*hzh)
!!end if
!!
!!! Initialize rho.
!!if (libxc_functionals_isgga()) then
!!    call razero(nrho, rho)
!!else
!!    ! There is no mpi_allreduce, therefore directly initialize to
!!    ! 10^-20 and not 10^-20/nproc.
!!    rho=1.d-20
!!    !call tenminustwenty(nrho, rho, nproc)
!!end if


! Check whether the communication has completed.
call mpi_barrier(mpi_comm_world, ierr)
call cpu_time(t1)
nfast=0
nsameproc=0
testLoop: do
    do jproc=0,nproc-1
        do korb=1,lin%comsr%noverlaps(jproc)
            if(lin%comsr%communComplete(korb,jproc)) cycle
            call mpi_test(lin%comsr%comarr(8,korb,jproc), sendComplete, stat, ierr)
            call mpi_test(lin%comsr%comarr(9,korb,jproc), receiveComplete, stat, ierr)
            ! Attention: mpi_test is a local function.
            if(sendComplete .and. receiveComplete) lin%comsr%communComplete(korb,jproc)=.true.
        end do
    end do
    ! If we made it until here, either all all the communication is
    ! complete or we better wait for each single orbital.
    exit testLoop
end do testLoop
call mpi_barrier(mpi_comm_world, ierr)
call cpu_time(t2)
time=t2-t1
if(iproc==0) write(*,'(a,es12.4)') 'time for test:',time

! Since mpi_test is a local function, check whether the communication has completed on all processes.
call mpi_barrier(mpi_comm_world, ierr)
call cpu_time(t1)
call mpiallred(lin%comsr%communComplete(1,0), nproc*maxval(lin%comsr%noverlaps), mpi_land, mpi_comm_world, ierr)
call mpi_barrier(mpi_comm_world, ierr)
call cpu_time(t2)
time=t2-t1
if(iproc==0) write(*,'(a,es12.4)') 'time for allreduce:',time



call mpi_barrier(mpi_comm_world, ierr)
call cpu_time(t1)
! Wait for the communications that have not completed yet
nslow=0
do jproc=0,nproc-1
    do korb=1,lin%comsr%noverlaps(jproc)
        if(lin%comsr%communComplete(korb,jproc)) then
            mpisource=lin%comsr%comarr(1,korb,jproc)
            mpidest=lin%comsr%comarr(5,korb,jproc)
            if(mpisource==mpidest) then
                nsameproc=nsameproc+1
            else
                nfast=nfast+1
            end if
            cycle
        end if
        !write(*,'(2(a,i0))') 'process ', iproc, ' is waiting for orbital ', korb
        nslow=nslow+1
        call mpi_wait(lin%comsr%comarr(8,korb,jproc), stat, ierr)
        call mpi_wait(lin%comsr%comarr(9,korb,jproc), stat, ierr)
        lin%comsr%communComplete(korb,jproc)=.true.
        lin%comsr%computComplete(korb,jproc)=.true.
    end do
end do
call mpi_barrier(mpi_comm_world, ierr)
call cpu_time(t2)
time=t2-t1
if(iproc==0) write(*,'(a,es12.4)') 'time for wait:',time
!call mpiallred(nreceives, 1, mpi_sum, mpi_comm_world, ierr)
!call mpiallred(nfast, 1, mpi_sum, mpi_comm_world, ierr)
!call mpiallred(nslow, 1, mpi_sum, mpi_comm_world, ierr)
!call mpiallred(nsameproc, 1, mpi_sum, mpi_comm_world, ierr)
if(iproc==0) write(*,'(1x,2(a,i0),a)') 'statistics: - ', nfast+nslow, ' point to point communications, of which ', &
                       nfast, ' could be overlapped with computation.'
if(iproc==0) write(*,'(1x,a,i0,a)') '            - ', nsameproc, ' copies on the same processor.'


do iorb=1,lin%comsr%noverlaps(iproc)
    if(.not. lin%comsr%communComplete(iorb,iproc)) then
        write(*,'(a,i0,a,i0,a)') 'ERROR: communication of orbital ', iorb, ' to process ', iproc, ' failed!'
        stop
    end if
    !!if(.not. lin%comsr%computComplete(iorb,iproc)) then
    !!    write(*,'(a,i0,a,i0,a)') 'ERROR: computation of orbital ', iorb, ' on process ', iproc, ' failed!'
    !!    stop
    !!end if
end do



! Now calculate the charge density. Each process calculates only one slice of the total charge density.
! Such a slice has the full extent in the x and y direction, but is limited in the z direction.
! The bounds of the slice are given by nscatterarr. To do so, each process has received all orbitals that
! extend into this slice. The number of these orbitals is given by lin%comsr%noverlaps(iproc).
call mpi_barrier(mpi_comm_world, ierr)
call cpu_time(t1)


! First determine the bounds of the auxiliary array

! Bounds of the slice in global coordinates.
is=nscatterarr(iproc,3)-14
ie=is+nscatterarr(iproc,1)-1



!allocate(comsr%auxarray(comsr%maxsize_auxarray,lin%comsr%noverlaps(iproc),lin%comsr%noverlaps(iproc), stat=istat)
!call memocc(istat, comsr%auxarray, 'comsr%auxarray', subname)




! Bounds of the slice in global coordinates.
is=nscatterarr(iproc,3)-14
ie=is+nscatterarr(iproc,1)-1

totalCharge=0.d0
ind=1
do iorb=1,lin%comsr%noverlaps(iproc)
    iiorb=lin%comsr%overlaps(iorb) !global index of orbital iorb
    ilr=lin%comsr%comarr(4,iorb,iproc) !localization region of orbital iorb
    istri=lin%comsr%comarr(6,iorb,iproc)-1 !starting index of orbital iorb in the receive buffer
    !do jorb=1,lin%comsr%noverlaps(iproc)
    do jorb=iorb,lin%comsr%noverlaps(iproc)
        jjorb=lin%comsr%overlaps(jorb) !global indes of orbital jorb
        jlr=lin%comsr%comarr(4,jorb,iproc) !localization region of orbital jorb
        istrj=lin%comsr%comarr(6,jorb,iproc)-1 !starting index of orbital jorb in the receive buffer
        ! Bounds of the overlap of orbital iorb and jorb in global coordinates.
        i1s=max(2*lin%lzd%llr(ilr)%ns1-14,2*lin%lzd%llr(jlr)%ns1-14)
        i1e=min(2*lin%lzd%llr(ilr)%ns1-14+lin%lzd%llr(ilr)%d%n1i-1,2*lin%lzd%llr(jlr)%ns1-14+lin%lzd%llr(jlr)%d%n1i-1)
        i2s=max(2*lin%lzd%llr(ilr)%ns2-14,2*lin%lzd%llr(jlr)%ns2-14)
        i2e=min(2*lin%lzd%llr(ilr)%ns2-14+lin%lzd%llr(ilr)%d%n2i-1,2*lin%lzd%llr(jlr)%ns2-14+lin%lzd%llr(jlr)%d%n2i-1)
        i3s=max(2*lin%lzd%llr(ilr)%ns3-14,2*lin%lzd%llr(jlr)%ns3-14,is)
        i3e=min(2*lin%lzd%llr(ilr)%ns3-14+lin%lzd%llr(ilr)%d%n3i-1,2*lin%lzd%llr(jlr)%ns3-14+lin%lzd%llr(jlr)%d%n3i-1,ie)
        !factorTimesDensKern = factor*densKern(iiorb,jjorb)
        ! Now loop over all points in the box in which the orbitals overlap.
        do i3=i3s,i3e !bounds in z direction
            !!i3d=i3-i3s+1 !z coordinate of orbital iorb with respect to the overlap box
            !!j3d=i3-i3s+1 !z coordinate of orbital jorb with respect to the overlap box
            i3d=i3-max(is,2*lin%lzd%llr(ilr)%ns3-14)+1 !z coordinate of orbital iorb with respect to the overlap box
            j3d=i3-max(is,2*lin%lzd%llr(jlr)%ns3-14)+1 !z coordinate of orbital jorb with respect to the overlap box
            indi3=(i3d-1)*lin%lzd%llr(ilr)%d%n2i*lin%lzd%llr(ilr)%d%n1i !z-part of the index of orbital iorb in the 1-dim receive buffer
            indj3=(j3d-1)*lin%lzd%llr(jlr)%d%n2i*lin%lzd%llr(jlr)%d%n1i !z-part of the index of orbital jorb in the 1-dim receive buffer
            indl3=(i3-is)*Glr%d%n2i*Glr%d%n1i !z-part of the index for which the charge density is beeing calculated
            do i2=i2s,i2e !bounds in y direction
                i2d=i2-2*lin%lzd%llr(ilr)%ns2 !y coordinate of orbital iorb with respect to the overlap box
                j2d=i2-2*lin%lzd%llr(jlr)%ns2 !y coordinate of orbital jorb with respect to the overlap box
                indi2=(i2d+15-1)*lin%lzd%llr(ilr)%d%n1i !y-part of the index of orbital iorb in the 1-dim receive buffer
                indj2=(j2d+15-1)*lin%lzd%llr(jlr)%d%n1i !y-part of the index of orbital jorb in the 1-dim receive buffer
                indl2=(i2+15-1)*Glr%d%n1i !y-part of the index for which the charge density is beeing calculated
                !!!! This is the old version.
                !!do i1=i1s,i1e
                !!    i1d=i1-2*lin%lzd%llr(ilr)%ns1
                !!    j1d=i1-2*lin%lzd%llr(jlr)%ns1
                !!    ! Now calculate the index in the boxes.
                !!    indri = indi3 + indi2 + i1d+15 + istri
                !!    indrj = indj3 + indj2 + j1d+15 + istrj
                !!    indLarge = indl3 + indl2 + i1+15
                !!    tt = factor*densKern(iiorb,jjorb)*lin%comsr%recvBuf(indri)*lin%comsr%recvBuf(indrj)
                !!    rho(indLarge) = rho(indLarge) + tt
                !!    totalCharge = totalCharge + tt
                !!end do
                ! #####################################################################
                ! This is the new version.
                m=mod(i1e-i1s+1,4)
                if(m/=0) then
                    ! The following five variables hold some intermediate results to speed up the code.
                    i1d0=-2*lin%lzd%llr(ilr)%ns1 
                    j1d0=-2*lin%lzd%llr(jlr)%ns1
                    indri0 = indi3 + indi2 + 15 + istri
                    indrj0 = indj3 + indj2 + 15 + istrj
                    indLarge0 = indl3 + indl2 + 15
                    do i1=i1s,i1s+m-1
                        i1d=i1d0+i1 !x coordinate of orbital iorb with respect to the overlap box
                        j1d=j1d0+i1 !x coordinate of orbital jorb with respect to the overlap box
                        indri = indri0 + i1d !index of orbital iorb in the 1-dim receive buffer
                        indrj = indrj0 + j1d !index of orbital jorb in the 1-dim receive buffer
                        indLarge = indLarge0 + i1 !index for which the charge density is beeing calculated
                        !tt = factorTimesDensKern*lin%comsr%recvBuf(indri)*lin%comsr%recvBuf(indrj)
                        lin%comsr%auxarray(ind)=lin%comsr%recvBuf(indri)*lin%comsr%recvBuf(indrj)
                        !rho(indLarge) = rho(indLarge) + tt !update the charge density at point indLarge
                        !totalCharge = totalCharge + tt !add the contribution to the total charge
                        ind=ind+1
                        !write(5000+iproc,*) indri, indrj
                    end do
                end if
                ! This is the same again, this time with unrolled loops.
                if(i1e-i1s+1>4) then
                    i1d0=-2*lin%lzd%llr(ilr)%ns1
                    j1d0=-2*lin%lzd%llr(jlr)%ns1
                    indri0 = indi3 + indi2 + 15 + istri
                    indrj0 = indj3 + indj2 + 15 + istrj
                    indLarge0 = indl3 + indl2 + 15
                    do i1=i1s+m,i1e,4
                        i1d=i1d0+i1
                        j1d=j1d0+i1
                        indri = indri0 + i1d
                        indrj = indrj0 + j1d
                        indLarge = indLarge0 + i1
                        !tt0 = factorTimesDensKern*lin%comsr%recvBuf(indri  )*lin%comsr%recvBuf(indrj  )
                        !tt1 = factorTimesDensKern*lin%comsr%recvBuf(indri+1)*lin%comsr%recvBuf(indrj+1)
                        !tt2 = factorTimesDensKern*lin%comsr%recvBuf(indri+2)*lin%comsr%recvBuf(indrj+2)
                        !tt3 = factorTimesDensKern*lin%comsr%recvBuf(indri+3)*lin%comsr%recvBuf(indrj+3)
                        lin%comsr%auxarray(ind  )=lin%comsr%recvBuf(indri  )*lin%comsr%recvBuf(indrj  )
                        lin%comsr%auxarray(ind+1)=lin%comsr%recvBuf(indri+1)*lin%comsr%recvBuf(indrj+1)
                        lin%comsr%auxarray(ind+2)=lin%comsr%recvBuf(indri+2)*lin%comsr%recvBuf(indrj+2)
                        lin%comsr%auxarray(ind+3)=lin%comsr%recvBuf(indri+3)*lin%comsr%recvBuf(indrj+3)
                        !rho(indLarge  ) = rho(indLarge  ) + tt0
                        !rho(indLarge+1) = rho(indLarge+1) + tt1
                        !rho(indLarge+2) = rho(indLarge+2) + tt2
                        !rho(indLarge+3) = rho(indLarge+3) + tt3
                        !totalCharge = totalCharge + tt0 + tt1 + tt2 + tt3
                        ind=ind+4
                        !write(5000+iproc,*) indri, indrj
                        !write(5000+iproc,*) indri+1, indrj+1
                        !write(5000+iproc,*) indri+2, indrj+2
                        !write(5000+iproc,*) indri+3, indrj+3
                    end do
                end if
            end do
        end do
    end do
end do
call mpi_barrier(mpi_comm_world, ierr)
call cpu_time(t2)
time=t2-t1
if(iproc==0) write(*,'(a,es12.4)') 'time for large loop:',time

!!call mpiallred(totalCharge, 1, mpi_sum, mpi_comm_world, ierr)
!!if(iproc==0) write(*,'(1x,a,es20.12)') 'done. TOTAL CHARGE = ', totalCharge*hxh*hyh*hzh
!!
!!iall=-product(shape(densKern))*kind(densKern)
!!deallocate(densKern, stat=istat)
!!call memocc(istat, iall, 'densKern', subname)


end subroutine sumrholinear_auxiliary






subroutine sumrholinear_withauxiliary(iproc, nproc, orbs, Glr, input, lin, coeff, nrho, rho, at, nscatterarr)
!
use module_base
use module_types
use libxc_functionals
use module_interfaces, exceptThisOne => sumrholinear_withauxiliary
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, nrho
type(orbitals_data),intent(in):: orbs
type(locreg_descriptors),intent(in):: Glr
type(input_variables),intent(in):: input
type(linearParameters),intent(inout):: lin
real(8),dimension(lin%lb%orbs%norb,orbs%norb),intent(in):: coeff
real(8),dimension(nrho),intent(out),target:: rho
type(atoms_data),intent(in):: at
integer, dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh

! Local variables
integer:: iorb, jorb, korb, istat, indLarge, i1, i2, i3, ilr, jlr
integer:: i1s, i1e, i2s, i2e, i3s, i3e, i1d, j1d, i2d, j2d, i3d, j3d, indri, indrj, ldim, iall, istr, istri, istrj
integer:: indi2, indi3, indj2, indj3, indl2, indl3, mpisource, mpidest, iiorb, jjorb
integer:: ierr, jproc, is, ie, nreceives, ind
integer:: nfast, nslow, nsameproc, m, i1d0, j1d0, indri0, indrj0, indLarge0
real(8):: tt, hxh, hyh, hzh, factor, totalCharge, tt0, tt1, tt2, tt3, factorTimesDensKern, t1, t2, time
real(8),dimension(:,:),allocatable:: densKern
character(len=*),parameter:: subname='sumrholinear_withauxiliary'
integer,dimension(mpi_status_size):: stat
logical:: sendComplete, receiveComplete


if(iproc==0) write(*,'(1x,a)',advance='no') 'Calculating charge density...'

!lin%comsr%communComplete=.false.
!lin%comsr%computComplete=.false.


! Allocate the density kernel.
allocate(densKern(lin%lb%orbs%norb,lin%lb%orbs%norb), stat=istat)
call memocc(istat, densKern, 'densKern', subname)

call mpi_barrier(mpi_comm_world, ierr)
call cpu_time(t1)
! Calculate the density kernel.
call dgemm('n', 't', lin%lb%orbs%norb, lin%lb%orbs%norb, orbs%norb, 1.d0, coeff(1,1), lin%lb%orbs%norb, &
     coeff(1,1), lin%lb%orbs%norb, 0.d0, densKern(1,1), lin%lb%orbs%norb)
call mpi_barrier(mpi_comm_world, ierr)
call cpu_time(t2)
time=t2-t1
if(iproc==0) write(*,'(a,es12.4)') 'time for kernel:',time


! Define some constant factors.
hxh=.5d0*input%hx
hyh=.5d0*input%hy
hzh=.5d0*input%hz
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


!!!! Check whether the communication has completed.
!!!call mpi_barrier(mpi_comm_world, ierr)
!!!call cpu_time(t1)
!!!nfast=0
!!!nsameproc=0
!!!testLoop: do
!!!    do jproc=0,nproc-1
!!!        do korb=1,lin%comsr%noverlaps(jproc)
!!!            if(lin%comsr%communComplete(korb,jproc)) cycle
!!!            call mpi_test(lin%comsr%comarr(8,korb,jproc), sendComplete, stat, ierr)
!!!            call mpi_test(lin%comsr%comarr(9,korb,jproc), receiveComplete, stat, ierr)
!!!            ! Attention: mpi_test is a local function.
!!!            if(sendComplete .and. receiveComplete) lin%comsr%communComplete(korb,jproc)=.true.
!!!        end do
!!!    end do
!!!    ! If we made it until here, either all all the communication is
!!!    ! complete or we better wait for each single orbital.
!!!    exit testLoop
!!!end do testLoop
!!!call mpi_barrier(mpi_comm_world, ierr)
!!!call cpu_time(t2)
!!!time=t2-t1
!!!if(iproc==0) write(*,'(a,es12.4)') 'time for test:',time
!!!
!!!! Since mpi_test is a local function, check whether the communication has completed on all processes.
!!!call mpi_barrier(mpi_comm_world, ierr)
!!!call cpu_time(t1)
!!!call mpiallred(lin%comsr%communComplete(1,0), nproc*maxval(lin%comsr%noverlaps), mpi_land, mpi_comm_world, ierr)
!!!call mpi_barrier(mpi_comm_world, ierr)
!!!call cpu_time(t2)
!!!time=t2-t1
!!!if(iproc==0) write(*,'(a,es12.4)') 'time for allreduce:',time
!!!
!!!
!!!
!!!call mpi_barrier(mpi_comm_world, ierr)
!!!call cpu_time(t1)
!!!! Wait for the communications that have not completed yet
!!!nslow=0
!!!do jproc=0,nproc-1
!!!    do korb=1,lin%comsr%noverlaps(jproc)
!!!        if(lin%comsr%communComplete(korb,jproc)) then
!!!            mpisource=lin%comsr%comarr(1,korb,jproc)
!!!            mpidest=lin%comsr%comarr(5,korb,jproc)
!!!            if(mpisource==mpidest) then
!!!                nsameproc=nsameproc+1
!!!            else
!!!                nfast=nfast+1
!!!            end if
!!!            cycle
!!!        end if
!!!        !write(*,'(2(a,i0))') 'process ', iproc, ' is waiting for orbital ', korb
!!!        nslow=nslow+1
!!!        call mpi_wait(lin%comsr%comarr(8,korb,jproc), stat, ierr)
!!!        call mpi_wait(lin%comsr%comarr(9,korb,jproc), stat, ierr)
!!!        lin%comsr%communComplete(korb,jproc)=.true.
!!!        lin%comsr%computComplete(korb,jproc)=.true.
!!!    end do
!!!end do
!!!call mpi_barrier(mpi_comm_world, ierr)
!!!call cpu_time(t2)
!!!time=t2-t1
!!!if(iproc==0) write(*,'(a,es12.4)') 'time for wait:',time
!!!!call mpiallred(nreceives, 1, mpi_sum, mpi_comm_world, ierr)
!!!!call mpiallred(nfast, 1, mpi_sum, mpi_comm_world, ierr)
!!!!call mpiallred(nslow, 1, mpi_sum, mpi_comm_world, ierr)
!!!!call mpiallred(nsameproc, 1, mpi_sum, mpi_comm_world, ierr)
!!!if(iproc==0) write(*,'(1x,2(a,i0),a)') 'statistics: - ', nfast+nslow, ' point to point communications, of which ', &
!!!                       nfast, ' could be overlapped with computation.'
!!!if(iproc==0) write(*,'(1x,a,i0,a)') '            - ', nsameproc, ' copies on the same processor.'
!!!
!!!
!!!do iorb=1,lin%comsr%noverlaps(iproc)
!!!    if(.not. lin%comsr%communComplete(iorb,iproc)) then
!!!        write(*,'(a,i0,a,i0,a)') 'ERROR: communication of orbital ', iorb, ' to process ', iproc, ' failed!'
!!!        stop
!!!    end if
!!!    !!if(.not. lin%comsr%computComplete(iorb,iproc)) then
!!!    !!    write(*,'(a,i0,a,i0,a)') 'ERROR: computation of orbital ', iorb, ' on process ', iproc, ' failed!'
!!!    !!    stop
!!!    !!end if
!!!end do



! Now calculate the charge density. Each process calculates only one slice of the total charge density.
! Such a slice has the full extent in the x and y direction, but is limited in the z direction.
! The bounds of the slice are given by nscatterarr. To do so, each process has received all orbitals that
! extend into this slice. The number of these orbitals is given by lin%comsr%noverlaps(iproc).
call mpi_barrier(mpi_comm_world, ierr)
call cpu_time(t1)


! Bounds of the slice in global coordinates.
is=nscatterarr(iproc,3)-14
ie=is+nscatterarr(iproc,1)-1

totalCharge=0.d0
do iorb=1,lin%comsr%noverlaps(iproc)
    iiorb=lin%comsr%overlaps(iorb) !global index of orbital iorb
    ilr=lin%comsr%comarr(4,iorb,iproc) !localization region of orbital iorb
    !istri=lin%comsr%comarr(6,iorb,iproc)-1 !starting index of orbital iorb in the receive buffer
    do jorb=1,lin%comsr%noverlaps(iproc)
        jjorb=lin%comsr%overlaps(jorb) !global indes of orbital jorb
        jlr=lin%comsr%comarr(4,jorb,iproc) !localization region of orbital jorb
        !istrj=lin%comsr%comarr(6,jorb,iproc)-1 !starting index of orbital jorb in the receive buffer
        ! Bounds of the overlap of orbital iorb and jorb in global coordinates.
        i1s=max(2*lin%lzd%llr(ilr)%ns1-14,2*lin%lzd%llr(jlr)%ns1-14)
        i1e=min(2*lin%lzd%llr(ilr)%ns1-14+lin%lzd%llr(ilr)%d%n1i-1,2*lin%lzd%llr(jlr)%ns1-14+lin%lzd%llr(jlr)%d%n1i-1)
        i2s=max(2*lin%lzd%llr(ilr)%ns2-14,2*lin%lzd%llr(jlr)%ns2-14)
        i2e=min(2*lin%lzd%llr(ilr)%ns2-14+lin%lzd%llr(ilr)%d%n2i-1,2*lin%lzd%llr(jlr)%ns2-14+lin%lzd%llr(jlr)%d%n2i-1)
        i3s=max(2*lin%lzd%llr(ilr)%ns3-14,2*lin%lzd%llr(jlr)%ns3-14,is)
        i3e=min(2*lin%lzd%llr(ilr)%ns3-14+lin%lzd%llr(ilr)%d%n3i-1,2*lin%lzd%llr(jlr)%ns3-14+lin%lzd%llr(jlr)%d%n3i-1,ie)
        factorTimesDensKern = factor*densKern(iiorb,jjorb)
        ! Now loop over all points in the box in which the orbitals overlap.
        if(jorb>=iorb) then
            ind=lin%comsr%startingindex(jorb,iorb)
        else
            ind=lin%comsr%startingindex(iorb,jorb)
        end if
        do i3=i3s,i3e !bounds in z direction
            !!i3d=i3-i3s+1 !z coordinate of orbital iorb with respect to the overlap box
            !!j3d=i3-i3s+1 !z coordinate of orbital jorb with respect to the overlap box
            !i3d=i3-max(is,2*lin%lzd%llr(ilr)%ns3-14)+1 !z coordinate of orbital iorb with respect to the overlap box
            !j3d=i3-max(is,2*lin%lzd%llr(jlr)%ns3-14)+1 !z coordinate of orbital jorb with respect to the overlap box
            !indi3=(i3d-1)*lin%lzd%llr(ilr)%d%n2i*lin%lzd%llr(ilr)%d%n1i !z-part of the index of orbital iorb in the 1-dim receive buffer
            !indj3=(j3d-1)*lin%lzd%llr(jlr)%d%n2i*lin%lzd%llr(jlr)%d%n1i !z-part of the index of orbital jorb in the 1-dim receive buffer
            indl3=(i3-is)*Glr%d%n2i*Glr%d%n1i !z-part of the index for which the charge density is beeing calculated
            do i2=i2s,i2e !bounds in y direction
                !i2d=i2-2*lin%lzd%llr(ilr)%ns2 !y coordinate of orbital iorb with respect to the overlap box
                !j2d=i2-2*lin%lzd%llr(jlr)%ns2 !y coordinate of orbital jorb with respect to the overlap box
                !indi2=(i2d+15-1)*lin%lzd%llr(ilr)%d%n1i !y-part of the index of orbital iorb in the 1-dim receive buffer
                !indj2=(j2d+15-1)*lin%lzd%llr(jlr)%d%n1i !y-part of the index of orbital jorb in the 1-dim receive buffer
                indl2=(i2+15-1)*Glr%d%n1i !y-part of the index for which the charge density is beeing calculated
                !!!! This is the old version.
                !!do i1=i1s,i1e
                !!    i1d=i1-2*lin%lzd%llr(ilr)%ns1
                !!    j1d=i1-2*lin%lzd%llr(jlr)%ns1
                !!    ! Now calculate the index in the boxes.
                !!    indri = indi3 + indi2 + i1d+15 + istri
                !!    indrj = indj3 + indj2 + j1d+15 + istrj
                !!    indLarge = indl3 + indl2 + i1+15
                !!    tt = factor*densKern(iiorb,jjorb)*lin%comsr%recvBuf(indri)*lin%comsr%recvBuf(indrj)
                !!    rho(indLarge) = rho(indLarge) + tt
                !!    totalCharge = totalCharge + tt
                !!end do
                ! #####################################################################
                ! This is the new version.
                m=mod(i1e-i1s+1,4)
                if(m/=0) then
                    ! The following five variables hold some intermediate results to speed up the code.
                    !i1d0=-2*lin%lzd%llr(ilr)%ns1 
                    !j1d0=-2*lin%lzd%llr(jlr)%ns1
                    !indri0 = indi3 + indi2 + 15 + istri
                    !indrj0 = indj3 + indj2 + 15 + istrj
                    indLarge0 = indl3 + indl2 + 15
                    do i1=i1s,i1s+m-1
                        !i1d=i1d0+i1 !x coordinate of orbital iorb with respect to the overlap box
                        !j1d=j1d0+i1 !x coordinate of orbital jorb with respect to the overlap box
                        !indri = indri0 + i1d !index of orbital iorb in the 1-dim receive buffer
                        !indrj = indrj0 + j1d !index of orbital jorb in the 1-dim receive buffer
                        indLarge = indLarge0 + i1 !index for which the charge density is beeing calculated
                        !tt = factorTimesDensKern*lin%comsr%recvBuf(indri)*lin%comsr%recvBuf(indrj)
                        tt = factorTimesDensKern*lin%comsr%auxarray(ind)
                        rho(indLarge) = rho(indLarge) + tt !update the charge density at point indLarge
                        totalCharge = totalCharge + tt !add the contribution to the total charge
                        ind=ind+1
                    end do
                end if
                ! This is the same again, this time with unrolled loops.
                if(i1e-i1s+1>4) then
                    !i1d0=-2*lin%lzd%llr(ilr)%ns1
                    !j1d0=-2*lin%lzd%llr(jlr)%ns1
                    !indri0 = indi3 + indi2 + 15 + istri
                    !indrj0 = indj3 + indj2 + 15 + istrj
                    indLarge0 = indl3 + indl2 + 15
                    do i1=i1s+m,i1e,4
                        !i1d=i1d0+i1
                        !j1d=j1d0+i1
                        !indri = indri0 + i1d
                        !indrj = indrj0 + j1d
                        indLarge = indLarge0 + i1
                        !!tt0 = factorTimesDensKern*lin%comsr%recvBuf(indri  )*lin%comsr%recvBuf(indrj  )
                        !!tt1 = factorTimesDensKern*lin%comsr%recvBuf(indri+1)*lin%comsr%recvBuf(indrj+1)
                        !!tt2 = factorTimesDensKern*lin%comsr%recvBuf(indri+2)*lin%comsr%recvBuf(indrj+2)
                        !!tt3 = factorTimesDensKern*lin%comsr%recvBuf(indri+3)*lin%comsr%recvBuf(indrj+3)
                        tt0 = factorTimesDensKern*lin%comsr%auxarray(ind  )
                        tt1 = factorTimesDensKern*lin%comsr%auxarray(ind+1)
                        tt2 = factorTimesDensKern*lin%comsr%auxarray(ind+2)
                        tt3 = factorTimesDensKern*lin%comsr%auxarray(ind+3)
                        rho(indLarge  ) = rho(indLarge  ) + tt0
                        rho(indLarge+1) = rho(indLarge+1) + tt1
                        rho(indLarge+2) = rho(indLarge+2) + tt2
                        rho(indLarge+3) = rho(indLarge+3) + tt3
                        totalCharge = totalCharge + tt0 + tt1 + tt2 + tt3
                        ind=ind+4
                    end do
                end if
            end do
        end do
    end do
end do
call mpi_barrier(mpi_comm_world, ierr)
call cpu_time(t2)
time=t2-t1
if(iproc==0) write(*,'(a,es12.4)') 'time for large loop:',time

call mpiallred(totalCharge, 1, mpi_sum, mpi_comm_world, ierr)
if(iproc==0) write(*,'(1x,a,es20.12)') 'done. TOTAL CHARGE = ', totalCharge*hxh*hyh*hzh

iall=-product(shape(densKern))*kind(densKern)
deallocate(densKern, stat=istat)
call memocc(istat, iall, 'densKern', subname)


end subroutine sumrholinear_withauxiliary









!> Initializes the parameters needed for the communication of the orbitals
!! when calculating the charge density.
!!
!! input arguments
!!  @param jproc        process to which the orbital shall be sent
!!  @param iorb         orbital that is to be sent
!!  @param istDest      the position on the MPI process to which it should be sent
!!  @param tag          communication tag
!!  @param lin          type containing the parameters for the linear scaling version
!! output arguments
!!  @param commsSumrho  contains the parameters
subroutine setCommunicationInformation(jproc, iorb, istDest, tag, lin, commsSumrho)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: jproc, iorb, istDest, tag
type(linearParameters),intent(in):: lin
integer,dimension(9),intent(out):: commsSumrho

! Local variables
integer:: mpisource, ist, jorb, jlr

! on which MPI process is the orbital that has to be sent to jproc
mpisource=lin%orbs%onWhichMPI(iorb)
commsSumrho(1)=mpisource

! starting index of the orbital on that MPI process
ist=1
do jorb=lin%orbs%isorb_par(mpisource)+1,iorb-1
    !jlr=lin%onWhichAtomAll(jorb)
    jlr=lin%orbs%inWhichLocreg(jorb)
    ist=ist+lin%lzd%llr(jlr)%wfd%nvctr_c+7*lin%lzd%llr(jlr)%wfd%nvctr_f
end do
commsSumrho(2)=ist

! amount of data to be sent
!jlr=lin%onWhichAtomAll(iorb)
jlr=lin%orbs%inWhichLocreg(iorb)
commsSumrho(3)=lin%lzd%llr(jlr)%wfd%nvctr_c+7*lin%lzd%llr(jlr)%wfd%nvctr_f

! localization region to which this orbital belongs to
!commsSumrho(4)=lin%onWhichAtomAll(iorb)
commsSumrho(4)=lin%orbs%inWhichLocreg(iorb)

! to which MPI process should this orbital be sent
commsSumrho(5)=jproc

! the position on the MPI process to which it should be sent
commsSumrho(6)=istDest

! the tag for this communication
commsSumrho(7)=tag

! commsSumrho(8): this entry is used a request for the mpi_isend.

! commsSumrho(9): this entry is used a request for the mpi_irecv.


end subroutine setCommunicationInformation




!> Initializes the parameters needed for the communication of the orbitals
!! when calculating the charge density.
!!
!! input arguments
!!  @param jproc        process to which the orbital shall be sent
!!  @param iorb         orbital that is to be sent
!!  @param istDest      the position on the MPI process to which it should be sent
!!  @param tag          communication tag
!!  @param lin          type containing the parameters for the linear scaling version
!! output arguments
!!  @param commsSumrho  contains the parameters
subroutine setCommunicationInformation2(jproc, iorb, is3ovrlp, n3ovrlp, istDest, tag, nlr, Llr, &
           onWhichAtomAll, orbs, commsSumrho)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: jproc, iorb, is3ovrlp, n3ovrlp, istDest, tag, nlr
type(locreg_descriptors),dimension(nlr),intent(in):: Llr
type(orbitals_data):: orbs
integer,dimension(orbs%norb),intent(in):: onWhichAtomAll
integer,dimension(9),intent(out):: commsSumrho

! Local variables
integer:: mpisource, ist, jorb, jlr

! on which MPI process is the orbital that has to be sent to jproc
mpisource=orbs%onWhichMPI(iorb)
commsSumrho(1)=mpisource

! starting index of the orbital on that MPI process
ist=1
do jorb=orbs%isorb_par(mpisource)+1,iorb-1
    jlr=onWhichAtomAll(jorb)
    !ist=ist+lin%lzd%llr(jlr)%wfd%nvctr_c+7*lin%lzd%llr(jlr)%wfd%nvctr_f
    ist = ist + Llr(jlr)%d%n1i*Llr(jlr)%d%n2i*Llr(jlr)%d%n3i
end do
jlr=onWhichAtomAll(iorb)
ist = ist + Llr(jlr)%d%n1i*Llr(jlr)%d%n2i*(is3ovrlp-1)
commsSumrho(2)=ist

! amount of data to be sent
jlr=onWhichAtomAll(iorb)
!commsSumrho(3)=lin%lzd%llr(jlr)%wfd%nvctr_c+7*lin%lzd%llr(jlr)%wfd%nvctr_f
commsSumrho(3)=Llr(jlr)%d%n1i*Llr(jlr)%d%n2i*n3ovrlp

! localization region to which this orbital belongs to
commsSumrho(4)=onWhichAtomAll(iorb)

! to which MPI process should this orbital be sent
commsSumrho(5)=jproc

! the position on the MPI process to which it should be sent
commsSumrho(6)=istDest

! the tag for this communication
commsSumrho(7)=tag

! commsSumrho(8): this entry is used as request for the mpi_isend.

! commsSumrho(9): this entry is used as request for the mpi_irecv.


end subroutine setCommunicationInformation2


!!subroutine setNscatterarr(iproc,nproc,datacode,atoms,n1i,n2i,n3i,ixc,nscatterarr)
!!  use module_base
!!  use module_types
!!  use Poisson_Solver
!!  use module_xc
!!  implicit none
!!  !Arguments
!!  character(len=1), intent(in) :: datacode
!!  integer, intent(in) :: iproc,nproc,ixc,n1i,n2i,n3i
!!  real(gp), intent(in) :: hxh,hyh,hzh
!!  type(atoms_data), intent(in) :: atoms
!!  integer, dimension(0:nproc-1,4), intent(out) :: nscatterarr
!!  !Local Variables
!!   integer :: jproc
!!
!!  if (datacode == 'D') then
!!     do jproc=0,iproc-1
!!        call PS_dim4allocation(atoms%geocode,datacode,jproc,nproc,n1i,n2i,n3i,ixc,n3d,n3p,n3pi,i3xcsh,i3s)
!!        nscatterarr(jproc,1)=n3d            !number of planes for the density
!!        nscatterarr(jproc,2)=n3p            !number of planes for the potential
!!        nscatterarr(jproc,3)=i3s+i3xcsh-1   !starting offset for the potential
!!        nscatterarr(jproc,4)=i3xcsh         !GGA XC shift between density and potential
!!     end do
!!     do jproc=iproc+1,nproc-1
!!        call PS_dim4allocation(atoms%geocode,datacode,jproc,nproc,n1i,n2i,n3i,ixc,n3d,n3p,n3pi,i3xcsh,i3s)
!!        nscatterarr(jproc,1)=n3d
!!        nscatterarr(jproc,2)=n3p
!!        nscatterarr(jproc,3)=i3s+i3xcsh-1
!!        nscatterarr(jproc,4)=i3xcsh
!!     end do
!!  end if
!!
!!  call PS_dim4allocation(atoms%geocode,datacode,iproc,nproc,n1i,n2i,n3i,ixc,n3d,n3p,n3pi,i3xcsh,i3s)
!!  nscatterarr(iproc,1)=n3d
!!  nscatterarr(iproc,2)=n3p
!!  nscatterarr(iproc,3)=i3s+i3xcsh-1
!!  nscatterarr(iproc,4)=i3xcsh
!!
!!end subroutine setNscatterarr
