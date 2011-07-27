!> @file 
!!    Calculate the electronic density (rho)
!! @author
!!   Copyright (C) 2007-2011 BigDFT group 
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 


!> Calculates the charge density by summing the square of all orbitals
!! Input: psi
!! Output: rho
subroutine sumrho(iproc,nproc,orbs,lr,hxh,hyh,hzh,psi,rho,&
     nscatterarr,nspin,GPU,symObj,irrzon,phnons,rhodsc)
  use module_base!, only: gp,dp,wp,ndebug,memocc
  use module_types
  use module_xc

  implicit none
  !Arguments
  integer, intent(in) :: iproc,nproc,nspin,symObj
  real(gp), intent(in) :: hxh,hyh,hzh
  type(rho_descriptors),intent(in) :: rhodsc
  type(orbitals_data), intent(in) :: orbs
  type(locreg_descriptors), intent(in) :: lr 
  integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%norbp*orbs%nspinor), intent(in) :: psi
  real(dp), dimension(max(lr%d%n1i*lr%d%n2i*nscatterarr(iproc,1),1),nspin), intent(out), target :: rho
  type(GPU_pointers), intent(inout) :: GPU
  integer, dimension(*), intent(in) :: irrzon
  real(dp), dimension(*), intent(in) :: phnons
  !Local variables
  character(len=*), parameter :: subname='sumrho'
  integer :: nrhotot,n3d,itmred
  integer :: nspinn
  integer :: i1,i2,i3,i3off,i3s,i,ispin,jproc,i_all,i_stat,ierr,j3,j3p,j
  real(dp) :: charge,tt
  real(dp), dimension(:,:), allocatable :: tmred
  real(dp), dimension(:,:), pointer :: rho_p
!!  real(dp), dimension(:,:), allocatable :: rho_p_OCL
!!  real(dp), dimension(:,:), allocatable :: psi_OCL

  !integer :: ncount0,ncount1,ncount2,ncount3,ncountmpi0,ncountmpi1,ncount_max,ncount_rate
  real(gp),dimension(:,:),allocatable :: dprho_comp
  real(4) ,dimension(:,:),allocatable :: sprho_comp

  call timing(iproc,'Rho_comput    ','ON')

  if (iproc==0 .and. verbose >= 1) then
     write(*,'(1x,a)',advance='no')&
          'Calculation of charge density...'
  end if

  !components of the charge density
  if (orbs%nspinor ==4) then
     nspinn=4
  else
     nspinn=nspin
  end if

  !flag for toggling the REDUCE_SCATTER stategy (deprecated, icomm used instead of ixc value)
  !rsflag=.not. ((ixc >= 11 .and. ixc <= 16) .or. &
  !     & (ixc < 0 .and. module_xc_isgga()))
  
!  write(*,*) 'RSFLAG stuffs ',(ixc >= 11 .and. ixc <= 16),&
!             (ixc < 0 .and. module_xc_isgga()), have_mpi2,rsflag

  !calculate dimensions of the complete array to be allocated before the reduction procedure
  if (rhodsc%icomm==1) then
     nrhotot=0
     do jproc=0,nproc-1
        nrhotot=nrhotot+nscatterarr(jproc,1)
     end do
  else
     nrhotot=lr%d%n3i
  end if

  if (nproc > 1) then
     !write(*,*) 'iproc,rhoarray dim', iproc, lr%d%n1i*lr%d%n2i*nrhotot,nspinn+ndebug
     allocate(rho_p(lr%d%n1i*lr%d%n2i*nrhotot,nspinn+ndebug),stat=i_stat)
     call memocc(i_stat,rho_p,'rho_p',subname)
  else
     rho_p => rho
  end if

  !call system_clock(ncount1,ncount_rate,ncount_max)
!!$  if (OCLconv) then
!!$     allocate(rho_p_OCL(max(nrho,1),nspin),stat=i_stat)
!!$     call memocc(i_stat,rho_p_OCL,'rho_p_OCL',subname)
!!$     allocate(psi_OCL((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),orbs%nspinor*orbs%norbp),stat=i_stat)
!!$     call memocc(i_stat,psi_OCL,'psi_OCL',subname)
!!$     call dcopy((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp, psi, 1, psi_OCL, 1)
!!$  end if


  !switch between GPU/CPU treatment of the density
  if (GPUconv) then
     call local_partial_density_GPU(iproc,nproc,orbs,nrhotot,lr,hxh,hyh,hzh,nspin,psi,rho_p,GPU)
  else if (OCLconv) then
     call local_partial_density_OCL(iproc,nproc,orbs,nrhotot,lr,hxh,hyh,hzh,nspin,psi,rho_p,GPU)
  else
     !initialize the rho array at 10^-20 instead of zero, due to the invcb ABINIT routine
     !otherwise use libXC routine
     call xc_init_rho(lr%d%n1i*lr%d%n2i*nrhotot*nspinn,rho_p,nproc)

     !for each of the orbitals treated by the processor build the partial densities
     call local_partial_density(iproc,nproc,(rhodsc%icomm==1),nscatterarr,&
          nrhotot,lr,hxh,hyh,hzh,nspin,orbs,psi,rho_p)
  end if

!!  if (OCLconv) then
!!     call local_partial_density_OCL(iproc,nproc,orbs,nrhotot,lr,hxh,hyh,hzh,nspin,psi_OCL,rho_p_OCL,GPU)
!!     maxdiff=0.0_wp
!!     do i=1,max(nrho,1)
!!       do j=1,nspin
!!        maxdiff=max(maxdiff,abs(rho_p(i,j)-rho_p_OCL(i,j)))
!!       end do
!!     end do
!!     print *,''
!!     print *,'maxdiff',maxdiff
!!     i_all=-product(shape(rho_p_OCL))*kind(rho_p_OCL)
!!     deallocate(rho_p_OCL,stat=i_stat)
!!     call memocc(i_stat,i_all,'rho_p_OCL',subname)
!!     i_all=-product(shape(psi_OCL))*kind(psi_OCL)
!!     deallocate(psi_OCL,stat=i_stat)
!!     call memocc(i_stat,i_all,'psi_OCL',subname)
!!  end if

  ! Symmetrise density, TODO...
  !after validation this point can be deplaced after the allreduce such as to reduce the number of operations
  if (symObj >= 0) then
     call symmetrise_density(0,1,lr%d%n1i,lr%d%n2i,lr%d%n3i,nscatterarr,nspin,&
          lr%d%n1i*lr%d%n2i*lr%d%n3i,&
          rho_p,symObj,irrzon,phnons)
  end if

  !write(*,*) 'iproc,TIMING:SR1',iproc,real(ncount1-ncount0)/real(ncount_rate)
  !the density must be communicated to meet the shape of the poisson solver
  if (nproc > 1) then
     call timing(iproc,'Rho_comput    ','OF')
     call timing(iproc,'Rho_commun    ','ON')
     !write(*,*) 'rsflag',rsflag
     !communication strategy for the density

     !LDA case (icomm==1)
     if (rhodsc%icomm==1) then
        do ispin=1,nspin
           call MPI_REDUCE_SCATTER(rho_p(1,ispin),rho(1,ispin),&
                lr%d%n1i*lr%d%n2i*nscatterarr(:,1),&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
           !write(*,*) 'LDA: MRC called'
        end do

     ! splitted single-double precision communication (icomm=2)
     else if (rhodsc%icomm==2) then
        !        if (rho_compress .and. rhodsc%geocode.eq.'F') then
        !if (rhodsc%geocode == 'F') then
        !call system_clock(ncount0,ncount_rate,ncount_max)
        !write(*,*) 'geocode=F, compress rho called'
        allocate(sprho_comp(rhodsc%sp_size,nspin),stat=i_stat)
        call memocc(i_stat,sprho_comp,'sprho_comp',subname)
        allocate(dprho_comp(rhodsc%dp_size,nspin),stat=i_stat)
        call memocc(i_stat,dprho_comp,'dprho_comp',subname)
        call compress_rho(iproc,nproc,rho_p,lr,nspin,rhodsc,sprho_comp,dprho_comp)
        !call system_clock(ncount1,ncount_rate,ncount_max)
        !write(*,*) 'TIMING:ARED1',real(ncount1-ncount0)/real(ncount_rate)
        call mpiallred(sprho_comp(1,1),rhodsc%sp_size*nspin,MPI_SUM,MPI_COMM_WORLD,ierr)
        call mpiallred(dprho_comp(1,1),rhodsc%dp_size*nspin,MPI_SUM,MPI_COMM_WORLD,ierr)
        !call system_clock(ncount2,ncount_rate,ncount_max)
        !write(*,*) 'TIMING:ARED2',real(ncount2-ncount1)/real(ncount_rate)
        i3s=nscatterarr(iproc,3)-nscatterarr(iproc,4)
        n3d=nscatterarr(iproc,1)
        call uncompress_rho(iproc,nproc,sprho_comp,dprho_comp,&
             lr,nspin,rhodsc,rho_p,i3s,n3d)
        !call system_clock(ncount3,ncount_rate,ncount_max)
        !write(*,*) 'TIMING:ARED3',real(ncount3-ncount2)/real(ncount_rate)
        i_all=-product(shape(sprho_comp))*kind(sprho_comp)
        deallocate(sprho_comp,stat=i_stat)
        call memocc(i_stat,i_all,'sprho_comp',subname)
        i_all=-product(shape(dprho_comp))*kind(dprho_comp)
        deallocate(dprho_comp,stat=i_stat)
        call memocc(i_stat,i_all,'dprho_comp',subname)

     !naive communication (unsplitted GGA case) (icomm=0)
     else if (rhodsc%icomm==0) then
        call mpiallred(rho_p(1,1),lr%d%n1i*lr%d%n2i*lr%d%n3i*nspin,&
             MPI_SUM,MPI_COMM_WORLD,ierr)
     else
        STOP 'DENSITY COMMUNICATION KEY UNVALID' 
     endif
     call timing(iproc,'Rho_commun    ','OF')
     call timing(iproc,'Rho_comput    ','ON')
     if (rhodsc%icomm /= 1) then
        !treatment which includes the periodic GGA
        !the density should meet the poisson solver distribution
        i3s=nscatterarr(iproc,3)-nscatterarr(iproc,4)
        n3d=nscatterarr(iproc,1)
        do ispin=1,nspin
           do i3=1,n3d
              j3=i3+i3s
              j3p=modulo(j3-1,lr%d%n3i)+1
              do i2=1,lr%d%n2i
                 do i1=1,lr%d%n1i
                    i=i1+(i2-1)*lr%d%n1i+lr%d%n1i*lr%d%n2i*(i3-1)
                    j=i1+(i2-1)*lr%d%n1i+lr%d%n1i*lr%d%n2i*(j3p-1)
                    rho(i,ispin)=rho_p(j,ispin)
                 end do
              end do
           end do
        end do
     end if
  end if

  ! Check
  tt=0.d0
  i3off=lr%d%n1i*lr%d%n2i*nscatterarr(iproc,4)

  !allocation of the magnetic density orientation array
  if (nproc > 1) then
     itmred=2
  else
     itmred=1
  end if

  !use this check also for the magnetic density orientation
  allocate(tmred(nspin+1,itmred+ndebug),stat=i_stat)
  call memocc(i_stat,tmred,'tmred',subname)

  tmred(nspin+1,itmred)=0.0_dp
  do ispin=1,nspin!n
     tmred(ispin,itmred)=0.0_dp
     do i=1,lr%d%n1i*lr%d%n2i*nscatterarr(iproc,2)
!!        tt=tt+rho(i+i3off,ispin)
        tmred(ispin,itmred)=tmred(ispin,itmred)+rho(i+i3off,ispin)
        !temporary check for debugging purposes
!!        if (rho(i+i3off,ispin)/rho(i+i3off,ispin) /= 1.d0) then
!!           print *,iproc,'error in density construction',rho(i+i3off,ispin)
!!        end if
     enddo
     tmred(nspin+1,itmred)=tmred(nspin+1,itmred)+tmred(ispin,itmred)
  end do

  if (nproc > 1) then
     i_all=-product(shape(rho_p))*kind(rho_p)
     deallocate(rho_p,stat=i_stat)
     call memocc(i_stat,i_all,'rho_p',subname)

     call timing(iproc,'Rho_comput    ','OF')
     call timing(iproc,'Rho_commun    ','ON')
!!     call MPI_REDUCE(tt,charge,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     call MPI_REDUCE(tmred(1,2),tmred(1,1),nspin+1,mpidtypd,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     call timing(iproc,'Rho_commun    ','OF')
     call timing(iproc,'Rho_comput    ','ON')
  else
     !useless, only for completeness
     nullify(rho_p)

     !charge=tt
  endif

  !write the results
  if (iproc == 0 .and. verbose >= 1) then
     if(nspin==4) then
        charge=tmred(1,1)
        tt=sqrt(tmred(2,1)**2+tmred(3,1)**2+tmred(4,1)**2)
     else
        charge=0._dp
        do ispin=1,nspin
           charge=charge+tmred(ispin,1)
        end do
     end if
     write(*,'(1x,a,f21.12)')&
          'done. Total electronic charge=',real(charge,gp)*hxh*hyh*hzh
     if(nspin == 4 .and. tt > 0._dp)&
          write(*,'(a,5f10.4)')'  Magnetic density orientation:',&
          (tmred(ispin,1)/tmred(1,1),ispin=2,nspin)
  end if
  
  i_all=-product(shape(tmred))*kind(tmred)
  deallocate(tmred,stat=i_stat)
  call memocc(i_stat,i_all,'tmred',subname)
  call timing(iproc,'Rho_comput    ','OF')

END SUBROUTINE sumrho


!> Here starts the routine for building partial density inside the localisation region
!! This routine should be treated as a building-block for the linear scaling code
subroutine local_partial_density(iproc,nproc,rsflag,nscatterarr,&
     nrhotot,lr,hxh,hyh,hzh,nspin,orbs,psi,rho_p)
  use module_base
  use module_types
  use module_interfaces
  implicit none
  logical, intent(in) :: rsflag
  integer, intent(in) :: iproc,nproc,nrhotot
  integer, intent(in) :: nspin
  real(gp), intent(in) :: hxh,hyh,hzh
  type(orbitals_data), intent(in) :: orbs
  type(locreg_descriptors), intent(in) :: lr
  integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(in) :: psi
  real(dp), dimension(lr%d%n1i,lr%d%n2i,nrhotot,max(nspin,orbs%nspinor)), intent(inout) :: rho_p
  !local variables
  character(len=*), parameter :: subname='local_partial_density'
  integer :: iorb,i_stat,i_all
  integer :: oidx,sidx,nspinn,npsir,ncomplex
  real(gp) :: hfac,spinval
  type(workarr_sumrho) :: w
  real(wp), dimension(:,:), allocatable :: psir

  call initialize_work_arrays_sumrho(lr,w)

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

  allocate(psir(lr%d%n1i*lr%d%n2i*lr%d%n3i,npsir+ndebug),stat=i_stat)
  call memocc(i_stat,psir,'psir',subname)
  !initialisation
  !print *,iproc,'there'
  if (lr%geocode == 'F') then
     call razero(lr%d%n1i*lr%d%n2i*lr%d%n3i*npsir,psir)
  end if

  do iorb=1,orbs%norbp
     !print *,'norbp',orbs%norbp,orbs%norb,orbs%nkpts,orbs%kwgts,orbs%iokpt,orbs%occup
     hfac=orbs%kwgts(orbs%iokpt(iorb))*(orbs%occup(orbs%isorb+iorb)/(hxh*hyh*hzh))
     spinval=orbs%spinsgn(orbs%isorb+iorb)

     if (hfac /= 0.d0) then

        !sum for complex function case, npsir=1 in that case
        do oidx=0,ncomplex

           do sidx=1,npsir
              call daub_to_isf(lr,w,psi(1,oidx+sidx,iorb),psir(1,sidx))
           end do

           !print *,'iorb,nrm',iorb,&
           !nrm2(lr%d%n1i*lr%d%n2i*lr%d%n3i*npsir,psir(1,1),1)

           select case(lr%geocode)
           case('F')

              call partial_density_free(rsflag,nproc,lr%d%n1i,lr%d%n2i,lr%d%n3i,&
                   npsir,nspinn,nrhotot,&
                   hfac,nscatterarr,spinval,psir,rho_p,lr%bounds%ibyyzz_r)

           case('P')

              call partial_density(rsflag,nproc,lr%d%n1i,lr%d%n2i,lr%d%n3i,&
                   npsir,nspinn,nrhotot,&
                   hfac,nscatterarr,spinval,psir,rho_p)

           case('S')

              call partial_density(rsflag,nproc,lr%d%n1i,lr%d%n2i,lr%d%n3i,&
                   npsir,nspinn,nrhotot,&
                   hfac,nscatterarr,spinval,psir,rho_p)

           end select

        end do
     end if
  enddo

  i_all=-product(shape(psir))*kind(psir)
  deallocate(psir,stat=i_stat)
  call memocc(i_stat,i_all,'psir',subname)
  call deallocate_work_arrays_sumrho(w)

END SUBROUTINE local_partial_density


!> Do partial density
subroutine partial_density(rsflag,nproc,n1i,n2i,n3i,npsir,nspinn,nrhotot,&
     hfac,nscatterarr,spinsgn,psir,rho_p)
  use module_base
  use module_types
  implicit none
  logical, intent(in) :: rsflag
  integer, intent(in) :: nproc,n1i,n2i,n3i,nrhotot,nspinn,npsir
  real(gp), intent(in) :: hfac,spinsgn
  integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr
  real(wp), dimension(n1i,n2i,n3i,npsir), intent(in) :: psir
  real(dp), dimension(n1i,n2i,nrhotot,nspinn), intent(inout) :: rho_p
  !local variables
  integer :: i3s,jproc,i3off,n3d,isjmp,i1,i2,i3,i1s,i1e,j3,i3sg
  real(gp) :: hfac2
  real(dp) :: psisq,p1,p2,p3,p4,r1,r2,r3,r4
!!!  integer :: ithread,nthread,omp_get_thread_num,omp_get_num_threads
  !sum different slices by taking into account the overlap
  i3sg=0
!$omp parallel default(private) shared(n1i,nproc,rsflag,nspinn,nscatterarr,spinsgn) &
!$omp shared(n2i,npsir,hfac,psir,rho_p,n3i,i3sg)
  i3s=0
  hfac2=2.0_gp*hfac
!!!  ithread=omp_get_thread_num()
!!!  nthread=omp_get_num_threads()

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
!!  if(mod(i3s,nthread) .eq. ithread) then
        !$omp do
        do i2=1,n2i
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
!!  end if
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

END SUBROUTINE partial_density


!> Do partial density for isolated system
subroutine partial_density_free(rsflag,nproc,n1i,n2i,n3i,npsir,nspinn,nrhotot,&
     hfac,nscatterarr,spinsgn,psir,rho_p,&
     ibyyzz_r) 
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
END SUBROUTINE partial_density_free


!> Symmetrise density
subroutine symmetrise_density(iproc,nproc,n1i,n2i,n3i,nscatterarr,nspin,nrho,rho,&
     symObj,irrzon,phnons)
  use module_base!, only: gp,dp,wp,ndebug,memocc
  use ab6_symmetry

  implicit none
  integer, intent(in) :: iproc,nproc,nrho,symObj,nspin, n1i, n2i, n3i
  integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  real(dp), dimension(n1i,n2i,n3i,nspin), intent(inout) :: rho
  integer, dimension(n1i*n2i*n3i,2,1), intent(in) :: irrzon 
  real(dp), dimension(2,n1i*n2i*n3i,1), intent(in) :: phnons 
  !local variables
  character(len=*), parameter :: subname='symmetrise_density'
  integer :: errno, ispden, nsym_used, nSym, isym, imagn, r2,i_stat,i_all,inzee,isign
  integer :: nd2, izone_max, numpt, izone, rep, nup, iup, ind, j, j1, j2, j3,i1,i2,i3
  real(dp) :: rhosu1, rhosu2
  real(dp), dimension(:,:), allocatable :: rhosu12
  real(dp), dimension(:,:,:,:,:), allocatable :: rhog
  integer, pointer :: sym(:,:,:)
  integer, pointer :: symAfm(:)
  real(gp), pointer :: transNon(:,:)

  call ab6_symmetry_get_matrices_p(symObj, nSym, sym, transNon, symAfm, errno)
  if (nSym == 1) return

!!$  ! Array sizes for the real-to-complex FFT: note that n1(there)=n1(here)+1
!!$  ! and the same for n2,n3. Not needed for the moment
!!$  call dimensions_fft(n1,n2,n3,nd1,nd2,nd3,n1f,n3f,n1b,n3b,nd1f,nd3f,nd1b,nd3b)

  !use this check also for the magnetic density orientation
  allocate(rhog(2,n1i+1,n2i+1,n3i+1,2+ndebug),stat=i_stat)
  call memocc(i_stat,rhog,'rhog',subname)


  ! Here imagn doesn't change since nspin == nsppol in BigDFT.
  imagn = 1

  !  Treat either full density, spin-up density or magnetization
  !  Note the decrease of ispden to the value 1, in order to finish
  !  with rhog of the total density (and not the spin-up density or magnetization)
  do ispden=nspin,1,-1

     !    Prepare the density to be symmetrized, in the reciprocal space
     nsym_used=0
     do isym=1,nSym
        if(symAfm(isym)==1)nsym_used=nsym_used+1
     end do

     !    rhor -fft-> rhog    (rhog is used as work space)
     !    Note : it should be possible to reuse rhog in the antiferromagnetic case
     !    this would avoid one FFT
     ! fft the input array x:
     do i3=0,n3i-1
        do i2=0,n2i-1
           do i1=0,n1i-1
              rhog(1,i1+1,i2+1,i3+1,1)=rho(i1+1,i2+1,i3+1,ispden)
              rhog(2,i1+1,i2+1,i3+1,1)=0.d0
           enddo
        enddo
     enddo

     inzee=1
     isign=-1

     call fft(n1i,n2i,n3i,n1i+1,n2i+1,n3i+1,rhog,isign,inzee)

!!     work(:)=rho(:,ispden)
!!     call fourdp(cplex,rhog,work,-1,mpi_enreg,nfft,ngfft,paral_kgb,0)

     !this should be modified via the nscatter array
     !for the moment we can put nproc=1
     nd2=n2i/nproc

     !    The following is only valid for total, up or dn density
     !    -------------------------------------------------------

     !    Get maxvalue of izone
     izone_max=count(irrzon(:,2,imagn)>0)
     allocate(rhosu12(2,izone_max+ndebug),stat=i_stat)
     call memocc(i_stat,rhosu12,'rhosu12',subname)

     numpt=0
     do izone=1,nrho

        !      Get repetition number
        rep=irrzon(izone,2,imagn)
        if(rep==0)exit

        !      Compute number of unique points in this symm class:
        nup=nsym_used/rep

        !      Accumulate charge over equivalent points
        rhosu1=0._dp
        rhosu2=0._dp
        do iup=1,nup
           ind=irrzon(iup+numpt,1,imagn)
           j=ind-1
           j1=modulo(j,n1i)
           j2=modulo(j/n1i,n2i)
           j3=j/(n1i*n2i)
           r2=modulo(j2,nd2)
           !here we should insert the condition that the planes should belong to iproc
           if(modulo(j/n1i,n2i)/nd2==iproc) then ! this ind is to be treated by me_fft
              ind=n1i*(nd2*j3+r2)+j1+1 !this is ind in the current proc
              rhosu1=rhosu1+rhog(1,j1+1,r2+1,j3+1,inzee)*phnons(1,iup+numpt,imagn)&
                   -rhog(2,j1+1,r2+1,j3+1,inzee)*phnons(2,iup+numpt,imagn)
              rhosu2=rhosu2+rhog(2,j1+1,r2+1,j3+1,inzee)*phnons(1,iup+numpt,imagn)&
                   +rhog(1,j1+1,r2+1,j3+1,inzee)*phnons(2,iup+numpt,imagn)
           end if

        end do
        rhosu1=rhosu1/real(nup,dp)
        rhosu2=rhosu2/real(nup,dp)
        rhosu12(1,izone)=rhosu1
        rhosu12(2,izone)=rhosu2
        !      Keep index of how many points have been considered:
        numpt=numpt+nup

        !      End loop over izone
     end do
     !reduction of the rho dimension to be discussed
     !call mpiallred(rhosu12(1,1),2*izone_max,MPI_SUM,MPI_COMM_WORLD,ierr)

     !    Reduction in case of FFT parallelization
!!     if(mpi_enreg%mode_para=='b')then
!!        old_paral_level=mpi_enreg%paral_level
!!        mpi_enreg%paral_level=3
!!        spaceComm=mpi_enreg%comm_fft
!!        call xsum_mpi(rhosu1_arr,spaceComm,ier)
!!        call xsum_mpi(rhosu2_arr,spaceComm,ier)
!!        mpi_enreg%paral_level=old_paral_level
!!     end if

     !    Now symmetrize the density
     numpt=0
     do izone=1,nrho

        !      Get repetition number
        rep=irrzon(izone,2,imagn)
        if(rep==0)exit

        !      Compute number of unique points in this symm class:
        nup=nsym_used/rep

        !      Define symmetrized rho(G) at equivalent points:
        do iup=1,nup
           ind=irrzon(iup+numpt,1,imagn)
           !        decompose ind-1=n1(n2 j3+ j2)+j1
           j=ind-1
           j1=modulo(j,n1i)
           j2=modulo(j/n1i,n2i)
           j3=j/(n1i*n2i)
           r2=modulo(j2,nd2)
           if(modulo(j/n1i,n2i)/nd2==iproc) then ! this ind is to be treated by me_fft
              !          ind in the proc ind-1=n1(nd2 j3+ r2)+j1
              ind=n1i*(nd2*j3+r2)+j1+1 !this is ind in the current proc
              rhog(1,j1+1,r2+1,j3+1,inzee)=rhosu12(1,izone)*phnons(1,iup+numpt,imagn)&
                   +rhosu12(2,izone)*phnons(2,iup+numpt,imagn)
              rhog(2,j1+1,r2+1,j3+1,inzee)=rhosu12(2,izone)*phnons(1,iup+numpt,imagn)&
                   -rhosu12(1,izone)*phnons(2,iup+numpt,imagn)
           end if
        end do

        !      Keep index of how many points have been considered:
        numpt=numpt+nup

        !      End loop over izone
     end do

     i_all=-product(shape(rhosu12))*kind(rhosu12)
     deallocate(rhosu12,stat=i_stat)
     call memocc(i_stat,i_all,'rhosu12',subname)

     !    Pull out full or spin up density, now symmetrized
!!     call fourdp(cplex,rhog,work,1,mpi_enreg,nfft,ngfft,paral_kgb,0)

     isign=1
     call fft(n1i,n2i,n3i,n1i+1,n2i+1,n3i+1,rhog,isign,inzee)

     do i3=0,n3i-1
        do i2=0,n2i-1
           do i1=0,n1i-1
              !correct the density in case it has negative values
              rho(i1+1,i2+1,i3+1,ispden)=max(rhog(1,i1+1,i2+1,i3+1,inzee)/real(n1i*n2i*n3i,dp),1.d-20)
           enddo
        enddo
     enddo
     !divide by the number of grid points
     !rho(:,ispden)=work(:)

  end do ! ispden

  i_all=-product(shape(rhog))*kind(rhog)
  deallocate(rhog,stat=i_stat)
  call memocc(i_stat,i_all,'rhog',subname)

END SUBROUTINE symmetrise_density


!> INPUT  
!!    rho_p: the partial rho array of the current proc
!!    spkey,dpkey: keys for coarse and fine regions
!! OUTPUT 
!!    sprho_comp, dprho_comp: compressed arrays of rho in single and double 
!!    precision
subroutine compress_rho(iproc,nprocs,rho_p,lr,nspin,rhodsc,sprho_comp,dprho_comp)
  use module_base
  use module_types
  implicit none
  type(locreg_descriptors), intent(in) :: lr 
  type(rho_descriptors),intent(in) :: rhodsc
  integer,intent(in) :: iproc,nspin,nprocs
  real(gp),dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,nspin),intent(in) :: rho_p
  real(gp),dimension(rhodsc%dp_size,nspin),intent(out) :: dprho_comp
  real(4),dimension(rhodsc%sp_size,nspin),intent(out) :: sprho_comp
  integer :: irho,jrho,iseg,ispin

  do ispin=1,nspin
    !$omp parallel default(shared) private(irho,jrho)
    !$omp do
    do iseg=1,rhodsc%n_csegs
      jrho=rhodsc%cseg_b(iseg)
      do irho=rhodsc%spkey(iseg,1),rhodsc%spkey(iseg,2)
        sprho_comp(jrho,ispin)=rho_p(irho,ispin)
        jrho=jrho+1
      enddo
    enddo
    !$omp enddo
    !$omp do
    do iseg=1,rhodsc%n_fsegs
      jrho=rhodsc%fseg_b(iseg)
      do irho=rhodsc%dpkey(iseg,1),rhodsc%dpkey(iseg,2)
        dprho_comp(jrho,ispin)=rho_p(irho,ispin)
        jrho=jrho+1
      enddo
    enddo
    !$omp enddo
    !$omp end parallel
  enddo
end subroutine compress_rho

! restore the necessary rho planes for using in GGA
subroutine uncompress_rho(iproc,nproc,sprho_comp,dprho_comp,&
    lr,nspin,rhodsc,rho_uncomp,i3s,n3d)
  use module_base
  use module_types
  implicit none
  integer,intent(in) :: iproc,nspin,nproc,i3s,n3d
  type(locreg_descriptors), intent(in) :: lr 
  type(rho_descriptors),intent(in) :: rhodsc
  real(gp),dimension(rhodsc%dp_size,nspin),intent(in) :: dprho_comp
  real(4), dimension(rhodsc%sp_size,nspin),intent(in) :: sprho_comp
  real(gp),dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,nspin),intent(out) :: rho_uncomp
  !local variables
  integer :: irho,jrho,ispin,iseg,ibegin,iend,j3p,imax
  logical :: overlap

  ! starting and endding point of rho array of interest
  imax=lr%d%n1i*lr%d%n2i*lr%d%n3i
  j3p=modulo(i3s,lr%d%n3i)+1
  ibegin=(j3p-1)*lr%d%n1i*lr%d%n2i+1
  j3p=modulo(i3s+n3d-1,lr%d%n3i)+1
  iend= j3p*lr%d%n1i*lr%d%n2i

  ! background 
  do ispin=1,nspin
    do irho=ibegin,iend
      rho_uncomp(irho,ispin)=1.d-20
    enddo
  enddo
  ! single precision region
  do ispin=1,nspin
    do iseg=1,rhodsc%n_csegs
      jrho=rhodsc%cseg_b(iseg)
      call is_overlap(ibegin,iend,rhodsc%spkey(iseg,1),rhodsc%spkey(iseg,2),overlap)
      if (overlap) then
        do irho=rhodsc%spkey(iseg,1),rhodsc%spkey(iseg,2)
           rho_uncomp(irho,ispin)=dble(sprho_comp(jrho,ispin))
          jrho=jrho+1
        enddo
      endif
    enddo
  enddo
  ! double precision region
  do ispin=1,nspin
    do iseg=1,rhodsc%n_fsegs
      jrho=rhodsc%fseg_b(iseg)
      call is_overlap(ibegin,iend,rhodsc%dpkey(iseg,1),rhodsc%dpkey(iseg,2),overlap)
      if (overlap) then
        do irho=rhodsc%dpkey(iseg,1),rhodsc%dpkey(iseg,2)
          rho_uncomp(irho,ispin)=dprho_comp(jrho,ispin)
          jrho=jrho+1
        enddo
      endif
    enddo
  enddo
end subroutine uncompress_rho

! restore the necessary rho planes for using in GGA
subroutine uncompress_rho_old(iproc,nproc,sprho_comp,dprho_comp,&
    lr,nspin,rhodsc,rho_uncomp,i3s,n3d)
  use module_base
  use module_types
  implicit none
  integer,intent(in) :: iproc,nspin,nproc,i3s,n3d
  type(locreg_descriptors), intent(in) :: lr 
  type(rho_descriptors),intent(in) :: rhodsc
  real(gp),dimension(rhodsc%dp_size,nspin),intent(in) :: dprho_comp
  real(4), dimension(rhodsc%sp_size,nspin),intent(in) :: sprho_comp
  real(gp),dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,nspin),intent(out) :: rho_uncomp
  !local variables
  integer :: irho,jrho,ispin,iseg,ibegin,iend,j3p,imax
  logical :: overflow,overlap

  ! starting and endding point of rho array of interest
  imax=lr%d%n1i*lr%d%n2i*lr%d%n3i
  j3p=modulo(i3s,lr%d%n3i)+1
  ibegin=(j3p-1)*lr%d%n1i*lr%d%n2i+1
  j3p=modulo(i3s+n3d-1,lr%d%n3i)+1
  iend= j3p*lr%d%n1i*lr%d%n2i
  if (ibegin.lt.iend) then 
    overflow=.false.
  else
    overflow=.true.
  endif
  ! background 
  do ispin=1,nspin
    do irho=ibegin,iend
      rho_uncomp(irho,ispin)=1.d-20
    enddo
  enddo
  ! single precision region
  do ispin=1,nspin
    do iseg=1,rhodsc%n_csegs
      jrho=rhodsc%cseg_b(iseg)
      if (overflow) then
        call is_overlap(ibegin,imax,rhodsc%spkey(iseg,1),rhodsc%spkey(iseg,2),overlap)
        if (overlap) then
          do irho=rhodsc%spkey(iseg,1),rhodsc%spkey(iseg,2)
             rho_uncomp(irho,ispin)=dble(sprho_comp(jrho,ispin))
            jrho=jrho+1
          enddo
        endif
        call is_overlap(1,iend,rhodsc%spkey(iseg,1),rhodsc%spkey(iseg,2),overlap)
        if (overlap) then
          do irho=rhodsc%spkey(iseg,1),rhodsc%spkey(iseg,2)
             rho_uncomp(irho,ispin)=dble(sprho_comp(jrho,ispin))
            jrho=jrho+1
          enddo
        endif
      else
        call is_overlap(ibegin,iend,rhodsc%spkey(iseg,1),rhodsc%spkey(iseg,2),overlap)
        if (overlap) then
          do irho=rhodsc%spkey(iseg,1),rhodsc%spkey(iseg,2)
             rho_uncomp(irho,ispin)=dble(sprho_comp(jrho,ispin))
            jrho=jrho+1
          enddo
        endif
      endif
    enddo
  enddo
  ! double precision region
  do ispin=1,nspin
    do iseg=1,rhodsc%n_fsegs
      jrho=rhodsc%fseg_b(iseg)
      if (overflow) then
        call is_overlap(ibegin,imax,rhodsc%dpkey(iseg,1),rhodsc%dpkey(iseg,2),overlap)
        if (overlap) then
          do irho=rhodsc%dpkey(iseg,1),rhodsc%dpkey(iseg,2)
            rho_uncomp(irho,ispin)=dprho_comp(jrho,ispin)
            jrho=jrho+1
          enddo
        endif
        call is_overlap(1,iend,rhodsc%dpkey(iseg,1),rhodsc%dpkey(iseg,2),overlap)
        if (overlap) then
          do irho=rhodsc%dpkey(iseg,1),rhodsc%dpkey(iseg,2)
            rho_uncomp(irho,ispin)=dprho_comp(jrho,ispin)
            jrho=jrho+1
          enddo
        endif
      else
        call is_overlap(ibegin,iend,rhodsc%dpkey(iseg,1),rhodsc%dpkey(iseg,2),overlap)
        if (overlap) then
          do irho=rhodsc%dpkey(iseg,1),rhodsc%dpkey(iseg,2)
            rho_uncomp(irho,ispin)=dprho_comp(jrho,ispin)
            jrho=jrho+1
          enddo
        endif
      endif
    enddo
  enddo
end subroutine uncompress_rho_old

subroutine rho_segkey(iproc,at,rxyz,crmult,frmult,radii_cf,&
    n1,n2,n3,n1i,n2i,n3i,hxh,hyh,hzh,nspin,rhodsc,iprint)
  use module_base
  use module_types
  implicit none
  integer,intent(in) :: n1,n2,n3,n1i,n2i,n3i,iproc,nspin
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(gp), intent(in) :: crmult,frmult,hxh,hyh,hzh
  real(gp), dimension(at%ntypes,3), intent(in) :: radii_cf
  logical,intent(in) :: iprint
  type(rho_descriptors),intent(inout) :: rhodsc
  !local variable
  integer :: i1,i2,i3,iseg,irho,i_stat,iat
  integer :: reg_c,reg_l
  integer(4),dimension(n1i*n2i*n3i) :: reg
  integer,dimension(n1i*n2i*n3i,2) :: dpkey,spkey
  integer :: n_fsegs,n_csegs
  character(len=*), parameter :: subname='rhokey'
  integer :: nbx,nby,nbz,nl1,nl2,nl3
  real(gp) :: dpmult,dsq,dsq_cr,dsq_fr,spadd
  integer :: i1min,i1max,i2min,i2max,i3min,i3max,nrhomin,nrhomax
  integer :: i1fmin,i1fmax,i2fmin,i2fmax,i3fmin,i3fmax
  integer :: i1cmin,i1cmax,i2cmin,i2cmax,i3cmin,i3cmax,csegstot,fsegstot,corx,cory,corz
  !integer :: ncount0,ncount1,ncount2,ncount3,ncount4,ncount_rate,ncount_max

  rhodsc%geocode=at%geocode
  !paprmeter to adjust the single precision and double precision regions
  spadd=10._gp
  dpmult=1.0_gp

  !call system_clock(ncount0,ncount_rate,ncount_max)

  ! calculate the corrections of the grid when transforming from 
  ! n1,n2,n3 to n1i, n2i, n3i
  call gridcorrection(nbx,nby,nbz,nl1,nl2,nl3,at%geocode)

  corx=nl1+nbx+1
  cory=nl2+nby+1
  corz=nl3+nbz+1 

  ! set the boundaries of the regions that we will determine the segment structure
  call get_boxbound(at,rxyz,radii_cf,crmult,frmult,hxh,hyh,hzh,spadd,dpmult,&
    n1i,n2i,n3i,corx,cory,corz,i1min,i1max,i2min,i2max,i3min,i3max)

  nrhomin = (i3min-1)*n1i*n2i+(i2min-1)*n1i+i1min
  nrhomax = (i3max-1)*n1i*n2i+(i2max-1)*n1i+i1max
  
  n_fsegs=0
  n_csegs=0

  do irho=nrhomin,nrhomax
    reg(irho)=0
  enddo
  !call system_clock(ncount1,ncount_rate,ncount_max)
  !write(*,*) 'TIMING:RHOKEY1',real(ncount1-ncount0)/real(ncount_rate)

  !$omp parallel default(none)&
  !$omp private(iat,irho,i1cmin,i1cmax,i2cmin,i2cmax,i3cmin)&
  !$omp private(i3cmax,i1fmin,i1fmax,i2fmin,i2fmax,i3fmin)&
  !$omp private(i3fmax,dsq,dsq_cr,dsq_fr,i1,i2,i3)&
  !$omp shared(at,rxyz,radii_cf,crmult,frmult,hxh,hyh,hzh)&
  !$omp shared(spadd,dpmult,n1i,n2i,n3i,corx,cory,corz,reg)
  !$omp do schedule(static,1)
  do iat=1,at%nat
    ! determine the bounds around the current atom within which
    ! search for single precision points done
    call get_atbound(iat,at,rxyz,radii_cf,crmult,frmult,hxh,&
      hyh,hzh,spadd,dpmult,n1i,n2i,n3i,corx,cory,corz,&
      i1cmin,i1cmax,i2cmin,i2cmax,i3cmin,i3cmax,&
      i1fmin,i1fmax,i2fmin,i2fmax,i3fmin,i3fmax)
    dsq_cr=(radii_cf(at%iatype(iat),1)*crmult+hxh*spadd)**2
    do i1=i1cmin,i1cmax
      do i2=i2cmin,i2cmax
        do i3=i3cmin,i3cmax
          dsq=(rxyz(1,iat)-(i1-corx)*hxh)**2+&
              (rxyz(2,iat)-(i2-cory)*hyh)**2+&
              (rxyz(3,iat)-(i3-corz)*hzh)**2          
          if(dsq.lt.dsq_cr) then
            irho = (i3-1)*n1i*n2i+(i2-1)*n1i+i1
            reg(irho)=1
          endif
        enddo
      enddo
    enddo
  enddo
  !$omp enddo

  !$omp do schedule(static,1)
  do iat=1,at%nat
    ! determine the bounds around the current atom within which
    ! search for double precision points done
    call get_atbound(iat,at,rxyz,radii_cf,crmult,frmult,hxh,hyh,hzh,&
      spadd,dpmult,n1i,n2i,n3i,corx,cory,corz,i1cmin,&
      i1cmax,i2cmin,i2cmax,i3cmin,i3cmax,i1fmin,i1fmax,i2fmin,&
      i2fmax,i3fmin,i3fmax)
    dsq_fr=((radii_cf(at%iatype(iat),2)*frmult*dpmult)**2)
    do i1=i1fmin,i1fmax
      do i2=i2fmin,i2fmax
        do i3=i3fmin,i3fmax
          dsq=(rxyz(1,iat)-(i1-corx)*hxh)**2+&
              (rxyz(2,iat)-(i2-cory)*hyh)**2+&
              (rxyz(3,iat)-(i3-corz)*hzh)**2          
          if(dsq.lt.dsq_fr) then
            irho = (i3-1)*n1i*n2i+(i2-1)*n1i+i1
            reg(irho)=2
          endif
        enddo
      enddo
    enddo
  enddo
  !$omp enddo
  !$omp end parallel
  !call system_clock(ncount2,ncount_rate,ncount_max)
  !write(*,*) 'TIMING:RHOKEY2',real(ncount2-ncount1)/real(ncount_rate)

  do irho=nrhomin,nrhomax
    if (irho.eq.nrhomin) then
      reg_c=reg(irho)
      select case (reg_c)
        case (2)
          n_fsegs=n_fsegs+1
          dpkey(n_fsegs,1)=irho
        case (1)
          n_csegs=n_csegs+1
          spkey(n_csegs,1)=irho
      end select
    else
      reg_c=reg(irho)
      reg_l=reg(irho-1)
      select case (reg_c)
        case (2)
          select case (reg_l)
            case (1)
              spkey(n_csegs,2)=irho-1
              n_fsegs=n_fsegs+1
              dpkey(n_fsegs,1)=irho
              if (irho.eq.nrhomax) then
                dpkey(n_fsegs,2)=irho
              endif
            case (2)
              if (irho.eq.nrhomax) then
                dpkey(n_fsegs,2)=irho
              endif
            case (0)
              n_fsegs=n_fsegs+1
              dpkey(n_fsegs,1)=irho
              if (irho.eq.nrhomax) then
                dpkey(n_fsegs,2)=irho
              endif
          end select
        case (1)
          select case (reg_l)
            case (2)
              dpkey(n_fsegs,2)=irho-1
              n_csegs=n_csegs+1
              spkey(n_csegs,1)=irho
              if (irho.eq.nrhomax) then
                spkey(n_csegs,2)=irho
              endif
            case (1)
              if (irho.eq.nrhomax) then
                spkey(n_csegs,2)=irho
              endif
            case (0)
              n_csegs=n_csegs+1
              spkey(n_csegs,1)=irho
              if (irho.eq.nrhomax) then
                spkey(n_csegs,2)=irho
              endif
          end select
        case (0)
          select case (reg_l)
            case (2)
              dpkey(n_fsegs,2)=irho-1
            case (1)
              spkey(n_csegs,2)=irho-1
          end select
      end select
    endif
  enddo

  !call system_clock(ncount3,ncount_rate,ncount_max)
  !write(*,*) 'TIMING:RHOKEY3',real(ncount3-ncount2)/real(ncount_rate)

  rhodsc%sp_size=0
  rhodsc%dp_size=0
  do iseg=1,n_csegs
    rhodsc%sp_size=rhodsc%sp_size+spkey(iseg,2)-spkey(iseg,1)+1
  enddo
  do iseg=1,n_fsegs
    rhodsc%dp_size=rhodsc%dp_size+dpkey(iseg,2)-dpkey(iseg,1)+1
  enddo

  rhodsc%n_fsegs=n_fsegs
  rhodsc%n_csegs=n_csegs

  allocate(rhodsc%dpkey(n_fsegs,2),stat=i_stat)
  call memocc(i_stat,rhodsc%dpkey,'dpkey',subname)
  allocate(rhodsc%spkey(n_csegs,2),stat=i_stat)
  call memocc(i_stat,rhodsc%spkey,'spkey',subname)
  allocate(rhodsc%cseg_b(n_csegs),stat=i_stat)
  call memocc(i_stat,rhodsc%cseg_b,'csegb',subname)
  allocate(rhodsc%fseg_b(n_fsegs),stat=i_stat)
  call memocc(i_stat,rhodsc%fseg_b,'fsegb',subname)

  csegstot=1
  fsegstot=1
  do iseg=1,n_fsegs
    rhodsc%fseg_b(iseg)=fsegstot
    rhodsc%dpkey(iseg,1)=dpkey(iseg,1)
    rhodsc%dpkey(iseg,2)=dpkey(iseg,2)
    fsegstot=fsegstot+dpkey(iseg,2)-dpkey(iseg,1)+1
  enddo
  do iseg=1,n_csegs
    rhodsc%cseg_b(iseg)=csegstot
    rhodsc%spkey(iseg,1)=spkey(iseg,1)
    rhodsc%spkey(iseg,2)=spkey(iseg,2)
    csegstot=csegstot+spkey(iseg,2)-spkey(iseg,1)+1
  enddo

  if (iprint) then
   if (iproc.eq.0) then
     !open(unit=1001,file='csegs.dat',status='unknown')
     !write(1001,*) n_csegs
     !do iseg=1,n_csegs
     ! write(1001,'(3(I5,1X))') iseg,rhodsc%spkey(iseg,1),rhodsc%spkey(iseg,2)
     !enddo
     !close(1001)
     !open(unit=1001,file='fsegs.dat',status='unknown')
     !write(1001,*) n_fsegs
     !do iseg=1,n_fsegs
     ! write(1001,'(3(I5,1X))') iseg,rhodsc%dpkey(iseg,1),rhodsc%dpkey(iseg,2)
     !enddo
     !close(1001)
     open(unit=1001,file='rhoseg.inf',status='unknown')
       write(1001,*) 'INFORMATION FROM RHO COMPRESSION PROCEDURE'
       write(1001,*) '----------------------------------------------------------------------'
       write(1001,'(1X,A,1X,I2,1X,2(F6.2,1X))') &
                     'nspin,spadd,dpmult                                 :',nspin,spadd,dpmult
       write(1001,*) 'spadd and dpmult can be found in "/src/sumrho.f90"'
       write(1001,*) '----------------------------------------------------------------------'
       write(1001,*) 'Number of single precision data points              :',&
                     rhodsc%sp_size*nspin
       write(1001,*) 'Number of double precision data points              :',&
                     rhodsc%dp_size*nspin
       write(1001,*) 'Total number of data points                         :',&
                     n1i*n2i*n3i*nspin
       write(1001,'(1X,A,1X,F4.2)') &
                     'Estimated compression ratio, number of data points  :',&
                     real(rhodsc%sp_size*nspin  +rhodsc%dp_size*nspin)/(n1i*n2i*n3i*nspin)
       write(1001,'(1X,A,1X,F4.2)') &
                     'Estimated compression ratio, data volume to be sent :',&
                     real(rhodsc%sp_size*nspin/2+rhodsc%dp_size*nspin)/(n1i*n2i*n3i*nspin)
       write(1001,*) ''
     close(1001)
   endif
 endif
 !call system_clock(ncount4,ncount_rate,ncount_max)
 !write(*,*) 'TIMING:RHOKEY4',real(ncount4-ncount3)/real(ncount_rate)
 !write(*,*) 'TIMING:RHOKEYA',real(ncount4-ncount0)/real(ncount_rate)
end subroutine rho_segkey

subroutine gridcorrection(nbx,nby,nbz,nl1,nl2,nl3,geocode)
  implicit none
  character(len=1),intent(in) :: geocode
  integer,intent(out) :: nbx,nby,nbz,nl1,nl2,nl3

  !conditions for periodicity in the three directions
  !value of the buffer in the x and z direction
  if (geocode /= 'F') then
     nl1=1
     nl3=1
     nbx = 1
     nbz = 1
  else
     nl1=15
     nl3=15
     nbx = 0
     nbz = 0
  end if
  !value of the buffer in the y direction
  if (geocode == 'P') then
     nl2=1
     nby = 1
  else
     nl2=15
     nby = 0
  end if
end subroutine gridcorrection

subroutine get_boxbound(at,rxyz,radii_cf,crmult,frmult,hxh,hyh,hzh,spadd,dpmult,&
    n1i,n2i,n3i,corx,cory,corz,i1min,i1max,i2min,i2max,i3min,i3max)
  use module_base
  use module_types
  implicit none
  real(gp),intent(in) :: dpmult,spadd
  integer,intent(in) :: n1i,n2i,n3i,corx,cory,corz
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(gp), intent(in) :: crmult,frmult,hxh,hyh,hzh
  real(gp), dimension(at%ntypes,3), intent(in) :: radii_cf
  integer,intent(out) :: i1min,i1max,i2min,i2max,i3min,i3max
  integer,dimension(at%nat,6) :: crbound,frbound
  real(gp) :: sprad,dprad
  integer :: iat
 
  i1min=n1i
  i2min=n2i
  i3min=n3i
  i1max=1
  i2max=1
  i3max=1
  ! setup up the cubes that determines the single and double precision segments
  do iat=1,at%nat
    sprad=radii_cf(at%iatype(iat),1)*crmult+hxh*spadd
    dprad=dpmult*frmult*radii_cf(at%iatype(iat),2)
    crbound(iat,1)=floor((rxyz(1,iat)-sprad)/hxh)+corx
    crbound(iat,2)=floor((rxyz(1,iat)+sprad)/hxh)+corx
    crbound(iat,3)=floor((rxyz(2,iat)-sprad)/hyh)+cory
    crbound(iat,4)=floor((rxyz(2,iat)+sprad)/hyh)+cory
    crbound(iat,5)=floor((rxyz(3,iat)-sprad)/hzh)+corz
    crbound(iat,6)=floor((rxyz(3,iat)+sprad)/hzh)+corz
    frbound(iat,1)=floor((rxyz(1,iat)-dprad)/hxh)+corx
    frbound(iat,2)=floor((rxyz(1,iat)+dprad)/hxh)+corx
    frbound(iat,3)=floor((rxyz(2,iat)-dprad)/hyh)+cory
    frbound(iat,4)=floor((rxyz(2,iat)+dprad)/hyh)+cory
    frbound(iat,5)=floor((rxyz(3,iat)-dprad)/hzh)+corz
    frbound(iat,6)=floor((rxyz(3,iat)+dprad)/hzh)+corz

    crbound(iat,1)=max(crbound(iat,1),1)
    crbound(iat,3)=max(crbound(iat,3),1)
    crbound(iat,5)=max(crbound(iat,5),1)
    frbound(iat,1)=max(frbound(iat,1),1)
    frbound(iat,3)=max(frbound(iat,3),1)
    frbound(iat,5)=max(frbound(iat,5),1)

    crbound(iat,2)=min(crbound(iat,2),n1i)
    crbound(iat,4)=min(crbound(iat,4),n2i)
    crbound(iat,6)=min(crbound(iat,6),n3i)
    frbound(iat,2)=min(frbound(iat,2),n1i)
    frbound(iat,4)=min(frbound(iat,4),n2i)
    frbound(iat,6)=min(frbound(iat,6),n3i)

    i1min=min(i1min,crbound(iat,1))
    i2min=min(i2min,crbound(iat,3))
    i3min=min(i3min,crbound(iat,5))
    i1max=max(i1max,crbound(iat,2))
    i2max=max(i2max,crbound(iat,4))
    i3max=max(i3max,crbound(iat,6))
  enddo
end subroutine get_boxbound

subroutine get_atbound(iat,at,rxyz,radii_cf,crmult,frmult,hxh,hyh,hzh,&
    spadd,dpmult,n1i,n2i,n3i,corx,cory,corz,i1cmin,&
    i1cmax,i2cmin,i2cmax,i3cmin,i3cmax,i1fmin,i1fmax,i2fmin,&
    i2fmax,i3fmin,i3fmax)
  use module_base
  use module_types
  implicit none
  real(gp),intent(in) :: dpmult,spadd
  integer,intent(in) :: n1i,n2i,n3i,iat,corx,cory,corz
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(gp), intent(in) :: crmult,frmult,hxh,hyh,hzh
  real(gp), dimension(at%ntypes,3), intent(in) :: radii_cf
  integer,intent(out) :: i1cmin,i1cmax,i2cmin,i2cmax,i3cmin,i3cmax
  integer,intent(out) :: i1fmin,i1fmax,i2fmin,i2fmax,i3fmin,i3fmax
  real(gp) :: sprad,dprad

  sprad=radii_cf(at%iatype(iat),1)*crmult+hxh*spadd
  dprad=dpmult*frmult*radii_cf(at%iatype(iat),2)

  i1cmin=floor((rxyz(1,iat)-sprad)/hxh)+corx
  i1cmax=floor((rxyz(1,iat)+sprad)/hxh)+corx
  i2cmin=floor((rxyz(2,iat)-sprad)/hyh)+cory
  i2cmax=floor((rxyz(2,iat)+sprad)/hyh)+cory
  i3cmin=floor((rxyz(3,iat)-sprad)/hzh)+corz
  i3cmax=floor((rxyz(3,iat)+sprad)/hzh)+corz

  i1fmin=floor((rxyz(1,iat)-dprad)/hxh)+corx
  i1fmax=floor((rxyz(1,iat)+dprad)/hxh)+corx
  i2fmin=floor((rxyz(2,iat)-dprad)/hyh)+cory
  i2fmax=floor((rxyz(2,iat)+dprad)/hyh)+cory
  i3fmin=floor((rxyz(3,iat)-dprad)/hzh)+corz
  i3fmax=floor((rxyz(3,iat)+dprad)/hzh)+corz

  i1cmin=max(i1cmin,1)
  i2cmin=max(i2cmin,1)
  i3cmin=max(i3cmin,1)
  i1fmin=max(i1fmin,1)
  i2fmin=max(i2fmin,1)
  i3fmin=max(i3fmin,1)

  i1cmax=min(i1cmax,n1i)
  i2cmax=min(i2cmax,n2i)
  i3cmax=min(i3cmax,n3i)
  i1fmax=min(i1fmax,n1i)
  i2fmax=min(i2fmax,n2i)
  i3fmax=min(i3fmax,n3i)

end subroutine get_atbound

subroutine is_overlap(a,b,x,y,overlap)
  integer,intent(in) :: a,b,x,y
  logical,intent(out) :: overlap

  if (((a.le.x).and.(x.le.b)).or.&
      ((a.le.y).and.(y.le.b))) then
    overlap=.true.
  else
    overlap=.false.
  endif
end subroutine is_overlap
