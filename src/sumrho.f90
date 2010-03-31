!!****f* BigDFT/sumrho
!! FUNCTION
!!    Calculate the electronic density (rho)
!! COPYRIGHT
!!    Copyright (C) 2007-2009 CEA, UNIBAS
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!! SOURCE
!!
subroutine sumrho(iproc,nproc,orbs,lr,ixc,hxh,hyh,hzh,psi,rho,nrho,&
     & nscatterarr,nspin,GPU,symObj,irrzon,phnons)
  ! Calculates the charge density by summing the square of all orbitals
  ! Input: psi
  ! Output: rho
  use module_base!, only: gp,dp,wp,ndebug,memocc
  use module_types
  use libxc_functionals

  implicit none
  !Arguments
  integer, intent(in) :: iproc,nproc,nrho,nspin,ixc,symObj
  real(gp), intent(in) :: hxh,hyh,hzh
  type(orbitals_data), intent(in) :: orbs
  type(locreg_descriptors), intent(in) :: lr 
  integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%norbp*orbs%nspinor), intent(in) :: psi
  real(dp), dimension(max(nrho,1),nspin), intent(out), target :: rho
  type(GPU_pointers), intent(inout) :: GPU
  integer, dimension(*), intent(in) :: irrzon
  real(dp), dimension(*), intent(in) :: phnons
  !Local variables
  character(len=*), parameter :: subname='sumrho'
  logical :: rsflag
  integer :: nw1,nw2,nrhotot,n3d,nxc,nxf,itmred
  integer :: ind1,ind2,ind3,ind1s,ind2s,ind3s,oidx,sidx,nspinn
  integer :: i00,i0,i1,i2,i3,i3off,i3s,isjmp,i,j,ispin,iorb,jproc,i_all,i_stat,ierr,j3,j3p
  real(dp) :: charge,tt,maxdiff
  real(dp), dimension(:,:), allocatable :: tmred
  real(dp), dimension(:,:), pointer :: rho_p
  real(dp), dimension(:,:), allocatable :: rho_p_OCL
  real(dp), dimension(:,:), allocatable :: psi_OCL
  integer, dimension(3) :: periodic

!  real(kind=8) :: stream_ptr

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

  !flag for toggling the REDUCE_SCATTER stategy
  rsflag=.not. ((ixc >= 11 .and. ixc <= 16) .or. &
       & (ixc < 0 .and. libxc_functionals_isgga())) .and. .not. have_mpi2

  !calculate dimensions of the complete array to be allocated before the reduction procedure
  if (rsflag) then
     nrhotot=0
     do jproc=0,nproc-1
        nrhotot=nrhotot+nscatterarr(jproc,1)
     end do
  else
     nrhotot=lr%d%n3i
  end if

  if (nproc > 1) then
     allocate(rho_p(lr%d%n1i*lr%d%n2i*nrhotot,nspinn+ndebug),stat=i_stat)
     call memocc(i_stat,rho_p,'rho_p',subname)
  else
     rho_p => rho
  end if

  OCLconv=.true.
  call ocl_create_gpu_context(GPU%context)
  call ocl_create_command_queue(GPU%queue,GPU%context)
  call ocl_build_kernels(GPU%context)
  call init_event_list
  if (lr%geocode /= 'F') then
    periodic(1) = 1
  else
    periodic(1) = 0
  endif
  if (lr%geocode == 'P') then
    periodic(2) = 1
  else
    periodic(2) = 0
  endif 
  if (lr%geocode /= 'F') then
    periodic(3) = 1
  else
    periodic(3) = 0
  endif
  call allocate_data_OCL(lr%d%n1,lr%d%n2,lr%d%n3,periodic,orbs%nspinor,hxh*2.0,hyh*2.0,hzh*2.0,lr%wfd,orbs,GPU)

  if (OCLconv) then
     allocate(rho_p_OCL(max(nrho,1),nspin),stat=i_stat)
     call memocc(i_stat,rho_p_OCL,'rho_p_OCL',subname)
     allocate(psi_OCL((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),orbs%nspinor*orbs%norbp),stat=i_stat)
     call memocc(i_stat,psi_OCL,'psi_OCL',subname)
     call dcopy((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp, psi, 1, psi_OCL, 1)
  end if

  !switch between GPU/CPU treatment of the density
  if (GPUconv) then
     call local_partial_density_GPU(iproc,nproc,orbs,nrhotot,lr,hxh,hyh,hzh,nspin,psi,rho_p,GPU)
  else
     !initialize the rho array at 10^-20 instead of zero, due to the invcb ABINIT routine
     !otherwise use libXC routine
     if (libxc_functionals_isgga()) then
        call razero(lr%d%n1i*lr%d%n2i*nrhotot*nspinn,rho_p)
     else
        call tenminustwenty(lr%d%n1i*lr%d%n2i*nrhotot*nspinn,rho_p,nproc)
     end if

     !for each of the orbitals treated by the processor build the partial densities
     call local_partial_density(iproc,nproc,rsflag,nscatterarr,&
          nrhotot,lr,hxh,hyh,hzh,nspin,orbs,psi,rho_p)
  end if

  if (OCLconv) then
     call local_partial_density_OCL(iproc,nproc,orbs,nrhotot,lr,hxh,hyh,hzh,nspin,psi_OCL,rho_p_OCL,GPU)
     maxdiff=0.0_wp
     do i=1,max(nrho,1)
       do j=1,nspin
        maxdiff=max(maxdiff,abs(rho_p(i,j)-rho_p_OCL(i,j)))
       end do
     end do
     print *,''
     print *,'maxdiff',maxdiff
     i_all=-product(shape(rho_p_OCL))*kind(rho_p_OCL)
     deallocate(rho_p_OCL,stat=i_stat)
     call memocc(i_stat,i_all,'rho_p_OCL',subname)
     i_all=-product(shape(psi_OCL))*kind(psi_OCL)
     deallocate(psi_OCL,stat=i_stat)
     call memocc(i_stat,i_all,'psi_OCL',subname)
     call free_gpu_OCL(GPU,orbs%norbp)
     call ocl_clean(GPU%queue,GPU%context)
     OCLconv=.false.
  end if

  ! Symmetrise density, TODO...
  !after validation this point can be deplaced after the allreduce such as to reduce the number of operations
  if (symObj >= 0 .and. .false.) then
     call symmetrise_density(0,1,lr%d%n1i,lr%d%n2i,lr%d%n3i,nscatterarr,nspin,lr%d%n1i*lr%d%n2i*lr%d%n3i,&
          rho_p,symObj,irrzon,phnons)
  end if

  !the density must be communicated to meet the shape of the poisson solver
  if (nproc > 1) then
     call timing(iproc,'Rho_comput    ','OF')
     call timing(iproc,'Rho_commun    ','ON')
     if (rsflag) then
        do ispin=1,nspin
          call MPI_REDUCE_SCATTER(rho_p(1,ispin),rho(1,ispin),&
               lr%d%n1i*lr%d%n2i*nscatterarr(:,1),&
               MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        end do
     else
         call MPI_ALLREDUCE(MPI_IN_PLACE,rho_p,lr%d%n1i*lr%d%n2i*lr%d%n3i*nspin,&
              MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
         !stop 'rsflag active in sumrho.f90, check MPI2 implementation'
     end if
     call timing(iproc,'Rho_commun    ','OF')
     call timing(iproc,'Rho_comput    ','ON')
     if (.not. rsflag) then
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

end subroutine sumrho
!!***


!!****f* BigDFT/local_partial_density
!! FUNCTION
!!   Here starts the routine for building partial density inside the localisation region
!!   This routine should be treated as a building-block for the linear scaling code
!! SOURCE
!!
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

end subroutine local_partial_density
!!***


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
  integer :: i3s,jproc,i3off,n3d,isjmp,i1,i2,i3,i1s,i1e,j3,i3sg,ithread,nthread
  real(gp) :: hfac2
  real(dp) :: psisq,p1,p2,p3,p4,r1,r2,r3,r4
  integer :: omp_get_thread_num,omp_get_num_threads
  !sum different slices by taking into account the overlap
  i3sg=0
!$omp parallel default(private) shared(n1i,nproc,rsflag,nspinn,nscatterarr,spinsgn) &
!$omp shared(n2i,npsir,hfac,psir,rho_p,n3i,i3sg)
  i3s=0
  hfac2=2.0_gp*hfac
!$  ithread=omp_get_thread_num()
!$  nthread=omp_get_num_threads()

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
!$  if(mod(i3s,nthread) .eq. ithread) then

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
!$  end if

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

end subroutine partial_density



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
  integer :: i3s,jproc,i3off,n3d,isjmp,i1,i2,i3,i1s,i1e,j3,i3sg,ithread,nthread
  real(gp) :: hfac2
  real(dp) :: psisq,p1,p2,p3,p4,r1,r2,r3,r4
  integer :: omp_get_thread_num,omp_get_num_threads
  !sum different slices by taking into account the overlap
  i3sg=0
!$omp parallel default(private) shared(n1i,nproc,rsflag,nspinn,nscatterarr,spinsgn) &
!$omp shared(n2i,npsir,hfac,psir,rho_p,n3i,i3sg,ibyyzz_r)
  i3s=0
!$   ithread=omp_get_thread_num()
!$   nthread=omp_get_num_threads()
  hfac2=2.0_gp*hfac

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
!$    if(mod(i3s,nthread) .eq. ithread) then

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
!$    end if

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

end subroutine partial_density_free

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

!!$     work(:)=rho(:,ispden)
!!$     call fourdp(cplex,rhog,work,-1,mpi_enreg,nfft,ngfft,paral_kgb,0)

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
!!$     if(mpi_enreg%mode_para=='b')then
!!$        old_paral_level=mpi_enreg%paral_level
!!$        mpi_enreg%paral_level=3
!!$        spaceComm=mpi_enreg%comm_fft
!!$        call xsum_mpi(rhosu1_arr,spaceComm,ier)
!!$        call xsum_mpi(rhosu2_arr,spaceComm,ier)
!!$        mpi_enreg%paral_level=old_paral_level
!!$     end if

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
!!$           if(modulo(j/n1i,n2i)/nd2==iproc) then ! this ind is to be treated by me_fft
!!$              !          ind in the proc ind-1=n1(nd2 j3+ r2)+j1
!!$              ind=n1i*(nd2*j3+r2)+j1+1 !this is ind in the current proc
!!$              rhog(1,j1+1,r2+1,j3+1,inzee)=rhosu12(1,izone)*phnons(1,iup+numpt,imagn)&
!!$                   +rhosu12(2,izone)*phnons(2,iup+numpt,imagn)
!!$              rhog(2,j1+1,r2+1,j3+1,inzee)=rhosu12(2,izone)*phnons(1,iup+numpt,imagn)&
!!$                   -rhosu12(1,izone)*phnons(2,iup+numpt,imagn)
!!$           end if
        end do

        !      Keep index of how many points have been considered:
        numpt=numpt+nup

        !      End loop over izone
     end do

     i_all=-product(shape(rhosu12))*kind(rhosu12)
     deallocate(rhosu12,stat=i_stat)
     call memocc(i_stat,i_all,'rhosu12',subname)

     !    Pull out full or spin up density, now symmetrized
!!$     call fourdp(cplex,rhog,work,1,mpi_enreg,nfft,ngfft,paral_kgb,0)

     isign=1
     call fft(n1i,n2i,n3i,n1i+1,n2i+1,n3i+1,rhog,isign,inzee)

     do i3=0,n3i-1
        do i2=0,n2i-1
           do i1=0,n1i-1
              rho(i1+1,i2+1,i3+1,ispden)=rhog(1,i1+1,i2+1,i3+1,inzee)/real(n1i*n2i*n3i,dp)
           enddo
        enddo
     enddo
     !divide by the number of grid points
     !rho(:,ispden)=work(:)

  end do ! ispden

  i_all=-product(shape(rhog))*kind(rhog)
  deallocate(rhog,stat=i_stat)
  call memocc(i_stat,i_all,'rhog',subname)

end subroutine symmetrise_density
