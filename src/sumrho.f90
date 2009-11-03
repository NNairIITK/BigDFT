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
subroutine sumrho(iproc,nproc,orbs,lr,ixc,hxh,hyh,hzh,psi,rho,nrho,nscatterarr,nspin,GPU)
  ! Calculates the charge density by summing the square of all orbitals
  ! Input: psi
  ! Output: rho
  use module_base!, only: gp,dp,wp,ndebug,memocc
  use module_types
  use libxc_functionals

  implicit none
  integer, intent(in) :: iproc,nproc,nrho,nspin,ixc
  real(gp), intent(in) :: hxh,hyh,hzh
  type(orbitals_data), intent(in) :: orbs
  type(locreg_descriptors), intent(in) :: lr 
  integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%norbp*orbs%nspinor), intent(in) :: psi
  real(dp), dimension(max(nrho,1),nspin), intent(out), target :: rho
  type(GPU_pointers), intent(inout) :: GPU
  !local variables
  character(len=*), parameter :: subname='sumrho'
  logical :: rsflag
  integer :: nw1,nw2,nrhotot,n3d,n1i,n2i,n3i,nxc,nxf,itmred
  integer :: ind1,ind2,ind3,ind1s,ind2s,ind3s,oidx,sidx,nspinn
  integer :: i00,i0,i1,i2,i3,i3off,i3s,isjmp,i,ispin,iorb,jproc,i_all,i_stat,ierr,j3,j3p,j
  real(dp) :: charge,tt
  real(dp), dimension(:,:), allocatable :: tmred
  real(dp), dimension(:,:), pointer :: rho_p

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

  !switch between GPU/CPU treatment of the density
  if (GPUconv) then
     call local_partial_density_GPU(iproc,nproc,orbs,nrhotot,lr,hxh,hyh,hzh,nspin,psi,rho_p,GPU)
  else
     !initialize the rho array at 10^-20 instead of zero, due to the invcb ABINIT routine
     !call razero(lr%d%n1i*lr%d%n2i*nrhotot*nspinn,rho_p)
     call tenminustwenty(lr%d%n1i*lr%d%n2i*nrhotot*nspinn,rho_p,nproc)

     !for each of the orbitals treated by the processor build the partial densities
     call local_partial_density(iproc,nproc,rsflag,nscatterarr,&
          nrhotot,lr,hxh,hyh,hzh,nspin,orbs,psi,rho_p)
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

