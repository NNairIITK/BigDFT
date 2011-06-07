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
subroutine sumrhoLinear(iproc,nproc,nlr,orbs,Glr,Llr,ixc,hxh,hyh,hzh,psi,rho,nrho,&
     & nscatterarr,nspin,GPU,symObj,irrzon,phnons)
  use module_base!, only: gp,dp,wp,ndebug,memocc
  use module_types
  use libxc_functionals

  implicit none
  !Arguments
  integer, intent(in) :: iproc,nproc,nlr,nrho,nspin,ixc,symObj
  real(gp), intent(in) :: hxh,hyh,hzh
  type(orbitals_data), intent(in) :: orbs
  type(locreg_descriptors),dimension(nlr),intent(in):: LLr
  type(locreg_descriptors), intent(in) :: Glr 
  integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  real(wp), dimension(orbs%npsidim), intent(in) :: psi
  real(dp), dimension(max(nrho,1),nspin), intent(out), target :: rho
  type(GPU_pointers), intent(inout) :: GPU
  integer, dimension(*), intent(in) :: irrzon
  real(dp), dimension(*), intent(in) :: phnons
  !Local variables
  character(len=*), parameter :: subname='sumrho'
  logical :: rsflag
  integer :: nrhotot,n3d,itmred
  integer :: nspinn
  integer :: i1,i2,i3,i3off,i3s,i,ispin,jproc,i_all,i_stat,ierr,j3,j3p,j
  real(dp) :: charge,tt
  real(dp), dimension(:,:), allocatable :: tmred
  real(dp), dimension(:,:), pointer :: rho_p
!!  real(dp), dimension(:,:), allocatable :: rho_p_OCL
!!  real(dp), dimension(:,:), allocatable :: psi_OCL

!  integer :: ncount0,ncount1,ncount2,ncount3,ncountmpi0,ncountmpi1,ncount_max,ncount_rate
!  real(kind=8) :: stream_ptr

  !call system_clock(ncount0,ncount_rate,ncount_max)

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
       & (ixc < 0 .and. libxc_functionals_isgga()))
  
!  write(*,*) 'RSFLAG stuffs ',(ixc >= 11 .and. ixc <= 16),&
!             (ixc < 0 .and. libxc_functionals_isgga()), have_mpi2,rsflag

  !calculate dimensions of the complete array to be allocated before the reduction procedure
  !if (rsflag) then
  !   nrhotot=0
  !   do jproc=0,nproc-1
  !      nrhotot=nrhotot+nscatterarr(jproc,1)
  !   end do
  !else
  !   nrhotot=Glr%d%n3i
  !end if

  !if (nproc > 1) then
  !   allocate(rho_p(Glr%d%n1i*Glr%d%n2i*nrhotot,nspinn+ndebug),stat=i_stat)
  !   call memocc(i_stat,rho_p,'rho_p',subname)
  !else
  !   rho_p => rho
  !end if

  !call system_clock(ncount1,ncount_rate,ncount_max)
!!$  if (OCLconv) then
!!$     allocate(rho_p_OCL(max(nrho,1),nspin),stat=i_stat)
!!$     call memocc(i_stat,rho_p_OCL,'rho_p_OCL',subname)
!!$     allocate(psi_OCL((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f),orbs%nspinor*orbs%norbp),stat=i_stat)
!!$     call memocc(i_stat,psi_OCL,'psi_OCL',subname)
!!$     call dcopy((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp, psi, 1, psi_OCL, 1)
!!$  end if



  !switch between GPU/CPU treatment of the density
  if (GPUconv) then
     stop 'local_partial_density_GPU not implemented!'
     call local_partial_density_GPU(iproc,nproc,orbs,nrhotot,Glr,hxh,hyh,hzh,nspin,psi,rho_p,GPU)
  else if (OCLconv) then
     stop 'local_partial_density_OCL not implemented'
     call local_partial_density_OCL(iproc,nproc,orbs,nrhotot,Glr,hxh,hyh,hzh,nspin,psi,rho_p,GPU)
  else
     !initialize the rho array at 10^-20 instead of zero, due to the invcb ABINIT routine
     !otherwise use libXC routine
     !if (libxc_functionals_isgga()) then
     !   call razero(Glr%d%n1i*Glr%d%n2i*nrhotot*nspinn,rho_p)
     !else
     !   call tenminustwenty(Glr%d%n1i*Glr%d%n2i*nrhotot*nspinn,rho_p,nproc)
     !end if

     !for each of the orbitals treated by the processor build the partial densities
     !call system_clock(ncount2,ncount_rate,ncount_max)
     call local_partial_densityLinear(iproc,nproc,nlr,rsflag,nscatterarr,&
          nrhotot,Glr,Llr,nrho,rho,hxh,hyh,hzh,nspin,orbs,psi)
     !call system_clock(ncount3,ncount_rate,ncount_max)
  end if

!!$  if (OCLconv) then
!!$     call local_partial_density_OCL(iproc,nproc,orbs,nrhotot,Glr,hxh,hyh,hzh,nspin,psi_OCL,rho_p_OCL,GPU)
!!$     maxdiff=0.0_wp
!!$     do i=1,max(nrho,1)
!!$       do j=1,nspin
!!$        maxdiff=max(maxdiff,abs(rho_p(i,j)-rho_p_OCL(i,j)))
!!$       end do
!!$     end do
!!$     print *,''
!!$     print *,'maxdiff',maxdiff
!!$     i_all=-product(shape(rho_p_OCL))*kind(rho_p_OCL)
!!$     deallocate(rho_p_OCL,stat=i_stat)
!!$     call memocc(i_stat,i_all,'rho_p_OCL',subname)
!!$     i_all=-product(shape(psi_OCL))*kind(psi_OCL)
!!$     deallocate(psi_OCL,stat=i_stat)
!!$     call memocc(i_stat,i_all,'psi_OCL',subname)
!!$  end if

  ! Symmetrise density, TODO...
  !after validation this point can be deplaced after the allreduce such as to reduce the number of operations
  if (symObj >= 0) then
     call symmetrise_density(0,1,Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,nscatterarr,nspin,&
          Glr%d%n1i*Glr%d%n2i*Glr%d%n3i,&
          rho_p,symObj,irrzon,phnons)
  end if

  !write(*,*) 'iproc,TIMING:SR1',iproc,real(ncount1-ncount0)/real(ncount_rate)
  !the density must be communicated to meet the shape of the poisson solver
  if (nproc > 1) then
     call timing(iproc,'Rho_comput    ','OF')
     call timing(iproc,'Rho_commun    ','ON')
     if (rsflag) then
        do ispin=1,nspin
          call MPI_REDUCE_SCATTER(rho_p(1,ispin),rho(1,ispin),&
               Glr%d%n1i*Glr%d%n2i*nscatterarr(:,1),&
               MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        end do
     else
       !  call system_clock(ncountmpi0,ncount_rate,ncount_max)
       !  call MPI_ALLREDUCE(MPI_IN_PLACE,rho_p,Glr%d%n1i*Glr%d%n2i*Glr%d%n3i*nspin,&
       !       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
         call mpiallred(rho_p(1,1),Glr%d%n1i*Glr%d%n2i*Glr%d%n3i*nspin,&
               MPI_SUM,MPI_COMM_WORLD,ierr)
       !  call system_clock(ncountmpi1,ncount_rate,ncount_max)
       !  write(*,'(A,4(1X,f6.3))') 'TIMING:ARED'
       !             real(ncount1-ncount0)/real(ncount_rate),&
       !             real(ncount2-ncount1)/real(ncount_rate),&
       !             real(ncount3-ncount2)/real(ncount_rate),&
       !             real(ncountmpi1-ncountmpi0)/real(ncount_rate)
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
              j3p=modulo(j3-1,Glr%d%n3i)+1
              do i2=1,Glr%d%n2i
                 do i1=1,Glr%d%n1i
                    i=i1+(i2-1)*Glr%d%n1i+Glr%d%n1i*Glr%d%n2i*(i3-1)
                    j=i1+(i2-1)*Glr%d%n1i+Glr%d%n1i*Glr%d%n2i*(j3p-1)
                    rho(i,ispin)=rho_p(j,ispin)
                 end do
              end do
           end do
        end do
     end if
  end if
  !call system_clock(ncount2,ncount_rate,ncount_max)
  !write(*,*) 'TIMING:SR2',real(ncount2-ncount1)/real(ncount_rate)
  ! Check
  tt=0.d0
  i3off=Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,4)

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
     do i=1,Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2)
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

  !call system_clock(ncount3,ncount_rate,ncount_max)
  !write(*,*) 'TIMING:SR3',real(ncount3-ncount2)/real(ncount_rate)
  !write(*,*) 'TIMING:SR',real(ncount3-ncount0)/real(ncount_rate)
END SUBROUTINE sumrhoLinear


!>   Here starts the routine for building partial density inside the localisation region
!!   This routine should be treated as a building-block for the linear scaling code
subroutine local_partial_densityLinear(iproc,nproc,nlr,rsflag,nscatterarr,&
     nrhotot,Glr,Llr,nrho,rho,hxh,hyh,hzh,nspin,orbs,psi)
  use module_base
  use module_types
  use module_interfaces
  use libxc_functionals
  implicit none
  logical, intent(in) :: rsflag
  integer, intent(in) :: iproc,nproc,nlr,nrho
  integer,intent(inout):: nrhotot
  integer, intent(in) :: nspin
  real(dp),dimension(max(nrho,1),nspin),intent(out):: rho
  real(gp), intent(in) :: hxh,hyh,hzh
  type(orbitals_data), intent(in) :: orbs
  type(locreg_descriptors), intent(in) :: Glr
  type(locreg_descriptors),dimension(nlr),intent(in) :: Llr
  integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  real(wp), dimension(orbs%npsidim), intent(in) :: psi
  !local variables
  character(len=*), parameter :: subname='local_partial_density'
  integer :: iorb,i_stat,i_all,ii, ind, indSmall, indLarge
  integer :: oidx,sidx,nspinn,npsir,ncomplex, i1, i2, i3, ilr, ispin
  real(gp) :: hfac,spinval
  type(workarr_sumrho) :: w
  real(wp), dimension(:,:), allocatable :: psir
  real(dp), dimension(:),allocatable :: rho_p
  real(8):: dnrm2

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

  !initialisation
  !print *,iproc,'there'
  
  ind=1
!  write(*,*) 'starting locRegLoop...'
  
 !initialize rho
  if (libxc_functionals_isgga()) then
      call razero(max(nrho,1)*nspin,rho)
  else
      call tenminustwenty(max(nrho,1)*nspin,rho,nproc)
  end if


  locRegLoop: do ilr=1,nlr
      call initialize_work_arrays_sumrho(Llr(ilr),w)
      allocate(rho_p(Llr(ilr)%d%n1i*Llr(ilr)%d%n2i*Llr(ilr)%d%n3i*nspinn), stat=i_stat)
      call memocc(i_stat,rho_p,'rho_p',subname)
      allocate(psir(Llr(ilr)%d%n1i*Llr(ilr)%d%n2i*Llr(ilr)%d%n3i,npsir+ndebug),stat=i_stat)
      call memocc(i_stat,psir,'psir',subname)
      
      if (Llr(ilr)%geocode == 'F') then
          call razero(Llr(ilr)%d%n1i*Llr(ilr)%d%n2i*Llr(ilr)%d%n3i*npsir,psir)
      end if
      nrhotot=Llr(ilr)%d%n3i
      orbitalsLoop: do iorb=1,orbs%norbp
         if(orbs%inwhichlocreg(iorb) /= ilr) cycle
         if (libxc_functionals_isgga()) then
             call razero(Llr(ilr)%d%n1i*Llr(ilr)%d%n2i*Llr(ilr)%d%n3i*nspinn,rho_p)
         else
             call tenminustwenty(Llr(ilr)%d%n1i*Llr(ilr)%d%n2i*Llr(ilr)%d%n3i*nspinn,rho_p,nproc)
         end if

         !rho_p=1.d-20

         !print *,'norbp',orbs%norbp,orbs%norb,orbs%nkpts,orbs%kwgts,orbs%iokpt,orbs%occup
         hfac=orbs%kwgts(orbs%iokpt(iorb))*(orbs%occup(orbs%isorb+iorb)/(hxh*hyh*hzh))
         spinval=orbs%spinsgn(orbs%isorb+iorb)

         if (hfac /= 0.d0) then

            !sum for complex function case, npsir=1 in that case
            do oidx=0,ncomplex

               do sidx=1,npsir
                  !call daub_to_isf(Glr,w,psi(1,oidx+sidx,iorb),psir(1,sidx))
                  call daub_to_isf(Llr(ilr),w,psi(ind),psir(1,sidx))
                  ind=ind+Llr(ilr)%wfd%nvctr_c+7*Llr(ilr)%wfd%nvctr_f
               end do

               !print *,'iorb,nrm',iorb,&
               !nrm2(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i*npsir,psir(1,1),1)

               select case(Llr(ilr)%geocode)
               case('F')
                  call partial_density_linear(rsflag,nproc,Llr(ilr)%d%n1i,Llr(ilr)%d%n2i,Llr(ilr)%d%n3i,&
                       npsir,nspinn,nrhotot,&
                       hfac,nscatterarr,spinval,psir,rho_p,Llr(ilr)%bounds%ibyyzz_r)

               case('P')

                  call partial_density_linear(rsflag,nproc,Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,&
                       npsir,nspinn,nrhotot,&
                       hfac,nscatterarr,spinval,psir,rho_p,Llr(ilr)%bounds%ibyyzz_r)

               case('S')

                  call partial_density_linear(rsflag,nproc,Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,&
                       npsir,nspinn,nrhotot,&
                       hfac,nscatterarr,spinval,psir,rho_p,Llr(ilr)%bounds%ibyyzz_r)

               end select

               ! Copy rho_p to the correct place in rho
               indSmall=0
               do ispin=1,nspinn
                   do i3=1,Llr(ilr)%d%n3i
                       do i2=1,Llr(ilr)%d%n2i
                           do i1=1,Llr(ilr)%d%n1i
                               ! indSmall is the index in the currect localization region
                               indSmall=indSmall+1
                               ! indLarge is the index in the whole box. 
                               indLarge=(Llr(ilr)%nsi3+i3-1)*Glr%d%n2i*Glr%d%n1i +&
                                   (Llr(ilr)%nsi2+i2-1)*Glr%d%n1i + Llr(ilr)%nsi1+i1
                       !        print *,'indLarge, Llr(ilr)%nsi',indLarge, Llr(ilr)%nsi1,Llr(ilr)%nsi2,Llr(ilr)%nsi3
                               rho(indLarge,ispin)=rho(indLarge,ispin)+rho_p(indSmall)
                           end do
                       end do
                   end do
               end do

            end do
         end if
      end do orbitalsLoop
      i_all=-product(shape(rho_p))*kind(rho_p)
      deallocate(rho_p,stat=i_stat)
      call memocc(i_stat,i_all,'rho_p',subname)
      i_all=-product(shape(psir))*kind(psir)
      deallocate(psir,stat=i_stat)
      call memocc(i_stat,i_all,'psir',subname)
      call deallocate_work_arrays_sumrho(w)
  end do locRegLoop

  !i_all=-product(shape(psir))*kind(psir)
  !deallocate(psir,stat=i_stat)
  !call memocc(i_stat,i_all,'psir',subname)
  

END SUBROUTINE local_partial_densityLinear

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

!  loop_xc_overlap: do jproc=0,nproc-1
     !case for REDUCE_SCATTER approach, not used for GGA since it enlarges the 
     !communication buffer
!     if (rsflag) then
!        i3off=nscatterarr(jproc,3)-nscatterarr(jproc,4)
!        n3d=nscatterarr(jproc,1)
!        if (n3d==0) exit loop_xc_overlap
!     else
        i3off=0
        n3d=n3i
!     end if
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
!     if (.not. rsflag) exit loop_xc_overlap !the whole range is already done
!  end do loop_xc_overlap
!$omp end parallel

  if (i3sg /= nrhotot) then
     write(*,'(1x,a,i0,1x,i0)')'ERROR: problem with rho_p: i3s,nrhotot,',i3sg,nrhotot
     stop
  end if

!  call system_clock(ncount1,ncount_rate,ncount_max)
!  write(*,*) 'TIMING:PDF',real(ncount1-ncount0)/real(ncount_rate)
END SUBROUTINE partial_density_linear


!> @file
!!  Application of the Hamiltonian + orthonormalize constraints
!! @author
!!    Copyright (C) 2007-2011 CEA
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Application of the Hamiltonian
subroutine LinearHamiltonianApplication(input,iproc,nproc,at,Lzd,hx,hy,hz,rxyz,&
     proj,ngatherarr,pot,psi,hpsi,&
     ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU,radii_cf,pkernel,orbsocc,psirocc)
  use module_base
  use module_types
  use libxc_functionals
  implicit none
  integer, intent(in) :: iproc,nproc,nspin
  real(gp), intent(in) :: hx,hy,hz
  type(atoms_data), intent(in) :: at
  type(input_variables), intent(in) :: input
  type(linear_zone_descriptors),intent(inout) :: Lzd
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(wp), dimension(Lzd%Gnlpspd%nprojel), intent(in) :: proj
  real(wp), dimension(Lzd%orbs%npsidim), intent(in) :: psi
  real(wp), dimension(:), pointer :: pot
  real(gp), intent(out) :: ekin_sum,epot_sum,eexctX,eproj_sum
  real(wp), target, dimension(Lzd%orbs%npsidim), intent(out) :: hpsi
  type(GPU_pointers), intent(inout) :: GPU
  real(gp), dimension(at%ntypes,3+ndebug), intent(in) :: radii_cf
  real(dp), dimension(*), optional :: pkernel
  type(orbitals_data), intent(in), optional :: orbsocc
  real(wp), dimension(:), pointer, optional :: psirocc
  !local variables
  real(gp) :: tmp_ekin_sum,tmp_epot_sum,tmp_eproj_sum
  real(gp), dimension(2,Lzd%orbs%norbp) :: ekin
  real(gp), dimension(2,Lzd%orbs%norbp) :: epot
  real(wp), dimension(:), pointer :: hpsi2
  character(len=*), parameter :: subname='HamiltonianApplication'
  logical :: exctX,op2p
  integer :: i_all,i_stat,ierr,iorb,n3p,ispot,istart_c,iat
  integer :: istart_ck,isorb,ieorb,ikpt,ispsi_k,nspinor,ispsi
  integer :: ilr,dimwf,ind,size_Lpot,size_pot
  real(dp),dimension(:),pointer:: Lpot
!OCL  integer, dimension(3) :: periodic
!OCL  real(wp) :: maxdiff
!OCL  real(gp) :: eproj,ek_fake,ep_fake
  real(gp), dimension(3,2) :: wrkallred
!OCL  real(wp), dimension(:), allocatable :: hpsi_OCL


  !check if the potential has been associated
  if (.not. associated(pot)) then
     if (iproc ==0) then
        write(*,*)' ERROR, HamiltonianApplication, potential not associated!'

        stop
     end if
  end if

  !initialise exact exchange energy 
  op2p=(eexctX == -99.0_gp)
  eexctX=0.0_gp

  exctX = libxc_functionals_exctXfac() /= 0.0_gp

  ! Allocate the nonlocal descriptors for the locregs
  allocate(Lzd%Lnlpspd(Lzd%nlr),stat=i_stat)   

  !initialize accumulators
  ekin_sum = 0.0_gp
  epot_sum = 0.0_gp
  eproj_sum= 0.0_gp
  ind = 1
  do ilr= 1, Lzd%nlr
 
     allocate(Lpot(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*nspin), stat=i_stat)
     call memocc(i_stat,Lpot,'Lpot',subname)
 
     ! replace orbse%norbp by Localnorb for HamiltonianApplication
     Lzd%orbs%norbp = Lzd%Llr(ilr)%Localnorb*nspin
 
     !determine the dimension of the potential array (copied from full_local_potential)
     if (exctX) then
        size_pot=Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*nspin + &
         max(max(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*Lzd%orbs%norbp,ngatherarr(0,1)*Lzd%orbs%norb),1) !part which refers to exact exchange
        size_Lpot=Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*nspin + &
           max(max(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*Lzd%orbs%norbp,&
           ngatherarr(0,1)*Lzd%orbs%norb),1) !CHECK THIS...DOES NOT WORK YET
     else
        size_pot=Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*nspin
        size_Lpot = Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*nspin
     end if
 
     ! Cut the potential into locreg pieces
     call global_to_local(Lzd%Glr,Lzd%Llr(ilr),nspin,size_pot,size_Lpot,pot,Lpot)

     ! Set some quantities: ispot=shift for potential, dimwf=dimension of wavefunction
     ispot=Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*nspin+1
     dimwf=(Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*Lzd%Llr(ilr)%Localnorb*&
           Lzd%orbs%nspinor*nspin

     ! EXACT EXCHANGE NOT TESTED: SHOULD CHECK IF EVERYTHING IF FINE
     !fill the rest of the potential with the exact-exchange terms
     if (present(pkernel) .and. exctX) then
        n3p=ngatherarr(iproc,1)/(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i)
        !exact exchange for virtual orbitals (needs psirocc)
   
        !here we have to add the round part
        if (present(psirocc) .and. present(orbsocc)) then
           call exact_exchange_potential_virt(iproc,nproc,at%geocode,nspin,&
                Lzd%Llr(ilr),orbsocc,Lzd%orbs,ngatherarr(0,1),n3p,&
                0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,pkernel,psirocc,psi(ind:ind+dimwf-1),Lpot)
           eexctX = 0._gp
        else
   !!$        call exact_exchange_potential_round(iproc,nproc,at%geocode,nspin,lr,orbs,&
   !!$             0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,pkernel,psi,pot(ispot),eexctX)
   
           !here the condition for the scheme should be chosen
           if (.not. op2p) then
              call exact_exchange_potential(iproc,nproc,at%geocode,nspin,&
                   Lzd%Llr(ilr),Lzd%orbs,ngatherarr(0,1),n3p,&
                   0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,pkernel,psi(ind:ind+dimwf-1),Lpot,eexctX)
           else
              !the psi should be transformed in real space
              call exact_exchange_potential_round(iproc,nproc,at%geocode,nspin,Lzd%Llr(ilr),Lzd%orbs,&
                   0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,pkernel,psi(ind:ind+dimwf-1),Lpot,eexctX)
   
           end if
        end if
     else
        eexctX = 0._gp
        !print *,'iproc,eexctX',iproc,eexctX
     end if

!     call timing(iproc,'ApplyLocPotKin','ON')

     !apply the local hamiltonian for each of the orbitals
     !given to each processor
     !pot=0.d0
     !psi=1.d0
     !switch between GPU/CPU treatment
!     do i=1,(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp
!          call random_number(psi(i))
!     end do

     if(OCLconv .and. ASYNCconv) then
       allocate(hpsi2((Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*Lzd%orbs%nspinor*Lzd%orbs%norbp),stat=i_stat)
       call memocc(i_stat,hpsi2,'hpsi2',subname)
       hpsi(:)=0.0
     else
       hpsi2 => hpsi
     end if
     if (GPUconv) then
        call local_hamiltonian_GPU(iproc,Lzd%orbs,Lzd%Llr(ilr),hx,hy,hz,nspin,Lpot,psi(ind:ind+dimwf-1),&
             hpsi(ind:ind+dimwf-1),tmp_ekin_sum,tmp_epot_sum,GPU)
     else if (OCLconv) then
        call local_hamiltonian_OCL(iproc,Lzd%orbs,Lzd%Llr(ilr),hx,hy,hz,nspin,Lpot,psi(ind:ind+dimwf-1),&
             hpsi2,tmp_ekin_sum,tmp_epot_sum,GPU,ekin,epot)
     else
        call local_hamiltonian(iproc,Lzd%orbs,Lzd%Llr(ilr),hx,hy,hz,nspin,Lpot,psi(ind:ind+dimwf-1),&
             hpsi(ind:ind+dimwf-1),tmp_ekin_sum,tmp_epot_sum)
     end if

     ekin_sum = ekin_sum + tmp_ekin_sum
     epot_sum = epot_sum + tmp_epot_sum

  !test part to check the results wrt OCL convolutions
!!$  if (OCLconv) then
!!$     allocate(hpsi_OCL((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp+ndebug),stat=i_stat)
!!$     call memocc(i_stat,hpsi_OCL,'hpsi_OCL',subname)
!!$     print *,'fulllocam',GPU%full_locham
!!$     call local_hamiltonian_OCL(iproc,orbs,at%geocode,lr,hx,hy,hz,nspin,pot,psi,hpsi,ek_fake,ep_fake,GPU)
!!$     maxdiff=0.0_wp
!!$     do i=1,(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp
!!$        maxdiff=max(maxdiff,abs(hpsi(i)-hpsi_OCL(i)))
!!$     end do
!!$     print *,'maxdiff',maxdiff
!!$     print *,'ekin_diff',abs(ek_fake-ekin_sum)
!!$     print *,'epot_diff',abs(ep_fake-epot_sum)
!!$     i_all=-product(shape(hpsi_OCL))*kind(hpsi_OCL)
!!$     deallocate(hpsi_OCL,stat=i_stat)
!!$     call memocc(i_stat,i_all,'hpsi_OCL',subname)
!!$  end if

!     call timing(iproc,'ApplyLocPotKin','OF')

  !  apply all PSP projectors for all orbitals belonging to iproc
!     call timing(iproc,'ApplyProj     ','ON')

  !here the localisation region should be changed, temporary only for cubic approach
     eproj_sum=0.0_gp

  ! CUBIC STUFF
  !apply the projectors following the strategy (On-the-fly calculation or not)
!!!  if (DistProjApply .and. .not.present(Lzd)) then
!!!     call applyprojectorsonthefly(iproc,orbs,at,lr,&
!!!          rxyz,hx,hy,hz,lr%wfd,nlpspd,proj,psi,hpsi,eproj_sum)
!!!  else if(orbs%norbp > 0 .and. .not.present(Lzd)) then
!!!     !apply the projectors  k-point of the processor
!!!     !starting k-point
!!!     ikpt=orbs%iokpt(1)
!!!     istart_ck=1
!!!     ispsi_k=1
!!!     loop_kpt: do
!!!
!!!        call orbs_in_kpt(ikpt,orbs,isorb,ieorb,nspinor)
!!!
!!!        ! loop over all my orbitals
!!!        ispsi=ispsi_k
!!!        do iorb=isorb,ieorb
!!!           istart_c=istart_ck
!!!           do iat=1,at%nat
!!!              call apply_atproj_iorb(iat,iorb,istart_c,at,orbs,lr%wfd,nlpspd,&
!!!                   proj,psi(ispsi),hpsi(ispsi),eproj_sum)
!!!           end do
!!!           ispsi=ispsi+(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*nspinor
!!!        end do
!!!        istart_ck=istart_c
!!!        if (ieorb == orbs%norbp) exit loop_kpt
!!!        ikpt=ikpt+1
!!!        ispsi_k=ispsi
!!!     end do loop_kpt
!!!
!!!     if (istart_ck-1 /= nlpspd%nprojel) stop 'incorrect once-and-for-all psp application'
!!!     if (ispsi-1 /= (lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp) stop 'incorrect V_nl psi application'
!!!
!!!  END OF CUBIC STUFF

     if(Lzd%orbs%norbp > 0 ) then
        ! allocate projflg
        allocate(Lzd%Llr(ilr)%projflg(at%nat),stat=i_stat)
        call memocc(i_stat,Lzd%Llr(ilr)%projflg,'Lzd%Llr(ilr)%projflg',subname)

        ! Make the local non-linear pseudopotentials descriptors
        call nlpspd_to_locreg(input,iproc,Lzd%Glr,Lzd%Llr(ilr),rxyz,at,Lzd%orbs,&
      &      radii_cf,input%frmult,input%frmult,input%hx,input%hy,input%hz,Lzd%Gnlpspd,Lzd%Lnlpspd(ilr),Lzd%Llr(ilr)%projflg)
   
        call apply_local_projectors(at,hx,hy,hz,Lzd%Llr(ilr),Lzd%Lnlpspd(ilr),proj,Lzd%orbs,&
                 Lzd%Llr(ilr)%projflg,psi(ind:ind+dimwf-1),rxyz,hpsi(ind:ind+dimwf-1),tmp_eproj_sum)

        eproj_sum = eproj_sum + tmp_eproj_sum
     end if
     ind = ind + dimwf

     ! deallocate Lpot
     call free_full_potential(nproc,Lpot,subname)

  end do
! END LINEAR MODIFICATIONS

  ! local potential and kinetic energy for all orbitals belonging to iproc
  if (iproc==0 .and. verbose > 1) then
     write(*,'(1x,a)',advance='no')&
          'Hamiltonian application...'
  end if

  if(OCLconv .and. ASYNCconv) then
    call finish_hamiltonian_OCL(Lzd%orbs,ekin_sum,epot_sum,GPU,ekin,epot)
    call daxpy(size(hpsi), 1.0_wp, hpsi2(1), 1, hpsi(1),1)
    i_all=-product(shape(hpsi2))*kind(hpsi2)
    deallocate(hpsi2,stat=i_stat)
    call memocc(i_stat,i_all,'hpsi2',subname)
  endif

!  call timing(iproc,'ApplyProj     ','OF')

  !energies reduction
  if (nproc > 1) then
     wrkallred(1,2)=ekin_sum
     wrkallred(2,2)=epot_sum
     wrkallred(3,2)=eproj_sum
     call MPI_ALLREDUCE(wrkallred(1,2),wrkallred(1,1),3,&
          mpidtypg,MPI_SUM,MPI_COMM_WORLD,ierr)
     ekin_sum=wrkallred(1,1)
     epot_sum=wrkallred(2,1)
     eproj_sum=wrkallred(3,1)
  endif



  !up to this point, the value of the potential energy is 
  !only taking into account the local potential part
  !whereas it should consider also the value coming from the 
  !exact exchange operator (twice the exact exchange energy)
  if (exctX) epot_sum=epot_sum+2.0_gp*eexctX

END SUBROUTINE LinearHamiltonianApplication


subroutine LinearDiagHam(iproc,etol,Lzd,orbs,nspin,Lhpsi,Lpsi,psit,orbsv)
  use module_base
  use module_types
  use module_interfaces
  implicit none
  integer, intent(in) :: iproc                                       !> current processor
  integer, intent(in) :: nspin                                       !> number of spin in spin polarized calculation 
  real(gp),intent(in) :: etol                                        !> Tolerance on energy for solve_eigensystem
  type(linear_zone_descriptors) :: Lzd                               !> Information about the locregs
  type(orbitals_data), intent(in) :: orbs                            !> description of orbitals after the input guess (less orbitals)
  type(orbitals_data), optional, intent(in) :: orbsv                 !> description of virtual orbitals (for davidson?)
  real(wp),dimension(Lzd%orbs%npsidim),intent(in):: Lhpsi            !> All the local H|Psi>
  real(wp),dimension(Lzd%orbs%npsidim),intent(in):: Lpsi             !> All the local |Psi>
  real(wp),dimension(orbs%npsidim),intent(inout):: psit                !> Eigenvectors
  !Local variables
  integer :: ilr,ilr2,ikptp,ikpt                             !> loop integers
  integer :: i_stat,i_all                                    !> error handling for allocation/deallocation and memocc
  integer :: ndim_hamovr                                     !> dimension of Global Hamiltonian/overlap matrix
  integer :: psishift1,psishift2                             !> shifting index of wavefunctions for locreg(1) and locreg(2)
  integer :: psidim1,psidim2                                 !> dimension of wavefunctions in locreg(1) and locreg(2)
  integer :: firstrow,lastrow                                !> index of first and last row of hamovr calculted locally
  integer :: firstcol,lastcol                                !> index of first and last col of hamovr calculted locally
  integer :: isovrlp                                         !> number of overlap between locreg(1) and locreg(2) [more then one because of periodic BC]
  integer :: dim_Lhamovr                                     !> dimension of the local, to locreg(1), hamiltonian/overlap matrix
  integer :: norbi_max                                       !> Maximum number of orbitals
  integer :: natsceff,ispsi,norbsc,ldim                      !> 
  integer :: nvctrp,norbtot,ispsie,ispsiv                    !>
  integer :: totshift                                        !> Total shift for the wavefunction when passing from local_to_Global
  integer :: Gpsidim                                         !> Global wavefunction dimension
  character(len=*),parameter :: subname='Linear_DiagHam'     ! name of subroutine
  integer, dimension(:,:), allocatable :: norbgrp            !>
  real(wp), dimension(:), pointer :: psivirt                 !> 
  real(wp),dimension(:,:,:),allocatable :: hamovr            !> Global Hamiltonian/overlap matrix
  integer,dimension(5) :: sizes                              !> array with the sizes for the reshaping of hamovr
  real(wp), dimension(:), pointer :: psi                     !> global quantities for now   SHOULD BE REPLACED
  real(wp),dimension(:,:,:,:,:),allocatable :: work1, work2  !> working arrays for hamovr
  real(wp),dimension(:,:,:),allocatable :: Lhamovr           !> Local,to locreg(1), hamiltonian/overlap matrix



  if (iproc == 0 .and. verbose > 1) write(*,'(1x,a)',advance='no')&
        'Overlap Matrix...'

  ! number of orbitals, dimension and allocation of Global hamiltonian/overlap matrix
  norbi_max=max(Lzd%orbs%norbu,Lzd%orbs%norbd)
  ndim_hamovr = norbi_max**2
  allocate(hamovr(nspin*ndim_hamovr,2,Lzd%orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,hamovr,'hamovr',subname)

  ! reshape hamovr for easy assignation
  allocate(work1(norbi_max, norbi_max, nspin,2,Lzd%orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,work1,'work1',subname)
  sizes = (/ norbi_max, norbi_max, nspin, 2, Lzd%orbs%nkpts+ndebug /)
  work1 = reshape(hamovr,sizes)

  psishift1 = 1
  firstrow  = 1
  lastrow   = 0
  ! The loop on ilr gives the row indexes, the loop on ilr2 gives the column indexes
  do ilr = 1, Lzd%nlr
     firstcol = 1
     lastcol  = 0
     psishift2 = 1
     lastrow  = lastrow  + Lzd%Llr(ilr)%Localnorb
     psidim1 = (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*Lzd%Llr(ilr)%Localnorb*orbs%nspinor
     do ilr2 = 1,Lzd%nlr

        psidim2 = (Lzd%Llr(ilr2)%wfd%nvctr_c+7*Lzd%Llr(ilr2)%wfd%nvctr_f)*Lzd%Llr(ilr2)%Localnorb*Lzd%orbs%nspinor

        call get_number_of_overlap_region(ilr,ilr2,Lzd%Glr,isovrlp,Lzd%Llr,Lzd%nlr)!,outofzone)
        
        ! If no overlap, increment the index of Lpsi and overlap matrix then cycle
        if(isovrlp == 0)then
           psishift2 = psishift2 + psidim2*nspin
           firstcol = firstcol + Lzd%Llr(ilr2)%Localnorb
           lastcol  = lastcol  + Lzd%Llr(ilr2)%Localnorb
           cycle
        end if
        ! dimensions and allocation of Local hamiltonian/overlap matrix
        dim_Lhamovr = Lzd%Llr(ilr)%Localnorb * Lzd%Llr(ilr2)%Localnorb
        allocate(Lhamovr(nspin*dim_Lhamovr,2,Lzd%orbs%nkpts+ndebug),stat=i_stat)
        call memocc(i_stat,Lhamovr,'Lhamovr',subname)

     ! In this routine, we begin by calculating the hamiltonian/overlap matrix between two locregs.
        call overlap_matrix_between_locreg(ilr,ilr2,isovrlp,Lzd%nlr,nspin,psidim1,psidim2,psishift1,&
           psishift2,Lzd%orbs%npsidim,Lzd%orbs,Lzd%Glr,Lzd%Llr,Lpsi,Lhpsi,dim_Lhamovr,Lhamovr)   !deleted outofzone,Localnorb

     ! update the shift for second wavefunction
        psishift2 = psishift2 + psidim2*nspin

     ! reshape the hamiltonian/overlap matrix for easy assignations
       allocate(work2(Lzd%Llr(ilr)%Localnorb,Lzd%Llr(ilr2)%Localnorb,nspin,2,Lzd%orbs%nkpts+ndebug),stat=i_stat)
       call memocc(i_stat,work2,'work2',subname)
       sizes = (/ Lzd%Llr(ilr)%Localnorb, Lzd%Llr(ilr2)%Localnorb, nspin, 2, Lzd%orbs%nkpts+ndebug /)
       work2 = reshape(Lhamovr,sizes)

     ! Assign the calculated values inside global matrix (for truly O(N) this should be replaced) 
       lastcol  = lastcol  + Lzd%Llr(ilr2)%Localnorb
       work1(firstrow:lastrow,firstcol:lastcol,:,:,:) = work2(:,:,:,:,:)

     ! deallocate this instance of Lhamovr
        i_all=-product(shape(work2))*kind(work2)
        deallocate(work2,stat=i_stat)
        call memocc(i_stat,i_all,'work2',subname)

        i_all=-product(shape(Lhamovr))*kind(Lhamovr)
        deallocate(Lhamovr,stat=i_stat)
        call memocc(i_stat,i_all,'Lhamovr',subname)

     ! update indexes
       firstcol = firstcol + Lzd%Llr(ilr2)%Localnorb
     end do
     ! increment the shift of wavefunctions
     psishift1 = psishift1 + psidim1*nspin
     firstrow = firstrow + Lzd%Llr(ilr)%Localnorb
  end do

  ! reshape back to original shape
  hamovr = reshape(work1,(/ nspin*ndim_hamovr, 2, Lzd%orbs%nkpts+ndebug /))
  i_all=-product(shape(work1))*kind(work1)
  deallocate(work1,stat=i_stat)
  call memocc(i_stat,i_all,'work1',subname)

! DEBUG
!  print *,'hamovr, ham:',hamovr(:,1,:)
!  print *,'hamovr, ovr:',hamovr(:,2,:)
! END DEBUG

  ! Don't need Lhpsi anymore
!  i_all=-product(shape(Lhpsi))*kind(Lhpsi)
!  deallocate(Lhpsi,stat=i_stat)
!  call memocc(i_stat,i_all,'Lhpsi',subname)

  ! Now solve the eigensystem: H |Lpsi> = epsilon S |Lpsi>  
  if(iproc==0) write(*,'(1x,a)') 'Direct diagonalization...'

!  call timing(iproc, 'Input_comput', 'ON')

  ! SET SOME VARIABLE FOR NOW (NO SEMICORE)
  ispsi=1
  natsceff = 0
  allocate(norbgrp(1,nspin+ndebug),stat=i_stat)
  call memocc(i_stat,norbgrp,'norbgrp',subname)
  norbsc=0
  norbgrp(1,1)=Lzd%orbs%norbu
  if (nspin == 2) norbgrp(1,2)=Lzd%orbs%norbd
  !it is important that the k-points repartition of the inputguess orbitals
  !coincides with the one of the SCF orbitals
  ! NOTE: USE ORBS OR ORBSE??
  do ikptp=1,Lzd%orbs%nkptsp
     ikpt=Lzd%orbs%iskpts+ikptp!orbs%ikptsp(ikptp)
     call solve_eigensystem(iproc,orbs%norb,orbs%norbu,orbs%norbd,norbi_max,&
          ndim_hamovr,natsceff,nspin,Lzd%orbs%nspinor,etol,norbgrp,hamovr(1,1,ikpt),&
          orbs%eval((ikpt-1)*orbs%norb+1))
  end do

  if (iproc == 0 .and. verbose > 1) write(*,'(1x,a)',advance='no')'Building orthogonal Wavefunctions...'


! FOR NOW, just transform the Lpsi to psi in global region.
  Gpsidim = (Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*Lzd%orbs%norb*Lzd%orbs%nspinor
  allocate(psi(Gpsidim+ndebug),stat=i_stat)
  call memocc(i_stat,psi,'psi',subname)
  call razero(Gpsidim+ndebug,psi)

! WATCH OUT, does not work for nspinor > 1
  psishift1 = 1
  totshift = 0
  do ilr = 1,Lzd%nlr
     ldim = (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*Lzd%Llr(ilr)%Localnorb*Lzd%orbs%nspinor*nspin
     call Lpsi_to_global(Lzd%Glr,Gpsidim,Lzd%Llr(ilr),Lpsi(psishift1:psishift1+ldim-1),&
          ldim,Lzd%Llr(ilr)%Localnorb,Lzd%orbs%nspinor,nspin,totshift,psi)
     psishift1 = psishift1 + ldim
     totshift = (Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*Lzd%Llr(ilr)%Localnorb*Lzd%orbs%nspinor
  end do

!  allocate(psit(orbs%npsidim+ndebug),stat=i_stat)
!  call memocc(i_stat,psit,'psit',subname)

  ispsi=1
  ispsie=1
  ispsiv=1
  norbtot = Lzd%orbs%norb
  do ikptp=1,Lzd%orbs%nkptsp
     ikpt=Lzd%orbs%iskpts+ikptp!orbsu%ikptsp(ikptp)\

!    nvctrp is not a simple quantity anymore has it depends on the locregs (can be different for every locreg)
!    for an O(N) code, should change these routines.
!    FOR NOW, just transform the Lpsi to psi in global region.
     nvctrp=Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f
     if (nvctrp == 0) cycle
     call build_eigenvectors(iproc,orbs%norbu,orbs%norbd,orbs%norb,norbtot,nvctrp,&
          natsceff,nspin,Lzd%orbs%nspinor,Lzd%orbs%nspinor,ndim_hamovr,norbgrp,hamovr(1,1,ikpt),&
          psi(ispsie:),psit(ispsi:))
     ispsi=ispsi+nvctrp*orbs%norb*Lzd%orbs%nspinor
     ispsie=ispsie+nvctrp*norbtot*Lzd%orbs%nspinor
     if (present(orbsv)) ispsiv=ispsiv+nvctrp*orbsv%norb*Lzd%orbs%nspinor
  end do


  !if(nproc==1.and.nspinor==4) call psitransspi(nvctrp,norbu+norbd,psit,.false.)
  if (iproc == 0 .and. verbose > 1) write(*,'(1x,a)') 'done.'

  !deallocate psi
  i_all=-product(shape(psi))*kind(psi)
  deallocate(psi,stat=i_stat)
  call memocc(i_stat,i_all,'psi',subname)

end subroutine LinearDiagHam

