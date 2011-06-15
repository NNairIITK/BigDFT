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
  integer :: iorb,i_stat,i_all,ii, ind, indSmall, indLarge, orbtot
  integer :: oidx,sidx,nspinn,npsir,ncomplex, i1, i2, i3, ilr, ispin
  integer :: nspincomp
  real(gp) :: hfac,spinval
  type(workarr_sumrho) :: w
  integer, allocatable :: inthisLocreg(:)
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
  nspincomp = 1
  if (nspin > 1) nspincomp = 2

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
      allocate(inthisLocreg(Llr(ilr)%localnorb*nspincomp),stat=i_stat)
      call memocc(i_stat,inthisLocreg,'inthisLocreg',subname)
         
      orbtot = 0      
      do iorb=1,orbs%norbp
         if (orbs%inWhichLocreg(iorb) == ilr) then
            orbtot = orbtot+1
            inthisLocreg(orbtot) = iorb
         end if
      end do
      if (orbtot .ne. Llr(ilr)%localnorb*nspincomp) then
         write(*,*) 'Problem in local_hamiltonian_Linear, orbtot=',orbtot,'is not equal to localnorb=',Llr(ilr)%localnorb*nspincomp
         
         stop
      end if

      orbitalsLoop: do ii=1,orbtot

         iorb = inthisLocreg(ii)   !using ii and iorb to identify the orbitals because in linear case, the ordering is different
                                   !orbitals are now orderer by locreg. So, iorb is the old numbering (i.e. in Global region)
                                   !while ii is it's numbering in the locreg.

         if (libxc_functionals_isgga()) then
             call razero(Llr(ilr)%d%n1i*Llr(ilr)%d%n2i*Llr(ilr)%d%n3i*nspinn,rho_p)
         else
             call tenminustwenty(Llr(ilr)%d%n1i*Llr(ilr)%d%n2i*Llr(ilr)%d%n3i*nspinn,rho_p,nproc)
         end if

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
                       if(Llr(ilr)%nsi3 + i3 - 1 < 0) cycle   !throwing away the extra buffer of the locreg, related to the change of boundary conditions
                       if(Llr(ilr)%nsi3 + i3  > Glr%d%n3i) cycle
                       do i2=1,Llr(ilr)%d%n2i
                           if(Llr(ilr)%nsi2+ i2 - 1 < 0) cycle !same
                           if(Llr(ilr)%nsi2 + i2  > Glr%d%n2i) cycle
                           do i1=1,Llr(ilr)%d%n1i
                               if(Llr(ilr)%nsi1+ i1 - 1 < 0) cycle ! same
                               if(Llr(ilr)%nsi1 + i1  > Glr%d%n1i) cycle
                               ! indSmall is the index in the currect localization region
                               indSmall=indSmall+1
                               ! indLarge is the index in the whole box. 
                               indLarge=(Llr(ilr)%nsi3+i3-1)*Glr%d%n2i*Glr%d%n1i +&
                                   (Llr(ilr)%nsi2+i2-1)*Glr%d%n1i + Llr(ilr)%nsi1+i1
                               rho(indLarge,ispin)=rho(indLarge,ispin)+rho_p(indSmall)
                           end do
                       end do
                   end do
               end do

            end do
         else
            ind=ind+(Llr(ilr)%wfd%nvctr_c+7*Llr(ilr)%wfd%nvctr_f)*max(ncomplex,1)*npsir
         end if
      end do orbitalsLoop
      i_all=-product(shape(rho_p))*kind(rho_p)
      deallocate(rho_p,stat=i_stat)
      call memocc(i_stat,i_all,'rho_p',subname)
      i_all=-product(shape(psir))*kind(psir)
      deallocate(psir,stat=i_stat)
      call memocc(i_stat,i_all,'psir',subname)
      call deallocate_work_arrays_sumrho(w)
      i_all=-product(shape(inthisLocreg))*kind(inthisLocreg)
      deallocate(inthisLocreg,stat=i_stat)
      call memocc(i_stat,i_all,'inthisLocreg',subname)
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
  integer :: tmp_norbp,nspincomp
  real(dp),dimension(:),pointer:: Lpot
  real(wp),dimension(:),allocatable :: hpsi_proj
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

  nspincomp = 1
  if (nspin > 1) then
     nspincomp = 2
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
 
 
     !determine the dimension of the potential array (copied from full_local_potential)
     if (exctX) then
        size_pot=Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*nspin + &
         max(max(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*Lzd%Llr(ilr)%Localnorb*nspin,ngatherarr(0,1)*Lzd%orbs%norb),1) !part which refers to exact exchange
        size_Lpot=Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*nspin + &
           max(max(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*Lzd%Llr(ilr)%Localnorb*nspin,&
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
           Lzd%orbs%nspinor*nspincomp

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
       allocate(hpsi2((Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*Lzd%orbs%nspinor*&
                Lzd%Llr(ilr)%Localnorb*nspin),stat=i_stat)
       call memocc(i_stat,hpsi2,'hpsi2',subname)
       hpsi(:)=0.0
     else
       hpsi2 => hpsi
     end if
     if (GPUconv) then  !does not work yet
        call local_hamiltonian_GPU(iproc,Lzd%orbs,Lzd%Llr(ilr),hx,hy,hz,nspin,Lpot,psi(ind:ind+dimwf-1),&
             hpsi(ind:ind+dimwf-1),tmp_ekin_sum,tmp_epot_sum,GPU,ilr)
     else if (OCLconv) then  ! does_not_work yet
        call local_hamiltonian_OCL(iproc,Lzd%orbs,Lzd%Llr(ilr),hx,hy,hz,nspin,Lpot,psi(ind:ind+dimwf-1),&
             hpsi2,tmp_ekin_sum,tmp_epot_sum,GPU,ekin,epot,ilr)
     else
        call local_hamiltonian_Linear(iproc,ilr,Lzd%orbs,Lzd%Llr(ilr),hx,hy,hz,nspin,Lpot,psi(ind:ind+dimwf-1),&
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
  !   eproj_sum=0.0_gp

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

     if(Lzd%orbs%norbp > 0) then
        !allocate
        if(ilr == 1) then
           allocate(hpsi_proj(Lzd%orbs%npsidim),stat=i_stat)
           call memocc(i_stat,hpsi_proj,'hpsi_proj',subname)
           hpsi_proj = 0.0_wp
        end if

        ! allocate projflg
        allocate(Lzd%Llr(ilr)%projflg(at%nat),stat=i_stat)
        call memocc(i_stat,Lzd%Llr(ilr)%projflg,'Lzd%Llr(ilr)%projflg',subname)

        ! Make the local non-linear pseudopotentials descriptors
        call nlpspd_to_locreg(input,iproc,Lzd%Glr,Lzd%Llr(ilr),rxyz,at,Lzd%orbs,&
      &      radii_cf,input%frmult,input%frmult,hx,hy,hz,Lzd%Gnlpspd,Lzd%Lnlpspd(ilr),Lzd%Llr(ilr)%projflg)

        call apply_local_projectors(ilr,nspin,at,hx,hy,hz,Lzd%Llr(ilr),Lzd%Lnlpspd(ilr),proj,Lzd%orbs,&
                 Lzd%Llr(ilr)%projflg,psi(ind:ind+dimwf-1),rxyz,hpsi(ind:ind+dimwf-1),eproj_sum)
        ! accumulate the new hpsi
        hpsi_proj(ind:ind+dimwf-1) = hpsi_proj(ind:ind+dimwf-1) + hpsi(ind:ind+dimwf-1)
     end if
     ind = ind + dimwf

     ! deallocate Lpot
     call free_full_potential(nproc,Lpot,subname)

  end do
! END LINEAR MODIFICATIONS

  ! Now that all is accumulated, rename hpsi_proj to hpsi
  hpsi = hpsi_proj

  !deallocate hpsi_proj
  i_all=-product(shape(hpsi_proj))*kind(hpsi_proj)
  deallocate(hpsi_proj,stat=i_stat)
  call memocc(i_stat,i_all,'hpsi_proj',subname)

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
     totshift = totshift + (Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*Lzd%Llr(ilr)%Localnorb*Lzd%orbs%nspinor
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
  end do


  !if(nproc==1.and.nspinor==4) call psitransspi(nvctrp,norbu+norbd,psit,.false.)
  if (iproc == 0 .and. verbose > 1) write(*,'(1x,a)') 'done.'

  !deallocate psi
  i_all=-product(shape(psi))*kind(psi)
  deallocate(psi,stat=i_stat)
  call memocc(i_stat,i_all,'psi',subname)

end subroutine LinearDiagHam


!> @file
!!  Routine to calculate the action of the hamiltonian
!! @author
!!   Copyright (C) 2005-2011 BigDFT group 
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 


!>   Calculate the action of the local hamiltonian on the orbitals
subroutine local_hamiltonian_Linear(iproc,ilr,orbs,lr,hx,hy,hz,&
     nspin,pot,psi,hpsi,ekin_sum,epot_sum)
  use module_base
  use module_types
  use module_interfaces
  use libxc_functionals
  implicit none
  integer, intent(in) :: iproc,nspin,ilr
  real(gp), intent(in) :: hx,hy,hz
  type(orbitals_data), intent(in) :: orbs
  type(locreg_descriptors), intent(in) :: lr
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*orbs%norbp), intent(in) :: psi
  real(wp), dimension(*) :: pot
  !real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i*nspin) :: pot
  real(gp), intent(out) :: ekin_sum,epot_sum
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*orbs%norbp), intent(out) :: hpsi
  !local variables
  character(len=*), parameter :: subname='local_hamiltonian_Linear'
  integer :: i_all,i_stat,iorb,npot,nsoffset,oidx,ispot
  integer :: ii,orbtot
  integer :: nspincomp
  integer,dimension(:),allocatable :: inthisLocreg
  real(wp) :: exctXcoeff
  real(gp) :: ekin,epot,kx,ky,kz,etest
  type(workarr_locham) :: wrk_lh
  real(wp), dimension(:,:), allocatable :: psir

  exctXcoeff=libxc_functionals_exctXfac()

  nspincomp = 1
  if(nspin > 1) then
    nspincomp = 2
  end if

  allocate(inthisLocreg(lr%localnorb*nspincomp),stat=i_stat)
  call memocc(i_stat,inthisLocreg,'inthisLocreg',subname)

  !initialise the work arrays
  call initialize_work_arrays_locham(lr,orbs%nspinor,wrk_lh)

  !components of the potential
  npot=orbs%nspinor
  if (orbs%nspinor == 2) npot=1

  ! Wavefunction in real space
  allocate(psir(lr%d%n1i*lr%d%n2i*lr%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
  call memocc(i_stat,psir,'psir',subname)

  call razero(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspinor,psir)

  ekin_sum=0.0_gp
  epot_sum=0.0_gp

  etest=0.0_gp
  
  orbtot = 0
  do iorb=1,orbs%norbp
     if (orbs%inWhichLocreg(iorb) == ilr) then
        orbtot = orbtot+1
        inthisLocreg(orbtot) = iorb
     end if
  end do 
  if (orbtot .ne. lr%localnorb*nspincomp) then
     write(*,*) 'Problem in local_hamiltonian_Linear, orbtot=',orbtot,'is not equal to localnorb=',lr%localnorb*nspincomp
     stop
  end if

  do ii=1,orbtot

     iorb = inthisLocreg(ii)   !using ii and iorb to identify the orbitals because in linear case, the ordering is different
                               !orbitals are now orderer by locreg. So, iorb is the old numbering (i.e. in Global region)
                               !while ii is it's numbering in the locreg.

     if(orbs%spinsgn(iorb+orbs%isorb)>0.0_gp .or. nspin == 1 .or. nspin == 4 ) then
        nsoffset=1
     else
        nsoffset=lr%d%n1i*lr%d%n2i*lr%d%n3i+1
     end if

     oidx=(ii-1)*orbs%nspinor+1

     !transform the wavefunction in Daubechies basis to the wavefunction in ISF basis
     !the psir wavefunction is given in the spinorial form
     call daub_to_isf_locham(orbs%nspinor,lr,wrk_lh,psi(1,oidx),psir)

     !ispot=1+lr%d%n1i*lr%d%n2i*lr%d%n3i*(nspin+iorb-1)
     !etest=etest+dot(lr%d%n1i*lr%d%n2i*lr%d%n3i,pot(ispot),1,psir(1,1),1)
     !print *,'epot, iorb,iproc,norbp',iproc,orbs%norbp,iorb,etest

     !apply the potential to the psir wavefunction and calculate potential energy
     select case(lr%geocode)
     case('F')

        call apply_potential(lr%d%n1,lr%d%n2,lr%d%n3,1,1,1,0,orbs%nspinor,npot,psir,&
             pot(nsoffset),epot,&
             lr%bounds%ibyyzz_r) !optional

     case('P')
        !here the hybrid BC act the same way
        call apply_potential(lr%d%n1,lr%d%n2,lr%d%n3,0,0,0,0,orbs%nspinor,npot,psir,&
             pot(nsoffset),epot)

     case('S')

        call apply_potential(lr%d%n1,lr%d%n2,lr%d%n3,0,1,0,0,orbs%nspinor,npot,psir,&
             pot(nsoffset),epot)
     end select

     !k-point values, if present
     kx=orbs%kpts(1,orbs%iokpt(iorb))
     ky=orbs%kpts(2,orbs%iokpt(iorb))
     kz=orbs%kpts(3,orbs%iokpt(iorb))

     if (exctXcoeff /= 0.0_gp) then
        ispot=1+lr%d%n1i*lr%d%n2i*lr%d%n3i*(nspincomp+ii-1)
        !add to the psir function the part of the potential coming from the exact exchange
        call axpy(lr%d%n1i*lr%d%n2i*lr%d%n3i,exctXcoeff,pot(ispot),1,psir(1,1),1)
     end if

     !apply the kinetic term, sum with the potential and transform back to Daubechies basis
     call isf_to_daub_kinetic(hx,hy,hz,kx,ky,kz,orbs%nspinor,lr,wrk_lh,&
          psir,hpsi(1,oidx),ekin)
!     print *,iorb, ekin+epot, epot
     ekin_sum=ekin_sum+orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*ekin
     epot_sum=epot_sum+orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*epot
  enddo

  !print *,'iproc,etest',etest

  !deallocations of work arrays
  i_all=-product(shape(psir))*kind(psir)
  deallocate(psir,stat=i_stat)
  call memocc(i_stat,i_all,'psir',subname)

  i_all=-product(shape(inthisLocreg))*kind(inthisLocreg)
  deallocate(inthisLocreg,stat=i_stat)
  call memocc(i_stat,i_all,'inthisLocreg',subname)

  call deallocate_work_arrays_locham(lr,wrk_lh)

END SUBROUTINE local_hamiltonian_Linear


subroutine sumrhoForLocalizedBasis(iproc, nproc, orbs, Glr, input, lin, coeff, phi, nrho, rho, at, rxyz, nscatterarr, &
            phibuff)
!
use module_base
use module_types
use libxc_functionals
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, nrho
type(orbitals_data),intent(in):: orbs
type(locreg_descriptors),intent(in):: Glr
type(input_variables),intent(in):: input
type(linearParameters),intent(inout):: lin
real(8),dimension(lin%orbs%norb,orbs%norb),intent(in):: coeff
real(8),dimension(lin%orbs%npsidim),intent(in):: phi
real(8),dimension(nrho),intent(out),target:: rho
type(atoms_data),intent(in):: at
real(8),dimension(3,at%nat),intent(in):: rxyz
integer, dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
real(8),dimension(lin%comsr%sizePhibuff):: phibuff

! Local variables
integer:: iorb, jorb, korb, istat, ind, ind1, ind2, indLarge, i1, i2, i3, ilr, jlr, i_f, iseg_f
integer:: i1s, i1e, i2s, i2e, i3s, i3e, i1d, j1d, i2d, j2d, i3d, j3d, indri, indrj, ldim, gdim, iall, istr, istri, istrj
integer:: indi2, indi3, indj2, indj3, indl2, indl3, mpisource, mpidest, orbitalsource, tag, lrsource, iiorb, jjorb
integer:: ist, klr, cnt, ierr, istrecv, ii, ilrprev, istSend, nSendRecv, nSendRecvCheck, requestStart, sizePhibuffr, kst
integer:: jproc, jprocprev, mpidestprev, is, ie, ioverlap, orbitaldest, sizePhibuff, kproc, nreceives, istsource, istdest, ncount
integer:: nsends, nfast, nslow, nsameproc
real(8):: tt, hxh, hyh, hzh, factor, totalCharge, tr, partialCharge
real(8),dimension(:),allocatable:: phir1, phir2, Lphi, Lphir1, lPhir2, rho2, phibuffr
real(8),dimension(:,:),allocatable:: densKern
real(8),dimension(:),pointer:: rhofull
type(workarr_sumrho):: w
character(len=*),parameter:: subname='sumrhoForLocalizedBasis'
real(8),dimension(0:3):: scal
integer,dimension(:),allocatable:: request, startInd, nSendRecvArr, requestCheck,fromWhichLocreg
!integer,dimension(:,:,:),allocatable:: commsSumrho
logical,dimension(:,:),allocatable:: communComplete, computComplete
integer,dimension(mpi_status_size):: stat
integer,dimension(mpi_status_size,2):: stats
logical:: quit, sendComplete, receiveComplete



!!! Determine which orbital are needed for each slice
!!if(iproc==0) then
!!    do jproc=0,nproc-1
!!        write(*,'(a,i5,4i9)') 'jproc, n3d, n3p, i3s+i3xcsh-1, i3xcsh', jproc, nscatterarr(jproc,1), nscatterarr(jproc,2), nscatterarr(jproc,3), nscatterarr(jproc,4)
!!    end do
!!end if

lin%comsr%communComplete=.false.
lin%comsr%computComplete=.false.


allocate(phibuffr(lin%comsr%sizePhibuffr), stat=istat)
call memocc(istat, phibuffr, 'phibuffr', subname)
phibuffr=0.d0




! Allocate the density kernel.
allocate(densKern(lin%orbs%norb,lin%orbs%norb), stat=istat)
call memocc(istat, densKern, 'densKern', subname)

! Calculate the density kernel.
do iorb=1,lin%orbs%norb
    do jorb=1,lin%orbs%norb
        tt=0.d0
        do korb=1,orbs%norb
            tt = tt + coeff(iorb,korb)*coeff(jorb,korb)
        end do
        densKern(iorb,jorb)=tt
        !write(*,'(a,2i6,es15.7)') 'iorb, jorb, densKern(iorb,jorb)', iorb, jorb, densKern(iorb,jorb)
    end do
end do


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





! Wait for the communications that have not completed yet
nfast=0
nsameproc=0
testLoop: do
    !do korb=1,nreceives
    do korb=1,lin%comsr%noverlaps(iproc)
        if(lin%comsr%communComplete(korb,iproc)) cycle
!        call mpi_test(lin%comsr%comarr(8,korb,iproc), sendComplete, stat, ierr)      !COMMENTED BY PB
!        call mpi_test(lin%comsr%comarr(9,korb,iproc), receiveComplete, stat, ierr)   !COMMENTED BY PB
        if(sendComplete .and. receiveComplete) lin%comsr%communComplete(korb,iproc)=.true.
        if(lin%comsr%communComplete(korb,iproc)) then
            !write(*,'(2(a,i0))') 'fast communication; process ', iproc, ' has received orbital ', korb
            mpisource=lin%comsr%comarr(1,korb,iproc)
            mpidest=lin%comsr%comarr(5,korb,iproc)
            if(mpisource/=mpidest) then
                nfast=nfast+1
            else
                nsameproc=nsameproc+1
            end if
            kst=lin%comsr%comarr(6,korb,iproc)
            klr=lin%comsr%comarr(4,korb,iproc)
            istr=lin%comsr%istrarr(korb)
            !write(*,'(a,i0,a,2i10)') 'process ', iproc, ' calls daub_to_isf with ', kst, istr
            call initialize_work_arrays_sumrho(lin%Llr(klr), w)
            call daub_to_isf(lin%Llr(klr), w, phibuff(kst), phibuffr(istr))
            call deallocate_work_arrays_sumrho(w)
            lin%comsr%computComplete(korb,iproc)=.true.
        end if
    end do
    ! If we made it until here, either all all the communication is
    ! complete or we better wait for each single orbital.
    exit testLoop
end do testLoop



! Wait for the communications that have not completed yet
!do korb=1,nreceives
nslow=0
do korb=1,lin%comsr%noverlaps(iproc)
    if(lin%comsr%communComplete(korb,iproc)) cycle
    !write(*,'(2(a,i0))') 'process ', iproc, ' is waiting for orbital ', korb
    nslow=nslow+1
!    call mpi_wait(lin%comsr%comarr(8,korb,iproc), stat, ierr)   !COMMENTED BY PB
!    call mpi_wait(lin%comsr%comarr(9,korb,iproc), stat, ierr)   !COMMENTED BY PB
    lin%comsr%communComplete(korb,iproc)=.true.
    !write(*,'(2(a,i0))') 'process ', iproc, ' has received orbital ', korb
    kst=lin%comsr%comarr(6,korb,iproc)
    klr=lin%comsr%comarr(4,korb,iproc)
    istr=lin%comsr%istrarr(korb)
    !write(*,'(a,i0,a,2i10,3(a,i0))') 'process ', iproc, ' call daub_to_isf with ', kst, istr, ', size(phibuffr)=',size(phibuffr), ' klr=', klr, ' size(phibuff)=', size(phibuff)
    call initialize_work_arrays_sumrho(lin%Llr(klr), w)
    call daub_to_isf(lin%Llr(klr), w, phibuff(kst), phibuffr(istr))
    call deallocate_work_arrays_sumrho(w)
    lin%comsr%computComplete(korb,iproc)=.true.
end do

call mpiallred(nreceives, 1, mpi_sum, mpi_comm_world, ierr)
call mpiallred(nfast, 1, mpi_sum, mpi_comm_world, ierr)
call mpiallred(nslow, 1, mpi_sum, mpi_comm_world, ierr)
call mpiallred(nsameproc, 1, mpi_sum, mpi_comm_world, ierr)
if(iproc==0) write(*,'(x,2(a,i0),a)') 'statistics: - ', nfast+nslow, ' point to point communications, of which ', &
                       nfast, ' could be overlapped with computation.'
if(iproc==0) write(*,'(x,a,i0,a)') '            - ', nsameproc, ' copies on the same processor.'


quit=.false.
do iorb=1,lin%comsr%noverlaps(iproc)
    if(.not. lin%comsr%communComplete(iorb,iproc)) then
        write(*,'(a,i0,a,i0,a)') 'ERROR: communication of orbital ', iorb, ' to process ', iproc, ' failed!'
        !quit=.true.
        stop
    end if
    if(.not. lin%comsr%computComplete(iorb,iproc)) then
        write(*,'(a,i0,a,i0,a)') 'ERROR: computation of orbital ', iorb, ' on process ', iproc, ' failed!'
        !quit=.true.
        stop
    end if
end do

!!do iall=1,lin%comsr%sizePhibuffr
!!  write(20000+iproc,*) iall, phibuffr(iall)
!!end do


! No calculate the charge density
! Bounds of slice in global coordinates
is=nscatterarr(iproc,3)-14
ie=is+nscatterarr(iproc,1)-1
istri=0
totalCharge=0.d0
do iorb=1,lin%comsr%noverlaps(iproc)
    iiorb=lin%comsr%overlaps(iorb)
    ilr=lin%comsr%comarr(4,iorb,iproc)
    istrj=0
    do jorb=1,lin%comsr%noverlaps(iproc)
        jjorb=lin%comsr%overlaps(jorb)
        jlr=lin%comsr%comarr(4,jorb,iproc)
        ! Bounds in global coordinates
        i1s=max(2*lin%Llr(ilr)%ns1-14,2*lin%Llr(jlr)%ns1-14)
        i1e=min(2*lin%Llr(ilr)%ns1-14+lin%Llr(ilr)%d%n1i-1,2*lin%Llr(jlr)%ns1-14+lin%Llr(jlr)%d%n1i-1)
        i2s=max(2*lin%Llr(ilr)%ns2-14,2*lin%Llr(jlr)%ns2-14)
        i2e=min(2*lin%Llr(ilr)%ns2-14+lin%Llr(ilr)%d%n2i-1,2*lin%Llr(jlr)%ns2-14+lin%Llr(jlr)%d%n2i-1)
        i3s=max(2*lin%Llr(ilr)%ns3-14,2*lin%Llr(jlr)%ns3-14,is)
        i3e=min(2*lin%Llr(ilr)%ns3-14+lin%Llr(ilr)%d%n3i-1,2*lin%Llr(jlr)%ns3-14+lin%Llr(jlr)%d%n3i-1,ie)
        do i3=i3s,i3e
            i3d=i3-2*lin%Llr(ilr)%ns3
            j3d=i3-2*lin%Llr(jlr)%ns3
            indi3=(i3d+15-1)*lin%Llr(ilr)%d%n2i*lin%Llr(ilr)%d%n1i
            indj3=(j3d+15-1)*lin%Llr(jlr)%d%n2i*lin%Llr(jlr)%d%n1i
            !indl3=(i3-is+15-1)*Glr%d%n2i*Glr%d%n1i
            indl3=(i3-is)*Glr%d%n2i*Glr%d%n1i
            do i2=i2s,i2e
                i2d=i2-2*lin%Llr(ilr)%ns2
                j2d=i2-2*lin%Llr(jlr)%ns2
                indi2=(i2d+15-1)*lin%Llr(ilr)%d%n1i
                indj2=(j2d+15-1)*lin%Llr(jlr)%d%n1i
                indl2=(i2+15-1)*Glr%d%n1i
                do i1=i1s,i1e
                    i1d=i1-2*lin%Llr(ilr)%ns1
                    j1d=i1-2*lin%Llr(jlr)%ns1
                    ! Now calculate the index in the boxes.
                    indri = indi3 + indi2 + i1d+15 + istri
                    indrj = indj3 + indj2 + j1d+15 + istrj
                    indLarge = indl3 + indl2 + i1+15
                    tt = factor*densKern(iorb,jorb)*phibuffr(indri)*phibuffr(indrj)
                    rho(indLarge) = rho(indLarge) + tt
                    !rhofull(indLarge) = rhofull(indLarge) + tt
                    totalCharge = totalCharge + tt
                end do
            end do
        end do
        istrj = istrj + lin%Llr(jlr)%d%n3i*lin%Llr(jlr)%d%n2i*lin%Llr(jlr)%d%n1i
    end do
    istri = istri + lin%Llr(ilr)%d%n3i*lin%Llr(ilr)%d%n2i*lin%Llr(ilr)%d%n1i
end do

!write(*,'(x,a,i5,es20.12)') 'iproc, TOTAL CHARGE = ', iproc, totalCharge*hxh*hyh*hzh
call mpiallred(totalCharge, 1, mpi_sum, mpi_comm_world, ierr)
if(iproc==0) write(*,'(x,a,es20.12)') 'TOTAL CHARGE = ', totalCharge*hxh*hyh*hzh

end subroutine sumrhoForLocalizedBasis




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
    jlr=lin%onWhichAtomAll(jorb)
    ist=ist+lin%Llr(jlr)%wfd%nvctr_c+7*lin%Llr(jlr)%wfd%nvctr_f
end do
commsSumrho(2)=ist

! amount of data to be sent
jlr=lin%onWhichAtomAll(iorb)
commsSumrho(3)=lin%Llr(jlr)%wfd%nvctr_c+7*lin%Llr(jlr)%wfd%nvctr_f

! localization region to which this orbital belongs to
commsSumrho(4)=lin%onWhichAtomAll(iorb)

! to which MPI process should this orbital be sent
commsSumrho(5)=jproc

! the position on the MPI process to which it should be sent
commsSumrho(6)=istDest

! the tag for this communication
commsSumrho(7)=tag

! commsSumrho(8): this entry is used a request for the mpi_isend.

! commsSumrho(9): this entry is used a request for the mpi_irecv.


end subroutine setCommunicationInformation
