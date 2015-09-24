!> @file
!!    Calculate the electronic density (rho)
!! @author
!!    Copyright (C) 2007-2012 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Selfconsistent potential is saved in rhopot, 
!! new arrays rho,pot for calculation of forces ground state electronic density
!! Potential from electronic charge density
subroutine density_and_hpot(dpbox,symObj,orbs,Lzd,pkernel,rhodsc,GPU,xc,psi,rho,vh,rho_ion,hstrten)
  use module_base
  use module_types
  use module_xc
  use module_interfaces, fake_name => density_and_hpot
  use Poisson_Solver, except_dp => dp, except_gp => gp
  implicit none
  !Arguments
  type(denspot_distribution), intent(in) :: dpbox
  type(rho_descriptors),intent(inout) :: rhodsc
  type(orbitals_data), intent(in) :: orbs
  type(local_zone_descriptors), intent(in) :: Lzd
  type(symmetry_data), intent(in) :: symObj
  type(coulomb_operator), intent(inout) :: pkernel
  type(xc_info), intent(in) :: xc
  real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
  type(GPU_pointers), intent(inout) :: GPU
  real(gp), dimension(6), intent(out) :: hstrten
  real(dp), dimension(:), pointer :: rho,vh
  real(dp), dimension(:,:,:,:), pointer :: rho_ion
  !local variables
  character(len=*), parameter :: subname='density_and_hpot'
  real(gp) :: ehart_fake
  real(dp), dimension(:,:), pointer :: rho_p

  if (.not. associated(rho)) then
     if (dpbox%ndimpot>0) then
        rho = f_malloc_ptr(dpbox%ndimpot*orbs%nspin,id='rho')
     else
        rho = f_malloc_ptr(1*orbs%nspin,id='rho')
     end if

     nullify(rho_p)
     call sumrho(dpbox,orbs,Lzd,GPU,symObj,rhodsc,xc,psi,rho_p)
     call communicate_density(dpbox,orbs%nspin,rhodsc,rho_p,rho,.false.)
  end if

  !Calculate the total density in the case of nspin==2
  if (orbs%nspin==2) then
     call axpy(dpbox%ndimpot,1.0_dp,rho(1+dpbox%ndimpot),1,rho(1),1)
  end if
  if (dpbox%ndimpot>0) then
     vh = f_malloc_ptr(dpbox%ndimpot,id='vh')
  else
     vh = f_malloc_ptr(1,id='vh')
  end if

  if (xc%id(1) /= XC_NO_HARTREE) then
     !Calculate electrostatic potential
     call vcopy(dpbox%ndimpot,rho(1),1,vh(1),1)
     !the ionic denisty is given in the case of the embedded solver
     if (pkernel%method /= 'VAC') then
        call H_potential('D',pkernel,vh,vh,ehart_fake,0.0_dp,.true.,stress_tensor=hstrten,rho_ion=rho_ion)
     else
        call H_potential('D',pkernel,vh,vh,ehart_fake,0.0_dp,.false.,stress_tensor=hstrten)
     end if
  else
     !Only to_zero vh
     vh=0.0_dp
     ehart_fake=0.0_dp
  end if

  !In principle symmetrization of the stress tensor is not needed since the density has been 
  !already symmetrized
  if (symObj%symObj >= 0 .and. pkernel%geocode=='P') &
       call symm_stress(hstrten,symObj%symObj)

END SUBROUTINE density_and_hpot


!> Calculates the charge density by summing the square of all orbitals
subroutine sumrho(dpbox,orbs,Lzd,GPU,symObj,rhodsc,xc,psi,rho_p,mapping)
   use module_base
   use module_types
   use module_xc
   use module_interfaces, except_this_one => sumrho
   use yaml_output
   implicit none
   !Arguments
   type(denspot_distribution), intent(in) :: dpbox
   type(rho_descriptors),intent(in) :: rhodsc
   type(orbitals_data), intent(in) :: orbs
   type(local_zone_descriptors), intent(in) :: Lzd
   type(symmetry_data), intent(in) :: symObj
   type(xc_info), intent(in) :: xc
   real(wp), dimension(orbs%npsidim_orbs), intent(in) :: psi
   real(dp), dimension(:,:), pointer :: rho_p
   !real(dp), dimension(max(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*nscatterarr(iproc,1),1),nspin), intent(out), target :: rho
   type(GPU_pointers), intent(inout) :: GPU
   integer,dimension(orbs%norb),intent(in),optional:: mapping
   !Local variables
   character(len=*), parameter :: subname='sumrho'
   logical :: writeout
   integer :: nspinn
   integer :: iorb
   integer,dimension(:),allocatable:: localmapping

   writeout=dpbox%mpi_env%iproc + dpbox%mpi_env%igroup==0 .and. verbose >= 1

   call timing(dpbox%mpi_env%iproc,'Rho_comput    ','ON')

   if (writeout) then
      call yaml_map('GPU acceleration',GPU%OCLconv)
   end if

   !components of the charge density
   if (orbs%nspinor ==4) then
      nspinn=4
   else
      nspinn=orbs%nspin
   end if

   if (associated(rho_p)) then
      call yaml_warning('(sumrho) rho_p already associated, exiting...')
      stop
   end if
   !print *,'here',Lzd%linear,present(mapping),dpbox%iproc_world
   !write(*,*) 'iproc,rhoarray dim', iproc, Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*nrhotot,nspinn
   rho_p = f_malloc_ptr((/ Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*rhodsc%nrhotot , nspinn /),id='rho_p')

   !switch between GPU/CPU treatment of the density
   !here also one might decide to save the value of psir and of its laplacian 
   if (GPU%OCLconv) then
      call local_partial_density_OCL(orbs,rhodsc%nrhotot,Lzd%Glr,&
           dpbox%hgrids(1),dpbox%hgrids(2),dpbox%hgrids(3),orbs%nspin,psi,rho_p,GPU)
   else if(Lzd%linear) then
       if(.not.present(mapping)) then
           if(dpbox%mpi_env%iproc + dpbox%mpi_env%igroup==0) then
              call yaml_newline()
              call yaml_warning('Using fake local mapping array. Check whether this is correct!')
           end if
           !write(*,'(1x,a)') &
           !    'WARNING: mapping is not present, using fake local mapping array. Check whether this is correct!'
           localmapping = f_malloc(orbs%norb,id='localmapping')
           do iorb=1,orbs%norb
               localmapping(iorb)=iorb
           end do
           call local_partial_densityLinear(dpbox%mpi_env%nproc,(rhodsc%icomm==1),dpbox%nscatterarr,rhodsc%nrhotot,&
                Lzd,dpbox%hgrids(1),dpbox%hgrids(2),dpbox%hgrids(3),xc,orbs%nspin,orbs,localmapping,psi,rho_p)
           call f_free(localmapping)
       else
           call local_partial_densityLinear(dpbox%mpi_env%nproc,(rhodsc%icomm==1),dpbox%nscatterarr,rhodsc%nrhotot,&
                Lzd,dpbox%hgrids(1),dpbox%hgrids(2),dpbox%hgrids(3),xc,orbs%nspin,orbs,mapping,psi,rho_p)
       end if
   else
      !initialize the rho array at 10^-20 instead of zero, due to the invcb ABINIT routine
      !otherwise use libXC routine
      call xc_init_rho(xc, Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*rhodsc%nrhotot*nspinn,rho_p,dpbox%mpi_env%nproc)

      !for each of the orbitals treated by the processor build the partial densities
      call local_partial_density(dpbox%mpi_env%nproc,(rhodsc%icomm==1),dpbox%nscatterarr,&
           rhodsc%nrhotot,Lzd%Glr,dpbox%hgrids(1),dpbox%hgrids(2),dpbox%hgrids(3),orbs%nspin,orbs,psi,rho_p)

   end if

   !after validation this point can be deplaced after the allreduce such as to reduce the number of operations
   !probably previous line is not suitable due to the fact that a extra communication would be needed
   if (symObj%symObj >= 0) then
      call symmetrise_density(0,1,Lzd%Glr%geocode,&
           dpbox%ndims(1),dpbox%ndims(2),dpbox%ndims(3),orbs%nspin,rho_p,symObj)
   end if
   call timing(dpbox%mpi_env%iproc,'Rho_comput    ','OF')

 END SUBROUTINE sumrho
 

!> Starting point for the communication routine of the density
subroutine communicate_density(dpbox,nspin,rhodsc,rho_p,rho,keep_rhop)
  use module_base
  use module_types
  use yaml_output
  implicit none
  logical, intent(in) :: keep_rhop !< preserves the total density in the rho_p array
  integer, intent(in) :: nspin
  type(denspot_distribution), intent(in) :: dpbox
  type(rho_descriptors),intent(in) :: rhodsc
  real(dp), dimension(:,:), pointer :: rho_p !< partial density in orbital distribution scheme
  real(dp), dimension(max(dpbox%ndims(1)*dpbox%ndims(2)*dpbox%n3d,1),nspin), intent(out) :: rho
  !local variables
  character(len=*), parameter :: subname='communicate_density'
  logical :: dump
  integer :: i1,i2,i3,i3off,i3s,i,ispin,ierr,j3,j3p,j,itmred,n3d,irho
  real(dp) :: charge,tt,rhotot_dbl
  real(dp), dimension(:,:), allocatable :: tmred
  !!  real(dp), dimension(:,:), allocatable :: rho_p_OCL
  !!  real(dp), dimension(:,:), allocatable :: psi_OCL
  !integer :: ncount0,ncount1,ncount2,ncount3,ncountmpi0,ncountmpi1,ncount_max,ncount_rate
  real(gp),dimension(:,:),allocatable :: dprho_comp
  real(4) ,dimension(:,:),allocatable :: sprho_comp
  
  dump=dpbox%mpi_env%iproc + dpbox%mpi_env%igroup == 0 .and. verbose >= 1

  !write(*,*) 'iproc,TIMING:SR1'!,dpbox%iproc_world,real(ncount1-ncount0)/real(ncount_rate)
  !the density must be communicated to meet the shape of the poisson solver
  if (dpbox%mpi_env%ngroup*dpbox%mpi_env%nproc> 1) then
     call timing(dpbox%mpi_env%iproc,'Rho_commun    ','ON')
     !write(*,*) 'rsflag',rsflag
     !communication strategy for the density
     !LDA case (icomm==1)
     if (rhodsc%icomm==1) then
        if (dump) call yaml_map('Rho Commun','RED_SCT')
        do ispin=1,nspin
            !call system_clock(ncount0,ncount_rate,ncount_max)
           call MPI_REDUCE_SCATTER(rho_p(1,ispin),rho(1,ispin),&
                &   dpbox%ngatherarr(0,3),mpidtypd,MPI_SUM,bigdft_mpi%mpi_comm,ierr)
            !call system_clock(ncount1,ncount_rate,ncount_max)
           !write(*,*) 'TIMING:LDA',ispin!,real(ncount1-ncount0)/real(ncount_rate)
        end do
        !call MPI_BARRIER(dpbox%mpi_comm,ierr)
        !print *,'iproc',dpbox%iproc
        ! splitted single-double precision communication (icomm=2)
     else if (rhodsc%icomm==2) then
        if (dump) call yaml_map('Rho Commun','ALLRED (Mix)')
         rhotot_dbl=0.0d0
         do ispin=1,nspin
           do irho=1, dpbox%ndims(1)*dpbox%ndims(2)*rhodsc%nrhotot
             rhotot_dbl=rhotot_dbl+rho_p(irho,ispin)*product(dpbox%hgrids)!hxh*hyh*hzh
           enddo
        enddo
        call mpiallred(rhotot_dbl,1,MPI_SUM,comm=bigdft_mpi%mpi_comm)

        !call system_clock(ncount0,ncount_rate,ncount_max)

        sprho_comp = f_malloc((/ rhodsc%sp_size, nspin /),id='sprho_comp')
        dprho_comp = f_malloc((/ rhodsc%dp_size, nspin /),id='dprho_comp')
        call compress_rho(rho_p,dpbox%ndimgrid,nspin,rhodsc,sprho_comp,dprho_comp)

        !call system_clock(ncount1,ncount_rate,ncount_max)
        !write(*,*) 'TIMING:ARED1',real(ncount1-ncount0)/real(ncount_rate)
        call mpiallred(sprho_comp,MPI_SUM,comm=bigdft_mpi%mpi_comm)
        call mpiallred(dprho_comp,MPI_SUM,comm=bigdft_mpi%mpi_comm)
        !call system_clock(ncount2,ncount_rate,ncount_max)
        !write(*,*) 'TIMING:ARED2',real(ncount2-ncount1)/real(ncount_rate)

        i3s=dpbox%nscatterarr(dpbox%mpi_env%iproc,3)-dpbox%nscatterarr(dpbox%mpi_env%iproc,4)
        n3d=dpbox%nscatterarr(dpbox%mpi_env%iproc,1)
        call uncompress_rho(sprho_comp,dprho_comp,dpbox%ndims,nspin,rhodsc,rho_p,i3s,n3d)

        !call system_clock(ncount3,ncount_rate,ncount_max)
        !write(*,*) 'TIMING:ARED3',real(ncount3-ncount2)/real(ncount_rate)
         !write(*,*) 'TIMING:MIX',real(ncount3-ncount0)/real(ncount_rate)

        call f_free(sprho_comp)
        call f_free(dprho_comp)

        !naive communication (unsplitted GGA case) (icomm=0)
     else if (rhodsc%icomm==0) then
        if (dump) call yaml_map('Rho Commun','ALLRED')
        call mpiallred(rho_p(1,1),dpbox%ndimgrid*nspin,&
             &   MPI_SUM,comm=bigdft_mpi%mpi_comm)
         !call system_clock(ncount1,ncount_rate,ncount_max)
         !write(*,*) 'TIMING:DBL',real(ncount1-ncount0)/real(ncount_rate)
     else
        STOP 'DENSITY COMMUNICATION KEY UNVALID' 
     endif
     call timing(dpbox%mpi_env%iproc,'Rho_commun    ','OF')
     call timing(dpbox%mpi_env%iproc,'Rho_comput    ','ON')
     if (rhodsc%icomm /= 1) then
        !treatment which includes the periodic GGA
        !the density should meet the poisson solver distribution
        i3s=dpbox%nscatterarr(dpbox%mpi_env%iproc,3)-dpbox%nscatterarr(dpbox%mpi_env%iproc,4)
        n3d=dpbox%nscatterarr(dpbox%mpi_env%iproc,1)
        do ispin=1,nspin
           do i3=1,n3d
              j3=i3+i3s
              j3p=modulo(j3-1,dpbox%ndims(3))+1
              do i2=1,dpbox%ndims(2)
                 do i1=1,dpbox%ndims(1)
                    i=i1+(i2-1)*dpbox%ndims(1)+dpbox%ndims(1)*dpbox%ndims(2)*(i3-1)
                    j=i1+(i2-1)*dpbox%ndims(1)+dpbox%ndims(1)*dpbox%ndims(2)*(j3p-1)
                    rho(i,ispin)=rho_p(j,ispin)
                 end do
              end do
           end do
        end do
     end if
  else
     call timing(dpbox%mpi_env%iproc,'Rho_comput    ','ON')
     call vcopy(dpbox%ndims(1)*dpbox%ndims(2)*dpbox%n3d*nspin,rho_p(1,1),1,&
          rho(1,1),1)
  end if
  !print *,'okhere',dpbox%iproc_world
  ! Check
  tt=0.d0
  i3off=dpbox%ndims(1)*dpbox%ndims(2)*dpbox%i3xcsh!Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*nscatterarr(iproc,4)

  !allocation of the magnetic density orientation array
  if (dpbox%mpi_env%nproc > 1) then
     itmred=2
  else
     itmred=1
  end if
  !use this check also for the magnetic density orientation
  tmred = f_malloc((/ nspin+1, itmred /),id='tmred')



  tmred(nspin+1,itmred)=0.0_dp
  do ispin=1,nspin!n
     tmred(ispin,itmred)=0.0_dp
     do i=1,dpbox%ndims(1)*dpbox%ndims(2)*dpbox%nscatterarr(dpbox%mpi_env%iproc,2)!dpbox%ndimpot!Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*nscatterarr(iproc,2)
        !!        tt=tt+rho(i+i3off,ispin)
        tmred(ispin,itmred)=tmred(ispin,itmred)+rho(i+i3off,ispin)
        !temporary check for debugging purposes
        !!        if (rho(i+i3off,ispin)/rho(i+i3off,ispin) /= 1.d0) then
        !!           print *,iproc,'error in density construction',rho(i+i3off,ispin)
        !!        end if
     enddo
     tmred(nspin+1,itmred)=tmred(nspin+1,itmred)+tmred(ispin,itmred)
  end do

  if (.not. keep_rhop) then
     call f_free_ptr(rho_p)
  end if
  
  !print *,'okhere1',dpbox%iproc_world
  !call mpi_barrier(bigdft_mpi%mpi_comm,ierr)
  if (dpbox%mpi_env%nproc > 1) then
     call timing(dpbox%mpi_env%iproc,'Rho_comput    ','OF')
     call timing(dpbox%mpi_env%iproc,'Rho_commun    ','ON')
     !!     call MPI_REDUCE(tt,charge,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,bigdft_mpi%mpi_comm,ierr)
     call MPI_REDUCE(tmred(1,2),tmred(1,1),nspin+1,mpidtypd,MPI_SUM,0,dpbox%mpi_env%mpi_comm,ierr)
     call timing(dpbox%mpi_env%iproc,'Rho_commun    ','OF')
     call timing(dpbox%mpi_env%iproc,'Rho_comput    ','ON')
  endif
  !print *,'okthere',dpbox%iproc_world
  !call mpi_barrier(bigdft_mpi%mpi_comm,ierr)
  !write the results
  if (dump) then
     if(nspin==4) then
        charge=tmred(1,1)
        tt=sqrt(tmred(2,1)**2+tmred(3,1)**2+tmred(4,1)**2)
     else
        charge=0._dp
        do ispin=1,nspin
           charge=charge+tmred(ispin,1)
        end do
     end if

     !yaml output
     call yaml_map('Total electronic charge',real(charge,gp)*product(dpbox%hgrids),fmt='(f21.12)')

     if (rhodsc%icomm==2 .and. dpbox%mpi_env%nproc*dpbox%mpi_env%ngroup > 1) then
        call yaml_map('Electronic charge changed by rho compression',&
             abs(rhotot_dbl-real(charge,gp)*product(dpbox%hgrids)),fmt='(1pe21.12)')
     endif
     if(nspin == 4 .and. tt > 0._dp) then
        call yaml_map('Magnetic density orientation',&
             (/(tmred(ispin,1)/tmred(1,1),ispin=2,nspin)/),fmt='(f10.4)')
     end if
     call yaml_newline()
  end if
!call yaml_mapping_close()

  call f_free(tmred)

  call timing(dpbox%mpi_env%iproc,'Rho_comput    ','OF')

end subroutine communicate_density

!> Here starts the routine for building partial density inside the localisation region
!! This routine should be treated as a building-block for the linear scaling code
subroutine local_partial_density(nproc,rsflag,nscatterarr,&
      &   nrhotot,lr,hxh,hyh,hzh,nspin,orbs,psi,rho_p)
   use module_base
   use module_types
   use module_interfaces
   use locreg_operations
   implicit none
   logical, intent(in) :: rsflag
   integer, intent(in) :: nproc,nrhotot
   integer, intent(in) :: nspin
   real(gp), intent(in) :: hxh,hyh,hzh
   type(orbitals_data), intent(in) :: orbs
   type(locreg_descriptors), intent(in) :: lr
   integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
   real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(in) :: psi
   real(dp), dimension(lr%d%n1i,lr%d%n2i,nrhotot,max(nspin,orbs%nspinor)), intent(inout) :: rho_p
   !local variables
   character(len=*), parameter :: subname='local_partial_density'
   integer :: iorb !n(c) i1,i2,i3,ii 
   integer :: oidx,sidx,nspinn,npsir,ncomplex
   real(gp) :: hfac,spinval
   type(workarr_sumrho) :: w
   real(wp), dimension(:,:), allocatable :: psir

   call initialize_work_arrays_sumrho(1,[lr],.true.,w)

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

   psir = f_malloc((/ lr%d%n1i*lr%d%n2i*lr%d%n3i, npsir /),id='psir')
   !initialisation
   !print *,iproc,'there'
   if (lr%geocode == 'F') call f_zero(psir)

   do iorb=1,orbs%norbp
      !print *,'norbp',orbs%norbp,orbs%norb,orbs%nkpts,orbs%kwgts,orbs%iokpt,orbs%occup
      hfac=orbs%kwgts(orbs%iokpt(iorb))*(orbs%occup(orbs%isorb+iorb)/(hxh*hyh*hzh))
      spinval=orbs%spinsgn(orbs%isorb+iorb)

      if (hfac /= 0.d0) then

         !sum for complex function case, npsir=1 in that case
         do oidx=0,ncomplex

!!$      print *,'iorb+orbs%isorb,SUMRHO',iorb+orbs%isorb,&
!!$                sum(psi(:,1,iorb))

            do sidx=1,npsir
               call daub_to_isf(lr,w,psi(1,oidx+sidx,iorb),psir(1,sidx))
            end do

            !print *,'iorb,nrm',iorb,npsir,&
            !     nrm2(lr%d%n1i*lr%d%n2i*lr%d%n3i*npsir,psir(1,1),1)

            select case(lr%geocode)
            case('F')

               call partial_density_free(rsflag,nproc,lr%d%n1i,lr%d%n2i,lr%d%n3i,&
                  &   npsir,nspinn,nrhotot,&
                  &   hfac,nscatterarr,spinval,psir,rho_p,lr%bounds%ibyyzz_r)

            case('P')

               call partial_density(rsflag,nproc,lr%d%n1i,lr%d%n2i,lr%d%n3i,&
                  &   npsir,nspinn,nrhotot,&
                  &   hfac,nscatterarr,spinval,psir,rho_p)

            case('S')

               call partial_density(rsflag,nproc,lr%d%n1i,lr%d%n2i,lr%d%n3i,&
                  &   npsir,nspinn,nrhotot,&
                  &   hfac,nscatterarr,spinval,psir,rho_p)

            end select

         end do
      end if

      !print *,'iorb,nrmBBB',iorb,npsir,&
      !     nrm2(lr%d%n1i*lr%d%n2i*nrhotot*max(nspin,orbs%nspinor),rho_p(1,1,1,1),1)


   enddo
   

   call f_free(psir)
   call deallocate_work_arrays_sumrho(w)

END SUBROUTINE local_partial_density


subroutine partial_density(rsflag,nproc,n1i,n2i,n3i,npsir,nspinn,nrhotot,&
      &   hfac,nscatterarr,spinsgn,psir,rho_p)
   use module_base
   use module_types
   use yaml_output
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
   !$omp parallel default(none) &
   !$omp private(jproc,i1,i3s,i1s,i1e,i3off,n3d,i3,j3,isjmp,psisq,p1,p2,p3,p4,r1,r2,r3,r4,i2) &
   !$omp shared(n1i,nproc,rsflag,nspinn,nscatterarr,spinsgn) &
   !$omp shared(n2i,npsir,hfac,hfac2,psir,rho_p,n3i,i3sg)
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
      call yaml_warning('Problem with rho_p, i3s=' // trim(yaml_toa(i3sg)) // &
           & 'nrhotot=' // trim(yaml_toa(nrhotot)))
      !write(*,'(1x,a,i0,1x,i0)')'ERROR: problem with rho_p: i3s,nrhotot,',i3sg,nrhotot
      stop
   end if

END SUBROUTINE partial_density


subroutine partial_density_free(rsflag,nproc,n1i,n2i,n3i,npsir,nspinn,nrhotot,&
      &   hfac,nscatterarr,spinsgn,psir,rho_p,&
      &   ibyyzz_r) 
   use module_base
   use module_types
   use yaml_output
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
                  !write(93000+iorb,'(4i8,es20.10)') i1,i2,i3,isjmp,rho_p(i1,i2,i3s,isjmp)
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
      call yaml_warning('Problem with rho_p, i3s=' // trim(yaml_toa(i3sg)) // &
           & 'nrhotot=' // trim(yaml_toa(nrhotot)))
      !write(*,'(1x,a,i0,1x,i0)')'ERROR: problem with rho_p: i3s,nrhotot,',i3sg,nrhotot
      stop
   end if

   !  call system_clock(ncount1,ncount_rate,ncount_max)
   !  write(*,*) 'TIMING:PDF',real(ncount1-ncount0)/real(ncount_rate)
END SUBROUTINE partial_density_free


!> Symmetrise the density using the symmetry operation
subroutine symmetrise_density(iproc,nproc,geocode,n1i,n2i,n3i,nspin,rho,& !n(c) nscatterarr (arg:6)
     sym)
  use module_base!, only: gp,dp,wp,ndebug,memocc
  use module_types
  use m_ab6_symmetry
  use yaml_output, only: yaml_warning
  implicit none
  !Arguments
  integer, intent(in) :: iproc,nproc,nspin, n1i, n2i, n3i
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  !n(c) integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  real(dp), dimension(n1i,n2i,n3i,nspin), intent(inout) :: rho
  type(symmetry_data), intent(in) :: sym
  !local variables
  character(len=*), parameter :: subname='symmetrise_density'
  integer :: errno, ispden, nsym_used, nSym, isym, imagn, r2,inzee,isign, n2i_eff
  integer :: nd2, izone_max, numpt, izone, rep, nup, iup, ind, j, j1, j2, j3,i1,i2,i3, i2_eff
  real(dp) :: rhosu1, rhosu2
  real(dp), dimension(:,:), allocatable :: rhosu12
  real(dp), dimension(:,:,:,:,:), allocatable :: rhog
  integer, parameter :: ncache = 4 * 1024
  real(dp), dimension(:,:,:), allocatable :: zw
  integer, pointer :: symRel(:,:,:)
  integer, pointer :: symAfm(:)
  real(gp), pointer :: transNon(:,:)

  call symmetry_get_matrices_p(sym%symObj, nSym, symRel, transNon, symAfm, errno)
  if (nSym == 1) return
  if (geocode == 'F') then
     !call yaml_warning('The symmetrization of the density is not implemented for the isolated systems')
     return
  end if

!!$  ! Array sizes for the real-to-complex FFT: note that n1(there)=n1(here)+1
!!$  ! and the same for n2,n3. Not needed for the moment
!!$  call dimensions_fft(n1,n2,n3,nd1,nd2,nd3,n1f,n3f,n1b,n3b,nd1f,nd3f,nd1b,nd3b)
  n2i_eff = n2i
  if (geocode == "S") then
     n2i_eff = 1
     zw = f_malloc((/ 2, ncache/4, 2 /),id='zw')
     !use this check also for the magnetic density orientation
     rhog = f_malloc((/ 2, n1i+1, 1, n3i+1, 2 /),id='rhog')
  else
     !use this check also for the magnetic density orientation
     rhog = f_malloc((/ 2, n1i+1, n2i+1, n3i+1, 2 /),id='rhog')
  end if


  ! Here imagn doesn't change since nspin == nsppol in BigDFT.
  imagn = 1

  !  Treat either full density, spin-up density or magnetization
  !  Note the decrease of ispden to the value 1, in order to finish
  !  with rhog of the total density (and not the spin-up density or magnetization)
  do ispden=nspin,1,-1

     do i2_eff = 0, (n2i - n2i_eff), 1
        imagn = 1 + i2_eff

        !    Prepare the density to be symmetrized, in the reciprocal space
        nsym_used=0
        do isym=1,nSym
           if(symAfm(isym)==1)nsym_used=nsym_used+1
        end do

        !    rhor -fft-> rhog    (rhog is used as work space)
        !    Note : it should be possible to reuse rhog in the antiferromagnetic case
        !    this would avoid one FFT
        ! fft the input array x:
        rhog=0.0_gp !put to avoid fpe in the FFT
        do i3=0,n3i-1
           do i2=0,n2i_eff-1
              do i1=0,n1i-1
                 rhog(1,i1+1,i2+1,i3+1,1)=rho(i1+1,i2_eff+i2+1,i3+1,ispden)
                 rhog(2,i1+1,i2+1,i3+1,1)=0.d0
              enddo
           enddo
        enddo

        inzee=1
        isign=-1

        if (geocode /= "S") then
           call fft(n1i,n2i_eff,n3i,n1i+1,n2i_eff+1,n3i+1,rhog,isign,inzee)
        else
           call fft2d(n1i,n3i,n1i+1,n3i+1,rhog,isign,inzee,zw,ncache)
        end if

        !!     work(:)=rho(:,ispden)
        !!     call fourdp(cplex,rhog,work,-1,mpi_enreg,nfft,ngfft,paral_kgb,0)

        !this should be modified via the nscatter array
        !for the moment we can put nproc=1
        nd2=n2i_eff/nproc

        !    The following is only valid for total, up or dn density
        !    -------------------------------------------------------

        !    Get maxvalue of izone
        izone_max=count(sym%irrzon(:,2,imagn)>0)
        rhosu12 = f_malloc((/ 2, izone_max /),id='rhosu12')

        numpt=0
        do izone=1,n1i*n2i_eff*n3i

           !      Get repetition number
           rep=sym%irrzon(izone,2,imagn)
           if(rep==0)exit

           !      Compute number of unique points in this symm class:
           nup=nsym_used/rep

           !      Accumulate charge over equivalent points
           rhosu1=0._dp
           rhosu2=0._dp
           do iup=1,nup
              ind=sym%irrzon(iup+numpt,1,imagn)
              j=ind-1
              j1=modulo(j,n1i)
              j2=modulo(j/n1i,n2i_eff)
              j3=j/(n1i*n2i_eff)
              r2=modulo(j2,nd2)
              !here we should insert the condition that the planes should belong to iproc
              if(modulo(j/n1i,n2i_eff)/nd2==iproc) then ! this ind is to be treated by me_fft
                 ind=n1i*(nd2*j3+r2)+j1+1 !this is ind in the current proc
                 rhosu1=rhosu1+rhog(1,j1+1,r2+1,j3+1,inzee)*sym%phnons(1,iup+numpt,imagn)&
                      &   -rhog(2,j1+1,r2+1,j3+1,inzee)*sym%phnons(2,iup+numpt,imagn)
                 rhosu2=rhosu2+rhog(2,j1+1,r2+1,j3+1,inzee)*sym%phnons(1,iup+numpt,imagn)&
                      &   +rhog(1,j1+1,r2+1,j3+1,inzee)*sym%phnons(2,iup+numpt,imagn)
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
        !call mpiallred(rhosu12(1,1),2*izone_max,MPI_SUM,bigdft_mpi%mpi_comm,ierr)

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
        do izone=1,n1i*n2i_eff*n3i

           !      Get repetition number
           rep=sym%irrzon(izone,2,imagn)
           if(rep==0)exit

           !      Compute number of unique points in this symm class:
           nup=nsym_used/rep

           !      Define symmetrized rho(G) at equivalent points:
           do iup=1,nup
              ind=sym%irrzon(iup+numpt,1,imagn)
              !        decompose ind-1=n1(n2 j3+ j2)+j1
              j=ind-1
              j1=modulo(j,n1i)
              j2=modulo(j/n1i,n2i_eff)
              j3=j/(n1i*n2i_eff)
              r2=modulo(j2,nd2)
              if(modulo(j/n1i,n2i_eff)/nd2==iproc) then ! this ind is to be treated by me_fft
                 !          ind in the proc ind-1=n1(nd2 j3+ r2)+j1
                 ind=n1i*(nd2*j3+r2)+j1+1 !this is ind in the current proc
                 rhog(1,j1+1,r2+1,j3+1,inzee)=rhosu12(1,izone)*sym%phnons(1,iup+numpt,imagn)&
                      &   +rhosu12(2,izone)*sym%phnons(2,iup+numpt,imagn)
                 rhog(2,j1+1,r2+1,j3+1,inzee)=rhosu12(2,izone)*sym%phnons(1,iup+numpt,imagn)&
                      &   -rhosu12(1,izone)*sym%phnons(2,iup+numpt,imagn)
              end if
           end do

           !      Keep index of how many points have been considered:
           numpt=numpt+nup

           !      End loop over izone
        end do

        call f_free(rhosu12)

        !    Pull out full or spin up density, now symmetrized
        !!     call fourdp(cplex,rhog,work,1,mpi_enreg,nfft,ngfft,paral_kgb,0)

        isign=1
        if (geocode /= "S") then
           call fft(n1i,n2i_eff,n3i,n1i+1,n2i_eff+1,n3i+1,rhog,isign,inzee)
        else
           call fft2d(n1i,n3i,n1i+1,n3i+1,rhog,isign,inzee,zw,ncache)
        end if

        do i3=0,n3i-1
           do i2=0,n2i_eff-1
              do i1=0,n1i-1
                 !correct the density in case it has negative values
                 rho(i1+1,i2_eff+i2+1,i3+1,ispden)=max(rhog(1,i1+1,i2+1,i3+1,inzee)/real(n1i*n2i_eff*n3i,dp),1.d-20)
              enddo
           enddo
        enddo
        !divide by the number of grid points
        !rho(:,ispden)=work(:)

     end do

  end do ! ispden

  call f_free(rhog)
  if (geocode == "S") then
     call f_free(zw)
  end if

END SUBROUTINE symmetrise_density


!> Compress the electronic density
subroutine compress_rho(rho_p,ndimgrid,nspin,rhodsc,sprho_comp,dprho_comp)
   use module_base
   use module_types
   implicit none
   type(rho_descriptors),intent(in) :: rhodsc
   integer,intent(in) :: nspin,ndimgrid
   real(gp),dimension(ndimgrid,nspin),intent(in) :: rho_p                 !< the partial rho array of the current proc
   real(gp),dimension(rhodsc%dp_size,nspin),intent(out) :: dprho_comp     !< compressed arrays of rho in double 
   real(kind=4),dimension(rhodsc%sp_size,nspin),intent(out) :: sprho_comp !< compressed arrays of rho in single 
   integer :: irho,jrho,iseg,ispin

   do ispin=1,nspin
      !$omp parallel default(shared) private(irho,jrho)
      !$omp do
      do iseg=1,rhodsc%n_csegs
         jrho=rhodsc%cseg_b(iseg)
         do irho=rhodsc%spkey(iseg,1),rhodsc%spkey(iseg,2)
            sprho_comp(jrho,ispin)=real(rho_p(irho,ispin),kind=4)
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
END SUBROUTINE compress_rho


!> Restore the necessary rho planes for using in GGA
subroutine uncompress_rho(sprho_comp,dprho_comp,ndims,nspin,rhodsc,rho_uncomp,i3s,n3d)
   use module_base
   use module_types
   implicit none
   integer,intent(in) :: nspin,i3s,n3d
   type(rho_descriptors),intent(in) :: rhodsc
   integer, dimension(3), intent(in) :: ndims
   real(gp),dimension(rhodsc%dp_size,nspin),intent(in) :: dprho_comp
   real(4), dimension(rhodsc%sp_size,nspin),intent(in) :: sprho_comp
   real(gp),dimension(ndims(1)*ndims(2)*ndims(3),nspin), intent(out) :: rho_uncomp
   !local variables
   integer :: irho,jrho,ispin,iseg,ibegin,iend,j3p !n(c) imax
   logical :: overlap

   ! starting and endding point of rho array of interest
   !n(c) imax=ndims(1)*ndims(2)*ndims(3)
   j3p=modulo(i3s,ndims(3))+1
   ibegin=(j3p-1)*ndims(1)*ndims(2)+1
   j3p=modulo(i3s+n3d-1,ndims(3))+1
   iend= j3p*ndims(1)*ndims(2)

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
END SUBROUTINE uncompress_rho


!> Restore the necessary rho planes for using in GGA
subroutine uncompress_rho_old(sprho_comp,dprho_comp,&
      &   lr,nspin,rhodsc,rho_uncomp,i3s,n3d)
   use module_base
   use module_types
   implicit none
   integer,intent(in) :: nspin,i3s,n3d
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
END SUBROUTINE uncompress_rho_old


subroutine rho_segkey(iproc,at,rxyz,crmult,frmult,&
      &   n1i,n2i,n3i,hxh,hyh,hzh,nspin,rhodsc,iprint)
   use module_base
   use module_types
   implicit none
   integer,intent(in) :: n1i,n2i,n3i,iproc,nspin
   type(atoms_data), intent(in) :: at
   real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
   real(gp), intent(in) :: crmult,frmult,hxh,hyh,hzh
   !real(gp), dimension(at%astruct%ntypes,3), intent(in) :: radii_cf
   logical,intent(in) :: iprint
   type(rho_descriptors),intent(inout) :: rhodsc
   !local variables
   real(gp), parameter :: epsilon=1.e-10_gp
   integer :: i1,i2,i3,iseg,irho,iat !n(c) ispin, i_all,jrho, nseg
   integer :: reg_c,reg_l
   !these give stack overflow!!!!
   !integer, dimension(n1i*n2i*n3i) :: reg
   !integer,dimension(n1i*n2i*n3i,2) :: dpkey,spkey
   integer :: n_fsegs,n_csegs
   character(len=*), parameter :: subname='rhokey'
   integer :: nbx,nby,nbz,nl1,nl2,nl3,nat
   real(gp) :: dpmult,dsq,spadd
   integer :: i1min,i1max,i2min,i2max,i3min,i3max,nrhomin,nrhomax
   integer,dimension(at%astruct%nat) :: i1fmin,i1fmax,i2fmin,i2fmax,i3fmin,i3fmax
   integer,dimension(at%astruct%nat) :: i1cmin,i1cmax,i2cmin,i2cmax,i3cmin,i3cmax!,dsq_cr,dsq_fr
   real(gp), dimension(at%astruct%nat) :: dsq_cr,dsq_fr
   integer :: csegstot,fsegstot,corx,cory,corz,ithread,nthreads
   integer, dimension(:), allocatable :: reg
   integer, dimension(:,:), allocatable :: dpkey,spkey
   !integer :: ncount0,ncount1,ncount2,ncount3,ncount4,ncount_rate,ncount_max, i_stat
   !$ integer :: omp_get_thread_num,omp_get_num_threads

   reg = f_malloc(n1i*n2i*n3i,id='reg')
   dpkey = f_malloc((/ n1i*n2i*n3i, 2 /),id='dpkey')
   spkey = f_malloc((/ n1i*n2i*n3i, 2 /),id='spkey')


   ithread=0
   nthreads=1

   rhodsc%geocode=at%astruct%geocode
   nat=at%astruct%nat

   !parameter to adjust the single precision and double precision regions
   spadd=5.0_gp
   dpmult=1.0_gp


   ! calculate the corrections of the grid when transforming from 
   ! n1,n2,n3 to n1i, n2i, n3i
   call gridcorrection(nbx,nby,nbz,nl1,nl2,nl3,at%astruct%geocode)

   corx=nl1+nbx+1
   cory=nl2+nby+1
   corz=nl3+nbz+1 

   ! set the boundaries of the regions that we will determine the segment structure
   call get_boxbound(at,rxyz,crmult,frmult,hxh,hyh,hzh,spadd,dpmult,&
      &   n1i,n2i,n3i,corx,cory,corz,i1min,i1max,i2min,i2max,i3min,i3max)

   nrhomin = (i3min-1)*n1i*n2i+(i2min-1)*n1i+i1min
   nrhomax = (i3max-1)*n1i*n2i+(i2max-1)*n1i+i1max

   n_fsegs=0
   n_csegs=0

   !$omp parallel default(none)&
   !$omp private(irho,i1,i2,i3,ithread,nthreads)&
   !$omp shared(i1min,i1max,i2min,i2max,i3min,i3max,reg,n1i,n2i,n3i)
   !$ ithread = omp_get_thread_num()
   !$ nthreads = omp_get_num_threads()
   !LG: the initialization have to be done in all the array
      do i3=1,n3i!i3min,i3max
      if (mod(i3,nthreads).eq.ithread) then
        do i2=1,n2i!i2min,i2max
          do i1=1,n1i!i1min,i1max
            irho = (i3-1)*n1i*n2i+(i2-1)*n1i+i1
            reg(irho)=0
          enddo
        enddo
      endif
      enddo
   !$omp end parallel

   do iat=1,nat
      call get_atbound(iat,at,rxyz,crmult,frmult,hxh,&
         &   hyh,hzh,spadd,dpmult,n1i,n2i,n3i,corx,cory,corz,&
         &   i1cmin(iat),i1cmax(iat),i2cmin(iat),i2cmax(iat),i3cmin(iat),i3cmax(iat),&
         &   i1fmin(iat),i1fmax(iat),i2fmin(iat),i2fmax(iat),i3fmin(iat),i3fmax(iat))

      dsq_cr(iat)=(at%radii_cf(at%astruct%iatype(iat),1)*crmult+hxh*spadd)**2
      dsq_fr(iat)=(at%radii_cf(at%astruct%iatype(iat),2)*frmult*dpmult)**2
   enddo

   !$omp parallel default(none)&
   !$omp private(iat,irho,dsq,i1,i2,i3,ithread,nthreads)&
   !$omp shared(nat,rxyz,hxh,hyh,hzh,dsq_cr,dsq_fr,reg)&
   !$omp shared(i1cmin,i1cmax,i2cmin,i2cmax,i3cmin,i3cmax)&
   !$omp shared(n1i,n2i,n3i,corx,cory,corz,iproc)
   !$ ithread = omp_get_thread_num()
   !$ nthreads = omp_get_num_threads()

   do iat=1,nat
      do i3=i3cmin(iat),i3cmax(iat)
      if (mod(i3,nthreads).eq.ithread) then
         do i2=i2cmin(iat),i2cmax(iat)
             do i1=i1cmin(iat),i1cmax(iat)
               dsq=(rxyz(1,iat)-real(i1-corx,gp)*hxh)**2+&
                  &   (rxyz(2,iat)-real(i2-cory,gp)*hyh)**2+&
                  &   (rxyz(3,iat)-real(i3-corz,gp)*hzh)**2
               irho = (i3-1)*n1i*n2i+(i2-1)*n1i+i1
               if(dsq-epsilon.lt.dsq_cr(iat)) then
                  reg(irho)=1
               endif
               !write(17+iproc,*)i1,i2,i3,dsq,dsq_cr(iat),irho,reg(irho)
            enddo
         enddo
      endif
      enddo
   enddo
   !$omp end parallel

   !$omp parallel default(none)&
   !$omp private(iat,irho,dsq,i1,i2,i3,ithread,nthreads)&
   !$omp shared(nat,rxyz,hxh,hyh,hzh,dsq_cr,dsq_fr,reg)&
   !$omp shared(i1fmin,i1fmax,i2fmin,i2fmax,i3fmin,i3fmax,iproc)&
   !$omp shared(n1i,n2i,n3i,corx,cory,corz)
   !$ ithread = omp_get_thread_num()
   !$ nthreads = omp_get_num_threads()
   do iat=1,nat
      do i3=i3fmin(iat),i3fmax(iat)
      if (mod(i3,nthreads).eq.ithread) then
         do i2=i2fmin(iat),i2fmax(iat)
            do i1=i1fmin(iat),i1fmax(iat)
               dsq=(rxyz(1,iat)-real(i1-corx,gp)*hxh)**2+&
                  &   (rxyz(2,iat)-real(i2-cory,gp)*hyh)**2+&
                  &   (rxyz(3,iat)-real(i3-corz,gp)*hzh)**2          
               irho = (i3-1)*n1i*n2i+(i2-1)*n1i+i1
               !write(17+iproc,*)i1,i2,i3,dsq,dsq_fr(iat)
               if(dsq-epsilon.lt.dsq_fr(iat)) then
                  reg(irho)=2
               endif
               !write(17+iproc,*)i1,i2,i3,dsq,dsq_fr(iat),irho,reg(irho)
            enddo
         enddo
      endif
      enddo
   enddo
   !$omp end parallel
   !print *,'here',nrhomin,nrhomax,iproc,spadd,dpmult
   !write(17+iproc,*)dsq_fr,dsq_cr,i1fmin,i1fmax,i2fmin,i2fmax,i3fmin,i3fmax
   !write(17+iproc,*)i1cmin,i1cmax,i2cmin,i2cmax,i3cmin,i3cmax,nrhomin,nrhomax,spadd
   !write(17+iproc,*)corx,cory,corz,i1min,i1max,i2min,i2max,i3min,i3max
   !write(17+iproc,*)i1min,i1max,i2min,i2max,i3min,i3max
   !do iat=1,nat
   !   write(17+iproc,*)i1cmin(iat),i1cmax(iat),i2cmin(iat),i2cmax(iat),i3cmin(iat),i3cmax(iat)
   !   write(17+iproc,*)i1fmin(iat),i1fmax(iat),i2fmin(iat),i2fmax(iat),i3fmin(iat),i3fmax(iat)
   !end do
   !do irho=nrhomin,nrhomax
   !   write(17+iproc,*)irho,reg(irho)
   !end do
   !call MPI_BARRIER(bigdft_mpi%mpi_comm,i_stat)
   !stop
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

   rhodsc%dpkey = f_malloc_ptr((/ n_fsegs, 2 /),id='rhodsc%dpkey')
   rhodsc%spkey = f_malloc_ptr((/ n_csegs, 2 /),id='rhodsc%spkey')
   rhodsc%cseg_b = f_malloc_ptr(n_csegs,id='rhodsc%cseg_b')
   rhodsc%fseg_b = f_malloc_ptr(n_fsegs,id='rhodsc%fseg_b')

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

   call f_free(reg)
   call f_free(dpkey)
   call f_free(spkey)


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
            &   'nspin,spadd,dpmult                                 :',nspin,spadd,dpmult
         write(1001,*) 'spadd and dpmult can be found in "/src/sumrho.f90"'
         write(1001,*) '----------------------------------------------------------------------'
         write(1001,*) 'Number of single precision data points              :',&
            &   rhodsc%sp_size*nspin
         write(1001,*) 'Number of double precision data points              :',&
            &   rhodsc%dp_size*nspin
         write(1001,*) 'Total number of data points                         :',&
            &   n1i*n2i*n3i*nspin
         write(1001,'(1X,A,1X,F5.2)') &
            &   'Estimated compression ratio, number of data points  :',&
            &   real(rhodsc%sp_size*nspin  +rhodsc%dp_size*nspin)/(n1i*n2i*n3i*nspin)
         write(1001,'(1X,A,1X,F5.2)') &
            &   'Estimated compression ratio, data volume to be sent :',&
            &   real(rhodsc%sp_size*nspin/2+rhodsc%dp_size*nspin)/(n1i*n2i*n3i*nspin)
         write(1001,*) ''
         close(1001)
      endif
   endif
   !call system_clock(ncount4,ncount_rate,ncount_max)
   !write(*,*) 'TIMING:RHOKEY4',real(ncount4-ncount3)/real(ncount_rate)
   !write(*,*) 'TIMING:RHOKEYA',real(ncount4-ncount0)/real(ncount_rate)
END SUBROUTINE rho_segkey


subroutine gridcorrection(nbx,nby,nbz,nl1,nl2,nl3,geocode)
   implicit none
   character(len=1),intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
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
END SUBROUTINE gridcorrection


subroutine get_boxbound(at,rxyz,crmult,frmult,hxh,hyh,hzh,spadd,dpmult,&
      &   n1i,n2i,n3i,corx,cory,corz,i1min,i1max,i2min,i2max,i3min,i3max)
   use module_base
   use module_types
   implicit none
   real(gp),intent(in) :: dpmult,spadd
   integer,intent(in) :: n1i,n2i,n3i,corx,cory,corz
   type(atoms_data), intent(in) :: at
   real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
   real(gp), intent(in) :: crmult,frmult,hxh,hyh,hzh
   !real(gp), dimension(at%astruct%ntypes,3), intent(in) :: radii_cf
   integer,intent(out) :: i1min,i1max,i2min,i2max,i3min,i3max
   integer,dimension(at%astruct%nat,6) :: crbound,frbound
   real(gp) :: sprad,dprad
   integer :: iat

   i1min=n1i
   i2min=n2i
   i3min=n3i
   i1max=1
   i2max=1
   i3max=1

!   write(*,*) 'i1min,i2min:',i1min,i2min,i3min,i1max,i2max,i3max

!write (*,*) 'hxh,hyh,hzh',hxh,hyh,hzh

   ! setup up the cubes that determines the single and double precision segments
   do iat=1,at%astruct%nat
      sprad=at%radii_cf(at%astruct%iatype(iat),1)*crmult+hxh*spadd
      dprad=at%radii_cf(at%astruct%iatype(iat),2)*frmult*dpmult
      crbound(iat,1)=floor((rxyz(1,iat)-sprad)/hxh)+corx
!write(*,*) 'crbound(iat,1)',crbound(iat,1)
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
!write(*,*) 'iat=',iat,i1min,i1max,i2min,i2max,i3min,i3max
   enddo
END SUBROUTINE get_boxbound


subroutine get_atbound(iat,at,rxyz,crmult,frmult,hxh,hyh,hzh,&
      &   spadd,dpmult,n1i,n2i,n3i,corx,cory,corz,&
      &   i1cmin,i1cmax,i2cmin,i2cmax,i3cmin,i3cmax,&
      &   i1fmin,i1fmax,i2fmin,i2fmax,i3fmin,i3fmax)
   use module_base
   use module_types
   implicit none
   real(gp),intent(in) :: dpmult,spadd
   integer,intent(in) :: n1i,n2i,n3i,iat,corx,cory,corz
   type(atoms_data), intent(in) :: at
   real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
   real(gp), intent(in) :: crmult,frmult,hxh,hyh,hzh
   !real(gp), dimension(at%astruct%ntypes,3), intent(in) :: radii_cf
   integer,intent(out) :: i1cmin,i1cmax,i2cmin,i2cmax,i3cmin,i3cmax
   integer,intent(out) :: i1fmin,i1fmax,i2fmin,i2fmax,i3fmin,i3fmax
   real(gp) :: sprad,dprad

   sprad=at%radii_cf(at%astruct%iatype(iat),1)*crmult+hxh*spadd
   dprad=dpmult*frmult*at%radii_cf(at%astruct%iatype(iat),2)

   i1cmin=floor((rxyz(1,iat)-sprad)/hxh)+corx+1
   i1cmax=floor((rxyz(1,iat)+sprad)/hxh)+corx-1
   i2cmin=floor((rxyz(2,iat)-sprad)/hyh)+cory+1
   i2cmax=floor((rxyz(2,iat)+sprad)/hyh)+cory-1
   i3cmin=floor((rxyz(3,iat)-sprad)/hzh)+corz+1
   i3cmax=floor((rxyz(3,iat)+sprad)/hzh)+corz-1

   i1fmin=floor((rxyz(1,iat)-dprad)/hxh)+corx+1
   i1fmax=floor((rxyz(1,iat)+dprad)/hxh)+corx-1
   i2fmin=floor((rxyz(2,iat)-dprad)/hyh)+cory+1
   i2fmax=floor((rxyz(2,iat)+dprad)/hyh)+cory-1
   i3fmin=floor((rxyz(3,iat)-dprad)/hzh)+corz+1
   i3fmax=floor((rxyz(3,iat)+dprad)/hzh)+corz-1

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

END SUBROUTINE get_atbound


subroutine is_overlap(a,b,x,y,overlap)
   implicit none
   integer,intent(in) :: a,b,x,y
   logical,intent(out) :: overlap

   if (((a.le.x).and.(x.le.b)).or.&
      &   ((a.le.y).and.(y.le.b))) then
   overlap=.true.
else
   overlap=.false.
endif
END SUBROUTINE is_overlap
