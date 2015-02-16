!> @file
!!  Routine to calculate the action of the hamiltonian
!! @author
!!   Copyright (C) 2005-2013 BigDFT group 
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 


!> Calculate the action of the local hamiltonian on the orbitals
subroutine local_hamiltonian(iproc,nproc,npsidim_orbs,orbs,Lzd,hx,hy,hz,&
     ipotmethod,confdatarr,pot,psi,hpsi,pkernel,xc,alphaSIC,ekin_sum,epot_sum,eSIC_DC,&
     dpbox,potential,comgp)
  use module_base
  use module_types
  use module_interfaces, except_this_one => local_hamiltonian
  use module_xc
  implicit none
  !Arguments
  integer, intent(in) :: iproc,nproc,npsidim_orbs
  integer, intent(in) :: ipotmethod !< Method which has to be chosen for applying the potential to the wavefunctions in the real space form:
                                    !! 0 is the traditional potential application
                                    !! 1 is the application of the exact exchange (which has to be precomputed and stored in the potential array)
                                    !! 2 is the application of the Perdew-Zunger SIC
                                    !! 3 is the application of the Non-Koopman's correction SIC
  real(gp), intent(in) :: hx,hy,hz,alphaSIC
  type(orbitals_data), intent(in) :: orbs
  type(local_zone_descriptors), intent(in) :: Lzd
  type(confpot_data), dimension(orbs%norbp), intent(in) :: confdatarr
  type(xc_info), intent(in) :: xc
  real(wp), dimension(npsidim_orbs), intent(in) :: psi              !< This dimension will be modified
  real(wp), dimension(:),pointer :: pot                             !< the potential, with the dimension compatible with the ipotmethod flag
  real(gp), intent(out) :: ekin_sum,epot_sum,eSIC_DC
  real(wp), dimension(npsidim_orbs), intent(inout) :: hpsi
  type(coulomb_operator), intent(in) :: pkernel                     !< the PSolver kernel which should be associated for the SIC schemes
  type(denspot_distribution),intent(in),optional :: dpbox
  real(wp), dimension(*), intent(in), optional, target :: potential !< Distributed potential. Might contain the density for the SIC treatments
  type(p2pComms),intent(inout), optional:: comgp
  !!real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i*nspin) :: pot
  !!real(wp), dimension(max(dpbox%ndimrhopot,orbs%nspin)), intent(in), optional, target :: potential !< Distributed potential. Might contain the density for the SIC treatments
  !local variables
  character(len=*), parameter :: subname='local_hamiltonian'
  logical :: dosome
  integer :: iorb,npot,ispot,ispsi,ilr,ilr_orb!,jproc,ierr
  real(wp) :: exctXcoeff
  real(gp) :: ekin,epot,kx,ky,kz,eSICi,eSIC_DCi !n(c) etest
  type(workarr_locham) :: wrk_lh
  real(wp), dimension(:,:), allocatable :: vsicpsir
  real(wp), dimension(:,:), allocatable :: psir
  !!write(*,*) 'condition',(present(dpbox) .and. present(potential) .and. present(comgp))

  call f_routine(id='local_hamiltonian')

  epot=0.d0
  ekin=0.d0

  !some checks
  exctXcoeff=xc_exctXfac(xc)

  if (exctXcoeff /= 0.0_gp .neqv. ipotmethod ==1) then
     if (iproc==0) write(*,*)&
          'ERROR (local_hamiltonian): potential method not compatible with exact exchange'
     stop
  end if

  if (.not.(associated(pkernel%kernel) .and. alphaSIC /=0.0_gp) .and. ipotmethod == 2) then
     if (iproc==0) write(*,*)&
          'ERROR (local_hamiltonian): potential method not compatible with SIC'
     stop
  end if

  ekin_sum=0.0_gp
  epot_sum=0.0_gp
  eSIC_DC=0.0_gp
!!$do jproc=0,nproc-1
!!$call MPI_BARRIER(bigdft_mpi%mpi_comm,ierr)
!!$if (jproc==iproc) then
  !loop on the localisation regions (so to create one work array set per lr)
  loop_lr: do ilr=1,Lzd%nlr
    !check if this localisation region is used by one of the orbitals
    dosome=.false.
    do iorb=1,orbs%norbp
      dosome = (orbs%inwhichlocreg(iorb+orbs%isorb) == ilr)
      if (dosome) exit
    end do
    if (.not. dosome) cycle loop_lr
      
    !components of the potential (four or one, depending on the spin)
    npot=orbs%nspinor
    if (orbs%nspinor == 2) npot=1
   
    ! Wavefunction in real space
    psir = f_malloc0((/ Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i, orbs%nspinor /),id='psir')

    call initialize_work_arrays_locham(1,Lzd%Llr(ilr),orbs%nspinor,.true.,wrk_lh)  
  
    ! wavefunction after application of the self-interaction potential
    if (ipotmethod == 2 .or. ipotmethod == 3) then
      vsicpsir = f_malloc((/ Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i, orbs%nspinor /),id='vsicpsir')
    end if

    ispsi=1
    loop_orbs: do iorb=1,orbs%norbp
      ilr_orb=orbs%inwhichlocreg(iorb+orbs%isorb)
      if (ilr_orb /= ilr) then
        ispsi=ispsi+&
             (Lzd%Llr(ilr_orb)%wfd%nvctr_c+7*Lzd%Llr(ilr_orb)%wfd%nvctr_f)*orbs%nspinor
        cycle loop_orbs
      end if
      
!!$      print *,'iorb+orbs%isorb,BEFORE',iorb+orbs%isorb,&
!!$                sum(psi(ispsi:&
!!$                ispsi+(Lzd%Llr(ilr_orb)%wfd%nvctr_c+7*Lzd%Llr(ilr_orb)%wfd%nvctr_f)*orbs%nspinor-1))

        
      call daub_to_isf_locham(orbs%nspinor,Lzd%Llr(ilr),wrk_lh,psi(ispsi),psir(1,1))

      !calculate the ODP, to be added to VPsi array
   
      !Perdew-Zunger SIC scheme
      eSIC_DCi=0.0_gp
      if (ipotmethod == 2) then
         !in this scheme the application of the potential is already done
         call PZ_SIC_potential(iorb,Lzd%Llr(ilr),orbs,xc,&
              0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,pkernel,psir,vsicpsir,eSICi,eSIC_DCi)
      !NonKoopmans' correction scheme
      else if (ipotmethod == 3) then 
         !in this scheme first we have calculated the potential then we apply it
         call vcopy(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,&
              psir(1,1),1,vsicpsir(1,1),1)
         !for the moment the ODP is supposed to be valid only with one lr
         call psir_to_vpsi(npot,orbs%nspinor,Lzd%Llr(ilr),&
              pot(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspin+&
              (iorb-1)*Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor+1),&
              vsicpsir,eSICi)
      end if
   
      call psir_to_vpsi(npot,orbs%nspinor,Lzd%Llr(ilr),&
           pot(orbs%ispot(iorb)),psir(1,1),epot,confdata=confdatarr(iorb))

!!$      !this ispot has to be better defined inside denspot structure
!!$      print *,'orbs, epot',orbs%isorb+iorb,epot,orbs%ispot(iorb),&
!!$           sum(pot(orbs%ispot(iorb):orbs%ispot(iorb):Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i-1))
   
      !ODP treatment (valid only for the nlr=1 case)
      if (ipotmethod==1) then !Exact Exchange
         ispot=1+Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*(orbs%nspin+iorb-1)
         !add to the psir function the part of the potential coming from the exact exchange
         call axpy(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,exctXcoeff,pot(ispot),1,psir(1,1),1)
      else if (ipotmethod == 2) then !PZ scheme
         !subtract the sic potential from the vpsi function
         call axpy(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,-alphaSIC,vsicpsir(1,1),1,psir(1,1),1)
         !add the SIC correction to the potential energy
         epot=epot-alphaSIC*eSICi
         !accumulate the Double-Counted SIC energy
         eSIC_DC=eSIC_DC+alphaSIC*eSIC_DCi
      else if (ipotmethod == 3) then !NK scheme
         !add the sic potential from the vpsi function
         call axpy(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,alphaSIC,vsicpsir(1,1),1,psir(1,1),1)
         epot=epot+alphaSIC*eSICi
         !accumulate the Double-Counted SIC energy
         eSIC_DC=eSIC_DC+alphaSIC*orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*eSICi
      end if
   
      !apply the kinetic term, sum with the potential and transform back to Daubechies basis
      !k-point values, if present
      kx=orbs%kpts(1,orbs%iokpt(iorb))
      ky=orbs%kpts(2,orbs%iokpt(iorb))
      kz=orbs%kpts(3,orbs%iokpt(iorb))
 
      call isf_to_daub_kinetic(hx,hy,hz,kx,ky,kz,orbs%nspinor,Lzd%Llr(ilr),wrk_lh,&
           psir(1,1),hpsi(ispsi),ekin)

!!$      print *,'iorb+orbs%isorb,AFTER',iorb+orbs%isorb,&
!!$                sum(hpsi(ispsi:&
!!$                ispsi+(Lzd%Llr(ilr_orb)%wfd%nvctr_c+7*Lzd%Llr(ilr_orb)%wfd%nvctr_f)*orbs%nspinor-1))
    

      ekin_sum=ekin_sum+orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*ekin
      epot_sum=epot_sum+orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*epot
           !print *,'iorb+orbs%isorb',iorb+orbs%isorb,ekin,epot
      ispsi=ispsi+&
           (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor
      !print *,'iorb,epot',orbs%isorb+iorb,epot

    enddo loop_orbs
   
    !deallocations of work arrays
    call f_free(psir)

    if (ipotmethod == 2 .or. ipotmethod ==3) then
       call f_free(vsicpsir)
    end if
    call deallocate_work_arrays_locham(wrk_lh)
   
  end do loop_lr
!!$end if
!!$end do

  call f_release_routine()

END SUBROUTINE local_hamiltonian

!> Calculate the action of the local potential on the orbitals
!! @param ipotmethod Indicates the method which has to be chosen for applying the potential to the wavefunctions in the 
!!                   real space form:
!!                   0 is the traditional potential application
!!                   1 is the application of the exact exchange (which has to be precomputed and stored in the potential array)
!!                   2 is the application of the Perdew-Zunger SIC
!!                   3 is the application of the Non-Koopman's correction SIC
subroutine psi_to_vlocpsi(iproc,npsidim_orbs,orbs,Lzd,&
     ipotmethod,confdatarr,pot,psi,vpsi,pkernel,xc,alphaSIC,epot_sum,evSIC,vpsi_noconf,econf_sum)
  use module_base
  use module_types
  use module_interfaces, except_this_one => psi_to_vlocpsi
  use module_xc
  implicit none
  integer, intent(in) :: iproc,ipotmethod,npsidim_orbs
  real(gp), intent(in) :: alphaSIC
  type(xc_info), intent(in) :: xc
  type(orbitals_data), intent(in) :: orbs
  type(local_zone_descriptors), intent(in) :: Lzd
  type(confpot_data), dimension(orbs%norbp), intent(in) :: confdatarr
  real(wp), dimension(npsidim_orbs), intent(in) :: psi !this dimension will be modified
  real(wp), dimension(*) :: pot !< the potential, with the dimension compatible with the ipotmethod flag
  real(gp), intent(out) :: epot_sum,evSIC
  real(wp), dimension(npsidim_orbs), intent(inout) :: vpsi
  type(coulomb_operator), intent(in) :: pkernel !< the PSolver kernel which should be associated for the SIC schemes
  real(wp), dimension(npsidim_orbs), intent(inout),optional :: vpsi_noconf
  real(gp),intent(out),optional :: econf_sum
  !local variables
  character(len=*), parameter :: subname='psi_to_vlocpsi'
  logical :: dosome
  integer :: iorb,npot,ispot,ispsi,ilr,ilr_orb,nbox,nvctr,ispinor
  real(wp) :: exctXcoeff
  real(gp) :: epot,eSICi,eSIC_DCi,econf !n(c) etest
  type(workarr_sumrho) :: w
  real(wp), dimension(:,:), allocatable :: psir,vsicpsir,psir_noconf

  call f_routine(id='psi_to_vlocpsi')

  !some checks
  exctXcoeff=xc_exctXfac(xc)

  if (exctXcoeff /= 0.0_gp .neqv. ipotmethod ==1) then
     if (iproc==0) write(*,*)&
          'ERROR (local_hamiltonian): potential method not compatible with exact exchange'
     stop
  end if

  if (.not.(associated(pkernel%kernel) .and. alphaSIC /=0.0_gp) .and. ipotmethod == 2) then
     if (iproc==0) write(*,*)&
          'ERROR (local_hamiltonian): potential method not compatible with SIC'
     stop
  end if

  epot_sum=0.0_gp
  evSIC=0.0_gp
  if (present(econf_sum)) then
      econf_sum=0.0_gp
  end if

  call initialize_work_arrays_sumrho(lzd%nlr,lzd%llr,.true.,w)

  !loop on the localisation regions (so to create one work array set per lr)
  loop_lr: do ilr=1,Lzd%nlr
     !check if this localisation region is used by one of the orbitals
     dosome=.false.
     do iorb=1,orbs%norbp
        dosome = (orbs%inwhichlocreg(iorb+orbs%isorb) == ilr)
        if (dosome) then
            exit
        end if
     end do
     if (.not. dosome) cycle loop_lr

     !initialise the work arrays
     call initialize_work_arrays_sumrho(1,lzd%llr(ilr),.false.,w)

     !box elements size
     nbox=Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i

     !components of the potential (four or one, depending on the spin)
     npot=orbs%nspinor
     if (orbs%nspinor == 2) npot=1

     ! Wavefunction in real space
     psir = f_malloc0((/ nbox, orbs%nspinor /),id='psir')

     if (present(vpsi_noconf)) then
         psir_noconf = f_malloc((/ nbox, orbs%nspinor /),id='psir_noconf')
     end if

     !call to_zero(nbox*orbs%nspinor,psir(1,1))

     ! wavefunction after application of the self-interaction potential
     if (ipotmethod == 2 .or. ipotmethod == 3) then
        vsicpsir = f_malloc((/ nbox, orbs%nspinor /),id='vsicpsir')
     end if

  !n(c) etest=0.0_gp

  ispsi=1
  loop_orbs: do iorb=1,orbs%norbp
     ilr_orb=orbs%inwhichlocreg(iorb+orbs%isorb)
     nvctr=Lzd%Llr(ilr_orb)%wfd%nvctr_c+7*Lzd%Llr(ilr_orb)%wfd%nvctr_f
     if (ilr_orb /= ilr) then
        ispsi=ispsi+nvctr*orbs%nspinor
        cycle loop_orbs
     end if
     
     !transform the wavefunction in Daubechies basis to the wavefunction in ISF basis
     !the psir wavefunction is given in the spinorial form
     do ispinor=1,orbs%nspinor
        call daub_to_isf(Lzd%Llr(ilr),w,psi(ispsi+nvctr*(ispinor-1)),psir(1,ispinor))
     end do

     !ispot=1+lr%d%n1i*lr%d%n2i*lr%d%n3i*(nspin+iorb-1)
     !etest=etest+dot(lr%d%n1i*lr%d%n2i*lr%d%n3i,pot(ispot),1,psir(1,1),1)
     !print *,'epot, iorb,iproc,norbp',iproc,orbs%norbp,iorb

     !calculate the ODP, to be added to VPsi array

     !Perdew-Zunger SIC scheme
     eSIC_DCi=0.0_gp
     if (ipotmethod == 2) then
        !in this scheme the application of the potential is already done
        call PZ_SIC_potential(iorb,Lzd%Llr(ilr),orbs,xc,&
             0.5_gp*Lzd%hgrids(1),0.5_gp*Lzd%hgrids(2),0.5_gp*Lzd%hgrids(3),&
             pkernel,psir,vsicpsir,eSICi,eSIC_DCi)
     !NonKoopmans' correction scheme
     else if (ipotmethod == 3) then 
        !in this scheme first we have calculated the potential then we apply it
        call vcopy(nbox*orbs%nspinor,psir(1,1),1,vsicpsir(1,1),1)
        !for the moment the ODP is supposed to be valid only with one lr
        call psir_to_vpsi(npot,orbs%nspinor,Lzd%Llr(ilr),&
             pot(nbox*(orbs%nspin+(iorb-1)*orbs%nspinor)+1),&
             vsicpsir,eSICi)
     end if

     !apply the potential to the psir wavefunction and calculate potential energy
     if (present(vpsi_noconf)) then
         if (.not.present(econf_sum)) then
             stop 'ERROR: econf must be present when psir_noconf is present'
         end if
         call vcopy(nbox*orbs%nspinor, psir(1,1), 1, psir_noconf(1,1), 1)
         call psir_to_vpsi(npot,orbs%nspinor,Lzd%Llr(ilr),&
              pot(orbs%ispot(iorb)),psir,epot,confdata=confdatarr(iorb),vpsir_noconf=psir_noconf,econf=econf)
     else
         call psir_to_vpsi(npot,orbs%nspinor,Lzd%Llr(ilr),&
              pot(orbs%ispot(iorb)),psir,epot,confdata=confdatarr(iorb))
     end if
     !!do i_stat=1,lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i
     !!    write(1000+ilr_orb,*) orbs%ispot(iorb)+i_stat-1, pot(orbs%ispot(iorb)+i_stat-1)
     !!end do
     !this ispot has to be better defined inside denspot structure

     !ODP treatment (valid only for the nlr=1 case)
     if (ipotmethod==1) then !Exact Exchange
        ispot=1+nbox*(orbs%nspin+iorb-1)
        !add to the psir function the part of the potential coming from the exact exchange
        call axpy(nbox,exctXcoeff,pot(ispot),1,psir(1,1),1)
     else if (ipotmethod == 2) then !PZ scheme
        !subtract the sic potential from the vpsi function
        call axpy(nbox*orbs%nspinor,-alphaSIC,vsicpsir(1,1),1,psir(1,1),1)
        !add the SIC correction to the potential energy
        epot=epot-alphaSIC*eSICi
        !accumulate the Double-Counted SIC energy
        evSIC=evSIC+alphaSIC*eSIC_DCi
     else if (ipotmethod == 3) then !NK scheme
        !add the sic potential from the vpsi function
        call axpy(nbox*orbs%nspinor,alphaSIC,vsicpsir(1,1),1,psir(1,1),1)
        epot=epot+alphaSIC*eSICi
        !accumulate the Double-Counted SIC energy
        evSIC=evSIC+alphaSIC*orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*eSICi
     end if

     do ispinor=1,orbs%nspinor
         call isf_to_daub(Lzd%Llr(ilr),w,psir(1,ispinor),vpsi(ispsi+nvctr*(ispinor-1)))
         if (present(vpsi_noconf)) then
             call isf_to_daub(Lzd%Llr(ilr),w,psir_noconf(1,ispinor),vpsi_noconf(ispsi+nvctr*(ispinor-1)))
         end if
     end do

     epot_sum=epot_sum+orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*epot
     if (present(econf_sum)) then
         econf_sum=econf_sum+orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*econf
     end if
     ispsi=ispsi+nvctr*orbs%nspinor
  enddo loop_orbs

  !deallocations of work arrays
  call f_free(psir)
  if (present(vpsi_noconf)) then
      call f_free(psir_noconf)
  end if
  if (ipotmethod == 2 .or. ipotmethod ==3) then
     call f_free(vsicpsir)
  end if

end do loop_lr

call deallocate_work_arrays_sumrho(w)

call f_release_routine()

END SUBROUTINE psi_to_vlocpsi


subroutine psi_to_kinpsi(iproc,npsidim_orbs,orbs,lzd,psi,hpsi,ekin_sum)
  use module_base
  use module_types
  use module_interfaces, except_this_one => psi_to_kinpsi
  implicit none
  integer, intent(in) :: iproc,npsidim_orbs
  type(orbitals_data), intent(in) :: orbs
  type(local_zone_descriptors), intent(in) :: Lzd
  real(wp), dimension(npsidim_orbs), intent(in) :: psi
  real(gp), intent(out) :: ekin_sum
  real(wp), dimension(npsidim_orbs), intent(inout) :: hpsi

  !local variables
  character(len=*), parameter :: subname='psi_to_kinpsi'
  logical :: dosome
  integer :: iorb,ispsi,ilr,ilr_orb
  real(gp) :: ekin
  type(workarr_locham) :: wrk_lh
  real(wp), dimension(:,:), allocatable :: psir
  real(gp) :: kx,ky,kz


  ekin=0.d0
  ekin_sum=0.0_gp

  call initialize_work_arrays_locham(lzd%nlr,lzd%llr,orbs%nspinor,.true.,wrk_lh)  

  !loop on the localisation regions (so to create one work array set per lr)
  loop_lr: do ilr=1,Lzd%nlr
    !check if this localisation region is used by one of the orbitals
    dosome=.false.
    do iorb=1,orbs%norbp
      dosome = (orbs%inwhichlocreg(iorb+orbs%isorb) == ilr)
      if (dosome) then
          exit
      end if
    end do
    if (.not. dosome) cycle loop_lr
   
    ! Wavefunction in real space
    psir = f_malloc0((/ Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i, orbs%nspinor /),id='psir')

    !initialise the work arrays
    call initialize_work_arrays_locham(1,Lzd%Llr(ilr),orbs%nspinor,.false.,wrk_lh)  

   
    ispsi=1
    loop_orbs: do iorb=1,orbs%norbp
      ilr_orb=orbs%inwhichlocreg(iorb+orbs%isorb)
      if (ilr_orb /= ilr) then
        ispsi=ispsi+&
             (Lzd%Llr(ilr_orb)%wfd%nvctr_c+7*Lzd%Llr(ilr_orb)%wfd%nvctr_f)*orbs%nspinor
        cycle loop_orbs
      end if
        
      !call daub_to_isf_locham(orbs%nspinor,Lzd%Llr(ilr),wrk_lh,psi(ispsi),psir(1,1))

      kx=orbs%kpts(1,orbs%iokpt(iorb))
      ky=orbs%kpts(2,orbs%iokpt(iorb))
      kz=orbs%kpts(3,orbs%iokpt(iorb))

      !call isf_to_daub_kinetic(lzd%hgrids(1),lzd%hgrids(2),lzd%hgrids(3),kx,ky,kz,orbs%nspinor,Lzd%Llr(ilr),wrk_lh,&
      !      psir(1,1),hpsi(ispsi),ekin)
      call psi_to_tpsi(lzd%hgrids,orbs%kpts(1,orbs%iokpt(iorb)),orbs%nspinor,&
           Lzd%Llr(ilr),psi(ispsi),wrk_lh,hpsi(ispsi),ekin)
   
      ekin_sum=ekin_sum+orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*ekin

      ispsi=ispsi+&
           (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor

    enddo loop_orbs

    call f_free(psir)

  end do loop_lr

  call deallocate_work_arrays_locham(wrk_lh)


end subroutine psi_to_kinpsi



!> apply the potential to the psir wavefunction and calculate potential energy
subroutine psir_to_vpsi(npot,nspinor,lr,pot,vpsir,epot,confdata,vpsir_noconf,econf)
  use module_base
  use module_types
  use module_interfaces, except_this_one => psir_to_vpsi
  implicit none
  integer, intent(in) :: npot,nspinor
  type(locreg_descriptors), intent(in) :: lr !< localization region of the wavefunction
  real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,npot), intent(in) :: pot
  real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,nspinor), intent(inout) :: vpsir
  real(gp), intent(out) :: epot
  type(confpot_data), intent(in), optional :: confdata !< data for the confining potential
  real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,nspinor), intent(inout), optional :: vpsir_noconf !< wavefunction with  the potential without confinement applied
  real(gp), intent(out),optional :: econf !< confinement energy
  !local variables
  logical :: confining
  integer, dimension(3) :: ishift !temporary variable in view of wavefunction creation

  call f_routine(id='psir_to_vpsi')

  !write(*,'(a,a4,2l5)') 'in psir_to_vpsi: lr%geocode, present(vpsir_noconf), present(econf)', lr%geocode, present(vpsir_noconf), present(econf)

  epot=0.0_gp
  ishift=(/0,0,0/)
  confining=present(confdata)
  if (confining) confining= (confdata%potorder /=0)

  if (confining) then
     if (lr%geocode == 'F') then
        if (present(vpsir_noconf)) then
            if (.not.present(econf)) stop 'ERROR: econf must be present when vpsir_noconf is present!'
            !call apply_potential_lr(lr%d%n1i,lr%d%n2i,lr%d%n3i,&
            call apply_potential_lr_conf_noconf(lr%d%n1i,lr%d%n2i,lr%d%n3i,&
                 lr%d%n1i,lr%d%n2i,lr%d%n3i,&
                 ishift,lr%d%n2,lr%d%n3,&
                 nspinor,npot,vpsir,pot,epot,&
                 confdata,lr%bounds%ibyyzz_r,vpsir_noconf,econf)
                 !confdata=confdata,ibyyzz_r=lr%bounds%ibyyzz_r,psir_noconf=vpsir_noconf,econf=econf)
        else
           !call apply_potential_lr(lr%d%n1i,lr%d%n2i,lr%d%n3i,&
           call apply_potential_lr_conf(lr%d%n1i,lr%d%n2i,lr%d%n3i,&
                 lr%d%n1i,lr%d%n2i,lr%d%n3i,&
                 ishift,lr%d%n2,lr%d%n3,&
                 nspinor,npot,vpsir,pot,epot,&
                 confdata,lr%bounds%ibyyzz_r)
                 !confdata=confdata,ibyyzz_r=lr%bounds%ibyyzz_r)
        end if
     else
        !!!if (present(vpsir_noconf)) then
        !!!if (.not.present(econf)) stop 'ERROR: econf must be present when vpsir_noconf is present!'
        !!!!call apply_potential_lr(lr%d%n1i,lr%d%n2i,lr%d%n3i,&
        !!!    call apply_potential_lr_conf_noconf_nobounds(lr%d%n1i,lr%d%n2i,lr%d%n3i,&
        !!!         lr%d%n1i,lr%d%n2i,lr%d%n3i,&
        !!!         ishift,lr%d%n2,lr%d%n3,&
        !!!         nspinor,npot,vpsir,pot,epot,&
        !!!         confdata,vpsir_noconf,econf)
        !!!         !confdata=confdata)
        !!!else
            call apply_potential_lr_conf_nobounds(lr%d%n1i,lr%d%n2i,lr%d%n3i,&
                 lr%d%n1i,lr%d%n2i,lr%d%n3i,&
                 ishift,lr%d%n2,lr%d%n3,&
                 nspinor,npot,vpsir,pot,epot,&
                 confdata)
                 !confdata=confdata)
        !!! end if
     end if

  else
     
     if (lr%geocode == 'F') then
        !call apply_potential_lr(lr%d%n1i,lr%d%n2i,lr%d%n3i,&
        call apply_potential_lr_bounds(lr%d%n1i,lr%d%n2i,lr%d%n3i,&
             lr%d%n1i,lr%d%n2i,lr%d%n3i,&
             ishift,lr%d%n2,lr%d%n3,&
             nspinor,npot,vpsir,pot,epot,&
             lr%bounds%ibyyzz_r)
        !     ibyyzz_r=lr%bounds%ibyyzz_r)
     else
        !call apply_potential_lr(lr%d%n1i,lr%d%n2i,lr%d%n3i,&
        call apply_potential_lr_nobounds(lr%d%n1i,lr%d%n2i,lr%d%n3i,&
             lr%d%n1i,lr%d%n2i,lr%d%n3i,&
             ishift,lr%d%n2,lr%d%n3,&
             nspinor,npot,vpsir,pot,epot)
     end if
  end if

  call f_release_routine()

end subroutine psir_to_vpsi


!>   Transpose the wavefunction into a real and imaginary part to be treated with k-points
!!   to be used only when nspinor=2 or 4
!!   here the dimensions are n1->n1+1
subroutine transpose_for_kpoints(nspinor,n1,n2,n3,x,ww,direct)
  use module_base
  implicit none
  logical, intent(in) :: direct
  integer, intent(in) :: nspinor,n1,n2,n3
  real(wp), dimension(nspinor*n1*n2*n3), intent(inout) :: x,ww
  !local variables
  integer :: i1,i2,i3,idx,id,id2,id3,isd,ispinor,it

  !k-points also admitted in non-collinear case
  if (direct) then
     do ispinor=1,nspinor/2
        isd=(ispinor-1)*2*n1*n2*n3
        do idx=1,2
           do i3=1,n3
              id3=(i3-1)*n1*n2
              do i2=1,n2
                 id2=(i2-1)*n1
                 do i1=1,n1
                    id=i1+id2+id3+(idx-1)*n1*n2*n3+isd
                    it=idx+2*(i1-1)+2*id2+2*id3+isd
                    ww(it)=x(id)
                 end do
              end do
           end do
        end do
     end do
  else
     do ispinor=1,nspinor/2
        isd=(ispinor-1)*2*n1*n2*n3
        do idx=1,2
           do i3=1,n3
              id3=(i3-1)*n1*n2
              do i2=1,n2
                 id2=(i2-1)*n1
                 do i1=1,n1
                    id=i1+id2+id3+(idx-1)*n1*n2*n3+isd
                    it=idx+2*(i1-1)+2*id2+2*id3+isd
                    ww(id)=x(it)
                 end do
              end do
           end do
        end do
     end do
  end if
  
  !for mixed precision code it should be changed
  call vcopy(nspinor*n1*n2*n3,ww(1),1,x(1),1)
END SUBROUTINE transpose_for_kpoints


!>   routine for applying the local potentials
!!   supports the non-collinear case, the buffer for tails and different Boundary Conditions
!!   Optimal also for the complex wavefuntion case
!!   might generalize the buffers for two different localization regions, provided that the potential lies in a bigger region
subroutine apply_potential(n1,n2,n3,nl1,nl2,nl3,nbuf,nspinor,npot,psir,pot,epot,&
     ibyyzz_r) !optional
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nl1,nl2,nl3,nbuf,nspinor,npot
  real(wp), dimension(-14*nl1:2*n1+1+15*nl1,-14*nl2:2*n2+1+15*nl2,-14*nl3:2*n3+1+15*nl3,nspinor), intent(inout) :: psir
  real(wp), dimension(-14*nl1:2*n1+1+15*nl1-4*nbuf,-14*nl2:2*n2+1+15*nl2-4*nbuf,&
       -14*nl3:2*n3+1+15*nl3-4*nbuf,npot), intent(in) :: pot
  integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r
  real(gp), intent(out) :: epot
  !local variables
  integer :: i1,i2,i3,i1s,i1e,ispinor
  real(wp) :: tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42,tt
  real(wp) :: psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4
  real(gp) :: epot_p
 
  !the Tail treatment is allowed only in the Free BC case
  if (nbuf /= 0 .and. nl1*nl2*nl3 == 0) stop 'NONSENSE: nbuf/=0 only for Free BC'

  epot=0.0_wp
!$omp parallel default(private)&
!$omp shared(pot,psir,n1,n2,n3,epot,ibyyzz_r,nl1,nl2,nl3,nbuf,nspinor)
  !case without bounds
  i1s=-14*nl1
  i1e=2*n1+1+15*nl1
  epot_p=0._gp
!$omp do
  do i3=-14*nl3,2*n3+1+15*nl3
     if (i3 >= -14+2*nbuf .and. i3 <= 2*n3+16-2*nbuf) then !check for the nbuf case
        do i2=-14*nl2,2*n2+1+15*nl2
           if (i2 >= -14+2*nbuf .and. i2 <= 2*n2+16-2*nbuf) then !check for the nbuf case
              !this if statement is inserted here for avoiding code duplication
              !it is to be seen whether the code results to be too much unoptimised
              if (present(ibyyzz_r)) then
                 !in this case we are surely in Free BC
                 !the min is to avoid to calculate for no bounds
                 do i1=-14+2*nbuf,min(ibyyzz_r(1,i2,i3),ibyyzz_r(2,i2,i3))-14-1
                    psir(i1,i2,i3,:)=0.0_wp
                 enddo
                 i1s=max(ibyyzz_r(1,i2,i3)-14,-14+2*nbuf)
                 i1e=min(ibyyzz_r(2,i2,i3)-14,2*n1+16-2*nbuf)
              end if
              !here we put the branchments wrt to the spin
              if (nspinor == 4) then
                 do i1=i1s,i1e
                    !wavefunctions
                    psir1=psir(i1,i2,i3,1)
                    psir2=psir(i1,i2,i3,2)
                    psir3=psir(i1,i2,i3,3)
                    psir4=psir(i1,i2,i3,4)
                    !potentials
                    pot1=pot(i1-2*nbuf,i2-2*nbuf,i3-2*nbuf,1)
                    pot2=pot(i1-2*nbuf,i2-2*nbuf,i3-2*nbuf,2)
                    pot3=pot(i1-2*nbuf,i2-2*nbuf,i3-2*nbuf,3)
                    pot4=pot(i1-2*nbuf,i2-2*nbuf,i3-2*nbuf,4)

                    !diagonal terms
                    tt11=pot1*psir1 !p1
                    tt22=pot1*psir2 !p2
                    tt33=pot4*psir3 !p3
                    tt44=pot4*psir4 !p4
                    !Rab*Rb
                    tt13=pot2*psir3 !p1
                    !Iab*Ib
                    tt14=pot3*psir4 !p1
                    !Rab*Ib
                    tt23=pot2*psir4 !p2
                    !Iab*Rb
                    tt24=pot3*psir3 !p2
                    !Rab*Ra
                    tt31=pot2*psir1 !p3
                    !Iab*Ia
                    tt32=pot3*psir2 !p3
                    !Rab*Ia
                    tt41=pot2*psir2 !p4
                    !Iab*Ra
                    tt42=pot3*psir1 !p4

                    ! Change epot later
                    epot_p=epot_p+tt11*psir1+tt22*psir2+tt33*psir3+tt44*psir4+&
                         2.0_gp*tt31*psir3-2.0_gp*tt42*psir4+2.0_gp*tt41*psir4+2.0_gp*tt32*psir3

                    !wavefunction update
                    !p1=h1p1+h2p3-h3p4
                    !p2=h1p2+h2p4+h3p3
                    !p3=h2p1+h3p2+h4p3
                    !p4=h2p2-h3p1+h4p4
                    psir(i1,i2,i3,1)=tt11+tt13-tt14
                    psir(i1,i2,i3,2)=tt22+tt23+tt24
                    psir(i1,i2,i3,3)=tt33+tt31+tt32
                    psir(i1,i2,i3,4)=tt44+tt41-tt42
                 end do
              else
                 do ispinor=1,nspinor
                    do i1=i1s,i1e
                       !the local potential is always real
                       tt=pot(i1-2*nbuf,i2-2*nbuf,i3-2*nbuf,1)*psir(i1,i2,i3,ispinor)
                       epot_p=epot_p+real(tt*psir(i1,i2,i3,ispinor),gp)
                       psir(i1,i2,i3,ispinor)=tt
                    end do
                 end do
              end if
              
              if (present(ibyyzz_r)) then
                 !the max is to avoid the calculation for no bounds
                 do i1=max(ibyyzz_r(1,i2,i3),ibyyzz_r(2,i2,i3))-14+1,2*n1+16-2*nbuf
                    psir(i1,i2,i3,:)=0.0_wp
                 enddo
              end if

           else
              do i1=-14,2*n1+16
                 psir(i1,i2,i3,:)=0.0_wp
              enddo
           endif
        enddo
     else
        do i2=-14,2*n2+16
           do i1=-14,2*n1+16
              psir(i1,i2,i3,:)=0.0_wp
           enddo
        enddo
     endif
  enddo
!$omp end do

!$omp critical
  epot=epot+epot_p
!$omp end critical

!$omp end parallel

END SUBROUTINE apply_potential

!>   Calculate on-the fly each projector for each atom, then applies the projectors 
!!   to all distributed orbitals
subroutine applyprojectorsonthefly(iproc,orbs,at,lr,&
     rxyz,hx,hy,hz,wfd,nlpsp,psi,hpsi,eproj_sum,&
     paw)
  use module_base
  use module_types
  use yaml_output
  use psp_projectors, only: PSPCODE_PAW
  implicit none
  integer, intent(in) :: iproc
  real(gp), intent(in) :: hx,hy,hz
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(DFT_PSP_projectors), intent(inout) :: nlpsp
  type(locreg_descriptors),intent(in) :: lr
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  real(wp), dimension((wfd%nvctr_c+7*wfd%nvctr_f)*orbs%nspinor*orbs%norbp), intent(in) :: psi
  real(wp), dimension((wfd%nvctr_c+7*wfd%nvctr_f)*orbs%nspinor*orbs%norbp), intent(inout) :: hpsi
  real(gp), intent(out) :: eproj_sum
  type(paw_objects),optional,intent(inout)::paw
  !local variables
  integer :: iat,nwarnings,iproj,iorb
  integer :: iatype
  integer :: istart_c,idir,isorb,ieorb,ikpt,nspinor,ispsi_k,ispsi

  call f_routine(id='applyprojectorsonthefly')
  
  !put idir=0, no derivative
  idir=0
  nwarnings=0
  eproj_sum=0.0_gp

  !quick return if no orbitals on this processor
  if (orbs%norbp == 0) then
     return
  end if

  !apply the projectors on the fly for each k-point of the processor
  !starting k-point
  ikpt=orbs%iokpt(1)

  ispsi_k=1
  loop_kpt: do

     call orbs_in_kpt(ikpt,orbs,isorb,ieorb,nspinor)

     !this may not work for non-collinear cases
     iproj=0
     do iat=1,at%astruct%nat
        iatype=at%astruct%iatype(iat)
        istart_c=1

        call atom_projector(nlpsp, iatype, iat, at%astruct%atomnames(iatype), &
             & at%astruct%geocode, idir, lr, hx, hy, hz, &
             & orbs%kpts(1,ikpt), orbs%kpts(2,ikpt), orbs%kpts(3,ikpt), &
             & istart_c, iproj, nwarnings)

        !apply the projector to all the orbitals belonging to the processor
        ispsi=ispsi_k
        do iorb=isorb,ieorb
           istart_c=1
           if(at%npspcode(iatype) == PSPCODE_PAW) then
           !    PAW case:
              call apply_atproj_iorb_paw(iat,iorb,istart_c,&
                   at,orbs,wfd,nlpsp,&
                   psi(ispsi),hpsi(ispsi),paw%spsi(ispsi),eproj_sum,paw)
           else
           !    HGH or GTH case:
              call apply_atproj_iorb_new(iat,iorb,istart_c,nlpsp%nprojel,&
                   at,orbs,wfd,nlpsp%pspd(iat)%plr,nlpsp%proj,&
                   psi(ispsi),hpsi(ispsi),eproj_sum)
           end if
           ispsi=ispsi+(wfd%nvctr_c+7*wfd%nvctr_f)*nspinor
        end do
     end do
     if (iproj /= nlpsp%nproj) then
        call yaml_warning('Incorrect number of projectors created')
        stop
     end if
     if (ieorb == orbs%norbp) exit loop_kpt
     ikpt=ikpt+1
     ispsi_k=ispsi
  end do loop_kpt

  if (nwarnings /= 0 .and. iproc == 0 .and. nlpsp%nproj /=0 .and. idir == 0) then
     call yaml_map('Calculating wavelets expansion of projectors, found warnings',nwarnings,fmt='(i0)')
     if (nwarnings /= 0) then
        call yaml_newline()
        call yaml_warning('Projectors too rough: Consider modifying hgrid and/or the localisation radii.')
        !write(*,'(1x,a,i0,a)') 'found ',nwarnings,' warnings.'
        !write(*,'(1x,a)') 'Some projectors may be too rough.'
        !write(*,'(1x,a,f6.3)') 'Consider the possibility of modifying hgrid and/or the localisation radii.'
     end if
  end if

  call f_release_routine()

END SUBROUTINE applyprojectorsonthefly


!!$!>   Applies the projector associated on a given atom on a corresponding orbital
!!$subroutine apply_atproj_iorb(iat,iorb,istart_c,at,orbs,wfd,nlpspd,proj,psi,hpsi,eproj)
!!$  use module_base
!!$  use module_types
!!$  implicit none
!!$  integer, intent(in) :: iat,iorb
!!$  type(atoms_data), intent(in) :: at
!!$  type(orbitals_data), intent(in) :: orbs
!!$  type(wavefunctions_descriptors), intent(in) :: wfd
!!$  type(nonlocal_psp_descriptors), intent(in) :: nlpspd
!!$  real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
!!$  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,orbs%nspinor), intent(in) :: psi
!!$  integer, intent(inout) :: istart_c !< address of the starting point of the projector in proj array
!!$  real(gp), intent(inout) :: eproj
!!$  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,orbs%nspinor), intent(inout) :: hpsi
!!$  !Local variables
!!$  integer :: ispinor,ityp,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,jseg_c,l,i,istart_c_i,ncplx
!!$  real(gp) :: eproj_spinor
!!$
!!$  !complex functions or not
!!$  !this should be decided as a function of the orbital
!!$  !features of the k-point ikpt
!!$  call ncplx_kpt(orbs%iokpt(iorb),orbs,ncplx)
!!$
!!$  istart_c_i=istart_c
!!$  do ispinor=1,orbs%nspinor,ncplx
!!$     eproj_spinor=0.0_gp
!!$     if (ispinor >= 2) istart_c=istart_c_i
!!$     ityp=at%astruct%iatype(iat)
!!$     mbvctr_c=nlpspd%nvctr_p(2*iat-1)-nlpspd%nvctr_p(2*iat-2)
!!$     mbvctr_f=nlpspd%nvctr_p(2*iat  )-nlpspd%nvctr_p(2*iat-1)
!!$     
!!$     mbseg_c=nlpspd%nseg_p(2*iat-1)-nlpspd%nseg_p(2*iat-2)
!!$     mbseg_f=nlpspd%nseg_p(2*iat  )-nlpspd%nseg_p(2*iat-1)
!!$     jseg_c=nlpspd%nseg_p(2*iat-2)+1
!!$     !GTH and HGH pseudopotentials
!!$     do l=1,4
!!$        do i=1,3
!!$           if (at%psppar(l,i,ityp) /= 0.0_gp) then
!!$              call applyprojector(ncplx,l,i,at%psppar(0,0,ityp),at%npspcode(ityp),&
!!$                   wfd%nvctr_c,wfd%nvctr_f,wfd%nseg_c,wfd%nseg_f,wfd%keyv,wfd%keyg,&
!!$                   mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
!!$                   nlpspd%keyv_p(jseg_c),nlpspd%keyg_p(1,jseg_c),proj(istart_c),&
!!$                   psi(1,ispinor),hpsi(1,ispinor),eproj_spinor)
!!$              istart_c=istart_c+(mbvctr_c+7*mbvctr_f)*(2*l-1)*ncplx
!!$           end if
!!$        enddo
!!$     enddo
!!$     eproj=eproj+&
!!$          orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*eproj_spinor
!!$  end do
!!$END SUBROUTINE apply_atproj_iorb


!> Build the Hifj matrix for PSP
subroutine build_hgh_hij_matrix(npspcode,psppar,hij)
  use module_base, only: gp
  use psp_projectors, only: PSPCODE_GTH, PSPCODE_HGH, PSPCODE_HGH_K, PSPCODE_HGH_K_NLCC, PSPCODE_PAW
  implicit none
  !Arguments
  integer, intent(in) :: npspcode
  real(gp), dimension(0:4,0:6), intent(in) :: psppar
  real(gp), dimension(3,3,4), intent(out) :: hij
  !Local variables
  integer :: l,i,j
  real(gp), dimension(2,2,3) :: offdiagarr

  !enter the coefficients for the off-diagonal terms (HGH case, npspcode=PSPCODE_HGH)
  offdiagarr(1,1,1)=-0.5_gp*sqrt(3._gp/5._gp)
  offdiagarr(2,1,1)=-0.5_gp*sqrt(100._gp/63._gp)
  offdiagarr(1,2,1)=0.5_gp*sqrt(5._gp/21._gp)
  offdiagarr(2,2,1)=0.0_gp !never used
  offdiagarr(1,1,2)=-0.5_gp*sqrt(5._gp/7._gp)  
  offdiagarr(2,1,2)=-7._gp/3._gp*sqrt(1._gp/11._gp)
  offdiagarr(1,2,2)=1._gp/6._gp*sqrt(35._gp/11._gp)
  offdiagarr(2,2,2)=0.0_gp !never used
  offdiagarr(1,1,3)=-0.5_gp*sqrt(7._gp/9._gp)
  offdiagarr(2,1,3)=-9._gp*sqrt(1._gp/143._gp)
  offdiagarr(1,2,3)=0.5_gp*sqrt(63._gp/143._gp)
  offdiagarr(2,2,3)=0.0_gp !never used

!  call to_zero(3*3*4,hij(1,1,1))
  hij=0.0_gp

  do l=1,4
     !term for all npspcodes
     loop_diag: do i=1,3
        hij(i,i,l)=psppar(l,i) !diagonal term
        if ((npspcode == PSPCODE_HGH .and. l/=4 .and. i/=3) .or. &
             ((npspcode == PSPCODE_HGH_K .or. npspcode == PSPCODE_HGH_K_NLCC) .and. i/=3)) then !HGH(-K) case, offdiagonal terms
           loop_offdiag: do j=i+1,3
              if (psppar(l,j) == 0.0_gp) exit loop_offdiag
              !offdiagonal HGH term
              if (npspcode == PSPCODE_HGH) then !traditional HGH convention
                 hij(i,j,l)=offdiagarr(i,j-i,l)*psppar(l,j)
              else !HGH-K convention
                 hij(i,j,l)=psppar(l,i+j+1)
              end if
              hij(j,i,l)=hij(i,j,l) !symmetrization
           end do loop_offdiag
        end if
     end do loop_diag
  end do
  
end subroutine build_hgh_hij_matrix


!> Apply the PSP projectors
subroutine applyprojector(ncplx,l,i,psppar,npspcode,&
     nvctr_c,nvctr_f,nseg_c,nseg_f,keyv,keyg,&
     mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p,keyg_p,proj,psi,hpsi,eproj)
  use module_base, only: gp,wp,dp
  use psp_projectors, only: PSPCODE_GTH, PSPCODE_HGH, PSPCODE_HGH_K, PSPCODE_HGH_K_NLCC, PSPCODE_PAW
  implicit none
  integer, intent(in) :: i,l,npspcode,ncplx
  integer, intent(in) :: nvctr_c,nvctr_f,nseg_c,nseg_f,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  integer, dimension(mbseg_c+mbseg_f), intent(in) :: keyv_p
  integer, dimension(2,mbseg_c+mbseg_f), intent(in) :: keyg_p
  real(wp), dimension(*), intent(in) :: proj
  real(gp), dimension(0:4,0:6), intent(in) :: psppar
  real(wp), dimension(nvctr_c+7*nvctr_f,ncplx), intent(in) :: psi
  real(gp), intent(inout) :: eproj
  real(wp), dimension(nvctr_c+7*nvctr_f,ncplx), intent(inout) :: hpsi
  !local variables
  integer :: j,m,istart_c,istart_c_i,istart_c_j,icplx
  real(dp), dimension(2) :: scpr,scprp,scpr_i,scprp_i,scpr_j,scprp_j
  real(gp), dimension(2,2,3) :: offdiagarr
  real(gp) :: hij

  !enter the coefficients for the off-diagonal terms (HGH case, npspcode=3)
  offdiagarr(1,1,1)=-0.5_gp*sqrt(3._gp/5._gp)
  offdiagarr(2,1,1)=-0.5_gp*sqrt(100._gp/63._gp)
  offdiagarr(1,2,1)=0.5_gp*sqrt(5._gp/21._gp)
  offdiagarr(2,2,1)=0.0_gp !never used
  offdiagarr(1,1,2)=-0.5_gp*sqrt(5._gp/7._gp)  
  offdiagarr(2,1,2)=-7._gp/3._gp*sqrt(1._gp/11._gp)
  offdiagarr(1,2,2)=1._gp/6._gp*sqrt(35._gp/11._gp)
  offdiagarr(2,2,2)=0.0_gp !never used
  offdiagarr(1,1,3)=-0.5_gp*sqrt(7._gp/9._gp)
  offdiagarr(2,1,3)=-9._gp*sqrt(1._gp/143._gp)
  offdiagarr(1,2,3)=0.5_gp*sqrt(63._gp/143._gp)
  offdiagarr(2,2,3)=0.0_gp !never used

  istart_c=1
  !start of the routine for projectors application
  do m=1,2*l-1

     call wpdot_wrap(ncplx,  &
          nvctr_c,nvctr_f,nseg_c,nseg_f,keyv,keyg,psi,  &
          mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p,keyg_p,proj(istart_c),scpr)
  
     do icplx=1,ncplx
        scprp(icplx)=scpr(icplx)*real(psppar(l,i),dp)
        eproj=eproj+real(scprp(icplx),gp)*real(scpr(icplx),gp)
     end do

     call waxpy_wrap(ncplx,scprp,&
          mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p,keyg_p,proj(istart_c),&
          nvctr_c,nvctr_f,nseg_c,nseg_f,keyv,keyg,hpsi)

     !print *,'scprp,m,l,i',scprp,m,l,i

     istart_c=istart_c+(mbvctr_c+7*mbvctr_f)*ncplx
  enddo
  if ((npspcode == PSPCODE_HGH .and. l/=4 .and. i/=3) .or. &
       ((npspcode == PSPCODE_HGH_K .or. npspcode == PSPCODE_HGH_K_NLCC ).and. i/=3)) then !HGH(-K) case, offdiagonal terms
     loop_j: do j=i+1,3
        if (psppar(l,j) == 0.0_gp) exit loop_j

        !offdiagonal HGH term
        if (npspcode == PSPCODE_HGH) then !traditional HGH convention
           hij=offdiagarr(i,j-i,l)*psppar(l,j)
        else !HGH-K convention
           hij=psppar(l,i+j+1)
        end if

        !starting addresses of the projectors
        istart_c_i=istart_c-(2*l-1)*(mbvctr_c+7*mbvctr_f)*ncplx
        istart_c_j=istart_c_i+(j-i)*(2*l-1)*(mbvctr_c+7*mbvctr_f)*ncplx
        do m=1,2*l-1
           call wpdot_wrap(ncplx,nvctr_c,nvctr_f,nseg_c,nseg_f,keyv,keyg,psi,  &
                mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p,keyg_p,proj(istart_c_j),scpr_j)

           call wpdot_wrap(ncplx,nvctr_c,nvctr_f,nseg_c,nseg_f,keyv,keyg,psi,  &
                mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p,keyg_p,proj(istart_c_i),scpr_i)

           do icplx=1,ncplx
              scprp_j(icplx)=scpr_j(icplx)*hij
              scprp_i(icplx)=scpr_i(icplx)*hij
              !scpr_i*h_ij*scpr_j+scpr_j*h_ij*scpr_i
              eproj=eproj+2._gp*hij*real(scpr_j(icplx),gp)*real(scpr_i(icplx),gp)
           end do

           !|hpsi>=|hpsi>+h_ij (<p_i|psi>|p_j>+<p_j|psi>|p_i>)
           call waxpy_wrap(ncplx,scprp_j,&
                mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p,keyg_p,&
                proj(istart_c_i),&
                nvctr_c,nvctr_f,nseg_c,nseg_f,keyv,keyg,hpsi)

           call waxpy_wrap(ncplx,scprp_i,&
                mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
                keyv_p,keyg_p,proj(istart_c_j),&
                nvctr_c,nvctr_f,nseg_c,nseg_f,keyv,keyg,hpsi)

           istart_c_j=istart_c_j+(mbvctr_c+7*mbvctr_f)*ncplx
           istart_c_i=istart_c_i+(mbvctr_c+7*mbvctr_f)*ncplx
        enddo
     end do loop_j
  end if
END SUBROUTINE applyprojector


!sij_opt=1 : obtain hpsi
!sij_opt=2 : obtain spsi
!sij_opt=3 : obtain hpsi and spsi
!
subroutine applyprojector_paw(ncplx,istart_c,iat,&
     nvctr_c,nvctr_f,nseg_c,nseg_f,keyv,keyg,&
     mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p,keyg_p,proj,&
     psi,hpsi,spsi,eproj,proj_G,paw_ij,&
     lmn2_size,cprj_out,sij_opt,sij)
  use module_base
  use module_types
  use gaussians, only:gaussian_basis_new, gaussian_basis_iter, gaussian_iter_start, gaussian_iter_next_shell
  use m_paw_ij, only: paw_ij_type
  use m_pawcprj, only: pawcprj_type
  implicit none
  integer,parameter::nspinor=1  !not yet implemented
  integer, intent(inout)::istart_c
  integer, intent(in) :: ncplx,lmn2_size,sij_opt,iat
  integer, intent(in) :: nvctr_c,nvctr_f,nseg_c,nseg_f,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  integer, dimension(mbseg_c+mbseg_f), intent(in) :: keyv_p
  integer, dimension(2,mbseg_c+mbseg_f), intent(in) :: keyg_p
  real(wp), dimension(*), intent(in) :: proj
  real(wp), dimension(nvctr_c+7*nvctr_f,ncplx), intent(in) :: psi
  type(gaussian_basis_new),intent(in)::proj_G
  type(paw_ij_type),intent(in)::paw_ij
  !type(cprj_objects),dimension(1,nspinor),intent(out)::cprj_out
  type(pawcprj_type),intent(inout)::cprj_out
  real(gp), intent(out) :: eproj
  real(wp), dimension(nvctr_c+7*nvctr_f,ncplx), intent(inout) :: hpsi
  real(wp), dimension(nvctr_c+7*nvctr_f,ncplx), intent(inout) :: spsi
  real(wp), dimension(lmn2_size),intent(in)::sij
  !local variables
  character(len=*),parameter::subname='applyprojector_paw'
  integer :: ilmn,jlmn,klmn,j0lmn,ispinor
  integer :: klmnc,i_m,j_m,iaux
  integer :: istart_j,icplx
  real(gp)::eproj_i
  real(gp)::ddot
  real(dp), dimension(2) :: scpr
  real(gp) :: dij
  !real(wp), dimension(:,:), allocatable :: cprj_i
  real(wp), dimension(:,:), allocatable :: cprj,dprj !scalar products with the projectors (always assumed to be complex and spinorial)
  integer :: proj_count
  type(gaussian_basis_iter) :: iter, iter2, iter0

! change: keyv_p by nlpspd%keyv_p(jseg_c),&

!
  proj_count= paw_ij%lmn_size
  cprj = f_malloc((/ nspinor*ncplx, proj_count /),id='cprj')
  dprj = f_malloc((/ nspinor*ncplx, proj_count /),id='dprj')

  !cprj_out(1,1:nspinor)%cp(1:ncplx,1:proj_count)=0.0_wp
  eproj=0.0_gp
  cprj=0.0_wp
  !call to_zero(4*7*3*4,cprj(1,1,1,1))

  !Use special subroutines for these number of projectors
  !if (proj_count.eq.4 .or. proj_count.eq.5 .or. proj_count.eq.8 .or. proj_count.eq.13 &
  !    .or. proj_count.eq.14 .or. proj_count.eq.18 .or. proj_count.eq.19 &
  !    .or. proj_count.eq.20 .or. proj_count.eq.22) then

  !  allocate(cprj_i(proj_count,ncplx))
  !
  !  !loop over all the components of the wavefunction
  !  do ispinor=1,orbs%nspinor,ncplx
  !               call wpdot_wrap1(ncplx,  &
  !                    wfd%nvctr_c,wfd%nvctr_f,wfd%nseg_c,wfd%nseg_f,&
  !                    wfd%keyv,wfd%keyglob,&
  !                    psi(1,ispinor), &
  !                    mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
  !                    plr%wfd%keyv,&!nlpspd%keyv_p(jseg_c),&
  !                    plr%wfd%keyglob,&!nlpspd%keyg_p(1,jseg_c),&
  !                    proj(istart_c),&
  !                    cprj_i,proj_count)
  !
  !    i_proj=1
  !    do i_shell=1,proj_G%nshltot
  !      i_l=proj_G%nam(i_shell)
  !      do i_m=1,2*i_l-1
  !        do icplx=1,ncplx
  !         cprj(ispinor+icplx-1,m,i,l) = cprj_i(i_proj,icplx)
  !        enddo
  !        i_proj=i_proj+1
  !      end do
  !    end do
  !  end do

  !  deallocate(cprj_i)
  !else !use standart subroutine for projector application

  call gaussian_iter_start(proj_G, iat, iter0)

! Get cprj:
  jlmn=0
  istart_j=istart_c
  iter = iter0
  do
     if (.not. gaussian_iter_next_shell(proj_G, iter)) exit
     do j_m=1,2*iter%l-1
        jlmn=jlmn+1
        !loop over all the components of the wavefunction
        do ispinor=1,nspinor,ncplx
           call wpdot_wrap(ncplx,  &
                nvctr_c,nvctr_f,nseg_c,nseg_f,&
                keyv,keyg,&
                psi(1,ispinor), &
                mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
                keyv_p,keyg_p,&
                proj(istart_j),&
                cprj(ispinor,jlmn))
        end do !ispinor
        istart_j=istart_j+(mbvctr_c+7*mbvctr_f)*ncplx
     end do !i_m
  end do !j_shell
  !
  !copy cprj into cprj_out
  !
!  do ispinor=1,nspinor
!    do icplx=1,ncplx
!      jspinor=jspinor+1
!      cprj_out(1,ispinor)%cp(icplx,:)=cprj(jspinor,:) 
!    end do
!  end do
  do ispinor=1,nspinor
     cprj_out%cp(ispinor,:)=cprj(ispinor,:)
     !write(*,*)'applyprojector_paw: erase me: l1212 cprj=',cprj_out%cp(ispinor,:)
  end do
!

  if(sij_opt==1 .or. sij_opt==3) then
  !CALCULATE |H|PSI>
     dprj=0.0_wp
     iaux=paw_ij%cplex_dij*paw_ij%lmn2_size
     !call calculate_dprj(paw_ij%dij,iaux,paw_ij%ndij)
     call calculate_dprj(paw_ij%dij(:,1),iaux)
     !
     !apply non-local operator
     istart_j=istart_c
     call apply_non_local_operator(hpsi,nvctr_c+7*nvctr_f,ncplx,istart_j)
     eproj=eproj+eproj_i
     !
     !DEBUG: calculate <PSI|H|PSI>, only for 1 orbital and ncplx=1
     do ispinor=1,nspinor,ncplx
        scpr(1)=ddot(nvctr_c+7*nvctr_f,psi(istart_c,ispinor),1,hpsi(istart_c,ispinor),1)
        !write(*,*)'erase me: applyprojector_paw l1231'
        !write(*,*)'<psi|H|psi>= ',scpr(1:ncplx)
     end do !ispinor
  end if
  if(sij_opt==2 .or. sij_opt==3) then
  !CALCULATE |S-I|PSI>
     dprj=0.0_wp
     !Pending: check if it works  for cplex_dij=2
     iaux=paw_ij%cplex_dij*paw_ij%lmn2_size
     !DEBUG
     !write(*,*)'erase me, applyprojector_paw, l1241 sij=',sij(1:iaux)
     !END DEBUG
     call calculate_dprj(sij(1:iaux),iaux)
     !
     !apply non-local operator
     istart_j=istart_c
     call apply_non_local_operator(spsi,nvctr_c+7*nvctr_f,ncplx,istart_j)
     !
     !DEBUG: calculate <PSI|S-I|PSI>, only for 1 orbital and ncplx=1
     !do ispinor=1,nspinor,ncplx
     !   scpr(1)=ddot(nvctr_c+7*nvctr_f,psi(istart_c,ispinor),1,spsi(istart_c,ispinor),1)
     !   write(*,*)'erase me: applyprojector_paw l1260'
     !   write(*,*)'<psi|S-I|psi>= ',scpr(1:ncplx)
     !end do !ispinor
  end if

  !update istart_c, note that we only used istart_j above.
  istart_c=istart_j

  call f_free(cprj)
  call f_free(dprj)

  contains

!real(8) function ddot(n,A,l1,B,l2)
!  implicit none
!  integer, intent(in)::n,l1,l2
!  real(8),intent(in),dimension(n)::A,B
!  real(8)::scpr
!  integer::i

!  scpr=0.00_dp
!  do i=1,n
!   ddot=ddot+A(i)*B(i)
!  end do
!end function ddot

  subroutine calculate_dprj(kij,dim1)
     !Here we calculate:
     !dprj(i)= sum_{j} dij <p_j|psi>
  
     implicit none
     integer,intent(in)::dim1
     real(wp),dimension(dim1),intent(in)::kij 

     !apply the matrix of the coefficients on the cprj array
     jlmn=0
     iter = iter0
     do
        if (.not. gaussian_iter_next_shell(proj_G, iter)) exit
        do j_m=1,2*iter%l-1
           jlmn=jlmn+1
           j0lmn=jlmn*(jlmn-1)/2
           !Diagonal components
           klmn=j0lmn+jlmn;klmnc=paw_ij%cplex_dij*(klmn-1)
           !case of cplex_dij pending
           !dij=paw_ij%dij(klmn,1)
           dij=kij(klmn)
           !write(*,'(2(i2,x),f18.10,i2)')j_shell,j_m,dij,klmn
           do ispinor=1,nspinor !real matrix
              dprj(ispinor,jlmn)=dprj(ispinor,jlmn)+&
              dij*cprj(ispinor,jlmn)
           end do
           !Off-diagonal components
           ilmn=0
           iter2 = iter0
           do
              if (.not. gaussian_iter_next_shell(proj_G, iter2)) exit
              do i_m=1,2*iter2%l-1
                 if(i_m>=j_m .and. iter%ishell==iter2%ishell) cycle
                 ilmn=ilmn+1
                 klmn=j0lmn+ilmn;klmnc=paw_ij%cplex_dij*(klmn-1)
                 dij=paw_ij%dij(klmn,1)
                 !write(*,'(4(i2,x),f18.10,i2)')j_shell,j_m,i_shell,i_m,dij,klmn
                 do ispinor=1,nspinor !real matrix
                     dprj(ispinor,jlmn)=dprj(ispinor,jlmn)+&
                         dij*cprj(ispinor,ilmn)
                     dprj(ispinor,ilmn)=dprj(ispinor,ilmn)+&
                         dij*cprj(ispinor,jlmn)
                 end do
              end do
           end do
        end do
     end do
  end subroutine calculate_dprj 

  subroutine apply_non_local_operator(apham,dim1,dim2,istart)

     implicit none
     integer,intent(in)::dim1,dim2
     integer,intent(inout)::istart
     real(wp),dimension(dim1,dim2), intent(inout) :: apham
     
 
     !build a single array via daxpy for the projectors
     !apply the non-local operator on the wavefunction
     !for the moment use the traditional waxpy instead of daxpy, for test purposes
     eproj_i=0.0_gp
     jlmn=0
     iter = iter0
     do
        if (.not. gaussian_iter_next_shell(proj_G, iter)) exit
        do j_m=1,2*iter%l-1
           jlmn=jlmn+1
           do ispinor=1,nspinor,ncplx
              do icplx=1,ncplx
                 eproj_i=eproj_i+dprj(ispinor+icplx-1,jlmn)*cprj(ispinor+icplx-1,jlmn)
              end do
              call waxpy_wrap(ncplx,dprj(ispinor,jlmn),&
                   mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
                   keyv_p,&!nlpspd%keyv_p(jseg_c),&
                   keyg_p,&!nlpspd%keyg_p(1,jseg_c),&
                   proj(istart),&
                   nvctr_c,nvctr_f,nseg_c,nseg_f,&
                   keyv,keyg,&
                   apham(1,ispinor))
           end do
           istart=istart+(mbvctr_c+7*mbvctr_f)*ncplx
        end do
     end do
  end subroutine apply_non_local_operator



end subroutine applyprojector_paw


!> Applies the projector associated on a given atom on a corresponding orbital
!! uses a generic representation of the projector to generalize the form of the projector  
subroutine apply_atproj_iorb_new(iat,iorb,istart_c,nprojel,at,orbs,wfd,&
     plr,proj,&
     psi,hpsi,eproj)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iat,iorb,nprojel
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(locreg_descriptors), intent(in) :: plr
  !type(nonlocal_psp_descriptors), intent(in) :: nlpspd
  real(wp), dimension(nprojel), intent(in) :: proj
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,orbs%nspinor), intent(in) :: psi
  integer, intent(inout) :: istart_c !< address of the starting point of the projector in proj array
  real(gp), intent(inout) :: eproj
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,orbs%nspinor), intent(inout) :: hpsi
  !local variables
  !character(len=*), parameter :: subname='apply_atproj_iorb'
  integer :: ispinor,ityp,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,l,i,istart_c_i,ncplx,m,j,icplx
  real(gp) :: eproj_i
  real(wp), dimension(4,7,3,4) :: cproj,dproj !<scalar products with the projectors (always assumed to be complex and spinorial)
  real(gp), dimension(3,3,4) :: hij_hgh 
!!$  integer :: jseg_c
!!$  real(wp), dimension(:,:), allocatable :: wproj !work array for the application of the projectors
  real(wp), dimension(:,:), allocatable :: cproj_i
  integer :: proj_count, i_proj
 

  !parameter for the descriptors of the projectors
  ityp=at%astruct%iatype(iat)

  call plr_segs_and_vctrs(plr,mbseg_c,mbseg_f,mbvctr_c,mbvctr_f)
 
  !complex functions or not
  !this should be decided as a function of the orbital
  !features of the k-point ikpt
  call ncplx_kpt(orbs%iokpt(iorb),orbs,ncplx)

  !build the matrix of the pseudopotential
  call build_hgh_hij_matrix(at%npspcode(ityp),at%psppar(0,0,ityp),hij_hgh)

!!$  allocate(wproj(mbvctr_c+7*mbvctr_f,ncplx+ndebug),stat=i_stat)
!!$  call memocc(i_stat,wproj,'wproj',subname)

  !calculate the scalar product with all the projectors of the atom
  !call to_zero(4*7*3*4,cproj(1,1,1,1))
  cproj=0.0_wp

  proj_count = 0
  !count over all the channels
  do l=1,4
     !loop over all the projectors of the channel
     do i=1,3
        !loop over all the components of the projector
        if (at%psppar(l,i,ityp) /= 0.0_gp) then
           do m=1,2*l-1
              proj_count=proj_count+1
           end do
        end if
     end do
  end do

  !Use special subroutines for these number of projectors
!  if (proj_count.eq.4 .or. proj_count.eq.5 .or. proj_count.eq.8 .or. proj_count.eq.13 &
!      .or. proj_count.eq.14 .or. proj_count.eq.18 .or. proj_count.eq.19 &
!      .or. proj_count.eq.20 .or. proj_count.eq.22) then
     !more elegant?
  if (any(proj_count==(/4,5,8,13,14,18,19,20,22/))) then

    cproj_i = f_malloc((/ proj_count, ncplx /),id='cproj_i')

    !loop over all the components of the wavefunction
    do ispinor=1,orbs%nspinor,ncplx
                 call wpdot_wrap1(ncplx,  &
                      wfd%nvctr_c,wfd%nvctr_f,wfd%nseg_c,wfd%nseg_f,&
                      wfd%keyvglob,wfd%keyglob,&
                      psi(1,ispinor), &
                      mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
                      plr%wfd%keyvglob,&
                      plr%wfd%keyglob,&
                      proj(istart_c),&
                      cproj_i,proj_count)

      i_proj=1
      do l=1,4
       !loop over all the projectors of the channel
       do i=1,3
        !loop over all the components of the projector
        if (at%psppar(l,i,ityp) /= 0.0_gp) then
           do m=1,2*l-1
             do icplx=1,ncplx
              cproj(ispinor+icplx-1,m,i,l) = cproj_i(i_proj,icplx)
             enddo
              i_proj=i_proj+1
           end do
        end if
       end do
      end do

    end do

    call f_free(cproj_i)

    !print *,'iorb,cproj',iorb,sum(cproj)

  else  ! use standard subroutine for projector application

    !index for performing the calculation with all the projectors
    istart_c_i=istart_c
    !loop over all the channels (from s to f)
    do l=1,4
       !loop over all the projectors of the channel
       do i=1,3
          !loop over all the components of the projector
          if (at%psppar(l,i,ityp) /= 0.0_gp) then
             do m=1,2*l-1
              !loop over all the components of the wavefunction
                do ispinor=1,orbs%nspinor,ncplx
                   call wpdot_wrap(ncplx,  &
                        wfd%nvctr_c,wfd%nvctr_f,wfd%nseg_c,wfd%nseg_f,&
                        wfd%keyvglob,wfd%keyglob,&
                        psi(1,ispinor), &
                        mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
                        plr%wfd%keyvglob,&!nlpspd%keyv_p(jseg_c),&
                        plr%wfd%keyglob,&!nlpspd%keyg_p(1,jseg_c),&
                        proj(istart_c_i),&
                        cproj(ispinor,m,i,l))
                end do
                istart_c_i=istart_c_i+(mbvctr_c+7*mbvctr_f)*ncplx
             end do
          end if
       end do
    end do

  endif

  !apply the matrix of the coefficients on the cproj array
  !call to_zero(4*7*3*4,dproj(1,1,1,1))
  dproj=0.0_wp
  do l=1,4 !diagonal in l
     do i=1,3
        do j=1,3
           do m=1,2*l-1 !diagonal in m
              do ispinor=1,orbs%nspinor !real matrix
                 dproj(ispinor,m,i,l)=dproj(ispinor,m,i,l)+&
                     hij_hgh(i,j,l)*cproj(ispinor,m,j,l)
              end do
           end do
        end do
     end do
  end do
  istart_c_i=istart_c
  !build a single array via daxpy for the projectors
  !apply the non-local operator on the wavefunction
  !for the moment use the traditional waxpy instead of daxpy, for test purposes
  eproj_i=0.0_gp
!!$  call nanosec(itsc0) 
  do l=1,4
     !loop over all the projectors of the channel
     do i=1,3
        !loop over all the components of the projector
        if (at%psppar(l,i,ityp) /= 0.0_gp) then
           do m=1,2*l-1
              !loop over all the components of the wavefunction
              do ispinor=1,orbs%nspinor,ncplx
                 
                 do icplx=1,ncplx
                    eproj_i=eproj_i+dproj(ispinor+icplx-1,m,i,l)*cproj(ispinor+icplx-1,m,i,l)
                 end do

                 call waxpy_wrap(ncplx,dproj(ispinor,m,i,l),&
                      mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
                      plr%wfd%keyvglob,&!nlpspd%keyv_p(jseg_c),&
                      plr%wfd%keyglob,&!nlpspd%keyg_p(1,jseg_c),&
                      proj(istart_c),&
                      wfd%nvctr_c,wfd%nvctr_f,wfd%nseg_c,wfd%nseg_f,&
                      wfd%keyvglob,wfd%keyglob,&
                      hpsi(1,ispinor))
              end do
              istart_c=istart_c+(mbvctr_c+7*mbvctr_f)*ncplx
           end do
        end if
     end do
  end do
!!$  call nanosec(itsc1)
!!$  print *,'normal time',real(itsc1-itsc0,gp)*1.e-9_gp

  eproj=eproj+&
       orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*eproj_i
!!$  istart_c=istart_c_i
!!$  call nanosec(itsc0) 
!!$  do l=1,4
!!$     !loop over all the projectors of the channel
!!$     do i=1,3
!!$        !loop over all the components of the projector
!!$        if (at%psppar(l,i,ityp) /= 0.0_gp) then
!!$           do m=1,2*l-1
!!$              !loop over all the components of the wavefunction
!!$              do ispinor=1,orbs%nspinor,ncplx
!!$                 
!!$                 do icplx=1,ncplx
!!$                    eproj_i=eproj_i+dproj(ispinor+icplx-1,m,i,l)*cproj(ispinor+icplx-1,m,i,l)
!!$                 end do
!!$
!!$                 call axpy((mbvctr_c+7*mbvctr_f)*ncplx,1.0_gp,proj(istart_c),1,wproj(1,1),1)
!!$
!!$              end do
!!$              istart_c=istart_c+(mbvctr_c+7*mbvctr_f)*ncplx
!!$           end do
!!$        end if
!!$     end do
!!$  end do
!!$  call nanosec(itsc1)
!!$  print *,'daxpy time',real(itsc1-itsc0,gp)*1.e-9_gp
!!$
!!$  i_all=-product(shape(wproj))*kind(wproj)
!!$  deallocate(wproj,stat=i_stat)
!!$  call memocc(i_stat,i_all,'wproj',subname)

END SUBROUTINE apply_atproj_iorb_new


!> Applies the projector associated on a given atom on a corresponding orbital
!! uses a generic representation of the projector to generalize the form of the projector  
subroutine apply_atproj_iorb_paw(iat,iorb,istart_c,at,orbs,wfd,&
     nlpsp,psi,hpsi,spsi,eproj,paw)
  use module_base
  use module_types
  use gaussians, only: gaussian_basis
  implicit none
  integer, intent(in) :: iat,iorb
  integer, intent(inout)::istart_c
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(DFT_PSP_projectors), intent(in) :: nlpsp
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,orbs%nspinor), intent(in) :: psi
  real(gp), intent(inout) :: eproj
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,orbs%nspinor), intent(inout) :: hpsi
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,orbs%nspinor), intent(inout) :: spsi
  type(paw_objects),intent(inout)::paw
  !local variables
  !character(len=*), parameter :: subname='apply_atproj_iorb'
  integer :: sij_opt
  integer :: ncplx
  integer :: ityp,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f
  real(gp) :: eproj_i

  !parameter for the descriptors of the projectors
  ityp=at%astruct%iatype(iat)

  call plr_segs_and_vctrs(nlpsp%pspd(iat)%plr,mbseg_c,mbseg_f,mbvctr_c,mbvctr_f)
 
  !complex functions or not
  !this should be decided as a function of the orbital
  !features of the k-point ikpt
  call ncplx_kpt(orbs%iokpt(iorb),orbs,ncplx)


  !calculate the scalar product with all the projectors of the atom
  !index for performing the calculation with all the projectors
  
  sij_opt=3 !get hpsi and spsi

   call applyprojector_paw(ncplx,istart_c,iat,&
        wfd%nvctr_c,wfd%nvctr_f,wfd%nseg_c,wfd%nseg_f,wfd%keyvglob,wfd%keyglob,&
        mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
        nlpsp%pspd(iat)%plr%wfd%keyvglob,& !nlpspd%keyv_p(jseg_c),
        nlpsp%pspd(iat)%plr%wfd%keyglob,& !nlpspd%keyg_p(1,jseg_c),&
        nlpsp%proj,&
        psi,hpsi,spsi,eproj_i,nlpsp%proj_G,paw%paw_ij(iat),&
        at%pawtab(ityp)%lmn2_size,paw%cprj(iat,iorb),&
        sij_opt,at%pawtab(ityp)%sij)  

  !DEBUG
  !do ii=1,wfd%nvctr_c+7*wfd%nvctr_f
  !   write(400,*)ii,hpsi(ii,1)
  !   write(401,*)ii,paw%spsi(ispsi+ii-1)
  !end do
  !DEBUG
  eproj=eproj+&
        &orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*eproj_i

end subroutine apply_atproj_iorb_paw


!> Find the starting and ending orbital for kpoint ikpt, and the corresponding nspinor
subroutine orbs_in_kpt(ikpt,orbs,isorb,ieorb,nspinor)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: ikpt
  type(orbitals_data), intent(in) :: orbs
  integer, intent(out) :: isorb,ieorb,nspinor
  !local variables
  integer :: iorb

  !disable starting and ending points for the case no orbitals on a given processor
  if (orbs%norbp == 0) then
     isorb=1
     ieorb=0
  end if

  !find starting orbital
  isorb=1 !default if orbs%norbp==0
  do iorb=1,orbs%norbp
     if (orbs%iokpt(iorb)==ikpt) then
        isorb=iorb
        exit
     end if
  end do

  !find ending orbital
  ieorb=0 !default if orbs%norbp==0
  do iorb=orbs%norbp,1,-1
     if (orbs%iokpt(iorb)==ikpt) then
        ieorb=iorb
        exit
     end if
  end do

  !nspinor for this k-point
  nspinor=orbs%nspinor

END SUBROUTINE orbs_in_kpt


!>   Determine whether the k-point is complex of real
!!   Find the starting and ending orbital for kpoint ikpt, and the corresponding nspinor
subroutine ncplx_kpt(ikpt,orbs,ncplx)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: ikpt
  type(orbitals_data), intent(in) :: orbs
  integer, intent(out) :: ncplx
  !local variables
  real(gp) :: kx,ky,kz

  !features of the k-point ikpt
  kx=orbs%kpts(1,ikpt)
  ky=orbs%kpts(2,ikpt)
  kz=orbs%kpts(3,ikpt)

  !evaluate the complexity of the k-point
  if (kx**2 + ky**2 + kz**2 == 0.0_gp) then
     ncplx=1
  else
     ncplx=2
  end if

END SUBROUTINE ncplx_kpt


!!!> Calculate the action of the local hamiltonian on the orbitals
!!subroutine local_hamiltonianParabola(iproc,orbs,lr,hx,hy,hz,&
!!     nspin,pot,psi,hpsi,ekin_sum,epot_sum, nat, rxyz, onWhichAtom, at)
!!  use module_base
!!  use module_types
!!  use module_interfaces
!!  use libxc_functionals
!!  implicit none
!!  integer, intent(in) :: iproc,nspin
!!  real(gp), intent(in) :: hx,hy,hz
!!  type(orbitals_data), intent(in) :: orbs
!!  type(locreg_descriptors), intent(in) :: lr
!!  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*orbs%norbp), intent(in) :: psi
!!  real(wp), dimension(*) :: pot
!!  !real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i*nspin) :: pot
!!  real(gp), intent(out) :: ekin_sum,epot_sum
!!  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*orbs%norbp), intent(out) :: hpsi
!!integer:: nat
!!real(8),dimension(3,nat):: rxyz 
!!integer,dimension(orbs%norbp):: onWhichAtom
!!type(atoms_data), intent(in) :: at
!!  !local variables
!!  character(len=*), parameter :: subname='local_hamiltonian'
!!  integer :: i_all,i_stat,iorb,npot,nsoffset,oidx,ispot
!!  real(wp) :: exctXcoeff
!!  real(gp) :: ekin,epot,kx,ky,kz,etest
!!  type(workarr_locham) :: wrk_lh
!!  real(wp), dimension(:,:), allocatable :: psir
!!real(8):: hxh, hyh, hzh
!!real(8),dimension(3):: rxyzShifted
!!
!!  exctXcoeff=libxc_functionals_exctXfac()
!!
!!  !initialise the work arrays
!!  call initialize_work_arrays_locham(lr,orbs%nspinor,wrk_lh)  
!!
!!  !components of the potential
!!  npot=orbs%nspinor
!!  if (orbs%nspinor == 2) npot=1
!!
!!  ! Wavefunction in real space
!!  allocate(psir(lr%d%n1i*lr%d%n2i*lr%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
!!  call memocc(i_stat,psir,'psir',subname)
!!
!!  call to_zero(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspinor,psir)
!!
!!  ekin_sum=0.0_gp
!!  epot_sum=0.0_gp
!!
!!  etest=0.0_gp
!!
!!  hxh=hx*.5d0
!!  hyh=hy*.5d0
!!  hzh=hz*.5d0
!!
!!  do iorb=1,orbs%norbp
!!
!!
!!     if(orbs%spinsgn(iorb+orbs%isorb)>0.0_gp .or. nspin == 1 .or. nspin == 4 ) then
!!        nsoffset=1
!!     else
!!        nsoffset=lr%d%n1i*lr%d%n2i*lr%d%n3i+1
!!     end if
!!
!!     oidx=(iorb-1)*orbs%nspinor+1
!!
!!     !transform the wavefunction in Daubechies basis to the wavefunction in ISF basis
!!     !the psir wavefunction is given in the spinorial form
!!     call daub_to_isf_locham(orbs%nspinor,lr,wrk_lh,psi(1,oidx),psir)
!!
!!     !ispot=1+lr%d%n1i*lr%d%n2i*lr%d%n3i*(nspin+iorb-1)
!!     !etest=etest+dot(lr%d%n1i*lr%d%n2i*lr%d%n3i,pot(ispot),1,psir(1,1),1)
!!     !print *,'epot, iorb,iproc,norbp',iproc,orbs%norbp,iorb,etest
!!
!!     !apply the potential to the psir wavefunction and calculate potential energy
!!     select case(lr%geocode)
!!     case('F')
!!
!!   ! ATTENTION: WITH SHIFTED PARABOLA
!!    rxyzShifted(1)=rxyz(1,onWhichAtom(iorb))+orbs%parabolaShift(1,iorb)
!!    rxyzShifted(2)=rxyz(2,onWhichAtom(iorb))+orbs%parabolaShift(2,iorb)
!!    rxyzShifted(3)=rxyz(3,onWhichAtom(iorb))+orbs%parabolaShift(3,iorb)
!!        call apply_potentialParabola(lr%d%n1,lr%d%n2,lr%d%n3,1,1,1,0,orbs%nspinor,npot,psir,&
!!             pot(nsoffset),epot, rxyzShifted, hxh, hyh, hzh, orbs%parabPrefacArr(at%astruct%iatype(onWhichAtom(iorb))), orbs%power, &
!!             lr%bounds%ibyyzz_r) !optional
!!   ! THIS WAS THE ORIGINAL
!!        !call apply_potentialParabola(lr%d%n1,lr%d%n2,lr%d%n3,1,1,1,0,orbs%nspinor,npot,psir,&
!!        !     pot(nsoffset),epot, rxyz(1,onWhichAtom(iorb)), hxh, hyh, hzh, orbs%parabPrefacArr(at%astruct%iatype(onWhichAtom(iorb))),  &
!!        !     lr%bounds%ibyyzz_r) !optional
!!
!!        !call apply_potentialParabola(lr%d%n1,lr%d%n2,lr%d%n3,1,1,1,0,orbs%nspinor,npot,psir,&
!!        !     pot(nsoffset),epot, rxyz(1,onWhichAtom(iorb)), hxh, hyh, hzh, orbs%parabPrefac,  &
!!        !     lr%bounds%ibyyzz_r) !optional
!!        !call apply_potential(lr%d%n1,lr%d%n2,lr%d%n3,1,1,1,0,orbs%nspinor,npot,psir,&
!!        !     pot(nsoffset),epot,&
!!        !     lr%bounds%ibyyzz_r) !optional
!!          
!!     case('P') 
!!        !here the hybrid BC act the same way
!!        call apply_potential(lr%d%n1,lr%d%n2,lr%d%n3,0,0,0,0,orbs%nspinor,npot,psir,&
!!             pot(nsoffset),epot)
!!
!!     case('S')
!!
!!        call apply_potential(lr%d%n1,lr%d%n2,lr%d%n3,0,1,0,0,orbs%nspinor,npot,psir,&
!!             pot(nsoffset),epot)
!!     end select
!!
!!     !k-point values, if present
!!     kx=orbs%kpts(1,orbs%iokpt(iorb))
!!     ky=orbs%kpts(2,orbs%iokpt(iorb))
!!     kz=orbs%kpts(3,orbs%iokpt(iorb))
!!
!!     if (exctXcoeff /= 0.0_gp) then
!!        ispot=1+lr%d%n1i*lr%d%n2i*lr%d%n3i*(nspin+iorb-1)
!!        !add to the psir function the part of the potential coming from the exact exchange
!!        call axpy(lr%d%n1i*lr%d%n2i*lr%d%n3i,exctXcoeff,pot(ispot),1,psir(1,1),1)
!!     end if
!!
!!     !apply the kinetic term, sum with the potential and transform back to Daubechies basis
!!     call isf_to_daub_kinetic(hx,hy,hz,kx,ky,kz,orbs%nspinor,lr,wrk_lh,&
!!          psir,hpsi(1,oidx),ekin)
!!
!!     ekin_sum=ekin_sum+orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*ekin
!!     epot_sum=epot_sum+orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*epot
!!
!!  enddo
!!
!!  !print *,'iproc,etest',etest
!!
!!  !deallocations of work arrays
!!  i_all=-product(shape(psir))*kind(psir)
!!  deallocate(psir,stat=i_stat)
!!  call memocc(i_stat,i_all,'psir',subname)
!!
!!  call deallocate_work_arrays_locham(lr,wrk_lh)
!!
!!END SUBROUTINE local_hamiltonianParabola
!!
!!
!!!> routine for applying the local potentials
!!!! supports the non-collinear case, the buffer for tails and different Boundary Conditions
!!!! Optimal also for the complex wavefuntion case
!!subroutine apply_potentialParabola(n1,n2,n3,nl1,nl2,nl3,nbuf,nspinor,npot,psir,pot,epot, rxyzParab, &
!!     hxh, hyh, hzh, parabPrefac, power, &
!!     ibyyzz_r) !optional
!!  use module_base
!!  implicit none
!!  integer, intent(in) :: n1,n2,n3,nl1,nl2,nl3,nbuf,nspinor,npot
!!  real(wp), dimension(-14*nl1:2*n1+1+15*nl1,-14*nl2:2*n2+1+15*nl2,-14*nl3:2*n3+1+15*nl3,nspinor), intent(inout) :: psir
!!  real(wp), dimension(-14*nl1:2*n1+1+15*nl1-4*nbuf,-14*nl2:2*n2+1+15*nl2-4*nbuf,&
!!       -14*nl3:2*n3+1+15*nl3-4*nbuf,npot), intent(in) :: pot
!!  integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r
!!  real(gp), intent(out) :: epot
!!real(8),dimension(3):: rxyzParab
!!real(8):: hxh, hyh, hzh, parabPrefac
!!integer:: power
!!  !local variables
!!  integer :: i1,i2,i3,i1s,i1e,ispinor
!!  real(wp) :: tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42,tt
!!  real(wp) :: psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4
!!  real(gp) :: epot_p
!!real(8):: hxh2, hyh2, hzh2
!!integer:: ix0, iy0, iz0, istat, ipot
!!real(8),dimension(:,:,:,:),allocatable:: potCopy
!!
!!!write(*,'(a,2i7)') '-14*nl3, 2*n3+1+15*nl3', -14*nl3, 2*n3+1+15*nl3
!!!write(*,'(a,2i7)') '-14*nl2, 2*n2+1+15*nl2', -14*nl2, 2*n2+1+15*nl2
!!!write(*,'(a,2i7)') 'max(ibyyzz_r(1,i2,i3)-14,-14+2*nbuf), min(ibyyzz_r(2,i2,i3)-14,2*n1+16-2*nbuf)', max(ibyyzz_r(1,i2,i3)-14,-14+2*nbuf), min(ibyyzz_r(2,i2,i3)-14,2*n1+16-2*nbuf)
!!allocate(potCopy(-14*nl1:2*n1+1+15*nl1-4*nbuf,-14*nl2:2*n2+1+15*nl2-4*nbuf,&
!!       -14*nl3:2*n3+1+15*nl3-4*nbuf,npot), stat=istat)
!!potCopy=0.d0
!!  
!!  !the Tail treatment is allowed only in the Free BC case
!!  if (nbuf /= 0 .and. nl1*nl2*nl3 == 0) stop 'NONSENSE: nbuf/=0 only for Free BC'
!!
!!  epot=0.0_wp
!!  ! Copy the potential
!!!write(*,*) 'size(pot)', size(pot)
!!  !potCopy=pot
!!  do ipot=1,npot
!!      do i3=-14*nl3,2*n3+1+15*nl3-4*nbuf
!!          do i2=-14*nl2,2*n2+1+15*nl2-4*nbuf
!!              do i1=-14*nl1,2*n1+1+15*nl1-4*nbuf
!! !write(1,'(a,4i8)') 'i1, i2, i3, ipot', i1, i2, i3, ipot
!!                  potCopy(i1,i2,i3,ipot)=pot(i1,i2,i3,ipot)
!!              end do
!!          end do
!!      end do
!!  end do
!!  !potCopy(-14*nl1:2*n1+1+15*nl1-4*nbuf,-14*nl2:2*n2+1+15*nl2-4*nbuf,&
!!  !     -14*nl3:2*n3+1+15*nl3-4*nbuf,1:npot) &
!!  !  =pot(-14*nl1:2*n1+1+15*nl1-4*nbuf,-14*nl2:2*n2+1+15*nl2-4*nbuf,&
!!  !     -14*nl3:2*n3+1+15*nl3-4*nbuf,1:npot)
!!
!!   ix0=nint(rxyzParab(1)/hxh)
!!   iy0=nint(rxyzParab(2)/hyh)
!!   iz0=nint(rxyzParab(3)/hzh)
!!   hxh2=hxh**2
!!   hyh2=hyh**2
!!   hzh2=hzh**2
!!
!!!$omp parallel default(private)&
!!!$omp shared(pot,psir,n1,n2,n3,epot,ibyyzz_r,nl1,nl2,nl3,nbuf,nspinor)
!!  !case without bounds
!!  i1s=-14*nl1
!!  i1e=2*n1+1+15*nl1
!!  epot_p=0._gp
!!!$omp do
!!  do i3=-14*nl3,2*n3+1+15*nl3
!!     if (i3 >= -14+2*nbuf .and. i3 <= 2*n3+16-2*nbuf) then !check for the nbuf case
!!        do i2=-14*nl2,2*n2+1+15*nl2
!!           if (i2 >= -14+2*nbuf .and. i2 <= 2*n2+16-2*nbuf) then !check for the nbuf case
!!              !this if statement is inserted here for avoiding code duplication
!!              !it is to be seen whether the code results to be too much unoptimised
!!              if (present(ibyyzz_r)) then
!!                 !in this case we are surely in Free BC
!!                 !the min is to avoid to calculate for no bounds
!!                 do i1=-14+2*nbuf,min(ibyyzz_r(1,i2,i3),ibyyzz_r(2,i2,i3))-14-1
!!                    psir(i1,i2,i3,:)=0.0_wp
!!                 enddo
!!                 i1s=max(ibyyzz_r(1,i2,i3)-14,-14+2*nbuf)
!!                 i1e=min(ibyyzz_r(2,i2,i3)-14,2*n1+16-2*nbuf)
!!              end if
!!              
!!              !here we put the branchments wrt to the spin
!!              if (nspinor == 4) then
!!                 do i1=i1s,i1e
!!                    !wavefunctions
!!                    psir1=psir(i1,i2,i3,1)
!!                    psir2=psir(i1,i2,i3,2)
!!                    psir3=psir(i1,i2,i3,3)
!!                    psir4=psir(i1,i2,i3,4)
!!                    !potentials
!!                    pot1=pot(i1-2*nbuf,i2-2*nbuf,i3-2*nbuf,1)
!!                    pot2=pot(i1-2*nbuf,i2-2*nbuf,i3-2*nbuf,2)
!!                    pot3=pot(i1-2*nbuf,i2-2*nbuf,i3-2*nbuf,3)
!!                    pot4=pot(i1-2*nbuf,i2-2*nbuf,i3-2*nbuf,4)
!!
!!                    !diagonal terms
!!                    tt11=pot1*psir1 !p1
!!                    tt22=pot1*psir2 !p2
!!                    tt33=pot4*psir3 !p3
!!                    tt44=pot4*psir4 !p4
!!                    !Rab*Rb
!!                    tt13=pot2*psir3 !p1
!!                    !Iab*Ib
!!                    tt14=pot3*psir4 !p1
!!                    !Rab*Ib
!!                    tt23=pot2*psir4 !p2
!!                    !Iab*Rb
!!                    tt24=pot3*psir3 !p2
!!                    !Rab*Ra
!!                    tt31=pot2*psir1 !p3
!!                    !Iab*Ia
!!                    tt32=pot3*psir2 !p3
!!                    !Rab*Ia
!!                    tt41=pot2*psir2 !p4
!!                    !Iab*Ra
!!                    tt42=pot3*psir1 !p4
!!
!!                    ! Change epot later
!!                    epot_p=epot_p+tt11*psir1+tt22*psir2+tt33*psir3+tt44*psir4+&
!!                         2.0_gp*tt31*psir3-2.0_gp*tt42*psir4+2.0_gp*tt41*psir4+2.0_gp*tt32*psir3
!!
!!                    !wavefunction update
!!                    !p1=h1p1+h2p3-h3p4
!!                    !p2=h1p2+h2p4+h3p3
!!                    !p3=h2p1+h3p2+h4p3
!!                    !p4=h2p2-h3p1+h4p4
!!                    psir(i1,i2,i3,1)=tt11+tt13-tt14
!!                    psir(i1,i2,i3,2)=tt22+tt23+tt24
!!                    psir(i1,i2,i3,3)=tt33+tt31+tt32
!!                    psir(i1,i2,i3,4)=tt44+tt41-tt42
!!                 end do
!!              else
!!                 do ispinor=1,nspinor
!!                    do i1=i1s,i1e
!!                       !the local potential is always real
!!                       ! Add the parabola to the potential
!!                       !tt=hxh**2*dble(i1-ix0)**2 + hyh**2*dble(i2-iy0)**2 + hzh**2*dble(i3-iz0)**2
!!                       !tt=hxh2*dble(i1-ix0)**2 + hyh2*dble(i2-iy0)**2 + hzh2*dble(i3-iz0)**2
!!                       tt=(hxh*dble(i1)-rxyzParab(1))**2 + (hyh*dble(i2)-rxyzParab(2))**2 + (hzh*dble(i3)-rxyzParab(3))**2
!!                       if(power==2) then
!!                           tt=parabPrefac*tt
!!                       else if(power==4) then
!!                           tt=parabPrefac*tt**2
!!                       else
!!                           write(*,'(a,i0)') "'power' must be 2 or 4, but we found ", power
!!                           stop
!!                       end if
!!                       potCopy(i1-2*nbuf,i2-2*nbuf,i3-2*nbuf,1)=potCopy(i1-2*nbuf,i2-2*nbuf,i3-2*nbuf,1)+tt
!!                       tt=potCopy(i1-2*nbuf,i2-2*nbuf,i3-2*nbuf,1)*psir(i1,i2,i3,ispinor)
!!                       epot_p=epot_p+real(tt*psir(i1,i2,i3,ispinor),gp)
!!                       psir(i1,i2,i3,ispinor)=tt
!!                    end do
!!                 end do
!!              end if
!!              
!!              if (present(ibyyzz_r)) then
!!                 !the max is to avoid the calculation for no bounds
!!                 do i1=max(ibyyzz_r(1,i2,i3),ibyyzz_r(2,i2,i3))-14+1,2*n1+16-2*nbuf
!!                    psir(i1,i2,i3,:)=0.0_wp
!!                 enddo
!!              end if
!!
!!           else
!!              do i1=-14,2*n1+16
!!                 psir(i1,i2,i3,:)=0.0_wp
!!              enddo
!!           endif
!!        enddo
!!     else
!!        do i2=-14,2*n2+16
!!           do i1=-14,2*n1+16
!!              psir(i1,i2,i3,:)=0.0_wp
!!           enddo
!!        enddo
!!     endif
!!  enddo
!!!$omp end do
!!
!!!$omp critical
!!  epot=epot+epot_p
!!!$omp end critical
!!
!!!$omp end parallel
!!
!!END SUBROUTINE apply_potentialParabola
