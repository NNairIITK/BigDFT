subroutine linearScaling(iproc,nproc,Glr,orbs,comms,tmb,tmbder,at,input,hx,hy,hz,&
     rxyz,fion,fdisp,denspot,rhopotold,nlpspd,proj,GPU,&
     energs,scpot,psi,energy)
use module_base
use module_types
use module_interfaces, exceptThisOne => linearScaling
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(locreg_descriptors),intent(in) :: Glr
type(orbitals_data),intent(inout):: orbs
type(communications_arrays),intent(in) :: comms
type(atoms_data),intent(inout):: at
type(input_variables),intent(in):: input
real(8),dimension(3,at%nat),intent(inout):: rxyz
real(8),dimension(3,at%nat),intent(in):: fion, fdisp
type(DFT_local_fields), intent(inout) :: denspot
real(gp), dimension(:), intent(inout) :: rhopotold
type(nonlocal_psp_descriptors),intent(in):: nlpspd
real(wp),dimension(nlpspd%nprojel),intent(inout):: proj
type(GPU_pointers),intent(in out):: GPU
type(energy_terms),intent(inout) :: energs
real(gp),intent(in):: hx,hy,hz
logical,intent(in):: scpot
real(8),dimension(:),pointer,intent(out):: psi
real(gp), dimension(:), pointer :: rho,pot
real(8),intent(out):: energy
type(DFT_wavefunction),intent(inout),target:: tmb
type(DFT_wavefunction),intent(inout),target:: tmbder

type(linear_scaling_control_variables):: lscv
real(8):: pnrm,trace,fnrm_tmb
integer:: infoCoeff,istat,iall,it_scc,ilr,tag,itout,iorb,ist,iiorb,ncnt,p2p_tag,scf_mode,info_scf, nit_highaccur
character(len=*),parameter:: subname='linearScaling'
real(8),dimension(:),pointer :: psit
real(8),dimension(:),allocatable:: rhopotold_out
real(8):: energyold, energyDiff, energyoldout
type(mixrhopotDIISParameters):: mixdiis
type(localizedDIISParameters):: ldiis, ldiis_coeff
type(DFT_wavefunction),pointer:: tmbmix
logical:: check_whether_derivatives_to_be_used,coeffs_copied, first_time_with_der,calculate_overlap_matrix
integer:: jorb, jjorb, iiat
real(8),dimension(:,:),allocatable:: density_kernel, overlapmatrix
!FOR DEBUG ONLY
!integer,dimension(:),allocatable:: debugarr
real(8),dimension(:),allocatable :: locrad_tmp, eval
type(DFT_wavefunction):: tmblarge, tmblargeder, tmblarge2
real(8),dimension(:,:),allocatable:: locregCenter, locregCenterTemp, kernel, Umat
real(8),dimension(:),pointer:: lhphilarge, lhphilargeold, lphilargeold, lhphilargeder, lhphilargeoldder, lphilargeoldder
real(8),dimension(:),pointer:: lhphilarge2, lhphilarge2old, lphilarge2old
integer,dimension(:),allocatable:: onwhichatom_reference, inwhichlocreg_reference


  if(iproc==0) then
      write(*,'(1x,a)') repeat('*',84)
      write(*,'(1x,a)') '****************************** LINEAR SCALING VERSION ******************************'
  end if

  ! Initialize everything related to the linear scaling version ###########################################################
  ! Initialize the tags for the p2p communication

  ! Now all initializations are done ######################################################################################


     !!!! Allocate the transposed TMBs
     !!!allocate(tmb%psit_c(tmb%collcom%ndimind_c), stat=istat)
     !!!call memocc(istat, tmb%psit_c, 'tmb%psit_c', subname)
     !!!allocate(tmb%psit_f(7*tmb%collcom%ndimind_f), stat=istat)
     !!!call memocc(istat, tmb%psit_f, 'tmb%psit_f', subname)
     !!!allocate(overlapmatrix(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
     !!!call memocc(istat, overlapmatrix, 'overlapmatrix', subname)
     !!!! Deallocate the transposed TMBs
     !!!iall=-product(shape(tmb%psit_c))*kind(tmb%psit_c)
     !!!deallocate(tmb%psit_c, stat=istat)
     !!!call memocc(istat, iall, 'tmb%psit_c', subname)
     !!!iall=-product(shape(tmb%psit_f))*kind(tmb%psit_f)
     !!!deallocate(tmb%psit_f, stat=istat)
     !!!call memocc(istat, iall, 'tmb%psit_f', subname)
     !!!iall=-product(shape(overlapmatrix))*kind(overlapmatrix)
     !!!deallocate(overlapmatrix, stat=istat)
     !!!call memocc(istat, iall, 'overlapmatrix', subname)
     !!!!!call random_seed()
     !!!!!call random_number(tmb%psi)
     !!!call memocc(istat, density_kernel, 'density_kernel', subname)
     !!!call calculate_density_kernel(iproc, nproc, tmb%orbs%norb, orbs%norb, orbs%norbp, orbs%isorb, &
     !!!     tmb%wfnmd%ld_coeff, tmb%wfnmd%coeff, density_kernel)
     !!!call sumrhoForLocalizedBasis2(iproc, nproc, &
     !!!iall = -product(shape(density_kernel))*kind(density_kernel)
     !!!deallocate(density_kernel,stat=istat)
     !!!call memocc(istat,iall,'density_kernel',subname)
! DEBUG (SEE IF HAMILTONIAN IS GOOD)
!!     call allocateCommunicationsBuffersPotential(tmb%comgp, subname)
!!     call post_p2p_communication(iproc, nproc, denspot%dpbox%ndimpot, denspot%rhov, &
!!          tmb%comgp%nrecvbuf, tmb%comgp%recvbuf, tmb%comgp)
!!     call full_local_potential(iproc,nproc,tmb%orbs,tmb%lzd,2,&
!!          denspot%dpbox,denspot%rhov,denspot%pot_work,tmb%comgp)
!!     call deallocateCommunicationsBuffersPotential(tmb%comgp, subname)
!!     allocate(confdatarr(tmb%orbs%norbp))
!!     call define_confinement_data(confdatarr,tmb%orbs,rxyz,at,hx,hy,hz,input%lin%confpotorder,&
!!          input%lin%potentialprefac_lowaccuracy,tmb%lzd,tmb%orbs%onwhichatom)
!!     allocate(lhchi(max(tmb%orbs%npsidim_orbs,tmb%orbs%npsidim_comp)),stat=istat)
!!     call memocc(istat, lhchi, 'lhchi', subname)
!!     call to_zero(max(tmb%orbs%npsidim_orbs,tmb%orbs%npsidim_comp),lhchi(1))
!!     allocate(tmb%lzd%doHamAppl(tmb%lzd%nlr))
!!     tmb%lzd%doHamAppl = .true.
!!     call LocalHamiltonianApplication(iproc,nproc,at,tmb%orbs,&
!!          tmb%lzd,confdatarr,denspot%dpbox%ngatherarr,denspot%pot_work,tmb%psi,lhchi(1),&
!!          energs,input%SIC,GPU,.false.,&
!!          pkernel=denspot%pkernelseq)
!!     call NonLocalHamiltonianApplication(iproc,at,tmb%orbs,&
!!          rxyz,proj,tmb%lzd,nlpspd,tmb%psi,lhchi(1),energs%eproj)
!!    call total_energies(energs,1)
!!    print *,'ebs,ekin,epot,eproj',energs%ebs,energs%ekin,energs%epot,energs%eproj
!!    print *,'eKS,ehart,exc,evxc',energs%ebs-energs%eh+energs%exc-energs%evxc-eexctX+eion+edisp,energs%eh,energs%exc,energs%evxc
!!    print *,'Wavefunction coefficients:'
!!    do iall = 1, orbs%norb
!!       do istat = 1, tmb%orbs%norb
!!          print *,iall,istat,tmb%wfnmd%coeff(istat,iall)
!!       end do
!!    end do
!!    ind1 = 1
!!    print *,'Hamiltonian matrix:'
!!    do iall = 1, tmb%lzd%nlr
!!       ind2 = 1
!!       do istat = 1, tmb%lzd%nlr
!!          call wpdot_wrap(1,tmb%lzd%Llr(iall)%wfd%nvctr_c,tmb%Lzd%Llr(iall)%wfd%nvctr_f,tmb%Lzd%Llr(iall)%wfd%nseg_c,&
!!               tmb%Lzd%Llr(iall)%wfd%nseg_f,tmb%Lzd%Llr(iall)%wfd%keyvglob,tmb%Lzd%Llr(iall)%wfd%keyglob,tmb%psi(ind1),&
!!               tmb%Lzd%Llr(istat)%wfd%nvctr_c,tmb%Lzd%Llr(istat)%wfd%nvctr_f,tmb%Lzd%Llr(istat)%wfd%nseg_c,&
!!               tmb%Lzd%Llr(istat)%wfd%nseg_f,tmb%Lzd%Llr(istat)%wfd%keyvglob,tmb%Lzd%Llr(istat)%wfd%keyglob,tmb%psi(ind2),overlap)
!!          call wpdot_wrap(1,tmb%lzd%Llr(iall)%wfd%nvctr_c,tmb%Lzd%Llr(iall)%wfd%nvctr_f,tmb%Lzd%Llr(iall)%wfd%nseg_c,&
!!               tmb%Lzd%Llr(iall)%wfd%nseg_f,tmb%Lzd%Llr(iall)%wfd%keyvglob,tmb%Lzd%Llr(iall)%wfd%keyglob,tmb%psi(ind1),&
!!               tmb%Lzd%Llr(istat)%wfd%nvctr_c,tmb%Lzd%Llr(istat)%wfd%nvctr_f,tmb%Lzd%Llr(istat)%wfd%nseg_c,&
!!               tmb%Lzd%Llr(istat)%wfd%nseg_f,tmb%Lzd%Llr(istat)%wfd%keyvglob,tmb%Lzd%Llr(istat)%wfd%keyglob,lhchi(ind2),hamil)
!!          ind2 = ind2 + tmb%Lzd%Llr(istat)%wfd%nvctr_c + 7*tmb%Lzd%Llr(istat)%wfd%nvctr_f
!!          print *,iall,istat,overlap,hamil
!!       end do
!!       ind1 = ind1+tmb%lzd%Llr(iall)%wfd%nvctr_c+7*tmb%Lzd%Llr(iall)%wfd%nvctr_f
!!    end do
!!stop
!END DEBUG






  !! Now one could calculate the charge density like this. It is not done since we would in this way overwrite
  !! the potential from the input guess.     
  !call allocateCommunicationbufferSumrho(iproc, with_auxarray, tmb%comsr, subname)
  !call communicate_basis_for_density(iproc, nproc, tmb%lzd, tmb%orbs, tmb%psi, tmb%comsr)
  !call sumrhoForLocalizedBasis2(iproc, nproc, orbs%norb, tmb%lzd, input, hx, hy, hz, tmb%orbs, tmb%comsr, &
  !     tmb%wfnmd%ld_coeff, tmb%wfnmd%coeff, Glr%d%n1i*Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov, at,denspot%dpbox%nscatterarr)
  !call deallocateCommunicationbufferSumrho(tmb%comsr, subname)

  allocate(lscv%locrad(tmb%lzd%nlr), stat=istat)
  call memocc(istat, lscv%locrad, 'lscv%locrad', subname)


  ! Allocate the old charge density (used to calculate the variation in the charge density)
  allocate(rhopotold_out(max(glr%d%n1i*glr%d%n2i*denspot%dpbox%n3p,1)*input%nspin), stat=istat)
  call memocc(istat, rhopotold_out, 'rhopotold_out', subname)

  if(input%lin%nItInguess>0) then
      tmb%wfnmd%bs%communicate_phi_for_lsumrho=.true.
      tmb%wfnmd%bs%target_function=TARGET_FUNCTION_IS_TRACE

      do ilr=1,tmb%lzd%nlr
          lscv%locrad(ilr)=max(input%lin%locrad_lowaccuracy(ilr),tmb%lzd%llr(ilr)%locrad)
      end do

      !if(trim(input%lin%mixingMethod)=='dens') then
      if(input%lin%scf_mode==LINEAR_MIXDENS_SIMPLE) then
          rhopotold_out=rhopotold
      end if

      !if(trim(input%lin%mixingMethod)=='pot') then
      if(input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE) then
          rhopotold_out=denspot%rhov
      end if

      ! Copy the current potential
      !if(trim(input%lin%mixingMethod)=='pot') then
      if(input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE) then
           call dcopy(max(Glr%d%n1i*Glr%d%n2i*denspot%dpbox%n3p,1)*input%nspin, denspot%rhov(1), 1, rhopotOld(1), 1)
      end if
  end if



  ! Allocate the communications buffers needed for the communications of teh potential and
  ! post the messages. This will send to each process the part of the potential that this process
  ! needs for the application of the Hamlitonian to all orbitals on that process.
  call allocateCommunicationsBuffersPotential(tmb%comgp, subname)
  call allocateCommunicationsBuffersPotential(tmbder%comgp, subname)

  ! Initialize the DIIS mixing of the potential if required.
  if(input%lin%mixHist_lowaccuracy>0) then
      call initializeMixrhopotDIIS(input%lin%mixHist_lowaccuracy, denspot%dpbox%ndimpot, mixdiis)
  end if

  ! Flag that indicates that the basis functions shall be improved in the following.
  tmb%wfnmd%bs%update_phi=.true.
  pnrm=1.d100
  lscv%pnrm_out=1.d100
  energyold=0.d0
  energyoldout=0.d0
  lscv%reduce_convergence_tolerance=.false.
  tmb%wfnmd%bs%target_function=TARGET_FUNCTION_IS_TRACE
  lscv%lowaccur_converged=.false.
  lscv%info_basis_functions=-1
  lscv%idecrease=0
  lscv%increase_locreg=0.d0
  lscv%decrease_factor_total=1.d10 !initialize to some large value
  lscv%ifail=0
  lscv%enlarge_locreg=.false.

  ! tmbmix is the types we use for the mixing. It will point to either tmb if we don't use the derivatives
  ! or to tmbder if we use the derivatives.
  if(input%lin%useDerivativeBasisFunctions) then
      tmbmix => tmbder
  else
      tmbmix => tmb
  end if

  ! Check whether it is possible to have variable localization regions or not.
  if(tmb%wfnmd%bs%nit_unitary_loop==-1 .and. tmb%wfnmd%bs%locreg_enlargement==1.d0) then
      lscv%variable_locregs=.false.
  else
      lscv%variable_locregs=.true.
  end if


  !!! Allocate the communication arrays for the calculation of the charge density.
  call allocateCommunicationbufferSumrho(iproc, tmb%comsr, subname)
  call allocateCommunicationbufferSumrho(iproc, tmbder%comsr, subname)

  ! Initialize DIIS...
  !!if(.not.lscv%lowaccur_converged) then
      call initializeDIIS(input%lin%DIISHistMax, tmb%lzd, tmb%orbs, tmb%orbs%norb, ldiis)
      ldiis%DIISHistMin=input%lin%DIISHistMin
      ldiis%DIISHistMax=input%lin%DIISHistMax
      ldiis%alphaSD=input%lin%alphaSD
      ldiis%alphaDIIS=input%lin%alphaDIIS
      ldiis%icountSDSatur=0
      ldiis%icountSwitch=0
      ldiis%icountDIISFailureTot=0
      ldiis%icountDIISFailureCons=0
      ldiis%is=0
      ldiis%switchSD=.false.
      ldiis%trmin=1.d100
      ldiis%trold=1.d100
  !!end if


  ! just to be sure...
  nullify(tmb%psit_c)
  nullify(tmb%psit_f)
  nullify(tmbder%psit_c)
  nullify(tmbder%psit_f)


  allocate(eval(tmb%orbs%norb), stat=istat)
  call memocc(istat, eval, 'eval', subname)
  call vcopy(tmb%orbs%norb, tmb%orbs%eval(1), 1, eval(1), 1)

  ! This is the main outer loop. Each iteration of this loop consists of a first loop in which the basis functions
  ! are optimized and a consecutive loop in which the density is mixed.
  coeffs_copied=.false.
  first_time_with_der=.false.
  nit_highaccur=0



  outerLoop: do itout=1,input%lin%nit_lowaccuracy+input%lin%nit_highaccuracy
      !!if(iproc==0) write(*,*) 'START LOOP: ldiis%hphiHist(1)',ldiis%hphiHist(1)


      ! First to some initialization and determine the value of some control parameters.

      !!! Initialize DIIS...
      !!! Keep the history in the high accuracy case.
      !!if(.not.lscv%lowaccur_converged) then
      !!    if(itout>1) call deallocateDIIS(ldiis)
      !!    call initializeDIIS(input%lin%DIISHistMax, tmb%lzd, tmb%orbs, tmb%orbs%norb, ldiis)
      !!    ldiis%DIISHistMin=input%lin%DIISHistMin
      !!    ldiis%DIISHistMax=input%lin%DIISHistMax
      !!    ldiis%alphaSD=input%lin%alphaSD
      !!    ldiis%alphaDIIS=input%lin%alphaDIIS
      !!    ldiis%icountSDSatur=0
      !!    ldiis%icountSwitch=0
      !!    ldiis%icountDIISFailureTot=0
      !!    ldiis%icountDIISFailureCons=0
      !!    ldiis%is=0
      !!    ldiis%switchSD=.false.
      !!    ldiis%trmin=1.d100
      !!    ldiis%trold=1.d100
      !!end if

      ! The basis functions shall be optimized
      tmb%wfnmd%bs%update_phi=.true.

      ! Convergence criterion for the self consistency loop
      lscv%self_consistent=input%lin%convCritMix


      ! Check whether the low accuracy part (i.e. with strong confining potential) has converged.
      call check_whether_lowaccuracy_converged(itout, input, lscv)

      ! Check whether the derivatives shall be used or not.
      lscv%withder=check_whether_derivatives_to_be_used(input, itout, lscv)

      if(lscv%withder .and. lscv%lowaccur_converged .and. .not.coeffs_copied) then
          !!tmbder%wfnmd%coeff=0.d0
          call to_zero(tmbder%orbs%norb*orbs%norb, tmbder%wfnmd%coeff(1,1))
          do iorb=1,orbs%norb
              jjorb=0
              do jorb=1,tmbder%orbs%norb,4
                  jjorb=jjorb+1
                  tmbder%wfnmd%coeff(jorb,iorb)=tmb%wfnmd%coeff(jjorb,iorb)
              end do
          end do
          coeffs_copied=.true.
      end if

      

      ! Set all remaining variables that we need for the optimizations of the basis functions and the mixing.
      call set_optimization_variables(input, at, tmb%orbs, tmb%lzd%nlr, tmb%orbs%onwhichatom, &
           tmb%confdatarr, tmb%wfnmd, lscv)
      call set_optimization_variables(input, at, tmbder%orbs, tmb%lzd%nlr, tmbder%orbs%onwhichatom, &
           tmbder%confdatarr, tmbder%wfnmd, lscv)

      if(lscv%lowaccur_converged) nit_highaccur=nit_highaccur+1
      if(nit_highaccur==1) lscv%enlarge_locreg=.true.


      !!if(iproc==0) write(*,*) 'MIDDLE 1: ldiis%hphiHist(1)',ldiis%hphiHist(1)
      ! Adjust the confining potential if required.
      call adjust_locregs_and_confinement(iproc, nproc, hx, hy, hz, &
           input, tmb, tmbder, denspot, ldiis, lscv)
      !!if(iproc==0) write(*,*) 'MIDDLE 2: ldiis%hphiHist(1)',ldiis%hphiHist(1)

      ! Some special treatement if we are in the high accuracy part
      call adjust_DIIS_for_high_accuracy(input, tmb, denspot, ldiis, mixdiis, lscv)
      !!if(lscv%exit_outer_loop) exit outerLoop

      !!if(iproc==0) write(*,*) 'MIDDLE: ldiis%hphiHist(1)',ldiis%hphiHist(1)

      if(lscv%withder) then
          call initialize_DIIS_coeff(3, tmbder, orbs, ldiis_coeff)
      else
          call initialize_DIIS_coeff(3, tmb, orbs, ldiis_coeff)
      end if

      ! Now all initializations are done...



      allocate(locregCenter(3,tmb%lzd%nlr), stat=istat)
      call memocc(istat, locregCenter, 'locregCenter', subname)
      allocate(locrad_tmp(tmb%lzd%nlr), stat=istat)
      call memocc(istat, locrad_tmp, 'locrad_tmp', subname)
      !!$$allocate(onwhichatom_reference(tmb%orbs%norb), stat=istat)
      !!$$call memocc(istat, onwhichatom_reference, 'onwhichatom_reference', subname)
      !!$$allocate(inwhichlocreg_reference(tmb%orbs%norb), stat=istat)
      !!$$call memocc(istat, inwhichlocreg_reference, 'inwhichlocreg_reference', subname)
      !!$$allocate(locregCenterTemp(3,tmb%lzd%nlr), stat=istat)
      !!$$call memocc(istat, locregCenterTemp, 'locregCenterTemp', subname)
      !!$$allocate(kernel(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
      !!$$call memocc(istat, kernel, 'kernel', subname)
      !!$$allocate(Umat(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
      !!$$call memocc(istat, Umat, 'Umat', subname)


!!$$      !## TEST ###################################################
!!$$
!!$$  call vcopy(tmb%orbs%norb, tmb%orbs%inwhichlocreg(1), 1, inwhichlocreg_reference(1), 1)
!!$$
!!$$
!!$$  ! Initialize largestructures if required
!!$$      do iorb=1,tmb%orbs%norb
!!$$          ilr=tmb%orbs%inwhichlocreg(iorb)
!!$$          locregCenter(:,ilr)=tmb%lzd%llr(ilr)%locregCenter
!!$$      end do
!!$$      locregCenterTemp=locregCenter
!!$$      do ilr=1,tmb%lzd%nlr
!!$$          locrad_tmp(ilr)=tmb%lzd%llr(ilr)%locrad+4.d0
!!$$      end do
!!$$      call update_locreg(iproc, nproc, tmb%lzd%nlr, locrad_tmp, inwhichlocreg_reference, locregCenter, tmb%lzd%glr, &
!!$$           .false., denspot%dpbox%nscatterarr, tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3), &
!!$$           tmb%orbs, tmblarge2%lzd, tmblarge2%orbs, tmblarge2%op, tmblarge2%comon, &
!!$$           tmblarge2%comgp, tmblarge2%comsr, tmblarge2%mad, tmblarge2%collcom)
!!$$      call update_ldiis_arrays(tmblarge2, subname, ldiis)
!!$$      call allocate_auxiliary_basis_function(tmblarge2%orbs%npsidim_orbs, subname, tmblarge2%psi, &
!!$$           lhphilarge2, lhphilarge2old, lphilarge2old)
!!$$      call copy_basis_performance_options(tmb%wfnmd%bpo, tmblarge2%wfnmd%bpo, subname)
!!$$      call copy_orthon_data(tmb%orthpar, tmblarge2%orthpar, subname)
!!$$      tmblarge2%wfnmd%nphi=tmblarge2%orbs%npsidim_orbs
!!$$      call vcopy(tmb%orbs%norb, tmb%orbs%onwhichatom(1), 1, onwhichatom_reference(1), 1)
!!$$      !call vcopy(tmb%orbs%norb, onwhichatom_reference(1), 1, tmblarge2%orbs%onwhichatom(1), 1)
!!$$
!!$$      call small_to_large_locreg(iproc, nproc, tmb%lzd, tmblarge2%lzd, tmb%orbs, tmblarge2%orbs, tmb%psi, tmblarge2%psi)
!!$$      allocate(tmblarge2%confdatarr(tmblarge2%orbs%norbp), stat=istat)
!!$$      call vcopy(tmb%orbs%norb, tmb%orbs%onwhichatom(1), 1, tmblarge2%orbs%onwhichatom(1), 1)
!!$$      if(.not.lscv%lowaccur_converged) then
!!$$          call define_confinement_data(tmblarge2%confdatarr,tmblarge2%orbs,rxyz,at,&
!!$$               tmblarge2%lzd%hgrids(1),tmblarge2%lzd%hgrids(2),tmblarge2%lzd%hgrids(3),&
!!$$               input%lin%ConfPotOrder,input%lin%potentialPrefac_lowaccuracy,tmblarge2%lzd,tmblarge2%orbs%onwhichatom)
!!$$      else
!!$$          call define_confinement_data(tmblarge2%confdatarr,tmblarge2%orbs,rxyz,at,&
!!$$               tmblarge2%lzd%hgrids(1),tmblarge2%lzd%hgrids(2),tmblarge2%lzd%hgrids(3),&
!!$$               input%lin%ConfPotOrder,input%lin%potentialPrefac_highaccuracy,tmblarge2%lzd,tmblarge2%orbs%onwhichatom)
!!$$      end if
!!$$
!!$$      call MLWFnew(iproc, nproc, tmblarge2%lzd, tmblarge2%orbs, at, tmblarge2%op, &
!!$$           tmblarge2%comon, tmblarge2%mad, rxyz, 0, kernel, &
!!$$           tmblarge2%confdatarr, tmb%lzd%hgrids(1), locregCenterTemp, 3.d0, tmblarge2%psi, Umat, locregCenter)
!!$$      deallocate(tmblarge2%confdatarr, stat=istat)
!!$$
!!$$      call check_locregCenters(iproc, tmb%lzd, locregCenter, tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3))
!!$$
!!$$      call vcopy(tmb%orbs%norb, tmb%orbs%onwhichatom(1), 1, onwhichatom_reference(1), 1)
!!$$      do ilr=1,tmb%lzd%nlr
!!$$          locrad_tmp(ilr)=tmb%lzd%llr(ilr)%locrad
!!$$      end do
!!$$      call destroy_new_locregs(iproc, nproc, tmb)
!!$$      call update_locreg(iproc, nproc, tmb%lzd%nlr, locrad_tmp, inwhichlocreg_reference, locregCenter, tmblarge2%lzd%glr, &
!!$$           .false., denspot%dpbox%nscatterarr, tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3), &
!!$$           tmblarge2%orbs, tmb%lzd, tmb%orbs, tmb%op, tmb%comon, &
!!$$           tmb%comgp, tmb%comsr, tmb%mad, tmb%collcom)
!!$$      call update_ldiis_arrays(tmb, subname, ldiis)
!!$$
!!$$      allocate(tmb%psi(tmb%orbs%npsidim_orbs), stat=istat)
!!$$      call memocc(istat, tmb%psi, 'tmb%psi', subname)
!!$$
!!$$      !!call update_auxiliary_basis_function(subname, tmb%orbs%npsidim_orbs, tmb%psi, lhphi, lhphiold, lphiold)
!!$$      call copy_basis_performance_options(tmblarge2%wfnmd%bpo, tmb%wfnmd%bpo, subname)
!!$$      call copy_orthon_data(tmblarge2%orthpar, tmb%orthpar, subname)
!!$$      call vcopy(tmb%orbs%norb, onwhichatom_reference(1), 1, tmb%orbs%onwhichatom(1), 1)
!!$$      tmb%wfnmd%nphi=tmb%orbs%npsidim_orbs
!!$$
!!$$      call large_to_small_locreg(iproc, nproc, tmb%lzd, tmblarge2%lzd, tmb%orbs, tmblarge2%orbs, tmblarge2%psi, tmb%psi)
!!$$      call destroy_new_locregs(iproc, nproc, tmblarge2)
!!$$
!!$$      iall=-product(shape(onwhichatom_reference))*kind(onwhichatom_reference)
!!$$      deallocate(onwhichatom_reference, stat=istat)
!!$$      call memocc(istat, iall, 'onwhichatom_reference', subname)
!!$$      iall=-product(shape(inwhichlocreg_reference))*kind(inwhichlocreg_reference)
!!$$      deallocate(inwhichlocreg_reference, stat=istat)
!!$$      call memocc(istat, iall, 'inwhichlocreg_reference', subname)
!!$$      iall=-product(shape(locregCenterTemp))*kind(locregCenterTemp)
!!$$      deallocate(locregCenterTemp, stat=istat)
!!$$      call memocc(istat, iall, 'locregCenterTemp', subname)
!!$$      iall=-product(shape(kernel))*kind(kernel)
!!$$      deallocate(kernel, stat=istat)
!!$$      call memocc(istat, iall, 'kernel', subname)
!!$$      iall=-product(shape(Umat))*kind(Umat)
!!$$      deallocate(Umat, stat=istat)
!!$$      call memocc(istat, iall, 'Umat', subname)
!!$$
!!$$      !## END TEST ###################################################








      do iorb=1,tmb%orbs%norb
          ilr=tmb%orbs%inwhichlocreg(iorb)
          locregCenter(:,ilr)=tmb%lzd%llr(ilr)%locregCenter
      end do
      do ilr=1,tmb%lzd%nlr
          locrad_tmp(ilr)=tmb%lzd%llr(ilr)%locrad+8.d0*tmb%lzd%hgrids(1)
      end do

      call update_locreg(iproc, nproc, tmb%lzd%nlr, locrad_tmp, tmb%orbs%inwhichlocreg, locregCenter, tmb%lzd%glr, &
           .false., denspot%dpbox%nscatterarr, tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3), &
           tmb%orbs, tmblarge%lzd, tmblarge%orbs, tmblarge%op, tmblarge%comon, &
           tmblarge%comgp, tmblarge%comsr, tmblarge%mad, tmblarge%collcom)
      call allocate_auxiliary_basis_function(max(tmblarge%orbs%npsidim_comp,tmblarge%orbs%npsidim_orbs), subname, &
           tmblarge%psi, lhphilarge, lhphilargeold, lphilargeold)
      call copy_basis_performance_options(tmb%wfnmd%bpo, tmblarge%wfnmd%bpo, subname)
      call copy_orthon_data(tmb%orthpar, tmblarge%orthpar, subname)
      tmblarge%wfnmd%nphi=tmblarge%orbs%npsidim_orbs
      tmblarge%can_use_transposed=.false.
      nullify(tmblarge%psit_c)
      nullify(tmblarge%psit_f)
      allocate(tmblarge%confdatarr(tmblarge%orbs%norbp), stat=istat)
      !call memocc(istat, tmblarge%confdatarr, 'tmblarge%confdatarr', subname)
      ! copy onwhichatom... maybe to be done somewhere else
      call vcopy(tmb%orbs%norb, tmb%orbs%onwhichatom(1), 1, tmblarge%orbs%onwhichatom(1), 1)
      if(.not.lscv%lowaccur_converged) then
          call define_confinement_data(tmblarge%confdatarr,tmblarge%orbs,rxyz,at,&
               tmblarge%lzd%hgrids(1),tmblarge%lzd%hgrids(2),tmblarge%lzd%hgrids(3),&
               input%lin%ConfPotOrder,input%lin%potentialPrefac_lowaccuracy,tmblarge%lzd,tmblarge%orbs%onwhichatom)
      else
          call define_confinement_data(tmblarge%confdatarr,tmblarge%orbs,rxyz,at,&
               tmblarge%lzd%hgrids(1),tmblarge%lzd%hgrids(2),tmblarge%lzd%hgrids(3),&
               input%lin%ConfPotOrder,input%lin%potentialPrefac_highaccuracy,tmblarge%lzd,tmblarge%orbs%onwhichatom)
      end if
      !write(*,*) 'tmb%confdatarr(1)%ioffset(:), tmblarge%confdatarr(1)%ioffset(:)',tmb%confdatarr(1)%ioffset(:), tmblarge%confdatarr(1)%ioffset(:)

      ! take the eigenvalues from the input guess for the preconditioning
      call vcopy(tmb%orbs%norb, eval(1), 1, tmblarge%orbs%eval(1), 1)




      if(lscv%withder) then
          call update_locreg(iproc, nproc, tmb%lzd%nlr, locrad_tmp, tmbder%orbs%inwhichlocreg, locregCenter, tmb%lzd%glr, &
               .false., denspot%dpbox%nscatterarr, tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3), &
               tmbder%orbs, tmblargeder%lzd, tmblargeder%orbs, tmblargeder%op, tmblargeder%comon, &
               tmblargeder%comgp, tmblargeder%comsr, tmblargeder%mad, tmblargeder%collcom)
          call allocate_auxiliary_basis_function(max(tmblargeder%orbs%npsidim_comp,tmblargeder%orbs%npsidim_orbs), subname, &
               tmblargeder%psi, lhphilargeder, lhphilargeoldder, lphilargeoldder)
          call copy_basis_performance_options(tmbder%wfnmd%bpo, tmblargeder%wfnmd%bpo, subname)
          call copy_orthon_data(tmbder%orthpar, tmblargeder%orthpar, subname)
          tmblargeder%wfnmd%nphi=tmblargeder%orbs%npsidim_orbs
          tmblargeder%can_use_transposed=.false.
          nullify(tmblargeder%psit_c)
          nullify(tmblargeder%psit_f)
          allocate(tmblargeder%confdatarr(tmblargeder%orbs%norbp), stat=istat)
          !call memocc(istat, tmblargeder%confdatarr, 'tmblargeder%confdatarr', subname)
          ! copy onwhichatom... maybe to be done somewhere else
          call vcopy(tmb%orbs%norb, tmbder%orbs%onwhichatom(1), 1, tmblargeder%orbs%onwhichatom(1), 1)
          if(.not.lscv%lowaccur_converged) then
              call define_confinement_data(tmblargeder%confdatarr,tmblargeder%orbs,rxyz,at,&
                   tmblargeder%lzd%hgrids(1),tmblargeder%lzd%hgrids(2),tmblargeder%lzd%hgrids(3),&
                   input%lin%ConfPotOrder,input%lin%potentialPrefac_lowaccuracy,tmblargeder%lzd,tmblargeder%orbs%onwhichatom)
          else
              call define_confinement_data(tmblargeder%confdatarr,tmblargeder%orbs,rxyz,at,&
                   tmblargeder%lzd%hgrids(1),tmblargeder%lzd%hgrids(2),tmblargeder%lzd%hgrids(3),&
                   input%lin%ConfPotOrder,input%lin%potentialPrefac_highaccuracy,tmblargeder%lzd,tmblargeder%orbs%onwhichatom)
          end if
      end if

      if(itout==1) then
          ! Orthonormalize the TMBs
          ! just to be sure...
          tmb%can_use_transposed=.false.
          nullify(tmb%psit_c)
          nullify(tmb%psit_f)
          call orthonormalizeLocalized(iproc, nproc, tmb%orthpar%methTransformOverlap, tmb%orthpar%nItOrtho, &
               tmb%orbs, tmb%op, tmb%comon, tmb%lzd, &
               tmb%mad, tmb%collcom, tmb%orthpar, tmb%wfnmd%bpo, tmb%psi, tmb%psit_c, tmb%psit_f, &
               tmb%can_use_transposed)
      end if


      ! The self consistency cycle. Here we try to get a self consistent density/potential.
      ! In the first lscv%nit_scc_when_optimizing iteration, the basis functions are optimized, whereas in the remaining
      ! iteration the basis functions are fixed.
      do it_scc=1,lscv%nit_scc

          !!if(lscv%withder .and. .not.first_time_with_der) then
          !!    first_time_with_der=.true.
          !!    !scf_mode=LINEAR_MIXDENS_SIMPLE
          !!    scf_mode=input%lin%scf_mode
          !!    call transform_coeffs_to_derivatives(iproc, nproc, orbs, tmb%lzd, tmb, tmbder)
          !!else
          !!    scf_mode=input%lin%scf_mode
          !!end if
          scf_mode=input%lin%scf_mode

          ! Do not update the TMB if it_scc>lscv%nit_scc_when_optimizing
          if(it_scc>lscv%nit_scc_when_optimizing) tmb%wfnmd%bs%update_phi=.false.


          !!call post_p2p_communication(iproc, nproc, denspot%dpbox%ndimpot, denspot%rhov, &
          !!     tmb%comgp%nrecvbuf, tmb%comgp%recvbuf, tmb%comgp)
          if(lscv%withder) then
              call post_p2p_communication(iproc, nproc, denspot%dpbox%ndimpot, denspot%rhov, &
                   tmbder%comgp%nrecvbuf, tmbder%comgp%recvbuf, tmbder%comgp)
          end if


         ! Improve the trace minimizing orbitals.
          if(tmb%wfnmd%bs%update_phi) then
              if(tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY) then
                  do iorb=1,orbs%norb
                      call dcopy(tmb%orbs%norb, tmb%wfnmd%coeff_proj(1,iorb), 1, tmb%wfnmd%coeff(1,iorb), 1)
                  end do
              end if



              call getLocalizedBasis(iproc,nproc,at,orbs,rxyz,denspot,GPU,trace,fnrm_tmb,lscv%info_basis_functions,&
                  nlpspd,proj,ldiis,input%SIC,lscv%locrad,tmb, tmblarge, lhphilarge, lhphilargeold, lphilargeold)
              tmb%can_use_transposed=.false. !since basis functions have changed...
              tmbder%can_use_transposed=.false. !since basis functions have changed...
              !allocate(denspot%pot_work(tmblarge%lzd%ndimpotisf+ndebug),stat=istat)
              !call memocc(istat,denspot%pot_work,'denspot%pot_work',subname)



              tmb%wfnmd%nphi=tmb%orbs%npsidim_orbs
              !reset counter for optimization of coefficients (otherwise step size will be decreases...)
              tmb%wfnmd%it_coeff_opt=0
              tmbder%wfnmd%it_coeff_opt=0
              tmb%wfnmd%alpha_coeff=.2d0 !reset to default value
              tmbder%wfnmd%alpha_coeff=.2d0 !reset to default value

              !!write(*,*) 'nit_highaccur',nit_highaccur
              if(nit_highaccur<=1) then
                  call deallocateDIIS(ldiis)
                  call initializeDIIS(input%lin%DIISHistMax, tmb%lzd, tmb%orbs, tmb%orbs%norb, ldiis)
                  ldiis%DIISHistMin=input%lin%DIISHistMin
                  ldiis%DIISHistMax=input%lin%DIISHistMax
                  ldiis%alphaSD=input%lin%alphaSD
                  ldiis%alphaDIIS=input%lin%alphaDIIS
                  ldiis%icountSDSatur=0
                  ldiis%icountSwitch=0
                  ldiis%icountDIISFailureTot=0
                  ldiis%icountDIISFailureCons=0
                  ldiis%is=0
                  ldiis%switchSD=.false.
                  ldiis%trmin=1.d100
                  ldiis%trold=1.d100
              else
                  ! Since the potential changes, the values of ldiis%trmin should be reset.
                  ldiis%switchSD=.false.
                  ldiis%trmin=1.d100
                  ldiis%trold=1.d100
                  ldiis%icountSDSatur=0
                  ldiis%icountSwitch=0
                  ldiis%icountDIISFailureTot=0
                  ldiis%icountDIISFailureCons=0
              end if
          end if

          ! Initialize DIIS...
          ! Keep the history in the high accuracy case.
          !if(.not.lscv%lowaccur_converged) then

          if((lscv%locreg_increased .or. (lscv%variable_locregs .and. tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY)) &
              .and. tmb%wfnmd%bs%update_phi) then
              ! Redefine some quantities if the localization region has changed.
              if(lscv%withder) then
                  call redefine_locregs_quantities(iproc, nproc, hx, hy, hz, tmb%lzd%llr(:)%locrad, &
                       .false., tmb%lzd, tmb, tmbder, denspot)
                  call post_p2p_communication(iproc, nproc, denspot%dpbox%ndimpot, denspot%rhov, &
                       tmbder%comgp%nrecvbuf, tmbder%comgp%recvbuf, tmbder%comgp)
              end if
          end if

          ! Decide whether we have to use the derivatives or not.
          if(lscv%withder) then
              tmbmix => tmbder
          else
              tmbmix => tmb
          end if

          ! Build the derivatives if required.
          if(tmb%wfnmd%bs%update_phi .or. it_scc==0) then
              if(tmbmix%wfnmd%bs%use_derivative_basis) then
                  if((lscv%locreg_increased .or. &
                      (lscv%variable_locregs .and. tmb%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY)) &
                      .and. tmb%wfnmd%bs%update_phi) then
                      call deallocate_p2pComms(tmbder%comrp, subname)
                      call initializeRepartitionOrbitals(iproc, nproc, tag, tmb%orbs, tmbder%orbs, tmb%lzd, tmbder%comrp)
                      tmbmix => tmbder
                  end if
                  if(iproc==0) write(*,'(1x,a)',advance='no') 'calculating derivative basis functions...'
                  call getDerivativeBasisFunctions(iproc,nproc,hx,tmb%lzd,tmb%orbs,tmbmix%orbs,tmbmix%comrp,&
                       max(tmb%orbs%npsidim_orbs,tmb%orbs%npsidim_comp),tmb%psi,tmbmix%psi)
                  if(iproc==0) write(*,'(a)') 'done.'
                  !! TEST ###############################################################################################
                  !write(*,*) 'test: orthonormalize derivatives'
                  !!call orthonormalizeLocalized(iproc, nproc, tmb%orthpar%methTransformOverlap, tmb%orthpar%nItOrtho, &
                  !!     tmbder%orbs, tmbder%op, tmbder%comon, tmb%lzd, &
                  !!     tmbder%mad, tmbder%collcom, tmbder%orthpar, tmbder%wfnmd%bpo, tmbder%psi, tmbder%psit_c, tmbder%psit_f, &
                  !!     tmbder%can_use_transposed)
                  !!if(tmbder%can_use_transposed) then
                  !!    ! This is not optimal, these quantities will be recalculated...
                  !!    iall = -product(shape(tmbder%psit_c))*kind(tmbder%psit_c)
                  !!    deallocate(tmbder%psit_c,stat=istat)
                  !!    call memocc(istat,iall,'tmbder%psit_c',subname)
                  !!    iall = -product(shape(tmbder%psit_f))*kind(tmbder%psit_f)
                  !!    deallocate(tmbder%psit_f,stat=istat)
                  !!    call memocc(istat,iall,'tmbder%psit_f',subname)
                  !!end if
                  !! END TEST ###########################################################################################
              else
                  call dcopy(tmb%wfnmd%nphi, tmb%psi(1), 1, tmbmix%psi(1), 1)
              end if

              !!! Allocate the transposed TMBs
              !!allocate(tmbmix%psit_c(tmbmix%collcom%ndimind_c), stat=istat)
              !!call memocc(istat, tmbmix%psit_c, 'tmbmix%psit_c', subname)
              !!allocate(tmbmix%psit_f(7*tmbmix%collcom%ndimind_f), stat=istat)
              !!call memocc(istat, tmbmix%psit_f, 'tmbmix%psit_f', subname)
              allocate(overlapmatrix(tmbmix%orbs%norb,tmbmix%orbs%norb), stat=istat)
              call memocc(istat, overlapmatrix, 'overlapmatrix', subname)


          end if

          ! Only communicate the TMB for sumrho if required (i.e. only if the TMB were optimized).
          if(it_scc<=lscv%nit_scc_when_optimizing) then
              tmbmix%wfnmd%bs%communicate_phi_for_lsumrho=.true.
              calculate_overlap_matrix=.true.
          else
              tmbmix%wfnmd%bs%communicate_phi_for_lsumrho=.false.
              calculate_overlap_matrix=.false.
          end if

          if(lscv%withder .and. .not.first_time_with_der) then
              first_time_with_der=.true.
              !scf_mode=LINEAR_MIXDENS_SIMPLE
              scf_mode=input%lin%scf_mode
              call transform_coeffs_to_derivatives(iproc, nproc, orbs, tmb%lzd, tmb, tmbder)
              tmb%wfnmd%alpha_coeff=1.d-2
              tmbder%wfnmd%alpha_coeff=1.d-2
          else
              scf_mode=input%lin%scf_mode
          end if


          allocate(density_kernel(tmbmix%orbs%norb,tmbmix%orbs%norb), stat=istat)
          call memocc(istat, density_kernel, 'density_kernel', subname)

          ! Calculate the coefficients
          if(.not.lscv%withder) then
              call get_coeff(iproc,nproc,scf_mode,tmb%lzd,orbs,at,rxyz,denspot,GPU,infoCoeff,energs%ebs,nlpspd,proj,&
                   tmbmix%wfnmd%bpo%blocksize_pdsyev,tmbder%wfnmd%bpo%nproc_pdsyev,&
                   hx,hy,hz,input%SIC,tmbmix,tmb,pnrm,density_kernel,overlapmatrix,calculate_overlap_matrix,&
                   tmblarge, lhphilarge, lhphilargeold, lphilargeold, ldiis_coeff)
          else
              call get_coeff(iproc,nproc,scf_mode,tmb%lzd,orbs,at,rxyz,denspot,GPU,infoCoeff,energs%ebs,nlpspd,proj,&
                   tmbmix%wfnmd%bpo%blocksize_pdsyev,tmbder%wfnmd%bpo%nproc_pdsyev,&
                   hx,hy,hz,input%SIC,tmbmix,tmb,pnrm,density_kernel,overlapmatrix,calculate_overlap_matrix,&
                   tmblargeder, lhphilargeder, lhphilargeoldder, lphilargeoldder, ldiis_coeff)
          end if



          ! Calculate the total energy.
          energy=energs%ebs-energs%eh+energs%exc-energs%evxc-energs%eexctX+energs%eion+energs%edisp
          !write(34,*) energy,energs%ebs,energs%eh,energs%exc,energs%evxc,energs%eexctX,energs%eion,energs%edisp
          energyDiff=energy-energyold
          energyold=energy
!DEBUG
if(iproc==0)then
print *,'ebs,eh,exc,evxc,eexctX,eion,edisp',energs%ebs,energs%eh,energs%exc,energs%evxc,energs%eexctX,energs%eion,energs%edisp
end if
!END DEBUG


          ! Calculate the charge density.
          call sumrhoForLocalizedBasis2(iproc, nproc, &
               tmb%lzd, input, hx, hy ,hz, tmbmix%orbs, tmbmix%comsr, &
               density_kernel, Glr%d%n1i*Glr%d%n2i*denspot%dpbox%n3d, &
               denspot%rhov, at, denspot%dpbox%nscatterarr)

          iall = -product(shape(density_kernel))*kind(density_kernel)
          deallocate(density_kernel,stat=istat)
          call memocc(istat,iall,'density_kernel',subname)

          ! Mix the density.
          !if(trim(input%lin%mixingMethod)=='dens') then
          if(input%lin%scf_mode==LINEAR_MIXDENS_SIMPLE) then
           lscv%compare_outer_loop = pnrm<lscv%self_consistent .or. it_scc==lscv%nit_scc
           call mix_main(iproc, nproc, lscv%mix_hist, lscv%compare_outer_loop, input, glr, lscv%alpha_mix, &
                denspot, mixdiis, rhopotold, rhopotold_out, pnrm, lscv%pnrm_out)
          end if

          ! Calculate the new potential.
          if(iproc==0) write(*,'(1x,a)') '---------------------------------------------------------------- Updating potential.'
          call updatePotential(input%ixc,input%nspin,denspot,energs%eh,energs%exc,energs%evxc)

          ! Mix the potential
          !if(trim(input%lin%mixingMethod)=='pot') then
          if(input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE) then
           lscv%compare_outer_loop = pnrm<lscv%self_consistent .or. it_scc==lscv%nit_scc
           call mix_main(iproc, nproc, lscv%mix_hist, lscv%compare_outer_loop, input, glr, lscv%alpha_mix, &
                denspot, mixdiis, rhopotold, rhopotold_out, pnrm, lscv%pnrm_out)
          end if


          ! Make sure that the previous communication is complete (only do that if this check
          ! for completeness has not been done in get_coeff)
          if(tmbmix%wfnmd%bs%use_derivative_basis .and. .not.tmb%wfnmd%bs%update_phi) then
              call wait_p2p_communication(iproc, nproc, tmb%comgp)
          end if
          if(lscv%withder) then
              call wait_p2p_communication(iproc, nproc, tmbder%comgp)
          end if

          ! Write some informations.
          call printSummary(iproc, it_scc, lscv%info_basis_functions, &
               infoCoeff, pnrm, energy, energyDiff, input%lin%scf_mode)
          if(pnrm<lscv%self_consistent) then
              info_scf=it_scc
              lscv%reduce_convergence_tolerance=.true.
              exit
          else
              info_scf=-1
              lscv%reduce_convergence_tolerance=.false.
          end if

          if(it_scc<lscv%nit_scc_when_optimizing) then
              ! Deallocate the transposed TMBs
              if(tmbmix%can_use_transposed) then
                  iall=-product(shape(tmbmix%psit_c))*kind(tmbmix%psit_c)
                  deallocate(tmbmix%psit_c, stat=istat)
                  call memocc(istat, iall, 'tmbmix%psit_c', subname)
                  iall=-product(shape(tmbmix%psit_f))*kind(tmbmix%psit_f)
                  deallocate(tmbmix%psit_f, stat=istat)
                  call memocc(istat, iall, 'tmbmix%psit_f', subname)
              end if
              write(*,*) 'deallocating overlapmatrix'
              iall=-product(shape(overlapmatrix))*kind(overlapmatrix)
              deallocate(overlapmatrix, stat=istat)
              call memocc(istat, iall, 'overlapmatrix', subname)
          end if

      end do



    call destroy_new_locregs(iproc, nproc, tmblarge)
    call deallocate_auxiliary_basis_function(subname, tmblarge%psi, lhphilarge, lhphilargeold, lphilargeold)
    if(tmblarge%can_use_transposed) then
        iall=-product(shape(tmblarge%psit_c))*kind(tmblarge%psit_c)
        deallocate(tmblarge%psit_c, stat=istat)
        call memocc(istat, iall, 'tmblarge%psit_c', subname)
        iall=-product(shape(tmblarge%psit_f))*kind(tmblarge%psit_f)
        deallocate(tmblarge%psit_f, stat=istat)
        call memocc(istat, iall, 'tmblarge%psit_f', subname)
    end if
    deallocate(tmblarge%confdatarr, stat=istat)


    if(lscv%withder) then
        call destroy_new_locregs(iproc, nproc, tmblargeder)
        call deallocate_auxiliary_basis_function(subname, tmblargeder%psi, lhphilargeder, lhphilargeoldder, lphilargeoldder)
        if(tmblargeder%can_use_transposed) then
            iall=-product(shape(tmblargeder%psit_c))*kind(tmblargeder%psit_c)
            deallocate(tmblargeder%psit_c, stat=istat)
            call memocc(istat, iall, 'tmblargeder%psit_c', subname)
            iall=-product(shape(tmblargeder%psit_f))*kind(tmblargeder%psit_f)
            deallocate(tmblargeder%psit_f, stat=istat)
            call memocc(istat, iall, 'tmblargeder%psit_f', subname)
        end if
        deallocate(tmblargeder%confdatarr, stat=istat)
    end if

    

      iall=-product(shape(locregCenter))*kind(locregCenter)
      deallocate(locregCenter, stat=istat)
      call memocc(istat, iall, 'locregCenter', subname)
      iall=-product(shape(locrad_tmp))*kind(locrad_tmp)
      deallocate(locrad_tmp, stat=istat)
      call memocc(istat, iall, 'locrad_tmp', subname)



      call deallocateDIIS(ldiis_coeff)


!! call set_optimization_variables(input, at, tmb%orbs, tmb%lzd%nlr, tmb%orbs%onwhichatom, &
!!      tmb%confdatarr, tmb%wfnmd, lscv)
!!
!!
!!
!!
!!write(*,*) 'allocated(overlapmatrix)',allocated(overlapmatrix)
!!      write(*,*) 'tmb%confdatarr(1)%prefac', tmb%confdatarr(1)%prefac
!!!!call set_optimization_variables(input, at, tmb%orbs, tmb%lzd%nlr, tmb%orbs%onwhichatom, &
!!!!     tmb%confdatarr, tmb%wfnmd, lscv)
      iall=-product(shape(overlapmatrix))*kind(overlapmatrix)
      deallocate(overlapmatrix, stat=istat)
      call memocc(istat, iall, 'overlapmatrix', subname)


! TEST
      !!do iorb=1,tmb%orbs%norbp
      !!    !!ilr=tmb%orbs%inwhichlocreg(tmb%orbs%isorb+iorb)
      !!    !!iiat=tmb%orbs%onwhichatom(tmb%orbs%isorb+iorb)
      !!!!    tmb%confdatarr(iorb)%prefac=input%lin%potentialPrefac_lowaccuracy(at%iatype(iiat))
      !!    !tmb%confdatarr(iorb)%prefac=1.d0
      !!    write(*,*) 'tmb%confdatarr(iorb)%prefac', tmb%confdatarr(iorb)%prefac
      !!end do
      !!tmb%wfnmd%bs%target_function=TARGET_FUNCTION_IS_TRACE
      !!tmb%wfnmd%bs%nit_basis_optimization=input%lin%nItBasis_lowaccuracy
      !!tmb%wfnmd%bs%conv_crit=input%lin%convCrit_lowaccuracy
      !!lscv%nit_scc=input%lin%nitSCCWhenOptimizing_lowaccuracy+input%lin%nitSCCWhenFixed_lowaccuracy
      !!lscv%nit_scc_when_optimizing=input%lin%nitSCCWhenOptimizing_lowaccuracy
      !!lscv%mix_hist=input%lin%mixHist_lowaccuracy
      !!do ilr=1,tmb%lzd%nlr
      !!    lscv%locrad(ilr)=input%lin%locrad_lowaccuracy(ilr)
      !!end do
      !!if(tmb%wfnmd%bs%update_phi) then
      !!    lscv%alpha_mix=input%lin%alphaMixWhenOptimizing_lowaccuracy
      !!else
      !!    lscv%alpha_mix=input%lin%alphaMixWhenFixed_lowaccuracy
      !!end if




!      ! Deallocate the transposed TMBs
!      write(*,*) 'associated(tmbmix%psit_c)',associated(tmbmix%psit_c)
!      write(*,*) 'associated(tmbmix%psit_f)',associated(tmbmix%psit_f)
!      write(*,*) 'tmbmix%can_use_transposed',tmbmix%can_use_transposed
!      write(*,*) 'associated(tmb%confdatarr)',associated(tmb%confdatarr)
      if(tmbmix%can_use_transposed) then
          iall=-product(shape(tmbmix%psit_c))*kind(tmbmix%psit_c)
          deallocate(tmbmix%psit_c, stat=istat)
          call memocc(istat, iall, 'tmbmix%psit_c', subname)
          iall=-product(shape(tmbmix%psit_f))*kind(tmbmix%psit_f)
          deallocate(tmbmix%psit_f, stat=istat)
          call memocc(istat, iall, 'tmbmix%psit_f', subname)
      end if
  !write(*,*) 'allocated(overlapmatrix)',allocated(overlapmatrix)
!!call set_optimization_variables(input, at, tmb%orbs, tmb%lzd%nlr, tmb%orbs%onwhichatom, &
!!     tmb%confdatarr, tmb%wfnmd, lscv)




      ! Print out values related to two iterations of the outer loop.
      if(iproc==0) then
          write(*,'(3x,a,7es18.10)') 'ebs, ehart, eexcu, vexcu, eexctX, eion, edisp', &
              energs%ebs, energs%eh, energs%exc, energs%evxc, energs%eexctX, energs%eion, energs%edisp
          !if(trim(input%lin%mixingMethod)=='dens') then
          if(input%lin%scf_mode==LINEAR_MIXDENS_SIMPLE) then
             if (.not. lscv%lowaccur_converged) then
                 write(*,'(3x,a,3x,i0,es11.2,es27.17,es14.4)')&
                      'itoutL, Delta DENSOUT, energy, energyDiff', itout, lscv%pnrm_out, energy, &
                      energy-energyoldout
             else
                 write(*,'(3x,a,3x,i0,es11.2,es27.17,es14.4)')&
                      'itoutH, Delta DENSOUT, energy, energyDiff', itout, lscv%pnrm_out, energy, &
                      energy-energyoldout
             end if
          else if(input%lin%scf_mode==LINEAR_MIXPOT_SIMPLE) then
             if (.not. lscv%lowaccur_converged) then
                 write(*,'(3x,a,3x,i0,es11.2,es27.17,es14.4)')&
                      'itoutH, Delta POTOUT, energy energyDiff', itout, lscv%pnrm_out, energy, energy-energyoldout
             else
                 write(*,'(3x,a,3x,i0,es11.2,es27.17,es14.4)')&
                      'itoutL, Delta POTOUT, energy energyDiff', itout, lscv%pnrm_out, energy, energy-energyoldout
             end if
          else if(input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) then
             if (.not. lscv%lowaccur_converged) then
                 write(*,'(3x,a,3x,i0,es11.2,es27.17,es14.4)')&
                      'itoutH, fnrm coeff, energy energyDiff', itout, lscv%pnrm_out, energy, energy-energyoldout
             else
                 write(*,'(3x,a,3x,i0,es11.2,es27.17,es14.4)')&
                      'itoutL, fnrm coeff, energy energyDiff', itout, lscv%pnrm_out, energy, energy-energyoldout
             end if
          end if
      end if
      call print_info(iproc, itout, lscv%info_basis_functions, info_scf, input%lin%scf_mode, tmb%wfnmd%bs%target_function, &
           fnrm_tmb, pnrm, trace, energy, energy-energyoldout)


      energyoldout=energy

      !!! Deallocate DIIS structures.
      !!call deallocateDIIS(ldiis)


      call check_for_exit(input, lscv)
      if(lscv%exit_outer_loop) exit outerLoop




  end do outerLoop

  ! Deallocate DIIS structures.
  call deallocateDIIS(ldiis)

  call deallocateCommunicationbufferSumrho(tmb%comsr, subname)
  call deallocateCommunicationbufferSumrho(tmbder%comsr, subname)

  call wait_p2p_communication(iproc, nproc, tmb%comgp)
  call deallocateCommunicationsBuffersPotential(tmb%comgp, subname)
  if(tmbder%wfnmd%bs%use_derivative_basis) then
     call wait_p2p_communication(iproc, nproc, tmbder%comgp)
     call deallocateCommunicationsBuffersPotential(tmbder%comgp, subname)
  end if

  iall=-product(shape(rhopotold_out))*kind(rhopotold_out)
  deallocate(rhopotold_out, stat=istat)
  call memocc(istat, iall, 'rhopotold_out', subname)

  if(input%lin%mixHist_highaccuracy>0) then
      call deallocateMixrhopotDIIS(mixdiis)
  end if

  !Write the linear wavefunctions to file if asked
  if(input%lin%plotBasisFunctions /= WF_FORMAT_NONE) then
    call writemywaves_linear(iproc,trim(input%dir_output) // 'minBasis',input%lin%plotBasisFunctions,tmb%Lzd,&
       tmbmix%orbs,orbs%norb,hx,hy,hz,at,rxyz,tmbmix%psi,tmbmix%wfnmd%coeff)
   end if

  ! Allocate the communication buffers for the calculation of the charge density.
  call allocateCommunicationbufferSumrho(iproc, tmbmix%comsr, subname)
  call communicate_basis_for_density(iproc, nproc, tmb%lzd, tmbmix%orbs, tmbmix%psi, tmbmix%comsr)
  allocate(density_kernel(tmbmix%orbs%norb,tmbmix%orbs%norb), stat=istat)
  call memocc(istat, density_kernel, 'density_kernel', subname)
  call calculate_density_kernel(iproc, nproc, tmbmix%orbs%norb, orbs%norb, orbs%norbp, orbs%isorb, &
       tmbmix%wfnmd%ld_coeff, tmbmix%wfnmd%coeff, density_kernel)
  call sumrhoForLocalizedBasis2(iproc, nproc, tmb%lzd, input, hx, hy, hz, &
       tmbmix%orbs, tmbmix%comsr, density_kernel, Glr%d%n1i*Glr%d%n2i*denspot%dpbox%n3d, &
       denspot%rhov, at,denspot%dpbox%nscatterarr)
  iall = -product(shape(density_kernel))*kind(density_kernel)
  deallocate(density_kernel,stat=istat)
  call memocc(istat,iall,'density_kernel',subname)

  call deallocateCommunicationbufferSumrho(tmbmix%comsr, subname)

  ! Build global orbitals psi (the physical ones).
  if(input%lin%transformToGlobal) then
     if(nproc>1) then
        allocate(psit(max(orbs%npsidim_orbs,orbs%npsidim_comp)), stat=istat)
        call memocc(istat, psit, 'psit', subname)
     else
        psit => psi
     end if
     call transformToGlobal(iproc, nproc, tmb%lzd, tmbmix%orbs, orbs, comms, input, tmbmix%wfnmd%ld_coeff, &
          tmbmix%wfnmd%coeff, tmbmix%psi, psi, psit)
     if(nproc>1) then
        iall=-product(shape(psit))*kind(psit)
        deallocate(psit, stat=istat)
        call memocc(istat, iall, 'psit', subname)
     else
        nullify(psit)
     end if
  end if


  nullify(rho,pot)

  iall=-product(shape(lscv%locrad))*kind(lscv%locrad)
  deallocate(lscv%locrad, stat=istat)
  call memocc(istat, iall, 'lscv%locrad', subname)

  iall=-product(shape(eval))*kind(eval)
  deallocate(eval, stat=istat)
  call memocc(istat, iall, 'eval', subname)

  call timing(iproc,'WFN_OPT','PR')

end subroutine linearScaling





subroutine printSummary(iproc, itSCC, infoBasisFunctions, infoCoeff, pnrm, energy, energyDiff, scf_mode)
use module_base
use module_types
!
! Purpose:
! ========
!   Print a short summary of some values calculated during the last iteration in the self
!   consistency cycle.
! 
! Calling arguments:
! ==================
!   Input arguments
!   ---------------
!
implicit none

! Calling arguments
integer,intent(in):: iproc, itSCC, infoBasisFunctions, infoCoeff, scf_mode
real(8),intent(in):: pnrm, energy, energyDiff

  if(iproc==0) then
      write(*,'(1x,a)') repeat('+',92 + int(log(real(itSCC))/log(10.)))
      write(*,'(1x,a,i0,a)') 'at iteration ', itSCC, ' of the density optimization:'
      !!if(infoCoeff<0) then
      !!    write(*,'(3x,a)') '- WARNING: coefficients not converged!'
      !!else if(infoCoeff>0) then
      !!    write(*,'(3x,a,i0,a)') '- coefficients converged in ', infoCoeff, ' iterations.'
      if(scf_mode==LINEAR_DIRECT_MINIMIZATION) then
          write(*,'(3x,a)') 'coefficients obtained by direct minimization.'
      else
          write(*,'(3x,a)') 'coefficients obtained by diagonalization.'
      end if
      !!end if
      !if(mixingMethod=='dens') then
      if(scf_mode==LINEAR_MIXDENS_SIMPLE) then
          write(*,'(3x,a,3x,i0,es11.2,es27.17,es14.4)') 'it, Delta DENS, energy, energyDiff', itSCC, pnrm, energy, energyDiff
      !else if(mixingMethod=='pot') then
      else if(scf_mode==LINEAR_MIXPOT_SIMPLE) then
          write(*,'(3x,a,3x,i0,es11.2,es27.17,es14.4)') 'it, Delta POT, energy, energyDiff', itSCC, pnrm, energy, energyDiff
      else if(scf_mode==LINEAR_DIRECT_MINIMIZATION) then
          write(*,'(3x,a,3x,i0,es11.2,es27.17,es14.4)') 'it, fnrm coeff, energy, energyDiff', itSCC, pnrm, energy, energyDiff
      end if
      write(*,'(1x,a)') repeat('+',92 + int(log(real(itSCC))/log(10.)))
  end if

end subroutine printSummary



subroutine print_info(iproc, itout, info_tmb, info_coeff, scf_mode, target_function, &
           fnrm_tmb, pnrm, value_tmb, energy, energyDiff)
use module_base
use module_types
!
! Purpose:
! ========
!   Print a short summary of some values calculated during the last iteration in the self
!   consistency cycle.
! 
! Calling arguments:
! ==================
!   Input arguments
!   ---------------
!
implicit none

! Calling arguments
integer,intent(in):: iproc, itout, info_tmb, info_coeff, scf_mode, target_function
real(8),intent(in):: fnrm_tmb, pnrm, value_tmb, energy, energyDiff

  if(iproc==0) then
      write(*,'(1x,a)') repeat('#',92 + int(log(real(itout))/log(10.)))
      write(*,'(1x,a,i0,a)') 'at iteration ', itout, ' of the outer loop:'
      write(*,'(3x,a)') '> basis functions optimization:'
      if(target_function==TARGET_FUNCTION_IS_TRACE) then
          write(*,'(5x,a)') '- target function is trace'
      else if(target_function==TARGET_FUNCTION_IS_ENERGY) then
          write(*,'(5x,a)') '- target function is energy'
      end if
      if(info_tmb<0) then
          write(*,'(5x,a)') '- WARNING: basis functions not converged!'
      else
          write(*,'(5x,a,i0,a)') '- basis functions converged in ', info_tmb, ' iterations.'
      end if
      write(*,'(5x,a,es15.6,2x,es10.2)') 'Final values: target function, fnrm', value_tmb, fnrm_tmb
      write(*,'(3x,a)') '> density optimization:'
      if(scf_mode==LINEAR_DIRECT_MINIMIZATION) then
          write(*,'(5x,a)') '- using direct minimization.'
      else
          write(*,'(5x,a)') '- using diagonalization / mixing.'
      end if
      if(info_coeff<0) then
          write(*,'(5x,a)') '- WARNING: density optimization not converged!'
      else
          write(*,'(5x,a,i0,a)') '- density optimization converged in ', info_coeff, ' iterations.'
      end if
      if(scf_mode==LINEAR_MIXDENS_SIMPLE) then
          write(*,'(5x,a,3x,i0,es12.2,es27.17)') 'FINAL values: it, Delta DENS, energy', itout, pnrm, energy
      else if(scf_mode==LINEAR_MIXPOT_SIMPLE) then
          write(*,'(5x,a,3x,i0,es12.2,es27.17)') 'FINAL values: it, Delta POT, energy', itout, pnrm, energy
      else if(scf_mode==LINEAR_DIRECT_MINIMIZATION) then
          write(*,'(5x,a,3x,i0,es12.2,es27.17)') 'FINAL values: it, fnrm coeff, energy', itout, pnrm, energy
      end if
      write(*,'(3x,a,es14.6)') '> energy difference to last iteration:', energyDiff
      write(*,'(1x,a)') repeat('#',92 + int(log(real(itout))/log(10.)))
  end if

end subroutine print_info



subroutine transformToGlobal(iproc,nproc,lzd,lorbs,orbs,comms,input,ld_coeff,coeff,lphi,psi,psit)
use module_base
use module_types
use module_interfaces, exceptThisOne => transformToGlobal
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, ld_coeff
type(local_zone_descriptors),intent(in):: lzd
type(orbitals_data),intent(in):: lorbs, orbs
type(communications_arrays):: comms
type(input_variables),intent(in):: input
real(8),dimension(ld_coeff,orbs%norb),intent(in):: coeff
real(8),dimension(lorbs%npsidim_orbs),intent(inout):: lphi
real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),target,intent(out):: psi
real(8),dimension(:),pointer,intent(out):: psit

! Local variables
integer:: ind1, ind2, istat, iall, iorb, ilr, ldim, gdim, nvctrp
real(8),dimension(:),pointer:: phiWork
real(8),dimension(:),allocatable:: phi
character(len=*),parameter:: subname='transformToGlobal'
type(orbitals_data):: gorbs
type(communications_arrays):: gcomms


  call nullify_orbitals_data(gorbs)
  call copy_orbitals_data(lorbs, gorbs, subname)
  call orbitals_communicators(iproc,nproc,lzd%glr,gorbs,gcomms)

  allocate(phi(max(gorbs%npsidim_orbs,gorbs%npsidim_comp)+ndebug), stat=istat)
  call memocc(istat, phi, 'phi', subname)
  allocate(phiWork(max(size(phi),size(psi))), stat=istat)
  call memocc(istat, phiWork, 'phiWork', subname)

  ind1=1
  ind2=1
  if (max(gorbs%npsidim_orbs,gorbs%npsidim_comp) > 0) &
       call to_zero(max(gorbs%npsidim_orbs,gorbs%npsidim_comp),phi(1))

  do iorb=1,lorbs%norbp
      ilr = lorbs%inWhichLocreg(lorbs%isorb+iorb)
      ldim=lzd%Llr(ilr)%wfd%nvctr_c+7*lzd%Llr(ilr)%wfd%nvctr_f
      gdim=lzd%Glr%wfd%nvctr_c+7*lzd%Glr%wfd%nvctr_f
      call Lpsi_to_global2(iproc,nproc,ldim,gdim,lorbs%norb,lorbs%nspinor,input%nspin,lzd%Glr,&
           lzd%Llr(ilr),lphi(ind2),phi(ind1))
      ind1=ind1+lzd%Glr%wfd%nvctr_c+7*lzd%Glr%wfd%nvctr_f
      ind2=ind2+lzd%Llr(ilr)%wfd%nvctr_c+7*lzd%Llr(ilr)%wfd%nvctr_f
  end do
  call transpose_v(iproc, nproc, lorbs, lzd%Glr%wfd, gcomms, phi, work=phiWork)


  if(iproc==0) then
      write(*,'(1x,a)', advance='no') '------------------------------------- Building linear combinations... '
  end if
  ! Build the extended orbital psi as a linear combination of localized basis functions phi. for real O(N)
  ! this has to replaced, but at the moment it is still needed.
  !call buildWavefunctionModified(iproc, nproc, orbs, gorbs, comms, gcomms, phi, psi, coeff)
  nvctrp=sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor
  call dgemm('n', 'n', nvctrp, orbs%norb, lorbs%norb, 1.d0, phi(1), nvctrp, coeff(1,1), &
       lorbs%norb, 0.d0, psi(1), nvctrp)


  if(nproc>1) then
      call dcopy(orbs%npsidim_comp, psi, 1, psit, 1)
  else
      psit => psi
  end if

  call untranspose_v(iproc, nproc, orbs, lzd%Glr%wfd, comms, psi, work=phiWork)

  if(iproc==0) write(*,'(a)') 'done.'


  iall=-product(shape(phi))*kind(phi)
  deallocate(phi, stat=istat)
  call memocc(istat, iall, 'phi', subname)
  iall=-product(shape(phiWork))*kind(phiWork)
  deallocate(phiWork, stat=istat)
  call memocc(istat, iall, 'phiWork', subname)

  call deallocate_orbitals_data(gorbs, subname)
  call deallocate_communications_arrays(gcomms, subname)

end subroutine transformToGlobal



subroutine set_optimization_variables(input, at, lorbs, nlr, onwhichatom, confdatarr, wfnmd, lscv)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: nlr
  type(orbitals_data),intent(in):: lorbs
  type(input_variables),intent(in):: input
  type(atoms_data),intent(in):: at
  integer,dimension(lorbs%norb),intent(in):: onwhichatom
  type(confpot_data),dimension(lorbs%norbp),intent(inout):: confdatarr
  type(wfn_metadata),intent(inout):: wfnmd
  type(linear_scaling_control_variables),intent(inout):: lscv

  ! Local variables
  integer:: iorb, ilr, iiat

  if(lscv%lowaccur_converged) then

      do iorb=1,lorbs%norbp
          ilr=lorbs%inwhichlocreg(lorbs%isorb+iorb)
          iiat=onwhichatom(lorbs%isorb+iorb)
          confdatarr(iorb)%prefac=input%lin%potentialPrefac_highaccuracy(at%iatype(iiat))
      end do
      wfnmd%bs%target_function=TARGET_FUNCTION_IS_ENERGY
      wfnmd%bs%nit_basis_optimization=input%lin%nItBasis_highaccuracy
      wfnmd%bs%conv_crit=input%lin%convCrit_highaccuracy
      lscv%nit_scc=input%lin%nitSCCWhenOptimizing_highaccuracy+input%lin%nitSCCWhenFixed_highaccuracy
      lscv%nit_scc_when_optimizing=input%lin%nitSCCWhenOptimizing_highaccuracy
      lscv%mix_hist=input%lin%mixHist_highaccuracy
      do ilr=1,nlr
          lscv%locrad(ilr)=input%lin%locrad_highaccuracy(ilr)
      end do
      if(wfnmd%bs%update_phi) then
          lscv%alpha_mix=input%lin%alphaMixWhenOptimizing_highaccuracy
      else
          lscv%alpha_mix=input%lin%alphaMixWhenFixed_highaccuracy
      end if
      !!if(.not.lscv%enlarge_locreg) lscv%enlarge_locreg=.true.


  else

      do iorb=1,lorbs%norbp
          ilr=lorbs%inwhichlocreg(lorbs%isorb+iorb)
          iiat=onwhichatom(lorbs%isorb+iorb)
          confdatarr(iorb)%prefac=input%lin%potentialPrefac_lowaccuracy(at%iatype(iiat))
      end do
      wfnmd%bs%target_function=TARGET_FUNCTION_IS_TRACE
      wfnmd%bs%nit_basis_optimization=input%lin%nItBasis_lowaccuracy
      wfnmd%bs%conv_crit=input%lin%convCrit_lowaccuracy
      lscv%nit_scc=input%lin%nitSCCWhenOptimizing_lowaccuracy+input%lin%nitSCCWhenFixed_lowaccuracy
      lscv%nit_scc_when_optimizing=input%lin%nitSCCWhenOptimizing_lowaccuracy
      lscv%mix_hist=input%lin%mixHist_lowaccuracy
      do ilr=1,nlr
          lscv%locrad(ilr)=input%lin%locrad_lowaccuracy(ilr)
      end do
      if(wfnmd%bs%update_phi) then
          lscv%alpha_mix=input%lin%alphaMixWhenOptimizing_lowaccuracy
      else
          lscv%alpha_mix=input%lin%alphaMixWhenFixed_lowaccuracy
      end if

  end if

end subroutine set_optimization_variables



subroutine adjust_locregs_and_confinement(iproc, nproc, hx, hy, hz, &
           input, tmb, tmbder, denspot, ldiis, lscv)
  use module_base
  use module_types
  use module_interfaces, except_this_one => adjust_locregs_and_confinement
  implicit none
  
  ! Calling argument
  integer,intent(in):: iproc, nproc
  real(8),intent(in):: hx, hy, hz
  type(input_variables),intent(in):: input
  type(DFT_wavefunction),intent(inout):: tmb, tmbder
  type(DFT_local_fields),intent(inout) :: denspot
  type(localizedDIISParameters),intent(inout):: ldiis
  type(linear_scaling_control_variables),intent(inout):: lscv

  ! Local variables
  integer :: ilr
  logical:: redefine_derivatives, change
  character(len=*),parameter:: subname='adjust_locregs_and_confinement'

  if(tmb%wfnmd%bs%confinement_decrease_mode==DECREASE_ABRUPT) then
      lscv%decrease_factor_total=1.d0
  else if(tmb%wfnmd%bs%confinement_decrease_mode==DECREASE_LINEAR) then
      if(lscv%info_basis_functions>0) then
          lscv%idecrease=lscv%idecrease+1
          lscv%ifail=0
      else
          lscv%ifail=lscv%ifail+1
      end if
      lscv%decrease_factor_total=1.d0-dble(lscv%idecrease)*input%lin%decrease_step
  end if
  if(tmbder%wfnmd%bs%target_function==TARGET_FUNCTION_IS_ENERGY) lscv%decrease_factor_total=1.d0
  if(iproc==0) write(*,'(1x,a,f6.2,a)') 'Reduce the confining potential to ', &
      100.d0*lscv%decrease_factor_total,'% of its initial value.'
  tmb%confdatarr(:)%prefac=lscv%decrease_factor_total*tmb%confdatarr(:)%prefac


  lscv%locreg_increased=.false.
  redefine_derivatives=.false.
  if(lscv%ifail>=input%lin%increase_locrad_after .and. .not.lscv%lowaccur_converged .and. &
     input%lin%locrad_increase_amount>0.d0) then
      lscv%increase_locreg=lscv%increase_locreg+input%lin%locrad_increase_amount
      !lscv%increase_locreg=lscv%increase_locreg+0.d0
      if(iproc==0) then
          write(*,'(1x,a)') 'It seems that the convergence criterion can not be reached with this localization radius.'
          write(*,'(1x,a,f6.2)') 'The localization radius is increased by totally',lscv%increase_locreg
      end if
      lscv%ifail=0
      lscv%locrad=lscv%locrad+lscv%increase_locreg
      if(lscv%withder) then
          redefine_derivatives=.true.
      end if
      lscv%locreg_increased=.true.
  end if

  !redefine_derivatives=.false.
  if(lscv%lowaccur_converged .and. lscv%enlarge_locreg) then
      change = .false.
      do ilr = 1, tmb%lzd%nlr
         if(input%lin%locrad_highaccuracy(ilr) /= input%lin%locrad_lowaccuracy(ilr)) then
             change = .true.
             exit
         end if
      end do
      if(change) then
         if(iproc==0) then
             write(*,'(1x,a)') 'Increasing the localization radius for the high accuracy part.'
         end if
         lscv%locreg_increased=.true.
      end if
      lscv%enlarge_locreg=.false. !flag to indicate that the locregs should not be increased any more in the following iterations
  end if
  if(lscv%locreg_increased) then
      call redefine_locregs_quantities(iproc, nproc, hx, hy, hz, lscv%locrad, .true., tmb%lzd, tmb, tmb, denspot, ldiis)
  end if
  if(redefine_derivatives) then
      call redefine_locregs_quantities(iproc, nproc, hx, hy, hz, tmb%lzd%llr(:)%locrad, .false., tmb%lzd, tmb, tmbder, denspot)
  end if

end subroutine adjust_locregs_and_confinement



subroutine adjust_DIIS_for_high_accuracy(input, tmb, denspot, ldiis, mixdiis, lscv)
  use module_base
  use module_types
  use module_interfaces, except_this_one => adjust_DIIS_for_high_accuracy
  implicit none
  
  ! Calling arguments
  type(input_variables),intent(in):: input
  type(DFT_wavefunction),intent(in):: tmb
  type(DFT_local_fields),intent(inout) :: denspot
  type(localizedDIISParameters),intent(inout):: ldiis
  type(mixrhopotDIISParameters),intent(inout):: mixdiis
  type(linear_scaling_control_variables),intent(inout):: lscv
  
  !!lscv%exit_outer_loop=.false.
  
  if(lscv%lowaccur_converged) then
      !!lscv%nit_highaccuracy=lscv%nit_highaccuracy+1
      if(lscv%nit_highaccuracy==input%lin%nit_highaccuracy+1) then
          ! Deallocate DIIS structures.
          !!call deallocateDIIS(ldiis)
          !!lscv%exit_outer_loop=.true.
      end if
      ! only use steepest descent if the localization regions may change
      if(input%lin%nItInnerLoop/=-1 .or. tmb%wfnmd%bs%locreg_enlargement/=1.d0) then
          ldiis%isx=0
      end if
  
      if(input%lin%mixHist_lowaccuracy==0 .and. input%lin%mixHist_highaccuracy>0) then
          call initializeMixrhopotDIIS(input%lin%mixHist_highaccuracy, denspot%dpbox%ndimpot, mixdiis)
      else if(input%lin%mixHist_lowaccuracy>0 .and. input%lin%mixHist_highaccuracy==0) then
          call deallocateMixrhopotDIIS(mixdiis)
      end if
  end if
  
end subroutine adjust_DIIS_for_high_accuracy


subroutine check_for_exit(input, lscv)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(input_variables),intent(in):: input
  type(linear_scaling_control_variables),intent(inout):: lscv

  lscv%exit_outer_loop=.false.
  
  if(lscv%lowaccur_converged) then
      lscv%nit_highaccuracy=lscv%nit_highaccuracy+1
      if(lscv%nit_highaccuracy==input%lin%nit_highaccuracy) then
          lscv%exit_outer_loop=.true.
      else if (lscv%pnrm_out<input%lin%highaccuracy_converged) then !lr408
          lscv%exit_outer_loop=.true.
      end if
  end if

end subroutine check_for_exit


function check_whether_derivatives_to_be_used(input, itout, lscv)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  type(input_variables),intent(in):: input
  integer,intent(in):: itout
  type(linear_scaling_control_variables),intent(in):: lscv
  logical:: check_whether_derivatives_to_be_used

  ! Local variables
  logical:: withder

  if(input%lin%mixedmode) then
      if( (.not.lscv%lowaccur_converged .and. &
           (itout==input%lin%nit_lowaccuracy+1 .or. lscv%pnrm_out<input%lin%lowaccuray_converged) ) &
          .or. lscv%lowaccur_converged ) then
          withder=.true.
      else
          withder=.false.
      end if
  else
      if(input%lin%useDerivativeBasisFunctions) then
          withder=.true.
      else
          withder=.false.
      end if
  end if
  check_whether_derivatives_to_be_used=withder

end function check_whether_derivatives_to_be_used



subroutine check_whether_lowaccuracy_converged(itout, input, lscv)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: itout
  type(input_variables),intent(in):: input
  type(linear_scaling_control_variables),intent(inout):: lscv
  
  if(.not.lscv%lowaccur_converged .and. &
     (itout==input%lin%nit_lowaccuracy+1 .or. lscv%pnrm_out<input%lin%lowaccuray_converged .or. &
      lscv%decrease_factor_total<1.d0-input%lin%decrease_amount)) then
      lscv%lowaccur_converged=.true.
      lscv%nit_highaccuracy=0
  end if 

end subroutine check_whether_lowaccuracy_converged
