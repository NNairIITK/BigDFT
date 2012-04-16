!!!!> Application of the Hamiltonian
!!!subroutine HamiltonianApplicationConfinement2(input,iproc,nproc,at,Lzd,orbs,lin,hx,hy,hz,rxyz,&
!!!     ngatherarr,ndimpot,pot,psi,hpsi,&
!!!     ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU,radii_cf, comgp, onWhichAtomp, withConfinement, energyReductionFlag, &
!!!     doNotCalculate, pkernel,orbsocc,psirocc)
!!!  use module_base
!!!  use module_types
!!!  use libxc_functionals
!!!  use module_interfaces, exceptThisOne => HamiltonianApplicationConfinement2
!!!  implicit none
!!!  integer, intent(in) :: iproc,nproc,nspin,ndimpot
!!!  real(gp), intent(in) :: hx,hy,hz
!!!  type(atoms_data), intent(in) :: at
!!!  type(input_variables), intent(in) :: input
!!!  type(local_zone_descriptors),intent(inout) :: Lzd
!!!  type(orbitals_data),intent(in):: orbs
!!!  type(linearParameters),intent(in):: lin
!!!  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
!!!  real(gp), dimension(3,at%nat), intent(in) :: rxyz
!!!  !real(wp), dimension(lin%Lorbs%npsidim), intent(in) :: psi
!!!  real(wp), dimension(orbs%npsidim), intent(in) :: psi
!!!  real(wp), dimension(max(ndimpot,1)*nspin), intent(in) :: pot
!!!  !real(wp), dimension(:), pointer :: pot
!!!  real(gp), intent(out) :: ekin_sum,epot_sum,eexctX,eproj_sum
!!!  !real(wp), target, dimension(lin%Lorbs%npsidim), intent(out) :: hpsi
!!!  real(wp), target, dimension(orbs%npsidim), intent(out) :: hpsi
!!!  type(GPU_pointers), intent(inout) :: GPU
!!!  real(gp), dimension(at%ntypes,3+ndebug), intent(in) :: radii_cf
!!!  type(p2pCommsgatherPot), intent(in):: comgp
!!!  integer,dimension(orbs%norbp),intent(in):: onWhichAtomp
!!!  logical,intent(in):: withConfinement
!!!  logical,intent(in):: energyReductionFlag
!!!  logical,dimension(lzd%nlr),intent(in),optional:: doNotCalculate
!!!  real(dp), dimension(*), optional :: pkernel
!!!  type(orbitals_data), intent(in), optional :: orbsocc
!!!  real(wp), dimension(:), pointer, optional :: psirocc
!!!  !local variables
!!!  real(gp) :: tmp_ekin_sum,tmp_epot_sum,tmp_eproj_sum
!!!  real(gp), dimension(2,orbs%norbp) :: ekin
!!!  real(gp), dimension(2,orbs%norbp) :: epot
!!!  real(wp), dimension(:), pointer :: hpsi2
!!!  character(len=*), parameter :: subname='LinearHamiltonianApplicationConfinement2'
!!!  logical :: exctX,op2p
!!!  integer :: i_all,i_stat,ierr,iorb,n3p,ispot,istart_c,iat, i3s, i3e, ind1, ind2, ldim, gdim, jlr
!!!  integer :: istart_ck,isorb,ieorb,ikpt,ispsi_k,nspinor,ispsi, localnorb
!!!  integer :: ilr,dimwf,ind,size_Lpot,size_pot
!!!  integer :: tmp_norbp, istorb
!!!  real(dp),dimension(:),pointer:: Lpot
!!!  real(wp),dimension(:),allocatable :: hpsi_proj
!!!!OCL  integer, dimension(3) :: periodic
!!!!OCL  real(wp) :: maxdiff
!!!!OCL  real(gp) :: eproj,ek_fake,ep_fake
!!!  real(gp), dimension(3,2) :: wrkallred
!!!!OCL  real(wp), dimension(:), allocatable :: hpsi_OCL
!!!
!!!  !!do i_all=1,size(psi)
!!!  !!    write(7100+iproc,*) psi(i_all)
!!!  !!    !write(7110+iproc,*) psi(i_all)
!!!  !!end do
!!!
!!!
!!!  !initialise exact exchange energy 
!!!  op2p=(eexctX == -99.0_gp)
!!!  eexctX=0.0_gp
!!!
!!!  exctX = libxc_functionals_exctXfac() /= 0.0_gp
!!!
!!!  ! Allocate the nonlocal descriptors for the locregs
!!!  allocate(Lzd%Lnlpspd(Lzd%nlr),stat=i_stat)   
!!!  do ilr=1,Lzd%nlr
!!!      call nullify_nonlocal_psp_descriptors(Lzd%Lnlpspd(ilr))
!!!  end do
!!!
!!!  !initialize accumulators
!!!  ekin_sum = 0.0_gp
!!!  epot_sum = 0.0_gp
!!!  eproj_sum= 0.0_gp
!!!  ind = 1
!!!  istorb=1
!!!  !write(*,'(a,i6,4x,100i4)') 'iproc, orbs%inWhichLocregp', iproc, orbs%inWhichLocregp
!!!  do ilr= 1, Lzd%nlr
!!!     if(lin%useDerivativeBasisFunctions) then
!!!         ! Set localnorb if we use the derivative basis functions.
!!!         localnorb=0
!!!         do iorb=1,orbs%norbp
!!!             if(orbs%inWhichLocregp(iorb)==ilr) then
!!!                 localnorb = localnorb+1
!!!             end if
!!!         end do
!!!     else
!!!         localnorb=lzd%Llr(ilr)%localnorb
!!!     end if
!!!     !write(*,'(a,3i7)') 'iproc, ilr, localnorb', iproc, ilr, localnorb
!!!     ! Cycle if the process does not have any orbitals belonging
!!!     ! to this localization region.
!!!     !write(*,'(a,3i8)') 'iproc, ilr, Lzd%Llr(ilr)%Localnorb', iproc, ilr, Lzd%Llr(ilr)%Localnorb
!!!     !if(Lzd%Llr(ilr)%Localnorb == 0) then
!!!     if(localnorb == 0) then
!!!         !write(*,'(a,i0,a,i0)') 'process ',iproc,' cycles for ilr=',ilr
!!!         cycle
!!!     end if
!!!
!!!     ! If no calculation for this localization region is required, only increase the index
!!!     if(present(doNotCalculate)) then
!!!         if(doNotCalculate(ilr)) then
!!!             !dimwf=(Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*Lzd%Llr(ilr)%Localnorb*orbs%nspinor*nspin
!!!             dimwf=(Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*localnorb*orbs%nspinor*nspin
!!!             ind = ind + dimwf
!!!             !write(*,'(a,i0,a,i0)') 'process ',iproc,' cycles for locreg ',ilr
!!!             cycle
!!!         end if
!!!     end if
!!!
!!!     allocate(Lpot(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*nspin), stat=i_stat)
!!!     call memocc(i_stat,Lpot,'Lpot',subname)
!!! 
!!! 
!!!     !determine the dimension of the potential array (copied from full_local_potential)
!!!     if (exctX) then
!!!        stop 'exctX not yet implemented!'
!!!        size_pot=Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*nspin + &
!!!         max(max(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*Lzd%Llr(ilr)%Localnorb*nspin,ngatherarr(0,1)*orbs%norb),1) !part which refers to exact exchange
!!!        size_Lpot=Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*nspin + &
!!!           max(max(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*Lzd%Llr(ilr)%Localnorb*nspin,&
!!!           ngatherarr(0,1)*orbs%norb),1) !CHECK THIS...DOES NOT WORK YET
!!!     else
!!!        !size_pot=Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*nspin
!!!        !size_pot=Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*nscatterarr(iproc,2)*nspin
!!!        size_Lpot = Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*nspin
!!!     end if
!!! 
!!!     ! Extract the part of the potential which is needed for the current localization region.
!!!     i3s=lzd%Llr(ilr)%nsi3-comgp%ise3(1,iproc)+2 ! starting index of localized potential with respect to total potential in comgp%recvBuf
!!!     i3e=lzd%Llr(ilr)%nsi3+lzd%Llr(ilr)%d%n3i-comgp%ise3(1,iproc)+1 ! ending index of localized potential with respect to total potential in comgp%recvBuf
!!!     !write(*,'(a,2i4,3i7)') 'iproc, ilr, lzd%Llr(ilr)%nsi3, comgp%ise3(1,iproc), lzd%Llr(ilr)%d%n3i', iproc, ilr, lzd%Llr(ilr)%nsi3, comgp%ise3(1,iproc), lzd%Llr(ilr)%d%n3i
!!!     if(i3e-i3s+1 /= Lzd%Llr(ilr)%d%n3i) then
!!!         write(*,'(a,i0,3x,i0)') 'ERROR: i3e-i3s+1 /= Lzd%Llr(ilr)%d%n3i', i3e-i3s+1, Lzd%Llr(ilr)%d%n3i
!!!         stop
!!!     end if
!!!     !write(*,'(a,i4,3x,2i8,i15)') 'iproc, i3s, i3e, (i3e-i3s+1)*lzd%glr%d%n2i*lzd%glr%d%n1i', iproc, i3s, i3e, (i3e-i3s+1)*lzd%glr%d%n2i*lzd%glr%d%n1i
!!!     call global_to_local_parallel(lzd%Glr, lzd%Llr(ilr), nspin, ndimpot, size_Lpot, pot, Lpot, i3s, i3e)
!!!
!!!     ! Set some quantities: ispot=shift for potential, dimwf=dimension of wavefunction
!!!     ispot=Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*nspin+1
!!!     !dimwf=(Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*Lzd%Llr(ilr)%Localnorb*orbs%nspinor*nspin
!!!     dimwf=(Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*localnorb*orbs%nspinor*nspin
!!!
!!!     ! EXACT EXCHANGE NOT TESTED: SHOULD CHECK IF EVERYTHING IF FINE
!!!     !fill the rest of the potential with the exact-exchange terms
!!!     if (present(pkernel) .and. exctX) then
!!!        stop 'exctx not yet implemented!'
!!!        n3p=ngatherarr(iproc,1)/(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i)
!!!        !exact exchange for virtual orbitals (needs psirocc)
!!!   
!!!        !here we have to add the round part
!!!        if (present(psirocc) .and. present(orbsocc)) then
!!!           call exact_exchange_potential_virt(iproc,nproc,at%geocode,nspin,&
!!!                Lzd%Llr(ilr),orbsocc,orbs,ngatherarr(0,1),n3p,&
!!!                0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,pkernel,psirocc,psi(ind:ind+dimwf-1),Lpot)
!!!           eexctX = 0._gp
!!!        else
!!!   !!$        call exact_exchange_potential_round(iproc,nproc,at%geocode,nspin,lr,orbs,&
!!!   !!$             0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,pkernel,psi,pot(ispot),eexctX)
!!!   
!!!           !here the condition for the scheme should be chosen
!!!           if (.not. op2p) then
!!!              call exact_exchange_potential(iproc,nproc,at%geocode,nspin,&
!!!                   Lzd%Llr(ilr),orbs,ngatherarr(0,1),n3p,&
!!!                   0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,pkernel,psi(ind:ind+dimwf-1),Lpot,eexctX)
!!!           else
!!!              !the psi should be transformed in real space
!!!              call exact_exchange_potential_round(iproc,nproc,at%geocode,nspin,Lzd%Llr(ilr),orbs,&
!!!                   0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,pkernel,psi(ind:ind+dimwf-1),Lpot,eexctX)
!!!   
!!!           end if
!!!        end if
!!!     else
!!!        eexctX = 0._gp
!!!        !print *,'iproc,eexctX',iproc,eexctX
!!!     end if
!!!
!!!!     call timing(iproc,'ApplyLocPotKin','ON')
!!!
!!!     !apply the local hamiltonian for each of the orbitals
!!!     !given to each processor
!!!     !pot=0.d0
!!!     !psi=1.d0
!!!     !switch between GPU/CPU treatment
!!!!     do i=1,(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp
!!!!          call random_number(psi(i))
!!!!     end do
!!!
!!!     if(OCLconv .and. ASYNCconv) then
!!!       !allocate(hpsi2((Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor*&
!!!       !         Lzd%Llr(ilr)%Localnorb*nspin),stat=i_stat)
!!!       allocate(hpsi2((Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor*&
!!!                localnorb*nspin),stat=i_stat)
!!!       call memocc(i_stat,hpsi2,'hpsi2',subname)
!!!       hpsi(:)=0.0
!!!     else
!!!       hpsi2 => hpsi
!!!     end if
!!!     if (GPUconv) then  !does not work yet
!!!        call local_hamiltonian_GPU(iproc,orbs,Lzd%Llr(ilr),hx,hy,hz,nspin,Lpot,psi(ind:ind+dimwf-1),&
!!!             hpsi(ind:ind+dimwf-1),tmp_ekin_sum,tmp_epot_sum,GPU,ilr)
!!!     else if (OCLconv) then  ! does_not_work yet
!!!        call local_hamiltonian_OCL(iproc,orbs,Lzd%Llr(ilr),hx,hy,hz,nspin,Lpot,psi(ind:ind+dimwf-1),&
!!!             hpsi2,tmp_ekin_sum,tmp_epot_sum,GPU,ekin,epot,ilr)
!!!     else
!!!        call local_hamiltonian_LinearConfinement(iproc, nproc, ilr, orbs, lzd%Llr(ilr), localnorb, hx, hy, hz, &
!!!              nspin, size_Lpot, Lpot, psi(ind), hpsi(ind), tmp_ekin_sum, tmp_epot_sum, lin, at, rxyz, onWhichAtomp, withConfinement)
!!!        !!do i_stat=ind,ind+dimwf-1
!!!        !!    !write(7300+iproc,*) hpsi(i_stat)
!!!        !!    write(7310+iproc,*) hpsi(i_stat)
!!!        !!end do
!!!        !!do i_all=ind,ind+dimwf-1
!!!        !!    !!write(600+iproc,*) i_all, psi(i_all)
!!!        !!    !!write(610+iproc,*) i_all, hpsi(i_all)
!!!        !!    !!if((psi(i_all)==0.d0 .and. hpsi(i_all)/=0.d0) .or. (hpsi(i_all)==0.d0 .and. psi(i_all)/=0.d0))  then
!!!        !!    !!    write(*,*) 'ERRORi hamapp: iproc, i_all', iproc, i_all
!!!        !!    !!end if
!!!        !!end do
!!!     end if
!!!
!!!     ekin_sum = ekin_sum + tmp_ekin_sum
!!!     epot_sum = epot_sum + tmp_epot_sum
!!!
!!!  !test part to check the results wrt OCL convolutions
!!!!!$  if (OCLconv) then
!!!!!$     allocate(hpsi_OCL((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp+ndebug),stat=i_stat)
!!!!!$     call memocc(i_stat,hpsi_OCL,'hpsi_OCL',subname)
!!!!!$     print *,'fulllocam',GPU%full_locham
!!!!!$     call local_hamiltonian_OCL(iproc,orbs,at%geocode,lr,hx,hy,hz,nspin,pot,psi,hpsi,ek_fake,ep_fake,GPU)
!!!!!$     maxdiff=0.0_wp
!!!!!$     do i=1,(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp
!!!!!$        maxdiff=max(maxdiff,abs(hpsi(i)-hpsi_OCL(i)))
!!!!!$     end do
!!!!!$     print *,'maxdiff',maxdiff
!!!!!$     print *,'ekin_diff',abs(ek_fake-ekin_sum)
!!!!!$     print *,'epot_diff',abs(ep_fake-epot_sum)
!!!!!$     i_all=-product(shape(hpsi_OCL))*kind(hpsi_OCL)
!!!!!$     deallocate(hpsi_OCL,stat=i_stat)
!!!!!$     call memocc(i_stat,i_all,'hpsi_OCL',subname)
!!!!!$  end if
!!!
!!!!     call timing(iproc,'ApplyLocPotKin','OF')
!!!
!!!  !  apply all PSP projectors for all orbitals belonging to iproc
!!!!     call timing(iproc,'ApplyProj     ','ON')
!!!
!!!  !here the localisation region should be changed, temporary only for cubic approach
!!!  !   eproj_sum=0.0_gp
!!!
!!!  ! CUBIC STUFF
!!!  !apply the projectors following the strategy (On-the-fly calculation or not)
!!!!!!  if (DistProjApply .and. .not.present(Lzd)) then
!!!!!!     call applyprojectorsonthefly(iproc,orbs,at,lr,&
!!!!!!          rxyz,hx,hy,hz,lr%wfd,nlpspd,proj,psi,hpsi,eproj_sum)
!!!!!!  else if(orbs%norbp > 0 .and. .not.present(Lzd)) then
!!!!!!     !apply the projectors  k-point of the processor
!!!!!!     !starting k-point
!!!!!!     ikpt=orbs%iokpt(1)
!!!!!!     istart_ck=1
!!!!!!     ispsi_k=1
!!!!!!     loop_kpt: do
!!!!!!
!!!!!!        call orbs_in_kpt(ikpt,orbs,isorb,ieorb,nspinor)
!!!!!!
!!!!!!        ! loop over all my orbitals
!!!!!!        ispsi=ispsi_k
!!!!!!        do iorb=isorb,ieorb
!!!!!!           istart_c=istart_ck
!!!!!!           do iat=1,at%nat
!!!!!!              call apply_atproj_iorb(iat,iorb,istart_c,at,orbs,lr%wfd,nlpspd,&
!!!!!!                   proj,psi(ispsi),hpsi(ispsi),eproj_sum)
!!!!!!           end do
!!!!!!           ispsi=ispsi+(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*nspinor
!!!!!!        end do
!!!!!!        istart_ck=istart_c
!!!!!!        if (ieorb == orbs%norbp) exit loop_kpt
!!!!!!        ikpt=ikpt+1
!!!!!!        ispsi_k=ispsi
!!!!!!     end do loop_kpt
!!!!!!
!!!!!!     if (istart_ck-1 /= nlpspd%nprojel) stop 'incorrect once-and-for-all psp application'
!!!!!!     if (ispsi-1 /= (lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp) stop 'incorrect V_nl psi application'
!!!!!!
!!!!!!  END OF CUBIC STUFF
!!!
!!!     if(orbs%norbp > 0) then
!!!        !allocate
!!!        !if(ilr == 1) then
!!!        if(.not.allocated(hpsi_proj)) then
!!!           allocate(hpsi_proj(orbs%npsidim),stat=i_stat)
!!!           call memocc(i_stat,hpsi_proj,'hpsi_proj',subname)
!!!           hpsi_proj = 0.0_wp
!!!        end if
!!!
!!!        ! allocate projflg
!!!        allocate(Lzd%Llr(ilr)%projflg(at%nat),stat=i_stat)
!!!        call memocc(i_stat,Lzd%Llr(ilr)%projflg,'Lzd%Llr(ilr)%projflg',subname)
!!!
!!!        ! Make the local non-linear pseudopotentials descriptors
!!!        call nlpspd_to_locreg(input,iproc,Lzd%Glr,Lzd%Llr(ilr),rxyz,at,orbs,&
!!!      &      radii_cf,input%frmult,input%frmult,input%hx,input%hy,input%hz,Lzd%Gnlpspd,Lzd%Lnlpspd(ilr),Lzd%Llr(ilr)%projflg)
!!!
!!!        ! proj is declared with intent in, but in apply_local_projectors it has intent out. Therefore
!!!        ! copy it to projCopy. Maybe not needed, since it is not used anywhere else in this subroutine
!!!        !call apply_local_projectors(ilr,nspin,at,hx,hy,hz,Lzd%Llr(ilr),Lzd%Lnlpspd(ilr),proj,orbs,&
!!!        !         Lzd%Llr(ilr)%projflg,psi(ind:ind+dimwf-1),rxyz,hpsi(ind:ind+dimwf-1),eproj_sum)
!!!        !allocate(projCopy(Lzd%Lnlpspd(ilr)%nprojel), stat=i_stat)
!!!        !call dcopy(Lzd%Lnlpspd(ilr)%nprojel, proj, 1, projCopy, 1)
!!!        !call apply_local_projectors(ilr,nspin,at,hx,hy,hz,Lzd%Llr(ilr),Lzd%Lnlpspd(ilr),projCopy,orbs,&
!!!        !         Lzd%Llr(ilr)%projflg,psi(ind:ind+dimwf-1),rxyz,hpsi(ind:ind+dimwf-1),eproj_sum)
!!!                                                                    
!!!        call apply_local_projectors2(ilr,iproc,localnorb,nspin,at,hx,hy,hz,Lzd%Llr(ilr),Lzd%Lnlpspd(ilr),orbs,&
!!!                 Lzd%Llr(ilr)%projflg,psi(ind),rxyz,hpsi(ind),eproj_sum)
!!!        !!do i_stat=ind,ind+dimwf-1
!!!        !!    !write(7200+iproc,*) psi(i_stat)
!!!        !!    write(7210+iproc,*) psi(i_stat)
!!!        !!end do
!!!
!!!        !deallocate(projCopy, stat=i_stat)
!!!        ! accumulate the new hpsi
!!!        hpsi_proj(ind:ind+dimwf-1) = hpsi_proj(ind:ind+dimwf-1) + hpsi(ind:ind+dimwf-1)
!!!
!!!        ! Deallocate projflg
!!!        i_all=-product(shape(Lzd%Llr(ilr)%projflg))*kind(Lzd%Llr(ilr)%projflg)
!!!        deallocate(Lzd%Llr(ilr)%projflg,stat=i_stat)
!!!        call memocc(i_stat,i_all,'Lzd%Llr(ilr)%projflg',subname)
!!!     end if
!!!     ind = ind + dimwf
!!!
!!!     ! deallocate Lpot
!!!     call free_full_potential(nproc,Lpot,subname)
!!!
!!!  end do
!!!! END LINEAR MODIFICATIONS
!!!
!!!
!!!  ! Now that all is accumulated, rename hpsi_proj to hpsi
!!!  ! Do this only if hpsi_proj is allcoated, i.e. if it has some content. It is not
!!!  ! allocated if we cycled all the times.
!!!  if(allocated(hpsi_proj)) then
!!!      hpsi = hpsi_proj
!!!
!!!      !deallocate hpsi_proj
!!!      i_all=-product(shape(hpsi_proj))*kind(hpsi_proj)
!!!      deallocate(hpsi_proj,stat=i_stat)
!!!      call memocc(i_stat,i_all,'hpsi_proj',subname)
!!!  end if
!!!
!!!
!!!  !! local potential and kinetic energy for all orbitals belonging to iproc
!!!  !if (iproc==0 .and. verbose > 1) then
!!!  !   write(*,'(1x,a)',advance='no')&
!!!  !        'Hamiltonian application...'
!!!  !end if
!!!
!!!  if(OCLconv .and. ASYNCconv) then
!!!    call finish_hamiltonian_OCL(orbs,ekin_sum,epot_sum,GPU,ekin,epot)
!!!    call daxpy(size(hpsi), 1.0_wp, hpsi2(1), 1, hpsi(1),1)
!!!    i_all=-product(shape(hpsi2))*kind(hpsi2)
!!!    deallocate(hpsi2,stat=i_stat)
!!!    call memocc(i_stat,i_all,'hpsi2',subname)
!!!  endif
!!!
!!!!  call timing(iproc,'ApplyProj     ','OF')
!!!
!!!
!!!  if(energyReductionFlag) then
!!!      !energies reduction
!!!      if (nproc > 1) then
!!!         wrkallred(1,2)=ekin_sum
!!!         wrkallred(2,2)=epot_sum
!!!         wrkallred(3,2)=eproj_sum
!!!         call MPI_ALLREDUCE(wrkallred(1,2),wrkallred(1,1),3,&
!!!              mpidtypg,MPI_SUM,MPI_COMM_WORLD,ierr)
!!!         ekin_sum=wrkallred(1,1)
!!!         epot_sum=wrkallred(2,1)
!!!         eproj_sum=wrkallred(3,1)
!!!      endif
!!!  else
!!!      ! Do not sum up the energies, but set them to default values
!!!      ekin_sum=0.d0
!!!      epot_sum=0.d0
!!!      eproj_sum=0.d0
!!!  end if
!!!          
!!!
!!!
!!!
!!!  !up to this point, the value of the potential energy is 
!!!  !only taking into account the local potential part
!!!  !whereas it should consider also the value coming from the 
!!!  !exact exchange operator (twice the exact exchange energy)
!!!  if (exctX) epot_sum=epot_sum+2.0_gp*eexctX
!!!
!!!  do ilr=1,Lzd%nlr
!!!      call deallocate_nonlocal_psp_descriptors(Lzd%Lnlpspd(ilr), subname)
!!!  end do
!!!
!!!
!!!  !!do i_all=1,size(hpsi)
!!!  !!    !write(7000+iproc,*) hpsi(i_all)
!!!  !!    write(7010+iproc,*) hpsi(i_all)
!!!  !!end do
!!!
!!!END SUBROUTINE HamiltonianApplicationConfinement2



!!!>   Calculate the action of the local hamiltonian on the orbitals
!!subroutine local_hamiltonian_LinearConfinement(iproc, nproc, ilr, orbs, lr, norb, hx, hy, hz, &
!!     nspin, ndimpot, pot, psi, hpsi, ekin_sum, epot_sum, lin, at, rxyz, onWhichAtomp, withConfinement)
!!  use module_base
!!  use module_types
!!  use module_interfaces, exceptThisOne => local_hamiltonian_LinearConfinement
!!  use libxc_functionals
!!  implicit none
!!  integer, intent(in) :: iproc, nproc, nspin, ilr, norb, ndimpot
!!  real(gp), intent(in) :: hx, hy, hz
!!  type(orbitals_data), intent(in) :: orbs
!!  type(locreg_descriptors), intent(in) :: lr
!!  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*norb), intent(in) :: psi
!!  real(wp), dimension(ndimpot) :: pot
!!  real(gp), intent(out) :: ekin_sum,epot_sum
!!  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*norb), intent(out) :: hpsi
!!  type(linearParameters),intent(in):: lin
!!  type(atoms_data),intent(in):: at
!!  real(8),dimension(3,at%nat),intent(in):: rxyz
!!  integer,dimension(orbs%norbp),intent(in):: onWhichAtomp
!!  logical,intent(in):: withConfinement
!!  !local variables
!!  character(len=*), parameter :: subname='local_hamiltonian_Linear'
!!  integer :: i_all,i_stat,iorb,npot,nsoffset,oidx,ispot
!!  integer :: ii,orbtot
!!  !integer,dimension(lr%localnorb*nspin) :: inthisLocreg
!!  integer,dimension(norb*nspin) :: inthisLocreg
!!  real(wp) :: exctXcoeff
!!  real(gp) :: ekin,epot,kx,ky,kz,etest, hxh, hyh, hzh
!!  type(workarr_locham) :: wrk_lh
!!  real(wp), dimension(:,:), allocatable :: psir
!!integer:: i, j, jj
!!
!!  hxh=.5d0*hx
!!  hyh=.5d0*hy
!!  hzh=.5d0*hz
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
!!  call razero(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspinor,psir)
!!
!!  ekin_sum=0.0_gp
!!  epot_sum=0.0_gp
!!
!!  etest=0.0_gp
!!  
!!  orbtot = 0
!!  do iorb=1,orbs%norbp
!!     !if (orbs%inWhichLocreg(iorb) == ilr) then
!!     if (orbs%inWhichLocregp(iorb) == ilr) then
!!        orbtot = orbtot+1
!!        inthisLocreg(orbtot) = iorb
!!     end if
!!  end do 
!!  
!!  !if (orbtot .ne. lr%localnorb*nspin) then
!!  if (orbtot .ne. norb*nspin) then
!!     write(*,'(3(a,i0))') 'process ',iproc, ': Problem in local_hamiltonian_Linear, orbtot=',orbtot,&
!!     ' is not equal to localnorb=',lr%localnorb*nspin
!!     stop
!!  end if
!!  !write(*,'(4(a,i0))') 'iproc ',iproc,' handles ',orbtot,' orbitals in locreg ',ilr,'. Data per orbital=',lr%wfd%nvctr_c+7*lr%wfd%nvctr_f
!!
!!
!!  do ii=1,orbtot
!!
!!     iorb = inthisLocreg(ii)   !using ii and iorb to identify the orbitals because in linear case, the ordering is different
!!                               !orbitals are now orderer by locreg. So, iorb is the old numbering (i.e. in Global region)
!!                               !while ii is it's numbering in the locreg.
!!
!!     if(orbs%spinsgn(iorb+orbs%isorb)>0.0_gp .or. nspin == 1 .or. nspin == 4 ) then
!!        nsoffset=1
!!     else
!!        nsoffset=lr%d%n1i*lr%d%n2i*lr%d%n3i+1
!!     end if
!!
!!     oidx=(ii-1)*orbs%nspinor+1
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
!!        if(withConfinement) then
!!            call apply_potentialConfinement2(iproc, lr%d%n1,lr%d%n2,lr%d%n3,1,1,1,0,orbs%nspinor,npot,psir,&
!!                 pot(nsoffset),epot, rxyz(1,onWhichAtomp(iorb)), hxh, hyh, hzh, &
!!                 lin%potentialprefac(at%iatype(onWhichAtomp(iorb))), lin%confpotorder, &
!!                 lr%nsi1, lr%nsi2, lr%nsi3, &
!!                 lr%bounds%ibyyzz_r) !optional
!!        else
!!            !!call apply_potentialConfinement2(lr%d%n1,lr%d%n2,lr%d%n3,1,1,1,0,orbs%nspinor,npot,psir,&
!!            !!     pot(nsoffset),epot, rxyz(1,onWhichAtomp(iorb)), hxh, hyh, hzh, &
!!            !!     0.d0, lin%confpotorder, &
!!            !!     lr%nsi1, lr%nsi2, lr%nsi3, &
!!            !!     lr%bounds%ibyyzz_r) !optional
!!            call apply_potential(lr%d%n1,lr%d%n2,lr%d%n3,1,1,1,0,orbs%nspinor,npot,psir,&
!!                pot(nsoffset),epot,&
!!                lr%bounds%ibyyzz_r) !optional
!!        end if
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
!!        ispot=1+lr%d%n1i*lr%d%n2i*lr%d%n3i*(nspin+ii-1)
!!        !add to the psir function the part of the potential coming from the exact exchange
!!        call axpy(lr%d%n1i*lr%d%n2i*lr%d%n3i,exctXcoeff,pot(ispot),1,psir(1,1),1)
!!     end if
!!
!!     !apply the kinetic term, sum with the potential and transform back to Daubechies basis
!!     call isf_to_daub_kinetic(hx,hy,hz,kx,ky,kz,orbs%nspinor,lr,wrk_lh,&
!!          psir,hpsi(1,oidx),ekin)
!!!     print *,iorb, ekin+epot, epot
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
!!
!!END SUBROUTINE local_hamiltonian_LinearConfinement


!!$!> @file
!!$!!  Application of the Hamiltonian + orthonormalize constraints
!!$!! @author
!!$!!    Copyright (C) 2007-2011 CEA
!!$!!    This file is distributed under the terms of the
!!$!!    GNU General Public License, see ~/COPYING file
!!$!!    or http://www.gnu.org/copyleft/gpl.txt .
!!$!!    For the list of contributors, see ~/AUTHORS 
!!$
!!$
!!$!> Application of the Hamiltonian
!!$subroutine LinearHamiltonianApplication(input,iproc,nproc,at,Lzd,orbs,hx,hy,hz,rxyz,&
!!$     proj,ngatherarr,pot,psi,Lhpsi,&
!!$     ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU,radii_cf,pkernel,orbsocc,psirocc)
!!$  use module_base
!!$  use module_types
!!$  use libxc_functionals
!!$  use module_interfaces, exceptThisOne => LinearHamiltonianApplication
!!$  implicit none
!!$  integer, intent(in) :: iproc,nproc,nspin
!!$  real(gp), intent(in) :: hx,hy,hz
!!$  type(atoms_data), intent(in) :: at
!!$  type(input_variables), intent(in) :: input
!!$  type(local_zone_descriptors),intent(inout) :: Lzd
!!$  type(orbitals_data),intent(in) :: orbs
!!$  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
!!$  real(gp), dimension(3,at%nat), intent(in) :: rxyz
!!$  real(wp), dimension(Lzd%Gnlpspd%nprojel), intent(in) :: proj
!!$  real(wp), dimension(Lzd%Lpsidimtot), intent(in) :: psi
!!$  real(wp), dimension(:), pointer :: pot
!!$  real(gp), intent(out) :: ekin_sum,epot_sum,eexctX,eproj_sum
!!$  real(wp), target, dimension(Lzd%Lpsidimtot), intent(out) :: Lhpsi
!!$  type(GPU_pointers), intent(inout) :: GPU
!!$  real(gp), dimension(at%ntypes,3+ndebug), intent(in) :: radii_cf
!!$  real(dp), dimension(*), optional :: pkernel
!!$  type(orbitals_data), intent(in), optional :: orbsocc
!!$  real(wp), dimension(:), pointer, optional :: psirocc
!!$  !local variables
!!$  real(gp) :: tmp_ekin_sum,tmp_epot_sum,tmp_eproj_sum
!!$  real(gp), dimension(2,orbs%norbp) :: ekin
!!$  real(gp), dimension(2,orbs%norbp) :: epot
!!$  real(wp), dimension(:), pointer :: hpsi2
!!$  character(len=*), parameter :: subname='LinearHamiltonianApplication'
!!$  logical :: exctX,op2p
!!$  integer :: i_all,i_stat,ierr,iorb,n3p,ispot,istart_c,iat
!!$  integer :: istart_ck,isorb,ieorb,ikpt,ispsi_k,nspinor,ispsi
!!$  integer :: ii,ilr,dimwf,ind,size_Lpot,size_pot,size_potxc
!!$  integer :: tmp_norbp,iels
!!$  real(wp),dimension(:),pointer :: Lpot
!!$  real(wp),dimension(:),allocatable :: potxc
!!$  real(wp),dimension(:),allocatable :: hpsi_proj
!!$!OCL  integer, dimension(3) :: periodic
!!$!OCL  real(wp) :: maxdiff
!!$!OCL  real(gp) :: eproj,ek_fake,ep_fake
!!$  real(gp), dimension(3,2) :: wrkallred
!!$!OCL  real(wp), dimension(:), allocatable :: hpsi_OCL
!!$
!!$  !check if the potential has been associated
!!$  if (.not. associated(pot)) then
!!$     if (iproc ==0) then
!!$        write(*,*)' ERROR, LinearHamiltonianApplication, potential not associated!'
!!$
!!$        stop
!!$     end if
!!$  end if
!!$
!!$  !initialise exact exchange energy 
!!$  op2p=(eexctX == -99.0_gp)
!!$  eexctX=0.0_gp
!!$
!!$  exctX = libxc_functionals_exctXfac() /= 0.0_gp
!!$
!!$  ! Allocate the nonlocal descriptors for the locregs
!!$  !allocate(Lzd%Lnlpspd(Lzd%nlr),stat=i_stat)
!!$
!!$  ! NOW do Exact Exchange potential
!!$  ! Still using the cubic approach using the whole simulation box
!!$  ! and total orbitals
!!$    size_potxc = 0
!!$    if (exctX) then
!!$       size_potxc = max(max(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*orbs%norb,ngatherarr(0,1)*orbs%norb),1)
!!$       allocate(potxc(size_potxc+ndebug),stat=i_stat)
!!$       if (present(pkernel) .and. present(orbsocc) .and. present(psirocc)) then
!!$          call cubic_exact_exchange(iproc,nproc,nspin,Lzd%Lpsidimtot,size_potxc,hx,hy,hz,Lzd%Glr,orbs,&
!!$               ngatherarr,psi,potxc,eexctX,pkernel,orbsocc,psirocc)
!!$       else if(present(pkernel)) then
!!$           call cubic_exact_exchange(iproc,nproc,nspin,Lzd%Lpsidimtot,size_potxc,hx,hy,hz,Lzd%Glr,orbs,&
!!$               ngatherarr,psi,potxc,eexctX,pkernel=pkernel)
!!$       else 
!!$           call cubic_exact_exchange(iproc,nproc,nspin,Lzd%Lpsidimtot,size_potxc,hx,hy,hz,Lzd%Glr,orbs,&
!!$               ngatherarr,psi,potxc,eexctX)
!!$       end if
!!$    end if
!!$
!!$  !initialize accumulators
!!$  ekin_sum = 0.0_gp
!!$  epot_sum = 0.0_gp
!!$  eproj_sum= 0.0_gp
!!$  ind = 1
!!$  do ii= 1, orbs%norbp
!!$     ilr = orbs%inwhichlocreg(ii+orbs%isorb)
!!$
!!$     !determine the dimension of the potential array (copied from full_local_potential)
!!$     ! For now, using the whole set of orbitals in the Glr (could diminish to the Llr??)
!!$     if (exctX) then
!!$        size_pot=Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*nspin + size_potxc
!!$        size_Lpot=Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*nspin + size_potxc 
!!$     else
!!$        size_pot=Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*nspin
!!$        size_Lpot = Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*nspin
!!$     end if
!!$
!!$     allocate(Lpot(size_Lpot+ndebug), stat=i_stat)
!!$     call memocc(i_stat,Lpot,'Lpot',subname)
!!$     call razero(size_Lpot,Lpot)
!!$
!!$     ! Cut the potential into locreg pieces
!!$     call global_to_local(Lzd%Glr,Lzd%Llr(ilr),nspin,size_pot,size_Lpot,pot,Lpot)
!!$  
!!$     ! Set some quantities: ispot=shift for potential, dimwf=dimension of wavefunction
!!$     ispot=Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*nspin+1
!!$     dimwf = (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor
!!$
!!$     ! fill in the exact-exchange part
!!$     if (size_potxc > 0) then
!!$        do iels = 1, size_potxc
!!$           Lpot(ispot+iels-1) = potxc(iels)
!!$        end do
!!$     end if
!!$
!!$     call timing(iproc,'ApplyLocPotKin','ON')
!!$     if(OCLconv .and. ASYNCconv) then
!!$       allocate(hpsi2((Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor*&
!!$                orbs%norb),stat=i_stat)
!!$       call memocc(i_stat,hpsi2,'hpsi2',subname)
!!$       Lhpsi(:)=0.0
!!$     else
!!$       hpsi2 => Lhpsi
!!$     end if
!!$     if (GPUconv) then  !does not work yet
!!$        call local_hamiltonian_GPU(iproc,orbs,Lzd%Llr(ilr),hx,hy,hz,nspin,Lpot,psi(ind),&
!!$             Lhpsi(ind),tmp_ekin_sum,tmp_epot_sum,GPU,ilr)
!!$     else if (OCLconv) then  ! does_not_work yet
!!$        call local_hamiltonian_OCL(iproc,orbs,Lzd%Llr(ilr),hx,hy,hz,nspin,Lpot,psi(ind),&
!!$             hpsi2,tmp_ekin_sum,tmp_epot_sum,GPU,ekin,epot,ilr)
!!$     else
!!$        call local_hamiltonian_Linear(iproc,ii,orbs,Lzd%Llr(ilr),hx,hy,hz,nspin,Lpot,psi(ind),&
!!$             Lhpsi(ind),tmp_ekin_sum,tmp_epot_sum)
!!$     end if
!!$     call timing(iproc,'ApplyLocPotKin','OF')
!!$
!!$     ekin_sum = ekin_sum + tmp_ekin_sum
!!$     epot_sum = epot_sum + tmp_epot_sum
!!$
!!$     if(orbs%norbp > 0) then
!!$!        call timing(iproc,'create_nlpspd ','ON')
!!$        !allocate
!!$        allocate(hpsi_proj((Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor),stat=i_stat)
!!$        call memocc(i_stat,hpsi_proj,'hpsi_proj',subname)
!!$        hpsi_proj = 0.0_wp
!!$
!!$!        ! allocate projflg
!!$!        allocate(Lzd%Llr(ilr)%projflg(at%nat),stat=i_stat)
!!$!        call memocc(i_stat,Lzd%Llr(ilr)%projflg,'Lzd%Llr(ilr)%projflg',subname)
!!$!        
!!$!        ! Make the local non-linear pseudopotentials descriptors
!!$!        call nlpspd_to_locreg(input,iproc,Lzd%Glr,Lzd%Llr(ilr),rxyz,at,orbs,&
!!$!      &      radii_cf,input%frmult,input%frmult,hx,hy,hz,Lzd%Gnlpspd,Lzd%Lnlpspd(ilr),Lzd%Llr(ilr)%projflg)
!!$!        call timing(iproc,'create_nlpspd ','OF')
!!$
!!$        call timing(iproc,'ApplyProj     ','ON')
!!$        call apply_local_projectors(ii,iproc,nspin,at,hx,hy,hz,Lzd%Llr(ilr),Lzd%Lnlpspd(ilr),&
!!$                 orbs,Lzd%Llr(ilr)%projflg,psi(ind),rxyz,hpsi_proj(1),eproj_sum)
!!$        ! accumulate the new hpsi
!!$        Lhpsi(ind:ind+dimwf-1) = Lhpsi(ind:ind+dimwf-1) + hpsi_proj(1:dimwf)
!!$
!!$       !deallocate hpsi_proj
!!$       i_all=-product(shape(hpsi_proj))*kind(hpsi_proj)
!!$       deallocate(hpsi_proj,stat=i_stat)
!!$       call memocc(i_stat,i_all,'hpsi_proj',subname)
!!$       call timing(iproc,'ApplyProj     ','OF')
!!$
!!$     end if
!!$     ind = ind + dimwf
!!$
!!$     ! deallocate Lpot
!!$     deallocate(Lpot,stat=i_stat)
!!$     !call free_full_potential(nproc,Lpot,subname)
!!$  end do
!!$! END LINEAR MODIFICATIONS
!!$
!!$  ! local potential and kinetic energy for all orbitals belonging to iproc
!!$  if (iproc==0 .and. verbose > 1) then
!!$     write(*,'(1x,a)',advance='no')&
!!$          'Hamiltonian application...'
!!$  end if
!!$
!!$  if(OCLconv .and. ASYNCconv) then
!!$    call finish_hamiltonian_OCL(orbs,ekin_sum,epot_sum,GPU,ekin,epot)
!!$    call daxpy(size(Lhpsi), 1.0_wp, hpsi2(1), 1, Lhpsi(1),1)
!!$    i_all=-product(shape(hpsi2))*kind(hpsi2)
!!$    deallocate(hpsi2,stat=i_stat)
!!$    call memocc(i_stat,i_all,'hpsi2',subname)
!!$  endif
!!$
!!$!  call timing(iproc,'ApplyProj     ','OF')
!!$
!!$  !energies reduction
!!$  if (nproc > 1) then
!!$     wrkallred(1,2)=ekin_sum
!!$     wrkallred(2,2)=epot_sum
!!$     wrkallred(3,2)=eproj_sum
!!$     call MPI_ALLREDUCE(wrkallred(1,2),wrkallred(1,1),3,&
!!$          mpidtypg,MPI_SUM,MPI_COMM_WORLD,ierr)
!!$     ekin_sum=wrkallred(1,1)
!!$     epot_sum=wrkallred(2,1)
!!$     eproj_sum=wrkallred(3,1)
!!$  endif
!!$
!!$  !up to this point, the value of the potential energy is 
!!$  !only taking into account the local potential part
!!$  !whereas it should consider also the value coming from the 
!!$  !exact exchange operator (twice the exact exchange energy)
!!$
!!$  if (exctX) epot_sum=epot_sum+2.0_gp*eexctX
!!$
!!$END SUBROUTINE LinearHamiltonianApplication


