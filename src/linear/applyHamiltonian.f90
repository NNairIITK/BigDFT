subroutine apply_potentialConfinement2(n1,n2,n3,nl1,nl2,nl3,nbuf,nspinor,npot,psir,pot,epot, &
     rxyzConfinement, hxh, hyh, hzh, potentialPrefac, confPotOrder, offsetx, offsety, offsetz, &
     ibyyzz_r) !optional
!
! Purpose:
! ========
!   Routine for applying the local potentials. Supports the non-collinear case, the buffer for tails 
!   and different Boundary Conditions. Optimal also for the complex wavefuntion case.
!   The potential includes also the confinement potential.
!
! Calling arguments:
! ==================
!   Input arguments:
!   ----------------
!     n1
!     n2
!     n3
!     nl1
!     nl2
!     nl3
!     nbuf                       ???
!     nspinor                    real or complex?
!     npot
!     pot                        the potential to be applied. The confinement potential
!                                  will be applied on the fly.
!     ibyyzz_r (optional)
!     rxyzConfinement            the center for the confinement potential
!     hxh                        the grid spacing in x direction divided by 2
!     hyh                        the grid spacing in y direction divided by 2
!     hzh                        the grid spacing in z direction divided by 2
!     potentialPrefac            the prefactor for the confinement potential
!   Input / Output arguments:
!   -------------------------
!     psir                       wave function om real space grid
!   Output arguments:
!   -----------------
!     epot                       potential energy
!
use module_base
implicit none
integer, intent(in) :: n1,n2,n3,nl1,nl2,nl3,nbuf,nspinor,npot, confPotOrder, offsetx, offsety, offsetz
real(wp), dimension(-14*nl1:2*n1+1+15*nl1,-14*nl2:2*n2+1+15*nl2,-14*nl3:2*n3+1+15*nl3,nspinor), intent(inout) :: psir
real(wp), dimension(-14*nl1:2*n1+1+15*nl1-4*nbuf,-14*nl2:2*n2+1+15*nl2-4*nbuf,&
     -14*nl3:2*n3+1+15*nl3-4*nbuf,npot), intent(in) :: pot
integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r
real(gp), intent(out) :: epot
real(8),dimension(3),intent(in):: rxyzConfinement
real(8),intent(in):: hxh, hyh, hzh, potentialPrefac
!local variables
integer :: i1,i2,i3,i1s,i1e,ispinor, order
real(wp) :: tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42,tt
real(wp) :: psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4
real(gp) :: epot_p


  !the Tail treatment is allowed only in the Free BC case
  if (nbuf /= 0 .and. nl1*nl2*nl3 == 0) stop 'NONSENSE: nbuf/=0 only for Free BC'

  ! The order of the cofinement potential (we take order divided by two, 
  ! since later we calculate (r**2)**order.
  if(confPotOrder==2) then
      ! parabolic potential
      order=1
  else if(confPotOrder==4) then
      ! quartic potential
      order=2
  else if(confPotOrder==6) then
      ! sextic potential
      order=3
  end if
  
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
                       ! Add the quartic confinement potential to the potential.
                       !tt=(hxh*dble(i1)-rxyzConfinement(1))**2 + (hyh*dble(i2)-rxyzConfinement(2))**2 + &
                       !    (hzh*dble(i3)-rxyzConfinement(3))**2
                       tt=(hxh*dble(i1+offsetx)-rxyzConfinement(1))**2 + (hyh*dble(i2+offsety)-rxyzConfinement(2))**2 + &
                           (hzh*dble(i3+offsetz)-rxyzConfinement(3))**2
                       tt=potentialPrefac*tt**order
                       tt=pot(i1-2*nbuf,i2-2*nbuf,i3-2*nbuf,1)+tt
                       tt=tt*psir(i1,i2,i3,ispinor)
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

END SUBROUTINE apply_potentialConfinement2
!!***




!> Application of the Hamiltonian
subroutine HamiltonianApplicationConfinement2(input,iproc,nproc,at,Lzd,lin,hx,hy,hz,rxyz,&
     proj,ngatherarr,ndimpot,pot,psi,hpsi,&
     ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU,radii_cf, comgp, onWhichAtomp, withConfinement, &
     pkernel,orbsocc,psirocc)
  use module_base
  use module_types
  use libxc_functionals
  use module_interfaces, exceptThisOne => HamiltonianApplicationConfinement2
  implicit none
  integer, intent(in) :: iproc,nproc,nspin,ndimpot
  real(gp), intent(in) :: hx,hy,hz
  type(atoms_data), intent(in) :: at
  type(input_variables), intent(in) :: input
  type(linear_zone_descriptors),intent(inout) :: Lzd
  type(linearParameters),intent(in):: lin
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(wp), dimension(Lzd%Gnlpspd%nprojel), intent(in) :: proj
  real(wp), dimension(lin%Lorbs%npsidim), intent(in) :: psi
  real(wp), dimension(max(ndimpot,1)*nspin), intent(in) :: pot
  !real(wp), dimension(:), pointer :: pot
  real(gp), intent(out) :: ekin_sum,epot_sum,eexctX,eproj_sum
  real(wp), target, dimension(lin%Lorbs%npsidim), intent(out) :: hpsi
  type(GPU_pointers), intent(inout) :: GPU
  real(gp), dimension(at%ntypes,3+ndebug), intent(in) :: radii_cf
  type(p2pCommsgatherPot), intent(in):: comgp
  integer,dimension(lzd%orbs%norbp),intent(in):: onWhichAtomp
  logical,intent(in):: withConfinement
  real(dp), dimension(*), optional :: pkernel
  type(orbitals_data), intent(in), optional :: orbsocc
  real(wp), dimension(:), pointer, optional :: psirocc
  !local variables
  real(gp) :: tmp_ekin_sum,tmp_epot_sum,tmp_eproj_sum
  real(gp), dimension(2,Lzd%orbs%norbp) :: ekin
  real(gp), dimension(2,Lzd%orbs%norbp) :: epot
  real(wp), dimension(:), pointer :: hpsi2, temparr
  character(len=*), parameter :: subname='LinearHamiltonianApplicationConfinement2'
  logical :: exctX,op2p
  integer :: i_all,i_stat,ierr,iorb,n3p,ispot,istart_c,iat, i3s, i3e, ind1, ind2, ldim, gdim, jlr
  integer :: istart_ck,isorb,ieorb,ikpt,ispsi_k,nspinor,ispsi
  integer :: ilr,dimwf,ind,size_Lpot,size_pot
  integer :: tmp_norbp, istorb
  real(dp),dimension(:),pointer:: Lpot
  real(wp),dimension(:),allocatable :: hpsi_proj
  real(wp),dimension(:),allocatable:: projCopy
!OCL  integer, dimension(3) :: periodic
!OCL  real(wp) :: maxdiff
!OCL  real(gp) :: eproj,ek_fake,ep_fake
  real(gp), dimension(3,2) :: wrkallred
!OCL  real(wp), dimension(:), allocatable :: hpsi_OCL



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
  istorb=1
  do ilr= 1, Lzd%nlr
     ! Cycle if the process does not have any orbitals belonging
     ! to this localization region.
     !write(*,'(a,3i8)') 'iproc, ilr, Lzd%Llr(ilr)%Localnorb', iproc, ilr, Lzd%Llr(ilr)%Localnorb
     if(Lzd%Llr(ilr)%Localnorb == 0) cycle
 
     allocate(Lpot(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*nspin), stat=i_stat)
     call memocc(i_stat,Lpot,'Lpot',subname)
 
 
     !determine the dimension of the potential array (copied from full_local_potential)
     if (exctX) then
        stop 'exctX not yet implemented!'
        size_pot=Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*nspin + &
         max(max(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*Lzd%Llr(ilr)%Localnorb*nspin,ngatherarr(0,1)*Lzd%orbs%norb),1) !part which refers to exact exchange
        size_Lpot=Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*nspin + &
           max(max(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*Lzd%Llr(ilr)%Localnorb*nspin,&
           ngatherarr(0,1)*Lzd%orbs%norb),1) !CHECK THIS...DOES NOT WORK YET
     else
        !size_pot=Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*nspin
        !size_pot=Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*nscatterarr(iproc,2)*nspin
        size_Lpot = Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*nspin
     end if
 
     ! Extract the part of the potential which is needed for the current localization region.
     i3s=lzd%Llr(ilr)%nsi3-comgp%ise3(1,iproc)+2 ! starting index of localized potential with respect to total potential in comgp%recvBuf
     i3e=lzd%Llr(ilr)%nsi3+lzd%Llr(ilr)%d%n3i-comgp%ise3(1,iproc)+1 ! ending index of localized potential with respect to total potential in comgp%recvBuf
     if(i3e-i3s+1 /= Lzd%Llr(ilr)%d%n3i) then
         write(*,'(a,i0,3x,i0)') 'ERROR: i3e-i3s+1 /= Lzd%Llr(ilr)%d%n3i', i3e-i3s+1, Lzd%Llr(ilr)%d%n3i
         stop
     end if
     call global_to_local_parallel(lzd%Glr, lzd%Llr(ilr), nspin, ndimpot, size_Lpot, pot, Lpot, i3s, i3e)

     ! Set some quantities: ispot=shift for potential, dimwf=dimension of wavefunction
     ispot=Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*nspin+1
     dimwf=(Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*Lzd%Llr(ilr)%Localnorb*&
           Lzd%orbs%nspinor*nspin

     ! EXACT EXCHANGE NOT TESTED: SHOULD CHECK IF EVERYTHING IF FINE
     !fill the rest of the potential with the exact-exchange terms
     if (present(pkernel) .and. exctX) then
        stop 'exctx not yet implemented!'
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
        call local_hamiltonian_LinearConfinement(iproc, nproc, ilr, lzd%orbs, lzd%Llr(ilr), lzd%Llr(ilr)%localnorb, hx, hy, hz, &
              nspin, size_Lpot, Lpot, psi(ind), hpsi(ind), tmp_ekin_sum, tmp_epot_sum, lin, at, rxyz, onWhichAtomp, withConfinement)
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
        !if(ilr == 1) then
        if(.not.allocated(hpsi_proj)) then
           allocate(hpsi_proj(Lzd%orbs%npsidim),stat=i_stat)
           call memocc(i_stat,hpsi_proj,'hpsi_proj',subname)
           hpsi_proj = 0.0_wp
        end if

        ! allocate projflg
        allocate(Lzd%Llr(ilr)%projflg(at%nat),stat=i_stat)
        call memocc(i_stat,Lzd%Llr(ilr)%projflg,'Lzd%Llr(ilr)%projflg',subname)

        ! Make the local non-linear pseudopotentials descriptors
        call nlpspd_to_locreg(input,iproc,Lzd%Glr,Lzd%Llr(ilr),rxyz,at,Lzd%orbs,&
      &      radii_cf,input%frmult,input%frmult,input%hx,input%hy,input%hz,Lzd%Gnlpspd,Lzd%Lnlpspd(ilr),Lzd%Llr(ilr)%projflg)

        ! proj is declared with intent in, but in apply_local_projectors it has intent out. Therefore
        ! copy it to projCopy. Maybe not needed, since it is not used anywhere else in this subroutine
        !call apply_local_projectors(ilr,nspin,at,hx,hy,hz,Lzd%Llr(ilr),Lzd%Lnlpspd(ilr),proj,Lzd%orbs,&
        !         Lzd%Llr(ilr)%projflg,psi(ind:ind+dimwf-1),rxyz,hpsi(ind:ind+dimwf-1),eproj_sum)
        !allocate(projCopy(Lzd%Lnlpspd(ilr)%nprojel), stat=i_stat)
        !call dcopy(Lzd%Lnlpspd(ilr)%nprojel, proj, 1, projCopy, 1)
        !call apply_local_projectors(ilr,nspin,at,hx,hy,hz,Lzd%Llr(ilr),Lzd%Lnlpspd(ilr),projCopy,Lzd%orbs,&
        !         Lzd%Llr(ilr)%projflg,psi(ind:ind+dimwf-1),rxyz,hpsi(ind:ind+dimwf-1),eproj_sum)
                                   
        call apply_local_projectors(ilr,nspin,at,hx,hy,hz,Lzd%Llr(ilr),Lzd%Lnlpspd(ilr),Lzd%orbs,Lzd%orbs,&
                 Lzd%Llr(ilr)%projflg,psi(ind:ind+dimwf-1),rxyz,hpsi(ind:ind+dimwf-1),eproj_sum)
        !deallocate(projCopy, stat=i_stat)
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

  !! local potential and kinetic energy for all orbitals belonging to iproc
  !if (iproc==0 .and. verbose > 1) then
  !   write(*,'(1x,a)',advance='no')&
  !        'Hamiltonian application...'
  !end if

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

END SUBROUTINE HamiltonianApplicationConfinement2



!>   Calculate the action of the local hamiltonian on the orbitals
subroutine local_hamiltonian_LinearConfinement(iproc, nproc, ilr, orbs, lr, norb, hx, hy, hz, &
     nspin, ndimpot, pot, psi, hpsi, ekin_sum, epot_sum, lin, at, rxyz, onWhichAtomp, withConfinement)
  use module_base
  use module_types
  use module_interfaces, exceptThisOne => local_hamiltonian_LinearConfinement
  use libxc_functionals
  implicit none
  integer, intent(in) :: iproc, nproc, nspin, ilr, norb, ndimpot
  real(gp), intent(in) :: hx, hy, hz
  type(orbitals_data), intent(in) :: orbs
  type(locreg_descriptors), intent(in) :: lr
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*norb), intent(in) :: psi
  real(wp), dimension(ndimpot) :: pot
  real(gp), intent(out) :: ekin_sum,epot_sum
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*norb), intent(out) :: hpsi
  type(linearParameters),intent(in):: lin
  type(atoms_data),intent(in):: at
  real(8),dimension(3,at%nat),intent(in):: rxyz
  integer,dimension(orbs%norbp),intent(in):: onWhichAtomp
  logical,intent(in):: withConfinement
  !local variables
  character(len=*), parameter :: subname='local_hamiltonian_Linear'
  integer :: i_all,i_stat,iorb,npot,nsoffset,oidx,ispot
  integer :: ii,orbtot
  integer,dimension(lr%localnorb*nspin) :: inthisLocreg
  real(wp) :: exctXcoeff
  real(gp) :: ekin,epot,kx,ky,kz,etest, hxh, hyh, hzh
  type(workarr_locham) :: wrk_lh
  real(wp), dimension(:,:), allocatable :: psir

  hxh=.5d0*hx
  hyh=.5d0*hy
  hzh=.5d0*hz

  exctXcoeff=libxc_functionals_exctXfac()

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
  
  if (orbtot .ne. lr%localnorb*nspin) then
     write(*,'(3(a,i0))') 'process ',iproc, ': Problem in local_hamiltonian_Linear, orbtot=',orbtot,&
     ' is not equal to localnorb=',lr%localnorb*nspin
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

        if(withConfinement) then
            call apply_potentialConfinement2(lr%d%n1,lr%d%n2,lr%d%n3,1,1,1,0,orbs%nspinor,npot,psir,&
                 pot(nsoffset),epot, rxyz(1,onWhichAtomp(iorb)), hxh, hyh, hzh, &
                 lin%potentialprefac(at%iatype(onWhichAtomp(iorb))), lin%confpotorder, &
                 lr%nsi1, lr%nsi2, lr%nsi3, &
                 lr%bounds%ibyyzz_r) !optional
        else
            !!call apply_potentialConfinement2(lr%d%n1,lr%d%n2,lr%d%n3,1,1,1,0,orbs%nspinor,npot,psir,&
            !!     pot(nsoffset),epot, rxyz(1,onWhichAtomp(iorb)), hxh, hyh, hzh, &
            !!     0.d0, lin%confpotorder, &
            !!     lr%nsi1, lr%nsi2, lr%nsi3, &
            !!     lr%bounds%ibyyzz_r) !optional
            call apply_potential(lr%d%n1,lr%d%n2,lr%d%n3,1,1,1,0,orbs%nspinor,npot,psir,&
                pot(nsoffset),epot,&
                lr%bounds%ibyyzz_r) !optional
        end if

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
        ispot=1+lr%d%n1i*lr%d%n2i*lr%d%n3i*(nspin+ii-1)
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

  call deallocate_work_arrays_locham(lr,wrk_lh)

END SUBROUTINE local_hamiltonian_LinearConfinement

