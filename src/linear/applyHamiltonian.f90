!!!!subroutine apply_potentialConfinement2(iproc, n1,n2,n3,nl1,nl2,nl3,nbuf,nspinor,npot,psir,pot,epot, &
!!!!     rxyzConfinement, hxh, hyh, hzh, potentialPrefac, confPotOrder, offsetx, offsety, offsetz, &
!!!!     ibyyzz_r) !optional
!!!!!
!!!!! Purpose:
!!!!! ========
!!!!!   Routine for applying the local potentials. Supports the non-collinear case, the buffer for tails 
!!!!!   and different Boundary Conditions. Optimal also for the complex wavefuntion case.
!!!!!   The potential includes also the confinement potential.
!!!!!
!!!!! Calling arguments:
!!!!! ==================
!!!!!   Input arguments:
!!!!!   ----------------
!!!!!     n1
!!!!!     n2
!!!!!     n3
!!!!!     nl1
!!!!!     nl2
!!!!!     nl3
!!!!!     nbuf                       ???
!!!!!     nspinor                    real or complex?
!!!!!     npot
!!!!!     pot                        the potential to be applied. The confinement potential
!!!!!                                  will be applied on the fly.
!!!!!     ibyyzz_r (optional)
!!!!!     rxyzConfinement            the center for the confinement potential
!!!!!     hxh                        the grid spacing in x direction divided by 2
!!!!!     hyh                        the grid spacing in y direction divided by 2
!!!!!     hzh                        the grid spacing in z direction divided by 2
!!!!!     potentialPrefac            the prefactor for the confinement potential
!!!!!   Input / Output arguments:
!!!!!   -------------------------
!!!!!     psir                       wave function om real space grid
!!!!!   Output arguments:
!!!!!   -----------------
!!!!!     epot                       potential energy
!!!!!
!!!!use module_base
!!!!implicit none
!!!!integer, intent(in) :: iproc, n1,n2,n3,nl1,nl2,nl3,nbuf,nspinor,npot, confPotOrder, offsetx, offsety, offsetz
!!!!real(wp), dimension(-14*nl1:2*n1+1+15*nl1,-14*nl2:2*n2+1+15*nl2,-14*nl3:2*n3+1+15*nl3,nspinor), intent(inout) :: psir
!!!!real(wp), dimension(-14*nl1:2*n1+1+15*nl1-4*nbuf,-14*nl2:2*n2+1+15*nl2-4*nbuf,&
!!!!     -14*nl3:2*n3+1+15*nl3-4*nbuf,npot), intent(in) :: pot
!!!!integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r
!!!!real(gp), intent(out) :: epot
!!!!real(8),dimension(3),intent(in):: rxyzConfinement
!!!!real(8),intent(in):: hxh, hyh, hzh, potentialPrefac
!!!!!local variables
!!!!integer :: i1,i2,i3,i1s,i1e,ispinor, order
!!!!real(wp) :: tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42,tt
!!!!real(wp) :: psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4
!!!!real(gp) :: epot_p
!!!!
!!!!
!!!!
!!!!
!!!!  !the Tail treatment is allowed only in the Free BC case
!!!!  if (nbuf /= 0 .and. nl1*nl2*nl3 == 0) stop 'NONSENSE: nbuf/=0 only for Free BC'
!!!!
!!!!  ! The order of the cofinement potential (we take order divided by two, 
!!!!  ! since later we calculate (r**2)**order.
!!!!  if(confPotOrder==2) then
!!!!      ! parabolic potential
!!!!      order=1
!!!!  else if(confPotOrder==4) then
!!!!      ! quartic potential
!!!!      order=2
!!!!  else if(confPotOrder==6) then
!!!!      ! sextic potential
!!!!      order=3
!!!!  end if
!!!!  
!!!!  epot=0.0_wp
!!!!
!!!!!!!$omp parallel default(private)&
!!!!!!!$omp shared(pot,psir,n1,n2,n3,epot,ibyyzz_r,nl1,nl2,nl3,nbuf,nspinor)
!!!!  !case without bounds
!!!!  i1s=-14*nl1
!!!!  i1e=2*n1+1+15*nl1
!!!!  epot_p=0._gp
!!!!!!!$omp do
!!!!!write(*,*) 'iproc, -14*nl3,2*n3+1+15*nl3', iproc, -14*nl3,2*n3+1+15*nl3
!!!!!write(*,*) 'iproc, -14*nl2,2*n2+1+15*nl2', iproc, -14*nl2,2*n2+1+15*nl2
!!!!!write(*,'(a,i5,3es14.6)') 'iproc, confinement center (on grid): ', iproc, rxyzConfinement(1)/hxh, rxyzConfinement(2)/hyh, rxyzConfinement(3)/hzh
!!!!  do i3=-14*nl3,2*n3+1+15*nl3
!!!!     if (i3 >= -14+2*nbuf .and. i3 <= 2*n3+16-2*nbuf) then !check for the nbuf case
!!!!        do i2=-14*nl2,2*n2+1+15*nl2
!!!!           if (i2 >= -14+2*nbuf .and. i2 <= 2*n2+16-2*nbuf) then !check for the nbuf case
!!!!              !this if statement is inserted here for avoiding code duplication
!!!!              !it is to be seen whether the code results to be too much unoptimised
!!!!              if (present(ibyyzz_r)) then
!!!!                 !in this case we are surely in Free BC
!!!!                 !the min is to avoid to calculate for no bounds
!!!!                 do i1=-14+2*nbuf,min(ibyyzz_r(1,i2,i3),ibyyzz_r(2,i2,i3))-14-1
!!!!                    psir(i1,i2,i3,:)=0.0_wp
!!!!                 enddo
!!!!                 i1s=max(ibyyzz_r(1,i2,i3)-14,-14+2*nbuf)
!!!!                 i1e=min(ibyyzz_r(2,i2,i3)-14,2*n1+16-2*nbuf)
!!!!              end if
!!!!              !write(*,'(a,5i8)') 'iproc, i1, i2, i1s, i1e', iproc, i1, i2, i1s, i1e
!!!!              
!!!!              !here we put the branchments wrt to the spin
!!!!              if (nspinor == 4) then
!!!!                 do i1=i1s,i1e
!!!!                    !wavefunctions
!!!!                    psir1=psir(i1,i2,i3,1)
!!!!                    psir2=psir(i1,i2,i3,2)
!!!!                    psir3=psir(i1,i2,i3,3)
!!!!                    psir4=psir(i1,i2,i3,4)
!!!!                    !potentials
!!!!                    pot1=pot(i1-2*nbuf,i2-2*nbuf,i3-2*nbuf,1)
!!!!                    pot2=pot(i1-2*nbuf,i2-2*nbuf,i3-2*nbuf,2)
!!!!                    pot3=pot(i1-2*nbuf,i2-2*nbuf,i3-2*nbuf,3)
!!!!                    pot4=pot(i1-2*nbuf,i2-2*nbuf,i3-2*nbuf,4)
!!!!
!!!!                    !diagonal terms
!!!!                    tt11=pot1*psir1 !p1
!!!!                    tt22=pot1*psir2 !p2
!!!!                    tt33=pot4*psir3 !p3
!!!!                    tt44=pot4*psir4 !p4
!!!!                    !Rab*Rb
!!!!                    tt13=pot2*psir3 !p1
!!!!                    !Iab*Ib
!!!!                    tt14=pot3*psir4 !p1
!!!!                    !Rab*Ib
!!!!                    tt23=pot2*psir4 !p2
!!!!                    !Iab*Rb
!!!!                    tt24=pot3*psir3 !p2
!!!!                    !Rab*Ra
!!!!                    tt31=pot2*psir1 !p3
!!!!                    !Iab*Ia
!!!!                    tt32=pot3*psir2 !p3
!!!!                    !Rab*Ia
!!!!                    tt41=pot2*psir2 !p4
!!!!                    !Iab*Ra
!!!!                    tt42=pot3*psir1 !p4
!!!!
!!!!                    ! Change epot later
!!!!                    epot_p=epot_p+tt11*psir1+tt22*psir2+tt33*psir3+tt44*psir4+&
!!!!                         2.0_gp*tt31*psir3-2.0_gp*tt42*psir4+2.0_gp*tt41*psir4+2.0_gp*tt32*psir3
!!!!
!!!!                    !wavefunction update
!!!!                    !p1=h1p1+h2p3-h3p4
!!!!                    !p2=h1p2+h2p4+h3p3
!!!!                    !p3=h2p1+h3p2+h4p3
!!!!                    !p4=h2p2-h3p1+h4p4
!!!!                    psir(i1,i2,i3,1)=tt11+tt13-tt14
!!!!                    psir(i1,i2,i3,2)=tt22+tt23+tt24
!!!!                    psir(i1,i2,i3,3)=tt33+tt31+tt32
!!!!                    psir(i1,i2,i3,4)=tt44+tt41-tt42
!!!!                 end do
!!!!              else
!!!!                 do ispinor=1,nspinor
!!!!                    do i1=i1s,i1e
!!!!                       !the local potential is always real
!!!!                       ! Add the quartic confinement potential to the potential.
!!!!                        tt=(hxh*dble(i1+offsetx)-rxyzConfinement(1))**2 + (hyh*dble(i2+offsety)-rxyzConfinement(2))**2 + &
!!!!                            (hzh*dble(i3+offsetz)-rxyzConfinement(3))**2
!!!!                        !!!!! EXPERIMENTAL ########################
!!!!                        !!tt=sqrt(tt)
!!!!                        !!tt=max(tt-3.d0,0.d0)
!!!!                        !!tt=tt**2
!!!!                        !!!########################################
!!!!
!!!!
!!!!                       !!! New trial
!!!!                       !!tt=(hxh*dble(i1+offsetx)-rxyzConfinement(1))**2 + (hyh*dble(i2+offsety)-rxyzConfinement(2))**2 + &
!!!!                       !!   (hzh*dble(i3+offsetz)-rxyzConfinement(3))**2
!!!!                       !!tt=sqrt(tt)
!!!!                       !!tt=tt/5.d0
!!!!                       !!tt=tt**2
!!!!                       !!tt=5.d-2*(exp(tt)-1.d0)
!!!!
!!!!
!!!!                       tt=potentialPrefac*tt**order
!!!!
!!!!
!!!!                       !!tt=(hxh*dble(i1+offsetx)-rxyzConfinement(1))**2 + (hyh*dble(i2+offsety)-rxyzConfinement(2))**2 + &
!!!!                       !!    (hzh*dble(i3+offsetz)-rxyzConfinement(3))**2
!!!!                       !!tt=.5d0*potentialPrefac*tt**2+.5d0*.01d0*potentialPrefac*tt**3
!!!!                       tt=pot(i1-2*nbuf,i2-2*nbuf,i3-2*nbuf,1)+tt
!!!!                       tt=tt*psir(i1,i2,i3,ispinor)
!!!!                       epot_p=epot_p+real(tt*psir(i1,i2,i3,ispinor),gp)
!!!!                       psir(i1,i2,i3,ispinor)=tt
!!!!                    end do
!!!!                 end do
!!!!              end if
!!!!              
!!!!              if (present(ibyyzz_r)) then
!!!!                 !the max is to avoid the calculation for no bounds
!!!!                 do i1=max(ibyyzz_r(1,i2,i3),ibyyzz_r(2,i2,i3))-14+1,2*n1+16-2*nbuf
!!!!                    psir(i1,i2,i3,:)=0.0_wp
!!!!                 enddo
!!!!              end if
!!!!
!!!!           else
!!!!              do i1=-14,2*n1+16
!!!!                 psir(i1,i2,i3,:)=0.0_wp
!!!!              enddo
!!!!           endif
!!!!        enddo
!!!!     else
!!!!        do i2=-14,2*n2+16
!!!!           do i1=-14,2*n1+16
!!!!              psir(i1,i2,i3,:)=0.0_wp
!!!!           enddo
!!!!        enddo
!!!!     endif
!!!!  enddo
!!!!!!!$omp end do
!!!!
!!!!
!!!!!!!$omp critical
!!!!  epot=epot+epot_p
!!!!!!!$omp end critical
!!!!
!!!!!!!$omp end parallel
!!!!
!!!!END SUBROUTINE apply_potentialConfinement2
!!!!!!***


!> @file
!!  Routine to calculate the action of the hamiltonian
!! @author
!!   Copyright (C) 2005-2011 BigDFT group 
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 


!>   Calculate the action of the local hamiltonian on the orbitals
subroutine local_hamiltonian_Linear(iproc,iorb,orbs,lr,hx,hy,hz,&
     nspin,pot,psi,hpsi,ekin_sum,epot_sum)
  use module_base
  use module_types
  use module_interfaces
  use libxc_functionals
  implicit none
  integer, intent(in) :: iproc,nspin,iorb
  real(gp), intent(in) :: hx,hy,hz
  type(orbitals_data), intent(in) :: orbs
  type(locreg_descriptors), intent(in) :: lr
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor), intent(in) :: psi
  real(wp), dimension(*) :: pot
  !real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i*nspin) :: pot
  real(gp), intent(out) :: ekin_sum,epot_sum
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor), intent(inout) :: hpsi
  !local variables
  character(len=*), parameter :: subname='local_hamiltonian_Linear'
  integer :: i_all,i_stat,npot,nsoffset,oidx,ispot
  real(wp) :: exctXcoeff
  real(gp) :: ekin,epot,kx,ky,kz,etest
  type(workarr_locham) :: wrk_lh
  real(wp), dimension(:,:), allocatable :: psir

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

  if(orbs%spinsgn(iorb+orbs%isorb)>0.0_gp .or. nspin == 1 .or. nspin == 4 ) then
     nsoffset=1
  else
     nsoffset=lr%d%n1i*lr%d%n2i*lr%d%n3i+1
  end if

!  oidx=(iorb-1)*orbs%nspinor+1
  oidx = 1
  !transform the wavefunction in Daubechies basis to the wavefunction in ISF basis
  !the psir wavefunction is given in the spinorial form
  call daub_to_isf_locham(orbs%nspinor,lr,wrk_lh,psi(1,1),psir)

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
     ispot=1+lr%d%n1i*lr%d%n2i*lr%d%n3i*(orbs%nspin+iorb-1)
     !add to the psir function the part of the potential coming from the exact exchange
     call axpy(lr%d%n1i*lr%d%n2i*lr%d%n3i,exctXcoeff,pot(ispot),1,psir(1,oidx),1)
  end if

  !apply the kinetic term, sum with the potential and transform back to Daubechies basis
  call isf_to_daub_kinetic(hx,hy,hz,kx,ky,kz,orbs%nspinor,lr,wrk_lh,&
       psir,hpsi(1,oidx),ekin)

  ekin_sum=ekin_sum+orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*ekin
  epot_sum=epot_sum+orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*epot

  !deallocations of work arrays
  i_all=-product(shape(psir))*kind(psir)
  deallocate(psir,stat=i_stat)
  call memocc(i_stat,i_all,'psir',subname)

  call deallocate_work_arrays_locham(lr,wrk_lh)

END SUBROUTINE local_hamiltonian_Linear


subroutine cubic_exact_exchange(iproc,nproc,nspin,npsidim,size_potxc,hx,hy,hz,Glr,orbs,&
           ngatherarr,psi,potxc,eexctX,pkernel,orbsocc,psirocc)
  use module_base
  use module_types
  use module_xc 
  implicit none
  integer, intent(in) :: iproc,nproc,nspin,npsidim,size_potxc
  real(gp), intent(in) :: hx,hy,hz
  type(locreg_descriptors) :: Glr
  type(orbitals_data) :: orbs
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
  real(wp), dimension(npsidim), intent(in) :: psi
  real(wp), dimension(size_potxc),intent(out) :: potxc
  real(gp), intent(out) :: eexctX
  real(dp), dimension(*), optional :: pkernel
  type(orbitals_data), intent(in), optional :: orbsocc
  real(wp), dimension(:), pointer, optional :: psirocc
!  local variables
  character(len=*), parameter :: subname='cubic_exact_exchange'
  logical :: exctX,op2p
  integer :: n3p

  !initialise exact exchange energy 
  op2p=(eexctX ==  UNINITIALIZED(1.0_gp))
  eexctX=0.0_gp

  exctX = xc_exctXfac() /= 0.0_gp

  !fill the rest of the potential with the exact-exchange terms
  if (present(pkernel) .and. exctX ) then
     n3p=ngatherarr(iproc,1)/(Glr%d%n1i*Glr%d%n2i)

     !exact exchange for virtual orbitals (needs psirocc)
     !here we have to add the round part
     if (present(psirocc) .and. present(orbsocc)) then
        call exact_exchange_potential_virt(iproc,nproc,Glr%geocode,nspin,&
             Glr,orbsocc,orbs,ngatherarr(0,1),n3p,&
             0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,pkernel,psirocc,psi,potxc)
        eexctX = 0._gp
     else
 !!$        call exact_exchange_potential_round(iproc,nproc,at%geocode,nspin,lr,orbs,&
 !!$             0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,pkernel,psi,pot(ispot),eexctX)

        !here the condition for the scheme should be chosen
        if (.not. op2p) then
           call exact_exchange_potential(iproc,nproc,Glr%geocode,nspin,&
                Glr,orbs,ngatherarr(0,1),n3p,&
                0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,pkernel,psi,potxc,eexctX)
        else
           !the psi should be transformed in real space
           call exact_exchange_potential_round(iproc,nproc,Glr%geocode,nspin,&
                Glr,orbs,0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,pkernel,&
                psi,potxc,eexctX)

        end if
     end if
  else
     eexctX = 0._gp
     !print *,'iproc,eexctX',iproc,eexctX
  end if

end subroutine cubic_exact_exchange








subroutine apply_confinement(iproc, n1, n2, n3, nl1, nl2, nl3, nbuf, nspinor, psir, &
     rxyzConfinement, hxh, hyh, hzh, potentialPrefac, confPotOrder, offsetx, offsety, offsetz, &
     ibyyzz_r) !optional
use module_base
implicit none
integer, intent(in) :: iproc, n1,n2,n3,nl1,nl2,nl3,nbuf,nspinor, confPotOrder, offsetx, offsety, offsetz
real(wp), dimension(-14*nl1:2*n1+1+15*nl1,-14*nl2:2*n2+1+15*nl2,-14*nl3:2*n3+1+15*nl3,nspinor), intent(inout) :: psir
integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r
real(8),dimension(3),intent(in):: rxyzConfinement
real(8),intent(in):: hxh, hyh, hzh, potentialPrefac
!local variables
integer :: i1,i2,i3,i1s,i1e,ispinor, order
real(wp) :: tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42,tt
real(wp) :: psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4
real(gp) :: epot_p, epot


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
  

!!!$omp parallel default(private)&
!!!$omp shared(psir,n1,n2,n3,epot,ibyyzz_r,nl1,nl2,nl3,nbuf,nspinor)
  !case without bounds
  i1s=-14*nl1
  i1e=2*n1+1+15*nl1
  !epot_p=0._gp
!!!$omp do
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
              !write(*,'(a,5i8)') 'iproc, i1, i2, i1s, i1e', iproc, i1, i2, i1s, i1e
              
              !here we put the branchments wrt to the spin
              if (nspinor == 4) then
                 stop 'this part is not yet implemented'
              else
                 do ispinor=1,nspinor
                    do i1=i1s,i1e
                       ! THIS IS CORRECT #################################################################
                       tt=(hxh*dble(i1+offsetx)-rxyzConfinement(1))**2 + (hyh*dble(i2+offsety)-rxyzConfinement(2))**2 + &
                          (hzh*dble(i3+offsetz)-rxyzConfinement(3))**2
                       !!!!! EXPERIMENTAL ########################
                       !!tt=sqrt(tt)
                       !!tt=max(tt-3.d0,0.d0)
                       !!tt=tt**2
                       !!!!########################################
                       tt=potentialPrefac*tt**order

                       !!! New trial
                       !!tt=(hxh*dble(i1+offsetx)-rxyzConfinement(1))**2 + (hyh*dble(i2+offsety)-rxyzConfinement(2))**2 + &
                       !!   (hzh*dble(i3+offsetz)-rxyzConfinement(3))**2
                       !!tt=sqrt(tt)
                       !!tt=tt/5.d0
                       !!tt=tt**2
                       !!tt=5.d-2*(exp(tt)-1.d0)




                       ! #################################################################################
                       !!tt=(hxh*dble(i1+offsetx)-rxyzConfinement(1))**2 + (hyh*dble(i2+offsety)-rxyzConfinement(2))**2 + &
                       !!   (hzh*dble(i3+offsetz)-rxyzConfinement(3))**2
                       !!tt=(hxh*dble(i1+offsetx)-rxyzConfinement(1))**4 - 2.d0*2.d0*(hxh*dble(i1+offsetx)-rxyzConfinement(1))**2 + 2.d0**2
                       !!tt=tt*potentialPrefac
                       !!if((hxh*dble(i1+offsetx)-rxyzConfinement(1))**2 + (hyh*dble(i2+offsety)-rxyzConfinement(2))**2 + (hzh*dble(i3+offsetz)-rxyzConfinement(3))**2 > 2.d0) then
                       !!else
                       !!    tt=0.d0
                       !!end if
                       !!tt=(hxh*dble(i1+offsetx)-rxyzConfinement(1))**2 + (hyh*dble(i2+offsety)-rxyzConfinement(2))**2 + &
                       !!   (hzh*dble(i3+offsetz)-rxyzConfinement(3))**2
                       !!tt=potentialPrefac*tt**3
                       !!tt = (hxh*dble(i1+offsetx)-rxyzConfinement(1))**6 + (hyh*dble(i2+offsety)-rxyzConfinement(2))**6 + &
                       !!     3.d0*(hxh*dble(i1+offsetx)-rxyzConfinement(1))**4*(hyh*dble(i2+offsety)-rxyzConfinement(2))**2 + &
                       !!     3.d0*(hxh*dble(i1+offsetx)-rxyzConfinement(1))**2*(hyh*dble(i2+offsety)-rxyzConfinement(2))**4
                       !!tt=potentialPrefac*tt
                       tt=tt*psir(i1,i2,i3,ispinor)
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
!!!$omp end do

!!!$omp end parallel

END SUBROUTINE apply_confinement
!!***
