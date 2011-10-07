!> @file
!!  Routine to calculate the action of the hamiltonian
!! @author
!!   Copyright (C) 2005-2011 BigDFT group 
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 


!>   Calculate the action of the local hamiltonian on the orbitals
!! @param ipotmethod Indicates the method which has to be chosen for applying the potential to the wavefunctions in the 
!!                   real space form:
!!                   0 is the traditional potential application
!!                   1 is the application of the exact exchange (which has to be precomputed and stored in the potential array)
!!                   2 is the application of the Perdew-Zunger SIC
!!                   3 is the application of the Non-Koopman's correction SIC
subroutine local_hamiltonian(iproc,orbs,lr,hx,hy,hz,&
     ipotmethod,pot,psi,hpsi,pkernel,ixc,alphaSIC,ekin_sum,epot_sum,eSIC_DC)
  use module_base
  use module_types
  use module_interfaces, except_this_one => local_hamiltonian
  use module_xc
  implicit none
  integer, intent(in) :: iproc,ipotmethod,ixc
  real(gp), intent(in) :: hx,hy,hz,alphaSIC
  type(orbitals_data), intent(in) :: orbs
  type(locreg_descriptors), intent(in) :: lr
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*orbs%norbp), intent(in) :: psi
  real(wp), dimension(*) :: pot !< the potential, with the dimension compatible with the ipotmethod flag
  !real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i*nspin) :: pot
  real(gp), intent(out) :: ekin_sum,epot_sum,eSIC_DC
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*orbs%norbp), intent(out) :: hpsi
  real(dp), dimension(:), pointer :: pkernel !< the PSolver kernel which should be associated for the SIC schemes
  !local variables
  character(len=*), parameter :: subname='local_hamiltonian'
  integer :: i_all,i_stat,iorb,npot,nsoffset,oidx,ispot
  real(wp) :: exctXcoeff
  real(gp) :: ekin,epot,kx,ky,kz,etest,eSICi,eSIC_DCi
  type(workarr_locham) :: wrk_lh
  real(wp), dimension(:,:), allocatable :: psir,vsicpsir

  !some checks
  exctXcoeff=xc_exctXfac()

  if (exctXcoeff /= 0.0_gp .neqv. ipotmethod ==1) then
     if (iproc==0) write(*,*)&
          'ERROR (local_hamiltonian): potential method not compatible with exact exchange'
     stop
  end if

  if (.not.(associated(pkernel) .and. alphaSIC /=0.0_gp) .and. ipotmethod == 2) then
     if (iproc==0) write(*,*)&
          'ERROR (local_hamiltonian): potential method not compatible with SIC'
     stop
  end if

  !initialise the work arrays
  call initialize_work_arrays_locham(lr,orbs%nspinor,wrk_lh)  

  !components of the potential
  npot=orbs%nspinor
  if (orbs%nspinor == 2) npot=1

  ! Wavefunction in real space
  allocate(psir(lr%d%n1i*lr%d%n2i*lr%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
  call memocc(i_stat,psir,'psir',subname)

  call razero(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspinor,psir)

  if (ipotmethod == 2 .or. ipotmethod == 3) then! wavefunction after application of the self-interaction potential
     allocate(vsicpsir(lr%d%n1i*lr%d%n2i*lr%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
     call memocc(i_stat,vsicpsir,'vsicpsir',subname)
  end if


  ekin_sum=0.0_gp
  epot_sum=0.0_gp
  eSIC_DC=0.0_gp

  etest=0.0_gp

  do iorb=1,orbs%norbp
     !this section should be replaced with the ispot array, calculated in fill_local_potential (or maybe before)

     if(orbs%spinsgn(iorb+orbs%isorb)>0.0_gp .or. orbs%nspin == 1 .or. orbs%nspin == 4 ) then
        nsoffset=1
     else
        nsoffset=lr%d%n1i*lr%d%n2i*lr%d%n3i+1
     end if

     oidx=(iorb-1)*orbs%nspinor+1

     !transform the wavefunction in Daubechies basis to the wavefunction in ISF basis
     !the psir wavefunction is given in the spinorial form
     call daub_to_isf_locham(orbs%nspinor,lr,wrk_lh,psi(1,oidx),psir)
     !ispot=1+lr%d%n1i*lr%d%n2i*lr%d%n3i*(nspin+iorb-1)
     !etest=etest+dot(lr%d%n1i*lr%d%n2i*lr%d%n3i,pot(ispot),1,psir(1,1),1)
     !print *,'epot, iorb,iproc,norbp',iproc,orbs%norbp,iorb,etest

     !Perdew-Zunger SIC scheme
     eSIC_DCi=0.0_gp
     if (ipotmethod == 2) then
        call PZ_SIC_potential(iorb,lr,orbs,ixc,0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,pkernel,psir,vsicpsir,eSICi,eSIC_DCi)
     else if (ipotmethod == 3) then !NonKoopmans' correction scheme
        call vcopy(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspinor,psir(1,1),1,vsicpsir(1,1),1)
        call psir_to_vpsi(npot,orbs%nspinor,lr,&
             pot(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspin+(iorb-1)*lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspinor+1),vsicpsir,eSICi)
     end if

     !apply the potential to the psir wavefunction and calculate potential energy
     call psir_to_vpsi(npot,orbs%nspinor,lr,pot(nsoffset),psir,epot)

!!$     select case(lr%geocode)
!!$     case('F')
!!$
!!$        call apply_potential(lr%d%n1,lr%d%n2,lr%d%n3,1,1,1,0,orbs%nspinor,npot,psir,&
!!$             pot(nsoffset),epot,&
!!$             lr%bounds%ibyyzz_r) !optional
!!$          
!!$     case('P') 
!!$        !here the hybrid BC act the same way
!!$        call apply_potential(lr%d%n1,lr%d%n2,lr%d%n3,0,0,0,0,orbs%nspinor,npot,psir,&
!!$             pot(nsoffset),epot)
!!$
!!$     case('S')
!!$
!!$        call apply_potential(lr%d%n1,lr%d%n2,lr%d%n3,0,1,0,0,orbs%nspinor,npot,psir,&
!!$             pot(nsoffset),epot)
!!$     end select

     if (ipotmethod==1) then
        ispot=1+lr%d%n1i*lr%d%n2i*lr%d%n3i*(orbs%nspin+iorb-1)
        !add to the psir function the part of the potential coming from the exact exchange
        call axpy(lr%d%n1i*lr%d%n2i*lr%d%n3i,exctXcoeff,pot(ispot),1,psir(1,1),1)
     else if (ipotmethod == 2) then
        !subtract the sic potential from the vpsi function
        call axpy(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspinor,-alphaSIC,vsicpsir(1,1),1,psir(1,1),1)
        !add the SIC caorrection to the potential energy
        epot=epot-alphaSIC*eSICi
        !accumulate the Double-Counted SIC energy
        eSIC_DC=eSIC_DC+alphaSIC*eSIC_DCi
     else if (ipotmethod == 3) then
        !add the sic potential from the vpsi function
        call axpy(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspinor,alphaSIC,vsicpsir(1,1),1,psir(1,1),1)
        epot=epot+alphaSIC*eSICi
        !accumulate the Double-Counted SIC energy
        eSIC_DC=eSIC_DC+alphaSIC*orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*eSICi
     end if

     !apply the kinetic term, sum with the potential and transform back to Daubechies basis
     !k-point values, if present
     kx=orbs%kpts(1,orbs%iokpt(iorb))
     ky=orbs%kpts(2,orbs%iokpt(iorb))
     kz=orbs%kpts(3,orbs%iokpt(iorb))
     call isf_to_daub_kinetic(hx,hy,hz,kx,ky,kz,orbs%nspinor,lr,wrk_lh,&
          psir,hpsi(1,oidx),ekin)

     !print *,iorb, ekin, epot,orbs%kwgts(orbs%iokpt(iorb))
     ekin_sum=ekin_sum+orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*ekin
     epot_sum=epot_sum+orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*epot
  enddo

  !deallocations of work arrays
  i_all=-product(shape(psir))*kind(psir)
  deallocate(psir,stat=i_stat)
  call memocc(i_stat,i_all,'psir',subname)
  if (ipotmethod == 2 .or. ipotmethod ==3) then
     i_all=-product(shape(vsicpsir))*kind(vsicpsir)
     deallocate(vsicpsir,stat=i_stat)
     call memocc(i_stat,i_all,'vsicpsir',subname)
  end if
  call deallocate_work_arrays_locham(lr,wrk_lh)

END SUBROUTINE local_hamiltonian


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
  call dcopy(nspinor*n1*n2*n3,ww,1,x,1)
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


subroutine realspace(ibyyzz_r,pot,psir,epot,n1,n2,n3)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in) :: ibyyzz_r
  real(wp), dimension(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16), intent(in) :: pot
  real(wp), dimension(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16), intent(inout) :: psir
  real(wp), intent(out) :: epot
  !local variables
  real(wp) :: tt
  integer :: i1,i2,i3

  epot=0.0_wp
  do i3=-14,2*n3+16
     do i2=-14,2*n2+16
        do i1=max(ibyyzz_r(1,i2,i3)-14,-14),min(ibyyzz_r(2,i2,i3)-14,2*n1+16)
           tt=pot(i1,i2,i3)*psir(i1,i2,i3)
           epot=epot+tt*psir(i1,i2,i3)
           psir(i1,i2,i3)=tt
        enddo
     enddo
  enddo

END SUBROUTINE realspace


subroutine realspace_nbuf(ibyyzz_r,pot,psir,epot,nb1,nb2,nb3,nbuf)
  implicit none
  !Arguments
  integer,intent(in)::nb1,nb2,nb3,nbuf
  integer,intent(in)::ibyyzz_r(2,-14:2*nb2+16,-14:2*nb3+16)
  real(kind=8),intent(in)::pot(-14:2*nb1+16-4*nbuf,-14:2*nb2+16-4*nbuf,-14:2*nb3+16-4*nbuf)
  real(kind=8),intent(inout)::psir(-14:2*nb1+16,-14:2*nb2+16,-14:2*nb3+16)
  real(kind=8),intent(out)::epot
  !Local variables
  real(kind=8) :: tt
  integer :: i1,i2,i3

  epot=0.d0
  do i3=-14,2*nb3+16
     if (i3 >= -14+2*nbuf .and. i3 <= 2*nb3+16-2*nbuf) then
        do i2=-14,2*nb2+16
           if (i2 >= -14+2*nbuf .and. i2 <= 2*nb2+16-2*nbuf) then
              do i1=-14+2*nbuf,ibyyzz_r(1,i2,i3)-14-1
                 psir(i1,i2,i3)=0.d0
              enddo
              do i1=max(ibyyzz_r(1,i2,i3)-14,-14+2*nbuf),min(ibyyzz_r(2,i2,i3)-14,2*nb1+16-2*nbuf)
                 tt=pot(i1-2*nbuf,i2-2*nbuf,i3-2*nbuf)*psir(i1,i2,i3)
                 epot=epot+tt*psir(i1,i2,i3)
                 psir(i1,i2,i3)=tt
              enddo
              do i1=ibyyzz_r(2,i2,i3)-14+1,2*nb1+16-2*nbuf
                 psir(i1,i2,i3)=0.d0
              enddo
           else
              do i1=-14,2*nb1+16
                 psir(i1,i2,i3)=0.d0
              enddo
           endif
        enddo
     else
        do i2=-14,2*nb2+16
           do i1=-14,2*nb1+16
              psir(i1,i2,i3)=0.d0
           enddo
        enddo
     endif
  enddo

END SUBROUTINE realspace_nbuf


subroutine realspaceINOUT(ibyyzz_r,pot,psirIN,psirOUT,epot,n1,n2,n3)
  implicit none
  integer,intent(in)::n1,n2,n3
  integer,intent(in)::ibyyzz_r(2,-14:2*n2+16,-14:2*n3+16)

  real(kind=8),intent(in)::pot(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16)
  real(kind=8),intent(in)::psirIN(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16)
 real(kind=8),intent(out)::psirOUT(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16)

  real(kind=8),intent(out)::epot
  real(kind=8) tt
  integer i1,i2,i3

  epot=0.d0
  do i3=-14,2*n3+16
     do i2=-14,2*n2+16
        do i1=max(ibyyzz_r(1,i2,i3)-14,-14),min(ibyyzz_r(2,i2,i3)-14,2*n1+16)
           tt=pot(i1,i2,i3)*psirIN(i1,i2,i3)
           epot=epot+tt*psirIN(i1,i2,i3)
           psirOUT(i1,i2,i3)=psirOUT(i1,i2,i3)+tt
        enddo
     enddo
  enddo

END SUBROUTINE realspaceINOUT


subroutine realspaceINOUT_nbuf(ibyyzz_r,pot,psirIN,psirOUT,epot,nb1,nb2,nb3,nbuf)
  implicit none
  !Arguments
  integer,intent(in) :: nb1,nb2,nb3,nbuf
  integer,intent(in) :: ibyyzz_r(2,-14:2*nb2+16,-14:2*nb3+16)
  real(kind=8),intent(in) :: pot(-14:2*nb1+16-4*nbuf,-14:2*nb2+16-4*nbuf,-14:2*nb3+16-4*nbuf)
  real(kind=8),intent(in) :: psirIN(-14:2*nb1+16,-14:2*nb2+16,-14:2*nb3+16)
  real(kind=8),intent(out) :: psirOUT(-14:2*nb1+16,-14:2*nb2+16,-14:2*nb3+16)
  real(kind=8),intent(out) :: epot
  !Local variables
  real(kind=8) :: tt
  integer :: i1,i2,i3

  epot=0.d0
  do i3=-14,2*nb3+16
     if (i3.ge.-14+2*nbuf .and. i3.le.2*nb3+16-2*nbuf) then
        do i2=-14,2*nb2+16
           if (i2.ge.-14+2*nbuf .and. i2.le.2*nb2+16-2*nbuf) then
              do i1=-14+2*nbuf,ibyyzz_r(1,i2,i3)-14-1
                 psirOUT(i1,i2,i3)=0.d0
              enddo
              do i1=max(ibyyzz_r(1,i2,i3)-14,-14+2*nbuf),min(ibyyzz_r(2,i2,i3)-14,2*nb1+16-2*nbuf)
                 tt=pot(i1-2*nbuf,i2-2*nbuf,i3-2*nbuf)*psirIN(i1,i2,i3)
                 epot=epot+tt*psirIN(i1,i2,i3)
                 psirOUT(i1,i2,i3)=tt
              enddo
              do i1=ibyyzz_r(2,i2,i3)-14+1,2*nb1+16-2*nbuf
                 psirOUT(i1,i2,i3)=0.d0
              enddo
           else
              do i1=-14,2*nb1+16
                 psirOUT(i1,i2,i3)=0.d0
              enddo
           endif
        enddo
     else
        do i2=-14,2*nb2+16
           do i1=-14,2*nb1+16
              psirOUT(i1,i2,i3)=0.d0
           enddo
        enddo
     endif
  enddo

END SUBROUTINE realspaceINOUT_nbuf


subroutine realspaceINPLACE(ibyyzz_r,pot,psir,epot,n1,n2,n3)
  implicit none
  integer,intent(in)::n1,n2,n3
  integer,intent(in)::ibyyzz_r(2,-14:2*n2+16,-14:2*n3+16)

  real(kind=8),intent(in)::pot(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16,4)
  real(kind=8),intent(inout)::psir(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16,4)

  real(kind=8),intent(out)::epot
  real(kind=8) tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42
  integer i1,i2,i3

  epot=0.d0
  do i3=-14,2*n3+16
     do i2=-14,2*n2+16
        do i1=max(ibyyzz_r(1,i2,i3)-14,-14),min(ibyyzz_r(2,i2,i3)-14,2*n1+16)
           !diagonal terms
           tt11=pot(i1,i2,i3,1)*psir(i1,i2,i3,1) !p1
           tt22=pot(i1,i2,i3,1)*psir(i1,i2,i3,2) !p2
           tt33=pot(i1,i2,i3,4)*psir(i1,i2,i3,3) !p3
           tt44=pot(i1,i2,i3,4)*psir(i1,i2,i3,4) !p4
           !Rab*Rb
           tt13=pot(i1,i2,i3,2)*psir(i1,i2,i3,3) !p1
           !Iab*Ib
           tt14=pot(i1,i2,i3,3)*psir(i1,i2,i3,4) !p1
           !Rab*Ib
           tt23=pot(i1,i2,i3,2)*psir(i1,i2,i3,4) !p2
           !Iab*Rb
           tt24=pot(i1,i2,i3,3)*psir(i1,i2,i3,3) !p2
           !Rab*Ra
           tt31=pot(i1,i2,i3,2)*psir(i1,i2,i3,1) !p3
           !Iab*Ia
           tt32=pot(i1,i2,i3,3)*psir(i1,i2,i3,2) !p3
           !Rab*Ia
           tt41=pot(i1,i2,i3,2)*psir(i1,i2,i3,2) !p4
           !Iab*Ra
           tt42=pot(i1,i2,i3,3)*psir(i1,i2,i3,1) !p4
           ! Change epot later
           epot=epot+tt11*psir(i1,i2,i3,1)+tt22*psir(i1,i2,i3,2)+tt33*psir(i1,i2,i3,3)+tt44*psir(i1,i2,i3,4)+&
                2.0d0*tt31*psir(i1,i2,i3,3)-2.0d0*tt42*psir(i1,i2,i3,4)+2.0d0*tt41*psir(i1,i2,i3,4)+2.0d0*tt32*psir(i1,i2,i3,3)
!p1=h1p1+h2p3-h3p4
!p2=h1p2+h2p4+h3p3
!p3=h2p1+h3p2+h4p3
!p4=h2p2-h3p1+h4p4
           psir(i1,i2,i3,1)=tt11+tt13-tt14
           psir(i1,i2,i3,2)=tt22+tt23+tt24
           psir(i1,i2,i3,3)=tt33+tt31+tt32
           psir(i1,i2,i3,4)=tt44+tt41-tt42
        enddo
     enddo
  enddo

END SUBROUTINE realspaceINPLACE


!>   Calculate on-the fly each projector for each atom, then applies the projectors 
!!   to all distributed orbitals
subroutine applyprojectorsonthefly(iproc,orbs,at,n1,n2,n3,&
     rxyz,hx,hy,hz,wfd,nlpspd,proj,psi,hpsi,eproj_sum)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,n1,n2,n3
  real(gp), intent(in) :: hx,hy,hz
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(nonlocal_psp_descriptors), intent(in) :: nlpspd
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(wp), dimension((wfd%nvctr_c+7*wfd%nvctr_f)*orbs%nspinor*orbs%norbp), intent(in) :: psi
  real(wp), dimension((wfd%nvctr_c+7*wfd%nvctr_f)*orbs%nspinor*orbs%norbp), intent(inout) :: hpsi
  real(gp), intent(out) :: eproj_sum
  real(wp), dimension(nlpspd%nprojel), intent(out) :: proj
  !local variables
  integer :: iat,nwarnings,iproj,iorb
  integer :: istart_c,idir,isorb,ieorb,ikpt,nspinor,ispsi_k,ispsi
  
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
     do iat=1,at%nat
        istart_c=1
        call atom_projector(ikpt,iat,idir,istart_c,iproj,&
             n1,n2,n3,hx,hy,hz,rxyz,at,orbs,nlpspd,proj,nwarnings)

        !apply the projector to all the orbitals belonging to the processor
        ispsi=ispsi_k
        do iorb=isorb,ieorb
           istart_c=1
           call apply_atproj_iorb(iat,iorb,istart_c,at,orbs,wfd,nlpspd,&
                proj,psi(ispsi),hpsi(ispsi),eproj_sum)
           ispsi=ispsi+(wfd%nvctr_c+7*wfd%nvctr_f)*nspinor
        end do

     end do
     if (iproj /= nlpspd%nproj) stop 'incorrect number of projectors created'
     if (ieorb == orbs%norbp) exit loop_kpt
     ikpt=ikpt+1
     ispsi_k=ispsi
  end do loop_kpt

  if (iproc == 0 .and. nlpspd%nproj /=0 .and. idir == 0) then
     if (nwarnings == 0) then
     else
        write(*,'(1x,a,i0,a)')'found ',nwarnings,' warnings.'
        write(*,'(1x,a)')'Some projectors may be too rough.'
        write(*,'(1x,a,f6.3)')&
             'Consider the possibility of reducing hgrid for having a more accurate run.'
     end if
  end if

END SUBROUTINE applyprojectorsonthefly


!>   Applies the projector associated on a given atom on a corresponding orbital
subroutine apply_atproj_iorb(iat,iorb,istart_c,at,orbs,wfd,nlpspd,proj,psi,hpsi,eproj)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iat,iorb
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(nonlocal_psp_descriptors), intent(in) :: nlpspd
  real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,orbs%nspinor), intent(in) :: psi
  integer, intent(inout) :: istart_c
  real(gp), intent(inout) :: eproj
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,orbs%nspinor), intent(inout) :: hpsi
  !local variables
  integer :: ispinor,ityp,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,jseg_c,l,i,istart_c_i,ncplx
  real(gp) :: eproj_spinor

  !complex functions or not
  !this should be decided as a function of the orbital
  !features of the k-point ikpt
  call ncplx_kpt(orbs%iokpt(iorb),orbs,ncplx)

  istart_c_i=istart_c
  do ispinor=1,orbs%nspinor,ncplx
     eproj_spinor=0.0_gp
     if (ispinor >= 2) istart_c=istart_c_i
     ityp=at%iatype(iat)
     mbvctr_c=nlpspd%nvctr_p(2*iat-1)-nlpspd%nvctr_p(2*iat-2)
     mbvctr_f=nlpspd%nvctr_p(2*iat  )-nlpspd%nvctr_p(2*iat-1)
     
     mbseg_c=nlpspd%nseg_p(2*iat-1)-nlpspd%nseg_p(2*iat-2)
     mbseg_f=nlpspd%nseg_p(2*iat  )-nlpspd%nseg_p(2*iat-1)
     jseg_c=nlpspd%nseg_p(2*iat-2)+1
     !GTH and HGH pseudopotentials
     do l=1,4
        do i=1,3
           if (at%psppar(l,i,ityp) /= 0.0_gp) then
              call applyprojector(ncplx,l,i,at%psppar(0,0,ityp),at%npspcode(ityp),&
                   wfd%nvctr_c,wfd%nvctr_f,wfd%nseg_c,wfd%nseg_f,wfd%keyv,wfd%keyg,&
                   mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
                   nlpspd%keyv_p(jseg_c),nlpspd%keyg_p(1,jseg_c),proj(istart_c),&
                   psi(1,ispinor),hpsi(1,ispinor),eproj_spinor)
              istart_c=istart_c+(mbvctr_c+7*mbvctr_f)*(2*l-1)*ncplx
           end if
        enddo
     enddo
     eproj=eproj+&
          orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*eproj_spinor
  end do
END SUBROUTINE apply_atproj_iorb


subroutine applyprojector(ncplx,l,i,psppar,npspcode,&
     nvctr_c,nvctr_f,nseg_c,nseg_f,keyv,keyg,&
     mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p,keyg_p,proj,psi,hpsi,eproj)
  use module_base
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
  if ((npspcode == 3 .and. l/=4 .and. i/=3) .or. &
       (npspcode == 10 .and. i/=3)) then !HGH(-K) case, offdiagonal terms
     loop_j: do j=i+1,3
        if (psppar(l,j) == 0.0_gp) exit loop_j

        !offdiagonal HGH term
        if (npspcode == 3) then !traditional HGH convention
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


subroutine applyprojector_old(l,i,psppar,npspcode,&
     nvctr_c,nvctr_f,nseg_c,nseg_f,keyv,keyg,&
     mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p,keyg_p,proj,psi,hpsi,eproj)
  use module_base
  implicit none
  integer, intent(in) :: i,l,npspcode
  integer, intent(in) :: nvctr_c,nvctr_f,nseg_c,nseg_f,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  integer, dimension(mbseg_c+mbseg_f), intent(in) :: keyv_p
  integer, dimension(2,mbseg_c+mbseg_f), intent(in) :: keyg_p
  real(wp), dimension(*), intent(in) :: proj
  real(gp), dimension(0:4,0:6), intent(in) :: psppar
  real(wp), dimension(nvctr_c+7*nvctr_f), intent(in) :: psi
  real(gp), intent(inout) :: eproj
  real(wp), dimension(nvctr_c+7*nvctr_f), intent(inout) :: hpsi
  !local variables
  integer :: j,m,istart_c,istart_f,istart_c_i,istart_c_j,istart_f_i,istart_f_j
  real(dp) :: scpr,scprp,scpr_i,scprp_i,scpr_j,scprp_j
  real(gp) :: offdiagcoeff
  real(gp), dimension(2,2,3) :: offdiagarr

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
     istart_f=istart_c+mbvctr_c

     call wpdot_wrap(1,  &
          nvctr_c,nvctr_f,nseg_c,nseg_f,keyv,keyg,psi,  &
          mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p,keyg_p,proj(istart_c),scpr)
  
     scprp=scpr*real(psppar(l,i),dp)
     eproj=eproj+real(scprp,gp)*real(scpr,gp)

     call waxpy_wrap(1,scprp,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p,keyg_p,proj(istart_c),&
          nvctr_c,nvctr_f,nseg_c,nseg_f,keyv,keyg,hpsi)

     !print *,'scprp,m,l,i',scprp,m,l,i

     istart_c=istart_f+7*mbvctr_f
  enddo
  if (npspcode == 3 .and. l/=4 .and. i/=3) then !HGH case, offdiagonal terms
     loop_j: do j=i+1,3
        if (psppar(l,j) == 0.0_gp) exit loop_j
        !calculate the coefficients for the off-diagonal terms
        if (l==1) then
           if (i==1) then
              if (j==2) offdiagcoeff=-0.5_gp*sqrt(3._gp/5._gp)
              if (j==3) offdiagcoeff=0.5_gp*sqrt(5._gp/21._gp)
           else
              offdiagcoeff=-0.5_gp*sqrt(100._gp/63._gp)
           end if
        else if (l==2) then
           if (i==1) then
              if (j==2) offdiagcoeff=-0.5_gp*sqrt(5._gp/7._gp)
              if (j==3) offdiagcoeff=1._gp/6._gp*sqrt(35._gp/11._gp)
           else
              offdiagcoeff=-7._gp/3._gp*sqrt(1._gp/11._gp)
           end if
        else if (l==3) then
           if (i==1) then
              if (j==2) offdiagcoeff=-0.5_gp*sqrt(7._gp/9._gp)
              if (j==3) offdiagcoeff=0.5_gp*sqrt(63._gp/143._gp)
           else
              offdiagcoeff=-9._gp*sqrt(1._gp/143._gp)
           end if
        end if
        istart_c_i=istart_c-(2*l-1)*(mbvctr_c+7*mbvctr_f)
        istart_c_j=istart_c_i+(j-i)*(2*l-1)*(mbvctr_c+7*mbvctr_f)
        do m=1,2*l-1
           !starting addresses of the projectors
           istart_f_j=istart_c_j+mbvctr_c
           istart_f_i=istart_c_i+mbvctr_c
           call wpdot_wrap(1,nvctr_c,nvctr_f,nseg_c,nseg_f,keyv,keyg,psi,  &
                mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p,keyg_p,proj(istart_c_j),scpr_j)

           call wpdot_wrap(1,nvctr_c,nvctr_f,nseg_c,nseg_f,keyv,keyg,psi,  &
                mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p,keyg_p,proj(istart_c_i),scpr_i)

           scprp_j=scpr_j*real(offdiagcoeff*psppar(l,j),dp)
           scprp_i=scpr_i*real(offdiagcoeff*psppar(l,j),dp)
           !scpr_i*h_ij*scpr_j+scpr_j*h_ij*scpr_i
           eproj=eproj+2._gp*real(scpr_j,gp)*&
                offdiagcoeff*psppar(l,j)*real(scpr_i,gp)

           !|hpsi>=|hpsi>+h_ij (<p_i|psi>|p_j>+<p_j|psi>|p_i>)
           call waxpy_wrap(1,scprp_j,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p,keyg_p,&
                proj(istart_c_i),&
                nvctr_c,nvctr_f,nseg_c,nseg_f,keyv,keyg,hpsi)

           call waxpy_wrap(1,scprp_i,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
                keyv_p,keyg_p,proj(istart_c_j),&
                nvctr_c,nvctr_f,nseg_c,nseg_f,keyv,keyg,hpsi)

           istart_c_j=istart_f_j+7*mbvctr_f
           istart_c_i=istart_f_i+7*mbvctr_f
        enddo
     end do loop_j
  else if (npspcode == 10 .and. i/=3) then !HGH-K case, offdiagonal terms
     loop_jK: do j=i+1,3
        if (psppar(l,j) .eq. 0._gp) exit loop_jK
        istart_c_i=istart_c-(2*l-1)*(mbvctr_c+7*mbvctr_f)
        istart_c_j=istart_c_i+(j-i)*(2*l-1)*(mbvctr_c+7*mbvctr_f)
        do m=1,2*l-1
           !starting addresses of the projectors
           istart_f_j=istart_c_j+mbvctr_c
           istart_f_i=istart_c_i+mbvctr_c
           call wpdot_wrap(1,nvctr_c,nvctr_f,nseg_c,nseg_f,keyv,keyg,psi,&
                mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p,keyg_p,proj(istart_c_j),scpr_j)

           call wpdot_wrap(1,nvctr_c,nvctr_f,nseg_c,nseg_f,keyv,keyg,psi,&
                mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p,keyg_p,proj(istart_c_i),scpr_i)

           !scpr_i*h_ij*scpr_j+scpr_j*h_ij*scpr_i (with symmetric h_ij)
           eproj=eproj+2._gp*real(scpr_i,gp)*psppar(l,i+j+1)*real(scpr_j,gp)
           scprp_j=scpr_j*real(psppar(l,i+j+1),dp)
           scprp_i=scpr_i*real(psppar(l,i+j+1),dp)

           !|hpsi>=|hpsi>+h_ij (<p_i|psi>|p_j>+<p_j|psi>|p_i>)
           call waxpy_wrap(1,scprp_j,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p,keyg_p,&
                proj(istart_c_i),&
                nvctr_c,nvctr_f,nseg_c,nseg_f,keyv,keyg,hpsi)

           call waxpy_wrap(1,scprp_i,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keyv_p,keyg_p,&
                proj(istart_c_j),&
                nvctr_c,nvctr_f,nseg_c,nseg_f,keyv,keyg,hpsi)

           istart_c_j=istart_f_j+7*mbvctr_f
           istart_c_i=istart_f_i+7*mbvctr_f
        enddo
     end do loop_jK
  end if
END SUBROUTINE applyprojector_old


!>   Find the starting and ending orbital for kpoint ikpt, and the corresponding nspinor
subroutine orbs_in_kpt(ikpt,orbs,isorb,ieorb,nspinor)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: ikpt
  type(orbitals_data), intent(in) :: orbs
  integer, intent(out) :: isorb,ieorb,nspinor
  !local variables
  integer :: iorb

  !find starting orbital
  do iorb=1,orbs%norbp
     if (orbs%iokpt(iorb)==ikpt) then
        isorb=iorb
        exit
     end if
  end do

  !find ending orbital
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

