!!subroutine check_cutoff(iproc, nproc, orbs, lzd, hx, hy, hz, locrad, confdatarr, psi)
!!use module_base
!!use module_types
!!use module_interfaces, except_this_one => apply_position_operators
!!implicit none
!!
!!! Calling arguments
!!integer,intent(in):: iproc, nproc
!!type(orbitals_data),intent(in):: orbs
!!type(local_zone_descriptors),intent(in):: lzd
!!real(8),intent(in):: hx, hy, hz, locrad
!!type(confpot_data),dimension(orbs%norbp),intent(in):: confdatarr
!!real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(in):: psi
!!
!!! Local variables
!!integer:: oidx, iorb, ilr, npot, icenter, i_stat, i_all, ist_c, ist_f, ist, iiorb, iall, ierr
!!real(8):: hxh, hyh, hzh, ddot, tt, t1, t2, time, weight_in, weight_out
!!real(8),dimension(:,:),allocatable:: psir, psirx, psiry, psirz
!!type(workarr_sumrho):: work_sr
!!real(8),dimension(0:3),parameter:: scal=1.d0
!!real(8),dimension(:,:,:),allocatable:: ypsitemp_c
!!real(8),dimension(:,:,:,:),allocatable:: ypsitemp_f
!!character(len=*),parameter:: subname='apply_orbitaldependent_potential'
!!integer, dimension(3) :: ishift !temporary variable in view of wavefunction creation
!!
!!  ishift=(/0,0,0/)
!!
!!  oidx = 0
!!  do iorb=1,orbs%norbp
!!     ilr = orbs%inwhichlocreg(iorb+orbs%isorb)
!!
!!  
!!     !initialise the work arrays
!!     call initialize_work_arrays_sumrho(lzd%llr(ilr), work_sr)
!!
!!     ! Wavefunction in real space
!!     allocate(psir(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
!!     call memocc(i_stat,psir,'psir',subname)
!!     call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psir)
!!
!!     call daub_to_isf(lzd%llr(ilr), work_sr, psi(1+oidx), psir)
!!     !apply the potential to the psir wavefunction and calculate potential energy
!!     !icenter=confinementCenter(iorb)
!!     !components of the potential
!!     npot=orbs%nspinor
!!     if (orbs%nspinor == 2) npot=1
!!
!!     write(*,'(a,2i8,3es16.7)') 'iproc, iorb, confdatarr(iorb)%rxyzConf', iproc, iorb, confdatarr(iorb)%rxyzConf
!!
!!     call get_cutoff_weight(lzd%llr(ilr)%d%n1i, lzd%llr(ilr)%d%n2i, lzd%llr(ilr)%d%n3i, &
!!                             lzd%llr(ilr)%d%n1i, lzd%llr(ilr)%d%n2i, lzd%llr(ilr)%d%n3i, &
!!                             ishift, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, orbs%nspinor, &
!!                             psir, locrad, weight_in, weight_out, &
!!                             confdatarr(iorb), lzd%llr(ilr)%bounds%ibyyzz_r) !optional
!!
!!     write(*,'(a,2i8,3es16.6)') 'iproc, iorb, weight_in, weight_out, ratio', &
!!         iproc, iorb, weight_in, weight_out, weight_in/(weight_in+weight_out)
!!
!!     i_all=-product(shape(psir))*kind(psir)
!!     deallocate(psir,stat=i_stat)
!!     call memocc(i_stat,i_all,'psir',subname)
!!
!!     call deallocate_work_arrays_sumrho(work_sr)
!!
!!     oidx = oidx + (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor
!!
!!  enddo
!!
!!
!!end subroutine check_cutoff
!!
!!
!!
!!subroutine get_cutoff_weight(n1i,n2i,n3i,n1ip,n2ip,n3ip,ishift,n2,n3,nspinor,psir,&
!!     cutoff, weight_in, weight_out, &
!!     confdata,ibyyzz_r) !optional
!!  use module_base
!!  use module_types
!!  implicit none
!!  integer, intent(in) :: n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,nspinor
!!  integer, dimension(3), intent(in) :: ishift !<offset of potential box in wfn box coords.
!!  real(wp), dimension(n1i,n2i,n3i,nspinor), intent(in) :: psir !< real-space wfn in lr
!!  real(8),intent(in):: cutoff
!!  real(8),intent(out):: weight_in, weight_out
!!  type(confpot_data), intent(in), optional :: confdata !< data for the confining potential
!!  integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r !< bounds in lr
!!  !local variables
!!  integer :: i1,i2,i3,ispinor,i1s,i1e,i2s,i2e,i3s,i3e,i1st,i1et
!!  real(wp) :: tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42,tt
!!  real(wp) :: psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4
!!  real(wp):: ttx, tty, ttz, cutoff2
!!
!!  ! Square of the cutoff radius
!!  cutoff2=cutoff**2
!!
!!  ! Initialize return values
!!  weight_in=0.d0
!!  weight_out=0.d0
!!
!!
!!  !loop on wavefunction
!!  !calculate the limits in all the directions
!!  !regions in which both the potential and wavefunctions are defined
!!  i3s=max(1,ishift(3)+1)
!!  i3e=min(n3i,n3ip+ishift(3))
!!  i2s=max(1,ishift(2)+1)
!!  i2e=min(n2i,n2ip+ishift(2))
!!  i1s=max(1,ishift(1)+1)
!!  i1e=min(n1i,n1ip+ishift(1))
!!
!!
!!  !$omp parallel default(none)&
!!  !$omp shared(psir,n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,ibyyzz_r,nspinor)&
!!  !$omp shared(i1s,i1e,i2s,i2e,i3s,i3e,ishift,confdata,weight_in,weight_out)&
!!  !$omp private(ispinor,i1,i2,i3,i1st,i1et)&
!!  !$omp private(tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42,tt)&
!!  !$omp private(psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4,ttx,tty,ttz,cutoff2)
!!
!!!!$  !$omp parallel default(private)&
!!!!$  !$omp shared(pot,psir,n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,ibyyzz_r,nspinor)&
!!!!$  !$omp shared(i1s,i1e,i2s,i2e,i3s,i3e,ishift)
!!  !case without bounds
!!
!!
!!  !!!put to zero the external part of psir if the potential is more little than the wavefunction
!!  !!!first part of the array
!!  !!do ispinor=1,nspinor
!!  !!   !$omp do 
!!  !!   do i3=1,i3s-1
!!  !!      do i2=1,n2i
!!  !!         do i1=1,n1i
!!  !!           psirx(i1,i2,i3,ispinor)=0.0_wp 
!!  !!           psiry(i1,i2,i3,ispinor)=0.0_wp 
!!  !!           psirz(i1,i2,i3,ispinor)=0.0_wp 
!!  !!         end do
!!  !!      end do
!!  !!   end do
!!  !!   !$omp end do
!!  !!end do
!!
!!  !!!central part of the array
!!  !!do ispinor=1,nspinor
!!  !!   !$omp do 
!!  !!   do i3=i3s,i3e
!!
!!  !!      !first part
!!  !!      do i2=1,i2s-1
!!  !!         do i1=1,n1i
!!  !!            psirx(i1,i2,i3,ispinor)=0.0_wp 
!!  !!            psiry(i1,i2,i3,ispinor)=0.0_wp 
!!  !!            psirz(i1,i2,i3,ispinor)=0.0_wp 
!!  !!         end do
!!  !!      end do
!!  !!      !central part
!!  !!      do i2=i2s,i2e
!!  !!         do i1=1,i1s-1
!!  !!            psirx(i1,i2,i3,ispinor)=0.0_wp 
!!  !!            psiry(i1,i2,i3,ispinor)=0.0_wp 
!!  !!            psirz(i1,i2,i3,ispinor)=0.0_wp 
!!  !!         end do
!!  !!         do i1=i1e+1,n1i
!!  !!            psirx(i1,i2,i3,ispinor)=0.0_wp 
!!  !!            psiry(i1,i2,i3,ispinor)=0.0_wp 
!!  !!            psirz(i1,i2,i3,ispinor)=0.0_wp 
!!  !!         end do
!!  !!      end do
!!  !!      !last part
!!  !!      do i2=i2e+1,n2i
!!  !!         do i1=1,n1i
!!  !!            psirx(i1,i2,i3,ispinor)=0.0_wp 
!!  !!            psiry(i1,i2,i3,ispinor)=0.0_wp 
!!  !!            psirz(i1,i2,i3,ispinor)=0.0_wp 
!!  !!         end do
!!  !!      end do
!!
!!  !!   end do
!!  !!   !$omp end do
!!  !!end do
!!
!!
!!  !!!last part of the array
!!  !!do ispinor=1,nspinor
!!  !!   !$omp do 
!!  !!   do i3=i3e+1,n3i
!!  !!      do i2=1,n2i
!!  !!         do i1=1,n1i
!!  !!            psirx(i1,i2,i3,ispinor)=0.0_wp 
!!  !!            psiry(i1,i2,i3,ispinor)=0.0_wp 
!!  !!            psirz(i1,i2,i3,ispinor)=0.0_wp 
!!  !!         end do
!!  !!      end do
!!  !!   end do
!!  !!   !$omp end do
!!  !!end do
!!
!!
!!  !important part of the array
!!  if (nspinor==4) then
!!      stop 'not yet implemented for nspinor==4!'
!!     !!!$omp do
!!     !!do i3=i3s,i3e
!!     !!   do i2=i2s,i2e
!!     !!      !thanks to the optional argument the conditional is done at compile time
!!     !!      if (present(ibyyzz_r)) then
!!     !!         i1st=max(i1s,ibyyzz_r(1,i2-15,i3-15)+1) !in bounds coordinates
!!     !!         i1et=min(i1e,ibyyzz_r(2,i2-15,i3-15)+1) !in bounds coordinates
!!     !!      else
!!     !!         i1st=i1s
!!     !!         i1et=i1e
!!     !!      end if
!!     !!      !no need of setting up to zero values outside wavefunction bounds
!!     !!      do i1=i1st,i1et
!!     !!         !wavefunctions
!!     !!         psir1=psir(i1,i2,i3,1)
!!     !!         psir2=psir(i1,i2,i3,2)
!!     !!         psir3=psir(i1,i2,i3,3)
!!     !!         psir4=psir(i1,i2,i3,4)
!!     !!         !potentials + confining term
!!     !!         pot1=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),1)+cp(i1,i2,i3)
!!     !!         pot2=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),2)+cp(i1,i2,i3)
!!     !!         pot3=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),3)+cp(i1,i2,i3)
!!     !!         pot4=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),4)+cp(i1,i2,i3)
!!
!!     !!         !diagonal terms
!!     !!         tt11=pot1*psir1 !p1
!!     !!         tt22=pot1*psir2 !p2
!!     !!         tt33=pot4*psir3 !p3
!!     !!         tt44=pot4*psir4 !p4
!!     !!         !Rab*Rb
!!     !!         tt13=pot2*psir3 !p1
!!     !!         !Iab*Ib
!!     !!         tt14=pot3*psir4 !p1
!!     !!         !Rab*Ib
!!     !!         tt23=pot2*psir4 !p2
!!     !!         !Iab*Rb
!!     !!         tt24=pot3*psir3 !p2
!!     !!         !Rab*Ra
!!     !!         tt31=pot2*psir1 !p3
!!     !!         !Iab*Ia
!!     !!         tt32=pot3*psir2 !p3
!!     !!         !Rab*Ia
!!     !!         tt41=pot2*psir2 !p4
!!     !!         !Iab*Ra
!!     !!         tt42=pot3*psir1 !p4
!!
!!     !!         !value of the potential energy
!!     !!         epot_p=epot_p+tt11*psir1+tt22*psir2+tt33*psir3+tt44*psir4+&
!!     !!              2.0_gp*tt31*psir3-2.0_gp*tt42*psir4+2.0_gp*tt41*psir4+2.0_gp*tt32*psir3
!!
!!     !!         !wavefunction update
!!     !!         !p1=h1p1+h2p3-h3p4
!!     !!         !p2=h1p2+h2p4+h3p3
!!     !!         !p3=h2p1+h3p2+h4p3
!!     !!         !p4=h2p2-h3p1+h4p4
!!     !!         psir(i1,i2,i3,1)=tt11+tt13-tt14
!!     !!         psir(i1,i2,i3,2)=tt22+tt23+tt24
!!     !!         psir(i1,i2,i3,3)=tt33+tt31+tt32
!!     !!         psir(i1,i2,i3,4)=tt44+tt41-tt42
!!     !!      end do
!!     !!   end do
!!     !!end do
!!     !!!$omp end do
!!
!!  else !case with nspinor /=4
!!     do ispinor=1,nspinor
!!        !$omp do
!!        do i3=i3s,i3e
!!           do i2=i2s,i2e
!!              !thanks to the optional argument the conditional is done at compile time
!!              if (present(ibyyzz_r)) then
!!                 i1st=max(i1s,ibyyzz_r(1,i2-15,i3-15)+1) !in bounds coordinates
!!                 i1et=min(i1e,ibyyzz_r(2,i2-15,i3-15)+1) !in bounds coordinates
!!              else
!!                 i1st=i1s
!!                 i1et=i1e
!!              end if
!!              !no need of setting up to zero values outside wavefunction bounds
!!              do i1=i1st,i1et
!!                 psir1=psir(i1,i2,i3,ispinor)
!!                 !the local potential is always real (npot=1) + confining term
!!                 !!pot1=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),1)+cp(i1,i2,i3)
!!                 ttx=(confdata%hh(1)*real(i1+confdata%ioffset(1),wp)-confdata%rxyzConf(1))**2
!!                 tty=(confdata%hh(2)*real(i2+confdata%ioffset(2),wp)-confdata%rxyzConf(2))**2
!!                 ttz=(confdata%hh(3)*real(i3+confdata%ioffset(3),wp)-confdata%rxyzConf(3))**2
!!                 tt=ttx+tty+ttz
!!                 !write(1000,*) tt, cutoff2, psir1**2
!!                 if(tt>cutoff2) then
!!                     weight_out=weight_out+psir1**2
!!                 else
!!                     weight_in=weight_in+psir1**2
!!                 end if
!!              end do
!!           end do
!!        end do
!!        !$omp end do
!!     end do
!!  end if
!!  
!!  
!!  !$omp end parallel
!!
!!
!!END SUBROUTINE get_cutoff_weight
