!> @file
!!  ensemble of routines needed to apply the local potential to the real-space wavefunction
!! @author
!!   Copyright (C) 2005-2013 BigDFT group 
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 


!> Routine for applying the local potential
!! Support the adding of a confining potential and the localisation region of the potential
subroutine apply_potential_lr_conf_noconf(n1i,n2i,n3i,n1ip,n2ip,n3ip,ishift,n2,n3,&
     nspinor,npot,psir,pot,epot,confdata,ibyyzz_r,psir_noconf,econf)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,nspinor,npot
  integer, dimension(3), intent(in) :: ishift !<offset of potential box in wfn box coords.
  real(wp), dimension(n1i,n2i,n3i,nspinor), intent(inout) :: psir !< real-space wfn in lr
  real(wp), dimension(n1ip,n2ip,n3ip,npot), intent(in) :: pot !< real-space pot in lrb
  type(confpot_data), intent(in) :: confdata !< data for the confining potential
  integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in) :: ibyyzz_r !< bounds in lr
  real(gp), intent(out) :: epot
  real(wp),dimension(n1i,n2i,n3i,nspinor),intent(inout) :: psir_noconf !< real-space wfn in lr where only the potential (without confinement) will be applied
  real(gp), intent(out) :: econf
  !local variables
  integer :: i1,i2,i3,ispinor,i1s,i1e,i2s,i2e,i3s,i3e,i1st,i1et,ii1,ii2,ii3,potorder_half, m
  real(wp) :: x0, x1, x2, x3, r02, r12, r22, r32, cp0, cp1, cp2, cp3, psir01, psir11, psir21, psir31
  real(wp) :: pot01, pot11, pot21, pot31, tt011, tt111, tt211, tt311 
  real(wp) :: pot01_noconf, pot11_noconf, pot21_noconf, pot31_noconf
  real(wp) :: tt011_noconf, tt111_noconf, tt211_noconf, tt311_noconf
  integer :: ii01, ii11, ii21, ii31
  real(wp) :: tt11,tt11_noconf,cp,r2,z2,y2
  real(wp) :: psir1,pot1,pot1_noconf,x,y,z
  real(gp) :: epot_p,econf_p!,ierr

  if (f_err_raise(nspinor==4,'nspinor=4 not supported with noconf')) return

  epot=0.0_wp
  econf=0.d0

  !loop on wavefunction
  !calculate the limits in all the directions
  !regions in which both the potential and wavefunctions are defined
  i3s=max(1,ishift(3)+1)
  i3e=min(n3i,n3ip+ishift(3))
  i2s=max(1,ishift(2)+1)
  i2e=min(n2i,n2ip+ishift(2))
  i1s=max(1,ishift(1)+1)
  i1e=min(n1i,n1ip+ishift(1))

  potorder_half=confdata%potorder/2

  !$omp parallel default(private)&
  !$omp shared(pot,psir,n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,epot,ibyyzz_r,nspinor)&
  !$omp shared(i1s,i1e,i2s,i2e,i3s,i3e,ishift,psir_noconf,econf,confdata,potorder_half)&
  !$omp private(ispinor,i1,i2,i3,epot_p,i1st,i1et,pot1_noconf,tt11_noconf,econf_p)&
  !$omp private(tt11,psir1,pot1,ii1,ii2,ii3)

  !case without bounds
  epot_p=0._gp
  econf_p=0._gp

  !put to zero the external part of psir if the potential is more little than the wavefunction
  !first part of the array
  do ispinor=1,nspinor
     !$omp do 
     do i3=1,i3s-1
        do i2=1,n2i
           do i1=1,n1i
             psir(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
     end do
     !$omp end do
  end do

  !central part of the array
  do ispinor=1,nspinor
     !$omp do 
     do i3=i3s,i3e

        !first part
        do i2=1,i2s-1
           do i1=1,n1i
              psir(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
        !central part
        do i2=i2s,i2e
           do i1=1,i1s-1
              psir(i1,i2,i3,ispinor)=0.0_wp 
           end do
           do i1=i1e+1,n1i
              psir(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
        !last part
        do i2=i2e+1,n2i
           do i1=1,n1i
              psir(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do

     end do
     !$omp end do
  end do

  !last part of the array
  do ispinor=1,nspinor
     !$omp do 
     do i3=i3e+1,n3i
        do i2=1,n2i
           do i1=1,n1i
              psir(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
     end do
     !$omp end do
  end do

  !important part of the array
  do ispinor=1,nspinor
     !$omp do reduction(+:epot,econf)
     do i3=i3s,i3e
        z=confdata%hh(3)*real(i3+confdata%ioffset(3),wp)-confdata%rxyzConf(3)
        z2=z**2
        ii3=i3-ishift(3)
        do i2=i2s,i2e
           y=confdata%hh(2)*real(i2+confdata%ioffset(2),wp)-confdata%rxyzConf(2)
           y2=y**2
           ii2=i2-ishift(2)

           i1st=max(i1s,ibyyzz_r(1,i2-15,i3-15)+1) !in bounds coordinates
           i1et=min(i1e,ibyyzz_r(2,i2-15,i3-15)+1) !in bounds coordinates

           !no need of setting up to zero values outside wavefunction bounds
           m=mod(i1et-i1st+1,4)
           !m=i1et-i1st+1
           !m=0
           if (m/=0) then
               do i1=i1st,i1st+m-1
                  x=confdata%hh(1)*real(i1+confdata%ioffset(1),wp)-confdata%rxyzConf(1)
                  r2=x**2+y2+z2
                  cp=confdata%prefac*r2**potorder_half

                  psir1=psir(i1,i2,i3,ispinor)
                  !the local potential is always real (npot=1) + confining term
                  ii1=i1-ishift(1)
                  pot1=pot(ii1,ii2,ii3,1)+cp
                  tt11=pot1*psir1

                  pot1_noconf=pot(ii1,ii2,ii3,1)
                  tt11_noconf=pot1_noconf*psir1
                  psir_noconf(i1,i2,i3,ispinor) = tt11_noconf
                  econf=econf+real(cp*psir1*psir1,wp)

                  epot=epot+real(tt11*psir1,wp)
                  psir(i1,i2,i3,ispinor)=tt11
               end do
           end if
           do i1=i1st+m,i1et,4
              x0=confdata%hh(1)*real(i1+0+confdata%ioffset(1),wp)-confdata%rxyzConf(1)
              x1=confdata%hh(1)*real(i1+1+confdata%ioffset(1),wp)-confdata%rxyzConf(1)
              x2=confdata%hh(1)*real(i1+2+confdata%ioffset(1),wp)-confdata%rxyzConf(1)
              x3=confdata%hh(1)*real(i1+3+confdata%ioffset(1),wp)-confdata%rxyzConf(1)
              r02=x0**2+y2+z2
              r12=x1**2+y2+z2
              r22=x2**2+y2+z2
              r32=x3**2+y2+z2
              cp0=confdata%prefac*r02**potorder_half
              cp1=confdata%prefac*r12**potorder_half
              cp2=confdata%prefac*r22**potorder_half
              cp3=confdata%prefac*r32**potorder_half

              psir01=psir(i1+0,i2,i3,ispinor)
              psir11=psir(i1+1,i2,i3,ispinor)
              psir21=psir(i1+2,i2,i3,ispinor)
              psir31=psir(i1+3,i2,i3,ispinor)
              !the local potential is always real (npot=1) + confining term
              ii01=i1+0-ishift(1)
              ii11=i1+1-ishift(1)
              ii21=i1+2-ishift(1)
              ii31=i1+3-ishift(1)
              pot01=pot(ii01,ii2,ii3,1)+cp0
              pot11=pot(ii11,ii2,ii3,1)+cp1
              pot21=pot(ii21,ii2,ii3,1)+cp2
              pot31=pot(ii31,ii2,ii3,1)+cp3
              tt011=pot01*psir01
              tt111=pot11*psir11
              tt211=pot21*psir21
              tt311=pot31*psir31

              pot01_noconf=pot(ii01,ii2,ii3,1)
              pot11_noconf=pot(ii11,ii2,ii3,1)
              pot21_noconf=pot(ii21,ii2,ii3,1)
              pot31_noconf=pot(ii31,ii2,ii3,1)
              tt011_noconf=pot01_noconf*psir01
              tt111_noconf=pot11_noconf*psir11
              tt211_noconf=pot21_noconf*psir21
              tt311_noconf=pot31_noconf*psir31

              psir_noconf(i1+0,i2,i3,ispinor) = tt011_noconf
              psir_noconf(i1+1,i2,i3,ispinor) = tt111_noconf
              psir_noconf(i1+2,i2,i3,ispinor) = tt211_noconf
              psir_noconf(i1+3,i2,i3,ispinor) = tt311_noconf

              econf=econf+real(cp0*psir01*psir01,wp)
              econf=econf+real(cp1*psir11*psir11,wp)
              econf=econf+real(cp2*psir21*psir21,wp)
              econf=econf+real(cp3*psir31*psir31,wp)

              epot=epot+real(tt011*psir01,wp)
              epot=epot+real(tt111*psir11,wp)
              epot=epot+real(tt211*psir21,wp)
              epot=epot+real(tt311*psir31,wp)
              psir(i1+0,i2,i3,ispinor)=tt011
              psir(i1+1,i2,i3,ispinor)=tt111
              psir(i1+2,i2,i3,ispinor)=tt211
              psir(i1+3,i2,i3,ispinor)=tt311
           end do
           !!do i1=i1st,i1et
           !!   x=confdata%hh(1)*real(i1+confdata%ioffset(1),wp)-confdata%rxyzConf(1)
           !!   r2=x**2+y2+z2
           !!   cp=confdata%prefac*r2**potorder_half

           !!   psir1=psir(i1,i2,i3,ispinor)
           !!   !the local potential is always real (npot=1) + confining term
           !!   ii1=i1-ishift(1)
           !!   pot1=pot(ii1,ii2,ii3,1)+cp
           !!   tt11=pot1*psir1

           !!   pot1_noconf=pot(ii1,ii2,ii3,1)
           !!   tt11_noconf=pot1_noconf*psir1
           !!   psir_noconf(i1,i2,i3,ispinor) = tt11_noconf
           !!   econf=econf+real(cp*psir1*psir1,wp)

           !!   epot=epot+real(tt11*psir1,wp)
           !!   psir(i1,i2,i3,ispinor)=tt11
           !!end do
        end do
     end do
     !$omp end do
  end do
  
  !!!$omp critical
  !!epot=epot+epot_p
  !!
  !!econf=econf+econf_p
  !!!$omp end critical
  
  !$omp end parallel

END SUBROUTINE apply_potential_lr_conf_noconf

!> routine for applying the local potential
!! Support the adding of a confining potential and the localisation region of the potential
subroutine apply_potential_lr_conf(n1i,n2i,n3i,n1ip,n2ip,n3ip,ishift,n2,n3,&
     nspinor,npot,psir,pot,epot,confdata,ibyyzz_r)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,nspinor,npot
  integer, dimension(3), intent(in) :: ishift !<offset of potential box in wfn box coords.
  real(wp), dimension(n1i,n2i,n3i,nspinor), intent(inout) :: psir !< real-space wfn in lr
  real(wp), dimension(n1ip,n2ip,n3ip,npot), intent(in) :: pot !< real-space pot in lrb
  type(confpot_data), intent(in) :: confdata !< data for the confining potential
  integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in) :: ibyyzz_r !< bounds in lr
  real(gp), intent(out) :: epot
  !local variables
  integer :: i1,i2,i3,ispinor,i1s,i1e,i2s,i2e,i3s,i3e,i1st,i1et,ii1,ii2,ii3,potorder_half
  real(wp) :: tt11,cp,r2,z2,y2
  real(wp) :: psir1,pot1,x,y,z
  real(gp) :: epot_p!,ierr

  if (f_err_raise(nspinor==4,&
       'nspinor=4 not supported with confining potential')) return

  epot=0.0_wp

  !loop on wavefunction
  !calculate the limits in all the directions
  !regions in which both the potential and wavefunctions are defined
  i3s=max(1,ishift(3)+1)
  i3e=min(n3i,n3ip+ishift(3))
  i2s=max(1,ishift(2)+1)
  i2e=min(n2i,n2ip+ishift(2))
  i1s=max(1,ishift(1)+1)
  i1e=min(n1i,n1ip+ishift(1))

  potorder_half=confdata%potorder/2

  !$omp parallel default(private)&
  !$omp shared(pot,psir,n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,epot,ibyyzz_r,nspinor)&
  !$omp shared(i1s,i1e,i2s,i2e,i3s,i3e,ishift,confdata,potorder_half)&
  !$omp private(ispinor,i1,i2,i3,epot_p,i1st,i1et)&
  !$omp private(tt11,psir1,pot1,ii1,ii2,ii3)

  !case without bounds
  epot_p=0._gp

  !put to zero the external part of psir if the potential is more little than the wavefunction
  !first part of the array
  do ispinor=1,nspinor
     !$omp do 
     do i3=1,i3s-1
        do i2=1,n2i
           do i1=1,n1i
             psir(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
     end do
     !$omp end do
  end do

  !central part of the array
  do ispinor=1,nspinor
     !$omp do 
     do i3=i3s,i3e

        !first part
        do i2=1,i2s-1
           do i1=1,n1i
              psir(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
        !central part
        do i2=i2s,i2e
           do i1=1,i1s-1
              psir(i1,i2,i3,ispinor)=0.0_wp 
           end do
           do i1=i1e+1,n1i
              psir(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
        !last part
        do i2=i2e+1,n2i
           do i1=1,n1i
              psir(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do

     end do
     !$omp end do
  end do

  !last part of the array
  do ispinor=1,nspinor
     !$omp do 
     do i3=i3e+1,n3i
        do i2=1,n2i
           do i1=1,n1i
              psir(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
     end do
     !$omp end do
  end do

  !important part of the array
  do ispinor=1,nspinor
     !$omp do reduction(+:epot)
     do i3=i3s,i3e
        z=confdata%hh(3)*real(i3+confdata%ioffset(3),wp)-confdata%rxyzConf(3)
        z2=z**2
        ii3=i3-ishift(3)
        do i2=i2s,i2e
           y=confdata%hh(2)*real(i2+confdata%ioffset(2),wp)-confdata%rxyzConf(2)
           y2=y**2
           ii2=i2-ishift(2)

           i1st=max(i1s,ibyyzz_r(1,i2-15,i3-15)+1) !in bounds coordinates
           i1et=min(i1e,ibyyzz_r(2,i2-15,i3-15)+1) !in bounds coordinates

           !no need of setting up to zero values outside wavefunction bounds
           do i1=i1st,i1et
              x=confdata%hh(1)*real(i1+confdata%ioffset(1),wp)-confdata%rxyzConf(1)
              r2=x**2+y2+z2
              cp=confdata%prefac*r2**potorder_half

              psir1=psir(i1,i2,i3,ispinor)
              !the local potential is always real (npot=1) + confining term
              ii1=i1-ishift(1)
              pot1=pot(ii1,ii2,ii3,1)+cp
              tt11=pot1*psir1

              epot=epot+real(tt11*psir1,wp)
              psir(i1,i2,i3,ispinor)=tt11
           end do
        end do
     end do
     !$omp end do
  end do
  
  !!!$omp critical
  !!epot=epot+epot_p
  !!!$omp end critical
  
  !$omp end parallel

END SUBROUTINE apply_potential_lr_conf


!> Routine for applying the local potential
!! Support the adding of a confining potential and the localisation region of the potential
subroutine apply_potential_lr_conf_nobounds(n1i,n2i,n3i,n1ip,n2ip,n3ip,ishift,n2,n3,nspinor,npot,psir,pot,epot,confdata)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,nspinor,npot
  integer, dimension(3), intent(in) :: ishift !<offset of potential box in wfn box coords.
  real(wp), dimension(n1i,n2i,n3i,nspinor), intent(inout) :: psir !< real-space wfn in lr
  real(wp), dimension(n1ip,n2ip,n3ip,npot), intent(in) :: pot !< real-space pot in lrb
  type(confpot_data), intent(in) :: confdata !< data for the confining potential
  real(gp), intent(out) :: epot
  !local variables
  integer :: i1,i2,i3,ispinor,i1s,i1e,i2s,i2e,i3s,i3e,ii1,ii2,ii3,potorder_half
  real(wp) :: tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42,r2
  real(wp) :: psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4,x,y,z,cp,z2,y2
  real(gp) :: epot_p!,ierr

  epot=0.0_wp

  !loop on wavefunction
  !calculate the limits in all the directions
  !regions in which both the potential and wavefunctions are defined
  i3s=max(1,ishift(3)+1)
  i3e=min(n3i,n3ip+ishift(3))
  i2s=max(1,ishift(2)+1)
  i2e=min(n2i,n2ip+ishift(2))
  i1s=max(1,ishift(1)+1)
  i1e=min(n1i,n1ip+ishift(1))

  potorder_half=confdata%potorder/2

  !$omp parallel default(private)&
  !$omp shared(pot,psir,n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,epot,nspinor)&
  !$omp shared(i1s,i1e,i2s,i2e,i3s,i3e,ishift,confdata,potorder_half)&
  !$omp private(ispinor,i1,i2,i3,epot_p)&
  !$omp private(tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42)&
  !$omp private(psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4,ii1,ii2,ii3)

  !case without bounds

  epot_p=0._gp

  !put to zero the external part of psir if the potential is more little than the wavefunction
  !first part of the array
  do ispinor=1,nspinor
     !$omp do 
     do i3=1,i3s-1
        do i2=1,n2i
           do i1=1,n1i
             psir(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
     end do
     !$omp end do
  end do

  !central part of the array
  do ispinor=1,nspinor
     !$omp do 
     do i3=i3s,i3e

        !first part
        do i2=1,i2s-1
           do i1=1,n1i
              psir(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
        !central part
        do i2=i2s,i2e
           do i1=1,i1s-1
              psir(i1,i2,i3,ispinor)=0.0_wp 
           end do
           do i1=i1e+1,n1i
              psir(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
        !last part
        do i2=i2e+1,n2i
           do i1=1,n1i
              psir(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
     end do
     !$omp end do
  end do

  !last part of the array
  do ispinor=1,nspinor
     !$omp do 
     do i3=i3e+1,n3i
        do i2=1,n2i
           do i1=1,n1i
              psir(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
     end do
     !$omp end do
  end do

  !important part of the array
  if (nspinor==4) then
     !$omp do reduction(+:epot)
     do i3=i3s,i3e
        z=confdata%hh(3)*real(i3+confdata%ioffset(3),wp)-confdata%rxyzConf(3)
        z2=z**2
        do i2=i2s,i2e
           y=confdata%hh(2)*real(i2+confdata%ioffset(2),wp)-confdata%rxyzConf(2)
           y2=y**2
           !thanks to the optional argument the conditional is done at compile time

           !no need of setting up to zero values outside wavefunction bounds
           do i1=i1s,i1e
              x=confdata%hh(1)*real(i1+confdata%ioffset(1),wp)-confdata%rxyzConf(1)
              r2=x**2+y2+z2
              cp=confdata%prefac*r2**potorder_half

              !wavefunctions
              psir1=psir(i1,i2,i3,1)
              psir2=psir(i1,i2,i3,2)
              psir3=psir(i1,i2,i3,3)
              psir4=psir(i1,i2,i3,4)
              !potentials + confining term
              pot1=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),1)+cp
              pot2=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),2)+cp
              pot3=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),3)+cp
              pot4=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),4)+cp

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

              !value of the potential energy
              epot=epot+tt11*psir1+tt22*psir2+tt33*psir3+tt44*psir4+&
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
        end do
     end do
     !$omp end do

  else !case with nspinor /=4
     do ispinor=1,nspinor
        !$omp do reduction(+:epot)
        do i3=i3s,i3e
           z=confdata%hh(3)*real(i3+confdata%ioffset(3),wp)-confdata%rxyzConf(3)
           z2=z**2
           ii3=i3-ishift(3)
           do i2=i2s,i2e
              y=confdata%hh(2)*real(i2+confdata%ioffset(2),wp)-confdata%rxyzConf(2)
              y2=y**2
              ii2=i2-ishift(2)
              !thanks to the optional argument the conditional is done at compile time
              !no need of setting up to zero values outside wavefunction bounds
              do i1=i1s,i1e
                 x=confdata%hh(1)*real(i1+confdata%ioffset(1),wp)-confdata%rxyzConf(1)
                 r2=x**2+y2+z2
                 cp=confdata%prefac*r2**potorder_half

                 psir1=psir(i1,i2,i3,ispinor)
                 !the local potential is always real (npot=1) + confining term
                 ii1=i1-ishift(1)
                 pot1=pot(ii1,ii2,ii3,1)+cp
                 tt11=pot1*psir1

                 epot=epot+real(tt11*psir1,wp)
                 psir(i1,i2,i3,ispinor)=tt11
              end do
           end do
        end do
        !$omp end do
     end do
  end if
  
  !!!$omp critical
  !!epot=epot+epot_p
  !!!$omp end critical
  
  !$omp end parallel

END SUBROUTINE apply_potential_lr_conf_nobounds

!>   routine for applying the local potential
!! Support the adding of a confining potential and the localisation region of the potential
subroutine apply_potential_lr_nobounds(n1i,n2i,n3i,n1ip,n2ip,n3ip,ishift,n2,n3,nspinor,npot,psir,pot,epot)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,nspinor,npot
  integer, dimension(3), intent(in) :: ishift !<offset of potential box in wfn box coords.
  real(wp), dimension(n1i,n2i,n3i,nspinor), intent(inout) :: psir !< real-space wfn in lr
  real(wp), dimension(n1ip,n2ip,n3ip,npot), intent(in) :: pot !< real-space pot in lrb
  real(gp), intent(out) :: epot
  !local variables
  integer :: i1,i2,i3,ispinor,i1s,i1e,i2s,i2e,i3s,i3e,ii1,ii2,ii3
  real(wp) :: tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42
  real(wp) :: psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4
  real(gp) :: epot_p!,ierr

  epot=0.0_wp

  !loop on wavefunction
  !calculate the limits in all the directions
  !regions in which both the potential and wavefunctions are defined
  i3s=max(1,ishift(3)+1)
  i3e=min(n3i,n3ip+ishift(3))
  i2s=max(1,ishift(2)+1)
  i2e=min(n2i,n2ip+ishift(2))
  i1s=max(1,ishift(1)+1)
  i1e=min(n1i,n1ip+ishift(1))

  !$omp parallel default(private)&
  !$omp shared(pot,psir,n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,epot,nspinor)&
  !$omp shared(i1s,i1e,i2s,i2e,i3s,i3e,ishift)&
  !$omp private(ispinor,i1,i2,i3,epot_p)&
  !$omp private(tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42)&
  !$omp private(psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4,ii1,ii2,ii3)

  !case without bounds

  epot_p=0._gp

  !put to zero the external part of psir if the potential is more little than the wavefunction
  !first part of the array
  do ispinor=1,nspinor
     !$omp do 
     do i3=1,i3s-1
        do i2=1,n2i
           do i1=1,n1i
             psir(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
     end do
     !$omp end do
  end do

  !central part of the array
  do ispinor=1,nspinor
     !$omp do 
     do i3=i3s,i3e

        !first part
        do i2=1,i2s-1
           do i1=1,n1i
              psir(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
        !central part
        do i2=i2s,i2e
           do i1=1,i1s-1
              psir(i1,i2,i3,ispinor)=0.0_wp 
           end do
           do i1=i1e+1,n1i
              psir(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
        !last part
        do i2=i2e+1,n2i
           do i1=1,n1i
              psir(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
     end do
     !$omp end do
  end do

  !last part of the array
  do ispinor=1,nspinor
     !$omp do 
     do i3=i3e+1,n3i
        do i2=1,n2i
           do i1=1,n1i
              psir(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
     end do
     !$omp end do
  end do

  !important part of the array
  if (nspinor==4) then
     !$omp do reduction(+:epot)
     do i3=i3s,i3e
        do i2=i2s,i2e
           !thanks to the optional argument the conditional is done at compile time

           !no need of setting up to zero values outside wavefunction bounds
           do i1=i1s,i1e

              !wavefunctions
              psir1=psir(i1,i2,i3,1)
              psir2=psir(i1,i2,i3,2)
              psir3=psir(i1,i2,i3,3)
              psir4=psir(i1,i2,i3,4)
              !potentials + confining term
              pot1=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),1)
              pot2=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),2)
              pot3=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),3)
              pot4=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),4)

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

              !value of the potential energy
              epot=epot+tt11*psir1+tt22*psir2+tt33*psir3+tt44*psir4+&
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
        end do
     end do
     !$omp end do

  else !case with nspinor /=4
     do ispinor=1,nspinor
        !$omp do reduction(+:epot)
        do i3=i3s,i3e
           ii3=i3-ishift(3)
           do i2=i2s,i2e
              ii2=i2-ishift(2)
              !thanks to the optional argument the conditional is done at compile time
              !no need of setting up to zero values outside wavefunction bounds
              do i1=i1s,i1e

                 psir1=psir(i1,i2,i3,ispinor)
                 !the local potential is always real (npot=1) + confining term
                 ii1=i1-ishift(1)
                 pot1=pot(ii1,ii2,ii3,1)
                 tt11=pot1*psir1

                 epot=epot+real(tt11*psir1,wp)
                 psir(i1,i2,i3,ispinor)=tt11
              end do
           end do
        end do
        !$omp end do
     end do
  end if
  
  !!!$omp critical
  !!epot=epot+epot_p
  !!!$omp end critical
  
  !$omp end parallel

END SUBROUTINE apply_potential_lr_nobounds

!>   routine for applying the local potential
!! Support the adding of a confining potential and the localisation region of the potential
subroutine apply_potential_lr_bounds(n1i,n2i,n3i,n1ip,n2ip,n3ip,ishift,n2,n3,nspinor,npot,psir,pot,epot,ibyyzz_r)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,nspinor,npot
  integer, dimension(3), intent(in) :: ishift !<offset of potential box in wfn box coords.
  real(wp), dimension(n1i,n2i,n3i,nspinor), intent(inout) :: psir !< real-space wfn in lr
  real(wp), dimension(n1ip,n2ip,n3ip,npot), intent(in) :: pot !< real-space pot in lrb
  integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in) :: ibyyzz_r !< bounds in lr
  real(gp), intent(out) :: epot
  !local variables
  integer :: i1,i2,i3,ispinor,i1s,i1e,i2s,i2e,i3s,i3e,ii1,ii2,ii3,i1st,i1et
  integer :: ii01, ii11, ii21, ii31, m
  real(wp) :: psir01, psir11, psir21, psir31, pot01, pot11, pot21, pot31, tt011, tt111, tt211, tt311
  real(wp) :: tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42
  real(wp) :: psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4
  real(gp) :: epot_p!,ierr

  epot=0.0_wp

  !loop on wavefunction
  !calculate the limits in all the directions
  !regions in which both the potential and wavefunctions are defined
  i3s=max(1,ishift(3)+1)
  i3e=min(n3i,n3ip+ishift(3))
  i2s=max(1,ishift(2)+1)
  i2e=min(n2i,n2ip+ishift(2))
  i1s=max(1,ishift(1)+1)
  i1e=min(n1i,n1ip+ishift(1))

  !$omp parallel default(private)&
  !$omp shared(pot,psir,n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,epot,ibyyzz_r,nspinor)&
  !$omp shared(i1s,i1e,i2s,i2e,i3s,i3e,ishift)&
  !$omp private(ispinor,i1,i2,i3,epot_p)&
  !$omp private(tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42)&
  !$omp private(psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4,ii1,ii2,ii3)

  !case without bounds

  epot_p=0._gp

  !put to zero the external part of psir if the potential is more little than the wavefunction
  !first part of the array
  do ispinor=1,nspinor
     !$omp do 
     do i3=1,i3s-1
        do i2=1,n2i
           do i1=1,n1i
             psir(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
     end do
     !$omp end do
  end do

  !central part of the array
  do ispinor=1,nspinor
     !$omp do 
     do i3=i3s,i3e

        !first part
        do i2=1,i2s-1
           do i1=1,n1i
              psir(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
        !central part
        do i2=i2s,i2e
           do i1=1,i1s-1
              psir(i1,i2,i3,ispinor)=0.0_wp 
           end do
           do i1=i1e+1,n1i
              psir(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
        !last part
        do i2=i2e+1,n2i
           do i1=1,n1i
              psir(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
     end do
     !$omp end do
  end do

  !last part of the array
  do ispinor=1,nspinor
     !$omp do 
     do i3=i3e+1,n3i
        do i2=1,n2i
           do i1=1,n1i
              psir(i1,i2,i3,ispinor)=0.0_wp 
           end do
        end do
     end do
     !$omp end do
  end do

  !important part of the array
  if (nspinor==4) then
     !$omp do reduction(+:epot)
     do i3=i3s,i3e
        do i2=i2s,i2e
           !thanks to the optional argument the conditional is done at compile time
           i1st=max(i1s,ibyyzz_r(1,i2-15,i3-15)+1) !in bounds coordinates
           i1et=min(i1e,ibyyzz_r(2,i2-15,i3-15)+1) !in bounds coordinates

           !no need of setting up to zero values outside wavefunction bounds
           do i1=i1st,i1et

              !wavefunctions
              psir1=psir(i1,i2,i3,1)
              psir2=psir(i1,i2,i3,2)
              psir3=psir(i1,i2,i3,3)
              psir4=psir(i1,i2,i3,4)
              !potentials + confining term
              pot1=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),1)
              pot2=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),2)
              pot3=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),3)
              pot4=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),4)

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

              !value of the potential energy
              epot=epot+tt11*psir1+tt22*psir2+tt33*psir3+tt44*psir4+&
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
        end do
     end do
     !$omp end do

  else !case with nspinor /=4
     do ispinor=1,nspinor
        !$omp do reduction(+:epot)
        do i3=i3s,i3e
           ii3=i3-ishift(3)
           do i2=i2s,i2e
              ii2=i2-ishift(2)
              i1st=max(i1s,ibyyzz_r(1,i2-15,i3-15)+1) !in bounds coordinates
              i1et=min(i1e,ibyyzz_r(2,i2-15,i3-15)+1) !in bounds coordinates
              !no need of setting up to zero values outside wavefunction bounds
              m=mod(i1et-i1st+1,4)
              if (m/=0) then
                  do i1=i1st,i1et

                     psir1=psir(i1,i2,i3,ispinor)
                     !the local potential is always real (npot=1) + confining term
                     ii1=i1-ishift(1)
                     pot1=pot(ii1,ii2,ii3,1)
                     tt11=pot1*psir1

                     epot_p=epot_p+real(tt11*psir1,wp)
                     psir(i1,i2,i3,ispinor)=tt11
                  end do
              end if
              do i1=i1st,i1st+m-1

                 psir01=psir(i1+0,i2,i3,ispinor)
                 psir11=psir(i1+1,i2,i3,ispinor)
                 psir21=psir(i1+2,i2,i3,ispinor)
                 psir31=psir(i1+3,i2,i3,ispinor)
                 !the local potential is always real (npot=1) + confining term
                 ii01=i1+0-ishift(1)
                 ii11=i1+1-ishift(1)
                 ii21=i1+2-ishift(1)
                 ii31=i1+3-ishift(1)
                 pot01=pot(ii01,ii2,ii3,1)
                 pot11=pot(ii11,ii2,ii3,1)
                 pot21=pot(ii21,ii2,ii3,1)
                 pot31=pot(ii31,ii2,ii3,1)
                 tt011=pot01*psir01
                 tt111=pot11*psir11
                 tt211=pot21*psir21
                 tt311=pot31*psir31

                 epot_p=epot_p+real(tt011*psir01,wp)
                 epot_p=epot_p+real(tt111*psir11,wp)
                 epot_p=epot_p+real(tt211*psir21,wp)
                 epot_p=epot_p+real(tt311*psir31,wp)
                 psir(i1+0,i2,i3,ispinor)=tt011
                 psir(i1+1,i2,i3,ispinor)=tt111
                 psir(i1+2,i2,i3,ispinor)=tt211
                 psir(i1+3,i2,i3,ispinor)=tt311
              end do
           end do
        end do
        !$omp end do
     end do
  end if
  
  !!!$omp critical
  !!epot=epot+epot_p
  !!!$omp end critical
  
  !$omp end parallel

END SUBROUTINE apply_potential_lr_bounds
