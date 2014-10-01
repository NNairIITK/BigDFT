!>   routine for applying the local potential
!! Support the adding of a confining potential and the localisation region of the potential
subroutine apply_potential_lr(n1i,n2i,n3i,n1ip,n2ip,n3ip,ishift,n2,n3,nspinor,npot,psir,pot,epot,&
     confdata,ibyyzz_r,psir_noconf,econf) !optional
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,nspinor,npot
  integer, dimension(3), intent(in) :: ishift !<offset of potential box in wfn box coords.
  real(wp), dimension(n1i,n2i,n3i,nspinor), intent(inout) :: psir !< real-space wfn in lr
  real(wp), dimension(n1ip,n2ip,n3ip,npot), intent(in) :: pot !< real-space pot in lrb
  type(confpot_data), intent(in), optional, target :: confdata !< data for the confining potential
  integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r !< bounds in lr
  real(gp), intent(out) :: epot
  real(wp),dimension(n1i,n2i,n3i,nspinor),intent(inout),optional :: psir_noconf !< real-space wfn in lr where only the potential (without confinement) will be applied
  real(gp), intent(out),optional :: econf
  !local variables
  integer :: i1,i2,i3,ispinor,i1s,i1e,i2s,i2e,i3s,i3e,i1st,i1et,ii1,ii2,ii3
  real(wp) :: tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42,tt11_noconf,ttt!,r2
  real(wp) :: psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4,pot1_noconf
  real(gp) :: epot_p,econf_p!,ierr
  !real(gp), dimension(3) :: hh,rxyzconf
  !integer, dimension(3) :: ioffset
  real(kind=8),dimension(:),pointer :: hh, rxyzConf
  integer,dimension(:),pointer :: ioffset
  real(kind=8),pointer :: prefac
  integer,pointer :: potorder
  real(kind=8),target :: zero_dble
  integer,target :: zero_int


  zero_dble=0.d0
  zero_int=0
  if (present(confdata)) then
      hh => confdata%hh
      ioffset => confdata%ioffset
      rxyzConf => confdata%rxyzConf
      prefac => confdata%prefac
      potorder => confdata%potorder
  else
     !these arrays are not deallocated!!!
      allocate(hh(3), rxyzConf(3), ioffset(3))
      hh=0.d0
      ioffset=0
      rxyzConf=0.d0
      prefac => zero_dble
      potorder => zero_int
  end if
  epot=0.0_wp

  if (present(econf)) then
      econf=0.d0
  end if

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
  !$omp shared(i1s,i1e,i2s,i2e,i3s,i3e,ishift,psir_noconf,econf)&
  !$omp shared(hh,ioffset,rxyzConf,prefac,potorder)&
  !$omp private(ispinor,i1,i2,i3,epot_p,i1st,i1et,pot1_noconf,tt11_noconf,econf_p)&
  !$omp private(tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42)&
  !$omp private(psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4,ii1,ii2,ii3,ttt)


!!$  !$omp parallel default(private)&
!!$  !$omp shared(pot,psir,n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,epot,ibyyzz_r,nspinor)&
!!$  !$omp shared(i1s,i1e,i2s,i2e,i3s,i3e,ishift)
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
  if (nspinor==4) then
     !$omp do
     do i3=i3s,i3e
        do i2=i2s,i2e
           !thanks to the optional argument the conditional is done at compile time
           if (present(ibyyzz_r)) then
              i1st=max(i1s,ibyyzz_r(1,i2-15,i3-15)+1) !in bounds coordinates
              i1et=min(i1e,ibyyzz_r(2,i2-15,i3-15)+1) !in bounds coordinates
           else
              i1st=i1s
              i1et=i1e
           end if
           !no need of setting up to zero values outside wavefunction bounds
           do i1=i1st,i1et
              !wavefunctions
              psir1=psir(i1,i2,i3,1)
              psir2=psir(i1,i2,i3,2)
              psir3=psir(i1,i2,i3,3)
              psir4=psir(i1,i2,i3,4)
              !potentials + confining term
              pot1=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),1)+cp(i1,i2,i3)
              pot2=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),2)+cp(i1,i2,i3)
              pot3=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),3)+cp(i1,i2,i3)
              pot4=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),4)+cp(i1,i2,i3)

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
        end do
     end do
     !$omp end do

  else !case with nspinor /=4
     do ispinor=1,nspinor
        !$omp do
        do i3=i3s,i3e
           do i2=i2s,i2e
              !thanks to the optional argument the conditional is done at compile time
              if (present(ibyyzz_r)) then
                 i1st=max(i1s,ibyyzz_r(1,i2-15,i3-15)+1) !in bounds coordinates
                 i1et=min(i1e,ibyyzz_r(2,i2-15,i3-15)+1) !in bounds coordinates
              else
                 i1st=i1s
                 i1et=i1e
              end if
              !no need of setting up to zero values outside wavefunction bounds
              do i1=i1st,i1et
                 psir1=psir(i1,i2,i3,ispinor)
                 !the local potential is always real (npot=1) + confining term
                 ii1=i1-ishift(1)
                 ii2=i2-ishift(2)
                 ii3=i3-ishift(3)
                 ttt=cp(i1,i2,i3)
                 pot1=pot(ii1,ii2,ii3,1)+ttt
                 tt11=pot1*psir1
                 if (present(psir_noconf)) then
                     pot1_noconf=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),1)
                     tt11_noconf=pot1_noconf*psir1
                     psir_noconf(i1,i2,i3,ispinor) = tt11_noconf
                 end if
                 if (present(econf)) then
                     econf_p=econf_p+real(ttt*psir1*psir1,wp)
                 end if

                 epot_p=epot_p+real(tt11*psir1,wp)
                 psir(i1,i2,i3,ispinor)=tt11
              end do
           end do
        end do
        !$omp end do
     end do
  end if

  
  !$omp critical
  epot=epot+epot_p
  if (present(econf)) then
      econf=econf+econf_p
  end if
  !$omp end critical
  
  !$omp end parallel

  if (.not. present(confdata)) then
     deallocate(hh,ioffset,rxyzConf)
  end if
contains
  
  !inline the definition of the confining potential
  real(wp) function cp(i1,i2,i3)
    implicit none
    integer, intent(in) :: i1,i2,i3
    !local variables
    real(wp) :: r2
    !to be sure that the conditional is executed at compile time
    if (present(confdata)) then
       r2=(hh(1)*real(i1+ioffset(1),wp)-rxyzConf(1))**2 +&
            (hh(2)*real(i2+ioffset(2),wp)-rxyzConf(2))**2 +&
            (hh(3)*real(i3+ioffset(3),wp)-rxyzConf(3))**2 
       cp=prefac*r2**(potorder/2)
    else
       cp=0.0_wp
    end if

  end function cp

END SUBROUTINE apply_potential_lr

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
