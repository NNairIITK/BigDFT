!> @file 
!!   Routines to do special convolution of linear toolbox
!! @author
!!   Copyright (C) 2011-2012 BigDFT group 
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 
subroutine ConvolQuartic4(iproc, nproc, n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3,  &
           hgrid, offsetx, offsety, offsetz, ibyz_c, ibxz_c, ibxy_c, ibyz_f, ibxz_f, ibxy_f, &
           rxyzConf, potentialPrefac, with_kinetic, cprecr, maxdim, &
           xx_c, xx_f1, xx_f, xy_c, xy_f2, xy_f,  xz_c, xz_f4, xz_f, &
           aeff0array, beff0array, ceff0array, eeff0array, &
           aeff0_2array, beff0_2array, ceff0_2array, eeff0_2array, &
           aeff0_2auxarray, beff0_2auxarray, ceff0_2auxarray, eeff0_2auxarray, &
           xya_c, xyb_c, xyc_c, xye_c, xza_c, xzb_c, xzc_c, xze_c, &
           yza_c, yzb_c, yzc_c, yze_c, xya_f, xyb_f, xyc_f, xye_f, &
           xza_f, xzb_f, xzc_f, xze_f, yza_f, yzb_f, yzc_f, yze_f, &
           aeff0, aeff1, aeff2, aeff3, beff0, beff1, beff2, beff3, &
           ceff0, ceff1, ceff2, ceff3, eeff0, eeff1, eeff2, eeff3, &
           aeff0_2, aeff1_2, aeff2_2, aeff3_2, beff0_2, beff1_2, beff2_2, beff3_2, &
           ceff0_2, ceff1_2, ceff2_2, ceff3_2, eeff0_2, eeff1_2, eeff2_2, eeff3_2, & 
           y_c, y_f)

  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3, offsetx, offsety, offsetz, maxdim
  real(gp),intent(in) :: hgrid, potentialPrefac, cprecr
  logical,intent(in) :: with_kinetic
  real(8),dimension(3) :: rxyzConf
  integer,dimension(2,0:n2,0:n3), intent(in) :: ibyz_c,ibyz_f
  integer,dimension(2,0:n1,0:n3), intent(in) :: ibxz_c,ibxz_f
  integer,dimension(2,0:n1,0:n2), intent(in) :: ibxy_c,ibxy_f
  real(wp),dimension(0:n1,0:n2,0:n3),intent(in) :: xx_c
  real(wp),dimension(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3),intent(in) :: xx_f1
  real(wp),dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3),intent(in) :: xx_f
  real(wp),dimension(0:n2,0:n1,0:n3),intent(in) :: xy_c
  real(wp),dimension(nfl2:nfu2,nfl1:nfu1,nfl3:nfu3),intent(in) :: xy_f2
  real(wp),dimension(7,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3),intent(in) :: xy_f
  real(wp),dimension(0:n3,0:n1,0:n2),intent(in) :: xz_c
  real(wp),dimension(nfl3:nfu3,nfl1:nfu1,nfl2:nfu2),intent(in) :: xz_f4
  real(wp),dimension(7,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2),intent(in) :: xz_f
  real(wp),dimension(-17:17,0:maxdim),intent(in):: aeff0array
  real(wp),dimension(-17:17,0:maxdim),intent(in):: beff0array
  real(wp),dimension(-17:17,0:maxdim),intent(in):: ceff0array
  real(wp),dimension(-14:14,0:maxdim),intent(in):: eeff0array
  real(wp),dimension(-17:17,0:maxdim),intent(in):: aeff0_2array
  real(wp),dimension(-17:17,0:maxdim),intent(in):: beff0_2array
  real(wp),dimension(-17:17,0:maxdim),intent(in):: ceff0_2array
  real(wp),dimension(-14:14,0:maxdim),intent(in):: eeff0_2array
  real(wp),dimension(-17:17,0:maxdim),intent(in):: aeff0_2auxarray
  real(wp),dimension(-17:17,0:maxdim),intent(in):: beff0_2auxarray
  real(wp),dimension(-17:17,0:maxdim),intent(in):: ceff0_2auxarray
  real(wp),dimension(-17:17,0:maxdim),intent(in):: eeff0_2auxarray
  real(wp),dimension(0:n2,0:n1,0:n3):: xya_c, xyb_c, xyc_c, xye_c
  real(wp),dimension(0:n3,0:n1,0:n2):: xza_c, xzb_c, xzc_c, xze_c, yza_c, yzb_c, yzc_c, yze_c
  real(wp),dimension(3,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3):: xya_f
  real(wp),dimension(4,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3):: xyb_f
  real(wp),dimension(3,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3):: xyc_f
  real(wp),dimension(4,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3):: xye_f
  real(wp),dimension(3,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2):: xza_f
  real(wp),dimension(4,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2):: xzb_f
  real(wp),dimension(3,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2):: xzc_f
  real(wp),dimension(4,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2):: xze_f
  real(wp),dimension(3,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2):: yza_f
  real(wp),dimension(4,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2):: yzb_f
  real(wp),dimension(3,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2):: yzc_f
  real(wp),dimension(4,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2):: yze_f
  real(wp),dimension(35):: aeff0, aeff1, aeff2, aeff3, beff0, beff1, beff2, beff3, ceff0, ceff1, ceff2, ceff3
  real(wp),dimension(29):: eeff0, eeff1, eeff2, eeff3
  real(wp),dimension(35):: aeff0_2, aeff1_2, aeff2_2, aeff3_2, beff0_2, beff1_2, beff2_2, beff3_2
  real(wp),dimension(35):: ceff0_2, ceff1_2, ceff2_2, ceff3_2
  real(wp),dimension(29):: eeff0_2, eeff1_2, eeff2_2, eeff3_2
  real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: y_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(out) :: y_f

  ! Local variables
  integer,parameter :: lowfil=-14,lupfil=14
  integer :: i,t,i1,i2,i3, icur,istart,iend,l, istat, iall
  real(wp) :: dyi,dyi0,dyi1,dyi2,dyi3,t112,t121,t122,t212,t221,t222,t211
  real(wp) :: tt112, tt121, tt122, tt212, tt221, tt222, tt211, tt0
  real(kind=8) :: x0, y0, z0
  real(kind=8) :: tt10, tt11, tt12, tt13
  real(kind=8) :: tt20, tt21, tt22, tt23
  real(kind=8) :: tt30, tt31, tt32, tt33
  real(kind=8) :: tt40, tt41, tt42, tt43
  real(kind=8) :: tt50, tt51, tt52, tt53
  real(kind=8) :: tt60, tt61, tt62, tt63
  real(kind=8) :: tt70
  real(kind=8) :: tt0a0, tt0a1, tt0a2, tt0a3
  real(kind=8) :: tt0b0, tt0b1, tt0b2, tt0b3
  real(kind=8) :: tt0c0, tt0c1, tt0c2, tt0c3
  real(kind=8) :: tt0e0, tt0e1, tt0e2, tt0e3
  real(kind=8) :: tt1a0, tt1b0, tt1c0, tt1e0                     
  real(kind=8) :: tt2a0, tt2b0, tt2c0, tt2e0                     
  real(kind=8) :: tt3a0, tt3b0, tt3c0, tt3e0                     
  real(kind=8) :: tt4a0, tt4b0, tt4c0, tt4e0                     
  real(kind=8) :: tt5a0, tt5b0, tt5c0, tt5e0                     
  real(kind=8) :: tt6a0, tt6b0, tt6c0, tt6e0                     
  real(kind=8) :: tt7a0, tt7b0, tt7c0, tt7e0                     
  logical:: with_confpot
  character(len=*),parameter :: subname='ConvolQuartic4'

  real(8)::t2,t1


  call timing(iproc,'convolQuartic ','ON')


  ! Flag indicating whether a confining quartic potential is used
  with_confpot=(potentialPrefac/=0.d0)

 

    do i1=0,n1
        x0=hgrid*(i1+offsetx)-rxyzConf(1)
        if(.not. with_kinetic) then
            call getFilterQuartic(potentialPrefac, hgrid, x0, aeff0array(lowfil,i1), 'a')
            call getFilterQuartic(potentialPrefac, hgrid, x0, beff0array(lowfil,i1), 'b')
            call getFilterQuartic(potentialPrefac, hgrid, x0, ceff0array(lowfil,i1), 'c')
            call getFilterQuartic(potentialPrefac, hgrid, x0, eeff0array(lowfil,i1), 'e')
        else
            call getEffectiveFilterQuartic(potentialPrefac, hgrid, x0, aeff0array(lowfil,i1), 'a')
            call getEffectiveFilterQuartic(potentialPrefac, hgrid, x0, beff0array(lowfil,i1), 'b')
            call getEffectiveFilterQuartic(potentialPrefac, hgrid, x0, ceff0array(lowfil,i1), 'c')
            call getEffectiveFilterQuartic(potentialPrefac, hgrid, x0, eeff0array(lowfil,i1), 'e')
        end if

        call getFilterQuadratic(1.d0, hgrid, x0, aeff0_2auxarray(lowfil,i1), 'a')
        call getFilterQuadratic(1.d0, hgrid, x0, beff0_2auxarray(lowfil,i1), 'b')
        call getFilterQuadratic(1.d0, hgrid, x0, ceff0_2auxarray(lowfil,i1), 'c')
        call getFilterQuadratic(1.d0, hgrid, x0, eeff0_2auxarray(lowfil,i1), 'e')
    end do

  !$omp parallel default(private) &
  !$omp shared(hgrid,offsetx,offsety,offsetz,rxyzConf,with_kinetic,potentialPrefac,with_confpot,cprecr) &
  !$omp shared(aeff0_2auxarray,beff0_2auxarray,ceff0_2auxarray,eeff0_2auxarray,aeff0array,beff0array,ceff0array,eeff0array)&
  !$omp shared(aeff0_2array,beff0_2array,ceff0_2array,eeff0_2array)&
  !$omp shared(nfu1,nfu2,nfu3,n1,n2,n3,nfl1,nfl2,nfl3)&
  !$omp shared(xya_c,xyb_c,xyc_c,xye_c,xza_c,xzb_c,xzc_c,xze_c,yza_c,yzb_c,yzc_c,yze_c)&
  !$omp shared(xya_f,xyb_f,xyc_f,xye_f,xza_f,xzb_f,xzc_f,xze_f,yza_f,yzb_f,yzc_f,yze_f)&
  !$omp shared(ibxy_c,ibxy_f,ibxz_c,ibyz_c,ibxz_f,ibyz_f,xx_c,xx_f,xx_f1,xy_c,xy_f,xz_f,xy_f2,xz_f4,xz_c)&
  !$omp shared(y_c,y_f)

  !$omp do 

    do i3=0,n3
       do i2=0,n2
          if (ibyz_c(2,i2,i3)-ibyz_c(1,i2,i3).ge.4) then
             do i1=ibyz_c(1,i2,i3),ibyz_c(2,i2,i3)-4,4
                dyi0=0.0_wp ; dyi1=0.0_wp ; dyi2=0.0_wp ; dyi3=0.0_wp 
                tt0a0=0.d0 ; tt0a1=0.d0 ; tt0a2=0.d0 ; tt0a3=0.d0
                tt0b0=0.d0 ; tt0b1=0.d0 ; tt0b2=0.d0 ; tt0b3=0.d0
                tt0c0=0.d0 ; tt0c1=0.d0 ; tt0c2=0.d0 ; tt0c3=0.d0
                tt0e0=0.d0 ; tt0e1=0.d0 ; tt0e2=0.d0 ; tt0e3=0.d0
  
                do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_c(2,i2,i3))
                   dyi0=dyi0 + xx_c(t,i2,i3)*aeff0array(t-i1-0,i1+0)
                   dyi1=dyi1 + xx_c(t,i2,i3)*aeff0array(t-i1-1,i1+1)
                   dyi2=dyi2 + xx_c(t,i2,i3)*aeff0array(t-i1-2,i1+2)
                   dyi3=dyi3 + xx_c(t,i2,i3)*aeff0array(t-i1-3,i1+3)
                end do
                y_c(i1+0,i2,i3)=dyi0+cprecr*xx_c(i1+0,i2,i3)
                y_c(i1+1,i2,i3)=dyi1+cprecr*xx_c(i1+1,i2,i3)
                y_c(i1+2,i2,i3)=dyi2+cprecr*xx_c(i1+2,i2,i3)
                y_c(i1+3,i2,i3)=dyi3+cprecr*xx_c(i1+3,i2,i3)

                ! sss coefficients
                if(with_confpot) then
                   do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_c(2,i2,i3))
                      tt0a0=tt0a0 + xx_c(t,i2,i3)*aeff0_2auxarray(t-i1-0,i1+0)
                      tt0a1=tt0a1 + xx_c(t,i2,i3)*aeff0_2auxarray(t-i1-1,i1+1)
                      tt0a2=tt0a2 + xx_c(t,i2,i3)*aeff0_2auxarray(t-i1-2,i1+2)
                      tt0a3=tt0a3 + xx_c(t,i2,i3)*aeff0_2auxarray(t-i1-3,i1+3)

                      tt0b0=tt0b0 + xx_c(t,i2,i3)*beff0_2auxarray(t-i1-0,i1+0)
                      tt0b1=tt0b1 + xx_c(t,i2,i3)*beff0_2auxarray(t-i1-1,i1+1)
                      tt0b2=tt0b2 + xx_c(t,i2,i3)*beff0_2auxarray(t-i1-2,i1+2)
                      tt0b3=tt0b3 + xx_c(t,i2,i3)*beff0_2auxarray(t-i1-3,i1+3)

                      tt0c0=tt0c0 + xx_c(t,i2,i3)*ceff0_2auxarray(t-i1-0,i1+0)
                      tt0c1=tt0c1 + xx_c(t,i2,i3)*ceff0_2auxarray(t-i1-1,i1+1)
                      tt0c2=tt0c2 + xx_c(t,i2,i3)*ceff0_2auxarray(t-i1-2,i1+2)
                      tt0c3=tt0c3 + xx_c(t,i2,i3)*ceff0_2auxarray(t-i1-3,i1+3)

                      tt0e0=tt0e0 + xx_c(t,i2,i3)*eeff0_2auxarray(t-i1-0,i1+0)
                      tt0e1=tt0e1 + xx_c(t,i2,i3)*eeff0_2auxarray(t-i1-1,i1+1)
                      tt0e2=tt0e2 + xx_c(t,i2,i3)*eeff0_2auxarray(t-i1-2,i1+2)
                      tt0e3=tt0e3 + xx_c(t,i2,i3)*eeff0_2auxarray(t-i1-3,i1+3)
                   end do

                   xya_c(i2,i1+0,i3)=tt0a0
                   xya_c(i2,i1+1,i3)=tt0a1
                   xya_c(i2,i1+2,i3)=tt0a2
                   xya_c(i2,i1+3,i3)=tt0a3
                   xza_c(i3,i1+0,i2)=tt0a0
                   xza_c(i3,i1+1,i2)=tt0a1
                   xza_c(i3,i1+2,i2)=tt0a2
                   xza_c(i3,i1+3,i2)=tt0a3

                   xyb_c(i2,i1+0,i3)=tt0b0
                   xyb_c(i2,i1+1,i3)=tt0b1
                   xyb_c(i2,i1+2,i3)=tt0b2
                   xyb_c(i2,i1+3,i3)=tt0b3
                   xzb_c(i3,i1+0,i2)=tt0b0
                   xzb_c(i3,i1+1,i2)=tt0b1
                   xzb_c(i3,i1+2,i2)=tt0b2
                   xzb_c(i3,i1+3,i2)=tt0b3

                   xyc_c(i2,i1+0,i3)=tt0c0
                   xyc_c(i2,i1+1,i3)=tt0c1
                   xyc_c(i2,i1+2,i3)=tt0c2
                   xyc_c(i2,i1+3,i3)=tt0c3
                   xzc_c(i3,i1+0,i2)=tt0c0
                   xzc_c(i3,i1+1,i2)=tt0c1
                   xzc_c(i3,i1+2,i2)=tt0c2
                   xzc_c(i3,i1+3,i2)=tt0c3

                   xye_c(i2,i1+0,i3)=tt0e0
                   xye_c(i2,i1+1,i3)=tt0e1
                   xye_c(i2,i1+2,i3)=tt0e2
                   xye_c(i2,i1+3,i3)=tt0e3
                   xze_c(i3,i1+0,i2)=tt0e0
                   xze_c(i3,i1+1,i2)=tt0e1
                   xze_c(i3,i1+2,i2)=tt0e2
                   xze_c(i3,i1+3,i2)=tt0e3
               end if

             enddo
             icur=i1
          else
             icur=ibyz_c(1,i2,i3)
          endif
  
          do i1=icur,ibyz_c(2,i2,i3)
             dyi=0.0_wp ; tt0a0=0.d0 ; tt0b0=0.d0 ; tt0c0=0.d0 ; tt0e0=0.d0
             do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_c(2,i2,i3))
                dyi=dyi + xx_c(t,i2,i3)*aeff0array(t-i1,i1)
             end do
             y_c(i1,i2,i3)=dyi+cprecr*xx_c(i1,i2,i3)

             if(with_confpot) then
                ! sss coefficients
                do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_c(2,i2,i3))
                   tt0a0=tt0a0 + xx_c(t,i2,i3)*aeff0_2auxarray(t-i1,i1)
                   tt0b0=tt0b0 + xx_c(t,i2,i3)*beff0_2auxarray(t-i1,i1)
                   tt0c0=tt0c0 + xx_c(t,i2,i3)*ceff0_2auxarray(t-i1,i1)
                   tt0e0=tt0e0 + xx_c(t,i2,i3)*eeff0_2auxarray(t-i1,i1)

                   xya_c(i2,i1,i3)=tt0a0
                   xza_c(i3,i1,i2)=tt0a0

                   xyb_c(i2,i1,i3)=tt0b0
                   xzb_c(i3,i1,i2)=tt0b0

                   xyc_c(i2,i1,i3)=tt0c0
                   xzc_c(i3,i1,i2)=tt0c0

                   xye_c(i2,i1,i3)=tt0e0
                   xze_c(i3,i1,i2)=tt0e0
                enddo
            end if

          enddo
  
          istart=max(ibyz_c(1,i2,i3),ibyz_f(1,i2,i3)-lupfil)
          iend=min(ibyz_c(2,i2,i3),ibyz_f(2,i2,i3)-lowfil)
  
          if (istart-iend.ge.4) then
             do i1=istart,iend-4,4
                dyi0=0.0_wp ; dyi1=0.0_wp ; dyi2=0.0_wp ; dyi3=0.0_wp
                do t=max(ibyz_f(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_f(2,i2,i3))
                   dyi0=dyi0 + xx_f1(t,i2,i3)*beff0array(t-i1-0,i1+0)
                   dyi1=dyi1 + xx_f1(t,i2,i3)*beff0array(t-i1-1,i1+1)
                   dyi2=dyi2 + xx_f1(t,i2,i3)*beff0array(t-i1-2,i1+2)
                   dyi3=dyi3 + xx_f1(t,i2,i3)*beff0array(t-i1-3,i1+3)
                enddo
                y_c(i1+0,i2,i3)=y_c(i1+0,i2,i3)+dyi0
                y_c(i1+1,i2,i3)=y_c(i1+1,i2,i3)+dyi1
                y_c(i1+2,i2,i3)=y_c(i1+2,i2,i3)+dyi2
                y_c(i1+3,i2,i3)=y_c(i1+3,i2,i3)+dyi3
             enddo
             istart=i1
          endif
  
          do i1=istart,iend
             dyi=0.0_wp ; tt0=0.0_wp
             do t=max(ibyz_f(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_f(2,i2,i3))
                dyi=dyi + xx_f1(t,i2,i3)*beff0array(t-i1,i1)
             enddo
             y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
          enddo
  
           if (ibyz_c(2,i2,i3)-ibyz_c(1,i2,i3).ge.4) then
             do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)-4,4
                dyi0=0.0_wp ; dyi1=0.0_wp ; dyi2=0.0_wp ; dyi3=0.0_wp 
                do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_c(2,i2,i3))
                   dyi0=dyi0 + xx_c(t,i2,i3)*ceff0array(t-i1-0,i1+0)
                   dyi1=dyi1 + xx_c(t,i2,i3)*ceff0array(t-i1-1,i1+1)
                   dyi2=dyi2 + xx_c(t,i2,i3)*ceff0array(t-i1-2,i1+2)
                   dyi3=dyi3 + xx_c(t,i2,i3)*ceff0array(t-i1-3,i1+3)
                enddo
                y_f(1,i1+0,i2,i3)=dyi0
                y_f(1,i1+1,i2,i3)=dyi1
                y_f(1,i1+2,i2,i3)=dyi2
                y_f(1,i1+3,i2,i3)=dyi3
             enddo
             icur=i1
          else
             icur=ibyz_f(1,i2,i3)
          endif
          do i1=icur,ibyz_f(2,i2,i3)
             dyi=0.0_wp ; tt0=0.0_wp 
             do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_c(2,i2,i3))
                dyi=dyi + xx_c(t,i2,i3)*ceff0array(t-i1,i1)
             enddo
             y_f(1,i1,i2,i3)=dyi
          enddo
       enddo
    enddo
    !$omp end do
    

  
  
  
    ! wavelet part
  
    !$omp do 
    do i3=nfl3,nfu3
       do i2=nfl2,nfu2
          do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
             t112=0.0_wp;t121=0.0_wp;t122=0.0_wp;t212=0.0_wp;t221=0.0_wp;t222=0.0_wp;t211=0.0_wp 
             tt112=0.0_wp;tt121=0.0_wp;tt122=0.0_wp;tt212=0.0_wp;tt221=0.0_wp;tt222=0.0_wp;tt211=0.0_wp 
             tt1a0=0.d0 ; tt1b0=0.d0 ; tt1c0=0.d0 ; tt1e0=0.d0
             tt2a0=0.d0 ; tt2b0=0.d0 ; tt2c0=0.d0 ; tt2e0=0.d0
             tt3a0=0.d0 ; tt3b0=0.d0 ; tt3c0=0.d0 ; tt3e0=0.d0
             tt4a0=0.d0 ; tt4b0=0.d0 ; tt4c0=0.d0 ; tt4e0=0.d0
             tt5a0=0.d0 ; tt5b0=0.d0 ; tt5c0=0.d0 ; tt5e0=0.d0
             tt6a0=0.d0 ; tt6b0=0.d0 ; tt6c0=0.d0 ; tt6e0=0.d0
             tt7a0=0.d0 ; tt7b0=0.d0 ; tt7c0=0.d0 ; tt7e0=0.d0
             do l=max(nfl1-i1,lowfil),min(lupfil,nfu1-i1)
                t112=t112 + xx_f(4,i1+l,i2,i3)*aeff0array(l,i1) + xx_f(5,i1+l,i2,i3)*beff0array(l,i1)
                t121=t121 + xx_f(2,i1+l,i2,i3)*aeff0array(l,i1) + xx_f(3,i1+l,i2,i3)*beff0array(l,i1)
                t122=t122 + xx_f(6,i1+l,i2,i3)*aeff0array(l,i1) + xx_f(7,i1+l,i2,i3)*beff0array(l,i1)
                t212=t212 + xx_f(4,i1+l,i2,i3)*ceff0array(l,i1) + xx_f(5,i1+l,i2,i3)*eeff0array(l,i1)
                t221=t221 + xx_f(2,i1+l,i2,i3)*ceff0array(l,i1) + xx_f(3,i1+l,i2,i3)*eeff0array(l,i1)
                t222=t222 + xx_f(6,i1+l,i2,i3)*ceff0array(l,i1) + xx_f(7,i1+l,i2,i3)*eeff0array(l,i1)
                t211=t211 + xx_f(1,i1+l,i2,i3)*eeff0array(l,i1)
             end do

             y_f(4,i1,i2,i3)=t112+cprecr*xx_f(4,i1,i2,i3)
             y_f(2,i1,i2,i3)=t121+cprecr*xx_f(2,i1,i2,i3)
             y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211+cprecr*xx_f(1,i1,i2,i3)
             y_f(6,i1,i2,i3)=t122+cprecr*xx_f(6,i1,i2,i3)
             y_f(5,i1,i2,i3)=t212+cprecr*xx_f(5,i1,i2,i3)
             y_f(3,i1,i2,i3)=t221+cprecr*xx_f(3,i1,i2,i3)
             y_f(7,i1,i2,i3)=t222+cprecr*xx_f(7,i1,i2,i3)

             if(with_confpot) then
                do l=max(nfl1-i1,lowfil),min(lupfil,nfu1-i1)
                   ! dss coefficients
                   tt1b0=tt1b0 + xx_f(1,i1+l,i2,i3)*beff0_2auxarray(l,i1)
                   tt1e0=tt1e0 + xx_f(1,i1+l,i2,i3)*eeff0_2auxarray(l,i1)
                   ! sds coefficients
                   tt2a0=tt2a0 + xx_f(2,i1+l,i2,i3)*aeff0_2auxarray(l,i1)
                   tt2c0=tt2c0 + xx_f(2,i1+l,i2,i3)*ceff0_2auxarray(l,i1)
                   ! dds coefficients
                   tt3b0=tt3b0 + xx_f(3,i1+l,i2,i3)*beff0_2auxarray(l,i1)
                   tt3e0=tt3e0 + xx_f(3,i1+l,i2,i3)*eeff0_2auxarray(l,i1)
                   ! ssd coefficients
                   tt4a0=tt4a0 + xx_f(4,i1+l,i2,i3)*aeff0_2auxarray(l,i1)
                   tt4c0=tt4c0 + xx_f(4,i1+l,i2,i3)*ceff0_2auxarray(l,i1)
                   ! dsd coefficients
                   tt5b0=tt5b0 + xx_f(5,i1+l,i2,i3)*beff0_2auxarray(l,i1)
                   tt5e0=tt5e0 + xx_f(5,i1+l,i2,i3)*eeff0_2auxarray(l,i1)
                   ! sdd coefficients
                   tt6a0=tt6a0 + xx_f(6,i1+l,i2,i3)*aeff0_2auxarray(l,i1)
                   tt6c0=tt6c0 + xx_f(6,i1+l,i2,i3)*ceff0_2auxarray(l,i1)
                   ! ddd coefficients
                   tt7b0=tt7b0 + xx_f(7,i1+l,i2,i3)*beff0_2auxarray(l,i1)
                   tt7e0=tt7e0 + xx_f(7,i1+l,i2,i3)*eeff0_2auxarray(l,i1)
                enddo

                ! dss coefficients
                xyb_f(1,i2,i1,i3)=tt1b0
                xye_f(1,i2,i1,i3)=tt1e0
                xzb_f(1,i3,i1,i2)=tt1b0
                xze_f(1,i3,i1,i2)=tt1e0
                ! sds coefficients
                xya_f(1,i2,i1,i3)=tt2a0
                xyc_f(1,i2,i1,i3)=tt2c0
                xza_f(1,i3,i1,i2)=tt2a0
                xzc_f(1,i3,i1,i2)=tt2c0
                ! dds coefficients
                !xyb_f3(i2,i1,i3)=tt3b0
                xyb_f(2,i2,i1,i3)=tt3b0
                xye_f(2,i2,i1,i3)=tt3e0
                xzb_f(2,i3,i1,i2)=tt3b0
                xze_f(2,i3,i1,i2)=tt3e0
                ! ssd coefficients
                xya_f(2,i2,i1,i3)=tt4a0
                xyc_f(2,i2,i1,i3)=tt4c0
                xza_f(2,i3,i1,i2)=tt4a0
                xzc_f(2,i3,i1,i2)=tt4c0
                ! dsd coefficients
                xyb_f(3,i2,i1,i3)=tt5b0
                xye_f(3,i2,i1,i3)=tt5e0
                xzb_f(3,i3,i1,i2)=tt5b0
                xze_f(3,i3,i1,i2)=tt5e0
                ! sdd coefficients
                xya_f(3,i2,i1,i3)=tt6a0
                xyc_f(3,i2,i1,i3)=tt6c0
                xza_f(3,i3,i1,i2)=tt6a0
                xzc_f(3,i3,i1,i2)=tt6c0
                ! sdd coefficients
                xyb_f(4,i2,i1,i3)=tt7b0
                xye_f(4,i2,i1,i3)=tt7e0
                xzb_f(4,i3,i1,i2)=tt7b0
                xze_f(4,i3,i1,i2)=tt7e0
             end if
          enddo
       enddo
    enddo
    !$omp enddo
    

    !$omp single
    do i2=0,n2
        y0=hgrid*(i2+offsety)-rxyzConf(2)
        if(.not. with_kinetic) then
            call getFilterQuartic(potentialPrefac, hgrid, y0, aeff0array(lowfil,i2), 'a')
            call getFilterQuartic(potentialPrefac, hgrid, y0, beff0array(lowfil,i2), 'b')
            call getFilterQuartic(potentialPrefac, hgrid, y0, ceff0array(lowfil,i2), 'c')
            call getFilterQuartic(potentialPrefac, hgrid, y0, eeff0array(lowfil,i2), 'e')
        else
            call getEffectiveFilterQuartic(potentialPrefac, hgrid, y0, aeff0array(lowfil,i2), 'a')
            call getEffectiveFilterQuartic(potentialPrefac, hgrid, y0, beff0array(lowfil,i2), 'b')
            call getEffectiveFilterQuartic(potentialPrefac, hgrid, y0, ceff0array(lowfil,i2), 'c')
            call getEffectiveFilterQuartic(potentialPrefac, hgrid, y0, eeff0array(lowfil,i2), 'e')
        end if

        call getFilterQuadratic(potentialPrefac, hgrid, y0, aeff0_2array(lowfil,i2), 'a')
        call getFilterQuadratic(potentialPrefac, hgrid, y0, beff0_2array(lowfil,i2), 'b')
        call getFilterQuadratic(potentialPrefac, hgrid, y0, ceff0_2array(lowfil,i2), 'c')
        call getFilterQuadratic(potentialPrefac, hgrid, y0, eeff0_2array(lowfil,i2), 'e')

        call getFilterQuadratic(1.d0, hgrid, y0, aeff0_2auxarray(lowfil,i2), 'a')
        call getFilterQuadratic(1.d0, hgrid, y0, beff0_2auxarray(lowfil,i2), 'b')
        call getFilterQuadratic(1.d0, hgrid, y0, ceff0_2auxarray(lowfil,i2), 'c')
        call getFilterQuadratic(1.d0, hgrid, y0, eeff0_2auxarray(lowfil,i2), 'e')
    end do
    !$omp end single
  
    ! + (1/2) d^2/dy^2
    !$omp do 
    do i3=0,n3
       do i1=0,n1
          if (ibxz_c(2,i1,i3)-ibxz_c(1,i1,i3).ge.4) then
             do i2=ibxz_c(1,i1,i3),ibxz_c(2,i1,i3)-4,4
                dyi0=0.0_wp ;  dyi1=0.0_wp ;  dyi2=0.0_wp ;  dyi3=0.0_wp 
                tt0a0=0.d0 ; tt0a1=0.d0 ; tt0a2=0.d0 ; tt0a3=0.d0
                tt0b0=0.d0 ; tt0b1=0.d0 ; tt0b2=0.d0 ; tt0b3=0.d0
                tt0c0=0.d0 ; tt0c1=0.d0 ; tt0c2=0.d0 ; tt0c3=0.d0
                tt0e0=0.d0 ; tt0e1=0.d0 ; tt0e2=0.d0 ; tt0e3=0.d0
                if(with_confpot) then
                   do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_c(2,i1,i3))
                      dyi0=dyi0 + xy_c(t,i1,i3)*aeff0array(t-i2-0,i2+0) + 2.d0*xya_c(t,i1,i3)*aeff0_2array(t-i2-0,i2+0)
                      dyi1=dyi1 + xy_c(t,i1,i3)*aeff0array(t-i2-1,i2+1) + 2.d0*xya_c(t,i1,i3)*aeff0_2array(t-i2-1,i2+1)
                      dyi2=dyi2 + xy_c(t,i1,i3)*aeff0array(t-i2-2,i2+2) + 2.d0*xya_c(t,i1,i3)*aeff0_2array(t-i2-2,i2+2)
                      dyi3=dyi3 + xy_c(t,i1,i3)*aeff0array(t-i2-3,i2+3) + 2.d0*xya_c(t,i1,i3)*aeff0_2array(t-i2-3,i2+3)
                   end do
                else
                   do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_c(2,i1,i3))
                      dyi0=dyi0 + xy_c(t,i1,i3)*aeff0array(t-i2-0,i2+0)
                      dyi1=dyi1 + xy_c(t,i1,i3)*aeff0array(t-i2-1,i2+1)
                      dyi2=dyi2 + xy_c(t,i1,i3)*aeff0array(t-i2-2,i2+2)
                      dyi3=dyi3 + xy_c(t,i1,i3)*aeff0array(t-i2-3,i2+3)
                   end do
                end if
                y_c(i1,i2+0,i3)=y_c(i1,i2+0,i3)+dyi0
                y_c(i1,i2+1,i3)=y_c(i1,i2+1,i3)+dyi1
                y_c(i1,i2+2,i3)=y_c(i1,i2+2,i3)+dyi2
                y_c(i1,i2+3,i3)=y_c(i1,i2+3,i3)+dyi3

                if(with_confpot) then
                   do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_c(2,i1,i3))
                      ! sss coefficients
                      tt0a0=tt0a0 + xy_c(t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)
                      tt0a1=tt0a1 + xy_c(t,i1,i3)*aeff0_2auxarray(t-i2-1,i2+1)
                      tt0a2=tt0a2 + xy_c(t,i1,i3)*aeff0_2auxarray(t-i2-2,i2+2)
                      tt0a3=tt0a3 + xy_c(t,i1,i3)*aeff0_2auxarray(t-i2-3,i2+3)

                      tt0b0=tt0b0 + xy_c(t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)
                      tt0b1=tt0b1 + xy_c(t,i1,i3)*beff0_2auxarray(t-i2-1,i2+1)
                      tt0b2=tt0b2 + xy_c(t,i1,i3)*beff0_2auxarray(t-i2-2,i2+2)
                      tt0b3=tt0b3 + xy_c(t,i1,i3)*beff0_2auxarray(t-i2-3,i2+3)

                      tt0c0=tt0c0 + xy_c(t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)
                      tt0c1=tt0c1 + xy_c(t,i1,i3)*ceff0_2auxarray(t-i2-1,i2+1)
                      tt0c2=tt0c2 + xy_c(t,i1,i3)*ceff0_2auxarray(t-i2-2,i2+2)
                      tt0c3=tt0c3 + xy_c(t,i1,i3)*ceff0_2auxarray(t-i2-3,i2+3)

                      tt0e0=tt0e0 + xy_c(t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)
                      tt0e1=tt0e1 + xy_c(t,i1,i3)*eeff0_2auxarray(t-i2-1,i2+1)
                      tt0e2=tt0e2 + xy_c(t,i1,i3)*eeff0_2auxarray(t-i2-2,i2+2)
                      tt0e3=tt0e3 + xy_c(t,i1,i3)*eeff0_2auxarray(t-i2-3,i2+3)
                   enddo

                   yza_c(i3,i1,i2+0)=tt0a0
                   yza_c(i3,i1,i2+1)=tt0a1
                   yza_c(i3,i1,i2+2)=tt0a2
                   yza_c(i3,i1,i2+3)=tt0a3
                             
                   yzb_c(i3,i1,i2+0)=tt0b0
                   yzb_c(i3,i1,i2+1)=tt0b1
                   yzb_c(i3,i1,i2+2)=tt0b2
                   yzb_c(i3,i1,i2+3)=tt0b3
                             
                   yzc_c(i3,i1,i2+0)=tt0c0
                   yzc_c(i3,i1,i2+1)=tt0c1
                   yzc_c(i3,i1,i2+2)=tt0c2
                   yzc_c(i3,i1,i2+3)=tt0c3
                             
                   yze_c(i3,i1,i2+0)=tt0e0
                   yze_c(i3,i1,i2+1)=tt0e1
                   yze_c(i3,i1,i2+2)=tt0e2
                   yze_c(i3,i1,i2+3)=tt0e3
                end if
             enddo
             icur=i2
          else
             icur=ibxz_c(1,i1,i3)
          endif
  
          do i2=icur,ibxz_c(2,i1,i3)
             dyi=0.0_wp ; tt0a0=0.d0 ; tt0b0=0.d0 ; tt0c0=0.d0 ; tt0e0=0.d0
             if(with_confpot) then
                do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
                   dyi=dyi + xy_c(t,i1,i3)*aeff0array(t-i2,i2) + 2.d0*xya_c(t,i1,i3)*aeff0_2array(t-i2,i2)
                end do
             else
                do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
                   dyi=dyi + xy_c(t,i1,i3)*aeff0array(t-i2,i2)
                end do
             end if
             y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi

             if(with_confpot) then
                do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
                   tt0a0=tt0a0 + xy_c(t,i1,i3)*aeff0_2auxarray(t-i2-0,i2)

                   tt0b0=tt0b0 + xy_c(t,i1,i3)*beff0_2auxarray(t-i2-0,i2)

                   tt0c0=tt0c0 + xy_c(t,i1,i3)*ceff0_2auxarray(t-i2-0,i2)

                   tt0e0=tt0e0 + xy_c(t,i1,i3)*eeff0_2auxarray(t-i2-0,i2)
                enddo

                yza_c(i3,i1,i2)=tt0a0

                yzb_c(i3,i1,i2)=tt0b0

                yzc_c(i3,i1,i2)=tt0c0

                yze_c(i3,i1,i2)=tt0e0
             end if
          enddo
          istart=max(ibxz_c(1,i1,i3),ibxz_f(1,i1,i3)-lupfil)
          iend= min(ibxz_c(2,i1,i3),ibxz_f(2,i1,i3)-lowfil)
  
          if (istart-iend.ge.4) then
             do i2=istart,iend-4,4
                dyi0=0.0_wp ; dyi1=0.0_wp ; dyi2=0.0_wp ; dyi3=0.0_wp
                if(with_confpot) then
                   do t=max(ibxz_f(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_f(2,i1,i3))
                      dyi0=dyi0 + xy_f2(t,i1,i3)*beff0array(t-i2-0,i2+0) + 2.d0*xyb_f(1,t,i1,i3)*aeff0_2array(t-i2-0,i2+0) + &
                                  2.d0*(xya_f(1,t,i1,i3)+xyb_f(2,t,i1,i3))*beff0_2array(t-i2-0,i2+0)
                      dyi1=dyi1 + xy_f2(t,i1,i3)*beff0array(t-i2-1,i2+1) + 2.d0*xyb_f(1,t,i1,i3)*aeff0_2array(t-i2-1,i2+1) + &
                                  2.d0*(xya_f(1,t,i1,i3)+xyb_f(2,t,i1,i3))*beff0_2array(t-i2-1,i2+1)
                      dyi2=dyi2 + xy_f2(t,i1,i3)*beff0array(t-i2-2,i2+2) + 2.d0*xyb_f(1,t,i1,i3)*aeff0_2array(t-i2-2,i2+2) + &
                                  2.d0*(xya_f(1,t,i1,i3)+xyb_f(2,t,i1,i3))*beff0_2array(t-i2-2,i2+2)
                      dyi3=dyi3 + xy_f2(t,i1,i3)*beff0array(t-i2-3,i2+3) + 2.d0*xyb_f(1,t,i1,i3)*aeff0_2array(t-i2-3,i2+3) + &
                                  2.d0*(xya_f(1,t,i1,i3)+xyb_f(2,t,i1,i3))*beff0_2array(t-i2-3,i2+3)
                   enddo
                else
                   do t=max(ibxz_f(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_f(2,i1,i3))
                      dyi0=dyi0 + xy_f2(t,i1,i3)*beff0array(t-i2-0,i2+0)
                      dyi1=dyi1 + xy_f2(t,i1,i3)*beff0array(t-i2-1,i2+1)
                      dyi2=dyi2 + xy_f2(t,i1,i3)*beff0array(t-i2-2,i2+2)
                      dyi3=dyi3 + xy_f2(t,i1,i3)*beff0array(t-i2-3,i2+3)
                   enddo
                end if
                y_c(i1,i2+0,i3)=y_c(i1,i2+0,i3)+dyi0
                y_c(i1,i2+1,i3)=y_c(i1,i2+1,i3)+dyi1
                y_c(i1,i2+2,i3)=y_c(i1,i2+2,i3)+dyi2
                y_c(i1,i2+3,i3)=y_c(i1,i2+3,i3)+dyi3
             enddo
             istart=i2
          endif
  
          do i2=istart,iend
             dyi0=0.0_wp
             if(with_confpot) then
                do t=max(ibxz_f(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_f(2,i1,i3))
                   dyi0=dyi0 + xy_f2(t,i1,i3)*beff0array(t-i2-0,i2) + 2.d0*xyb_f(1,t,i1,i3)*aeff0_2array(t-i2-0,i2) + &
                               2.d0*(xya_f(1,t,i1,i3)+xyb_f(2,t,i1,i3))*beff0_2array(t-i2-0,i2)
                enddo
             else
                do t=max(ibxz_f(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_f(2,i1,i3))
                   dyi0=dyi0 + xy_f2(t,i1,i3)*beff0array(t-i2-0,i2)
                enddo
             end if
             y_c(i1,i2+0,i3)=y_c(i1,i2+0,i3)+dyi0
          enddo
  
           if (ibxz_f(2,i1,i3)-ibxz_f(1,i1,i3).ge.4) then
             do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)-4,4
                dyi0=0.0_wp ; dyi1=0.0_wp ; dyi2=0.0_wp ; dyi3=0.0_wp 
                tt10=0.0_wp ; tt11=0.0_wp ; tt12=0.0_wp ; tt13=0.0_wp 
                tt20=0.0_wp ; tt21=0.0_wp ; tt22=0.0_wp ; tt23=0.0_wp 
                tt30=0.0_wp ; tt31=0.0_wp ; tt32=0.0_wp ; tt33=0.0_wp 
                do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_c(2,i1,i3))
                   dyi0=dyi0 + xy_c(t,i1,i3)*ceff0array(t-i2-0,i2+0)
                   dyi1=dyi1 + xy_c(t,i1,i3)*ceff0array(t-i2-1,i2+1)
                   dyi2=dyi2 + xy_c(t,i1,i3)*ceff0array(t-i2-2,i2+2)
                   dyi3=dyi3 + xy_c(t,i1,i3)*ceff0array(t-i2-3,i2+3)
                end do
                y_f(2,i1,i2+0,i3)=y_f(2,i1,i2+0,i3)+dyi0
                y_f(2,i1,i2+1,i3)=y_f(2,i1,i2+1,i3)+dyi1
                y_f(2,i1,i2+2,i3)=y_f(2,i1,i2+2,i3)+dyi2
                y_f(2,i1,i2+3,i3)=y_f(2,i1,i2+3,i3)+dyi3
  
                if(with_confpot) then
                   do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_c(2,i1,i3))
                      tt10=tt10 + 2.d0*xyc_c(t,i1,i3)*aeff0_2array(t-i2-0,i2+0)
                      tt11=tt11 + 2.d0*xyc_c(t,i1,i3)*aeff0_2array(t-i2-1,i2+1)
                      tt12=tt12 + 2.d0*xyc_c(t,i1,i3)*aeff0_2array(t-i2-2,i2+2)
                      tt13=tt13 + 2.d0*xyc_c(t,i1,i3)*aeff0_2array(t-i2-3,i2+3)
  
                      tt20=tt20 + 2.d0*xya_c(t,i1,i3)*ceff0_2array(t-i2-0,i2+0)
                      tt21=tt21 + 2.d0*xya_c(t,i1,i3)*ceff0_2array(t-i2-1,i2+1)
                      tt22=tt22 + 2.d0*xya_c(t,i1,i3)*ceff0_2array(t-i2-2,i2+2)
                      tt23=tt23 + 2.d0*xya_c(t,i1,i3)*ceff0_2array(t-i2-3,i2+3)
  
                      tt30=tt30 + 2.d0*xyc_c(t,i1,i3)*ceff0_2array(t-i2-0,i2+0)
                      tt31=tt31 + 2.d0*xyc_c(t,i1,i3)*ceff0_2array(t-i2-1,i2+1)
                      tt32=tt32 + 2.d0*xyc_c(t,i1,i3)*ceff0_2array(t-i2-2,i2+2)
                      tt33=tt33 + 2.d0*xyc_c(t,i1,i3)*ceff0_2array(t-i2-3,i2+3)
                   enddo
  
                   y_f(1,i1,i2+0,i3)=y_f(1,i1,i2+0,i3)+tt10
                   y_f(1,i1,i2+1,i3)=y_f(1,i1,i2+1,i3)+tt11
                   y_f(1,i1,i2+2,i3)=y_f(1,i1,i2+2,i3)+tt12
                   y_f(1,i1,i2+3,i3)=y_f(1,i1,i2+3,i3)+tt13
  
                   y_f(2,i1,i2+0,i3)=y_f(2,i1,i2+0,i3)+tt20
                   y_f(2,i1,i2+1,i3)=y_f(2,i1,i2+1,i3)+tt21
                   y_f(2,i1,i2+2,i3)=y_f(2,i1,i2+2,i3)+tt22
                   y_f(2,i1,i2+3,i3)=y_f(2,i1,i2+3,i3)+tt23
  
                   y_f(3,i1,i2+0,i3)=y_f(3,i1,i2+0,i3)+tt30
                   y_f(3,i1,i2+1,i3)=y_f(3,i1,i2+1,i3)+tt31
                   y_f(3,i1,i2+2,i3)=y_f(3,i1,i2+2,i3)+tt32
                   y_f(3,i1,i2+3,i3)=y_f(3,i1,i2+3,i3)+tt33
                end if
             enddo
             icur=i2
          else
             icur=ibxz_f(1,i1,i3)
          endif
  
          do i2=icur,ibxz_f(2,i1,i3)
             dyi0=0.0_wp ; tt10=0.0_wp ; tt20=0.0_wp ; tt30=0.0_wp 
             do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
                dyi0=dyi0 + xy_c(t,i1,i3)*ceff0array(t-i2-0,i2)
             end do
             y_f(2,i1,i2+0,i3)=y_f(2,i1,i2+0,i3)+dyi0

             if(with_confpot) then
                do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
                   tt10=tt10 + 2.d0*xyc_c(t,i1,i3)*aeff0_2array(t-i2-0,i2)

                   tt20=tt20 + 2.d0*xya_c(t,i1,i3)*ceff0_2array(t-i2-0,i2)

                   tt30=tt30 + 2.d0*xyc_c(t,i1,i3)*ceff0_2array(t-i2-0,i2)
                enddo
                y_f(1,i1,i2+0,i3)=y_f(1,i1,i2+0,i3)+tt10
                y_f(2,i1,i2+0,i3)=y_f(2,i1,i2+0,i3)+tt20
                y_f(3,i1,i2+0,i3)=y_f(3,i1,i2+0,i3)+tt30
             end if
          enddo
       enddo
    enddo
    !$omp end do
    


    ! wavelet part
  
    !$omp do 
    do i3=nfl3,nfu3
       do i1=nfl1,nfu1
          do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
             ! Get the effective filters for the y dimension
             tt10 = 0.d0 ; tt20 = 0.d0 ; tt30 = 0.d0 ; tt40 = 0.d0 ; tt50 = 0.d0 ; tt60 = 0.d0 ; tt70 = 0.d0
             tt1a0=0.d0 ; tt1b0=0.d0 ; tt1c0=0.d0 ; tt1e0=0.d0
             tt2a0=0.d0 ; tt2b0=0.d0 ; tt2c0=0.d0 ; tt2e0=0.d0
             tt3a0=0.d0 ; tt3b0=0.d0 ; tt3c0=0.d0 ; tt3e0=0.d0
             tt4a0=0.d0 ; tt4b0=0.d0 ; tt4c0=0.d0 ; tt4e0=0.d0
             tt5a0=0.d0 ; tt5b0=0.d0 ; tt5c0=0.d0 ; tt5e0=0.d0
             tt6a0=0.d0 ; tt6b0=0.d0 ; tt6c0=0.d0 ; tt6e0=0.d0
             tt7a0=0.d0 ; tt7b0=0.d0 ; tt7c0=0.d0 ; tt7e0=0.d0
             if(with_confpot) then
                do l=max(nfl2-i2,lowfil),min(lupfil,nfu2-i2)
                   tt10 = tt10 + xy_f(1,i2+l,i1,i3)*aeff0array(l,i2) + xy_f(3,i2+l,i1,i3)*beff0array(l,i2) + &
                                 2.d0*xye_f(1,i2+l,i1,i3)*aeff0_2array(l,i2) + &
                                 2.d0*(xyc_f(1,i2+l,i1,i3)+xye_f(2,i2+l,i1,i3))*beff0_2array(l,i2)
  
                   tt20 = tt20 + xy_f(2,i2+l,i1,i3)*eeff0array(l,i2) +                                      &
                                 2.d0*xyb_f(1,i2+l,i1,i3)*ceff0_2array(l,i2) + &
                                 2.d0*(xya_f(1,i2+l,i1,i3)+xyb_f(2,i2+l,i1,i3))*eeff0_2array(l,i2)
  
                   tt30 = tt30 + xy_f(1,i2+l,i1,i3)*ceff0array(l,i2) + xy_f(3,i2+l,i1,i3)*eeff0array(l,i2) + &
                                 2.d0*xye_f(1,i2+l,i1,i3)*ceff0_2array(l,i2) + &
                                 2.d0*(xyc_f(1,i2+l,i1,i3)+xye_f(2,i2+l,i1,i3))*eeff0_2array(l,i2)
  
                   tt40 = tt40 + xy_f(4,i2+l,i1,i3)*aeff0array(l,i2) + xy_f(6,i2+l,i1,i3)*beff0array(l,i2) + &
                                 2.d0*(xya_f(2,i2+l,i1,i3)+xyb_f(3,i2+l,i1,i3))*aeff0_2array(l,i2) + &
                                 2.d0*(xya_f(3,i2+l,i1,i3)+xyb_f(4,i2+l,i1,i3))*beff0_2array(l,i2)
  
                   tt50 = tt50 + xy_f(5,i2+l,i1,i3)*aeff0array(l,i2) + xy_f(7,i2+l,i1,i3)*beff0array(l,i2) + &
                                 2.d0*(xyc_f(2,i2+l,i1,i3)+xye_f(3,i2+l,i1,i3))*aeff0_2array(l,i2) + &
                                 2.d0*(xyc_f(3,i2+l,i1,i3)+xye_f(4,i2+l,i1,i3))*beff0_2array(l,i2)
                   
                   tt60 = tt60 + xy_f(4,i2+l,i1,i3)*ceff0array(l,i2) + xy_f(6,i2+l,i1,i3)*eeff0array(l,i2) + &
                                 2.d0*(xya_f(2,i2+l,i1,i3)+xyb_f(3,i2+l,i1,i3))*ceff0_2array(l,i2) + &
                                 2.d0*(xya_f(3,i2+l,i1,i3)+xyb_f(4,i2+l,i1,i3))*eeff0_2array(l,i2)
  
                   tt70 = tt70 + xy_f(5,i2+l,i1,i3)*ceff0array(l,i2) + xy_f(7,i2+l,i1,i3)*eeff0array(l,i2) + &
                                 2.d0*(xyc_f(2,i2+l,i1,i3)+xye_f(3,i2+l,i1,i3))*ceff0_2array(l,i2) + &
                                 2.d0*(xyc_f(3,i2+l,i1,i3)+xye_f(4,i2+l,i1,i3))*eeff0_2array(l,i2)
                end do
             else
                do l=max(nfl2-i2,lowfil),min(lupfil,nfu2-i2)
                   tt10 = tt10 + xy_f(1,i2+l,i1,i3)*aeff0array(l,i2) + xy_f(3,i2+l,i1,i3)*beff0array(l,i2)
  
                   tt20 = tt20 + xy_f(2,i2+l,i1,i3)*eeff0array(l,i2)
  
                   tt30 = tt30 + xy_f(1,i2+l,i1,i3)*ceff0array(l,i2) + xy_f(3,i2+l,i1,i3)*eeff0array(l,i2)
  
                   tt40 = tt40 + xy_f(4,i2+l,i1,i3)*aeff0array(l,i2) + xy_f(6,i2+l,i1,i3)*beff0array(l,i2)
  
                   tt50 = tt50 + xy_f(5,i2+l,i1,i3)*aeff0array(l,i2) + xy_f(7,i2+l,i1,i3)*beff0array(l,i2)
                   
                   tt60 = tt60 + xy_f(4,i2+l,i1,i3)*ceff0array(l,i2) + xy_f(6,i2+l,i1,i3)*eeff0array(l,i2)
  
                   tt70 = tt70 + xy_f(5,i2+l,i1,i3)*ceff0array(l,i2) + xy_f(7,i2+l,i1,i3)*eeff0array(l,i2)
                end do
             end if
             y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+tt40
             y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+tt20
             y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+tt10
             y_f(6,i1,i2,i3)=y_f(6,i1,i2,i3)+tt60
             y_f(5,i1,i2,i3)=y_f(5,i1,i2,i3)+tt50
             y_f(3,i1,i2,i3)=y_f(3,i1,i2,i3)+tt30
             y_f(7,i1,i2,i3)=y_f(7,i1,i2,i3)+tt70

             if(with_confpot) then
                do l=max(nfl2-i2,lowfil),min(lupfil,nfu2-i2)
                   ! dss coefficients
                   tt1a0=tt1a0 + xy_f(1,i2+l,i1,i3)*aeff0_2auxarray(l,i2)
                   tt1c0=tt1c0 + xy_f(1,i2+l,i1,i3)*ceff0_2auxarray(l,i2)
                   ! sds coefficients
                   tt2b0=tt2b0 + xy_f(2,i2+l,i1,i3)*beff0_2auxarray(l,i2)
                   tt2e0=tt2e0 + xy_f(2,i2+l,i1,i3)*eeff0_2auxarray(l,i2)
                   ! dds coefficients
                   tt3b0=tt3b0 + xy_f(3,i2+l,i1,i3)*beff0_2auxarray(l,i2)
                   tt3e0=tt3e0 + xy_f(3,i2+l,i1,i3)*eeff0_2auxarray(l,i2)
                   ! ssd coefficients
                   tt4a0=tt4a0 + xy_f(4,i2+l,i1,i3)*aeff0_2auxarray(l,i2)
                   tt4c0=tt4c0 + xy_f(4,i2+l,i1,i3)*ceff0_2auxarray(l,i2)
                   ! dsd coefficients
                   tt5a0=tt5a0 + xy_f(5,i2+l,i1,i3)*aeff0_2auxarray(l,i2)
                   tt5c0=tt5c0 + xy_f(5,i2+l,i1,i3)*ceff0_2auxarray(l,i2)
                   ! sdd coefficients
                   tt6b0=tt6b0 + xy_f(6,i2+l,i1,i3)*beff0_2auxarray(l,i2)
                   tt6e0=tt6e0 + xy_f(6,i2+l,i1,i3)*eeff0_2auxarray(l,i2)
                   ! ddd coefficients
                   tt7b0=tt7b0 + xy_f(7,i2+l,i1,i3)*beff0_2auxarray(l,i2)
                   tt7e0=tt7e0 + xy_f(7,i2+l,i1,i3)*eeff0_2auxarray(l,i2)
                enddo

                ! dss coefficients
                yza_f(1,i3,i1,i2)=tt1a0
                yzc_f(1,i3,i1,i2)=tt1c0
                ! sds coefficients
                yzb_f(1,i3,i1,i2)=tt2b0
                yze_f(1,i3,i1,i2)=tt2e0
                ! dds coefficients
                yzb_f(2,i3,i1,i2)=tt3b0
                yze_f(2,i3,i1,i2)=tt3e0
                ! ssd coefficients
                yza_f(2,i3,i1,i2)=tt4a0
                yzc_f(2,i3,i1,i2)=tt4c0
                ! dsd coefficients
                yza_f(3,i3,i1,i2)=tt5a0
                yzc_f(3,i3,i1,i2)=tt5c0
                ! sdd coefficients
                yzb_f(3,i3,i1,i2)=tt6b0
                yze_f(3,i3,i1,i2)=tt6e0
                ! sdd coefficients
                yzb_f(4,i3,i1,i2)=tt7b0
                yze_f(4,i3,i1,i2)=tt7e0
            end if
          enddo
       enddo
    enddo
    !$omp enddo


  
    !$omp single
    do i3=0,n3
        z0=hgrid*(i3+offsetz)-rxyzConf(3)
        if(.not. with_kinetic) then
            call getFilterQuartic(potentialPrefac, hgrid, z0, aeff0array(lowfil,i3), 'a')
            call getFilterQuartic(potentialPrefac, hgrid, z0, beff0array(lowfil,i3), 'b')
            call getFilterQuartic(potentialPrefac, hgrid, z0, ceff0array(lowfil,i3), 'c')
            call getFilterQuartic(potentialPrefac, hgrid, z0, eeff0array(lowfil,i3), 'e')
        else
            call getEffectiveFilterQuartic(potentialPrefac, hgrid, z0, aeff0array(lowfil,i3), 'a')
            call getEffectiveFilterQuartic(potentialPrefac, hgrid, z0, beff0array(lowfil,i3), 'b')
            call getEffectiveFilterQuartic(potentialPrefac, hgrid, z0, ceff0array(lowfil,i3), 'c')
            call getEffectiveFilterQuartic(potentialPrefac, hgrid, z0, eeff0array(lowfil,i3), 'e')
        end if

        call getFilterQuadratic(potentialPrefac, hgrid, z0, aeff0_2array(lowfil,i3), 'a')
        call getFilterQuadratic(potentialPrefac, hgrid, z0, beff0_2array(lowfil,i3), 'b')
        call getFilterQuadratic(potentialPrefac, hgrid, z0, ceff0_2array(lowfil,i3), 'c')
        call getFilterQuadratic(potentialPrefac, hgrid, z0, eeff0_2array(lowfil,i3), 'e')
    end do
    !$omp end single

  ! + (1/2) d^2/dz^2
  !$omp do 
  do i2=0,n2
     do i1=0,n1
        if (ibxy_c(2,i1,i2)-ibxy_c(1,i1,i2).ge.4) then
           do i3=ibxy_c(1,i1,i2),ibxy_c(2,i1,i2)-4,4
              dyi0=0.0_wp ; dyi1=0.0_wp ; dyi2=0.0_wp ; dyi3=0.0_wp 
              if(with_confpot) then
                 do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_c(2,i1,i2))
                    dyi0=dyi0 + xz_c(t,i1,i2)*aeff0array(t-i3-0,i3+0) + &
                                2.d0*(xza_c(t,i1,i2)+yza_c(t,i1,i2))*aeff0_2array(t-i3-0,i3+0)
                    dyi1=dyi1 + xz_c(t,i1,i2)*aeff0array(t-i3-1,i3+1) + &
                                2.d0*(xza_c(t,i1,i2)+yza_c(t,i1,i2))*aeff0_2array(t-i3-1,i3+1)
                    dyi2=dyi2 + xz_c(t,i1,i2)*aeff0array(t-i3-2,i3+2) + &
                                2.d0*(xza_c(t,i1,i2)+yza_c(t,i1,i2))*aeff0_2array(t-i3-2,i3+2)
                    dyi3=dyi3 + xz_c(t,i1,i2)*aeff0array(t-i3-3,i3+3) + &
                                2.d0*(xza_c(t,i1,i2)+yza_c(t,i1,i2))*aeff0_2array(t-i3-3,i3+3)
                 enddo
              else
                 do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_c(2,i1,i2))
                    dyi0=dyi0 + xz_c(t,i1,i2)*aeff0array(t-i3-0,i3+0)
                    dyi1=dyi1 + xz_c(t,i1,i2)*aeff0array(t-i3-1,i3+1)
                    dyi2=dyi2 + xz_c(t,i1,i2)*aeff0array(t-i3-2,i3+2)
                    dyi3=dyi3 + xz_c(t,i1,i2)*aeff0array(t-i3-3,i3+3)
                 enddo
              end if
              y_c(i1,i2,i3+0)=y_c(i1,i2,i3+0)+dyi0
              y_c(i1,i2,i3+1)=y_c(i1,i2,i3+1)+dyi1
              y_c(i1,i2,i3+2)=y_c(i1,i2,i3+2)+dyi2
              y_c(i1,i2,i3+3)=y_c(i1,i2,i3+3)+dyi3
           enddo
           icur=i3
        else
           icur=ibxy_c(1,i1,i2)
        endif

        do i3=icur,ibxy_c(2,i1,i2)
           dyi0=0.0_wp
           if(with_confpot) then
              do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_c(2,i1,i2))
                 dyi0=dyi0 + xz_c(t,i1,i2)*aeff0array(t-i3-0,i3) + &
                             2.d0*(xza_c(t,i1,i2)+yza_c(t,i1,i2))*aeff0_2array(t-i3-0,i3+0)
              enddo
           else
              do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_c(2,i1,i2))
                 dyi0=dyi0 + xz_c(t,i1,i2)*aeff0array(t-i3-0,i3)
              enddo
           end if
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi0
        enddo

        istart=max(ibxy_c(1,i1,i2),ibxy_f(1,i1,i2)-lupfil)
        iend=min(ibxy_c(2,i1,i2),ibxy_f(2,i1,i2)-lowfil)

        if (istart-iend.ge.4) then
           do i3=istart,iend-4,4
              dyi0=0.0_wp ; dyi1=0.0_wp ; dyi2=0.0_wp ; dyi3=0.0_wp
              if(with_confpot) then
                 do t=max(ibxy_f(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_f(2,i1,i2))
                    dyi0 = dyi0 + xz_f4(t,i1,i2)*beff0array(t-i3-0,i3+0) + &
                                  2.d0*(xzb_f(1,t,i1,i2)+yzb_f(1,t,i1,i2))*aeff0_2array(t-i3-0,i3+0) + &
                                  2.d0*(xza_f(2,t,i1,i2)+xzb_f(3,t,i1,i2)+yza_f(2,t,i1,i2)+yzb_f(3,t,i1,i2))&
                                      *beff0_2array(t-i3-0,i3+0)
                    dyi1 = dyi1 + xz_f4(t,i1,i2)*beff0array(t-i3-1,i3+1) + &
                                  2.d0*(xzb_f(1,t,i1,i2)+yzb_f(1,t,i1,i2))*aeff0_2array(t-i3-1,i3+1) + &
                                  2.d0*(xza_f(2,t,i1,i2)+xzb_f(3,t,i1,i2)+yza_f(2,t,i1,i2)+yzb_f(3,t,i1,i2))&
                                      *beff0_2array(t-i3-1,i3+1)
                    dyi2 = dyi2 + xz_f4(t,i1,i2)*beff0array(t-i3-2,i3+2) + &
                                  2.d0*(xzb_f(1,t,i1,i2)+yzb_f(1,t,i1,i2))*aeff0_2array(t-i3-2,i3+2) + &
                                  2.d0*(xza_f(2,t,i1,i2)+xzb_f(3,t,i1,i2)+yza_f(2,t,i1,i2)+yzb_f(3,t,i1,i2))&
                                      *beff0_2array(t-i3-2,i3+2)
                    dyi3 = dyi3 + xz_f4(t,i1,i2)*beff0array(t-i3-3,i3+3) + &
                                  2.d0*(xzb_f(1,t,i1,i2)+yzb_f(1,t,i1,i2))*aeff0_2array(t-i3-3,i3+3) + &
                                  2.d0*(xza_f(2,t,i1,i2)+xzb_f(3,t,i1,i2)+yza_f(2,t,i1,i2)+yzb_f(3,t,i1,i2))&
                                      *beff0_2array(t-i3-3,i3+3)
                 enddo
              else
                 do t=max(ibxy_f(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_f(2,i1,i2))
                    dyi0 = dyi0 + xz_f4(t,i1,i2)*beff0array(t-i3-0,i3+0)
                    dyi1 = dyi1 + xz_f4(t,i1,i2)*beff0array(t-i3-1,i3+1)
                    dyi2 = dyi2 + xz_f4(t,i1,i2)*beff0array(t-i3-2,i3+2)
                    dyi3 = dyi3 + xz_f4(t,i1,i2)*beff0array(t-i3-3,i3+3)
                 enddo
              end if
              y_c(i1,i2,i3+0)=y_c(i1,i2,i3+0)+dyi0
              y_c(i1,i2,i3+1)=y_c(i1,i2,i3+1)+dyi1
              y_c(i1,i2,i3+2)=y_c(i1,i2,i3+2)+dyi2
              y_c(i1,i2,i3+3)=y_c(i1,i2,i3+3)+dyi3
           enddo
           istart=i2
        endif

        do i3=istart,iend
           dyi0=0.0_wp ; tt0=0.0_wp
           if(with_confpot) then
              do t=max(ibxy_f(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_f(2,i1,i2))
                 dyi0=dyi0 + xz_f4(t,i1,i2)*beff0array(t-i3-0,i3) + &
                             2.d0*(xzb_f(1,t,i1,i2)+yzb_f(1,t,i1,i2))*aeff0_2array(t-i3-0,i3) + &
                             2.d0*(xza_f(2,t,i1,i2)+xzb_f(3,t,i1,i2)+yza_f(2,t,i1,i2)+yzb_f(3,t,i1,i2))&
                                 *beff0_2array(t-i3-0,i3)
              enddo
           else
              do t=max(ibxy_f(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_f(2,i1,i2))
                 dyi0=dyi0 + xz_f4(t,i1,i2)*beff0array(t-i3-0,i3)
              enddo
           end if
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi0
        enddo

         if (ibxy_f(2,i1,i2)-ibxy_f(1,i1,i2).ge.4) then
           do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)-4,4
              dyi0=0.0_wp ; dyi1=0.0_wp ; dyi2=0.0_wp ; dyi3=0.0_wp 
              tt10 = 0.d0 ; tt11 = 0.d0 ; tt12 = 0.d0 ; tt13 = 0.d0
              tt40 = 0.d0 ; tt41 = 0.d0 ; tt42 = 0.d0 ; tt43 = 0.d0
              tt50 = 0.d0 ; tt51 = 0.d0 ; tt52 = 0.d0 ; tt53 = 0.d0
              tt20 = 0.d0 ; tt21 = 0.d0 ; tt22 = 0.d0 ; tt23 = 0.d0
              tt60 = 0.d0 ; tt61 = 0.d0 ; tt62 = 0.d0 ; tt63 = 0.d0
              do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_c(2,i1,i2))
                 dyi0=dyi0 + xz_c(t,i1,i2)*ceff0array(t-i3-0,i3+0)
                 dyi1=dyi1 + xz_c(t,i1,i2)*ceff0array(t-i3-1,i3+1)
                 dyi2=dyi2 + xz_c(t,i1,i2)*ceff0array(t-i3-2,i3+2)
                 dyi3=dyi3 + xz_c(t,i1,i2)*ceff0array(t-i3-3,i3+3)
              end do

              if(with_confpot) then
                 do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_c(2,i1,i2))
                    tt10 = tt10 + 2.d0*xzc_c(t,i1,i2)*aeff0_2array(t-i3-0,i3+0)
                    tt11 = tt11 + 2.d0*xzc_c(t,i1,i2)*aeff0_2array(t-i3-1,i3+1)
                    tt12 = tt12 + 2.d0*xzc_c(t,i1,i2)*aeff0_2array(t-i3-2,i3+2)
                    tt13 = tt13 + 2.d0*xzc_c(t,i1,i2)*aeff0_2array(t-i3-3,i3+3)

                    tt40 = tt40 + 2.d0*xza_c(t,i1,i2)*ceff0_2array(t-i3-0,i3+0)
                    tt41 = tt41 + 2.d0*xza_c(t,i1,i2)*ceff0_2array(t-i3-1,i3+1)
                    tt42 = tt42 + 2.d0*xza_c(t,i1,i2)*ceff0_2array(t-i3-2,i3+2)
                    tt43 = tt43 + 2.d0*xza_c(t,i1,i2)*ceff0_2array(t-i3-3,i3+3)

                    tt50 = tt50 + 2.d0*xzc_c(t,i1,i2)*ceff0_2array(t-i3-0,i3+0)
                    tt51 = tt51 + 2.d0*xzc_c(t,i1,i2)*ceff0_2array(t-i3-1,i3+1)
                    tt52 = tt52 + 2.d0*xzc_c(t,i1,i2)*ceff0_2array(t-i3-2,i3+2)
                    tt53 = tt53 + 2.d0*xzc_c(t,i1,i2)*ceff0_2array(t-i3-3,i3+3)

                    tt20 = tt20 + 2.d0*yzc_c(t,i1,i2)*aeff0_2array(t-i3-0,i3+0)
                    tt21 = tt21 + 2.d0*yzc_c(t,i1,i2)*aeff0_2array(t-i3-1,i3+1)
                    tt22 = tt22 + 2.d0*yzc_c(t,i1,i2)*aeff0_2array(t-i3-2,i3+2)
                    tt23 = tt23 + 2.d0*yzc_c(t,i1,i2)*aeff0_2array(t-i3-3,i3+3)

                    tt40 = tt40 + 2.d0*yza_c(t,i1,i2)*ceff0_2array(t-i3-0,i3+0)
                    tt41 = tt41 + 2.d0*yza_c(t,i1,i2)*ceff0_2array(t-i3-1,i3+1)
                    tt42 = tt42 + 2.d0*yza_c(t,i1,i2)*ceff0_2array(t-i3-2,i3+2)
                    tt43 = tt43 + 2.d0*yza_c(t,i1,i2)*ceff0_2array(t-i3-3,i3+3)

                    tt60 = tt60 + 2.d0*yzc_c(t,i1,i2)*ceff0_2array(t-i3-0,i3+0)
                    tt61 = tt61 + 2.d0*yzc_c(t,i1,i2)*ceff0_2array(t-i3-1,i3+1)
                    tt62 = tt62 + 2.d0*yzc_c(t,i1,i2)*ceff0_2array(t-i3-2,i3+2)
                    tt63 = tt63 + 2.d0*yzc_c(t,i1,i2)*ceff0_2array(t-i3-3,i3+3)
                 enddo
                 y_f(1,i1,i2,i3+0) = y_f(1,i1,i2,i3+0) + tt10
                 y_f(1,i1,i2,i3+1) = y_f(1,i1,i2,i3+1) + tt11
                 y_f(1,i1,i2,i3+2) = y_f(1,i1,i2,i3+2) + tt12
                 y_f(1,i1,i2,i3+3) = y_f(1,i1,i2,i3+3) + tt13

                 y_f(5,i1,i2,i3+0) = y_f(5,i1,i2,i3+0) + tt50
                 y_f(5,i1,i2,i3+1) = y_f(5,i1,i2,i3+1) + tt51
                 y_f(5,i1,i2,i3+2) = y_f(5,i1,i2,i3+2) + tt52
                 y_f(5,i1,i2,i3+3) = y_f(5,i1,i2,i3+3) + tt53

                 y_f(2,i1,i2,i3+0) = y_f(2,i1,i2,i3+0) + tt20
                 y_f(2,i1,i2,i3+1) = y_f(2,i1,i2,i3+1) + tt21
                 y_f(2,i1,i2,i3+2) = y_f(2,i1,i2,i3+2) + tt22
                 y_f(2,i1,i2,i3+3) = y_f(2,i1,i2,i3+3) + tt23

                 y_f(6,i1,i2,i3+0) = y_f(6,i1,i2,i3+0) + tt60
                 y_f(6,i1,i2,i3+1) = y_f(6,i1,i2,i3+1) + tt61
                 y_f(6,i1,i2,i3+2) = y_f(6,i1,i2,i3+2) + tt62
                 y_f(6,i1,i2,i3+3) = y_f(6,i1,i2,i3+3) + tt63
                 end if
                 y_f(4,i1,i2,i3+0) = y_f(4,i1,i2,i3+0) + dyi0 + tt40
                 y_f(4,i1,i2,i3+1) = y_f(4,i1,i2,i3+1) + dyi1 + tt41
                 y_f(4,i1,i2,i3+2) = y_f(4,i1,i2,i3+2) + dyi2 + tt42
                 y_f(4,i1,i2,i3+3) = y_f(4,i1,i2,i3+3) + dyi3 + tt43
           enddo
           icur=i3
        else
           icur=ibxy_f(1,i1,i2)
        endif

        do i3=icur,ibxy_f(2,i1,i2)
           dyi0=0.0_wp ; tt10 = 0.d0 ; tt40 = 0.d0 ; tt50 = 0.d0 ; tt20 = 0.d0 ; tt60 = 0.d0
           do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_c(2,i1,i2))
              dyi0=dyi0 + xz_c(t,i1,i2)*ceff0array(t-i3-0,i3)
           end do

           if(with_confpot) then
              do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_c(2,i1,i2))
                 tt10 = tt10 + 2.d0*xzc_c(t,i1,i2)*aeff0_2array(t-i3-0,i3)

                 tt40 = tt40 + 2.d0*xza_c(t,i1,i2)*ceff0_2array(t-i3-0,i3)

                 tt50 = tt50 + 2.d0*xzc_c(t,i1,i2)*ceff0_2array(t-i3-0,i3)

                 tt20 = tt20 + 2.d0*yzc_c(t,i1,i2)*aeff0_2array(t-i3-0,i3)

                 tt40 = tt40 + 2.d0*yza_c(t,i1,i2)*ceff0_2array(t-i3-0,i3)

                 tt60 = tt60 + 2.d0*yzc_c(t,i1,i2)*ceff0_2array(t-i3-0,i3)
              enddo
           end if
           y_f(4,i1,i2,i3+0) = y_f(4,i1,i2,i3+0) + dyi0 + tt40
           y_f(1,i1,i2,i3+0) = y_f(1,i1,i2,i3+0) + tt10
           y_f(5,i1,i2,i3+0) = y_f(5,i1,i2,i3+0) + tt50
           y_f(2,i1,i2,i3+0) = y_f(2,i1,i2,i3+0) + tt20
           y_f(6,i1,i2,i3+0) = y_f(6,i1,i2,i3+0) + tt60
        enddo
     enddo
  enddo
  !$omp enddo


  ! wavelet part

  !$omp do 
  do i2=nfl2,nfu2
     do i1=nfl1,nfu1
        do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
           tt10 = 0.d0 ; tt20 = 0.d0 ; tt30 = 0.d0 ; tt40 = 0.d0 ; tt50 = 0.d0 ; tt60 = 0.d0 ; tt70 = 0.d0

           if(with_confpot) then
              do l=max(nfl3-i3,lowfil),min(lupfil,nfu3-i3)
                 tt10 = tt10 + xz_f(1,i3+l,i1,i2)*aeff0array(l,i3) + xz_f(5,i3+l,i1,i2)*beff0array(l,i3) + &
                               2.d0*                    xze_f(1,i3+l,i1,i2) *aeff0_2array(l,i3) + &
                               2.d0*(xzc_f(2,i3+l,i1,i2)+xze_f(3,i3+l,i1,i2))*beff0_2array(l,i3) + &
                               2.d0*(yza_f(1,i3+l,i1,i2)+yzb_f(2,i3+l,i1,i2))*aeff0_2array(l,i3) + &
                               2.d0*(yza_f(3,i3+l,i1,i2)+yzb_f(4,i3+l,i1,i2))*beff0_2array(l,i3)

                 tt20 = tt20 + xz_f(2,i3+l,i1,i2)*aeff0array(l,i3) + xz_f(6,i3+l,i1,i2)*beff0array(l,i3) + &
                               2.d0*(xza_f(1,i3+l,i1,i2)+xzb_f(2,i3+l,i1,i2))*aeff0_2array(l,i3) + &
                               2.d0*(xza_f(3,i3+l,i1,i2)+xzb_f(4,i3+l,i1,i2))*beff0_2array(l,i3) + &
                               2.d0*                    yze_f(1,i3+l,i1,i2) *aeff0_2array(l,i3) + &
                               2.d0*(yzc_f(2,i3+l,i1,i2)+yze_f(3,i3+l,i1,i2))*beff0_2array(l,i3)

                 tt30 = tt30 + xz_f(3,i3+l,i1,i2)*aeff0array(l,i3) + xz_f(7,i3+l,i1,i2)*beff0array(l,i3) + &
                               2.d0*(xzc_f(1,i3+l,i1,i2)+xze_f(2,i3+l,i1,i2))*aeff0_2array(l,i3) + &
                               2.d0*(xzc_f(3,i3+l,i1,i2)+xze_f(4,i3+l,i1,i2))*beff0_2array(l,i3) + &
                               2.d0*(yzc_f(1,i3+l,i1,i2)+yze_f(2,i3+l,i1,i2))*aeff0_2array(l,i3) + &
                               2.d0*(yzc_f(3,i3+l,i1,i2)+yze_f(4,i3+l,i1,i2))*beff0_2array(l,i3)

                 tt40 = tt40 + xz_f(4,i3+l,i1,i2)*eeff0array(l,i3)                                      + &
                               2.d0*                    xzb_f(1,i3+l,i1,i2) *ceff0_2array(l,i3) + &
                               2.d0*(xza_f(2,i3+l,i1,i2)+xzb_f(3,i3+l,i1,i2))*eeff0_2array(l,i3) + &
                               2.d0*                    yzb_f(1,i3+l,i1,i2) *ceff0_2array(l,i3) + &
                               2.d0*(yza_f(2,i3+l,i1,i2)+yzb_f(3,i3+l,i1,i2))*eeff0_2array(l,i3)

                 tt50 = tt50 + xz_f(1,i3+l,i1,i2)*ceff0array(l,i3) + xz_f(5,i3+l,i1,i2)*eeff0array(l,i3) + &
                               2.d0*                    xze_f(1,i3+l,i1,i2) *ceff0_2array(l,i3) + &
                               2.d0*(xzc_f(2,i3+l,i1,i2)+xze_f(3,i3+l,i1,i2))*eeff0_2array(l,i3) + &
                               2.d0*(yza_f(1,i3+l,i1,i2)+yzb_f(2,i3+l,i1,i2))*ceff0_2array(l,i3) + &
                               2.d0*(yza_f(3,i3+l,i1,i2)+yzb_f(4,i3+l,i1,i2))*eeff0_2array(l,i3)

                 tt60 = tt60 + xz_f(2,i3+l,i1,i2)*ceff0array(l,i3) + xz_f(6,i3+l,i1,i2)*eeff0array(l,i3) + &
                               2.d0*(xza_f(1,i3+l,i1,i2)+xzb_f(2,i3+l,i1,i2))*ceff0_2array(l,i3) + &
                               2.d0*(xza_f(3,i3+l,i1,i2)+xzb_f(4,i3+l,i1,i2))*eeff0_2array(l,i3) + &
                               2.d0*                    yze_f(1,i3+l,i1,i2) *ceff0_2array(l,i3) + &
                               2.d0*(yzc_f(2,i3+l,i1,i2)+yze_f(3,i3+l,i1,i2))*eeff0_2array(l,i3)

                 tt70 = tt70 + xz_f(3,i3+l,i1,i2)*ceff0array(l,i3) + xz_f(7,i3+l,i1,i2)*eeff0array(l,i3) + &
                               2.d0*(xzc_f(1,i3+l,i1,i2)+xze_f(2,i3+l,i1,i2))*ceff0_2array(l,i3) + &
                               2.d0*(xzc_f(3,i3+l,i1,i2)+xze_f(4,i3+l,i1,i2))*eeff0_2array(l,i3) + &
                               2.d0*(yzc_f(1,i3+l,i1,i2)+yze_f(2,i3+l,i1,i2))*ceff0_2array(l,i3) + &
                               2.d0*(yzc_f(3,i3+l,i1,i2)+yze_f(4,i3+l,i1,i2))*eeff0_2array(l,i3)
              enddo
           else
              do l=max(nfl3-i3,lowfil),min(lupfil,nfu3-i3)
                 tt10 = tt10 + xz_f(1,i3+l,i1,i2)*aeff0array(l,i3) + xz_f(5,i3+l,i1,i2)*beff0array(l,i3)

                 tt20 = tt20 + xz_f(2,i3+l,i1,i2)*aeff0array(l,i3) + xz_f(6,i3+l,i1,i2)*beff0array(l,i3)

                 tt30 = tt30 + xz_f(3,i3+l,i1,i2)*aeff0array(l,i3) + xz_f(7,i3+l,i1,i2)*beff0array(l,i3)

                 tt40 = tt40 + xz_f(4,i3+l,i1,i2)*eeff0array(l,i3) 

                 tt50 = tt50 + xz_f(1,i3+l,i1,i2)*ceff0array(l,i3) + xz_f(5,i3+l,i1,i2)*eeff0array(l,i3)

                 tt60 = tt60 + xz_f(2,i3+l,i1,i2)*ceff0array(l,i3) + xz_f(6,i3+l,i1,i2)*eeff0array(l,i3)

                 tt70 = tt70 + xz_f(3,i3+l,i1,i2)*ceff0array(l,i3) + xz_f(7,i3+l,i1,i2)*eeff0array(l,i3)
              enddo
           end if
           y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+tt40
           y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+tt20
           y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+tt10
           y_f(6,i1,i2,i3)=y_f(6,i1,i2,i3)+tt60
           y_f(5,i1,i2,i3)=y_f(5,i1,i2,i3)+tt50
           y_f(3,i1,i2,i3)=y_f(3,i1,i2,i3)+tt30
           y_f(7,i1,i2,i3)=y_f(7,i1,i2,i3)+tt70
        enddo
     enddo
  enddo
  !$omp enddo
  !$omp end parallel

  call timing(iproc,'convolQuartic ','OF')


END SUBROUTINE ConvolQuartic4




subroutine createDerivativeBasis(n1,n2,n3, &
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
     hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,&
     w_c, w_f, w_f1, w_f2, w_f3, x_c, x_f, y_c, y_f, z_c, z_f)
     
  use module_base
  use filterModule
  implicit none

  ! Calling arguments
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  real(gp), intent(in) :: hgrid
  integer, dimension(2,0:n2,0:n3), intent(in) :: ibyz_c,ibyz_f
  integer, dimension(2,0:n1,0:n3), intent(in) :: ibxz_c,ibxz_f
  integer, dimension(2,0:n1,0:n2), intent(in) :: ibxy_c,ibxy_f
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: w_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in) :: w_f
  real(wp), dimension(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in) :: w_f1
  real(wp), dimension(nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), intent(in) :: w_f2
  real(wp), dimension(nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), intent(in) :: w_f3
  real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: x_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(out) :: x_f
  real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: y_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(out) :: y_f
  real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: z_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(out) :: z_f
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  integer :: i,t,i1,i2,i3
  integer :: icur,istart,iend,l
  real(wp) :: dyi,dyi0,dyi1,dyi2,dyi3,t112,t121,t122,t212,t221,t222,t211
  real(wp), dimension(-3+lowfil:lupfil+3) :: ad1_ext
  real(wp), dimension(-3+lowfil:lupfil+3) :: bd1_ext
  real(wp), dimension(-3+lowfil:lupfil+3) :: cd1_ext


! Copy the filters to the 'extended filters', i.e. add some zeros.
! This seems to be required since we use loop unrolling.
ad1_ext=0.d0
bd1_ext=0.d0
cd1_ext=0.d0

do i=lowfil,lupfil
    ad1_ext(i)=ad1(i)
    bd1_ext(i)=bd1(i)
    cd1_ext(i)=cd1(i)
end do


!$omp parallel default(private) &
!$omp shared(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3) &
!$omp shared(ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,w_c,w_f,y_c,y_f,x_c,x_f,z_c,z_f)& 
!$omp shared(w_f1,w_f2,w_f3,ad1_ext,bd1_ext,cd1_ext)

  ! x direction
  !$omp do
  do i3=0,n3
     do i2=0,n2
        if (ibyz_c(2,i2,i3)-ibyz_c(1,i2,i3).ge.4) then
           do i1=ibyz_c(1,i2,i3),ibyz_c(2,i2,i3)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_c(2,i2,i3))
                 dyi0=dyi0 + w_c(t,i2,i3)*ad1_ext(t-i1-0)
                 dyi1=dyi1 + w_c(t,i2,i3)*ad1_ext(t-i1-1)
                 dyi2=dyi2 + w_c(t,i2,i3)*ad1_ext(t-i1-2)
                 dyi3=dyi3 + w_c(t,i2,i3)*ad1_ext(t-i1-3)
              enddo
              x_c(i1+0,i2,i3)=dyi0
              x_c(i1+1,i2,i3)=dyi1
              x_c(i1+2,i2,i3)=dyi2
              x_c(i1+3,i2,i3)=dyi3
           enddo
           icur=i1
        else
           icur=ibyz_c(1,i2,i3)
        endif

        do i1=icur,ibyz_c(2,i2,i3)
           dyi=0.0_wp 
           !! Get the effective a-filters for the x dimension
           do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_c(2,i2,i3))
              dyi=dyi + w_c(t,i2,i3)*ad1_ext(t-i1)
           enddo
           x_c(i1,i2,i3)=dyi
        enddo

        istart=max(ibyz_c(1,i2,i3),ibyz_f(1,i2,i3)-lupfil)
        iend=min(ibyz_c(2,i2,i3),ibyz_f(2,i2,i3)-lowfil)

        if (istart-iend.ge.4) then
           do i1=istart,iend-4,4
              dyi0=0.0_wp
              dyi1=0.0_wp
              dyi2=0.0_wp
              dyi3=0.0_wp
              do t=max(ibyz_f(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_f(2,i2,i3))
                 dyi0=dyi0 + w_f1(t,i2,i3)*bd1_ext(t-i1-0)
                 dyi1=dyi1 + w_f1(t,i2,i3)*bd1_ext(t-i1-1)
                 dyi2=dyi2 + w_f1(t,i2,i3)*bd1_ext(t-i1-2)
                 dyi3=dyi3 + w_f1(t,i2,i3)*bd1_ext(t-i1-3)
              enddo
              x_c(i1+0,i2,i3)=x_c(i1+0,i2,i3)+dyi0
              x_c(i1+1,i2,i3)=x_c(i1+1,i2,i3)+dyi1
              x_c(i1+2,i2,i3)=x_c(i1+2,i2,i3)+dyi2
              x_c(i1+3,i2,i3)=x_c(i1+3,i2,i3)+dyi3
           enddo
           istart=i1
        endif

        do i1=istart,iend
           dyi=0.0_wp
           do t=max(ibyz_f(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_f(2,i2,i3))
              dyi=dyi + w_f1(t,i2,i3)*bd1_ext(t-i1)
           enddo
           x_c(i1,i2,i3)=x_c(i1,i2,i3)+dyi
        enddo

         if (ibyz_c(2,i2,i3)-ibyz_c(1,i2,i3).ge.4) then
           do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_c(2,i2,i3))
                 dyi0=dyi0 + w_c(t,i2,i3)*cd1_ext(t-i1-0)
                 dyi1=dyi1 + w_c(t,i2,i3)*cd1_ext(t-i1-1)
                 dyi2=dyi2 + w_c(t,i2,i3)*cd1_ext(t-i1-2)
                 dyi3=dyi3 + w_c(t,i2,i3)*cd1_ext(t-i1-3)
              enddo
              x_f(1,i1+0,i2,i3)=dyi0
              x_f(1,i1+1,i2,i3)=dyi1
              x_f(1,i1+2,i2,i3)=dyi2
              x_f(1,i1+3,i2,i3)=dyi3
           enddo
           icur=i1
        else
           icur=ibyz_f(1,i2,i3)
        endif
        do i1=icur,ibyz_f(2,i2,i3)
           dyi=0.0_wp 
           do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_c(2,i2,i3))
              dyi=dyi + w_c(t,i2,i3)*cd1_ext(t-i1)
           enddo
           x_f(1,i1,i2,i3)=dyi
        enddo
     enddo
  enddo
  !$omp enddo
  

  ! y direction
  !$omp do
  do i3=0,n3
     do i1=0,n1
        if (ibxz_c(2,i1,i3)-ibxz_c(1,i1,i3).ge.4) then
           do i2=ibxz_c(1,i1,i3),ibxz_c(2,i1,i3)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_c(2,i1,i3))
                 dyi0=dyi0 + w_c(i1,t,i3)*ad1_ext(t-i2-0)
                 dyi1=dyi1 + w_c(i1,t,i3)*ad1_ext(t-i2-1)
                 dyi2=dyi2 + w_c(i1,t,i3)*ad1_ext(t-i2-2)
                 dyi3=dyi3 + w_c(i1,t,i3)*ad1_ext(t-i2-3)
              enddo
              y_c(i1,i2+0,i3)=dyi0
              y_c(i1,i2+1,i3)=dyi1
              y_c(i1,i2+2,i3)=dyi2
              y_c(i1,i2+3,i3)=dyi3
           enddo
           icur=i2
        else
           icur=ibxz_c(1,i1,i3)
        endif

        do i2=icur,ibxz_c(2,i1,i3)
           dyi=0.0_wp 
           do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
              dyi=dyi + w_c(i1,t,i3)*ad1_ext(t-i2)
           enddo
           y_c(i1,i2,i3)=dyi
        enddo
        istart=max(ibxz_c(1,i1,i3),ibxz_f(1,i1,i3)-lupfil)
        iend= min(ibxz_c(2,i1,i3),ibxz_f(2,i1,i3)-lowfil)

        if (istart-iend.ge.4) then
           do i2=istart,iend-4,4
              dyi0=0.0_wp
              dyi1=0.0_wp
              dyi2=0.0_wp
              dyi3=0.0_wp
              do t=max(ibxz_f(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_f(2,i1,i3))
                 dyi0=dyi0 + w_f2(t,i1,i3)*bd1_ext(t-i2-0)
                 dyi1=dyi1 + w_f2(t,i1,i3)*bd1_ext(t-i2-1)
                 dyi2=dyi2 + w_f2(t,i1,i3)*bd1_ext(t-i2-2)
                 dyi3=dyi3 + w_f2(t,i1,i3)*bd1_ext(t-i2-3)
              enddo
              y_c(i1,i2+0,i3)=y_c(i1,i2+0,i3)+dyi0
              y_c(i1,i2+1,i3)=y_c(i1,i2+1,i3)+dyi1
              y_c(i1,i2+2,i3)=y_c(i1,i2+2,i3)+dyi2
              y_c(i1,i2+3,i3)=y_c(i1,i2+3,i3)+dyi3
           enddo
           istart=i2
        endif

        do i2=istart,iend
           dyi=0.0_wp
           do t=max(ibxz_f(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_f(2,i1,i3))
              dyi=dyi + w_f2(t,i1,i3)*bd1_ext(t-i2)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
        enddo

         if (ibxz_f(2,i1,i3)-ibxz_f(1,i1,i3).ge.4) then
           do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_c(2,i1,i3))
                 dyi0=dyi0 + w_c(i1,t,i3)*cd1_ext(t-i2-0)
                 dyi1=dyi1 + w_c(i1,t,i3)*cd1_ext(t-i2-1)
                 dyi2=dyi2 + w_c(i1,t,i3)*cd1_ext(t-i2-2)
                 dyi3=dyi3 + w_c(i1,t,i3)*cd1_ext(t-i2-3)
              enddo
              y_f(2,i1,i2+0,i3)=dyi0
              y_f(2,i1,i2+1,i3)=dyi1
              y_f(2,i1,i2+2,i3)=dyi2
              y_f(2,i1,i2+3,i3)=dyi3
           enddo
           icur=i2
        else
           icur=ibxz_f(1,i1,i3)
        endif

        do i2=icur,ibxz_f(2,i1,i3)
           dyi=0.0_wp 
           do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
              dyi=dyi + w_c(i1,t,i3)*cd1_ext(t-i2)
           enddo
           y_f(2,i1,i2,i3)=dyi
        enddo
     enddo
  enddo
  !$omp enddo


  ! z direction
  !$omp do
  do i2=0,n2
     do i1=0,n1
        if (ibxy_c(2,i1,i2)-ibxy_c(1,i1,i2).ge.4) then
           do i3=ibxy_c(1,i1,i2),ibxy_c(2,i1,i2)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_c(2,i1,i2))
                 dyi0=dyi0 + w_c(i1,i2,t)*ad1_ext(t-i3-0)
                 dyi1=dyi1 + w_c(i1,i2,t)*ad1_ext(t-i3-1)
                 dyi2=dyi2 + w_c(i1,i2,t)*ad1_ext(t-i3-2)
                 dyi3=dyi3 + w_c(i1,i2,t)*ad1_ext(t-i3-3)
              enddo
              z_c(i1,i2,i3+0)=dyi0
              z_c(i1,i2,i3+1)=dyi1
              z_c(i1,i2,i3+2)=dyi2
              z_c(i1,i2,i3+3)=dyi3
           enddo
           icur=i3
        else
           icur=ibxy_c(1,i1,i2)
        endif

        do i3=icur,ibxy_c(2,i1,i2)
           dyi=0.0_wp
           do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_c(2,i1,i2))
              dyi=dyi + w_c(i1,i2,t)*ad1_ext(t-i3)
           enddo
           z_c(i1,i2,i3)=dyi
        enddo
        istart=max(ibxy_c(1,i1,i2),ibxy_f(1,i1,i2)-lupfil)
        iend=min(ibxy_c(2,i1,i2),ibxy_f(2,i1,i2)-lowfil)

        if (istart-iend.ge.4) then
           do i3=istart,iend-4,4
              dyi0=0.0_wp
              dyi1=0.0_wp
              dyi2=0.0_wp
              dyi3=0.0_wp
              do t=max(ibxy_f(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_f(2,i1,i2))
                 dyi0=dyi0 + w_f3(t,i1,i2)*bd1_ext(t-i3-0)
                 dyi1=dyi1 + w_f3(t,i1,i2)*bd1_ext(t-i3-1)
                 dyi2=dyi2 + w_f3(t,i1,i2)*bd1_ext(t-i3-2)
                 dyi3=dyi3 + w_f3(t,i1,i2)*bd1_ext(t-i3-3)
              enddo
              z_c(i1,i2,i3+0)=z_c(i1,i2,i3+0)+dyi0
              z_c(i1,i2,i3+1)=z_c(i1,i2,i3+1)+dyi1
              z_c(i1,i2,i3+2)=z_c(i1,i2,i3+2)+dyi2
              z_c(i1,i2,i3+3)=z_c(i1,i2,i3+3)+dyi3
           enddo
           istart=i2
        endif

        do i3=istart,iend
           dyi=0.0_wp
           do t=max(ibxy_f(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_f(2,i1,i2))
              dyi=dyi + w_f3(t,i1,i2)*bd1_ext(t-i3)
           enddo
           z_c(i1,i2,i3)=z_c(i1,i2,i3)+dyi
        enddo

         if (ibxy_f(2,i1,i2)-ibxy_f(1,i1,i2).ge.4) then
           do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_c(2,i1,i2))
                 dyi0=dyi0 + w_c(i1,i2,t)*cd1_ext(t-i3-0)
                 dyi1=dyi1 + w_c(i1,i2,t)*cd1_ext(t-i3-1)
                 dyi2=dyi2 + w_c(i1,i2,t)*cd1_ext(t-i3-2)
                 dyi3=dyi3 + w_c(i1,i2,t)*cd1_ext(t-i3-3)
              enddo
              z_f(4,i1,i2,i3+0)=dyi0
              z_f(4,i1,i2,i3+1)=dyi1
              z_f(4,i1,i2,i3+2)=dyi2
              z_f(4,i1,i2,i3+3)=dyi3
           enddo
           icur=i3
        else
           icur=ibxy_f(1,i1,i2)
        endif

        do i3=icur,ibxy_f(2,i1,i2)
           dyi=0.0_wp 
           do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_c(2,i1,i2))
              dyi=dyi + w_c(i1,i2,t)*cd1_ext(t-i3)
           enddo
           z_f(4,i1,i2,i3)=dyi
        enddo
     enddo
  enddo
  !$omp enddo

  ! wavelet part

  ! x direction
  !$omp do
  do i3=nfl3,nfu3
     do i2=nfl2,nfu2
        do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
           t112=0.0_wp;t121=0.0_wp;t122=0.0_wp;t212=0.0_wp;t221=0.0_wp;t222=0.0_wp;t211=0.0_wp 
           do l=max(nfl1-i1,lowfil),min(lupfil,nfu1-i1)
              t121=t121 + w_f(2,i1+l,i2,i3)*ad1_ext(l) + w_f(3,i1+l,i2,i3)*bd1_ext(l)
              t221=t221 + w_f(2,i1+l,i2,i3)*cd1_ext(l) + w_f(3,i1+l,i2,i3)*ed1(l)
           enddo
           x_f(4,i1,i2,i3)=t112
           x_f(2,i1,i2,i3)=t121
           x_f(1,i1,i2,i3)=x_f(1,i1,i2,i3)+t211
           x_f(6,i1,i2,i3)=t122
           x_f(5,i1,i2,i3)=t212
           x_f(3,i1,i2,i3)=t221
           x_f(7,i1,i2,i3)=t222
        enddo
     enddo
  enddo
  !$omp enddo


  ! y direction
  !$omp do
  do i3=nfl3,nfu3
     do i1=nfl1,nfu1
        do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
           t112=0.0_wp;t121=0.0_wp;t122=0.0_wp;t212=0.0_wp;t221=0.0_wp;t222=0.0_wp;t211=0.0_wp 
           do l=max(nfl2-i2,lowfil),min(lupfil,nfu2-i2)
              t112=t112 + w_f(4,i1,i2+l,i3)*ad1_ext(l) + w_f(6,i1,i2+l,i3)*bd1_ext(l)
              t211=t211 + w_f(1,i1,i2+l,i3)*ad1_ext(l) + w_f(3,i1,i2+l,i3)*bd1_ext(l)
              t122=t122 + w_f(4,i1,i2+l,i3)*cd1_ext(l) + w_f(6,i1,i2+l,i3)*ed1(l)
              t212=t212 + w_f(5,i1,i2+l,i3)*ad1_ext(l) + w_f(7,i1,i2+l,i3)*bd1_ext(l)
              t221=t221 + w_f(1,i1,i2+l,i3)*cd1_ext(l) + w_f(3,i1,i2+l,i3)*ed1(l)
              t222=t222 + w_f(5,i1,i2+l,i3)*cd1_ext(l) + w_f(7,i1,i2+l,i3)*ed1(l)
              t121=t121 + w_f(2,i1,i2+l,i3)*ed1(l)
           enddo
           y_f(4,i1,i2,i3)=t112
           y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+t121
           y_f(1,i1,i2,i3)=t211
           y_f(6,i1,i2,i3)=t122
           y_f(5,i1,i2,i3)=t212
           y_f(3,i1,i2,i3)=t221
           y_f(7,i1,i2,i3)=t222
        enddo
     enddo
  enddo
  !$omp enddo


  ! z direction
  !$omp do
  do i2=nfl2,nfu2
     do i1=nfl1,nfu1
        do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
           t112=0.0_wp;t121=0.0_wp;t122=0.0_wp;t212=0.0_wp;t221=0.0_wp;t222=0.0_wp;t211=0.0_wp 
           do l=max(nfl3-i3,lowfil),min(lupfil,nfu3-i3)
              t121=t121 + w_f(2,i1,i2,i3+l)*ad1_ext(l) + w_f(6,i1,i2,i3+l)*bd1_ext(l)
              t211=t211 + w_f(1,i1,i2,i3+l)*ad1_ext(l) + w_f(5,i1,i2,i3+l)*bd1_ext(l)
              t122=t122 + w_f(2,i1,i2,i3+l)*cd1_ext(l) + w_f(6,i1,i2,i3+l)*ed1(l)
              t212=t212 + w_f(1,i1,i2,i3+l)*cd1_ext(l) + w_f(5,i1,i2,i3+l)*ed1(l)
              t221=t221 + w_f(3,i1,i2,i3+l)*ad1_ext(l) + w_f(7,i1,i2,i3+l)*bd1_ext(l)
              t222=t222 + w_f(3,i1,i2,i3+l)*cd1_ext(l) + w_f(7,i1,i2,i3+l)*ed1(l)
              t112=t112 + w_f(4,i1,i2,i3+l)*ed1(l)
           enddo
           z_f(4,i1,i2,i3)=z_f(4,i1,i2,i3)+t112
           z_f(2,i1,i2,i3)=t121
           z_f(1,i1,i2,i3)=t211
           z_f(6,i1,i2,i3)=t122
           z_f(5,i1,i2,i3)=t212
           z_f(3,i1,i2,i3)=t221
           z_f(7,i1,i2,i3)=t222

        enddo
     enddo
  enddo
  !$omp enddo

  !$omp end parallel

END SUBROUTINE createDerivativeBasis
