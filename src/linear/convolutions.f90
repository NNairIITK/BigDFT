subroutine ConvolQuartic4(iproc, nproc, n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3,  &
           hgrid, offsetx, offsety, offsetz, ibyz_c, ibxz_c, ibxy_c, ibyz_f, ibxz_f, ibxy_f, &
           rxyzConf, potentialPrefac, with_kinetic, cprecr, &
           xx_c, xx_f1, xx_f, xy_c, xy_f2, xy_f,  xz_c, xz_f4, xz_f, &
           y_c, y_f)

  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc, n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3, offsetx, offsety, offsetz
  real(gp),intent(in):: hgrid, potentialPrefac, cprecr
  logical,intent(in):: with_kinetic
  real(8),dimension(3):: rxyzConf
  integer,dimension(2,0:n2,0:n3), intent(in) :: ibyz_c,ibyz_f
  integer,dimension(2,0:n1,0:n3), intent(in) :: ibxz_c,ibxz_f
  integer,dimension(2,0:n1,0:n2), intent(in) :: ibxy_c,ibxy_f
  real(wp),dimension(0:n1,0:n2,0:n3),intent(in):: xx_c
  real(wp),dimension(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3),intent(in):: xx_f1
  real(wp),dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3),intent(in):: xx_f
  real(wp),dimension(0:n2,0:n1,0:n3),intent(in):: xy_c
  real(wp),dimension(nfl2:nfu2,nfl1:nfu1,nfl3:nfu3),intent(in):: xy_f2
  real(wp),dimension(7,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3),intent(in):: xy_f
  real(wp),dimension(0:n3,0:n1,0:n2),intent(in):: xz_c
  real(wp),dimension(nfl3:nfu3,nfl1:nfu1,nfl2:nfu2),intent(in):: xz_f4
  real(wp),dimension(7,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2),intent(in):: xz_f
  real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: y_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(out) :: y_f

  ! Local variables
  integer,parameter:: lowfil=-14,lupfil=14
  integer:: i,t,i1,i2,i3, icur,istart,iend,l, istat, iall
  real(wp):: dyi,dyi0,dyi1,dyi2,dyi3,t112,t121,t122,t212,t221,t222,t211
  real(wp):: tt112, tt121, tt122, tt212, tt221, tt222, tt211, tt0
  real(wp),dimension(-3+lowfil:lupfil+3):: a, aeff0, aeff1, aeff2, aeff3
  real(wp),dimension(-3+lowfil:lupfil+3):: b, beff0, beff1, beff2, beff3
  real(wp),dimension(-3+lowfil:lupfil+3):: c, ceff0, ceff1, ceff2, ceff3
  real(wp),dimension(lowfil:lupfil):: e, eeff0, eeff1, eeff2, eeff3
  real(wp),dimension(-3+lowfil:lupfil+3):: aeff0_2, aeff1_2, aeff2_2, aeff3_2
  real(wp),dimension(-3+lowfil:lupfil+3):: beff0_2, beff1_2, beff2_2, beff3_2
  real(wp),dimension(-3+lowfil:lupfil+3):: ceff0_2, ceff1_2, ceff2_2, ceff3_2
  real(wp),dimension(lowfil:lupfil):: eeff0_2, eeff1_2, eeff2_2, eeff3_2
  real(wp),dimension(:,:),allocatable:: aeff0array, beff0array, ceff0array, eeff0array
  real(wp),dimension(:,:),allocatable:: aeff0_2array, beff0_2array, ceff0_2array, eeff0_2array
  real(wp),dimension(:,:),allocatable:: aeff0_2auxarray, beff0_2auxarray, ceff0_2auxarray, eeff0_2auxarray
  real(wp),dimension(:,:,:),allocatable:: xya_c, xyb_c, xyc_c, xye_c, xza_c, xzb_c, xzc_c, xze_c, yza_c, yzb_c, yzc_c, yze_c
  real(wp),dimension(:,:,:,:),allocatable:: xya_f, xyb_f, xyc_f, xye_f
  real(wp),dimension(:,:,:,:),allocatable:: xza_f, xzb_f, xzc_f, xze_f
  real(wp),dimension(:,:,:,:),allocatable:: yza_f, yzb_f, yzc_f, yze_f
  real(8):: x0, y0, z0
  real(8):: x1, y1, z1
  real(8):: x2, y2, z2
  real(8):: x3, y3, z3
  real(8):: tt00, tt01, tt02, tt03
  real(8):: tt10, tt11, tt12, tt13
  real(8):: tt20, tt21, tt22, tt23
  real(8):: tt30, tt31, tt32, tt33
  real(8):: tt40, tt41, tt42, tt43
  real(8):: tt50, tt51, tt52, tt53
  real(8):: tt60, tt61, tt62, tt63
  real(8):: tt70, tt71, tt72, tt73
  real(8):: tt0a0, tt0a1, tt0a2, tt0a3
  real(8):: tt0b0, tt0b1, tt0b2, tt0b3
  real(8):: tt0c0, tt0c1, tt0c2, tt0c3
  real(8):: tt0e0, tt0e1, tt0e2, tt0e3
  real(8):: tt1a0, tt1a1, tt1a2, tt1a3
  real(8):: tt1b0, tt1b1, tt1b2, tt1b3
  real(8):: tt1c0, tt1c1, tt1c2, tt1c3
  real(8):: tt1e0, tt1e1, tt1e2, tt1e3
  real(8):: tt2a0, tt2a1, tt2a2, tt2a3
  real(8):: tt2b0, tt2b1, tt2b2, tt2b3
  real(8):: tt2c0, tt2c1, tt2c2, tt2c3
  real(8):: tt2e0, tt2e1, tt2e2, tt2e3
  real(8):: tt3a0, tt3a1, tt3a2, tt3a3
  real(8):: tt3b0, tt3b1, tt3b2, tt3b3
  real(8):: tt3c0, tt3c1, tt3c2, tt3c3
  real(8):: tt3e0, tt3e1, tt3e2, tt3e3
  real(8):: tt4a0, tt4a1, tt4a2, tt4a3
  real(8):: tt4b0, tt4b1, tt4b2, tt4b3
  real(8):: tt4c0, tt4c1, tt4c2, tt4c3
  real(8):: tt4e0, tt4e1, tt4e2, tt4e3
  real(8):: tt5a0, tt5a1, tt5a2, tt5a3
  real(8):: tt5b0, tt5b1, tt5b2, tt5b3
  real(8):: tt5c0, tt5c1, tt5c2, tt5c3
  real(8):: tt5e0, tt5e1, tt5e2, tt5e3
  real(8):: tt6a0, tt6a1, tt6a2, tt6a3
  real(8):: tt6b0, tt6b1, tt6b2, tt6b3
  real(8):: tt6c0, tt6c1, tt6c2, tt6c3
  real(8):: tt6e0, tt6e1, tt6e2, tt6e3
  real(8):: tt7a0, tt7a1, tt7a2, tt7a3
  real(8):: tt7b0, tt7b1, tt7b2, tt7b3
  real(8):: tt7c0, tt7c1, tt7c2, tt7c3
  real(8):: tt7e0, tt7e1, tt7e2, tt7e3
  character(len=*),parameter:: subname='ConvolQuartic4'


call timing(iproc,'convolQuartic ','ON')

! Allocate all arrays

i=max(n1,n2,n3)
allocate(aeff0array(-3+lowfil:lupfil+3,0:i), stat=istat)
call memocc(istat, aeff0array, 'aeff0array', subname)
allocate(beff0array(-3+lowfil:lupfil+3,0:i), stat=istat)
call memocc(istat, beff0array, 'beff0array', subname)
allocate(ceff0array(-3+lowfil:lupfil+3,0:i), stat=istat)
call memocc(istat, ceff0array, 'ceff0array', subname)
allocate(eeff0array(lowfil:lupfil,0:i), stat=istat)
call memocc(istat, eeff0array, 'eeff0array', subname)
aeff0array=0.d0
beff0array=0.d0
ceff0array=0.d0
eeff0array=0.d0

allocate(aeff0_2array(-3+lowfil:lupfil+3,0:i), stat=istat)
call memocc(istat, aeff0_2array, 'aeff0_2array', subname)
allocate(beff0_2array(-3+lowfil:lupfil+3,0:i), stat=istat)
call memocc(istat, beff0_2array, 'beff0_2array', subname)
allocate(ceff0_2array(-3+lowfil:lupfil+3,0:i), stat=istat)
call memocc(istat, ceff0_2array, 'ceff0_2array', subname)
allocate(eeff0_2array(lowfil:lupfil,0:i), stat=istat)
call memocc(istat, eeff0_2array, 'eeff0_2array', subname)
aeff0_2array=0.d0
beff0_2array=0.d0
ceff0_2array=0.d0
eeff0_2array=0.d0

allocate(aeff0_2auxarray(-3+lowfil:lupfil+3,0:i), stat=istat)
call memocc(istat, aeff0_2auxarray, 'aeff0_2auxarray', subname)
allocate(beff0_2auxarray(-3+lowfil:lupfil+3,0:i), stat=istat)
call memocc(istat, beff0_2auxarray, 'beff0_2auxarray', subname)
allocate(ceff0_2auxarray(-3+lowfil:lupfil+3,0:i), stat=istat)
call memocc(istat, ceff0_2auxarray, 'ceff0_2auxarray', subname)
allocate(eeff0_2auxarray(-3+lowfil:lupfil+3,0:i), stat=istat)
call memocc(istat, eeff0_2auxarray, 'eeff0_2auxarray', subname)
aeff0_2auxarray=0.d0
beff0_2auxarray=0.d0
ceff0_2auxarray=0.d0
eeff0_2auxarray=0.d0

allocate(xya_c(0:n2,0:n1,0:n3), stat=istat)
call memocc(istat, xya_c, 'xya_c', subname)
allocate(xyb_c(0:n2,0:n1,0:n3), stat=istat)
call memocc(istat, xyb_c, 'xyb_c', subname)
allocate(xyc_c(0:n2,0:n1,0:n3), stat=istat)
call memocc(istat, xyc_c, 'xyc_c', subname)
allocate(xye_c(0:n2,0:n1,0:n3), stat=istat)
call memocc(istat, xye_c, 'xye_c', subname)
xya_c=0.d0
xyb_c=0.d0
xyc_c=0.d0
xye_c=0.d0

allocate(xza_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, xza_c, 'xza_c', subname)
allocate(xzb_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, xzb_c, 'xzb_c', subname)
allocate(xzc_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, xzc_c, 'xzc_c', subname)
allocate(xze_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, xze_c, 'xze_c', subname)
xza_c=0.d0
xzb_c=0.d0
xzc_c=0.d0
xze_c=0.d0

allocate(yza_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, yza_c, 'yza_c', subname)
allocate(yzb_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, yzb_c, 'yzb_c', subname)
allocate(yzc_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, yzc_c, 'yzc_c', subname)
allocate(yze_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, yze_c, 'yze_c', subname)
yza_c=0.d0
yzb_c=0.d0
yzc_c=0.d0
yze_c=0.d0

allocate(xya_f(3,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), stat=istat)
call memocc(istat, xya_f, 'xya_f', subname)
allocate(xyb_f(4,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), stat=istat)
call memocc(istat, xyb_f, 'xyb_f', subname)
allocate(xyc_f(3,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), stat=istat)
call memocc(istat, xyc_f, 'xyc_f', subname)
allocate(xye_f(4,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), stat=istat)
call memocc(istat, xye_f, 'xye_f', subname)
xya_f=0.d0
xyb_f=0.d0
xyc_f=0.d0
xye_f=0.d0

allocate(xza_f(3,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, xza_f, 'xza_f', subname)
allocate(xzb_f(4,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, xzb_f, 'xzb_f', subname)
allocate(xzc_f(3,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, xzc_f, 'xzc_f', subname)
allocate(xze_f(4,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, xze_f, 'xze_f', subname)
xza_f=0.d0
xzb_f=0.d0
xzc_f=0.d0
xze_f=0.d0

allocate(yza_f(3,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, yza_f, 'yza_f', subname)
allocate(yzb_f(4,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, yzb_f, 'yzb_f', subname)
allocate(yzc_f(3,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, yzc_f, 'yzc_f', subname)
allocate(yze_f(4,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, yze_f, 'yze_f', subname)
yza_f=0.d0
yzb_f=0.d0
yzc_f=0.d0
yze_f=0.d0

aeff0=0.d0 ; beff0=0.d0 ; ceff0=0.d0 ; eeff0=0.0
aeff1=0.d0 ; beff1=0.d0 ; ceff1=0.d0 ; eeff1=0.0
aeff2=0.d0 ; beff2=0.d0 ; ceff2=0.d0 ; eeff2=0.0
aeff3=0.d0 ; beff3=0.d0 ; ceff3=0.d0 ; eeff3=0.0

aeff0_2=0.d0 ; beff0_2=0.d0 ; ceff0_2=0.d0 ; eeff0_2=0.0
aeff1_2=0.d0 ; beff1_2=0.d0 ; ceff1_2=0.d0 ; eeff1_2=0.0
aeff2_2=0.d0 ; beff2_2=0.d0 ; ceff2_2=0.d0 ; eeff2_2=0.0
aeff3_2=0.d0 ; beff3_2=0.d0 ; ceff3_2=0.d0 ; eeff3_2=0.0


  !!$!$omp parallel default(private) &
  !!$!$omp shared(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3) &
  !!$!$omp shared(cprecr,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,x_c,x_f,y_c,y_f)& 
  !!$!$omp shared(x_f1,x_f2,x_f3,a,b,c,e)
    !!!$omp do  

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
    do i3=0,n3
       do i2=0,n2
          if (ibyz_c(2,i2,i3)-ibyz_c(1,i2,i3).ge.4) then
             do i1=ibyz_c(1,i2,i3),ibyz_c(2,i2,i3)-4,4
                dyi0=0.0_wp 
                dyi1=0.0_wp 
                dyi2=0.0_wp 
                dyi3=0.0_wp 

                tt0a0=0.d0
                tt0a1=0.d0
                tt0a2=0.d0
                tt0a3=0.d0

                tt0b0=0.d0
                tt0b1=0.d0
                tt0b2=0.d0
                tt0b3=0.d0

                tt0c0=0.d0
                tt0c1=0.d0
                tt0c2=0.d0
                tt0c3=0.d0

                tt0e0=0.d0
                tt0e1=0.d0
                tt0e2=0.d0
                tt0e3=0.d0
  
                do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_c(2,i2,i3))
                   dyi0=dyi0 + xx_c(t,i2,i3)*aeff0array(t-i1-0,i1+0)
                   dyi1=dyi1 + xx_c(t,i2,i3)*aeff0array(t-i1-1,i1+1)
                   dyi2=dyi2 + xx_c(t,i2,i3)*aeff0array(t-i1-2,i1+2)
                   dyi3=dyi3 + xx_c(t,i2,i3)*aeff0array(t-i1-3,i1+3)

                   ! sss coefficients
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
                enddo
                y_c(i1+0,i2,i3)=dyi0+cprecr*xx_c(i1+0,i2,i3)
                y_c(i1+1,i2,i3)=dyi1+cprecr*xx_c(i1+1,i2,i3)
                y_c(i1+2,i2,i3)=dyi2+cprecr*xx_c(i1+2,i2,i3)
                y_c(i1+3,i2,i3)=dyi3+cprecr*xx_c(i1+3,i2,i3)

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
             enddo
             icur=i1
          else
             icur=ibyz_c(1,i2,i3)
          endif
  
          do i1=icur,ibyz_c(2,i2,i3)
             dyi=0.0_wp 
             tt0a0=0.d0
             tt0b0=0.d0
             tt0c0=0.d0
             tt0e0=0.d0
             ! Get the effective a-filters for the x dimension
             do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_c(2,i2,i3))
                dyi=dyi + xx_c(t,i2,i3)*aeff0array(t-i1,i1)
                ! sss coefficients
                tt0a0=tt0a0 + xx_c(t,i2,i3)*aeff0_2auxarray(t-i1,i1)
                tt0b0=tt0b0 + xx_c(t,i2,i3)*beff0_2auxarray(t-i1,i1)
                tt0c0=tt0c0 + xx_c(t,i2,i3)*ceff0_2auxarray(t-i1,i1)
                tt0e0=tt0e0 + xx_c(t,i2,i3)*eeff0_2auxarray(t-i1,i1)
             enddo
             y_c(i1,i2,i3)=dyi+cprecr*xx_c(i1,i2,i3)

             xya_c(i2,i1,i3)=tt0a0
             xza_c(i3,i1,i2)=tt0a0

             xyb_c(i2,i1,i3)=tt0b0
             xzb_c(i3,i1,i2)=tt0b0

             xyc_c(i2,i1,i3)=tt0c0
             xzc_c(i3,i1,i2)=tt0c0

             xyc_c(i2,i1,i3)=tt0c0
             xzc_c(i3,i1,i2)=tt0c0
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
             dyi=0.0_wp
             tt0=0.0_wp
             do t=max(ibyz_f(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_f(2,i2,i3))
                dyi=dyi + xx_f1(t,i2,i3)*beff0array(t-i1,i1)
             enddo
             y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
          enddo
  
           if (ibyz_c(2,i2,i3)-ibyz_c(1,i2,i3).ge.4) then
             do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)-4,4
                dyi0=0.0_wp 
                dyi1=0.0_wp 
                dyi2=0.0_wp 
                dyi3=0.0_wp 
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
             dyi=0.0_wp 
             tt0=0.0_wp 
             do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_c(2,i2,i3))
                dyi=dyi + xx_c(t,i2,i3)*ceff0array(t-i1,i1)
             enddo
             y_f(1,i1,i2,i3)=dyi
          enddo
       enddo
    enddo
    !!!$omp enddo
  
  
  
  
    ! wavelet part
  
    !!!$omp do
    do i3=nfl3,nfu3
       do i2=nfl2,nfu2
          do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
             t112=0.0_wp;t121=0.0_wp;t122=0.0_wp;t212=0.0_wp;t221=0.0_wp;t222=0.0_wp;t211=0.0_wp 
             tt112=0.0_wp;tt121=0.0_wp;tt122=0.0_wp;tt212=0.0_wp;tt221=0.0_wp;tt222=0.0_wp;tt211=0.0_wp 
             tt1a0=0.d0
             tt1b0=0.d0
             tt1c0=0.d0
             tt1e0=0.d0

             tt2a0=0.d0
             tt2b0=0.d0
             tt2c0=0.d0
             tt2e0=0.d0

             tt3a0=0.d0
             tt3b0=0.d0
             tt3c0=0.d0
             tt3e0=0.d0

             tt4a0=0.d0
             tt4b0=0.d0
             tt4c0=0.d0
             tt4e0=0.d0

             tt5a0=0.d0
             tt5b0=0.d0
             tt5c0=0.d0
             tt5e0=0.d0

             tt6a0=0.d0
             tt6b0=0.d0
             tt6c0=0.d0
             tt6e0=0.d0

             tt7a0=0.d0
             tt7b0=0.d0
             tt7c0=0.d0
             tt7e0=0.d0
             do l=max(nfl1-i1,lowfil),min(lupfil,nfu1-i1)
                t112=t112 + xx_f(4,i1+l,i2,i3)*aeff0array(l,i1) + xx_f(5,i1+l,i2,i3)*beff0array(l,i1)
                t121=t121 + xx_f(2,i1+l,i2,i3)*aeff0array(l,i1) + xx_f(3,i1+l,i2,i3)*beff0array(l,i1)
                t122=t122 + xx_f(6,i1+l,i2,i3)*aeff0array(l,i1) + xx_f(7,i1+l,i2,i3)*beff0array(l,i1)
                t212=t212 + xx_f(4,i1+l,i2,i3)*ceff0array(l,i1) + xx_f(5,i1+l,i2,i3)*eeff0array(l,i1)
                t221=t221 + xx_f(2,i1+l,i2,i3)*ceff0array(l,i1) + xx_f(3,i1+l,i2,i3)*eeff0array(l,i1)
                t222=t222 + xx_f(6,i1+l,i2,i3)*ceff0array(l,i1) + xx_f(7,i1+l,i2,i3)*eeff0array(l,i1)
                t211=t211 + xx_f(1,i1+l,i2,i3)*eeff0array(l,i1)
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
             y_f(4,i1,i2,i3)=t112+cprecr*xx_f(4,i1,i2,i3)
             y_f(2,i1,i2,i3)=t121+cprecr*xx_f(2,i1,i2,i3)
             y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211+cprecr*xx_f(1,i1,i2,i3)
             y_f(6,i1,i2,i3)=t122+cprecr*xx_f(6,i1,i2,i3)
             y_f(5,i1,i2,i3)=t212+cprecr*xx_f(5,i1,i2,i3)
             y_f(3,i1,i2,i3)=t221+cprecr*xx_f(3,i1,i2,i3)
             y_f(7,i1,i2,i3)=t222+cprecr*xx_f(7,i1,i2,i3)
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
          enddo
       enddo
    enddo
    !!!$omp enddo


  
  

  
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
  
  
    ! + (1/2) d^2/dy^2
    !!!$omp do
    do i3=0,n3
       do i1=0,n1
          if (ibxz_c(2,i1,i3)-ibxz_c(1,i1,i3).ge.4) then
             do i2=ibxz_c(1,i1,i3),ibxz_c(2,i1,i3)-4,4
                dyi0=0.0_wp 
                dyi1=0.0_wp 
                dyi2=0.0_wp 
                dyi3=0.0_wp 

                tt0a0=0.d0
                tt0a1=0.d0
                tt0a2=0.d0
                tt0a3=0.d0

                tt0b0=0.d0
                tt0b1=0.d0
                tt0b2=0.d0
                tt0b3=0.d0

                tt0c0=0.d0
                tt0c1=0.d0
                tt0c2=0.d0
                tt0c3=0.d0

                tt0e0=0.d0
                tt0e1=0.d0
                tt0e2=0.d0
                tt0e3=0.d0
                do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_c(2,i1,i3))
                   dyi0=dyi0 + xy_c(t,i1,i3)*aeff0array(t-i2-0,i2+0) + 2.d0*xya_c(t,i1,i3)*aeff0_2array(t-i2-0,i2+0)
                   dyi1=dyi1 + xy_c(t,i1,i3)*aeff0array(t-i2-1,i2+1) + 2.d0*xya_c(t,i1,i3)*aeff0_2array(t-i2-1,i2+1)
                   dyi2=dyi2 + xy_c(t,i1,i3)*aeff0array(t-i2-2,i2+2) + 2.d0*xya_c(t,i1,i3)*aeff0_2array(t-i2-2,i2+2)
                   dyi3=dyi3 + xy_c(t,i1,i3)*aeff0array(t-i2-3,i2+3) + 2.d0*xya_c(t,i1,i3)*aeff0_2array(t-i2-3,i2+3)

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
                y_c(i1,i2+0,i3)=y_c(i1,i2+0,i3)+dyi0
                y_c(i1,i2+1,i3)=y_c(i1,i2+1,i3)+dyi1
                y_c(i1,i2+2,i3)=y_c(i1,i2+2,i3)+dyi2
                y_c(i1,i2+3,i3)=y_c(i1,i2+3,i3)+dyi3

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
             enddo
             icur=i2
          else
             icur=ibxz_c(1,i1,i3)
          endif
  
          do i2=icur,ibxz_c(2,i1,i3)
             dyi=0.0_wp 
             tt0a0=0.d0
             tt0b0=0.d0
             tt0c0=0.d0
             tt0e0=0.d0
             do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
                dyi=dyi + xy_c(t,i1,i3)*aeff0array(t-i2,i2) + 2.d0*xya_c(t,i1,i3)*aeff0_2array(t-i2,i2)

                tt0a0=tt0a0 + xy_c(t,i1,i3)*aeff0_2auxarray(t-i2-0,i2)

                tt0b0=tt0b0 + xy_c(t,i1,i3)*beff0_2auxarray(t-i2-0,i2)

                tt0c0=tt0c0 + xy_c(t,i1,i3)*ceff0_2auxarray(t-i2-0,i2)

                tt0e0=tt0e0 + xy_c(t,i1,i3)*eeff0_2auxarray(t-i2-0,i2)
             enddo
             y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi

             yza_c(i3,i1,i2)=tt0a0

             yzb_c(i3,i1,i2)=tt0b0

             yzc_c(i3,i1,i2)=tt0c0

             yze_c(i3,i1,i2)=tt0e0
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
                   dyi0=dyi0 + xy_f2(t,i1,i3)*beff0array(t-i2-0,i2+0) + 2.d0*xyb_f(1,t,i1,i3)*aeff0_2array(t-i2-0,i2+0) + &
                               2.d0*(xya_f(1,t,i1,i3)+xyb_f(2,t,i1,i3))*beff0_2array(t-i2-0,i2+0)
                   dyi1=dyi1 + xy_f2(t,i1,i3)*beff0array(t-i2-1,i2+1) + 2.d0*xyb_f(1,t,i1,i3)*aeff0_2array(t-i2-1,i2+1) + &
                               2.d0*(xya_f(1,t,i1,i3)+xyb_f(2,t,i1,i3))*beff0_2array(t-i2-1,i2+1)
                   dyi2=dyi2 + xy_f2(t,i1,i3)*beff0array(t-i2-2,i2+2) + 2.d0*xyb_f(1,t,i1,i3)*aeff0_2array(t-i2-2,i2+2) + &
                               2.d0*(xya_f(1,t,i1,i3)+xyb_f(2,t,i1,i3))*beff0_2array(t-i2-2,i2+2)
                   dyi3=dyi3 + xy_f2(t,i1,i3)*beff0array(t-i2-3,i2+3) + 2.d0*xyb_f(1,t,i1,i3)*aeff0_2array(t-i2-3,i2+3) + &
                               2.d0*(xya_f(1,t,i1,i3)+xyb_f(2,t,i1,i3))*beff0_2array(t-i2-3,i2+3)
                enddo
                y_c(i1,i2+0,i3)=y_c(i1,i2+0,i3)+dyi0
                y_c(i1,i2+1,i3)=y_c(i1,i2+1,i3)+dyi1
                y_c(i1,i2+2,i3)=y_c(i1,i2+2,i3)+dyi2
                y_c(i1,i2+3,i3)=y_c(i1,i2+3,i3)+dyi3
             enddo
             istart=i2
          endif
  
          do i2=istart,iend
             dyi0=0.0_wp
             do t=max(ibxz_f(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_f(2,i1,i3))
                dyi0=dyi0 + xy_f2(t,i1,i3)*beff0array(t-i2-0,i2) + 2.d0*xyb_f(1,t,i1,i3)*aeff0_2array(t-i2-0,i2) + &
                            2.d0*(xya_f(1,t,i1,i3)+xyb_f(2,t,i1,i3))*beff0_2array(t-i2-0,i2)
             enddo
             y_c(i1,i2+0,i3)=y_c(i1,i2+0,i3)+dyi0
          enddo
  
           if (ibxz_f(2,i1,i3)-ibxz_f(1,i1,i3).ge.4) then
             do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)-4,4
                dyi0=0.0_wp 
                dyi1=0.0_wp 
                dyi2=0.0_wp 
                dyi3=0.0_wp 
                tt10=0.0_wp 
                tt11=0.0_wp 
                tt12=0.0_wp 
                tt13=0.0_wp 
                tt20=0.0_wp 
                tt21=0.0_wp 
                tt22=0.0_wp 
                tt23=0.0_wp 
                tt30=0.0_wp 
                tt31=0.0_wp 
                tt32=0.0_wp 
                tt33=0.0_wp 
                do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_c(2,i1,i3))
                   dyi0=dyi0 + xy_c(t,i1,i3)*ceff0array(t-i2-0,i2+0)
                   dyi1=dyi1 + xy_c(t,i1,i3)*ceff0array(t-i2-1,i2+1)
                   dyi2=dyi2 + xy_c(t,i1,i3)*ceff0array(t-i2-2,i2+2)
                   dyi3=dyi3 + xy_c(t,i1,i3)*ceff0array(t-i2-3,i2+3)
  
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
                y_f(2,i1,i2+0,i3)=y_f(2,i1,i2+0,i3)+dyi0
                y_f(2,i1,i2+1,i3)=y_f(2,i1,i2+1,i3)+dyi1
                y_f(2,i1,i2+2,i3)=y_f(2,i1,i2+2,i3)+dyi2
                y_f(2,i1,i2+3,i3)=y_f(2,i1,i2+3,i3)+dyi3
  
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
             enddo
             icur=i2
          else
             icur=ibxz_f(1,i1,i3)
          endif
  
          do i2=icur,ibxz_f(2,i1,i3)
             dyi0=0.0_wp 
             tt10=0.0_wp 
             tt20=0.0_wp 
             tt30=0.0_wp 
             do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
                dyi0=dyi0 + xy_c(t,i1,i3)*ceff0array(t-i2-0,i2)

                tt10=tt10 + 2.d0*xyc_c(t,i1,i3)*aeff0_2array(t-i2-0,i2)

                tt20=tt20 + 2.d0*xya_c(t,i1,i3)*ceff0_2array(t-i2-0,i2)

                tt30=tt30 + 2.d0*xyc_c(t,i1,i3)*ceff0_2array(t-i2-0,i2)
             enddo
             y_f(2,i1,i2+0,i3)=y_f(2,i1,i2+0,i3)+dyi0
             y_f(1,i1,i2+0,i3)=y_f(1,i1,i2+0,i3)+tt10
             y_f(2,i1,i2+0,i3)=y_f(2,i1,i2+0,i3)+tt20
             y_f(3,i1,i2+0,i3)=y_f(3,i1,i2+0,i3)+tt30
          enddo
       enddo
    enddo
    !!!$omp enddo
  
  
    ! wavelet part
  
    !!!$omp do
    do i3=nfl3,nfu3
       do i1=nfl1,nfu1
          do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
             ! Get the effective filters for the y dimension
             tt10 = 0.d0
             tt20 = 0.d0
             tt30 = 0.d0
             tt40 = 0.d0
             tt50 = 0.d0
             tt60 = 0.d0
             tt70 = 0.d0

             tt1a0=0.d0
             tt1b0=0.d0
             tt1c0=0.d0
             tt1e0=0.d0

             tt2a0=0.d0
             tt2b0=0.d0
             tt2c0=0.d0
             tt2e0=0.d0

             tt3a0=0.d0
             tt3b0=0.d0
             tt3c0=0.d0
             tt3e0=0.d0

             tt4a0=0.d0
             tt4b0=0.d0
             tt4c0=0.d0
             tt4e0=0.d0

             tt5a0=0.d0
             tt5b0=0.d0
             tt5c0=0.d0
             tt5e0=0.d0

             tt6a0=0.d0
             tt6b0=0.d0
             tt6c0=0.d0
             tt6e0=0.d0

             tt7a0=0.d0
             tt7b0=0.d0
             tt7c0=0.d0
             tt7e0=0.d0
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
             y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+tt40
             y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+tt20
             y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+tt10
             y_f(6,i1,i2,i3)=y_f(6,i1,i2,i3)+tt60
             y_f(5,i1,i2,i3)=y_f(5,i1,i2,i3)+tt50
             y_f(3,i1,i2,i3)=y_f(3,i1,i2,i3)+tt30
             y_f(7,i1,i2,i3)=y_f(7,i1,i2,i3)+tt70

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
          enddo
       enddo
    enddo
    !!!$omp enddo




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


  ! + (1/2) d^2/dz^2

  !!!$omp do
  do i2=0,n2
     do i1=0,n1
        if (ibxy_c(2,i1,i2)-ibxy_c(1,i1,i2).ge.4) then
           do i3=ibxy_c(1,i1,i2),ibxy_c(2,i1,i2)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_c(2,i1,i2))
                 dyi0=dyi0 + xz_c(t,i1,i2)*aeff0array(t-i3-0,i3+0) + 2.d0*(xza_c(t,i1,i2)+yza_c(t,i1,i2))*aeff0_2array(t-i3-0,i3+0)
                 dyi1=dyi1 + xz_c(t,i1,i2)*aeff0array(t-i3-1,i3+1) + 2.d0*(xza_c(t,i1,i2)+yza_c(t,i1,i2))*aeff0_2array(t-i3-1,i3+1)
                 dyi2=dyi2 + xz_c(t,i1,i2)*aeff0array(t-i3-2,i3+2) + 2.d0*(xza_c(t,i1,i2)+yza_c(t,i1,i2))*aeff0_2array(t-i3-2,i3+2)
                 dyi3=dyi3 + xz_c(t,i1,i2)*aeff0array(t-i3-3,i3+3) + 2.d0*(xza_c(t,i1,i2)+yza_c(t,i1,i2))*aeff0_2array(t-i3-3,i3+3)
              enddo
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
           do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_c(2,i1,i2))
              dyi0=dyi0 + xz_c(t,i1,i2)*aeff0array(t-i3-0,i3) + 2.d0*(xza_c(t,i1,i2)+yza_c(t,i1,i2))*aeff0_2array(t-i3-0,i3+0)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi0
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
                 dyi0 = dyi0 + xz_f4(t,i1,i2)*beff0array(t-i3-0,i3+0) + &
                               2.d0*(xzb_f(1,t,i1,i2)+yzb_f(1,t,i1,i2))*aeff0_2array(t-i3-0,i3+0) + &
                               2.d0*(xza_f(2,t,i1,i2)+xzb_f(3,t,i1,i2)+yza_f(2,t,i1,i2)+yzb_f(3,t,i1,i2))*beff0_2array(t-i3-0,i3+0)
                 dyi1 = dyi1 + xz_f4(t,i1,i2)*beff0array(t-i3-1,i3+1) + &
                               2.d0*(xzb_f(1,t,i1,i2)+yzb_f(1,t,i1,i2))*aeff0_2array(t-i3-1,i3+1) + &
                               2.d0*(xza_f(2,t,i1,i2)+xzb_f(3,t,i1,i2)+yza_f(2,t,i1,i2)+yzb_f(3,t,i1,i2))*beff0_2array(t-i3-1,i3+1)
                 dyi2 = dyi2 + xz_f4(t,i1,i2)*beff0array(t-i3-2,i3+2) + &
                               2.d0*(xzb_f(1,t,i1,i2)+yzb_f(1,t,i1,i2))*aeff0_2array(t-i3-2,i3+2) + &
                               2.d0*(xza_f(2,t,i1,i2)+xzb_f(3,t,i1,i2)+yza_f(2,t,i1,i2)+yzb_f(3,t,i1,i2))*beff0_2array(t-i3-2,i3+2)
                 dyi3 = dyi3 + xz_f4(t,i1,i2)*beff0array(t-i3-3,i3+3) + &
                               2.d0*(xzb_f(1,t,i1,i2)+yzb_f(1,t,i1,i2))*aeff0_2array(t-i3-3,i3+3) + &
                               2.d0*(xza_f(2,t,i1,i2)+xzb_f(3,t,i1,i2)+yza_f(2,t,i1,i2)+yzb_f(3,t,i1,i2))*beff0_2array(t-i3-3,i3+3)
              enddo
              y_c(i1,i2,i3+0)=y_c(i1,i2,i3+0)+dyi0
              y_c(i1,i2,i3+1)=y_c(i1,i2,i3+1)+dyi1
              y_c(i1,i2,i3+2)=y_c(i1,i2,i3+2)+dyi2
              y_c(i1,i2,i3+3)=y_c(i1,i2,i3+3)+dyi3
           enddo
           istart=i2
        endif

        do i3=istart,iend
           dyi0=0.0_wp
           tt0=0.0_wp
           do t=max(ibxy_f(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_f(2,i1,i2))
              dyi0=dyi0 + xz_f4(t,i1,i2)*beff0array(t-i3-0,i3) + &
                          2.d0*(xzb_f(1,t,i1,i2)+yzb_f(1,t,i1,i2))*aeff0_2array(t-i3-0,i3) + &
                          2.d0*(xza_f(2,t,i1,i2)+xzb_f(3,t,i1,i2)+yza_f(2,t,i1,i2)+yzb_f(3,t,i1,i2))*beff0_2array(t-i3-0,i3)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi0
        enddo

         if (ibxy_f(2,i1,i2)-ibxy_f(1,i1,i2).ge.4) then
           do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              tt10 = 0.d0
              tt11 = 0.d0
              tt12 = 0.d0
              tt13 = 0.d0
              tt40 = 0.d0
              tt41 = 0.d0
              tt42 = 0.d0
              tt43 = 0.d0
              tt50 = 0.d0
              tt51 = 0.d0
              tt52 = 0.d0
              tt53 = 0.d0
              tt20 = 0.d0
              tt21 = 0.d0
              tt22 = 0.d0
              tt23 = 0.d0
              tt60 = 0.d0
              tt61 = 0.d0
              tt62 = 0.d0
              tt63 = 0.d0
              do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_c(2,i1,i2))
                 dyi0=dyi0 + xz_c(t,i1,i2)*ceff0array(t-i3-0,i3+0)
                 dyi1=dyi1 + xz_c(t,i1,i2)*ceff0array(t-i3-1,i3+1)
                 dyi2=dyi2 + xz_c(t,i1,i2)*ceff0array(t-i3-2,i3+2)
                 dyi3=dyi3 + xz_c(t,i1,i2)*ceff0array(t-i3-3,i3+3)

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
              y_f(4,i1,i2,i3+0) = y_f(4,i1,i2,i3+0) + dyi0 + tt40
              y_f(4,i1,i2,i3+1) = y_f(4,i1,i2,i3+1) + dyi1 + tt41
              y_f(4,i1,i2,i3+2) = y_f(4,i1,i2,i3+2) + dyi2 + tt42
              y_f(4,i1,i2,i3+3) = y_f(4,i1,i2,i3+3) + dyi3 + tt43

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
           enddo
           icur=i3
        else
           icur=ibxy_f(1,i1,i2)
        endif

        do i3=icur,ibxy_f(2,i1,i2)
           dyi0=0.0_wp 
           tt10 = 0.d0
           tt40 = 0.d0
           tt50 = 0.d0
           tt20 = 0.d0
           tt60 = 0.d0
           do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_c(2,i1,i2))
              dyi0=dyi0 + xz_c(t,i1,i2)*ceff0array(t-i3-0,i3)

              tt10 = tt10 + 2.d0*xzc_c(t,i1,i2)*aeff0_2array(t-i3-0,i3)

              tt40 = tt40 + 2.d0*xza_c(t,i1,i2)*ceff0_2array(t-i3-0,i3)

              tt50 = tt50 + 2.d0*xzc_c(t,i1,i2)*ceff0_2array(t-i3-0,i3)

              tt20 = tt20 + 2.d0*yzc_c(t,i1,i2)*aeff0_2array(t-i3-0,i3)

              tt40 = tt40 + 2.d0*yza_c(t,i1,i2)*ceff0_2array(t-i3-0,i3)

              tt60 = tt60 + 2.d0*yzc_c(t,i1,i2)*ceff0_2array(t-i3-0,i3)
           enddo
           y_f(4,i1,i2,i3+0) = y_f(4,i1,i2,i3+0) + dyi0 + tt40
           y_f(1,i1,i2,i3+0) = y_f(1,i1,i2,i3+0) + tt10
           y_f(5,i1,i2,i3+0) = y_f(5,i1,i2,i3+0) + tt50
           y_f(2,i1,i2,i3+0) = y_f(2,i1,i2,i3+0) + tt20
           y_f(6,i1,i2,i3+0) = y_f(6,i1,i2,i3+0) + tt60
        enddo
     enddo
  enddo
  !!!$omp enddo


  ! wavelet part

  !!!$omp do
  do i2=nfl2,nfu2
     do i1=nfl1,nfu1
        do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
           tt10 = 0.d0
           tt20 = 0.d0
           tt30 = 0.d0
           tt40 = 0.d0
           tt50 = 0.d0
           tt60 = 0.d0
           tt70 = 0.d0
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
  !!!$omp enddo


  iall=-product(shape(aeff0array))*kind(aeff0array)
  deallocate(aeff0array, stat=istat)
  call memocc(istat, iall, 'aeff0array', subname)

  iall=-product(shape(beff0array))*kind(beff0array)
  deallocate(beff0array, stat=istat)
  call memocc(istat, iall, 'beff0array', subname)

  iall=-product(shape(ceff0array))*kind(ceff0array)
  deallocate(ceff0array, stat=istat)
  call memocc(istat, iall, 'ceff0array', subname)

  iall=-product(shape(eeff0array))*kind(eeff0array)
  deallocate(eeff0array, stat=istat)
  call memocc(istat, iall, 'eeff0array', subname)


  iall=-product(shape(aeff0_2array))*kind(aeff0_2array)
  deallocate(aeff0_2array, stat=istat)
  call memocc(istat, iall, 'aeff0_2array', subname)

  iall=-product(shape(beff0_2array))*kind(beff0_2array)
  deallocate(beff0_2array, stat=istat)
  call memocc(istat, iall, 'beff0_2array', subname)

  iall=-product(shape(ceff0_2array))*kind(ceff0_2array)
  deallocate(ceff0_2array, stat=istat)
  call memocc(istat, iall, 'ceff0_2array', subname)

  iall=-product(shape(eeff0_2array))*kind(eeff0_2array)
  deallocate(eeff0_2array, stat=istat)
  call memocc(istat, iall, 'eeff0_2array', subname)


  iall=-product(shape(aeff0_2auxarray))*kind(aeff0_2auxarray)
  deallocate(aeff0_2auxarray, stat=istat)
  call memocc(istat, iall, 'aeff0_2auxarray', subname)

  iall=-product(shape(beff0_2auxarray))*kind(beff0_2auxarray)
  deallocate(beff0_2auxarray, stat=istat)
  call memocc(istat, iall, 'beff0_2auxarray', subname)

  iall=-product(shape(ceff0_2auxarray))*kind(ceff0_2auxarray)
  deallocate(ceff0_2auxarray, stat=istat)
  call memocc(istat, iall, 'ceff0_2auxarray', subname)

  iall=-product(shape(eeff0_2auxarray))*kind(eeff0_2auxarray)
  deallocate(eeff0_2auxarray, stat=istat)
  call memocc(istat, iall, 'eeff0_2auxarray', subname)


  iall=-product(shape(xya_c))*kind(xya_c)
  deallocate(xya_c, stat=istat)
  call memocc(istat, iall, 'xya_c', subname)

  iall=-product(shape(xyb_c))*kind(xyb_c)
  deallocate(xyb_c, stat=istat)
  call memocc(istat, iall, 'xyb_c', subname)

  iall=-product(shape(xyc_c))*kind(xyc_c)
  deallocate(xyc_c, stat=istat)
  call memocc(istat, iall, 'xyc_c', subname)

  iall=-product(shape(xye_c))*kind(xye_c)
  deallocate(xye_c, stat=istat)
  call memocc(istat, iall, 'xye_c', subname)



  iall=-product(shape(xza_c))*kind(xza_c)
  deallocate(xza_c, stat=istat)
  call memocc(istat, iall, 'xza_c', subname)

  iall=-product(shape(xzb_c))*kind(xzb_c)
  deallocate(xzb_c, stat=istat)
  call memocc(istat, iall, 'xzb_c', subname)

  iall=-product(shape(xzc_c))*kind(xzc_c)
  deallocate(xzc_c, stat=istat)
  call memocc(istat, iall, 'xzc_c', subname)

  iall=-product(shape(xze_c))*kind(xze_c)
  deallocate(xze_c, stat=istat)
  call memocc(istat, iall, 'xze_c', subname)



  iall=-product(shape(yza_c))*kind(yza_c)
  deallocate(yza_c, stat=istat)
  call memocc(istat, iall, 'yza_c', subname)

  iall=-product(shape(yzb_c))*kind(yzb_c)
  deallocate(yzb_c, stat=istat)
  call memocc(istat, iall, 'yzb_c', subname)

  iall=-product(shape(yzc_c))*kind(yzc_c)
  deallocate(yzc_c, stat=istat)
  call memocc(istat, iall, 'yzc_c', subname)

  iall=-product(shape(yze_c))*kind(yze_c)
  deallocate(yze_c, stat=istat)
  call memocc(istat, iall, 'yze_c', subname)


  iall=-product(shape(xya_f))*kind(xya_f)
  deallocate(xya_f, stat=istat)
  call memocc(istat, iall, 'xya_f', subname)
  iall=-product(shape(xyb_f))*kind(xyb_f)
  deallocate(xyb_f, stat=istat)
  call memocc(istat, iall, 'xyb_f', subname)
  iall=-product(shape(xyc_f))*kind(xyc_f)
  deallocate(xyc_f, stat=istat)
  call memocc(istat, iall, 'xyc_f', subname)
  iall=-product(shape(xye_f))*kind(xye_f)
  deallocate(xye_f, stat=istat)
  call memocc(istat, iall, 'yze_f7', subname)

  iall=-product(shape(xza_f))*kind(xza_f)
  deallocate(xza_f, stat=istat)
  call memocc(istat, iall, 'xza_f', subname)
  iall=-product(shape(xzb_f))*kind(xzb_f)
  deallocate(xzb_f, stat=istat)
  call memocc(istat, iall, 'xzb_f', subname)
  iall=-product(shape(xzc_f))*kind(xzc_f)
  deallocate(xzc_f, stat=istat)
  call memocc(istat, iall, 'xzc_f', subname)
  iall=-product(shape(xze_f))*kind(xze_f)
  deallocate(xze_f, stat=istat)
  call memocc(istat, iall, 'zze_f7', subname)

  iall=-product(shape(yza_f))*kind(yza_f)
  deallocate(yza_f, stat=istat)
  call memocc(istat, iall, 'yza_f', subname)
  iall=-product(shape(yzb_f))*kind(yzb_f)
  deallocate(yzb_f, stat=istat)
  call memocc(istat, iall, 'yzb_f', subname)
  iall=-product(shape(yzc_f))*kind(yzc_f)
  deallocate(yzc_f, stat=istat)
  call memocc(istat, iall, 'yzc_f', subname)
  iall=-product(shape(yze_f))*kind(yze_f)
  deallocate(yze_f, stat=istat)
  call memocc(istat, iall, 'zze_f7', subname)


  call timing(iproc,'convolQuartic ','OF')


END SUBROUTINE ConvolQuartic4






!>  Applies the following operation: 
!!  y = [((r-r0)^6)]*x
!!! WARNING: MIGHT CONTAIN A BUG!
subroutine ConvolSextic(n1, n2, n3, &
     nfl1, nfu1, nfl2, nfu2, nfl3, nfu3,  &
     hgrid, offsetx, offsety, offsetz, &
     ibyz_c, ibxz_c, ibxy_c, ibyz_f, ibxz_f, ibxy_f, &
     rxyzConf, potentialPrefac, withKinetic, cprecr, &
     xx_c, xx_f1, xx_f, &
     xy_c, xy_f2, xy_f, &
     xz_c, xz_f4, xz_f, &
     y_c, y_f)

  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, offsetx, offsety, offsetz
  real(gp), intent(in) :: hgrid, potentialPrefac, cprecr
  logical,intent(in):: withKinetic
  real(8),dimension(3):: rxyzConf
  integer, dimension(2,0:n2,0:n3), intent(in) :: ibyz_c,ibyz_f
  integer, dimension(2,0:n1,0:n3), intent(in) :: ibxz_c,ibxz_f
  integer, dimension(2,0:n1,0:n2), intent(in) :: ibxy_c,ibxy_f
  real(wp),dimension(0:n1,0:n2,0:n3),intent(in):: xx_c
  real(wp),dimension(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3),intent(in):: xx_f1
  real(wp),dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3),intent(in):: xx_f
  real(wp),dimension(0:n2,0:n1,0:n3),intent(in):: xy_c
  real(wp),dimension(nfl2:nfu2,nfl1:nfu1,nfl3:nfu3),intent(in):: xy_f2
  real(wp),dimension(7,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3),intent(in):: xy_f
  real(wp),dimension(0:n3,0:n1,0:n2),intent(in):: xz_c
  real(wp),dimension(nfl3:nfu3,nfl1:nfu1,nfl2:nfu2),intent(in):: xz_f4
  real(wp),dimension(7,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2),intent(in):: xz_f
  real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: y_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(out) :: y_f
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  !logical :: firstcall=.true. 
  !integer, save :: mflop1,mflop2,mflop3,nflop1,nflop2,nflop3
  !integer :: ncount1,ncount_rate,ncount_max,ncount2,ncount3,ncount4,ncount5,ncount6
  integer :: i,t,i1,i2,i3
  integer :: icur,istart,iend,l
  real(wp) :: scale,dyi,dyi0,dyi1,dyi2,dyi3,t112,t121,t122,t212,t221,t222,t211
  real(wp):: tt112, tt121, tt122, tt212, tt221, tt222, tt211
  real(wp):: tt0, tt1, tt2, tt3, tt4, tt5, tt6, tt7
  real(wp), dimension(-3+lowfil:lupfil+3) :: a, aeff0, aeff1, aeff2, aeff3
  real(wp), dimension(-3+lowfil:lupfil+3) :: b, beff0, beff1, beff2, beff3
  real(wp), dimension(-3+lowfil:lupfil+3) :: c, ceff0, ceff1, ceff2, ceff3
  real(wp), dimension(lowfil:lupfil) :: e, eeff0, eeff1, eeff2, eeff3
  real(wp), dimension(-3+lowfil:lupfil+3) :: aeff0_2, aeff1_2, aeff2_2, aeff3_2
  real(wp), dimension(-3+lowfil:lupfil+3) :: beff0_2, beff1_2, beff2_2, beff3_2
  real(wp), dimension(-3+lowfil:lupfil+3) :: ceff0_2, ceff1_2, ceff2_2, ceff3_2
  real(wp),dimension(:,:),allocatable:: aeff0array, beff0array, ceff0array, eeff0array
  real(wp),dimension(:,:),allocatable:: aeff0_2array, beff0_2array, ceff0_2array, eeff0_2array
  real(wp),dimension(:,:),allocatable:: aeff0_4array, beff0_4array, ceff0_4array, eeff0_4array
  real(wp),dimension(:,:),allocatable:: aeff0_2auxarray, beff0_2auxarray, ceff0_2auxarray, eeff0_2auxarray
  real(wp),dimension(:,:),allocatable:: aeff0_4auxarray, beff0_4auxarray, ceff0_4auxarray, eeff0_4auxarray
  real(wp), dimension(lowfil:lupfil) :: eeff0_2, eeff1_2, eeff2_2, eeff3_2
  real(wp),dimension(:,:,:),allocatable:: x2ya_c, x2yb_c, x2yc_c, x2ye_c, x2za_c, x2zb_c, x2zc_c, x2ze_c, &
                                          y2za_c, y2zb_c, y2zc_c, y2ze_c
  real(wp),dimension(:,:,:),allocatable:: x4ya_c, x4yb_c, x4yc_c, x4ye_c, x4za_c, x4zb_c, x4zc_c, x4ze_c, &
                                          y4za_c, y4zb_c, y4zc_c, y4ze_c
  real(wp),dimension(:,:,:,:),allocatable:: x2ya_f, x2yb_f, x2yc_f, x2ye_f, x2za_f, x2zb_f, x2zc_f, x2ze_f, &
                                            y2za_f, y2zb_f, y2zc_f, y2ze_f
  real(wp),dimension(:,:,:,:),allocatable:: x4ya_f, x4yb_f, x4yc_f, x4ye_f, x4za_f, x4zb_f, x4zc_f, x4ze_f, &
                                            y4za_f, y4zb_f, y4zc_f, y4ze_f
  real(wp),dimension(:,:,:,:),allocatable:: x2ya_f2, x2yb_f2, x2yc_f2, x2ye_f2
  real(wp),dimension(:,:,:),allocatable:: x2y2aa_c, x2y2ab_c, x2y2ac_c, x2y2ae_c
  real(wp),dimension(:,:,:),allocatable:: x2y2ba_c, x2y2bb_c, x2y2bc_c, x2y2be_c
  real(wp),dimension(:,:,:),allocatable:: x2y2ca_c, x2y2cb_c, x2y2cc_c, x2y2ce_c
  real(wp),dimension(:,:,:),allocatable:: x2y2ea_c, x2y2eb_c, x2y2ec_c, x2y2ee_c
  real(wp),dimension(:,:,:,:),allocatable:: x2y2aa_f, x2y2ab_f, x2y2ac_f, x2y2ae_f
  real(wp),dimension(:,:,:,:),allocatable:: x2y2ba_f, x2y2bb_f, x2y2bc_f, x2y2be_f
  real(wp),dimension(:,:,:,:),allocatable:: x2y2ca_f, x2y2cb_f, x2y2cc_f, x2y2ce_f
  real(wp),dimension(:,:,:,:),allocatable:: x2y2ea_f, x2y2eb_f, x2y2ec_f, x2y2ee_f

  !!real(wp),dimension(:,:,:),allocatable:: xy_c, xz_c
  !!real(wp),dimension(:,:,:),allocatable:: xy_f1, xy_f2, xy_f3, xy_f4, xy_f5, xy_f6, xy_f7
  !!real(wp),dimension(:,:,:),allocatable:: xz_f1, xz_f2, xz_f3, xz_f4, xz_f5, xz_f6, xz_f7
real(8):: x0, y0, z0
integer:: ii, istat, iall
character(len=*),parameter:: subname='ConvolSextic'
real(8):: tt00, tt01, tt02, tt03
real(8):: tt10, tt11, tt12, tt13
real(8):: tt20, tt21, tt22, tt23
real(8):: tt30, tt31, tt32, tt33
real(8):: tt40, tt41, tt42, tt43
real(8):: tt50, tt51, tt52, tt53
real(8):: tt60, tt61, tt62, tt63
real(8):: tt70, tt71, tt72, tt73
real(8):: tt0a0_2, tt0a1_2, tt0a2_2, tt0a3_2, tt0a0_4, tt0a1_4, tt0a2_4, tt0a3_4
real(8):: tt0b0_2, tt0b1_2, tt0b2_2, tt0b3_2, tt0b0_4, tt0b1_4, tt0b2_4, tt0b3_4
real(8):: tt0c0_2, tt0c1_2, tt0c2_2, tt0c3_2, tt0c0_4, tt0c1_4, tt0c2_4, tt0c3_4
real(8):: tt0e0_2, tt0e1_2, tt0e2_2, tt0e3_2, tt0e0_4, tt0e1_4, tt0e2_4, tt0e3_4
real(8):: tt1a0_2, tt1a1_2, tt1a2_2, tt1a3_2, tt1a0_4, tt1a1_4, tt1a2_4, tt1a3_4
real(8):: tt1b0_2, tt1b1_2, tt1b2_2, tt1b3_2, tt1b0_4, tt1b1_4, tt1b2_4, tt1b3_4
real(8):: tt1c0_2, tt1c1_2, tt1c2_2, tt1c3_2, tt1c0_4, tt1c1_4, tt1c2_4, tt1c3_4
real(8):: tt1e0_2, tt1e1_2, tt1e2_2, tt1e3_2, tt1e0_4, tt1e1_4, tt1e2_4, tt1e3_4
real(8):: tt2a0_2, tt2a1_2, tt2a2_2, tt2a3_2, tt2a0_4, tt2a1_4, tt2a2_4, tt2a3_4
real(8):: tt2b0_2, tt2b1_2, tt2b2_2, tt2b3_2, tt2b0_4, tt2b1_4, tt2b2_4, tt2b3_4
real(8):: tt2c0_2, tt2c1_2, tt2c2_2, tt2c3_2, tt2c0_4, tt2c1_4, tt2c2_4, tt2c3_4
real(8):: tt2e0_2, tt2e1_2, tt2e2_2, tt2e3_2, tt2e0_4, tt2e1_4, tt2e2_4, tt2e3_4
real(8):: tt3a0_2, tt3a1_2, tt3a2_2, tt3a3_2, tt3a0_4, tt3a1_4, tt3a2_4, tt3a3_4
real(8):: tt3b0_2, tt3b1_2, tt3b2_2, tt3b3_2, tt3b0_4, tt3b1_4, tt3b2_4, tt3b3_4
real(8):: tt3c0_2, tt3c1_2, tt3c2_2, tt3c3_2, tt3c0_4, tt3c1_4, tt3c2_4, tt3c3_4
real(8):: tt3e0_2, tt3e1_2, tt3e2_2, tt3e3_2, tt3e0_4, tt3e1_4, tt3e2_4, tt3e3_4
real(8):: tt4a0_2, tt4a1_2, tt4a2_2, tt4a3_2, tt4a0_4, tt4a1_4, tt4a2_4, tt4a3_4
real(8):: tt4b0_2, tt4b1_2, tt4b2_2, tt4b3_2, tt4b0_4, tt4b1_4, tt4b2_4, tt4b3_4
real(8):: tt4c0_2, tt4c1_2, tt4c2_2, tt4c3_2, tt4c0_4, tt4c1_4, tt4c2_4, tt4c3_4
real(8):: tt4e0_2, tt4e1_2, tt4e2_2, tt4e3_2, tt4e0_4, tt4e1_4, tt4e2_4, tt4e3_4
real(8):: tt5a0_2, tt5a1_2, tt5a2_2, tt5a3_2, tt5a0_4, tt5a1_4, tt5a2_4, tt5a3_4
real(8):: tt5b0_2, tt5b1_2, tt5b2_2, tt5b3_2, tt5b0_4, tt5b1_4, tt5b2_4, tt5b3_4
real(8):: tt5c0_2, tt5c1_2, tt5c2_2, tt5c3_2, tt5c0_4, tt5c1_4, tt5c2_4, tt5c3_4
real(8):: tt5e0_2, tt5e1_2, tt5e2_2, tt5e3_2, tt5e0_4, tt5e1_4, tt5e2_4, tt5e3_4
real(8):: tt6a0_2, tt6a1_2, tt6a2_2, tt6a3_2, tt6a0_4, tt6a1_4, tt6a2_4, tt6a3_4
real(8):: tt6b0_2, tt6b1_2, tt6b2_2, tt6b3_2, tt6b0_4, tt6b1_4, tt6b2_4, tt6b3_4
real(8):: tt6c0_2, tt6c1_2, tt6c2_2, tt6c3_2, tt6c0_4, tt6c1_4, tt6c2_4, tt6c3_4
real(8):: tt6e0_2, tt6e1_2, tt6e2_2, tt6e3_2, tt6e0_4, tt6e1_4, tt6e2_4, tt6e3_4
real(8):: tt7a0_2, tt7a1_2, tt7a2_2, tt7a3_2, tt7a0_4, tt7a1_4, tt7a2_4, tt7a3_4
real(8):: tt7b0_2, tt7b1_2, tt7b2_2, tt7b3_2, tt7b0_4
real(8):: tt7c0_2, tt7c1_2, tt7c2_2, tt7c3_2, tt7c0_4
real(8):: tt7e0_2,                            tt7e0_4
real(8):: ttaa0, ttab0, ttac0, ttae0
real(8):: ttaa1, ttab1, ttac1, ttae1
real(8):: ttaa2, ttab2, ttac2, ttae2
real(8):: ttaa3, ttab3, ttac3, ttae3
real(8):: ttaa4, ttab4, ttac4, ttae4
real(8):: ttaa5, ttab5, ttac5, ttae5
real(8):: ttaa6, ttab6, ttac6, ttae6
real(8):: ttaa7, ttab7, ttac7, ttae7
real(8):: ttba0, ttbb0, ttbc0, ttbe0
real(8):: ttba1, ttbb1, ttbc1, ttbe1
real(8):: ttba2, ttbb2, ttbc2, ttbe2
real(8):: ttba3, ttbb3, ttbc3, ttbe3
real(8):: ttba4, ttbb4, ttbc4, ttbe4
real(8):: ttba5, ttbb5, ttbc5, ttbe5
real(8):: ttba6, ttbb6, ttbc6, ttbe6
real(8):: ttba7, ttbb7, ttbc7, ttbe7
real(8):: ttca0, ttcb0, ttcc0, ttce0
real(8):: ttca1, ttcb1, ttcc1, ttce1
real(8):: ttca2, ttcb2, ttcc2, ttce2
real(8):: ttca3, ttcb3, ttcc3, ttce3
real(8):: ttca4, ttcb4, ttcc4, ttce4
real(8):: ttca5, ttcb5, ttcc5, ttce5
real(8):: ttca6, ttcb6, ttcc6, ttce6
real(8):: ttca7, ttcb7, ttcc7, ttce7
real(8):: ttea0, tteb0, ttec0, ttee0
real(8):: ttea1, tteb1, ttec1, ttee1
real(8):: ttea2, tteb2, ttec2, ttee2
real(8):: ttea3, tteb3, ttec3, ttee3
real(8):: ttea4, tteb4, ttec4, ttee4
real(8):: ttea5, tteb5, ttec5, ttee5
real(8):: ttea6, tteb6, ttec6, ttee6
real(8):: ttea7, tteb7, ttec7, ttee7
real(8):: tt1a0
real(8):: tt1b0
real(8):: tt1c0
real(8):: tt1e0
real(8):: tt2a0
real(8):: tt2b0
real(8):: tt2c0
real(8):: tt2e0
real(8):: tt3a0
real(8):: tt3b0
real(8):: tt3c0
real(8):: tt3e0
real(8):: tt4a0
real(8):: tt4b0
real(8):: tt4c0
real(8):: tt4e0
real(8):: tt5a0
real(8):: tt5b0
real(8):: tt5c0
real(8):: tt5e0
real(8):: tt6a0
real(8):: tt6b0
real(8):: tt6c0
real(8):: tt6e0
real(8):: tt7a0
real(8):: tt7b0
real(8):: tt7c0
real(8):: tt7e0

integer:: it=1!debug
real(8):: t1, t2, time, t3, t4, time2


  scale=-.5_wp/real(hgrid**2,wp)

  !---------------------------------------------------------------------------
  ! second derivative filters for Daubechies 16
  !  <phi|D^2|phi_i>
  a(0)=   -3.5536922899131901941296809374_wp*scale
  a(1)=    2.2191465938911163898794546405_wp*scale
  a(2)=   -0.6156141465570069496314853949_wp*scale
  a(3)=    0.2371780582153805636239247476_wp*scale
  a(4)=   -0.0822663999742123340987663521_wp*scale
  a(5)=    0.02207029188482255523789911295638968409_wp*scale
  a(6)=   -0.409765689342633823899327051188315485e-2_wp*scale
  a(7)=    0.45167920287502235349480037639758496e-3_wp*scale
  a(8)=   -0.2398228524507599670405555359023135e-4_wp*scale
  a(9)=    2.0904234952920365957922889447361e-6_wp*scale
  a(10)=  -3.7230763047369275848791496973044e-7_wp*scale
  a(11)=  -1.05857055496741470373494132287e-8_wp*scale
  a(12)=  -5.813879830282540547959250667e-11_wp*scale
  a(13)=   2.70800493626319438269856689037647576e-13_wp*scale
  a(14)=  -6.924474940639200152025730585882e-18_wp*scale

  a(15)=0.0_wp
  a(16)=0.0_wp 
  a(17)=0.0_wp
  
  do i=1,14+3
     a(-i)=a(i)
  enddo
  !  <phi|D^2|psi_i>
  c(-17)=0.0_wp
  c(-16)=0.0_wp
  c(-15)=0.0_wp
  
  c(-14)=     -3.869102413147656535541850057188e-18_wp*scale
  c(-13)=      1.5130616560866154733900029272077362e-13_wp*scale
  c(-12)=     -3.2264702314010525539061647271983988409e-11_wp*scale
  c(-11)=     -5.96264938781402337319841002642e-9_wp*scale
  c(-10)=     -2.1656830629214041470164889350342e-7_wp*scale
  c(-9 )=      8.7969704055286288323596890609625e-7_wp*scale
  c(-8 )=     -0.00001133456724516819987751818232711775_wp*scale
  c(-7 )=      0.00021710795484646138591610188464622454_wp*scale
  c(-6 )=     -0.0021356291838797986414312219042358542_wp*scale
  c(-5 )=      0.00713761218453631422925717625758502986_wp*scale
  c(-4 )=     -0.0284696165863973422636410524436931061_wp*scale
  c(-3 )=      0.14327329352510759457155821037742893841_wp*scale
  c(-2 )=     -0.42498050943780130143385739554118569733_wp*scale
  c(-1 )=      0.65703074007121357894896358254040272157_wp*scale
  c( 0 )=     -0.42081655293724308770919536332797729898_wp*scale
  c( 1 )=     -0.21716117505137104371463587747283267899_wp*scale
  c( 2 )=      0.63457035267892488185929915286969303251_wp*scale
  c( 3 )=     -0.53298223962800395684936080758073568406_wp*scale
  c( 4 )=      0.23370490631751294307619384973520033236_wp*scale
  c( 5 )=     -0.05657736973328755112051544344507997075_wp*scale
  c( 6 )=      0.0080872029411844780634067667008050127_wp*scale
  c( 7 )=     -0.00093423623304808664741804536808932984_wp*scale
  c( 8 )=      0.00005075807947289728306309081261461095_wp*scale
  c( 9 )=     -4.62561497463184262755416490048242e-6_wp*scale
  c( 10)=      6.3919128513793415587294752371778e-7_wp*scale
  c( 11)=      1.87909235155149902916133888931e-8_wp*scale
  c( 12)=      1.04757345962781829480207861447155543883e-10_wp*scale
  c( 13)=     -4.84665690596158959648731537084025836e-13_wp*scale
  c( 14)=      1.2392629629188986192855777620877e-17_wp*scale

  c(15)=0.0_wp
  c(16)=0.0_wp
  c(17)=0.0_wp
  !  <psi|D^2|phi_i>
  do i=-14-3,14+3
     b(i)=c(-i)
  enddo
  !<psi|D^2|psi_i>
  e(0)=   -24.875846029392331358907766562_wp*scale
  e(1)=   -7.1440597663471719869313377994_wp*scale
  e(2)=   -0.04251705323669172315864542163525830944_wp*scale
  e(3)=   -0.26995931336279126953587091167128839196_wp*scale
  e(4)=    0.08207454169225172612513390763444496516_wp*scale
  e(5)=   -0.02207327034586634477996701627614752761_wp*scale
  e(6)=    0.00409765642831595181639002667514310145_wp*scale
  e(7)=   -0.00045167920287507774929432548999880117_wp*scale
  e(8)=    0.00002398228524507599670405555359023135_wp*scale
  e(9)=   -2.0904234952920365957922889447361e-6_wp*scale
  e(10)=   3.7230763047369275848791496973044e-7_wp*scale
  e(11)=   1.05857055496741470373494132287e-8_wp*scale
  e(12)=   5.8138798302825405479592506674648873655e-11_wp*scale
  e(13)=  -2.70800493626319438269856689037647576e-13_wp*scale
  e(14)=   6.924474940639200152025730585882e-18_wp*scale
  do i=1,14
     e(-i)=e(i)
  enddo



! Allocate all arrays

i=max(n1,n2,n3)
allocate(aeff0array(-3+lowfil:lupfil+3,0:i), stat=istat)
call memocc(istat, aeff0array, 'aeff0array', subname)
allocate(beff0array(-3+lowfil:lupfil+3,0:i), stat=istat)
call memocc(istat, beff0array, 'beff0array', subname)
allocate(ceff0array(-3+lowfil:lupfil+3,0:i), stat=istat)
call memocc(istat, ceff0array, 'ceff0array', subname)
allocate(eeff0array(lowfil:lupfil,0:i), stat=istat)
call memocc(istat, eeff0array, 'eeff0array', subname)
aeff0array=0.d0
beff0array=0.d0
ceff0array=0.d0
eeff0array=0.d0

allocate(aeff0_2array(-3+lowfil:lupfil+3,0:i), stat=istat)
call memocc(istat, aeff0_2array, 'aeff0_2array', subname)
allocate(beff0_2array(-3+lowfil:lupfil+3,0:i), stat=istat)
call memocc(istat, beff0_2array, 'beff0_2array', subname)
allocate(ceff0_2array(-3+lowfil:lupfil+3,0:i), stat=istat)
call memocc(istat, ceff0_2array, 'ceff0_2array', subname)
allocate(eeff0_2array(lowfil:lupfil,0:i), stat=istat)
call memocc(istat, eeff0_2array, 'eeff0_2array', subname)
aeff0_2array=0.d0
beff0_2array=0.d0
ceff0_2array=0.d0
eeff0_2array=0.d0


allocate(aeff0_4array(-3+lowfil:lupfil+3,0:i), stat=istat)
call memocc(istat, aeff0_4array, 'aeff0_4array', subname)
allocate(beff0_4array(-3+lowfil:lupfil+3,0:i), stat=istat)
call memocc(istat, beff0_4array, 'beff0_4array', subname)
allocate(ceff0_4array(-3+lowfil:lupfil+3,0:i), stat=istat)
call memocc(istat, ceff0_4array, 'ceff0_4array', subname)
allocate(eeff0_4array(lowfil:lupfil,0:i), stat=istat)
call memocc(istat, eeff0_4array, 'eeff0_4array', subname)
aeff0_4array=0.d0
beff0_4array=0.d0
ceff0_4array=0.d0
eeff0_4array=0.d0


allocate(aeff0_2auxarray(-3+lowfil:lupfil+3,0:i), stat=istat)
call memocc(istat, aeff0_2auxarray, 'aeff0_2auxarray', subname)
allocate(beff0_2auxarray(-3+lowfil:lupfil+3,0:i), stat=istat)
call memocc(istat, beff0_2auxarray, 'beff0_2auxarray', subname)
allocate(ceff0_2auxarray(-3+lowfil:lupfil+3,0:i), stat=istat)
call memocc(istat, ceff0_2auxarray, 'ceff0_2auxarray', subname)
allocate(eeff0_2auxarray(-3+lowfil:lupfil+3,0:i), stat=istat)
call memocc(istat, eeff0_2auxarray, 'eeff0_2auxarray', subname)
aeff0_2auxarray=0.d0
beff0_2auxarray=0.d0
ceff0_2auxarray=0.d0
eeff0_2auxarray=0.d0


allocate(aeff0_4auxarray(-3+lowfil:lupfil+3,0:i), stat=istat)
call memocc(istat, aeff0_4auxarray, 'aeff0_4auxarray', subname)
allocate(beff0_4auxarray(-3+lowfil:lupfil+3,0:i), stat=istat)
call memocc(istat, beff0_4auxarray, 'beff0_4auxarray', subname)
allocate(ceff0_4auxarray(-3+lowfil:lupfil+3,0:i), stat=istat)
call memocc(istat, ceff0_4auxarray, 'ceff0_4auxarray', subname)
allocate(eeff0_4auxarray(-3+lowfil:lupfil+3,0:i), stat=istat)
call memocc(istat, eeff0_4auxarray, 'eeff0_4auxarray', subname)
aeff0_4auxarray=0.d0
beff0_4auxarray=0.d0
ceff0_4auxarray=0.d0
eeff0_4auxarray=0.d0


allocate(x2ya_c(0:n2,0:n1,0:n3), stat=istat)
call memocc(istat, x2ya_c, 'x2ya_c', subname)
allocate(x2yb_c(0:n2,0:n1,0:n3), stat=istat)
call memocc(istat, x2yb_c, 'x2yb_c', subname)
allocate(x2yc_c(0:n2,0:n1,0:n3), stat=istat)
call memocc(istat, x2yc_c, 'x2yc_c', subname)
allocate(x2ye_c(0:n2,0:n1,0:n3), stat=istat)
call memocc(istat, x2ye_c, 'x2ye_c', subname)
x2ya_c=0.d0
x2yb_c=0.d0
x2yc_c=0.d0
x2ye_c=0.d0







allocate(x2za_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, x2za_c, 'x2za_c', subname)
allocate(x2zb_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, x2zb_c, 'x2zb_c', subname)
allocate(x2zc_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, x2zc_c, 'x2zc_c', subname)
allocate(x2ze_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, x2ze_c, 'x2ze_c', subname)
x2za_c=0.d0
x2zb_c=0.d0
x2zc_c=0.d0
x2ze_c=0.d0





allocate(y2za_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, y2za_c, 'y2za_c', subname)
allocate(y2zb_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, y2zb_c, 'y2zb_c', subname)
allocate(y2zc_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, y2zc_c, 'y2zc_c', subname)
allocate(y2ze_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, y2ze_c, 'y2ze_c', subname)
y2za_c=0.d0
y2zb_c=0.d0
y2zc_c=0.d0
y2ze_c=0.d0





allocate(x4ya_c(0:n2,0:n1,0:n3), stat=istat)
call memocc(istat, x4ya_c, 'x4ya_c', subname)
allocate(x4yb_c(0:n2,0:n1,0:n3), stat=istat)
call memocc(istat, x4yb_c, 'x4yb_c', subname)
allocate(x4yc_c(0:n2,0:n1,0:n3), stat=istat)
call memocc(istat, x4yc_c, 'x4yc_c', subname)
allocate(x4ye_c(0:n2,0:n1,0:n3), stat=istat)
call memocc(istat, x4ye_c, 'x4ye_c', subname)
x4ya_c=0.d0
x4yb_c=0.d0
x4yc_c=0.d0
x4ye_c=0.d0

allocate(x4za_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, x4za_c, 'x4za_c', subname)
allocate(x4zb_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, x4zb_c, 'x4zb_c', subname)
allocate(x4zc_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, x4zc_c, 'x4zc_c', subname)
allocate(x4ze_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, x4ze_c, 'x4ze_c', subname)
x4za_c=0.d0
x4zb_c=0.d0
x4zc_c=0.d0
x4ze_c=0.d0

allocate(y4za_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, y4za_c, 'y4za_c', subname)
allocate(y4zb_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, y4zb_c, 'y4zb_c', subname)
allocate(y4zc_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, y4zc_c, 'y4zc_c', subname)
allocate(y4ze_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, y4ze_c, 'y4ze_c', subname)
y4za_c=0.d0
y4zb_c=0.d0
y4zc_c=0.d0
y4ze_c=0.d0



allocate(x2y2aa_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, x2y2aa_c, 'x2y2aa_c', subname)
allocate(x2y2ab_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, x2y2ab_c, 'x2y2ab_c', subname)
allocate(x2y2ac_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, x2y2ac_c, 'x2y2ac_c', subname)
allocate(x2y2ae_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, x2y2ae_c, 'x2y2ae_c', subname)
x2y2aa_c=0.d0
x2y2ab_c=0.d0
x2y2ac_c=0.d0
x2y2ae_c=0.d0

allocate(x2y2ba_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, x2y2ba_c, 'x2y2ba_c', subname)
allocate(x2y2bb_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, x2y2bb_c, 'x2y2bb_c', subname)
allocate(x2y2bc_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, x2y2bc_c, 'x2y2bc_c', subname)
allocate(x2y2be_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, x2y2be_c, 'x2y2be_c', subname)
x2y2ba_c=0.d0
x2y2bb_c=0.d0
x2y2bc_c=0.d0
x2y2be_c=0.d0

allocate(x2y2ca_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, x2y2ca_c, 'x2y2ca_c', subname)
allocate(x2y2cb_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, x2y2cb_c, 'x2y2cb_c', subname)
allocate(x2y2cc_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, x2y2cc_c, 'x2y2cc_c', subname)
allocate(x2y2ce_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, x2y2ce_c, 'x2y2ce_c', subname)
x2y2ca_c=0.d0
x2y2cb_c=0.d0
x2y2cc_c=0.d0
x2y2ce_c=0.d0

allocate(x2y2ea_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, x2y2ea_c, 'x2y2ea_c', subname)
allocate(x2y2eb_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, x2y2eb_c, 'x2y2eb_c', subname)
allocate(x2y2ec_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, x2y2ec_c, 'x2y2ec_c', subname)
allocate(x2y2ee_c(0:n3,0:n1,0:n2), stat=istat)
call memocc(istat, x2y2ee_c, 'x2y2ee_c', subname)
x2y2ea_c=0.d0
x2y2eb_c=0.d0
x2y2ec_c=0.d0
x2y2ee_c=0.d0



allocate(x2y2aa_f(7,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, x2y2aa_f, 'x2y2aa_f', subname)
allocate(x2y2ab_f(7,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, x2y2ab_f, 'x2y2ab_f', subname)
allocate(x2y2ac_f(7,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, x2y2ac_f, 'x2y2ac_f', subname)
allocate(x2y2ae_f(7,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, x2y2ae_f, 'x2y2ae_f', subname)
x2y2aa_f=0.d0
x2y2ab_f=0.d0
x2y2ac_f=0.d0
x2y2ae_f=0.d0

allocate(x2y2ba_f(7,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, x2y2ba_f, 'x2y2ba_f', subname)
allocate(x2y2bb_f(7,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, x2y2bb_f, 'x2y2bb_f', subname)
allocate(x2y2bc_f(7,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, x2y2bc_f, 'x2y2bc_f', subname)
allocate(x2y2be_f(7,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, x2y2be_f, 'x2y2be_f', subname)
x2y2ba_f=0.d0
x2y2bb_f=0.d0
x2y2bc_f=0.d0
x2y2be_f=0.d0

allocate(x2y2ca_f(7,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, x2y2ca_f, 'x2y2ca_f', subname)
allocate(x2y2cb_f(7,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, x2y2cb_f, 'x2y2cb_f', subname)
allocate(x2y2cc_f(7,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, x2y2cc_f, 'x2y2cc_f', subname)
allocate(x2y2ce_f(7,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, x2y2ce_f, 'x2y2ce_f', subname)
x2y2ca_f=0.d0
x2y2cb_f=0.d0
x2y2cc_f=0.d0
x2y2ce_f=0.d0


allocate(x2y2ea_f(7,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, x2y2ea_f, 'x2y2ea_f', subname)
allocate(x2y2eb_f(7,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, x2y2eb_f, 'x2y2eb_f', subname)
allocate(x2y2ec_f(7,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, x2y2ec_f, 'x2y2ec_f', subname)
allocate(x2y2ee_f(7,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, x2y2ee_f, 'x2y2ee_f', subname)
x2y2ea_f=0.d0
x2y2eb_f=0.d0
x2y2ec_f=0.d0
x2y2ee_f=0.d0



allocate(x2ya_f(3,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), stat=istat)
call memocc(istat, x2ya_f, 'x2ya_f', subname)
allocate(x2yb_f(4,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), stat=istat)
call memocc(istat, x2yb_f, 'x2yb_f', subname)
allocate(x2yc_f(3,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), stat=istat)
call memocc(istat, x2yc_f, 'x2yc_f', subname)
allocate(x2ye_f(4,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), stat=istat)
call memocc(istat, x2ye_f, 'x2ye_f', subname)
x2ya_f=0.d0
x2yb_f=0.d0
x2yc_f=0.d0
x2ye_f=0.d0


allocate(x2za_f(3,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, x2za_f, 'x2za_f', subname)
allocate(x2zb_f(4,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, x2zb_f, 'x2zb_f', subname)
allocate(x2zc_f(3,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, x2zc_f, 'x2zc_f', subname)
allocate(x2ze_f(4,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, x2ze_f, 'x2ze_f', subname)
x2za_f=0.d0
x2zb_f=0.d0
x2zc_f=0.d0
x2ze_f=0.d0


allocate(y2za_f(3,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, y2za_f, 'y2za_f', subname)
allocate(y2zb_f(4,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, y2zb_f, 'y2zb_f', subname)
allocate(y2zc_f(3,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, y2zc_f, 'y2zc_f', subname)
allocate(y2ze_f(4,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, y2ze_f, 'y2ze_f', subname)
y2za_f=0.d0
y2zb_f=0.d0
y2zc_f=0.d0
y2ze_f=0.d0



allocate(x2ya_f2(7,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), stat=istat)
call memocc(istat, x2ya_f2, 'x2ya_f2', subname)
allocate(x2yb_f2(7,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), stat=istat)
call memocc(istat, x2yb_f2, 'x2yb_f2', subname)
allocate(x2yc_f2(7,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), stat=istat)
call memocc(istat, x2yc_f2, 'x2yc_f2', subname)
allocate(x2ye_f2(7,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), stat=istat)
call memocc(istat, x2ye_f2, 'x2ye_f2', subname)
x2ya_f2=0.d0
x2yb_f2=0.d0
x2yc_f2=0.d0
x2ye_f2=0.d0




allocate(x4ya_f(3,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), stat=istat)
call memocc(istat, x4ya_f, 'x4ya_f', subname)
allocate(x4yb_f(4,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), stat=istat)
call memocc(istat, x4yb_f, 'x4yb_f', subname)
allocate(x4yc_f(3,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), stat=istat)
call memocc(istat, x4yc_f, 'x4yc_f', subname)
allocate(x4ye_f(4,nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), stat=istat)
call memocc(istat, x4ye_f, 'x4ye_f', subname)
x4ya_f=0.d0
x4yb_f=0.d0
x4yc_f=0.d0
x4ye_f=0.d0


allocate(x4za_f(3,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, x4za_f, 'x4za_f', subname)
allocate(x4zb_f(4,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, x4zb_f, 'x4zb_f', subname)
allocate(x4zc_f(3,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, x4zc_f, 'x4zc_f', subname)
allocate(x4ze_f(4,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, x4ze_f, 'x4ze_f', subname)
x4za_f=0.d0
x4zb_f=0.d0
x4zc_f=0.d0
x4ze_f=0.d0


allocate(y4za_f(3,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, y4za_f, 'y4za_f', subname)
allocate(y4zb_f(4,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, y4zb_f, 'y4zb_f', subname)
allocate(y4zc_f(3,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, y4zc_f, 'y4zc_f', subname)
allocate(y4ze_f(4,nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), stat=istat)
call memocc(istat, y4ze_f, 'y4ze_f', subname)
y4za_f=0.d0
y4zb_f=0.d0
y4zc_f=0.d0
y4ze_f=0.d0








aeff0=0.d0 ; beff0=0.d0 ; ceff0=0.d0 ; eeff0=0.0
aeff1=0.d0 ; beff1=0.d0 ; ceff1=0.d0 ; eeff1=0.0
aeff2=0.d0 ; beff2=0.d0 ; ceff2=0.d0 ; eeff2=0.0
aeff3=0.d0 ; beff3=0.d0 ; ceff3=0.d0 ; eeff3=0.0

aeff0_2=0.d0 ; beff0_2=0.d0 ; ceff0_2=0.d0 ; eeff0_2=0.0
aeff1_2=0.d0 ; beff1_2=0.d0 ; ceff1_2=0.d0 ; eeff1_2=0.0
aeff2_2=0.d0 ; beff2_2=0.d0 ; ceff2_2=0.d0 ; eeff2_2=0.0
aeff3_2=0.d0 ; beff3_2=0.d0 ; ceff3_2=0.d0 ; eeff3_2=0.0


  !!$!$omp parallel default(private) &
  !!$!$omp shared(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3) &
  !!$!$omp shared(cprecr,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,x_c,x_f,y_c,y_f)& 
  !!$!$omp shared(x_f1,x_f2,x_f3,a,b,c,e)
    !!!$omp do  

    do i1=0,n1
        x0=hgrid*(i1+offsetx)-rxyzConf(1)
        if(.not. WithKinetic) then
            call getFilterSextic(it, potentialPrefac, hgrid, x0, aeff0array(lowfil,i1), 'a')
            call getFilterSextic(it, potentialPrefac, hgrid, x0, beff0array(lowfil,i1), 'b')
            call getFilterSextic(it, potentialPrefac, hgrid, x0, ceff0array(lowfil,i1), 'c')
            call getFilterSextic(it, potentialPrefac, hgrid, x0, eeff0array(lowfil,i1), 'e')
        else
            call getEffectiveFilterSextic(it, potentialPrefac, hgrid, x0, aeff0array(lowfil,i1), 'a')
            call getEffectiveFilterSextic(it, potentialPrefac, hgrid, x0, beff0array(lowfil,i1), 'b')
            call getEffectiveFilterSextic(it, potentialPrefac, hgrid, x0, ceff0array(lowfil,i1), 'c')
            call getEffectiveFilterSextic(it, potentialPrefac, hgrid, x0, eeff0array(lowfil,i1), 'e')
        end if

        call getFilterQuartic(1.d0, hgrid, x0, aeff0_4auxarray(lowfil,i1), 'a')
        call getFilterQuartic(1.d0, hgrid, x0, beff0_4auxarray(lowfil,i1), 'b')
        call getFilterQuartic(1.d0, hgrid, x0, ceff0_4auxarray(lowfil,i1), 'c')
        call getFilterQuartic(1.d0, hgrid, x0, eeff0_4auxarray(lowfil,i1), 'e')

        call getFilterQuadratic(1.d0, hgrid, x0, aeff0_2auxarray(lowfil,i1), 'a')
        call getFilterQuadratic(1.d0, hgrid, x0, beff0_2auxarray(lowfil,i1), 'b')
        call getFilterQuadratic(1.d0, hgrid, x0, ceff0_2auxarray(lowfil,i1), 'c')
        call getFilterQuadratic(1.d0, hgrid, x0, eeff0_2auxarray(lowfil,i1), 'e')
    end do
!t1=mpi_wtime()
    do i3=0,n3
       do i2=0,n2
          if (ibyz_c(2,i2,i3)-ibyz_c(1,i2,i3).ge.4) then
             do i1=ibyz_c(1,i2,i3),ibyz_c(2,i2,i3)-4,4
                dyi0=0.0_wp 
                dyi1=0.0_wp 
                dyi2=0.0_wp 
                dyi3=0.0_wp 

                tt0a0_2=0.d0  ; tt0a0_4=0.d0
                tt0a1_2=0.d0  ; tt0a1_4=0.d0
                tt0a2_2=0.d0  ; tt0a2_4=0.d0
                tt0a3_2=0.d0  ; tt0a3_4=0.d0

                tt0b0_2=0.d0  ; tt0b0_4=0.d0
                tt0b1_2=0.d0  ; tt0b1_4=0.d0
                tt0b2_2=0.d0  ; tt0b2_4=0.d0
                tt0b3_2=0.d0  ; tt0b3_4=0.d0

                tt0c0_2=0.d0  ; tt0c0_4=0.d0
                tt0c1_2=0.d0  ; tt0c1_4=0.d0
                tt0c2_2=0.d0  ; tt0c2_4=0.d0
                tt0c3_2=0.d0  ; tt0c3_4=0.d0

                tt0e0_2=0.d0  ; tt0e0_4=0.d0
                tt0e1_2=0.d0  ; tt0e1_4=0.d0
                tt0e2_2=0.d0  ; tt0e2_4=0.d0
                tt0e3_2=0.d0  ; tt0e3_4=0.d0
  
                do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1+3,ibyz_c(2,i2,i3))
                   dyi0=dyi0 + xx_c(t,i2,i3)*aeff0array(t-i1-0,i1+0)
                   dyi1=dyi1 + xx_c(t,i2,i3)*aeff0array(t-i1-1,i1+1)
                   dyi2=dyi2 + xx_c(t,i2,i3)*aeff0array(t-i1-2,i1+2)
                   dyi3=dyi3 + xx_c(t,i2,i3)*aeff0array(t-i1-3,i1+3)

                   ! sss coefficients
                   tt0a0_2=tt0a0_2 + xx_c(t,i2,i3)*aeff0_2auxarray(t-i1-0,i1+0)
                   tt0a1_2=tt0a1_2 + xx_c(t,i2,i3)*aeff0_2auxarray(t-i1-1,i1+1)
                   tt0a2_2=tt0a2_2 + xx_c(t,i2,i3)*aeff0_2auxarray(t-i1-2,i1+2)
                   tt0a3_2=tt0a3_2 + xx_c(t,i2,i3)*aeff0_2auxarray(t-i1-3,i1+3)

                   tt0b0_2=tt0b0_2 + xx_c(t,i2,i3)*beff0_2auxarray(t-i1-0,i1+0)
                   tt0b1_2=tt0b1_2 + xx_c(t,i2,i3)*beff0_2auxarray(t-i1-1,i1+1)
                   tt0b2_2=tt0b2_2 + xx_c(t,i2,i3)*beff0_2auxarray(t-i1-2,i1+2)
                   tt0b3_2=tt0b3_2 + xx_c(t,i2,i3)*beff0_2auxarray(t-i1-3,i1+3)

                   tt0c0_2=tt0c0_2 + xx_c(t,i2,i3)*ceff0_2auxarray(t-i1-0,i1+0)
                   tt0c1_2=tt0c1_2 + xx_c(t,i2,i3)*ceff0_2auxarray(t-i1-1,i1+1)
                   tt0c2_2=tt0c2_2 + xx_c(t,i2,i3)*ceff0_2auxarray(t-i1-2,i1+2)
                   tt0c3_2=tt0c3_2 + xx_c(t,i2,i3)*ceff0_2auxarray(t-i1-3,i1+3)

                   tt0e0_2=tt0e0_2 + xx_c(t,i2,i3)*eeff0_2auxarray(t-i1-0,i1+0)
                   tt0e1_2=tt0e1_2 + xx_c(t,i2,i3)*eeff0_2auxarray(t-i1-1,i1+1)
                   tt0e2_2=tt0e2_2 + xx_c(t,i2,i3)*eeff0_2auxarray(t-i1-2,i1+2)
                   tt0e3_2=tt0e3_2 + xx_c(t,i2,i3)*eeff0_2auxarray(t-i1-3,i1+3)

                   tt0a0_4=tt0a0_4 + xx_c(t,i2,i3)*aeff0_4auxarray(t-i1-0,i1+0)
                   tt0a1_4=tt0a1_4 + xx_c(t,i2,i3)*aeff0_4auxarray(t-i1-1,i1+1)
                   tt0a2_4=tt0a2_4 + xx_c(t,i2,i3)*aeff0_4auxarray(t-i1-2,i1+2)
                   tt0a3_4=tt0a3_4 + xx_c(t,i2,i3)*aeff0_4auxarray(t-i1-3,i1+3)

                   tt0b0_4=tt0b0_4 + xx_c(t,i2,i3)*beff0_4auxarray(t-i1-0,i1+0)
                   tt0b1_4=tt0b1_4 + xx_c(t,i2,i3)*beff0_4auxarray(t-i1-1,i1+1)
                   tt0b2_4=tt0b2_4 + xx_c(t,i2,i3)*beff0_4auxarray(t-i1-2,i1+2)
                   tt0b3_4=tt0b3_4 + xx_c(t,i2,i3)*beff0_4auxarray(t-i1-3,i1+3)

                   tt0c0_4=tt0c0_4 + xx_c(t,i2,i3)*ceff0_4auxarray(t-i1-0,i1+0)
                   tt0c1_4=tt0c1_4 + xx_c(t,i2,i3)*ceff0_4auxarray(t-i1-1,i1+1)
                   tt0c2_4=tt0c2_4 + xx_c(t,i2,i3)*ceff0_4auxarray(t-i1-2,i1+2)
                   tt0c3_4=tt0c3_4 + xx_c(t,i2,i3)*ceff0_4auxarray(t-i1-3,i1+3)

                   tt0e0_4=tt0e0_4 + xx_c(t,i2,i3)*eeff0_4auxarray(t-i1-0,i1+0)
                   tt0e1_4=tt0e1_4 + xx_c(t,i2,i3)*eeff0_4auxarray(t-i1-1,i1+1)
                   tt0e2_4=tt0e2_4 + xx_c(t,i2,i3)*eeff0_4auxarray(t-i1-2,i1+2)
                   tt0e3_4=tt0e3_4 + xx_c(t,i2,i3)*eeff0_4auxarray(t-i1-3,i1+3)
                enddo
                y_c(i1+0,i2,i3)=dyi0+cprecr*xx_c(i1+0,i2,i3)
                y_c(i1+1,i2,i3)=dyi1+cprecr*xx_c(i1+1,i2,i3)
                y_c(i1+2,i2,i3)=dyi2+cprecr*xx_c(i1+2,i2,i3)
                y_c(i1+3,i2,i3)=dyi3+cprecr*xx_c(i1+3,i2,i3)

                x2ya_c(i2,i1+0,i3)=tt0a0_2
                x2ya_c(i2,i1+1,i3)=tt0a1_2
                x2ya_c(i2,i1+2,i3)=tt0a2_2
                x2ya_c(i2,i1+3,i3)=tt0a3_2
                x2za_c(i3,i1+0,i2)=tt0a0_2
                x2za_c(i3,i1+1,i2)=tt0a1_2
                x2za_c(i3,i1+2,i2)=tt0a2_2
                x2za_c(i3,i1+3,i2)=tt0a3_2

                x2yb_c(i2,i1+0,i3)=tt0b0_2
                x2yb_c(i2,i1+1,i3)=tt0b1_2
                x2yb_c(i2,i1+2,i3)=tt0b2_2
                x2yb_c(i2,i1+3,i3)=tt0b3_2
                x2zb_c(i3,i1+0,i2)=tt0b0_2
                x2zb_c(i3,i1+1,i2)=tt0b1_2
                x2zb_c(i3,i1+2,i2)=tt0b2_2
                x2zb_c(i3,i1+3,i2)=tt0b3_2

                x2yc_c(i2,i1+0,i3)=tt0c0_2
                x2yc_c(i2,i1+1,i3)=tt0c1_2
                x2yc_c(i2,i1+2,i3)=tt0c2_2
                x2yc_c(i2,i1+3,i3)=tt0c3_2
                x2zc_c(i3,i1+0,i2)=tt0c0_2
                x2zc_c(i3,i1+1,i2)=tt0c1_2
                x2zc_c(i3,i1+2,i2)=tt0c2_2
                x2zc_c(i3,i1+3,i2)=tt0c3_2

                x2ye_c(i2,i1+0,i3)=tt0e0_2
                x2ye_c(i2,i1+1,i3)=tt0e1_2
                x2ye_c(i2,i1+2,i3)=tt0e2_2
                x2ye_c(i2,i1+3,i3)=tt0e3_2
                x2ze_c(i3,i1+0,i2)=tt0e0_2
                x2ze_c(i3,i1+1,i2)=tt0e1_2
                x2ze_c(i3,i1+2,i2)=tt0e2_2
                x2ze_c(i3,i1+3,i2)=tt0e3_2

                x4ya_c(i2,i1+0,i3)=tt0a0_4
                x4ya_c(i2,i1+1,i3)=tt0a1_4
                x4ya_c(i2,i1+2,i3)=tt0a2_4
                x4ya_c(i2,i1+3,i3)=tt0a3_4
                x4za_c(i3,i1+0,i2)=tt0a0_4
                x4za_c(i3,i1+1,i2)=tt0a1_4
                x4za_c(i3,i1+2,i2)=tt0a2_4
                x4za_c(i3,i1+3,i2)=tt0a3_4

                x4yb_c(i2,i1+0,i3)=tt0b0_4
                x4yb_c(i2,i1+1,i3)=tt0b1_4
                x4yb_c(i2,i1+2,i3)=tt0b2_4
                x4yb_c(i2,i1+3,i3)=tt0b3_4
                x4zb_c(i3,i1+0,i2)=tt0b0_4
                x4zb_c(i3,i1+1,i2)=tt0b1_4
                x4zb_c(i3,i1+2,i2)=tt0b2_4
                x4zb_c(i3,i1+3,i2)=tt0b3_4

                x4yc_c(i2,i1+0,i3)=tt0c0_4
                x4yc_c(i2,i1+1,i3)=tt0c1_4
                x4yc_c(i2,i1+2,i3)=tt0c2_4
                x4yc_c(i2,i1+3,i3)=tt0c3_4
                x4zc_c(i3,i1+0,i2)=tt0c0_4
                x4zc_c(i3,i1+1,i2)=tt0c1_4
                x4zc_c(i3,i1+2,i2)=tt0c2_4
                x4zc_c(i3,i1+3,i2)=tt0c3_4

                x4ye_c(i2,i1+0,i3)=tt0e0_4
                x4ye_c(i2,i1+1,i3)=tt0e1_4
                x4ye_c(i2,i1+2,i3)=tt0e2_4
                x4ye_c(i2,i1+3,i3)=tt0e3_4
                x4ze_c(i3,i1+0,i2)=tt0e0_4
                x4ze_c(i3,i1+1,i2)=tt0e1_4
                x4ze_c(i3,i1+2,i2)=tt0e2_4
                x4ze_c(i3,i1+3,i2)=tt0e3_4
             enddo
             icur=i1
          else
             icur=ibyz_c(1,i2,i3)
          endif
  
          do i1=icur,ibyz_c(2,i2,i3)
             dyi=0.0_wp 
             tt0a0_2=0.d0 ; tt0a0_4=0.d0
             tt0b0_2=0.d0 ; tt0b0_4=0.d0
             tt0c0_2=0.d0 ; tt0c0_4=0.d0
             tt0e0_2=0.d0 ; tt0e0_4=0.d0
             ! Get the effective a-filters for the x dimension
             do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_c(2,i2,i3))
                dyi=dyi + xx_c(t,i2,i3)*aeff0array(t-i1,i1)
                ! sss coefficients
                tt0a0_2=tt0a0_2 + xx_c(t,i2,i3)*aeff0_2auxarray(t-i1,i1)
                tt0b0_2=tt0b0_2 + xx_c(t,i2,i3)*beff0_2auxarray(t-i1,i1)
                tt0c0_2=tt0c0_2 + xx_c(t,i2,i3)*ceff0_2auxarray(t-i1,i1)
                tt0e0_2=tt0e0_2 + xx_c(t,i2,i3)*eeff0_2auxarray(t-i1,i1)
                tt0a0_4=tt0a0_4 + xx_c(t,i2,i3)*aeff0_4auxarray(t-i1,i1)
                tt0b0_4=tt0b0_4 + xx_c(t,i2,i3)*beff0_4auxarray(t-i1,i1)
                tt0c0_4=tt0c0_4 + xx_c(t,i2,i3)*ceff0_4auxarray(t-i1,i1)
                tt0e0_4=tt0e0_4 + xx_c(t,i2,i3)*eeff0_4auxarray(t-i1,i1)
             enddo
             y_c(i1,i2,i3)=dyi+cprecr*xx_c(i1,i2,i3)

             x2ya_c(i2,i1,i3)=tt0a0_2
             x2za_c(i3,i1,i2)=tt0a0_2

             x2yb_c(i2,i1,i3)=tt0b0_2
             x2zb_c(i3,i1,i2)=tt0b0_2

             x2yc_c(i2,i1,i3)=tt0c0_2
             x2zc_c(i3,i1,i2)=tt0c0_2

             x2yc_c(i2,i1,i3)=tt0c0_2
             x2zc_c(i3,i1,i2)=tt0c0_2

             x4ya_c(i2,i1,i3)=tt0a0_4
             x4za_c(i3,i1,i2)=tt0a0_4

             x4yb_c(i2,i1,i3)=tt0b0_4
             x4zb_c(i3,i1,i2)=tt0b0_4

             x4yc_c(i2,i1,i3)=tt0c0_4
             x4zc_c(i3,i1,i2)=tt0c0_4

             x4yc_c(i2,i1,i3)=tt0c0_4
             x4zc_c(i3,i1,i2)=tt0c0_4
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
             dyi=0.0_wp
             tt0=0.0_wp
             do t=max(ibyz_f(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_f(2,i2,i3))
                dyi=dyi + xx_f1(t,i2,i3)*beff0array(t-i1,i1)
             enddo
             y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi
          enddo
  
           if (ibyz_c(2,i2,i3)-ibyz_c(1,i2,i3).ge.4) then
             do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)-4,4
                dyi0=0.0_wp 
                dyi1=0.0_wp 
                dyi2=0.0_wp 
                dyi3=0.0_wp 
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
             dyi=0.0_wp 
             tt0=0.0_wp 
             do t=max(ibyz_c(1,i2,i3),lowfil+i1),min(lupfil+i1,ibyz_c(2,i2,i3))
                dyi=dyi + xx_c(t,i2,i3)*ceff0array(t-i1,i1)
             enddo
             y_f(1,i1,i2,i3)=dyi
          enddo
       enddo
    enddo
    !!!$omp enddo
  
  
  
  
    ! wavelet part
  
    !!!$omp do
    do i3=nfl3,nfu3
       do i2=nfl2,nfu2
          do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
             t112=0.0_wp;t121=0.0_wp;t122=0.0_wp;t212=0.0_wp;t221=0.0_wp;t222=0.0_wp;t211=0.0_wp 
             tt112=0.0_wp;tt121=0.0_wp;tt122=0.0_wp;tt212=0.0_wp;tt221=0.0_wp;tt222=0.0_wp;tt211=0.0_wp 
             tt1a0_2=0.d0 ; tt1a0_4=0.d0
             tt1b0_2=0.d0 ; tt1b0_4=0.d0
             tt1c0_2=0.d0 ; tt1c0_4=0.d0
             tt1e0_2=0.d0 ; tt1e0_4=0.d0

             tt2a0_2=0.d0 ; tt2a0_4=0.d0
             tt2b0_2=0.d0 ; tt2b0_4=0.d0
             tt2c0_2=0.d0 ; tt2c0_4=0.d0
             tt2e0_2=0.d0 ; tt2e0_4=0.d0

             tt3a0_2=0.d0 ; tt3a0_4=0.d0
             tt3b0_2=0.d0 ; tt3b0_4=0.d0
             tt3c0_2=0.d0 ; tt3c0_4=0.d0
             tt3e0_2=0.d0 ; tt3e0_4=0.d0

             tt4a0_2=0.d0 ; tt4a0_4=0.d0
             tt4b0_2=0.d0 ; tt4b0_4=0.d0
             tt4c0_2=0.d0 ; tt4c0_4=0.d0
             tt4e0_2=0.d0 ; tt4e0_4=0.d0

             tt5a0_2=0.d0 ; tt5a0_4=0.d0
             tt5b0_2=0.d0 ; tt5b0_4=0.d0
             tt5c0_2=0.d0 ; tt5c0_4=0.d0
             tt5e0_2=0.d0 ; tt5e0_4=0.d0

             tt6a0_2=0.d0 ; tt6a0_4=0.d0
             tt6b0_2=0.d0 ; tt6b0_4=0.d0
             tt6c0_2=0.d0 ; tt6c0_4=0.d0
             tt6e0_2=0.d0 ; tt6e0_4=0.d0

             tt7a0_2=0.d0 ; tt7a0_4=0.d0
             tt7b0_2=0.d0 ; tt7b0_4=0.d0
             tt7c0_2=0.d0 ; tt7c0_4=0.d0
             tt7e0_2=0.d0 ; tt7e0_4=0.d0
             do l=max(nfl1-i1,lowfil),min(lupfil,nfu1-i1)
                t112=t112 + xx_f(4,i1+l,i2,i3)*aeff0array(l,i1) + xx_f(5,i1+l,i2,i3)*beff0array(l,i1)
                t121=t121 + xx_f(2,i1+l,i2,i3)*aeff0array(l,i1) + xx_f(3,i1+l,i2,i3)*beff0array(l,i1)
                t122=t122 + xx_f(6,i1+l,i2,i3)*aeff0array(l,i1) + xx_f(7,i1+l,i2,i3)*beff0array(l,i1)
                t212=t212 + xx_f(4,i1+l,i2,i3)*ceff0array(l,i1) + xx_f(5,i1+l,i2,i3)*eeff0array(l,i1)
                t221=t221 + xx_f(2,i1+l,i2,i3)*ceff0array(l,i1) + xx_f(3,i1+l,i2,i3)*eeff0array(l,i1)
                t222=t222 + xx_f(6,i1+l,i2,i3)*ceff0array(l,i1) + xx_f(7,i1+l,i2,i3)*eeff0array(l,i1)
                t211=t211 + xx_f(1,i1+l,i2,i3)*eeff0array(l,i1)
                ! dss coefficients
                tt1b0_2=tt1b0_2 + xx_f(1,i1+l,i2,i3)*beff0_2auxarray(l,i1)
                tt1e0_2=tt1e0_2 + xx_f(1,i1+l,i2,i3)*eeff0_2auxarray(l,i1)
                tt1b0_4=tt1b0_4 + xx_f(1,i1+l,i2,i3)*beff0_4auxarray(l,i1)
                tt1e0_4=tt1e0_4 + xx_f(1,i1+l,i2,i3)*eeff0_4auxarray(l,i1)
                ! sds coefficients
                tt2a0_2=tt2a0_2 + xx_f(2,i1+l,i2,i3)*aeff0_2auxarray(l,i1)
                tt2c0_2=tt2c0_2 + xx_f(2,i1+l,i2,i3)*ceff0_2auxarray(l,i1)
                tt2a0_4=tt2a0_4 + xx_f(2,i1+l,i2,i3)*aeff0_4auxarray(l,i1)
                tt2c0_4=tt2c0_4 + xx_f(2,i1+l,i2,i3)*ceff0_4auxarray(l,i1)
                ! dds coefficients
                tt3b0_2=tt3b0_2 + xx_f(3,i1+l,i2,i3)*beff0_2auxarray(l,i1)
                tt3e0_2=tt3e0_2 + xx_f(3,i1+l,i2,i3)*eeff0_2auxarray(l,i1)
                tt3b0_4=tt3b0_4 + xx_f(3,i1+l,i2,i3)*beff0_4auxarray(l,i1)
                tt3e0_4=tt3e0_4 + xx_f(3,i1+l,i2,i3)*eeff0_4auxarray(l,i1)
                ! ssd coefficients
                tt4a0_2=tt4a0_2 + xx_f(4,i1+l,i2,i3)*aeff0_2auxarray(l,i1)
                tt4c0_2=tt4c0_2 + xx_f(4,i1+l,i2,i3)*ceff0_2auxarray(l,i1)
                tt4a0_4=tt4a0_4 + xx_f(4,i1+l,i2,i3)*aeff0_4auxarray(l,i1)
                tt4c0_4=tt4c0_4 + xx_f(4,i1+l,i2,i3)*ceff0_4auxarray(l,i1)
                ! dsd coefficients
                tt5b0_2=tt5b0_2 + xx_f(5,i1+l,i2,i3)*beff0_2auxarray(l,i1)
                tt5e0_2=tt5e0_2 + xx_f(5,i1+l,i2,i3)*eeff0_2auxarray(l,i1)
                tt5b0_4=tt5b0_4 + xx_f(5,i1+l,i2,i3)*beff0_4auxarray(l,i1)
                tt5e0_4=tt5e0_4 + xx_f(5,i1+l,i2,i3)*eeff0_4auxarray(l,i1)
                ! sdd coefficients
                tt6a0_2=tt6a0_2 + xx_f(6,i1+l,i2,i3)*aeff0_2auxarray(l,i1)
                tt6c0_2=tt6c0_2 + xx_f(6,i1+l,i2,i3)*ceff0_2auxarray(l,i1)
                tt6a0_4=tt6a0_4 + xx_f(6,i1+l,i2,i3)*aeff0_4auxarray(l,i1)
                tt6c0_4=tt6c0_4 + xx_f(6,i1+l,i2,i3)*ceff0_4auxarray(l,i1)
                ! ddd coefficients
                tt7b0_2=tt7b0_2 + xx_f(7,i1+l,i2,i3)*beff0_2auxarray(l,i1)
                tt7e0_2=tt7e0_2 + xx_f(7,i1+l,i2,i3)*eeff0_2auxarray(l,i1)
                tt7b0_4=tt7b0_4 + xx_f(7,i1+l,i2,i3)*beff0_4auxarray(l,i1)
                tt7e0_4=tt7e0_4 + xx_f(7,i1+l,i2,i3)*eeff0_4auxarray(l,i1)
             enddo
             y_f(4,i1,i2,i3)=t112+cprecr*xx_f(4,i1,i2,i3)
             y_f(2,i1,i2,i3)=t121+cprecr*xx_f(2,i1,i2,i3)
             y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+t211+cprecr*xx_f(1,i1,i2,i3)
             y_f(6,i1,i2,i3)=t122+cprecr*xx_f(6,i1,i2,i3)
             y_f(5,i1,i2,i3)=t212+cprecr*xx_f(5,i1,i2,i3)
             y_f(3,i1,i2,i3)=t221+cprecr*xx_f(3,i1,i2,i3)
             y_f(7,i1,i2,i3)=t222+cprecr*xx_f(7,i1,i2,i3)
             ! dss coefficients
             x2ya_f2(1,i2,i1,i3)=tt1a0_2
             x2yb_f2(1,i2,i1,i3)=tt1b0_2
             x2yc_f2(1,i2,i1,i3)=tt1c0_2
             x2ye_f2(1,i2,i1,i3)=tt1e0_2
             x2yb_f(1,i2,i1,i3)=tt1b0_2
             x2ye_f(1,i2,i1,i3)=tt1e0_2
             x2zb_f(1,i3,i1,i2)=tt1b0_2
             x2ze_f(1,i3,i1,i2)=tt1e0_2
             x4yb_f(1,i2,i1,i3)=tt1b0_4
             x4ye_f(1,i2,i1,i3)=tt1e0_4
             x4zb_f(1,i3,i1,i2)=tt1b0_4
             x4ze_f(1,i3,i1,i2)=tt1e0_4
             ! sds coefficients
             x2ya_f2(2,i2,i1,i3)=tt2a0_2
             x2yb_f2(2,i2,i1,i3)=tt2b0_2
             x2yc_f2(2,i2,i1,i3)=tt2c0_2
             x2ye_f2(2,i2,i1,i3)=tt2e0_2
             x2ya_f(1,i2,i1,i3)=tt2a0_2
             x2yc_f(1,i2,i1,i3)=tt2c0_2
             x2za_f(1,i3,i1,i2)=tt2a0_2
             x2zc_f(1,i3,i1,i2)=tt2c0_2
             x4ya_f(1,i2,i1,i3)=tt2a0_4
             x4yc_f(1,i2,i1,i3)=tt2c0_4
             x4za_f(1,i3,i1,i2)=tt2a0_4
             x4zc_f(1,i3,i1,i2)=tt2c0_4
             ! dds coefficients
             x2ya_f2(3,i2,i1,i3)=tt3a0_2
             x2yb_f2(3,i2,i1,i3)=tt3b0_2
             x2yc_f2(3,i2,i1,i3)=tt3c0_2
             x2ye_f2(3,i2,i1,i3)=tt3e0_2
             x2yb_f(2,i2,i1,i3)=tt3b0_2
             x2ye_f(2,i2,i1,i3)=tt3e0_2
             x2zb_f(2,i3,i1,i2)=tt3b0_2
             x2ze_f(2,i3,i1,i2)=tt3e0_2
             x4yb_f(2,i2,i1,i3)=tt3b0_4
             x4ye_f(2,i2,i1,i3)=tt3e0_4
             x4zb_f(2,i3,i1,i2)=tt3b0_4
             x4ze_f(2,i3,i1,i2)=tt3e0_4
             ! ssd coefficients
             x2ya_f2(4,i2,i1,i3)=tt4a0_2
             x2yb_f2(4,i2,i1,i3)=tt4b0_2
             x2yc_f2(4,i2,i1,i3)=tt4c0_2
             x2ye_f2(4,i2,i1,i3)=tt4e0_2
             x2ya_f(2,i2,i1,i3)=tt4a0_2
             x2yc_f(2,i2,i1,i3)=tt4c0_2
             x2za_f(2,i3,i1,i2)=tt4a0_2
             x2zc_f(2,i3,i1,i2)=tt4c0_2
             x4ya_f(2,i2,i1,i3)=tt4a0_4
             x4yc_f(2,i2,i1,i3)=tt4c0_4
             x4za_f(2,i3,i1,i2)=tt4a0_4
             x4zc_f(2,i3,i1,i2)=tt4c0_4
             ! dsd coefficients
             x2ya_f2(5,i2,i1,i3)=tt5a0_2
             x2yb_f2(5,i2,i1,i3)=tt5b0_2
             x2yc_f2(5,i2,i1,i3)=tt5c0_2
             x2ye_f2(5,i2,i1,i3)=tt5e0_2
             x2yb_f(3,i2,i1,i3)=tt5b0_2
             x2ye_f(3,i2,i1,i3)=tt5e0_2
             x2zb_f(3,i3,i1,i2)=tt5b0_2
             x2ze_f(3,i3,i1,i2)=tt5e0_2
             x4yb_f(3,i2,i1,i3)=tt5b0_4
             x4ye_f(3,i2,i1,i3)=tt5e0_4
             x4zb_f(3,i3,i1,i2)=tt5b0_4
             x4ze_f(3,i3,i1,i2)=tt5e0_4
             ! sdd coefficients
             x2ya_f2(6,i2,i1,i3)=tt6a0_2
             x2yb_f2(6,i2,i1,i3)=tt6b0_2
             x2yc_f2(6,i2,i1,i3)=tt6c0_2
             x2ye_f2(6,i2,i1,i3)=tt6e0_2
             x2ya_f(3,i2,i1,i3)=tt6a0_2
             x2yc_f(3,i2,i1,i3)=tt6c0_2
             x2za_f(3,i3,i1,i2)=tt6a0_2
             x2zc_f(3,i3,i1,i2)=tt6c0_2
             x4ya_f(3,i2,i1,i3)=tt6a0_4
             x4yc_f(3,i2,i1,i3)=tt6c0_4
             x4za_f(3,i3,i1,i2)=tt6a0_4
             x4zc_f(3,i3,i1,i2)=tt6c0_4
             ! sdd coefficients
             x2ya_f2(7,i2,i1,i3)=tt7a0_2
             x2yb_f2(7,i2,i1,i3)=tt7b0_2
             x2yc_f2(7,i2,i1,i3)=tt7c0_2
             x2ye_f2(7,i2,i1,i3)=tt7e0_2
             x2yb_f(4,i2,i1,i3)=tt7b0_2
             x2ye_f(4,i2,i1,i3)=tt7e0_2
             x2zb_f(4,i3,i1,i2)=tt7b0_2
             x2ze_f(4,i3,i1,i2)=tt7e0_2
             x4yb_f(4,i2,i1,i3)=tt7b0_4
             x4ye_f(4,i2,i1,i3)=tt7e0_4
             x4zb_f(4,i3,i1,i2)=tt7b0_4
             x4ze_f(4,i3,i1,i2)=tt7e0_4
          enddo
       enddo
    enddo
    !!!$omp enddo


  
  

  
    do i2=0,n2
        y0=hgrid*(i2+offsety)-rxyzConf(2)
        if(.not. withKinetic) then
            call getFilterSextic(it, potentialPrefac, hgrid, y0, aeff0array(lowfil,i2), 'a')
            call getFilterSextic(it, potentialPrefac, hgrid, y0, beff0array(lowfil,i2), 'b')
            call getFilterSextic(it, potentialPrefac, hgrid, y0, ceff0array(lowfil,i2), 'c')
            call getFilterSextic(it, potentialPrefac, hgrid, y0, eeff0array(lowfil,i2), 'e')
        else
            call getEffectiveFilterSextic(it, potentialPrefac, hgrid, y0, aeff0array(lowfil,i2), 'a')
            call getEffectiveFilterSextic(it, potentialPrefac, hgrid, y0, beff0array(lowfil,i2), 'b')
            call getEffectiveFilterSextic(it, potentialPrefac, hgrid, y0, ceff0array(lowfil,i2), 'c')
            call getEffectiveFilterSextic(it, potentialPrefac, hgrid, y0, eeff0array(lowfil,i2), 'e')
        end if

        call getFilterQuartic(potentialPrefac, hgrid, y0, aeff0_4array(lowfil,i2), 'a')
        call getFilterQuartic(potentialPrefac, hgrid, y0, beff0_4array(lowfil,i2), 'b')
        call getFilterQuartic(potentialPrefac, hgrid, y0, ceff0_4array(lowfil,i2), 'c')
        call getFilterQuartic(potentialPrefac, hgrid, y0, eeff0_4array(lowfil,i2), 'e')

        call getFilterQuadratic(potentialPrefac, hgrid, y0, aeff0_2array(lowfil,i2), 'a')
        call getFilterQuadratic(potentialPrefac, hgrid, y0, beff0_2array(lowfil,i2), 'b')
        call getFilterQuadratic(potentialPrefac, hgrid, y0, ceff0_2array(lowfil,i2), 'c')
        call getFilterQuadratic(potentialPrefac, hgrid, y0, eeff0_2array(lowfil,i2), 'e')

        call getFilterQuartic(1.d0, hgrid, y0, aeff0_4auxarray(lowfil,i2), 'a')
        call getFilterQuartic(1.d0, hgrid, y0, beff0_4auxarray(lowfil,i2), 'b')
        call getFilterQuartic(1.d0, hgrid, y0, ceff0_4auxarray(lowfil,i2), 'c')
        call getFilterQuartic(1.d0, hgrid, y0, eeff0_4auxarray(lowfil,i2), 'e')

        call getFilterQuadratic(1.d0, hgrid, y0, aeff0_2auxarray(lowfil,i2), 'a')
        call getFilterQuadratic(1.d0, hgrid, y0, beff0_2auxarray(lowfil,i2), 'b')
        call getFilterQuadratic(1.d0, hgrid, y0, ceff0_2auxarray(lowfil,i2), 'c')
        call getFilterQuadratic(1.d0, hgrid, y0, eeff0_2auxarray(lowfil,i2), 'e')
    end do
  
  
    ! + (1/2) d^2/dy^2
    !!!$omp do
    do i3=0,n3
       do i1=0,n1
          if (ibxz_c(2,i1,i3)-ibxz_c(1,i1,i3).ge.4) then
             do i2=ibxz_c(1,i1,i3),ibxz_c(2,i1,i3)-4,4
                dyi0=0.0_wp 
                dyi1=0.0_wp 
                dyi2=0.0_wp 
                dyi3=0.0_wp 

                tt0a0_2=0.d0 ; tt0a0_4=0.d0
                tt0a1_2=0.d0 ; tt0a1_4=0.d0
                tt0a2_2=0.d0 ; tt0a2_4=0.d0
                tt0a3_2=0.d0 ; tt0a3_4=0.d0

                tt0b0_2=0.d0 ; tt0b0_4=0.d0
                tt0b1_2=0.d0 ; tt0b1_4=0.d0
                tt0b2_2=0.d0 ; tt0b2_4=0.d0
                tt0b3_2=0.d0 ; tt0b3_4=0.d0

                tt0c0_2=0.d0 ; tt0c0_4=0.d0
                tt0c1_2=0.d0 ; tt0c1_4=0.d0
                tt0c2_2=0.d0 ; tt0c2_4=0.d0
                tt0c3_2=0.d0 ; tt0c3_4=0.d0

                tt0e0_2=0.d0 ; tt0e0_4=0.d0
                tt0e1_2=0.d0 ; tt0e1_4=0.d0
                tt0e2_2=0.d0 ; tt0e2_4=0.d0
                tt0e3_2=0.d0 ; tt0e3_4=0.d0

                ttaa0=0.d0 ; ttba0=0.d0 ; ttca0=0.d0 ; ttea0=0.d0
                ttaa1=0.d0 ; ttba1=0.d0 ; ttca1=0.d0 ; ttea1=0.d0
                ttaa2=0.d0 ; ttba2=0.d0 ; ttca2=0.d0 ; ttea2=0.d0
                ttaa3=0.d0 ; ttba3=0.d0 ; ttca3=0.d0 ; ttea3=0.d0

                ttab0=0.d0 ; ttbb0=0.d0 ; ttcb0=0.d0 ; tteb0=0.d0
                ttab1=0.d0 ; ttbb1=0.d0 ; ttcb1=0.d0 ; tteb1=0.d0
                ttab2=0.d0 ; ttbb2=0.d0 ; ttcb2=0.d0 ; tteb2=0.d0
                ttab3=0.d0 ; ttbb3=0.d0 ; ttcb3=0.d0 ; tteb3=0.d0

                ttac0=0.d0 ; ttbc0=0.d0 ; ttcc0=0.d0 ; ttec0=0.d0
                ttac1=0.d0 ; ttbc1=0.d0 ; ttcc1=0.d0 ; ttec1=0.d0
                ttac2=0.d0 ; ttbc2=0.d0 ; ttcc2=0.d0 ; ttec2=0.d0
                ttac3=0.d0 ; ttbc3=0.d0 ; ttcc3=0.d0 ; ttec3=0.d0

                ttae0=0.d0 ; ttbe0=0.d0 ; ttce0=0.d0 ; ttee0=0.d0
                ttae1=0.d0 ; ttbe1=0.d0 ; ttce1=0.d0 ; ttee1=0.d0
                ttae2=0.d0 ; ttbe2=0.d0 ; ttce2=0.d0 ; ttee2=0.d0
                ttae3=0.d0 ; ttbe3=0.d0 ; ttce3=0.d0 ; ttee3=0.d0

                   
                do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_c(2,i1,i3))
                   dyi0=dyi0 + xy_c(t,i1,i3)*aeff0array(t-i2-0,i2+0) + 3.d0*x4ya_c(t,i1,i3)*aeff0_2array(t-i2-0,i2+0) + &
                               3.d0*x2ya_c(t,i1,i3)*aeff0_4array(t-i2-0,i2+0)
                   dyi1=dyi1 + xy_c(t,i1,i3)*aeff0array(t-i2-1,i2+1) + 3.d0*x4ya_c(t,i1,i3)*aeff0_2array(t-i2-1,i2+1) + &
                               3.d0*x2ya_c(t,i1,i3)*aeff0_4array(t-i2-1,i2+1)
                   dyi2=dyi2 + xy_c(t,i1,i3)*aeff0array(t-i2-2,i2+2) + 3.d0*x4ya_c(t,i1,i3)*aeff0_2array(t-i2-2,i2+2) + &
                               3.d0*x2ya_c(t,i1,i3)*aeff0_4array(t-i2-2,i2+2)
                   dyi3=dyi3 + xy_c(t,i1,i3)*aeff0array(t-i2-3,i2+3) + 3.d0*x4ya_c(t,i1,i3)*aeff0_2array(t-i2-3,i2+3) + &
                               3.d0*x2ya_c(t,i1,i3)*aeff0_4array(t-i2-3,i2+3)

                   ! sss coefficients
                   tt0a0_2=tt0a0_2 + xy_c(t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)
                   tt0a1_2=tt0a1_2 + xy_c(t,i1,i3)*aeff0_2auxarray(t-i2-1,i2+1)
                   tt0a2_2=tt0a2_2 + xy_c(t,i1,i3)*aeff0_2auxarray(t-i2-2,i2+2)
                   tt0a3_2=tt0a3_2 + xy_c(t,i1,i3)*aeff0_2auxarray(t-i2-3,i2+3)

                   tt0b0_2=tt0b0_2 + xy_c(t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)
                   tt0b1_2=tt0b1_2 + xy_c(t,i1,i3)*beff0_2auxarray(t-i2-1,i2+1)
                   tt0b2_2=tt0b2_2 + xy_c(t,i1,i3)*beff0_2auxarray(t-i2-2,i2+2)
                   tt0b3_2=tt0b3_2 + xy_c(t,i1,i3)*beff0_2auxarray(t-i2-3,i2+3)

                   tt0c0_2=tt0c0_2 + xy_c(t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)
                   tt0c1_2=tt0c1_2 + xy_c(t,i1,i3)*ceff0_2auxarray(t-i2-1,i2+1)
                   tt0c2_2=tt0c2_2 + xy_c(t,i1,i3)*ceff0_2auxarray(t-i2-2,i2+2)
                   tt0c3_2=tt0c3_2 + xy_c(t,i1,i3)*ceff0_2auxarray(t-i2-3,i2+3)

                   tt0e0_2=tt0e0_2 + xy_c(t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)
                   tt0e1_2=tt0e1_2 + xy_c(t,i1,i3)*eeff0_2auxarray(t-i2-1,i2+1)
                   tt0e2_2=tt0e2_2 + xy_c(t,i1,i3)*eeff0_2auxarray(t-i2-2,i2+2)
                   tt0e3_2=tt0e3_2 + xy_c(t,i1,i3)*eeff0_2auxarray(t-i2-3,i2+3)

                   tt0a0_4=tt0a0_4 + xy_c(t,i1,i3)*aeff0_4auxarray(t-i2-0,i2+0)
                   tt0a1_4=tt0a1_4 + xy_c(t,i1,i3)*aeff0_4auxarray(t-i2-1,i2+1)
                   tt0a2_4=tt0a2_4 + xy_c(t,i1,i3)*aeff0_4auxarray(t-i2-2,i2+2)
                   tt0a3_4=tt0a3_4 + xy_c(t,i1,i3)*aeff0_4auxarray(t-i2-3,i2+3)

                   tt0b0_4=tt0b0_4 + xy_c(t,i1,i3)*beff0_4auxarray(t-i2-0,i2+0)
                   tt0b1_4=tt0b1_4 + xy_c(t,i1,i3)*beff0_4auxarray(t-i2-1,i2+1)
                   tt0b2_4=tt0b2_4 + xy_c(t,i1,i3)*beff0_4auxarray(t-i2-2,i2+2)
                   tt0b3_4=tt0b3_4 + xy_c(t,i1,i3)*beff0_4auxarray(t-i2-3,i2+3)

                   tt0c0_4=tt0c0_4 + xy_c(t,i1,i3)*ceff0_4auxarray(t-i2-0,i2+0)
                   tt0c1_4=tt0c1_4 + xy_c(t,i1,i3)*ceff0_4auxarray(t-i2-1,i2+1)
                   tt0c2_4=tt0c2_4 + xy_c(t,i1,i3)*ceff0_4auxarray(t-i2-2,i2+2)
                   tt0c3_4=tt0c3_4 + xy_c(t,i1,i3)*ceff0_4auxarray(t-i2-3,i2+3)

                   tt0e0_4=tt0e0_4 + xy_c(t,i1,i3)*eeff0_4auxarray(t-i2-0,i2+0)
                   tt0e1_4=tt0e1_4 + xy_c(t,i1,i3)*eeff0_4auxarray(t-i2-1,i2+1)
                   tt0e2_4=tt0e2_4 + xy_c(t,i1,i3)*eeff0_4auxarray(t-i2-2,i2+2)
                   tt0e3_4=tt0e3_4 + xy_c(t,i1,i3)*eeff0_4auxarray(t-i2-3,i2+3)

                   ! x^2y^2, a coefficients
                   ttaa0=ttaa0 + x2ya_c(t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)
                   ttaa1=ttaa1 + x2ya_c(t,i1,i3)*aeff0_2auxarray(t-i2-1,i2+1)
                   ttaa2=ttaa2 + x2ya_c(t,i1,i3)*aeff0_2auxarray(t-i2-2,i2+2)
                   ttaa3=ttaa3 + x2ya_c(t,i1,i3)*aeff0_2auxarray(t-i2-3,i2+3)

                   ttab0=ttab0 + x2yb_c(t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)
                   ttab1=ttab1 + x2yb_c(t,i1,i3)*aeff0_2auxarray(t-i2-1,i2+1)
                   ttab2=ttab2 + x2yb_c(t,i1,i3)*aeff0_2auxarray(t-i2-2,i2+2)
                   ttab3=ttab3 + x2yb_c(t,i1,i3)*aeff0_2auxarray(t-i2-3,i2+3)

                   ttac0=ttac0 + x2yc_c(t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)
                   ttac1=ttac1 + x2yc_c(t,i1,i3)*aeff0_2auxarray(t-i2-1,i2+1)
                   ttac2=ttac2 + x2yc_c(t,i1,i3)*aeff0_2auxarray(t-i2-2,i2+2)
                   ttac3=ttac3 + x2yc_c(t,i1,i3)*aeff0_2auxarray(t-i2-3,i2+3)

                   ttae0=ttae0 + x2ye_c(t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)
                   ttae1=ttae1 + x2ye_c(t,i1,i3)*aeff0_2auxarray(t-i2-1,i2+1)
                   ttae2=ttae2 + x2ye_c(t,i1,i3)*aeff0_2auxarray(t-i2-2,i2+2)
                   ttae3=ttae3 + x2ye_c(t,i1,i3)*aeff0_2auxarray(t-i2-3,i2+3)

                   ! x^2y^2, b coefficients
                   ttba0=ttba0 + x2ya_c(t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)
                   ttba1=ttba1 + x2ya_c(t,i1,i3)*beff0_2auxarray(t-i2-1,i2+1)
                   ttba2=ttba2 + x2ya_c(t,i1,i3)*beff0_2auxarray(t-i2-2,i2+2)
                   ttba3=ttba3 + x2ya_c(t,i1,i3)*beff0_2auxarray(t-i2-3,i2+3)

                   ttbb0=ttbb0 + x2yb_c(t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)
                   ttbb1=ttbb1 + x2yb_c(t,i1,i3)*beff0_2auxarray(t-i2-1,i2+1)
                   ttbb2=ttbb2 + x2yb_c(t,i1,i3)*beff0_2auxarray(t-i2-2,i2+2)
                   ttbb3=ttbb3 + x2yb_c(t,i1,i3)*beff0_2auxarray(t-i2-3,i2+3)

                   ttbc0=ttbc0 + x2yc_c(t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)
                   ttbc1=ttbc1 + x2yc_c(t,i1,i3)*beff0_2auxarray(t-i2-1,i2+1)
                   ttbc2=ttbc2 + x2yc_c(t,i1,i3)*beff0_2auxarray(t-i2-2,i2+2)
                   ttbc3=ttbc3 + x2yc_c(t,i1,i3)*beff0_2auxarray(t-i2-3,i2+3)

                   ttbe0=ttbe0 + x2ye_c(t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)
                   ttbe1=ttbe1 + x2ye_c(t,i1,i3)*beff0_2auxarray(t-i2-1,i2+1)
                   ttbe2=ttbe2 + x2ye_c(t,i1,i3)*beff0_2auxarray(t-i2-2,i2+2)
                   ttbe3=ttbe3 + x2ye_c(t,i1,i3)*beff0_2auxarray(t-i2-3,i2+3)

                   ! x^2y^2, c coefficients
                   ttca0=ttca0 + x2ya_c(t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)
                   ttca1=ttca1 + x2ya_c(t,i1,i3)*ceff0_2auxarray(t-i2-1,i2+1)
                   ttca2=ttca2 + x2ya_c(t,i1,i3)*ceff0_2auxarray(t-i2-2,i2+2)
                   ttca3=ttca3 + x2ya_c(t,i1,i3)*ceff0_2auxarray(t-i2-3,i2+3)

                   ttcb0=ttcb0 + x2yb_c(t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)
                   ttcb1=ttcb1 + x2yb_c(t,i1,i3)*ceff0_2auxarray(t-i2-1,i2+1)
                   ttcb2=ttcb2 + x2yb_c(t,i1,i3)*ceff0_2auxarray(t-i2-2,i2+2)
                   ttcb3=ttcb3 + x2yb_c(t,i1,i3)*ceff0_2auxarray(t-i2-3,i2+3)

                   ttcc0=ttcc0 + x2yc_c(t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)
                   ttcc1=ttcc1 + x2yc_c(t,i1,i3)*ceff0_2auxarray(t-i2-1,i2+1)
                   ttcc2=ttcc2 + x2yc_c(t,i1,i3)*ceff0_2auxarray(t-i2-2,i2+2)
                   ttcc3=ttcc3 + x2yc_c(t,i1,i3)*ceff0_2auxarray(t-i2-3,i2+3)

                   ttce0=ttce0 + x2ye_c(t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)
                   ttce1=ttce1 + x2ye_c(t,i1,i3)*ceff0_2auxarray(t-i2-1,i2+1)
                   ttce2=ttce2 + x2ye_c(t,i1,i3)*ceff0_2auxarray(t-i2-2,i2+2)
                   ttce3=ttce3 + x2ye_c(t,i1,i3)*ceff0_2auxarray(t-i2-3,i2+3)

                   ! x^2y^2, e coefficients
                   ttea0=ttea0 + x2ya_c(t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)
                   ttea1=ttea1 + x2ya_c(t,i1,i3)*eeff0_2auxarray(t-i2-1,i2+1)
                   ttea2=ttea2 + x2ya_c(t,i1,i3)*eeff0_2auxarray(t-i2-2,i2+2)
                   ttea3=ttea3 + x2ya_c(t,i1,i3)*eeff0_2auxarray(t-i2-3,i2+3)

                   tteb0=tteb0 + x2yb_c(t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)
                   tteb1=tteb1 + x2yb_c(t,i1,i3)*eeff0_2auxarray(t-i2-1,i2+1)
                   tteb2=tteb2 + x2yb_c(t,i1,i3)*eeff0_2auxarray(t-i2-2,i2+2)
                   tteb3=tteb3 + x2yb_c(t,i1,i3)*eeff0_2auxarray(t-i2-3,i2+3)

                   ttec0=ttec0 + x2yc_c(t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)
                   ttec1=ttec1 + x2yc_c(t,i1,i3)*eeff0_2auxarray(t-i2-1,i2+1)
                   ttec2=ttec2 + x2yc_c(t,i1,i3)*eeff0_2auxarray(t-i2-2,i2+2)
                   ttec3=ttec3 + x2yc_c(t,i1,i3)*eeff0_2auxarray(t-i2-3,i2+3)

                   ttee0=ttee0 + x2ye_c(t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)
                   ttee1=ttee1 + x2ye_c(t,i1,i3)*eeff0_2auxarray(t-i2-1,i2+1)
                   ttee2=ttee2 + x2ye_c(t,i1,i3)*eeff0_2auxarray(t-i2-2,i2+2)
                   ttee3=ttee3 + x2ye_c(t,i1,i3)*eeff0_2auxarray(t-i2-3,i2+3)
                enddo
                y_c(i1,i2+0,i3)=y_c(i1,i2+0,i3)+dyi0
                y_c(i1,i2+1,i3)=y_c(i1,i2+1,i3)+dyi1
                y_c(i1,i2+2,i3)=y_c(i1,i2+2,i3)+dyi2
                y_c(i1,i2+3,i3)=y_c(i1,i2+3,i3)+dyi3

                y2za_c(i3,i1,i2+0)=tt0a0_2
                y2za_c(i3,i1,i2+1)=tt0a1_2
                y2za_c(i3,i1,i2+2)=tt0a2_2
                y2za_c(i3,i1,i2+3)=tt0a3_2
                          
                y2zb_c(i3,i1,i2+0)=tt0b0_2
                y2zb_c(i3,i1,i2+1)=tt0b1_2
                y2zb_c(i3,i1,i2+2)=tt0b2_2
                y2zb_c(i3,i1,i2+3)=tt0b3_2
                          
                y2zc_c(i3,i1,i2+0)=tt0c0_2
                y2zc_c(i3,i1,i2+1)=tt0c1_2
                y2zc_c(i3,i1,i2+2)=tt0c2_2
                y2zc_c(i3,i1,i2+3)=tt0c3_2
                          
                y2ze_c(i3,i1,i2+0)=tt0e0_2
                y2ze_c(i3,i1,i2+1)=tt0e1_2
                y2ze_c(i3,i1,i2+2)=tt0e2_2
                y2ze_c(i3,i1,i2+3)=tt0e3_2

                y4za_c(i3,i1,i2+0)=tt0a0_4
                y4za_c(i3,i1,i2+1)=tt0a1_4
                y4za_c(i3,i1,i2+2)=tt0a2_4
                y4za_c(i3,i1,i2+3)=tt0a3_4
                          
                y4zb_c(i3,i1,i2+0)=tt0b0_4
                y4zb_c(i3,i1,i2+1)=tt0b1_4
                y4zb_c(i3,i1,i2+2)=tt0b2_4
                y4zb_c(i3,i1,i2+3)=tt0b3_4
                          
                y4zc_c(i3,i1,i2+0)=tt0c0_4
                y4zc_c(i3,i1,i2+1)=tt0c1_4
                y4zc_c(i3,i1,i2+2)=tt0c2_4
                y4zc_c(i3,i1,i2+3)=tt0c3_4
                          
                y4ze_c(i3,i1,i2+0)=tt0e0_4
                y4ze_c(i3,i1,i2+1)=tt0e1_4
                y4ze_c(i3,i1,i2+2)=tt0e2_4
                y4ze_c(i3,i1,i2+3)=tt0e3_4

                ! x^2y^2, a coefficients
                x2y2aa_c(i3,i1,i2+0)=ttaa0
                x2y2aa_c(i3,i1,i2+1)=ttaa1
                x2y2aa_c(i3,i1,i2+2)=ttaa2
                x2y2aa_c(i3,i1,i2+3)=ttaa3

                x2y2ab_c(i3,i1,i2+0)=ttab0
                x2y2ab_c(i3,i1,i2+1)=ttab1
                x2y2ab_c(i3,i1,i2+2)=ttab2
                x2y2ab_c(i3,i1,i2+3)=ttab3

                x2y2ac_c(i3,i1,i2+0)=ttac0
                x2y2ac_c(i3,i1,i2+1)=ttac1
                x2y2ac_c(i3,i1,i2+2)=ttac2
                x2y2ac_c(i3,i1,i2+3)=ttac3

                x2y2ae_c(i3,i1,i2+0)=ttae0
                x2y2ae_c(i3,i1,i2+1)=ttae1
                x2y2ae_c(i3,i1,i2+2)=ttae2
                x2y2ae_c(i3,i1,i2+3)=ttae3

                ! x^2y^2, b coefficients
                x2y2ba_c(i3,i1,i2+0)=ttba0
                x2y2ba_c(i3,i1,i2+1)=ttba1
                x2y2ba_c(i3,i1,i2+2)=ttba2
                x2y2ba_c(i3,i1,i2+3)=ttba3

                x2y2bb_c(i3,i1,i2+0)=ttbb0
                x2y2bb_c(i3,i1,i2+1)=ttbb1
                x2y2bb_c(i3,i1,i2+2)=ttbb2
                x2y2bb_c(i3,i1,i2+3)=ttbb3

                x2y2bc_c(i3,i1,i2+0)=ttbc0
                x2y2bc_c(i3,i1,i2+1)=ttbc1
                x2y2bc_c(i3,i1,i2+2)=ttbc2
                x2y2bc_c(i3,i1,i2+3)=ttbc3

                x2y2be_c(i3,i1,i2+0)=ttbe0
                x2y2be_c(i3,i1,i2+1)=ttbe1
                x2y2be_c(i3,i1,i2+2)=ttbe2
                x2y2be_c(i3,i1,i2+3)=ttbe3

                ! x^2y^2, c coefficients
                x2y2ca_c(i3,i1,i2+0)=ttca0
                x2y2ca_c(i3,i1,i2+1)=ttca1
                x2y2ca_c(i3,i1,i2+2)=ttca2
                x2y2ca_c(i3,i1,i2+3)=ttca3

                x2y2cb_c(i3,i1,i2+0)=ttcb0
                x2y2cb_c(i3,i1,i2+1)=ttcb1
                x2y2cb_c(i3,i1,i2+2)=ttcb2
                x2y2cb_c(i3,i1,i2+3)=ttcb3

                x2y2cc_c(i3,i1,i2+0)=ttcc0
                x2y2cc_c(i3,i1,i2+1)=ttcc1
                x2y2cc_c(i3,i1,i2+2)=ttcc2
                x2y2cc_c(i3,i1,i2+3)=ttcc3

                x2y2ce_c(i3,i1,i2+0)=ttce0
                x2y2ce_c(i3,i1,i2+1)=ttce1
                x2y2ce_c(i3,i1,i2+2)=ttce2
                x2y2ce_c(i3,i1,i2+3)=ttce3

                ! x^2y^2, e coefficients
                x2y2ea_c(i3,i1,i2+0)=ttea0
                x2y2ea_c(i3,i1,i2+1)=ttea1
                x2y2ea_c(i3,i1,i2+2)=ttea2
                x2y2ea_c(i3,i1,i2+3)=ttea3

                x2y2eb_c(i3,i1,i2+0)=tteb0
                x2y2eb_c(i3,i1,i2+1)=tteb1
                x2y2eb_c(i3,i1,i2+2)=tteb2
                x2y2eb_c(i3,i1,i2+3)=tteb3

                x2y2ec_c(i3,i1,i2+0)=ttec0
                x2y2ec_c(i3,i1,i2+1)=ttec1
                x2y2ec_c(i3,i1,i2+2)=ttec2
                x2y2ec_c(i3,i1,i2+3)=ttec3

                x2y2ee_c(i3,i1,i2+0)=ttee0
                x2y2ee_c(i3,i1,i2+1)=ttee1
                x2y2ee_c(i3,i1,i2+2)=ttee2
                x2y2ee_c(i3,i1,i2+3)=ttee3
             enddo
             icur=i2
          else
             icur=ibxz_c(1,i1,i3)
          endif
  
          do i2=icur,ibxz_c(2,i1,i3)
             dyi=0.0_wp 

             tt0a0_2=0.d0 ; tt0a0_4=0.d0

             tt0b0_2=0.d0 ; tt0b0_4=0.d0

             tt0c0_2=0.d0 ; tt0c0_4=0.d0

             tt0e0_2=0.d0 ; tt0e0_4=0.d0

             ttaa0=0.d0 ; ttba0=0.d0 ; ttca0=0.d0 ; ttea0=0.d0

             ttab0=0.d0 ; ttbb0=0.d0 ; ttcb0=0.d0 ; tteb0=0.d0

             ttac0=0.d0 ; ttbc0=0.d0 ; ttcc0=0.d0 ; ttec0=0.d0

             ttae0=0.d0 ; ttbe0=0.d0 ; ttce0=0.d0 ; ttee0=0.d0
             do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
                dyi=dyi + xy_c(t,i1,i3)*aeff0array(t-i2-0,i2+0) + 3.d0*x4ya_c(t,i1,i3)*aeff0_2array(t-i2-0,i2+0) + &
                          3.d0*x2ya_c(t,i1,i3)*aeff0_4array(t-i2-0,i2+0)

                ! sss coefficients
                tt0a0_2=tt0a0_2 + xy_c(t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)

                tt0b0_2=tt0b0_2 + xy_c(t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)

                tt0c0_2=tt0c0_2 + xy_c(t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)

                tt0e0_2=tt0e0_2 + xy_c(t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)

                tt0a0_4=tt0a0_4 + xy_c(t,i1,i3)*aeff0_4auxarray(t-i2-0,i2+0)

                tt0b0_4=tt0b0_4 + xy_c(t,i1,i3)*beff0_4auxarray(t-i2-0,i2+0)

                tt0c0_4=tt0c0_4 + xy_c(t,i1,i3)*ceff0_4auxarray(t-i2-0,i2+0)

                tt0e0_4=tt0e0_4 + xy_c(t,i1,i3)*eeff0_4auxarray(t-i2-0,i2+0)

                ! x^2y^2, a coefficients
                ttaa0=ttaa0 + x2ya_c(t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)

                ttab0=ttab0 + x2yb_c(t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)

                ttac0=ttac0 + x2yc_c(t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)

                ttae0=ttae0 + x2ye_c(t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)

                ! x^2y^2, b coefficients
                ttba0=ttba0 + x2ya_c(t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)

                ttbb0=ttbb0 + x2yb_c(t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)

                ttbc0=ttbc0 + x2yc_c(t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)

                ttbe0=ttbe0 + x2ye_c(t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)

                ! x^2y^2, c coefficients
                ttca0=ttca0 + x2ya_c(t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)

                ttcb0=ttcb0 + x2yb_c(t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)

                ttcc0=ttcc0 + x2yc_c(t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)

                ttce0=ttce0 + x2ye_c(t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)

                ! x^2y^2, e coefficients
                ttea0=ttea0 + x2ya_c(t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)

                tteb0=tteb0 + x2yb_c(t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)

                ttec0=ttec0 + x2yc_c(t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)

                ttee0=ttee0 + x2ye_c(t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)


             enddo
             y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi


             y2za_c(i3,i1,i2+0)=tt0a0_2
                       
             y2zb_c(i3,i1,i2+0)=tt0b0_2
                       
             y2zc_c(i3,i1,i2+0)=tt0c0_2
                       
             y2ze_c(i3,i1,i2+0)=tt0e0_2

             y4za_c(i3,i1,i2+0)=tt0a0_4
                       
             y4zb_c(i3,i1,i2+0)=tt0b0_4
                       
             y4zc_c(i3,i1,i2+0)=tt0c0_4
                       
             y4ze_c(i3,i1,i2+0)=tt0e0_4

             ! x^2y^2, a coefficients
             x2y2aa_c(i3,i1,i2+0)=ttaa0

             x2y2ab_c(i3,i1,i2+0)=ttab0

             x2y2ac_c(i3,i1,i2+0)=ttac0

             x2y2ae_c(i3,i1,i2+0)=ttae0

             ! x^2y^2, b coefficients
             x2y2ba_c(i3,i1,i2+0)=ttba0

             x2y2bb_c(i3,i1,i2+0)=ttbb0

             x2y2bc_c(i3,i1,i2+0)=ttbc0

             x2y2be_c(i3,i1,i2+0)=ttbe0

             ! x^2y^2, c coefficients
             x2y2ca_c(i3,i1,i2+0)=ttca0

             x2y2cb_c(i3,i1,i2+0)=ttcb0

             x2y2cc_c(i3,i1,i2+0)=ttcc0

             x2y2ce_c(i3,i1,i2+0)=ttce0

             ! x^2y^2, e coefficients
             x2y2ea_c(i3,i1,i2+0)=ttea0

             x2y2eb_c(i3,i1,i2+0)=tteb0

             x2y2ec_c(i3,i1,i2+0)=ttec0

             x2y2ee_c(i3,i1,i2+0)=ttee0
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
                   dyi0=dyi0 + xy_f2(t,i1,i3)*beff0array(t-i2-0,i2+0) + &
                               3.d0*x4yb_f(1,t,i1,i3)*aeff0_2array(t-i2-0,i2+0) + &
                               3.d0*(x4ya_f(1,t,i1,i3)+x4yb_f(2,t,i1,i3))*beff0_2array(t-i2-0,i2+0) + &
                               3.d0*x2yb_f(1,t,i1,i3)*aeff0_4array(t-i2-0,i2+0) + &
                               3.d0*(x2ya_f(1,t,i1,i3)+x2yb_f(2,t,i1,i3))*beff0_4array(t-i2-0,i2+0)
                   dyi1=dyi1 + xy_f2(t,i1,i3)*beff0array(t-i2-1,i2+1) + &
                               3.d0*x4yb_f(1,t,i1,i3)*aeff0_2array(t-i2-1,i2+1) + &
                               3.d0*(x4ya_f(1,t,i1,i3)+x4yb_f(2,t,i1,i3))*beff0_2array(t-i2-1,i2+1) + &
                               3.d0*x2yb_f(1,t,i1,i3)*aeff0_4array(t-i2-1,i2+1) + &
                               3.d0*(x2ya_f(1,t,i1,i3)+x2yb_f(2,t,i1,i3))*beff0_4array(t-i2-1,i2+1)
                   dyi2=dyi2 + xy_f2(t,i1,i3)*beff0array(t-i2-2,i2+2) + &
                               3.d0*x4yb_f(1,t,i1,i3)*aeff0_2array(t-i2-2,i2+2) + &
                               3.d0*(x4ya_f(1,t,i1,i3)+x4yb_f(2,t,i1,i3))*beff0_2array(t-i2-2,i2+2) + &
                               3.d0*x2yb_f(1,t,i1,i3)*aeff0_4array(t-i2-2,i2+2) + &
                               3.d0*(x2ya_f(1,t,i1,i3)+x2yb_f(2,t,i1,i3))*beff0_4array(t-i2-2,i2+2)
                   dyi3=dyi3 + xy_f2(t,i1,i3)*beff0array(t-i2-3,i2+3) + &
                               3.d0*x4yb_f(1,t,i1,i3)*aeff0_2array(t-i2-3,i2+3) + &
                               3.d0*(x4ya_f(1,t,i1,i3)+x4yb_f(2,t,i1,i3))*beff0_2array(t-i2-3,i2+3) + &
                               3.d0*x2yb_f(1,t,i1,i3)*aeff0_4array(t-i2-3,i2+3) + &
                               3.d0*(x2ya_f(1,t,i1,i3)+x2yb_f(2,t,i1,i3))*beff0_4array(t-i2-3,i2+3)
                enddo
                y_c(i1,i2+0,i3)=y_c(i1,i2+0,i3)+dyi0
                y_c(i1,i2+1,i3)=y_c(i1,i2+1,i3)+dyi1
                y_c(i1,i2+2,i3)=y_c(i1,i2+2,i3)+dyi2
                y_c(i1,i2+3,i3)=y_c(i1,i2+3,i3)+dyi3
             enddo
             istart=i2
          endif
  
          do i2=istart,iend
             dyi0=0.0_wp
             do t=max(ibxz_f(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_f(2,i1,i3))
                dyi0=dyi0 + xy_f2(t,i1,i3)*beff0array(t-i2-0,i2+0) + &
                            3.d0*x4yb_f(1,t,i1,i3)*aeff0_2array(t-i2-0,i2+0) + &
                            3.d0*(x4ya_f(1,t,i1,i3)+x4yb_f(2,t,i1,i3))*beff0_2array(t-i2-0,i2+0) + &
                            3.d0*x2yb_f(1,t,i1,i3)*aeff0_4array(t-i2-0,i2+0) + &
                            3.d0*(x2ya_f(1,t,i1,i3)+x2yb_f(2,t,i1,i3))*beff0_4array(t-i2-0,i2+0)
             enddo
             y_c(i1,i2+0,i3)=y_c(i1,i2+0,i3)+dyi0
          enddo
  
           if (ibxz_f(2,i1,i3)-ibxz_f(1,i1,i3).ge.4) then
             do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)-4,4
                dyi0=0.0_wp 
                dyi1=0.0_wp 
                dyi2=0.0_wp 
                dyi3=0.0_wp 
                tt10=0.0_wp 
                tt11=0.0_wp 
                tt12=0.0_wp 
                tt13=0.0_wp 
                tt20=0.0_wp 
                tt21=0.0_wp 
                tt22=0.0_wp 
                tt23=0.0_wp 
                tt30=0.0_wp 
                tt31=0.0_wp 
                tt32=0.0_wp 
                tt33=0.0_wp 
                do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2+3,ibxz_c(2,i1,i3))
                   dyi0=dyi0 + xy_c(t,i1,i3)*ceff0array(t-i2-0,i2+0)
                   dyi1=dyi1 + xy_c(t,i1,i3)*ceff0array(t-i2-1,i2+1)
                   dyi2=dyi2 + xy_c(t,i1,i3)*ceff0array(t-i2-2,i2+2)
                   dyi3=dyi3 + xy_c(t,i1,i3)*ceff0array(t-i2-3,i2+3)
  
                   tt10=tt10 + 3.d0*( x4yc_c(t,i1,i3)*aeff0_2array(t-i2-0,i2+0) + x2yc_c(t,i1,i3)*aeff0_4array(t-i2-0,i2+0) )
                   tt11=tt11 + 3.d0*( x4yc_c(t,i1,i3)*aeff0_2array(t-i2-1,i2+1) + x2yc_c(t,i1,i3)*aeff0_4array(t-i2-1,i2+1) )
                   tt12=tt12 + 3.d0*( x4yc_c(t,i1,i3)*aeff0_2array(t-i2-2,i2+2) + x2yc_c(t,i1,i3)*aeff0_4array(t-i2-2,i2+2) )
                   tt13=tt13 + 3.d0*( x4yc_c(t,i1,i3)*aeff0_2array(t-i2-3,i2+3) + x2yc_c(t,i1,i3)*aeff0_4array(t-i2-3,i2+3) )

                   tt20=tt20 + 3.d0*( x4ya_c(t,i1,i3)*ceff0_2array(t-i2-0,i2+0) + x2ya_c(t,i1,i3)*ceff0_4array(t-i2-0,i2+0) )
                   tt21=tt21 + 3.d0*( x4ya_c(t,i1,i3)*ceff0_2array(t-i2-1,i2+1) + x2ya_c(t,i1,i3)*ceff0_4array(t-i2-1,i2+1) )
                   tt22=tt22 + 3.d0*( x4ya_c(t,i1,i3)*ceff0_2array(t-i2-2,i2+2) + x2ya_c(t,i1,i3)*ceff0_4array(t-i2-2,i2+2) )
                   tt23=tt23 + 3.d0*( x4ya_c(t,i1,i3)*ceff0_2array(t-i2-3,i2+3) + x2ya_c(t,i1,i3)*ceff0_4array(t-i2-3,i2+3) )

                   tt30=tt30 + 3.d0*( x4yc_c(t,i1,i3)*ceff0_2array(t-i2-0,i2+0) + x2yc_c(t,i1,i3)*ceff0_4array(t-i2-0,i2+0) )
                   tt31=tt31 + 3.d0*( x4yc_c(t,i1,i3)*ceff0_2array(t-i2-1,i2+1) + x2yc_c(t,i1,i3)*ceff0_4array(t-i2-1,i2+1) )
                   tt32=tt32 + 3.d0*( x4yc_c(t,i1,i3)*ceff0_2array(t-i2-2,i2+2) + x2yc_c(t,i1,i3)*ceff0_4array(t-i2-2,i2+2) )
                   tt33=tt33 + 3.d0*( x4yc_c(t,i1,i3)*ceff0_2array(t-i2-3,i2+3) + x2yc_c(t,i1,i3)*ceff0_4array(t-i2-3,i2+3) )
                enddo
                y_f(2,i1,i2+0,i3)=y_f(2,i1,i2+0,i3)+dyi0
                y_f(2,i1,i2+1,i3)=y_f(2,i1,i2+1,i3)+dyi1
                y_f(2,i1,i2+2,i3)=y_f(2,i1,i2+2,i3)+dyi2
                y_f(2,i1,i2+3,i3)=y_f(2,i1,i2+3,i3)+dyi3
  
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
             enddo
             icur=i2
          else
             icur=ibxz_f(1,i1,i3)
          endif
  
          do i2=icur,ibxz_f(2,i1,i3)
             dyi0=0.0_wp 
             tt10=0.0_wp 
             tt20=0.0_wp 
             tt30=0.0_wp 
             do t=max(ibxz_c(1,i1,i3),lowfil+i2),min(lupfil+i2,ibxz_c(2,i1,i3))
                dyi0=dyi0 + xy_c(t,i1,i3)*ceff0array(t-i2-0,i2)

                tt10=tt10 + 3.d0*( x4yc_c(t,i1,i3)*aeff0_2array(t-i2-0,i2+0) + x2yc_c(t,i1,i3)*aeff0_4array(t-i2-0,i2+0) )

                tt20=tt20 + 3.d0*( x4ya_c(t,i1,i3)*ceff0_2array(t-i2-0,i2+0) + x2ya_c(t,i1,i3)*ceff0_4array(t-i2-0,i2+0) )

                tt30=tt30 + 3.d0*( x4yc_c(t,i1,i3)*ceff0_2array(t-i2-0,i2+0) + x2yc_c(t,i1,i3)*ceff0_4array(t-i2-0,i2+0) )
             enddo
             y_f(2,i1,i2+0,i3)=y_f(2,i1,i2+0,i3)+dyi0
             y_f(1,i1,i2+0,i3)=y_f(1,i1,i2+0,i3)+tt10
             y_f(2,i1,i2+0,i3)=y_f(2,i1,i2+0,i3)+tt20
             y_f(3,i1,i2+0,i3)=y_f(3,i1,i2+0,i3)+tt30
          enddo
       enddo
    enddo
    !!!$omp enddo
  
  
    ! wavelet part
  
    !!!$omp do
    do i3=nfl3,nfu3
       do i1=nfl1,nfu1
          do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
             ! Get the effective filters for the y dimension
             tt10 = 0.d0
             tt20 = 0.d0
             tt30 = 0.d0
             tt40 = 0.d0
             tt50 = 0.d0
             tt60 = 0.d0
             tt70 = 0.d0

             tt1a0=0.d0
             tt1b0=0.d0
             tt1c0=0.d0
             tt1e0=0.d0

             tt2a0=0.d0
             tt2b0=0.d0
             tt2c0=0.d0
             tt2e0=0.d0

             tt3a0=0.d0
             tt3b0=0.d0
             tt3c0=0.d0
             tt3e0=0.d0

             tt4a0=0.d0
             tt4b0=0.d0
             tt4c0=0.d0
             tt4e0=0.d0

             tt5a0=0.d0
             tt5b0=0.d0
             tt5c0=0.d0
             tt5e0=0.d0

             tt6a0=0.d0
             tt6b0=0.d0
             tt6c0=0.d0
             tt6e0=0.d0

             tt7a0=0.d0
             tt7b0=0.d0
             tt7c0=0.d0
             tt7e0=0.d0

                ttaa1=0.d0 ; ttab1=0.d0 ; ttac1=0.d0 ; ttae1=0.d0
                ttaa2=0.d0 ; ttab2=0.d0 ; ttac2=0.d0 ; ttae2=0.d0
                ttaa3=0.d0 ; ttab3=0.d0 ; ttac3=0.d0 ; ttae3=0.d0
                ttaa4=0.d0 ; ttab4=0.d0 ; ttac4=0.d0 ; ttae4=0.d0
                ttaa5=0.d0 ; ttab5=0.d0 ; ttac5=0.d0 ; ttae5=0.d0
                ttaa6=0.d0 ; ttab6=0.d0 ; ttac6=0.d0 ; ttae6=0.d0
                ttaa7=0.d0 ; ttab7=0.d0 ; ttac7=0.d0 ; ttae7=0.d0

                ttba1=0.d0 ; ttbb1=0.d0 ; ttbc1=0.d0 ; ttbe1=0.d0
                ttba2=0.d0 ; ttbb2=0.d0 ; ttbc2=0.d0 ; ttbe2=0.d0
                ttba3=0.d0 ; ttbb3=0.d0 ; ttbc3=0.d0 ; ttbe3=0.d0
                ttba4=0.d0 ; ttbb4=0.d0 ; ttbc4=0.d0 ; ttbe4=0.d0
                ttba5=0.d0 ; ttbb5=0.d0 ; ttbc5=0.d0 ; ttbe5=0.d0
                ttba6=0.d0 ; ttbb6=0.d0 ; ttbc6=0.d0 ; ttbe6=0.d0
                ttba7=0.d0 ; ttbb7=0.d0 ; ttbc7=0.d0 ; ttbe7=0.d0

                ttca1=0.d0 ; ttcb1=0.d0 ; ttcc1=0.d0 ; ttce1=0.d0
                ttca2=0.d0 ; ttcb2=0.d0 ; ttcc2=0.d0 ; ttce2=0.d0
                ttca3=0.d0 ; ttcb3=0.d0 ; ttcc3=0.d0 ; ttce3=0.d0
                ttca4=0.d0 ; ttcb4=0.d0 ; ttcc4=0.d0 ; ttce4=0.d0
                ttca5=0.d0 ; ttcb5=0.d0 ; ttcc5=0.d0 ; ttce5=0.d0
                ttca6=0.d0 ; ttcb6=0.d0 ; ttcc6=0.d0 ; ttce6=0.d0
                ttca7=0.d0 ; ttcb7=0.d0 ; ttcc7=0.d0 ; ttce7=0.d0

                ttea1=0.d0 ; tteb1=0.d0 ; ttec1=0.d0 ; ttee1=0.d0
                ttea2=0.d0 ; tteb2=0.d0 ; ttec2=0.d0 ; ttee2=0.d0
                ttea3=0.d0 ; tteb3=0.d0 ; ttec3=0.d0 ; ttee3=0.d0
                ttea4=0.d0 ; tteb4=0.d0 ; ttec4=0.d0 ; ttee4=0.d0
                ttea5=0.d0 ; tteb5=0.d0 ; ttec5=0.d0 ; ttee5=0.d0
                ttea6=0.d0 ; tteb6=0.d0 ; ttec6=0.d0 ; ttee6=0.d0
                ttea7=0.d0 ; tteb7=0.d0 ; ttec7=0.d0 ; ttee7=0.d0
                
                
             do l=max(nfl2-i2,lowfil),min(lupfil,nfu2-i2)
                tt10 = tt10 + xy_f(1,i2+l,i1,i3)*aeff0array(l,i2) + xy_f(3,i2+l,i1,i3)*beff0array(l,i2) + &
                              3.d0* x4ye_f(1,i2+l,i1,i3)*                      aeff0_2array(l,i2) + &
                              3.d0*(x4yc_f(1,i2+l,i1,i3)+x4ye_f(2,i2+l,i1,i3))*beff0_2array(l,i2) + &
                              3.d0* x2ye_f(1,i2+l,i1,i3)*                      aeff0_4array(l,i2) + &
                              3.d0*(x2yc_f(1,i2+l,i1,i3)+x2ye_f(2,i2+l,i1,i3))*beff0_4array(l,i2)
  
                tt20 = tt20 + xy_f(2,i2+l,i1,i3)*eeff0array(l,i2) +                                      &
                              3.d0* x4yb_f(1,i2+l,i1,i3)*                      ceff0_2array(l,i2) + &
                              3.d0*(x4ya_f(1,i2+l,i1,i3)+x4yb_f(2,i2+l,i1,i3))*eeff0_2array(l,i2) + &
                              3.d0* x2yb_f(1,i2+l,i1,i3)*                      ceff0_4array(l,i2) + &
                              3.d0*(x2ya_f(1,i2+l,i1,i3)+x2yb_f(2,i2+l,i1,i3))*eeff0_4array(l,i2)
  
                tt30 = tt30 + xy_f(1,i2+l,i1,i3)*ceff0array(l,i2) + xy_f(3,i2+l,i1,i3)*eeff0array(l,i2) + &
                              3.d0* x4ye_f(1,i2+l,i1,i3)*                      ceff0_2array(l,i2) + &
                              3.d0*(x4yc_f(1,i2+l,i1,i3)+x4ye_f(2,i2+l,i1,i3))*eeff0_2array(l,i2) + &
                              3.d0* x2ye_f(1,i2+l,i1,i3)*                      ceff0_4array(l,i2) + &
                              3.d0*(x2yc_f(1,i2+l,i1,i3)+x2ye_f(2,i2+l,i1,i3))*eeff0_4array(l,i2)
  
                tt40 = tt40 + xy_f(4,i2+l,i1,i3)*aeff0array(l,i2) + xy_f(6,i2+l,i1,i3)*beff0array(l,i2) + &
                              3.d0*(x4ya_f(2,i2+l,i1,i3)+x4yb_f(3,i2+l,i1,i3))*aeff0_2array(l,i2) + &
                              3.d0*(x4ya_f(3,i2+l,i1,i3)+x4yb_f(4,i2+l,i1,i3))*beff0_2array(l,i2) + &
                              3.d0*(x2ya_f(2,i2+l,i1,i3)+x2yb_f(3,i2+l,i1,i3))*aeff0_4array(l,i2) + &
                              3.d0*(x2ya_f(3,i2+l,i1,i3)+x2yb_f(4,i2+l,i1,i3))*beff0_4array(l,i2)
  
                tt50 = tt50 + xy_f(5,i2+l,i1,i3)*aeff0array(l,i2) + xy_f(7,i2+l,i1,i3)*beff0array(l,i2) + &
                              3.d0*(x4yc_f(2,i2+l,i1,i3)+x4ye_f(3,i2+l,i1,i3))*aeff0_2array(l,i2) + &
                              3.d0*(x4yc_f(3,i2+l,i1,i3)+x4ye_f(4,i2+l,i1,i3))*beff0_2array(l,i2) + &
                              3.d0*(x2yc_f(2,i2+l,i1,i3)+x2ye_f(3,i2+l,i1,i3))*aeff0_4array(l,i2) + &
                              3.d0*(x2yc_f(3,i2+l,i1,i3)+x2ye_f(4,i2+l,i1,i3))*beff0_4array(l,i2)
                
                tt60 = tt60 + xy_f(4,i2+l,i1,i3)*ceff0array(l,i2) + xy_f(6,i2+l,i1,i3)*eeff0array(l,i2) + &
                              3.d0*(x4ya_f(2,i2+l,i1,i3)+x4yb_f(3,i2+l,i1,i3))*ceff0_2array(l,i2) + &
                              3.d0*(x4ya_f(3,i2+l,i1,i3)+x4yb_f(4,i2+l,i1,i3))*eeff0_2array(l,i2) + &
                              3.d0*(x2ya_f(2,i2+l,i1,i3)+x2yb_f(3,i2+l,i1,i3))*ceff0_4array(l,i2) + &
                              3.d0*(x2ya_f(3,i2+l,i1,i3)+x2yb_f(4,i2+l,i1,i3))*eeff0_4array(l,i2)
  
                tt70 = tt70 + xy_f(5,i2+l,i1,i3)*ceff0array(l,i2) + xy_f(7,i2+l,i1,i3)*eeff0array(l,i2) + &
                              3.d0*(x4yc_f(2,i2+l,i1,i3)+x4ye_f(3,i2+l,i1,i3))*ceff0_2array(l,i2) + &
                              3.d0*(x4yc_f(3,i2+l,i1,i3)+x4ye_f(4,i2+l,i1,i3))*eeff0_2array(l,i2) + &
                              3.d0*(x2yc_f(2,i2+l,i1,i3)+x2ye_f(3,i2+l,i1,i3))*ceff0_4array(l,i2) + &
                              3.d0*(x2yc_f(3,i2+l,i1,i3)+x2ye_f(4,i2+l,i1,i3))*eeff0_4array(l,i2)

                ! dss coefficients
                tt1a0_2=tt1a0_2 + xy_f(1,i2+l,i1,i3)*aeff0_2auxarray(l,i2)
                tt1c0_2=tt1c0_2 + xy_f(1,i2+l,i1,i3)*ceff0_2auxarray(l,i2)
                tt1a0_4=tt1a0_4 + xy_f(1,i2+l,i1,i3)*aeff0_4auxarray(l,i2)
                tt1c0_4=tt1c0_4 + xy_f(1,i2+l,i1,i3)*ceff0_4auxarray(l,i2)
                ! sds coefficients
                tt2b0_2=tt2b0_2 + xy_f(2,i2+l,i1,i3)*beff0_2auxarray(l,i2)
                tt2e0_2=tt2e0_2 + xy_f(2,i2+l,i1,i3)*eeff0_2auxarray(l,i2)
                tt2b0_4=tt2b0_4 + xy_f(2,i2+l,i1,i3)*beff0_4auxarray(l,i2)
                tt2e0_4=tt2e0_4 + xy_f(2,i2+l,i1,i3)*eeff0_4auxarray(l,i2)
                ! dds coefficients
                tt3b0_2=tt3b0_2 + xy_f(3,i2+l,i1,i3)*beff0_2auxarray(l,i2)
                tt3e0_2=tt3e0_2 + xy_f(3,i2+l,i1,i3)*eeff0_2auxarray(l,i2)
                tt3b0_4=tt3b0_4 + xy_f(3,i2+l,i1,i3)*beff0_4auxarray(l,i2)
                tt3e0_4=tt3e0_4 + xy_f(3,i2+l,i1,i3)*eeff0_4auxarray(l,i2)
                ! ssd coefficients
                tt4a0_2=tt4a0_2 + xy_f(4,i2+l,i1,i3)*aeff0_2auxarray(l,i2)
                tt4c0_2=tt4c0_2 + xy_f(4,i2+l,i1,i3)*ceff0_2auxarray(l,i2)
                tt4a0_4=tt4a0_4 + xy_f(4,i2+l,i1,i3)*aeff0_4auxarray(l,i2)
                tt4c0_4=tt4c0_4 + xy_f(4,i2+l,i1,i3)*ceff0_4auxarray(l,i2)
                ! dsd coefficients
                tt5a0_2=tt5a0_2 + xy_f(5,i2+l,i1,i3)*aeff0_2auxarray(l,i2)
                tt5c0_2=tt5c0_2 + xy_f(5,i2+l,i1,i3)*ceff0_2auxarray(l,i2)
                tt5a0_4=tt5a0_4 + xy_f(5,i2+l,i1,i3)*aeff0_4auxarray(l,i2)
                tt5c0_4=tt5c0_4 + xy_f(5,i2+l,i1,i3)*ceff0_4auxarray(l,i2)
                ! sdd coefficients
                tt6b0_2=tt6b0_2 + xy_f(6,i2+l,i1,i3)*beff0_2auxarray(l,i2)
                tt6e0_2=tt6e0_2 + xy_f(6,i2+l,i1,i3)*eeff0_2auxarray(l,i2)
                tt6b0_4=tt6b0_4 + xy_f(6,i2+l,i1,i3)*beff0_4auxarray(l,i2)
                tt6e0_4=tt6e0_4 + xy_f(6,i2+l,i1,i3)*eeff0_4auxarray(l,i2)
                ! ddd coefficients
                tt7b0_2=tt7b0_2 + xy_f(7,i2+l,i1,i3)*beff0_2auxarray(l,i2)
                tt7b0_2=tt7b0_2 + xy_f(7,i2+l,i1,i3)*beff0_2auxarray(l,i2)
                tt7e0_4=tt7e0_4 + xy_f(7,i2+l,i1,i3)*eeff0_4auxarray(l,i2)
                tt7e0_4=tt7e0_4 + xy_f(7,i2+l,i1,i3)*eeff0_4auxarray(l,i2)


                ! x^2y^2, a coefficients
                ttaa1=ttaa1 + x2ya_f2(1,t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)
                ttaa2=ttaa2 + x2ya_f2(2,t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)
                ttaa3=ttaa3 + x2ya_f2(3,t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)
                ttaa4=ttaa4 + x2ya_f2(4,t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)
                ttaa5=ttaa5 + x2ya_f2(5,t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)
                ttaa6=ttaa6 + x2ya_f2(6,t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)
                ttaa7=ttaa7 + x2ya_f2(7,t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)

                ttab1=ttab1 + x2yb_f2(1,t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)
                ttab2=ttab2 + x2yb_f2(2,t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)
                ttab3=ttab3 + x2yb_f2(3,t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)
                ttab4=ttab4 + x2yb_f2(4,t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)
                ttab5=ttab5 + x2yb_f2(5,t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)
                ttab6=ttab6 + x2yb_f2(6,t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)
                ttab7=ttab7 + x2yb_f2(7,t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)

                ttac1=ttac1 + x2yc_f2(1,t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)
                ttac2=ttac2 + x2yc_f2(2,t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)
                ttac3=ttac3 + x2yc_f2(3,t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)
                ttac4=ttac4 + x2yc_f2(4,t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)
                ttac5=ttac5 + x2yc_f2(5,t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)
                ttac6=ttac6 + x2yc_f2(6,t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)
                ttac7=ttac7 + x2yc_f2(7,t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)

                ttae1=ttae1 + x2ye_f2(1,t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)
                ttae2=ttae2 + x2ye_f2(2,t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)
                ttae3=ttae3 + x2ye_f2(3,t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)
                ttae4=ttae4 + x2ye_f2(4,t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)
                ttae5=ttae5 + x2ye_f2(5,t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)
                ttae6=ttae6 + x2ye_f2(6,t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)
                ttae7=ttae7 + x2ye_f2(7,t,i1,i3)*aeff0_2auxarray(t-i2-0,i2+0)

                ! x^2y^2, b coefficients
                ttba1=ttba1 + x2ya_f2(1,t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)
                ttba2=ttba2 + x2ya_f2(2,t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)
                ttba3=ttba3 + x2ya_f2(3,t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)
                ttba4=ttba4 + x2ya_f2(4,t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)
                ttba5=ttba5 + x2ya_f2(5,t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)
                ttba6=ttba6 + x2ya_f2(6,t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)
                ttba7=ttba7 + x2ya_f2(7,t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)

                ttbb1=ttbb1 + x2yb_f2(1,t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)
                ttbb2=ttbb2 + x2yb_f2(2,t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)
                ttbb3=ttbb3 + x2yb_f2(3,t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)
                ttbb4=ttbb4 + x2yb_f2(4,t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)
                ttbb5=ttbb5 + x2yb_f2(5,t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)
                ttbb6=ttbb6 + x2yb_f2(6,t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)
                ttbb7=ttbb7 + x2yb_f2(7,t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)

                ttbc1=ttbc1 + x2yc_f2(1,t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)
                ttbc2=ttbc2 + x2yc_f2(2,t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)
                ttbc3=ttbc3 + x2yc_f2(3,t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)
                ttbc4=ttbc4 + x2yc_f2(4,t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)
                ttbc5=ttbc5 + x2yc_f2(5,t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)
                ttbc6=ttbc6 + x2yc_f2(6,t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)
                ttbc7=ttbc7 + x2yc_f2(7,t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)

                ttbe1=ttbe1 + x2ye_f2(1,t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)
                ttbe2=ttbe2 + x2ye_f2(2,t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)
                ttbe3=ttbe3 + x2ye_f2(3,t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)
                ttbe4=ttbe4 + x2ye_f2(4,t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)
                ttbe5=ttbe5 + x2ye_f2(5,t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)
                ttbe6=ttbe6 + x2ye_f2(6,t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)
                ttbe7=ttbe7 + x2ye_f2(7,t,i1,i3)*beff0_2auxarray(t-i2-0,i2+0)

                ! x^2y^2, c coefficients
                ttca1=ttca1 + x2ya_f2(1,t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)
                ttca2=ttca2 + x2ya_f2(2,t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)
                ttca3=ttca3 + x2ya_f2(3,t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)
                ttca4=ttca4 + x2ya_f2(4,t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)
                ttca5=ttca5 + x2ya_f2(5,t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)
                ttca6=ttca6 + x2ya_f2(6,t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)
                ttca7=ttca7 + x2ya_f2(7,t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)

                ttcb1=ttcb1 + x2yb_f2(1,t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)
                ttcb2=ttcb2 + x2yb_f2(2,t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)
                ttcb3=ttcb3 + x2yb_f2(3,t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)
                ttcb4=ttcb4 + x2yb_f2(4,t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)
                ttcb5=ttcb5 + x2yb_f2(5,t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)
                ttcb6=ttcb6 + x2yb_f2(6,t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)
                ttcb7=ttcb7 + x2yb_f2(7,t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)

                ttcc1=ttcc1 + x2yc_f2(1,t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)
                ttcc2=ttcc2 + x2yc_f2(2,t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)
                ttcc3=ttcc3 + x2yc_f2(3,t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)
                ttcc4=ttcc4 + x2yc_f2(4,t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)
                ttcc5=ttcc5 + x2yc_f2(5,t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)
                ttcc6=ttcc6 + x2yc_f2(6,t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)
                ttcc7=ttcc7 + x2yc_f2(7,t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)

                ttce1=ttce1 + x2ye_f2(1,t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)
                ttce2=ttce2 + x2ye_f2(2,t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)
                ttce3=ttce3 + x2ye_f2(3,t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)
                ttce4=ttce4 + x2ye_f2(4,t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)
                ttce5=ttce5 + x2ye_f2(5,t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)
                ttce6=ttce6 + x2ye_f2(6,t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)
                ttce7=ttce7 + x2ye_f2(7,t,i1,i3)*ceff0_2auxarray(t-i2-0,i2+0)

                ! x^2y^2, e coefficients
                ttea1=ttea1 + x2ya_f2(1,t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)
                ttea2=ttea2 + x2ya_f2(2,t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)
                ttea3=ttea3 + x2ya_f2(3,t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)
                ttea4=ttea4 + x2ya_f2(4,t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)
                ttea5=ttea5 + x2ya_f2(5,t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)
                ttea6=ttea6 + x2ya_f2(6,t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)
                ttea7=ttea7 + x2ya_f2(7,t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)

                tteb1=tteb1 + x2yb_f2(1,t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)
                tteb2=tteb2 + x2yb_f2(2,t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)
                tteb3=tteb3 + x2yb_f2(3,t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)
                tteb4=tteb4 + x2yb_f2(4,t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)
                tteb5=tteb5 + x2yb_f2(5,t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)
                tteb6=tteb6 + x2yb_f2(6,t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)
                tteb7=tteb7 + x2yb_f2(7,t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)

                ttec1=ttec1 + x2yc_f2(1,t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)
                ttec2=ttec2 + x2yc_f2(2,t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)
                ttec3=ttec3 + x2yc_f2(3,t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)
                ttec4=ttec4 + x2yc_f2(4,t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)
                ttec5=ttec5 + x2yc_f2(5,t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)
                ttec6=ttec6 + x2yc_f2(6,t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)
                ttec7=ttec7 + x2yc_f2(7,t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)

                ttee1=ttee1 + x2ye_f2(1,t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)
                ttee2=ttee2 + x2ye_f2(2,t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)
                ttee3=ttee3 + x2ye_f2(3,t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)
                ttee4=ttee4 + x2ye_f2(4,t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)
                ttee5=ttee5 + x2ye_f2(5,t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)
                ttee6=ttee6 + x2ye_f2(6,t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)
                ttee7=ttee7 + x2ye_f2(7,t,i1,i3)*eeff0_2auxarray(t-i2-0,i2+0)
             enddo
             y_f(4,i1,i2,i3)=y_f(4,i1,i2,i3)+tt40
             y_f(2,i1,i2,i3)=y_f(2,i1,i2,i3)+tt20
             y_f(1,i1,i2,i3)=y_f(1,i1,i2,i3)+tt10
             y_f(6,i1,i2,i3)=y_f(6,i1,i2,i3)+tt60
             y_f(5,i1,i2,i3)=y_f(5,i1,i2,i3)+tt50
             y_f(3,i1,i2,i3)=y_f(3,i1,i2,i3)+tt30
             y_f(7,i1,i2,i3)=y_f(7,i1,i2,i3)+tt70

             ! dss coefficients
             y4za_f(1,i3,i1,i2)=tt1a0_4
             y4zc_f(1,i3,i1,i2)=tt1c0_4
             y2za_f(1,i3,i1,i2)=tt1a0_2
             y2zc_f(1,i3,i1,i2)=tt1c0_2
             ! sds coefficients
             y4zb_f(1,i3,i1,i2)=tt2b0_4
             y4ze_f(1,i3,i1,i2)=tt2e0_4
             y2zb_f(1,i3,i1,i2)=tt2b0_2
             y2ze_f(1,i3,i1,i2)=tt2e0_2
             ! dds coefficients
             y4zb_f(2,i3,i1,i2)=tt3b0_4
             y4ze_f(2,i3,i1,i2)=tt3e0_4
             y2zb_f(2,i3,i1,i2)=tt3b0_2
             y2ze_f(2,i3,i1,i2)=tt3e0_2
             ! ssd coefficients
             y4za_f(2,i3,i1,i2)=tt4a0_4
             y4zc_f(2,i3,i1,i2)=tt4c0_4
             y2za_f(2,i3,i1,i2)=tt4a0_2
             y2zc_f(2,i3,i1,i2)=tt4c0_2
             ! dsd coefficients
             y4za_f(3,i3,i1,i2)=tt5a0_4
             y4zc_f(3,i3,i1,i2)=tt5c0_4
             y2za_f(3,i3,i1,i2)=tt5a0_2
             y2zc_f(3,i3,i1,i2)=tt5c0_2
             ! sdd coefficients
             y4zb_f(3,i3,i1,i2)=tt6b0_4
             y4ze_f(3,i3,i1,i2)=tt6e0_4
             y2zb_f(3,i3,i1,i2)=tt6b0_2
             y2ze_f(3,i3,i1,i2)=tt6e0_2
             ! sdd coefficients
             y4zb_f(4,i3,i1,i2)=tt7b0_4
             y4ze_f(4,i3,i1,i2)=tt7e0_4
             y2zb_f(4,i3,i1,i2)=tt7b0_2
             y2ze_f(4,i3,i1,i2)=tt7e0_2




                ! x^2y^2, a coefficients
                x2y2aa_f(1,i3,i1,i2)=ttaa1
                x2y2aa_f(2,i3,i1,i2)=ttaa2
                x2y2aa_f(3,i3,i1,i2)=ttaa3
                x2y2aa_f(4,i3,i1,i2)=ttaa4
                x2y2aa_f(5,i3,i1,i2)=ttaa5
                x2y2aa_f(6,i3,i1,i2)=ttaa6
                x2y2aa_f(7,i3,i1,i2)=ttaa7

                x2y2ab_f(1,i3,i1,i2)=ttab1
                x2y2ab_f(2,i3,i1,i2)=ttab2
                x2y2ab_f(3,i3,i1,i2)=ttab3
                x2y2ab_f(4,i3,i1,i2)=ttab4
                x2y2ab_f(5,i3,i1,i2)=ttab5
                x2y2ab_f(6,i3,i1,i2)=ttab6
                x2y2ab_f(7,i3,i1,i2)=ttab7

                x2y2ac_f(1,i3,i1,i2)=ttac1
                x2y2ac_f(2,i3,i1,i2)=ttac2
                x2y2ac_f(3,i3,i1,i2)=ttac3
                x2y2ac_f(4,i3,i1,i2)=ttac4
                x2y2ac_f(5,i3,i1,i2)=ttac5
                x2y2ac_f(6,i3,i1,i2)=ttac6
                x2y2ac_f(7,i3,i1,i2)=ttac7

                x2y2ae_f(1,i3,i1,i2)=ttae1
                x2y2ae_f(2,i3,i1,i2)=ttae2
                x2y2ae_f(3,i3,i1,i2)=ttae3
                x2y2ae_f(4,i3,i1,i2)=ttae4
                x2y2ae_f(5,i3,i1,i2)=ttae5
                x2y2ae_f(6,i3,i1,i2)=ttae6
                x2y2ae_f(7,i3,i1,i2)=ttae7

                x2y2ba_f(2,i3,i1,i2)=ttba1
                x2y2ba_f(3,i3,i1,i2)=ttba2
                x2y2ba_f(4,i3,i1,i2)=ttba3
                x2y2ba_f(5,i3,i1,i2)=ttba4
                x2y2ba_f(6,i3,i1,i2)=ttba5
                x2y2ba_f(7,i3,i1,i2)=ttba6
                x2y2ba_f(7,i3,i1,i2)=ttba7

                x2y2bb_f(1,i3,i1,i2)=ttbb1
                x2y2bb_f(2,i3,i1,i2)=ttbb2
                x2y2bb_f(3,i3,i1,i2)=ttbb3
                x2y2bb_f(4,i3,i1,i2)=ttbb4
                x2y2bb_f(5,i3,i1,i2)=ttbb5
                x2y2bb_f(6,i3,i1,i2)=ttbb6
                x2y2bb_f(7,i3,i1,i2)=ttbb7

                x2y2bc_f(1,i3,i1,i2)=ttbc1
                x2y2bc_f(2,i3,i1,i2)=ttbc2
                x2y2bc_f(3,i3,i1,i2)=ttbc3
                x2y2bc_f(4,i3,i1,i2)=ttbc4
                x2y2bc_f(5,i3,i1,i2)=ttbc5
                x2y2bc_f(6,i3,i1,i2)=ttbc6
                x2y2bc_f(7,i3,i1,i2)=ttbc7

                x2y2be_f(1,i3,i1,i2)=ttbe1
                x2y2be_f(2,i3,i1,i2)=ttbe2
                x2y2be_f(3,i3,i1,i2)=ttbe3
                x2y2be_f(4,i3,i1,i2)=ttbe4
                x2y2be_f(5,i3,i1,i2)=ttbe5
                x2y2be_f(6,i3,i1,i2)=ttbe6
                x2y2be_f(7,i3,i1,i2)=ttbe7

                x2y2ca_f(1,i3,i1,i2)=ttca1
                x2y2ca_f(2,i3,i1,i2)=ttca2
                x2y2ca_f(3,i3,i1,i2)=ttca3
                x2y2ca_f(4,i3,i1,i2)=ttca4
                x2y2ca_f(5,i3,i1,i2)=ttca5
                x2y2ca_f(6,i3,i1,i2)=ttca6
                x2y2ca_f(7,i3,i1,i2)=ttca7

                x2y2cb_f(1,i3,i1,i2)=ttcb1
                x2y2cb_f(2,i3,i1,i2)=ttcb2
                x2y2cb_f(3,i3,i1,i2)=ttcb3
                x2y2cb_f(4,i3,i1,i2)=ttcb4
                x2y2cb_f(5,i3,i1,i2)=ttcb5
                x2y2cb_f(6,i3,i1,i2)=ttcb6
                x2y2cb_f(7,i3,i1,i2)=ttcb7

                x2y2cc_f(1,i3,i1,i2)=ttcc1
                x2y2cc_f(2,i3,i1,i2)=ttcc2
                x2y2cc_f(3,i3,i1,i2)=ttcc3
                x2y2cc_f(4,i3,i1,i2)=ttcc4
                x2y2cc_f(5,i3,i1,i2)=ttcc5
                x2y2cc_f(6,i3,i1,i2)=ttcc6
                x2y2cc_f(7,i3,i1,i2)=ttcc7

                x2y2ce_f(1,i3,i1,i2)=ttce1
                x2y2ce_f(2,i3,i1,i2)=ttce2
                x2y2ce_f(3,i3,i1,i2)=ttce3
                x2y2ce_f(4,i3,i1,i2)=ttce4
                x2y2ce_f(5,i3,i1,i2)=ttce5
                x2y2ce_f(6,i3,i1,i2)=ttce6
                x2y2ce_f(7,i3,i1,i2)=ttce7

                x2y2ea_f(1,i3,i1,i2)=ttea1
                x2y2ea_f(2,i3,i1,i2)=ttea2
                x2y2ea_f(3,i3,i1,i2)=ttea3
                x2y2ea_f(4,i3,i1,i2)=ttea4
                x2y2ea_f(5,i3,i1,i2)=ttea5
                x2y2ea_f(6,i3,i1,i2)=ttea6
                x2y2ea_f(7,i3,i1,i2)=ttea7

                x2y2eb_f(1,i3,i1,i2)=tteb1
                x2y2eb_f(2,i3,i1,i2)=tteb2
                x2y2eb_f(3,i3,i1,i2)=tteb3
                x2y2eb_f(4,i3,i1,i2)=tteb4
                x2y2eb_f(5,i3,i1,i2)=tteb5
                x2y2eb_f(6,i3,i1,i2)=tteb6
                x2y2eb_f(7,i3,i1,i2)=tteb7

                x2y2ec_f(1,i3,i1,i2)=ttec1
                x2y2ec_f(2,i3,i1,i2)=ttec2
                x2y2ec_f(3,i3,i1,i2)=ttec3
                x2y2ec_f(4,i3,i1,i2)=ttec4
                x2y2ec_f(5,i3,i1,i2)=ttec5
                x2y2ec_f(6,i3,i1,i2)=ttec6
                x2y2ec_f(7,i3,i1,i2)=ttec7

                x2y2ee_f(1,i3,i1,i2)=ttee1
                x2y2ee_f(2,i3,i1,i2)=ttee2
                x2y2ee_f(3,i3,i1,i2)=ttee3
                x2y2ee_f(4,i3,i1,i2)=ttee4
                x2y2ee_f(5,i3,i1,i2)=ttee5
                x2y2ee_f(6,i3,i1,i2)=ttee6
                x2y2ee_f(7,i3,i1,i2)=ttee7
          enddo
       enddo
    enddo
    !!!$omp enddo




    do i3=0,n3
        z0=hgrid*(i3+offsetz)-rxyzConf(3)
        if(.not. withKinetic) then
            call getFilterSextic(it, potentialPrefac, hgrid, z0, aeff0array(lowfil,i3), 'a')
            call getFilterSextic(it, potentialPrefac, hgrid, z0, beff0array(lowfil,i3), 'b')
            call getFilterSextic(it, potentialPrefac, hgrid, z0, ceff0array(lowfil,i3), 'c')
            call getFilterSextic(it, potentialPrefac, hgrid, z0, eeff0array(lowfil,i3), 'e')
        else
            call getEffectiveFilterSextic(it, potentialPrefac, hgrid, z0, aeff0array(lowfil,i3), 'a')
            call getEffectiveFilterSextic(it, potentialPrefac, hgrid, z0, beff0array(lowfil,i3), 'b')
            call getEffectiveFilterSextic(it, potentialPrefac, hgrid, z0, ceff0array(lowfil,i3), 'c')
            call getEffectiveFilterSextic(it, potentialPrefac, hgrid, z0, eeff0array(lowfil,i3), 'e')
        end if

        call getFilterQuartic(potentialPrefac, hgrid, z0, aeff0_4array(lowfil,i3), 'a')
        call getFilterQuartic(potentialPrefac, hgrid, z0, beff0_4array(lowfil,i3), 'b')
        call getFilterQuartic(potentialPrefac, hgrid, z0, ceff0_4array(lowfil,i3), 'c')
        call getFilterQuartic(potentialPrefac, hgrid, z0, eeff0_4array(lowfil,i3), 'e')

        call getFilterQuadratic(potentialPrefac, hgrid, z0, aeff0_2array(lowfil,i3), 'a')
        call getFilterQuadratic(potentialPrefac, hgrid, z0, beff0_2array(lowfil,i3), 'b')
        call getFilterQuadratic(potentialPrefac, hgrid, z0, ceff0_2array(lowfil,i3), 'c')
        call getFilterQuadratic(potentialPrefac, hgrid, z0, eeff0_2array(lowfil,i3), 'e')
    end do


  ! + (1/2) d^2/dz^2

  !!!$omp do
  do i2=0,n2
     do i1=0,n1
        if (ibxy_c(2,i1,i2)-ibxy_c(1,i1,i2).ge.4) then
           do i3=ibxy_c(1,i1,i2),ibxy_c(2,i1,i2)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_c(2,i1,i2))
                 dyi0=dyi0 + xz_c(t,i1,i2)*aeff0array(t-i3-0,i3+0) + &
                             3.d0*( (x4za_c(t,i1,i2)+y4za_c(t,i1,i2))*aeff0_2array(t-i3-0,i3+0) + &
                             (x2za_c(t,i1,i2)+y2za_c(t,i1,i2))*aeff0_4array(t-i3-0,i3+0) ) + &
                             6.d0*x2y2aa_c(t,i1,i2)*aeff0_2array(t-i3-0,i3+0)
                 dyi1=dyi1 + xz_c(t,i1,i2)*aeff0array(t-i3-1,i3+1) + &
                             3.d0*( (x4za_c(t,i1,i2)+y4za_c(t,i1,i2))*aeff0_2array(t-i3-1,i3+1) + &
                             (x2za_c(t,i1,i2)+y2za_c(t,i1,i2))*aeff0_4array(t-i3-1,i3+1) ) + &
                             6.d0*x2y2aa_c(t,i1,i2)*aeff0_2array(t-i3-1,i3+1)
                 dyi2=dyi2 + xz_c(t,i1,i2)*aeff0array(t-i3-2,i3+2) + &
                             3.d0*( (x4za_c(t,i1,i2)+y4za_c(t,i1,i2))*aeff0_2array(t-i3-2,i3+2) + &
                             (x2za_c(t,i1,i2)+y2za_c(t,i1,i2))*aeff0_4array(t-i3-2,i3+2) ) + &
                             6.d0*x2y2aa_c(t,i1,i2)*aeff0_2array(t-i3-2,i3+2)
                 dyi3=dyi3 + xz_c(t,i1,i2)*aeff0array(t-i3-3,i3+3) + &
                             3.d0*( (x4za_c(t,i1,i2)+y4za_c(t,i1,i2))*aeff0_2array(t-i3-3,i3+3) + &
                             (x2za_c(t,i1,i2)+y2za_c(t,i1,i2))*aeff0_4array(t-i3-3,i3+3) ) + &
                             6.d0*x2y2aa_c(t,i1,i2)*aeff0_2array(t-i3-3,i3+3)
              enddo
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
           do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_c(2,i1,i2))
              dyi0=dyi0 + xz_c(t,i1,i2)*aeff0array(t-i3-0,i3+0) + &
                          3.d0*( (x4za_c(t,i1,i2)+y4za_c(t,i1,i2))*aeff0_2array(t-i3-0,i3+0) + &
                          (x2za_c(t,i1,i2)+y2za_c(t,i1,i2))*aeff0_4array(t-i3-0,i3+0) ) + & 
                          6.d0*x2y2aa_c(t,i1,i2)*aeff0_2array(t-i3-0,i3+0)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi0
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
                 dyi0 = dyi0 + xz_f4(t,i1,i2)*beff0array(t-i3-0,i3+0) + &
                               3.d0*(x4zb_f(1,t,i1,i2)+y4zb_f(1,t,i1,i2)                                    )* &
                                 aeff0_2array(t-i3-0,i3+0) + &
                               3.d0*(x4za_f(2,t,i1,i2)+x4zb_f(3,t,i1,i2)+y4za_f(2,t,i1,i2)+y4zb_f(3,t,i1,i2))* &
                                 beff0_2array(t-i3-0,i3+0) + &
                               3.d0*(x2zb_f(1,t,i1,i2)+y2zb_f(1,t,i1,i2)                                    )* &
                                 aeff0_4array(t-i3-0,i3+0) + &
                               3.d0*(x2za_f(2,t,i1,i2)+x2zb_f(3,t,i1,i2)+y2za_f(2,t,i1,i2)+y2zb_f(3,t,i1,i2))* &
                                 beff0_4array(t-i3-0,i3+0) + &
                               6.d0*(                    x2y2ba_f(1,t,i1,i2)+x2y2ab_f(2,t,i1,i2)+x2y2bb_f(3,t,i1,i2))* &
                                 aeff0_2array(t-i3-0,i3+0) + &
                               6.d0*(x2y2aa_f(4,t,i1,i2)+x2y2ba_f(5,t,i1,i2)+x2y2ab_f(6,t,i1,i2)+x2y2bb_f(7,t,i1,i2))* &
                                 beff0_2array(t-i3-0,i3+0)
                 dyi1 = dyi1 + xz_f4(t,i1,i2)*beff0array(t-i3-1,i3+1) + &
                               3.d0*(x4zb_f(1,t,i1,i2)+y4zb_f(1,t,i1,i2)                                    )* &
                                 aeff0_2array(t-i3-1,i3+1) + &
                               3.d0*(x4za_f(2,t,i1,i2)+x4zb_f(3,t,i1,i2)+y4za_f(2,t,i1,i2)+y4zb_f(3,t,i1,i2))* &
                                 beff0_2array(t-i3-1,i3+1) + &
                               3.d0*(x2zb_f(1,t,i1,i2)+y2zb_f(1,t,i1,i2)                                    )* &
                                 aeff0_4array(t-i3-1,i3+1) + &
                               3.d0*(x2za_f(2,t,i1,i2)+x2zb_f(3,t,i1,i2)+y2za_f(2,t,i1,i2)+y2zb_f(3,t,i1,i2))* &
                                 beff0_4array(t-i3-1,i3+1) + &
                               6.d0*(                    x2y2ba_f(1,t,i1,i2)+x2y2ab_f(2,t,i1,i2)+x2y2bb_f(3,t,i1,i2))* &
                                 aeff0_2array(t-i3-1,i3+1) + &
                               6.d0*(x2y2aa_f(4,t,i1,i2)+x2y2ba_f(5,t,i1,i2)+x2y2ab_f(6,t,i1,i2)+x2y2bb_f(7,t,i1,i2))* &
                                 beff0_2array(t-i3-1,i3+1)
                 dyi2 = dyi2 + xz_f4(t,i1,i2)*beff0array(t-i3-2,i3+2) + &
                               3.d0*(x4zb_f(1,t,i1,i2)+y4zb_f(1,t,i1,i2)                                    )* &
                                 aeff0_2array(t-i3-2,i3+2) + &
                               3.d0*(x4za_f(2,t,i1,i2)+x4zb_f(3,t,i1,i2)+y4za_f(2,t,i1,i2)+y4zb_f(3,t,i1,i2))* &
                                 beff0_2array(t-i3-2,i3+2) + &
                               3.d0*(x2zb_f(1,t,i1,i2)+y2zb_f(1,t,i1,i2)                                    )* &
                                 aeff0_4array(t-i3-2,i3+2) + &
                               3.d0*(x2za_f(2,t,i1,i2)+x2zb_f(3,t,i1,i2)+y2za_f(2,t,i1,i2)+y2zb_f(3,t,i1,i2))* &
                                 beff0_4array(t-i3-2,i3+2) + &
                               6.d0*(                    x2y2ba_f(1,t,i1,i2)+x2y2ab_f(2,t,i1,i2)+x2y2bb_f(3,t,i1,i2))* &
                                 aeff0_2array(t-i3-2,i3+2) + &
                               6.d0*(x2y2aa_f(4,t,i1,i2)+x2y2ba_f(5,t,i1,i2)+x2y2ab_f(6,t,i1,i2)+x2y2bb_f(7,t,i1,i2))* &
                                 beff0_2array(t-i3-2,i3+2)
                 dyi3 = dyi3 + xz_f4(t,i1,i2)*beff0array(t-i3-3,i3+3) + &
                               3.d0*(x4zb_f(1,t,i1,i2)+y4zb_f(1,t,i1,i2)                                    )* &
                                 aeff0_2array(t-i3-3,i3+3) + &
                               3.d0*(x4za_f(2,t,i1,i2)+x4zb_f(3,t,i1,i2)+y4za_f(2,t,i1,i2)+y4zb_f(3,t,i1,i2))* &
                                 beff0_2array(t-i3-3,i3+3) + &
                               3.d0*(x2zb_f(1,t,i1,i2)+y2zb_f(1,t,i1,i2)                                    )* &
                                 aeff0_4array(t-i3-3,i3+3) + &
                               3.d0*(x2za_f(2,t,i1,i2)+x2zb_f(3,t,i1,i2)+y2za_f(2,t,i1,i2)+y2zb_f(3,t,i1,i2))* &
                                 beff0_4array(t-i3-3,i3+3) + &
                               6.d0*(                    x2y2ba_f(1,t,i1,i2)+x2y2ab_f(2,t,i1,i2)+x2y2bb_f(3,t,i1,i2))* &
                                 aeff0_2array(t-i3-3,i3+3) + &
                               6.d0*(x2y2aa_f(4,t,i1,i2)+x2y2ba_f(5,t,i1,i2)+x2y2ab_f(6,t,i1,i2)+x2y2bb_f(7,t,i1,i2))* &
                                 beff0_2array(t-i3-3,i3+3)
              enddo
              y_c(i1,i2,i3+0)=y_c(i1,i2,i3+0)+dyi0
              y_c(i1,i2,i3+1)=y_c(i1,i2,i3+1)+dyi1
              y_c(i1,i2,i3+2)=y_c(i1,i2,i3+2)+dyi2
              y_c(i1,i2,i3+3)=y_c(i1,i2,i3+3)+dyi3
           enddo
           istart=i2
        endif

        do i3=istart,iend
           dyi0=0.0_wp
           tt0=0.0_wp
           do t=max(ibxy_f(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_f(2,i1,i2))
              dyi0 = dyi0 + xz_f4(t,i1,i2)*beff0array(t-i3-0,i3+0) + &
                            3.d0*(x4zb_f(1,t,i1,i2)+y4zb_f(1,t,i1,i2)                                    )* &
                              aeff0_2array(t-i3-0,i3+0) + &
                            3.d0*(x4za_f(2,t,i1,i2)+x4zb_f(3,t,i1,i2)+y4za_f(2,t,i1,i2)+y4zb_f(3,t,i1,i2))* &
                              beff0_2array(t-i3-0,i3+0) + &
                            3.d0*(x2zb_f(1,t,i1,i2)+y2zb_f(1,t,i1,i2)                                    )* &
                              aeff0_4array(t-i3-0,i3+0) + &
                            3.d0*(x2za_f(2,t,i1,i2)+x2zb_f(3,t,i1,i2)+y2za_f(2,t,i1,i2)+y2zb_f(3,t,i1,i2))* &
                              beff0_4array(t-i3-0,i3+0) + &
                            6.d0*(                    x2y2ba_f(1,t,i1,i2)+x2y2ab_f(2,t,i1,i2)+x2y2bb_f(3,t,i1,i2))* &
                              aeff0_2array(t-i3-0,i3+0) + &
                            6.d0*(x2y2aa_f(4,t,i1,i2)+x2y2ba_f(5,t,i1,i2)+x2y2ab_f(6,t,i1,i2)+x2y2bb_f(7,t,i1,i2))* &
                              beff0_2array(t-i3-0,i3+0)
           enddo
           y_c(i1,i2,i3)=y_c(i1,i2,i3)+dyi0
        enddo

         if (ibxy_f(2,i1,i2)-ibxy_f(1,i1,i2).ge.4) then
           do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)-4,4
              dyi0=0.0_wp 
              dyi1=0.0_wp 
              dyi2=0.0_wp 
              dyi3=0.0_wp 
              tt10 = 0.d0
              tt11 = 0.d0
              tt12 = 0.d0
              tt13 = 0.d0
              tt40 = 0.d0
              tt41 = 0.d0
              tt42 = 0.d0
              tt43 = 0.d0
              tt50 = 0.d0
              tt51 = 0.d0
              tt52 = 0.d0
              tt53 = 0.d0
              tt20 = 0.d0
              tt21 = 0.d0
              tt22 = 0.d0
              tt23 = 0.d0
              tt60 = 0.d0
              tt61 = 0.d0
              tt62 = 0.d0
              tt63 = 0.d0
              tt30 = 0.d0
              tt31 = 0.d0
              tt32 = 0.d0
              tt33 = 0.d0
              tt70 = 0.d0
              tt71 = 0.d0
              tt72 = 0.d0
              tt73 = 0.d0
              do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3+3,ibxy_c(2,i1,i2))
                 dyi0=dyi0 + xz_c(t,i1,i2)*ceff0array(t-i3-0,i3+0)
                 dyi1=dyi1 + xz_c(t,i1,i2)*ceff0array(t-i3-1,i3+1)
                 dyi2=dyi2 + xz_c(t,i1,i2)*ceff0array(t-i3-2,i3+2)
                 dyi3=dyi3 + xz_c(t,i1,i2)*ceff0array(t-i3-3,i3+3)

                 tt10 = tt10 + 3.d0*( x4zc_c(t,i1,i2)*aeff0_2array(t-i3-0,i3+0) + x2zc_c(t,i1,i2)*aeff0_4array(t-i3-0,i3+0) ) + &
                   6.d0*x2y2ca_c(t,i1,i2)*aeff0_2array(t-i3-0,i3+0)
                 tt11 = tt11 + 3.d0*( x4zc_c(t,i1,i2)*aeff0_2array(t-i3-1,i3+1) + x2zc_c(t,i1,i2)*aeff0_4array(t-i3-1,i3+1) ) + &
                   6.d0*x2y2ca_c(t,i1,i2)*aeff0_2array(t-i3-1,i3+1)
                 tt12 = tt12 + 3.d0*( x4zc_c(t,i1,i2)*aeff0_2array(t-i3-2,i3+2) + x2zc_c(t,i1,i2)*aeff0_4array(t-i3-2,i3+2) ) + &
                   6.d0*x2y2ca_c(t,i1,i2)*aeff0_2array(t-i3-2,i3+2)
                 tt13 = tt13 + 3.d0*( x4zc_c(t,i1,i2)*aeff0_2array(t-i3-3,i3+3) + x2zc_c(t,i1,i2)*aeff0_4array(t-i3-3,i3+3) ) + &
                   6.d0*x2y2ca_c(t,i1,i2)*aeff0_2array(t-i3-3,i3+3)

                 tt40 = tt40 + 3.d0*( x4za_c(t,i1,i2)*ceff0_2array(t-i3-0,i3+0) + x2za_c(t,i1,i2)*ceff0_4array(t-i3-0,i3+0) ) + &
                   6.d0*x2y2aa_c(t,i1,i2)*ceff0_2array(t-i3-0,i3+0)
                 tt41 = tt41 + 3.d0*( x4za_c(t,i1,i2)*ceff0_2array(t-i3-1,i3+1) + x2za_c(t,i1,i2)*ceff0_4array(t-i3-1,i3+1) ) + &
                   6.d0*x2y2aa_c(t,i1,i2)*ceff0_2array(t-i3-1,i3+1)
                 tt42 = tt42 + 3.d0*( x4za_c(t,i1,i2)*ceff0_2array(t-i3-2,i3+2) + x2za_c(t,i1,i2)*ceff0_4array(t-i3-2,i3+2) ) + &
                   6.d0*x2y2aa_c(t,i1,i2)*ceff0_2array(t-i3-2,i3+2)
                 tt43 = tt43 + 3.d0*( x4za_c(t,i1,i2)*ceff0_2array(t-i3-3,i3+3) + x2za_c(t,i1,i2)*ceff0_4array(t-i3-3,i3+3) ) + &
                   6.d0*x2y2aa_c(t,i1,i2)*ceff0_2array(t-i3-3,i3+3)

                 tt50 = tt50 + 3.d0*( x4zc_c(t,i1,i2)*ceff0_2array(t-i3-0,i3+0) + x2zc_c(t,i1,i2)*ceff0_4array(t-i3-0,i3+0) ) + &
                   6.d0*x2y2ca_c(t,i1,i2)*ceff0_2array(t-i3-0,i3+0)
                 tt51 = tt51 + 3.d0*( x4zc_c(t,i1,i2)*ceff0_2array(t-i3-1,i3+1) + x2zc_c(t,i1,i2)*ceff0_4array(t-i3-1,i3+1) ) + &
                   6.d0*x2y2ca_c(t,i1,i2)*ceff0_2array(t-i3-1,i3+1)
                 tt52 = tt52 + 3.d0*( x4zc_c(t,i1,i2)*ceff0_2array(t-i3-2,i3+2) + x2zc_c(t,i1,i2)*ceff0_4array(t-i3-2,i3+2) ) + &
                   6.d0*x2y2ca_c(t,i1,i2)*ceff0_2array(t-i3-2,i3+2)
                 tt53 = tt53 + 3.d0*( x4zc_c(t,i1,i2)*ceff0_2array(t-i3-3,i3+3) + x2zc_c(t,i1,i2)*ceff0_4array(t-i3-3,i3+3) ) + &
                   6.d0*x2y2ca_c(t,i1,i2)*ceff0_2array(t-i3-3,i3+3)

                 tt20 = tt20 + 3.d0*( y4zc_c(t,i1,i2)*aeff0_2array(t-i3-0,i3+0) + y2zc_c(t,i1,i2)*aeff0_4array(t-i3-0,i3+0) ) + &
                   6.d0*x2y2ac_c(t,i1,i2)*aeff0_2array(t-i3-0,i3+0)
                 tt21 = tt21 + 3.d0*( y4zc_c(t,i1,i2)*aeff0_2array(t-i3-1,i3+1) + y2zc_c(t,i1,i2)*aeff0_4array(t-i3-1,i3+1) ) + &
                   6.d0*x2y2ac_c(t,i1,i2)*aeff0_2array(t-i3-1,i3+1)
                 tt22 = tt22 + 3.d0*( y4zc_c(t,i1,i2)*aeff0_2array(t-i3-2,i3+2) + y2zc_c(t,i1,i2)*aeff0_4array(t-i3-2,i3+2) ) + &
                   6.d0*x2y2ac_c(t,i1,i2)*aeff0_2array(t-i3-2,i3+2)
                 tt23 = tt23 + 3.d0*( y4zc_c(t,i1,i2)*aeff0_2array(t-i3-3,i3+3) + y2zc_c(t,i1,i2)*aeff0_4array(t-i3-3,i3+3) ) + &
                   6.d0*x2y2ac_c(t,i1,i2)*aeff0_2array(t-i3-3,i3+3)

                 tt40 = tt40 + 3.d0*( y4za_c(t,i1,i2)*ceff0_2array(t-i3-0,i3+0) + y2za_c(t,i1,i2)*ceff0_4array(t-i3-0,i3+0) )
                 tt41 = tt41 + 3.d0*( y4za_c(t,i1,i2)*ceff0_2array(t-i3-1,i3+1) + y2za_c(t,i1,i2)*ceff0_4array(t-i3-1,i3+1) )
                 tt42 = tt42 + 3.d0*( y4za_c(t,i1,i2)*ceff0_2array(t-i3-2,i3+2) + y2za_c(t,i1,i2)*ceff0_4array(t-i3-2,i3+2) )
                 tt43 = tt43 + 3.d0*( y4za_c(t,i1,i2)*ceff0_2array(t-i3-3,i3+3) + y2za_c(t,i1,i2)*ceff0_4array(t-i3-3,i3+3) )

                 tt60 = tt60 + 3.d0*( y4zc_c(t,i1,i2)*ceff0_2array(t-i3-0,i3+0) + y2zc_c(t,i1,i2)*ceff0_4array(t-i3-0,i3+0) ) + &
                   6.d0*x2y2ac_c(t,i1,i2)*ceff0_2array(t-i3-0,i3+0)
                 tt61 = tt61 + 3.d0*( y4zc_c(t,i1,i2)*ceff0_2array(t-i3-1,i3+1) + y2zc_c(t,i1,i2)*ceff0_4array(t-i3-1,i3+1) ) + &
                   6.d0*x2y2ac_c(t,i1,i2)*ceff0_2array(t-i3-1,i3+1)
                 tt62 = tt62 + 3.d0*( y4zc_c(t,i1,i2)*ceff0_2array(t-i3-2,i3+2) + y2zc_c(t,i1,i2)*ceff0_4array(t-i3-2,i3+2) ) + &
                   6.d0*x2y2ac_c(t,i1,i2)*ceff0_2array(t-i3-2,i3+2)
                 tt63 = tt63 + 3.d0*( y4zc_c(t,i1,i2)*ceff0_2array(t-i3-3,i3+3) + y2zc_c(t,i1,i2)*ceff0_4array(t-i3-3,i3+3) ) + &
                   6.d0*x2y2ac_c(t,i1,i2)*ceff0_2array(t-i3-3,i3+3)

                 tt30 = tt30 + 6.d0*x2y2cc_c(t,i1,i2)*aeff0_2array(t-i3-0,i3+0)
                 tt31 = tt31 + 6.d0*x2y2cc_c(t,i1,i2)*aeff0_2array(t-i3-1,i3+1)
                 tt32 = tt32 + 6.d0*x2y2cc_c(t,i1,i2)*aeff0_2array(t-i3-2,i3+2)
                 tt33 = tt33 + 6.d0*x2y2cc_c(t,i1,i2)*aeff0_2array(t-i3-3,i3+3)

                 tt70 = tt70 + 6.d0*x2y2cc_c(t,i1,i2)*ceff0_2array(t-i3-0,i3+0)
                 tt71 = tt71 + 6.d0*x2y2cc_c(t,i1,i2)*ceff0_2array(t-i3-1,i3+1)
                 tt72 = tt72 + 6.d0*x2y2cc_c(t,i1,i2)*ceff0_2array(t-i3-2,i3+2)
                 tt73 = tt73 + 6.d0*x2y2cc_c(t,i1,i2)*ceff0_2array(t-i3-3,i3+3)

              enddo
              y_f(4,i1,i2,i3+0) = y_f(4,i1,i2,i3+0) + dyi0 + tt40
              y_f(4,i1,i2,i3+1) = y_f(4,i1,i2,i3+1) + dyi1 + tt41
              y_f(4,i1,i2,i3+2) = y_f(4,i1,i2,i3+2) + dyi2 + tt42
              y_f(4,i1,i2,i3+3) = y_f(4,i1,i2,i3+3) + dyi3 + tt43

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

              y_f(3,i1,i2,i3+0) = y_f(3,i1,i2,i3+0) + tt30
              y_f(3,i1,i2,i3+1) = y_f(3,i1,i2,i3+1) + tt31
              y_f(3,i1,i2,i3+2) = y_f(3,i1,i2,i3+2) + tt32
              y_f(3,i1,i2,i3+3) = y_f(3,i1,i2,i3+3) + tt33

              y_f(7,i1,i2,i3+0) = y_f(7,i1,i2,i3+0) + tt70
              y_f(7,i1,i2,i3+1) = y_f(7,i1,i2,i3+1) + tt71
              y_f(7,i1,i2,i3+2) = y_f(7,i1,i2,i3+2) + tt72
              y_f(7,i1,i2,i3+3) = y_f(7,i1,i2,i3+3) + tt73
           enddo
           icur=i3
        else
           icur=ibxy_f(1,i1,i2)
        endif

        do i3=icur,ibxy_f(2,i1,i2)
           dyi0=0.0_wp 
           tt10 = 0.d0
           tt40 = 0.d0
           tt50 = 0.d0
           tt20 = 0.d0
           tt60 = 0.d0
           do t=max(ibxy_c(1,i1,i2),lowfil+i3),min(lupfil+i3,ibxy_c(2,i1,i2))
              dyi0=dyi0 + xz_c(t,i1,i2)*ceff0array(t-i3-0,i3)

              tt10 = tt10 + 3.d0*( x4zc_c(t,i1,i2)*aeff0_2array(t-i3-0,i3+0) + x2zc_c(t,i1,i2)*aeff0_4array(t-i3-0,i3+0) ) + &
                6.d0*x2y2ca_c(t,i1,i2)*aeff0_2array(t-i3-0,i3+0)

              tt40 = tt40 + 3.d0*( x4za_c(t,i1,i2)*ceff0_2array(t-i3-0,i3+0) + x2za_c(t,i1,i2)*ceff0_4array(t-i3-0,i3+0) ) + &
                6.d0*x2y2aa_c(t,i1,i2)*ceff0_2array(t-i3-0,i3+0)

              tt50 = tt50 + 3.d0*( x4zc_c(t,i1,i2)*ceff0_2array(t-i3-0,i3+0) + x2zc_c(t,i1,i2)*ceff0_4array(t-i3-0,i3+0) ) + &
                6.d0*x2y2ca_c(t,i1,i2)*ceff0_2array(t-i3-0,i3+0)

              tt20 = tt20 + 3.d0*( y4zc_c(t,i1,i2)*aeff0_2array(t-i3-0,i3+0) + y2zc_c(t,i1,i2)*aeff0_4array(t-i3-0,i3+0) ) + &
                6.d0*x2y2ac_c(t,i1,i2)*aeff0_2array(t-i3-0,i3+0)

              tt40 = tt40 + 3.d0*( y4za_c(t,i1,i2)*ceff0_2array(t-i3-0,i3+0) + y2za_c(t,i1,i2)*ceff0_4array(t-i3-0,i3+0) )

              tt60 = tt60 + 3.d0*( y4zc_c(t,i1,i2)*ceff0_2array(t-i3-0,i3+0) + y2zc_c(t,i1,i2)*ceff0_4array(t-i3-0,i3+0) ) + &
                6.d0*x2y2ac_c(t,i1,i2)*ceff0_2array(t-i3-0,i3+0)

              tt30 = tt30 + 6.d0*x2y2cc_c(t,i1,i2)*aeff0_2array(t-i3-0,i3+0)

              tt70 = tt70 + 6.d0*x2y2cc_c(t,i1,i2)*ceff0_2array(t-i3-0,i3+0)

           enddo
           y_f(4,i1,i2,i3+0) = y_f(4,i1,i2,i3+0) + dyi0 + tt40
           y_f(1,i1,i2,i3+0) = y_f(1,i1,i2,i3+0) + tt10
           y_f(5,i1,i2,i3+0) = y_f(5,i1,i2,i3+0) + tt50
           y_f(2,i1,i2,i3+0) = y_f(2,i1,i2,i3+0) + tt20
           y_f(6,i1,i2,i3+0) = y_f(6,i1,i2,i3+0) + tt60
           y_f(3,i1,i2,i3+0) = y_f(3,i1,i2,i3+0) + tt30
           y_f(7,i1,i2,i3+0) = y_f(7,i1,i2,i3+0) + tt70
        enddo
     enddo
  enddo
  !!!$omp enddo


  ! wavelet part

  !!!$omp do
  do i2=nfl2,nfu2
     do i1=nfl1,nfu1
        do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
           tt10 = 0.d0
           tt20 = 0.d0
           tt30 = 0.d0
           tt40 = 0.d0
           tt50 = 0.d0
           tt60 = 0.d0
           tt70 = 0.d0
           do l=max(nfl3-i3,lowfil),min(lupfil,nfu3-i3)

              tt10 = tt10 + xz_f(1,i3+l,i1,i2)*aeff0array(l,i3) + xz_f(5,i3+l,i1,i2)*beff0array(l,i3) + &
                            3.d0*                      x4ze_f(1,i3+l,i1,i2) *aeff0_2array(l,i3) + &
                            3.d0*                      x2ze_f(1,i3+l,i1,i2) *aeff0_4array(l,i3) + &
                            3.d0*(x4zc_f(2,i3+l,i1,i2)+x4ze_f(3,i3+l,i1,i2))*beff0_2array(l,i3) + &
                            3.d0*(x2zc_f(2,i3+l,i1,i2)+x4ze_f(3,i3+l,i1,i2))*beff0_4array(l,i3) + &
                            3.d0*(y4za_f(1,i3+l,i1,i2)+y4zb_f(2,i3+l,i1,i2))*aeff0_2array(l,i3) + &
                            3.d0*(y2za_f(1,i3+l,i1,i2)+y2zb_f(2,i3+l,i1,i2))*aeff0_4array(l,i3) + &
                            3.d0*(y4za_f(3,i3+l,i1,i2)+y4zb_f(4,i3+l,i1,i2))*beff0_2array(l,i3) + &
                            3.d0*(y2za_f(3,i3+l,i1,i2)+y2zb_f(4,i3+l,i1,i2))*beff0_4array(l,i3) + &
                            6.d0*(                       x2y2ea_f(1,i3+l,i1,i2)+x2y2cb_f(1,i3+l,i1,i2)+x2y2eb_f(1,i3+l,i1,i2))* &
                              aeff0_2array(l,i3) + &
                            6.d0*(x2y2ca_f(4,i3+l,i1,i2)+x2y2ea_f(5,i3+l,i1,i2)+x2y2cb_f(6,i3+l,i1,i2)+x2y2eb_f(7,i3+l,i1,i2))* &
                              beff0_2array(l,i3) 

              tt20 = tt20 + xz_f(2,i3+l,i1,i2)*aeff0array(l,i3) + xz_f(6,i3+l,i1,i2)*beff0array(l,i3) + &
                            3.d0*(x4za_f(1,i3+l,i1,i2)+x4zb_f(2,i3+l,i1,i2))*aeff0_2array(l,i3) + &
                            3.d0*(x2za_f(1,i3+l,i1,i2)+x2zb_f(2,i3+l,i1,i2))*aeff0_4array(l,i3) + &
                            3.d0*(x4za_f(3,i3+l,i1,i2)+x4zb_f(4,i3+l,i1,i2))*beff0_2array(l,i3) + &
                            3.d0*(x2za_f(3,i3+l,i1,i2)+x2zb_f(4,i3+l,i1,i2))*beff0_4array(l,i3) + &
                            3.d0*                      y4ze_f(1,i3+l,i1,i2) *aeff0_2array(l,i3) + &
                            3.d0*                      y2ze_f(1,i3+l,i1,i2) *aeff0_4array(l,i3) + &
                            3.d0*(y4zc_f(2,i3+l,i1,i2)+y4ze_f(3,i3+l,i1,i2))*beff0_2array(l,i3) + &
                            3.d0*(y2zc_f(2,i3+l,i1,i2)+y2ze_f(3,i3+l,i1,i2))*beff0_4array(l,i3) + &
                            6.d0*(                       x2y2bc_f(1,i3+l,i1,i2)+x2y2ae_f(1,i3+l,i1,i2)+x2y2be_f(1,i3+l,i1,i2))* &
                              aeff0_2array(l,i3) + &
                            6.d0*(x2y2ac_f(4,i3+l,i1,i2)+x2y2bc_f(5,i3+l,i1,i2)+x2y2ae_f(6,i3+l,i1,i2)+x2y2be_f(7,i3+l,i1,i2))* &
                              beff0_2array(l,i3) 

              tt30 = tt30 + xz_f(3,i3+l,i1,i2)*aeff0array(l,i3) + xz_f(7,i3+l,i1,i2)*beff0array(l,i3) + &
                            3.d0*(x4zc_f(1,i3+l,i1,i2)+x4ze_f(2,i3+l,i1,i2))*aeff0_2array(l,i3) + &
                            3.d0*(x2zc_f(1,i3+l,i1,i2)+x2ze_f(2,i3+l,i1,i2))*aeff0_4array(l,i3) + &
                            3.d0*(x4zc_f(3,i3+l,i1,i2)+x4ze_f(4,i3+l,i1,i2))*beff0_2array(l,i3) + &
                            3.d0*(x2zc_f(3,i3+l,i1,i2)+x2ze_f(4,i3+l,i1,i2))*beff0_4array(l,i3) + &
                            3.d0*(y4zc_f(1,i3+l,i1,i2)+y4ze_f(2,i3+l,i1,i2))*aeff0_2array(l,i3) + &
                            3.d0*(y2zc_f(1,i3+l,i1,i2)+y2ze_f(2,i3+l,i1,i2))*aeff0_4array(l,i3) + &
                            3.d0*(y4zc_f(3,i3+l,i1,i2)+y4ze_f(4,i3+l,i1,i2))*beff0_2array(l,i3) + &
                            3.d0*(y2zc_f(3,i3+l,i1,i2)+y2ze_f(4,i3+l,i1,i2))*beff0_4array(l,i3) + &
                            6.d0*(                       x2y2ec_f(1,i3+l,i1,i2)+x2y2ce_f(1,i3+l,i1,i2)+x2y2ee_f(1,i3+l,i1,i2))* &
                              aeff0_2array(l,i3) + &
                            6.d0*(x2y2cc_f(4,i3+l,i1,i2)+x2y2ec_f(5,i3+l,i1,i2)+x2y2ce_f(6,i3+l,i1,i2)+x2y2ee_f(7,i3+l,i1,i2))* &
                              beff0_2array(l,i3) 

              tt40 = tt40 + xz_f(4,i3+l,i1,i2)*eeff0array(l,i3)                                      + &
                            3.d0*                      x4zb_f(1,i3+l,i1,i2) *ceff0_2array(l,i3) + &
                            3.d0*                      x2zb_f(1,i3+l,i1,i2) *ceff0_4array(l,i3) + &
                            3.d0*(x4za_f(2,i3+l,i1,i2)+x4zb_f(3,i3+l,i1,i2))*eeff0_2array(l,i3) + &
                            3.d0*(x2za_f(2,i3+l,i1,i2)+x2zb_f(3,i3+l,i1,i2))*eeff0_4array(l,i3) + &
                            3.d0*                      y4zb_f(1,i3+l,i1,i2) *ceff0_2array(l,i3) + &
                            3.d0*                      y2zb_f(1,i3+l,i1,i2) *ceff0_4array(l,i3) + &
                            3.d0*(y4za_f(2,i3+l,i1,i2)+y4zb_f(3,i3+l,i1,i2))*eeff0_2array(l,i3) + &
                            3.d0*(y2za_f(2,i3+l,i1,i2)+y2zb_f(3,i3+l,i1,i2))*eeff0_4array(l,i3) + &
                            6.d0*(                       x2y2ba_f(1,i3+l,i1,i2)+x2y2ab_f(1,i3+l,i1,i2)+x2y2bb_f(1,i3+l,i1,i2))* &
                              ceff0_2array(l,i3) + &
                            6.d0*(x2y2aa_f(4,i3+l,i1,i2)+x2y2ba_f(5,i3+l,i1,i2)+x2y2ab_f(6,i3+l,i1,i2)+x2y2bb_f(7,i3+l,i1,i2))* &
                              eeff0_2array(l,i3) 

              tt50 = tt50 + xz_f(1,i3+l,i1,i2)*ceff0array(l,i3) + xz_f(5,i3+l,i1,i2)*eeff0array(l,i3) + &
                            3.d0*                      x4ze_f(1,i3+l,i1,i2) *ceff0_2array(l,i3) + &
                            3.d0*                      x2ze_f(1,i3+l,i1,i2) *ceff0_4array(l,i3) + &
                            3.d0*(x4zc_f(2,i3+l,i1,i2)+x4ze_f(3,i3+l,i1,i2))*eeff0_2array(l,i3) + &
                            3.d0*(x2zc_f(2,i3+l,i1,i2)+x2ze_f(3,i3+l,i1,i2))*eeff0_4array(l,i3) + &
                            3.d0*(y4za_f(1,i3+l,i1,i2)+y4zb_f(2,i3+l,i1,i2))*ceff0_2array(l,i3) + &
                            3.d0*(y2za_f(1,i3+l,i1,i2)+y2zb_f(2,i3+l,i1,i2))*ceff0_4array(l,i3) + &
                            3.d0*(y4za_f(3,i3+l,i1,i2)+y4zb_f(4,i3+l,i1,i2))*eeff0_2array(l,i3) + &
                            3.d0*(y2za_f(3,i3+l,i1,i2)+y2zb_f(4,i3+l,i1,i2))*eeff0_4array(l,i3) + &
                            6.d0*(                       x2y2ea_f(1,i3+l,i1,i2)+x2y2cb_f(1,i3+l,i1,i2)+x2y2eb_f(1,i3+l,i1,i2))* &
                              ceff0_2array(l,i3) + &
                            6.d0*(x2y2ca_f(4,i3+l,i1,i2)+x2y2ea_f(5,i3+l,i1,i2)+x2y2cb_f(6,i3+l,i1,i2)+x2y2eb_f(7,i3+l,i1,i2))* &
                              eeff0_2array(l,i3) 

              tt60 = tt60 + xz_f(2,i3+l,i1,i2)*ceff0array(l,i3) + xz_f(6,i3+l,i1,i2)*eeff0array(l,i3) + &
                            3.d0*(x4za_f(1,i3+l,i1,i2)+x4zb_f(2,i3+l,i1,i2))*ceff0_2array(l,i3) + &
                            3.d0*(x2za_f(1,i3+l,i1,i2)+x2zb_f(2,i3+l,i1,i2))*ceff0_4array(l,i3) + &
                            3.d0*(x4za_f(3,i3+l,i1,i2)+x4zb_f(4,i3+l,i1,i2))*eeff0_2array(l,i3) + &
                            3.d0*(x2za_f(3,i3+l,i1,i2)+x2zb_f(4,i3+l,i1,i2))*eeff0_4array(l,i3) + &
                            3.d0*                      y4ze_f(1,i3+l,i1,i2) *ceff0_2array(l,i3) + &
                            3.d0*                      y2ze_f(1,i3+l,i1,i2) *ceff0_4array(l,i3) + &
                            3.d0*(y4zc_f(2,i3+l,i1,i2)+y4ze_f(3,i3+l,i1,i2))*eeff0_2array(l,i3) + &
                            3.d0*(y2zc_f(2,i3+l,i1,i2)+y2ze_f(3,i3+l,i1,i2))*eeff0_4array(l,i3) + &
                            6.d0*(                       x2y2bc_f(1,i3+l,i1,i2)+x2y2ae_f(1,i3+l,i1,i2)+x2y2be_f(1,i3+l,i1,i2))* &
                              ceff0_2array(l,i3) + &
                            6.d0*(x2y2ac_f(4,i3+l,i1,i2)+x2y2bc_f(5,i3+l,i1,i2)+x2y2ae_f(6,i3+l,i1,i2)+x2y2be_f(7,i3+l,i1,i2))* &
                              eeff0_2array(l,i3) 

              tt70 = tt70 + xz_f(3,i3+l,i1,i2)*ceff0array(l,i3) + xz_f(7,i3+l,i1,i2)*eeff0array(l,i3) + &
                            3.d0*(x4zc_f(1,i3+l,i1,i2)+x4ze_f(2,i3+l,i1,i2))*ceff0_2array(l,i3) + &
                            3.d0*(x2zc_f(1,i3+l,i1,i2)+x2ze_f(2,i3+l,i1,i2))*ceff0_4array(l,i3) + &
                            3.d0*(x4zc_f(3,i3+l,i1,i2)+x4ze_f(4,i3+l,i1,i2))*eeff0_2array(l,i3) + &
                            3.d0*(x2zc_f(3,i3+l,i1,i2)+x2ze_f(4,i3+l,i1,i2))*eeff0_4array(l,i3) + &
                            3.d0*(y4zc_f(1,i3+l,i1,i2)+y4ze_f(2,i3+l,i1,i2))*ceff0_2array(l,i3) + &
                            3.d0*(y2zc_f(1,i3+l,i1,i2)+y2ze_f(2,i3+l,i1,i2))*ceff0_4array(l,i3) + &
                            3.d0*(y4zc_f(3,i3+l,i1,i2)+y4ze_f(4,i3+l,i1,i2))*eeff0_2array(l,i3) + &
                            3.d0*(y2zc_f(3,i3+l,i1,i2)+y2ze_f(4,i3+l,i1,i2))*eeff0_4array(l,i3) + &
                            6.d0*(                       x2y2ec_f(1,i3+l,i1,i2)+x2y2ce_f(1,i3+l,i1,i2)+x2y2ee_f(1,i3+l,i1,i2))* &
                              ceff0_2array(l,i3) + &
                            6.d0*(x2y2cc_f(4,i3+l,i1,i2)+x2y2ec_f(5,i3+l,i1,i2)+x2y2ce_f(6,i3+l,i1,i2)+x2y2ee_f(7,i3+l,i1,i2))* &
                              eeff0_2array(l,i3) 
           enddo
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
  !!!$omp enddo






  iall=-product(shape(aeff0array))*kind(aeff0array)
  deallocate(aeff0array, stat=istat)
  call memocc(istat, iall, 'aeff0array', subname)

  iall=-product(shape(beff0array))*kind(beff0array)
  deallocate(beff0array, stat=istat)
  call memocc(istat, iall, 'beff0array', subname)

  iall=-product(shape(ceff0array))*kind(ceff0array)
  deallocate(ceff0array, stat=istat)
  call memocc(istat, iall, 'ceff0array', subname)

  iall=-product(shape(eeff0array))*kind(eeff0array)
  deallocate(eeff0array, stat=istat)
  call memocc(istat, iall, 'eeff0array', subname)


  iall=-product(shape(aeff0_2array))*kind(aeff0_2array)
  deallocate(aeff0_2array, stat=istat)
  call memocc(istat, iall, 'aeff0_2array', subname)

  iall=-product(shape(beff0_2array))*kind(beff0_2array)
  deallocate(beff0_2array, stat=istat)
  call memocc(istat, iall, 'beff0_2array', subname)

  iall=-product(shape(ceff0_2array))*kind(ceff0_2array)
  deallocate(ceff0_2array, stat=istat)
  call memocc(istat, iall, 'ceff0_2array', subname)

  iall=-product(shape(eeff0_2array))*kind(eeff0_2array)
  deallocate(eeff0_2array, stat=istat)
  call memocc(istat, iall, 'eeff0_2array', subname)


  iall=-product(shape(aeff0_4array))*kind(aeff0_4array)
  deallocate(aeff0_4array, stat=istat)
  call memocc(istat, iall, 'aeff0_4array', subname)

  iall=-product(shape(beff0_4array))*kind(beff0_4array)
  deallocate(beff0_4array, stat=istat)
  call memocc(istat, iall, 'beff0_4array', subname)

  iall=-product(shape(ceff0_4array))*kind(ceff0_4array)
  deallocate(ceff0_4array, stat=istat)
  call memocc(istat, iall, 'ceff0_4array', subname)

  iall=-product(shape(eeff0_4array))*kind(eeff0_4array)
  deallocate(eeff0_4array, stat=istat)
  call memocc(istat, iall, 'eeff0_4array', subname)


  iall=-product(shape(aeff0_2auxarray))*kind(aeff0_2auxarray)
  deallocate(aeff0_2auxarray, stat=istat)
  call memocc(istat, iall, 'aeff0_2auxarray', subname)

  iall=-product(shape(beff0_2auxarray))*kind(beff0_2auxarray)
  deallocate(beff0_2auxarray, stat=istat)
  call memocc(istat, iall, 'beff0_2auxarray', subname)

  iall=-product(shape(ceff0_2auxarray))*kind(ceff0_2auxarray)
  deallocate(ceff0_2auxarray, stat=istat)
  call memocc(istat, iall, 'ceff0_2auxarray', subname)

  iall=-product(shape(eeff0_2auxarray))*kind(eeff0_2auxarray)
  deallocate(eeff0_2auxarray, stat=istat)
  call memocc(istat, iall, 'eeff0_2auxarray', subname)


  iall=-product(shape(aeff0_4auxarray))*kind(aeff0_4auxarray)
  deallocate(aeff0_4auxarray, stat=istat)
  call memocc(istat, iall, 'aeff0_4auxarray', subname)

  iall=-product(shape(beff0_4auxarray))*kind(beff0_4auxarray)
  deallocate(beff0_4auxarray, stat=istat)
  call memocc(istat, iall, 'beff0_4auxarray', subname)

  iall=-product(shape(ceff0_4auxarray))*kind(ceff0_4auxarray)
  deallocate(ceff0_4auxarray, stat=istat)
  call memocc(istat, iall, 'ceff0_4auxarray', subname)

  iall=-product(shape(eeff0_4auxarray))*kind(eeff0_4auxarray)
  deallocate(eeff0_4auxarray, stat=istat)
  call memocc(istat, iall, 'eeff0_4auxarray', subname)


  iall=-product(shape(x2ya_c))*kind(x2ya_c)
  deallocate(x2ya_c, stat=istat)
  call memocc(istat, iall, 'x2ya_c', subname)

  iall=-product(shape(x2yb_c))*kind(x2yb_c)
  deallocate(x2yb_c, stat=istat)
  call memocc(istat, iall, 'x2yb_c', subname)

  iall=-product(shape(x2yc_c))*kind(x2yc_c)
  deallocate(x2yc_c, stat=istat)
  call memocc(istat, iall, 'x2yc_c', subname)

  iall=-product(shape(x2ye_c))*kind(x2ye_c)
  deallocate(x2ye_c, stat=istat)
  call memocc(istat, iall, 'x2ye_c', subname)



  iall=-product(shape(x2za_c))*kind(x2za_c)
  deallocate(x2za_c, stat=istat)
  call memocc(istat, iall, 'x2za_c', subname)

  iall=-product(shape(x2zb_c))*kind(x2zb_c)
  deallocate(x2zb_c, stat=istat)
  call memocc(istat, iall, 'x2zb_c', subname)

  iall=-product(shape(x2zc_c))*kind(x2zc_c)
  deallocate(x2zc_c, stat=istat)
  call memocc(istat, iall, 'x2zc_c', subname)

  iall=-product(shape(x2ze_c))*kind(x2ze_c)
  deallocate(x2ze_c, stat=istat)
  call memocc(istat, iall, 'x2ze_c', subname)



  iall=-product(shape(y2za_c))*kind(y2za_c)
  deallocate(y2za_c, stat=istat)
  call memocc(istat, iall, 'y2za_c', subname)

  iall=-product(shape(y2zb_c))*kind(y2zb_c)
  deallocate(y2zb_c, stat=istat)
  call memocc(istat, iall, 'y2zb_c', subname)

  iall=-product(shape(y2zc_c))*kind(y2zc_c)
  deallocate(y2zc_c, stat=istat)
  call memocc(istat, iall, 'y2zc_c', subname)

  iall=-product(shape(y2ze_c))*kind(y2ze_c)
  deallocate(y2ze_c, stat=istat)
  call memocc(istat, iall, 'y2ze_c', subname)



  iall=-product(shape(x4ya_c))*kind(x4ya_c)
  deallocate(x4ya_c, stat=istat)
  call memocc(istat, iall, 'x4ya_c', subname)

  iall=-product(shape(x4yb_c))*kind(x4yb_c)
  deallocate(x4yb_c, stat=istat)
  call memocc(istat, iall, 'x4yb_c', subname)

  iall=-product(shape(x4yc_c))*kind(x4yc_c)
  deallocate(x4yc_c, stat=istat)
  call memocc(istat, iall, 'x4yc_c', subname)

  iall=-product(shape(x4ye_c))*kind(x4ye_c)
  deallocate(x4ye_c, stat=istat)
  call memocc(istat, iall, 'x4ye_c', subname)



  iall=-product(shape(x4za_c))*kind(x4za_c)
  deallocate(x4za_c, stat=istat)
  call memocc(istat, iall, 'x4za_c', subname)

  iall=-product(shape(x4zb_c))*kind(x4zb_c)
  deallocate(x4zb_c, stat=istat)
  call memocc(istat, iall, 'x4zb_c', subname)

  iall=-product(shape(x4zc_c))*kind(x4zc_c)
  deallocate(x4zc_c, stat=istat)
  call memocc(istat, iall, 'x4zc_c', subname)

  iall=-product(shape(x4ze_c))*kind(x4ze_c)
  deallocate(x4ze_c, stat=istat)
  call memocc(istat, iall, 'x4ze_c', subname)



  iall=-product(shape(y4za_c))*kind(y4za_c)
  deallocate(y4za_c, stat=istat)
  call memocc(istat, iall, 'y4za_c', subname)

  iall=-product(shape(y4zb_c))*kind(y4zb_c)
  deallocate(y4zb_c, stat=istat)
  call memocc(istat, iall, 'y4zb_c', subname)

  iall=-product(shape(y4zc_c))*kind(y4zc_c)
  deallocate(y4zc_c, stat=istat)
  call memocc(istat, iall, 'y4zc_c', subname)

  iall=-product(shape(y4ze_c))*kind(y4ze_c)
  deallocate(y4ze_c, stat=istat)
  call memocc(istat, iall, 'y4ze_c', subname)



  iall=-product(shape(x2y2aa_c))*kind(x2y2aa_c)
  deallocate(x2y2aa_c, stat=istat)
  call memocc(istat, iall, 'x2y2aa_c', subname)

  iall=-product(shape(x2y2ab_c))*kind(x2y2ab_c)
  deallocate(x2y2ab_c, stat=istat)
  call memocc(istat, iall, 'x2y2ab_c', subname)

  iall=-product(shape(x2y2ac_c))*kind(x2y2ac_c)
  deallocate(x2y2ac_c, stat=istat)
  call memocc(istat, iall, 'x2y2ac_c', subname)

  iall=-product(shape(x2y2ae_c))*kind(x2y2ae_c)
  deallocate(x2y2ae_c, stat=istat)
  call memocc(istat, iall, 'x2y2ae_c', subname)


  iall=-product(shape(x2y2ba_c))*kind(x2y2ba_c)
  deallocate(x2y2ba_c, stat=istat)
  call memocc(istat, iall, 'x2y2ba_c', subname)

  iall=-product(shape(x2y2bb_c))*kind(x2y2bb_c)
  deallocate(x2y2bb_c, stat=istat)
  call memocc(istat, iall, 'x2y2bb_c', subname)

  iall=-product(shape(x2y2bc_c))*kind(x2y2bc_c)
  deallocate(x2y2bc_c, stat=istat)
  call memocc(istat, iall, 'x2y2bc_c', subname)

  iall=-product(shape(x2y2be_c))*kind(x2y2be_c)
  deallocate(x2y2be_c, stat=istat)
  call memocc(istat, iall, 'x2y2be_c', subname)


  iall=-product(shape(x2y2ca_c))*kind(x2y2ca_c)
  deallocate(x2y2ca_c, stat=istat)
  call memocc(istat, iall, 'x2y2ca_c', subname)

  iall=-product(shape(x2y2cb_c))*kind(x2y2cb_c)
  deallocate(x2y2cb_c, stat=istat)
  call memocc(istat, iall, 'x2y2cb_c', subname)

  iall=-product(shape(x2y2cc_c))*kind(x2y2cc_c)
  deallocate(x2y2cc_c, stat=istat)
  call memocc(istat, iall, 'x2y2cc_c', subname)

  iall=-product(shape(x2y2ce_c))*kind(x2y2ce_c)
  deallocate(x2y2ce_c, stat=istat)
  call memocc(istat, iall, 'x2y2ce_c', subname)


  iall=-product(shape(x2y2ea_c))*kind(x2y2ea_c)
  deallocate(x2y2ea_c, stat=istat)
  call memocc(istat, iall, 'x2y2ea_c', subname)

  iall=-product(shape(x2y2eb_c))*kind(x2y2eb_c)
  deallocate(x2y2eb_c, stat=istat)
  call memocc(istat, iall, 'x2y2eb_c', subname)

  iall=-product(shape(x2y2ec_c))*kind(x2y2ec_c)
  deallocate(x2y2ec_c, stat=istat)
  call memocc(istat, iall, 'x2y2ec_c', subname)

  iall=-product(shape(x2y2ee_c))*kind(x2y2ee_c)
  deallocate(x2y2ee_c, stat=istat)
  call memocc(istat, iall, 'x2y2ee_c', subname)




  iall=-product(shape(x2y2aa_f))*kind(x2y2aa_f)
  deallocate(x2y2aa_f, stat=istat)
  call memocc(istat, iall, 'x2y2aa_f', subname)

  iall=-product(shape(x2y2ab_f))*kind(x2y2ab_f)
  deallocate(x2y2ab_f, stat=istat)
  call memocc(istat, iall, 'x2y2ab_f', subname)

  iall=-product(shape(x2y2ac_f))*kind(x2y2ac_f)
  deallocate(x2y2ac_f, stat=istat)
  call memocc(istat, iall, 'x2y2ac_f', subname)

  iall=-product(shape(x2y2ae_f))*kind(x2y2ae_f)
  deallocate(x2y2ae_f, stat=istat)
  call memocc(istat, iall, 'x2y2ae_f', subname)


  iall=-product(shape(x2y2ba_f))*kind(x2y2ba_f)
  deallocate(x2y2ba_f, stat=istat)
  call memocc(istat, iall, 'x2y2ba_f', subname)

  iall=-product(shape(x2y2bb_f))*kind(x2y2bb_f)
  deallocate(x2y2bb_f, stat=istat)
  call memocc(istat, iall, 'x2y2bb_f', subname)

  iall=-product(shape(x2y2bc_f))*kind(x2y2bc_f)
  deallocate(x2y2bc_f, stat=istat)
  call memocc(istat, iall, 'x2y2bc_f', subname)

  iall=-product(shape(x2y2be_f))*kind(x2y2be_f)
  deallocate(x2y2be_f, stat=istat)
  call memocc(istat, iall, 'x2y2be_f', subname)


  iall=-product(shape(x2y2ca_f))*kind(x2y2ca_f)
  deallocate(x2y2ca_f, stat=istat)
  call memocc(istat, iall, 'x2y2ca_f', subname)

  iall=-product(shape(x2y2cb_f))*kind(x2y2cb_f)
  deallocate(x2y2cb_f, stat=istat)
  call memocc(istat, iall, 'x2y2cb_f', subname)

  iall=-product(shape(x2y2cc_f))*kind(x2y2cc_f)
  deallocate(x2y2cc_f, stat=istat)
  call memocc(istat, iall, 'x2y2cc_f', subname)

  iall=-product(shape(x2y2ce_f))*kind(x2y2ce_f)
  deallocate(x2y2ce_f, stat=istat)
  call memocc(istat, iall, 'x2y2ce_f', subname)


  iall=-product(shape(x2y2ea_f))*kind(x2y2ea_f)
  deallocate(x2y2ea_f, stat=istat)
  call memocc(istat, iall, 'x2y2ea_f', subname)

  iall=-product(shape(x2y2eb_f))*kind(x2y2eb_f)
  deallocate(x2y2eb_f, stat=istat)
  call memocc(istat, iall, 'x2y2eb_f', subname)

  iall=-product(shape(x2y2ec_f))*kind(x2y2ec_f)
  deallocate(x2y2ec_f, stat=istat)
  call memocc(istat, iall, 'x2y2ec_f', subname)

  iall=-product(shape(x2y2ee_f))*kind(x2y2ee_f)
  deallocate(x2y2ee_f, stat=istat)
  call memocc(istat, iall, 'x2y2ee_f', subname)










  iall=-product(shape(x2ya_f))*kind(x2ya_f)
  deallocate(x2ya_f, stat=istat)
  call memocc(istat, iall, 'x2ya_f', subname)
  iall=-product(shape(x2yb_f))*kind(x2yb_f)
  deallocate(x2yb_f, stat=istat)
  call memocc(istat, iall, 'x2yb_f', subname)
  iall=-product(shape(x2yc_f))*kind(x2yc_f)
  deallocate(x2yc_f, stat=istat)
  call memocc(istat, iall, 'x2yc_f', subname)
  iall=-product(shape(x2ye_f))*kind(x2ye_f)
  deallocate(x2ye_f, stat=istat)
  call memocc(istat, iall, 'yze_f7', subname)

  iall=-product(shape(x2za_f))*kind(x2za_f)
  deallocate(x2za_f, stat=istat)
  call memocc(istat, iall, 'x2za_f', subname)
  iall=-product(shape(x2zb_f))*kind(x2zb_f)
  deallocate(x2zb_f, stat=istat)
  call memocc(istat, iall, 'x2zb_f', subname)
  iall=-product(shape(x2zc_f))*kind(x2zc_f)
  deallocate(x2zc_f, stat=istat)
  call memocc(istat, iall, 'x2zc_f', subname)
  iall=-product(shape(x2ze_f))*kind(x2ze_f)
  deallocate(x2ze_f, stat=istat)
  call memocc(istat, iall, 'zze_f7', subname)

  iall=-product(shape(y2za_f))*kind(y2za_f)
  deallocate(y2za_f, stat=istat)
  call memocc(istat, iall, 'y2za_f', subname)
  iall=-product(shape(y2zb_f))*kind(y2zb_f)
  deallocate(y2zb_f, stat=istat)
  call memocc(istat, iall, 'y2zb_f', subname)
  iall=-product(shape(y2zc_f))*kind(y2zc_f)
  deallocate(y2zc_f, stat=istat)
  call memocc(istat, iall, 'y2zc_f', subname)
  iall=-product(shape(y2ze_f))*kind(y2ze_f)
  deallocate(y2ze_f, stat=istat)
  call memocc(istat, iall, 'zze_f7', subname)



  iall=-product(shape(x2ya_f2))*kind(x2ya_f2)
  deallocate(x2ya_f2, stat=istat)
  call memocc(istat, iall, 'x2ya_f2', subname)
  iall=-product(shape(x2yb_f2))*kind(x2yb_f2)
  deallocate(x2yb_f2, stat=istat)
  call memocc(istat, iall, 'x2yb_f2', subname)
  iall=-product(shape(x2yc_f2))*kind(x2yc_f2)
  deallocate(x2yc_f2, stat=istat)
  call memocc(istat, iall, 'x2yc_f2', subname)
  iall=-product(shape(x2ye_f2))*kind(x2ye_f2)
  deallocate(x2ye_f2, stat=istat)
  call memocc(istat, iall, 'yze_f27', subname)





  iall=-product(shape(x4ya_f))*kind(x4ya_f)
  deallocate(x4ya_f, stat=istat)
  call memocc(istat, iall, 'x4ya_f', subname)
  iall=-product(shape(x4yb_f))*kind(x4yb_f)
  deallocate(x4yb_f, stat=istat)
  call memocc(istat, iall, 'x4yb_f', subname)
  iall=-product(shape(x4yc_f))*kind(x4yc_f)
  deallocate(x4yc_f, stat=istat)
  call memocc(istat, iall, 'x4yc_f', subname)
  iall=-product(shape(x4ye_f))*kind(x4ye_f)
  deallocate(x4ye_f, stat=istat)
  call memocc(istat, iall, 'yze_f7', subname)

  iall=-product(shape(x4za_f))*kind(x4za_f)
  deallocate(x4za_f, stat=istat)
  call memocc(istat, iall, 'x4za_f', subname)
  iall=-product(shape(x4zb_f))*kind(x4zb_f)
  deallocate(x4zb_f, stat=istat)
  call memocc(istat, iall, 'x4zb_f', subname)
  iall=-product(shape(x4zc_f))*kind(x4zc_f)
  deallocate(x4zc_f, stat=istat)
  call memocc(istat, iall, 'x4zc_f', subname)
  iall=-product(shape(x4ze_f))*kind(x4ze_f)
  deallocate(x4ze_f, stat=istat)
  call memocc(istat, iall, 'zze_f7', subname)

  iall=-product(shape(y4za_f))*kind(y4za_f)
  deallocate(y4za_f, stat=istat)
  call memocc(istat, iall, 'y4za_f', subname)
  iall=-product(shape(y4zb_f))*kind(y4zb_f)
  deallocate(y4zb_f, stat=istat)
  call memocc(istat, iall, 'y4zb_f', subname)
  iall=-product(shape(y4zc_f))*kind(y4zc_f)
  deallocate(y4zc_f, stat=istat)
  call memocc(istat, iall, 'y4zc_f', subname)
  iall=-product(shape(y4ze_f))*kind(y4ze_f)
  deallocate(y4ze_f, stat=istat)
  call memocc(istat, iall, 'zze_f7', subname)




END SUBROUTINE ConvolSextic













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


!!!$omp parallel default(private) &
!!!$omp shared(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3) &
!!!$omp shared(ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,w_c,w_f,y_c,y_f)& 
!!!$omp shared(w_f1,w_f2,w_f3,ad1_ext,bd1_ext,cd1_ext)

  ! x direction
  !!!$omp do  
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
  !!!$omp enddo
  

  ! y direction
  !!!$omp do
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
  !!!$omp enddo


  ! z direction
  !!!$omp do
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
  !!!$omp enddo


  

  ! wavelet part

  ! x direction
  !!!$omp do
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
  !!!$omp enddo


  ! y direction
  !!!$omp do
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
  !!!$omp enddo


  ! z direction
  !!!$omp do
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
  !!!$omp enddo

  !!!$omp end parallel


END SUBROUTINE createDerivativeBasis





