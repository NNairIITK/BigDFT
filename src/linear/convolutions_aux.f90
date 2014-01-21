!> @file
!! Convolutions for linear version (with quartic potentials)
!! @author
!!    Copyright (C) 2011-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
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
  real(wp), dimension(lowfil:lupfil) :: ed1_ext


! Copy the filters to the 'extended filters', i.e. add some zeros.
! This seems to be required since we use loop unrolling.
ad1_ext=0.d0
bd1_ext=0.d0
cd1_ext=0.d0
ed1_ext=0.d0

do i=lowfil,lupfil
    ad1_ext(i)=ad1(i)
    bd1_ext(i)=bd1(i)
    cd1_ext(i)=cd1(i)
    ed1_ext(i)=ed1(i) !this is only needed due to OpenMP, since it seems I cannot declare ed1 (which is in filterModule) as shared
end do


!$omp parallel default(none) &
!$omp shared(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3) &
!$omp shared(ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,w_c,w_f,y_c,y_f,x_c,x_f,z_c,z_f)& 
!$omp shared(w_f1,w_f2,w_f3,ad1_ext,bd1_ext,cd1_ext,ed1_ext)&
!$omp private(i,t,i1,i2,i3,icur,istart,iend,l)&
!$omp private(dyi,dyi0,dyi1,dyi2,dyi3,t112,t121,t122,t212,t221,t222,t211)

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

        if (iend-istart.ge.4) then
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

        if (iend-istart.ge.4) then
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

        if (iend-istart.ge.4) then
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
           istart=i3
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
              t221=t221 + w_f(2,i1+l,i2,i3)*cd1_ext(l) + w_f(3,i1+l,i2,i3)*ed1_ext(l)
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
              t122=t122 + w_f(4,i1,i2+l,i3)*cd1_ext(l) + w_f(6,i1,i2+l,i3)*ed1_ext(l)
              t212=t212 + w_f(5,i1,i2+l,i3)*ad1_ext(l) + w_f(7,i1,i2+l,i3)*bd1_ext(l)
              t221=t221 + w_f(1,i1,i2+l,i3)*cd1_ext(l) + w_f(3,i1,i2+l,i3)*ed1_ext(l)
              t222=t222 + w_f(5,i1,i2+l,i3)*cd1_ext(l) + w_f(7,i1,i2+l,i3)*ed1_ext(l)
              t121=t121 + w_f(2,i1,i2+l,i3)*ed1_ext(l)
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
              t122=t122 + w_f(2,i1,i2,i3+l)*cd1_ext(l) + w_f(6,i1,i2,i3+l)*ed1_ext(l)
              t212=t212 + w_f(1,i1,i2,i3+l)*cd1_ext(l) + w_f(5,i1,i2,i3+l)*ed1_ext(l)
              t221=t221 + w_f(3,i1,i2,i3+l)*ad1_ext(l) + w_f(7,i1,i2,i3+l)*bd1_ext(l)
              t222=t222 + w_f(3,i1,i2,i3+l)*cd1_ext(l) + w_f(7,i1,i2,i3+l)*ed1_ext(l)
              t112=t112 + w_f(4,i1,i2,i3+l)*ed1_ext(l)
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


subroutine createDerivativeX(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
     hgrid,ibyz_c,ibyz_f,w_c, w_f, w_f1, x_c, x_f)

  use module_base
  use filterModule
  implicit none

  ! Calling arguments
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  real(gp), intent(in) :: hgrid
  integer, dimension(2,0:n2,0:n3), intent(in) :: ibyz_c,ibyz_f
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: w_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in) :: w_f
  real(wp), dimension(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in) :: w_f1
  real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: x_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(out) :: x_f
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  integer :: i,t,i1,i2,i3
  integer :: icur,istart,iend,l
  real(wp) :: dyi,dyi0,dyi1,dyi2,dyi3,t112,t121,t122,t212,t221,t222,t211
  real(wp), dimension(-3+lowfil:lupfil+3) :: ad1_ext
  real(wp), dimension(-3+lowfil:lupfil+3) :: bd1_ext
  real(wp), dimension(-3+lowfil:lupfil+3) :: cd1_ext
  real(wp), dimension(lowfil:lupfil) :: ed1_ext


! Copy the filters to the 'extended filters', i.e. add some zeros.
! This seems to be required since we use loop unrolling.
ad1_ext=0.d0
bd1_ext=0.d0
cd1_ext=0.d0

do i=lowfil,lupfil
    ad1_ext(i)=ad1(i)
    bd1_ext(i)=bd1(i)
    cd1_ext(i)=cd1(i)
    ed1_ext(i)=ed1(i) !this is only needed due to OpenMP, since it seems I cannot declare ed1 (which is in filterModule) as shared
end do


!$omp parallel default(none) &
!$omp shared(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3) &
!$omp shared(ibyz_c,ibyz_f,w_c,w_f,x_c,x_f)& 
!$omp shared(w_f1,ad1_ext,bd1_ext,cd1_ext,ed1_ext)&
!$omp private(i,t,i1,i2,i3,icur,istart,iend,l)&
!$omp private(dyi,dyi0,dyi1,dyi2,dyi3,t112,t121,t122,t212,t221,t222,t211)

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

        if (iend-istart.ge.4) then
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

   ! wavelet part

  ! x direction
  !$omp do
  do i3=nfl3,nfu3
     do i2=nfl2,nfu2
        do i1=ibyz_f(1,i2,i3),ibyz_f(2,i2,i3)
           t112=0.0_wp;t121=0.0_wp;t122=0.0_wp;t212=0.0_wp;t221=0.0_wp;t222=0.0_wp;t211=0.0_wp
           do l=max(nfl1-i1,lowfil),min(lupfil,nfu1-i1)
              t121=t121 + w_f(2,i1+l,i2,i3)*ad1_ext(l) + w_f(3,i1+l,i2,i3)*bd1_ext(l)
              t221=t221 + w_f(2,i1+l,i2,i3)*cd1_ext(l) + w_f(3,i1+l,i2,i3)*ed1_ext(l)
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

!$omp end parallel

END SUBROUTINE createDerivativeX


subroutine createDerivativeY(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
     hgrid,ibxz_c,ibxz_f,w_c, w_f, w_f2, y_c, y_f)

  use module_base
  use filterModule
  implicit none

  ! Calling arguments
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  real(gp), intent(in) :: hgrid
  integer, dimension(2,0:n1,0:n3), intent(in) :: ibxz_c,ibxz_f
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: w_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in) :: w_f
  real(wp), dimension(nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), intent(in) :: w_f2
  real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: y_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(out) :: y_f
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  integer :: i,t,i1,i2,i3
  integer :: icur,istart,iend,l
  real(wp) :: dyi,dyi0,dyi1,dyi2,dyi3,t112,t121,t122,t212,t221,t222,t211
  real(wp), dimension(-3+lowfil:lupfil+3) :: ad1_ext
  real(wp), dimension(-3+lowfil:lupfil+3) :: bd1_ext
  real(wp), dimension(-3+lowfil:lupfil+3) :: cd1_ext
  real(wp), dimension(lowfil:lupfil) :: ed1_ext


! Copy the filters to the 'extended filters', i.e. add some zeros.
! This seems to be required since we use loop unrolling.
ad1_ext=0.d0
bd1_ext=0.d0
cd1_ext=0.d0

do i=lowfil,lupfil
    ad1_ext(i)=ad1(i)
    bd1_ext(i)=bd1(i)
    cd1_ext(i)=cd1(i)
    ed1_ext(i)=ed1(i) !this is only needed due to OpenMP, since it seems I cannot declare ed1 (which is in filterModule) as shared
end do


!$omp parallel default(none) &
!$omp shared(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3) &
!$omp shared(ibxz_c,ibxz_f,w_c,w_f,y_c,y_f)& 
!$omp shared(w_f2,ad1_ext,bd1_ext,cd1_ext,ed1_ext)&
!$omp private(i,t,i1,i2,i3,icur,istart,iend,l)&
!$omp private(dyi,dyi0,dyi1,dyi2,dyi3,t112,t121,t122,t212,t221,t222,t211)

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

        if (iend-istart.ge.4) then
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

  ! y direction
  !$omp do
  do i3=nfl3,nfu3
     do i1=nfl1,nfu1
        do i2=ibxz_f(1,i1,i3),ibxz_f(2,i1,i3)
           t112=0.0_wp;t121=0.0_wp;t122=0.0_wp;t212=0.0_wp;t221=0.0_wp;t222=0.0_wp;t211=0.0_wp
           do l=max(nfl2-i2,lowfil),min(lupfil,nfu2-i2)
              t112=t112 + w_f(4,i1,i2+l,i3)*ad1_ext(l) + w_f(6,i1,i2+l,i3)*bd1_ext(l)
              t211=t211 + w_f(1,i1,i2+l,i3)*ad1_ext(l) + w_f(3,i1,i2+l,i3)*bd1_ext(l)
              t122=t122 + w_f(4,i1,i2+l,i3)*cd1_ext(l) + w_f(6,i1,i2+l,i3)*ed1_ext(l)
              t212=t212 + w_f(5,i1,i2+l,i3)*ad1_ext(l) + w_f(7,i1,i2+l,i3)*bd1_ext(l)
              t221=t221 + w_f(1,i1,i2+l,i3)*cd1_ext(l) + w_f(3,i1,i2+l,i3)*ed1_ext(l)
              t222=t222 + w_f(5,i1,i2+l,i3)*cd1_ext(l) + w_f(7,i1,i2+l,i3)*ed1_ext(l)
              t121=t121 + w_f(2,i1,i2+l,i3)*ed1_ext(l)
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
!$omp end parallel

END SUBROUTINE createDerivativeY


subroutine createDerivativeZ(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
     hgrid,ibxy_c,ibxy_f,w_c, w_f, w_f3, z_c, z_f)

  use module_base
  use filterModule
  implicit none

  ! Calling arguments
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  real(gp), intent(in) :: hgrid
  integer, dimension(2,0:n1,0:n2), intent(in) :: ibxy_c,ibxy_f
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: w_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in) :: w_f
  real(wp), dimension(nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), intent(in) :: w_f3
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
  real(wp), dimension(lowfil:lupfil) :: ed1_ext


! Copy the filters to the 'extended filters', i.e. add some zeros.
! This seems to be required since we use loop unrolling.
ad1_ext=0.d0
bd1_ext=0.d0
cd1_ext=0.d0

do i=lowfil,lupfil
    ad1_ext(i)=ad1(i)
    bd1_ext(i)=bd1(i)
    cd1_ext(i)=cd1(i)
    ed1_ext(i)=ed1(i) !this is only needed due to OpenMP, since it seems I cannot declare ed1 (which is in filterModule) as shared
end do

!$omp parallel default(none) &
!$omp shared(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3) &
!$omp shared(ibxy_c,ibxy_f,w_c,w_f,z_c,z_f)& 
!$omp shared(w_f3,ad1_ext,bd1_ext,cd1_ext,ed1_ext)&
!$omp private(i,t,i1,i2,i3,icur,istart,iend,l)&
!$omp private(dyi,dyi0,dyi1,dyi2,dyi3,t112,t121,t122,t212,t221,t222,t211)

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

        if (iend-istart.ge.4) then
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
           istart=i3
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

  ! z direction
  !$omp do
  do i2=nfl2,nfu2
     do i1=nfl1,nfu1
        do i3=ibxy_f(1,i1,i2),ibxy_f(2,i1,i2)
           t112=0.0_wp;t121=0.0_wp;t122=0.0_wp;t212=0.0_wp;t221=0.0_wp;t222=0.0_wp;t211=0.0_wp
           do l=max(nfl3-i3,lowfil),min(lupfil,nfu3-i3)
              t121=t121 + w_f(2,i1,i2,i3+l)*ad1_ext(l) + w_f(6,i1,i2,i3+l)*bd1_ext(l)
              t211=t211 + w_f(1,i1,i2,i3+l)*ad1_ext(l) + w_f(5,i1,i2,i3+l)*bd1_ext(l)
              t122=t122 + w_f(2,i1,i2,i3+l)*cd1_ext(l) + w_f(6,i1,i2,i3+l)*ed1_ext(l)
              t212=t212 + w_f(1,i1,i2,i3+l)*cd1_ext(l) + w_f(5,i1,i2,i3+l)*ed1_ext(l)
              t221=t221 + w_f(3,i1,i2,i3+l)*ad1_ext(l) + w_f(7,i1,i2,i3+l)*bd1_ext(l)
              t222=t222 + w_f(3,i1,i2,i3+l)*cd1_ext(l) + w_f(7,i1,i2,i3+l)*ed1_ext(l)
              t112=t112 + w_f(4,i1,i2,i3+l)*ed1_ext(l)
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


END SUBROUTINE createDerivativeZ


!> Calculates the effective filter for the operator [kineticEnergy + (x-x0)^4].
subroutine getEffectiveFilterQuartic(parabPrefac,hgrid, x0, eff, filterCode)
use filterModule
implicit none

! Calling arguments
real(kind=8),intent(in) :: parabPrefac
real(kind=8),intent(in) :: hgrid                  !< Grid spacing
real(kind=8),intent(in) :: x0                     !< The center of the parabolic potential (x-x0)^2
real(kind=8),dimension(lb:ub),intent(out) :: eff  !< The effective filter
character(len=*), intent(in) :: filterCode

! Local variables
integer :: i
real(kind=8) :: fac, fac2, prefac1, hgrid2, hgrid3, x02, x03
prefac1=-.5d0/hgrid**2
fac=parabPrefac
fac2=parabPrefac*hgrid
hgrid2=hgrid**2
hgrid3=hgrid**3
x02=x0**2
x03=x0**3
! Determine which filter we have to calculate

if(filterCode=='a')then
    do i=lb,ub
        eff(i)=prefac1*a(i) + fac2*( hgrid3*a4(i) + 4*hgrid2*x0*a3(i) + 6*hgrid*x02*a2(i) + 4*x03*a1(i))
    end do
    eff(0)=eff(0)+fac*x0**4
elseif(filterCode=='b') then
    do i=lb,ub
        eff(i)=prefac1*b(i) + fac2*( hgrid3*b4(i) + 4*hgrid2*x0*b3(i) + 6*hgrid*x02*b2(i) + 4*x03*b1(i))
    end do
elseif(filterCode=='c') then
    do i=lb,ub
        eff(i)=prefac1*c(i) + fac2*( hgrid3*c4(i) + 4*hgrid2*x0*c3(i) + 6*hgrid*x02*c2(i) + 4*x03*c1(i))
    end do
elseif(filterCode=='e') then
    do i=lb,ub
        eff(i)=prefac1*e(i) + fac2*( hgrid3*e4(i) + 4*hgrid2*x0*e3(i) + 6*hgrid*x02*e2(i) + 4*x03*e1(i))
    end do
    eff(0)=eff(0)+fac*x0**4
else
    write(*,*) "ERROR: allowed values for 'filterCode' are 'a', 'b', 'c', 'e', whereas we found ", trim(filterCode)
    stop
end if

end subroutine getEffectiveFilterQuartic


!> Calculates the effective filter for the operator (x-x0)^4
pure subroutine getFilterQuartic(parabPrefac,hgrid, x0, eff, filterCode)

  use filterModule
  implicit none

  ! Calling arguments
  ! Calling arguments
  real(kind=8),intent(in) :: parabPrefac
  real(kind=8),intent(in) :: hgrid                  !< Grid spacing
  real(kind=8),intent(in) :: x0                     !< The center of the quartic potential (x-x0)^4
  real(kind=8),dimension(lb:ub),intent(out) :: eff  !< The effective filter
  character(len=*), intent(in) :: filterCode

  ! Local variables
  integer :: i
  real(kind=8) :: fac, fac2,  hgrid2, hgrid3, x02, x03
  real(kind=8) :: scale
  scale=1.d0
  !scale=1.d-1
  !scale=0.d-1
  !scale=5.d-2
  !fac=dble(max(100-int(dble(it)/2.d0),1))*parabPrefac
  !fac2=dble(max(100-int(dble(it)/2.d0),1))*parabPrefac*hgrid
  fac=parabPrefac*scale
  fac2=parabPrefac*hgrid*scale
  hgrid2=hgrid**2
  hgrid3=hgrid**3
  x02=x0**2
  x03=x0**3
  ! Determine which filter we have to calculate

  select case(filtercode)
  case('a')
     !if(filterCode=='a') then
     do i=lb,ub
        !eff(i)=prefac1*a(i) + fac2*(hgrid*a2(i)+2*x0*a1(i))
        eff(i) = fac2*( hgrid3*a4(i) + 4*hgrid2*x0*a3(i) + 6*hgrid*x02*a2(i) + 4*x03*a1(i))
     end do
     !eff(0)=eff(0)+fac*x0**2
     eff(0)=eff(0)+fac*x0**4
  case('b')
  !elseif(filterCode=='b') then
     do i=lb,ub
        !eff(i)=prefac1*b(i) + fac2*(hgrid*b2(i)+2*x0*b1(i))
        eff(i) = fac2*( hgrid3*b4(i) + 4*hgrid2*x0*b3(i) + 6*hgrid*x02*b2(i) + 4*x03*b1(i))
     end do
  case('c')
  !elseif(filterCode=='c') then
     do i=lb,ub
        !eff(i)=prefac1*c(i) + fac2*(hgrid*c2(i)+2*x0*c1(i))
        eff(i) = fac2*( hgrid3*c4(i) + 4*hgrid2*x0*c3(i) + 6*hgrid*x02*c2(i) + 4*x03*c1(i))
     end do
  case('e')
  !elseif(filterCode=='e') then
     do i=lb,ub
        !eff(i)=prefac1*e(i) + fac2*(hgrid*e2(i)+2*x0*e1(i))
        eff(i) = fac2*( hgrid3*e4(i) + 4*hgrid2*x0*e3(i) + 6*hgrid*x02*e2(i) + 4*x03*e1(i))
     end do
     !eff(0)=eff(0)+fac*x0**2
     eff(0)=eff(0)+fac*x0**4
  !else
  !   write(*,*) "ERROR: allowed values for 'filterCode' are 'a', 'b', 'c', 'e', whereas we found ", trim(filterCode)
  !   stop
  !end if
  end select

end subroutine getFilterQuartic


!> Calculates the effective filter for the operator (x-x0)^2
subroutine getFilterQuadratic(parabPrefac,hgrid, x0, eff, filterCode)
use filterModule
implicit none

! Calling arguments
real(kind=8),intent(in) :: parabPrefac
real(kind=8),intent(in) :: hgrid                  !< Grid spacing
real(kind=8),intent(in) :: x0                     !< The center of the parabolic potential (x-x0)^2
real(kind=8),dimension(lb:ub),intent(out) :: eff  !< The effective filter
character(len=*), intent(in) :: filterCode


! Local variables
integer :: i
real(kind=8) :: fac, fac2, hgrid2, hgrid3, x02, x03
real(kind=8) :: scale
scale=1.d0
!scale=1.d-1
!scale=0.d-1
!scale=5.d-2
!fac=dble(max(100-int(dble(it)/2.d0),1))*parabPrefac
!fac2=dble(max(100-int(dble(it)/2.d0),1))*parabPrefac*hgrid
fac=parabPrefac*scale
fac2=parabPrefac*hgrid*scale
hgrid2=hgrid**2
hgrid3=hgrid**3
x02=x0**2
x03=x0**3
! Determine which filter we have to calculate

if(filterCode=='a') then
    do i=lb,ub
        eff(i) = fac2*( hgrid*a2(i) + 2.d0*x0*a1(i) )
        !eff(i) = fac2*( hgrid3*a4(i) + 4*hgrid2*x0*a3(i) + 6*hgrid*x02*a2(i) + 4*x03*a1(i))
    end do
    eff(0)=eff(0)+fac*x0**2
    !eff(0)=eff(0)+fac*x0**4
elseif(filterCode=='b') then
    do i=lb,ub
        eff(i) = fac2*( hgrid*b2(i) + 2.d0*x0*b1(i) )
        !eff(i) = fac2*( hgrid3*b4(i) + 4*hgrid2*x0*b3(i) + 6*hgrid*x02*b2(i) + 4*x03*b1(i))
    end do
elseif(filterCode=='c') then
    do i=lb,ub
        eff(i) = fac2*( hgrid*c2(i) + 2.d0*x0*c1(i) )
        !eff(i) = fac2*( hgrid3*c4(i) + 4*hgrid2*x0*c3(i) + 6*hgrid*x02*c2(i) + 4*x03*c1(i))
    end do
elseif(filterCode=='e') then
    do i=lb,ub
        eff(i) = fac2*( hgrid*e2(i) + 2.d0*x0*e1(i) )
        !eff(i) = fac2*( hgrid3*e4(i) + 4*hgrid2*x0*e3(i) + 6*hgrid*x02*e2(i) + 4*x03*e1(i))
    end do
    eff(0)=eff(0)+fac*x0**2
    !eff(0)=eff(0)+fac*x0**4
else
    write(*,*) "ERROR: allowed values for 'filterCode' are 'a', 'b', 'c', 'e', whereas we found ", trim(filterCode)
    stop
end if

end subroutine getFilterQuadratic


!> Expands the compressed wavefunction in vector form (psi_c,psi_f) into the psig format
subroutine uncompress_for_quartic_convolutions(n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3,  & 
     mseg_c, mvctr_c, keyg_c, keyv_c,  & 
     mseg_f, mvctr_f, keyg_f, keyv_f,  & 
     scal, psi_c, psi_f, &
     work)
  use module_base
  use module_types
  implicit none
  integer,intent(in) :: n1, n2, n3, nfl1, nfu1, nfl2, nfu2, nfl3, nfu3, mseg_c, mvctr_c, mseg_f, mvctr_f
  integer,dimension(mseg_c),intent(in) :: keyv_c
  integer,dimension(mseg_f),intent(in) :: keyv_f
  integer,dimension(2,mseg_c),intent(in) :: keyg_c
  integer,dimension(2,mseg_f),intent(in) :: keyg_f
  real(wp),dimension(0:3),intent(in) :: scal
  real(wp),dimension(mvctr_c),intent(in) :: psi_c
  real(wp),dimension(7,mvctr_f),intent(in) :: psi_f
  type(workarrays_quartic_convolutions),intent(inout) :: work
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

  !$omp parallel default(private) &
  !$omp shared(scal) &
  !$omp shared(psi_c,psi_f,keyv_c,keyg_c,keyv_f,keyg_f,n1,n2,n3,mseg_c,mseg_f,work)
  !!! coarse part
  !$omp do
  do iseg=1,mseg_c
     jj=keyv_c(iseg)
     j0=keyg_c(1,iseg)
     j1=keyg_c(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        work%xx_c(i,i2,i3)=psi_c(i-i0+jj)*scal(0)
        work%xy_c(i2,i,i3)=psi_c(i-i0+jj)*scal(0)
        work%xz_c(i3,i,i2)=psi_c(i-i0+jj)*scal(0)
     enddo
  enddo
  !$omp enddo
  !!! fine part
  !$omp do
  do iseg=1,mseg_f
     jj=keyv_f(iseg)
     j0=keyg_f(1,iseg)
     j1=keyg_f(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        work%xx_f1(i,i2,i3)=psi_f(1,i-i0+jj)*scal(1)
        work%xx_f(1,i,i2,i3)=psi_f(1,i-i0+jj)*scal(1)
        work%xy_f(1,i2,i,i3)=psi_f(1,i-i0+jj)*scal(1)
        work%xz_f(1,i3,i,i2)=psi_f(1,i-i0+jj)*scal(1)

        work%xy_f2(i2,i,i3)=psi_f(2,i-i0+jj)*scal(1)
        work%xx_f(2,i,i2,i3)=psi_f(2,i-i0+jj)*scal(1)
        work%xy_f(2,i2,i,i3)=psi_f(2,i-i0+jj)*scal(1)
        work%xz_f(2,i3,i,i2)=psi_f(2,i-i0+jj)*scal(1)

        work%xx_f(3,i,i2,i3)=psi_f(3,i-i0+jj)*scal(2)
        work%xy_f(3,i2,i,i3)=psi_f(3,i-i0+jj)*scal(2)
        work%xz_f(3,i3,i,i2)=psi_f(3,i-i0+jj)*scal(2)

        work%xz_f4(i3,i,i2)=psi_f(4,i-i0+jj)*scal(1)
        work%xx_f(4,i,i2,i3)=psi_f(4,i-i0+jj)*scal(1)
        work%xy_f(4,i2,i,i3)=psi_f(4,i-i0+jj)*scal(1)
        work%xz_f(4,i3,i,i2)=psi_f(4,i-i0+jj)*scal(1)

        work%xx_f(5,i,i2,i3)=psi_f(5,i-i0+jj)*scal(2)
        work%xy_f(5,i2,i,i3)=psi_f(5,i-i0+jj)*scal(2)
        work%xz_f(5,i3,i,i2)=psi_f(5,i-i0+jj)*scal(2)

        work%xx_f(6,i,i2,i3)=psi_f(6,i-i0+jj)*scal(2)
        work%xy_f(6,i2,i,i3)=psi_f(6,i-i0+jj)*scal(2)
        work%xz_f(6,i3,i,i2)=psi_f(6,i-i0+jj)*scal(2)

        work%xx_f(7,i,i2,i3)=psi_f(7,i-i0+jj)*scal(3)
        work%xy_f(7,i2,i,i3)=psi_f(7,i-i0+jj)*scal(3)
        work%xz_f(7,i3,i,i2)=psi_f(7,i-i0+jj)*scal(3)
     enddo
  enddo
 !$omp enddo
 !$omp end parallel
 

END SUBROUTINE uncompress_for_quartic_convolutions
