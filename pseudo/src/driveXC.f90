!> @file
!! Driver of the XC potential for pseudo program
!! @author
!!    Alex Willand, under the supervision of Stefan Goedecker
!!    gpu accelerated routines by Raffael Widmer
!!    parts of this program were based on the fitting program by matthias krack
!!    http://cvs.berlios.de/cgi-bin/viewcvs.cgi/cp2k/potentials/goedecker/pseudo/v2.2/


!> Integrates LDA or GGA energy and potential terms.
!! Derived from the routine ggaenergy_15.f.
!! This version uses libXC via the bindings in xcfunction.f90.
!! The routine is split into four cases for each LDA and GGA
!! in closed shell or spin polarized calculations.
!! In contrast, the original ggaenergy_15.f always did the entire
!! GGA procedure of calculating gradients and their corrections.
subroutine driveXC(nspol,nrad,r,rw,rd,rho,enexc,pot,eps)
   use libxcModule
   implicit none
   !Arguments
   integer, intent(in) :: nspol                             !< Number of spin components
   integer, intent(in) :: nrad                              !< Number of grid points
   real(kind=8), dimension(nrad), intent(in) :: r,rw,rd     !< Grid points (normal, strange1 and strange2)
   real(kind=8), dimension(nrad,nspol), intent(in) :: rho   !< Electronic density
   real(kind=8), dimension(nrad,nspol), intent(out) :: pot  !< XC potential
   real(kind=8), dimension(nrad), intent(out) :: eps        !< XC energy density
   real(kind=8), intent(out) :: enexc                       !< XC energy

   if (libxc_functionals_isgga()) then
      select case(nspol)
         case (1)
            call driveGGAsimple(nspol,nrad,rw,rd,rho,enexc,pot,eps)
         case (2)
            call driveGGApolarized(nspol,nrad,rw,rd,rho,enexc,pot,eps)
      end select
   else
      select case(nspol)
         case (1)
            call driveLDAsimple(nspol,nrad,r,rho,enexc,pot,eps)
         case (2)
            call driveLDApolarized(nspol,nrad,rw,rho,enexc,pot,eps)
      end select
   end if

end subroutine driveXC


!> Greatly simplified version of ggaenergy_15 for LDA XC
subroutine driveLDAsimple(nspol,nrad,r,rho,enexc,pot,eps)
   use libxcModule
   implicit real*8 (a-h,o-z)
   !Arguments
   integer, intent(in) :: nspol                             !< Number of spin components
   integer, intent(in) :: nrad                              !< Number of grid points
   real(kind=8), dimension(nrad), intent(in) :: r           !< Grid points
   real(kind=8), dimension(nrad,nspol), intent(in) :: rho   !< Electronic density
   real(kind=8), dimension(nrad,nspol), intent(out) :: pot  !< XC potential
   real(kind=8), dimension(nrad), intent(out) :: eps        !< XC energy density
   real(kind=8), intent(out) :: enexc                       !< XC energy
   !Local variables
   real(kind=8), dimension(1) :: Vxc,dEdg, grad,rhoj        !< Arguments to XCfunction
   real(kind=8) :: Exc                                      !< Arguments to XCfunction

   if (nspol /= 1) then
      write(*,*) 'nspol=',nspol
      stop 'driveLDAsimple: nspol should be equal to 1!'
   end if
   enexc=0.d0
   pot=0d0
   eps=0d0
   grad=0d0

   do j=1,nrad
      rhoj=rho(j,:)
      call xcfunction(nspol,rhoj,grad,Exc,Vxc,dEdg)
      enexc=enexc+Exc*r(j)*rhoj(1)
      eps(j)=eps(j)+Exc
      pot(j,1)=pot(j,1)+Vxc(1)*r(j)
   enddo


   do j=2,nrad
      pot(j,1)=pot(j,1)/r(j)
   enddo
   pot(1,1)=pot(2,1)

end subroutine driveLDAsimple


!> The simple LDA XC driver generalized for two spin channels
!! @warning This routine is called with nspol == 2 !
subroutine driveLDApolarized(nspol,nrad,rw,rho,enexc,pot,eps)
   use libxcModule
   implicit real*8 (a-h,o-z)
   !Arguments
   integer, intent(in) :: nspol                             !< Number of spin components
   integer, intent(in) :: nrad                              !< Number of grid points
   real(kind=8), dimension(nrad), intent(in) :: rw          !< Integration grid points
   real(kind=8), dimension(nrad,nspol), intent(in) :: rho   !< Electronic density
   real(kind=8), dimension(nrad,nspol), intent(out) :: pot  !< XC potential
   real(kind=8), dimension(nrad), intent(out) :: eps        !< XC energy density
   real(kind=8), intent(out) :: enexc                       !< XC energy
   !Local variables
   real(kind=8), dimension(2) :: Vxc,dEdg,grad,rhoj         !< Arguments to XCfunction
   real(kind=8) :: Exc                                      !< Arguments to XCfunction

   if (nspol /= 2) then
      write(*,*) 'nspol=',nspol
      stop 'driveLDApolarized: nspol should be equal to 2!'
   end if
   enexc=0.d0
   pot=0d0
   eps=0d0
   grad=0d0

   do j=1,nrad
      rhoj = rho(j,:)
      call xcfunction(nspol,rhoj,grad,Exc,Vxc,dEdg)
      enexc=enexc+Exc*rw(j)*(rhoj(1)+rhoj(2))
      eps(j)=eps(j)+Exc
      pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
      pot(j,2)=pot(j,2)+Vxc(2)*rw(j)
   enddo

   do j=2,nrad
      pot(j,1)=pot(j,1)/rw(j)
      pot(j,2)=pot(j,2)/rw(j)
   enddo
   pot(1,:)=pot(2,:)

end subroutine driveLDApolarized


!> This one is very close to the original version of ggaenergy_15 
subroutine driveGGAsimple(nspol,nrad,rw,rd,rho,enexc,pot,eps)
   use libxcModule
   implicit real*8 (a-h,o-z)
   !Arguments
   integer, intent(in) :: nspol                             !< Number of spin components
   integer, intent(in) :: nrad                              !< Number of grid points
   real(kind=8), dimension(nrad), intent(in) :: rw,rd       !< Grid points (up, down)
   real(kind=8), dimension(nrad,nspol), intent(in) :: rho   !< Electronic density
   real(kind=8), dimension(nrad,nspol), intent(out) :: pot  !< XC potential
   real(kind=8), dimension(nrad), intent(out) :: eps        !< XC energy density
   real(kind=8), intent(out) :: enexc                       !< XC energy
   !Local variables
   real(kind=8), dimension(1) :: Vxc,dEdg, grad,rhoj        !< Arguments to XCfunction
   real(kind=8) :: Exc                                      !< Arguments to XCfunction
   real(kind=8) :: c(-8:8), sign(1)                         !< Variables for derivatives 

   if (nspol /= 1) then
      write(*,*) 'nspol=',nspol
      stop 'driveGGAsimple: nspol should be equal to 1!'
   end if
   enexc=0.d0
   pot=0d0
   eps=0d0
   grad=0d0

   j=1

   c(0)=-2.717857142857143d0
   c(1)=8.d0
   c(2)=-14.d0
   c(3)=18.66666666666667d0
   c(4)=-17.5d0
   c(5)=11.2d0
   c(6)=-4.666666666666667d0
   c(7)=1.142857142857143d0
   c(8)=-0.125d0

   grad=0.d0
   do i=-0,8
      grad(1)=grad(1)+c(i)*rho(j+i,1)
   enddo
   if (grad(1).ge.0.d0) then 
      sign(1)=rd(j)
   else
      sign(1)=-rd(j)
   endif
   grad(1)=sign(1)*grad(1)
   rhoj=rho(j,:)
   call xcfunction(nspol,rhoj,grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*rho(j,1)
   eps(j)=eps(j)+Exc
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   do i=-0,8
      pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   enddo

   j=2
    c(-1)=-0.1111111111111111d0
    c(0)=-1.717857142857143d0
    c(1)=4.d0
    c(2)=-4.666666666666667d0
    c(3)=4.666666666666667d0
    c(4)=-3.5d0
    c(5)=1.866666666666667d0
    c(6)=-0.6666666666666666d0
    c(7)=0.1428571428571428d0
    c(8)=-0.01388888888888889d0
   grad=0.d0
   do i=-1,8
      grad(1)=grad(1)+c(i)*rho(j+i,1)
   enddo
   if (grad(1).ge.0.d0) then 
      sign(1)=rd(j)
   else
      sign(1)=-rd(j)
   endif
   grad(1)=sign(1)*grad(1)
   rhoj=rho(j,:)
   call xcfunction(nspol,rhoj,grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*rho(j,1)
   eps(j)=eps(j)+Exc
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   do i=-1,8
      pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   enddo

   j=3
    c(-2)=0.01111111111111111d0
    c(-1)=-0.2222222222222222d0
    c(0)=-1.217857142857143d0
    c(1)=2.666666666666666d0
    c(2)=-2.333333333333333d0
    c(3)=1.866666666666667d0
    c(4)=-1.166666666666667d0
    c(5)=0.5333333333333333d0
    c(6)=-0.1666666666666666d0
    c(7)=0.03174603174603174d0
    c(8)=-0.2777777777777778d-2
   grad=0.d0
   do i=-2,8
      grad(1)=grad(1)+c(i)*rho(j+i,1)
   enddo
   if (grad(1).ge.0.d0) then 
      sign(1)=rd(j)
   else
      sign(1)=-rd(j)
   endif
   grad(1)=sign(1)*grad(1)
   rhoj=rho(j,:)
   call xcfunction(nspol,rhoj,grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*rho(j,1)
   eps(j)=eps(j)+Exc
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   do i=-2,8
      pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   enddo

   j=4
    c(-3)=-0.202020202020202d-2
    c(-2)=0.03333333333333333d0
    c(-1)=-0.3333333333333333d0
    c(0)=-0.88452380952381d0
    c(1)=2.d0
    c(2)=-1.4d0
    c(3)=0.933333333333333d0
    c(4)=-0.5d0
    c(5)=0.2d0
    c(6)=-0.05555555555555556d0
    c(7)=0.952380952380952d-2
    c(8)=-0.7575757575757577d-3
   grad=0.d0
   do i=-3,8
      grad(1)=grad(1)+c(i)*rho(j+i,1)
   enddo
   if (grad(1).ge.0.d0) then 
      sign(1)=rd(j)
   else
      sign(1)=-rd(j)
   endif
   grad(1)=sign(1)*grad(1)
   rhoj=rho(j,:)
   call xcfunction(nspol,rhoj,grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*rho(j,1)
   eps(j)=eps(j)+Exc
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   do i=-3,8
      pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   enddo

   j=5
    c(-4)=0.5050505050505051d-3
    c(-3)=-0.808080808080808d-2
    c(-2)=0.06666666666666666d0
    c(-1)=-0.4444444444444445d0
    c(0)=-0.6345238095238095d0
    c(1)=1.6d0
    c(2)=-0.933333333333333d0
    c(3)=0.5333333333333333d0
    c(4)=-0.25d0
    c(5)=0.0888888888888889d0
    c(6)=-0.02222222222222222d0
    c(7)=0.3463203463203463d-2
    c(8)=-0.2525252525252525d-3
   grad=0.d0
   do i=-4,8
      grad(1)=grad(1)+c(i)*rho(j+i,1)
   enddo
   if (grad(1).ge.0.d0) then 
      sign(1)=rd(j)
   else
      sign(1)=-rd(j)
   endif
   grad(1)=sign(1)*grad(1)
   rhoj=rho(j,:)
   call xcfunction(nspol,rhoj,grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*rho(j,1)
   eps(j)=eps(j)+Exc
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   do i=-4,8
      pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   enddo

   j=6
    c(-5)=-0.1554001554001554d-3
    c(-4)=0.2525252525252525d-2
    c(-3)=-0.0202020202020202d0
    c(-2)=0.1111111111111111d0
    c(-1)=-0.5555555555555556d0
    c(0)=-0.4345238095238095d0
    c(1)=1.333333333333333d0
    c(2)=-0.6666666666666666d0
    c(3)=0.3333333333333333d0
    c(4)=-0.1388888888888889d0
    c(5)=0.04444444444444445d0
    c(6)=-0.0101010101010101d0
    c(7)=0.1443001443001443d-2
    c(8)=-0.971250971250971d-4
   grad=0.d0
   do i=-5,8
      grad(1)=grad(1)+c(i)*rho(j+i,1)
   enddo
   if (grad(1).ge.0.d0) then 
      sign(1)=rd(j)
   else
      sign(1)=-rd(j)
   endif
   grad(1)=sign(1)*grad(1)
   rhoj=rho(j,:)
   call xcfunction(nspol,rhoj,grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*rho(j,1)
   eps(j)=eps(j)+Exc
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   do i=-5,8
      pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   enddo

   j=7
    c(-6)=0.555000555000555d-4
    c(-5)=-0.932400932400932d-3
    c(-4)=0.7575757575757577d-2
    c(-3)=-0.04040404040404041d0
    c(-2)=0.1666666666666666d0
    c(-1)=-0.6666666666666666d0
    c(0)=-0.2678571428571428d0
    c(1)=1.142857142857143d0
    c(2)=-0.5d0
    c(3)=0.2222222222222222d0
    c(4)=-0.0833333333333333d0
    c(5)=0.02424242424242424d0
    c(6)=-0.5050505050505051d-2
    c(7)=0.6660006660006659d-3
    c(8)=-0.4162504162504162d-4
   grad=0.d0
   do i=-6,8
      grad(1)=grad(1)+c(i)*rho(j+i,1)
   enddo
   if (grad(1).ge.0.d0) then 
      sign(1)=rd(j)
   else
      sign(1)=-rd(j)
   endif
   grad(1)=sign(1)*grad(1)
   rhoj=rho(j,:)
   call xcfunction(nspol,rhoj,grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*rho(j,1)
   eps(j)=eps(j)+Exc
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   do i=-6,8
      pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   enddo

   j=8
    c(-7)=-0.222000222000222d-4
    c(-6)=0.3885003885003884d-3
    c(-5)=-0.3263403263403263d-2
    c(-4)=0.01767676767676768d0
    c(-3)=-0.07070707070707071d0
    c(-2)=0.2333333333333333d0
    c(-1)=-0.7777777777777778d0
    c(0)=-0.125d0
    c(1)=1.d0
    c(2)=-0.3888888888888889d0
    c(3)=0.1555555555555556d0
    c(4)=-0.05303030303030303d0
    c(5)=0.01414141414141414d0
    c(6)=-0.2719502719502719d-2
    c(7)=0.3330003330003329d-3
    c(8)=-0.1942501942501942d-4
   grad=0.d0
   do i=-7,8
      grad(1)=grad(1)+c(i)*rho(j+i,1)
   enddo
   if (grad(1).ge.0.d0) then 
      sign(1)=rd(j)
   else
      sign(1)=-rd(j)
   endif
   grad(1)=sign(1)*grad(1)
   rhoj=rho(j,:)
   call xcfunction(nspol,rhoj,grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*rho(j,1)
   eps(j)=eps(j)+Exc
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   do i=-7,8
      pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   enddo


    c(-8)=9.71250971250971d-6
    c(-7)=-0.1776001776001776d-3
    c(-6)=0.1554001554001554d-2
    c(-5)=-0.87024087024087d-2
    c(-4)=0.3535353535353535d-1
    c(-3)=-0.1131313131313131d0
    c(-2)=0.3111111111111111d0
    c(-1)=-0.888888888888889d0
    c(0)=0.d0
    c(1)=0.888888888888889d0
    c(2)=-0.3111111111111111d0
    c(3)=0.1131313131313131d0
    c(4)=-0.3535353535353535d-1
    c(5)=0.87024087024087d-2
    c(6)=-0.1554001554001554d-2
    c(7)=0.1776001776001776d-3
    c(8)=-9.71250971250971d-6
   do 100,j=9,nrad-8
   grad=0.d0
   do i=-8,8
      grad(1)=grad(1)+c(i)*rho(j+i,1)
   enddo
   if (grad(1).ge.0.d0) then 
      sign(1)=rd(j)
   else
      sign(1)=-rd(j)
   endif
   grad(1)=sign(1)*grad(1)
   rhoj=rho(j,:)
   call xcfunction(nspol,rhoj,grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*rho(j,1)
   eps(j)=eps(j)+Exc
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   do i=-8,8
      pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   enddo
100   continue

   j=nrad-7
   grad=0.d0
   do i=-8,7
      grad(1)=grad(1)+c(i)*rho(j+i,1)
   enddo
   if (grad(1).ge.0.d0) then 
      sign(1)=rd(j)
   else
      sign(1)=-rd(j)
   endif
   grad(1)=sign(1)*grad(1)
   rhoj=rho(j,:)
   call xcfunction(nspol,rhoj,grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*rho(j,1)
   eps(j)=eps(j)+Exc
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   do i=-8,7
      pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   enddo

   j=nrad-6
   grad=0.d0
   do i=-8,6
      grad(1)=grad(1)+c(i)*rho(j+i,1)
   enddo
   if (grad(1).ge.0.d0) then 
      sign(1)=rd(j)
   else
      sign(1)=-rd(j)
   endif
   grad(1)=sign(1)*grad(1)
   rhoj=rho(j,:)
   call xcfunction(nspol,rhoj,grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*rho(j,1)
   eps(j)=eps(j)+Exc
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   do i=-8,6
      pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   enddo

   j=nrad-5
   grad=0.d0
   do i=-8,5
      grad(1)=grad(1)+c(i)*rho(j+i,1)
   enddo
   if (grad(1).ge.0.d0) then 
      sign(1)=rd(j)
   else
      sign(1)=-rd(j)
   endif
   grad(1)=sign(1)*grad(1)
   rhoj=rho(j,:)
   call xcfunction(nspol,rhoj,grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*rho(j,1)
   eps(j)=eps(j)+Exc
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   do i=-8,5
      pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   enddo

   j=nrad-4
   grad=0.d0
   do i=-8,4
      grad(1)=grad(1)+c(i)*rho(j+i,1)
   enddo
   if (grad(1).ge.0.d0) then 
      sign(1)=rd(j)
   else
      sign(1)=-rd(j)
   endif
   grad(1)=sign(1)*grad(1)
   rhoj=rho(j,:)
   call xcfunction(nspol,rhoj,grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*rho(j,1)
   eps(j)=eps(j)+Exc
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   do i=-8,4
      pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   enddo

   j=nrad-3
   grad=0.d0
   do i=-8,3
      grad(1)=grad(1)+c(i)*rho(j+i,1)
   enddo
   if (grad(1).ge.0.d0) then 
      sign(1)=rd(j)
   else
      sign(1)=-rd(j)
   endif
   grad(1)=sign(1)*grad(1)
   rhoj=rho(j,:)
   call xcfunction(nspol,rhoj,grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*rho(j,1)
   eps(j)=eps(j)+Exc
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   do i=-8,3
      pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   enddo

   j=nrad-2
   grad=0.d0
   do i=-8,2
      grad(1)=grad(1)+c(i)*rho(j+i,1)
   enddo
   if (grad(1).ge.0.d0) then 
      sign(1)=rd(j)
   else
      sign(1)=-rd(j)
   endif
   grad(1)=sign(1)*grad(1)
   rhoj=rho(j,:)
   call xcfunction(nspol,rhoj,grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*rho(j,1)
   eps(j)=eps(j)+Exc
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   do i=-8,2
      pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   enddo

   j=nrad-1
   grad=0.d0
   do i=-8,1
      grad(1)=grad(1)+c(i)*rho(j+i,1)
   enddo
   if (grad(1).ge.0.d0) then 
      sign(1)=rd(j)
   else
      sign(1)=-rd(j)
   endif
   grad(1)=sign(1)*grad(1)
   rhoj=rho(j,:)
   call xcfunction(nspol,rhoj,grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*rho(j,1)
   eps(j)=eps(j)+Exc
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   do i=-8,1
      pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   enddo

   j=nrad-0
   grad=0.d0
   do i=-8,0
      grad(1)=grad(1)+c(i)*rho(j+i,1)
   enddo
   if (grad(1).ge.0.d0) then 
      sign(1)=rd(j)
   else
      sign(1)=-rd(j)
   endif
   grad(1)=sign(1)*grad(1)
   rhoj=rho(j,:)
   call xcfunction(nspol,rhoj,grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*rho(j,1)
   eps(j)=eps(j)+Exc
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   do i=-8,0
      pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   enddo

   do j=2,nrad
      pot(j,1)=pot(j,1)/rw(j)
   enddo
   pot(1,1)=pot(2,1)

end subroutine driveGGAsimple


!> ggaenergy_15 generalized to colinear spin polarization
subroutine driveGGApolarized(nspol,nrad,rw,rd,rho,enexc,pot,eps)
   use libxcModule
   implicit real*8 (a-h,o-z)
   !Arguments
   integer, intent(in) :: nspol                             !< Number of spin components
   integer, intent(in) :: nrad                              !< Number of grid points
   real(kind=8), dimension(nrad), intent(in) :: rw,rd       !< Grid points (strange1 and strange2)
   real(kind=8), dimension(nrad,nspol), intent(in) :: rho   !< Electronic density
   real(kind=8), dimension(nrad,nspol), intent(out) :: pot  !< XC potential
   real(kind=8), dimension(nrad), intent(out) :: eps        !< XC energy density
   real(kind=8), intent(out) :: enexc                       !< XC energy
   !Local variables
   real(kind=8), dimension(2) :: Vxc,dEdg, grad,rhoj        !< Arguments to XCfunction
   real(kind=8) :: Exc                                      !< Arguments to XCfunction
   real(kind=8) :: c(-8:8), sign(2)                         !< Variables for derivatives 

   if (nspol /= 2) then
      write(*,*) 'nspol=',nspol
      stop 'driveGGApolarized: nspol should be equal to 2!'
   end if
   enexc=0.d0
   pot=0d0
   eps=0d0
   grad=0d0

   j=1
    c(0)=-2.717857142857143d0
    c(1)=8.d0
    c(2)=-14.d0
    c(3)=18.66666666666667d0
    c(4)=-17.5d0
    c(5)=11.2d0
    c(6)=-4.666666666666667d0
    c(7)=1.142857142857143d0
    c(8)=-0.125d0
   grad=0.d0
   do i=-0,8
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   grad(2)=grad(2)+c(i)*rho(j+i,2)
   enddo
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   if (grad(2).ge.0.d0) then 
   sign(2)=rd(j)
   else
   sign(2)=-rd(j)
   endif
   grad=sign*grad
   rhoj=rho(j,:)
   call xcfunction(nspol,rhoj,grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*(rho(j,1)+rho(j,2))
   eps(j)=eps(j)+Exc
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   pot(j,2)=pot(j,2)+Vxc(2)*rw(j)
   do i=-0,8
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   pot(j+i,2)=pot(j+i,2)+(sign(2)*c(i)*dEdg(2))*rw(j)
   enddo

   j=2
    c(-1)=-0.1111111111111111d0
    c(0)=-1.717857142857143d0
    c(1)=4.d0
    c(2)=-4.666666666666667d0
    c(3)=4.666666666666667d0
    c(4)=-3.5d0
    c(5)=1.866666666666667d0
    c(6)=-0.6666666666666666d0
    c(7)=0.1428571428571428d0
    c(8)=-0.01388888888888889d0
   grad=0.d0
   do i=-1,8
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   grad(2)=grad(2)+c(i)*rho(j+i,2)
   enddo
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   if (grad(2).ge.0.d0) then 
   sign(2)=rd(j)
   else
   sign(2)=-rd(j)
   endif
   grad=sign*grad
   rhoj=rho(j,:)
   call xcfunction(nspol,rhoj,grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*(rho(j,1)+rho(j,2))
   eps(j)=eps(j)+Exc
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   pot(j,2)=pot(j,2)+Vxc(2)*rw(j)
   do i=-1,8
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   pot(j+i,2)=pot(j+i,2)+(sign(2)*c(i)*dEdg(2))*rw(j)
   enddo

   j=3
    c(-2)=0.01111111111111111d0
    c(-1)=-0.2222222222222222d0
    c(0)=-1.217857142857143d0
    c(1)=2.666666666666666d0
    c(2)=-2.333333333333333d0
    c(3)=1.866666666666667d0
    c(4)=-1.166666666666667d0
    c(5)=0.5333333333333333d0
    c(6)=-0.1666666666666666d0
    c(7)=0.03174603174603174d0
    c(8)=-0.2777777777777778d-2
   grad=0.d0
   do i=-2,8
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   grad(2)=grad(2)+c(i)*rho(j+i,2)
   enddo
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   if (grad(2).ge.0.d0) then 
   sign(2)=rd(j)
   else
   sign(2)=-rd(j)
   endif
   grad=sign*grad
   rhoj=rho(j,:)
   call xcfunction(nspol,rhoj,grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*(rho(j,1)+rho(j,2))
   eps(j)=eps(j)+Exc
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   pot(j,2)=pot(j,2)+Vxc(2)*rw(j)
   do i=-2,8
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   pot(j+i,2)=pot(j+i,2)+(sign(2)*c(i)*dEdg(2))*rw(j)
   enddo

   j=4
    c(-3)=-0.202020202020202d-2
    c(-2)=0.03333333333333333d0
    c(-1)=-0.3333333333333333d0
    c(0)=-0.88452380952381d0
    c(1)=2.d0
    c(2)=-1.4d0
    c(3)=0.933333333333333d0
    c(4)=-0.5d0
    c(5)=0.2d0
    c(6)=-0.05555555555555556d0
    c(7)=0.952380952380952d-2
    c(8)=-0.7575757575757577d-3
   grad=0.d0
   do i=-3,8
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   grad(2)=grad(2)+c(i)*rho(j+i,2)
   enddo
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   if (grad(2).ge.0.d0) then 
   sign(2)=rd(j)
   else
   sign(2)=-rd(j)
   endif
   grad=sign*grad
   rhoj=rho(j,:)
   call xcfunction(nspol,rhoj,grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*(rho(j,1)+rho(j,2))
   eps(j)=eps(j)+Exc
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   pot(j,2)=pot(j,2)+Vxc(2)*rw(j)
   do i=-3,8
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   pot(j+i,2)=pot(j+i,2)+(sign(2)*c(i)*dEdg(2))*rw(j)
   enddo

   j=5
    c(-4)=0.5050505050505051d-3
    c(-3)=-0.808080808080808d-2
    c(-2)=0.06666666666666666d0
    c(-1)=-0.4444444444444445d0
    c(0)=-0.6345238095238095d0
    c(1)=1.6d0
    c(2)=-0.933333333333333d0
    c(3)=0.5333333333333333d0
    c(4)=-0.25d0
    c(5)=0.0888888888888889d0
    c(6)=-0.02222222222222222d0
    c(7)=0.3463203463203463d-2
    c(8)=-0.2525252525252525d-3
   grad=0.d0
   do i=-4,8
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   grad(2)=grad(2)+c(i)*rho(j+i,2)
   enddo
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   if (grad(2).ge.0.d0) then 
   sign(2)=rd(j)
   else
   sign(2)=-rd(j)
   endif
   grad=sign*grad
   rhoj=rho(j,:)
   call xcfunction(nspol,rhoj,grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*(rho(j,1)+rho(j,2))
   eps(j)=eps(j)+Exc
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   pot(j,2)=pot(j,2)+Vxc(2)*rw(j)
   do i=-4,8
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   pot(j+i,2)=pot(j+i,2)+(sign(2)*c(i)*dEdg(2))*rw(j)
   enddo

   j=6
    c(-5)=-0.1554001554001554d-3
    c(-4)=0.2525252525252525d-2
    c(-3)=-0.0202020202020202d0
    c(-2)=0.1111111111111111d0
    c(-1)=-0.5555555555555556d0
    c(0)=-0.4345238095238095d0
    c(1)=1.333333333333333d0
    c(2)=-0.6666666666666666d0
    c(3)=0.3333333333333333d0
    c(4)=-0.1388888888888889d0
    c(5)=0.04444444444444445d0
    c(6)=-0.0101010101010101d0
    c(7)=0.1443001443001443d-2
    c(8)=-0.971250971250971d-4
   grad=0.d0
   do i=-5,8
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   grad(2)=grad(2)+c(i)*rho(j+i,2)
   enddo
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   if (grad(2).ge.0.d0) then 
   sign(2)=rd(j)
   else
   sign(2)=-rd(j)
   endif
   grad=sign*grad
   rhoj=rho(j,:)
   call xcfunction(nspol,rhoj,grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*(rho(j,1)+rho(j,2))
   eps(j)=eps(j)+Exc
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   pot(j,2)=pot(j,2)+Vxc(2)*rw(j)
   do i=-5,8
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   pot(j+i,2)=pot(j+i,2)+(sign(2)*c(i)*dEdg(2))*rw(j)
   enddo

   j=7
    c(-6)=0.555000555000555d-4
    c(-5)=-0.932400932400932d-3
    c(-4)=0.7575757575757577d-2
    c(-3)=-0.04040404040404041d0
    c(-2)=0.1666666666666666d0
    c(-1)=-0.6666666666666666d0
    c(0)=-0.2678571428571428d0
    c(1)=1.142857142857143d0
    c(2)=-0.5d0
    c(3)=0.2222222222222222d0
    c(4)=-0.0833333333333333d0
    c(5)=0.02424242424242424d0
    c(6)=-0.5050505050505051d-2
    c(7)=0.6660006660006659d-3
    c(8)=-0.4162504162504162d-4
   grad=0.d0
   do i=-6,8
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   grad(2)=grad(2)+c(i)*rho(j+i,2)
   enddo
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   if (grad(2).ge.0.d0) then 
   sign(2)=rd(j)
   else
   sign(2)=-rd(j)
   endif
   grad=sign*grad
   rhoj=rho(j,:)
   call xcfunction(nspol,rhoj,grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*(rho(j,1)+rho(j,2))
   eps(j)=eps(j)+Exc
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   pot(j,2)=pot(j,2)+Vxc(2)*rw(j)
   do i=-6,8
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   pot(j+i,2)=pot(j+i,2)+(sign(2)*c(i)*dEdg(2))*rw(j)
   enddo

   j=8
    c(-7)=-0.222000222000222d-4
    c(-6)=0.3885003885003884d-3
    c(-5)=-0.3263403263403263d-2
    c(-4)=0.01767676767676768d0
    c(-3)=-0.07070707070707071d0
    c(-2)=0.2333333333333333d0
    c(-1)=-0.7777777777777778d0
    c(0)=-0.125d0
    c(1)=1.d0
    c(2)=-0.3888888888888889d0
    c(3)=0.1555555555555556d0
    c(4)=-0.05303030303030303d0
    c(5)=0.01414141414141414d0
    c(6)=-0.2719502719502719d-2
    c(7)=0.3330003330003329d-3
    c(8)=-0.1942501942501942d-4
   grad=0.d0
   do i=-7,8
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   grad(2)=grad(2)+c(i)*rho(j+i,2)
   enddo
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   if (grad(2).ge.0.d0) then 
   sign(2)=rd(j)
   else
   sign(2)=-rd(j)
   endif
   grad=sign*grad
   rhoj=rho(j,:)
   call xcfunction(nspol,rhoj,grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*(rho(j,1)+rho(j,2))
   eps(j)=eps(j)+Exc
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   pot(j,2)=pot(j,2)+Vxc(2)*rw(j)
   do i=-7,8
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   pot(j+i,2)=pot(j+i,2)+(sign(2)*c(i)*dEdg(2))*rw(j)
   enddo


    c(-8)=9.71250971250971d-6
    c(-7)=-0.1776001776001776d-3
    c(-6)=0.1554001554001554d-2
    c(-5)=-0.87024087024087d-2
    c(-4)=0.3535353535353535d-1
    c(-3)=-0.1131313131313131d0
    c(-2)=0.3111111111111111d0
    c(-1)=-0.888888888888889d0
    c(0)=0.d0
    c(1)=0.888888888888889d0
    c(2)=-0.3111111111111111d0
    c(3)=0.1131313131313131d0
    c(4)=-0.3535353535353535d-1
    c(5)=0.87024087024087d-2
    c(6)=-0.1554001554001554d-2
    c(7)=0.1776001776001776d-3
    c(8)=-9.71250971250971d-6
   do 100,j=9,nrad-8
   grad=0.d0
   do i=-8,8
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   grad(2)=grad(2)+c(i)*rho(j+i,2)
   enddo
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   if (grad(2).ge.0.d0) then 
   sign(2)=rd(j)
   else
   sign(2)=-rd(j)
   endif
   grad=sign*grad
   rhoj=rho(j,:)
   call xcfunction(nspol,rhoj,grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*(rho(j,1)+rho(j,2))
   eps(j)=eps(j)+Exc
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   pot(j,2)=pot(j,2)+Vxc(2)*rw(j)
   do i=-8,8
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   pot(j+i,2)=pot(j+i,2)+(sign(2)*c(i)*dEdg(2))*rw(j)
   enddo
100   continue

   j=nrad-7
   grad=0.d0
   do i=-8,7
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   grad(2)=grad(2)+c(i)*rho(j+i,2)
   enddo
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   if (grad(2).ge.0.d0) then 
   sign(2)=rd(j)
   else
   sign(2)=-rd(j)
   endif
   grad=sign*grad
   rhoj=rho(j,:)
   call xcfunction(nspol,rhoj,grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*(rho(j,1)+rho(j,2))
   eps(j)=eps(j)+Exc
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   pot(j,2)=pot(j,2)+Vxc(2)*rw(j)
   do i=-8,7
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   pot(j+i,2)=pot(j+i,2)+(sign(2)*c(i)*dEdg(2))*rw(j)
   enddo

   j=nrad-6
   grad=0.d0
   do i=-8,6
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   grad(2)=grad(2)+c(i)*rho(j+i,2)
   enddo
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   if (grad(2).ge.0.d0) then 
   sign(2)=rd(j)
   else
   sign(2)=-rd(j)
   endif
   grad=sign*grad
   rhoj=rho(j,:)
   call xcfunction(nspol,rhoj,grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*(rho(j,1)+rho(j,2))
   eps(j)=eps(j)+Exc
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   pot(j,2)=pot(j,2)+Vxc(2)*rw(j)
   do i=-8,6
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   pot(j+i,2)=pot(j+i,2)+(sign(2)*c(i)*dEdg(2))*rw(j)
   enddo

   j=nrad-5
   grad=0.d0
   do i=-8,5
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   grad(2)=grad(2)+c(i)*rho(j+i,2)
   enddo
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   if (grad(2).ge.0.d0) then 
   sign(2)=rd(j)
   else
   sign(2)=-rd(j)
   endif
   grad=sign*grad
   rhoj=rho(j,:)
   call xcfunction(nspol,rhoj,grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*(rho(j,1)+rho(j,2))
   eps(j)=eps(j)+Exc
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   pot(j,2)=pot(j,2)+Vxc(2)*rw(j)
   do i=-8,5
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   pot(j+i,2)=pot(j+i,2)+(sign(2)*c(i)*dEdg(2))*rw(j)
   enddo

   j=nrad-4
   grad=0.d0
   do i=-8,4
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   grad(2)=grad(2)+c(i)*rho(j+i,2)
   enddo
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   if (grad(2).ge.0.d0) then 
   sign(2)=rd(j)
   else
   sign(2)=-rd(j)
   endif
   grad=sign*grad
   rhoj=rho(j,:)
   call xcfunction(nspol,rhoj,grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*(rho(j,1)+rho(j,2))
   eps(j)=eps(j)+Exc
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   pot(j,2)=pot(j,2)+Vxc(2)*rw(j)
   do i=-8,4
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   pot(j+i,2)=pot(j+i,2)+(sign(2)*c(i)*dEdg(2))*rw(j)
   enddo

   j=nrad-3
   grad=0.d0
   do i=-8,3
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   grad(2)=grad(2)+c(i)*rho(j+i,2)
   enddo
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   if (grad(2).ge.0.d0) then 
   sign(2)=rd(j)
   else
   sign(2)=-rd(j)
   endif
   grad=sign*grad
   rhoj=rho(j,:)
   call xcfunction(nspol,rhoj,grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*(rho(j,1)+rho(j,2))
   eps(j)=eps(j)+Exc
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   pot(j,2)=pot(j,2)+Vxc(2)*rw(j)
   do i=-8,3
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   pot(j+i,2)=pot(j+i,2)+(sign(2)*c(i)*dEdg(2))*rw(j)
   enddo

   j=nrad-2
   grad=0.d0
   do i=-8,2
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   grad(2)=grad(2)+c(i)*rho(j+i,2)
   enddo
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   if (grad(2).ge.0.d0) then 
   sign(2)=rd(j)
   else
   sign(2)=-rd(j)
   endif
   grad=sign*grad
   rhoj=rho(j,:)
   call xcfunction(nspol,rhoj,grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*(rho(j,1)+rho(j,2))
   eps(j)=eps(j)+Exc
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   pot(j,2)=pot(j,2)+Vxc(2)*rw(j)
   do i=-8,2
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   pot(j+i,2)=pot(j+i,2)+(sign(2)*c(i)*dEdg(2))*rw(j)
   enddo

   j=nrad-1
   grad=0.d0
   do i=-8,1
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   grad(2)=grad(2)+c(i)*rho(j+i,2)
   enddo
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   if (grad(2).ge.0.d0) then 
   sign(2)=rd(j)
   else
   sign(2)=-rd(j)
   endif
   grad=sign*grad
   rhoj=rho(j,:)
   call xcfunction(nspol,rhoj,grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*(rho(j,1)+rho(j,2))
   eps(j)=eps(j)+Exc
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   pot(j,2)=pot(j,2)+Vxc(2)*rw(j)
   do i=-8,1
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   pot(j+i,2)=pot(j+i,2)+(sign(2)*c(i)*dEdg(2))*rw(j)
   enddo

   j=nrad-0
   grad=0.d0
   do i=-8,0
   grad(1)=grad(1)+c(i)*rho(j+i,1)
   grad(2)=grad(2)+c(i)*rho(j+i,2)
   enddo
   if (grad(1).ge.0.d0) then 
   sign(1)=rd(j)
   else
   sign(1)=-rd(j)
   endif
   if (grad(2).ge.0.d0) then 
   sign(2)=rd(j)
   else
   sign(2)=-rd(j)
   endif
   grad=sign*grad
   rhoj=rho(j,:)
   call xcfunction(nspol,rhoj,grad,Exc,Vxc,dEdg)
   enexc=enexc+Exc*rw(j)*(rho(j,1)+rho(j,2))
   eps(j)=eps(j)+Exc
   pot(j,1)=pot(j,1)+Vxc(1)*rw(j)
   pot(j,2)=pot(j,2)+Vxc(2)*rw(j)
   do i=-8,0
   pot(j+i,1)=pot(j+i,1)+(sign(1)*c(i)*dEdg(1))*rw(j)
   pot(j+i,2)=pot(j+i,2)+(sign(2)*c(i)*dEdg(2))*rw(j)
   enddo

   do j=2,nrad
      pot(j,1)=pot(j,1)/rw(j)
      pot(j,2)=pot(j,2)/rw(j)
   enddo
   pot(1,:)=pot(2,:)

end subroutine driveGGApolarized
