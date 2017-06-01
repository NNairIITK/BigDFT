!> @file
!! Part of the pseudo program (pseudopotential generation)
!! @author
!!    Alex Willand, under the supervision of Stefan Goedecker
!!    gpu accelerated routines by Raffael Widmer
!!    parts of this program were based on the fitting program by Matthias Krack
!!    http://cvs.berlios.de/cgi-bin/viewcvs.cgi/cp2k/potentials/goedecker/pseudo/v2.2/
!!
!!    Copyright (C) 2010-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Calculation of the total energy (pseudo program)
subroutine etot(verbose,nspol,  &
        noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,nsmx,  &
        aeval,  &
        rprb,zion,rloc,gpot,r_l,r_l2,hsep,  &
        xp,ud,nint,ng,ngmx,psi,rho,pp1,pp2,pp3,  &
        xcgrd,excgrd,rhogrd,rhocore,occup,rr,rw,  &
        expxpr,etotal)
   
   implicit none
   real(kind=8), parameter :: fourpi = 16.d0*atan(1.d0)
   !Arguments
   logical, intent(in) :: verbose
   integer, intent(in) :: nspol,noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,nsmx
   integer, intent(in) :: ng,nint,ngmx
   real(kind=8), intent(in) :: rprb,zion,rloc
   real(kind=8), dimension(noccmx,lmx,nsmx), intent(in) :: aeval,occup
   real(kind=8), dimension(4), intent(in) :: gpot
   real(kind=8), dimension(nint,nspol), intent(in) :: xcgrd, rhocore
   real(kind=8), dimension(nint,nspol), intent(inout) :: rhogrd
   real(kind=8), dimension(0:ng,lmax+1), intent(in) :: pp1,pp2,pp3
   real(kind=8), dimension(lmx), intent(in) :: r_l,r_l2
   real(kind=8), dimension(6,lpmx,nsmx), intent(in) :: hsep
   real(kind=8), intent(out) :: etotal
   !Local variables
   real(kind=8), dimension(0:ng,nint) :: expxpr
   real(kind=8), dimension(0:ng) :: xp
   real(kind=8), dimension(nint,((ng+1)*(ng+2))/2,lcx+1) :: ud
   real(kind=8), dimension(0:ngmx,noccmx,lmx,nsmx) :: psi
   real(kind=8), dimension(((ng+1)*(ng+2))/2,lmax+1,nspol) :: rho
   real(kind=8), dimension(nsmx) :: scpr1,scpr2,scpr3
   real(kind=8), dimension(nint) :: vexgrd,excgrd,vhgrd,rr,rw
   real(kind=8) :: dd,denfull,eext,eh,ehart,eigsum,ekin,enl,epot,exc
   real(kind=8) :: psigrd,r,rkin,rnrm1,rnrm2,rnrm3,sep,vxc,zz
   real(kind=8) :: ddot,wave2
   integer :: i,iocc,ispin,ispol,k,l,ll
   
   real(kind=8), external :: gamma
   
   exc  = 0.0d0
   eext = 0.0d0
   ekin = 0.0d0
   vxc  = 0.0d0
   ehart= 0.0d0
   eigsum=0.0d0
   enl  = 0.0d0
   
   !     vxc integral is with rho, exc integral with rho+rhocore
   
   !     calc exc and vxc first
   !     rhogrd already contains the core charge at this point in case of a nlcc.
   !     note: vxcgrd as returned by ggaenergy was already multifplied by rw/4pi in gatom
   !     use the same loop to strip off the core charge 
   !     and to add up the charges from  both channels.
   !     we do not need polarization nor core charge for the other energy terms!
   do ispol=1,nspol
      do k=1,nint
         exc   = exc   +    excgrd(k)      *rhogrd(k,ispol)*rw(k)
         rhogrd(k,ispol)=rhogrd(k,ispol)-rhocore(k,ispol)        
         vxc   = vxc   +     xcgrd(k,ispol)*rhogrd(k,ispol)
      end do
   end do
   
   
   if(nspol.eq.2)then
      do k=1,((ng+1)*(ng+2))/2
         do l=1,lcx+1
            rho(k,l,1)=rho(k,l,1)+rho(k,l,2)
         end do
      end do
      do k=1,nint
         rhogrd(k,1)=rhogrd(k,1)+rhogrd(k,2)
      end do
   end if
   
   ! calc hartree potential
   
   
   !     first and possibly only spin channel:
   call dgemv('n',nint,((ng+1)*(ng+2))/2*(lcx+1),1.d0,ud,nint,  &
        rho(1,1,1),1,0.d0,vhgrd,1)
   
   
   
   
   !   calc eext and ehart
   do k=1,nint
      r=rr(k)
      vexgrd(k)=.5d0*(r/rprb**2)**2-zion*derf(r/(sqrt(2.d0)*rloc))/r   &
           + exp(-.5d0*(r/rloc)**2)*  &
           ( gpot(1) + gpot(2)*(r/rloc)**2 + gpot(3)*(r/rloc)**4 +   &
           gpot(4)*(r/rloc)**6 )
      !        write(88,*)rr(k),vhgrd(k),rhogrd(k,1)+rhogrd(k,2),rw(k)
      eext  = eext  +      vexgrd(k)*rhogrd(k,1)*rw(k)
      ehart = ehart + 0.5d0*vhgrd(k)*rhogrd(k,1)*rw(k)
      
   enddo
   
   !     calc the kinetic energy term and the separable part of the pseudopotential
   
   
   vxc=vxc*fourpi
   do ll=0,lmax
      do ispin=1,max(min(2*ll+1,nspin),nspol)
         if (ll.le.lpx-1) then
            !       note: here sigma= r_l2(l) is used for the pp2 and pp3 projcetors
            rnrm1=1.d0/sqrt(.5d0*gamma(ll+1.5d0)*r_l(ll+1)**(2*ll+3))
            rnrm2=1.d0/sqrt(.5d0*gamma(ll+3.5d0)*r_l2(ll+1)**(2*ll+7))
            rnrm3=1.d0/sqrt(.5d0*gamma(ll+5.5d0)*r_l2(ll+1)**(2*ll+11))
         endif
         do iocc=1,noccmax
            zz = occup(iocc,ll+1,ispin)
            if (zz.gt.1.0d-8  ) then 
               eigsum=eigsum + aeval(iocc,ll+1,ispin) *zz
               !     separabel part
               if (ll.le.lpx-1) then
                  scpr1(ispin)=ddot(ng+1,psi(0,iocc,ll+1,ispin),  &
                       1,pp1(0,ll+1),1)
                  scpr2(ispin)=ddot(ng+1,psi(0,iocc,ll+1,ispin),  &
                       1,pp2(0,ll+1),1)
                  scpr3(ispin)=ddot(ng+1,psi(0,iocc,ll+1,ispin),  &
                       1,pp3(0,ll+1),1)
               endif
               dd=ekin
               do k=1,nint
                  r=rr(k)
                  !     wavefunction on grid
                  psigrd = wave2(ng,ll,psi(0,iocc,ll+1,ispin),  &
                       expxpr,r,k,nint)
                  !     kinetic energy
                  rkin=0.d0
                  do i=0,ng
                     rkin=rkin+psi(i,iocc,ll+1,ispin)*(xp(i)*  &
                          (3.d0+2.d0*ll-2.d0*xp(i)*r**2)  &
                          *exp(-xp(i)*r**2))
                  enddo
                  !                    below line seems to be nonsense!
                  !                    rkin=rkin*r**ll
                  ekin = ekin+rkin*(psigrd*r**ll)*zz*rw(k)/fourpi
                  !     separabel part
                  if (ll.le.lpx-1) then
                     sep = (scpr1(ispin)*hsep(1,ll+1,ispin)   &
                          + scpr2(ispin)*hsep(2,ll+1,ispin)   &
                          + scpr3(ispin)*hsep(4,ll+1,ispin))  &
                          *rnrm1*r**ll*exp(-.5d0*(r/r_l(ll+1))**2)+  &
                          (scpr1(ispin)*hsep(2,ll+1,ispin)   &
                          + scpr2(ispin)*hsep(3,ll+1,ispin)   &
                          + scpr3(ispin)*hsep(5,ll+1,ispin))  &
                          *rnrm2*r**(ll+2)  &
                          *exp(-.5d0*(r/r_l2(ll+1))**2)     &
                          +(scpr1(ispin)*hsep(4,ll+1,ispin)   &
                          + scpr2(ispin)*hsep(5,ll+1,ispin)   &
                          + scpr3(ispin)*hsep(6,ll+1,ispin))  &
                          *rnrm3*r**(ll+4)  &
                          *exp(-.5d0*(r/r_l2(ll+1))**2)
                     enl = enl+sep*(psigrd*r**ll)*zz*rw(k)/fourpi
                     
                  else
                     sep=0.d0
                  endif
               enddo
               !             write(*,'(a,f20.12,3i3,f10.3)')
               !    :       'DEBUG: ekin component',(ekin-dd)/zz,ll,ispin,iocc,zz
            endif
         enddo
      enddo
   enddo
   eh     =          ehart + vxc - exc - enl
   etotal = eigsum - ehart - vxc + exc
   epot   = etotal - ekin  + eh !  -.5*enl
   if(verbose)then ! condition in previos versions was energ
      write(6,*) 
      write(6,*)' pseudo atom energies'
      write(6,*)' ____________________'
      write(6,*) 
      write(6,'(a,f16.10)')' kinetic energy            =',ekin
      write(6,'(a,f16.10)')' potential energy          =',epot
      write(6,'(a,f16.10)')' non local energy          =',enl
      write(6,'(a,f16.10)')' hartree energy with xc    =',eh
      write(6,'(a,f16.10)')' exchange + corr energy    =',exc
      write(6,'(a,f16.10)')' vxc    correction         =',vxc
      write(6,'(a)')' -------------------------------------------'
      write(6,'(a,f16.10)')' el-el  interaction energy =',ehart
      write(6,'(a,f16.10)')' external energy           =',eext
      write(6,'(a,f16.10)')' sum of eigenvalues        =',eigsum
      write(6,'(a)')' -------------------------------------------'
      write(6,'(a,f16.10)')' total energy              =',etotal
      denfull=0.0d0
      write(6,*)
      do k=1,nint
         denfull=denfull+rhogrd(k,1)*rw(k)
      enddo
      write(6,'(a,e10.4)')' atomic+electronic charge ',zion-denfull
      write(6,*)
   end if
   return
end
