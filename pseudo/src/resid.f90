!> @file
!! atomic program for generating and optimizing hgh pseudo-potentials.
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


!> Spin polarized treatment if nspol == 2
subroutine resid(nspol, &
   noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,nsmx, &
   aeval,res, &
   hsep, &
   ud,nint,ng,ngmx,psi,rho,pp1,pp2,pp3, &
   potgrd,pexgrd,xcgrd,rr,rw, &
   ppr1,ppr2,ppr3,aux1,aux2, &
   expxpr)

   implicit none
   !Arguments
   integer, intent(in) :: nspol,noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,nsmx
   integer, intent(in) :: nint,ng,ngmx
   real(kind=8), dimension(noccmx,lmx,nsmx), intent(in) :: aeval
   real(kind=8), dimension(noccmx,lmx,nsmx), intent(inout) :: res
   real(kind=8), dimension(6,lpmx,nsmx), intent(in) :: hsep
   real(kind=8), dimension(nint,((ng+1)*(ng+2))/2,lcx+1), intent(in) :: ud
   real(kind=8), dimension(0:ngmx,noccmx,lmx,nsmx), intent(in) :: psi
   real(kind=8), dimension(((ng+1)*(ng+2))/2,lmax+1,nspol), intent(in) :: rho
   real(kind=8), dimension(0:ng,lmax+1), intent(in) :: pp1,pp2,pp3
   real(kind=8), dimension(nint), intent(out) :: potgrd
   real(kind=8), dimension(nint), intent(in) :: pexgrd
   real(kind=8), dimension(nint,nspol), intent(in) :: xcgrd !< Vxc (ij, ispin) ; xcgrd is actually a waste of memory
   real(kind=8), dimension(nint), intent(in) :: rr,rw
   real(kind=8), dimension(nint,lmax), intent(in) :: ppr1,ppr2,ppr3
   real(kind=8), dimension(nint), intent(in) :: aux1
   real(kind=8), dimension(nint,0:ng,lmax+1), intent(in) :: aux2
   real(kind=8), dimension(0:ng,nint), intent(in) :: expxpr
   !Local variables
   real(kind=8), parameter :: fourpi = 16.d0*atan(1.d0)

   logical :: nothing
   real(kind=8), dimension(nsmx) :: scpr1,scpr2,scpr3
   real(kind=8) :: psigrd,r,rkin,sep,tt
   real(kind=8), external :: gamma,ddot,wave2
   integer :: i,iocc,ispin,k,ll

   nothing = .true.
   do ll=0,lmax
      do ispin=1,max(nspol,min(2*ll+1,nspin))
         do iocc=1,noccmax
            if (res(iocc,ll+1,ispin).ne.-1.0d0) nothing=.false.
         enddo
      enddo
   enddo
   if (nothing)  return
   
   !   external and exc potential on grid 
   do k=1,nint
      r=rr(k)
      ! some older comment line?
      ! potgrd(k)=pexgrd(k) + fourpi*xcgrd(k)/rw(k)
      ! nonpolarized version was
      ! potgrd(k)=pexgrd(k) + aux1(k)*xcgrd(k)
      potgrd(k)=pexgrd(k) 
      ! the missing XC term is added later to the temp variable tt
   enddo

   ! add hartree potential
   ! potgrd(k)=+ud(k,i,j,l+1)*rho(i,j,l+1)
   call DGEMV('N',nint,((ng+1)*(ng+2))/2*(lcx+1),1.d0,ud,nint, rho(1,1,1),1,1.d0,potgrd,1)

   ! polarized: both channels of rho are needed for the Hartree term
   if (nspol.eq.2) &
      call DGEMV('N',nint,((ng+1)*(ng+2))/2*(lcx+1),1.d0,ud,nint, &
                rho(1,1,2),1,1.d0,potgrd,1)

   do ll=0,lmax
      ! if nspol=2, s orbitals have two spin states, too
      do ispin=1,max(min(2*ll+1,nspin),nspol)
         do iocc=1,noccmax
            if (res(iocc,ll+1,ispin).ne.-1.0d0) then
               ! separable part
               if (ll.le.lpx-1) then
                  scpr1(ispin)=DDOT(ng+1,psi(0,iocc,ll+1,ispin), &
                       1,pp1(0,ll+1),1)
                  scpr2(ispin)=DDOT(ng+1,psi(0,iocc,ll+1,ispin), &
                       1,pp2(0,ll+1),1)
                  scpr3(ispin)=DDOT(ng+1,psi(0,iocc,ll+1,ispin), &
                       1,pp3(0,ll+1),1)
               endif
               res(iocc,ll+1,ispin)=0.d0
               do k=1,nint
                  r=rr(k)
                  ! wavefunction on grid
                  psigrd= wave2(ng,ll,psi(0,iocc,ll+1,ispin), &
                       expxpr,r,k,nint)
                  ! kinetic energy
                  rkin=0.d0
                  do i=0,ng
                     rkin=rkin+psi(i,iocc,ll+1,ispin) &
                          *aux2(k,i,ll+1)
                  enddo
                  rkin=rkin*r**ll
                  ! separable part
                  if (ll.le.lpx-1) then
                     sep = (scpr1(ispin)*hsep(1,ll+1,ispin) &
                          + scpr2(ispin)*hsep(2,ll+1,ispin) &
                          + scpr3(ispin)*hsep(4,ll+1,ispin)) &
                          *ppr1(k,ll+1)+ &
                          (scpr1(ispin)*hsep(2,ll+1,ispin) & 
                          + scpr2(ispin)*hsep(3,ll+1,ispin) & 
                          + scpr3(ispin)*hsep(5,ll+1,ispin)) &
                          *ppr2(k,ll+1) &
                          +(scpr1(ispin)*hsep(4,ll+1,ispin) & 
                          + scpr2(ispin)*hsep(5,ll+1,ispin) &
                          + scpr3(ispin)*hsep(6,ll+1,ispin)) &
                          *ppr3(k,ll+1)
                  else
                     sep=0.d0
                  endif
                  ! residue
                  ! the term from VXC can be spin polarized,
                  ! therefore  the index min(nspol, ispin)
                  tt=aux1(k)*xcgrd(k,min(nspol,ispin))
                  ! add external and Ha potential, KS ev
                  tt=tt+potgrd(k)-aeval(iocc,ll+1,ispin)
                  ! times psi, plus kinetic term and separable part
                  tt=psigrd*tt+rkin+sep
                  ! write(22,*)'DEBUG resid:iocc,ll,ispin,tt,rkin,sep', iocc,ll,ispin, &
                  !    (tt-rkin-sep)*rw(k) ,rkin*rw(k),sep*rw(k)
                  res(iocc,ll+1,ispin)=res(iocc,ll+1,ispin)+tt*tt*rw(k)/fourpi
               enddo
            else
               res(iocc,ll+1,ispin)=0.0d0
            endif
         enddo
      enddo
   enddo
      
end subroutine resid
