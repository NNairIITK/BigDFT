!> @file
!! Paw generation using dirac equation (?) in pseudo program
!! @author
!!    Alex Willand, under the supervision of Stefan Goedecker
!!    gpu accelerated routines by Raffael Widmer
!!    parts of this program were based on the fitting program by matthias krack
!!    http://cvs.berlios.de/cgi-bin/viewcvs.cgi/cp2k/potentials/goedecker/pseudo/v2.2/
!!
!!    Copyright (C) 2010-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> difnrl integrates the Schroedinger equation
!! if finds the eigenvalue ev, the wavefunction ar
!! and the derivative br = d(ar)/dr
      subroutine difnrl(v,ar,br,&
       nr,a,b,r,rab,np,lo,znuc,&
       ev)
      implicit real*8 (a-h,o-z)
      !> Tolerance
      real(kind=8), parameter :: etol=-1.d-7
      real(kind=8), parameter :: tol=1.0d-14
!
      dimension v(nr),ar(nr),br(nr),r(nr),rab(nr) 

!
!    Arrays added to gain speed.
!
      dimension rabrlo(5),rlp(5),rab2(nr),fa(nr),fb(nr)
!
!------Machine dependent parameter-
!------Require exp(-2*expzer) to be within the range of the machine
! IBM
      expzer = 3.7D2
!Iris     expzer = 3.7E2
!Apollo   expzer = 3.7D2
!Sun      expzer = 3.7D2
!Vax      expzer = 44.D0
!ray      expzer =  2.8E3
!
!     for numerical stability:
!
      expzer = expzer/2
!
!      integration coefficients
!
       abc1 = 1901.D0/720.D0
       abc2 = -1387.D0/360.D0
       abc3 = 109.D0/30.D0
       abc4 = -637.D0/360.D0
       abc5 = 251.D0/720.D0
       amc0 = 251.D0/720.D0
       amc1 = 323.D0/360.D0
       amc2 = -11.D0/30.D0
       amc3 = 53.D0/360.D0
       amc4 = -19.D0/720.D0
      itmax = 100
      lp = lo+1
      ar(1) = 0.0d0
      if (lo == 0) then
        br(1) = b*a
      else
        br(1) = 0.0d0
      endif
      do 1 j=2,nr
        ar(j) = 0.0d0
 1    continue
      do 2 j=2,nr
        br(j) =0.0d0
 2    continue
      do 4 j=2,5
        rlp(j)=r(j)**lp
 4    continue
      do 5 j=2,5
        rabrlo(j)=rab(j)*r(j)**lo
 5    continue
      do 6 j=1,nr
        rab2(j)=rab(j)*rab(j)
 6    continue
!
!   set underflow trap
!
      juflow=1
      do 42 j=2,nr
        if (lp*abs(log(r(j))) .ge. expzer/2) juflow = j
 42   continue
!
!   determine effective charge and vzero for startup of
!   outward integration
!   ar = r**(l+1) * (1 + aa r + bb r**2 + ... )
!   aa = -znuc / lp     bb = (-2 znuc aa + v(0) - e)/(4 l + 6)
!
      zeff=znuc
      aa = -zeff/lp
      vzero = -2*zeff*aa

      var0 = 0.0d0
      if (lo == 0) var0=-2*zeff
      if (lo == 1) var0=2.0d0
      emax = 0.0d0
      emin = -200000.0d0


      if (ev .gt. emax) ev = emax
      !! write(6,15) ev,nodes
 10   continue
       !! if (itmax .lt. 2) write(6,15) ev,nodes
       !! 15   format(' ev =',1pe18.10,' nodes =',i2)

      if (itmax == 0) return
      if (ev .gt. 0.0) then
        write(6,1000)iorb
        stop 'difnrl one'
      endif
 1000 format(//,' error in difnrl - ev(',i2,&
       ') greater then v(infinty)')
!
!   find practical infinity ninf and classical turning
!   point nctp for orbital
!

      icount=0
 20   continue
      icount=icount+1
      do j=nr,2,-1
        temp = v(j) -ev
        if (temp .lt. 0.0) temp = 0.0d0
        if (r(j)*sqrt(temp) .lt. expzer) goto 23
      end do
 23   continue
      ninf=j
      nctp = ninf - 5
      do j=2,ninf-5
        if (v(j) .lt. ev) nctp = j
      end do
      if (ev .ge. etol*10) nctp=ninf-5
      if (ev .ge. etol) ev=0.0d0
      if (nctp .le. 6) then
        ev = 0.9d0*ev
        if (icount .gt. 100) then
          write(*,*)
          write(6,1010)
          stop 'difnrl two'
        endif
        goto 20
      endif
 1010 format(//,'error in difnrl - cannot find the classical turning point ')
!
!   outward integration from 1 to nctp
!   startup
!
      bb = (vzero-ev)/(4*lp+2)
      do 35 j=2,5
        ar(j) = rlp(j) * (1+(aa+bb*r(j))*r(j))
        br(j) = rabrlo(j) * (lp+(aa*(lp+1)+bb*(lp+2)*r(j))*r(j))
 35   continue
!
!    Predictor-corrector array added.
!
      fa(1) = br(1)
      fb(1) = b*br(1) + rab2(1)*var0
      fa(2) = br(2)
      fb(2) = b*br(2) + rab2(2)*(v(2)-ev )*ar(2)
      fa(3) = br(3)
      fb(3) = b*br(3) + rab2(3)*(v(3)-ev )*ar(3)
      fa(4) = br(4)
      fb(4) = b*br(4) + rab2(4)*(v(4)-ev )*ar(4)
      fa(5) = br(5)
      fb(5) = b*br(5) + rab2(5)*(v(5)-ev )*ar(5)
!
!   intergration loop
!
      nodes = 0
      do 40 j=6,nctp
!
!   predictor (Adams-Bashforth)
!
        j1=j-1
        j2=j-2
        j3=j-3
        j4=j-4
        j5=j-5
        vev=v(j)-ev
        arp = ar(j1) + abc1*fa(j1)+abc2*fa(j2)+abc3*fa(j3)+&
         abc4*fa(j4)+abc5*fa(j5)
        brp = br(j1) + abc1*fb(j1)+abc2*fb(j2)+abc3*fb(j3)+&
         abc4*fb(j4)+abc5*fb(j5)
        fb1 = b*brp + rab2(j)*vev*arp
!
!   corrector (Adams-Moulton)
!
        arc = ar(j1) + amc0*brp+amc1*fa(j1)+amc2*fa(j2)+&
         amc3*fa(j3)+amc4*fa(j4)
        brc = br(j1) + amc0*fb1+amc1*fb(j1)+amc2*fb(j2)+&
         amc3*fb(j3)+amc4*fb(j4)
        fb0 = b*brc + rab2(j)*vev*arc
!
!   error reduction step
!
        ar(j) = arc + amc0*(brc-brp)
        br(j) = brc + amc0*(fb0-fb1)
        fa(j) = br(j)
        fb(j) = b*br(j) + rab2(j)*vev*ar(j)
!
!   count nodes - if no underflow
!
        if(j.gt.juflow.and.ar(j)*ar(j-1).lt.0.0)nodes=nodes+1
 40   continue
!
      arctp = ar(nctp)
      brctp = br(nctp)
!
!   end outward integration
!
!   if number of nodes correct, start inward integration
!   else modify energy stepwise and try again
!
      if (nodes /= np-lo-1) then
!     c.hartwig
        !! write(6,*) 'nodes,ev',nodes,ev
        if (nodes .lt. np-lo-1) then
!
!  too few nodes; increase ev
!
          if (ev .gt. emin) emin = ev
          ev = ev - ev/10
        else
!
!  too many nodes; decrease ev
!
          if (ev .lt. emax) emax = ev
          ev = ev + ev/10
        endif
        itmax = itmax-1
        goto 10
      endif
!
!   inward integration from ninf to nctp
!   startup
!
      do 71 j=ninf,ninf-4,-1
        alf = v(j) - ev
        if (alf .lt. 0.0) alf = 0.0d0
        alf = sqrt(alf)
        ar(j) = exp(-alf*r(j))
        br(j) = -rab(j)*alf*ar(j)
 71   continue
!
!    Array for predictor-corrector added.
!
      fa(ninf) = br(ninf)
      fb(ninf) = b*br(ninf) + rab2(ninf)*&
       (v(ninf)-ev)*ar(ninf)
      ninf1 = ninf - 1
      fa(ninf1) = br(ninf1)
      fb(ninf1) = b*br(ninf1) + rab2(ninf1)*&
             (v(ninf1)-ev)*ar(ninf1)
      ninf2 = ninf - 2
      fa(ninf2) = br(ninf2)
      fb(ninf2) = b*br(ninf2) + rab2(ninf2)*&
             (v(ninf2)-ev)*ar(ninf2)
      ninf3 = ninf - 3
      fa(ninf3) = br(ninf3)
      fb(ninf3) = b*br(ninf3) + rab2(ninf3)*&
             (v(ninf3)-ev)*ar(ninf3)
      ninf4 = ninf - 4
      fa(ninf4) = br(ninf4)
      fb(ninf4) = b*br(ninf4) + rab2(ninf4)*&
             (v(ninf4)-ev)*ar(ninf4)
!
!   integration loop
!
      istop = ninf - nctp
      if (istop .lt. 5) goto 222
      do 80 j=ninf-5,nctp,-1
!
!   predictor (Adams-Bashforth)
!
        j1 = j + 1
        j2 = j + 2
        j3 = j + 3
        j4 = j + 4
        j5 = j + 5
        vev = v(j)-ev
        arp = ar(j1) - (abc1*fa(j1)+abc2*fa(j2)+abc3*fa(j3)+&
         abc4*fa(j4)+abc5*fa(j5))
        brp = br(j1) - (abc1*fb(j1)+abc2*fb(j2)+abc3*fb(j3)+&
         abc4*fb(j4)+abc5*fb(j5))
        fb0 = b*brp + rab2(j)*vev*arp
!
!   corrector (Adams-Moulton)
!
        arc = ar(j1) - (amc0*brp+amc1*fa(j1)+amc2*fa(j2)+&
         amc3*fa(j3)+amc4*fa(j4))
        brc = br(j1) - (amc0*fb0+amc1*fb(j1)+amc2*fb(j2)+&
         amc3*fb(j3)+amc4*fb(j4))
!
        fb1 = b*brc + rab2(j)*vev*arc
!
!   error reduction step
!
        ar(j) = arc - amc0*(brc-brp)
        br(j) = brc - amc0*(fb1-fb0)
        fa(j) = br(j)
        fb(j) = b*br(j) + rab2(j)*vev*ar(j)
 80   continue
!
!   end inward integration
!
!   rescale ar and br outside nctp to match ar(nctp) from
!   outward integration
!
  222 factor = arctp/ar(nctp)
      do 90 j=nctp,ninf
        ar(j) = factor * ar(j)
        br(j) = factor * br(j)
 90   continue
!
!   find normalizing factor
!
      factor = 0.0d0
      ll = 4
      do 100 j=2,ninf
        factor = factor + ll*ar(j)*ar(j)*rab(j)
        ll = 6 - ll
 100  continue
      factor = factor / 3
!
!   modify eigenvalue ev
!
      dev = arctp * (brctp-br(nctp)) / (factor * rab(nctp))
      if (5*abs(dev) .gt. -ev) dev=sign(ev,dev)/5
      itmax = itmax-1
      evold = ev
      ev = ev + dev
      if (ev .gt. emax) ev = (evold + emax) / 2
      if (ev .lt. emin) ev = (evold + emin) / 2
      if (abs(dev) .gt. tol*(1-ev)) goto 10
!
!   normalize wavefunction and change br from d(ar)/dj to d(ar)/dr
!
      factor = 1 / sqrt(factor)
      do 110 j=1,ninf
        ar(j) = factor*ar(j)
        br(j) = factor*br(j) / rab(j)
 110  continue
      end
!
!      *****************************************************************
!
 
