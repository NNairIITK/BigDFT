       subroutine velect(iter,iconv,iXC,ispp,nspol,ifcore,
     1 nrmax,nr,a,b,r,rab,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zsh,rsh,zel,zcore,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
       implicit double precision(a-h,o-z)
c
c      velect generates the electronic output potential from
c      the electron charge density.
c      the ionic part is added in dsolve.
c
       dimension r(nr),rab(nr),
     2 no(norb),lo(norb),so(norb),zo(norb),
     3 cdd(nr),cdu(nr),cdc(nr),
     4 viod(lmax,nr),viou(lmax,nr),vid(nr),viu(nr),vod(nr),vou(nr),
     5 etot(10),ev(norb),ek(norb),ep(norb)
       dimension vtemp(1000)
       character*2 ispp*1,nameat,itype
       integer:: iXC
c
      parameter ( mesh = 2000 )
       dimension y(mesh),yp(mesh),ypp(mesh),w(3*mesh),s1(mesh),s2(mesh)
       common  y,yp,ypp,w,s1,s2
c
c      for use in routine atomwr:
       parameter (ntitle = 40)
       character*40 text(ntitle)
       character irel*3, xccore*4, cdtyp*2

c     c.hartwig
c     convention for spol inded as in dsolv: spin down=1 and spin up=2
      dimension rho(nr,nspol),excgrd(nr),vxcgrd(nr,nspol)
      dimension rw(10000),rd(10000)
      common /intgrd/ rw,rd
      INCLUDE 'func.inc'

c
c
       pi = 4*atan(1.D0)
c
c      fit cd/r by splines
c
       y(1) = 0.D0
       do 10 i=2,nr
       y(i) = (cdd(i)+cdu(i))/r(i)
c      below test output proofs cdd cdu are charge densities
c      write(22,'(i4,3f20.8)')i,r(i),cdu(i),cdd(i)
       if (ifcore .eq. 2) y(i) = y(i) + cdc(i)/r(i)
 10    continue
       isx = 0
       a1 = 0.D0
       an = 0.D0
       b1 = 0.D0
       bn = 0.D0
       call splift(r,y,yp,ypp,nr,w,ierr,isx,a1,b1,an,bn)
c
c      compute the integrals of cd/r and cd from
c      r(1)=0 to r(i)
c
       xlo = 0.D0
       call spliq(r,y,yp,ypp,nr,xlo,r,nr,s2,ierr)
c      s2 ==    ans(i) = integral from xlo to xup(i)
       do 20 i=1,nr
       ypp(i) = r(i)*ypp(i) + 2*yp(i)
       yp(i)  = r(i)*yp(i)  + y(i)
       y(i)   = r(i)*y(i)
 20    continue
       call spliq(r,y,yp,ypp,nr,xlo,r,nr,s1,ierr)
c
c      check normalization
c
       xnorm = 0.D0
       if (zel .ne. 0.D0) xnorm = zel/s1(nr)


c      let us try this
       if (iter .gt. 3 .and. abs(zel-s1(nr)) .gt. 0.01) then
         if (zel .lt. s1(nr)+1.0 ) then
           write(6,24) iter,xnorm
 24    format(/,' warning *** charge density rescaled in',
     1 ' velect',/,' iteration number',i4,3x,
     2 'scaling factor =',f6.3,/)
         else
           xnorm=.99d0*xnorm
           write(6,25) iter,xnorm
 25    format(/,' warning *** charge density partially rescaled in',
     1 ' velect',/,' iteration number',i4,3x,
     2 'scaling factor =',f6.3,/)
         endif
       endif


c      rather than:
c      if (iter .gt. 0 .and. abs(zel-s1(nr)) .gt. 0.01D0)
c    1 write(6,25) iter,xnorm
c25    format(/,46h warning *** charge density rescaled in velect,
c    1 /,17h iteration number,i4,3x,16hscaling factor =,g10.3,/)



c
c      compute new hartree potential
c      renormalize the charge density
c
       do 30 i=2,nr
c      modification needed?
c      looks like at this point, no difference btw up and dwn V is needed
c      OR is this exactly the unscreening effect we are MISSING?
       vod(i) = 2 * xnorm*(s1(i)/r(i) + s2(nr) - s2(i))
       vou(i) = vod(i)    
       cdd(i) = xnorm*cdd(i)*(zel+1)/2/zel
       cdu(i) = xnorm*cdu(i)
c      brute guess for unscreening 
c      vou(i) = vou(i)  -cdu(i)/r(i)  
c      vod(i) = vod(i)  -cdd(i)/r(i)  
 30    continue
c
       if (iconv .ne. 1) goto 50
c
c      compute hartree contribution to total energy
c      does not look spin polarized yet
c
       ehart = 0.D0
       ll = 4
       do 40 i=2,nr
       ehart = ehart + ll * (cdd(i)+cdu(i)) * vod(i) * rab(i)
c      ehart = ehart + ll * (cdd(i)*vod(i)+cdu(i)*vod(i))* rab(i)
c      ^^ identical if vod=vou, which is the case for now

       ll = 6 - ll
 40    continue
       ehart = ehart / 6
c
c      find derivatives of the charge density
c
       do 45 i=2,nr
c      ??????????????????????????????????????????
c      probably moved to ggaenergy17
 45    continue
c
c      store the atomic Coulomb (ionic + Hartree) potential on file
c
c      first construct the total potential, store in array vtemp:
c
c       ifile = 2
c       irectp = 31
       do 300 l = 1, 3
c       do 310 i = 1, nr
c         vtemp(i) = viod(l,i) + vod(i) * r(i)
c310    continue
c       cdtyp = ' '
c       itype = ' '
c       call atomwr
c     +  (ifile,irectp,nameat,iXC,irel,xccore,zcore,norb,text,
c     +   nr,aa,bb,r,nql,delql,nqnl,delqnl,numnl,
c     +   itype,cdtyp,0,(l-1),mode,vtemp)
c
c       if (ispp .ne. ' ') then
c         do 220 i = 1, nr
c           vtemp(i) = viou(l,i) + vou(i) * r(i)
c220      continue
c         cdtyp = 'up'
c         call atomwr
c     +    (ifile,irectp,nameat,iXC,irel,xccore,zcore,norb,text,
c     +     nr,aa,bb,r,nql,delql,nqnl,delqnl,numnl,
c     +     itype,cdtyp,0,(l-1),mode,vtemp)
c        endif
300     continue

c        goto 50
c
c      add exchange and correlation
c

 50     continue

ccccccccccccccccccccccccccccccccccccccccccccccccc
c here was the functional specification section c
ccccccccccccccccccccccccccccccccccccccccccccccccc


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C NEED TO PUT IN THE TREATMENT OF SPIN POLARIZATION FROM LIBXC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       do i=2,nr
c         so how do we define this line now:
c         rho(i)=(cdd(i)+cdu(i))/4.d0/pi/r(i)**2


c         CAREFUL: this looks clumsy       
          if(nspol==2)then 
             rho(i,1)=(cdd(i))/2.d0/pi/r(i)**2
             rho(i,2)=(cdu(i))/2.d0/pi/r(i)**2
          else
             rho(i,1)=(cdd(i)+cdu(i))/4.d0/pi/r(i)**2
          end if
c         what about the core charge cdc?
       enddo
c      some : added here
       rho(1,:)=rho(2,:)-(rho(3,:)-rho(2,:))*r(2)/(r(3)-r(2)) 
c
cCMK   this should avoid problems with XC-functionals
cCMK   with kinks (BLYP,....)
c      if (iter.lt.30) then
c        mfxcx=0
c        mfxcc=9
c        mgcc=0
c        mgcx=0
c      else if (iter.eq.30) then
c        write(6,*) 'Switching from LDA to the requested functional'
c      endif
c
c     hutter's routine
c       call evxc(nr,r,rho,vxcgrd,excgrd)
c     goedecker's routine
       call ggaenergy_15(nspol,nr,rw,rd,rho,enexc,vxcgrd,excgrd)
c                rho and vxcgr are now of dimension  (ng,nspol)
c                                        
      write(17,*)'atom iter,E',iter,enexc
c
c     c.hartwig modified integration
       exc=0.d0
       vxc=0.d0
       do i=1,nr
c     need energy/potential in ryd
c     this section was and is a bit misleading
c     let us keep this style for now
          if(nspol==1)then
             exct = 2.d0*excgrd(i)*rho(i,1)
          else
             exct = 1.d0*excgrd(i)*(rho(i,1)+rho(i,2))
          end if 

          if(nspol==1)then
            vxcd =  2.d0*vxcgrd(i,1)
            vxcu=vxcd
          else  
c ??????????????????????????
            vxcd =  2.d0*vxcgrd(i,1)
            vxcu =  2.d0*vxcgrd(i,2)
          end if
          if(nspol==1)then
            rhodw=rho(i,1)/2.d0
            rhoup=rhodw
          else
c ??????????????????????????
            rhodw=rho(i,1)/1.d0
            rhoup=rho(i,2)/1.d0
          end if
          vod(i) = vod(i) + vxcd
          vou(i) = vou(i) + vxcu
          vxc = vxc + (vxcd*rhodw + vxcu*rhoup) * rw(i)
          exc = exc + exct * rw(i)
          write(18,*)vxc, vxcd
       enddo
       etot(4) = ehart
       etot(5) = vxc
       etot(7) = exc
       return
       end
c
