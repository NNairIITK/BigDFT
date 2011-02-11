
      subroutine gatom(nspol,energ,verbose,
     :     noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,nsmx,
     :     occup,aeval,chrg,dhrg,ehrg,res,wght,wfnode,psir0,
     :     rcov,rprb,rcore,zcore,znuc,zion,rloc,gpot,r_l,hsep,
     :     vh,xp,rmt,rmtg,ud,nint,ng,ngmx,psi,igrad,
     :     rr,rw,rd,ntime,itertot,etotal)


c          updated version with spin polarization and libXC support
c          ________________________________________________________


c          Some notes about input variables:

c          nspol are the spin channels for the XC function.
c          nspin are spin components of orbitals with l > 0.
c          nspin may differ from nspol in the relativistic case.
c          the XC treatment and gaussian basis are initialized in the calling routine.
c          energ requests a total energy calculation.
c          verbose requests detailed output after the wfn is converged.
c          Above two logical variables have been moved from a common block.

      
c optimized version of gatom with level 2 BLAS calls
      implicit real*8 (a-h,o-z)
      dimension occup(noccmx,lmx,nsmx),
     :     aeval(noccmx,lmx,nsmx),chrg(noccmx,lmx,nsmx),
     :     dhrg(noccmx,lmx,nsmx),ehrg(noccmx,lmx,nsmx),
     :     res(noccmx,lmx,nsmx),wght(noccmx,lmx,nsmx,8),
     :     wfnode(noccmx,lmx,nsmx,3),
     :     gpot(4),r_l(lmx),hsep(6,lpmx,nsmx),
     :     vh(((ng+1)*(ng+2))/2,lcx+1,((ng+1)*(ng+2))/2,lmax+1),
     :     xp(0:ng), rmt(nint,((ng+1)*(ng+2))/2,lmax+1),
     :     rmtg(nint,((ng+1)*(ng+2))/2,lmax+1),
     :     ud(nint,((ng+1)*(ng+2))/2,lcx+1),
     :     psi(0:ngmx,noccmx,lmx,nsmx)
c     local arrays
c     notice nspin and nspol differ
      dimension  hh(0:ng,0:ng,lmax+1,nspin),ss(0:ng,0:ng,lmax+1),
     :     hht(0:ng,0:ng),sst(0:ng,0:ng),hhsc(((ng+1)*(ng+2))/2,lmax+1),
c          new work array for XC, used only in the polarized case
c    :     hhxc(((ng+1)*(ng+2))/2,lmax+1,2*(nspol-1)),
     :     hhxc(((ng+1)*(ng+2))/2,lmax+1,nspin),
     :     eval(0:ng),evec(0:ng,0:ng),pp1(0:ng,lpx+1),
     :     pp2(0:ng,lpx+1),pp3(0:ng,lpx+1),potgrd(nint),
     :     rho(((ng+1)*(ng+2))/2,lmax+1,nspol),
     :     rhoold(((ng+1)*(ng+2))/2,lmax+1,nspol),excgrd(nint),
     :     vxcgrd(nint,nspol),rr(nint),rw(nint),rd(nint),pexgrd(nint),
     :     ppr1(nint,lmax+1),ppr2(nint,lmax+1),ppr3(nint,lmax+1),
     :     aux1(nint),aux2(nint,0:ng,lmax+1),
     :     expxpr(0:ng,nint), tts(nspol)

      dimension rhogrd(nint,nspol),drhogrd(nint,nspol),
     :                             rhocore(nint,nspol)



      dimension y1(nint),y2(nint),y3(nint),
     :     rlist(0:nint),drlist(0:nint),ddrlist(0:nint)

       character*10 is(2)
      real*8 gamma
      logical energ,igrad, verbose, pol

      external gamma,wave,dwave,ddwave
      save nscf,nscfo,delta,odelta
      fourpi = 16.d0*atan(1.d0)
c      print*,'entered gatom, nspin and nspol are',nspin,nspol
c      is(1)= 'so=+0.5'
c      is(2)= 'so=-0.5'
c      print*,'rcov,rprb,zion,rloc,gpot'
c      print*,rcov,rprb,zion,rloc
c      print*,gpot(1)
c      print*,gpot(2)
c      print*,gpot(3)
c      print*,gpot(4)
c      print*,'lpx,lmax,lcx,noccmax,nspin'
c      print*,lpx,lmax,lcx,noccmax,nspin
c      print*,'ng,nint:',ng,nint
c            do l=0,lpx
c               write(6,*) 'l=',l
c               write(6,'(f7.3,t8,6e11.3,t76,a)') r_l(l+1),
c     :              (hsep(i,l+1,1),i=1,6),'r_l(),hsep(), '//is(1)
c               if (l.gt.0 .and. nspin.eq.2)
c     :              write(6,'(t8,6e11.3,t76,a)')
c     :              (hsep(i,l+1,2),i=1,6),'       hsep(), '//is(2)
c            enddo
c       print*,'xp',xp
cc----------------------------------------------end of debug output




      if (ntime .eq. 0) then
c     c. hartwig    modified scf mixing
c     initialize variables for first iteration rmix
         delta= 0.25d0
         odelta = 0.0d0
         nscf = 0
         nscfo=2**10
c     initial guess for wavefunction
         do l=0,lmax
            do iocc=1,noccmax
               do ispin=1,nspin
                  do i=0,ng
                     psi(i,iocc,l+1,ispin)=1.d-8
                  enddo
               enddo
            enddo
         enddo
      endif
      ntime=ntime+1

c      just to be sure
       rhocore=0d0


       if(zcore>0d0)then
c      ***************************
c      *Nonlinear Core Correction*
c      ***************************

c      set up a simplest model charge distrubution
c      to take into account a frozen core charge for
c      the evaluation of exchange correlation terms. 
c      For now, use EITHER one  gaussian of the form

c      rhocore(r) = zcore *  exp(-0.5*(r/rcore)**2)
c                            /(sqrt(2pi)*rcore)**3

c      OR a simple exponential form:

c      rhocore(r) = zcore *  exp(-2*r/rcore)*(1+r/rcore)
c                                /(4pi rcore**3)

c      Notice that zcore may be a fraction of the 
c      true core charge. Both zcore and rcore are
c      fitting parameters.


Cc----------------------------------------------------------------
c          A few TEST lines to check the impact of the 
Cc         exact core charge, linear interpolation
cc         from a dumpfile created by atom.f
C
C          open(22,file='ae.core.dens.plt')
C          read(22,*)
C          read(22,*)
C          rless=0d0
C          ibig=0
C          rhocore=0d0
C          tt=1d0/(16d0*atan(1d0))
C          do k=1,nint
C             r=rr(k)
Cc            find two grid points with rr(k) in between
C             do
C                 ibig=ibig+1
C                 read(22,*) rbig, rhobig
Cc                write(127,*)k,ibig,rbig,rhobig
C                 if(rbig>r)exit
C                 rless=rbig
C                 psiless=psibig
C             end do
C             rhocore(k,1)=(znuc-zion)*
C     :         ((rbig-r)*rholess+(r-rless)*rhobig)
C     :                   /(rbig-rless)
Cc            divide by 4pi r**2
C             rhocore(k,1)=rhocore(k,1)*tt/(r**2)  
C            
C             if( psibig**2+psiless**2<1e-20 ) exit
C          end do
C          close(22)
Cc----------------------------------------------------------------


c         Set up the core charge using rcore and zcore

c         Let us use the sign of rcore to choose from the
c         two different analytic forms we try here.
          
          if(rcore>0d0)then
c            use a single gaussian with sigma=rcore and charge zcore

c            normalization factor of one 3d gaussian
             tt=zcore* (sqrt(0.125d0/atan(1d0))/rcore)**3
c              =zcore   *  1/sqrt(2pi)**3  *  1/rcore**3
             do k= 1,nint
                rhocore(k,1)=tt*exp(-0.5d0*(rr(k)/rcore)**2)
             end do
          else
c               if rcore is negative use an
c               experimental analytic form:
c                         exponential decay
c                         smooth at r=0
             tt =  zcore/(16d0*atan(1d0)*abs(rcore)**3)
c               =  zcore/(4pi rcore**3) 
             do k= 1,nint
                r=rr(k)
                rhocore(k,1)=tt*(1d0+2d0*r/rcore)*exp(2d0*r/rcore)
             end do
          end if

Cc         a plot for debuging
C           tt=0d0
C           write(33,*) '# rcore zcore',rcore,zcore
C           do k=1,nint
C              tt=tt+rhocore(k,1)*rw(k)
C              write(33,*) rr(k),rhocore(k,1),tt
C           end do
          if(nspol==2)then
c           split the charge equally among the two channels
c           even though the core charge should not be polarized,
c           it is stored in two spin channels for now.
            rhocore(:,1)=rhocore(:,1)*.5d0
            rhocore(:,2)=rhocore(:,1)
          end if
        end if

c       done setting up frozen core charge for nonlinear core correction
c***********************************************************************


c set up all quantities that are not density dependent
c
      do l=lcx+1,lmax
c        no spin down s orbitals in the r case: nspin=2, nspol=1, 2l+1=1
         do ispin=1,max(min(2*l+1,nspin),nspol)
            if (occup(1,l+1,ispin).gt.1.d-10) stop 'lcx too small'
         enddo
      enddo

c projectors
      do l=0,lpx
c     gaussians
         gml1=sqrt( gamma(l+1.5d0) / (2.d0*r_l(l+1)**(2*l+3)) )
         gml2=sqrt( gamma(l+3.5d0) / (2.d0*r_l(l+1)**(2*l+7)) )
     1        /(l+2.5d0)
         gml3=sqrt( gamma(l+5.5d0) / (2.d0*r_l(l+1)**(2*l+11)) )
     1        /((l+3.5d0)*(l+4.5d0))
         tt=1.d0/(2.d0*r_l(l+1)**2)
         do i=0,ng
            ttt=1.d0/(xp(i)+tt)
            pp1(i,l+1)=gml1*(sqrt(ttt)**(2*l+3))
            pp2(i,l+1)=gml2*ttt*(sqrt(ttt)**(2*l+3))
            pp3(i,l+1)=gml3*ttt**2*(sqrt(ttt)**(2*l+3))
         enddo
c     radial grid
        rnrm1=1.d0/sqrt(.5d0*gamma(l+1.5d0)*r_l(l+1)**(2*l+3))
        rnrm2=1.d0/sqrt(.5d0*gamma(l+3.5d0)*r_l(l+1)**(2*l+7))
        rnrm3=1.d0/sqrt(.5d0*gamma(l+5.5d0)*r_l(l+1)**(2*l+11))
        do k=1,nint
           r=rr(k)
           ppr1(k,l+1)=rnrm1*r**l*exp(-.5d0*(r/r_l(l+1))**2)
           ppr2(k,l+1)=rnrm2*r**(l+2)*exp(-.5d0*(r/r_l(l+1))**2)
           ppr3(k,l+1)=rnrm3*r**(l+4)*exp(-.5d0*(r/r_l(l+1))**2)
        enddo
      enddo

c   external potential on grid
      do k=1,nint
         r=rr(k)
         pexgrd(k)=.5d0*(r/rprb**2)**2-zion*Derf(r/(sqrt(2.d0)*rloc))/r
     1        + exp(-.5d0*(r/rloc)**2)*
     1        ( gpot(1) + gpot(2)*(r/rloc)**2 + gpot(3)*(r/rloc)**4 +
     1        gpot(4)*(r/rloc)**6 )
      enddo

c     store exp(-xp(i)*r**2) in expxpr()
      do k=1,nint
         r=rr(k)
         do i=0,ng
            expxpr(i,k)= exp(-xp(i)*r**2)
         enddo
      enddo

c     auxillary grids for resid:
      fourpi = 16.d0*atan(1.d0)
      do k=1,nint
         r=rr(k)
         aux1(k)=fourpi/rw(k)
         do ll=0,lmax
            do i=0,ng
               aux2(k,i,ll+1)=(xp(i)*(3.d0+2.d0
     :              *ll-2.d0*xp(i)*r**2)*expxpr(i,k))
            enddo
         enddo
      enddo
c
c set up charge independent part of hamiltonian
c
      do l=0,lmax
c        no spin down s orbitals in the r case: nspin=2, nspol=1, 2l+1=1
         do ispin=1,max(min(2*l+1,nspin),nspol)
            gml=0.5d0*gamma(l+0.5d0)
c

c     lower triangles only
c
            do j=0,ng
               do i=j,ng
                  hhij=0.0d0
                  d=xp(i)+xp(j)
                  sxp=1.d0/d
                  const=gml*sqrt(sxp)**(2*l+1)
c     overlap
                  ss(i,j,l+1)=const*sxp*(l+.5d0)
c     kinetic energy
                  hhij=.5d0*const*sxp**2*(3.d0*xp(i)*xp(j)+
     1                 l*(6.d0*xp(i)*xp(j)-xp(i)**2-xp(j)**2) -
     1                 l**2*(xp(i)-xp(j))**2  )+ .5d0*l*(l+1.d0)*const
c     potential energy from parabolic potential
                  hhij=hhij+
     1                 .5d0*const*sxp**2*(l+.5d0)*(l+1.5d0)/rprb**4
c     hartree potential from ionic core charge
                  tt=sqrt(1.d0+2.d0*rloc**2*d)
                  if (l.eq.0) then
                     hhij=hhij-zion/(2.d0*d*tt)
                  else if (l.eq.1) then
                     hhij=hhij-zion*
     &                    (1.d0 + 3.d0*rloc**2*d)/(2.d0*d**2*tt**3)
                  else if (l.eq.2) then
                     hhij=hhij-zion*
     &                    (2.d0+10.d0*rloc**2*d+15.d0*rloc**4*d**2)
     &                    /(2.d0*d**3*tt**5)
                  else if (l.eq.3) then
                     hhij=hhij-zion*3.d0*
     &                    (2.d0+14.d0*rloc**2*d+35.d0*rloc**4*d**2
     &                    +35.d0*rloc**6*d**3)/(2.d0*d**4*tt**7)
                  else
                     stop 'l too big'
                  endif
c     potential from repulsive gauss potential
                  tt=rloc**2/(.5d0+d*rloc**2)

                  pw1=1.5d0+dble(l)
                  pw2=2.5d0+dble(l)
                  pw3=3.5d0+dble(l)
                  pw4=4.5d0+dble(l)
                  hhij=hhij+gpot(1)*.5d0*gamma(pw1)*tt**pw1
     1                 + (gpot(2)/rloc**2)*.5d0*gamma(pw2)*tt**pw2
     1                 + (gpot(3)/rloc**4)*.5d0*gamma(pw3)*tt**pw3
     1                 + (gpot(4)/rloc**6)*.5d0*gamma(pw4)*tt**pw4
c     separabel terms
                  if (l.le.lpx) then
                     hhij = hhij
     1                    + pp1(i,l+1)*hsep(1,l+1,ispin)*pp1(j,l+1)
     1                    + pp1(i,l+1)*hsep(2,l+1,ispin)*pp2(j,l+1)
     1                    + pp2(i,l+1)*hsep(2,l+1,ispin)*pp1(j,l+1)
     1                    + pp2(i,l+1)*hsep(3,l+1,ispin)*pp2(j,l+1)
     1                    + pp1(i,l+1)*hsep(4,l+1,ispin)*pp3(j,l+1)
     1                    + pp3(i,l+1)*hsep(4,l+1,ispin)*pp1(j,l+1)
     1                    + pp2(i,l+1)*hsep(5,l+1,ispin)*pp3(j,l+1)
     1                    + pp3(i,l+1)*hsep(5,l+1,ispin)*pp2(j,l+1)
     1                    + pp3(i,l+1)*hsep(6,l+1,ispin)*pp3(j,l+1)
                  endif
                  hh(i,j,l+1,ispin)=hhij
               enddo
            enddo
         enddo
      enddo
c     hhxc is kept constant at zero in the unpolarized case
      hhxc=0d0

c
c finished setup of hh()
c


c initial charge and ev
      do l=0,lmax
         ij=0
         do j=0,ng
            do i=j,ng
               ij=ij+1
               rho(ij,l+1,:)=0.d0
            enddo
         enddo
      enddo
      evsum=1.d30



cccccccccccccccccccccccccccccccc
c     begin the SCF cycles     c
cccccccccccccccccccccccccccccccc


      do it=1,200
         evsumold=evsum
         evsum=0.d0
c
c     coefficients of charge density
c
c     the index for spin polarization shall be
         isp=1
c     in the unpolarized case, then
c     rho(i,l,1)  holds the total charge,
c     while in the polarized case nspol=2
c     rho(i,l,isp) are the two spin channels

         do l=0,lmax
            ij=0
            do j=0,ng
               do i=j,ng
                  ij=ij+1
                  rhoold(ij,l+1,:)=rho(ij,l+1,:)
                  rho(ij,l+1,:)=0.d0
               enddo
            enddo
         enddo
         do l=0,lmax
c           no spin down s orbitals in the r case: nspin=2, nspol=1, 2l+1=1
            do ispin=1,max(min(2*l+1,nspin),nspol)
c              isp is always one in the unpolarized case
c              it determines which spin index of rho is addressed
               isp=min(ispin,nspol)
               do iocc=1,noccmax
                  if (occup(iocc,l+1,ispin).gt.1.d-10) then
                     ij=0
                     do j=0,ng
                        i=j
                        ij=ij+1
                        rho(ij,l+1,isp)=rho(ij,l+1,isp) +
     :                       psi(i,iocc,l+1,ispin)
     :                       *psi(j,iocc,l+1,ispin)
     :                       *occup(iocc,l+1,ispin)
                        do i=j+1,ng
                           ij=ij+1
                           rho(ij,l+1,isp)=rho(ij,l+1,isp) +
     1                          psi(i,iocc,l+1,ispin)
     1                          *psi(j,iocc,l+1,ispin)
     1                          *(2.d0*occup(iocc,l+1,ispin))
                        enddo
                     enddo
                  endif
               enddo
            enddo
         enddo

c determine optimal value for rmix
c     for minimum number of scf iterations
         if ( mod(ntime,20).eq.0 .and. it.eq.2) then
            tttt = delta
            if (nscf .lt. nscfo) then
               if (delta .gt. odelta) then
                  delta = delta + 0.05d0
               else
                  delta = delta - 0.05d0
               endif
            else
               if (delta .gt. odelta) then
                  delta = delta - 0.05d0
               else
                  delta = delta + 0.05d0
               endif
            endif
            delta = max(0.1d0,delta)
            delta = min(0.9d0,delta)
            odelta = tttt
            nscfo = nscf
            nscf = 0
         endif
c     F90 intrinsic
         call random_number(rmix)
         rmix = delta + (.5d0-delta/2.d0)*rmix
c     Intel (ifc)
c        rmix = delta + (.5d0-delta/2.d0)*dble(rand(0.0d0))
c     IBM/DEC/PGI
c        rmix = delta + (.5d0-delta/2.d0)*dble(rand())
c     CRAY
c        rmix = delta + (.5d0-delta/2.d0)*ranf()
c     rmix = delta
         if (it.eq.1) rmix=1.d0
         do l=0,lmax
            ij=0
            do j=0,ng
               do i=j,ng
                  ij=ij+1
c                 the : is over nspol components, 1 or 2
                  tts=rmix*rho(ij,l+1,:) + (1.d0-rmix)*rhoold(ij,l+1,:)
                  rho(ij,l+1,:)=tts
               enddo
            enddo
         enddo
c
c     calc. gradient only if xc-func. with gradient-corrections
c     rho on grid ij=1,nint:
c     rhogrd(k) =+ rmt(k,i,j,l+1)*rho(i,j,l+1)/(4*pi)
c     corresponding gradient:
c     drhogrd(k)=+ rmtg(k,i,j,l+1)*rho(i,j,l+1)/(4*pi)





         tt=1.d0/(16.d0*atan(1.d0))
         call DGEMV('N',nint,((ng+1)*(ng+2))/2*(lcx+1),
     &        tt,rmt,nint,rho(:,:,1),1,0.d0,rhogrd(:,1),1)
c        try yo keep it simple. Same procedure for down gradient
         if(nspol==2)then
         call DGEMV('N',nint,((ng+1)*(ng+2))/2*(lcx+1),
     &        tt,rmt,nint,rho(:,:,2),1,0.d0,rhogrd(:,2),1)
c        there was a factor of two to take into account in atom.f
         rhogrd=rhogrd*1.0d0
         end if
c     for ggaenergy15, we don't need the gradient, as that driver
c     provides the derivative by finite differences on the radial grid.
c     Therefore, calculation of drhogrid is commented out and not
c     generalized to the spin polarized case.

c        if(igrad) call DGEMV('N',nint,((ng+1)*(ng+2))/2*(lcx+1),
c    &           tt,rmtg,nint,rho,1,0.d0,drhogrd,1)



c      The Nonlinear Core Correction (NLCC)
c      is added to the charge density
c      prior to calling the XC drivers
       if(zcore>0d0) rhogrd = rhogrd + rhocore



c     hutter
c      call evxc(nint,rr,rhogrd,drhogrd,vxcgrd,excgrd)
c     goedecker
c     libXC wrapper
c     call ggaenergy_15(nspol,nint,rw,rd,rhogrd,enexc,vxcgrd,excgrd)
c     call ggaenergy_15(nspol,nint,rr,rw,rd,rhogrd,enexc,vxcgrd,excgrd)
      call driveXC(nspol,nint,rr,rw,rd,rhogrd,enexc,vxcgrd,excgrd)
c        multiply with dr*r^2 to speed up calculation of matrix elements
c        open(11,file='rhogrd')
         do k=1,nint
            vxcgrd(k,:)=vxcgrd(k,:)*rw(k)/fourpi
         enddo
c        close(11)

c      IMPORTANT: since the real space representation
c      rhogrd will be used again only after recomputing
c      it from the gaussian representation rho,
c      there is no need to subtract the core charge here.
c      EXCEPTION: rhogrd will be passed to etot, which
c      needs the core charge to calculate energy terms
c      from XC and then subtracts rhocore from rhogrd
c      before computing the other energy functionals.

c
c     charge dependent part of hamiltonian
c     hartree potential from valence charge distribution

c     if nspol=2, the up and down potential differs due to Vxc 

C     do 4982,lp=0,lcx
C     do 4982,jp=0,ng
C     do 4982,ip=0,ng
C     hhsc(i,j,l+1) =+ vh(ip,jp,lp+1,i,j,l+1)*rho(ip,jp,lp+1)

         call DGEMV('T',(lcx+1)*((ng+1)*(ng+2))/2,
     &        (lmax+1)*((ng+1)*(ng+2))/2,1.d0,
     &        vh,(lcx+1)*((ng+1)*(ng+2))/2,rho(:,:,1),1,0.d0,hhsc,1)

c     if nspol=2, add the hartree potential from the 2nd charge channel
c     Note: It seems easier to add the charges first and then do one dgemv.
         if(nspol==2) 
     &   call DGEMV('T',(lcx+1)*((ng+1)*(ng+2))/2,
     &        (lmax+1)*((ng+1)*(ng+2))/2,1.d0,
     &        vh,(lcx+1)*((ng+1)*(ng+2))/2,rho(:,:,2),1,1.d0,hhsc,1)
c                                                  ^    ^ 

c     potential from XC libraries 

C     do 8049,k=1,nint
C     8049 hhsc(i,j,l+1) =+ vxcgrd(k)*rmt(k,i,j,l+1)

c     MODIFICATION: if spin polarized, add this term to hhxc, not hhsc.
c                   hxc is a spin polarized matrix only for that purpose.
         hhxc=0d0
         if(nspol==1)then
              call DGEMV('T',nint,(lmax+1)*((ng+1)*(ng+2))/2,1.0d0,
     &        rmt,nint,vxcgrd(:,1),1,1.d0,hhsc,1)
         else
              call DGEMV('T',nint,(lmax+1)*((ng+1)*(ng+2))/2,1.0d0,
     &        rmt,nint,vxcgrd(:,1),1,0.d0,hhxc(:,:,1),1)
              call DGEMV('T',nint,(lmax+1)*((ng+1)*(ng+2))/2,1.0d0,
     &        rmt,nint,vxcgrd(:,2),1,1.d0,hhxc(:,:,2),1)

c        spin polarized XC term end
         end if


c     DIAGONALIZE
         do l=0,lmax
            do ispin=1,max(min(2*l+1,nspin), nspol)
C     LAPACK
               ij=0
               do j=0,ng
                  do i=j,ng
                     ij=ij+1
                     hht(i,j)=hh(i,j,l+1,ispin)+hhsc(ij,l+1)
     :                       +hhxc(ij,l+1,ispin) 
                     sst(i,j)=ss(i,j,l+1)
                  enddo
               enddo
c     IBM/DEC
               call DSYGV(1,'V','L',ng+1,hht,ng+1,sst,ng+1,
     1              eval,evec,(ng+1)**2,info)
c     the routine DSYGV is also included in sub_lapack.f
c     CRAY:
c     call SSYGV(1,'V','L',ng+1,hht,ng+1,sst,ng+1,
c     1       eval,evec,(ng+1)**2,info)
               if (info.ne.0) write(6,*) 'LAPACK',info
               do iocc=0,noccmax-1
                  do i=0,ng
                     evec(i,iocc)=hht(i,iocc)
                  enddo
               enddo
C     end LAPACK
               do iocc=1,noccmax
c                 write(6,*)'DEBUG: E(iocc,ispin,l,it)',
c    :                       eval(iocc-1),iocc-1,ispin,l,it
                  evsum=evsum+eval(iocc-1)
                  aeval(iocc,l+1,ispin)=eval(iocc-1)
                  do i=0,ng
                     psi(i,iocc,l+1,ispin)=evec(i,iocc-1)
                  enddo
               enddo
c     write(6,*) 'eval',l
c     55         format(5(e14.7))
c     write(6,55) eval
c     write(6,*) 'evec',l
c     do i=0,ng
c     33            format(10(e9.2))
c     write(6,33) (evec(i,iocc),iocc=0,noccmax-1)
c     enddo
           enddo
         enddo
         tt=abs(evsum-evsumold)
c        write(6,*)'DEBUG: residue=',tt,it
         if (tt.lt.1.d-8) goto 3000
      enddo
      write(6,*) 'WARNING: NO SC CONVERGENCE',tt
 3000 continue

cccccccccccccccccccccccccccccccc
c    end of the SCF cycles     c
cccccccccccccccccccccccccccccccc

c     write(*,*)'DEBUG: KS eigenvalues',aeval
      itertot=itertot+it
      nscf = nscf +it
      call resid(nspol,
     :     noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,nsmx,
     :     aeval,res,
     :     hsep,
     :     ud,nint,ng,ngmx,psi,rho,pp1,pp2,pp3,
     :     potgrd,pexgrd,vxcgrd,rr,rw,ppr1,ppr2,ppr3,aux1,aux2,
     :     expxpr)
c     etot evaluates Ehartree using rhogrd,
c     so get rid of the frozen core charge first
c     !!!!!!!!!!! dangerous: we want the rhocharge for the integration of Exc and Vxc!
c     what about this: Keep it for now here... ehart will have an offset that SHOULD
c     not depend on the occupation numbers and excitation energies will finally work.
c     if (energ.and. zcore>0d0) rhogrd=rhogrd-rhocore
      if (energ) call etot(verbose,nspol,
     :     noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,nsmx,
     :     aeval,
     :     rprb,zion,rloc,gpot,r_l,hsep,
     :     xp,ud,nint,ng,ngmx,psi,rho,pp1,pp2,pp3,
     :     vxcgrd,excgrd,rhogrd,rhocore,occup,rr,rw,
     :     expxpr,etotal)
c
c     charge up to radius rcov or infinity
c

c     MODIFICATION: Can one integrate up to a different rcov for semicore states?
c     problem: loop over l is implicit here, one can not tell easily where nl<occmax(l)

      if (lmax.gt.3) stop 'cannot calculate chrg'
      do l=0,lmax
c        no spin down s orbitals in the r case: nspin=2, nspol=1, 2l+1=1
         do ispin=1,max(min(2*l+1,nspin), nspol)
            do iocc=1,noccmax
               chrg(iocc,l+1,ispin)=0.d0
               dhrg(iocc,l+1,ispin)=0.d0
               ehrg(iocc,l+1,ispin)=0.d0
            enddo
         enddo
      enddo
      do ispin=1,max(min(2*l+1,nspin),nspol)
c     here, l=lmax+1, so do ispin=1,2 if lmax>0 and nspin=2
         do 3711,iocc=1,noccmax
c        loop over all nl(l)
            do 3762,j=0,ng
               do 3762,i=0,ng
                  d=xp(i)+xp(j)
                  sd=sqrt(d)
                  terf=Derf(sd*rcov)
                  texp=exp(-d*rcov**2)
                  tt0=0.4431134627263791d0*terf/sd**3-0.5d0*rcov*texp/d
                  tt1=0.6646701940895686d0*terf/sd**5 +
     &                 (-0.75d0*rcov*texp - 0.5d0*d*rcov**3*texp)/d**2
                  chrg(iocc,1,ispin)=chrg(iocc,1,ispin) +
     1                 psi(i,iocc,1,ispin)*psi(j,iocc,1,ispin)*tt0
c     integrate up to rcov
c               dhrg(iocc,1,ispin)=dhrg(iocc,1,ispin) +
c     1              psi(i,iocc,1,ispin)*psi(j,iocc,1,ispin)*tt1
c     integrate up to inf
                  dhrg(iocc,1,ispin)=dhrg(iocc,1,ispin) +
     :                 psi(i,iocc,1,ispin)*psi(j,iocc,1,ispin)
     :                 *0.6646701940895686d0/sd**5
                  ehrg(iocc,1,ispin)=ehrg(iocc,1,ispin) +
     :                 psi(i,iocc,1,ispin)*psi(j,iocc,1,ispin)
     :                 *1.66167548522392d0/sd**7
                  if (lmax.eq.0) goto 3762

                  tt2=1.661675485223921d0*terf/sd**7 +
     &                 (-1.875d0*rcov*texp-1.25d0*d*rcov**3*texp-
     &                 0.5d0*d**2*rcov**5*texp)/d**3
                  chrg(iocc,2,ispin)=chrg(iocc,2,ispin) +
     1                 psi(i,iocc,2,ispin)*psi(j,iocc,2,ispin)*tt1
c     integrate up to rcov
c               dhrg(iocc,2,ispin)=dhrg(iocc,2,ispin) +
c     1              psi(i,iocc,2,ispin)*psi(j,iocc,2,ispin)*tt2
c     integrate up to inf
                  dhrg(iocc,2,ispin)=dhrg(iocc,2,ispin) +
     1                 psi(i,iocc,2,ispin)*psi(j,iocc,2,ispin)
     2                 *1.661675485223921d0/sd**7
                  ehrg(iocc,2,ispin)=ehrg(iocc,2,ispin) +
     1                 psi(i,iocc,2,ispin)*psi(j,iocc,2,ispin)
     2                 *5.815864198283725d0/sd**9
                  if (lmax.eq.1) goto 3762

                  tt3=5.815864198283725d0*terf/sd**9 +
     &                 (-6.5625d0*rcov*texp-4.375d0*d*rcov**3*texp-
     &                 1.75d0*d**2*rcov**5*texp -
     &                 0.5d0*d**3*rcov**7*texp)/d**4
                  chrg(iocc,3,ispin)=chrg(iocc,3,ispin) +
     1                 psi(i,iocc,3,ispin)*psi(j,iocc,3,ispin)*tt2
c     integrate up to rcov
c               dhrg(iocc,3,ispin)=dhrg(iocc,3,ispin) +
c     1              psi(i,iocc,3,ispin)*psi(j,iocc,3,ispin)*tt3
c     integrate up to inf
                  dhrg(iocc,3,ispin)=dhrg(iocc,3,ispin) +
     1                 psi(i,iocc,3,ispin)*psi(j,iocc,3,ispin)
     2                 *5.815864198283725d0/sd**9
                  ehrg(iocc,3,ispin)=ehrg(iocc,3,ispin) +
     1                 psi(i,iocc,3,ispin)*psi(j,iocc,3,ispin)
     2                 *26.17138889227676d0/sd**11

                  if (lmax.eq.2) goto 3762

                  chrg(iocc,4,ispin)=chrg(iocc,4,ispin) +
     1                 psi(i,iocc,4,ispin)*psi(j,iocc,4,ispin)*tt3
c     integrate up to rcov
c                  tt4=26.17138889227676d0*terf/sd**11+(-29.53125d0*
c     :                 rcov*texp-19.6875d0*d*rcov**3*texp-7.875d0*d**2
c     :                 *rcov**5*texp-2.25d0*d**3*rcov**7*texp-
c     &                 0.5d0*d**4*rcov**9*texp)/d**5
c               dhrg(iocc,4,ispin)=dhrg(iocc,4,ispin) +
c     1              psi(i,iocc,4,ispin)*psi(j,iocc,4,ispin)*tt4
c     integrate up to inf
                  dhrg(iocc,4,ispin)=dhrg(iocc,4,ispin) +
     1                 psi(i,iocc,4,ispin)*psi(j,iocc,4,ispin)
     2                 *26.17138889227676d0/sd**11
                  ehrg(iocc,4,ispin)=dhrg(iocc,4,ispin) +
     1                 psi(i,iocc,4,ispin)*psi(j,iocc,4,ispin)
     2                 *143.9426389075222d0/sd**13

 3762       continue
 3711    continue
      enddo

CCc     NEW loop: Recalculate chrg(l,nl,s) for each 'semicore'
CCc               using a different value than rcov. Condition:
CCc               higher occupied orbital with same l exists. 
CCc               For now, pick rsemi = rcov/2; later: read from psp.par
CC      rsemi=0.5d0*rcov
CC      do l=0,lmax-1
CC         do ispin=1,min(2*l+1,nspin)
CC            do iocc=1,noccmax
CCc              check if the next higher orbital of this l channel is occupied
CC               if(occup(iocc,l+2,nspin)<1d-4) exit
CCc              do the charge integral up to rsemi
CC               do j=0,ng
CC                  do i=0,ng
CC                     d=xp(i)+xp(j)
CC                     sd=sqrt(d)
CC                     terf=Derf(sd*rsemi)
CC                     texp=exp(-d*rsemi**2)
CC                     if(l==0) tt=0.4431134627263791d0*terf
CC     :                           /sd**3-0.5d0*rsemi*texp/d
CC                     if(l==1) tt=0.6646701940895686d0*terf/sd**5 +
CC     &                   (-0.75d0*rcov*texp - 0.5d0*d*rcov**3*texp)/d**2
CC                     chrg(iocc,l+1,ispin)=chrg(iocc,l+1,ispin) +
CC     1                    psi(i,iocc,l+1,ispin)*psi(j,iocc,l+1,ispin)*tt
CC                  end do
CC               end do
CC            end do
CC         end do
CC      end do
CC
CC               


c
c     value at origin
c
      psir0=0.d0
      do i=0,ng
         psir0=psir0+psi(i,1,1,1)
      enddo
      psir0=psir0**2

c     node locations of psi*r
c     n-1 nodes allowed!
c     search only for nodes if the corresponding weights are <> zero
c     to avoid bumpy wavefunctions: no node of the first derivative of
c     the pseudowavefunction*r  and only one node
c     of the second derivative  up to the rmax of the lowest valence state
c
      tol =1.0d-12
c initialize ALL elements of wfnode to zero, also unused ones. 
c this is only to test whether this helps to get rid of a bug
      wfnode=0d0
      do l=0,lmax
c        no spin down s orbitals in the r case: nspin=2, nspol=1, 2l+1=1
         do ispin=1, min(min(2*l+1,nspin),nspol)
            do nocc=1,noccmax
               if ( (wght(nocc,l+1,ispin,6).ne.0.d0)
     :              .or.  (wght(nocc,l+1,ispin,7).ne.0.d0)
     :              .or.  (wght(nocc,l+1,ispin,8).ne.0.d0) ) then

c       print*,'node search, l, nocc:',l,' ',nocc
                  wfnode(nocc,l+1,ispin,1)=0.d0
                  wfnode(nocc,l+1,ispin,2)=0.d0
                  wfnode(nocc,l+1,ispin,3)=0.d0
                  nnode=0
                  ndnode=0
                  nddnode=0
                  rnode=0.d0
                  rrnode=0.0d0
                  rrdnode=0.d0
                  dnode=0.d0
                  ddnode=0.0d0
c     find outer max of psi, search from ~10 bohr down
                  call detnp(nint,rr,10.0d0,kout)
                  ttrmax=rr(kout)
                  ra=ttrmax
                  ttmax= dabs(wave2(ng,l,psi(0,nocc,l+1,ispin),
     :                 expxpr,ra,kout,nint))
c      print*,'ttmax=',ttmax
                  do k=kout,1, -1
                     ra= rr(k)
                     ttpsi= dabs(wave2(ng,l,psi(0,nocc,l+1,ispin),
     :                 expxpr,ra,k,nint))
                     if ( ttpsi .gt. ttmax
     :                    .and. ttpsi .gt. 1.0d-4 ) then
                        ttmax=ttpsi
                        ttrmax=ra
                     endif
                  if (ttpsi.lt.ttmax .and. ttpsi.gt.1.0d-4) goto 3456
                  enddo
 3456             continue
c     search up to 90% of rmax
                  ttrmax=max(0.90d0*ttrmax,rr(1))
                  call detnp(nint,rr,ttrmax,kout)
                  ttrmax=rr(kout)
c       print*,'search up to ',ttrmax,ttmax
c     calc wavefunction and it's first two derivatives on the grid
c
                  do k=1,kout
                     call wave3(ng,l,xp,psi(0,nocc,l+1,ispin),
     :                    expxpr,rr(k),k,nint,y1(k),y2(k),y3(k))
                  enddo

                     do k = 2,kout
c     nodes of wavefunction
                     if (y1(k)*y1(k-1).lt.0.d0) then
                        nnode = nnode +1
                        x1=rr(k-1)
                        x2=rr(k)
                        rrnode = zbrent(wave,ng,ngmx,l,lmx,xp,psi,
     :                       nocc,noccmx,ispin,nsmx,
     :                       X1,X2,TOL)
                        if (nnode .ge.nocc) then
                           rnode=rnode+rrnode
c                          print*,'found rnode at:',rrnode
                        endif
                        rlist(nnode)=rrnode
                     endif
c     nodes of first derivative
                     if (y2(k)*y2(k-1).lt.0.d0) then
                        ndnode = ndnode +1
                        x1=rr(k-1)
                        x2=rr(k)
                        rrnode = zbrent(dwave,ng,ngmx,l,lmx,xp,psi,
     :                       nocc,noccmx,ispin,nsmx,
     :                       X1,X2,TOL)
                        if (ndnode .ge.nocc) then
                           dnode=dnode+rrnode
c                        print*,'found dnode at:',rrnode
                        endif
                        drlist(ndnode)=rrnode
                     endif
c     second derivative test:
                     if (y3(k)*y3(k-1).lt.0.d0) then
                        nddnode = nddnode + 1
                        x1=rr(k-1)
                        x2=rr(k)
                        rrnode = zbrent(ddwave,ng,ngmx,l,lmx,xp,psi,
     :                       nocc,noccmx,ispin,nsmx,
     :                       X1,X2,TOL)
c     only add the lowest node! (this one shoud dissapear)
                        if (nddnode .ge. nocc +1 ) then
                           ddnode = ddnode + rrdnode
c                          print*,'found ddnode at:',rrnode
                        else
                           rrdnode=rrnode
                        endif
                        ddrlist(nddnode)=rrnode
                     endif
                  enddo

c     print*,'rnode,dnode,ddnode',rnode,dnode,ddnode,nnode

c     new version: use integral of the relevant functions between the nodes
c     not the node-locations!
c     calc. necessary integrals:
                  sum1=0.0d0
                  sum2=0.0d0
                  sum3=0.0d0
c     rnodes:
                  do i=nnode+1-nocc,1,-2
                     aa=Wwav(ng,l,xp,psi(0,nocc,l+1,ispin),rlist(i))
     :                 -Wwav(ng,l,xp,psi(0,nocc,l+1,ispin),rlist(i-1))
c                    WE sometimes get NaNs.... do some test here to trace
c                    if(.not. (aa>0d0.or.aa<=0d0))
c    :               write(6,*)'wwav=NaN for l,s,nl,node',l,ispin,nocc,i
                     sum1 = sum1+aa
                  enddo
c     dnodes
                  do i=ndnode+1-nocc,1,-2
                     aa=wave(ng,l,xp,psi(0,nocc,l+1,ispin),drlist(i))
     :                 -wave(ng,l,xp,psi(0,nocc,l+1,ispin),drlist(i-1))
                     sum2 = sum2+aa
                  enddo
c     ddnodes
                  do i=nddnode+1-nocc,1,-2
                     aa=dwave(ng,l,xp,psi(0,nocc,l+1,ispin),ddrlist(i))
     :                -dwave(ng,l,xp,psi(0,nocc,l+1,ispin),ddrlist(i-1))
                     sum3 = sum3+aa
                  enddo
c     old version for nodes as used in the paper:
c                  wfnode(nocc,l+1,ispin,1)=rnode
c                  wfnode(nocc,l+1,ispin,2)=dnode
c                  wfnode(nocc,l+1,ispin,3)=ddnode
c     new version, using the integrals of the function between the nodes
                  wfnode(nocc,l+1,ispin,1)=sum1
                  wfnode(nocc,l+1,ispin,2)=sum2
                  wfnode(nocc,l+1,ispin,3)=sum3
               endif
            enddo
         enddo
      enddo

c     print*,'leave gatom'

      end

