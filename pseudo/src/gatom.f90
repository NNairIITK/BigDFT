!> @file
!! Generate atomic electronic configuration
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


!> Updated version with spin polarization and libxc support
!! some notes about input variables:
!! nspol are the spin channels for the xc function.
!! nspin are spin components of orbitals with l > 0.
!! nspin may differ from nspol in the relativistic case.
!! the xc treatment and gaussian basis are initialized in the calling routine.
!! energ requests a total energy calculation.
!! verbose requests detailed output after the wfn is converged.
!! above two logical variables have been moved from a common block.
!! @note a note about the r_l2 variable: 
!! this is a feature to test gaussian type projectors with two
!! different length scales r_l and r_l2.
subroutine gatom(energ, verbose)
   
   use pseudovars
   use gatomvars

   implicit none

   !Arguments
   logical, intent(in) :: energ, verbose
   !Local variables
   real(kind=8), parameter :: fourpi = 16.d0*atan(1.d0)
   real(kind=8), dimension(0:ng,0:ng,lmax+1,nspin) :: hh
   real(kind=8), dimension(0:ng,0:ng,lmax+1) :: ss
   real(kind=8), dimension(0:ng,0:ng) :: hht,sst
   real(kind=8), dimension(((ng+1)*(ng+2))/2,lmax+1) :: hhsc
   real(kind=8), dimension(((ng+1)*(ng+2))/2,lmax+1,nspin) :: hhxc !< only for the spin polarized case
   real(kind=8), dimension(((ng+1)*(ng+2))/2,lmax+1,nspol) :: rho, rhoold
   real(kind=8), dimension(0:ng) :: eval
   real(kind=8), dimension(0:ng,0:ng) :: evec
   real(kind=8), dimension(0:ng,lpx) :: pp1,pp2,pp3
   real(kind=8), dimension(nint,lmax+1) :: ppr1,ppr2,ppr3
   real(kind=8), dimension(nint,0:ng,lmax+1) :: aux2
   real(kind=8), dimension(0:ng,nint) :: expxpr
   real(kind=8), dimension(nspol) :: tts
   real(kind=8), dimension(nint,nspol) :: vxcgrd, rhogrd, rhocore
   real(kind=8), dimension(nint) :: potgrd, excgrd, pexgrd, aux1, y1,y2,y3
   real(kind=8), dimension(0:nint) :: rlist,drlist,ddrlist

   real(kind=8) :: aa,const,ddnode,d,dnode,ddrnode,drnode
   real(kind=8) :: enexc,evsum,evsumold,gml,gml1,gml2,gml3,hhij
   real(kind=8) :: pw1,pw2,pw3,pw4,r,r2,ra,rmix,rnode,rnrm1,rnrm2,rnrm3,rrdnode,rrnode,sd
   real(kind=8) :: sum1,sum2,sum3,sxp,terf,texp,tol,tt,tt0,tt1,tt2,tt3,ttmax,ttpsi,ttrmax
   real(kind=8) :: ttt,ttt2,tttt,x1,x2
   integer :: l,iocc,ispin,i,k,ij,info,isp,it,j,kout,ll,lq,nddnode,ndnode,nnode,nocc

   real(kind=8), external :: gamma, wave, dwave, ddwave, wave2,zbrent,wwav

   real(kind=8), save :: delta, odelta
   integer, save :: nscf, nscfo

   if (ntime .eq. 0) then
      !     c. hartwig    modified scf mixing
      !     initialize variables for first iteration rmix
      delta= 0.25d0
      odelta = 0.0d0
      nscf = 0
      nscfo=2**10
      !     initial guess for wavefunction
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
   
   
   !***********************************************************************
   
   ! just to be sure
   rhocore=0d0
   
   if(rcore>0d0)then
      !      ***************************
      !      *nonlinear core correction*
      !      ***************************
      
      !      set up a simplest model charge distrubution
      !      to take into account a frozen core charge for
      !      the evaluation of exchange correlation terms. 
      
      !      to be consistent with the convention in bigdft
      !      let us use a factor of four pi for the nlcc here.
      
      !      careful: in bigdft, the polynomial part is
      !               not scaled by rcore. this results
      !               easily in large coeffs and is not
      !               optimal for fitting. therefore,
      !               we have slightly different conventions.
      !               the corresponding transformation is
      !               done once in pseudo for i/o of nlcc.
      !
      
      do k= 1,nint
         r2=(rr(k)/rcore)**2
         rhocore(k,1) =  &
              exp(-.5d0*r2) /fourpi  *(  &
              gcore(1)  &
              + gcore(2)*r2   &
              + gcore(3)*r2**2   &
              + gcore(4)*r2**3 )
         !           write(17,*)rr(k),rhocore(k,1)
      end do
      
      if(nspol==2)then
         !           split the charge equally among the two channels
         !           even though the core charge should not be polarized,
         !           it is stored in two spin channels for further testing.
         do k= 1,nint
            rhocore(k,1)=rhocore(k,1)*.5d0
            rhocore(k,2)=rhocore(k,1)
         end do
      end if
   end if
   
   !***********************************************************************
   
   
   ! set up all quantities that are not density dependent
   !
   do l=lcx+1,lmax
      !        no spin down s orbitals in the r case: nspin=2, nspol=1, 2l+1=1
      do ispin=1,max(min(2*l+1,nspin),nspol)
         if (occup(1,l+1,ispin).gt.1.d-10) stop 'lcx too small'
      enddo
   enddo
   
   ! projectors
   
   !     r_l2 can be  disabled here:
   !     r_l2 = r_l
   !     it is used instead of r_l for 2nd or 3rd order projectors
   !     which enter the seperable part for any hij other than  h11.
   
   do l=1,lpx
      !        l is a positive index, lq the quantum number
      lq=l-1
      !     gaussians
      gml1=sqrt( gamma(lq+1.5d0) / (2.d0*r_l(l)**(2*lq+3)) )  
      gml2=sqrt( gamma(lq+3.5d0) / (2.d0*r_l2(l)**(2*lq+7)) )  &
           /(lq+2.5d0)
      gml3=sqrt( gamma(lq+5.5d0) / (2.d0*r_l2(l)**(2*lq+11)) )  &
           /((lq+3.5d0)*(lq+4.5d0))
      tt =1.d0/(2.d0*r_l(l)**2)
      tt2=1.d0/(2.d0*r_l2(l)**2)
      do i=0,ng
         ttt =1.d0/(xp(i)+tt)
         ttt2=1.d0/(xp(i)+tt2)
         pp1(i,l)=gml1*(sqrt(ttt)**(2*lq+3))
         pp2(i,l)=gml2*ttt2*(sqrt(ttt2)**(2*lq+3))
         pp3(i,l)=gml3*ttt2**2*(sqrt(ttt2)**(2*lq+3))
      enddo
      !     radial grid
      rnrm1=1.d0/sqrt(.5d0*gamma(lq+1.5d0)*r_l(l)**(2*lq+3))
      rnrm2=1.d0/sqrt(.5d0*gamma(lq+3.5d0)*r_l2(l)**(2*lq+7))
      rnrm3=1.d0/sqrt(.5d0*gamma(lq+5.5d0)*r_l2(l)**(2*lq+11))
      do k=1,nint
         r=rr(k)
         ppr1(k,l)=rnrm1*r**lq    *exp(-.5d0*(r/r_l(l))**2)
         ppr2(k,l)=rnrm2*r**(lq+2)*exp(-.5d0*(r/r_l2(l))**2)
         ppr3(k,l)=rnrm3*r**(lq+4)*exp(-.5d0*(r/r_l2(l))**2)
      enddo
   enddo
   
   
   
   
   !   external potential on grid
   do k=1,nint
      r=rr(k)
      pexgrd(k)=.5d0*(r/rprb**2)**2-zion*derf(r/(sqrt(2.d0)*rloc))/r  &
           + exp(-.5d0*(r/rloc)**2)*  &
           ( gpot(1) + gpot(2)*(r/rloc)**2 + gpot(3)*(r/rloc)**4 +  &
           gpot(4)*(r/rloc)**6 )
   enddo
   
   !     store exp(-xp(i)*r**2) in expxpr()
   do k=1,nint
      r=rr(k)
      do i=0,ng
         expxpr(i,k)= exp(-xp(i)*r**2)
      enddo
   enddo
   
   !     auxillary grids for resid:
   do k=1,nint
      r=rr(k)
      aux1(k)=fourpi/rw(k)
      do ll=0,lmax
         do i=0,ng
            aux2(k,i,ll+1)=(xp(i)*(3.d0+2.d0  &
                 *ll-2.d0*xp(i)*r**2)*expxpr(i,k))
         enddo
      enddo
   enddo
   !
   ! set up charge independent part of hamiltonian
   !
   do l=1,lmax+1
      !        l is a positive index, lq the quantum number
      lq=l-1
      !        no spin down s orbitals in the r case: nspin=2, nspol=1, 2l+1=1
      do ispin=1,max(min(2*lq+1,nspin),nspol)
         gml=0.5d0*gamma(lq+0.5d0)
         !
         
         !     lower triangles only
         !
         do j=0,ng
            do i=j,ng
               hhij=0.0d0
               d=xp(i)+xp(j)
               sxp=1.d0/d
               const=gml*sqrt(sxp)**(2*lq+1)
               !     overlap
               ss(i,j,l)=const*sxp*(lq+.5d0)
               !     kinetic energy
               hhij=.5d0*const*sxp**2*(3.d0*xp(i)*xp(j)+  &
                    lq*(6.d0*xp(i)*xp(j)-xp(i)**2-xp(j)**2) -  &
                    lq**2*(xp(i)-xp(j))**2  )+ .5d0*lq*(lq+1.d0)*const
               !     potential energy from parabolic potential
               hhij=hhij+  &
                    .5d0*const*sxp**2*(lq+.5d0)*(lq+1.5d0)/rprb**4
               !     hartree potential from ionic core charge
               tt=sqrt(1.d0+2.d0*rloc**2*d)
               if (lq.eq.0) then
                  hhij=hhij-zion/(2.d0*d*tt)
               else if (lq.eq.1) then
                  hhij=hhij-zion*  &
                       (1.d0 + 3.d0*rloc**2*d)/(2.d0*d**2*tt**3)
               else if (lq.eq.2) then
                  hhij=hhij-zion*  &
                       (2.d0+10.d0*rloc**2*d+15.d0*rloc**4*d**2)  &
                       /(2.d0*d**3*tt**5)
               else if (lq.eq.3) then
                  hhij=hhij-zion*3.d0*  &
                       (2.d0+14.d0*rloc**2*d+35.d0*rloc**4*d**2  &
                       +35.d0*rloc**6*d**3)/(2.d0*d**4*tt**7)
               else
                  stop 'l too big'
               endif
               !     potential from repulsive gauss potential
               tt=rloc**2/(.5d0+d*rloc**2)
               
               pw1=1.5d0+dble(lq)
               pw2=2.5d0+dble(lq)
               pw3=3.5d0+dble(lq)
               pw4=4.5d0+dble(lq)
               hhij=hhij+gpot(1)*.5d0*gamma(pw1)*tt**pw1  &
                    + (gpot(2)/rloc**2)*.5d0*gamma(pw2)*tt**pw2  &
                    + (gpot(3)/rloc**4)*.5d0*gamma(pw3)*tt**pw3  &
                    + (gpot(4)/rloc**6)*.5d0*gamma(pw4)*tt**pw4
               !     separabel terms
               if (l.le.lpx) then
                  hhij = hhij  &
                       + pp1(i,l)*hsep(1,l,ispin)*pp1(j,l)  &
                       + pp1(i,l)*hsep(2,l,ispin)*pp2(j,l)  &
                       + pp2(i,l)*hsep(2,l,ispin)*pp1(j,l)  &
                       + pp2(i,l)*hsep(3,l,ispin)*pp2(j,l)  &
                       + pp1(i,l)*hsep(4,l,ispin)*pp3(j,l)  &
                       + pp3(i,l)*hsep(4,l,ispin)*pp1(j,l)  &
                       + pp2(i,l)*hsep(5,l,ispin)*pp3(j,l)  &
                       + pp3(i,l)*hsep(5,l,ispin)*pp2(j,l)  &
                       + pp3(i,l)*hsep(6,l,ispin)*pp3(j,l) 
                  
                  
               endif
               hh(i,j,l,ispin)=hhij
            enddo
         enddo
      enddo
   enddo
   !     hhxc is kept constant at zero in the unpolarized case
   hhxc=0d0
   
   !
   ! finished setup of hh()
   !
   
   
   ! initial charge and ev
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
   
   
   
   !ccccccccccccccccccccccccccccccc
   !     begin the scf cycles     c
   !ccccccccccccccccccccccccccccccc
   
   
   do it=1,200
      evsumold=evsum
      evsum=0.d0
      !
      !     coefficients of charge density
      !
      !     the index for spin polarization shall be
      isp=1
      !     in the unpolarized case, then
      !     rho(i,l,1)  holds the total charge,
      !     while in the polarized case nspol=2
      !     rho(i,l,isp) are the two spin channels
      
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
         !           no spin down s orbitals in the r case: nspin=2, nspol=1, 2l+1=1
         do ispin=1,max(min(2*l+1,nspin),nspol)
            !              isp is always one in the unpolarized case
            !              it determines which spin index of rho is addressed
            isp=min(ispin,nspol)
            do iocc=1,noccmax
               if (occup(iocc,l+1,ispin).gt.1.d-10) then
                  ij=0
                  do j=0,ng
                     i=j
                     ij=ij+1
                     rho(ij,l+1,isp)=rho(ij,l+1,isp) +  &
                          psi(i,iocc,l+1,ispin)  &
                          *psi(j,iocc,l+1,ispin)  &
                          *occup(iocc,l+1,ispin)
                     do i=j+1,ng
                        ij=ij+1
                        rho(ij,l+1,isp)=rho(ij,l+1,isp) +  &
                             psi(i,iocc,l+1,ispin)  &
                             *psi(j,iocc,l+1,ispin)  &
                             *(2.d0*occup(iocc,l+1,ispin))
                     enddo
                  enddo
               endif
            enddo
         enddo
      enddo
      
      ! determine optimal value for rmix
      !     for minimum number of scf iterations
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
      !     f90 intrinsic
      call random_number(rmix)
      rmix = delta + (.5d0-delta/2.d0)*rmix
      !     intel (ifc)
      !        rmix = delta + (.5d0-delta/2.d0)*dble(rand(0.0d0))
      !     ibm/dec/pgi
      !        rmix = delta + (.5d0-delta/2.d0)*dble(rand())
      !     cray
      !        rmix = delta + (.5d0-delta/2.d0)*ranf()
      !     rmix = delta
      if (it.eq.1) rmix=1.d0
      do l=0,lmax
         ij=0
         do j=0,ng
            do i=j,ng
               ij=ij+1
               !                 the : is over nspol components, 1 or 2
               tts=rmix*rho(ij,l+1,:) + (1.d0-rmix)*rhoold(ij,l+1,:)
               rho(ij,l+1,:)=tts
            enddo
         enddo
      enddo
      !
      !     calc. gradient only if xc-func. with gradient-corrections
      !     rho on grid ij=1,nint:
      !     rhogrd(k) =+ rmt(k,i,j,l+1)*rho(i,j,l+1)/(4*pi)
      !     corresponding gradient:
      !     drhogrd(k)=+ rmtg(k,i,j,l+1)*rho(i,j,l+1)/(4*pi)
      
      
      
      
      
      tt=1.d0/(16.d0*atan(1.d0))
      call dgemv('n',nint,((ng+1)*(ng+2))/2*(lcx+1),  &
           tt,rmt,nint,rho(:,:,1),1,0.d0,rhogrd(:,1),1)
      !        try yo keep it simple. same procedure for spin down charge
      if(nspol==2)then
         call dgemv('n',nint,((ng+1)*(ng+2))/2*(lcx+1),  &
              tt,rmt,nint,rho(:,:,2),1,0.d0,rhogrd(:,2),1)
      end if
      !     for ggaenergy15, we don't need the gradient, as that driver
      !     provides the derivative by finite differences on the radial grid.
      !     therefore, calculation of drhogrid is commented out and not
      !     generalized to the spin polarized case.
      
      !        if(igrad) call dgemv('n',nint,((ng+1)*(ng+2))/2*(lcx+1),
      !    &           tt,rmtg,nint,rho,1,0.d0,drhogrd,1)
      
      
      
      !      the nonlinear core correction (nlcc)
      !      is added to the charge density
      !      prior to calling the xc drivers
      if(rcore>0d0)then
         do k=1,nint
            rhogrd(k,1) = rhogrd(k,1) + rhocore(k,1)
         end do
      end if
      
      if(rcore>0d0.and.nspol==2)then
         do k=1,nint
            rhogrd(k,2) = rhogrd(k,2) + rhocore(k,2)
         end do
      end if
      
      
      !     hutter
      !      call evxc(nint,rr,rhogrd,drhogrd,vxcgrd,excgrd)
      !     goedecker
      !     libxc wrapper
      !     call ggaenergy_15(nspol,nint,rw,rd,rhogrd,enexc,vxcgrd,excgrd)
      !     call ggaenergy_15(nspol,nint,rr,rw,rd,rhogrd,enexc,vxcgrd,excgrd)
      call drivexc(nspol,nint,rr,rw,rd,rhogrd,enexc,vxcgrd,excgrd)
      !        multiply with dr*r^2 to speed up calculation of matrix elements
      !        open(11,file='rhogrd')
      do k=1,nint
         vxcgrd(k,:)=vxcgrd(k,:)*rw(k)/fourpi
      enddo
      !        close(11)
      
      !      important: since the real space representation
      !      rhogrd will be used again only after recomputing
      !      it from the gaussian representation rho,
      !      there is no need to subtract the core charge here.
      
      !      exception: rhogrd will be passed to etot, 
      !      which needs the core charge to calculate 
      !      the exc energy term and then subtracts 
      !      rhocore from rhogrd to compute the other
      !      energy functionals (including vxc).
      
      !
      !     charge dependent part of hamiltonian
      !     hartree potential from valence charge distribution
      
      !     if nspol=2, the up and down potential differs due to vxc 
      
      !     do 4982,lp=0,lcx
      !     do 4982,jp=0,ng
      !     do 4982,ip=0,ng
      !     hhsc(i,j,l+1) =+ vh(ip,jp,lp+1,i,j,l+1)*rho(ip,jp,lp+1)
      
      call dgemv('t',(lcx+1)*((ng+1)*(ng+2))/2,  &
           (lmax+1)*((ng+1)*(ng+2))/2,1.d0,  &
           vh,(lcx+1)*((ng+1)*(ng+2))/2,rho(:,:,1),1,0.d0,hhsc,1)
      
      !     if nspol=2, add the hartree potential from the 2nd charge channel
      !     note: it seems easier to add the charges first and then do one dgemv.
      if(nspol==2)   &
           call dgemv('t',(lcx+1)*((ng+1)*(ng+2))/2,  &
           (lmax+1)*((ng+1)*(ng+2))/2,1.d0,  &
           vh,(lcx+1)*((ng+1)*(ng+2))/2,rho(:,:,2),1,1.d0,hhsc,1)
      !                                                  ^    ^ 
      
      !     potential from xc libraries 
      
      !     do 8049,k=1,nint
      !     8049 hhsc(i,j,l+1) =+ vxcgrd(k)*rmt(k,i,j,l+1)
      
      !     modification: if spin polarized, add this term to hhxc, not hhsc.
      !                   hxc is a spin polarized matrix only for that purpose.
      hhxc=0d0
      if(nspol==1)then
         call dgemv('t',nint,(lmax+1)*((ng+1)*(ng+2))/2,1.0d0,  &
              rmt,nint,vxcgrd(:,1),1,1.d0,hhsc,1)
      else
         call dgemv('t',nint,(lmax+1)*((ng+1)*(ng+2))/2,1.0d0,  &
              rmt,nint,vxcgrd(:,1),1,0.d0,hhxc(:,:,1),1)
         call dgemv('t',nint,(lmax+1)*((ng+1)*(ng+2))/2,1.0d0,  &
              rmt,nint,vxcgrd(:,2),1,1.d0,hhxc(:,:,2),1)
         
         !        spin polarized xc term end
      end if
      
      
      !     diagonalize
      do l=0,lmax
         do ispin=1,max(min(2*l+1,nspin), nspol)
            !     lapack
            ij=0
            do j=0,ng
               do i=j,ng
                  ij=ij+1
                  hht(i,j)=hh(i,j,l+1,ispin)+hhsc(ij,l+1)  &
                       +hhxc(ij,l+1,ispin) 
                  sst(i,j)=ss(i,j,l+1)
               enddo
            enddo
            !     ibm/dec
            call dsygv(1,'v','l',ng+1,hht,ng+1,sst,ng+1,  &
                 eval,evec,(ng+1)**2,info)
            !     the routine dsygv is also included in sub_lapack.f
            !     cray:
            !     call ssygv(1,'v','l',ng+1,hht,ng+1,sst,ng+1,
            !     1       eval,evec,(ng+1)**2,info)
            if (info.ne.0) write(6,*) 'lapack',info
            do iocc=0,noccmax-1
               do i=0,ng
                  evec(i,iocc)=hht(i,iocc)
               enddo
            enddo
            !     end lapack
            do iocc=1,noccmax
               !                 write(6,*)'debug: e(iocc,ispin,l,it)',
               !    :                       eval(iocc-1),iocc-1,ispin,l,it
               evsum=evsum+eval(iocc-1)
               aeval(iocc,l+1,ispin)=eval(iocc-1)
               do i=0,ng
                  psi(i,iocc,l+1,ispin)=evec(i,iocc-1)
               enddo
            enddo
            !     write(6,*) 'eval',l
            !     55         format(5(e14.7))
            !     write(6,55) eval
            !     write(6,*) 'evec',l
            !     do i=0,ng
            !     33            format(10(e9.2))
            !     write(6,33) (evec(i,iocc),iocc=0,noccmax-1)
            !     enddo
         enddo
      enddo
      tt=abs(evsum-evsumold)
      !        write(6,*)'debug: residue=',tt,it
      if (tt.lt.1.d-8) exit
   enddo
   if(tt>1d-8) write(6,*) 'warning: no sc convergence',tt
   
   !ccccccccccccccccccccccccccccccc
   !    end of the scf cycles     c
   !ccccccccccccccccccccccccccccccc
   
   !     write(*,*)'debug: ks eigenvalues',aeval
   itertot=itertot+it
   nscf = nscf +it
   call resid(nspol,  &
        noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,nsmx,  &
        aeval,res,  &
        hsep,  &
        ud,nint,ng,ngmx,psi,rho,pp1,pp2,pp3,  &
        potgrd,pexgrd,vxcgrd,rr,rw,ppr1,ppr2,ppr3,aux1,aux2,  &
        expxpr)
   !     etot evaluates ehartree using rhogrd,
   if (energ) call etot(verbose,nspol,  &
        noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,nsmx,  &
        aeval,  &
        rprb,zion,rloc,gpot,r_l,r_l2,hsep,  &
        xp,ud,nint,ng,ngmx,psi,rho,pp1,pp2,pp3,  &
        vxcgrd,excgrd,rhogrd,rhocore,occup,rr,rw,  &
        expxpr,etotal)
   !
   !     charge up to radius rcov or infinity
   !
   
   !     modification: can one integrate up to a different rcov for semicore states?
   !     problem: loop over l is implicit here, one can not tell easily where nl<occmax(l)
   
   if (lmax.gt.3) stop 'cannot calculate chrg'
   do l=0,lmax
      !        no spin down s orbitals in the r case: nspin=2, nspol=1, 2l+1=1
      do ispin=1,max(min(2*l+1,nspin), nspol)
         do iocc=1,noccmax
            chrg(iocc,l+1,ispin)=0.d0
            dhrg(iocc,l+1,ispin)=0.d0
            ehrg(iocc,l+1,ispin)=0.d0
         enddo
      enddo
   enddo
   do ispin=1,max(min(2*l+1,nspin),nspol)
      !     here, l=lmax+1, so do ispin=1,2 if lmax>0 and nspin=2
      do iocc=1,noccmax
         !        loop over all nl(l)
         loopij: do j=0,ng
            do i=0,ng
               d=xp(i)+xp(j)
               sd=sqrt(d)
               terf=derf(sd*rcov)
               texp=exp(-d*rcov**2)
               tt0=0.4431134627263791d0*terf/sd**3-0.5d0*rcov*texp/d
               tt1=0.6646701940895686d0*terf/sd**5 +  &  
                    (-0.75d0*rcov*texp - 0.5d0*d*rcov**3*texp)/d**2
               chrg(iocc,1,ispin)=chrg(iocc,1,ispin) +  &
                    psi(i,iocc,1,ispin)*psi(j,iocc,1,ispin)*tt0
               !     integrate up to rcov
               !               dhrg(iocc,1,ispin)=dhrg(iocc,1,ispin) +
               !     1              psi(i,iocc,1,ispin)*psi(j,iocc,1,ispin)*tt1
               !     integrate up to inf
               dhrg(iocc,1,ispin)=dhrg(iocc,1,ispin) +  &
                    psi(i,iocc,1,ispin)*psi(j,iocc,1,ispin)  &
                    *0.6646701940895686d0/sd**5
               ehrg(iocc,1,ispin)=ehrg(iocc,1,ispin) +  &
                    psi(i,iocc,1,ispin)*psi(j,iocc,1,ispin)  &
                    *1.66167548522392d0/sd**7
               if (lmax.eq.0) exit loopij
               
               tt2=1.661675485223921d0*terf/sd**7 +  &
                    (-1.875d0*rcov*texp-1.25d0*d*rcov**3*texp-  &
                    0.5d0*d**2*rcov**5*texp)/d**3
               chrg(iocc,2,ispin)=chrg(iocc,2,ispin) +  &
                    psi(i,iocc,2,ispin)*psi(j,iocc,2,ispin)*tt1
               !     integrate up to rcov
               !               dhrg(iocc,2,ispin)=dhrg(iocc,2,ispin) +
               !     1              psi(i,iocc,2,ispin)*psi(j,iocc,2,ispin)*tt2
               !     integrate up to inf
               dhrg(iocc,2,ispin)=dhrg(iocc,2,ispin) +  &
                    psi(i,iocc,2,ispin)*psi(j,iocc,2,ispin)  &
                    *1.661675485223921d0/sd**7
               ehrg(iocc,2,ispin)=ehrg(iocc,2,ispin) +  &
                    psi(i,iocc,2,ispin)*psi(j,iocc,2,ispin)  &
                    *5.815864198283725d0/sd**9
               if (lmax.eq.1) exit loopij
               
               tt3=5.815864198283725d0*terf/sd**9 +  &
                    (-6.5625d0*rcov*texp-4.375d0*d*rcov**3*texp-  &
                    1.75d0*d**2*rcov**5*texp -  &
                    0.5d0*d**3*rcov**7*texp)/d**4
               chrg(iocc,3,ispin)=chrg(iocc,3,ispin) +  &
                    psi(i,iocc,3,ispin)*psi(j,iocc,3,ispin)*tt2
               !     integrate up to rcov
               !               dhrg(iocc,3,ispin)=dhrg(iocc,3,ispin) +
               !     1              psi(i,iocc,3,ispin)*psi(j,iocc,3,ispin)*tt3
               !     integrate up to inf
               dhrg(iocc,3,ispin)=dhrg(iocc,3,ispin) +  &
                    psi(i,iocc,3,ispin)*psi(j,iocc,3,ispin)  &
                    *5.815864198283725d0/sd**9
               ehrg(iocc,3,ispin)=ehrg(iocc,3,ispin) +  &
                    psi(i,iocc,3,ispin)*psi(j,iocc,3,ispin)  &
                    *26.17138889227676d0/sd**11
               
               if (lmax.eq.2) exit loopij
               
               chrg(iocc,4,ispin)=chrg(iocc,4,ispin) +  &
                    psi(i,iocc,4,ispin)*psi(j,iocc,4,ispin)*tt3
               !     integrate up to rcov
               !                  tt4=26.17138889227676d0*terf/sd**11+(-29.53125d0*
               !     :                 rcov*texp-19.6875d0*d*rcov**3*texp-7.875d0*d**2
               !     :                 *rcov**5*texp-2.25d0*d**3*rcov**7*texp-  &
                    !     &                 0.5d0*d**4*rcov**9*texp)/d**5
               !               dhrg(iocc,4,ispin)=dhrg(iocc,4,ispin) +
               !     1              psi(i,iocc,4,ispin)*psi(j,iocc,4,ispin)*tt4
               !     integrate up to inf
               dhrg(iocc,4,ispin)=dhrg(iocc,4,ispin) +  &
                    psi(i,iocc,4,ispin)*psi(j,iocc,4,ispin)  &
                    *26.17138889227676d0/sd**11
               ehrg(iocc,4,ispin)=dhrg(iocc,4,ispin) +  &
                    psi(i,iocc,4,ispin)*psi(j,iocc,4,ispin)  &
                    *143.9426389075222d0/sd**13
               
            end do
         end do loopij
      end do
   enddo
   
   !
   !     value at origin
   !
   psir0=0.d0
   do i=0,ng
      psir0=psir0+psi(i,1,1,1)
   enddo
   psir0=psir0**2
   
   !     node locations of psi*r
   !     n-1 nodes allowed!
   !     search only for nodes if the corresponding weights are <> zero
   !     to avoid bumpy wavefunctions: no node of the first derivative of
   !     the pseudowavefunction*r  and only one node
   !     of the second derivative  up to the rmax of the lowest valence state
   !
   tol =1.0d-12
   ! initialize all elements of wfnode to zero, also unused ones. 
   ! this is only to test whether this helps to get rid of a bug
   wfnode=0d0
   do l=0,lmax
      !        no spin down s orbitals in the r case: nspin=2, nspol=1, 2l+1=1
      do ispin=1, min(min(2*l+1,nspin),nspol)
         do nocc=1,noccmax
            if ( (wght(nocc,l+1,ispin,6).ne.0.d0)  &
                 .or.  (wght(nocc,l+1,ispin,7).ne.0.d0)  &
                 .or.  (wght(nocc,l+1,ispin,8).ne.0.d0) ) then
               
               !       print*,'node search, l, nocc:',l,' ',nocc
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
               !     find outer max of psi, search from ~10 bohr down
               call detnp(nint,rr,10.0d0,kout)
               ttrmax=rr(kout)
               ra=ttrmax
               ttmax= dabs(wave2(ng,l,psi(0,nocc,l+1,ispin),  &
                    expxpr,ra,kout,nint))
               !      print*,'ttmax=',ttmax
               do k=kout,1, -1
                  ra= rr(k)
                  ttpsi= dabs(wave2(ng,l,psi(0,nocc,l+1,ispin),  &
                       expxpr,ra,k,nint))
                  if ( ttpsi .gt. ttmax  &
                       .and. ttpsi .gt. 1.0d-4 ) then
                     ttmax=ttpsi
                     ttrmax=ra
                  endif
                  if (ttpsi.lt.ttmax .and. ttpsi.gt.1.0d-4) exit
               enddo
               !     search up to 90% of rmax
               ttrmax=max(0.90d0*ttrmax,rr(1))
               call detnp(nint,rr,ttrmax,kout)
               ttrmax=rr(kout)
               !       print*,'search up to ',ttrmax,ttmax
               !     calc wavefunction and it's first two derivatives on the grid
               !
               do k=1,kout
                  call wave3(ng,l,xp,psi(0,nocc,l+1,ispin),  &
                       expxpr,rr(k),k,nint,y1(k),y2(k),y3(k))
               enddo
               
               do k = 2,kout
                  !     nodes of wavefunction
                  if (y1(k)*y1(k-1).lt.0.d0) then
                     nnode = nnode +1
                     x1=rr(k-1)
                     x2=rr(k)
                     rrnode = zbrent(wave,ng,ngmx,l,lmx,xp,psi,  &
                          nocc,noccmx,ispin,nsmx,  &
                          x1,x2,tol)
                     if (nnode .ge.nocc) then
                        rnode=rnode+rrnode
                        !                          print*,'found rnode at:',rrnode
                     endif
                     rlist(nnode)=rrnode
                  endif
                  !     nodes of first derivative
                  if (y2(k)*y2(k-1).lt.0.d0) then
                     ndnode = ndnode +1
                     x1=rr(k-1)
                     x2=rr(k)
                     rrnode = zbrent(dwave,ng,ngmx,l,lmx,xp,psi,  &
                          nocc,noccmx,ispin,nsmx,  &
                          x1,x2,tol)
                     if (ndnode .ge.nocc) then
                        dnode=dnode+rrnode
                        !                        print*,'found dnode at:',rrnode
                     endif
                     drlist(ndnode)=rrnode
                  endif
                  !     second derivative test:
                  if (y3(k)*y3(k-1).lt.0.d0) then
                     nddnode = nddnode + 1
                     x1=rr(k-1)
                     x2=rr(k)
                     rrnode = zbrent(ddwave,ng,ngmx,l,lmx,xp,psi,  &
                          nocc,noccmx,ispin,nsmx,  &
                          x1,x2,tol)
                     !     only add the lowest node! (this one shoud dissapear)
                     if (nddnode .ge. nocc +1 ) then
                        ddnode = ddnode + rrdnode
                        !                          print*,'found ddnode at:',rrnode
                     else
                        rrdnode=rrnode
                     endif
                     ddrlist(nddnode)=rrnode
                  endif
               enddo
               
               !     print*,'rnode,dnode,ddnode',rnode,dnode,ddnode,nnode
               
               !     new version: use integral of the relevant functions between the nodes
               !     not the node-locations!
               !     calc. necessary integrals:
               sum1=0.0d0
               sum2=0.0d0
               sum3=0.0d0
               !     rnodes:
               do i=nnode+1-nocc,1,-2
                  aa=wwav(ng,l,xp,psi(0,nocc,l+1,ispin),rlist(i))  &
                       -wwav(ng,l,xp,psi(0,nocc,l+1,ispin),rlist(i-1))
                  sum1 = sum1+aa
               enddo
               !     dnodes
               do i=ndnode+1-nocc,1,-2
                  aa=wave(ng,l,xp,psi(0,nocc,l+1,ispin),drlist(i))  &
                       -wave(ng,l,xp,psi(0,nocc,l+1,ispin),drlist(i-1))
                  sum2 = sum2+aa
               enddo
               !     ddnodes
               do i=nddnode+1-nocc,1,-2
                  aa=dwave(ng,l,xp,psi(0,nocc,l+1,ispin),ddrlist(i))  &
                       -dwave(ng,l,xp,psi(0,nocc,l+1,ispin),ddrlist(i-1))
                  !                    this test line is quite slow, for debuging purposes
                  sum3 = sum3+aa
               enddo
               !     old version for nodes as used in the paper:
               !                  wfnode(nocc,l+1,ispin,1)=rnode
               !                  wfnode(nocc,l+1,ispin,2)=dnode
               !                  wfnode(nocc,l+1,ispin,3)=ddnode
               !     new version, using the integrals of the function between the nodes
               wfnode(nocc,l+1,ispin,1)=sum1
               wfnode(nocc,l+1,ispin,2)=sum2
               wfnode(nocc,l+1,ispin,3)=sum3
               
               !                 some lines for bugfixing of wfnode 
               if(.not. sum1*sum1+sum2*sum2+sum3*sum3 >= 0d0)then
                  !                   let us use this condition as an isnan function
                  !                   that should work with all fortran compilers.
                  !                   indeed nan was observed sometimes for penalty terms from nodes
                  write(6,*)'ouch! some node integral is nan for'
                  write(6,*)'l=',l,' s=',ispin,' n-ncore(l)=',nocc
                  write(6,*)'(nan?)   node=',sum1,'   rnode=',rnode
                  write(6,*)'(nan?)  dnode=',sum2,'  drnode=',drnode
                  write(6,*)'(nan?) ddnode=',sum3,' ddrnode=',ddrnode
                  if(.not. sum1*sum1 >=0d0 )wfnode(nocc,l+1,ispin,1)=0d0
                  if(.not. sum2*sum2 >=0d0 )wfnode(nocc,l+1,ispin,2)=0d0
                  if(.not. sum3*sum3 >=0d0 )wfnode(nocc,l+1,ispin,3)=0d0
               end if
            endif
         enddo
      enddo
   enddo
   
   !     print*,'leave gatom'
   
end subroutine gatom

