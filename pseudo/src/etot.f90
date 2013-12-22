      subroutine etot(verbose,nspol,  &
           noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,nsmx,  &
           aeval,  &
           rprb,zion,rloc,gpot,r_l,r_l2,hsep,  &
           xp,ud,nint,ng,ngmx,psi,rho,pp1,pp2,pp3,  &
           xcgrd,excgrd,rhogrd,rhocore,occup,rr,rw,  &
           expxpr,ekin_orb,dertwo,etotal)

      implicit real*8 (a-h,o-z)
      logical  energ, verbose
      dimension aeval(noccmx,lmx,nsmx),ekin_orb(noccmx,lmx,nsmx),  &
           gpot(4),r_l(lmx),r_l2(lmx),hsep(6,lpmx,nsmx),  &
           xp(0:ng),  ud(nint,((ng+1)*(ng+2))/2,lcx+1),   &
           psi(0:ngmx,noccmx,lmx,nsmx),  &
           rho(((ng+1)*(ng+2))/2,lmax+1,nspol),  &
           pp1(0:ng,lmax+1),pp2(0:ng,lmax+1),pp3(0:ng,lmax+1),  &
           xcgrd(nint,nspol), rhogrd(nint,nspol),rhocore(nint,nspol),  &
           occup(noccmx,lmx,nsmx),expxpr(0:ng,nint),  & 
           hhk(0:ng,0:ng) ! automatic array
      dimension scpr1(nsmx),scpr2(nsmx),scpr3(nsmx)

      dimension vexgrd(nint),excgrd(nint),vhgrd(nint)
      dimension rr(nint),rw(nint)
      real(kind=8) ehart
      external gamma


      fourpi = 16.d0*atan(1.d0)
      exc  = 0.0d0
      eext = 0.0d0
      !ekin_w = 0.0d0
      dertwo = 0.0d0
      vxc  = 0.0d0
      ehart= 0.0d0
      eigsum=0.0d0
      enl  = 0.0d0

!     VXC integral is with rho, EXC integral with rho+rhocore

!     calc Exc and Vxc first
!     rhogrd already contains the core charge at this point in case of a NLCC.
!     note: vxcgrd as returned by ggaenergy was already multifplied by rw/4pi in gatom
!     use the same loop to strip off the core charge 
!     and to add up the charges from  both channels.
!     We do not need polarization nor core charge for the other energy terms!
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
      call DGEMV('N',nint,((ng+1)*(ng+2))/2*(lcx+1),1.d0,ud,nint,  &
                   rho(1,1,1),1,0.d0,vhgrd,1)


      

!   calc eext and ehart
      do k=1,nint
         r=rr(k)
         vexgrd(k)=.5d0*(r/rprb**2)**2-zion*Derf(r/(sqrt(2.d0)*rloc))/r   &
              + exp(-.5d0*(r/rloc)**2)*  &
              ( gpot(1) + gpot(2)*(r/rloc)**2 + gpot(3)*(r/rloc)**4 +   &
              gpot(4)*(r/rloc)**6 )
!        write(88,*)rr(k),vhgrd(k),rhogrd(k,1)+rhogrd(k,2),rw(k)
         eext  = eext  +      vexgrd(k)*rhogrd(k,1)*rw(k)
         ehart = ehart + 0.5d0*vhgrd(k)*rhogrd(k,1)*rw(k)

      enddo

!     calc the kinetic energy term and the separable part of the pseudopotential

      ekin_orb=0.d0
      vxc=vxc*fourpi
      do ll=0,lmax
!     kinetic energy matrix
            gml=0.5d0*gamma(ll+0.5d0)
            do j=0,ng
               do i=j,ng
                  d=xp(i)+xp(j)
                  sxp=1.d0/d
                  const=gml*sqrt(sxp)**(2*ll+1)
                  hhij=.5d0*const*sxp**2*(3.d0*xp(i)*xp(j)+  &
                       ll*(6.d0*xp(i)*xp(j)-xp(i)**2-xp(j)**2) -  &
                       ll**2*(xp(i)-xp(j))**2  )+ .5d0*(ll*(ll+1))*const
                  hhk(i,j)=hhij
                  hhk(j,i)=hhij
               enddo
            enddo      

         do ispin=1,max(min(2*ll+1,nspin),nspol)
            if (ll.le.lpx-1) then
!       Note: Here sigma= r_l2(l) is used for the pp2 and pp3 projcetors
        rnrm1=1.d0/sqrt(.5d0*gamma(ll+1.5d0)*r_l(ll+1)**(2*ll+3))
        rnrm2=1.d0/sqrt(.5d0*gamma(ll+3.5d0)*r_l2(ll+1)**(2*ll+7))
        rnrm3=1.d0/sqrt(.5d0*gamma(ll+5.5d0)*r_l2(ll+1)**(2*ll+11))
            endif
            do iocc=1,noccmax
               zz = occup(iocc,ll+1,ispin)
               if (zz.gt.1.0d-8  ) then 
                  eigsum=eigsum + aeval(iocc,ll+1,ispin) *zz

!     kinetic energy
              tee=0.d0
              do i=0,ng
                te=0.d0
                do j=0,ng
                te=te+hhk(i,j)*psi(j,iocc,ll+1,ispin)
                enddo
                tee=tee+te*psi(i,iocc,ll+1,ispin)
              enddo
              ekin_orb(iocc,ll+1,ispin)=tee
        
!     separabel part
                  if (ll.le.lpx-1) then
                     scpr1(ispin)=DDOT(ng+1,psi(0,iocc,ll+1,ispin),  &
                          1,pp1(0,ll+1),1)
                     scpr2(ispin)=DDOT(ng+1,psi(0,iocc,ll+1,ispin),  &
                          1,pp2(0,ll+1),1)
                     scpr3(ispin)=DDOT(ng+1,psi(0,iocc,ll+1,ispin),  &
                          1,pp3(0,ll+1),1)
                  endif

                  psigrdold=wave2(ng,ll,psi(0,iocc,ll+1,ispin),expxpr,r,1,nint)
                  psigrd=psigrdold
                  rold=-rr(1)
                  r=0.d0
                  der2=0.d0
                  do k=1,nint
                     roldold=rold
                     rold=r
                     r=rr(k)
!     wavefunction on grid
                     psigrdoldold=psigrdold
                     psigrdold=psigrd
                     psigrd=wave2(ng,ll,psi(0,iocc,ll+1,ispin),expxpr,r,k,nint)
!     kinetic energy
                     d2=(psigrd-psigrdold)/(r-rold)
                     d2old=(psigrdold-psigrdoldold)/(rold-roldold)
                     if (k.ge.3) der2=der2+exp(-2.d0*rold)*.5d0*((d2-d2old)**2)/(r-roldold)
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
                  dertwo=dertwo+der2**2
                  if (verbose) write(6,'(a,i3,i3,i3,e10.3)') &
                    ' iocc,ll,ispin,der2 ',iocc,ll,ispin,der2
               endif
            enddo
         enddo
      enddo
      eh     =          ehart + vxc - exc - enl
      etotal = eigsum - ehart - vxc + exc
     
      if(verbose)then ! condition in previos versions was energ
         write(6,*) 
         write(6,*)' Pseudo atom energies'
         write(6,*)' ____________________'
         write(6,*) 
         write(6,'(a,e26.17)') ' Sec. Dervative Softness   =',dertwo
         write(6,'(a,f16.10)')' non local energy          =',enl
         write(6,'(a,f16.10)')' hartree energy with XC    =',eh
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
         write(6,'(a,e11.4)')' atomic+electronic charge ',zion-denfull
         write(6,*)
      end if
      !write(*,*)'From etot.f90 ekin=',ekin !Santanu
      return
      end
