

      subroutine resid(nspol, &
           noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,nsmx, &
           aeval,res, &
           hsep, &
           ud,nint,ng,ngmx,psi,rho,pp1,pp2,pp3, &
           potgrd,pexgrd,xcgrd,rr,rw, &
           ppr1,ppr2,ppr3,aux1,aux2, &
           expxpr)

!          NEW
!          spin polarized treatment if nspol == 2

      implicit real*8 (a-h,o-z)
      dimension aeval(noccmx,lmx,nsmx),res(noccmx,lmx,nsmx), &
           hsep(6,lpmx,nsmx), &
           ud(nint,((ng+1)*(ng+2))/2,lcx+1), & 
           psi(0:ngmx,noccmx,lmx,nsmx), &
!          rho (ij,l+1, ispin)
           rho(((ng+1)*(ng+2))/2,lmax+1,nspol), &
           pp1(0:ng,lmax+1),pp2(0:ng,lmax+1),pp3(0:ng,lmax+1), &
           potgrd(nint),pexgrd(nint), &
!          Vxc (ij, ispin) ; xcgrid is actually a waste of memory
           xcgrd(nint,nspol),rr(nint),rw(nint), &
           ppr1(nint,lmax),ppr2(nint,lmax),ppr3(nint,lmax), &
           aux1(nint),aux2(nint,0:ng,lmax+1), &
           expxpr(0:ng,nint)

      dimension scpr1(nsmx),scpr2(nsmx),scpr3(nsmx)
      logical nothing
      external gamma

      nothing = .true.
      do ll=0,lmax
         do ispin=1,max(nspol,min(2*ll+1,nspin))
            do iocc=1,noccmax
               if (res(iocc,ll+1,ispin).ne.-1.0d0) nothing=.false.
            enddo
         enddo
      enddo
      if (nothing)  return
     
      fourpi = 16.d0*atan(1.d0)
!   external and exc potential on grid 
      do k=1,nint
         r=rr(k)
!        some older comment line?
!        potgrd(k)=pexgrd(k) + fourpi*xcgrd(k)/rw(k)

!        nonpolarized version was
!        potgrd(k)=pexgrd(k) + aux1(k)*xcgrd(k)

         potgrd(k)=pexgrd(k) 
!        the missing XC term is added later to the temp variable tt
         
      enddo

!
! add hartree potential

!        potgrd(k)=+ud(k,i,j,l+1)*rho(i,j,l+1)
      call DGEMV('N',nint,((ng+1)*(ng+2))/2*(lcx+1),1.d0,ud,nint, &
                   rho(1,1,1),1,1.d0,potgrd,1)
!     polarized: both channels of rho are needed for the Hartree term
      if(nspol.eq.2) &
      call DGEMV('N',nint,((ng+1)*(ng+2))/2*(lcx+1),1.d0,ud,nint, &
                   rho(1,1,2),1,1.d0,potgrd,1)
      do ll=0,lmax
!        if nspol=2, s orbitals have two spin states, too
         do ispin=1,max(min(2*ll+1,nspin),nspol)
            do iocc=1,noccmax
               if (res(iocc,ll+1,ispin).ne.-1.0d0) then
!     separabel part
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
!     wavefunction on grid
                     psigrd= wave2(ng,ll,psi(0,iocc,ll+1,ispin), &
                          expxpr,r,k,nint)
!     kinetic energy
                     rkin=0.d0
                     do i=0,ng
                        rkin=rkin+psi(i,iocc,ll+1,ispin) &
                             *aux2(k,i,ll+1)
                     enddo
                     rkin=rkin*r**ll
!     separabel part
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
! resdidue
!                    the term from VXC can be spin polarized,
!                    therefore  the index min(nspol, ispin)
                     tt=aux1(k)*xcgrd(k,min(nspol,ispin))
!                    add external and Ha potential, KS ev
                     tt=tt+potgrd(k)-aeval(iocc,ll+1,ispin)
!                    times psi, plus kinetic term and separable part
                     tt=psigrd*tt+rkin+sep
!                    write(22,*)'DEBUG resid:iocc,ll,ispin,tt,rkin,sep', &
!                                            iocc,ll,ispin, &
!                           (tt-rkin-sep)*rw(k) ,rkin*rw(k),sep*rw(k)
                     res(iocc,ll+1,ispin)=res(iocc,ll+1,ispin) &
                          +tt*tt*rw(k)/fourpi
                  enddo
               else
                  res(iocc,ll+1,ispin)=0.0d0
               endif
            enddo
         enddo
      enddo
      return
      end
