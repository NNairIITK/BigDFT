       subroutine vionic(itype,iXC,ifcore,
     1 nrmax,nr,a,b,r,rab,rprb,lmax,
     2 nameat,norb,ncore,no,lo,so,zo,
     3 znuc,zsh,rsh,zel,zcore,cdd,cdu,cdc,
     4 viod,viou,vid,viu,vod,vou,
     5 etot,ev,ek,ep)
       implicit double precision(a-h,o-z)
c
c      vionic sets up the ionic potential
c      note that vio is the ionic potential times r
c
       dimension r(*),rab(*),
     2 no(norb),lo(norb),so(norb),zo(norb),
     3 cdd(*),cdu(*),cdc(*),
     4 viod(lmax,*),viou(lmax,*),vid(*),viu(*),vod(*),vou(*),
     5 etot(10),ev(norb),ek(norb),ep(norb)
       character*2 itype,nameat,icalc,cdtyp
       integer:: iXC
c
       dimension iray(6)
       character namef*6,iray*8,
     1 namet*2,icorrt*2,mcore*4,irel*3
c.....files
      common /files/iinput,iout,in290,in213,istore,iunit7,iunit8,istruc,
     +               ivnlkk,isumry,ikpts
c
c      2*znuc part (Rydberg units)
c
       ifcore = 0
       do 10 i=1,lmax
       do 12 j=1,nrmax
c  c.hartwig  add confining potential
          viod(i,j) = -2.0d0*(znuc -.5d0*(r(j)/rprb**2)**2*r(j))
          viou(i,j) = -2.0d0*(znuc -.5d0*(r(j)/rprb**2)**2*r(j))
c         viod(i,j) = -2.0*(       -.5d0*(r(j)/rprb**2)**2*r(j))
c         viou(i,j) = -2.0*(       -.5d0*(r(j)/rprb**2)**2*r(j))
c
c     c.hartwig  shift potential to avoid positive eigenvalues
c     and convergence problems
          vshift=-15.0d0*r(j)
          viod(i,j) = viod(i,j)+vshift
          viou(i,j) = viou(i,j)+vshift
 12    continue
 10   continue
c
c      add potential from shell charge
c
 105   if (zsh .eq. 0.D0) return
       do 110 i=1,lmax
       do 110 j=1,nr
       if (r(j) .ge. rsh) viod(i,j) = viod(i,j) - 2*zsh
       if (r(j) .ge. rsh) viou(i,j) = viou(i,j) - 2*zsh
       if (r(j) .lt. rsh) viod(i,j) = viod(i,j) - 2*zsh*r(i)/rsh
       if (r(j) .lt. rsh) viou(i,j) = viou(i,j) - 2*zsh*r(i)/rsh
 110   continue
       return
       end
