        subroutine wave3(ng,ll,xp,psi,expxpr,r,k,nint,w,dw,ddw)
        implicit real*8 (a-h,o-z)
        dimension psi(0:ng),xp(0:ng),expxpr(0:ng,nint)

        w=0.d0
        dw=0.d0
        ddw=0.d0
        byr  = 1.d0/r
        byrr = byr/r
        rpll = xpown(r,ll)
        uu1  = -2.0d0*r*rpll
        uu2  = rpll*ll*byr
        uu3  = -2.0d0*rpll
        uu4  = +4.d0*rpll*r*r
        uu5  = -4.d0*rpll*ll
        uu6  = rpll*ll*(ll-1)*byrr
        do i=0,ng
           tt1=expxpr(i,k)
           tt2= psi(i)*tt1
           w=w + tt2
           tt3=tt2*xp(i)
           dw = dw +uu1*tt3 +uu2*tt2 
           tt4=tt3*xp(i)
           ddw=ddw +uu3*tt3+uu4*tt4+uu5*tt3+uu6*tt2
        enddo
!     cwh
        if (ll.gt.0) w=w*rpll
        return
        end



