        real*8 function xpown(r,n)
c     calc. r^n, n=int
c     this routine avoids poroblem with then CRAY
        implicit real*8 (a-h,o-z)
        if (n.eq.0) then
           xpown=1.d0
           return
        endif
        if (r.eq.0.d0) then
           xpown=0.d0
           return
        endif
        xpown=1.d0
        do i=1,abs(n)
           xpown=xpown*r
        enddo
        if (n.lt.0) xpown=1.d0/xpown
        return
        end



