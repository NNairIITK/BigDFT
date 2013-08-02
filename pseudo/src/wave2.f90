        real*8 function wave2(ng,ll,psi,expxpr,r,k,nint)
        implicit real*8 (a-h,o-z)
        dimension psi(0:ng),expxpr(0:ng,nint)

        wave2=DDOT(ng+1,psi,1,expxpr(0,k),1)
!     cwh
        if (ll.gt.0) wave2=wave2*r**ll
        return
        end

