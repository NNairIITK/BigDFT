
      real*8 function gamma(x)
c
c     calc. gamma(x), fuer halb/ganzzahlige positive x
c
      implicit real*8 (a-h,o-z)

c      print*,'entered gamma with x=',x
      if (x.le.0) then 
         write(6,*) 'stopped in gamma x<=0'
         stop
      endif
      if ( (x-int(x)).eq.0.d0) then
         n=int(x)
         if (n.eq.0) then 
            write(6,*) 'stopped in gamma x=0'
            stop
         else
            gam=1
            do i=1,n-1
               gam=gam*i
            enddo
         endif
      else
         xx=x-0.5d0
         if ( (xx-int(xx)).ne.0.d0) then
            write(6,*) 'stopped in gamma x<>n+1/2'
            write(6,*) 'x=',x
            stop
         endif
         n=int(xx)
         sqrtpi=sqrt(4.d0*atan(1.0d0))
         gam=sqrtpi
         do i=1,2*n-1,2
            gam=gam*i
         enddo
         if (n.gt.0) gam=gam/2**n
      endif
      gamma = gam
      return
      end
