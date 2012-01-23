!!!  zero is already defined in atom.splines.f90
!!!  and atom.splines.f90 has been added to pseudo source
!!! because it is used by pawpatch.  AM 22/9/2011

c$$$        subroutine zero(n,x)
c$$$        implicit real*8 (a-h,o-z)
c$$$        dimension x(n)
c$$$        do 10,i=1,n-1,2
c$$$        x(i)=0.d0
c$$$        x(i+1)=0.d0
c$$$10      continue
c$$$        istart=i
c$$$        do 11,i=istart,n
c$$$        x(i)=0.d0
c$$$11      continue
c$$$        return
c$$$        end

