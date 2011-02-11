        subroutine zero(n,x)
        implicit real*8 (a-h,o-z)
        dimension x(n)
        do 10,i=1,n-1,2
        x(i)=0.d0
        x(i+1)=0.d0
10      continue
        istart=i
        do 11,i=istart,n
        x(i)=0.d0
11      continue
        return
        end

