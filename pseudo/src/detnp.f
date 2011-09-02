

        subroutine detnp(nn,r,rad0,npsp)
C       DETERMINE NPSP
        implicit real*8 (a-h,o-z)
        dimension r(nn)
        rmin=1.d10
        do 9237,i=1,nn
        if (abs(r(i)-rad0).lt.rmin) then
        rmin=abs(r(i)-rad0)
        npsp=i
        endif
9237    continue
        return
        end
