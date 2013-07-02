subroutine radgrid(nrad,rr,rw,rd,a_grd,b_grd,rmax)
   implicit real*8 (a-h,o-z)
   dimension rr(nrad),rw(nrad),rd(nrad)
   ! generates logarithmic grid
   ! rr() radial grid: rr(i)=a_grd*b_grd**(i-1)-c_grd
   ! rw() weights for radial integration (dr/di)
   ! rd() di/dr
   !            ^ /4pir?
   
   fourpi=16.d0*atan(1.d0)
   do i=1,nrad
      rr(i)=a_grd*exp(b_grd*(i-1))
      rw(i)=b_grd*rr(i)
      rd(i)=1.d0/rw(i)
      rw(i)=rw(i)*fourpi*rr(i)**2
      if (rr(i).gt.rmax) exit
   end do
   if(rr(min(i,nrad))<rmax) then
      !write(6,*)'rmax too large, stopped in rradgrid'
      !stop
      write(6,*)'WARNING: The largest distance on the radial grid'
      write(6,*)'         may be too short. Consider to either'
      write(6,*)'         raise ng or to lower rprb.'
   end if
   nrad=i-1
   ! modify weights at en point for improved accuracy
   rw(1)=rw(1)*17.d0/48.d0
   rw(2)=rw(2)*59.d0/48.d0
   rw(3)=rw(3)*43.d0/48.d0
   rw(4)=rw(4)*49.d0/48.d0
   
   return
end subroutine radgrid
