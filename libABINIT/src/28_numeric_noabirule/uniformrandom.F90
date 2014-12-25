 function uniformrandom(seed) 
! Returns a uniform random deviate between 0.0 and 1.0.
! Set seed to any value < 0 to initialize or reinitialize sequence.
! Parameters are chosen from integer overflow=2**23 (conservative).
! For some documentation, see Numerical Recipes, 1986, p196.

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_16_hideleave
!End of the abilint section

 implicit none
!Input/Output
 double precision :: uniformrandom
 integer :: seed
!Local variables
 integer, parameter :: im1=11979,ia1= 430,ic1=2531
 integer, parameter :: im2= 6655,ia2= 936,ic2=1399
 integer, parameter :: im3= 6075,ia3=1366,ic3=1283
 integer, save :: init=0
 integer, save :: ii1,ii2,ii3
 integer :: kk 
 double precision :: im1inv,im2inv
 double precision, save :: table(97)

 im1inv=1.0d0/im1 ; im2inv=1.0d0/im2

!Initialize on first call or when seed<0:
 if (seed<0.or.init==0) then
   seed=-abs(seed) 

!  First generator
   ii1=mod(ic1-seed,im1)
   ii1=mod(ia1*ii1+ic1,im1) 
!  Second generator
   ii2=mod(ii1,im2)
   ii1=mod(ia1*ii1+ic1,im1) 
!  Third generator
   ii3=mod(ii1,im3) 

!  Fill table
   do kk=1,97
     ii1=mod(ia1*ii1+ic1,im1)
     ii2=mod(ia2*ii2+ic2,im2) 
     table(kk)=(dble(ii1)+dble(ii2)*im2inv)*im1inv
   enddo

   init=1 ; seed=1
 end if 

!Third generator gives index
 ii3=mod(ia3*ii3+ic3,im3) 
 kk=1+(97*ii3)/im3
 if (kk<1.or.kk>97) then
   write(06, '(a,2i10,a)' ) ' trouble in uniformrandom; ii3,kk=',ii3,kk,' =>stop'
   call abi_leave_new('COLL')
 end if 
 uniformrandom=table(kk) 

!Replace old value, based on generators 1 and 2
 ii1=mod(ia1*ii1+ic1,im1)
 ii2=mod(ia2*ii2+ic2,im2)
 table(kk)=(dble(ii1)+dble(ii2)*im2inv)*im1inv

 return
 end
