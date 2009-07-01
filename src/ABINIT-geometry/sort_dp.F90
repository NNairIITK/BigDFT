 subroutine sort_dp(n,list,iperm,tol)

! Sort double precision array list(n) into ascending numerical order using Heapsort
! algorithm, while making corresponding rearrangement of the integer
! array iperm. Consider that two double precision numbers
! within tolerance tol are equal.

 implicit none
 integer n,iperm(n)
 double precision list(n)
 integer l,ir,iap,i,j
 double precision ap,tol

 if (n==1) then

! Accomodate case of array of length 1: already sorted!
  return

 else if (n<1) then

! Should not call with n<1
  write(06,1000) n
  1000  format(/,' sort_dp has been called with array length n=',i12,/, &
&  ' having a value less than 1.  This is not allowed.')
  stop

 else ! n>1

! Conduct the usual sort

  l=n/2+1
  ir=n

  do   ! Infinite do-loop

   if (l>1) then

    l=l-1
    ap=list(l)
    iap=iperm(l)

   else ! l<=1

    ap=list(ir)
    iap=iperm(ir)
    list(ir)=list(1)
    iperm(ir)=iperm(1)
    ir=ir-1

    if (ir==1) then
     list(1)=ap
     iperm(1)=iap
     exit   ! This is the end of this algorithm
    end if

   end if ! l>1

   i=l
   j=l+l

   do while (j<=ir) 
    if (j<ir) then
     if ( list(j)<list(j+1)-tol .or.  &
&        (list(j)<list(j+1)+tol.and.iperm(j)<iperm(j+1))) j=j+1
    endif
    if (ap<list(j)-tol .or. (ap<list(j)+tol.and.iap<iperm(j))) then
     list(i)=list(j)
     iperm(i)=iperm(j)
     i=j
     j=j+j
    else
     j=ir+1
    end if
   enddo

   list(i)=ap
   iperm(i)=iap

  enddo ! End infinite do-loop

 end if ! n>1

 end
