 subroutine abi_sort_int(n,list,iperm)

! Sort integer array list(n) into ascending numerical order using Heapsort
! algorithm, while making corresponding rearrangement of the integer
! array iperm. 

 implicit none
 integer n,l,ir,i,j,ip,ipp
 integer list(n),iperm(n)

 if (n==1) then

! Accomodate case of array of length 1: already sorted!
  return

 else if (n<1) then

! Should not call with n<1
  write(06,1000) n
  1000  format(/,' abi_sort_int has been called with array length n=',i12,/, &
&  ' having a value less than 1.  This is not allowed.')
  stop

 else ! n>1

! Conduct the usual sort

  l=n/2+1
  ir=n

  do   ! Infinite do-loop
 
   if (l>1) then

    l=l-1
    ip=list(l)
    ipp=iperm(l)

   else

    ip=list(ir)
    ipp=iperm(ir)
    list(ir)=list(1)
    iperm(ir)=iperm(1)
    ir=ir-1

    if (ir==1) then
     list(1)=ip
     iperm(1)=ipp
     exit   ! This is the end of this algorithm
    end if

   end if ! l>1

   i=l
   j=l+l

   do while (j<=ir)
    if (j<ir) then
     if (list(j).lt.list(j+1)) j=j+1
    end if
    if (ip.lt.list(j)) then
     list(i)=list(j)
     iperm(i)=iperm(j)
     i=j
     j=j+j
    else
     j=ir+1
    end if
   enddo

   list(i)=ip
   iperm(i)=ipp

  enddo ! End infinite do-loop

 end if ! n>1

 end
