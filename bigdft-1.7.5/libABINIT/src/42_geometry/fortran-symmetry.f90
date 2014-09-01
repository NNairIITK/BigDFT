!try to interface with C
program symm
   use iso_c_binding
   implicit none
   integer :: nat =2
   integer, dimension(2) :: typat
   real(kind=8), dimension(3,2) :: xat
   type(c_ptr) :: cpointer
   typat(1) = 1
   typat(2) = 1
   xat(:,:) = 0.d0
   xat(1,2) = 1.d0
   xat(2,2) = 1.d0
   call find_symmetries(nat,typat,xat,cpointer)
end program symm
