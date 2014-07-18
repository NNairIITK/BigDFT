module copy_utils
  use module_base
  implicit none

  private

  interface allocate_and_copy
    module procedure allocate_and_copy_i1, allocate_and_copy_i2, allocate_and_copy_d1, allocate_and_copy_d2
  end interface allocate_and_copy

  public :: allocate_and_copy


  contains

    subroutine allocate_and_copy_i1(array_in, array_out, id)
      implicit none
      ! Calling arguments
      integer,dimension(:),pointer,intent(in) :: array_in
      integer,dimension(:),pointer,intent(out) :: array_out
      character(len=*),intent(in) :: id
      ! Local variables
      integer :: iis1, iie1, ncount

      if (associated(array_out)) then
          call f_free_ptr(array_out)
      end if
      !array_out=f_malloc_ptr(src_ptr=array_in,id=id)
      if (associated(array_in)) then
          !array_out=f_malloc_ptr(lbounds=lbound(array_in),ubounds=ubound(array_in),id=id)
          !call f_memcpy(src=array_in,dest=array_out)
          iis1=lbound(array_in,1)
          iie1=ubound(array_in,1)
          ncount=iie1-iis1+1
          if (ncount>0) then
              array_out = f_malloc_ptr(iis1.to.iie1,id=id)
              call vcopy(ncount, array_in(iis1), 1, array_out(iis1), 1)
          else
              nullify(array_out)
          end if
      end if

    end subroutine allocate_and_copy_i1


    subroutine allocate_and_copy_i2(array_in, array_out, id)
      implicit none
      ! Calling arguments
      integer,dimension(:,:),pointer,intent(in) :: array_in
      integer,dimension(:,:),pointer,intent(out) :: array_out
      character(len=*),intent(in) :: id
      ! Local variables
      integer :: iis1, iie1, iis2, iie2, ncount

      if (associated(array_out)) then
          call f_free_ptr(array_out)
      end if
      if (associated(array_in)) then
          iis1=lbound(array_in,1)
          iie1=ubound(array_in,1)
          iis2=lbound(array_in,2)
          iie2=ubound(array_in,2)
          ncount=(iie1-iis1+1)*(iie2-iis2+1)
          if (ncount>0) then
              array_out = f_malloc_ptr((/iis1.to.iie1,iis2.to.iie2/),id=id)
              call vcopy(ncount, array_in(iis1,iis2), 1, array_out(iis1,iis2), 1)
          else
              nullify(array_out)
          end if
      end if

    end subroutine allocate_and_copy_i2


    subroutine allocate_and_copy_d1(array_in, array_out, id)
      implicit none
      ! Calling arguments
      double precision,dimension(:),pointer,intent(in) :: array_in
      double precision,dimension(:),pointer,intent(out) :: array_out
      character(len=*),intent(in) :: id
      ! Local variables
      integer :: iis1, iie1, ncount

      if (associated(array_out)) then
          call f_free_ptr(array_out)
      end if
      if (associated(array_in)) then
          iis1=lbound(array_in,1)
          iie1=ubound(array_in,1)
          ncount=iie1-iis1+1
          if (ncount>0) then
              array_out = f_malloc_ptr(iis1.to.iie1,id=id)
              call vcopy(ncount, array_in(iis1), 1, array_out(iis1), 1)
          else
              nullify(array_out)
          end if
      end if

    end subroutine allocate_and_copy_d1


    subroutine allocate_and_copy_d2(array_in, array_out, id)
      implicit none
      ! Calling arguments
      double precision,dimension(:,:),pointer,intent(in) :: array_in
      double precision,dimension(:,:),pointer,intent(out) :: array_out
      character(len=*),intent(in) :: id
      ! Local variables
      integer :: iis1, iie1, iis2, iie2, ncount

      if (associated(array_out)) then
          call f_free_ptr(array_out)
      end if
      if (associated(array_in)) then
          iis1=lbound(array_in,1)
          iie1=ubound(array_in,1)
          iis2=lbound(array_in,2)
          iie2=ubound(array_in,2)
          ncount=(iie1-iis1+1)*(iie2-iis2+1)
          if (ncount>0) then
              array_out = f_malloc_ptr((/iis1.to.iie1,iis2.to.iie2/),id=id)
              call vcopy(ncount, array_in(iis1,iis2), 1, array_out(iis1,iis2), 1)
          else
              nullify(array_out)
          end if
      end if

    end subroutine allocate_and_copy_d2

end module copy_utils
