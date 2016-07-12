  if(m%shared .eqv. .true.) then 
    p = smpi_shared_malloc(product(m%shape(1:(m%rank-1)))*(m%shape(m%rank)+ndebug)*f_sizeof(d),&
trim(m%array_id)//char(0), int8(m%rank))
    call c_f_pointer(p, array, m%shape(1:m%rank))
    ierror=0
  include 'allocate-inc.f90'
  return
  end if
