  if (present(id)) then
     lgt=min(len(id),namelen)
     m%array_id(1:lgt)=id(1:lgt)
  end if
  if (present(routine_id)) then
     lgt=min(len(routine_id),namelen)
     m%routine_id(1:lgt)=routine_id(1:lgt)
  else
     m%routine_id=present_routine
  end if

  if(present(try)) m%try=try
