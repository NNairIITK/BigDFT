subroutine kswfn_free_scf_data(KSwfn, freePsit)
  use module_base
  use module_types
  use m_profiling
  implicit none
  type(DFT_wavefunction), intent(inout) :: KSwfn
  logical, intent(in) :: freePsit
  
  character(len = *), parameter :: subname = "kswfn_free_scf_data"
  integer :: i_all, i_stat

  ! Clean KSwfn parts only needed in the SCF loop.
  call deallocate_diis_objects(KSwfn%diis,subname)
  i_all=-product(shape(KSwfn%hpsi))*kind(KSwfn%hpsi)
  deallocate(KSwfn%hpsi,stat=i_stat)
  call memocc(i_stat,i_all,'hpsi',subname)
  if (freePsit) then
     i_all=-product(shape(KSwfn%psit))*kind(KSwfn%psit)
     deallocate(KSwfn%psit,stat=i_stat)
     call memocc(i_stat,i_all,'psit',subname)
  else
     nullify(KSwfn%psit)
  end if
end subroutine kswfn_free_scf_data
