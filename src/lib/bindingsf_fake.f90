subroutine denspot_emit_rhov(c_obj, istep)
  implicit none
  double precision, intent(in) :: c_obj
  integer, intent(in) :: istep
end subroutine denspot_emit_rhov

subroutine wf_emit_psi(c_obj, istep)
  implicit none
  double precision, intent(in) :: c_obj
  integer, intent(in) :: istep
end subroutine wf_emit_psi

subroutine energs_emit(c_obj, kind)
  use module_types
  implicit none
  double precision, intent(in) :: c_obj
  integer, intent(in) :: kind
END SUBROUTINE energs_emit
