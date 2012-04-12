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

subroutine wf_new_from_fortran(c_obj, wf)
  use module_types
  implicit none
  type(DFT_wavefunction), intent(inout) :: wf
  double precision, intent(out) :: c_obj
END SUBROUTINE wf_new_from_fortran

subroutine bigdft_wf_free(c_obj)
  implicit none
  double precision, intent(in) :: c_obj
END SUBROUTINE bigdft_wf_free

subroutine bigdft_signals_start(c_obj)
  implicit none
  double precision, intent(in) :: c_obj
END SUBROUTINE bigdft_signals_start

subroutine bigdft_signals_stop(c_obj)
  implicit none
  double precision, intent(in) :: c_obj
END SUBROUTINE bigdft_signals_stop
