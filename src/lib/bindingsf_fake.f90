subroutine wf_emit_psi(c_obj, istep)
  implicit none
  double precision, intent(in) :: c_obj
  integer, intent(in) :: istep
end subroutine wf_emit_psi
subroutine wf_copy_from_fortran(c_obj, radii, crmult, frmult)
  use module_base
  implicit none
  double precision, intent(in) :: c_obj
  real(gp), intent(in) :: crmult, frmult
  real(gp), intent(in) :: radii(:,:)
end subroutine wf_copy_from_fortran

subroutine energs_emit(c_obj, kind)
  use module_types
  implicit none
  double precision, intent(in) :: c_obj
  integer, intent(in) :: kind
END SUBROUTINE energs_emit

subroutine localfields_emit_v_ext(c_obj)
  use module_types
  implicit none
  double precision, intent(in) :: c_obj
END SUBROUTINE localfields_emit_v_ext
subroutine localfields_emit_rhov(c_obj, iter)
  use module_types
  implicit none
  double precision, intent(in) :: c_obj
  integer, intent(in) :: iter
END SUBROUTINE localfields_emit_rhov

subroutine optloop_emit_iter(optloop, id, energs, iproc, nproc)
  use module_types
  implicit none
  type(DFT_optimization_loop), intent(inout) :: optloop
  type(energy_terms), intent(in) :: energs
  integer, intent(in) :: id, iproc, nproc
END SUBROUTINE optloop_emit_iter
subroutine optloop_emit_done(optloop, id, energs, iproc, nproc)
  use module_types
  implicit none
  type(DFT_optimization_loop), intent(inout) :: optloop
  type(energy_terms), intent(in) :: energs
  integer, intent(in) :: id, iproc, nproc
END SUBROUTINE optloop_emit_done

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

subroutine bigdft_signals_init(c_obj, type, domain, ln)
  implicit none
  double precision, intent(out) :: c_obj
  integer, intent(in) :: type
  integer, intent(in) :: ln
  character(len = ln), intent(in) :: domain
END SUBROUTINE bigdft_signals_init
subroutine bigdft_signals_free(c_obj)
  implicit none
  double precision, intent(in) :: c_obj
END SUBROUTINE bigdft_signals_free
subroutine bigdft_signals_start(c_obj, timeout)
  implicit none
  double precision, intent(in) :: c_obj
  integer, intent(in) :: timeout
END SUBROUTINE bigdft_signals_start
subroutine bigdft_signals_stop(c_obj)
  implicit none
  double precision, intent(in) :: c_obj
END SUBROUTINE bigdft_signals_stop
subroutine bigdft_signals_rm_denspot(c_obj)
  implicit none
  double precision, intent(in) :: c_obj
END SUBROUTINE bigdft_signals_rm_denspot
subroutine bigdft_signals_rm_wf(c_obj)
  implicit none
  double precision, intent(in) :: c_obj
END SUBROUTINE bigdft_signals_rm_wf
subroutine bigdft_signals_rm_energs(c_obj)
  implicit none
  double precision, intent(in) :: c_obj
END SUBROUTINE bigdft_signals_rm_energs
subroutine bigdft_signals_rm_optloop(c_obj)
  implicit none
  double precision, intent(in) :: c_obj
END SUBROUTINE bigdft_signals_rm_optloop
subroutine bigdft_signals_add_denspot(loop, c_obj)
  implicit none
  double precision, intent(in) :: c_obj, loop
END SUBROUTINE bigdft_signals_add_denspot
subroutine bigdft_signals_add_wf(loop, c_obj)
  implicit none
  double precision, intent(in) :: c_obj, loop
END SUBROUTINE bigdft_signals_add_wf
subroutine bigdft_signals_add_energs(loop, c_obj)
  implicit none
  double precision, intent(in) :: c_obj, loop
END SUBROUTINE bigdft_signals_add_energs
subroutine bigdft_signals_add_optloop(loop, c_obj)
  implicit none
  double precision, intent(in) :: c_obj, loop
END SUBROUTINE bigdft_signals_add_optloop

subroutine localfields_free_wrapper(c_obj)
  implicit none
  double precision, intent(in) :: c_obj
END SUBROUTINE localfields_free_wrapper
subroutine energs_free_wrapper(c_obj)
  implicit none
  double precision, intent(in) :: c_obj
END SUBROUTINE energs_free_wrapper
subroutine wf_free_wrapper(c_obj)
  implicit none
  double precision, intent(in) :: c_obj
END SUBROUTINE wf_free_wrapper
subroutine optloop_free_wrapper(c_obj)
  implicit none
  double precision, intent(in) :: c_obj
END SUBROUTINE optloop_free_wrapper

subroutine wf_new_wrapper(c_obj, f_st)
  use module_types
  implicit none
  double precision, intent(out) :: c_obj
  type(DFT_wavefunction), intent(in) :: f_st
END SUBROUTINE wf_new_wrapper
subroutine energs_new_wrapper(c_obj, f_st)
  use module_types
  implicit none
  double precision, intent(out) :: c_obj
  type(energy_terms), intent(in) :: f_st
END SUBROUTINE energs_new_wrapper
subroutine localfields_new_wrapper(c_obj, f_st)
  use module_types
  implicit none
  double precision, intent(out) :: c_obj
  type(DFT_local_fields), intent(in) :: f_st
END SUBROUTINE localfields_new_wrapper
subroutine optloop_new_wrapper(c_obj, f_st)
  use module_types
  implicit none
  double precision, intent(out) :: c_obj
  type(DFT_optimization_loop), intent(in) :: f_st
END SUBROUTINE optloop_new_wrapper
