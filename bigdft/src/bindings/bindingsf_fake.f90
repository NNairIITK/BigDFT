!> @file
!! Bindings for BigDFT
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


subroutine wf_emit_psi(c_obj, istep)
  implicit none
  integer(kind = 8), intent(in) :: c_obj
  integer, intent(in) :: istep
end subroutine wf_emit_psi
subroutine wf_emit_hpsi(c_obj, istep)
  implicit none
  integer(kind = 8), intent(in) :: c_obj
  integer, intent(in) :: istep
end subroutine wf_emit_hpsi
subroutine wf_emit_lzd(c_obj)
  implicit none
  integer(kind = 8), intent(in) :: c_obj
end subroutine wf_emit_lzd
subroutine wf_copy_from_fortran(c_obj, radii, crmult, frmult)
  use module_base
  implicit none
  integer(kind = 8), intent(in) :: c_obj
  real(gp), intent(in) :: crmult, frmult
  real(gp), intent(in) :: radii(:,:)
end subroutine wf_copy_from_fortran

subroutine localfields_emit_v_ext(c_obj)
  use module_types
  implicit none
  integer(kind = 8), intent(in) :: c_obj
END SUBROUTINE localfields_emit_v_ext
subroutine localfields_emit_rhov(c_obj, iter)
  use module_types
  implicit none
  integer(kind = 8), intent(in) :: c_obj
  integer, intent(in) :: iter
END SUBROUTINE localfields_emit_rhov

subroutine optloop_emit(optloop_c_obj, id, energs_c_obj)
  use module_types
  implicit none
  integer, intent(in) :: id
  integer(kind = 8), intent(in) :: optloop_c_obj, energs_c_obj
END SUBROUTINE optloop_emit

subroutine bigdft_signals_init(c_obj, type, domain, ln)
  implicit none
  integer(kind = 8), intent(out) :: c_obj
  integer, intent(in) :: type
  integer, intent(in) :: ln
  character(len = ln), intent(in) :: domain
END SUBROUTINE bigdft_signals_init
subroutine bigdft_signals_free(c_obj)
  implicit none
  integer(kind = 8), intent(in) :: c_obj
END SUBROUTINE bigdft_signals_free
subroutine bigdft_signals_start(c_obj, timeout)
  implicit none
  integer(kind = 8), intent(in) :: c_obj
  integer, intent(in) :: timeout
END SUBROUTINE bigdft_signals_start
subroutine bigdft_signals_stop(c_obj)
  implicit none
  integer(kind = 8), intent(in) :: c_obj
END SUBROUTINE bigdft_signals_stop
subroutine bigdft_signals_rm_denspot(c_obj)
  implicit none
  integer(kind = 8), intent(in) :: c_obj
END SUBROUTINE bigdft_signals_rm_denspot
subroutine bigdft_signals_rm_wf(c_obj)
  implicit none
  integer(kind = 8), intent(in) :: c_obj
END SUBROUTINE bigdft_signals_rm_wf
subroutine bigdft_signals_rm_energs(c_obj)
  implicit none
  integer(kind = 8), intent(in) :: c_obj
END SUBROUTINE bigdft_signals_rm_energs
subroutine bigdft_signals_rm_optloop(c_obj)
  implicit none
  integer(kind = 8), intent(in) :: c_obj
END SUBROUTINE bigdft_signals_rm_optloop
subroutine bigdft_signals_add_denspot(loop, c_obj)
  implicit none
  integer(kind = 8), intent(in) :: c_obj, loop
END SUBROUTINE bigdft_signals_add_denspot
subroutine bigdft_signals_add_wf(loop, c_obj)
  implicit none
  integer(kind = 8), intent(in) :: c_obj, loop
END SUBROUTINE bigdft_signals_add_wf
subroutine bigdft_signals_add_energs(loop, c_obj)
  implicit none
  integer(kind = 8), intent(in) :: c_obj, loop
END SUBROUTINE bigdft_signals_add_energs
subroutine bigdft_signals_add_optloop(loop, c_obj)
  implicit none
  integer(kind = 8), intent(in) :: c_obj, loop
END SUBROUTINE bigdft_signals_add_optloop

subroutine localfields_free_wrapper(c_obj)
  implicit none
  integer(kind = 8), intent(in) :: c_obj
END SUBROUTINE localfields_free_wrapper
subroutine wf_free_wrapper(c_obj)
  implicit none
  integer(kind = 8), intent(in) :: c_obj
END SUBROUTINE wf_free_wrapper
subroutine optloop_free_wrapper(c_obj)
  implicit none
  integer(kind = 8), intent(in) :: c_obj
END SUBROUTINE optloop_free_wrapper

subroutine wf_new_wrapper(c_obj, f_st)
  use module_types
  implicit none
  integer(kind = 8), intent(out) :: c_obj
  type(DFT_wavefunction), intent(in) :: f_st
END SUBROUTINE wf_new_wrapper
subroutine localfields_new_wrapper(c_obj, f_st)
  use module_types
  implicit none
  integer(kind = 8), intent(out) :: c_obj
  type(DFT_local_fields), intent(in) :: f_st
END SUBROUTINE localfields_new_wrapper
subroutine optloop_new_wrapper(c_obj, f_st)
  use module_types
  implicit none
  integer(kind = 8), intent(out) :: c_obj
  type(DFT_optimization_loop), intent(in) :: f_st
END SUBROUTINE optloop_new_wrapper
