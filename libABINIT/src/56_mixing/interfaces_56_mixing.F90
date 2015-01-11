!!****m* ABINIT/interfaces_56_mixing
!! NAME
!! interfaces_56_mixing
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/56_mixing
!!
!! COPYRIGHT
!! Copyright (C) 2010-2011 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! THIS FILE IS GENERATED AUTOMATICALLY BY abilint.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module interfaces_56_mixing

 implicit none

interface
 subroutine aprxdr(cplex,choice,dedv_mix,dedv_new,dedv_old,&  
  &  f_atm,f_fftgr,i_rhor2,i_vresid,moved_atm_inside,&  
  &  natom,nfft,nfftot,nspden,n_fftgr,rhor,ucvol,xred,fdot,user_data)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: choice
  integer,intent(in) :: cplex
  integer,intent(in) :: i_rhor2
  integer,intent(in) :: moved_atm_inside
  integer,intent(in) :: n_fftgr
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftot
  integer,intent(in) :: nspden
  real(dp),intent(out) :: dedv_mix
  real(dp),intent(out) :: dedv_new
  real(dp),intent(out) :: dedv_old
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: i_vresid(3)
  integer, intent(in) :: user_data(:)
  real(dp),intent(in) :: f_atm(3,natom,n_fftgr)
  real(dp),intent(in) :: f_fftgr(cplex*nfft,nspden,n_fftgr)
  real(dp),intent(in) :: rhor(cplex*nfft,nspden)
  real(dp),intent(in) :: xred(3,natom)

  interface
     function fdot(x,y,cplex,nfft,nspden,opt_denpot,user_data)
       integer, intent(in) :: cplex,nfft,nspden,opt_denpot
       double precision, intent(in) :: x(*), y(*)
       integer, intent(in) :: user_data(:)

       double precision :: fdot
     end function fdot
  end interface
 end subroutine aprxdr
end interface

interface
 subroutine dotprodm_v(cplex,cpldot,dot,index1,index2,mpi_comm,mpi_summarize,&  
  &  mult1,mult2,nfft,npot1,npot2,nspden,opt_storage,potarr1,potarr2)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: cpldot
  integer,intent(in) :: cplex
  integer,intent(in) :: index1
  integer,intent(in) :: index2
  integer,intent(in) :: mpi_comm
  integer,intent(in) :: mult1
  integer,intent(in) :: mult2
  integer,intent(in) :: nfft
  integer,intent(in) :: npot1
  integer,intent(in) :: npot2
  integer,intent(in) :: nspden
  integer,intent(in) :: opt_storage
  logical, intent(in) :: mpi_summarize
  real(dp),intent(out) :: dot(cpldot,mult1,mult2)
  real(dp),intent(in) :: potarr1(cplex*nfft,nspden,npot1)
  real(dp),intent(in) :: potarr2(cplex*nfft,nspden,npot2)
 end subroutine dotprodm_v
end interface

interface
 subroutine dotprodm_vn(cplex,cpldot,denarr,dot,id,ip,mpi_comm, mpi_summarize,multd,multp,&  
  &  nden,nfft,npot,nspden,potarr)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: cpldot
  integer,intent(in) :: cplex
  integer,intent(in) :: id
  integer,intent(in) :: ip
  integer,intent(in) :: mpi_comm
  integer,intent(in) :: multd
  integer,intent(in) :: multp
  integer,intent(in) :: nden
  integer,intent(in) :: nfft
  integer,intent(in) :: npot
  integer,intent(in) :: nspden
  logical, intent(in) :: mpi_summarize
  real(dp),intent(in) :: denarr(cplex*nfft,nspden,nden)
  real(dp),intent(out) :: dot(cpldot,multp,multd)
  real(dp),intent(in) :: potarr(cplex*nfft,nspden,npot)
 end subroutine dotprodm_vn
end interface

interface
 subroutine findminscf(choice,dedv_1,dedv_2,dedv_predict,&  
  &  d2edv2_1,d2edv2_2,d2edv2_predict,&  
  &  etotal_1,etotal_2,etotal_predict,&  
  &  lambda_1,lambda_2,lambda_predict,errid,errmess)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: choice
  integer,intent(out) :: errid
  real(dp),intent(out) :: d2edv2_1
  real(dp),intent(out) :: d2edv2_2
  real(dp),intent(out) :: d2edv2_predict
  real(dp),intent(in) :: dedv_1
  real(dp),intent(in) :: dedv_2
  real(dp),intent(out) :: dedv_predict
  character(len=500), intent(out) :: errmess
  real(dp),intent(in) :: etotal_1
  real(dp),intent(in) :: etotal_2
  real(dp),intent(out) :: etotal_predict
  real(dp),intent(in) :: lambda_1
  real(dp),intent(in) :: lambda_2
  real(dp),intent(out) :: lambda_predict
 end subroutine findminscf
end interface

interface
 subroutine scfcge(cplex,dbl_nnsclo,dtn_pc,etotal,f_atm,&  
  &  f_fftgr,initialized,iscf,isecur,istep,&  
  &  i_rhor,i_vresid,i_vrespc,moved_atm_inside,&  
  &  natom,nfft,nfftot,nspden,n_fftgr,n_index,opt_denpot,response,rhor,ucvol,vtrial,xred,&
  & fnrm,fdot,user_data,errid,errmess)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: cplex
  integer,intent(out) :: dbl_nnsclo
  integer,intent(out) :: errid
  integer,intent(in) :: initialized
  integer,intent(in) :: iscf
  integer,intent(in) :: isecur
  integer,intent(in) :: istep
  integer,intent(in) :: moved_atm_inside
  integer,intent(in) :: n_fftgr
  integer,intent(in) :: n_index
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nfftot
  integer,intent(in) :: nspden
  integer,intent(in) :: opt_denpot
  integer,intent(in) :: response
  integer, intent(in) :: user_data(:)
  character(len = 500), intent(out) :: errmess
  real(dp),intent(in) :: etotal
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: dtn_pc(3,natom)
  real(dp),intent(inout) :: f_atm(3,natom,n_fftgr)
  real(dp),intent(inout) :: f_fftgr(cplex*nfft,nspden,n_fftgr)
  integer,intent(inout) :: i_rhor(n_index)
  integer,intent(inout) :: i_vresid(n_index)
  integer,intent(inout) :: i_vrespc(n_index)
  real(dp),intent(in) :: rhor(cplex*nfft,nspden)
  real(dp),intent(inout) :: vtrial(cplex*nfft,nspden)
  real(dp),intent(inout) :: xred(3,natom)

  interface
     function fdot(x,y,cplex,nfft,nspden,opt_denpot,user_data)
       integer, intent(in) :: cplex,nfft,nspden,opt_denpot
       double precision, intent(in) :: x(*), y(*)
       integer, intent(in) :: user_data(:)

       double precision :: fdot
     end function fdot

     function fnrm(x,cplex,nfft,nspden,opt_denpot,user_data)
       integer, intent(in) :: cplex,nfft,nspden,opt_denpot
       double precision, intent(in) :: x(*)
       integer, intent(in) :: user_data(:)

       double precision :: fnrm
     end function fnrm
  end interface
 end subroutine scfcge
end interface

interface
 subroutine scfeig(istep,nfft,nspden,vrespc,vtrial,vtrial0,work,errid,errmess)
  use abi_defs_basis
  implicit none
  integer,intent(out) :: errid
  integer,intent(in) :: istep
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  character(len = 500), intent(out) :: errmess
  real(dp),intent(inout) :: vrespc(nfft,nspden)
  real(dp), intent(inout) :: vtrial(nfft,nspden)
  real(dp),intent(inout) :: vtrial0(nfft,nspden)
  real(dp),intent(inout) :: work(nfft,nspden,2)
 end subroutine scfeig
end interface

interface
 subroutine scfopt(cplex,f_fftgr,f_paw,iscf,istep,i_vrespc,i_vtrial,&  
  &  nfft,npawmix,nspden,n_fftgr,&  
  &  n_index,opt_denpot,pawoptmix,usepaw,vpaw,vresid,vtrial,fnrm,fdot,user_data,errid,errmess)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: cplex
  integer,intent(out) :: errid
  integer,intent(in) :: iscf
  integer,intent(in) :: istep
  integer,intent(in) :: n_fftgr
  integer,intent(in) :: n_index
  integer,intent(in) :: nfft
  integer,intent(in) :: npawmix
  integer,intent(in) :: nspden
  integer,intent(in) :: opt_denpot
  integer,intent(in) :: pawoptmix
  integer,intent(in) :: usepaw
  character(len = 500), intent(out) :: errmess
  real(dp),intent(out) :: vresid
  real(dp),intent(inout) :: f_fftgr(cplex*nfft,nspden,n_fftgr)
  real(dp),intent(inout) :: f_paw(npawmix,n_fftgr*usepaw)
  integer, intent(in) :: user_data(:)
  integer,intent(inout) :: i_vrespc(n_index)
  integer,intent(inout) :: i_vtrial(n_index)
  real(dp),intent(inout) :: vpaw(npawmix*usepaw)
  real(dp),intent(inout) :: vtrial(cplex*nfft,nspden)
  interface
     function fdot(x,y,cplex,nfft,nspden,opt_denpot,user_data)
       integer, intent(in) :: cplex,nfft,nspden,opt_denpot
       double precision, intent(in) :: x(*), y(*)
       integer, intent(in) :: user_data(:)

       double precision :: fdot
     end function fdot

     function fnrm(x,cplex,nfft,nspden,opt_denpot,user_data)
       integer, intent(in) :: cplex,nfft,nspden,opt_denpot
       double precision, intent(in) :: x(*)
       integer, intent(in) :: user_data(:)

       double precision :: fnrm
     end function fnrm
  end interface
 end subroutine scfopt
end interface

interface
 subroutine sqnormm_v(cplex,index,mpi_comm, mpi_summarize,mult,nfft,norm2,npot,nspden,opt_storage,potarr)
  use abi_defs_basis
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: index
  integer,intent(in) :: mpi_comm
  integer,intent(in) :: mult
  integer,intent(in) :: nfft
  integer,intent(in) :: npot
  integer,intent(in) :: nspden
  integer,intent(in) :: opt_storage
  logical, intent(in) :: mpi_summarize
  real(dp),intent(out) :: norm2(mult)
  real(dp),intent(in) :: potarr(cplex*nfft,nspden,npot)
 end subroutine sqnormm_v
end interface

end module interfaces_56_mixing
!!***
