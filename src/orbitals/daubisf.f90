!> @file
!!  Daubechies to Interpolation scaling functions routines
!! @author
!!    Copyright (C) 2010-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

!> Transform a daubechies function in compressed form to a function in real space via
!! the Magic Filter operation
!! do this for a single component (spinorial and/or complex)
subroutine daub_to_isf(lr,w,psi,psir)
  use module_base
  use locregs
  use locreg_operations
  implicit none
  type(locreg_descriptors), intent(in) :: lr
  type(workarr_sumrho), intent(inout) :: w
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f), intent(in) :: psi
  real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i), intent(out) :: psir
  !local variables
  integer :: i,i_f,iseg_f
  real(wp), dimension(0:3) :: scal

  call f_routine(id='daub_to_isf')

  i_f=min(lr%wfd%nvctr_f,1)
  iseg_f=min(lr%wfd%nseg_f,1)

  do i=0,3
     scal(i)=1.0_wp
  enddo

  select case(lr%geocode)
  case('F')
     call f_zero(psir)

     call uncompress_forstandard_short(lr%d%n1,lr%d%n2,lr%d%n3,&
          lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3,&
          lr%wfd%nseg_c,lr%wfd%nvctr_c,lr%wfd%keygloc(1,1),lr%wfd%keyvloc(1),  & 
          lr%wfd%nseg_f,lr%wfd%nvctr_f,&
          lr%wfd%keygloc(1,lr%wfd%nseg_c+iseg_f),lr%wfd%keyvloc(lr%wfd%nseg_c+iseg_f), &
          scal,psi(1),psi(lr%wfd%nvctr_c+i_f),&
          w%x_c,w%x_f)

     call comb_grow_all(lr%d%n1,lr%d%n2,lr%d%n3,&
          lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3,&
          w%w1,w%w2,w%x_c,w%x_f,  & 
          psir,lr%bounds%kb%ibyz_c,&
          lr%bounds%gb%ibzxx_c,lr%bounds%gb%ibxxyy_c,&
          lr%bounds%gb%ibyz_ff,lr%bounds%gb%ibzxx_f,lr%bounds%gb%ibxxyy_f)

  case('P')
     if (lr%hybrid_on) then
        ! hybrid case
        call uncompress_per_f_short(lr%d%n1,lr%d%n2,lr%d%n3,lr%wfd%nseg_c,&
             lr%wfd%nvctr_c,lr%wfd%keygloc(1,1),lr%wfd%keyvloc(1),& 
             lr%wfd%nseg_f,lr%wfd%nvctr_f,&
             lr%wfd%keygloc(1,lr%wfd%nseg_c+iseg_f),lr%wfd%keyvloc(lr%wfd%nseg_c+iseg_f), &
             psi(1),psi(lr%wfd%nvctr_c+i_f),&
             w%x_c,w%x_f,&
             lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3)

        call comb_grow_all_hybrid(lr%d%n1,lr%d%n2,lr%d%n3,&
             lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3,&
             w%nw1,w%nw2,w%w1,w%w2,w%x_c,w%x_f,psir,lr%bounds%gb)
     else
        call uncompress_per(lr%d%n1,lr%d%n2,lr%d%n3,lr%wfd%nseg_c,&
             lr%wfd%nvctr_c,lr%wfd%keygloc(1,1),lr%wfd%keyvloc(1),&
             lr%wfd%nseg_f,lr%wfd%nvctr_f,&
             lr%wfd%keygloc(1,lr%wfd%nseg_c+iseg_f),lr%wfd%keyvloc(lr%wfd%nseg_c+iseg_f),&
             psi,psi(lr%wfd%nvctr_c+i_f),&
             w%x_c,psir)

        call convolut_magic_n_per_self(2*lr%d%n1+1,2*lr%d%n2+1,2*lr%d%n3+1,&
             w%x_c,psir) 
     endif

  case('S')
     call uncompress_slab(lr%d%n1,lr%d%n2,lr%d%n3,lr%wfd%nseg_c,lr%wfd%nvctr_c,&
          lr%wfd%keygloc(1,1),lr%wfd%keyvloc(1),&
          lr%wfd%nseg_f,lr%wfd%nvctr_f,lr%wfd%keygloc(1,lr%wfd%nseg_c+iseg_f),&
          lr%wfd%keyvloc(lr%wfd%nseg_c+iseg_f), &
          psi(1),psi(lr%wfd%nvctr_c+i_f),w%x_c,&
          psir)

     call convolut_magic_n_slab_self(2*lr%d%n1+1,2*lr%d%n2+15,2*lr%d%n3+1,w%x_c,&
          psir) 

  end select

  call f_release_routine()


END SUBROUTINE daub_to_isf

!> Transforms a wavefunction written in real space basis into a 
!! wavefunction in Daubechies form
!! does the job for all supported BC
!! @warning: 
!!  - the psir is destroyed for some BCs (slab and periodic)
!!  - psi must already be initialized (to zero) before entering this routine
subroutine isf_to_daub(lr,w,psir,psi)
  use module_defs, only: wp
  use locregs
  use locreg_operations
  implicit none
  type(locreg_descriptors), intent(in) :: lr
  type(workarr_sumrho), intent(inout) :: w
  real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i), intent(inout) :: psir
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f), intent(inout) :: psi
  !local variables
  integer :: i,i_f,iseg_f,isegf,ipsif
  real(wp), dimension(0:3) :: scal

  do i=0,3
     scal(i)=1.0_wp
  enddo

  !starting point for the fine degrees, to avoid boundary problems
  i_f=min(1,lr%wfd%nvctr_f)
  iseg_f=min(1,lr%wfd%nseg_f)
  ipsif=lr%wfd%nvctr_c+i_f
  isegf=lr%wfd%nseg_c+iseg_f

  select case(lr%geocode)
  case('F')
     call comb_shrink(lr%d%n1,lr%d%n2,lr%d%n3,&
          lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3,&
          w%w1,w%w2,psir,&
          lr%bounds%kb%ibxy_c,lr%bounds%sb%ibzzx_c,lr%bounds%sb%ibyyzz_c,&
          lr%bounds%sb%ibxy_ff,lr%bounds%sb%ibzzx_f,lr%bounds%sb%ibyyzz_f,&
          w%x_c,w%x_f)

     call compress_and_accumulate_standard(lr%d,lr%wfd,&
          lr%wfd%keyvloc(1),lr%wfd%keyvloc(isegf),&
          lr%wfd%keygloc(1,1),lr%wfd%keygloc(1,isegf),&
          w%x_c,w%x_f,psi(1),psi(ipsif))

  case('S')

     call convolut_magic_t_slab_self(2*lr%d%n1+1,2*lr%d%n2+15,2*lr%d%n3+1,&
          psir(1),w%x_c(1))

     call analyse_slab_self(lr%d%n1,lr%d%n2,lr%d%n3,&
          w%x_c,psir(1))
     call compress_and_accumulate_mixed(lr%d,lr%wfd,&
          lr%wfd%keyvloc(1),lr%wfd%keyvloc(isegf),&
          lr%wfd%keygloc(1,1),lr%wfd%keygloc(1,isegf),&
          psir(1),psi(1),psi(ipsif))

  case('P')

     if (lr%hybrid_on) then

        call comb_shrink_hyb(lr%d%n1,lr%d%n2,lr%d%n3,&
             lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3,&
             w%w2,w%w1,psir(1),w%x_c(1),w%x_f(1),lr%bounds%sb)
        call compress_and_accumulate_standard(lr%d,lr%wfd,&
             lr%wfd%keyvloc(1),lr%wfd%keyvloc(isegf),&
             lr%wfd%keygloc(1,1),lr%wfd%keygloc(1,isegf),&
             w%x_c(1),w%x_f(1),psi(1),psi(ipsif))

     else

        call convolut_magic_t_per_self(2*lr%d%n1+1,2*lr%d%n2+1,2*lr%d%n3+1,&
             psir(1),w%x_c(1))

        call analyse_per_self(lr%d%n1,lr%d%n2,lr%d%n3,&
             w%x_c(1),psir(1))
        call compress_and_accumulate_mixed(lr%d,lr%wfd,&
             lr%wfd%keyvloc(1),lr%wfd%keyvloc(isegf),&
             lr%wfd%keygloc(1,1),lr%wfd%keygloc(1,isegf),&
             psir(1),psi(1),psi(ipsif))

     end if

  end select

END SUBROUTINE isf_to_daub

subroutine daub_to_isf_locham(nspinor,lr,w,psi,psir)
  use module_base, only: wp,f_zero
  use locregs
  use locreg_operations
  implicit none
  integer, intent(in) :: nspinor
  type(locreg_descriptors), intent(in) :: lr
  type(workarr_locham), intent(inout) :: w
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,nspinor), intent(in) :: psi
  real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,nspinor), intent(out) :: psir
  !local variables
  integer :: idx,i,i_f,iseg_f
  real(wp), dimension(0:3) :: scal


  do i=0,3
     scal(i)=1.0_wp
  enddo

  !starting point for the fine degrees, to avoid boundary problems
  i_f=min(1,lr%wfd%nvctr_f)
  iseg_f=min(1,lr%wfd%nseg_f)

  !call f_zero((2*n1+31)*(2*n2+31)*(2*n3+31)*nspinor,psir)
  !call MPI_COMM_RANK(bigdft_mpi%mpi_comm,iproc,ierr)
  select case(lr%geocode)
  case('F')
     call f_zero(psir)
     !call timing(iproc,'CrtDescriptors','ON') !temporary
     do idx=1,nspinor  
        call uncompress_forstandard(lr%d%n1,lr%d%n2,lr%d%n3,&
             lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3,  & 
             lr%wfd%nseg_c,lr%wfd%nvctr_c,&
             lr%wfd%keygloc(1,1),lr%wfd%keyvloc(1),  & 
             lr%wfd%nseg_f,lr%wfd%nvctr_f,&
             lr%wfd%keygloc(1,lr%wfd%nseg_c+iseg_f),lr%wfd%keyvloc(lr%wfd%nseg_c+iseg_f),   &
             scal,psi(1,idx),psi(lr%wfd%nvctr_c+i_f,idx),  &
             w%x_c(1,idx),w%x_f(1,idx),&
             w%x_f1(1,idx),w%x_f2(1,idx),w%x_f3(1,idx))
        !call timing(iproc,'CrtDescriptors','OF') !temporary
        !call timing(iproc,'CrtLocPot     ','ON') !temporary
        call comb_grow_all(lr%d%n1,lr%d%n2,lr%d%n3,&
             lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3,&
             w%w1,w%w2,w%x_c(1,idx),w%x_f(1,idx), & 
             psir(1,idx),lr%bounds%kb%ibyz_c,lr%bounds%gb%ibzxx_c,&
             lr%bounds%gb%ibxxyy_c,lr%bounds%gb%ibyz_ff,&
             lr%bounds%gb%ibzxx_f,lr%bounds%gb%ibxxyy_f)
        !call timing(iproc,'CrtLocPot     ','OF') !temporary
     end do

  case('S')

     do idx=1,nspinor
        call uncompress_slab(lr%d%n1,lr%d%n2,lr%d%n3,&
             lr%wfd%nseg_c,lr%wfd%nvctr_c,&
             lr%wfd%keygloc(1,1),lr%wfd%keyvloc(1),   &
             lr%wfd%nseg_f,lr%wfd%nvctr_f,&
             lr%wfd%keygloc(1,lr%wfd%nseg_c+iseg_f),lr%wfd%keyvloc(lr%wfd%nseg_c+iseg_f),   &
             psi(1,idx),psi(lr%wfd%nvctr_c+i_f,idx),w%x_c(1,idx),psir(1,idx))

        call convolut_magic_n_slab(2*lr%d%n1+1,2*lr%d%n2+15,2*lr%d%n3+1,w%x_c(1,idx),psir(1,idx),w%y_c(1,idx)) 

     end do

  case('P')

     if (lr%hybrid_on) then

        do idx=1,nspinor
           call uncompress_per_f(lr%d%n1,lr%d%n2,lr%d%n3,&
                lr%wfd%nseg_c,lr%wfd%nvctr_c,&
                lr%wfd%keygloc(1,1),lr%wfd%keyvloc(1),   &
                lr%wfd%nseg_f,lr%wfd%nvctr_f,&
                lr%wfd%keygloc(1,lr%wfd%nseg_c+iseg_f),lr%wfd%keyvloc(lr%wfd%nseg_c+iseg_f),   &
                psi(1,idx),psi(lr%wfd%nvctr_c+i_f,idx),w%x_c(1,idx),w%x_f(1,idx),&
                w%x_f1(1,idx),w%x_f2(1,idx),w%x_f3(1,idx),&
                lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3)

           call comb_grow_all_hybrid(lr%d%n1,lr%d%n2,lr%d%n3,&
                lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3,&
                w%nw1,w%nw2,&
                w%w1,w%w2,w%x_c(1,idx),w%x_f(1,idx),psir(1,idx),lr%bounds%gb)
        end do

     else

        do idx=1,nspinor
           call uncompress_per(lr%d%n1,lr%d%n2,lr%d%n3,&
                lr%wfd%nseg_c,lr%wfd%nvctr_c,&
                lr%wfd%keygloc(1,1),lr%wfd%keyvloc(1),   &
                lr%wfd%nseg_f,lr%wfd%nvctr_f,&
                lr%wfd%keygloc(1,lr%wfd%nseg_c+iseg_f),lr%wfd%keyvloc(lr%wfd%nseg_c+iseg_f),   &
                psi(1,idx),psi(lr%wfd%nvctr_c+i_f,idx),w%x_c(1,idx),psir(1,idx))

           call convolut_magic_n_per(2*lr%d%n1+1,2*lr%d%n2+1,2*lr%d%n3+1,w%x_c(1,idx),psir(1,idx),w%y_c(1,idx)) 
        end do

     end if
  end select

END SUBROUTINE daub_to_isf_locham


subroutine psig_to_psir_free(n1,n2,n3,work,psig_psir)
 implicit none
 integer, intent(in) :: n1,n2,n3 !< dimensions in the daubechies grid
 real(kind=8),dimension((2*n1+31)*(2*n2+31)*(2*n3+31)), intent(inout):: work !< enlarged buffer 
 real(kind=8),dimension((2*n1+31)*(2*n2+31)*(2*n3+31)), intent(inout):: psig_psir  !< final result, containing psig data but big enough to contain psir (used as work array)

 call synthese_free_self(n1,n2,n3,psig_psir,work)
 call convolut_magic_n_free_self(2*n1+15,2*n2+15,2*n3+15,work,psig_psir)

end subroutine psig_to_psir_free

