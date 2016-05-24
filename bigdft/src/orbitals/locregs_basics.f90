!> @file
!! Basic operations relative to the localization regions
!! @author
!!    Copyright (C) 2007-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
subroutine psi_to_tpsi(hgrids,kptv,nspinor,lr,psi,w,hpsi,ekin,k_strten)
  use module_base
  use locregs, only: locreg_descriptors
  use locreg_operations, only: workarr_locham
  implicit none
  integer, intent(in) :: nspinor
  real(gp), dimension(3), intent(in) :: hgrids,kptv
  type(locreg_descriptors), intent(in) :: lr
  type(workarr_locham), intent(inout) :: w
  real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,nspinor), intent(in) :: psi
  real(gp), intent(out) :: ekin
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,nspinor), intent(inout) :: hpsi
  real(wp), dimension(6), optional :: k_strten
  !Local variables
  logical, parameter :: transpose=.false.
  logical :: usekpts
  integer :: idx,i,i_f,iseg_f,ipsif,isegf
  real(gp) :: ekino
  real(wp), dimension(0:3) :: scal
  real(gp), dimension(3) :: hgridh
  real(wp), dimension(6) :: kstrten,kstrteno


  !control whether the k points are to be used
  !real k-point different from Gamma still not implemented
  usekpts = nrm2(3,kptv(1),1) > 0.0_gp .or. nspinor == 2

  hgridh=.5_gp*hgrids

  do i=0,3
     scal(i)=1.0_wp
  enddo

  !starting point for the fine degrees, to avoid boundary problems
  i_f=min(1,lr%wfd%nvctr_f)
  iseg_f=min(1,lr%wfd%nseg_f)
  ipsif=lr%wfd%nvctr_c+i_f
  isegf=lr%wfd%nseg_c+iseg_f

  !call MPI_COMM_RANK(bigdft_mpi%mpi_comm,iproc,ierr)
  ekin=0.0_gp

  kstrten=0.0_wp
  select case(lr%geocode)
  case('F')

     !here kpoints cannot be used (for the moment, to be activated for the 
     !localisation region scheme
     if (usekpts) stop 'K points not allowed for Free BC locham'

     do idx=1,nspinor
        call uncompress_forstandard(lr%d%n1,lr%d%n2,lr%d%n3,&
             lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3,  & 
             lr%wfd%nseg_c,lr%wfd%nvctr_c,&
             lr%wfd%keygloc(1,1),lr%wfd%keyvloc(1),  & 
             lr%wfd%nseg_f,lr%wfd%nvctr_f,&
             lr%wfd%keygloc(1,isegf),lr%wfd%keyvloc(isegf),   &
             scal,psi(1,idx),psi(ipsif,idx),  &
             w%x_c(1,idx),w%x_f(1,idx),&
             w%x_f1(1,idx),w%x_f2(1,idx),w%x_f3(1,idx))

        call f_zero(w%nyc,w%y_c(1,idx))
        call f_zero(w%nyf,w%y_f(1,idx))

        call ConvolkineticT(lr%d%n1,lr%d%n2,lr%d%n3,&
             lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3,  &
             hgrids(1),hgrids(2),hgrids(3), &        !here the grid spacings are supposed to be equal. SM: not any more
             lr%bounds%kb%ibyz_c,lr%bounds%kb%ibxz_c,lr%bounds%kb%ibxy_c,&
             lr%bounds%kb%ibyz_f,lr%bounds%kb%ibxz_f,lr%bounds%kb%ibxy_f, &
             w%x_c(1,idx),w%x_f(1,idx),&
             w%y_c(1,idx),w%y_f(1,idx),ekino, &
             w%x_f1(1,idx),w%x_f2(1,idx),w%x_f3(1,idx),111)
        ekin=ekin+ekino

        !new compression routine in standard form
        call compress_and_accumulate_standard(lr%d,lr%wfd,&
             lr%wfd%keyvloc(1),lr%wfd%keyvloc(isegf),&
             lr%wfd%keygloc(1,1),lr%wfd%keygloc(1,isegf),&
             w%y_c(1,idx),w%y_f(1,idx),&
             hpsi(1,idx),hpsi(ipsif,idx))

     end do

  case('S')

     if (usekpts) then
        !first calculate the proper arrays then transpose them before passing to the
        !proper routine
        do idx=1,nspinor
           call uncompress_slab(lr%d%n1,lr%d%n2,lr%d%n3,&
                lr%wfd%nseg_c,lr%wfd%nvctr_c,&
                lr%wfd%keygloc(1,1),lr%wfd%keyvloc(1),   &
                lr%wfd%nseg_f,lr%wfd%nvctr_f,&
                lr%wfd%keygloc(1,isegf),lr%wfd%keyvloc(isegf),   &
                psi(1,idx),psi(ipsif,idx),w%x_c(1,idx),w%y_c(1,idx))
        end do

        !Transposition of the work arrays (use y_c as workspace)
        call transpose_for_kpoints(nspinor,2*lr%d%n1+2,2*lr%d%n2+31,2*lr%d%n3+2,&
             w%x_c,w%y_c,.true.)
        call f_zero(nspinor*w%nyc,w%y_c(1,1))

        ! compute the kinetic part and add  it to psi_out
        ! the kinetic energy is calculated at the same time
        ! do this thing for both components of the spinors
        do idx=1,nspinor,2
           call convolut_kinetic_slab_T_k(2*lr%d%n1+1,2*lr%d%n2+15,2*lr%d%n3+1,&
                hgridh,w%x_c(1,idx),w%y_c(1,idx),ekino,kptv(1),kptv(2),kptv(3))
           ekin=ekin+ekino        
        end do

        !re-Transposition of the work arrays (use x_c as workspace)
        call transpose_for_kpoints(nspinor,2*lr%d%n1+2,2*lr%d%n2+31,2*lr%d%n3+2,&
             w%y_c,w%x_c,.false.)

        do idx=1,nspinor
           !new compression routine in mixed form
           call analyse_slab_self(lr%d%n1,lr%d%n2,lr%d%n3,&
                w%y_c(1,idx),w%x_c(1,idx))
           call compress_and_accumulate_mixed(lr%d,lr%wfd,&
                lr%wfd%keyvloc(1),lr%wfd%keyvloc(isegf),&
                lr%wfd%keygloc(1,1),lr%wfd%keygloc(1,isegf),&
                w%x_c(1,idx),hpsi(1,idx),hpsi(ipsif,idx))

        end do

     else
        do idx=1,nspinor
           call uncompress_slab(lr%d%n1,lr%d%n2,lr%d%n3,&
                lr%wfd%nseg_c,lr%wfd%nvctr_c,&
                lr%wfd%keygloc(1,1),lr%wfd%keyvloc(1),   &
                lr%wfd%nseg_f,lr%wfd%nvctr_f,&
                lr%wfd%keygloc(1,isegf),lr%wfd%keyvloc(isegf),   &
                psi(1,idx),psi(ipsif,idx),w%x_c(1,idx),w%y_c(1,idx))

           call f_zero(w%nyc,w%y_c(1,idx))
           ! compute the kinetic part and add  it to psi_out
           ! the kinetic energy is calculated at the same time
           call convolut_kinetic_slab_T(2*lr%d%n1+1,2*lr%d%n2+15,2*lr%d%n3+1,&
                hgridh,w%x_c(1,idx),w%y_c(1,idx),ekino)
           ekin=ekin+ekino

           !new compression routine in mixed form
           call analyse_slab_self(lr%d%n1,lr%d%n2,lr%d%n3,&
                w%y_c(1,idx),w%x_c(1,idx))
           call compress_and_accumulate_mixed(lr%d,lr%wfd,&
                lr%wfd%keyvloc(1),lr%wfd%keyvloc(isegf),&
                lr%wfd%keygloc(1,1),lr%wfd%keygloc(1,isegf),&
                w%x_c(1,idx),hpsi(1,idx),hpsi(ipsif,idx))
        end do
     end if

  case('P')

     if (lr%hybrid_on) then

        !here kpoints cannot be used, such BC are used in general to mimic the Free BC
        if (usekpts) stop 'K points not allowed for hybrid BC locham'

        !here the grid spacing is not halved
        hgridh=hgrids
        do idx=1,nspinor
           call uncompress_per_f(lr%d%n1,lr%d%n2,lr%d%n3,&
                lr%wfd%nseg_c,lr%wfd%nvctr_c,&
                lr%wfd%keygloc(1,1),lr%wfd%keyvloc(1),   &
                lr%wfd%nseg_f,lr%wfd%nvctr_f,&
                lr%wfd%keygloc(1,isegf),lr%wfd%keyvloc(isegf),   &
                psi(1,idx),psi(ipsif,idx),w%x_c(1,idx),w%x_f(1,idx),&
                w%x_f1(1,idx),w%x_f2(1,idx),w%x_f3(1,idx),&
                lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3)

           call f_zero(w%nyc,w%y_c(1,idx))
           call f_zero(w%nyf,w%y_f(1,idx))

           call convolut_kinetic_hyb_T(lr%d%n1,lr%d%n2,lr%d%n3, &
                lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3,  &
                hgridh,w%x_c(1,idx),w%x_f(1,idx),w%y_c(1,idx),w%y_f(1,idx),kstrteno,&
                w%x_f1(1,idx),w%x_f2(1,idx),w%x_f3(1,idx),lr%bounds%kb%ibyz_f,&
                lr%bounds%kb%ibxz_f,lr%bounds%kb%ibxy_f)
           kstrten=kstrten+kstrteno
           !ekin=ekin+ekino

           call compress_and_accumulate_standard(lr%d,lr%wfd,&
                lr%wfd%keyvloc(1),lr%wfd%keyvloc(isegf),&
                lr%wfd%keygloc(1,1),lr%wfd%keygloc(1,isegf),&
                w%y_c(1,idx),w%y_f(1,idx),hpsi(1,idx),hpsi(ipsif,idx))

        end do
     else

        if (usekpts) then

           do idx=1,nspinor
              call uncompress_per(lr%d%n1,lr%d%n2,lr%d%n3,&
                   lr%wfd%nseg_c,lr%wfd%nvctr_c,&
                   lr%wfd%keygloc(1,1),lr%wfd%keyvloc(1),   &
                   lr%wfd%nseg_f,lr%wfd%nvctr_f,&
                   lr%wfd%keygloc(1,isegf),lr%wfd%keyvloc(isegf),   &
                   psi(1,idx),psi(ipsif,idx),w%x_c(1,idx),w%y_c(1,idx))
           end do

           if (transpose) then
              !Transposition of the work arrays (use psir as workspace)
              call transpose_for_kpoints(nspinor,2*lr%d%n1+2,2*lr%d%n2+2,2*lr%d%n3+2,&
                   w%x_c,w%y_c,.true.)

              call f_zero(w%y_c)
              ! compute the kinetic part and add  it to psi_out
              ! the kinetic energy is calculated at the same time
              do idx=1,nspinor,2
                 !print *,'AAA',2*lr%d%n1+1,2*lr%d%n2+1,2*lr%d%n3+1,hgridh

                 call convolut_kinetic_per_T_k(2*lr%d%n1+1,2*lr%d%n2+1,2*lr%d%n3+1,&
                      hgridh,w%x_c(1,idx),w%y_c(1,idx),kstrteno,kptv(1),kptv(2),kptv(3))
                 kstrten=kstrten+kstrteno
                 !ekin=ekin+ekino
              end do

              !Transposition of the work arrays (use psir as workspace)
              call transpose_for_kpoints(nspinor,2*lr%d%n1+2,2*lr%d%n2+2,2*lr%d%n3+2,&
                   w%y_c,w%x_c,.false.)

           else
              call f_zero(w%y_c)
              do idx=1,nspinor,2
                 call convolut_kinetic_per_T_k_notranspose(2*lr%d%n1+1,2*lr%d%n2+1,2*lr%d%n3+1,&
                      hgridh,w%x_c(1,idx),w%y_c(1,idx),kstrteno,kptv(1),kptv(2),kptv(3))
                 kstrten=kstrten+kstrteno
              end do
           end if

           do idx=1,nspinor

              call analyse_per_self(lr%d%n1,lr%d%n2,lr%d%n3,&
                   w%y_c(1,idx),w%x_c(1,idx))
              call compress_and_accumulate_mixed(lr%d,lr%wfd,&
                   lr%wfd%keyvloc(1),lr%wfd%keyvloc(isegf),&
                   lr%wfd%keygloc(1,1),lr%wfd%keygloc(1,isegf),&
                   w%x_c(1,idx),hpsi(1,idx),hpsi(ipsif,idx))

           end do
        else
           !first calculate the proper arrays then transpose them before passing to the
           !proper routine
           do idx=1,nspinor
              call uncompress_per(lr%d%n1,lr%d%n2,lr%d%n3,&
                   lr%wfd%nseg_c,lr%wfd%nvctr_c,&
                   lr%wfd%keygloc(1,1),lr%wfd%keyvloc(1),   &
                   lr%wfd%nseg_f,lr%wfd%nvctr_f,&
                   lr%wfd%keygloc(1,isegf),lr%wfd%keyvloc(isegf),   &
                   psi(1,idx),psi(ipsif,idx),w%x_c(1,idx),w%y_c(1,idx))

              call f_zero(w%nyc,w%y_c(1,idx))
              ! compute the kinetic part and add  it to psi_out
              ! the kinetic energy is calculated at the same time
              call convolut_kinetic_per_t(2*lr%d%n1+1,2*lr%d%n2+1,2*lr%d%n3+1,&
                   hgridh,w%x_c(1,idx),w%y_c(1,idx),kstrteno)
              kstrten=kstrten+kstrteno

              call analyse_per_self(lr%d%n1,lr%d%n2,lr%d%n3,&
                   w%y_c(1,idx),w%x_c(1,idx))
              call compress_and_accumulate_mixed(lr%d,lr%wfd,&
                   lr%wfd%keyvloc(1),lr%wfd%keyvloc(isegf),&
                   lr%wfd%keygloc(1,1),lr%wfd%keygloc(1,isegf),&
                   w%x_c(1,idx),hpsi(1,idx),hpsi(ipsif,idx))

           end do
        end if

     end if
     ekin=ekin+kstrten(1)+kstrten(2)+kstrten(3)
     if (present(k_strten)) k_strten=kstrten 

  end select

END SUBROUTINE psi_to_tpsi


!> In 3d,            
!! Applies the magic filter transposed, then analysis wavelet transformation.
!! The size of the data is forced to shrink
!! The input array y is not overwritten
subroutine comb_shrink_hyb(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w1,w2,y,xc,xf,sb)
  use module_defs, only: wp
  use locregs, only: shrink_bounds
  implicit none
  type(shrink_bounds),intent(in):: sb
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  real(wp), dimension(0:2*n1+1,0:2*n2+1,0:2*n3+1), intent(in) :: y
  real(wp), dimension(max(2*(2*n2+2)*(2*n3+2)*(nfu1-nfl1+1),&
       (2*n2+2)*(2*n3+2)*(n1+1))), intent(inout) :: w1
  real(wp), dimension(max(4*(2*n3+2)*(nfu1-nfl1+1)*(nfu2-nfl2+1),&
       (2*n3+2)*(n1+1)*(n2+1))), intent(inout) :: w2
  real(wp), dimension(0:n1,0:n2,0:n3), intent(inout) :: xc
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(inout) :: xf

  integer nt

  !perform the combined transform    

  call comb_shrink_hyb_c(n1,n2,n3,w1,w2,y,xc)

  ! I1,I2,I3 -> I2,I3,i1
  nt=(2*n2+2)*(2*n3+2)
  call comb_rot_shrink_hyb_1_ib(nt,n1,nfl1,nfu1,y,w1,sb%ibyyzz_f)

  ! I2,I3,i1 -> I3,i1,i2
  nt=(2*n3+2)*(nfu1-nfl1+1)
  call comb_rot_shrink_hyb_2_ib(nt,w1,w2,nfl2,nfu2,n2,sb%ibzzx_f)

  ! I3,i1,i2 -> i1,i2,i3
  nt=(nfu1-nfl1+1)*(nfu2-nfl2+1)
  call comb_rot_shrink_hyb_3_ib(nt,w2,xf,nfl3,nfu3,n3,sb%ibxy_ff)

END SUBROUTINE comb_shrink_hyb

subroutine comb_grow_all_hybrid(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nw1,nw2&
     ,w1,w2,xc,xf,y,gb)
  use module_defs, only: wp
  use locregs, only: grow_bounds
  implicit none
  type(grow_bounds),intent(in):: gb
  integer,intent(in)::n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nw1,nw2
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: xc
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in) :: xf
  real(wp), dimension(nw1), intent(inout) :: w1 !work
  real(wp), dimension(nw2), intent(inout) :: w2 ! work
  real(wp), dimension(0:2*n1+1,0:2*n2+1,0:2*n3+1), intent(out) :: y

  call comb_grow_c_simple(n1,n2,n3,w1,w2,xc,y)

  call comb_rot_grow_ib_1(n1      ,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,xf,w1,gb%ibyz_ff,gb%ibzxx_f)
  call comb_rot_grow_ib_2(n1,n2   ,          nfl2,nfu2,nfl3,nfu3,w1,w2,gb%ibzxx_f,gb%ibxxyy_f)
  call comb_rot_grow_ib_3(n1,n2,n3,                    nfl3,nfu3,w2,y,gb%ibxxyy_f)

END SUBROUTINE comb_grow_all_hybrid
