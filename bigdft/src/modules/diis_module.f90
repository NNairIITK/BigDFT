!> @file
!! DIIS and Steepest descent routines
!! @author
!!    Copyright (C) 2007-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Module handling DIIS (Direct Inversion in the Iterative Subspace) and SD (Steepest Descent) optimization
module diis_sd_optimization
  use module_base

  !> Control objects of DIIS procedure
  type, public :: DIIS_ctrl
    logical :: switchSD     !< Switch to Steepest Descent if .true.
     integer :: idiistol    !< Number of iterations when the energy is increasing
     real(gp) :: energy_min !< Minimal energy during the iterated process
     real(gp) :: energy_old !< Previous value already fulfilled
     real(gp) :: energy     !< Current value of energy
     real(gp) :: alpha_max  !< Maximal value of alpha (for step size with SD)
  end type DIIS_ctrl

  !> Contains the arguments needed for the diis procedure
  type, public :: DIIS_obj
     type(DIIS_ctrl) :: ctrl
     integer :: mids !< Size of the current DIIS history (or matrix) <= idsx
     integer :: ids  !< Iteration number
     integer :: idsx !< History of the diis (also if idiistol > idsx switch to SD)
     real(tp), dimension(:), pointer :: psidst   !< History of the given vectors (psi)
     real(tp), dimension(:), pointer :: hpsidst  !< History of the corresponding hpsi
     real(tp), dimension(:,:,:), pointer :: ads  !< DIIS matrix
     real(gp) :: alpha_coeff !< Mixing coefficient
  end type DIIS_obj

  private

  public :: diis_set,diis_free,diis_update_psi,diis_update_errors,diis_step

contains

!!$  function DIIS_ctrl_init() result(ctrl)
!!$    
!!$  end function DIIS_ctrl_init


  !> Allocate diis objects
  subroutine DIIS_set(idsx,alphaSD,ndim_psi,ngrpp,diis) !n(m)
    use module_base
    use module_types
    implicit none
    integer, intent(in) :: idsx,ndim_psi,ngrpp
    real(gp), intent(in) :: alphaSD
    type(DIIS_obj), intent(inout) :: diis

    diis%ids=0
    diis%mids=1
    !local variable for the diis history
    diis%idsx=idsx
    diis%psidst =f_malloc_ptr(ndim_psi*idsx,id='psidst')
    diis%hpsidst=f_malloc_ptr(ndim_psi*idsx,id='hpsidst')
    diis%ads    =f_malloc0_ptr((/idsx+1,idsx+1,ngrpp/),id='ads')

    !initialize scalar variables
    !diis initialisation variables
    diis%alpha_coeff=alphaSD
!!$    diis%alpha_max=alphadiis
!!$    diis%energy=1.d10
!!$    !minimum value of the energy during the minimisation procedure
!!$    diis%energy_min=1.d10
!!$    !previous value already fulfilled
!!$    diis%energy_old=diis%energy
!!$    !logical control variable for switch DIIS-SD
!!$    diis%switchSD=.false.
  END SUBROUTINE DIIS_set


  !> De-Allocate diis objects
  subroutine DIIS_free(diis)
    use module_base
    use module_types
    implicit none
    type(DIIS_obj), intent(inout) :: diis

    call f_free_ptr(diis%psidst)
    call f_free_ptr(diis%hpsidst)
    call f_free_ptr(diis%ads)

  END SUBROUTINE DIIS_free


  !> Fill the DIIS matrices with the error of the previous step
  subroutine DIIS_update_errors(ngrp,isgrp,ngrpp,ncomp_grp,ndim_psi,psi,hpsi,diis)
    use yaml_strings, only: yaml_toa
    implicit none 
    integer, intent(in) :: ngrp,isgrp,ngrpp
    integer, intent(in) :: ndim_psi !< should be greater or equal to sum(ncomp_grp(isgrp+1:isgrp+ngrpp)
    integer, dimension(ngrp), intent(in) :: ncomp_grp !< number of components per group
    real(wp), dimension(ndim_psi), intent(in) :: psi,hpsi
    type(DIIS_obj), intent(inout) :: diis
    !local variables
    integer :: ispsi,ispsidst,ncomp,igrpp,igrp

    if (f_err_raise(sum(ncomp_grp(isgrp+1:isgrp+ngrpp)) > ndim_psi,&
         'Size inconsistency in DIIS_update_errors, '//&
         trim(yaml_toa(sum(ncomp_grp(isgrp+1:isgrp+ngrpp))))//' > '//&
         trim(yaml_toa(ndim_psi)))) return

    ispsi=1
    ispsidst=1
    do igrpp=1,ngrpp
       igrp=isgrp+igrpp
       ncomp=ncomp_grp(igrp)
       if (ncomp == 0) cycle
       !here we can choose to store the DIIS arrays with single precision
       !psidst=psit
       call vcopy(ncomp,psi(ispsi),1,diis%psidst(ispsidst+ncomp*(diis%mids-1)),1)
       !hpsidst=hpsi
       call vcopy(ncomp,hpsi(ispsi),1,diis%hpsidst(ispsidst+ncomp*(diis%mids-1)),1)
       ispsi=ispsi+ncomp
       ispsidst=ispsidst+ncomp*diis%idsx
    end do

  end subroutine DIIS_update_errors


  !> Fill the psi array with the DIIS combination of previous errors
  subroutine DIIS_update_psi(ngrp,isgrp,ngrpp,ncomp_grp,ndim_psi,psi,diis)
    use yaml_strings, only: yaml_toa
    implicit none 
    integer, intent(in) :: ngrp,isgrp,ngrpp
    integer, intent(in) :: ndim_psi !< should be greater or equal to sum(ncomp_grp(isgrp+1:isgrp+ngrpp)
    integer, dimension(ngrp), intent(in) :: ncomp_grp !< number of components per group
    real(wp), dimension(ndim_psi), intent(inout) :: psi
    type(DIIS_obj), intent(inout) :: diis
    !local variables
    integer :: ispsi,ispsidst,ncomp,igrpp,igrp

    !update the psit array with the difference stored in the psidst work array
    ispsi=1
    ispsidst=1
    do igrpp=1,ngrpp
       igrp=isgrp+igrpp
       ncomp=ncomp_grp(igrp)
       if (ncomp == 0) cycle
       call axpy(ncomp,1.0_dp,&
            diis%psidst(ispsidst+(mod(diis%ids,diis%idsx))*ncomp),1,&
            psi(ispsi),1)
       ispsi=ispsi+ncomp
       ispsidst=ispsidst+ncomp*diis%idsx
    end do

  end subroutine DIIS_update_psi


  !> Calculates the DIIS extrapolated solution psit in the ids-th DIIS step 
  !! using the previous iteration points psidst and the associated error 
  !! vectors (preconditioned gradients) hpsidst
  subroutine diis_step(iproc,nproc,ngrp,isgrp,ngrpp,igrpproc,ncomp_grp,diis)
    use module_types
    use yaml_strings, only: yaml_toa
    implicit none
    ! Arguments
    integer, intent(in) :: nproc,iproc,ngrp,isgrp,ngrpp
    integer, dimension(ngrp), intent(in) :: igrpproc !<array which associate each group to only one iproc for broadcasting
    integer, dimension(ngrp), intent(in) :: ncomp_grp !< number of components per group
    type(DIIS_obj), intent(inout) :: diis
    ! Local variables
    character(len=*), parameter :: subname='diisstp'
    integer :: i,j,ist,jst,mi,info,jj,mj
    integer :: ispsi,ispsidst,ncomp,iacc_add,igrpp,igrp
    real(tp) :: psicoeff
    integer, dimension(:), allocatable :: ipiv
    real(tp), dimension(:,:), allocatable :: adsw
    real(tp), dimension(:,:), allocatable :: rds

    call f_routine(id=subname)

    ipiv=f_malloc(diis%idsx+1,id='ipiv')
    rds=f_malloc0((/diis%idsx+1,ngrp/),id='rds')
    adsw=f_malloc0((/diis%idsx+1,diis%idsx+1/),id='adsw')

    ispsidst=1
    do igrpp=1,ngrpp
       igrp=isgrp+igrpp
       ncomp=ncomp_grp(igrp)
       if (ncomp == 0) cycle
       ! set up DIIS matrix (upper triangle)
       if (diis%ids > diis%idsx) then
          ! shift left up matrix
          do i=1,diis%idsx-1
             do j=1,i
                diis%ads(j,i,igrpp)=diis%ads(j+1,i+1,igrpp)
             end do
          end do
       end if

       ! calculate new line, use rds as work array for summation
       ist=max(1,diis%ids-diis%idsx+1)
       do i=ist,diis%ids
          mi=mod(i-1,diis%idsx)+1
          !useful in care of more than one group
          rds(i-ist+1,igrp)=0.0_tp
          rds(i-ist+1,igrp)=rds(i-ist+1,igrp)+dot(ncomp,&
               diis%hpsidst(ispsidst+(diis%mids-1)*ncomp),1,&
               diis%hpsidst(ispsidst+(mi-1)*ncomp),1)
       end do
       ispsidst=ispsidst+ncomp*diis%idsx
    end do
    if (nproc > 1) then
       call mpiallred(rds,MPI_SUM,comm=bigdft_mpi%mpi_comm)
       if (f_err_raise(f_err_check(err_name='ERR_MPI_WRAPPERS'),&
            'Error in allreduce operation, '//subname,BIGDFT_MPI_ERROR)) then
          call free_and_exit()
          return
       end if
    endif

    ispsi=1
    ispsidst=1
    do igrpp=1,ngrpp
       igrp=isgrp+igrpp
       ncomp=ncomp_grp(igrp)
       if (ncomp == 0) cycle
          !update the matrix of the DIIS errors
          do i=1,min(diis%ids,diis%idsx)
             diis%ads(i,min(diis%idsx,diis%ids),igrpp)=rds(i,igrp)
          end do

          ! copy to work array, right hand side, boundary elements
          do j=1,min(diis%idsx,diis%ids)
             adsw(j,min(diis%idsx,diis%ids)+1)=0.0_tp
             adsw(j,min(diis%idsx,diis%ids)+1)=1.0_tp
             rds(j,igrp)=0.0_tp
             do i=j,min(diis%idsx,diis%ids)
                adsw(j,i)=diis%ads(j,i,igrpp)
             end do
          end do
          adsw(min(diis%idsx,diis%ids)+1,min(diis%idsx,diis%ids)+1)=0.0_tp
          rds(min(diis%idsx,diis%ids)+1,igrp)=1.0_tp

          !make the matrix symmetric (hermitian) to use DGESV (ZGESV) (no work array, more stable)
          do j=1,min(diis%idsx,diis%ids)+1
             do i=1,min(diis%idsx,diis%ids)+1
                adsw(i,j)=adsw(j,i)
             end do
          end do
          !if(iproc==0)  write(6,*) 'DIIS matrix'
          !do i=1,min(diis%idsx,ids)+1
          !  if(iproc==0)  write(6,'(i3,12(1x,e9.2))') iproc,(ads(i,j,2),j=1,min(diis%idsx,ids)+1),rds(i)
          !enddo
          if (diis%ids > 1) then
             ! solve linear system, supposing it is general. More stable, no need of work array
             call gesv(min(diis%idsx,diis%ids)+1,1,adsw(1,1),diis%idsx+1,  & 
                  ipiv(1),rds(1,igrp),diis%idsx+1,info)
             if (info /= 0) then
                call f_err_throw('Error in GESV operation, info='//trim(yaml_toa(info))//&
                     ' Size='//trim(yaml_toa(min(diis%idsx,diis%ids)+1)),BIGDFT_LINALG_ERROR)
                call free_and_exit()
                return
             end if
          else
             rds(1,igrp)=1.0_tp
          endif

          !change the approach and fill only the difference between the original psit and the updated one
          jst=max(1,diis%ids-diis%idsx+1)
          !use the array which will be erased in the next step as the work array
          iacc_add=ispsidst+(mod(diis%ids,diis%idsx))*ncomp
          if (diis%ids < diis%idsx) then
             call f_zero(ncomp,diis%psidst(iacc_add))
          end if

          jj=0
          do j=jst,diis%ids
             jj=jj+1
             mj=mod(j-1,diis%idsx)+1

             !correct the coefficient for the daxpy in psi for the first and the last cycle
             if ((j==jst .and. diis%ids >= diis%idsx) .or. j==diis%ids) then
                psicoeff=rds(jj,igrp)-1.0_tp
             else
                psicoeff=rds(jj,igrp)
             end if
             !use axpy for updating the array (can be done for all the orbitals in the group)
             !the last step is the update with psi
             call axpy(ncomp,psicoeff,&
                  diis%psidst(ispsidst+(mj-1)*ncomp),1,&
                  diis%psidst(iacc_add),1)
             !this will work only if the errors are written in the same precision
             call axpy(ncomp,-rds(jj,igrp),&
                  diis%hpsidst(ispsidst+(mj-1)*ncomp),1,&
                  diis%psidst(iacc_add),1)
          end do
       ispsi=ispsi+ncomp
       ispsidst=ispsidst+ncomp*diis%idsx
    end do

    ! Output to screen, depending on policy.
    if (verbose >= 10) then
       call broadcast_kpt_objects(nproc, ngrp, (diis%idsx+1), rds, igrpproc)
    end if
    if (iproc == 0) then 
       call write_diis_weights(1,diis%idsx,1,ngrp,min(diis%idsx,diis%ids),rds)
    endif

    call free_and_exit()

    contains

      subroutine free_and_exit()
        implicit none
        
        call f_free(ipiv)
        call f_free(rds)
        call f_free(adsw)
        call f_release_routine()

      end subroutine free_and_exit

  END SUBROUTINE diis_step

  

  !> Temporary routine to test the diis_step procedure
  subroutine DIIS_obj_fill(diis_old,diis)
    use module_types
    implicit none
    type(diis_objects), intent(in) :: diis_old
    type(DIIS_obj), intent(out) :: diis
    
    !fill the integers
    diis%ids  =diis_old%ids 
    diis%mids =diis_old%mids
    diis%idsx =diis_old%idsx
    
    diis%ads=f_malloc_ptr((/diis%idsx+1,diis%idsx+1,size(diis_old%ads,5)/),id='ads')
    
    call vcopy(product(shape(diis%ads)),diis_old%ads(1,1,1,1,1,1),1,diis%ads(1,1,1),1)

    !associate pointers
    diis%psidst=>diis_old%psidst
    diis%hpsidst=>diis_old%hpsidst

  end subroutine DIIS_obj_fill


  subroutine DIIS_obj_release(diis,diis_old)
    use module_types
    implicit none
    type(DIIS_obj), intent(inout) :: diis
    type(diis_objects), intent(inout) :: diis_old
    
    !Fill the integers
    diis_old%ids  =diis%ids 
    diis_old%mids =diis%mids
    diis_old%idsx =diis%idsx
       
    call vcopy(product(shape(diis%ads)),diis%ads(1,1,1),1,diis_old%ads(1,1,1,1,1,1),1)

    !Nullify pointers
    nullify(diis%psidst)
    nullify(diis%hpsidst)

    call f_free_ptr(diis%ads)

  end subroutine DIIS_obj_release

end module diis_sd_optimization


subroutine diis_opt(iproc,nproc,ngrp,isgrp,ngrpp,igrpproc,ncomp_grp,ndim_psi,psi,hpsi,diis)
    use module_defs, only: wp
    use diis_sd_optimization
    implicit none
    ! Arguments
    integer, intent(in) :: nproc,iproc,ngrp,isgrp,ngrpp
    integer, dimension(ngrp), intent(in) :: igrpproc  !< Array which associate each group to only one iproc for broadcasting
    integer, dimension(ngrp), intent(in) :: ncomp_grp !< Number of components per group
    type(DIIS_obj), intent(inout) :: diis
    integer, intent(in) :: ndim_psi                   !< Should be greater or equal to sum(ncomp_grp(isgrp+1:isgrp+ngrpp)
    real(wp), dimension(ndim_psi), intent(inout) :: psi
    real(wp), dimension(ndim_psi), intent(in) :: hpsi

    call DIIS_update_errors(ngrp,isgrp,ngrpp,ncomp_grp,ndim_psi,psi,hpsi,diis)

    call diis_step(iproc,nproc,ngrp,isgrp,ngrpp,igrpproc,ncomp_grp,diis)

    call DIIS_update_psi(ngrp,isgrp,ngrpp,ncomp_grp,ndim_psi,psi,diis)

end subroutine diis_opt
