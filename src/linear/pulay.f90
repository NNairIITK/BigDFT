!> @file
!! Pulay correction calculation for linear version
!! @author
!!    Copyright (C) 2007-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


subroutine pulay_correction_new(iproc, nproc, tmb, orbs, at, fpulay)
  use module_base
  use module_types
  use module_interfaces, except_this_one => pulay_correction_new
  use yaml_output
  use communications_base, only: TRANSPOSE_FULL
  use communications, only: transpose_localized
  use sparsematrix_base, only: sparsematrix_malloc_ptr, DENSE_FULL, assignment(=), &
                               sparsematrix_malloc, SPARSE_FULL
  use sparsematrix, only: compress_matrix, uncompress_matrix, gather_matrix_from_taskgroups_inplace
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(DFT_wavefunction),intent(inout) :: tmb
  type(orbitals_data),intent(in) :: orbs
  type(atoms_data),intent(in) :: at
  real(kind=8),dimension(3,at%astruct%nat),intent(out) :: fpulay

  ! Local variables
  integer :: iat, isize, iorb, jorb, korb, idir, iiorb, ierr, num_points, num_points_tot
  real(kind=8),dimension(:,:),allocatable :: phi_delta, energykernel, tempmat, phi_delta_large
  real(kind=8),dimension(:),allocatable :: hphit_c, hphit_f, denskern_tmp, delta_phit_c, delta_phit_f, tmparr
  real(kind=8) :: tt

  call timing(iproc,'new_pulay_corr','ON') 

  call f_routine(id='pulay_correction_new')

  phi_delta=f_malloc0((/tmb%npsidim_orbs,3/),id='phi_delta')
  ! Get the values of the support functions on the boundary of the localization region
  call extract_boundary(tmb, phi_delta, num_points, num_points_tot)


  ! calculate the "energy kernel"
  energykernel=f_malloc0((/tmb%orbs%norb,tmb%orbs%norb/),id='energykernel')
  do iorb=1,tmb%orbs%norbp
      iiorb=tmb%orbs%isorb+iorb
      do jorb=1,tmb%orbs%norb
          do korb=1,orbs%norb
              energykernel(jorb,iiorb) = energykernel(jorb,iiorb) &
                                      + tmb%coeff(jorb,korb)*tmb%coeff(iiorb,korb)*tmb%orbs%eval(korb)
          end do
      end do
  end do

  if (nproc > 1) then
     call mpiallred(energykernel(1,1), tmb%orbs%norb**2, mpi_sum, bigdft_mpi%mpi_comm)
  end if

  ! calculate the overlap matrix
  !!if(.not.associated(tmb%psit_c)) then
  !!    isize=sum(tmb%collcom%nrecvcounts_c)
  !!    tmb%psit_c=f_malloc_ptr(isize,id='tmb%psit_c')
  !!end if
  !!if(.not.associated(tmb%psit_f)) then
  !!    isize=7*sum(tmb%collcom%nrecvcounts_f)
  !!    tmb%psit_f=f_malloc_ptr(isize,id=' tmb%psit_f')
  !!end if
  call transpose_localized(iproc, nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
       TRANSPOSE_FULL, tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd)
  call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%collcom, tmb%psit_c, &
       tmb%psit_c, tmb%psit_f, tmb%psit_f, tmb%linmat%s, tmb%linmat%ovrlp_)
  call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%s, tmb%linmat%ovrlp_)
  ! This can then be deleted if the transition to the new type has been completed.
  !tmb%linmat%ovrlp%matrix_compr=tmb%linmat%ovrlp_%matrix_compr


  !!call f_free_ptr(tmb%psit_c)
  !!call f_free_ptr(tmb%psit_f)


  ! Construct the array chi
  call construct_chi()
  call f_free(energykernel)


  ! transform phi_delta to the shamop region
  phi_delta_large=f_malloc((/tmb%ham_descr%npsidim_orbs,3/),id='phi_delta_large')
  call small_to_large_locreg(iproc, tmb%npsidim_orbs, tmb%ham_descr%npsidim_orbs, tmb%lzd, tmb%ham_descr%lzd, &
                      tmb%orbs, phi_delta(1,1), phi_delta_large(1,1))
  call small_to_large_locreg(iproc, tmb%npsidim_orbs, tmb%ham_descr%npsidim_orbs, tmb%lzd, tmb%ham_descr%lzd, &
                      tmb%orbs, phi_delta(1,2), phi_delta_large(1,2))
  call small_to_large_locreg(iproc, tmb%npsidim_orbs, tmb%ham_descr%npsidim_orbs, tmb%lzd, tmb%ham_descr%lzd, &
                      tmb%orbs, phi_delta(1,3), phi_delta_large(1,3))

  isize=sum(tmb%ham_descr%collcom%nrecvcounts_c)
  delta_phit_c=f_malloc(isize,id='delta_phit_c')
  isize=7*sum(tmb%ham_descr%collcom%nrecvcounts_f)
  delta_phit_f=f_malloc(isize,id='delta_phit_f')
  !fpulay=f_malloc((/3,at%astruct%nat/),id='fpulay')
  tmb%linmat%kernel_%matrix = sparsematrix_malloc_ptr(tmb%linmat%l, iaction=DENSE_FULL, id='tmb%linmat%kernel_%matrix')
  call uncompress_matrix(iproc,tmb%linmat%l, inmat=tmb%linmat%kernel_%matrix_compr, outmat=tmb%linmat%kernel_%matrix)
  call to_zero(3*at%astruct%nat, fpulay(1,1))
  do idir=1,3
      ! calculate the overlap matrix among hphi and phi_delta_large
      call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
           TRANSPOSE_FULL, phi_delta_large(1,idir), delta_phit_c, delta_phit_f, tmb%ham_descr%lzd)
      call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%ham_descr%collcom, &
           hphit_c, delta_phit_c, hphit_f, delta_phit_f, tmb%linmat%m, tmb%linmat%ham_)
      call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%m, tmb%linmat%ham_)
      ! This can then be deleted if the transition to the new type has been completed.
      !tmb%linmat%ham%matrix_compr=tmb%linmat%ham_%matrix_compr

      tmb%linmat%ham_%matrix = sparsematrix_malloc_ptr(tmb%linmat%m, iaction=DENSE_FULL, id='tmb%linmat%ham_%matrix')
      call uncompress_matrix(iproc, tmb%linmat%m, &
           inmat=tmb%linmat%ham_%matrix_compr, outmat=tmb%linmat%ham_%matrix)

      do iorb=1,tmb%orbs%norbp
          iiorb=tmb%orbs%isorb+iorb
          iat=tmb%orbs%onwhichatom(iiorb)
          tt=0.d0
          do jorb=1,tmb%orbs%norb
              tt = tt -2.d0*tmb%linmat%kernel_%matrix(jorb,iiorb,1)*tmb%linmat%ham_%matrix(jorb,iiorb,1)
              !if (iproc==0) write(*,*) 'kern, ovrlp', tmb%linmat%denskern%matrix(jorb,iiorb), tmb%linmat%ham%matrix(iiorb,jorb)
          end do  
          fpulay(idir,iat)=fpulay(idir,iat)+tt
      end do
      call f_free_ptr(tmb%linmat%ham_%matrix)
  end do

  if (nproc > 1) then
     call mpiallred(fpulay(1,1), 3*at%astruct%nat, mpi_sum, bigdft_mpi%mpi_comm)
  end if
  call f_free_ptr(tmb%linmat%kernel_%matrix)


  if(iproc==0) then
       call yaml_comment('new Pulay correction',hfill='-')
       call yaml_sequence_open('Pulay forces (Ha/Bohr)')
          do iat=1,at%astruct%nat
             call yaml_sequence(advance='no')
             call yaml_mapping_open(flow=.true.)
             call yaml_map(trim(at%astruct%atomnames(at%astruct%iatype(iat))),fpulay(1:3,iat),fmt='(1es20.12)')
             call yaml_mapping_close(advance='no')
             call yaml_comment(trim(yaml_toa(iat,fmt='(i4.4)')))
          end do
          call yaml_sequence_close()
  end if


    
  call f_free(phi_delta)
  call f_free(phi_delta_large)
  call f_free(hphit_c)
  call f_free(hphit_f)
  call f_free(delta_phit_c)
  call f_free(delta_phit_f)
  
  call f_release_routine()

  call timing(iproc,'new_pulay_corr','OF') 

  contains

    subroutine construct_chi()

      tempmat=f_malloc((/tmb%orbs%norb,tmb%orbs%norb/),id='tempmat')
      tmb%linmat%ovrlp_%matrix = sparsematrix_malloc_ptr(tmb%linmat%s, iaction=DENSE_FULL, id='tmb%linmat%ovrlp_%matrix')
      call uncompress_matrix(iproc, tmb%linmat%s, inmat=tmb%linmat%ovrlp_%matrix_compr, outmat=tmb%linmat%ovrlp_%matrix)
      call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norb, tmb%orbs%norb, 1.d0, &
                 tmb%linmat%ovrlp_%matrix, tmb%orbs%norb, energykernel, tmb%orbs%norb, &
                 0.d0, tempmat, tmb%orbs%norb)
      call f_free_ptr(tmb%linmat%ovrlp_%matrix)
      isize=sum(tmb%ham_descr%collcom%nrecvcounts_c)
      hphit_c=f_malloc(isize,id='hphit_c')
      isize=7*sum(tmb%ham_descr%collcom%nrecvcounts_f)
      hphit_f=f_malloc(isize,id='hphit_f')
      call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
                               TRANSPOSE_FULL, tmb%hpsi, hphit_c, hphit_f, tmb%ham_descr%lzd)
      denskern_tmp=f_malloc(tmb%linmat%l%nvctr,id='denskern_tmp')
      denskern_tmp=tmb%linmat%kernel_%matrix_compr
      tmb%linmat%kernel_%matrix = sparsematrix_malloc_ptr(tmb%linmat%l, iaction=DENSE_FULL, id='tmb%linmat%kernel_%matrix')
      tmb%linmat%kernel_%matrix(:,:,1)=tempmat
      call compress_matrix(iproc, tmb%linmat%l, inmat=tmb%linmat%kernel_%matrix, outmat=tmb%linmat%kernel_%matrix_compr)
      call f_free_ptr(tmb%linmat%kernel_%matrix)

      tmparr = sparsematrix_malloc(tmb%linmat%l,iaction=SPARSE_FULL,id='tmparr')
      call vcopy(tmb%linmat%l%nvctr, tmb%linmat%kernel_%matrix_compr(1), 1, tmparr(1), 1)
      call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%l, tmb%linmat%kernel_)
      call build_linear_combination_transposed(tmb%ham_descr%collcom, &
           tmb%linmat%l, tmb%linmat%kernel_, tmb%ham_descr%psit_c, tmb%ham_descr%psit_f, &
           .false., hphit_c, hphit_f, iproc)
      call vcopy(tmb%linmat%l%nvctr, tmparr(1), 1, tmb%linmat%kernel_%matrix_compr(1), 1)
      call f_free(tmparr)

      tmb%linmat%kernel_%matrix_compr=denskern_tmp
    
      call f_free(tempmat)
      call f_free(denskern_tmp)
    
    end subroutine construct_chi
 
end subroutine pulay_correction_new




subroutine extract_boundary(tmb, phi_delta, numpoints, numpoints_tot)
  use module_base
  use module_types
  use module_interfaces
  use yaml_output
  implicit none

  ! Calling arguments
  type(DFT_wavefunction),intent(in) :: tmb
  real(kind=8),dimension(tmb%npsidim_orbs,3),intent(out) :: phi_delta
  integer, intent(out) :: numpoints, numpoints_tot

  ! Local variables
  integer :: ishift, iorb, iiorb, ilr, iseg, jj_prev, j0_prev, j1_prev, ii_prev, i3_prev, i2_prev, i1_prev, i0_prev
  integer :: jj, j0, j1, ii, i3, i2, i1, i0, jorb, korb, istat, idir, iat, i, isize
  logical,dimension(:,:,:),allocatable :: boundaryarray
  real(kind=8) :: dist, crit, xsign, ysign, zsign

  call f_routine(id='extract_boundary')
  numpoints=0
  ! First copy the boundary elements of the first array to a temporary array,
  ! filling the remaining part with zeros.
  call to_zero(3*tmb%npsidim_orbs,phi_delta(1,1))
  ishift=0
  do iorb=1,tmb%orbs%norbp
      iiorb=tmb%orbs%isorb+iorb
      ilr=tmb%orbs%inwhichlocreg(iiorb)
      boundaryarray=f_malloc((/0.to.tmb%lzd%llr(ilr)%d%n1,0.to.tmb%lzd%llr(ilr)%d%n2,0.to.tmb%lzd%llr(ilr)%d%n3/), &
                             id='boundaryarray')
      boundaryarray=.false.


      ! coarse part
      do iseg=1,tmb%lzd%llr(ilr)%wfd%nseg_c

          ! indizes of the last element of the previous segment
          if (iseg>1) then
              jj_prev=tmb%lzd%llr(ilr)%wfd%keyvloc(iseg-1)
              j0_prev=tmb%lzd%llr(ilr)%wfd%keygloc(1,iseg-1)
              j1_prev=tmb%lzd%llr(ilr)%wfd%keygloc(2,iseg-1)
              ii_prev=j0_prev-1
              i3_prev=ii_prev/((tmb%lzd%llr(ilr)%d%n1+1)*(tmb%lzd%llr(ilr)%d%n2+1))
              ii_prev=ii_prev-i3_prev*(tmb%lzd%llr(ilr)%d%n1+1)*(tmb%lzd%llr(ilr)%d%n2+1)
              i2_prev=ii_prev/(tmb%lzd%llr(ilr)%d%n1+1)
              i0_prev=ii_prev-i2_prev*(tmb%lzd%llr(ilr)%d%n1+1)
              i1_prev=i0_prev+j1_prev-j0_prev
          else
              !just some large values
              i1_prev=-99999
              i2_prev=-99999
              i3_prev=-99999
          end if

          ! indizes of the first element of the current segment
          jj=tmb%lzd%llr(ilr)%wfd%keyvloc(iseg)
          j0=tmb%lzd%llr(ilr)%wfd%keygloc(1,iseg)
          j1=tmb%lzd%llr(ilr)%wfd%keygloc(2,iseg)
          ii=j0-1
          i3=ii/((tmb%lzd%llr(ilr)%d%n1+1)*(tmb%lzd%llr(ilr)%d%n2+1))
          ii=ii-i3*(tmb%lzd%llr(ilr)%d%n1+1)*(tmb%lzd%llr(ilr)%d%n2+1)
          i2=ii/(tmb%lzd%llr(ilr)%d%n1+1)
          i0=ii-i2*(tmb%lzd%llr(ilr)%d%n1+1)
          i1=i0

          if (i2/=i2_prev .or. i3/=i3_prev) then
              ! Segment starts on a new line, i.e. this is a boundary element.
              ! Furthermore the last element of the previous segment must be a
              ! boundary element as well. However this is only true if the
              ! distance from the locreg center corresponds to the cutoff
              ! radius, otherwise it is only a boundary element of the global
              ! grid and not of the localization region.
              if (iseg>1) then
                  ! this makes only sense if we are not in the first segment
                  dist= sqrt(((tmb%lzd%llr(ilr)%ns1+i1_prev)*tmb%lzd%hgrids(1)-tmb%lzd%llr(ilr)%locregcenter(1))**2 &
                            +((tmb%lzd%llr(ilr)%ns2+i2_prev)*tmb%lzd%hgrids(2)-tmb%lzd%llr(ilr)%locregcenter(2))**2 &
                            +((tmb%lzd%llr(ilr)%ns3+i3_prev)*tmb%lzd%hgrids(3)-tmb%lzd%llr(ilr)%locregcenter(3))**2)
                  crit=tmb%lzd%llr(ilr)%locrad-sqrt(tmb%lzd%hgrids(1)**2+tmb%lzd%hgrids(2)**2+tmb%lzd%hgrids(3)**2)
                  if (dist>=crit) then
                      ! boundary element of the locreg
                      xsign=((tmb%lzd%llr(ilr)%ns1+i1_prev)*tmb%lzd%hgrids(1)-tmb%lzd%llr(ilr)%locregcenter(1))/dist
                      ysign=((tmb%lzd%llr(ilr)%ns2+i2_prev)*tmb%lzd%hgrids(2)-tmb%lzd%llr(ilr)%locregcenter(2))/dist
                      zsign=((tmb%lzd%llr(ilr)%ns3+i3_prev)*tmb%lzd%hgrids(3)-tmb%lzd%llr(ilr)%locregcenter(3))/dist
                      !xsign=1.d0 ; ysign=1.d0 ; zsign=1.d0
                      phi_delta(ishift+jj-1,1)=xsign*tmb%psi(ishift+jj-1)
                      phi_delta(ishift+jj-1,2)=ysign*tmb%psi(ishift+jj-1)
                      phi_delta(ishift+jj-1,3)=zsign*tmb%psi(ishift+jj-1)
                      boundaryarray(i1_prev,i2_prev,i3_prev)=.true.
                      numpoints=numpoints+1
                  end if
                  numpoints_tot=numpoints_tot+1
              end if
              dist= sqrt(((tmb%lzd%llr(ilr)%ns1+i1)*tmb%lzd%hgrids(1)-tmb%lzd%llr(ilr)%locregcenter(1))**2 &
                        +((tmb%lzd%llr(ilr)%ns2+i2)*tmb%lzd%hgrids(2)-tmb%lzd%llr(ilr)%locregcenter(2))**2 &
                        +((tmb%lzd%llr(ilr)%ns3+i3)*tmb%lzd%hgrids(3)-tmb%lzd%llr(ilr)%locregcenter(3))**2)
              crit=tmb%lzd%llr(ilr)%locrad-sqrt(tmb%lzd%hgrids(1)**2+tmb%lzd%hgrids(2)**2+tmb%lzd%hgrids(3)**2)
              if (dist>=crit) then
                  xsign=((tmb%lzd%llr(ilr)%ns1+i1)*tmb%lzd%hgrids(1)-tmb%lzd%llr(ilr)%locregcenter(1))/dist
                  ysign=((tmb%lzd%llr(ilr)%ns2+i2)*tmb%lzd%hgrids(2)-tmb%lzd%llr(ilr)%locregcenter(2))/dist
                  zsign=((tmb%lzd%llr(ilr)%ns3+i3)*tmb%lzd%hgrids(3)-tmb%lzd%llr(ilr)%locregcenter(3))/dist
                  !xsign=1.d0 ; ysign=1.d0 ; zsign=1.d0
                  phi_delta(ishift+jj,1)=xsign*tmb%psi(ishift+jj)
                  phi_delta(ishift+jj,2)=ysign*tmb%psi(ishift+jj)
                  phi_delta(ishift+jj,3)=zsign*tmb%psi(ishift+jj)
                  boundaryarray(i1,i2,i3)=.true.
                  numpoints=numpoints+1
              end if
              numpoints_tot=numpoints_tot+1
          end if
      end do

      ! fine part
      do iseg=tmb%lzd%llr(ilr)%wfd%nseg_c+1,tmb%lzd%llr(ilr)%wfd%nseg_c+tmb%lzd%llr(ilr)%wfd%nseg_f
          jj=tmb%lzd%llr(ilr)%wfd%keyvloc(iseg)
          j0=tmb%lzd%llr(ilr)%wfd%keygloc(1,iseg)
          j1=tmb%lzd%llr(ilr)%wfd%keygloc(2,iseg)
          ii=j0-1
          i3=ii/((tmb%lzd%llr(ilr)%d%n1+1)*(tmb%lzd%llr(ilr)%d%n2+1))
          ii=ii-i3*(tmb%lzd%llr(ilr)%d%n1+1)*(tmb%lzd%llr(ilr)%d%n2+1)
          i2=ii/(tmb%lzd%llr(ilr)%d%n1+1)
          i0=ii-i2*(tmb%lzd%llr(ilr)%d%n1+1)
          i1=i0+j1-j0
          ! The segments goes now from i0 to i1
          ! Check the beginnig of the segment. If it was a boundary element of
          ! the coarse grid, copy its content also for the fine part.
          if (boundaryarray(i0,i2,i3)) then
              dist= sqrt(((tmb%lzd%llr(ilr)%ns1+i0)*tmb%lzd%hgrids(1)-tmb%lzd%llr(ilr)%locregcenter(1))**2 &
                        +((tmb%lzd%llr(ilr)%ns2+i2)*tmb%lzd%hgrids(2)-tmb%lzd%llr(ilr)%locregcenter(2))**2 &
                        +((tmb%lzd%llr(ilr)%ns3+i3)*tmb%lzd%hgrids(3)-tmb%lzd%llr(ilr)%locregcenter(3))**2)
              xsign=((tmb%lzd%llr(ilr)%ns1+i0)*tmb%lzd%hgrids(1)-tmb%lzd%llr(ilr)%locregcenter(1))/dist
              ysign=((tmb%lzd%llr(ilr)%ns2+i2)*tmb%lzd%hgrids(2)-tmb%lzd%llr(ilr)%locregcenter(2))/dist
              zsign=((tmb%lzd%llr(ilr)%ns3+i3)*tmb%lzd%hgrids(3)-tmb%lzd%llr(ilr)%locregcenter(3))/dist
              !xsign=1.d0 ; ysign=1.d0 ; zsign=1.d0
              ! x direction
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+1,1) = &
                  xsign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+1)
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+2,1) = &
                  xsign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+2)
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+3,1) = &
                  xsign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+3)
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+4,1) = &
                  xsign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+4)
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+5,1) = &
                  xsign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+5)
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+6,1) = &
                  xsign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+6)
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+7,1) = &
                  xsign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+7)
              ! y direction
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+1,2) = &
                  ysign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+1)
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+2,2) = &
                  ysign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+2)
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+3,2) = &
                  ysign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+3)
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+4,2) = &
                  ysign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+4)
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+5,2) = &
                  ysign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+5)
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+6,2) = &
                  ysign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+6)
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+7,2) = &
                  ysign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+7)
              ! z direction
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+1,3) = &
                  zsign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+1)
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+2,3) = &
                  zsign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+2)
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+3,3) = &
                  zsign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+3)
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+4,3) = &
                  zsign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+4)
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+5,3) = &
                  zsign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+5)
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+6,3) = &
                  zsign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+6)
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+7,3) = &
                  zsign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(jj-1)+7)
              numpoints=numpoints+7
          end if
          numpoints_tot=numpoints_tot+1
          ! Check the end of the segment. If it was a boundary element of
          ! the coarse grid, copy its content also for the fine part.
          if (boundaryarray(i1,i2,i3)) then
              dist= sqrt(((tmb%lzd%llr(ilr)%ns1+i1)*tmb%lzd%hgrids(1)-tmb%lzd%llr(ilr)%locregcenter(1))**2 &
                        +((tmb%lzd%llr(ilr)%ns2+i2)*tmb%lzd%hgrids(2)-tmb%lzd%llr(ilr)%locregcenter(2))**2 &
                        +((tmb%lzd%llr(ilr)%ns3+i3)*tmb%lzd%hgrids(3)-tmb%lzd%llr(ilr)%locregcenter(3))**2)
              xsign=((tmb%lzd%llr(ilr)%ns1+i1)*tmb%lzd%hgrids(1)-tmb%lzd%llr(ilr)%locregcenter(1))/dist
              ysign=((tmb%lzd%llr(ilr)%ns2+i2)*tmb%lzd%hgrids(2)-tmb%lzd%llr(ilr)%locregcenter(2))/dist
              zsign=((tmb%lzd%llr(ilr)%ns3+i3)*tmb%lzd%hgrids(3)-tmb%lzd%llr(ilr)%locregcenter(3))/dist
              !xsign=1.d0 ; ysign=1.d0 ; zsign=1.d0
              ! x direction
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+1,1) = &
                  xsign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+1)
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+2,1) = &
                  xsign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+2)
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+3,1) = &
                  xsign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+3)
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+4,1) = &
                  xsign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+4)
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+5,1) = &
                  xsign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+5)
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+6,1) = &
                  xsign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+6)
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+7,1) = &
                  xsign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+7)
              ! y direction
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+1,2) = &
                  ysign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+1)
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+2,2) = &
                  ysign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+2)
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+3,2) = &
                  ysign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+3)
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+4,2) = &
                  ysign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+4)
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+5,2) = &
                  ysign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+5)
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+6,2) = &
                  ysign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+6)
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+7,2) = &
                  ysign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+7)
              ! z direction
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+1,3) = &
                  zsign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+1)
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+2,3) = &
                  zsign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+2)
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+3,3) = &
                  zsign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+3)
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+4,3) = &
                  zsign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+4)
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+5,3) = &
                  zsign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+5)
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+6,3) = &
                  zsign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+6)
              phi_delta(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+7,3) = &
                  zsign*tmb%psi(ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*(i1-i0+jj-1)+7)
              numpoints=numpoints+7
          end if
          numpoints_tot=numpoints_tot+1
      end do
      
      ishift=ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
      call f_free(boundaryarray)
  end do

call f_release_routine()

end subroutine extract_boundary


subroutine pulay_correction(iproc, nproc, orbs, at, rxyz, nlpsp, SIC, denspot, GPU, tmb, fpulay)
  use module_base
  use module_types
  use module_interfaces, except_this_one => pulay_correction
  use yaml_output
  use communications_base, only: TRANSPOSE_FULL
  use communications, only: transpose_localized, start_onesided_communication
  use sparsematrix_base, only: sparse_matrix, sparse_matrix_null, deallocate_sparse_matrix, &
                               matrices_null, allocate_matrices, deallocate_matrices, &
                               sparsematrix_malloc_ptr, SPARSE_FULL, assignment(=)
  use sparsematrix, only: gather_matrix_from_taskgroups_inplace
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(orbitals_data),intent(in) :: orbs
  type(atoms_data),intent(in) :: at
  real(kind=8),dimension(3,at%astruct%nat),intent(in) :: rxyz
  type(DFT_PSP_projectors),intent(inout) :: nlpsp
  type(SIC_data),intent(in) :: SIC
  type(DFT_local_fields), intent(inout) :: denspot
  type(GPU_pointers),intent(inout) :: GPU
  type(DFT_wavefunction),intent(inout) :: tmb
  real(kind=8),dimension(3,at%astruct%nat),intent(out) :: fpulay

  ! Local variables
  integer:: istat, iall, ierr, iialpha, jorb
  integer:: iorb, ii, iseg, isegstart, isegend, is, ie, ishift, ispin
  integer:: jat, jdir, ibeta
  !!integer :: ialpha, iat, iiorb
  real(kind=8) :: kernel, ekernel
  real(kind=8),dimension(:),allocatable :: lhphilarge, psit_c, psit_f, hpsit_c, hpsit_f, lpsit_c, lpsit_f
  type(sparse_matrix) :: dovrlp(3), dham(3)
  type(matrices) :: dovrlp_(3), dham_(3)
  type(energy_terms) :: energs
  type(confpot_data),dimension(:),allocatable :: confdatarrtmp
  character(len=*),parameter :: subname='pulay_correction'
  type(matrices) :: ham_

  call f_routine(id='pulay_correction')
  energs=energy_terms_null()
  
  ! Begin by updating the Hpsi
  call local_potential_dimensions(iproc,tmb%ham_descr%lzd,tmb%orbs,denspot%xc,denspot%dpbox%ngatherarr(0,1))

  lhphilarge = f_malloc0(tmb%ham_descr%npsidim_orbs,id='lhphilarge')

  !!call post_p2p_communication(iproc, nproc, denspot%dpbox%ndimpot, denspot%rhov, &
  !!     tmb%ham_descr%comgp%nrecvbuf, tmb%ham_descr%comgp%recvbuf, tmb%ham_descr%comgp, tmb%ham_descr%lzd)
  call start_onesided_communication(iproc, nproc, denspot%dpbox%ndims(1), denspot%dpbox%ndims(2), &
       max(denspot%dpbox%nscatterarr(:,2),1), denspot%rhov, &
       tmb%ham_descr%comgp%nspin*tmb%ham_descr%comgp%nrecvbuf, tmb%ham_descr%comgp%recvbuf, tmb%ham_descr%comgp, tmb%ham_descr%lzd)

  allocate(confdatarrtmp(tmb%orbs%norbp))
  call default_confinement_data(confdatarrtmp,tmb%orbs%norbp)


  call NonLocalHamiltonianApplication(iproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,&
       tmb%ham_descr%lzd,nlpsp,tmb%ham_descr%psi,lhphilarge,energs%eproj,tmb%paw)

  ! only kinetic because waiting for communications
  call LocalHamiltonianApplication(iproc,nproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,&
       tmb%ham_descr%lzd,confdatarrtmp,denspot%dpbox%ngatherarr,denspot%pot_work,tmb%ham_descr%psi,lhphilarge,&
       energs,SIC,GPU,3,denspot%xc,pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,&
       & potential=denspot%rhov,comgp=tmb%ham_descr%comgp)
  call full_local_potential(iproc,nproc,tmb%orbs,tmb%ham_descr%lzd,2,denspot%dpbox,&
       & denspot%xc,denspot%rhov,denspot%pot_work,tmb%ham_descr%comgp)
  ! only potential
  call LocalHamiltonianApplication(iproc,nproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,&
       tmb%ham_descr%lzd,confdatarrtmp,denspot%dpbox%ngatherarr,denspot%pot_work,tmb%ham_descr%psi,lhphilarge,&
       energs,SIC,GPU,2,denspot%xc,pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,&
       & potential=denspot%rhov,comgp=tmb%ham_descr%comgp)

  call timing(iproc,'glsynchham1','ON') !lr408t
  call SynchronizeHamiltonianApplication(nproc,tmb%ham_descr%npsidim_orbs,&
       & tmb%orbs,tmb%ham_descr%lzd,GPU,denspot%xc,lhphilarge,&
       energs%ekin,energs%epot,energs%eproj,energs%evsic,energs%eexctX)
  call timing(iproc,'glsynchham1','OF') !lr408t
  deallocate(confdatarrtmp)
  

  ! Now transpose the psi and hpsi
  lpsit_c = f_malloc(tmb%ham_descr%collcom%ndimind_c,id='lpsit_c')
  lpsit_f = f_malloc(7*tmb%ham_descr%collcom%ndimind_f,id='lpsit_f')
  hpsit_c = f_malloc(tmb%ham_descr%collcom%ndimind_c,id='hpsit_c')
  hpsit_f = f_malloc(7*tmb%ham_descr%collcom%ndimind_f,id='hpsit_f')
  psit_c = f_malloc(tmb%ham_descr%collcom%ndimind_c,id='psit_c')
  psit_f = f_malloc(7*tmb%ham_descr%collcom%ndimind_f,id='psit_f')


  call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
       TRANSPOSE_FULL, tmb%ham_descr%psi, lpsit_c, lpsit_f, tmb%ham_descr%lzd)

  call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
       TRANSPOSE_FULL, lhphilarge, hpsit_c, hpsit_f, tmb%ham_descr%lzd)

  !now build the derivative and related matrices <dPhi_a | H | Phi_b> and <dPhi_a | Phi_b>

  ! DOVRLP AND DHAM SHOULD HAVE DIFFERENT SPARSITIES, BUT TO MAKE LIFE EASIER KEEPING THEM THE SAME FOR NOW
  ! also array of structure a bit inelegant at the moment

  ham_ = matrices_null()
  call allocate_matrices(tmb%linmat%m, allocate_full=.false., &
       matname='ham_', mat=ham_)
  do jdir = 1, 3
    !call nullify_sparse_matrix(dovrlp(jdir))
    !call nullify_sparse_matrix(dham(jdir))
    dovrlp(jdir)=sparse_matrix_null()
    dovrlp_(jdir)=matrices_null()
    dham(jdir)=sparse_matrix_null()
    dham_(jdir)=matrices_null()
    !call sparse_copy_pattern(tmb%linmat%m,dovrlp(jdir),iproc,subname) 
    !call sparse_copy_pattern(tmb%linmat%m,dham(jdir),iproc,subname)
    call copy_sparse_matrix(tmb%linmat%m,dovrlp(jdir))
    call copy_sparse_matrix(tmb%linmat%m,dham(jdir))
    !dham_(jdir)%matrix_compr=f_malloc_ptr(dham(jdir)%nvctr,id='dham(jdir)%matrix_compr')
    !dovrlp_(jdir)%matrix_compr=f_malloc_ptr(dovrlp(jdir)%nvctr,id='dovrlp(jdir)%matrix_compr')
    dham_(jdir)%matrix_compr=sparsematrix_malloc_ptr(dham(jdir),iaction=SPARSE_FULL,id='dham(jdir)%matrix_compr')
    dovrlp_(jdir)%matrix_compr=sparsematrix_malloc_ptr(dovrlp(jdir),iaction=SPARSE_FULL,id='dovrlp(jdir)%matrix_compr')

    call get_derivative(jdir, tmb%ham_descr%npsidim_orbs, tmb%ham_descr%lzd%hgrids(1), tmb%orbs, &
         tmb%ham_descr%lzd, tmb%ham_descr%psi, lhphilarge)

    call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
         TRANSPOSE_FULL, lhphilarge, psit_c, psit_f, tmb%ham_descr%lzd)

    call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%ham_descr%collcom,&
         psit_c, lpsit_c, psit_f, lpsit_f, tmb%linmat%m, ham_)
    call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%m, tmb%linmat%ham_)
    ! This can then be deleted if the transition to the new type has been completed.
    dovrlp_(jdir)%matrix_compr=ham_%matrix_compr

    call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%ham_descr%collcom,&
         psit_c, hpsit_c, psit_f, hpsit_f, tmb%linmat%m, ham_)
    call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%m, tmb%linmat%ham_)
    ! This can then be deleted if the transition to the new type has been completed.
    dham_(jdir)%matrix_compr=ham_%matrix_compr
  end do
  call deallocate_matrices(ham_)


  !DEBUG
  !!print *,'iproc,tmb%orbs%norbp',iproc,tmb%orbs%norbp
  !!if(iproc==0)then
  !!do iorb = 1, tmb%orbs%norb
  !!   do iiorb=1,tmb%orbs%norb
  !!      !print *,'Hamiltonian of derivative: ',iorb, iiorb, (matrix(iorb,iiorb,jdir),jdir=1,3)
  !!      print *,'Overlap of derivative: ',iorb, iiorb, (dovrlp(iorb,iiorb,jdir),jdir=1,3)
  !!   end do
  !!end do
  !!end if
  !!!Check if derivatives are orthogonal to functions
  !!if(iproc==0)then
  !!  do iorb = 1, tmbder%orbs%norb
  !!     !print *,'overlap of derivative: ',iorb, (dovrlp(iorb,iiorb),iiorb=1,tmb%orbs%norb)
  !!     do iiorb=1,tmbder%orbs%norb
  !!         write(*,*) iorb, iiorb, dovrlp(iorb,iiorb)
  !!     end do
  !!  end do
  !!end if
  !END DEBUG

   ! needs generalizing if dovrlp and dham are to have different structures
   call to_zero(3*at%astruct%nat, fpulay(1,1))
   do jdir=1,3
     !do ialpha=1,tmb%orbs%norb
     do ispin=1,tmb%linmat%m%nspin
         ishift=(ispin-1)*tmb%linmat%m%nvctr
         if (tmb%linmat%m%nfvctrp>0) then
             isegstart=dham(jdir)%istsegline(tmb%linmat%m%isfvctr_par(iproc)+1)
             if (tmb%linmat%m%isfvctr+tmb%linmat%m%nfvctrp<tmb%linmat%m%nfvctr) then
                 isegend=dham(jdir)%istsegline(tmb%linmat%m%isfvctr_par(iproc+1)+1)-1
             else
                 isegend=dham(jdir)%nseg
             end if
             do iseg=isegstart,isegend
                  ii=dham(jdir)%keyv(iseg)-1
                  do jorb=dham(jdir)%keyg(1,1,iseg),dham(jdir)%keyg(2,1,iseg)
                      ii=ii+1
                      iialpha = dham(jdir)%keyg(1,2,iseg)
                      ibeta = jorb
                      jat=tmb%orbs%onwhichatom(iialpha)
                      kernel = 0.d0
                      ekernel= 0.d0
                      if (ispin==1) then
                          is=1
                          ie=orbs%norbu
                      else
                          is=orbs%norbu+1
                          ie=orbs%norb
                      end if
                      do iorb=is,ie
                          kernel  = kernel+orbs%occup(iorb)*tmb%coeff(iialpha,iorb)*tmb%coeff(ibeta,iorb)
                          !!ekernel = ekernel+tmb%orbs%eval(iorb)*orbs%occup(iorb) &
                          !!     *tmb%coeff(iialpha,iorb)*tmb%coeff(ibeta,iorb) 
                          ekernel = ekernel+orbs%eval(iorb)*orbs%occup(iorb) &
                               *tmb%coeff(iialpha,iorb)*tmb%coeff(ibeta,iorb) 
                      end do
                      fpulay(jdir,jat)=fpulay(jdir,jat)+&
                             2.0_gp*(kernel*dham_(jdir)%matrix_compr(ishift+ii)-ekernel*dovrlp_(jdir)%matrix_compr(ishift+ii))
                  end do
             end do
         end if
      end do
   end do 

   if (nproc > 1) then
      call mpiallred(fpulay(1,1), 3*at%astruct%nat, mpi_sum, bigdft_mpi%mpi_comm)
   end if

  if(iproc==0) then
       !!do jat=1,at%astruct%nat
       !!    write(*,'(a,i5,3es16.6)') 'iat, fpulay', jat, fpulay(1:3,jat)
       !!end do
       call yaml_comment('Pulay Correction',hfill='-')
       call yaml_sequence_open('Pulay Forces (Ha/Bohr)')
          do jat=1,at%astruct%nat
             call yaml_sequence(advance='no')
             call yaml_mapping_open(flow=.true.)
             call yaml_map(trim(at%astruct%atomnames(at%astruct%iatype(jat))),fpulay(1:3,jat),fmt='(1es20.12)')
             call yaml_mapping_close(advance='no')
             call yaml_comment(trim(yaml_toa(jat,fmt='(i4.4)')))
          end do
          call yaml_sequence_close()
  end if


  call f_free(psit_c)
  call f_free(psit_f)
  call f_free(hpsit_c)
  call f_free(hpsit_f)
  call f_free(lpsit_c)
  call f_free(lpsit_f)
  call f_free(lhphilarge)

  call f_free_ptr(denspot%pot_work)


  do jdir=1,3
     call deallocate_sparse_matrix(dovrlp(jdir))
     call deallocate_sparse_matrix(dham(jdir))
     call deallocate_matrices(dovrlp_(jdir))
     call deallocate_matrices(dham_(jdir))
  end do

  !!if(iproc==0) write(*,'(1x,a)') 'done.'

  call f_release_routine()

end subroutine pulay_correction
