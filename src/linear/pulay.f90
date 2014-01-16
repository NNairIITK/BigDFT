subroutine pulay_correction_new(iproc, nproc, tmb, orbs, at, fpulay)
  use module_base
  use module_types
  use module_interfaces, except_this_one => pulay_correction_new
  use yaml_output
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
  real(kind=8),dimension(:),allocatable :: hphit_c, hphit_f, denskern_tmp, delta_phit_c, delta_phit_f
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
  call mpiallred(energykernel(1,1), tmb%orbs%norb**2, mpi_sum, bigdft_mpi%mpi_comm, ierr)

  ! calculate the overlap matrix
  if(.not.associated(tmb%psit_c)) then
      isize=sum(tmb%collcom%nrecvcounts_c)
      tmb%psit_c=f_malloc_ptr(isize,id='tmb%psit_c')
  end if
  if(.not.associated(tmb%psit_f)) then
      isize=7*sum(tmb%collcom%nrecvcounts_f)
      tmb%psit_f=f_malloc_ptr(isize,id=' tmb%psit_f')
  end if
  call transpose_localized(iproc, nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
       tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd)
  call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%collcom, tmb%psit_c, &
       tmb%psit_c, tmb%psit_f, tmb%psit_f, tmb%linmat%ovrlp)

  call f_free_ptr(tmb%psit_c)
  call f_free_ptr(tmb%psit_f)


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
  tmb%linmat%denskern%matrix=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='tmb%linmat%denskern%matrix')
  call uncompressMatrix(iproc,tmb%linmat%denskern)
  call to_zero(3*at%astruct%nat, fpulay(1,1))
  do idir=1,3
      ! calculate the overlap matrix among hphi and phi_delta_large
      call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
           phi_delta_large(1,idir), delta_phit_c, delta_phit_f, tmb%ham_descr%lzd)
      call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%ham_descr%collcom, &
           hphit_c, delta_phit_c, hphit_f, delta_phit_f, tmb%linmat%ham)
      tmb%linmat%ham%matrix=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='tmb%linmat%ham%matrix')
      call uncompressMatrix(iproc,tmb%linmat%ham)

      do iorb=1,tmb%orbs%norbp
          iiorb=tmb%orbs%isorb+iorb
          iat=tmb%orbs%onwhichatom(iiorb)
          tt=0.d0
          do jorb=1,tmb%orbs%norb
              tt = tt -2.d0*tmb%linmat%denskern%matrix(jorb,iiorb)*tmb%linmat%ham%matrix(jorb,iiorb)
              !if (iproc==0) write(*,*) 'kern, ovrlp', tmb%linmat%denskern%matrix(jorb,iiorb), tmb%linmat%ham%matrix(iiorb,jorb)
          end do  
          fpulay(idir,iat)=fpulay(idir,iat)+tt
      end do
      call f_free_ptr(tmb%linmat%ham%matrix)
  end do
  call mpiallred(fpulay(1,1), 3*at%astruct%nat, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  call f_free_ptr(tmb%linmat%denskern%matrix)

  if(iproc==0) then
       call yaml_comment('new Pulay correction',hfill='-')
       call yaml_open_sequence('Pulay forces (Ha/Bohr)')
          do iat=1,at%astruct%nat
             call yaml_sequence(advance='no')
             call yaml_open_map(flow=.true.)
             call yaml_map(trim(at%astruct%atomnames(at%astruct%iatype(iat))),fpulay(1:3,iat),fmt='(1es20.12)')
             call yaml_close_map(advance='no')
             call yaml_comment(trim(yaml_toa(iat,fmt='(i4.4)')))
          end do
          call yaml_close_sequence()
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
      tmb%linmat%ovrlp%matrix=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='tmb%linmat%ovrlp%matrix')
      call uncompressMatrix(iproc,tmb%linmat%ovrlp)
      call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norb, tmb%orbs%norb, 1.d0, &
                 tmb%linmat%ovrlp%matrix, tmb%orbs%norb, energykernel, tmb%orbs%norb, &
                 0.d0, tempmat, tmb%orbs%norb)
      call f_free_ptr(tmb%linmat%ovrlp%matrix)
      isize=sum(tmb%ham_descr%collcom%nrecvcounts_c)
      hphit_c=f_malloc(isize,id='hphit_c')
      isize=7*sum(tmb%ham_descr%collcom%nrecvcounts_f)
      hphit_f=f_malloc(isize,id='hphit_f')
      call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
                               tmb%hpsi, hphit_c, hphit_f, tmb%ham_descr%lzd)
      isize=size(tmb%linmat%denskern%matrix_compr)
      denskern_tmp=f_malloc(isize,id='denskern_tmp')
      denskern_tmp=tmb%linmat%denskern%matrix_compr
      tmb%linmat%denskern%matrix=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='tmb%linmat%denskern%matrix')
      tmb%linmat%denskern%matrix=tempmat
      call compress_matrix_for_allreduce(iproc,tmb%linmat%denskern)
      call f_free_ptr(tmb%linmat%denskern%matrix)
      call build_linear_combination_transposed(tmb%ham_descr%collcom, &
           tmb%linmat%denskern, tmb%ham_descr%psit_c, tmb%ham_descr%psit_f, .false., hphit_c, hphit_f, iproc)
      tmb%linmat%denskern%matrix_compr=denskern_tmp
    
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





subroutine pulay_correction(iproc, nproc, orbs, at, rxyz, nlpspd, proj, SIC, denspot, GPU, tmb, fpulay)
  use module_base
  use module_types
  use module_interfaces, except_this_one => pulay_correction
  use yaml_output
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(orbitals_data),intent(in) :: orbs
  type(atoms_data),intent(in) :: at
  real(kind=8),dimension(3,at%astruct%nat),intent(in) :: rxyz
  type(nonlocal_psp_descriptors),intent(in) :: nlpspd
  real(wp),dimension(nlpspd%nprojel),intent(inout) :: proj
  type(SIC_data),intent(in) :: SIC
  type(DFT_local_fields), intent(inout) :: denspot
  type(GPU_pointers),intent(inout) :: GPU
  type(DFT_wavefunction),intent(inout) :: tmb
  real(kind=8),dimension(3,at%astruct%nat),intent(out) :: fpulay

  ! Local variables
  integer:: istat, iall, ierr, iialpha, jorb
  integer:: iorb, ii, iseg, isegstart, isegend
  integer:: jat, jdir, ibeta
  !!integer :: ialpha, iat, iiorb
  real(kind=8) :: kernel, ekernel
  real(kind=8),dimension(:),allocatable :: lhphilarge, psit_c, psit_f, hpsit_c, hpsit_f, lpsit_c, lpsit_f
  type(sparseMatrix) :: dovrlp(3), dham(3)
  type(energy_terms) :: energs
  type(confpot_data),dimension(:),allocatable :: confdatarrtmp
  character(len=*),parameter :: subname='pulay_correction'

  ! Begin by updating the Hpsi
  call local_potential_dimensions(tmb%ham_descr%lzd,tmb%orbs,denspot%dpbox%ngatherarr(0,1))

  allocate(lhphilarge(tmb%ham_descr%npsidim_orbs), stat=istat)
  call memocc(istat, lhphilarge, 'lhphilarge', subname)
  call to_zero(tmb%ham_descr%npsidim_orbs,lhphilarge(1))

  !!call post_p2p_communication(iproc, nproc, denspot%dpbox%ndimpot, denspot%rhov, &
  !!     tmb%ham_descr%comgp%nrecvbuf, tmb%ham_descr%comgp%recvbuf, tmb%ham_descr%comgp, tmb%ham_descr%lzd)
  call start_onesided_communication(iproc, nproc, denspot%dpbox%ndimpot, denspot%rhov, &
       tmb%ham_descr%comgp%nrecvbuf, tmb%ham_descr%comgp%recvbuf, tmb%ham_descr%comgp, tmb%ham_descr%lzd)

  allocate(confdatarrtmp(tmb%orbs%norbp))
  call default_confinement_data(confdatarrtmp,tmb%orbs%norbp)


  call NonLocalHamiltonianApplication(iproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,rxyz,&
       proj,tmb%ham_descr%lzd,nlpspd,tmb%ham_descr%psi,lhphilarge,energs%eproj)

  ! only kinetic because waiting for communications
  call LocalHamiltonianApplication(iproc,nproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,&
       tmb%ham_descr%lzd,confdatarrtmp,denspot%dpbox%ngatherarr,denspot%pot_work,tmb%ham_descr%psi,lhphilarge,&
       energs,SIC,GPU,3,pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,potential=denspot%rhov,comgp=tmb%ham_descr%comgp)
  call full_local_potential(iproc,nproc,tmb%orbs,tmb%ham_descr%lzd,2,denspot%dpbox,denspot%rhov,denspot%pot_work, &
       tmb%ham_descr%comgp)
  ! only potential
  call LocalHamiltonianApplication(iproc,nproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,&
       tmb%ham_descr%lzd,confdatarrtmp,denspot%dpbox%ngatherarr,denspot%pot_work,tmb%ham_descr%psi,lhphilarge,&
       energs,SIC,GPU,2,pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,potential=denspot%rhov,comgp=tmb%ham_descr%comgp)

  call timing(iproc,'glsynchham1','ON') !lr408t
  call SynchronizeHamiltonianApplication(nproc,tmb%ham_descr%npsidim_orbs,tmb%orbs,tmb%ham_descr%lzd,GPU,lhphilarge,&
       energs%ekin,energs%epot,energs%eproj,energs%evsic,energs%eexctX)
  call timing(iproc,'glsynchham1','OF') !lr408t
  deallocate(confdatarrtmp)
  

  ! Now transpose the psi and hpsi
  allocate(lpsit_c(tmb%ham_descr%collcom%ndimind_c))
  call memocc(istat, lpsit_c, 'lpsit_c', subname)
  allocate(lpsit_f(7*tmb%ham_descr%collcom%ndimind_f))
  call memocc(istat, lpsit_f, 'lpsit_f', subname)
  allocate(hpsit_c(tmb%ham_descr%collcom%ndimind_c))
  call memocc(istat, hpsit_c, 'hpsit_c', subname)
  allocate(hpsit_f(7*tmb%ham_descr%collcom%ndimind_f))
  call memocc(istat, hpsit_f, 'hpsit_f', subname)
  allocate(psit_c(tmb%ham_descr%collcom%ndimind_c))
  call memocc(istat, psit_c, 'psit_c', subname)
  allocate(psit_f(7*tmb%ham_descr%collcom%ndimind_f))
  call memocc(istat, psit_f, 'psit_f', subname)

  call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
       tmb%ham_descr%psi, lpsit_c, lpsit_f, tmb%ham_descr%lzd)

  call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
       lhphilarge, hpsit_c, hpsit_f, tmb%ham_descr%lzd)

  !now build the derivative and related matrices <dPhi_a | H | Phi_b> and <dPhi_a | Phi_b>

  ! DOVRLP AND DHAM SHOULD HAVE DIFFERENT SPARSITIES, BUT TO MAKE LIFE EASIER KEEPING THEM THE SAME FOR NOW
  ! also array of structure a bit inelegant at the moment
  do jdir = 1, 3
    call nullify_sparsematrix(dovrlp(jdir))
    call nullify_sparsematrix(dham(jdir))
    call sparse_copy_pattern(tmb%linmat%ham,dovrlp(jdir),iproc,subname) 
    call sparse_copy_pattern(tmb%linmat%ham,dham(jdir),iproc,subname)
    allocate(dham(jdir)%matrix_compr(dham(jdir)%nvctr), stat=istat)
    call memocc(istat, dham(jdir)%matrix_compr, 'dham%matrix_compr', subname)
    allocate(dovrlp(jdir)%matrix_compr(dovrlp(jdir)%nvctr), stat=istat)
    call memocc(istat, dovrlp(jdir)%matrix_compr, 'dovrlp%matrix_compr', subname)

    call get_derivative(jdir, tmb%ham_descr%npsidim_orbs, tmb%ham_descr%lzd%hgrids(1), tmb%orbs, &
         tmb%ham_descr%lzd, tmb%ham_descr%psi, lhphilarge)

    call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
         lhphilarge, psit_c, psit_f, tmb%ham_descr%lzd)

    call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%ham_descr%collcom,&
         psit_c, lpsit_c, psit_f, lpsit_f, dovrlp(jdir))

    call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%ham_descr%collcom,&
         psit_c, hpsit_c, psit_f, hpsit_f, dham(jdir))
  end do


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
     if (tmb%orbs%norbp>0) then
         isegstart=dham(jdir)%istsegline(tmb%orbs%isorb_par(iproc)+1)
         if (tmb%orbs%isorb+tmb%orbs%norbp<tmb%orbs%norb) then
             isegend=dham(jdir)%istsegline(tmb%orbs%isorb_par(iproc+1)+1)-1
         else
             isegend=dham(jdir)%nseg
         end if
         do iseg=isegstart,isegend
              ii=dham(jdir)%keyv(iseg)-1
              do jorb=dham(jdir)%keyg(1,iseg),dham(jdir)%keyg(2,iseg)
                  ii=ii+1
                  iialpha = (jorb-1)/tmb%orbs%norb + 1
                  ibeta = jorb - (iialpha-1)*tmb%orbs%norb
                  jat=tmb%orbs%onwhichatom(iialpha)
                  kernel = 0.d0
                  ekernel= 0.d0
                  do iorb=1,orbs%norb
                      kernel  = kernel+orbs%occup(iorb)*tmb%coeff(iialpha,iorb)*tmb%coeff(ibeta,iorb)
                      ekernel = ekernel+tmb%orbs%eval(iorb)*orbs%occup(iorb) &
                           *tmb%coeff(iialpha,iorb)*tmb%coeff(ibeta,iorb) 
                  end do
                  fpulay(jdir,jat)=fpulay(jdir,jat)+&
                         2.0_gp*(kernel*dham(jdir)%matrix_compr(ii)-ekernel*dovrlp(jdir)%matrix_compr(ii))
              end do
         end do
     end if
   end do 

   call mpiallred(fpulay(1,1), 3*at%astruct%nat, mpi_sum, bigdft_mpi%mpi_comm, ierr)

  if(iproc==0) then
       !!do jat=1,at%astruct%nat
       !!    write(*,'(a,i5,3es16.6)') 'iat, fpulay', jat, fpulay(1:3,jat)
       !!end do
       call yaml_comment('Pulay Correction',hfill='-')
       call yaml_open_sequence('Pulay Forces (Ha/Bohr)')
          do jat=1,at%astruct%nat
             call yaml_sequence(advance='no')
             call yaml_open_map(flow=.true.)
             call yaml_map(trim(at%astruct%atomnames(at%astruct%iatype(jat))),fpulay(1:3,jat),fmt='(1es20.12)')
             call yaml_close_map(advance='no')
             call yaml_comment(trim(yaml_toa(jat,fmt='(i4.4)')))
          end do
          call yaml_close_sequence()
  end if

  iall=-product(shape(psit_c))*kind(psit_c)
  deallocate(psit_c, stat=istat)
  call memocc(istat, iall, 'psit_c', subname)
  iall=-product(shape(psit_f))*kind(psit_f)
  deallocate(psit_f, stat=istat)
  call memocc(istat, iall, 'psit_f', subname)
  iall=-product(shape(hpsit_c))*kind(hpsit_c)
  deallocate(hpsit_c, stat=istat)
  call memocc(istat, iall, 'hpsit_c', subname)
  iall=-product(shape(hpsit_f))*kind(hpsit_f)
  deallocate(hpsit_f, stat=istat)
  call memocc(istat, iall, 'hpsit_f', subname)
  iall=-product(shape(lpsit_c))*kind(lpsit_c)
  deallocate(lpsit_c, stat=istat)
  call memocc(istat, iall, 'lpsit_c', subname)
  iall=-product(shape(lpsit_f))*kind(lpsit_f)
  deallocate(lpsit_f, stat=istat)
  call memocc(istat, iall, 'lpsit_f', subname)

  iall=-product(shape(lhphilarge))*kind(lhphilarge)
  deallocate(lhphilarge, stat=istat)
  call memocc(istat, iall, 'lhphilarge', subname)

  iall=-product(shape(denspot%pot_work))*kind(denspot%pot_work)
  deallocate(denspot%pot_work, stat=istat)
  call memocc(istat, iall, 'denspot%pot_work', subname)

  do jdir=1,3
     call deallocate_sparseMatrix(dovrlp(jdir),subname)
     call deallocate_sparseMatrix(dham(jdir),subname)
  end do

  !!if(iproc==0) write(*,'(1x,a)') 'done.'

end subroutine pulay_correction

