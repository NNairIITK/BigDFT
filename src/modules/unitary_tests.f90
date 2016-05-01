!> @file
!!  Unitary tests for the communication subroutines (potential and density)
!! @author
!!    Copyright (C) 2015-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Module containing unitary tests for the potential and density communications
module unitary_tests
  use module_base
  implicit none

  private

  !> Public routines
  public :: check_communication_potential
  public :: check_communication_sumrho
  public :: check_communications_locreg


  contains


    !> Perform the communication needed for the potential and verify that the results is as expected
    subroutine check_communication_potential(iproc,nproc,denspot,tmb)
      use module_base
      use module_types
      use yaml_output
      use dictionaries, only: f_err_throw
      use communications, only: start_onesided_communication
      use rhopotential, only: full_local_potential
      implicit none
      integer,intent(in) :: iproc,nproc
      type(DFT_wavefunction), intent(inout) :: tmb
      type(DFT_local_fields), intent(inout) :: denspot
      !local variables
      logical :: dosome, abort, wrong
      integer :: i1,i2,i3,ind,n3p,ilr,iorb,ilr_orb,n2i,n1i,numtot,ishift,ispin
      integer :: i1s, i1e, i2s, i2e, i3s, i3e, ii1, ii2, ii3, ierr, jproc
      !integer :: ierr
      real(dp) :: maxdiff,sumdiff,testval
      real(dp),parameter :: tol_calculation_mean=1.d-12
      real(dp),parameter :: tol_calculation_max=1.d-10
      character(len=200), parameter :: subname='check_communication_potential'
      integer,dimension(:),allocatable :: n3p_withmax
    
      call timing(bigdft_mpi%iproc,'check_pot','ON')
    
      !assign constants
      i3s=denspot%dpbox%nscatterarr(bigdft_mpi%iproc,3)+1 !< starting point of the planes in the z direction
      n3p=denspot%dpbox%nscatterarr(bigdft_mpi%iproc,2) !< number of planes for the potential
      n2i=denspot%dpbox%mesh%ndims(2) !< size of the global domain in y direction
      n1i=denspot%dpbox%mesh%ndims(1) !< size of the global domain in x direction
    
      !fill the values of the rhov array
      ind=0
      do ispin=1,denspot%dpbox%nrhodim
          ishift=(ispin-1)*n1i*n2i*n3p
          do i3=i3s,i3s+n3p-1
             do i2=1,n2i
                do i1=1,n1i
                   ind=ind+1
                   !denspot%rhov(ind)=real(ishift+i1+(i2-1)*n1i+(i3-1)*n1i*n2i,dp)
                   denspot%rhov(ind)=real((-1)**(ispin+1)*(i1+(i2-1)*n1i+(i3-1)*n1i*n2i),dp)
                   !write(500,'(es16.8)') denspot%rhov(ind)
                end do
             end do
          end do
      end do

     
      !!write(*,'(a,3i12)') 'iproc, denspot%dpbox%ndimpot*denspot%dpbox%nrhodim, size(denspot%rhov)', iproc, denspot%dpbox%ndimpot*denspot%dpbox%nrhodim, size(denspot%rhov)
    
      !calculate the dimensions and communication of the potential element with mpi_get
      !max(denspot%dpbox%nscatterarr(:,2),1)
      !denspot%dpbox%nscatterarr(:,2) creates a temporary array, to be avoided
      call local_potential_dimensions(iproc,tmb%ham_descr%lzd,tmb%orbs,denspot%xc,denspot%dpbox%ngatherarr(0,1))
      n3p_withmax = f_malloc(0.to.nproc-1,id='n3p_withmax')
      do jproc=0,nproc-1
          n3p_withmax(jproc) = max(denspot%dpbox%nscatterarr(jproc,2),1)
      end do
      call start_onesided_communication(bigdft_mpi%iproc, bigdft_mpi%nproc, &
           denspot%dpbox%mesh%ndims(1), denspot%dpbox%mesh%ndims(2), n3p_withmax, denspot%rhov, &
           tmb%ham_descr%comgp%nspin*tmb%ham_descr%comgp%nrecvbuf, tmb%ham_descr%comgp%recvbuf, &
           tmb%ham_descr%comgp, tmb%ham_descr%lzd)
      call f_free(n3p_withmax)
    
      !check the fetching of the potential element, destroy the MPI window, results in pot_work
      !!write(*,*) 'kind(2)',kind(2)
      call full_local_potential(bigdft_mpi%iproc,bigdft_mpi%nproc,tmb%orbs,tmb%ham_descr%lzd,&
           2,denspot%dpbox,denspot%xc,denspot%rhov,denspot%pot_work,tmb%ham_descr%comgp)
    
      !!do ind=1,tmb%ham_descr%comgp%nspin*tmb%ham_descr%comgp%nrecvbuf
      !!    write(5200+iproc,'(a,i10,es16.7)') 'ind, val', ind, tmb%ham_descr%comgp%recvbuf(ind)
      !!end do
    
    
      maxdiff=0.0_dp
      sumdiff=0.0_dp
      numtot=0
      loop_lr: do ilr=1,tmb%ham_descr%Lzd%nlr
         !check if this localisation region is used by one of the orbitals
         dosome=.false.
         do iorb=1,tmb%orbs%norbp
            dosome = (tmb%orbs%inwhichlocreg(iorb+tmb%orbs%isorb) == ilr)
            if (dosome) exit
         end do
         if (.not. dosome) cycle loop_lr
    
         loop_orbs: do iorb=1,tmb%orbs%norbp
            if (tmb%orbs%spinsgn(iorb+tmb%orbs%isorb)>0.d0) then
                ispin=1
            else
                ispin=2
            end if
            ilr_orb=tmb%orbs%inwhichlocreg(iorb+tmb%orbs%isorb)
            if (ilr_orb /= ilr) cycle loop_orbs
    
            ind=tmb%orbs%ispot(iorb)-1
            i3s=tmb%ham_descr%Lzd%Llr(ilr)%nsi3+1
            i3e=i3s+tmb%ham_descr%Lzd%Llr(ilr)%d%n3i-1
            i2s=tmb%ham_descr%Lzd%Llr(ilr)%nsi2+1
            i2e=i2s+tmb%ham_descr%Lzd%Llr(ilr)%d%n2i-1
            i1s=tmb%ham_descr%Lzd%Llr(ilr)%nsi1+1
            i1e=i1s+tmb%ham_descr%Lzd%Llr(ilr)%d%n1i-1
            !do i3=1,tmb%ham_descr%Lzd%Llr(ilr)%d%n3i
            do i3=i3s,i3e
               ii3=modulo(i3-1,tmb%ham_descr%Lzd%glr%d%n3i)+1
               !do i2=1,tmb%ham_descr%Lzd%Llr(ilr)%d%n2i
               do i2=i2s,i2e
                  ii2=modulo(i2-1,tmb%ham_descr%Lzd%glr%d%n2i)+1
                  !do i1=1,tmb%ham_descr%Lzd%Llr(ilr)%d%n1i
                  do i1=i1s,i1e
                     ii1=modulo(i1-1,tmb%ham_descr%Lzd%glr%d%n1i)+1
                     ind=ind+1
                     !!testval=real(i1+tmb%ham_descr%Lzd%Llr(ilr)%nsi1+&
                     !!     (i2+tmb%ham_descr%Lzd%Llr(ilr)%nsi2-1)*n1i+&
                     !!     (i3+tmb%ham_descr%Lzd%Llr(ilr)%nsi3-1)*n1i*n2i,dp)
                     !testval=real((-1)**(ispin+1)*(i1+tmb%ham_descr%Lzd%Llr(ilr)%nsi1+&
                     !     (i2+tmb%ham_descr%Lzd%Llr(ilr)%nsi2-1)*n1i+&
                     !     (i3+tmb%ham_descr%Lzd%Llr(ilr)%nsi3-1)*n1i*n2i),dp)
                     testval=real((-1)**(ispin+1)*(ii1+(ii2-1)*n1i+(ii3-1)*n1i*n2i),dp)
                     !if (iproc==0) write(*,'(a,5i8,2es14.3)') 'ispin, i1, i2, i3, ind, val, ref', ispin, i1, i2, i3, ind, denspot%pot_work(ind), testval
                     !write(2000+iproc,'(a,6i8,2es14.3)') 'ispin, ilr, i1, i2, i3, ind, val, ref', &
                     !    ispin, ilr, i1, i2, i3, ind, denspot%pot_work(ind), testval
                     testval=abs(denspot%pot_work(ind)-testval)
                     maxdiff=max(maxdiff,testval)
                     sumdiff=sumdiff+testval
                     numtot=numtot+1
                  end do
               end do
            end do
    
         end do loop_orbs
    
      end do loop_lr
    
      if (numtot>0) sumdiff = sumdiff/numtot
    
      ! Reduce the results
      if (bigdft_mpi%nproc>1) then
          call mpiallred(sumdiff, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
          call mpiallred(maxdiff, 1, mpi_max, comm=bigdft_mpi%mpi_comm)
      end if
        
      ! Get mean value for the sum
      sumdiff=sqrt(sumdiff)
    
      if (bigdft_mpi%iproc==0) call yaml_mapping_open('Checking operations for potential communication')    
      ! Print the results
      if (bigdft_mpi%iproc==0) then
          call yaml_map('Tolerance for the following test',tol_calculation_mean,fmt='(1es25.18)')
          if (sumdiff>tol_calculation_mean) then
             call yaml_warning('CALCULATION ERROR: total difference of '//trim(yaml_toa(sumdiff,fmt='(1es25.18)')))
          else
             call yaml_map('calculation check, error sum', sumdiff,fmt='(1es25.18)')
          end if
          call yaml_map('Tolerance for the following test',tol_calculation_max,fmt='(1es25.18)')
          if (maxdiff>tol_calculation_max) then
             call yaml_warning('CALCULATION ERROR: max difference of '//trim(yaml_toa(maxdiff,fmt='(1es25.18)')))
          else
             call yaml_map('calculation check, error max', maxdiff,fmt='(1es25.18)')
          end if
      end if
      if (bigdft_mpi%iproc==0) call yaml_mapping_close()
    
      wrong = (sumdiff>tol_calculation_mean .or. maxdiff>tol_calculation_max)
      abort = .false.
      if (wrong) then
          call yaml_warning('The communication of the potential is not correct for this setup, check communication routines')
          !abort = .true.
          !!call f_err_throw('The communication of the potential is not correct for this setup, check communication routines',&
          !!        err_name='BIGDFT_MPI_ERROR')
      end if
      if (abort) call MPI_ABORT(bigdft_mpi%mpi_comm,10,ierr)
    
      call f_free_ptr(denspot%pot_work)
    
      nullify(denspot%pot_work)
    
      call timing(bigdft_mpi%iproc,'check_pot','OF')
    
      !!call mpi_finalize(ierr)
      !!stop
    
    end subroutine check_communication_potential
    
    
    subroutine check_communication_sumrho(iproc, nproc, orbs, lzd, collcom_sr, denspot, denskern, denskern_, check_sumrho)
      use module_base
      use module_types
      use yaml_output
      use locregs, only: check_whether_bounds_overlap
      use communications, only: transpose_switch_psir, transpose_communicate_psir, transpose_unswitch_psirt
      use sparsematrix_base, only: sparse_matrix, matrices
      use sparsematrix_init, only: matrixindex_in_compressed
      use rhopotential, only: sumrho_for_TMBs
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(local_zone_descriptors),intent(in) :: lzd
      type(orbitals_data),intent(in) :: orbs
      type(comms_linear),intent(inout) :: collcom_sr
      type(DFT_local_fields),intent(in) :: denspot
      type(sparse_matrix),intent(in) :: denskern
      type(matrices),intent(inout) :: denskern_
      integer,intent(in) :: check_sumrho
    
      ! Local variables
      integer :: ist, iorb, iiorb, ilr, i, ii, iixyz, nxyz, ipt, i0, ierr, jproc
      integer :: i1, i2, i3, is1, is2, is3, ie1, ie2, ie3, ii3s, ii3e, nmax, jj, j, ind, ikernel
      !integer :: iz, iy, ix, iix, iiy, iiz, iim
      integer :: iorbmin, iorbmax, jorb, ispin,ishift
      integer :: n2, n3, ii1, ii2, ii3
      real(kind=8) :: maxdiff, sumdiff, tt, tti, ttj, hxh, hyh, hzh, factor, ref_value
      real(kind=8) :: diff
      real(kind=8),dimension(:),allocatable :: psir, psirwork, psirtwork, rho, rho_check
      integer,dimension(:,:,:),allocatable :: weight
      integer,dimension(:,:,:,:),allocatable :: orbital_id
      integer,dimension(:),allocatable :: istarr
      integer,dimension(:,:),allocatable :: matrixindex_in_compressed_auxilliary
      real(kind=8),parameter :: tol_transpose=1.d-14
      real(kind=8),parameter :: tol_calculation_mean=1.d-12
      real(kind=8),parameter :: tol_calculation_max=1.d-10
      character(len=*), parameter :: subname='check_sumrho'
      logical :: rho_negative
    
      call timing(iproc,'check_sumrho','ON')
    
      if (iproc==0) call yaml_mapping_open('Checking operations for sumrho')
    
      call f_routine(id='check_communication_sumrho')
    
      ! Allocate all the main arrays arrays
      psir=f_malloc(collcom_sr%ndimpsi_c,id='psir') !direct array
    
      psirwork=f_malloc(collcom_sr%ndimpsi_c,id='psirwork') !direct workarray
    
      psirtwork=f_malloc(collcom_sr%ndimind_c,id='psirtwork') !transposed workarray
    
      ! Size of global box
      nxyz=lzd%glr%d%n1i*lzd%glr%d%n2i*lzd%glr%d%n3i
    
      call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
      call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
    
      ! Fill the direct array with a recognizable pattern
      ist=0
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          ilr=orbs%inWhichLocreg(iiorb)
          !!$omp parallel default(none) &
          !!$omp shared(orbs, lzd, psir, iorb, iiorb, ilr, ist, nxyz) &
          !!$omp private(i, ii, iz, iy, ix, iix, iiy, iiz, iixyz)
          !!$omp do
          !!do i=1,lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i
          !!    ! coordinates within locreg
          !!    ii=i-1
          !!    iz=ii/(lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i)+1
          !!    ii=ii-(iz-1)*lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i
          !!    iy=ii/lzd%llr(ilr)%d%n1i+1
          !!    ix=ii-(iy-1)*lzd%llr(ilr)%d%n1i+1
          !!    ! coordinates within global region
          !!    iix=ix+lzd%llr(ilr)%nsi1
          !!    iiy=iy+lzd%llr(ilr)%nsi2
          !!    iiz=iz+lzd%llr(ilr)%nsi3
          !!    iixyz=(iiz-1)*lzd%glr%d%n1i*lzd%glr%d%n2i+(iiy-1)*lzd%glr%d%n1i+iix
          !!    ! assign unique value
          !!    !psir(ist+i)=dble(iiorb)*orbs%spinsgn(iiorb)!test_value_sumrho(iiorb,iixyz,nxyz)
          !!    psir(ist+i)=test_value_sumrho(iiorb,iixyz,nxyz)
          !!end do
          !!$omp end do
          !!$omp end parallel
          is3 = modulo(1+lzd%llr(ilr)%nsi3-1,lzd%glr%d%n3i)+1
          ie3 = is3+lzd%llr(ilr)%d%n3i-1
          is2 = modulo(1+lzd%llr(ilr)%nsi2-1,lzd%glr%d%n2i)+1
          ie2 = is2+lzd%llr(ilr)%d%n2i-1
          is1 = modulo(1+lzd%llr(ilr)%nsi1-1,lzd%glr%d%n1i)+1
          ie1 = is1+lzd%llr(ilr)%d%n1i-1
          i = 0
          do i3=is3,ie3
              ii3 = modulo(i3-1,lzd%glr%d%n3i)+1
              n3 = (ii3-1)*lzd%glr%d%n1i*lzd%glr%d%n2i
              do i2=is2,ie2
                  ii2 = modulo(i2-1,lzd%glr%d%n2i)+1
                  n2 = (ii2-1)*lzd%glr%d%n1i
                  do i1=is1,ie1
                      ii1 = modulo(i1-1,lzd%glr%d%n1i)+1
                      iixyz = n3 + n2 + ii1
                      i = i + 1
                      psir(ist+i)=test_value_sumrho(iiorb,iixyz,nxyz)
                  end do
              end do
          end do
          ist = ist + lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i
      end do
      if(ist/=collcom_sr%ndimpsi_c) then
          write(*,'(a,i0,a)') 'ERROR on process ',iproc,' : ist/=collcom_sr%ndimpsi_c'
          stop
      end if
    
      call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
      call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
    
      ! Rearrange data
      call transpose_switch_psir(collcom_sr, psir, psirwork)
    
      !!! TEST #################################
      !!do ipt=1,collcom_sr%ndimpsi_c
      !!if (iproc==0) write(*,'(a,i9,2f13.2)') 'psir, psirwork: ipt, vals', ipt, psir(ipt), psirwork(ipt)
      !!end do
      !!! END TEST #############################
    
      ! direct array not needed anymore
      call f_free(psir)
    
      ! Communicate the data
      call transpose_communicate_psir(iproc, nproc, collcom_sr, psirwork, psirtwork)
      !!do ipt=1,size(psirwork)
      !!    write(800+iproc,*) psirwork(ipt)
      !!end do
      !!do ipt=1,size(psirtwork)
      !!    write(810+iproc,*) psirtwork(ipt)
      !!end do
      !! TEST #################################
      !do ipt=1,collcom_sr%ndimind_c
      !    if (iproc==0) write(*,'(a,i9,2f13.2)') 'psirtwork, ipt, val', ipt, psirtwork(ipt)
      !end do
      !! END TEST #############################
    
      ! Direct workarray not needed anymore
      call f_free(psirwork)
    
      ! Rearrange array
      !LG: WARNING: it is bad practice to consider collcom_sr as intent(in)
      !and collcom_sr%psit_c and intent(out) or intent(inout)!!!
      call transpose_unswitch_psirt(collcom_sr, psirtwork, collcom_sr%psit_c)
    
      !! TEST #################################
      !do ipt=1,collcom_sr%ndimind_c
      !    if (iproc==0) write(*,'(a,i9,2f13.2)') 'psirtwork, collcom_sr%psit_c, ipt, vals', ipt, psirtwork(ipt), collcom_sr%psit_c(ipt)
      !end do
      !! END TEST #############################
    
    
      !!! TEST #################################
      !!do ispin=1,denskern%nspin
      !!    do ipt=1,collcom_sr%nptsp_c
      !!        ii=collcom_sr%norb_per_gridpoint_c(ipt)
      !!        i0=collcom_sr%isptsp_c(ipt)+(ispin-1)*collcom_sr%ndimind_c/2
      !!        do i=1,ii
      !!            if (iproc==0) write(3000+iproc,'(a,4i9,f11.2)') 'ipt, i0, i, npg, val', ipt, i0, i, ii, collcom_sr%psit_c(i0+i)
      !!        end do
      !!    end do
      !!end do
      !!! END TEST #############################
    
      ! Transposed workarray not needed anymore
      call f_free(psirtwork)
    
      call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
      call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
    
      ! Check the layout of the transposed data
      maxdiff=0.d0
      sumdiff=0.d0
    
      ! Get the starting point of each MPI task
      istarr=f_malloc((/0.to.nproc-1/),id='istarr')
      istarr=0
      istarr(iproc)=collcom_sr%nptsp_c
    
      if (nproc > 1) then
         call mpiallred(istarr(0), nproc, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if
    
      ist=0
      do jproc=0,iproc-1
          ist=ist+istarr(jproc)
      end do
      call f_free(istarr)
      
      ! Iterate through all the transposed values and check whether they are correct
      do ispin=1,denskern%nspin
          !$omp parallel default(none) &
          !$omp shared(collcom_sr, ist, nxyz, maxdiff, sumdiff, ispin) &
          !$omp private(ipt, ii, i0, iixyz, i, iiorb, tt, ref_value, diff)
          !$omp do reduction(+:sumdiff) reduction(max:maxdiff)
          do ipt=1,collcom_sr%nptsp_c
              ii=collcom_sr%norb_per_gridpoint_c(ipt)
              !i0=collcom_sr%isptsp_c(ipt)
              i0=collcom_sr%isptsp_c(ipt)+(ispin-1)*collcom_sr%ndimind_c/2
              iixyz=ist+ipt
              do i=1,ii
                  iiorb=collcom_sr%indexrecvorbital_c(i0+i)
                  tt=collcom_sr%psit_c(i0+i)
                  ref_value=test_value_sumrho(iiorb,iixyz,nxyz)
                  !write(*,'(a,4i9,2f11.2)') 'ipt, i0, i, npg, val, ref', ipt, i0, i, ii, tt, ref_value
                  diff=abs(tt-ref_value)
                  if (diff>maxdiff) maxdiff=diff
                  sumdiff=sumdiff+diff**2
              end do
          end do
          !$omp end do
          !$omp end parallel
      end do
    
    
      ! Reduce the results
      if (nproc>1) then
          call mpiallred(sumdiff, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
          call mpiallred(maxdiff, 1, mpi_max, comm=bigdft_mpi%mpi_comm)
      end if
      call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
      call mpi_barrier(bigdft_mpi%mpi_comm, ierr)
    
      ! Get mean value for the sum
      sumdiff = sumdiff/(lzd%glr%d%n1i*lzd%glr%d%n2i*lzd%glr%d%n3i)
      sumdiff=sqrt(sumdiff)
    
      ! Print the results
      if (iproc==0) then
          call yaml_map('Tolerance for the following test',tol_transpose,fmt='(1es25.18)')
          if (sumdiff>tol_transpose) then
             !call yaml_warning('TRANSPOSITION ERROR: mean difference of '//trim(yaml_toa(sumdiff,fmt='(1es25.18)')))
             call f_err_throw('TRANSPOSITION ERROR: mean difference of '//trim(yaml_toa(sumdiff,fmt='(1es25.18)')),&
                  err_name='BIGDFT_MPI_ERROR')
          else
             call yaml_map('transposition check, mean error ', sumdiff,fmt='(1es25.18)')
          end if
          if (maxdiff>tol_transpose) then
             !call yaml_warning('TRANSPOSITION ERROR: max difference of '//trim(yaml_toa(maxdiff,fmt='(1es25.18)')))
             call f_err_throw('TRANSPOSITION ERROR: max difference of '//trim(yaml_toa(maxdiff,fmt='(1es25.18)')), &
                  err_name='BIGDFT_MPI_ERROR')
          else
             call yaml_map('transposition check, max error ', maxdiff,fmt='(1es25.18)')
          end if
      end if
    
    
      ! Now comes the full check.. Do it depending on the value of check_sumrho
      if (check_sumrho==2) then
      
          ! Now simulate the calculation of the charge density. Take the same reproducable
          ! values as above. In this way the charge density can be calculated without
          ! the communication.
          
          ! First determine how many orbitals one has for each grid point in the current slice
          ii3s=denspot%dpbox%nscatterarr(iproc,3)-denspot%dpbox%nscatterarr(iproc,4)+1
          ii3e=denspot%dpbox%nscatterarr(iproc,3)-denspot%dpbox%nscatterarr(iproc,4)+denspot%dpbox%nscatterarr(iproc,1)
          weight=f_malloc0((/1.to.lzd%glr%d%n1i,1.to.lzd%glr%d%n2i,ii3s.to.ii3e/),id='weight')
    
          !if (denspot%dpbox%nscatterarr(iproc,1)>0) then
          !    call f_zero(lzd%glr%d%n1i*lzd%glr%d%n2i*denspot%dpbox%nscatterarr(iproc,1), weight(1,1,ii3s))
          !end if
    
          do i3=ii3s,ii3e
              ii3=modulo(i3-1,lzd%glr%d%n3i)+1
              do iorb=1,orbs%norbu
                  ilr=orbs%inwhichlocreg(iorb)
                  !is3=1+lzd%Llr(ilr)%nsi3
                  !ie3=lzd%Llr(ilr)%nsi3+lzd%llr(ilr)%d%n3i
                  is3=modulo(1+lzd%Llr(ilr)%nsi3-1,lzd%glr%d%n3i)+1
                  !ie3=is3+lzd%llr(ilr)%d%n3i-1
                  ie3=modulo(lzd%llr(ilr)%nsi3+lzd%llr(ilr)%d%n3i-1,lzd%glr%d%n3i)+1
                  !if (is3>i3 .or. i3>ie3) cycle
                  !if (.not.check_whether_bounds_overlap(is3,ie3,i3,i3)) cycle
                  if (.not.check_whether_bounds_overlap(is3,ie3,ii3,ii3)) cycle
                  !is1=1+lzd%Llr(ilr)%nsi1
                  !ie1=lzd%Llr(ilr)%nsi1+lzd%llr(ilr)%d%n1i
                  is1=modulo(1+lzd%Llr(ilr)%nsi1-1,lzd%glr%d%n1i)+1
                  ie1=is1+lzd%llr(ilr)%d%n1i-1
                  !is2=1+lzd%Llr(ilr)%nsi2
                  !ie2=lzd%Llr(ilr)%nsi2+lzd%llr(ilr)%d%n2i
                  is2=modulo(1+lzd%Llr(ilr)%nsi2-1,lzd%glr%d%n2i)+1
                  ie2=is2+lzd%llr(ilr)%d%n2i-1
                  !$omp parallel default(none) &
                  !$omp shared(is2, ie2, is1, ie1, weight, i3, lzd) &
                  !$omp private(i2, i1, ii2, ii1) 
                  !$omp do
                  do i2=is2,ie2
                      ii2=modulo(i2-1,lzd%glr%d%n2i)+1
                      do i1=is1,ie1
                          ii1=modulo(i1-1,lzd%glr%d%n1i)+1
                          weight(ii1,ii2,i3) = weight(ii1,ii2,i3)+1
                      end do
                  end do
                  !$omp end do
                  !$omp end parallel
              end do
          end do
        
          ! The array orbital_id contains the IDs of the orbitals touching a given gridpoint
          nmax=maxval(weight)
    
          orbital_id=f_malloc((/nmax,lzd%glr%d%n1i,lzd%glr%d%n2i,ii3e-ii3s+1/),lbounds=(/1,1,1,ii3s/),id='orbital_id')
    
          !if (denspot%dpbox%nscatterarr(iproc,1)>0) then
          !    call f_zero(lzd%glr%d%n1i*lzd%glr%d%n2i*denspot%dpbox%nscatterarr(iproc,1), weight(1,1,ii3s))
          !end if
          call f_zero(weight)
          iorbmin=1000000000
          iorbmax=-1000000000
          do i3=ii3s,ii3e
              ii3=modulo(i3-1,lzd%glr%d%n3i)+1
              do iorb=1,orbs%norbu
                  ilr=orbs%inwhichlocreg(iorb)
                  !is3=1+lzd%Llr(ilr)%nsi3
                  !ie3=lzd%Llr(ilr)%nsi3+lzd%llr(ilr)%d%n3i
                  is3=modulo(1+lzd%Llr(ilr)%nsi3-1,lzd%glr%d%n3i)+1
                  ie3=modulo(lzd%llr(ilr)%nsi3+lzd%llr(ilr)%d%n3i-1,lzd%glr%d%n3i)+1
                  !if (is3>i3 .or. i3>ie3) cycle
                  !if (.not.check_whether_bounds_overlap(is3,ie3,i3,i3)) cycle
                  if (.not.check_whether_bounds_overlap(is3,ie3,ii3,ii3)) cycle
                  !is1=1+lzd%Llr(ilr)%nsi1
                  !ie1=lzd%Llr(ilr)%nsi1+lzd%llr(ilr)%d%n1i
                  !is2=1+lzd%Llr(ilr)%nsi2
                  !ie2=lzd%Llr(ilr)%nsi2+lzd%llr(ilr)%d%n2i
                  is1=modulo(1+lzd%Llr(ilr)%nsi1-1,lzd%glr%d%n1i)+1
                  ie1=is1+lzd%llr(ilr)%d%n1i-1
                  is2=modulo(1+lzd%Llr(ilr)%nsi2-1,lzd%glr%d%n2i)+1
                  ie2=is2+lzd%llr(ilr)%d%n2i-1
                  !$omp parallel default(none) &
                  !$omp shared(is2, ie2, is1, ie1, weight, orbital_id, i3, iorb, iorbmin, iorbmax, lzd) &
                  !$omp private(i2, i1, ii2, ii1, jj)
                  !$omp do reduction(min:iorbmin) reduction(max:iorbmax)
                  do i2=is2,ie2
                      ii2=modulo(i2-1,lzd%glr%d%n2i)+1
                      do i1=is1,ie1
                          ii1=modulo(i1-1,lzd%glr%d%n1i)+1
                          jj=weight(ii1,ii2,i3)+1
                          !weight(i1,i2,i3) = weight(i1,i2,i3)+1
                          !orbital_id(weight(i1,i2,i3),i1,i2,i3) = iorb
                          orbital_id(jj,ii1,ii2,i3) = iorb
                          if (iorb<iorbmin) iorbmin=iorb
                          if (iorb>iorbmax) iorbmax=iorb
                          weight(ii1,ii2,i3)=jj
                      end do
                  end do
                  !$omp end do
                  !$omp end parallel
              end do
          end do
        
          ! Make sure that the bounds are okay for all processes
          if (iorbmin>iorbmax) then
              iorbmin=1
              iorbmax=1
          end if
        
        
          ! Now calculate the charge density. Of course this is only possible since the
          ! value of each gridpoint is given by the special pattern and therefore always known.
        
          ! First fill the kernel with some numbers.
          do i=1,denskern%nvctrp_tg*denskern%nspin
              denskern_%matrix_compr(i)=sine_taylor(real(denskern%nvctr*denskern%nspin-i+denskern%isvctrp_tg+1,kind=8))
              !denskern_%matrix_compr(i)=sine_taylor(real(mod(denskern%nspin*denskern%nvctr-i+1-1,denskern%nvctr)+1,kind=8))
              !write(660+iproc,'(a,2i8,2es13.5)') 'i, mod(denskern%nspin*denskern%nvctr-i+1-1,denskern%nvctr)+1, arg, val', &
              !     i, mod(denskern%nspin*denskern%nvctr-i+1-1,denskern%nvctr)+1, real(mod(denskern%nspin*denskern%nvctr-i+1-1,denskern%nvctr)+1,kind=8), denskern_%matrix_compr(i)
          end do
        
          hxh=.5d0*lzd%hgrids(1)
          hyh=.5d0*lzd%hgrids(2)
          hzh=.5d0*lzd%hgrids(3)
          factor=1.d0/(hxh*hyh*hzh)
        
          ! Use an auxilliary array to store the indices of the kernel in the compressed
          ! format. The usual denskern%matrixindex_in_compressed_fortransposed can not be used 
          ! since we are not in the transposed distribution.
          matrixindex_in_compressed_auxilliary=f_malloc((/iorbmin.to.iorbmax,iorbmin.to.iorbmax /), &
              id='matrixindex_in_compressed_auxilliary')
    
          !$omp parallel default(none) &
          !$omp shared(iorbmin, iorbmax, matrixindex_in_compressed_auxilliary, denskern) &
          !$omp private(iorb, jorb)
          !$omp do
          do iorb=iorbmin,iorbmax
              do jorb=iorbmin,iorbmax
                  matrixindex_in_compressed_auxilliary(jorb,iorb)=matrixindex_in_compressed(denskern, jorb, iorb)
              end do
          end do
          !$omp end do
          !$omp end parallel
        
          ! Now calculate the charge density and store the result in rho_check
          rho_check=f_malloc(max(lzd%glr%d%n1i*lzd%glr%d%n2i*(ii3e-ii3s+1)*denskern%nspin,1),id='rho_check')
          !!write(*,*) 'iproc, ii3s, ii3e', iproc, ii3s, ii3e
          do ispin=1,denskern%nspin
              ishift=(ispin-1)*lzd%glr%d%n1i*lzd%glr%d%n2i*(ii3e-ii3s+1)
              !!$omp parallel default (none) &
              !!$omp private (i3, i2, i1, ii3, iixyz, ind, tt, i,j, ii, tti, ikernel, jj, ttj) &
              !!$omp shared (ii3s, ii3e, lzd, weight, orbital_id, denskern, denskern_, rho_check) &
              !!$omp shared (nxyz, factor, matrixindex_in_compressed_auxilliary, ispin, orbs, ishift)
              do i3=ii3s,ii3e
                  ii3=modulo(i3-1,lzd%glr%d%n3i)+1
                  !!$omp do
                  do i2=1,lzd%glr%d%n2i
                      do i1=1,lzd%glr%d%n1i
                          !iixyz=(i3-1)*lzd%glr%d%n1i*lzd%glr%d%n2i+(i2-1)*lzd%glr%d%n1i+i1
                          iixyz=(ii3-1)*lzd%glr%d%n1i*lzd%glr%d%n2i+(i2-1)*lzd%glr%d%n1i+i1
                          ind=(i3-ii3s)*lzd%glr%d%n1i*lzd%glr%d%n2i+(i2-1)*lzd%glr%d%n1i+i1+ishift
                          tt=1.d-20
                          do i=1,weight(i1,i2,i3) !the number of orbitals touching this grid point
                              ii=orbital_id(i,i1,i2,i3)+(ispin-1)*orbs%norbu
                              !!iim=mod(ii-1,orbs%norbu)+1 !orbital number regardless of the spin
                              !!ispin=(ii-1)/orbs%norbu+1 !integer division to get the spin (1 for spin up (or non polarized), 2 for spin down)
                              tti=test_value_sumrho(ii,iixyz,nxyz)
                              !ikernel=matrixindex_in_compressed_auxilliary(ii,ii)
                              ikernel=matrixindex_in_compressed(denskern,ii,ii)-denskern%isvctrp_tg
                              tt=tt+denskern_%matrix_compr(ikernel)*tti*tti
                              do j=i+1,weight(i1,i2,i3)
                                  jj=orbital_id(j,i1,i2,i3)+(ispin-1)*orbs%norbu
                                  !ikernel=matrixindex_in_compressed_auxilliary(jj,ii)
                                  ikernel=matrixindex_in_compressed(denskern,jj,ii)-denskern%isvctrp_tg
                                  if (ikernel==0) cycle
                                  ttj=test_value_sumrho(jj,iixyz,nxyz)
                                  tt=tt+2.d0*denskern_%matrix_compr(ikernel)*tti*ttj
                                  !!if (mod(ind-1,lzd%glr%d%n1i*lzd%glr%d%n2i*(ii3e-ii3s+1))+1==865737) then
                                  !!    write(6500,'(a,6i8,3es13.5)') 'ind, i, j, ii, jj, ikernel, tti, ttj, valk', &
                                  !!        ind, i, j, ii, jj, ikernel, tti, ttj, denskern_%matrix_compr(ikernel)
                                  !!end if
                              end do
                          end do
                          tt=tt*factor
                          rho_check(ind)=tt
                          !!write(2100+iproc,'(a,4i9,es18.8)') 'ind, i1, i2, ii3, rho_check(ind)',ind, i1, i2, ii3, rho_check(ind)
                      end do
                  end do
                  !!$omp end do
              end do
              !!$omp end parallel
          end do
        
          call f_free(matrixindex_in_compressed_auxilliary)
        
          ! Now calculate the charge density in the transposed way using the standard routine
          rho=f_malloc(max(lzd%glr%d%n1i*lzd%glr%d%n2i*(ii3e-ii3s+1)*denskern%nspin,1),id='rho')
          !denskern_%matrix_compr = denskern%matrix_compr
          call sumrho_for_TMBs(iproc, nproc, lzd%hgrids(1), lzd%hgrids(2), lzd%hgrids(3), collcom_sr, denskern, denskern_, &
               denspot%dpbox%ndimrhopot, rho, rho_negative, .false.)
        
          ! Determine the difference between the two versions
          sumdiff=0.d0
          maxdiff=0.d0
          ii=0
          do ispin=1,denskern%nspin
              ii=(ispin-1)*lzd%glr%d%n1i*lzd%glr%d%n2i*(ii3e-ii3s+1)
              !!$omp parallel default(none) shared(ii,lzd,ii3e,ii3s,rho,rho_check,sumdiff,maxdiff) private(i,tt)
              !!$omp do reduction(+:sumdiff) reduction(max:maxdiff) 
              do i=1,lzd%glr%d%n1i*lzd%glr%d%n2i*(ii3e-ii3s+1)
                  !ii=ii+1
                  tt=abs(rho(ii+i)-rho_check(ii+i))
                  !!write(2000+iproc,'(a,2i9,4es18.8)') 'i,ii,rho(ii+i),rho_check(ii+1),diff,sumdiff',i,ii,rho(ii+i),rho_check(ii+i), tt, sumdiff
                  sumdiff = sumdiff + tt**2
                  if (tt>maxdiff) maxdiff=tt
              end do
              !!$omp end do
              !!$omp end parallel
          end do
        
          ! Reduce the results
          if (nproc>1) then
              call mpiallred(sumdiff, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
              call mpiallred(maxdiff, 1, mpi_max, comm=bigdft_mpi%mpi_comm)
          end if
        
          ! Get mean value for the sum
          sumdiff = sumdiff/(lzd%glr%d%n1i*lzd%glr%d%n2i*lzd%glr%d%n3i*denskern%nspin)
          sumdiff=sqrt(sumdiff)
        
          ! Print the results
          if (iproc==0) then
              call yaml_map('Tolerance for the following test',tol_calculation_mean,fmt='(1es25.18)')
              if (sumdiff>tol_calculation_mean) then
                 !call yaml_warning('CALCULATION ERROR: total difference of '//trim(yaml_toa(sumdiff,fmt='(1es25.18)')))
                 call f_err_throw('CALCULATION ERROR: total difference of '//trim(yaml_toa(sumdiff,fmt='(1es25.18)')), &
                      err_name='BIGDFT_RUNTIME_ERROR')
              else
                 call yaml_map('calculation check, error sum', sumdiff,fmt='(1es25.18)')
              end if
              call yaml_map('Tolerance for the following test',tol_calculation_max,fmt='(1es25.18)')
              if (sumdiff>tol_calculation_max) then
                 !call yaml_warning('CALCULATION ERROR: max difference of '//trim(yaml_toa(maxdiff,fmt='(1es25.18)')))
                 call f_err_throw('CALCULATION ERROR: max difference of '//trim(yaml_toa(maxdiff,fmt='(1es25.18)')), &
                      err_name='BIGDFT_RUNTIME_ERROR')
              else
                 call yaml_map('calculation check, error max', maxdiff,fmt='(1es25.18)')
              end if
          end if
        
        
          call f_free(weight)
          call f_free(orbital_id)
          call f_free(rho_check)
          call f_free(rho)
    
      end if
    
      if (iproc==0) call yaml_mapping_close()
    
      call timing(iproc,'check_sumrho','OF')
    
      call f_release_routine()
    
      contains
    
        function test_value_sumrho(i, j, n)
          implicit none
    
          ! Calling arguments
          integer,intent(in) :: i, j, n
          real(kind=8) :: test_value_sumrho
    
          ! Local variables
          real(kind=8) :: ri, rj, rn
          real(kind=8),parameter :: fac=1.d-8
    
          ri=real(i,kind=8)
          rj=real(j,kind=8)
          rn=real(n,kind=8)
          !test_value=fac*real((i-1)*n+j,dp)
          !test_value_sumrho=fac*(ri-1.d0)*rn+rj
          test_value_sumrho=sine_taylor((ri-1.d0)*rn)*cosine_taylor(rj)
          !test_value_sumrho=0.d0
    
    
        end function test_value_sumrho
    
        function sine_taylor(xx)
          implicit none
    
          ! Calling arguments
          real(kind=8),intent(in) :: xx
          real(kind=8) :: sine_taylor
    
          ! Local variables
          real(kind=8) :: x, x2, x3, x5, x7, x9, x11, x13, x15
          real(kind=8),parameter :: pi=3.14159265358979323846d0
          real(kind=8),parameter :: pi2=6.28318530717958647693d0
          real(kind=8),parameter :: inv6=1.66666666666666666667d-1
          real(kind=8),parameter :: inv120=8.33333333333333333333d-3
          real(kind=8),parameter :: inv5040=1.98412698412698412698d-4
          real(kind=8),parameter :: inv362880=2.75573192239858906526d-6
          real(kind=8),parameter :: inv39916800=2.50521083854417187751d-8
          real(kind=8),parameter :: inv6227020800=1.60590438368216145994d-10
          real(kind=8),parameter :: inv1307674368000=7.6471637318198164759d-13
    
          ! The Taylor approximation is most accurate around 0, so shift by pi to be centered around this point.
          ! This first part is equivalent to x=mod(xx,pi2)-pi
          x=xx/pi2
          x=real(int(x,kind=8),kind=8)*pi2
          x=xx-x-pi
    
          x2=x*x
          x3=x2*x
          x5=x3*x2
          x7=x5*x2
          x9=x7*x2
          x11=x9*x2
          x13=x11*x2
          x15=x13*x2
    
          ! Calculate the value
          sine_taylor = x - x3*inv6 + x5*inv120 - x7*inv5040 + x9*inv362880 &
                        - x11*inv39916800 + x13*inv6227020800 - x15*inv1307674368000
    
          ! Undo the shift of pi, which corresponds to a multiplication with -1
          sine_taylor=-1.d0*sine_taylor
    
        end function sine_taylor
    
        function cosine_taylor(xx)
          implicit none
    
          ! Calling arguments
          real(kind=8),intent(in) :: xx
          real(kind=8) :: cosine_taylor
    
          ! Local variables
          real(kind=8) :: x, x2, x4, x6, x8, x10, x12, x14
          real(kind=8),parameter :: pi=3.14159265358979323846d0
          real(kind=8),parameter :: pi2=6.28318530717958647693d0
          real(kind=8),parameter :: inv2=5.d-1
          real(kind=8),parameter :: inv24=4.16666666666666666667d-2
          real(kind=8),parameter :: inv720=1.38888888888888888889d-3
          real(kind=8),parameter :: inv40320=2.48015873015873015873d-5
          real(kind=8),parameter :: inv3628800=2.75573192239858906526d-7
          real(kind=8),parameter :: inv479001600=2.08767569878680989792d-9
          real(kind=8),parameter :: inv87178291200=1.14707455977297247139d-11
    
          ! The Taylor approximation is most accurate around 0, so shift by pi to be centered around this point.
          ! This first part is equivalent to x=mod(xx,pi2)-pi
          x=xx/pi2
          x=real(int(x,kind=8),kind=8)*pi2
          x=xx-x-pi
    
          x2=x*x
          x4=x2*x2
          x6=x4*x2
          x8=x6*x2
          x10=x8*x2
          x12=x10*x2
          x14=x12*x2
    
          ! Calculate the value
          cosine_taylor = 1.d0 - x2*inv2 + x4*inv24 - x6*inv720 + x8*inv40320 &
                          - x10*inv3628800 + x12*inv479001600 - x14*inv87178291200
    
          ! Undo the shift of pi, which corresponds to a multiplication with -1
          cosine_taylor=-1.d0*cosine_taylor
    
        end function cosine_taylor
    
    end subroutine check_communication_sumrho


    subroutine check_communications_locreg(iproc,nproc,orbs,nspin,Lzd,collcom,smat,mat,npsidim_orbs,npsidim_comp,check_overlap)
       use module_base!, only: wp, bigdft_mpi, mpi_sum, mpi_max, mpiallred
       use module_types, only: orbitals_data, local_zone_descriptors, linear_matrices
       use yaml_output
       use communications_base, only: comms_linear, TRANSPOSE_FULL
       use communications, only: transpose_localized, untranspose_localized
       use sparsematrix_base, only : sparse_matrix, matrices, DENSE_PARALLEL
       use sparsematrix, only : compress_matrix_distributed_wrapper, gather_matrix_from_taskgroups_inplace
       use transposed_operations, only: calculate_overlap_transposed
       use locregs, only: check_overlap_cubic_periodic
       use locreg_operations, only: Lpsi_to_global2
       !use dynamic_memory
       implicit none
       integer, intent(in) :: iproc,nproc,nspin,check_overlap
       type(orbitals_data), intent(in) :: orbs
       type(local_zone_descriptors), intent(in) :: lzd
       type(comms_linear), intent(in) :: collcom
       type(sparse_matrix),intent(in) :: smat
       type(matrices),intent(inout) :: mat
       integer, intent(in) :: npsidim_orbs, npsidim_comp
       !local variables
       character(len=*), parameter :: subname='check_communications'
       integer, parameter :: ilog=6
       integer :: i,ispinor,iorb,indspin,ikptsp
       integer :: ikpt,ierr,i0,ifine,ii,iiorb,ipt,jorb,indorb_tmp
       integer :: ispin
       !integer :: icomp
       !!$integer :: ipsi,ipsic,ipsif,ipsiworkc,ipsiworkf,jcomp,jkpt
       real(wp) :: psival,maxdiff,tt
       real(wp), dimension(:), allocatable :: psi,psit_c,psit_f
       real(wp), dimension(:,:), allocatable :: checksum
       real(wp) :: epsilon,tol
       logical :: abort, isoverlap
       integer :: jjorb, ilr, jlr, ldim, gdim, jjspin, ist, niorb, njorb, jjjorb, is, ie
       !integer :: iispin
       real(kind=8),dimension(:),allocatable :: psii, psij, psiig, psijg, mat_compr
       real(kind=8),dimension(:,:),allocatable :: matp
       real(kind=8) :: ddot
    
       call f_routine(id='check_communications_locreg')
    
       if (check_overlap > 0) then
    
           !allocate the "wavefunction" and fill it, and also the workspace
           psi = f_malloc(max(npsidim_orbs, npsidim_comp),id='psi')
           psit_c = f_malloc(sum(collcom%nrecvcounts_c),id='psit_c')
           psit_f = f_malloc(7*sum(collcom%nrecvcounts_f),id='psit_f')
           checksum = f_malloc0((/ orbs%norb*orbs%nspinor, 2 /),id='checksum')
           if (orbs%norbp>0) then
              tol=1.e-10*real(npsidim_orbs,wp)/real(orbs%norbp,wp)
           else
              tol=0.0_wp
           end if
        
           do iorb=1,orbs%norbp
              ikpt=(orbs%isorb+iorb-1)/orbs%norb+1
              indorb_tmp=ind_orb(iorb)
              do ispinor=1,orbs%nspinor
                 indspin=(ispinor-1)*nvctr_orb(iorb)+indorb_tmp
                 !checksum(orbs%isorb+iorb+(ispinor-1)*orbs%nspinor,1)=0.0_wp
                 tt=0.0_wp
                 do i=1,nvctr_orb(iorb)
                    !vali=real(i,wp)/512.0_wp  ! *1.d-5
                    call test_value_locreg(ikpt,orbs%isorb+iorb-(ikpt-1)*orbs%norb,ispinor,i,psival)
                  !psival=dble(mod(iorb+orbs%isorb-1,orbs%norbu)+1)*orbs%spinsgn(iorb+orbs%isorb)  
                    !psi(i+indspin+ind_orb(iorb))=psival!(valorb+vali)*(-1)**(ispinor-1)
                    !psi(i+indspin)=dble(iorb+orbs%isorb)*orbs%spinsgn(iorb+orbs%isorb)!psival!(valorb+vali)*(-1)**(ispinor-1)
                    psi(i+indspin)=psival!(valorb+vali)*(-1)**(ispinor-1)
                    tt=tt+psival
                    !checksum(orbs%isorb+iorb+(ispinor-1)*orbs%nspinor,1)=&
                    !     checksum(orbs%isorb+iorb+(ispinor-1)*orbs%nspinor,1)+psival
                 end do
                 checksum(orbs%isorb+iorb+(ispinor-1)*orbs%nspinor,1)=tt
              end do
           end do
        
           call transpose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, TRANSPOSE_FULL, &
                psi, psit_c, psit_f, lzd)
           !!do i=1,size(psit_c)
           !!    write(7000+iproc,*) i, psit_c(i)
           !!end do
           !!do i=1,size(psit_f)
           !!    write(7100+iproc,*) i, psit_f(i)
           !!end do
           
           !check the results of the transposed wavefunction
           maxdiff=0.0_wp
           if (iproc==0) call yaml_map('Number of coarse and fine DoF (MasterMPI task)',&
                (/collcom%nptsp_c,collcom%nptsp_f/),fmt='(i8)')
        
           do ikptsp=1,1!orbs%nkptsp !should be one for the moment
              ikpt=orbs%iskpts+ikptsp!orbs%ikptsp(ikptsp)
              ispinor=1 !for the (long?) moment
              !icomp=1
              do ispin=1,nspin
                 if (collcom%nptsp_c>0) then
                    do ipt=1,collcom%nptsp_c 
                       ii=collcom%norb_per_gridpoint_c(ipt)
                       i0 = collcom%isptsp_c(ipt) + (ispin-1)*collcom%ndimind_c/nspin
                       do i=1,ii
                          iiorb=collcom%indexrecvorbital_c(i0+i)
                          !write(5000+iproc,'(a,4i8,es16.6)') 'ispin, ipt, ii, i0+i, psit_c(i0+i)', ispin, ipt, ii, i0+i, psit_c(i0+i)
           !!$               !here a function which determin the address after mpi_alltoall
           !!$               !procedure should be called
           !!$               ipsitworkc=collcom%iexpand_c(icomp)
           !!$               !ipsiglob=collcom%nrecvdspls_c(iproc)+1+(ipsitworkc-1)*sum(
           !!$               ipsic=collcom%isendbuf_c(ipsiworkc)
           !!$               ipsi=ipsic
           !!$               do jorb=1,iiorb-1
           !!$                  ipsi=ipsi-nvctr_c_orb(jorb)
           !!$               end do
           !!$               call test_value_locreg(ikpt,iiorb-(ikpt-1)*orbs%norb,ispinor,&
           !!$                    ipsi,psival)
           !!$               indspin=(ispinor-1)*nvctr_orb(iiorb)
           !!$               maxdiff=max(abs(psit_c(i0+i)-psival),maxdiff)
                          checksum(iiorb,2)=checksum(iiorb,2)+psit_c(i0+i)
                          !icomp=icomp+1
                       end do
                    end do
                 end if
                 !icomp=1
                 if (collcom%nptsp_f>0) then
                    do ipt=1,collcom%nptsp_f 
                       ii=collcom%norb_per_gridpoint_f(ipt) 
                       i0 = collcom%isptsp_f(ipt) + (ispin-1)*collcom%ndimind_f/nspin
                       do i=1,ii
                          iiorb=collcom%indexrecvorbital_f(i0+i)
           !!$               ipsitworkf=collcom%iexpand_f(icomp)
           !!$               ipsif=collcom%isendbuf_f(ipsiworkf)
           !!$               ipsi=ipsif
           !!$               do jorb=1,iiorb-1
           !!$                  ipsi=ipsi-nvctr_f_orb(jorb)
           !!$               end do
                          tt=0.d0
                          do ifine=1,7
           !!$                  call test_value_locreg(ikpt,iiorb-(ikpt-1)*orbs%norb,ispinor,&
           !!$                       nvctr_c_orb(iiorb)+7*(ipsi-1)+ifine,psival) 
           !!$                  tt=abs(psit_f(7*(i0+i-1)+ifine)-psival)
           !!$                  if (tt > maxdiff) then
           !!$                     maxdiff=tt
           !!$                     !call wrong_components(psival,jkpt,jorb,jcomp)
           !!$                  end if
                             !checksum(iiorb,2)=checksum(iiorb,2)+psit_f(7*(i0+i-1)+ifine)
                             tt=tt+psit_f(7*(i0+i-1)+ifine)
                          end do
                          checksum(iiorb,2)=checksum(iiorb,2)+tt
                          !icomp=icomp+1
                       end do
                    end do
                 end if
              end do
           end do
        !!$
           if (iproc==0) then
              call yaml_map('Tolerances for this check',&
                (/tol,real(orbs%norb,wp)*epsilon(1.0_wp)/),fmt='(1pe25.17)')
           end if
        
           if (nproc > 1) then
              !call MPI_BARRIER(bigdft_mpi%mpi_comm, ierr)
              call mpiallred(checksum(1,1), 2*orbs%norb*orbs%nspinor, MPI_SUM, comm=bigdft_mpi%mpi_comm)
           end if
        
           !if (iproc==0) then
              maxdiff=0.0_wp
              do jorb=1,orbs%norb*orbs%nspinor
                 tt=abs(checksum(jorb,1)-checksum(jorb,2))
                 if (tt > maxdiff) then
                    maxdiff=tt
                    if (maxdiff > tol) then 
                       call yaml_warning('ERROR of checksum for orbital'//trim(yaml_toa(jorb))//&
                            ': difference of '//trim(yaml_toa(tt,fmt='(1pe12.5)')))
                    end if
                 end if
              end do
           !end if
           if (iproc==0) call yaml_map('Maxdiff for transpose (checksum)',&
                maxdiff,fmt='(1pe25.17)')
        
        
           abort = .false.
           if (abs(maxdiff) >tol) then
              call yaml_warning('ERROR (Transposition): process'//trim(yaml_toa(iproc))//&
                   ' found an error of:'//trim(yaml_toa(maxdiff,fmt='(1pe15.7)')))
              !call yaml_map('Some wrong results in',(/jkpt,jorb,jcomp/),fmt='(i8)')
              !abort=.true.
           end if
        
           if (abort) call MPI_ABORT(bigdft_mpi%mpi_comm,10,ierr)
        
        
           !@NEW: check the calculation of the overlap matrices #############
           if (check_overlap > 1) then
               call calculate_overlap_transposed(iproc, nproc, orbs, collcom, psit_c, &
                    psit_c, psit_f, psit_f, smat, mat)
               !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, smat, mat)
               !!do i=1,smat%nvctr*nspin
               !!    write(6000+iproc,'(a,2i8,es16.7)') 'i, mod(i-1,nvctr)+1, val', i, mod(i-1,smat%nvctr)+1, mat%matrix_compr(i)
               !!end do
               ! Alternative calculation of the overlap matrix
               gdim=lzd%glr%wfd%nvctr_c+7*lzd%glr%wfd%nvctr_f
               psiig = f_malloc(gdim,id='psiig')
               psijg = f_malloc(gdim,id='psijg')
               matp = f_malloc((/smat%nfvctr,smat%nfvctrp/),id='matp')
               mat_compr = f_malloc(smat%nvctr*smat%nspin,id='mat_compr')
               do ispin=1,smat%nspin
                   niorb=0
                   njorb=0
                   !not possible to iterate over norbp since the distributions over the MPI tasks might be incompatible with smat%nfvctrp
                   is=(ispin-1)*orbs%norbu+orbs%isorbu+1
                   ie=(ispin-1)*orbs%norbu+orbs%isorbu+orbs%norbup
                   !is=(ispin-1)*smat%nfvctr+smat%isfvctr+1
                   !is=(ispin-1)*smat%nfvctr+smat%isfvctr+smat%nfvctrp
                   do iiorb=is,ie
                       !iiorb=orbs%isorb+iorb
                       !if (orbs%spinsgn(iiorb)>0) then
                       !    iispin=1
                       !else
                       !    iispin=2
                       !end if
                       !if (iispin/=ispin) cycle
                       niorb=niorb+1
                       ilr=orbs%inwhichlocreg(iiorb)
                       ldim=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
                       psii = f_malloc(ldim,id='psii')
                       do i=1,ldim
                          call test_value_locreg(1,iiorb,1,i,psival)
                          psii(i)=psival
                       end do
                       call f_zero(psiig)
                       call Lpsi_to_global2(iproc, ldim, gdim, orbs%norb, orbs%nspinor, 1, lzd%glr, &
                                            lzd%llr(ilr), psii, psiig)
                       do jjorb=1,orbs%norb
                           if (orbs%spinsgn(jjorb)>0) then
                               jjspin=1
                           else
                               jjspin=2
                           end if
                           if (jjspin/=ispin) cycle
                           njorb=njorb+1
                           jjjorb=mod(jjorb-1,smat%nfvctr)+1 !index regardless of the spin
                           jlr=orbs%inwhichlocreg(jjorb)
                           ! check if there is an overlap, else cycle
                           call check_overlap_cubic_periodic(lzd%glr, lzd%llr(ilr), lzd%llr(jlr), isoverlap)
                           if (.not.isoverlap) then
                               matp(jjjorb,niorb)=0.d0
                               cycle
                           end if
                           ldim=lzd%llr(jlr)%wfd%nvctr_c+7*lzd%llr(jlr)%wfd%nvctr_f
                           psij = f_malloc(ldim,id='psij')
                           do i=1,ldim
                              call test_value_locreg(1,jjorb,1,i,psival)
                              psij(i)=psival
                           end do
                           call f_zero(psijg)
                           !!write(4200+iproc,'(a,2i8,l4)') 'iproc, jlr, associated(lzd%llr(jlr)%wfd%keygloc)', iproc, jlr, associated(lzd%llr(jlr)%wfd%keygloc)
                           call Lpsi_to_global2(iproc, ldim, gdim, orbs%norb, orbs%nspinor, 1, lzd%glr, &
                                                lzd%llr(jlr), psij, psijg)
                           matp(jjjorb,niorb)=ddot(gdim, psiig, 1, psijg, 1)
                           call f_free(psij)
                       end do
                       call f_free(psii)
                   end do
                   !ist=(ispin-1)*smat%nvctr+smat%isvctrp_tg+1
                   ist=(ispin-1)*smat%nvctrp_tg+smat%isvctrp_tg+1
                   call compress_matrix_distributed_wrapper(iproc, nproc, smat, DENSE_PARALLEL, &
                        matp, mat_compr(ist:))
               end do
               maxdiff=0.d0
               call f_free(psiig)
               call f_free(psijg)
               call f_free(matp)
               !write(*,'(3(a,i0))') 'task ',iproc,' checks the values from ',1+smat%isvctrp_tg,' to ',smat%nvctrp_tg+smat%isvctrp_tg
               !write(*,'(3(a,i0))') 'task ',iproc,' checks the values from ',smat%istartend_local(1),' to ',smat%istartend_local(2)
               !do i=1,smat%nvctrp_tg
               do i=smat%istartend_local(1),smat%istartend_local(2)
                   !maxdiff=max(abs(mat_compr(i+smat%isvctrp_tg)-mat%matrix_compr(i)),maxdiff)
                   maxdiff=max(abs(mat_compr(i)-mat%matrix_compr(i-smat%isvctrp_tg)),maxdiff)
                   !write(8000+iproc,'(a,i7,2es15.5)') 'i, mat_compr(i), mat%matrix_compr(i)', &
                   !    i, mat_compr(i), mat%matrix_compr(i)
               end do
               if (nproc>1) then
                   call mpiallred(maxdiff, 1, mpi_max, comm=bigdft_mpi%mpi_comm)
               end if
               call f_free(mat_compr)
               if (iproc==0) call yaml_map('Maxdiff for overlap calculation',maxdiff,fmt='(1es25.17)')
           end if
           !@END NEW ########################################################
        
        
           call untranspose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, &
                TRANSPOSE_FULL, psit_c, psit_f, psi, lzd)
        
           maxdiff=0.0_wp
           do iorb=1,orbs%norbp
              ikpt=(orbs%isorb+iorb-1)/orbs%norb+1
              do ispinor=1,orbs%nspinor
                 indspin=(ispinor-1)*nvctr_orb(iorb)
                 do i=1,nvctr_orb(iorb)
                    call test_value_locreg(ikpt,orbs%isorb+iorb-(ikpt-1)*orbs%norb,ispinor,i,psival)
                    maxdiff=max(abs(psi(i+indspin+ind_orb(iorb))-psival),maxdiff)
                 end do
              end do
           end do
        
        
           abort = .false.
           if (abs(maxdiff) > real(orbs%norb,wp)*epsilon(1.0_wp)) then
              call yaml_warning('ERROR (Inverse Transposition): process'//trim(yaml_toa(iproc))//&
                   ' found an error of:'//trim(yaml_toa(maxdiff,fmt='(1pe15.7)')))
              !abort = .true.
           end if
        
       if (abort) call MPI_ABORT(bigdft_mpi%mpi_comm,11,ierr)
        
           if (nproc > 1) then
              !call MPI_BARRIER(bigdft_mpi%mpi_comm, ierr)
              call mpiallred(maxdiff, 1, MPI_MAX, comm=bigdft_mpi%mpi_comm)
           end if
        
           if (iproc==0) call yaml_map('Maxdiff for untranspose',maxdiff,fmt='(1pe25.17)')
        
           call f_free(psi)
           call f_free(psit_c)
           call f_free(psit_f)
           call f_free(checksum)
        
       end if
    
       call f_release_routine()
    
     contains
       
    
       function ind_orb(iorb)
         implicit none
         integer, intent(in) :: iorb
         integer :: ind_orb
         !local variables
         integer :: jorb
         ind_orb=0
         do jorb=1,iorb-1
            ind_orb=ind_orb+nvctr_orb(jorb)
         end do
       end function ind_orb
    
       function nvctr_orb(iorb)
         implicit none
         integer, intent(in) :: iorb
         integer :: nvctr_orb
         !local variables
         integer :: jlr
    
         jlr = orbs%inwhichlocreg(iorb+orbs%isorb)
         nvctr_orb=(Lzd%Llr(jlr)%wfd%nvctr_c+7*Lzd%Llr(jlr)%wfd%nvctr_f)
         
       end function nvctr_orb
    
       function nvctr_c_orb(iorb)
         implicit none
         integer, intent(in) :: iorb
         integer :: nvctr_c_orb
         !local variables
         integer :: jlr
    
         jlr = orbs%inwhichlocreg(iorb+orbs%isorb)
         nvctr_c_orb=Lzd%Llr(jlr)%wfd%nvctr_c
         
       end function nvctr_c_orb
    
       function nvctr_f_orb(iorb)
         implicit none
         integer, intent(in) :: iorb
         integer :: nvctr_f_orb
         !local variables
         integer :: jlr
    
         jlr = orbs%inwhichlocreg(iorb+orbs%isorb)
         nvctr_f_orb=Lzd%Llr(jlr)%wfd%nvctr_f
         
       end function nvctr_f_orb
    
    
       !> define a value for the wavefunction which is dependent of the indices
       subroutine test_value_locreg(ikpt,iorb,ispinor,icomp,val)
         use module_base
         implicit none
         integer, intent(in) :: ikpt,icomp,iorb,ispinor
         real(wp), intent(out) :: val
         !local variables
         real(wp) :: valkpt,valorb,vali
    
         ! recognizable pattern, for debugging
         valkpt=real(10**ilog*(ikpt-1),wp)!real(512*ikpt,wp)
         valorb=real(iorb,wp)+valkpt
         vali=real(icomp,wp)*10.0_wp**(-ilog)  !real(icomp,wp)/512.0_wp  ! *1.d-5
         val=(valorb+vali)*(-1)**(ispinor-1)
    
       END SUBROUTINE test_value_locreg
    
       !>determine the components which were not communicated correctly
       !! works only with the recognizable pattern of test function
       subroutine wrong_components_locreg(psival,ikpt,iorb,icomp)
         use module_base
         implicit none
         real(wp), intent(in) :: psival
         integer, intent(out) :: ikpt,iorb,icomp
    
         icomp=nint((psival-real(floor(psival),wp))*10.0_wp**ilog)
         ikpt=floor(psival)/(10**ilog)
         iorb=floor(psival)-(ikpt-1)*(10**ilog)
    
       end subroutine wrong_components_locreg
    
    
    
     END SUBROUTINE check_communications_locreg

end module unitary_tests
