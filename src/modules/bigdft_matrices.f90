module bigdft_matrices
  implicit none

  private

  public :: check_local_matrix_extents
  public :: get_modulo_array
  public :: init_matrixindex_in_compressed_fortransposed


  contains

    subroutine check_local_matrix_extents(iproc, nproc, collcom, collcom_sr, smmd, smat, &
               ind_min, ind_max, ind_trans_min, ind_trans_max, irow, icol)
          use module_base
          use sparsematrix_base, only: sparse_matrix, sparse_matrix_metadata
          use sparsematrix_init, only: matrixindex_in_compressed, check_matmul_layout, &
                                       check_compress_distributed_layout
          use communications_base, only: comms_linear
          implicit none
    
          ! Caling arguments
          integer,intent(in) :: iproc, nproc
          type(comms_linear),intent(in) :: collcom, collcom_sr
          type(sparse_matrix_metadata),intent(in) :: smmd
          type(sparse_matrix),intent(in) :: smat
          integer,intent(out) :: ind_min, ind_max, ind_trans_min, ind_trans_max
          integer,dimension(2),intent(out) :: irow, icol
    
          ! Local variables
          integer :: i, ii_ref, iorb, jorb, ii, iseg
          logical :: found

          real(kind=4) :: tr0, tr1, trt0, trt1
          real(kind=8) :: time0, time1, time2, time3, time4, time5, ttime
          logical, parameter :: extra_timing=.false.
          integer,dimension(:),pointer :: moduloarray
                        
    
          call timing(iproc,'matrix_extents','ON')
          if (extra_timing) call cpu_time(trt0)  

          ind_min = smat%nvctr
          ind_max = 0

          if (extra_timing) call cpu_time(tr0)
          ! The operations done in the transposed wavefunction layout
          !call check_transposed_layout()
          call get_modulo_array(smat, moduloarray)
          call find_minmax_transposed(smat%matrixindex_in_compressed_fortransposed,collcom,smat%nfvctr,moduloarray,ind_min,ind_max)
          ind_trans_min = ind_min
          ind_trans_max = ind_max

          !write(*,'(a,2i8)') 'after check_transposed_layout: ind_min, ind_max', ind_min, ind_max
          if (extra_timing) call cpu_time(tr1)
          if (extra_timing) time0=real(tr1-tr0,kind=8)    


          ! Now check the compress_distributed layout
          !call check_compress_distributed_layout()
          call check_compress_distributed_layout(smat,ind_min,ind_max)

          !write(*,'(a,2i8)') 'after check_compress_distributed_layout: ind_min, ind_max', ind_min, ind_max
          if (extra_timing) call cpu_time(tr0)
          if (extra_timing) time1=real(tr0-tr1,kind=8)        

          ! Now check the matrix matrix multiplications layout
          if (smat%smatmul_initialized) then
              call check_matmul_layout(smat%smmm%nseq,smat%smmm%indices_extract_sequential,ind_min,ind_max)
          end if
          !write(*,'(a,2i8)') 'after check_matmul_layout: ind_min, ind_max', ind_min, ind_max
          if (extra_timing) call cpu_time(tr1)
          if (extra_timing) time2=real(tr1-tr0,kind=8)    
    
          ! Now check the sumrho operations
          !call check_sumrho_layout()
          call check_sumrho_layout(collcom_sr,smat%nfvctr,moduloarray,smat%matrixindex_in_compressed_fortransposed,ind_min,ind_max)

          call f_free_ptr(moduloarray)
          !write(*,'(a,2i8)') 'after check_sumrho_layout: ind_min, ind_max', ind_min, ind_max
          if (extra_timing) call cpu_time(tr0)
          if (extra_timing) time3=real(tr0-tr1,kind=8)    
    
          ! Now check the pseudo-exact orthonormalization during the input guess
          !call check_ortho_inguess()
          call check_ortho_inguess(smat,ind_min,ind_max)
          !write(*,'(a,2i8)') 'after check_ortho_inguess: ind_min, ind_max', ind_min, ind_max
          if (extra_timing) call cpu_time(tr1)
          if (extra_timing) time4=real(tr1-tr0,kind=8)        

          ! Now check the submatrix extraction for the projector charge analysis
          call check_projector_charge_analysis(iproc, nproc, smmd, smat, ind_min, ind_max)

          !!write(*,'(a,3i8)') 'after check_local_matrix_extents: iproc, ind_min, ind_max', iproc, ind_min, ind_max

          ! Get the global indices of ind_min and ind_max
          do i=1,2
              if (i==1) then
                  ii_ref = ind_min
              else
                  ii_ref = ind_max
              end if
              ! Search the indices iorb,jorb corresponding to ii_ref
              found=.false.

              ! not sure if OpenMP is really worth it here
              !$omp parallel default(none) &
              !$omp private(iseg,ii,iorb,jorb) &
              !$omp shared(smat,ii_ref,irow,icol,found,i)
              !$omp do
              outloop: do iseg=1,smat%nseg
                  if (.not. found) then
                     iorb = smat%keyg(1,2,iseg)
                     do jorb=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
                         ii = matrixindex_in_compressed(smat, jorb, iorb)
                         !if (iproc==0) write(*,'(a,5i9)') 'i, ii_ref, ii, iorb, jorb', i, ii_ref, ii, iorb, jorb
                         if (ii==ii_ref) then
                             irow(i) = jorb
                             icol(i) = iorb
                             !exit outloop
                             !SM: I think one should do this within a critical section since it is shared, just to be sure...
                             !$omp critical
                             found=.true.
                             !$omp end critical
                         end if
                     end do
                  end if
              end do outloop
              !$omp end do
              !$omp end parallel

          end do
          if (extra_timing) call cpu_time(tr0)
          if (extra_timing) time5=real(tr0-tr1,kind=8)    
          if (extra_timing) call cpu_time(trt1)  
          if (extra_timing) ttime=real(trt1-trt0,kind=8)  

          if (extra_timing.and.iproc==0) print*,'matextent',time0,time1,time2,time3,time4,time5,&
               time0+time1+time2+time3+time4+time5,ttime

          call timing(iproc,'matrix_extents','OF')    
    
    end subroutine check_local_matrix_extents


    subroutine find_minmax_transposed(matrixindex_in_compressed_fortransposed,collcom,nfvctr,moduloarray,ind_min,ind_max)
      use communications_base, only: comms_linear
      implicit none
      integer, intent(in) :: nfvctr
      type(comms_linear),intent(in) :: collcom
      integer, dimension(:,:), intent(in) :: matrixindex_in_compressed_fortransposed
      integer, dimension(nfvctr), intent(in) :: moduloarray
      integer, intent(inout) :: ind_min,ind_max
      !local variables
      integer :: ipt, ii, i0, i, i0i, iiorb, j, i0j, jjorb, ind, iorb, jorb

      !$omp parallel default(none) &
      !$omp private(ipt,ii,i0,i,i0i,iiorb,iorb,j,i0j,jjorb,jorb,ind) &
      !$omp shared(collcom,moduloarray,ind_min,ind_max,matrixindex_in_compressed_fortransposed,nfvctr)
      !$omp do reduction(min: ind_min) reduction(max: ind_max)
      do ipt=1,collcom%nptsp_c
         ii=collcom%norb_per_gridpoint_c(ipt)
         i0 = collcom%isptsp_c(ipt)
         do i=1,ii
            i0i=i0+i
            iiorb=collcom%indexrecvorbital_c(i0i)
            iorb=moduloarray(iiorb)
            do j=1,ii
               i0j=i0+j
               jjorb=collcom%indexrecvorbital_c(i0j)
               jorb=moduloarray(jjorb)
               !ind = smat%matrixindex_in_compressed_fortransposed(jorb,iorb)
               ind = matrixindex_in_compressed_fortransposed(jorb,iorb)
               ind_min = min(ind_min,ind)
               ind_max = max(ind_max,ind)
            end do
         end do
      end do
      !$omp end do

      !$omp do reduction(min: ind_min) reduction(max: ind_max)
      do ipt=1,collcom%nptsp_f
         ii=collcom%norb_per_gridpoint_f(ipt)
         i0 = collcom%isptsp_f(ipt)
         do i=1,ii
            i0i=i0+i
            iiorb=collcom%indexrecvorbital_f(i0i)
            iorb=moduloarray(iiorb)
            do j=1,ii
               i0j=i0+j
               jjorb=collcom%indexrecvorbital_f(i0j)
               jorb=moduloarray(jjorb)
               !ind = smat%matrixindex_in_compressed_fortransposed(jorb,iorb)
               ind = matrixindex_in_compressed_fortransposed(jorb,iorb)
               ind_min = min(ind_min,ind)
               ind_max = max(ind_max,ind)
            end do
         end do
      end do
      !$omp end do
      !$omp end parallel

    end subroutine find_minmax_transposed




    !> Copied from projector_for_charge_analysis and extract_matrix
    subroutine check_projector_charge_analysis(iproc, nproc, smmd, smat, ind_min, ind_max)
      use module_base
      use sparsematrix_base, only: sparse_matrix, sparse_matrix_metadata
      use sparsematrix_init, only: matrixindex_in_compressed
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(sparse_matrix_metadata),intent(in) :: smmd
      type(sparse_matrix),intent(in) :: smat
      integer,intent(inout) :: ind_min, ind_max

      integer :: ii, natp, jj, isat, kat, iatold, kkat, i, iat, j, ind
      integer,dimension(:,:),allocatable :: orbs_atom_id
      integer,dimension(:),allocatable :: neighbor_id
      integer,parameter :: ntmb_max = 16 !maximal number of TMBs per atom

      ! Parallelization over the number of atoms
      ii = smmd%nat/nproc
      natp = ii
      jj = smmd%nat - nproc*natp
      if (iproc<jj) then
          natp = natp + 1
      end if
      isat = (iproc)*ii + min(iproc,jj)

      orbs_atom_id = f_malloc0((/0.to.ntmb_max,1.to.natp/),id='orbs_atom_id')
      do i=1,smat%nfvctr
          kkat = smmd%on_which_atom(i)
          if (kkat>isat .and. kkat<=isat+natp) then
              kat = kkat - isat
              orbs_atom_id(0,kat) = orbs_atom_id(0,kat) + 1
              orbs_atom_id(orbs_atom_id(0,kat),kat) = i
          end if
      end do

      neighbor_id = f_malloc(0.to.smat%nfvctr,id='neighbor_id')
      do kat=1,natp
          ! Determine the "neighbors"
          iatold = 0
          kkat = kat + isat
          neighbor_id(0) = 0
          do ii=1,orbs_atom_id(0,kat)
              i = orbs_atom_id(ii,kat)
              do j=1,smat%nfvctr
                  ind =  matrixindex_in_compressed(smat, j, i)
                  if (ind/=0) then
                     neighbor_id(0) = neighbor_id(0) + 1
                     neighbor_id(neighbor_id(0)) = j
                  end if
              end do
          end do

          ! Determine the size of the matrix needed
          do ii=1,neighbor_id(0)
              i = neighbor_id(ii)
              do jj=1,neighbor_id(0)
                  j = neighbor_id(jj)
                  ind =  matrixindex_in_compressed(smat, j, i)
                  if (ind>0) then
                      ind_min = min(ind_min,ind)
                      ind_max = max(ind_max,ind)
                  end if
              end do
          end do
      end do

      call f_free(orbs_atom_id)
      call f_free(neighbor_id)

    end subroutine check_projector_charge_analysis


    subroutine check_sumrho_layout(collcom_sr,nfvctr,moduloarray,matrixindex_in_compressed_fortransposed,ind_min,ind_max)
      use communications_base, only: comms_linear
      use sparsematrix_init, only: matrixindex_in_compressed
      implicit none
      integer, intent(in) :: nfvctr
      type(comms_linear),intent(in) :: collcom_sr
      integer, dimension(:,:), intent(in) :: matrixindex_in_compressed_fortransposed
      integer, dimension(nfvctr), intent(in) :: moduloarray
      integer, intent(inout) :: ind_min,ind_max
      !local variables
      integer :: ipt, ii, i0, i, iiorb, ind, iorb

      !$omp parallel default(none) &
      !$omp private(ipt,ii,i0,iiorb,iorb,ind,i) &
      !$omp shared(collcom_sr,moduloarray,matrixindex_in_compressed_fortransposed,ind_min,ind_max)
      !$omp do reduction(min: ind_min) reduction(max: ind_max)
      do ipt=1,collcom_sr%nptsp_c
         ii=collcom_sr%norb_per_gridpoint_c(ipt)
         i0=collcom_sr%isptsp_c(ipt)
         do i=1,ii
            iiorb=collcom_sr%indexrecvorbital_c(i0+i)
            iorb=moduloarray(iiorb)
            !ind=smat%matrixindex_in_compressed_fortransposed(iiorb,iiorb)
            ind=matrixindex_in_compressed_fortransposed(iorb,iorb)
            !ind=get_transposed_index(smat,iiorb,iiorb)
            ind_min = min(ind_min,ind)
            ind_max = max(ind_max,ind)
         end do
      end do
      !$omp end do
      !$omp end parallel
    end subroutine check_sumrho_layout

    subroutine check_ortho_inguess(smat,ind_min,ind_max)
      use module_base
      use sparsematrix_base, only: sparse_matrix
      use sparsematrix_init, only: matrixindex_in_compressed
      implicit none
      type(sparse_matrix),intent(in) :: smat
      integer, intent(inout) :: ind_min,ind_max
      !local variables
      integer :: iorb, iiorb, isegstart, isegend, iseg, j, i, jorb, korb, ind, nthread, ithread
      logical, dimension(:,:), allocatable :: in_neighborhood
      !$ integer :: omp_get_max_threads, omp_get_thread_num

      !call f_routine(id='check_ortho_inguess')

      ! Allocate the array for all threads to avoid that it has to be declared private
      nthread = 1
      !$ nthread = omp_get_max_threads()
      in_neighborhood = f_malloc((/1.to.smat%nfvctr,0.to.nthread-1/),id='in_neighborhood')

      ithread = 0
      !$omp parallel default(none) &
      !$omp shared(smat, in_neighborhood, ind_min, ind_max) &
      !$omp private(iorb, iiorb, isegstart, isegend, iseg, j, jorb, korb, ind,i) &
      !$omp firstprivate(ithread)
      !$omp do reduction(min: ind_min) reduction(max: ind_max)
      do iorb=1,smat%nfvctrp
         !$ ithread = omp_get_thread_num()

         iiorb = smat%isfvctr + iorb
         isegstart = smat%istsegline(iiorb)
         isegend = smat%istsegline(iiorb) + smat%nsegline(iiorb) -1
         in_neighborhood(:,ithread) = .false.
         do iseg=isegstart,isegend
            ! A segment is always on one line, therefore no double loop
            j = smat%keyg(1,2,iseg)
            do i=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
               in_neighborhood(i,ithread) = .true.
            end do
         end do

         do jorb=1,smat%nfvctr
            if (.not.in_neighborhood(jorb,ithread)) cycle
            do korb=1,smat%nfvctr
               if (.not.in_neighborhood(korb,ithread)) cycle
               ind = matrixindex_in_compressed(smat,korb,jorb)
               if (ind>0) then
                  ind_min = min(ind_min,ind)
                  ind_max = max(ind_max,ind)
               end if
            end do
         end do

      end do
      !$omp end do
      !$omp end parallel

      call f_free(in_neighborhood)


      !call f_release_routine()

    end subroutine check_ortho_inguess


    subroutine get_modulo_array(smat, moduloarray)
      use module_base
      use sparsematrix_base, only: sparse_matrix
      implicit none
      ! Calling arguments
      type(sparse_matrix),intent(in) :: smat
      integer,dimension(:),pointer :: moduloarray
      ! Local variables
      integer :: i
      moduloarray = f_malloc_ptr(smat%nfvctr,id='moduloarray')
      !$omp parallel default(none) &
      !$omp shared(moduloarray,smat) &
      !$omp private(i)
      !$omp do
      do i=1,smat%nfvctr
          moduloarray(i) = modulo(i-smat%offset_matrixindex_in_compressed_fortransposed,smat%nfvctr)+1
      end do
      !$omp end do
      !$omp end parallel
    end subroutine get_modulo_array


    subroutine init_matrixindex_in_compressed_fortransposed(iproc, nproc, collcom, collcom_shamop, &
               collcom_sr, sparsemat)
      use module_base
      use sparsematrix_base, only: sparse_matrix
      use communications_base, only: comms_linear
      use sparsematrix_init, only: matrixindex_in_compressed
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(comms_linear),intent(in) :: collcom, collcom_shamop, collcom_sr
      type(sparse_matrix),intent(inout) :: sparsemat
      
      ! Local variables
      integer :: iorb, jorb, istat, imin, imax, nmiddle, imax_old, imin_old, iiorb, jjorb
      integer :: ii, imin_new, imax_new, i, nlen, j
      !integer :: kproc,jproc,jjorbold,jjorb,isend,irecv,ilr,ijorb,iiorb,ind,ierr, irow, irowold, iseg
      !integer :: compressed_index
    !  integer,dimension(:,:),allocatable :: sendbuf, requests, iminmaxarr
      character(len=*),parameter :: subname='init_sparse_matrix'
    
      call f_routine(id='init_matrixindex_in_compressed_fortransposed')
      call timing(iproc,'init_matrCompr','ON')
    
      ! for the calculation of overlaps and the charge density
      !imin=minval(collcom%indexrecvorbital_c)
      !imin=min(imin,minval(collcom%indexrecvorbital_f))
      !imin=min(imin,minval(collcom_shamop%indexrecvorbital_c))
      !imin=min(imin,minval(collcom_shamop%indexrecvorbital_f))
      !imin=min(imin,minval(collcom_sr%indexrecvorbital_c))
      !imax=maxval(collcom%indexrecvorbital_c)
      !imax=max(imax,maxval(collcom%indexrecvorbital_f))
      !imax=max(imax,maxval(collcom_shamop%indexrecvorbital_c))
      !imax=max(imax,maxval(collcom_shamop%indexrecvorbital_f))
      !imax=max(imax,maxval(collcom_sr%indexrecvorbital_c))
    
      nmiddle = sparsemat%nfvctr/2 + 1
    
      imin_old = huge(1)
      imax_old = 0
      imin_new = huge(1)
      imax_new = 0
      do i=1,size(collcom%indexrecvorbital_c)
          ii = mod(collcom%indexrecvorbital_c(i)-1,sparsemat%nfvctr)+1
          imin_old = min(imin_old,ii)
          imax_old = max(imax_old,ii)
          if (ii>nmiddle) then
              imin_new = min(imin_new,ii)
          else
              imax_new = max(imax_new,ii+sparsemat%nfvctr)
          end if
      end do
      do i=1,size(collcom%indexrecvorbital_f)
          ii = mod(collcom%indexrecvorbital_f(i)-1,sparsemat%nfvctr)+1
          imin_old = min(imin_old,ii)
          imax_old = max(imax_old,ii)
          if (ii>nmiddle) then
              imin_new = min(imin_new,ii)
          else
              imax_new = max(imax_new,ii+sparsemat%nfvctr)
          end if
      end do
      do i=1,size(collcom_shamop%indexrecvorbital_c)
          ii = mod(collcom_shamop%indexrecvorbital_c(i)-1,sparsemat%nfvctr)+1
          imin_old = min(imin_old,ii)
          imax_old = max(imax_old,ii)
          if (ii>nmiddle) then
              imin_new = min(imin_new,ii)
          else
              imax_new = max(imax_new,ii+sparsemat%nfvctr)
          end if
      end do
      do i=1,size(collcom_shamop%indexrecvorbital_f)
          ii = mod(collcom_shamop%indexrecvorbital_f(i)-1,sparsemat%nfvctr)+1
          imin_old = min(imin_old,ii)
          imax_old = max(imax_old,ii)
          if (ii>nmiddle) then
              imin_new = min(imin_new,ii)
          else
              imax_new = max(imax_new,ii+sparsemat%nfvctr)
          end if
      end do
      do i=1,size(collcom_sr%indexrecvorbital_c)
          ii = mod(collcom_sr%indexrecvorbital_c(i)-1,sparsemat%nfvctr)+1
          imin_old = min(imin_old,ii)
          imax_old = max(imax_old,ii)
          if (ii>nmiddle) then
              imin_new = min(imin_new,ii)
          else
              imax_new = max(imax_new,ii+sparsemat%nfvctr)
          end if
      end do
    
    
      !!write(*,*) 'iproc, imin_old, imax_old', iproc, imin_old, imax_old
      !!write(*,*) 'iproc, imin_new, imax_new', iproc, imin_new, imax_new
    
      !! values regardless of the spin
      !imin=mod(imin-1,sparsemat%nfvctr)+1
      !imax=mod(imax-1,sparsemat%nfvctr)+1
    
    
      ! Determine with which size the array should be allocated
      if (imax_new-imin_new<0) then
          ! everything in either first or second half
          imin = imin_old
          imax = imax_old
          !sparsemat%offset_matrixindex_in_compressed_fortransposed = 1
      else
          ! in both half
          if (imax_old-imin_old>imax_new-imin_new) then
              ! wrap around
              imin = imin_new
              imax = imax_new
              !sparsemat%offset_matrixindex_in_compressed_fortransposed = imin_new
          else
              ! no wrap around
              imin = imin_old
              imax = imax_old
              !sparsemat%offset_matrixindex_in_compressed_fortransposed = 1
          end if
      end if
    
      !!! Check
      !!if (sparsemat%offset_matrixindex_in_compressed_fortransposed<sparsemat%nfvctr/2+1) then
      !!    stop 'sparsemat%offset_matrixindex_in_compressed_fortransposed<sparsemat%nfvctr/2+1'
      !!end if
    
      nlen = imax - imin + 1
      sparsemat%offset_matrixindex_in_compressed_fortransposed = imin
      !!write(*,*) 'iproc, imin, imax, nlen', iproc, imin, imax, nlen
    
      !!! This is a temporary solution for spin polarized systems
      !!imax=min(imax,orbs%norbu)
    
    
    
      !!allocate(sparsemat%matrixindex_in_compressed_fortransposed(imin:imax,imin:imax), stat=istat)
      !!call memocc(istat, sparsemat%matrixindex_in_compressed_fortransposed, &
      !sparsemat%matrixindex_in_compressed_fortransposed=f_malloc_ptr((/imin.to.imax,imin.to.imax/),&
      !    id='sparsemat%matrixindex_in_compressed_fortransposed')
      sparsemat%matrixindex_in_compressed_fortransposed=f_malloc_ptr((/nlen,nlen/),&
          id='sparsemat%matrixindex_in_compressed_fortransposed')
    
      !$omp parallel do default(private) shared(sparsemat,imin,imax)
      do iorb=imin,imax
          i = iorb - imin + 1
          do jorb=imin,imax
              j = jorb - imin + 1
              !@ii=(jorb-1)*sparsemat%nfvctr+iorb
              !@ispin=(ii-1)/sparsemat%nfvctr+1 !integer division to get the spin (1 for spin up (or non polarized), 2 for spin down)
              !@iiorb=mod(iorb-1,sparsemat%nfvctr)+1 !orbital number regardless of the spin
              !@jjorb=mod(jorb-1,sparsemat%nfvctr)+1 !orbital number regardless of the spin
              !sparsemat%matrixindex_in_compressed_fortransposed(iorb,jorb)=compressed_index(iiorb,jjorb,orbs%norbu,sparsemat)
              iiorb = mod(iorb-1,sparsemat%nfvctr)+1
              jjorb = mod(jorb-1,sparsemat%nfvctr)+1
              !sparsemat%matrixindex_in_compressed_fortransposed(iorb,jorb)=matrixindex_in_compressed(sparsemat, iiorb, jjorb)
              sparsemat%matrixindex_in_compressed_fortransposed(i,j)=matrixindex_in_compressed(sparsemat, iiorb, jjorb)
              !sendbuf(jorb,iorb)=compressed_index(jorb,iiorb,orbs%norb,sparsemat)
              !sendbuf(iorb,jorb)=compressed_index(iiorb,jorb,orbs%norb,sparsemat)
          end do
      end do
      !$omp end parallel do
    
      !@! Add the spin shift (i.e. the index is in the spin polarized matrix which is at the end)
      !@if (ispin==2) then
      !@    matrixindex_in_compressed = matrixindex_in_compressed + sparsemat%nvctr
      !@end if
    
      call timing(iproc,'init_matrCompr','OF')
      call f_release_routine()
    
    end subroutine init_matrixindex_in_compressed_fortransposed

end module bigdft_matrices
