module transposed_operations

  implicit none

  private

  public :: calculate_overlap_transposed
  public :: build_linear_combination_transposed
  public :: normalize_transposed
  public :: init_matrixindex_in_compressed_fortransposed


  contains


    subroutine calculate_overlap_transposed(iproc, nproc, orbs, collcom, &
               psit_c1, psit_c2, psit_f1, psit_f2, smat, ovrlp)
      use module_base
      use module_types, only: orbitals_data
      use communications_base, only: comms_linear
      use sparsematrix_base, only: sparse_matrix, matrices
      use sparsematrix_init, only: get_modulo_array
      use sparsematrix, only: synchronize_matrix_taskgroups
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(orbitals_data),intent(in) :: orbs
      type(comms_linear),intent(in) :: collcom
      real(kind=8),dimension(collcom%ndimind_c),intent(in) :: psit_c1, psit_c2
      real(kind=8),dimension(7*collcom%ndimind_f),intent(in) :: psit_f1, psit_f2
      type(sparse_matrix),intent(inout) :: smat
      type(matrices),intent(inout) :: ovrlp
    
      ! Local variables
      integer :: i0, ipt, ii, iiorb, j, jjorb, i, ierr, istat, m, tid, norb, nthreads, ispin, ishift_mat
      integer :: istart, iend, orb_rest, ind0, ind1, ind2, ind3, ind4, ind5, ind6, i07i, i07j, i0i, i0j
      integer :: jjorb0, jjorb1, jjorb2, jjorb3, jjorb4, jjorb5, jjorb6
      integer :: jorb0, jorb1, jorb2, jorb3, jorb4, jorb5, jorb6
      real(kind=8) :: tt00, tt01, tt02, tt03, tt04, tt05, tt06
      real(kind=8) :: tt10, tt11, tt12, tt13, tt14, tt15, tt16
      real(kind=8) :: tt20, tt21, tt22, tt23, tt24, tt25, tt26
      real(kind=8) :: tt30, tt31, tt32, tt33, tt34, tt35, tt36
      real(kind=8) :: tt40, tt41, tt42, tt43, tt44, tt45, tt46
      real(kind=8) :: tt50, tt51, tt52, tt53, tt54, tt55, tt56
      real(kind=8) :: tt60, tt61, tt62, tt63, tt64, tt65, tt66
      integer,dimension(:),allocatable :: n
      !$ integer  :: omp_get_thread_num,omp_get_max_threads
      integer(kind=8) :: totops
      integer :: avops, ops, opsn
      integer, allocatable, dimension(:) :: numops
      integer,dimension(:),pointer :: moduloarray
      logical :: ifnd, jfnd
      integer :: iorb, jorb, imat, iseg, iorb_shift, itg, iitg, ist_send, ist_recv, ncount, ishift
      real(kind=8) :: res
      integer,dimension(:),allocatable :: request
      real(kind=8),dimension(:),allocatable :: recvbuf
      integer,dimension(2) :: irowcol
      integer,parameter :: GLOBAL_MATRIX=101, SUBMATRIX=102
      integer,parameter :: data_strategy=SUBMATRIX!GLOBAL_MATRIX
    
      call timing(iproc,'ovrlptransComp','ON')
    
      call f_routine(id='calculate_overlap_transposed')

      call get_modulo_array(smat, moduloarray)
    
      call f_zero(smat%nvctrp_tg*smat%nspin, ovrlp%matrix_compr(1))
    
      ! WARNING: METHOD 2 NOT EXTENSIVELY TESTED
      method_if: if (collcom%imethod_overlap==2) then
    
          stop 'collcom%imethod_overlap=2 is deprecated'
    
       !!   !!!iicnt=0
       !!   !!!sm_it = iterator(collcom)
       !!   !!!do while(valid(sm_it))
       !!   !!!icnt=icnt+1
       !!   !!!call get_position(sm_it,shift=ind0(:),iorb=i0i,jorb=i0j)
       !!   !!!call ge_orbitals(sm_it,iiorb,ijorb)
       !!   !!!call get_ind0(iiorb+(norb)*jjorb,ind0)
       !!   !!!ovrlp%matrix_compr(ind0) = ovrlp%matrix_compr(ind0) + psit_c1(i0i)*psit_c2(i0j)
       !!   !!!call next(sm_it)
       !!   !!!end do
       !!   iorb_shift=(ispin-1)*smat%nfvctr    
       !!   do ispin=1,smat%nspin
       !!     do iseg=1,smat%nseg
       !!       imat=smat%keyv(iseg)
       !!       ! A segment is always on one line, therefore no double loop
       !!       do j=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
       !!         !call get_orbs(smat,i,iorb,jorb) !lookup on work array of size smat%nvctr 
       !!         !iorb=smat%orb_from_index(2,imat)
       !!         !jorb=smat%orb_from_index(1,imat)
       !!         irowcol = orb_from_index(smat, j)
       !!         ovrlp%matrix_compr(imat)=0.0_wp
       !!   
       !!         do ipt=1,collcom%nptsp_c
       !!           ii=collcom%norb_per_gridpoint_c(ipt)
       !!           i0 = collcom%isptsp_c(ipt) + (ispin-1)*collcom%ndimind_c/smat%nspin
       !!           ifnd=.false.
       !!           jfnd=.false.
       !!           do i=1,ii
       !!             iiorb=collcom%indexrecvorbital_c(i0+i) - iorb_shift
       !!             !iiorb=mod(iiorb-1,smat%nfvctr)+1
       !!             if (iiorb == irowcol(1)) then        
       !!                ifnd=.true.
       !!                i0i=i0+i
       !!                !i0i=collcom%iextract_c(i0+i)
       !!             end if 
       !!             if (iiorb == irowcol(2)) then
       !!                 jfnd=.true.
       !!                 i0j=i0+i
       !!                 !i0j=collcom%iextract_c(i0+i)
       !!             end if
       !!             if (.not. (jfnd .and. ifnd)) cycle
       !!             ovrlp%matrix_compr(imat) = ovrlp%matrix_compr(imat) +  psit_c1(i0i)*psit_c2(i0j)
       !!             if (jfnd .and. ifnd) exit
       !!           end do
       !!         end do
       !!   
       !!         do ipt=1,collcom%nptsp_f
       !!           ii=collcom%norb_per_gridpoint_f(ipt)
       !!           i0 = collcom%isptsp_f(ipt) + (ispin-1)*collcom%ndimind_f/smat%nspin
       !!           ifnd=.false.
       !!           jfnd=.false.
       !!           do i=1,ii
       !!             iiorb=collcom%indexrecvorbital_f(i0+i) - iorb_shift
       !!             !iiorb=mod(iiorb-1,smat%nfvctr)+1
       !!             if (iiorb == irowcol(1)) then        
       !!                ifnd=.true.
       !!                i0i=i0+i
       !!                !i0i=collcom%iextract_f(i0+i)
       !!             end if 
       !!             if (iiorb == irowcol(2)) then
       !!                 jfnd=.true.
       !!                 i0j=i0+i
       !!                 !i0j=collcom%iextract_f(i0+i)
       !!             end if
       !!             if (.not. (jfnd .and. ifnd)) cycle
       !!             i07i=7*i0i
       !!             i07j=7*i0j
       !!             ovrlp%matrix_compr(imat) = ovrlp%matrix_compr(imat) +  psit_f1(i07i-6)*psit_f2(i07j-6)
       !!             ovrlp%matrix_compr(imat) = ovrlp%matrix_compr(imat) +  psit_f1(i07i-5)*psit_f2(i07j-5)
       !!             ovrlp%matrix_compr(imat) = ovrlp%matrix_compr(imat) +  psit_f1(i07i-4)*psit_f2(i07j-4)
       !!             ovrlp%matrix_compr(imat) = ovrlp%matrix_compr(imat) +  psit_f1(i07i-3)*psit_f2(i07j-3)
       !!             ovrlp%matrix_compr(imat) = ovrlp%matrix_compr(imat) +  psit_f1(i07i-2)*psit_f2(i07j-2)
       !!             ovrlp%matrix_compr(imat) = ovrlp%matrix_compr(imat) +  psit_f1(i07i-1)*psit_f2(i07j-1)
       !!             ovrlp%matrix_compr(imat) = ovrlp%matrix_compr(imat) +  psit_f1(i07i-0)*psit_f2(i07j-0)
       !!             if (jfnd .and. ifnd) exit
       !!           end do
       !!         end do
       !!         imat=imat+1
       !!        end do 
       !!      end do    
       !!   end do
    
      else if (collcom%imethod_overlap==1) then method_if
    
          !only optimized for spin=1 for now
          ispin=1
          nthreads=1
          !$ nthreads = OMP_GET_max_threads()
          n = f_malloc(nthreads,id='n')
          iorb_shift=(ispin-1)*smat%nfvctr
          ! calculate number of operations for better load balancing of OpenMP
          if (nthreads>1) then
             numops = f_malloc(orbs%norb,id='numops')
             !coarse
             numops=0
             do ipt=1,collcom%nptsp_c
                ii=collcom%norb_per_gridpoint_c(ipt)
                i0 = collcom%isptsp_c(ipt) + (ispin-1)*collcom%ndimind_c/smat%nspin
                do i=1,ii
                   i0i=i0+i
                   iiorb=collcom%indexrecvorbital_c(i0+i) - iorb_shift
                   numops(iiorb)=numops(iiorb)+ii
                end do
             end do
             totops=sum(int(numops,kind=8))
             avops=nint(dble(totops)/dble(nthreads))
             jjorb=1
             do i=1,nthreads
                res=dble(nthreads-i)
                ops=0
                do j=jjorb,orbs%norb
                   opsn=ops+numops(j)
                   if (opsn>=avops) then
                      if ((opsn-avops)<(avops-ops)) then
                         n(i)=j
                         jjorb=j+1
                         totops=totops-int(opsn,kind=8)
                      else
                         n(i)=j-1
                         jjorb=j
                         totops=totops-int(ops,kind=8)
                      end if
                      exit
                   end if
                   ops=opsn
                end do
                if (res /=0.d0) avops=nint(dble(totops)/res)
             end do
             call f_free(numops)
          end if 
    
          n(nthreads)=orbs%norb
        
    
          !$omp parallel default(none) &
          !$omp shared(collcom, smat, ovrlp, psit_c1, psit_c2, psit_f1, psit_f2, n, moduloarray) &
          !$omp private(tid, ispin, iend, istart, ipt, ii, i0, i, iiorb, m, j, i0j, jjorb, ishift_mat, iorb_shift, ind0) &
          !$omp private(jjorb0, jjorb1, ind1, jjorb2, ind2, jjorb3, ind3, jjorb4, ind4, jjorb5, ind5, jjorb6, ind6) &
          !$omp private(i0i, i07i, i07j, tt06, tt05, tt04, tt03, tt02, tt01, tt00) &
          !$omp private(tt16, tt15, tt14, tt13, tt12, tt11, tt10) & 
          !$omp private(tt26, tt25, tt24, tt23, tt22, tt21, tt20) &
          !$omp private(tt36, tt35, tt34, tt33, tt32, tt31, tt30) &
          !$omp private(tt46, tt45, tt44, tt43, tt42, tt41, tt40) &
          !$omp private(tt56, tt55, tt54, tt53, tt52, tt51, tt50) &
          !$omp private(tt66, tt65, tt64, tt63, tt62, tt61, tt60) &
          !$omp private(iorb, jorb, jorb0, jorb1, jorb2, jorb3, jorb4, jorb5, jorb6)
          tid=0
          !$ tid = OMP_GET_THREAD_NUM()
          iend=n(tid+1)
          if (tid==0) then
             istart=1
          else
             istart=n(tid)+1
          end if
        
    
          !SM: check if the modulo operations take a lot of time. If so, try to use an
          !auxiliary array with shifted bounds in order to access smat%matrixindex_in_compressed_fortransposed
          spin_loop: do ispin=1,smat%nspin
        
              ishift_mat=(ispin-1)*smat%nvctrp_tg-smat%isvctrp_tg
              iorb_shift=(ispin-1)*smat%nfvctr
              if (collcom%nptsp_c>0) then
        
                  do ipt=1,collcom%nptsp_c 
                      ii=collcom%norb_per_gridpoint_c(ipt) 
                      i0 = collcom%isptsp_c(ipt) + (ispin-1)*collcom%ndimind_c/smat%nspin
                      do i=1,ii
                          i0i=i0+i
                          iiorb=collcom%indexrecvorbital_c(i0i) - iorb_shift
                          iorb=moduloarray(iiorb)
                          !iiorb=mod(iiorb-1,smat%nfvctr)+1
                          if(iiorb < istart .or. iiorb > iend) cycle
                          m=mod(ii,7)
                          if(m/=0) then
                              do j=1,m
                                  i0j=i0+j
                                  jjorb=collcom%indexrecvorbital_c(i0j) - iorb_shift
                                  jorb=moduloarray(jjorb)
                                  !jjorb=mod(jjorb-1,smat%nfvctr)+1
                                  !ind0 = ishift_mat + smat%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
                                  !ind0 = smat%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
                                  ind0 = smat%matrixindex_in_compressed_fortransposed(jorb,iorb)
                                  !ind0 = get_transposed_index(smat,jjorb,iiorb)
                                  ind0=ind0+ishift_mat
                                  !if (ind0>=smat%nvctr-smat%nfvctr .and.  ind0<=smat%nvctr) then
                                  !    write(*,'(a,3i9)') 'iiorb, jjorb, ind0', iiorb, jjorb, ind0
                                  !end if
                                  !!write(880,'(a,5i8,es14.6)') 'ispin, ipt, i, ind0, i0j, val', ispin, ipt, i, ind0, i0j, psit_c1(i0i)
                                  ovrlp%matrix_compr(ind0) = ovrlp%matrix_compr(ind0) + psit_c1(i0i)*psit_c2(i0j)
                              end do
                          end if
                          do j=m+1,ii,7
                              i0j=i0+j
        
                              jjorb0=collcom%indexrecvorbital_c(i0j+0) - iorb_shift
                              jorb0=moduloarray(jjorb0)
                              !jjorb0=mod(jjorb0-1,smat%nfvctr)+1
                              !ind0 = ishift_mat + smat%matrixindex_in_compressed_fortransposed(jjorb0,iiorb)
                              !ind0 = smat%matrixindex_in_compressed_fortransposed(jjorb0,iiorb)
                              ind0 = smat%matrixindex_in_compressed_fortransposed(jorb0,iorb)
                              !ind0 = get_transposed_index(smat,jjorb0,iiorb)
                              ind0=ind0+ishift_mat
                              ovrlp%matrix_compr(ind0) = ovrlp%matrix_compr(ind0) + psit_c1(i0i)*psit_c2(i0j+0)
        
                              jjorb1=collcom%indexrecvorbital_c(i0j+1) - iorb_shift
                              jorb1=moduloarray(jjorb1)
                              !jjorb1=mod(jjorb1-1,smat%nfvctr)+1
                              !ind1 = ishift_mat + smat%matrixindex_in_compressed_fortransposed(jjorb1,iiorb)
                              !ind1 = smat%matrixindex_in_compressed_fortransposed(jjorb1,iiorb)
                              ind1 = smat%matrixindex_in_compressed_fortransposed(jorb1,iorb)
                              !ind1 = get_transposed_index(smat,jjorb1,iiorb)
                              ind1=ind1+ishift_mat
                              ovrlp%matrix_compr(ind1) = ovrlp%matrix_compr(ind1) + psit_c1(i0i)*psit_c2(i0j+1)
        
                              jjorb2=collcom%indexrecvorbital_c(i0j+2) - iorb_shift
                              jorb2=moduloarray(jjorb2)
                              !jjorb2=mod(jjorb2-1,smat%nfvctr)+1
                              !ind2 = ishift_mat + smat%matrixindex_in_compressed_fortransposed(jjorb2,iiorb)
                              !ind2 = smat%matrixindex_in_compressed_fortransposed(jjorb2,iiorb)
                              ind2 = smat%matrixindex_in_compressed_fortransposed(jorb2,iorb)
                              !ind2 = get_transposed_index(smat,jjorb2,iiorb)
                              ind2=ind2+ishift_mat
                              ovrlp%matrix_compr(ind2) = ovrlp%matrix_compr(ind2) + psit_c1(i0i)*psit_c2(i0j+2)
        
                              jjorb3=collcom%indexrecvorbital_c(i0j+3) - iorb_shift
                              jorb3=moduloarray(jjorb3)
                              !jjorb3=mod(jjorb3-1,smat%nfvctr)+1
                              !ind3 = ishift_mat + smat%matrixindex_in_compressed_fortransposed(jjorb3,iiorb)
                              !ind3 = smat%matrixindex_in_compressed_fortransposed(jjorb3,iiorb)
                              ind3 = smat%matrixindex_in_compressed_fortransposed(jorb3,iorb)
                              !ind3 = get_transposed_index(smat,jjorb3,iiorb)
                              ind3=ind3+ishift_mat
                              ovrlp%matrix_compr(ind3) = ovrlp%matrix_compr(ind3) + psit_c1(i0i)*psit_c2(i0j+3)
        
                              jjorb4=collcom%indexrecvorbital_c(i0j+4) - iorb_shift
                              jorb4=moduloarray(jjorb4)
                              !jjorb4=mod(jjorb4-1,smat%nfvctr)+1
                              !ind4 = ishift_mat + smat%matrixindex_in_compressed_fortransposed(jjorb4,iiorb)
                              !ind4 = smat%matrixindex_in_compressed_fortransposed(jjorb4,iiorb)
                              ind4 = smat%matrixindex_in_compressed_fortransposed(jorb4,iorb)
                              !ind4 = get_transposed_index(smat,jjorb4,iiorb)
                              ind4=ind4+ishift_mat
                              ovrlp%matrix_compr(ind4) = ovrlp%matrix_compr(ind4) + psit_c1(i0i)*psit_c2(i0j+4)
        
                              jjorb5=collcom%indexrecvorbital_c(i0j+5) - iorb_shift
                              jorb5=moduloarray(jjorb5)
                              !jjorb5=mod(jjorb5-1,smat%nfvctr)+1
                              !ind5 = ishift_mat + smat%matrixindex_in_compressed_fortransposed(jjorb5,iiorb)
                              !ind5 = smat%matrixindex_in_compressed_fortransposed(jjorb5,iiorb)
                              ind5 = smat%matrixindex_in_compressed_fortransposed(jorb5,iorb)
                              !ind5 = get_transposed_index(smat,jjorb5,iiorb)
                              ind5=ind5+ishift_mat
                              ovrlp%matrix_compr(ind5) = ovrlp%matrix_compr(ind5) + psit_c1(i0i)*psit_c2(i0j+5)
        
                              jjorb6=collcom%indexrecvorbital_c(i0j+6) - iorb_shift
                              jorb6=moduloarray(jjorb6)
                              !jjorb6=mod(jjorb6-1,smat%nfvctr)+1
                              !ind6 = ishift_mat + smat%matrixindex_in_compressed_fortransposed(jjorb6,iiorb)
                              !ind6 = smat%matrixindex_in_compressed_fortransposed(jjorb6,iiorb)
                              ind6 = smat%matrixindex_in_compressed_fortransposed(jorb6,iorb)
                              !ind6 = get_transposed_index(smat,jjorb6,iiorb)
                              ind6=ind6+ishift_mat
                              ovrlp%matrix_compr(ind6) = ovrlp%matrix_compr(ind6) + psit_c1(i0i)*psit_c2(i0j+6)
        
                          end do
                      end do
                  end do
              end if
          end do spin_loop
          !$omp end parallel
    
    
          !recalculate best OpenMP load balancing for fine grid points - not necessarily the same as coarse
          !only optimized for spin=1 for now
          ispin=1
          nthreads=1
          !$  nthreads = OMP_GET_max_threads()
          iorb_shift=(ispin-1)*smat%nfvctr
          ! calculate number of operations for better load balancing of OpenMP
          if (nthreads>1) then
             numops = f_malloc(orbs%norb,id='numops')
             numops=0
             !fine
             do ipt=1,collcom%nptsp_f
                ii=collcom%norb_per_gridpoint_f(ipt)
                i0 = collcom%isptsp_f(ipt) + (ispin-1)*collcom%ndimind_f/smat%nspin
                do i=1,ii
                   i0i=i0+i
                   iiorb=collcom%indexrecvorbital_f(i0i) - iorb_shift
                   numops(iiorb)=numops(iiorb)+ii  !*7
                end do
             end do
             totops=sum(int(numops,kind=8))
             avops=nint(totops/dble(nthreads))
             jjorb=1
             do i=1,nthreads
                res=dble(nthreads-i)
                ops=0
                do j=jjorb,orbs%norb
                   opsn=ops+numops(j)
                   if (opsn>=avops) then
                      if ((opsn-avops)<(avops-ops)) then
                         n(i)=j
                         jjorb=j+1
                         totops=totops-int(opsn,kind=8)
                      else
                         n(i)=j-1
                         jjorb=j
                         totops=totops-int(ops,kind=8)
                      end if
                      exit
                   end if
                   ops=opsn
                end do
                if (res /= 0.d0) avops=nint(dble(totops)/res)
             end do
             call f_free(numops)
          end if    
    
          n(nthreads)=orbs%norb
        
    
          !$omp parallel default(none) &
          !$omp shared(collcom, smat, ovrlp, psit_c1, psit_c2, psit_f1, psit_f2, n, moduloarray) &
          !$omp private(tid, ispin, iend, istart, ipt, ii, i0, i, iiorb, m, j, i0j, jjorb, ishift_mat, iorb_shift, ind0) &
          !$omp private(jjorb0, jjorb1, ind1, jjorb2, ind2, jjorb3, ind3, jjorb4, ind4, jjorb5, ind5, jjorb6, ind6) &
          !$omp private(i0i, i07i, i07j, tt06, tt05, tt04, tt03, tt02, tt01, tt00) &
          !$omp private(tt16, tt15, tt14, tt13, tt12, tt11, tt10) & 
          !$omp private(tt26, tt25, tt24, tt23, tt22, tt21, tt20) &
          !$omp private(tt36, tt35, tt34, tt33, tt32, tt31, tt30) &
          !$omp private(tt46, tt45, tt44, tt43, tt42, tt41, tt40) &
          !$omp private(tt56, tt55, tt54, tt53, tt52, tt51, tt50) &
          !$omp private(tt66, tt65, tt64, tt63, tt62, tt61, tt60) &
          !$omp private(iorb, jorb, jorb0, jorb1, jorb2, jorb3, jorb4, jorb5, jorb6)
          tid=0
          !$ tid = OMP_GET_THREAD_NUM()
          iend=n(tid+1)
          if (tid==0) then
             istart=1
          else
             istart=n(tid)+1
          end if
        
          spin_loopf: do ispin=1,smat%nspin
        
              ishift_mat=(ispin-1)*smat%nvctr-smat%isvctrp_tg
              iorb_shift=(ispin-1)*smat%nfvctr
    
              if (collcom%nptsp_f>0) then
                  do ipt=1,collcom%nptsp_f 
                      ii=collcom%norb_per_gridpoint_f(ipt) 
                      i0 = collcom%isptsp_f(ipt) + (ispin-1)*collcom%ndimind_f/smat%nspin
                      do i=1,ii
                          i0i=i0+i
                          iiorb=collcom%indexrecvorbital_f(i0i) - iorb_shift
                          iorb=moduloarray(iiorb)
                          !iiorb=mod(iiorb-1,smat%nfvctr)+1
                          if(iiorb < istart .or. iiorb > iend) cycle
                          i07i=7*i0i
                          m=mod(ii,7)
                          if(m/=0) then
                              do j=1,m
                                  i0j=i0+j
                                  i07j=7*i0j
                                  jjorb0=collcom%indexrecvorbital_f(i0j) - iorb_shift
                                  jorb0=moduloarray(jjorb0)
                                  !jjorb0=mod(jjorb0-1,smat%nfvctr)+1
                                  !ind0 = ishift_mat + smat%matrixindex_in_compressed_fortransposed(jjorb0,iiorb)
                                  !ind0 = smat%matrixindex_in_compressed_fortransposed(jjorb0,iiorb)
                                  ind0 = smat%matrixindex_in_compressed_fortransposed(jorb0,iorb)
                                  !ind0 = get_transposed_index(smat,jjorb0,iiorb)
                                  ind0=ind0+ishift_mat
                                  tt06 = psit_f1(i07i-6)*psit_f2(i07j-6)
                                  tt05 = psit_f1(i07i-5)*psit_f2(i07j-5)
                                  tt04 = psit_f1(i07i-4)*psit_f2(i07j-4)
                                  tt03 = psit_f1(i07i-3)*psit_f2(i07j-3)
                                  tt02 = psit_f1(i07i-2)*psit_f2(i07j-2)
                                  tt01 = psit_f1(i07i-1)*psit_f2(i07j-1)
                                  tt00 = psit_f1(i07i-0)*psit_f2(i07j-0)
        
                                  ovrlp%matrix_compr(ind0) = ovrlp%matrix_compr(ind0) &
                                      + tt06 + tt05 + tt04 + tt03 + tt02 + tt01 + tt00
                              end do
                          end if
                          do j=m+1,ii,7
                              i0j=i0+j
                              i07j=7*i0j
                              jjorb0=collcom%indexrecvorbital_f(i0j+0) - iorb_shift
                              jorb0=moduloarray(jjorb0)
                              !jjorb0=mod(jjorb0-1,smat%nfvctr)+1
                              !ind0 = ishift_mat + smat%matrixindex_in_compressed_fortransposed(jjorb0,iiorb)
                              !ind0 = smat%matrixindex_in_compressed_fortransposed(jjorb0,iiorb)
                              ind0 = smat%matrixindex_in_compressed_fortransposed(jorb0,iorb)
                              !ind0 = get_transposed_index(smat,jjorb0,iiorb)
                              ind0=ind0+ishift_mat
                              tt06 = psit_f1(i07i-6)*psit_f2(i07j-6)
                              tt05 = psit_f1(i07i-5)*psit_f2(i07j-5)
                              tt04 = psit_f1(i07i-4)*psit_f2(i07j-4)
                              tt03 = psit_f1(i07i-3)*psit_f2(i07j-3)
                              tt02 = psit_f1(i07i-2)*psit_f2(i07j-2)
                              tt01 = psit_f1(i07i-1)*psit_f2(i07j-1)
                              tt00 = psit_f1(i07i-0)*psit_f2(i07j-0)
                              ovrlp%matrix_compr(ind0) = ovrlp%matrix_compr(ind0) + tt06 + tt05 + tt04 + tt03 + tt02 + tt01 + tt00
        
                              jjorb1=collcom%indexrecvorbital_f(i0j+1) - iorb_shift
                              jorb1=moduloarray(jjorb1)
                              !jjorb1=mod(jjorb1-1,smat%nfvctr)+1
                              !ind1 = ishift_mat + smat%matrixindex_in_compressed_fortransposed(jjorb1,iiorb)
                              !ind1 = smat%matrixindex_in_compressed_fortransposed(jjorb1,iiorb)
                              ind1 = smat%matrixindex_in_compressed_fortransposed(jorb1,iorb)
                              !ind1 = get_transposed_index(smat,jjorb1,iiorb)
                              ind1=ind1+ishift_mat
                              tt16 = psit_f1(i07i-6)*psit_f2(i07j+1) !+1*7-6
                              tt15 = psit_f1(i07i-5)*psit_f2(i07j+2) !+1*7-5
                              tt14 = psit_f1(i07i-4)*psit_f2(i07j+3) !+1*7-4
                              tt13 = psit_f1(i07i-3)*psit_f2(i07j+4) !+1*7-3
                              tt12 = psit_f1(i07i-2)*psit_f2(i07j+5) !+1*7-2
                              tt11 = psit_f1(i07i-1)*psit_f2(i07j+6) !+1*7-1
                              tt10 = psit_f1(i07i-0)*psit_f2(i07j+7) !+1*7-0
                              ovrlp%matrix_compr(ind1) = ovrlp%matrix_compr(ind1) + tt16 + tt15 + tt14 + tt13 + tt12 + tt11 + tt10
        
                              jjorb2=collcom%indexrecvorbital_f(i0j+2) - iorb_shift
                              jorb2=moduloarray(jjorb2)
                              !jjorb2=mod(jjorb2-1,smat%nfvctr)+1
                              !ind2 = ishift_mat + smat%matrixindex_in_compressed_fortransposed(jjorb2,iiorb)
                              !ind2 = smat%matrixindex_in_compressed_fortransposed(jjorb2,iiorb)
                              ind2 = smat%matrixindex_in_compressed_fortransposed(jorb2,iorb)
                              !ind2 = get_transposed_index(smat,jjorb2,iiorb)
                              ind2=ind2+ishift_mat
                              tt26 = psit_f1(i07i-6)*psit_f2(i07j+8) !+2*7-6
                              tt25 = psit_f1(i07i-5)*psit_f2(i07j+9) !+2*7-5
                              tt24 = psit_f1(i07i-4)*psit_f2(i07j+10) !+2*7-4
                              tt23 = psit_f1(i07i-3)*psit_f2(i07j+11) !+2*7-3
                              tt22 = psit_f1(i07i-2)*psit_f2(i07j+12) !+2*7-2
                              tt21 = psit_f1(i07i-1)*psit_f2(i07j+13) !+2*7-1
                              tt20 = psit_f1(i07i-0)*psit_f2(i07j+14) !+2*7-0
                              ovrlp%matrix_compr(ind2) = ovrlp%matrix_compr(ind2) + tt26 + tt25 + tt24 + tt23 + tt22 + tt21 + tt20
        
                              jjorb3=collcom%indexrecvorbital_f(i0j+3) - iorb_shift
                              jorb3=moduloarray(jjorb3)
                              !jjorb3=mod(jjorb3-1,smat%nfvctr)+1
                              !ind3 = ishift_mat + smat%matrixindex_in_compressed_fortransposed(jjorb3,iiorb)
                              !ind3 = smat%matrixindex_in_compressed_fortransposed(jjorb3,iiorb)
                              ind3 = smat%matrixindex_in_compressed_fortransposed(jorb3,iorb)
                              !ind3 = get_transposed_index(smat,jjorb3,iiorb)
                              ind3=ind3+ishift_mat
                              tt36 = psit_f1(i07i-6)*psit_f2(i07j+15) !+3*7-6
                              tt35 = psit_f1(i07i-5)*psit_f2(i07j+16) !+3*7-5
                              tt34 = psit_f1(i07i-4)*psit_f2(i07j+17) !+3*7-4
                              tt33 = psit_f1(i07i-3)*psit_f2(i07j+18) !+3*7-3
                              tt32 = psit_f1(i07i-2)*psit_f2(i07j+19) !+3*7-2
                              tt31 = psit_f1(i07i-1)*psit_f2(i07j+20) !+3*7-1
                              tt30 = psit_f1(i07i-0)*psit_f2(i07j+21) !+3*7-0
                              ovrlp%matrix_compr(ind3) = ovrlp%matrix_compr(ind3) + tt36 + tt35 + tt34 + tt33 + tt32 + tt31 + tt30
        
                              jjorb4=collcom%indexrecvorbital_f(i0j+4) - iorb_shift
                              jorb4=moduloarray(jjorb4)
                              !jjorb4=mod(jjorb4-1,smat%nfvctr)+1
                              !ind4 = ishift_mat + smat%matrixindex_in_compressed_fortransposed(jjorb4,iiorb)
                              !ind4 = smat%matrixindex_in_compressed_fortransposed(jjorb4,iiorb)
                              ind4 = smat%matrixindex_in_compressed_fortransposed(jorb4,iorb)
                              !ind4 = get_transposed_index(smat,jjorb4,iiorb)
                              ind4=ind4+ishift_mat
                              tt46 = psit_f1(i07i-6)*psit_f2(i07j+22) !+4*7-6
                              tt45 = psit_f1(i07i-5)*psit_f2(i07j+23) !+4*7-5
                              tt44 = psit_f1(i07i-4)*psit_f2(i07j+24) !+4*7-4
                              tt43 = psit_f1(i07i-3)*psit_f2(i07j+25) !+4*7-3
                              tt42 = psit_f1(i07i-2)*psit_f2(i07j+26) !+4*7-2
                              tt41 = psit_f1(i07i-1)*psit_f2(i07j+27) !+4*7-1
                              tt40 = psit_f1(i07i-0)*psit_f2(i07j+28) !+4*7-0
                              ovrlp%matrix_compr(ind4) = ovrlp%matrix_compr(ind4) + tt46 + tt45 + tt44 + tt43 + tt42 + tt41 + tt40
        
                              jjorb5=collcom%indexrecvorbital_f(i0j+5) - iorb_shift
                              jorb5=moduloarray(jjorb5)
                              !jjorb5=mod(jjorb5-1,smat%nfvctr)+1
                              !ind5 = ishift_mat + smat%matrixindex_in_compressed_fortransposed(jjorb5,iiorb)
                              !ind5 = smat%matrixindex_in_compressed_fortransposed(jjorb5,iiorb)
                              ind5 = smat%matrixindex_in_compressed_fortransposed(jorb5,iorb)
                              !ind5 = get_transposed_index(smat,jjorb5,iiorb)
                              ind5=ind5+ishift_mat
                              tt56 = psit_f1(i07i-6)*psit_f2(i07j+29) !+5*7-6
                              tt55 = psit_f1(i07i-5)*psit_f2(i07j+30) !+5*7-5
                              tt54 = psit_f1(i07i-4)*psit_f2(i07j+31) !+5*7-4
                              tt53 = psit_f1(i07i-3)*psit_f2(i07j+32) !+5*7-3
                              tt52 = psit_f1(i07i-2)*psit_f2(i07j+33) !+5*7-2
                              tt51 = psit_f1(i07i-1)*psit_f2(i07j+34) !+5*7-1
                              tt50 = psit_f1(i07i-0)*psit_f2(i07j+35) !+5*7-0
                              ovrlp%matrix_compr(ind5) = ovrlp%matrix_compr(ind5) + tt56 + tt55 + tt54 + tt53 + tt52 + tt51 + tt50
        
                              jjorb6=collcom%indexrecvorbital_f(i0j+6) - iorb_shift
                              jorb6=moduloarray(jjorb6)
                              !jjorb6=mod(jjorb6-1,smat%nfvctr)+1
                              !ind6 = ishift_mat + smat%matrixindex_in_compressed_fortransposed(jjorb6,iiorb)
                              !ind6 = smat%matrixindex_in_compressed_fortransposed(jjorb6,iiorb)
                              ind6 = smat%matrixindex_in_compressed_fortransposed(jorb6,iorb)
                              !ind6 = get_transposed_index(smat,jjorb6,iiorb)
                              ind6=ind6+ishift_mat
                              tt66 = psit_f1(i07i-6)*psit_f2(i07j+36) !+6*7-6
                              tt65 = psit_f1(i07i-5)*psit_f2(i07j+37) !+6*7-5
                              tt64 = psit_f1(i07i-4)*psit_f2(i07j+38) !+6*7-4
                              tt63 = psit_f1(i07i-3)*psit_f2(i07j+39) !+6*7-3
                              tt62 = psit_f1(i07i-2)*psit_f2(i07j+40) !+6*7-2
                              tt61 = psit_f1(i07i-1)*psit_f2(i07j+41) !+6*7-1
                              tt60 = psit_f1(i07i-0)*psit_f2(i07j+42) !+6*7-0
                              ovrlp%matrix_compr(ind6) = ovrlp%matrix_compr(ind6) + tt66 + tt65 + tt64 + tt63 + tt62 + tt61 + tt60
                          end do
                      end do
                  end do
              end if
        
          end do spin_loopf
          !$omp end parallel
    
          call f_free(n)
          call f_free_ptr(moduloarray)
    
      else method_if
          stop 'wrong value of imethod_if'
      end if method_if
    
      call timing(iproc,'ovrlptransComp','OF')
    
      call timing(iproc,'ovrlptransComm','ON')

      if (data_strategy==GLOBAL_MATRIX) then
          if(nproc > 1) then
              call mpiallred(ovrlp%matrix_compr(1), smat%nvctr*smat%nspin, mpi_sum,comm=bigdft_mpi%mpi_comm)
          end if
    
      else if (data_strategy==SUBMATRIX) then
    
          call synchronize_matrix_taskgroups(iproc, nproc, smat, ovrlp)
    
          !!if (nproc>1) then
          !!    request = f_malloc(smat%ntaskgroupp,id='request')
          !!    ncount = 0
          !!    do itg=1,smat%ntaskgroupp
          !!        iitg = smat%taskgroupid(itg)
          !!        ncount = ncount + smat%taskgroup_startend(2,1,iitg)-smat%taskgroup_startend(1,1,iitg)+1
          !!    end do
          !!    recvbuf = f_malloc(ncount,id='recvbuf')
          !!    do ispin=1,smat%nspin
          !!        ishift = (ispin-1)*smat%nvctr
    
          !!        ncount = 0
          !!        do itg=1,smat%ntaskgroupp
          !!            iitg = smat%taskgroupid(itg)
          !!            ist_send = smat%taskgroup_startend(1,1,iitg)
          !!            ist_recv = ncount + 1
          !!            ncount = smat%taskgroup_startend(2,1,iitg)-smat%taskgroup_startend(1,1,iitg)+1
          !!            !!call mpi_iallreduce(ovrlp%matrix_compr(ist_send), recvbuf(ist_recv), ncount, &
          !!            !!     mpi_double_precision, mpi_sum, smat%mpi_groups(iitg)%mpi_comm, request(itg), ierr)
          !!            if (nproc>1) then
          !!                call mpiiallred(ovrlp%matrix_compr(ishift+ist_send), recvbuf(ist_recv), ncount, &
          !!                     mpi_sum, smat%mpi_groups(iitg)%mpi_comm, request(itg))
          !!            else
          !!                call vcopy(ncount, ovrlp%matrix_compr(ishift+ist_send), 1, recvbuf(ist_recv), 1)
          !!            end if
          !!        end do
          !!        if (nproc>1) then
          !!            call mpiwaitall(smat%ntaskgroupp, request)
          !!        end if
          !!        ncount = 0
          !!        do itg=1,smat%ntaskgroupp
          !!            iitg = smat%taskgroupid(itg)
          !!            ist_send = smat%taskgroup_startend(1,1,iitg)
          !!            ist_recv = ncount + 1
          !!            ncount = smat%taskgroup_startend(2,1,iitg)-smat%taskgroup_startend(1,1,iitg)+1
          !!            call vcopy(ncount, recvbuf(ist_recv), 1, ovrlp%matrix_compr(ishift+ist_send), 1)
          !!        end do
          !!    end do
          !!    call f_free(request)
          !!    call f_free(recvbuf)
          !!end if
      else
          stop 'calculate_overlap_transposed: wrong data_strategy'
      end if
    
      ! Indicate the matrix is the "original one" (and not the inverse etc.)
      ovrlp%power = 0.d0
    
    
      smat%can_use_dense=.false.
    
      call f_release_routine()
      call timing(iproc,'ovrlptransComm','OF')

      !contains

        !function get_transposed_index(jorb,iorb) res(ind)
        !    integer,intent(in) :: jorb, iorb
        !    integer :: ind
        !    integer :: jjorb,iiorb
        !    ! If iorb is smaller than the offset, add a periodic shift
        !    if (iorb<smat%offset_matrixindex_in_compressed_fortransposed) then
        !        iiorb = iorb + smat%nfvctr
        !    else
        !        iiorb = iorb
        !    end if
        !    if (jorb<smat%offset_matrixindex_in_compressed_fortransposed) then
        !        jjorb = jorb + smat%nfvctr
        !    else
        !        jjorb = jorb
        !    end if
        !    ind = smat%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
        !end function get_transposed_index
    
    end subroutine calculate_overlap_transposed
    
    
    
    subroutine normalize_transposed(iproc, nproc, orbs, nspin, collcom, psit_c, psit_f, norm)
      use module_base
      use module_types, only: orbitals_data
      use communications_base, only: comms_linear
      implicit none
      
      ! Calling arguments
      integer,intent(in):: iproc, nproc, nspin
      type(orbitals_data),intent(in):: orbs
      type(comms_linear),intent(in):: collcom
      real(8),dimension(collcom%ndimind_c),intent(inout):: psit_c
      real(8),dimension(7*collcom%ndimind_f),intent(inout):: psit_f
      real(8),dimension(orbs%norb),intent(out):: norm
      
      ! Local variables
      integer:: i0, ipt, ii, iiorb, i, ierr, iorb, i07i, i0i, ispin
    
      call timing(iproc,'norm_trans','ON')
    
      call f_zero(norm)
    
    
      spin_loop: do ispin=1,nspin
    
          !$omp parallel default(private) &
          !$omp shared(collcom, norm, psit_c,psit_f,orbs,ispin,nspin)
          if (collcom%nptsp_c>0) then
              !$omp do reduction(+:norm)
              do ipt=1,collcom%nptsp_c 
                  ii=collcom%norb_per_gridpoint_c(ipt)
                  i0 = collcom%isptsp_c(ipt) + (ispin-1)*collcom%ndimind_c/nspin
                  do i=1,ii
                      i0i=i0+i
                      iiorb=collcom%indexrecvorbital_c(i0i)
                      !!write(720,'(a,6i8,es13.5)') 'ipt, ispin, i0, i, i0i, iiorb, psit_c(i0i)', &
                      !!    ipt, ispin, i0, i, i0i, iiorb, psit_c(i0i)
                      norm(iiorb)=norm(iiorb)+psit_c(i0i)**2
                  end do
              end do
              !$omp end do
          end if
    
          if (collcom%nptsp_f>0) then
              !$omp do reduction(+:norm)
              do ipt=1,collcom%nptsp_f 
                  ii=collcom%norb_per_gridpoint_f(ipt) 
                  i0 = collcom%isptsp_f(ipt) + (ispin-1)*collcom%ndimind_f/nspin
                  do i=1,ii
                      i0i=i0+i
                      i07i=7*i0i
                      iiorb=collcom%indexrecvorbital_f(i0i)
                      norm(iiorb)=norm(iiorb)+psit_f(i07i-6)**2
                      norm(iiorb)=norm(iiorb)+psit_f(i07i-5)**2
                      norm(iiorb)=norm(iiorb)+psit_f(i07i-4)**2
                      norm(iiorb)=norm(iiorb)+psit_f(i07i-3)**2
                      norm(iiorb)=norm(iiorb)+psit_f(i07i-2)**2
                      norm(iiorb)=norm(iiorb)+psit_f(i07i-1)**2
                      norm(iiorb)=norm(iiorb)+psit_f(i07i-0)**2
                  end do
              end do
              !$omp end do
          end if
          !$omp end parallel
    
      end do spin_loop
      
      if(nproc>1) then
          call mpiallred(norm, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if
    
      do iorb=1,orbs%norb
         norm(iorb)=1.d0/sqrt(norm(iorb))
      end do
    
      spin_loop2: do ispin=1,nspin
    
          !$omp parallel default(private) shared(norm,orbs,collcom,psit_c,psit_f,ispin,nspin)  
          if (collcom%nptsp_c>0) then
              !$omp do
              do ipt=1,collcom%nptsp_c 
                  ii=collcom%norb_per_gridpoint_c(ipt)
                  i0=collcom%isptsp_c(ipt) + (ispin-1)*collcom%ndimind_c/nspin
                  do i=1,ii
                      i0i=i0+i
                      iiorb=collcom%indexrecvorbital_c(i0i)
                      psit_c(i0i)=psit_c(i0i)*norm(iiorb)
                  end do 
              end do
              !$omp end do
          end if
          if (collcom%nptsp_f>0) then
              !$omp do
              do ipt=1,collcom%nptsp_f 
                  ii=collcom%norb_per_gridpoint_f(ipt)
                  i0 = collcom%isptsp_f(ipt) + (ispin-1)*collcom%ndimind_f/nspin
                  do i=1,ii
                      i0i=i0+i
                      i07i=7*i0i
                      iiorb=collcom%indexrecvorbital_f(i0i)
                      psit_f(i07i-6)=psit_f(i07i-6)*norm(iiorb)
                      psit_f(i07i-5)=psit_f(i07i-5)*norm(iiorb)
                      psit_f(i07i-4)=psit_f(i07i-4)*norm(iiorb)
                      psit_f(i07i-3)=psit_f(i07i-3)*norm(iiorb)
                      psit_f(i07i-2)=psit_f(i07i-2)*norm(iiorb)
                      psit_f(i07i-1)=psit_f(i07i-1)*norm(iiorb)
                      psit_f(i07i-0)=psit_f(i07i-0)*norm(iiorb)
                  end do
              end do
              !$omp end do
          end if
          !$omp end parallel
    
      end do spin_loop2
    
      call timing(iproc,'norm_trans','OF')
    
    end subroutine normalize_transposed
    
    
    subroutine build_linear_combination_transposed(collcom, sparsemat, mat, psitwork_c, psitwork_f, &
         reset, psit_c, psit_f, iproc)
      use module_base
      use communications_base, only: comms_linear
      use sparsematrix_base, only: sparse_matrix, matrices
      use sparsematrix_init, only: get_modulo_array
      implicit none
      
      ! Calling arguments
      type(sparse_matrix),intent(in) :: sparsemat
      type(matrices),intent(in) :: mat
      type(comms_linear),intent(in) :: collcom
      real(kind=8),dimension(collcom%ndimind_c),intent(in) :: psitwork_c
      real(kind=8),dimension(7*collcom%ndimind_f),intent(in) :: psitwork_f
      logical,intent(in) :: reset
      real(kind=8),dimension(collcom%ndimind_c),intent(inout) :: psit_c
      real(kind=8),dimension(7*collcom%ndimind_f),intent(inout) :: psit_f
      integer, intent(in) :: iproc
      ! Local variables
      integer :: i0, ipt, ii, j, iiorb, jjorb, i, m, ind0, ind1, ind2, ind3, i0i, i0j, i07i, i07j, iorb_shift
      integer :: jorb0, jorb1, jorb2, jorb3, jorb4, jorb5, jorb6, iorb, jorb
      integer :: ind4, ind5, ind6, jjorb0, jjorb1, jjorb2, jjorb3, jjorb4, jjorb5, jjorb6, ispin, ishift_mat
      integer,dimension(:),pointer :: moduloarray
      real(kind=8) :: tt0, tt1, tt2, tt3, tt4, tt5, tt6
      real(kind=8) :: tt00, tt01, tt02, tt03, tt04, tt05, tt06
      real(kind=8) :: tt10, tt11, tt12, tt13, tt14, tt15, tt16
      real(kind=8) :: tt20, tt21, tt22, tt23, tt24, tt25, tt26
      real(kind=8) :: tt30, tt31, tt32, tt33, tt34, tt35, tt36
      real(kind=8) :: tt40, tt41, tt42, tt43, tt44, tt45, tt46
      real(kind=8) :: tt50, tt51, tt52, tt53, tt54, tt55, tt56
      real(kind=8) :: tt60, tt61, tt62, tt63, tt64, tt65, tt66
    
      call f_routine(id='build_linear_combination_transposed')
      call timing(iproc,'lincombtrans  ','ON')

      call get_modulo_array(sparsemat, moduloarray)

      if(reset) then
         call f_zero(psit_c)
         call f_zero(psit_f)
      end if
    
    
      !SM: check if the modulo operations take a lot of time. If so, try to use an
      !auxiliary array with shifted bounds in order to access smat%matrixindex_in_compressed_fortransposed
    
      spin_loop: do ispin=1,sparsemat%nspin
    
          ishift_mat=(ispin-1)*sparsemat%nvctr-sparsemat%isvctrp_tg
          iorb_shift=(ispin-1)*sparsemat%nfvctr
    
          !$omp parallel default(private) &
          !$omp shared(collcom, psit_c, psitwork_c, psit_f, psitwork_f, sparsemat) &
          !$omp shared(mat, ispin, ishift_mat, iorb_shift, moduloarray)
    
          !$omp do schedule(static,1)
           do ipt=1,collcom%nptsp_c 
              ii=collcom%norb_per_gridpoint_c(ipt) 
              i0 = collcom%isptsp_c(ipt) + (ispin-1)*collcom%ndimind_c/sparsemat%nspin
              do i=1,ii
                  i0i=i0+i
                  iiorb=collcom%indexrecvorbital_c(i0i) - iorb_shift
                  iorb=moduloarray(iiorb)
                  !iiorb=mod(iiorb-1,sparsemat%nfvctr)+1
                  m=mod(ii,7)
                  tt0=0.d0 ; tt1=0.d0 ; tt2=0.d0 ; tt3=0.d0 ; tt4=0.d0 ; tt5=0.d0 ; tt6=0.d0
                  if(m/=0) then
                      do j=1,m
                          i0j=i0+j
                          jjorb=collcom%indexrecvorbital_c(i0j) - iorb_shift
                          jorb=moduloarray(jjorb)
                          !jjorb=mod(jjorb-1,sparsemat%nfvctr)+1
                          !ind0 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
                          ind0 = sparsemat%matrixindex_in_compressed_fortransposed(jorb,iorb)
                          !ind0 = get_transposed_index(sparsemat,jjorb,iiorb)
                          ind0=ind0+ishift_mat
                          tt0=tt0+mat%matrix_compr(ind0)*psitwork_c(i0j)
                      end do
                  end if
                  do j=m+1,ii,7
                      i0j=i0+j
    
                      jjorb0=collcom%indexrecvorbital_c(i0j+0) - iorb_shift
                      jorb0=moduloarray(jjorb0)
                      !jjorb0=mod(jjorb0-1,sparsemat%nfvctr)+1
                      !ind0 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb0,iiorb)
                      ind0 = sparsemat%matrixindex_in_compressed_fortransposed(jorb0,iorb)
                      !ind0 = get_transposed_index(sparsemat,jjorb0,iiorb)
                      ind0=ind0+ishift_mat
                      tt0=tt0+mat%matrix_compr(ind0)*psitwork_c(i0j+0)
    
                      jjorb1=collcom%indexrecvorbital_c(i0j+1) - iorb_shift
                      jorb1=moduloarray(jjorb1)
                      !jjorb1=mod(jjorb1-1,sparsemat%nfvctr)+1
                      !ind1 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb1,iiorb)
                      ind1 = sparsemat%matrixindex_in_compressed_fortransposed(jorb1,iorb)
                      !ind1 = get_transposed_index(sparsemat,jjorb1,iiorb)
                      ind1=ind1+ishift_mat
                      tt1=tt1+mat%matrix_compr(ind1)*psitwork_c(i0j+1)
    
                      jjorb2=collcom%indexrecvorbital_c(i0j+2) - iorb_shift
                      jorb2=moduloarray(jjorb2)
                      !jjorb2=mod(jjorb2-1,sparsemat%nfvctr)+1
                      !ind2 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb2,iiorb)
                      ind2 = sparsemat%matrixindex_in_compressed_fortransposed(jorb2,iorb)
                      !ind2 = get_transposed_index(sparsemat,jjorb2,iiorb)
                      ind2=ind2+ishift_mat
                      tt2=tt2+mat%matrix_compr(ind2)*psitwork_c(i0j+2)
    
                      jjorb3=collcom%indexrecvorbital_c(i0j+3) - iorb_shift
                      jorb3=moduloarray(jjorb3)
                      !jjorb3=mod(jjorb3-1,sparsemat%nfvctr)+1
                      !ind3 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb3,iiorb)
                      ind3 = sparsemat%matrixindex_in_compressed_fortransposed(jorb3,iorb)
                      !ind3 = get_transposed_index(sparsemat,jjorb3,iiorb)
                      ind3=ind3+ishift_mat
                      tt3=tt3+mat%matrix_compr(ind3)*psitwork_c(i0j+3)
    
                      jjorb4=collcom%indexrecvorbital_c(i0j+4) - iorb_shift
                      jorb4=moduloarray(jjorb4)
                      !jjorb4=mod(jjorb4-1,sparsemat%nfvctr)+1
                      !ind4 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb4,iiorb)
                      ind4 = sparsemat%matrixindex_in_compressed_fortransposed(jorb4,iorb)
                      !ind4 = get_transposed_index(sparsemat,jjorb4,iiorb)
                      ind4=ind4+ishift_mat
                      tt4=tt4+mat%matrix_compr(ind4)*psitwork_c(i0j+4)
    
                      jjorb5=collcom%indexrecvorbital_c(i0j+5) - iorb_shift
                      jorb5=moduloarray(jjorb5)
                      !jjorb5=mod(jjorb5-1,sparsemat%nfvctr)+1
                      !ind5 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb5,iiorb)
                      ind5 = sparsemat%matrixindex_in_compressed_fortransposed(jorb5,iorb)
                      !ind5 = get_transposed_index(sparsemat,jjorb5,iiorb)
                      ind5=ind5+ishift_mat
                      tt5=tt5+mat%matrix_compr(ind5)*psitwork_c(i0j+5)
    
                      jjorb6=collcom%indexrecvorbital_c(i0j+6) - iorb_shift
                      jorb6=moduloarray(jjorb6)
                      !jjorb6=mod(jjorb6-1,sparsemat%nfvctr)+1
                      !ind6 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb6,iiorb)
                      ind6 = sparsemat%matrixindex_in_compressed_fortransposed(jorb6,iorb)
                      !ind6 = get_transposed_index(sparsemat,jjorb6,iiorb)
                      ind6=ind6+ishift_mat
                      tt6=tt6+mat%matrix_compr(ind6)*psitwork_c(i0j+6)
                  end do
                  psit_c(i0i)=psit_c(i0i)+tt0+tt1+tt2+tt3+tt4+tt5+tt6
              end do
          end do
          !$omp end do
    
          !$omp do schedule(static,1)
          do ipt=1,collcom%nptsp_f 
              ii=collcom%norb_per_gridpoint_f(ipt) 
              i0 = collcom%isptsp_f(ipt) + (ispin-1)*collcom%ndimind_f/sparsemat%nspin
              do i=1,ii
                  i0i=i0+i
                  i07i=7*i0i
                  iiorb=collcom%indexrecvorbital_f(i0i) - iorb_shift
                  iorb=moduloarray(iiorb)
                  !iiorb=mod(iiorb-1,sparsemat%nfvctr)+1
                  m=mod(ii,7)
                  tt00=0.d0 ; tt01=0.d0 ; tt02=0.d0 ; tt03=0.d0 ; tt04=0.d0 ; tt05=0.d0 ; tt06=0.d0
                  tt10=0.d0 ; tt11=0.d0 ; tt12=0.d0 ; tt13=0.d0 ; tt14=0.d0 ; tt15=0.d0 ; tt16=0.d0
                  tt20=0.d0 ; tt21=0.d0 ; tt22=0.d0 ; tt23=0.d0 ; tt24=0.d0 ; tt25=0.d0 ; tt26=0.d0
                  tt30=0.d0 ; tt31=0.d0 ; tt32=0.d0 ; tt33=0.d0 ; tt34=0.d0 ; tt35=0.d0 ; tt36=0.d0
                  tt40=0.d0 ; tt41=0.d0 ; tt42=0.d0 ; tt43=0.d0 ; tt44=0.d0 ; tt45=0.d0 ; tt46=0.d0
                  tt50=0.d0 ; tt51=0.d0 ; tt52=0.d0 ; tt53=0.d0 ; tt54=0.d0 ; tt55=0.d0 ; tt56=0.d0
                  tt60=0.d0 ; tt61=0.d0 ; tt62=0.d0 ; tt63=0.d0 ; tt64=0.d0 ; tt65=0.d0 ; tt66=0.d0
                  if(m/=0) then
                      do j=1,m
                          i0j=i0+j
                          i07j=7*i0j
                          jjorb=collcom%indexrecvorbital_f(i0j) - iorb_shift
                          jorb=moduloarray(jjorb)
                          !jjorb=mod(jjorb-1,sparsemat%nfvctr)+1
                          !ind0 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
                          ind0 = sparsemat%matrixindex_in_compressed_fortransposed(jorb,iorb)
                          !ind0 = get_transposed_index(sparsemat,jjorb,iiorb)
                          ind0=ind0+ishift_mat
                          tt06 = tt06 + mat%matrix_compr(ind0)*psitwork_f(i07j-6)
                          tt05 = tt05 + mat%matrix_compr(ind0)*psitwork_f(i07j-5)
                          tt04 = tt04 + mat%matrix_compr(ind0)*psitwork_f(i07j-4)
                          tt03 = tt03 + mat%matrix_compr(ind0)*psitwork_f(i07j-3)
                          tt02 = tt02 + mat%matrix_compr(ind0)*psitwork_f(i07j-2)
                          tt01 = tt01 + mat%matrix_compr(ind0)*psitwork_f(i07j-1)
                          tt00 = tt00 + mat%matrix_compr(ind0)*psitwork_f(i07j-0)
                      end do
                  end if
                  do j=m+1,ii,7
                      i0j=i0+j
                      i07j=7*i0j
                      jjorb0=collcom%indexrecvorbital_f(i0j+0) - iorb_shift
                      jorb0=moduloarray(jjorb0)
                      !jjorb0=mod(jjorb0-1,sparsemat%nfvctr)+1
                      !ind0 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb0,iiorb)
                      ind0 = sparsemat%matrixindex_in_compressed_fortransposed(jorb0,iorb)
                      !ind0 = get_transposed_index(sparsemat,jjorb0,iiorb)
                      ind0=ind0+ishift_mat
                      tt06 = tt06 + mat%matrix_compr(ind0)*psitwork_f(i07j-6)
                      tt05 = tt05 + mat%matrix_compr(ind0)*psitwork_f(i07j-5)
                      tt04 = tt04 + mat%matrix_compr(ind0)*psitwork_f(i07j-4)
                      tt03 = tt03 + mat%matrix_compr(ind0)*psitwork_f(i07j-3)
                      tt02 = tt02 + mat%matrix_compr(ind0)*psitwork_f(i07j-2)
                      tt01 = tt01 + mat%matrix_compr(ind0)*psitwork_f(i07j-1)
                      tt00 = tt00 + mat%matrix_compr(ind0)*psitwork_f(i07j-0)
    
                      jjorb1=collcom%indexrecvorbital_f(i0j+1) - iorb_shift
                      jorb1=moduloarray(jjorb1)
                      !jjorb1=mod(jjorb1-1,sparsemat%nfvctr)+1
                      !ind1 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb1,iiorb)
                      ind1 = sparsemat%matrixindex_in_compressed_fortransposed(jorb1,iorb)
                      !ind1 = get_transposed_index(sparsemat,jjorb1,iiorb)
                      ind1=ind1+ishift_mat
                      tt16 = tt16 + mat%matrix_compr(ind1)*psitwork_f(i07j+1) !+1*7-6
                      tt15 = tt15 + mat%matrix_compr(ind1)*psitwork_f(i07j+2) !+1*7-5
                      tt14 = tt14 + mat%matrix_compr(ind1)*psitwork_f(i07j+3) !+1*7-4
                      tt13 = tt13 + mat%matrix_compr(ind1)*psitwork_f(i07j+4) !+1*7-3
                      tt12 = tt12 + mat%matrix_compr(ind1)*psitwork_f(i07j+5) !+1*7-2
                      tt11 = tt11 + mat%matrix_compr(ind1)*psitwork_f(i07j+6) !+1*7-1
                      tt10 = tt10 + mat%matrix_compr(ind1)*psitwork_f(i07j+7) !+1*7-0
    
                      jjorb2=collcom%indexrecvorbital_f(i0j+2) - iorb_shift
                      jorb2=moduloarray(jjorb2)
                      !jjorb2=mod(jjorb2-1,sparsemat%nfvctr)+1
                      !ind2 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb2,iiorb)
                      ind2 = sparsemat%matrixindex_in_compressed_fortransposed(jorb2,iorb)
                      !ind2 = get_transposed_index(sparsemat,jjorb2,iiorb)
                      ind2=ind2+ishift_mat
                      tt26 = tt26 + mat%matrix_compr(ind2)*psitwork_f(i07j+8) !+2*7-6
                      tt25 = tt25 + mat%matrix_compr(ind2)*psitwork_f(i07j+9) !+2*7-5
                      tt24 = tt24 + mat%matrix_compr(ind2)*psitwork_f(i07j+10) !+2*7-4
                      tt23 = tt23 + mat%matrix_compr(ind2)*psitwork_f(i07j+11) !+2*7-3
                      tt22 = tt22 + mat%matrix_compr(ind2)*psitwork_f(i07j+12) !+2*7-2
                      tt21 = tt21 + mat%matrix_compr(ind2)*psitwork_f(i07j+13) !+2*7-1
                      tt20 = tt20 + mat%matrix_compr(ind2)*psitwork_f(i07j+14) !+2*7-0
    
                      jjorb3=collcom%indexrecvorbital_f(i0j+3) - iorb_shift
                      jorb3=moduloarray(jjorb3)
                      !jjorb3=mod(jjorb3-1,sparsemat%nfvctr)+1
                      !ind3 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb3,iiorb)
                      ind3 = sparsemat%matrixindex_in_compressed_fortransposed(jorb3,iorb)
                      !ind3 = get_transposed_index(sparsemat,jjorb3,iiorb)
                      ind3=ind3+ishift_mat
                      tt36 = tt36 + mat%matrix_compr(ind3)*psitwork_f(i07j+15) !+3*7-6
                      tt35 = tt35 + mat%matrix_compr(ind3)*psitwork_f(i07j+16) !+3*7-5
                      tt34 = tt34 + mat%matrix_compr(ind3)*psitwork_f(i07j+17) !+3*7-4
                      tt33 = tt33 + mat%matrix_compr(ind3)*psitwork_f(i07j+18) !+3*7-3
                      tt32 = tt32 + mat%matrix_compr(ind3)*psitwork_f(i07j+19) !+3*7-2
                      tt31 = tt31 + mat%matrix_compr(ind3)*psitwork_f(i07j+20) !+3*7-1
                      tt30 = tt30 + mat%matrix_compr(ind3)*psitwork_f(i07j+21) !+3*7-0
    
                      jjorb4=collcom%indexrecvorbital_f(i0j+4) - iorb_shift
                      jorb4=moduloarray(jjorb4)
                      !jjorb4=mod(jjorb4-1,sparsemat%nfvctr)+1
                      !ind4 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb4,iiorb)
                      ind4 = sparsemat%matrixindex_in_compressed_fortransposed(jorb4,iorb)
                      !ind4 = get_transposed_index(sparsemat,jjorb4,iiorb)
                      ind4=ind4+ishift_mat
                      tt46 = tt46 + mat%matrix_compr(ind4)*psitwork_f(i07j+22) !+4*7-6
                      tt45 = tt45 + mat%matrix_compr(ind4)*psitwork_f(i07j+23) !+4*7-5
                      tt44 = tt44 + mat%matrix_compr(ind4)*psitwork_f(i07j+24) !+4*7-4
                      tt43 = tt43 + mat%matrix_compr(ind4)*psitwork_f(i07j+25) !+4*7-3
                      tt42 = tt42 + mat%matrix_compr(ind4)*psitwork_f(i07j+26) !+4*7-2
                      tt41 = tt41 + mat%matrix_compr(ind4)*psitwork_f(i07j+27) !+4*7-1
                      tt40 = tt40 + mat%matrix_compr(ind4)*psitwork_f(i07j+28) !+4*7-0
    
                      jjorb5=collcom%indexrecvorbital_f(i0j+5) - iorb_shift
                      jorb5=moduloarray(jjorb5)
                      !jjorb5=mod(jjorb5-1,sparsemat%nfvctr)+1
                      !ind5 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb5,iiorb)
                      ind5 = sparsemat%matrixindex_in_compressed_fortransposed(jorb5,iorb)
                      !ind5 = get_transposed_index(sparsemat,jjorb5,iiorb)
                      ind5=ind5+ishift_mat
                      tt56 = tt56 + mat%matrix_compr(ind5)*psitwork_f(i07j+29) !+5*7-6
                      tt55 = tt55 + mat%matrix_compr(ind5)*psitwork_f(i07j+30) !+5*7-5
                      tt54 = tt54 + mat%matrix_compr(ind5)*psitwork_f(i07j+31) !+5*7-4
                      tt53 = tt53 + mat%matrix_compr(ind5)*psitwork_f(i07j+32) !+5*7-3
                      tt52 = tt52 + mat%matrix_compr(ind5)*psitwork_f(i07j+33) !+5*7-2
                      tt51 = tt51 + mat%matrix_compr(ind5)*psitwork_f(i07j+34) !+5*7-1
                      tt50 = tt50 + mat%matrix_compr(ind5)*psitwork_f(i07j+35) !+5*7-0
    
                      jjorb6=collcom%indexrecvorbital_f(i0j+6) - iorb_shift
                      jorb6=moduloarray(jjorb6)
                      !jjorb6=mod(jjorb6-1,sparsemat%nfvctr)+1
                      !ind6 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb6,iiorb)
                      ind6 = sparsemat%matrixindex_in_compressed_fortransposed(jorb6,iorb)
                      !ind6 = get_transposed_index(sparsemat,jjorb6,iiorb)
                      ind6=ind6+ishift_mat
                      tt66 = tt66 + mat%matrix_compr(ind6)*psitwork_f(i07j+36) !+6*7-6
                      tt65 = tt65 + mat%matrix_compr(ind6)*psitwork_f(i07j+37) !+6*7-5
                      tt64 = tt64 + mat%matrix_compr(ind6)*psitwork_f(i07j+38) !+6*7-4
                      tt63 = tt63 + mat%matrix_compr(ind6)*psitwork_f(i07j+39) !+6*7-3
                      tt62 = tt62 + mat%matrix_compr(ind6)*psitwork_f(i07j+40) !+6*7-2
                      tt61 = tt61 + mat%matrix_compr(ind6)*psitwork_f(i07j+41) !+6*7-1
                      tt60 = tt60 + mat%matrix_compr(ind6)*psitwork_f(i07j+42) !+6*7-0
                  end do
                  psit_f(i07i-6) = psit_f(i07i-6) + tt06 + tt16 + tt26 + tt36 + tt46 + tt56 + tt66
                  psit_f(i07i-5) = psit_f(i07i-5) + tt05 + tt15 + tt25 + tt35 + tt45 + tt55 + tt65
                  psit_f(i07i-4) = psit_f(i07i-4) + tt04 + tt14 + tt24 + tt34 + tt44 + tt54 + tt64
                  psit_f(i07i-3) = psit_f(i07i-3) + tt03 + tt13 + tt23 + tt33 + tt43 + tt53 + tt63
                  psit_f(i07i-2) = psit_f(i07i-2) + tt02 + tt12 + tt22 + tt32 + tt42 + tt52 + tt62
                  psit_f(i07i-1) = psit_f(i07i-1) + tt01 + tt11 + tt21 + tt31 + tt41 + tt51 + tt61
                  psit_f(i07i-0) = psit_f(i07i-0) + tt00 + tt10 + tt20 + tt30 + tt40 + tt50 + tt60
              end do  
          end do
          !$omp end do
          !$omp end parallel
    
      end do spin_loop
    
      call f_free_ptr(moduloarray)

      call f_release_routine()
      call timing(iproc,'lincombtrans  ','OF')

      !contains

        !function get_transposed_index(jorb,iorb) res(ind)
        !    integer,intent(in) :: jorb, iorb
        !    integer :: ind
        !    integer :: jjorb,iiorb
        !    ! If iorb is smaller than the offset, add a periodic shift
        !    if (iorb<sparsemat%offset_matrixindex_in_compressed_fortransposed) then
        !        iiorb = iorb + sparsemat%nfvctr
        !    else
        !        iiorb = iorb
        !    end if
        !    if (jorb<sparsemat%offset_matrixindex_in_compressed_fortransposed) then
        !        jjorb = jorb + sparsemat%nfvctr
        !    else
        !        jjorb = jorb
        !    end if
        !    ind = sparsemat%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
        !end function get_transposed_index
    
    end subroutine build_linear_combination_transposed


    subroutine init_matrixindex_in_compressed_fortransposed(iproc, nproc, orbs, collcom, collcom_shamop, &
               collcom_sr, sparsemat)
      use module_base
      use module_types, only: orbitals_data
      use communications_base, only: comms_linear
      !use module_interfaces, except_this_one => init_matrixindex_in_compressed_fortransposed
      use sparsematrix_base, only: sparse_matrix
      use sparsematrix_init, only: matrixindex_in_compressed!compressed_index
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(orbitals_data),intent(in) :: orbs
      type(comms_linear),intent(in) :: collcom, collcom_shamop, collcom_sr
      type(sparse_matrix), intent(inout) :: sparsemat
      
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
    
      !$omp parallel do default(private) shared(sparsemat,orbs,imin,imax)
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



end module transposed_operations