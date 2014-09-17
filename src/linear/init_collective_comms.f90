!> @file
!! Intialization of the collective communications for the linear version
!! @author
!!    Copyright (C) 2011-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS



subroutine check_communications_locreg(iproc,nproc,orbs,nspin,Lzd,collcom,smat,mat,npsidim_orbs,npsidim_comp,check_overlap)
   use module_base!, only: wp, bigdft_mpi, mpi_sum, mpi_max, mpiallred
   use module_types, only: orbitals_data, local_zone_descriptors, linear_matrices
   use yaml_output
   use communications_base, only: comms_linear
   use communications, only: transpose_localized, untranspose_localized
   use sparsematrix_base, only : sparse_matrix, matrices, DENSE_PARALLEL
   use sparsematrix, only : compress_matrix_distributed
   !use dynamic_memory
   implicit none
   integer, intent(in) :: iproc,nproc,nspin,check_overlap
   type(orbitals_data), intent(in) :: orbs
   type(local_zone_descriptors), intent(in) :: lzd
   type(comms_linear), intent(in) :: collcom
   type(sparse_matrix),intent(inout) :: smat
   type(matrices),intent(inout) :: mat
   integer, intent(in) :: npsidim_orbs, npsidim_comp
   !local variables
   character(len=*), parameter :: subname='check_communications'
   integer, parameter :: ilog=6
   integer :: i,ispinor,iorb,indspin,i_stat,i_all,ikptsp
   integer :: ikpt,ierr,i0,ifine,ii,iiorb,ipt,jorb,indorb_tmp
   integer :: icomp,ispin
   !!$integer :: ipsi,ipsic,ipsif,ipsiworkc,ipsiworkf,jcomp,jkpt
   real(wp) :: psival,maxdiff,tt
   real(wp), dimension(:), allocatable :: psi,psit_c,psit_f
   real(wp), dimension(:,:), allocatable :: checksum
   real(wp) :: epsilon,tol
   logical :: abort, isoverlap
   integer :: jjorb, ilr, jlr, ldim, gdim, iispin, jjspin, ist, niorb, njorb, jjjorb, is, ie
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
    
       call transpose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, psi, psit_c, psit_f, lzd)
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
          call MPI_BARRIER(bigdft_mpi%mpi_comm, ierr)
          call mpiallred(checksum(1,1),2*orbs%norb*orbs%nspinor,MPI_SUM,bigdft_mpi%mpi_comm)
       end if
    
       if (iproc==0) then
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
       end if
       if (iproc==0) call yaml_map('Maxdiff for transpose (checksum)',&
            maxdiff,fmt='(1pe25.17)')
    
    
       abort = .false.
       if (abs(maxdiff) >tol) then
          call yaml_comment('ERROR (Transposition): process'//trim(yaml_toa(iproc))//&
               ' found an error of:'//trim(yaml_toa(maxdiff,fmt='(1pe15.7)')))
          !call yaml_map('Some wrong results in',(/jkpt,jorb,jcomp/),fmt='(i8)')
          abort=.true.
       end if
    
       if (abort) call MPI_ABORT(bigdft_mpi%mpi_comm,10,ierr)
    
    
       !@NEW: check the calculation of the overlap matrices #############
       if (check_overlap > 1) then
           call calculate_overlap_transposed(iproc, nproc, orbs, collcom, psit_c, &
                psit_c, psit_f, psit_f, smat, mat)
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
                   call to_zero(gdim, psiig(1))
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
                       call to_zero(gdim, psijg(1))
                       !!write(4200+iproc,'(a,2i8,l4)') 'iproc, jlr, associated(lzd%llr(jlr)%wfd%keygloc)', iproc, jlr, associated(lzd%llr(jlr)%wfd%keygloc)
                       call Lpsi_to_global2(iproc, ldim, gdim, orbs%norb, orbs%nspinor, 1, lzd%glr, &
                                            lzd%llr(jlr), psij, psijg)
                       matp(jjjorb,niorb)=ddot(gdim, psiig, 1, psijg, 1)
                       call f_free(psij)
                   end do
                   call f_free(psii)
               end do
               ist=(ispin-1)*smat%nvctr+1
               call compress_matrix_distributed(iproc, nproc, smat, DENSE_PARALLEL, matp, mat_compr(ist))
           end do
           maxdiff=0.d0
           call f_free(psiig)
           call f_free(psijg)
           call f_free(matp)
           do i=1,smat%nvctr
               maxdiff=max(abs(mat_compr(i)-mat%matrix_compr(i)),maxdiff)
               !write(8000+iproc,'(a,i7,2es15.5)') 'i, mat_compr(i), mat%matrix_compr(i)', &
               !    i, mat_compr(i), mat%matrix_compr(i)
           end do
           call f_free(mat_compr)
           if (iproc==0) call yaml_map('Maxdiff for overlap calculation',maxdiff,fmt='(1es25.17)')
       end if
       !@END NEW ########################################################
    
    
       call untranspose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, psit_c, psit_f, psi, lzd)
    
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
          call yaml_comment('ERROR (Inverse Transposition): process'//trim(yaml_toa(iproc))//&
               ' found an error of:'//trim(yaml_toa(maxdiff,fmt='(1pe15.7)')))
          abort = .true.
       end if
    
       if (abort) call MPI_ABORT(bigdft_mpi%mpi_comm,ierr)
    
       if (nproc > 1) then
          call MPI_BARRIER(bigdft_mpi%mpi_comm, ierr)
          call mpiallred(maxdiff,1,MPI_MAX,bigdft_mpi%mpi_comm)
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


subroutine calculate_overlap_transposed(iproc, nproc, orbs, collcom, &
           psit_c1, psit_c2, psit_f1, psit_f2, smat, ovrlp)
  use module_base
  use module_types
  use sparsematrix_base, only: sparse_matrix
  use sparsematrix, only : orb_from_index
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
  real(kind=8) :: tt00, tt01, tt02, tt03, tt04, tt05, tt06
  real(kind=8) :: tt10, tt11, tt12, tt13, tt14, tt15, tt16
  real(kind=8) :: tt20, tt21, tt22, tt23, tt24, tt25, tt26
  real(kind=8) :: tt30, tt31, tt32, tt33, tt34, tt35, tt36
  real(kind=8) :: tt40, tt41, tt42, tt43, tt44, tt45, tt46
  real(kind=8) :: tt50, tt51, tt52, tt53, tt54, tt55, tt56
  real(kind=8) :: tt60, tt61, tt62, tt63, tt64, tt65, tt66
  integer,dimension(:),allocatable :: n
  !$ integer  :: omp_get_thread_num,omp_get_max_threads
  real(kind=8) :: totops
  integer :: avops, ops, opsn
  integer, allocatable, dimension(:) :: numops
  logical :: ifnd, jfnd
  integer :: iorb, jorb, imat, iseg, iorb_shift
  integer,dimension(2) :: irowcol

  call timing(iproc,'ovrlptransComp','ON') !lr408t

  call f_routine(id='calculate_overlap_transposed')

  call to_zero(smat%nvctr*smat%nspin, ovrlp%matrix_compr(1))

  ! WARNING: METHOD 2 NOT EXTENSIVELY TESTED
  method_if: if (collcom%imethod_overlap==2) then

      !!!iicnt=0
      !!!sm_it = iterator(collcom)
      !!!do while(valid(sm_it))
      !!!icnt=icnt+1
      !!!call get_position(sm_it,shift=ind0(:),iorb=i0i,jorb=i0j)
      !!!call ge_orbitals(sm_it,iiorb,ijorb)
      !!!call get_ind0(iiorb+(norb)*jjorb,ind0)
      !!!ovrlp%matrix_compr(ind0) = ovrlp%matrix_compr(ind0) + psit_c1(i0i)*psit_c2(i0j)
      !!!call next(sm_it)
      !!!end do
      iorb_shift=(ispin-1)*smat%nfvctr    
      do ispin=1,smat%nspin
        do iseg=1,smat%nseg
          imat=smat%keyv(iseg)
          do j=smat%keyg(1,iseg),smat%keyg(2,iseg)
            !call get_orbs(smat,i,iorb,jorb) !lookup on work array of size smat%nvctr 
            !iorb=smat%orb_from_index(2,imat)
            !jorb=smat%orb_from_index(1,imat)
            irowcol = orb_from_index(smat, j)
            ovrlp%matrix_compr(imat)=0.0_wp
      
            do ipt=1,collcom%nptsp_c
              ii=collcom%norb_per_gridpoint_c(ipt)
              i0 = collcom%isptsp_c(ipt) + (ispin-1)*collcom%ndimind_c/smat%nspin
              ifnd=.false.
              jfnd=.false.
              do i=1,ii
                iiorb=collcom%indexrecvorbital_c(i0+i) - iorb_shift
                !iiorb=mod(iiorb-1,smat%nfvctr)+1
                if (iiorb == irowcol(1)) then        
                   ifnd=.true.
                   i0i=i0+i
                   !i0i=collcom%iextract_c(i0+i)
                end if 
                if (iiorb == irowcol(2)) then
                    jfnd=.true.
                    i0j=i0+i
                    !i0j=collcom%iextract_c(i0+i)
                end if
                if (.not. (jfnd .and. ifnd)) cycle
                ovrlp%matrix_compr(imat) = ovrlp%matrix_compr(imat) +  psit_c1(i0i)*psit_c2(i0j)
                if (jfnd .and. ifnd) exit
              end do
            end do
      
            do ipt=1,collcom%nptsp_f
              ii=collcom%norb_per_gridpoint_f(ipt)
              i0 = collcom%isptsp_f(ipt) + (ispin-1)*collcom%ndimind_f/smat%nspin
              ifnd=.false.
              jfnd=.false.
              do i=1,ii
                iiorb=collcom%indexrecvorbital_f(i0+i) - iorb_shift
                !iiorb=mod(iiorb-1,smat%nfvctr)+1
                if (iiorb == irowcol(1)) then        
                   ifnd=.true.
                   i0i=i0+i
                   !i0i=collcom%iextract_f(i0+i)
                end if 
                if (iiorb == irowcol(2)) then
                    jfnd=.true.
                    i0j=i0+i
                    !i0j=collcom%iextract_f(i0+i)
                end if
                if (.not. (jfnd .and. ifnd)) cycle
                i07i=7*i0i
                i07j=7*i0j
                ovrlp%matrix_compr(imat) = ovrlp%matrix_compr(imat) +  psit_f1(i07i-6)*psit_f2(i07j-6)
                ovrlp%matrix_compr(imat) = ovrlp%matrix_compr(imat) +  psit_f1(i07i-5)*psit_f2(i07j-5)
                ovrlp%matrix_compr(imat) = ovrlp%matrix_compr(imat) +  psit_f1(i07i-4)*psit_f2(i07j-4)
                ovrlp%matrix_compr(imat) = ovrlp%matrix_compr(imat) +  psit_f1(i07i-3)*psit_f2(i07j-3)
                ovrlp%matrix_compr(imat) = ovrlp%matrix_compr(imat) +  psit_f1(i07i-2)*psit_f2(i07j-2)
                ovrlp%matrix_compr(imat) = ovrlp%matrix_compr(imat) +  psit_f1(i07i-1)*psit_f2(i07j-1)
                ovrlp%matrix_compr(imat) = ovrlp%matrix_compr(imat) +  psit_f1(i07i-0)*psit_f2(i07j-0)
                if (jfnd .and. ifnd) exit
              end do
            end do
            imat=imat+1
           end do 
         end do    
      end do

  else if (collcom%imethod_overlap==1) then method_if

      !only optimized for spin=1 for now
      ispin=1
      nthreads=1
      !$  nthreads = OMP_GET_max_threads()
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
         totops=sum(numops)
         avops=nint(totops/dble(nthreads))
         jjorb=1
         do i=1,nthreads
            ops=0
            do j=jjorb,orbs%norb
               opsn=ops+numops(j)
               if (opsn>=avops) then
                  if ((opsn-avops)<(avops-ops)) then
                     n(i)=j
                     jjorb=j+1
                     totops=totops-opsn
                  else
                     n(i)=j-1
                     jjorb=j
                     totops=totops-ops
                  end if
                  exit
               end if
               ops=opsn
            end do
            if (i/=nthreads) then
               avops=nint(totops/dble(nthreads-i))
            end if
         end do
         call f_free(numops)
      end if

      n(nthreads)=orbs%norb
    

      !$omp parallel default(none) &
      !$omp shared(collcom, smat, ovrlp, psit_c1, psit_c2, psit_f1, psit_f2, n) &
      !$omp private(tid, ispin, iend, istart, ipt, ii, i0, i, iiorb, m, j, i0j, jjorb, ishift_mat, iorb_shift, ind0) &
      !$omp private(jjorb0, jjorb1, ind1, jjorb2, ind2, jjorb3, ind3, jjorb4, ind4, jjorb5, ind5, jjorb6, ind6) &
      !$omp private(i0i, i07i, i07j, tt06, tt05, tt04, tt03, tt02, tt01, tt00) &
      !$omp private(tt16, tt15, tt14, tt13, tt12, tt11, tt10) & 
      !$omp private(tt26, tt25, tt24, tt23, tt22, tt21, tt20) &
      !$omp private(tt36, tt35, tt34, tt33, tt32, tt31, tt30) &
      !$omp private(tt46, tt45, tt44, tt43, tt42, tt41, tt40) &
      !$omp private(tt56, tt55, tt54, tt53, tt52, tt51, tt50) &
      !$omp private(tt66, tt65, tt64, tt63, tt62, tt61, tt60)
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
    
          ishift_mat=(ispin-1)*smat%nvctr
          iorb_shift=(ispin-1)*smat%nfvctr
          if (collcom%nptsp_c>0) then
    
              do ipt=1,collcom%nptsp_c 
                  ii=collcom%norb_per_gridpoint_c(ipt) 
                  i0 = collcom%isptsp_c(ipt) + (ispin-1)*collcom%ndimind_c/smat%nspin
                  do i=1,ii
                      i0i=i0+i
                      iiorb=collcom%indexrecvorbital_c(i0i) - iorb_shift
                      !iiorb=mod(iiorb-1,smat%nfvctr)+1
                      if(iiorb < istart .or. iiorb > iend) cycle
                      m=mod(ii,7)
                      if(m/=0) then
                          do j=1,m
                              i0j=i0+j
                              jjorb=collcom%indexrecvorbital_c(i0j) - iorb_shift
                              !jjorb=mod(jjorb-1,smat%nfvctr)+1
                              !ind0 = ishift_mat + smat%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
                              ind0 = smat%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
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
                          !jjorb0=mod(jjorb0-1,smat%nfvctr)+1
                          !ind0 = ishift_mat + smat%matrixindex_in_compressed_fortransposed(jjorb0,iiorb)
                          ind0 = smat%matrixindex_in_compressed_fortransposed(jjorb0,iiorb)
                          ind0=ind0+ishift_mat
                          ovrlp%matrix_compr(ind0) = ovrlp%matrix_compr(ind0) + psit_c1(i0i)*psit_c2(i0j+0)
    
                          jjorb1=collcom%indexrecvorbital_c(i0j+1) - iorb_shift
                          !jjorb1=mod(jjorb1-1,smat%nfvctr)+1
                          !ind1 = ishift_mat + smat%matrixindex_in_compressed_fortransposed(jjorb1,iiorb)
                          ind1 = smat%matrixindex_in_compressed_fortransposed(jjorb1,iiorb)
                          ind1=ind1+ishift_mat
                          ovrlp%matrix_compr(ind1) = ovrlp%matrix_compr(ind1) + psit_c1(i0i)*psit_c2(i0j+1)
    
                          jjorb2=collcom%indexrecvorbital_c(i0j+2) - iorb_shift
                          !jjorb2=mod(jjorb2-1,smat%nfvctr)+1
                          !ind2 = ishift_mat + smat%matrixindex_in_compressed_fortransposed(jjorb2,iiorb)
                          ind2 = smat%matrixindex_in_compressed_fortransposed(jjorb2,iiorb)
                          ind2=ind2+ishift_mat
                          ovrlp%matrix_compr(ind2) = ovrlp%matrix_compr(ind2) + psit_c1(i0i)*psit_c2(i0j+2)
    
                          jjorb3=collcom%indexrecvorbital_c(i0j+3) - iorb_shift
                          !jjorb3=mod(jjorb3-1,smat%nfvctr)+1
                          !ind3 = ishift_mat + smat%matrixindex_in_compressed_fortransposed(jjorb3,iiorb)
                          ind3 = smat%matrixindex_in_compressed_fortransposed(jjorb3,iiorb)
                          ind3=ind3+ishift_mat
                          ovrlp%matrix_compr(ind3) = ovrlp%matrix_compr(ind3) + psit_c1(i0i)*psit_c2(i0j+3)
    
                          jjorb4=collcom%indexrecvorbital_c(i0j+4) - iorb_shift
                          !jjorb4=mod(jjorb4-1,smat%nfvctr)+1
                          !ind4 = ishift_mat + smat%matrixindex_in_compressed_fortransposed(jjorb4,iiorb)
                          ind4 = smat%matrixindex_in_compressed_fortransposed(jjorb4,iiorb)
                          ind4=ind4+ishift_mat
                          ovrlp%matrix_compr(ind4) = ovrlp%matrix_compr(ind4) + psit_c1(i0i)*psit_c2(i0j+4)
    
                          jjorb5=collcom%indexrecvorbital_c(i0j+5) - iorb_shift
                          !jjorb5=mod(jjorb5-1,smat%nfvctr)+1
                          !ind5 = ishift_mat + smat%matrixindex_in_compressed_fortransposed(jjorb5,iiorb)
                          ind5 = smat%matrixindex_in_compressed_fortransposed(jjorb5,iiorb)
                          ind5=ind5+ishift_mat
                          ovrlp%matrix_compr(ind5) = ovrlp%matrix_compr(ind5) + psit_c1(i0i)*psit_c2(i0j+5)
    
                          jjorb6=collcom%indexrecvorbital_c(i0j+6) - iorb_shift
                          !jjorb6=mod(jjorb6-1,smat%nfvctr)+1
                          !ind6 = ishift_mat + smat%matrixindex_in_compressed_fortransposed(jjorb6,iiorb)
                          ind6 = smat%matrixindex_in_compressed_fortransposed(jjorb6,iiorb)
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
         totops=sum(numops)
         avops=nint(totops/dble(nthreads))
         jjorb=1
         do i=1,nthreads
            ops=0
            do j=jjorb,orbs%norb
               opsn=ops+numops(j)
               if (opsn>=avops) then
                  if ((opsn-avops)<(avops-ops)) then
                     n(i)=j
                     jjorb=j+1
                     totops=totops-opsn
                  else
                     n(i)=j-1
                     jjorb=j
                     totops=totops-ops
                  end if
                  exit
               end if
               ops=opsn
            end do
            if (i/=nthreads) then
               avops=nint(totops/dble(nthreads-i))
            end if
         end do
         call f_free(numops)
      end if    

      n(nthreads)=orbs%norb
    

      !$omp parallel default(none) &
      !$omp shared(collcom, smat, ovrlp, psit_c1, psit_c2, psit_f1, psit_f2, n) &
      !$omp private(tid, ispin, iend, istart, ipt, ii, i0, i, iiorb, m, j, i0j, jjorb, ishift_mat, iorb_shift, ind0) &
      !$omp private(jjorb0, jjorb1, ind1, jjorb2, ind2, jjorb3, ind3, jjorb4, ind4, jjorb5, ind5, jjorb6, ind6) &
      !$omp private(i0i, i07i, i07j, tt06, tt05, tt04, tt03, tt02, tt01, tt00) &
      !$omp private(tt16, tt15, tt14, tt13, tt12, tt11, tt10) & 
      !$omp private(tt26, tt25, tt24, tt23, tt22, tt21, tt20) &
      !$omp private(tt36, tt35, tt34, tt33, tt32, tt31, tt30) &
      !$omp private(tt46, tt45, tt44, tt43, tt42, tt41, tt40) &
      !$omp private(tt56, tt55, tt54, tt53, tt52, tt51, tt50) &
      !$omp private(tt66, tt65, tt64, tt63, tt62, tt61, tt60)
      tid=0
      !$ tid = OMP_GET_THREAD_NUM()
      iend=n(tid+1)
      if (tid==0) then
         istart=1
      else
         istart=n(tid)+1
      end if
    
      spin_loopf: do ispin=1,smat%nspin
    
          ishift_mat=(ispin-1)*smat%nvctr
          iorb_shift=(ispin-1)*smat%nfvctr

          if (collcom%nptsp_f>0) then
              do ipt=1,collcom%nptsp_f 
                  ii=collcom%norb_per_gridpoint_f(ipt) 
                  i0 = collcom%isptsp_f(ipt) + (ispin-1)*collcom%ndimind_f/smat%nspin
                  do i=1,ii
                      i0i=i0+i
                      iiorb=collcom%indexrecvorbital_f(i0i) - iorb_shift
                      !iiorb=mod(iiorb-1,smat%nfvctr)+1
                      if(iiorb < istart .or. iiorb > iend) cycle
                      i07i=7*i0i
                      m=mod(ii,7)
                      if(m/=0) then
                          do j=1,m
                              i0j=i0+j
                              i07j=7*i0j
                              jjorb0=collcom%indexrecvorbital_f(i0j) - iorb_shift
                              !jjorb0=mod(jjorb0-1,smat%nfvctr)+1
                              !ind0 = ishift_mat + smat%matrixindex_in_compressed_fortransposed(jjorb0,iiorb)
                              ind0 = smat%matrixindex_in_compressed_fortransposed(jjorb0,iiorb)
                              ind0=ind0+ishift_mat
                              tt06 = psit_f1(i07i-6)*psit_f2(i07j-6)
                              tt05 = psit_f1(i07i-5)*psit_f2(i07j-5)
                              tt04 = psit_f1(i07i-4)*psit_f2(i07j-4)
                              tt03 = psit_f1(i07i-3)*psit_f2(i07j-3)
                              tt02 = psit_f1(i07i-2)*psit_f2(i07j-2)
                              tt01 = psit_f1(i07i-1)*psit_f2(i07j-1)
                              tt00 = psit_f1(i07i-0)*psit_f2(i07j-0)
    
                              ovrlp%matrix_compr(ind0) = ovrlp%matrix_compr(ind0) + tt06 + tt05 + tt04 + tt03 + tt02 + tt01 + tt00
                          end do
                      end if
                      do j=m+1,ii,7
                          i0j=i0+j
                          i07j=7*i0j
                          jjorb0=collcom%indexrecvorbital_f(i0j+0) - iorb_shift
                          !jjorb0=mod(jjorb0-1,smat%nfvctr)+1
                          !ind0 = ishift_mat + smat%matrixindex_in_compressed_fortransposed(jjorb0,iiorb)
                          ind0 = smat%matrixindex_in_compressed_fortransposed(jjorb0,iiorb)
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
                          !jjorb1=mod(jjorb1-1,smat%nfvctr)+1
                          !ind1 = ishift_mat + smat%matrixindex_in_compressed_fortransposed(jjorb1,iiorb)
                          ind1 = smat%matrixindex_in_compressed_fortransposed(jjorb1,iiorb)
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
                          !jjorb2=mod(jjorb2-1,smat%nfvctr)+1
                          !ind2 = ishift_mat + smat%matrixindex_in_compressed_fortransposed(jjorb2,iiorb)
                          ind2 = smat%matrixindex_in_compressed_fortransposed(jjorb2,iiorb)
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
                          !jjorb3=mod(jjorb3-1,smat%nfvctr)+1
                          !ind3 = ishift_mat + smat%matrixindex_in_compressed_fortransposed(jjorb3,iiorb)
                          ind3 = smat%matrixindex_in_compressed_fortransposed(jjorb3,iiorb)
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
                          !jjorb4=mod(jjorb4-1,smat%nfvctr)+1
                          !ind4 = ishift_mat + smat%matrixindex_in_compressed_fortransposed(jjorb4,iiorb)
                          ind4 = smat%matrixindex_in_compressed_fortransposed(jjorb4,iiorb)
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
                          !jjorb5=mod(jjorb5-1,smat%nfvctr)+1
                          !ind5 = ishift_mat + smat%matrixindex_in_compressed_fortransposed(jjorb5,iiorb)
                          ind5 = smat%matrixindex_in_compressed_fortransposed(jjorb5,iiorb)
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
                          !jjorb6=mod(jjorb6-1,smat%nfvctr)+1
                          !ind6 = ishift_mat + smat%matrixindex_in_compressed_fortransposed(jjorb6,iiorb)
                          ind6 = smat%matrixindex_in_compressed_fortransposed(jjorb6,iiorb)
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

  else method_if
      stop 'wrong value of imethod_if'
  end if method_if

  call timing(iproc,'ovrlptransComp','OF') !lr408t

  call timing(iproc,'ovrlptransComm','ON') !lr408t

  if(nproc > 1) then
      call mpiallred(ovrlp%matrix_compr(1), smat%nvctr*smat%nspin, mpi_sum, bigdft_mpi%mpi_comm)
  end if


  smat%can_use_dense=.false.

  call f_release_routine()
  call timing(iproc,'ovrlptransComm','OF') !lr408t

end subroutine calculate_overlap_transposed



subroutine calculate_pulay_overlap(iproc, nproc, orbs1, orbs2, collcom1, collcom2, psit_c1, psit_c2, psit_f1, psit_f2, ovrlp)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(orbitals_data),intent(in) :: orbs1, orbs2
  type(comms_linear),intent(in) :: collcom1, collcom2
  real(kind=8),dimension(collcom1%ndimind_c),intent(in) :: psit_c1
  real(kind=8),dimension(collcom2%ndimind_c),intent(in) :: psit_c2
  real(kind=8),dimension(7*collcom1%ndimind_f),intent(in) :: psit_f1
  real(kind=8),dimension(7*collcom2%ndimind_f),intent(in) :: psit_f2
  real(kind=8),dimension(orbs1%norb,orbs2%norb),intent(out) :: ovrlp
  
  ! Local variables
  integer :: i0, j0, ipt, ii, iiorb, j, jj, jjorb, i, ierr  

  call timing(iproc,'ovrlptransComp','ON') !lr408t
  call to_zero(orbs1%norb*orbs2%norb, ovrlp(1,1))
  if(collcom1%nptsp_c/=collcom2%nptsp_c) then
      write(*,'(a,i0,a)') 'ERROR on process ',iproc,': collcom1%nptsp_c/=collcom2%nptsp_c'
      stop
  end if
  if(collcom1%nptsp_f/=collcom2%nptsp_f) then
      write(*,'(a,i0,a)') 'ERROR on process ',iproc,': collcom1%nptsp_f/=collcom2%nptsp_f'
      stop
  end if

  i0=0
  j0=0
  do ipt=1,collcom1%nptsp_c 
      ii=collcom1%norb_per_gridpoint_c(ipt)
      jj=collcom2%norb_per_gridpoint_c(ipt)
      do i=1,ii
          iiorb=collcom1%indexrecvorbital_c(i0+i)
          do j=1,jj
              jjorb=collcom2%indexrecvorbital_c(j0+j)
              ovrlp(iiorb,jjorb)=ovrlp(iiorb,jjorb)+psit_c1(i0+i)*psit_c2(j0+j)
          end do
      end do
      i0=i0+ii
      j0=j0+jj
  end do

  i0=0
  j0=0
  do ipt=1,collcom1%nptsp_f 
      ii=collcom1%norb_per_gridpoint_f(ipt)
      jj=collcom2%norb_per_gridpoint_f(ipt)
      do i=1,ii
          iiorb=collcom1%indexrecvorbital_f(i0+i)
          do j=1,jj
              jjorb=collcom2%indexrecvorbital_f(j0+j)
              ovrlp(iiorb,jjorb)=ovrlp(iiorb,jjorb)+psit_f1(7*(i0+i)-6)*psit_f2(7*(j0+j)-6)
              ovrlp(iiorb,jjorb)=ovrlp(iiorb,jjorb)+psit_f1(7*(i0+i)-5)*psit_f2(7*(j0+j)-5)
              ovrlp(iiorb,jjorb)=ovrlp(iiorb,jjorb)+psit_f1(7*(i0+i)-4)*psit_f2(7*(j0+j)-4)
              ovrlp(iiorb,jjorb)=ovrlp(iiorb,jjorb)+psit_f1(7*(i0+i)-3)*psit_f2(7*(j0+j)-3)
              ovrlp(iiorb,jjorb)=ovrlp(iiorb,jjorb)+psit_f1(7*(i0+i)-2)*psit_f2(7*(j0+j)-2)
              ovrlp(iiorb,jjorb)=ovrlp(iiorb,jjorb)+psit_f1(7*(i0+i)-1)*psit_f2(7*(j0+j)-1)
              ovrlp(iiorb,jjorb)=ovrlp(iiorb,jjorb)+psit_f1(7*(i0+i)-0)*psit_f2(7*(j0+j)-0)
          end do
      end do
      i0=i0+ii
      j0=j0+jj
  end do

  call timing(iproc,'ovrlptransComp','OF') !lr408t

  call timing(iproc,'ovrlptransComm','ON') !lr408t

  if(nproc>1) then
      call mpiallred(ovrlp(1,1), orbs1%norb*orbs2%norb, mpi_sum, bigdft_mpi%mpi_comm)
  end if
  call timing(iproc,'ovrlptransComm','OF') !lr408t
end subroutine calculate_pulay_overlap

subroutine build_linear_combination_transposed(collcom, sparsemat, mat, psitwork_c, psitwork_f, &
     reset, psit_c, psit_f, iproc)
  use module_base
  use module_types
  use sparsematrix_base, only: sparse_matrix
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
  integer :: ind4, ind5, ind6, jjorb0, jjorb1, jjorb2, jjorb3, jjorb4, jjorb5, jjorb6, ispin, ishift_mat
  real(kind=8) :: tt0, tt1, tt2, tt3, tt4, tt5, tt6
  real(kind=8) :: tt00, tt01, tt02, tt03, tt04, tt05, tt06
  real(kind=8) :: tt10, tt11, tt12, tt13, tt14, tt15, tt16
  real(kind=8) :: tt20, tt21, tt22, tt23, tt24, tt25, tt26
  real(kind=8) :: tt30, tt31, tt32, tt33, tt34, tt35, tt36
  real(kind=8) :: tt40, tt41, tt42, tt43, tt44, tt45, tt46
  real(kind=8) :: tt50, tt51, tt52, tt53, tt54, tt55, tt56
  real(kind=8) :: tt60, tt61, tt62, tt63, tt64, tt65, tt66

  call f_routine(id='build_linear_combination_transposed')
  call timing(iproc,'lincombtrans  ','ON') !lr408t
  if(reset) then
      if(collcom%ndimind_c>0) call to_zero(collcom%ndimind_c, psit_c(1))
      if(collcom%ndimind_f>0) call to_zero(7*collcom%ndimind_f, psit_f(1))
  end if


  !SM: check if the modulo operations take a lot of time. If so, try to use an
  !auxiliary array with shifted bounds in order to access smat%matrixindex_in_compressed_fortransposed

  spin_loop: do ispin=1,sparsemat%nspin

      ishift_mat=(ispin-1)*sparsemat%nvctr
      iorb_shift=(ispin-1)*sparsemat%nfvctr

      !$omp parallel default(private) &
      !$omp shared(collcom, psit_c, psitwork_c, psit_f, psitwork_f, sparsemat, mat, ispin, ishift_mat, iorb_shift)

      !$omp do schedule(static,1)
       do ipt=1,collcom%nptsp_c 
          ii=collcom%norb_per_gridpoint_c(ipt) 
          i0 = collcom%isptsp_c(ipt) + (ispin-1)*collcom%ndimind_c/sparsemat%nspin
          do i=1,ii
              i0i=i0+i
              iiorb=collcom%indexrecvorbital_c(i0i) - iorb_shift
              !iiorb=mod(iiorb-1,sparsemat%nfvctr)+1
              m=mod(ii,7)
              tt0=0.d0 ; tt1=0.d0 ; tt2=0.d0 ; tt3=0.d0 ; tt4=0.d0 ; tt5=0.d0 ; tt6=0.d0
              if(m/=0) then
                  do j=1,m
                      i0j=i0+j
                      jjorb=collcom%indexrecvorbital_c(i0j) - iorb_shift
                      !jjorb=mod(jjorb-1,sparsemat%nfvctr)+1
                      ind0 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
                      ind0=ind0+ishift_mat
                      tt0=tt0+mat%matrix_compr(ind0)*psitwork_c(i0j)
                  end do
              end if
              do j=m+1,ii,7
                  i0j=i0+j

                  jjorb0=collcom%indexrecvorbital_c(i0j+0) - iorb_shift
                  !jjorb0=mod(jjorb0-1,sparsemat%nfvctr)+1
                  ind0 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb0,iiorb)
                  ind0=ind0+ishift_mat
                  tt0=tt0+mat%matrix_compr(ind0)*psitwork_c(i0j+0)

                  jjorb1=collcom%indexrecvorbital_c(i0j+1) - iorb_shift
                  !jjorb1=mod(jjorb1-1,sparsemat%nfvctr)+1
                  ind1 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb1,iiorb)
                  ind1=ind1+ishift_mat
                  tt1=tt1+mat%matrix_compr(ind1)*psitwork_c(i0j+1)

                  jjorb2=collcom%indexrecvorbital_c(i0j+2) - iorb_shift
                  !jjorb2=mod(jjorb2-1,sparsemat%nfvctr)+1
                  ind2 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb2,iiorb)
                  ind2=ind2+ishift_mat
                  tt2=tt2+mat%matrix_compr(ind2)*psitwork_c(i0j+2)

                  jjorb3=collcom%indexrecvorbital_c(i0j+3) - iorb_shift
                  !jjorb3=mod(jjorb3-1,sparsemat%nfvctr)+1
                  ind3 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb3,iiorb)
                  ind3=ind3+ishift_mat
                  tt3=tt3+mat%matrix_compr(ind3)*psitwork_c(i0j+3)

                  jjorb4=collcom%indexrecvorbital_c(i0j+4) - iorb_shift
                  !jjorb4=mod(jjorb4-1,sparsemat%nfvctr)+1
                  ind4 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb4,iiorb)
                  ind4=ind4+ishift_mat
                  tt4=tt4+mat%matrix_compr(ind4)*psitwork_c(i0j+4)

                  jjorb5=collcom%indexrecvorbital_c(i0j+5) - iorb_shift
                  !jjorb5=mod(jjorb5-1,sparsemat%nfvctr)+1
                  ind5 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb5,iiorb)
                  ind5=ind5+ishift_mat
                  tt5=tt5+mat%matrix_compr(ind5)*psitwork_c(i0j+5)

                  jjorb6=collcom%indexrecvorbital_c(i0j+6) - iorb_shift
                  !jjorb6=mod(jjorb6-1,sparsemat%nfvctr)+1
                  ind6 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb6,iiorb)
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
                      !jjorb=mod(jjorb-1,sparsemat%nfvctr)+1
                      ind0 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
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
                  !jjorb0=mod(jjorb0-1,sparsemat%nfvctr)+1
                  ind0 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb0,iiorb)
                  ind0=ind0+ishift_mat
                  tt06 = tt06 + mat%matrix_compr(ind0)*psitwork_f(i07j-6)
                  tt05 = tt05 + mat%matrix_compr(ind0)*psitwork_f(i07j-5)
                  tt04 = tt04 + mat%matrix_compr(ind0)*psitwork_f(i07j-4)
                  tt03 = tt03 + mat%matrix_compr(ind0)*psitwork_f(i07j-3)
                  tt02 = tt02 + mat%matrix_compr(ind0)*psitwork_f(i07j-2)
                  tt01 = tt01 + mat%matrix_compr(ind0)*psitwork_f(i07j-1)
                  tt00 = tt00 + mat%matrix_compr(ind0)*psitwork_f(i07j-0)

                  jjorb1=collcom%indexrecvorbital_f(i0j+1) - iorb_shift
                  !jjorb1=mod(jjorb1-1,sparsemat%nfvctr)+1
                  ind1 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb1,iiorb)
                  ind1=ind1+ishift_mat
                  tt16 = tt16 + mat%matrix_compr(ind1)*psitwork_f(i07j+1) !+1*7-6
                  tt15 = tt15 + mat%matrix_compr(ind1)*psitwork_f(i07j+2) !+1*7-5
                  tt14 = tt14 + mat%matrix_compr(ind1)*psitwork_f(i07j+3) !+1*7-4
                  tt13 = tt13 + mat%matrix_compr(ind1)*psitwork_f(i07j+4) !+1*7-3
                  tt12 = tt12 + mat%matrix_compr(ind1)*psitwork_f(i07j+5) !+1*7-2
                  tt11 = tt11 + mat%matrix_compr(ind1)*psitwork_f(i07j+6) !+1*7-1
                  tt10 = tt10 + mat%matrix_compr(ind1)*psitwork_f(i07j+7) !+1*7-0

                  jjorb2=collcom%indexrecvorbital_f(i0j+2) - iorb_shift
                  !jjorb2=mod(jjorb2-1,sparsemat%nfvctr)+1
                  ind2 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb2,iiorb)
                  ind2=ind2+ishift_mat
                  tt26 = tt26 + mat%matrix_compr(ind2)*psitwork_f(i07j+8) !+2*7-6
                  tt25 = tt25 + mat%matrix_compr(ind2)*psitwork_f(i07j+9) !+2*7-5
                  tt24 = tt24 + mat%matrix_compr(ind2)*psitwork_f(i07j+10) !+2*7-4
                  tt23 = tt23 + mat%matrix_compr(ind2)*psitwork_f(i07j+11) !+2*7-3
                  tt22 = tt22 + mat%matrix_compr(ind2)*psitwork_f(i07j+12) !+2*7-2
                  tt21 = tt21 + mat%matrix_compr(ind2)*psitwork_f(i07j+13) !+2*7-1
                  tt20 = tt20 + mat%matrix_compr(ind2)*psitwork_f(i07j+14) !+2*7-0

                  jjorb3=collcom%indexrecvorbital_f(i0j+3) - iorb_shift
                  !jjorb3=mod(jjorb3-1,sparsemat%nfvctr)+1
                  ind3 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb3,iiorb)
                  ind3=ind3+ishift_mat
                  tt36 = tt36 + mat%matrix_compr(ind3)*psitwork_f(i07j+15) !+3*7-6
                  tt35 = tt35 + mat%matrix_compr(ind3)*psitwork_f(i07j+16) !+3*7-5
                  tt34 = tt34 + mat%matrix_compr(ind3)*psitwork_f(i07j+17) !+3*7-4
                  tt33 = tt33 + mat%matrix_compr(ind3)*psitwork_f(i07j+18) !+3*7-3
                  tt32 = tt32 + mat%matrix_compr(ind3)*psitwork_f(i07j+19) !+3*7-2
                  tt31 = tt31 + mat%matrix_compr(ind3)*psitwork_f(i07j+20) !+3*7-1
                  tt30 = tt30 + mat%matrix_compr(ind3)*psitwork_f(i07j+21) !+3*7-0

                  jjorb4=collcom%indexrecvorbital_f(i0j+4) - iorb_shift
                  !jjorb4=mod(jjorb4-1,sparsemat%nfvctr)+1
                  ind4 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb4,iiorb)
                  ind4=ind4+ishift_mat
                  tt46 = tt46 + mat%matrix_compr(ind4)*psitwork_f(i07j+22) !+4*7-6
                  tt45 = tt45 + mat%matrix_compr(ind4)*psitwork_f(i07j+23) !+4*7-5
                  tt44 = tt44 + mat%matrix_compr(ind4)*psitwork_f(i07j+24) !+4*7-4
                  tt43 = tt43 + mat%matrix_compr(ind4)*psitwork_f(i07j+25) !+4*7-3
                  tt42 = tt42 + mat%matrix_compr(ind4)*psitwork_f(i07j+26) !+4*7-2
                  tt41 = tt41 + mat%matrix_compr(ind4)*psitwork_f(i07j+27) !+4*7-1
                  tt40 = tt40 + mat%matrix_compr(ind4)*psitwork_f(i07j+28) !+4*7-0

                  jjorb5=collcom%indexrecvorbital_f(i0j+5) - iorb_shift
                  !jjorb5=mod(jjorb5-1,sparsemat%nfvctr)+1
                  ind5 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb5,iiorb)
                  ind5=ind5+ishift_mat
                  tt56 = tt56 + mat%matrix_compr(ind5)*psitwork_f(i07j+29) !+5*7-6
                  tt55 = tt55 + mat%matrix_compr(ind5)*psitwork_f(i07j+30) !+5*7-5
                  tt54 = tt54 + mat%matrix_compr(ind5)*psitwork_f(i07j+31) !+5*7-4
                  tt53 = tt53 + mat%matrix_compr(ind5)*psitwork_f(i07j+32) !+5*7-3
                  tt52 = tt52 + mat%matrix_compr(ind5)*psitwork_f(i07j+33) !+5*7-2
                  tt51 = tt51 + mat%matrix_compr(ind5)*psitwork_f(i07j+34) !+5*7-1
                  tt50 = tt50 + mat%matrix_compr(ind5)*psitwork_f(i07j+35) !+5*7-0

                  jjorb6=collcom%indexrecvorbital_f(i0j+6) - iorb_shift
                  !jjorb6=mod(jjorb6-1,sparsemat%nfvctr)+1
                  ind6 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb6,iiorb)
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

  call f_release_routine()
  call timing(iproc,'lincombtrans  ','OF') !lr408t

end subroutine build_linear_combination_transposed




subroutine check_grid_point_from_boxes(i1, i2, i3, lr, overlap_possible)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: i1, i2, i3
  type(locreg_descriptors),intent(in) :: lr  
  logical,intent(out) :: overlap_possible

  ! Local variables
  logical :: ovrlpx, ovrlpy, ovrlpz
  
  ovrlpx = (i1>=lr%ns1 .and. i1<=lr%ns1+lr%d%n1)
  ovrlpy = (i2>=lr%ns2 .and. i2<=lr%ns2+lr%d%n2)
  ovrlpz = (i3>=lr%ns3 .and. i3<=lr%ns3+lr%d%n3)
  if(ovrlpx .and. ovrlpy .and. ovrlpz) then
      overlap_possible=.true.
  else
      overlap_possible=.true.
  end if

end subroutine check_grid_point_from_boxes


!!subroutine get_reverse_indices(n, indices, reverse_indices)
!!  use module_base
!!  implicit none
!!  
!!  ! Calling arguments
!!  integer,intent(in) :: n
!!  integer,dimension(n),intent(in) :: indices
!!  integer,dimension(n),intent(out) :: reverse_indices
!!
!!  ! Local variables
!!  integer :: i, j, m, j0, j1, j2, j3
!!
!!  !$omp parallel default(private) &
!!  !$omp shared(n, m, indices, reverse_indices)
!!
!!  m=mod(n,4)
!!  if (m/=0) then
!!      do i=1,m
!!          j=indices(i)
!!          reverse_indices(j)=i
!!      end do
!!  end if
!!
!!  !$omp do
!!  do i=m+1,n,4
!!      j0=indices(i+0)
!!      reverse_indices(j0)=i+0
!!      j1=indices(i+1)
!!      reverse_indices(j1)=i+1
!!      j2=indices(i+2)
!!      reverse_indices(j2)=i+2
!!      j3=indices(i+3)
!!      reverse_indices(j3)=i+3
!!  end do
!!  !$omp end do
!!
!!  !$omp end parallel
!!
!!  !!do i=1,n
!!  !!    j=indices(i)
!!  !!    reverse_indices(j)=i
!!  !!end do
!!
!!end subroutine get_reverse_indices


subroutine normalize_transposed(iproc, nproc, orbs, nspin, collcom, psit_c, psit_f, norm)
  use module_base
  use module_types
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

  call to_zero(orbs%norb, norm(1))


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
      call mpiallred(norm(1), orbs%norb, mpi_sum, bigdft_mpi%mpi_comm)
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



subroutine init_matrixindex_in_compressed_fortransposed(iproc, nproc, orbs, collcom, collcom_shamop, &
           collcom_sr, sparsemat)
  use module_base
  use module_types
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
  integer :: iorb, jorb, istat, imin, imax
  !integer :: kproc,jproc,jjorbold,jjorb,isend,irecv,ilr,ijorb,iiorb,ind,ierr, irow, irowold, iseg
  !integer :: compressed_index
!  integer,dimension(:,:),allocatable :: sendbuf, requests, iminmaxarr
  character(len=*),parameter :: subname='init_sparse_matrix'


  ! for the calculation of overlaps and the charge density
  imin=minval(collcom%indexrecvorbital_c)
  imin=min(imin,minval(collcom%indexrecvorbital_f))
  imin=min(imin,minval(collcom_shamop%indexrecvorbital_c))
  imin=min(imin,minval(collcom_shamop%indexrecvorbital_f))
  imin=min(imin,minval(collcom_sr%indexrecvorbital_c))
  imax=maxval(collcom%indexrecvorbital_c)
  imax=max(imax,maxval(collcom%indexrecvorbital_f))
  imax=max(imax,maxval(collcom_shamop%indexrecvorbital_c))
  imax=max(imax,maxval(collcom_shamop%indexrecvorbital_f))
  imax=max(imax,maxval(collcom_sr%indexrecvorbital_c))

  ! values regardless of the spin
  imin=mod(imin-1,sparsemat%nfvctr)+1
  imax=mod(imax-1,sparsemat%nfvctr)+1
  

  !!! This is a temporary solution for spin polarized systems
  !!imax=min(imax,orbs%norbu)



  !!allocate(sparsemat%matrixindex_in_compressed_fortransposed(imin:imax,imin:imax), stat=istat)
  !!call memocc(istat, sparsemat%matrixindex_in_compressed_fortransposed, &
  sparsemat%matrixindex_in_compressed_fortransposed=f_malloc_ptr((/imin.to.imax,imin.to.imax/),&
      id='sparsemat%matrixindex_in_compressed_fortransposed')

  !$omp parallel do default(private) shared(sparsemat,orbs,imin,imax)
  do iorb=imin,imax
      do jorb=imin,imax
          !@ii=(jorb-1)*sparsemat%nfvctr+iorb
          !@ispin=(ii-1)/sparsemat%nfvctr+1 !integer division to get the spin (1 for spin up (or non polarized), 2 for spin down)
          !@iiorb=mod(iorb-1,sparsemat%nfvctr)+1 !orbital number regardless of the spin
          !@jjorb=mod(jorb-1,sparsemat%nfvctr)+1 !orbital number regardless of the spin
          !sparsemat%matrixindex_in_compressed_fortransposed(iorb,jorb)=compressed_index(iiorb,jjorb,orbs%norbu,sparsemat)
          sparsemat%matrixindex_in_compressed_fortransposed(iorb,jorb)=matrixindex_in_compressed(sparsemat, iorb, jorb)
          !sendbuf(jorb,iorb)=compressed_index(jorb,iiorb,orbs%norb,sparsemat)
          !sendbuf(iorb,jorb)=compressed_index(iiorb,jorb,orbs%norb,sparsemat)
      end do
  end do
  !$omp end parallel do

  !@! Add the spin shift (i.e. the index is in the spin polarized matrix which is at the end)
  !@if (ispin==2) then
  !@    matrixindex_in_compressed = matrixindex_in_compressed + sparsemat%nvctr
  !@end if


end subroutine init_matrixindex_in_compressed_fortransposed
