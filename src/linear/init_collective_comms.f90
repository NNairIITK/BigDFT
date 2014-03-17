!> @file
!! Intialization of the collective communications for the linear version
!! @author
!!    Copyright (C) 2011-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS



subroutine check_communications_locreg(iproc,nproc,orbs,Lzd,collcom,npsidim_orbs,npsidim_comp)
   use module_base
   use module_types
   use yaml_output
   use communications, only: transpose_localized, untranspose_localized
   implicit none
   integer, intent(in) :: iproc,nproc
   type(orbitals_data), intent(in) :: orbs
   type(local_zone_descriptors), intent(in) :: lzd
   type(comms_linear), intent(in) :: collcom
   integer, intent(in) :: npsidim_orbs, npsidim_comp
   !local variables
   character(len=*), parameter :: subname='check_communications'
   integer, parameter :: ilog=6
   integer :: i,ispinor,iorb,indspin,i_stat,i_all,ikptsp
   integer :: ikpt,ierr,i0,ifine,ii,iiorb,ipt,jorb
   integer :: icomp
   !!$integer :: ipsi,ipsic,ipsif,ipsiworkc,ipsiworkf,jcomp,jkpt
   real(wp) :: psival,maxdiff,tt
   real(wp), dimension(:), allocatable :: psi,psit_c,psit_f
   real(wp), dimension(:,:), allocatable :: checksum
   real(wp) :: epsilon,tol
   logical :: abort

   !allocate the "wavefunction" and fill it, and also the workspace
   allocate(psi(max(npsidim_orbs,npsidim_comp)+ndebug),stat=i_stat)
   call memocc(i_stat,psi,'psi',subname)
   allocate(psit_c(sum(collcom%nrecvcounts_c)), stat=i_stat)
   call memocc(i_stat, psit_c, 'psit_c', subname)
   allocate(psit_f(7*sum(collcom%nrecvcounts_f)), stat=i_stat)
   call memocc(i_stat, psit_f, 'psit_f', subname)
   allocate(checksum(orbs%norb*orbs%nspinor,2), stat=i_stat)
   call memocc(i_stat, checksum, 'checksum', subname)
   if (orbs%norbp>0) then
      tol=1.e-10*real(npsidim_orbs,wp)/real(orbs%norbp,wp)
   else
      tol=0.0_wp
   end if

   checksum(:,:)=0.0_wp
   do iorb=1,orbs%norbp
      ikpt=(orbs%isorb+iorb-1)/orbs%norb+1
      do ispinor=1,orbs%nspinor
         indspin=(ispinor-1)*nvctr_orb(iorb)
         checksum(orbs%isorb+iorb+(ispinor-1)*orbs%nspinor,1)=0.0_wp
         do i=1,nvctr_orb(iorb)
            !vali=real(i,wp)/512.0_wp  ! *1.d-5
            call test_value_locreg(ikpt,orbs%isorb+iorb-(ikpt-1)*orbs%norb,ispinor,i,psival)
            psi(i+indspin+ind_orb(iorb))=psival!(valorb+vali)*(-1)**(ispinor-1)
            checksum(orbs%isorb+iorb+(ispinor-1)*orbs%nspinor,1)=&
                 checksum(orbs%isorb+iorb+(ispinor-1)*orbs%nspinor,1)+psival
         end do
      end do
   end do

   call transpose_localized(iproc, nproc, npsidim_orbs, orbs, collcom, psi, psit_c, psit_f, lzd)
   
   !check the results of the transposed wavefunction
   maxdiff=0.0_wp
   if (iproc==0) call yaml_map('Number of coarse and fine DoF (MasterMPI task)',&
        (/collcom%nptsp_c,collcom%nptsp_f/),fmt='(i8)')

   do ikptsp=1,1!orbs%nkptsp !should be one for the moment
      ikpt=orbs%iskpts+ikptsp!orbs%ikptsp(ikptsp)
      ispinor=1 !for the (long?) moment
      icomp=1
      if (collcom%nptsp_c>0) then
         do ipt=1,collcom%nptsp_c 
            ii=collcom%norb_per_gridpoint_c(ipt)
            i0 = collcom%isptsp_c(ipt) 
            do i=1,ii
               iiorb=collcom%indexrecvorbital_c(i0+i)
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
               icomp=icomp+1
            end do
         end do
      end if
      icomp=1
      if (collcom%nptsp_f>0) then
         do ipt=1,collcom%nptsp_f 
            ii=collcom%norb_per_gridpoint_f(ipt) 
            i0 = collcom%isptsp_f(ipt) 
            do i=1,ii
               iiorb=collcom%indexrecvorbital_f(i0+i)
!!$               ipsitworkf=collcom%iexpand_f(icomp)
!!$               ipsif=collcom%isendbuf_f(ipsiworkf)
!!$               ipsi=ipsif
!!$               do jorb=1,iiorb-1
!!$                  ipsi=ipsi-nvctr_f_orb(jorb)
!!$               end do

               do ifine=1,7
!!$                  call test_value_locreg(ikpt,iiorb-(ikpt-1)*orbs%norb,ispinor,&
!!$                       nvctr_c_orb(iiorb)+7*(ipsi-1)+ifine,psival) 
!!$                  tt=abs(psit_f(7*(i0+i-1)+ifine)-psival)
!!$                  if (tt > maxdiff) then
!!$                     maxdiff=tt
!!$                     !call wrong_components(psival,jkpt,jorb,jcomp)
!!$                  end if
                  checksum(iiorb,2)=checksum(iiorb,2)+psit_f(7*(i0+i-1)+ifine)
               end do
               icomp=icomp+1
            end do
         end do
      end if
   end do
!!$
   if (iproc==0) call yaml_map('Tolerances for this check',&
        (/tol,real(orbs%norb,wp)*epsilon(1.0_wp)/),fmt='(1pe25.17)')

   call MPI_BARRIER(bigdft_mpi%mpi_comm, ierr)
   call mpiallred(checksum(1,1),2*orbs%norb*orbs%nspinor,MPI_SUM,&
        bigdft_mpi%mpi_comm,ierr)
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

   if (abort) call MPI_ABORT(bigdft_mpi%mpi_comm,ierr)


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
   call MPI_BARRIER(bigdft_mpi%mpi_comm, ierr)

   call mpiallred(maxdiff,1,MPI_MAX,bigdft_mpi%mpi_comm,ierr)
   if (iproc==0) call yaml_map('Maxdiff for untranspose',maxdiff,fmt='(1pe25.17)')

   i_all=-product(shape(psi))*kind(psi)
   deallocate(psi,stat=i_stat)
   call memocc(i_stat,i_all,'psi',subname)
   i_all=-product(shape(psit_c))*kind(psit_c)
   deallocate(psit_c, stat=i_stat)
   call memocc(i_stat, i_all, 'psit_c', subname)
   i_all=-product(shape(psit_f))*kind(psit_f)
   deallocate(psit_f, stat=i_stat)
   call memocc(i_stat, i_all, 'psit_f', subname)
   i_all=-product(shape(checksum))*kind(checksum)
   deallocate(checksum, stat=i_stat)
   call memocc(i_stat, i_all, 'checksum', subname)



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
           psit_c1, psit_c2, psit_f1, psit_f2, ovrlp)
  use module_base
  use module_types
  use sparsematrix_base, only: sparseMatrix
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(orbitals_data),intent(in) :: orbs
  type(comms_linear),intent(in) :: collcom
  real(kind=8),dimension(collcom%ndimind_c),intent(in) :: psit_c1, psit_c2
  real(kind=8),dimension(7*collcom%ndimind_f),intent(in) :: psit_f1, psit_f2
  type(sparse_matrix),intent(inout) :: ovrlp

  ! Local variables
  integer :: i0, ipt, ii, iiorb, j, jjorb, i, ierr, istat, m, tid, norb, nthreads
  integer :: istart, iend, orb_rest, ind0, ind1, ind2, ind3
  integer,dimension(:),allocatable :: n
  !$ integer  :: omp_get_thread_num,omp_get_max_threads

  call timing(iproc,'ovrlptransComp','ON') !lr408t

  call to_zero(ovrlp%nvctr, ovrlp%matrix_compr(1))

  nthreads=1
  !$  nthreads = OMP_GET_max_threads()

  allocate(n(nthreads),stat=istat)

  norb = orbs%norb / nthreads
  
  orb_rest = orbs%norb-norb*nthreads

  n = norb

  do i = 1,orb_rest
     n(i) = n(i) + 1
  end do

  !$omp parallel default(private) &
  !$omp shared(collcom, ovrlp, psit_c1, psit_c2, psit_f1, psit_f2, n)
  tid=0
  !$ tid = OMP_GET_THREAD_NUM()
  istart = 1
  iend = 0

  do i = 0,tid-1
     istart = istart + n(i+1) 
  end do

  do i = 0,tid
     iend = iend + n(i+1)
  end do

 ! write(*,*) '#########################################',tid,istart,iend

  if (collcom%nptsp_c>0) then
  
      do ipt=1,collcom%nptsp_c 
          ii=collcom%norb_per_gridpoint_c(ipt) 
          i0 = collcom%isptsp_c(ipt)

          do i=1,ii
              iiorb=collcom%indexrecvorbital_c(i0+i)
              if(iiorb < istart .or. iiorb > iend) cycle

              m=mod(ii,4)
              if(m/=0) then
                  do j=1,m
                      jjorb=collcom%indexrecvorbital_c(i0+j)
                      ind0 = ovrlp%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
                      !ind0 = collcom%matrixindex_in_compressed(jjorb,iiorb)
                      ovrlp%matrix_compr(ind0) = ovrlp%matrix_compr(ind0) + psit_c1(i0+i)*psit_c2(i0+j)
                  end do
              end if
              do j=m+1,ii,4

                  jjorb=collcom%indexrecvorbital_c(i0+j+0)
                  ind0 = ovrlp%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
                  !ind0 = collcom%matrixindex_in_compressed(jjorb,iiorb)
                  ovrlp%matrix_compr(ind0) = ovrlp%matrix_compr(ind0) + psit_c1(i0+i)*psit_c2(i0+j+0)

                  jjorb=collcom%indexrecvorbital_c(i0+j+1)
                  ind1 = ovrlp%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
                  !ind1 = collcom%matrixindex_in_compressed(jjorb,iiorb)
                  ovrlp%matrix_compr(ind1) = ovrlp%matrix_compr(ind1) + psit_c1(i0+i)*psit_c2(i0+j+1)

                  jjorb=collcom%indexrecvorbital_c(i0+j+2)
                  ind2 = ovrlp%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
                  !ind2 = collcom%matrixindex_in_compressed(jjorb,iiorb)
                  ovrlp%matrix_compr(ind2) = ovrlp%matrix_compr(ind2) + psit_c1(i0+i)*psit_c2(i0+j+2)

                  jjorb=collcom%indexrecvorbital_c(i0+j+3)
                  ind3 = ovrlp%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
                  !ind3 = collcom%matrixindex_in_compressed(jjorb,iiorb)
                  ovrlp%matrix_compr(ind3) = ovrlp%matrix_compr(ind3) + psit_c1(i0+i)*psit_c2(i0+j+3)

              end do
          end do
      end do

  end if
  

  if (collcom%nptsp_f>0) then
   
      do ipt=1,collcom%nptsp_f 
          ii=collcom%norb_per_gridpoint_f(ipt) 
          i0 = collcom%isptsp_f(ipt)

          do i=1,ii

              iiorb=collcom%indexrecvorbital_f(i0+i)
              if(iiorb < istart .or. iiorb > iend) cycle

              do j=1,ii
                  jjorb=collcom%indexrecvorbital_f(i0+j)
                  ind0 = ovrlp%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
                  !ind0 = collcom%matrixindex_in_compressed(jjorb,iiorb)
                  ovrlp%matrix_compr(ind0) = ovrlp%matrix_compr(ind0) + psit_f1(7*(i0+i)-6)*psit_f2(7*(i0+j)-6)
                  ovrlp%matrix_compr(ind0) = ovrlp%matrix_compr(ind0) + psit_f1(7*(i0+i)-5)*psit_f2(7*(i0+j)-5)
                  ovrlp%matrix_compr(ind0) = ovrlp%matrix_compr(ind0) + psit_f1(7*(i0+i)-4)*psit_f2(7*(i0+j)-4)
                  ovrlp%matrix_compr(ind0) = ovrlp%matrix_compr(ind0) + psit_f1(7*(i0+i)-3)*psit_f2(7*(i0+j)-3)
                  ovrlp%matrix_compr(ind0) = ovrlp%matrix_compr(ind0) + psit_f1(7*(i0+i)-2)*psit_f2(7*(i0+j)-2)
                  ovrlp%matrix_compr(ind0) = ovrlp%matrix_compr(ind0) + psit_f1(7*(i0+i)-1)*psit_f2(7*(i0+j)-1)
                  ovrlp%matrix_compr(ind0) = ovrlp%matrix_compr(ind0) + psit_f1(7*(i0+i)-0)*psit_f2(7*(i0+j)-0)
              end do
          end do
      end do
  end if

  !$omp end parallel

  call timing(iproc,'ovrlptransComp','OF') !lr408t

  call timing(iproc,'ovrlptransComm','ON') !lr408t

  if(nproc>1) then
      call mpiallred(ovrlp%matrix_compr(1), ovrlp%nvctr, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  end if

  ovrlp%can_use_dense=.false.
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
      call mpiallred(ovrlp(1,1), orbs1%norb*orbs2%norb, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  end if
  call timing(iproc,'ovrlptransComm','OF') !lr408t
end subroutine calculate_pulay_overlap

subroutine build_linear_combination_transposed(collcom, sparsemat, psitwork_c, psitwork_f, &
     reset, psit_c, psit_f, iproc)
  use module_base
  use module_types
  use sparsematrix_base, only: sparseMatrix
  implicit none
  
  ! Calling arguments
  type(sparse_matrix),intent(in) :: sparsemat
  type(comms_linear),intent(in) :: collcom
  real(kind=8),dimension(collcom%ndimind_c),intent(in) :: psitwork_c
  real(kind=8),dimension(7*collcom%ndimind_f),intent(in) :: psitwork_f
  logical,intent(in) :: reset
  real(kind=8),dimension(collcom%ndimind_c),intent(inout) :: psit_c
  real(kind=8),dimension(7*collcom%ndimind_f),intent(inout) :: psit_f
  integer, intent(in) :: iproc
  ! Local variables
  integer :: i0, ipt, ii, j, iiorb, jjorb, i, m, ind0, ind1, ind2, ind3

  call timing(iproc,'lincombtrans  ','ON') !lr408t
  if(reset) then
      if(collcom%ndimind_c>0) call to_zero(collcom%ndimind_c, psit_c(1))
      if(collcom%ndimind_f>0) call to_zero(7*collcom%ndimind_f, psit_f(1))
  end if

 
  !$omp parallel default(private) &
  !$omp shared(collcom, psit_c, psitwork_c, psit_f, psitwork_f, sparsemat)

  !!write(*,'(a,i4,4i8)') 'iproc, lbound, ubound, minval, maxval',&
  !!iproc, lbound(sparsemat%matrixindex_in_compressed_fortransposed,2),&
  !!ubound(sparsemat%matrixindex_in_compressed_fortransposed,2),&
  !!minval(collcom%indexrecvorbital_c),maxval(collcom%indexrecvorbital_c)

  !$omp do
   do ipt=1,collcom%nptsp_c 
      ii=collcom%norb_per_gridpoint_c(ipt) 
      i0 = collcom%isptsp_c(ipt)
      do i=1,ii
          iiorb=collcom%indexrecvorbital_c(i0+i)
          m=mod(ii,4)
          if(m/=0) then
              do j=1,m
                  jjorb=collcom%indexrecvorbital_c(i0+j)
                  ind0 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
                  !ind0 = collcom%matrixindex_in_compressed(jjorb,iiorb)
                  !write(41,*) jjorb, iiorb, sparsemat%matrixindex_in_compressed(jjorb,iiorb)
                  !write(42,*) jjorb, iiorb, collcom%matrixindex_in_compressed(jjorb,iiorb)
                  psit_c(i0+i)=psit_c(i0+i)+sparsemat%matrix_compr(ind0)*psitwork_c(i0+j)
              end do
          end if
          do j=m+1,ii,4
              jjorb=collcom%indexrecvorbital_c(i0+j+0)
              ind0 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
              !ind0 = collcom%matrixindex_in_compressed(jjorb,iiorb)
              psit_c(i0+i)=psit_c(i0+i)+sparsemat%matrix_compr(ind0)*psitwork_c(i0+j+0)

              jjorb=collcom%indexrecvorbital_c(i0+j+1)
              ind1 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
              !ind1 = collcom%matrixindex_in_compressed(jjorb,iiorb)
              psit_c(i0+i)=psit_c(i0+i)+sparsemat%matrix_compr(ind1)*psitwork_c(i0+j+1)

              jjorb=collcom%indexrecvorbital_c(i0+j+2)
              ind2 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
              !ind2 = collcom%matrixindex_in_compressed(jjorb,iiorb)
              psit_c(i0+i)=psit_c(i0+i)+sparsemat%matrix_compr(ind2)*psitwork_c(i0+j+2)

              jjorb=collcom%indexrecvorbital_c(i0+j+3)
              ind3 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
              !ind3 = collcom%matrixindex_in_compressed(jjorb,iiorb)
              psit_c(i0+i)=psit_c(i0+i)+sparsemat%matrix_compr(ind3)*psitwork_c(i0+j+3)

          end do
      end do
  end do
  !$omp end do

  !$omp do
  do ipt=1,collcom%nptsp_f 
      ii=collcom%norb_per_gridpoint_f(ipt) 
       i0 = collcom%isptsp_f(ipt)
      do i=1,ii
          iiorb=collcom%indexrecvorbital_f(i0+i)
          do j=1,ii
              jjorb=collcom%indexrecvorbital_f(i0+j)
              ind0 = sparsemat%matrixindex_in_compressed_fortransposed(jjorb,iiorb)
              !ind0 = collcom%matrixindex_in_compressed(jjorb,iiorb)
              psit_f(7*(i0+i)-6) = psit_f(7*(i0+i)-6) + sparsemat%matrix_compr(ind0)*psitwork_f(7*(i0+j)-6)
              psit_f(7*(i0+i)-5) = psit_f(7*(i0+i)-5) + sparsemat%matrix_compr(ind0)*psitwork_f(7*(i0+j)-5)
              psit_f(7*(i0+i)-4) = psit_f(7*(i0+i)-4) + sparsemat%matrix_compr(ind0)*psitwork_f(7*(i0+j)-4)
              psit_f(7*(i0+i)-3) = psit_f(7*(i0+i)-3) + sparsemat%matrix_compr(ind0)*psitwork_f(7*(i0+j)-3)
              psit_f(7*(i0+i)-2) = psit_f(7*(i0+i)-2) + sparsemat%matrix_compr(ind0)*psitwork_f(7*(i0+j)-2)
              psit_f(7*(i0+i)-1) = psit_f(7*(i0+i)-1) + sparsemat%matrix_compr(ind0)*psitwork_f(7*(i0+j)-1)
              psit_f(7*(i0+i)-0) = psit_f(7*(i0+i)-0) + sparsemat%matrix_compr(ind0)*psitwork_f(7*(i0+j)-0)
          end do
      end do
     
  end do
  !$omp end do
  !$omp end parallel

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


subroutine normalize_transposed(iproc, nproc, orbs, collcom, psit_c, psit_f, norm)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(orbitals_data),intent(in):: orbs
  type(comms_linear),intent(in):: collcom
  real(8),dimension(collcom%ndimind_c),intent(inout):: psit_c
  real(8),dimension(7*collcom%ndimind_f),intent(inout):: psit_f
  real(8),dimension(orbs%norb),intent(out):: norm
  
  ! Local variables
  integer:: i0, ipt, ii, iiorb, i, ierr, iorb

  call timing(iproc,'norm_trans','ON')

  call to_zero(orbs%norb, norm(1))

  !$omp parallel default(private) &
  !$omp shared(collcom, norm, psit_c,psit_f,orbs)

  if (collcom%nptsp_c>0) then
      !$omp do reduction(+:norm)
      do ipt=1,collcom%nptsp_c 
          ii=collcom%norb_per_gridpoint_c(ipt)
          i0 = collcom%isptsp_c(ipt) 
          do i=1,ii
              iiorb=collcom%indexrecvorbital_c(i0+i)
              norm(iiorb)=norm(iiorb)+psit_c(i0+i)**2
          end do
      end do
      !$omp end do
  end if

  if (collcom%nptsp_f>0) then
      !$omp do reduction(+:norm)
      do ipt=1,collcom%nptsp_f 
          ii=collcom%norb_per_gridpoint_f(ipt) 
          i0 = collcom%isptsp_f(ipt) 
          do i=1,ii
              iiorb=collcom%indexrecvorbital_f(i0+i)
              norm(iiorb)=norm(iiorb)+psit_f(7*(i0+i)-6)**2
              norm(iiorb)=norm(iiorb)+psit_f(7*(i0+i)-5)**2
              norm(iiorb)=norm(iiorb)+psit_f(7*(i0+i)-4)**2
              norm(iiorb)=norm(iiorb)+psit_f(7*(i0+i)-3)**2
              norm(iiorb)=norm(iiorb)+psit_f(7*(i0+i)-2)**2
              norm(iiorb)=norm(iiorb)+psit_f(7*(i0+i)-1)**2
              norm(iiorb)=norm(iiorb)+psit_f(7*(i0+i)-0)**2
          end do
      end do
      !$omp end do
  end if

  !$omp end parallel
  
  if(nproc>1) then
      call mpiallred(norm(1), orbs%norb, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  end if

  do iorb=1,orbs%norb
     norm(iorb)=1.d0/sqrt(norm(iorb))
  end do

  !$omp parallel default(private) shared(norm,orbs,collcom,psit_c,psit_f)  
  !$omp do
  do ipt=1,collcom%nptsp_c 
      ii=collcom%norb_per_gridpoint_c(ipt)
      i0=collcom%isptsp_c(ipt)
      do i=1,ii
          iiorb=collcom%indexrecvorbital_c(i0+i)
          psit_c(i0+i)=psit_c(i0+i)*norm(iiorb)
      end do
    
  end do
  !$omp end do
  !$omp do
  do ipt=1,collcom%nptsp_f 
      ii=collcom%norb_per_gridpoint_f(ipt)
      i0 = collcom%isptsp_f(ipt) 
      do i=1,ii
          iiorb=collcom%indexrecvorbital_f(i0+i)
          psit_f(7*(i0+i)-6)=psit_f(7*(i0+i)-6)*norm(iiorb)
          psit_f(7*(i0+i)-5)=psit_f(7*(i0+i)-5)*norm(iiorb)
          psit_f(7*(i0+i)-4)=psit_f(7*(i0+i)-4)*norm(iiorb)
          psit_f(7*(i0+i)-3)=psit_f(7*(i0+i)-3)*norm(iiorb)
          psit_f(7*(i0+i)-2)=psit_f(7*(i0+i)-2)*norm(iiorb)
          psit_f(7*(i0+i)-1)=psit_f(7*(i0+i)-1)*norm(iiorb)
          psit_f(7*(i0+i)-0)=psit_f(7*(i0+i)-0)*norm(iiorb)
      end do
  
  end do
!$omp end do
!$omp end parallel

  call timing(iproc,'norm_trans','OF')

end subroutine normalize_transposed



subroutine init_matrixindex_in_compressed_fortransposed(iproc, nproc, orbs, collcom, collcom_shamop, &
           collcom_sr, sparsemat)
  use module_base
  use module_types
  use module_interfaces, except_this_one => init_matrixindex_in_compressed_fortransposed
  use sparsematrix_base, only: sparseMatrix
  use sparsematrix_init, only: compressed_index
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
  character(len=*),parameter :: subname='initSparseMatrix'


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

  !!allocate(sparsemat%matrixindex_in_compressed_fortransposed(imin:imax,imin:imax), stat=istat)
  !!call memocc(istat, sparsemat%matrixindex_in_compressed_fortransposed, &
  sparsemat%matrixindex_in_compressed_fortransposed=f_malloc_ptr((/imin.to.imax,imin.to.imax/),&
      id='sparsemat%matrixindex_in_compressed_fortransposed')


  do iorb=imin,imax
      do jorb=imin,imax
          sparsemat%matrixindex_in_compressed_fortransposed(iorb,jorb)=compressed_index(iorb,jorb,orbs%norb,sparsemat)
          !sendbuf(jorb,iorb)=compressed_index(jorb,iiorb,orbs%norb,sparsemat)
          !sendbuf(iorb,jorb)=compressed_index(iiorb,jorb,orbs%norb,sparsemat)
      end do
  end do

end subroutine init_matrixindex_in_compressed_fortransposed
