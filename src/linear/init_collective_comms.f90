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
   use communications_base, only: comms_linear, TRANSPOSE_FULL
   use communications, only: transpose_localized, untranspose_localized
   use sparsematrix_base, only : sparse_matrix, matrices, DENSE_PARALLEL
   use sparsematrix, only : compress_matrix_distributed, gather_matrix_from_taskgroups_inplace
   use transposed_operations, only: calculate_overlap_transposed
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
               ist=(ispin-1)*smat%nvctr+smat%isvctrp_tg+1
               call compress_matrix_distributed(iproc, nproc, smat, DENSE_PARALLEL, &
                    matp, mat_compr(ist:))
           end do
           maxdiff=0.d0
           call f_free(psiig)
           call f_free(psijg)
           call f_free(matp)
           do i=1,smat%nvctrp_tg
               maxdiff=max(abs(mat_compr(i+smat%isvctrp_tg)-mat%matrix_compr(i)),maxdiff)
               !write(8000+iproc,'(a,i7,2es15.5)') 'i, mat_compr(i), mat%matrix_compr(i)', &
               !    i, mat_compr(i), mat%matrix_compr(i)
           end do
           if (nproc>1) then
               call mpiallred(maxdiff, 1, mpi_max, bigdft_mpi%mpi_comm)
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
          call yaml_comment('ERROR (Inverse Transposition): process'//trim(yaml_toa(iproc))//&
               ' found an error of:'//trim(yaml_toa(maxdiff,fmt='(1pe15.7)')))
          abort = .true.
       end if
    
   if (abort) call MPI_ABORT(bigdft_mpi%mpi_comm,11,ierr)
    
       if (nproc > 1) then
          !call MPI_BARRIER(bigdft_mpi%mpi_comm, ierr)
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
  call f_zero(ovrlp)
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

  call f_routine(id='init_matrixindex_in_compressed_fortransposed')


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

  call f_release_routine()

end subroutine init_matrixindex_in_compressed_fortransposed


subroutine synchronize_matrix_taskgroups(iproc, nproc, smat, mat)
  use module_base
  use sparsematrix_base, only: sparse_matrix, matrices
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(sparse_matrix),intent(in) :: smat
  type(matrices),intent(in) :: mat

  ! Local variables
  integer :: ncount, itg, iitg, ispin, ishift, ist_send, ist_recv
  integer,dimension(:),allocatable :: request
  real(kind=8),dimension(:),allocatable :: recvbuf

  if (nproc>1) then
      request = f_malloc(smat%ntaskgroupp,id='request')
      ncount = 0
      do itg=1,smat%ntaskgroupp
          iitg = smat%taskgroupid(itg)
          ncount = ncount + smat%taskgroup_startend(2,1,iitg)-smat%taskgroup_startend(1,1,iitg)+1
      end do
      recvbuf = f_malloc(ncount,id='recvbuf')
      do ispin=1,smat%nspin
          ishift = (ispin-1)*smat%nvctrp_tg

          ncount = 0
          do itg=1,smat%ntaskgroupp
              iitg = smat%taskgroupid(itg)
              ist_send = smat%taskgroup_startend(1,1,iitg) - smat%isvctrp_tg
              ist_recv = ncount + 1
              ncount = smat%taskgroup_startend(2,1,iitg)-smat%taskgroup_startend(1,1,iitg)+1
              !!call mpi_iallreduce(mat%matrix_compr(ist_send), recvbuf(ist_recv), ncount, &
              !!     mpi_double_precision, mpi_sum, smat%mpi_groups(iitg)%mpi_comm, request(itg), ierr)
              if (nproc>1) then
                  call mpiiallred(mat%matrix_compr(ishift+ist_send), recvbuf(ist_recv), ncount, &
                       mpi_sum, smat%mpi_groups(iitg)%mpi_comm, request(itg))
              else
                  call vcopy(ncount, mat%matrix_compr(ishift+ist_send), 1, recvbuf(ist_recv), 1)
              end if
          end do
          if (nproc>1) then
              call mpiwaitall(smat%ntaskgroupp, request)
          end if
          ncount = 0
          do itg=1,smat%ntaskgroupp
              iitg = smat%taskgroupid(itg)
              ist_send = smat%taskgroup_startend(1,1,iitg) - smat%isvctrp_tg
              ist_recv = ncount + 1
              ncount = smat%taskgroup_startend(2,1,iitg)-smat%taskgroup_startend(1,1,iitg)+1
              !call vcopy(ncount, recvbuf(ist_recv), 1, mat%matrix_compr(ishift+ist_send), 1)
              call dcopy(ncount, recvbuf(ist_recv), 1, mat%matrix_compr(ishift+ist_send), 1)
          end do
      end do
      call f_free(request)
      call f_free(recvbuf)
  end if
end subroutine synchronize_matrix_taskgroups
