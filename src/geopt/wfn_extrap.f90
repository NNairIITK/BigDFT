!> @file
!!  Routines for wavefunction extrapolation for BOMD
!! @author
!!    Nisanth Nair
!!    Copyright (C) 2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
MODULE wfn_extrap
use module_base
use module_types
use yaml_output
use bigdft_run
use module_interfaces
use communications, only: transpose_v, untranspose_v
use communications_base, only: comms_cubic
implicit none

CONTAINS


!> Rotate psi_in wavefunction by using the overlap (ovrlp) with the ref. basis 
SUBROUTINE rotate_wavefunction(orbs,comms,lzd,nproc,iproc,nspin,norbTot,ndim_ovrlp, &
                               ovrlp,alpha,psi_in,psi_out)

integer :: nproc
type(orbitals_data), intent(in) :: orbs
type(local_zone_descriptors) :: lzd
type(comms_cubic), intent(in) :: comms
integer  :: nspin, iproc
real(wp) :: alpha


integer,dimension(orbs%nspin)        :: norbTot
integer,dimension(nspin,0:orbs%nkpts):: ndim_ovrlp

!real(wp), dimension(orbs%npsidim_comp)  :: psi_in
real(wp), dimension(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f,orbs%nspinor*orbs%norbp) :: psi_in
real(wp), dimension(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f,orbs%nspinor*orbs%norbp) :: psi_out
!real(wp), dimension(orbs%npsidim_comp), intent(inout) :: psi_out
real(kind=8),dimension(ndim_ovrlp(nspin,orbs%nkpts)):: ovrlp

integer :: ikptp, ikpt, ist, ispin, nvctrp, norb, norbs, ncomp, nspinor, ii, jj

real(wp), dimension(:,:) , allocatable :: psit,psit_out,psi_out2
real(wp), dimension(:) , allocatable :: work

real(wp) :: norm1, norm2

!work=f_malloc_ptr(max(orbs%npsidim_orbs,orbs%npsidim_comp),id='work')
!psit = f_malloc_ptr(max(orbs%npsidim_orbs,orbs%npsidim_comp),id='psit')
!psit_out = f_malloc_ptr(max(orbs%npsidim_orbs,orbs%npsidim_comp),id='psit_out')
!psi_out2 = f_malloc_ptr(max(orbs%npsidim_orbs,orbs%npsidim_comp),id='psi_out2')
work=f_malloc(max(orbs%npsidim_orbs,orbs%npsidim_comp),id='work')
psit=f_malloc((/Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f,orbs%nspinor*orbs%norbp/),id='psit')
psit_out=f_malloc((/Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f,orbs%nspinor*orbs%norbp/),id='psit_out')
psi_out2=f_malloc((/Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f,orbs%nspinor*orbs%norbp/),id='psi_out2')

call transpose_v(iproc,nproc,orbs,lzd%glr%wfd,comms,psi_in(1,1),&
                 work(1),out_add=psit(1,1))
!call transpose_v(iproc,nproc,orbs,lzd%glr%wfd,comms,psi_out(1),&
!                 work(1),out_add=psit_out(1))

call f_zero(psit_out) 
call f_zero(psi_out2) 

!if(iproc.eq.0)print *, "NNdbg: size of psi_in=",size(psi_in), "size psit=",size(psit)
ist=1
do ikptp=1,orbs%nkptsp
   ikpt=orbs%iskpts+ikptp
   do ispin=1,nspin
!      print *, "ist =", ist
      call orbitals_and_components(iproc,ikpt,ispin,orbs,comms,&
                                   nvctrp,norb,norbs,ncomp,nspinor)
!        call dgemm('n','n',nvctrp,norb,norb,alpha,psit(ist),nvctrp, &
!                   ovrlp(ndim_ovrlp(ispin,ikpt-1)+1),norb,1.d0,psit_out(ist),nvctrp)

!TODO Dgemm is currently in testing mode:
!Nndbg the following is only for testing.. alpha->1.0 and beta-> 0.0 instead of
!alpha, and 1.d0, respectively.
        call gemm('n','n',nvctrp,norb,norb,1.d0,psit(1,1),nvctrp, &
                   ovrlp(ndim_ovrlp(ispin,ikpt-1)+1),norb,0.d0,psit_out(1,1),nvctrp)

        call untranspose_v(iproc,nproc,orbs,lzd%glr%wfd,comms,psit_out(1,1),&
                    work(1),out_add=psi_out2(1,1)) !NN: TODO should be outside this loop?

!!         do ii=1,(lzd%Glr%wfd%nvctr_c+7*lzd%Glr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp
!         do jj=1,orbs%nspinor*orbs%norbp
!           norm1=0.0_wp
!           norm2=0.0_wp
!           do ii=1,Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f
!             norm1=norm1+psi_out(ii,jj)*psi_out(ii,jj)
!             norm2=norm2+psi_out2(ii,jj)*psi_out2(ii,jj)
!           end do
!           print *, "NNdbg: <Psi_i|Psi_i> OF psi_out (in): =", norm1, " jj=",jj,"in iproc=",iproc
!           print *, "NNdbg: <Psi_i|Psi_i> OF psi_out2 (local): =", norm2, " jj=",jj,"in iproc=",iproc
!         end do
!         print *, "NNdbg: psi_out := psi_out +", alpha , "x   psi_out2"
         do jj=1,orbs%nspinor*orbs%norbp
!           do ii=1,comms%nvctr_par(iproc,0)
            do ii=1,Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f
!              if(ii.lt.5.and.iproc.eq.0)print *, "NNdbg: psi_out(in) =",psi_out(ii,jj), " psi_out2 (local)=",psi_out2(ii,jj), " psi_in=",psi_in(ii,jj)
              psi_out(ii,jj)=psi_out(ii,jj)+alpha*psi_out2(ii,jj)
!              if(ii.lt.5.and.iproc.eq.0)print *, "NNdbg: psi_out(out) =",psi_out(ii,jj), " iproc=",iproc
            end do
         end do
!         do jj=1,orbs%nspinor*orbs%norbp
!           norm1=0.0_wp
!!           do ii=1,comms%nvctr_par(iproc,0)
!           do ii=1,Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f
!             norm1=norm1+psi_out(ii,jj)*psi_out(ii,jj)
!           end do
!           print *, "NNdbg: <Psi_i|Psi_i> OF psi_out (added): =", norm1, " ii=", jj," in iproc=",iproc
!         end do
!         do jj=1,orbs%nspinor*orbs%norbp
!           norm1=0.0_wp
!           do ii=1,Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f
!             !psi_out2(ii,jj)=psi_out2(ii,jj)+psi_out2(ii,jj)
!             psi_out2(ii,jj)=2.0_wp*psi_out2(ii,jj)
!             norm1=norm1+psi_out2(ii,jj)*psi_out2(ii,jj)
!           end do
!           print *, "NNdbg: <Psi_i|Psi_i> OF psi_out2+psi_out2 (added): =", norm1, " ii=", jj," in iproc=",iproc
!         end do

!      else
!        call zgemm('n','n',nvctrp,norb,norb,alpha,psit(ist),nvctrp, &
!                   ovrlp(ndim_ovrlp(ispin,ikpt-1)+1,1),norb,1.d0,psit_out(ist))
!      end if
      ist=ist+nvctrp*norbTot(ispin)*nspinor 
   end do
end do

!if(iproc==0)print *, "calling daxpy: alpha =",alpha
!call axpy(orbs%npsidim_comp,alpha,psi_out2(1),1,psi_out(1),1)
!do ist=1,orbs%npsidim_comp



!call untranspose_v(iproc,nproc,orbs,lzd%glr%wfd,comms,psit_out(1),&
!                    work(1),out_add=psi_out(1))

call f_free(psit)
call f_free(psit_out)
call f_free(work)
call f_free(psi_out2)

END SUBROUTINE rotate_wavefunction

SUBROUTINE extrapolation_coeff(istep,wfn_history,cc)
  integer :: istep, wfn_history
  real(wp), dimension(0:5), intent(out) :: cc

  real(wp), dimension(2), parameter :: c1 = (/-1.0_wp,2.0_wp/) 
  real(wp), dimension(3), parameter :: c2 = (/0.5_wp,-2.0_wp,2.5_wp/) 
  real(wp), dimension(5), parameter :: c4 = (/0.07142857_wp,-0.57142857_wp,1.92857143_wp,-3.42857143_wp,3.00000000_wp/)

  cc(0)=1._wp
  cc(1:5)=0._wp
  SELECT CASE (wfn_history)
  CASE(1)
    if(istep==1)then
      cc(0)=c1(1)
      cc(1)=c1(2)
    else if(istep==0)then
      cc(0)=c1(2)
      cc(1)=c1(1)
    end if
  CASE(2)
    if(istep==2)then
      cc(0)=c2(1)
      cc(1)=c2(2)
      cc(2)=c2(3)
    else if(istep==0)then
      cc(0)=c2(3)
      cc(1)=c2(1)
      cc(2)=c2(2)
    else if(istep==1)then
      cc(0)=c2(2)
      cc(1)=c2(3)
      cc(2)=c2(1)
    end if
  CASE(4)
    if(istep==4)then
      cc(0)=c4(1)
      cc(1)=c4(2)
      cc(2)=c4(3)
      cc(3)=c4(4)
      cc(4)=c4(5)
    else if(istep==0)then
      cc(0)=c4(5)
      cc(1)=c4(1)
      cc(2)=c4(2)
      cc(3)=c4(3)
      cc(4)=c4(4)
    else if(istep==1)then
      cc(0)=c4(4)
      cc(1)=c4(5)
      cc(2)=c4(1)
      cc(3)=c4(2)
      cc(4)=c4(3)
    else if(istep==2)then
      cc(0)=c4(3)
      cc(1)=c4(4)
      cc(2)=c4(5)
      cc(3)=c4(1)
      cc(4)=c4(2)
    else if(istep==3)then
      cc(0)=c4(2)
      cc(1)=c4(3)
      cc(2)=c4(4)
      cc(3)=c4(5)
      cc(4)=c4(1)
    end if
  CASE default
     call f_err_throw('Wavefunction extrapolation order is not implemented', &
                       err_name='BIGDFT_INPUT_VARIABLES_ERROR')
  END SELECT

END SUBROUTINE extrapolation_coeff

!NNdbg This is a debugging subroutine. will be removed after the full
!developement
SUBROUTINE myoverlap(iproc,nproc,lzd,orbs,comms,psi,phi)
!use module_base
!use module_types
!use yaml_output
!use bigdft_run
!use module_interfaces
!use communications, only: transpose_v
implicit none


INTEGER, INTENT(in) :: iproc, nproc
!, iwfn, nwfn 
!type(atoms_data), intent(in) :: atoms
type(local_zone_descriptors) :: lzd
type(orbitals_data), intent(in) :: orbs
type(comms_cubic), intent(in) :: comms
!type(old_wavefunction), intent(inout) :: oldpsis(nwfn)
!real(wp), dimension(orbs%npsidim_comp), intent(inout) :: psi
!real(wp), dimension(orbs%npsidim_comp), intent(inout) :: phi
!Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f,orbs%nspinor*orbs%norbp
real(wp), dimension((Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp), intent(inout) :: psi
real(wp), dimension((Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp), intent(inout) :: phi
real(wp), dimension(:) , pointer :: phit,psit,work_array

!real(gp) :: rxyz(3,atoms%astruct%nat)
logical, parameter :: debug_flag=.false.
integer, save :: icall = 0
integer :: ispin, nspin, jwfn, nsize
!real(gp) :: c(nwfn)
!real(wp), allocatable :: psi_istep(:), psi_nstep(:)

integer, dimension(:,:), allocatable :: ndim_ovrlp
real(wp), dimension(:), allocatable :: ovrlp
integer,dimension(:),allocatable:: norbArr

icall=icall+1


nspin=1
if(orbs%norbd>0)nspin=2

! ndim_ovrlp describes the shape of the overlap matrix.
ndim_ovrlp = f_malloc((/ 1.to.nspin, 0.to.orbs%nkpts /),id='ndim_ovrlp')
!print *, "NNdbg: orbs%nkpts= ", orbs%nkpts
!print *, "NNdbg: nspin= ", nspin
!print *, "NNdbg: size(ndim_ovrlp)= ", size(ndim_ovrlp)
call dimension_ovrlp(nspin,orbs,ndim_ovrlp)
!print *, "NNdbg: in myoverlap: dimension_ovrlp=",ndim_ovrlp(1, orbs%nkpts)

! Allocate norbArr which contains the number of up and down orbitals.
norbArr = f_malloc(nspin,id='norbArr')
do ispin=1,nspin
  if(ispin==1) norbArr(ispin)=orbs%norbu
  if(ispin==2) norbArr(ispin)=orbs%norbd
end do

!allocate transposed principal wavefunction
work_array =f_malloc_ptr(max(orbs%npsidim_orbs,orbs%npsidim_comp),id='work_array')
psit = f_malloc_ptr(max(orbs%npsidim_orbs,orbs%npsidim_comp),id='psit')
phit = f_malloc_ptr(max(orbs%npsidim_orbs,orbs%npsidim_comp),id='phit')
call transpose_v(iproc,nproc,orbs,lzd%glr%wfd,comms,psi(1),&
          &   work_array(1),out_add=psit(1))
call transpose_v(iproc,nproc,orbs,lzd%glr%wfd,comms,phi(1),&
          &   work_array(1),out_add=phit(1))

! Allocate overlap matrix
ovrlp = f_malloc(ndim_ovrlp(nspin, orbs%nkpts),id='ovrlp')

if(iproc.eq.0 .and. debug_flag)print *, "!!!!NNdbg: calling my getoverlap!!!!",iproc

call Overlap_PhiPsi(iproc,nproc,nspin,norbArr(1),orbs,comms,&
     phit,psit,ndim_ovrlp,ovrlp,norbArr,1,1)

if(iproc.eq.0 .and. debug_flag)then
  print *, "!!!!NNdbg: my overlap", ovrlp(1:ndim_ovrlp(nspin,orbs%nkpts))
end if

if(iproc.eq.0 .and. debug_flag)print *, "done getoverlap",iproc


!deallocate(psi_istep)
!deallocate(psi_nstep)
call f_free(ndim_ovrlp)
call f_free(norbArr)
call f_free(ovrlp)
call f_free_ptr(psit)
call f_free_ptr(phit)
call f_free_ptr(work_array)
return
END SUBROUTINE myOverlap


!ovrlp=<Psit|Psit> Note: psi have to be transposed
!NOTE: this routine is not used at the moment
subroutine Overlap_PsiPsi(iproc,nproc,nspin,norbIn,orbs,comms,&
     psi,ndim_ovrlp,ovrlp,norbTot,block1,ispinIn)

!  use module_base
!  use module_types
!  implicit none

  ! Calling arguments
  integer,intent(in):: iproc,nproc,nspin,norbIn,block1,ispinIn
  type(orbitals_data),intent(in):: orbs
  type(comms_cubic),intent(in) :: comms
  real(wp),dimension(comms%nvctr_par(iproc,0)*orbs%nspinor*orbs%norb),intent(in) :: psi
  integer,dimension(nspin,0:orbs%nkpts),intent(in):: ndim_ovrlp
  real(wp),dimension(ndim_ovrlp(nspin,orbs%nkpts)),intent(out):: ovrlp
  integer,dimension(nspin),intent(in):: norbTot

  ! Local variables
  logical, parameter :: debug_flag=.false.
  integer:: ispsi,ikptp,ikpt,ispin,nspinor,ncomp,norbs,nvctrp,norb



  ! Set the whole overlap matrix to zero. This is necessary since each process treats only a part
  ! of the matrix.
  call f_zero(ovrlp)

!NNdbg
!  print *, "NNdbg: go Syrk", iproc
!  print *, "NNdbg: check sizes: comms%nvctr_par(iproc,0)", comms%nvctr_par(iproc,0), iproc
!  print *, "NNdbg: check sizes: orbs%nspinor", orbs%nspinor, iproc
!  print *, "NNdbg: check sizes: orbs%norb",orbs%norb, iproc
!  print *, "NNdbg: check sizes: ndim_ovrlp(nspin,orbs%nkpts) =",ndim_ovrlp(nspin,orbs%nkpts), iproc
!NNdbg

  ispsi=1
  ! First make a loop over the k points handled by this process.
  do ikptp=1,orbs%nkptsp
     ! ikpt is the index of the k point.
     ikpt=orbs%iskpts+ikptp

     ! Now make also a loop over spin up/down.
     do ispin=1,nspin

        ! This subroutine gives essentially back nvctrp, i.e. the length of the vectors for which the overlap
        ! matrix shall be calculated. In addition it sets the value of nspinor to orbs%nspinor.
        call orbitals_and_components(iproc,ikpt,ispin,orbs,comms,&
             nvctrp,norb,norbs,ncomp,nspinor)
!NNdbg
!  print *, "NNdbg: orb and comp done", iproc
!NNdbg
        ! The subroutine also overwrite the variable norb with the total number of orbitals.
        ! However we want to keep the value of norbIn (since we treat only a part of the orbitals).
        norb=norbIn
        ! Put the starting index to the right place. The current block of vector starts at the block1-th vector.
        ispsi=ispsi+nvctrp*(block1-1)*nspinor
        if(ispin==ispinIn) then
           if (nvctrp == 0) cycle

           ! Now calclulate one part of the overlap matrix. The starting index of this part is given by ndim_ovrlp(ispin,ikpt-1)+1.
           if(nspinor==1) then
              call syrk('L','T',norb,nvctrp,1.0_wp,psi(ispsi),max(1,nvctrp),&
                   0.0_wp,ovrlp(ndim_ovrlp(ispin,ikpt-1)+1),norb)
!NNdbg
!       print *, "NNdbg: syrk done ", iproc
!NNdbg
           else
              call herk('L','C',norb,ncomp*nvctrp,1.0_wp,psi(ispsi),&
                   max(1,ncomp*nvctrp),0.0_wp,ovrlp(ndim_ovrlp(ispin,ikpt-1)+1),norb)
           end if
        end if
        ! Move the starting indices to the end of the actual k point. This is necessary since nvctrp is
        ! different for the next k point and we cannot jump directly to the starting indices of our block for 
        ! the next k point.
        ispsi=ispsi+nvctrp*(norbTot(ispin)-block1+1)*nspinor
     end do
  end do

!  print *, "size overlap in", iproc," is =",size(ovrlp)
!  print *, "ndim_ovrlp(nspin,nkpt) in", iproc," is =",ndim_ovrlp(nspin,orbs%nkpts)

  !call MPI_BARRIER(bigdft_mpi%mpi_comm,ierr)
  !print *,'here',iproc
!NNdbg
  !call MPI_BARRIER(bigdft_mpi%mpi_comm,ispsi)
!  print *, "NNdbg: loop done", iproc
!NNdbg

  if (nproc > 1) then
     call mpiallred(ovrlp,MPI_SUM,comm=bigdft_mpi%mpi_comm)
  end if
!NNdbg
  if(iproc.eq.0  .and. debug_flag)print *, "overlap-inside", ovrlp(1:ndim_ovrlp(nspin,orbs%nkpts))
!  print *, "NNdbg: getOverlap done", iproc
!NNdbg

  ! Now each processors knows all the overlap matrices for each k-point
  ! even if it does not handle it.
  ! This is somehow redundant but it is one way of reducing the number of communications
  ! without defining group of processors.

END SUBROUTINE Overlap_PsiPsi



!ovrlp=<Phit|Psit> Note: phi and psi have to be transposed
subroutine Overlap_PhiPsi(iproc,nproc,nspin,norbIn,orbs,comms,&
     phi,psi,ndim_ovrlp,ovrlp,norbTot,block1,ispinIn)

!  use module_base
!  use module_types
!  use communications_base, only: comms_cubic
!  implicit none

  ! Calling arguments
  integer,intent(in):: iproc,nproc,nspin,norbIn,block1,ispinIn
  type(orbitals_data),intent(in):: orbs
  type(comms_cubic),intent(in) :: comms
  real(wp),dimension(comms%nvctr_par(iproc,0)*orbs%nspinor*orbs%norb),intent(in) :: psi
  real(wp),dimension(comms%nvctr_par(iproc,0)*orbs%nspinor*orbs%norb),intent(in) :: phi
  integer,dimension(nspin,0:orbs%nkpts),intent(in):: ndim_ovrlp
  real(wp),dimension(ndim_ovrlp(nspin,orbs%nkpts)),intent(out):: ovrlp
  integer,dimension(nspin),intent(in):: norbTot

  ! Local variables
  integer:: ispsi,ikptp,ikpt,ispin,nspinor,ncomp,norbs,nvctrp,norb


  ! Set the whole overlap matrix to zero. This is necessary since each process treats only a part
  ! of the matrix.
  call f_zero(ovrlp)

!NNdbg
  !print *, "NNdbg 2: go Syrk", iproc
  !print *, "NNdbg 2: check sizes: comms%nvctr_par(iproc,0)", comms%nvctr_par(iproc,0), iproc
  !print *, "NNdbg 2: check sizes: orbs%nspinor", orbs%nspinor, iproc
  !print *, "NNdbg 2: check sizes: orbs%norb",orbs%norb, iproc
  !print *, "NNdbg 2: check sizes: ndim_ovrlp(nspin,orbs%nkpts) =",ndim_ovrlp(nspin,orbs%nkpts), iproc
!NNdbg

  ispsi=1
  do ikptp=1,orbs%nkptsp
     ikpt=orbs%iskpts+ikptp
     do ispin=1,nspin
        call orbitals_and_components(iproc,ikpt,ispin,orbs,comms,&
             nvctrp,norb,norbs,ncomp,nspinor)
!NNdbg
!  print *, "NNdbg 2: orb and comp done", iproc
!NNdbg
        norb=norbIn
        ! Put the starting index to the right place. The current block of vector starts at the block1-th vector.
        ispsi=ispsi+nvctrp*(block1-1)*nspinor
        if(ispin==ispinIn) then
           if (nvctrp == 0) cycle
           if(nspinor==1) then
              call gemm('t','n',norb,norb,ncomp*nvctrp,1.0_wp,psi(ispsi),&
                   ncomp*nvctrp,phi(ispsi),ncomp*nvctrp,0.d0,ovrlp(ndim_ovrlp(ispin,ikpt-1)+1),norb)
!NNdbg
!              print *, "NNdbg 2: gemm done ", iproc
!NNdbg
           else !NNdbg: TODO: not tested
              call c_gemm('c','n',norb,norb,ncomp*nvctrp,(1.0_wp,0.0_wp),psi(ispsi),&
                   ncomp*nvctrp,phi(ispsi),ncomp*nvctrp,(0.d0,0.d0),ovrlp(ndim_ovrlp(ispin,ikpt-1)+1),norb)
           end if
        end if
        ispsi=ispsi+nvctrp*(norbTot(ispin)-block1+1)*nspinor
     end do
  end do

  if (nproc > 1) then
     call mpiallred(ovrlp,MPI_SUM,comm=bigdft_mpi%mpi_comm)
  end if
!NNdbg
!  if(iproc.eq.0)print *, "overlap-inside 2", ovrlp(1:ndim_ovrlp(nspin,orbs%nkpts))
!  print *, "NNdbg 2: getOverlap 2 done", iproc
!NNdbg

  ! Now each processors knows all the overlap matrices for each k-point
  ! even if it does not handle it.
  ! This is somehow redundant but it is one way of reducing the number of communications
  ! without defining group of processors.

END SUBROUTINE Overlap_PhiPsi


!> Compute overlap matrix ovrlp=<Phi|Psi>
!Input wfns Phi and Psi should NOT be transposed
!This is a wrapper for Overlap_PhiPsi routine
SUBROUTINE myoverlap2(iproc,nproc,lzd,orbs,comms,nspin,norbArr,ndim_ovrlp,psi,phi,ovrlp)
implicit none

INTEGER, INTENT(in) :: iproc, nproc
type(local_zone_descriptors) :: lzd
type(orbitals_data), intent(in) :: orbs
type(comms_cubic), intent(in) :: comms
real(wp), dimension((Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp), intent(inout) :: psi
real(wp), dimension((Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp), intent(inout) :: phi
integer :: nspin
integer,dimension(nspin,0:orbs%nkpts),intent(in):: ndim_ovrlp
real(wp),dimension(ndim_ovrlp(nspin,orbs%nkpts)),intent(out):: ovrlp
!local variables
logical, parameter :: debug_flag=.true.
integer,dimension(orbs%nspin)        :: norbArr

real(wp), dimension(:) , pointer :: psit,phit,work_array

integer :: ispin, jwfn, nsize

!allocate transposed principal wavefunction
work_array =f_malloc_ptr(max(orbs%npsidim_orbs,orbs%npsidim_comp),id='work_array')
psit = f_malloc0_ptr(max(orbs%npsidim_orbs,orbs%npsidim_comp),id='psit')
phit = f_malloc0_ptr(max(orbs%npsidim_orbs,orbs%npsidim_comp),id='phit')
call transpose_v(iproc,nproc,orbs,lzd%glr%wfd,comms,psi(1),&
          &   work_array(1),out_add=psit(1))
call transpose_v(iproc,nproc,orbs,lzd%glr%wfd,comms,phi(1),&
          &   work_array(1),out_add=phit(1))
!

!if(iproc.eq.0)print *, "!!!!NNdbg: calling my getoverlap2!!!!",iproc
!
call Overlap_PhiPsi(iproc,nproc,nspin,norbArr(1),orbs,comms,&
     phit,psit,ndim_ovrlp,ovrlp,norbArr,1,1)

if(iproc==0 .and. debug_flag)then
  print *, "!!!!NNdbg: overlap", ovrlp(1:ndim_ovrlp(nspin,orbs%nkpts))
end if

!if(iproc.eq.0)print *, "done getoverlap",iproc

call f_free_ptr(psit)
call f_free_ptr(phit)
call f_free_ptr(work_array)
END SUBROUTINE myOverlap2



!> Normalize psi
SUBROUTINE normalize_wavefunction(orbs,comms,lzd,nproc,iproc,psi)

integer :: nproc
type(orbitals_data), intent(in) :: orbs
type(local_zone_descriptors) :: lzd
type(comms_cubic), intent(in) :: comms
integer  :: nspin, iproc
real(wp) :: alpha


real(wp), dimension(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f,orbs%nspinor*orbs%norbp) :: psi

integer :: ikptp, ikpt, ist, ispin, nvctrp, norb, norbs, ncomp, nspinor, ii, jj

real(wp) :: norm1, norm2

!TODO cleanup using dnrm2()
do jj=1,orbs%nspinor*orbs%norbp
  norm1=0.0_wp
  do ii=1,Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f
    norm1=norm1+psi(ii,jj)*psi(ii,jj)
  end do
  norm1=1.0_wp/sqrt(norm1)
  call vscal(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f,norm1,psi(1,jj),1)
end do
END SUBROUTINE normalize_wavefunction

end module wfn_extrap
