!!****f* BigDFT/davidson
!! AUTHOR
!!   Alexander Willand
!! DESCRIPTION
!!   Davidsons method for iterative diagonalization of virtual Kohn Sham orbitals
!!   under orthogonality constraints to occupied orbitals psi. The nvirt input
!!   variable gives the number of unoccupied orbitals for which the exit criterion
!!   for the gradients norm holds. nvirte = norbe - norb >= nvirt is the number of
!!   virtual orbitals processed by the method. The dimension of the subspace for
!!   diagonalization is 2*nvirte = n2virt
!!                                                                   Alex Willand
!!   Algorithm
!!   _________
!!   (parallel)
    
    
!!   (transpose psi, v is already transposed)
!!   orthogonality of v to psi
!!   orthogonalize v
!!   (retranspose v)
!!   Hamilton(v) --> hv
!!   transpose v and hv
!!   Rayleigh quotients  e
!!   do
!!      gradients g= e*v -hv
!!      exit condition gnrm
!!      orthogonality of g to psi
!!      (retranspose g)
!!      preconditioning of g
!!      (transpose g again)
!!      orthogonality of g to psi
!!      (retranspose g)
!!      Hamilton(g) --> hg
!!      (transpose g and hg)
!!      subspace matrices H and S
!!      DSYGV(H,e,S)  --> H
!!      update v with eigenvectors H
!!      orthogonality of v to psi
!!      orthogonalize v
!!      (retranspose v)
!!      Hamilton(v) --> hv
!!      (transpose v and hv)
!!   end do
!!   (retranspose v and psi)
!!
!! COPYRIGHT
!!    Copyright (C) 2007-2010 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! SOURCE
!!
subroutine davidson(iproc,nproc,n1i,n2i,in,at,&
     orbs,orbsv,nvirt,lr,comms,&
     hx,hy,hz,rxyz,rhopot,i3xcsh,n3p,nlpspd,proj,pkernel,psi,v,ngatherarr,GPU)
  use module_base
  use module_types
  use module_interfaces, except_this_one => davidson
  use libxc_functionals
  implicit none
  integer, intent(in) :: iproc,nproc,n1i,n2i
  integer, intent(in) :: i3xcsh
  integer, intent(in) :: nvirt,n3p
  type(input_variables), intent(in) :: in
  type(atoms_data), intent(in) :: at
  type(nonlocal_psp_descriptors), intent(in) :: nlpspd
  type(locreg_descriptors), intent(in) :: lr 
  type(orbitals_data), intent(in) :: orbs
  type(communications_arrays), intent(in) :: comms
  real(gp), intent(in) :: hx,hy,hz
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
  real(dp), dimension(*), intent(in) :: pkernel,rhopot
  type(orbitals_data), intent(inout) :: orbsv
  type(GPU_pointers), intent(inout) :: GPU
  !this is a Fortran 95 standard, should be avoided (it is a pity IMHO)
  !real(kind=8), dimension(:,:,:,:), allocatable :: rhopot 
  real(wp), dimension(:), pointer :: psi,v!=psivirt(nvctrp,nvirtep*nproc) 
                        !v, that is psivirt, is transposed on input and direct on output
  !local variables
  character(len=*), parameter :: subname='davidson'
  character(len=10) :: comment
  character(len=11) :: orbname,denname
  logical :: msg,exctX,occorbs !extended output
  integer :: ierr,i_stat,i_all,iorb,jorb,iter,nwork,ind,norb,nspinor
  integer :: ise,j,ispsi,ikpt,ikptp,nvctrp,ncplx,ncomp,norbs,ispin,ish1,ish2,nspin
  real(gp) :: tt,gnrm,epot_sum,eexctX,ekin_sum,eproj_sum,gnrm_fake
  type(communications_arrays) :: commsv
  integer, dimension(:,:), allocatable :: ndimovrlp
  real(wp), dimension(:), allocatable :: work,work_rp,hamovr
  real(wp), dimension(:), allocatable :: hv,g,hg,ew,psirocc
  real(wp), dimension(:,:,:), allocatable :: e
  real(wp), dimension(:), pointer :: psiw
!OCL  integer, dimension(3) :: periodic

  !logical flag which control to othogonalise wrt the occupied orbitals or not
  if (orbs%nkpts /= orbsv%nkpts) then
     occorbs=.false.
  else
     occorbs=.true.
     do ikpt = 1, orbs%nkpts
        if (abs(maxval(orbs%kpts(:,ikpt) - orbsv%kpts(:,ikpt))) > 1.d-6) then
           occorbs=.false.
           exit
        end if
     end do
  end if

  !in the GPU case, the wavefunction should be copied to the card 
  !at each HamiltonianApplication
  !rebind the GPU pointers to the orbsv structure
  if (GPUconv) then
     call free_gpu(GPU,orbs%norbp)
     call prepare_gpu_for_locham(lr%d%n1,lr%d%n2,lr%d%n3,in%nspin,&
          hx,hy,hz,lr%wfd,orbsv,GPU)
  else if (OCLconv) then
     call free_gpu_OCL(GPU,orbs%norbp)    
     call allocate_data_OCL(lr%d%n1,lr%d%n2,lr%d%n3,at%geocode,&
          in%nspin,hx,hy,hz,lr%wfd,orbsv,GPU)
  end if
 
  GPU%full_locham=.true.
  !verify whether the calculation of the exact exchange term
  !should be preformed
  exctX = libxc_functionals_exctXfac() /= 0.0_gp

  !last index of e and hamovr are for mpi_allreduce. 
  !e (eigenvalues) is also used as 2 work arrays
  
  msg=verbose > 2 .and. iproc ==0! no extended output
  !msg =(iproc==0)!extended output

  if(iproc==0)write(*,'(1x,a)')"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  if(iproc==0)write(*,'(1x,a)')"Iterative subspace diagonalization of virtual orbitals."

  !if(msg)write(*,*)'shape(v)',shape(v),'size(v)',size(v)


  !before transposition, create the array of the occupied
  !wavefunctions in real space, for exact exchange calculations
  if (exctX) then
     allocate(psirocc(max(max(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%norbp,&
          ngatherarr(0,1)*orbs%norb),1)+ndebug),stat=i_stat)
     call memocc(i_stat,psirocc,'psirocc',subname)

     call prepare_psirocc(iproc,nproc,lr,orbs,n3p,ngatherarr(0,1),psi,psirocc)
  end if

  !n2virt=2*orbsv%norb! the dimension of the subspace

  !disassociate work array for transposition in serial
  if (nproc > 1) then
     allocate(psiw(orbs%npsidim+ndebug),stat=i_stat)
     call memocc(i_stat,psiw,'psiw',subname)
  else
     psiw => null()
  endif

  !transpose the wavefunction psi 
  call transpose_v(iproc,nproc,orbs,lr%wfd,comms,psi,work=psiw)

  if (nproc > 1) then
     i_all=-product(shape(psiw))*kind(psiw)
     deallocate(psiw,stat=i_stat)
     call memocc(i_stat,i_all,'psiw',subname)
  end if

  !allocate communications arrays for virtual orbitals
  call orbitals_communicators(iproc,nproc,lr,orbsv,commsv)  

  !prepare the v array starting from a set of gaussians
  call psivirt_from_gaussians(iproc,nproc,at,orbsv,lr,commsv,rxyz,hx,hy,hz,in%nspin,v)

  if(iproc==0)write(*,'(1x,a)',advance="no")"Orthogonality to occupied psi..."
  !project v such that they are orthogonal to all occupied psi
  !Orthogonalize before and afterwards.

  !here nvirte=orbsv%norb
  !     nvirtep=orbsv%norbp

  !this is the same also in serial
  call orthogonalize(iproc,nproc,orbsv,commsv,lr%wfd,v)

  if (occorbs) then
     call orthon_virt_occup(iproc,nproc,orbs,orbsv,comms,commsv,psi,v,msg)
     !and orthonormalize them using "gram schmidt"  (should conserve orthogonality to psi)
     call orthogonalize(iproc,nproc,orbsv,commsv,lr%wfd,v)
  end if

  !retranspose v
  if(nproc > 1)then
     !reallocate the work array with the good size
     allocate(psiw(orbsv%npsidim+ndebug),stat=i_stat)
     call memocc(i_stat,psiw,'psiw',subname)
  end if

  call untranspose_v(iproc,nproc,orbsv,lr%wfd,commsv,v,work=psiw)

  ! 1st Hamilton application on psivirt
  if(iproc==0)write(*,'(1x,a)',advance="no")"done. first "

  allocate(hv(orbsv%npsidim+ndebug),stat=i_stat)
  call memocc(i_stat,hv,'hv',subname)
  
  call HamiltonianApplication(iproc,nproc,at,orbsv,hx,hy,hz,rxyz,&
       nlpspd,proj,lr,ngatherarr,n1i*n2i*n3p,&
       rhopot,v,hv,ekin_sum,epot_sum,eexctX,eproj_sum,in%nspin,GPU,&
       pkernel,orbs,psirocc) ! optional arguments

  !if(iproc==0)write(*,'(1x,a)',advance="no")"done. Rayleigh quotients..."

  allocate(e(orbsv%norb,orbsv%nkpts,2+ndebug),stat=i_stat)
  call memocc(i_stat,e,'e',subname)

  !transpose  v and hv
  call transpose_v(iproc,nproc,orbsv,lr%wfd,commsv,v,work=psiw)
  call transpose_v(iproc,nproc,orbsv,lr%wfd,commsv,hv,work=psiw)

  call timing(iproc,'Davidson      ','ON')
  !Timing excludes transposition, hamilton application and preconditioning
  call razero(orbsv%norb*2*orbsv%nkpts,e)
  ! Rayleigh quotients.

  !probably this loop can be rewritten using GEMMs
  ispsi=1
  do ikptp=1,orbsv%nkptsp
     ikpt=orbsv%iskpts+ikptp!orbsv%ikptsp(ikptp)
     nvctrp=commsv%nvctr_par(iproc,ikptp)
     if (nvctrp == 0) cycle

     nspinor=orbsv%nspinor
     do iorb=1,orbsv%norb ! temporary variables 
        !for complex wavefunctions the diagonal is always real
        e(iorb,ikpt,1)= dot(nvctrp*nspinor,v(ispsi+nvctrp*nspinor*(iorb-1)),1,&
             hv(ispsi+nvctrp*nspinor*(iorb-1)),1)          != <psi|H|psi> 
        e(iorb,ikpt,2)= nrm2(nvctrp*nspinor,v(ispsi+nvctrp*nspinor*(iorb-1)),1)**2   != <psi|psi> 
     end do
     ispsi=ispsi+nvctrp*orbsv%norb*orbsv%nspinor
  end do

  if(nproc > 1)then
     !sum up the contributions of nproc sets with 
     !commsv%nvctr_par(iproc,1) wavelet coefficients each
     !call MPI_ALLREDUCE(e(1,1,1,2),e(1,1,1,1),2*orbsv%norb*orbsv%nkpts,&
     !     mpidtypw,MPI_SUM,MPI_COMM_WORLD,ierr)
     call mpiallred(e(1,1,1),2*orbsv%norb*orbsv%nkpts,MPI_SUM,MPI_COMM_WORLD,ierr)

  end if

  if(iproc==0)write(*,'(1x,a)')"done."
  if(iproc==0)write(*,'(1x,a)')"     sqnorm                Rayleigh quotient"
 
  do ikpt=1,orbsv%nkpts
     do iorb=1,orbsv%norb
        !e(:,1,1) = <psi|H|psi> / <psi|psi>
        e(iorb,ikpt,1)=e(iorb,ikpt,1)/e(iorb,ikpt,2)
        if(iproc==0) write(*,'(1x,i3,2(1x,1pe21.14))')&
             iorb,e(iorb,ikpt,2),e(iorb,ikpt,1)
     end do
  end do

!if(msg)then
!write(*,*)"******** transposed v,hv 1st elements"
!do iorb=1,10
!  write(*,*)v(iorb,1),hv(iorb,1)
!end do
!write(*,*)"**********"
!end if

  !calculate the dimension of the overlap matrix for each k-point
  if (orbsv%norbd > 0) then
     nspin=2
  else
     nspin=1
  end if

  !number of components for the overlap matrix in wp-kind real numbers

  allocate(ndimovrlp(nspin,0:orbsv%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,ndimovrlp,'ndimovrlp',subname)

  call dimension_ovrlp(nspin,orbsv,ndimovrlp)
  
  !the dimension should be chosen with the max between k-points
  !allocate(hamovr(n2virt,n2virt,2,orbsv%nkpts+ndebug),stat=i_stat)
  allocate(hamovr(8*ndimovrlp(nspin,orbsv%nkpts)+ndebug),stat=i_stat)
  call memocc(i_stat,hamovr,'hamovr',subname)

  !put to zero all the k-points which are not needed
  call razero(8*ndimovrlp(nspin,orbsv%nkpts),hamovr)

  nwork=max(10,16*orbsv%norb)
  allocate(work(nwork+ndebug),stat=i_stat)
  call memocc(i_stat,work,'work',subname)
  
  allocate(ew(2*orbsv%norb+ndebug),stat=i_stat)
  call memocc(i_stat,ew,'ew',subname)

  allocate(orbsv%eval(orbsv%norb*orbsv%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,orbsv%eval,'eval',subname)


  !itermax=... use the input variable instead
  iter=1
  davidson_loop: do 

     if(iproc==0)write(*,'(1x,a,i3)')&
     "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~iter",iter
     if(msg) write(*,'(1x,a)')"squared norm of the (nvirt) gradients"

     allocate(g(orbsv%npsidim+ndebug),stat=i_stat)
     call memocc(i_stat,g,'g',subname)

     call dcopy(orbsv%npsidim,hv,1,g,1)! don't overwrite hv

     call razero(orbsv%norb*orbsv%nkpts,e(1,1,2))
     !also these operations are presumably GEMMs
     !here we should add the ncomp term for non-collinear case
     ispsi=1
     do ikptp=1,orbsv%nkptsp
        ikpt=orbsv%iskpts+ikptp!orbsv%ikptsp(ikptp)
        nvctrp=commsv%nvctr_par(iproc,ikptp)
        if (nvctrp == 0) cycle

        nspinor=orbsv%nspinor
        do iorb=1,orbsv%norb
           !gradient = hv-e*v
           call axpy(nvctrp*nspinor,-e(iorb,ikpt,1),v(ispsi+nvctrp*nspinor*(iorb-1)),1,&
                g(ispsi+nvctrp*nspinor*(iorb-1)),1)

           !local contribution to the square norm
           e(iorb,ikpt,2)= nrm2(nvctrp*nspinor,g(ispsi+nvctrp*nspinor*(iorb-1)),1)**2
        end do
        ispsi=ispsi+nvctrp*orbsv%norb*orbsv%nspinor
     end do

     if(nproc > 1)then
        !sum up the contributions of nproc sets with nvctrp wavelet coefficients each
        call mpiallred(e(1,1,2),orbsv%norb*orbsv%nkpts,MPI_SUM,MPI_COMM_WORLD,ierr)
     end if

     gnrm=0._dp
     do ikpt=1,orbsv%nkpts
        do iorb=1,nvirt
           tt=real(e(iorb,ikpt,2)*orbsv%kwgts(ikpt),dp)
           if(msg)write(*,'(1x,i3,1x,1pe21.14)')iorb,tt
           gnrm=gnrm+tt
        end do
        if (nspin == 2) then
           do iorb=1,nvirt
              tt=real(e(iorb+orbsv%norbu,ikpt,2)*orbsv%kwgts(ikpt),dp)
              if(msg)write(*,'(1x,i3,1x,1pe21.14)')iorb+orbsv%norbu,tt
              gnrm=gnrm+tt
           end do
        end if
     end do
     !should we divide by nvrit or by orbsv%norb?
     gnrm=dsqrt(gnrm/real(orbsv%norb,dp))

     if(iproc == 0)write(*,'(1x,a,2(1x,1pe12.5))')&
          "|gradient|=gnrm and exit criterion ",gnrm,in%gnrm_cv
     if(gnrm < in%gnrm_cv) then
        i_all=-product(shape(g))*kind(g)
        deallocate(g,stat=i_stat)
        call memocc(i_stat,i_all,'g',subname)
        exit davidson_loop! iteration loop
     end if
     call timing(iproc,'Davidson      ','OF')

     if(iproc==0)write(*,'(1x,a)',advance="no")&
          "Orthogonality of gradients to occupied psi..."

     !project g such that they are orthogonal to all occupied psi. 
     !Gradients do not need orthogonality.
     if (occorbs) then
        call orthon_virt_occup(iproc,nproc,orbs,orbsv,comms,commsv,psi,g,msg)
     end if

     call timing(iproc,'Davidson      ','ON')
     if(iproc==0)write(*,'(1x,a)',advance="no")"done."
     if(msg)write(*,'(1x,a)')"squared norm of all gradients after projection"

     call razero(orbsv%norb*orbsv%nkpts,e(1,1,2))
     ispsi=1
     do ikptp=1,orbsv%nkptsp
        ikpt=orbsv%iskpts+ikptp!orbsv%ikptsp(ikptp)
        nvctrp=commsv%nvctr_par(iproc,ikptp)
        if (nvctrp == 0) cycle

        nspinor=orbsv%nspinor
        do iorb=1,orbsv%norb
           e(iorb,ikpt,2)= nrm2(nvctrp*nspinor,g(ispsi+nvctrp*nspinor*(iorb-1)),1)**2
        end do

        ispsi=ispsi+nvctrp*orbsv%norb*orbsv%nspinor
     end do

     if(nproc > 1)then
        !sum up the contributions of nproc sets with nvctrp wavelet coefficients each
        call mpiallred(e(1,1,2),orbsv%norb*orbsv%nkpts,MPI_SUM,MPI_COMM_WORLD,ierr)
     end if

     gnrm=0._dp
     do ikpt=1,orbsv%nkpts
        do iorb=1,nvirt
           tt=real(e(iorb,ikpt,2)*orbsv%kwgts(ikpt),dp)
           if(msg)write(*,'(1x,i3,1x,1pe21.14)')iorb,tt
           gnrm=gnrm+tt
        end do
        if (nspin == 2) then
           do iorb=1,nvirt
              tt=real(e(iorb+orbsv%norbu,ikpt,2)*orbsv%kwgts(ikpt),dp)
              if(msg)write(*,'(1x,i3,1x,1pe21.14)')iorb,tt
              gnrm=gnrm+tt
           end do
        end if
     end do
     gnrm=sqrt(gnrm/real(orbsv%norb,dp))

     if(msg)write(*,'(1x,a,2(1x,1pe21.14))')"gnrm of all ",gnrm

     if (iproc==0)write(*,'(1x,a)',advance='no')'Preconditioning...'

     call timing(iproc,'Davidson      ','OF')

     !retranspose the gradient g 
     call untranspose_v(iproc,nproc,orbsv,lr%wfd,commsv,g,work=psiw)

     ! Here the gradients norm could be calculated in the direct form instead,
     ! as it is done in hpsiortho before preconditioning. 
     ! However, this should not make a difference and is not really simpler 

     call timing(iproc,'Precondition  ','ON')

     !we fill the values of the eval for the orbitals used in the preconditioner
     do ikpt=1,orbsv%nkpts
        do iorb=1,orbsv%norb
           orbsv%eval(iorb+(ikpt-1)*orbsv%norb)=e(iorb,ikpt,1)
        end do
     end do
     !we use for preconditioning the eval from the lowest value of the KS wavefunctions
     !WARNING this pointer association is nor coherent: 
     !        it will give errors if orbs%norb < orbsv%norb
!     orbsv%eval=>orbs%eval
!     if (orbs%norb < orbsv%norb) then
!        write(*,*)'ERROR: too many virtual orbitals demanded, cannot proceed.'
!        write(*,*)'       change the eigenvalues array to fix this problem'
!        stop
!     end if

     call preconditionall(iproc,nproc,orbsv,lr,hx,hy,hz,in%ncong,g,gnrm_fake)

     call timing(iproc,'Precondition  ','OF')
     if (iproc==0)write(*,'(1x,a)')'done.'

     if(iproc==0)write(*,'(1x,a)',advance="no")&
                 "Orthogonality of preconditioned gradients to occupied psi..."

     !transpose  g 
     call transpose_v(iproc,nproc,orbsv,lr%wfd,commsv,g,work=psiw)

     if (occorbs) then
        !project g such that they are orthogonal to all occupied psi
        call orthon_virt_occup(iproc,nproc,orbs,orbsv,comms,commsv,psi,g,msg)
     end if
     !retranspose the gradient g
     call untranspose_v(iproc,nproc,orbsv,lr%wfd,commsv,g,work=psiw)

     if(iproc==0)write(*,'(1x,a)')"done."

     allocate(hg(orbsv%npsidim+ndebug),stat=i_stat)
     call memocc(i_stat,hg,'hg',subname)

     call HamiltonianApplication(iproc,nproc,at,orbsv,hx,hy,hz,rxyz,&
          nlpspd,proj,lr,ngatherarr,n1i*n2i*n3p,&
          rhopot,g,hg,ekin_sum,epot_sum,eexctX,eproj_sum,in%nspin,GPU,&
          pkernel,orbs,psirocc) ! optional argument

     !transpose  g and hg
     call transpose_v(iproc,nproc,orbsv,lr%wfd,commsv,g,work=psiw)
     call transpose_v(iproc,nproc,orbsv,lr%wfd,commsv,hg,work=psiw)

     call timing(iproc,'Davidson      ','ON')
     if(iproc==0)write(*,'(1x,a)',advance="no")"done."


     if(msg)write(*,'(1x,a)')"Norm of all preconditioned gradients"

     call razero(orbsv%norb*orbsv%nkpts,e(1,1,2))

     ispsi=1
     do ikptp=1,orbsv%nkptsp
        ikpt=orbsv%iskpts+ikptp!orbsv%ikptsp(ikptp)
        nvctrp=commsv%nvctr_par(iproc,ikptp)
        if (nvctrp == 0) cycle

        nspinor=orbsv%nspinor
        do iorb=1,orbsv%norb
           e(iorb,ikpt,2)=nrm2(nvctrp*nspinor,g(ispsi+nvctrp*nspinor*(iorb-1)),1)**2
        end do
        ispsi=ispsi+nvctrp*orbsv%norb*orbsv%nspinor
     end do
     
     if(nproc > 1)then
        !sum up the contributions of nproc sets with nvctrp wavelet coefficients each
        call mpiallred(e(1,1,2),orbsv%norb*orbsv%nkpts,MPI_SUM,MPI_COMM_WORLD,ierr)
     end if

     gnrm=0.0_dp
     do ikpt=1,orbsv%nkpts
        do iorb=1,orbsv%norb
           tt=real(e(iorb,ikpt,2)*orbsv%kwgts(ikpt),dp)
           if(msg)write(*,'(1x,i3,1x,1pe21.14)')iorb,tt
           gnrm=gnrm+tt
        end do
     end do
     gnrm=sqrt(gnrm/real(orbsv%norb,dp))

     if(msg)write(*,'(1x,a,2(1x,1pe21.14))')"gnrm of all",gnrm

     if(iproc==0)write(*,'(1x,a)',advance="no")"Expanding subspace matrices..."

     !                 <vi | hvj>      <vi | hgj-n>                   <vi | vj>      <vi | gj-n>
     ! hamovr(i,j,1)=                               ;  hamovr(i,j,2)=  
     !                 <gi-n | hvj>  <gi-n | hgj-n>                   <gi-n | vj>  <gi-n | gj-n>
     !put to zero all the k-points which are not needed
     call razero(8*ndimovrlp(nspin,orbsv%nkpts),hamovr)


     ! store upper triangular part of these matrices only
     ! therefore, element (iorb+nvirte,jorb) is transposed to (j,nvirt+iorb)
     ispsi=1
     do ikptp=1,orbsv%nkptsp
        ikpt=orbsv%iskpts+ikptp!orbsv%ikptsp(ikptp)

        do ispin=1,nspin

           call orbitals_and_components(iproc,ikptp,ispin,orbsv,commsv,&
                nvctrp,norb,norbs,ncomp,nspinor)
           if (nvctrp == 0) cycle
           if (nspinor > 1) then
              ncplx=2
           else
              ncplx=1
           end if
!print *,iproc,ikpt,ispin,norb,nspinor,ncplx,nvctrp,8*ndimovrlp(ispin,ikpt-1)+1,8*ndimovrlp(nspin,orbsv%nkpts)
           call Davidson_subspace_hamovr(norb,nspinor,ncplx,nvctrp,&
                hamovr(8*ndimovrlp(ispin,ikpt-1)+1),&
                v(ispsi),g(ispsi),hv(ispsi),hg(ispsi))

           ispsi=ispsi+nvctrp*norb*nspinor
        end do
     end do

     i_all=-product(shape(hg))*kind(hg)
     deallocate(hg,stat=i_stat)
     call memocc(i_stat,i_all,'hg',subname)

     if(nproc > 1)then
        !sum up the contributions of nproc sets with nvctrp wavelet coefficients each
        call mpiallred(hamovr(1),8*ndimovrlp(nspin,orbsv%nkpts),&
             MPI_SUM,MPI_COMM_WORLD,ierr)
     end if

     if(iproc==0)write(*,'(1x,a)')"done."

     ispsi=1
     do ikptp=1,orbsv%nkptsp
        ikpt=orbsv%iskpts+ikptp!orbsv%ikptsp(ikptp)

        if(msg .or. (iproc==0 .and. ikpt == 0)) write(*,'(1x,a)',advance='no')"Diagonalization..."

        do ispin=1,nspin
           call orbitals_and_components(iproc,ikptp,ispin,orbsv,commsv,&
                nvctrp,norb,norbs,ncomp,nspinor)
           if (nvctrp == 0) cycle

           if (nspinor /= 1) then
              allocate(work_rp(6*norb+1+ndebug),stat=i_stat)
              call memocc(i_stat,work_rp,'work_rp',subname)
           end if

           ish1=8*ndimovrlp(ispin,ikpt-1)+1
           ish2=8*ndimovrlp(ispin,ikpt-1)+4*norbs*norb+1

           if(msg)then
              write(*,*)"subspace matrices, upper triangular (diagonal elements first)"
              write(*,'(1x)')
              write(*,*)"subspace H "
              do iorb=1,2*norbs
                 write(*,'(100(1x,1pe12.5))')&
                      (hamovr(ish1-1+iorb+(jorb-1)*2*norbs),jorb=1,2*norb)
                 write(*,*)
              end do
              write(*,*)"subspace S"
              write(*,*)
              do iorb=1,2*norbs
                 write(*,'(100(1x,1pe12.5))')&
                      (hamovr(ish2-1+iorb+(jorb-1)*2*norbs),jorb=1,2*norb)
                 write(*,*)
              end do
           end if

           if (nspinor == 1) then
              call sygv(1,'V','U',2*norb,hamovr(ish1),2*norb,hamovr(ish2),2*norb,&
                   ew(1),work(1),nwork,i_stat)! Lapack GEVP
              if (i_stat /= 0) write(*,*) &
                   'Error in SYGV on process ',iproc,', infocode ', i_stat

           else
              call hegv(1,'V','U',2*norb,hamovr(ish1),2*norb,hamovr(ish2),2*norb,&
                   ew(1),work(1),nwork,work_rp(1),i_stat)! Lapack GEVP

              if (i_stat /= 0) write(*,*) &
                   'Error in HEGV on process ',iproc,', infocode ', i_stat

           end if

           if (ispin==1) ise=0
           do iorb=1,norb
              e(iorb+ise,ikpt,1)=ew(iorb)
              e(iorb+ise,ikpt,2)=ew(iorb+norb)
           end do
           ise=norb

           if (nspinor /= 1) then
              i_all=-product(shape(work_rp))*kind(work_rp)
              deallocate(work_rp,stat=i_stat)
              call memocc(i_stat,i_all,'work_rp',subname)
           end if

           if(msg)then
              write(*,'(1x,a)')'    e(update)           e(not used)'
              do iorb=1,orbsv%norb
                 write(*,'(1x,i3,2(1pe21.14))')iorb, (e(iorb,ikpt,j),j=1,2)
              end do
              write(*,*)
              write(*,*)"and the eigenvectors are"
              write(*,*)
              do iorb=1,2*norbs
                 write(*,'(100(1x,1pe12.5))')&
                      (hamovr(ish1-1+iorb+(jorb-1)*2*norbs),jorb=1,2*norb)
                 write(*,*)
              end do
           end if


           if(msg .or. (iproc==0 .and. ikpt == 0))write(*,'(1x,a)',advance="no")&
                "done. Update v with eigenvectors..."

!!$     !Update v, that is the wavefunction, using the eigenvectors stored in hamovr(:,:,1)
!!$     !Lets say we have 4 quarters top/bottom left/right, then
!!$     !v = matmul(v, hamovr(topleft)  ) + matmul(g, hamovr(bottomleft)  )     needed    
!!$     !g=  matmul(v, hamovr(topright) ) + matmul(g, hamovr(bottomright) ) not needed
!!$     !use hv as work arrray
!!$
!!$     ispsi=1
!!$     do ikptp=1,orbsv%nkptsp
!!$        ikpt=orbsv%iskpts+ikptp!orbsv%ikptsp(ikptp)
!!$        nvctrp=commsv%nvctr_par(iproc,ikptp)
!!$
!!$        do jorb=1,orbsv%norb! v to update
!!$           call razero(nvctrp,hv(ispsi+nvctrp*(jorb-1)))
!!$           do iorb=1,orbsv%norb ! sum over v and g
!!$              tt=hamovr(iorb,jorb,ikpt,1)
!!$              call axpy(nvctrp,tt,v(ispsi+nvctrp*(iorb-1)),1,&
!!$                   hv(ispsi+nvctrp*(jorb-1)),1)
!!$              tt=hamovr(iorb+orbsv%norb,jorb,ikpt,1)
!!$              call axpy(nvctrp,tt,g(ispsi+nvctrp*(iorb-1)),1,&
!!$                   hv(ispsi+nvctrp*(jorb-1)),1)
!!$           enddo
!!$        enddo
!!$
!!$        call dcopy(nvctrp*orbsv%norb,hv(ispsi),1,v(ispsi),1)
!!$
!!$        ispsi=ispsi+nvctrp*orbsv%norb*orbsv%nspinor
!!$     end do

           call update_psivirt(norb,nspinor,ncplx,nvctrp,&
                hamovr(ish1),v(ispsi),g(ispsi),hv(ispsi))

           ispsi=ispsi+nvctrp*norb*nspinor

           if(msg .or. (iproc==0 .and. ikpt == 0)) then
              if (nspin ==1) then
                 write(*,'(1x,a)')'done. The refined eigenvalues are'
                 do iorb=1,nvirt
                    write(*,'(1x,i3,2(1pe21.14))')iorb,(e(iorb,ikpt,j),j=1,2)
                 end do
              else if (ispin == 2) then
                 write(*,'(1x,a)')'done. The refined eigenvalues are'
                 do iorb=1,nvirt
                    write(*,'(1x,i3,4(1pe21.14))')&
                         iorb,(e(iorb,ikpt,j),j=1,2),(e(iorb+orbsv%norbu,ikpt,j),j=1,2)
                 end do
              end if
           end if
           
        end do
        
     end do

     i_all=-product(shape(g))*kind(g)
     deallocate(g,stat=i_stat)
     call memocc(i_stat,i_all,'g',subname)

     if(iproc==0)write(*,'(1x,a)')"done."
     if(iproc==0)write(*,'(1x,a)',advance="no")"Orthogonality to occupied psi..."
     !project v such that they are orthogonal to all occupied psi
     !Orthogonalize before and afterwards.

     call timing(iproc,'Davidson      ','OF')

     !these routines should work both in parallel or in serial
     call orthogonalize(iproc,nproc,orbsv,commsv,lr%wfd,v)

     if (occorbs) then
        call orthon_virt_occup(iproc,nproc,orbs,orbsv,comms,commsv,psi,v,msg)
     !and orthonormalize them using "gram schmidt"  (should conserve orthogonality to psi)
        call orthogonalize(iproc,nproc,orbsv,commsv,lr%wfd,v)
     end if

     !retranspose v
     call untranspose_v(iproc,nproc,orbsv,lr%wfd,commsv,v,work=psiw)

     ! Hamilton application on v
     if(iproc==0)write(*,'(1x,a)',advance="no")"done."

     call HamiltonianApplication(iproc,nproc,at,orbsv,hx,hy,hz,rxyz,&
          nlpspd,proj,lr,ngatherarr,n1i*n2i*n3p,&
          rhopot,v,hv,ekin_sum,epot_sum,eexctX,eproj_sum,in%nspin,GPU,&
          pkernel,orbs,psirocc) !optional arguments

     !transpose  v and hv
     call transpose_v(iproc,nproc,orbsv,lr%wfd,commsv,v,work=psiw)
     call transpose_v(iproc,nproc,orbsv,lr%wfd,commsv,hv,work=psiw)

     if(iproc==0 .and. verbose > 1) write(*,'(1x,a)')"done. "
     call timing(iproc,'Davidson      ','ON')
     iter=iter+1
     if(iter>in%itermax+100)then !an input variable should be put
        if(iproc==0)write(*,'(1x,a)')&
             'No convergence within the allowed number of minimization steps (itermax + 100)'
        exit davidson_loop
     end if

  end do davidson_loop

  i_all=-product(shape(ndimovrlp))*kind(ndimovrlp)
  deallocate(ndimovrlp,stat=i_stat)
  call memocc(i_stat,i_all,'ndimovrlp',subname)

  i_all=-product(shape(hamovr))*kind(hamovr)
  deallocate(hamovr,stat=i_stat)
  call memocc(i_stat,i_all,'hamovr',subname)

  i_all=-product(shape(work))*kind(work)
  deallocate(work,stat=i_stat)
  call memocc(i_stat,i_all,'work',subname)

  i_all=-product(shape(ew))*kind(ew)
  deallocate(ew,stat=i_stat)
  call memocc(i_stat,i_all,'ew',subname)

  !deallocate real array of wavefunctions
  if(exctX)then
     i_all=-product(shape(psirocc))*kind(psirocc)
     deallocate(psirocc,stat=i_stat)
     call memocc(i_stat,i_all,'psirocc',subname)
  end if


  if(iter <=in%itermax) then
     if(iproc==0)write(*,'(1x,a,i3,a)')&
          'Davidsons method: Convergence after ',iter-1,' iterations.'
  end if
  !finalize: Retranspose, deallocate

  ! Send all eigenvalues to all procs.
  call broadcast_kpt_objects(nproc, orbsv%nkpts, orbsv%norb, e(1,1,1), orbsv%ikptproc)

  if(iproc==0)then
     if (nspin==1) then
        write(*,'(1x,a)')'Complete list of energy eigenvalues'
        do ikpt=1,orbsv%nkpts
           write(*,"(1x,A,I3.3,A,3F12.6)") &
                & "Kpt #", ikpt, " BZ coord. = ", orbsv%kpts(:, ikpt)
           do iorb=1,orbs%norb
              write(*,'(1x,a,i4,a,1x,1pe21.14)') 'e_occupied(',iorb,')=',&
                   orbs%eval(iorb+(ikpt-1)*orbs%norb)
           end do
           write(*,'(1x,a,1pe21.14,a,0pf8.4,a)')&
                'HOMO LUMO gap   =',e(1,ikpt,1)-orbs%eval(orbs%norb+(ikpt-1)*orbs%norb),&
                ' (',ha2ev*(e(1,ikpt,1)-orbs%eval(orbs%norb+(ikpt-1)*orbs%norb)),&
                ' eV)'
           do iorb=1,orbsv%norb
              write(*,'(1x,a,i4,a,1x,1pe21.14)') 'e_virtual(',iorb,')=',e(iorb,ikpt,1)
           end do
        end do
     else
        do ikpt=1,orbsv%nkpts
           write(*,'(1x,a)')'Complete list of energy eigenvalues'
           do iorb=1,min(orbs%norbu,orbs%norbd)
              jorb=orbs%norbu+iorb
              write(*,'(1x,a,i4,a,1x,1pe21.14,14x,a,i4,a,1x,1pe21.14)') &
                   'e_occ(',iorb,',u)=',orbs%eval(iorb+(ikpt-1)*orbs%norb),&
                   'e_occ(',iorb,',d)=',orbs%eval(jorb+(ikpt-1)*orbs%norb)
           end do
           if (orbs%norbu > orbs%norbd) then
              do iorb=orbs%norbd+1,orbs%norbu
                 write(*,'(1x,a,i4,a,1x,1pe21.14)') &
                      'e_occ(',iorb,',u)=',orbs%eval(iorb+(ikpt-1)*orbs%norb)
              end do
           else if (orbs%norbd > orbs%norbu) then
              do iorb=2*orbs%norbu+1,orbs%norbu+orbs%norbd
                 write(*,'(50x,a,i4,a,1x,1pe21.14)') &
                      'e_occ(',iorb-orbs%norbu,',d)=',orbs%eval(iorb+(ikpt-1)*orbs%norb)
              end do
           end if
           write(*,'(1x,a,1x,1pe21.14,a,0pf8.4,a,a,1x,1pe21.14,a,0pf8.4,a)') &
                'HOMO LUMO gap, u =',&
                e(1,ikpt,1)-orbs%eval(orbs%norbu+(ikpt-1)*orbs%norb),&
                ' (',ha2ev*(e(1,ikpt,1)-orbs%eval(orbs%norbu+(ikpt-1)*orbs%norb)),&
                ' eV)',&
                ',d =',e(orbsv%norbu+1,ikpt,1)-orbs%eval(orbs%norb+(ikpt-1)*orbs%norb),&
                ' (',&
                ha2ev*(e(orbsv%norbu+1,ikpt,1)-orbs%eval(orbs%norb+(ikpt-1)*orbs%norb)),&
                ' eV)'
           do iorb=1,min(orbsv%norbu,orbsv%norbd)
              jorb=orbsv%norbu+iorb
              write(*,'(1x,a,i4,a,1x,1pe21.14,14x,a,i4,a,1x,1pe21.14)') &
                   'e_vrt(',iorb,',u)=',e(iorb,ikpt,1),&
                   'e_vrt(',iorb,',d)=',e(jorb,ikpt,1)
           end do
           if (orbsv%norbu > orbsv%norbd) then
              do iorb=orbsv%norbd+1,orbsv%norbu
                 write(*,'(1x,a,i4,a,1x,1pe21.14)') &
                      'e_vrt(',iorb,',u)=',e(iorb,ikpt,1)
              end do
           else if (orbsv%norbd > orbsv%norbu) then
              do iorb=2*orbsv%norbu+1,orbsv%norbu+orbsv%norbd
                 write(*,'(50x,a,i4,a,1x,1pe21.14)') &
                      'e_vrt(',iorb-orbsv%norbu,',d)=',e(iorb,ikpt,1)
              end do
           end if
        end do
     end if
     end if

  call timing(iproc,'Davidson      ','OF')

  !retranspose v and psi
  call untranspose_v(iproc,nproc,orbsv,lr%wfd,commsv,v,work=psiw)

  !resize work array before final transposition
  if(nproc > 1)then
     i_all=-product(shape(psiw))*kind(psiw)
     deallocate(psiw,stat=i_stat)
     call memocc(i_stat,i_all,'psiw',subname)

     allocate(psiw(orbs%npsidim+ndebug),stat=i_stat)
     call memocc(i_stat,psiw,'psiw',subname)
  end if

  call untranspose_v(iproc,nproc,orbs,lr%wfd,comms,psi,work=psiw)

  if(nproc > 1) then
     i_all=-product(shape(psiw))*kind(psiw)
     deallocate(psiw,stat=i_stat)
     call memocc(i_stat,i_all,'psiw',subname)
  end if

  i_all=-product(shape(hv))*kind(hv)
  deallocate(hv,stat=i_stat)
  call memocc(i_stat,i_all,'hv',subname)

  call deallocate_comms(commsv,subname)

  ! PLOTTING

  !plot the converged wavefunctions in the different orbitals.
  !nplot is the requested total of orbitals to plot, where
  !states near the HOMO/LUMO gap are given higher priority.
  !Occupied orbitals are only plotted when nplot>nvirt,
  !otherwise a comment is given in the out file.

  if(abs(in%nplot)>orbs%norb+nvirt)then
     if(iproc==0)write(*,'(1x,A,i3)')&
          "WARNING: More plots requested than orbitals calculated." 
  end if

  !add a modulo operator to get rid of the particular k-point
  do iorb=1,orbsv%norbp!requested: nvirt of nvirte orbitals
     if(modulo(iorb+orbsv%isorb-1,orbsv%norb)+1 > abs(in%nplot))then
        if(iproc == 0 .and. abs(in%nplot) > 0)write(*,'(A)')&
             'WARNING: No plots of occupied orbitals requested.'
        exit 
     end if
     ind=1+(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*(iorb-1)
     !plot the orbital and the density
     write(orbname,'(A,i3.3)')'virtual',iorb+orbsv%isorb
     write(denname,'(A,i3.3)')'denvirt',iorb+orbsv%isorb
     write(comment,'(1pe10.3)')e(modulo(iorb+orbsv%isorb-1,orbsv%norb)+1,orbsv%iokpt(iorb),1)
     !choose the way of plotting the wavefunctions
     if (in%nplot > 0) then
        call plot_wf('POT',orbname,1,at,lr,hx,hy,hz,rxyz,v(ind:),comment)
        call plot_wf('POT',denname,2,at,lr,hx,hy,hz,rxyz,v(ind:),comment)
     else if (in%nplot < 0) then
        call plot_wf('CUBE',orbname,1,at,lr,hx,hy,hz,rxyz,v(ind:),comment)
        call plot_wf('CUBE',denname,2,at,lr,hx,hy,hz,rxyz,v(ind:),comment)
     end if
  end do

  do iorb=orbs%norbp,1,-1 ! sweep over highest occupied orbitals
     if(modulo(orbs%norb-iorb-orbs%isorb-1,orbs%norb)+1 > abs(in%nplot)) then
        exit! we have written nplot pot files
     end if
     !address
     ind=1+(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*(iorb-1)
     write(orbname,'(A,i3.3)')'orbital',iorb+orbs%isorb
     write(denname,'(A,i3.3)')'densocc',iorb+orbs%isorb
     write(comment,'(1pe10.3)')orbs%eval(iorb+orbs%isorb)
     !choose the way of plotting the wavefunctions
     if (in%nplot > 0) then
        call plot_wf('POT',orbname,1,at,lr,hx,hy,hz,rxyz,psi(ind:),comment)
        call plot_wf('POT',denname,2,at,lr,hx,hy,hz,rxyz,psi(ind:),comment)
     else if (in%nplot < 0) then
        call plot_wf('CUBE',orbname,1,at,lr,hx,hy,hz,rxyz,psi(ind:),comment)
        call plot_wf('CUBE',denname,2,at,lr,hx,hy,hz,rxyz,psi(ind:),comment)
     end if

  end do
  ! END OF PLOTTING

  i_all=-product(shape(e))*kind(e)
  deallocate(e,stat=i_stat)
  call memocc(i_stat,i_all,'e',subname)

  if (GPUconv) then
     call free_gpu(GPU,orbsv%norbp)
  else if (OCLconv) then
     call free_gpu_OCL(GPU,orbsv%norbp)
  end if


END SUBROUTINE davidson
!!***


!!****f* BigDFT/Davidson_subspace_hamovr
!! FUNCTION
!!   Generate upper triangular matrix in the subspace of Davidson algorithm
!! SOURCE
!!
subroutine Davidson_subspace_hamovr(norb,nspinor,ncplx,nvctrp,hamovr,v,g,hv,hg)
  use module_base
  implicit none
  integer, intent(in) :: norb,nvctrp,nspinor,ncplx
  real(wp), dimension(nspinor*nvctrp*norb), intent(in) :: v,g,hv,hg
  real(wp), dimension(ncplx,2*norb,2*norb,2), intent(out) :: hamovr
  !local variables
  character(len=*), parameter :: subname='Davidson_subspace_hamovr'
  integer :: iorb,jorb,icplx,ncomp

  if (nspinor == 4) then
     ncomp=2
  else
     ncomp=1
  end if

  !                 <vi | hvj>      <vi | hgj-n>                   <vi | vj>      <vi | gj-n>
  ! hamovr(i,j,1)=                               ;  hamovr(i,j,2)=  
  !                 <gi-n | hvj>  <gi-n | hgj-n>                   <gi-n | vj>  <gi-n | gj-n>

  !  do iorb=1,norb
  !     do jorb=iorb,norb
  !        hamovr(1,iorb,jorb,1)=               dot(nvctrp,v(1,iorb),1,hv(1,jorb),1)
  !        hamovr(1,jorb,iorb+norb,1)=        dot(nvctrp,g(1,iorb),1,hv(1,jorb),1)
  !        hamovr(1,iorb,jorb+norb,1)=        dot(nvctrp,v(1,iorb),1,hg(1,jorb),1)
  !        hamovr(1,iorb+norb,jorb+norb,1)= dot(nvctrp,g(1,iorb),1,hg(1,jorb),1)
  !               
  !        hamovr(1,iorb,jorb,2)=               dot(nvctrp,v(1,iorb),1, v(1,jorb),1)
  !        hamovr(1,jorb,iorb+norb,2)=       dot(nvctrp,g(1,iorb),1, v(1,jorb),1)
  !        hamovr(1,iorb,jorb+norb,2)=        dot(nvctrp,v(1,iorb),1, g(1,jorb),1)
  !        hamovr(1,iorb+norb,jorb+norb,2)= dot(nvctrp,g(1,iorb),1, g(1,jorb),1)
  !     enddo
  !  enddo
  !
  !use lapack operations to generalise the calculation to different nspinor

  !4 gemm + 2 dsyrk operations

  !<vi | hvj> 
  if(nspinor==1) then
     call gemm('T','N',norb,norb,nvctrp,1.0_wp,v(1),&
          max(1,nvctrp),hv(1),max(1,nvctrp),0.0_wp,&
          hamovr(1,1,1,1),2*norb)
  else
     call c_gemm('C','N',norb,norb,ncomp*nvctrp,(1.0_wp,0.0_wp),v(1),&
          max(1,ncomp*nvctrp), &
          hv(1),max(1,ncomp*nvctrp),(0.0_wp,0.0_wp),&
          hamovr(1,1,1,1),2*norb)
  end if

  !<gi | hvj> 
  if(nspinor==1) then
     call gemm('T','N',norb,norb,nvctrp,1.0_wp,g(1),&
          max(1,nvctrp),hv(1),max(1,nvctrp),0.0_wp,&
          hamovr(1,norb+1,1,1),2*norb)
  else
     call c_gemm('C','N',norb,norb,ncomp*nvctrp,(1.0_wp,0.0_wp),g(1),&
          max(1,ncomp*nvctrp), &
          hv(1),max(1,ncomp*nvctrp),(0.0_wp,0.0_wp),&
          hamovr(1,norb+1,1,1),2*norb)
  end if

  !<gi | hgj>
  if(nspinor==1) then
     call gemm('T','N',norb,norb,nvctrp,1.0_wp,g(1),&
          max(1,nvctrp),hg(1),max(1,nvctrp),0.0_wp,&
          hamovr(1,norb+1,norb+1,1),2*norb)
  else
     call c_gemm('C','N',norb,norb,ncomp*nvctrp,(1.0_wp,0.0_wp),g(1),&
          max(1,ncomp*nvctrp), &
          hg(1),max(1,ncomp*nvctrp),(0.0_wp,0.0_wp),&
          hamovr(1,norb+1,norb+1,1),2*norb)
  end if

  !<vi | vj> 
  if(nspinor==1) then
     call syrk('U','T',norb,nvctrp,1.0_wp,v(1),max(1,nvctrp),&
          0.0_wp,hamovr(1,1,1,2),2*norb)
  else
     call herk('U','C',norb,ncomp*nvctrp,1.0_wp,v(1),max(1,ncomp*nvctrp),&
          0.0_wp,hamovr(1,1,1,2),2*norb)
  end if

  !<gi | vj> => hsub(:,:,:,5)
  if(nspinor==1) then
     call gemm('T','N',norb,norb,nvctrp,1.0_wp,g(1),&
          max(1,nvctrp),v(1),max(1,nvctrp),0.0_wp,&
          hamovr(1,norb+1,1,2),2*norb)
  else
     call c_gemm('C','N',norb,norb,ncomp*nvctrp,(1.0_wp,0.0_wp),g(1),&
          max(1,ncomp*nvctrp), &
          v(1),max(1,ncomp*nvctrp),(0.0_wp,0.0_wp),&
          hamovr(1,norb+1,1,2),2*norb)
  end if

  !<gi | gj>
  if(nspinor==1) then
     call syrk('U','T',norb,nvctrp,1.0_wp,g(1),max(1,nvctrp),&
          0.0_wp,hamovr(1,norb+1,norb+1,2),2*norb)
  else
     call herk('U','C',norb,ncomp*nvctrp,1.0_wp,g(1),max(1,ncomp*nvctrp),&
          0.0_wp,hamovr(1,norb+1,norb+1,2),2*norb)
  end if


  !fill hamovr (Upper triangular)
  do jorb=1,norb
     do iorb=1,norb
        do icplx=1,ncplx
           hamovr(icplx,iorb,jorb+norb,1) = (-1)**(icplx-1)*hamovr(icplx,jorb+norb,iorb,1)  
           hamovr(icplx,iorb,jorb+norb,2) = (-1)**(icplx-1)*hamovr(icplx,jorb+norb,iorb,2)  
        end do
     enddo
  enddo

END SUBROUTINE Davidson_subspace_hamovr
!!***


subroutine update_psivirt(norb,nspinor,ncplx,nvctrp,hamovr,v,g,work)
  use module_base
  implicit none
  integer, intent(in) :: norb,nvctrp,nspinor,ncplx
  real(wp), dimension(nspinor*nvctrp*norb), intent(in) :: g
  real(wp), dimension(nspinor*nvctrp*norb), intent(inout) :: v
  real(wp), dimension(nspinor*nvctrp*norb), intent(inout) :: work
  real(wp), dimension(ncplx,2*norb,2*norb), intent(out) :: hamovr
  !local variables
  character(len=*), parameter :: subname='update_psivirt'
  integer :: ncomp

  if (nspinor == 4) then
     ncomp=2
  else
     ncomp=1
  end if

  !Update v, that is the wavefunction, using eigenvectors stored in hamovr(:,:,1)
  !Lets say we have 4 quarters top/bottom left/right, then
  !v = matmul(v, hamovr(topleft)  ) + matmul(g, hamovr(bottomleft)  )  needed    
  !g=  matmul(v, hamovr(topright) ) + matmul(g, hamovr(bottomright) ) not needed
  !use hv as work arrray

  !Note: The previous data layout allowed level 3 BLAS
  !call DGEMM('N','N',nvctrp,nvirte,n2virt,1.d0,v(1,1),nvctrp,hamovr(1,1,1),n2virt,0.d0,hv(1,1),nvctrp)
  !    dimensions    =m      =n   =k          m,k        k,n                   m,n             
  !call DCOPY(nvctrp*nvirte,hv(1,1),1,v(1,1),1)

  if(nspinor==1) then
     call gemm('N','N',nvctrp,norb,norb,1.0_wp,v(1),&
          max(1,nvctrp),hamovr(1,1,1),max(1,2*norb),0.0_wp,&
          work(1),nvctrp)
     call gemm('N','N',nvctrp,norb,norb,1.0_wp,g(1),&
          max(1,nvctrp),hamovr(1,norb+1,1),max(1,2*norb),1.0_wp,&
          work(1),nvctrp)

  else
     call c_gemm('N','N',ncomp*nvctrp,norb,norb,(1.0_wp,0.0_wp),v(1),&
          max(1,ncomp*nvctrp),hamovr(1,1,1),max(1,2*norb),(0.0_wp,0.0_wp),&
          work(1),ncomp*nvctrp)
     call c_gemm('N','N',ncomp*nvctrp,norb,norb,(1.0_wp,0.0_wp),g(1),&
          max(1,ncomp*nvctrp),hamovr(1,norb+1,1),max(1,2*norb),(1.0_wp,0.0_wp),&
          work(1),ncomp*nvctrp)
  end if

  call dcopy(nspinor*nvctrp*norb,work(1),1,v(1),1)

END SUBROUTINE update_psivirt

subroutine psivirt_from_gaussians(iproc,nproc,at,orbs,lr,comms,rxyz,hx,hy,hz,nspin,psivirt)
  use module_base
  use module_types
  use module_interfaces
  implicit none
  integer, intent(in) :: iproc,nproc,nspin
  real(gp), intent(in) :: hx,hy,hz
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  type(locreg_descriptors), intent(in) :: lr
  type(communications_arrays), intent(in) :: comms
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(wp), dimension(orbs%npsidim), intent(out) :: psivirt
  !local variables
  character(len=*), parameter :: subname='psivirt_from_gaussians'
  integer :: iorb,icoeff,i_all,i_stat,jproc
  real(kind=4) :: tt
  real(wp), dimension(:,:), allocatable :: gaucoeffs
  type(gaussian_basis) :: G
  real(wp), dimension(:), pointer :: gbd_occ,psiw


  !initialise some coefficients in the gaussian basis
  !nullify the G%rxyz pointer
  nullify(G%rxyz)
  !extract the gaussian basis from the pseudowavefunctions
  !use a better basis than the input guess
  call gaussian_pswf_basis(31,iproc,nspin,at,rxyz,G,gbd_occ)

  allocate(gaucoeffs(G%ncoeff,orbs%norbp*orbs%nspinor+ndebug),stat=i_stat)
  call memocc(i_stat,gaucoeffs,'gaucoeffs',subname)

  !fill randomly the gaussian coefficients for the orbitals considered
  do iorb=1,orbs%norbp*orbs%nspinor
     do icoeff=1,G%ncoeff
        !be sure to call always a different random number
        do jproc=0,iproc-1
           call random_number(tt)
        end do
        call random_number(tt)
        gaucoeffs(icoeff,iorb)=real(tt,wp)
        do jproc=iproc+1,nproc-1
           call random_number(tt)
        end do
     end do
  end do

  !othogonalise the gaussian basis (wrong with k-points)
  !call gaussian_orthogonality(iproc,nproc,norb,norbp,G,coeffs)

  call gaussians_to_wavelets_new(iproc,nproc,lr,orbs,hx,hy,hz,G,&
       gaucoeffs,psivirt)

  !deallocate the gaussian basis descriptors
  call deallocate_gwf(G,subname)

  !deallocate gaussian array
  i_all=-product(shape(gaucoeffs))*kind(gaucoeffs)
  deallocate(gaucoeffs,stat=i_stat)
  call memocc(i_stat,i_all,'gaucoeffs',subname)
  i_all=-product(shape(gbd_occ))*kind(gbd_occ)
  deallocate(gbd_occ,stat=i_stat)
  call memocc(i_stat,i_all,'gbd_occ',subname)


  !transpose v
  if(nproc > 1)then
     !reallocate the work array with the good size
     allocate(psiw(orbs%npsidim+ndebug),stat=i_stat)
     call memocc(i_stat,psiw,'psiw',subname)
  end if

  !transpose the wavefunction in wavelet basis
  call transpose_v(iproc,nproc,orbs,lr%wfd,comms,psivirt,work=psiw)

  if(nproc > 1)then
     i_all=-product(shape(psiw))*kind(psiw)
     deallocate(psiw,stat=i_stat)
     call memocc(i_stat,i_all,'psiw',subname)
  end if

  
END SUBROUTINE psivirt_from_gaussians

