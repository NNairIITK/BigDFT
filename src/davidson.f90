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
!!    Copyright (C) 2007-2009 CEA, UNIBAS
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! SOURCE
!!
subroutine davidson(iproc,nproc,n1i,n2i,n3i,in,at,cpmult,fpmult,radii_cf,&
     orbs,orbsv,nvirt,lr,comms,&
     hx,hy,hz,rxyz,rhopot,i3xcsh,n3p,nlpspd,proj,pkernel,psi,v,ngatherarr)
  use module_base
  use module_types
  use module_interfaces, except_this_one => davidson
  implicit none
  integer, intent(in) :: iproc,nproc,n1i,n2i,n3i
  integer, intent(in) :: i3xcsh
  integer, intent(in) :: nvirt,n3p
  type(input_variables), intent(in) :: in
  type(atoms_data), intent(in) :: at
  type(nonlocal_psp_descriptors), intent(in) :: nlpspd
  type(locreg_descriptors), intent(in) :: lr 
  type(orbitals_data), intent(in) :: orbs
  type(communications_arrays), intent(in) :: comms
  real(gp), dimension(at%ntypes,3), intent(in) :: radii_cf  
  real(gp), intent(in) :: hx,hy,hz,cpmult,fpmult
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
  real(dp), dimension(*), intent(in) :: pkernel,rhopot
  type(orbitals_data), intent(inout) :: orbsv
  !this is a Fortran 95 standard, should be avoided (it is a pity IMHO)
  !real(kind=8), dimension(:,:,:,:), allocatable :: rhopot 
  real(wp), dimension(:), pointer :: psi,v!=psivirt(nvctrp,nvirtep*nproc) 
                        !v, that is psivirt, is transposed on input and direct on output
  !local variables
  character(len=*), parameter :: subname='davidson'
  character(len=10) :: comment
  character(len=11) :: orbname
  logical :: msg !extended output
  integer :: n2virt,n2virtp,ierr,i_stat,i_all,iorb,jorb,iter,nwork,ind,i1,i2!<-last 3 for debug
  integer :: ise,ish,jnd,j,ispsi,ikpt,ikptp,nvctrp
  real(kind=8) :: tt,gnrm,eks,eexcu,vexcu,epot_sum,ekin_sum,ehart,eproj_sum,etol,gnrm_fake
  type(communications_arrays) :: commsv
  type(GPU_pointers) :: GPU !added for interface compatibility, not working here
  real(wp), dimension(:), allocatable :: work
  real(wp), dimension(:), allocatable :: hv,g,hg,ew
  real(wp), dimension(:,:,:,:), allocatable :: e,hamovr
  real(wp), dimension(:), pointer :: psiw

  !last index of e and hamovr are for mpi_allreduce. 
  !e (eigenvalues) is also used as 2 work arrays
  
  msg=verbose > 2!.false.! no extended output
  !msg =(iproc==0)!extended output

  if(iproc==0)write(*,'(1x,a)')"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  if(iproc==0)write(*,'(1x,a)')"Iterative subspace diagonalization of virtual orbitals."

  !if(msg)write(*,*)'shape(v)',shape(v),'size(v)',size(v)

  !dimensions for allocations of matrices
  if (nproc > 1) then
     ise=2
     ish=4
  else
     ise=1
     ish=2
  end if

  n2virt=2*orbsv%norb! the dimension of the subspace
  !to be modified for complex wavefunctions

  if (nproc == 1) then
   
  end if
 
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
  !call allocate_comms(nproc,orbsv,commsv,subname)
  call orbitals_communicators(iproc,nproc,lr,orbsv,commsv)  

  i_all=-product(shape(orbsv%norb_par))*kind(orbsv%norb_par)
  deallocate(orbsv%norb_par,stat=i_stat)
  call memocc(i_stat,i_all,'orbsv%norb_par',subname)


  if(iproc==0)write(*,'(1x,a)',advance="no")"Orthogonality to occupied psi..."
  !project v such that they are orthogonal to all occupied psi
  !Orthogonalize before and afterwards.

  !here nvirte=orbsv%norb
  !     nvirtep=orbsv%norbp

  !this is the same also in serial
  call orthogonalize(iproc,nproc,orbsv,commsv,lr%wfd,v)
  !call orthon_p(iproc,nproc,orbsv%norb,commsv%nvctr_par(iproc,1),lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,v,1)
  call orthon_virt_occup(iproc,nproc,orbs,orbsv,comms,commsv,psi,v,msg)

!!$  call orthoconvirt_p(iproc,nproc,orbs%norbu,orbsv%norb,comms%nvctr_par(iproc,1),psi,v,msg)
!!$  if(orbs%norbd > 0) then
!!$     call orthoconvirt_p(iproc,nproc,orbs%norbd,&
!!$       orbsv%norb,comms%nvctr_par(iproc,1),psi(1+comms%nvctr_par(iproc,1)*orbs%norbu),v,msg)
!!$  end if
  !and orthonormalize them using "gram schmidt"  (should conserve orthogonality to psi)
  call orthogonalize(iproc,nproc,orbsv,commsv,lr%wfd,v)
  !call orthon_p(iproc,nproc,orbsv%norb,commsv%nvctr_par(iproc,1),lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,v,1)

  !retranspose v
  if(nproc > 1)then
     !reallocate the work array with the good sizeXS
     allocate(psiw(orbsv%npsidim+ndebug),stat=i_stat)
     call memocc(i_stat,psiw,'psiw',subname)
  end if

  call untranspose_v(iproc,nproc,orbsv,lr%wfd,commsv,v,work=psiw)

  ! 1st Hamilton application on psivirt
  if(iproc==0)write(*,'(1x,a)',advance="no")"done. first "

  allocate(hv(orbsv%npsidim+ndebug),stat=i_stat)
  call memocc(i_stat,hv,'hv',subname)
  
  call HamiltonianApplication(iproc,nproc,at,orbsv,hx,hy,hz,rxyz,cpmult,fpmult,radii_cf,&
       nlpspd,proj,lr,ngatherarr,n1i*n2i*n3p,&
       rhopot(1+i3xcsh*n1i*n2i),v,hv,ekin_sum,epot_sum,eproj_sum,in%nspin,GPU)

  !if(iproc==0)write(*,'(1x,a)',advance="no")"done. Rayleigh quotients..."

  allocate(e(orbsv%norb,orbsv%nkpts,2,ise+ndebug),stat=i_stat)
  call memocc(i_stat,e,'e',subname)

  !transpose  v and hv
  call transpose_v(iproc,nproc,orbsv,lr%wfd,commsv,v,work=psiw)
  call transpose_v(iproc,nproc,orbsv,lr%wfd,commsv,hv,work=psiw)

  call timing(iproc,'Davidson      ','ON')
  !Timing excludes transposition, hamilton application and preconditioning
  call razero(orbsv%norb*2*orbsv%nkpts*ise,e)
  ! Rayleigh quotients.


  !probably this loop can be rewritten using GEMMs
  ispsi=1
  do ikptp=1,orbsv%nkptsp
     ikpt=ikptp+orbsv%iskpts
     nvctrp=commsv%nvctr_par(iproc,ikptp)
     do iorb=1,orbsv%norb ! temporary variables 
        e(iorb,ikpt,1,ise)= dot(nvctrp,v(ispsi+nvctrp*(iorb-1)),1,&
             hv(ispsi+nvctrp*(iorb-1)),1)          != <psi|H|psi> 
        e(iorb,ikpt,2,ise)= nrm2(nvctrp,v(ispsi+nvctrp*(iorb-1)),1)**2   != <psi|psi> 
     end do
     ispsi=ispsi+nvctrp*orbsv%norb*orbsv%nspinor
  end do

  if(nproc > 1)then
     !sum up the contributions of nproc sets with commsv%nvctr_par(iproc,1) wavelet coefficients each
     call MPI_ALLREDUCE(e(1,1,1,2),e(1,1,1,1),2*orbsv%norb*orbsv%nkpts,&
          mpidtypw,MPI_SUM,MPI_COMM_WORLD,ierr)
  end if



  if(iproc==0)write(*,'(1x,a)')"done."
  if(iproc==0)write(*,'(1x,a)')"     sqnorm                Rayleigh quotient"
 
  do ikpt=1,orbsv%nkpts
     do iorb=1,orbsv%norb
        !e(:,1,1) = <psi|H|psi> / <psi|psi>
        e(iorb,ikpt,1,1)=e(iorb,ikpt,1,1)/e(iorb,ikpt,2,1)
        if(iproc==0) write(*,'(1x,i3,2(1x,1pe21.14))')iorb, e(iorb,ikpt,2,1), e(iorb,ikpt,1,1)
     end do
  end do

!if(msg)then
!write(*,*)"******** transposed v,hv 1st elements"
!do iorb=1,10
!  write(*,*)v(iorb,1),hv(iorb,1)
!end do
!write(*,*)"**********"
!end if
  

  !itermax=... use the input variable instead
  do iter=1,in%itermax
     if(iproc==0)write(*,'(1x,a,i3)')&
     "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~iter",iter
     if(msg) write(*,'(1x,a)')"squared norm of the (nvirt) gradients"

     allocate(g(orbsv%npsidim+ndebug),stat=i_stat)
     call memocc(i_stat,g,'g',subname)

     call dcopy(orbsv%npsidim,hv,1,g,1)! don't overwrite hv

     call razero(orbsv%norb*orbsv%nkpts,e(1,1,2,ise))
     !also these operations are presumably GEMMs
     ispsi=1
     do ikptp=1,orbsv%nkptsp
        ikpt=ikptp+orbsv%iskpts
        nvctrp=commsv%nvctr_par(iproc,ikptp)

        do iorb=1,orbsv%norb
           call axpy(nvctrp,-e(iorb,ikpt,1,1),v(ispsi+nvctrp*(iorb-1)),1,&
                g(ispsi+nvctrp*(iorb-1)),1)
           !gradient = hv-e*v
           e(iorb,ikpt,2,ise)= nrm2(nvctrp,g(ispsi+nvctrp*(iorb-1)),1)**2
           !local contribution to the square norm
        end do

        ispsi=ispsi+nvctrp*orbsv%norb*orbsv%nspinor
     end do

     if(nproc > 1)then
        !sum up the contributions of nproc sets with nvctrp wavelet coefficients each
        call MPI_ALLREDUCE(e(1,1,2,2),e(1,1,2,1),orbsv%norb*orbsv%nkpts,&
             mpidtypw,MPI_SUM,MPI_COMM_WORLD,ierr)
     end if

     gnrm=0._dp
     do ikpt=1,orbsv%nkpts
        do iorb=1,nvirt
           tt=real(e(iorb,ikpt,2,1)*orbs%kwgts(ikpt),dp)
           if(msg)write(*,'(1x,i3,1x,1pe21.14)')iorb,tt
           gnrm=gnrm+tt
        end do
     end do
     gnrm=dsqrt(gnrm/real(orbsv%norb,dp))

     if(iproc == 0)write(*,'(1x,a,2(1x,1pe12.5))')&
          "|gradient|=gnrm and exit criterion ",gnrm,in%gnrm_cv
     if(gnrm < in%gnrm_cv) then
        i_all=-product(shape(g))*kind(g)
        deallocate(g,stat=i_stat)
        call memocc(i_stat,i_all,'g',subname)
        exit ! iteration loop
     end if
     call timing(iproc,'Davidson      ','OF')

     if(iproc==0)write(*,'(1x,a)',advance="no")&
          "Orthogonality of gradients to occupied psi..."

     !project g such that they are orthogonal to all occupied psi. 
     !Gradients do not need orthogonality.
     call orthon_virt_occup(iproc,nproc,orbs,orbsv,comms,commsv,psi,g,msg)
!!$     call  orthoconvirt_p(iproc,nproc,orbs%norbu,orbsv%norb,comms%nvctr_par(iproc,1),psi,g,msg)
!!$     if(orbs%norbd > 0 ) then
!!$        call orthoconvirt_p(iproc,nproc,orbs%norbd,&
!!$          orbsv%norb,comms%nvctr_par(iproc,1),psi(1+comms%nvctr_par(iproc,1)*orbs%norbu),g,msg)
!!$     end if

     call timing(iproc,'Davidson      ','ON')
     if(iproc==0)write(*,'(1x,a)',advance="no")"done."
     if(msg)write(*,'(1x,a)')"squared norm of all gradients after projection"

     call razero(orbsv%norb*orbsv%nkpts,e(1,1,2,ise))
     ispsi=1
     do ikptp=1,orbsv%nkptsp
        ikpt=ikptp+orbsv%iskpts
        nvctrp=commsv%nvctr_par(iproc,ikptp)

        do iorb=1,orbsv%norb
           e(iorb,ikpt,2,ise)= nrm2(nvctrp,g(ispsi+nvctrp*(iorb-1)),1)**2
        end do

        ispsi=ispsi+nvctrp*orbsv%norb*orbsv%nspinor
     end do

     if(nproc > 1)then
        !sum up the contributions of nproc sets with nvctrp wavelet coefficients each
        call MPI_ALLREDUCE(e(1,1,2,2),e(1,1,2,1),orbsv%norb*orbsv%nkpts,&
             mpidtypw,MPI_SUM,MPI_COMM_WORLD,ierr)
     end if

     gnrm=0._dp
     do ikpt=1,orbsv%nkpts
        do iorb=1,nvirt
           tt=real(e(iorb,ikpt,2,1)*orbs%kwgts(ikpt),dp)
           if(msg)write(*,'(1x,i3,1x,1pe21.14)')iorb,tt
           gnrm=gnrm+tt
        end do
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

     !we use for preconditioning the eval from the lowest value of the KS wavefunctions
     !WARNING this pointer association is nor coherent: it will give errors if orbs%norb < orbsv%norb
     orbsv%eval=>orbs%eval
     if (orbs%norb < orbsv%norb) then
        write(*,*)'ERROR: too many virtual orbitals demanded, cannot proceed.'
        write(*,*)'       chenge the eigenvalues array to fix this problem'
        stop
     end if
     call preconditionall(iproc,nproc,orbsv,lr,hx,hy,hz,in%ncong,g,gnrm_fake)

     call timing(iproc,'Precondition  ','OF')
     if (iproc==0)write(*,'(1x,a)')'done.'

     if(iproc==0)write(*,'(1x,a)',advance="no")&
                 "Orthogonality of preconditioned gradients to occupied psi..."

     !transpose  g 
     call transpose_v(iproc,nproc,orbsv,lr%wfd,commsv,g,work=psiw)

     !project g such that they are orthogonal to all occupied psi
     call orthon_virt_occup(iproc,nproc,orbs,orbsv,comms,commsv,psi,g,msg)
!!$     call  orthoconvirt_p(iproc,nproc,orbs%norbu,orbsv%norb,comms%nvctr_par(iproc,1),psi,g,msg)
!!$     if(orbs%norbd > 0) then 
!!$        call orthoconvirt_p(iproc,nproc,orbs%norbd,orbsv%norb,&
!!$          comms%nvctr_par(iproc,1),psi(1+comms%nvctr_par(iproc,1)*orbs%norbu),g,msg)
!!$     end if
     !retranspose the gradient g
     call untranspose_v(iproc,nproc,orbsv,lr%wfd,commsv,g,work=psiw)

     if(iproc==0)write(*,'(1x,a)')"done."

     allocate(hg(orbsv%npsidim+ndebug),stat=i_stat)
     call memocc(i_stat,hg,'hg',subname)

     call HamiltonianApplication(iproc,nproc,at,orbsv,hx,hy,hz,rxyz,cpmult,fpmult,&
          radii_cf,nlpspd,proj,lr,ngatherarr,n1i*n2i*n3p,&
          rhopot(1+i3xcsh*n1i*n2i),g,hg,ekin_sum,epot_sum,eproj_sum,in%nspin,GPU)

     !transpose  g and hg
     call transpose_v(iproc,nproc,orbsv,lr%wfd,commsv,g,work=psiw)
     call transpose_v(iproc,nproc,orbsv,lr%wfd,commsv,hg,work=psiw)

     call timing(iproc,'Davidson      ','ON')
     if(iproc==0)write(*,'(1x,a)',advance="no")"done."


     if(msg)write(*,'(1x,a)')"Norm of all preconditioned gradients"

     call razero(orbsv%norb*orbsv%nkpts,e(1,1,2,ise))

     ispsi=1
     do ikptp=1,orbsv%nkptsp
        ikpt=ikptp+orbsv%iskpts
        nvctrp=commsv%nvctr_par(iproc,ikptp)

        do iorb=1,orbsv%norb
           e(iorb,ikpt,2,ise)=nrm2(nvctrp,g(ispsi+nvctrp*(iorb-1)),1)**2
        end do

        ispsi=ispsi+nvctrp*orbsv%norb*orbsv%nspinor
     end do
     
     if(nproc > 1)then
        !sum up the contributions of nproc sets with nvctrp wavelet coefficients each
        call MPI_ALLREDUCE(e(1,1,2,2),e(1,1,2,1),orbsv%norb*orbsv%nkpts,&
             mpidtypw,MPI_SUM,MPI_COMM_WORLD,ierr)
     end if

     gnrm=0.0_dp
     do ikpt=1,orbsv%nkpts
        do iorb=1,orbsv%norb
           tt=real(e(iorb,ikpt,2,1)*orbs%kwgts(ikpt),dp)
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

     !the dimension should be chosen with the max between k-points
     allocate(hamovr(n2virt,n2virt,orbsv%nkpts,ish+ndebug),stat=i_stat)
     call memocc(i_stat,hamovr,'hamovr',subname)

     ! store upper triangular part of these matrices only
     ! therefore, element (iorb+nvirte,jorb) is transposed to (j,nvirt+iorb)
     ispsi=1
     do ikptp=1,orbsv%nkptsp
        ikpt=ikptp+orbsv%iskpts
        nvctrp=commsv%nvctr_par(iproc,ikptp)

        do iorb=1,orbsv%norb
           do jorb=iorb,orbsv%norb!or 1,nvirte 
              ind=ispsi+nvctrp*(iorb-1)
              jnd=ispsi+nvctrp*(jorb-1)
              hamovr(iorb,jorb,ikpt,ish-1)=               dot(nvctrp,v(ind),1,hv(jnd),1)
              hamovr(jorb,iorb+orbsv%norb,ikpt,ish-1)=        dot(nvctrp,g(ind),1,hv(jnd),1)
              !=hamovr(iorb+orbsv%norb,jorb,ikpt,ish-1)=        dot(nvctrp,g(ind),1,hv(jnd),1)
              hamovr(iorb,jorb+orbsv%norb,ikpt,ish-1)=        dot(nvctrp,v(ind),1,hg(jnd),1)
              hamovr(iorb+orbsv%norb,jorb+orbsv%norb,ikpt,ish-1)= dot(nvctrp,g(ind),1,hg(jnd),1)

              hamovr(iorb,jorb,ikpt,ish)=               dot(nvctrp,v(ind),1, v(jnd),1)
              hamovr(jorb,iorb+orbsv%norb,ikpt,ish)=       dot(nvctrp,g(ind),1, v(jnd),1)
              !=hamovr(iorb+orbsv%norb,jorb,ikpt,ish)=        dot(nvctrp,g(ind),1, v(jnd),1)
              hamovr(iorb,jorb+orbsv%norb,ikpt,ish)=        dot(nvctrp,v(ind),1, g(jnd),1)
              hamovr(iorb+orbsv%norb,jorb+orbsv%norb,ikpt,ish)= dot(nvctrp,g(ind),1, g(jnd),1)
           enddo
        enddo
        ispsi=ispsi+nvctrp*orbsv%norb*orbsv%nspinor
     end do


     !Note: The previous data layout allowed level 3 BLAS
!    call DGEMM('T','N',nvirte,nvirte,nvctrp,1.d0,v(1,1),nvctrp,&
!         hv(1,1),nvctrp,0.d0,hamovr(1,1,ish-1),n2virt)
!    call DSYRK('U','T',n2virt,nvctrp,1.d0,v(1,1),nvctrp,0.d0,hamovr(1,1,ish),n2virt)!upper


     i_all=-product(shape(hg))*kind(hg)
     deallocate(hg,stat=i_stat)
     call memocc(i_stat,i_all,'hg',subname)

     if(nproc > 1)then
        !sum up the contributions of nproc sets with nvctrp wavelet coefficients each
        call MPI_ALLREDUCE(hamovr(1,1,1,3),hamovr(1,1,1,1),2*n2virt**2*orbsv%nkpts,&
             mpidtypw,MPI_SUM,MPI_COMM_WORLD,ierr)
     end if
     
     if(msg)then
        do ikpt=1,orbsv%nkpts
           write(*,*)"subspace matrices, upper triangular (diagonal elements first)"
           write(*,'(1x)')
           write(*,*)"subspace H "
           do iorb=1,n2virt
              write(*,*)hamovr(iorb,iorb:n2virt,ikpt,1)
              write(*,*)
           end do
           write(*,*)"subspace S"
           write(*,*)
           do iorb=1,n2virt
              write(*,*)hamovr(iorb,iorb:n2virt,ikpt,2)
              write(*,*)
           end do
        end do
     end if

     if(iproc==0)write(*,'(1x,a)')"done."
     if(iproc==0)write(*,'(1x,a)',advance='no')"Diagonalization..."

     nwork=max(10,4*n2virt)
     allocate(work(nwork+ndebug),stat=i_stat)
     call memocc(i_stat,work,'work',subname)

     allocate(ew(n2virt+ndebug),stat=i_stat)
     call memocc(i_stat,ew,'ew',subname)


     do ikptp=1,orbsv%nkptsp
        ikpt=ikptp+orbsv%iskpts
        
        call sygv(1,'V','U',n2virt,hamovr(1,1,ikpt,1),n2virt,hamovr(1,1,ikpt,2),n2virt,&
             ew(1),work(1),nwork,i_stat)! Lapack GEVP

        do iorb=1,orbsv%norb
           e(iorb,ikpt,1,1)=ew(iorb)
           e(iorb,ikpt,2,1)=ew(iorb+orbsv%norb)
        end do

     end do

     if (i_stat /= 0) write(*,*) 'Error in DSYGV on process ',iproc,', infocode ', i_stat

     i_all=-product(shape(work))*kind(work)
     deallocate(work,stat=i_stat)
     call memocc(i_stat,i_all,'work',subname)

     i_all=-product(shape(ew))*kind(ew)
     deallocate(ew,stat=i_stat)
     call memocc(i_stat,i_all,'ew',subname)


     if(iproc == 0)write(*,'(1x,a)')'done. The refined eigenvalues are'
     if(msg)then
     write(*,'(1x,a)')'    e(update)           e(not used)'
     do ikptp=1,orbsv%nkptsp
        ikpt=ikptp+orbsv%iskpts
        do iorb=1,orbsv%norb
          write(*,'(1x,i3,2(1pe21.14))')iorb, (e(iorb,ikpt,j,1),j=1,2)
        end do
          write(*,*)
        write(*,*)"and the eigenvectors are"
        write(*,*)
        do iorb=1,n2virt
          write(*,*)hamovr(iorb,:,ikpt,1)!iorb:n2virt,1)
          write(*,*)
        end do
     end do
     else
     do ikptp=1,orbsv%nkptsp
        ikpt=ikptp+orbsv%iskpts
        do iorb=1,nvirt
          if(iproc==0)write(*,'(1x,i3,2(1pe21.14))')iorb,(e(iorb,ikpt,j,1),j=1,2)
        end do
     end do
     end if
     if(iproc==0)write(*,'(1x,a)',advance="no")"Update v with eigenvectors..."

     !Update v, that is the wavefunction, using the eigenvectors stored in hamovr(:,:,1)
     !Lets say we have 4 quarters top/bottom left/right, then
     !v = matmul(v, hamovr(topleft)  ) + matmul(g, hamovr(bottomleft)  )     needed    
     !g=  matmul(v, hamovr(topright) ) + matmul(g, hamovr(bottomright) ) not needed
     !use hv as work arrray

     ispsi=1
     do ikptp=1,orbsv%nkptsp
        ikpt=ikptp+orbsv%iskpts
        nvctrp=commsv%nvctr_par(iproc,ikptp)

        do jorb=1,orbsv%norb! v to update
           call razero(nvctrp,hv(ispsi+nvctrp*(jorb-1)))
           do iorb=1,orbsv%norb ! sum over v and g
              tt=hamovr(iorb,jorb,ikpt,1)
              call axpy(nvctrp,tt,v(ispsi+nvctrp*(iorb-1)),1,&
                   hv(ispsi+nvctrp*(jorb-1)),1)
              tt=hamovr(iorb+orbsv%norb,jorb,ikpt,1)
              call axpy(nvctrp,tt,g(ispsi+nvctrp*(iorb-1)),1,&
                   hv(ispsi+nvctrp*(jorb-1)),1)
           enddo
        enddo

        call dcopy(nvctrp*orbsv%norb,hv(ispsi),1,v(ispsi),1)

        ispsi=ispsi+nvctrp*orbsv%norb*orbsv%nspinor
     end do


     !Note: The previous data layout allowed level 3 BLAS
     !call DGEMM('N','N',nvctrp,nvirte,n2virt,1.d0,v(1,1),nvctrp,hamovr(1,1,1),n2virt,0.d0,hv(1,1),nvctrp)
     !    dimensions    =m      =n   =k          m,k        k,n                   m,n             
     !call DCOPY(nvctrp*nvirte,hv(1,1),1,v(1,1),1)

     i_all=-product(shape(g))*kind(g)
     deallocate(g,stat=i_stat)
     call memocc(i_stat,i_all,'g',subname)

     i_all=-product(shape(hamovr))*kind(hamovr)
     deallocate(hamovr,stat=i_stat)
     call memocc(i_stat,i_all,'hamovr',subname)

     if(iproc==0)write(*,'(1x,a)')"done."
     if(iproc==0)write(*,'(1x,a)',advance="no")"Orthogonality to occupied psi..."
     !project v such that they are orthogonal to all occupied psi
     !Orthogonalize before and afterwards.

     call timing(iproc,'Davidson      ','OF')

     !these routines should work both in parallel or in serial
     call orthogonalize(iproc,nproc,orbsv,commsv,lr%wfd,v)
     !call  orthon_p(iproc,nproc,orbsv%norb,commsv%nvctr_par(iproc,1),lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,v,1)
     call orthon_virt_occup(iproc,nproc,orbs,orbsv,comms,commsv,psi,v,msg)
!!$     call  orthoconvirt_p(iproc,nproc,orbs%norbu,orbsv%norb,comms%nvctr_par(iproc,1),psi,v,msg)
!!$     if(orbs%norbd > 0) then
!!$        call orthoconvirt_p(iproc,nproc,orbs%norbd,orbsv%norb,&
!!$          comms%nvctr_par(iproc,1),psi(1+comms%nvctr_par(iproc,1)*orbs%norbu),v,msg)
!!$     end if
     !and orthonormalize them using "gram schmidt"  (should conserve orthogonality to psi)
     call orthogonalize(iproc,nproc,orbsv,commsv,lr%wfd,v)
     !call  orthon_p(iproc,nproc,orbsv%norb,commsv%nvctr_par(iproc,1),lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,v,1)

     !retranspose v
     call untranspose_v(iproc,nproc,orbsv,lr%wfd,commsv,v,work=psiw)
 
     ! Hamilton application on v
     if(iproc==0)write(*,'(1x,a)',advance="no")"done."
  
     call HamiltonianApplication(iproc,nproc,at,orbsv,hx,hy,hz,rxyz,cpmult,fpmult,&
          radii_cf,nlpspd,proj,lr,ngatherarr,n1i*n2i*n3p,&
          rhopot(1+i3xcsh*n1i*n2i),v,hv,ekin_sum,epot_sum,eproj_sum,in%nspin,GPU)

     !transpose  v and hv
     call transpose_v(iproc,nproc,orbsv,lr%wfd,commsv,v,work=psiw)
     call transpose_v(iproc,nproc,orbsv,lr%wfd,commsv,hv,work=psiw)

     if(iproc==0)write(*,'(1x,a)')"done. "
     call timing(iproc,'Davidson      ','ON')
  end do! davidson iterations

  if(iter>in%itermax)then
     if(iproc==0)write(*,'(1x,a)')'No convergence within the allowed number of minimization steps'
  else
     if(iproc==0)write(*,'(1x,a,i3,a)')'Davidsons method: Convergence after ',iter-1,' iterations.'
  end if
  !finalize: Retranspose, deallocate

  if(iproc==0)then
     do ikpt=1,orbsv%nkpts
        write(*,'(1x,a)')'Complete list of energy eigenvalues'
        do iorb=1,orbs%norb
           write(*,'(1x,a,i4,a,1x,1pe21.14)') 'e_occupied(',iorb,')=',orbs%eval(iorb+(ikpt-1)*orbs%norb)
        end do
        write(*,'(1x,a,1pe21.14)')&
                'HOMO LUMO gap   =',e(1,ikpt,1,1)-orbs%eval(orbs%norb+(ikpt-1)*orbs%norb)
        if(orbs%norbd > 0)write(*,'(1x,a,1pe21.14)')&
                '    and (spin up)',e(1,ikpt,1,1)-orbs%eval(orbs%norbu+(ikpt-1)*orbs%norb) !to be recontrolled
        do iorb=1,orbsv%norb
           write(*,'(1x,a,i4,a,1x,1pe21.14)') 'e_virtual(',iorb,')=',e(iorb,ikpt,1,1)
        end do 
     end do
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

  i_all=-product(shape(orbsv%occup))*kind(orbsv%occup)
  deallocate(orbsv%occup,stat=i_stat)
  call memocc(i_stat,i_all,'orbsv%occup',subname)
  i_all=-product(shape(orbsv%spinsgn))*kind(orbsv%spinsgn)
  deallocate(orbsv%spinsgn,stat=i_stat)
  call memocc(i_stat,i_all,'orbsv%spinsgn',subname)
  i_all=-product(shape(orbsv%kpts))*kind(orbsv%kpts)
  deallocate(orbsv%kpts,stat=i_stat)
  call memocc(i_stat,i_all,'orbsv%kpts',subname)
  i_all=-product(shape(orbsv%kwgts))*kind(orbsv%kwgts)
  deallocate(orbsv%kwgts,stat=i_stat)
  call memocc(i_stat,i_all,'orbsv%kwgts',subname)


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
     write(orbname,'(A,i3.3)')'virtual',iorb+orbsv%isorb
     write(comment,'(1pe10.3)')e(modulo(iorb+orbsv%isorb-1,orbsv%norb)+1,orbsv%iokpt(iorb),1,1)
     !choose the way of plotting the wavefunctions
     if (in%nplot > 0) then
        call plot_wf('POT',orbname,at,lr,hx,hy,hz,rxyz,v(ind:),comment)
     else if (in%nplot < 0) then
        call plot_wf('CUBE',orbname,at,lr,hx,hy,hz,rxyz,v(ind:),comment)
     end if
  end do

  do iorb=orbs%norbp,1,-1 ! sweep over highest occupied orbitals
     if(modulo(orbs%norb-iorb-orbs%isorb-1+nvirt-1,orbs%norb)+1 > abs(in%nplot))exit! we have written nplot pot files
     !address
     ind=1+(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*(iorb-1)
     write(orbname,'(A,i3.3)')'orbital',iorb+orbs%isorb
     write(comment,'(1pe10.3)')orbs%eval(iorb+orbs%isorb)
     !choose the way of plotting the wavefunctions
     if (in%nplot > 0) then
        call plot_wf('POT',orbname,at,lr,hx,hy,hz,rxyz,psi(ind:),comment)
     else if (in%nplot < 0) then
        call plot_wf('CUBE',orbname,at,lr,hx,hy,hz,rxyz,psi(ind:),comment)
     end if

  end do
  ! END OF PLOTTING

  i_all=-product(shape(e))*kind(e)
  deallocate(e,stat=i_stat)
  call memocc(i_stat,i_all,'e',subname)
  i_all=-product(shape(orbsv%iokpt))*kind(orbsv%iokpt)
  deallocate(orbsv%iokpt,stat=i_stat)
  call memocc(i_stat,i_all,'orbsv%iokpt',subname)



END SUBROUTINE davidson
!!***


!!****f* BigDFT/orthoconvirt
!! DESCRIPTION
!!   Makes sure all psivirt/gradients are othogonal to the occupied states psi
!!   This routine is almost the same as orthoconstraint. Only differences:
!!   hpsi(:,norb) -->  psivirt(:,nvirte) , therefore different dimensions.
!!
!! WARNING
!!   Orthogonality to spin polarized channels is achieved in two calls,
!    because up and down orbitals of psi are not orthogonal.
!! SOURCE
!! 
subroutine orthoconvirt(norb,nvirte,nvctrp,psi,hpsi,msg)
  use module_base
  implicit none! real(kind=8) (a-h,o-z)
  integer::norb,nvirte,nvctrp,i_all,i_stat,iorb,jorb,iproc
  logical, parameter :: parallel=.false.
  real(8):: psi(nvctrp,norb),hpsi(nvctrp,nvirte)!,occup(norb)
  real(8), allocatable :: alag(:,:,:)
  real(8)::scprsum,tt
  character(len=*), parameter :: subname='orthoconvirt'
  logical::msg

  iproc=0
  call timing(iproc,'LagrM_comput  ','ON')

  allocate(alag(norb,nvirte,2+ndebug),stat=i_stat)
  call memocc(i_stat,alag,'alag',subname)

  if (nvctrp == 0) then
     call razero(norb*nvirte*2,alag)
  end if
  !     alag(jorb,iorb,2)=+psi(k,jorb)*hpsi(k,iorb)

  call DGEMM('T','N',norb,nvirte,nvctrp,1.d0,psi(1,1),nvctrp,hpsi(1,1),nvctrp,0.d0,alag(1,1,1),norb)

  if(msg)write(*,'(1x,a)')'scalar products are'
  if(msg)write(*,'(1x,a)')'iocc ivirt       value'!                  zero if<1d-12'

  scprsum=0.0_dp
  do iorb=1,norb
   do jorb=1,nvirte
     tt=alag(iorb,jorb,1)
     if(msg)write(*,'(1x,2i3,1pe21.14)')iorb,jorb,tt
     scprsum=scprsum+tt**2
     !if(abs(tt)<1d-12)alag(iorb,jorb,1)=0d0 
     !if(msg)write(*,'(i5,1x,i5,7x,2(1pe21.14,1x))')iorb,jorb,tt,alag(iorb,jorb,1)
   end do
  enddo
  scprsum=dsqrt(scprsum/dble(norb)/dble(nvirte))
  if(msg)write(*,'(1x,a,1pe21.14)')'sqrt sum squares is',scprsum
  if(msg)write(*,'(1x)')
  ! hpsi(k,iorb)=-psi(k,jorb)*alag(jorb,iorb,1)
  !if(maxval(alag(:,:,1))>0d0)&
  call DGEMM('N','N',nvctrp,nvirte,norb,&
             -1.d0,psi(1,1),nvctrp,alag(1,1,1),norb,1.d0,hpsi(1,1),nvctrp)

  i_all=-product(shape(alag))*kind(alag)
  deallocate(alag,stat=i_stat)
  call memocc(i_stat,i_all,'alag',subname)

  call timing(iproc,'LagrM_comput  ','OF')

END SUBROUTINE orthoconvirt
!!***


!!****f* BigDFT/orthoconvirt_p
!! DESCRIPTION
!!   Makes sure all psivirt/gradients are othogonal to the occupied states psi.
!!   This routine is almost the same as orthoconstraint_p. Difference:
!!   hpsi(:,norb) -->  psivirt(:,nvirte) , therefore rectangular alag.
!! 
!! WARNING
!!   Orthogonality to spin polarized channels is achieved in two calls,
!!   because up and down orbitals of psi are not orthogonal.
!! SOURCE
!!
subroutine orthoconvirt_p(iproc,nproc,norb,nvirte,nvctrp,psi,hpsi,msg)
  use module_base
  implicit none
  integer, intent(in) :: norb,nvirte,nvctrp,iproc,nproc
  real(wp), dimension(nvctrp,norb), intent(in) :: psi
  real(wp), dimension(nvctrp,nvirte), intent(out) :: hpsi
  !local variables
  character(len=*), parameter :: subname='orthoconvirt_p'
  logical :: msg
  integer :: i_all,i_stat,ierr,iorb,jorb,istart
  real(wp), dimension(:,:,:), allocatable :: alag
  real(wp) :: scprsum,tt

  istart=1
  if (nproc > 1) istart=2

  call timing(iproc,'LagrM_comput  ','ON')

  allocate(alag(norb,nvirte,istart+ndebug),stat=i_stat)
  call memocc(i_stat,alag,'alag',subname)

  !     alag(jorb,iorb,2)=+psi(k,jorb)*hpsi(k,iorb)
  call DGEMM('T','N',norb,nvirte,nvctrp,1.d0,psi(1,1),nvctrp,hpsi(1,1),nvctrp,&
       0.d0,alag(1,1,istart),norb)

  if (nproc > 1) then
     call timing(iproc,'LagrM_comput  ','OF')
     call timing(iproc,'LagrM_commun  ','ON')
     call MPI_ALLREDUCE(alag(1,1,2),alag(1,1,1),norb*nvirte,&
          mpidtypw,MPI_SUM,MPI_COMM_WORLD,ierr)
     call timing(iproc,'LagrM_commun  ','OF')
     call timing(iproc,'LagrM_comput  ','ON')
  end if

  if(msg)write(*,'(1x,a)')'scalar products are'
  if(msg)write(*,'(1x,a)')'iocc  ivirt       value'!               zero if<1d-12'

  scprsum=0.0_dp
  do iorb=1,norb
   do jorb=1,nvirte
     tt=alag(iorb,jorb,1)
     if(msg) write(*,'(1x,2i3,1pe21.14)')iorb,jorb,tt
     scprsum=scprsum+tt**2
     !if(abs(tt)<1d-12)alag(iorb,jorb,1)=0d0
     !if(msg)write(*,'(2(i3),7x,2(1pe21.14))')iorb,jorb,tt,alag(iorb,jorb,1)
   end do
  enddo
  scprsum=dsqrt(scprsum/dble(norb)/dble(nvirte))
  if(msg)write(*,'(1x,a,1pe21.14)')'sqrt sum squares is',scprsum
  if(msg)write(*,'(1x)')
  !hpsi(k,iorb)=-psi(k,jorb)*alag(jorb,iorb,1)
  !if(maxval(alag(:,:,1))>0d0)
  call DGEMM('N','N',nvctrp,nvirte,norb,&
       -1.d0,psi(1,1),nvctrp,alag(1,1,1),norb,1.d0,hpsi(1,1),nvctrp)

  i_all=-product(shape(alag))*kind(alag)
  deallocate(alag,stat=i_stat)
  call memocc(i_stat,i_all,'alag',subname)

  call timing(iproc,'LagrM_comput  ','OF')

END SUBROUTINE orthoconvirt_p
!!***
