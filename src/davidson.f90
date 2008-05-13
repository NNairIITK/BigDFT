! Davidsons method for iterative diagonalization of virtual Kohn Sham orbitals
! under orthogonality constraints to occupied orbitals psi. The nvirt input
! variable gives the number of unoccupied orbitals for which the exit criterion
! for the gradients norm holds. nvirte = norbe - norb >= nvirt is the number of
! virtual orbitals processed by the method. The dimension of the subspace for
! diagonalization is 2*nvirte = n2virt
!                                                                 Alex Willand
! Algorithm
! _________
! (parallel)


! (tanspose psi, v is already transposed)
! orthogonality of v to psi
! orthogonalize v
! (retranspose v)
! Hamilton(v) --> hv
! transpose v and hv
! Rayleigh quotients  e
! do
!    gradients g= e*v -hv
!    exit condition gnrm
!    orthogonality of g to psi
!    (retranspose g)
!    preconditioning of g
!    (transpose g again)
!    orthogonality of g to psi
!    (retranspose g)
!    Hamilton(g) --> hg
!    (transpose g and hg)
!    subspace matrices H and S
!    DSYGV(H,e,S)  --> H
!    update v with eigenvectors H
!    orthogonality of v to psi
!    orthogonalize v
!    (retranspose v)
!    Hamilton(v) --> hv
!    (transpose v and hv)
! end do
! (retranspose v and psi) 

subroutine davidson(geocode,iproc,nproc,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,n1i,n2i,n3i,at,&
           norb,norbu,norbp,nvirte,nvirtep,nvirt,gnrm_cv,nplot,n1,n2,n3,nvctrp,&
           hx,hy,hz,rxyz,rhopot,occup,i3xcsh,n3p,itermax,wfd,bounds,nlpspd,proj,  & 
           pkernel,ixc,psi,v,eval,ncong,nscatterarr,ngatherarr)
  use module_base
  use module_types
  use module_interfaces, except_this_one => davidson
  implicit none
  include 'mpif.h'
  type(atoms_data), intent(in) :: at
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(nonlocal_psp_descriptors), intent(in) :: nlpspd
  type(convolutions_bounds), intent(in) :: bounds
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: iproc,nproc,norb,norbp,n1,n2,n3,ixc,n1i,n2i,n3i
  integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,i3xcsh,nvctrp,norbu
  integer, intent(in) :: nvirte,nvirtep,nvirt,ncong,n3p,itermax,nplot
  real(gp), dimension(norb), intent(in) :: occup
  real(dp), intent(in) :: gnrm_cv
  real(gp), intent(in) :: hx,hy,hz!convergence criterion for gradients
  integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
  real(dp), dimension(*), intent(in) :: pkernel,rhopot
  !this is a Fortran 95 standard, should be avoided (it is a pity IMHO)
  !real(kind=8), dimension(:,:,:,:), allocatable :: rhopot 
  real(wp), dimension(norb), intent(in) :: eval
  real(wp), dimension(:,:), pointer :: psi,v!=psivirt(nvctrp,nvirtep*nproc) 
                        !v, that is psivirt, is transposed on input and direct on output
  !local variables
  character(len=*), parameter :: subname='davidson'
  character(len=10)::orbname
  logical :: msg !extended output
  integer :: n2virt,n2virtp,ierr,i_stat,i_all,iorb,jorb,iter,nwork,ind,i1,i2!<-last 3 for debug
  integer :: ise,ish
  real(kind=8), external :: ddot,dnrm2
  real(kind=8) :: tt,gnrm,eks,eexcu,vexcu,epot_sum,ekin_sum,ehart,eproj_sum,etol,gnrm_fake
  real(gp), dimension(:), allocatable :: ones
  real(kind=8), dimension(:), allocatable :: work
  real(wp), dimension(:,:), allocatable :: hv,g,hg
  real(wp), dimension(:,:), pointer :: psiw
  real(dp), dimension(:,:,:), allocatable :: e,hamovr
  !last index of e and hamovr are for mpi_allraduce. e (eigenvalues) is also used as 2 work arrays

  !this should be avoided, it creates a copy of all the bounds pointers, occupying more memory
  !better to define a set of pointers kbounds and associate them with a procedure
  !kbounds => bounds%kb
!!$  type(kinetic_bounds):: kbounds  ! this copy is for readability when preoonditioning is called
!!$  kbounds=bounds%kb
  
  msg=.false.! no extended output
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

  n2virt=nvirte*2! the dimension of the subspace


  if (nproc == 1) then
   
  end if
 
  !disassociate work array for transposition in serial
  if (nproc > 1) then

     allocate(psiw(nvctrp,norbp*nproc+ndebug),stat=i_stat)
     call memocc(i_stat,psiw,'psiw',subname)
  else
     psiw => null()
  endif

  !transpose the wavefunction psi 
  call transpose(iproc,nproc,norb,norbp,1,wfd,nvctrp,psi,work=psiw)

  if (nproc > 1) then
     i_all=-product(shape(psiw))*kind(psiw)
     deallocate(psiw,stat=i_stat)
     call memocc(i_stat,i_all,'psiw',subname)
  end if


  if(iproc==0)write(*,'(1x,a)',advance="no")"Orthogonality to occupied psi..."
  !project v such that they are orthogonal to all occupied psi
  !Orthogonalize before and afterwards.

  !this is the same also in serial
  call  orthon_p(iproc,nproc,nvirte,nvctrp,wfd%nvctr_c+7*wfd%nvctr_f,v,1)
  call  orthoconvirt_p(iproc,nproc,norbu,nvirte,nvctrp,psi,v,msg)
  if(norbu<norb)call  orthoconvirt_p(iproc,nproc,norb-norbu,nvirte,nvctrp,&
       psi(1,norbu+1),v,msg)
  !and orthonormalize them using "gram schmidt"  (should conserve orthogonality to psi)
  call  orthon_p(iproc,nproc,nvirte,nvctrp,wfd%nvctr_c+7*wfd%nvctr_f,v,1)


  !retranspose v
  if(nproc > 1)then
     !reallocate the work array with the good sizeXS
     allocate(psiw(nvctrp,nvirtep*nproc+ndebug),stat=i_stat)
     call memocc(i_stat,psiw,'psiw',subname)
  end if

  call untranspose(iproc,nproc,nvirte,nvirtep,1,wfd,nvctrp,v,work=psiw)


  !DEBUG output of initial guess
!  do iorb=iproc*nvirtep+1,min((iproc+1)*nvirtep,nvirt)
!     !adress
!     if (parallel) then
!        ind=1+(wfd%nvctr_c+7*wfd%nvctr_f)*(iorb-iproc*nvirtep-1)
!        i1=mod(ind-1,nvctrp)+1
!        i2=(ind-i1)/nvctrp+1
!     else
!        i1=1
!        i2=iorb-iproc*norbp
!     end if
!
!     if(iorb>5)then
!        exit
!     end if
!     if(iproc==0)write(*,*)'nvrt total is',wfd%nvctr_c+7*wfd%nvctr_f
!     write(orbname,'(A,i3.3)')'vguess-',iorb
!     write(*,*)'1-norm for v guess',iorb,1d0-dnrm2(wfd%nvctr_c+7*wfd%nvctr_f,v(1,iorb-iproc*norbp),1)
!     call plot_wf(orbname,n1,n2,n3,hgrid,wfd%nseg_c,wfd%nvctr_c,wfd%keyg,wfd%keyv,wfd%nseg_f,wfd%nvctr_f,  &
!          rxyz(1,1),rxyz(2,1),rxyz(3,1),v(i1,i2),&
!          bounds%kb%ibyz_c,bounds%gb%ibzxx_c,bounds%gb%ibxxyy_c,bounds%gb%ibyz_ff,bounds%gb%ibzxx_f,&
!          bounds%gb%ibxxyy_f,bounds%ibyyzz_r,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3)
!  end do
 
  ! 1st Hamilton application on psivirt
  if(iproc==0)write(*,'(1x,a)',advance="no")"done. first "

  allocate(hv(nvctrp,nproc*nvirtep+ndebug),stat=i_stat)
  call memocc(i_stat,hv,'hv',subname)
  
  allocate(ones(nvirtep*nproc+ndebug),stat=i_stat)
  call memocc(i_stat,ones,'ones',subname)
  ones(:)=1.0_gp ! a placeholder for the arrays of occupation and spins

!debug
!write(*,*)"shape(rhopot)",shape(rhopot)
!write(*,*)"start adress: ",1,1,1+i3xcsh,1
!write(*,*)"start adress: ",1+(2*n1+31)*(2*n2+31)*nscatterarr(iproc,4)

  call HamiltonianApplication(geocode,iproc,nproc,at,hx,hy,hz,&
       nvirte,nvirtep,ones,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
       wfd,bounds,nlpspd,proj,ngatherarr,n1i*n2i*n3p,&
       rhopot(1+i3xcsh*n1i*n2i),v,hv,ekin_sum,epot_sum,eproj_sum,1,1,ones)

  !if(iproc==0)write(*,'(1x,a)',advance="no")"done. Rayleigh quotients..."

  allocate(e(nvirte,2,ise+ndebug),stat=i_stat)
  call memocc(i_stat,e,'e',subname)

  !transpose  v and hv
  call transpose(iproc,nproc,nvirte,nvirtep,1,wfd,nvctrp,v,work=psiw)
  call transpose(iproc,nproc,nvirte,nvirtep,1,wfd,nvctrp,hv,work=psiw)

  call timing(iproc,'Davidson      ','ON')
  !Timing excludes transposition, hamilton application and preconditioning

  ! Rayleigh quotients.
  do iorb=1,nvirte ! temporary variables 
     e(iorb,1,ise)= ddot(nvctrp,v(1,iorb),1,hv(1,iorb),1)          != <psi|H|psi> 
     e(iorb,2,ise)= dnrm2(nvctrp,v(1,iorb),1)**2                    != <psi|psi> 
  end do

  if(nproc > 1)then
     !sum up the contributions of nproc sets with nvctrp wavelet coefficients each
     call MPI_ALLREDUCE(e(1,1,2),e(1,1,1),2*nvirte,&
          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  end if

  if(iproc==0)write(*,'(1x,a)')"done."
  if(iproc==0)write(*,'(1x,a)')"     sqnorm                Reyleygh quotient"
 
  do iorb=1,nvirte
     !e(:,1,1) = <psi|H|psi> / <psi|psi>
     e(iorb,1,1)=e(iorb,1,1)/e(iorb,2,1)
     if(iproc==0)write(*,'(1x,i3,2(1x,1pe21.14))')iorb, e(iorb,2,1), e(iorb,1,1)
  end do

!if(msg)then
!write(*,*)"******** transposed v,hv 1st elements"
!do iorb=1,10
!  write(*,*)v(iorb,1),hv(iorb,1)
!end do
!write(*,*)"**********"
!end if
  
  !itermax=... use the input variable instead
  do iter=1,itermax
     if(iproc==0)write(*,'(1x,a,i3)')&
     "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~iter",iter
     if(msg) write(*,'(1x,a)')"squared norm of the (nvirt) gradients"

     allocate(g(nvctrp,nproc*nvirtep+ndebug),stat=i_stat)
     call memocc(i_stat,g,'g',subname)

     call DCOPY(nvctrp*nvirtep*nproc,hv(1,1),1,g(1,1),1)! don't overwrite hv
     do iorb=1,nvirte
         call DAXPY(nvctrp,-e(iorb,1,1),v(1,iorb),1,g(1,iorb),1)
         !gradient = hv-e*v
         e(iorb,2,ise)=dnrm2(nvctrp,g(1,iorb),1)**2
         !local contribution to the square norm
     end do

     if(nproc > 1)then
        !sum up the contributions of nproc sets with nvctrp wavelet coefficients each
        call MPI_ALLREDUCE(e(1,2,2),e(1,2,1),nvirte,&
             MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     end if

     gnrm=0_dp
     do iorb=1,nvirt
         tt=e(iorb,2,1)
         if(msg)write(*,'(1x,i3,1x,1pe21.14)')iorb,tt
         gnrm=gnrm+tt
     end do
     gnrm=dsqrt(gnrm/dble(nvirte))

     if(iproc == 0)write(*,'(1x,a,2(1x,1pe12.5))')&
          "|gradient|=gnrm and exit criterion ",gnrm,gnrm_cv
     if(gnrm < gnrm_cv) then
        i_all=-product(shape(g))*kind(g)
        deallocate(g,stat=i_stat)
        call memocc(i_stat,i_all,'g',subname)
        exit ! iteration loop
     end if
     call timing(iproc,'Davidson      ','OF')

     if(iproc==0)write(*,'(1x,a)',advance="no")"Orthogonality of gradients to occupied psi..."

     !project g such that they are orthogonal to all occupied psi. 
     !Gradients do not need orthogonality.
     call  orthoconvirt_p(iproc,nproc,norbu,nvirte,nvctrp,psi,g,msg)
     if(norbu<norb)call  orthoconvirt_p(iproc,nproc,norb-norbu,nvirte,nvctrp,&
          psi(1,norbu+1),g,msg)

     call timing(iproc,'Davidson      ','ON')
     if(iproc==0)write(*,'(1x,a)',advance="no")"done."
     if(msg)write(*,'(1x,a)')"squared norm of all gradients after projection"

     do iorb=1,nvirte
         e(iorb,2,ise)=dnrm2(nvctrp,g(1,iorb),1)**2
     end do

     if(nproc > 1)then
        !sum up the contributions of nproc sets with nvctrp wavelet coefficients each
        call MPI_ALLREDUCE(e(1,2,2),e(1,2,1),nvirte,&
             MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     end if

     gnrm=0_dp

     do iorb=1,nvirte
        tt=e(iorb,2,1)
        if(msg)write(*,'(1x,i3,1x,1pe21.14)')iorb,tt
        gnrm=gnrm+tt
     end do
     gnrm=dsqrt(gnrm/dble(nvirte))

     if(msg)write(*,'(1x,a,2(1x,1pe21.14))')"gnrm of all ",gnrm

     if (iproc==0)write(*,'(1x,a)',advance='no')'Preconditioning...'

     call timing(iproc,'Davidson      ','OF')

     !retranspose the gradient g 
     call untranspose(iproc,nproc,nvirte,nvirtep,1,wfd,nvctrp,g,work=psiw)

     ! Here the gradients norm could be calculated in the direct form instead,
     ! as it is done in hpsiortho before preconditioning. 
     ! However, this should not make a difference and is not really simpler 

     call timing(iproc,'Precondition  ','ON')

     !LG: why we use for preconditioning the eval form the initial input guess?
     call preconditionall(geocode,iproc,nproc,nvirte,nvirtep,&
          n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,hx,hy,hz, &
          ncong,1,wfd,eval,bounds%kb,g,gnrm_fake)

     call timing(iproc,'Precondition  ','OF')
     if (iproc==0)write(*,'(1x,a)')'done.'

     if(iproc==0)write(*,'(1x,a)',advance="no")&
                 "Orthogonality of preconditioned gradients to occupied psi..."

     !transpose  g 
     call transpose(iproc,nproc,nvirte,nvirtep,1,wfd,nvctrp,g,work=psiw)

     !project g such that they are orthogonal to all occupied psi
     call  orthoconvirt_p(iproc,nproc,norbu,nvirte,nvctrp,psi,g,msg)
     if(norbu<norb)call  orthoconvirt_p(iproc,nproc,norb-norbu,nvirte,nvctrp,&
                                           psi(1,norbu+1),g,msg)
     !retranspose the gradient g
     call untranspose(iproc,nproc,nvirte,nvirtep,1,wfd,nvctrp,g,work=psiw)

     if(iproc==0)write(*,'(1x,a)')"done."

     allocate(hg(nvctrp,nvirtep*nproc+ndebug),stat=i_stat)
     call memocc(i_stat,hg,'hg',subname)

     call HamiltonianApplication(geocode,iproc,nproc,at,hx,hy,hz,&
          nvirte,nvirtep,ones,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
          wfd,bounds,nlpspd,proj,ngatherarr,n1i*n2i*n3p,&
          rhopot(1+i3xcsh*n1i*n2i),g,hg,ekin_sum,epot_sum,eproj_sum,1,1,ones)

                              !ixcs
!  and the syntax from init, wfn_diag
!  call HamiltonianApplication(parallel,datacode,iproc,nproc,nat,ntypes,iatype,hgrid,&
!       psppar,npspcode,2*nvirte,2*nvirtep,ones(1),n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
!       wfd,bounds,nlpspd,proj,&
!       ngatherarr,nscatterarr(iproc,2),rhopot(1+(2*n1+31)*(2*n2+31)*nscatterarr(iproc,4)),&
!       g,hg,ekin_sum,epot_sum,eproj_sum,1,ones(1))

     !transpose  g and hg
     call transpose(iproc,nproc,nvirte,nvirtep,1,wfd,nvctrp,g,work=psiw)
     call transpose(iproc,nproc,nvirte,nvirtep,1,wfd,nvctrp,hg,work=psiw)

     call timing(iproc,'Davidson      ','ON')
     if(iproc==0)write(*,'(1x,a)',advance="no")"done."


     if(msg)write(*,'(1x,a)')"Norm of all preconditioned gradients"
     do iorb=1,nvirte
         e(iorb,2,ise)=dnrm2(nvctrp,g(1,iorb),1)**2
     end do
     if(nproc > 1)then
        !sum up the contributions of nproc sets with nvctrp wavelet coefficients each
        call MPI_ALLREDUCE(e(1,2,2),e(1,2,1),nvirte,&
             MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     end if
     gnrm=0.0_dp
     do iorb=1,nvirte
         tt=e(iorb,2,1)
         if(msg)write(*,'(1x,i3,1x,1pe21.14)')iorb,tt
         gnrm=gnrm+tt
     end do
     gnrm=dsqrt(gnrm/dble(nvirte))
     if(msg)write(*,'(1x,a,2(1x,1pe21.14))')"gnrm of all",gnrm

     if(iproc==0)write(*,'(1x,a)',advance="no")"Expanding subspace matrices..."

     !                 <vi | hvj>      <vi | hgj-n>                   <vi | vj>      <vi | gj-n>
     ! hamovr(i,j,1)=                               ;  hamovr(i,j,2)=  
     !                 <gi-n | hvj>  <gi-n | hgj-n>                   <gi-n | vj>  <gi-n | gj-n>

     allocate(hamovr(n2virt,n2virt,ish+ndebug),stat=i_stat)
     call memocc(i_stat,hamovr,'hamovr',subname)

     ! store upper triangular part of these matrices only
     ! therefore, element (iorb+nvirte,jorb) is transposed to (j,nvirt+iorb)
     do iorb=1,nvirte
        do jorb=iorb,nvirte!or 1,nvirte 
           hamovr(iorb,jorb,ish-1)=               ddot(nvctrp,v(1,iorb),1,hv(1,jorb),1)
           hamovr(jorb,iorb+nvirte,ish-1)=        ddot(nvctrp,g(1,iorb),1,hv(1,jorb),1)
           !=hamovr(iorb+nvirte,jorb,ish-1)=        ddot(nvctrp,g(1,iorb),1,hv(1,jorb),1)
           hamovr(iorb,jorb+nvirte,ish-1)=        ddot(nvctrp,v(1,iorb),1,hg(1,jorb),1)
           hamovr(iorb+nvirte,jorb+nvirte,ish-1)= ddot(nvctrp,g(1,iorb),1,hg(1,jorb),1)

           hamovr(iorb,jorb,ish)=               ddot(nvctrp,v(1,iorb),1, v(1,jorb),1)
           hamovr(jorb,iorb+nvirte,ish)=       ddot(nvctrp,g(1,iorb),1, v(1,jorb),1)
           !=hamovr(iorb+nvirte,jorb,ish)=        ddot(nvctrp,g(1,iorb),1, v(1,jorb),1)
           hamovr(iorb,jorb+nvirte,ish)=        ddot(nvctrp,v(1,iorb),1, g(1,jorb),1)
           hamovr(iorb+nvirte,jorb+nvirte,ish)= ddot(nvctrp,g(1,iorb),1, g(1,jorb),1)
        enddo
     enddo

     !Note: The previous data layout allowed level 3 BLAS
!    call DGEMM('T','N',nvirte,nvirte,nvctrp,1.d0,v(1,1),nvctrp,&
!         hv(1,1),nvctrp,0.d0,hamovr(1,1,ish-1),n2virt)
!    call DSYRK('U','T',n2virt,nvctrp,1.d0,v(1,1),nvctrp,0.d0,hamovr(1,1,ish),n2virt)!upper


     i_all=-product(shape(hg))*kind(hg)
     deallocate(hg,stat=i_stat)
     call memocc(i_stat,i_all,'hg',subname)

     if(nproc > 1)then
        !sum up the contributions of nproc sets with nvctrp wavelet coefficients each
        call MPI_ALLREDUCE(hamovr(1,1,3),hamovr(1,1,1),2*n2virt**2,&
             MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     end if
     
     if(msg)then
       write(*,*)"supscpace matrices, upper triangular (diagonal elements first)"
       write(*,'(1x)')
       write(*,*)"supscpace H "
       do iorb=1,n2virt
         write(*,*)hamovr(iorb,iorb:n2virt,1)
         write(*,*)
       end do
       write(*,*)"supscpace S"
       write(*,*)
       do iorb=1,n2virt
         write(*,*)hamovr(iorb,iorb:n2virt,2)
         write(*,*)
       end do
     end if

     if(iproc==0)write(*,'(1x,a)')"done."
     if(iproc==0)write(*,'(1x,a)',advance='no')"Diagonalization..."

     nwork=max(10,4*n2virt)
     allocate(work(nwork+ndebug),stat=i_stat)
     call memocc(i_stat,work,'work',subname)

     call DSYGV(1,'V','U',n2virt,hamovr(1,1,1),n2virt,hamovr(1,1,2),n2virt,&
          e(1,1,1),work,nwork,i_stat)! Lapack GEVP

     if (i_stat.ne.0) write(*,*) 'Error in DSYGV on process ',iproc,', infocode ', i_stat

     i_all=-product(shape(work))*kind(work)
     deallocate(work,stat=i_stat)
     call memocc(i_stat,i_all,'work',subname)

     if(iproc==0)write(*,'(1x,a)')'done. The refined eigenvalues are'
     if(msg)then
     write(*,'(1x,a)')'    e(update)           e(not used)'
        do iorb=1,nvirte
          write(*,'(1x,i3,2(1x,g20.14))')iorb, e(iorb,:,1)
        end do
        write(*,*)
        write(*,*)"and the eigenvectors are"
        write(*,*)
        do iorb=1,n2virt
          write(*,*)hamovr(iorb,:,1)!iorb:n2virt,1)
          write(*,*)
        end do
     else
        do iorb=1,nvirt
          if(iproc==0)write(*,'(1x,i3,2(1x,g20.14))')iorb, e(iorb,:,1)
        end do
     end if
     if(iproc==0)write(*,'(1x,a)',advance="no")"Update v with eigenvectors..."

     !Update v, that is the wavefunction, using the eigenvectors stored in hamovr(:,:,1)
     !Lets say we have 4 quarters top/bottom left/right, then
     !v = matmul(v, hamovr(topleft)  ) + matmul(g, hamovr(bottomleft)  )     needed    
     !g=  matmul(v, hamovr(topright) ) + matmul(g, hamovr(bottomright) ) not needed
     !use hv as work arrray

     do jorb=1,nvirte! v to update
        call razero(nvctrp,hv(1,jorb))
        do iorb=1,nvirte ! sum over v and g
           tt=hamovr(iorb,jorb,1)
           call daxpy(nvctrp,tt,v(1,iorb),1,hv(1,jorb),1)
           tt=hamovr(iorb+nvirte,jorb,1)
           call daxpy(nvctrp,tt,g(1,iorb),1,hv(1,jorb),1)
        enddo
     enddo

     call DCOPY(nvctrp*nvirtep*nproc,hv(1,1),1,v(1,1),1)

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
     call  orthon_p(iproc,nproc,nvirte,nvctrp,wfd%nvctr_c+7*wfd%nvctr_f,v,1)
     call  orthoconvirt_p(iproc,nproc,norbu,nvirte,nvctrp,psi,v,msg)
     if(norbu<norb)call  orthoconvirt_p(iproc,nproc,norb-norbu,nvirte,nvctrp,&
          psi(1,norbu+1),v,msg)
     !and orthonormalize them using "gram schmidt"  (should conserve orthogonality to psi)
     call  orthon_p(iproc,nproc,nvirte,nvctrp,wfd%nvctr_c+7*wfd%nvctr_f,v,1)

     !retranspose v
     call untranspose(iproc,nproc,nvirte,nvirtep,1,wfd,nvctrp,v,work=psiw)
 
     ! Hamilton application on v
     if(iproc==0)write(*,'(1x,a)',advance="no")"done."
  
     call HamiltonianApplication(geocode,iproc,nproc,at,hx,hy,hz,&
          nvirte,nvirtep,ones,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
          wfd,bounds,nlpspd,proj,ngatherarr,n1i*n2i*n3p,&
          rhopot(1+i3xcsh*n1i*n2i),v,hv,ekin_sum,epot_sum,eproj_sum,1,1,ones)

     !transpose  v and hv
     call transpose(iproc,nproc,nvirte,nvirtep,1,wfd,nvctrp,v,work=psiw)
     call transpose(iproc,nproc,nvirte,nvirtep,1,wfd,nvctrp,hv,work=psiw)

     if(iproc==0)write(*,'(1x,a)')"done. "
     call timing(iproc,'Davidson      ','ON')
  end do! davidson iterations

  if(iter>itermax)then
     if(iproc==0)write(*,'(1x,a)')'No convergence within the allowed number of minimization steps'
  else
     if(iproc==0)write(*,'(1x,a,i3,a)')'Davidsons method: Convergence after ',iter-1,' iterations.'
  end if
  !finalize: Retranspose, deallocate

  if(iproc==0)then
        write(*,'(1x,a)')'Complete list of energy eigenvalues'
        do iorb=1,norb
           write(*,'(1x,a,i4,a,1x,1pe21.14)') 'e_occupied(',iorb,')=',eval(iorb)
        end do 
        write(*,'(1x,a,1pe21.14)')&
                'HOMO LUMO gap   =',e(1,1,1)-eval(norb)
        if(norbu<norb)write(*,'(1x,a,1pe21.14)')&
                '    and (spin up)',e(1,1,1)-eval(norbu)! right?
        do iorb=1,nvirte
           write(*,'(1x,a,i4,a,1x,1pe21.14)') 'e_virtual(',iorb,')=',e(iorb,1,1)
        end do 
  end if

  call timing(iproc,'Davidson      ','OF')

  !retranspose v and psi
  call untranspose(iproc,nproc,nvirte,nvirtep,1,wfd,nvctrp,v,work=psiw)

  !resize work array before final transposition
  if(nproc > 1)then
     i_all=-product(shape(psiw))*kind(psiw)
     deallocate(psiw,stat=i_stat)
     call memocc(i_stat,i_all,'psiw',subname)

     allocate(psiw(nvctrp,norbp*nproc+ndebug),stat=i_stat)
     call memocc(i_stat,psiw,'psiw',subname)
  end if

  call untranspose(iproc,nproc,norb,norbp,1,wfd,nvctrp,psi,work=psiw)

  if(nproc > 1) then
     i_all=-product(shape(psiw))*kind(psiw)
     deallocate(psiw,stat=i_stat)
     call memocc(i_stat,i_all,'psiw',subname)
  end if

  i_all=-product(shape(hv))*kind(hv)
  deallocate(hv,stat=i_stat)
  call memocc(i_stat,i_all,'hv',subname)

  i_all=-product(shape(e))*kind(e)
  deallocate(e,stat=i_stat)
  call memocc(i_stat,i_all,'e',subname)

  i_all=-product(shape(ones))*kind(ones)
  deallocate(ones,stat=i_stat)
  call memocc(i_stat,i_all,'ones',subname)

  !THIS PART WAS INSIDE CLUSTER ROUTINE

  !plot the converged wavefunctions in the different orbitals.
  !nplot is the requested total of orbitals to plot, where
  !states near the HOMO/LUMO gap are given higher priority.
  !Occupied orbitals are only plotted when nplot>nvirt,
  !otherwise a comment is given in the out file.
  
  if(nplot>norb+nvirt)then
     if(iproc==0)write(*,'(1x,A,i3)')&
          "WARNING: More plots requested than orbitals calculated." 
  end if

  do iorb=iproc*nvirtep+1,min((iproc+1)*nvirtep,nvirt)!requested: nvirt of nvirte orbitals
     if(iorb>nplot)then
        if(iproc==0.and.nplot>0)write(*,'(A)')&
             'WARNING: No plots of occupied orbitals requested.'
        exit 
     end if
     !calculate the address to start from, since psivirt and
     !psi are allocated in the transposed way. Physical shape
     !is (nvtrp,norbp*nproc), logical shape is (nvctr_tot,norbp).
     if (nproc > 1) then
        ind=1+(wfd%nvctr_c+7*wfd%nvctr_f)*(iorb-iproc*nvirtep-1)
        i1=mod(ind-1,nvctrp)+1
        i2=(ind-i1)/nvctrp+1
     else
        i1=1
        i2=iorb-iproc*norbp
     end if

     write(orbname,'(A,i3.3)')'virtual',iorb
     call plot_wf(orbname,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,hx,rxyz(1,1),rxyz(2,1),rxyz(3,1),wfd,&
          bounds,v(i1:,i2:))
  end do

  do iorb=min((iproc+1)*norbp,norb),iproc*norbp+1,-1 ! sweep over highest occupied orbitals
     if(norb-iorb+1+nvirt>nplot)exit! we have written nplot pot files
     !adress
     if (nproc > 1) then
        ind=1+(wfd%nvctr_c+7*wfd%nvctr_f)*(iorb-iproc*norbp-1)
        i1=mod(ind-1,nvctrp)+1
        i2=(ind-i1)/nvctrp+1
     else
        i1=1
        i2=iorb-iproc*norbp
     end if

     write(orbname,'(A,i3.3)')'orbital',iorb
     call plot_wf(orbname,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,hx,rxyz(1,1),rxyz(2,1),rxyz(3,1),wfd,&
          bounds,psi(i1:,i2:))
  end do
  ! END OF PLOTTING


end subroutine davidson




subroutine orthoconvirt(norb,nvirte,nvctrp,psi,hpsi,msg)
  !Makes sure all psivirt/gradients are othogonal to the occupied states psi
  !This routine is almost the same as orthoconstraint. Only differences:
  !hpsi(:,norb) -->  psivirt(:,nvirte) , therefore different dimensions.

  !Note: Orthogonality to spin polarized channels is achieved in two calls,
  !because up and down orbitals of psi are not orthogonal.
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


!makes sure all psivirt/gradients are othogonal to the occupied states psi.
!This routine is almost the same as orthoconstraint_p. Difference:
!hpsi(:,norb) -->  psivirt(:,nvirte) , therefore rectangular alag.

!Note: Orthogonality to spin polarized channels is achieved in two calls,
!because up and down orbitals of psi are not orthogonal.
subroutine orthoconvirt_p(iproc,nproc,norb,nvirte,nvctrp,psi,hpsi,msg)
  use module_base
  implicit none
  integer, intent(in) :: norb,nvirte,nvctrp,iproc,nproc
  real(wp), dimension(nvctrp,norb), intent(in) :: psi
  real(wp), dimension(nvctrp,nvirte), intent(out) :: hpsi
  !local variables
  include 'mpif.h'
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
          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
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
