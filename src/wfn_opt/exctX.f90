!!$!!****f* BigDFT/exact_exchange_potential
!!$!! FUNCTION
!!$!!    Calculate the exact exchange potential
!!$!!
!!$!! COPYRIGHT
!!$!!    Copyright (C) 2009-2010 BigDFT group 
!!$!!    This file is distributed under the terms of the
!!$!!    GNU General Public License, see ~/COPYING file
!!$!!    or http://www.gnu.org/copyleft/gpl.txt .
!!$!!    For the list of contributors, see ~/AUTHORS 
!!$!!
!!$!! AUTHOR
!!$!!    Luigi Genovese
!!$!!
!!$!! SOURCE
!!$!! 
!!$!second version, with MPI non-blocking communications
!!$subroutine exact_exchange_potential_round(iproc,nproc,geocode,nspin,lr,orbs,&
!!$     hxh,hyh,hzh,pkernelseq,psi,psir,eexctX)
!!$  use module_base
!!$  use module_types
!!$  use Poisson_Solver
!!$  use libxc_functionals
!!$  implicit none
!!$  character(len=1), intent(in) :: geocode
!!$  integer, intent(in) :: iproc,nproc,nspin
!!$  real(gp), intent(in) :: hxh,hyh,hzh
!!$  type(locreg_descriptors), intent(in) :: lr
!!$  type(orbitals_data), intent(in) :: orbs
!!$  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(in) :: psi
!!$  real(dp), dimension(*), intent(in) :: pkernelseq !this kernel should be built in sequential
!!$  real(gp), intent(out) :: eexctX
!!$  !real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,orbs%norbp), intent(out) :: XCi
!!$  !local variables
!!$  character(len=*), parameter :: subname='exact_exchange_potential_round'
!!$  integer :: i_all,i_stat,ierr,ispinor,ispsiw,ispin,norb,jpsend,jprecv,ishift,mpirequest,mpirequestr
!!$  integer :: i1,i2,i3,iorb,iorbs,jorb,jorbs,ispsir,ind3,ind2,ind1i,ind1j,jproc,igran,ngran,norbpm
!!$  real(gp) :: ehart,zero,hfac,exctXfac,sign,sfac,hfaci,hfacj,kerneloff,hfactor
!!$  type(workarr_sumrho) :: w
!!$  integer, dimension(MPI_STATUS_SIZE) :: istatus
!!$  integer, dimension(:), allocatable :: mpirequests,jsorb
!!$  real(wp), dimension(:,:,:,:), allocatable :: rp_ij
!!$  real(wp), dimension(:,:,:,:,:), allocatable :: psir
!!$  real(wp), dimension(:,:,:,:,:,:), allocatable :: psisr
!!$
!!$  !call timing(iproc,'Exchangecorr  ','ON')
!!$
!!$  exctXfac = libxc_functionals_exctXfac()
!!$
!!$  eexctX=0.0_gp
!!$
!!$  if (nspin==2) then
!!$     sfac=1.0_gp
!!$  else 
!!$     sfac=0.5_gp
!!$  end if
!!$  hfactor=1/(hxh*hyh*hzh)
!!$  !orbital quantities
!!$  allocate(jsorb(0:nproc-1+ndebug),stat=i_stat)
!!$  call memocc(i_stat,jsorb,'jsorb',subname)
!!$  !maximum value of orbitals
!!$  norbpm=maxval(orbs%norb_par(:))
!!$  !starting orbital for any processor
!!$  jsorb(0)=0
!!$  do jproc=1,nproc-1
!!$     jsorb(jproc)=jsorb(jproc-1)+orbs%norb_par(jproc-1)
!!$  end do
!!$
!!$  call initialize_work_arrays_sumrho(lr,w)
!!$  
!!$  !allocate arrays for the procedure
!!$  !partial densities with always granularity one
!!$  allocate(rp_ij(lr%d%n1i,lr%d%n2i,lr%d%n3i,orbs%norbp*(orbs%norbp+1)/2+ndebug),stat=i_stat)
!!$  call memocc(i_stat,rp_ij,'rp_ij',subname)
!!$  !final result of the exctX potential
!!$  allocate(psir(lr%d%n1i,lr%d%n2i,lr%d%n3i,orbs%nspinor,orbs%norbp+ndebug),stat=i_stat)
!!$  call memocc(i_stat,psir,'psir',subname)
!!$  !array for round-robin communication
!!$  allocate(psisr(lr%d%n1i,lr%d%n2i,lr%d%n3i,orbs%nspinor,norbpm,2+ndebug),stat=i_stat)
!!$  call memocc(i_stat,psisr,'psisr',subname)
!!$  !array for mpi_requests
!!$  allocate(mpirequests(0:nproc-1+ndebug),stat=i_stat)
!!$  call memocc(i_stat,mpirequests,'mpirequests',subname)
!!$
!!$  if (geocode == 'F') then
!!$     call razero(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspinor*orbs%norbp,psir)
!!$  end if
!!$
!!$  !uncompress the wavefunction in the real grid 
!!$  !and copy them on the psir array
!!$  ispsiw=1
!!$  do iorb=1,orbs%norbp
!!$     do ispinor=1,orbs%nspinor
!!$        call daub_to_isf(lr,w,psi(1,ispinor,iorb),psir(1,1,1,ispinor,iorb))
!!$     end do
!!$  end do
!!$  call deallocate_work_arrays_sumrho(w)
!!$
!!$  !copy the psir array in the sending space
!!$  call dcopy(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspinor*orbs%norbp,psir,1,psisr,1)
!!$
!!$  !initialise the sending arrays to iproc
!!$  psisr=real(iproc,wp)
!!$
!!$  !for any of the processes send the wavefunctions to the next one (torus topology)
!!$  do ishift=0,nproc/2
!!$
!!$     jprecv=modulo(iproc+ishift+1,nproc)
!!$     jpsend=modulo(iproc-ishift,nproc)
!!$     
!!$     !the tag is associated to the receiving process
!!$     !first of all, send the psi for the next processor
!!$     call MPI_ISEND(psisr(1,1,1,1,1,1),lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspinor*orbs%norb_par(iproc),&
!!$          mpidtypw,jprecv,jprecv,MPI_COMM_WORLD,mpirequest,ierr)
!!$
!!$     !build the partial densities between the orbitals of the first processor and the one which has been 
!!$     !communicated
!!$     if (ishift == 0 ) then
!!$        jorbs=iorb
!!$     else
!!$        jorbs=1
!!$     end if
!!$     ind=0
!!$     do iorb=1,orbs%norbp
!!$        do jorb=jorbs,orbs%norb_par(jpsend)
!!$           !if (orbs%spinsgn(iorb+orbs%isorb) == orbs%spinsgn(jorb+orbs%isorb)) then
!!$              ind=ind+1
!!$              !calculate partial density (real functions), no spin-polarisation
!!$              do i3=1,lr%d%n3i
!!$                 do i2=1,lr%d%n2i
!!$                    do i1=1,lr%d%n1i
!!$                       rp_ij(i1,i2,i3,ind)=hfactor*psir(i1,i2,i3,1,iorb)*psisr(i1,i2,i3,1,jorb,2)
!!$                    end do
!!$                 end do
!!$              end do
!!$           !end if
!!$        end do
!!$     end do
!!$
!!$     !now the wavefunctions can be received, psisr not anymore necessary
!!$     !receive the wavefunctions to avoid deadlocks
!!$     call MPI_IRECV(psisr(1,1,1,1,1,2),lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspinor*orbs%norb_par(jpsend),&
!!$          mpidtypw,jpsend,iproc,MPI_COMM_WORLD,mpirequestr,ierr)
!!$
!!$     
!!$     !the poisson solver can be executed for all the partial densities
!!$     ind=0
!!$     do iorb=1,orbs%norbp
!!$        do jorb=jorbs,orbs%norbp
!!$           !if (orbs%spinsgn(iorb+orbs%isorb) == orbs%spinsgn(jorb+orbs%isorb)) then
!!$              ind=ind+1
!!$              !if (iproc == 0 .and. verbose > 1) then
!!$              !   write(*,*)'Exact exchange calculation: spin, orbitals:',ispin,iorb,jorb
!!$              !end if
!!$              !poisson solver in sequential mode, thus no communications
!!$              !call H_potential(geocode,'D',0,1,lr%d%n1i,lr%d%n2i,lr%d%n3i,hxh,hyh,hzh,&
!!$              !     rp_ij(1,1,1,ind),pkernelseq,rp_ij,ehart,0.0_dp,.false.,&
!!$              !     quiet='YES')
!!$              !hartree energy
!!$              !this factor is only valid with one k-point
!!$              !can be easily generalised to the k-point case
!!$              !here the occupation numbers have to be properly calculated
!!$              hfac=sfac*orbs%occup(iorb)*orbs%occup(jorb)
!!$              eexctX=eexctX+2.0_gp*hfac*real(ehart,gp)
!!$           end if
!!$        end do
!!$     end do
!!$
!!$     !now the partial potentials are known, they have to be used for building the potential 
!!$     !either in the present or in the sending processor
!!$     !fill the exctX potential with the results coming from partial potential
!!$     ind=0
!!$     do iorb=1,orbs%norbp
!!$        do jorb=iorb+1,orbs%norbp
!!$           if (orbs%spinsgn(iorb+orbs%isorb) == orbs%spinsgn(jorb+orbs%isorb)) then
!!$              ind=ind+1
!!$              !this factor is only valid with one k-point
!!$              !we have to correct with the kwgts if we want more than one k-point
!!$              hfaci=-sfac*orbs%occup(jorb)
!!$              hfacj=-sfac*orbs%occup(iorb)
!!$
!!$              !accumulate the results for each of the wavefunctions concerned
!!$              do i3=1,lr%d%n3i
!!$                 do i2=1,lr%d%n2i
!!$                    do i1=1,lr%d%n1i
!!$                       psir(i1,i2,i3,1,iorb)=psir(i1,i2,i3,1,iorb)+&
!!$                            hfaci*rp_ij(i1,i2,i3,ind)*psiw(i1,i2,i3,1,jorb)
!!$                       psir(i1,i2,i3,1,jorb)=psir(i1,i2,i3,1,jorb)+&
!!$                            hfacj*rp_ij(i1,i2,i3,ind)*psiw(i1,i2,i3,1,iorb)
!!$                    end do
!!$                 end do
!!$              end do
!!$           end if
!!$        end do
!!$     end do
!!$

!!$
!!$     !the exact exchange energy is half the Hartree energy (which already has another half)
!!$     eexctX=-exctXfac*eexctX
!!$
!!$     if (iproc == 0) write(*,'(a,1x,1pe18.11)')'Exact Exchange Energy:',eexctX
!!$     
!!$     
!!$     !wait until the transfer has completed (data have been received)
!!$     call MPI_WAIT(mpirequestr,istatus)     
!!$     !verify the difference with the wavefunctions
!!$     hfac=0.0_wp
!!$     loop_check : do iorb=1,orbs%norb_par(jpsend)
!!$        do i3=1,lr%d%n3i
!!$           do i2=1,lr%d%n2i
!!$              do i1=1,lr%d%n1i
!!$                 hfac=max(hfac,psisr(i1,i2,i3,1,iorb,2)-real(jpsend,wp))
!!$                 if (hfac/=0.0_wp) exit loop_check
!!$              end do
!!$           end do
!!$        end do
!!$     end do loop_check
!!$     
!!$     if (hfac /= 0.0_wp) then
!!$        print *,'error from ',jpsend,' to ',iproc,' in: ',i1,i2,i3,iorb,' :',hfac,psisr(i1,i2,i3,1,iorb,2)
!!$     else
!!$        print *,'maxdiff,iproc',iproc,ishift,jpsend,jprecv,hfac
!!$     end if
!!$  end do
!!$
!!$  i_all=-product(shape(rp_ij))*kind(rp_ij)
!!$  deallocate(rp_ij,stat=i_stat)
!!$  call memocc(i_stat,i_all,'rp_ij',subname)
!!$  
!!$  i_all=-product(shape(psir))*kind(psir)
!!$  deallocate(psir,stat=i_stat)
!!$  call memocc(i_stat,i_all,'psir',subname)
!!$  i_all=-product(shape(psisr))*kind(psisr)
!!$  deallocate(psisr,stat=i_stat)
!!$  call memocc(i_stat,i_all,'psisr',subname)
!!$  i_all=-product(shape(mpirequests))*kind(mpirequests)
!!$  deallocate(mpirequests,stat=i_stat)
!!$  call memocc(i_stat,i_all,'mpirequests',subname)
!!$  i_all=-product(shape(jsorb))*kind(jsorb)
!!$  deallocate(jsorb,stat=i_stat)
!!$  call memocc(i_stat,i_all,'jsorb',subname)
!!$
!!$
!!$stop
!!$  
!!$  !build the partial densities for the poisson solver, calculate the partial potential
!!$  !and accumulate the result
!!$
!!$
!!$end subroutine exact_exchange_potential_round

subroutine exact_exchange_potential(iproc,nproc,geocode,nspin,lr,orbs,n3parr,n3p,&
     hxh,hyh,hzh,pkernel,psi,psir,eexctX)
  use module_base
  use module_types
  use Poisson_Solver
  use libxc_functionals
  implicit none
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: iproc,nproc,n3p,nspin
  real(gp), intent(in) :: hxh,hyh,hzh
  type(locreg_descriptors), intent(in) :: lr
  type(orbitals_data), intent(in) :: orbs
  integer, dimension(0:nproc-1), intent(in) :: n3parr
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(in) :: psi
  real(dp), dimension(*), intent(in) :: pkernel
  real(gp), intent(out) :: eexctX
  real(wp), dimension(max(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%norbp,n3parr(0)*orbs%norb)), intent(out) :: psir
  !local variables
  character(len=*), parameter :: subname='exact_exchange_potential'
  integer :: i_all,i_stat,ierr,ispinor,ispsiw,ispin,norb
  integer :: i1,i2,i3p,iorb,iorbs,jorb,jorbs,ispsir,ind3,ind2,ind1i,ind1j,jproc,igran,ngran
  real(gp) :: ehart,hfac,exctXfac,sign,sfac,hfaci,hfacj
  type(workarr_sumrho) :: w
  integer, dimension(:,:), allocatable :: ncommarr
  real(wp), dimension(:), allocatable :: psiw
  real(wp), dimension(:,:,:,:), allocatable :: rp_ij

  !call timing(iproc,'Exchangecorr  ','ON')

  exctXfac = libxc_functionals_exctXfac()

  eexctX=0.0_gp

  call initialize_work_arrays_sumrho(lr,w)
  
  !save the value of the previous offset of the kernel
  !kerneloff=pkernel(1)
  !put to szero the offset to subtract the energy of the momentum
  !pkernel(1)=0.0_dp

  !the granularity of the calculation is set by ngran
  !for the moment it is irrelevant but if the poisson solver is modified
  !we may increase this value
  ngran=1

  !partial densities with a given granularity
  allocate(rp_ij(lr%d%n1i,lr%d%n2i,n3p,ngran+ndebug),stat=i_stat)
  call memocc(i_stat,rp_ij,'rp_ij',subname)
  allocate(psiw(max(max(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%norbp,n3parr(0)*orbs%norb),1)+ndebug),stat=i_stat)
  call memocc(i_stat,psiw,'psiw',subname)

  if (geocode == 'F') then
     call razero(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%norbp,psiw)
  end if

  !uncompress the wavefunction in the real grid
  !and switch the values of the function
  ispinor=1
  ispsiw=1
  do iorb=1,orbs%norbp
     call daub_to_isf(lr,w,psi(1,ispinor,iorb),psiw(ispsiw))
     ispsir=1+(iorb-1)*n3parr(0)
     do jproc=0,nproc-1
        !write(*,'(a,1x,8(i10))'),'iproc,jproc',iproc,jproc,iorb,orbs%norbp,ispsir,ispsiw,&
        !     lr%d%n1i*lr%d%n2i*max(lr%d%n3i*orbs%norbp,n3p*orbs%norb),n3parr(jproc)
        call dcopy(n3parr(jproc),psiw(ispsiw),1,psir(ispsir),1)
        ispsiw=ispsiw+n3parr(jproc)
        if (jproc /= nproc-1) then
           do jorb=iorb,orbs%norbp
              ispsir=ispsir+n3parr(jproc)
           end do
           do jorb=1,iorb-1
              ispsir=ispsir+n3parr(jproc+1)
           end do
        end if
     end do
  end do
  call deallocate_work_arrays_sumrho(w)

  !communicate them between processors
  if (nproc > 1) then
     !arrays for the communication between processors
     !valid only for one k-point for the moment
     !and only real functions (nspinor=1)
     !this distribution is in principle valid also for k-points

     allocate(ncommarr(0:nproc-1,4+ndebug),stat=i_stat)
     call memocc(i_stat,ncommarr,'ncommarr',subname)

     !count array for orbitals => components
     do jproc=0,nproc-1
        ncommarr(jproc,1)=n3parr(jproc)*orbs%norb_par(iproc)
     end do
     !displacement array for orbitals => components
     ncommarr(0,2)=0
     do jproc=1,nproc-1
        ncommarr(jproc,2)=ncommarr(jproc-1,2)+ncommarr(jproc-1,1)
     end do
     !count array for components => orbitals
     do jproc=0,nproc-1
        ncommarr(jproc,3)=n3parr(iproc)*orbs%norb_par(jproc)
     end do
     !displacement array for components => orbitals
     ncommarr(0,4)=0
     do jproc=1,nproc-1
        ncommarr(jproc,4)=ncommarr(jproc-1,4)+ncommarr(jproc-1,3)
     end do

     call MPI_ALLTOALLV(psir,ncommarr(0,1),ncommarr(0,2),mpidtypw, &
          psiw,ncommarr(0,3),ncommarr(0,4),mpidtypw,MPI_COMM_WORLD,ierr)

  else
     call dcopy(lr%d%n1i*lr%d%n2i*n3p*orbs%norb,psir,1,psiw,1)
  end if

  !this is the array of the actions of the X potential on psi
  call razero(lr%d%n1i*lr%d%n2i*n3p*orbs%norb,psir)

  !build the partial densities for the poisson solver, calculate the partial potential
  !and accumulate the result
  !do it for different spins
  !for non spin-polarised systems there is a factor of two
  !non-collinear spin not yet implemented
  if (nspin==2) then
     sfac=1.0_gp
  else 
     sfac=0.5_gp
  end if

  !number of orbitals, all quantum numbers
  norb=orbs%norb!*orbs%nkpts, for the future
  iorb=1
  jorb=1
  
  do ispin=1,nspin
     if (ispin==1) then
        iorb=1
        jorb=1
        norb=orbs%norbu
        sign=1.0_gp
     else
        iorb=orbs%norbu+1
        jorb=orbs%norbu+1
        norb=orbs%norb
        sign=-1.0_gp
     end if
     orbital_loop: do
        iorbs=iorb
        jorbs=jorb
        hfac=1/(hxh*hyh*hzh)
        do igran=1,ngran
           if (iorb > norb) exit orbital_loop
           if (orbs%spinsgn(iorb) == sign .and. orbs%spinsgn(jorb) == sign) then
              !calculate partial density (real functions), no spin-polarisation
              do i3p=1,n3p
                 ind3=(i3p-1)*lr%d%n1i*lr%d%n2i
                 do i2=1,lr%d%n2i
                    ind2=(i2-1)*lr%d%n1i+ind3
                    do i1=1,lr%d%n1i
                       ind1i=i1+ind2+(iorb-1)*lr%d%n1i*lr%d%n2i*n3p
                       ind1j=i1+ind2+(jorb-1)*lr%d%n1i*lr%d%n2i*n3p
                       rp_ij(i1,i2,i3p,igran)=hfac*psiw(ind1i)*psiw(ind1j)
                    end do
                 end do
              end do
           end if
           jorb=jorb+1
           if (jorb > norb) then
              iorb=iorb+1
              jorb=iorb
           end if
        end do
        jorb=jorbs
        iorb=iorbs
        do igran=1,ngran
           if (iorb > norb) exit orbital_loop
           if (orbs%spinsgn(iorb) == sign .and. orbs%spinsgn(jorb) == sign) then
              !this factor is only valid with one k-point
              !can be easily generalised to the k-point case
              hfac=sfac*orbs%occup(iorb)*orbs%occup(jorb)

              !print *,'test',iproc,iorb,jorb,sum(rp_ij(:,:,:,igran))
              !partial exchange term for each partial density
              if (iproc == 0 .and. verbose > 1) then
                 write(*,*)'Exact exchange calculation: spin, orbitals:',ispin,iorb,jorb
              end if
              call H_potential(geocode,'D',iproc,nproc,&
                   lr%d%n1i,lr%d%n2i,lr%d%n3i,hxh,hyh,hzh,&
                   rp_ij(1,1,1,igran),pkernel,rp_ij,ehart,0.0_dp,.false.,&
                   quiet='YES')

!!$              call PSolver(geocode,'D',iproc,nproc,lr%d%n1i,lr%d%n2i,lr%d%n3i,&
!!$                   0,hxh,hyh,hzh,rp_ij(1,1,1,igran),pkernel,rp_ij,ehart,zero,zero,&
!!$                   0.d0,.false.,1,quiet='YES')
              if (iorb==jorb) then
                 eexctX=eexctX+hfac*real(ehart,gp)
              else
                 eexctX=eexctX+2.0_gp*hfac*real(ehart,gp)
              end if
              !print *,'PSOLVER,ehart,iproc',iproc,ehart,hfac
           end if
           jorb=jorb+1
           if (jorb > norb) then
              iorb=iorb+1
              jorb=iorb
           end if
        end do
        jorb=jorbs
        iorb=iorbs
        do igran=1,ngran
           if (iorb > norb) exit orbital_loop
           if (orbs%spinsgn(iorb) == sign .and. orbs%spinsgn(jorb) == sign) then
              !this factor is only valid with one k-point
              !we have to correct with the kwgts if we want more than one k-point
              hfaci=-sfac*orbs%occup(jorb)
              hfacj=-sfac*orbs%occup(iorb)

              if (iorb /= jorb) then
                 !accumulate the results for each of the wavefunctions concerned
                 do i3p=1,n3p
                    ind3=(i3p-1)*lr%d%n1i*lr%d%n2i
                    do i2=1,lr%d%n2i
                       ind2=(i2-1)*lr%d%n1i+ind3
                       do i1=1,lr%d%n1i
                          ind1i=i1+ind2+(iorb-1)*lr%d%n1i*lr%d%n2i*n3p
                          ind1j=i1+ind2+(jorb-1)*lr%d%n1i*lr%d%n2i*n3p
                          psir(ind1i)=psir(ind1i)+hfaci*rp_ij(i1,i2,i3p,igran)*psiw(ind1j)
                          psir(ind1j)=psir(ind1j)+hfacj*rp_ij(i1,i2,i3p,igran)*psiw(ind1i)
                       end do
                    end do
                 end do
              else
                 !accumulate the results for each of the wavefunctions concerned
                 do i3p=1,n3p
                    ind3=(i3p-1)*lr%d%n1i*lr%d%n2i
                    do i2=1,lr%d%n2i
                       ind2=(i2-1)*lr%d%n1i+ind3
                       do i1=1,lr%d%n1i
                          ind1i=i1+ind2+(iorb-1)*lr%d%n1i*lr%d%n2i*n3p
                          ind1j=i1+ind2+(jorb-1)*lr%d%n1i*lr%d%n2i*n3p
                          psir(ind1i)=psir(ind1i)+hfacj*rp_ij(i1,i2,i3p,igran)*psiw(ind1j)
                       end do
                    end do
                 end do
              end if
           end if
           jorb=jorb+1
           if (jorb > norb) then
              iorb=iorb+1
              jorb=iorb
           end if
        end do
     end do orbital_loop
  end do

  !the exact exchange energy is half the Hartree energy (which already has another half)
  eexctX=-exctXfac*eexctX

  if (iproc == 0) write(*,'(a,1x,1pe18.11)')'Exact Exchange Energy:',eexctX

  !assign the potential for each function
  if (nproc > 1) then
     !call dcopy(lr%d%n1i*lr%d%n2i*n3p*orbs%norb,psir,1,psirt,1)
     !recommunicate the values in the psir array
     call MPI_ALLTOALLV(psir,ncommarr(0,3),ncommarr(0,4),mpidtypw, &
          psiw,ncommarr(0,1),ncommarr(0,2),mpidtypw,MPI_COMM_WORLD,ierr)
     !redress the potential
     ispsiw=1
     do iorb=1,orbs%norbp
        ispsir=1+(iorb-1)*n3parr(0)
        do jproc=0,nproc-1
           call dcopy(n3parr(jproc),psiw(ispsir),1,psir(ispsiw),1)
           ispsiw=ispsiw+n3parr(jproc)
           if (jproc /= nproc-1) then
              do jorb=iorb,orbs%norbp
                 ispsir=ispsir+n3parr(jproc)
              end do
              do jorb=1,iorb-1
                 ispsir=ispsir+n3parr(jproc+1)
              end do
           end if
        end do
     end do
  end if

  i_all=-product(shape(rp_ij))*kind(rp_ij)
  deallocate(rp_ij,stat=i_stat)
  call memocc(i_stat,i_all,'rp_ij',subname)
  
  i_all=-product(shape(psiw))*kind(psiw)
  deallocate(psiw,stat=i_stat)
  call memocc(i_stat,i_all,'psiw',subname)


  if (nproc > 1) then
     i_all=-product(shape(ncommarr))*kind(ncommarr)
     deallocate(ncommarr,stat=i_stat)
     call memocc(i_stat,i_all,'ncommarr',subname)
  end if

  !call timing(iproc,'Exchangecorr  ','OF')

END SUBROUTINE exact_exchange_potential
!!***


!!****f* BigDFT/prepare_psirocc
!! SOURCE
!! 
subroutine prepare_psirocc(iproc,nproc,lr,orbsocc,n3p,n3parr,psiocc,psirocc)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,nproc,n3p
  type(locreg_descriptors), intent(in) :: lr
  type(orbitals_data), intent(in) :: orbsocc
  integer, dimension(0:nproc-1), intent(in) :: n3parr
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbsocc%nspinor,orbsocc%norbp), intent(in) :: psiocc
  real(wp), dimension(max(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbsocc%norbp,n3parr(0)*orbsocc%norb)), intent(out) :: psirocc
  !local variables
  character(len=*), parameter :: subname='prepare_psirocc'
  integer :: i_all,i_stat,ierr,ispinor,ispsiw
  integer :: iorb,jorb,ispsir,jproc
  type(workarr_sumrho) :: w
  integer, dimension(:,:), allocatable :: ncommocc
  real(wp), dimension(:), allocatable :: psiwocc

  call initialize_work_arrays_sumrho(lr,w)

  allocate(psiwocc(max(max(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbsocc%norbp,n3parr(0)*orbsocc%norb),1)+ndebug),stat=i_stat)
  call memocc(i_stat,psiwocc,'psiwocc',subname)

  if (lr%geocode == 'F') then
     call razero(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbsocc%norbp,psirocc)
  end if

  !uncompress the wavefunction in the real grid
  !and switch the values of the function

  !occupied orbitals
  ispinor=1
  ispsiw=1
  do iorb=1,orbsocc%norbp
     call daub_to_isf(lr,w,psiocc(1,ispinor,iorb),psirocc(ispsiw))
     ispsir=1+(iorb-1)*n3parr(0)
     do jproc=0,nproc-1
        !write(*,'(a,1x,8(i10))'),'iproc,jproc',iproc,jproc,iorb,orbs%norbp,ispsir,ispsiw,&
        !     lr%d%n1i*lr%d%n2i*max(lr%d%n3i*orbs%norbp,n3p*orbs%norb),n3parr(jproc)
        call dcopy(n3parr(jproc),psirocc(ispsiw),1,psiwocc(ispsir),1)
        ispsiw=ispsiw+n3parr(jproc)
        if (jproc /= nproc-1) then
           do jorb=iorb,orbsocc%norbp
              ispsir=ispsir+n3parr(jproc)
           end do
           do jorb=1,iorb-1
              ispsir=ispsir+n3parr(jproc+1)
           end do
        end if
     end do
  end do
  call deallocate_work_arrays_sumrho(w)

  !communicate them between processors
  !occupied orbitals
  if (nproc > 1) then
     !arrays for the communication between processors
     !valid only for one k-point for the moment
     !and only real functions (nspinor=1)
     !this distribution is in principle valid also for k-points
     allocate(ncommocc(0:nproc-1,4+ndebug),stat=i_stat)
     call memocc(i_stat,ncommocc,'ncommocc',subname)

     !count occay for orbitals => components
     do jproc=0,nproc-1
        ncommocc(jproc,1)=n3parr(jproc)*orbsocc%norb_par(iproc)
     end do
     !displacement array for orbitals => components
     ncommocc(0,2)=0
     do jproc=1,nproc-1
        ncommocc(jproc,2)=ncommocc(jproc-1,2)+ncommocc(jproc-1,1)
     end do
     !count occay for components => orbitals
     do jproc=0,nproc-1
        ncommocc(jproc,3)=n3parr(iproc)*orbsocc%norb_par(jproc)
     end do
     !displacement array for components => orbitals
     ncommocc(0,4)=0
     do jproc=1,nproc-1
        ncommocc(jproc,4)=ncommocc(jproc-1,4)+ncommocc(jproc-1,3)
     end do

     call MPI_ALLTOALLV(psiwocc,ncommocc(0,1),ncommocc(0,2),mpidtypw, &
          psirocc,ncommocc(0,3),ncommocc(0,4),mpidtypw,MPI_COMM_WORLD,ierr)
  else
     call dcopy(lr%d%n1i*lr%d%n2i*n3p*orbsocc%norb,psiwocc,1,psirocc,1)
  end if
  i_all=-product(shape(psiwocc))*kind(psiwocc)
  deallocate(psiwocc,stat=i_stat)
  call memocc(i_stat,i_all,'psiwocc',subname)

  if (nproc > 1) then
     i_all=-product(shape(ncommocc))*kind(ncommocc)
     deallocate(ncommocc,stat=i_stat)
     call memocc(i_stat,i_all,'ncommocc',subname)
  end if

END SUBROUTINE prepare_psirocc
!!***


!!****f* BigDFT/exact_exchange_potential_virt
!! FUNCTION
!!   Calculate the exact exchange potential only on virtual orbitals
!!   by knowing the occupied orbitals and their distribution
!!   both sets of orbitals are to be 
!! SOURCE
!! 
subroutine exact_exchange_potential_virt(iproc,nproc,geocode,nspin,lr,orbsocc,orbsvirt,n3parr,n3p,&
     hxh,hyh,hzh,pkernel,psirocc,psivirt,psirvirt)
  use module_base
  use module_types
  use Poisson_Solver
  implicit none
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: iproc,nproc,n3p,nspin
  real(gp), intent(in) :: hxh,hyh,hzh
  type(locreg_descriptors), intent(in) :: lr
  type(orbitals_data), intent(in) :: orbsocc,orbsvirt
  integer, dimension(0:nproc-1), intent(in) :: n3parr
  real(wp), dimension(max(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbsocc%norbp,n3parr(0)*orbsocc%norb)), intent(in) :: psirocc
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbsvirt%nspinor,orbsvirt%norbp), intent(in) :: psivirt
  real(dp), dimension(*), intent(in) :: pkernel
  real(wp), dimension(max(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbsvirt%norbp,n3parr(0)*orbsvirt%norb)), intent(out) :: psirvirt
  !local variables
  character(len=*), parameter :: subname='exact_exchange_potential_virt'
  integer :: i_all,i_stat,ierr,ispinor,ispsiw,ispin,norbocc,norbvirt
  integer :: i1,i2,i3p,iorb,iorbs,jorb,jorbs,ispsir,ind3,ind2,ind1i,ind1j,jproc,igran,ngran
  real(gp) :: ehart,hfac,sign,sfac,hfaci
  type(workarr_sumrho) :: w
  integer, dimension(:,:), allocatable :: ncommvirt
  real(wp), dimension(:), allocatable :: psiwvirt
  real(wp), dimension(:,:,:,:), allocatable :: rp_ij

  !call timing(iproc,'Exchangecorr  ','ON')

  call initialize_work_arrays_sumrho(lr,w)
  
  !the granularity of the calculation is set by ngran
  !for the moment it is irrelevant but if the poisson solver is modified
  !we may increase this value
  ngran=1

  !partial densities with a given granularity
  allocate(rp_ij(lr%d%n1i,lr%d%n2i,n3p,ngran+ndebug),stat=i_stat)
  call memocc(i_stat,rp_ij,'rp_ij',subname)
  allocate(psiwvirt(max(max(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbsvirt%norbp,n3parr(0)*orbsvirt%norb),1)+ndebug),stat=i_stat)
  call memocc(i_stat,psiwvirt,'psiwvirt',subname)

  if (geocode == 'F') then
     call razero(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbsvirt%norbp,psiwvirt)
  end if

  !uncompress the wavefunction in the real grid
  !and switch the values of the function
  !unoccupied orbitals
  ispinor=1
  ispsiw=1
  do iorb=1,orbsvirt%norbp
     call daub_to_isf(lr,w,psivirt(1,ispinor,iorb),psiwvirt(ispsiw))
     ispsir=1+(iorb-1)*n3parr(0)
     do jproc=0,nproc-1
        !write(*,'(a,1x,8(i10))'),'iproc,jproc',iproc,jproc,iorb,orbs%norbp,ispsir,ispsiw,&
        !     lr%d%n1i*lr%d%n2i*max(lr%d%n3i*orbs%norbp,n3p*orbs%norb),n3parr(jproc)
        call dcopy(n3parr(jproc),psiwvirt(ispsiw),1,psirvirt(ispsir),1)
        ispsiw=ispsiw+n3parr(jproc)
        if (jproc /= nproc-1) then
           do jorb=iorb,orbsvirt%norbp
              ispsir=ispsir+n3parr(jproc)
           end do
           do jorb=1,iorb-1
              ispsir=ispsir+n3parr(jproc+1)
           end do
        end if
     end do
  end do

  call deallocate_work_arrays_sumrho(w)

  !communicate them between processors
  !unoccupied orbitals
  if (nproc > 1) then
     !arrays for the communication between processors
     !valid only for one k-point for the moment
     !and only real functions (nspinor=1)
     !this distribution is in principle valid also for k-points

     allocate(ncommvirt(0:nproc-1,4+ndebug),stat=i_stat)
     call memocc(i_stat,ncommvirt,'ncommvirt',subname)

     !count array for orbitals => components
     do jproc=0,nproc-1
        ncommvirt(jproc,1)=n3parr(jproc)*orbsvirt%norb_par(iproc)
     end do
     !displacement array for orbitals => components
     ncommvirt(0,2)=0
     do jproc=1,nproc-1
        ncommvirt(jproc,2)=ncommvirt(jproc-1,2)+ncommvirt(jproc-1,1)
     end do
     !count array for components => orbitals
     do jproc=0,nproc-1
        ncommvirt(jproc,3)=n3parr(iproc)*orbsvirt%norb_par(jproc)
     end do
     !displacement array for components => orbitals
     ncommvirt(0,4)=0
     do jproc=1,nproc-1
        ncommvirt(jproc,4)=ncommvirt(jproc-1,4)+ncommvirt(jproc-1,3)
     end do

     call MPI_ALLTOALLV(psirvirt,ncommvirt(0,1),ncommvirt(0,2),mpidtypw, &
          psiwvirt,ncommvirt(0,3),ncommvirt(0,4),mpidtypw,MPI_COMM_WORLD,ierr)

  else
     call dcopy(lr%d%n1i*lr%d%n2i*n3p*orbsvirt%norb,psirvirt,1,psiwvirt,1)
  end if

  !this is the array of the actions of the X potential on psi
  call razero(lr%d%n1i*lr%d%n2i*n3p*orbsvirt%norb,psirvirt)

  !build the partial densities for the poisson solver, calculate the partial potential
  !and accumulate the result
  !do it for different spins
  !for non spin-polarised systems there is a factor of two
  !non-collinear spin not yet implemented
  if (nspin==2) then
     sfac=1.0_gp
  else 
     sfac=0.5_gp
  end if

  !number of orbitals, all quantum numbers
  norbocc=orbsocc%norb!*orbs%nkpts, for the future
  norbvirt=orbsvirt%norb
  iorb=1 !index on virtual orbitals
  jorb=1 !index on occupied orbitals
  
  do ispin=1,nspin
     if (ispin==1) then
        iorb=1
        jorb=1
        norbocc=orbsocc%norbu
        norbvirt=orbsvirt%norbu
        sign=1.0_gp
     else
        iorb=orbsvirt%norbu+1
        jorb=orbsocc%norbu+1
        norbocc=orbsocc%norb
        norbvirt=orbsvirt%norb
        sign=-1.0_gp
     end if
     orbital_loop: do
        iorbs=iorb
        jorbs=jorb
        hfac=1/(hxh*hyh*hzh)
        do igran=1,ngran
           if (iorb > norbvirt) exit orbital_loop
           if (orbsvirt%spinsgn(iorb) == sign .and. orbsocc%spinsgn(jorb) == sign) then
              !calculate partial density (real functions), no spin-polarisation
              do i3p=1,n3p
                 ind3=(i3p-1)*lr%d%n1i*lr%d%n2i
                 do i2=1,lr%d%n2i
                    ind2=(i2-1)*lr%d%n1i+ind3
                    do i1=1,lr%d%n1i
                       ind1i=i1+ind2+(iorb-1)*lr%d%n1i*lr%d%n2i*n3p
                       ind1j=i1+ind2+(jorb-1)*lr%d%n1i*lr%d%n2i*n3p
                       rp_ij(i1,i2,i3p,igran)=hfac*psiwvirt(ind1i)*psirocc(ind1j)
                    end do
                 end do
              end do
           end if
           jorb=jorb+1
           if (jorb > norbocc) then
              iorb=iorb+1
              jorb=1
           end if
        end do
        jorb=jorbs
        iorb=iorbs
        do igran=1,ngran
           if (iorb > norbvirt) exit orbital_loop
           if (orbsvirt%spinsgn(iorb) == sign .and. orbsocc%spinsgn(jorb) == sign) then

              !partial exchange term for each partial density
              if (iproc == 0 .and. verbose > 1) then
                 write(*,*)'Exact exchange calculation: spin, orbitals:',ispin,iorb,jorb
              end if
              call H_potential(geocode,'D',iproc,nproc,&
                   lr%d%n1i,lr%d%n2i,lr%d%n3i,hxh,hyh,hzh,&
                   rp_ij(1,1,1,igran),pkernel,rp_ij,ehart,0.0_dp,.false.,&
                   quiet='YES')

!!$              call PSolver(geocode,'D',iproc,nproc,lr%d%n1i,lr%d%n2i,lr%d%n3i,&
!!$                   0,hxh,hyh,hzh,rp_ij(1,1,1,igran),pkernel,rp_ij,ehart,zero,zero,&
!!$                   0.d0,.false.,1,quiet='YES')

           end if
           jorb=jorb+1
           if (jorb > norbocc) then
              iorb=iorb+1
              jorb=1
           end if
        end do
        jorb=jorbs
        iorb=iorbs
        do igran=1,ngran
           if (iorb > norbvirt) exit orbital_loop
           if (orbsvirt%spinsgn(iorb) == sign .and. orbsocc%spinsgn(jorb) == sign) then
              !this factor is only valid with one k-point
              !we have to correct with the kwgts if we want more than one k-point
              hfaci=-sfac*orbsocc%occup(jorb)

              !accumulate the results for each of the wavefunctions concerned
              do i3p=1,n3p
                 ind3=(i3p-1)*lr%d%n1i*lr%d%n2i
                 do i2=1,lr%d%n2i
                    ind2=(i2-1)*lr%d%n1i+ind3
                    do i1=1,lr%d%n1i
                       ind1i=i1+ind2+(iorb-1)*lr%d%n1i*lr%d%n2i*n3p
                       ind1j=i1+ind2+(jorb-1)*lr%d%n1i*lr%d%n2i*n3p
                       psirvirt(ind1i)=psirvirt(ind1i)+hfaci*rp_ij(i1,i2,i3p,igran)*psirocc(ind1j)
                    end do
                 end do
              end do
           end if
           jorb=jorb+1
           if (jorb > norbocc) then
              iorb=iorb+1
              jorb=1
           end if
        end do
     end do orbital_loop
  end do

  !assign the potential for each function
  if (nproc > 1) then
     !call dcopy(lr%d%n1i*lr%d%n2i*n3p*orbs%norb,psir,1,psirt,1)
     !recommunicate the values in the psir array
     call MPI_ALLTOALLV(psirvirt,ncommvirt(0,3),ncommvirt(0,4),mpidtypw, &
          psiwvirt,ncommvirt(0,1),ncommvirt(0,2),mpidtypw,MPI_COMM_WORLD,ierr)
     !redress the potential
     ispsiw=1
     do iorb=1,orbsvirt%norbp
        ispsir=1+(iorb-1)*n3parr(0)
        do jproc=0,nproc-1
           call dcopy(n3parr(jproc),psiwvirt(ispsir),1,psirvirt(ispsiw),1)
           ispsiw=ispsiw+n3parr(jproc)
           if (jproc /= nproc-1) then
              do jorb=iorb,orbsvirt%norbp
                 ispsir=ispsir+n3parr(jproc)
              end do
              do jorb=1,iorb-1
                 ispsir=ispsir+n3parr(jproc+1)
              end do
           end if
        end do
     end do
  end if

  i_all=-product(shape(rp_ij))*kind(rp_ij)
  deallocate(rp_ij,stat=i_stat)
  call memocc(i_stat,i_all,'rp_ij',subname)

  i_all=-product(shape(psiwvirt))*kind(psiwvirt)
  deallocate(psiwvirt,stat=i_stat)
  call memocc(i_stat,i_all,'psiwvirt',subname)


  if (nproc > 1) then
     i_all=-product(shape(ncommvirt))*kind(ncommvirt)
     deallocate(ncommvirt,stat=i_stat)
     call memocc(i_stat,i_all,'ncommvirt',subname)
  end if

end subroutine exact_exchange_potential_virt
!!***
