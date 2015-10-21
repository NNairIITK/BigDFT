module forces_linear

  implicit none
  private

  ! Public routines
  public :: nonlocal_forces_linear
  public :: local_hamiltonian_stress_linear
  public :: erf_stress

  contains

    subroutine erf_stress(at,rxyz,hxh,hyh,hzh,n1i,n2i,n3i,n3p,iproc,nproc,ngatherarr,rho,tens)
      use module_base
      use module_types
      use module_fft_sg
      implicit none
      !passed var
      type(atoms_data), intent(in) :: at
      real(gp), dimension(3,at%astruct%nat), target, intent(in) :: rxyz
      real(gp), intent(in) :: hxh,hyh,hzh
      integer,intent(in) :: n1i,n2i,n3i,n3p,iproc,nproc
      real(kind=8), dimension(n1i*n2i*max(n3p,1)), intent(in), target :: rho
      integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
      real(dp),dimension(6), intent(out) :: tens
      !local var
      character(len=*), parameter :: subname='erf_stress'
      real(kind=8),allocatable :: rhog(:,:,:,:,:)
      real(kind=8),dimension(:),pointer :: rhor
      integer :: ierr
      real(kind=8) :: p(3),g2,rloc,setv,fac
      real(kind=8) :: rx,ry,rz,sfr,sfi,rhore,rhoim
      real(kind=8) :: potg,potg2
      real(kind=8) :: Zion
      integer :: iat,ityp,iout,iiout,nout,noutp,ii,isout,jproc
      integer :: j1,j2,j3,i1,i2,i3,inzee,ind
      real(kind=8) :: routp
      integer,dimension(:),allocatable :: nout_par

      call f_routine(id='erf_stress')

      !write(*,*) 'iproc,n3i,n3p',iproc,n3i,n3p
      !write(*,*) 'iproc',iproc, ngatherarr(iproc-1,1),ngatherarr(iproc-1,2)

      if (nproc > 1) then
         rhor = f_malloc_ptr(n1i*n2i*n3i,id='rhor')
         call MPI_ALLGATHERV(rho(1),ngatherarr(iproc,1),&
              &   mpidtypw,rhor(1),ngatherarr(0,1),&
              ngatherarr(0,2),mpidtypw,bigdft_mpi%mpi_comm,ierr)
      else
         rhor => rho        
      end if

      !  pi = 4.0_gp*atan(1.0_gp)
      rhog = f_malloc((/2,n1i+1,n2i+1,n3i+1,2/),id='rhog')
      tens(:)=0.0_dp ; p(:)=0.0_dp

      ! calculate total rho(G)
      rhog=0.0_dp
      do i3=1,n3i
         do i2=1,n2i
            do i1=1,n1i
               ind=i1+(i2-1)*n1i+(i3-1)*n1i*n2i
               rhog(1,i1,i2,i3,1)=rhor(ind)
               rhog(2,i1,i2,i3,1)=0.d0
            end do
         end do
      end do

      ! DO FFT OF DENSITY: Rho(r) -FFT-> Rho(G)
      inzee=1
      call FFT(n1i,n2i,n3i,n1i+1,n2i+1,n3i+1,rhog,-1,inzee)   
      rhog=rhog/real(n1i*n2i*n3i,kind=8)

      tens=0.0_dp

      ! Collapse the loops over nat and n3i
      nout = at%astruct%nat*n3i
      routp = real(nout,kind=8)/real(nproc,kind=8)
      noutp = floor(routp)
      ii = nout - noutp*nproc
      if (iproc<ii) noutp = noutp + 1
      nout_par = f_malloc0(0.to.nproc-1,id='nout_par')
      nout_par(iproc) = noutp
      if (nproc>1) then
         !call mpiallred(nout_par(0), nproc, mpi_sum, comm=bigdft_mpi%mpi_comm)
         call mpiallred(nout_par,op=mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if
      if (sum(nout_par)/=nout) then
         call f_err_throw('wrong partition of the outer loop',err_name='BIGDT_RUNTIME_ERROR')
      end if
      isout = 0
      do jproc=0,iproc-1
         isout = isout + nout_par(jproc)
      end do

      !do iat=1,at%astruct%nat                          ! SUM OVER ATOMS
      do iout=1,noutp!at%astruct%nat*n3i

         iiout = isout + iout

         iat = (iiout-1)/n3i + 1
         ityp=at%astruct%iatype(iat)                ! ityp
         rloc=at%psppar(0,0,ityp)           ! take corresp. r_loc
         Zion = real(at%nelpsp(ityp),kind=8)

         !coordinates of new center (atom)
         rx=real(rxyz(1,iat),kind=8) 
         ry=real(rxyz(2,iat),kind=8) 
         rz=real(rxyz(3,iat),kind=8)

         potg=0.0_dp

         i3 = modulo(iiout-1,n3i) + 1
         j3=i3-(i3/(n3i/2+2))*n3i-1
         p(3)=real(j3,dp)/(n3i*hzh)         

         !$omp parallel default(none) &
         !$omp shared(n3i, inzee, i3, n2i, hyh, n1i, hxh, rx, ry, rz, rhog, Zion, rloc, hzh, tens) &
         !$omp private(i2, j2, i1, j1, sfr, sfi, rhore, rhoim, g2, fac, setv, potg, potg2) &
         !$omp firstprivate(p)
         !$omp do reduction(+:tens)
         do i2=1,n2i
            j2=i2-(i2/(n2i/2+2))*n2i-1
            p(2)=real(j2,dp)/(n2i*hyh)    

            do i1=1,n1i
               j1=i1-(i1/(n1i/2+2))*n1i-1
               p(1)=real(j1,dp)/(n1i*hxh)

               ! calculate structural factor exp(-2*pi*i*R*G)
               sfr=cos(-2.0_gp*pi*(rx*p(1)+ry*p(2)+rz*p(3)))
               sfi=sin(-2.0_gp*pi*(rx*p(1)+ry*p(2)+rz*p(3)))

               ! multiply density by shift /structural/ factor 
               rhore=rhog(1,i1,i2,i3,inzee)*sfr
               rhoim=rhog(2,i1,i2,i3,inzee)*sfi

               ! use analytical expression for rhog^el
               ! multiply by  Kernel 1/pi/G^2 in reciprocal space

               g2=real(p(1)**2+p(2)**2+p(3)**2,kind=8)

               !set = rhog^el (analytic)
               fac=(Zion/rloc**3.0_gp)/sqrt(2.0_gp*pi)/(2.0_gp*pi)
               fac=fac/real(n1i*hxh*n2i*hyh*n3i*hzh,kind=8)                    !Division by Volume
               setv=((sqrt(pi*2.0_gp*rloc**2.0_gp))**3)*fac*exp(-pi*pi*g2*2.0_gp*rloc**2.0_gp)

               if (g2 /= 0.d0) then

                  potg = -setv/(pi*g2)  ! V^el(G)
                  potg2 = (setv/pi)*((real(1.d0,kind=8)/g2**2.d0)&
                       +real(pi*pi*2.d0*rloc**2.d0/g2,kind=8))

                  !STRESS TENSOR
                  tens(1)=tens(1)-(rhore+rhoim)*(potg2*2.0_gp*p(1)*p(1)+potg)
                  tens(2)=tens(2)-(rhore+rhoim)*(potg2*2.0_gp*p(2)*p(2)+potg)
                  tens(3)=tens(3)-(rhore+rhoim)*(potg2*2.0_gp*p(3)*p(3)+potg)
                  tens(6)=tens(6)-(rhore+rhoim)*(potg2*2.0_gp*p(1)*p(2))
                  tens(5)=tens(5)-(rhore+rhoim)*(potg2*2.0_gp*p(1)*p(3))
                  tens(4)=tens(4)-(rhore+rhoim)*(potg2*2.0_gp*p(2)*p(3))

               end if  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! g2 /=0
            end do !i1
         end do !i2
         !$omp end do
         !$omp end parallel

      end do !iat -atoms

!!$  !not needed as the stress tensor is reduced afterwards
!!$  if (nproc>1) then
!!$      call mpiallred(tens(1), 6, mpi_sum, comm=bigdft_mpi%mpi_comm)
!!$  end if


      call f_free(nout_par)

      if (nproc>1) then
         call f_free_ptr(rhor)
      else
         nullify(rhor)
      end if
      call f_free(rhog)

      call f_release_routine()

    END SUBROUTINE erf_stress


    !> Calculates the nonlocal forces on all atoms arising from the wavefunctions 
    !! belonging to iproc and adds them to the force array
    !! recalculate the projectors at the end if refill flag is .true.
    subroutine nonlocal_forces_linear(iproc,nproc,npsidim_orbs,lr,hx,hy,hz,at,rxyz,&
         orbs,nlpsp,lzd,phi,denskern,denskern_mat,fsep,refill,calculate_strten,strten)
      use module_base
      use module_types
      use sparsematrix_base, only: sparse_matrix, matrices, sparsematrix_malloc, assignment(=), SPARSE_FULL
      use sparsematrix, only: gather_matrix_from_taskgroups
      use psp_projectors, only: projector_has_overlap
      use public_enums, only: PSPCODE_HGH,PSPCODE_HGH_K,PSPCODE_HGH_K_NLCC,&
           PSPCODE_PAW
    
      use yaml_output
      use locregs, only: check_whether_bounds_overlap
      implicit none
      !Arguments-------------
      type(atoms_data), intent(in) :: at
      type(local_zone_descriptors), intent(in) :: lzd
      type(DFT_PSP_projectors), intent(inout) :: nlpsp
      logical, intent(in) :: refill,calculate_strten
      integer, intent(in) :: iproc, nproc, npsidim_orbs
      real(gp), intent(in) :: hx,hy,hz
      type(locreg_descriptors) :: lr
      type(orbitals_data), intent(in) :: orbs
      real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
      real(wp), dimension(npsidim_orbs), intent(in) :: phi
      type(sparse_matrix),intent(in) :: denskern
      type(matrices),intent(inout) :: denskern_mat
      real(gp), dimension(3,at%astruct%nat), intent(inout) :: fsep
      real(gp), dimension(6), intent(out) :: strten
      !local variables--------------
      integer :: istart_c,iproj,iat,ityp,i,j,l,m,iorbout,iiorb,ilr
      integer :: mbseg_c,mbseg_f,jseg_c,jseg_f,ind,iseg,jjorb,ispin
      integer :: mbvctr_c,mbvctr_f,iorb,nwarnings,nspinor,ispinor,jorbd,ncount,ist_send
      real(gp) :: offdiagcoeff,hij,sp0,spi,sp0i,sp0j,spj,Enl,vol
      !real(gp) :: orbfac,strc
      real(gp) :: strc
      integer :: idir,ncplx,icplx,isorb,ikpt,ieorb,istart_ck,ispsi_k,ispsi,jorb,jproc,ii,ist,ierr,iiat,iiiat
      real(gp), dimension(2,2,3) :: offdiagarr
      !real(gp), dimension(:,:), allocatable :: fxyz_orb
      real(dp), dimension(:,:,:,:,:,:,:), allocatable :: scalprod
      !real(gp), dimension(6) :: sab
      integer,dimension(:),allocatable :: nat_par, isat_par!, sendcounts, recvcounts, senddspls, recvdspls
      integer,dimension(:),allocatable :: is_supfun_per_atom, supfun_per_atom, scalprod_send_lookup
      integer,dimension(:),pointer :: scalprod_lookup
      !integer,dimension(:,:),allocatable :: iat_startend
      !real(dp),dimension(:,:,:,:,:,:,:),allocatable :: scalprod_sendbuf
      real(dp),dimension(:,:,:,:,:,:),allocatable :: scalprod_sendbuf_new
      real(dp),dimension(:,:,:,:,:,:),pointer :: scalprod_new
      real(dp),dimension(:),allocatable :: scalprod_recvbuf
      integer :: ndir!=9 !3 for forces, 9 for forces and stresses
      real(kind=8),dimension(:),allocatable :: denskern_gathered
      !integer,dimension(:,:),allocatable :: iorbminmax, iatminmax 
      integer :: iorbmin, jorbmin, iorbmax, jorbmax, nscalprod_send, jat
      integer :: i1s, i1e, j1s, j1e, i2s, i2e, j2s, j2e, i3s, i3e, j3s, j3e
      integer :: nat_per_iteration, isat, natp, iat_out, nat_out, norbp_max
      integer :: l_max, i_max, m_max
      real(kind=8),dimension(:,:),allocatable :: fxyz_orb
    
      !integer,parameter :: MAX_SIZE=268435456 !max size of the array scalprod, in elements
    
      !integer :: ldim, gdim
      !real(8),dimension(:),allocatable :: phiglobal
      !real(8),dimension(2,0:9,7,3,4,at%astruct%nat,orbsglobal%norb) :: scalprodglobal
      !scalprodglobal=0.d0
      !!allocate(phiglobal(lzd%glr%wfd%nvctr_c+7*lzd%glr%wfd%nvctr_f))
    
      real(kind=4) :: tr0, tr1, trt0, trt1
      real(kind=8) :: time0, time1, time2, time3, time4, time5, time6, time7, ttime
      logical, parameter :: extra_timing=.false.
    
    
    
      call f_routine(id='nonlocal_forces_linear')
    
      if (extra_timing) call cpu_time(trt0)
    
    
      !fxyz_orb = f_malloc0((/ 3, at%astruct%nat /),id='fxyz_orb')
    
      ! Gather together the entire density kernel
      denskern_gathered = sparsematrix_malloc(denskern,iaction=SPARSE_FULL,id='denskern_gathered')
      call gather_matrix_from_taskgroups(iproc, nproc, denskern, denskern_mat%matrix_compr, denskern_gathered)
    
      isat = 1
    
    
      natp = at%astruct%nat
    
      ! Determine how many atoms each MPI task will handle
      nat_par = f_malloc(0.to.nproc-1,id='nat_par')
      isat_par = f_malloc(0.to.nproc-1,id='isat_par')
      ii=natp/nproc
      nat_par(0:nproc-1)=ii
      ii=natp-ii*nproc
      do i=0,ii-1
         nat_par(i)=nat_par(i)+1
      end do
      isat_par(0)=0
      do jproc=1,nproc-1
         isat_par(jproc)=isat_par(jproc-1)+nat_par(jproc-1)
      end do
    
    
      ! Number of support functions having an overlap with the projector of a given atom
      supfun_per_atom = f_malloc0(at%astruct%nat,id='supfun_per_atom')
      is_supfun_per_atom = f_malloc0(at%astruct%nat,id='is_supfun_per_atom')
    
    
      call f_zero(strten) 
    
    
      !always put complex scalprod
      !also nspinor for the moment is the biggest as possible
    
    
      Enl=0._gp
      !strten=0.d0
      vol=real(at%astruct%cell_dim(1)*at%astruct%cell_dim(2)*at%astruct%cell_dim(3),gp)
      !sab=0.d0
    
      !calculate the coefficients for the off-diagonal terms
      do l=1,3
         do i=1,2
            do j=i+1,3
               offdiagcoeff=0.0_gp
               if (l==1) then
                  if (i==1) then
                     if (j==2) offdiagcoeff=-0.5_gp*sqrt(3._gp/5._gp)
                     if (j==3) offdiagcoeff=0.5_gp*sqrt(5._gp/21._gp)
                  else
                     offdiagcoeff=-0.5_gp*sqrt(100._gp/63._gp)
                  end if
               else if (l==2) then
                  if (i==1) then
                     if (j==2) offdiagcoeff=-0.5_gp*sqrt(5._gp/7._gp)
                     if (j==3) offdiagcoeff=1._gp/6._gp*sqrt(35._gp/11._gp)
                  else
                     offdiagcoeff=-7._gp/3._gp*sqrt(1._gp/11._gp)
                  end if
               else if (l==3) then
                  if (i==1) then
                     if (j==2) offdiagcoeff=-0.5_gp*sqrt(7._gp/9._gp)
                     if (j==3) offdiagcoeff=0.5_gp*sqrt(63._gp/143._gp)
                  else
                     offdiagcoeff=-9._gp*sqrt(1._gp/143._gp)
                  end if
               end if
               offdiagarr(i,j-i,l)=offdiagcoeff
            end do
         end do
      end do
    
    
      ! Determine the maximal value of l and i
      l_max = 1
      i_max = 1
      do iat=1,natp
         iiat = iat+isat-1
         ityp=at%astruct%iatype(iiat)
         do l=1,4
            do i=1,3
               if (at%psppar(l,i,ityp) /= 0.0_gp) then
                  l_max = max(l,l_max)
                  i_max = max(i,i_max)
               end if
            end do
         end do
      end do
      m_max = 2*l_max-1
    
    
      ! Determine the size of the array scalprod_sendbuf (indicated by iat_startend)
      call determine_dimension_scalprod(calculate_strten, natp, isat, at, lzd, nlpsp, &
               orbs, supfun_per_atom, ndir, nscalprod_send)
      scalprod_sendbuf_new = f_malloc0((/ 1.to.2, 0.to.ndir, 1.to.m_max, 1.to.i_max, 1.to.l_max, &
           1.to.max(nscalprod_send,1) /),id='scalprod_sendbuf_new')
      scalprod_send_lookup = f_malloc(max(nscalprod_send,1), id='scalprod_send_lookup')
    
      is_supfun_per_atom(1) = 0
      do jat=2,at%astruct%nat
         is_supfun_per_atom(jat) = is_supfun_per_atom(jat-1) + supfun_per_atom(jat-1)
      end do
    
      ! Calculate the values of scalprod
    
      if (extra_timing) call cpu_time(tr0)
      !call calculate_scalprod()
      call  calculate_scalprod(iproc, natp, isat, ndir, i_max, l_max, m_max, npsidim_orbs, orbs, &
            lzd, nlpsp, at, lr, hx, hy, hz, is_supfun_per_atom, phi, &
            size(scalprod_send_lookup), scalprod_send_lookup, scalprod_sendbuf_new)
      if (extra_timing) call cpu_time(tr1)
      if (extra_timing) time0=real(tr1-tr0,kind=8)
    
    
      ! Communicate scalprod
      !call transpose_scalprod()
      call transpose_scalprod(iproc, nproc, at, nat_par, isat_par, &
               ndir, i_max, l_max, m_max, scalprod_sendbuf_new, &
               supfun_per_atom, is_supfun_per_atom, &
               size(scalprod_send_lookup), scalprod_send_lookup, scalprod_new, scalprod_lookup)
      if (extra_timing) call cpu_time(tr0)
      if (extra_timing) time1=real(tr0-tr1,kind=8)
    
      !allocate the temporary array
      fxyz_orb = f_malloc0((/3,nat_par(iproc)/),id='fxyz_orb')
    
      !call calculate_forces(nat_par(iproc),fxyz_orb)
      call calculate_forces_kernel(iproc, nproc, ndir, nat_par(iproc), isat, i_max, l_max, m_max, isat_par, orbs, &
               denskern, at, offdiagarr, denskern_gathered, scalprod_new, supfun_per_atom, &
               is_supfun_per_atom, calculate_strten, vol, size(scalprod_lookup), scalprod_lookup, &
               fxyz_orb, Enl, strten, fsep)
      call f_free(fxyz_orb)
    
      if (extra_timing) call cpu_time(tr1)
      if (extra_timing) time2=real(tr1-tr0,kind=8)
    
    
    
      call f_free(is_supfun_per_atom)
      call f_free(supfun_per_atom)
      call f_free(scalprod_sendbuf_new)
      call f_free(scalprod_send_lookup)
      call f_free_ptr(scalprod_lookup)
    
    
      !!call mpiallred(Enl,1,mpi_sum, bigdft_mpi%mpi_comm)
      !!if (bigdft_mpi%iproc==0) call yaml_map('Enl',Enl)
    
      call f_free_ptr(scalprod_new)
      call f_free(nat_par)
      call f_free(isat_par)
      !call f_free(fxyz_orb)
      call f_free(denskern_gathered)
    
      if (extra_timing) call cpu_time(trt1)
      if (extra_timing) ttime=real(trt1-trt0,kind=8)
    
      if (extra_timing.and.iproc==0) print*,'nonloc (scal, trans, calc):',time0,time1,time2,time0+time1+time2,ttime
    
    
    
      call f_release_routine()
    
    
    END SUBROUTINE nonlocal_forces_linear
    
    
    
    subroutine determine_dimension_scalprod(calculate_strten, natp, isat, at, lzd, nlpsp, &
               orbs, supfun_per_atom, ndir, nscalprod_send)
      use module_base
      use module_types, only: atoms_data, local_zone_descriptors, &
                              DFT_PSP_projectors, orbitals_data
      use psp_projectors, only: projector_has_overlap
      implicit none
    
      ! Calling arguments
      logical,intent(in) :: calculate_strten
      integer,intent(in) :: natp, isat
      type(atoms_data),intent(in) :: at
      type(local_zone_descriptors),intent(in) :: lzd
      type(DFT_PSP_projectors), intent(in) :: nlpsp
      type(orbitals_data),intent(in) :: orbs
      integer,dimension(at%astruct%nat),intent(inout) :: supfun_per_atom 
      integer,intent(out) :: ndir, nscalprod_send
      !logical :: projector_has_overlap

      ! Local variables
      integer :: ikpt, ii, isorb, ieorb, iorb, iiorb, ilr, nspinor, iat, iiat, ityp 

      call f_routine(id='determine_dimension_scalprod')

      if (calculate_strten) then
         ndir=9
      else
         ndir=3
      end if
    
      nscalprod_send = 0
      norbp_if: if (orbs%norbp>0) then
    
         !look for the strategy of projectors application
         if (DistProjApply) then
            !apply the projectors on the fly for each k-point of the processor
            !starting k-point
            ikpt=orbs%iokpt(1)
            loop_kptD: do
    
               call orbs_in_kpt(ikpt,orbs,isorb,ieorb,nspinor)
    
    
               do iat=1,natp
                  iiat = iat+isat-1
    
                  ityp=at%astruct%iatype(iiat)
                  ii = 0
                  do iorb=isorb,ieorb
                     iiorb=orbs%isorb+iorb
                     ilr=orbs%inwhichlocreg(iiorb)
                     ! Check whether there is an overlap between projector and support functions
                     if (.not.projector_has_overlap(iat, ilr, lzd%llr(ilr), lzd%glr, nlpsp)) then
                        cycle 
                     else
                        ii = ii +1
                     end if
                  end do
                  nscalprod_send = nscalprod_send + ii
                  supfun_per_atom(iat) = supfun_per_atom(iat) + ii
    
               end do
    
               if (ieorb == orbs%norbp) exit loop_kptD
               ikpt=ikpt+1
            end do loop_kptD
    
         else
    
            stop 'carefully test this section...'
    !calculate all the scalar products for each direction and each orbitals
    
            !!   !apply the projectors  k-point of the processor
            !!   !starting k-point
            !!   ikpt=orbs%iokpt(1)
            !!   loop_kpt: do
    
            !!      call orbs_in_kpt(ikpt,orbs,isorb,ieorb,nspinor)
    
            !!      do iorb=isorb,ieorb
            !!         iiorb=orbs%isorb+iorb
            !!         ilr=orbs%inwhichlocreg(iiorb)
            !!         ! Quick check
            !!         if (lzd%llr(ilr)%ns1>nlpsp%pspd(iiat)%plr%ns1+nlpsp%pspd(iiat)%plr%d%n1 .or. &
            !!             nlpsp%pspd(iiat)%plr%ns1>lzd%llr(ilr)%ns1+lzd%llr(ilr)%d%n1 .or. &
            !!             lzd%llr(ilr)%ns2>nlpsp%pspd(iiat)%plr%ns2+nlpsp%pspd(iiat)%plr%d%n2 .or. &
            !!             nlpsp%pspd(iiat)%plr%ns2>lzd%llr(ilr)%ns2+lzd%llr(ilr)%d%n2 .or. &
            !!             lzd%llr(ilr)%ns3>nlpsp%pspd(iiat)%plr%ns3+nlpsp%pspd(iiat)%plr%d%n3 .or. &
            !!             nlpsp%pspd(iiat)%plr%ns3>lzd%llr(ilr)%ns3+lzd%llr(ilr)%d%n3) then
            !!             cycle 
            !!         else
            !!             iat_startend(1,iproc) = min(iat_startend(1,iproc),iat)
            !!             iat_startend(2,iproc) = max(iat_startend(2,iproc),iat)
            !!         end if
            !!      end do
            !!      if (ieorb == orbs%norbp) exit loop_kpt
            !!      ikpt=ikpt+1
            !!   end do loop_kpt
    
         end if
    
      else norbp_if
    
         !!iat_startend(1,iproc) = 1
         !!iat_startend(2,iproc) = 1
    
      end if norbp_if
    
      !if (nproc>1) then
      !   call mpiallred(iat_startend, mpi_sum, bigdft_mpi%mpi_comm)
      !end if
    
      call f_release_routine()
    
    end subroutine determine_dimension_scalprod
    
    
    subroutine calculate_scalprod(iproc, natp, isat, ndir, i_max, l_max, m_max, npsidim_orbs, orbs, &
               lzd, nlpsp, at, lr, hx, hy, hz, is_supfun_per_atom, phi, &
               nscalprod_send, scalprod_send_lookup, scalprod_sendbuf_new)
      use module_base
      use module_types, only: orbitals_data, local_zone_descriptors, &
                              atoms_data, DFT_PSP_projectors
      use psp_projectors, only: projector_has_overlap
      use locregs, only: locreg_descriptors
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, natp, isat, ndir, i_max, l_max, m_max, npsidim_orbs
      integer,intent(in) :: nscalprod_send
      type(orbitals_data),intent(in) :: orbs
      type(local_zone_descriptors),intent(in) :: lzd
      type(DFT_PSP_projectors),intent(in) :: nlpsp
      type(atoms_data),intent(in) :: at
      type(locreg_descriptors),intent(in) :: lr
      real(kind=8),intent(in) :: hx, hy, hz
      integer,dimension(at%astruct%nat),intent(in) :: is_supfun_per_atom
      real(wp), dimension(npsidim_orbs),intent(in) :: phi
      integer,dimension(nscalprod_send),intent(inout) :: scalprod_send_lookup
      real(kind=8),dimension(1:2,0:ndir,1:m_max,1:i_max,1:l_max,1:nscalprod_send),intent(inout) ::  scalprod_sendbuf_new
    
      ! Local variables
      integer :: ikpt, ispsi, ispsi_k, jorb, jorbd, iii, nwarnings, iproj, nspinor
      integer :: iat, isorb, ieorb, iorb, iiorb, ilr, ityp, istart_c, i, m, l, ispinor
      integer :: mbseg_c, mbseg_f, mbvctr_c, mbvctr_f, jseg_c, jseg_f, idir, ncplx, iiat
      logical :: increase
      integer,dimension(:),allocatable :: is_supfun_per_atom_tmp
      real(kind=8) :: scpr
      logical :: need_proj
      real(kind=4) :: tr0, tr1, trt0, trt1
      real(kind=8) :: time0, time1, time2, time3, time4, time5, time6, time7, ttime
      logical, parameter :: extra_timing=.false.

    
      call f_routine(id='calculate_scalprod')
    
      time0=0.0d0
      time1=0.0d0
    
      if (extra_timing) call cpu_time(trt0)
    
      is_supfun_per_atom_tmp = f_malloc(at%astruct%nat,id='is_supfun_per_atom_tmp')
    
      norbp_if: if (orbs%norbp>0) then
    
         !look for the strategy of projectors application
         if (DistProjApply) then
            !apply the projectors on the fly for each k-point of the processor
            !starting k-point
            ikpt=orbs%iokpt(1)
            ispsi_k=1
            jorb=0
            loop_kptD: do
    
               call orbs_in_kpt(ikpt,orbs,isorb,ieorb,nspinor)
    
               call ncplx_kpt(ikpt,orbs,ncplx)
    
               nwarnings=0 !not used, simply initialised 
               iproj=0 !should be equal to four times nproj at the end
               jorbd=jorb
               do iat=1,natp
                  iiat = iat+isat-1
    
                  !first check if we actually need this projector (i.e. if any of our orbs overlap with it)
                  need_proj=.false.
                  do iorb=isorb,ieorb
                     iiorb=orbs%isorb+iorb
                     ilr=orbs%inwhichlocreg(iiorb)
                     ! Check whether there is an overlap between projector and support functions
                     if (projector_has_overlap(iat, ilr, lzd%llr(ilr), lzd%glr, nlpsp)) then
                        need_proj=.true.
                        exit
                     end if
                  end do
    
                  if (.not. need_proj) cycle
    
    
                  call plr_segs_and_vctrs(nlpsp%pspd(iiat)%plr,&
                       mbseg_c,mbseg_f,mbvctr_c,mbvctr_f)
                  jseg_c=1
                  jseg_f=1
    
    
                  do idir=0,ndir
    
                     call vcopy(at%astruct%nat, is_supfun_per_atom(1), 1, is_supfun_per_atom_tmp(1), 1)
    
                     ityp=at%astruct%iatype(iiat)
                     !calculate projectors
                     istart_c=1
    
    
                     if (extra_timing) call cpu_time(tr0)
                     call atom_projector(nlpsp, ityp, iiat, at%astruct%atomnames(ityp), &
                          & at%astruct%geocode, idir, lr, hx, hy, hz, &
                          & orbs%kpts(1,ikpt), orbs%kpts(2,ikpt), orbs%kpts(3,ikpt), &
                          & istart_c, iproj, nwarnings)
                     if (extra_timing) call cpu_time(tr1)
                     if (extra_timing) time0=time0+real(tr1-tr0,kind=8)
                     if (extra_timing) call cpu_time(tr0)
                     !calculate the contribution for each orbital
                     !here the nspinor contribution should be adjusted
                     ! loop over all my orbitals
                     ispsi=ispsi_k
                     jorb=jorbd
                     do iorb=isorb,ieorb
                        iiorb=orbs%isorb+iorb
                        ilr=orbs%inwhichlocreg(iiorb)
                        ! Check whether there is an overlap between projector and support functions
                        if (.not.projector_has_overlap(iat, ilr, lzd%llr(ilr), lzd%glr, nlpsp)) then
                           jorb=jorb+1
                           ispsi=ispsi+(lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f)*ncplx
                           cycle 
                        end if
                        increase = .true.
                        do ispinor=1,nspinor,ncplx
                           jorb=jorb+1
                           istart_c=1
                           do l=1,l_max!4
                              do i=1,i_max!3
                                 if (at%psppar(l,i,ityp) /= 0.0_gp) then
                                    do m=1,2*l-1
                                       call wpdot_wrap(ncplx,&
                                            lzd%llr(ilr)%wfd%nvctr_c,lzd%llr(ilr)%wfd%nvctr_f,&
                                            lzd%llr(ilr)%wfd%nseg_c,lzd%llr(ilr)%wfd%nseg_f,&
                                            lzd%llr(ilr)%wfd%keyvglob,lzd%llr(ilr)%wfd%keyglob,phi(ispsi),&
                                            mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
                                            nlpsp%pspd(iiat)%plr%wfd%keyvglob(jseg_c),&
                                            nlpsp%pspd(iiat)%plr%wfd%keyglob(1,jseg_c),&
                                            nlpsp%proj(istart_c),&
                                            scpr)
                                       !if (scpr/=0.d0) then
                                       ! SM: In principle it would be sufficient to update only is_supfun_per_atom_tmp
                                       ! and then put iii =  is_supfun_per_atom_tmp(iat) after the if, but this
                                       ! causes a crash with gfortran
                                       iii = is_supfun_per_atom_tmp(iat)
                                       if (increase) then
                                          is_supfun_per_atom_tmp(iat) = is_supfun_per_atom_tmp(iat)+1
                                          iii = iii + 1
                                          increase = .false.
                                       end if
                                       scalprod_sendbuf_new(1,idir,m,i,l,iii) = scpr
                                       scalprod_send_lookup(iii) = iiorb
                                       !else
                                       !    stop 'scalprod should not be zero'
                                       !end if
                                       istart_c=istart_c+(mbvctr_c+7*mbvctr_f)*ncplx
                                    end do
                                 end if
                              end do
                           end do
                           ispsi=ispsi+(lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f)*ncplx
                        end do
                     end do
                     if (extra_timing) call cpu_time(tr1)
                     if (extra_timing) time1=time1+real(tr1-tr0,kind=8)
    
                     if (istart_c-1  > nlpsp%nprojel) stop '2:applyprojectors'
                  end do
               end do
    
               if (ieorb == orbs%norbp) exit loop_kptD
               ikpt=ikpt+1
               ispsi_k=ispsi
            end do loop_kptD
    
         else
            stop "shouldn't enter this section"
         end if
    
      end if norbp_if
    
    
      call f_free(is_supfun_per_atom_tmp)
    
      if (extra_timing) call cpu_time(trt1)
      if (extra_timing) ttime=real(trt1-trt0,kind=8)
    
      if (extra_timing.and.iproc==0) print*,'calcscalprod (proj, wpdot):',time0,time1,time0+time1,ttime
    
      call f_release_routine()
    
    end subroutine calculate_scalprod
    
    
    
    subroutine transpose_scalprod(iproc, nproc, at, nat_par, isat_par, &
               ndir, i_max, l_max, m_max, scalprod_sendbuf_new, &
               supfun_per_atom, is_supfun_per_atom, &
               nscalprod_send, scalprod_send_lookup, scalprod_new, scalprod_lookup)
      use module_base
      use module_types, only: atoms_data
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, ndir, i_max, l_max, m_max, nscalprod_send
      integer,dimension(0:nproc-1),intent(in) :: nat_par, isat_par
      type(atoms_data),intent(in) :: at
      real(kind=8),dimension(1:2,0:ndir,1:m_max,1:i_max,1:l_max,1:nscalprod_send),intent(in) ::  scalprod_sendbuf_new
      integer,dimension(at%astruct%nat),intent(inout) :: supfun_per_atom 
      integer,dimension(at%astruct%nat),intent(inout) :: is_supfun_per_atom
      integer,dimension(nscalprod_send),intent(inout) :: scalprod_send_lookup
      real(kind=8),dimension(:,:,:,:,:,:),pointer :: scalprod_new
      integer,dimension(:),pointer :: scalprod_lookup
    
      ! Local variables
      integer :: window, iatmin, iatmax, is, ie, nat_on_task, ist_recv, nsize, size_of_double, jjat
      integer :: isrc, idst, ncount, nel, iat, ii, iiat, jat, jproc, iorb, nscalprod_recv
      integer,dimension(:),allocatable :: datatypes, is_supfun_per_atom_tmp, scalprod_lookup_recvbuf
      integer,dimension(:),allocatable :: nsendcounts, nsenddspls, nrecvcounts, nrecvdspls, supfun_per_atom_recv
      integer,dimension(:),allocatable :: nsendcounts_tmp, nsenddspls_tmp, nrecvcounts_tmp, nrecvdspls_tmp
      real(kind=8),dimension(:),allocatable :: scalprod_recvbuf
    
      call f_routine(id='transpose_scalprod')
    
    
    
      ! Prepare communication arrays for alltoallv
      nsendcounts = f_malloc0(0.to.nproc-1,id='nsendcounts')
      nsenddspls = f_malloc(0.to.nproc-1,id='nsenddspls')
      nrecvcounts = f_malloc(0.to.nproc-1,id='nrecvcounts')
      nrecvdspls = f_malloc(0.to.nproc-1,id='nrecvdspls')
      do jproc=0,nproc-1
         do jat=1,nat_par(jproc)
            jjat = isat_par(jproc) + jat
            nsendcounts(jproc) = nsendcounts(jproc) + supfun_per_atom(jjat)
         end do
      end do
      nsenddspls(0) = 0
      do jproc=1,nproc-1
         nsenddspls(jproc) = nsenddspls(jproc-1) + nsendcounts(jproc-1)
      end do
    
    
      ! Communicate the send counts and receive displacements using another alltoallv
      nsendcounts_tmp = f_malloc(0.to.nproc-1,id='nsendcounts_tmp')
      nsenddspls_tmp = f_malloc(0.to.nproc-1,id='nsenddspls_tmp')
      nrecvcounts_tmp = f_malloc(0.to.nproc-1,id='nrecvcounts_tmp')
      nrecvdspls_tmp = f_malloc(0.to.nproc-1,id='nrecvdspls_tmp')
      do jproc=0,nproc-1
         nsendcounts_tmp(jproc) = 1
         nsenddspls_tmp(jproc) = jproc
         nrecvcounts_tmp(jproc) = 1
         nrecvdspls_tmp(jproc) = jproc
      end do
      if (nproc>1) then
         call mpialltoallv(nsendcounts(0), nsendcounts_tmp, nsenddspls_tmp, &
              nrecvcounts(0), nrecvcounts_tmp, nrecvdspls_tmp, bigdft_mpi%mpi_comm)
      else
         call f_memcpy(n=nsendcounts_tmp(0), src=nsendcounts(0), dest=nrecvcounts(0))
      end if
      nrecvdspls(0) = 0
      do jproc=1,nproc-1
         nrecvdspls(jproc) = nrecvdspls(jproc-1) + nrecvcounts(jproc-1)
      end do
    
      ! Communicate the number of scalprods per atom, which will be needed as a switch
      supfun_per_atom_recv = f_malloc(max(nat_par(iproc),1)*nproc,id='supfun_per_atom_recv')
      do jproc=0,nproc-1
         nsendcounts_tmp(jproc) = nat_par(jproc)
         nsenddspls_tmp(jproc) = isat_par(jproc)
         nrecvcounts_tmp(jproc) = nat_par(iproc)
         nrecvdspls_tmp(jproc) = jproc*nat_par(iproc)
      end do
      if (nproc>1) then
         call mpialltoallv(supfun_per_atom(1), nsendcounts_tmp, nsenddspls_tmp, &
              supfun_per_atom_recv(1), nrecvcounts_tmp, nrecvdspls_tmp, bigdft_mpi%mpi_comm)
      else
         call f_memcpy(n=nsendcounts_tmp(0), src=supfun_per_atom(1), dest=supfun_per_atom_recv(1))
      end if
    
      ! Determine the size of the receive buffer
      nscalprod_recv = sum(nrecvcounts)
      scalprod_new = f_malloc_ptr((/ 1.to.2, 0.to.ndir, 1.to.m_max, 1.to.i_max, 1.to.l_max, &
           1.to.max(nscalprod_recv,1) /),id='scalprod_new')
      scalprod_recvbuf = f_malloc(2*(ndir+1)*m_max*i_max*l_max*max(nscalprod_recv,1),id='scalprod_recvbuf')
      scalprod_lookup_recvbuf = f_malloc(max(nscalprod_recv,1), id='scalprod_send_lookup_recvbuf')
      scalprod_lookup = f_malloc_ptr(max(nscalprod_recv,1), id='scalprod_lookup')
    
      ! Communicate the lookup array
      if (nproc>1) then
         call mpialltoallv(scalprod_send_lookup(1), nsendcounts, nsenddspls, &
              scalprod_lookup_recvbuf(1), nrecvcounts, nrecvdspls, bigdft_mpi%mpi_comm)
      else
         call f_memcpy(n=nsendcounts(0), src=scalprod_send_lookup(1), dest=scalprod_lookup_recvbuf(1))
      end if
    
      ! Communicate the scalprods
      ncount = 2*(ndir+1)*m_max*i_max*l_max
      nsendcounts(:) = nsendcounts(:)*ncount
      nsenddspls(:) = nsenddspls(:)*ncount
      nrecvcounts(:) = nrecvcounts(:)*ncount
      nrecvdspls(:) = nrecvdspls(:)*ncount
      if (nproc>1) then
         call mpialltoallv(scalprod_sendbuf_new(1,0,1,1,1,1), nsendcounts, nsenddspls, &
              scalprod_recvbuf(1), nrecvcounts, nrecvdspls, bigdft_mpi%mpi_comm)
      else
         call f_memcpy(n=nsendcounts(0), src=scalprod_sendbuf_new(1,0,1,1,1,1), dest=scalprod_recvbuf(1))
      end if
    
    
    
      ! Now the value of supfun_per_atom can be summed up, since each task has all scalprods for a given atom.
      ! In principle this is only necessary for the atoms handled by a given task, but do it for the moment for all...
      if (nproc>1) then
         call mpiallred(supfun_per_atom, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if
      ! The starting points have to be recalculated.
      is_supfun_per_atom(1) = 0
      do jat=2,at%astruct%nat
         is_supfun_per_atom(jat) = is_supfun_per_atom(jat-1) + supfun_per_atom(jat-1)
      end do
    
    
      ! Rearrange the elements
      is_supfun_per_atom_tmp = f_malloc(at%astruct%nat,id='is_supfun_per_atom_tmp')
      !ncount = 2*(ndir+1)*m_max*i_max*l_max
      call vcopy(at%astruct%nat, is_supfun_per_atom(1), 1, is_supfun_per_atom_tmp(1), 1)
      ii = 1
      isrc = 1
      do jproc=0,nproc-1
         do iat=1,nat_par(iproc)
            iiat = isat_par(iproc) + iat
            nel = supfun_per_atom_recv(ii)
            idst = is_supfun_per_atom_tmp(iiat) - is_supfun_per_atom(isat_par(iproc)+1) + 1
            if (nel>0) then
               call vcopy(nel*ncount, scalprod_recvbuf((isrc-1)*ncount+1), 1, scalprod_new(1,0,1,1,1,idst), 1)
               call vcopy(nel, scalprod_lookup_recvbuf(isrc), 1, scalprod_lookup(idst), 1)
            end if
            is_supfun_per_atom_tmp(iiat) = is_supfun_per_atom_tmp(iiat) + nel
            isrc = isrc + nel
            ii = ii + 1
         end do
      end do
    
      call f_free(is_supfun_per_atom_tmp)
    
    
      do iat=1,nat_par(iproc)
         iiat = isat_par(iproc) + iat
         do iorb=1,supfun_per_atom(iiat)
            ii = is_supfun_per_atom(iiat)+iorb
         end do
      end do
    
      call f_free(scalprod_recvbuf)
      call f_free(nsendcounts)
      call f_free(nsenddspls)
      call f_free(nrecvcounts)
      call f_free(nrecvdspls)
      call f_free(nsendcounts_tmp)
      call f_free(nsenddspls_tmp)
      call f_free(nrecvcounts_tmp)
      call f_free(nrecvdspls_tmp)
      call f_free(scalprod_lookup_recvbuf)
      call f_free(supfun_per_atom_recv)
    
      call f_release_routine()
    
    end subroutine transpose_scalprod
    
    
    
    subroutine calculate_forces_kernel(iproc, nproc, ndir, ntmp, isat, i_max, l_max, m_max, isat_par, orbs, &
               denskern, at, offdiagarr, denskern_gathered, scalprod_new, supfun_per_atom, &
               is_supfun_per_atom, calculate_strten, vol, nscalprod_recv, scalprod_lookup, &
               fxyz_orb, Enl, strten, fsep)
      use module_base
      use module_types, only: atoms_data, orbitals_data
      use sparsematrix_base, only: sparse_matrix
      use sparsematrix_init, only: matrixindex_in_compressed
      use public_enums, only: PSPCODE_HGH, PSPCODE_HGH_K, PSPCODE_HGH_K_NLCC
      implicit none
    
      ! Calling arguments
      integer, intent(in) :: iproc, nproc, ndir, ntmp, isat, i_max, l_max, m_max, nscalprod_recv
      integer,dimension(0:nproc-1),intent(in) :: isat_par
      type(orbitals_data),intent(in) :: orbs
      type(sparse_matrix),intent(in) :: denskern
      type(atoms_data),intent(in) :: at
      real(gp),dimension(2,2,3),intent(in) :: offdiagarr
      real(kind=8),dimension(denskern%nvctr*denskern%nspin),intent(in) :: denskern_gathered
      real(kind=8),dimension(1:2,0:ndir,1:m_max,1:i_max,1:l_max,1:nscalprod_recv),intent(in) :: scalprod_new
      integer,dimension(at%astruct%nat),intent(in) :: supfun_per_atom 
      integer,dimension(at%astruct%nat),intent(in) :: is_supfun_per_atom
      logical,intent(in) :: calculate_strten
      real(gp),intent(in) :: vol
      integer,dimension(nscalprod_recv),intent(in) :: scalprod_lookup
      real(gp), dimension(3,ntmp), intent(inout) :: fxyz_orb
      real(gp),intent(inout) :: Enl
      real(gp),dimension(6),intent(inout) :: strten
      real(gp),dimension(3,at%astruct%nat),intent(inout) :: fsep
    
      ! Local variables
      integer :: jj, iispin, jjspin, ikpt, isorb, ieorb, iorb, ispinor, ncplx, ityp, idir, ispsi, iat
      integer :: iiat, iiorb, ind, ispin, j, nspinor, jorb, i, l, m, jjorb, ii, icplx
      !real(kind=8),dimension(:,:),allocatable :: fxyz_orb
      !real(kind=8),dimension(:),allocatable :: sab, strten_loc
      real(kind=8),dimension(6) :: sab, strten_loc
      real(kind=8) :: tt, tt1, sp0, spi, spj
      real(gp) :: hij, sp0i, sp0j, strc
    
      call f_routine(id='calculate_forces')
    
      !fxyz_orb = f_malloc0((/3,nat_par(iproc)/),id='fxyz_orb')
      !sab = f_malloc0(6,id='sab')
      !strten_loc = f_malloc(6,id='strten_loc')
    
      natp_if2: if (ntmp>0) then
    
         !apply the projectors  k-point of the processor
         !starting k-point
         ikpt=orbs%iokpt(1)
         jorb=0
         loop_kptF2: do
    
            call orbs_in_kpt(ikpt,orbs,isorb,ieorb,nspinor)
    
            call ncplx_kpt(ikpt,orbs,ncplx)
            strten_loc(:) = 0.d0
    
            ! Do the OMP loop over supfun_per_atom, as ntmp is typically rather small
    
            !$omp parallel default(none) &
            !$omp shared(denskern, ntmp, iproc, isat_par, at, supfun_per_atom, is_supfun_per_atom) &
            !$omp shared(scalprod_lookup, l_max, i_max, scalprod_new, fxyz_orb, denskern_gathered) &
            !$omp shared(offdiagarr, strten, strten_loc, vol, Enl, nspinor,ncplx,ndir,calculate_strten) &
            !$omp private(ispin, iat, iiat, ityp, iorb, ii, iiorb, jorb, jj, jjorb, ind, sab, ispinor) &
            !$omp private(l, i, m, icplx, sp0, idir, spi, strc, j, hij, sp0i, sp0j, spj, iispin, jjspin, tt, tt1)
            spin_loop2: do ispin=1,denskern%nspin
    
    
               do iat=1,ntmp
                  iiat=isat_par(iproc)+iat
                  ityp=at%astruct%iatype(iiat)
                  !$omp do reduction(+:fxyz_orb,strten_loc,Enl)
                  do iorb=1,supfun_per_atom(iiat)
                     ii = is_supfun_per_atom(iiat) - is_supfun_per_atom(isat_par(iproc)+1) + iorb
                     iiorb = scalprod_lookup(ii)
                     iispin = (iiorb-1)/denskern%nfvctr + 1
                     if (iispin/=ispin) cycle
                     !  if (ispin==2) then
                     !      ! spin shift
                     !      iiorb = iiorb + denskern%nfvctr
                     !  end if
                     do jorb=1,supfun_per_atom(iiat)
                        jj = is_supfun_per_atom(iiat) - is_supfun_per_atom(isat_par(iproc)+1) + jorb
                        jjorb = scalprod_lookup(jj)
                        jjspin = (jjorb-1)/denskern%nfvctr + 1
                        if (jjspin/=ispin) cycle
                        !     if (ispin==2) then
                        !         !spin shift
                        !         jjorb = jjorb + denskern%nfvctr
                        !     end if
                        ind = matrixindex_in_compressed(denskern, jjorb, iiorb)
                        if (ind==0) cycle
                        sab=0.0_gp
                        ! Loop over all projectors
                        do ispinor=1,nspinor,ncplx
                           if (denskern_gathered(ind)==0.d0) cycle
                           do l=1,l_max!4
                              do i=1,i_max!3
                                 if (at%psppar(l,i,ityp) /= 0.0_gp) then
                                    tt=denskern_gathered(ind)*at%psppar(l,i,ityp)
                                    do m=1,2*l-1
                                       do icplx=1,ncplx
                                          ! scalar product with the derivatives in all the directions
                                          sp0=real(scalprod_new(icplx,0,m,i,l,ii),gp)
                                          tt1=tt*sp0
                                          do idir=1,3
                                             spi=real(scalprod_new(icplx,idir,m,i,l,jj),gp)
                                             fxyz_orb(idir,iat)=fxyz_orb(idir,iat)+&
                                                  tt1*spi
                                          end do
                                          spi=real(scalprod_new(icplx,0,m,i,l,jj),gp)
                                          Enl=Enl+tt1*spi
                                          do idir=4,ndir !for stress
                                             strc=real(scalprod_new(icplx,idir,m,i,l,jj),gp)
                                             sab(idir-3) = sab(idir-3)+&   
                                                  tt1*2.0_gp*strc
                                          end do
                                       end do
                                    end do
                                 end if
                              end do
                           end do
                           !HGH case, offdiagonal terms
                           if (at%npspcode(ityp) == PSPCODE_HGH .or. &
                                at%npspcode(ityp) == PSPCODE_HGH_K .or. &
                                at%npspcode(ityp) == PSPCODE_HGH_K_NLCC) then
                              do l=1,3!min(l_max,3) !no offdiagoanl terms for l=4 in HGH-K case
                                 do i=1,2!min(i_max,2)
                                    if (at%psppar(l,i,ityp) /= 0.0_gp) then 
                                       loop_j2: do j=i+1,3
                                          if (at%psppar(l,j,ityp) == 0.0_gp) exit loop_j2
                                          !offdiagonal HGH term
                                          if (at%npspcode(ityp) == PSPCODE_HGH) then !traditional HGH convention
                                             hij=offdiagarr(i,j-i,l)*at%psppar(l,j,ityp)
                                          else !HGH-K convention
                                             hij=at%psppar(l,i+j+1,ityp)
                                          end if
                                          tt=denskern_gathered(ind)*hij
                                          do m=1,2*l-1
                                             !F_t= 2.0*h_ij (<D_tp_i|psi><psi|p_j>+<p_i|psi><psi|D_tp_j>)
                                             !(the factor two is below)
                                             do icplx=1,ncplx
                                                sp0i=real(scalprod_new(icplx,0,m,i,l,ii),gp)
                                                sp0j=real(scalprod_new(icplx,0,m,j,l,ii),gp)
                                                do idir=1,3
                                                   spi=real(scalprod_new(icplx,idir,m,i,l,jj),gp)
                                                   spj=real(scalprod_new(icplx,idir,m,j,l,jj),gp)
                                                   fxyz_orb(idir,iat)=fxyz_orb(idir,iat)+&
                                                        tt*(sp0j*spi+spj*sp0i)
                                                end do
                                                spi=real(scalprod_new(icplx,0,m,i,l,jj),gp)
                                                spj=real(scalprod_new(icplx,0,m,j,l,jj),gp)
                                                Enl=Enl+denskern_gathered(ind)*(sp0i*spj+sp0j*spi)*hij
                                                do idir=4,ndir
                                                   spi=real(scalprod_new(icplx,idir,m,i,l,jj),gp)
                                                   spj=real(scalprod_new(icplx,idir,m,j,l,jj),gp)
                                                   sab(idir-3)=sab(idir-3)+&   
                                                        2.0_gp*tt*(sp0j*spi+sp0i*spj)!&
                                                end do
                                             end do
                                          end do
                                       end do loop_j2
                                    end if
                                 end do
                              end do
                           end if
                        end do
    
    
                        !seq: strten(1:6) =  11 22 33 23 13 12 
                        !!tmparr(1)=tmparr(1)+1.d0!sab(1)/vol 
                        !!tmparr(2)=tmparr(2)+1.d0!sab(2)/vol 
                        !!tmparr(3)=tmparr(3)+1.d0!sab(3)/vol 
                        !!tmparr(4)=tmparr(4)+1.d0!sab(5)/vol
                        !!tmparr(5)=tmparr(5)+1.d0!sab(6)/vol
                        !!tmparr(6)=tmparr(6)+1.d0!sab(4)/vol
                        if (calculate_strten) then
                           strten_loc(1)=strten_loc(1)+sab(1)/vol 
                           strten_loc(2)=strten_loc(2)+sab(2)/vol 
                           strten_loc(3)=strten_loc(3)+sab(3)/vol 
                           strten_loc(4)=strten_loc(4)+sab(5)/vol
                           strten_loc(5)=strten_loc(5)+sab(6)/vol
                           strten_loc(6)=strten_loc(6)+sab(4)/vol
                        end if
                     end do
                  end do
                  !$omp end do
               end do
            end do spin_loop2
            !$omp end parallel
            do iat=1,ntmp
               iiat=isat_par(iproc)+iat+isat-1
               fsep(1,iiat)=fsep(1,iiat)+2.d0*fxyz_orb(1,iat)
               fsep(2,iiat)=fsep(2,iiat)+2.d0*fxyz_orb(2,iat)
               fsep(3,iiat)=fsep(3,iiat)+2.d0*fxyz_orb(3,iat)
            end do
            if (calculate_strten) strten(:) = strten(:) + strten_loc(:)
            if (ieorb == orbs%norbp) exit loop_kptF2
            ikpt=ikpt+1
            !ispsi_k=ispsi
         end do loop_kptF2
    
      end if natp_if2
    
      !call f_free(sab)
      !call f_free(strten_loc)
      !call f_free(fxyz_orb)
    
      if (calculate_strten) then
         !Adding Enl to the diagonal components of strten after loop over kpts is finished...
         do i=1,3
            strten(i)=strten(i)+Enl/vol
         end do
      end if
    
    
      call f_release_routine()
    
    end subroutine calculate_forces_kernel


    subroutine local_hamiltonian_stress_linear(iproc, nproc, orbs, lzd, hx, hy, hz, npsidim, &
         psi, &!psit_c, psit_f, &
         collcom, msmat, mmat, lsmat, lmat, tens)
      use module_base
      use module_types
      !use module_xc
      use sparsematrix_base, only: sparse_matrix, matrices
      use sparsematrix, only: trace_sparse
      use communications_base, only: TRANSPOSE_FULL
      use communications, only: transpose_localized
      use transposed_operations, only: calculate_overlap_transposed
      use locreg_operations
      implicit none
      integer,intent(in) :: iproc, nproc
      type(orbitals_data), intent(in) :: orbs
      type(local_zone_descriptors), intent(in) :: lzd
      type(comms_linear),intent(in) :: collcom
      type(sparse_matrix),intent(in) :: lsmat
      type(sparse_matrix),intent(inout) :: msmat
      type(matrices),intent(inout) :: mmat, lmat
      real(gp), intent(in) :: hx,hy,hz
      integer,intent(in) :: npsidim
      real(wp), dimension(npsidim), intent(in) :: psi
      !real(kind=8),dimension(collcom%ndimind_c),intent(inout) :: psit_c
      !real(kind=8),dimension(7*collcom%ndimind_f),intent(inout) :: psit_f
      real(gp), intent(inout) :: tens(6)
      !local variables
      real(gp) :: ekin_sum,epot_sum
      character(len=*), parameter :: subname='local_hamiltonian_stress'
      integer :: iorb, npot, i_f, iseg_f, iiorb, ilr, ist, idir, iidim
      !real(wp) :: kinstr(6)
      real(gp) :: ekin,kx,ky,kz,etest, tt, ekino
      type(workarr_locham) :: w
      real(wp), dimension(:,:), allocatable :: psir
      real(kind=8),dimension(0:3) :: scal=1.d0
      real(kind=8),dimension(3) :: hgridh
      real(kind=8),dimension(:,:),allocatable :: hpsit_c, hpsit_f, hpsi
      real(kind=8),dimension(:),allocatable :: psit_c, psit_f
    
      call f_routine(id='local_hamiltonian_stress')
    
      !@ NEW ####################################################
    
      hgridh(1)=hx*.5_gp
      hgridh(2)=hy*.5_gp
      hgridh(3)=hz*.5_gp
    
      hpsit_c = f_malloc((/collcom%ndimind_c,3/),id='hpsit_c')
      hpsit_f = f_malloc((/7*collcom%ndimind_f,3/),id='hpsit_f')
      hpsi = f_malloc((/npsidim,3/),id='hpsi')
      psit_c = f_malloc(collcom%ndimind_c,id='psit_c')
      psit_f = f_malloc(7*collcom%ndimind_f,id='psit_f')
    
      call initialize_work_arrays_locham(lzd%nlr, lzd%llr, orbs%nspinor, .true., w)
      ist = 1
      do iorb=1,orbs%norbp
         iiorb = orbs%isorb + iorb
         ilr = orbs%inwhichlocreg(iiorb)
    
         i_f=min(1,lzd%llr(ilr)%wfd%nvctr_f)
         iseg_f=min(1,lzd%llr(ilr)%wfd%nseg_f)
    
    
         !!psir = f_malloc0((/lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i,orbs%nspinor/),id='psir')
         call initialize_work_arrays_locham(1, [Lzd%Llr(ilr)], orbs%nspinor, .false., w)
    
         !!call daub_to_isf_locham(orbs%nspinor, lzd%llr(ilr), wrk_lh, psi(ist), psir)
         call uncompress_forstandard(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
              lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
              lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, &
              lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, &
              lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, &
              lzd%llr(ilr)%wfd%keygloc(1,1), lzd%llr(ilr)%wfd%keyvloc(1), &
              lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
              lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+iseg_f), &
              lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+iseg_f), &
              scal, psi(ist), psi(ist+lzd%llr(ilr)%wfd%nvctr_c+i_f-1), &
              w%x_c(1,1), w%x_f(1,1), w%x_f1(1,1), w%x_f2(1,1), w%x_f3(1,1))
    
         do idir=1,3
            iidim = 10**(3-idir)
            call f_zero(w%y_c)
            call f_zero(w%y_f)
            !!call convolut_kinetic_per_T_1D(2*lzd%llr(ilr)%d%n1+1, 2*lzd%llr(ilr)%d%n2+1, 2*lzd%llr(ilr)%d%n3+1, &
            !!     hgridh(idir), idir, w%x_c(1,1), w%y_c(1,1))
            call ConvolkineticT(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
                 lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
                 lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, &
                 lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, &
                 hx, hy, hz, &      !here the grid spacings are supposed to be equal.  SM: not any more
                 lzd%llr(ilr)%bounds%kb%ibyz_c, lzd%llr(ilr)%bounds%kb%ibxz_c, lzd%llr(ilr)%bounds%kb%ibxy_c, &
                 lzd%llr(ilr)%bounds%kb%ibyz_f, lzd%llr(ilr)%bounds%kb%ibxz_f, lzd%llr(ilr)%bounds%kb%ibxy_f, &
                 w%x_c(1,1), w%x_f(1,1), &
                 w%y_c(1,1), w%y_f(1,1), ekino, &
                 w%x_f1(1,1), w%x_f2(1,1), w%x_f3(1,1),iidim)
            call compress_forstandard(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
                 lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
                 lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, &
                 lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, &
                 lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, &
                 lzd%llr(ilr)%wfd%keygloc(1,1), lzd%llr(ilr)%wfd%keyvloc(1), &
                 lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
                 lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+iseg_f), &
                 lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+iseg_f), &
                 scal, w%y_c(1,1), w%y_f(1,1), hpsi(ist,idir), hpsi(ist+lzd%llr(ilr)%wfd%nvctr_c+i_f-1,idir))
         end do
    
         ist = ist + lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f
      end do
    
      call transpose_localized(iproc, nproc, npsidim, orbs, collcom, TRANSPOSE_FULL, &
           psi, psit_c, psit_f, lzd)
      do idir=1,3
         call transpose_localized(iproc, nproc, npsidim, orbs, collcom, TRANSPOSE_FULL, &
              hpsi(1,idir), hpsit_c, hpsit_f, lzd)
         call calculate_overlap_transposed(iproc, nproc, orbs, collcom, &
              psit_c, hpsit_c, psit_f, hpsit_f, msmat, mmat)
         tt = trace_sparse(iproc, nproc, orbs, msmat, lsmat, mmat%matrix_compr, lmat%matrix_compr, 1)
         !tens(idir) = tens(idir) + -8.0_gp/(hx*hy*hz)/real(lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i,gp)*tt
         tens(idir) = tens(idir) - 2.0_gp*8.0_gp/(hx*hy*hz)/real(lzd%glr%d%n1i*lzd%glr%d%n2i*lzd%glr%d%n3i,gp)*tt
         tens(idir) = tens(idir)/real(nproc,kind=8) !divide by nproc since an allreduce will follow
      end do
    
      call deallocate_work_arrays_locham(w)
      call f_free(hpsit_c)
      call f_free(hpsit_f)
      call f_free(psit_c)
      call f_free(psit_f)
      call f_free(hpsi)
    
      !@ END NEW ################################################
    
    !!!initialise the work arrays
      !!call initialize_work_arrays_locham(1,lr,orbs%nspinor,.true.,wrk_lh)  
    
      !!tens=0.d0
    
    !!!components of the potential
      !!npot=orbs%nspinor
      !!if (orbs%nspinor == 2) npot=1
    
      !!hpsi = f_malloc((/ lr%wfd%nvctr_c+7*lr%wfd%nvctr_f , orbs%nspinor*orbs%norbp /),id='hpsi')
      !!hpsi=0.0_wp
    !!! Wavefunction in real space
      !!psir = f_malloc0((/ lr%d%n1i*lr%d%n2i*lr%d%n3i, orbs%nspinor /),id='psir')
    !!!call to_zero(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%nspinor,psir)
    
    
    
      !!ekin_sum=0.0_gp
      !!epot_sum=0.0_gp
    
      !!etest=0.0_gp
    
    
      !!do iorb=1,orbs%norbp
      !!   kinstr=0._wp
      !!   oidx=(iorb-1)*orbs%nspinor+1
    
      !!   call daub_to_isf_locham(orbs%nspinor,lr,wrk_lh,psi(1,oidx),psir)
    
      !!   kx=orbs%kpts(1,orbs%iokpt(iorb))
      !!   ky=orbs%kpts(2,orbs%iokpt(iorb))
      !!   kz=orbs%kpts(3,orbs%iokpt(iorb))
    
      !!   call isf_to_daub_kinetic(hx,hy,hz,kx,ky,kz,orbs%nspinor,lr,wrk_lh,&
      !!        psir,hpsi(1,oidx),ekin,kinstr)
    
      !!   kinstr = -kinstr*8.0_gp/(hx*hy*hz)/real(lr%d%n1i*lr%d%n2i*lr%d%n3i,gp)
      !!   tens=tens+kinstr*2.0_gp*&
      !!   orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)
    
      !!end do !loop over orbitals: finished
    
    !!!deallocations of work arrays
      !!call f_free(psir)
    
      !!call f_free(hpsi)
    
      !!call deallocate_work_arrays_locham(wrk_lh)
    
    
      call f_release_routine()
    
    END SUBROUTINE local_hamiltonian_stress_linear

end module forces_linear
