!> @file
!!   Routines to precondition wavefunctions
!! @author
!!    Copyright (C) 2005-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
 

!> Calls the preconditioner for each orbital treated by the processor
subroutine preconditionall(orbs,lr,hx,hy,hz,ncong,hpsi,gnrm,gnrm_zero)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: ncong
  real(gp), intent(in) :: hx,hy,hz
  type(locreg_descriptors), intent(in) :: lr
  type(orbitals_data), intent(in) :: orbs
  real(dp), intent(out) :: gnrm,gnrm_zero
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(inout) :: hpsi
  !local variables
  integer :: iorb,inds,ncplx,ikpt,jorb
  real(wp) :: cprecr,scpr,evalmax,eval_zero
  real(gp) :: kx,ky,kz

  ! Preconditions all orbitals belonging to iproc
  !and calculate the norm of the residue

  ! norm of gradient
  gnrm=0.0_dp
  !norm of gradient of unoccupied orbitals
  gnrm_zero=0.0_dp


  !commented out, never used
!   evalmax=orbs%eval(orbs%isorb+1)
!   do iorb=1,orbs%norbp
!     evalmax=max(orbs%eval(orbs%isorb+iorb),evalmax)
!   enddo
!   call MPI_ALLREDUCE(evalmax,eval_zero,1,mpidtypd,&
!        MPI_MAX,bigdft_mpi%mpi_comm,ierr)


  if (orbs%norbp >0) ikpt=orbs%iokpt(1)
  do iorb=1,orbs%norbp
     !if it is the first orbital or the k-point has changed calculate the max
     if (orbs%iokpt(iorb) /= ikpt .or. iorb == 1) then
        !the eval array contains all the values
        !take the max for all k-points
        !one may think to take the max per k-point
        evalmax=orbs%eval((orbs%iokpt(iorb)-1)*orbs%norb+1)
        do jorb=1,orbs%norb
           evalmax=max(orbs%eval((orbs%iokpt(iorb)-1)*orbs%norb+jorb),evalmax)
        enddo
        eval_zero=evalmax
        ikpt=orbs%iokpt(iorb)
     end if

     !indo=(iorb-1)*nspinor+1
     !loop over the spinorial components
     !k-point values, if present
     kx=orbs%kpts(1,orbs%iokpt(iorb))
     ky=orbs%kpts(2,orbs%iokpt(iorb))
     kz=orbs%kpts(3,orbs%iokpt(iorb))
!       print *, iorb, orbs%kpts(1,orbs%iokpt(iorb)), orbs%kpts(2,orbs%iokpt(iorb)), orbs%kpts(3,orbs%iokpt(iorb))

     !real k-point different from Gamma still not implemented
     if (kx**2+ky**2+kz**2 > 0.0_gp .or. orbs%nspinor==2 ) then
        ncplx=2
     else
        ncplx=1
     end if

     do inds=1,orbs%nspinor,ncplx

        !the nrm2 function can be replaced here by ddot
        scpr=nrm2(ncplx*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),hpsi(1,inds,iorb),1)
        if (orbs%occup(orbs%isorb+iorb) == 0.0_gp) then
           gnrm_zero=gnrm_zero+orbs%kwgts(orbs%iokpt(iorb))*scpr**2
        else
           !write(17,*)'iorb,gnrm',orbs%isorb+iorb,scpr**2
           gnrm=gnrm+orbs%kwgts(orbs%iokpt(iorb))*scpr**2
        end if
        !write(*,*) 'preconditionall: verbosity ',verbose
!           write(*,*)'iorb,gnrm',orbs%isorb+iorb,scpr**2

       if (scpr /= 0.0_wp) then

          call cprecr_from_eval(lr%geocode,eval_zero,orbs%eval(orbs%isorb+iorb),cprecr)          
           !cases with no CG iterations, diagonal preconditioning
           !for Free BC it is incorporated in the standard procedure
           if (ncong == 0 .and. lr%geocode /= 'F') then
              select case(lr%geocode)
              case('F')
              case('S')
                 call prec_fft_slab(lr%d%n1,lr%d%n2,lr%d%n3, &
                      lr%wfd%nseg_c,lr%wfd%nvctr_c,lr%wfd%nseg_f,&
                      lr%wfd%nvctr_f,lr%wfd%keygloc,lr%wfd%keyvloc, &
                      cprecr,hx,hy,hz,hpsi(1,inds,iorb))
              case('P')
                 call prec_fft(lr%d%n1,lr%d%n2,lr%d%n3, &
                      lr%wfd%nseg_c,lr%wfd%nvctr_c,lr%wfd%nseg_f,lr%wfd%nvctr_f,&
                      lr%wfd%keygloc,lr%wfd%keyvloc, &
                      cprecr,hx,hy,hz,hpsi(1,inds,iorb))
              end select

           else !normal preconditioner
              
              call precondition_residue(lr,ncplx,ncong,cprecr,&
                   hx,hy,hz,kx,ky,kz,hpsi(1,inds,iorb))

           end if

        end if

!     print *,iorb,inds,dot(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f, hpsi(1,inds,iorb),1,hpsi(1,inds,iorb),1)
!     print *,iorb,inds+1,dot(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f, hpsi(1,inds+1,iorb),1,hpsi(1,inds+1,iorb),1)
     end do
  enddo

END SUBROUTINE preconditionall


!> Generalized for the Linearscaling code
subroutine preconditionall2(iproc,nproc,orbs,Lzd,hx,hy,hz,ncong,npsidim,hpsi,confdatarr,gnrm,gnrm_zero,&
                            linear_precond_convol_workarrays, linear_precond_workarrays)
  use module_base
  use module_types
  use module_interfaces, except_this_one => preconditionall2
  use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
  use yaml_output
  implicit none
  integer, intent(in) :: iproc,nproc,ncong,npsidim
  real(gp), intent(in) :: hx,hy,hz
  type(local_zone_descriptors), intent(in) :: Lzd
  type(orbitals_data), intent(in) :: orbs
  real(dp), intent(out) :: gnrm,gnrm_zero
  real(wp), dimension(npsidim), intent(inout) :: hpsi
  type(confpot_data), dimension(orbs%norbp), intent(in) :: confdatarr !< used in the linear scaling but also for the cubic case
  type(workarrays_quartic_convolutions),dimension(orbs%norbp),intent(inout),optional :: linear_precond_convol_workarrays !< convolution workarrays for the linear case
  type(workarr_precond),dimension(orbs%norbp),intent(inout),optional :: linear_precond_workarrays !< workarrays for the linear case
  !local variables
  character(len=*), parameter :: subname='preconditionall2'
  integer :: iorb,inds,ncplx,ikpt,jorb,ist,ilr,ierr,jproc
  real(wp) :: cprecr,scpr,evalmax,eval_zero,gnrm_orb
  real(gp) :: kx,ky,kz
!!$  integer :: i_stat,i_all,ispinor,nbox
!!$  logical, parameter :: newp=.true.
!!$  real(gp) :: eh_fake,monop
!!$  real(wp), dimension(:), allocatable :: hpsir
!!$  type(coulomb_operator) :: G_Helmholtz
!!$  real(wp), dimension(:,:), allocatable :: gnrm_per_orb
!!$  type(workarr_sumrho) :: w
  integer, dimension(:,:), allocatable :: ncntdsp
  real(wp), dimension(:), allocatable :: gnrms,gnrmp
  !debug
!!$  type(atoms_data) atoms_fake
!!$  integer :: iter=0
!!$  iter=iter+1

  call f_routine(id=subname)

  ! Preconditions all orbitals belonging to iproc
  !and calculate the norm of the residue
  ! norm of gradient
  gnrm=0.0_dp
  !norm of gradient of unoccupied orbitals
  gnrm_zero=0.0_dp

  !prepare the arrays for the 
  if (verbose >= 3) then
     gnrmp = f_malloc(max(orbs%norbp, 1),id='gnrmp')
  end if

!!$  if (newp) then
!!$     ilr=1
!!$     hpsir=f_malloc(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,id='hpsir',routine_id=subname)
!!$     call initialize_work_arrays_sumrho(1,Lzd%Llr(ilr),.true.,w)
!!$  end if
  !if (iproc.eq. 0 .and. verbose.ge.3) write(*,*) ' '
  ist = 0
  if (orbs%norbp >0) ikpt=orbs%iokpt(1)
  do iorb=1,orbs%norbp
     ilr = orbs%inwhichlocreg(iorb+orbs%isorb)
     !if it is the first orbital or the k-point has changed calculate the max
     if (orbs%iokpt(iorb) /= ikpt .or. iorb == 1) then
        !the eval array contains all the values
        !take the max for all k-points
        !one may think to take the max per k-point
        evalmax=orbs%eval((orbs%iokpt(iorb)-1)*orbs%norb+1)
        do jorb=1,orbs%norb
           evalmax=max(orbs%eval((orbs%iokpt(iorb)-1)*orbs%norb+jorb),evalmax)
        enddo
        eval_zero=evalmax
        ikpt=orbs%iokpt(iorb)
     end if
     !print *,'iorb,eval,evalmax',iorb+orbs%isorb,orbs%eval(iorb+orbs%isorb),eval_zero
     !indo=(iorb-1)*nspinor+1
     !loop over the spinorial components
     !k-point values, if present
     kx=orbs%kpts(1,orbs%iokpt(iorb))
     ky=orbs%kpts(2,orbs%iokpt(iorb))
     kz=orbs%kpts(3,orbs%iokpt(iorb))
!       print *, iorb, orbs%kpts(1,orbs%iokpt(iorb)), orbs%kpts(2,orbs%iokpt(iorb)), orbs%kpts(3,orbs%iokpt(iorb))

     !real k-point different from Gamma still not implemented
     if (kx**2+ky**2+kz**2 > 0.0_gp .or. orbs%nspinor==2) then
        ncplx=2
     else
        ncplx=1
     end if

     gnrm_orb=0.0_wp
     do inds=1,orbs%nspinor,ncplx

        !the nrm2 function can be replaced here by ddot
        scpr=nrm2(ncplx*(Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f),hpsi(1+ist),1)
        if (orbs%occup(orbs%isorb+iorb) == 0.0_gp) then
           gnrm_zero=gnrm_zero+orbs%kwgts(orbs%iokpt(iorb))*scpr**2
        else
           !write(*,*)'iorb,gnrm',orbs%isorb+iorb,scpr**2,ilr
           gnrm=gnrm+orbs%kwgts(orbs%iokpt(iorb))*scpr**2
        end if
        if (verbose >= 3) then
           gnrm_orb=gnrm_orb+scpr
           if (inds+ncplx-1==orbs%nspinor) gnrmp(iorb)=gnrm_orb
           !write(*,*) 'iorb,gnrm,ilr',orbs%isorb+iorb,scpr,ilr,gnrm_orb,iproc
        end if

!!$        print *,'plotting gradient iorb,initial',iorb,&
!!$             nrm2(Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f,hpsi(1+ist),1)
!!$        atoms_fake=atoms_null()
!!$        atoms_fake%nat=0
!!$        atoms_fake%geocode='F'
!!$        atoms_fake%alat1=Lzd%hgrids(1)*(Lzd%Llr(ilr)%d%n1+1)
!!$        atoms_fake%alat2=Lzd%hgrids(2)*(Lzd%Llr(ilr)%d%n2+1)
!!$        atoms_fake%alat3=Lzd%hgrids(3)*(Lzd%Llr(ilr)%d%n3+1)
!!$        call plot_wf('G'//trim(adjustl(yaml_toa(iorb)))//'-'//trim(adjustl(yaml_toa(iter))),1,&
!!$             atoms_fake,1.0_gp,Lzd%llr(ilr),&
!!$             Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),(/0.0_gp,0.0_gp,0.0_gp/),hpsi(1+ist))
!!$
!!$        if (newp) then
!!$           call cprecr_from_eval(Lzd%Llr(ilr)%geocode,eval_zero,orbs%eval(orbs%isorb+iorb),cprecr)
!!$           !alternative preconditioner
!!$           call vscal(ncplx*(Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f),0.5_gp/pi_param,hpsi(1+ist),1)
!!$           call daub_to_isf(Lzd%Llr(ilr),w,hpsi(1+ist),hpsir(1))
!!$           !sequential kernel
!!$           G_Helmholtz=pkernel_init(.false.,0,1,0,Lzd%Llr(ilr)%geocode,&
!!$                (/Lzd%Llr(ilr)%d%n1i,Lzd%Llr(ilr)%d%n2i,Lzd%Llr(ilr)%d%n3i/),&
!!$                0.5_gp*Lzd%hgrids,16,mu0_screening=sqrt(2.0_gp*abs(cprecr)))
!!$           
!!$           call pkernel_set(G_Helmholtz,.true.)
!!$           
!!$           !apply it to the gradient to smooth it
!!$           call H_potential('D',G_Helmholtz,hpsir(1),hpsir(1),eh_fake,0.d0,.false.)
!!$           call pkernel_free(G_Helmholtz,subname)
!!$           !convert the gradient back to the locreg
!!$           call isf_to_daub(Lzd%Llr(ilr),w,hpsir(1),hpsi(1+ist))
!!$           
!!$        end if

        
       if (scpr /= 0.0_wp) then
          call cprecr_from_eval(Lzd%Llr(ilr)%geocode,eval_zero,orbs%eval(orbs%isorb+iorb),cprecr)
           !cases with no CG iterations, diagonal preconditioning
           !for Free BC it is incorporated in the standard procedure
           if (ncong == 0 .and. Lzd%Llr(ilr)%geocode /= 'F') then
              select case(Lzd%Llr(ilr)%geocode)
              case('F')
              case('S')
                 call prec_fft_slab(Lzd%Llr(ilr)%d%n1,Lzd%Llr(ilr)%d%n2,Lzd%Llr(ilr)%d%n3, &
                      Lzd%Llr(ilr)%wfd%nseg_c,Lzd%Llr(ilr)%wfd%nvctr_c,Lzd%Llr(ilr)%wfd%nseg_f,&
                      Lzd%Llr(ilr)%wfd%nvctr_f,Lzd%Llr(ilr)%wfd%keygloc,Lzd%Llr(ilr)%wfd%keyvloc, &
                      cprecr,hx,hy,hz,hpsi(1+ist))
              case('P')
                 call prec_fft(Lzd%Llr(ilr)%d%n1,Lzd%Llr(ilr)%d%n2,Lzd%Llr(ilr)%d%n3, &
                      Lzd%Llr(ilr)%wfd%nseg_c,Lzd%Llr(ilr)%wfd%nvctr_c,&
                      Lzd%Llr(ilr)%wfd%nseg_f,Lzd%Llr(ilr)%wfd%nvctr_f,&
                      Lzd%Llr(ilr)%wfd%keygloc,Lzd%Llr(ilr)%wfd%keyvloc, &
                      cprecr,hx,hy,hz,hpsi(1+ist))
              end select

           else !normal preconditioner
              !case active only in the linear scaling case
              if(confdatarr(iorb)%prefac > 0.0_gp .or. confdatarr(iorb)%potorder > 0)then
              !   call yaml_map('Localizing preconditioner factor',confdatarr(iorb)%prefac)
              !!write(1000+orbs%isorb+iorb,'(a,2i4,3x,2i8)') 'id, ilr, centerx, startx', orbs%isorb+iorb, ilr, nint(Lzd%Llr(ilr)%locregCenter(1)/hx), lzd%llr(ilr)%ns1
              !!write(1000+orbs%isorb+iorb,'(a,2i4,3x,2i8)') 'id, ilr, centery, starty', orbs%isorb+iorb, ilr, nint(Lzd%Llr(ilr)%locregCenter(2)/hy), lzd%llr(ilr)%ns2
              !!write(1000+orbs%isorb+iorb,'(a,2i4,3x,2i8)') 'id, ilr, centerz, startz', orbs%isorb+iorb, ilr, nint(Lzd%Llr(ilr)%locregCenter(3)/hz), lzd%llr(ilr)%ns3

              if (.not.present(linear_precond_convol_workarrays)) then
                  call f_err_throw("linear_precond_convol_workarrays must be present when calling the linear preconditioner", &
                                   err_name='BIGDFT_RUNTIME_ERROR')
              end if
              if (.not.present(linear_precond_workarrays)) then
                  call f_err_throw("linear_precond_workarrays must be present when calling the linear preconditioner", &
                                   err_name='BIGDFT_RUNTIME_ERROR')
              end if
                 call solvePrecondEquation(iproc,nproc,Lzd%Llr(ilr),ncplx,ncong,&
                      cprecr,&
                      hx,hy,hz,kx,ky,kz,hpsi(1+ist),&
                      Lzd%Llr(ilr)%locregCenter, orbs,&
                      confdatarr(iorb)%prefac,&
                      confdatarr(iorb)%potorder,&
                      linear_precond_convol_workarrays(iorb),linear_precond_workarrays(iorb))

!                 call solvePrecondEquation(Lzd%Llr(ilr),ncplx,ncong,cprecr,&
!                   hx,hy,hz,kx,ky,kz,hpsi(1+ist), rxyz(1,ilr), orbs,&                         !here should change rxyz to be center of Locreg
!                   potentialPrefac(ilr), confPotOrder, 1)                         ! should depend on locreg not atom type? 'it' is commented in lower routines, so put 1
              else
                 call precondition_residue(Lzd%Llr(ilr),ncplx,ncong,cprecr,&
                      hx,hy,hz,kx,ky,kz,hpsi(1+ist))
              end if
           end if

        end if

!!$        print *,'iorb,newgradient',iorb,scpr,nrm2(ncplx*(Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f),hpsi(1+ist),1),eh_fake
!!$        call plot_wf('GPnew'//trim(adjustl(yaml_toa(iorb)))//'-'//trim(adjustl(yaml_toa(iter))),1,&
!!$             atoms_fake,1.0_gp,Lzd%llr(ilr),&
!!$             Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),(/0.0_gp,0.0_gp,0.0_gp/),hpsi(1+ist))


       ist = ist + (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*ncplx
!     print *,iorb,inds,dot(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f, hpsi(1,inds,iorb),1,hpsi(1,inds,iorb),1)
!     print *,iorb,inds+1,dot(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f, hpsi(1,inds+1,iorb),1,hpsi(1,inds+1,iorb),1)
    end do

 enddo

!!$ if (newp) then
!!$    call deallocate_work_arrays_sumrho(w)
!!$    call f_free(hpsir)
!!$ end if
  !gather the results of the gnrm per orbital in the case of high verbosity
  if (verbose >= 3) then
     gnrms = f_malloc(orbs%norb*orbs%nkpts,id='gnrms')
     !prepare displacements arrays
     ncntdsp = f_malloc((/ nproc, 2 /),id='ncntdsp')
     ncntdsp(1,2)=0
     ncntdsp(1,1)=orbs%norb_par(0,0)
     do jproc=1,nproc-1
        ncntdsp(jproc+1,2)=ncntdsp(jproc,2)+ncntdsp(jproc,1)
        ncntdsp(jproc+1,1)=orbs%norb_par(jproc,0)
     end do
     call to_zero(orbs%norb*orbs%nkpts,gnrms(1))
     !root mpi task collects the data
     if (nproc > 1) then
        call MPI_GATHERV(gnrmp(1),orbs%norbp,mpidtypw,gnrms(1),ncntdsp(1,1),&
             ncntdsp(1,2),mpidtypw,0,bigdft_mpi%mpi_comm,ierr)
     else
        call vcopy(orbs%norb*orbs%nkpts,gnrmp(1),1,gnrms(1),1)
     end if

     !if (iproc ==0) print *,'ciao',gnrmp,orbs%nspinor

     !write the values per orbitals
     if (iproc ==0) call write_gnrms(orbs%nkpts,orbs%norb,gnrms)


     call f_free(ncntdsp)
     call f_free(gnrms)
     call f_free(gnrmp)
  end if

  call f_release_routine()

END SUBROUTINE preconditionall2


! > This function has been created also for the GPU-ported routines
subroutine cprecr_from_eval(geocode,eval_zero,eval,cprecr)
  use module_base
  implicit none
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  real(gp), intent(in) :: eval,eval_zero
  real(gp), intent(out) :: cprecr

  select case(geocode)
  case('F')
     cprecr=sqrt(.2d0**2+min(0.d0,eval)**2)
  case('S')
     cprecr=sqrt(0.2d0**2+(eval-eval_zero)**2)
  case('P')
     cprecr=sqrt(0.2d0**2+(eval-eval_zero)**2)
  end select

END SUBROUTINE cprecr_from_eval


!> Routine used for the k-points, eventually to be used for all cases
subroutine precondition_residue(lr,ncplx,ncong,cprecr,&
     hx,hy,hz,kx,ky,kz,x)
  use module_base
  use module_types
  ! Solves (KE+cprecr*I)*xx=yy by conjugate gradient method
  ! x is the right hand side on input and the solution on output
  implicit none
  integer, intent(in) :: ncong,ncplx
  real(gp), intent(in) :: hx,hy,hz,cprecr,kx,ky,kz
  type(locreg_descriptors), intent(in) :: lr
  real(wp), dimension((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*ncplx), intent(inout) :: x
  ! local variables
  character(len=*), parameter :: subname='precondition_residue'
  real(gp), dimension(0:7) :: scal
  real(wp) :: rmr_old,rmr_new,alpha,beta
  integer :: icong
  type(workarr_precond) :: w
  real(wp), dimension(:), allocatable :: b,r,d

  !arrays for the CG procedure
  b = f_malloc(ncplx*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),id='b')
  r = f_malloc(ncplx*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),id='r')
  d = f_malloc(ncplx*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),id='d')

  call allocate_work_arrays(lr%geocode,lr%hybrid_on,ncplx,lr%d,w)

  call precondition_preconditioner(lr,ncplx,hx,hy,hz,scal,cprecr,w,x,b)

  call precond_locham(ncplx,lr,hx,hy,hz,kx,ky,kz,cprecr,x,d,w,scal)

  rmr_new=dot(ncplx*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),d(1),1,d(1),1)
  !write(*,*)'debug1',rmr_new

  !this operation should be rewritten in a better way
  r=b-d ! r=b-Ax

  call calculate_rmr_new(lr%geocode,lr%hybrid_on,ncplx,lr%wfd,scal,r,d,rmr_new)
  !stands for
  !d=r
  !rmr_new=dot_product(r,r)


  do icong=1,ncong 
!     write(*,*)'hello',icong,rmr_new

     call precond_locham(ncplx,lr,hx,hy,hz,kx,ky,kz,cprecr,d,b,w,scal)! b:=Ad

     !in the complex case these objects are to be supposed real
     alpha=rmr_new/dot(ncplx*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),d(1),1,b(1),1)

     call axpy(ncplx*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),alpha,d(1),1,x(1),1)
     call axpy(ncplx*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),-alpha,b(1),1,r(1),1)

     if (icong==ncong) exit

     rmr_old=rmr_new

     call calculate_rmr_new(lr%geocode,lr%hybrid_on,ncplx,lr%wfd,scal,r,b,rmr_new)

     beta=rmr_new/rmr_old
!print *,'beta.icong',icong,beta
     d=b+beta*d
    
  enddo

  call finalise_precond_residue(lr%geocode,lr%hybrid_on,ncplx,lr%wfd,scal,x)

  !write(*,*)'debug2',dot(ncplx*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),x(1),1,x(1),1)


  call f_free(b)
  call f_free(r)
  call f_free(d)

  call deallocate_work_arrays(lr%geocode,lr%hybrid_on,ncplx,w)

END SUBROUTINE precondition_residue


subroutine finalise_precond_residue(geocode,hybrid_on,ncplx,wfd,scal,x)
  use module_base
  use module_types
  implicit none
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  logical, intent(in) :: hybrid_on
  integer, intent(in) :: ncplx
  type(wavefunctions_descriptors), intent(in) :: wfd
  real(gp), dimension(0:7), intent(in) :: scal
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,ncplx), intent(inout) :: x
  !local variables
  integer :: idx

  call f_routine(id='finalise_precond_residue')

  if (geocode == 'F') then
     do idx=1,ncplx
        call wscalv_wrap(wfd%nvctr_c,wfd%nvctr_f,scal,x(1,idx))
     end do
  else if ((geocode == 'P' .and. .not. hybrid_on) .or. geocode == 'S') then
     do idx=1,ncplx
        ! x=D^{-1/2}x'
        call wscal_per_self(wfd%nvctr_c,wfd%nvctr_f,scal,x(1,idx),&
             x(wfd%nvctr_c+min(1,wfd%nvctr_f),idx))
        !	write(30,*) x
        !	stop
     end do
  else
  end if

  call f_release_routine()

END SUBROUTINE finalise_precond_residue


subroutine calculate_rmr_new(geocode,hybrid_on,ncplx,wfd,scal,r,b,rmr_new)
  use module_base
  use module_types
  implicit none
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  logical, intent(in) :: hybrid_on
  integer, intent(in) :: ncplx
  type(wavefunctions_descriptors), intent(in) :: wfd
  real(gp), dimension(0:7), intent(in) :: scal
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,ncplx), intent(in) :: r
  real(wp), intent(out) :: rmr_new
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,ncplx), intent(out) :: b
  !local variables
  logical :: noscal
  integer :: idx

  call f_routine(id='calculate_rmr_new')

  noscal = ((geocode == 'P' .and. .not. hybrid_on) .or. &
       geocode == 'F' .or. geocode == 'S')

  if (noscal) then
     call vcopy(ncplx*(wfd%nvctr_c+7*wfd%nvctr_f),r(1,1),1,b(1,1),1) 
     rmr_new=dot(ncplx*(wfd%nvctr_c+7*wfd%nvctr_f),r(1,1),1,r(1,1),1)
  else 
     do idx=1,ncplx
        call wscal_per(wfd%nvctr_c,wfd%nvctr_f,scal,r(1,idx),&
             r(wfd%nvctr_c+min(1,wfd%nvctr_f),idx),&
             b(1,idx),b(wfd%nvctr_c+min(1,wfd%nvctr_f),idx))
     end do
     rmr_new=dot(ncplx*(wfd%nvctr_c+7*wfd%nvctr_f),r(1,1),1,b(1,1),1)
  end if

  call f_release_routine()

END SUBROUTINE calculate_rmr_new


subroutine precondition_preconditioner(lr,ncplx,hx,hy,hz,scal,cprecr,w,x,b)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: ncplx
  real(gp), intent(in) :: hx,hy,hz,cprecr
  type(locreg_descriptors), intent(in) :: lr
  type(workarr_precond), intent(inout) :: w
  real(gp), dimension(0:7), intent(out) :: scal
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,ncplx), intent(inout) ::  x
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,ncplx), intent(out) ::  b
  !local variables
  logical, parameter :: inguess_on=.true.
  !       wavelet and scaling function second derivative filters
  real(wp), parameter :: b2=24.8758460293923314_wp, a2=3.55369228991319019_wp
  integer :: nd1,nd2,nd3,idx
  integer :: n1f,n3f,n1b,n3b,nd1f,nd3f,nd1b,nd3b 
  real(gp) :: fac
  real(wp) :: fac_h,h0,h1,h2,h3

  call f_routine(id='precondition_preconditioner')
    
  if (lr%geocode == 'F') then
     !using hx instead of hgrid for isolated bc
     fac_h=1.0_wp/real(hx,wp)**2
     h0=    1.5_wp*a2*fac_h
     h1=(a2+b2*.5_wp)*fac_h
     h2=(a2*.5_wp+b2)*fac_h
     h3=    1.5_wp*b2*fac_h

     scal(0)=sqrt(1.0_wp/(h0+cprecr)) 
     scal(1)=sqrt(1.0_wp/(h1+cprecr)) 
     scal(2)=sqrt(1.0_wp/(h2+cprecr)) 
     scal(3)=sqrt(1.0_wp/(h3+cprecr))

     do idx=1,ncplx
        if (inguess_on) then
           !the right hand side is temporarily stored in the rpsi array
           !rpsi=hpsi           
           call vcopy(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,x(1,idx),1,b(1,idx),1) 
           !          and preconditioned with d^{-1/2} as usual:
           call wscalv_wrap(lr%wfd%nvctr_c,lr%wfd%nvctr_f,scal,b(1,idx))
           !hpsi is now diagonally preconditioned with alexey's old preconditioner;
           !inside the diagonal preconditioner a factor of d^{1/2} was added
           !to make the overall factor d^{-1/2} again

           call prec_diag(lr%d%n1,lr%d%n2,lr%d%n3,hx,lr%wfd%nseg_c,&
                lr%wfd%nvctr_c,lr%wfd%nvctr_f,&
                lr%wfd%keygloc,lr%wfd%keyvloc,&
                x(1,idx),x(lr%wfd%nvctr_c+min(1,lr%wfd%nvctr_f),idx),cprecr,scal,a2,b2)

        else
           !assume as input guess x=y
           !hpsi is preconditioned with d^{-1/2} as usual
           call wscalv_wrap(lr%wfd%nvctr_c,lr%wfd%nvctr_f,scal,x(1,idx))

           !b=x
           call vcopy(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,x(1,idx),1,b(1,idx),1) 
        endif
     end do

     !initalize to zero the work arrays, probably not needed
     call to_zero((lr%d%nfu1-lr%d%nfl1+1)*(lr%d%nfu2-lr%d%nfl2+1)*(lr%d%nfu3-lr%d%nfl3+1),&
          w%x_f1(1))
     call to_zero((lr%d%nfu1-lr%d%nfl1+1)*(lr%d%nfu2-lr%d%nfl2+1)*(lr%d%nfu3-lr%d%nfl3+1),&
          w%x_f2(1))
     call to_zero((lr%d%nfu1-lr%d%nfl1+1)*(lr%d%nfu2-lr%d%nfl2+1)*(lr%d%nfu3-lr%d%nfl3+1),&
          w%x_f3(1))
     call to_zero((lr%d%n1+1)*(lr%d%n2+1)*(lr%d%n3+1),w%xpsig_c(0,0,0))
     call to_zero(7*(lr%d%nfu1-lr%d%nfl1+1)*(lr%d%nfu2-lr%d%nfl2+1)*(lr%d%nfu3-lr%d%nfl3+1),&
          w%xpsig_f(1,lr%d%nfl1,lr%d%nfl2,lr%d%nfl3))

     call to_zero((lr%d%n1+1)*(lr%d%n2+1)*(lr%d%n3+1),w%ypsig_c(0,0,0))
     call to_zero(7*(lr%d%nfu1-lr%d%nfl1+1)*(lr%d%nfu2-lr%d%nfl2+1)*(lr%d%nfu3-lr%d%nfl3+1),&
          w%ypsig_f(1,lr%d%nfl1,lr%d%nfl2,lr%d%nfl3))

  else if (lr%geocode == 'P') then

     call dimensions_fft(lr%d%n1,lr%d%n2,lr%d%n3,&
          nd1,nd2,nd3,n1f,n3f,n1b,n3b,nd1f,nd3f,nd1b,nd3b)

     if (ncplx /=2 .and. .not. lr%hybrid_on) then
        call prepare_sdc(lr%d%n1,lr%d%n2,lr%d%n3,&
          w%modul1,w%modul2,w%modul3,w%af,w%bf,w%cf,w%ef,hx,hy,hz)
     end if
     !	initializes the wavelet scaling coefficients	
     call wscal_init_per(scal,hx,hy,hz,cprecr)


     if (lr%hybrid_on) then
        do idx=1,ncplx
           !b=x
           call vcopy(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,x(1,idx),1,b(1,idx),1) 
           
           call prec_fft_fast(lr%d%n1,lr%d%n2,lr%d%n3,&
                lr%wfd%nseg_c,lr%wfd%nvctr_c,lr%wfd%nseg_f,lr%wfd%nvctr_f,&
                lr%wfd%keygloc,lr%wfd%keyvloc, &
                cprecr,hx,hy,hz,x(1,idx),&
                w%kern_k1,w%kern_k2,w%kern_k3,w%z1,w%z3,w%x_c,&
                nd1,nd2,nd3,n1f,n1b,n3f,n3b,nd1f,nd1b,nd3f,nd3b)
        end do
     else
        ! Array sizes for the real-to-complex FFT: note that n1(there)=n1(here)+1
        ! and the same for lr%d%n2,n3.

        do idx=1,ncplx
           !	scale the r.h.s. that is also the scaled input guess :
           !	b'=D^{-1/2}b
           call wscal_per_self(lr%wfd%nvctr_c,lr%wfd%nvctr_f,scal,&
                x(1,idx),x(lr%wfd%nvctr_c+1,idx))
           !b=x
           call vcopy(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,x(1,idx),1,b(1,idx),1) 

           !if GPU is swithced on and there is no call to GPU preconditioner
           !do not do the FFT preconditioning (not valid anymore)
           if (.not. GPUconv .or. .true.) then
              !	compute the input guess x via a Fourier transform in a cubic box.
              !	Arrays psifscf and ww serve as work arrays for the Fourier
              fac=1.0_gp/scal(0)**2
              call prec_fft_c(lr%d%n1,lr%d%n2,lr%d%n3,lr%wfd%nseg_c,&
                   lr%wfd%nvctr_c,lr%wfd%nseg_f,lr%wfd%nvctr_f,&
                   lr%wfd%keygloc,lr%wfd%keyvloc, &
                   cprecr,hx,hy,hz,x(1,idx),&
                   w%psifscf(1),w%psifscf(lr%d%n1+2),&
                   w%psifscf(lr%d%n1+lr%d%n2+3),w%ww(1),w%ww(nd1b*nd2*nd3*4+1),&
                   w%ww(nd1b*nd2*nd3*4+nd1*nd2*nd3f*4+1),&
                   nd1,nd2,nd3,n1f,n1b,n3f,n3b,nd1f,nd1b,nd3f,nd3b,fac)
           end if
        end do
     end if


  else if (lr%geocode == 'S') then

     if (ncplx == 1) then
        call prepare_sdc_slab(lr%d%n1,lr%d%n3,w%modul1,w%modul3,&
          w%af,w%bf,w%cf,w%ef,hx,hy,hz)
     end if
    
     !	initializes the wavelet scaling coefficients	
     call wscal_init_per(scal,hx,hy,hz,cprecr)
    
     do idx=1,ncplx

        !recently added
        !	scale the r.h.s. that is also the scaled input guess :
        !	b'=D^{-1/2}b
        call wscal_per_self(lr%wfd%nvctr_c,lr%wfd%nvctr_f,scal,&
             x(1,idx),x(lr%wfd%nvctr_c+1,idx))
        !end of that

        !b=x
        call vcopy(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,x(1,idx),1,b(1,idx),1) 
        
        !	compute the input guess x via a Fourier transform in a cubic box.
        !	Arrays psifscf and ww serve as work arrays for the Fourier
        call prec_fft_slab_fast(lr%d%n1,lr%d%n2,lr%d%n3,lr%wfd%nseg_c,lr%wfd%nvctr_c,&
             lr%wfd%nseg_f,lr%wfd%nvctr_f,lr%wfd%keygloc,lr%wfd%keyvloc, &
             cprecr,hx,hy,hz,x(1,idx),&
             w%psifscf(1),w%psifscf(lr%d%n1+2),w%ww(1),&
             w%ww(2*((lr%d%n1+1)/2+1)*(lr%d%n2+1)*(lr%d%n3+1)+1))

        !we will probably have to rescale x by fac=1.0_gp/scal(0)**2
        call dscal(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,1.0_gp/scal(0)**2,x(1,idx),1)
        
     end do

  end if

  call f_release_routine()
  
END SUBROUTINE precondition_preconditioner


subroutine allocate_work_arrays(geocode,hybrid_on,ncplx,d,w)
  use module_base
  use module_types
  implicit none
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  logical, intent(in) :: hybrid_on
  integer, intent(in) :: ncplx
  type(grid_dimensions), intent(in) :: d
  type(workarr_precond), intent(out) :: w
  !local variables
  character(len=*), parameter :: subname='allocate_work_arrays'
  integer, parameter :: lowfil=-14,lupfil=14
  integer :: nd1,nd2,nd3
  integer :: n1f,n3f,n1b,n3b,nd1f,nd3f,nd1b,nd3b
  integer :: nf

  if (geocode == 'F') then

     nf=(d%nfu1-d%nfl1+1)*(d%nfu2-d%nfl2+1)*(d%nfu3-d%nfl3+1)
     !allocate work arrays
     w%xpsig_c = f_malloc_ptr((/ 0.to.d%n1, 0.to.d%n2, 0.to.d%n3 /),id='w%xpsig_c')
     w%xpsig_f = f_malloc_ptr((/ 1.to.7, d%nfl1.to.d%nfu1, d%nfl2.to.d%nfu2, d%nfl3.to.d%nfu3 /),id='w%xpsig_f')
     w%ypsig_c = f_malloc_ptr((/ 0.to.d%n1, 0.to.d%n2, 0.to.d%n3 /),id='w%ypsig_c')
     w%ypsig_f = f_malloc_ptr((/ 1.to.7, d%nfl1.to.d%nfu1, d%nfl2.to.d%nfu2, d%nfl3.to.d%nfu3 /),id='w%ypsig_f')

     w%x_f1 = f_malloc_ptr(nf,id='w%x_f1')
     w%x_f2 = f_malloc_ptr(nf,id='w%x_f2')
     w%x_f3 = f_malloc_ptr(nf,id='w%x_f3')
    
  else if (geocode == 'P') then
     
     if (hybrid_on) then
          
        call dimensions_fft(d%n1,d%n2,d%n3,&
             nd1,nd2,nd3,n1f,n3f,n1b,n3b,nd1f,nd3f,nd1b,nd3b)

        nf=(d%nfu1-d%nfl1+1)*(d%nfu2-d%nfl2+1)*(d%nfu3-d%nfl3+1)

        w%kern_k1 = f_malloc_ptr(0.to.d%n1,id='w%kern_k1')
        w%kern_k2 = f_malloc_ptr(0.to.d%n2,id='w%kern_k2')
        w%kern_k3 = f_malloc_ptr(0.to.d%n3,id='w%kern_k3')
        w%z1 = f_malloc_ptr((/ 2, nd1b, nd2, nd3, 2 /),id='w%z1')
        w%z3 = f_malloc_ptr((/ 2, nd1, nd2, nd3f, 2 /),id='w%z3')
        w%x_c = f_malloc_ptr((/ 0.to.d%n1, 0.to.d%n2, 0.to.d%n3 /),id='w%x_c')

        w%x_f = f_malloc_ptr((/ 1.to.7, d%nfl1.to.d%nfu1, d%nfl2.to.d%nfu2, d%nfl3.to.d%nfu3 /),id='w%x_f')
        w%x_f1 = f_malloc_ptr(nf,id='w%x_f1')
        w%x_f2 = f_malloc_ptr(nf,id='w%x_f2')
        w%x_f3 = f_malloc_ptr(nf,id='w%x_f3')
        w%y_f = f_malloc_ptr((/ 1.to.7, d%nfl1.to.d%nfu1, d%nfl2.to.d%nfu2, d%nfl3.to.d%nfu3 /),id='w%y_f')
        w%ypsig_c = f_malloc_ptr((/ 0.to.d%n1, 0.to.d%n2, 0.to.d%n3 /),id='w%ypsig_c')


     else 

        if (ncplx == 1) then
           !periodic, not k-points
           w%modul1 = f_malloc_ptr(lowfil.to.d%n1+lupfil,id='w%modul1')
           w%modul2 = f_malloc_ptr(lowfil.to.d%n2+lupfil,id='w%modul2')
           w%modul3 = f_malloc_ptr(lowfil.to.d%n3+lupfil,id='w%modul3')
           w%af = f_malloc_ptr((/ lowfil.to.lupfil, 1.to.3 /),id='w%af')
           w%bf = f_malloc_ptr((/ lowfil.to.lupfil, 1.to.3 /),id='w%bf')
           w%cf = f_malloc_ptr((/ lowfil.to.lupfil, 1.to.3 /),id='w%cf')
           w%ef = f_malloc_ptr((/ lowfil.to.lupfil, 1.to.3 /),id='w%ef')
        end if

        w%psifscf = f_malloc_ptr(ncplx*(2*d%n1+2)*(2*d%n2+2)*(2*d%n3+2),id='w%psifscf')
        w%ww = f_malloc_ptr(ncplx*(2*d%n1+2)*(2*d%n2+2)*(2*d%n3+2),id='w%ww')

     end if

  else if (geocode == 'S') then

     if (ncplx == 1) then
        w%modul1 = f_malloc_ptr(lowfil.to.d%n1+lupfil,id='w%modul1')
        w%modul3 = f_malloc_ptr(lowfil.to.d%n3+lupfil,id='w%modul3')
        w%af = f_malloc_ptr((/ lowfil.to.lupfil, 1.to.3 /),id='w%af')
        w%bf = f_malloc_ptr((/ lowfil.to.lupfil, 1.to.3 /),id='w%bf')
        w%cf = f_malloc_ptr((/ lowfil.to.lupfil, 1.to.3 /),id='w%cf')
        w%ef = f_malloc_ptr((/ lowfil.to.lupfil, 1.to.3 /),id='w%ef')
     end if
        
     w%psifscf = f_malloc_ptr(ncplx*(2*d%n1+2)*(2*d%n2+16)*(2*d%n3+2),id='w%psifscf')
     w%ww = f_malloc_ptr(ncplx*(2*d%n1+2)*(2*d%n2+16)*(2*d%n3+2),id='w%ww')

  end if

END SUBROUTINE allocate_work_arrays


subroutine memspace_work_arrays_precond(geocode,hybrid_on,ncplx,d,memwork)
  use module_base
  use module_types
  implicit none
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  logical, intent(in) :: hybrid_on
  integer, intent(in) :: ncplx
  type(grid_dimensions), intent(in) :: d
  integer(kind=8), intent(out) :: memwork
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  integer :: nd1,nd2,nd3
  integer :: n1f,n3f,n1b,n3b,nd1f,nd3f,nd1b,nd3b
  integer :: nf


  if (geocode == 'F') then

     nf=(d%nfu1-d%nfl1+1)*(d%nfu2-d%nfl2+1)*(d%nfu3-d%nfl3+1)

     memwork=2*(d%n1+1)*(d%n2+1)*(d%n3+1)+2*7*(d%nfu1-d%nfl1+1)*(d%nfu2-d%nfl2+1)*(d%nfu3-d%nfl3+1)+3*nf
     
    
  else if (geocode == 'P') then
     
     if (hybrid_on) then
          
        call dimensions_fft(d%n1,d%n2,d%n3,&
             nd1,nd2,nd3,n1f,n3f,n1b,n3b,nd1f,nd3f,nd1b,nd3b)

        nf=(d%nfu1-d%nfl1+1)*(d%nfu2-d%nfl2+1)*(d%nfu3-d%nfl3+1)

        memwork=(d%n1+1)+(d%n2+1)+(d%n3+1)+2*nd1b*nd2*nd3*2+2*nd1*nd2*nd3f*2+&
             (d%n1+1)*(d%n2+1)*(d%n3+1)+2*7*(d%nfu1-d%nfl1+1)*(d%nfu2-d%nfl2+1)*(d%nfu3-d%nfl3+1)+3*nf

     else 

        memwork=0
        if (ncplx == 1) then
           memwork=d%n1+d%n2+d%n3+15*(lupfil-lowfil+1)
        end if
        memwork=memwork+2*ncplx*(2*d%n1+2)*(2*d%n2+2)*(2*d%n3+2)

     end if

  else if (geocode == 'S') then

     memwork=0
     if (ncplx == 1) then
        memwork=d%n1+d%n3+14*(lupfil-lowfil+1)
     end if
     memwork=memwork+2*ncplx*(2*d%n1+2)*(2*d%n2+16)*(2*d%n3+2)
  end if

END SUBROUTINE memspace_work_arrays_precond


subroutine deallocate_work_arrays(geocode,hybrid_on,ncplx,w)
  use module_base
  use module_types
  implicit none
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  logical, intent(in) :: hybrid_on
  integer, intent(in) :: ncplx
  type(workarr_precond), intent(inout) :: w
  !local variables
  character(len=*), parameter :: subname='deallocate_work_arrays'

  if (geocode == 'F') then

     call f_free_ptr(w%xpsig_c)
     call f_free_ptr(w%ypsig_c)
     call f_free_ptr(w%xpsig_f)
     call f_free_ptr(w%ypsig_f)
     call f_free_ptr(w%x_f1)
     call f_free_ptr(w%x_f2)
     call f_free_ptr(w%x_f3)

  else if ((geocode == 'P' .and. .not. hybrid_on) .or. geocode == 'S') then

     if (ncplx == 1) then
        call f_free_ptr(w%modul1)
        if (geocode /= 'S') then
           call f_free_ptr(w%modul2)
        end if
        call f_free_ptr(w%modul3)
        call f_free_ptr(w%af)
        call f_free_ptr(w%bf)
        call f_free_ptr(w%cf)
        call f_free_ptr(w%ef)
     end if

     call f_free_ptr(w%psifscf)
     call f_free_ptr(w%ww)

  else if (geocode == 'P' .and. hybrid_on) then

     call f_free_ptr(w%z1)
     call f_free_ptr(w%z3)
     call f_free_ptr(w%kern_k1)
     call f_free_ptr(w%kern_k2)
     call f_free_ptr(w%kern_k3)
     call f_free_ptr(w%x_c)
     call f_free_ptr(w%x_f)
     call f_free_ptr(w%x_f1)
     call f_free_ptr(w%x_f2)
     call f_free_ptr(w%x_f3)
     call f_free_ptr(w%y_f)
     call f_free_ptr(w%ypsig_c)


  end if

END SUBROUTINE deallocate_work_arrays


subroutine precond_locham(ncplx,lr,hx,hy,hz,kx,ky,kz,&
     cprecr,x,y,w,scal)! y:=Ax
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: ncplx
  real(gp), intent(in) :: hx,hy,hz,cprecr,kx,ky,kz
  type(locreg_descriptors), intent(in) :: lr
  real(gp), dimension(0:7), intent(in) :: scal
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,ncplx), intent(in) ::  x
  type(workarr_precond), intent(inout) :: w
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,ncplx), intent(out) ::  y
  !local variables
  logical :: sseprecond=.false.
  integer :: idx,nf,isegf,ipsif

  isegf=lr%wfd%nseg_c+min(1,lr%wfd%nseg_f)
  ipsif=lr%wfd%nvctr_c+min(1,lr%wfd%nvctr_f)

  if (lr%geocode == 'F') then
     do idx=1,ncplx

        if (sseprecond) then
           call uncompress_standard_scal(lr%d,lr%wfd,scal,&
                lr%wfd%keyvloc(1),lr%wfd%keyvloc(isegf),&
                lr%wfd%keygloc(1,1),lr%wfd%keygloc(1,isegf),&
                x(1,idx),x(ipsif,idx),&
                w%xpsig_c,w%xpsig_f)
!commented out, not working correctly        
!!$           call Convolkinetic_SSE(lr%d%n1,lr%d%n2,lr%d%n3, &
!!$                lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3,  &
!!$                cprecr,hx,&
!!$                lr%bounds%kb%ibyz_c,lr%bounds%kb%ibxz_c,lr%bounds%kb%ibxy_c,&
!!$                lr%bounds%kb%ibyz_f,lr%bounds%kb%ibxz_f,lr%bounds%kb%ibxy_f,&
!!$                w%xpsig_c,w%xpsig_f,w%ypsig_c,w%ypsig_f)
           
           call compress_standard_scal(lr%d,lr%wfd,scal,&
                lr%wfd%keyvloc(1),lr%wfd%keyvloc(isegf),&
                lr%wfd%keygloc(1,1),lr%wfd%keygloc(1,isegf),&
                w%ypsig_c,w%ypsig_f,&
                y(1,idx),y(ipsif,idx))

        else

           call calc_grad_reza(lr%d%n1,lr%d%n2,lr%d%n3,&
                lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3, &
                lr%wfd%nseg_c,lr%wfd%nvctr_c,lr%wfd%keygloc,lr%wfd%keyvloc,&
                lr%wfd%nseg_f,lr%wfd%nvctr_f,&
                lr%wfd%keygloc(1,lr%wfd%nseg_c+min(1,lr%wfd%nseg_f)),&
                lr%wfd%keyvloc(lr%wfd%nseg_c+min(1,lr%wfd%nseg_f)), &
                scal,cprecr,hx,&
                lr%bounds%kb%ibyz_c,lr%bounds%kb%ibxz_c,lr%bounds%kb%ibxy_c,&
                lr%bounds%kb%ibyz_f,lr%bounds%kb%ibxz_f,lr%bounds%kb%ibxy_f,&
                x(1,idx),x(lr%wfd%nvctr_c+min(1,lr%wfd%nvctr_f),idx),&
                y(1,idx),y(lr%wfd%nvctr_c+min(1,lr%wfd%nvctr_f),idx),&
                w%xpsig_c,w%xpsig_f,w%ypsig_c,w%ypsig_f,&
                w%x_f1,w%x_f2,w%x_f3)
        end if
     end do
  else if (lr%geocode == 'P') then
     if (lr%hybrid_on) then

        nf=(lr%d%nfu1-lr%d%nfl1+1)*(lr%d%nfu2-lr%d%nfl2+1)*(lr%d%nfu3-lr%d%nfl3+1)
        do idx=1,ncplx
           call apply_hp_hyb(lr%d%n1,lr%d%n2,lr%d%n3,&
                lr%wfd%nseg_c,lr%wfd%nvctr_c,lr%wfd%nseg_f,lr%wfd%nvctr_f,&
                lr%wfd%keygloc,lr%wfd%keyvloc, &
                cprecr,hx,hy,hz,x(1,idx),y(1,idx),&
                w%x_f,w%x_c,w%x_f1,w%x_f2,w%x_f3,w%y_f,w%ypsig_c,&
                lr%d%nfl1,lr%d%nfl2,lr%d%nfl3,lr%d%nfu1,lr%d%nfu2,lr%d%nfu3,nf,&
                lr%bounds%kb%ibyz_f,lr%bounds%kb%ibxz_f,lr%bounds%kb%ibxy_f)
        end do
     else
        if (ncplx == 1) then
           call apply_hp_scal(lr%d%n1,lr%d%n2,lr%d%n3,&
                lr%wfd%nseg_c,lr%wfd%nvctr_c,lr%wfd%nseg_f,&
                lr%wfd%nvctr_f,lr%wfd%keygloc,lr%wfd%keyvloc, &
                cprecr,x,y,w%psifscf,w%ww,w%modul1,w%modul2,w%modul3,&
                w%af,w%bf,w%cf,w%ef,scal) 
        else
           call apply_hp_per_k(lr%d%n1,lr%d%n2,lr%d%n3,&
                lr%wfd%nseg_c,lr%wfd%nvctr_c,lr%wfd%nseg_f,&
                lr%wfd%nvctr_f,lr%wfd%keygloc,lr%wfd%keyvloc, &
                !cprecr,hx,hy,hz,0.0_gp,0.0_gp,0.0_gp,x,y,w%psifscf,w%ww,scal) 
                cprecr,hx,hy,hz,kx,ky,kz,x,y,w%psifscf,w%ww,scal) 
        end if
     end if
  else if (lr%geocode == 'S') then
     if (ncplx == 1) then
        call apply_hp_slab_sd_scal(lr%d%n1,lr%d%n2,lr%d%n3,&
             lr%wfd%nseg_c,lr%wfd%nvctr_c,lr%wfd%nseg_f,&
             lr%wfd%nvctr_f,lr%wfd%keygloc,lr%wfd%keyvloc, &
             cprecr,x,y,w%psifscf,w%ww,w%modul1,w%modul3,&
             w%af,w%bf,w%cf,w%ef,scal)
     else
        call apply_hp_slab_k(lr%d%n1,lr%d%n2,lr%d%n3,&
             lr%wfd%nseg_c,lr%wfd%nvctr_c,lr%wfd%nseg_f,&
             lr%wfd%nvctr_f,lr%wfd%keygloc,lr%wfd%keyvloc, &
             cprecr,hx,hy,hz,kx,ky,kz,x,y,w%psifscf,w%ww,scal) 

     end if
   end if
END SUBROUTINE precond_locham


!> ypsi = @f$(1/2) \nabla^2 xpsi + cprecr xpsi@f$
subroutine calc_grad_reza(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
     nseg_c,nvctr_c,keyg_c,keyv_c,nseg_f,nvctr_f,keyg_f,keyv_f, &
     scal,cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,&
     xpsi_c,xpsi_f,ypsi_c,ypsi_f,&
     xpsig_c,xpsig_f,ypsig_c,ypsig_f,x_f1,x_f2,x_f3)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  integer, intent(in) :: nseg_c,nvctr_c,nseg_f,nvctr_f
  real(wp), intent(in) :: cprecr
  real(gp), intent(in) :: hgrid
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(nseg_f), intent(in) :: keyv_f
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  integer, dimension(2,nseg_f), intent(in) :: keyg_f
  integer, dimension(2,0:n2,0:n3), intent(in) :: ibyz_c,ibyz_f
  integer, dimension(2,0:n1,0:n3), intent(in) :: ibxz_c,ibxz_f
  integer, dimension(2,0:n1,0:n2), intent(in) :: ibxy_c,ibxy_f
  real(wp), dimension(0:3), intent(in) :: scal
  real(wp), dimension(nvctr_c), intent(in) :: xpsi_c
  real(wp), dimension(7,nvctr_f), intent(in) :: xpsi_f
  real(wp), dimension(nvctr_c), intent(out) :: ypsi_c
  real(wp), dimension(7,nvctr_f), intent(out) :: ypsi_f
  real(wp), dimension(0:n1,0:n2,0:n3), intent(inout) :: xpsig_c,ypsig_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(inout) :: xpsig_f,ypsig_f
  real(wp), dimension(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(inout) :: x_f1
  real(wp), dimension(nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), intent(inout) :: x_f2
  real(wp), dimension(nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), intent(inout) :: x_f3

  call uncompress_forstandard(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
       nseg_c,nvctr_c,keyg_c,keyv_c,  & 
       nseg_f,nvctr_f,keyg_f,keyv_f,  & 
       scal,xpsi_c,xpsi_f,xpsig_c,xpsig_f,x_f1,x_f2,x_f3)

!!  ypsig_c=xpsig_c
!!  ypsig_f=xpsig_f
  call Convolkinetic(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
       cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,xpsig_c,&
       xpsig_f,ypsig_c,ypsig_f,x_f1,x_f2,x_f3)

  call compress_forstandard(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
       nseg_c,nvctr_c,keyg_c,keyv_c,  & 
       nseg_f,nvctr_f,keyg_f,keyv_f,  & 
       scal,ypsig_c,ypsig_f,ypsi_c,ypsi_f)

END SUBROUTINE calc_grad_reza


subroutine prec_diag(n1,n2,n3,hgrid,nseg_c,nvctr_c,nvctr_f,&
     keyg_c,keyv_c,hpsi_c,hpsi_f,c,scal,a2,b2)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nseg_c,nvctr_c,nvctr_f
  real(wp), intent(in) :: c,a2,b2
  real(gp), intent(in) :: hgrid
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  real(wp), dimension(0:3), intent(in) :: scal
  real(wp), dimension(nvctr_c), intent(inout) :: hpsi_c
  real(wp), dimension(7,nvctr_f), intent(inout) :: hpsi_f
  !local variables
  character(len=*), parameter :: subname='prec_diag'
  real(gp), parameter ::atomic_length=2.0_gp,fac_len=2.0_gp
  integer :: num_trans,n2_nt,nd1,nd2,nd3,iseg,jj,j0,ii,i3,i2,i
  integer :: nn1,nn2,nn3,nnn1,nnn2,nnn3,i0,i1,j1
  real(wp) :: h0,h1,h2,h3,fac_h
  real(wp), dimension(:,:,:), allocatable :: hpsip

  !      number of sweeps in wavelet transformation
  !      the biggest scaling function step: atomic_length*fac_len
  !      (not just atomic_length, because so it is better in practice) 
  num_trans=nint(log(atomic_length*fac_len/hgrid)/log(2.0_gp))
  n2_nt=2**num_trans
  !write(*,'(1x,a)') 'number of wavelet transforms (sweeps)',num_trans

  ! find right leading dimensions for array

  !       nd1+1 is the multiple of n2_n
  !       which is closest to n1+1 from above. 
  nd1=ceiling( real(n1+1,kind=8)/real(n2_nt,kind=8)) *n2_nt-1
  !       the same for nd2,nd3.
  nd2=ceiling( real(n2+1,kind=8)/real(n2_nt,kind=8)) *n2_nt-1
  nd3=ceiling( real(n3+1,kind=8)/real(n2_nt,kind=8)) *n2_nt-1

  !write(*,'(3(1x,a,i0))')'nd1=',nd1,'nd2=',nd2,'nd3=',nd3

  hpsip = f_malloc((/ 0.to.nd1, 0.to.nd2, 0.to.nd3 /),id='hpsip')

  hpsip=0.0_wp

  ! coarse part
  !$omp parallel default(shared)&
  !$omp private(iseg,jj,j0,j1,ii,i3,i2,i1,i,i0)
  !$omp do !!!!schedule(static,1)
  do iseg=1,nseg_c
     jj=keyv_c(iseg)
     j0=keyg_c(1,iseg)
     j1=keyg_c(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        hpsip(i,i2,i3)=hpsi_c(i-i0+jj)
     enddo
  enddo
  !$omp enddo
  !$omp end parallel

  fac_h=real(1.0_gp/((hgrid*real(n2_nt,gp))**2),wp)

  h0=    1.5_wp*a2*fac_h
  h1=(a2+b2*.5d0)*fac_h
  h2=(a2*.5_wp+b2)*fac_h
  h3=    1.5_wp*b2*fac_h

  !       forward transform the coarse scaling functions num_trans times
  call ana_repeated_per(nd1,nd2,nd3,hpsip,num_trans,nn1,nn2,nn3) 

  nnn1=nn1
  nnn2=nn2
  nnn3=nn3 

  !       diagonally precondition the resulting coarse wavelets
  call precond_proper(nd1,nd2,nd3,hpsip,num_trans,nnn1,nnn2,nnn3,h0,h1,h2,h3,c)

  hpsip=hpsip/scal(0) ! apply (wscal)^(-1)

  !       backward transform the coarse scaling functions num_trans times
  call syn_repeated_per(nd1,nd2,nd3,hpsip,num_trans,nn1,nn2,nn3)

  !       diagonally precondition the fine wavelets
  !$omp parallel default(shared)&
  !$omp private(i)
  !$omp do !!!!schedule(static,1)
  do i=1,nvctr_f
     hpsi_f(1,i)=hpsi_f(1,i)*scal(1)
     hpsi_f(2,i)=hpsi_f(2,i)*scal(1)
     hpsi_f(4,i)=hpsi_f(4,i)*scal(1)

     hpsi_f(3,i)=hpsi_f(3,i)*scal(2)
     hpsi_f(5,i)=hpsi_f(5,i)*scal(2)
     hpsi_f(6,i)=hpsi_f(6,i)*scal(2)

     hpsi_f(7,i)=hpsi_f(7,i)*scal(3)
  enddo
  !$omp enddo
  !$omp end parallel

  ! coarse part
  !$omp parallel default(shared)&
  !$omp private(iseg,jj,j0,j1,ii,i3,i2,i1,i0,i)
  !$omp do !!!!schedule(static,1)
  do iseg=1,nseg_c
     jj=keyv_c(iseg)
     j0=keyg_c(1,iseg)
     j1=keyg_c(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        hpsi_c(i-i0+jj)=hpsip(i,i2,i3)
     enddo
  enddo
  !$omp enddo
  !$omp end parallel

  call f_free(hpsip)

END SUBROUTINE prec_diag


subroutine precond_proper(nd1,nd2,nd3,x,num_trans,n1,n2,n3,h0,h1,h2,h3,eps)
  use module_base
  implicit none
  integer, intent(in) :: nd1,nd2,nd3,num_trans
  integer, intent(inout) :: n1,n2,n3
  real(wp), intent(in) :: eps,h0
  real(wp), intent(inout) :: h1,h2,h3
  real(wp), dimension(0:nd1,0:nd2,0:nd3), intent(inout) :: x
  !local variables
  integer :: i_trans,n1p,n2p,n3p,n1pp,n2pp,n3pp,i1,i2,i3,i1p,i2p,i3p
  real(wp) :: f0,f1,f2,f3


  do i_trans=1,num_trans
     n1p=2*(n1+1)-1
     n2p=2*(n2+1)-1
     n3p=2*(n3+1)-1

     if (n1p.gt.nd1) stop 'n1 beyond borders'
     if (n2p.gt.nd2) stop 'n2 beyond borders'
     if (n3p.gt.nd3) stop 'n3 beyond borders'

     n1pp=n1+1
     n2pp=n2+1
     n3pp=n3+1

     f1=1.0_wp/(h1+eps)
     f2=1.0_wp/(h2+eps)
     f3=1.0_wp/(h3+eps)       

     if (i_trans == 1) then 

        f0=1.d0/(h0+eps)

     !$omp parallel default(shared)&   !*
     !$omp private(i3,i3p,i2,i2p,i1,i1p)
     !$omp do !!!!schedule(static,1)
        do i3=0,n3
           i3p=i3+n3pp
           do i2=0,n2
              i2p=i2+n2pp
              do i1=0,n1
                 i1p=i1+n1pp

                 x(i1,i2,i3)=x(i1,i2,i3)*f0

                 x(i1p,i2,i3)=x(i1p,i2,i3)*f1
                 x(i1,i2p,i3)=x(i1,i2p,i3)*f1
                 x(i1,i2,i3p)=x(i1,i2,i3p)*f1

                 x(i1p,i2p,i3)=x(i1p,i2p,i3)*f2
                 x(i1,i2p,i3p)=x(i1,i2p,i3p)*f2
                 x(i1p,i2,i3p)=x(i1p,i2,i3p)*f2

                 x(i1p,i2p,i3p)=x(i1p,i2p,i3p)*f3

              enddo
           enddo
        enddo
     !$omp enddo
     !$omp end parallel

     else

     !$omp parallel default(shared)&   !*
     !$omp private(i3,i3p,i2,i2p,i1,i1p)
     !$omp do !!!!schedule(static,1)
        do i3=0,n3
           i3p=i3+n3pp
           do i2=0,n2
              i2p=i2+n2pp
              do i1=0,n1
                 i1p=i1+n1pp

                 x(i1p,i2,i3)=x(i1p,i2,i3)*f1
                 x(i1,i2p,i3)=x(i1,i2p,i3)*f1
                 x(i1,i2,i3p)=x(i1,i2,i3p)*f1

                 x(i1p,i2p,i3)=x(i1p,i2p,i3)*f2
                 x(i1,i2p,i3p)=x(i1,i2p,i3p)*f2
                 x(i1p,i2,i3p)=x(i1p,i2,i3p)*f2

                 x(i1p,i2p,i3p)=x(i1p,i2p,i3p)*f3

              enddo
           enddo
        enddo
     !$omp enddo
     !$omp end parallel

     endif

     n1=n1p
     n2=n2p
     n3=n3p

     h1=h1*4.0_wp
     h2=h2*4.0_wp
     h3=h3*4.0_wp

  enddo

END SUBROUTINE precond_proper


!> Solves (KE+cprecr*I)*xx=yy by conjugate gradient method
!! hpsi is the right hand side on input and the solution on output
!!
!! The input guess consists of diagonal preconditioning of the original gradient.
!! In contrast to older version, not only the wavelet part and the scfunction
!! part are multiplied by different factors, but the scfunction part is 
!! subjected to wavelet analysis with periodic boundaries. Then the wavelets
!! on different scales are multiplied by different factors and backward wavelet 
!! transformed to scaling functions.
!!
!! The new input guess is turned on if the parameter INGUESS_ON
!! has value .TRUE.
!! @warning
!!  This routine is sensitive in OpenMP versus the number of threads.
subroutine precong(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, &
     nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
     ncong,cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,hpsi)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  integer, intent(in) :: nseg_c,nvctr_c,nseg_f,nvctr_f,ncong
  real(gp), intent(in) :: hgrid
  real(dp), intent(in) :: cprecr
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  integer, dimension(2,0:n2,0:n3), intent(in) :: ibyz_c,ibyz_f
  integer, dimension(2,0:n1,0:n3), intent(in) :: ibxz_c,ibxz_f
  integer, dimension(2,0:n1,0:n2), intent(in) :: ibxy_c,ibxy_f
  real(wp), dimension(nvctr_c+7*nvctr_f), intent(inout) :: hpsi
  !local variables
  character(len=*), parameter :: subname='precong'
  logical, parameter :: inguess_on=.true.
  !       wavelet and scaling function second derivative filters
  real(wp), parameter :: b2=24.8758460293923314_wp, a2=3.55369228991319019_wp
  integer :: i,icong
  real(wp) :: fac_h,h0,h1,h2,h3,tt,alpha1,alpha2,alpha,beta1,beta2,beta,aa1,aa2
  real(wp), dimension(0:3) :: scal
  real(wp), dimension(:), allocatable :: rpsi,ppsi,wpsi
  real(wp), dimension(:,:,:,:), allocatable :: xpsig_f,ypsig_f
  real(wp), dimension(:,:,:), allocatable :: xpsig_c,ypsig_c,x_f1,x_f2,x_f3

  rpsi = f_malloc(nvctr_c+7*nvctr_f,id='rpsi')
  ppsi = f_malloc(nvctr_c+7*nvctr_f,id='ppsi')
  wpsi = f_malloc(nvctr_c+7*nvctr_f,id='wpsi')

!!  !array of initial wavefunction
!!  allocate(spsi(nvctr_c+7*nvctr_f),stat=i_stat)
!!  call memocc(i_stat,spsi,'spsi',subname)
!!  do i=1,nvctr_c+7*nvctr_f
!!     spsi(i)=hpsi(i)
!!  enddo

  fac_h=1.0_wp/real(hgrid,wp)**2
  h0=    1.5_wp*a2*fac_h
  h1=(a2+b2*.5_wp)*fac_h
  h2=(a2*.5_wp+b2)*fac_h
  h3=    1.5_wp*b2*fac_h

  scal(0)=sqrt(1.0_wp/(h0+cprecr)) 
  scal(1)=sqrt(1.0_wp/(h1+cprecr)) 
  scal(2)=sqrt(1.0_wp/(h2+cprecr)) 
  scal(3)=sqrt(1.0_wp/(h3+cprecr))

  if (inguess_on) then
     !          the right hand side is temporarily stored in the rpsi array
     !rpsi=hpsi           
     call vcopy(nvctr_c+7*nvctr_f,hpsi(1),1,rpsi(1),1) 
     !          and preconditioned with d^{-1/2} as usual:
     call  wscalv(nvctr_c,nvctr_f,scal,rpsi,rpsi(nvctr_c+1))

     !          hpsi is now diagonally preconditioned with alexey's old preconditioner;
     !          inside the diagonal preconditioner a factor of d^{1/2} was added
     !          to make the overall factor d^{-1/2} again
     call prec_diag(n1,n2,n3,hgrid,nseg_c,nvctr_c,nvctr_f,&
          keyg,keyv,hpsi,hpsi(nvctr_c+1),cprecr,scal,a2,b2)
  else
     !          assume as input guess x=y
     !          hpsi is preconditioned with d^{-1/2} as usual
     call  wscalv(nvctr_c,nvctr_f,scal,hpsi,hpsi(nvctr_c+1))
  endif

  !allocate work arrays
  xpsig_c = f_malloc((/ 0.to.n1, 0.to.n2, 0.to.n3 /),id='xpsig_c')
  xpsig_f = f_malloc((/ 1.to.7, nfl1.to.nfu1, nfl2.to.nfu2, nfl3.to.nfu3 /),id='xpsig_f')
  ypsig_c = f_malloc((/ 0.to.n1, 0.to.n2, 0.to.n3 /),id='ypsig_c')
  ypsig_f = f_malloc((/ 1.to.7, nfl1.to.nfu1, nfl2.to.nfu2, nfl3.to.nfu3 /),id='ypsig_f')

  x_f1 = f_malloc((/ nfl1.to.nfu1, nfl2.to.nfu2, nfl3.to.nfu3 /),id='x_f1')
  x_f2 = f_malloc((/ nfl2.to.nfu2, nfl1.to.nfu1, nfl3.to.nfu3 /),id='x_f2')
  x_f3 = f_malloc((/ nfl3.to.nfu3, nfl1.to.nfu1, nfl2.to.nfu2 /),id='x_f3')
  
  !initalize to zero the work arrays, probably not needed
  call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1),x_f1(nfl1,nfl2,nfl3))
  call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1),x_f2(nfl2,nfl1,nfl3))
  call to_zero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1),x_f3(nfl3,nfl1,nfl2))

  call to_zero((n1+1)*(n2+1)*(n3+1),xpsig_c(0,0,0))
  call to_zero(7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1),xpsig_f(1,nfl1,nfl2,nfl3))

  call to_zero((n1+1)*(n2+1)*(n3+1),ypsig_c(0,0,0))
  call to_zero(7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1),ypsig_f(1,nfl1,nfl2,nfl3))
  
  call calc_grad_reza(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, &
       nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), &
       scal,cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,hpsi,&
       hpsi(nvctr_c+1),wpsi,wpsi(nvctr_c+1),&
       xpsig_c,xpsig_f,ypsig_c,ypsig_f,&
       x_f1,x_f2,x_f3)


  IF (INGUESS_ON) THEN 
     !$omp parallel default(shared)&   !*
     !$omp private(i,tt)
     !$omp do !!!!schedule(static,1)
     do i=1,nvctr_c+7*nvctr_f
        tt=wpsi(i)-rpsi(i)  ! rpsi instead of hpsi: alexey
        rpsi(i)=tt
        ppsi(i)=tt
     enddo
     !$omp enddo
     !$omp end parallel
  ELSE
     !$omp parallel default(shared)&   !*
     !$omp private(i,tt)
     !$omp do !!!!schedule(static,1)
     do i=1,nvctr_c+7*nvctr_f
        tt=wpsi(i)-hpsi(i)  ! normal
        rpsi(i)=tt
        ppsi(i)=tt
     enddo
     !$omp enddo
     !$omp end parallel
  ENDIF

  loop_precond: do icong=2,ncong

     call calc_grad_reza(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, &
          nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), &
          scal,cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,&
          ibxy_f,ppsi,ppsi(nvctr_c+1),wpsi,wpsi(nvctr_c+1),&
          xpsig_c,xpsig_f,ypsig_c,ypsig_f,&
          x_f1,x_f2,x_f3)

     alpha1=0.0_wp 
     alpha2=0.0_wp

     !!@warning
     !!This section is very sensitive versus the number of threads

     !$omp parallel default(shared)&   !*
     !$omp private(i,aa1,aa2)
     aa1=0.0_wp
     aa2=0.0_wp
     !$omp do !!!! schedule(static,1)
     do i=1,nvctr_c+7*nvctr_f
        aa1=aa1+rpsi(i)*rpsi(i)
        aa2=aa2+rpsi(i)*wpsi(i)
     enddo
     !$omp enddo

     !$omp critical
     alpha1=alpha1+aa1
     alpha2=alpha2+aa2
     !$omp end critical

     !$omp end parallel
     !write(*,*)icong,alpha1,alpha2

     !residues(icong)=alpha1
     alpha=alpha1/alpha2        

     !write(10+iorb,'(1x,i0,3(1x,1pe24.17))')icong,alpha1,alpha2,alpha

     !$omp parallel default(shared)&
     !$omp private(i)
     !$omp do !!!!schedule (static,1)
     do i=1,nvctr_c+7*nvctr_f
        hpsi(i)=hpsi(i)-alpha*ppsi(i)
        rpsi(i)=rpsi(i)-alpha*wpsi(i)
     end do
     !$omp enddo
     !$omp end parallel

     if (icong >= ncong) exit loop_precond

     beta1=0.0_wp 
     beta2=0.0_wp

     !$omp parallel default(shared)&
     !$omp private(i,aa1,aa2)
     aa1=0.0_wp
     aa2=0.0_wp
     !$omp do !!!! schedule (static,1)
     do i=1,nvctr_c+7*nvctr_f
        aa1=aa1+rpsi(i)*wpsi(i)
        aa2=aa2+ppsi(i)*wpsi(i)
     enddo
     !$omp enddo

     !$omp critical
     beta1=beta1+aa1
     beta2=beta2+aa2
     !$omp end critical

     !$omp end parallel

     beta=beta1/beta2        

     !omp parallel default(shared)&
     !omp private(i)
     !omp do schedule(static,1)
     do i=1,nvctr_c+7*nvctr_f
        ppsi(i)=rpsi(i)-beta*ppsi(i)
     end do
     !omp enddo
     !omp end parallel

  end do loop_precond

  !  D^{-1/2} times solution
  call wscalv(nvctr_c,nvctr_f,scal,hpsi,hpsi(nvctr_c+1))

  !write(*,'(i4,(100(1x,e8.2)))') iorb,(residues(icong),icong=2,ncong)

!!  ! check final residue of original equation
!!  do i=0,3
!!     scal(i)=1.d0
!!  enddo
!!
!!  call CALC_GRAD_REZA(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, &
!!       nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), &
!!       scal,cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,&
!!       ibxy_f,hpsi,hpsi(nvctr_c+1),wpsi,wpsi(nvctr_c+1),&
!!       xpsig_c,xpsig_f,ypsig_c,ypsig_f,&
!!       x_f1,x_f2,x_f3)
!!     
!!  tt=0.d0
!!  do i=1,nvctr_c+7*nvctr_f
!!     tt=tt+(wpsi(i)-spsi(i))**2
!!  enddo
!!  !write(*,'(1x,a,1x,i0,1x,1pe13.6)') 'Precond, final residue',iorb,sqrt(tt)
!!  i_all=-product(shape(spsi))*kind(spsi)
!!  deallocate(spsi,stat=i_stat)
!!  call memocc(i_stat,i_all,'spsi',subname)
!!  ! checkend

  call f_free(rpsi)
  call f_free(ppsi)
  call f_free(wpsi)


  call f_free(xpsig_c)

  call f_free(ypsig_c)

  call f_free(xpsig_f)

  call f_free(ypsig_f)

  call f_free(x_f1)

  call f_free(x_f2)

  call f_free(x_f3)
     
END SUBROUTINE precong
