!> @file
!!   Module used to calculate the exact exchange potential with
!! overlap point to point
!! @author
!!    Copyright (C) 2011-2013 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Module used to calculate op2p (overlap point-to-point) for the exact exchange
module module_exctx_op2p
  use module_base
  use module_types
  use overlap_point_to_point
  implicit none

  private

  !> Public routines
  public :: op2p_exctx_run,op2p_descriptors,orbs_to_attributes,op2p_exctx_set,op2p_exctx_clear

  integer :: ncalls,ncalltot
  real(gp) :: sfac,hfac,hxh,hyh,hzh,eexctX_int
  type(orbitals_data), pointer :: orbs
  type(coulomb_operator), pointer :: pkernel

contains

  !> Initialize internal variables for module
  subroutine OP2P_exctx_set(OP2P,orbs_ext,pkernel_ext)
    implicit none
    type(OP2P_descriptors), intent(inout) :: OP2P
    type(coulomb_operator), intent(inout), target :: pkernel_ext
    type(orbitals_data), intent(in), target :: orbs_ext

    pkernel => pkernel_ext
    orbs => orbs_ext

    eexctX_int=0.0_gp

    hfac=1.0_gp/(product(pkernel%hgrids))

  end subroutine OP2P_exctx_set

  subroutine OP2P_exctx_clear(OP2P,eexctX_ext)
    implicit none
    real(gp), intent(out) :: eexctX_ext
    type(OP2P_descriptors), intent(inout) :: OP2P
    !local variables
    character(len=*), parameter :: subname='OP2P_exctx_clear'
    !print *,'iproc,icount2',iproc,icount
    call free_OP2P_descriptors(OP2P)
    
    nullify(pkernel,orbs)

    !recuperate the exact exchange energy
    eexctX_ext=eexctX_int

  end subroutine OP2P_exctx_clear

  !> Internal operation which calculates the partial densities for any of the orbitals
  !! then uses these information to evaluate the exact exchange operator
  subroutine internal_exctx_operation(istep,iproc,igroup,remote_result,&
     isorb,jsorb,iorbs,jorbs,norbi,norbj,&
     nvctri,nvctrj,nvctri_results,nvctrj_results,&
     psir_i,psir_j,&
     dpsir_i,dpsir_j)
    use module_base
    use module_types
    use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
    implicit none
    logical, intent(in) :: remote_result
    integer, intent(in) :: istep,iproc,igroup,isorb,jsorb,iorbs,jorbs,norbi,norbj
    integer, intent(in) :: nvctri,nvctrj,nvctri_results,nvctrj_results
    real(wp), dimension(nvctri), intent(in) :: psir_i
    real(wp), dimension(nvctrj), intent(in) :: psir_j 
    real(wp), dimension(nvctri_results), intent(inout) :: dpsir_i
    real(wp), dimension(nvctrj_results), intent(inout) :: dpsir_j
    !local variables
    character(len=*), parameter :: subname='internal_exctx_operation'
    integer :: iorb,jorb,iorb_glob,jorb_glob,ierr,ndim,i
    real(gp) :: hfaci,hfacj,hfac2,ehart
    real(wp), dimension(:), allocatable :: rp_ij

    !it is assumed that all the products are done on the same dimension
    ndim=product(pkernel%ndims)

    rp_ij = f_malloc(ndim,id='rp_ij')
    !print *,'norbi,norbj,jorbs,iorbs,isorb,jsorb',norbi,norbj,jorbs,iorbs,isorb,jsorb
    !for the first step do only the upper triangular part
    do iorb=1,norbi
       iorb_glob=iorb+isorb!+iorbs
       hfacj=-sfac*orbs%occup(iorb+isorb)
       do jorb=1,norbj
          jorb_glob=jorb+jsorb!+jorbs
          !first cross-check whether the spin indices are the same
          if (orbs%spinsgn(iorb_glob) /= orbs%spinsgn(jorb_glob)) then
             write(*,*)'ERROR in partitioning the orbitals',&
                  iorb_glob,jorb_glob,igroup,jsorb,iproc
             call MPI_ABORT(bigdft_mpi%mpi_comm,ierr)
          end if
          hfaci=-sfac*orbs%occup(jorb_glob)
          !do it only for upper triangular results
          if (istep /= 0 .or. jorb_glob >= iorb_glob) then
             if (istep == 0 ) then
                do i=1,ndim
                   rp_ij(i)=hfac*psir_i(i+(iorb-1)*ndim)*psir_i(i+(jorb-1)*ndim)
                end do
             else
                do i=1,ndim
                   rp_ij(i)=hfac*psir_i(i+(iorb-1)*ndim)*psir_j(i+(jorb-1)*ndim)
                end do
             end if
             ncalls=ncalls+1                    

             call H_potential('D',pkernel,rp_ij,rp_ij,ehart,0.0_dp,.false.,&
                  quiet='YES')
             !print *,'done',iproc
             !Poisson solver in sequential (commented out for the moment
             !if (iproc == OP2P%iprocref .and. verbose > 1) then
             !   write(*,'(1x,a,i3,a2)')'Exact exchange calculation: ',&
             !        nint(real(ncalls,gp)/real(ncalltot,gp)*100.0_gp),' %'
             !end if

             !this factor is only valid with one k-point (weigth equal to one)
             !can be easily generalised to the k-point case
             hfac2=sfac*orbs%occup(iorb_glob)*orbs%occup(jorb_glob)

             !exact exchange energy
             if (iorb_glob == jorb_glob) then
                eexctX_int=eexctX_int+hfac2*real(ehart,gp)
             else
                !if the result has to be sent away
                if (remote_result .or. istep==0) then
                   eexctX_int=eexctX_int+2.0_gp*hfac2*real(ehart,gp)
                else !otherwise other processors are already calculating it
                   eexctX_int=eexctX_int+hfac2*real(ehart,gp)
                end if
             end if
             !accumulate the results for each of the wavefunctions concerned
             if (istep == 0) then
                do i=1,ndim
                   dpsir_i(i+(iorb-1)*ndim)=dpsir_i(i+(iorb-1)*ndim)+&
                        hfaci*rp_ij(i)*psir_i(i+(jorb-1)*ndim)
                end do
                if (jorb_glob /= iorb_glob) then
                   do i=1,ndim
                      dpsir_i(i+(jorb-1)*ndim)=dpsir_i(i+(jorb-1)*ndim)+&
                           hfacj*rp_ij(i)*psir_i(i+(iorb-1)*ndim)
                   end do
                   !write(100+iproc,*)jorb+jsorb,iorb+isorb,igrpr(igroup) 
                end if
             else
                do i=1,ndim
                   dpsir_i(i+(iorb-1)*ndim)=dpsir_i(i+(iorb-1)*ndim)+&
                        hfaci*rp_ij(i)*psir_j(i+(jorb-1)*ndim)
                end do
             end if
          end if

          !fill the set of the vector to be sent to the other processes
          !in the first step the results are self-contained
          if (remote_result) then
             !write(100+iproc,*)jorb+jsorb,iorb+isorb,igrpr(igroup) 
             do i=1,ndim
                dpsir_j(i+(jorb-1)*ndim)=dpsir_j(i+(jorb-1)*ndim)+&
                     hfacj*rp_ij(i)*psir_i(i+(iorb-1)*ndim)
             end do
          end if
       end do
    end do

    call f_free(rp_ij)

  end subroutine internal_exctx_operation


  subroutine OP2P_exctx_run(iproc,nproc,ngrid,norbp,psir,dpsir,OP2P)
    use module_base
    use module_types
    implicit none
    integer, intent(in) :: iproc,nproc,ngrid,norbp
    type(OP2P_descriptors), intent(in) :: OP2P
    real(wp), dimension(ngrid,norbp), intent(in) :: psir
    real(wp), dimension(ngrid,norbp), intent(out) :: dpsir
    !local variables
    integer :: ierr

    !here any processor will initialise the global communications arrays needed for executing the op2p
    !the operation is symmetric
 
    call OP2P_communication(iproc,nproc,OP2P,psir,dpsir,internal_exctx_operation)

    call MPI_BARRIER(bigdft_mpi%mpi_comm,ierr)

  end subroutine OP2P_exctx_run

  subroutine orbs_to_attributes(iproc,nproc,orbs,OP2P,ndim)
    implicit none
    integer, intent(in) :: ndim,iproc,nproc
    type(orbitals_data), intent(in) :: orbs
    type(OP2P_descriptors), intent(out) :: OP2P
    !local variables
    character(len=*), parameter :: subname='orbs_to_attributes'
    integer :: igroup,iorb
    integer, dimension(:,:), allocatable :: orbs_attributes
    

    if (orbs%nkpts /=1) stop 'OP2P not allowed with K-points so far' 

    !build the array of the groups needed for OP2P module
    orbs_attributes = f_malloc((/ orbs%norb, 3 /),id='orbs_attributes')
    !the rule for the objects are listed here
    !objects_attributes(:,1) <= group to which the object belongs
    !objects_attributes(:,2) <= size in number of elements of the object
    !objects_attributes(:,3) <= size in number of elements of the results of the operations associated to the object

    !the number of groups is associated to the spin value
    !wavefunctions which belong to the same group has to be contiguously distributed on the processor
    !this might poses problems once generalising this technique to different k-points, since the k-points 
    !should be considered altogether
    !A lookup array should be defined, once the Poisson Solver updated
    do iorb=1,orbs%norb
       if (orbs%spinsgn(iorb) == 1.0_gp) then
          igroup=1
       else
          igroup=2
       end if
       orbs_attributes(iorb,1)=igroup
       orbs_attributes(iorb,2)=ndim !wavefunction
       orbs_attributes(iorb,3)=ndim !results
    end do
    
    !print *,'test',orbs%norb_par(:,0),orbs%norb,orbs_attributes

    !initialize the descriptors for the communications
    !the ExcatX operator is symmetric
    call initialize_OP2P_descriptors(.true.,iproc,nproc,orbs%norb,orbs_attributes,&
         orbs%norb_par(:,0),OP2P)

    !call the module for executing the exctX routine

    call f_free(orbs_attributes)

    if (orbs%nspin==2) then
       sfac=1.0_gp
    else 
       sfac=0.5_gp
    end if


  end subroutine orbs_to_attributes


end module module_exctx_op2p


!> Routine which applies the op2p module to calculate the exact exchange
!! Defines the interface module in the same file
subroutine exact_exchange_potential_op2p(iproc,nproc,xc,lr,orbs,pkernel,psi,dpsir,eexctX)
  use module_base
  use module_types
  use module_xc
  use module_exctx_op2p
  implicit none
  integer, intent(in) :: iproc,nproc
  type(xc_info), intent(in) :: xc
  type(locreg_descriptors), intent(in) :: lr
  type(orbitals_data), intent(in) :: orbs
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(in) :: psi !> wavefunctions in wavelet form
  type(coulomb_operator), intent(inout) :: pkernel !> Poisson Solver kernel, sequential in this case
  real(gp), intent(out) :: eexctX !> exact exchange energy
  real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,orbs%norbp), intent(out) :: dpsir !>Fock operator applied on the real-space wavefunctions
  !local variables
  character(len=*), parameter :: subname='exact_exchange_potential_op2p'
  integer :: iorb
  real(gp) :: exctXfac
  type(workarr_sumrho) :: w
  type(OP2P_descriptors) :: OP2P
  real(wp), dimension(:,:), allocatable :: psir

  call initialize_work_arrays_sumrho(1,lr,.true.,w)
  psir = f_malloc((/ lr%d%n1i*lr%d%n2i*lr%d%n3i, orbs%norbp /),id='psir')

  call f_zero(dpsir)

  !uncompress the wavefunction in the real grid
  do iorb=1,orbs%norbp
     !here ispinor is equal to one
     call daub_to_isf(lr,w,psi(1,1,iorb),psir(1,iorb))
  end do

  call deallocate_work_arrays_sumrho(w)

  !create the communicators descriptors
  call orbs_to_attributes(iproc,nproc,orbs,OP2P,lr%d%n1i*lr%d%n2i*lr%d%n3i)

  !pass the values to the interface module
  call OP2P_exctx_set(OP2P,orbs,pkernel)

  !run the calculation
  call OP2P_exctx_run(iproc,nproc,lr%d%n1i*lr%d%n2i*lr%d%n3i,orbs%norbp,psir,dpsir,OP2P)

  !recuperate the exctX energy and purge allocated variables
  call OP2P_exctx_clear(OP2P,eexctX)

  if (nproc>1) call mpiallred(eexctX,1,MPI_SUM,comm=bigdft_mpi%mpi_comm)

  exctXfac = xc_exctXfac(xc)
  eexctX=-exctXfac*eexctX

  call f_free(psir)

end subroutine exact_exchange_potential_op2p

