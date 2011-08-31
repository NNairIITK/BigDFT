!> Routine which applies the op2p module to calculate the exact exchange
!! Defines the interface module in the same file
subroutine exact_exchange_potential_op2p(iproc,nproc,geocode,nspin,lr,orbs,&
     hxh,hyh,hzh,pkernel,psi,dpsir,eexctX)
  use module_base
  use module_types
  use module_xc
  use module_exctx_op2p
  implicit none
  character(len=1), intent(in) :: geocode !> the boundary conditions
  integer, intent(in) :: iproc,nproc,nspin 
  real(gp), intent(in) :: hxh,hyh,hzh
  type(locreg_descriptors), intent(in) :: lr
  type(orbitals_data), intent(in) :: orbs
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(in) :: psi !> wavefunctions in wavelet form
  real(dp), dimension(*), intent(in) :: pkernel !> Poisson Solver kernel, sequential in this case
  real(gp), intent(out) :: eexctX !> exact exchange energy
  real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,orbs%norbp), intent(out) :: dpsir !>Fock operator applied on the real-space wavefunctions
  !local variables
  character(len=*), parameter :: subname='exact_exchange_potential_op2p'
  integer :: i_stat,i_all
  type(workarr_sumrho) :: w

  call initialize_work_arrays_sumrho(lr,w)
  allocate(psir(lr%d%n1i*lr%d%n2i*lr%d%n3i,orbs%norbp+ndebug),stat=i_stat)
  call memocc(i_stat,psir,'psir',subname)

  call razero(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%norbp,psir)

  !uncompress the wavefunction in the real grid
  do iorb=1,orbs%norbp
     !here ispinor is equal to one
     call daub_to_isf(lr,w,psi(1,1,iorb),psir(1,iorb))
  end do

  call deallocate_work_arrays_sumrho(w)

  !build the array of the groups needed for OP2P module
  allocate(orbs_attributes(norb,3+ndebug),stat=i_stat)
  call memocc(i_stat,orbs_attributes,'orbs_attributes',subname)
  !the rule for the objects are listed here
  !objects_attributes(:,1) <= group to which the object belongs
  !objects_attributes(:,2) <= size in number of elements of the object
  !objects_attributes(:,3) <= size in number of elements of the results of the operations associated to the object
  
  !the number of groups is associated to the spin value
  !wavefunctions which belong to the same group has to be contiguously distributed on the processor
  !this might poses problems once generalising this technique to different k-points.
  !A lookup array should be defined, once the Poisson Solver updated
  do iorb=1,norb
     if (orbs%spinsgn(iorb) == 1.0_gp) then
        igroup=1
     else
        igroup=2
     end if
     objects_attributes(iorb,1)=igroup
     objects_attributes(iorb,2)=lr%d%n1i*lr%d%n2i*lr%d%n3i !wavefunction
     objects_attributes(iorb,3)=lr%d%n1i*lr%d%n2i*lr%d%n3i !results
  end do



  i_all=-product(shape(orbs_attributes))*kind(orbs_attributes)
  deallocate(orbs_attributes,stat=i_stat)
  call memocc(i_stat,i_all,'orbs_attributes',subname)

  i_all=-product(shape(psir))*kind(psir)
  deallocate(psir,stat=i_stat)
  call memocc(i_stat,i_all,'psir',subname)

end subroutine exact_exchange_potential_op2p

module module_exctx_op2p
  use module_base
  use module_types
  use op2p_module
  implicit none

  private

  public :: op2p_exctx_binding

  character(len=1) :: geocode
  real(gp) :: sfac,hfac,hxh,hyh,hzh,eexctX
  type(locreg_descriptors), pointer :: lr
  type(orbitals_data), pointer :: orbs
  real(dp), dimension(:), pointer :: pkernel

contains

  subroutine internal_exctx_operation(istep,iproc,igroup,remote_result,&
     isorb,jsorb,iorbs,jorbs,norbi,norbj,&
     nvctri,nvctrj,nvctri_results,nvctrj_results,&
     psir_i,psir_j,&
     dpsir_i,dpsir_j)
    use module_base
    use module_types
    use Poisson_Solver
    implicit none
    logical, intent(in) :: remote_result
    integer, intent(in) :: istep,iproc,igroup,isorb,jsorb,iorbs,jorbs,norbi,norbj
    integer, intent(in) :: nvctri,nvctrj,nvctri_results,nvctrj_results
    real(wp), dimension(nvctri), intent(in) :: psir_i
    real(wp), dimension(nvctrj), intent(in) :: psir_j 
    real(wp), dimension(nvctri_results), intent(inout) :: dpsir_i
    real(wp), dimension(nvctrj_results), intent(inout) :: dpsir_j
    !local variables
    integer :: iorb,jorb,i_all,i_stat,iorb_glob,jorb_glob
    real(gp) :: hfaci,hfacj,hfac2,ehart
    real(wp), dimension(:), allocatable :: rp_ij

    allocate(rp_ij(lr%d%n1i*lr%d%n2i*lr%d%n3i+ndebug),stat=i_stat)
    call memocc(i_stat,rp_ij,'rp_ij',subname)

    !for the first step do only the upper triangular part
    do iorb=1,norbi
       iorb_glob=iorb+iorbs+isorb
       hfacj=-sfac*orbs%occup(iorb+isorb)
       do jorb=1,norbj
          jorb_glob=jorb+jorbs+jsorb
          !first cross-check whether the spin indices are the same
          if (orbs%spinsgn(iorb_glob) /= orbs%spinsgn(jorb_glob)) then
             write(*,*)'ERROR in partitioning the orbitals',&
                  iorb_glob,jorb_glob,igroup,jsorb,iproc
             stop
          end if
          hfaci=-sfac*orbs%occup(jorb_glob)
          !do it only for upper triangular results
          if (istep /= 0 .or. jorb_glob >= iorb_glob) then
             if (istep == 0 ) then
                do i=1,lr%d%n1i*lr%d%n2i*lr%d%n3i
                   rp_ij(i)=hfac*psir_i(i,iorb)*psir_i(i,jorb)
                end do
             else
                do i=1,lr%d%n1i*lr%d%n2i*lr%d%n3i
                   rp_ij(i)=hfac*psir_i(i,iorb)*psir_j(i,jorb)
                end do
             end if
             ncalls=ncalls+1                    
             !Poisson solver in sequential
             if (iproc == OP2P%iprocref .and. verbose > 1) then
                write(*,'(1x,a,i3,a2)')'Exact exchange calculation: ',&
                     nint(real(ncalls,gp)/real(ncalltot,gp)*100.0_gp),' %'
             end if

             call H_potential(geocode,'D',0,1,&
                  lr%d%n1i,lr%d%n2i,lr%d%n3i,hxh,hyh,hzh,&
                  rp_ij,pkernel,rp_ij,ehart,0.0_dp,.false.,&
                  quiet='YES')

             !this factor is only valid with one k-point
             !can be easily generalised to the k-point case
             hfac2=sfac*orbs%occup(iorb_glob)*orbs%occup(jorb_glob)

             !exact exchange energy
             if (iorb_glob == jorb_glob) then
                eexctX=eexctX+hfac2*real(ehart,gp)
             else
                !if the result has to be sent away
                if (remote_results) then
                   eexctX=eexctX+2.0_gp*hfac2*real(ehart,gp)
                else !otherwise other processors are already calculating it
                   eexctX=eexctX+hfac2*real(ehart,gp)
                end if
             end if
             !accumulate the results for each of the wavefunctions concerned
             if (istep == 0) then
                do i=1,lr%d%n1i*lr%d%n2i*lr%d%n3i
                   dpsir_i(i,iorb)=dpsir_i(i,iorb)+&
                        hfaci*rp_ij(i)*psir_i(i,jorb)
                end do
                if (jorb_glob /= iorb_glob) then
                   do i=1,lr%d%n1i*lr%d%n2i*lr%d%n3i
                      dpsir_i(i,jorb)=dpsir_i(i,jorb)+&
                           hfacj*rp_ij(i)*psir_i(i,iorb)
                   end do
                   !write(100+iproc,*)jorb+jsorb,iorb+isorb,igrpr(igroup) 
                end if
             else
                do i=1,lr%d%n1i*lr%d%n2i*lr%d%n3i
                   dpsir_i(i,iorb)=dpsir_i(i,iorb)+&
                        hfaci*rp_ij(i)*psir_j(i,jorb)
                end do
             end if
          end if

          !fill the set of the vector to be sent to the other processes
          !in the first step the results are self-contained
          if (remote_results) then
             !write(100+iproc,*)jorb+jsorb,iorb+isorb,igrpr(igroup) 
             do i=1,lr%d%n1i*lr%d%n2i*lr%d%n3i
                dpsir_j(i,jorb)=dpsir_j(i,jorb)+&
                     hfacj*rp_ij(i)*psir_i(i,iorb)
             end do
          end if
       end do
    end do

    i_all=-product(shape(rp_ij))*kind(rp_ij)
    deallocate(rp_ij,stat=i_stat)
    call memocc(i_stat,i_all,'rp_ij',subname)

  end subroutine internal_exctx_operation


  subroutine op2p_exctx_binding(iproc,nproc,norb,orbs_attributes,norb_par
    use module_base
    use module_types
    implicit none


    if (nspin==2) then
       sfac=1.0_gp
       ngroup=2
    else 
       sfac=0.5_gp
       ngroup=1
    end if

    hfac=1.0_gp/(hxh*hyh*hzh)

    !here any processor will initialise the global communications arrays needed for executing the op2p
    !the operation is symmetric
    call initialize_OP2P_descriptors(.true.,iproc,nproc,norb,orbs_attributes,norb_par,OP2P)

    call OP2P_communication(iproc,nproc,OP2P,psir,dpsir,fake_operation)

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    !print *,'iproc,icount2',iproc,icount
    call free_OP2P_descriptors(OP2P,subname)

  end subroutine op2p_exctx_binding

  subroutine bind_op2p_eexctx
    use module_base
    use module_types
    implicit none
    
  end subroutine bind_op2p_eexctx

end module module_exctx_op2p