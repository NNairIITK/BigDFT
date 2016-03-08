!> @file
!!    Define routines for Lanczos diagonalization
!! @author
!!    Copyright (C) 2009-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Module defining the interfaces for routines which handle diagonalization
module lanczos_interface
  !use module_base
  use module_types
  use module_xc, only : xc_info
  use module_abscalc
  use communications_base, only: comms_cubic
  implicit none

  private

  !> Contains the arguments needed for the application of the hamiltonian
  type, public :: lanczos_args
     !arguments for the hamiltonian
     integer :: iproc,nproc,ndimpot,nspin, in_iat_absorber, Labsorber
     real(gp) :: hx,hy,hz
     type(energy_terms) :: energs
     !real(gp) :: ekin_sum,epot_sum,eexctX,eproj_sum,eSIC_DC
     type(atoms_data), pointer :: at
     type(orbitals_data), pointer :: orbs
     type(comms_cubic) :: comms
     type(DFT_PSP_projectors), pointer :: nlpsp
     type(local_zone_descriptors), pointer :: Lzd
     !type(gaussian_basis), pointer :: Gabsorber    !unused?
     type(SIC_data), pointer :: SIC
     type(xc_info), pointer :: xc
     integer, dimension(:,:), pointer :: ngatherarr 
     real(gp), dimension(:,:),  pointer :: rxyz
     !real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*orbs%norbp), pointer :: psi
     real(wp), dimension(:), pointer :: potential
     real(wp), dimension(:), pointer :: Gabs_coeffs
     !real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor*orbs%norbp) :: hpsi
     type(GPU_pointers), pointer :: GPU
     type(pcproj_data_type), pointer :: PPD
     type(pawproj_data_type), pointer :: PAWD
     ! removed from orbs, not sure if needed here or not
     !integer :: npsidim_orbs  !< Number of elements inside psi in the orbitals distribution scheme
     !integer :: npsidim_comp  !< Number of elements inside psi in the components distribution scheme
  end type lanczos_args

  !calculate the allocation dimensions
  public :: EP_inizializza,EP_initialize_start,EP_allocate_for_eigenprob,&
       &   get_EP_dim,set_EP_shift,&
       &   EP_mat_mult,EP_make_dummy_vectors, EP_normalizza,EP_copy,EP_scalare,&
       &   EP_copia_per_prova,&
       &   EP_set_all_random,EP_GramSchmidt,EP_add_from_vect_with_fact,&
       &   EP_Moltiplica, &
       &   EP_free,  EP_norma2_initialized_state ,EP_store_occupied_orbitals,EP_occprojections,&
       &   EP_multbyfact,EP_precondition,EP_Moltiplica4spectra,EP_ApplySinv,EP_ApplyS,&
       &   EP_scalare_multik,xabs_lanczos,xabs_cg,xabs_chebychev


  character(len=*), parameter :: subname='lanczos_interface'
  integer :: ierr
  type(lanczos_args), pointer :: ha

  real(kind=8), pointer :: matrix(:,:)
  real(wp), pointer :: Qvect   (:,:)
  real(wp), pointer :: dumQvect(:,:)
  real(wp), pointer :: occQvect(:)

  real(kind=8), pointer :: passed_matrix(:,:)

  real(kind=8) :: EP_shift
  integer :: EP_dim
  integer :: EP_dim_tot

  real(wp), dimension(:), pointer :: Qvect_tmp,wrk, wrk1, wrk2
  real(wp), dimension(:), pointer :: EP_norma2_initialized_state

  logical :: EP_doorthoocc
  integer :: EP_norb
  real(wp), dimension(:), pointer :: EP_occprojections

  real(wp), dimension(:,:,:,:,:,:), allocatable :: psi_gross
  logical , allocatable :: logrid(:,:,:)

contains

!!!  subroutine settapointer(p)
!!!    real(kind=8), intent(inout), TARGET :: p(:,:)
!!!    passed_matrix => p 
!!!    return
!!!  END SUBROUTINE settapointer

!!!  subroutine testapointer()
!!!    integer i,j
!!!    integer , pointer :: shape_a(:)
!!!    print *, " ecco il pointer " 
!!!    print *, " shape is " , shape(passed_matrix)
!!!    do i=1,  size(passed_matrix, 1)
!!!       print *, (passed_matrix(i,j), j=1, size(passed_matrix,2))
!!!    enddo
!!!    return
!!!  END SUBROUTINE testapointer



  integer function get_EP_dim()
    get_EP_dim=EP_dim_tot
    return
  END FUNCTION get_EP_dim


  subroutine set_EP_shift(shift)
    real(kind=8) shift
    EP_shift=shift
    return
  END SUBROUTINE set_EP_shift


  subroutine EP_inizializza(ha_actual)
    implicit none
    !Arguments
    type(lanczos_args), target :: ha_actual

    integer EP_dim_tot_touse

    ha=>ha_actual

    !! for the transposed representation
    EP_dim=sum( ha%comms%nvctr_par(ha%iproc,ha%orbs%iskpts+1:ha%orbs%iskpts+ ha%orbs%nkptsp )) *ha%orbs%nspinor

    !!   EP_dim_tot=(ha%Lzd%Glr%wfd%nvctr_c+7*ha%Lzd%Glr%wfd%nvctr_f)*ha%orbs%nspinor*ha%orbs%nkpts
    !! for the direct representation

    EP_dim_tot=(ha%Lzd%Glr%wfd%nvctr_c+7*ha%Lzd%Glr%wfd%nvctr_f)*ha%orbs%nspinor*ha%orbs%norbp

    EP_dim_tot_touse=max(EP_dim_tot,1)


!!$    if( (ha%Lzd%Glr%wfd%nvctr_c+7*ha%Lzd%Glr%wfd%nvctr_f)*ha%orbs%nkptsp /= &
!!$         sum(ha%comms%nvctr_par(:,1))  ) then
!!$       stop "array size inconsistency" 
!!$    endif

    if( (ha%Lzd%Glr%wfd%nvctr_c+7*ha%Lzd%Glr%wfd%nvctr_f) * ha%orbs%nkptsp  /= &
         &   sum(ha%comms%nvctr_par(:, ha%orbs%iskpts+1:ha%orbs%iskpts+ ha%orbs%nkptsp   ))  ) then
       stop "array size inconsistency" 
    endif
    !! arrays which are used in the direct rep after retrieval
    !! from transposed memory

    Qvect_tmp = f_malloc_ptr(EP_dim_tot_touse , id='Qvect_tmp')
    wrk = f_malloc_ptr(EP_dim_tot_touse,id='wrk')
    wrk1 = f_malloc_ptr(EP_dim_tot_touse,id='wrk1')
    wrk2 = f_malloc_ptr(EP_dim_tot_touse,id='wrk2')

    EP_shift=0.0

    EP_norb=0
    EP_doorthoocc=.false.

    EP_norma2_initialized_state = f_malloc_ptr(ha%orbs%norbp,id='EP_norma2_initialized_state')

    !added for nullification of the pointers in the lanczos_base module
    nullify(Qvect,dumQvect)

  END SUBROUTINE EP_inizializza


  subroutine EP_mat_mult( m,k ,  EV   ) ! usare i dumvectors a partire da -1
    implicit none
    !Arguments
    integer , intent(in):: m,k
    real(kind=8), intent (in) :: EV(1)

    call gemm('N','N', EP_dim , k, m  ,1.0_wp ,&
         &   Qvect(1,0)  , EP_dim ,&
         &   EV(1) ,m ,&
         &   0.0_wp ,  dumQvect(1,1) ,  EP_dim )

!!!
!!!    do i=0, k-1
!!!       do j=1, EP_dim
!!!          dumQvect(j,i+1)=0.0d0
!!!          do l=0, m-1 ! controllare il range
!!!             dumQvect(j,i+1)=dumQvect(j,i+1)+EV( l,i )*Qvect(  j,l ) 
!!!          enddo 
!!!       enddo
!!!    enddo

  END SUBROUTINE EP_mat_mult


  !> Allocate the wavefunctions in the transposed form, for lanczos
  subroutine EP_free(iproc)
    use yaml_output
    implicit none
    integer, intent(in) :: iproc

    if (iproc == 0)  call yaml_comment("DEALLOCATING")

    call f_free_ptr(Qvect)
    call f_free_ptr(dumQvect)
    call f_free_ptr(Qvect_tmp)
    call f_free_ptr(wrk)
    call f_free_ptr(wrk1)
    call f_free_ptr(wrk2)
    if(EP_doorthoocc) then
       call f_free_ptr(EP_occprojections)
    endif

    call f_free_ptr(EP_norma2_initialized_state)

  END SUBROUTINE EP_free


  subroutine  EP_allocate_for_eigenprob(nsteps)
    implicit none

    !Arguments
    integer, intent(in):: nsteps

    integer i, j

    if(associated(Qvect) ) then
       call f_free_ptr(Qvect)
    endif

    Qvect = f_malloc_ptr( (/1.to.EP_dim , 0.to.nsteps /),id='Qvect')
    do i=0,nsteps
       do j=1,EP_dim
          Qvect(j,i)=0.0_wp
       end do
    enddo
  END SUBROUTINE EP_allocate_for_eigenprob

  subroutine EP_make_dummy_vectors(nd)
    implicit none
    integer, intent(in):: nd

    if(associated(dumQvect) ) then
       call f_free_ptr(dumQvect)
    endif

    dumQvect = f_malloc_ptr( (/1.to.EP_dim , 0.to.nd+1/),id='dumQvect')
  END SUBROUTINE EP_make_dummy_vectors

  subroutine EP_store_occupied_orbitals( norb, Occ_psit )
    implicit none
    integer, intent(in) :: norb
    real(wp), dimension(:), pointer :: Occ_psit

    EP_doorthoocc=.true.
    EP_norb=norb
    occQvect=>Occ_psit
    EP_occprojections = f_malloc_ptr(EP_norb,id='EP_occprojections')

  END SUBROUTINE EP_store_occupied_orbitals

  subroutine EP_normalizza(j)
    implicit none
    integer, intent(in)::j
    ! ::::::::::::::::::::::::::::::::::::::
    if(J.ge.0) then
       call EP_normalizza_interno( Qvect(1:,j) )
    else
       call EP_normalizza_interno( dumQvect(1:,-j) )
    endif
    return 
  END SUBROUTINE EP_normalizza

  subroutine EP_multbyfact(j, fact)
    implicit none
    integer, intent(in)::j
    real(gp) fact
    ! ::::::::::::::::::::::::::::::::::::::
    if(J.ge.0) then
       call EP_multbyfact_interno( Qvect(1,j), fact )
    else
       call EP_multbyfact_interno( dumQvect(1,-j), fact )
    endif
    return 
  END SUBROUTINE EP_multbyfact

  subroutine EP_normalizza_interno(Q)
    implicit none
    !Arguments
    real(wp) :: Q(EP_dim)
    !Local variables
    real(wp) :: sump, sumtot

    sump = dot(EP_dim,Q(1),1 ,Q(1),1)
    sumtot=0

    if(ha%nproc/=1) then
       call MPI_Allreduce(sump,sumtot,1,mpidtypw, MPI_SUM,bigdft_mpi%mpi_comm ,ierr )
    else
       sumtot=sump
    endif
    sump=1.0_wp/sqrt(sumtot)

    call vscal(EP_dim,sump, Q(1), 1  ) 

  END SUBROUTINE EP_normalizza_interno

  subroutine EP_multbyfact_interno(Q, fact)
    implicit none
    !Arguments
    real(wp) :: Q
    real(gp) :: fact

    call vscal(EP_dim,fact, Q, 1  ) 

  END SUBROUTINE EP_multbyfact_interno

  subroutine EP_copy(i,j)
    implicit none
    !Arguments
    integer, intent(in) :: i,j

    if( i.ge.0 .and. j.ge.0) then
       call vcopy(EP_dim,Qvect(1,j),1,Qvect(1,i),1)
       !Qvect(:,i)=Qvect(:,j)
    else  if( i.lt.0 .and. j.ge.0) then
       call vcopy(EP_dim,Qvect(1,j),1,dumQvect(1,-i),1)
       !dumQvect(:,-i)=Qvect(:,j)
    else  if( i.ge.0 .and. j.lt.0) then
       call vcopy(EP_dim,dumQvect(1,-j),1,Qvect(1,i),1)
       !Qvect(:,i)=dumQvect(:,-j)
    else 
       call vcopy(EP_dim,dumQvect(1,-j),1,dumQvect(1,-i),1)
       !dumQvect(:,-i)=dumQvect(:,-j)
    endif

  END SUBROUTINE EP_copy

  subroutine   EP_scalare_multik(i,j, scalari)
    !use module_base
    implicit none
    !Arguments
    integer, intent(in) :: i,j
    real(wp) :: scalari(*)


    if( i.ge.0 .and. j.ge.0) then
       call EP_scalare_interna_multik(   Qvect(:,i)    ,   Qvect(:,j) , scalari )
    else  if( i.lt.0 .and. j.ge.0) then
       call EP_scalare_interna_multik(   dumQvect(:,-i),   Qvect(:,j) , scalari )
    else  if( i.ge.0 .and. j.lt.0) then
       call EP_scalare_interna_multik(   Qvect(:,i)    ,   dumQvect(:,-j) , scalari )
    else 
       call EP_scalare_interna_multik(   dumQvect(:,-i),   dumQvect(:,-j) , scalari )
    endif

  END SUBROUTINE EP_scalare_multik

  subroutine EP_scalare_interna_multik(a,b, scalari)
    implicit none
    !Arguments
    real(wp), intent(in):: a(*),b(*)
    real(wp) :: scalari(*)
    !Local variables
    real(wp) :: ovrlp_local ( ha%orbs%nkpts  )
    real(wp) :: ovrlp_global( ha%orbs%nkpts  )
    integer :: ik, ipos, ic, iorb
    call f_zero(ovrlp_local)
    call f_zero(ovrlp_global)
    if(ha%orbs%nspinor==1) then
       ipos=1
       do ik = 1, ha%orbs%nkptsp  !! this supposes norb=1 for chebychev
          do ic = 1, ha%comms%nvctr_par(ha%iproc,ha%orbs%iskpts+ik)
             ovrlp_local( ha%orbs%iskpts+ik)=ovrlp_local( ha%orbs%iskpts+ik)+a(ipos)*b(ipos)
             ipos=ipos+1
          end do
       end do
    else if (ha%orbs%nspinor ==2) then
       ipos=1
       do ik = 1, ha%orbs%nkptsp  !! this supposes norb=1 for chebychev
          do ic = 1, ha%comms%nvctr_par(ha%iproc,ha%orbs%iskpts+ik)
             ovrlp_local( ha%orbs%iskpts+ik)=ovrlp_local( ha%orbs%iskpts+ik)+a(ipos)*b(ipos)
             ovrlp_local( ha%orbs%iskpts+ik)=ovrlp_local( ha%orbs%iskpts+ik)+a(ipos+1)*b(ipos+1)
             !! for chebychev all scalar product should be real
             ipos=ipos+2
          end do
       end do
    else if (ha%orbs%nspinor ==4) then
       STOP "nspinor=4 not treated yet in EP_scalare_interna_multik"
    end if

    if(ha%nproc/=1) then
       call MPI_Allreduce(ovrlp_local,ovrlp_global,ha%orbs%nkpts ,mpidtypw, MPI_SUM,bigdft_mpi%mpi_comm ,ierr )
       do iorb = 1, ha%orbs%norbp !! this supposes norb=1 for chebychev
          scalari(iorb) = ovrlp_global(ha%orbs%iokpt(iorb ))
       end do
    else
       do ik=1, ha%orbs%norbp 
          scalari(ik)=ovrlp_local(ik)
       end do
    endif
  END SUBROUTINE EP_scalare_interna_multik





  real(kind=8) function   EP_scalare(i,j)
    implicit none
    !Arguments
    integer, intent(in) :: i,j

    if( i.ge.0 .and. j.ge.0) then
       EP_scalare = EP_scalare_interna(Qvect(1,i), Qvect(1,j) )
    else  if( i.lt.0 .and. j.ge.0) then
       EP_scalare = EP_scalare_interna(dumQvect(1,-i), Qvect(1,j) )
    else  if( i.ge.0 .and. j.lt.0) then
       EP_scalare = EP_scalare_interna(Qvect(1,i), dumQvect(1,-j) )
    else 
       EP_scalare = EP_scalare_interna(dumQvect(1,-i), dumQvect(1,-j) )
    endif
  END FUNCTION EP_scalare


  real(kind=8) function EP_scalare_interna(a,b)
    implicit none
    !Arguments
    real(kind=8), intent(in):: a,b
    !Local variables
    real(wp) :: sump, sumtot


    sump = dot(EP_dim,a,1 ,b,1)
    sumtot=0


    if(ha%nproc/=1) then
       call MPI_Allreduce(sump,sumtot,1,mpidtypw, MPI_SUM,bigdft_mpi%mpi_comm ,ierr )
    else
       sumtot=sump
    endif
    EP_scalare_interna=sumtot

  END FUNCTION EP_scalare_interna

  !  real(kind=8) function EP_scalare_interna(a,b)
  !    implicit none
  !    real(kind=8), intent(in):: a(EP_dim), b(EP_dim)
  !    ! ::::::::::::::::::::::::::::::::::::::::::::::
  !    integer i,j
  !
  !    real(wp) sump, sumtot
  !
  !    sump = dot(EP_dim,a(1),1 ,b(1),1)
  !    sumtot=0
  !
  !    if(ha%nproc/=1) then
  !       call MPI_Allreduce(sump,sumtot,1,mpidtypw, MPI_SUM,bigdft_mpi%mpi_comm ,ierr )
  !    else
  !       sumtot=sump
  !    endif
  !    EP_scalare_interna=sumtot
  !
  !  END FUNCTION EP_scalare_interna
  !

  subroutine  EP_copia_per_prova(psi)
    use communications, only: transpose_v
    !Arguments
    real(wp), dimension(*), target :: psi ! per testare happlication
    !Local variables
    integer :: i

    if( ha%nproc/=1) then
       call transpose_v(ha%iproc,ha%nproc,ha%orbs,ha%lzd%glr%wfd,ha%comms,&
            &   psi(1),wrk(1),out_add=Qvect(1,0))  
    else
       do i=1, EP_dim_tot
          Qvect(i,0)= psi(i)
       enddo
    endif

  END SUBROUTINE EP_copia_per_prova

  subroutine EP_initialize_start()
    use communications, only: transpose_v
    !Local variables
    integer :: i

    if ( ha%at%paw_NofL(ha%at%astruct%iatype(ha%in_iat_absorber )).gt.0  ) then
       if (ha%iproc == 0) write(*,*) "USING PTILDES TO BUILD INITIAL WAVE"
       !!if(EP_dim_tot.gt.0) then
       call f_zero(Qvect_tmp)
       call applyPAWprojectors(ha%orbs,ha%at,&
            &   ha%hx,ha%hy,ha%hz,ha%Lzd%Glr,ha%PAWD,Qvect_tmp,Qvect_tmp, ha%at%paw_S_matrices, &
            &   .true. ,    ha%in_iat_absorber, ha%Labsorber+1, &
            &   ha%Gabs_coeffs               ) 
       !!end if
    else
       STOP " ha%at%paw_NofL(ha%at%astruct%iatype(ha%in_iat_absorber )).gt.0  is false" 
!!$       Note G%psiat and G%xp have now 2 dimenstions.
!!$       call gaussians_to_wavelets_nonorm(ha%iproc,ha%nproc,ha%Lzd%Glr%geocode,ha%orbs,ha%Lzd%Glr%d,&
!!$            ha%hx,ha%hy,ha%hz,ha%Lzd%Glr%wfd,EP_Gabsorber,ha%Gabs_coeffs,Qvect_tmp )
    endif

!!$
!!$    inquire(file='DOWRITEDOWNINITIALORBITAL', exist=exist)
!!$    if(exist) then
!!$       if(  sum( ha%at%paw_NofL ).gt.0 ) then
!!$          if(  associated( ha%PAWD) ) then
!!$             call f_zero(EP_dim_tot  ,  wrk  )
!!$             call applyPAWprojectors(ha%orbs,ha%at,&
!!$                  ha%hx,ha%hy,ha%hz,ha%Lzd%Glr,ha%PAWD,Qvect_tmp,wrk, ha%at%paw_S_matrices, &
!!$                  .false.)      
!!$             do i=1, EP_dim_tot
!!$                Qvect_tmp(i) = Qvect_tmp(i) +wrk(i)
!!$             end do
!!$          endif
!!$       endif
!!$       
!!$       write(orbname,'(A)')'in_orb_'
!!$       call plot_wf_cube(orbname,ha%at,  ha%Lzd%Glr,ha%hx,ha%hy,ha%hz,ha%rxyz,Qvect_tmp ,"initial orbital" ) ! solo spinore 1!!$    endif

    if (.not. associated(Qvect)) then
       write(*,*)'ERROR: initialization vector not allocated!'
       stop
    end if

    if( ha%nproc/=1) then
       call transpose_v(ha%iproc,ha%nproc,ha%orbs,ha%lzd%glr%wfd,ha%comms,&
            &   Qvect_tmp(1),wrk(1),out_add=Qvect(1,0))
    else
       do i=1, EP_dim_tot
          Qvect(i,0)= Qvect_tmp(i)
       enddo
    endif



    call EP_scalare_multik(0,0, EP_norma2_initialized_state )


  END SUBROUTINE EP_initialize_start


  subroutine EP_set_all_random(i)
    implicit none
    !Arguments
    integer, intent(in) :: i
    !Local variables
    if ( i.ge.0 ) then
       call EP_set_random_interna(Qvect(1:,i))
    else 
       call EP_set_random_interna(dumQvect(1:,-i))
    endif

  END SUBROUTINE EP_set_all_random


  subroutine EP_set_random_interna(Q)
    implicit none
    !Arguments
    real(wp) :: Q(EP_dim)
    !Local variables
    real(4) :: rand
    integer :: i
    do i=1,EP_dim
       call random_number(rand)
       Q(i)=real(rand,wp)
    enddo
  END SUBROUTINE EP_set_random_interna


  subroutine EP_GramSchmidt(i,n)
    implicit none
    !Arguments
    integer, intent(in) :: i,n

    if( i.ge.0 ) then
       call EP_GramSchmidt_interna(Qvect(1:,i),n)
    else
       call EP_GramSchmidt_interna(dumQvect(1:,-i),n)
    endif

  END SUBROUTINE EP_GramSchmidt


  subroutine EP_GramSchmidt_interna(Q,n)
    implicit none
    !Arguments
    integer, intent(in) :: n
    real(wp) Q(EP_dim)
    !Local variables
    integer :: i, volta
    real(wp) :: scals(0:n-1), scalstot(0:n-1), occscals(1:EP_norb), occscalstot(1:EP_norb)

    if(ha%iproc.eq.0) then
       print *, " starting GramSchmidt ", n
    endif

    do volta=1,2
       call gemm('T','N', n , 1, EP_dim ,1.0_wp ,  Qvect(1,0)  , EP_dim ,Q(1) ,EP_dim , 0.0_wp , scals(0) ,  n )
       if(ha%nproc/=1) then
          call MPI_Allreduce(scals(0) ,scalstot(0) , n ,mpidtypw, MPI_SUM,bigdft_mpi%mpi_comm ,ierr )
       else
          do i=0,n-1
             scalstot(i)=scals(i)
          enddo
       endif
       do i=0, n-1
          call axpy(EP_dim, -scalstot(i)  ,   Qvect(1,i)  , 1,  Q(1) , 1)
       enddo

       if ( EP_doorthoocc) then
          call gemm('T','N', EP_norb , 1, EP_dim ,1.0_wp ,  occQvect(1)  , EP_dim ,Q(1) ,EP_dim , 0.0_wp , occscals(1) ,  EP_norb )
          if(ha%nproc/=1) then
             call MPI_Allreduce(occscals(1) ,occscalstot(1) , EP_norb  ,mpidtypw, MPI_SUM,bigdft_mpi%mpi_comm ,ierr )
          else
             do i=1,EP_norb
                occscalstot(i)=occscals(i)
             enddo
          endif
          do i=1, EP_norb
             call axpy(EP_dim, -occscalstot(i)  ,   occQvect(1+i*EP_dim)  , 1,  Q(1) , 1)
          enddo
       endif

    enddo

    if(ha%iproc.eq.0) then
       print *, "GramSchmidt done "
    endif

  END SUBROUTINE EP_GramSchmidt_interna


  !>   Hits the input array x with the kernel @f$((-1/2\Delta+C)_{ij})^{-1}@f$
  subroutine hit_with_kernel_spectra(x,z1,z3,kern_k1,kern_k2,kern_k3,n1,n2,n3,nd1,nd2,nd3,&
       &   n1f,n1b,n3f,n3b,nd1f,nd1b,nd3f,nd3b,ene, gamma )
    !use module_base
    implicit none
    ! Arguments
    integer, intent(in) :: n1,n2,n3,nd1,nd2,nd3
    integer, intent(in) :: n1f,n1b,n3f,n3b,nd1f,nd1b,nd3f,nd3b
    real(gp), intent(in) :: kern_k1(n1)
    real(gp), intent(in) :: kern_k2(n2)
    real(gp), intent(in) :: kern_k3(n3)
    real(kind=8), intent(in):: ene, gamma

    real(wp),intent(inout) :: x(n1,n2,n3)! input/output
    !Local variables
    real(wp) :: z1(2,nd1b,nd2,nd3,2)! work array
    real(wp) :: z3(2,nd1,nd2,nd3f,2)! work array
    real(gp) :: tt
    integer :: i1,i2,i3,inzee

    ! fft the input array x:

    call FFT_for(n1,n2,n3,n1f,n3f,nd1,nd2,nd3,nd1f,nd3f,x,z1,z3,inzee)

    ! hit the Fourier transform of x with the kernel. At the same time, transform the array
    ! from the form z3 (where only half of values of i3 are stored)
    ! to the form z1   (where only half of values of i1 are stored)
    ! The latter thing could be done separately by the subroutine z3_to_z1 that is contained
    ! in FFT_back, but then the code would be slower.

    !$omp parallel default (none) &
    !$omp private(i1,i2,i3,tt) &
    !$omp shared(z1,z3,kern_k1,kern_k2,kern_k3,n1b,n3f,inzee,n1,n2,n3,ene,gamma)

    ! i3=1: then z1 is contained in z3 
    !$omp do 
    do i2=1,n2
       do i1=1,n1b
          ! tt=1._gp/(kern_k1(i1)+kern_k2(i2)+kern_k3(1)+c)
          tt= 1._gp/((kern_k1(i1)+kern_k2(i2)+kern_k3(1)-ene)**2+gamma**2)
          z1(1,i1,i2,1,inzee)=z3(1,i1,i2,1,inzee)*tt
          z1(2,i1,i2,1,inzee)=z3(2,i1,i2,1,inzee)*tt
       enddo
    enddo
    !$omp enddo

    !$omp do
    do i3=2,n3f
       ! i2=1
       ! i1=1
       ! tt=1._gp/(kern_k1(1)+kern_k2(1)+kern_k3(i3)+c)
       tt= 1._gp/((kern_k1(1)+kern_k2(1)+kern_k3(i3)-ene)**2+gamma**2)
       z1(1,1,1,i3,inzee)=z3(1,1,1,i3,inzee)*tt
       z1(2,1,1,i3,inzee)=z3(2,1,1,i3,inzee)*tt

       z1(1,1,1,n3+2-i3,inzee)=z3(1,1,1,i3,inzee)*tt
       z1(2,1,1,n3+2-i3,inzee)=-z3(2,1,1,i3,inzee)*tt

       ! i2=1
       do i1=2,n1b  
          ! tt=1._gp/(kern_k1(i1)+kern_k2(1)+kern_k3(i3)+c)
          tt= 1._gp/((kern_k1(i1)+kern_k2(1)+kern_k3(i3)-ene)**2+gamma**2)
          z1(1,i1,1,i3,inzee)=z3(1,i1,1,i3,inzee)*tt
          z1(2,i1,1,i3,inzee)=z3(2,i1,1,i3,inzee)*tt

          z1(1,i1,1,n3+2-i3,inzee)= z3(1,n1+2-i1,1,i3,inzee)*tt
          z1(2,i1,1,n3+2-i3,inzee)=-z3(2,n1+2-i1,1,i3,inzee)*tt
       enddo

       do i2=2,n2
          ! i1=1
          ! tt=1._gp/(kern_k1(1)+kern_k2(i2)+kern_k3(i3)+c)
          tt= 1._gp/((kern_k1(1)+kern_k2(i2)+kern_k3(i3)-ene)**2+gamma**2)
          z1(1,1,i2,i3,inzee)=z3(1,1,i2,i3,inzee)*tt
          z1(2,1,i2,i3,inzee)=z3(2,1,i2,i3,inzee)*tt

          z1(1,1,i2,n3+2-i3,inzee)= z3(1,1,n2+2-i2,i3,inzee)*tt
          z1(2,1,i2,n3+2-i3,inzee)=-z3(2,1,n2+2-i2,i3,inzee)*tt

          do i1=2,n1b
             ! tt=1._gp/(kern_k1(i1)+kern_k2(i2)+kern_k3(i3)+c)
             tt= 1._gp/((kern_k1(i1)+kern_k2(i2)+kern_k3(i3)-ene)**2+gamma**2)
             z1(1,i1,i2,i3,inzee)=z3(1,i1,i2,i3,inzee)*tt
             z1(2,i1,i2,i3,inzee)=z3(2,i1,i2,i3,inzee)*tt

             z1(1,i1,i2,n3+2-i3,inzee)= z3(1,n1+2-i1,n2+2-i2,i3,inzee)*tt
             z1(2,i1,i2,n3+2-i3,inzee)=-z3(2,n1+2-i1,n2+2-i2,i3,inzee)*tt
          enddo
       enddo
    enddo
    !$omp enddo

    !$omp end parallel

    call FFT_back(n1,n2,n3,n1b,n3f,n3b,nd1,nd2,nd3,nd1b,nd3f,nd3b,x,z1,z3,inzee)

  END SUBROUTINE hit_with_kernel_spectra



  !> Multiplies a wavefunction psi_c,psi_f (in vector form) with a scaling vector (scal)
  subroutine wscal_f_spectra(mvctr_f,psi_f,hx,hy,hz,ene, gamma)
    ! multiplies a wavefunction psi_c,psi_f (in vector form) with a scaling vector (scal)
    !use module_base
    implicit none
    integer,intent(in)::mvctr_f
    real(gp),intent(in)::hx,hy,hz, ene, gamma

    real(wp)::psi_f(7,mvctr_f)
    real(gp)::scal(7),hh(3)
    !WAVELET AND SCALING FUNCTION SECOND DERIVATIVE FILTERS, diagonal elements
    real(gp),PARAMETER::B2=24.8758460293923314_gp,A2=3.55369228991319019_gp

    integer i

    hh(1)=.5_gp/hx**2
    hh(2)=.5_gp/hy**2
    hh(3)=.5_gp/hz**2

    scal(1)=1._gp/((b2*hh(1)+a2*hh(2)+a2*hh(3)-ene)**2+gamma*gamma)       !  2 1 1
    scal(2)=1._gp/((a2*hh(1)+b2*hh(2)+a2*hh(3)-ene)**2+gamma*gamma)       !  1 2 1
    scal(3)=1._gp/((b2*hh(1)+b2*hh(2)+a2*hh(3)-ene)**2+gamma*gamma)       !  2 2 1
    scal(4)=1._gp/((a2*hh(1)+a2*hh(2)+b2*hh(3)-ene)**2+gamma*gamma)       !  1 1 2
    scal(5)=1._gp/((b2*hh(1)+a2*hh(2)+b2*hh(3)-ene)**2+gamma*gamma)       !  2 1 2
    scal(6)=1._gp/((a2*hh(1)+b2*hh(2)+b2*hh(3)-ene)**2+gamma*gamma)       !  1 2 2
    scal(7)=1._gp/((b2*hh(1)+b2*hh(2)+b2*hh(3)-ene)**2+gamma*gamma)       !  2 2 2

    do i=1,mvctr_f
       psi_f(1,i)=psi_f(1,i)*scal(1)       !  2 1 1
       psi_f(2,i)=psi_f(2,i)*scal(2)       !  1 2 1
       psi_f(3,i)=psi_f(3,i)*scal(3)       !  2 2 1
       psi_f(4,i)=psi_f(4,i)*scal(4)       !  1 1 2
       psi_f(5,i)=psi_f(5,i)*scal(5)       !  2 1 2
       psi_f(6,i)=psi_f(6,i)*scal(6)       !  1 2 2
       psi_f(7,i)=psi_f(7,i)*scal(7)       !  2 2 2
    enddo

  END SUBROUTINE wscal_f_spectra




  !> Solves ((KE-ene)**2+gamma**2*I)*xx=yy by FFT in a cubic box 
  !! x_c is the right hand side on input and the solution on output
  !! This version uses work arrays kern_k1-kern_k3 and z allocated elsewhere
  subroutine prec_fft_fast_spectra(n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
       &   ene, gamma,hx,hy,hz,hpsi,&
       &   kern_k1,kern_k2,kern_k3,z1,z3,x_c,&
       &   nd1,nd2,nd3,n1f,n1b,n3f,n3b,nd1f,nd1b,nd3f,nd3b)
    !use module_base
    implicit none 
    integer, intent(in) :: n1,n2,n3
    integer,intent(in)::nd1,nd2,nd3
    integer,intent(in)::n1f,n3f,n1b,n3b,nd1f,nd3f,nd1b,nd3b
    integer, intent(in) :: nseg_c,nvctr_c,nseg_f,nvctr_f
    real(gp), intent(in) :: hx,hy,hz,ene, gamma
    integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
    integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
    ! real(wp), intent(inout) ::  hpsi(*) 
    real(wp), intent(inout) ::  hpsi(nvctr_c+7*nvctr_f) 

    !work arrays
    real(gp):: kern_k1(0:n1),kern_k2(0:n2),kern_k3(0:n3)
    real(wp),dimension(0:n1,0:n2,0:n3):: x_c! in and out of Fourier preconditioning
    real(wp)::z1(2,nd1b,nd2,nd3,2)! work array
    real(wp)::z3(2,nd1,nd2,nd3f,2)! work array 

    x_c=0.0_wp
    z1=0.0_wp
    z3=0.0_wp

    if (nvctr_f > 0) then
       call wscal_f_spectra(nvctr_f,hpsi(nvctr_c+1),hx,hy,hz,ene, gamma)
    end if

    call make_kernel(n1,hx,kern_k1)
    call make_kernel(n2,hy,kern_k2)
    call make_kernel(n3,hz,kern_k3)

    call uncompress_c(hpsi,x_c,keyg(1,1),keyv(1),nseg_c,nvctr_c,n1,n2,n3)

    call  hit_with_kernel_spectra(x_c,z1,z3,kern_k1,kern_k2,kern_k3,n1+1,n2+1,n3+1,nd1,nd2,nd3,&
         &   n1f,n1b,n3f,n3b,nd1f,nd1b,nd3f,nd3b,ene, gamma)

    call compress_c(hpsi,x_c,keyg(1,1),keyv(1),nseg_c,nvctr_c,n1,n2,n3)

  END SUBROUTINE prec_fft_fast_spectra


  subroutine EP_precondition(p,i, ene, gamma)
    !use module_base
    use communications, only: transpose_v, untranspose_v
    use locreg_operations
    !Arguments
    implicit none
    integer, intent(in) :: p,i
    real(gp) ene, gamma
    !Local variables
    integer :: k
!!$ real(wp), parameter :: b2=24.8758460293923314d0,a2=3.55369228991319019d0
    integer :: nd1,nd2,nd3,n1f,n3f,n1b,n3b,nd1f,nd3f,nd1b,nd3b 
    type(workarr_precond) :: w
    logical :: dopcproj
    dopcproj=.true.


    call allocate_work_arrays('P',.true.,1,ha%Lzd%Glr%d,w)

!!$
!!$    hh(1)=.5_wp/ha%hx**2
!!$    hh(2)=.5_wp/ha%hy**2
!!$    hh(3)=.5_wp/ha%hz**2
!!$    
!!$    scal(0)=1._wp/((a2*hh(1)+a2*hh(2)+a2*hh(3)-ene)**2+gamma*gamma)       !  1 1 1
!!$    scal(1)=1._wp/((b2*hh(1)+a2*hh(2)+a2*hh(3)-ene)**2+gamma*gamma)       !  2 1 1
!!$    scal(2)=1._wp/((a2*hh(1)+b2*hh(2)+a2*hh(3)-ene)**2+gamma*gamma)       !  1 2 1
!!$    scal(3)=1._wp/((b2*hh(1)+b2*hh(2)+a2*hh(3)-ene)**2+gamma*gamma)       !  2 2 1
!!$    scal(4)=1._wp/((a2*hh(1)+a2*hh(2)+b2*hh(3)-ene)**2+gamma*gamma)       !  1 1 2
!!$    scal(5)=1._wp/((b2*hh(1)+a2*hh(2)+b2*hh(3)-ene)**2+gamma*gamma)       !  2 1 2
!!$    scal(6)=1._wp/((a2*hh(1)+b2*hh(2)+b2*hh(3)-ene)**2+gamma*gamma)       !  1 2 2
!!$    scal(7)=1._wp/((b2*hh(1)+b2*hh(2)+b2*hh(3)-ene)**2+gamma*gamma)       !  2 2 2
!!$    

    call dimensions_fft(ha%Lzd%Glr%d%n1,ha%Lzd%Glr%d%n2,ha%Lzd%Glr%d%n3,&
         &   nd1,nd2,nd3,n1f,n3f,n1b,n3b,nd1f,nd3f,nd1b,nd3b)

    if( ha%nproc > 1) then
       if(i>=0) then
          call untranspose_v(ha%iproc,ha%nproc,ha%orbs,ha%Lzd%Glr%wfd,ha%comms,&
               &   Qvect(1,i), wrk(1),out_add=Qvect_tmp(1) )  
       else
          call untranspose_v(ha%iproc,ha%nproc,ha%orbs,ha%Lzd%Glr%wfd,ha%comms,&
               &   dumQvect(1,-i), wrk(1),out_add=Qvect_tmp(1) )  
       endif
    else
       if(i>=0) then
          do k=1, EP_dim_tot
             Qvect_tmp(k)=  Qvect(k,i)
          enddo
       else
          do k=1, EP_dim_tot
             Qvect_tmp(k)=  dumQvect(k,-i)
          enddo
       endif
    endif

!!$
!!$    do k=1,ha%lr%wfd%nvctr_c
!!$       wrk(k)=Qvect_tmp(k)*scal(0)           !  1 1 1
!!$    enddo
!!$    
!!$    do k=1,         ha%lr%wfd%nvctr_f
!!$       ind=ha%lr%wfd%nvctr_c+7*(k-1)
!!$       wrk(ind+1)=Qvect_tmp(ind+1)*scal(1)       !  2 1 1
!!$       wrk(ind+2)=Qvect_tmp(ind+2)*scal(2)       !  1 2 1
!!$       wrk(ind+3)=Qvect_tmp(ind+3)*scal(3)       !  2 2 1
!!$       wrk(ind+4)=Qvect_tmp(ind+4)*scal(4)       !  1 1 2
!!$       wrk(ind+5)=Qvect_tmp(ind+5)*scal(5)       !  2 1 2
!!$       wrk(ind+6)=Qvect_tmp(ind+6)*scal(6)       !  1 2 2
!!$       wrk(ind+7)=Qvect_tmp(ind+7)*scal(7)       !  2 2 2
!!$    enddo

    call vcopy(EP_dim_tot, Qvect_tmp(1),1,wrk(1),1) 


    if( dopcproj) then

       call f_zero(wrk1)
       ha%PPD%iproj_to_factor(:) =  1.0_gp
       call applyPCprojectors(ha%orbs,ha%at,ha%hx,ha%hy,ha%hz,&
            &   ha%Lzd%Glr,ha%PPD,wrk,wrk1 )


       call f_zero(wrk2)
!!$       do k=1, ha%PPD%mprojtot
!!$          print *, ha%PPD%iproj_to_ene(k), ha%PPD%iproj_to_l(k)
!!$       end do
       ha%PPD%iproj_to_factor = ( ha%PPD%iproj_to_ene - ene   )
       ha%PPD%iproj_to_factor =1.0_gp/(ha%PPD%iproj_to_factor**2 + gamma*gamma)

       call applyPCprojectors(ha%orbs,ha%at,ha%hx,ha%hy,ha%hz,&
            &   ha%Lzd%Glr,ha%PPD,wrk,wrk2 )
!!$
!!$

       call axpy(EP_dim_tot, -1.0_gp  ,  wrk1(1)   , 1,  wrk(1) , 1)

    endif


    call prec_fft_fast_spectra(ha%Lzd%Glr%d%n1,ha%Lzd%Glr%d%n2,ha%Lzd%Glr%d%n3,&
         &   ha%Lzd%Glr%wfd%nseg_c,ha%Lzd%Glr%wfd%nvctr_c,ha%Lzd%Glr%wfd%nseg_f,ha%Lzd%Glr%wfd%nvctr_f,&
         &   ha%Lzd%Glr%wfd%keygloc,ha%Lzd%Glr%wfd%keyvloc, &
         &   ene, gamma,ha%hx,ha%hy,ha%hz,wrk(1:),&
         &   w%kern_k1,w%kern_k2,w%kern_k3,w%z1,w%z3,w%x_c,&
         &   nd1,nd2,nd3,n1f,n1b,n3f,n3b,nd1f,nd1b,nd3f,nd3b)



    if( dopcproj) then

       call f_zero(wrk1)
       ha%PPD%iproj_to_factor(:) =  1.0_gp
       call applyPCprojectors(ha%orbs,ha%at,ha%hx,ha%hy,ha%hz,&
            &   ha%Lzd%Glr,ha%PPD,wrk,wrk1)


       call axpy(EP_dim_tot, -1.0_gp  ,  wrk1(1)   , 1,  wrk(1) , 1)
       call axpy(EP_dim_tot, +1.0_gp  ,  wrk2(1)   , 1,  wrk(1) , 1)
    endif


    if(  ha%iproc ==0 ) write(*,*)" done "



    if(p<0) then
       if(  ha%nproc/=1) then
          call transpose_v(ha%iproc,ha%nproc,ha%orbs,ha%lzd%glr%wfd,ha%comms,&
               &   wrk(1),Qvect_tmp(1),out_add=dumQvect(1,-p))  
       else
          do k=1, EP_dim_tot
             dumQvect(k,-p) =  wrk(k)
          enddo
       endif
    else
       if(  ha%nproc/=1) then
          call transpose_v(ha%iproc,ha%nproc,ha%orbs,ha%lzd%glr%wfd,ha%comms,&
               &   wrk(1),Qvect_tmp(1),out_add=Qvect(1,p))  
       else
          do k=1, EP_dim_tot
             Qvect(k,p) =  wrk(k)
          enddo
       endif
    endif

    call deallocate_work_arrays('P',.true.,1,w)

  END SUBROUTINE EP_precondition


  subroutine EP_Moltiplica4spectra(p,i, ene, gamma)
    use module_interfaces, only: FullHamiltonianApplication
    use gaussians, only: gaussian_basis
    use communications, only: transpose_v, untranspose_v
    use locreg_operations, only: confpot_data
    !Arguments
    implicit none
    integer, intent(in) :: p,i
    real(gp) :: ene, gamma
    !Local variables
    integer :: k
    type(confpot_data), dimension(ha%orbs%norbp) :: confdatarr
    type(paw_objects) :: paw

    if( ha%nproc > 1) then
       if(i>=0) then
          call untranspose_v(ha%iproc,ha%nproc,ha%orbs,ha%Lzd%Glr%wfd,ha%comms,&
               &   Qvect(1,i),wrk(1),out_add=Qvect_tmp(1) )  
       else
          call untranspose_v(ha%iproc,ha%nproc,ha%orbs,ha%Lzd%Glr%wfd,ha%comms,&
               &   dumQvect(1,-i),wrk(1),out_add=Qvect_tmp(1) )  
       endif
    else

       if(i>=0) then
          do k=1, EP_dim_tot
             Qvect_tmp(k)=  Qvect(k,i)
          enddo
       else
          do k=1, EP_dim_tot
             Qvect_tmp(k)=  dumQvect(k,-i)
          enddo
       endif
    endif

    call default_confinement_data(confdatarr,ha%orbs%norbp)

    paw%usepaw = .false.
    call FullHamiltonianApplication(ha%iproc,ha%nproc,ha%at,ha%orbs,&
         ha%Lzd,ha%nlpsp,confdatarr,ha%ngatherarr,ha%potential,Qvect_tmp,wrk,paw,&
         ha%energs,ha%SIC,ha%GPU,ha%xc)

    call axpy(EP_dim_tot, -ene  ,  Qvect_tmp(1)   , 1,  wrk(1) , 1)
    call vcopy(EP_dim_tot,wrk(1),1,Qvect_tmp(1),1)
    !Qvect_tmp   =  wrk

    call FullHamiltonianApplication(ha%iproc,ha%nproc,ha%at,ha%orbs,&
         ha%Lzd,ha%nlpsp,confdatarr,ha%ngatherarr,ha%potential,Qvect_tmp,wrk,paw,&
         ha%energs,ha%SIC,ha%GPU,ha%xc)

    call axpy(EP_dim_tot, -ene  ,  Qvect_tmp(1)   , 1,  wrk(1) , 1)

    if(  ha%iproc ==0 ) write(*,*)" done "

    if(p<0) then
       if(  ha%nproc/=1) then
          call transpose_v(ha%iproc,ha%nproc,ha%orbs,ha%lzd%glr%wfd,ha%comms,&
               &   wrk(1),Qvect_tmp(1),out_add=dumQvect(1,-p))  
       else
          do k=1, EP_dim_tot
             dumQvect(k,-p) =  wrk(k)
          enddo
       endif

       if(i>=0) then

          do k=1, EP_dim_tot
             dumQvect(k,-p)=Qvect(k,p)+gamma*gamma*Qvect(k,i)
          end do
       else
          do k=1, EP_dim_tot
             dumQvect(k,-p)=Qvect(k,p)+gamma*gamma*dumQvect(k,-i)
          end do
       endif

    else
       if(  ha%nproc/=1) then
          call transpose_v(ha%iproc,ha%nproc,ha%orbs,ha%lzd%glr%wfd,ha%comms,&
               &   wrk(1),Qvect_tmp(1),Qvect(1,p))  
       else
          do k=1, EP_dim_tot
             Qvect(k,p) =  wrk(k)
          enddo
       endif

       if(i>=0) then
          do k=1, EP_dim_tot
             Qvect(k,p)=Qvect(k,p)+gamma*gamma*Qvect(k,i)
          end do
       else
          do k=1, EP_dim_tot
             Qvect(k,p)=Qvect(k,p)+gamma*gamma*dumQvect(k,-i)
          end do
       endif

    endif

  END SUBROUTINE EP_Moltiplica4spectra


  subroutine EP_Moltiplica(p,i)
    use module_interfaces, only: FullHamiltonianApplication
    use communications, only: transpose_v, untranspose_v
    use locreg_operations, only: confpot_data
    !Arguments
    implicit none
    integer, intent(in) :: p,i
    !Local variables
    integer :: k
    type(confpot_data), dimension(ha%orbs%norbp) :: confdatarr
    type(paw_objects) :: paw

    if( ha%nproc > 1) then
       if(i>=0) then
          call untranspose_v(ha%iproc,ha%nproc,ha%orbs,ha%Lzd%Glr%wfd,ha%comms,&
               &   Qvect(1,i),wrk(1),out_add=Qvect_tmp(1) )  
       else
          call untranspose_v(ha%iproc,ha%nproc,ha%orbs,ha%Lzd%Glr%wfd,ha%comms,&
               &   dumQvect(1,-i),wrk(1),out_add=Qvect_tmp(1) )  
       endif
    else
       if(i>=0) then
          do k=1, EP_dim_tot
             Qvect_tmp(k)=  Qvect(k,i)
          enddo
       else
          do k=1, EP_dim_tot
             Qvect_tmp(k)=  dumQvect(k,-i)
          enddo
       endif
    endif

    if(ha%iproc==0 .and. .false. ) then
       !here the projector should be applied
       call lowpass_projector(ha%Lzd%Glr%d%n1,ha%Lzd%Glr%d%n2,ha%Lzd%Glr%d%n3,ha%Lzd%Glr%wfd%nvctr_c, Qvect_tmp )
    endif

    call default_confinement_data(confdatarr,ha%orbs%norbp)

    paw%usepaw = .false.
    call FullHamiltonianApplication(ha%iproc,ha%nproc,ha%at,ha%orbs,&
         ha%Lzd,ha%nlpsp,confdatarr,ha%ngatherarr,ha%potential,Qvect_tmp,wrk,paw,&
         ha%energs,ha%SIC,ha%GPU,ha%xc)

    if(  ha%iproc ==0 ) write(*,*)" done "


    if(  sum( ha%at%paw_NofL ).gt.0 ) then
       if(associated( ha%PAWD) ) then
          call applyPAWprojectors(ha%orbs,ha%at,&
               &   ha%hx,ha%hy,ha%hz,ha%Lzd%Glr,ha%PAWD,Qvect_tmp,wrk,ha%at%paw_H_matrices, .false.  )

!!$          call f_zero(EP_dim_tot  ,  wrk1  )
!!$          call applyPAWprojectors(ha%orbs,ha%at,&
!!$               ha%hx,ha%hy,ha%hz,ha%Lzd%Glr,ha%PAWD,wrk,wrk1, ha%at%paw_S_matrices )
!!$          do k=1, EP_dim_tot
!!$             wrk(k)  =   wrk(k)+  wrk1(k)
!!$          enddo

       endif
    endif



    if(ha%iproc==0 .and. .false.) then
       !here the projector should be applied
       print *,  "here the projector should be applied   ha%iproc   ",  ha%iproc 
       print *, ha%Lzd%Glr%d%n1,ha%Lzd%Glr%d%n2,ha%Lzd%Glr%d%n3
       call lowpass_projector(ha%Lzd%Glr%d%n1,ha%Lzd%Glr%d%n2,ha%Lzd%Glr%d%n3,ha%Lzd%Glr%wfd%nvctr_c,wrk)
       print *, " projector OK ha%iproc ",  ha%iproc  
    endif

    if(ha%iproc.eq.0) then
       if(EP_shift /= 0 ) then
          call axpy(EP_dim_tot, EP_shift  ,  Qvect_tmp(1)   , 1,  wrk(1) , 1)
       endif
    endif


    if(p<0) then
       if(  ha%nproc/=1) then
          call transpose_v(ha%iproc,ha%nproc,ha%orbs,ha%lzd%glr%wfd,ha%comms,&
               &   wrk(1),Qvect_tmp(1),out_add=dumQvect(1,-p))  
       else
          do k=1, EP_dim_tot
             dumQvect(k,-p) =  wrk(k)
          enddo
       endif
    else
       if(  ha%nproc/=1) then
          call transpose_v(ha%iproc,ha%nproc,ha%orbs,ha%lzd%glr%wfd,ha%comms,&
               &   wrk(1),Qvect_tmp(1),out_add=Qvect(1,p))  
       else
          do k=1, EP_dim_tot
             Qvect(k,p) =  wrk(k)
          enddo
       endif
    endif
    return 
  END SUBROUTINE EP_Moltiplica

  subroutine EP_ApplySinv(p,i)
    use communications, only: transpose_v, untranspose_v
    !Arguments
    implicit none
    integer, intent(in) :: p,i
    !Local variables
    integer :: k

    if( ha%nproc > 1) then
       if(i>=0) then
          call untranspose_v(ha%iproc,ha%nproc,ha%orbs,ha%Lzd%Glr%wfd,ha%comms,&
               &   Qvect(1,i),wrk(1),out_add=Qvect_tmp(1) )  
       else
          call untranspose_v(ha%iproc,ha%nproc,ha%orbs,ha%Lzd%Glr%wfd,ha%comms,&
               &   dumQvect(1,-i),wrk(1),out_add=Qvect_tmp(1) )  
       endif
    else
       if(i>=0) then
          do k=1, EP_dim_tot
             Qvect_tmp(k)=  Qvect(k,i)
          enddo
       else
          do k=1, EP_dim_tot
             Qvect_tmp(k)=  dumQvect(k,-i)
          enddo
       endif
    endif
    call f_zero(wrk)
    if(  sum( ha%at%paw_NofL ).gt.0 ) then
       if(  associated( ha%PAWD) ) then
          call applyPAWprojectors(ha%orbs,ha%at,&
               &   ha%hx,ha%hy,ha%hz,ha%Lzd%Glr,ha%PAWD,Qvect_tmp,wrk, ha%at%paw_Sm1_matrices, &
               &   .false.)      
       endif
    endif

    do k=1, EP_dim_tot
       wrk(k)=Qvect_tmp(k)+wrk(k)
    enddo



    if(p<0) then
       if(  ha%nproc/=1) then
          call transpose_v(ha%iproc,ha%nproc,ha%orbs,ha%lzd%glr%wfd,ha%comms,&
               &   wrk(1),Qvect_tmp(1),out_add=dumQvect(1,-p))  
       else
          do k=1, EP_dim_tot
             dumQvect(k,-p) =  wrk(k)
          enddo
       endif
    else
       if(  ha%nproc/=1) then
          call transpose_v(ha%iproc,ha%nproc,ha%orbs,ha%lzd%glr%wfd,ha%comms,&
               &   wrk(1),Qvect_tmp(1),out_add=Qvect(1,p))  
       else
          do k=1, EP_dim_tot
             Qvect(k,p) =  wrk(k)
          enddo
       endif
    endif
    return 
  END SUBROUTINE EP_ApplySinv




  subroutine EP_ApplyS(p,i)
    use communications, only: transpose_v, untranspose_v
    !Arguments
    implicit none
    integer, intent(in) :: p,i
    !Local variables
    integer :: k

    if( ha%nproc > 1) then
       if(i>=0) then
          call untranspose_v(ha%iproc,ha%nproc,ha%orbs,ha%Lzd%Glr%wfd,ha%comms,&
               &   Qvect(1,i),wrk(1),out_add=Qvect_tmp(1) )  
       else
          call untranspose_v(ha%iproc,ha%nproc,ha%orbs,ha%Lzd%Glr%wfd,ha%comms,&
               &   dumQvect(1,-i),wrk(1),out_add=Qvect_tmp(1) )  
       endif
    else
       if(i>=0) then
          do k=1, EP_dim_tot
             Qvect_tmp(k)=  Qvect(k,i)
          enddo
       else
          do k=1, EP_dim_tot
             Qvect_tmp(k)=  dumQvect(k,-i)
          enddo
       endif
    endif
    call f_zero(wrk)
    if(  sum( ha%at%paw_NofL ).gt.0 ) then
       if(  associated( ha%PAWD) ) then
          call applyPAWprojectors(ha%orbs,ha%at,&
               &   ha%hx,ha%hy,ha%hz,ha%Lzd%Glr,ha%PAWD,Qvect_tmp,wrk, ha%at%paw_S_matrices, &
               &   .false. )      
       endif
    endif

    do k=1, EP_dim_tot
       wrk(k)=Qvect_tmp(k)+wrk(k)
    enddo

    if(p<0) then
       if(  ha%nproc/=1) then
          call transpose_v(ha%iproc,ha%nproc,ha%orbs,ha%lzd%glr%wfd,ha%comms,&
               &   wrk(1),Qvect_tmp(1),out_add=dumQvect(1,-p))  
       else
          do k=1, EP_dim_tot
             dumQvect(k,-p) =  wrk(k)
          enddo
       endif
    else
       if(  ha%nproc/=1) then
          call transpose_v(ha%iproc,ha%nproc,ha%orbs,ha%lzd%glr%wfd,ha%comms,&
               &   wrk(1),Qvect_tmp(1),out_add=Qvect(1,p))  
       else
          do k=1, EP_dim_tot
             Qvect(k,p) =  wrk(k)
          enddo
       endif
    endif
    return 
  END SUBROUTINE EP_ApplyS



  subroutine EP_add_from_vect_with_fact(p,i, fact)
    implicit none
    !Arguments
    integer, intent(in) :: p,i
    real(kind=8), intent(in) :: fact

    if(i.ge.0  .and. p<0) then
       call axpy(EP_dim, fact,   Qvect(1,i)  , 1,  dumQvect(1,-p)  , 1)
    elseif(i.ge.0  .and. p>=0) then
       call axpy(EP_dim, fact,   Qvect(1,i)  , 1,  Qvect(1,p)  , 1)
    elseif(i<0  .and. p<0) then
       call axpy(EP_dim, fact,   dumQvect(1,-i)  , 1,  dumQvect(1,-p)  , 1)
    else
       call axpy(EP_dim, fact,   dumQvect(1,-i)  , 1,  Qvect(1, p)  , 1)
    endif

!!!    do k=1,EP_dim
!!!       dumQvect(k,-p) = dumQvect(k,-p)+fact*Qvect(k,i)
!!!    enddo

  END SUBROUTINE EP_add_from_vect_with_fact



  subroutine gaussians_to_wavelets_nonorm(iproc,nproc,geocode,orbs,grid,hx,hy,hz,wfd,G,wfn_gau,psi)
    !use module_base
    use module_types
    use gaussians
    implicit none
    character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
    integer, intent(in) :: iproc,nproc
    real(gp), intent(in) :: hx,hy,hz
    type(grid_dimensions), intent(in) :: grid
    type(wavefunctions_descriptors), intent(in) :: wfd
    type(orbitals_data), intent(in) :: orbs
    type(gaussian_basis), intent(in) :: G
    real(wp), dimension(G%ncoeff,orbs%nspinor,orbs%norbp), intent(in) :: wfn_gau
    real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(out) :: psi

    !local variables
    character(len=*), parameter :: subname='gaussians_to_wavelets'
    integer, parameter :: nterm_max=3
    logical :: maycalc
    integer :: ishell,iexpo,icoeff,iat,isat,ng,l,m,iorb,jorb,nterm,ierr,ispinor
    real(dp) :: normdev,tt,scpr,totnorm
    real(gp) :: rx,ry,rz
    integer, dimension(nterm_max) :: lx,ly,lz
    real(gp), dimension(nterm_max) :: fac_arr
    real(wp), dimension(:), allocatable :: tpsi

    if(iproc == 0 .and. verbose > 1) write(*,'(1x,a)',advance='no')'Writing wavefunctions in wavelet form '
    
    tpsi = f_malloc(wfd%nvctr_c+7*wfd%nvctr_f,id='tpsi')

    !initialize the wavefunction
    call f_zero(psi)
    !this can be changed to be passed only once to all the gaussian basis
    !eks=0.d0
    !loop over the atoms
    ishell=0
    iexpo=1
    icoeff=1



    do iat=1,G%nat
       rx=G%rxyz(1,iat)
       ry=G%rxyz(2,iat)
       rz=G%rxyz(3,iat)

       !loop over the number of shells of the atom type
       do isat=1,G%nshell(iat)
          ishell=ishell+1
          !the degree of contraction of the basis function
          !is the same as the ng value of the createAtomicOrbitals routine
          ng=G%ndoc(ishell)
          !angular momentum of the basis set(shifted for compatibility with BigDFT routines
          l=G%nam(ishell)
          !print *,iproc,iat,ishell,G%nam(ishell),G%nshell(iat)
          !multiply the values of the gaussian contraction times the orbital coefficient




          do m=1,2*l-1
             call calc_coeff_inguess(l,m,nterm_max,nterm,lx,ly,lz,fac_arr)
             !control whether the basis element may be
             !contribute to some of the orbital of the processor
             maycalc=.false.
             loop_calc: do iorb=1,orbs%norb
                if (orbs%isorb < iorb .and. iorb <= orbs%isorb+orbs%norbp) then
                   jorb=iorb-orbs%isorb
                   do ispinor=1,orbs%nspinor
                      if (wfn_gau(icoeff,ispinor,jorb) /= 0.0_wp) then
                         maycalc=.true.
                         exit loop_calc
                      end if
                   end do
                end if
             end do loop_calc
             if (maycalc) then
                call crtonewave(geocode,grid%n1,grid%n2,grid%n3,ng,nterm,lx,ly,lz,fac_arr,&
                     &   G%xp(1,iexpo),G%psiat(1,iexpo),&
                     &   rx,ry,rz,hx,hy,hz,&
                     &   0,grid%n1,0,grid%n2,0,grid%n3,&
                     &   grid%nfl1,grid%nfu1,grid%nfl2,grid%nfu2,grid%nfl3,grid%nfu3,  & 
                     wfd%nseg_c,wfd%nvctr_c,wfd%keygloc,wfd%keyvloc,wfd%nseg_f,wfd%nvctr_f,&
                     &   wfd%keygloc(1,wfd%nseg_c+1),wfd%keyvloc(wfd%nseg_c+1),&
                     &   tpsi(1),tpsi(wfd%nvctr_c+1))
             end if
             !sum the result inside the orbital wavefunction
             !loop over the orbitals
             do iorb=1,orbs%norb
                if (orbs%isorb < iorb .and. iorb <= orbs%isorb+orbs%norbp) then
                   jorb=iorb-orbs%isorb
                   do ispinor=1,orbs%nspinor
                      call axpy(wfd%nvctr_c+7*wfd%nvctr_f,wfn_gau(icoeff,ispinor,jorb),&
                           &   tpsi(1),1,psi(1,ispinor,jorb),1)
                   end do
                end if
             end do
             icoeff=icoeff+1
          end do
          iexpo=iexpo+ng
       end do
       if (iproc == 0 .and. verbose > 1) then
          write(*,'(a)',advance='no') &
               &   repeat('.',(iat*40)/G%nat-((iat-1)*40)/G%nat)
       end if
    end do

    if (iproc==0)    call gaudim_check(iexpo,icoeff,ishell,G%nexpo,G%ncoeff,G%nshltot)

    !if (iproc ==0  .and. verbose > 1) write(*,'(1x,a)')'done.'
    !renormalize the orbitals
    !calculate the deviation from 1 of the orbital norm
    normdev=0.0_dp
    tt=0.0_dp

    do iorb=1,orbs%norb
       if (orbs%isorb < iorb .and. iorb <= orbs%isorb+orbs%norbp) then
          jorb=iorb-orbs%isorb
          totnorm=0.0_dp
          do ispinor=1,orbs%nspinor !to be verified in case of nspinor=4
             call wnrm(wfd%nvctr_c,wfd%nvctr_f,psi(1,ispinor,jorb),psi(wfd%nvctr_c+1,ispinor,jorb),scpr) 
             totnorm=totnorm+scpr
          end do
          !! print *, " TOTNORM QUI " , totnorm 
          totnorm=1.0_wp
          do ispinor=1,orbs%nspinor !to be verified in case of nspinor=4
             call wscal(wfd%nvctr_c,wfd%nvctr_f,real(1.0_dp/sqrt(totnorm),wp),psi(1,ispinor,jorb),psi(wfd%nvctr_c+1,ispinor,jorb))
          end do
          !write(*,'(1x,a,i5,1pe14.7)')'norm of orbital ',iorb,totnorm
          tt=max(tt,abs(1.0_dp-totnorm))
       end if
    end do
    if (nproc > 1) then
       call MPI_REDUCE(tt,normdev,1,mpidtypd,MPI_MAX,0,bigdft_mpi%mpi_comm,ierr)
    else
       normdev=tt
    end if
    if (iproc ==0) write(*,'(1x,a,1pe12.2)')&
         &   'Deviation from normalization of the imported orbitals',normdev

    call f_free(tpsi)

  END SUBROUTINE gaussians_to_wavelets_nonorm


  subroutine lowpass_projector(n1,n2,n3,nvctr,psi)
    !use module_base
    implicit none
    integer, intent(in) :: nvctr,n1,n2,n3
    real(wp), dimension(nvctr), intent(inout) :: psi
    !local variables
    character(len=*), parameter :: subname='projector'
    real(gp) rx,ry,rz
    real(gp) rxc,ryc,rzc
    integer ia, iat
    logical isgross
    real(gp) radius_gross
    integer ix,iy,iz
    real(gp) r
    integer countgross



    !test whether the coarse grid fills the whole box
    if (nvctr /= (n1+1)*(n2+1)*(n3+1)) then
       write(*,*)' ERROR: projector nont implemented for non-filled coarse grids'
       stop
    end if

    !for the moment the dimensions should be odd
    if (mod(n1,2) /=1 .or. mod(n2,2) /=1 .or. mod(n3,2) /=1) then
       write(*,*)' ERROR: the dimensions should be odd'
       if (ha%nproc > 1) call MPI_FINALIZE(ierr)
       stop
    end if


    if(.not. allocated(psi_gross)) then 
       psi_gross = f_malloc((/ 0.to.n1/2, 1.to.2, 0.to.n2/2, 1.to.2, 0.to.n3/2, 1.to.2 /),id='psi_gross')
       logrid = f_malloc((/ 0.to.n1/2, 0.to.n2/2, 0.to.n3/2 /),id='logrid')
    endif

    call analyse_per_self(n1/2,n2/2,n3/2,psi,psi_gross)

!!$
!!$
!!$  ! coarse grid quantities
!!$  call fill_logrid(ha%at%astruct%geocode,n1/2,n2/2,n3/2,0,n1/2,0,n2/2,0,n3/2,0,ha%at%astruct%nat,&
!!$       ha%at%astruct%ntypes,ha%at%astruct%iatype,ha%rxyz,ha%radii_cf(1,1),8.0_gp,2.0_gp*ha%hx,2.0_gp*ha%hy,2.0_gp*ha%hz,logrid)

!!$
    do ia =1, ha%at%astruct%nat
       iat=ha%at%astruct%iatype(ia)
       radius_gross = ha%at%radii_cf(  iat ,1 )*8
       if(ha%iproc==0) print *, " for atomo " , ia , " fine degrees of freedo,g will be thrown beyond radius =  " ,  radius_gross
    enddo
    countgross=0
    do ix=1,n1/2
       do iy=1,n2/2
          do iz = 1,n3/2
             rx = ha%hx*(ix-1)           *2.0  
             ry = ha%hy*(iy-1)           *2.0  
             rz = ha%hz*(iz-1)           *2.0  
             isgross=.true.
             do ia =1, ha%at%astruct%nat
                iat=ha%at%astruct%iatype(ia)
                radius_gross = ha%at%radii_cf(  iat ,1 )*8
                rxc = rx  -  ha%rxyz(1,ia )
                ryc = ry  -  ha%rxyz(2,ia )
                rzc = rz  -  ha%rxyz(3,ia )
                r  = sqrt( rxc*rxc+ryc*ryc+rzc*rzc)
                if( r < radius_gross) isgross=.false.
                if(.not. isgross) exit
             enddo
             if(isgross) then
                countgross=countgross+1
                psi_gross(ix,2 ,iy,1 , iz,1    )=0.0_wp
!!$
                psi_gross(ix,2 ,iy,2 , iz,2    )=0.0_wp
                psi_gross(ix,1 ,iy,2 , iz,2    )=0.0_wp
                psi_gross(ix,2 ,iy,1 , iz,2    )=0.0_wp
                psi_gross(ix,2 ,iy,2 , iz,1    )=0.0_wp
                psi_gross(ix,1 ,iy,1 , iz,2    )=0.0_wp
                psi_gross(ix,1 ,iy,2 , iz,1    )=0.0_wp
!!$
             endif
          enddo
       enddo
    enddo


!!$   do iz=1,n3/2
!!$      do iy=1,n2/2
!!$         do ix = 1,n1/2
!!$            if(.not. logrid(ix,iy,iz)) then
!!$               psi_gross(ix,2 ,iy,1 , iz,1    )=0.0_wp
!!$            endif
!!$         enddo
!!$      end do
!!$      do iy=1,n2/2
!!$         do ix = 1,n1/2
!!$            if(.not. logrid(ix,iy,iz)) then
!!$               psi_gross(ix,1 ,iy,2 , iz,1    )=0.0_wp
!!$               psi_gross(ix,2 ,iy,2 , iz,1    )=0.0_wp
!!$            endif
!!$         enddo
!!$      end do
!!$   end do
!!$   do iz=1,n3/2
!!$      do iy=1,n2/2
!!$         do ix = 1,n1/2
!!$            if(.not. logrid(ix,iy,iz)) then
!!$               psi_gross(ix,1 ,iy,1 , iz,2    )=0.0_wp
!!$               psi_gross(ix,2 ,iy,1 , iz,2    )=0.0_wp
!!$            endif
!!$         enddo
!!$      end do
!!$      do iy=1,n2/2
!!$         do ix = 1,n1/2
!!$            if(.not. logrid(ix,iy,iz)) then
!!$               psi_gross(ix,1 ,iy,2 , iz,2    )=0.0_wp
!!$               psi_gross(ix,2 ,iy,2 , iz,2    )=0.0_wp
!!$            endif
!!$         enddo
!!$      end do
!!$   enddo


    !here the projector operation

!!$   psi_gross=0.0_wp

    call synthese_per_self(n1/2,n2/2,n3/2,psi_gross,psi)

  END SUBROUTINE lowpass_projector

  !> Lanczos diagonalization
  subroutine xabs_lanczos(iproc,nproc,at,hx,hy,hz,rxyz,&
       nlpsp,Lzd,dpcom,potential,&
       energs,xc,nspin,GPU,in_iat_absorber,&
       in , PAWD , orbs )! add to interface
    !use module_base
    use module_dpbox, only: denspot_distribution
    use module_types
    use module_xc
    use lanczos_base
    use module_interfaces, only: free_full_potential
    use communications_init, only: orbitals_communicators
    use communications_base, only: deallocate_comms
    use rhopotential, only: full_local_potential
    implicit none
    integer, intent(in) :: iproc,nproc,nspin
    real(gp), intent(in) :: hx,hy,hz
    type(atoms_data), intent(in), target :: at
    type(DFT_PSP_projectors), intent(in), target :: nlpsp
    type(local_zone_descriptors), intent(inout), target :: Lzd
    type(denspot_distribution), intent(in), target :: dpcom
    type(xc_info), intent(in), target :: xc
    real(gp), dimension(3,at%astruct%nat), intent(in), target :: rxyz
    real(wp), dimension(max(dpcom%ndimpot,1),nspin), target :: potential
    type(energy_terms), intent(inout) :: energs
    type(GPU_pointers), intent(inout) , target :: GPU
    integer, intent(in) :: in_iat_absorber

    type(input_variables),intent(in), target :: in
    type(pawproj_data_type), target ::PAWD
    type(orbitals_data), intent(inout), target :: orbs

    !local variables
    character(len=*), parameter :: subname='lanczos'
    type(lanczos_args) :: ha
    integer :: i

    real(wp), pointer :: Gabs_coeffs(:)
    real(wp), dimension(:), pointer :: pot

    if(iproc==0) write(*,*) " IN ROUTINE LANCZOS "

    GPU%full_locham=.true.
    if (GPU%OCLconv) then
       call allocate_data_OCL(Lzd%Glr%d%n1,Lzd%Glr%d%n2,Lzd%Glr%d%n3,at%astruct%geocode,&
            &   in%nspin,Lzd%Glr%wfd,orbs,GPU)
       if (iproc == 0) write(*,*) 'GPU data allocated'
    end if


    orbs%eval = f_malloc_ptr(orbs%norb,id='orbs%eval')
    orbs%occup(1:orbs%norb)=1.0_gp
    orbs%spinsgn(1:orbs%norb)=1.0_gp
    orbs%eval(1:orbs%norb)=1.0_gp
    !call allocate_comms(nproc,ha%comms,subname)
    call orbitals_communicators(iproc,nproc,Lzd%Glr,orbs,ha%comms)  

    call local_potential_dimensions(iproc,Lzd,orbs,xc,dpcom%ngatherarr(0,1))

    Gabs_coeffs = f_malloc_ptr(2*in%L_absorber+1,id='Gabs_coeffs')


    if(   at%paw_NofL( at%astruct%iatype(   in_iat_absorber ) ) .gt. 0   ) then     
       Gabs_coeffs(:)=in%Gabs_coeffs(:)
    else     
       print * ," You are asking for a spectra for atom " , in_iat_absorber
       print *, " but at%paw_NofL( at%astruct%iatype(   in_iat_absorber ) )=0 " 
       print *, " this mean that the pseudopotential file is not pawpatched. "
       print *, " You'll have to generated the patch with pseudo"
       STOP     
    endif

    call full_local_potential(iproc,nproc,orbs,Lzd,0,dpcom,xc,potential,pot)

!!$   call full_local_potential(iproc,nproc,ndimpot,Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i,&
!!$        in%nspin,Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*in%nspin,0,&
!!$        orbs,Lzd,0,ngatherarr,potential,pot)

    ha%in_iat_absorber=in_iat_absorber
    ha%Labsorber  = in%L_absorber
    ha%iproc=iproc
    ha%nproc=nproc
    ha%at=>at !!
    ha%hx=hx 
    ha%hy=hy
    ha%hz=hz
    ha%rxyz=>rxyz

    ha%nlpsp=>nlpsp !!
    ha%Lzd=>Lzd !!!
    ha%ngatherarr=>dpcom%ngatherarr
    ha%ndimpot=dpcom%ndimpot
    ha%potential=>pot
    ha%energs = energs
    ha%nspin=nspin
    ha%GPU=>GPU !!
    ha%Gabs_coeffs=>Gabs_coeffs
    ha%PAWD=> PAWD
    ha%SIC=>in%SIC
    ha%xc=>xc
    ha%orbs=>orbs

    call EP_inizializza(ha) 

    if(.true.) then
       LB_nsteps =in%nsteps
       call LB_allocate_for_lanczos( )
       call EP_allocate_for_eigenprob(LB_nsteps)
       call EP_make_dummy_vectors(10)


       call LB_passeggia(0,LB_nsteps,      get_EP_dim, EP_initialize_start , EP_normalizza,&
            &   EP_Moltiplica, EP_GramSchmidt ,EP_set_all_random, EP_copy,   EP_mat_mult, &
            &   EP_scalare,EP_add_from_vect_with_fact     )


       if(ha%iproc==0) then
          open(unit=22,file="alphabeta")
          write(22,*) LB_nsteps, EP_norma2_initialized_state
          print *, " alpha and beta from lanczos "
          WRITE(*,'(I5,1000(1ES23.16))')    LB_nsteps, EP_norma2_initialized_state
          do i=0, LB_nsteps-1
             write(22,*) LB_alpha(i), LB_beta(i)
             WRITE(*,'(2ES23.16)')  LB_alpha(i), LB_beta(i)
          enddo

          close(unit=22)
       endif

       call LB_de_allocate_for_lanczos( )

    endif

    call deallocate_comms(ha%comms)

    call EP_free(ha%iproc)

    if (GPU%OCLconv) then
       call free_gpu_OCL(GPU,orbs,in%nspin)
    end if

    call f_free_ptr(orbs%eval)

    call f_free_ptr(Gabs_coeffs)

    call free_full_potential(dpcom%mpi_env%nproc,0,xc,pot)

  END SUBROUTINE xabs_lanczos


  !> Chebychev polynomials to calculate the density of states
  subroutine xabs_chebychev(iproc,nproc,at,hx,hy,hz,rxyz,&
       nlpsp,Lzd,dpcom,potential,&
       energs,xc,nspin,GPU,in_iat_absorber,in, PAWD , orbs  )

    !use module_base
    use module_dpbox, only: denspot_distribution
    use module_types
    use module_xc
    use lanczos_base
    ! per togliere il bug 
    use module_interfaces, only: free_full_potential
    use communications_init, only: orbitals_communicators
    use communications_base, only: deallocate_comms
    use rhopotential, only: full_local_potential

    implicit none
    integer  :: iproc,nproc,nspin
    real(gp)  :: hx,hy,hz
    type(atoms_data), target :: at
    type(DFT_PSP_projectors), target :: nlpsp
    type(local_zone_descriptors), target :: Lzd
    type(denspot_distribution), intent(in), target :: dpcom
    type(xc_info), intent(in), target :: xc
    real(gp), dimension(3,at%astruct%nat), target :: rxyz
    real(wp), dimension(max(dpcom%ndimpot,1),nspin), target :: potential
    type(energy_terms), intent(inout) :: energs
    type(GPU_pointers), intent(inout) , target :: GPU
    integer, intent(in) :: in_iat_absorber 

    type(input_variables),intent(in), target :: in
    type(pawproj_data_type), target ::PAWD
    type(orbitals_data), intent(inout), target :: orbs

    !Local variables
    character(len=*), parameter :: subname='chebychev'
    type(lanczos_args) :: ha
    integer :: i

    real(wp), dimension(:), pointer  :: pot

    real(gp) :: eval_min, eval_max, fact_cheb, cheb_shift
!    real(gp) :: Pi
    logical:: dopaw
    real(gp) :: GetBottom

    if (iproc==0) print *, " IN ROUTINE  chebychev  "

!    Pi=acos(-1.0_gp)

    GPU%full_locham=.true.

    if (GPU%OCLconv) then
       call allocate_data_OCL(Lzd%Glr%d%n1,Lzd%Glr%d%n2,Lzd%Glr%d%n3,at%astruct%geocode,&
            &   in%nspin,Lzd%Glr%wfd,orbs,GPU)
       if (iproc == 0) write(*,*)&
            &   'GPU data allocated'
    end if

    orbs%eval = f_malloc_ptr(orbs%norb*orbs%nkpts,id='orbs%eval')
    orbs%occup(1:orbs%norb*orbs%nkpts )=1.0_gp
    orbs%spinsgn(1:orbs%norb*orbs%nkpts )=1.0_gp
    orbs%eval(1:orbs%norb*orbs%nkpts )=1.0_gp

    call orbitals_communicators(iproc,nproc,Lzd%Glr,orbs,ha%comms)  

    call local_potential_dimensions(iproc,Lzd,orbs,xc,dpcom%ngatherarr(0,1))

    if(   at%paw_NofL( at%astruct%iatype(   in_iat_absorber ) ) .gt. 0   ) then     
    else
       print * ," You are asking for a spactra for atom " , in_iat_absorber
       print *, " but at%paw_NofL( at%astruct%iatype(   in_iat_absorber ) )=0 " 
       print *, " this mean that the pseudopotential file is not pawpatched. "
       print *, " You'll have to generated the patch with pseudo"
       STOP     
    endif

    call full_local_potential(iproc,nproc,orbs,Lzd,0,dpcom,xc,potential,pot)
!!$   call full_local_potential(iproc,nproc,ndimpot,Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i,&
!!$        in%nspin,Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*in%nspin,0,&
!!$        orbs,Lzd,0,ngatherarr,potential,pot)

    !associate hamapp_arg pointers
    ha%in_iat_absorber=in_iat_absorber
    ha%Labsorber  = in%L_absorber
    ha%iproc=iproc
    ha%nproc=nproc
    ha%at=>at !!
    ha%hx=hx 
    ha%hy=hy
    ha%hz=hz
    ha%rxyz=>rxyz
    ha%nlpsp=>nlpsp !!
    ha%Lzd=>Lzd !!!
    ha%ngatherarr=>dpcom%ngatherarr
    ha%ndimpot=dpcom%ndimpot
    ha%potential=>pot
    ha%energs = energs
    ha%nspin=nspin
    ha%GPU=>GPU !!
    ha%Gabs_coeffs=>in%Gabs_coeffs
    ha%PAWD=> PAWD
    ha%SIC=>in%SIC
    ha%orbs=>orbs
    ha%xc=>xc

    call EP_inizializza(ha)  

!!$  if(.false.) then
!!$
!!$     ! trova il valore massimo 
!!$     shift =-0.0
!!$     tol   =1.0D-8
!!$     accontentati_di=1
!!$     
!!$     cercacount = LB_cerca( 10, shift, tol, set_EP_shift, EP_allocate_for_eigenprob, EP_make_dummy_vectors, &
!!$          get_EP_dim, EP_initialize_start , EP_normalizza, EP_Moltiplica, EP_GramSchmidt, &
!!$          EP_set_all_random, EP_copy , EP_mat_mult,  EP_scalare,EP_add_from_vect_with_fact,accontentati_di)
!!$     
!!$     if(iproc==0) then
!!$        print *, " maximal eigenvalues " 
!!$        print *, LB_eval
!!$     endif
!!$     eval_max = LB_eval(0)
!!$     
!!$     ! trova il valore minimo 
!!$     shift =-10000
!!$     
!!$     
!!$     cercacount = LB_cerca( 10, shift, tol, set_EP_shift, EP_allocate_for_eigenprob, EP_make_dummy_vectors, &
!!$          get_EP_dim, EP_initialize_start , EP_normalizza, EP_Moltiplica, EP_GramSchmidt, &
!!$          EP_set_all_random, EP_copy , EP_mat_mult,  EP_scalare,EP_add_from_vect_with_fact,accontentati_di)
!!$     
!!$     if(iproc==0) then
!!$        print *, " minima eigenvalues" 
!!$        print *, LB_eval
!!$     endif
!!$     eval_min = LB_eval(0)+10000
!!$
!!$  else

    eval_min = GetBottom(  at, nspin)-1.0 - in%abscalc_bottomshift
    eval_max = 4.0*Pi*Pi*(1.0/hx/hx + 1.0/hy/hy + 1.0/hz/hz  )/2.0*1.1 +2

    !   endif

    cheb_shift=0.5*(eval_min+ eval_max) 
    fact_cheb = (2-0.0001)/(eval_max-eval_min)

    LB_nsteps = in%nsteps
    LB_norbp  = orbs%norbp
    LB_nproc=nproc
    LB_iproc=iproc

    call LB_allocate_for_chebychev( )
    call EP_allocate_for_eigenprob(6) 
    call EP_make_dummy_vectors(4)


    !! call set_EP_shift(-cheb_shift) 
    call set_EP_shift(0.0_gp) 

    if( sum( ha%at%paw_NofL ).gt.0 ) then
       dopaw=.true.
    else
       dopaw=.false.
    endif
    if(ha%iproc==0) then
       print *, "weigths ", orbs%kwgts
    endif

    call LB_passeggia_Chebychev (LB_nsteps, cheb_shift,  fact_cheb,     get_EP_dim, EP_initialize_start , EP_normalizza,&
         &   EP_Moltiplica, EP_GramSchmidt ,EP_set_all_random, EP_copy,   EP_mat_mult, &
         &   EP_scalare_multik,EP_add_from_vect_with_fact  , EP_multbyfact, EP_ApplySinv, EP_ApplyS, dopaw, &
         &   in%abscalc_S_do_cg,  in%abscalc_Sinv_do_cg, in%xabs_res_prefix, orbs%nkpts , orbs%norb_par, &
         &   orbs%kwgts(   orbs%iskpts+1    :    orbs%iskpts + orbs%norb_par(ha%iproc,0))   )

    if(ha%iproc==0) then
       print *, "coefficients from Chebychev "
       WRITE(*,'(I5,2ES23.16)')   2*LB_nsteps, cheb_shift,  fact_cheb
       print *,"... " 
       do i=0, 2*LB_nsteps-1
          if(i>2*LB_nsteps-1 -10) then
             WRITE(*,'(1ES23.16)')   LB_alpha_cheb(1, i)
          endif
       enddo
    endif

    call free_full_potential(dpcom%mpi_env%nproc,0,xc,pot)
    nullify(ha%potential)


    !deallocate communication and orbitals descriptors
    call deallocate_comms(ha%comms)

    call EP_free(ha%iproc)
    call  LB_de_allocate_for_cheb( )

!!$ this free is already executed by bigdft
!!$
    if (GPU%OCLconv) then
       call free_gpu_OCL(GPU,orbs,in%nspin)
    end if

    call f_free_ptr(orbs%eval)


!!$  i_all=-product(shape(Gabsorber%nshell))*kind(Gabsorber%nshell)
!!$  deallocate(Gabsorber%nshell,stat=i_stat)
!!$  call memocc(i_stat,i_all,'Gabsorber%nshell',subname)
!!$
!!$  i_all=-product(shape(Gabsorber%nam))*kind(Gabsorber%nam)
!!$  deallocate(Gabsorber%nam,stat=i_stat)
!!$  call memocc(i_stat,i_all,'Gabsorber%nam',subname)
!!$
!!$  i_all=-product(shape(Gabsorber%ndoc))*kind(Gabsorber%ndoc)
!!$  deallocate(Gabsorber%ndoc,stat=i_stat)
!!$  call memocc(i_stat,i_all,'Gabsorber%ndoc',subname)
!!$
!!$  i_all=-product(shape(Gabsorber%xp))*kind(Gabsorber%xp)
!!$  deallocate(Gabsorber%xp,stat=i_stat)
!!$  call memocc(i_stat,i_all,'Gabsorber%xp',subname)
!!$
!!$  i_all=-product(shape(Gabsorber%psiat))*kind(Gabsorber%psiat)
!!$  deallocate(Gabsorber%psiat,stat=i_stat)
!!$  call memocc(i_stat,i_all,'Gabsorber%psiat',subname)
!!$
!!$  if( associated(Gabsorber%rxyz)) then
!!$     i_all=-product(shape(Gabsorber%rxyz))*kind(Gabsorber%rxyz)
!!$     deallocate(Gabsorber%rxyz,stat=i_stat)
!!$     call memocc(i_stat,i_all,'Gabsorber%rxyz',subname)
!!$  endif

  END SUBROUTINE xabs_chebychev


  !> Finds the spectra solving  (H-omega)x=b
  subroutine xabs_cg(iproc,nproc,at,hx,hy,hz,rxyz,&
       nlpsp,Lzd,dpcom,potential,&
       energs,xc,nspin,GPU,in_iat_absorber,&
       in , rhoXanes, PAWD , PPD, orbs )
    !use module_base
    use module_dpbox, only: denspot_distribution
    use module_types
    use lanczos_base
    use module_xc
    use communications_base, only: deallocate_comms
    ! per togliere il bug 
    use module_interfaces, only: free_full_potential
    use communications_init, only: orbitals_communicators
    use rhopotential, only: full_local_potential

    implicit none

    integer  :: iproc,nproc,nspin
    real(gp)  :: hx,hy,hz
    type(atoms_data), target :: at
    type(DFT_PSP_projectors), target :: nlpsp
    type(local_zone_descriptors), target :: Lzd
    type(pcproj_data_type), target ::PPD
    type(denspot_distribution), intent(in), target :: dpcom
    type(xc_info), intent(in), target :: xc
    real(gp), dimension(3,at%astruct%nat), target :: rxyz
    real(wp), dimension(max(dpcom%ndimpot,1),nspin), target :: potential
    real(wp), dimension(max(dpcom%ndimpot,1),nspin), target :: rhoXanes
    type(energy_terms), intent(inout) :: energs
    type(GPU_pointers), intent(inout) , target :: GPU
    integer, intent(in) :: in_iat_absorber
    type(pawproj_data_type), target ::PAWD
    type(input_variables),intent(in), target :: in
    type(orbitals_data), intent(inout), target :: orbs

    !local variables
    character(len=*), parameter :: subname='xabs_cg'
    type(lanczos_args) :: ha
    integer :: i,j
    real(gp) Ene,gamma,  res 

    real(wp),   pointer  :: Gabs_coeffs(:)
    real(wp), dimension(:), pointer  :: pot

    logical:: useold
    real(gp) , pointer ::potentialclone(:,:)

    if( iand( in%potshortcut,16)>0) then
       potentialclone = f_malloc_ptr((/ max(dpcom%ndimpot, 1), nspin /),id='potentialclone')
       potentialclone=potential
    endif

    if(iproc==0) print *, " IN ROUTINE xabs_cg "

    GPU%full_locham=.true.
    if (GPU%OCLconv) then
       call allocate_data_OCL(Lzd%Glr%d%n1,Lzd%Glr%d%n2,Lzd%Glr%d%n3,at%astruct%geocode,&
            &   in%nspin,Lzd%Glr%wfd,orbs,GPU)
       if (iproc == 0) write(*,*)&
            &   'GPU data allocated'
    end if

    orbs%eval = f_malloc_ptr(orbs%norb,id='orbs%eval')

    orbs%occup(1:orbs%norb)=1.0_gp
    orbs%spinsgn(1:orbs%norb)=1.0_gp
    orbs%eval(1:orbs%norb)=1.0_gp
    !call allocate_comms(nproc,ha%comms,subname)
    call orbitals_communicators(iproc,nproc,Lzd%Glr,orbs,ha%comms)  

    call local_potential_dimensions(iproc,Lzd,orbs,xc,dpcom%ngatherarr(0,1))

    Gabs_coeffs = f_malloc_ptr(2*in%L_absorber+1,id='Gabs_coeffs')

    if(   at%paw_NofL( at%astruct%iatype(   in_iat_absorber ) ) .gt. 0   ) then     
       Gabs_coeffs(:)=in%Gabs_coeffs(:)
    else
       print * ," You are asking for a spactra for atom " , in_iat_absorber
       print *, " but at%paw_NofL( at%astruct%iatype(   in_iat_absorber ) )=0 " 
       print *, " this mean that the pseudopotential file is not pawpatched. "
       print *, " You'll have to generated the patch with pseudo"
       STOP     
    endif

    call full_local_potential(iproc,nproc,orbs,Lzd,0,dpcom,xc,potential,pot)
!!$   call full_local_potential(iproc,nproc,ndimpot,Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i,&
!!$        in%nspin,Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*in%nspin,0,&
!!$        orbs,Lzd,0,ngatherarr,potential,pot)

    ha%in_iat_absorber=in_iat_absorber
    ha%Labsorber  = in%L_absorber
    ha%iproc=iproc
    ha%nproc=nproc
    ha%at=>at !!
    ha%hx=hx 
    ha%hy=hy
    ha%hz=hz
    ha%rxyz=>rxyz

    ha%nlpsp=>nlpsp !!
    ha%Lzd=>Lzd !!!
    ha%ngatherarr=>dpcom%ngatherarr
    ha%ndimpot=dpcom%ndimpot
    ha%potential=>pot
    ha%energs = energs
    ha%nspin=nspin
    ha%GPU=>GPU !!
    ha%Gabs_coeffs=>Gabs_coeffs
    ha%PAWD=> PAWD 
    ha%PPD=> PPD
    ha%SIC=>in%SIC
    ha%xc=>xc
    ha%orbs=>orbs

    call EP_inizializza(ha) 

    if(.true.) then
       LB_nsteps =in%nsteps
       call LB_allocate_for_lanczos( )
       call EP_allocate_for_eigenprob(10)
       call EP_make_dummy_vectors(10)

       do i=0,0

          ene = 0.22_gp + i*0.03_gp

          if( iand( in%potshortcut,16)>0) then
             potential=potentialclone
             do j=1, dpcom%ndimpot

                if( mod(j-1,100)==0) then
                   print *, " dirac_hara punto",j
                endif
                call dirac_hara (rhoXanes(j,1), ene , potential(j,1))
             enddo
          endif

          gamma = 0.03_gp

          if(i==0) then
             useold=.false.
          else
             useold=.true.
          endif

          res =  LB_cg(    get_EP_dim, EP_initialize_start , EP_normalizza,&
               &   EP_Moltiplica4spectra,  EP_copy,  &
               &   EP_scalare,EP_add_from_vect_with_fact   , EP_multbyfact  ,EP_precondition, Ene, gamma, 1.0D-2, useold )

          print *, ene, res
          open(unit=22,file="cgspectra.dat", position="append")
          write(22,*) ene, res
          close(unit=22)

       enddo

       call LB_de_allocate_for_lanczos( )

    endif

    call deallocate_comms(ha%comms)

    call EP_free(ha%iproc)

    if (GPU%OCLconv) then
       call free_gpu_OCL(GPU,orbs,in%nspin)
    end if

    call f_free_ptr(orbs%eval)

    call f_free_ptr(Gabs_coeffs)

    if( iand( in%potshortcut,16)>0) then
       call f_free_ptr(potentialclone)
    endif

    call free_full_potential(dpcom%mpi_env%nproc,0,xc,pot)
    nullify(ha%potential)

  END SUBROUTINE xabs_cg

END MODULE lanczos_interface
