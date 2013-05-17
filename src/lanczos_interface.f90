!> @file
!!    Define routines for Lanczos diagonalization
!! @author
!!    Copyright (C) 2009-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Interface for routines which handle diagonalization
module lanczos_interface
   use module_base
   use module_types
   implicit none

   private

   !calculate the allocation dimensions
   public :: EP_inizializza,EP_initialize_start,EP_allocate_for_eigenprob,&
      &   get_EP_dim,set_EP_shift,&
      &   EP_mat_mult,EP_make_dummy_vectors, EP_normalizza,EP_copy,EP_scalare,&
      &   EP_copia_per_prova,&
      &   EP_set_all_random,EP_GramSchmidt,EP_add_from_vect_with_fact,&
      &   EP_Moltiplica, &
      &   EP_free,  EP_norma2_initialized_state ,EP_store_occupied_orbitals,EP_occprojections,&
      &   EP_multbyfact,EP_precondition,EP_Moltiplica4spectra,EP_ApplySinv,EP_ApplyS,&
      &   EP_scalare_multik


   character(len=*), parameter :: subname='lanczos_interface'
   integer :: i_all,i_stat,ierr
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

allocate(Qvect_tmp( EP_dim_tot_touse  +ndebug) , stat=i_stat) 
call memocc(i_stat,Qvect_tmp,'Qvect_tmp',subname)

allocate(wrk  ( EP_dim_tot_touse   +ndebug) , stat=i_stat )
call memocc(i_stat,wrk,'wrk',subname)
allocate(wrk1 (  EP_dim_tot_touse    +ndebug) , stat=i_stat )
call memocc(i_stat,wrk1,'wrk1',subname)
allocate(wrk2 (  EP_dim_tot_touse    +ndebug) , stat=i_stat )
call memocc(i_stat,wrk2,'wrk2',subname)

EP_shift=0.0

EP_norb=0
EP_doorthoocc=.false.

allocate( EP_norma2_initialized_state( ha%orbs%norbp+ndebug ) )
call memocc(i_stat,EP_norma2_initialized_state,'EP_norma2_initialized_state' , subname)

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

     i_all=-product(shape(Qvect))*kind(Qvect)
     deallocate(Qvect,stat=i_stat)
     call memocc(i_stat,i_all,'Qvect',subname)

     i_all=-product(shape(dumQvect))*kind(dumQvect)
     deallocate(dumQvect,stat=i_stat)
     call memocc(i_stat,i_all,'dumQvect',subname)

     i_all=-product(shape(Qvect_tmp))*kind(Qvect_tmp)
     deallocate(Qvect_tmp,stat=i_stat)
     call memocc(i_stat,i_all,'Qvect_tmp',subname)

     i_all=-product(shape(wrk))*kind(wrk)
     deallocate(wrk,stat=i_stat)
     call memocc(i_stat,i_all,'wrk',subname)

     i_all=-product(shape(wrk1))*kind(wrk1)
     deallocate(wrk1,stat=i_stat)
     call memocc(i_stat,i_all,'wrk1',subname)

     i_all=-product(shape(wrk2))*kind(wrk2)
     deallocate(wrk2,stat=i_stat)
     call memocc(i_stat,i_all,'wrk2',subname)

     if(EP_doorthoocc) then
        i_all=-product(shape(EP_occprojections))*kind(EP_occprojections)
        deallocate(EP_occprojections,stat=i_stat)
        call memocc(i_stat,i_all,'EP_occprojections',subname)
     endif

     i_all=-product(shape(EP_norma2_initialized_state))*kind(EP_norma2_initialized_state)
     deallocate(EP_norma2_initialized_state,stat=i_stat)
     call memocc(i_stat,i_all,'EP_norma2_initialized_state',subname)

  END SUBROUTINE EP_free


  subroutine  EP_allocate_for_eigenprob(nsteps)
     implicit none

     !Arguments
     integer, intent(in):: nsteps

     integer i, j

     if(associated(Qvect) ) then
        i_all=-product(shape(Qvect))*kind(Qvect)

        deallocate(Qvect,stat=i_stat)
        call memocc(i_stat,i_all,'Qvect',subname)
     endif

     allocate(Qvect(EP_dim,0:nsteps+ndebug),&
        &   stat=i_stat)
     call memocc(i_stat,Qvect,'Qvect',subname)
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
        i_all=-product(shape(dumQvect))*kind(dumQvect)
        deallocate(dumQvect,stat=i_stat)
        call memocc(i_stat,i_all,'dumQvect',subname)
     endif

     allocate(dumQvect(EP_dim,0:nd+1+ndebug),&
        &   stat=i_stat)
     call memocc(i_stat,dumQvect,'dumQvect',subname)
  END SUBROUTINE EP_make_dummy_vectors

  subroutine EP_store_occupied_orbitals( norb, Occ_psit )
     implicit none
     integer, intent(in) :: norb
     real(wp), dimension(:), pointer :: Occ_psit

     EP_doorthoocc=.true.
     EP_norb=norb
     occQvect=>Occ_psit
     allocate(EP_occprojections(EP_norb+ndebug) , stat=i_stat )
     call memocc(i_stat,EP_occprojections,'EP_occprojections',subname)

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
        call dcopy(EP_dim,Qvect(1,j),1,Qvect(1,i),1)
        !Qvect(:,i)=Qvect(:,j)
     else  if( i.lt.0 .and. j.ge.0) then
        call dcopy(EP_dim,Qvect(1,j),1,dumQvect(1,-i),1)
        !dumQvect(:,-i)=Qvect(:,j)
     else  if( i.ge.0 .and. j.lt.0) then
        call dcopy(EP_dim,dumQvect(1,-j),1,Qvect(1,i),1)
        !Qvect(:,i)=dumQvect(:,-j)
     else 
        call dcopy(EP_dim,dumQvect(1,-j),1,dumQvect(1,-i),1)
        !dumQvect(:,-i)=dumQvect(:,-j)
     endif

  END SUBROUTINE EP_copy


  subroutine   EP_scalare_multik(i,j, scalari)
     use module_base
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
     call razero(ha%orbs%nkpts,ovrlp_local)
     call razero(ha%orbs%nkpts,ovrlp_global)
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
     use module_interfaces
     !Arguments
     real(wp), dimension(*), target :: psi ! per testare happlication
     !Local variables
     integer :: i

     if( ha%nproc/=1) then
        call transpose_v(ha%iproc,ha%nproc,ha%orbs,ha%Lzd%Glr%wfd,ha%comms,&
           &   psi(1),work=wrk,outadd=Qvect(1,0))  
     else
        do i=1, EP_dim_tot
           Qvect(i,0)= psi(i)
        enddo
     endif

  END SUBROUTINE EP_copia_per_prova

  subroutine EP_initialize_start()
     use module_interfaces
     !Local variables
     integer :: i

     if ( ha%at%paw_NofL(ha%at%iatype(ha%in_iat_absorber )).gt.0  ) then
        if (ha%iproc == 0) write(*,*) "USING PTILDES TO BUILD INITIAL WAVE"
        !!if(EP_dim_tot.gt.0) then
        call razero(EP_dim_tot  ,  Qvect_tmp  )
        call applyPAWprojectors(ha%orbs,ha%at,&
           &   ha%hx,ha%hy,ha%hz,ha%Lzd%Glr,ha%PAWD,Qvect_tmp,Qvect_tmp, ha%at%paw_S_matrices, &
           &   .true. ,    ha%in_iat_absorber, ha%Labsorber+1, &
           &   ha%Gabs_coeffs               ) 
        !!end if
     else
        STOP " ha%at%paw_NofL(ha%at%iatype(ha%in_iat_absorber )).gt.0  is false" 
        !!$       Note G%psiat and G%xp have now 2 dimenstions.
        !!$       call gaussians_to_wavelets_nonorm(ha%iproc,ha%nproc,ha%Lzd%Glr%geocode,ha%orbs,ha%Lzd%Glr%d,&
        !!$            ha%hx,ha%hy,ha%hz,ha%Lzd%Glr%wfd,EP_Gabsorber,ha%Gabs_coeffs,Qvect_tmp )
     endif

     !!$
     !!$    inquire(file='DOWRITEDOWNINITIALORBITAL', exist=exist)
     !!$    if(exist) then
     !!$       if(  sum( ha%at%paw_NofL ).gt.0 ) then
     !!$          if(  associated( ha%PAWD) ) then
     !!$             call razero(EP_dim_tot  ,  wrk  )
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
        call transpose_v(ha%iproc,ha%nproc,ha%orbs,ha%Lzd%Glr%wfd,ha%comms,&
           &   Qvect_tmp,work=wrk,outadd=Qvect(1,0))  
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
     use module_base
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

     !$omp parallel default (private) shared(z1,z3,kern_k1,kern_k2,kern_k3)&
     !$omp & shared(n1b,n3f,inzee,n1,n2,n3,ene,gamma)

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
     use module_base
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
     use module_base
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
     use module_interfaces
     use module_base
     !Arguments
     implicit none
     integer, intent(in) :: p,i
     real(gp) ene, gamma
     !Local variables
     integer :: k
     real(wp), parameter :: b2=24.8758460293923314d0,a2=3.55369228991319019d0
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
              &   Qvect(1:,i), work=wrk,outadd= Qvect_tmp(1) )  
        else
           call untranspose_v(ha%iproc,ha%nproc,ha%orbs,ha%Lzd%Glr%wfd,ha%comms,&
              &   dumQvect(1:,-i), work=wrk,outadd= Qvect_tmp(1) )  
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

     call dcopy(EP_dim_tot, Qvect_tmp(1),1,wrk(1),1) 


     if( dopcproj) then

        call razero(EP_dim_tot  ,  wrk1  )
        ha%PPD%iproj_to_factor(:) =  1.0_gp
        call applyPCprojectors(ha%orbs,ha%at,ha%hx,ha%hy,ha%hz,&
           &   ha%Lzd%Glr,ha%PPD,wrk,wrk1 )


        call razero(EP_dim_tot  ,  wrk2  )
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

        call razero(EP_dim_tot  ,  wrk1  )
        ha%PPD%iproj_to_factor(:) =  1.0_gp
        call applyPCprojectors(ha%orbs,ha%at,ha%hx,ha%hy,ha%hz,&
           &   ha%Lzd%Glr,ha%PPD,wrk,wrk1)


        call axpy(EP_dim_tot, -1.0_gp  ,  wrk1(1)   , 1,  wrk(1) , 1)
        call axpy(EP_dim_tot, +1.0_gp  ,  wrk2(1)   , 1,  wrk(1) , 1)
     endif


     if(  ha%iproc ==0 ) write(*,*)" done "



     if(p<0) then
        if(  ha%nproc/=1) then
           call transpose_v(ha%iproc,ha%nproc,ha%orbs,ha%Lzd%Glr%wfd,ha%comms,&
              &   wrk , work= Qvect_tmp ,outadd=dumQvect(1,-p))  
        else
           do k=1, EP_dim_tot
              dumQvect(k,-p) =  wrk(k)
           enddo
        endif
     else
        if(  ha%nproc/=1) then
           call transpose_v(ha%iproc,ha%nproc,ha%orbs,ha%Lzd%Glr%wfd,ha%comms,&
              &   wrk , work= Qvect_tmp ,outadd=Qvect(1,p))  
        else
           do k=1, EP_dim_tot
              Qvect(k,p) =  wrk(k)
           enddo
        endif
     endif

     call deallocate_work_arrays('P',.true.,1,w)

  END SUBROUTINE EP_precondition


  subroutine EP_Moltiplica4spectra(p,i, ene, gamma)
     use module_interfaces
     use gaussians, only: gaussian_basis
     !Arguments
     implicit none
     integer, intent(in) :: p,i
     real(gp) :: ene, gamma
     !Local variables
     integer :: k
     type(confpot_data), dimension(ha%orbs%norbp) :: confdatarr

     if( ha%nproc > 1) then
        if(i>=0) then
           call untranspose_v(ha%iproc,ha%nproc,ha%orbs,ha%Lzd%Glr%wfd,ha%comms,&
              &   Qvect(1:,i), work=wrk,outadd= Qvect_tmp(1) )  
        else
           call untranspose_v(ha%iproc,ha%nproc,ha%orbs,ha%Lzd%Glr%wfd,ha%comms,&
              &   dumQvect(1:,-i), work=wrk,outadd= Qvect_tmp(1) )  
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

     call FullHamiltonianApplication(ha%iproc,ha%nproc,ha%at,ha%orbs,ha%rxyz,&
          ha%proj,ha%Lzd,ha%nlpspd,confdatarr,ha%ngatherarr,ha%potential,Qvect_tmp,wrk,&
          ha%energs,ha%SIC,ha%GPU)

     call axpy(EP_dim_tot, -ene  ,  Qvect_tmp(1)   , 1,  wrk(1) , 1)
     call vcopy(EP_dim_tot,wrk(1),1,Qvect_tmp(1),1)
     !Qvect_tmp   =  wrk

     call FullHamiltonianApplication(ha%iproc,ha%nproc,ha%at,ha%orbs,ha%rxyz,&
          ha%proj,ha%Lzd,ha%nlpspd,confdatarr,ha%ngatherarr,ha%potential,Qvect_tmp,wrk,&
          ha%energs,ha%SIC,ha%GPU)

     call axpy(EP_dim_tot, -ene  ,  Qvect_tmp(1)   , 1,  wrk(1) , 1)

     if(  ha%iproc ==0 ) write(*,*)" done "

     if(p<0) then
        if(  ha%nproc/=1) then
           call transpose_v(ha%iproc,ha%nproc,ha%orbs,ha%Lzd%Glr%wfd,ha%comms,&
              &   wrk , work= Qvect_tmp ,outadd=dumQvect(1,-p))  
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
           call transpose_v(ha%iproc,ha%nproc,ha%orbs,ha%Lzd%Glr%wfd,ha%comms,&
              &   wrk , work= Qvect_tmp ,outadd=Qvect(1,p))  
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
     use module_interfaces
     !Arguments
     implicit none
     integer, intent(in) :: p,i
     !Local variables
     integer :: k
     type(confpot_data), dimension(ha%orbs%norbp) :: confdatarr

     if( ha%nproc > 1) then
        if(i>=0) then
           call untranspose_v(ha%iproc,ha%nproc,ha%orbs,ha%Lzd%Glr%wfd,ha%comms,&
              &   Qvect(1:,i), work=wrk,outadd= Qvect_tmp(1) )  
        else
           call untranspose_v(ha%iproc,ha%nproc,ha%orbs,ha%Lzd%Glr%wfd,ha%comms,&
              &   dumQvect(1:,-i), work=wrk,outadd= Qvect_tmp(1) )  
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

     call FullHamiltonianApplication(ha%iproc,ha%nproc,ha%at,ha%orbs,ha%rxyz,&
          ha%proj,ha%Lzd,ha%nlpspd,confdatarr,ha%ngatherarr,ha%potential,Qvect_tmp,wrk,&
          ha%energs,ha%SIC,ha%GPU)

     if(  ha%iproc ==0 ) write(*,*)" done "


     if(  sum( ha%at%paw_NofL ).gt.0 ) then
        if(associated( ha%PAWD) ) then
           call applyPAWprojectors(ha%orbs,ha%at,&
              &   ha%hx,ha%hy,ha%hz,ha%Lzd%Glr,ha%PAWD,Qvect_tmp,wrk,ha%at%paw_H_matrices, .false.  )

           !!$          call razero(EP_dim_tot  ,  wrk1  )
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
           call transpose_v(ha%iproc,ha%nproc,ha%orbs,ha%Lzd%Glr%wfd,ha%comms,&
              &   wrk , work= Qvect_tmp ,outadd=dumQvect(1,-p))  
        else
           do k=1, EP_dim_tot
              dumQvect(k,-p) =  wrk(k)
           enddo
        endif
     else
        if(  ha%nproc/=1) then
           call transpose_v(ha%iproc,ha%nproc,ha%orbs,ha%Lzd%Glr%wfd,ha%comms,&
              &   wrk , work= Qvect_tmp ,outadd=Qvect(1,p))  
        else
           do k=1, EP_dim_tot
              Qvect(k,p) =  wrk(k)
           enddo
        endif
     endif
     return 
  END SUBROUTINE EP_Moltiplica

  subroutine EP_ApplySinv(p,i)
     use module_interfaces
     !Arguments
     implicit none
     integer, intent(in) :: p,i
     !Local variables
     integer :: k

     if( ha%nproc > 1) then
        if(i>=0) then
           call untranspose_v(ha%iproc,ha%nproc,ha%orbs,ha%Lzd%Glr%wfd,ha%comms,&
              &   Qvect(1:,i), work=wrk,outadd= Qvect_tmp(1) )  
        else
           call untranspose_v(ha%iproc,ha%nproc,ha%orbs,ha%Lzd%Glr%wfd,ha%comms,&
              &   dumQvect(1:,-i), work=wrk,outadd= Qvect_tmp(1) )  
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
     call razero(EP_dim_tot  ,  wrk  )
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
           call transpose_v(ha%iproc,ha%nproc,ha%orbs,ha%Lzd%Glr%wfd,ha%comms,&
              &   wrk , work= Qvect_tmp ,outadd=dumQvect(1,-p))  
        else
           do k=1, EP_dim_tot
              dumQvect(k,-p) =  wrk(k)
           enddo
        endif
     else
        if(  ha%nproc/=1) then
           call transpose_v(ha%iproc,ha%nproc,ha%orbs,ha%Lzd%Glr%wfd,ha%comms,&
              &   wrk , work= Qvect_tmp ,outadd=Qvect(1,p))  
        else
           do k=1, EP_dim_tot
              Qvect(k,p) =  wrk(k)
           enddo
        endif
     endif
     return 
  END SUBROUTINE EP_ApplySinv




  subroutine EP_ApplyS(p,i)
     use module_interfaces
     !Arguments
     implicit none
     integer, intent(in) :: p,i
     !Local variables
     integer :: k

     if( ha%nproc > 1) then
        if(i>=0) then
           call untranspose_v(ha%iproc,ha%nproc,ha%orbs,ha%Lzd%Glr%wfd,ha%comms,&
              &   Qvect(1:,i), work=wrk,outadd= Qvect_tmp(1) )  
        else
           call untranspose_v(ha%iproc,ha%nproc,ha%orbs,ha%Lzd%Glr%wfd,ha%comms,&
              &   dumQvect(1:,-i), work=wrk,outadd= Qvect_tmp(1) )  
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
     call razero(EP_dim_tot  ,  wrk  )
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
           call transpose_v(ha%iproc,ha%nproc,ha%orbs,ha%Lzd%Glr%wfd,ha%comms,&
              &   wrk , work= Qvect_tmp ,outadd=dumQvect(1,-p))  
        else
           do k=1, EP_dim_tot
              dumQvect(k,-p) =  wrk(k)
           enddo
        endif
     else
        if(  ha%nproc/=1) then
           call transpose_v(ha%iproc,ha%nproc,ha%orbs,ha%Lzd%Glr%wfd,ha%comms,&
              &   wrk , work= Qvect_tmp ,outadd=Qvect(1,p))  
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
     use module_base
     use module_types
     use gaussians
     implicit none
     character(len=1), intent(in) :: geocode
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
     integer :: i_stat,i_all,ishell,iexpo,icoeff,iat,isat,ng,l,m,iorb,jorb,nterm,ierr,ispinor
     real(dp) :: normdev,tt,scpr,totnorm
     real(gp) :: rx,ry,rz
     integer, dimension(nterm_max) :: lx,ly,lz
     real(gp), dimension(nterm_max) :: fac_arr
     real(wp), dimension(:), allocatable :: tpsi

     if(iproc == 0 .and. verbose > 1) write(*,'(1x,a)',advance='no')'Writing wavefunctions in wavelet form '

     allocate(tpsi(wfd%nvctr_c+7*wfd%nvctr_f+ndebug),stat=i_stat)
     call memocc(i_stat,tpsi,'tpsi',subname)

     !initialize the wavefunction
     call razero((wfd%nvctr_c+7*wfd%nvctr_f)*orbs%norbp*orbs%nspinor,psi)
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

     i_all=-product(shape(tpsi))*kind(tpsi)
     deallocate(tpsi,stat=i_stat)
     call memocc(i_stat,i_all,'tpsi',subname)

  END SUBROUTINE gaussians_to_wavelets_nonorm


  subroutine lowpass_projector(n1,n2,n3,nvctr,psi)
     use module_base
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
        allocate(psi_gross(0:n1/2,2,0:n2/2,2,0:n3/2,2+ndebug),stat=i_stat)
        call memocc(i_stat,psi_gross,'psi_gross',subname)
        allocate(logrid(0:n1/2,0:n2/2,0:n3/2+ndebug),stat=i_stat)
        call memocc(i_stat,logrid,'logrid',subname)
     endif

     call analyse_per_self(n1/2,n2/2,n3/2,psi,psi_gross)

     !!$
     !!$
     !!$  ! coarse grid quantities
     !!$  call fill_logrid(ha%at%geocode,n1/2,n2/2,n3/2,0,n1/2,0,n2/2,0,n3/2,0,ha%at%nat,&
     !!$       ha%at%ntypes,ha%at%iatype,ha%rxyz,ha%radii_cf(1,1),8.0_gp,2.0_gp*ha%hx,2.0_gp*ha%hy,2.0_gp*ha%hz,logrid)

     !!$
     do ia =1, ha%at%nat
        iat=ha%at%iatype(ia)
        radius_gross = ha%radii_cf(  iat ,1 )*8
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
              do ia =1, ha%at%nat
                 iat=ha%at%iatype(ia)
                 radius_gross = ha%radii_cf(  iat ,1 )*8
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

END MODULE lanczos_interface


subroutine applyPAWprojectors(orbs,at,&
      &   hx,hy,hz,Glr,PAWD,psi,hpsi,  paw_matrix, dosuperposition , &
      &   sup_iatom, sup_l, sup_arraym)
   use module_base
   use module_types
   use module_interfaces, except_this_one => applyPAWprojectors

   implicit none

   type(orbitals_data), intent(inout) :: orbs
   type(atoms_data) :: at
   real(gp), intent(in) :: hx,hy,hz
   type(locreg_descriptors), intent(in) :: Glr
   type(pawproj_data_type) ::PAWD
   real(wp), dimension(:), pointer :: psi, hpsi, paw_matrix
   logical ::  dosuperposition
   integer , optional :: sup_iatom, sup_l
   real(wp) , dimension(:), pointer, optional :: sup_arraym 
   ! local variables
   character(len=*), parameter :: subname='applyPAWprojectors'
   integer :: ikpt, istart_ck, ispsi_k, isorb,ieorb, ispsi, iproj, istart_c,&
      &   mproj, mdone, ispinor, istart_c_i, mbvctr_c, mbvctr_f, mbseg_c, mbseg_f, &
      &   jseg_c, iproj_old, iorb, ncplx, l, jorb, lsign, ncplx_global
   real(gp) :: eproj_spinor
   real(gp) :: kx,ky,kz

   integer , parameter :: dotbuffersize = 1000
   real(dp)  :: dotbuffer(dotbuffersize), dotbufferbis(dotbuffersize)
   integer :: ibuffer, ichannel, nchannels, imatrix
   logical :: lfound_sup
   integer :: iat, old_istart_c, iatat , m, nspinor

   if (orbs%norbp.gt.0) then


      !apply the projectors  k-point of the processor
      !starting k-point
      ikpt=orbs%iokpt(1)
      istart_ck=1
      ispsi_k=1
      imatrix=1

      !!$ check that the coarse wavelets cover the whole box
      if(Glr%wfd%nvctr_c .ne. &
           ((Glr%d%n1+1)*(Glr%d%n2+1)*(Glr%d%n3+1))) then
         print *," WARNING : coarse wavelets dont cover the whole box "
      endif

      if(dosuperposition) then 
         lfound_sup=.false.
      endif
      ncplx_global=min(orbs%nspinor,2)

      loop_kpt: do

         call orbs_in_kpt(ikpt,orbs,isorb,ieorb,nspinor)

         ! loop over all my orbitals
         do iorb=isorb,ieorb
            istart_c=istart_ck
            iproj=1
            iat=0

            do iatat=1, at%nat
               if (  at%paw_NofL(at%iatype(iatat)).gt.0  ) then
                  iat=iat+1
                  istart_c_i=istart_c
                  iproj_old=iproj
                  ispsi=ispsi_k
                  !!!! do iorb=isorb,ieorb


                  mproj= PAWD%ilr_to_mproj(iat)

                  if( ikpt .ne. orbs%iokpt(iorb) ) then
                     STOP " ikpt .ne. orbs%iokpt(iorb) in applypawprojectors " 
                  end if
                  kx=orbs%kpts(1,ikpt)
                  ky=orbs%kpts(2,ikpt)
                  kz=orbs%kpts(3,ikpt)
                  call ncplx_kpt(orbs%iokpt(iorb),orbs,ncplx)

                  do ispinor=1,orbs%nspinor,ncplx_global
                     eproj_spinor=0.0_gp
                     if (ispinor >= 2) istart_c=istart_c_i
                     call plr_segs_and_vctrs(PAWD%paw_nlpspd%plr(iat),&
                          mbseg_c,mbseg_f,mbvctr_c,mbvctr_f)
                     jseg_c=1
!!$                     mbvctr_c=PAWD%paw_nlpspd%nvctr_p(2*iat-1)-PAWD%paw_nlpspd%nvctr_p(2*iat-2)
!!$                     mbvctr_f=PAWD%paw_nlpspd%nvctr_p(2*iat  )-PAWD%paw_nlpspd%nvctr_p(2*iat-1)
!!$                     mbseg_c=PAWD%paw_nlpspd%nseg_p(2*iat-1)-PAWD%paw_nlpspd%nseg_p(2*iat-2)
!!$                     mbseg_f=PAWD%paw_nlpspd%nseg_p(2*iat  )-PAWD%paw_nlpspd%nseg_p(2*iat-1)
!!$                     jseg_c=PAWD%paw_nlpspd%nseg_p(2*iat-2)+1
                     mdone=0
                     iproj=iproj_old
                     if(mproj>0) then
                        if(  PAWD%DistProjApply) then
                           jorb=1
                           do while(jorb<=PAWD%G%ncoeff .and. PAWD%iorbtolr(jorb)/= iat) 
                              jorb=jorb+1
                           end do
                           if(jorb<PAWD%G%ncoeff) then
                              call fillPawProjOnTheFly(PAWD,Glr,iat,hx,hy,hz,&
                                 &   kx,ky,kz,&
                                 &   jorb,istart_c,at%geocode,at,iatat ) 
                           endif
                        end if
                     endif

                     do while(mdone< mproj)
                        lsign = PAWD%iproj_to_l(iproj)
                        l=abs(lsign)
                        ibuffer = 0
                        nchannels =  PAWD% iproj_to_paw_nchannels(iproj)
                        imatrix=PAWD%iprojto_imatrixbeg(iproj)
                        !!$
                        !!$                    if(.not. dosuperposition) then
                        !!$                       print *, "applying paw for l= ", l,&
                        !!$                            "  primo elemento ", paw_matrix(PAWD%iprojto_imatrixbeg(iproj))
                        !!$                    end if
                        old_istart_c=istart_c
                        do ichannel=1, nchannels
                           do m=1,2*l-1
                              ibuffer=ibuffer+1
                              if(ibuffer.gt.dotbuffersize ) then
                                 STOP 'ibuffer.gt.dotbuffersize'
                              end if

                              if( .not. dosuperposition .and. lsign>0 ) then
                                 call wpdot_wrap(ncplx,  &
                                      Glr%wfd%nvctr_c,Glr%wfd%nvctr_f,&
                                      Glr%wfd%nseg_c,Glr%wfd%nseg_f,&
                                      Glr%wfd%keyvglob(1),Glr%wfd%keyglob(1,1),&
                                      psi(ispsi+&
                                      (ispinor-1)*(orbs%npsidim_orbs/orbs%nspinor)),&
                                      mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
!!$                                      PAWD%paw_nlpspd%keyv_p(jseg_c),&
!!$                                      PAWD%paw_nlpspd%keyg_p(1,jseg_c),&
                                      PAWD%paw_nlpspd%plr(iat)%wfd%keyvglob(jseg_c),&
                                      PAWD%paw_nlpspd%plr(iat)%wfd%keyglob(1,jseg_c),&
                                      PAWD%paw_proj(istart_c),&
                                      dotbuffer( ibuffer ) )
                              end if
                              ibuffer=ibuffer + (ncplx-1)

                              !!$                          !! TTTTTTTTTTTTTTTTTTTt TEST TTTTTTTTTTTTTTTTTTT
                              !!$                          call wpdot_wrap(ncplx,  &
                              !!$                               mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,PAWD%paw_nlpspd%keyv_p(jseg_c),&
                              !!$                               PAWD%paw_nlpspd%keyg_p(1,jseg_c),PAWD%paw_proj(istart_c),& 
                              !!$                               mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,PAWD%paw_nlpspd%keyv_p(jseg_c),&
                              !!$                               PAWD%paw_nlpspd%keyg_p(1,jseg_c),&
                              !!$                               PAWD%paw_proj(istart_c),&
                              !!$                               eproj_spinor)
                              !!$                          print *, "TEST:  THE PROJECTOR ichannel = ", ichannel, " m=",m, " HAS SQUARED MODULUS  " ,eproj_spinor 

                              !!$                          !! plot -------------------------------------------------------------
                              !!$                          Plr%d%n1 = Glr%d%n1
                              !!$                          Plr%d%n2 = Glr%d%n2
                              !!$                          Plr%d%n3 = Glr%d%n3
                              !!$                          Plr%geocode = at%geocode
                              !!$                          Plr%wfd%nvctr_c  =PAWD%paw_nlpspd%nvctr_p(2*iat-1)-PAWD%paw_nlpspd%nvctr_p(2*iat-2)
                              !!$                          Plr%wfd%nvctr_f  =PAWD%paw_nlpspd%nvctr_p(2*iat  )-PAWD%paw_nlpspd%nvctr_p(2*iat-1)
                              !!$                          Plr%wfd%nseg_c   =PAWD%paw_nlpspd%nseg_p(2*iat-1)-PAWD%paw_nlpspd%nseg_p(2*iat-2)
                              !!$                          Plr%wfd%nseg_f   =PAWD%paw_nlpspd%nseg_p(2*iat  )-PAWD%paw_nlpspd%nseg_p(2*iat-1)
                              !!$                          call allocate_wfd(Plr%wfd,subname)
                              !!$                          Plr%wfd%keyv(:)  = PAWD%paw_nlpspd%keyv_p( PAWD%paw_nlpspd%nseg_p(2*iat-2)+1:&
                              !!$                               PAWD%paw_nlpspd%nseg_p(2*iat)   )
                              !!$                          Plr%wfd%keyg(1:2, :)  = PAWD%paw_nlpspd%keyg_p( 1:2,  PAWD%paw_nlpspd%nseg_p(2*iat-2)+1:&
                              !!$                               PAWD%paw_nlpspd%nseg_p(2*iat)   )
                              !!$
                              !!$                          !! ---------------  use this to plot projectors
                              !!$                          write(orbname,'(A,i4.4)')'paw_',iproj
                              !!$                          Plr%bounds = Glr%bounds
                              !!$                          Plr%d          = Glr%d
                              !!$                          call plot_wf_cube(orbname,at,Plr,hx,hy,hz,rxyz, PAWD%paw_proj(istart_c) ,"1234567890" ) 
                              !!$                          !! END plot ----------------------------------------------------------


                              !! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

                              istart_c=istart_c+(mbvctr_c+7*mbvctr_f)*ncplx
                              iproj=iproj+1
                              mdone=mdone+1
                           end do
                        end do

                        !!$                    call DGEMM('N','N', nchannels ,(2*l-1)  , nchannels  ,&
                        !!$                         1.0d0 , paw_matrix(imatrix) , nchannels ,&
                        !!$                         dotbuffer  ,(2*l-1), 0.0D0 , dotbufferbis  ,(2*l-1))


                        if( .not. dosuperposition) then
                           if(lsign>0) then
                              call DGEMM('N','N',(2*l-1)*ncplx  , nchannels , nchannels  ,&
                                 &   1.0d0 ,dotbuffer , (2*l-1)*ncplx ,&
                                 &   paw_matrix(imatrix)  ,nchannels , 0.0D0 , dotbufferbis  ,(2*l-1)*ncplx )
                           else
                              dotbufferbis=0.0_wp
                           endif
                        else
                           !print *,'here',nchannels
                           if( sup_iatom .eq. iatat .and. (-sup_l) .eq. lsign ) then
                              do ichannel=1, nchannels
                                 do m=1,2*l-1
                                    dotbufferbis((ichannel-1)*(2*l-1)*ncplx+m*ncplx           ) = 0.0_gp ! keep this before
                                    dotbufferbis((ichannel-1)*(2*l-1)*ncplx+m*ncplx -(ncplx-1)) = sup_arraym(m)
                                 end do
                              enddo
                              lfound_sup=.true.
                           else
                              do ichannel=1, nchannels
                                 do m=1,2*l-1
                                    dotbufferbis((ichannel-1)*(2*l-1)*ncplx+m*ncplx           ) = 0.0_gp 
                                    dotbufferbis((ichannel-1)*(2*l-1)*ncplx+m*ncplx -(ncplx-1)) = 0.0_gp
                                 end do
                              enddo
                           endif
                        endif


                        ibuffer=0
                        iproj    =  iproj  - nchannels * ( 2*l-1 )
                        istart_c = old_istart_c

                        do ichannel=1, nchannels
                           do m=1,2*l-1
                              ibuffer=ibuffer+1

                              call waxpy_wrap(ncplx,dotbufferbis( ibuffer ) ,&
                                   mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
!!$                                 &   PAWD%paw_nlpspd%keyv_p(jseg_c),PAWD%paw_nlpspd%keyg_p(1,jseg_c),&
                                   PAWD%paw_nlpspd%plr(iat)%wfd%keyvglob(jseg_c),&
                                   PAWD%paw_nlpspd%plr(iat)%wfd%keyglob(1,jseg_c),&
                                   PAWD%paw_proj(istart_c),&
                                   Glr%wfd%nvctr_c,Glr%wfd%nvctr_f,Glr%wfd%nseg_c,Glr%wfd%nseg_f,&
                                   Glr%wfd%keyvglob(1),Glr%wfd%keyglob(1,1),&
                                   hpsi(ispsi+(ispinor-1)*(orbs%npsidim_orbs/orbs%nspinor)  )&
                                   )


                              istart_c=istart_c+(mbvctr_c+7*mbvctr_f)*ncplx
                              iproj=iproj+1
                              ibuffer=ibuffer + (ncplx-1)
                           end do
                        end do
                     end do

                     mdone=0
                     !!$ iproj=iproj_old
                     istart_c=istart_c_i
                     ispsi=ispsi+(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*nspinor

                  end do

                  if( PAWD%DistProjApply ) then
                     istart_c=1
                  else
                     istart_c=istart_c+(mbvctr_c+7*mbvctr_f)*mproj*ncplx
                  endif

               end if
            end do

            ispsi_k=ispsi
         end do
         istart_ck=istart_c

         if(  dosuperposition ) then
            if(.not. lfound_sup) then
               print *, " initial state not found in routine ",subname
               STOP 
            endif
         endif


         if (ieorb == orbs%norbp) exit loop_kpt
         ikpt=ikpt+1


      end do loop_kpt
   end if
END SUBROUTINE applyPAWprojectors


subroutine applyPCprojectors(orbs,at,&
      &   hx,hy,hz,Glr,PPD,psi,hpsi, dotest)

   use module_base
   use module_types
   use module_interfaces, except_this_one => applyPCprojectors

   implicit none

   type(orbitals_data), intent(inout) :: orbs
   type(atoms_data) :: at
   real(gp), intent(in) :: hx,hy,hz
   type(locreg_descriptors), intent(in) :: Glr
   type(pcproj_data_type) ::PPD
   real(wp), dimension(:), pointer :: psi, hpsi
   logical, optional :: dotest


   ! local variables
   character(len=*), parameter :: subname='applyPCprojectors'
   character(len=11) :: orbname
   type(locreg_descriptors) :: Plr
   integer :: ikpt, istart_ck, ispsi_k, isorb,ieorb, ispsi, iproj, istart_c,&
      &   mproj, mdone, ispinor, istart_c_i, mbvctr_c, mbvctr_f, mbseg_c, mbseg_f, &
      &   jseg_c, iproj_old, iorb, ncplx, l, i, jorb
   real(gp) eproj_spinor,psppar_aux(0:4, 0:6)
   integer :: iat , nspinor

   !apply the projectors  k-point of the processor
   !starting k-point
   ikpt=orbs%iokpt(1)
   istart_ck=1
   ispsi_k=1
   loop_kpt: do

      call orbs_in_kpt(ikpt,orbs,isorb,ieorb,nspinor)

      ! loop over all my orbitals
      istart_c=1
      iproj=1

      do iat=1,at%nat
         istart_c_i=istart_c
         iproj_old=iproj
         ispsi=ispsi_k
         do iorb=isorb,ieorb

            mproj= PPD%ilr_to_mproj(iat)

            call ncplx_kpt(orbs%iokpt(iorb),orbs,ncplx)


            do ispinor=1,orbs%nspinor,ncplx
               eproj_spinor=0.0_gp

               if (ispinor >= 2) istart_c=istart_c_i

               call plr_segs_and_vctrs(PPD%pc_nlpspd%plr(iat),&
                    mbseg_c,mbseg_f,mbvctr_c,mbvctr_f)
               jseg_c=1

!!$               mbvctr_c=PPD%pc_nlpspd%nvctr_p(2*iat-1)-PPD%pc_nlpspd%nvctr_p(2*iat-2)
!!$               mbvctr_f=PPD%pc_nlpspd%nvctr_p(2*iat  )-PPD%pc_nlpspd%nvctr_p(2*iat-1)
!!$
!!$               mbseg_c=PPD%pc_nlpspd%nseg_p(2*iat-1)-PPD%pc_nlpspd%nseg_p(2*iat-2)
!!$               mbseg_f=PPD%pc_nlpspd%nseg_p(2*iat  )-PPD%pc_nlpspd%nseg_p(2*iat-1)
!!$               jseg_c=PPD%pc_nlpspd%nseg_p(2*iat-2)+1

               mdone=0
               iproj=iproj_old

               if(mproj>0) then
                  if(  PPD%DistProjApply) then
                     jorb=1
                     do while( jorb<=PPD%G%ncoeff         .and. PPD%iorbtolr(jorb)/= iat) 
                        jorb=jorb+1
                     end do
                     if(jorb<PPD%G%ncoeff) then
                        call fillPcProjOnTheFly(PPD,Glr,iat,at,hx,hy,hz,&
                             jorb,PPD%ecut_pc,istart_c) 
                     endif
                  end if
               endif


               do while(mdone< mproj)

                  l = PPD%iproj_to_l(iproj)

                  i=1
                  psppar_aux=0.0_gp
                  !! psppar_aux(l,i)=1.0_gp/PPD%iproj_to_ene(iproj)
                  !! psppar_aux(l,i)=1.0_gp  ! *iorb
                  psppar_aux(l,i)=PPD%iproj_to_factor(iproj)  


                  call applyprojector(ncplx,l,i, psppar_aux(0,0), 2 ,&
                       Glr%wfd%nvctr_c,Glr%wfd%nvctr_f, Glr%wfd%nseg_c,&
                       Glr%wfd%nseg_f,&
                       Glr%wfd%keyvglob(1),Glr%wfd%keyglob(1,1),&
                       mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
                       PPD%pc_nlpspd%plr(iat)%wfd%keyvglob(jseg_c),&
                       PPD%pc_nlpspd%plr(iat)%wfd%keyglob(1,jseg_c),&
!!$                       PPD%pc_nlpspd%keyv_p(jseg_c),PPD%pc_nlpspd%keyg_p(1,jseg_c),&
                       PPD%pc_proj(istart_c),&
                       psi(ispsi+ (ispinor-1)*(orbs%npsidim_orbs/orbs%nspinor)  ),&
                       hpsi(ispsi+(ispinor-1)*(orbs%npsidim_orbs/orbs%nspinor)  ),&
                       eproj_spinor)


                  if(iorb==1) then         
                     if( present(dotest) ) then
                        eproj_spinor=0.0_gp
!!$                        call wpdot_wrap(ncplx,  &
!!$                             mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
!!$                             PPD%pc_nlpspd%keyv_p(jseg_c),&
!!$                             PPD%pc_nlpspd%keyg_p(1,jseg_c),&
!!$                             PPD%pc_proj(istart_c),& 
!!$                             mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,&
!!$                             PPD%pc_nlpspd%keyv_p(jseg_c),&
!!$                             PPD%pc_nlpspd%keyg_p(1,jseg_c),&
!!$                             PPD%pc_proj(istart_c),&
!!$                             eproj_spinor)
                        print *, " IL PROIETTORE HA MODULO QUADRO  " ,eproj_spinor 
                        if(dotest) then
                           !! ---------------  use this to plot projectors
                           write(orbname,'(A,i4.4)')'pc_',iproj                     
                           Plr%d%n1=Glr%d%n1
                           Plr%d%n2=Glr%d%n2
                           Plr%d%n3=Glr%d%n3
                           Plr%geocode=at%geocode
                           
                           call plr_segs_and_vctrs(PPD%pc_nlpspd%plr(iat),&
                                Plr%wfd%nseg_c,Plr%wfd%nseg_f,&
                                Plr%wfd%nvctr_c,Plr%wfd%nvctr_f)                  
!!$                           Plr%wfd%nvctr_c  =PPD%pc_nlpspd%nvctr_p(2*iat-1)-PPD%pc_nlpspd%nvctr_p(2*iat-2)
!!$                           Plr%wfd%nvctr_f  =PPD%pc_nlpspd%nvctr_p(2*iat  )-PPD%pc_nlpspd%nvctr_p(2*iat-1)
!!$                           Plr%wfd%nseg_c   =PPD%pc_nlpspd%nseg_p(2*iat-1 )-PPD%pc_nlpspd%nseg_p(2*iat-2)
!!$                           Plr%wfd%nseg_f   =PPD%pc_nlpspd%nseg_p(2*iat  ) -PPD%pc_nlpspd%nseg_p(2*iat-1)
                           call allocate_wfd(Plr%wfd,subname)
!!$                           Plr%wfd%keyv(:)=PPD%pc_nlpspd%keyv_p(PPD%pc_nlpspd%nseg_p(2*iat-2)+1:&
!!$                              &   PPD%pc_nlpspd%nseg_p(2*iat)   )
!!$                           Plr%wfd%keyg(1:2, :)  = PPD%pc_nlpspd%keyg_p( 1:2,  PPD%pc_nlpspd%nseg_p(2*iat-2)+1:&
!!$                              &   PPD%pc_nlpspd%nseg_p(2*iat)   )
                           Plr%bounds = Glr%bounds
                           Plr%d          = Glr%d                    
                           !! call plot_wf_cube(orbname,at,Plr,hx,hy,hz,rxyz, PPD%pc_proj(istart_c) ,"1234567890" ) 
                           call deallocate_wfd(Plr%wfd,subname)
                        endif
                     endif

                  endif
                  istart_c=istart_c+(mbvctr_c+7*mbvctr_f)*(2*l-1)*ncplx
                  iproj=iproj+(2*l-1)
                  mdone=mdone+(2*l-1)
               end do
            end  do
            istart_c=istart_c_i

            if( present(dotest) ) then
               if(dotest) then
                  eproj_spinor=0.0_gp
!!$                  call wpdot_wrap(ncplx,  &
!!$                     &   Glr%wfd%nvctr_c,Glr%wfd%nvctr_f, Glr%wfd%nseg_c, Glr%wfd%nseg_f,&
!!$                     &   Glr%wfd%keyv(1),Glr%wfd%keyg(1,1), hpsi(ispsi + 0 ), &
!!$                     &   Glr%wfd%nvctr_c,Glr%wfd%nvctr_f, Glr%wfd%nseg_c, Glr%wfd%nseg_f,&
!!$                     &   Glr%wfd%keyv(1),Glr%wfd%keyg(1,1), hpsi(ispsi + 0 ),  &
!!$                     &   eproj_spinor)
                  print *, "hpsi  HA MODULO QUADRO  " ,eproj_spinor 
                  eproj_spinor=0.0_gp
!!$                  call wpdot_wrap(ncplx,  &
!!$                     &   Glr%wfd%nvctr_c,Glr%wfd%nvctr_f, Glr%wfd%nseg_c, Glr%wfd%nseg_f,&
!!$                     &   Glr%wfd%keyv(1),Glr%wfd%keyg(1,1), psi(ispsi + 0 ), &
!!$                     &   Glr%wfd%nvctr_c,Glr%wfd%nvctr_f, Glr%wfd%nseg_c, Glr%wfd%nseg_f,&
!!$                     &   Glr%wfd%keyv(1),Glr%wfd%keyg(1,1), psi(ispsi + 0 ),  &
!!$                     &   eproj_spinor)
                  print *, "psi  HA MODULO QUADRO  " ,eproj_spinor         
                  !! CECCARE IPROJ = mproj tot, istart_c=nelproj 
                  write(orbname,'(A,i4.4)')'pcorb_',iorb
                  !! call plot_wf_cube(orbname,at,Glr,hx,hy,hz,rxyz,hpsi(ispsi + 0 ),"dopoprec.." ) ! solo spinore 1
               end if
            endif

            ispsi=ispsi+(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*nspinor

         end do


         if( PPD%DistProjApply ) then
            istart_c=1
         else
            istart_c=istart_c+(mbvctr_c+7*mbvctr_f)*mproj
         endif

      end do

      !! istart_ck=istart_c  non si incrementa

      if (ieorb == orbs%norbp) exit loop_kpt
      ikpt=ikpt+1
      ispsi_k=ispsi
   end do loop_kpt

END SUBROUTINE applyPCprojectors  
