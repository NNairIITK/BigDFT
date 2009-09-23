module lanczos_interface
  use module_base
  use module_types
  implicit none

  private

  !calculate the allocation dimensions
  public :: EP_inizializza,EP_initialize_start,EP_allocate_for_eigenprob, get_EP_dim, set_EP_shift,&
       EP_mat_mult,EP_make_dummy_vectors,  EP_normalizza,EP_copy,EP_scalare, EP_copia_per_prova, &
       EP_set_all_random, EP_GramSchmidt, EP_add_from_vect_with_fact, EP_Moltiplica, EP_memorizza_stato, &
       EP_free, EP_initialize_start_0,  EP_norma2_initialized_state , EP_store_occupied_orbitals, EP_occprojections

 
  character(len=*), parameter :: subname='lanczos_interface'
  integer :: i_all,i_stat,ierr
  type(lanczos_args), pointer :: ha

  real(8), pointer :: matrix(:,:)
  real(wp), pointer :: Qvect   (:,:)
  real(wp), pointer :: dumQvect(:,:)
  real(wp), pointer :: occQvect(:)

  real(8), pointer :: passed_matrix(:,:)

  real(8) EP_shift
  integer EP_dim
  integer EP_dim_tot

  type(gaussian_basis), pointer :: EP_Gabsorber

  real(wp), dimension(:), pointer :: Qvect_tmp,wrk
  real(8)  EP_norma2_initialized_state

  logical EP_doorthoocc
  integer EP_norb
  real(wp), dimension(:), pointer :: EP_occprojections
contains

!!$  subroutine settapointer(p)
!!$    real(8), intent(inout), TARGET :: p(:,:)
!!$    passed_matrix => p 
!!$    return
!!$  end subroutine settapointer

!!$  subroutine testapointer()
!!$    integer i,j
!!$    integer , pointer :: shape_a(:)
!!$    print *, " ecco il pointer " 
!!$    print *, " shape is " , shape(passed_matrix)
!!$    do i=1,  size(passed_matrix, 1)
!!$       print *, (passed_matrix(i,j), j=1, size(passed_matrix,2))
!!$    enddo
!!$    return
!!$  end subroutine testapointer

  

  integer function get_EP_dim()
    get_EP_dim=EP_dim_tot
    return
  end function get_EP_dim

  subroutine set_EP_shift(shift)
    real(8) shift
    EP_shift=shift
    return
  end subroutine set_EP_shift

  subroutine EP_inizializza(ha_actual)
    implicit none
    type(lanczos_args), target :: ha_actual
    
    ha=>ha_actual

    EP_dim=ha%comms%nvctr_par(ha%iproc,1)*ha%orbs%nspinor

    EP_dim_tot=(ha%lr%wfd%nvctr_c+7*ha%lr%wfd%nvctr_f)*ha%orbs%nspinor

    if( (ha%lr%wfd%nvctr_c+7*ha%lr%wfd%nvctr_f) /= &
         sum(ha%comms%nvctr_par(:,1))  ) then
       stop "array size inconsistency" 
    endif
    
    allocate(Qvect_tmp(ha%orbs%norbp*ha%orbs%nspinor*(ha%lr%wfd%nvctr_c+7*ha%lr%wfd%nvctr_f  +ndebug)) ) 
    allocate(wrk      (ha%orbs%norbp*ha%orbs%nspinor*(ha%lr%wfd%nvctr_c+7*ha%lr%wfd%nvctr_f  +ndebug)) )
 
    call memocc(i_stat,Qvect_tmp,'Qvect_tmp',subname)
    call memocc(i_stat,wrk,'wrk',subname)

    EP_shift=0.0

    EP_norb=0
    EP_doorthoocc=.false.

  end subroutine EP_inizializza
  

  subroutine EP_mat_mult( m,k ,  EV   ) ! usare i dumvectors a partire da -1
    implicit none
    
    integer , intent(in):: m,k
    real(8), intent (in) :: EV(0:m-1,0:k-1)
! ::::::::::::::::::::::::::::::::::::
    integer i,j,l

    call gemm('N','N', EP_dim , k, m  ,1.0_wp ,&
         Qvect(1,0)  , EP_dim ,&
         EV(0 ,0) ,m ,&
         0.0_wp ,  dumQvect(1,1) ,  EP_dim )


!!$
!!$    do i=0, k-1
!!$       do j=1, EP_dim
!!$          dumQvect(j,i+1)=0.0d0
!!$          do l=0, m-1 ! controllare il range
!!$             dumQvect(j,i+1)=dumQvect(j,i+1)+EV( l,i )*Qvect(  j,l ) 
!!$          enddo 
!!$       enddo
!!$    enddo


    return 
  end subroutine EP_mat_mult
  

  !allocate the wavefunctions in the transposed form, for lancsoz


  subroutine EP_free()
  

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



    if(EP_doorthoocc) then
       i_all=-product(shape(EP_occprojections))*kind(EP_occprojections)
       deallocate(EP_occprojections,stat=i_stat)
       call memocc(i_stat,i_all,'EP_occprojections',subname)
    endif

  end subroutine EP_free


  subroutine  EP_allocate_for_eigenprob(nsteps)
    implicit none
    integer, intent(in):: nsteps
! :::::::::::::::::::::::
    allocate(Qvect(EP_dim,0:nsteps+ndebug),&
         stat=i_stat)
    call memocc(i_stat,Qvect,'Qvect',subname)
  end subroutine EP_allocate_for_eigenprob


  subroutine EP_make_dummy_vectors(nd)
    implicit none
    integer, intent(in):: nd
! :::::::::::::::::::::::    
    allocate(dumQvect(EP_dim,0:nd+1+ndebug),&
         stat=i_stat)
    call memocc(i_stat,dumQvect,'dumQvect',subname)
  end subroutine EP_make_dummy_vectors

  subroutine EP_store_occupied_orbitals(iproc, nproc, norb, Occ_psit )
    implicit none
    integer, intent(IN) :: iproc, nproc, norb
    real(wp), dimension(:), pointer :: Occ_psit

    ! :::::::::::::::::::::::    
    EP_doorthoocc=.true.
    EP_norb=norb
    occQvect=>Occ_psit
    allocate(EP_occprojections(EP_norb) , stat=i_stat )
    call memocc(i_stat,EP_occprojections,'EP_occprojections',subname)

  end subroutine EP_store_occupied_orbitals


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
  end subroutine EP_normalizza
    


  subroutine EP_normalizza_interno(Q)
    implicit none
    real(wp) Q(EP_dim)
! :::::::::::::::::::::::
    integer i,j,l
    real(wp) sump, sumtot

    sump = dot(EP_dim,Q(1),1 ,Q(1),1)
    sumtot=0


    if(ha%nproc/=1) then
       call MPI_Allreduce(sump,sumtot,1,mpidtypw, MPI_SUM,MPI_COMM_WORLD ,ierr )
    else
       sumtot=sump
    endif
    sump=1.0_wp/sqrt(sumtot)

    call vscal(EP_dim,sump, Q(1), 1  ) 

  end subroutine EP_normalizza_interno


  subroutine EP_copy(i,j)
    implicit none
    integer, intent(in) :: i,j

    integer k
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

  end subroutine EP_copy
  

  real(8) function   EP_scalare(i,j)
    implicit none
    integer, intent(in) :: i,j
    if( i.ge.0 .and. j.ge.0) then
       EP_scalare = EP_scalare_interna(Qvect(1:,i), Qvect(1:,j) )
    else  if( i.lt.0 .and. j.ge.0) then
       EP_scalare = EP_scalare_interna(dumQvect(1:,-i), Qvect(1:,j) )
    else  if( i.ge.0 .and. j.lt.0) then
       EP_scalare = EP_scalare_interna(Qvect(1:,i), dumQvect(1:,-j) )
    else 
       EP_scalare = EP_scalare_interna(dumQvect(1:,-i), dumQvect(1:,-j) )
    endif
    
    return 
  end function EP_scalare
  
 
  real(8) function EP_scalare_interna(a,b)
    implicit none
    real(8), intent(in):: a(EP_dim), b(EP_dim)
    ! ::::::::::::::::::::::::::::::::::::::::::::::
    integer i,j

    real(wp) sump, sumtot

    sump = dot(EP_dim,a(1),1 ,b(1),1)
    sumtot=0

    if(ha%nproc/=1) then
       call MPI_Allreduce(sump,sumtot,1,mpidtypw, MPI_SUM,MPI_COMM_WORLD ,ierr )
    else
       sumtot=sump
    endif
    EP_scalare_interna=sumtot

  end function EP_scalare_interna
    
  subroutine  EP_copia_per_prova(psi)
    use module_interfaces
    real(wp), dimension(*), target :: psi ! per testare happlication

    real(8) sumdum
    integer i

    if( ha%nproc/=1) then
       call transpose_v(ha%iproc,ha%nproc,ha%orbs,ha%lr%wfd,ha%comms,&
            psi(1),work=wrk,outadd=Qvect(1,0))  
    else
       do i=1, EP_dim_tot
          Qvect(i,0)= psi(i)
       enddo
    endif
    

  end subroutine EP_copia_per_prova
  
  subroutine EP_memorizza_stato(Gabsorber)
    type(gaussian_basis), target :: Gabsorber
  
    EP_Gabsorber => Gabsorber
    
  end subroutine EP_memorizza_stato


  subroutine EP_initialize_start_0( Gabsorber)
    use module_interfaces
    type(gaussian_basis) :: Gabsorber
    integer :: i

    call gaussians_to_wavelets_nonorm(ha%iproc,ha%nproc,ha%lr%geocode,ha%orbs,ha%lr%d,&
         ha%hx,ha%hy,ha%hz,ha%lr%wfd,Gabsorber,ha%Gabs_coeffs,Qvect_tmp )
       
    
    if (.not. associated(Qvect)) then
       write(*,*)'ERROR: initialization vector not allocated!'
       stop
    end if
    
    if( ha%nproc/=1) then
       call transpose_v(ha%iproc,ha%nproc,ha%orbs,ha%lr%wfd,ha%comms,&
            Qvect_tmp,work=wrk,outadd=Qvect(1,0))  
    else
       do i=1, EP_dim_tot
          Qvect(i,0)= Qvect_tmp(i)
       enddo
    endif
    


   end subroutine EP_initialize_start_0

  subroutine EP_initialize_start()
    use module_interfaces
    integer :: i,k, volta
    real(wp), pointer ::  scals(:), scalstot(:)

!!$    if( associated(EP_Gabsorber) ) then
!!$    if( ha%iproc == 0  ) then
       ! print *, " inizializzo da Gabsorber " 
!!$        print *, EP_Gabsorber%rxyz(1,1)
       call gaussians_to_wavelets_nonorm(ha%iproc,ha%nproc,ha%lr%geocode,ha%orbs,ha%lr%d,&
            ha%hx,ha%hy,ha%hz,ha%lr%wfd,EP_Gabsorber,ha%Gabs_coeffs,Qvect_tmp )
 
!!$
!!$       if(ha%iproc.eq.0) then
!!$          call plot_wf('testa',ha%lr,ha%hx,ha%hy,ha%hz,ha%rxyz(1,1),ha%rxyz(2,1),ha%rxyz(3,1),&
!!$               Qvect_tmp,'comm')
!!$       endif
!!$    endif
    if (.not. associated(Qvect)) then
       write(*,*)'ERROR: initialization vector not allocated!'
       stop
    end if
    

    if( ha%nproc/=1) then
       call transpose_v(ha%iproc,ha%nproc,ha%orbs,ha%lr%wfd,ha%comms,&
            Qvect_tmp,work=wrk,outadd=Qvect(1,0))  
    else
       do i=1, EP_dim_tot
          Qvect(i,0)= Qvect_tmp(i)
       enddo
    endif
    
    
    
    

    if(EP_doorthoocc) then
       print *, " doorthoocc " 
       allocate(scals( EP_norb) ,stat=i_stat)
       call memocc(i_stat,scals,'scals',subname)
       allocate(scalstot( EP_norb) ,stat=i_stat)
       call memocc(i_stat,scalstot,'scalstot',subname)
       
       print *, " EP_dim", EP_dim , " EP_norb  ",EP_norb,  ha%iproc
       print *," associated(occQvect)    " ,  associated(occQvect)  ,  ha%iproc
       print *, "  shape(occQvect )  " , shape(occQvect )  ,  ha%iproc
       

       do volta=1,2 
          call gemm('T','N', EP_norb , 1, EP_dim ,1.0_wp ,  occQvect(1)  ,&
               EP_dim ,  Qvect(1,0) ,EP_dim , 0.0_wp , scals(1) ,  EP_norb )
          
          if(ha%nproc/=1) then
             call MPI_Allreduce(scals(1) ,scalstot(1) ,EP_norb  ,mpidtypw, MPI_SUM,MPI_COMM_WORLD ,ierr )
          else
             do i=1,EP_norb
                scalstot(i)=scals(i)
             enddo
          endif
          do i=1, EP_norb
             call axpy(EP_dim, -scalstot(i)  ,   occQvect(1+i*EP_dim)  , 1, Qvect(1,0)  , 1)
             
             if(volta==1) then
                EP_occprojections(i)= scalstot(i)
             endif
          enddo
       enddo
       
       print *, 2

       i_all=-product(shape(scals))*kind(scals)
       deallocate(scals,stat=i_stat)
       call memocc(i_stat,i_all,'scals',subname)
       
       
       i_all=-product(shape(scalstot))*kind(scalstot)
       deallocate(scalstot,stat=i_stat)
       call memocc(i_stat,i_all,'scalstot',subname)
    endif







!!$    call untranspose_v(ha%iproc,ha%nproc,ha%orbs%norbp,ha%orbs%nspinor,ha%lr%wfd,ha%comms,&
!!$               Qvect(1,0), work=wrk,outadd= Qvect_tmp(1) )  
!!$
!!$    if(ha%iproc.eq.0) then
!!$       call plot_wf('testadopo',ha%lr,ha%hx,ha%hy,ha%hz,ha%rxyz(1,1),ha%rxyz(2,1),ha%rxyz(3,1),&
!!$            Qvect_tmp,'comm')
!!$    endif
!!$ else
!!$    call EP_set_all_random(0)
!!$ endif

    EP_norma2_initialized_state= EP_scalare(0,0)

end subroutine EP_initialize_start
     
     subroutine EP_set_all_random(i)
    implicit none
    integer, intent(in) :: i
    ! ::::::::::::::::::::::::::::::
    if( i.ge.0 ) then
       call EP_set_random_interna(Qvect(1:,i))
    else 
       call EP_set_random_interna(dumQvect(1:,-i))
    endif

  end subroutine EP_set_all_random


  subroutine EP_set_random_interna(Q)
    implicit none
    real(wp) Q(EP_dim)

    ! :::::::::::::::::::::::
    real(4) rand
    integer i
    do i=1,EP_dim
       call random_number(rand)
       Q(i)=real(rand,wp)
    enddo
  end subroutine EP_set_random_interna
  

  subroutine EP_GramSchmidt(i,n)
    implicit none
    integer, intent(in) :: i,n
    ! ::::::::::::::::::::::::::::::

    if( i.ge.0 ) then
       call EP_GramSchmidt_interna(Qvect(1:,i),n)
    else
       call EP_GramSchmidt_interna(dumQvect(1:,-i),n)
    endif
    return
  end subroutine EP_GramSchmidt

  subroutine EP_GramSchmidt_interna(Q,n)
    implicit none
    integer, intent(in) :: n
    real(wp) Q(EP_dim)
    ! ::::::::::::::::::::::::::::::
    integer i,k, volta
    real(wp) scals(0:n-1), scalstot(0:n-1), occscals(1:EP_norb), occscalstot(1:EP_norb)

    if(ha%iproc.eq.0) then
       print *, " inizio GramSchmidt ", n
    endif

    do volta=1,2
       call gemm('T','N', n , 1, EP_dim ,1.0_wp ,  Qvect(1,0)  , EP_dim ,Q(1) ,EP_dim , 0.0_wp , scals(0) ,  n )
       if(ha%nproc/=1) then
          call MPI_Allreduce(scals(0) ,scalstot(0) , n ,mpidtypw, MPI_SUM,MPI_COMM_WORLD ,ierr )
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
             call MPI_Allreduce(occscals(1) ,occscalstot(1) , EP_norb  ,mpidtypw, MPI_SUM,MPI_COMM_WORLD ,ierr )
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
       print *, " finito  GramSchmidt "
    endif

    return
  end subroutine EP_GramSchmidt_interna
 


  subroutine EP_Moltiplica(p,i)
    use module_interfaces
    implicit none
    integer, intent(in) :: p,i
! :::::::::::::::::::::::::::
    integer k, nprocbidon
    real(8) sumdum
 



   if(p.ge.0) then
       print *, " problema Moltiplica solo con dumvect per result "
       stop
    endif
    if(i.lt.0) then
       print *, " problema Moltiplica solo con vect per input "
       stop
    endif
    
    if( ha%nproc > 1) then
       call untranspose_v(ha%iproc,ha%nproc,ha%orbs,ha%lr%wfd,ha%comms,&
            Qvect(1:,i), work=wrk,outadd= Qvect_tmp(1) )  
    else
       do k=1, EP_dim_tot
          Qvect_tmp(k)=  Qvect(k,i)
       enddo
    endif

 

    if(  ha%iproc ==0 ) print *, "chiamo hamiltonian "
    call HamiltonianApplication(ha%iproc,ha%nproc,ha%at,ha%orbs,ha%hx,ha%hy,ha%hz,&
         ha%rxyz,ha%cpmult,ha%fpmult,ha%radii_cf,&
         ha%nlpspd,ha%proj,ha%lr,ha%ngatherarr,            &
         ha%ndimpot, &
         ha%potential,  Qvect_tmp    ,  wrk   ,ha%ekin_sum,ha%epot_sum,ha%eproj_sum,1,ha%GPU)
    if(  ha%iproc ==0 ) print *, " done "

   if(ha%iproc.eq.0) then
      if(EP_shift /= 0 ) then
         call axpy(EP_dim_tot, EP_shift  ,  Qvect_tmp(1)   , 1,  wrk(1) , 1)
      endif
   endif

   if(  ha%nproc/=1) then
      call transpose_v(ha%iproc,ha%nproc,ha%orbs,ha%lr%wfd,ha%comms,&
           wrk , work= Qvect_tmp ,outadd=dumQvect(1,-p))  
   else
      do k=1, EP_dim_tot
         dumQvect(k,-p) =  wrk(k)
      enddo
   endif
   
   return 
 end subroutine EP_Moltiplica
 
 
 
 subroutine EP_add_from_vect_with_fact(p,i, fact)
   implicit none
   integer, intent(in) :: p,i
   real(8), intent(in) :: fact
   ! :::::::::::::::::::::::::::
   integer k
   
   if(p.ge.0) then
      print *, " problema Moltiplica solo con dumvect per result "
      stop
   endif
   if(i.lt.0) then
      print *, " problema Moltiplica solo con vect per input "
      stop
   endif
   
   call axpy(EP_dim, fact,   Qvect(1,i)  , 1,  dumQvect(1,-p)  , 1)
   
!!$    do k=1,EP_dim
!!$       dumQvect(k,-p) = dumQvect(k,-p)+fact*Qvect(k,i)
!!$    enddo
   
   
   return
 end subroutine EP_add_from_vect_with_fact
 
 
 
 subroutine gaussians_to_wavelets_nonorm(iproc,nproc,geocode,orbs,grid,hx,hy,hz,wfd,G,wfn_gau,psi)
   use module_base
   use module_types
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
   logical :: myorbital,maycalc
   integer :: i_stat,i_all,ishell,iexpo,icoeff,iat,isat,ng,l,m,iorb,jorb,i,nterm,ierr,ig,ispinor
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
                    G%xp(iexpo),G%psiat(iexpo),&
                    rx,ry,rz,hx,hy,hz,&
                    0,grid%n1,0,grid%n2,0,grid%n3,&
                    grid%nfl1,grid%nfu1,grid%nfl2,grid%nfu2,grid%nfl3,grid%nfu3,  & 
                    wfd%nseg_c,wfd%nvctr_c,wfd%keyg,wfd%keyv,wfd%nseg_f,wfd%nvctr_f,&
                    wfd%keyg(1,wfd%nseg_c+1),wfd%keyv(wfd%nseg_c+1),&
                    tpsi(1),tpsi(wfd%nvctr_c+1))
            end if
            !sum the result inside the orbital wavefunction
            !loop over the orbitals
            do iorb=1,orbs%norb
               if (orbs%isorb < iorb .and. iorb <= orbs%isorb+orbs%norbp) then
                  jorb=iorb-orbs%isorb
                  do ispinor=1,orbs%nspinor
                     call axpy(wfd%nvctr_c+7*wfd%nvctr_f,wfn_gau(icoeff,ispinor,jorb),&
                          tpsi(1),1,psi(1,ispinor,jorb),1)
                  end do
               end if
            end do
            icoeff=icoeff+1
         end do
         iexpo=iexpo+ng
      end do
      if (iproc == 0 .and. verbose > 1) then
         write(*,'(a)',advance='no') &
              repeat('.',(iat*40)/G%nat-((iat-1)*40)/G%nat)
      end if
   end do

   call gaudim_check(iexpo,icoeff,ishell,G%nexpo,G%ncoeff,G%nshltot)

   if (iproc ==0  .and. verbose > 1) write(*,'(1x,a)')'done.'
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
         totnorm=1.0_wp
         do ispinor=1,orbs%nspinor !to be verified in case of nspinor=4
            call wscal(wfd%nvctr_c,wfd%nvctr_f,real(1.0_dp/sqrt(totnorm),wp),psi(1,ispinor,jorb),psi(wfd%nvctr_c+1,ispinor,jorb))
         end do
         !write(*,'(1x,a,i5,1pe14.7)')'norm of orbital ',iorb,totnorm
         tt=max(tt,abs(1.0_dp-totnorm))
      end if
   end do
   if (nproc > 1) then
      call MPI_REDUCE(tt,normdev,1,mpidtypd,MPI_MAX,0,MPI_COMM_WORLD,ierr)
   else
      normdev=tt
   end if
   if (iproc ==0) write(*,'(1x,a,1pe12.2)')&
        'Deviation from normalization of the imported orbitals',normdev

   i_all=-product(shape(tpsi))*kind(tpsi)
   deallocate(tpsi,stat=i_stat)
   call memocc(i_stat,i_all,'tpsi',subname)

 end subroutine gaussians_to_wavelets_nonorm



 
 
end module lanczos_interface
