module lanczos_interface
  use module_base
  use module_types
  implicit none

  private

  !calculate the allocation dimensions
  public :: inizializza,EP_initialize_start,EP_allocate_for_eigenprob

 
  character(len=*), parameter :: subname='lanczos_interface'
  integer :: i_all,i_stat,ierr
  type(lanczos_args), pointer :: ha

  real(8), pointer :: matrix(:,:)
  real(wp), pointer :: Qvect   (:,:)
  real(wp), pointer :: dumQvect(:,:)

  real(8), pointer :: passed_matrix(:,:)

  real(8) EP_shift
  integer EP_dim
  integer EP_dim_tot

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

  subroutine inizializza(ha_actual)
    implicit none
    type(lanczos_args), target :: ha_actual
    
    ha=>ha_actual

    EP_dim=ha%comms%nvctr_par(ha%iproc)*ha%orbs%nspinor

    EP_dim_tot=(ha%lr%wfd%nvctr_c+7*ha%lr%wfd%nvctr_f)*ha%orbs%nspinor

    if( (ha%lr%wfd%nvctr_c+7*ha%lr%wfd%nvctr_f) /= &
         sum(ha%comms%nvctr_par(:))  ) then
       stop "array size inconsistency" 
    endif

  end subroutine inizializza
  

  subroutine EP_mat_mult( m,k ,  EV   ) ! usare i dumvectors a partire da -1
    implicit none
    integer , intent(in):: m,k
    real(8), intent (in) :: EV(0:m-1,0:k-1)
! ::::::::::::::::::::::::::::::::::::
    integer i,j,l

    stop 'matmult'

    do i=0, k-1
       do j=1, EP_dim
          dumQvect(j,i+1)=0.0d0
          do l=0, m-1 ! controllare il range
             dumQvect(j,i+1)=dumQvect(j,i+1)+EV( l,i )*Qvect(  j,l ) 
          enddo 
       enddo
    enddo


    return 
  end subroutine EP_mat_mult
  

  !allocate the wavefunctions in the transposed form, for lancsoz
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

    call MPI_Allreduce(sump,sumtot,1,mpidtypw, MPI_SUM,MPI_COMM_WORLD ,ierr )

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

    call MPI_Allreduce(sump,sumtot,1,mpidtypw, MPI_SUM,MPI_COMM_WORLD ,ierr )

    EP_scalare_interna=sumtot

  end function EP_scalare_interna
    
  subroutine EP_initialize_start(Gabsorber)
    integer :: i
    type(gaussian_basis) :: Gabsorber
    if (.not. associated(Qvect)) then
       write(*,*)'ERROR: initialization vector not allocated!'
       stop
    end if

    call gaussians_to_wavelets(ha%iproc,ha%nproc,ha%lr%geocode,ha%orbs,ha%lr%d,&
         ha%hx,ha%hy,ha%hz,ha%lr%wfd,Gabsorber,ha%Gabs_coeffs,Qvect(1,0))

    call plot_wf('test',ha%lr,ha%hx,ha%hy,ha%hz,ha%rxyz(1,1),ha%rxyz(2,1),ha%rxyz(3,1),&
         Qvect,'')

    ! plot Qvect(1,0)

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


    return 
  end subroutine EP_set_all_random


  subroutine EP_set_random_interna(Q)
    implicit none
    real(8) Q(EP_dim)

    stop 'random interna'
!!$    ! :::::::::::::::::::::::
!!$    real(8) drand
!!$    external drand
!!$    integer i
!!$    do i=1,EP_dim
!!$       Q(i)=drand(0)
!!$       Q(i)=(-1+10.0D0/i)**i
!!$
!!$    enddo
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
    real(8) Q(EP_dim)
    ! ::::::::::::::::::::::::::::::
    integer i,k, volta
    real(8) scal



    do volta=1,4       
       do i=0, n-1
          scal=0.0d0
          do k=1, EP_dim
             scal=scal+ Q(k)*Qvect(k,i)
          enddo         

          do k=1, EP_dim
             Q(k)=Q(k)-scal*Qvect(k,i)

          enddo
       enddo
    enddo
    return
  end subroutine EP_GramSchmidt_interna



  subroutine EP_Moltiplica(p,i)
    implicit none
    integer, intent(in) :: p,i
! :::::::::::::::::::::::::::
    integer k
    real(8) alpha, beta, sum
    if(p.ge.0) then
       print *, " problema Moltiplica solo con dumvect per result "
       stop
    endif
    if(i.lt.0) then
       print *, " problema Moltiplica solo con vect per input "
       stop
    endif
    
    alpha=1.0
    beta=0.0d0



    call dgemm('N','N', EP_dim , 1 ,EP_dim,alpha, matrix ,EP_dim,   Qvect(1,i) ,EP_dim,beta, dumQvect(1,-p) , EP_dim)
    


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
    do k=1,EP_dim
       dumQvect(k,-p) = dumQvect(k,-p)+fact*Qvect(k,i)
    enddo
    return
  end subroutine EP_add_from_vect_with_fact


end module lanczos_interface
