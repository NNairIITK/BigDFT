

module lanczos_base
  use module_base
  implicit none
  private

  public ::LB_alpha , LB_beta, LB_eval, LB_shift, LB_nsteps, LB_allocate_for_lanczos, &
       LB_de_allocate_for_lanczos, LB_cerca, LB_passeggia, LB_allocate_for_chebychev,&
        LB_iproc, LB_nproc, LB_passeggia_Chebychev
  
  character(len=*), parameter :: subname='lanczos_base'


  real(8), pointer :: LB_alpha(:)
  real(8), pointer :: LB_beta(:)
  real(8), pointer :: omega(:,:)
  real(8), pointer :: evect(:,:)
  real(8), pointer :: LB_eval(:)
  real(8), pointer :: oldalpha(:)
  real(8), pointer :: diagwork(:)
  real(8) LB_shift
  real(8) LANCZOS_tol
  integer LB_nsteps
  integer hasoldalpha
  integer i_stat,i_all
  integer LB_iproc, LB_nproc
contains

  subroutine LB_allocate_for_chebychev( )
   if(associated(LB_alpha )) then
      call  LB_de_allocate_for_lanczos( )
   endif

   allocate(LB_alpha(0: 3*LB_nsteps ) , stat=i_stat)
   call memocc(i_stat,LB_alpha,'LB_alpha',subname)

   allocate(LB_beta(0: 1 ) , stat=i_stat)
   call memocc(i_stat,LB_beta,'LB_beta',subname)

   allocate(omega( 0:1,  0:1      )  , stat=i_stat)
   call memocc(i_stat,omega,'omega',subname)

   
   allocate(evect( 0:1,  0:1    )   , stat=i_stat)
   call memocc(i_stat,evect,'evect',subname)

   allocate(LB_eval(0: 1 )   , stat=i_stat)
   call memocc(i_stat,LB_eval,'LB_eval',subname)

   allocate(diagwork( 0:1*(3+1)   )  , stat=i_stat)
   call memocc(i_stat,diagwork,'diagwork',subname)

   allocate(oldalpha (0: 1 )        , stat=i_stat  )
   call memocc(i_stat,oldalpha,'oldalpha',subname)

  end subroutine LB_allocate_for_chebychev


  subroutine LB_allocate_for_lanczos( )
   if(associated(LB_alpha )) then
      call  LB_de_allocate_for_lanczos( )
   endif

   allocate(LB_alpha(0: LB_nsteps ) , stat=i_stat)
   call memocc(i_stat,LB_alpha,'LB_alpha',subname)

   allocate(LB_beta(0: LB_nsteps-1 ) , stat=i_stat)
   call memocc(i_stat,LB_beta,'LB_beta',subname)

   allocate(omega( 0:LB_nsteps,  0:LB_nsteps      )  , stat=i_stat)
   call memocc(i_stat,omega,'omega',subname)

   
   allocate(evect( 0:LB_nsteps-1,  0:LB_nsteps-1    )   , stat=i_stat)
   call memocc(i_stat,evect,'evect',subname)

   allocate(LB_eval(0: LB_nsteps-1 )   , stat=i_stat)
   call memocc(i_stat,LB_eval,'LB_eval',subname)

   allocate(diagwork( 0:LB_nsteps*(3+LB_nsteps)   )  , stat=i_stat)
   call memocc(i_stat,diagwork,'diagwork',subname)

   allocate(oldalpha (0: LB_nsteps )        , stat=i_stat  )
   call memocc(i_stat,oldalpha,'oldalpha',subname)

   

   omega(:,:)=0.0D0

 end subroutine LB_allocate_for_lanczos

  subroutine LB_de_allocate_for_lanczos( )

    i_all=-product(shape(LB_alpha))*kind(LB_alpha)
    deallocate(LB_alpha)
    call memocc(i_stat,i_all,'LB_alpha',subname)


    i_all=-product(shape(LB_beta))*kind(LB_beta)
    deallocate(LB_beta)
    call memocc(i_stat,i_all,'LB_beta',subname)


    i_all=-product(shape(omega))*kind(omega)
    deallocate(omega )
    call memocc(i_stat,i_all,'omega',subname)

    i_all=-product(shape(evect))*kind(evect)
    deallocate(evect  )
    call memocc(i_stat,i_all,'evect',subname)

    i_all=-product(shape(LB_eval))*kind(LB_eval)
    deallocate(LB_eval  )
    call memocc(i_stat,i_all,'LB_eval',subname)

    i_all=-product(shape(diagwork))*kind(diagwork)
    deallocate(diagwork    )
    call memocc(i_stat,i_all,'diagwork',subname)

    i_all=-product(shape(oldalpha))*kind(oldalpha)
    deallocate(oldalpha          )
    call memocc(i_stat,i_all,'oldalpha',subname)

  end subroutine LB_de_allocate_for_lanczos
   
  integer function converged(m)
    implicit none
    integer, intent(in):: m
! ::::::::::::::::::::::::::::::::::::::::::
    integer j,i,k
    real(8) dum, dumvect(0:LB_nsteps )
    do j=1, m-1
       do i=0, m-j-1
          if( abs(LB_eval(i) ) .lt. abs(LB_eval(i+1) ) ) then
             dum = LB_eval(i)
             LB_eval(i)= LB_eval(i+1)
             LB_eval(i+1)=dum
             do k=0, LB_nsteps-1
                dumvect(k)   = evect(k,i)
                evect(k,  i) = evect(k,i+1)
                evect(k,i+1)  = dumvect(k)
             enddo
          endif
       enddo
    enddo
    

    converged=0
    
    if( hasoldalpha.eq.1) then
       do while(converged.lt.m)
          if(LB_iproc==0) print *, LB_eval(converged),oldalpha(converged),  LANCZOS_tol
          if (abs(LB_eval(converged)-oldalpha(converged))/abs(oldalpha(converged)-LB_shift)  .gt. LANCZOS_tol) then
             exit
          endif
          converged=converged+1
       enddo
    else
       hasoldalpha=1
    endif

   
    oldalpha(:m-1)=LB_eval(:m-1)

  end function converged



  subroutine diago(k,m)
    implicit none
    integer, intent(in):: k,m
! :::::::::::::::::::::::::::::::::::::::::::
    integer i,j
    integer info, lwork
    
    evect(:,:)=0.0D0
    do i=0, m-1
       evect(i,i)=LB_alpha(i)
    enddo
    do i=k,m-2
       evect(i,i+1)=LB_beta(i)
    enddo
    evect(0:k-1,k)=LB_beta(0:k-1)


!!$    if(LB_iproc==0) then
!!$       print *, " in diag "
!!$       print *, "  evect  " , evect
!!$       print *, " LB_eval  " , LB_eval
!!$    endif
    lwork = LB_nsteps*(3+LB_nsteps)

    call DSYEV( 'V', 'U', m, evect(0,0) , m, LB_eval(0), diagwork(0), LWORK, INFO )

    if(LB_iproc==0) then
       print *, " evals  " 
       print *, LB_eval-LB_shift
    endif

!!$    if(LB_nproc>1) call MPI_finalize()
!!$    stop

    if(info .ne.0) then
       print *, " problem with dsyev"
       stop
    endif
    return 
  end subroutine diago

  integer function  LB_cerca( nd, shift, tol, set_EP_shift, EP_allocate_for_eigenprob,&
       EP_make_dummy_vectors, get_EP_dim, EP_initialize_start , EP_normalizza,&
       EP_Moltiplica, EP_GramSchmidt ,EP_set_all_random, EP_copy , EP_mat_mult, &
       EP_scalare ,EP_add_from_vect_with_fact, accontentati_di_n)

    implicit none

    integer, intent(IN) :: nd
    real(8), intent(IN) :: shift
    real(8), intent(IN) :: tol
    interface
       subroutine set_EP_shift(shift)
         real(8) shift
       end subroutine set_EP_shift

       subroutine  EP_allocate_for_eigenprob(nsteps)
         implicit none
         integer, intent(in):: nsteps
       end subroutine
       subroutine EP_make_dummy_vectors(nd)
         implicit none
         integer, intent(in):: nd
         ! :::::::::::::::::::::::    
       end subroutine
       integer function get_EP_dim()
       end function

       subroutine EP_initialize_start()
       end subroutine 
       subroutine EP_normalizza(i)
         integer i
       end subroutine 
       subroutine EP_Moltiplica(i,j)
         integer i,j
       end subroutine 
       real(8) function EP_scalare(i,j)
         integer i,j
       end function 
       subroutine EP_add_from_vect_with_fact( i, j  ,   a )
         integer i,j
         real(8) a
       end subroutine 
       subroutine EP_GramSchmidt(i,j)
         integer i,j
       end subroutine 
       subroutine EP_set_all_random(i)
         integer i
       end subroutine 
       subroutine EP_copy(i,j)
         integer i,j
       end subroutine
       subroutine EP_mat_mult(m,k ,  EV )
         integer  m,k
         real(8)  EV(1 )
       end subroutine        

    end interface

    integer , optional :: accontentati_di_n

! :::::::::::::::::::::::::::::::::::::::::::::::::::
    integer k,nc,m
    integer n_converged_cercati
    real(8) dum

    LANCZOS_tol=tol
    call set_EP_shift(shift)
    
    LB_shift=shift
    
    k = get_EP_dim()

    m = min(4*nd, get_EP_dim())

    LB_nsteps = m
    hasoldalpha=0

    call LB_allocate_for_lanczos( )
    call EP_allocate_for_eigenprob(LB_nsteps)
    call EP_make_dummy_vectors(LB_nsteps)



    k=0
    nc=0
 

    call LB_passeggia(k,m,      get_EP_dim, EP_initialize_start , EP_normalizza,&
         EP_Moltiplica, EP_GramSchmidt ,EP_set_all_random, EP_copy,   EP_mat_mult, &
         EP_scalare,EP_add_from_vect_with_fact     )

    n_converged_cercati = nd
    if(present(accontentati_di_n)) then
       n_converged_cercati=accontentati_di_n
    endif

    do while(nc.lt.n_converged_cercati)

       call diago(k,m)



       nc=converged(m)

       if ( k.gt.0) then
          if( notzero( k, LANCZOS_tol).eq.0 ) then
             exit
          endif
       endif


       if ( (nc+2*nd) .ge. m) then 
          k=m-1
       else
          k=nc+2*nd
       endif
  

       

       call ricipolla(k,m, EP_normalizza ,EP_copy, EP_mat_mult)


       call LB_passeggia(k,m , get_EP_dim, EP_initialize_start , EP_normalizza, &
            EP_Moltiplica, EP_GramSchmidt ,EP_set_all_random, EP_copy ,&
            EP_mat_mult,   EP_scalare,EP_add_from_vect_with_fact)
    enddo
    if( m.eq.get_EP_dim()) then
       LB_cerca=m
    else
       LB_cerca= k
    endif
    return 
  end function LB_cerca

  integer function notzero(k, tol)
    implicit none
    integer, intent(in) :: k
    real(8), intent(in) :: tol
    !::::::::::::::::::::::::::::::`
    integer i
    notzero=0
    do i=0,k-1
       if( abs(LB_beta(i)).gt.tol) then
          notzero=notzero+1
       endif
    enddo
    return 
  end function notzero


  subroutine ricipolla(k,m, EP_normalizza ,EP_copy, EP_mat_mult )
    implicit none
    integer, intent(in):: k,m
    interface
       subroutine EP_normalizza(i)
         integer i
       end subroutine EP_normalizza
       subroutine EP_copy(i,j)
         integer i,j
       end subroutine EP_copy
       subroutine EP_mat_mult(m,k ,  EV )
         integer  m,k
         real(8)  EV(1 )
       end subroutine EP_mat_mult
    end interface
    
! :::::::::::::::::::::::::::::::::::::
    integer i
    real(8), pointer :: dumomega(:,:), dumomega2(:,:)
    real(8) acoeff, bcoeff
    integer lda, ldb, ldc




    LB_alpha(0:k-1)=LB_eval(0:k-1)

    LB_beta(0:k-1) = LB_beta(m-1)*evect(m-1, 0:k-1)



    call EP_mat_mult(m,k,  evect(0:, 0:  )   ) 


    do i=0, k-1
       call EP_normalizza(-i-1)
       
       call EP_copy(i, -i-1 )
    enddo


    call EP_copy(k, m )
    
  
    allocate(dumomega(0:m-1,0:m-1))
    allocate(dumomega2(0:m-1,0:m-1))
    

    acoeff=1.0D0
    bcoeff =0.0D0
    


    call dgemm('N','N',m,m,m,acoeff,omega(0,0) ,LB_nsteps+1,  evect(0,0) ,LB_nsteps,bcoeff,dumomega(0,0),m)

    omega(0:k,k)=dumomega(0:k,k)
    omega(k,0:k)=dumomega(0:k,k)

    call dgemm('T','N',m,m,m,acoeff,evect(0,0) ,LB_nsteps,  dumomega (0,0) ,m,bcoeff,dumomega2(0,0),m)
    omega(0:k-1,0:k-1)= dumomega2 (0:k-1,0:k-1)


    deallocate(dumomega)
    deallocate( dumomega2)
    

  end subroutine ricipolla



  subroutine LB_passeggia( k, m, get_EP_dim, EP_initialize_start , EP_normalizza, &
       EP_Moltiplica, EP_GramSchmidt ,EP_set_all_random, EP_copy,  EP_mat_mult,&
       EP_scalare,EP_add_from_vect_with_fact )

    implicit none
    integer, intent(in)::k
    integer, intent(in)::m
    ! ::::::::::::::::::::::::::::::::::::::::
    real(8) sn,eu,eusn,eq
    real(8) max0
    integer pvect
    integer i,l,j, w
    interface
       subroutine EP_mat_mult(m,k ,  EV )
         integer  m,k
         real(8)  EV(1 )
       end subroutine  
       integer function get_EP_dim()
       end function
       subroutine EP_initialize_start()
       end subroutine 
       subroutine EP_normalizza(i)
         integer i
       end subroutine 
       subroutine EP_Moltiplica(i,j)
         integer i,j
       end subroutine 
       real(8) function EP_scalare(i,j)
         integer i,j
       end function 
       subroutine EP_add_from_vect_with_fact( i, j  ,   a )
         integer i,j
         real(8) a
       end subroutine 
       subroutine EP_GramSchmidt(i,j)
         integer i,j
       end subroutine 
       subroutine EP_set_all_random(i)
         integer i
       end subroutine 
       subroutine EP_copy(i,j)
         integer i,j
       end subroutine
    end interface
    integer p
    real(8) add
    real(8) condition
    integer ipa

    p =-1
        
    sn   = sqrt(get_EP_dim()*1.0D0 )
    eu   = 1.1D-15
    eusn = eu*sn
    eq   = sqrt(eu)

    if (k.eq.0) then
       call EP_initialize_start()
       call EP_normalizza(0)
    endif
    
    DO i=k, m-1

       if(LB_iproc==0 .and. modulo( m-1-i   ,100)==0 ) then
          open(unit=22,file="alphabeta_p")
          write(22,*) i, 1.0
          do ipa=0, i-1
             write(22,*) LB_alpha(ipa), LB_beta(ipa)
          enddo
          close(unit=22)
       end if

       call EP_Moltiplica(p, i )
       
       LB_alpha(i) = EP_scalare(p , i)
       
       if (i.eq.k) then

          call EP_add_from_vect_with_fact( p, k  ,   -LB_alpha(k) )
          do l=0,k-1

             call EP_add_from_vect_with_fact( p, l  ,   -LB_beta(l) )
          enddo
       else

          call EP_add_from_vect_with_fact( p, i    ,   -LB_alpha(i)   )

          call EP_add_from_vect_with_fact( p, i-1  ,   -LB_beta (i-1) )
       endif
       
       

       LB_beta(i)=sqrt(EP_scalare(p,p))


       
       omega(i,i)=1.
       
       max0 = 0.0


       
       if( LB_beta(i).ne.0.0) then
          do j=0,i
             omega(i+1,j) = eusn

             if( j.lt.k) then
                


                add = 2 * eusn + abs( LB_alpha(j)-LB_alpha(i)) * abs(omega(i,j))
                
                if(i.ne.k) then
                   add = add+  LB_beta(j)*abs(omega(i,k))
                endif
                if ( i.gt.0 .and. j.ne.(i-1)) then
                   add = add +  LB_beta(i-1)*abs(omega(i-1,j) )
                endif
                omega(i+1,j) = omega(i+1,j)+ add / LB_beta(i)
                
                
             else if (j.eq.k)  then                        
                add = 2 * eusn + abs(LB_alpha(j)-LB_alpha(i))* abs( omega(i,j) )
                do w=0,k-1 
                   add = add + LB_beta(w)* abs( omega(i,w) )
                enddo
                if (i.ne.(k+1)) then
                   add  =  add +  LB_beta(k)*abs( omega(i,k+1) )
                endif
                
                if(  i.gt.0 .and. i.ne.(k+1)) then
                   add = add +  LB_beta(i-1)* abs(omega(i-1,k))
                endif
                
                omega(i+1,j)  = omega(i+1,j) + add / LB_beta(i)
                
             else if( j.lt.i) then 
                
                add = 2 * eusn + abs(LB_alpha(j)- LB_alpha(i))  * abs(omega(i,j) )
                
                if( i.ne.(j+1)) then
                   add =  add +  LB_beta(j) * abs( omega(i,j+1) ) 
                endif
                if (i.gt.0 .and. j.gt.0) then
                   add = add +  LB_beta(j-1)*abs( omega(i-1,j-1))
                endif
                if(  i.gt.0 .and. i.ne.(j+1)) then
                   add = add +   LB_beta(i-1)*abs(omega(i-1,j))
                endif
                omega(i+1,j)  = omega(i+1,j) +  add / LB_beta(i)
                
             else
                
                add = eusn
                
                if (i.gt.0) then 
                   add = add +  LB_beta(i-1)*abs( omega(i,i-1))
                endif
                
                omega(i+1,j)  = omega(i+1,j) +  add / LB_beta(i)
             endif
             


             omega(j,i+1) = omega(i+1,j)
             
             max0 = max0+ omega(i+1,j)**2
          enddo
       endif
       


       
       if ( LB_beta(i).eq.0.0 .or.  max0.gt.eu  ) then

          if (i.gt.0) then

             call EP_GramSchmidt(i,i)
             

             call EP_normalizza(i)
             

             call EP_Moltiplica(p, i)
             
             LB_alpha(i) = EP_scalare(p , i)
          endif
          

          call EP_GramSchmidt(p,i+1)
          
          LB_beta(i) = sqrt(EP_scalare(p,p))
          

          call EP_normalizza(p)
          if (i.gt.0) then
             condition = eu * sqrt(get_EP_dim() * ( LB_alpha(i)**2+ LB_beta(i-1)**2))
          else
             condition = eu * sqrt(get_EP_dim() * (LB_alpha(i)**2))
          endif
          
          if ( LB_beta(i).lt. condition) then
             
             if(LB_iproc==0) print *, "starting with a perpendicular vector"
             
             LB_beta(i)=0.0D0
             
             call EP_set_all_random(p)
             

             call EP_GramSchmidt(p,i+1)


             call EP_normalizza(p)

             
          endif
          
          do l=0,i-1
             omega(i,l)=eusn
             omega(l,i)=eusn
          enddo


          do l=0,i
             omega(i+1,l)=eusn
             omega(l,i+1)=eusn
          enddo


       else

          call EP_normalizza(p)
       endif

       call EP_copy(i+1, p)
       ! emergenza
!!!       call EP_GramSchmidt(i+1,i+1)
!!!       call EP_normalizza(i+1)


       
!!!
!!!
!!!       print *, "LB_beta(",i,")=" , LB_beta(i)
!!!       ! controlla hermitianicita
!!!       call EP_Moltiplica(p, i+1 )
!!!       print *, "alla rovescio " , EP_scalare(p,i)- LB_beta(i)
!!!       
!!!       if (i.gt.2) then
!!!          print *, "con due prima " , EP_scalare(i-1,i+1)
!!!
!!!       endif


    enddo
    
  end subroutine LB_passeggia


  subroutine CalcolaSpettroChebychev( cheb_shift, fact_cheb,   Nu ) 
    real(gp) cheb_shift, fact_cheb
    integer Nu

    real(gp), pointer :: Xs(:), res(:), cfftreal(:), cfftimag(:)
    complex(kind=8), pointer :: alphas(:), expn(:)
    complex(kind=8) ctemp

    integer Nbar
    integer i,n
    real(gp) Pi
    character(200)  filename, saux
    real(gp) fact



    write(filename,'(a,i0)') "cheb_spectra_" , Nu
    print *, " writing spectra to " , filename 

    Pi=acos(-1.0)
    print *, Pi
    Nbar =1
    do while(Nbar<Nu) 
       Nbar=Nbar*2
    enddo
    Nbar=Nbar*2   


    
    allocate(Xs(0:Nbar-1) , stat=i_stat)
    call memocc(i_stat,Xs,'Xs',subname)
    
    allocate(res(0:Nbar-1) , stat=i_stat)
    call memocc(i_stat,res,'res',subname)
        
    allocate(alphas(0:Nbar-1) , stat=i_stat)
    !! call memocc(i_stat,alphas,'alphas',subname)
    
    allocate(expn(0:Nbar-1) , stat=i_stat)


    allocate(cfftreal(0:2*Nbar-1) , stat=i_stat)
    call memocc(i_stat,cfftreal,'cfftreal',subname)
  
    allocate(cfftimag(0:2*Nbar-1) , stat=i_stat)
    call memocc(i_stat,cfftimag,'cfftimag',subname)
    


    cfftreal=0
    cfftimag=0
    do i=0,Nbar-1
       Xs(i)=cos( (Nbar-1 - i +0.5_gp)*Pi/Nbar)
!!$       res(i)=LB_alpha(0)
!!$       alphas(i) = exp(   ((0.0_gp,1.0_gp) * ( Nbar-1 - i +0.5_gp)) *Pi/Nbar     ) 
!!$       expn(i)=(1.0_gp,0)
!!$

!!$
!!$       exp(   ((0.0_gp,1.0_gp) * (  - i )) *Pi/Nbar     ) 
!!$       exp(   ((0.0_gp,1.0_gp) * ( Nbar-1  +0.5_gp)) *Pi/Nbar     ) = -exp(   ((0.0_gp,1.0_gp) * ( -0.5_gp)) *Pi/Nbar     )
       
       if(i<Nu) then
          fact = 1.0/(Nu+1)*( (Nu-i +1.0)*cos(pi* i /(Nu+1.0))+ &
               sin(pi* i /(Nu+1.0))* cos(pi/(Nu+1))/sin(pi/(Nu+1))    )
          ctemp= LB_alpha(i)*fact*exp(   ((0.0_gp,1.0_gp) * ( Nbar -0.5_gp)) *i*Pi/Nbar     )
          cfftreal(i)=REAL(ctemp)
          cfftimag(i)=AIMAG(ctemp)
       else
          cfftreal(i)=0
          cfftimag(i)=0
       endif
    enddo


    call FFT842(0, 2*Nbar, cfftreal(0:), cfftimag(0:) )
    cfftreal=cfftreal-LB_alpha(0)*0.5

!!$    do n=1,Nu-1
!!$       expn(:)=expn(:)*alphas(:)
!!$       fact = 1.0/(Nu+1)*( (Nu-n +1.0)*cos(pi* n /(Nu+1))+ &
!!$       sin(pi* n /(Nu+1.0))* cos(pi/(Nu+1))/sin(pi/(Nu+1))    )
!!$       res(:)=res(:)+2*REAL(expn(:))*LB_alpha(n)*fact
!!$    enddo

    print *, " done " 
    
!!$    res =res/Pi/sqrt(1-Xs*Xs)
    cfftreal(0:Nbar-1) =2*cfftreal(0:Nbar-1)/Pi/sqrt(1-Xs*Xs)
    
    
    open(unit=22,file=filename)
    do i=0, Nbar-1
       write(22,*) Xs(i) / fact_cheb + cheb_shift , cfftreal(i)
    enddo
    close(unit=22)

    i_all=-product(shape(Xs))*kind(Xs)
    deallocate(Xs)
    call memocc(i_stat,i_all,'Xs',subname)
    
    i_all=-product(shape(res))*kind(res)
    deallocate(res)
    call memocc(i_stat,i_all,'res',subname)
    
    i_all=-product(shape(alphas))*kind(alphas)
    deallocate(alphas)
    !! call memocc(i_stat,i_all,'alphas',subname)
    
    i_all=-product(shape(expn))*kind(expn)
    deallocate(expn)
    !! call memocc(i_stat,i_all,'expn',subname)

    i_all=-product(shape(cfftreal))*kind(cfftreal)
    deallocate(cfftreal)
    call memocc(i_stat,i_all,'cfftreal',subname)
    
    i_all=-product(shape(cfftimag))*kind(cfftimag)
    deallocate(cfftimag)
    call memocc(i_stat,i_all,'cfftimag',subname)
    



  end subroutine CalcolaSpettroChebychev
  
  
  subroutine LB_passeggia_Chebychev (  m, cheb_shift,  fact_cheb,  get_EP_dim,  EP_initialize_start , EP_normalizza, &
       EP_Moltiplica, EP_GramSchmidt ,EP_set_all_random, EP_copy,  EP_mat_mult,&
       EP_scalare,EP_add_from_vect_with_fact , EP_multbyfact)
    use module_base
    implicit none
    integer, intent(in)::m
    ! ::::::::::::::::::::::::::::::::::::::::
    real(8) sn,eu,eusn,eq
    real(8) max0
    real(gp) cheb_shift,  fact_cheb
    integer pvect
    integer i,l,j, w
    interface
       subroutine EP_mat_mult(m,k ,  EV )
         integer  m,k
         real(8)  EV(1 )
       end subroutine  
       integer function get_EP_dim()
       end function
       subroutine EP_initialize_start()
       end subroutine 
       subroutine EP_normalizza(i)
         integer i
       end subroutine 
       subroutine EP_Moltiplica(i,j)
         integer i,j
       end subroutine 
       real(8) function EP_scalare(i,j)
         integer i,j
       end function 
       subroutine EP_add_from_vect_with_fact( i, j  ,   a )
         integer i,j
         real(8) a
       end subroutine 
       subroutine EP_GramSchmidt(i,j)
         integer i,j
       end subroutine 
       subroutine EP_set_all_random(i)
         integer i
       end subroutine 
       subroutine EP_copy(i,j)
         integer i,j
       end subroutine
       subroutine EP_multbyfact(j, fact)
         use module_base
         implicit none
         integer, intent(in)::j
         real(gp) fact
         ! ::::::::::::::::::::::::::::::::::::::

       end subroutine EP_multbyfact
       

    end interface
    integer p
    real(8) add
    real(8) condition
    integer ipa
    integer tmp1, attuale, precedente

    tmp1 = 1
    attuale= 2
    precedente= 3

    call EP_initialize_start()
    call EP_copy(attuale,0)
    call EP_copy(precedente,0)
    LB_alpha(0)=EP_scalare(attuale,attuale)
    
    DO i=0, m-1
       if(LB_iproc==0 .and. modulo( m-1-i   ,100)==0 ) then
          open(unit=22,file="coeffs_chebychev")
          write(22,*) 2*i+1,  cheb_shift, fact_cheb
          do ipa=0, 2*i
             write(22,*) LB_alpha(ipa)
          enddo
          close(unit=22)
          call CalcolaSpettroChebychev( cheb_shift, fact_cheb,   2*i+1 ) 
      endif
       
       call EP_Moltiplica(tmp1, attuale )
       call EP_multbyfact(tmp1,fact_cheb)
       if(i==0) then
          LB_alpha(1)=EP_scalare(tmp1,attuale)
       else
          call EP_multbyfact(tmp1,2.0_gp)
          call EP_add_from_vect_with_fact(tmp1,precedente,-1.0_gp ) 
       endif
       
       
       LB_alpha(2*i+1) = 2*EP_scalare(tmp1,attuale)-LB_alpha(1)
       LB_alpha(2*i+2) = 2*EP_scalare(tmp1,tmp1)-LB_alpha(0)
       
       call EP_copy(precedente,attuale)
       call EP_copy(attuale,tmp1)
       
    enddo
    
  end subroutine LB_passeggia_Chebychev




!
!-----------------------------------------------------------------------
! SUBROUTINE:  FFT842
! FAST FOURIER TRANSFORM FOR N=2**M
! COMPLEX INPUT
!-----------------------------------------------------------------------
!
      SUBROUTINE FFT842(IN, N, X, Y)
!
! THIS PROGRAM REPLACES THE VECTOR Z=X+IY BY ITS  FINITE
! DISCRETE, COMPLEX FOURIER TRANSFORM IF IN=0.  THE INVERSE TRANSFORM
! IS CALCULATED FOR IN=1.  IT PERFORMS AS MANY BASE
! 8 ITERATIONS AS POSSIBLE AND THEN FINISHES WITH A BASE 4 ITERATION
! OR A BASE 2 ITERATION IF NEEDED.
!
! THE SUBROUTINE IS CALLED AS SUBROUTINE FFT842 (IN,N,X,Y).
! THE INTEGER N (A POWER OF 2), THE N REAL LOCATION ARRAY X, AND
! THE N REAL LOCATION ARRAY Y MUST BE SUPPLIED TO THE SUBROUTINE.
!

      IMPLICIT REAL*8(A-H, O-Z), INTEGER(I-N)

      DIMENSION X(*), Y(*), L(15)
      COMMON /CON2/ PI2, P7
      EQUIVALENCE (L15,L(1)), (L14,L(2)), (L13,L(3)), (L12,L(4)),&
         (L11,L(5)), (L10,L(6)), (L9,L(7)), (L8,L(8)), (L7,L(9)),&
         (L6,L(10)), (L5,L(11)), (L4,L(12)), (L3,L(13)), (L2,L(14)),&
         (L1,L(15))

      PI2 = 8.*ATAN(1.)
      P7 = 1./SQRT(2.)
      DO 10 I=1,15
        M = I
        NT = 2**I
        IF (N.EQ.NT) GO TO 20
  10  CONTINUE
      WRITE (*,9999)
9999  FORMAT (35H N IS NOT A POWER OF TWO FOR FFT842)
      STOP
  20  N2POW = M
      NTHPO = N
      FN = NTHPO
      IF (IN.EQ.1) GO TO 40
      DO 30 I=1,NTHPO
        Y(I) = -Y(I)
  30  CONTINUE
  40  N8POW = N2POW/3
      IF (N8POW.EQ.0) GO TO 60
!
! RADIX 8 PASSES,IF ANY.
!
      DO 50 IPASS=1,N8POW
        NXTLT = 2**(N2POW-3*IPASS)
        LENGT = 8*NXTLT
        CALL R8TX(NXTLT, NTHPO, LENGT, X(1), X(NXTLT+1), X(2*NXTLT+1),&
           X(3*NXTLT+1), X(4*NXTLT+1), X(5*NXTLT+1), X(6*NXTLT+1),&
           X(7*NXTLT+1), Y(1), Y(NXTLT+1), Y(2*NXTLT+1), Y(3*NXTLT+1),&
           Y(4*NXTLT+1), Y(5*NXTLT+1), Y(6*NXTLT+1), Y(7*NXTLT+1))
  50  CONTINUE
!
! IS THERE A FOUR FACTOR LEFT
!
  60  IF (N2POW-3*N8POW-1) 90, 70, 80
!
! GO THROUGH THE BASE 2 ITERATION
!
!
  70  CALL R2TX(NTHPO, X(1), X(2), Y(1), Y(2))
      GO TO 90
!
! GO THROUGH THE BASE 4 ITERATION
!
  80  CALL R4TX(NTHPO, X(1), X(2), X(3), X(4), Y(1), Y(2), Y(3), Y(4))
!
  90  DO 110 J=1,15
        L(J) = 1
        IF (J-N2POW) 100, 100, 110
 100    L(J) = 2**(N2POW+1-J)
 110  CONTINUE
      IJ = 1
      DO 130 J1=1,L1
      DO 130 J2=J1,L2,L1
      DO 130 J3=J2,L3,L2
      DO 130 J4=J3,L4,L3
      DO 130 J5=J4,L5,L4
      DO 130 J6=J5,L6,L5
      DO 130 J7=J6,L7,L6
      DO 130 J8=J7,L8,L7
      DO 130 J9=J8,L9,L8
      DO 130 J10=J9,L10,L9
      DO 130 J11=J10,L11,L10
      DO 130 J12=J11,L12,L11
      DO 130 J13=J12,L13,L12
      DO 130 J14=J13,L14,L13
      DO 130 JI=J14,L15,L14
        IF (IJ-JI) 120, 130, 130
 120    R = X(IJ)
        X(IJ) = X(JI)
        X(JI) = R
        FI = Y(IJ)
        Y(IJ) = Y(JI)
        Y(JI) = FI
 130    IJ = IJ + 1
      IF (IN.EQ.1) GO TO 150
      DO 140 I=1,NTHPO
        Y(I) = -Y(I)
 140  CONTINUE
      GO TO 170
 150  DO 160 I=1,NTHPO
        X(I) = X(I)/FN
        Y(I) = Y(I)/FN
 160  CONTINUE
 170  RETURN
     END SUBROUTINE FFT842





!
!-----------------------------------------------------------------------
! SUBROUTINE:  R2TX
! RADIX 2 ITERATION SUBROUTINE
!-----------------------------------------------------------------------
!
      SUBROUTINE R2TX(NTHPO, CR0, CR1, CI0, CI1)
      IMPLICIT REAL*8(A-H, O-Z), INTEGER(I-N)
      DIMENSION CR0(*), CR1(*), CI0(*), CI1(*)
      DO 10 K=1,NTHPO,2
        R1 = CR0(K) + CR1(K)
        CR1(K) = CR0(K) - CR1(K)
        CR0(K) = R1
        FI1 = CI0(K) + CI1(K)
        CI1(K) = CI0(K) - CI1(K)
        CI0(K) = FI1
  10  CONTINUE
      RETURN
      END  SUBROUTINE R2TX
!
!-----------------------------------------------------------------------
! SUBROUTINE:  R4TX
! RADIX 4 ITERATION SUBROUTINE
!-----------------------------------------------------------------------
!
      SUBROUTINE R4TX(NTHPO, CR0, CR1, CR2, CR3, CI0, CI1, CI2, CI3)
      IMPLICIT REAL*8(A-H, O-Z), INTEGER(I-N)
      DIMENSION CR0(*), CR1(*), CR2(*), CR3(*), CI0(*), CI1(*), CI2(*),&
         CI3(*)
      DO 10 K=1,NTHPO,4
        R1 = CR0(K) + CR2(K)
        R2 = CR0(K) - CR2(K)
        R3 = CR1(K) + CR3(K)
        R4 = CR1(K) - CR3(K)
        FI1 = CI0(K) + CI2(K)
        FI2 = CI0(K) - CI2(K)
        FI3 = CI1(K) + CI3(K)
        FI4 = CI1(K) - CI3(K)
        CR0(K) = R1 + R3
        CI0(K) = FI1 + FI3
        CR1(K) = R1 - R3
        CI1(K) = FI1 - FI3
        CR2(K) = R2 - FI4
        CI2(K) = FI2 + R4
        CR3(K) = R2 + FI4
        CI3(K) = FI2 - R4
  10  CONTINUE
      RETURN
      END  SUBROUTINE R4TX
!
!-----------------------------------------------------------------------
! SUBROUTINE:  R8TX
! RADIX 8 ITERATION SUBROUTINE
!-----------------------------------------------------------------------
!
      SUBROUTINE R8TX(NXTLT, NTHPO, LENGT, CR0, CR1, CR2, CR3, CR4,&
         CR5, CR6, CR7, CI0, CI1, CI2, CI3, CI4, CI5, CI6, CI7)
      IMPLICIT REAL*8(A-H, O-Z), INTEGER(I-N)
      DIMENSION CR0(*), CR1(*), CR2(*), CR3(*), CR4(*), CR5(*), CR6(*),&
         CR7(*), CI1(*), CI2(*), CI3(*), CI4(*), CI5(*), CI6(*),&
         CI7(*), CI0(*)
      COMMON /CON2/ PI2, P7
!
      SCALE = PI2/FLOAT(LENGT)
      DO 30 J=1,NXTLT
        ARG = FLOAT(J-1)*SCALE
        C1 = COS(ARG)
        S1 = SIN(ARG)
        C2 = C1**2 - S1**2
        S2 = C1*S1 + C1*S1
        C3 = C1*C2 - S1*S2
        S3 = C2*S1 + S2*C1
        C4 = C2**2 - S2**2
        S4 = C2*S2 + C2*S2
        C5 = C2*C3 - S2*S3
        S5 = C3*S2 + S3*C2
        C6 = C3**2 - S3**2
        S6 = C3*S3 + C3*S3
        C7 = C3*C4 - S3*S4
        S7 = C4*S3 + S4*C3
        DO 20 K=J,NTHPO,LENGT
          AR0 = CR0(K) + CR4(K)
          AR1 = CR1(K) + CR5(K)
          AR2 = CR2(K) + CR6(K)
          AR3 = CR3(K) + CR7(K)
          AR4 = CR0(K) - CR4(K)
          AR5 = CR1(K) - CR5(K)
          AR6 = CR2(K) - CR6(K)
          AR7 = CR3(K) - CR7(K)
          AI0 = CI0(K) + CI4(K)
          AI1 = CI1(K) + CI5(K)
          AI2 = CI2(K) + CI6(K)
          AI3 = CI3(K) + CI7(K)
          AI4 = CI0(K) - CI4(K)
          AI5 = CI1(K) - CI5(K)
          AI6 = CI2(K) - CI6(K)
          AI7 = CI3(K) - CI7(K)
          BR0 = AR0 + AR2
          BR1 = AR1 + AR3
          BR2 = AR0 - AR2
          BR3 = AR1 - AR3
          BR4 = AR4 - AI6
          BR5 = AR5 - AI7
          BR6 = AR4 + AI6
          BR7 = AR5 + AI7
          BI0 = AI0 + AI2
          BI1 = AI1 + AI3
          BI2 = AI0 - AI2
          BI3 = AI1 - AI3
          BI4 = AI4 + AR6
          BI5 = AI5 + AR7
          BI6 = AI4 - AR6
          BI7 = AI5 - AR7
          CR0(K) = BR0 + BR1
          CI0(K) = BI0 + BI1
          IF (J.LE.1) GO TO 10
          CR1(K) = C4*(BR0-BR1) - S4*(BI0-BI1)
          CI1(K) = C4*(BI0-BI1) + S4*(BR0-BR1)
          CR2(K) = C2*(BR2-BI3) - S2*(BI2+BR3)
          CI2(K) = C2*(BI2+BR3) + S2*(BR2-BI3)
          CR3(K) = C6*(BR2+BI3) - S6*(BI2-BR3)
          CI3(K) = C6*(BI2-BR3) + S6*(BR2+BI3)
          TR = P7*(BR5-BI5)
          TI = P7*(BR5+BI5)
          CR4(K) = C1*(BR4+TR) - S1*(BI4+TI)
          CI4(K) = C1*(BI4+TI) + S1*(BR4+TR)
          CR5(K) = C5*(BR4-TR) - S5*(BI4-TI)
          CI5(K) = C5*(BI4-TI) + S5*(BR4-TR)
          TR = -P7*(BR7+BI7)
          TI = P7*(BR7-BI7)
          CR6(K) = C3*(BR6+TR) - S3*(BI6+TI)
          CI6(K) = C3*(BI6+TI) + S3*(BR6+TR)
          CR7(K) = C7*(BR6-TR) - S7*(BI6-TI)
          CI7(K) = C7*(BI6-TI) + S7*(BR6-TR)
          GO TO 20
  10      CR1(K) = BR0 - BR1
          CI1(K) = BI0 - BI1
          CR2(K) = BR2 - BI3
          CI2(K) = BI2 + BR3
          CR3(K) = BR2 + BI3
          CI3(K) = BI2 - BR3
          TR = P7*(BR5-BI5)
          TI = P7*(BR5+BI5)
          CR4(K) = BR4 + TR
          CI4(K) = BI4 + TI
          CR5(K) = BR4 - TR
          CI5(K) = BI4 - TI
          TR = -P7*(BR7+BI7)
          TI = P7*(BR7-BI7)
          CR6(K) = BR6 + TR
          CI6(K) = BI6 + TI
          CR7(K) = BR6 - TR
          CI7(K) = BI6 - TI
  20    CONTINUE
  30  CONTINUE
      RETURN
      END SUBROUTINE R8TX








  
end module lanczos_base


