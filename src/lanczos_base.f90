

module lanczos_base
  implicit none
  private

  public ::LB_alpha , LB_beta, LB_eval, LB_shift, LB_nsteps, LB_allocate_for_lanczos, &
       LB_de_allocate_for_lanczos, LB_cerca, LB_passeggia
  


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
  

contains

  subroutine LB_allocate_for_lanczos( )
    allocate(LB_alpha(0: LB_nsteps ))
    allocate(LB_beta(0: LB_nsteps-1 ))
    allocate(omega( 0:LB_nsteps,  0:LB_nsteps      ) )

    allocate(evect( 0:LB_nsteps-1,  0:LB_nsteps-1    )  )
    allocate(LB_eval(0: LB_nsteps-1 )  )
    allocate(diagwork( 0:LB_nsteps*(3+LB_nsteps)   ) )
    allocate(oldalpha (0: LB_nsteps-1 )         )


    omega(:,:)=0.0D0

  end subroutine LB_allocate_for_lanczos


  subroutine LB_de_allocate_for_lanczos( )
    deallocate(LB_alpha)
    deallocate(LB_beta)
    deallocate(omega )

    deallocate(evect  )
    deallocate(LB_eval  )
    deallocate(diagwork    )
    deallocate(oldalpha          )
  end subroutine LB_de_allocate_for_lanczos
   
  integer function converged(m)
    implicit none
    integer, intent(in):: m
! ::::::::::::::::::::::::::::::::::::::::::
    integer j,i,k
    real(8) dum, dumvect(LB_nsteps )
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
          if (abs(LB_eval(converged)-oldalpha(converged))/abs(oldalpha(converged))  .gt. LANCZOS_tol) then
             exit
          endif
          converged=converged+1
       enddo
    else
       hasoldalpha=1
    endif
    oldalpha(:m-1)=LB_eval(:m-1)
    
    return
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




    lwork = LB_nsteps*(3+LB_nsteps)

    call DSYEV( 'V', 'U', m, evect, m, LB_eval, diagwork, LWORK, INFO )

    
    print *, " evals  " 
    print *, LB_eval-LB_shift


    if(info .ne.0) then
       print *, " problema con dsyev"
       stop
    endif
    return 
  end subroutine diago

  integer function  LB_cerca( nd, shift, tol, set_EP_shift, EP_allocate_for_eigenprob,&
       EP_make_dummy_vectors, get_EP_dim, EP_initialize_start , EP_normalizza,&
       EP_Moltiplica, EP_GramSchmidt ,EP_set_all_random, EP_copy , EP_mat_mult, &
       EP_scalare ,EP_add_from_vect_with_fact)

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
         real(8)  EV(0:m-1,0:k-1 )
       end subroutine        

    end interface
! :::::::::::::::::::::::::::::::::::::::::::::::::::
    integer k,nc,m
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



    do while(nc.lt.nd)

       call diago(k,m)
       nc=converged(m)
       if ( k.gt.0) then
          if( notzero( LB_beta ,k, LANCZOS_tol).eq.0 ) then
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

  integer function notzero(LB_beta, k, tol)
    implicit none
    real(8), intent(in) ::LB_beta(0: LB_nsteps-1 )
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
         real(8)  EV(0:m-1,0:k-1 )
       end subroutine EP_mat_mult
    end interface
    
! :::::::::::::::::::::::::::::::::::::
    integer i
    real(8), pointer :: dumomega(:,:), dumomega2(:,:)
    real(8) acoeff, bcoeff
    integer lda, ldb, ldc



    LB_alpha(0:k-1)=LB_eval(0:k-1)

    LB_beta(0:k-1) = LB_beta(m-1)*evect(m-1, 0:k-1)




    call EP_mat_mult(m,k,  evect(0:m-1, 0:k-1  )   ) ! usare i dumvectors a partire da -1


    

    do i=0, k-1
       call EP_normalizza(-i-1)
       
       call EP_copy(i, -i-1 )
    enddo


    call EP_copy(k, m )
    
  
    allocate(dumomega(0:m-1,0:m-1))
    allocate(dumomega2(0:m-1,0:m-1))
    

    acoeff=1.0D0
    bcoeff =0.0D0
    


    call dgemm('N','N',m,m,m,acoeff,omega ,LB_nsteps+1,  evect ,LB_nsteps,bcoeff,dumomega,m)

    omega(0:k,k)=dumomega(0:k,k)
    omega(k,0:k)=dumomega(0:k,k)

    call dgemm('T','N',m,m,m,acoeff,evect ,LB_nsteps,  dumomega  ,m,bcoeff,dumomega2,m)
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
         real(8)  EV(0:m-1,0:k-1 )
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

       open(unit=22,file="alphabeta_p")
       write(22,*) i, 1.0
       do ipa=0, i-1
          write(22,*) LB_alpha(ipa), LB_beta(ipa)
       enddo
       close(unit=22)


       call EP_Moltiplica(p, i )
       
       
       LB_alpha(i) = EP_scalare(p , i)
       
       ! print *, LB_alpha(i)

       
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

       ! print *, "LB_beta(",i,")=" , LB_beta(i)

       
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
             
             print *, "starting with a vector perpendicular"
             
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
!!$       ! emergenza
       call EP_GramSchmidt(i+1,i+1)
       call EP_normalizza(i+1)


       
!!$
!!$
!!$       print *, "LB_beta(",i,")=" , LB_beta(i)
!!$       ! controlla hermitianicita
!!$       call EP_Moltiplica(p, i+1 )
!!$       print *, "alla rovescio " , EP_scalare(p,i)- LB_beta(i)
!!$       
!!$       if (i.gt.2) then
!!$          print *, "con due prima " , EP_scalare(i-1,i+1)
!!$
!!$       endif


    enddo
    
  end subroutine LB_passeggia
   
end module lanczos_base


