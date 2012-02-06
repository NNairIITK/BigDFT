!> @file 
!!   Diagonalization schme (module)
!! @author
!!   Copyright (C) 2010-2011 BigDFT group 
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 


!> Module to handle diagonalization scheme
module lanczos_base
   use module_base
   implicit none
   private

   public ::LB_alpha ,LB_alpha_cheb , LB_beta, LB_eval, LB_shift, LB_nsteps, LB_allocate_for_lanczos, &
      &   LB_de_allocate_for_lanczos, LB_cerca, LB_passeggia, LB_allocate_for_chebychev,&
      &   LB_iproc, LB_nproc, LB_passeggia_Chebychev, LB_cg, LB_norbp , LB_de_allocate_for_cheb

   character(len=*), parameter :: subname='lanczos_base'

   real(kind=8), pointer :: LB_alpha_cheb(:,:)
   real(kind=8), pointer :: LB_alpha(:)
   real(kind=8), pointer :: LB_beta(:)
   real(kind=8), pointer :: omega(:,:)
   real(kind=8), pointer :: evect(:,:)
   real(kind=8), pointer :: LB_eval(:)
   real(kind=8), pointer :: oldalpha(:)
   real(kind=8), pointer :: diagwork(:)
   real(kind=8) :: LB_shift
   real(kind=8) :: LANCZOS_tol
   integer :: LB_nsteps
   integer :: hasoldalpha
   integer :: i_stat,i_all
   integer :: LB_iproc, LB_nproc, LB_norbp 

   contains

   subroutine LB_allocate_for_chebychev( )
      !not associated at initialisation
      !if(associated(LB_alpha_cheb )) then
      !   call  LB_de_allocate_for_cheb ( )
      !endif

      if(LB_norbp.gt.0) then
         allocate(LB_alpha_cheb( LB_norbp  ,   0: 3*LB_nsteps+ndebug ) , stat=i_stat)
      else
         allocate(LB_alpha_cheb( LB_norbp+1  ,   0: 3*LB_nsteps+ndebug ) , stat=i_stat)      
      endif

      call memocc(i_stat,LB_alpha_cheb,'LB_alpha_cheb',subname)

      !    allocate(LB_alpha(0: 3*LB_nsteps+ndebug ) , stat=i_stat)
      !    call memocc(i_stat,LB_alpha,'LB_alpha',subname)

      !    allocate(LB_beta(0: 1+ndebug ) , stat=i_stat)
      !    call memocc(i_stat,LB_beta,'LB_beta',subname)

      !    allocate(omega( 0:1,  0:1    +ndebug  )  , stat=i_stat)
      !    call memocc(i_stat,omega,'omega',subname)


      !    allocate(evect( 0:1,  0:1   +ndebug )   , stat=i_stat)
      !    call memocc(i_stat,evect,'evect',subname)

      !    allocate(LB_eval(0: 1+ndebug )   , stat=i_stat)
      !    call memocc(i_stat,LB_eval,'LB_eval',subname)

      !    allocate(diagwork( 0:1*(3+1)+ndebug   )  , stat=i_stat)
      !    call memocc(i_stat,diagwork,'diagwork',subname)

      !    allocate(oldalpha (0: 1 +ndebug )        , stat=i_stat  )
      !    call memocc(i_stat,oldalpha,'oldalpha',subname)

   END SUBROUTINE LB_allocate_for_chebychev


   subroutine LB_allocate_for_lanczos( )
      !if(associated(LB_alpha )) then
      !   call  LB_de_allocate_for_lanczos( )
      !endif

      allocate(LB_alpha(0: LB_nsteps+ndebug ) , stat=i_stat)
      call memocc(i_stat,LB_alpha,'LB_alpha',subname)

      allocate(LB_beta(0: LB_nsteps-1+ndebug ) , stat=i_stat)
      call memocc(i_stat,LB_beta,'LB_beta',subname)

      allocate(omega( 0:LB_nsteps,  0:LB_nsteps     +ndebug )  , stat=i_stat)
      call memocc(i_stat,omega,'omega',subname)


      allocate(evect( 0:LB_nsteps-1,  0:LB_nsteps-1   +ndebug )   , stat=i_stat)
      call memocc(i_stat,evect,'evect',subname)

      allocate(LB_eval(0: LB_nsteps-1+ndebug )   , stat=i_stat)
      call memocc(i_stat,LB_eval,'LB_eval',subname)

      allocate(diagwork( 0:LB_nsteps*(3+LB_nsteps+ndebug)   )  , stat=i_stat)
      call memocc(i_stat,diagwork,'diagwork',subname)

      allocate(oldalpha (0: LB_nsteps+ndebug )        , stat=i_stat  )
      call memocc(i_stat,oldalpha,'oldalpha',subname)



      omega(:,:)=0.0D0

   END SUBROUTINE LB_allocate_for_lanczos

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

   END SUBROUTINE LB_de_allocate_for_lanczos

   subroutine LB_de_allocate_for_cheb( )

      i_all=-product(shape(LB_alpha_cheb))*kind(LB_alpha_cheb)
      deallocate(LB_alpha_cheb)
      call memocc(i_stat,i_all,'LB_alpha_cheb',subname)

   END SUBROUTINE LB_de_allocate_for_cheb

   function converged(m)
      implicit none
      integer converged
      integer, intent(in):: m

      integer j,i,k
      real(kind=8) dum, dumvect(0:LB_nsteps )
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

   END FUNCTION converged


   subroutine diago(k,m)
      implicit none
      integer, intent(in) :: k,m
      integer :: i
      integer :: info, lwork

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
   END SUBROUTINE diago

   function  LB_cerca( nd, shift, tol, set_EP_shift, EP_allocate_for_eigenprob,&
         &   EP_make_dummy_vectors, get_EP_dim, EP_initialize_start , EP_normalizza,&
         &   EP_Moltiplica, EP_GramSchmidt ,EP_set_all_random, EP_copy , EP_mat_mult, &
         &   EP_scalare ,EP_add_from_vect_with_fact, accontentati_di_n)

      implicit none

      integer :: LB_cerca
      integer, intent(IN) :: nd
      real(kind=8), intent(IN) :: shift
      real(kind=8), intent(IN) :: tol
      interface
         subroutine set_EP_shift(shift)
            real(kind=8) shift
         END SUBROUTINE set_EP_shift

         subroutine  EP_allocate_for_eigenprob(nsteps)
            implicit none
            integer, intent(in):: nsteps
         END SUBROUTINE
         subroutine EP_make_dummy_vectors(nd)
            implicit none
            integer, intent(in):: nd
            ! :::::::::::::::::::::::    
         END SUBROUTINE
         function get_EP_dim()
            implicit none
            integer :: get_EP_dim
         END FUNCTION

         subroutine EP_initialize_start()
         END SUBROUTINE 
         subroutine EP_normalizza(i)
            integer, intent(in) :: i
         END SUBROUTINE 
         subroutine EP_Moltiplica(i,j)
            integer, intent(in) :: i,j
         END SUBROUTINE 
         function EP_scalare(i,j)
            real(kind=8) :: EP_scalare
            integer, intent(in) :: i,j
         END FUNCTION 
         subroutine EP_add_from_vect_with_fact( i, j  ,   a )
            integer, intent(in) :: i,j
            real(kind=8), intent(in) :: a
         END SUBROUTINE 
         subroutine EP_GramSchmidt(i,j)
            integer, intent(in) :: i,j
         END SUBROUTINE 
         subroutine EP_set_all_random(i)
            integer, intent(in) :: i
         END SUBROUTINE 
         subroutine EP_copy(i,j)
            integer, intent(in) :: i,j
         END SUBROUTINE
         subroutine EP_mat_mult(m,k ,  EV )
            integer, intent(in):: m,k
            real(kind=8), intent(in):: EV(1 )
         END SUBROUTINE        

      end interface

      integer , optional :: accontentati_di_n

      integer :: k,nc,m
      integer :: n_converged_cercati

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
         &   EP_Moltiplica, EP_GramSchmidt ,EP_set_all_random, EP_copy,   EP_mat_mult, &
         &   EP_scalare,EP_add_from_vect_with_fact     )

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
            &   EP_Moltiplica, EP_GramSchmidt ,EP_set_all_random, EP_copy ,&
            &   EP_mat_mult,   EP_scalare,EP_add_from_vect_with_fact)
      enddo
      if( m.eq.get_EP_dim()) then
         LB_cerca=m
      else
         LB_cerca= k
      endif
      return 
   END FUNCTION LB_cerca

   function notzero(k, tol)
      implicit none
      integer :: notzero
      integer, intent(in) :: k
      real(kind=8), intent(in) :: tol
      !::::::::::::::::::::::::::::::`
      integer i
      notzero=0
      do i=0,k-1
         if( abs(LB_beta(i)).gt.tol) then
            notzero=notzero+1
         endif
      enddo
      return 
   END FUNCTION notzero


   subroutine ricipolla(k,m, EP_normalizza ,EP_copy, EP_mat_mult )
      implicit none
      integer, intent(in):: k,m
      interface
         subroutine EP_normalizza(i)
            integer, intent(in) :: i
         END SUBROUTINE EP_normalizza
         subroutine EP_copy(i,j)
            integer, intent(in) :: i,j
         END SUBROUTINE EP_copy
         subroutine EP_mat_mult(m,k ,  EV )
            integer, intent(in)::  m,k
            real(kind=8), intent(in)::  EV(1 )
         END SUBROUTINE EP_mat_mult
      end interface

      integer :: i
      real(kind=8), pointer :: dumomega(:,:), dumomega2(:,:)
      real(kind=8) :: acoeff, bcoeff

      LB_alpha(0:k-1)=LB_eval(0:k-1)

      LB_beta(0:k-1) = LB_beta(m-1)*evect(m-1, 0:k-1)

      call EP_mat_mult(m,k,  evect(0:, 0:  )   ) 

      do i=0, k-1
         call EP_normalizza(-i-1)

         call EP_copy(i, -i-1 )
      enddo

      call EP_copy(k, m )

      allocate(dumomega(0:m-1,0:m-1+ndebug))
      allocate(dumomega2(0:m-1,0:m-1+ndebug))

      acoeff=1.0D0
      bcoeff =0.0D0

      call dgemm('N','N',m,m,m,acoeff,omega(0,0) ,LB_nsteps+1,  evect(0,0) ,LB_nsteps,bcoeff,dumomega(0,0),m)

      omega(0:k,k)=dumomega(0:k,k)
      omega(k,0:k)=dumomega(0:k,k)

      call dgemm('T','N',m,m,m,acoeff,evect(0,0) ,LB_nsteps,  dumomega (0,0) ,m,bcoeff,dumomega2(0,0),m)
      omega(0:k-1,0:k-1)= dumomega2 (0:k-1,0:k-1)


      deallocate(dumomega)
      deallocate( dumomega2)


   END SUBROUTINE ricipolla



   subroutine LB_passeggia( k, m, get_EP_dim, EP_initialize_start , EP_normalizza, &
         &   EP_Moltiplica, EP_GramSchmidt ,EP_set_all_random, EP_copy,  EP_mat_mult,&
         &   EP_scalare,EP_add_from_vect_with_fact )

      implicit none
      integer, intent(in) :: k
      integer, intent(in) :: m
      real(kind=8) :: sn,eu,eusn,eq
      real(kind=8) :: max0
      integer :: i,l,j, w

      interface
         subroutine EP_mat_mult(m,k ,  EV )
            integer, intent(in)::  m,k
            real(kind=8), intent(in)::  EV(1 )
         END SUBROUTINE  
         function get_EP_dim()
            implicit none
            integer :: get_EP_dim
         END FUNCTION
         subroutine EP_initialize_start()
         END SUBROUTINE 
         subroutine EP_normalizza(i)
            integer, intent(in) :: i
         END SUBROUTINE 
         subroutine EP_Moltiplica(i,j)
            integer, intent(in) :: i,j
         END SUBROUTINE 
         function EP_scalare(i,j)
            implicit none
            real(kind=8) :: EP_scalare
            integer, intent(in) :: i,j
         END FUNCTION 
         subroutine EP_add_from_vect_with_fact( i, j  ,   a )
            integer, intent(in) :: i,j
            real(kind=8), intent(in) :: a
         END SUBROUTINE 
         subroutine EP_GramSchmidt(i,j)
            integer, intent(in) :: i,j
         END SUBROUTINE 
         subroutine EP_set_all_random(i)
            integer, intent(in) :: i
         END SUBROUTINE 
         subroutine EP_copy(i,j)
            integer, intent(in) :: i,j
         END SUBROUTINE
      end interface

      integer :: p
      real(kind=8) :: add
      real(kind=8) :: condition
      integer :: ipa

      p = -1

      sn   = sqrt(get_EP_dim()*1.0D0 )
      eu   = 1.1D-15
      eusn = eu*sn
      eq   = sqrt(eu)

      if (k == 0) then
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

   END SUBROUTINE LB_passeggia


   subroutine CalcolaSpettroChebychev( cheb_shift, fact_cheb,   Nu, xabs_res_prefix, norbp , wgts) 
      real(gp) cheb_shift, fact_cheb
      integer Nu, norbp
      character(len=*) :: xabs_res_prefix
      real(gp) wgts(norbp)


      real(gp), pointer :: Xs(:), res(:), cfftreal(:),result(:), cfftimag(:), zinout(:,:,:)
      complex(kind=8), pointer :: alphas(:), expn(:)
      complex(kind=8) :: ctemp
      integer :: inzee

      integer :: Nbar
      integer :: i, iorb, ierr
      real(gp) :: Pi
      character(len=1000) :: filename
      real(gp) :: fact


      Pi=acos(-1.0_gp)
      Nbar =1
      do while(Nbar<Nu) 
         Nbar=Nbar*2
      enddo
      Nbar=Nbar*2   

      allocate(Xs(0:Nbar-1+ndebug) , stat=i_stat)
      call memocc(i_stat,Xs,'Xs',subname)

      allocate(res(0:Nbar-1+ndebug) , stat=i_stat)
      call memocc(i_stat,res,'res',subname)


      !! memocc does not work in complex
      allocate(alphas(0:Nbar-1+ndebug) , stat=i_stat)
      !! call memocc(i_stat,alphas,'alphas',subname)

      !! memocc does not work in complex
      allocate(expn(0:Nbar-1+ndebug) , stat=i_stat)


      allocate(cfftreal(0:2*Nbar-1+ndebug) , stat=i_stat)
      call memocc(i_stat,cfftreal,'cfftreal',subname)

      allocate(result(0:2*Nbar-1+ndebug) , stat=i_stat)
      call memocc(i_stat,result,'result',subname)


      allocate(cfftimag(0:2*Nbar-1+ndebug) , stat=i_stat)
      call memocc(i_stat,cfftimag,'cfftimag',subname)

      allocate(zinout (2,2*Nbar,2+ndebug) , stat=i_stat)
      call memocc(i_stat,zinout,'zinout',subname)

      result=0.0_gp

      do iorb=1, norbp

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
                  &   sin(pi* i /(Nu+1.0))* cos(pi/(Nu+1))/sin(pi/(Nu+1))    )
               ctemp= LB_alpha_cheb(iorb, i)*fact*exp(   ((0.0_gp,1.0_gp) * ( Nbar -0.5_gp)) *i*Pi/Nbar     )
               cfftreal(i)=REAL(ctemp) *wgts(iorb)
               cfftimag(i)=AIMAG(ctemp) *wgts(iorb)
            else
               cfftreal(i)=0
               cfftimag(i)=0
            endif
         enddo

         zinout=0.0

         zinout(1,1:2*Nbar,1) = cfftreal(:)
         zinout(2,1:2*Nbar,1) = cfftimag(:)


         call dcopy(2*Nbar,cfftreal(0),1,zinout(1,1,1),2)
         call dcopy(2*Nbar,cfftimag(0),1,zinout(2,1,1),2)
         !zinout(1,1:2*Nbar,1) = cfftreal(:)
         !zinout(2,1:2*Nbar,1) = cfftimag(:)


         call fft_1d_ctoc(1 ,1, 2*Nbar ,zinout(1,1,1) ,inzee)

         call dcopy(2*Nbar,zinout(1,1,inzee),2,cfftreal(0),1)
         call dcopy(2*Nbar,zinout(2,1,inzee),2,cfftimag(0),1)
         !cfftreal(:)  =   zinout(1,1:2*Nbar,inzee)   
         !cfftimag(:)  =   zinout(2,1:2*Nbar,inzee)    

         cfftreal=cfftreal-(LB_alpha_cheb(iorb, 0)*0.5) *wgts(iorb)

         !!$    do n=1,Nu-1
         !!$       expn(:)=expn(:)*alphas(:)
         !!$       fact = 1.0/(Nu+1)*( (Nu-n +1.0)*cos(pi* n /(Nu+1))+ &
               !!$       sin(pi* n /(Nu+1.0))* cos(pi/(Nu+1))/sin(pi/(Nu+1))    )
         !!$       res(:)=res(:)+2*REAL(expn(:))*LB_alpha(n)*fact
         !!$    enddo

         print *, " done " 

         !!$    res =res/Pi/sqrt(1-Xs*Xs)

         cfftreal(0:Nbar-1) =2*cfftreal(0:Nbar-1)/Pi/sqrt(1-Xs*Xs)*fact_cheb


         result=result+cfftreal


         write(filename,'(a,a,I0,a,I0)') trim(xabs_res_prefix),  "cheb_spectra_proc_", LB_iproc,"_orbp_", iorb
         write(*, '(a,1x, a)') " writing spectra to " , trim(filename) 
         open(unit=22,file=filename)
         do i=0, Nbar-1
            write(22,*) Xs(i) / fact_cheb + cheb_shift , cfftreal(i)/wgts(iorb)
         enddo
         close(unit=22)

      end do

      cfftreal=0

      if(LB_nproc/=1) then
         call MPI_REDUCE(result(0),cfftreal(0), 2*Nbar ,mpidtypd,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      else
         cfftreal=result
      endif


      if(LB_iproc==0) then
         write(filename,'(a,a,I0)') trim(xabs_res_prefix),  "cheb_spectra_", Nu+1
         write(*, '(a,1x, a)') " writing spectra to " , trim(filename) 

         open(unit=22,file=filename)
         do i=0, Nbar-1
            write(22,*) Xs(i) / fact_cheb + cheb_shift , cfftreal(i)
         enddo
         close(unit=22)
      endif

      i_all=-product(shape(Xs))*kind(Xs)
      deallocate(Xs)
      call memocc(i_stat,i_all,'Xs',subname)

      i_all=-product(shape(res))*kind(res)
      deallocate(res)
      call memocc(i_stat,i_all,'res',subname)

      !! memocc does not work in complex
      i_all=-product(shape(alphas))*kind(alphas)
      deallocate(alphas)
      !! call memocc(i_stat,i_all,'alphas',subname)


      !! memocc does not work in complex
      i_all=-product(shape(expn))*kind(expn)
      deallocate(expn)
      !! call memocc(i_stat,i_all,'expn',subname)

      i_all=-product(shape(cfftreal))*kind(cfftreal)
      deallocate(cfftreal)
      call memocc(i_stat,i_all,'cfftreal',subname)

      i_all=-product(shape(result))*kind(result)
      deallocate(result)
      call memocc(i_stat,i_all,'result',subname)

      i_all=-product(shape(cfftimag))*kind(cfftimag)
      deallocate(cfftimag)
      call memocc(i_stat,i_all,'cfftimag',subname)


      i_all=-product(shape(zinout))*kind(zinout)
      deallocate(zinout)
      call memocc(i_stat,i_all,'zinout',subname)

   END SUBROUTINE CalcolaSpettroChebychev


   subroutine LB_passeggia_Chebychev (  m, cheb_shift,  fact_cheb,  get_EP_dim,  EP_initialize_start , EP_normalizza, &
         &   EP_Moltiplica, EP_GramSchmidt ,EP_set_all_random, EP_copy,  EP_mat_mult,&
         &   EP_scalare_multik,EP_add_from_vect_with_fact , EP_multbyfact, EP_ApplySinv, EP_ApplyS, dopaw, &
         &   S_do_cg ,Sinv_do_cg , xabs_res_prefix, nkpts, norb_par, wgts)
      use module_base
      implicit none
      integer, intent(in) :: m
      real(gp) :: cheb_shift, fact_cheb
      logical dopaw
      character(len=*) :: xabs_res_prefix
      integer, intent(in) :: nkpts, norb_par(0:LB_nproc-1)
      real(gp) wgts(LB_nproc )
      integer ierr

      interface
         subroutine EP_mat_mult(m,k ,  EV )
            integer, intent(in)::  m,k
            real(kind=8), intent(in)::  EV(1 )
         END SUBROUTINE  
         function get_EP_dim()
            implicit none
            integer :: get_EP_dim
         END FUNCTION
         subroutine EP_initialize_start()
         END SUBROUTINE 
         subroutine EP_normalizza(i)
            integer, intent(in) :: i
         END SUBROUTINE 
         subroutine EP_Moltiplica(i,j)
            integer, intent(in) :: i,j
         END SUBROUTINE 
         function EP_scalare(i,j)
            implicit none
            real(kind=8) :: EP_scalare
            integer, intent(in) :: i,j
         END FUNCTION

         subroutine EP_scalare_multik(i,j, scalari)
            use module_base
            integer, intent(in) :: i,j
            real(wp) :: scalari(*)
         END SUBROUTINE EP_scalare_multik

         subroutine EP_add_from_vect_with_fact( i, j  ,   a )
            integer, intent(in) :: i,j
            real(kind=8), intent(in) :: a
         END SUBROUTINE 

         subroutine EP_GramSchmidt(i,j)
            integer, intent(in) :: i,j
         END SUBROUTINE 

         subroutine EP_set_all_random(i)
            integer, intent(in) :: i
         END SUBROUTINE 

         subroutine EP_copy(i,j)
            integer, intent(in) :: i,j
         END SUBROUTINE

         subroutine EP_multbyfact(j, fact)
            use module_base
            implicit none
            integer, intent(in) :: j
            real(gp) :: fact
         END SUBROUTINE EP_multbyfact

         subroutine EP_ApplySinv(p,i)
            implicit none
            integer, intent(in) :: p,i
         END SUBROUTINE EP_ApplySinv

         subroutine EP_ApplyS(p,i)
            implicit none
            integer, intent(in) :: p,i
         END SUBROUTINE EP_ApplyS

      end interface

      integer :: i
      integer :: tmp1, attuale, precedente
      integer :: Stmp1, Sattuale, Sprecedente, dspl( 0:LB_nproc-1)
      logical :: S_do_cg ,Sinv_do_cg 
      real(gp) tol
      real(wp) alphacollect(nkpts)

      if(LB_iproc == 0) print *, "LB_nproc= ", LB_nproc

      if(LB_nproc/=1) then
         dspl(0)=0
         do i=1, LB_nproc-1
            dspl(i ) = dspl(i-1)+norb_par(i-1)
         end do
      endif

      attuale= 1
      Sattuale= 2
      tmp1 = 3
      precedente= 4

      Stmp1 = 5
      Sprecedente= 6

      call EP_initialize_start()
      if (LB_iproc == 0) write(*,*) " EP_initialize_start    OK "
      call EP_copy(attuale,0)
      call EP_copy(precedente,0)
      if(dopaw) then
         tol=1.0D-20

         if(   S_do_cg ) then
            !!$          call  LB_generic_cg(    get_EP_dim,EP_normalizza,&
               !!$               EP_ApplySinv ,  EP_copy,  &
               !!$               EP_scalare,EP_add_from_vect_with_fact    , EP_multbyfact,  tol, Sattuale, attuale)
            STOP " S_do_cg temporarily disactivate beacause chebychev is potentially multik " 


            !!$       call  EP_ApplySinv( tmp1 , Sattuale)
            !!$       call EP_add_from_vect_with_fact( tmp1  , attuale ,-1.0_gp ) 
            !!$       print *, " residuo " , EP_scalare(tmp1,tmp1) 
            !!$       STOP "  !! qua fare un test rimoltiplicando per sinv "

         else
            call  EP_ApplyS(Sattuale   ,attuale )
         endif

         call EP_scalare_multik(attuale,Sattuale, LB_alpha_cheb(:,0)  ) 
      else
         call EP_scalare_multik(attuale,attuale, LB_alpha_cheb(:,0) )
      endif

      do i=0, m-1
         if( modulo( m-1-i   ,100)==0 ) then
            !!$          open(unit=22,file=(trim(xabs_res_prefix)//"coeffs_chebychev"))
            !!$          write(22,*) 2*i+1,  cheb_shift, fact_cheb
            !!$          do ipa=0, 2*i
            !!$             write(22,*) LB_alpha(ipa)
            !!$          enddo
            !!$          close(unit=22)
            call CalcolaSpettroChebychev( cheb_shift,&
               &   fact_cheb,   2*i+1 , xabs_res_prefix, &
               &   norb_par(LB_iproc), wgts) 
         endif

         call EP_Moltiplica(tmp1, attuale )
         if(dopaw) then
            call EP_add_from_vect_with_fact( tmp1 , Sattuale, -cheb_shift ) 
         else
            call EP_add_from_vect_with_fact(  tmp1 ,attuale  ,-cheb_shift )           
         endif


         call EP_multbyfact(tmp1,fact_cheb)
         if(i==0) then
            !!$          LB_alpha(1)=EP_scalare(tmp1,attuale)  
            call EP_scalare_multik(tmp1,attuale, LB_alpha_cheb(:, 1))  
            if(dopaw) then
               call EP_copy(Stmp1,tmp1)
               if(   Sinv_do_cg ) then
                  !!$                call  LB_generic_cg(    get_EP_dim,EP_normalizza,&
!!$                     EP_ApplyS ,  EP_copy,  &
!!$                     EP_scalare,EP_add_from_vect_with_fact    , EP_multbyfact,  tol,tmp1 ,Stmp1 )  
                  STOP " Sinv_do_cg temporarily disactivate beacaause chebychev is potentially multik" 

               else
                  call EP_ApplySinv(tmp1, Stmp1 )             
               endif
            end if
         else
            if(dopaw) then
               call EP_copy(Stmp1,tmp1)
               call EP_multbyfact(Stmp1,2.0_gp)
               call EP_add_from_vect_with_fact( Stmp1 , Sprecedente,-1.0_gp ) 

               if(   Sinv_do_cg ) then
                  !!$                call  LB_generic_cg(    get_EP_dim,EP_normalizza,&
                     !!$                     EP_ApplyS ,  EP_copy,  &
                     !!$                     EP_scalare,EP_add_from_vect_with_fact    , EP_multbyfact,  tol, tmp1,Stmp1 )
                  STOP " Sinv_do_cg temporarily disactivate beacaause chebychev is potentially multik " 
               else
                  call EP_ApplySinv(tmp1, Stmp1 )             
               endif
            else
               call EP_multbyfact(tmp1,2.0_gp)
               call EP_add_from_vect_with_fact( tmp1 , precedente,-1.0_gp )              
            end if
         endif
         if(dopaw) then
            if(i.gt.0) then
               call  EP_scalare_multik(tmp1,Sattuale,  LB_alpha_cheb(:,2*i+1)  )
               LB_alpha_cheb(:,2*i+1)= 2.0_wp * LB_alpha_cheb(:,2*i+1)  - LB_alpha_cheb(:, 1)
               !!$ LB_alpha(2*i+1) = 2*EP_scalare(tmp1,Sattuale)-LB_alpha(1)
            endif
            call EP_scalare_multik(tmp1,Stmp1, LB_alpha_cheb(:, 2*i+2) )
            LB_alpha_cheb(:, 2*i+2) = 2_wp*LB_alpha_cheb(:, 2*i+2)  -LB_alpha_cheb(:,0)
            !!$ LB_alpha(2*i+2) = 2*EP_scalare(tmp1,Stmp1)-LB_alpha(0)
         else

            call  EP_scalare_multik(tmp1, attuale,  LB_alpha_cheb(:,2*i+1)  )
            LB_alpha_cheb(:,2*i+1)= 2.0_wp * LB_alpha_cheb(:,2*i+1)  - LB_alpha_cheb(:, 1)

            !!$ LB_alpha(2*i+1) = 2*EP_scalare(tmp1,attuale)-LB_alpha(1)


            call EP_scalare_multik(tmp1,tmp1, LB_alpha_cheb(:, 2*i+2))
            LB_alpha_cheb(:, 2*i+2) = 2_wp*LB_alpha_cheb(:, 2*i+2)  -LB_alpha_cheb(:,0)

            !!$ LB_alpha(2*i+2) = 2*EP_scalare(tmp1,tmp1)-LB_alpha(0)
         end if

         if(LB_nproc/=1) then
            call MPI_Gatherv( LB_alpha_cheb(1,2*i+1) , norb_par(LB_iproc),mpidtypw, &
               &   alphacollect(1),norb_par(0), dspl(0), mpidtypw, 0, MPI_COMM_WORLD, ierr)
            call MPI_Gatherv( LB_alpha_cheb(1,2*i+2) , norb_par(LB_iproc),mpidtypw, &
               &   alphacollect(1),norb_par(0), dspl(0), mpidtypw, 0, MPI_COMM_WORLD, ierr)
         end if
         if(LB_iproc == 0) then
            write (*,"(A)",advance="no")   "new cheb coeffs "
            write( *,"(1x,10(d20.10,1x))")   LB_alpha_cheb(1:LB_norbp ,2*i+1) , LB_alpha_cheb(1:LB_norbp,2*i+2)
         endif
         call EP_copy(precedente,attuale)
         call EP_copy(attuale,tmp1)

         if(dopaw) then
            call EP_copy(Sprecedente,Sattuale)
            call EP_copy(Sattuale,Stmp1)
         end if
      enddo

   END SUBROUTINE LB_passeggia_Chebychev

   function  LB_cg(    get_EP_dim, EP_initialize_start , EP_normalizza,&
         &   EP_Moltiplica4spectra,  EP_copy,  &
         &   EP_scalare,EP_add_from_vect_with_fact    , EP_multbyfact ,EP_precondition , Ene, gamma, tol, useold )


      use module_base
      implicit none
      real(gp) :: LB_cg
      real(gp) :: ene, gamma, tol
      logical :: useold

      interface
         integer function get_EP_dim()
      END FUNCTION
      subroutine EP_initialize_start()
      END SUBROUTINE 
      subroutine EP_normalizza(i)
         integer , intent(in)  :: i
      END SUBROUTINE 
      subroutine EP_Moltiplica4spectra(i,j, ene, gamma)
         use module_base
         integer , intent(in)  :: i,j
         real(gp) :: ene, gamma
      END SUBROUTINE 
      function EP_scalare(i,j)
         implicit none
         real(kind=8) :: EP_scalare
         integer , intent(in)  :: i,j
      END FUNCTION 
      subroutine EP_add_from_vect_with_fact( i, j  ,   a )
         integer , intent(in)  :: i,j
         real(kind=8) , intent(in)  :: a
      END SUBROUTINE 
      subroutine EP_copy(i,j)
         integer , intent(in)  :: i,j
      END SUBROUTINE
      subroutine EP_multbyfact(j, fact)
         use module_base
         implicit none
         integer, intent(in) :: j
         real(gp) :: fact
         ! ::::::::::::::::::::::::::::::::::::::
      END SUBROUTINE EP_multbyfact
      subroutine EP_precondition(p,i, ene, gamma)
         use module_interfaces
         !Arguments
         use module_base
         implicit none
         integer, intent(in) :: p,i
         real(gp) ene, gamma
      END SUBROUTINE EP_precondition

   end interface

   real(gp)::rho, beta, err, alpha, err0, rhoold, maxerr
   integer :: k
   integer :: b, x, r, p, z, Ap, remember


   b=0
   x=1
   r=2
   p=3
   z=4
   Ap=5
   remember=6


   call EP_initialize_start()

   err0= EP_scalare(b,b)

   if(.not. useold) then
      call EP_copy(r,b)
      call EP_precondition(z,r, ene, gamma)
      call EP_multbyfact(x,0.0_gp)
   else
      call EP_copy(x,remember)
      call EP_Moltiplica4spectra(r, x, ene, gamma)
      call EP_multbyfact(r,-1.0_gp)
      call EP_add_from_vect_with_fact(r,b,1.0_gp)
      call EP_precondition(z,r, ene, gamma)
   endif

   beta = 0.0
   rho = EP_scalare(z,r)
   err= rho
   k = 0

   print *, " initial error " , err0
   call EP_copy(p,z)

   maxerr=err0

   do while(k<500 .and. err/maxerr>tol) 
      if(err>maxerr) maxerr=err
      call EP_Moltiplica4spectra(Ap, p, ene, gamma)

      alpha = rho/ EP_scalare( p,Ap)



      call EP_add_from_vect_with_fact(x,p,alpha)
      call EP_add_from_vect_with_fact(r, Ap,-alpha)

      ! precondiziona qui 
      call EP_precondition(z,r, ene, gamma)
      ! call EP_copy(z,r)

      rhoold=rho
      rho = EP_scalare(z,r)

      beta= rho/rhoold

      call EP_multbyfact(p,beta)

      call EP_add_from_vect_with_fact( p,z,1.0_gp  )

      err = EP_scalare(r,r)
      !! write(*,'(1x,A30,ES13.2,ES13.4)') " CG ERR ,  result :  ", err,   EP_scalare(x,b)
      k = k+1
   enddo
   write(*,'(1x,A30,ES13.2,ES13.4)') " CG ERR ,  result :  ", err,   EP_scalare(x,b)
   LB_cg= EP_scalare(x,b)
   call EP_copy(remember,x)

END FUNCTION LB_cg




subroutine  LB_generic_cg(    get_EP_dim,  EP_normalizza,&
      &   EP_gen_Moltiplica,  EP_copy,  &
      &   EP_scalare,EP_add_from_vect_with_fact    , EP_multbyfact,  tol, x,b )


   use module_base
   implicit none
   real(gp) tol

   interface
      function get_EP_dim()
         implicit none
         integer :: get_EP_dim
      END FUNCTION
      subroutine EP_normalizza(i)
         integer , intent(in):: i
      END SUBROUTINE 
      subroutine EP_gen_Moltiplica(i,j)
         use module_base
         integer, intent(in) :: i,j
         real(gp) :: ene, gamma
      END SUBROUTINE EP_gen_Moltiplica
      function EP_scalare(i,j)
         real(kind=8) :: EP_scalare
         integer , intent(in) :: i,j
      END FUNCTION 
      subroutine EP_add_from_vect_with_fact( i, j  ,   a )
         integer , intent(in) :: i,j
         real(kind=8) , intent(in)  :: a
      END SUBROUTINE 
      subroutine EP_copy(i,j)
         integer , intent(in) :: i,j
      END SUBROUTINE
      subroutine EP_multbyfact(j, fact)
         use module_base
         implicit none
         integer, intent(in) :: j
         real(gp) :: fact
         ! ::::::::::::::::::::::::::::::::::::::
      END SUBROUTINE EP_multbyfact


   end interface

   real(gp)::rho, beta, err, alpha, err0, rhoold, maxerr
   integer :: k
   integer :: b, x, r, p, Ap

   !!  b=0
   !!  x=1
   r=-1
   p=-2
   Ap=-3

   err0= EP_scalare(b,b)

   call EP_copy(r,b)
   call EP_multbyfact(x,0.0_gp)

   beta = 0.0
   rho = EP_scalare(r,r)
   err= rho
   k = 0

   print *, " initial error " , err0
   call EP_copy(p,r)

   maxerr=err0

   do while(k<500 .and. err/maxerr>tol) 
      if(err>maxerr) maxerr=err
      call EP_gen_Moltiplica(Ap, p)

      alpha = rho/ EP_scalare( p,Ap)

      call EP_add_from_vect_with_fact(x,p,alpha)
      call EP_add_from_vect_with_fact(r, Ap,-alpha)

      rhoold=rho
      rho = EP_scalare(r,r)

      beta= rho/rhoold

      call EP_multbyfact(p,beta)

      call EP_add_from_vect_with_fact( p,r,1.0_gp  )

      err = rho
      print * , " err " , err
      k = k+1
   enddo
   write(*,'(1x,A30,ES13.2)') " generic CG ERR ,  result :  ", err


   return  

END SUBROUTINE LB_generic_cg

END MODULE lanczos_base
