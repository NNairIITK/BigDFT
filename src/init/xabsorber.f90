


subroutine find_pfproj( Nsol,Ngrid,rgrid,  psi1s, psigrid, real_start, psigrid_pseudo)
  use module_base
  implicit none
 
  integer, intent(in) :: Nsol,Ngrid, real_start
  real(gp), intent(in) :: psi1s(Ngrid), rgrid(Ngrid)
  real(gp), intent(inout) :: psigrid(Ngrid,Nsol),psigrid_pseudo(Ngrid,Nsol)
  ! ------------------------------------------------------------------------------------
  real(gp) dumgrid(Ngrid), coeffs(Nsol), dumgrid2(Ngrid), mass, mass_pseudo
  integer segno(Nsol), segno_pseudo(Nsol)
  integer i,j,k
  do i=1, Nsol
     do k=1, Ngrid
        dumgrid(k)=psigrid(k,i)*psi1s(k)
     enddo
     call integrate(dumgrid, dumgrid2, rgrid, Ngrid)
     coeffs(i)=dumgrid2(Ngrid)

     mass=0.0D0
     mass_pseudo=0.0D0

     do k=1, Ngrid
        if( abs(psigrid(k,i))>mass) mass= abs(psigrid(k,i))
        if( abs(psigrid_pseudo(k,i))>mass_pseudo) mass_pseudo= abs(psigrid_pseudo(k,i))
     enddo
     

     do k=Ngrid, 1,-1
        if( abs(psigrid(k,i))>mass*0.01) then
           if( psigrid(k,i).gt.0.0D0) then
              segno(i)=1.0
           else
              segno(i)=-1
           endif
           exit
        endif
     enddo


     do k=Ngrid, 1,-1
        if( abs(psigrid_pseudo(k,i))>mass_pseudo*0.01) then
           if( psigrid_pseudo(k,i).gt.0.0D0) then
              segno_pseudo(i)=1
           else
              segno_pseudo(i)=-1
           endif
           exit
        endif
     enddo

  enddo

  call  DGEMM('N','N',Ngrid ,1,   Nsol,1.0d0 ,psigrid , Ngrid ,coeffs ,Nsol, 0.0D0 , dumgrid , Ngrid)


  do i=1,Nsol
     coeffs(i)=coeffs(i)*segno(i)*segno_pseudo(i)
  enddo

  call  DGEMM('N','N',Ngrid ,1,   Nsol-real_start+1  ,1.0d0 ,psigrid_pseudo , Ngrid ,&
       coeffs(real_start) ,Nsol-real_start+1, 0.0D0 , dumgrid2 , Ngrid)

  psigrid(:,1)=dumgrid
  psigrid_pseudo(:,1)=dumgrid2

  return
end subroutine find_pfproj
  


subroutine     find_Scoeffs_grid( ng,  expo, Ngrid, rgrid, psi1s , gcoeffs , l )
  use module_base
  implicit none
 
  integer, intent(in) :: ng,  l, Ngrid
  real(gp), intent(in) :: rgrid(Ngrid), psi1s(Ngrid)
  real(gp), dimension(ng+1), intent(in) :: expo
  real(gp), dimension(0:ng), intent(out) :: gcoeffs

  ! -----------------------------------------------------------
  real(gp)  ::   Soverlap(0:ng,0:ng)
  real(gp)  ::   Score(0:ng), dumgrid(Ngrid), dumgrid2(Ngrid)
  integer :: i,j,k,n,m, iw, INFO, LWORK, nord
  real(gp) :: a1,b1,a2,b2, A,B, spi, pi
  real(gp) :: W(0:ng), WORK(3*(ng+1)*(ng+1))
  real(gp)::  sum, totalpow, gin, gim, gip, ggg, gamma



  lwork= 3*(ng+1)*(ng+1)

  spi=1.772453850905516_gp
  pi=spi*spi

  totalpow=2.0_gp + 2*l  
  do i=0,ng
     b1 = 1.0/2.0/expo(i+1)/expo(i+1)
     do j=0,ng
        Soverlap(i,j)=0.0
        b2 = 1.0/2.0/expo(j+1)/expo(j+1)
        
        B=b1+b2
        
        ggg= gamma( (1.0_gp+totalpow)/2.0_gp   )
        Soverlap(i,j)=0.5_gp* ggg * B**(-(1.0D0+totalpow)/2 )
     enddo
     do k=1, Ngrid
        dumgrid(k)= psi1s(k)*exp( -b1*rgrid(k)*rgrid(k) ) *rgrid(k)**( l+1 ) 
     enddo
     call integrate(dumgrid, dumgrid2, rgrid, Ngrid)
     Score(i) = dumgrid2(Ngrid)
  enddo
  

  call DSYEV( 'V', 'L', ng+1, Soverlap(0,0) , ng+1, W(0), WORK, LWORK, INFO )
  



  gcoeffs(:)=0.0
  do n=0, ng
     if( W(n).gt. 1.0e-13*W(ng)    ) then
        sum=0.0
        do i=0, ng
           sum =  sum+ Soverlap(i,n)*Score(i)
        enddo
        sum=sum/W(n)
        do i=0,ng
           gcoeffs(i)=gcoeffs(i)+Soverlap(i,n)*sum
        enddo
     endif
  enddo
  

  return 
end subroutine find_Scoeffs_grid




subroutine   dump_1gauwf_on_radgrid(prefix, ng , expo,psi   ,lpow   )
  use module_base
  implicit none
 
  character(*) , intent(in) ::  prefix
  integer, intent(in) :: ng,lpow
  real(gp), dimension(ng+1), intent(in) :: expo
  real(gp), dimension(0:ng), intent(in) :: psi

  ! local
  integer, parameter :: n_int=100
  character(200)  filename
  integer l,i, iocc,ig
  real(8) r,sum



  write(filename,'(a)') prefix


  open(unit=22,file=filename)
  do i=1, 2000
     r=0.01*i
     sum=0.0
     do ig = 0,ng
        sum=sum+psi(ig)*exp( -r*r/2.0/expo(ig+1)/expo(ig+1) )
     enddo
     write(22,*) r, sum*(r**lpow)
  enddo
  close(unit=22)

return 
end subroutine dump_1gauwf_on_radgrid

real(gp)  function  value_at_r(r, ng , expo,psi     )
  use module_base, only: gp

  implicit none
 
  real(gp) , intent(in) ::  r
  integer, intent(in) :: ng
  real(gp), dimension(ng+1), intent(in) :: expo
  real(gp), dimension(0:ng), intent(in) :: psi

  ! local
  integer, parameter :: n_int=100

  integer l,i, iocc,ig
  real(gp) sum



  sum=0.0
  do ig = 0,ng
     sum=sum+psi(ig)*exp( -r*r/2.0/expo(ig+1)/expo(ig+1) )
  enddo
  value_at_r=sum

  return 
end function value_at_r


subroutine dump_gauwf_on_radgrid(prefix  ,ng ,noccmax , lmax , expo,psi,aeval, occup     )
  use module_base, only: gp
  implicit none
 
  character(*) , intent(in) ::  prefix
  integer, intent(in) :: ng,noccmax, lmax
  real(gp), dimension(ng+1), intent(in) :: expo
  real(gp), dimension(0:ng,noccmax,lmax+1), intent(in) :: psi
  real(gp), dimension(noccmax,lmax+1  ), intent(in) ::  aeval,occup

  ! local
  integer, parameter :: n_int=100
  character(200)  filename
  integer l,i,k, iocc,ig
  real(8) r,sum

  do i=1,noccmax
     do l=0,lmax
  
        write(filename,'(a,a1,i1,a1,i1)') prefix,'_',i,'_',l
        
        open(unit=22,file=filename)
        do k=1 ,2000
           r=0.01*k
           sum=0.0
           do ig = 0,ng
              sum=sum+psi(ig,i,l+1)*exp( -r*r/2.0/expo(ig+1)/expo(ig+1) )
           enddo
           write(22,*) r, sum
        enddo
        close(unit=22)

        write(filename,'(a,a1,i1,a1,i1,a)') prefix,'_',i,'_',l,'_coeff'
        open(unit=22,file=filename)
        do ig = 0,ng
           write(22,*) psi(ig,i,l+1)
        enddo
        

        close(unit=22)

        


     enddo
  enddo
return
end subroutine dump_gauwf_on_radgrid



subroutine abs_generator_modified(iproc,izatom,ielpsp,psppar,npspcode,ng, noccmax, lmax ,expo,psi, aeval, occup, psp_modifier, &
     Nsol, Labs, Ngrid,Egrid,  rgrid , psigrid  )
  use module_base, only: gp, memocc
  implicit none
  integer, intent(in) :: iproc,izatom,ielpsp,ng,npspcode,noccmax, lmax, Nsol, labs, Ngrid
  real(gp), dimension(0:4,0:6), intent(in) :: psppar
  integer, intent(in) :: psp_modifier
  

  real(gp), dimension(ng+1), intent(out) :: expo

  integer, parameter :: n_int=100

  real(gp), dimension(0:ng,noccmax,lmax+1), intent(out) :: psi, Egrid(Nsol), rgrid(Ngrid), psigrid(Ngrid,Nsol  )
  real(gp), dimension(noccmax,lmax+1  ), intent(out) ::  aeval,occup
  
  !local variables
  character(len=*), parameter :: subname='iguess_generator'
  character(len=27) :: string 
  character(len=2) :: symbol
  real(gp), parameter :: fact=4.0_gp
  integer, dimension(6,4) :: neleconf
  real(gp), dimension(3) :: gpot
  real(gp), dimension(6) :: ott
  real(gp), dimension(noccmax,lmax+1) ::chrg,res
  real(gp), dimension(:), allocatable :: xp,alps
  real(gp), dimension(:,:), allocatable :: vh,hsep,ofdcoef

  real(gp), dimension(:,:,:,:), allocatable :: rmt
  logical :: exists
  integer :: lpx,ncount,nsccode,mxpl,mxchg
  integer :: l,i,j,iocc,il,lwrite,i_all,i_stat
  real(gp) :: alpz,alpl,rcov,rprb,zion,rij,a,a0,a0in,tt,ehomo
  real(gp) value_at_r
  integer :: igrid

  !filename = 'psppar.'//trim(atomname)

  lpx=0
  if (psp_modifier.ne.0) then
     lpx = 0
  else
     lpx_determination: do i=1,4
        if (psppar(i,0) == 0.0_gp) then
           exit lpx_determination
        else
           lpx=i-1
        end if
     end do lpx_determination
  endif
  allocate(alps(lpx+1),stat=i_stat)
  call memocc(i_stat,alps,'alps',subname)
  allocate(hsep(6,lpx+1),stat=i_stat)
  call memocc(i_stat,hsep,'hsep',subname)

  !assignation of radii and coefficients of the local part


  if (psp_modifier.ne.0) then
     alpz=0.001_gp
     alpl=alpz
     alps(1:lpx+1)=0.0_gp
     gpot(1:3)=0.0_gp
     zion= izatom 
  else
     alpz=psppar(0,0)
     alpl=psppar(0,0)
     alps(1:lpx+1)=psppar(1:lpx+1,0)
     gpot(1:3)=psppar(0,1:3)
     zion=real(ielpsp,gp)
  endif

  !assignation of the coefficents for the nondiagonal terms
  if (npspcode == 2) then !GTH case
     do l=1,lpx+1
        hsep(1,l)=psppar(l,1)
        hsep(2,l)=0.0_gp
        hsep(3,l)=psppar(l,2)
        hsep(4,l)=0.0_gp
        hsep(5,l)=0.0_gp
        hsep(6,l)=psppar(l,3)
     end do
  else if (npspcode == 3) then !HGH case
     allocate(ofdcoef(3,4),stat=i_stat)
     call memocc(i_stat,ofdcoef,'ofdcoef',subname)

     ofdcoef(1,1)=-0.5_gp*sqrt(3._gp/5._gp) !h2
     ofdcoef(2,1)=0.5_gp*sqrt(5._gp/21._gp) !h4
     ofdcoef(3,1)=-0.5_gp*sqrt(100.0_gp/63._gp) !h5

     ofdcoef(1,2)=-0.5_gp*sqrt(5._gp/7._gp) !h2
     ofdcoef(2,2)=1._gp/6._gp*sqrt(35._gp/11._gp) !h4
     ofdcoef(3,2)=-7._gp/3._gp*sqrt(1._gp/11._gp) !h5

     ofdcoef(1,3)=-0.5_gp*sqrt(7._gp/9._gp) !h2
     ofdcoef(2,3)=0.5_gp*sqrt(63._gp/143._gp) !h4
     ofdcoef(3,3)=-9._gp*sqrt(1._gp/143._gp) !h5

     ofdcoef(1,4)=0.0_gp !h2
     ofdcoef(2,4)=0.0_gp !h4
     ofdcoef(3,4)=0.0_gp !h5

     !define the values of hsep starting from the pseudopotential file
     do l=1,lpx+1
        hsep(1,l)=psppar(l,1)
        hsep(2,l)=psppar(l,2)*ofdcoef(1,l)
        hsep(3,l)=psppar(l,2)
        hsep(4,l)=psppar(l,3)*ofdcoef(2,l)
        hsep(5,l)=psppar(l,3)*ofdcoef(3,l)
        hsep(6,l)=psppar(l,3)
     end do
     i_all=-product(shape(ofdcoef))*kind(ofdcoef)
     deallocate(ofdcoef,stat=i_stat)
     call memocc(i_stat,i_all,'ofdcoef',subname)
  else if (npspcode == 10) then !HGH-K case
     do l=1,lpx+1
        hsep(1,l)=psppar(l,1) !h11
        hsep(2,l)=psppar(l,4) !h12
        hsep(3,l)=psppar(l,2) !h22
        hsep(4,l)=psppar(l,5) !h13
        hsep(5,l)=psppar(l,6) !h23
        hsep(6,l)=psppar(l,3) !h33
     end do
  end if

  !Now the treatment of the occupation number

  if(psp_modifier.ne.0) then
     call modified_eleconf(izatom,ielpsp,symbol,rcov,rprb,ehomo,neleconf,nsccode,mxpl,mxchg)
  else
     call eleconf(izatom,ielpsp,symbol,rcov,rprb,ehomo,neleconf,nsccode,mxpl,mxchg)
  endif

  occup(:,:)=0.000000000_gp
   do l=0,lmax
     iocc=0
     do i=1,6
        ott(i)=real(neleconf(i,l+1),gp)
        if (ott(i) > 0.0_gp) then
           iocc=iocc+1
            if (iocc > noccmax) stop 'iguess_generator: noccmax too small'
           occup(iocc,l+1)=ott(i)
        endif
     end do

  end do



  !allocate arrays for the gatom routine
  allocate(vh(4*(ng+1)**2,4*(ng+1)**2),stat=i_stat)
  call memocc(i_stat,vh,'vh',subname)

  allocate(xp(0:ng),stat=i_stat)
  call memocc(i_stat,xp,'xp',subname)
  allocate(rmt(n_int,0:ng,0:ng,lmax+1),stat=i_stat)
  call memocc(i_stat,rmt,'rmt',subname)


  !can be switched on for debugging
  !if (iproc.eq.0) write(*,'(1x,a,a7,a9,i3,i3,a9,i3,f5.2)')&
  !     'Input Guess Generation for atom',trim(atomname),&
  !     'Z,Zion=',izatom,ielpsp,'ng,rprb=',ng+1,rprb

  rij=3._gp
  ! exponents of gaussians
  a0in=alpz
  a0=a0in/rij
  !       tt=sqrt(sqrt(2._gp))
  tt=2._gp**.3_gp
  do i=0,ng
     a=a0*tt**i
     xp(i)=.5_gp/a**2
  end do

  ! initial guess
  do l=0,lmax
     do iocc=1,noccmax
        do i=0,ng
           psi(i,iocc,l+1)=0.0_gp
        end do
     end do
  end do


  call crtvh(ng,lmax,xp,vh,rprb,fact,n_int,rmt)



!!$  call gatom(rcov,rprb,lmax,lpx,noccmax,occup,&
!!$       zion,alpz,gpot,alpl,hsep,alps,vh,xp,rmt,fact,n_int,&
!!$       aeval,ng,psi,res,chrg)



  call gatom_modified(rcov,rprb,lmax,lpx,noccmax,occup,&
                 zion,alpz,gpot,alpl,hsep,alps,vh,xp,rmt,fact,n_int,&
                 aeval,ng,psi,res,chrg,&
                 Nsol, Labs, Ngrid,Egrid,  rgrid , psigrid )




  !post-treatment of the inguess data
  do i=1,ng+1
     expo(i)=sqrt(0.5_gp/xp(i-1))
  end do


  do l=0,lmax
     do iocc=1,noccmax
        if( value_at_r(rprb, ng , expo,psi(0,iocc,l+1)).lt.0.0     ) then
           do i=0,ng
              psi(i,iocc,l+1)=-psi(i,iocc,l+1)
           enddo
        endif
     enddo
  enddo
  

  i_all=-product(shape(vh))*kind(vh)
  deallocate(vh,stat=i_stat)
  call memocc(i_stat,i_all,'vh',subname)
  i_all=-product(shape(psi))*kind(psi)

  i_all=-product(shape(xp))*kind(xp)
  deallocate(xp,stat=i_stat)
  call memocc(i_stat,i_all,'xp',subname)
  i_all=-product(shape(rmt))*kind(rmt)
  deallocate(rmt,stat=i_stat)
  call memocc(i_stat,i_all,'rmt',subname)
  i_all=-product(shape(hsep))*kind(hsep)
  deallocate(hsep,stat=i_stat)
  call memocc(i_stat,i_all,'hsep',subname)
  i_all=-product(shape(alps))*kind(alps)
  deallocate(alps,stat=i_stat)
  call memocc(i_stat,i_all,'alps',subname)

END SUBROUTINE abs_generator_modified





subroutine integrate(f,fint,x,Nx)
  use module_base, only: gp
  implicit none
  real(gp), intent(in) :: f(0:Nx-1), x(0:Nx-1)
  real(gp), intent(out):: fint(0:Nx-1)
  integer, intent(in):: Nx
  ! -------------------------
  real(gp) sfin, sgro
  integer i
  real(gp) d,h
  real(gp) p,q,r
  real(gp) fp,fq,fr
  
  
  
  if(Nx.lt.3) then
     print *, " error, Nx<3 in routine integrate "
     stop
  endif
  
  
  sfin=0.0_gp
  sgro=0.0_gp
  fint(0)=0.0_gp
  
  do i=2, Nx-1,2
     sgro=sgro+ (f(i-2)+ f(i))*(x(i)-x(i-2))
     sfin=sfin+(f(i-2) + f(i-1))*(x(i-1)-x(i-2))
     sfin=sfin+(f(i-1)+f(i))*(x(i)-x(i-1))
     fint(i)=(4*sfin-sgro)/6.0_gp
  enddo
  
  
  
  do i=1, Nx-1,2
     
     if( i .lt. Nx-1) then
        d=x(i)-x(i-1)
        h=x(i+1)-x(i)
        
        p=d*(2*d+3*h)/6/(d+h)
        q=d*(d+3*h)/6./h
        r=-d*d*d/6./h/(d+h)
        
        fp=f(i-1)
        fq=f(i)
        fr=f(i+1)
        
        fint(i)=fint(i-1)+p*fp+q*fq+r*fr
     else
        h=x(i-1)-x(i-2)
        d=x(i)-x(i-1)
        p=-d*d*d/6./h/(d+h)
        q=d*(d+3*h)/6./h
        r=d*(2*d+3*h)/6/(d+h)
        
        fp=f(i-2)
        fq=f(i-1)
        fr=f(i)
        fint(i)=fint(i-1)+p*fp+q*fq+r*fr
     endif
     
  enddo
  return
end subroutine integrate


real(gp) function pow(x,n)
  use module_base, only: gp,wp
  real(gp) x
  integer n
  pow=x**n
  return
end function pow

real(gp) function phase( E, N, rgrid, V, nonloc, y,l,normalize,Z, onlyout)
  use module_base, only: gp,wp
  implicit none
  integer N, normalize, onlyout
  real(gp) E,rgrid(N),V(N), nonloc(N), Z, y(N),l
  ! ----------------------------------------------
  integer ii, i,j
  real(gp)  ypi
  
  integer yNcross, yflag
  real(gp)  dh,dl,dh2,dl2,dp,dl3,dhdl2,dh2dl,dh3,add1,add2,deno,num
  real(gp)  Gh,Gl,G0,Nlc,Nla,Nlb
  real(gp)  ca,cb,ep,em,y1,y2
  
  real(gp) Ga,Gb,Gc
  real(gp) ya,yb,yc
  real(gp) r, PI
  
  real(gp) fact,norm, func(N), funcdum(N), pow
  integer count
  
 



  PI=4*atan(1.0d0)
  
  phase=0
  yNcross=0
  yflag=0
  
  ! -- Calcul du point intermediaire où le V est minimum
  if( onlyout.eq.1 ) then
     ii=N-1  
  else
     ii=0
     do i=1 , N-1
        if(V(i+1).ge. V(i) )  then
           ii=i
           if(i.eq.1)  then
              ii=10
              if(ii.ge.N-3) then
                 print *, "  attention !!! if(I>=n-3) in phase";
                 stop
              endif
           endif
           exit
        endif
     enddo
     
     if(ii.eq.0) then 
        ii=N-10;
        print *, " attention !!!I=N-1 in phase  "
        print *, " l est " ,  l
        stop
     endif
     
     
     !  if(I<100) I=100;
     
  endif
  


  ! print * , " propagation de 1 a " , ii
  
  ! //------------ Propagation de  1 a I  ----------------------
  
  do i=1,2
     if(rgrid(i).eq.0.0 .and. l.ne.-1) then
        y(i)=0
     else if (l.eq.-4) then
        y(i)=0.0_gp
     else 
        y(i) = exp((l+1)*log(rgrid(i)))
     endif

  enddo
  
  
  do i=2,ii
     dh = rgrid(i+1)-rgrid(i)
     dl = rgrid(i)-rgrid(i-1)
     G0 = 2*(V(i)  -E)
     Nlb=nonloc(i)
     if(dabs(G0*dh*dl).lt.1.0)  then
        dh2= dh*dh
        dl2= dl*dl
        dp = dh*dl
        
        Gh = 2*(V(i+1)-E)
        Gl = 2*(V(i-1)-E)
        
        Nlc=nonloc(i+1)
        Nla=nonloc(i-1)
        
        R  = dh2 + dl2 + dp 
        
        dh3=dh2*dh
        dhdl2=dl*dp
        dh2dl=dh*dp
        dl3=dl*dl2
        
        ! // ********************************************************************************
        ! // ** 17 luglio 1998. Aggiunta della parte per il termine non locale.
        ! // ** Vedere file numerov.mathematica
        ! // ** ( h^3*(-Nla  + Nlb ) + l^3*(Nlb  - Nlc) + 
        ! // **      h^2*l*(Nla+4*Nlbd+Nlc)+h*l^2*(Nla+4*Nlb+Nlc))/l*(12-(h^2+h*l-l^2)*Gc))
        
        deno=(dl*(1._gp-Gh*(R-2*dl2)/12._gp))
        
        add1=( y(i)*(dh+(dl*Gh*(R-2*dl2)+G0*(dh+dl)*(R+2*dh*dl))/12) &
             -y(i-1)*dh*(1.-Gl*(R-2*dh2)/12.) &
             )/deno;
        
        
        add2=( dh3*(-Nla  + Nlb ) + dl3*(Nlb  - Nlc) + &
             dh2dl*(Nla  + 4*Nlb + Nlc ) + dhdl2*(Nla + 4*Nlb + Nlc) )/deno/6
        
        y(i+1) = y(i)+add1+add2


     else
        print *, " needs taking smaller steps in the grid for the schroedinger equation , in phase"
        stop
     endif
     
     
     if(dabs(y(i+1)).gt.1.0D40 .and. l.ne.-4) then 
        fact= 1./dabs(y(i+1)) 
        count =0
        do j=i+1,1,-1
           if(y(j).eq.0) then
              count=count+1
           else 
              count =0
           endif
           if(count .gt. 2) exit
           y(j)=y(j)*fact 
        enddo
     endif
     if( (i+1) .le.ii ) then !  // LA PHASE EST CALCULE EN I
        if((y(i)*y(i+1)).lt.0) yNcross=yNcross+1
        if(yflag.eq.1 .and. (y(i)*y(i+1)).eq.0)then
           yflag=1-yflag
           yNcross=yNcross+1
        endif
     endif
  enddo
  
  if( y(ii).eq.0.0) then 
     print *, " y[I] == 0.0 dans Schroedinger , in phase"
     stop
  endif
  
  !  // ypI=(y[I+1]*dl2-y[I-1]*dh2+y[I]*(dh2-dl2))/((dh+dl)*dl*dh); // vecchia maniera di calcolare yp
  !  // sostituita con la maniera seguente, presa dal metodo numerov (vedi numerov.nb)
  
  
  
  
  i=ii
  dh = rgrid(i+1)-rgrid(i)
  dl = rgrid(i  )-rgrid(i-1)
  
  Nlc=nonloc(i+1)
  Nla=nonloc(i-1)
  Nlb=nonloc(i)
  
  Gb = 2*(V(i)  -E);
  if(dabs(G0*dh*dl).gt.1.0) then
     print *, " problem with  fabs(G0*dh*dl)>1.0 in calculation  of yp in I , in phase"
     stop
  end if

  Gc = 2*(V(i+1)-E)
  Ga = 2*(V(i-1)-E)
  
  ya=y(ii-1)
  yb=y(ii)
  yc=y(ii+1)
  
  ypI=(dh*(Ga*pow(dl,2)*(-12*dh + Gc*pow(dh,3) - 6*dl + &
       Gc*pow(dh,2)*dl) - &
       6 *(-12*dh + Gc*pow(dh,3) - 12*dl + 2*Gc*pow(dh,2)*dl))*ya + &
       (dh + dl)*(Gb*pow(dl,2)* &
       (-24*dh + 2*Gc*pow(dh,3) - 6*dl + 3*Gc*pow(dh,2)*dl) + &
       6 *(-12*dh + Gc*pow(dh,3) + Gc*pow(dh,2)*dl - &
       Gc*dh*pow(dl,2) + Gc*pow(dl,3)))*yb)/&
       (6.*dh*dl*(dh + dl)*(-12 + Gc*(pow(dh,2) + dh*dl - pow(dl,2)))) 
  
  ! // *****************************************************************************
  ! // ** Termine aggiuntivo
  ! // **
  ! // 6*l^3*Nlc + 
  ! // Nla*(-12*h^2*l - 6*h*l^2 + h^3*l*(h + l)*Gc) + 
  ! // Nlb*(-24*h^2*l - 30*h*l^2 - 6*l^3 +      2*h^3*l*(h + l)*Gc + 3*h^2*l^2*(h + l)*Gc)
  ! // / al denominatore
  ! //
  ! //  6*h*(h + l)*(-12 + (h^2 + h*l - l^2)*Gc)
  
  deno= 6*dh*(dh + dl)*(-12 + (dh*dh + dh*dl - dl*dl)*Gc)
	    
  num=  6*dl*dl*dl*Nlc&
       +Nla*(-12*dh*dh*dl-6*dh*dl*dl+dh*dh*dh*dl*(dh + dl)*Gc)&
       +Nlb*(-24*dh*dh*dl - 30*dh*dl*dl - &
       6*dl*dl*dl+2*dh*dh*dh*dl*(dh + dl)*Gc +3*dh*dh*dl*dl*(dh+dl)*Gc)

  ypI=ypI+num/deno
	    
  
  if( onlyout.eq.1) then
     ypI = ypI+ dh*(  y(ii)*Gb  +y(ii+1)*Gc )/2 
     phase=ypI
     return
  endif


	  
	  
  if(dabs(y(ii)) .gt. dabs( ypI ) ) then
     phase=-atan(ypI/y(ii))
  else
     r = y(ii)/ypI ;
     if(  r .gt.0.  ) then
        phase=-( PI/2.0 - atan(r) )
     else
        phase=- ( - PI/2. - atan( r ));
     endif
  endif



  
  if( dabs(y(ii) ) .gt. dabs(ypI)) then
     y(ii+1)=y(ii+1)/y(ii)
     do i=1, ii
        y(i)=y(i)/y(ii)
     enddo
  else
     do i=1, ii+1
        y(i)=y(i)/ypI
     enddo
  endif
  
  ! //------------ Propagation de   I   à  rinf --------------------------
  
  do i=N,N-1,-1 
     y(i)=N-i    ! //y[i]=exp(-sqrt(-2*E)*r[i]);
  enddo
  
  do i=N-1, ii,-1
     dh = rgrid(i)-rgrid(i-1)
     dl = rgrid(i+1)-rgrid(i)
     G0 = 2*(V(i)  -E)
     Nlb=nonloc(i)
     
     if(dabs(G0*dh*dl).lt.1.0) then
        dh2= dh*dh
        dl2=dl*dl
        dp = dh*dl
        
        Gh = 2*(V(i-1)-E)
        Gl = 2*(V(i+1)-E)
        
        Nlc=nonloc(i-1)
        Nla=nonloc(i+1)
		  
        R  = dh2 + dl2 + dp 
		  
        dh3=dh2*dh
        dhdl2=dl*dp
        dh2dl=dh*dp
        dl3=dl*dl2
        
        ! // ********************************************************************************
        ! // ** 17 luglio 1998. Aggiunta della parte per il termine non locale.
        ! // ** Vedere file numerov.mathematica
        ! // ** ( h^3*(-Nla  + Nlb ) + l^3*(Nlb  - Nlc) + 
        ! // **      h^2*l*(Nla+4*Nlbd+Nlc)+h*l^2*(Nla+4*Nlb+Nlc))/l*(12-(h^2+h*l-l^2)*Gc))
        
        deno=(dl*(1.-Gh*(R-2*dl2)/12.))
        
        add1=( y(i)*(dh+(dl*Gh*(R-2*dl2)+G0*(dh+dl)*(R+2*dh*dl))/12)&
             -y(i+1)*dh*(1.-Gl*(R-2*dh2)/12.) &
             )/deno
        
        add2=( dh3*(-Nla  + Nlb ) + dl3*(Nlb  - Nlc) + &
             dh2dl*(Nla  + 4*Nlb + Nlc ) + dhdl2*(Nla + 4*Nlb + Nlc) )/deno
		  
        y(i-1) = y(i)+add1+add2
        
        ! /*  
        ! y[i-1] = y[i]+( y[i]*(dh+(dl*Gh*(R-2*dl2)+G0*(dh+dl)*(R+2*dh*dl))/12)
        ! -y[i+1]*dh*(1.-Gl*(R-2*dh2)/12.)
        ! )/
        ! (dl*(1.-Gh*(R-2*dl2)/12.));
        ! */
     else
        print *, "needs taking smaller steps in the grid for the schroedinger equation  , in phase "
        stop
     endif
 
     

     if(dabs(y(i-1)).gt.1.0D8) then 
        fact= 1./dabs(y(i-1)) 
        count =0
        do j=i-1, N 
           if(y(j).eq.0) then
              count=count+1
           else 
              count =0
           endif
           if(count.ne.0) exit
           y(j)=y(j)*fact 
        enddo
     endif

     if( (i-1) .ge.ii ) then!  // LA PHASE EST CALCULE EN I
        if((y(i)*y(i-1) ).lt.0) yNcross=yNcross+1
        if(yflag.eq.1 .and. (y(i)*y(i-1)).eq.0) then
           yflag=1-yflag
           yNcross=yNcross+1
        endif
     endif
  enddo
  !  // ypI=(y[I+1]*dh2-y[I-1]*dl2+y[I]*(dl2-dh2))/((dh+dl)*dl*dh);
  !  // sostituita con la maniera seguente, presa dal metodo numerov (vedi numerov.nb)
  ! {

  i=ii
  dh = rgrid(i+1)-rgrid(i)
  dl = rgrid(i  )-rgrid(i-1)
  
  Nlc=nonloc(i+1)
  Nla=nonloc(i-1)
  Nlb=nonloc(i)
  
  
  Gb = 2*(V(i)  -E)
  if(dabs(G0*dh*dl).gt.1.0) then
     print *, " problem with fabs(G0*dh*dl)>1.0 in calculation di yp in I , in phase"
     stop
  endif
  Gc = 2*(V(i+1)-E)
  Ga = 2*(V(i-1)-E)

  ya=y(ii-1)
  yb=y(ii)
  yc=y(ii+1)
	    
  ypI=(dh*(Ga*pow(dl,2)*(-12*dh + Gc*pow(dh,3) - 6*dl + &
       Gc*pow(dh,2)*dl) - &
       6*(-12*dh + Gc*pow(dh,3) - 12*dl + 2*Gc*pow(dh,2)*dl))*ya + &
       (dh + dl)*(Gb*pow(dl,2)*&
       (-24*dh + 2*Gc*pow(dh,3) - 6*dl + 3*Gc*pow(dh,2)*dl) + &
       6*(-12*dh + Gc*pow(dh,3) + Gc*pow(dh,2)*dl - &
       Gc*dh*pow(dl,2) + Gc*pow(dl,3)))*yb)/&
       (6.*dh*dl*(dh + dl)*(-12 + Gc*(pow(dh,2) + dh*dl - pow(dl,2)))) 
	    
  ! // *****************************************************************************
  !	    // ** Termine aggiuntivo
  !	    // **
  !	    // 6*l^3*Nlc + 
  !	    // Nla*(-12*h^2*l - 6*h*l^2 + h^3*l*(h + l)*Gc) + 
  !	    // Nlb*(-24*h^2*l - 30*h*l^2 - 6*l^3 +      2*h^3*l*(h + l)*Gc + 3*h^2*l^2*(h + l)*Gc)
  !	    // / al denominatore
  !	    //
  !	    //  6*h*(h + l)*(-12 + (h^2 + h*l - l^2)*Gc)
	    
  deno= 6*dh*(dh + dl)*(-12 + (dh*dh + dh*dl - dl*dl)*Gc)
	    
  num=  6*dl*dl*dl*Nlc &
       +Nla*(-12*dh*dh*dl-6*dh*dl*dl+dh*dh*dh*dl*(dh + dl)*Gc)&
       +Nlb*(-24*dh*dh*dl - 30*dh*dl*dl - &
       6*dl*dl*dl+2*dh*dh*dh*dl*(dh + dl)*Gc +3*dh*dh*dl*dl*(dh+dl)*Gc&
       )
	    
  ypI=ypI+num/deno
	    
  !}
	
  phase=phase+PI*yNcross 
	

 
  if(dabs(y(ii)) .gt. dabs( ypI ) ) Then
     phase =phase+atan(ypI/y(ii))
  else
     r = y(ii)/ypI 
     if(  r .gt.0.  ) then
        phase =phase+( PI/2. - atan(r) )
     else
        phase = phase+ ( - PI/2. - atan( r ))
     endif
  endif
	  
	  
  if( dabs(y(ii)) .gt. dabs(ypI)) then
     y(ii-1)=y(ii-1)/y(ii)
     do i=N,ii,-1
        y(i)=y(i)/y(ii)
     enddo
  else
     do i=N,ii-1,-1
        y(i)=y(i)/ypI
     enddo
  endif
	  
	  
  ! //------------------------------------------------------------------------------
  if(normalize.eq.1) then
     do i=1,N
        func(i)=y(i)*y(i)
     enddo
     
     call integrate(func,funcdum,rgrid,N)
     norm=sqrt(funcdum(N) )
     do i=1,N
        y(i)=y(i)/norm
     enddo
  endif
  return 
end function phase
	
 

subroutine schro(E, r,  V,nonloc, y, NGRID, nsol, l,  Z)
  use module_base, only: gp,wp
  implicit none
  integer, intent(IN) :: nsol,ngrid
  real(gp), intent(IN) ::  r(ngrid),v(ngrid), nonloc(ngrid)
  real(gp), intent(OUT) ::  E
  real(gp), intent(OUT) ::  y(ngrid)
  real(gp), intent(IN) ::  l,Z
  ! -------------------------------
  real(gp) Elow, Ehigh, Eguess
  real(gp) pathh, pathl, fase
  integer i, igrid
  real(gp) Phase, but
  real(gp) PI



  PI=4.0*atan(1.0)
  Elow=0;

  do i=1, NGRID
      if(V(i) .lt. Elow) Elow=V(i)
   enddo
  if( Elow .lt. -Z*Z) Elow = -Z*Z
  Elow=-Z*Z

  Ehigh = 0.0;

  

  


  
  pathh = Phase(Ehigh,NGRID,r,v,nonloc,y,  l ,0, Z,0);
  pathl = Phase(Elow ,NGRID, r,v,nonloc,y,  l ,0,  Z,0);
 

!!$  print *, Ehigh, pathh
!!$  print *, Elow, pathl

  but= PI*(nsol-l-1)

!!$  print *, " Z " , Z



  if( pathl.gt.but)  then
     print *, " pathl>but " 
     print *, " Elow " , Elow
     print *, " now exiting , routine schro" 
     stop
  endif
  
      
  if(but .gt. pathh) then
     Ehigh = (but+1)*(but+1)/r(NGRID-1)/r(NGRID-1)
     do while( but .gt. Phase(Ehigh ,NGRID, r,V,nonloc,y,  l , 0, Z,0 )	 )
         Ehigh =2*Ehigh;
      enddo
   endif
   

 
   do while(1.gt.0) 

      Eguess = (Elow+Ehigh)/2


      fase=  Phase(Eguess ,NGRID, r,V,nonloc,y,  l , 0, Z ,0)
      



      if( fase.gt.but) then
         Ehigh=Eguess
      else
         Elow=Eguess
      endif

      if(dabs(but-fase).lt.1.0e-8) exit
   enddo

   
   fase  = Phase(Eguess,NGRID, r,v,nonloc,y, l,1 ,  Z  ,0)
   E=Eguess;



   return
 end subroutine schro











subroutine gatom_modified(rcov,rprb,lmax,lpx,noccmax,occup,&
                 zion,alpz,gpot,alpl,hsep,alps,vh,xp,rmt,fact,nintp,&
                 aeval,ng,psi,res,chrg,&
                 Nsol, Labs, Ngrid,Egrid,  rgrid , psigrid )
  use module_base, only: gp
  implicit real(gp) (a-h,o-z)
  logical :: noproj
  integer, parameter :: n_int=100
  dimension psi(0:ng,noccmax,lmax+1),aeval(noccmax,lmax+1),&
       hh(0:ng,0:ng),ss(0:ng,0:ng),eval(0:ng),evec(0:ng,0:ng),&
       aux(2*ng+2),&
       gpot(3),hsep(6,lpx+1),rmt(n_int,0:ng,0:ng,lmax+1),&
       pp1(0:ng,lpx+1),pp2(0:ng,lpx+1),pp3(0:ng,lpx+1),alps(lpx+1),&
       potgrd(n_int),&
       rho(0:ng,0:ng,lmax+1),rhoold(0:ng,0:ng,lmax+1),xcgrd(n_int),&
       occup(noccmax,lmax+1),chrg(noccmax,lmax+1),&
       vh(0:ng,0:ng,4,0:ng,0:ng,4),&
       res(noccmax,lmax+1),xp(0:ng),& 
       rgrid(Ngrid), psigrid(Ngrid, Nsol),psigrid_naked(Ngrid,Nsol), projgrid(Ngrid,3), &
       rhogrid(Ngrid), potgrid(Ngrid), &
       vxcgrid(Ngrid), &
       dumgrid1(Ngrid),dumgrid2(Ngrid),  &
       Egrid(nsol), ppgrid(Nsol,3), work(nsol*nsol*2), &
       H(Nsol, Nsol)

  if (nintp.ne.n_int) stop 'n_int><nintp'


  do l=0,lmax
     if (occup(1,l+1).gt.0._gp) lcx=l
  end do
  !write(6,*) 'lcx',lcx
 
  noproj=.true.
  do l=1,lpx+1
     noproj = noproj .and. (alps(l) .eq. 0._gp)
  end do


  

! projectors, just in case
  if (.not. noproj) then
     do l=0,lpx
        gml1=sqrt( gamma(real(l,gp)+1.5_gp) / (2._gp*alps(l+1)**(2*l+3)) )
        gml2=sqrt( gamma(real(l,gp)+3.5_gp) / (2._gp*alps(l+1)**(2*l+7)) )&
            /(real(l,gp)+2.5_gp)
        gml3=sqrt( gamma(real(l,gp)+5.5_gp) / (2._gp*alps(l+1)**(2*l+11)) )&
            /((real(l,gp)+3.5_gp)*(real(l,gp)+4.5_gp))
        tt=1._gp/(2._gp*alps(l+1)**2)
        do i=0,ng
           ttt=1._gp/(xp(i)+tt)
           pp1(i,l+1)=gml1*(sqrt(ttt)**(2*l+3))
           pp2(i,l+1)=gml2*ttt*(sqrt(ttt)**(2*l+3))
           pp3(i,l+1)=gml3*ttt**2*(sqrt(ttt)**(2*l+3))
        end do
        
        if(l.eq.Labs) then
           do iorder=1,3
              do igrid=1,Ngrid
                 r=rgrid(igrid)
                 projgrid(igrid,iorder) = (  exp( -r*r*tt)*(r**l *r) )* r**(2*(i-1))
              end do
              
              do igrid=1, Ngrid
                 dumgrid1(igrid) =   projgrid(igrid,iorder)*projgrid(igrid,iorder)
              enddo
              call integrate( dumgrid1, dumgrid2, rgrid, Ngrid)
              do igrid=1, Ngrid
                 projgrid(igrid,iorder) = projgrid(igrid,iorder) / sqrt( dumgrid2(Ngrid) ) 
              enddo
           enddo
        endif
     end do
  else
     pp1(:,:)=0._gp
     pp2(:,:)=0._gp
     pp3(:,:)=0._gp
     projgrid(:,:)=0.0_gp
  end if

  do l=0,lmax
     do j=0,ng
        do i=0,ng
           rho(i,j,l+1)=0._gp
        end do
     end do
  end do


  do igrid=1,Ngrid
     rhogrid(igrid)=0.0_gp
  enddo

  evsum=1.d30
  big_loop: do it=1,50
     evsumold=evsum
     evsum=0._gp
     
! coefficients of charge density
     do l=0,lmax
        do j=0,ng
           do i=0,ng
              rhoold(i,j,l+1)=rho(i,j,l+1)
              rho(i,j,l+1)=0._gp        
           end do
        end do
     end do

     do l=0,lmax
        do iocc=1,noccmax
           if (occup(iocc,l+1).gt.0._gp) then
              do j=0,ng
                 do i=0,ng
                    rho(i,j,l+1)=rho(i,j,l+1) + &
                         psi(i,iocc,l+1)*psi(j,iocc,l+1)*occup(iocc,l+1)
                 end do
              end do
           end if
        end do

       
     end do

     do igrid=1,Ngrid
        r=rgrid(igrid)
        rhogrid(igrid)=0.0_gp
        do l=0, lmax
           do iocc=1,noccmax
              if (occup(iocc,l+1).gt.0._gp) then
                 dum=0.0_gp
                 do j=0,ng
                    dum=dum + psi(j,iocc,l+1)*exp(-r*r *   xp(j) )* r**(l+1)
                 end do
                 rhogrid(igrid)=rhogrid(igrid)+dum*dum*occup(iocc,l+1)
              end if
           end do
        enddo
        dum = rhogrid(igrid)/r/r *0.07957747154594768_gp
        vxcgrid(igrid)=emuxc(dum)
     enddo
     
  
     rmix=.5_gp
     if (it.eq.1) rmix=1._gp
     do l=0,lmax
        do j=0,ng
           do i=0,ng
              tt=rmix*rho(i,j,l+1) + (1._gp-rmix)*rhoold(i,j,l+1)
              rho(i,j,l+1)=tt
           end do
        end do
     end do

! XC potential on grid
!        do k=1,n_int
!           xcgrd(k)=0._gp
!        end do
!        do l=0,lmax
!           do j=0,ng
!              do i=0,ng
!                 do k=1,n_int
!                    xcgrd(k)=xcgrd(k)+rmt(k,i,j,l+1)*rho(i,j,l+1)
!                 end do
!              end do
!           end do
!        end do
     call DGEMV('N',n_int,(lcx+1)*(ng+1)**2,1._gp,&
                rmt,n_int,rho,1,0._gp,xcgrd,1)

     dr=fact*rprb/real(n_int,gp)
     do k=1,n_int
        r=(real(k,gp)-.5_gp)*dr
! divide by 4 pi
        tt=xcgrd(k)*0.07957747154594768_gp
! multiply with r^2 to speed up calculation of matrix elements
        xcgrd(k)=emuxc(tt)*r**2
     end do

     do igrid=1, Ngrid
        r=rgrid(igrid)
        potgrid(igrid) =0.5_gp*r*r  /    rprb**4 
        potgrid(igrid) = potgrid(igrid) - zion/r * derf( r/alpz/sqrt(2.0)   )
        rr = r/alpz
        potgrid(igrid) = potgrid(igrid) + exp(-0.5_gp * rr**2 )*( gpot(1)+gpot(2)*rr**2 + gpot(3)*rr**4 )
     enddo
! poisson per potgrid
     call integrate( rhogrid, dumgrid1, rgrid, Ngrid)
     do igrid=1, Ngrid
        potgrid(igrid)=potgrid(igrid)+dumgrid1(igrid)/rgrid(igrid)
     enddo

     do igrid=1, Ngrid
        dumgrid1(igrid) = rhogrid(igrid)/rgrid(igrid)
     enddo
     call integrate( dumgrid1, dumgrid2, rgrid, Ngrid)
     do igrid=1, Ngrid
        potgrid(igrid)=potgrid(igrid)-dumgrid2(igrid) + dumgrid2(Ngrid)
     enddo
! ------------------------------------------------------------------------------------



 

     loop_l: do l=0,lmax
        gml=.5_gp*gamma(.5_gp+real(l,gp))

!  lower triangles only
        loop_i: do i=0,ng
           loop_j: do j=0,i
              d=xp(i)+xp(j)
              sxp=1._gp/d
              const=gml*sqrt(sxp)**(2*l+1)
! overlap
              ss(i,j)=const*sxp*(real(l,gp)+.5_gp)
! kinetic energy
              hh(i,j)=.5_gp*const*sxp**2* ( 3._gp*xp(i)*xp(j) +&
                   real(l,gp)*(6._gp*xp(i)*xp(j)-xp(i)**2-xp(j)**2) -&
                   real(l,gp)**2*(xp(i)-xp(j))**2  ) + .5_gp*real(l,gp)*(real(l,gp)+1._gp)*const
! potential energy from parabolic potential
              hh(i,j)=hh(i,j) +&
                   .5_gp*const*sxp**2*(real(l,gp)+.5_gp)*(real(l,gp)+1.5_gp)/rprb**4 
! hartree potential from ionic core charge
              tt=sqrt(1._gp+2._gp*alpz**2*d)
              if (l.eq.0) then
                 hh(i,j)=hh(i,j) -zion/(2._gp*d*tt)
              else if (l.eq.1) then
                 hh(i,j)=hh(i,j) -zion* &
                      (1._gp + 3._gp*alpz**2*d)/(2._gp*d**2*tt**3)
              else if (l.eq.2) then
                 hh(i,j)=hh(i,j) -zion* &
                      (2._gp + 10._gp*alpz**2*d + 15._gp*alpz**4*d**2)/(2._gp*d**3*tt**5)
              else if (l.eq.3) then
                 hh(i,j)=hh(i,j) -zion*3._gp* &
                      (2._gp +14._gp*alpz**2*d +35._gp*alpz**4*d**2 +35._gp*alpz**6*d**3)/&
                      (2._gp*d**4*tt**7)
              else 
                 stop 'l too big'
              end if
! potential from repulsive gauss potential
              tt=alpl**2/(.5_gp+d*alpl**2)
              if (1.eq.1) then
              hh(i,j)=hh(i,j)+ gpot(1)*.5_gp*gamma(1.5_gp+real(l,gp))*tt**(1.5_gp+real(l,gp))&
                   + (gpot(2)/alpl**2)*.5_gp*gamma(2.5_gp+real(l,gp))*tt**(2.5_gp+real(l,gp))&
                   + (gpot(3)/alpl**4)*.5_gp*gamma(3.5_gp+real(l,gp))*tt**(3.5_gp+real(l,gp))
           endif
! separable terms
              if (1.eq.1 .and. l.le.lpx) then
                 hh(i,j)=hh(i,j) + pp1(i,l+1)*hsep(1,l+1)*pp1(j,l+1)&
                      + pp1(i,l+1)*hsep(2,l+1)*pp2(j,l+1)&
                      + pp2(i,l+1)*hsep(2,l+1)*pp1(j,l+1)&
                      + pp2(i,l+1)*hsep(3,l+1)*pp2(j,l+1)&
                      + pp1(i,l+1)*hsep(4,l+1)*pp3(j,l+1)&
                      + pp3(i,l+1)*hsep(4,l+1)*pp1(j,l+1)&
                      + pp2(i,l+1)*hsep(5,l+1)*pp3(j,l+1)&
                      + pp3(i,l+1)*hsep(5,l+1)*pp2(j,l+1)&
                      + pp3(i,l+1)*hsep(6,l+1)*pp3(j,l+1)
              end if
! hartree potential from valence charge distribution
!              tt=0._gp
!              do lp=0,lcx
!                 do jp=0,ng
!                    do ip=0,ng
!                       tt=tt + vh(ip,jp,lp+1,i,j,l+1)*rho(ip,jp,lp+1)
!                    end do
!                 end do
!              end do
              tt=DDOT((lcx+1)*(ng+1)**2,vh(0,0,1,i,j,l+1),1,rho(0,0,1),1)
              hh(i,j)=hh(i,j) + tt
! potential from XC potential
              dr=fact*rprb/real(n_int,gp)
!              tt=0._gp
!              do k=1,n_int
!                 tt=tt+xcgrd(k)*rmt(k,i,j,l+1)
!              end do
              tt=DDOT(n_int,rmt(1,i,j,l+1),1,xcgrd(1),1)
              hh(i,j)=hh(i,j)+tt*dr
           end do loop_j
        end do loop_i

! ESSL
!        call DSYGV(1,hh,ng+1,ss,ng+1,eval,evec,ng+1,ng+1,aux,2*ng+2)
! LAPACK
        call DSYGV(1,'V','L',ng+1,hh,ng+1,ss,ng+1,eval,evec,(ng+1)**2,info)
        if (info.ne.0) write(6,*) 'LAPACK',info
        do iocc=0,noccmax-1
           do i=0,ng
              evec(i,iocc)=hh(i,iocc)
           end do
        end do
! end LAPACK
        do iocc=1,noccmax
           evsum=evsum+eval(iocc-1)
           aeval(iocc,l+1)=eval(iocc-1)
           do i=0,ng
              psi(i,iocc,l+1)=evec(i,iocc-1)
           end do
        end do
!        write(6,*) 'eval',l
!55      format(5(e14.7))
!        write(6,55) eval 
!        write(6,*) 'diff eval'
!        write(6,55) (eval(i)-eval(i-1),i=1,ng)
!        write(6,*) 'evec',l
!33      format(10(e9.2))
!        do i=0,ng
!           write(6,33) (evec(i,iocc),iocc=0,noccmax-1)
!        end do

     end do loop_l

     tt=abs(evsum-evsumold)
!     write(6,*) 'evdiff',it,tt
     if (tt.lt.1.e-12_gp) then
         exit big_loop
     end if
  end do big_loop
! End of the big loop




  dumgrid1(:)=0.0_gp
  do igrid=1, ngrid
     r=rgrid(igrid)
     potgrid(igrid)=potgrid(igrid)+ 0.5_gp*labs*(labs+1.0_gp)/r/r
  enddo
  


  do isol=1,nsol

     call schro(Egrid(isol) , rgrid ,  potgrid , dumgrid1, psigrid_naked(:,isol) , ngrid , isol+labs , labs*1.0_gp ,  zion)

  enddo
  

  H(:,:)=0.0D0
  do i=1,Nsol
     H(i,i)=Egrid(i)
     do iproj=1,3
        do igrid=1,Ngrid
           dumgrid1(igrid)=psigrid_naked(igrid,i)*projgrid(igrid,iproj)
        enddo
        call integrate(dumgrid1,dumgrid2,rgrid,Ngrid)
        ppgrid(i,iproj)=dumgrid2(Ngrid)
     enddo
  enddo

  do i=1,Nsol
     do j=1, Nsol
        if ( labs.le.lpx) then
           H(i,j)=H(i,j)+ ppgrid(i,1)*hsep(1,labs+1)*ppgrid(j,1)&
                + ppgrid(i,1)*hsep(2,labs+1)*ppgrid(j,2)&
                + ppgrid(i,2)*hsep(2,labs+1)*ppgrid(j,1)&
                + ppgrid(i,2)*hsep(3,labs+1)*ppgrid(j,2)&
                + ppgrid(i,1)*hsep(4,labs+1)*ppgrid(j,3)&
                + ppgrid(i,3)*hsep(4,labs+1)*ppgrid(j,1)&
                + ppgrid(i,2)*hsep(5,labs+1)*ppgrid(j,3)&
                + ppgrid(i,3)*hsep(5,labs+1)*ppgrid(j,2)&
                + ppgrid(i,3)*hsep(6,labs+1)*ppgrid(j,3)
        endif
        do igrid=1,Ngrid
           dumgrid1(igrid)=psigrid_naked(igrid,i)*psigrid_naked(igrid,j)*vxcgrid(igrid)
        enddo
        call integrate(dumgrid1,dumgrid2,rgrid,Ngrid)
        H(i,j)=H(i,j)+dumgrid2(Ngrid)
     enddo
  enddo
  



  call DSYEV('V','U', Nsol, H, Nsol,Egrid , WORK, Nsol*Nsol*2, INFO)
  
  call  DGEMM('N','N',Ngrid ,Nsol,   Nsol,1.0d0 ,psigrid_naked, Ngrid ,H,Nsol, 0.0D0 , psigrid , Ngrid)

  call resid(lmax,lpx,noccmax,rprb,xp,aeval,psi,rho,ng,res,&
             zion,alpz,alpl,gpot,pp1,pp2,pp3,alps,hsep,fact,n_int,&
             potgrd,xcgrd,noproj)





! charge up to radius rcov
  if (lmax.gt.3) stop 'cannot calculate chrg'
  do l=0,lmax
     do iocc=1,noccmax
        chrg(iocc,l+1)=0._gp
     end do
  end do

  do iocc=1,noccmax
     do j=0,ng
        do i=0,ng
           d=xp(i)+xp(j)
           sd=sqrt(d)
           terf=derf(sd*rcov) 
           texp=exp(-d*rcov**2)

           tt=0.4431134627263791_gp*terf/sd**3 - 0.5_gp*rcov*texp/d
           chrg(iocc,1)=chrg(iocc,1) + psi(i,iocc,1)*psi(j,iocc,1)*tt
           if (lmax.eq.0) then
              cycle
           end if
           tt=0.6646701940895686_gp*terf/sd**5 + &
              (-0.75_gp*rcov*texp - 0.5_gp*d*rcov**3*texp)/d**2
           chrg(iocc,2)=chrg(iocc,2) + psi(i,iocc,2)*psi(j,iocc,2)*tt
           if (lmax.eq.1) then
               cycle
           end if
           tt=1.661675485223921_gp*terf/sd**7 + &
              (-1.875_gp*rcov*texp-1.25_gp*d*rcov**3*texp-.5_gp*d**2*rcov**5*texp) &
              /d**3
           chrg(iocc,3)=chrg(iocc,3) + psi(i,iocc,3)*psi(j,iocc,3)*tt
           if (lmax.eq.2) then
              cycle
           end if
           tt=5.815864198283725_gp*terf/sd**9 + &
              (-6.5625_gp*rcov*texp - 4.375_gp*d*rcov**3*texp - &
              1.75_gp*d**2*rcov**5*texp - .5_gp*d**3*rcov**7*texp)/d**4
           chrg(iocc,4)=chrg(iocc,4) + psi(i,iocc,4)*psi(j,iocc,4)*tt
        end do
     end do
  end do



! ------------------------------------------------
  


! -----------------------------------------------




! writing lines suppressed
!!$        write(66,*)  lmax+1
!!$        write(66,*) ' #LINETYPE{1324}' 
!!$        write(66,*) ' $' 
!!$  do l=0,lmax
!!$           write(66,*) ' 161'
!!$     r=0._gp
!!$     do
!!$        tt= wave(ng,l,xp,psi(0,1,l+1),r)
!!$              write(66,*) r,tt
!!$        r=r+.025_gp
!!$        if(r > 4.00001_gp) exit
!!$     end do
!!$  end do
! writing lines suppressed
!!$        write(67,*) min(lmax+1,3)
!!$        write(67,*) ' #LINETYPE{132}'
!!$        write(67,*) ' #TITLE{FOURIER}' 
!!$        write(67,*) ' $'
  dr=6.28_gp/rprb/200._gp
!!$        write(67,*) ' 200'
  rk=0._gp
  loop_rk1: do 
     tt=0._gp
     do i=0,ng
        texp=exp(-.25_gp*rk**2/xp(i))
!        texp=exp(-.5_gp*energy/xp(i))
        sd=sqrt(xp(i))
        tt=tt+psi(i,1,1)*0.4431134627263791_gp*texp/sd**3
     end do
!!$           write(67,*) rk,tt
     rk=rk+dr
     if(rk > 6.28_gp/rprb-.5_gp*dr) exit loop_rk1
  end do loop_rk1
  if (lmax.ge.1) then
!!$           write(67,*) ' 200'
     rk=0._gp
     loop_rk2: do 
        tt=0._gp
        do i=0,ng
           texp=exp(-.25_gp*rk**2/xp(i))
           sd=sqrt(xp(i))
           tt=tt+psi(i,1,2)*0.2215567313631895_gp*rk*texp/sd**5
        end do
!!$              write(67,*) rk,tt
        rk=rk+dr
        if (rk > 6.28_gp/rprb-.5_gp*dr) exit loop_rk2
     end do loop_rk2
  end if
  if (lmax.ge.2) then
!!$           write(67,*) ' 200'
     rk=0._gp
     do 
        tt=0._gp
        do i=0,ng
           texp=exp(-.25_gp*rk**2/xp(i))
           sd=sqrt(xp(i))
           tt=tt+psi(i,1,3)*0.1107783656815948_gp*rk**2*texp/sd**7
        end do
!!$              write(67,*) rk,tt
        rk=rk+dr
        if (rk > 6.28_gp/rprb-.5_gp*dr) exit
     end do
  end if

END SUBROUTINE gatom_modified


subroutine GetExcitedOrbitalAsG( in_iat_absorber ,Gabsorber, atoms, rxyz, nproc, iproc, dump_functions,Gabs_coeffs)

  use module_base
  use module_types
  use module_interfaces

  implicit none
  integer, parameter :: abs_final_L = 1
  integer, intent(in) :: in_iat_absorber, nproc, iproc, dump_functions
  type(gaussian_basis) , intent(out) :: Gabsorber
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
  real(wp), dimension(2*abs_final_L+1), intent(out) :: Gabs_coeffs
  
  ! -----------------------------------------------------------
  

  integer :: ity, ng , noccmax, lmax, j, ierr, i_all
  real(gp) , pointer :: expo(:), psi(:,:,:), aeval(:,:), occup(:,:), gcoeffs(:)
  integer :: psp_modifier
  integer :: ig, nord, iocc, iexpo

  integer ::  abs_initial_L
  integer, parameter :: Norder=4
  real(gp) :: Scoeffs(Norder)
  integer :: iocc_for_j(  Norder )
  real(gp) :: ene_for_j(  Norder )
  real(gp) , parameter :: sphere_radius=3.0
  real(gp) , pointer:: psi1s(:) 
  integer :: iw
  integer :: real_start
  real(gp) :: cradius
  
  integer ::  Nsol , Ngrid, igrid

  integer :: ng_fine
  real(gp), pointer :: expo_fine(:)

  real(gp), pointer :: Egrid(:) ,  rgrid(:) , psigrid (:,:) , Egrid_pseudo(:) ,  psigrid_pseudo (:,:) 
  integer i_stat
  character(len=*), parameter :: subname='GetExcitedOrbitalAsG'


  ! if (in_iat_absorber.ne.0) then

  ity = atoms%iatype(in_iat_absorber)
  ng  = 30
  noccmax = 5 
  lmax=3
  
  ng_fine= 200
  
  Nsol=50
  Ngrid=10000
  
  cradius=4.0
  

  
  allocate(expo_fine(ng_fine  +ndebug ), stat=i_stat)
  call memocc(i_stat,expo_fine,'expo_fine',subname)
  
  allocate(expo(ng +ndebug  ), stat=i_stat)
  call memocc(i_stat,expo,'expo',subname)
  
  allocate(psi ( 0:ng-1  ,noccmax,lmax+1+ndebug ), stat=i_stat)
  call memocc(i_stat,psi,'psi',subname)
  
  
  allocate(gcoeffs ( 0:ng_fine-1  +ndebug ), stat=i_stat)
  call memocc(i_stat,gcoeffs,'gcoeffs',subname)
  
  allocate(aeval ( noccmax  ,lmax+1+ndebug ), stat=i_stat)
  call memocc(i_stat,aeval,'aeval',subname)
  
  allocate(occup ( noccmax  ,lmax+1+ndebug ), stat=i_stat)
  call memocc(i_stat,occup,'occup',subname)
  
  allocate( Egrid(Nsol +ndebug ), stat=i_stat)
  call memocc(i_stat,Egrid,'Egrid',subname)
  
  allocate( rgrid(Ngrid +ndebug ), stat=i_stat)
  call memocc(i_stat,rgrid,'rgrid',subname)
  
  allocate( psigrid(Ngrid  , Nsol +ndebug ), stat=i_stat)
  call memocc(i_stat,psigrid,'psigrid',subname)
  
  allocate( Egrid_pseudo(Nsol +ndebug ), stat=i_stat)
  call memocc(i_stat,Egrid_pseudo,'Egrid_pseudo',subname)
  
  allocate( psigrid_pseudo(Ngrid  , Nsol +ndebug), stat=i_stat)
  call memocc(i_stat,psigrid_pseudo,'psigrid_pseudo',subname)

  allocate(psi1s( Ngrid +ndebug ), stat=i_stat)
  call memocc(i_stat,psi1s,'psi1s',subname)
   
  do igrid=1, Ngrid
     rgrid(igrid) = igrid*1.0_gp/Ngrid * cradius
  enddo
  
  

  abs_initial_L = 0
  
 
  
  if(iproc.eq.0)   print * , " routine GetExcitedOrbitalAsG  , calculate pseudo " 
  psp_modifier=0
  call abs_generator_modified(iproc,atoms%nzatom(ity), atoms%nelpsp(ity),atoms%psppar(0,0,ity),&
       atoms%npspcode(ity),ng-1 ,noccmax , lmax , expo,psi,aeval, occup , psp_modifier , &
       Nsol, abs_final_L , Ngrid,Egrid_pseudo,  rgrid , psigrid_pseudo  )

  

  if(iproc.eq.0 .and. dump_functions.eq.1) then
     if(psp_modifier.eq.0) then
        call dump_gauwf_on_radgrid("pseudo_wf_radgrid", ng-1 ,noccmax , lmax , expo,psi,aeval, occup     )
     else
        call dump_gauwf_on_radgrid("real_wf_radgrid", ng-1 ,noccmax , lmax , expo,psi,aeval, occup     )
     endif
  endif
  

  if(iproc.eq.0 .and. dump_functions.eq.1) then
     open(unit=22,file='numerov_pseudo.dat')
     do igrid=1, Ngrid
        write(22,'(200(f20.10,1x))') rgrid(igrid), (psigrid_pseudo(igrid,j ), j=1,Nsol)   
     enddo
     close(unit=22)
  endif



  if(iproc.eq.0) print * , " routine GetExcitedOrbitalAsG  , generate  atom to  extract 1S " 
  psp_modifier=1
  
  call abs_generator_modified(iproc,atoms%nzatom(ity), atoms%nelpsp(ity),atoms%psppar(0,0,ity),&
       atoms%npspcode(ity),ng-1 ,noccmax , lmax , expo,psi,aeval, occup , psp_modifier  , &
       Nsol, abs_initial_L , Ngrid,Egrid,  rgrid , psigrid  )
  
  !! retrieve 1 s *r 

  
  do igrid=1,Ngrid
     psi1s(igrid) =   psigrid(igrid,1)  *rgrid(igrid)
  enddo


 
  if(iproc.eq.0) print * , " routine GetExcitedOrbitalAsG  , calculate  atom with  real-pot " 
  
  psp_modifier=1
  call abs_generator_modified(iproc,atoms%nzatom(ity), atoms%nelpsp(ity),atoms%psppar(0,0,ity),&
       atoms%npspcode(ity),ng-1 ,noccmax , lmax , expo,psi,aeval, occup , psp_modifier ,  &
       Nsol, abs_final_L , Ngrid,Egrid,  rgrid , psigrid  )

  

  if(iproc.eq.0 .and. dump_functions.eq.1) then
     open(unit=22,file='numerov.dat')
     do igrid=1, Ngrid
        write(22,'(200(f20.10,1x))') rgrid(igrid), (psigrid(igrid,j ), j=1,Nsol)   
     enddo
     close(unit=22)
  endif

  if(iproc.eq.0 .and. dump_functions.eq.1) then
     if(psp_modifier.eq.0) then
        call dump_gauwf_on_radgrid("pseudo_wf_radgrid", ng-1 ,noccmax , lmax , expo,psi,aeval, occup     )
     else
        call dump_gauwf_on_radgrid("real_wf_radgrid", ng-1 ,noccmax , lmax , expo,psi,aeval, occup     )
     endif
  endif
  
  real_start=-1
  do iocc=1, Nsol
     if((Egrid(iocc)+0.1).ge.Egrid_pseudo(1)) then
        real_start = iocc
        exit
     endif
  enddo
  
  if(real_start.eq.-1) then
     print *, " routine GetExcitedOrbitalAsG  ,  not found  eigenvalues of the real potential problem which could correspond "
     print *, " to the first one of the pseudo one , abs_final_L=", abs_final_L
     if (nproc > 1) call MPI_FINALIZE(ierr)
     stop
  endif
  
  if(iproc.eq.0) then 
     print *, " routine GetExcitedOrbitalAsG  ,  comparaison between  energies real and  pseudo "
     do iocc=1, Nsol
        if(iocc.lt.real_start) then
           print *,  iocc, Egrid(iocc) 
        else
           print *,  iocc, Egrid(iocc) , Egrid_pseudo(iocc-real_start +1)
        endif
     enddo
  endif
  
  if(iproc.eq.0) print *, " routine GetExcitedOrbitalAsG  , PROJECT 1s*r on pseudos "
  call find_pfproj( Nsol,Ngrid, rgrid, psi1s, psigrid, real_start, psigrid_pseudo)
  
  if(iproc.eq.0 .and. dump_functions.eq.1) then 
     open(unit=22,file='projres.dat')
     do igrid=1, Ngrid
        write(22,'(200(f20.10,1x))') rgrid(igrid),  psi1s(igrid), psigrid(igrid,1), psigrid_pseudo(igrid,1)
     enddo
     close(unit=22)
  endif
  
  
  
  expo_fine(1)=             expo(1)/3.0
  expo_fine(ng_fine) =       cradius*2.0      
  
  do ig=1, ng_fine
     expo_fine(ig) = exp( ( log(expo_fine(1))*(ng_fine-ig) +log(expo_fine(ng_fine))*(ig-1) )/(ng_fine-1))
  enddo
  
  call find_Scoeffs_grid(ng_fine-1,  expo_fine, Ngrid, rgrid, psigrid_pseudo(:,1)  , gcoeffs , abs_final_L  )
  

  if(iproc.eq.0 .and. dump_functions.eq.1)  call  dump_1gauwf_on_radgrid("proje_gau_proje_pseudo.dat",&
       ng_fine-1 , expo_fine,gcoeffs   ,   abs_final_L +1  )
  
 
  
  ! ----------------- Gabsorber --------------------------------------------------------
  Gabsorber%nat = 1
  allocate(Gabsorber%rxyz(3,Gabsorber%nat  ))
  allocate(Gabsorber%nshell(Gabsorber%nat ),stat=i_stat)

  Gabsorber%rxyz(:,1) = rxyz(:,in_iat_absorber )

  Gabsorber%nshell(1)=1
  Gabsorber%nshltot  =1

  allocate(Gabsorber%ndoc(1),stat=i_stat)
  allocate(Gabsorber%nam (1),stat=i_stat)


  call memocc(i_stat,Gabsorber%rxyz ,'Gabsorber%rxyz',subname)
  call memocc(i_stat,Gabsorber%nshell,'Gabsorber%nshell',subname)

  call memocc(i_stat,Gabsorber%ndoc,'Gabsorber%ndoc',subname)
  call memocc(i_stat,Gabsorber%nam,'Gabsorber%nam',subname)

  Gabsorber%nexpo=0
  Gabsorber%ncoeff=0


  Gabsorber%ndoc(1)  =  ng_fine
  Gabsorber%nam(1)   =  abs_final_l+1
  Gabsorber%nexpo         =  Gabsorber%nexpo+ ng_fine 
  Gabsorber%ncoeff        =  Gabsorber%ncoeff+2* abs_final_l + 1

  allocate(Gabsorber%psiat(Gabsorber%nexpo),stat=i_stat)
  allocate(Gabsorber%xp(Gabsorber%nexpo),stat=i_stat)

  call memocc(i_stat,Gabsorber%psiat , 'Gabsorber%psiat',subname)
  call memocc(i_stat,Gabsorber%xp    , 'Gabsorber%xp',subname)

  iexpo=0
  do ig=1,Gabsorber%ndoc(1)
     iexpo=iexpo+1
     Gabsorber%psiat(iexpo)=gcoeffs (ig-1)
     Gabsorber%xp(iexpo)=expo_fine(ig)
  end do

  print *,'expo',shape(expo),ng_fine,expo(:)

  ! -------------------------------------------------------------------------

  !fill the polarisation, hard coded for the moment
  Gabs_coeffs(1)=0.0_wp
  Gabs_coeffs(2)=0.0_wp
  Gabs_coeffs(3)=1.0_wp

  i_all=-product(shape(psigrid_pseudo))*kind(psigrid_pseudo)
  deallocate(psigrid_pseudo,stat=i_stat)
  call memocc(i_stat,i_all,'psigrid_pseudo',subname)

  
  i_all=-product(shape(Egrid))*kind(Egrid)
  deallocate(Egrid,stat=i_stat)
  call memocc(i_stat,i_all,'Egrid',subname)

  i_all=-product(shape(psigrid))*kind(psigrid)
  deallocate(psigrid,stat=i_stat)
  call memocc(i_stat,i_all,'psigrid',subname)
  
 
  i_all=-product(shape(occup))*kind(occup)
  deallocate(occup,stat=i_stat)
  call memocc(i_stat,i_all,'occup',subname)
  
  i_all=-product(shape(aeval))*kind(aeval)
  deallocate(aeval,stat=i_stat)
  call memocc(i_stat,i_all,'aeval',subname)

  i_all=-product(shape(gcoeffs))*kind(gcoeffs)
  deallocate(gcoeffs,stat=i_stat)
  call memocc(i_stat,i_all,'gcoeffs',subname) 
  
  i_all=-product(shape(psi))*kind(psi)
  deallocate(psi,stat=i_stat)
  call memocc(i_stat,i_all,'psi',subname) 
  
  i_all=-product(shape(expo))*kind(expo)
  deallocate(expo,stat=i_stat)
  call memocc(i_stat,i_all,'expo',subname)
    
  i_all=-product(shape(expo_fine))*kind(expo_fine)
  deallocate(expo_fine,stat=i_stat)
  call memocc(i_stat,i_all,'expo_fine',subname)
  
  i_all=-product(shape(Egrid_pseudo))*kind(Egrid_pseudo)
  deallocate(Egrid_pseudo,stat=i_stat)
  call memocc(i_stat,i_all,'Egrid_pseudo',subname)
  
  i_all=-product(shape(rgrid))*kind(rgrid)
  deallocate(rgrid,stat=i_stat)
  call memocc(i_stat,i_all,'rgrid',subname)
  
  i_all=-product(shape(psi1s))*kind(psi1s)
  deallocate(psi1s,stat=i_stat)
  call memocc(i_stat,i_all,'psi1s',subname)
  

  return
end subroutine GetExcitedOrbitalAsG
